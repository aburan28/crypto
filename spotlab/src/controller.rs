//! The controller: submit jobs, and reconcile every job's state with what is
//! running. Only the controller writes `record.json` and `result.json`, so a
//! single reconciler at a time keeps them consistent; a lease in the store
//! keeps a second controller (another laptop, a cron job) from acting at the
//! same time.

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;

use crate::provider::{
    aws::Aws, gcp::Gcp, local::Local, sort_offers, InstState, LaunchRequest, Offer, Provider,
};
use crate::record::*;
use crate::spec::{Config, JobSpec, Placement, Resources, JOB_SCHEMA, RESULT_SCHEMA};
use crate::store::Store;
use crate::util::{now, sha256_file, sha256_hex};

const BOOTSTRAP: &str = include_str!("bootstrap.sh");
const LEASE_KEY: &str = "controller/lease.json";

pub struct Controller {
    pub cfg: Config,
    pub store: Store,
    pub providers: Vec<Box<dyn Provider>>,
    pub owner: String,
}

#[derive(Serialize, Deserialize)]
struct Lease {
    owner: String,
    expires: f64,
}

impl Controller {
    pub fn new(cfg: Config) -> Result<Controller> {
        let store = Store::open(&cfg.store)?;
        let mut providers: Vec<Box<dyn Provider>> = Vec::new();
        if let Some(c) = &cfg.aws {
            providers.push(Box::new(Aws { cfg: c.clone() }));
        }
        if let Some(c) = &cfg.gcp {
            providers.push(Box::new(Gcp { cfg: c.clone() }));
        }
        if let Some(c) = &cfg.local {
            providers.push(Box::new(Local { cfg: c.clone() }));
        }
        // The lease belongs to the host; processes on one host are serialised by a file lock.
        let mut host = [0u8; 256];
        unsafe { libc::gethostname(host.as_mut_ptr() as *mut libc::c_char, host.len()) };
        let host = String::from_utf8_lossy(&host)
            .trim_end_matches('\0')
            .to_string();
        let owner = format!("{host}:{}", unsafe { libc::getuid() });
        Ok(Controller {
            cfg,
            store,
            providers,
            owner,
        })
    }

    fn provider(&self, name: &str) -> Option<&dyn Provider> {
        self.providers
            .iter()
            .find(|p| p.name() == name)
            .map(|b| b.as_ref())
    }

    // ---------------------------------------------------------------- submit

    /// Validate, upload local inputs, store the spec and a queued record.
    pub fn submit(&self, mut spec: JobSpec, base_dir: &Path) -> Result<String> {
        spec.validate()?;
        for i in spec.inputs.iter_mut() {
            if let Some(lf) = i.local_file.take() {
                let p = if Path::new(&lf).is_absolute() {
                    lf.clone().into()
                } else {
                    base_dir.join(&lf)
                };
                let sha = sha256_file(&p).with_context(|| format!("input {}", p.display()))?;
                let key = format!("blobs/{sha}");
                if !self.store.exists(&key)? {
                    self.store.put_file(&key, &p)?;
                }
                if !i.executable {
                    use std::os::unix::fs::PermissionsExt;
                    i.executable = std::fs::metadata(&p)?.permissions().mode() & 0o111 != 0;
                }
                i.store_key = Some(key);
                i.sha256 = Some(sha);
            } else if let Some(c) = &i.content {
                i.sha256 = Some(sha256_hex(c.as_bytes()));
            }
        }
        let (canonical, spec_sha) = spec.canonical()?;
        let job_id = format!(
            "j{}",
            &sha256_hex(format!("{spec_sha}{}{}", now(), std::process::id()).as_bytes())[..12]
        );
        self.store
            .put(&job_key(&job_id, "spec.json"), canonical.as_bytes())?;
        let t = now();
        let rec = JobRecord {
            schema: JOB_SCHEMA.into(),
            job_id: job_id.clone(),
            name: spec.name.clone(),
            spec_sha256: spec_sha,
            state: JobState::Queued,
            created_at: t,
            updated_at: t,
            attempt: 0,
            attempts: vec![],
            reason: Some("submitted".into()),
            labels: spec.labels.clone(),
        };
        self.store
            .put_json(&job_key(&job_id, "record.json"), &rec)?;
        Ok(job_id)
    }

    pub fn cancel(&self, job_id: &str) -> Result<()> {
        self.record(job_id)?;
        self.store.put(
            &job_key(job_id, "cancel"),
            format!("{}\n", now()).as_bytes(),
        )
    }

    pub fn record(&self, job_id: &str) -> Result<JobRecord> {
        self.store
            .get_json(&job_key(job_id, "record.json"))?
            .with_context(|| format!("no job {job_id}"))
    }

    pub fn spec(&self, job_id: &str) -> Result<JobSpec> {
        self.store
            .get_json(&job_key(job_id, "spec.json"))?
            .with_context(|| format!("job {job_id} has no spec"))
    }

    pub fn job_ids(&self) -> Result<Vec<String>> {
        let mut ids: Vec<String> = self
            .store
            .list("jobs/")?
            .iter()
            .filter_map(|k| {
                k.strip_prefix("jobs/")?
                    .strip_suffix("/record.json")
                    .map(str::to_string)
            })
            .filter(|id| !id.contains('/'))
            .collect();
        ids.sort();
        ids.dedup();
        Ok(ids)
    }

    // ---------------------------------------------------------------- offers

    /// Every offer that fits, across the allowed providers, cheapest first.
    /// Errors from one provider are returned beside the offers of the others.
    pub fn offers(&self, r: &Resources, p: &Placement) -> (Vec<Offer>, Vec<String>) {
        let mut all = Vec::new();
        let mut errors = Vec::new();
        for prov in &self.providers {
            if !p.providers.is_empty() && !p.providers.iter().any(|n| n == prov.name()) {
                continue;
            }
            match prov.offers(r, p) {
                Ok(v) => all.extend(v),
                Err(e) => errors.push(format!("{}: {e:#}", prov.name())),
            }
        }
        if let Some(max) = p.max_hourly_usd {
            all.retain(|o| o.hourly_usd <= max);
        }
        sort_offers(&mut all);
        (all, errors)
    }

    // ------------------------------------------------------------- reconcile

    fn take_lease(&self, ttl: f64) -> Result<bool> {
        if let Some(l) = self.store.get_json::<Lease>(LEASE_KEY)? {
            if l.owner != self.owner && l.expires > now() {
                return Ok(false);
            }
        }
        self.store.put_json(
            LEASE_KEY,
            &Lease {
                owner: self.owner.clone(),
                expires: now() + ttl,
            },
        )?;
        Ok(true)
    }

    /// One pass over every job. Returns a line per action taken.
    pub fn release_lease(&self) -> Result<()> {
        if let Some(l) = self.store.get_json::<Lease>(LEASE_KEY)? {
            if l.owner == self.owner {
                self.store.put_json(
                    LEASE_KEY,
                    &Lease {
                        owner: self.owner.clone(),
                        expires: now(),
                    },
                )?;
            }
        }
        Ok(())
    }

    pub fn reconcile(&self, lease_ttl_s: f64) -> Result<Vec<String>> {
        let _host_lock = HostLock::acquire(&self.store.url())?;
        if !self.take_lease(lease_ttl_s)? {
            return Ok(vec!["another controller holds the lease; skipped".into()]);
        }
        let mut events = Vec::new();
        let mut records: Vec<JobRecord> = Vec::new();
        for id in self.job_ids()? {
            match self.record(&id) {
                Ok(r) if !r.state.terminal() => records.push(r),
                Ok(_) => {}
                Err(e) => events.push(format!("{id}: unreadable record: {e:#}")),
            }
        }
        records.sort_by(|a, b| {
            a.created_at
                .partial_cmp(&b.created_at)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Pass 1: settle running attempts (frees budget before launching).
        for rec in records.iter_mut() {
            if rec.state == JobState::Running {
                match self.settle(rec) {
                    Ok(Some(ev)) => events.push(format!("{}: {ev}", rec.job_id)),
                    Ok(None) => {}
                    Err(e) => events.push(format!("{}: {e:#}", rec.job_id)),
                }
            }
        }
        // Pass 2: cancellations and launches, oldest job first, within the budget.
        let t = now();
        let mut running = 0usize;
        let mut hourly = 0.0;
        for r in records.iter().filter(|r| r.state == JobState::Running) {
            running += 1;
            hourly += r.current().map(|a| a.hourly_usd).unwrap_or(0.0);
        }
        let mut offer_cache: BTreeMap<String, (Vec<Offer>, Vec<String>)> = BTreeMap::new();
        for rec in records.iter_mut() {
            if rec.state.terminal() {
                continue;
            }
            if self.store.exists(&job_key(&rec.job_id, "cancel"))? {
                if rec.state == JobState::Running {
                    if let Some(a) = rec.current().cloned() {
                        if let Some(p) = self.provider(&a.instance.provider) {
                            p.terminate(&a.instance).ok();
                        }
                    }
                    if let Some(a) = rec.current_mut() {
                        a.ended_at.get_or_insert(t);
                        a.end_reason.get_or_insert("cancelled".into());
                    }
                }
                rec.state = JobState::Cancelled;
                rec.reason = Some("cancelled".into());
                self.finish(rec, None)?;
                events.push(format!("{}: cancelled", rec.job_id));
                continue;
            }
            if rec.state != JobState::Queued {
                continue;
            }
            let spec = self.spec(&rec.job_id)?;
            if let Some(why) = self.exhausted(rec, &spec, t) {
                rec.state = JobState::Dead;
                rec.reason = Some(why.clone());
                self.finish(rec, None)?;
                events.push(format!("{}: dead: {why}", rec.job_id));
                continue;
            }
            let key = serde_json::to_string(&(&spec.resources, &spec.placement))?;
            let (offers, errors) = offer_cache
                .entry(key)
                .or_insert_with(|| self.offers(&spec.resources, &spec.placement))
                .clone();
            match self.launch(rec, &spec, &offers, running, hourly) {
                Ok(Some(price)) => {
                    running += 1;
                    hourly += price;
                    let a = rec.current().expect("just launched");
                    events.push(format!(
                        "{}: attempt {} launched on {} {} in {} at ${:.4}/h",
                        rec.job_id,
                        rec.attempt,
                        a.instance.provider,
                        a.instance.instance_type,
                        a.instance.zone,
                        price
                    ));
                }
                Ok(None) => {
                    let mut why = rec.reason.clone().unwrap_or_default();
                    if !errors.is_empty() {
                        why = format!("{why}; provider errors: {}", errors.join("; "));
                        rec.reason = Some(why.clone());
                        rec.updated_at = now();
                        self.save(rec)?;
                    }
                    events.push(format!("{}: waiting: {why}", rec.job_id));
                }
                Err(e) => events.push(format!("{}: launch failed: {e:#}", rec.job_id)),
            }
        }
        Ok(events)
    }

    fn save(&self, rec: &JobRecord) -> Result<()> {
        self.store
            .put_json(&job_key(&rec.job_id, "record.json"), rec)
    }

    fn exhausted(&self, rec: &JobRecord, spec: &JobSpec, t: f64) -> Option<String> {
        if rec.attempt >= spec.limits.max_attempts {
            return Some(format!("{} attempts used", rec.attempt));
        }
        if rec.total_hours(t) >= spec.limits.max_total_hours {
            return Some(format!(
                "{:.2} h used of {}",
                rec.total_hours(t),
                spec.limits.max_total_hours
            ));
        }
        if let Some(max) = spec.limits.max_cost_usd {
            if rec.cost_usd(t) >= max {
                return Some(format!("${:.2} spent of ${max:.2}", rec.cost_usd(t)));
            }
        }
        None
    }

    /// Look at a running attempt; move the job on if the attempt ended.
    fn settle(&self, rec: &mut JobRecord) -> Result<Option<String>> {
        let t = now();
        let Some(cur) = rec.current().cloned() else {
            rec.state = JobState::Queued;
            self.save(rec)?;
            return Ok(Some("running without an attempt; requeued".into()));
        };
        let spec = self.spec(&rec.job_id)?;
        let prov = self.provider(&cur.instance.provider);
        let exit: Option<ExitDoc> =
            self.store
                .get_json(&attempt_key(&rec.job_id, rec.attempt, "exit.json"))?;
        if let Some(x) = exit {
            if let Some(p) = prov {
                // The bootstrap powers the instance off itself; this is belt and braces.
                p.terminate(&cur.instance).ok();
            }
            let a = rec.current_mut().expect("current");
            a.ended_at = Some(x.ended_at.max(a.launched_at));
            a.end_reason = Some(x.reason.clone());
            let ev = match x.reason.as_str() {
                "completed" => {
                    let ok = x.exit_code == Some(0);
                    rec.state = if ok {
                        JobState::Succeeded
                    } else {
                        JobState::Failed
                    };
                    rec.reason = Some(match (x.exit_code, x.signal) {
                        (Some(c), _) => format!("command exited {c}"),
                        (None, Some(s)) => format!("command killed by signal {s}"),
                        _ => "command ended".into(),
                    });
                    self.finish(rec, Some(&x))?;
                    format!(
                        "attempt {} completed: {}",
                        rec.attempt,
                        rec.reason.clone().unwrap_or_default()
                    )
                }
                "preempted" => {
                    rec.state = JobState::Queued;
                    rec.reason = Some(format!(
                        "attempt {} preempted; resumes from checkpoint {}",
                        rec.attempt,
                        x.last_checkpoint_seq
                            .map(|s| s.to_string())
                            .unwrap_or_else(|| "none (starts over)".into())
                    ));
                    self.save(rec)?;
                    rec.reason.clone().unwrap()
                }
                other => {
                    rec.state = JobState::Queued;
                    rec.reason = Some(format!(
                        "attempt {} ended {other}: {}",
                        rec.attempt,
                        x.error.clone().unwrap_or_default()
                    ));
                    self.save(rec)?;
                    rec.reason.clone().unwrap()
                }
            };
            return Ok(Some(ev));
        }

        let state = match prov {
            Some(p) => p.state(&cur.instance)?,
            None => InstState::Gone(format!("provider {} not configured", cur.instance.provider)),
        };
        match state {
            InstState::Gone(why) => {
                let a = rec.current_mut().expect("current");
                a.ended_at = Some(t);
                a.end_reason = Some(format!("lost: {why}"));
                rec.state = JobState::Queued;
                rec.reason = Some(format!(
                    "attempt {} instance gone without an exit record ({why}); resuming",
                    rec.attempt
                ));
                self.save(rec)?;
                Ok(Some(rec.reason.clone().unwrap()))
            }
            InstState::Pending | InstState::Running => {
                let hb: Option<Heartbeat> = self.store.get_json(&attempt_key(
                    &rec.job_id,
                    rec.attempt,
                    "heartbeat.json",
                ))?;
                let last = hb.map(|h| h.at).unwrap_or(cur.launched_at);
                let over_budget = self.exhausted(rec, &spec, t);
                if t - last > spec.limits.heartbeat_stale_s as f64 || over_budget.is_some() {
                    if let Some(p) = prov {
                        p.terminate(&cur.instance)?;
                    }
                    let a = rec.current_mut().expect("current");
                    a.ended_at = Some(t);
                    match over_budget {
                        Some(why) => {
                            a.end_reason = Some("deadline".into());
                            rec.state = JobState::Dead;
                            rec.reason = Some(why);
                            self.finish(rec, None)?;
                        }
                        None => {
                            a.end_reason = Some("stale".into());
                            rec.state = JobState::Queued;
                            rec.reason = Some(format!(
                                "attempt {} silent for {:.0} s; terminated, resuming",
                                rec.attempt,
                                t - last
                            ));
                            self.save(rec)?;
                        }
                    }
                    return Ok(Some(rec.reason.clone().unwrap()));
                }
                Ok(None)
            }
        }
    }

    /// Launch the next attempt on the cheapest offer that accepts it.
    /// `Ok(Some(price))` when launched, `Ok(None)` when it has to wait.
    fn launch(
        &self,
        rec: &mut JobRecord,
        spec: &JobSpec,
        offers: &[Offer],
        running: usize,
        hourly: f64,
    ) -> Result<Option<f64>> {
        if offers.is_empty() {
            rec.reason = Some("no offer fits the resources, placement and price cap (GCP and local entries need a spot_hourly_usd)".into());
            return Ok(None);
        }
        if running >= self.cfg.budget.max_running {
            rec.reason = Some(format!("budget: {running} instances already running"));
            return Ok(None);
        }
        let attempt = rec.attempt + 1;
        let remaining_h = (spec.limits.max_total_hours - rec.total_hours(now())).max(0.05);
        let mut tried = Vec::new();
        for offer in offers.iter().take(6) {
            if hourly + offer.hourly_usd > self.cfg.budget.max_hourly_usd {
                tried.push(format!(
                    "{} {}: over the ${}/h budget",
                    offer.provider, offer.instance_type, self.cfg.budget.max_hourly_usd
                ));
                continue;
            }
            if offer.provider != "local"
                && !self.store.exists(&format!("bin/spotlab-{}", offer.arch))?
            {
                tried.push(format!(
                    "{}: no agent binary for {} (run `spotlab publish-agent`)",
                    offer.instance_type, offer.arch
                ));
                continue;
            }
            let Some(p) = self.provider(&offer.provider) else {
                continue;
            };
            let req = LaunchRequest {
                name: format!("spotlab-{}-a{attempt}", rec.job_id),
                job_id: rec.job_id.clone(),
                attempt,
                bootstrap: self.bootstrap(rec, offer, attempt, remaining_h),
                max_run_s: (remaining_h * 3600.0).ceil() as u64,
                max_hourly_usd: spec.placement.max_hourly_usd,
                store_url: self.store.url(),
            };
            // Record the attempt first: the agent checks that it is current.
            let prev = rec.clone();
            rec.attempt = attempt;
            rec.state = JobState::Running;
            rec.reason = Some(format!(
                "launching on {} {}",
                offer.provider, offer.instance_type
            ));
            rec.updated_at = now();
            self.save(rec)?;
            match p.launch(offer, &req) {
                Ok(inst) => {
                    rec.attempts.push(AttemptRecord {
                        attempt,
                        instance: inst,
                        hourly_usd: offer.hourly_usd,
                        launched_at: now(),
                        ended_at: None,
                        end_reason: None,
                    });
                    rec.reason = None;
                    rec.updated_at = now();
                    self.save(rec)?;
                    return Ok(Some(offer.hourly_usd));
                }
                Err(e) => {
                    *rec = prev;
                    self.save(rec)?;
                    tried.push(format!(
                        "{} {} in {}: {e:#}",
                        offer.provider, offer.instance_type, offer.zone
                    ));
                }
            }
        }
        rec.reason = Some(format!("no launch succeeded: {}", tried.join("; ")));
        rec.updated_at = now();
        self.save(rec)?;
        Ok(None)
    }

    fn bootstrap(&self, rec: &JobRecord, offer: &Offer, attempt: u32, remaining_h: f64) -> String {
        let extra = match offer.provider.as_str() {
            "aws" => self
                .cfg
                .aws
                .as_ref()
                .map(|c| c.extra_bootstrap.join("\n"))
                .unwrap_or_default(),
            "gcp" => self
                .cfg
                .gcp
                .as_ref()
                .map(|c| c.extra_bootstrap.join("\n"))
                .unwrap_or_default(),
            _ => String::new(),
        };
        BOOTSTRAP
            .replace("@STORE@", &self.store.url())
            .replace("@JOB@", &rec.job_id)
            .replace("@ATTEMPT@", &attempt.to_string())
            .replace("@PROVIDER@", &offer.provider)
            .replace(
                "@DEADLINE_MIN@",
                &((remaining_h * 60.0).ceil() as u64 + 5).to_string(),
            )
            .replace("@EXTRA@", &extra)
    }

    /// Write the result document and the final record.
    fn finish(&self, rec: &mut JobRecord, last: Option<&ExitDoc>) -> Result<()> {
        let t = now();
        rec.updated_at = t;
        self.save(rec)?;
        let mut attempts = Vec::new();
        let mut resumed = false;
        for a in &rec.attempts {
            let exit: Option<ExitDoc> = self
                .store
                .get_json(&attempt_key(&rec.job_id, a.attempt, "exit.json"))
                .unwrap_or(None);
            resumed |= exit.as_ref().and_then(|x| x.restored_seq).is_some();
            attempts.push(AttemptSummary {
                attempt: a.attempt,
                instance: a.instance.clone(),
                hourly_usd: a.hourly_usd,
                launched_at: a.launched_at,
                ended_at: a.ended_at,
                end_reason: a.end_reason.clone(),
                exit,
            });
        }
        let preemptions = rec
            .attempts
            .iter()
            .filter(|a| matches!(a.end_reason.as_deref(), Some(r) if r == "preempted" || r == "stale" || r.starts_with("lost")))
            .count() as u32;
        let doc = ResultDoc {
            schema: RESULT_SCHEMA.into(),
            job_id: rec.job_id.clone(),
            name: rec.name.clone(),
            spec_sha256: rec.spec_sha256.clone(),
            status: rec.state.as_str().into(),
            exit_code: last.and_then(|x| x.exit_code),
            reason: rec.reason.clone(),
            attempts,
            resumed,
            preemptions,
            total_hours: rec.total_hours(t),
            cost_usd_estimate: rec.cost_usd(t),
            outputs: last.map(|x| x.outputs.clone()).unwrap_or_default(),
            evidence: EVIDENCE.into(),
            labels: rec.labels.clone(),
        };
        let key = job_key(&rec.job_id, "result.json");
        if self.store.exists(&key)? {
            bail!("result for {} already written", rec.job_id);
        }
        self.store.put_json(&key, &doc)
    }
}

/// An exclusive `flock` held for one reconcile pass, so two processes on the
/// same host (a CLI call and an MCP server's loop) never interleave.
struct HostLock(libc::c_int);

impl HostLock {
    fn acquire(store_url: &str) -> Result<HostLock> {
        let dir = std::env::temp_dir().join("spotlab-locks");
        std::fs::create_dir_all(&dir)?;
        let path = dir.join(format!("{}.lock", &sha256_hex(store_url.as_bytes())[..16]));
        let c = std::ffi::CString::new(path.as_os_str().as_encoded_bytes())?;
        let fd = unsafe {
            libc::open(
                c.as_ptr(),
                libc::O_CREAT | libc::O_RDWR | libc::O_CLOEXEC,
                0o600,
            )
        };
        if fd < 0 {
            bail!(
                "open {}: {}",
                path.display(),
                std::io::Error::last_os_error()
            );
        }
        if unsafe { libc::flock(fd, libc::LOCK_EX) } != 0 {
            unsafe { libc::close(fd) };
            bail!(
                "flock {}: {}",
                path.display(),
                std::io::Error::last_os_error()
            );
        }
        Ok(HostLock(fd))
    }
}

impl Drop for HostLock {
    fn drop(&mut self) {
        unsafe {
            libc::flock(self.0, libc::LOCK_UN);
            libc::close(self.0);
        }
    }
}
