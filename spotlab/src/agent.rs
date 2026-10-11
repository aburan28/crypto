//! The agent: runs one attempt of one job on the instance it was launched on.
//!
//! It restores the newest checkpoint, starts the command, uploads the
//! checkpoint directory whenever it changes (at most every
//! `checkpoint.interval_s`), and polls the provider's interruption notice.
//! On a notice it sends the command SIGTERM, waits `grace_s`, uploads the
//! final checkpoint and records the attempt as preempted. The controller then
//! launches the next attempt somewhere else, which resumes from that upload.
//!
//! The command's side of the contract: keep resumable state in
//! `$SPOTLAB_CHECKPOINT_DIR`, replace files there atomically (write, then
//! rename), restore from it at start-up when it is not empty, and save on
//! SIGTERM. Final results go in `$SPOTLAB_OUTPUT_DIR`.

use anyhow::{Context, Result};
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

use crate::record::*;
use crate::spec::JobSpec;
use crate::store::Store;
use crate::util::{now, run, s, sha256_file};

pub struct AgentArgs {
    pub store: String,
    pub job: String,
    pub attempt: u32,
    pub provider: String,
    pub workdir: Option<PathBuf>,
    pub preempt_file: Option<PathBuf>,
}

static TERM: AtomicBool = AtomicBool::new(false);

extern "C" fn on_term(_: libc::c_int) {
    TERM.store(true, Ordering::SeqCst);
}

/// What the interruption check found.
#[derive(Debug, PartialEq)]
enum Notice {
    None,
    /// The provider will reclaim the instance soon: stop and save.
    Interrupt(String),
    /// Elevated risk (AWS rebalance recommendation): save now, keep running.
    Rebalance,
}

struct NoticeSource {
    provider: String,
    preempt_file: Option<PathBuf>,
    aws_token: Option<(String, Instant)>,
}

impl NoticeSource {
    fn poll(&mut self) -> Notice {
        if TERM.load(Ordering::SeqCst) {
            return Notice::Interrupt("agent received SIGTERM (instance shutting down)".into());
        }
        if let Some(f) = &self.preempt_file {
            if f.exists() {
                return Notice::Interrupt("preempt file".into());
            }
        }
        match self.provider.as_str() {
            "aws" => self.poll_aws(),
            "gcp" => {
                let out = run(
                    &s(&[
                        "curl",
                        "-s",
                        "-m",
                        "2",
                        "-H",
                        "Metadata-Flavor: Google",
                        "http://metadata.google.internal/computeMetadata/v1/instance/preempted",
                    ]),
                    None,
                );
                match out {
                    Ok(o) if o.ok() && o.stdout_str().trim() == "TRUE" => {
                        Notice::Interrupt("gce preempted".into())
                    }
                    _ => Notice::None,
                }
            }
            _ => Notice::None,
        }
    }

    fn poll_aws(&mut self) -> Notice {
        let fresh =
            matches!(&self.aws_token, Some((_, t)) if t.elapsed() < Duration::from_secs(240));
        if !fresh {
            let tok = run(
                &s(&[
                    "curl",
                    "-s",
                    "-m",
                    "2",
                    "-X",
                    "PUT",
                    "http://169.254.169.254/latest/api/token",
                    "-H",
                    "X-aws-ec2-metadata-token-ttl-seconds: 300",
                ]),
                None,
            );
            self.aws_token = match tok {
                Ok(o) if o.ok() && !o.stdout.is_empty() => Some((o.stdout_str(), Instant::now())),
                _ => None,
            };
        }
        let Some((tok, _)) = &self.aws_token else {
            return Notice::None;
        };
        let get = |path: &str| -> bool {
            let hdr = format!("X-aws-ec2-metadata-token: {tok}");
            let url = format!("http://169.254.169.254/latest/meta-data/{path}");
            matches!(
                run(&s(&["curl", "-s", "-m", "2", "-o", "/dev/null", "-w", "%{http_code}", "-H", &hdr, &url]), None),
                Ok(o) if o.stdout_str().trim() == "200"
            )
        };
        if get("spot/instance-action") {
            Notice::Interrupt("ec2 spot interruption notice".into())
        } else if get("events/recommendations/rebalance") {
            Notice::Rebalance
        } else {
            Notice::None
        }
    }
}

struct Attempt {
    store: Store,
    spec: JobSpec,
    job: String,
    attempt: u32,
    ckpt_dir: PathBuf,
    out_dir: PathBuf,
    log_dir: PathBuf,
    started: f64,
    restored_seq: Option<u64>,
    last_seq: Option<u64>,
    last_digest: Option<String>,
}

impl Attempt {
    fn heartbeat(&self, phase: &str) {
        let hb = Heartbeat {
            attempt: self.attempt,
            at: now(),
            phase: phase.into(),
            elapsed_s: now() - self.started,
            last_checkpoint_seq: self.last_seq,
        };
        if let Err(e) = self
            .store
            .put_json(&attempt_key(&self.job, self.attempt, "heartbeat.json"), &hb)
        {
            eprintln!("spotlab agent: heartbeat failed: {e:#}");
        }
    }

    /// False when this attempt may no longer write: the controller moved on or the job was cancelled.
    fn still_current(&self) -> bool {
        match self
            .store
            .get_json::<JobRecord>(&job_key(&self.job, "record.json"))
        {
            Ok(Some(r)) => {
                r.attempt == self.attempt
                    && !r.state.terminal()
                    && !self
                        .store
                        .exists(&job_key(&self.job, "cancel"))
                        .unwrap_or(false)
            }
            // A store hiccup is not a reason to kill a long run.
            _ => true,
        }
    }

    /// Upload the checkpoint directory if its contents changed since the last upload or restore.
    fn sync(&mut self) -> Result<bool> {
        let (files, digest) = hash_dir(&self.ckpt_dir)?;
        if files.is_empty() || self.last_digest.as_deref() == Some(digest.as_str()) {
            return Ok(false);
        }
        if !self.still_current() {
            return Ok(false);
        }
        let existing = checkpoint_seqs(&self.store, &self.job)?
            .last()
            .copied()
            .unwrap_or(0);
        let seq = existing.max(self.last_seq.unwrap_or(0)) + 1;
        let m = upload_checkpoint(&self.store, &self.job, seq, self.attempt, &self.ckpt_dir)?;
        self.last_seq = Some(seq);
        self.last_digest = Some(m.digest);
        prune_checkpoints(&self.store, &self.job, self.spec.checkpoint.keep.max(1))?;
        eprintln!(
            "spotlab agent: checkpoint {seq} uploaded ({} files)",
            m.files.len()
        );
        Ok(true)
    }

    fn upload_logs(&self) {
        for name in ["stdout.log", "stderr.log"] {
            let p = self.log_dir.join(name);
            if p.exists() {
                if let Err(e) = self
                    .store
                    .put_file(&attempt_key(&self.job, self.attempt, name), &p)
                {
                    eprintln!("spotlab agent: upload {name}: {e:#}");
                }
            }
        }
    }

    fn upload_outputs(&self) -> Result<Vec<FileEntry>> {
        let (files, _) = hash_dir(&self.out_dir)?;
        for f in &files {
            self.store.put_file(
                &job_key(&self.job, &format!("output/{}", f.path)),
                &self.out_dir.join(&f.path),
            )?;
        }
        Ok(files)
    }
}

pub fn run_agent(a: AgentArgs) -> Result<()> {
    unsafe {
        libc::signal(libc::SIGTERM, on_term as *const () as libc::sighandler_t);
        libc::signal(libc::SIGINT, on_term as *const () as libc::sighandler_t);
    }
    let store = Store::open(&a.store)?;
    let exit_key = attempt_key(&a.job, a.attempt, "exit.json");
    if store.exists(&exit_key)? {
        eprintln!(
            "spotlab agent: attempt {} already ended; nothing to do",
            a.attempt
        );
        return Ok(());
    }
    let spec: JobSpec = store
        .get_json(&job_key(&a.job, "spec.json"))?
        .with_context(|| format!("job {} has no spec", a.job))?;
    let workdir = a
        .workdir
        .clone()
        .unwrap_or_else(|| PathBuf::from(format!("/var/lib/spotlab/{}/a{:04}", a.job, a.attempt)));
    let mut at = Attempt {
        store: store.clone(),
        spec,
        job: a.job.clone(),
        attempt: a.attempt,
        ckpt_dir: workdir.join("checkpoint"),
        out_dir: workdir.join("output"),
        log_dir: workdir.join("logs"),
        started: now(),
        restored_seq: None,
        last_seq: None,
        last_digest: None,
    };
    let work = workdir.join("work");
    for d in [&work, &at.ckpt_dir, &at.out_dir, &at.log_dir] {
        std::fs::create_dir_all(d)?;
    }

    let end = |at: &Attempt,
               reason: &str,
               code: Option<i32>,
               sig: Option<i32>,
               ru: Rusage,
               err: Option<String>,
               outputs: Vec<FileEntry>| {
        let doc = ExitDoc {
            attempt: at.attempt,
            reason: reason.into(),
            exit_code: code,
            signal: sig,
            started_at: at.started,
            ended_at: now(),
            wall_s: now() - at.started,
            user_cpu_s: ru.user,
            sys_cpu_s: ru.sys,
            max_rss_kb: ru.maxrss_kb,
            restored_seq: at.restored_seq,
            last_checkpoint_seq: at.last_seq,
            error: err,
            outputs,
        };
        at.upload_logs();
        if let Err(e) = at
            .store
            .put_json(&attempt_key(&at.job, at.attempt, "exit.json"), &doc)
        {
            eprintln!("spotlab agent: could not write exit.json: {e:#}");
        }
    };

    if !at.still_current() {
        end(
            &at,
            "fenced",
            None,
            None,
            Rusage::default(),
            Some("not the current attempt".into()),
            vec![],
        );
        return Ok(());
    }
    at.heartbeat("fetch");
    if let Err(e) = prepare(&mut at, &work) {
        end(
            &at,
            "error",
            None,
            None,
            Rusage::default(),
            Some(format!("{e:#}")),
            vec![],
        );
        return Err(e);
    }

    // Start the command in its own process group so signals reach all of it.
    let cwd = match &at.spec.command.cwd {
        Some(c) => work.join(c),
        None => work.clone(),
    };
    let mut cmd = Command::new(&at.spec.command.argv[0]);
    cmd.args(&at.spec.command.argv[1..])
        .current_dir(&cwd)
        .envs(&at.spec.command.env)
        .env("SPOTLAB_JOB_ID", &at.job)
        .env("SPOTLAB_ATTEMPT", at.attempt.to_string())
        .env("SPOTLAB_CHECKPOINT_DIR", &at.ckpt_dir)
        .env("SPOTLAB_OUTPUT_DIR", &at.out_dir)
        .env(
            "SPOTLAB_RESUMED",
            if at.restored_seq.is_some() { "1" } else { "0" },
        )
        .env(
            "SPOTLAB_RESTORED_SEQ",
            at.restored_seq.map(|s| s.to_string()).unwrap_or_default(),
        )
        .stdin(Stdio::null())
        .stdout(std::fs::File::create(at.log_dir.join("stdout.log"))?)
        .stderr(std::fs::File::create(at.log_dir.join("stderr.log"))?)
        .process_group(0);
    let child = match cmd.spawn() {
        Ok(c) => c,
        Err(e) => {
            let msg = format!("spawn {:?}: {e}", at.spec.command.argv[0]);
            end(
                &at,
                "error",
                None,
                None,
                Rusage::default(),
                Some(msg.clone()),
                vec![],
            );
            anyhow::bail!(msg);
        }
    };
    let pid = child.id() as i32;
    at.heartbeat("running");

    let mut notices = NoticeSource {
        provider: a.provider.clone(),
        preempt_file: a.preempt_file.clone(),
        aws_token: None,
    };
    let interval = Duration::from_secs(at.spec.checkpoint.interval_s);
    let mut last_sync = Instant::now();
    let mut last_hb = Instant::now();
    let mut last_fence = Instant::now();
    let mut last_poll = Instant::now() - Duration::from_secs(10);
    let mut rebalance_saved = false;

    loop {
        if let Some((st, ru)) = wait_child(pid, false) {
            // The command finished on its own.
            if let Err(e) = at.sync() {
                eprintln!("spotlab agent: final checkpoint: {e:#}");
            }
            let outputs = at.upload_outputs().unwrap_or_else(|e| {
                eprintln!("spotlab agent: outputs: {e:#}");
                vec![]
            });
            end(&at, "completed", st.code, st.signal, ru, None, outputs);
            return Ok(());
        }
        if last_poll.elapsed() >= Duration::from_secs(2) {
            last_poll = Instant::now();
            match notices.poll() {
                Notice::Interrupt(why) => {
                    eprintln!("spotlab agent: {why}; stopping the command");
                    at.heartbeat("preempting");
                    let (st, ru) = stop(pid, Duration::from_secs(at.spec.checkpoint.grace_s));
                    if let Err(e) = at.sync() {
                        eprintln!("spotlab agent: checkpoint on preemption failed: {e:#}");
                    }
                    end(&at, "preempted", st.code, st.signal, ru, Some(why), vec![]);
                    return Ok(());
                }
                Notice::Rebalance if !rebalance_saved => {
                    rebalance_saved = true;
                    eprintln!("spotlab agent: rebalance recommendation; checkpointing now");
                    at.sync().ok();
                    last_sync = Instant::now();
                }
                _ => {}
            }
        }
        if last_fence.elapsed() >= Duration::from_secs(30) {
            last_fence = Instant::now();
            if !at.still_current() {
                eprintln!("spotlab agent: superseded or cancelled; stopping");
                let (st, ru) = stop(pid, Duration::from_secs(at.spec.checkpoint.grace_s));
                end(
                    &at,
                    "fenced",
                    st.code,
                    st.signal,
                    ru,
                    Some("superseded or cancelled".into()),
                    vec![],
                );
                return Ok(());
            }
        }
        if last_sync.elapsed() >= interval {
            last_sync = Instant::now();
            if let Err(e) = at.sync() {
                eprintln!("spotlab agent: checkpoint upload failed (will retry): {e:#}");
            }
        }
        if last_hb.elapsed() >= Duration::from_secs(30) {
            last_hb = Instant::now();
            at.heartbeat("running");
        }
        std::thread::sleep(Duration::from_millis(250));
    }
}

/// Inputs into the working directory; newest checkpoint into the checkpoint directory.
fn prepare(at: &mut Attempt, work: &Path) -> Result<()> {
    for i in at.spec.inputs.clone() {
        let dst = work.join(crate::util::safe_rel(&i.path)?);
        if let Some(d) = dst.parent() {
            std::fs::create_dir_all(d)?;
        }
        if let Some(c) = &i.content {
            std::fs::write(&dst, c)?;
        } else if let Some(k) = &i.store_key {
            if !at.store.get_file(k, &dst)? {
                anyhow::bail!("input {}: {} not in the store", i.path, k);
            }
        } else {
            anyhow::bail!("input {}: local_file was not uploaded at submit", i.path);
        }
        if let Some(want) = &i.sha256 {
            let got = sha256_file(&dst)?;
            if &got != want {
                anyhow::bail!("input {}: sha256 {got} does not match {want}", i.path);
            }
        }
        if i.executable {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(0o755))?;
        }
    }
    if let Some(&seq) = checkpoint_seqs(&at.store, &at.job)?.last() {
        let m = restore_checkpoint(&at.store, &at.job, seq, &at.ckpt_dir)?;
        eprintln!(
            "spotlab agent: restored checkpoint {seq} ({} files, from attempt {})",
            m.files.len(),
            m.attempt
        );
        at.restored_seq = Some(seq);
        at.last_seq = Some(seq);
        at.last_digest = Some(m.digest);
    }
    Ok(())
}

#[derive(Default, Clone, Copy)]
struct Rusage {
    user: f64,
    sys: f64,
    maxrss_kb: i64,
}

struct WaitStatus {
    code: Option<i32>,
    signal: Option<i32>,
}

/// Reap the command with `wait4`, so the rusage is the command's own and not
/// that of the agent's other children (curl, the cloud CLIs).
fn wait_child(pid: i32, block: bool) -> Option<(WaitStatus, Rusage)> {
    let mut status = 0;
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    let flags = if block { 0 } else { libc::WNOHANG };
    let r = unsafe { libc::wait4(pid, &mut status, flags, &mut ru) };
    if r != pid {
        return None;
    }
    let tv = |t: libc::timeval| t.tv_sec as f64 + t.tv_usec as f64 / 1e6;
    let st = if libc::WIFEXITED(status) {
        WaitStatus {
            code: Some(libc::WEXITSTATUS(status)),
            signal: None,
        }
    } else if libc::WIFSIGNALED(status) {
        WaitStatus {
            code: None,
            signal: Some(libc::WTERMSIG(status)),
        }
    } else {
        WaitStatus {
            code: None,
            signal: None,
        }
    };
    Some((
        st,
        Rusage {
            user: tv(ru.ru_utime),
            sys: tv(ru.ru_stime),
            maxrss_kb: ru.ru_maxrss as i64,
        },
    ))
}

/// SIGTERM the command's process group, wait up to `grace`, then SIGKILL.
fn stop(pid: i32, grace: Duration) -> (WaitStatus, Rusage) {
    unsafe {
        libc::kill(-pid, libc::SIGTERM);
    }
    let t0 = Instant::now();
    while t0.elapsed() < grace {
        if let Some(r) = wait_child(pid, false) {
            return r;
        }
        std::thread::sleep(Duration::from_millis(100));
    }
    unsafe {
        libc::kill(-pid, libc::SIGKILL);
    }
    wait_child(pid, true).unwrap_or((
        WaitStatus {
            code: None,
            signal: Some(libc::SIGKILL),
        },
        Rusage::default(),
    ))
}
