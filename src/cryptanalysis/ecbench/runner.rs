//! The session runner.
//!
//! ```text
//! spec.json ──ecbench run──► <out>/
//!                              spec.json      byte copy of the spec
//!                              plan.json      arms, workloads, every execution in order
//!                              host.json      the host capsule
//!                              session.json   identity, reservation, preflight, hashes, status
//!                              records.jsonl  one sealed record per execution, in order
//!                              exec/<seq>.stderr   a child's stderr, when it wrote any
//! ```
//!
//! Each execution is a fresh process (`ecbench exec`), so one-time
//! costs — page faults, lazy statics, a cold allocator — land on every
//! execution alike instead of on whichever arm ran first.  The child is
//! pinned and its memory bound between fork and exec, with a pinned
//! environment; the runner waits for it blocked (no polling on any CPU),
//! a watchdog thread kills its process group at the timeout, and the
//! runner verifies the answer in its own process.  A record is written
//! and flushed for every execution, warm-ups, errors and timeouts
//! included, before the next one starts.
//!
//! The output directory is claimed with an atomic `mkdir` while the lock
//! is held, so two sessions queued for one directory cannot both write
//! it.  An interruption (`signals`) or an error leaves the session marked
//! `interrupted`, every evicted thread restored, and no child running.

use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant};

use serde::{Deserialize, Serialize};
use serde_json::json;

use crate::cryptanalysis::ecbench::canonical::{sha256_hex, short_id};
use crate::cryptanalysis::ecbench::host::{capture, HostCapsule};
use crate::cryptanalysis::ecbench::isolation::{
    allowed_cpus, conditions, evict, leave_reserved, plan_cpus, preflight, ticks_per_second,
    topology, BenchLock, CpuPlan, CpuRequest, Eviction, EvictionSummary, Preflight, Thresholds,
};
use crate::cryptanalysis::ecbench::methods::MethodSpec;
use crate::cryptanalysis::ecbench::record::{
    floor_s, grade, Boundaries, ChildInput, ChildOutput, Cost, GradeInput, IsolationRecord,
    Outcome, Record, Timing, CHILD_INPUT_SCHEMA, GRADING_VERSION, RECORD_SCHEMA, UNIT,
};
use crate::cryptanalysis::ecbench::signals;
use crate::cryptanalysis::ecbench::spec::{plan, Execution, Plan, Spec};
use crate::cryptanalysis::ecbench::workload::{CurveSpec, Workload};

pub const SESSION_SCHEMA: &str = "ecbench.session/v1";

/// Where and how a session runs; recorded whole in `session.json`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RunOptions {
    pub cpus: CpuRequest,
    pub lock_path: String,
    pub wait_for_lock: bool,
    /// Record a refused preflight and continue (every run then stays
    /// below L2), instead of refusing to start.
    pub allow_busy: bool,
    pub thresholds: Thresholds,
    /// The `ecbench` binary each measured child execs.  On Linux this is
    /// `/proc/self/exe`, so a rebuild during the session cannot swap the
    /// code under the hash recorded at its start.
    pub exe: PathBuf,
}

fn default_hz() -> u64 {
    100
}
fn grading_v1() -> u32 {
    1
}

/// The session document.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Session {
    pub schema: String,
    pub session_id: String,
    pub label: String,
    pub spec_id: String,
    pub spec_sha256: String,
    pub env_class_id: String,
    pub host_sha256: String,
    pub plan_sha256: String,
    pub binary_sha256: Option<String>,
    pub git_commit: Option<String>,
    pub git_dirty: Option<bool>,
    pub options: RunOptions,
    pub cpu_plan: Option<CpuPlan>,
    pub preflight: Option<Preflight>,
    pub eviction: Option<EvictionSummary>,
    /// The runner shared the run CPU (it had nowhere else to go).
    #[serde(default)]
    pub runner_on_run_cpu: bool,
    /// The CPUs the runner moved itself to.
    #[serde(default)]
    pub runner_cpus: Option<Vec<u32>>,
    /// `_SC_CLK_TCK` on the host, which grading divides ticks by.
    #[serde(default = "default_hz")]
    pub ticks_per_second: u64,
    /// The rules the records' levels were graded under; an audit regrades
    /// a session only under the rules it was graded with.
    #[serde(default = "grading_v1")]
    pub grading_version: u32,
    /// SHA-256 of the sorted `k=v` lines of the child environment.
    pub child_env_sha256: String,
    pub child_env: BTreeMap<String, String>,
    pub started_unix_ms: u128,
    pub finished_unix_ms: Option<u128>,
    /// `running`, `complete` or `interrupted`.
    pub status: String,
    pub executions_planned: u64,
    pub records_written: u64,
    pub records_sha256: Option<String>,
    pub status_counts: BTreeMap<String, u64>,
}

/// The plan as written to `plan.json`.
#[derive(Serialize, Deserialize)]
pub struct PlanDoc {
    pub schema: String,
    pub spec_id: String,
    pub arms: Vec<crate::cryptanalysis::ecbench::spec::Arm>,
    pub workloads: Vec<Workload>,
    pub executions: Vec<Execution>,
}

fn now_ms() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

/// The child's environment: nothing inherited but `PATH` and `HOME`,
/// thread counts and locale pinned.
fn child_env() -> BTreeMap<String, String> {
    let mut env = BTreeMap::new();
    for k in ["PATH", "HOME"] {
        if let Ok(v) = std::env::var(k) {
            env.insert(k.to_string(), v);
        }
    }
    for (k, v) in [
        ("LC_ALL", "C"),
        ("RAYON_NUM_THREADS", "1"),
        ("OMP_NUM_THREADS", "1"),
        ("RUST_BACKTRACE", "0"),
        ("RUST_MIN_STACK", "8388608"),
    ] {
        env.insert(k.to_string(), v.to_string());
    }
    env
}

fn write_json(path: &Path, v: &impl Serialize) -> Result<String, String> {
    let text = serde_json::to_string_pretty(v).map_err(|e| e.to_string())? + "\n";
    std::fs::write(path, &text).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(sha256_hex(text.as_bytes()))
}

/// What the runner observed of one finished child.
struct Reaped {
    output: Option<ChildOutput>,
    parse_error: Option<String>,
    stderr: Vec<u8>,
    timing: Timing,
    exit: Option<String>,
    timed_out: bool,
}

/// Spawn, pin, feed, wait, reap.
fn run_child(
    opts: &RunOptions,
    plan: Option<&CpuPlan>,
    env: &BTreeMap<String, String>,
    input: &ChildInput,
    timeout: Duration,
) -> Result<Reaped, String> {
    let mut cmd = Command::new(&opts.exe);
    cmd.arg("exec")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .env_clear()
        .envs(env);
    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        // `ps` shows `ecbench exec`, not `/proc/self/exe exec`.
        cmd.arg0("ecbench");
        let pin = plan.map(|p| (p.run_cpu, p.node));
        let multi_node = topology().nodes.len() > 1;
        let parent = std::process::id();
        // SAFETY: only async-signal-safe system calls between fork and exec.
        unsafe {
            cmd.pre_exec(move || {
                if libc::setpgid(0, 0) != 0 {
                    return Err(std::io::Error::last_os_error());
                }
                signals::prepare_child(parent)?;
                #[cfg(target_os = "linux")]
                if let Some((cpu, node)) = pin {
                    crate::cryptanalysis::ecbench::isolation::pin_in_child(
                        cpu,
                        if multi_node { node } else { None },
                    )?;
                }
                #[cfg(not(target_os = "linux"))]
                let _ = (pin, multi_node);
                Ok(())
            });
        }
    }
    let payload = serde_json::to_vec(input).map_err(|e| e.to_string())?;
    let started = Instant::now();
    let mut child = cmd
        .spawn()
        .map_err(|e| format!("cannot exec {}: {e}", opts.exe.display()))?;
    let pid = child.id() as i32;
    signals::set_child(pid);
    {
        let mut stdin = child.stdin.take().ok_or("no child stdin")?;
        let _ = stdin.write_all(&payload);
    }
    let mut stdout = child.stdout.take().ok_or("no child stdout")?;
    let mut stderr = child.stderr.take().ok_or("no child stderr")?;
    let out_reader = std::thread::spawn(move || {
        let mut buf = Vec::new();
        let _ = stdout.read_to_end(&mut buf);
        buf
    });
    let err_reader = std::thread::spawn(move || {
        let mut buf = Vec::new();
        let _ = stderr.read_to_end(&mut buf);
        buf
    });
    // The watchdog sleeps until the deadline or until told the child is
    // done; it costs no CPU while the child runs.
    let done = Arc::new((Mutex::new(false), Condvar::new()));
    let fired = Arc::new(Mutex::new(false));
    let watchdog = {
        let done = Arc::clone(&done);
        let fired = Arc::clone(&fired);
        std::thread::spawn(move || {
            let (lock, cv) = &*done;
            let guard = lock.lock().unwrap_or_else(|e| e.into_inner());
            let (guard, _) = cv
                .wait_timeout_while(guard, timeout, |finished| !*finished)
                .unwrap_or_else(|e| e.into_inner());
            if !*guard {
                *fired.lock().unwrap_or_else(|e| e.into_inner()) = true;
                #[cfg(unix)]
                // SAFETY: signalling the child's own process group.
                unsafe {
                    libc::kill(-pid, libc::SIGKILL);
                }
            }
        })
    };
    let mut timing = Timing::default();
    let mut exit = None;
    #[cfg(unix)]
    {
        // SAFETY: waiting for our own child; status and rusage are filled.
        unsafe {
            let mut status: libc::c_int = 0;
            let mut ru: libc::rusage = std::mem::zeroed();
            let rc = loop {
                let rc = libc::wait4(pid, &mut status, 0, &mut ru);
                if rc == -1
                    && std::io::Error::last_os_error().kind() == std::io::ErrorKind::Interrupted
                {
                    continue;
                }
                break rc;
            };
            timing.process_wall_ns = started.elapsed().as_nanos() as u64;
            if rc == pid {
                let tv =
                    |t: libc::timeval| t.tv_sec as u64 * 1_000_000_000 + t.tv_usec as u64 * 1000;
                timing.user_ns = tv(ru.ru_utime);
                timing.sys_ns = tv(ru.ru_stime);
                // Linux reports kibibytes, macOS bytes.
                timing.max_rss_kib = if cfg!(target_os = "macos") {
                    ru.ru_maxrss as u64 / 1024
                } else {
                    ru.ru_maxrss as u64
                };
                timing.minor_faults = ru.ru_minflt as u64;
                timing.major_faults = ru.ru_majflt as u64;
                timing.voluntary_switches = ru.ru_nvcsw as u64;
                timing.involuntary_switches = ru.ru_nivcsw as u64;
                exit = Some(if libc::WIFEXITED(status) {
                    format!("exit {}", libc::WEXITSTATUS(status))
                } else if libc::WIFSIGNALED(status) {
                    format!("signal {}", libc::WTERMSIG(status))
                } else {
                    format!("status {status}")
                });
            }
        }
    }
    #[cfg(not(unix))]
    {
        let st = child.wait().map_err(|e| e.to_string())?;
        timing.process_wall_ns = started.elapsed().as_nanos() as u64;
        exit = Some(format!("{st}"));
    }
    signals::clear_child();
    {
        let (lock, cv) = &*done;
        *lock.lock().unwrap_or_else(|e| e.into_inner()) = true;
        cv.notify_all();
    }
    let _ = watchdog.join();
    // The child is reaped; std must not wait for it again.
    std::mem::forget(child);
    let out = out_reader.join().unwrap_or_default();
    let err = err_reader.join().unwrap_or_default();
    let timed_out = *fired.lock().unwrap_or_else(|e| e.into_inner());
    let (output, parse_error) = match serde_json::from_slice::<ChildOutput>(&out) {
        Ok(o) => (Some(o), None),
        Err(e) if !timed_out => (None, Some(format!("child output: {e}"))),
        Err(_) => (None, None),
    };
    if let Some(o) = &output {
        timing.schedstat_solve = o.schedstat_solve;
        timing.hw_solve = o.hw_solve.clone();
        timing.solve_wall_ns = o.report.as_ref().map(|r| r.solve_wall_ns);
    }
    Ok(Reaped {
        output,
        parse_error,
        stderr: err,
        timing,
        exit,
        timed_out,
    })
}

fn interrupted_error() -> Option<String> {
    signals::interrupted().map(|sig| format!("interrupted by signal {sig}"))
}

/// Size and modification time of the binary: a cheap check, where
/// `/proc/self/exe` is not available, that it was not rebuilt mid-session.
fn exe_fingerprint(exe: &Path) -> Option<(u64, std::time::SystemTime)> {
    let m = std::fs::metadata(exe).ok()?;
    Some((m.len(), m.modified().ok()?))
}

/// Run a whole session.  Refuses an existing output directory; refuses
/// to start on a busy host unless `allow_busy`.
pub fn run_session(
    spec_text: &str,
    out: &Path,
    opts: RunOptions,
    mut progress: impl FnMut(&str),
) -> Result<Session, String> {
    signals::install();
    // A fast refusal; the authoritative one is the atomic mkdir below.
    if out.exists() {
        return Err(format!(
            "{} exists; a session never overwrites another",
            out.display()
        ));
    }
    let spec = Spec::from_json(spec_text)?;
    let plan_: Plan = plan(spec)?;
    let capsule: HostCapsule = capture()?;
    let topo = capsule.stable.topology.clone();
    let allowed = allowed_cpus();
    let cpu_plan = plan_cpus(&opts.cpus, &topo, allowed.as_deref())?;

    // Lock, check and preflight before anything is written.
    let lock = BenchLock::acquire(&opts.lock_path, opts.wait_for_lock)?;
    let mut runner_on_run_cpu = false;
    let mut runner_cpus = None;
    // Restores on drop: a failure or panic from here on puts the host back.
    let mut eviction: Option<Eviction> = None;
    let mut pre = None;
    if let Some(p) = &cpu_plan {
        // The runner leaves the reservation.  Under an inherited
        // placement (an isolab cgroup) there may be nothing outside it;
        // the runner then takes an idle sibling, where it sleeps while
        // the child runs, and only as a last resort the run CPU itself.
        match leave_reserved(&p.reserved).or_else(|_| leave_reserved(&[p.run_cpu])) {
            Ok(cpus) => runner_cpus = Some(cpus),
            Err(_) => runner_on_run_cpu = true,
        }
        eviction = Some(evict(&p.reserved));
        let pf = preflight(p, &opts.thresholds);
        if let Some(e) = interrupted_error() {
            return Err(e);
        }
        if !pf.quiet && !opts.allow_busy {
            return Err(format!(
                "host is not quiet ({}); wait, or pass --allow-busy to record runs below L2",
                pf.reasons.join("; ")
            ));
        }
        pre = Some(pf);
    }
    // Claim the directory atomically, under the lock: a session queued
    // behind another for the same directory fails here instead of
    // writing over it.
    if let Some(parent) = out.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    std::fs::create_dir(out).map_err(|e| {
        if e.kind() == std::io::ErrorKind::AlreadyExists {
            format!(
                "{} exists; a session never overwrites another",
                out.display()
            )
        } else {
            format!("{}: {e}", out.display())
        }
    })?;
    let result = execute_session(
        out,
        &plan_,
        &capsule,
        cpu_plan.clone(),
        pre,
        eviction.as_ref().map(|e| e.summary.clone()),
        runner_on_run_cpu,
        runner_cpus,
        spec_text,
        opts,
        &mut progress,
    );
    // Put the host back before the session is declared finished either
    // way; the lock is released last.
    drop(eviction);
    let result = match result {
        Ok(mut session) => {
            session.status = "complete".into();
            session.finished_unix_ms = Some(now_ms());
            write_json(&out.join("session.json"), &session)?;
            Ok(session)
        }
        Err(e) => {
            // Mark what was written as interrupted; never delete it.
            if let Ok(text) = std::fs::read_to_string(out.join("session.json")) {
                if let Ok(mut s) = serde_json::from_str::<Session>(&text) {
                    s.status = "interrupted".into();
                    s.finished_unix_ms = Some(now_ms());
                    let _ = write_json(&out.join("session.json"), &s);
                }
            }
            Err(e)
        }
    };
    drop(lock);
    result
}

#[allow(clippy::too_many_arguments)]
fn execute_session(
    out: &Path,
    p: &Plan,
    capsule: &HostCapsule,
    cpu_plan: Option<CpuPlan>,
    pre: Option<Preflight>,
    eviction: Option<EvictionSummary>,
    runner_on_run_cpu: bool,
    runner_cpus: Option<Vec<u32>>,
    spec_text: &str,
    opts: RunOptions,
    progress: &mut impl FnMut(&str),
) -> Result<Session, String> {
    std::fs::create_dir(out.join("exec")).map_err(|e| format!("{}: {e}", out.display()))?;
    std::fs::write(out.join("spec.json"), spec_text).map_err(|e| e.to_string())?;
    let host_sha256 = write_json(&out.join("host.json"), capsule)?;
    let plan_doc = PlanDoc {
        schema: "ecbench.plan/v1".into(),
        spec_id: p.spec_id.clone(),
        arms: p.arms.clone(),
        workloads: p.workloads.clone(),
        executions: p.executions.clone(),
    };
    let plan_sha256 = write_json(&out.join("plan.json"), &plan_doc)?;
    let env = child_env();
    let env_lines: Vec<String> = env.iter().map(|(k, v)| format!("{k}={v}")).collect();
    let started = now_ms();
    let (session_id, _) = short_id(
        "ECBS1",
        &json!({
            "spec_sha256": p.spec_sha256,
            "host_sha256": host_sha256,
            "plan_sha256": plan_sha256,
            "started_unix_ms": started.to_string(),
        }),
    )?;
    let hz = ticks_per_second();
    let mut session = Session {
        schema: SESSION_SCHEMA.into(),
        session_id: session_id.clone(),
        label: p.spec.label.clone(),
        spec_id: p.spec_id.clone(),
        spec_sha256: p.spec_sha256.clone(),
        env_class_id: capsule.env_class_id.clone(),
        host_sha256,
        plan_sha256,
        binary_sha256: capsule.build.binary_sha256.clone(),
        git_commit: capsule.build.git_commit.clone(),
        git_dirty: capsule.build.git_dirty,
        options: opts.clone(),
        cpu_plan: cpu_plan.clone(),
        preflight: pre.clone(),
        eviction: eviction.clone(),
        runner_on_run_cpu,
        runner_cpus,
        ticks_per_second: hz,
        grading_version: GRADING_VERSION,
        child_env_sha256: sha256_hex(env_lines.join("\n").as_bytes()),
        child_env: env.clone(),
        started_unix_ms: started,
        finished_unix_ms: None,
        status: "running".into(),
        executions_planned: p.executions.len() as u64,
        records_written: 0,
        records_sha256: None,
        status_counts: BTreeMap::new(),
    };
    write_json(&out.join("session.json"), &session)?;

    let mut records = std::fs::OpenOptions::new()
        .create_new(true)
        .append(true)
        .open(out.join("records.jsonl"))
        .map_err(|e| e.to_string())?;
    let mut reps: BTreeMap<(String, String), u32> = BTreeMap::new();
    let timeout = Duration::from_secs(p.spec.measurement.timeout_seconds.max(1));
    let observed = cpu_plan
        .as_ref()
        .map(|c| c.observed_cpus())
        .unwrap_or_default();
    let exe_at_start = exe_fingerprint(&opts.exe);
    let result = (|| -> Result<(), String> {
        for ex in &p.executions {
            if let Some(e) = interrupted_error() {
                return Err(e);
            }
            if !opts.exe.starts_with("/proc") && exe_fingerprint(&opts.exe) != exe_at_start {
                return Err(format!(
                    "{} changed during the session; its records would not match the recorded binary hash",
                    opts.exe.display()
                ));
            }
            let arm = &p.arms[ex.arm];
            let w = &p.workloads[ex.workload];
            let inst = p.instance(ex.workload);
            let input = ChildInput {
                schema: CHILD_INPUT_SCHEMA.into(),
                curve: CurveSpec::explicit(inst).unwrap_or_else(|| w.curve_spec.clone()),
                target_seed: w.target_seed,
                target_index: w.target_index,
                expected_workload_id: w.workload_id.clone(),
                method: MethodSpec {
                    id: arm.method.id.clone(),
                    params: arm.method.params.clone(),
                },
                expected_method_id: arm.method.method_id.clone(),
                algorithm_seed: ex.algorithm_seed,
            };
            let before = conditions(&observed);
            let reaped = run_child(&opts, cpu_plan.as_ref(), &env, &input, timeout)?;
            let after = conditions(&observed);
            // n counts every execution of this (method, workload) in the
            // session, warm-ups and A/A arms included, so a run id is
            // unique in its session and n ≥ 1 as the convention requires.
            let rep_n = {
                let n = reps
                    .entry((arm.method.method_id.clone(), w.workload_id.clone()))
                    .or_insert(0);
                *n += 1;
                *n
            };
            let mut rec = assemble(
                &session_id,
                ex,
                rep_n,
                arm,
                w,
                inst,
                capsule,
                &reaped,
                before,
                after,
                &GradeContext {
                    plan: cpu_plan.as_ref(),
                    preflight: pre.as_ref(),
                    eviction: eviction.as_ref(),
                    runner_on_run_cpu,
                    thresholds: &opts.thresholds,
                    hz: hz as f64,
                },
            );
            let interrupted = interrupted_error();
            if let Some(e) = &interrupted {
                // The child was killed by the interruption: say so.
                rec.outcome.error = Some(match rec.outcome.error.take() {
                    Some(prev) => format!("{prev}; session {e}"),
                    None => format!("session {e}"),
                });
            }
            if !reaped.stderr.is_empty() {
                let _ = std::fs::write(
                    out.join("exec").join(format!("{:06}.stderr", ex.seq)),
                    &reaped.stderr,
                );
            }
            // Sealed over the exact bytes written (see `Record::seal`).
            let line = rec.seal() + "\n";
            records
                .write_all(line.as_bytes())
                .map_err(|e| e.to_string())?;
            records.sync_data().map_err(|e| e.to_string())?;
            *session
                .status_counts
                .entry(rec.outcome.status.clone())
                .or_insert(0) += 1;
            session.records_written += 1;
            progress(&format!(
                "{:>4}/{} {:<14} {} {} {:>9} S={} {}",
                ex.seq + 1,
                p.executions.len(),
                arm.name,
                w.curve.slug,
                w.workload_id,
                rec.outcome.status,
                rec.cost
                    .s
                    .map(|s| format!("{s:.3}"))
                    .unwrap_or_else(|| "-".into()),
                rec.isolation.level.name(),
            ));
            if let Some(e) = interrupted {
                return Err(e);
            }
        }
        Ok(())
    })();
    drop(records);
    // The counts and the hash so far, whatever happened: an interrupted
    // session keeps what it measured.
    if let Ok(body) = std::fs::read(out.join("records.jsonl")) {
        session.records_sha256 = Some(sha256_hex(&body));
    }
    write_json(&out.join("session.json"), &session)?;
    result.map(|()| session)
}

/// What grading needs from the session, beside the record's own
/// observations.
pub struct GradeContext<'a> {
    pub plan: Option<&'a CpuPlan>,
    pub preflight: Option<&'a Preflight>,
    pub eviction: Option<&'a EvictionSummary>,
    pub runner_on_run_cpu: bool,
    pub thresholds: &'a Thresholds,
    pub hz: f64,
}

/// Build the record for one execution from what the runner observed.
#[allow(clippy::too_many_arguments)]
fn assemble(
    session_id: &str,
    ex: &Execution,
    rep_n: u32,
    arm: &crate::cryptanalysis::ecbench::spec::Arm,
    w: &Workload,
    inst: &crate::cryptanalysis::ecbench::workload::Instance,
    capsule: &HostCapsule,
    reaped: &Reaped,
    before: crate::cryptanalysis::ecbench::isolation::Conditions,
    after: crate::cryptanalysis::ecbench::isolation::Conditions,
    ctx: &GradeContext,
) -> Record {
    let report = reaped.output.as_ref().and_then(|o| o.report.clone());
    let child_error = reaped
        .output
        .as_ref()
        .and_then(|o| o.error.clone())
        .or_else(|| reaped.parse_error.clone());
    let exited_cleanly = reaped.exit.as_deref() == Some("exit 0");
    let recovered = report.as_ref().and_then(|r| r.recovered);
    // Verification in the runner's own process: [k]G against the target.
    let matches_target = recovered.map(|k| inst.mul_generator_hex(k).as_ref() == Some(&w.target));
    let matches_planted = recovered.map(|k| k == w.planted);
    let status = if reaped.timed_out {
        "timeout"
    } else if !exited_cleanly && report.is_none() {
        "crashed"
    } else if child_error.is_some() || report.is_none() {
        "error"
    } else if matches_target == Some(true) && matches_planted == Some(true) {
        "verified"
    } else if recovered.is_some() {
        "wrong_answer"
    } else {
        "exhausted"
    };
    let sqrt_r = (w.curve.r as f64).sqrt();
    let floor = floor_s(w.curve.automorphisms_available);
    let s = report.as_ref().map(|r| r.total_gae / sqrt_r);
    let placement = reaped
        .output
        .as_ref()
        .map(|o| o.placement.clone())
        .unwrap_or_default();
    let g = grade(&GradeInput {
        plan: ctx.plan,
        preflight: ctx.preflight,
        eviction: ctx.eviction,
        runner_on_run_cpu: ctx.runner_on_run_cpu,
        capsule,
        thresholds: ctx.thresholds,
        placement: &placement,
        before: &before,
        after: &after,
        timing: &reaped.timing,
        hz: ctx.hz,
    });
    Record {
        schema: RECORD_SCHEMA.into(),
        record_id: String::new(),
        session_id: session_id.into(),
        run_id: format!("{}{}R{}", arm.method.method_id, w.workload_id, rep_n),
        seq: ex.seq,
        round: ex.round,
        warmup: ex.warmup,
        arm: arm.name.clone(),
        role: serde_json::to_value(arm.role)
            .ok()
            .and_then(|v| v.as_str().map(String::from))
            .unwrap_or_default(),
        method: arm.method.clone(),
        workload: w.clone(),
        algorithm_seed: ex.algorithm_seed.to_string(),
        outcome: Outcome {
            status: status.into(),
            recovered: recovered.map(|k| k.to_string()),
            matches_target,
            matches_planted,
            error: child_error,
            exit: reaped.exit.clone(),
            stderr_sha256: (!reaped.stderr.is_empty()).then(|| sha256_hex(&reaped.stderr)),
        },
        cost: Cost {
            unit: UNIT.into(),
            total_gae: report.as_ref().map(|r| r.total_gae),
            s,
            lower_bound: report
                .as_ref()
                .map(|r| !r.unpriced.is_empty())
                .unwrap_or(false),
            unpriced: report
                .as_ref()
                .map(|r| r.unpriced.clone())
                .unwrap_or_default(),
            deterministic: report.as_ref().map(|r| r.deterministic).unwrap_or(false),
            nondeterminism: report
                .as_ref()
                .map(|r| r.nondeterminism.clone())
                .unwrap_or_default(),
        },
        boundaries: Boundaries {
            automorphisms_available: w.curve.automorphisms_available,
            automorphisms_used: report.as_ref().map(|r| r.automorphisms_used),
            floor_s: floor,
            ratio_to_floor: s.map(|s| s / floor),
        },
        phases: report
            .as_ref()
            .map(|r| r.phases.clone())
            .unwrap_or_default(),
        counters: report
            .as_ref()
            .map(|r| r.counters.clone())
            .unwrap_or_default(),
        factor_base: report.as_ref().and_then(|r| r.factor_base.clone()),
        detail: report
            .as_ref()
            .map(|r| r.detail.clone())
            .unwrap_or(serde_json::Value::Null),
        online: report.as_ref().and_then(|r| r.online.clone()),
        time: reaped.timing.clone(),
        isolation: IsolationRecord {
            level: g.level,
            blockers: g.blockers,
            run_cpu: ctx.plan.map(|c| c.run_cpu),
            reserved: ctx.plan.map(|c| c.reserved.clone()).unwrap_or_default(),
            node: ctx.plan.and_then(|c| c.node),
            placement,
            foreign_busy_ticks: g.foreign_busy_ticks,
            steal_ticks: g.steal_ticks,
            idle_sibling_busy_ticks: g.idle_sibling_busy_ticks,
            before,
            after,
        },
        env_class_id: capsule.env_class_id.clone(),
        binary_sha256: capsule.build.binary_sha256.clone(),
    }
}

/// Read a session directory's records.
pub fn read_records(dir: &Path) -> Result<Vec<Record>, String> {
    Ok(read_record_lines(dir)?
        .into_iter()
        .map(|(_, r)| r)
        .collect())
}

/// Each record with the raw line it was read from, for the seal check.
pub fn read_record_lines(dir: &Path) -> Result<Vec<(String, Record)>, String> {
    let text = std::fs::read_to_string(dir.join("records.jsonl"))
        .map_err(|e| format!("{}: {e}", dir.join("records.jsonl").display()))?;
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .enumerate()
        .map(|(i, l)| {
            serde_json::from_str(l)
                .map(|r| (l.to_string(), r))
                .map_err(|e| format!("records.jsonl line {}: {e}", i + 1))
        })
        .collect()
}

pub fn read_session(dir: &Path) -> Result<Session, String> {
    let text = std::fs::read_to_string(dir.join("session.json"))
        .map_err(|e| format!("{}: {e}", dir.join("session.json").display()))?;
    serde_json::from_str(&text).map_err(|e| format!("session.json: {e}"))
}
