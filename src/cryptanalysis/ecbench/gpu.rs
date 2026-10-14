//! A bounded, provenance-preserving front door to `gpu/ecc/bench`.
//!
//! Stage benchmarks remain separate from the matched-workload GPU receipt.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use serde::{Deserialize, Serialize};

use super::canonical::sha256_hex;
use super::record::Record;
use super::workload::{Instance, Workload};

pub const GPU_STAGE_SCHEMA: &str = "ecbench.gpu_stage/v1";
pub const GPU_MATCH_SCHEMA: &str = "ecbench.gpu_match/v1";

#[derive(Clone, Debug)]
pub struct GpuOptions {
    pub binary: PathBuf,
    pub out: PathBuf,
    pub mode: String,
    pub args: Vec<String>,
    pub timeout_seconds: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CommandReceipt {
    pub argv: Vec<String>,
    pub started_unix_ms: u128,
    pub elapsed_ms: u128,
    #[serde(default)]
    pub elapsed_ns: u128,
    #[serde(default)]
    pub child_pid: u32,
    pub exit_code: Option<i32>,
    pub timed_out: bool,
    pub stdout_sha256: String,
    pub stderr_sha256: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GpuStage {
    pub schema: String,
    pub kind: String,
    pub binary: String,
    pub binary_sha256: String,
    pub cuda_visible_devices: Option<String>,
    pub mode: String,
    pub args: Vec<String>,
    pub timeout_seconds: u64,
    pub selftest: CommandReceipt,
    pub benchmark: Option<CommandReceipt>,
    pub benchmark_output_recognized: Option<bool>,
    pub status: String,
    pub limitation: String,
}

fn unix_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn run_one(
    binary: &Path,
    args: &[String],
    out: &Path,
    stem: &str,
    timeout: u64,
) -> Result<(CommandReceipt, String), String> {
    let mut cmd = Command::new(binary);
    cmd.args(args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        // Keep children of the CUDA driver in one group for the watchdog.
        unsafe {
            cmd.pre_exec(|| {
                if libc::setpgid(0, 0) == 0 {
                    Ok(())
                } else {
                    Err(std::io::Error::last_os_error())
                }
            });
        }
    }
    let started_unix_ms = unix_ms();
    let start = Instant::now();
    let child = match cmd.spawn() {
        Ok(child) => child,
        Err(e) => {
            let message = format!("{}: {e}", binary.display());
            std::fs::write(out.join(format!("{stem}.stdout")), b"").map_err(|e| e.to_string())?;
            std::fs::write(out.join(format!("{stem}.stderr")), message.as_bytes())
                .map_err(|e| e.to_string())?;
            return Ok((
                CommandReceipt {
                    argv: std::iter::once(binary.display().to_string())
                        .chain(args.iter().cloned())
                        .collect(),
                    started_unix_ms,
                    elapsed_ms: start.elapsed().as_millis(),
                    elapsed_ns: start.elapsed().as_nanos(),
                    child_pid: 0,
                    exit_code: None,
                    timed_out: false,
                    stdout_sha256: sha256_hex(b""),
                    stderr_sha256: sha256_hex(message.as_bytes()),
                },
                String::new(),
            ));
        }
    };
    let pid = child.id() as i32;
    let done = Arc::new((Mutex::new(false), Condvar::new()));
    let timed_out = Arc::new(Mutex::new(false));
    let watchdog = {
        let done = Arc::clone(&done);
        let timed_out = Arc::clone(&timed_out);
        std::thread::spawn(move || {
            let (lock, cv) = &*done;
            let guard = lock.lock().unwrap_or_else(|e| e.into_inner());
            let (guard, _) = cv
                .wait_timeout_while(guard, Duration::from_secs(timeout), |finished| !*finished)
                .unwrap_or_else(|e| e.into_inner());
            if !*guard {
                *timed_out.lock().unwrap_or_else(|e| e.into_inner()) = true;
                #[cfg(unix)]
                unsafe {
                    libc::kill(-pid, libc::SIGKILL);
                }
                #[cfg(not(unix))]
                let _ = pid;
            }
        })
    };
    let output = child.wait_with_output().map_err(|e| e.to_string());
    {
        let (lock, cv) = &*done;
        *lock.lock().unwrap_or_else(|e| e.into_inner()) = true;
        cv.notify_one();
    }
    watchdog
        .join()
        .map_err(|_| "GPU watchdog panicked".to_string())?;
    let output = output?;
    std::fs::write(out.join(format!("{stem}.stdout")), &output.stdout)
        .map_err(|e| e.to_string())?;
    std::fs::write(out.join(format!("{stem}.stderr")), &output.stderr)
        .map_err(|e| e.to_string())?;
    let receipt = CommandReceipt {
        argv: std::iter::once(binary.display().to_string())
            .chain(args.iter().cloned())
            .collect(),
        started_unix_ms,
        elapsed_ms: start.elapsed().as_millis(),
        elapsed_ns: start.elapsed().as_nanos(),
        child_pid: pid as u32,
        exit_code: output.status.code(),
        timed_out: *timed_out.lock().unwrap_or_else(|e| e.into_inner()),
        stdout_sha256: sha256_hex(&output.stdout),
        stderr_sha256: sha256_hex(&output.stderr),
    };
    Ok((
        receipt,
        String::from_utf8_lossy(&output.stdout).into_owned(),
    ))
}

/// Run the CUDA binary's device self-test, then one bounded stage benchmark.
/// The output directory is new and never overwritten, even on failure.
pub fn run(opts: &GpuOptions) -> Result<GpuStage, String> {
    #[cfg(not(unix))]
    return Err("the GPU stage watchdog currently requires a Unix process group".into());
    if !["field", "mul", "rho", "bsgs"].contains(&opts.mode.as_str()) {
        return Err("GPU mode must be field, mul, rho, or bsgs".into());
    }
    if opts.timeout_seconds == 0 {
        return Err("GPU timeout must be positive".into());
    }
    if opts.args.iter().any(|a| a == "--solve") {
        return Err("--solve is excluded from this bounded stage path; use the GPU solver directly under its own experiment protocol".into());
    }
    let binary = std::fs::canonicalize(&opts.binary)
        .map_err(|e| format!("{}: {e}", opts.binary.display()))?;
    let bytes = std::fs::read(&binary).map_err(|e| format!("{}: {e}", binary.display()))?;
    std::fs::create_dir(&opts.out)
        .map_err(|e| format!("{}: {e}; use a new output directory", opts.out.display()))?;
    // Execute a private byte copy so a concurrent rebuild cannot change
    // the binary between the device self-test and the measured stage.
    let snapshot = opts.out.join("gpu-bench.binary");
    std::fs::copy(&binary, &snapshot).map_err(|e| format!("{}: {e}", snapshot.display()))?;
    if sha256_hex(&std::fs::read(&snapshot).map_err(|e| e.to_string())?) != sha256_hex(&bytes) {
        return Err("the CUDA binary changed while being copied".into());
    }
    let (selftest, output) = run_one(
        &snapshot,
        &["selftest".into()],
        &opts.out,
        "selftest",
        opts.timeout_seconds,
    )?;
    let ok = selftest.exit_code == Some(0) && !selftest.timed_out && output.contains("selftest OK");
    let (benchmark, benchmark_output_recognized) = if ok {
        let args = std::iter::once(opts.mode.clone())
            .chain(opts.args.iter().cloned())
            .collect::<Vec<_>>();
        let (receipt, stdout) = run_one(
            &snapshot,
            &args,
            &opts.out,
            "benchmark",
            opts.timeout_seconds,
        )?;
        let marker = match opts.mode.as_str() {
            "field" => "mul (throughput):",
            "mul" => "variable-base scalar mul:",
            "rho" => "Gstep/s",
            "bsgs" => "targets solved and verified",
            _ => unreachable!(),
        };
        (Some(receipt), Some(stdout.contains(marker)))
    } else {
        (None, None)
    };
    let status = if !ok {
        "selftest_failed"
    } else if benchmark
        .as_ref()
        .is_some_and(|r| r.exit_code == Some(0) && !r.timed_out)
        && benchmark_output_recognized == Some(true)
    {
        "stage_completed"
    } else {
        "stage_failed"
    };
    let receipt = GpuStage {
        schema: GPU_STAGE_SCHEMA.into(),
        kind: "gpu_ecc_stage".into(),
        binary: binary.display().to_string(),
        binary_sha256: sha256_hex(&bytes),
        cuda_visible_devices: std::env::var("CUDA_VISIBLE_DEVICES").ok(),
        mode: opts.mode.clone(),
        args: opts.args.clone(),
        timeout_seconds: opts.timeout_seconds,
        selftest,
        benchmark,
        benchmark_output_recognized,
        status: status.into(),
        limitation: "Fixed-curve GPU stage output: no shared EC Bench workload identity, external answer replay, matched CPU baseline, GPU contention sampling, or admitted end-to-end comparison".into(),
    };
    let body = serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())? + "\n";
    std::fs::write(opts.out.join("gpu-stage.json"), body).map_err(|e| e.to_string())?;
    Ok(receipt)
}

#[derive(Clone, Debug)]
pub struct GpuMatchOptions {
    pub binary: PathBuf,
    pub workload: PathBuf,
    pub out: PathBuf,
    pub timeout_seconds: u64,
    pub telemetry_binary: PathBuf,
    pub cpu_record: Option<PathBuf>,
    pub args: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GpuSample {
    pub unix_ms: u128,
    pub completed_unix_ms: u128,
    pub gpu_csv: Option<String>,
    pub pcie_traffic_csv: Option<String>,
    #[serde(default)]
    pub pcie_error: Option<String>,
    pub processes_csv: Option<String>,
    pub error: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GpuMatch {
    pub schema: String,
    pub workload_id: String,
    pub workload_sha256: String,
    pub binary_sha256: String,
    pub cuda_visible_devices: Option<String>,
    pub cpu_record_id: Option<String>,
    pub cpu_record_sha256: Option<String>,
    pub selftest: CommandReceipt,
    pub solve: Option<CommandReceipt>,
    pub recovered: Option<String>,
    pub independently_verified: bool,
    pub baby_steps: Option<u64>,
    pub giant_steps: Option<u64>,
    pub seed_ops: Option<u64>,
    pub device_seconds: Option<f64>,
    pub host_seconds: Option<f64>,
    pub setup_seconds: Option<f64>,
    pub candidate_verification_seconds: Option<f64>,
    pub transfer_seconds: Option<f64>,
    pub h2d_bytes: Option<u64>,
    pub d2h_bytes: Option<u64>,
    pub d2d_bytes: Option<u64>,
    pub telemetry: Vec<GpuSample>,
    pub telemetry_binary: String,
    pub wall_admission: String,
    pub wall_blockers: Vec<String>,
    pub status: String,
    pub cost_scope: String,
}

fn sample_gpu(binary: &Path) -> GpuSample {
    let started = unix_ms();
    let query = |arg: &str| -> Result<String, String> {
        let o = Command::new(binary)
            .args([arg, "--format=csv,noheader,nounits"])
            .stdin(Stdio::null())
            .output()
            .map_err(|e| e.to_string())?;
        if !o.status.success() {
            let detail = if o.stderr.is_empty() {
                &o.stdout
            } else {
                &o.stderr
            };
            return Err(format!(
                "{}: {}",
                arg,
                String::from_utf8_lossy(detail).trim()
            ));
        }
        String::from_utf8(o.stdout).map_err(|e| e.to_string())
    };
    let gpu = query("--query-gpu=index,uuid,clocks.current.sm,clocks.current.memory,utilization.gpu,utilization.memory,pcie.link.gen.current,pcie.link.width.current,memory.used,memory.total");
    let pcie = query("--query-gpu=pcie.rx_util,pcie.tx_util");
    let processes = query("--query-compute-apps=pid,gpu_uuid");
    GpuSample {
        unix_ms: started,
        completed_unix_ms: unix_ms(),
        gpu_csv: gpu.as_ref().ok().cloned(),
        pcie_traffic_csv: pcie.as_ref().ok().cloned(),
        pcie_error: pcie.as_ref().err().cloned(),
        processes_csv: processes.as_ref().ok().cloned(),
        error: gpu.err().or_else(|| processes.as_ref().err().cloned()),
    }
}

fn fields(line: &str, marker: &str) -> Option<std::collections::BTreeMap<String, String>> {
    let mut found = None;
    for row in line.lines().filter(|row| row.starts_with(marker)) {
        if found.is_some() {
            return None;
        }
        let pairs = row.split_whitespace().skip(1);
        let mut map = std::collections::BTreeMap::new();
        for pair in pairs {
            let (k, v) = pair.split_once('=')?;
            if map.insert(k.to_string(), v.to_string()).is_some() {
                return None;
            }
        }
        found = Some(map);
    }
    found
}

fn field_u64(map: &std::collections::BTreeMap<String, String>, key: &str) -> Option<u64> {
    map.get(key)?.parse().ok()
}

fn parse_scalar_hex(value: &str) -> Option<u128> {
    let digits = value.strip_prefix("0x")?.trim_start_matches('0');
    u128::from_str_radix(if digits.is_empty() { "0" } else { digits }, 16).ok()
}

fn parse_coordinate_hex(value: &str) -> Result<u64, String> {
    let digits = value
        .strip_prefix("0x")
        .or_else(|| value.strip_prefix("0X"))
        .unwrap_or(value);
    u64::from_str_radix(digits, 16)
        .map_err(|e| format!("invalid workload target coordinate {value:?}: {e}"))
}

fn gpu_admission(
    samples: &[GpuSample],
    child_pid: u32,
    started: u128,
    elapsed: u128,
) -> Vec<String> {
    let mut blockers = Vec::new();
    if samples.len() < 2 {
        blockers.push("fewer than two device samples".into());
    }
    let ended = started.saturating_add(elapsed);
    if let Some(before) = samples.first().and_then(|s| s.gpu_csv.as_ref()) {
        if let Some(row) = before.lines().next() {
            let cols: Vec<_> = row.split(',').map(str::trim).collect();
            if cols
                .get(4)
                .and_then(|v| v.parse::<f64>().ok())
                .is_some_and(|v| v > 20.0)
                || cols
                    .get(5)
                    .and_then(|v| v.parse::<f64>().ok())
                    .is_some_and(|v| v > 20.0)
            {
                blockers.push("GPU activity exceeded 20 percent before solve".into());
            }
            if let (Some(used), Some(total)) = (
                cols.get(8).and_then(|v| v.parse::<f64>().ok()),
                cols.get(9).and_then(|v| v.parse::<f64>().ok()),
            ) {
                if total > 0.0 && used / total > 0.90 {
                    blockers.push("GPU memory exceeded 90 percent before solve".into());
                }
            }
        }
    }
    let mut clocks = Vec::new();
    let mut identities = std::collections::BTreeSet::new();
    for s in samples {
        if let Some(e) = &s.error {
            blockers.push(format!("telemetry failed: {e}"));
            continue;
        }
        let Some(gpu) = &s.gpu_csv else {
            blockers.push("device telemetry missing".into());
            continue;
        };
        let rows: Vec<_> = gpu.lines().collect();
        if rows.len() != 1 {
            blockers.push("visible GPU is ambiguous".into());
            continue;
        }
        let cols: Vec<_> = rows[0].split(',').map(str::trim).collect();
        if cols.len() != 10 || cols.iter().any(|c| c.is_empty() || *c == "[Not Supported]") {
            blockers.push("device clocks, activity, or PCIe state unavailable".into());
            continue;
        }
        identities.insert(cols[1].to_string());
        if s.pcie_traffic_csv.as_ref().is_none_or(|v| {
            let parts: Vec<_> = v.trim().split(',').map(str::trim).collect();
            parts.len() != 2 || parts.iter().any(|p| p.parse::<f64>().is_err())
        }) {
            blockers.push(match &s.pcie_error {
                Some(e) => format!("PCIe traffic counters unavailable: {e}"),
                None => "PCIe traffic counters unavailable".into(),
            });
        }
        match cols[2].parse::<f64>() {
            Ok(v) if v > 0.0 => {
                if s.unix_ms <= ended && s.completed_unix_ms >= started {
                    clocks.push(v);
                }
            }
            _ => blockers.push("SM clock is invalid".into()),
        }
        if cols[3].parse::<f64>().is_err()
            || cols[4].parse::<f64>().is_err()
            || cols[5].parse::<f64>().is_err()
            || cols[6].parse::<u32>().is_err()
            || cols[7].parse::<u32>().is_err()
            || !matches!((cols[8].parse::<f64>(), cols[9].parse::<f64>()),
                (Ok(used), Ok(total)) if used >= 0.0 && total > 0.0 && used <= total)
        {
            blockers
                .push("device memory capacity, clock, activity, or PCIe link is invalid".into());
        }
        let Some(apps) = &s.processes_csv else {
            blockers.push("GPU process list unavailable".into());
            continue;
        };
        for row in apps.lines() {
            if row.trim().is_empty() || row.contains("No running processes") {
                continue;
            }
            match row.split(',').next().unwrap_or("").trim().parse::<u32>() {
                Ok(pid) if pid == child_pid => {}
                Ok(pid) => blockers.push(format!("competing GPU process {pid}")),
                Err(_) => blockers.push("GPU process list could not be parsed".into()),
            }
        }
    }
    if identities.len() != 1 {
        blockers.push("GPU identity changed or is unavailable".into());
    }
    if clocks.len() < 2 {
        blockers.push("fewer than two in-run SM clock samples".into());
    }
    if let (Some(min), Some(max)) = (
        clocks.iter().copied().reduce(f64::min),
        clocks.iter().copied().reduce(f64::max),
    ) {
        if max > min * 1.10 {
            blockers.push("SM clock varied by more than 10 percent".into());
        }
    }
    blockers.sort();
    blockers.dedup();
    blockers
}

/// Run the existing CUDA BSGS path on a frozen EC Bench point. The CUDA
/// process never receives the planted scalar. The result is checked with
/// EC Bench's independent group implementation and all raw output is kept.
pub fn run_match(opts: &GpuMatchOptions) -> Result<GpuMatch, String> {
    if opts.timeout_seconds == 0 {
        return Err("GPU timeout must be positive".into());
    }
    let raw_workload = std::fs::read(&opts.workload).map_err(|e| e.to_string())?;
    let workload: Workload = serde_json::from_slice(&raw_workload).map_err(|e| e.to_string())?;
    if workload.recompute_id()? != workload.workload_id {
        return Err("workload identity does not reproduce".into());
    }
    let (rebuilt, inst) = Workload::build(
        &workload.curve_spec,
        workload.target_seed,
        workload.target_index,
        workload.kind(),
    )?;
    if rebuilt != workload {
        return Err("workload does not rebuild from its declared seed and curve".into());
    }
    let Instance::Prime(prime) = &inst else {
        return Err("CUDA BSGS currently accepts one-word prime curves".into());
    };
    if prime.r > (1u64 << 63) - 1 {
        return Err("CUDA BSGS interval exceeds 63 bits".into());
    }
    if opts.args.len() % 2 != 0
        || opts.args.chunks_exact(2).any(|pair| {
            !["--walks", "--iters", "--w", "--neg", "--variant"].contains(&pair[0].as_str())
        })
    {
        return Err(
            "GPU tuning args must be value pairs: --walks, --iters, --w, --neg, or --variant"
                .into(),
        );
    }
    let target = workload.target.clone();
    let x = parse_coordinate_hex(&target[0])?;
    let y = parse_coordinate_hex(&target[1])?;
    let mut cpu_id = None;
    let mut cpu_sha = None;
    let mut cpu_raw = None;
    if let Some(path) = &opts.cpu_record {
        let bytes = std::fs::read(path).map_err(|e| e.to_string())?;
        let rec: Record = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
        if !super::record::check_line_seal(
            std::str::from_utf8(&bytes)
                .map_err(|e| e.to_string())?
                .trim_end(),
            &rec.record_id,
        ) {
            return Err("CPU record line seal is invalid".into());
        }
        if rec.workload.workload_id != workload.workload_id || !rec.counts() {
            return Err("CPU record must be a verified measurement of this workload".into());
        }
        cpu_id = Some(rec.record_id);
        cpu_sha = Some(sha256_hex(&bytes));
        cpu_raw = Some(bytes);
    }
    let binary = std::fs::canonicalize(&opts.binary).map_err(|e| e.to_string())?;
    let bytes = std::fs::read(&binary).map_err(|e| e.to_string())?;
    std::fs::create_dir(&opts.out)
        .map_err(|e| format!("{}: {e}; use a new output directory", opts.out.display()))?;
    std::fs::write(opts.out.join("workload.json"), &raw_workload).map_err(|e| e.to_string())?;
    if let Some(raw) = cpu_raw {
        std::fs::write(opts.out.join("cpu-record.json"), raw).map_err(|e| e.to_string())?;
    }
    let snapshot = opts.out.join("gpu-bench.binary");
    std::fs::copy(&binary, &snapshot).map_err(|e| e.to_string())?;
    if sha256_hex(&std::fs::read(&snapshot).map_err(|e| e.to_string())?) != sha256_hex(&bytes) {
        return Err("CUDA binary changed while being copied".into());
    }
    let (selftest, output) = run_one(
        &snapshot,
        &["selftest".into()],
        &opts.out,
        "selftest",
        opts.timeout_seconds,
    )?;
    let mut receipt = GpuMatch {
        schema: GPU_MATCH_SCHEMA.into(),
        workload_id: workload.workload_id.clone(),
        workload_sha256: sha256_hex(&raw_workload),
        binary_sha256: sha256_hex(&bytes),
        cuda_visible_devices: std::env::var("CUDA_VISIBLE_DEVICES").ok(),
        cpu_record_id: cpu_id,
        cpu_record_sha256: cpu_sha,
        selftest,
        solve: None,
        recovered: None,
        independently_verified: false,
        baby_steps: None,
        giant_steps: None,
        seed_ops: None,
        device_seconds: None,
        host_seconds: None,
        setup_seconds: None,
        candidate_verification_seconds: None,
        transfer_seconds: None,
        h2d_bytes: None,
        d2h_bytes: None,
        d2d_bytes: None,
        telemetry: Vec::new(),
        telemetry_binary: opts.telemetry_binary.display().to_string(),
        wall_admission: "refused".into(),
        wall_blockers: Vec::new(),
        status: "selftest_failed".into(),
        cost_scope: "solve.elapsed_ns includes process startup, allocation, transfers, table build, device kernels, candidate recovery and teardown. host_seconds runs from BSGS entry through the answer, excluding process startup and teardown. CUDA event time covers measured device phases only. GPU steps are not EC Bench GAE.".into(),
    };
    if receipt.selftest.exit_code == Some(0)
        && !receipt.selftest.timed_out
        && output.contains("selftest OK")
    {
        // Keep the self-test's rolling utilization window out of preflight.
        std::thread::sleep(Duration::from_secs(2));
        let stop = Arc::new(AtomicBool::new(false));
        let samples = Arc::new(Mutex::new(vec![sample_gpu(&opts.telemetry_binary)]));
        let monitor = {
            let stop = Arc::clone(&stop);
            let samples = Arc::clone(&samples);
            let telemetry_binary = opts.telemetry_binary.clone();
            std::thread::spawn(move || {
                while !stop.load(Ordering::Relaxed) {
                    let sample = sample_gpu(&telemetry_binary);
                    samples
                        .lock()
                        .unwrap_or_else(|e| e.into_inner())
                        .push(sample);
                    std::thread::sleep(Duration::from_millis(200));
                }
            })
        };
        let args = std::iter::once("bsgs".to_string())
            .chain([
                "--ecbench-target-x".into(),
                format!("0x{x:x}"),
                "--ecbench-target-y".into(),
                format!("0x{y:x}"),
                "--targets".into(),
                "1".into(),
            ])
            .chain(opts.args.iter().cloned())
            .collect::<Vec<String>>();
        let result = run_one(&snapshot, &args, &opts.out, "solve", opts.timeout_seconds);
        stop.store(true, Ordering::Relaxed);
        monitor
            .join()
            .map_err(|_| "GPU telemetry thread panicked".to_string())?;
        samples
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .push(sample_gpu(&opts.telemetry_binary));
        receipt.telemetry = samples.lock().unwrap_or_else(|e| e.into_inner()).clone();
        let (solve, stdout) = result?;
        receipt.wall_blockers = gpu_admission(
            &receipt.telemetry,
            solve.child_pid,
            solve.started_unix_ms,
            solve.elapsed_ms,
        );
        if receipt.wall_blockers.is_empty() {
            receipt.wall_admission = "admitted".into();
        }
        let curve = fields(&stdout, "ECBENCH_CURVE");
        let curve_matches = curve.as_ref().is_some_and(|v| {
            field_u64(v, "p") == Some(prime.curve.p)
                && field_u64(v, "a") == Some(prime.curve.a)
                && field_u64(v, "b") == Some(prime.curve.b)
                && field_u64(v, "r") == Some(prime.r)
                && field_u64(v, "gx") == Some(prime.generator.0)
                && field_u64(v, "gy") == Some(prime.generator.1)
                && field_u64(v, "cofactor") == Some(prime.cofactor)
        });
        let result = fields(&stdout, "ECBENCH_RESULT");
        if let Some(v) = &result {
            receipt.baby_steps = field_u64(v, "baby_steps");
            receipt.giant_steps = field_u64(v, "giant_steps");
            receipt.seed_ops = field_u64(v, "seed_ops");
            receipt.device_seconds = v.get("device_seconds").and_then(|s| s.parse().ok());
            receipt.host_seconds = v.get("host_seconds").and_then(|s| s.parse().ok());
            receipt.setup_seconds = v.get("setup_seconds").and_then(|s| s.parse().ok());
            receipt.candidate_verification_seconds = v
                .get("candidate_verification_seconds")
                .and_then(|s| s.parse().ok());
            receipt.transfer_seconds = v.get("transfer_seconds").and_then(|s| s.parse().ok());
            receipt.h2d_bytes = field_u64(v, "h2d_bytes");
            receipt.d2h_bytes = field_u64(v, "d2h_bytes");
            receipt.d2d_bytes = field_u64(v, "d2d_bytes");
            receipt.recovered = v.get("k").cloned();
        }
        let k = receipt.recovered.as_deref().and_then(parse_scalar_hex);
        receipt.independently_verified = k.is_some_and(|k| {
            k < u128::from(prime.r) && inst.mul_generator_hex(k).as_ref() == Some(&target)
        });
        receipt.status = if solve.exit_code != Some(0) || solve.timed_out {
            "solve_failed"
        } else if !curve_matches {
            "curve_mismatch"
        } else if !receipt.independently_verified
            || receipt.baby_steps.is_none()
            || receipt.giant_steps.is_none()
            || receipt.seed_ops.is_none()
            || !receipt
                .device_seconds
                .is_some_and(|v| v.is_finite() && v > 0.0)
            || !receipt
                .host_seconds
                .is_some_and(|v| v.is_finite() && v > 0.0)
            || !receipt
                .setup_seconds
                .is_some_and(|v| v.is_finite() && v >= 0.0)
            || !receipt
                .candidate_verification_seconds
                .is_some_and(|v| v.is_finite() && v >= 0.0)
            || !receipt
                .transfer_seconds
                .is_some_and(|v| v.is_finite() && v >= 0.0)
            || receipt.h2d_bytes.is_none()
            || receipt.d2h_bytes.is_none()
            || receipt.d2d_bytes.is_none()
        {
            "unverified"
        } else {
            "verified"
        }
        .into();
        if receipt.status != "verified" {
            receipt.wall_admission = "refused".into();
            receipt.wall_blockers.push("solve not verified".into());
            receipt.wall_blockers.sort();
            receipt.wall_blockers.dedup();
        }
        receipt.solve = Some(solve);
    } else {
        receipt.wall_blockers.push("device selftest failed".into());
    }
    std::fs::write(
        opts.out.join("gpu-match.json"),
        serde_json::to_vec_pretty(&receipt).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    Ok(receipt)
}

/// Re-derive the workload, binary and raw-stream bindings and replay the
/// answer with EC Bench arithmetic. This reads the frozen directory only.
pub fn audit_match(dir: &Path) -> Result<Vec<String>, String> {
    let raw = std::fs::read(dir.join("gpu-match.json")).map_err(|e| e.to_string())?;
    let receipt: GpuMatch = serde_json::from_slice(&raw).map_err(|e| e.to_string())?;
    let mut issues = Vec::new();
    if receipt.schema != GPU_MATCH_SCHEMA {
        issues.push("receipt schema mismatch".into());
    }
    let workload_raw = std::fs::read(dir.join("workload.json")).map_err(|e| e.to_string())?;
    let workload: Workload = serde_json::from_slice(&workload_raw).map_err(|e| e.to_string())?;
    if sha256_hex(&workload_raw) != receipt.workload_sha256
        || workload.workload_id != receipt.workload_id
        || workload.recompute_id()? != workload.workload_id
    {
        issues.push("workload binding mismatch".into());
    }
    let (rebuilt, inst) = Workload::build(
        &workload.curve_spec,
        workload.target_seed,
        workload.target_index,
        workload.kind(),
    )?;
    if rebuilt != workload {
        issues.push("workload reconstruction mismatch".into());
    }
    let binary = std::fs::read(dir.join("gpu-bench.binary")).map_err(|e| e.to_string())?;
    if sha256_hex(&binary) != receipt.binary_sha256 {
        issues.push("binary snapshot hash mismatch".into());
    }
    let check_command = |name: &str,
                         command: &CommandReceipt,
                         issues: &mut Vec<String>|
     -> Result<String, String> {
        let stdout =
            std::fs::read(dir.join(format!("{name}.stdout"))).map_err(|e| e.to_string())?;
        let stderr =
            std::fs::read(dir.join(format!("{name}.stderr"))).map_err(|e| e.to_string())?;
        if sha256_hex(&stdout) != command.stdout_sha256
            || sha256_hex(&stderr) != command.stderr_sha256
        {
            issues.push(format!("{name} raw stream hash mismatch"));
        }
        String::from_utf8(stdout).map_err(|e| e.to_string())
    };
    let selftest = check_command("selftest", &receipt.selftest, &mut issues)?;
    if receipt.selftest.argv.get(1).map(String::as_str) != Some("selftest")
        || receipt
            .selftest
            .argv
            .first()
            .and_then(|s| Path::new(s).file_name())
            != Some(std::ffi::OsStr::new("gpu-bench.binary"))
    {
        issues.push("selftest command mismatch".into());
    }
    let selftest_ok = receipt.selftest.exit_code == Some(0)
        && !receipt.selftest.timed_out
        && selftest.contains("selftest OK");
    if selftest_ok == (receipt.status == "selftest_failed") {
        issues.push("selftest status does not reproduce".into());
    }
    if let Some(command) = &receipt.solve {
        let stdout = check_command("solve", command, &mut issues)?;
        let expected_x = format!("0x{:x}", parse_coordinate_hex(&workload.target[0])?);
        let expected_y = format!("0x{:x}", parse_coordinate_hex(&workload.target[1])?);
        let argv = &command.argv;
        if argv.first().and_then(|s| Path::new(s).file_name())
            != Some(std::ffi::OsStr::new("gpu-bench.binary"))
            || argv.get(1).map(String::as_str) != Some("bsgs")
            || argv.get(2).map(String::as_str) != Some("--ecbench-target-x")
            || argv.get(3).map(String::as_str) != Some(expected_x.as_str())
            || argv.get(4).map(String::as_str) != Some("--ecbench-target-y")
            || argv.get(5).map(String::as_str) != Some(expected_y.as_str())
        {
            issues.push("solve command is not bound to the workload target".into());
        }
        let curve = fields(&stdout, "ECBENCH_CURVE");
        let curve_ok = if let (Some(v), Instance::Prime(prime)) = (&curve, &inst) {
            field_u64(v, "p") == Some(prime.curve.p)
                && field_u64(v, "a") == Some(prime.curve.a)
                && field_u64(v, "b") == Some(prime.curve.b)
                && field_u64(v, "r") == Some(prime.r)
                && field_u64(v, "gx") == Some(prime.generator.0)
                && field_u64(v, "gy") == Some(prime.generator.1)
                && field_u64(v, "cofactor") == Some(prime.cofactor)
        } else {
            false
        };
        if !curve_ok && receipt.status == "verified" {
            issues.push("compiled curve does not match workload".into());
        }
        let result = fields(&stdout, "ECBENCH_RESULT");
        let k = result
            .as_ref()
            .and_then(|v| v.get("k"))
            .and_then(|v| parse_scalar_hex(v));
        let verified = k.is_some_and(|k| {
            k < workload.curve.r && inst.mul_generator_hex(k).as_ref() == Some(&workload.target)
        });
        if verified != receipt.independently_verified || (!verified && receipt.status == "verified")
        {
            issues.push("independent scalar replay failed".into());
        }
        if result.as_ref().and_then(|v| v.get("k")) != receipt.recovered.as_ref() {
            issues.push("recovered scalar differs from raw output".into());
        }
        if result.as_ref().and_then(|v| field_u64(v, "baby_steps")) != receipt.baby_steps
            || result.as_ref().and_then(|v| field_u64(v, "giant_steps")) != receipt.giant_steps
            || result.as_ref().and_then(|v| field_u64(v, "seed_ops")) != receipt.seed_ops
            || result.as_ref().and_then(|v| field_u64(v, "h2d_bytes")) != receipt.h2d_bytes
            || result.as_ref().and_then(|v| field_u64(v, "d2h_bytes")) != receipt.d2h_bytes
            || result.as_ref().and_then(|v| field_u64(v, "d2d_bytes")) != receipt.d2d_bytes
        {
            issues.push("cost counters differ from raw output".into());
        }
        for (key, value) in [
            ("device_seconds", receipt.device_seconds),
            ("host_seconds", receipt.host_seconds),
            ("setup_seconds", receipt.setup_seconds),
            (
                "candidate_verification_seconds",
                receipt.candidate_verification_seconds,
            ),
            ("transfer_seconds", receipt.transfer_seconds),
        ] {
            let replay = result
                .as_ref()
                .and_then(|v| v.get(key))
                .and_then(|s| s.parse::<f64>().ok());
            if replay != value {
                issues.push(format!("{key} differs from raw output"));
            }
        }
        let mut blockers = gpu_admission(
            &receipt.telemetry,
            command.child_pid,
            command.started_unix_ms,
            command.elapsed_ms,
        );
        if receipt.status != "verified" {
            blockers.push("solve not verified".into());
            blockers.sort();
            blockers.dedup();
        }
        if blockers != receipt.wall_blockers
            || (blockers.is_empty() && receipt.wall_admission != "admitted")
            || (!blockers.is_empty() && receipt.wall_admission != "refused")
        {
            issues.push("GPU contention admission does not regrade".into());
        }
    } else if receipt.status == "verified" {
        issues.push("verified receipt lacks solve command".into());
    }
    if let Some(expected) = &receipt.cpu_record_sha256 {
        let raw = std::fs::read(dir.join("cpu-record.json")).map_err(|e| e.to_string())?;
        let rec: Record = serde_json::from_slice(&raw).map_err(|e| e.to_string())?;
        if sha256_hex(&raw) != *expected
            || receipt.cpu_record_id.as_ref() != Some(&rec.record_id)
            || rec.workload.workload_id != workload.workload_id
            || !rec.counts()
            || !super::record::check_line_seal(
                std::str::from_utf8(&raw)
                    .map_err(|e| e.to_string())?
                    .trim_end(),
                &rec.record_id,
            )
        {
            issues.push("CPU record binding mismatch".into());
        }
    }
    Ok(issues)
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;
    use crate::cryptanalysis::ecbench::workload::{CurveSpec, TargetKind};
    use std::os::unix::fs::PermissionsExt;

    #[test]
    fn preexisting_device_memory_pressure_refuses_wall_time() {
        let sample = GpuSample {
            unix_ms: 20,
            completed_unix_ms: 21,
            gpu_csv: Some("0, GPU-test, 1500, 5000, 0, 0, 4, 16, 950, 1000".into()),
            pcie_traffic_csv: Some("1, 1".into()),
            pcie_error: None,
            processes_csv: Some(String::new()),
            error: None,
        };
        let blockers = gpu_admission(&[sample.clone(), sample], 42, 10, 100);
        assert!(blockers
            .iter()
            .any(|reason| reason.contains("memory exceeded 90 percent")));
    }

    #[test]
    fn unavailable_pcie_counters_do_not_erase_gpu_identity_or_clock_samples() {
        let sample = GpuSample {
            unix_ms: 20,
            completed_unix_ms: 21,
            gpu_csv: Some("0, GPU-test, 1500, 5000, 0, 0, 4, 16, 100, 1000".into()),
            pcie_traffic_csv: None,
            pcie_error: Some("field unsupported".into()),
            processes_csv: Some(String::new()),
            error: None,
        };
        let blockers = gpu_admission(&[sample.clone(), sample], 42, 10, 100);
        assert_eq!(
            blockers,
            vec!["PCIe traffic counters unavailable: field unsupported".to_string()]
        );
    }

    #[test]
    fn selftest_gates_stage_and_outputs_are_bound() {
        let root = std::env::temp_dir().join(format!(
            "ecbench-gpu-stage-test-{}-{}",
            std::process::id(),
            unix_ms()
        ));
        std::fs::create_dir(&root).unwrap();
        let driver = root.join("driver.sh");
        std::fs::write(&driver, "#!/bin/sh\nif [ \"$1\" = selftest ]; then echo 'selftest OK'; else echo '1.0 Gstep/s'; fi\n").unwrap();
        std::fs::set_permissions(&driver, std::fs::Permissions::from_mode(0o700)).unwrap();
        let out = root.join("ok");
        let r = run(&GpuOptions {
            binary: driver.clone(),
            out: out.clone(),
            mode: "rho".into(),
            args: vec!["--iters".into(), "8".into()],
            timeout_seconds: 3,
        })
        .unwrap();
        assert_eq!(r.status, "stage_completed");
        assert_eq!(
            r.benchmark.as_ref().unwrap().argv[1..],
            ["rho", "--iters", "8"]
        );
        let raw = std::fs::read(out.join("benchmark.stdout")).unwrap();
        assert_eq!(
            r.benchmark.as_ref().unwrap().stdout_sha256,
            sha256_hex(&raw)
        );
        assert!(run(&GpuOptions {
            binary: driver.clone(),
            out,
            mode: "rho".into(),
            args: vec![],
            timeout_seconds: 3
        })
        .is_err());

        std::fs::write(&driver, "#!/bin/sh\nif [ \"$1\" = selftest ]; then echo 'SELFTEST FAILED'; exit 1; fi\necho should-not-run\n").unwrap();
        let bad = root.join("bad");
        let r = run(&GpuOptions {
            binary: driver,
            out: bad.clone(),
            mode: "rho".into(),
            args: vec![],
            timeout_seconds: 3,
        })
        .unwrap();
        assert_eq!(r.status, "selftest_failed");
        assert!(r.benchmark.is_none());
        assert!(!bad.join("benchmark.stdout").exists());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn matched_receipt_replays_and_detects_raw_tampering() {
        let root = std::env::temp_dir().join(format!(
            "ecbench-gpu-match-test-{}-{}",
            std::process::id(),
            unix_ms()
        ));
        std::fs::create_dir(&root).unwrap();
        let (workload, inst) = Workload::build(
            &CurveSpec::PrimeSearch { bits: 12, seed: 17 },
            9,
            0,
            TargetKind::Planted,
        )
        .unwrap();
        let Instance::Prime(prime) = inst else {
            unreachable!()
        };
        let workload_path = root.join("workload.json");
        std::fs::write(&workload_path, serde_json::to_vec(&workload).unwrap()).unwrap();
        let driver = root.join("driver.sh");
        std::fs::write(&driver, format!(
            "#!/bin/sh\nif [ \"$1\" = selftest ]; then echo 'selftest OK'; exit 0; fi\nsleep 1\necho 'ECBENCH_CURVE p={} a={} b={} r={} gx={} gy={} cofactor={}'\necho 'ECBENCH_RESULT k=0x{:064x} baby_steps=100 giant_steps=50 seed_ops=10 device_seconds=0.1 host_seconds=0.2 setup_seconds=0.02 candidate_verification_seconds=0.01 transfer_seconds=0.01 h2d_bytes=64 d2h_bytes=64 d2d_bytes=0'\n",
            prime.curve.p, prime.curve.a, prime.curve.b, prime.r,
            prime.generator.0, prime.generator.1, prime.cofactor,
            workload.planted.unwrap()
        )).unwrap();
        std::fs::set_permissions(&driver, std::fs::Permissions::from_mode(0o700)).unwrap();
        let telemetry = root.join("telemetry.sh");
        std::fs::write(&telemetry, "#!/bin/sh\ncase \"$1\" in\n  --query-gpu=pcie.rx_util*) echo '1, 1';;\n  --query-gpu=*) echo '0, GPU-test, 1500, 5000, 20, 10, 4, 16, 100, 1000';;\n  *) echo 'No running processes found';;\nesac\n").unwrap();
        std::fs::set_permissions(&telemetry, std::fs::Permissions::from_mode(0o700)).unwrap();
        let out = root.join("receipt");
        let receipt = run_match(&GpuMatchOptions {
            binary: driver,
            workload: workload_path,
            out: out.clone(),
            timeout_seconds: 5,
            telemetry_binary: telemetry,
            cpu_record: None,
            args: vec!["--iters".into(), "8".into()],
        })
        .unwrap();
        assert_eq!(receipt.status, "verified");
        assert_eq!(receipt.wall_admission, "admitted");
        assert!(audit_match(&out).unwrap().is_empty());
        std::fs::write(out.join("solve.stdout"), "tampered\n").unwrap();
        assert!(!audit_match(&out).unwrap().is_empty());
        let bad_driver = root.join("driver.sh");
        std::fs::write(&bad_driver, "#!/bin/sh\necho 'SELFTEST FAILED'\nexit 1\n").unwrap();
        let bad = root.join("failed-receipt");
        let failed = run_match(&GpuMatchOptions {
            binary: bad_driver,
            workload: root.join("workload.json"),
            out: bad.clone(),
            timeout_seconds: 5,
            telemetry_binary: root.join("telemetry.sh"),
            cpu_record: None,
            args: Vec::new(),
        })
        .unwrap();
        assert_eq!(failed.status, "selftest_failed");
        assert!(audit_match(&bad).unwrap().is_empty());
        std::fs::remove_dir_all(root).unwrap();
    }
}
