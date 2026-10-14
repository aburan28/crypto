//! A bounded, provenance-preserving front door to `gpu/ecc/bench`.
//!
//! These are GPU stage benchmarks on that binary's own fixed curves and
//! targets. They do not enter the matched-workload `ecbench run` records.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use serde::{Deserialize, Serialize};

use super::canonical::sha256_hex;

pub const GPU_STAGE_SCHEMA: &str = "ecbench.gpu_stage/v1";

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
    let child = cmd
        .spawn()
        .map_err(|e| format!("{}: {e}", binary.display()))?;
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

#[cfg(all(test, unix))]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

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
}
