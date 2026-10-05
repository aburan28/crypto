//! Frozen native process driver for the n37 six-sum / signed-Frobenius-rho panel.
//! Run only after the three producer/replay examples have been built in release mode.

use crypto_lib::hash::sha256;
use serde_json::{json, Map, Value};
use std::fs::{self, File, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

const CORPUS: &str = "compact-disjoint-cold-v2-n37-L1024-b01-20261001";
const POINTS: &str =
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b01.points.jsonl";
const NOTE: &str = "research/notes/ecc2k130/n37_native_residual_rho_20261005";
const WALL_CAP: u64 = 180;
const RSS_CAP: u64 = 4 * 1024 * 1024 * 1024;
const SOURCES: &[(&str, &str)] = &[
    (
        "examples/n37_native_m6_residual.rs",
        "086a6fbdd12e2b104ec18f809d4048453d747ad76c5a74f695befd08d4efffa8",
    ),
    (
        "examples/n37_native_m6_residual_replay.rs",
        "68a4ce03029b82821fdabb9854a14559dde243bd5ac0f409de884df032f7394a",
    ),
    (
        "examples/koblitz_rho_batch_ks_v3.rs",
        "98b6e8a27d821ebfc7ae716410184dffe89aadf27dfe36bcac8ce834c39cf04c",
    ),
    (
        "Cargo.toml",
        "f88ac8c9527bc2403cc550b8b00aaa778630911a4c9d1de47f67dbbff055d6a2",
    ),
    (
        "Cargo.lock",
        "b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365",
    ),
    (
        POINTS,
        "78553fdff5ae66521d3a1052978962258e78d48c92285df6e1df962034b43717",
    ),
];

fn digest(path: &Path) -> io::Result<String> {
    Ok(hex::encode(sha256(&fs::read(path)?)))
}

fn version(program: &str) -> io::Result<String> {
    let result = Command::new(program).arg("--version").output()?;
    if !result.status.success() {
        return Err(io::Error::other(format!("{program} --version failed")));
    }
    Ok(String::from_utf8_lossy(&result.stdout).trim().to_owned())
}

fn loadavg() -> [f64; 3] {
    let mut values = [0.0; 3];
    let count = unsafe { libc::getloadavg(values.as_mut_ptr(), 3) };
    if count == 3 {
        values
    } else {
        [f64::NAN; 3]
    }
}

fn timeval_ns(value: libc::timeval) -> u64 {
    value.tv_sec as u64 * 1_000_000_000 + value.tv_usec as u64 * 1_000
}

fn run_arm(
    root: &Path,
    out: &Path,
    block: usize,
    position: usize,
    method: &str,
) -> io::Result<Value> {
    let stem = format!("b{block}_{position}_{method}");
    let stdout_name = format!("{stem}.stdout");
    let stderr_name = format!("{stem}.stderr");
    let result_name = if method == "ic" {
        format!("{stem}.result.json")
    } else {
        stdout_name.clone()
    };
    let executable = root
        .join("target/release/examples")
        .join(if method == "ic" {
            "n37_native_m6_residual"
        } else {
            "koblitz_rho_batch_ks_v3"
        });
    let points = root.join(POINTS);
    let args = if method == "ic" {
        vec![out.join(&result_name).display().to_string()]
    } else {
        vec![
            "37".into(),
            "0".into(),
            "signed_frobenius".into(),
            "1024".into(),
            "2026100110101".into(),
        ]
    };
    let mut environment = Map::new();
    if method == "rho" {
        environment.insert(
            "KIC_RHO_POINT_INPUT".into(),
            json!(points.display().to_string()),
        );
        environment.insert("KIC_RHO_BATCH_CORPUS".into(), json!(CORPUS));
        environment.insert("KIC_RHO_CANON_BACKEND".into(), json!("normal_basis"));
        environment.insert("KIC_RHO_DP_BITS".into(), json!("8"));
    }
    let mut command = Command::new(&executable);
    command
        .args(&args)
        .current_dir(root)
        .stdout(Stdio::from(File::create(out.join(&stdout_name))?))
        .stderr(Stdio::from(File::create(out.join(&stderr_name))?));
    for (key, value) in &environment {
        command.env(key, value.as_str().expect("string environment"));
    }
    let before = loadavg();
    let start = Instant::now();
    let child = command.spawn()?;
    let pid = child.id() as libc::pid_t;
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    let (exit_code, timed_out) = loop {
        let mut status = 0;
        let waited = unsafe { libc::wait4(pid, &mut status, libc::WNOHANG, &mut usage) };
        if waited == pid {
            let code = if libc::WIFEXITED(status) {
                libc::WEXITSTATUS(status)
            } else if libc::WIFSIGNALED(status) {
                -libc::WTERMSIG(status)
            } else {
                -255
            };
            break (code, false);
        }
        if waited < 0 {
            let error = io::Error::last_os_error();
            if error.raw_os_error() == Some(libc::EINTR) {
                continue;
            }
            return Err(error);
        }
        if start.elapsed().as_secs_f64() > WALL_CAP as f64 {
            unsafe {
                libc::kill(pid, libc::SIGKILL);
            }
            let waited = unsafe { libc::wait4(pid, &mut status, 0, &mut usage) };
            if waited != pid {
                return Err(io::Error::last_os_error());
            }
            let code = if libc::WIFSIGNALED(status) {
                -libc::WTERMSIG(status)
            } else {
                -255
            };
            break (code, true);
        }
        thread::sleep(Duration::from_millis(10));
    };
    let wall_ns = start.elapsed().as_nanos() as u64;
    let peak = if cfg!(target_os = "macos") {
        usage.ru_maxrss as u64
    } else {
        usage.ru_maxrss as u64 * 1024
    };
    let result = out.join(&result_name);
    let stdout = out.join(&stdout_name);
    let stderr = out.join(&stderr_name);
    let mut full_command = vec![executable.display().to_string()];
    full_command.extend(args);
    Ok(json!({
        "schema": "n37-native-batch-rho-arm-v2", "block": block, "position": position,
        "method": method, "command": full_command, "environment": environment,
        "exit_code": exit_code, "timed_out": timed_out, "wall_ns": wall_ns,
        "user_ns": timeval_ns(usage.ru_utime), "system_ns": timeval_ns(usage.ru_stime),
        "peak_rss_bytes": peak, "loadavg_before": before, "loadavg_after": loadavg(),
        "stdout": stdout_name, "stdout_sha256": digest(&stdout)?,
        "stderr": stderr_name, "stderr_sha256": digest(&stderr)?,
        "result": if result.exists() { Some(result_name) } else { None },
        "result_sha256": if result.exists() { Some(digest(&result)?) } else { None },
        "resource_accepted": exit_code == 0 && !timed_out && peak <= RSS_CAP,
    }))
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args_os();
    let _program = args.next();
    let out = PathBuf::from(
        args.next()
            .ok_or("usage: n37_native_batch_rho_panel OUT_DIR")?,
    );
    if args.next().is_some() || out.exists() {
        return Err("expected one new output directory".into());
    }
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let note = root.join(NOTE);
    let mut source_sha = Map::new();
    for (path, expected) in SOURCES {
        let actual = digest(&root.join(path))?;
        if actual != *expected {
            return Err(format!("source/input hash changed: {path}").into());
        }
        source_sha.insert((*path).into(), json!(actual));
    }
    let mut binaries = Map::new();
    for (key, name) in [
        ("ic", "n37_native_m6_residual"),
        ("ic_replay", "n37_native_m6_residual_replay"),
        ("rho", "koblitz_rho_batch_ks_v3"),
    ] {
        binaries.insert(
            key.into(),
            json!(digest(&root.join("target/release/examples").join(name))?),
        );
    }
    let schedule: Vec<Vec<&str>> = (0..5)
        .map(|block| {
            if block % 2 == 0 {
                vec!["ic", "rho", "rho", "ic"]
            } else {
                vec!["rho", "ic", "ic", "rho"]
            }
        })
        .collect();
    let mut manifest = json!({
        "schema": "n37-native-batch-rho-panel-v2", "status": "running",
        "schedule": schedule, "wall_cap_seconds": WALL_CAP, "rss_cap_bytes": RSS_CAP,
        "host": {"os": std::env::consts::OS, "architecture": std::env::consts::ARCH,
                 "rustc": version("rustc")?, "cargo": version("cargo")?,
                 "isolation": "unverified_unisolated_local_host"},
        "source_sha256": source_sha, "protocol_sha256": digest(&note.join("PROTOCOL.md"))?,
        "runner_sha256": digest(&root.join("examples/n37_native_batch_rho_panel.rs"))?,
        "binary_sha256": binaries, "point_file_sha256": digest(&root.join(POINTS))?,
        "fixture_read_by_producers": false,
    });
    fs::create_dir(&out)?;
    fs::write(
        out.join("manifest.json"),
        format!("{}\n", serde_json::to_string_pretty(&manifest)?),
    )?;
    let mut stream = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(out.join("arms.jsonl"))?;
    for (block, order) in schedule.iter().enumerate() {
        for (position, method) in order.iter().enumerate() {
            let row = run_arm(root, &out, block, position, method)?;
            writeln!(stream, "{}", serde_json::to_string(&row)?)?;
            stream.sync_data()?;
            println!(
                "{}",
                json!({"block": block, "position": position,
                "method": method, "exit_code": row["exit_code"],
                "timed_out": row["timed_out"], "wall_ns": row["wall_ns"],
                "resource_accepted": row["resource_accepted"]})
            );
        }
    }
    manifest["status"] = json!("completed_schedule");
    fs::write(
        out.join("manifest.json"),
        format!("{}\n", serde_json::to_string_pretty(&manifest)?),
    )?;
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("panel failed: {error}");
        std::process::exit(1);
    }
}
