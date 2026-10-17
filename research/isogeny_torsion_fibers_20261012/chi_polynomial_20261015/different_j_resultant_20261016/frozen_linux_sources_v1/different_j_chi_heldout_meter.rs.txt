//! Whole-child CPU, wall and RSS meter for a frozen held-out chi run.

use std::fs::{self, File};
use std::os::unix::process::{CommandExt, ExitStatusExt};
use std::process::{Command, Stdio};
use std::time::Instant;

use serde_json::json;

const LIMIT_BYTES: u64 = 512 * 1024 * 1024;

fn timeval_us(t: libc::timeval) -> u64 {
    t.tv_sec as u64 * 1_000_000 + t.tv_usec as u64
}

fn children() -> libc::rusage {
    let mut raw = std::mem::MaybeUninit::<libc::rusage>::zeroed();
    assert_eq!(unsafe { libc::getrusage(libc::RUSAGE_CHILDREN, raw.as_mut_ptr()) }, 0);
    unsafe { raw.assume_init() }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 10,
        "usage: meter BINARY PUBLIC ROUND ARM POSITION STDOUT STDERR RECEIPT LIMIT_SECONDS");
    let stdout = File::create(&args[6]).unwrap();
    let stderr = File::create(&args[7]).unwrap();
    let wall_limit: u64 = args[9].parse().unwrap();
    assert_eq!(wall_limit, 300);
    let start = Instant::now();
    let mut command = Command::new(&args[1]);
    command
        .arg("run-heldout").arg(&args[2]).arg(&args[3]).arg(&args[4]).arg(&args[5])
        .stdout(Stdio::from(stdout)).stderr(Stdio::from(stderr));
    unsafe {
        command.pre_exec(|| {
            let limit = libc::rlimit {
                rlim_cur: LIMIT_BYTES as libc::rlim_t,
                rlim_max: LIMIT_BYTES as libc::rlim_t,
            };
            if !cfg!(target_os = "macos") && libc::setrlimit(libc::RLIMIT_AS, &limit) != 0 {
                return Err(std::io::Error::last_os_error());
            }
            Ok(())
        });
    }
    let mut child = command.spawn().unwrap();
    let (status, timed_out) = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break (status, false);
        }
        if start.elapsed().as_secs() >= wall_limit {
            child.kill().unwrap();
            break (child.wait().unwrap(), true);
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    };
    let elapsed = start.elapsed();
    let usage = children();
    let rss = if cfg!(target_os = "macos") {
        usage.ru_maxrss as u64
    } else {
        usage.ru_maxrss as u64 * 1024
    };
    let record = json!({
        "schema":"different_j_chi_heldout_meter/v1",
        "round":args[3].parse::<u64>().unwrap(),"arm":args[4],
        "position":args[5].parse::<usize>().unwrap(),
        "exit_code":status.code(),"signal":status.signal(),
        "wall_ns":elapsed.as_nanos(),"wall_limit_ns":wall_limit*1_000_000_000,
        "memory_limit_bytes":LIMIT_BYTES,
        "rlimit_as_applied":!cfg!(target_os = "macos"),
        "rlimit_data_applied":false,
        "timed_out":timed_out,
        "child_user_cpu_us":timeval_us(usage.ru_utime),
        "child_system_cpu_us":timeval_us(usage.ru_stime),
        "child_peak_rss_bytes":rss,
        "stdout_bytes":fs::metadata(&args[6]).unwrap().len(),
        "stderr_bytes":fs::metadata(&args[7]).unwrap().len()
    });
    fs::write(&args[8], format!("{record}\n")).unwrap();
    println!("{record}");
    assert!(!timed_out && rss <= LIMIT_BYTES);
    if !status.success() {
        std::process::exit(status.code().unwrap_or(128 + status.signal().unwrap_or(0)));
    }
}
