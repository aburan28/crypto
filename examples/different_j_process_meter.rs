//! Measures one complete native batch solver child, including startup and stdout.

use std::fs::{self, File};
use std::os::unix::process::ExitStatusExt;
use std::process::{Command, Stdio};
use std::time::Instant;

use serde_json::json;

fn timeval_us(time: libc::timeval) -> u64 {
    (time.tv_sec as u64) * 1_000_000 + time.tv_usec as u64
}

fn child_usage() -> libc::rusage {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::zeroed();
    assert_eq!(
        unsafe { libc::getrusage(libc::RUSAGE_CHILDREN, usage.as_mut_ptr()) },
        0
    );
    unsafe { usage.assume_init() }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        8,
        "usage: different_j_process_meter BINARY PUBLIC ROUND ARM STDOUT STDERR METER.json"
    );
    let stdout = File::create(&args[5]).unwrap();
    let stderr = File::create(&args[6]).unwrap();
    let wall = Instant::now();
    let status = Command::new(&args[1])
        .arg("batch-run")
        .arg(&args[2])
        .arg(&args[3])
        .arg(&args[4])
        .stdout(Stdio::from(stdout))
        .stderr(Stdio::from(stderr))
        .status()
        .unwrap();
    let wall_ns = wall.elapsed().as_nanos() as u64;
    let usage = child_usage();
    let user_cpu_us = timeval_us(usage.ru_utime);
    let system_cpu_us = timeval_us(usage.ru_stime);
    let peak_rss_bytes = if cfg!(target_os = "macos") {
        usage.ru_maxrss as u64
    } else {
        usage.ru_maxrss as u64 * 1024
    };
    let receipt = json!({"schema":"different_j_child_process_meter/v1",
        "round":args[3].parse::<u64>().unwrap(),"arm":args[4],
        "exit_code":status.code(),"signal":status.signal(),
        "wall_ns":wall_ns,"child_user_cpu_us":user_cpu_us,
        "child_system_cpu_us":system_cpu_us,
        "child_cpu_us":user_cpu_us+system_cpu_us,
        "child_peak_rss_bytes":peak_rss_bytes,
        "stdout_bytes":fs::metadata(&args[5]).unwrap().len(),
        "stderr_bytes":fs::metadata(&args[6]).unwrap().len()});
    fs::write(&args[7], serde_json::to_vec(&receipt).unwrap()).unwrap();
    println!("{receipt}");
    if !status.success() {
        std::process::exit(status.code().unwrap_or(128 + status.signal().unwrap_or(0)));
    }
}
