//! Bounded native driver; geometry and point-count stages have separate receipts.
use std::{
    env,
    fs::{self, File},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::{Duration, Instant},
};
#[allow(dead_code)]
#[path = "../../../src/hash/sha256.rs"]
mod sha256;
fn digest(path: &Path) -> String {
    let mut h = sha256::Sha256::new();
    h.update(&fs::read(path).unwrap());
    h.finalize().iter().map(|x| format!("{x:02x}")).collect()
}
fn main() {
    let args: Vec<_> = env::args().skip(1).collect();
    assert_eq!(args.len(), 2, "run_controls STUDY_DIR OUTPUT_DIR");
    let study = fs::canonicalize(&args[0]).unwrap();
    let out = PathBuf::from(&args[1]);
    fs::create_dir_all(&out).unwrap();
    let script = study.join("control.gp");
    assert!(!out.join("receipts.tsv").exists(), "preserve earlier run");
    fs::write(out.join("source_freeze.txt"),format!("GP source sha256={}\nNative driver sha256={}\nGP=/opt/homebrew/bin/gp\nstack_initial_bytes=268435456\ngeometry_cap_seconds=30\ncount_cap_seconds=120\nper_curve_seed=20261009+curve_number\n",digest(&script),digest(&study.join("run_controls.rs")))).unwrap();
    let cases = [
        (257u64, 3u64, 8usize),
        (509, 3, 8),
        (1009, 3, 8),
        (2003, 3, 4),
        (8191, 3, 2),
        (65537, 3, 1),
        (1048583, 3, 1),
        (4294967311, 3, 1),
        (4398046511119, 3, 1),
        (13, 5, 4),
        (257, 5, 2),
        (1009, 5, 1),
        (65537, 5, 1),
        (13, 7, 4),
        (257, 7, 1),
        (1009, 7, 1),
        (65537, 7, 1),
    ];
    let mut receipts=String::from("p\tbase_degree\todd_degree\tmode\tsamples\tseed\tcap_seconds\tstatus\texit_code\twall_ms\tstdout_sha256\tstderr_sha256\n");
    // Finish all geometric verification before the separately bounded counts.
    for mode in ["geometry", "count"] {
        for &(p, n, geometry_samples) in &cases {
            let samples = if mode == "geometry" {
                geometry_samples
            } else {
                1
            };
            let cap = if mode == "geometry" { 30 } else { 120 };
            let stdout = out.join(format!("p{p}_n{n}_{mode}.stdout"));
            let stderr = out.join(format!("p{p}_n{n}_{mode}.stderr"));
            assert!(!stdout.exists() && !stderr.exists());
            let start = Instant::now();
            let mut child = Command::new("/opt/homebrew/bin/gp")
                .args(["-q", "-f", "-s", "256M"])
                .arg(&script)
                .env("ISO1_P", p.to_string())
                .env("ISO1_A", "2")
                .env("ISO1_N", n.to_string())
                .env("ISO1_SAMPLES", samples.to_string())
                .env("ISO1_SEED", "20261009")
                .env("ISO1_MODE", mode)
                .stdout(Stdio::from(File::create(&stdout).unwrap()))
                .stderr(Stdio::from(File::create(&stderr).unwrap()))
                .spawn()
                .unwrap();
            let (code, timed_out) = loop {
                if let Some(s) = child.try_wait().unwrap() {
                    break (s.code().unwrap_or(-1), false);
                }
                if start.elapsed() >= Duration::from_secs(cap) {
                    child.kill().unwrap();
                    let s = child.wait().unwrap();
                    break (s.code().unwrap_or(-1), true);
                }
                std::thread::sleep(Duration::from_millis(100));
            };
            let text = fs::read_to_string(&stdout).unwrap();
            let marker = format!("COMPLETE|{p}|2|{n}|{samples}");
            let prefix = if mode == "geometry" {
                "GEOMETRY|"
            } else {
                "COUNT|"
            };
            let records = text.lines().filter(|l| l.starts_with(prefix)).count();
            let status = if timed_out {
                "TIMEOUT"
            } else if code == 0
                && fs::metadata(&stderr).unwrap().len() == 0
                && text.lines().any(|l| l == marker)
                && records == samples
            {
                "COMPLETE"
            } else {
                "FAILED"
            };
            let line = format!(
                "{p}\t2\t{n}\t{mode}\t{samples}\t20261009\t{cap}\t{status}\t{code}\t{}\t{}\t{}\n",
                start.elapsed().as_millis(),
                digest(&stdout),
                digest(&stderr)
            );
            print!("{line}");
            receipts.push_str(&line);
            fs::write(out.join("receipts.tsv"), &receipts).unwrap();
        }
    }
}
