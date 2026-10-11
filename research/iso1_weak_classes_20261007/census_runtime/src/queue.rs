//! Sequential census queue with frozen hashes, independent controls, and storage guards.
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[allow(dead_code)]
#[path = "../../../../src/hash/sha256.rs"]
mod sha256_source;

fn hash(path: &Path) -> String {
    let mut file = File::open(path).unwrap();
    let mut state = sha256_source::Sha256::new();
    let mut buffer = [0u8; 65536];
    loop {
        let n = file.read(&mut buffer).unwrap();
        if n == 0 {
            break;
        }
        state.update(&buffer[..n]);
    }
    state
        .finalize()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
fn free_bytes(path: &Path) -> u64 {
    let result = Command::new("df").arg("-Pk").arg(path).output().unwrap();
    assert!(result.status.success());
    let text = String::from_utf8(result.stdout).unwrap();
    text.lines()
        .last()
        .unwrap()
        .split_whitespace()
        .nth(3)
        .unwrap()
        .parse::<u64>()
        .unwrap()
        * 1024
}
fn required_bytes(p: u64) -> u64 {
    1024 * 1024 * 1024 + 2 * (p.pow(3) + 1) * 140
}
fn event(log: &mut File, p: u64, status: &str, details: &str) {
    let epoch = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs();
    writeln!(log, "{epoch}\t{p}\t{status}\t{details}").unwrap();
    log.sync_all().unwrap();
    println!("p={p} status={status} {details}");
    std::io::stdout().flush().unwrap();
}
fn run(cmd: &mut Command, stdout: &Path, stderr: &Path) -> i32 {
    let out = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(stdout)
        .unwrap();
    let err = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(stderr)
        .unwrap();
    cmd.stdout(Stdio::from(out))
        .stderr(Stdio::from(err))
        .status()
        .unwrap()
        .code()
        .unwrap_or(-1)
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert!(
        args.len() >= 2,
        "usage: iso1_census_queue OUTPUT_DIR [--after-pid PID CENSUS.csv CONTROL.csv] PRIME..."
    );
    let first_prime = if args.get(1).is_some_and(|x| x == "--after-pid") {
        assert!(args.len() >= 6);
        5
    } else {
        1
    };
    let out = PathBuf::from(&args[0]);
    fs::create_dir_all(&out).unwrap();
    let out = fs::canonicalize(out).unwrap();
    let manifest_dir = std::env::var_os("ISO1_QUEUE_RUNTIME_MANIFEST_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")));
    let manifest = manifest_dir.as_path();
    let study = manifest.parent().unwrap();
    let root = manifest.ancestors().nth(3).unwrap();
    let bin = fs::canonicalize(std::env::current_exe().unwrap()).unwrap();
    let census = bin.parent().unwrap().join("iso1_class_census");
    let audit = bin.parent().unwrap().join("iso1_census_audit");
    let gp = study.join("gp_large_prime_control.gp");
    let binary_hash = hash(&census);
    let sources = [
        root.join("examples/iso1_class_census.rs"),
        manifest.join("src/field_kernel.rs"),
        manifest.join("src/point_count_kernel.rs"),
        manifest.join("src/scalar_kernel.rs"),
        manifest.join("Cargo.lock"),
        gp.clone(),
    ];
    let mut log = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(out.join("queue_status.tsv"))
        .unwrap();
    let mut freeze = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(out.join("source_freeze.txt"))
        .unwrap();
    writeln!(freeze,"census_binary={} sha256={}\naudit_binary={} sha256={}\nqueue_binary={} sha256={}\nRAYON_NUM_THREADS=4\nHost: correctness campaign; wall time is resource accounting.",census.display(),binary_hash,audit.display(),hash(&audit),bin.display(),hash(&bin)).unwrap();
    for path in sources {
        writeln!(freeze, "{} {}", hash(&path), path.display()).unwrap();
    }
    for program in [
        vec!["uname", "-a"],
        vec!["rustc", "--version"],
        vec!["gp", "--version"],
    ] {
        let output = Command::new(program[0])
            .args(&program[1..])
            .output()
            .unwrap();
        writeln!(
            freeze,
            "{} {}{}",
            program.join(" "),
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        )
        .unwrap();
    }
    freeze.sync_all().unwrap();
    if first_prime == 5 {
        let pid: u32 = args[2].parse().unwrap();
        event(
            &mut log,
            53,
            "WAITING_FOR_PRIOR_CENSUS",
            &format!("pid={pid} census={}", args[3]),
        );
        loop {
            let proc = Command::new("ps")
                .args(["-p", &pid.to_string(), "-o", "command="])
                .output()
                .unwrap();
            if !String::from_utf8_lossy(&proc.stdout).contains("iso1_class_census 53") {
                break;
            }
            std::thread::sleep(std::time::Duration::from_secs(5));
        }
        let prior_control = fs::read_to_string(&args[4]).unwrap();
        assert_eq!(prior_control.lines().count(), 5000);
        let ac = run(
            Command::new(&audit)
                .arg(&args[3])
                .arg("samples")
                .arg(&args[4]),
            &out.join("prior_p53_validation.txt"),
            &out.join("prior_p53_validation.stderr"),
        );
        if ac != 0 {
            event(
                &mut log,
                53,
                "FAILED_PRIOR_CENSUS_AUDIT",
                &format!("exit_code={ac}"),
            );
            return;
        }
        event(
            &mut log,
            53,
            "PRIOR_CENSUS_VALIDATED",
            "5000 independent controls agree",
        );
    }
    let mut seen = std::collections::BTreeSet::new();
    for arg in &args[first_prime..] {
        let p = arg.parse::<u64>().unwrap();
        assert!(p > 3 && p <= 199 && (2..).take_while(|d| d * d <= p).all(|d| p % d != 0));
        assert!(seen.insert(p), "duplicate prime");
        let free = free_bytes(&out);
        let required = required_bytes(p);
        if free < required {
            event(
                &mut log,
                p,
                "STOPPED_STORAGE_PREFLIGHT",
                &format!("available_bytes={free} required_bytes={required}"),
            );
            return;
        }
        assert_eq!(
            hash(&census),
            binary_hash,
            "census binary changed after freeze"
        );
        let csv = out.join(format!("p{p}_absolute.csv"));
        assert!(
            !csv.exists() && !out.join(format!("p{p}_absolute.csv.zst")).exists(),
            "preserve existing evidence"
        );
        event(
            &mut log,
            p,
            "RUNNING",
            &format!(
                "expected_point_counts={} available_bytes={free}",
                (p.pow(4) + 3 * p * p + 8) / 12
            ),
        );
        let start = Instant::now();
        let code = run(
            Command::new(&census)
                .arg(p.to_string())
                .arg(&csv)
                .args(["--derive-twists", "--absolute-orbit-quotient"])
                .env("RAYON_NUM_THREADS", "4"),
            &out.join(format!("p{p}_census.stdout")),
            &out.join(format!("p{p}_census.stderr")),
        );
        let mut receipt = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(out.join(format!("p{p}_receipt.txt")))
            .unwrap();
        writeln!(receipt,"prime={p}\ncommand=RAYON_NUM_THREADS=4 {} {p} {} --derive-twists --absolute-orbit-quotient\nexit_code={code}\nresource_wall_seconds={:.3}\ncensus_binary_sha256={binary_hash}",census.display(),csv.display(),start.elapsed().as_secs_f64()).unwrap();
        if code != 0 {
            event(&mut log, p, "FAILED_CENSUS", &format!("exit_code={code}"));
            return;
        }
        assert_eq!(hash(&census), binary_hash);
        writeln!(receipt, "csv_sha256={}", hash(&csv)).unwrap();
        receipt.sync_all().unwrap();
        let control = out.join(format!("gp_p{p}_control.csv"));
        let gp_error = out.join(format!("gp_p{p}_control.stderr"));
        let gp_code = run(
            Command::new("gp")
                .args(["-q", "-f"])
                .arg(&gp)
                .env("ISO1_P", p.to_string()),
            &control,
            &gp_error,
        );
        let lines = fs::read_to_string(&control).unwrap().lines().count();
        writeln!(
            receipt,
            "gp_exit_code={gp_code}\ngp_samples={lines}\ngp_sha256={}",
            hash(&control)
        )
        .unwrap();
        if gp_code != 0 || fs::metadata(gp_error).unwrap().len() != 0 || lines != 5000 {
            event(
                &mut log,
                p,
                "FAILED_GP_CONTROL",
                &format!("exit_code={gp_code} samples={lines}"),
            );
            return;
        }
        let ac = run(
            Command::new(&audit).arg(&csv).arg("samples").arg(&control),
            &out.join(format!("p{p}_validation.txt")),
            &out.join(format!("p{p}_validation.stderr")),
        );
        writeln!(receipt, "audit_exit_code={ac}").unwrap();
        receipt.sync_all().unwrap();
        if ac != 0 {
            event(&mut log, p, "FAILED_AUDIT", &format!("exit_code={ac}"));
            return;
        }
        let compressed = out.join(format!("p{p}_absolute.csv.zst"));
        let compression = Command::new("zstd")
            .args(["-q", "-6"])
            .arg(&csv)
            .arg("-o")
            .arg(&compressed)
            .status()
            .unwrap();
        if !compression.success() {
            event(&mut log, p, "FAILED_COMPRESSION", "raw CSV preserved");
            return;
        }
        let archive_test = Command::new("zstd")
            .args(["-q", "-t"])
            .arg(&compressed)
            .status()
            .unwrap();
        if !archive_test.success() {
            event(&mut log, p, "FAILED_ARCHIVE_TEST", "raw CSV preserved");
            return;
        }
        let mut decoder = Command::new("zstd")
            .args(["-q", "-d", "-c"])
            .arg(&compressed)
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        let mut state = sha256_source::Sha256::new();
        let mut buffer = [0u8; 65536];
        let mut stream = decoder.stdout.take().unwrap();
        loop {
            let n = stream.read(&mut buffer).unwrap();
            if n == 0 {
                break;
            }
            state.update(&buffer[..n]);
        }
        let replay_hash: String = state
            .finalize()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        assert!(decoder.wait().unwrap().success());
        assert_eq!(replay_hash, hash(&csv), "compressed replay differs");
        writeln!(receipt,"archive_sha256={}\narchive_replay_sha256={replay_hash}\nstatus=COMPLETE_DIAGNOSTIC_WITH_GP_CONTROL",hash(&compressed)).unwrap();
        receipt.sync_all().unwrap();
        fs::remove_file(&csv).unwrap();
        event(
            &mut log,
            p,
            "COMPLETE_DIAGNOSTIC_WITH_GP_CONTROL",
            &format!(
                "compressed_csv={} replay_sha256={replay_hash}",
                compressed.display()
            ),
        );
    }
    event(
        &mut log,
        0,
        "QUEUE_COMPLETE",
        "all requested primes in this queue completed",
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn storage_guard_covers_a_large_raw_csv_and_reserve() {
        assert!(required_bytes(199) > 3_200_000_000);
        assert!(required_bytes(53) > 1024 * 1024 * 1024);
    }
}
