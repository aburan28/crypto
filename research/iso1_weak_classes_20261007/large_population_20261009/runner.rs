//! Native bounded process driver. Preserve every attempted cell and raw output.
use crypto_lib::hash::sha256::sha256;
use std::{
    env,
    fs::{self, File},
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
    thread,
    time::{Duration, Instant},
};
fn hash(p: &Path) -> String {
    hex::encode(sha256(&fs::read(p).unwrap()))
}
const FIELDS: [(u64, u64, usize); 10] = [
    (4294967291, 3, 64),
    (17179869143, 3, 64),
    (68719476731, 3, 64),
    (274877906899, 3, 64),
    (1099511627689, 3, 64),
    (4398046511093, 3, 64),
    (602233, 5, 32),
    (38543917, 5, 32),
    (13421, 7, 32),
    (262139, 7, 32),
];
#[derive(Clone)]
struct Job {
    field: usize,
    p: u64,
    n: u64,
    sample: usize,
    seed: u64,
    mode: &'static str,
}
fn main() {
    let args: Vec<_> = env::args().skip(1).collect();
    assert_eq!(
        args.len(),
        3,
        "iso1_population_run STUDY_DIR OUTPUT_DIR small|pilot|population|controls"
    );
    let study = fs::canonicalize(&args[0]).unwrap();
    let out = PathBuf::from(&args[1]);
    let phase = &args[2];
    let mut jobs = Vec::new();
    if phase == "small" {
        for i in 1..=64 {
            jobs.push(Job {
                field: 0,
                p: 7,
                n: 3,
                sample: i,
                seed: 202610090000 + i as u64,
                mode: "population",
            });
        }
        for i in 1..=4 {
            jobs.push(Job {
                field: 0,
                p: 7,
                n: 3,
                sample: i,
                seed: 202610099000 + i as u64,
                mode: "control",
            });
        }
    } else {
        assert!(["pilot", "population", "controls"].contains(&phase.as_str()));
        for (f, &(p, n, samples)) in FIELDS.iter().enumerate() {
            let mode = if phase == "controls" {
                "control"
            } else {
                "population"
            };
            let count = if phase == "population" { samples } else { 1 };
            for sample in 1..=count {
                let seed = if mode == "control" {
                    202610099000 + 10000 * (f as u64 + 1)
                } else {
                    202610090000 + 10000 * (f as u64 + 1) + sample as u64
                };
                jobs.push(Job {
                    field: f + 1,
                    p,
                    n,
                    sample,
                    seed,
                    mode,
                });
            }
        }
    }
    fs::create_dir_all(&out).unwrap();
    assert!(!out.join("receipts.tsv").exists(), "preserve previous run");
    let script = study.join("sample.gp");
    let exe = fs::canonicalize(env::current_exe().unwrap()).unwrap();
    let gp = Path::new("/opt/homebrew/bin/gp");
    let revision = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(&study)
        .output()
        .unwrap();
    let version = Command::new(gp).arg("--version-short").output().unwrap();
    fs::write(out.join("source_freeze.txt"),format!(
        "revision={}protocol_sha256={}\ngp_script_sha256={}\nnative_driver_sha256={}\nrunner_binary_sha256={}\ngp_binary_sha256={}\ngp_version={}phase={}\nworkers=4\nparent_cap_seconds=240\nsearch_cap_ms=90000\nvertex_cap=256\nstack_initial_bytes=268435456\npopulation_seed=202610090000+10000*field_number+sample_number\n",
        String::from_utf8(revision.stdout).unwrap(),hash(&study.join("PROTOCOL.md")),hash(&script),hash(&study.join("runner.rs")),hash(&exe),hash(gp),String::from_utf8(version.stdout).unwrap(),phase)).unwrap();
    let mut receipt = File::create(out.join("receipts.tsv")).unwrap();
    writeln!(receipt,"field\tp\tn\tsample\tseed\tmode\tstatus\texit_code\twall_ms\tstdout_sha256\tstderr_sha256\tstdout_file").unwrap();
    let receipts = Arc::new(Mutex::new(receipt));
    let jobs = Arc::new(jobs);
    let cursor = Arc::new(AtomicUsize::new(0));
    let completed = Arc::new(AtomicUsize::new(0));
    let start = Instant::now();
    let mut handles = Vec::new();
    for _ in 0..4 {
        let jobs = jobs.clone();
        let cursor = cursor.clone();
        let receipts = receipts.clone();
        let completed = completed.clone();
        let out = out.clone();
        let script = script.clone();
        handles.push(thread::spawn(move || loop {
            let j = cursor.fetch_add(1, Ordering::SeqCst);
            if j >= jobs.len() {
                break;
            }
            let job = &jobs[j];
            let name = format!(
                "f{:02}_p{}_n{}_{}_s{:03}",
                job.field, job.p, job.n, job.mode, job.sample
            );
            let stdout = out.join(format!("{name}.stdout"));
            let stderr = out.join(format!("{name}.stderr"));
            assert!(!stdout.exists() && !stderr.exists());
            let now = Instant::now();
            let mut child = Command::new("/opt/homebrew/bin/gp")
                .args(["-q", "-f", "-s", "256M"])
                .arg(&script)
                .env("ISO1_P", job.p.to_string())
                .env("ISO1_N", job.n.to_string())
                .env("ISO1_SEED", job.seed.to_string())
                .env("ISO1_MODE", job.mode)
                .env("ISO1_BUDGET", "256")
                .stdout(Stdio::from(File::create(&stdout).unwrap()))
                .stderr(Stdio::from(File::create(&stderr).unwrap()))
                .spawn()
                .unwrap();
            let (code, timeout) = loop {
                if let Some(code) = child.try_wait().unwrap() {
                    break (code.code().unwrap_or(-1), false);
                }
                if now.elapsed() >= Duration::from_secs(240) {
                    child.kill().unwrap();
                    break (child.wait().unwrap().code().unwrap_or(-1), true);
                }
                thread::sleep(Duration::from_millis(100));
            };
            let raw = fs::read_to_string(&stdout).unwrap();
            let marker = format!("COMPLETE|{}|{}|{}|{}", job.p, job.n, job.seed, job.mode);
            let ok = code == 0
                && fs::metadata(&stderr).unwrap().len() == 0
                && raw.lines().any(|l| l == marker)
                && raw.lines().filter(|l| l.starts_with("SAMPLE|")).count() == 1
                && raw.lines().filter(|l| l.starts_with("SEARCH|")).count() == 1;
            let status = if timeout {
                "TIMEOUT"
            } else if ok {
                "COMPLETE"
            } else {
                "FAILED"
            };
            let line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{name}.stdout\n",
                job.field,
                job.p,
                job.n,
                job.sample,
                job.seed,
                job.mode,
                status,
                code,
                now.elapsed().as_millis(),
                hash(&stdout),
                hash(&stderr)
            );
            let mut lock = receipts.lock().unwrap();
            lock.write_all(line.as_bytes()).unwrap();
            lock.flush().unwrap();
            let done = completed.fetch_add(1, Ordering::SeqCst) + 1;
            if status != "COMPLETE" || done % 16 == 0 || done == jobs.len() {
                print!("completed={done}/{} {line}", jobs.len());
            }
        }));
    }
    for h in handles {
        h.join().unwrap();
    }
    fs::write(
        out.join("process_complete.txt"),
        format!(
            "attempts={}\nwall_ms={}\n",
            jobs.len(),
            start.elapsed().as_millis()
        ),
    )
    .unwrap();
}
