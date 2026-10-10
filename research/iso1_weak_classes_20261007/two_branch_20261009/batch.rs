//! Four bounded workers, one process and immutable receipt per independent source.
use crypto_lib::hash::sha256::sha256;
use serde_json::json;
use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
    thread,
    time::Instant,
};
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 3, "STUDY INPUTS OUTPUT_DIR");
    let study = fs::canonicalize(&a[0]).unwrap();
    let inputs = fs::canonicalize(&a[1]).unwrap();
    let out = PathBuf::from(&a[2]);
    assert!(!out.exists(), "preserve previous batch");
    fs::create_dir_all(&out).unwrap();
    let out = fs::canonicalize(out).unwrap();
    let source = fs::read_to_string(&inputs).unwrap();
    let jobs: Vec<_> = source
        .lines()
        .filter(|s| s.starts_with("construct_source("))
        .map(str::to_string)
        .collect();
    assert!(!jobs.is_empty());
    let exe = std::env::current_exe()
        .unwrap()
        .with_file_name("iso1_two_branch_run");
    assert!(exe.exists());
    fs::write(out.join("batch_freeze.json"),serde_json::to_vec_pretty(&json!({"schema":"iso1.two_branch_batch/v1","input_sha256":hex::encode(sha256(source.as_bytes())),"sources":jobs.len(),"workers":4,"per_source_parent_cap_seconds":240,"search_cap_ms":90000,"vertex_cap":256,"batch_source_sha256":hex::encode(sha256(&fs::read(study.join("batch.rs")).unwrap())),"batch_binary_sha256":hex::encode(sha256(&fs::read(std::env::current_exe().unwrap()).unwrap()))})).unwrap()).unwrap();
    let count = jobs.len();
    let jobs = Arc::new(jobs);
    let cursor = Arc::new(AtomicUsize::new(0));
    let done = Arc::new(AtomicUsize::new(0));
    let failed = Arc::new(AtomicUsize::new(0));
    let start = Instant::now();
    let mut handles = Vec::new();
    for _ in 0..4 {
        let jobs = jobs.clone();
        let cursor = cursor.clone();
        let done = done.clone();
        let failed = failed.clone();
        let study = study.clone();
        let out = out.clone();
        let exe = exe.clone();
        handles.push(thread::spawn(move || loop {
            let i = cursor.fetch_add(1, Ordering::Relaxed);
            if i >= jobs.len() {
                break;
            }
            let cell = out.join(format!("source{:04}", i + 1));
            fs::create_dir(&cell).unwrap();
            let input = cell.join("input.gp");
            fs::write(
                &input,
                format!("{}\nprint(\"CONSTRUCTION_COMPLETE|sources=1\");\n", jobs[i]),
            )
            .unwrap();
            let result = Command::new(&exe)
                .arg(&study)
                .arg(&cell)
                .args(["construct.gp", "240"])
                .arg(format!("ISO1_INPUTS={}", input.display()))
                .output()
                .unwrap();
            fs::write(cell.join("driver_stdout.txt"), &result.stdout).unwrap();
            fs::write(cell.join("driver_stderr.txt"), &result.stderr).unwrap();
            if !result.status.success() {
                failed.fetch_add(1, Ordering::Relaxed);
            }
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            if d % 32 == 0 || d == jobs.len() {
                println!(
                    "completed={d}/{} failed={}",
                    jobs.len(),
                    failed.load(Ordering::Relaxed)
                );
            }
        }));
    }
    for h in handles {
        h.join().unwrap();
    }
    let result = json!({"sources":count,"completed":done.load(Ordering::Relaxed),"failed":failed.load(Ordering::Relaxed),"wall_ms":start.elapsed().as_millis().to_string(),"scope":"shared-host resource diagnostic"});
    fs::write(
        out.join("batch_receipt.json"),
        serde_json::to_vec_pretty(&result).unwrap(),
    )
    .unwrap();
    println!("{}", result);
}
