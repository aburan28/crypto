use serde_json::Value;
use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicU64, Ordering},
};

static NEXT: AtomicU64 = AtomicU64::new(0);
fn directory() -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "ic-fixed-test-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir(&path).unwrap();
    path
}
fn run(path: &PathBuf, stage: &str, pairs: u32) -> (bool, Value) {
    let output = Command::new(env!("CARGO_BIN_EXE_ic"))
        .args(["fixed", "--params"])
        .arg(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("docs/ic/params/k0n9-fixed.json"))
        .arg("--dir")
        .arg(path)
        .args([
            "--stage",
            stage,
            "--pair-budget",
            &pairs.to_string(),
            "--attempts",
            "256",
            "--query-seconds",
            "2",
            "--json",
        ])
        .output()
        .unwrap();
    let report = serde_json::from_slice(&output.stdout)
        .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&output.stderr)));
    (output.status.success(), report)
}

#[test]
fn fixed_cli_preserves_incomplete_then_resumes_and_reuses() {
    let path = directory();
    let (success, partial) = run(&path, "pairs", 7);
    assert!(!success);
    assert_eq!(partial["status"], "pair_budget");
    let (success, complete) = run(&path, "all", 10000);
    assert!(success, "{complete}");
    assert_eq!(complete["status"], "complete");
    assert!(complete["targets"]
        .as_object()
        .unwrap()
        .values()
        .all(|v| v["verified"] == true));
    let (success, replay) = run(&path, "all", 0);
    assert!(success, "{replay}");
    assert!(replay["reuse"]["new_attempts"].is_null());
    assert!(replay["reuse"]["pair_candidates_built"].is_null());
    assert_eq!(replay["relations_saved"], complete["relations_saved"]);
    fs::remove_dir_all(path).unwrap();
}

#[test]
fn fixed_cli_reports_missing_interpreter() {
    let path = directory();
    let output = Command::new(env!("CARGO_BIN_EXE_ic"))
        .args([
            "fixed",
            "--params",
            "unused.json",
            "--python",
            "/nonexistent/ic-python",
            "--dir",
        ])
        .arg(&path)
        .arg("--json")
        .output()
        .unwrap();
    assert!(!output.status.success());
    let report: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["status"], "error");
    assert!(report["message"]
        .as_str()
        .unwrap()
        .contains("could not run fixed workflow"));
    fs::remove_dir_all(path).unwrap();
}
