//! End to end through the `ecbench` binary: plan, run, verify with
//! replays, compare, and load the database SQL, on a session small
//! enough for CI.  Runs with `--cpus none` so it works on any host; the
//! Linux CI job (`.github/workflows/ecbench.yml`) runs a pinned session.

use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

fn exe() -> &'static str {
    env!("CARGO_BIN_EXE_ecbench")
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("ecbench-test-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn ecbench(args: &[&str]) -> (bool, String, String) {
    let out = Command::new(exe())
        .args(args)
        .output()
        .expect("ecbench runs");
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

const SPEC: &str = r#"{
  "schema": "ecbench.spec/v1",
  "label": "integration test",
  "workloads": {
    "curves": [{"kind": "prime_search", "bits": 16, "seed": 59297},
               {"kind": "koblitz", "a": 0, "n": 13}],
    "targets_per_curve": 2,
    "target_seed": 7
  },
  "arms": [
    {"name": "rho", "role": "reference", "method": {"id": "rho.negation"}},
    {"name": "bsgs", "role": "candidate", "method": {"id": "bsgs.negation"}},
    {"name": "kangaroo", "role": "candidate", "method": {"id": "kangaroo.vow"}},
    {"name": "rho-aa", "role": "control", "method": {"id": "rho.negation"}}
  ],
  "measurement": {"rounds": 2, "warmup": 1, "seed": 3, "isolation_required": "L0", "timeout_seconds": 60}
}"#;

fn lock(dir: &Path) -> String {
    // A private lock, so the test never waits on (or blocks) a real
    // benchmark holding the shared one.
    dir.join("bench.lock").display().to_string()
}

#[test]
fn a_session_runs_verifies_compares_and_loads() {
    let dir = scratch("session");
    let spec = dir.join("spec.json");
    std::fs::write(&spec, SPEC).unwrap();
    let out = dir.join("s");
    let (ok, _, err) = ecbench(&[
        "run",
        "--spec",
        spec.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--cpus",
        "none",
        "--lock",
        &lock(&dir),
        "--quiet",
    ]);
    assert!(ok, "run failed: {err}");

    // Every execution has a sealed record, all verified.
    let session: Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("session.json")).unwrap()).unwrap();
    assert_eq!(session["status"], "complete");
    assert_eq!(session["records_written"], 3 * 4 * 4); // rounds × workloads × arms
    assert_eq!(session["status_counts"]["verified"], 48);

    // A second run into the same directory is refused.
    let (ok, _, err) = ecbench(&[
        "run",
        "--spec",
        spec.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--cpus",
        "none",
        "--lock",
        &lock(&dir),
        "--quiet",
    ]);
    assert!(!ok && err.contains("never overwrites"), "{err}");

    // The audit re-derives everything and replays runs exactly.
    let (ok, stdout, err) = ecbench(&[
        "verify",
        "--dir",
        out.to_str().unwrap(),
        "--replay",
        "6",
        "--exit-code",
    ]);
    assert!(ok, "verify failed: {err}");
    let audit: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(audit["ok"], true, "{stdout}");
    assert_eq!(audit["replays"].as_array().unwrap().len(), 6);

    // Same method, same seeds: the A/A arm's counts are identical.
    let (ok, stdout, err) = ecbench(&[
        "compare",
        "--dir",
        out.to_str().unwrap(),
        "--a",
        "rho",
        "--b",
        "rho-aa",
        "--json",
    ]);
    assert!(ok, "{err}");
    let c: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(c["ops"]["ratio_b_over_a"], 1.0);
    assert_eq!(c["curves"].as_array().unwrap().len(), 2);
    // At L0 a wall-clock figure is never admitted, whatever was required.
    let (ok, stdout, _) = ecbench(&[
        "compare",
        "--dir",
        out.to_str().unwrap(),
        "--a",
        "rho",
        "--b",
        "bsgs",
        "--json",
        "--save",
    ]);
    assert!(ok);
    let c: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(c["ops"]["status"], "ok");
    assert!(c["ops"]["ci95"].is_array());

    // Tampering with a record is caught.
    let tampered = dir.join("t");
    copy_dir(&out, &tampered);
    let recs = std::fs::read_to_string(tampered.join("records.jsonl")).unwrap();
    std::fs::write(
        tampered.join("records.jsonl"),
        recs.replacen("\"status\":\"verified\"", "\"status\":\"wrong_answer\"", 1),
    )
    .unwrap();
    let (ok, stdout, _) = ecbench(&["verify", "--dir", tampered.to_str().unwrap(), "--exit-code"]);
    assert!(!ok);
    let audit: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(audit["ok"], false);

    // The SQL loads: schema first, every statement idempotent.
    let (ok, sql, err) = ecbench(&["db", "sql", out.to_str().unwrap()]);
    assert!(ok, "{err}");
    assert!(sql.starts_with("-- ecbench database schema"));
    assert!(sql.contains("INSERT INTO runs VALUES"));
    assert!(sql.contains("ON CONFLICT (record_id) DO NOTHING"));
    assert_eq!(sql.matches("INSERT INTO runs VALUES").count(), 48);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn bad_specs_are_refused_before_anything_runs() {
    let dir = scratch("refuse");
    let bad = SPEC.replace(
        "\"id\": \"bsgs.negation\"",
        "\"id\": \"rho.signed_frobenius\"",
    );
    let spec = dir.join("spec.json");
    std::fs::write(&spec, bad).unwrap();
    let out = dir.join("s");
    let (ok, _, err) = ecbench(&[
        "run",
        "--spec",
        spec.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--cpus",
        "none",
        "--lock",
        &lock(&dir),
    ]);
    assert!(!ok && err.contains("Koblitz curves only"), "{err}");
    assert!(!out.exists(), "a refused spec must leave nothing behind");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_factor_base_dump_carries_its_points() {
    let (ok, stdout, err) = ecbench(&[
        "fb",
        "--curve",
        r#"{"kind":"prime_search","bits":16,"seed":59297}"#,
        "--factor-base",
        "prime-abscissa:size=8",
    ]);
    assert!(ok, "{err}");
    let d: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(d["schema"], "ecbench.factor_base_dump/v1");
    assert_eq!(d["factor_base"]["columns"], 8);
    assert_eq!(
        d["points"].as_array().unwrap().len() as u64,
        d["factor_base"]["signed_points"].as_u64().unwrap()
    );
    assert!(d["factor_base"]["fb_id"]
        .as_str()
        .unwrap()
        .starts_with("FB1h"));
}

fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for e in std::fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let p = e.path();
        if p.is_dir() {
            copy_dir(&p, &to.join(e.file_name()));
        } else {
            std::fs::copy(&p, to.join(e.file_name())).unwrap();
        }
    }
}
