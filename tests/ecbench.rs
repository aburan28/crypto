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

fn run_spec(dir: &Path, out: &Path, extra: &[&str]) -> (bool, String) {
    let spec = dir.join("spec.json");
    if !spec.exists() {
        std::fs::write(&spec, SPEC).unwrap();
    }
    let mut args = vec![
        "run",
        "--spec",
        spec.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--cpus",
        "none",
        "--lock",
        Box::leak(lock(dir).into_boxed_str()),
        "--quiet",
    ];
    args.extend_from_slice(extra);
    let (ok, _, err) = ecbench(&args);
    (ok, err)
}

#[test]
fn two_sessions_queued_for_one_directory_cannot_both_write_it() {
    // Both start before either holds the lock; the second waits, then must
    // fail on the atomic mkdir instead of writing over the first.
    let dir = scratch("race");
    std::fs::write(dir.join("spec.json"), SPEC).unwrap();
    let out = dir.join("s");
    let spawn = || {
        Command::new(exe())
            .args([
                "run",
                "--spec",
                dir.join("spec.json").to_str().unwrap(),
                "--out",
                out.to_str().unwrap(),
                "--cpus",
                "none",
                "--lock",
                &lock(&dir),
                "--wait",
                "--quiet",
            ])
            .spawn()
            .unwrap()
    };
    let (mut a, mut b) = (spawn(), spawn());
    let (sa, sb) = (a.wait().unwrap(), b.wait().unwrap());
    assert!(
        sa.success() != sb.success(),
        "exactly one of the two queued sessions may succeed"
    );
    let (ok, stdout, err) = ecbench(&["verify", "--dir", out.to_str().unwrap(), "--exit-code"]);
    assert!(ok, "the surviving session must audit clean: {err} {stdout}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_candidate_session_compares_against_a_baseline_session() {
    // Two sessions of one spec (a baseline and a candidate binary, here
    // the same one): pairs match by (workload, round) across sessions, the
    // ratio is exactly 1 and the interval exists.
    let dir = scratch("cross");
    let (a, b) = (dir.join("a"), dir.join("b"));
    assert!(run_spec(&dir, &a, &[]).0);
    assert!(run_spec(&dir, &b, &[]).0);
    let (ok, stdout, err) = ecbench(&[
        "compare",
        "--dir",
        a.to_str().unwrap(),
        "--a",
        "bsgs",
        "--b-dir",
        b.to_str().unwrap(),
        "--b",
        "bsgs",
        "--json",
    ]);
    assert!(ok, "{err}");
    let c: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(c["ops"]["ratio_b_over_a"], 1.0);
    assert_eq!(c["ops"]["ci95"], serde_json::json!([1.0, 1.0]));
    assert_eq!(c["ops"]["same_seeds"], true);
    assert_eq!(c["ops"]["status"], "ok");
    assert_eq!(c["wall"]["status"], "refused");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_resealed_hand_edited_figure_is_caught() {
    use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
    use crypto_lib::cryptanalysis::ecbench::record::Record;
    use crypto_lib::cryptanalysis::ecbench::runner::Session;
    let dir = scratch("forge");
    let out = dir.join("s");
    assert!(run_spec(&dir, &out, &[]).0);
    // Halve one verified run's S, reseal the record and rehash the file:
    // every hash agrees, but S no longer follows from the counts.
    let text = std::fs::read_to_string(out.join("records.jsonl")).unwrap();
    let mut lines: Vec<String> = text.lines().map(String::from).collect();
    let i = lines
        .iter()
        .position(|l| l.contains("\"warmup\":false") && l.contains("\"status\":\"verified\""))
        .unwrap();
    let mut rec: Record = serde_json::from_str(&lines[i]).unwrap();
    rec.cost.s = rec.cost.s.map(|s| s / 2.0);
    lines[i] = rec.seal();
    let body = lines.join("\n") + "\n";
    std::fs::write(out.join("records.jsonl"), &body).unwrap();
    let mut session: Session =
        serde_json::from_str(&std::fs::read_to_string(out.join("session.json")).unwrap()).unwrap();
    session.records_sha256 = Some(sha256_hex(body.as_bytes()));
    std::fs::write(
        out.join("session.json"),
        serde_json::to_string_pretty(&session).unwrap() + "\n",
    )
    .unwrap();
    let (ok, stdout, _) = ecbench(&["verify", "--dir", out.to_str().unwrap(), "--exit-code"]);
    assert!(!ok);
    assert!(
        stdout.contains("does not follow from the record"),
        "{stdout}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_database_refuses_a_slug_with_another_identity() {
    let Ok(probe) = Command::new("sqlite3").arg("-version").output() else {
        eprintln!("sqlite3 not installed; skipped");
        return;
    };
    if !probe.status.success() {
        return;
    }
    let dir = scratch("db");
    let out = dir.join("s");
    assert!(run_spec(&dir, &out, &[]).0);
    let (ok, sql, err) = ecbench(&["db", "sql", out.to_str().unwrap()]);
    assert!(ok, "{err}");
    let db = dir.join("t.db");
    let load = |sql: &str| {
        let mut child = Command::new("sqlite3")
            .args(["-bail", db.to_str().unwrap()])
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .unwrap();
        use std::io::Write;
        child
            .stdin
            .take()
            .unwrap()
            .write_all(sql.as_bytes())
            .unwrap();
        child.wait_with_output().unwrap().status.success()
    };
    assert!(load(&sql), "the first load must succeed");
    assert!(
        load(&sql),
        "loading the same session again must change nothing and succeed"
    );
    // The same slug arriving on another ICV1 string must be an error.
    let slug = "icv1-fp16-t295-8d3c3165";
    assert!(sql.contains(slug));
    let forged = format!(
        "INSERT INTO curves VALUES ('{slug}', 'ICV1:forged', 'prime', NULL, 16, '1', '1', '1', 1.0, 2, 1.0, 1) ON CONFLICT (slug) DO NOTHING;"
    );
    assert!(
        !load(&forged),
        "a conflicting identity must not be silently dropped"
    );
    let _ = std::fs::remove_dir_all(&dir);
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
