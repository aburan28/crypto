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

#[test]
fn online_windows_split_into_exclusive_claim_phases() {
    // The strong single-target reference and an index-calculus pipeline on
    // a registered Koblitz curve: every run verifies, every record carries
    // an online window whose exclusive phases sum to its wall time under
    // the claim schema's names, and the hardware counts are recorded or
    // their absence explained.
    let dir = scratch("online");
    let spec = r#"{
      "schema": "ecbench.spec/v1",
      "label": "online window test",
      "workloads": {"curves": [{"kind": "koblitz", "a": 1, "n": 17}],
                    "targets_per_curve": 2, "target_seed": 11},
      "arms": [
        {"name": "rho-strong", "role": "reference", "method": {"id": "rho.signed_frobenius_strong"}},
        {"name": "ic", "role": "candidate", "method": {"id": "ic.pipeline", "params": {
          "factor_base": "koblitz-orbit:divisor=0;1", "oracle": "mitm-frobenius:m=2"}}}
      ],
      "measurement": {"rounds": 2, "warmup": 0, "seed": 5, "isolation_required": "L0", "timeout_seconds": 120}
    }"#;
    std::fs::write(dir.join("spec.json"), spec).unwrap();
    let out = dir.join("s");
    let (ok, err) = run_spec(&dir, &out, &[]);
    assert!(ok, "{err}");
    let text = std::fs::read_to_string(out.join("records.jsonl")).unwrap();
    let mut saw_ic = false;
    for line in text.lines() {
        let r: Value = serde_json::from_str(line).unwrap();
        assert_eq!(r["outcome"]["status"], "verified", "{line}");
        let w = &r["online"];
        assert!(w.is_object(), "no online window: {line}");
        let wall = w["wall_ns"].as_u64().unwrap();
        let sum: u64 = w["phases_ns"]
            .as_object()
            .unwrap()
            .values()
            .map(|v| v.as_u64().unwrap())
            .sum();
        assert_eq!(sum, wall, "exclusive phases must sum to the window");
        let stages: Vec<&str> = w["included_stages"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        if r["arm"] == "ic" {
            saw_ic = true;
            for s in [
                "target_query",
                "target_PDP",
                "target_relation_check",
                "target_descent",
                "target_recovery_check",
            ] {
                assert!(stages.contains(&s), "{s} missing from {stages:?}");
                assert!(w["phases_ns"][s].is_u64(), "{s} has no time");
            }
            assert!(w["phases_ns"]["target_PDP"].as_u64().unwrap() > 0);
        } else {
            assert_eq!(stages, ["walk", "collision", "recovery_check"]);
            assert_eq!(r["method"]["id"], "rho.signed_frobenius_strong");
        }
        let hw = &r["time"]["hw_solve"];
        assert!(
            hw["instructions"].is_u64() || hw["error"].is_string(),
            "hardware counts neither read nor explained: {hw}"
        );
    }
    assert!(saw_ic);
    let (ok, _, err) = ecbench(&[
        "verify",
        "--dir",
        out.to_str().unwrap(),
        "--replay-all",
        "--exit-code",
    ]);
    assert!(ok, "{err}");
    // These targets are planted: no claim can be built from them.
    let w = first_workload(&dir.join("spec.json"));
    let (ok, _, err) = ecbench(&[
        "claim",
        "build",
        "--dir",
        out.to_str().unwrap(),
        "--ic",
        "ic",
        "--rho",
        "rho-strong",
        "--workload",
        &w,
    ]);
    assert!(!ok && err.contains("planted"), "{err}");
    let _ = std::fs::remove_dir_all(&dir);
}

fn first_workload(spec: &Path) -> String {
    let (ok, out, err) = ecbench(&["plan", "--spec", spec.to_str().unwrap(), "--json"]);
    assert!(ok, "{err}");
    let plan: Value = serde_json::from_str(&out).unwrap();
    plan["workloads"][0]["workload_id"]
        .as_str()
        .unwrap()
        .to_string()
}

#[test]
fn a_public_target_session_yields_a_checked_vs_rho_claim() {
    // One public target, the strong rho reference, an IC pipeline and the
    // plain signed-Frobenius rho: the claim builds from the IC and strong
    // rho runs, the checker fails it on the independent replay alone, a
    // receipt from the session's own host class is refused, and a receipt
    // from another host class (simulated by its recorded class) passes.
    let dir = scratch("claim");
    let spec = r#"{
      "schema": "ecbench.spec/v1",
      "label": "claim test",
      "workloads": {"curves": [{"kind": "koblitz", "a": 1, "n": 17}],
                    "targets_per_curve": 1, "target_seed": 11, "target_kind": "public"},
      "arms": [
        {"name": "rho-strong", "role": "reference", "method": {"id": "rho.signed_frobenius_strong"}},
        {"name": "rho", "role": "baseline", "method": {"id": "rho.signed_frobenius"}},
        {"name": "ic", "role": "candidate", "method": {"id": "ic.pipeline", "params": {
          "factor_base": "koblitz-orbit:divisor=0;1", "oracle": "mitm-frobenius:m=2"}}}
      ],
      "measurement": {"rounds": 1, "warmup": 0, "seed": 5, "isolation_required": "L0", "timeout_seconds": 120}
    }"#;
    std::fs::write(dir.join("spec.json"), spec).unwrap();
    let out = dir.join("s");
    let (ok, err) = run_spec(&dir, &out, &[]);
    assert!(ok, "{err}");
    let w = first_workload(&dir.join("spec.json"));
    let s = out.to_str().unwrap();
    let claim_path = dir.join("claim.json");
    let build = |rho: &str, extra: &[&str]| {
        let mut args = vec![
            "claim",
            "build",
            "--dir",
            s,
            "--ic",
            "ic",
            "--rho",
            rho,
            "--workload",
            &w,
            "--out",
            claim_path.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        ecbench(&args)
    };
    let check = || {
        let (ok, out, err) = ecbench(&["claim", "check", "--report", claim_path.to_str().unwrap()]);
        (ok, serde_json::from_str::<Value>(&out).unwrap(), err)
    };

    // The plain walk is not the strong reference.
    let (ok, _, err) = build("rho", &[]);
    assert!(
        !ok && err.contains("strong single-target reference"),
        "{err}"
    );

    // Without an independent replay: built, and failed on exactly that.
    let (ok, _, err) = build("rho-strong", &[]);
    assert!(ok, "{err}");
    let claim: Value =
        serde_json::from_str(&std::fs::read_to_string(&claim_path).unwrap()).unwrap();
    assert!(claim["candidate_id"]
        .as_str()
        .unwrap()
        .starts_with("IC1N17Ckb1fb442PDP2mitmfrobeniusRCwalkLAincrementalgaussTDjointISO0h"));
    assert_eq!(claim["target_count"], 1);
    assert_eq!(claim["ic_target_hash"], claim["rho_target_hash"]);
    assert_eq!(
        claim["ic_resource_envelope"],
        claim["rho_resource_envelope"]
    );
    let (ok, c, _) = check();
    assert!(!ok);
    assert_eq!(c["status"], "FAIL");
    assert_eq!(
        c["validation_errors"],
        serde_json::json!([
            "independent_validation must be true",
            "ic_replay_certificate_sha256 must be a full lowercase SHA-256 digest",
            "rho_replay_certificate_sha256 must be a full lowercase SHA-256 digest"
        ])
    );
    assert_eq!(
        c["missing_stage_fields"],
        serde_json::json!([
            "ic_replay_certificate_sha256",
            "rho_replay_certificate_sha256"
        ])
    );
    assert_eq!(c["missing_global_provenance"], serde_json::json!([]));

    // A replay on the session's own host class is not independent.
    let receipt = dir.join("receipt.json");
    let (ok, _, err) = ecbench(&[
        "verify",
        "--dir",
        s,
        "--replay-all",
        "--exit-code",
        "--out",
        receipt.to_str().unwrap(),
    ]);
    assert!(ok, "{err}");
    let (ok, _, err) = build(
        "rho-strong",
        &[
            "--independent-receipt",
            receipt.to_str().unwrap(),
            "--pointer",
            "test",
        ],
    );
    assert!(!ok && err.contains("own host class"), "{err}");

    // The same receipt as another host class would record it.  The
    // receipt is the auditor's word: a claim cites its digest and pointer
    // so a reader can fetch it and rerun the audit.
    let mut r: Value = serde_json::from_str(&std::fs::read_to_string(&receipt).unwrap()).unwrap();
    r["auditor_env_class_id"] = Value::from("ECBENV2h000000000000");
    let elsewhere = dir.join("receipt-elsewhere.json");
    std::fs::write(&elsewhere, serde_json::to_string_pretty(&r).unwrap()).unwrap();
    let (ok, _, err) = build(
        "rho-strong",
        &[
            "--independent-receipt",
            elsewhere.to_str().unwrap(),
            "--pointer",
            "test",
            "--exit-code",
        ],
    );
    assert!(ok, "{err}");
    let (ok, c, err) = check();
    assert!(ok, "{c} {err}");
    assert_eq!(c["status"], "PASS");

    // The claim loads into the database beside its session, checked again.
    let (ok, sql, err) = ecbench(&["db", "sql", s, claim_path.to_str().unwrap()]);
    assert!(ok, "{err}");
    assert!(
        sql.contains("INSERT OR REPLACE INTO claims VALUES ("),
        "{sql}"
    );
    assert!(sql.contains("'PASS'"));

    // A receipt for other bytes is refused.
    let mut r2 = r.clone();
    r2["files"]["records.jsonl"] = Value::from("0".repeat(64));
    std::fs::write(&elsewhere, serde_json::to_string_pretty(&r2).unwrap()).unwrap();
    let (ok, _, err) = build(
        "rho-strong",
        &[
            "--independent-receipt",
            elsewhere.to_str().unwrap(),
            "--pointer",
            "test",
        ],
    );
    assert!(!ok && err.contains("records.jsonl"), "{err}");
    let _ = std::fs::remove_dir_all(&dir);
}

fn sameness_spec(curve: &str) -> String {
    format!(
        r#"{{
      "schema": "ecbench.spec/v1",
      "label": "wide sameness",
      "workloads": {{"curves": [{curve}], "targets_per_curve": 2, "target_seed": 3}},
      "arms": [
        {{"name": "rho-frob", "role": "reference", "method": {{"id": "rho.signed_frobenius"}}}},
        {{"name": "budget", "role": "candidate", "method": {{"id": "rho.signed_frobenius_budget", "params": {{"steps": "20000"}}}}}}
      ],
      "measurement": {{"rounds": 1, "warmup": 0, "seed": 9, "isolation_required": "L0", "timeout_seconds": 120}}
    }}"#
    )
}

fn records_of(session: &Path) -> Vec<Value> {
    std::fs::read_to_string(session.join("records.jsonl"))
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}

#[test]
fn a_two_word_session_of_a_word_size_curve_matches_the_narrow_session() {
    // The same curve built narrow (`koblitz`) and wide (`koblitz_explicit`
    // from the narrow plan's facts): the same workload ids, and for every
    // arm the same status, answer, cost, counters and phases, record for
    // record.  The walk is the reference on both; the budgeted walk
    // exhausts on both.
    let dir = scratch("wide-narrow");
    std::fs::write(
        dir.join("spec.json"),
        sameness_spec(r#"{"kind": "koblitz", "a": 0, "n": 41}"#),
    )
    .unwrap();
    let (ok, out, err) = ecbench(&[
        "plan",
        "--spec",
        dir.join("spec.json").to_str().unwrap(),
        "--json",
    ]);
    assert!(ok, "{err}");
    let plan: Value = serde_json::from_str(&out).unwrap();
    let c = &plan["workloads"][0]["curve"];
    let wide_curve = format!(
        r#"{{"kind": "koblitz_explicit", "a": 0, "n": 41, "modulus": "{}", "r": "{}", "cofactor": "{}", "gx": "{}", "gy": "{}"}}"#,
        c["modulus"].as_str().unwrap(),
        c["r"],
        c["cofactor"],
        c["generator"][0].as_str().unwrap(),
        c["generator"][1].as_str().unwrap()
    );
    let wide = scratch("wide-wide");
    std::fs::write(wide.join("spec.json"), sameness_spec(&wide_curve)).unwrap();
    let (ok, err) = run_spec(&dir, &dir.join("s"), &[]);
    assert!(ok, "{err}");
    let (ok, err) = run_spec(&wide, &wide.join("s"), &[]);
    assert!(ok, "{err}");
    let (a, b) = (records_of(&dir.join("s")), records_of(&wide.join("s")));
    assert_eq!(a.len(), 4);
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        let tag = format!("{} on {}", x["arm"], x["workload"]["workload_id"]);
        assert_eq!(
            x["workload"]["workload_id"], y["workload"]["workload_id"],
            "{tag}"
        );
        assert_eq!(x["arm"], y["arm"], "{tag}");
        assert_eq!(x["outcome"]["status"], y["outcome"]["status"], "{tag}");
        assert_eq!(
            x["outcome"]["recovered"], y["outcome"]["recovered"],
            "{tag}"
        );
        assert_eq!(x["cost"]["total_gae"], y["cost"]["total_gae"], "{tag}");
        assert!(
            x.get("counters").is_some() && x.get("phases").is_some(),
            "{x}"
        );
        assert_eq!(x["counters"], y["counters"], "{tag}");
        assert_eq!(x["phases"], y["phases"], "{tag}");
    }
    assert!(a
        .iter()
        .filter(|r| r["arm"] == "rho-frob")
        .all(|r| r["outcome"]["status"] == "verified"));
    assert!(a
        .iter()
        .filter(|r| r["arm"] == "budget")
        .all(|r| r["outcome"]["status"] == "exhausted"));
    for d in [&dir, &wide] {
        let (ok, _, err) = ecbench(&[
            "verify",
            "--dir",
            d.join("s").to_str().unwrap(),
            "--replay-all",
            "--exit-code",
        ]);
        assert!(ok, "{err}");
    }
    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::remove_dir_all(&wide);
}

#[test]
fn the_m83_gate_curve_plans_runs_a_budgeted_walk_and_replays() {
    // The committed step-rate spec: the registered gate slug with
    // two-word facts, every run `exhausted` at the budget with the walk's
    // shape recorded, exact replays, and a database load that carries the
    // subgroup order as text.
    let dir = scratch("m83");
    let spec =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/ecbench/specs/gate-m83-steprate.json");
    let (ok, out, err) = ecbench(&["plan", "--spec", spec.to_str().unwrap(), "--json"]);
    assert!(ok, "{err}");
    let plan: Value = serde_json::from_str(&out).unwrap();
    let c = &plan["workloads"][0]["curve"];
    assert_eq!(c["slug"], "icv1-f2m83-tm6151469093347-debefd74");
    assert_eq!(c["registered"], true);
    assert_eq!(c["r"], "2417851639230796216685689");
    assert_eq!(c["cofactor"], 4);
    assert_eq!(c["automorphisms_available"], 166);
    std::fs::copy(&spec, dir.join("spec.json")).unwrap();
    let out_dir = dir.join("s");
    let (ok, err) = run_spec(&dir, &out_dir, &[]);
    assert!(ok, "{err}");
    let records = records_of(&out_dir);
    assert!(records.len() >= 4, "{}", records.len());
    for r in &records {
        assert_eq!(r["outcome"]["status"], "exhausted", "{r}");
        assert!(r["counters"]["walk_operations"].as_u64().unwrap() >= 2_000_000);
        assert_eq!(r["counters"]["distinguished_point_bits"], 20);
        assert_eq!(r["counters"]["jumps"], 16);
        assert!(
            r["cost"]["s"].as_f64().unwrap() < 1e-4,
            "{}",
            r["cost"]["s"]
        );
        assert_eq!(r["workload"]["curve"]["r"], "2417851639230796216685689");
    }
    let audit = dir.join("audit.json");
    let (ok, _, err) = ecbench(&[
        "verify",
        "--dir",
        out_dir.to_str().unwrap(),
        "--replay-all",
        "--out",
        audit.to_str().unwrap(),
        "--exit-code",
    ]);
    assert!(ok, "{err}");
    let receipt: Value = serde_json::from_str(&std::fs::read_to_string(&audit).unwrap()).unwrap();
    assert_eq!(receipt["ok"], true);
    assert_eq!(
        receipt["replays"].as_array().unwrap().len(),
        records.iter().filter(|r| r["warmup"] == false).count()
    );
    assert!(receipt["replays"]
        .as_array()
        .unwrap()
        .iter()
        .all(|r| r["reproduced"] == true));
    let (ok, sql, err) = ecbench(&["db", "sql", out_dir.to_str().unwrap()]);
    assert!(ok, "{err}");
    assert!(
        sql.contains("'2417851639230796216685689'"),
        "{}",
        &sql[..2000]
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
