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

    // The resource vector includes every measured attempt, excludes the
    // warm-ups, and preserves native work in its original units.
    let (ok, stdout, err) = ecbench(&["resources", "--dir", out.to_str().unwrap()]);
    assert!(ok, "resources failed: {err}");
    let resources: Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(resources["schema"], "ecbench.resources/v1");
    let arms = resources["arms"].as_array().unwrap();
    assert_eq!(arms.len(), 4);
    for arm in arms {
        assert_eq!(arm["measured"], 8);
        assert_eq!(arm["verified"], 8);
        assert!(arm["process_wall_ns_sum"].as_u64().unwrap() > 0);
        assert!(arm["peak_rss_kib"].as_u64().unwrap() > 0);
        assert!(arm["method_counter_totals"].is_object());
    }

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
fn counted_ic_sat_replays_with_different_solver_wall_times() {
    let dir = scratch("counted-sat");
    let spec = dir.join("spec.json");
    std::fs::write(
        &spec,
        r#"{
  "schema":"ecbench.spec/v1",
  "label":"counted SAT replay regression",
  "workloads":{"curves":[{"kind":"koblitz","a":1,"n":17}],"targets_per_curve":2,"target_seed":41,"target_kind":"public"},
  "arms":[
    {"name":"rho","role":"reference","method":{"id":"rho.signed_frobenius_strong"}},
    {"name":"sat-cnf","role":"candidate","method":{"id":"ic.pipeline_counted","params":{"factor_base":"koblitz-orbit:divisor=0;1","oracle":"descent-algebraic:m=2","solver":"sat-cdcl:sat_conflict_budget=20000,sat_xor_encoding=cnf","max_trials":"20000"}}},
    {"name":"sat-xor","role":"candidate","method":{"id":"ic.pipeline_counted","params":{"factor_base":"koblitz-orbit:divisor=0;1","oracle":"descent-algebraic:m=2","solver":"sat-cdcl:sat_conflict_budget=20000,sat_xor_encoding=native","max_trials":"20000"}}}
  ],
  "measurement":{"rounds":1,"warmup":0,"seed":9123401,"isolation_required":"L0","timeout_seconds":60}
}"#,
    )
    .unwrap();
    let out = dir.join("session");
    let (ok, err) = run_spec(&dir, &out, &[]);
    assert!(ok, "{err}");
    let records: Vec<Value> = std::fs::read_to_string(out.join("records.jsonl"))
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(records.len(), 6);
    for r in records.iter().filter(|r| r["method"]["family"] == "ic") {
        assert_eq!(r["outcome"]["status"], "verified");
        assert_eq!(r["cost"]["lower_bound"], true);
        assert!(r["cost"]["unpriced"]
            .as_array()
            .unwrap()
            .iter()
            .any(|u| u == "solver_conflicts_uncharged"));
    }
    let receipt = dir.join("audit.json");
    let (ok, _, err) = ecbench(&[
        "verify",
        "--dir",
        out.to_str().unwrap(),
        "--replay-all",
        "--out",
        receipt.to_str().unwrap(),
        "--exit-code",
    ]);
    assert!(ok, "{err}");
    let audit: Value = serde_json::from_str(&std::fs::read_to_string(receipt).unwrap()).unwrap();
    assert_eq!(audit["replays"].as_array().unwrap().len(), 6);
    assert!(audit["replays"]
        .as_array()
        .unwrap()
        .iter()
        .all(|r| r["reproduced"] == true));
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

#[test]
fn pdp3_koblitz_f4_and_f6_verify_on_one_query_stream() {
    // #1333's inherited-F4 and F6-IC decomposers inside ic.pipeline on the
    // E_0 m = 13 curve: every run verifies; both arms, being exact solvers
    // on the same base and seed, issue the same queries and find the same
    // relations; F6-IC's gate additions are charged and the solver's word
    // XORs are counted but left out of S, which is therefore a lower bound.
    let dir = scratch("pdp3");
    let spec = r#"{
      "schema": "ecbench.spec/v1",
      "label": "pdp3-koblitz test",
      "workloads": {"curves": [{"kind": "koblitz", "a": 0, "n": 13}],
                    "targets_per_curve": 2, "target_seed": 7, "target_kind": "public"},
      "arms": [
        {"name": "ic-f4", "role": "baseline", "method": {"id": "ic.pipeline", "params": {
          "factor_base": "koblitz-standard-subspace:dimension=5",
          "oracle": "pdp3-koblitz:m=3,engine=inherited-f4,degree=3,node_budget=8192"}}},
        {"name": "ic-f6", "role": "candidate", "method": {"id": "ic.pipeline", "params": {
          "factor_base": "koblitz-standard-subspace:dimension=5",
          "oracle": "pdp3-koblitz:m=3,engine=f6-ic,degree=3,node_budget=8192"}}}
      ],
      "measurement": {"rounds": 1, "warmup": 0, "seed": 5, "isolation_required": "L0", "timeout_seconds": 120}
    }"#;
    std::fs::write(dir.join("spec.json"), spec).unwrap();
    let out = dir.join("s");
    let (ok, err) = run_spec(&dir, &out, &[]);
    assert!(ok, "{err}");
    let text = std::fs::read_to_string(out.join("records.jsonl")).unwrap();
    let mut by_workload: std::collections::BTreeMap<String, Vec<Value>> = Default::default();
    for line in text.lines() {
        let r: Value = serde_json::from_str(line).unwrap();
        assert_eq!(r["outcome"]["status"], "verified", "{line}");
        assert_eq!(r["cost"]["lower_bound"], true, "{line}");
        let unpriced = r["cost"]["unpriced"].to_string();
        assert!(unpriced.contains("word_XORs"), "{unpriced}");
        let extra = &r["solver"]["extra"];
        let geometry = extra["geometric_group_additions"].as_u64().unwrap();
        assert_eq!(extra["geometric_fallbacks"], 0, "{line}");
        if r["arm"] == "ic-f6" {
            assert!(geometry > 0, "F6-IC charged no geometry: {line}");
        } else {
            assert_eq!(geometry, 0, "{line}");
        }
        assert!(r["solver"]["ops"].as_u64().unwrap() > 0, "{line}");
        let w = r["workload"]["workload_id"].as_str().unwrap().to_string();
        by_workload.entry(w).or_default().push(r);
    }
    assert_eq!(by_workload.len(), 2);
    for (w, rs) in by_workload {
        assert_eq!(rs.len(), 2, "{w}");
        for key in ["targets_tried", "relations_found", "matrix_rank"] {
            assert_eq!(
                rs[0]["counters"][key], rs[1]["counters"][key],
                "{w}: the arms' query streams diverged at {key}"
            );
        }
    }
    let (ok, _, err) = ecbench(&[
        "verify",
        "--dir",
        out.to_str().unwrap(),
        "--replay-all",
        "--exit-code",
    ]);
    assert!(ok, "replay: {err}");
}

const BOUNDS_SPEC: &str = r#"{
  "schema": "ecbench.spec/v1",
  "label": "bounds integration test: four sizes",
  "workloads": {
    "curves": [{"kind": "prime_search", "bits": 12, "seed": 59297},
               {"kind": "prime_search", "bits": 14, "seed": 59297},
               {"kind": "prime_search", "bits": 16, "seed": 59297},
               {"kind": "prime_search", "bits": 18, "seed": 59297}],
    "targets_per_curve": 4,
    "target_seed": 11
  },
  "arms": [
    {"name": "rho-plain", "role": "reference", "method": {"id": "rho.plain"}},
    {"name": "rho-neg", "role": "candidate", "method": {"id": "rho.negation"}},
    {"name": "bsgs-neg", "role": "candidate", "method": {"id": "bsgs.negation"}},
    {"name": "rho-plain-aa", "role": "control", "method": {"id": "rho.plain"}}
  ],
  "measurement": {"rounds": 2, "warmup": 0, "seed": 5, "isolation_required": "L0", "timeout_seconds": 60}
}"#;

fn read_json(p: &Path) -> Value {
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

#[test]
fn bounds_a_frontier_and_challenge_verdicts_round_trip() {
    let dir = scratch("bounds");
    let s = |p: &Path| p.to_str().unwrap().to_string();
    let spec = dir.join("spec.json");
    std::fs::write(&spec, BOUNDS_SPEC).unwrap();
    let session = dir.join("s0");
    let (ok, _, err) = ecbench(&[
        "run",
        "--spec",
        &s(&spec),
        "--out",
        &s(&session),
        "--cpus",
        "none",
        "--lock",
        &lock(&dir),
        "--quiet",
    ]);
    assert!(ok, "run failed: {err}");

    // Bounds of three arms, sealed, re-derivable, and refused when touched.
    let records = dir.join("records");
    std::fs::create_dir_all(&records).unwrap();
    let mut ids = std::collections::BTreeMap::new();
    for arm in ["rho-plain", "rho-neg", "bsgs-neg"] {
        let out = records.join(format!("{arm}.json"));
        let (ok, _, err) = ecbench(&[
            "bound",
            "fit",
            "--dir",
            &s(&session),
            "--arm",
            arm,
            "--root",
            &s(&dir),
            "--out",
            &s(&out),
        ]);
        assert!(ok, "bound fit {arm}: {err}");
        let b = read_json(&out);
        assert_eq!(b["schema"], "ecbench.bound/v1");
        let id = b["bound_id"].as_str().unwrap().to_string();
        assert!(id.starts_with("ECBND1h") && id.len() == 7 + 12, "{id}");
        assert_eq!(b["domain"]["family"], "prime");
        assert_eq!(b["domain"]["tier"], "toy");
        assert_eq!(b["fit"]["sizes"], 4);
        assert_eq!(b["fit"]["scaling_claim"], true);
        assert_eq!(b["fit"]["declared_alpha"], 0.5);
        assert_eq!(b["level"], "exponent");
        assert_eq!(b["admissibility"]["status"], "admissible");
        assert_eq!(b["dimensions"]["memory"]["known"], true);
        assert_eq!(b["provenance"]["verified"], 32);
        assert_eq!(b["provenance"]["sessions"][0]["dir"], "s0");
        ids.insert(arm, id);
    }
    let (ok, _, err) = ecbench(&[
        "bound",
        "check",
        "--root",
        &s(&dir),
        "--record",
        &s(&records.join("rho-plain.json")),
        &s(&records.join("rho-neg.json")),
    ]);
    assert!(ok, "bound check: {err}");
    let touched = dir.join("touched.json");
    let text = std::fs::read_to_string(records.join("rho-neg.json")).unwrap();
    std::fs::write(
        &touched,
        text.replacen("\"verified\": 32", "\"verified\": 33", 1),
    )
    .unwrap();
    let (ok, _, err) = ecbench(&[
        "bound",
        "check",
        "--root",
        &s(&dir),
        "--record",
        &s(&touched),
    ]);
    assert!(!ok && err.contains("seal mismatch"), "{err}");

    // The frontier: BSGS leads on ops, the rho walks on memory; the page
    // is current right after it is built and stale once edited.
    let fjson = dir.join("frontier.json");
    let fmd = dir.join("FRONTIER.md");
    let (ok, _, err) = ecbench(&[
        "frontier",
        "build",
        "--bounds",
        &s(&records),
        "--out",
        &s(&fjson),
        "--markdown",
        &s(&fmd),
    ]);
    assert!(ok, "frontier build: {err}");
    let f = read_json(&fjson);
    assert_eq!(f["schema"], "ecbench.frontier/v1");
    assert_eq!(f["domains"].as_array().unwrap().len(), 1);
    let d = &f["domains"][0];
    assert_eq!(d["entries"].as_array().unwrap().len(), 3);
    assert_eq!(d["ops_leader"], ids["bsgs-neg"]);
    let frontier_ids: Vec<&str> = d["entries"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|e| e["is_frontier"] == true)
        .map(|e| e["bound_id"].as_str().unwrap())
        .collect();
    assert!(
        frontier_ids.contains(&ids["bsgs-neg"].as_str()),
        "{frontier_ids:?}"
    );
    let (ok, _, err) = ecbench(&[
        "frontier",
        "build",
        "--bounds",
        &s(&records),
        "--out",
        &s(&fjson),
        "--markdown",
        &s(&fmd),
        "--check",
    ]);
    assert!(ok, "frontier check: {err}");
    std::fs::write(&fmd, "edited\n").unwrap();
    let (ok, _, err) = ecbench(&[
        "frontier",
        "build",
        "--bounds",
        &s(&records),
        "--out",
        &s(&fjson),
        "--markdown",
        &s(&fmd),
        "--check",
    ]);
    assert!(!ok && err.contains("stale"), "{err}");

    // A challenge against rho.plain, holding the bound fitted above.
    let draft = dir.join("draft.json");
    std::fs::write(
        &draft,
        format!(
            r#"{{
  "label": "test: beat rho.plain",
  "domain": {{"problem": "ecdlp.single_target", "family": "prime", "target_kind": "planted",
             "unit": "ecbench.gae", "tier": "toy",
             "envelope": {{"targets": 1, "precomputation": "none", "threads": 1}}}},
  "incumbent": {{"bound_id": "{}", "method": {{"id": "rho.plain"}}}},
  "workloads": {{"curves": [{{"kind": "prime_search", "bits": 12, "seed": 59297}},
                           {{"kind": "prime_search", "bits": 14, "seed": 59297}},
                           {{"kind": "prime_search", "bits": 16, "seed": 59297}},
                           {{"kind": "prime_search", "bits": 18, "seed": 59297}}],
                "targets_per_curve": 4, "nonce": 7}},
  "measurement": {{"rounds": 2, "warmup": 0, "isolation_required": "L0", "timeout_seconds": 60}},
  "acceptance": {{"min_runs_per_size": 8}}
}}"#,
            ids["rho-plain"]
        ),
    )
    .unwrap();
    let challenge = dir.join("challenge.json");
    let (ok, _, err) = ecbench(&[
        "challenge",
        "seal",
        "--draft",
        &s(&draft),
        "--out",
        &s(&challenge),
    ]);
    assert!(ok, "challenge seal: {err}");
    let (ok, _, err) = ecbench(&["challenge", "check", "--file", &s(&challenge)]);
    assert!(ok, "challenge check: {err}");
    let c = read_json(&challenge);
    assert!(c["challenge_id"].as_str().unwrap().starts_with("ECCH1h"));
    assert_eq!(
        c["acceptance"]["axes"],
        serde_json::json!(["ops", "memory"])
    );

    // Epoch 1: BSGS with the negation map.  Fewer operations, a √r table:
    // a trade, and a new bound that names the incumbent's.
    let spec1 = dir.join("spec1.json");
    let (ok, _, err) = ecbench(&[
        "challenge",
        "spec",
        "--challenge",
        &s(&challenge),
        "--candidate",
        r#"{"id":"bsgs.negation"}"#,
        "--epoch",
        "1",
        "--out",
        &s(&spec1),
    ]);
    assert!(ok, "challenge spec: {err}");
    let s1 = dir.join("s1");
    let (ok, _, err) = ecbench(&[
        "run",
        "--spec",
        &s(&spec1),
        "--out",
        &s(&s1),
        "--cpus",
        "none",
        "--lock",
        &lock(&dir),
        "--quiet",
    ]);
    assert!(ok, "run s1: {err}");
    let v1 = dir.join("verdict1.json");
    let b1 = dir.join("bound1.json");
    let (ok, _, err) = ecbench(&[
        "challenge",
        "verdict",
        "--challenge",
        &s(&challenge),
        "--dir",
        &s(&s1),
        "--epoch",
        "1",
        "--replay-all",
        "--bounds",
        &s(&records),
        "--root",
        &s(&dir),
        "--out",
        &s(&v1),
        "--bound-out",
        &s(&b1),
        "--exit-code",
    ]);
    assert!(ok, "verdict 1: {err}");
    let v = read_json(&v1);
    assert_eq!(v["schema"], "ecbench.verdict/v1");
    assert_eq!(v["outcome"], "trade", "{}", v["statement"]);
    assert_eq!(v["advances_on"], serde_json::json!(["ops"]));
    assert_eq!(v["regresses_on"], serde_json::json!(["memory"]));
    assert_eq!(v["level_moved"], Value::Null);
    assert_eq!(v["session"]["spec_matches"], true);
    assert_eq!(v["audit"]["ok"], true);
    assert_eq!(v["audit"]["replay_all"], true);
    assert_eq!(v["audit"]["replays"], v["audit"]["replays_reproduced"]);
    assert_eq!(v["axes"][0]["axis"], "ops");
    assert_eq!(v["axes"][0]["pairs"], 32);
    assert_eq!(v["control"]["ratio"], 1.0);
    assert!(v["incumbent_drift"]["recorded_bound_id"] == ids["rho-plain"]);
    assert_eq!(v["fits"]["candidate"]["sizes"], 4);
    let nb = read_json(&b1);
    assert_eq!(nb["improves_on"], serde_json::json!([ids["rho-plain"]]));
    assert_eq!(nb["verdict_id"], v["verdict_id"]);
    assert_eq!(nb["method"]["id"], "bsgs.negation");

    // Epoch 2: the incumbent against itself.  Same seeds, same counts:
    // the ratio is exactly 1, the outcome `matches`, no bound is written.
    let spec2 = dir.join("spec2.json");
    assert!(
        ecbench(&[
            "challenge",
            "spec",
            "--challenge",
            &s(&challenge),
            "--candidate",
            r#"{"id":"rho.plain"}"#,
            "--epoch",
            "2",
            "--out",
            &s(&spec2),
        ])
        .0
    );
    let s2 = dir.join("s2");
    assert!(
        ecbench(&[
            "run",
            "--spec",
            &s(&spec2),
            "--out",
            &s(&s2),
            "--cpus",
            "none",
            "--lock",
            &lock(&dir),
            "--quiet",
        ])
        .0
    );
    let v2 = dir.join("verdict2.json");
    let b2 = dir.join("bound2.json");
    let (ok, _, err) = ecbench(&[
        "challenge",
        "verdict",
        "--challenge",
        &s(&challenge),
        "--dir",
        &s(&s2),
        "--epoch",
        "2",
        "--replay-all",
        "--root",
        &s(&dir),
        "--out",
        &s(&v2),
        "--bound-out",
        &s(&b2),
        "--exit-code",
    ]);
    assert!(ok, "verdict 2: {err}");
    let v = read_json(&v2);
    assert_eq!(v["outcome"], "matches", "{}", v["statement"]);
    assert_eq!(v["axes"][0]["ratio_candidate_over_incumbent"], 1.0);
    assert!(!b2.exists(), "a match writes no bound");

    // The wrong epoch is somebody else's spec: inadmissible, exit 1, and
    // the paired ratio is still reported for information.
    let v3 = dir.join("verdict3.json");
    let (ok, _, err) = ecbench(&[
        "challenge",
        "verdict",
        "--challenge",
        &s(&challenge),
        "--dir",
        &s(&s2),
        "--epoch",
        "3",
        "--root",
        &s(&dir),
        "--out",
        &s(&v3),
        "--exit-code",
    ]);
    assert!(!ok, "{err}");
    let v = read_json(&v3);
    assert_eq!(v["outcome"], "inadmissible");
    assert_eq!(v["session"]["spec_matches"], false);
    assert!(v["reasons"]
        .as_array()
        .unwrap()
        .iter()
        .any(|r| r.as_str().unwrap().contains("not the challenge's spec")));
    assert!(v["reasons"]
        .as_array()
        .unwrap()
        .iter()
        .any(|r| r.as_str().unwrap().contains("replay-all")));
    let _ = std::fs::remove_dir_all(&dir);
}

const FIELD_SPEC: &str = r#"{
  "schema": "ecbench.spec/v1",
  "label": "field operations integration test: four prime sizes",
  "workloads": {
    "curves": [{"kind": "prime_search", "bits": 12, "seed": 59297},
               {"kind": "prime_search", "bits": 14, "seed": 59297},
               {"kind": "prime_search", "bits": 16, "seed": 59297},
               {"kind": "prime_search", "bits": 18, "seed": 59297}],
    "targets_per_curve": 2,
    "target_seed": 13
  },
  "arms": [
    {"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
    {"name": "bsgs-neg", "role": "candidate", "method": {"id": "bsgs.negation"}}
  ],
  "measurement": {"rounds": 1, "warmup": 0, "seed": 9, "isolation_required": "L0", "timeout_seconds": 60}
}"#;

const KOBLITZ_FIELD_SPEC: &str = r#"{
  "schema": "ecbench.spec/v1",
  "label": "field operations integration test: a Koblitz curve counts none",
  "workloads": {
    "curves": [{"kind": "koblitz", "a": 0, "n": 13}],
    "targets_per_curve": 1,
    "target_seed": 13
  },
  "arms": [
    {"name": "bsgs-neg", "role": "reference", "method": {"id": "bsgs.negation"}}
  ],
  "measurement": {"rounds": 1, "warmup": 0, "seed": 9, "isolation_required": "L0", "timeout_seconds": 60}
}"#;

/// The primitive level, end to end: prime-field curves count the
/// multiplications, squarings and inversions behind their group
/// operations; the block rides on every record outside `counters` and
/// `phases`; a bound carries the three field axes with intervals; the
/// frontier page grows their columns; a challenge may name one, and the
/// verdict then decides on it and reports the other two.  A Koblitz curve
/// counts none and its record and bound say nothing, not zero.
#[test]
fn field_operations_are_counted_bounded_and_judged() {
    let dir = scratch("field-ops");
    let s = |p: &Path| p.to_str().unwrap().to_string();
    let spec = dir.join("spec.json");
    std::fs::write(&spec, FIELD_SPEC).unwrap();
    let session = dir.join("s0");
    let (ok, _, err) = ecbench(&[
        "run",
        "--spec",
        &s(&spec),
        "--out",
        &s(&session),
        "--cpus",
        "none",
        "--lock",
        &lock(&dir),
        "--quiet",
    ]);
    assert!(ok, "run failed: {err}");

    // Every record carries the block, and nothing the replay of a
    // committed record compares.  Every addition or doubling that did
    // field work did exactly one inversion, at least two multiplications
    // and at least one squaring (`PrimeCurve::FIELD_OPS_PER_ADD`,
    // `FIELD_OPS_PER_DOUBLE`); the special cases did none.  The group
    // operations are `total_gae` (additions plus doublings): the tuned
    // rho's `search` phase carries its gae without an operation split.
    let lines = std::fs::read_to_string(session.join("records.jsonl")).unwrap();
    let mut records = 0;
    for line in lines.lines().filter(|l| !l.trim().is_empty()) {
        let r: Value = serde_json::from_str(line).unwrap();
        assert_eq!(r["outcome"]["status"], "verified", "{}", r["run_id"]);
        let f = &r["field_ops"];
        let (muls, sqrs, invs) = (
            f["muls"].as_u64().unwrap(),
            f["sqrs"].as_u64().unwrap(),
            f["invs"].as_u64().unwrap(),
        );
        let group_ops = r["cost"]["total_gae"].as_f64().unwrap();
        assert!(
            invs > 0 && invs as f64 <= group_ops && muls >= 2 * invs && sqrs >= invs,
            "{}: {f} over {group_ops} group operations",
            r["run_id"]
        );
        assert!(r["counters"]
            .as_object()
            .unwrap()
            .keys()
            .all(|k| !k.contains("field")));
        assert!(r["phases"]
            .as_array()
            .unwrap()
            .iter()
            .all(|p| p.get("field_ops").is_none()));
        records += 1;
    }
    assert_eq!(records, 16);
    let (ok, _, err) = ecbench(&[
        "verify",
        "--dir",
        &s(&session),
        "--replay-all",
        "--exit-code",
    ]);
    assert!(ok, "replay: {err}");

    // Each arm's bound carries the three field axes, known, with the
    // two-stage interval the other axes have, and re-derives.
    let records_dir = dir.join("records");
    std::fs::create_dir_all(&records_dir).unwrap();
    let mut ids = std::collections::BTreeMap::new();
    for arm in ["rho-neg", "bsgs-neg"] {
        let out = records_dir.join(format!("{arm}.json"));
        let (ok, _, err) = ecbench(&[
            "bound",
            "fit",
            "--dir",
            &s(&session),
            "--arm",
            arm,
            "--root",
            &s(&dir),
            "--out",
            &s(&out),
        ]);
        assert!(ok, "bound fit {arm}: {err}");
        let b = read_json(&out);
        assert_eq!(b["fit"]["sizes"], 4);
        for axis in ["field_muls", "field_sqrs", "field_invs"] {
            let d = &b["dimensions"][axis];
            assert_eq!(d["known"], true, "{arm} {axis}: {d}");
            assert!(d["value"].as_f64().unwrap() > 0.0, "{arm} {axis}: {d}");
            assert_eq!(d["ci95"].as_array().unwrap().len(), 2, "{arm} {axis}: {d}");
            assert_eq!(d["lower_is_better"], true);
        }
        // Per `√r`: at most one inversion per group operation, so the
        // inversions axis is bounded by `S`, and two multiplications each.
        let dim = |axis: &str| b["dimensions"][axis]["value"].as_f64().unwrap();
        assert!(dim("field_invs") <= b["constant"]["s"]["value"].as_f64().unwrap());
        assert!(dim("field_muls") >= 2.0 * dim("field_invs"));
        assert!(dim("field_sqrs") >= dim("field_invs"));
        ids.insert(arm, b["bound_id"].as_str().unwrap().to_string());
    }
    let (ok, _, err) = ecbench(&[
        "bound",
        "check",
        "--root",
        &s(&dir),
        "--record",
        &s(&records_dir.join("rho-neg.json")),
        &s(&records_dir.join("bsgs-neg.json")),
    ]);
    assert!(ok, "bound check: {err}");

    // The frontier carries the axes and the page their columns; the axes
    // are in the vocabulary dominance may read.
    let fjson = dir.join("frontier.json");
    let fmd = dir.join("FRONTIER.md");
    let (ok, _, err) = ecbench(&[
        "frontier",
        "build",
        "--bounds",
        &s(&records_dir),
        "--out",
        &s(&fjson),
        "--markdown",
        &s(&fmd),
    ]);
    assert!(ok, "frontier build: {err}");
    let f = read_json(&fjson);
    let entries = f["domains"][0]["entries"].as_array().unwrap();
    assert_eq!(entries.len(), 2);
    for e in entries {
        for axis in ["field_muls", "field_sqrs", "field_invs"] {
            assert_eq!(e["axes"][axis]["known"], true, "{}", e["bound_id"]);
        }
    }
    let page = std::fs::read_to_string(&fmd).unwrap();
    assert!(page.contains("| field muls (/√r) | field sqrs (/√r) | field invs (/√r) |"));
    assert!(page.contains("the primitive level"));
    let (ok, _, err) = ecbench(&[
        "frontier",
        "build",
        "--bounds",
        &s(&records_dir),
        "--out",
        &s(&fjson),
        "--markdown",
        &s(&fmd),
        "--check",
    ]);
    assert!(ok, "frontier check: {err}");
    let (ok, _, err) = ecbench(&[
        "frontier",
        "build",
        "--bounds",
        &s(&records_dir),
        "--axes",
        "ops,field_sqrs",
        "--json",
    ]);
    assert!(ok, "frontier on a field axis: {err}");

    // A challenge naming `field_sqrs`: the verdict decides on it and
    // reports the other two field axes.
    let draft = dir.join("draft.json");
    std::fs::write(
        &draft,
        format!(
            r#"{{
  "label": "test: fewer squarings than rho.negation",
  "domain": {{"problem": "ecdlp.single_target", "family": "prime", "target_kind": "planted",
             "unit": "ecbench.gae", "tier": "toy",
             "envelope": {{"targets": 1, "precomputation": "none", "threads": 1}}}},
  "incumbent": {{"bound_id": "{}", "method": {{"id": "rho.negation"}}}},
  "workloads": {{"curves": [{{"kind": "prime_search", "bits": 12, "seed": 59297}},
                           {{"kind": "prime_search", "bits": 14, "seed": 59297}},
                           {{"kind": "prime_search", "bits": 16, "seed": 59297}},
                           {{"kind": "prime_search", "bits": 18, "seed": 59297}}],
                "targets_per_curve": 2, "nonce": 3}},
  "measurement": {{"rounds": 1, "warmup": 0, "isolation_required": "L0", "timeout_seconds": 60}},
  "acceptance": {{"axes": ["ops", "memory", "field_sqrs"], "min_runs_per_size": 2}}
}}"#,
            ids["rho-neg"]
        ),
    )
    .unwrap();
    let challenge = dir.join("challenge.json");
    let (ok, _, err) = ecbench(&[
        "challenge",
        "seal",
        "--draft",
        &s(&draft),
        "--out",
        &s(&challenge),
    ]);
    assert!(ok, "challenge seal: {err}");
    let (ok, _, err) = ecbench(&["challenge", "check", "--file", &s(&challenge)]);
    assert!(ok, "challenge check: {err}");
    let c = read_json(&challenge);
    assert_eq!(
        c["acceptance"]["axes"],
        serde_json::json!(["ops", "memory", "field_sqrs"])
    );
    let spec1 = dir.join("spec1.json");
    let (ok, _, err) = ecbench(&[
        "challenge",
        "spec",
        "--challenge",
        &s(&challenge),
        "--candidate",
        r#"{"id":"bsgs.negation"}"#,
        "--epoch",
        "1",
        "--out",
        &s(&spec1),
    ]);
    assert!(ok, "challenge spec: {err}");
    let s1 = dir.join("s1");
    let (ok, _, err) = ecbench(&[
        "run",
        "--spec",
        &s(&spec1),
        "--out",
        &s(&s1),
        "--cpus",
        "none",
        "--lock",
        &lock(&dir),
        "--quiet",
    ]);
    assert!(ok, "run s1: {err}");
    let v1 = dir.join("verdict1.json");
    let b1 = dir.join("bound1.json");
    let (_, _, err) = ecbench(&[
        "challenge",
        "verdict",
        "--challenge",
        &s(&challenge),
        "--dir",
        &s(&s1),
        "--epoch",
        "1",
        "--replay-all",
        "--bounds",
        &s(&records_dir),
        "--root",
        &s(&dir),
        "--out",
        &s(&v1),
        "--bound-out",
        &s(&b1),
    ]);
    assert!(v1.exists(), "verdict 1 was not written: {err}");
    let v = read_json(&v1);
    assert_eq!(v["schema"], "ecbench.verdict/v1");
    assert_ne!(v["outcome"], "inadmissible", "{}", v["statement"]);
    assert_eq!(v["audit"]["ok"], true);
    assert_eq!(v["audit"]["replays"], v["audit"]["replays_reproduced"]);
    assert_eq!(
        v["acceptance"]["axes"],
        serde_json::json!(["ops", "memory", "field_sqrs"])
    );
    let axes = v["axes"].as_array().unwrap();
    let axis = |name: &str| {
        axes.iter()
            .find(|a| a["axis"] == name)
            .unwrap_or_else(|| panic!("no axis {name} in {axes:?}"))
    };
    let ops_pairs = axis("ops")["pairs"].as_u64().unwrap();
    assert_eq!(ops_pairs, 8);
    for (name, decides) in [
        ("field_sqrs", true),
        ("field_muls", false),
        ("field_invs", false),
    ] {
        let a = axis(name);
        assert_eq!(a["known"], true, "{a}");
        assert_eq!(a["decides"], decides, "{a}");
        assert_eq!(a["pairs"], ops_pairs, "{a}");
        assert!(a["ratio_candidate_over_incumbent"].as_f64().unwrap() > 0.0);
        assert!(
            ["better", "worse", "indistinguishable"].contains(&a["verdict"].as_str().unwrap()),
            "{a}"
        );
    }
    // The axes named in the statement include the field ones, deciding or
    // reported; the level, if any, is one of the three.
    let statement = v["statement"].as_str().unwrap();
    assert!(statement.contains("field_sqrs"), "{statement}");
    assert!(
        statement.contains("field_muls") && statement.contains("reported, not deciding"),
        "{statement}"
    );
    match v["level_moved"].as_str() {
        None => assert_ne!(v["outcome"], "advances"),
        Some(l) => {
            assert_eq!(v["outcome"], "advances");
            assert!(["exponent", "constant", "primitive"].contains(&l), "{l}");
        }
    }
    if b1.exists() {
        let nb = read_json(&b1);
        assert_eq!(nb["dimensions"]["field_sqrs"]["known"], true);
        assert_eq!(nb["verdict_id"], v["verdict_id"]);
    }

    // A Koblitz curve does not count: no block on the record, no axis on
    // the bound.  Unknown is not zero.
    let kspec = dir.join("kspec.json");
    std::fs::write(&kspec, KOBLITZ_FIELD_SPEC).unwrap();
    let ksession = dir.join("k0");
    let (ok, _, err) = ecbench(&[
        "run",
        "--spec",
        &s(&kspec),
        "--out",
        &s(&ksession),
        "--cpus",
        "none",
        "--lock",
        &lock(&dir),
        "--quiet",
    ]);
    assert!(ok, "koblitz run failed: {err}");
    let lines = std::fs::read_to_string(ksession.join("records.jsonl")).unwrap();
    let krecords: Vec<Value> = lines
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(krecords.len(), 1);
    assert_eq!(krecords[0]["outcome"]["status"], "verified");
    assert!(krecords[0].get("field_ops").is_none(), "absent, never zero");
    let kout = dir.join("koblitz-bound.json");
    let (ok, _, err) = ecbench(&[
        "bound",
        "fit",
        "--dir",
        &s(&ksession),
        "--arm",
        "bsgs-neg",
        "--root",
        &s(&dir),
        "--out",
        &s(&kout),
    ]);
    assert!(ok, "koblitz bound fit: {err}");
    let kb = read_json(&kout);
    let dims: Vec<&String> = kb["dimensions"].as_object().unwrap().keys().collect();
    assert_eq!(dims, vec!["memory", "ops", "uncharged"], "{dims:?}");

    // A Koblitz challenge that names no field axis, judged on a session
    // that counted none: the verdict carries exactly the axes it always
    // carried, so a committed verdict re-derives under this binary.
    let kdraft = dir.join("kdraft.json");
    std::fs::write(
        &kdraft,
        r#"{
  "label": "test: koblitz, no field axes",
  "domain": {"problem": "ecdlp.single_target", "family": "koblitz", "target_kind": "planted",
             "unit": "ecbench.gae", "tier": "toy",
             "envelope": {"targets": 1, "precomputation": "none", "threads": 1}},
  "incumbent": {"bound_id": null, "method": {"id": "rho.negation"}},
  "workloads": {"curves": [{"kind": "koblitz", "a": 1, "n": 17},
                           {"kind": "koblitz", "a": 0, "n": 19},
                           {"kind": "koblitz", "a": 0, "n": 23},
                           {"kind": "koblitz", "a": 1, "n": 29}],
                "targets_per_curve": 2, "nonce": 5},
  "measurement": {"rounds": 1, "warmup": 0, "isolation_required": "L0", "timeout_seconds": 60},
  "acceptance": {"min_runs_per_size": 2}
}"#,
    )
    .unwrap();
    let kchallenge = dir.join("kchallenge.json");
    let (ok, _, err) = ecbench(&[
        "challenge",
        "seal",
        "--draft",
        &s(&kdraft),
        "--out",
        &s(&kchallenge),
    ]);
    assert!(ok, "koblitz challenge seal: {err}");
    let kspec1 = dir.join("kspec1.json");
    let (ok, _, err) = ecbench(&[
        "challenge",
        "spec",
        "--challenge",
        &s(&kchallenge),
        "--candidate",
        r#"{"id":"bsgs.negation"}"#,
        "--epoch",
        "1",
        "--out",
        &s(&kspec1),
    ]);
    assert!(ok, "koblitz challenge spec: {err}");
    let k1 = dir.join("k1");
    let (ok, _, err) = ecbench(&[
        "run",
        "--spec",
        &s(&kspec1),
        "--out",
        &s(&k1),
        "--cpus",
        "none",
        "--lock",
        &lock(&dir),
        "--quiet",
    ]);
    assert!(ok, "koblitz run k1: {err}");
    let kv1 = dir.join("kverdict1.json");
    let (_, _, err) = ecbench(&[
        "challenge",
        "verdict",
        "--challenge",
        &s(&kchallenge),
        "--dir",
        &s(&k1),
        "--epoch",
        "1",
        "--replay-all",
        "--root",
        &s(&dir),
        "--out",
        &s(&kv1),
    ]);
    assert!(kv1.exists(), "koblitz verdict was not written: {err}");
    let kv = read_json(&kv1);
    assert_ne!(kv["outcome"], "inadmissible", "{}", kv["statement"]);
    let names: Vec<&str> = kv["axes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|a| a["axis"].as_str().unwrap())
        .collect();
    assert_eq!(names, vec!["ops", "memory", "uncharged"], "{names:?}");
    assert!(!kv["statement"].as_str().unwrap().contains("field_"));
    let _ = std::fs::remove_dir_all(&dir);
}
