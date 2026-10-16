//! Fixed public-map controls: tracing changes evidence, never mathematical verdicts.
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicUsize, Ordering},
};
fn save(path: &Path, value: &Value) {
    fs::write(path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
}
fn read(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}
fn path(path: &Path) -> &str {
    path.to_str().unwrap()
}
fn call(root: &Path, args: &[&str], trace: Option<&Path>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_isogeny"));
    command
        .args(args)
        .current_dir(std::env::temp_dir())
        .env_remove("ISOGENY_REPLAY_TRACE");
    if let Some(trace) = trace {
        command.env("ISOGENY_REPLAY_TRACE", trace);
    }
    let output = command.output().unwrap();
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let receipt = root.join(format!(
        "command-{:03}.json",
        NEXT.fetch_add(1, Ordering::SeqCst)
    ));
    assert!(!receipt.exists());
    save(
        &receipt,
        &json!({"args":args, "trace":trace, "exit_code":output.status.code(), "stdout":String::from_utf8_lossy(&output.stdout), "stderr":String::from_utf8_lossy(&output.stderr)}),
    );
    output
}
fn result(output: Output, code: i32) -> Value {
    assert_eq!(
        output.status.code(),
        Some(code),
        "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}
fn clone_artifact(source: &Path, target: &Path) {
    fs::create_dir(target).unwrap();
    fs::create_dir(target.join("p192")).unwrap();
    for file in ["search.json", "p192/ell-13.json", "p192/ell-13.stderr.txt"] {
        fs::copy(source.join(file), target.join(file)).unwrap();
    }
}
fn alter_map(root: &Path, change: impl FnOnce(&mut Value)) {
    let file = root.join("p192/ell-13.json");
    let mut document = read(&file);
    change(&mut document);
    save(&file, &document);
    let mut summary = read(&root.join("search.json"));
    summary["attempts"][0]["receipt"]["stdout_sha256"] =
        json!(hex::encode(Sha256::digest(fs::read(file).unwrap())));
    save(&root.join("search.json"), &summary);
}
fn trace_records(path: &Path) -> Vec<Value> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect()
}
fn structure(records: &[Value]) -> Vec<Value> {
    records
        .iter()
        .map(|record| {
            json!([
                record["stage"],
                record["event"],
                record["stage_id"],
                record["parent_id"],
                record["status"]
            ])
        })
        .collect()
}
fn check_complete_trace(records: &[Value]) {
    let mut stack = vec![];
    let mut elapsed = 0;
    for record in records {
        assert_eq!(record["schema"], "isogeny-replay-stage/v1");
        let now = record["elapsed_ms_operational_only"].as_u64().unwrap();
        assert!(now >= elapsed);
        elapsed = now;
        match record["event"].as_str().unwrap() {
            "begin" => {
                assert_eq!(
                    record["parent_id"],
                    stack.last().cloned().unwrap_or(Value::Null)
                );
                stack.push(record["stage_id"].clone());
            }
            "end" => assert_eq!(stack.pop().unwrap(), record["stage_id"]),
            "terminal" => {
                assert!(stack.is_empty());
                assert_eq!(record["status"], "COMPLETED");
            }
            other => panic!("unexpected successful trace event {other}"),
        }
    }
    assert!(stack.is_empty());
    assert_eq!(records.last().unwrap()["event"], "terminal");
    for name in [
        "setup/loading/parsing",
        "squarefree_kernel_gcd",
        "selected_division_polynomial_torsion",
        "subgroup_inverse",
        "subgroup_modular_composition",
        "denominator_kernel_square",
        "numerator_denominator_gcd",
        "exact_rational_map_substitution",
        "target_generator_subgroup_order",
        "point_transport_19",
        "total",
    ] {
        assert!(
            records
                .iter()
                .any(|record| record["stage"] == name && record["event"] == "begin"),
            "missing {name}"
        );
    }
}
#[test]
fn tracing_preserves_fixed_map_and_adversarial_verdicts_and_refuses_bad_evidence_paths() {
    let retained = std::env::var_os("ISOGENY_TRACE_EVIDENCE");
    let root = retained.as_ref().map(PathBuf::from).unwrap_or_else(|| {
        std::env::temp_dir().join(format!(
            "isogeny-trace-controls-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    });
    fs::create_dir(&root).unwrap();
    let root = fs::canonicalize(root).unwrap();
    let generated = root.join("public-p192-13");
    let constructed = result(
        call(
            &root,
            &[
                "challenge",
                "run",
                "--curve",
                "p192",
                "--ell",
                "13",
                "--timeout",
                "60",
                "--out",
                path(&generated),
            ],
            None,
        ),
        0,
    );
    assert_eq!(constructed["status"], "ACCEPTED");
    let original = generated.join("construction");
    let off = root.join("valid-off");
    let on = root.join("valid-on");
    let trace = root.join("valid.jsonl");
    let no_trace = result(
        call(
            &root,
            &[
                "challenge",
                "judge",
                "--input",
                path(&original),
                "--out",
                path(&off),
            ],
            None,
        ),
        0,
    );
    let traced = result(
        call(
            &root,
            &[
                "challenge",
                "judge",
                "--input",
                path(&original),
                "--out",
                path(&on),
            ],
            Some(&trace),
        ),
        0,
    );
    assert_eq!(no_trace["score"], traced["score"]);
    assert_eq!(no_trace["source_icv1"], traced["source_icv1"]);
    assert_eq!(
        fs::read(off.join("verification.json")).unwrap(),
        fs::read(on.join("verification.json")).unwrap()
    );
    let records = trace_records(&trace);
    check_complete_trace(&records);
    let second = root.join("valid-second");
    let second_trace = root.join("valid-second.jsonl");
    result(
        call(
            &root,
            &[
                "challenge",
                "judge",
                "--input",
                path(&original),
                "--out",
                path(&second),
            ],
            Some(&second_trace),
        ),
        0,
    );
    assert_eq!(
        structure(&records),
        structure(&trace_records(&second_trace))
    );
    for kind in [
        "composite",
        "inflated-degree",
        "source",
        "kernel",
        "map",
        "compact-only",
        "traversal",
        "missing-hash",
        "duplicate-kernel",
    ] {
        let input = root.join(format!("bad-{kind}"));
        clone_artifact(&original, &input);
        match kind {
            "composite" => {
                let mut summary = read(&input.join("search.json"));
                summary["attempts"][0]["ell"] = json!(15);
                save(&input.join("search.json"), &summary);
            }
            "inflated-degree" => {
                alter_map(&input, |document| document["result"]["ell"] = json!(999983))
            }
            "source" => alter_map(&input, |document| document["curve"]["a"] = json!("1")),
            "kernel" => alter_map(&input, |document| {
                let v = &mut document["result"]["isogenies"][0]["kernel"][0];
                *v = if *v == "1" { json!("2") } else { json!("1") };
            }),
            "map" => alter_map(&input, |document| {
                let v = &mut document["result"]["isogenies"][0]["x_map"]["num"][0];
                *v = if *v == "1" { json!("2") } else { json!("1") };
            }),
            "compact-only" => alter_map(&input, |document| {
                document["result"]["isogenies"][0]["x_map"] = json!({"ideal":"unexpanded claim"})
            }),
            "duplicate-kernel" => alter_map(&input, |document| {
                document["result"]["isogenies"][1] = document["result"]["isogenies"][0].clone()
            }),
            "traversal" => {
                let mut summary = read(&input.join("search.json"));
                summary["attempts"][0]["curve"] = json!("../../p192");
                save(&input.join("search.json"), &summary);
            }
            "missing-hash" => {
                let mut summary = read(&input.join("search.json"));
                summary["attempts"][0]["receipt"]
                    .as_object_mut()
                    .unwrap()
                    .remove("stdout_sha256");
                save(&input.join("search.json"), &summary);
            }
            _ => unreachable!(),
        }
        let off = root.join(format!("{kind}-off"));
        let on = root.join(format!("{kind}-on"));
        let trace = root.join(format!("{kind}.jsonl"));
        let a = result(
            call(
                &root,
                &[
                    "challenge",
                    "judge",
                    "--input",
                    path(&input),
                    "--out",
                    path(&off),
                ],
                None,
            ),
            1,
        );
        let b = result(
            call(
                &root,
                &[
                    "challenge",
                    "judge",
                    "--input",
                    path(&input),
                    "--out",
                    path(&on),
                ],
                Some(&trace),
            ),
            1,
        );
        assert_eq!(a["status"], b["status"]);
        assert!(a["score"].is_null() && b["score"].is_null());
        if matches!(kind, "kernel" | "map") {
            let records = trace_records(&trace);
            assert_eq!(records.last().unwrap()["event"], "terminal");
            assert_eq!(records.last().unwrap()["status"], "ERROR");
            assert!(records.iter().any(|record| record["event"] == "error"));
        }
    }
    let malformed = root.join("malformed");
    fs::create_dir(&malformed).unwrap();
    fs::write(malformed.join("search.json"), b"{").unwrap();
    let trace = root.join("malformed.jsonl");
    let a = call(&root, &["verify", "--input", path(&malformed)], None);
    let b = call(
        &root,
        &["verify", "--input", path(&malformed)],
        Some(&trace),
    );
    assert_eq!(a.status.code(), b.status.code());
    assert!(!a.status.success());
    assert_eq!(a.stdout, b.stdout);
    let failed = trace_records(&trace);
    assert_eq!(failed.last().unwrap()["status"], "ERROR");
    assert!(!failed
        .iter()
        .any(|record| record["stage"] == "setup/loading/parsing" && record["event"] == "end"));
    let existing = root.join("existing.jsonl");
    fs::write(&existing, b"retained trace sentinel\n").unwrap();
    for trace in [
        &existing,
        &root.join("missing-parent/trace.jsonl"),
        Path::new("relative-trace-forbidden.jsonl"),
        &root,
    ] {
        let output = call(&root, &["verify", "--input", path(&original)], Some(trace));
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("requested replay trace"));
    }
    assert_eq!(fs::read(existing).unwrap(), b"retained trace sentinel\n");
    save(
        &root.join("control-summary.json"),
        &json!({"schema":"isogeny-replay-trace-controls/v1","public_curve":"P-192","prime_degree":"13","large_target_executed":false,"mathematical_verification_bytes_identical":true,"adversarial_pairs":9,"malformed_pair":1,"trace_structures_identical":true,"timing_comparison":false}),
    );
    if retained.is_none() {
        fs::remove_dir_all(root).unwrap();
    }
}
