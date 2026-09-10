use serde_json::{json, Value};
use std::{
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
fn path() -> PathBuf {
    std::env::temp_dir().join(format!(
        "ic-test-{}-{}.json",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}
fn command(args: &[&str]) -> (bool, Value) {
    let output = Command::new(env!("CARGO_BIN_EXE_ic"))
        .args(args)
        .arg("--json")
        .output()
        .unwrap();
    let value = serde_json::from_slice(&output.stdout)
        .unwrap_or_else(|_| panic!("not JSON: {:?} {:?}", output.stdout, output.stderr));
    (output.status.success(), value)
}
fn inspect(value: &Value) -> (bool, Value) {
    let file = path();
    std::fs::write(&file, serde_json::to_vec(value).unwrap()).unwrap();
    let result = command(&["inspect", "--file", file.to_str().unwrap()]);
    std::fs::remove_file(file).unwrap();
    result
}
fn generated() -> Value {
    let (ok, v) = command(&["generate", "--degree", "9", "--seed", "41"]);
    assert!(ok);
    v
}
#[test]
fn named_profiles_report_inspection_only_and_missing_representation() {
    for name in ["ecc2k-130", "ecc2k-95", "sect163k1", "secp256k1"] {
        let (ok, v) = command(&[name]);
        assert!(ok, "{name}: {v}");
        assert_eq!(v["operation"], "inspect");
        assert_eq!(v["capabilities"]["imported_target_solving"], false);
        assert_eq!(v["capabilities"]["factor_base_built"], false);
        assert!(v["software"]["binary_blake3"].as_str().unwrap().len() == 64);
        if name == "ecc2k-130" {
            assert_eq!(v["parameters"]["field"]["degree"], 131);
            assert_eq!(v["capabilities"]["point_checks_attempted"], false);
            assert!(v["checks"]
                .as_array()
                .unwrap()
                .iter()
                .any(|c| c["name"] == "field_representation" && c["status"] == "not_checked"));
        }
    }
}
#[test]
fn generated_documents_round_trip_and_reproduce() {
    let v = generated();
    assert_eq!(v, generated());
    let (ok, report) = inspect(&v);
    assert!(ok, "{report}");
    assert!(report["checks"]
        .as_array()
        .unwrap()
        .iter()
        .any(|c| c["name"] == "known_answer" && c["status"] == "pass"));
}
#[test]
fn malformed_binary_parameters_fail_without_false_success() {
    let original = generated();
    let mut reducible = original.clone();
    reducible["field"]["polynomial_terms"] = json!([9, 0]);
    let mut duplicate = original.clone();
    duplicate["field"]["polynomial_terms"] = json!([9, 1, 1, 0]);
    let mut noncanonical = original.clone();
    noncanonical["generator"]["x"] = json!("0x200");
    let mut wrong_order = original.clone();
    wrong_order["subgroup_order"] = json!("6");
    let mut wrong_point = original.clone();
    wrong_point["point"] = json!({"x":"0","y":"0"});
    let mut wrong_answer = original;
    wrong_answer["fixture"]["known_log"] = json!("0");
    for value in [
        reducible,
        duplicate,
        noncanonical,
        wrong_order,
        wrong_point,
        wrong_answer,
    ] {
        let (ok, report) = inspect(&value);
        assert!(!ok, "accepted {report}");
        assert_eq!(report["status"], "invalid");
    }
}
#[test]
fn custom_prime_curve_checks_points_and_known_answer() {
    let v = json!({"schema_version":1,"name":"prime-toy","field":{"kind":"prime","modulus":"97"},
        "a":"2","b":"3","subgroup_order":"5","cofactor":"20",
        "generator":{"x":"3","y":"6"},"point":{"x":"3","y":"91"},"fixture":{"known_log":"4","seed":0}});
    let (ok, report) = inspect(&v);
    assert!(ok, "{report}");
    assert!(report["checks"]
        .as_array()
        .unwrap()
        .iter()
        .any(|c| c["name"] == "known_answer" && c["status"] == "pass"));
    let mut bad = v;
    bad["field"]["modulus"] = json!("91");
    assert!(!inspect(&bad).0);
}
#[test]
fn schema_unknown_fields_and_ambiguous_abstract_coordinates_are_rejected() {
    let mut v = generated();
    v["unrecognized"] = json!(true);
    assert!(!inspect(&v).0);
    let (_, normal) = command(&["ecc2k-130"]);
    let mut p = normal["parameters"].clone();
    p["generator"] = json!({"x":"1","y":"1"});
    let (ok, report) = inspect(&p);
    assert!(!ok);
    assert_eq!(report["status"], "invalid");
}
#[test]
fn supported_larger_fixture_completes_and_reports_resources() {
    // Control accounting: one column per Frobenius orbit, fixed surplus,
    // serial trials — the historical 91 + 4 relations.
    let (ok, v) = command(&[
        "run",
        "--degree",
        "11",
        "--curve-a",
        "1",
        "--known-log",
        "53",
        "--solver",
        "enumerate",
        "--control",
        "--batch",
        "1",
    ]);
    assert!(ok, "{v}");
    assert_eq!(v["result"]["verified"], true);
    assert_eq!(v["result"]["recovered"], "53");
    assert_eq!(v["counts"]["relations"], 95);
    assert_eq!(v["mode"]["control"], true);
    assert_eq!(v["factor_base"]["spec"]["kind"], "factor");
    // Fast accounting: signed and cofactor-projected columns, and the
    // incremental solve stops as soon as the scalar is pinned.
    let (ok, fast) = command(&[
        "run",
        "--degree",
        "11",
        "--curve-a",
        "1",
        "--known-log",
        "53",
        "--solver",
        "enumerate",
        "--batch",
        "1",
    ]);
    assert!(ok, "{fast}");
    assert_eq!(fast["result"]["verified"], true);
    assert_eq!(fast["result"]["recovered"], "53");
    assert!(fast["counts"]["columns"].as_u64().unwrap() < v["counts"]["columns"].as_u64().unwrap());
    assert!(fast["counts"]["relations"].as_u64().unwrap() < 95);
    assert_eq!(fast["counts"]["inconsistent_relations"], 0);
    assert_eq!(fast["counts"]["verification_failures"], 0);
    assert_eq!(
        v["stages"].as_array().unwrap().last().unwrap()["status"],
        "pass"
    );
    #[cfg(any(target_os = "macos", target_os = "linux"))]
    {
        assert!(v["resources"]["peak_rss_bytes"].as_u64().unwrap() > 0);
        assert!(v["resources"]["cpu_seconds"].as_f64().unwrap() > 0.0);
    }
}
#[test]
fn synthetic_targets_are_reproducible_and_imports_are_not_run_arguments() {
    let args = [
        "run",
        "--random-target",
        "--seed",
        "17",
        "--solver",
        "enumerate",
    ];
    let (ok, a) = command(&args);
    let (ok2, b) = command(&args);
    assert!(ok && ok2);
    assert_eq!(a["parameters"], b["parameters"]);
    let status = Command::new(env!("CARGO_BIN_EXE_ic"))
        .args(["run", "--file", "anything.json"])
        .output()
        .unwrap()
        .status;
    assert!(!status.success());
    let status = Command::new(env!("CARGO_BIN_EXE_ic"))
        .args(["run", "--degree", "131"])
        .output()
        .unwrap()
        .status;
    assert!(!status.success());
}
#[test]
fn incomplete_comparisons_never_select_a_winner() {
    let (ok, v) = command(&[
        "compare",
        "--samples",
        "1",
        "--holdout",
        "1",
        "--max-trials",
        "1",
    ]);
    assert!(!ok);
    assert_eq!(v["status"], "inconclusive");
    assert!(v["selected_factor_index"].is_null());
    for c in v["candidates"].as_array().unwrap() {
        assert_eq!(c["eligible"], false);
        assert!(c["median_process_seconds"].is_null());
    }
}
#[test]
fn successful_comparison_has_same_training_inputs_and_separate_holdout() {
    let (ok, v) = command(&[
        "compare",
        "--samples",
        "2",
        "--holdout",
        "1",
        "--max-trials",
        "100",
    ]);
    assert!(ok, "{v}");
    let candidates = v["candidates"].as_array().unwrap();
    assert_eq!(candidates.len(), 2);
    assert_eq!(
        candidates[0]["training"][0]["parameters"],
        candidates[1]["training"][0]["parameters"]
    );
    let index = v["selected_factor_index"].as_u64().unwrap() as usize;
    assert_eq!(candidates[index]["eligible"], true);
    assert_eq!(v["holdout"][0]["result"]["verified"], true);
    assert_ne!(
        v["holdout"][0]["arguments"]["seed"],
        candidates[0]["training"][0]["arguments"]["seed"]
    );
}
#[test]
fn report_files_are_created_once_and_generated_reports_are_importable() {
    let file = path();
    let (ok, v) = command(&["generate", "--out", file.to_str().unwrap()]);
    assert!(ok);
    let bytes = std::fs::read(&file).unwrap();
    assert_eq!(serde_json::from_slice::<Value>(&bytes).unwrap(), v);
    let run = Command::new(env!("CARGO_BIN_EXE_ic"))
        .args(["generate", "--out", file.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!run.status.success());
    assert_eq!(std::fs::read(&file).unwrap(), bytes);
    assert!(command(&["inspect", "--file", file.to_str().unwrap()]).0);
    std::fs::remove_file(file).unwrap();
}

#[test]
fn point_only_input_reports_the_checks_that_were_attempted() {
    let value = json!({"schema_version":1,"name":"point-only","field":{"kind":"prime","modulus":"97"},
        "a":"2","b":"3","subgroup_order":"5","cofactor":"20","point":{"x":"3","y":"6"}});
    let (ok, report) = inspect(&value);
    assert!(ok, "{report}");
    assert_eq!(report["capabilities"]["point_checks_attempted"], true);
    assert!(report["checks"]
        .as_array()
        .unwrap()
        .iter()
        .any(|c| c["name"] == "generator" && c["status"] == "not_checked"));
    assert!(report["checks"]
        .as_array()
        .unwrap()
        .iter()
        .any(|c| c["name"] == "point_subgroup" && c["status"] == "pass"));
}

#[test]
fn search_finds_a_validated_base_where_the_legacy_family_yields_nothing() {
    // K_1 / GF(2^15) with factor index 0 collects no relation at all.
    let (ok, legacy) = command(&[
        "run",
        "--degree",
        "15",
        "--curve-a",
        "1",
        "--solver",
        "pair-table",
        "--max-trials",
        "500",
    ]);
    assert!(!ok);
    assert_eq!(legacy["status"], "incomplete");
    assert_eq!(legacy["counts"]["relations"], 0);

    let spec = path();
    let (ok, v) = command(&[
        "search",
        "--degree",
        "15",
        "--curve-a",
        "1",
        "--family",
        "divisor",
        "--max-dimension",
        "8",
        "--no-saturate",
        "--validate-top",
        "1",
        "--holdout",
        "1",
        "--spec-out",
        spec.to_str().unwrap(),
    ]);
    assert!(ok, "{v}");
    assert_eq!(v["status"], "complete");
    assert_eq!(v["exhaustive_targets"], true);
    assert_eq!(v["targets"], 210);
    let candidates = v["candidates"].as_array().unwrap();
    assert!(!candidates.is_empty());
    // Ranked by expected trials, finite ones first.
    let scores: Vec<f64> = candidates
        .iter()
        .map(|c| c["expected_trials"].as_f64().unwrap_or(f64::INFINITY))
        .collect();
    assert!(scores.windows(2).all(|w| w[0] <= w[1]));
    assert!(candidates
        .iter()
        .all(|c| c["spec"]["kind"] == "divisor" || c["spec"]["kind"] == "pruned"));
    assert_eq!(v["validation"]["runs"][0]["eligible"], true);
    assert_eq!(v["selected"]["degree"], 15);
    assert_eq!(v["selected"]["curve_a"], 1);

    // The saved recipe solves a fresh known-answer fixture end to end.
    let (ok, run) = command(&[
        "run",
        "--degree",
        "15",
        "--curve-a",
        "1",
        "--solver",
        "pair-table",
        "--factor-base",
        spec.to_str().unwrap(),
        "--known-log",
        "97",
    ]);
    assert!(ok, "{run}");
    assert_eq!(run["result"]["verified"], true);
    assert_eq!(run["result"]["recovered"], "97");
    assert_eq!(run["factor_base"]["spec"], v["selected"]["spec"]);
    assert_eq!(run["counts"]["inconsistent_relations"], 0);

    // A recipe is bound to its curve.
    let (ok, wrong) = command(&[
        "run",
        "--degree",
        "9",
        "--curve-a",
        "0",
        "--factor-base",
        spec.to_str().unwrap(),
    ]);
    assert!(!ok);
    assert_eq!(wrong["operation"], "error");
    // The recipe file is never overwritten.
    let (ok, _) = command(&[
        "search",
        "--degree",
        "15",
        "--curve-a",
        "1",
        "--family",
        "factor",
        "--validate-top",
        "0",
        "--spec-out",
        spec.to_str().unwrap(),
    ]);
    assert!(!ok);
    std::fs::remove_file(spec).unwrap();
}

#[test]
fn every_solver_recovers_the_same_scalar_with_three_summands() {
    for solver in ["pair-table", "enumerate", "groebner", "sat"] {
        let (ok, v) = command(&[
            "run",
            "--degree",
            "9",
            "--curve-a",
            "0",
            "--known-log",
            "77",
            "--summands",
            "3",
            "--solver",
            solver,
        ]);
        assert!(ok, "{solver}: {v}");
        assert_eq!(v["result"]["verified"], true, "{solver}");
        assert_eq!(v["result"]["recovered"], "77", "{solver}");
        assert_eq!(v["counts"]["sat_invalid_models"], 0);
        assert_eq!(v["counts"]["inconsistent_relations"], 0);
    }
}
