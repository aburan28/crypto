//! Integration tests for the `ecdlp-nist` binary: run the executable across
//! the fifteen NIST curves, the toy curves and a constructed anomalous
//! curve, and check the JSON it prints.

use std::process::Command;

use serde_json::Value;

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ecdlp-nist"))
}

fn run_json(args: &[&str]) -> Value {
    let out = bin()
        .args(args)
        .arg("--json")
        .output()
        .expect("run ecdlp-nist");
    assert!(
        out.status.success(),
        "ecdlp-nist {:?} failed: {}\n{}",
        args,
        String::from_utf8_lossy(&out.stderr),
        String::from_utf8_lossy(&out.stdout)
    );
    serde_json::from_slice(&out.stdout)
        .unwrap_or_else(|e| panic!("ecdlp-nist {args:?} did not print valid JSON: {e}"))
}

#[test]
fn list_covers_the_fifteen_nist_curves_and_finds_no_structural_weakness() {
    let v = run_json(&["list"]);
    assert_eq!(v["operation"], "list");
    assert_eq!(v["count"], 15);
    let curves = v["curves"].as_array().unwrap();
    let mut names: Vec<&str> = curves
        .iter()
        .map(|c| c["curve"]["name"].as_str().unwrap())
        .collect();
    names.sort_unstable();
    let mut expected = vec![
        "P-192", "P-224", "P-256", "P-384", "P-521", "K-163", "K-233", "K-283", "K-409", "K-571",
        "B-163", "B-233", "B-283", "B-409", "B-571",
    ];
    expected.sort_unstable();
    assert_eq!(names, expected);
    for c in curves {
        let a = &c["audit"];
        assert_eq!(a["anomalous"], false, "{}", c["curve"]["name"]);
        assert_eq!(a["order_is_prime"], true, "{}", c["curve"]["name"]);
        assert!(a["embedding_degree"].is_null(), "{}", c["curve"]["name"]);
        assert_eq!(a["recommended"], "rho");
    }
}

#[test]
fn audit_reports_frobenius_fold_on_koblitz() {
    let v = run_json(&["audit", "K-283"]);
    let a = &v["audit"];
    assert_eq!(a["family"], "koblitz");
    assert_eq!(a["rho_class_size"], 2 * 283);
    assert_eq!(a["koblitz_order_consistent"], true);
    assert_eq!(a["extension_degree_prime"], true);
    assert!(a["frobenius_eigenvalue"].is_string());
    let folded = a["rho_log2_expected"].as_f64().unwrap();
    let plain = a["rho_log2_expected_unfolded"].as_f64().unwrap();
    assert!((plain - folded - (2.0 * 283f64).log2() / 2.0).abs() < 1e-9);
}

#[test]
fn solve_recovers_planted_scalars_on_prime_koblitz_and_binary_curves() {
    for (curve, method) in [
        ("P-521", "bsgs"),
        ("K-571", "kangaroo"),
        ("B-163", "auto"),
        ("P-224", "kangaroo"),
    ] {
        let v = run_json(&[
            "solve",
            curve,
            "--secret",
            "0x1234567",
            "--interval-bits",
            "26",
            "--method",
            method,
            "--threads",
            "2",
        ]);
        assert_eq!(v["operation"], "solve");
        assert_eq!(v["planted_recovered"], true, "{curve} {method}: {v}");
        assert_eq!(v["report"]["verified"], true);
        assert_eq!(v["report"]["scalar"], "19088743");
    }
}

#[test]
fn whole_group_rho_on_toy_koblitz_and_refusal_on_real_curve() {
    let v = run_json(&[
        "solve",
        "toy-k23a1",
        "--random-bits",
        "21",
        "--threads",
        "2",
    ]);
    assert_eq!(v["report"]["method"], "rho");
    assert_eq!(v["planted_recovered"], true, "{v}");

    let out = bin()
        .args(["solve", "P-256", "--secret", "42", "--json"])
        .output()
        .unwrap();
    assert!(!out.status.success(), "infeasible rho must exit non-zero");
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert!(v["report"]["failure"]
        .as_str()
        .unwrap()
        .contains("feasibility"));
}

#[test]
fn smart_attack_breaks_a_constructed_anomalous_curve() {
    let v = run_json(&["smart", "--bits", "192", "--seed", "4"]);
    assert_eq!(v["operation"], "smart");
    assert_eq!(v["planted_recovered"], true, "{v}");
    assert_eq!(v["report"]["method"], "smart");
    assert_eq!(v["report"]["audit"]["anomalous"], true);
    let bits = v["curve"]["bits"].as_u64().unwrap();
    assert!((191..=194).contains(&bits), "{bits}");
}

#[test]
fn pohlig_hellman_runs_on_a_composite_order_toy() {
    let v = run_json(&["solve", "toy-k19a0-full", "--secret", "99999"]);
    assert_eq!(v["report"]["method"], "pohlig_hellman");
    assert_eq!(v["planted_recovered"], true, "{v}");
}
