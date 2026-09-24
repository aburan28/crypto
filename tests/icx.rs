//! Integration tests for the `icx` binary: they run the actual executable
//! across the whole standardized-curve catalog, which is the completion
//! criterion for the engine ("you can run this against every curve").

use std::process::Command;

fn icx() -> Command {
    Command::new(env!("CARGO_BIN_EXE_icx"))
}

/// Parse the top-level JSON object printed by an `icx --json` invocation.
fn run_json(args: &[&str]) -> serde_json::Value {
    let out = icx()
        .args(args)
        .arg("--json")
        .output()
        .expect("run icx");
    assert!(
        out.status.success(),
        "icx {:?} exited with failure: {}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_slice(&out.stdout)
        .unwrap_or_else(|e| panic!("icx {args:?} did not print valid JSON: {e}"))
}

#[test]
fn list_reports_every_family() {
    let v = run_json(&["list"]);
    assert_eq!(v["operation"], "list");
    let curves = v["curves"].as_array().expect("curves array");
    assert!(curves.len() >= 35, "catalog unexpectedly small: {}", curves.len());
    // Every named field family is represented.
    let families: std::collections::HashSet<&str> =
        curves.iter().filter_map(|c| c["family"].as_str()).collect();
    for fam in ["prime", "binary", "koblitz"] {
        assert!(families.contains(fam), "family {fam} missing from listing");
    }
    // A provenance block records the dispatched field backend.
    assert!(v["software"]["field_backend"].is_string());
}

#[test]
fn inspect_verifies_every_catalog_curve() {
    // The load-bearing gate: run `icx inspect` on EVERY curve the catalog
    // lists, and require each to verify. This is what "run icx against every
    // standardized curve" means for the inspect front door.
    let list = run_json(&["list"]);
    let curves = list["curves"].as_array().unwrap().clone();
    assert!(!curves.is_empty());
    let mut checked = 0;
    for c in &curves {
        let name = c["name"].as_str().unwrap();
        let rep = run_json(&["inspect", name]);
        assert_eq!(rep["status"], "checks_passed", "curve {name} failed inspect");
        assert_eq!(rep["verified"], true, "curve {name} not verified");
        checked += 1;
    }
    assert!(checked >= 35, "only {checked} curves checked");
}

#[test]
fn estimate_labels_prime_as_not_ic_relevant() {
    let v = run_json(&["estimate", "p256"]);
    assert_eq!(v["index_calculus"]["relevant"], false);
    // ~128-bit generic security.
    let sec = v["generic_attack"]["security_bits"].as_f64().unwrap();
    assert!((120.0..132.0).contains(&sec), "p256 security {sec}");
}

#[test]
fn estimate_labels_koblitz_as_ic_relevant() {
    let v = run_json(&["estimate", "sect163k1"]);
    assert_eq!(v["index_calculus"]["relevant"], true);
    assert_eq!(v["family"], "koblitz");
}

#[test]
fn aliases_resolve() {
    // NIST/SEC aliases must map to the canonical curve.
    for (alias, canonical) in [("nistp256", "p256"), ("k-163", "sect163k1")] {
        let v = run_json(&["inspect", alias]);
        assert_eq!(v["curve"]["name"], canonical, "alias {alias}");
    }
}

#[test]
fn unknown_curve_errors_cleanly() {
    let out = icx().args(["inspect", "no-such-curve"]).output().unwrap();
    assert!(!out.status.success());
}
