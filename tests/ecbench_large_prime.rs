use std::collections::BTreeMap;
use std::path::Path;
use std::process::Command;

use crypto_lib::cryptanalysis::ecbench::methods::{resolve, solve as ecbench_solve, MethodSpec};
use crypto_lib::cryptanalysis::ecbench::workload::{CurveSpec, Workload};
use crypto_lib::cryptanalysis::ecbench_large_prime::{
    import_family, load_manifest, load_manifest_with_source, solve, SolveConfig,
};

fn fixture() -> crypto_lib::cryptanalysis::ecbench_large_prime::ImportedFamily {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ecbench_manifest_small.json");
    let manifest = load_manifest(&path).unwrap();
    import_family(&manifest.families[0]).unwrap().0
}

#[test]
fn corpus_source_identity_is_recorded_from_exact_bytes() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ecbench_manifest_small.json");
    let loaded = load_manifest_with_source(&path).unwrap();
    assert_eq!(loaded.bytes, 568);
    assert_eq!(
        loaded.sha256,
        "235f49a6d07a5f7afb18a37f574ebc7ee1ae0b69dc9d4a33cd0b8ff15556b5b7"
    );
    assert_eq!(loaded.manifest.families.len(), 1);
}

#[test]
fn ic_report_binds_the_exact_manifest_bytes() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ecbench_manifest_small.json");
    let output = Command::new(env!("CARGO_BIN_EXE_ic"))
        .args(["--json", "large-prime", "--manifest"])
        .arg(&path)
        .arg("--validate-only")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["status"], "checks_passed");
    assert_eq!(report["manifest_bytes"], 568);
    assert_eq!(
        report["manifest_sha256"],
        "235f49a6d07a5f7afb18a37f574ebc7ee1ae0b69dc9d4a33cd0b8ff15556b5b7"
    );
    assert_eq!(report["manifest_schema"], "ecbench.curve-corpus/v1");
}

#[test]
fn corpus_smallest_imports_and_solves_with_n_minus_one() {
    let family = fixture();
    let report = solve(
        &family.instance,
        family.target,
        SolveConfig {
            small_dimension: family.small_dimension,
            envelope_dimension: family.envelope_dimension,
            summands: family.instance.n - 1,
            max_large_primes: 2,
            max_trials: 20_000,
            max_states: 100_000,
            seed: 7,
        },
    )
    .unwrap();
    assert_eq!(report.recovered, Some(family.expected));
    assert!(report.verified_fast && report.verified_independent);
    assert!(report.verification_performed);
    assert!(report.counters.relation_histogram[2] > 0);
}

#[test]
fn exact_state_cap_refuses_without_downscaling() {
    let family = fixture();
    let err = solve(
        &family.instance,
        family.target,
        SolveConfig {
            small_dimension: family.small_dimension,
            envelope_dimension: family.envelope_dimension,
            summands: family.instance.n - 1,
            max_large_primes: 2,
            max_trials: 1,
            max_states: 10,
            seed: 7,
        },
    )
    .unwrap_err();
    assert!(err.contains("exact 10-summand MITM"));
}

#[test]
fn native_ecbench_runs_large_prime_and_matched_rho_on_one_workload() {
    let curve = CurveSpec::BinaryExplicit {
        n: 11,
        modulus: 2323,
        a: 1,
        b: 1,
        group_order: 1982,
        r: 991,
        gx: 1193,
        gy: 1007,
    };
    let instance = curve.build().unwrap();
    let workload = Workload::on(&curve, &instance, 20261004, 0, Default::default()).unwrap();
    let method = resolve(&MethodSpec {
        id: "ic.large_prime".into(),
        params: BTreeMap::from([
            ("small_dimension".into(), "3".into()),
            ("envelope_dimension".into(), "4".into()),
            ("summands".into(), "n-1".into()),
            ("large_primes".into(), "2".into()),
            ("max_trials".into(), "20000".into()),
            ("max_states".into(), "100000".into()),
        ]),
    })
    .unwrap();
    let ic = ecbench_solve(&method, &instance, &workload.curve, &workload.target, 7).unwrap();
    assert_eq!(ic.recovered, workload.planted);
    assert!(!ic.detail["verification_performed"].as_bool().unwrap());
    assert!(ic.unpriced.contains(&"row_ops_uncharged".to_string()));

    let rho = resolve(&MethodSpec {
        id: "rho.signed_frobenius".into(),
        params: BTreeMap::new(),
    })
    .unwrap();
    let generic = ecbench_solve(&rho, &instance, &workload.curve, &workload.target, 7).unwrap();
    assert_eq!(generic.recovered, workload.planted);
}
