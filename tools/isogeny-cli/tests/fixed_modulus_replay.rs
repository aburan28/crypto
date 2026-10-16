//! Bounded exact controls and public cache-path admission, no timing claim.
#![allow(dead_code)]
#[path = "../src/verification_field.rs"]
mod verification_field;
mod field {
    pub type Field = crate::verification_field::PrimeField<9>;
    pub type Fe = crate::verification_field::Element<9>;
    pub const FIELD_BITS: usize = 576;
}
#[path = "../../../src/cryptanalysis/isogeny_walk/curve.rs"]
mod curve;
#[path = "../../../src/cryptanalysis/identity_certificate.rs"]
mod identity_certificate;
#[path = "../../../research/p192_p224_large_degree_isogenies_20261009/verification_kernel.rs"]
mod kernel;
#[path = "../src/verification_modpoly.rs"]
mod modpoly;
#[path = "../../../research/p192_p224_large_degree_isogenies_20261009/verification_poly.rs"]
mod poly;
mod hash {
    pub mod sha256 {
        pub fn sha256(bytes: &[u8]) -> [u8; 32] {
            use sha2::Digest;
            sha2::Sha256::digest(bytes).into()
        }
    }
}
#[cfg(test)]
#[path = "../../../src/cryptanalysis/modular_polynomial.rs"]
pub mod modular_polynomial;
#[cfg(test)]
pub mod cryptanalysis {
    pub use crate::modular_polynomial;
}
#[cfg(test)]
#[path = "../src/catalog.rs"]
mod catalog;
#[cfg(test)]
const REGISTRY: &str = include_str!("../../../docs/curves/registry.json");

use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::Command,
};
fn sha(b: &[u8]) -> String {
    hex::encode(Sha256::digest(b))
}
fn save(p: &Path, v: &Value) {
    let mut f = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(p)
        .unwrap();
    f.write_all(&serde_json::to_vec_pretty(v).unwrap()).unwrap();
}
fn read(p: &Path) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn root() -> (PathBuf, bool) {
    let retained = std::env::var_os("ISOGENY_FIXED_MODULUS_EVIDENCE");
    let p = retained.as_ref().map(PathBuf::from).unwrap_or_else(|| {
        std::env::temp_dir().join(format!(
            "fixed-modulus-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    });
    fs::create_dir(&p).unwrap();
    (fs::canonicalize(p).unwrap(), retained.is_some())
}
#[test]
fn p521_arithmetic_and_kernel_reference_controls() {
    let p = (num_bigint::BigUint::from(1u64) << 521usize) - num_bigint::BigUint::from(1u64);
    let f = field::Field::new(&p).unwrap();
    poly::fixed_modulus_control_suite(&f);
    kernel::fixed_modulus_kernel_controls(&f);
}
#[test]
fn newton_identity_native_issue_recheck_and_optional_export() {
    use identity_certificate::{self as idc, Expr as E, IdentityStatement};
    let b = E::var("B");
    let g = E::var("G");
    let bg = E::times(b, g);
    let error = E::plus(E::Const(1), E::minus(bg.clone()));
    let statement = IdentityStatement {
        name: "Newton reciprocal update residual".into(),
        source: "fixed-modulus-newton.md; scalar supplement to formal-series reduction proof"
            .into(),
        variables: vec!["B".into(), "G".into()],
        degree_bound: 4,
        lhs: E::times(bg.clone(), E::plus(E::Const(2), E::minus(bg))),
        rhs: E::plus(E::Const(1), E::minus(E::times(error.clone(), error))),
    };
    let certificate=idc::issue(statement,20261016,32,
        "native issuance and same-process checker plus domain-separated extra points; independent source review pending",None).unwrap();
    assert!(matches!(
        idc::check_with_extra(&certificate, Some((0x7265636970726f63, 32))).unwrap(),
        idc::Verdict::Accept { .. }
    ));
    if let Some(path) = std::env::var_os("ISOGENY_NEWTON_IDENTITY_OUTPUT") {
        save(
            Path::new(&path),
            &serde_json::to_value(&certificate).unwrap(),
        );
    }
}
#[test]
fn accepted_p224_1471_cache_control_and_rehashed_negatives() {
    let (root, retained) = root();
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let source = repo.join("research/p192_p224_large_degree_isogenies_20261009/kernel-p224-1471");
    let source_raw = fs::read(source.join("p224/ell-1471.json")).unwrap();
    assert_eq!(
        sha(&source_raw),
        "ecf1109dbff2a72000fd84d9939127db68e2287e55e4ec10d143229874c230c8"
    );
    let doc: Value = serde_json::from_slice(&source_raw).unwrap();
    let p = num_bigint::BigUint::parse_bytes(doc["curve"]["p"].as_str().unwrap().as_bytes(), 10)
        .unwrap();
    let f = field::Field::new(&p).unwrap();
    let h: poly::Poly = doc["result"]["isogenies"][0]["kernel"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| {
            f.from_big(
                &num_bigint::BigUint::parse_bytes(v.as_str().unwrap().as_bytes(), 10).unwrap(),
            )
        })
        .collect();
    let context = poly::FixedModulus::new(&f, &h).unwrap();
    assert_eq!(h.len(), 736);
    assert!(context.has_reciprocal());
    let product = poly::mul(&f, &h[..735].to_vec(), &h[..735].to_vec());
    assert!(product.len() >= h.len() && product.len() < 2 * h.len());
    assert_eq!(context.rem(&product), poly::divrem(&f, &product, &h).1);
    let binary = std::env::var_os("ISOGENY_CONTROL_BINARY")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_BIN_EXE_isogeny")));
    save(
        &root.join("cache-activation-diagnostic.json"),
        &json!({"schema":"fixed-modulus-diagnostic/v1",
        "kernel_length":h.len(),"stored_reciprocal":context.has_reciprocal(),"eligible_product_length":product.len(),
        "compiled_feature":cfg!(feature="cached-fixed-modulus-reduction"),"per_operation_instrumentation":false,
        "input_sha256":sha(&source_raw),"binary_sha256":sha(&fs::read(&binary).unwrap()),"timing_claim":false}),
    );
    let mut results = vec![];
    for kind in ["valid", "kernel", "numerator"] {
        let input = root.join(format!("{kind}-input"));
        fs::create_dir(&input).unwrap();
        fs::create_dir(input.join("p224")).unwrap();
        let mut document = doc.clone();
        if kind != "valid" {
            let v = if kind == "kernel" {
                &mut document["result"]["isogenies"][0]["kernel"][0]
            } else {
                &mut document["result"]["isogenies"][0]["x_map"]["num"][0]
            };
            *v = if *v == "1" { json!("2") } else { json!("1") };
        }
        let map = input.join("p224/ell-1471.json");
        if kind == "valid" {
            fs::write(&map, &source_raw).unwrap();
        } else {
            save(&map, &document);
        }
        fs::copy(
            source.join("p224/ell-1471.stderr.txt"),
            input.join("p224/ell-1471.stderr.txt"),
        )
        .unwrap();
        let mut summary = read(&source.join("search.json"));
        summary["attempts"][0]["receipt"]["stdout_sha256"] = json!(sha(&fs::read(&map).unwrap()));
        save(&input.join("search.json"), &summary);
        let out = root.join(format!("{kind}-judge"));
        let trace = root.join(format!("{kind}.jsonl"));
        let args = [
            "challenge",
            "judge",
            "--input",
            input.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
        ];
        let output = Command::new(&binary)
            .args(args)
            .env("ISOGENY_REPLAY_TRACE", &trace)
            .output()
            .unwrap();
        fs::write(root.join(format!("{kind}.stdout")), &output.stdout).unwrap();
        fs::write(root.join(format!("{kind}.stderr")), &output.stderr).unwrap();
        let result: Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(
            output.status.code(),
            Some(if kind == "valid" { 0 } else { 1 }),
            "{result}"
        );
        assert_eq!(
            result["status"],
            if kind == "valid" {
                "ACCEPTED"
            } else {
                "REJECTED"
            }
        );
        if kind == "valid" {
            let records: Vec<Value> = fs::read_to_string(&trace)
                .unwrap()
                .lines()
                .map(|s| serde_json::from_str(s).unwrap())
                .collect();
            let setup: Vec<_> = records
                .iter()
                .filter(|v| v["stage"] == "fixed_modulus_context_setup")
                .collect();
            assert_eq!(
                setup.len(),
                if cfg!(feature = "cached-fixed-modulus-reduction") {
                    2
                } else {
                    0
                }
            );
            if !setup.is_empty() {
                assert_eq!(setup[0]["event"], "begin");
                assert_eq!(setup[1]["event"], "end");
            }
            assert_eq!(records.last().unwrap()["status"], "COMPLETED");
        } else {
            assert!(result["score"].is_null());
        }
        results.push(
            json!({"case":kind,"exit":output.status.code(),"status":result["status"],"entry":result,
            "command":args,"input_sha256":sha(&fs::read(&map).unwrap())}),
        );
    }
    save(
        &root.join("controls.json"),
        &json!({"schema":"fixed-modulus-public-controls/v1","controls":results,
        "mathematical_exclusions":[],"large_target_executed":false}),
    );
    if !retained {
        fs::remove_dir_all(root).unwrap();
    }
}

#[test]
fn companion_waits_for_ack_and_marks_done_even_on_parse_failure() {
    use std::{
        process::Stdio,
        time::{Duration, Instant},
    };
    let root = std::env::temp_dir().join(format!(
        "fixed-modulus-worker-control-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&root).unwrap();
    let bench = env!("CARGO_BIN_EXE_fixed-modulus-bench");
    for bad_done in [false, true] {
        let dir = root.join(if bad_done {
            "writer-failure"
        } else {
            "parse-failure"
        });
        fs::create_dir(&dir).unwrap();
        let ready = dir.join("resources-ready");
        let mut child = Command::new(bench)
            .args([
                "fixed-modulus-bench",
                "__worker",
                "measure",
                "--input",
                dir.join("missing-dataset.json").to_str().unwrap(),
                "--out",
                dir.to_str().unwrap(),
                "--case",
                "control",
                "--rss-ready",
                ready.to_str().unwrap(),
            ])
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        std::thread::sleep(Duration::from_millis(50));
        assert!(
            child.try_wait().unwrap().is_none(),
            "worker parsed before RSS acknowledgement"
        );
        assert!(!ready.with_extension("done").exists());
        if bad_done {
            fs::write(ready.with_extension("done"), b"sentinel").unwrap();
        }
        fs::write(&ready, b"").unwrap();
        let start = Instant::now();
        while child.try_wait().unwrap().is_none() {
            if start.elapsed() > Duration::from_secs(5) {
                child.kill().unwrap();
                panic!("worker failed to terminate after acknowledgement");
            }
            std::thread::sleep(Duration::from_millis(5));
        }
        let output = child.wait_with_output().unwrap();
        assert_eq!(output.status.code(), Some(if bad_done { 2 } else { 1 }));
        assert!(ready.with_extension("done").exists());
        if bad_done {
            assert_eq!(fs::read(ready.with_extension("done")).unwrap(), b"sentinel");
        } else {
            assert_eq!(
                serde_json::from_slice::<Value>(&output.stdout).unwrap()["status"],
                "FAIL"
            );
        }
    }
    // Actual parent/child native shape must acknowledge before malformed parsing.
    let invalid = root.join("invalid.json");
    fs::write(&invalid, b"{").unwrap();
    let out = root.join("supervised-malformed");
    let output = Command::new(bench)
        .args([
            "measure",
            "--input",
            invalid.to_str().unwrap(),
            "--case",
            "control",
            "--out",
            out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    let receipt = read(&out.join("native.receipt.json"));
    assert_eq!(receipt["status"], "FAIL");
    assert_eq!(receipt["unavailable_memory_samples"], 0);
    assert!(
        receipt["sampled_peak_resident_bytes"]
            .as_str()
            .unwrap()
            .parse::<u64>()
            .unwrap()
            > 0
    );
    assert!(out.join("resources-ready").is_file() && out.join("resources-ready.done").is_file());
    fs::remove_dir_all(root).unwrap();
}
