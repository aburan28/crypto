//! Native replay and evidence sealing for the fixed F5 row-label screen.

use super::*;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("read JSON")).expect("parse JSON")
}

fn raw_lines(path: &Path) -> Vec<Value> {
    let input = fs::File::open(path).expect("raw file");
    BufReader::new(input)
        .lines()
        .map(|line| serde_json::from_str(&line.expect("raw line")).expect("raw JSON"))
        .collect()
}

fn compute_report(phase: &str, raw: &Path) -> Value {
    assert_eq!(sha256(PROTOCOL_BYTES), FROZEN_PROTOCOL_SHA256);
    let protocol: Protocol = serde_json::from_slice(PROTOCOL_BYTES).expect("protocol");
    let (split, seeds) = match phase {
        "discovery" => ("discovery", protocol.discovery_seeds),
        "holdout" => ("holdout", protocol.holdout_seeds),
        _ => panic!("phase"),
    };
    let raw_size = fs::metadata(raw).expect("raw size").len() as usize;
    assert!(raw_size <= protocol.evidence_cap_bytes);
    let lines = raw_lines(raw);
    let mut cursor = 0usize;
    let mut take = || {
        let next = lines
            .get(cursor)
            .unwrap_or_else(|| panic!("truncated at {cursor}"));
        cursor += 1;
        next
    };
    let expected_header = json!({
        "kind":"campaign", "phase":phase,
        "protocol_sha256":sha256(PROTOCOL_BYTES),
        "source_sha256":sha256(SOURCE_BYTES),
        "verifier_sha256":sha256(VERIFY_BYTES)
    });
    assert_eq!(take(), &expected_header);
    let mut cells = 0;
    let mut systems = 0;
    let mut selected_labels = 0u64;
    let mut criterion_word_xors = 0u64;
    let mut primary = Vec::new();
    let mut gate = true;
    for &n in &protocol.variables {
        for &seed in &seeds {
            for family in &protocol.families {
                for &batch in &protocol.batches {
                    let expected = cell_records(
                        n,
                        seed,
                        batch,
                        family,
                        protocol.selected_row_budget_per_system,
                        split,
                    );
                    for row in &expected {
                        assert_eq!(take(), row, "replay cell {n}/{seed}/{family}/{batch}");
                        if row["kind"] == "system" {
                            systems += 1;
                            selected_labels += row["signature"]["selected"].as_u64().unwrap();
                            criterion_word_xors +=
                                row["signature"]["criterion_word_xors"].as_u64().unwrap();
                        }
                    }
                    let summary = expected.last().unwrap();
                    if protocol.primary_variables.contains(&n) && batch == protocol.primary_batch {
                        let largest = summary["largest_class"].as_u64().unwrap() as usize;
                        let pass = 5 * largest >= 4 * batch;
                        assert_eq!(summary["gate"], pass);
                        gate &= pass;
                        primary.push(json!({
                            "n":n,"seed":seed,"family":family,"batch":batch,
                            "largest_class":largest,"distinct_signatures":summary["distinct_signatures"],
                            "adjacent_equal":summary["adjacent_equal"],
                            "base_equal":summary["base_equal"],"pass":pass
                        }));
                    }
                    cells += 1;
                }
            }
        }
    }
    assert_eq!(cursor, lines.len(), "trailing raw records");
    assert_eq!(cells, 72);
    assert_eq!(systems, 1008);
    assert_eq!(primary.len(), 18);
    assert_eq!(protocol.minimum_largest_signature_fraction, 0.8);
    json!({
        "schema_version":1,
        "phase":phase,
        "protocol_sha256":sha256(PROTOCOL_BYTES),
        "source_sha256":sha256(SOURCE_BYTES),
        "verifier_sha256":sha256(VERIFY_BYTES),
        "raw_sha256":sha256(&fs::read(raw).unwrap()),
        "raw_bytes":raw_size,
        "cells":cells,"systems":systems,
        "selected_labels":selected_labels,
        "criterion_word_xors":criterion_word_xors,
        "primary_cells":primary,
        "primary_passed":primary.iter().filter(|row| row["pass"] == true).count(),
        "gate_pass":gate,
        "complete_f5_cost":null,
        "full_ic_cost":null,
        "rho_ratio":null
    })
}

pub(super) fn verify_run(phase: &str, raw: &str, results: &str) {
    let output = Path::new(results);
    assert!(!output.exists(), "refuse to overwrite result");
    if phase == "holdout" {
        validate_binding(&output.with_file_name("binding.json"));
    }
    let report = compute_report(phase, Path::new(raw));
    fs::write(output, serde_json::to_vec_pretty(&report).unwrap()).expect("write result");
}

fn validate_binding(path: &Path) {
    let binding = read_json(path);
    assert_eq!(binding["protocol_sha256"], sha256(PROTOCOL_BYTES));
    assert_eq!(binding["source_sha256"], sha256(SOURCE_BYTES));
    assert_eq!(binding["verifier_sha256"], sha256(VERIFY_BYTES));
    for field in ["discovery_manifest_sha256", "discovery_results_sha256"] {
        assert_eq!(binding[field].as_str().unwrap().len(), 64);
    }
}

pub(super) fn seal(dir: &str) {
    let root = Path::new(dir);
    let manifest = root.join("manifest.json");
    assert!(!manifest.exists(), "refuse to overwrite manifest");
    let mut files = BTreeMap::new();
    for entry in fs::read_dir(root).expect("bundle directory") {
        let entry = entry.expect("directory member");
        let kind = entry.file_type().expect("file type");
        assert!(!kind.is_symlink(), "no symlinks");
        if !kind.is_file() {
            continue;
        }
        let name = entry.file_name().into_string().expect("UTF-8 member name");
        assert_ne!(name, "manifest.json");
        files.insert(name, sha256(&fs::read(entry.path()).expect("member bytes")));
    }
    for required in ["worker.rs", "verify.rs", "protocol.json", "raw.jsonl"] {
        assert!(files.contains_key(required), "missing {required}");
    }
    fs::write(
        manifest,
        serde_json::to_vec_pretty(&json!({
            "schema_version":1,"files":files
        }))
        .unwrap(),
    )
    .expect("write manifest");
}

pub(super) fn verify_bundle(dir: &str) -> Value {
    let root = Path::new(dir);
    let manifest = read_json(&root.join("manifest.json"));
    assert_eq!(manifest["schema_version"], 1);
    let files = manifest["files"].as_object().expect("file map");
    assert!(files.len() >= 8);
    for (name, digest) in files {
        assert!(!name.contains('/') && !name.contains(".."));
        assert_eq!(
            sha256(&fs::read(root.join(name)).expect("member")),
            digest.as_str().unwrap(),
            "member {name}"
        );
    }
    assert_eq!(fs::read(root.join("worker.rs")).unwrap(), SOURCE_BYTES);
    assert_eq!(fs::read(root.join("verify.rs")).unwrap(), VERIFY_BYTES);
    assert_eq!(
        fs::read(root.join("protocol.json")).unwrap(),
        PROTOCOL_BYTES
    );
    if root.join("failure.json").exists() {
        let failure = read_json(&root.join("failure.json"));
        assert_eq!(failure["performance_admitted"], false);
        return failure;
    }
    let old = read_json(&root.join("results.json"));
    let phase = old["phase"].as_str().expect("phase");
    if phase == "holdout" {
        validate_binding(&root.join("binding.json"));
    }
    let expected = compute_report(phase, &root.join("raw.jsonl"));
    assert_eq!(old, expected, "bundle replay differs");
    old
}

pub(super) fn check_discovery(bundle: &str, binding_path: &str) {
    let result = verify_bundle(bundle);
    assert_eq!(result["phase"], "discovery");
    assert_eq!(result["gate_pass"], true, "discovery screen failed");
    assert!(!Path::new(binding_path).exists());
    let root = Path::new(bundle);
    fs::write(
        binding_path,
        serde_json::to_vec_pretty(&json!({
            "protocol_sha256":sha256(PROTOCOL_BYTES),
            "source_sha256":sha256(SOURCE_BYTES),
            "verifier_sha256":sha256(VERIFY_BYTES),
            "discovery_manifest_sha256":sha256(&fs::read(root.join("manifest.json")).unwrap()),
            "discovery_results_sha256":sha256(&fs::read(root.join("results.json")).unwrap())
        }))
        .unwrap(),
    )
    .expect("write binding");
}

pub(super) fn failure(dir: &str, reason: &str) {
    let path = Path::new(dir).join("failure.json");
    assert!(!path.exists(), "refuse to overwrite failure");
    fs::write(
        path,
        serde_json::to_vec_pretty(&json!({
            "reason":reason,"performance_admitted":false,
            "complete_f5_cost":null,"full_ic_cost":null,"rho_ratio":null
        }))
        .unwrap(),
    )
    .expect("write failure");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn changed_signature_bit_is_detected() {
        let mut records = cell_records(12, 17, 2, "independent_affine", 8192, "development");
        let original = records[0].clone();
        records[0]["signature"]["selected_bits_hex"] = json!("00");
        assert_ne!(records[0], original);
        let original_digest = original["signature"]["selected_bits_sha256"]
            .as_str()
            .unwrap();
        assert_ne!(sha256(&hex::decode("00").unwrap()), original_digest);
    }

    #[test]
    fn sealed_failure_cannot_be_promoted() {
        let root = std::env::temp_dir().join(format!(
            "f5-signature-failure-{}-{}",
            std::process::id(),
            41
        ));
        assert!(!root.exists());
        fs::create_dir(&root).unwrap();
        fs::write(root.join("worker.rs"), SOURCE_BYTES).unwrap();
        fs::write(root.join("verify.rs"), VERIFY_BYTES).unwrap();
        fs::write(root.join("protocol.json"), PROTOCOL_BYTES).unwrap();
        fs::write(root.join("raw.jsonl"), b"").unwrap();
        for name in ["Cargo.toml", "Cargo.lock", "PROTOCOL.md", "README.md"] {
            fs::write(root.join(name), b"retained").unwrap();
        }
        failure(root.to_str().unwrap(), "producer-failed");
        seal(root.to_str().unwrap());
        let result = verify_bundle(root.to_str().unwrap());
        assert_eq!(result["performance_admitted"], false);
        fs::remove_dir_all(root).unwrap();
    }
}
