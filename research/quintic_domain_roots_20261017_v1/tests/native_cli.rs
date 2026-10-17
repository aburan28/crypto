use std::process::Command;

#[test]
fn tiny_domains_match_independent_hasse_oracle() {
    let output = Command::new(env!("CARGO_BIN_EXE_quintic-domain-roots"))
        .arg("controls").output().expect("native controls process");
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let result: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(result["status"], "passed");
    assert_eq!(result["independent_domain_root_cases"], 2816);
    assert_eq!(result["exact_quotient_square_vectors"], 1536);
}

#[test]
fn wide_implicit_domains_retain_every_planted_root_and_multiplicity() {
    let output = Command::new(env!("CARGO_BIN_EXE_quintic-domain-roots"))
        .arg("native-domain-root-unit-control").output().expect("native wide control");
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let result: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(result["status"], "passed");
    assert_eq!(result["evidence"]["field_count"], 6);
    assert_eq!(result["evidence"]["domain_enumerated"], false);
    assert_eq!(result["evidence"]["release_gate_claim"], false);
    let cases = result["evidence"]["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 36);
    assert!(cases.iter().all(|row| row["component_complete"] == true));
    assert!(cases.iter().all(|row| row["nodes"].as_u64().unwrap() <= 217));
}
