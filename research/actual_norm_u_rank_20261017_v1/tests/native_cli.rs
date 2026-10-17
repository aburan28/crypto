use std::process::Command;
#[test]
fn literal_norm_and_exhaustive_tiny_families_reject_missing_term() {
    let out=Command::new(env!("CARGO_BIN_EXE_actual-norm-u-rank")).arg("controls").output().unwrap();
    assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
    let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["status"],"passed");
    assert_eq!(v["literal_norm_cases"],2304);
    assert_eq!(v["tiny_family_memberships"],256);
    assert!(v["old_transfer_formula_counterexamples"].as_u64().unwrap()>0);
}
