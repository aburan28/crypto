use std::process::Command;
#[test]
fn repeated_pair_equivalence_retains_zero_and_nonrational_fibers() {
    let out=Command::new(env!("CARGO_BIN_EXE_actual-norm-repeated-pair")).arg("controls").output().unwrap();
    assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
    let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["status"],"passed");
    assert_eq!(v["literal_norm_cases"],2304);
    assert_eq!(v["tiny_family_memberships"],256);
    assert_eq!(v["repeated_pair_divisibility_controls"],144);
}
