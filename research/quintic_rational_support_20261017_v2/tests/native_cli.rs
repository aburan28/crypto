use std::process::Command;
#[test]
fn all_seven_exact_branches_match_independent_oracles() {
    let out=Command::new(env!("CARGO_BIN_EXE_quintic-rational-support")).arg("controls").output().unwrap();
    assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
    let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["status"],"passed");
    assert_eq!(v["support_oracle_cases"],3072);
    assert_eq!(v["rational_membership_cases"],616);
    assert_eq!(v["odd_prime_extension_rational_cases"],3400);
    for branch in ["squarefree","repeated_cubic","repeated_outside_reject","exceptional_quadratic","exceptional_linear_outside","repeated_zero_cubic","repeated_zero_outside"] {
        assert!(v["branches"][branch].as_u64().unwrap()>0);
    }
}
