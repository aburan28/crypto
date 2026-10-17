use std::process::Command;
#[test]
fn affine_cell_solutions_match_literal_exhaustive_membership() {
    let out=Command::new(env!("CARGO_BIN_EXE_quintic-fixed-u-cells")).arg("controls").output().unwrap();
    assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
    let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["status"],"passed");
    assert_eq!(v["cells"],150);
    assert_eq!(v["independent_parameter_memberships"],9600);
}
