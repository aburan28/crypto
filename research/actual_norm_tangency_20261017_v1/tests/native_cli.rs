use std::process::Command;
#[test]
fn closed_tangency_fibers_equal_exhaustive_norm_derivative_roots() {
 let out=Command::new(env!("CARGO_BIN_EXE_actual-norm-tangency")).arg("controls").output().unwrap();
 assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
 let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
 assert_eq!(v["status"],"passed");assert_eq!(v["tangency_exhaustive_parameter_checks"],15664);assert_eq!(v["artin_schreier_rhs_checks"],28);
}
