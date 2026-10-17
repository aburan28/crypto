use std::process::Command;
#[test]
fn full_norm_predicate_matches_literal_norm_for_every_tiny_assignment() {
 let out=Command::new(env!("CARGO_BIN_EXE_full-binary-norm-circuit")).arg("controls").output().unwrap();
 assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
 let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
 assert_eq!(v["status"],"passed");assert_eq!(v["independent_full_norm_assignments"],73216);assert_eq!(v["models"],8);
}
