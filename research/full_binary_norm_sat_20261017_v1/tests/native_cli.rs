use std::process::Command;
#[test]
fn norm_math_matches_independent_literal_oracle() {
 let out=Command::new(env!("CARGO_BIN_EXE_full-binary-norm-sat")).arg("controls").output().unwrap();assert!(out.status.success());let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();assert_eq!(v["independent_full_norm_assignments"],73216);
}
#[test]
#[ignore="requires the frozen local modified CMS library; explicitly run with ROOTFINDING_CMS_LIBRARY"]
fn native_backend_matches_exhaustive_tiny_root_set() {
 let p=std::env::var("ROOTFINDING_CMS_LIBRARY").unwrap();let out=Command::new(env!("CARGO_BIN_EXE_full-binary-norm-sat")).args(["native-controls",&p]).output().unwrap();assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();assert_eq!(v["native_tiny_primary_models"],28);assert_eq!(v["native_complete"],true);
}
