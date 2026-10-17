use std::process::Command;
#[test]
fn affine_forms_preserve_all_tiny_primary_assignments() {
 let out=Command::new(env!("CARGO_BIN_EXE_norm-xor-cones")).arg("controls").output().unwrap();assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();assert_eq!(v["independent_full_norm_assignments"],73216);assert_eq!(v["affine_cone_complete_primary_assignments"],16384);
}
