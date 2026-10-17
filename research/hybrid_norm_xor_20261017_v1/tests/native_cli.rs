#[test]
#[ignore = "requires the source-pinned native long-XOR wrapper"]
fn exact_hybrid_norm_model_sets() {
 let lib=std::env::var("ROOTFINDING_CMS_LIBRARY").unwrap();
 let out=std::process::Command::new(env!("CARGO_BIN_EXE_hybrid-norm-xor")).args(["controls",&lib]).output().unwrap();
 assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
 let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
 assert_eq!(v["status"],"passed");assert_eq!(v["hybrid_complete_primary_assignments"],116480);assert_eq!(v["native_complete_model_set_arms"],7);
}
