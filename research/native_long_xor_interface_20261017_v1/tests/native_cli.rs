#[test]
#[ignore = "requires the source-pinned modified native backend"]
fn native_long_xor_model_sets() {
 let lib=std::env::var("ROOTFINDING_CMS_LIBRARY").expect("source-pinned long-XOR wrapper path");
 let out=std::process::Command::new(env!("CARGO_BIN_EXE_native-long-xor-interface")).args(["integration-controls",&lib]).output().unwrap();
 assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
 let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
 assert_eq!(v["status"],"passed");assert_eq!(v["evidence"]["controls"]["complete_primary_models"],128);
}
