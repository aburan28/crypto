#[test]
#[ignore = "requires the source-pinned modified native wrapper"]
fn literal_all_models_and_budget_boundaries() {
 let lib=std::env::var("ROOTFINDING_CMS_LIBRARY").unwrap();
 let out=std::process::Command::new(env!("CARGO_BIN_EXE_hybrid-norm-search")).args(["controls",&lib]).output().unwrap();
 assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
 let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
 assert_eq!(v["status"],"passed");assert_eq!(v["independent_primary_assignments"],2048);assert_eq!(v["complete_native_sets"].as_array().unwrap().len(),6);
 assert_eq!(v["model_cap_is_inconclusive"],true);assert_eq!(v["expired_budget_is_inconclusive"],true);
}
