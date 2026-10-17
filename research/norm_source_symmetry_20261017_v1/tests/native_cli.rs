#[test]
fn restores_every_ordered_literal_norm_root() {
 let out=std::process::Command::new(env!("CARGO_BIN_EXE_norm-source-symmetry")).output().unwrap();
 assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
 let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
 assert_eq!(v["status"],"passed");assert_eq!(v["literal_primary_assignments"],34816);
 assert_eq!(v["uniform_division_by120_refuted"],true);assert_eq!(v["cases"].as_array().unwrap().len(),8);
}
