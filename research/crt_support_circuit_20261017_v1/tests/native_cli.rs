use std::process::Command;
#[test]
fn exact_crt_matches_literal_products_and_rejects_insufficient_degree() {
 let out=Command::new(env!("CARGO_BIN_EXE_crt-support-circuit")).arg("controls").output().unwrap();
 assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
 let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
 assert_eq!(v["status"],"passed");assert_eq!(v["independent_literal_coefficient_cases"],27648);assert_eq!(v["insufficient_CRT_degree_counterexample"],true);
}
