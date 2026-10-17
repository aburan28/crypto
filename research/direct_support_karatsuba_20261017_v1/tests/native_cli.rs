use std::process::Command;
#[test]
fn recursive_and_unbalanced_products_match_literal_convolution() {
 let out=Command::new(env!("CARGO_BIN_EXE_direct-support-karatsuba")).arg("controls").output().unwrap();
 assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
 let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
 assert_eq!(v["status"],"passed");assert_eq!(v["independent_literal_coefficient_cases"],27648);assert_eq!(v["independent_adaptive_multiplication_vectors"],288);
}
