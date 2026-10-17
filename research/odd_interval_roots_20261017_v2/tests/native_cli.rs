use std::process::Command;
#[test]
fn odd_roots_and_multiplicities_match_independent_hasse_oracle() {
 let out=Command::new(env!("CARGO_BIN_EXE_odd-interval-roots")).arg("controls").output().unwrap();
 assert!(out.status.success(),"{}",String::from_utf8_lossy(&out.stderr));
 let v:serde_json::Value=serde_json::from_slice(&out.stdout).unwrap();
 assert_eq!(v["status"],"passed");assert_eq!(v["monic_quintic_cases"],3624);
}
