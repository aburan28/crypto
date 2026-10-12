use std::path::PathBuf;
use std::process::Command;

use serde_json::Value;

fn scratch() -> PathBuf {
    let path = std::env::temp_dir().join(format!("p256-factor-base-test-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&path);
    std::fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn a_wide_dump_verifies_and_writes_repository_sql() {
    let dir = scratch();
    let dump = dir.join("p256.factor-base.json");
    let sql = dir.join("p256.factor-base.sql");
    let output = Command::new(env!("CARGO_BIN_EXE_p256_factor_base"))
        .args([
            "--curve",
            "icv1-fp256-t89188191154553853111372247798585809583-f188c491",
            "--factor-base",
            "dickson-torus:depth=3",
            "--out",
            dump.to_str().unwrap(),
            "--sql-out",
            sql.to_str().unwrap(),
            "--relation-length",
            "2",
            "--verify",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("(verified)"));

    let value: Value = serde_json::from_str(&std::fs::read_to_string(&dump).unwrap()).unwrap();
    assert_eq!(value["schema"], "ecbench.factor_base_dump/v1-wide");
    assert_eq!(value["curve"]["ec1"], "EC1P256Cp256h0523b774e066");
    assert_eq!(
        value["points"].as_array().unwrap().len() as u64,
        value["factor_base"]["signed_points"].as_u64().unwrap()
    );

    let sql = std::fs::read_to_string(sql).unwrap();
    assert!(sql.contains("INSERT INTO curves VALUES"));
    assert!(sql.contains("INSERT INTO curve_representations VALUES"));
    assert_eq!(
        sql.matches("INSERT INTO factor_base_points VALUES").count(),
        value["factor_base"]["signed_points"].as_u64().unwrap() as usize
    );
    let _ = std::fs::remove_dir_all(dir);
}
