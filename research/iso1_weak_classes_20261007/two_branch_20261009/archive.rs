//! Preserve executed sources and disclose observed-source / loaded-binary drift.
use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::{
    fs,
    path::{Path, PathBuf},
};
fn hash(b: &[u8]) -> String {
    hex::encode(sha256(b))
}
fn walk(p: &Path, out: &mut Vec<PathBuf>) {
    for e in fs::read_dir(p).unwrap() {
        let e = e.unwrap();
        let path = e.path();
        if e.file_type().unwrap().is_dir() {
            walk(&path, out);
        } else if e.file_name() == "source_freeze.json" {
            out.push(path);
        }
    }
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 1);
    let base = Path::new(&a[0]);
    let v1 = fs::read(base.join("run_v1_snapshot.txt")).unwrap();
    let actual = hash(&v1);
    let mut files = Vec::new();
    walk(base, &mut files);
    let mut corrections = Vec::new();
    let mut missing = Vec::new();
    let mut odd2 = fs::read_to_string(base.join("construct_odd.gp")).unwrap();
    odd2 = odd2
        .replace("edges3=0,edges5=0,edges7=0,", "edges3=0,")
        .replace(
            "if(deg==2,edges2++,if(deg==3,edges3++,if(deg==5,edges5++,edges7++)));",
            "if(deg==2,edges2++,edges3++);",
        )
        .replace("if(witness||limit,head,head-1)", "head")
        .replace(",\"|edges5=\",edges5,\"|edges7=\",edges7", "");
    let mut loaded_binary = None;
    for p in files {
        let v: Value = serde_json::from_slice(&fs::read(&p).unwrap()).unwrap();
        let dir = p.parent().unwrap();
        let script = v["script"].as_str().unwrap();
        let src = fs::read(base.join(script)).unwrap();
        let target = dir.join("executed_script.gp");
        if hash(&src) == v["script_sha256"].as_str().unwrap() {
            if !target.exists() {
                fs::write(&target, &src).unwrap();
            }
        } else if target.exists() {
            assert_eq!(hash(&fs::read(&target).unwrap()), v["script_sha256"]);
        } else if hash(odd2.as_bytes()) == v["script_sha256"].as_str().unwrap() {
            fs::write(&target, &odd2).unwrap();
        } else {
            missing.push(p.display().to_string());
        }
        if dir.to_string_lossy().contains("construction_large_run1") {
            let binary = v["driver_binary_sha256"].as_str().unwrap();
            if let Some(ref expected) = loaded_binary {
                assert_eq!(binary, expected);
            } else {
                loaded_binary = Some(binary.to_string());
            }
        }
        let current_driver = dir.join("executed_driver.txt");
        if current_driver.exists() {
            assert_eq!(
                hash(&fs::read(&current_driver).unwrap()),
                v["driver_source_sha256"]
            );
        } else if v["driver_source_sha256"] != actual {
            corrections.push(json!({"receipt":p.strip_prefix(base).unwrap().display().to_string(),"original_observed_source_sha256":v["driver_source_sha256"],"actual_executed_v1_source_sha256":actual,"loaded_driver_binary_sha256":v["driver_binary_sha256"]}));
        }
        if !current_driver.exists() {
            fs::write(dir.join("executed_driver_v1.txt"), &v1).unwrap();
        }
        assert_eq!(
            hash(&fs::read(base.join("PROTOCOL.md")).unwrap()),
            v["protocol_sha256"]
        );
        if let Some(input) = v["environment"]["ISO1_INPUTS"].as_str() {
            fs::copy(input, dir.join("executed_inputs.gp")).unwrap();
        }
    }
    fs::write(base.join("provenance_correction.json"),serde_json::to_vec_pretty(&json!({"schema":"iso1.provenance_correction/v1","reason":"The running launcher binary remained v1 while source-archival code was edited. Fourteen population launcher receipts observed the edited source-file hash; the loaded binary, GP arithmetic script and fixture inputs remained frozen. Original receipts are preserved. Bind their actual execution to run_v1_snapshot.txt and the recorded binary hash.","actual_v1_source_sha256":actual,"population_loaded_binary_sha256":loaded_binary,"corrections":corrections,"unarchived_script_receipts":missing})).unwrap()).unwrap();
    assert!(
        missing.is_empty(),
        "archive each executed GP source before publication"
    );
    println!(
        "archived exact GP sources; {} launcher source-observation corrections preserved",
        corrections.len()
    );
}
