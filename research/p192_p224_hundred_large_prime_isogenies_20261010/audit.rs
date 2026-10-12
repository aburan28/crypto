use serde_json::{json, Value};
use std::{fs, path::Path};
#[path = "../../src/hash/sha256.rs"]
#[allow(dead_code)]
mod sha256;
fn hash(b: &[u8]) -> String {
    sha256::sha256(b)
        .iter()
        .map(|x| format!("{x:02x}"))
        .collect()
}
fn main() {
    let repo = std::env::current_dir().unwrap();
    let study = repo.join("research/p192_p224_hundred_large_prime_isogenies_20261010");
    let m: Value = serde_json::from_slice(&fs::read(study.join("MANIFEST.json")).unwrap()).unwrap();
    let mut checked = 0usize;
    for (key, base) in [
        ("files", study.as_path()),
        ("canonical_outputs", repo.as_path()),
    ] {
        for r in m[key].as_array().unwrap() {
            let p = r["path"].as_str().unwrap();
            assert!(
                !Path::new(p).is_absolute()
                    && !Path::new(p)
                        .components()
                        .any(|c| matches!(c, std::path::Component::ParentDir))
            );
            let b = fs::read(base.join(p)).unwrap();
            assert_eq!(r["bytes"], b.len(), "{p}");
            assert_eq!(r["sha256"], hash(&b), "{p}");
            checked += 1;
        }
    }
    let s: Value = serde_json::from_slice(&fs::read(study.join("SUMMARY.json")).unwrap()).unwrap();
    assert_eq!(s["new_distinct_maps"], 100);
    assert_eq!(s["distinct_prime_degrees"], 2);
    let mut seen = std::collections::HashSet::new();
    for r in s["records"].as_array().unwrap() {
        assert!(seen.insert(format!("{}:{}", r["source_icv1"], r["certificate_sha256"])));
    }
    assert_eq!(seen.len(), 100);
    println!("{}",serde_json::to_string_pretty(&json!({"schema":"hundred-prime-artifact-audit/v1","status":"PASS","files_checked":checked,"distinct_maps":100,"distinct_prime_degrees":2,"manifest_sha256":hash(&fs::read(study.join("MANIFEST.json")).unwrap())})).unwrap());
}
