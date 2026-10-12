//! Bind the proof addition to its exact sources, receipts and frozen model input.
#[allow(dead_code)]
#[path = "../../src/hash/sha256.rs"]
mod sha256;
use std::{fs, path::Path};
const SEAL: &str = "EVIDENCE_SHA256.tsv";
fn digest(bytes: &[u8]) -> String {
    sha256::sha256(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
fn collect(base: &Path, at: &Path, rows: &mut Vec<(String, String)>) {
    for entry in fs::read_dir(at).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();
        let name = entry.file_name();
        let name = name.to_str().unwrap();
        if name == ".lake" || name == SEAL {
            continue;
        }
        let meta = fs::symlink_metadata(&path).unwrap();
        assert!(!meta.file_type().is_symlink(), "unexpected source symlink");
        if meta.is_dir() {
            collect(base, &path, rows);
        } else {
            let relative = path.strip_prefix(base).unwrap().to_str().unwrap();
            rows.push((relative.to_string(), digest(&fs::read(path).unwrap())));
        }
    }
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(
        args.len(),
        1,
        "--write or --check; run in the proof directory"
    );
    let input = "../../docs/curves/weak-families/models.json";
    assert_eq!(
        digest(&fs::read(input).unwrap()),
        "9738f6bd6c016732f88e4e7979ac09e9aa6fd16eba96220ec8781c0017b6055e",
        "frozen numerical inventory changed"
    );
    let mut rows = vec![];
    collect(Path::new("."), Path::new("."), &mut rows);
    for path in [
        "../../src/hash/sha256.rs",
        "../../.github/workflows/weak-curve-lean.yml",
    ] {
        rows.push((path.to_string(), digest(&fs::read(path).unwrap())));
    }
    rows.sort();
    let output: String = rows.iter().map(|(p, h)| format!("{h}\t{p}\n")).collect();
    match args[0].as_str() {
        "--write" => fs::write(SEAL, &output).unwrap(),
        "--check" => assert_eq!(
            fs::read_to_string(SEAL).unwrap(),
            output,
            "proof evidence changed"
        ),
        _ => panic!("unknown operation"),
    }
    println!("PROOF_EVIDENCE_VERIFIED|files={}", rows.len());
}
