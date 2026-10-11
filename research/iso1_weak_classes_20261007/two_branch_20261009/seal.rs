//! Seal retained evidence with the already validated native SHA-256 implementation.
use crypto_lib::hash::sha256::Sha256;
use std::{fs, io::Read, path::Path};

fn files(base: &Path, at: &Path, out: &mut Vec<String>) {
    for entry in fs::read_dir(at).unwrap() {
        let entry = entry.unwrap();
        let kind = entry.file_type().unwrap();
        assert!(!kind.is_symlink(), "evidence must contain no symlinks");
        let path = entry.path();
        if kind.is_dir() {
            files(base, &path, out);
        } else {
            let name = path.strip_prefix(base).unwrap().to_str().unwrap();
            if !["EVIDENCE_SHA256.tsv", "REPORT.log"].contains(&name) {
                assert!(!name.contains(['\n', '\t']));
                out.push(name.to_string());
            }
        }
    }
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert!(args.len() == 1 || (args.len() == 2 && args[1] == "--check"));
    let base = Path::new(&args[0]);
    let mut names = Vec::new();
    files(base, base, &mut names);
    names.sort();
    let mut manifest = String::from("sha256\tbytes\tpath\n");
    for name in &names {
        let mut file = fs::File::open(base.join(name)).unwrap();
        let mut hash = Sha256::new();
        let mut buffer = [0u8; 65536];
        let mut bytes = 0u64;
        loop {
            let n = file.read(&mut buffer).unwrap();
            if n == 0 {
                break;
            }
            hash.update(&buffer[..n]);
            bytes += n as u64;
        }
        manifest.push_str(&format!(
            "{}\t{bytes}\t{name}\n",
            hex::encode(hash.finalize())
        ));
    }
    let path = base.join("EVIDENCE_SHA256.tsv");
    if args.len() == 2 {
        assert_eq!(
            fs::read_to_string(path).unwrap(),
            manifest,
            "evidence seal mismatch"
        );
        println!("EVIDENCE_SEAL_VERIFIED|files={}", names.len());
    } else {
        fs::write(path, manifest).unwrap();
        println!("EVIDENCE_SEAL_CREATED|files={}", names.len());
    }
}
