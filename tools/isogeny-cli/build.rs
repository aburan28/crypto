use sha2::{Digest, Sha256};
use std::{env, fs, path::Path, process::Command};

fn git(repo: &Path, args: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(repo)
        .args(args)
        .output()
        .ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())
}

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = root.canonicalize().unwrap();
    let registry_path = root.join("docs/curves/registry.json");
    println!("cargo:rerun-if-changed={}", registry_path.display());
    let registry_bytes = fs::read(&registry_path).unwrap();
    let registry: serde_json::Value = serde_json::from_slice(&registry_bytes).unwrap();
    let rows: Vec<_> = registry["curves"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            let mut entry = serde_json::Map::new();
            for key in [
                "icv1",
                "slug",
                "family",
                "params",
                "order",
                "aliases",
                "standard_names",
                "representations",
            ] {
                entry.insert(key.to_owned(), row[key].clone());
            }
            serde_json::Value::Object(entry)
        })
        .collect();
    let index = serde_json::json!({"curves": rows});
    fs::write(
        Path::new(&env::var("OUT_DIR").unwrap()).join("execution_catalogue.json"),
        serde_json::to_vec(&index).unwrap(),
    )
    .unwrap();
    let checksum: String = Sha256::digest(&registry_bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    println!("cargo:rustc-env=ISOGENY_REGISTRY_SHA256={checksum}");
    let mut engines = String::new();
    for limbs in 1..=10 {
        engines.push_str(&format!("mod n{limbs} {{\nmod field {{ pub type Field=crate::verification_field::PrimeField<{limbs}>; pub type Fe=crate::verification_field::Element<{limbs}>; pub const FIELD_BITS:usize={limbs}*64; }}\n"));
        for (name, path) in [
            ("curve", "src/cryptanalysis/isogeny_walk/curve.rs"),
            ("modpoly", "tools/isogeny-cli/src/verification_modpoly.rs"),
            (
                "poly",
                "research/p192_p224_large_degree_isogenies_20261009/verification_poly.rs",
            ),
            (
                "kernel",
                "research/p192_p224_large_degree_isogenies_20261009/verification_kernel.rs",
            ),
            ("map_identity", "tools/isogeny-cli/src/verification_map.rs"),
            ("replay", "tools/isogeny-cli/src/replay.rs"),
        ] {
            let path = root.join(path);
            println!("cargo:rerun-if-changed={}", path.display());
            engines.push_str(&format!("#[path = {path:?}] mod {name};\n"));
        }
        engines.push_str("pub fn verify(out:&std::path::Path,write:bool)->serde_json::Value { replay::verify(out,write) }\n}\n");
    }
    engines.push_str("pub fn verify(out:&std::path::Path,write:bool)->serde_json::Value {\nlet summary:serde_json::Value=serde_json::from_slice(&std::fs::read(out.join(\"search.json\")).unwrap()).unwrap();\nlet item=&summary[\"attempts\"][0];\nlet data:serde_json::Value=serde_json::from_slice(&std::fs::read(out.join(item[\"curve\"].as_str().unwrap()).join(item[\"receipt\"][\"stdout\"].as_str().unwrap())).unwrap()).unwrap();\nlet p=crate::catalog::integer(&data[\"curve\"][\"p\"]);\nmatch (p.bits()+63)/64 {\n");
    for limbs in 1..=10 {
        engines.push_str(&format!("{limbs} => n{limbs}::verify(out,write),\n"));
    }
    engines.push_str("_ => panic!(\"unsupported verification field width\") } }\n");
    fs::write(
        Path::new(&env::var("OUT_DIR").unwrap()).join("verification_engines.rs"),
        engines,
    )
    .unwrap();
    println!("cargo:rerun-if-env-changed=CRYPTO_BUILD_GIT_COMMIT");
    println!("cargo:rerun-if-env-changed=CRYPTO_BUILD_GIT_DIRTY");
    for reference in [
        "HEAD".to_string(),
        git(&root, &["symbolic-ref", "HEAD"]).unwrap_or_else(|| "HEAD".into()),
    ] {
        if let Some(path) = git(&root, &["rev-parse", "--git-path", &reference]) {
            println!("cargo:rerun-if-changed={path}");
        }
    }
    let revision = env::var("CRYPTO_BUILD_GIT_COMMIT")
        .ok()
        .or_else(|| git(&root, &["rev-parse", "HEAD"]))
        .unwrap_or_else(|| "unknown".into());
    let dirty = env::var("CRYPTO_BUILD_GIT_DIRTY").ok().unwrap_or_else(|| {
        git(
            &root,
            &["status", "--porcelain", "--untracked-files=normal"],
        )
        .map(|s| (!s.is_empty()).to_string())
        .unwrap_or_else(|| "unknown".into())
    });
    println!("cargo:rustc-env=ISOGENY_GIT_COMMIT={revision}");
    println!("cargo:rustc-env=ISOGENY_GIT_DIRTY={dirty}");
}
