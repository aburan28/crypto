//! Build the canonical cover catalogue with unchanged native checker modules.
//! The narrow crate avoids unrelated root-library compilation failures.
#![allow(dead_code)]
extern crate self as crypto_lib;
#[path = "../../src/binary_ecc/f2m.rs"]
pub mod binary_ecc;
#[path = "../../src/ct_bignum.rs"]
pub mod ct_bignum;
#[path = "../../src/cryptanalysis/ecc2k130_guard.rs"]
pub mod ecc2k130_guard;
pub mod cryptanalysis {
    pub use crate::ecc2k130_guard;
}
#[path = "../../src/hash/sha256.rs"]
pub mod hash_sha256;
#[path = "../../src/utils/mod.rs"]
pub mod utils;
pub mod hash {
    pub use crate::hash_sha256 as sha256;
}
#[path = "../../src/bin/curve_cover_check/checker.rs"]
mod checker;
#[path = "../../src/bin/curve_cover_check/links.rs"]
mod links;
#[path = "../../src/bin/curve_cover_check/models.rs"]
mod models;

fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 && !(args.len() == 4 && args[3] == "--check") {
        return Err(
            "usage: catalogue-covers REGISTRY_JSON COVERS_JSON COVER_LINKS_YAML [--check]".into(),
        );
    }
    let bytes = std::fs::read(&args[0]).map_err(|e| e.to_string())?;
    let report = checker::catalog(&bytes)?;
    let registry = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    let graph = links::graph(&registry, &report)?;
    for (path, value) in [(&args[1], &report), (&args[2], &graph)] {
        let rendered = serde_json::to_string_pretty(value).map_err(|e| e.to_string())? + "\n";
        if args.len() == 4 {
            if std::fs::read_to_string(path).map_err(|e| e.to_string())? != rendered {
                return Err(format!("stale or altered catalogue output: {path}"));
            }
        } else {
            std::fs::write(path, rendered).map_err(|e| e.to_string())?;
        }
    }
    println!("{}", report["summary"]);
    if report["summary"]["invalid_input"].as_u64().unwrap_or(0) > 0 {
        return Err("invalid catalogue models; evidence retained".into());
    }
    Ok(())
}
