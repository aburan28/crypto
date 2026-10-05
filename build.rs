//! Build-time facts the crate's `cfg`s depend on.
//!
//! `has_f6_ic` is set when the `koblitz_index_calculus.rs` being compiled
//! defines PR #1333's `groebner_decompose_f6_ic`. The ecbench `pdp3-koblitz`
//! oracle (`src/cryptanalysis/ic_framework/pdp3_koblitz.rs`) calls it, and
//! the frozen-source replay workflows rebuild this crate with pre-F6
//! snapshots of that file
//! (`research/notes/ecc2k130/compact_frozen_source_replay_20260929`).
//! Keying the module on the function's presence keeps those builds
//! compiling without editing `Cargo.toml`, which other frozen evaluations
//! pin by hash. A sparse checkout without this script leaves the module out,
//! which is what those builds need.

fn main() {
    const SOURCE: &str = "src/cryptanalysis/koblitz_index_calculus.rs";
    println!("cargo::rustc-check-cfg=cfg(has_f6_ic)");
    println!("cargo::rerun-if-changed=build.rs");
    println!("cargo::rerun-if-changed={SOURCE}");
    let defines_f6 = std::fs::read_to_string(SOURCE)
        .map(|s| s.contains("pub fn groebner_decompose_f6_ic("))
        .unwrap_or(false);
    if defines_f6 {
        println!("cargo::rustc-cfg=has_f6_ic");
    }
}
