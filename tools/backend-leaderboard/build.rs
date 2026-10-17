use std::path::PathBuf;

fn main() {
    // One counter-to-unit vocabulary: compile the harness's literal declaration
    // without depending on its solver/curve implementation or root build.
    let root = PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR").unwrap()).join("../..");
    let path = root.join("src/cryptanalysis/ecbench/methods.rs");
    println!("cargo:rerun-if-changed={}", path.display());
    let text = std::fs::read_to_string(path).expect("read native pricing vocabulary");
    let start = text
        .find("pub const PRICED_NATIVE:")
        .expect("native pricing declaration");
    let end = start
        + text[start..]
            .find("];")
            .expect("end native pricing declaration")
        + 2;
    let out = PathBuf::from(std::env::var_os("OUT_DIR").unwrap()).join("priced_native.rs");
    std::fs::write(out, &text[start..end]).expect("write native pricing vocabulary");
}
