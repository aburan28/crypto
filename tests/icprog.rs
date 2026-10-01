//! The native harness against the frozen outputs of the scripts it
//! replaces (IC_TOOL_PROGRAM.md §10a).  A port replaces a script only while
//! it reproduces that script's frozen outputs from the same frozen inputs.

use std::path::Path;
use std::process::Command;

/// R03's `analysis.json` was written by the retired `analyse.py` from
/// R03's `runs.tar.xz`; `icprog analyse r03` must write the same bytes.
#[test]
fn icprog_reproduces_r03s_frozen_analysis_byte_for_byte() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let round = root.join("research/ic_tool_program/rounds/R03-curve-construction");
    let dir = std::env::temp_dir().join(format!("icprog-r03-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let tar = Command::new("tar")
        .arg("-xJf")
        .arg(round.join("runs.tar.xz"))
        .arg("-C")
        .arg(&dir)
        .status()
        .expect("tar runs");
    assert!(tar.success(), "tar -xJf runs.tar.xz failed");
    let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
        .args(["analyse", "r03", "--root"])
        .arg(root)
        .arg("--runs")
        .arg(dir.join("runs"))
        .output()
        .expect("icprog runs");
    std::fs::remove_dir_all(&dir).ok();
    assert!(
        out.status.success(),
        "icprog failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let frozen = std::fs::read(round.join("analysis.json")).unwrap();
    assert!(
        out.stdout == frozen,
        "icprog's R03 analysis differs from the frozen analysis.json"
    );
}

/// A run tree that is not one is refused, with the reason.
#[test]
fn icprog_refuses_a_missing_run_tree() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
        .args(["analyse", "r05", "--root"])
        .arg(root)
        .args(["--runs", "/nonexistent/run/tree"])
        .output()
        .expect("icprog runs");
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("is not a run tree"));
}
