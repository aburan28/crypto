//! The native harness against the frozen outputs of the scripts it
//! replaces (IC_TOOL_PROGRAM.md §10a).  A port replaces a script only while
//! it reproduces that script's frozen outputs from the same frozen inputs.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

fn round(dir: &str) -> PathBuf {
    root().join("research/ic_tool_program/rounds").join(dir)
}

/// A round's `runs.tar.xz`, unpacked into a directory of the test's own.
fn unpack(round: &Path, tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("icprog-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let tar = Command::new("tar")
        .arg("-xJf")
        .arg(round.join("runs.tar.xz"))
        .arg("-C")
        .arg(&dir)
        .status()
        .expect("tar runs");
    assert!(tar.success(), "tar -xJf runs.tar.xz failed");
    dir
}

fn analyse(name: &str, runs: &Path) -> Output {
    let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
        .args(["analyse", name, "--root"])
        .arg(root())
        .arg("--runs")
        .arg(runs)
        .output()
        .expect("icprog runs");
    assert!(
        out.status.success(),
        "icprog failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

/// R03's `analysis.json` was written by the retired `analyse.py` from
/// R03's `runs.tar.xz`; `icprog analyse r03` must write the same bytes.
#[test]
fn icprog_reproduces_r03s_frozen_analysis_byte_for_byte() {
    let round = round("R03-curve-construction");
    let dir = unpack(&round, "r03");
    let out = analyse("r03", &dir.join("runs"));
    std::fs::remove_dir_all(&dir).ok();
    let frozen = std::fs::read(round.join("analysis.json")).unwrap();
    assert!(
        out.stdout == frozen,
        "icprog's R03 analysis differs from the frozen analysis.json"
    );
}

/// R05's `analysis.json` is `icprog analyse r05`'s output from R05's
/// `runs.tar.xz`, and stays so.
#[test]
fn icprog_reproduces_r05s_committed_analysis_byte_for_byte() {
    let round = round("R05-presence-filter");
    let dir = unpack(&round, "r05");
    let out = analyse("r05", &dir.join("runs"));
    std::fs::remove_dir_all(&dir).ok();
    let committed = std::fs::read(round.join("analysis.json")).unwrap();
    assert!(
        out.stdout == committed,
        "icprog's R05 analysis differs from the committed analysis.json"
    );
}

/// The extension rule is re-tested from the runs: a record of extended
/// sizes that the rule's test does not give is a reason to reject.
#[test]
fn icprog_flags_an_extension_record_its_test_does_not_give() {
    let dir = unpack(&round("R05-presence-filter"), "r05-ext");
    let runs = dir.join("runs");
    std::fs::write(
        runs.join("extended.json"),
        "{\n \"suite\": [\n  \"icv1-f2m53-tm56619371-dac20a85\"\n ]\n}\n",
    )
    .unwrap();
    let out = analyse("r05", &runs);
    std::fs::remove_dir_all(&dir).ok();
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.contains("the extension record differs from the rule's test"));
    assert!(text.contains("\"accepted\": false"));
}

/// A run tree that is not one is refused, with the reason.
#[test]
fn icprog_refuses_a_missing_run_tree() {
    let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
        .args(["analyse", "r05", "--root"])
        .arg(root())
        .args(["--runs", "/nonexistent/run/tree"])
        .output()
        .expect("icprog runs");
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("is not a run tree"));
}
