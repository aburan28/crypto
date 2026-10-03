//! The native harness against the frozen outputs of the scripts it
//! replaces (IC_TOOL_PROGRAM.md §10a).  A port replaces a script only while
//! it reproduces that script's frozen outputs from the same frozen inputs.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn disclosed_sat_cli_writes_verified_evidence_once_and_preserves_it_on_rejection() {
    let dir = std::env::temp_dir().join(format!(
        "icprog-native-sat-cli-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&dir).unwrap();
    let report = dir.join("report.json");
    let run = || {
        Command::new(env!("CARGO_BIN_EXE_icprog"))
            .args(["sat-source-replay", "--root"])
            .arg(root())
            .arg("--out")
            .arg(&report)
            .output()
            .unwrap()
    };
    let first = run();
    assert!(
        first.status.success(),
        "{}",
        String::from_utf8_lossy(&first.stderr)
    );
    let bytes = std::fs::read(&report).unwrap();
    let value: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(value["status"], "PASS_NATIVE_RETAINED_SAT_SOURCE_REPLAY");
    assert_eq!(value["preparation_admission"]["rank"], 29);
    assert_eq!(value["native_solvers_executed"], 0);
    assert!(value["online_speedup"].is_null());
    let second = run();
    assert!(!second.status.success());
    assert!(String::from_utf8_lossy(&second.stderr).contains("create immutable"));
    assert_eq!(std::fs::read(&report).unwrap(), bytes);
    let rejected = dir.join("rejected.json");
    let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
        .args(["sat-source-replay", "--root"])
        .arg(&dir)
        .arg("--out")
        .arg(&rejected)
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(!rejected.exists());
    std::fs::remove_dir_all(dir).unwrap();
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

/// R05's pin was written by the declared `run.py pin`; with its record
/// removed, `icprog run r05 pin` must write the same bytes from the
/// candidate's outputs in the run tree (no binary runs: every output is
/// there).
#[test]
fn icprog_reproduces_r05s_pin_byte_for_byte() {
    let dir = unpack(&round("R05-presence-filter"), "r05-pin");
    let pin = dir.join("runs/pin/pin.json");
    let frozen = std::fs::read(&pin).unwrap();
    std::fs::remove_file(&pin).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
        .args(["run", "r05", "pin", "--root"])
        .arg(root())
        .arg("--runs")
        .arg(dir.join("runs"))
        .args(["--base", "/bin/true", "--cand", "/bin/true"])
        .args(["--isolate", "/bin/true"])
        .output()
        .expect("icprog runs");
    assert!(
        out.status.success(),
        "icprog failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let again = std::fs::read(&pin).unwrap();
    std::fs::remove_dir_all(&dir).ok();
    assert!(again == frozen, "the native pin differs from R05's");
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

/// R02's callgrind phase splits were written by the retired
/// `callgrind_phases.py` from R02's raw profiles, which R02's run tree
/// archives as `callgrind-parts.tar.xz`; `icprog callgrind-phases` must
/// write the same bytes.  It needs valgrind's `callgrind_annotate`, which
/// the harness workflow installs.
#[test]
#[ignore = "needs valgrind's callgrind_annotate; the ic programme harness workflow runs it"]
fn icprog_reproduces_r02s_callgrind_phases_byte_for_byte() {
    let dir = unpack(&round("R02-wide-tail-kernel"), "r02-callgrind");
    let cg = dir.join("runs/callgrind");
    let parts = dir.join("parts");
    std::fs::create_dir_all(&parts).unwrap();
    let tar = Command::new("tar")
        .arg("-xJf")
        .arg(cg.join("callgrind-parts.tar.xz"))
        .arg("-C")
        .arg(&parts)
        .status()
        .expect("tar runs");
    assert!(tar.success(), "tar -xJf callgrind-parts.tar.xz failed");
    let mut stems: Vec<String> = std::fs::read_dir(&parts)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".callgrind.out"))
        .collect();
    stems.sort();
    assert_eq!(stems.len(), 6, "R02 profiled two arms at three sizes");
    for stem in &stems {
        let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
            .arg("callgrind-phases")
            .arg(&parts)
            .arg(stem)
            .args(["--top", "40"])
            .output()
            .expect("icprog runs");
        assert!(
            out.status.success(),
            "icprog failed on {stem}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let tag = stem.trim_end_matches(".callgrind.out");
        let frozen = std::fs::read(cg.join(format!("{tag}.phases.json"))).unwrap();
        assert!(out.stdout == frozen, "{tag}: the phases differ from R02's");
    }
    std::fs::remove_dir_all(&dir).ok();
}

/// Every file under `dir`, as paths relative to it, sorted.
fn files(dir: &Path) -> Vec<PathBuf> {
    fn walk(base: &Path, dir: &Path, out: &mut Vec<PathBuf>) {
        for e in std::fs::read_dir(dir).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                walk(base, &p, out);
            } else {
                out.push(p.strip_prefix(base).unwrap().to_path_buf());
            }
        }
    }
    let mut out = Vec::new();
    walk(dir, dir, &mut out);
    out.sort();
    out
}

/// Ledger §23's claims, manifests and `analysis.json` were written by its
/// retired `claims.py` and `analyse.py` from its frozen run tree;
/// `icprog rule claims` and `icprog rule analyse` must write the same
/// bytes.  The one difference is by design: a replay record names its
/// checker, which is now `icprog`'s oracle and not `oracle.py`.
///
/// The claims hash the sources at the binary's build commit with `git
/// show`, so the commit's objects must be present; the harness workflow
/// fetches them.  The run writes into a scratch tree shaped like the
/// repository, so that every path a record holds is the frozen one.
#[test]
#[ignore = "needs §23's build commit 0bf67f16 in the object store; the ic programme harness workflow fetches it"]
fn icprog_reproduces_s23s_claims_and_analysis_byte_for_byte() {
    let s23 = root().join("research/ic_single_target_20260930");
    let scratch = std::env::temp_dir().join(format!("icprog-s23-{}", std::process::id()));
    let here = scratch.join("research/ic_single_target_20260930");
    for dir in [
        &here,
        &scratch.join("docs/ic"),
        &scratch.join("research/sat_factor_base_review_20260908/autolab"),
    ] {
        std::fs::create_dir_all(dir).unwrap();
    }
    std::os::unix::fs::symlink(s23.join("runs"), here.join("runs")).unwrap();
    std::os::unix::fs::symlink(
        root().join("research/ic_descent_20260930"),
        scratch.join("research/ic_descent_20260930"),
    )
    .unwrap();
    for f in [
        "research/ic_single_target_20260930/curve_ids.json",
        "research/ic_single_target_20260930/prediction.json",
        "docs/ic/boundary_targets.json",
        "research/sat_factor_base_review_20260908/autolab/protocol.json",
    ] {
        std::fs::copy(root().join(f), scratch.join(f)).unwrap();
    }
    let rule = |step: &str| {
        let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
            .args(["rule", step, "--root"])
            .arg(&scratch)
            .arg("--git")
            .arg(root())
            .output()
            .expect("icprog runs");
        assert!(
            out.status.success(),
            "icprog rule {step} failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        out.stdout
    };
    assert_eq!(
        String::from_utf8(rule("claims")).unwrap(),
        "408 of 408 rows pass the vs_rho check\n"
    );
    for sub in ["manifests", "claims"] {
        let (frozen, native) = (s23.join(sub), here.join(sub));
        let names = files(&frozen);
        assert_eq!(names, files(&native), "{sub}: the files differ from §23's");
        for name in &names {
            let mut got = std::fs::read_to_string(native.join(name)).unwrap();
            if name.to_string_lossy().ends_with(".replay.json") {
                got = got.replace(
                    "\"checker\": \"icprog's oracle (its own field arithmetic, not ic's)\"",
                    "\"checker\": \"oracle.py (Python field arithmetic)\"",
                );
            }
            let want = std::fs::read_to_string(frozen.join(name)).unwrap();
            assert!(got == want, "{sub}/{}: differs from §23's", name.display());
        }
    }
    let analysis = rule("analyse");
    std::fs::remove_dir_all(&scratch).ok();
    let frozen = std::fs::read(s23.join("analysis.json")).unwrap();
    assert!(
        analysis == frozen,
        "icprog's §23 analysis differs from the frozen analysis.json"
    );
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
