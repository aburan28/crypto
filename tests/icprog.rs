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
    // This is only a preflight placeholder: frozen outputs prevent any binary run.
    let true_bin = if cfg!(target_os = "macos") {
        "/usr/bin/true"
    } else {
        "/bin/true"
    };
    let dir = unpack(&round("R05-presence-filter"), "r05-pin");
    let pin = dir.join("runs/pin/pin.json");
    let frozen = std::fs::read(&pin).unwrap();
    std::fs::remove_file(&pin).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
        .args(["run", "r05", "pin", "--root"])
        .arg(root())
        .arg("--runs")
        .arg(dir.join("runs"))
        .args(["--base", true_bin, "--cand", true_bin])
        .args(["--isolate", true_bin])
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

fn holdouts(round_dir: &Path, sizes: &[&str], seeds: &str, first: &str, check: bool) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_icprog"));
    cmd.arg("holdouts").arg(round_dir).arg("--root").arg(root());
    for s in sizes {
        cmd.args(["--size", s]);
    }
    cmd.args(["--seeds", seeds, "--first-target", first]);
    if check {
        cmd.arg("--check");
    }
    cmd.output().expect("icprog runs")
}

/// R02b's and R05's holdouts were drawn by suite v1's `make_suite.py`
/// construction; `icprog holdouts --check` must re-derive every file and
/// each `SHA256SUMS` byte for byte, as it must every suite v1 row.
#[test]
fn icprog_reproduces_the_frozen_holdouts_and_suite_rows_byte_for_byte() {
    for (dir, sizes, seeds, first) in [
        (
            "R02b-wide-tail-retest",
            &["1,59", "0,61"][..],
            "206,207,208,209",
            "103",
        ),
        (
            "R05-presence-filter",
            &["0,53", "1,59", "0,61"][..],
            "210,211,212,213",
            "111",
        ),
    ] {
        let out = holdouts(&round(dir), sizes, seeds, first, true);
        assert!(
            out.status.success(),
            "{dir}: {}{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
    // Suite v1's own rows: seeds 201–204 from target 1, file for file.
    let tmp = std::env::temp_dir().join(format!("icprog-holdouts-{}", std::process::id()));
    let suite = root().join("research/ic_tool_program/suite/v1/params/S");
    for (a, n) in [
        (1, 19),
        (1, 23),
        (1, 45),
        (0, 37),
        (1, 43),
        (1, 47),
        (0, 57),
        (0, 41),
        (0, 53),
        (1, 59),
        (0, 61),
    ] {
        let dir = tmp.join(format!("k{a}n{n}"));
        let size = format!("{a},{n}");
        let out = holdouts(&dir, &[&size], "201,202,203,204", "1", false);
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let drawn = dir.join("holdouts");
        let slug = std::fs::read_dir(&drawn)
            .unwrap()
            .filter_map(|e| e.ok())
            .find(|e| e.path().is_dir())
            .expect("a slug directory")
            .path();
        for i in 1..=8u32 {
            let m = i.div_ceil(2);
            let ours = std::fs::read(slug.join(format!("M{m}-T{i}.json"))).unwrap();
            let frozen = std::fs::read(suite.join(format!("k{a}n{n}/M{m}-T{i:02}.json"))).unwrap();
            assert!(
                ours == frozen,
                "k{a}n{n} M{m}-T{i:02} differs from suite v1's"
            );
        }
        // Drawn once: a second draw into the same round refuses.
        let again = holdouts(&dir, &[&size], "201,202,203,204", "1", false);
        assert!(!again.status.success());
    }
    std::fs::remove_dir_all(&tmp).ok();
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

/// A stub `ic`: each case's first argument says what it does.
const STUB_IC: &str = r#"#!/bin/sh
case "$1" in
  report)
    printf '{"status": "%s", "build_commit": "abc123", "detail": {"n": 3, "items": [{"code": "x", "k": 1}, {"code": "y"}]}}\n' "$2" > "$3"
    exit 0 ;;
  refuse) echo "ic: refused: $2" >&2; exit 2 ;;
  panic) echo "thread 'main' panicked" >&2; exit 101 ;;
  hang) exec sleep 30 ;;
  cat) cat "$2" >&2; exit 0 ;;
esac
exit 3
"#;

/// The conformance runner (N4) on a stub `ic`, every rule of the scripts
/// it replaces (`conformance/run.py`, `v2/run.py`) exercised: the exit
/// rule, under which a panic never passes; the timeout; `stderr_contains`;
/// the report's `json_equals` (an absent key is `None`),
/// `json_equals_build_commit`, `json_paths`, `json_contains` and
/// `same_outputs_as`; a file copied with a key set, from `{here}`; and the
/// `until` and `supersedes` rules that choose the cases.
#[test]
fn icprog_conformance_runs_each_rule_on_a_stub_ic() {
    let dir = std::env::temp_dir().join(format!(
        "icprog-conformance-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let conf = dir.join("research/ic_tool_program/conformance");
    std::fs::create_dir_all(conf.join("v1")).unwrap();
    std::fs::create_dir_all(conf.join("v2/params")).unwrap();
    let report = "{tmp}/r.json";
    let v1 = serde_json::json!({"cases": [
        {"id": "C001", "argv": ["report", "complete", report], "timeout_s": 10,
         "expect": {"exit": "zero", "json_file": report,
                    "json_equals": {"status": "complete", "absent": null},
                    "json_equals_build_commit": true,
                    "json_paths": {"detail.n": 3.0, "detail.items.1.code": "y"},
                    "json_contains": {"detail.items": [{"code": "x"}]},
                    "same_outputs_as": {"argv": ["report", "complete", "{tmp}/s.json"],
                                        "json_file": "{tmp}/s.json",
                                        "paths": ["status", "detail.n"]}}},
        {"id": "C002", "argv": ["refuse", "bad-field"], "timeout_s": 10,
         "expect": {"exit": "nonzero", "stderr_contains": ["refused: bad-field"]}},
        {"id": "C003", "argv": ["panic"], "timeout_s": 10, "expect": {"exit": "nonzero"}},
        {"id": "C004", "argv": ["hang"], "timeout_s": 0.5, "expect": {"exit": "zero"}},
        {"id": "C005", "until": "B1", "argv": ["refuse", "x"], "timeout_s": 10,
         "expect": {"exit": 2}},
        {"id": "C006", "argv": ["report", "refused", report], "timeout_s": 10,
         "expect": {"exit": "zero", "json_file": report,
                    "json_equals": {"status": "complete"}}}
    ]});
    let v2 = serde_json::json!({"cases": [
        {"id": "C010", "step": "B1", "timeout_s": 10,
         "files": {"doc.json": {"copy": "{here}/base.json", "set": {"a.b": 5}}},
         "argv": ["cat", "{tmp}/doc.json"],
         "expect": {"exit": "zero", "stderr_contains": ["\"b\": 5", "\"c\": 2"]}},
        {"id": "C011", "step": "B1", "supersedes": "C002", "argv": ["refuse", "other"],
         "timeout_s": 10,
         "expect": {"exit": "nonzero", "stderr_contains": ["refused: other"]}}
    ]});
    std::fs::write(conf.join("v1/cases.json"), v1.to_string()).unwrap();
    std::fs::write(conf.join("v2/cases.json"), v2.to_string()).unwrap();
    std::fs::write(
        conf.join("v2/params/base.json"),
        r#"{"a": {"b": 1, "c": 2}}"#,
    )
    .unwrap();
    let ic = dir.join("ic");
    std::fs::write(&ic, STUB_IC).unwrap();
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&ic, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    let run = |steps: &str, commit: &str| -> (bool, serde_json::Value) {
        let out_file = dir.join(format!("report-{steps}.json"));
        let out = Command::new(env!("CARGO_BIN_EXE_icprog"))
            .arg("conformance")
            .arg("--ic")
            .arg(&ic)
            .args(["--steps", steps, "--build-commit", commit, "--root"])
            .arg(&dir)
            .arg("--out")
            .arg(&out_file)
            .output()
            .unwrap();
        let printed: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
        let written: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&out_file).unwrap()).unwrap();
        assert_eq!(printed, written);
        (out.status.success(), written)
    };
    let verdicts = |doc: &serde_json::Value| -> Vec<(String, bool, Vec<String>)> {
        doc["results"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| {
                let why = r["why"].as_array().unwrap();
                (
                    r["id"].as_str().unwrap().to_string(),
                    r["pass"].as_bool().unwrap(),
                    why.iter()
                        .map(|w| w.as_str().unwrap().to_string())
                        .collect(),
                )
            })
            .collect()
    };
    let case = |id: &str, pass: bool, why: &[&str]| {
        (
            id.to_string(),
            pass,
            why.iter().map(|w| w.to_string()).collect::<Vec<_>>(),
        )
    };

    let (all, b0) = run("B0", "abc123");
    assert!(!all);
    assert_eq!(b0["steps"], serde_json::json!(["B0"]));
    assert_eq!(
        (b0["cases"].as_u64(), b0["passed"].as_u64()),
        (Some(6), Some(3))
    );
    assert_eq!(
        verdicts(&b0),
        vec![
            case("C001", true, &[]),
            case("C002", true, &[]),
            case("C003", false, &["panicked (exit status 101)"]),
            case("C004", false, &["no exit within 0.5 s"]),
            case("C005", true, &[]),
            case(
                "C006",
                false,
                &[r#"status = "refused", expected "complete""#]
            ),
        ]
    );
    assert_eq!(b0["results"][2]["exit"], 101);
    assert_eq!(
        b0["results"][1]["stderr_tail"],
        serde_json::json!(["ic: refused: bad-field"])
    );

    // B1's cases join; C005 ends at B1 and C011 supersedes C002.  The
    // steps are reported in the scripts' order, whatever order they came in.
    let (_, b1) = run("B1,B0", "zzz");
    assert_eq!(b1["steps"], serde_json::json!(["B0", "B1"]));
    assert_eq!(
        verdicts(&b1),
        vec![
            case(
                "C001",
                false,
                &[r#"build_commit = "abc123", expected "zzz""#]
            ),
            case("C003", false, &["panicked (exit status 101)"]),
            case("C004", false, &["no exit within 0.5 s"]),
            case(
                "C006",
                false,
                &[r#"status = "refused", expected "complete""#]
            ),
            case("C010", true, &[]),
            case("C011", true, &[]),
        ]
    );
    let unknown = Command::new(env!("CARGO_BIN_EXE_icprog"))
        .arg("conformance")
        .arg("--ic")
        .arg(&ic)
        .args(["--steps", "B9", "--root"])
        .arg(&dir)
        .output()
        .unwrap();
    assert!(!unknown.status.success());
    assert!(String::from_utf8_lossy(&unknown.stderr).contains("unknown step \"B9\""));
    // A set whose files no longer match its `SHA256SUMS` runs nothing.
    std::fs::write(
        conf.join("v2/SHA256SUMS"),
        format!("{}  params/base.json\n", "0".repeat(64)),
    )
    .unwrap();
    let tampered = Command::new(env!("CARGO_BIN_EXE_icprog"))
        .arg("conformance")
        .arg("--ic")
        .arg(&ic)
        .args(["--steps", "B0", "--root"])
        .arg(&dir)
        .output()
        .unwrap();
    assert!(!tampered.status.success());
    assert!(String::from_utf8_lossy(&tampered.stderr).contains("does not match SHA256SUMS"));
    assert!(tampered.stdout.is_empty());
    std::fs::remove_dir_all(&dir).unwrap();
}
