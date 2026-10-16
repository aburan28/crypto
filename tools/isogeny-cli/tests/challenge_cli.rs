//! End-to-end admission checks use real public curves and the independent verifier.
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicUsize, Ordering},
};
fn call(args: &[&str]) -> Output {
    let exe = std::env::var_os("ISOGENY_CHALLENGE_BINARY")
        .unwrap_or_else(|| env!("CARGO_BIN_EXE_isogeny").into());
    let output = Command::new(&exe)
        .args(args)
        .env("RAYON_NUM_THREADS", "1")
        .output()
        .unwrap();
    if let Some(root) = std::env::var_os("ISOGENY_CHALLENGE_EVIDENCE") {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let file = PathBuf::from(root).join(format!(
            "command-{:03}.json",
            NEXT.fetch_add(1, Ordering::SeqCst)
        ));
        assert!(!file.exists(), "evidence must be fresh");
        save(
            &file,
            &json!({"executable":exe,"args":args,"exit_code":output.status.code(),
            "stdout":String::from_utf8_lossy(&output.stdout),
            "stderr":String::from_utf8_lossy(&output.stderr)}),
        );
    }
    output
}
fn record(output: Output, code: i32) -> Value {
    assert_eq!(
        output.status.code(),
        Some(code),
        "stderr: {}\nstdout: {}",
        String::from_utf8_lossy(&output.stderr),
        String::from_utf8_lossy(&output.stdout)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}
fn path(p: &Path) -> &str {
    p.to_str().unwrap()
}
fn read(p: &Path) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn save(p: &Path, j: &Value) {
    fs::write(p, serde_json::to_vec_pretty(j).unwrap()).unwrap();
}
fn sha(b: &[u8]) -> String {
    Sha256::digest(b)
        .iter()
        .map(|v| format!("{v:02x}"))
        .collect()
}
fn clone_artifact(source: &Path, out: &Path) -> PathBuf {
    fs::create_dir(out).unwrap();
    fs::create_dir(out.join("p192")).unwrap();
    for file in ["search.json", "p192/ell-13.json", "p192/ell-13.stderr.txt"] {
        fs::copy(source.join(file), out.join(file)).unwrap();
    }
    out.to_owned()
}
fn alter_map(root: &Path, change: impl FnOnce(&mut Value)) {
    let file = root.join("p192/ell-13.json");
    let mut doc = read(&file);
    change(&mut doc);
    save(&file, &doc);
    let mut summary = read(&root.join("search.json"));
    summary["attempts"][0]["receipt"]["stdout_sha256"] = json!(sha(&fs::read(file).unwrap()));
    save(&root.join("search.json"), &summary);
}
#[test]
fn native_challenge_replays_maps_rejects_forgery_and_ranks_curves_separately() {
    let evidence_path = std::env::var_os("ISOGENY_CHALLENGE_EVIDENCE");
    let root = evidence_path
        .as_ref()
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            std::env::temp_dir().join(format!(
                "isogeny-challenge-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ))
        });
    fs::create_dir(&root).unwrap();
    let policy = record(call(&["challenge", "policy"]), 0);
    assert!(policy["curve_count"].as_u64().unwrap() > 20);
    let a = root.join("p192");
    let p521 = root.join("p521");
    let first = record(
        call(&[
            "challenge",
            "run",
            "--curve",
            "p192",
            "--ell",
            "13",
            "--timeout",
            "60",
            "--out",
            path(&a),
        ]),
        0,
    );
    let other = record(
        call(&[
            "challenge",
            "run",
            "--curve",
            "P-521",
            "--ell",
            "7",
            "--timeout",
            "60",
            "--out",
            path(&p521),
        ]),
        0,
    );
    assert_eq!(first["status"], "ACCEPTED");
    assert_eq!(first["score"], "13");
    assert_eq!(other["score"], "7");
    assert_ne!(first["source_icv1"], other["source_icv1"]);
    for result in [&first, &other] {
        assert!(
            result["verification"]["sampled_peak_resident_bytes"]
                .as_str()
                .unwrap()
                .parse::<u64>()
                .unwrap()
                > 0
        );
    }
    let original = a.join("construction");
    let before = fs::read(original.join("search.json")).unwrap();
    save(
        &original.join("replay.json"),
        &json!({"status":"PASS","score":"999983"}),
    );
    let clean = root.join("fresh-judge");
    let checked = record(
        call(&[
            "challenge",
            "judge",
            "--input",
            path(&original),
            "--out",
            path(&clean),
        ]),
        0,
    );
    assert_eq!(checked["score"], "13");
    assert_eq!(fs::read(original.join("search.json")).unwrap(), before);
    let retained = fs::read(clean.join("entry.json")).unwrap();
    record(
        call(&[
            "challenge",
            "judge",
            "--input",
            path(&original),
            "--out",
            path(&clean),
        ]),
        2,
    );
    assert_eq!(fs::read(clean.join("entry.json")).unwrap(), retained);
    for kind in [
        "composite",
        "inflated-degree",
        "source",
        "kernel",
        "map",
        "compact-only",
        "traversal",
        "missing-hash",
        "duplicate-kernel",
    ] {
        let input = clone_artifact(&original, &root.join(format!("bad-{kind}")));
        match kind {
            "composite" => {
                let mut j = read(&input.join("search.json"));
                j["attempts"][0]["ell"] = json!(15);
                save(&input.join("search.json"), &j);
            }
            "inflated-degree" => alter_map(&input, |j| {
                j["result"]["ell"] = json!(999983);
            }),
            "source" => alter_map(&input, |j| {
                j["curve"]["a"] = json!("1");
            }),
            "kernel" => alter_map(&input, |j| {
                let c = &mut j["result"]["isogenies"][0]["kernel"][0];
                *c = if *c == "1" { json!("2") } else { json!("1") };
            }),
            "map" => alter_map(&input, |j| {
                let c = &mut j["result"]["isogenies"][0]["x_map"]["num"][0];
                *c = if *c == "1" { json!("2") } else { json!("1") };
            }),
            "compact-only" => alter_map(&input, |j| {
                j["result"]["isogenies"][0]["x_map"] = json!({"ideal":"claimed compact map"});
            }),
            "duplicate-kernel" => alter_map(&input, |j| {
                j["result"]["isogenies"][1] = j["result"]["isogenies"][0].clone();
            }),
            "traversal" => {
                let mut j = read(&input.join("search.json"));
                j["attempts"][0]["curve"] = json!("../../p192");
                save(&input.join("search.json"), &j);
            }
            "missing-hash" => {
                let mut j = read(&input.join("search.json"));
                j["attempts"][0]["receipt"]
                    .as_object_mut()
                    .unwrap()
                    .remove("stdout_sha256");
                save(&input.join("search.json"), &j);
            }
            _ => unreachable!(),
        }
        let out = root.join(format!("judged-{kind}"));
        let rejected = record(
            call(&[
                "challenge",
                "judge",
                "--input",
                path(&input),
                "--out",
                path(&out),
            ]),
            1,
        );
        assert_eq!(rejected["status"], "REJECTED", "{kind}");
        assert!(rejected["score"].is_null(), "{kind}");
        assert_eq!(read(&out.join("entry.json"))["status"], "REJECTED");
    }
    let fake = root.join("forged-score");
    fs::create_dir(&fake).unwrap();
    save(
        &fake.join("entry.json"),
        &json!({"status":"ACCEPTED","source_icv1":first["source_icv1"],"prime_degree":"999983"}),
    );
    let list = root.join("submissions.json");
    save(&list, &json!([original, p521.join("construction"), fake]));
    let board = root.join("board");
    let j = record(
        call(&[
            "challenge",
            "leaderboard",
            "--input",
            path(&list),
            "--out",
            path(&board),
        ]),
        0,
    );
    let rows = j["curves"].as_array().unwrap();
    assert_eq!(rows.len() as u64, policy["curve_count"].as_u64().unwrap());
    let winners: Vec<_> = rows.iter().filter(|r| r["status"] == "VERIFIED").collect();
    assert_eq!(winners.len(), 2);
    assert_eq!(
        winners
            .iter()
            .find(|r| r["source_icv1"] == first["source_icv1"])
            .unwrap()["best_prime_degree"],
        "13"
    );
    assert_eq!(
        winners
            .iter()
            .find(|r| r["source_icv1"] == other["source_icv1"])
            .unwrap()["best_prime_degree"],
        "7"
    );
    assert_eq!(j["attempts"][2]["status"], "REJECTED");
    let p = record(
        call(&[
            "challenge",
            "plan",
            "--curve",
            "p192",
            "--from",
            "10454",
            "--to",
            "12000",
        ]),
        0,
    );
    for c in p["curves"][0]["candidates"].as_array().unwrap() {
        assert_eq!(c["proof_status"], "CANDIDATE_ONLY");
        assert!(c["ell"].as_str().unwrap().parse::<u64>().unwrap() > 10453);
    }
    // Leave failed-attempt evidence intact if any assertion above fails.
    if evidence_path.is_none() {
        fs::remove_dir_all(root).unwrap();
    }
}
