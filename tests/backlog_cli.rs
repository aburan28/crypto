//! End to end through the `backlog` binary: the committed seed cards validate,
//! `new` scaffolds a card that validates as a draft, `cairn-post` builds an
//! objective draft from the card's funding (and refuses a card whose checker
//! is blocked), and — when `cairn` is on PATH — a lease race between two lab
//! replicas resolves to one holder, and publish → claim → exec → release →
//! status round-trips. Everything runs in temp directories; nothing touches a
//! real lab, node or remote.

#![cfg(unix)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_backlog")
}

fn repo() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("backlog-it-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn run(cmd: &mut Command) -> Output {
    let out = cmd.output().expect("spawn");
    if !out.status.success() {
        eprintln!(
            "{:?} exited {:?}\nstdout:\n{}\nstderr:\n{}",
            cmd,
            out.status.code(),
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
    out
}

fn stdout(o: &Output) -> String {
    String::from_utf8_lossy(&o.stdout).into_owned()
}

/// A scratch repository root with the schema, for `new` and fixtures.
fn scratch_root(name: &str) -> PathBuf {
    let root = tmp(name);
    fs::create_dir_all(root.join("backlog/schema")).unwrap();
    fs::copy(
        repo().join("backlog/schema/card.v1.json"),
        root.join("backlog/schema/card.v1.json"),
    )
    .unwrap();
    fs::create_dir_all(root.join("problems")).unwrap();
    root
}

#[test]
fn the_committed_cards_validate() {
    let out = run(Command::new(bin())
        .arg("--root")
        .arg(repo())
        .args(["validate", "--all", "--no-git"]));
    assert!(out.status.success(), "{}", stdout(&out));
    let text = stdout(&out);
    for id in [
        "XP-20261015-p192-isogeny-aux-search-t48",
        "XP-20261015-pkm-m3-n21-m4-n20",
        "XP-20261015-shape-g1-frozen-cells",
        "XP-20261015-kummer-p224-depth",
    ] {
        assert!(text.contains(&format!("{id}: ok")), "{text}");
    }
    // `list` orders by funding: lead 1 first, the P-192 search last.
    let list = stdout(&run(Command::new(bin())
        .arg("--root")
        .arg(repo())
        .arg("list")));
    let lines: Vec<&str> = list.lines().collect();
    assert!(
        lines
            .first()
            .unwrap()
            .starts_with("XP-20261015-pkm-m3-n21-m4-n20"),
        "{list}"
    );
    assert!(
        lines
            .last()
            .unwrap()
            .starts_with("XP-20261015-p192-isogeny-aux-search-t48"),
        "{list}"
    );
}

#[test]
fn new_scaffolds_a_draft_that_validates_and_refuses_to_leave_draft_with_todos() {
    let root = scratch_root("new");
    let out = run(Command::new(bin()).arg("--root").arg(&root).args([
        "new",
        "demo-search",
        "--date",
        "2026-10-15",
        "--title",
        "A demo",
    ]));
    assert!(out.status.success());
    let dir = root.join("backlog/XP-20261015-demo-search");
    assert!(dir.join("card.json").is_file() && dir.join("README.md").is_file());
    let v = run(Command::new(bin())
        .arg("--root")
        .arg(&root)
        .args(["validate", "XP-20261015-demo-search"]));
    assert!(v.status.success(), "{}", stdout(&v));
    assert!(stdout(&v).contains("warning:") && stdout(&v).contains("still TODO"));
    // The same card marked blocked fails: TODOs and an empty blocked_on.
    let path = dir.join("card.json");
    let text = fs::read_to_string(&path)
        .unwrap()
        .replace("\"status\": \"draft\"", "\"status\": \"blocked\"");
    fs::write(&path, text).unwrap();
    let v = run(Command::new(bin())
        .arg("--root")
        .arg(&root)
        .args(["validate", "--all"]));
    assert_eq!(v.status.code(), Some(1));
    assert!(
        stdout(&v).contains("names what blocks it"),
        "{}",
        stdout(&v)
    );
    // A second `new` with the same id is refused.
    let again = Command::new(bin())
        .arg("--root")
        .arg(&root)
        .args(["new", "demo-search", "--date", "2026-10-15"])
        .output()
        .unwrap();
    assert!(!again.status.success());
}

/// A minimal objective-mode card whose checker is available.
fn objective_card(root: &Path, id: &str) {
    let dir = root.join("backlog").join(id);
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("README.md"), "# demo\n").unwrap();
    let card = serde_json::json!({
        "schema": "backlog.card/v1", "id": id, "title": "Demo search", "created": "2026-10-15",
        "authors": ["test"], "status": "draft", "blocked_on": [], "links": {},
        "question": "Is there a witness?",
        "decision_value": "A witness found here would change which construction the next lane builds, and a bounded failure would be recorded with its exact seed space so nobody repeats it blindly.",
        "protocol": {"path": null, "blob": null},
        "code": {"repo": "crypto", "commit": null, "branch": "x", "build": []},
        "resources": {"cpus": 1, "ram_gb": 1, "disk_gb": 1, "gpu": null, "est_core_hours": 1, "est_wall_hours": null},
        "run": {"argv": ["solver", "{seed_lo}", "{seed_hi}"], "env": {}, "timeout_per_shard_s": 60},
        "shards": {"count": 8, "scheme": "seed_range", "params": {"seed_lo": "0", "range_size": "16"}, "stop": "first_success"},
        "checkpoint": {"file": null, "resume_flag": null},
        "outputs": {"result_schema": "x", "verify": ["check"], "land": "research/x"},
        "cairn": {"mode": "objective", "kind": "workspace", "goal": "GOAL-demo",
                  "checker": {"status": "available", "bin": "check", "argv": ["check"], "artifact_schema": {"type": "object"}},
                  "verifier": {"kind": "workspace", "base": "base.manifest.json", "base_sha256": "0000000000000000000000000000000000000000000000000000000000000000",
                               "editable_paths": ["claim"], "score_command": ["cargo", "run", "--offline", "-q"], "score_path": "score.json"},
                  "piecework": {"units": 8, "unit": "16 seeds"}, "stop": "first_success",
                  "funding": {"units": 800, "rationale": "test"}}
    });
    fs::write(
        dir.join("card.json"),
        serde_json::to_string_pretty(&card).unwrap(),
    )
    .unwrap();
}

#[test]
fn cairn_post_is_a_dry_run_by_default_and_refuses_a_blocked_checker() {
    let root = scratch_root("post");
    objective_card(&root, "XP-20261015-demo-post");
    let out = run(Command::new(bin())
        .arg("--root")
        .arg(&root)
        .args(["cairn-post", "XP-20261015-demo-post"]));
    assert!(out.status.success());
    let draft: serde_json::Value = serde_json::from_slice(&out.stdout).expect("the draft is JSON");
    assert_eq!(draft["reward"].as_u64(), Some(800));
    assert_eq!(draft["piecework"]["unit_price"].as_u64(), Some(100));
    assert_eq!(draft["goal"], "GOAL-demo");
    assert_eq!(draft["verifier"]["kind"], "workspace");
    assert!(String::from_utf8_lossy(&out.stderr).contains("nothing posted"));
    // --confirm without a node is refused before anything is signed.
    let c = Command::new(bin())
        .arg("--root")
        .arg(&root)
        .args(["cairn-post", "XP-20261015-demo-post", "--confirm"])
        .output()
        .unwrap();
    assert!(!c.status.success());
    // The committed P-192 card's checker is blocked: nothing to post.
    let p = Command::new(bin())
        .arg("--root")
        .arg(repo())
        .args(["cairn-post", "XP-20261015-p192-isogeny-aux-search-t48"])
        .output()
        .unwrap();
    assert!(!p.status.success());
    assert!(String::from_utf8_lossy(&p.stderr).contains("cairn-certificate-checker-python-only"));
    // A lab-task card is not an objective.
    let l = Command::new(bin())
        .arg("--root")
        .arg(repo())
        .args(["cairn-post", "XP-20261015-pkm-m3-n21-m4-n20"])
        .output()
        .unwrap();
    assert!(!l.status.success());
}

// ── with cairn ──────────────────────────────────────────────────────────

fn have_cairn() -> bool {
    let ok = Command::new("cairn")
        .args(["lab", "help"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if !ok {
        eprintln!("SKIPPED: `cairn` is not on PATH; the lab lease and publish tests need it (docs/CONDUCTOR.md: the cloud session hook installs it)");
    }
    ok
}

fn cairn(args: &[&str]) -> Output {
    run(Command::new("cairn").args(args))
}

fn public_key(identity: &Path) -> String {
    let v: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(identity).unwrap()).unwrap();
    v["public"].as_str().unwrap().to_string()
}

/// A lab at `w/A` with two writers, and its clone at `w/B`.
fn two_replicas(w: &Path) -> (PathBuf, PathBuf, PathBuf, PathBuf) {
    let (id1, id2) = (w.join("one.identity.json"), w.join("two.identity.json"));
    let (a, b) = (w.join("A"), w.join("B"));
    for id in [&id1, &id2] {
        assert!(cairn(&["lab", "identity", "--out", id.to_str().unwrap()])
            .status
            .success());
    }
    assert!(cairn(&[
        "lab",
        "--lab",
        a.to_str().unwrap(),
        "init",
        "--name",
        "backlog-test",
        "--identity",
        id1.to_str().unwrap()
    ])
    .status
    .success());
    let status: serde_json::Value = serde_json::from_slice(
        &cairn(&["lab", "--lab", a.to_str().unwrap(), "status", "--json"]).stdout,
    )
    .unwrap();
    let space = status["space"].as_str().unwrap().to_string();
    assert!(cairn(&[
        "lab",
        "--lab",
        a.to_str().unwrap(),
        "admit",
        &public_key(&id2),
        "--roles",
        "writer",
        "--identity",
        id1.to_str().unwrap()
    ])
    .status
    .success());
    assert!(cairn(&[
        "lab",
        "--lab",
        b.to_str().unwrap(),
        "clone",
        "--space",
        &space,
        "--dir",
        a.to_str().unwrap()
    ])
    .status
    .success());
    (a, b, id1, id2)
}

#[test]
fn concurrent_claims_resolve_to_one_holder_after_sync() {
    if !have_cairn() {
        return;
    }
    let w = tmp("race");
    let (a, b, id1, id2) = two_replicas(&w);
    let task = "XP-20261015-pkm-m3-n21-m4-n20/shard-0";
    // Both claim offline, neither having seen the other.
    assert!(cairn(&[
        "lab",
        "--lab",
        a.to_str().unwrap(),
        "claim",
        task,
        "--holder",
        "box-a",
        "--ttl",
        "600",
        "--identity",
        id1.to_str().unwrap()
    ])
    .status
    .success());
    assert!(cairn(&[
        "lab",
        "--lab",
        b.to_str().unwrap(),
        "claim",
        task,
        "--holder",
        "box-b",
        "--ttl",
        "600",
        "--identity",
        id2.to_str().unwrap()
    ])
    .status
    .success());
    assert!(cairn(&[
        "lab",
        "--lab",
        a.to_str().unwrap(),
        "sync",
        "--dir",
        b.to_str().unwrap()
    ])
    .status
    .success());
    let mut holders = Vec::new();
    for lab in [&a, &b] {
        let out = run(Command::new(bin())
            .arg("--root")
            .arg(repo())
            .args(["cairn-status", "XP-20261015-pkm-m3-n21-m4-n20", "--lab"])
            .arg(lab));
        assert!(out.status.success());
        let text = stdout(&out);
        let row = text
            .lines()
            .find(|l| l.starts_with(task))
            .expect("a row for the task")
            .to_string();
        let cols: Vec<&str> = row.split_whitespace().collect();
        assert_eq!(cols[1], "held", "{row}");
        assert_eq!(cols[3], "1", "one contended claim: {row}");
        holders.push(cols[2].to_string());
        // The other shard was never claimed.
        assert!(
            text.lines()
                .any(|l| l.starts_with("XP-20261015-pkm-m3-n21-m4-n20/shard-1")
                    && l.contains("unclaimed")),
            "{text}"
        );
    }
    assert_eq!(holders[0], holders[1], "both replicas agree on the holder");
    assert!(holders[0] == "box-a" || holders[0] == "box-b");
}

#[test]
fn publish_claim_exec_release_and_status_round_trip() {
    if !have_cairn() {
        return;
    }
    let w = tmp("publish");
    let (a, _b, id1, _id2) = two_replicas(&w);
    let lab = a.to_str().unwrap();
    let card = "XP-20261015-kummer-p224-depth";
    let publish = || {
        run(Command::new(bin()).arg("--root").arg(repo()).args([
            "cairn-publish",
            card,
            "--lab",
            lab,
            "--identity",
            id1.to_str().unwrap(),
        ]))
    };
    let first = stdout(&publish());
    assert!(
        first.contains(&format!("written   doc backlog/{card}.card")),
        "{first}"
    );
    assert!(
        first.contains(&format!("written   backlog/{card}/README.md")),
        "{first}"
    );
    // Publishing again signs nothing new.
    let second = stdout(&publish());
    assert!(!second.contains("written"), "{second}");
    let doc: serde_json::Value = serde_json::from_slice(
        &cairn(&[
            "lab",
            "--lab",
            lab,
            "doc",
            &format!("backlog/{card}"),
            "--json",
        ])
        .stdout,
    )
    .unwrap();
    assert_eq!(doc["card"]["id"], card);
    assert_eq!(doc["funding"].as_u64(), Some(400000));
    let readme = stdout(&cairn(&[
        "lab",
        "--lab",
        lab,
        "cat",
        &format!("backlog/{card}/README.md"),
    ]));
    assert_eq!(
        readme,
        fs::read_to_string(repo().join("backlog").join(card).join("README.md")).unwrap()
    );

    // Claim shard 0, run a stand-in command under it, release completed.
    let task = format!("{card}/shard-0");
    let claim = stdout(&cairn(&[
        "lab",
        "--lab",
        lab,
        "claim",
        &task,
        "--holder",
        "box-a",
        "--ttl",
        "600",
        "--identity",
        id1.to_str().unwrap(),
    ]));
    let claim_id = claim
        .split_whitespace()
        .nth(1)
        .expect("lease id")
        .to_string();
    let env = w.join("rootfs");
    fs::create_dir_all(&env).unwrap();
    fs::write(env.join("README"), "stand-in environment\n").unwrap();
    assert!(cairn(&[
        "lab",
        "--lab",
        lab,
        "env",
        "import",
        "stand-in",
        "--dir",
        env.to_str().unwrap(),
        "--identity",
        id1.to_str().unwrap()
    ])
    .status
    .success());
    let exec = cairn(&[
        "lab",
        "--lab",
        lab,
        "exec",
        "--env",
        "stand-in",
        "--identity",
        id1.to_str().unwrap(),
        "--sandbox",
        "none",
        "--task",
        &task,
        "--publish",
        &format!("runs/{card}/shard-0"),
        "--",
        "/bin/sh",
        "-c",
        "echo '{\"D\": 7}' > \"$CAIRN_LAB_OUT/rows.jsonl\"",
    ]);
    assert!(exec.status.success());
    assert!(cairn(&[
        "lab",
        "--lab",
        lab,
        "release",
        &claim_id,
        "--outcome",
        "completed",
        "--identity",
        id1.to_str().unwrap()
    ])
    .status
    .success());

    let out = stdout(&run(Command::new(bin()).arg("--root").arg(repo()).args([
        "cairn-status",
        card,
        "--lab",
        lab,
    ])));
    let row = out
        .lines()
        .find(|l| l.starts_with(&task))
        .expect("row")
        .to_string();
    assert!(
        row.contains("completed") && row.contains("exit=0") && row.contains("backend=none"),
        "{row}"
    );
    let other = out
        .lines()
        .find(|l| l.starts_with(&format!("{card}/shard-1")))
        .expect("row");
    assert!(
        other.contains("unclaimed") && other.trim_end().ends_with('-'),
        "{other}"
    );
}
