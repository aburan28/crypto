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
}

#[test]
fn list_orders_by_funding_then_id_including_unfunded_drafts() {
    // Use a fixed inventory: adding a committed card must not invalidate
    // the ordering test. Include zero funding and an equal-funding tie.
    let root = scratch_root("list-order");
    for (id, funding) in [
        ("XP-20261015-zero", 0),
        ("XP-20261015-tie-b", 800),
        ("XP-20261015-high", 1600),
        ("XP-20261015-tie-a", 800),
    ] {
        objective_card(&root, id);
        let path = root.join("backlog").join(id).join("card.json");
        let mut card: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();
        card["cairn"]["funding"]["units"] = serde_json::json!(funding);
        fs::write(path, serde_json::to_string_pretty(&card).unwrap()).unwrap();
    }
    let out = run(Command::new(bin()).arg("--root").arg(&root).arg("list"));
    assert!(out.status.success(), "{}", stdout(&out));
    let list = stdout(&out);
    let ids: Vec<&str> = list
        .lines()
        .map(|line| line.split_whitespace().next().unwrap())
        .collect();
    assert_eq!(
        ids,
        [
            "XP-20261015-high",
            "XP-20261015-tie-a",
            "XP-20261015-tie-b",
            "XP-20261015-zero",
        ],
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
fn cairn_post_is_a_dry_run_by_default_and_refuses_a_blocked_checker_or_a_lab_task() {
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
    // A card whose checker is blocked has nothing to post.
    let path = root.join("backlog/XP-20261015-demo-post/card.json");
    let mut card: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();
    card["cairn"]["checker"]["status"] = serde_json::json!("blocked");
    card["cairn"]["checker"]["blocked_by"] = serde_json::json!("pr:1 the verifier route");
    fs::write(&path, card.to_string()).unwrap();
    let p = Command::new(bin())
        .arg("--root")
        .arg(&root)
        .args(["cairn-post", "XP-20261015-demo-post"])
        .output()
        .unwrap();
    assert!(!p.status.success());
    assert!(String::from_utf8_lossy(&p.stderr).contains("pr:1 the verifier route"));
    // The P-192 card runs as lab tasks today: not an objective yet.
    let q = Command::new(bin())
        .arg("--root")
        .arg(repo())
        .args(["cairn-post", "XP-20261015-p192-isogeny-aux-search-t48"])
        .output()
        .unwrap();
    assert!(!q.status.success());
    assert!(String::from_utf8_lossy(&q.stderr).contains("lab-task"));
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
        // CI's lab job sets this, so a missing cairn there fails instead of skipping.
        assert!(
            std::env::var_os("BACKLOG_REQUIRE_CAIRN").is_none(),
            "BACKLOG_REQUIRE_CAIRN is set but `cairn` is not on PATH"
        );
        eprintln!("SKIPPED: `cairn` is not on PATH; the lab lease, publish and found tests need it (docs/CONDUCTOR.md: the cloud session hook installs it; CI's lab job installs a pinned build)");
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

// ── the native checker and the first-success rule ───────────────────────

/// A point of order exactly 2^t over F_p, by brute force with repeated
/// addition (independent of the checker's doubling chain).
fn small_witness(p: u64, t: u32) -> (u64, u64, u64, u64) {
    let inv = |v: u64| {
        let (mut r, mut b, mut e) = (1u64, v % p, p - 2);
        while e > 0 {
            if e & 1 == 1 {
                r = r * b % p;
            }
            b = b * b % p;
            e >>= 1;
        }
        r
    };
    let order = |a: u64, x: u64, y: u64| {
        let (mut cx, mut cy, mut n) = (x, y, 1u64);
        loop {
            let l = if cx == x {
                if (cy + y).is_multiple_of(p) {
                    return n + 1;
                }
                (3 * x % p * x % p + a) % p * inv(2 * y % p) % p
            } else {
                (y + p - cy) % p * inv((x + p - cx) % p) % p
            };
            let nx = (l * l % p + 2 * p - cx - x) % p;
            let ny = (l * ((cx + p - nx) % p) % p + p - cy) % p;
            cx = nx;
            cy = ny;
            n += 1;
        }
    };
    for a in 0..p {
        for b in 1..p {
            if (4 * a * a % p * a + 27 * b % p * b).is_multiple_of(p) {
                continue;
            }
            for x in 0..p {
                let r = (x * x % p * x + a * x + b) % p;
                if let Some(y) = (0..p).find(|y| y * y % p == r) {
                    if order(a, x, y) == 1 << t {
                        return (a, b, x, y);
                    }
                }
            }
        }
    }
    panic!("none")
}

fn write_artifact(path: &Path, p: u64, (a, b, x, y): (u64, u64, u64, u64), t: u32) {
    let v = serde_json::json!({"p": p.to_string(), "a": a.to_string(), "b": b.to_string(),
                               "Qx": x.to_string(), "Qy": y.to_string(), "t": t});
    fs::write(path, v.to_string()).unwrap();
}

fn checker() -> &'static str {
    env!("CARGO_BIN_EXE_backlog_check_aux2t")
}

#[test]
fn the_checker_binary_exit_codes() {
    let w = tmp("checker");
    let good = w.join("good.json");
    write_artifact(&good, 241, small_witness(241, 4), 4);
    let chk = |args: &[&str]| Command::new(checker()).args(args).output().unwrap();
    let g = good.to_str().unwrap();
    let out = chk(&[g, "--expect-t", "4", "--allow-test-prime"]);
    assert_eq!(out.status.code(), Some(0));
    let v: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["verdict"], "accept");
    assert_eq!(v["test_prime"], true);
    // --expect-t is required, and pins t.
    assert_eq!(chk(&[g, "--allow-test-prime"]).status.code(), Some(2));
    assert_eq!(
        chk(&[g, "--expect-t", "5", "--allow-test-prime"])
            .status
            .code(),
        Some(1)
    );
    assert_eq!(
        chk(&[g, "--expect-t", "x", "--allow-test-prime"])
            .status
            .code(),
        Some(2)
    );
    // The same curve is the wrong prime for the real check.
    let o = chk(&[g, "--expect-t", "4"]);
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stdout).contains("P-192"));
    // Malformed.
    fs::write(w.join("junk.json"), "{not json").unwrap();
    assert_eq!(
        chk(&[w.join("junk.json").to_str().unwrap(), "--expect-t", "4"])
            .status
            .code(),
        Some(2)
    );
    assert_eq!(chk(&[]).status.code(), Some(2));
}

/// Review regression (HIGH 1): a genuine point of order 2 over the P-192 prime,
/// claimed with t = 1, is rejected by the card's own checker argv.
#[test]
fn the_p192_cards_checker_pins_t_48() {
    let w = tmp("pin48");
    let art = w.join("t1.json");
    let p192 = "6277101735386680763835789423207666416083908700390324961279";
    fs::write(
        &art,
        serde_json::json!({"p": p192, "a": "1", "b": "0", "Qx": "0", "Qy": "0", "t": 1})
            .to_string(),
    )
    .unwrap();
    let card: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(
            repo().join("backlog/XP-20261015-p192-isogeny-aux-search-t48/card.json"),
        )
        .unwrap(),
    )
    .unwrap();
    for field in ["/cairn/checker/argv", "/outputs/verify"] {
        let argv: Vec<String> = card
            .pointer(field)
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .map(|s| {
                s.as_str()
                    .unwrap()
                    .replace("{artifact}", art.to_str().unwrap())
                    .replace("{out}/found/artifact.json", art.to_str().unwrap())
            })
            .collect();
        let out = Command::new(checker()).args(&argv[1..]).output().unwrap();
        assert_eq!(out.status.code(), Some(1), "{field}: {}", stdout(&out));
        assert!(
            stdout(&out).contains("t = 1, but this check is for t = 48"),
            "{}",
            stdout(&out)
        );
    }
}

#[test]
fn ids_are_checked_on_every_subcommand() {
    for args in [
        vec!["show", "../etc"],
        vec!["show", "XP-20261015-a/../b"],
        vec!["validate", "XP-20261015-x/y"],
        vec!["cairn-status", "../x", "--lab", "/nonexistent"],
        vec!["cairn-found", "/abs", "--lab", "/nonexistent"],
        vec!["cairn-post", "XP-2026-x"],
    ] {
        let o = Command::new(bin())
            .arg("--root")
            .arg(repo())
            .args(&args)
            .output()
            .unwrap();
        assert!(!o.status.success(), "{args:?}");
        let text = format!("{}{}", stdout(&o), String::from_utf8_lossy(&o.stderr));
        assert!(text.contains("is not a card id"), "{args:?}: {text}");
    }
}

#[test]
fn cairn_post_refuses_plain_http_nodes() {
    let root = scratch_root("node");
    objective_card(&root, "XP-20261015-demo-post");
    for (node, extra, ok_until_identity) in [
        ("http://example.com", vec![], false),
        ("http://localhost:9", vec![], false),
        ("ftp://x", vec![], false),
        ("https://", vec![], false),
        ("http://localhost:9", vec!["--allow-local-http"], true),
        ("https://node.example", vec![], true),
    ] {
        let o = Command::new(bin())
            .arg("--root")
            .arg(&root)
            .args([
                "cairn-post",
                "XP-20261015-demo-post",
                "--confirm",
                "--node",
                node,
            ])
            .args(&extra)
            .output()
            .unwrap();
        assert!(!o.status.success());
        let err = String::from_utf8_lossy(&o.stderr).to_string();
        assert_eq!(
            err.contains("--confirm needs --identity"),
            ok_until_identity,
            "{node}: {err}"
        );
    }
}

/// A draft test card with the P-192 card's shape on a small prime: the
/// checker argv carries --allow-test-prime and --expect-t 4.
fn test_found_card(root: &Path, id: &str) -> PathBuf {
    let src = repo().join("backlog/XP-20261015-p192-isogeny-aux-search-t48");
    let dir = root.join("backlog").join(id);
    fs::create_dir_all(&dir).unwrap();
    fs::copy(src.join("README.md"), dir.join("README.md")).unwrap();
    let mut card: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(src.join("card.json")).unwrap()).unwrap();
    card["id"] = serde_json::json!(id);
    card["status"] = serde_json::json!("draft");
    card["blocked_on"] = serde_json::json!([]);
    card["cairn"]["tasks"] = serde_json::json!([format!("{id}/shard-{{k}}")]);
    card["cairn"]["exec"]["publish"] = serde_json::json!(format!("runs/{id}"));
    card["cairn"]["checker"]["argv"] = serde_json::json!([
        "backlog_check_aux2t",
        "{artifact}",
        "--expect-t",
        "4",
        "--allow-test-prime"
    ]);
    card["outputs"]["verify"] = serde_json::json!([
        "backlog_check_aux2t",
        "{out}/found/artifact.json",
        "--expect-t",
        "4",
        "--allow-test-prime"
    ]);
    fs::write(dir.join("card.json"), card.to_string()).unwrap();
    dir
}

#[test]
fn the_stop_rule_fails_closed_on_a_missing_lab() {
    if !have_cairn() {
        return;
    }
    let o = Command::new(bin())
        .arg("--root")
        .arg(repo())
        .args([
            "cairn-found",
            "XP-20261015-p192-isogeny-aux-search-t48",
            "--lab",
            "/nonexistent/lab",
        ])
        .output()
        .unwrap();
    assert_eq!(
        o.status.code(),
        Some(2),
        "{}{}",
        stdout(&o),
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(String::from_utf8_lossy(&o.stderr).contains("no lab at"));
}

#[test]
fn every_seed_card_publishes_unchanged() {
    if !have_cairn() {
        return;
    }
    let w = tmp("publish-all");
    let (a, _b, id1, _id2) = two_replicas(&w);
    for entry in fs::read_dir(repo().join("backlog")).unwrap().flatten() {
        let id = entry.file_name().to_string_lossy().to_string();
        if !id.starts_with("XP-") {
            continue;
        }
        let o = run(Command::new(bin())
            .arg("--root")
            .arg(repo())
            .args(["cairn-publish", &id, "--lab"])
            .arg(&a)
            .arg("--identity")
            .arg(&id1));
        assert!(o.status.success(), "{id}");
        let doc: serde_json::Value = serde_json::from_slice(
            &cairn(&[
                "lab",
                "--lab",
                a.to_str().unwrap(),
                "doc",
                &format!("backlog/{id}"),
                "--json",
            ])
            .stdout,
        )
        .unwrap();
        let local: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(entry.path().join("card.json")).unwrap())
                .unwrap();
        assert_eq!(doc["card"], local, "{id} round-trips through the lab");
    }
}

#[test]
fn first_success_found_is_set_only_on_a_checked_published_artifact() {
    if !have_cairn() {
        return;
    }
    let w = tmp("found");
    let (a, _b, id1, _id2) = two_replicas(&w);
    let lab = a.to_str().unwrap();
    let ident = id1.to_str().unwrap();
    let root = scratch_root("found-root");
    let id = "XP-20261015-demo-found";
    let dir = test_found_card(&root, id);
    let backlog = |args: &[&str]| {
        Command::new(bin())
            .arg("--root")
            .arg(&root)
            .args(args)
            .output()
            .unwrap()
    };
    let say = |o: &Output| format!("{}{}", stdout(o), String::from_utf8_lossy(&o.stderr));

    let none = backlog(&["cairn-found", id, "--lab", lab]);
    assert_eq!(none.status.code(), Some(0), "{}", say(&none));
    assert!(stdout(&none).contains("not found"));
    assert!(
        backlog(&["cairn-publish", id, "--lab", lab, "--identity", ident])
            .status
            .success()
    );

    let task = format!("{id}/shard-7");
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
        ident,
    ]));
    let claim_id = claim.split_whitespace().nth(1).unwrap().to_string();
    let put = |path: &str, file: &Path| {
        assert!(cairn(&[
            "lab",
            "--lab",
            lab,
            "put",
            path,
            "--file",
            file.to_str().unwrap(),
            "--identity",
            ident
        ])
        .status
        .success());
    };
    let record = |lab_path: &str, claim: &str, holder: &str| {
        backlog(&[
            "cairn-found",
            id,
            "--lab",
            lab,
            "--lab-path",
            lab_path,
            "--claim",
            claim,
            "--holder",
            holder,
            "--identity",
            ident,
            "--checker-bin",
            checker(),
            "--allow-status",
        ])
    };
    let good = w.join("good.json");
    write_artifact(&good, 241, small_witness(241, 4), 4);
    let bad = w.join("bad.json");
    write_artifact(&bad, 241, small_witness(241, 4), 5);
    let good_path = format!("runs/{id}/shard-7/found/artifact.json");
    let bad_path = format!("runs/{id}/shard-7/found/wrong-t.json");
    put(&good_path, &good);
    put(&bad_path, &bad);

    // The checker sees the lab's bytes: a wrong t records nothing.
    assert_eq!(record(&bad_path, &claim_id, "box-a").status.code(), Some(1));
    // Paths outside runs/<ID>/, a lease that is not live for that holder, a missing artifact: errors.
    for (path, c, h) in [
        (format!("runs/{id}-x/a.json"), claim_id.as_str(), "box-a"),
        (format!("runs/{id}/../x/a.json"), claim_id.as_str(), "box-a"),
        (
            format!("backlog/{id}/README.md"),
            claim_id.as_str(),
            "box-a",
        ),
        (good_path.clone(), claim_id.as_str(), "box-b"),
        (good_path.clone(), "sha256:00", "box-a"),
        (
            format!("runs/{id}/shard-7/found/absent.json"),
            claim_id.as_str(),
            "box-a",
        ),
    ] {
        let o = record(&path, c, h);
        assert_eq!(o.status.code(), Some(2), "{path} {c} {h}: {}", say(&o));
    }
    // Without --allow-status a draft card records nothing.
    let o = backlog(&[
        "cairn-found",
        id,
        "--lab",
        lab,
        "--lab-path",
        &good_path,
        "--claim",
        &claim_id,
        "--holder",
        "box-a",
        "--identity",
        ident,
        "--checker-bin",
        checker(),
    ]);
    assert_eq!(o.status.code(), Some(2));
    assert!(say(&o).contains("only a ready or running card"));
    // A local card that differs from the published one is refused.
    let card_file = dir.join("card.json");
    let original = fs::read_to_string(&card_file).unwrap();
    fs::write(
        &card_file,
        original.replace("\"--expect-t\",\"4\"", "\"--expect-t\",\"3\""),
    )
    .unwrap();
    let o = record(&good_path, &claim_id, "box-a");
    assert_eq!(o.status.code(), Some(2));
    assert!(say(&o).contains("differs from doc"), "{}", say(&o));
    fs::write(&card_file, &original).unwrap();

    // The checked, published witness sets `found`; every runner then stops.
    let ok = record(&good_path, &claim_id, "box-a");
    assert_eq!(ok.status.code(), Some(0), "{}", say(&ok));
    let stop = backlog(&["cairn-found", id, "--lab", lab]);
    assert_eq!(stop.status.code(), Some(3));
    let text = stdout(&stop);
    assert!(
        text.contains("box-a")
            && text.contains("accept")
            && text.contains("checker_sha256")
            && text.contains(&task),
        "{text}"
    );
    assert!(cairn(&[
        "lab",
        "--lab",
        lab,
        "release",
        &claim_id,
        "--outcome",
        "completed",
        "--identity",
        ident
    ])
    .status
    .success());
    let st = stdout(&backlog(&["cairn-status", id, "--lab", lab]));
    assert!(
        st.contains("FOUND") && st.contains("1 of 65536 shards"),
        "{st}"
    );
    assert!(
        st.lines()
            .any(|l| l.starts_with(&task) && l.contains("completed")),
        "{st}"
    );
    assert_eq!(
        record(&good_path, &claim_id, "box-a").status.code(),
        Some(3)
    );
}

/// Review regression (MEDIUM 2b): two replicas record `found` concurrently.
/// cairn keeps both values as a register conflict; the stop rule must still
/// read it as found.
#[test]
fn a_conflicted_found_still_stops_every_runner() {
    if !have_cairn() {
        return;
    }
    let w = tmp("conflict");
    let (a, b, id1, id2) = two_replicas(&w);
    let root = scratch_root("conflict-root");
    let id = "XP-20261015-demo-conflict";
    test_found_card(&root, id);
    let doc = format!("backlog/{id}");
    assert!(cairn(&[
        "lab",
        "--lab",
        a.to_str().unwrap(),
        "doc",
        "set",
        &doc,
        "found",
        r#"{"holder":"box-a"}"#,
        "--identity",
        id1.to_str().unwrap()
    ])
    .status
    .success());
    assert!(cairn(&[
        "lab",
        "--lab",
        b.to_str().unwrap(),
        "doc",
        "set",
        &doc,
        "found",
        r#"{"holder":"box-b"}"#,
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
    for lab in [&a, &b] {
        let o = Command::new(bin())
            .arg("--root")
            .arg(&root)
            .args(["cairn-found", id, "--lab"])
            .arg(lab)
            .output()
            .unwrap();
        assert_eq!(o.status.code(), Some(3), "{}", stdout(&o));
        let text = stdout(&o);
        assert!(
            text.contains("conflict") && text.contains("box-a") && text.contains("box-b"),
            "{text}"
        );
        let st = stdout(
            &Command::new(bin())
                .arg("--root")
                .arg(&root)
                .args(["cairn-status", id, "--lab"])
                .arg(lab)
                .output()
                .unwrap(),
        );
        assert!(st.contains("FOUND (conflict: 2"), "{st}");
    }
}

/// Review regression (MEDIUM 2): a cairn outside the supported range is refused.
#[test]
fn an_unsupported_cairn_version_is_refused() {
    use std::os::unix::fs::PermissionsExt;
    let w = tmp("version");
    let fake = w.join("cairn");
    fs::write(&fake, "#!/bin/sh\necho 'cairn 7.9.0'\n").unwrap();
    fs::set_permissions(&fake, fs::Permissions::from_mode(0o755)).unwrap();
    let o = Command::new(bin())
        .env("CAIRN", &fake)
        .arg("--root")
        .arg(repo())
        .args([
            "cairn-found",
            "XP-20261015-p192-isogeny-aux-search-t48",
            "--lab",
        ])
        .arg(&w)
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("outside the supported range"));
}
