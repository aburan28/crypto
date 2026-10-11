//! End to end with the local provider: the real controller, agent, store
//! and checkpoint code, with "instances" that are agent processes.

use serde_json::Value;
use sha2::{Digest, Sha256};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

const BIN: &str = env!("CARGO_BIN_EXE_spotlab");

struct Lab {
    root: PathBuf,
    cfg: PathBuf,
}

impl Lab {
    fn new(tag: &str) -> Lab {
        let root = std::env::temp_dir().join(format!("spotlab-e2e-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).unwrap();
        let cfg = root.join("spotlab.json");
        let doc = serde_json::json!({
            "store": root.join("store"),
            "local": {
                "agent_binary": BIN,
                "work_root": root.join("work"),
                "catalog": [
                    {"instance_type": "small", "vcpus": 2, "memory_gb": 4.0, "arch": std::env::consts::ARCH, "spot_hourly_usd": 0.01},
                    {"instance_type": "big", "vcpus": 8, "memory_gb": 16.0, "arch": std::env::consts::ARCH, "spot_hourly_usd": 0.05}
                ]
            },
            "budget": {"max_running": 4, "max_hourly_usd": 1.0}
        });
        std::fs::write(&cfg, serde_json::to_vec_pretty(&doc).unwrap()).unwrap();
        Lab { root, cfg }
    }

    fn cli(&self, args: &[&str]) -> (bool, String) {
        let out = Command::new(BIN)
            .arg("--config")
            .arg(&self.cfg)
            .args(args)
            .output()
            .unwrap();
        let text = format!(
            "{}{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        (out.status.success(), text)
    }

    fn json(&self, args: &[&str]) -> Value {
        let (ok, text) = self.cli(args);
        assert!(ok, "{args:?} failed: {text}");
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("{args:?}: {e}: {text}"))
    }

    fn submit_demo(&self, steps: u64, step_us: u64) -> String {
        let (ok, text) = self.cli(&[
            "run",
            "--name",
            "demo",
            "--vcpus",
            "2",
            "--interval",
            "10",
            "--",
            BIN,
            "demo-job",
            "--to",
            &steps.to_string(),
            "--every",
            "500",
            "--step-us",
            &step_us.to_string(),
        ]);
        assert!(ok, "run failed: {text}");
        let id = text
            .lines()
            .find_map(|l| l.strip_prefix("submitted "))
            .expect(&text)
            .trim()
            .to_string();
        assert!(text.contains("attempt 1 launched on local small"), "{text}");
        id
    }

    fn state(&self, id: &str) -> Value {
        self.json(&["status", id])
    }

    /// Reconcile until `pred` holds on the status, or fail after `secs`.
    fn until(&self, id: &str, secs: u64, what: &str, pred: impl Fn(&Value) -> bool) -> Value {
        let t0 = Instant::now();
        loop {
            let (ok, text) = self.cli(&["reconcile"]);
            assert!(ok, "reconcile failed: {text}");
            let st = self.state(id);
            if pred(&st) {
                return st;
            }
            if t0.elapsed() > Duration::from_secs(secs) {
                panic!(
                    "timed out waiting for {what}: {st:#}\nagent logs:\n{}",
                    self.agent_logs()
                );
            }
            std::thread::sleep(Duration::from_millis(300));
        }
    }

    fn agent_logs(&self) -> String {
        let mut s = String::new();
        if let Ok(rd) = std::fs::read_dir(self.root.join("work")) {
            for e in rd.flatten() {
                if e.path().extension().is_some_and(|x| x == "log") {
                    s += &format!(
                        "== {}\n{}\n",
                        e.path().display(),
                        std::fs::read_to_string(e.path()).unwrap_or_default()
                    );
                }
            }
        }
        s
    }
}

fn expected_hash(n: u64) -> String {
    let mut h: [u8; 32] = Sha256::digest(b"spotlab").into();
    for i in 1..=n {
        let mut d = Sha256::new();
        d.update(h);
        d.update(i.to_le_bytes());
        h = d.finalize().into();
    }
    hex::encode(h)
}

fn running_with_progress(st: &Value, min_starts: u64) -> bool {
    st["record"]["state"] == "running"
        && st["record"]["attempt"].as_u64().unwrap_or(0) >= min_starts
}

fn read_json(p: &Path) -> Value {
    serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap()
}

#[test]
fn preempted_and_lost_attempts_resume_to_the_uninterrupted_answer() {
    let lab = Lab::new("resume");
    let steps = 60_000;
    let id = lab.submit_demo(steps, 100); // about 6 s of work per uninterrupted run

    // Attempt 1: let it work, then deliver an interruption notice.
    lab.until(&id, 30, "attempt 1 running", |s| {
        running_with_progress(s, 1)
    });
    std::thread::sleep(Duration::from_millis(1500));
    let (ok, text) = lab.cli(&["local-preempt", &id]);
    assert!(ok, "{text}");
    let st = lab.until(&id, 60, "attempt 2 launched after the preemption", |s| {
        running_with_progress(s, 2)
    });
    let a1 = &st["record"]["attempts"][0];
    assert_eq!(a1["end_reason"], "preempted", "{st:#}");
    assert!(
        !st["checkpoints"].as_array().unwrap().is_empty(),
        "the preemption must leave a checkpoint: {st:#}"
    );

    // Attempt 2: kill the agent and the command outright, as if the host vanished.
    std::thread::sleep(Duration::from_millis(1500));
    let pid: i32 = st["record"]["attempts"][1]["instance"]["instance_id"]
        .as_str()
        .unwrap()["pid-".len()..]
        .parse()
        .unwrap();
    unsafe {
        libc::kill(-pid, libc::SIGKILL);
        libc::kill(pid, libc::SIGKILL);
    }
    let st = lab.until(&id, 60, "attempt 3 after the lost instance", |s| {
        running_with_progress(s, 3)
    });
    assert!(
        st["record"]["attempts"][1]["end_reason"]
            .as_str()
            .unwrap()
            .starts_with("lost"),
        "{st:#}"
    );

    // Attempt 3 runs to the end.
    lab.until(&id, 120, "success", |s| s["record"]["state"] == "succeeded");
    let res = lab.json(&["result", &id]);
    assert_eq!(res["status"], "succeeded");
    assert_eq!(res["exit_code"], 0);
    assert_eq!(res["resumed"], true);
    assert_eq!(res["preemptions"], 2, "{res:#}");
    assert_eq!(res["attempts"].as_array().unwrap().len(), 3);
    assert!(res["evidence"]
        .as_str()
        .unwrap()
        .contains("not performance measurements"));

    let dest = lab.root.join("fetched");
    let f = lab.json(&["fetch", &id, "--dest", dest.to_str().unwrap()]);
    assert_eq!(f["verified"], true);
    let out = read_json(&dest.join("result.json"));
    assert_eq!(out["steps"], steps);
    assert_eq!(
        out["hash"],
        expected_hash(steps),
        "resumed chain differs from the uninterrupted one"
    );
    assert!(
        out["starts"].as_u64().unwrap() >= 2,
        "the command must have resumed: {out}"
    );

    // Nothing is left running.
    let ov = lab.json(&["overview"]);
    assert_eq!(ov["running"].as_array().unwrap().len(), 0);
    std::fs::remove_dir_all(&lab.root).ok();
}

#[test]
fn cancel_terminates_and_writes_a_result() {
    let lab = Lab::new("cancel");
    let id = lab.submit_demo(10_000_000, 100);
    lab.until(&id, 30, "running", |s| running_with_progress(s, 1));
    let (ok, text) = lab.cli(&["cancel", &id]);
    assert!(ok, "{text}");
    lab.until(&id, 30, "cancelled", |s| {
        s["record"]["state"] == "cancelled"
    });
    let res = lab.json(&["result", &id]);
    assert_eq!(res["status"], "cancelled");
    assert_eq!(
        lab.json(&["overview"])["running"].as_array().unwrap().len(),
        0
    );
    std::fs::remove_dir_all(&lab.root).ok();
}

#[test]
fn offers_respect_resources_and_price_cap() {
    let lab = Lab::new("offers");
    let o = lab.json(&["offers", "--vcpus", "4"]);
    let offers = o["offers"].as_array().unwrap();
    assert_eq!(offers.len(), 1);
    assert_eq!(offers[0]["instance_type"], "big");
    let o = lab.json(&["offers", "--vcpus", "1", "--max-hourly", "0.02"]);
    assert_eq!(o["offers"].as_array().unwrap().len(), 1);
    assert_eq!(o["offers"][0]["instance_type"], "small");
    std::fs::remove_dir_all(&lab.root).ok();
}

#[test]
fn mcp_lists_tools_and_answers_overview() {
    use std::io::{BufRead, BufReader, Write};
    let lab = Lab::new("mcp");
    let mut child = Command::new(BIN)
        .arg("--config")
        .arg(&lab.cfg)
        .args(["mcp", "--reconcile-every", "0"])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    let mut out = BufReader::new(child.stdout.take().unwrap());
    let mut ask = |v: Value| -> Value {
        writeln!(stdin, "{v}").unwrap();
        let mut line = String::new();
        out.read_line(&mut line).unwrap();
        serde_json::from_str(&line).unwrap()
    };
    let init = ask(
        serde_json::json!({"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18"}}),
    );
    assert_eq!(init["result"]["serverInfo"]["name"], "spotlab");
    let tools = ask(serde_json::json!({"jsonrpc":"2.0","id":2,"method":"tools/list"}));
    let names: Vec<&str> = tools["result"]["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t["name"].as_str().unwrap())
        .collect();
    for n in [
        "spotlab_overview",
        "spotlab_offers",
        "spotlab_submit",
        "spotlab_status",
        "spotlab_result",
        "spotlab_cancel",
    ] {
        assert!(names.contains(&n), "{names:?}");
    }
    let ov = ask(
        serde_json::json!({"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"spotlab_overview","arguments":{}}}),
    );
    assert_eq!(ov["result"]["isError"], false, "{ov:#}");
    let sub = ask(
        serde_json::json!({"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"spotlab_submit","arguments":{
        "name":"mcp-demo","argv":[BIN,"demo-job","--to","1000"],"vcpus":1}}}),
    );
    assert_eq!(sub["result"]["isError"], false, "{sub:#}");
    let body: Value =
        serde_json::from_str(sub["result"]["content"][0]["text"].as_str().unwrap()).unwrap();
    let id = body["job_id"].as_str().unwrap().to_string();
    assert!(
        body["events"][0].as_str().unwrap().contains("launched"),
        "{body:#}"
    );
    drop(stdin);
    child.wait().unwrap();
    lab.until(&id, 30, "mcp job success", |s| {
        s["record"]["state"] == "succeeded"
    });
    std::fs::remove_dir_all(&lab.root).ok();
}

#[test]
fn example_config_parses_and_unpriced_gcp_entries_are_not_offered() {
    let cfg = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/spotlab.json");
    let tmp = std::env::temp_dir().join(format!("spotlab-example-{}", std::process::id()));
    // Point the store at a local directory; everything else as shipped.
    let mut doc: Value = serde_json::from_slice(&std::fs::read(&cfg).unwrap()).unwrap();
    doc["store"] = Value::String(tmp.join("store").display().to_string());
    std::fs::create_dir_all(&tmp).unwrap();
    let p = tmp.join("spotlab.json");
    std::fs::write(&p, serde_json::to_vec(&doc).unwrap()).unwrap();
    let out = Command::new(BIN)
        .arg("--config")
        .arg(&p)
        .args(["offers", "--provider", "gcp"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["offers"].as_array().unwrap().len(), 0);
    std::fs::remove_dir_all(&tmp).ok();
}
