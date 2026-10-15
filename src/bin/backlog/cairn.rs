//! Thin shell-outs to `cairn lab`. Nothing here reimplements a lab operation:
//! every write is a `cairn lab` command signed with the caller's identity, and
//! every read is `cairn lab … --json`. `CAIRN` overrides the executable.

use std::path::Path;
use std::process::Command;

use serde_json::Value;

fn cairn() -> String {
    std::env::var("CAIRN").unwrap_or_else(|_| "cairn".to_string())
}

/// Run `cairn lab --lab DIR <args>` and return stdout, or the command's error.
pub fn lab(dir: &Path, args: &[&str]) -> Result<String, String> {
    let out = Command::new(cairn())
        .arg("lab")
        .arg("--lab")
        .arg(dir)
        .args(args)
        .output()
        .map_err(|e| {
            format!(
                "cannot run {}: {e} (is cairn on PATH? set CAIRN to its path)",
                cairn()
            )
        })?;
    if out.status.success() {
        Ok(String::from_utf8_lossy(&out.stdout).into_owned())
    } else {
        Err(format!(
            "cairn lab {} exited {}: {}",
            args.first().unwrap_or(&""),
            out.status
                .code()
                .map_or("on a signal".to_string(), |c| c.to_string()),
            String::from_utf8_lossy(&out.stderr).trim()
        ))
    }
}

/// What publishing would change, and the result of doing it.
#[derive(Debug, Default)]
pub struct Published {
    pub written: Vec<String>,
    pub unchanged: Vec<String>,
}

/// Write the card into the lab as doc `backlog/<ID>` (field `card`) and its
/// README/PROTOCOL as lab files under `backlog/<ID>/`. Unchanged content is
/// not rewritten, so publishing twice signs nothing new.
pub fn publish(
    lab_dir: &Path,
    identity: &Path,
    card_dir: &Path,
    id: &str,
    card: &Value,
) -> Result<Published, String> {
    let mut done = Published::default();
    let doc = format!("backlog/{id}");
    // A doc nobody has written yet answers "not found": publish creates it.
    let current: Value = match lab(lab_dir, &["doc", &doc, "--json"]) {
        Ok(text) => serde_json::from_str(&text).unwrap_or(Value::Null),
        Err(e) if e.contains("not found") => Value::Null,
        Err(e) => return Err(e),
    };
    let ident = identity.to_string_lossy().into_owned();
    if current.get("card") == Some(card) {
        done.unchanged.push(format!("doc {doc}.card"));
    } else {
        let json = serde_json::to_string(card).map_err(|e| e.to_string())?;
        lab(
            lab_dir,
            &["doc", "set", &doc, "card", &json, "--identity", &ident],
        )?;
        done.written.push(format!("doc {doc}.card"));
    }
    // Funding and status as their own fields, so a runner ordering the
    // space's work reads two small registers rather than every card.
    for (key, value) in [
        (
            "funding",
            card.pointer("/cairn/funding/units")
                .cloned()
                .unwrap_or(Value::Null),
        ),
        ("status", card.get("status").cloned().unwrap_or(Value::Null)),
    ] {
        if current.get(key) == Some(&value) {
            done.unchanged.push(format!("doc {doc}.{key}"));
        } else {
            lab(
                lab_dir,
                &[
                    "doc",
                    "set",
                    &doc,
                    key,
                    &value.to_string(),
                    "--identity",
                    &ident,
                ],
            )?;
            done.written.push(format!("doc {doc}.{key}"));
        }
    }
    for name in ["README.md", "PROTOCOL.md"] {
        let file = card_dir.join(name);
        if !file.is_file() {
            continue;
        }
        let path = format!("backlog/{id}/{name}");
        let bytes = std::fs::read(&file).map_err(|e| format!("{}: {e}", file.display()))?;
        let listed: Value =
            serde_json::from_str(&lab(lab_dir, &["ls", &path, "--json"])?).unwrap_or(Value::Null);
        let exists = listed.as_array().is_some_and(|a| {
            a.iter()
                .any(|e| e.get("path").and_then(Value::as_str) == Some(path.as_str()))
        });
        let same = exists
            && lab(lab_dir, &["cat", &path])
                .map(|s| s.as_bytes() == bytes.as_slice())
                .unwrap_or(false);
        if same {
            done.unchanged.push(path);
            continue;
        }
        let file_arg = file.to_string_lossy().into_owned();
        let mut args = vec![
            "put",
            path.as_str(),
            "--file",
            file_arg.as_str(),
            "--identity",
            ident.as_str(),
        ];
        if !exists {
            args.push("--create");
        }
        lab(lab_dir, &args)?;
        done.written.push(path);
    }
    Ok(done)
}

/// One row of `backlog cairn-status`.
#[derive(Debug)]
pub struct TaskRow {
    pub task: String,
    pub state: String,
    pub holder: String,
    pub contended: u64,
    pub expired: u64,
    pub released: u64,
    pub runs: Vec<String>,
}

/// Per-task lease state and run receipts for a card's tasks, from
/// `cairn lab tasks --json` and `cairn lab runs --json`.
pub fn status(lab_dir: &Path, tasks: &[String]) -> Result<Vec<TaskRow>, String> {
    let leases: Value = serde_json::from_str(&lab(lab_dir, &["tasks", "--json"])?)
        .map_err(|e| format!("cairn lab tasks --json is not JSON: {e}"))?;
    let runs: Value = serde_json::from_str(&lab(lab_dir, &["runs", "--json"])?)
        .map_err(|e| format!("cairn lab runs --json is not JSON: {e}"))?;
    let mut rows = Vec::new();
    for task in tasks {
        let lease = leases.get(task);
        let n = |k: &str| {
            lease
                .and_then(|l| l.get(k))
                .and_then(Value::as_u64)
                .unwrap_or(0)
        };
        let holder = lease
            .and_then(|l| l.get("holder"))
            .and_then(|h| h.get("holder"))
            .and_then(Value::as_str)
            .unwrap_or("-")
            .to_string();
        let completed = lease
            .and_then(|l| l.get("completed"))
            .and_then(Value::as_bool)
            .unwrap_or(false);
        let state = match lease {
            None => "unclaimed",
            Some(_) if completed => "completed",
            Some(_) if holder != "-" => "held",
            Some(_) if n("expired") > 0 => "expired",
            Some(_) => "released",
        }
        .to_string();
        let receipts = runs
            .as_array()
            .map(|a| {
                a.iter()
                    .filter(|r| {
                        r.get("receipt")
                            .and_then(|x| x.get("task"))
                            .and_then(Value::as_str)
                            == Some(task.as_str())
                    })
                    .map(|r| {
                        let rc = r.get("receipt");
                        let exit = rc
                            .and_then(|x| x.get("exit_status"))
                            .map_or("?".to_string(), |v| v.to_string());
                        let backend = rc
                            .and_then(|x| x.get("backend"))
                            .and_then(Value::as_str)
                            .unwrap_or("?");
                        let id = r.get("id").and_then(Value::as_str).unwrap_or("?");
                        format!("{} exit={exit} backend={backend}", short(id))
                    })
                    .collect()
            })
            .unwrap_or_default();
        rows.push(TaskRow {
            task: task.clone(),
            state,
            holder,
            contended: n("contended"),
            expired: n("expired"),
            released: n("released"),
            runs: receipts,
        });
    }
    Ok(rows)
}

fn short(id: &str) -> &str {
    let s = id.strip_prefix("sha256:").unwrap_or(id);
    &s[..s.len().min(12)]
}

/// The cairn objective draft for an objective-mode card: the record
/// `cairn sign-objective` signs and `POST /submit?kind=objective` admits.
/// The card's funding is the reward; with piecework the reward is the pool and
/// each novel unit pays `reward / units`.
pub fn objective_draft(card: &Value, created_at: &str) -> Result<Value, String> {
    let c = card.get("cairn").ok_or("the card has no cairn block")?;
    if c.get("mode").and_then(Value::as_str) != Some("objective") {
        return Err("only an objective-mode card is posted; a lab-task card's funding orders runners through the lab doc (backlog cairn-publish)".into());
    }
    if c.pointer("/checker/status").and_then(Value::as_str) != Some("available") {
        return Err(format!(
            "the card's checker is not available (blocked by {}); nothing can be posted until a native verifier runs under cairn",
            c.pointer("/checker/blocked_by").and_then(Value::as_str).unwrap_or("an unnamed blocker")
        ));
    }
    let verifier = match c.get("verifier") {
        Some(v @ Value::Object(_)) => v.clone(),
        _ => return Err("cairn.verifier is null: there is no verifier block to post".into()),
    };
    let reward = c
        .pointer("/funding/units")
        .and_then(Value::as_u64)
        .ok_or("cairn.funding.units is missing")?;
    let title = card.get("title").and_then(Value::as_str).unwrap_or("");
    let question = card.get("question").and_then(Value::as_str).unwrap_or("");
    let statement = c
        .get("statement")
        .and_then(Value::as_str)
        .map(String::from)
        .unwrap_or_else(|| format!("{title}. {question}"));
    let mut draft = serde_json::json!({
        "goal": c.get("goal").and_then(Value::as_str).ok_or("cairn.goal is missing")?,
        "statement": statement,
        "verifier": verifier,
        "reward": reward,
        "funder": "",
        "created_at": created_at,
    });
    if let Some(schema) = c.pointer("/checker/artifact_schema") {
        draft["artifact_schema"] = schema.clone();
    }
    if let Some(units) = c.pointer("/piecework/units").and_then(Value::as_u64) {
        let price = reward / units;
        if price == 0 {
            return Err(format!(
                "funding {reward} is below one unit per piece ({units} pieces)"
            ));
        }
        draft["piecework"] = serde_json::json!({"units": units, "unit_price": price});
    }
    Ok(draft)
}

/// Sign the draft with `cairn sign-objective` and send it to a node's
/// `POST /submit?kind=objective` (curl). Returns the node's answer.
pub fn post(draft: &Value, node: &str, identity: &Path, scratch: &Path) -> Result<String, String> {
    let draft_file = scratch.join("objective.draft.json");
    std::fs::write(
        &draft_file,
        serde_json::to_string(draft).map_err(|e| e.to_string())?,
    )
    .map_err(|e| format!("{}: {e}", draft_file.display()))?;
    let out = Command::new(cairn())
        .arg("sign-objective")
        .arg(&draft_file)
        .arg("--identity")
        .arg(identity)
        .output()
        .map_err(|e| format!("cannot run cairn sign-objective: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "cairn sign-objective failed: {}",
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }
    let signed_file = scratch.join("objective.signed.json");
    std::fs::write(&signed_file, &out.stdout).map_err(|e| e.to_string())?;
    let url = format!("{}/submit?kind=objective", node.trim_end_matches('/'));
    let resp = Command::new("curl")
        .args([
            "-sS",
            "--fail-with-body",
            "-X",
            "POST",
            "-H",
            "content-type: application/json",
            "--data-binary",
        ])
        .arg(format!("@{}", signed_file.display()))
        .arg(&url)
        .output()
        .map_err(|e| format!("cannot run curl: {e}"))?;
    let body = String::from_utf8_lossy(&resp.stdout).into_owned();
    if resp.status.success() {
        Ok(body)
    } else {
        Err(format!(
            "{url} refused the objective: {body} {}",
            String::from_utf8_lossy(&resp.stderr).trim()
        ))
    }
}
