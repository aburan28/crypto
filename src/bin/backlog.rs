//! `backlog` — the experiment backlog: cards in `backlog/<ID>/`, validated
//! here, executed and leased through `cairn lab`.
//!
//! ```text
//! backlog list [--status S]...            one line per card
//! backlog show <ID> [--shard K]           a card in full, with shard K's resolved argv
//! backlog new <slug> [--title T]          scaffold backlog/XP-<today>-<slug>/
//! backlog validate [<ID>... | --all]      the card rules (exit 1 on any error)
//! backlog cairn-publish <ID> --lab DIR --identity FILE
//!                                         the card as lab doc backlog/<ID>, README/PROTOCOL as lab files
//! backlog cairn-status <ID> --lab DIR     per-task lease holder/state and run receipts
//! backlog cairn-post <ID> [--confirm --node URL --identity FILE]
//!                                         an objective-mode card as a cairn objective (dry run by default)
//! backlog cairn-found <ID> --lab DIR [--lab-path P --claim C --holder H --identity I --checker-bin B]
//!                                         first_success: is `found` set (exit 3)? / check the lab's artifact and record it
//! ```
//!
//! The protocol is `backlog/README.md`.

#[path = "backlog/cairn.rs"]
mod cairn;
#[path = "backlog/card.rs"]
mod card;

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, Subcommand};
use serde_json::Value;

#[derive(Parser)]
#[command(
    name = "backlog",
    version,
    about = "Experiment cards: list, show, scaffold, validate, publish to a cairn lab"
)]
struct Cli {
    /// Repository root (default: the nearest ancestor of the working directory holding backlog/schema/).
    #[arg(long, global = true)]
    root: Option<PathBuf>,
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// One line per card, highest funding first: id, status, funding, title.
    List {
        /// Only these statuses (repeatable).
        #[arg(long)]
        status: Vec<String>,
    },
    /// Print a card, its blockers and shard K's resolved command.
    Show {
        id: String,
        #[arg(long, default_value_t = 0)]
        shard: u64,
    },
    /// Scaffold a draft card: backlog/XP-<date>-<slug>/{card.json,README.md}.
    New {
        slug: String,
        #[arg(long)]
        title: Option<String>,
        /// YYYY-MM-DD (default: today, UTC).
        #[arg(long)]
        date: Option<String>,
    },
    /// Check cards against backlog.card/v1. Exit 1 on any error.
    Validate {
        ids: Vec<String>,
        #[arg(long)]
        all: bool,
        /// Do not ask git whether a frozen protocol's blob is committed.
        #[arg(long)]
        no_git: bool,
        /// Also require every ready/running/done crypto card's code.commit to be
        /// an ancestor of this ref (CI passes origin/main). Locally only the
        /// 40-hex form is checked.
        #[arg(long)]
        merged_into: Option<String>,
    },
    /// Write the card into a cairn lab as doc backlog/<ID> and its README/PROTOCOL as files.
    CairnPublish {
        id: String,
        #[arg(long)]
        lab: PathBuf,
        #[arg(long)]
        identity: PathBuf,
    },
    /// Post an objective-mode card to a cairn node, funded with the card's
    /// funding. Prints the draft and posts nothing unless --confirm is given.
    CairnPost {
        id: String,
        #[arg(long)]
        node: Option<String>,
        #[arg(long)]
        identity: Option<PathBuf>,
        /// Print what would be posted (the default).
        #[arg(long)]
        dry_run: bool,
        /// Sign with --identity and send to --node. Funding reorders everyone's work: never implicit.
        #[arg(long, conflicts_with = "dry_run")]
        confirm: bool,
        /// Allow plain http to 127.0.0.1/localhost (a local node); otherwise https only.
        #[arg(long)]
        allow_local_http: bool,
    },
    /// Lease holder/state and run receipts for the card's cairn tasks.
    CairnStatus {
        id: String,
        #[arg(long)]
        lab: PathBuf,
    },
    /// first_success stop rule. Without --lab-path: is `found` set? (exit 3
    /// if so, or if it is conflicted: stop claiming; 0 if not; 2 on any
    /// error). With --lab-path: fetch the published artifact from the lab, run
    /// the pinned checker on those bytes, and on accept record `found`.
    CairnFound {
        id: String,
        #[arg(long)]
        lab: PathBuf,
        /// Where the run published the artifact: a lab path under runs/<ID>/.
        #[arg(long, requires_all = ["claim", "holder", "identity", "checker_bin"])]
        lab_path: Option<String>,
        /// The live lease the finding run holds on one of the card's tasks.
        #[arg(long)]
        claim: Option<String>,
        /// That lease's holder.
        #[arg(long)]
        holder: Option<String>,
        #[arg(long)]
        identity: Option<PathBuf>,
        /// The checker binary from the card's pinned build; its sha256 is recorded in `found`.
        #[arg(long)]
        checker_bin: Option<PathBuf>,
        /// Tests only: record on a card that is not ready or running.
        #[arg(long)]
        allow_status: bool,
    },
}

fn find_root(explicit: Option<PathBuf>) -> Result<PathBuf, String> {
    if let Some(r) = explicit {
        return Ok(r);
    }
    let mut at = std::env::current_dir().map_err(|e| e.to_string())?;
    loop {
        if at.join("backlog/schema/card.v1.json").is_file() {
            return Ok(at);
        }
        if !at.pop() {
            return Err(
                "no backlog/schema/card.v1.json above the working directory; pass --root".into(),
            );
        }
    }
}

fn card_dir(root: &Path, id: &str) -> PathBuf {
    root.join("backlog").join(id)
}

fn ctx(root: &Path, id: &str, check_git: bool) -> card::Ctx {
    card::Ctx {
        root: root.to_path_buf(),
        dir: card_dir(root, id),
        known_ids: card::card_ids(root),
        check_git,
    }
}

fn load(root: &Path, id: &str) -> Result<Value, String> {
    // An id names a directory: refuse anything that is not one (no `/`, no `..`).
    if !card::valid_id(id) {
        return Err(format!("{id:?} is not a card id (XP-<YYYYMMDD>-<slug>)"));
    }
    card::load(&card_dir(root, id).join("card.json"))
}

fn s<'a>(v: &'a Value, path: &str) -> &'a str {
    let mut at = v;
    for p in path.split('.') {
        match at.get(p) {
            Some(x) => at = x,
            None => return "",
        }
    }
    at.as_str().unwrap_or("")
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(cli) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("backlog: {e}");
            ExitCode::from(2)
        }
    }
}

fn run(cli: Cli) -> Result<ExitCode, String> {
    let root = find_root(cli.root)?;
    match cli.cmd {
        Cmd::List { status } => {
            let mut rows: Vec<(i64, String, String, String)> = Vec::new();
            for id in card::card_ids(&root) {
                let c = match load(&root, &id) {
                    Ok(c) => c,
                    Err(e) => {
                        rows.push((i64::MIN, id, "unreadable".into(), e));
                        continue;
                    }
                };
                let st = s(&c, "status").to_string();
                if !status.is_empty() && !status.contains(&st) {
                    continue;
                }
                let prio = card::funding(&c) as i64;
                rows.push((prio, id, st, s(&c, "title").to_string()));
            }
            rows.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
            for (prio, id, st, title) in rows {
                println!("{id:<48} {st:<9} {prio:>9}  {title}");
            }
            Ok(ExitCode::SUCCESS)
        }
        Cmd::Show { id, shard } => {
            let c = load(&root, &id)?;
            println!("{id}  [{}]  {}", s(&c, "status"), s(&c, "title"));
            if !s(&c, "status_note").is_empty() {
                println!("  note: {}", s(&c, "status_note"));
            }
            if let Some(b) = c.get("blocked_on").and_then(Value::as_array) {
                for l in b {
                    println!("  blocked on: {}", l.as_str().unwrap_or("?"));
                }
            }
            println!("  question: {}", s(&c, "question"));
            println!("  decision value: {}", s(&c, "decision_value"));
            let code = c.get("code").cloned().unwrap_or(Value::Null);
            let pin = if s(&code, "commit").is_empty() {
                format!("branch {}", s(&code, "branch"))
            } else {
                format!("commit {}", s(&code, "commit"))
            };
            println!("  code: {} at {pin}", s(&code, "repo"));
            let r = c.get("resources").cloned().unwrap_or(Value::Null);
            println!(
                "  resources: {} cpus, {} GB RAM, {} GB disk, est {} core-h / {} wall-h",
                r.get("cpus").map_or("?".into(), |v| v.to_string()),
                r.get("ram_gb").map_or("?".into(), |v| v.to_string()),
                r.get("disk_gb").map_or("?".into(), |v| v.to_string()),
                r.get("est_core_hours")
                    .map_or("?".into(), |v| v.to_string()),
                r.get("est_wall_hours")
                    .map_or("?".into(), |v| v.to_string()),
            );
            println!(
                "  shards: {} ({}, stop {})",
                c.pointer("/shards/count")
                    .map_or("?".into(), |v| v.to_string()),
                s(&c, "shards.scheme"),
                s(&c, "shards.stop")
            );
            let argv: Vec<String> = c
                .pointer("/run/argv")
                .and_then(Value::as_array)
                .map(|a| {
                    a.iter()
                        .filter_map(|x| x.as_str().map(String::from))
                        .collect()
                })
                .unwrap_or_default();
            match card::shard_vars(&c, shard, "<out>").and_then(|v| card::substitute(&argv, &v)) {
                Ok(cmd) => println!("  shard {shard}: {}", cmd.join(" ")),
                Err(e) => println!("  shard {shard}: cannot resolve: {e}"),
            }
            println!(
                "  cairn: {} — {}",
                s(&c, "cairn.mode"),
                s(&c, "cairn.notes")
            );
            println!("  lands in: {}", s(&c, "outputs.land"));
            println!(
                "  narrative: {}",
                card_dir(&root, &id).join("README.md").display()
            );
            Ok(ExitCode::SUCCESS)
        }
        Cmd::New { slug, title, date } => {
            let date = match date {
                Some(d) => d,
                None => today()?,
            };
            if !card::valid_date(&date) {
                return Err(format!("{date:?} is not YYYY-MM-DD"));
            }
            let id = format!("XP-{}-{slug}", date.replace('-', ""));
            if !card::valid_id(&id) {
                return Err(format!("{id:?} is not a valid card id (slug: lowercase letters, digits, single hyphens)"));
            }
            let dir = card_dir(&root, &id);
            if dir.exists() {
                return Err(format!("{} already exists", dir.display()));
            }
            std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
            let title = title.unwrap_or_else(|| "TODO: a one-line title".to_string());
            let card = scaffold(&id, &date, &title);
            std::fs::write(dir.join("card.json"), pretty(&card)).map_err(|e| e.to_string())?;
            std::fs::write(dir.join("README.md"), readme_template(&id, &title))
                .map_err(|e| e.to_string())?;
            println!("{}", dir.display());
            Ok(ExitCode::SUCCESS)
        }
        Cmd::Validate {
            ids,
            all,
            no_git,
            merged_into,
        } => {
            let ids: Vec<String> = if all || ids.is_empty() {
                card::card_ids(&root).into_iter().collect()
            } else {
                ids
            };
            if ids.is_empty() {
                println!("no cards under {}", root.join("backlog").display());
            }
            let mut failed = 0;
            for id in &ids {
                let issues = match load(&root, id) {
                    Ok(c) => {
                        let mut issues = card::validate(&c, &ctx(&root, id, !no_git));
                        let executable = card::EXECUTABLE.contains(&s(&c, "status"));
                        if let (Some(r), true, "crypto") =
                            (&merged_into, executable, s(&c, "code.repo"))
                        {
                            let commit = s(&c, "code.commit");
                            let merged = std::process::Command::new("git")
                                .arg("-C")
                                .arg(&root)
                                .args(["merge-base", "--is-ancestor", commit, r])
                                .status()
                                .is_ok_and(|st| st.success());
                            if !merged {
                                issues.push(card::Issue {
                                    error: true,
                                    field: "code.commit".into(),
                                    msg: format!(
                                        "{commit} is not an ancestor of {r}: pin a merged commit"
                                    ),
                                });
                            }
                        }
                        issues
                    }
                    Err(e) => vec![card::Issue {
                        error: true,
                        field: "card.json".into(),
                        msg: e,
                    }],
                };
                let errors = issues.iter().filter(|i| i.error).count();
                println!("{id}: {}", if errors == 0 { "ok" } else { "FAILED" });
                for i in &issues {
                    println!("  {i}");
                }
                if errors > 0 {
                    failed += 1;
                }
            }
            if failed > 0 {
                println!("{failed} of {} cards failed", ids.len());
                Ok(ExitCode::from(1))
            } else {
                Ok(ExitCode::SUCCESS)
            }
        }
        Cmd::CairnPublish { id, lab, identity } => {
            let c = load(&root, &id)?;
            cairn::check_version()?;
            let issues = card::validate(&c, &ctx(&root, &id, true));
            if let Some(e) = issues.iter().find(|i| i.error) {
                return Err(format!(
                    "{id} does not validate ({e}); run `backlog validate {id}`"
                ));
            }
            let done = cairn::publish(&lab, &identity, &card_dir(&root, &id), &id, &c)?;
            for w in &done.written {
                println!("written   {w}");
            }
            for u in &done.unchanged {
                println!("unchanged {u}");
            }
            Ok(ExitCode::SUCCESS)
        }
        Cmd::CairnPost {
            id,
            node,
            identity,
            dry_run: _,
            confirm,
            allow_local_http,
        } => {
            let c = load(&root, &id)?;
            let issues = card::validate(&c, &ctx(&root, &id, true));
            if let Some(e) = issues.iter().find(|i| i.error) {
                return Err(format!(
                    "{id} does not validate ({e}); run `backlog validate {id}`"
                ));
            }
            let draft = cairn::objective_draft(&c, &now_rfc3339()?)?;
            if !confirm {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&draft).map_err(|e| e.to_string())?
                );
                eprintln!(
                    "dry run: nothing posted. With --confirm --node URL --identity FILE this draft is signed \
                     (`cairn sign-objective`) and sent to URL/submit?kind=objective, funded with {} units.",
                    card::funding(&c)
                );
                return Ok(ExitCode::SUCCESS);
            }
            let node = node.ok_or("--confirm needs --node URL")?;
            cairn::check_node(&node, allow_local_http)?;
            let identity = identity.ok_or("--confirm needs --identity FILE")?;
            cairn::check_version()?;
            let scratch = cairn::private_tempdir("post")?;
            let answer = cairn::post(&draft, &node, &identity, &scratch);
            let _ = std::fs::remove_dir_all(&scratch);
            println!("{}", answer?);
            Ok(ExitCode::SUCCESS)
        }
        Cmd::CairnFound {
            id,
            lab,
            lab_path,
            claim,
            holder,
            identity,
            checker_bin,
            allow_status,
        } => {
            let c = load(&root, &id)?;
            cairn::check_version()?;
            let doc_name = format!("backlog/{id}");
            match cairn::doc_field(&lab, &doc_name, "found")? {
                cairn::Field::One(f) => {
                    println!("found: {f}");
                    return Ok(ExitCode::from(3));
                }
                cairn::Field::Conflict(values) => {
                    println!("found (conflict: {} concurrent records):", values.len());
                    for v in values {
                        println!("  {v}");
                    }
                    return Ok(ExitCode::from(3));
                }
                cairn::Field::Unset => {}
            }
            let Some(lab_path) = lab_path else {
                println!("not found");
                return Ok(ExitCode::SUCCESS);
            };
            // Only a valid, runnable card that is the one published in the lab.
            let issues = card::validate(&c, &ctx(&root, &id, true));
            if let Some(e) = issues.iter().find(|i| i.error) {
                return Err(format!("{id} does not validate ({e})"));
            }
            let status = s(&c, "status").to_string();
            if !allow_status && status != "ready" && status != "running" {
                return Err(format!(
                    "{id} is {status}: only a ready or running card records a find"
                ));
            }
            let published = cairn::doc_get(&lab, &doc_name)?;
            if published.get("card") != Some(&c) {
                return Err(format!(
                    "the local card differs from doc {doc_name}.card in the lab: run `backlog cairn-publish {id}` from the merged card first"
                ));
            }
            if !card::under_runs(&id, &lab_path, false) {
                return Err(format!(
                    "--lab-path {lab_path:?} is not a path below runs/{id}/"
                ));
            }
            let (claim, holder) = (claim.unwrap_or_default(), holder.unwrap_or_default());
            let count = c
                .pointer("/shards/count")
                .and_then(Value::as_u64)
                .unwrap_or(1);
            let mut task = None;
            for t in cairn::task_keys(&lab, &id, count)?
                .into_iter()
                .chain((count == 1).then(|| id.clone()))
            {
                if cairn::live_holder(&lab, &t)? == Some((claim.clone(), holder.clone())) {
                    task = Some(t);
                    break;
                }
            }
            let Some(task) = task else {
                return Err(format!(
                    "{claim} is not the live lease of {holder:?} on any task of {id}: record the find before releasing"
                ));
            };
            // The bytes the lab holds, in a private temp file.
            let bytes = cairn::lab_bytes(&lab, &["cat", &lab_path])?;
            let tmp = cairn::private_tempdir("found")?;
            let file = tmp.join("artifact.json");
            let result = (|| -> Result<ExitCode, String> {
                std::fs::write(&file, &bytes).map_err(|e| e.to_string())?;
                let checker_bin = checker_bin.ok_or("--checker-bin is required")?;
                let bin_bytes = std::fs::read(&checker_bin)
                    .map_err(|e| format!("{}: {e}", checker_bin.display()))?;
                let mut argv: Vec<String> = c
                    .pointer("/cairn/checker/argv")
                    .and_then(Value::as_array)
                    .map(|a| {
                        a.iter()
                            .filter_map(|x| x.as_str().map(String::from))
                            .collect()
                    })
                    .unwrap_or_default();
                if argv.is_empty()
                    || c.pointer("/cairn/checker/status").and_then(Value::as_str)
                        != Some("available")
                {
                    return Err(format!("{id} has no available checker"));
                }
                // argv[0] names the checker; the binary run is the pinned build's.
                argv[0] = checker_bin.to_string_lossy().into_owned();
                let mut vars = std::collections::BTreeMap::new();
                vars.insert("artifact".to_string(), file.to_string_lossy().into_owned());
                let argv = card::substitute(&argv, &vars)?;
                let out = std::process::Command::new(&argv[0])
                    .args(&argv[1..])
                    .output()
                    .map_err(|e| format!("cannot run the checker {}: {e}", argv[0]))?;
                let line = String::from_utf8_lossy(&out.stdout).trim().to_string();
                let verdict: Value =
                    serde_json::from_str(&line).unwrap_or(Value::String(line.clone()));
                let accepted = out.status.success()
                    && verdict.get("verdict").and_then(Value::as_str) == Some("accept");
                if !accepted {
                    println!("checker did not accept: {line}");
                    return Ok(ExitCode::from(1));
                }
                if verdict.get("test_prime").and_then(Value::as_bool) == Some(true)
                    && status != "draft"
                {
                    return Err("a test-prime verdict on a non-draft card is refused".into());
                }
                let found = serde_json::json!({
                    "artifact": lab_path,
                    "artifact_sha256": hex::encode(crypto_lib::hash::sha256::sha256(&bytes)),
                    "verdict": verdict,
                    "claim": claim,
                    "holder": holder,
                    "task": task,
                    "checker_sha256": hex::encode(crypto_lib::hash::sha256::sha256(&bin_bytes)),
                });
                let identity = identity.ok_or("--identity is required")?;
                cairn::doc_set(&lab, &identity, &doc_name, "found", &found)?;
                println!("recorded found in doc {doc_name}: {found}");
                Ok(ExitCode::SUCCESS)
            })();
            let _ = std::fs::remove_dir_all(&tmp);
            result
        }
        Cmd::CairnStatus { id, lab } => {
            let c = load(&root, &id)?;
            cairn::check_version()?;
            let mut tasks: Vec<String> = c
                .pointer("/cairn/tasks")
                .and_then(Value::as_array)
                .map(|a| {
                    a.iter()
                        .filter_map(|x| x.as_str().map(String::from))
                        .collect()
                })
                .unwrap_or_default();
            if tasks.is_empty() {
                tasks.push(id.clone());
            }
            // `<ID>/shard-{k}` stands for every shard: list the ones the lab has seen.
            let pattern = format!("{id}/shard-{{k}}");
            let total = c
                .pointer("/shards/count")
                .and_then(Value::as_u64)
                .unwrap_or(1);
            if tasks == [pattern] {
                let seen = cairn::task_keys(&lab, &id, total)?;
                println!(
                    "{} of {total} shards have lease history in the lab; the rest are unclaimed",
                    seen.len()
                );
                tasks = seen;
            }
            match cairn::doc_field(&lab, &format!("backlog/{id}"), "found")? {
                cairn::Field::One(f) => println!("FOUND (stop: every runner stops claiming): {f}"),
                cairn::Field::Conflict(values) => {
                    println!(
                        "FOUND (conflict: {} concurrent records; every runner stops claiming):",
                        values.len()
                    );
                    for v in values {
                        println!("  {v}");
                    }
                }
                cairn::Field::Unset => {}
            }
            let rows = cairn::status(&lab, &tasks)?;
            println!(
                "{:<52} {:<10} {:<24} {:>4} {:>4} {:>4}  runs",
                "task", "state", "holder", "cont", "exp", "rel"
            );
            for r in rows {
                println!(
                    "{:<52} {:<10} {:<24} {:>4} {:>4} {:>4}  {}",
                    r.task,
                    r.state,
                    r.holder,
                    r.contended,
                    r.expired,
                    r.released,
                    if r.runs.is_empty() {
                        "-".to_string()
                    } else {
                        r.runs.join("; ")
                    }
                );
            }
            Ok(ExitCode::SUCCESS)
        }
    }
}

fn pretty(v: &Value) -> String {
    let mut s = serde_json::to_string_pretty(v).unwrap_or_default();
    s.push('\n');
    s
}

fn now_secs() -> Result<i64, String> {
    Ok(std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|e| e.to_string())?
        .as_secs() as i64)
}

/// Now, UTC, as `YYYY-MM-DDTHH:MM:SS+00:00` (cairn's `created_at` form),
/// from one reading of the clock.
fn now_rfc3339() -> Result<String, String> {
    let secs = now_secs()?;
    let t = secs.rem_euclid(86_400);
    Ok(format!(
        "{}T{:02}:{:02}:{:02}+00:00",
        civil_date(secs),
        t / 3600,
        t % 3600 / 60,
        t % 60
    ))
}

/// Today's UTC date from the system clock (no date crate in this workspace).
fn today() -> Result<String, String> {
    Ok(civil_date(now_secs()?))
}

/// The UTC civil date of a Unix time (H. Hinnant's civil-from-days).
fn civil_date(secs: i64) -> String {
    let z = secs.div_euclid(86_400) + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}-{m:02}-{d:02}")
}

fn scaffold(id: &str, date: &str, title: &str) -> Value {
    serde_json::json!({
        "schema": card::SCHEMA,
        "id": id,
        "title": title,
        "created": date,
        "authors": ["TODO: person, role or lane (no email)"],
        "status": "draft",
        "status_note": "",
        "blocked_on": [],
        "links": {"lead": "TODO: encyclopedia section or note", "notes": [], "prs": []},
        "question": "TODO: the question, in one sentence",
        "decision_value": "TODO: one paragraph — which result would change which decision, and what each outcome leads to next. A run whose outcome changes nothing is not worth its cost.",
        "protocol": {"path": null, "blob": null},
        "code": {"repo": "crypto", "commit": null, "branch": "TODO: branch, until a commit is pinned", "build": [["cargo", "build", "--release", "--bin", "TODO"]]},
        "resources": {"cpus": 1, "ram_gb": 1, "disk_gb": 1, "gpu": null, "est_core_hours": null, "est_wall_hours": null, "notes": "TODO: the cost derivation lives in README.md"},
        "run": {"argv": ["target/release/TODO", "--out", "{out}"], "env": {"RAYON_NUM_THREADS": "1"}, "timeout_per_shard_s": 3600},
        "shards": {"count": 1, "scheme": "single", "params": {}, "stop": "all"},
        "checkpoint": {"file": null, "resume_flag": null},
        "outputs": {"result_schema": "backlog.shard-result/v1", "verify": [], "land": "research/TODO"},
        "cairn": {"mode": "lab-task", "space": "TBD", "tasks": [id], "ttl_s": 86400,
                  "env": {"name": "crypto-rust", "recipe": "TODO: how to build the environment tree from code.commit"},
                  "exec": {"cpus": 1, "memory_mb": 1024, "timeout_s": 3600, "sandbox": "auto", "publish": format!("runs/{id}")},
                  "notes": "",
                  "funding": {"units": 0, "rationale": "TODO: proposed by the author in the PR, approved by the reviewer who merges it"}},
        "owner_hint": null
    })
}

fn readme_template(id: &str, title: &str) -> String {
    format!(
        "# {id}: {title}\n\n\
         **Status:** draft. See `card.json`; the protocol is `backlog/README.md`.\n\n\
         ## Question\n\nTODO\n\n\
         ## Why it matters (decision value)\n\nTODO: which result changes which decision.\n\n\
         ## Prerequisites\n\nTODO: code, data, hardware; link each blocker (PR, problems/*.md, another XP card).\n\n\
         ## Expected cost\n\nTODO: the derivation, with every input cited and every extrapolation marked as one.\n\n\
         ## What \"done\" means\n\nTODO: the stop condition, what lands in research/, and the independent check.\n"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn civil_dates_and_one_clock_reading() {
        assert_eq!(civil_date(0), "1970-01-01");
        assert_eq!(civil_date(951_782_400), "2000-02-29");
        assert_eq!(civil_date(1_792_077_870), "2026-10-15");
        let t = now_rfc3339().unwrap();
        assert!(
            card::valid_date(&t[..10]) && t.ends_with("+00:00") && t.len() == 25,
            "{t}"
        );
    }

    #[test]
    fn private_tempdirs_are_unpredictable_and_owner_only() {
        use std::os::unix::fs::PermissionsExt;
        let a = cairn::private_tempdir("t").unwrap();
        let b = cairn::private_tempdir("t").unwrap();
        assert_ne!(a, b);
        assert_eq!(
            std::fs::metadata(&a).unwrap().permissions().mode() & 0o777,
            0o700
        );
        let _ = std::fs::remove_dir_all(&a);
        let _ = std::fs::remove_dir_all(&b);
    }

    #[test]
    fn node_urls() {
        assert!(cairn::check_node("https://node.example", false).is_ok());
        assert!(cairn::check_node("http://127.0.0.1:8080", true).is_ok());
        assert!(cairn::check_node("http://127.0.0.1:8080", false).is_err());
        assert!(cairn::check_node("http://127.0.0.1.evil.example", true).is_err());
        assert!(cairn::check_node("http://localhost@evil.example", true).is_err());
        assert!(cairn::check_node("-Ohttps://x", false).is_err());
    }
}
