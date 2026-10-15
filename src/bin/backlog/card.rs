//! The experiment card (`backlog.card/v1`): loading, validation and the
//! run-template substitution a card's `run.argv` relies on.
//!
//! The rules live here and in `backlog/schema/card.v1.json`; a unit test pins
//! the schema's required keys and enums to the constants below so the two
//! cannot drift. `backlog/README.md` is the prose authority.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

pub const SCHEMA: &str = "backlog.card/v1";

/// Every top-level key a card may carry. Required unless listed in `OPTIONAL_KEYS`.
pub const KEYS: &[&str] = &[
    "schema",
    "id",
    "title",
    "created",
    "authors",
    "status",
    "status_note",
    "blocked_on",
    "links",
    "question",
    "decision_value",
    "protocol",
    "code",
    "resources",
    "run",
    "shards",
    "checkpoint",
    "outputs",
    "cairn",
    "owner_hint",
];
pub const OPTIONAL_KEYS: &[&str] = &["status_note", "owner_hint"];

pub const STATUSES: &[&str] = &[
    "draft",
    "ready",
    "blocked",
    "running",
    "done",
    "failed",
    "abandoned",
];
/// Statuses that claim the card can be (or was) executed as written: they need
/// a frozen protocol, a pinned commit and a verifier.
pub const EXECUTABLE: &[&str] = &["ready", "running", "done"];
pub const SCHEMES: &[&str] = &["seed_range", "list", "single"];
pub const STOPS: &[&str] = &["all", "first_success"];
pub const PLACEHOLDERS: &[&str] = &["shard", "seed_lo", "seed_hi", "out", "checkpoint"];
pub const CAIRN_MODES: &[&str] = &["lab-task", "objective"];
pub const SANDBOXES: &[&str] = &["auto", "none"];
/// Prefixes a `blocked_on` entry may use.
pub const LINK_KINDS: &[&str] = &["pr", "branch", "problem", "xp", "hw", "url", "note"];

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Issue {
    pub error: bool,
    pub field: String,
    pub msg: String,
}

impl std::fmt::Display for Issue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let level = if self.error { "error" } else { "warning" };
        write!(f, "{level}: {}: {}", self.field, self.msg)
    }
}

/// What validation may consult beyond the card itself.
pub struct Ctx {
    /// Repository root (holds `backlog/` and `problems/`).
    pub root: PathBuf,
    /// The card's directory.
    pub dir: PathBuf,
    /// Every card id under `backlog/`, for `xp:` links.
    pub known_ids: BTreeSet<String>,
    /// Ask git whether a protocol blob is in the object database.
    pub check_git: bool,
}

struct V<'a> {
    issues: Vec<Issue>,
    card: &'a Value,
}

impl<'a> V<'a> {
    fn err(&mut self, field: &str, msg: impl Into<String>) {
        self.issues.push(Issue {
            error: true,
            field: field.to_string(),
            msg: msg.into(),
        });
    }
    fn warn(&mut self, field: &str, msg: impl Into<String>) {
        self.issues.push(Issue {
            error: false,
            field: field.to_string(),
            msg: msg.into(),
        });
    }
    fn get(&self, path: &str) -> Option<&'a Value> {
        let mut at = self.card;
        for part in path.split('.') {
            at = at.get(part)?;
        }
        Some(at)
    }
    fn string(&mut self, path: &str) -> Option<&'a str> {
        match self.get(path) {
            Some(Value::String(s)) if !s.trim().is_empty() => Some(s.as_str()),
            Some(Value::String(_)) => {
                self.err(path, "must be a non-empty string");
                None
            }
            Some(_) => {
                self.err(path, "must be a string");
                None
            }
            None => {
                self.err(path, "is required");
                None
            }
        }
    }
    fn opt_string(&mut self, path: &str) -> Option<&'a str> {
        match self.get(path) {
            None | Some(Value::Null) => None,
            Some(Value::String(s)) => Some(s.as_str()),
            Some(_) => {
                self.err(path, "must be a string or null");
                None
            }
        }
    }
    fn object(&mut self, path: &str) -> Option<&'a serde_json::Map<String, Value>> {
        match self.get(path) {
            Some(Value::Object(m)) => Some(m),
            Some(_) => {
                self.err(path, "must be an object");
                None
            }
            None => {
                self.err(path, "is required");
                None
            }
        }
    }
    fn uint(&mut self, path: &str, min: u64) -> Option<u64> {
        match self.get(path).map(|v| (v, v.as_u64())) {
            Some((_, Some(n))) if n >= min => Some(n),
            Some((_, Some(_))) => {
                self.err(path, format!("must be an integer >= {min}"));
                None
            }
            Some(_) => {
                self.err(path, "must be a non-negative integer");
                None
            }
            None => {
                self.err(path, "is required");
                None
            }
        }
    }
    fn number(&mut self, path: &str, nullable: bool) -> Option<f64> {
        match self.get(path) {
            Some(Value::Null) if nullable => None,
            Some(v) if v.as_f64().is_some() => {
                let x = v.as_f64().unwrap_or(-1.0);
                if x < 0.0 {
                    self.err(path, "must be >= 0");
                }
                Some(x)
            }
            Some(_) => {
                self.err(
                    path,
                    if nullable {
                        "must be a number or null"
                    } else {
                        "must be a number"
                    },
                );
                None
            }
            None => {
                self.err(path, "is required");
                None
            }
        }
    }
    fn strings(&mut self, path: &str, nonempty: bool) -> Option<Vec<&'a str>> {
        match self.get(path) {
            Some(Value::Array(a)) => {
                let mut out = Vec::with_capacity(a.len());
                for (i, v) in a.iter().enumerate() {
                    match v.as_str() {
                        Some(s) => out.push(s),
                        None => self.err(&format!("{path}[{i}]"), "must be a string"),
                    }
                }
                if nonempty && out.is_empty() {
                    self.err(path, "must not be empty");
                }
                Some(out)
            }
            Some(_) => {
                self.err(path, "must be an array of strings");
                None
            }
            None => {
                self.err(path, "is required");
                None
            }
        }
    }
    fn one_of(&mut self, path: &str, allowed: &[&str]) -> Option<&'a str> {
        let s = self.string(path)?;
        if allowed.contains(&s) {
            Some(s)
        } else {
            self.err(path, format!("{s:?} is not one of {allowed:?}"));
            None
        }
    }
}

pub fn is_hex(s: &str, len: usize) -> bool {
    s.len() == len
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

/// `XP-<YYYYMMDD>-<slug>`, slug lowercase alphanumerics and single hyphens.
pub fn valid_id(id: &str) -> bool {
    let Some(rest) = id.strip_prefix("XP-") else {
        return false;
    };
    let (date, slug) = match rest.split_once('-') {
        Some(p) => p,
        None => return false,
    };
    date.len() == 8
        && date.bytes().all(|b| b.is_ascii_digit())
        && valid_date(&format!("{}-{}-{}", &date[..4], &date[4..6], &date[6..]))
        && !slug.is_empty()
        && slug.len() <= 64
        && !slug.starts_with('-')
        && !slug.ends_with('-')
        && !slug.contains("--")
        && slug
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
}

/// `YYYY-MM-DD` with a month in 1..=12 and a day in 1..=31.
pub fn valid_date(s: &str) -> bool {
    let b = s.as_bytes();
    if b.len() != 10 || b[4] != b'-' || b[7] != b'-' {
        return false;
    }
    let num = |r: std::ops::Range<usize>| s[r].parse::<u32>().ok();
    match (num(0..4), num(5..7), num(8..10)) {
        (Some(y), Some(m), Some(d)) => y >= 2000 && (1..=12).contains(&m) && (1..=31).contains(&d),
        _ => false,
    }
}

/// The id git gives a file's bytes: SHA-1 over `blob <len>\0<bytes>`.
pub fn git_blob_sha(bytes: &[u8]) -> String {
    let mut buf = format!("blob {}\0", bytes.len()).into_bytes();
    buf.extend_from_slice(bytes);
    hex::encode(crypto_lib::hash::sha1::sha1(&buf))
}

/// Every `{name}` token in a template string.
pub fn placeholders(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut rest = s;
    while let Some(start) = rest.find('{') {
        let after = &rest[start + 1..];
        match after.find('}') {
            Some(end) => {
                out.push(after[..end].to_string());
                rest = &after[end + 1..];
            }
            None => break,
        }
    }
    out
}

/// Replace every `{name}` in `argv` from `vars`. An unknown or unset name is an
/// error, never left in place: a literal `{seed_lo}` reaching a solver would
/// run the wrong search silently.
pub fn substitute(argv: &[String], vars: &BTreeMap<String, String>) -> Result<Vec<String>, String> {
    argv.iter()
        .map(|arg| {
            let mut out = String::with_capacity(arg.len());
            let mut rest = arg.as_str();
            while let Some(start) = rest.find('{') {
                out.push_str(&rest[..start]);
                let after = &rest[start + 1..];
                let end = after
                    .find('}')
                    .ok_or_else(|| format!("unterminated placeholder in {arg:?}"))?;
                let name = &after[..end];
                let value = vars
                    .get(name)
                    .ok_or_else(|| format!("placeholder {{{name}}} in {arg:?} has no value"))?;
                out.push_str(value);
                rest = &after[end + 1..];
            }
            out.push_str(rest);
            Ok(out)
        })
        .collect()
}

/// The placeholder values for shard `k` of a card, with `out` the output
/// directory a runner chose.
pub fn shard_vars(card: &Value, k: u64, out: &str) -> Result<BTreeMap<String, String>, String> {
    let shards = card.get("shards").ok_or("card has no shards block")?;
    let count = shards
        .get("count")
        .and_then(Value::as_u64)
        .ok_or("shards.count must be an integer")?;
    if k >= count {
        return Err(format!("shard {k} is out of range: the card has {count}"));
    }
    let mut vars = BTreeMap::new();
    vars.insert("shard".to_string(), k.to_string());
    vars.insert("out".to_string(), out.to_string());
    if let Some(file) = card
        .get("checkpoint")
        .and_then(|c| c.get("file"))
        .and_then(Value::as_str)
    {
        vars.insert("checkpoint".to_string(), format!("{out}/{file}"));
    }
    match shards.get("scheme").and_then(Value::as_str) {
        Some("seed_range") => {
            let params = shards.get("params").ok_or("shards.params is required")?;
            let lo = decimal(params.get("seed_lo"))
                .ok_or("shards.params.seed_lo must be a decimal string")?;
            let size = decimal(params.get("range_size"))
                .ok_or("shards.params.range_size must be a decimal string")?;
            let start = lo
                .checked_add(size.checked_mul(u128::from(k)).ok_or("seed overflow")?)
                .ok_or("seed overflow")?;
            let end = start.checked_add(size).ok_or("seed overflow")?;
            vars.insert("seed_lo".to_string(), start.to_string());
            vars.insert("seed_hi".to_string(), end.to_string());
        }
        Some("list") => {
            let items = shards
                .get("params")
                .and_then(|p| p.get("items"))
                .and_then(Value::as_array)
                .ok_or("shards.params.items must be an array")?;
            if let Some(item) = items
                .get(k as usize)
                .and_then(|i| i.get("vars"))
                .and_then(Value::as_object)
            {
                for (key, v) in item {
                    if let Some(s) = v.as_str() {
                        vars.insert(key.clone(), s.to_string());
                    }
                }
            }
        }
        _ => {}
    }
    Ok(vars)
}

/// A u128 written as a decimal string (JSON numbers cannot carry 2^64+ exactly
/// in every reader, so seed bounds are strings).
fn decimal(v: Option<&Value>) -> Option<u128> {
    let s = v?.as_str()?;
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    s.parse().ok()
}

/// Every card id under `<root>/backlog/`.
pub fn card_ids(root: &Path) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    if let Ok(entries) = fs::read_dir(root.join("backlog")) {
        for e in entries.flatten() {
            let name = e.file_name().to_string_lossy().to_string();
            if name.starts_with("XP-") && e.path().join("card.json").is_file() {
                out.insert(name);
            }
        }
    }
    out
}

pub fn load(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    serde_json::from_str(&text).map_err(|e| format!("{}: not JSON: {e}", path.display()))
}

/// Validate one card. An empty list, or warnings only, means it passes.
pub fn validate(card: &Value, ctx: &Ctx) -> Vec<Issue> {
    let mut v = V {
        issues: Vec::new(),
        card,
    };
    let Some(top) = card.as_object() else {
        v.err("card", "must be a JSON object");
        return v.issues;
    };
    for key in top.keys() {
        if !KEYS.contains(&key.as_str()) {
            v.err(key, "is not a card field (typo?)");
        }
    }
    for key in KEYS {
        if !OPTIONAL_KEYS.contains(key) && !top.contains_key(*key) {
            v.err(key, "is required");
        }
    }
    if top.len() < 4 {
        return v.issues;
    }

    if v.get("schema").and_then(Value::as_str) != Some(SCHEMA) {
        v.err("schema", format!("must be {SCHEMA:?}"));
    }

    // Identity.
    let id = v.string("id");
    if let Some(id) = id {
        if !valid_id(id) {
            v.err("id", format!("{id:?} is not XP-<YYYYMMDD>-<slug>"));
        }
        let dir_name = ctx.dir.file_name().map(|n| n.to_string_lossy().to_string());
        if dir_name.as_deref() != Some(id) {
            v.err("id", format!("does not match its directory {dir_name:?}"));
        }
    }
    if let Some(created) = v.string("created") {
        if !valid_date(created) {
            v.err("created", "must be YYYY-MM-DD");
        } else if let Some(id) = id {
            if id.get(3..11) != Some(&created.replace('-', "")) {
                v.err("created", "must be the date in the id");
            }
        }
    }
    v.string("title");
    v.string("question");
    if let Some(dv) = v.string("decision_value") {
        if dv.split_whitespace().count() < 25 {
            v.err(
                "decision_value",
                "must be a paragraph (at least 25 words): which result changes which decision",
            );
        }
    }
    if let Some(authors) = v.strings("authors", true) {
        for a in authors {
            if a.contains('@') {
                v.err("authors", "no email addresses: name a person, role or lane");
            }
        }
    }
    if let Some(note) = v.opt_string("owner_hint") {
        if note.contains('@') {
            v.err("owner_hint", "no email addresses");
        }
    }
    // Status and what blocks it.
    let status = v.one_of("status", STATUSES);
    v.opt_string("status_note");
    let executable = status.is_some_and(|s| EXECUTABLE.contains(&s));
    if let Some(links) = v.strings("blocked_on", false) {
        if status == Some("blocked") && links.is_empty() {
            v.err("blocked_on", "a blocked card names what blocks it");
        }
        if executable && !links.is_empty() {
            v.err(
                "blocked_on",
                format!("a {} card cannot be blocked", status.unwrap_or("")),
            );
        }
        for (i, link) in links.iter().enumerate() {
            check_link(&mut v, &format!("blocked_on[{i}]"), link, ctx);
        }
    }
    if v.object("links").is_some() {
        if let Some(Value::Array(prs)) = v.get("links.prs") {
            if prs.iter().any(|p| p.as_u64().is_none()) {
                v.err("links.prs", "must be PR numbers");
            }
        }
    }

    // Protocol: frozen (committed, blob recorded) before the card is executable.
    if v.object("protocol").is_some() {
        let path = v.opt_string("protocol.path");
        let blob = v.opt_string("protocol.blob");
        if let Some(b) = blob {
            if !is_hex(b, 40) {
                v.err("protocol.blob", "must be a 40-hex git blob id");
            }
        }
        if executable && (path.is_none() || blob.is_none()) {
            v.err(
                "protocol",
                "a ready/running/done card needs a frozen PROTOCOL.md (path and git blob)",
            );
        }
        if let Some(p) = path {
            let file = ctx.dir.join(p);
            match fs::read(&file) {
                Ok(bytes) => {
                    let actual = git_blob_sha(&bytes);
                    match blob {
                        Some(b) if b != actual => v.err(
                            "protocol.blob",
                            format!("{p} has changed since it was frozen: blob is {actual}, card says {b}"),
                        ),
                        Some(b) if ctx.check_git && !git_has_blob(&ctx.root, b) => v.err(
                            "protocol.blob",
                            format!("{b} is not in git's object database: commit {p} before marking the card {}", status.unwrap_or("?")),
                        ),
                        None => v.warn("protocol.blob", format!("not frozen yet (current blob {actual})")),
                        _ => {}
                    }
                }
                Err(_) => v.err(
                    "protocol.path",
                    format!("{} does not exist", file.display()),
                ),
            }
        }
    }

    // Code.
    if v.object("code").is_some() {
        v.string("code.repo");
        let commit = v.opt_string("code.commit");
        let branch = v.opt_string("code.branch");
        match (commit, branch) {
            (Some(c), _) if !is_hex(c, 40) => v.err("code.commit", "must be a full 40-hex commit"),
            (None, None) => v.err("code", "names a commit or a branch"),
            _ => {}
        }
        if executable && commit.is_none() {
            v.err(
                "code.commit",
                "a ready/running/done card pins a full commit, not a branch",
            );
        }
        match v.get("code.build") {
            Some(Value::Array(steps)) => {
                for (i, step) in steps.iter().enumerate() {
                    if !step
                        .as_array()
                        .is_some_and(|a| !a.is_empty() && a.iter().all(Value::is_string))
                    {
                        v.err(
                            &format!("code.build[{i}]"),
                            "must be a non-empty argv (array of strings)",
                        );
                    }
                }
            }
            _ => v.err("code.build", "must be an array of argv arrays"),
        }
    }

    // Resources.
    if v.object("resources").is_some() {
        v.uint("resources.cpus", 1);
        if v.number("resources.ram_gb", false) == Some(0.0) {
            v.err("resources.ram_gb", "must be > 0");
        }
        v.number("resources.disk_gb", false);
        match v.get("resources.gpu") {
            Some(Value::Null) | Some(Value::String(_)) => {}
            _ => v.err("resources.gpu", "must be a string (what GPU) or null"),
        }
        let core = v.number("resources.est_core_hours", true);
        v.number("resources.est_wall_hours", true);
        v.opt_string("resources.notes");
        if executable && core.is_none() {
            v.err(
                "resources.est_core_hours",
                "a ready card states its expected cost",
            );
        }
    }

    // Shards, run template, checkpoint.
    let count = if v.object("shards").is_some() {
        let count = v.uint("shards.count", 1);
        let scheme = v.one_of("shards.scheme", SCHEMES);
        v.one_of("shards.stop", STOPS);
        let params = v.object("shards.params");
        match (scheme, count, params) {
            (Some("single"), Some(c), _) if c != 1 => {
                v.err("shards.count", "a single-shard card has count 1")
            }
            (Some("seed_range"), Some(c), Some(_)) => {
                let lo = decimal(v.get("shards.params.seed_lo"));
                let size = decimal(v.get("shards.params.range_size"));
                match (lo, size) {
                    (Some(lo), Some(size)) if size > 0 => {
                        let end = size
                            .checked_mul(u128::from(c))
                            .and_then(|span| span.checked_add(lo));
                        if end.is_none() {
                            v.err("shards.params", "seed_lo + count * range_size overflows u128");
                        }
                        if let (Some(end), Some(hi)) = (end, decimal(v.get("shards.params.seed_hi"))) {
                            if end != hi {
                                v.err("shards.params.seed_hi", format!("must equal seed_lo + count * range_size = {end}"));
                            }
                        }
                    }
                    _ => v.err(
                        "shards.params",
                        "a seed_range card gives seed_lo and a positive range_size as decimal strings",
                    ),
                }
            }
            (Some("list"), Some(c), Some(_)) => match v.get("shards.params.items") {
                Some(Value::Array(items)) if items.len() as u64 == c => {}
                Some(Value::Array(items)) => v.err(
                    "shards.params.items",
                    format!("has {} items for {c} shards", items.len()),
                ),
                _ => v.err(
                    "shards.params.items",
                    "a list card gives one item per shard",
                ),
            },
            _ => {}
        }
        count
    } else {
        None
    };
    let checkpoint_file = if v.object("checkpoint").is_some() {
        let file = v.opt_string("checkpoint.file");
        let flag = v.opt_string("checkpoint.resume_flag");
        if executable && (file.is_none() || flag.is_none()) {
            v.err(
                "checkpoint",
                "a ready card's shard runs are resumable: give file and resume_flag",
            );
        }
        file
    } else {
        None
    };
    if v.object("run").is_some() {
        let argv = v.strings("run.argv", true);
        match v.get("run.env") {
            Some(Value::Object(m)) if m.values().all(Value::is_string) => {}
            _ => v.err("run.env", "must be an object of strings"),
        }
        v.uint("run.timeout_per_shard_s", 1);
        if let (Some(argv), Some(count)) = (argv, count) {
            let argv: Vec<String> = argv.iter().map(|s| s.to_string()).collect();
            for arg in &argv {
                for name in placeholders(arg) {
                    let list_var = v
                        .get("shards.params.items")
                        .and_then(Value::as_array)
                        .is_some_and(|items| {
                            items
                                .iter()
                                .all(|it| it.get("vars").and_then(|x| x.get(&name)).is_some())
                        });
                    if !PLACEHOLDERS.contains(&name.as_str()) && !list_var {
                        v.err("run.argv", format!("unknown placeholder {{{name}}}; allowed {PLACEHOLDERS:?} or a vars key every list item defines"));
                    }
                    if name == "checkpoint" && checkpoint_file.is_none() {
                        v.err("run.argv", "{checkpoint} needs checkpoint.file");
                    }
                }
            }
            // Resolvable at both ends of the shard range.
            for k in [0, count.saturating_sub(1)] {
                match shard_vars(card, k, "OUT").and_then(|vars| substitute(&argv, &vars)) {
                    Ok(_) => {}
                    Err(e) => v.err("run.argv", format!("shard {k}: {e}")),
                }
            }
        }
    }

    // Outputs.
    if v.object("outputs").is_some() {
        v.string("outputs.result_schema");
        let verify = v.strings("outputs.verify", false);
        if executable && verify.is_some_and(|a| a.is_empty()) {
            v.err(
                "outputs.verify",
                "a ready card names an independent check of any success claim",
            );
        }
        if let Some(land) = v.string("outputs.land") {
            if !land.starts_with("research/") {
                v.err("outputs.land", "results land under research/");
            }
        }
    }

    // cairn: how the card is executed (docs: backlog/README.md).
    check_cairn(&mut v, id, count, ctx);

    // A scaffold's TODOs may stay while the card is a draft, and only then.
    let mut todos = Vec::new();
    find_todos(card, "", &mut todos);
    for field in todos {
        if status == Some("draft") {
            v.warn(&field, "still TODO");
        } else {
            v.err(&field, "still TODO: fill it in before leaving draft");
        }
    }

    // The narrative.
    if !ctx.dir.join("README.md").is_file() {
        v.err("README.md", "every card directory has a README.md");
    }
    v.issues
}

fn find_todos(v: &Value, at: &str, out: &mut Vec<String>) {
    match v {
        Value::String(s) if s.contains("TODO") => out.push(if at.is_empty() {
            "card".to_string()
        } else {
            at.to_string()
        }),
        Value::Array(a) => a
            .iter()
            .enumerate()
            .for_each(|(i, x)| find_todos(x, &format!("{at}[{i}]"), out)),
        Value::Object(m) => m.iter().for_each(|(k, x)| {
            find_todos(
                x,
                &if at.is_empty() {
                    k.clone()
                } else {
                    format!("{at}.{k}")
                },
                out,
            )
        }),
        _ => {}
    }
}

fn check_link(v: &mut V, field: &str, link: &str, ctx: &Ctx) {
    let Some((kind, rest)) = link.split_once(':') else {
        v.err(
            field,
            format!("{link:?} is not <kind>:<value>, kind one of {LINK_KINDS:?}"),
        );
        return;
    };
    let rest = rest.trim();
    match kind {
        "pr" => {
            if rest
                .split_whitespace()
                .next()
                .and_then(|n| n.parse::<u64>().ok())
                .is_none()
            {
                v.err(field, "pr: takes a PR number, then optional text");
            }
        }
        "branch" => {
            if rest.split_whitespace().next().is_none() {
                v.err(field, "branch: takes a branch name");
            }
        }
        "problem" => {
            let path = rest.split_whitespace().next().unwrap_or("");
            if !path.starts_with("problems/") || !ctx.root.join(path).is_file() {
                v.err(field, format!("{path:?} is not a file under problems/"));
            }
        }
        "xp" => {
            let id = rest.split_whitespace().next().unwrap_or("");
            if !ctx.known_ids.contains(id) {
                v.err(field, format!("{id:?} is not a card under backlog/"));
            }
        }
        "url" => {
            if !rest.starts_with("https://") {
                v.err(field, "url: takes an https URL");
            }
        }
        "hw" | "note" => {
            if rest.is_empty() {
                v.err(field, format!("{kind}: needs text"));
            }
        }
        _ => v.err(
            field,
            format!("unknown link kind {kind:?}; use one of {LINK_KINDS:?}"),
        ),
    }
}

/// `cairn.checker`: the native checker of a success claim.
fn check_checker(v: &mut V, ctx: &Ctx) {
    let st = v.one_of("cairn.checker.status", &["available", "blocked"]);
    v.string("cairn.checker.bin");
    v.strings("cairn.checker.argv", true);
    v.object("cairn.checker.artifact_schema");
    if st == Some("blocked") {
        match v.get("cairn.checker.blocked_by").and_then(Value::as_str) {
            Some(link) => check_link(v, "cairn.checker.blocked_by", link, ctx),
            None => v.err(
                "cairn.checker.blocked_by",
                "a blocked checker names its blocker",
            ),
        }
    }
}

fn check_cairn(v: &mut V, id: Option<&str>, count: Option<u64>, ctx: &Ctx) {
    if v.object("cairn").is_none() {
        return;
    }
    let mode = v.one_of("cairn.mode", CAIRN_MODES);
    v.opt_string("cairn.notes");
    match mode {
        Some("lab-task") => {
            v.string("cairn.space");
            if let Some(tasks) = v.strings("cairn.tasks", true) {
                let pattern = id.map(|id| format!("{id}/shard-{{k}}"));
                if tasks.len() > 1 && tasks.iter().any(|t| Some(t.to_string()) == pattern) {
                    v.err(
                        "cairn.tasks",
                        "<id>/shard-{k} stands for every shard and is the only entry",
                    );
                }
                for (i, t) in tasks.iter().enumerate() {
                    if Some(t.to_string()) == pattern && tasks.len() == 1 {
                        continue;
                    }
                    let ok = match (id, t.split_once("/shard-")) {
                        (Some(id), None) => *t == id && count == Some(1),
                        (Some(id), Some((base, k))) => {
                            base == id
                                && k.parse::<u64>().ok().zip(count).is_some_and(|(k, c)| k < c)
                        }
                        _ => false,
                    };
                    if !ok {
                        v.err(
                            &format!("cairn.tasks[{i}]"),
                            format!("{t:?} must be the card id (one shard) or <id>/shard-<k> with k < shards.count"),
                        );
                    }
                }
            }
            v.uint("cairn.ttl_s", 60);
            if v.object("cairn.env").is_some() {
                v.string("cairn.env.name");
                v.string("cairn.env.recipe");
            }
            if v.object("cairn.exec").is_some() {
                v.uint("cairn.exec.cpus", 1);
                v.uint("cairn.exec.memory_mb", 1);
                v.uint("cairn.exec.timeout_s", 1);
                let sandbox = v.one_of("cairn.exec.sandbox", SANDBOXES);
                if let Some(publish) = v.string("cairn.exec.publish") {
                    if let Some(id) = id {
                        if !publish.starts_with(&format!("runs/{id}")) {
                            v.err("cairn.exec.publish", format!("must start with runs/{id}"));
                        }
                    }
                }
                if sandbox == Some("none") && v.opt_string("cairn.exec.sandbox_reason").is_none() {
                    v.err(
                        "cairn.exec.sandbox_reason",
                        "sandbox \"none\" says why (the receipt records an unconfined run)",
                    );
                }
            }
            // first_success on lab tasks: the finder runs the card's checker
            // and records `found` in the lab doc (backlog/README.md), so the
            // card names an available native checker taking {artifact}.
            let first = v.get("shards.stop").and_then(Value::as_str) == Some("first_success");
            if first || v.get("cairn.checker").is_some() {
                if v.object("cairn.checker").is_some() {
                    check_checker(v, ctx);
                }
                if first {
                    if v.get("cairn.checker.status").and_then(Value::as_str) != Some("available") {
                        v.err(
                            "cairn.checker.status",
                            "a first_success lab-task card needs an available checker: it decides `found`",
                        );
                    }
                    let takes_artifact = v
                        .get("cairn.checker.argv")
                        .and_then(Value::as_array)
                        .is_some_and(|a| {
                            a.iter()
                                .any(|s| s.as_str().is_some_and(|s| s.contains("{artifact}")))
                        });
                    if !takes_artifact {
                        v.err("cairn.checker.argv", "takes the artifact as {artifact}");
                    }
                }
            }
        }
        Some("objective") => {
            if v.one_of("cairn.kind", &["certificate", "workspace"])
                .is_some()
            {
                // The verifier: native only (AGENTS.md "no Python").
                if v.object("cairn.checker").is_some() {
                    check_checker(v, ctx);
                }
            }
            if v.object("cairn.piecework").is_some() {
                v.uint("cairn.piecework.units", 1);
                v.string("cairn.piecework.unit");
            }
            if let Some(stop) = v.one_of("cairn.stop", STOPS) {
                if v.get("shards.stop").and_then(Value::as_str) != Some(stop) {
                    v.err("cairn.stop", "must agree with shards.stop");
                }
            }
            // What `backlog cairn-post` sends: the cairn goal, the statement a
            // stranger reads, and the verifier block exactly as cairn takes it
            // (null until a native route exists).
            v.string("cairn.goal");
            v.opt_string("cairn.statement");
            match v.get("cairn.verifier") {
                Some(Value::Null) => {
                    if v.get("cairn.checker.status").and_then(Value::as_str) == Some("available") {
                        v.err(
                            "cairn.verifier",
                            "an available checker needs the cairn verifier block it posts with",
                        );
                    }
                }
                Some(Value::Object(m)) => {
                    match m.get("kind").and_then(Value::as_str) {
                        Some(k) if Some(k) == v.get("cairn.kind").and_then(Value::as_str) => {}
                        _ => v.err("cairn.verifier.kind", "must equal cairn.kind"),
                    }
                    // AGENTS.md "no Python": no pinned interpreter script.
                    let text = Value::Object(m.clone()).to_string();
                    if text.contains(".py\"") || text.contains("python") {
                        v.err(
                            "cairn.verifier",
                            "pins Python; this repository's verifiers are native (AGENTS.md)",
                        );
                    }
                }
                _ => v.err(
                    "cairn.verifier",
                    "must be the cairn verifier object, or null",
                ),
            }
        }
        _ => {}
    }
    // Funding is priority: it orders every runner's work and, for an
    // objective, is the reward `backlog cairn-post` funds it with.
    if v.object("cairn.funding").is_some() {
        v.uint("cairn.funding.units", 0);
        v.string("cairn.funding.rationale");
    }
}

/// A card's funding units (0 when absent): the order runners take work in.
pub fn funding(card: &Value) -> u64 {
    card.pointer("/cairn/funding/units")
        .and_then(Value::as_u64)
        .unwrap_or(0)
}

fn git_has_blob(root: &Path, blob: &str) -> bool {
    Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["cat-file", "-t", blob])
        .output()
        .map(|o| o.status.success() && String::from_utf8_lossy(&o.stdout).trim() == "blob")
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn fixture(root: &Path, id: &str, card: &Value) -> Ctx {
        let dir = root.join("backlog").join(id);
        fs::create_dir_all(&dir).unwrap();
        fs::write(
            dir.join("card.json"),
            serde_json::to_string_pretty(card).unwrap(),
        )
        .unwrap();
        fs::write(dir.join("README.md"), "# x\n").unwrap();
        Ctx {
            root: root.to_path_buf(),
            dir,
            known_ids: card_ids(root),
            check_git: false,
        }
    }

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("backlog-unit-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    pub fn base(id: &str) -> Value {
        json!({
            "schema": SCHEMA, "id": id, "title": "t", "created": "2026-10-15",
            "authors": ["lane backlog"], "status": "draft", "blocked_on": [],
            "links": {"prs": [1]},
            "question": "q?",
            "decision_value": "If the search finds a curve within the budget we build the factor base at P-192 size; if it does not, we stop the tower line at this size and record the bounded negative result in problems/ with its exact limits.",
            "protocol": {"path": null, "blob": null},
            "code": {"repo": "crypto", "commit": null, "branch": "claude/x", "build": [["cargo", "build", "--release"]]},
            "resources": {"cpus": 1, "ram_gb": 1, "disk_gb": 1, "gpu": null, "est_core_hours": 1.0, "est_wall_hours": null, "notes": ""},
            "run": {"argv": ["bin", "--from", "{seed_lo}", "--to", "{seed_hi}", "--out", "{out}/r.json", "--ckpt", "{checkpoint}"], "env": {}, "timeout_per_shard_s": 60},
            "shards": {"count": 4, "scheme": "seed_range", "params": {"seed_lo": "0", "range_size": "1000"}, "stop": "first_success"},
            "checkpoint": {"file": "ckpt.bin", "resume_flag": "--resume"},
            "outputs": {"result_schema": "backlog.shard-result/v1", "verify": ["check", "{out}"], "land": "research/x"},
            "cairn": {"mode": "objective", "kind": "workspace", "goal": "GOAL-x", "statement": "find one",
                      "checker": {"status": "available", "bin": "chk", "argv": ["chk", "artifact.json"], "artifact_schema": {}},
                      "verifier": {"kind": "workspace", "base": "b.json", "base_sha256": "00", "editable_paths": ["claim"],
                                   "score_command": ["cargo", "run", "--offline", "-q"], "score_path": "score.json"},
                      "piecework": {"units": 4, "unit": "a seed range"}, "stop": "first_success",
                      "funding": {"units": 400, "rationale": "test"}}
        })
    }

    fn errors(issues: &[Issue]) -> Vec<String> {
        issues
            .iter()
            .filter(|i| i.error)
            .map(|i| i.to_string())
            .collect()
    }

    #[test]
    fn a_complete_draft_passes() {
        let root = tmp("ok");
        let id = "XP-20261015-unit-ok";
        let c = base(id);
        let ctx = fixture(&root, id, &c);
        assert_eq!(errors(&validate(&c, &ctx)), Vec::<String>::new());
    }

    #[test]
    fn required_fields_typos_and_ids_are_errors() {
        let root = tmp("req");
        let id = "XP-20261015-unit-req";
        let mut c = base(id);
        c.as_object_mut().unwrap().remove("decision_value");
        c["titel"] = json!("x");
        c["created"] = json!("2026-10-14");
        let ctx = fixture(&root, id, &c);
        let e = errors(&validate(&c, &ctx)).join("\n");
        assert!(e.contains("decision_value: is required"), "{e}");
        assert!(e.contains("titel: is not a card field"), "{e}");
        assert!(e.contains("created: must be the date in the id"), "{e}");
        assert!(
            !valid_id("XP-2026101-x") && !valid_id("XP-20261015-Bad") && !valid_id("XP-20261315-x")
        );
        assert!(valid_id("XP-20261015-p192-isogeny-aux-search-t48"));
    }

    #[test]
    fn blocked_needs_links_and_links_are_typed() {
        let root = tmp("blk");
        let id = "XP-20261015-unit-blk";
        let mut c = base(id);
        c["status"] = json!("blocked");
        let ctx = fixture(&root, id, &c);
        assert!(errors(&validate(&c, &ctx))
            .iter()
            .any(|e| e.contains("names what blocks it")));
        c["blocked_on"] = json!([
            "pr:1889 open",
            "branch:claude/x",
            "problem:problems/nope.md",
            "xp:XP-20261015-nope",
            "bogus"
        ]);
        let e = errors(&validate(&c, &ctx)).join("\n");
        assert!(
            e.contains("blocked_on[2]")
                && e.contains("blocked_on[3]")
                && e.contains("blocked_on[4]"),
            "{e}"
        );
        assert!(
            !e.contains("blocked_on[0]") && !e.contains("blocked_on[1]"),
            "{e}"
        );
        c["blocked_on"] = json!(["xp:XP-20261015-unit-blk"]);
        assert!(errors(&validate(&c, &ctx)).is_empty());
    }

    #[test]
    fn ready_needs_a_frozen_protocol_a_commit_and_a_verifier() {
        let root = tmp("ready");
        let id = "XP-20261015-unit-ready";
        let mut c = base(id);
        c["status"] = json!("ready");
        let ctx = fixture(&root, id, &c);
        let e = errors(&validate(&c, &ctx)).join("\n");
        assert!(
            e.contains("frozen PROTOCOL.md") && e.contains("code.commit"),
            "{e}"
        );
        fs::write(ctx.dir.join("PROTOCOL.md"), "frozen\n").unwrap();
        c["protocol"] = json!({"path": "PROTOCOL.md", "blob": git_blob_sha(b"frozen\n")});
        c["code"]["commit"] = json!("0123456789abcdef0123456789abcdef01234567");
        assert_eq!(errors(&validate(&c, &ctx)), Vec::<String>::new());
        // Editing the frozen protocol is caught.
        fs::write(ctx.dir.join("PROTOCOL.md"), "edited\n").unwrap();
        assert!(errors(&validate(&c, &ctx))
            .iter()
            .any(|e| e.contains("has changed since it was frozen")));
        c["outputs"]["verify"] = json!([]);
        assert!(errors(&validate(&c, &ctx))
            .iter()
            .any(|e| e.contains("independent check")));
    }

    #[test]
    fn git_blob_ids_match_git() {
        // `printf 'hello\n' | git hash-object --stdin`
        assert_eq!(
            git_blob_sha(b"hello\n"),
            "ce013625030ba8dba906f756967f9e9ca394464a"
        );
        assert_eq!(
            git_blob_sha(b""),
            "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391"
        );
    }

    #[test]
    fn placeholders_substitute_or_fail_loudly() {
        let c = base("XP-20261015-unit-sub");
        let argv: Vec<String> = c["run"]["argv"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_string())
            .collect();
        let vars = shard_vars(&c, 3, "/o").unwrap();
        assert_eq!(
            substitute(&argv, &vars).unwrap(),
            vec![
                "bin",
                "--from",
                "3000",
                "--to",
                "4000",
                "--out",
                "/o/r.json",
                "--ckpt",
                "/o/ckpt.bin"
            ]
        );
        assert!(shard_vars(&c, 4, "/o").is_err());
        let bad = vec!["x{nope}".to_string()];
        assert!(substitute(&bad, &vars).unwrap_err().contains("{nope}"));
        assert!(substitute(&["x{shard".to_string()], &vars).is_err());
        // Seeds past 2^64 stay exact.
        let mut big = c.clone();
        big["shards"]["params"]["seed_lo"] = json!("18446744073709551616");
        let vars = shard_vars(&big, 1, "o").unwrap();
        assert_eq!(vars["seed_lo"], "18446744073709552616");
    }

    #[test]
    fn unknown_or_unresolvable_placeholders_are_errors() {
        let root = tmp("ph");
        let id = "XP-20261015-unit-ph";
        let mut c = base(id);
        c["run"]["argv"] = json!(["bin", "{seed}", "{checkpoint}"]);
        c["checkpoint"]["file"] = Value::Null;
        let ctx = fixture(&root, id, &c);
        let e = errors(&validate(&c, &ctx)).join("\n");
        assert!(e.contains("unknown placeholder {seed}"), "{e}");
        assert!(e.contains("{checkpoint} needs checkpoint.file"), "{e}");
        // {seed_lo} on a single-shard card has no value.
        let mut s = base(id);
        s["shards"] = json!({"count": 1, "scheme": "single", "params": {}, "stop": "all"});
        s["cairn"]["piecework"]["units"] = json!(1);
        s["cairn"]["stop"] = json!("all");
        let e = errors(&validate(&s, &ctx)).join("\n");
        assert!(e.contains("{seed_lo}") && e.contains("has no value"), "{e}");
    }

    #[test]
    fn cairn_lab_task_mode() {
        let root = tmp("lab");
        let id = "XP-20261015-unit-lab";
        let mut c = base(id);
        c["shards"] = json!({"count": 2, "scheme": "list", "params": {"items": [{"name": "a", "vars": {"m": "3"}}, {"name": "b", "vars": {"m": "4"}}]}, "stop": "all"});
        c["run"]["argv"] = json!(["bin", "--m", "{m}", "--out", "{out}"]);
        c["cairn"] = json!({"mode": "lab-task", "space": "TBD", "tasks": [format!("{id}/shard-0"), format!("{id}/shard-1")], "ttl_s": 3600,
            "env": {"name": "crypto-rust", "recipe": "cargo build --release at code.commit"},
            "exec": {"cpus": 8, "memory_mb": 262144, "timeout_s": 600000, "sandbox": "none", "sandbox_reason": "memory above the jail", "publish": format!("runs/{id}")},
            "funding": {"units": 1000, "rationale": "the open exponent question"}});
        let ctx = fixture(&root, id, &c);
        assert_eq!(errors(&validate(&c, &ctx)), Vec::<String>::new());
        c["cairn"]["tasks"] = json!([format!("{id}/shard-2"), "XP-20261015-other"]);
        c["cairn"]["exec"]
            .as_object_mut()
            .unwrap()
            .remove("sandbox_reason");
        c["cairn"]["exec"]["publish"] = json!("elsewhere");
        let e = errors(&validate(&c, &ctx)).join("\n");
        assert!(
            e.contains("cairn.tasks[0]") && e.contains("cairn.tasks[1]"),
            "{e}"
        );
        assert!(
            e.contains("sandbox_reason") && e.contains("cairn.exec.publish"),
            "{e}"
        );
    }

    #[test]
    fn first_success_lab_tasks_need_a_checker_and_may_use_the_shard_pattern() {
        let root = tmp("first");
        let id = "XP-20261015-unit-first";
        let mut c = base(id);
        c["cairn"] = json!({"mode": "lab-task", "space": "TBD", "tasks": [format!("{id}/shard-{{k}}")], "ttl_s": 21600,
            "env": {"name": "e", "recipe": "r"},
            "exec": {"cpus": 1, "memory_mb": 1024, "timeout_s": 60, "sandbox": "auto", "publish": format!("runs/{id}")},
            "checker": {"status": "available", "bin": "chk", "argv": ["chk", "{artifact}"], "artifact_schema": {}},
            "funding": {"units": 1, "rationale": "r"}});
        let ctx = fixture(&root, id, &c);
        assert_eq!(errors(&validate(&c, &ctx)), Vec::<String>::new());
        // No checker, or one that does not take the artifact: refused.
        let mut d = c.clone();
        d["cairn"]["checker"]["argv"] = json!(["chk"]);
        assert!(errors(&validate(&d, &ctx))
            .iter()
            .any(|e| e.contains("{artifact}")));
        d["cairn"].as_object_mut().unwrap().remove("checker");
        assert!(errors(&validate(&d, &ctx))
            .iter()
            .any(|e| e.contains("needs an available checker")));
        // The pattern is the only entry.
        let mut e = c.clone();
        e["cairn"]["tasks"] = json!([format!("{id}/shard-{{k}}"), format!("{id}/shard-0")]);
        assert!(errors(&validate(&e, &ctx))
            .iter()
            .any(|x| x.contains("only entry")));
    }

    #[test]
    fn cairn_objective_mode() {
        let root = tmp("obj");
        let id = "XP-20261015-unit-obj";
        let mut c = base(id);
        let ctx = fixture(&root, id, &c);
        c["cairn"]["funding"] = json!(1000);
        c["cairn"]["stop"] = json!("all");
        c["cairn"]["checker"]["status"] = json!("blocked");
        let e = errors(&validate(&c, &ctx)).join("\n");
        assert!(
            e.contains("cairn.funding") && e.contains("cairn.stop") && e.contains("blocked_by"),
            "{e}"
        );
        // An available checker needs the verifier block it would be posted with.
        let mut d = base(id);
        d["cairn"]["verifier"] = Value::Null;
        assert!(errors(&validate(&d, &ctx))
            .iter()
            .any(|e| e.contains("cairn.verifier")));
        // A pinned Python checker is refused.
        let mut p = base(id);
        p["cairn"]["kind"] = json!("certificate");
        p["cairn"]["verifier"] = json!({"kind": "certificate", "checker": "chk.py", "checker_sha256": "00", "entrypoint": "check"});
        assert!(errors(&validate(&p, &ctx))
            .iter()
            .any(|e| e.contains("pins Python")));
        assert_eq!(funding(&base(id)), 400);
        c["cairn"]["mode"] = json!("python");
        assert!(errors(&validate(&c, &ctx))
            .iter()
            .any(|e| e.contains("cairn.mode")));
    }

    #[test]
    fn the_published_schema_agrees_with_the_validator() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let schema = load(&root.join("backlog/schema/card.v1.json")).unwrap();
        let required: BTreeSet<&str> = schema["required"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        let ours: BTreeSet<&str> = KEYS
            .iter()
            .copied()
            .filter(|k| !OPTIONAL_KEYS.contains(k))
            .collect();
        assert_eq!(required, ours);
        let props: BTreeSet<&str> = schema["properties"]
            .as_object()
            .unwrap()
            .keys()
            .map(String::as_str)
            .collect();
        assert_eq!(props, KEYS.iter().copied().collect());
        let st: Vec<&str> = schema["properties"]["status"]["enum"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(st, STATUSES);
        let modes: Vec<&str> = schema["properties"]["cairn"]["properties"]["mode"]["enum"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(modes, CAIRN_MODES);
    }
}
