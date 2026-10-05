//! The conformance runner (N4): `research/ic_tool_program/conformance/run.py`
//! and the case runner it loads, `v2/run.py`, natively.
//!
//! The cases are the steps' frozen `cases.json` files: `v1/` (B0's), `v2/`
//! (B1's), and every `v2-*/` in name order (the later steps').  Each set
//! with a `SHA256SUMS` (all but v1) is checked against it before any case
//! runs.  A case runs when its step is among the
//! steps named and its `until` step is not; a case another live case
//! `supersedes` does not run (design `schema-v2.md` §9, B3b's amendment 2).
//!
//! Every rule is the scripts': the files a case materialises (a copy,
//! optionally with dotted keys set, written as `json.dumps(doc, indent=1)`
//! writes it), the run with the case's environment and timeout, and the
//! expectations — the exit status (`zero`, `nonzero` or a number, and a
//! panic's 101 never passes), `stderr_contains`, and in the case's JSON
//! report `json_equals`, `json_equals_build_commit`, `json_paths`,
//! `json_contains` and `same_outputs_as`.  Values compare as Python's `==`
//! does.  Only the failure messages' wording differs from the scripts',
//! which printed Python's `repr` of each value; these print JSON.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use super::json::{self, obj, J};
use super::suite;

/// The steps, in the scripts' order (`conformance/run.py`'s `STEPS`).
pub const STEPS: [&str; 12] = [
    "B0", "B1", "B2", "B2b", "B3", "B3b", "B4", "B5", "B6", "B7", "B7a", "B7b",
];

/// A panicking Rust binary exits with this status.
const PANIC_STATUS: i32 = 101;

/// Each set's directory and its cases, as `run.py`'s `case_sets` reads
/// them: v1's (all step B0), v2's, then each `v2-*` with a `cases.json`.
fn case_sets(dir: &Path) -> Result<Vec<(PathBuf, Vec<J>)>, String> {
    let read = |d: &Path| -> Result<Vec<J>, String> {
        Ok(json::read(&d.join("cases.json"))?
            .at("cases")?
            .as_arr()
            .ok_or_else(|| format!("{}: `cases` is not a list", d.display()))?
            .to_vec())
    };
    let v1 = dir.join("v1");
    let mut sets = vec![(
        v1.clone(),
        read(&v1)?
            .into_iter()
            .map(|c| with_key(c, "step", J::Str("B0".into())))
            .collect(),
    )];
    sets.push((dir.join("v2"), read(&dir.join("v2"))?));
    let mut later: Vec<PathBuf> = std::fs::read_dir(dir)
        .map_err(|e| format!("{}: {e}", dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("v2-"))
                && p.join("cases.json").is_file()
        })
        .collect();
    later.sort();
    for d in later {
        let cases = read(&d)?;
        sets.push((d, cases));
    }
    Ok(sets)
}

/// `{**c, key: value}`: the key set in place, or appended.
fn with_key(case: J, key: &str, value: J) -> J {
    let J::Obj(mut kv) = case else { return case };
    match kv.iter_mut().find(|(k, _)| k == key) {
        Some((_, v)) => *v = value,
        None => kv.push((key.to_string(), value)),
    }
    J::Obj(kv)
}

fn str_of<'a>(case: &'a J, key: &str) -> Option<&'a str> {
    case.get(key).and_then(J::as_str)
}

/// The cases the named steps run, each with `{here}` made its own set's
/// `params/`, as `run.py`'s `cases_for` does it: in the case's text.
pub fn cases_for(dir: &Path, steps: &[String]) -> Result<Vec<J>, String> {
    let has = |s: Option<&str>| s.is_some_and(|s| steps.iter().any(|t| t == s));
    let live = |c: &J| has(str_of(c, "step")) && !has(str_of(c, "until"));
    let sets = case_sets(dir)?;
    let superseded: Vec<String> = sets
        .iter()
        .flat_map(|(_, cases)| cases.iter())
        .filter(|c| live(c))
        .filter_map(|c| str_of(c, "supersedes").map(str::to_string))
        .collect();
    let mut out = Vec::new();
    for (d, cases) in &sets {
        let here = d.join("params").to_string_lossy().into_owned();
        for c in cases {
            let id = str_of(c, "id").unwrap_or("");
            if live(c) && !superseded.iter().any(|s| s == id) {
                // `json.dumps(c).replace("{here}", …)`, then parsed again:
                // the path is inserted into the case's JSON text.
                let text =
                    json::dumps_compact(c, false, true).replace("{here}", &json_escape(&here));
                out.push(json::parse(&text)?);
            }
        }
    }
    Ok(out)
}

/// A path as it sits inside a JSON string (`json.dumps` escapes `\` and `"`).
fn json_escape(s: &str) -> String {
    let quoted = json::dumps(&J::Str(s.to_string()), 0);
    quoted[1..quoted.len() - 1].to_string()
}

/// What a case's strings may name.
struct Places {
    tmp: PathBuf,
    suite: PathBuf,
    cases: PathBuf,
}

impl Places {
    fn expand(&self, value: &str) -> String {
        value
            .replace("{tmp}", &self.tmp.to_string_lossy())
            .replace("{suite}", &self.suite.to_string_lossy())
            .replace("{cases}", &self.cases.to_string_lossy())
    }
}

/// `doc[a][b]…[last] = value`, through objects only, as `set_path` does.
fn set_path(doc: &mut J, dotted: &str, value: J) -> Result<(), String> {
    let keys: Vec<&str> = dotted.split('.').collect();
    let (last, parents) = keys.split_last().ok_or("an empty dotted path")?;
    let mut node = doc;
    for key in parents {
        let J::Obj(kv) = node else {
            return Err(format!("{dotted}: {key} is not inside an object"));
        };
        node = kv
            .iter_mut()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v)
            .ok_or_else(|| format!("{dotted}: no key {key:?}"))?;
    }
    let J::Obj(kv) = node else {
        return Err(format!("{dotted}: its parent is not an object"));
    };
    match kv.iter_mut().find(|(k, _)| k == last) {
        Some((_, v)) => *v = value,
        None => kv.push((last.to_string(), value)),
    }
    Ok(())
}

/// The value at a dotted path: list elements by index, object members by
/// key; `None` where `get_path` returns `MISSING`.
fn get_path<'a>(doc: &'a J, dotted: &str) -> Option<&'a J> {
    let mut node = doc;
    for key in dotted.split('.') {
        node = match node {
            J::Arr(items) if !key.is_empty() && key.bytes().all(|b| b.is_ascii_digit()) => {
                items.get(key.parse::<usize>().ok()?)?
            }
            J::Obj(_) => node.get(key)?,
            _ => return None,
        };
    }
    Some(node)
}

/// Whether `have` holds everything in `want`: objects key by key,
/// anything else by Python's equality.
fn contains(want: &J, have: &J) -> bool {
    match want {
        J::Obj(kv) => {
            matches!(have, J::Obj(_))
                && kv
                    .iter()
                    .all(|(k, v)| have.get(k).is_some_and(|h| contains(v, h)))
        }
        _ => json::py_eq(want, have),
    }
}

fn materialise(files: Option<&J>, at: &Places) -> Result<(), String> {
    for (rel, spec) in files.and_then(J::as_obj).unwrap_or(&[]) {
        let path = at.tmp.join(rel);
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
        }
        let text = if let Some(copy) = spec.get("copy").and_then(J::as_str) {
            let source = at.expand(copy);
            let text = std::fs::read_to_string(&source).map_err(|e| format!("{source}: {e}"))?;
            match spec.get("set").and_then(J::as_obj) {
                Some(sets) => {
                    let mut doc = json::parse(&text)?;
                    for (dotted, value) in sets {
                        set_path(&mut doc, dotted, value.clone())?;
                    }
                    json::dumps(&doc, 1) + "\n"
                }
                None => text,
            }
        } else {
            spec.get("text")
                .and_then(J::as_str)
                .ok_or_else(|| format!("file {rel}: neither `copy` nor `text`"))?
                .to_string()
        };
        std::fs::write(&path, text).map_err(|e| format!("{}: {e}", path.display()))?;
    }
    Ok(())
}

/// A finished run: its exit status (a signal is `-signal`, as Python
/// reports it) and its standard error.
struct Ran {
    status: i32,
    stderr: String,
}

/// `ic` with the case's arguments and environment, within its timeout;
/// `None` if it did not exit in time (it is then killed).
fn run(
    ic: &Path,
    argv: &[J],
    env: Option<&J>,
    at: &Places,
    timeout_s: f64,
) -> Result<Option<Ran>, String> {
    let mut cmd = Command::new(ic);
    for a in argv {
        cmd.arg(at.expand(a.as_str().ok_or("an argv entry is not a string")?));
    }
    for (k, v) in env.and_then(J::as_obj).unwrap_or(&[]) {
        cmd.env(k, v.as_str().ok_or("an env value is not a string")?);
    }
    let mut child = cmd
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("{}: {e}", ic.display()))?;
    // Both pipes drain on threads of their own, so a chatty child never
    // blocks on a full pipe while it is being waited for.
    let mut out = child.stdout.take().expect("piped stdout");
    let mut err = child.stderr.take().expect("piped stderr");
    let out_t = std::thread::spawn(move || {
        let mut sink = Vec::new();
        let _ = out.read_to_end(&mut sink);
    });
    let err_t = std::thread::spawn(move || {
        let mut text = Vec::new();
        let _ = err.read_to_end(&mut text);
        String::from_utf8_lossy(&text).into_owned()
    });
    let deadline = Instant::now() + Duration::from_secs_f64(timeout_s.max(0.0));
    let status = loop {
        if let Some(status) = child.try_wait().map_err(|e| e.to_string())? {
            break Some(status);
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            let _ = child.wait();
            break None;
        }
        std::thread::sleep(Duration::from_millis(5));
    };
    let _ = out_t.join();
    let stderr = err_t.join().unwrap_or_default();
    Ok(status.map(|s| Ran {
        status: s.code().unwrap_or_else(|| {
            use std::os::unix::process::ExitStatusExt;
            -s.signal().unwrap_or(0)
        }),
        stderr,
    }))
}

fn check_exit(want: Option<&J>, status: i32, failures: &mut Vec<String>) {
    if status == PANIC_STATUS {
        failures.push(format!("panicked (exit status {PANIC_STATUS})"));
    }
    match want {
        Some(J::Str(s)) if s == "zero" && status != 0 => {
            failures.push(format!("exit status {status}, expected 0"))
        }
        Some(J::Str(s)) if s == "nonzero" && status == 0 => {
            failures.push("exit status 0, expected a refusal".into())
        }
        Some(J::Int(n)) if *n != i128::from(status) => {
            failures.push(format!("exit status {status}, expected {n}"))
        }
        _ => {}
    }
}

fn load_report(path: &Path, failures: &mut Vec<String>) -> Option<J> {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    match std::fs::read_to_string(path)
        .map_err(|e| e.to_string())
        .and_then(|t| json::parse(&t))
    {
        Ok(doc) => Some(doc),
        Err(e) => {
            failures.push(format!("no JSON report at {name}: {e}"));
            None
        }
    }
}

/// A number as Python's `str` prints it, for a message.
fn num_text(v: &J) -> String {
    match v {
        J::Int(i) => i.to_string(),
        J::Float(x) => json::py_float(*x),
        other => json::dumps_line(other, false),
    }
}

fn shown(v: Option<&J>) -> String {
    v.map_or("None".into(), |v| json::dumps_compact(v, false, false))
}

fn check_report(doc: &J, expect: &J, build_commit: Option<&str>, failures: &mut Vec<String>) {
    for (key, want) in expect.get("json_equals").and_then(J::as_obj).unwrap_or(&[]) {
        // `doc.get(key) != want`: an absent key is `None`, which equals
        // only a `null` expectation.
        let have = doc.get(key);
        let equal = match have {
            Some(h) => json::py_eq(h, want),
            None => matches!(want, J::Null),
        };
        if !equal {
            failures.push(format!(
                "{key} = {}, expected {}",
                shown(have),
                shown(Some(want))
            ));
        }
    }
    if expect
        .get("json_equals_build_commit")
        .is_some_and(J::truthy)
    {
        if let Some(commit) = build_commit {
            if doc.get("build_commit").and_then(J::as_str) != Some(commit) {
                failures.push(format!(
                    "build_commit = {}, expected {}",
                    shown(doc.get("build_commit")),
                    shown(Some(&J::Str(commit.into())))
                ));
            }
        }
    }
    for (dotted, want) in expect.get("json_paths").and_then(J::as_obj).unwrap_or(&[]) {
        match get_path(doc, dotted) {
            None => failures.push(format!(
                "{dotted} is absent, expected {}",
                shown(Some(want))
            )),
            Some(have) if !json::py_eq(have, want) => failures.push(format!(
                "{dotted} = {}, expected {}",
                shown(Some(have)),
                shown(Some(want))
            )),
            _ => {}
        }
    }
    for (dotted, wants) in expect
        .get("json_contains")
        .and_then(J::as_obj)
        .unwrap_or(&[])
    {
        let Some(J::Arr(items)) = get_path(doc, dotted) else {
            failures.push(format!("{dotted} is not a list"));
            continue;
        };
        for want in wants.as_arr().unwrap_or(&[]) {
            if !items.iter().any(|item| contains(want, item)) {
                failures.push(format!(
                    "{dotted} has no element containing {}",
                    shown(Some(want))
                ));
            }
        }
    }
}

/// One case, as `v2/run.py`'s `run_case` runs it.
pub fn run_case(
    case: &J,
    ic: &Path,
    build_commit: Option<&str>,
    suite: &Path,
    cases_params: &Path,
) -> Result<J, String> {
    let id = str_of(case, "id").ok_or("a case has no id")?.to_string();
    let tmp = std::env::temp_dir().join(format!(
        "ic-conf-{id}-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0)
    ));
    std::fs::create_dir_all(&tmp).map_err(|e| format!("{}: {e}", tmp.display()))?;
    let at = Places {
        tmp: tmp.clone(),
        suite: suite.to_path_buf(),
        cases: cases_params.to_path_buf(),
    };
    let result = (|| -> Result<J, String> {
        let mut failures = Vec::new();
        materialise(case.get("files"), &at)?;
        let env = case.get("env");
        let timeout = case
            .at("timeout_s")?
            .as_f64()
            .ok_or("`timeout_s` is not a number")?;
        let argv = case.at("argv")?.as_arr().ok_or("`argv` is not a list")?;
        let Some(ran) = run(ic, argv, env, &at, timeout)? else {
            return Ok(obj([
                ("id", J::Str(id.clone())),
                ("pass", J::Bool(false)),
                (
                    "why",
                    J::Arr(vec![J::Str(format!(
                        "no exit within {} s",
                        num_text(case.at("timeout_s")?)
                    ))]),
                ),
            ]));
        };
        let expect = case.at("expect")?;
        check_exit(expect.get("exit"), ran.status, &mut failures);
        for needle in expect
            .get("stderr_contains")
            .and_then(J::as_arr)
            .unwrap_or(&[])
        {
            let needle = needle.as_str().unwrap_or("");
            if !ran.stderr.contains(needle) {
                failures.push(format!(
                    "stderr lacks {}",
                    shown(Some(&J::Str(needle.into())))
                ));
            }
        }
        let mut doc = None;
        if let Some(file) = expect.get("json_file").and_then(J::as_str) {
            doc = load_report(Path::new(&at.expand(file)), &mut failures);
            if let Some(d) = &doc {
                check_report(d, expect, build_commit, &mut failures);
            }
        }
        if let (Some(same), Some(d)) = (expect.get("same_outputs_as"), &doc) {
            let argv = same.at("argv")?.as_arr().ok_or("`argv` is not a list")?;
            match run(ic, argv, env, &at, timeout)? {
                None => failures.push(format!(
                    "the comparison run did not exit within {} s",
                    num_text(case.at("timeout_s")?)
                )),
                Some(other) => {
                    check_exit(Some(&J::Str("zero".into())), other.status, &mut failures);
                    let file = same.at("json_file")?.as_str().ok_or("`json_file`")?;
                    if let Some(reference) = load_report(Path::new(&at.expand(file)), &mut failures)
                    {
                        for dotted in same.at("paths")?.as_arr().unwrap_or(&[]) {
                            let dotted = dotted.as_str().unwrap_or("");
                            let (a, b) = (get_path(d, dotted), get_path(&reference, dotted));
                            let differ = match (a, b) {
                                (Some(a), Some(b)) => !json::py_eq(a, b),
                                (Some(_), None) => true,
                                (None, _) => true,
                            };
                            if differ {
                                failures
                                    .push(format!("{dotted} differs from the comparison run's"));
                            }
                        }
                    }
                }
            }
        }
        let tail: Vec<J> = {
            let lines: Vec<&str> = ran.stderr.trim().lines().collect();
            lines[lines.len().saturating_sub(3)..]
                .iter()
                .map(|l| J::Str((*l).to_string()))
                .collect()
        };
        Ok(J::Obj(vec![
            ("id".into(), J::Str(id.clone())),
            (
                "step".into(),
                J::Str(str_of(case, "step").unwrap_or("B0").to_string()),
            ),
            ("pass".into(), J::Bool(failures.is_empty())),
            (
                "why".into(),
                J::Arr(failures.into_iter().map(J::Str).collect()),
            ),
            ("exit".into(), J::Int(ran.status.into())),
            ("stderr_tail".into(), J::Arr(tail)),
        ]))
    })();
    let _ = std::fs::remove_dir_all(&tmp);
    result
}

/// `conformance/run.py --ic <ic> --steps <steps>`: every live case, and the
/// report (`binary`, `steps` in the steps' order, `cases`, `passed`,
/// `results`).  The second value is whether every case passed.
pub fn run_steps(
    programme: &Path,
    ic: &Path,
    steps: &[String],
    build_commit: Option<&str>,
) -> Result<(J, bool), String> {
    if let Some(bad) = steps.iter().find(|s| !STEPS.contains(&s.as_str())) {
        return Err(format!(
            "unknown step {bad:?}; the steps are {}",
            STEPS.join(", ")
        ));
    }
    let dir = programme.join("conformance");
    let suite = programme.join("suite").join("v1");
    let cases_params = dir.join("v2").join("params");
    // Each set with a `SHA256SUMS` is checked against it before any case
    // runs; v1, which has none, is pinned by the repository's history.
    for (set, _) in case_sets(&dir)? {
        if set.join("SHA256SUMS").is_file() {
            suite::check_sums(&set)?;
        }
    }
    let mut results = Vec::new();
    for case in cases_for(&dir, steps)? {
        results.push(run_case(&case, ic, build_commit, &suite, &cases_params)?);
    }
    let passed = results
        .iter()
        .filter(|r| r.get("pass").is_some_and(J::truthy))
        .count();
    let mut ordered: Vec<&str> = STEPS
        .iter()
        .copied()
        .filter(|s| steps.iter().any(|t| t == s))
        .collect();
    ordered.dedup();
    let all = passed == results.len();
    let report = obj([
        ("binary", J::Str(ic.to_string_lossy().into_owned())),
        (
            "steps",
            J::Arr(ordered.into_iter().map(|s| J::Str(s.into())).collect()),
        ),
        ("cases", J::Int(results.len() as i128)),
        ("passed", J::Int(passed as i128)),
        ("results", J::Arr(results)),
    ]);
    Ok((report, all))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc(text: &str) -> J {
        json::parse(text).unwrap()
    }

    #[test]
    fn paths_read_lists_by_index_and_objects_by_key() {
        let d = doc(r#"{"a": [{"b": 1}, {"b": 2}], "c": {"d": true}}"#);
        assert!(json::py_eq(get_path(&d, "a.1.b").unwrap(), &J::Int(2)));
        assert!(get_path(&d, "a.2.b").is_none());
        assert!(get_path(&d, "c.e").is_none());
        assert!(get_path(&d, "a.x").is_none());
        // Python's equality: `true == 1`.
        assert!(json::py_eq(get_path(&d, "c.d").unwrap(), &J::Int(1)));
    }

    #[test]
    fn set_path_replaces_in_place_and_appends_a_new_key() {
        let mut d = doc(r#"{"a": {"b": 1, "c": 2}}"#);
        set_path(&mut d, "a.b", J::Int(5)).unwrap();
        set_path(&mut d, "a.z", J::Int(9)).unwrap();
        assert_eq!(
            json::dumps_line(&d, false),
            r#"{"a": {"b": 5, "c": 2, "z": 9}}"#
        );
        assert!(set_path(&mut d, "x.y", J::Int(1)).is_err());
    }

    #[test]
    fn contains_is_a_subset_for_objects_and_equality_otherwise() {
        let have = doc(r#"{"code": "x", "detail": {"n": 3, "m": 4}}"#);
        assert!(contains(&doc(r#"{"code": "x"}"#), &have));
        assert!(contains(&doc(r#"{"detail": {"n": 3.0}}"#), &have));
        assert!(!contains(&doc(r#"{"detail": {"n": 4}}"#), &have));
        assert!(!contains(&doc(r#"{"other": 1}"#), &have));
        assert!(contains(&J::Int(3), &J::Float(3.0)));
    }

    #[test]
    fn the_exit_rule_never_passes_a_panic() {
        let mut f = Vec::new();
        check_exit(Some(&J::Str("nonzero".into())), PANIC_STATUS, &mut f);
        assert_eq!(f.len(), 1);
        f.clear();
        check_exit(Some(&J::Int(2)), 2, &mut f);
        check_exit(Some(&J::Str("zero".into())), 0, &mut f);
        assert!(f.is_empty());
        check_exit(Some(&J::Str("nonzero".into())), 0, &mut f);
        assert_eq!(f, ["exit status 0, expected a refusal"]);
    }

    #[test]
    fn json_equals_treats_an_absent_key_as_none() {
        let d = doc(r#"{"status": "complete"}"#);
        let mut f = Vec::new();
        check_report(
            &d,
            &doc(r#"{"json_equals": {"absent": null}}"#),
            None,
            &mut f,
        );
        assert!(f.is_empty());
        check_report(
            &d,
            &doc(r#"{"json_equals": {"status": "refused"}}"#),
            None,
            &mut f,
        );
        assert_eq!(f.len(), 1);
    }
}
