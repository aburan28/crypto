//! A round's run tree, as the programme's runner writes it
//! (IC_TOOL_PROGRAM.md §5).
//!
//! Every timed process is one `ic price --single-target` on one row.  Its
//! report is `<arm>/<row id>/r<k>.price.json`, its isolation record the
//! same stem with `.isolation.jsonl`.  A contended or failed process is
//! run again at most twice, as `r<k>-retry1` and `r<k>-retry2`; the first
//! clean, complete attempt is the row's figure, and every attempt stays.

use std::path::{Path, PathBuf};

use super::json::{self, J};
use super::stats;

pub const RETRIES: usize = 2;

/// `out` without `.price.json` (or `.json`).
pub fn stem(out: &Path) -> PathBuf {
    let name = out.file_name().and_then(|n| n.to_str()).unwrap_or("");
    for suffix in [".price.json", ".json"] {
        if let Some(s) = name.strip_suffix(suffix) {
            return out.with_file_name(s);
        }
    }
    out.to_path_buf()
}

pub fn record_path(out: &Path) -> PathBuf {
    let s = stem(out);
    let name = format!(
        "{}.isolation.jsonl",
        s.file_name().unwrap().to_string_lossy()
    );
    s.with_file_name(name)
}

/// The isolation tool's record of the attempt: its last line's `run`.
pub fn run_record(out: &Path) -> Result<Option<J>, String> {
    let rec = record_path(out);
    if !rec.exists() {
        return Ok(None);
    }
    let text = std::fs::read_to_string(&rec).map_err(|e| format!("{}: {e}", rec.display()))?;
    let last = text.lines().last().unwrap_or("");
    let line = json::parse(last).map_err(|e| format!("{}: {e}", rec.display()))?;
    Ok(line.get("run").cloned())
}

/// Ran to completion, with an exit status of 0, uncontended.
pub fn clean(out: &Path) -> Result<bool, String> {
    Ok(match run_record(out)? {
        Some(run) => {
            run.get("exit_status").and_then(J::as_i128) == Some(0)
                && !run.get("contended").is_some_and(J::truthy)
        }
        None => false,
    })
}

/// The report, or `{"status": "no report"}` when it is missing or unreadable.
pub fn load(out: &Path) -> J {
    json::read(out).unwrap_or_else(|_| J::Obj(vec![("status".into(), J::Str("no report".into()))]))
}

fn complete(rep: &J) -> bool {
    rep.get("status").and_then(J::as_str) == Some("complete")
}

pub fn attempt(out: &Path, k: usize) -> PathBuf {
    if k == 0 {
        out.to_path_buf()
    } else {
        let s = stem(out);
        let name = format!(
            "{}-retry{k}.price.json",
            s.file_name().unwrap().to_string_lossy()
        );
        s.with_file_name(name)
    }
}

/// The attempt that is the row's figure: the first clean, complete one,
/// or `out` itself when there is none.
pub fn figure_path(out: &Path) -> Result<PathBuf, String> {
    for k in 0..=RETRIES {
        let a = attempt(out, k);
        if a.exists() && clean(&a)? && complete(&load(&a)) {
            return Ok(a);
        }
    }
    Ok(out.to_path_buf())
}

/// The row's figure, or `None` when no attempt is clean and complete.
pub fn figure(out: &Path) -> Result<Option<J>, String> {
    let p = figure_path(out)?;
    if !p.exists() {
        return Ok(None);
    }
    let rep = json::read(&p)?;
    Ok((complete(&rep) && clean(&p)?).then_some(rep))
}

/// The rounds `arm` ran for `row_id`, from its first attempts' names.
pub fn rounds_of(d: &Path, arm: &str, row_id: &str) -> Vec<u32> {
    let dir = d.join(arm).join(row_id);
    let mut out: Vec<u32> = std::fs::read_dir(&dir)
        .into_iter()
        .flatten()
        .flatten()
        .filter_map(|e| {
            let name = e.file_name().to_string_lossy().into_owned();
            if !name.starts_with('r') || !name.ends_with(".price.json") || name.contains("-retry") {
                return None;
            }
            name[1..].split('.').next()?.parse().ok()
        })
        .collect();
    out.sort_unstable();
    out
}

/// Every isolation record under `d`, in a stable order.
fn records(d: &Path, out: &mut Vec<PathBuf>) {
    let mut entries: Vec<PathBuf> = std::fs::read_dir(d)
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| e.path())
        .collect();
    entries.sort();
    for p in entries {
        if p.is_dir() {
            records(&p, out);
        } else if p.to_string_lossy().ends_with(".isolation.jsonl") {
            out.push(p);
        }
    }
}

/// Processes, contended runs and failures under `d`, from every record.
pub fn accounting(d: &Path) -> Result<J, String> {
    let mut paths = Vec::new();
    records(d, &mut paths);
    let (mut processes, mut contended, mut failed) = (0i128, 0i128, 0i128);
    for rec in paths {
        let text = std::fs::read_to_string(&rec).map_err(|e| format!("{}: {e}", rec.display()))?;
        for line in text.lines() {
            let v = json::parse(line).map_err(|e| format!("{}: {e}", rec.display()))?;
            let run = v.at("run")?;
            processes += 1;
            contended += i128::from(run.at("contended")?.truthy());
            failed += i128::from(run.at("exit_status")?.as_i128() != Some(0));
        }
    }
    Ok(json::obj([
        ("processes", J::Int(processes)),
        ("contended", J::Int(contended)),
        ("failed", J::Int(failed)),
    ]))
}

// ── the measures a round compares (IC_TOOL_PROGRAM.md §5) ───────────

fn reps(rep: &J) -> Result<&[J], String> {
    rep.at("repetitions")?
        .as_arr()
        .ok_or_else(|| "`repetitions` is not a list".into())
}

fn num(v: &J, path: &[&str]) -> Result<f64, String> {
    let mut cur = v;
    for k in path {
        cur = cur.at(k)?;
    }
    cur.as_f64()
        .ok_or_else(|| format!("`{}` is not a number", path.join(".")))
}

fn median_of(rep: &J, f: impl Fn(&J) -> Result<f64, String>) -> Result<f64, String> {
    let xs = reps(rep)?.iter().map(f).collect::<Result<Vec<_>, _>>()?;
    Ok(stats::median(&xs))
}

/// The index calculus's cold time: set-up plus the online interval, per
/// repetition, at the median repetition.
pub fn ic_cold_ns(rep: &J) -> Result<f64, String> {
    median_of(rep, |r| {
        Ok(num(r, &["setup_ns"])? + num(r, &["ic_online", "wall_ns"])?)
    })
}

pub fn ic_online_ns(rep: &J) -> Result<f64, String> {
    median_of(rep, |r| num(r, &["ic_online", "wall_ns"]))
}

pub fn rho_online_ns(rep: &J) -> Result<f64, String> {
    median_of(rep, |r| num(r, &["rho_online", "wall_ns"]))
}

/// The median repetition's `key`, e.g. `setup_ns`.
pub fn median_rep(rep: &J, key: &str) -> Result<f64, String> {
    median_of(rep, |r| num(r, &[key]))
}

/// The median over repetitions of one set-up phase.
pub fn setup_phase_ns(rep: &J, phase: &str) -> Result<f64, String> {
    median_of(rep, |r| num(r, &["setup_phases_ns", phase]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn attempts_and_records_are_named_from_the_stem() {
        let out = Path::new("/x/base/row/r3.price.json");
        assert_eq!(
            record_path(out),
            Path::new("/x/base/row/r3.isolation.jsonl")
        );
        assert_eq!(attempt(out, 0), out);
        assert_eq!(
            attempt(out, 2),
            Path::new("/x/base/row/r3-retry2.price.json")
        );
        assert_eq!(
            record_path(&attempt(out, 1)),
            Path::new("/x/base/row/r3-retry1.isolation.jsonl")
        );
    }

    #[test]
    fn the_figure_is_the_first_clean_complete_attempt() {
        let dir = std::env::temp_dir().join(format!("icprog-runs-{}", std::process::id()));
        let row = dir.join("base").join("row");
        std::fs::create_dir_all(&row).unwrap();
        let write = |name: &str, text: &str| std::fs::write(row.join(name), text).unwrap();
        let rec = |contended: bool, exit: i32| {
            format!("{{\"run\": {{\"contended\": {contended}, \"exit_status\": {exit}}}}}\n")
        };
        let rep = |s: f64| {
            format!(
                "{{\"status\": \"complete\", \"repetitions\": [{{\"setup_ns\": {s}, \"ic_online\": {{\"wall_ns\": 10}}}}]}}"
            )
        };
        write("r1.price.json", &rep(100.0));
        write("r1.isolation.jsonl", &rec(true, 0));
        write("r1-retry1.price.json", &rep(200.0));
        write("r1-retry1.isolation.jsonl", &rec(false, 0));
        write("r2.price.json", &rep(300.0));
        write("r2.isolation.jsonl", &rec(false, 1));
        let r1 = figure(&row.join("r1.price.json")).unwrap().unwrap();
        assert_eq!(ic_cold_ns(&r1).unwrap(), 210.0);
        assert!(figure(&row.join("r2.price.json")).unwrap().is_none());
        assert_eq!(rounds_of(&dir, "base", "row"), vec![1, 2]);
        let acc = accounting(&dir).unwrap();
        assert_eq!(acc.get("processes"), Some(&J::Int(3)));
        assert_eq!(acc.get("contended"), Some(&J::Int(1)));
        assert_eq!(acc.get("failed"), Some(&J::Int(1)));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
