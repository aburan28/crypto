//! The pin (IC_TOOL_PROGRAM.md §5): a candidate's outputs on all 90 suite
//! rows against v0's from R01's profile pass, untimed.  This is the
//! rounds' declared `run.py pin` (R05's and R02b's are the same), native,
//! with the same `pin.json` byte for byte.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crypto_lib::cryptanalysis::curve_id;

use super::bench;
use super::json::{self, J};
use super::runs;
use super::suite::{self, Row};

/// Everything a pin compares: what the run decided, never how long it
/// took.
pub fn outputs(rep: &J) -> Vec<(&'static str, J)> {
    let get = |k: &str| rep.get(k).cloned().unwrap_or(J::Null);
    let scalar = |arm: &str| {
        rep.get("certificates")
            .and_then(|c| c.get(arm))
            .and_then(|c| c.get("scalar"))
            .cloned()
            .unwrap_or(J::Null)
    };
    vec![
        ("status", get("status")),
        ("counts", get("counts")),
        ("recovered_ic", scalar("ic")),
        ("recovered_rho", scalar("rho")),
        ("rho_counts", get("rho_counts")),
        ("all_verified", get("all_verified")),
        ("ic_and_rho_agree", get("ic_and_rho_agree")),
    ]
}

/// R01's run tree, extracted beside its archive on first use after its
/// SHA-256 is checked, as the declared runner did.
pub fn r01_runs(programme: &Path) -> Result<PathBuf, String> {
    let dir = programme.join("rounds").join("R01-baseline-v0");
    let runs = dir.join("runs");
    if runs.exists() {
        return Ok(runs);
    }
    let want = std::fs::read_to_string(dir.join("runs.tar.xz.sha256"))
        .map_err(|e| format!("R01's runs.tar.xz.sha256: {e}"))?;
    let want = want.split_whitespace().next().unwrap_or("");
    if bench::sha256_file(&dir.join("runs.tar.xz"))? != want {
        return Err("R01's runs.tar.xz does not match runs.tar.xz.sha256".into());
    }
    // Into a scratch directory first, so a half-extracted tree is never
    // taken for R01's.
    let tmp = dir.join(format!(".runs-extract-{}", std::process::id()));
    std::fs::create_dir_all(&tmp).map_err(|e| format!("{}: {e}", tmp.display()))?;
    let ok = Command::new("tar")
        .arg("-xJf")
        .arg(dir.join("runs.tar.xz"))
        .arg("-C")
        .arg(&tmp)
        .status()
        .map_err(|e| format!("tar: {e}"))?
        .success();
    if !ok {
        let _ = std::fs::remove_dir_all(&tmp);
        return Err("tar -xJf R01's runs.tar.xz failed".into());
    }
    let done = std::fs::rename(tmp.join("runs"), &runs);
    let _ = std::fs::remove_dir_all(&tmp);
    if done.is_err() && !runs.exists() {
        return Err(format!(
            "{}: could not move R01's runs into place",
            runs.display()
        ));
    }
    Ok(runs)
}

/// Counts only, under `taskset`: for pins, never for time.  An existing
/// output is read, not re-run.
fn untimed(binary: &Path, row: &Row, out: &Path) -> Result<J, String> {
    if !out.exists() {
        if let Some(parent) = out.parent() {
            std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
        }
        Command::new("taskset")
            .args(["-c", bench::CPUS])
            .args(bench::price_cmd(binary, row, out, &[])?)
            .env("RAYON_NUM_THREADS", "1")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .map_err(|e| format!("taskset: {e}"))?;
    }
    Ok(runs::load(out))
}

/// The candidate's pin against v0's outputs, written to `<runs>/pin/pin.json`
/// once; an existing record is returned as it is.
pub fn pin(programme: &Path, runs_dir: &Path, cand: &Path) -> Result<J, String> {
    let out = runs_dir.join("pin").join("pin.json");
    if out.exists() {
        return json::read(&out);
    }
    let r01 = r01_runs(programme)?;
    let mut result = Vec::new();
    let mut tiers = Vec::new();
    for tier in ["S", "smoke"] {
        tiers.extend(suite::rows(programme, tier)?.into_iter().map(|r| (tier, r)));
    }
    for (tier, row) in tiers {
        let suite_id = row.suite_id.clone().ok_or("a suite row has no suite id")?;
        let old = if tier == "S" {
            let first = r01
                .join("profile")
                .join("v0")
                .join(&suite_id)
                .join("r1.price.json");
            runs::load(&runs::figure_path(&first)?)
        } else {
            runs::load(&r01.join("smoke").join(format!("{suite_id}.price.json")))
        };
        let slug = suite::curve_slug(row.a, row.n)?;
        let v0_name = curve_id::resolve(old.get("curve").and_then(J::as_str).unwrap_or(""));
        let new = untimed(
            cand,
            &row,
            &runs_dir
                .join("pin")
                .join("candidate")
                .join(&row.id)
                .join("pin.price.json"),
        )?;
        let (a, b) = (outputs(&new), outputs(&old));
        let mut differs: Vec<&str> = a
            .iter()
            .zip(&b)
            .filter(|((_, x), (_, y))| !json::py_eq(x, y))
            .map(|((k, _), _)| *k)
            .collect();
        differs.sort_unstable();
        result.push(json::obj([
            ("row", J::Str(row.id.clone())),
            ("suite_id", J::Str(suite_id)),
            ("slug", J::Str(slug.clone())),
            ("equal", J::Bool(differs.is_empty())),
            (
                "differs_in",
                J::Arr(differs.iter().map(|k| J::Str((*k).into())).collect()),
            ),
            (
                "names_agree",
                J::Bool(
                    v0_name == Some(slug.as_str())
                        && new.get("curve").and_then(J::as_str) == Some(slug.as_str()),
                ),
            ),
        ]));
    }
    let all = |key: &str| {
        result
            .iter()
            .all(|e| e.get(key).is_some_and(|v| *v == J::Bool(true)))
    };
    let (held, names_agree) = (all("equal"), all("names_agree"));
    let doc = json::obj([
        ("rows", J::Arr(result)),
        ("held", J::Bool(held)),
        ("names_agree", J::Bool(names_agree)),
    ]);
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    std::fs::write(&out, json::dumps(&doc, 1) + "\n")
        .map_err(|e| format!("{}: {e}", out.display()))?;
    Ok(doc)
}

#[cfg(test)]
mod tests {
    use super::super::json::parse;
    use super::*;

    #[test]
    fn equality_is_pythons() {
        let eq = |a: &str, b: &str| json::py_eq(&parse(a).unwrap(), &parse(b).unwrap());
        assert!(eq(r#"{"a": 1, "b": [1, 2]}"#, r#"{"b": [1, 2.0], "a": 1}"#));
        assert!(eq("true", "1"));
        assert!(!eq(r#"{"a": 1}"#, r#"{"a": 1, "b": null}"#));
        assert!(!eq("[1, 2]", "[2, 1]"));
        assert!(!eq(r#""1""#, "1"));
        assert!(eq("null", "null"));
    }

    #[test]
    fn a_pin_reads_the_seven_outputs() {
        let rep =
            parse(r#"{"status": "complete", "certificates": {"ic": {"scalar": 5}}, "median": 1}"#)
                .unwrap();
        let keys: Vec<&str> = outputs(&rep).iter().map(|(k, _)| *k).collect();
        assert_eq!(
            keys,
            [
                "status",
                "counts",
                "recovered_ic",
                "recovered_rho",
                "rho_counts",
                "all_verified",
                "ic_and_rho_agree"
            ]
        );
        assert_eq!(outputs(&rep)[2].1, J::Int(5));
        assert_eq!(outputs(&rep)[3].1, J::Null);
    }
}
