//! Suite v1's rows (IC_TOOL_PROGRAM.md §4), and the names a run tree
//! files them under (AGENTS.md §11).

use std::path::{Path, PathBuf};

use crypto_lib::cryptanalysis::curve_id;

use super::json::{self, J};

/// One row of a comparison: a suite row, or a round's own holdout.
#[derive(Clone, Debug)]
pub struct Row {
    /// `<slug>/<recipe>-<target>`, the run tree's directory.
    pub id: String,
    pub a: u8,
    pub n: u32,
    /// The subgroup order.
    pub r: f64,
    pub recipe_seed: Option<i128>,
}

/// The programme's directory, `research/ic_tool_program`, under `root`.
pub fn programme(root: &Path) -> Result<PathBuf, String> {
    let p = root.join("research").join("ic_tool_program");
    if p.is_dir() {
        Ok(p)
    } else {
        Err(format!(
            "{} is not a checkout of this repository (no research/ic_tool_program); pass --root",
            root.display()
        ))
    }
}

/// The ICV1 slug of the suite's Koblitz curve with coefficient `a` over
/// `GF(2^n)`, from the registry.
pub fn curve_slug(a: u8, n: u32) -> Result<String, String> {
    curve_id::resolve(&format!("K_{a} / GF(2^{n})"))
        .map(str::to_string)
        .ok_or_else(|| format!("the registry has no Koblitz curve a={a} over GF(2^{n})"))
}

fn int(row: &J, key: &str) -> Result<i128, String> {
    row.at(key)?
        .as_i128()
        .ok_or_else(|| format!("suite row `{key}` is not an integer"))
}

/// The suite's rows of `tier`, keyed as a new run tree names them:
/// `<slug>/<recipe>-<target>` from the frozen `k<a>n<n>-<recipe>-<target>`.
pub fn rows(programme: &Path, tier: &str) -> Result<Vec<Row>, String> {
    let doc = json::read(&programme.join("suite").join("v1").join("SUITE.json"))?;
    let mut out = Vec::new();
    for row in doc.at("rows")?.as_arr().ok_or("`rows` is not a list")? {
        if row.at("tier")?.as_str() != Some(tier) {
            continue;
        }
        let (a, n) = (int(row, "a")? as u8, int(row, "n")? as u32);
        let frozen = row.at("id")?.as_str().ok_or("`id` is not a string")?;
        let (_, rest) = frozen
            .split_once('-')
            .ok_or_else(|| format!("suite id `{frozen}` has no recipe part"))?;
        out.push(Row {
            id: format!("{}/{rest}", curve_slug(a, n)?),
            a,
            n,
            r: int(row, "r")? as f64,
            recipe_seed: Some(int(row, "recipe_seed")?),
        });
    }
    Ok(out)
}

/// The rows grouped by size, the sizes in increasing `r`.
pub fn by_size(rows: &[Row]) -> Vec<((u8, u32), Vec<Row>)> {
    let mut out: Vec<((u8, u32), Vec<Row>)> = Vec::new();
    for r in rows {
        match out.iter_mut().find(|(k, _)| *k == (r.a, r.n)) {
            Some((_, v)) => v.push(r.clone()),
            None => out.push(((r.a, r.n), vec![r.clone()])),
        }
    }
    out.sort_by(|x, y| x.1[0].r.partial_cmp(&y.1[0].r).expect("finite r"));
    out
}

/// Every `<sha256>  <path>` line of a `SHA256SUMS` in `dir` holds.
pub fn check_sums(dir: &Path) -> Result<usize, String> {
    let sums = dir.join("SHA256SUMS");
    let text = std::fs::read_to_string(&sums).map_err(|e| format!("{}: {e}", sums.display()))?;
    let mut n = 0;
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let (want, name) = line
            .split_once("  ")
            .ok_or_else(|| format!("{}: malformed line `{line}`", sums.display()))?;
        let path = dir.join(name);
        let bytes = std::fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let got = hex::encode(crypto_lib::hash::sha256::sha256(&bytes));
        if got != want {
            return Err(format!("{} does not match SHA256SUMS", path.display()));
        }
        n += 1;
    }
    Ok(n)
}
