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
    /// The parameter file a timed process prices, absolute.
    pub params: Option<PathBuf>,
    /// The rho walk's seed for the row's single target.
    pub rho_seed: Option<i128>,
    /// The suite's frozen id, `k<a>n<n>-<recipe>-<target>`, for a suite row.
    pub suite_id: Option<String>,
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
        let params = row
            .at("params")?
            .as_str()
            .ok_or("`params` is not a string")?;
        out.push(Row {
            id: format!("{}/{rest}", curve_slug(a, n)?),
            a,
            n,
            r: int(row, "r")? as f64,
            recipe_seed: Some(int(row, "recipe_seed")?),
            params: Some(programme.join("suite").join("v1").join(params)),
            rho_seed: Some(int(row, "rho_seed")?),
            suite_id: Some(frozen.to_string()),
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

// ── fresh holdouts, by suite v1's own construction ──────────────────

/// Suite v1's target and rho-seed bases (`make_suite.py`).
const TARGET_SEED_BASE: i128 = 23000;

/// One holdout row's parameter file as suite v1's `make_suite.py` writes
/// it: §20's recipe rule (`research/ic_exponent_20260926/make_params.py`,
/// `params`) at the size's suite recipe, with the recipe seed and the
/// target free, renamed `-T<i>` and aimed at public hash target
/// `23000 + i`.  Python's `round` is half-to-even and its `int / int` is
/// correctly rounded; every operand here is an integer below `2^53`, so
/// `f64` reproduces both exactly.
fn holdout_text(a: u8, n: u32, r: i128, columns: i128, m: i128, seed: i128, i: i128) -> String {
    let n = i128::from(n);
    let points = 2 * n * columns;
    let window = ((points as f64) / 32.0).round_ties_even().max(1.0) as i128;
    let unit = ((2 * r) as f64 / (n * points * points) as f64)
        .round_ties_even()
        .clamp(16.0, 65536.0) as i128;
    let cap = (64 * ((2.02 * r as f64) / (points * points) as f64).ceil() as i128).max(100_000);
    let doc = J::Obj(vec![
        ("schema_version".into(), J::Int(1)),
        (
            "name".into(),
            J::Str(format!("k{a}n{n}-c{columns}-m{m}-s{seed}-T{i:02}")),
        ),
        (
            "curve".into(),
            json::obj([("degree", J::Int(n)), ("curve_a", J::Int(a.into()))]),
        ),
        ("summands".into(), J::Int(3)),
        ("descent_summands".into(), J::Int(m)),
        ("collection_window".into(), J::Int(window)),
        ("collection_aim".into(), J::Bool(true)),
        ("solver".into(), J::Str("pair_table".into())),
        ("seed".into(), J::Int(seed)),
        ("max_trials".into(), J::Int(cap.min((1 << 32) - 1))),
        (
            "collection".into(),
            json::obj([
                ("unit_trials", J::Int(unit)),
                ("units", J::Int(1)),
                ("max_units", J::Int(100_000)),
            ]),
        ),
        (
            "factor_base".into(),
            json::obj([
                ("mode", J::Str("spec".into())),
                (
                    "spec",
                    json::obj([
                        ("kind", J::Str("subgroup_orbits".into())),
                        ("points", J::Int(points)),
                        ("seed", J::Int(seed)),
                    ]),
                ),
            ]),
        ),
        (
            "targets".into(),
            J::Arr(vec![json::obj([(
                "public_hash_seed",
                J::Int(TARGET_SEED_BASE + i),
            )])]),
        ),
    ]);
    json::dumps(&doc, 1) + "\n"
}

/// A round's fresh holdouts, as `(relative path, text)` in `SHA256SUMS`
/// order: at each size, the suite's recipe (its `M1` row's columns and
/// descent summands) under each recipe seed, two targets to a seed,
/// numbered on from `first_target`.  Files are `<slug>/M<seed − 200>-T<i>.json`,
/// as R02b's and R05's were.
pub fn holdouts(
    programme: &Path,
    sizes: &[(u8, u32)],
    seeds: &[i128],
    first_target: i128,
) -> Result<Vec<(String, String)>, String> {
    let doc = json::read(&programme.join("suite").join("v1").join("SUITE.json"))?;
    let rows = doc.at("rows")?.as_arr().ok_or("`rows` is not a list")?;
    let mut out = Vec::new();
    for &(a, n) in sizes {
        let id = format!("k{a}n{n}-M1-T01");
        let m1 = rows
            .iter()
            .find(|r| r.get("id").and_then(J::as_str) == Some(id.as_str()))
            .ok_or_else(|| format!("suite v1 has no row {id}"))?;
        let slug = curve_slug(a, n)?;
        let (r, columns, m) = (
            int(m1, "r")?,
            int(m1, "columns")?,
            int(m1, "descent_summands")?,
        );
        for (j, &seed) in seeds.iter().enumerate() {
            for t in 0..2 {
                let i = first_target + 2 * j as i128 + t;
                out.push((
                    format!("{slug}/M{}-T{i}.json", seed - 200),
                    holdout_text(a, n, r, columns, m, seed, i),
                ));
            }
        }
    }
    out.sort_by(|x, y| x.0.cmp(&y.0));
    Ok(out)
}

/// The `SHA256SUMS` text for a set of files, in the order given.
pub fn sums_text(files: &[(String, String)]) -> String {
    files
        .iter()
        .map(|(path, text)| {
            format!(
                "{}  {path}\n",
                hex::encode(crypto_lib::hash::sha256::sha256(text.as_bytes()))
            )
        })
        .collect()
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
