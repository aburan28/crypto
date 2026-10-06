//! A speed round's figures and decision, from its run tree only
//! (IC_TOOL_PROGRAM.md §5 and each round's `PROTOCOL.md`).
//!
//! Every ratio is the base over the candidate, so a ratio above 1 means
//! the candidate is faster.  The machinery is shared; each round names
//! its rows, its measures and its acceptance rule.  `r03` reproduces
//! R03's frozen `analysis.json` from its frozen runs, which is how the
//! native analysis was checked against the one it replaces.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use super::bench::{self, Arm, Bench};
use super::callgrind;
use super::json::{self, obj, opt_f64, J};
use super::pin;
use super::runs;
use super::stats;
use super::suite::{self, Row};

const LEVEL: f64 = 0.95;

pub struct Ctx {
    pub programme: PathBuf,
    pub round_dir: PathBuf,
    pub runs: PathBuf,
}

type Measure = fn(&J) -> Result<f64, String>;

fn report(d: &Path, arm: &str, row: &str, k: u32) -> Result<Option<J>, String> {
    runs::figure(&d.join(arm).join(row).join(format!("r{k}.price.json")))
}

/// Every figure pair of two arms, by the first arm's rounds.
fn pairs_of(
    d: &Path,
    rows: &[Row],
    arms: (&str, &str),
) -> Result<Vec<(Option<J>, Option<J>)>, String> {
    let mut out = Vec::new();
    for r in rows {
        for k in runs::rounds_of(d, arms.0, &r.id) {
            out.push((report(d, arms.0, &r.id, k)?, report(d, arms.1, &r.id, k)?));
        }
    }
    Ok(out)
}

/// Every (base, candidate) figure pair, by the base's rounds.
fn pairs(d: &Path, rows: &[Row]) -> Result<Vec<(Option<J>, Option<J>)>, String> {
    pairs_of(d, rows, ("base", "cand"))
}

/// The paired ratios' interval, base over candidate, with the pairs that
/// had no figure.
fn paired(d: &Path, rows: &[Row], measure: Measure, half_width: bool) -> Result<J, String> {
    paired_of(d, rows, ("base", "cand"), measure, half_width)
}

/// The same for any two arms, the first over the second.
pub(crate) fn paired_of(
    d: &Path,
    rows: &[Row],
    arms: (&str, &str),
    measure: Measure,
    half_width: bool,
) -> Result<J, String> {
    let (mut ratios, mut missing) = (Vec::new(), 0i128);
    for (a, b) in pairs_of(d, rows, arms)? {
        match (a, b) {
            (Some(a), Some(b)) => ratios.push(measure(&a)? / measure(&b)?),
            _ => missing += 1,
        }
    }
    let ci = stats::geo_ci(&ratios, LEVEL);
    let mut kv = ci.fields();
    kv.push(("missing_pairs".into(), J::Int(missing)));
    if half_width {
        if let (Some(hi), Some(g)) = (ci.hi, ci.geomean) {
            kv.push(("half_width".into(), J::Float(hi / g - 1.0)));
        }
    }
    Ok(J::Obj(kv))
}

/// Every figure `arm` has for `rows`, row by row and round by round.
fn figures(d: &Path, arm: &str, rows: &[Row]) -> Result<Vec<J>, String> {
    let mut out = Vec::new();
    for r in rows {
        for k in runs::rounds_of(d, arm, &r.id) {
            if let Some(rep) = report(d, arm, &r.id, k)? {
                out.push(rep);
            }
        }
    }
    Ok(out)
}

/// `S` cold in v0's unit: per row the median over its rounds, then the
/// median over the size's rows.
fn s_cold(d: &Path, arm: &str, rows: &[Row], unit_ns: f64) -> Result<Option<f64>, String> {
    let mut per_row = Vec::new();
    for r in rows {
        let mut cold = Vec::new();
        for k in runs::rounds_of(d, arm, &r.id) {
            if let Some(rep) = report(d, arm, &r.id, k)? {
                cold.push(runs::ic_cold_ns(&rep)?);
            }
        }
        if !cold.is_empty() {
            per_row.push(stats::median(&cold));
        }
    }
    Ok((!per_row.is_empty()).then(|| stats::median(&per_row) / unit_ns / rows[0].r.sqrt()))
}

/// R01's A/A bands by slug: the cold interval of baseline against itself.
fn r01_aa(c: &Ctx) -> Result<HashMap<String, J>, String> {
    let doc = json::read(&c.programme.join("rounds/R01-baseline-v0/analysis.json"))?;
    let mut out = HashMap::new();
    for row in doc.at("aa")?.as_arr().ok_or("`aa` is not a list")? {
        let size = row.at("size")?.as_str().ok_or("`size` is not a string")?;
        // `k<a>n<n>`
        let a: u8 = size[1..2].parse().map_err(|_| format!("size `{size}`"))?;
        let n: u32 = size[3..].parse().map_err(|_| format!("size `{size}`"))?;
        out.insert(suite::curve_slug(a, n)?, row.at("cold")?.clone());
    }
    Ok(out)
}

/// v0's unit, nanoseconds per batched addition, by slug.
fn unit_v0(c: &Ctx) -> Result<HashMap<String, f64>, String> {
    let doc = json::read(&c.programme.join("baselines.json"))?;
    let v0 = doc
        .at("baselines")?
        .as_arr()
        .ok_or("`baselines` is not a list")?
        .iter()
        .find(|b| b.get("baseline").and_then(J::as_str) == Some("v0"))
        .ok_or("baselines.json has no v0")?;
    let mut out = HashMap::new();
    for s in v0.at("sizes")?.as_arr().ok_or("`sizes` is not a list")? {
        let alias = s.at("curve_alias")?.as_str().ok_or("`curve_alias`")?;
        let unit = s.at("unit_ns")?.as_f64().ok_or("`unit_ns`")?;
        out.insert(alias.to_string(), unit);
    }
    Ok(out)
}

pub(crate) fn round3(x: f64) -> f64 {
    format!("{x:.3}").parse().expect("a formatted float parses")
}

/// The A/A band's keys, and whether the row regresses beyond it.
pub(crate) fn aa_fields(kv: &mut Vec<(String, J)>, band: Option<&J>) -> Result<(), String> {
    let Some(band) = band.filter(|b| b.truthy()) else {
        return Ok(());
    };
    let cold_hi = kv
        .iter()
        .find(|(k, _)| k == "cold")
        .and_then(|(_, v)| v.get("hi"))
        .and_then(J::as_f64);
    let (Some(lo), Some(hi)) = (band.get("lo").and_then(J::as_f64), cold_hi) else {
        return Ok(());
    };
    kv.push((
        "aa_cold".into(),
        obj([
            ("geomean", band.at("geomean")?.clone()),
            ("lo", band.at("lo")?.clone()),
            ("hi", band.at("hi")?.clone()),
        ]),
    ));
    kv.push(("regresses_beyond_aa".into(), J::Bool(hi < lo)));
    Ok(())
}

pub(crate) fn size_head(slug: &str, a: u8, n: u32, rs: &[Row]) -> Vec<(String, J)> {
    vec![
        ("slug".into(), J::Str(slug.into())),
        ("a".into(), J::Int(a.into())),
        ("n".into(), J::Int(n.into())),
        ("log2_r".into(), J::Float(round3(rs[0].r.log2()))),
        ("rows".into(), J::Int(rs.len() as i128)),
    ]
}

/// The reasons a target fails its acceptance bound on one set.
fn bound_reasons(rows: &J, slug: &str, name: &str, accept_lo: f64, reasons: &mut Vec<J>) {
    let row = rows
        .as_arr()
        .unwrap_or(&[])
        .iter()
        .find(|r| r.get("slug").and_then(J::as_str) == Some(slug));
    let lo = row
        .and_then(|r| r.get("cold"))
        .and_then(|c| c.get("lo"))
        .and_then(J::as_f64);
    if lo.is_none_or(|lo| lo <= accept_lo) {
        let shown = lo.map_or("None".to_string(), json::py_float);
        reasons.push(J::Str(format!(
            "{slug} on the {name}: interval's lower end {shown} is not above {}",
            json::py_float(accept_lo)
        )));
    }
}

fn aa_reasons(suite: &J, reasons: &mut Vec<J>) {
    for r in suite.as_arr().unwrap_or(&[]) {
        if r.get("regresses_beyond_aa").is_some_and(J::truthy) {
            let slug = r.get("slug").and_then(J::as_str).unwrap_or("");
            reasons.push(J::Str(format!("{slug} regresses beyond its A/A band")));
        }
    }
}

pub(crate) fn pin_summary(pin: &Option<J>) -> Result<J, String> {
    Ok(match pin {
        Some(p) => obj([
            ("held", p.at("held")?.clone()),
            ("names_agree", p.at("names_agree")?.clone()),
            (
                "rows",
                J::Int(p.at("rows")?.as_arr().ok_or("`rows`")?.len() as i128),
            ),
        ]),
        None => J::Null,
    })
}

fn pin_flag(pin: &Option<J>, key: &str) -> bool {
    pin.as_ref().and_then(|p| p.get(key)).is_some_and(J::truthy)
}

fn accounting(c: &Ctx, steps: &[&str]) -> Result<J, String> {
    let mut kv = Vec::new();
    for step in steps {
        let d = c.runs.join(step);
        if d.exists() {
            kv.push((step.to_string(), runs::accounting(&d)?));
        }
    }
    Ok(J::Obj(kv))
}

pub(crate) fn opt(v: Option<J>) -> J {
    v.unwrap_or(J::Null)
}

// ── R03: the curve's construction at composite degrees ──────────────

pub mod r03 {
    use super::*;

    const TARGET: (u8, u32) = (0, 57);
    const ACCEPT_LO: f64 = 1.3;
    const COMPOSITE: [(u8, u32); 2] = [(1, 45), (0, 57)];

    fn construction(rep: &J) -> Result<f64, String> {
        runs::setup_phase_ns(rep, "setup")
    }

    fn setup(rep: &J) -> Result<f64, String> {
        runs::median_rep(rep, "setup_ns")
    }

    fn comparison(
        d: &Path,
        rows: &[Row],
        aa: &HashMap<String, J>,
        units: &HashMap<String, f64>,
    ) -> Result<J, String> {
        let mut out = Vec::new();
        for ((a, n), rs) in suite::by_size(rows) {
            let slug = suite::curve_slug(a, n)?;
            let mut kv = size_head(&slug, a, n, &rs);
            kv.push(("composite".into(), J::Bool(COMPOSITE.contains(&(a, n)))));
            kv.push(("cold".into(), paired(d, &rs, runs::ic_cold_ns, false)?));
            kv.push(("setup".into(), paired(d, &rs, setup, false)?));
            kv.push((
                "construction_stage".into(),
                paired(d, &rs, construction, false)?,
            ));
            let ps = pairs(d, &rs)?;
            let base_c = ps
                .iter()
                .filter_map(|(a, _)| a.as_ref().map(construction))
                .collect::<Result<Vec<_>, _>>()?;
            let cand_c = ps
                .iter()
                .filter_map(|(_, b)| b.as_ref().map(construction))
                .collect::<Result<Vec<_>, _>>()?;
            if !base_c.is_empty() && !cand_c.is_empty() {
                kv.push((
                    "construction_ms_median".into(),
                    obj([
                        ("base", J::Float(stats::median(&base_c) / 1e6)),
                        ("cand", J::Float(stats::median(&cand_c) / 1e6)),
                    ]),
                ));
            }
            if let Some(&unit) = units.get(&slug) {
                kv.push(("s_cold_base".into(), opt_f64(s_cold(d, "base", &rs, unit)?)));
                kv.push(("s_cold_cand".into(), opt_f64(s_cold(d, "cand", &rs, unit)?)));
            }
            aa_fields(&mut kv, aa.get(&slug))?;
            out.push(J::Obj(kv));
        }
        Ok(J::Arr(out))
    }

    pub fn analyse(c: &Ctx) -> Result<J, String> {
        let (aa, units) = (r01_aa(c)?, unit_v0(c)?);
        let rows: Vec<Row> = suite::rows(&c.programme, "S")?
            .into_iter()
            .filter(|r| COMPOSITE.contains(&(r.a, r.n)) || r.recipe_seed == Some(201))
            .collect();
        let mut holdout_rows = Vec::new();
        for (a, n) in COMPOSITE {
            let slug = suite::curve_slug(a, n)?;
            let r0 = rows
                .iter()
                .find(|r| (r.a, r.n) == (a, n))
                .ok_or("no suite row at a composite size")?;
            for i in [101, 102] {
                holdout_rows.push(Row {
                    id: format!("{slug}/M5-T{i}"),
                    a,
                    n,
                    r: r0.r,
                    recipe_seed: None,
                    params: None,
                    rho_seed: None,
                    suite_id: None,
                });
            }
        }
        let pin = json::read_opt(&c.runs.join("pin").join("pin.json"))?;
        let suite_j = comparison(&c.runs.join("compare"), &rows, &aa, &units)?;
        let held = comparison(
            &c.runs.join("holdout"),
            &holdout_rows,
            &HashMap::new(),
            &units,
        )?;
        let mut reasons = Vec::new();
        if !pin_flag(&pin, "held") {
            reasons.push(J::Str("an output differs from v0's".into()));
        }
        let slug = suite::curve_slug(TARGET.0, TARGET.1)?;
        bound_reasons(&suite_j, &slug, "suite", ACCEPT_LO, &mut reasons);
        bound_reasons(&held, &slug, "holdouts", ACCEPT_LO, &mut reasons);
        aa_reasons(&suite_j, &mut reasons);
        let decision = obj([
            ("accepted", J::Bool(reasons.is_empty())),
            ("reasons", J::Arr(reasons)),
        ]);
        Ok(obj([
            (
                "what_this_is",
                J::Str("R03, the curve's construction at composite degrees: every figure the README, the ledger and the scoreboard quote".into()),
            ),
            ("host", opt(json::read_opt(&c.runs.join("host.json"))?)),
            ("aa_source", opt(json::read_opt(&c.runs.join("aa-source.json"))?)),
            ("accounting", accounting(c, &["aa", "compare", "holdout"])?),
            ("pin", pin_summary(&pin)?),
            ("suite", suite_j),
            ("holdouts", held),
            ("decision", decision),
        ]))
    }
}

// ── speed rounds with R05's shape ───────────────────────────────────

/// A speed round with R05's shape, which R02b's protocol shares: `M1`'s
/// rows at every size and every suite row at the target sizes, eight frozen
/// holdout rows at each target size, five rounds, and rounds 6–10 for a set
/// whose interval is wider than 3%.
pub struct Spec {
    pub targets: &'static [(u8, u32)],
    pub holdouts: &'static [(i128, u32)],
    pub what_this_is: &'static str,
    /// The sizes its callgrind profiles cover, `M1`'s first target each;
    /// none for a round without them.
    pub callgrind: &'static [(u8, u32)],
    /// What those profiles decide.
    pub callgrind_role: CallgrindRole,
    /// What the target sizes must show.
    pub accept: Accept,
    /// The note `manifest-resumed` writes into `host-resumed.json`.
    pub resumed_note: &'static str,
}

/// What a round's target sizes must show to be accepted, beyond the pin
/// and the A/A bands that every round needs.
#[derive(Clone, Copy)]
pub enum Accept {
    /// The paired cold ratio's interval lies above this at every target
    /// size, on the suite rows and on the holdouts separately.
    Above(f64),
    /// No bound: a round that re-bases the programme on a newer commit is
    /// accepted on its pin and its A/A bands alone, and its figures are
    /// what the new baseline costs.
    Pinned,
}

/// What a round's callgrind profiles decide.
#[derive(Clone, Copy, PartialEq)]
pub enum CallgrindRole {
    /// R02b's control: its rule (the arithmetic kernels identical, the
    /// total within tolerance, the same logarithm) must hold.
    Control,
    /// Plan §5's deterministic cross-check: the instruction ratios are
    /// reported, and only a different logarithm decides.
    CrossCheck,
}

pub mod speed {
    use super::*;

    pub const ACCEPT_LO: f64 = 1.10;
    /// Suite v1's rho seeds: `0x230000 +` the target's number.
    pub const RHO_SEED_BASE: i128 = 0x230000;
    pub const ROUNDS: u32 = 5;
    pub const EXTENDED_ROUNDS: u32 = 10;
    pub const HALF_WIDTH_LIMIT: f64 = 0.03;

    fn collect(rep: &J) -> Result<f64, String> {
        runs::setup_phase_ns(rep, "collect")
    }

    fn per_summand(rep: &J) -> Result<f64, String> {
        let summands = rep
            .at("counts")?
            .at("pass")?
            .at("collect")?
            .at("summands_scanned")?
            .as_f64()
            .ok_or("`summands_scanned` is not a number")?;
        Ok(collect(rep)? / summands)
    }

    /// `M1`'s rows at every size: 22.
    pub fn m1_rows(c: &Ctx) -> Result<Vec<Row>, String> {
        Ok(suite::rows(&c.programme, "S")?
            .into_iter()
            .filter(|r| r.recipe_seed == Some(201))
            .collect())
    }

    /// `M1`'s 22 rows, and every suite row at the target sizes.
    pub fn suite_rows(c: &Ctx, s: &Spec) -> Result<Vec<Row>, String> {
        Ok(suite::rows(&c.programme, "S")?
            .into_iter()
            .filter(|r| r.recipe_seed == Some(201) || s.targets.contains(&(r.a, r.n)))
            .collect())
    }

    /// Eight fresh rows at each target size, from the frozen `holdouts/`,
    /// checked against its `SHA256SUMS`.
    pub fn holdout_rows(c: &Ctx, s: &Spec) -> Result<Vec<Row>, String> {
        let dir = c.round_dir.join("holdouts");
        suite::check_sums(&dir)?;
        let all = suite::rows(&c.programme, "S")?;
        let mut out = Vec::new();
        for &(a, n) in s.targets {
            let slug = suite::curve_slug(a, n)?;
            let r = all
                .iter()
                .find(|x| (x.a, x.n) == (a, n))
                .ok_or("no suite row at a target size")?
                .r;
            for &(seed, i) in s.holdouts {
                let name = format!("M{}-T{i}", seed - 200);
                let path = dir.join(&slug).join(format!("{name}.json"));
                if !path.exists() {
                    return Err(format!("{} is missing", path.display()));
                }
                out.push(Row {
                    id: format!("{slug}/{name}"),
                    a,
                    n,
                    r,
                    recipe_seed: Some(seed),
                    params: Some(path),
                    rho_seed: Some(RHO_SEED_BASE + i as i128),
                    suite_id: None,
                });
            }
        }
        Ok(out)
    }

    fn comparison(
        d: &Path,
        rows: &[Row],
        aa: &HashMap<String, J>,
        units: &HashMap<String, f64>,
    ) -> Result<J, String> {
        let mut out = Vec::new();
        for ((a, n), rs) in suite::by_size(rows) {
            let slug = suite::curve_slug(a, n)?;
            let mut kv = size_head(&slug, a, n, &rs);
            let rounds = rs
                .iter()
                .map(|r| runs::rounds_of(d, "base", &r.id).len())
                .max()
                .unwrap_or(0);
            kv.push(("rounds".into(), J::Int(rounds as i128)));
            kv.push(("cold".into(), paired(d, &rs, runs::ic_cold_ns, true)?));
            kv.push(("online".into(), paired(d, &rs, runs::ic_online_ns, true)?));
            kv.push(("collect_stage".into(), paired(d, &rs, collect, true)?));
            kv.push((
                "rho_online".into(),
                paired(d, &rs, runs::rho_online_ns, true)?,
            ));
            if let Some(&unit) = units.get(&slug) {
                kv.push(("s_cold_base".into(), opt_f64(s_cold(d, "base", &rs, unit)?)));
                kv.push(("s_cold_cand".into(), opt_f64(s_cold(d, "cand", &rs, unit)?)));
                let base = figures(d, "base", &rs)?;
                let cand = figures(d, "cand", &rs)?;
                if !base.is_empty() && !cand.is_empty() {
                    let med = |reps: &[J]| -> Result<f64, String> {
                        let xs = reps
                            .iter()
                            .map(per_summand)
                            .collect::<Result<Vec<_>, _>>()?;
                        Ok(stats::median(&xs) / unit)
                    };
                    kv.push((
                        "collect_units_per_summand".into(),
                        obj([
                            ("base", J::Float(med(&base)?)),
                            ("cand", J::Float(med(&cand)?)),
                        ]),
                    ));
                }
            }
            aa_fields(&mut kv, aa.get(&slug))?;
            out.push(J::Obj(kv));
        }
        Ok(J::Arr(out))
    }

    /// The round's own A/A, for a host that is not R01's: the base against
    /// a byte-identical copy (`A` and `A2`) on `M1`'s rows, size by size, as
    /// R01 measured its own.  Also the bands by slug.
    fn own_aa(c: &Ctx) -> Result<(J, HashMap<String, J>), String> {
        let d = c.runs.join("aa");
        let arms = ("A", "A2");
        let (mut out, mut bands) = (Vec::new(), HashMap::new());
        for ((a, n), rs) in suite::by_size(&m1_rows(c)?) {
            let slug = suite::curve_slug(a, n)?;
            let mut kv = size_head(&slug, a, n, &rs);
            let rounds = rs
                .iter()
                .map(|r| runs::rounds_of(&d, arms.0, &r.id).len())
                .max()
                .unwrap_or(0);
            kv.push(("rounds".into(), J::Int(rounds as i128)));
            let cold = paired_of(&d, &rs, arms, runs::ic_cold_ns, false)?;
            kv.push(("cold".into(), cold.clone()));
            kv.push((
                "online".into(),
                paired_of(&d, &rs, arms, runs::ic_online_ns, false)?,
            ));
            kv.push((
                "rho_online".into(),
                paired_of(&d, &rs, arms, runs::rho_online_ns, false)?,
            ));
            if cold.get("lo").is_some() {
                bands.insert(slug, cold);
            }
            out.push(J::Obj(kv));
        }
        Ok((J::Arr(out), bands))
    }

    /// Whether the round's host is R01's, as its `aa-source.json` records.
    fn host_matches_r01(c: &Ctx) -> Result<bool, String> {
        Ok(json::read_opt(&c.runs.join("aa-source.json"))?
            .and_then(|d| d.get("host_matches_r01").cloned())
            .is_none_or(|v| v.truthy()))
    }

    /// The round's figures and decision, from its run tree only.
    pub fn analyse(c: &Ctx, s: &Spec) -> Result<J, String> {
        let units = unit_v0(c)?;
        let same_host = host_matches_r01(c)?;
        let (own, aa) = if same_host {
            (None, r01_aa(c)?)
        } else {
            let (record, bands) = own_aa(c)?;
            (Some(record), bands)
        };
        let pin = json::read_opt(&c.runs.join("pin").join("pin.json"))?;
        let suite_rows = suite_rows(c, s)?;
        let suite_j = comparison(&c.runs.join("compare"), &suite_rows, &aa, &units)?;
        let holdouts = comparison(
            &c.runs.join("holdout"),
            &holdout_rows(c, s)?,
            &HashMap::new(),
            &units,
        )?;
        let control = if s.callgrind.is_empty() {
            None
        } else {
            Some(callgrind::control(&c.runs.join("callgrind"))?)
        };
        let mut reasons = Vec::new();
        if let Some(cg) = &control {
            match s.callgrind_role {
                CallgrindRole::Control => {
                    if !cg.get("control_held").is_some_and(J::truthy) {
                        reasons.push(J::Str("the callgrind control did not hold".into()));
                    }
                }
                CallgrindRole::CrossCheck => {
                    let pairs = cg.get("pairs").and_then(J::as_obj).unwrap_or(&[]);
                    if pairs.len() != s.callgrind.len() {
                        reasons.push(J::Str(format!(
                            "the callgrind cross-check has {} of its {} profile pairs",
                            pairs.len(),
                            s.callgrind.len()
                        )));
                    }
                    for (row, p) in pairs {
                        if !p.get("same_log").is_some_and(J::truthy) {
                            reasons.push(J::Str(format!(
                                "{row}: the arms under callgrind did not recover the same logarithm"
                            )));
                        }
                    }
                }
            }
        }
        if !pin_flag(&pin, "held") {
            reasons.push(J::Str("an output differs from v0's".into()));
        }
        if !pin_flag(&pin, "names_agree") {
            reasons.push(J::Str("a name disagrees with the registry".into()));
        }
        if let Accept::Above(accept_lo) = s.accept {
            for &(a, n) in s.targets {
                let slug = suite::curve_slug(a, n)?;
                bound_reasons(&suite_j, &slug, "suite", accept_lo, &mut reasons);
                bound_reasons(&holdouts, &slug, "holdouts", accept_lo, &mut reasons);
            }
        }
        aa_reasons(&suite_j, &mut reasons);
        if !same_host {
            for ((a, n), _) in suite::by_size(&suite_rows) {
                let slug = suite::curve_slug(a, n)?;
                if !aa.contains_key(&slug) {
                    reasons.push(J::Str(format!(
                        "{slug}: the host is not R01's, and the round's own A/A has no band there"
                    )));
                }
            }
        }
        let extension = extension_record(c, s)?;
        let recorded = json::read_opt(&c.runs.join("extended.json"))?;
        if let Some(rec) = &recorded {
            let listed = |e: &J| -> bool {
                let (set, slug) = (e.get("set"), e.get("slug"));
                rec.get(set.and_then(J::as_str).unwrap_or(""))
                    .and_then(J::as_arr)
                    .is_some_and(|v| v.iter().any(|s| Some(s) == slug))
            };
            let tests = extension.as_arr().unwrap_or(&[]);
            let agree = tests
                .iter()
                .all(|e| listed(e) == (e.get("extended") == Some(&J::Bool(true))));
            if !agree {
                reasons.push(J::Str(
                    "the extension record differs from the rule's test".into(),
                ));
            }
        }
        let decision = obj([
            ("accepted", J::Bool(reasons.is_empty())),
            ("reasons", J::Arr(reasons)),
            (
                "note",
                J::Str("the tests run with cargo before the timed steps and are recorded beside this analysis".into()),
            ),
        ]);
        let mut by_tool = Vec::new();
        for step in ["aa", "compare", "holdout"] {
            let d = c.runs.join(step);
            if step != "aa" || d.exists() {
                by_tool.push((step.to_string(), runs::by_tool(&d)?));
            }
        }
        let mut doc = vec![
            ("what_this_is".to_string(), J::Str(s.what_this_is.into())),
            (
                "host".into(),
                opt(json::read_opt(&c.runs.join("host.json"))?),
            ),
            (
                "host_resumed".into(),
                opt(json::read_opt(&c.runs.join("host-resumed.json"))?),
            ),
        ];
        let later = later_resumes(c)?;
        if !later.is_empty() {
            doc.push(("hosts_resumed_later".into(), J::Arr(later)));
        }
        doc.push((
            "aa_source".into(),
            opt(json::read_opt(&c.runs.join("aa-source.json"))?),
        ));
        if let Some(own) = own {
            doc.push(("aa".into(), own));
        }
        doc.extend([
            (
                "accounting".to_string(),
                accounting(c, &["aa", "compare", "holdout"])?,
            ),
            ("isolation_by_tool".into(), J::Obj(by_tool)),
            ("pin".into(), pin_summary(&pin)?),
            ("extended".into(), opt(recorded)),
            ("extension_test".into(), extension),
            ("suite".into(), suite_j),
            ("holdouts".into(), holdouts),
        ]);
        if let Some(cg) = control {
            doc.push(("callgrind".into(), callgrind_for_role(cg, s.callgrind_role)));
        }
        doc.push(("decision".into(), decision));
        Ok(J::Obj(doc))
    }

    /// The callgrind block as the round's role reads it.  A control keeps
    /// [`callgrind::control`]'s verdict and rule.  A cross-check (plan
    /// §5) decides one thing, the same logarithm from both arms, so its
    /// block says that and nothing about a control it is not.
    fn callgrind_for_role(cg: J, role: CallgrindRole) -> J {
        match role {
            CallgrindRole::Control => cg,
            CallgrindRole::CrossCheck => {
                let pairs = cg.get("pairs").cloned().unwrap_or(J::Obj(Vec::new()));
                let held = pairs.as_obj().is_some_and(|ps| {
                    !ps.is_empty()
                        && ps
                            .iter()
                            .all(|(_, p)| p.get("same_log").is_some_and(J::truthy))
                });
                obj([
                    ("pairs", pairs),
                    ("role", J::Str("cross-check".into())),
                    ("cross_check_held", J::Bool(held)),
                    (
                        "rule",
                        J::Str(
                            "both arms recover the same logarithm at every profiled row; the instruction \
                             counts and the functions that differ are reported, not tested \
                             (PROTOCOL.md, plan §5)"
                                .into(),
                        ),
                    ),
                ])
            }
        }
    }

    // ── the declared runs, natively (each PROTOCOL.md's "Rows") ──────

    fn status(rep: &J) -> Option<&str> {
        rep.get("status").and_then(J::as_str)
    }

    /// The cold ratios the extension rule reads: each round's figure, or
    /// its first attempt when none is clean, as the declared runner read them.
    fn cold_ratios(d: &Path, rs: &[Row], rounds: u32) -> Result<Vec<f64>, String> {
        let mut out = Vec::new();
        for r in rs {
            for k in 1..=rounds {
                let path = |arm: &str| d.join(arm).join(&r.id).join(format!("r{k}.price.json"));
                let a = runs::load(&runs::figure_path(&path("base"))?);
                let b = runs::load(&runs::figure_path(&path("cand"))?);
                if status(&a) == Some("complete") && status(&b) == Some("complete") {
                    out.push(runs::ic_cold_ns(&a)? / runs::ic_cold_ns(&b)?);
                }
            }
        }
        Ok(out)
    }

    fn sets(c: &Ctx, s: &Spec) -> Result<[(&'static str, PathBuf, Vec<Row>); 2], String> {
        Ok([
            ("suite", c.runs.join("compare"), suite_rows(c, s)?),
            ("holdouts", c.runs.join("holdout"), holdout_rows(c, s)?),
        ])
    }

    /// The extension rule's test, set by set at each target size: the
    /// cold interval over the first five rounds and its half-width.  The
    /// test reads the interval's width only.
    fn extension_tests(
        c: &Ctx,
        s: &Spec,
    ) -> Result<Vec<(&'static str, String, stats::GeoCi)>, String> {
        let mut out = Vec::new();
        for (name, d, rows) in sets(c, s)? {
            for &(a, n) in s.targets {
                let rs: Vec<Row> = rows
                    .iter()
                    .filter(|r| (r.a, r.n) == (a, n))
                    .cloned()
                    .collect();
                let ci = stats::geo_ci(&cold_ratios(&d, &rs, ROUNDS)?, LEVEL);
                out.push((name, suite::curve_slug(a, n)?, ci));
            }
        }
        Ok(out)
    }

    fn half_width(ci: &stats::GeoCi) -> Option<f64> {
        Some(ci.hi? / ci.geomean? - 1.0)
    }

    /// What the extension rule read, for the analysis: rounds 1–5 only,
    /// whatever ran after.
    fn extension_record(c: &Ctx, s: &Spec) -> Result<J, String> {
        let mut out = Vec::new();
        for (set, slug, ci) in extension_tests(c, s)? {
            let mut kv = vec![
                ("set".to_string(), J::Str(set.into())),
                ("slug".to_string(), J::Str(slug)),
                ("rounds".to_string(), J::Int(ROUNDS as i128)),
            ];
            kv.extend(ci.fields());
            kv.push(("half_width".into(), opt_f64(half_width(&ci))));
            kv.push((
                "extended".into(),
                J::Bool(half_width(&ci).is_some_and(|h| h > HALF_WIDTH_LIMIT)),
            ));
            out.push(J::Obj(kv));
        }
        Ok(J::Arr(out))
    }

    /// Rounds 6–10 for a set whose half-width at a target size exceeds 3%
    /// after five rounds.  The test reads the interval's width only.
    pub fn extend(c: &Ctx, s: &Spec, b: &Bench, arms: &[Arm]) -> Result<J, String> {
        let record = c.runs.join("extended.json");
        let sets = sets(c, s)?;
        let done = if record.exists() {
            json::read(&record)?
        } else {
            let mut kv: Vec<(String, J)> = Vec::new();
            for (name, slug, ci) in extension_tests(c, s)? {
                if half_width(&ci).is_some_and(|h| h > HALF_WIDTH_LIMIT) {
                    let slug = J::Str(slug);
                    match kv.iter_mut().find(|(k, _)| k == name) {
                        Some((_, J::Arr(v))) => v.push(slug),
                        _ => kv.push((name.to_string(), J::Arr(vec![slug]))),
                    }
                }
            }
            let doc = J::Obj(kv);
            std::fs::write(&record, json::dumps(&doc, 1) + "\n")
                .map_err(|e| format!("{}: {e}", record.display()))?;
            doc
        };
        for (name, d, rows) in &sets {
            let slugs: Vec<&str> = done
                .get(name)
                .and_then(J::as_arr)
                .unwrap_or(&[])
                .iter()
                .filter_map(J::as_str)
                .collect();
            let mut rs = Vec::new();
            for r in rows {
                if slugs.contains(&suite::curve_slug(r.a, r.n)?.as_str()) {
                    rs.push(r.clone());
                }
            }
            if !rs.is_empty() {
                b.interleave(arms, &rs, EXTENDED_ROUNDS, d, &[])?;
            }
        }
        Ok(done)
    }

    /// What `compare` and `holdout` would still run, in their order: the
    /// rounds with no first attempt yet, and the attempts owed a retry.
    pub fn plan(c: &Ctx, s: &Spec, arms: &[Arm]) -> Result<J, String> {
        let mut out = Vec::new();
        for (name, d, rows) in sets(c, s)? {
            let mut missing = Vec::new();
            for k in 1..=ROUNDS {
                let order: Vec<&Arm> = if k % 2 == 1 {
                    arms.iter().collect()
                } else {
                    arms.iter().rev().collect()
                };
                for r in &rows {
                    for arm in &order {
                        let first = d
                            .join(&arm.name)
                            .join(&r.id)
                            .join(format!("r{k}.price.json"));
                        let figure = runs::figure_path(&first)?;
                        let done = figure.exists()
                            && runs::clean(&figure)?
                            && status(&runs::load(&figure)) == Some("complete");
                        let tried = (0..=runs::RETRIES)
                            .filter(|&a| runs::attempt(&first, a).exists())
                            .count();
                        if !done && tried <= runs::RETRIES {
                            missing.push(J::Str(format!(
                                "r{k} {} {}{}",
                                r.id,
                                arm.name,
                                if tried > 0 { " (retry)" } else { "" }
                            )));
                        }
                    }
                }
            }
            out.push((
                name.to_string(),
                json::obj([
                    ("to_run", J::Int(missing.len() as i128)),
                    ("processes", J::Arr(missing)),
                ]),
            ));
        }
        Ok(J::Obj(out))
    }

    /// Each arm's binary, as a manifest names it: its file name, SHA-256 and,
    /// where given, the commit it was built from.
    pub(crate) fn binaries(arms: &[Arm], commits: &[Option<String>]) -> Result<J, String> {
        let mut kv = Vec::new();
        for (i, a) in arms.iter().enumerate() {
            let name = if a.name == "cand" {
                "candidate".to_string()
            } else {
                a.name.clone()
            };
            let mut fields = vec![
                (
                    "path_basename".to_string(),
                    J::Str(
                        a.binary
                            .file_name()
                            .map(|f| f.to_string_lossy().into_owned())
                            .unwrap_or_default(),
                    ),
                ),
                ("sha256".into(), J::Str(bench::sha256_file(&a.binary)?)),
            ];
            if let Some(Some(commit)) = commits.get(i) {
                fields.push(("built_from".into(), J::Str(commit.clone())));
            }
            kv.push((name, J::Obj(fields)));
        }
        Ok(J::Obj(kv))
    }

    /// The host a round starts on (`host.json`), and whether it is R01's
    /// (`aa-source.json`): the A/A bands are R01's if so, the round's own
    /// if not.  Both are written once.
    pub fn manifest(
        c: &Ctx,
        b: &Bench,
        arms: &[Arm],
        root: &Path,
        commits: &[Option<String>],
        note: &str,
    ) -> Result<J, String> {
        let doc = bench::host_manifest(
            &c.runs.join("host.json"),
            root,
            binaries(arms, commits)?,
            &b.isolate,
            note,
        )?;
        let out = c.runs.join("aa-source.json");
        if out.exists() {
            return json::read(&out);
        }
        let r01 = json::read(&pin::r01_runs(&c.programme)?.join("host.json"))?;
        let same = [
            "cpu_model",
            "cpu_flags_relevant",
            "logical_cores",
            "memory",
            "transparent_hugepage",
            "os",
        ]
        .iter()
        .all(|k| match (doc.get(k), r01.get(k)) {
            (Some(x), Some(y)) => json::py_eq(x, y),
            _ => false,
        });
        let record = obj([
            ("host_matches_r01", J::Bool(same)),
            (
                "aa",
                J::Str(if same { "R01's" } else { "the round's own" }.into()),
            ),
        ]);
        std::fs::write(&out, json::dumps(&record, 1) + "\n")
            .map_err(|e| format!("{}: {e}", out.display()))?;
        Ok(record)
    }

    /// The host a round resumed on after its container changed:
    /// `host-resumed.json` the first time, and `host-resumed-<k>.json`,
    /// `k = 2, 3, …`, at each later resume, so that every container a
    /// round ran on is on record.  None is ever overwritten.
    pub fn manifest_resumed(
        c: &Ctx,
        s: &Spec,
        b: &Bench,
        arms: &[Arm],
        root: &Path,
    ) -> Result<J, String> {
        let first = c.runs.join("host-resumed.json");
        let out = if first.exists() {
            (2..)
                .map(|k| c.runs.join(format!("host-resumed-{k}.json")))
                .find(|p| !p.exists())
                .expect("an unbounded range has a free name")
        } else {
            first
        };
        bench::host_manifest(&out, root, binaries(arms, &[])?, &b.isolate, s.resumed_note)
    }

    /// The hosts of every resume after the first, in order: empty for a
    /// round that resumed at most once.
    fn later_resumes(c: &Ctx) -> Result<Vec<J>, String> {
        let mut out = Vec::new();
        for k in 2.. {
            match json::read_opt(&c.runs.join(format!("host-resumed-{k}.json")))? {
                Some(doc) => out.push(doc),
                None => break,
            }
        }
        Ok(out)
    }

    /// The round's own A/A: the base against a byte-identical copy beside
    /// it, on `M1`'s 22 rows, five rounds, as R01's.
    pub fn aa(c: &Ctx, b: &Bench, arms: &[Arm]) -> Result<J, String> {
        let base = &arms
            .iter()
            .find(|a| a.name == "base")
            .ok_or("no base arm")?
            .binary;
        let name = base
            .file_name()
            .ok_or("the base arm has no file name")?
            .to_string_lossy();
        let copy = base.with_file_name(format!("{name}-aa-copy"));
        if !copy.exists() {
            std::fs::copy(base, &copy).map_err(|e| format!("{}: {e}", copy.display()))?;
        }
        if bench::sha256_file(&copy)? != bench::sha256_file(base)? {
            return Err(format!(
                "{} is not byte-identical to the base",
                copy.display()
            ));
        }
        let pair = [
            Arm {
                name: "A".into(),
                binary: base.clone(),
            },
            Arm {
                name: "A2".into(),
                binary: copy,
            },
        ];
        b.interleave(&pair, &m1_rows(c)?, ROUNDS, &c.runs.join("aa"), &[])?;
        Ok(J::Null)
    }

    /// The callgrind control's profiles: `ic workflow` under callgrind on
    /// `M1`'s first target at each of the round's sizes, both arms, untimed
    /// under `taskset`; then each profile's phases, natively.
    pub fn callgrind_step(c: &Ctx, s: &Spec, arms: &[Arm]) -> Result<J, String> {
        let d = c.runs.join("callgrind");
        std::fs::create_dir_all(&d).map_err(|e| format!("{}: {e}", d.display()))?;
        let rows = suite_rows(c, s)?;
        for arm in arms {
            for &(a, n) in s.callgrind {
                let row = rows
                    .iter()
                    .find(|r| (r.a, r.n) == (a, n) && r.id.ends_with("/M1-T01"))
                    .ok_or("no M1-T01 row at a callgrind size")?;
                let tag = format!("{}-{}", arm.name, row.id.replace('/', "-"));
                let stem = format!("{tag}.callgrind.out");
                let phases = d.join(format!("{tag}.phases.json"));
                if phases.exists() {
                    continue;
                }
                let params = row.params.as_ref().ok_or("a suite row has no parameters")?;
                let work = std::env::temp_dir()
                    .join(format!("icprog-workflow-{}-{tag}", std::process::id()));
                std::fs::create_dir_all(&work).map_err(|e| format!("{}: {e}", work.display()))?;
                let stdout = std::fs::File::create(d.join(format!("{tag}.workflow.json")))
                    .map_err(|e| format!("{tag}.workflow.json: {e}"))?;
                let stderr = std::fs::File::create(d.join(format!("{tag}.valgrind.log")))
                    .map_err(|e| format!("{tag}.valgrind.log: {e}"))?;
                std::process::Command::new("taskset")
                    .args([
                        "-c",
                        bench::CPUS,
                        "valgrind",
                        "--tool=callgrind",
                        "--cache-sim=yes",
                    ])
                    .args(["--I1=32768,8,64", "--D1=32768,8,64", "--LL=2097152,16,64"])
                    .arg(format!("--callgrind-out-file={}", d.join(&stem).display()))
                    .arg(&arm.binary)
                    .args(["--json", "workflow", "--params"])
                    .arg(params)
                    .arg("--dir")
                    .arg(&work)
                    .env("RAYON_NUM_THREADS", "1")
                    .stdout(stdout)
                    .stderr(stderr)
                    .status()
                    .map_err(|e| format!("valgrind: {e}"))?;
                let _ = std::fs::remove_dir_all(&work);
                let doc = callgrind::phases(&d, &stem, 60)?;
                std::fs::write(&phases, json::dumps(&doc, 1) + "\n")
                    .map_err(|e| format!("{}: {e}", phases.display()))?;
                println!("callgrind {tag}: done");
            }
        }
        Ok(J::Null)
    }

    /// The steps every round of this shape has: `plan`, `compare`,
    /// `holdout` and `extend`.
    pub fn run_common(
        c: &Ctx,
        s: &Spec,
        step: &str,
        b: &Bench,
        arms: &[Arm],
    ) -> Option<Result<J, String>> {
        Some(match step {
            "plan" => plan(c, s, arms),
            "compare" => suite_rows(c, s)
                .and_then(|rows| b.interleave(arms, &rows, ROUNDS, &c.runs.join("compare"), &[]))
                .map(|()| J::Null),
            "holdout" => holdout_rows(c, s)
                .and_then(|rows| b.interleave(arms, &rows, ROUNDS, &c.runs.join("holdout"), &[]))
                .map(|()| J::Null),
            "extend" => extend(c, s, b, arms),
            _ => return None,
        })
    }
}

// ── R05: the sharper presence filter ────────────────────────────────

pub mod r05 {
    use super::*;

    pub const SPEC: Spec = Spec {
        targets: &[(0, 53), (1, 59), (0, 61)],
        holdouts: &[
            (210, 111),
            (210, 112),
            (211, 113),
            (211, 114),
            (212, 115),
            (212, 116),
            (213, 117),
            (213, 118),
        ],
        what_this_is: "R05, the sharper presence filter: every figure the README, the ledger and the scoreboard quote",
        callgrind: &[],
        callgrind_role: CallgrindRole::Control,
        accept: Accept::Above(speed::ACCEPT_LO),
        resumed_note: "the host R05 resumed on after its container was rebuilt mid-holdouts",
    };

    pub fn analyse(c: &Ctx) -> Result<J, String> {
        speed::analyse(c, &SPEC)
    }

    /// One declared step on the native runner.  `manifest` ran under the
    /// declared runner at the round's start, and its record is in the run
    /// tree.  `pin` reads the run tree's record, or recomputes it from the
    /// candidate outputs there when the record is absent.
    pub fn run(c: &Ctx, step: &str, b: &Bench, arms: &[Arm], root: &Path) -> Result<J, String> {
        if let Some(done) = speed::run_common(c, &SPEC, step, b, arms) {
            return done;
        }
        match step {
            "manifest-resumed" => speed::manifest_resumed(c, &SPEC, b, arms, root),
            "pin" => pin::pin(&c.programme, &c.runs, &arms[1].binary),
            "manifest" => Err(
                "`manifest` ran under the declared runner at the round's start; its record is in the run tree"
                    .into(),
            ),
            other => Err(format!(
                "unknown step `{other}`; try plan, manifest-resumed, pin, compare, holdout or extend"
            )),
        }
    }
}

// ── R02b: the wide-tail kernel, re-tested ───────────────────────────

pub mod r02b {
    use super::*;

    pub const SPEC: Spec = Spec {
        targets: &[(1, 59), (0, 61)],
        holdouts: &[
            (206, 103),
            (206, 104),
            (207, 105),
            (207, 106),
            (208, 107),
            (208, 108),
            (209, 109),
            (209, 110),
        ],
        what_this_is: "R02b, the wide-tail kernel re-tested on fresh holdouts: every figure the README, the ledger and the scoreboard quote",
        callgrind: &[(0, 61), (1, 59), (0, 41)],
        callgrind_role: CallgrindRole::Control,
        accept: Accept::Above(speed::ACCEPT_LO),
        resumed_note: "the host R02b resumed on after its container changed",
    };

    pub fn analyse(c: &Ctx) -> Result<J, String> {
        speed::analyse(c, &SPEC)
    }

    /// One of R02b's declared steps, natively (amendment 2): `manifest`,
    /// `pin`, `aa`, `compare`, `holdout`, `extend`, `callgrind`, and `plan`
    /// and `manifest-resumed` as for R05.
    pub fn run(
        c: &Ctx,
        step: &str,
        b: &Bench,
        arms: &[Arm],
        root: &Path,
        commits: &[Option<String>],
    ) -> Result<J, String> {
        if let Some(done) = speed::run_common(c, &SPEC, step, b, arms) {
            return done;
        }
        match step {
            "manifest" => speed::manifest(
                c,
                b,
                arms,
                root,
                commits,
                "the host R02b ran on, at the round's start",
            ),
            "manifest-resumed" => speed::manifest_resumed(c, &SPEC, b, arms, root),
            "pin" => pin::pin(&c.programme, &c.runs, &arms[1].binary),
            "aa" => speed::aa(c, b, arms),
            "callgrind" => speed::callgrind_step(c, &SPEC, arms),
            other => Err(format!(
                "unknown step `{other}`; try plan, manifest, pin, aa, compare, holdout, extend, callgrind or manifest-resumed"
            )),
        }
    }
}

// ── R07: the programme re-based on main's head ──────────────────────

pub mod r07 {
    use super::*;

    pub const SPEC: Spec = Spec {
        targets: &[(0, 53), (1, 59), (0, 61)],
        holdouts: &[
            (214, 119),
            (214, 120),
            (215, 121),
            (215, 122),
            (216, 123),
            (216, 124),
            (217, 125),
            (217, 126),
        ],
        what_this_is: "R07, the programme re-based on main's head: every figure the README, the ledger and the scoreboard quote",
        callgrind: &[(0, 61), (0, 53)],
        callgrind_role: CallgrindRole::CrossCheck,
        accept: Accept::Pinned,
        resumed_note: "the host R07 resumed on after its container changed",
    };

    pub fn analyse(c: &Ctx) -> Result<J, String> {
        speed::analyse(c, &SPEC)
    }

    /// One of R07's declared steps: those R02b has, with the callgrind
    /// profiles as plan §5's cross-check rather than a control.
    pub fn run(
        c: &Ctx,
        step: &str,
        b: &Bench,
        arms: &[Arm],
        root: &Path,
        commits: &[Option<String>],
    ) -> Result<J, String> {
        if let Some(done) = speed::run_common(c, &SPEC, step, b, arms) {
            return done;
        }
        match step {
            "manifest" => speed::manifest(
                c,
                b,
                arms,
                root,
                commits,
                "the host R07 ran on, at the round's start",
            ),
            "manifest-resumed" => speed::manifest_resumed(c, &SPEC, b, arms, root),
            "pin" => pin::pin(&c.programme, &c.runs, &arms[1].binary),
            "aa" => speed::aa(c, b, arms),
            "callgrind" => speed::callgrind_step(c, &SPEC, arms),
            other => Err(format!(
                "unknown step `{other}`; try plan, manifest, pin, aa, compare, holdout, extend, callgrind or manifest-resumed"
            )),
        }
    }
}

// ── R06: the scan's canonical key by funnel shifts ───────────────────

pub mod r06 {
    use super::*;

    /// R06's declaration: the two sizes the exploration found the key
    /// worth most at, eight fresh holdouts at each (recipe seeds 218 to
    /// 221, `T127` to `T134`), and an interval above 1.03 at both, on the
    /// suite rows and the holdouts separately.  Valgrind cannot run the
    /// candidate's VBMI2 kernel, so there are no callgrind profiles.
    pub const SPEC: Spec = Spec {
        targets: &[(1, 59), (0, 61)],
        holdouts: &[
            (218, 127),
            (218, 128),
            (219, 129),
            (219, 130),
            (220, 131),
            (220, 132),
            (221, 133),
            (221, 134),
        ],
        what_this_is: "R06, the scan's canonical key by funnel shifts: every figure the README, the ledger and the scoreboard quote",
        callgrind: &[],
        callgrind_role: CallgrindRole::CrossCheck,
        accept: Accept::Above(1.03),
        resumed_note: "the host R06 resumed on after its container changed",
    };

    pub fn analyse(c: &Ctx) -> Result<J, String> {
        speed::analyse(c, &SPEC)
    }

    /// One of R06's declared steps: R07's, without the callgrind profiles.
    pub fn run(
        c: &Ctx,
        step: &str,
        b: &Bench,
        arms: &[Arm],
        root: &Path,
        commits: &[Option<String>],
    ) -> Result<J, String> {
        if let Some(done) = speed::run_common(c, &SPEC, step, b, arms) {
            return done;
        }
        match step {
            "manifest" => speed::manifest(
                c,
                b,
                arms,
                root,
                commits,
                "the host R06 ran on, at the round's start",
            ),
            "manifest-resumed" => speed::manifest_resumed(c, &SPEC, b, arms, root),
            "pin" => pin::pin(&c.programme, &c.runs, &arms[1].binary),
            "aa" => speed::aa(c, b, arms),
            other => Err(format!(
                "unknown step `{other}`; try plan, manifest, pin, aa, compare, holdout, extend or manifest-resumed"
            )),
        }
    }
}
