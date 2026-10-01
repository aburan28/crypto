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
use super::json::{self, obj, opt_f64, J};
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

/// Every (base, candidate) figure pair, by the base's rounds.
fn pairs(d: &Path, rows: &[Row]) -> Result<Vec<(Option<J>, Option<J>)>, String> {
    let mut out = Vec::new();
    for r in rows {
        for k in runs::rounds_of(d, "base", &r.id) {
            out.push((report(d, "base", &r.id, k)?, report(d, "cand", &r.id, k)?));
        }
    }
    Ok(out)
}

/// The paired ratios' interval, with the pairs that had no figure.
fn paired(d: &Path, rows: &[Row], measure: Measure, half_width: bool) -> Result<J, String> {
    let (mut ratios, mut missing) = (Vec::new(), 0i128);
    for (a, b) in pairs(d, rows)? {
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

fn round3(x: f64) -> f64 {
    format!("{x:.3}").parse().expect("a formatted float parses")
}

/// The A/A band's keys, and whether the row regresses beyond it.
fn aa_fields(kv: &mut Vec<(String, J)>, band: Option<&J>) -> Result<(), String> {
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

fn size_head(slug: &str, a: u8, n: u32, rs: &[Row]) -> Vec<(String, J)> {
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

fn pin_summary(pin: &Option<J>) -> Result<J, String> {
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

fn opt(v: Option<J>) -> J {
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

// ── R05: the sharper presence filter ────────────────────────────────

pub mod r05 {
    use super::*;

    pub const TARGETS: [(u8, u32); 3] = [(0, 53), (1, 59), (0, 61)];
    pub const HOLDOUTS: [(i128, u32); 8] = [
        (210, 111),
        (210, 112),
        (211, 113),
        (211, 114),
        (212, 115),
        (212, 116),
        (213, 117),
        (213, 118),
    ];
    const ACCEPT_LO: f64 = 1.10;
    /// Suite v1's rho seeds: `0x230000 +` the target's number.
    const RHO_SEED_BASE: i128 = 0x230000;
    const ROUNDS: u32 = 5;
    const EXTENDED_ROUNDS: u32 = 10;
    const HALF_WIDTH_LIMIT: f64 = 0.03;

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

    /// `M1`'s 22 rows, and every suite row at the three target sizes: 40.
    pub fn suite_rows(c: &Ctx) -> Result<Vec<Row>, String> {
        Ok(suite::rows(&c.programme, "S")?
            .into_iter()
            .filter(|r| r.recipe_seed == Some(201) || TARGETS.contains(&(r.a, r.n)))
            .collect())
    }

    /// Eight fresh rows at each target size, from the frozen `holdouts/`.
    pub fn holdout_rows(c: &Ctx) -> Result<Vec<Row>, String> {
        let dir = c.round_dir.join("holdouts");
        suite::check_sums(&dir)?;
        let all = suite::rows(&c.programme, "S")?;
        let mut out = Vec::new();
        for (a, n) in TARGETS {
            let slug = suite::curve_slug(a, n)?;
            let r = all
                .iter()
                .find(|x| (x.a, x.n) == (a, n))
                .ok_or("no suite row at a target size")?
                .r;
            for (seed, i) in HOLDOUTS {
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

    pub fn analyse(c: &Ctx) -> Result<J, String> {
        let (aa, units) = (r01_aa(c)?, unit_v0(c)?);
        let pin = json::read_opt(&c.runs.join("pin").join("pin.json"))?;
        let suite_j = comparison(&c.runs.join("compare"), &suite_rows(c)?, &aa, &units)?;
        let holdouts = comparison(
            &c.runs.join("holdout"),
            &holdout_rows(c)?,
            &HashMap::new(),
            &units,
        )?;
        let mut reasons = Vec::new();
        if !pin_flag(&pin, "held") {
            reasons.push(J::Str("an output differs from v0's".into()));
        }
        if !pin_flag(&pin, "names_agree") {
            reasons.push(J::Str("a name disagrees with the registry".into()));
        }
        for (a, n) in TARGETS {
            let slug = suite::curve_slug(a, n)?;
            bound_reasons(&suite_j, &slug, "suite", ACCEPT_LO, &mut reasons);
            bound_reasons(&holdouts, &slug, "holdouts", ACCEPT_LO, &mut reasons);
        }
        aa_reasons(&suite_j, &mut reasons);
        let decision = obj([
            ("accepted", J::Bool(reasons.is_empty())),
            ("reasons", J::Arr(reasons)),
            (
                "note",
                J::Str("the tests run with cargo before the timed steps and are recorded beside this analysis".into()),
            ),
        ]);
        Ok(obj([
            (
                "what_this_is",
                J::Str("R05, the sharper presence filter: every figure the README, the ledger and the scoreboard quote".into()),
            ),
            ("host", opt(json::read_opt(&c.runs.join("host.json"))?)),
            (
                "host_resumed",
                opt(json::read_opt(&c.runs.join("host-resumed.json"))?),
            ),
            ("aa_source", opt(json::read_opt(&c.runs.join("aa-source.json"))?)),
            ("accounting", accounting(c, &["compare", "holdout"])?),
            (
                "isolation_by_tool",
                J::Obj(vec![
                    ("compare".into(), runs::by_tool(&c.runs.join("compare"))?),
                    ("holdout".into(), runs::by_tool(&c.runs.join("holdout"))?),
                ]),
            ),
            ("pin", pin_summary(&pin)?),
            ("extended", opt(json::read_opt(&c.runs.join("extended.json"))?)),
            ("suite", suite_j),
            ("holdouts", holdouts),
            ("decision", decision),
        ]))
    }

    // ── the declared runs, natively (PROTOCOL.md "What runs") ────────

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

    fn sets(c: &Ctx) -> Result<[(&'static str, PathBuf, Vec<Row>); 2], String> {
        Ok([
            ("suite", c.runs.join("compare"), suite_rows(c)?),
            ("holdouts", c.runs.join("holdout"), holdout_rows(c)?),
        ])
    }

    /// Rounds 6–10 for a set whose half-width at a target size exceeds 3%
    /// after five rounds.  The test reads the interval's width only.
    pub fn extend(c: &Ctx, b: &Bench, arms: &[Arm]) -> Result<J, String> {
        let record = c.runs.join("extended.json");
        let sets = sets(c)?;
        let done = if record.exists() {
            json::read(&record)?
        } else {
            let mut kv: Vec<(String, J)> = Vec::new();
            for (name, d, rows) in &sets {
                for (a, n) in TARGETS {
                    let rs: Vec<Row> = rows
                        .iter()
                        .filter(|r| (r.a, r.n) == (a, n))
                        .cloned()
                        .collect();
                    let ci = stats::geo_ci(&cold_ratios(d, &rs, ROUNDS)?, LEVEL);
                    if let (Some(hi), Some(g)) = (ci.hi, ci.geomean) {
                        if hi / g - 1.0 > HALF_WIDTH_LIMIT {
                            let slug = J::Str(suite::curve_slug(a, n)?);
                            match kv.iter_mut().find(|(k, _)| k == name) {
                                Some((_, J::Arr(v))) => v.push(slug),
                                _ => kv.push((name.to_string(), J::Arr(vec![slug]))),
                            }
                        }
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
    pub fn plan(c: &Ctx, arms: &[Arm]) -> Result<J, String> {
        let mut out = Vec::new();
        for (name, d, rows) in sets(c)? {
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

    /// One declared step on the native runner.  `manifest` and `pin`
    /// ran under the declared runner at the round's start; their records
    /// are in the run tree and are not re-run.  `manifest-resumed` records
    /// the host the round resumed on, once.
    pub fn run(c: &Ctx, step: &str, b: &Bench, arms: &[Arm], root: &Path) -> Result<J, String> {
        match step {
            "plan" => plan(c, arms),
            "manifest-resumed" => {
                let binaries = J::Obj(
                    arms.iter()
                        .map(|a| {
                            Ok((
                                if a.name == "cand" {
                                    "candidate".to_string()
                                } else {
                                    a.name.clone()
                                },
                                json::obj([
                                    (
                                        "path_basename",
                                        J::Str(
                                            a.binary
                                                .file_name()
                                                .map(|f| f.to_string_lossy().into_owned())
                                                .unwrap_or_default(),
                                        ),
                                    ),
                                    ("sha256", J::Str(bench::sha256_file(&a.binary)?)),
                                ]),
                            ))
                        })
                        .collect::<Result<Vec<_>, String>>()?,
                );
                bench::host_manifest(
                    &c.runs.join("host-resumed.json"),
                    root,
                    binaries,
                    &b.isolate,
                    "the host R05 resumed on after its container was rebuilt mid-holdouts",
                )
            }
            "compare" => {
                b.interleave(arms, &suite_rows(c)?, ROUNDS, &c.runs.join("compare"), &[])?;
                Ok(J::Null)
            }
            "holdout" => {
                b.interleave(arms, &holdout_rows(c)?, ROUNDS, &c.runs.join("holdout"), &[])?;
                Ok(J::Null)
            }
            "extend" => extend(c, b, arms),
            "manifest" | "pin" => Err(format!(
                "`{step}` ran under the declared runner at the round's start; its record is in the run tree"
            )),
            other => Err(format!(
                "unknown step `{other}`; try plan, manifest-resumed, compare, holdout or extend"
            )),
        }
    }
}
