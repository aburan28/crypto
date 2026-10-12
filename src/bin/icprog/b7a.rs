//! B7a's own measurements, natively (`rounds/B7a-f1-sampled/run.py` and
//! `analyse.py`, ported, with B7a's amendment 3: the chain on main).
//!
//! - `f1` (measurement 5): F1 on `M1`'s 22 rows, three rounds, isolated.
//!   Each row's document is its v2 translation by B7a's arm (`bround
//!   translate`'s) with `method.fidelity` set to `F1`.
//! - `partial` (measurement 7): F1 forced to partial tables,
//!   `--f1-partial` at 1/2 and 1/4, on `M1`'s rows at the six largest
//!   sizes, one round.
//! - `analyse`: measurements 5–7, the falsification target and the
//!   decision, with measurements 2–4 from `bround`.  F0's figures are the
//!   chain's B7a arm, the processes measurement 4 runs under the chain.
//!
//! A contended or failed process is kept and run again, at most twice; the
//! first clean, extrapolated attempt is the process's figure.

use std::path::{Path, PathBuf};

use super::bench::Bench;
use super::bround;
use super::json::{self, obj, J};
use super::rounds::{self, speed};
use super::runs;
use super::stats;
use super::suite::{self, Row};

const F1_ROUNDS: u32 = 3;
const PARTIAL_FRACTIONS: [f64; 2] = [0.5, 0.25];
const LARGEST: usize = 6;
const TIGHT_SIZES: [u32; 3] = [53, 59, 61];
const TIGHT_BAND: (f64, f64) = (0.85, 1.18);
const INSIDE_NEEDED: usize = 9;

fn fraction_tag(g: f64) -> String {
    // Python's `f"g{g}"`: 0.5 and 0.25.
    format!("g{}", json::py_float(g))
}

/// The row's v2 translation at `fidelity: F1`, written once.
fn f1_document(r: &bround::Run, cand: &Path, row: &Row) -> Result<PathBuf, String> {
    let path = r
        .ctx
        .runs
        .join("f1")
        .join("docs")
        .join(format!("{}.json", row.id));
    if !path.exists() {
        let J::Obj(mut kv) = bround::translate_row(r, cand, row)? else {
            return Err(format!("{}: the translation is not one document", row.id));
        };
        for (k, v) in kv.iter_mut() {
            if k == "method" {
                if let J::Obj(m) = v {
                    match m.iter_mut().find(|(mk, _)| mk == "fidelity") {
                        Some((_, f)) => *f = J::Str("F1".into()),
                        None => m.push(("fidelity".into(), J::Str("F1".into()))),
                    }
                }
            }
        }
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
        }
        std::fs::write(&path, json::dumps(&J::Obj(kv), 1) + "\n")
            .map_err(|e| format!("{}: {e}", path.display()))?;
    }
    Ok(path)
}

/// The first clean, extrapolated attempt: the process's figure.
fn figure(out: &Path) -> Result<PathBuf, String> {
    for k in 0..=runs::RETRIES {
        let attempt = runs::attempt(out, k);
        if attempt.exists()
            && runs::clean(&attempt)?
            && runs::load(&attempt).get("status").and_then(J::as_str) == Some("extrapolated")
        {
            return Ok(attempt);
        }
    }
    Ok(out.to_path_buf())
}

/// The first clean run of up to three; every attempt is kept.
fn f1_price(
    b: &Bench,
    cand: &Path,
    doc: &Path,
    out: &Path,
    log_dir: &Path,
    extra: &[String],
) -> Result<J, String> {
    let mut rep = J::Null;
    for k in 0..=runs::RETRIES {
        let attempt = runs::attempt(out, k);
        if !attempt.exists() {
            if let Some(parent) = attempt.parent() {
                std::fs::create_dir_all(parent)
                    .map_err(|e| format!("{}: {e}", parent.display()))?;
            }
            let mut cmd = vec![
                cand.to_string_lossy().into_owned(),
                "price".into(),
                "--params".into(),
                doc.to_string_lossy().into_owned(),
                "--json".into(),
                "--out".into(),
                attempt.to_string_lossy().into_owned(),
            ];
            cmd.extend(extra.iter().cloned());
            b.launch(&cmd, &attempt, log_dir)?;
        }
        rep = runs::load(&attempt);
        if runs::clean(&attempt)? && rep.get("status").and_then(J::as_str) == Some("extrapolated") {
            return Ok(rep);
        }
    }
    Ok(rep)
}

fn note(rep: &J) -> String {
    rep.get("extrapolated")
        .and_then(|e| e.get("ic"))
        .and_then(|e| e.get("cold"))
        .and_then(|e| e.get("units"))
        .and_then(|u| u.get("median"))
        .and_then(J::as_f64)
        .map_or(String::new(), |m| format!("median {m:.4e} units"))
}

pub fn f1(r: &bround::Run, b: &Bench, cand: &Path) -> Result<J, String> {
    let d = r.ctx.runs.join("f1");
    for k in 1..=F1_ROUNDS {
        for row in speed::m1_rows(&r.ctx)? {
            let out = d
                .join(format!("r{k}"))
                .join(format!("{}.price.json", row.id));
            let doc = f1_document(r, cand, &row)?;
            let rep = f1_price(b, cand, &doc, &out, &d, &[])?;
            let status = rep.get("status").and_then(J::as_str).unwrap_or("None");
            println!("r{k} {} F1: {status} {}", row.id, note(&rep));
        }
    }
    Ok(J::Null)
}

/// `M1`'s rows at the six largest sizes.
fn largest_rows(c: &rounds::Ctx) -> Result<Vec<Row>, String> {
    let m1 = speed::m1_rows(c)?;
    let mut sizes: Vec<(u8, u32, f64)> = Vec::new();
    for row in &m1 {
        if !sizes.iter().any(|s| (s.0, s.1) == (row.a, row.n)) {
            sizes.push((row.a, row.n, row.r));
        }
    }
    sizes.sort_by(|x, y| x.2.partial_cmp(&y.2).expect("no NaN"));
    let keep: Vec<(u8, u32)> = sizes
        .iter()
        .rev()
        .take(LARGEST)
        .map(|s| (s.0, s.1))
        .collect();
    Ok(m1
        .into_iter()
        .filter(|r| keep.contains(&(r.a, r.n)))
        .collect())
}

pub fn partial(r: &bround::Run, b: &Bench, cand: &Path) -> Result<J, String> {
    let d = r.ctx.runs.join("partial");
    for g in PARTIAL_FRACTIONS {
        for row in largest_rows(&r.ctx)? {
            let out = d
                .join(fraction_tag(g))
                .join(format!("{}.price.json", row.id));
            let doc = f1_document(r, cand, &row)?;
            let extra = ["--f1-partial".to_string(), json::py_float(g)];
            let rep = f1_price(b, cand, &doc, &out, &d, &extra)?;
            let status = rep.get("status").and_then(J::as_str).unwrap_or("None");
            println!(
                "g={} {}: {status} {}",
                json::py_float(g),
                row.id,
                note(&rep)
            );
        }
    }
    Ok(J::Null)
}

// ── the analysis ───────────────────────────────────────────────────

fn probes(mu: f64, granule: f64) -> f64 {
    granule / -(-granule / mu).exp_m1()
}

fn num(v: &J, path: &[&str]) -> Option<f64> {
    let mut cur = v;
    for k in path {
        cur = cur.get(k)?;
    }
    cur.as_f64()
}

fn f0_cold_units(rep: &J) -> Option<f64> {
    let reps = rep.get("repetitions")?.as_arr()?;
    let mut xs = Vec::new();
    for r in reps {
        xs.push(
            (num(r, &["setup_ns"])? + num(r, &["ic_online", "wall_ns"])?) / num(r, &["unit_ns"])?,
        );
    }
    (!xs.is_empty()).then(|| stats::median(&xs))
}

fn wall(out: &Path) -> Result<Option<f64>, String> {
    Ok(runs::run_record(out)?.and_then(|r| r.get("wall_seconds").and_then(J::as_f64)))
}

/// The chain's B7a-arm processes for the row: measurement 4's under the
/// chain.
fn f0_reports(chain: &Path, row: &Row) -> Result<Vec<(PathBuf, J)>, String> {
    let mut out = Vec::new();
    for k in 1..=bround::ROUNDS {
        let p = runs::figure_path(
            &chain
                .join("B7a")
                .join(&row.id)
                .join(format!("r{k}.price.json")),
        )?;
        let rep = runs::load(&p);
        if runs::clean(&p)? && rep.get("status").and_then(J::as_str) == Some("complete") {
            out.push((p, rep));
        }
    }
    Ok(out)
}

fn f1_reports(runs_dir: &Path, row: &Row) -> Result<Vec<(PathBuf, J)>, String> {
    let mut out = Vec::new();
    for k in 1..=F1_ROUNDS {
        let p = figure(
            &runs_dir
                .join("f1")
                .join(format!("r{k}"))
                .join(format!("{}.price.json", row.id)),
        )?;
        let rep = runs::load(&p);
        if runs::clean(&p)? && rep.get("status").and_then(J::as_str) == Some("extrapolated") {
            out.push((p, rep));
        }
    }
    Ok(out)
}

fn med(values: &[Option<f64>]) -> Option<f64> {
    let v: Vec<f64> = values
        .iter()
        .flatten()
        .copied()
        .filter(|x| x.is_finite())
        .collect();
    (!v.is_empty()).then(|| stats::median(&v))
}

fn opt_f(x: Option<f64>) -> J {
    x.map_or(J::Null, J::Float)
}

fn measurement5(
    r: &bround::Run,
    chain: &Path,
    sizes: &[((u8, u32), Vec<Row>)],
) -> Result<Vec<J>, String> {
    let mut out = Vec::new();
    for ((a, n), rows) in sizes {
        let (mut f0, mut f1) = (Vec::new(), Vec::new());
        for row in rows {
            f0.extend(f0_reports(chain, row)?);
            f1.extend(f1_reports(&r.ctx.runs, row)?);
        }
        let f0_cold = med(&f0
            .iter()
            .map(|(_, rep)| f0_cold_units(rep))
            .collect::<Vec<_>>());
        let key = |k: &str| {
            med(&f1
                .iter()
                .map(|(_, rep)| num(rep, &["extrapolated", "ic", "cold", "units", k]))
                .collect::<Vec<_>>())
        };
        let keys = [
            "median",
            "expectation_lo",
            "expectation_hi",
            "predictive_lo",
            "predictive_hi",
        ];
        let bounds: Vec<(String, Option<f64>)> =
            keys.iter().map(|k| (k.to_string(), key(k))).collect();
        let get = |k: &str| bounds.iter().find(|(b, _)| b == k).and_then(|(_, v)| *v);
        let ratio = match (get("median"), f0_cold) {
            (Some(m), Some(f)) if f != 0.0 => Some(m / f),
            _ => None,
        };
        let inside = match (f0_cold, get("predictive_lo"), get("predictive_hi")) {
            (Some(f), Some(lo), Some(hi)) => lo <= f && f <= hi,
            _ => false,
        };
        let mut walls_f0 = Vec::new();
        for (p, _) in &f0 {
            walls_f0.push(wall(p)?);
        }
        let mut walls_f1 = Vec::new();
        for (p, _) in &f1 {
            walls_f1.push(wall(p)?);
        }
        let wall_ratio = match (med(&walls_f1), med(&walls_f0)) {
            (Some(x), Some(y)) if y != 0.0 => Some(x / y),
            _ => None,
        };
        let tight = if TIGHT_SIZES.contains(n) {
            ratio.map_or(J::Null, |q| J::Bool(TIGHT_BAND.0 <= q && q <= TIGHT_BAND.1))
        } else {
            J::Null
        };
        out.push(obj([
            ("slug", J::Str(suite::curve_slug(*a, *n)?)),
            ("a", J::Int((*a).into())),
            ("n", J::Int((*n).into())),
            ("log2_r", J::Float(rounds::round3(rows[0].r.log2()))),
            ("f0_processes", J::Int(f0.len() as i128)),
            ("f1_processes", J::Int(f1.len() as i128)),
            ("f0_cold_units", opt_f(f0_cold)),
            (
                "f1_cold_units",
                J::Obj(bounds.into_iter().map(|(k, v)| (k, opt_f(v))).collect()),
            ),
            ("ratio_f1_over_f0", opt_f(ratio)),
            ("f0_inside_predictive", J::Bool(inside)),
            ("wall_f1_over_f0", opt_f(wall_ratio)),
            ("in_tight_band", tight),
        ]));
    }
    Ok(out)
}

fn measurement6(
    r: &bround::Run,
    chain: &Path,
    sizes: &[((u8, u32), Vec<Row>)],
) -> Result<J, String> {
    let (mut rows_out, mut ns, mut kappas) = (Vec::new(), Vec::new(), Vec::new());
    for ((a, n), rows) in sizes {
        let (mut f0, mut f1) = (Vec::new(), Vec::new());
        for row in rows {
            f0.extend(f0_reports(chain, row)?);
            f1.extend(f1_reports(&r.ctx.runs, row)?);
        }
        let slug = suite::curve_slug(*a, *n)?;
        let Some((_, first)) = f1.first() else {
            rows_out.push(obj([
                ("slug", J::Str(slug)),
                ("n", J::Int((*n).into())),
                ("missing", J::Str("no F1 report".into())),
            ]));
            continue;
        };
        let collect = first.at("phases")?.at("ic")?.at("collect")?;
        let count = collect.at("count")?;
        let sample = collect.at("kappa_sample")?;
        let (mu, window) = (
            num(count, &["mu"]).ok_or("count.mu")?,
            num(count, &["window"]).ok_or("count.window")?,
        );
        let (mut rhos, mut ks) = (Vec::new(), Vec::new());
        for (_, rep) in &f0 {
            let pc = rep.at("counts")?.at("pass")?;
            rhos.push(
                num(pc, &["logs", "relations_accepted"])
                    .zip(num(pc, &["select", "columns"]))
                    .map(|(x, y)| x / y),
            );
            let found = num(pc, &["collect", "relations_collected"]).unwrap_or(0.0);
            if found != 0.0 {
                ks.push(
                    num(pc, &["collect", "summands_scanned"])
                        .map(|s| s / found / probes(mu, window)),
                );
            }
        }
        let kappa_f0 = med(&ks);
        if let Some(k) = kappa_f0 {
            ns.push(f64::from(*n));
            kappas.push(k);
        }
        rows_out.push(obj([
            ("slug", J::Str(slug)),
            ("n", J::Int((*n).into())),
            ("kappa_sample", sample.at("kappa")?.clone()),
            ("kappa_sample_interval", sample.at("interval")?.clone()),
            ("relations_found", sample.at("relations_found")?.clone()),
            ("f0_rho", opt_f(med(&rhos))),
            ("f0_kappa", opt_f(kappa_f0)),
            ("carried", collect.at("carried")?.clone()),
        ]));
    }
    Ok(obj([
        ("sizes", J::Arr(rows_out)),
        (
            "fit_f0_kappa_on_n",
            stats::fit(&ns, &kappas, 0.95).unwrap_or(J::Null),
        ),
    ]))
}

fn measurement7(r: &bround::Run, sizes: &[((u8, u32), Vec<Row>)]) -> Result<Vec<J>, String> {
    let mut out = Vec::new();
    let largest: Vec<(u8, u32)> = largest_rows(&r.ctx)?.iter().map(|x| (x.a, x.n)).collect();
    let per_pair = |ic: &J| -> Option<f64> {
        Some(num(ic, &["build", "ns", "build"])? / num(ic, &["build", "stored_pairs"])?)
    };
    let descent = |ic: &J| -> Option<f64> {
        let s = ic.get("descent")?.get("sample")?;
        let xs: Vec<f64> = s
            .get("ns_per_probe")
            .or_else(|| s.get("ns_per_summand"))?
            .as_arr()?
            .iter()
            .filter_map(J::as_f64)
            .collect();
        (!xs.is_empty()).then(|| stats::median(&xs))
    };
    for ((a, n), rows) in sizes {
        if !largest.contains(&(*a, *n)) {
            continue;
        }
        for g in PARTIAL_FRACTIONS {
            let names = [
                "build_per_pair",
                "scan_per_summand",
                "descent_per_probe",
                "keys_over_g2_full",
            ];
            let mut ratios: Vec<Vec<Option<f64>>> = vec![Vec::new(); 4];
            for row in rows {
                let full = f1_reports(&r.ctx.runs, row)?;
                let p = figure(
                    &r.ctx
                        .runs
                        .join("partial")
                        .join(fraction_tag(g))
                        .join(format!("{}.price.json", row.id)),
                )?;
                let part = runs::load(&p);
                let Some((_, full)) = full.first() else {
                    continue;
                };
                if !(runs::clean(&p)?
                    && part.get("status").and_then(J::as_str) == Some("extrapolated"))
                {
                    continue;
                }
                let (fic, pic) = (full.at("phases")?.at("ic")?, part.at("phases")?.at("ic")?);
                let div = |x: Option<f64>, y: Option<f64>| x.zip(y).map(|(x, y)| x / y);
                ratios[0].push(div(per_pair(pic), per_pair(fic)));
                ratios[1].push(div(
                    num(pic, &["collect", "sample", "ns_per_summand"]),
                    num(fic, &["collect", "sample", "ns_per_summand"]),
                ));
                ratios[2].push(div(descent(pic), descent(fic)));
                let fraction = num(pic, &["build", "partial_table", "fraction"]);
                ratios[3].push(
                    num(pic, &["build", "distinct_keys"])
                        .zip(fraction)
                        .zip(num(fic, &["build", "distinct_keys"]))
                        .map(|((pk, f), fk)| pk / (f * f * fk)),
                );
            }
            let mut kv = vec![
                ("slug".to_string(), J::Str(suite::curve_slug(*a, *n)?)),
                ("n".into(), J::Int((*n).into())),
                ("fraction_asked".into(), J::Float(g)),
            ];
            for (name, values) in names.iter().zip(&ratios) {
                kv.push((name.to_string(), opt_f(med(values))));
            }
            kv.push(("rows".into(), J::Int(ratios[0].len() as i128)));
            out.push(J::Obj(kv));
        }
    }
    Ok(out)
}

/// B7a's decision: measurements 2–4 from `bround` (its own run tree for
/// conformance, the pin and the translation; the chain for the timing),
/// and 5–7 here.
pub fn analyse(r: &bround::Run, chain: &bround::Run) -> Result<J, String> {
    let J::Obj(mut out) = bround::analyse(r)? else {
        return Err("bround's analysis is not an object".into());
    };
    let chain_doc = bround::analyse(chain)?;
    let timing = chain_doc
        .get("chain")
        .and_then(|c| c.get("steps"))
        .and_then(|s| s.get("B7a"))
        .cloned()
        .unwrap_or(J::Null);
    out.push(("timing_in_the_chain".into(), timing.clone()));
    let mut sizes = suite::by_size(&speed::m1_rows(&r.ctx)?);
    sizes.sort_by(|x, y| x.1[0].r.partial_cmp(&y.1[0].r).expect("no NaN"));
    let m5 = measurement5(r, &chain.ctx.runs.join("chain"), &sizes)?;
    let f1_dir = r.ctx.runs.join("f1");
    out.push((
        "f1_against_f0".into(),
        obj([
            (
                "what",
                J::Str(
                    "per size: F0's cold units (the chain's B7a arm, median of its processes) \
                     against F1's extrapolated cold units (median of its processes' figures)"
                        .into(),
                ),
            ),
            (
                "accounting",
                if f1_dir.exists() {
                    runs::accounting(&f1_dir)?
                } else {
                    J::Null
                },
            ),
            ("sizes", J::Arr(m5.clone())),
        ]),
    ));
    out.push((
        "count_check".into(),
        measurement6(r, &chain.ctx.runs.join("chain"), &sizes)?,
    ));
    out.push((
        "partial_tables".into(),
        if r.ctx.runs.join("partial").exists() {
            J::Arr(measurement7(r, &sizes)?)
        } else {
            J::Null
        },
    ));
    let truthy = |s: &J, k: &str| s.get(k).is_some_and(J::truthy);
    let inside = m5
        .iter()
        .filter(|s| truthy(s, "f0_inside_predictive"))
        .count();
    let tight: Vec<&J> = m5
        .iter()
        .filter(|s| {
            s.get("n")
                .and_then(J::as_i128)
                .is_some_and(|n| TIGHT_SIZES.iter().any(|t| i128::from(*t) == n))
        })
        .collect();
    let met = inside >= INSIDE_NEEDED
        && tight.len() == TIGHT_SIZES.len()
        && tight.iter().all(|s| truthy(s, "in_tight_band"));
    out.push((
        "falsification_target".into(),
        obj([
            ("inside_predictive", J::Int(inside as i128)),
            ("inside_needed", J::Int(INSIDE_NEEDED as i128)),
            ("sizes", J::Int(m5.len() as i128)),
            (
                "tight_band",
                J::Arr(vec![J::Float(TIGHT_BAND.0), J::Float(TIGHT_BAND.1)]),
            ),
            (
                "tight",
                J::Obj(
                    tight
                        .iter()
                        .map(|s| {
                            (
                                s.get("slug").and_then(J::as_str).unwrap_or("").to_string(),
                                s.get("ratio_f1_over_f0").cloned().unwrap_or(J::Null),
                            )
                        })
                        .collect(),
                ),
            ),
            ("met", J::Bool(met)),
        ]),
    ));
    let incomplete: Vec<J> = m5
        .iter()
        .filter(|s| {
            s.get("f0_processes").and_then(J::as_i128) == Some(0)
                || s.get("f1_processes").and_then(J::as_i128) == Some(0)
        })
        .filter_map(|s| s.get("slug").cloned())
        .collect();
    let decision = if !incomplete.is_empty() {
        obj([
            ("accepted", J::Bool(false)),
            ("complete", J::Bool(false)),
            (
                "reasons",
                J::Arr(vec![J::Str(format!(
                    "measurement 5 has no F0 or no F1 figure at {}",
                    json::dumps_line(&J::Arr(incomplete), false)
                ))]),
            ),
        ])
    } else {
        let mut reasons = Vec::new();
        let doc = J::Obj(out.clone());
        if let Some(failed) = doc
            .get("conformance")
            .and_then(|c| c.get("candidate"))
            .and_then(|c| c.get("failed"))
            .and_then(J::as_arr)
            .filter(|f| !f.is_empty())
        {
            reasons.push(J::Str(format!(
                "conformance cases failed: {}",
                json::dumps_line(&J::Arr(failed.to_vec()), false)
            )));
        }
        for name in ["pin", "translate"] {
            if let Some(section) = doc.get(name) {
                if !section.get("held").is_some_and(J::truthy) {
                    reasons.push(J::Str(format!(
                        "the {name} check failed on {}",
                        json::dumps_line(section.get("failing_rows").unwrap_or(&J::Null), false)
                    )));
                }
            }
        }
        if timing.get("any_regression").is_some_and(J::truthy) {
            reasons.push(J::Str("a size regresses beyond its A/A band at F0".into()));
        }
        if !met {
            reasons.push(J::Str("the falsification target is missed".into()));
        }
        obj([
            ("accepted", J::Bool(reasons.is_empty())),
            ("complete", J::Bool(true)),
            ("reasons", J::Arr(reasons)),
            (
                "note",
                J::Str(
                    "the tests (measurement 1) are run with cargo and recorded beside this analysis"
                        .into(),
                ),
            ),
        ])
    };
    out.push(("decision".into(), decision));
    Ok(J::Obj(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probes_is_the_expected_count_to_a_window() {
        // Far below the window the count is the window; far above it,
        // μ + w/2; at μ = w it is w / (1 − e^{−1}).
        assert!((probes(1e-9, 16.0) - 16.0).abs() < 1e-9);
        assert!((probes(1e9, 16.0) - (1e9 + 8.0)).abs() < 1e-3);
        assert!((probes(16.0, 16.0) - 16.0 / (1.0 - (-1.0f64).exp())).abs() < 1e-12);
    }

    #[test]
    fn the_fit_recovers_a_line_and_needs_three_points() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0];
        let ys: Vec<f64> = xs.iter().map(|x| 2.0 * x + 1.0).collect();
        let f = stats::fit(&xs, &ys, 0.95).unwrap();
        assert!((f.get("beta").and_then(J::as_f64).unwrap() - 2.0).abs() < 1e-12);
        assert!((f.get("alpha").and_then(J::as_f64).unwrap() - 1.0).abs() < 1e-12);
        assert!(stats::fit(&xs[..2], &ys[..2], 0.95).is_none());
    }

    #[test]
    fn the_fractions_name_their_directories_as_python_did() {
        assert_eq!(fraction_tag(0.5), "g0.5");
        assert_eq!(fraction_tag(0.25), "g0.25");
    }
}
