//! B3b's measurement 6, the two-word premium, natively: a stage
//! diagnostic, not a speed figure.  At `n = 47`, `53` and `61`, on `M1`'s
//! rows, three rounds, B3b's arm prices each row with `--kic-wide` and
//! without, the order alternating by round, isolated.  Each pipeline's own
//! reported costs give three prices a process:
//! - the build per stored pair: the `build` phase over `stored_pairs`;
//! - the scan per summand: the `collect` phase over `summands_scanned`;
//! - the descent per probe: the online interval's `target_pdp` phase (the
//!   descent's probe loop) over the descent's `trials_total`.
//!
//! Each price is the median over the process's repetitions, and each
//! size's premium is the median over its pairs of the two-word price over
//! the one-word one.  B7b reads F1 at `n = 83` with it.

use std::path::Path;

use super::bench::Bench;
use super::bround;
use super::json::{obj, J};
use super::rounds::{self, speed};
use super::runs;
use super::stats;
use super::suite::{self, Row};

const SIZES: [u32; 3] = [47, 53, 61];
const ROUNDS: u32 = 3;
const ARMS: [&str; 2] = ["one-word", "two-word"];

fn rows(r: &bround::Run) -> Result<Vec<Row>, String> {
    Ok(speed::m1_rows(&r.ctx)?
        .into_iter()
        .filter(|row| SIZES.contains(&row.n))
        .collect())
}

fn extra(arm: &str) -> Vec<String> {
    if arm == "two-word" {
        vec!["--kic-wide".into()]
    } else {
        Vec::new()
    }
}

pub fn run(r: &bround::Run, b: &Bench, cand: &Path) -> Result<J, String> {
    let d = r.ctx.runs.join("premium");
    for k in 1..=ROUNDS {
        let order: Vec<&str> = if k % 2 == 1 {
            ARMS.to_vec()
        } else {
            ARMS.iter().rev().copied().collect()
        };
        for row in rows(r)? {
            for arm in &order {
                let out = d.join(arm).join(&row.id).join(format!("r{k}.price.json"));
                let rep = b.price(cand, &row, &out, &d, &extra(arm))?;
                let status = rep.get("status").and_then(J::as_str).unwrap_or("None");
                println!("r{k} {} {arm}: {status}", row.id);
            }
        }
    }
    Ok(J::Null)
}

fn num(v: &J, path: &[&str]) -> Option<f64> {
    let mut cur = v;
    for k in path {
        cur = cur.get(k)?;
    }
    cur.as_f64()
}

/// A process's three prices, each the median over its repetitions.
fn prices(rep: &J) -> Option<[f64; 3]> {
    let pass = rep.get("counts")?.get("pass")?;
    let pairs = num(pass, &["build", "stored_pairs"])?;
    let summands = num(pass, &["collect", "summands_scanned"])?;
    let probes = num(pass, &["descent", "trials_total"])?;
    let reps = rep.get("repetitions")?.as_arr()?;
    let mut cols = [Vec::new(), Vec::new(), Vec::new()];
    for x in reps {
        cols[0].push(num(x, &["setup_phases_ns", "build"])? / pairs);
        cols[1].push(num(x, &["setup_phases_ns", "collect"])? / summands);
        cols[2].push(num(x, &["ic_online", "phases_ns", "target_pdp"])? / probes);
    }
    if cols[0].is_empty() {
        return None;
    }
    Some([
        stats::median(&cols[0]),
        stats::median(&cols[1]),
        stats::median(&cols[2]),
    ])
}

pub fn analyse(r: &bround::Run) -> Result<J, String> {
    let d = r.ctx.runs.join("premium");
    let names = ["build_per_pair", "scan_per_summand", "descent_per_probe"];
    let mut out = Vec::new();
    for ((a, n), rs) in suite::by_size(&rows(r)?) {
        let mut ratios: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        let mut missing = 0;
        for row in &rs {
            for k in 1..=ROUNDS {
                let load = |arm: &str| -> Result<Option<J>, String> {
                    runs::figure(&d.join(arm).join(&row.id).join(format!("r{k}.price.json")))
                };
                match (
                    load(ARMS[0])?.as_ref().and_then(prices),
                    load(ARMS[1])?.as_ref().and_then(prices),
                ) {
                    (Some(one), Some(two)) => {
                        for i in 0..3 {
                            ratios[i].push(two[i] / one[i]);
                        }
                    }
                    _ => missing += 1,
                }
            }
        }
        let mut kv = rounds::size_head(&suite::curve_slug(a, n)?, a, n, &rs);
        for (name, values) in names.iter().zip(&ratios) {
            kv.push((
                format!("two_word_over_one_word_{name}"),
                if values.is_empty() {
                    J::Null
                } else {
                    J::Float(stats::median(values))
                },
            ));
        }
        kv.push(("pairs".into(), J::Int(ratios[0].len() as i128)));
        kv.push(("missing_pairs".into(), J::Int(missing)));
        out.push(J::Obj(kv));
    }
    Ok(obj([
        (
            "what",
            J::Str(
                "the two-word pipeline's price over the one-word pipeline's, per stage, \
                 median over each size's pairs: a stage diagnostic, not a speed figure"
                    .into(),
            ),
        ),
        (
            "accounting",
            if d.exists() {
                runs::accounting(&d)?
            } else {
                J::Null
            },
        ),
        ("sizes", J::Arr(out)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::super::json;
    use super::*;

    #[test]
    fn prices_read_each_stage_over_its_count() {
        let rep = json::parse(
            r#"{"counts": {"pass": {"build": {"stored_pairs": 100},
                 "collect": {"summands_scanned": 1000}, "descent": {"trials_total": 10}}},
                "repetitions": [
                 {"setup_phases_ns": {"build": 300, "collect": 9000}, "ic_online": {"phases_ns": {"target_pdp": 50}}},
                 {"setup_phases_ns": {"build": 200, "collect": 8000}, "ic_online": {"phases_ns": {"target_pdp": 40}}},
                 {"setup_phases_ns": {"build": 400, "collect": 7000}, "ic_online": {"phases_ns": {"target_pdp": 60}}}]}"#,
        )
        .unwrap();
        assert_eq!(prices(&rep), Some([3.0, 8.0, 5.0]));
        assert_eq!(extra("two-word"), vec!["--kic-wide".to_string()]);
        assert!(extra("one-word").is_empty());
    }
}
