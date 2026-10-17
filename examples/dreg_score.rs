//! Score, compare and reference refutation-degree runs (`dreg_ladder` output).
//!
//! ```sh
//! # per-cell values and the registered reading, with the semi-regular reference
//! cargo run --release --example dreg_score -- score runs/*.jsonl
//! # draw-by-draw identity of a rerun against committed rows (all but `secs`)
//! cargo run --release --example dreg_score -- compare committed.jsonl rerun.jsonl
//! # the semi-regular degree of regularity of chained-system shapes
//! cargo run --release --example dreg_score -- semireg 4 5:2 6:2 7:2
//! # upper bounds on the degree-D Macaulay matrix shape
//! cargo run --release --example dreg_score -- shape 4 6:2:8 8:2:7
//! ```
//!
//! The native replacement for the degree-7 study's `score.py`, `compare.py`
//! and `semireg.py`, which stay in the repository as that study's record.
//!
//! **Reading of a cell** (the rule the refutation-degree studies register):
//! a resolved draw is its degree, exact; a full matrix without resolution
//! at `d_max` is a bound `≥ d_max + 1` and counts as that value; caps-hit
//! draws are excluded.  With at least three values, the median's middle
//! values give the reading: one degree `D` (all exact), `≥D` (all bounds),
//! or `a / b` when they differ.  Fewer than three values: not testable.

use crypto_lib::cryptanalysis::koblitz_bench::semi_regular_dreg;
use serde_json::Value;
use std::collections::BTreeMap;

/// The equation degrees of the `m`-summand chained system over GF(2^n):
/// `(m − 2)·n` cubic links' equations and `n` quadratic ones (the last link
/// closes on the known `x(R)`).
fn chain_degrees(m: usize, n: usize) -> Vec<u32> {
    std::iter::repeat_n(3, (m - 2) * n)
        .chain(std::iter::repeat_n(2, n))
        .collect()
}

fn chain_dreg(m: usize, n: usize, ell: usize) -> u32 {
    semi_regular_dreg((m - 2) * n + m * ell, &chain_degrees(m, n)).expect("reference exists")
}

/// (degree, exact) of one measured draw; `None` for caps-hit and satisfiable.
fn value(outcome: &Value) -> Option<(u64, bool)> {
    let degree = outcome.get("degree")?.as_u64()?;
    match outcome.get("kind")?.as_str()? {
        "resolved" => Some((degree, true)),
        "at_least" => Some((degree, false)),
        _ => None,
    }
}

fn show(v: (u64, bool)) -> String {
    if v.1 {
        v.0.to_string()
    } else {
        format!("≥{}", v.0)
    }
}

fn reading(values: &[(u64, bool)]) -> String {
    if values.len() < 3 {
        return "not testable".into();
    }
    let mut sorted = values.to_vec();
    sorted.sort();
    let k = sorted.len();
    let middle: Vec<(u64, bool)> = if k % 2 == 1 {
        vec![sorted[k / 2]]
    } else {
        vec![sorted[k / 2 - 1], sorted[k / 2]]
    };
    if middle.iter().all(|v| *v == middle[0]) {
        show(middle[0])
    } else {
        format!("{} / {}", show(middle[0]), show(middle[1]))
    }
}

fn lines(path: &str) -> Vec<Value> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{path}: {e}"))
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("{path}: {e}")))
        .filter(|v: &Value| v.get("control").is_none())
        .collect()
}

fn score(paths: &[String]) {
    // (m, n, ℓ) -> per-draw (unsat_index or order, value text, value, ffd)
    let mut cells: BTreeMap<(u64, u64, u64), Vec<(String, Option<(u64, bool)>, String)>> =
        BTreeMap::new();
    for path in paths {
        for v in lines(path) {
            let outcome = &v["outcome"];
            if outcome["kind"] == "satisfiable" {
                continue;
            }
            let key = (
                v.get("m").and_then(Value::as_u64).unwrap_or(3),
                v["n"].as_u64().expect("n"),
                v["ell"].as_u64().expect("ell"),
            );
            let val = value(outcome);
            let shown = val.map_or_else(|| outcome.to_string(), show);
            let ffd = v["ffd"].as_u64().map_or("-".into(), |f| f.to_string());
            cells.entry(key).or_default().push((shown, val, ffd));
        }
    }
    println!("| m | cell (n, ℓ) | N | S | values | reading | semi-regular D_reg | FFD |");
    println!("|--:|---|--:|--:|---|---|--:|---|");
    for ((m, n, ell), draws) in &cells {
        let (m, n, ell) = (*m as usize, *n as usize, *ell as usize);
        let vals: Vec<(u64, bool)> = draws.iter().filter_map(|d| d.1).collect();
        println!(
            "| {m} | ({n}, {ell}) | {} | {:+} | {} | {} | {} | {} |",
            (m - 2) * n + m * ell,
            n as i64 - (m * ell) as i64,
            draws
                .iter()
                .map(|d| d.0.as_str())
                .collect::<Vec<_>>()
                .join(" "),
            reading(&vals),
            chain_dreg(m, n, ell),
            draws
                .iter()
                .map(|d| d.2.as_str())
                .collect::<Vec<_>>()
                .join(" "),
        );
    }
}

fn compare(reference: &str, rerun: &str) {
    let key = |v: &Value| {
        [
            "draw",
            "v_basis",
            "x_r",
            "solutions",
            "outcome",
            "ffd",
            "n",
            "ell",
        ]
        .iter()
        .map(|k| v.get(*k).cloned().unwrap_or(Value::Null).to_string())
        .collect::<Vec<_>>()
    };
    let (want, got) = (lines(reference), lines(rerun));
    let mut mismatches = want.len().abs_diff(got.len());
    let mut measured = 0;
    for (a, b) in want.iter().zip(&got) {
        if key(a) != key(b) {
            mismatches += 1;
        }
        if a["outcome"]["kind"] != "satisfiable" {
            measured += 1;
        }
    }
    println!(
        "{rerun}: {} rows, {measured} measured, mismatches {mismatches}",
        want.len()
    );
    if mismatches > 0 {
        std::process::exit(1);
    }
}

fn semireg(m: usize, cells: &[String]) {
    println!("| m | cell (n, ℓ) | N | S | equations | semi-regular D_reg |");
    println!("|--:|---|--:|--:|---|--:|");
    for c in cells {
        let (n, ell) = c.split_once(':').expect("cell is n:ell");
        let (n, ell): (usize, usize) = (n.parse().expect("n"), ell.parse().expect("ell"));
        println!(
            "| {m} | ({n}, {ell}) | {} | {:+} | {} cubic + {n} quadratic | {} |",
            (m - 2) * n + m * ell,
            n as i64 - (m * ell) as i64,
            (m - 2) * n,
            chain_dreg(m, n, ell)
        );
    }
}

/// `C(N, ≤ d)`: square-free monomials of degree at most `d` in `N` unknowns.
fn monomials_up_to(n_vars: u128, d: i64) -> u128 {
    let mut total = 0u128;
    let mut binom = 1u128;
    for i in 0..=d.max(-1) {
        if i as u128 > n_vars {
            break;
        }
        total += binom;
        binom = binom * (n_vars - i as u128) / (i as u128 + 1);
    }
    total
}

/// Upper bounds on the degree-`D` Macaulay matrix of the `m`-summand chain:
/// every equation times every monomial of degree `≤ D − deg`, and every
/// monomial of degree `≤ D` as a column.
fn shape(m: usize, cells: &[String]) {
    println!("| m | cell (n, ℓ) | N | degree | rows | columns |");
    println!("|--:|---|--:|--:|--:|--:|");
    for c in cells {
        let p: Vec<usize> = c.split(':').map(|x| x.parse().expect("n:ℓ:D")).collect();
        let (n, ell, d) = (p[0], p[1], p[2] as i64);
        let n_vars = ((m - 2) * n + m * ell) as u128;
        let rows = ((m - 2) * n) as u128 * monomials_up_to(n_vars, d - 3)
            + n as u128 * monomials_up_to(n_vars, d - 2);
        println!(
            "| {m} | ({n}, {ell}) | {n_vars} | {d} | {rows} | {} |",
            monomials_up_to(n_vars, d)
        );
    }
}

/// Score `dreg_closure` output against PREREGISTRATION.md of
/// `research/dreg_closure_mutantxl_20261016`: per draw the closure degree
/// and its gap to the committed XL degree; per cell the Q1 verdict
/// (below / not below / mixed), or "exploratory" when `C(N, ≤ XL − 1)`
/// exceeds the column cap (the registered testability rule); then Q1's
/// study verdict, Q2's gap counts and Q3's per-ℓ medians at m = 3, and the
/// registered descriptive cost table.  `--satisfiable` output files are
/// tallied as the soundness control.
fn closure(paths: &[String]) {
    struct Draw {
        xl: u64,
        known: bool,
        closure: Option<u64>,
        measured_to: u64,
        // At the closure's refuting degree: columns, products, word ops.
        cost: Option<(u64, u64, u64)>,
    }
    let mut cells: BTreeMap<(u64, u64, u64), Vec<Draw>> = BTreeMap::new();
    let mut max_cols = 0u128;
    let (mut sat_draws, mut sat_unsound) = (0usize, 0usize);
    let mut sat_degree: BTreeMap<u64, usize> = BTreeMap::new();
    for path in paths {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let v: Value = serde_json::from_str(line).expect("json");
            if let Some(bad) = v["refuted_satisfiable"].as_bool() {
                sat_draws += 1;
                sat_unsound += bad as usize;
                *sat_degree
                    .entry(v["measured_to"].as_u64().unwrap())
                    .or_default() += 1;
                continue;
            }
            max_cols = max_cols.max(v["max_cols"].as_f64().unwrap() as u128);
            let key = (
                v["m"].as_u64().unwrap(),
                v["n"].as_u64().unwrap(),
                v["ell"].as_u64().unwrap(),
            );
            let o = &v["committed_outcome"];
            let closure = v["closure_degree"].as_u64();
            let cost = closure.and_then(|c| {
                v["scan"].as_array().unwrap().iter().find_map(|s| {
                    (s["degree"].as_u64() == Some(c)).then(|| {
                        (
                            s["columns"].as_u64().unwrap(),
                            s["closure_products"].as_u64().unwrap(),
                            s["closure_word_ops"].as_u64().unwrap(),
                        )
                    })
                })
            });
            cells.entry(key).or_default().push(Draw {
                xl: o["degree"].as_u64().unwrap(),
                known: o["kind"] == "resolved",
                closure,
                measured_to: v["measured_to"].as_u64().unwrap(),
                cost,
            });
        }
    }
    let med = |mut v: Vec<f64>| -> f64 {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let k = v.len();
        if k == 0 {
            f64::NAN
        } else if k % 2 == 1 {
            v[k / 2]
        } else {
            (v[k / 2 - 1] + v[k / 2]) / 2.0
        }
    };
    println!("| m | cell (n, ℓ) | N | committed XL | closure degree | gap XL − closure | Q1 |");
    println!("|--:|---|--:|---|---|---|---|");
    let (mut below, mut not_below, mut testable) = (0, 0, 0);
    let (mut gap1, mut gap2) = (0, 0);
    let mut by_ell: BTreeMap<u64, Vec<f64>> = BTreeMap::new();
    let mut cost_rows = Vec::new();
    for ((m, n, ell), draws) in &cells {
        let n_vars = (m - 2) * n + m * ell;
        // Registered rule: testable when the degree XL − 1 fits the cap.
        let xl_min = draws.iter().map(|d| d.xl).min().unwrap();
        let fits = monomials_up_to(n_vars as u128, xl_min as i64 - 1) <= max_cols;
        let k = draws
            .iter()
            .filter(|d| matches!(d.closure, Some(c) if c < d.xl))
            .count();
        let unmeasured = draws
            .iter()
            .any(|d| d.closure.is_none() && d.measured_to + 1 < d.xl - (!d.known) as u64);
        let verdict = if !fits {
            if k == draws.len() {
                "exploratory: below".to_string()
            } else {
                "exploratory".into()
            }
        } else if unmeasured {
            "untestable".into()
        } else {
            testable += 1;
            if k == draws.len() {
                below += 1;
                "below".into()
            } else if k == 0 {
                not_below += 1;
                "not below".into()
            } else {
                "mixed".into()
            }
        };
        for d in draws {
            if let Some(c) = d.closure {
                if d.known && c < d.xl {
                    if d.xl - c == 1 {
                        gap1 += 1
                    } else {
                        gap2 += 1
                    }
                }
                if *m == 3 {
                    by_ell.entry(*ell).or_default().push(c as f64);
                }
            }
        }
        let xl: Vec<String> = draws
            .iter()
            .map(|d| {
                if d.known {
                    d.xl.to_string()
                } else {
                    format!("≥{}", d.xl)
                }
            })
            .collect();
        let cl: Vec<String> = draws
            .iter()
            .map(|d| {
                d.closure
                    .map_or(format!(">{}", d.measured_to), |c| c.to_string())
            })
            .collect();
        let gaps: Vec<String> = draws
            .iter()
            .map(|d| match (d.known, d.closure) {
                (true, Some(c)) => (d.xl - c).to_string(),
                (false, Some(c)) => format!("≥{}", d.xl - c),
                _ => "–".into(),
            })
            .collect();
        println!(
            "| {m} | ({n}, {ell}) | {n_vars} | {} | {} | {} | {verdict} |",
            xl.join(" "),
            cl.join(" "),
            gaps.join(" ")
        );
        let lg = |d: u64| (monomials_up_to(n_vars as u128, d as i64) as f64).log2();
        let costed: Vec<&Draw> = draws.iter().filter(|d| d.cost.is_some()).collect();
        cost_rows.push(format!(
            "| {m} | ({n}, {ell}) | {n_vars} | {:.1} | {:.1} | {:.1} | {:.1} | {:.1} |",
            med(draws.iter().map(|d| d.xl as f64).collect()),
            med(draws.iter().map(|d| lg(d.xl)).collect()),
            med(costed.iter().map(|d| lg(d.closure.unwrap())).collect()),
            med(costed
                .iter()
                .map(|d| (d.cost.unwrap().1.max(1) as f64).log2())
                .collect()),
            med(costed
                .iter()
                .map(|d| (d.cost.unwrap().2.max(1) as f64).log2())
                .collect()),
        ));
    }
    let study = if 3 * below >= 2 * testable {
        "closure lowers"
    } else if 3 * not_below >= 2 * testable {
        "no effect"
    } else {
        "mixed"
    };
    println!();
    println!("Q1: {below} of {testable} testable cells below, {not_below} not below → {study}");
    println!("Q2: gaps among resolved-XL draws refuted below XL: {gap1} of 1, {gap2} of ≥ 2");
    let medians: Vec<String> = by_ell
        .iter()
        .map(|(ell, v)| format!("ℓ={ell}: {} (n={})", med(v.clone()), v.len()))
        .collect();
    println!(
        "Q3: m = 3 closure-degree medians by ℓ: {}",
        medians.join(", ")
    );
    if sat_draws > 0 {
        let to: Vec<String> = sat_degree
            .iter()
            .map(|(d, k)| format!("{k} to degree {d}"))
            .collect();
        println!(
            "Soundness control: {sat_draws} satisfiable draws ({}), {sat_unsound} refuted",
            to.join(", ")
        );
    }
    println!();
    println!("Descriptive cost, cell medians (log2; columns are C(N, ≤D)):");
    println!();
    println!("| m | cell (n, ℓ) | N | XL degree | log2 cols at XL | log2 cols at closure | log2 closure products | log2 closure word ops |");
    println!("|--:|---|--:|--:|--:|--:|--:|--:|");
    for r in cost_rows {
        println!("{r}");
    }
}

/// Score `research/dreg_closure_scale_20261016` (its PREREGISTRATION.md):
/// one `dreg_closure --to 4 --draw K` output file per draw.  An empty file
/// is a draw the memory or time limit killed (unmeasured, never evidence);
/// its cell is read from the file name (`fresh-m<M>-<n>-<ℓ>.d<K>` or
/// `committed-cell-<n>-<ℓ>-*`, which are `m = 3`).  Per draw: closure ≤ 4,
/// above 4 (degree 4 built and saturated without refuting), or unmeasured;
/// per cell flat / rises / unmeasured; the study verdict; Q5 (XL never
/// refutes at 4); Q6 (log₂ word ops ÷ log₂ columns at the refuting degree).
/// `--satisfiable` outputs are tallied as the soundness control.
fn scale(paths: &[String]) {
    #[derive(Default)]
    struct Cell {
        n_vars: u64,
        degrees: Vec<String>,
        flat: usize,
        above: usize,
        unmeasured: usize,
        xl_at_4: usize,
        ratios: Vec<f64>,
        secs: Vec<f64>,
    }
    let mut cells: BTreeMap<(u64, u64, u64), Cell> = BTreeMap::new();
    let (mut sat_draws, mut sat_bad, mut sat_unmeasured) = (0usize, 0usize, 0usize);
    for path in paths {
        let name = std::path::Path::new(path)
            .file_name()
            .and_then(|f| f.to_str())
            .unwrap_or("");
        let text = std::fs::read_to_string(path).unwrap_or_default();
        let line = text.lines().find(|l| !l.trim().is_empty());
        let Some(line) = line else {
            // A killed control is an unmeasured control, never evidence.
            if name.starts_with("control-") {
                sat_unmeasured += 1;
                continue;
            }
            // Killed: read the cell from the file name.
            let key = if let Some(rest) = name.strip_prefix("fresh-m") {
                let p: Vec<u64> = rest
                    .split(['-', '.'])
                    .take(3)
                    .map(|x| x.parse().expect("fresh-m<M>-<n>-<l>"))
                    .collect();
                (p[0], p[1], p[2])
            } else if let Some(rest) = name.strip_prefix("committed-cell-") {
                let p: Vec<u64> = rest
                    .split('-')
                    .take(2)
                    .map(|x| x.parse().expect("committed-cell-<n>-<l>"))
                    .collect();
                (3, p[0], p[1])
            } else {
                panic!("{path}: empty and unnamed");
            };
            let c = cells.entry(key).or_default();
            c.unmeasured += 1;
            c.degrees.push("killed".into());
            continue;
        };
        let v: Value = serde_json::from_str(line).expect("json");
        let scan = v["scan"].as_array().expect("scan");
        if let Some(bad) = v["refuted_satisfiable"].as_bool() {
            sat_draws += 1;
            sat_bad += bad as usize;
            continue;
        }
        let key = (
            v["m"].as_u64().unwrap(),
            v["n"].as_u64().unwrap(),
            v["ell"].as_u64().unwrap(),
        );
        let c = cells.entry(key).or_default();
        c.n_vars = v["n_vars"].as_u64().unwrap();
        let at4 = scan.iter().find(|s| s["degree"] == 4);
        c.xl_at_4 += at4.is_some_and(|s| s["xl_refuted"] == true) as usize;
        c.secs
            .push(scan.iter().map(|s| s["secs"].as_f64().unwrap()).sum());
        match v["closure_degree"].as_u64() {
            Some(d) if d <= 4 => {
                c.flat += 1;
                c.degrees.push(d.to_string());
                let s = scan.iter().find(|s| s["degree"] == d).unwrap();
                let cols = s["columns"].as_f64().unwrap();
                let ops = s["closure_word_ops"].as_f64().unwrap().max(1.0);
                c.ratios.push(ops.log2() / cols.log2());
            }
            _ if at4.is_some_and(|s| s["saturated"] == true) => {
                c.above += 1;
                c.degrees.push(">4".into());
            }
            _ => {
                c.unmeasured += 1;
                c.degrees.push(format!("to {}", v["measured_to"]));
            }
        }
    }
    println!("| m | (n, ℓ) | N | closure degree per draw | cell | XL refutes at 4 | log₂ ops ÷ log₂ cols | secs per draw |");
    println!("|--:|---|--:|---|---|--:|---|---|");
    let (mut flat_max, mut rise_min) = (None::<u64>, None::<u64>);
    let mut q5 = 0;
    let mut all_ratios = Vec::new();
    for ((m, n, ell), c) in &cells {
        let n_vars = if c.n_vars > 0 {
            c.n_vars
        } else {
            (m - 2) * n + m * ell
        };
        let verdict = if c.above > 0 {
            rise_min = Some(rise_min.map_or(n_vars, |r| r.min(n_vars)));
            "rises"
        } else if c.flat >= 3 {
            flat_max = Some(flat_max.map_or(n_vars, |f| f.max(n_vars)));
            "flat"
        } else {
            "unmeasured"
        };
        q5 += c.xl_at_4;
        all_ratios.extend(c.ratios.iter().copied());
        let r: Vec<String> = c.ratios.iter().map(|x| format!("{x:.2}")).collect();
        let t: Vec<String> = c.secs.iter().map(|x| format!("{x:.0}")).collect();
        println!(
            "| {m} | ({n}, {ell}) | {n_vars} | {} | {verdict} | {} | {} | {} |",
            c.degrees.join(" "),
            c.xl_at_4,
            r.join(" "),
            t.join(" ")
        );
    }
    println!();
    match (rise_min, flat_max) {
        (Some(y), _) => println!("Study: rises at N = {y}"),
        (None, Some(x)) => println!("Study: flat to N = {x}"),
        (None, None) => println!("Study: nothing measured"),
    }
    println!("Q5: draws where XL refutes at degree 4: {q5}");
    let (lo, hi) = all_ratios
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
    println!(
        "Q6: log₂ ops ÷ log₂ cols over {} refuting draws: {lo:.2}–{hi:.2} (registered band 2.0–2.4)",
        all_ratios.len()
    );
    println!(
        "Soundness control: {sat_draws} satisfiable draws, {sat_bad} refuted, {sat_unmeasured} unmeasured"
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("score") => score(&args[1..]),
        Some("compare") if args.len() == 3 => compare(&args[1], &args[2]),
        Some("semireg") if args.len() >= 3 => semireg(args[1].parse().expect("m"), &args[2..]),
        Some("shape") if args.len() >= 3 => shape(args[1].parse().expect("m"), &args[2..]),
        Some("closure") => closure(&args[1..]),
        Some("scale") => scale(&args[1..]),
        _ => {
            eprintln!(
                "usage: dreg_score score FILES… | compare REF RERUN | semireg M n:ℓ… | shape M n:ℓ:D… | closure FILES… | scale FILES…"
            );
            std::process::exit(2);
        }
    }
}
