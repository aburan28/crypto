//! Score, compare and reference refutation-degree runs (`dreg_ladder` output).
//!
//! ```sh
//! # per-cell values and the registered reading, with the semi-regular reference
//! cargo run --release --example dreg_score -- score runs/*.jsonl
//! # draw-by-draw identity of a rerun against committed rows (all but `secs`)
//! cargo run --release --example dreg_score -- compare committed.jsonl rerun.jsonl
//! # the semi-regular degree of regularity of chained-system shapes
//! cargo run --release --example dreg_score -- semireg 4 5:2 6:2 7:2
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

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("score") => score(&args[1..]),
        Some("compare") if args.len() == 3 => compare(&args[1], &args[2]),
        Some("semireg") if args.len() >= 3 => semireg(args[1].parse().expect("m"), &args[2..]),
        _ => {
            eprintln!("usage: dreg_score score FILES… | compare REF RERUN | semireg M n:ℓ…");
            std::process::exit(2);
        }
    }
}
