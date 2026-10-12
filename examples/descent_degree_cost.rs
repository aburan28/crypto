//! What the Weil-descent route costs at ECC2K-130 as a function of the
//! Macaulay degree: the native port of
//! `research/notes/ecc2k130/descent_degree_cost_20260930.py`, which stays in
//! the repository as the degree-7 study's record.
//!
//! ```sh
//! cargo run --release --example descent_degree_cost
//! ```
//!
//! Derived, not measured.  The cost model of
//! `RESEARCH_ECC2K130_DECOMPOSITION.md` §5.1 with the oracle priced per
//! target, from below, by the size of the degree-`D` Macaulay matrix of the
//! Weil-descended chained `S₃` system:
//!
//! ```text
//!     relations      |F| = 2^l                      (l = dim V, an integer)
//!     targets        m! * 2^{n-(m-1)l}, at least one per relation
//!     oracle         C(N, <=D)^w per target, N = (m-2) n + m l
//!     linear algebra m * 2^{2l}
//! ```
//!
//! `w = 1` prices one touch per column, which no elimination beats; `w = 2`
//! is the usual estimate at `ω = 2`.  Units: column touches against rho's
//! group operations, not converted.  Every row is minimised over integer
//! `l` in 4..70, as the Python did; the rows the Python printed come out
//! identical, and the `m = 4` semi-regular rows are new (the `m = 4` study,
//! `research/dreg_m4_chain_20261006/`).

use crypto_lib::cryptanalysis::koblitz_bench::semi_regular_dreg;
use std::collections::HashMap;

const N_FIELD: usize = 131;
const RHO: f64 = 60.81; // RESEARCH_ECC2K130_DECOMPOSITION.md §5.2, <-1> x <pi>

fn log2_add(a: f64, b: f64) -> f64 {
    a.max(b) + (1.0 + (-(a - b).abs()).exp2()).log2()
}

/// log2 C(N, <= d).
fn log2_monomials(n_vars: usize, d: u32) -> f64 {
    let (mut total, mut binom) = (0f64, 1f64);
    for i in 0..=d as usize {
        if i > n_vars {
            break;
        }
        total += binom;
        binom = binom * (n_vars - i) as f64 / (i + 1) as f64;
    }
    total.log2()
}

fn factorial(m: usize) -> f64 {
    (1..=m).map(|k| k as f64).product()
}

struct Cost {
    total: f64,
    targets: f64,
    oracle: f64,
    algebra: f64,
    unknowns: usize,
}

fn cost(m: usize, l: usize, degree: u32, w: u32) -> Cost {
    let n = N_FIELD;
    let unknowns = (m - 2) * n + m * l;
    let per_relation = (factorial(m).log2() + n as f64 - (m * l) as f64).max(0.0);
    let targets = l as f64 + per_relation;
    let oracle = if degree == 0 {
        0.0
    } else {
        w as f64 * log2_monomials(unknowns, degree)
    };
    let algebra = (m as f64).log2() + (2 * l) as f64;
    Cost {
        total: log2_add(targets + oracle, algebra),
        targets,
        oracle,
        algebra,
        unknowns,
    }
}

/// The equation degrees of the `m`-summand chain at `n = 131`.
fn chain_degrees(m: usize) -> Vec<u32> {
    std::iter::repeat_n(3, (m - 2) * N_FIELD)
        .chain(std::iter::repeat_n(2, N_FIELD))
        .collect()
}

fn main() {
    let n = N_FIELD;
    let mut dreg_cache: HashMap<(usize, usize), u32> = HashMap::new();
    let mut dreg = |m: usize, l: usize| -> u32 {
        *dreg_cache.entry((m, l)).or_insert_with(|| {
            semi_regular_dreg((m - 2) * n + m * l, &chain_degrees(m)).expect("reference exists")
        })
    };
    println!(
        "| m | degree D | w | best l | unknowns N | D | targets | oracle per target | linear algebra | total | vs rho |"
    );
    println!("|--:|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    for m in [3usize, 4] {
        let mut scenarios: Vec<(
            &str,
            Box<dyn Fn(usize, &mut dyn FnMut(usize, usize) -> u32) -> u32>,
        )> = vec![
            ("free oracle (floor)", Box::new(|_, _| 0)),
            ("2 (any Macaulay matrix)", Box::new(|_, _| 2)),
            ("6, constant", Box::new(|_, _| 6)),
            ("7, constant", Box::new(|_, _| 7)),
        ];
        scenarios.push((
            "semi-regular D_reg(N)",
            Box::new(move |l, dreg: &mut dyn FnMut(usize, usize) -> u32| dreg(m, l)),
        ));
        // The measured degree is flat in n at fixed l and rises with l
        // (RESEARCH_DREG_MEASUREMENT.md Results 5-10).  The two small-l
        // fits, extrapolated: m = 3 reads ceil(l/2) + 4 at l = 2..5 except
        // (4,3) and (7,4) at surplus -5 and half of (8,4)'s draws; m = 4
        // reads l + 4 at l = 1, 2 and >= 7 at l = 3.
        if m == 3 {
            scenarios.push((
                "ceil(l/2) + 4, the m = 3 fit",
                Box::new(|l, _| (l as u32).div_ceil(2) + 4),
            ));
        } else {
            scenarios.push(("l + 4, the m = 4 fit", Box::new(|l, _| l as u32 + 4)));
        }
        for (name, degree_of) in &scenarios {
            let ws: &[u32] = if name.starts_with("free") {
                &[1]
            } else {
                &[1, 2]
            };
            for &w in ws {
                let mut best: Option<(f64, usize)> = None;
                for l in 4..70 {
                    let t = cost(m, l, degree_of(l, &mut dreg), w).total;
                    if best.is_none_or(|(b, _)| t < b) {
                        best = Some((t, l));
                    }
                }
                let (_, l) = best.expect("a row");
                let d = degree_of(l, &mut dreg);
                let c = cost(m, l, d, w);
                println!(
                    "| {m} | {name} | {w} | {l} | {} | {} | 2^{:.2} | 2^{:.2} | 2^{:.2} | **2^{:.2}** | 2^{:+.2} |",
                    c.unknowns,
                    if d == 0 { "—".to_string() } else { d.to_string() },
                    c.targets,
                    c.oracle,
                    c.algebra,
                    c.total,
                    c.total - RHO
                );
            }
        }
    }
    // The table searches l in 4..70, as the 2026-09-30 model did, and the
    // fit rows' minima sit at l = 4.  Disclose what lies below that range:
    // the w = 1 cost of the fits at l = 1..4, with the target count, which
    // passes the 2^129 subgroup there.
    println!();
    for (m, fit) in [(3usize, "ceil(l/2) + 4"), (4, "l + 4")] {
        let rows: Vec<String> = (1..=4usize)
            .map(|l| {
                let d = if m == 3 {
                    (l as u32).div_ceil(2) + 4
                } else {
                    l as u32 + 4
                };
                let c = cost(m, l, d, 1);
                format!(
                    "l={l} (D={d}): 2^{:.2} (targets 2^{:.2})",
                    c.total, c.targets
                )
            })
            .collect();
        println!(
            "below the search range, m = {m}, {fit}, w = 1: {}",
            rows.join("; ")
        );
    }
    println!();
    for m in [3usize, 4] {
        let ls: [usize; 5] = [20, 26, 27, 33, 43];
        let shown: Vec<String> = ls
            .iter()
            .map(|&l| format!("l={l}: {}", dreg(m, l)))
            .collect();
        println!(
            "semi-regular D_reg of the m = {m} shape at n = 131: {}",
            shown.join(", ")
        );
    }
}
