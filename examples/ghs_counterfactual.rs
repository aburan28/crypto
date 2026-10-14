//! Composite-degree GHS counterfactual of
//! `research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md` §7 under the
//! Maurer–Menezes–Teske magic number, as JSON.
//!
//! ```text
//! cargo run --release --example ghs_counterfactual -- --out results.json
//! ```
//!
//! The JSON holds, for `N = 130` and `N = 131`, the table under the MMT rule and
//! under the pre-2026-10-13 rule (a regression check of the port against the
//! note), and the model cost against genus for `N = 130` (`l = 1`) with the rho
//! reference.  Every cost is a model figure for a hypothetical curve.

use crypto_lib::cryptanalysis::curve_construction::counterfactual::{
    counterfactual, log2_rho, Row, Rule,
};
use crypto_lib::cryptanalysis::curve_construction::window::{price, rational_model_points};
use serde::Serialize;

#[derive(Serialize)]
struct SweepPoint {
    genus: usize,
    smoothness_bound_b: usize,
    log2_cost: f64,
    log2_relations: f64,
    log2_linear_algebra: f64,
    log2_ratio_to_rho: f64,
    beats_rho: bool,
}

#[derive(Serialize)]
struct Output {
    schema: &'static str,
    command: String,
    model: &'static str,
    admissibility_rule: &'static str,
    tables: Vec<Table>,
    sweep_n130_l1: Vec<SweepPoint>,
    log2_rho_n130: f64,
}

#[derive(Serialize)]
struct Table {
    field_degree: u32,
    rule: Rule,
    rows: Vec<Row>,
}

fn main() {
    let mut out_path = None;
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--out" => out_path = args.next(),
            other => panic!("unknown argument {other}"),
        }
    }
    let mut tables = Vec::new();
    for n in [130u32, 131] {
        for rule in [Rule::Mmt, Rule::Legacy] {
            tables.push(Table {
                field_degree: n,
                rule,
                rows: counterfactual(n, rule),
            });
        }
    }
    let rho = log2_rho(130);
    let genera: Vec<usize> = (96..=320).step_by(8).chain([255, 256]).collect();
    let mut genera = genera;
    genera.sort_unstable();
    genera.dedup();
    let sweep = genera
        .into_iter()
        .map(|g| {
            let c = price(&rational_model_points(g), g, 200, "random-polynomial", rho);
            SweepPoint {
                genus: g,
                smoothness_bound_b: c.smoothness_bound_b,
                log2_cost: c.log2_total,
                log2_relations: c.log2_relations,
                log2_linear_algebra: c.log2_linear_algebra,
                log2_ratio_to_rho: c.log2_ratio_to_rho,
                beats_rho: c.beats_rho,
            }
        })
        .collect();
    let out = Output {
        schema: "ghs_counterfactual/v1",
        command: "cargo run --release --example ghs_counterfactual -- --out results.json".into(),
        model: "curve_construction::window::price, random-polynomial point counts #C(F_2^k) = 2^k + 1, smoothness bound b <= 200; fixed-genus g!q^(2-2/g) for l >= 2",
        admissibility_rule: "l * genus >= N - 2 (the note's rule, unchanged)",
        tables,
        sweep_n130_l1: sweep,
        log2_rho_n130: rho,
    };
    let json = serde_json::to_string_pretty(&out).expect("serialise");
    match out_path {
        Some(p) => std::fs::write(&p, json + "\n").expect("write"),
        None => println!("{json}"),
    }
}
