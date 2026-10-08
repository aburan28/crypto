//! **Does cofactor balance of the factor base predict relation yield?**
//!
//! For every row of a committed cost sweep, compute the exact
//! cofactor-cancellation count `C(E)` of the member's factor base
//! ([`cofactor_cancellation_count`]) and compare the zero-parameter
//! prediction `relations ≈ probes · 2C/r` with the relations measured.
//! The design and its pre-registered criteria are in
//! `research/notes/koblitz-isogeny/yield-spread-design-20261003.md`.
//!
//! Run: `cargo run --release --example koblitz_yield_predictor <n> <l> [rows.jsonl ...]`
//! With no `.jsonl` arguments the rows come from
//! `experiments/koblitz_isogeny_cost_sweep.json` (the sweep with that `n`, `l`).

use crypto_lib::cryptanalysis::koblitz_isogeny_cost::{
    cofactor_cancellation_by_part, cofactor_cancellation_count, field_for,
};
use rayon::prelude::*;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let n: u32 = args.first().and_then(|v| v.parse().ok()).expect("n");
    let l: u32 = args.get(1).and_then(|v| v.parse().ok()).expect("l");
    let files: Vec<&String> = args.iter().skip(2).collect();

    let mut rows: Vec<serde_json::Value> = Vec::new();
    if files.is_empty() {
        let text = std::fs::read_to_string("experiments/koblitz_isogeny_cost_sweep.json")
            .expect("sweep snapshot");
        let v: serde_json::Value = serde_json::from_str(&text).expect("json");
        for s in v["sweeps"].as_array().expect("sweeps") {
            if s["n"].as_u64() == Some(n as u64) && s["l"].as_u64() == Some(l as u64) {
                // The family lives on the sweep, not on its rows.
                for r in s["rows"].as_array().expect("rows") {
                    let mut r = r.clone();
                    r["a2"] = s["a2"].clone();
                    rows.push(r);
                }
            }
        }
    } else {
        for f in files {
            for line in std::fs::read_to_string(f).expect("rows file").lines() {
                let v: serde_json::Value = serde_json::from_str(line).expect("json line");
                if v.get("attempted").is_none() {
                    rows.push(v);
                }
            }
        }
    }
    assert!(!rows.is_empty(), "no rows for n={n}, l={l}");
    let irr = field_for(n).expect("field");

    let out: Vec<serde_json::Value> = rows
        .par_iter()
        .map(|r| {
            let a6 = r["a6"].as_u64().unwrap();
            let a2 = r["a2"].as_u64().expect("a2") as u8;
            let order = r["order"].as_u64().unwrap();
            let rr = r["r"].as_u64().unwrap();
            let (pts, c) =
                cofactor_cancellation_count(n, &irr, a2, a6, l, Some(order)).expect("member");
            let parts = cofactor_cancellation_by_part(n, &irr, a2, a6, l, Some(order)).expect("member");
            let probes = r["trials"].as_u64().unwrap() as f64;
            let expected = probes * 2.0 * c as f64 / rr as f64;
            serde_json::json!({
                "a6": a6, "factor_base_points": pts, "unknowns": r["unknowns"],
                "C": c, "expected": expected,
                "parts": parts.iter().map(|(q, cq)| serde_json::json!({"prime_power": q, "C": cq})).collect::<Vec<_>>(),
                "relations": r["relations"], "probes": probes,
            })
        })
        .collect();

    let obs: Vec<f64> = out
        .iter()
        .map(|o| o["relations"].as_f64().unwrap())
        .collect();
    let exp: Vec<f64> = out
        .iter()
        .map(|o| o["expected"].as_f64().unwrap())
        .collect();
    let size: Vec<f64> = out
        .iter()
        .map(|o| o["unknowns"].as_f64().unwrap().powi(2))
        .collect();
    let k = obs.len() as f64;
    let calibration = obs.iter().sum::<f64>() / exp.iter().sum::<f64>();
    // Pearson dispersion of obs about the calibrated expectation, against
    // Poisson: ~1 when the predictor accounts for everything but counting noise.
    let dispersion = obs
        .iter()
        .zip(&exp)
        .map(|(o, e)| (o - calibration * e).powi(2) / (calibration * e).max(1e-12))
        .sum::<f64>()
        / (k - 1.0);
    let corr = |a: &[f64], b: &[f64]| {
        let (ma, mb) = (a.iter().sum::<f64>() / k, b.iter().sum::<f64>() / k);
        let cov: f64 = a.iter().zip(b).map(|(x, y)| (x - ma) * (y - mb)).sum();
        let va: f64 = a.iter().map(|x| (x - ma).powi(2)).sum();
        let vb: f64 = b.iter().map(|y| (y - mb).powi(2)).sum();
        cov / (va * vb).sqrt()
    };
    let (c_pred, c_size) = (corr(&obs, &exp), corr(&obs, &size));
    println!("n = {n}, l = {l}: {} members", obs.len());
    println!("   calibration Σobs/Σexpected = {calibration:.3}");
    println!("   dispersion of obs about calibrated expectation = {dispersion:.2}   (Poisson ≈ 1)");
    println!("   corr(obs, expected) = {c_pred:.3}    corr(obs, |F|²) = {c_size:.3}");
    let doc = serde_json::json!({
        "n": n, "l": l, "members": out.len(), "calibration": calibration,
        "dispersion": dispersion, "corr_expected": c_pred, "corr_size": c_size, "rows": out,
    });
    let path = format!("experiments/koblitz_yield_predictor_n{n}_l{l}.json");
    std::fs::write(&path, serde_json::to_string_pretty(&doc).unwrap() + "\n").unwrap();
    println!("   wrote {path}");
}
