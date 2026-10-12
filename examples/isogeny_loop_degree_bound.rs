//! Class-group-free lower bound on the cost of any isogeny-loop endomorphism.
//!
//! A loop of horizontal isogenies of prime degrees l_1, ..., l_m is an
//! endomorphism `alpha` of degree `N(alpha) = l_1 ... l_m`, and costs
//! `sum c(l_i) >= KAPPA * log2 N(alpha)` where `KAPPA = min_l c(l)/log2 l`
//! is the cheapest cost per bit of degree in the cost model.  A non-scalar
//! `alpha = (A + B sqrt D)/2` in an order of discriminant `D` has `B != 0`,
//! so `N(alpha) = (A^2 + |D| B^2)/4 >= d_min(D) = floor((|D| + 1)/4)`, and
//! `d_min` only grows in sub-orders (`D = g^2 D_K`).  Hence every loop
//! endomorphism of every curve whose Frobenius has fundamental discriminant
//! `D_K` costs at least `KAPPA * log2 d_min(D_K)`, whatever the class group.
//! A 2-dimensional GLV split saves `c_dbl * ceil(log2 r / 2)`; when the bound
//! exceeds that, no loop can pay for itself under that cost model.
//!
//! Scope of the bound: separable isogenies of odd characteristic and the
//! separable steps in characteristic 2.  In characteristic 2 the Frobenius
//! `F` is an inseparable 2-isogeny costing two squarings, which is cheaper per
//! bit than any `KAPPA` here, so for binary curves the bound covers only
//! loops without `F` steps; loops that use `F` were covered by the earlier
//! exhaustive weighted search (`research/isogeny_loop_glv_20261010/`).
//!
//! Inputs are the frozen sweep records of `research/isogeny_loop_glv_20261010/`
//! (produced by the legacy Python sweeper; see that directory's README).  For
//! each curve this re-verifies, in native big-integer arithmetic, that
//! `t^2 - 4q = f^2 * D_K` with the recorded `D_K` and conductor `f`, then
//! prints the bound under three per-step cost models:
//!
//! * `projective`: the paper's `7.5 l` (Section 2.2), `KAPPA = 22.5/log2 3`;
//! * `affine`: Kohel affine `2 l + 2`, inversions free, `KAPPA = 8/log2 3`;
//! * `xonly`: the cheapest published x-only evaluation per step (3-isogeny
//!   `4M + 2S`, counted `S = M`), `KAPPA = 6/log2 3`.  Optimistic: it prices a
//!   GLV split that x-only arithmetic does not directly support.
//!
//! The saving uses the paper's doubling cost `c_dbl = 8 M` in all three rows,
//! so only the per-step model varies.
//!
//! ```bash
//! cargo run --release --example isogeny_loop_degree_bound -- \
//!     research/isogeny_loop_glv_20261010/group_*.json [--json out.json]
//! ```

use num_bigint::BigInt;
use num_traits::{Signed, ToPrimitive, Zero};
use serde_json::Value;
use std::str::FromStr;

fn big(v: &Value) -> Option<BigInt> {
    match v {
        Value::String(s) => BigInt::from_str(s).ok(),
        Value::Number(n) => BigInt::from_str(&n.to_string()).ok(),
        _ => None,
    }
}

fn log2(x: &BigInt) -> f64 {
    let x = x.abs();
    let bits = x.bits();
    if bits <= 1000 {
        let shift = bits.saturating_sub(64);
        let top = (&x >> shift).to_u64().unwrap() as f64;
        top.log2() + shift as f64
    } else {
        unreachable!()
    }
}

struct Model {
    name: &'static str,
    kappa: f64,
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let json_out = args
        .iter()
        .position(|a| a == "--json")
        .and_then(|i| args.get(i + 1))
        .cloned();
    let files: Vec<&String> = args
        .iter()
        .enumerate()
        .filter(|(i, a)| a.ends_with(".json") && (*i == 0 || args[i - 1] != "--json"))
        .map(|(_, a)| a)
        .collect();
    let l3 = 3f64.log2();
    let models = [
        Model {
            name: "projective",
            kappa: 22.5 / l3,
        },
        Model {
            name: "affine",
            kappa: 8.0 / l3,
        },
        Model {
            name: "xonly",
            kappa: 6.0 / l3,
        },
    ];
    let c_dbl = 8.0;
    println!("# Class-group-free lower bound on loop-endomorphism cost\n");
    println!(
        "KAPPA (M per bit of degree): projective {:.3}, affine {:.3}, xonly {:.3}; saving = {c_dbl} * ceil(log2 r / 2) M.\n",
        models[0].kappa, models[1].kappa, models[2].kappa
    );
    println!("| curve | family | log2 r | log2 d_min(D_K) | Delta = f^2 D_K verified | saving M | bound projective M (ratio) | bound affine M (ratio) | bound xonly M (ratio) | cheapest loop found, M (frozen; projective) |");
    println!("|:--|:--|--:|--:|:--|--:|--:|--:|--:|--:|");
    let mut rows = vec![];
    for f in files {
        let doc: Value =
            serde_json::from_str(&std::fs::read_to_string(f).expect("read")).expect("json");
        for res in doc["results"].as_array().expect("results") {
            let curve = &res["curve"];
            let fr = &res["frobenius"];
            let name = curve["name"].as_str().unwrap_or("?").to_string();
            let family = curve["family"].as_str().unwrap_or("?").to_string();
            let (Some(q), Some(t), Some(r), Some(dk), Some(fc)) = (
                big(&curve["q"]),
                big(&curve["trace"]),
                big(&curve["r"]),
                big(&fr["DK"]),
                big(&fr["conductor_f"]),
            ) else {
                println!("| {name} | {family} | — | — | record incomplete | | | | | |");
                continue;
            };
            let delta = &t * &t - BigInt::from(4) * &q;
            let verified = delta == &fc * &fc * &dk;
            let complete = fr["unfactored_cofactor"].is_null();
            let ad = dk.abs();
            let dmin: BigInt = if (&ad % 4u32).is_zero() {
                &ad / 4u32
            } else {
                (&ad + 1u32) / 4u32
            };
            let ldmin = if dmin > BigInt::from(1) {
                log2(&dmin)
            } else {
                0.0
            };
            let lr = r.bits();
            let saving = c_dbl * lr.div_ceil(2) as f64;
            let measured = res["orders"]
                .as_array()
                .and_then(|o| o.first())
                .and_then(|o| o["cheapest_non_frobenius_loop"]["cost_M"].as_f64());
            let mut cells = vec![];
            let mut bounds = serde_json::Map::new();
            for m in &models {
                let b = m.kappa * ldmin;
                cells.push(format!("{:.0} ({:.2})", b, b / saving));
                bounds.insert(
                    m.name.into(),
                    serde_json::json!({"bound_M": b, "ratio": b / saving}),
                );
            }
            // the bound must not exceed what was actually found (odd characteristic only:
            // binary searches used Frobenius steps, which the bound excludes)
            if let (Some(mc), false) = (measured, family == "binary" || family == "koblitz") {
                assert!(
                    models[0].kappa * ldmin <= mc + 1e-6,
                    "{name}: bound above a found loop"
                );
            }
            println!(
                "| {name} | {family} | {lr} | {:.1} | {}{} | {:.0} | {} | {} | {} | {} |",
                ldmin,
                if verified { "yes" } else { "NO" },
                if complete {
                    ""
                } else {
                    " (D_K assumes the recorded squarefree cofactor)"
                },
                saving,
                cells[0],
                cells[1],
                cells[2],
                measured
                    .map(|m| format!("{m:.1}"))
                    .unwrap_or_else(|| "not computed".into()),
            );
            rows.push(serde_json::json!({
                "curve": name, "family": family, "log2_r": lr, "log2_dmin": ldmin,
                "delta_verified": verified, "factorization_complete": complete,
                "saving_M": saving, "bounds": bounds, "cheapest_found_M": measured, "source": f,
            }));
        }
    }
    println!(
        "\nRatio >= 1: no loop endomorphism of that curve can cost less than the doublings a 2-dimensional split saves, under that per-step model."
    );
    if let Some(path) = json_out {
        let doc = serde_json::json!({
            "schema": "isogeny_loop_degree_bound/v1",
            "kappa": {"projective": models[0].kappa, "affine": models[1].kappa, "xonly": models[2].kappa},
            "c_dbl": c_dbl, "rows": rows,
        });
        std::fs::write(&path, serde_json::to_string_pretty(&doc).unwrap()).expect("write");
        eprintln!("wrote {path}");
    }
}
