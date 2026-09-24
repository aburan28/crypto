//! **Admissibility screen for chained Gröbner-stage cells.**  Decides no
//! target: for every `(a, n, factor index, m)` of a grid it reports whether
//! the curve and the factor base exist, whether `m` summands can close the
//! group identity at all (`FrobeniusFactorBase::m_can_decompose`, the check
//! `groebner_stage_bench` skips a cell on), and how many Boolean unknowns the
//! chained system has (`m·ℓ + (m − 2)·n`, at most 64).
//!
//! Companion to §5 of `research/notes/ecc2k130/RESEARCH_CHAIN_SPLIT_ORDER.md`,
//! which chose its supplementary holdout from this screen's output before any
//! arm ran on it.
//!
//! ```text
//! cargo run --release --example chain_ladder_screen > research/chain_split_order_20260924/screen.json
//! ```

use crypto_lib::cryptanalysis::koblitz_groebner::MAX_VARS;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_frobenius_factor_base, KoblitzCurve,
};

fn main() {
    let mut cells = Vec::new();
    for n in [9u32, 11, 13, 15, 17, 19, 23] {
        for a in [0u8, 1] {
            let Some(kc) = KoblitzCurve::new(a, n) else {
                cells.push(serde_json::json!({"a": a, "n": n, "curve": false}));
                continue;
            };
            for factor_index in 0..4usize {
                let Some(fb) = build_frobenius_factor_base(&kc, factor_index) else {
                    continue;
                };
                for m in [3usize, 4] {
                    let n_vars = m * fb.ell as usize + (m - 2) * n as usize;
                    cells.push(serde_json::json!({
                        "a": a,
                        "n": n,
                        "curve": true,
                        "factor_index": factor_index,
                        "ell": fb.ell,
                        "factor_base_points": fb.points.len(),
                        "m": m,
                        "n_vars": n_vars,
                        "fits": n_vars <= MAX_VARS,
                        "admissible": fb.m_can_decompose(&kc, m),
                    }));
                }
            }
        }
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::json!({
            "screen": "chained-cell admissibility; decides no target",
            "cells": cells,
        }))
        .unwrap()
    );
}
