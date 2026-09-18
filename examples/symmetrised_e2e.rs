//! X4: unchained symmetrised oracle through end-to-end collection.
//!
//! Frozen protocol in `RESEARCH_ECC2K130_ROUTE_TARGETS.md` (2026-09-18).
//! Do not change `m`, the divisor, `Q_enum`, or the seed after seeing a cell.
//!
//! ```bash
//! cargo run --release --example symmetrised_e2e -- --a 1 --n 7,9
//! cargo run --release --example symmetrised_e2e -- --a 1 --n 7,9,15 --control
//! ```
//!
//! Unit: `Λ = trials · C(|F_u|, m−1) / 2^n`. Ratio to the floor is
//! `Λ · n / m`. Algebraic reductions are a stage diagnostic until T3.
//! Wall-clock is a footnote.

use crypto_lib::cryptanalysis::koblitz_groebner::FieldStructure;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    enumerate_decompose, koblitz_index_calculus_dlp_with_factor_base, DecompositionStrategy,
    KoblitzCurve, KoblitzIcOptions,
};
use crypto_lib::cryptanalysis::koblitz_symmetrised::{
    prepare_symmetrised_attack, symmetrised_groebner_decompose,
};
use num_bigint::BigUint;
use serde_json::json;
use std::env;
use std::sync::Arc;
use std::time::Instant;

const M: usize = 3;
const SEED: u64 = 0x5EED_0004;
const SAMPLE: u32 = 12;

fn binom(n: u64, k: u64) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut acc = 1.0f64;
    for i in 0..k {
        acc *= (n - i) as f64;
        acc /= (i + 1) as f64;
    }
    acc
}

fn parse_ns(s: &str) -> Vec<u32> {
    s.split(',')
        .filter_map(|t| t.trim().parse().ok())
        .collect()
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut a = 1u8;
    let mut ns = vec![7u32, 9];
    let mut control = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--a" => {
                a = args.get(i + 1).and_then(|s| s.parse().ok()).unwrap_or(a);
                i += 2;
            }
            "--n" => {
                if let Some(s) = args.get(i + 1) {
                    ns = parse_ns(s);
                }
                i += 2;
            }
            "--control" => {
                control = true;
                i += 1;
            }
            other => {
                eprintln!("unknown argument: {other}");
                std::process::exit(2);
            }
        }
    }

    println!("X4 protocol m={M} seed={SEED:#x} a={a} n={ns:?} control={control}");
    println!(
        "{:>2} {:>3} {:>4} {:>5} {:>6} {:>8} {:>8} {:>8} {:>8} {:>10} {:>6} {}",
        "n", "ell", "|Fu|", "adm", "agree", "trials", "C(F,2)", "Lambda", "ratio", "reductions",
        "ok", "note"
    );

    let mut rows = Vec::new();
    for n in ns {
        let Some(kc) = KoblitzCurve::new(a, n) else {
            println!("{n:>2} skip (no curve)");
            continue;
        };
        let setup = Instant::now();
        let Some((fu, view)) = prepare_symmetrised_attack(&kc, M) else {
            println!("{n:>2} skip (|F_u|<3 or ell=1 or no view)");
            rows.push(json!({
                "a": a,
                "n": n,
                "skip": true,
                "reason": "prepare_symmetrised_attack returned None",
            }));
            continue;
        };
        let setup_ns = setup.elapsed().as_nanos();
        let admissible = view.m_can_decompose(&kc, M);
        let fu_len = fu.points.len();
        let c_enum = binom(fu_len as u64, (M - 1) as u64);
        let index_of = view.index_map();
        let g = kc.generator().clone();
        let field = FieldStructure::new(kc.n, &kc.curve.irreducible);

        let mut agreed = 0usize;
        let mut completed = 0usize;
        let mut disagreed = 0usize;
        let mut sample_found = 0usize;
        if admissible {
            for k in 1..=SAMPLE {
                let target = kc.mul(&g, &BigUint::from(k));
                let enumerated = enumerate_decompose(&kc, &view, &index_of, &target, M);
                match symmetrised_groebner_decompose(
                    &kc, &fu, &field, &target, M, Default::default(), 50_000,
                ) {
                    Some(o) if o.complete => {
                        completed += 1;
                        if o.relation.is_some() {
                            sample_found += 1;
                        }
                        if o.relation.is_some() == enumerated.is_some() {
                            agreed += 1;
                        } else {
                            disagreed += 1;
                        }
                    }
                    _ => {}
                }
            }
        }

        let d = BigUint::from(29u32);
        let q = kc.mul(&g, &d);
        let fu_ell = fu.ell;
        let fu = Arc::new(fu);
        let mut opts = KoblitzIcOptions {
            m: M,
            strategy: DecompositionStrategy::Symmetrised,
            node_budget: 50_000,
            extra_relations: 4,
            max_trials: 20_000,
            allow_direct_relation: false,
            stop_on_verified_rank: true,
            collapse_negation: true,
            seed: SEED,
            symmetrised_fb: Some(fu),
            ..KoblitzIcOptions::default()
        };
        let run = Instant::now();
        let report = if admissible {
            koblitz_index_calculus_dlp_with_factor_base(&kc, &q, &view, &opts)
        } else {
            None
        };
        let wall_ns = run.elapsed().as_nanos();
        let verified = report
            .as_ref()
            .and_then(|r| r.log.as_ref())
            .is_some_and(|log| kc.mul(&g, log) == q && *log == d);
        let trials = report.as_ref().map(|r| r.trials).unwrap_or(0);
        let lambda = trials as f64 * c_enum / 2f64.powi(n as i32);
        let ratio = lambda * (n as f64) / (M as f64);
        let reductions = report.as_ref().map(|r| r.reductions).unwrap_or(0);
        let la_ns = report.as_ref().map(|r| r.linear_algebra_ns).unwrap_or(0);
        let ok = if verified { "yes" } else { "no" };
        let agree = if disagreed == 0 && completed > 0 {
            format!("{agreed}/{completed}")
        } else {
            format!("{agreed}/{completed} d={disagreed}")
        };
        let note = if !admissible {
            "m not cofactor-admissible"
        } else if disagreed > 0 {
            "Enumerate disagreement"
        } else if !verified {
            "log missing or wrong"
        } else {
            "verified"
        };
        println!(
            "{n:>2} {ell:>3} {fu:>4} {adm:>5} {agree:>6} {trials:>8} {c:>8.1} {lambda:>8.4} {ratio:>8.4} {reductions:>10} {ok:>6} {note}",
            ell = fu_ell,
            fu = fu_len,
            adm = if admissible { "yes" } else { "no" },
            c = c_enum,
        );

        let mut row = json!({
            "a": a,
            "n": n,
            "m": M,
            "ell": fu_ell,
            "view_ell": view.ell,
            "sample_found": sample_found,
            "fu": fu_len,
            "orbits": view.unknowns(),
            "admissible": admissible,
            "setup_ns": setup_ns,
            "wall_ns": wall_ns,
            "agree_completed": completed,
            "agree_ok": agreed,
            "agree_disagree": disagreed,
            "trials": trials,
            "c_enum": c_enum,
            "lambda": lambda,
            "ratio": ratio,
            "reductions": reductions,
            "infeasible_branches": report.as_ref().map(|r| r.infeasible_branches),
            "linear_algebra_ns": la_ns,
            "relation_collection_ns": report.as_ref().map(|r| r.relation_collection_ns),
            "rank_checks": report.as_ref().map(|r| r.rank_checks),
            "relations": report.as_ref().map(|r| r.relations),
            "verified": verified,
            "direct_relation": report.as_ref().map(|r| r.direct_relation),
            "note": note,
            "unit": "lambda = trials * C(|F_u|, m-1) / 2^n",
            "floor": "m/n",
        });

        if control && admissible {
            opts.strategy = DecompositionStrategy::Enumerate;
            opts.symmetrised_fb = None;
            let enum_run = Instant::now();
            let enum_report = koblitz_index_calculus_dlp_with_factor_base(&kc, &q, &view, &opts);
            let enum_verified = enum_report
                .as_ref()
                .and_then(|r| r.log.as_ref())
                .is_some_and(|log| kc.mul(&g, log) == q && *log == d);
            let enum_trials = enum_report.as_ref().map(|r| r.trials).unwrap_or(0);
            let enum_lambda = enum_trials as f64 * c_enum / 2f64.powi(n as i32);
            row["enumerate"] = json!({
                "trials": enum_trials,
                "lambda": enum_lambda,
                "ratio": enum_lambda * (n as f64) / (M as f64),
                "verified": enum_verified,
                "wall_ns": enum_run.elapsed().as_nanos(),
                "linear_algebra_ns": enum_report.as_ref().map(|r| r.linear_algebra_ns),
            });
        }
        rows.push(row);
    }
    println!("{}", json!({ "protocol": "X4-2026-09-18", "rows": rows }));
}
