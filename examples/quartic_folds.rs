//! **E18 of `RESEARCH_GLV_INVARIANT_FACTOR_BASES.md` (§9.1): folds of the
//! base on the `k = 4` linear algebra, against the matched rho.**
//!
//! One meet-in-the-middle relation stream per curve feeds the control
//! (`⟨−1⟩`, §11.19's matrix) and the folded arm (`ψ` on `j = 0`, `τ_T` on
//! a curve with rational 2-torsion); each is filtered, solved by
//! Wiedemann and verified; rho runs unfolded, negation-folded and (on
//! `j = 0`) `ψ`-folded on one code path.
//!
//! ```text
//! cargo run --release --example quartic_folds -- --family psi --sizes 277,541,769,1033 --seeds 2 --rho-runs 16 --json experiments/26_quartic_folds_psi.json
//! cargo run --release --example quartic_folds -- --family tau --sizes 269,521,769,1033 --seeds 2 --rho-runs 16 --json experiments/26_quartic_folds_tau.json
//! # --first-seed k runs seeds k, k+1, …: one curve per invocation survives a restart of the host
//! cargo run --release --example quartic_folds -- --family tau --sizes 1033 --first-seed 1 --seeds 1 --rho-runs 16 --json experiments/26_quartic_folds_tau_1033s1.json
//! ```
//!
//! Tables: `cargo run --release --example glv_invariant_experiment_tables
//! -- experiments/26_quartic_folds_psi.json experiments/26_quartic_folds_tau.json`.

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::quartic_folds::{run_k4_folds, Family4, RhoFold};
use serde_json::Value;

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut family = Family4::Psi;
    let mut sizes: Vec<u64> = vec![277, 541, 769, 1033];
    let mut seeds = 2u64;
    let mut first_seed = 1u64;
    let mut rho_runs = 16usize;
    let mut json: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--family" => {
                i += 1;
                family = match args[i].as_str() {
                    "psi" => Family4::Psi,
                    "tau" => Family4::Translation,
                    other => panic!("unknown family {other}; psi or tau"),
                };
            }
            "--sizes" => {
                i += 1;
                sizes = args[i]
                    .split(',')
                    .map(|v| v.parse().expect("--sizes"))
                    .collect();
            }
            "--first-seed" => {
                i += 1;
                first_seed = args[i].parse().expect("--first-seed");
            }
            "--seeds" => {
                i += 1;
                seeds = args[i].parse().expect("--seeds");
            }
            "--rho-runs" => {
                i += 1;
                rho_runs = args[i].parse().expect("--rho-runs");
            }
            "--json" => {
                i += 1;
                json = Some(args[i].clone());
            }
            "--probe" => {
                // Which of the sizes have an instance of the family (seed 1).
                for &p in &sizes {
                    let ok = match family {
                        Family4::Psi => {
                            crypto_lib::cryptanalysis::quartic_folds::generate_psi_instance(p, 1)
                                .map(|i| (i.n, i.cofactor))
                        }
                        Family4::Translation => {
                            crypto_lib::cryptanalysis::quartic_folds::generate_translation_instance(
                                p, 1, 4,
                            )
                            .map(|i| (i.n, i.cofactor))
                        }
                    };
                    println!("probe {} p={p}: {ok:?}", family.name());
                }
                return;
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    let mut rows: Vec<Value> = Vec::new();
    for &p in &sizes {
        for seed in first_seed..first_seed + seeds {
            let r = match run_k4_folds(family, p, seed, rho_runs) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("e18 {} p={p} seed={seed}: {e}", family.name());
                    continue;
                }
            };
            let walk = |f: RhoFold| -> (f64, f64, bool) {
                let w: Vec<_> = r.rho.iter().filter(|w| w.fold == f).collect();
                let muls = w.iter().map(|w| w.fp_muls as f64).sum::<f64>() / w.len().max(1) as f64;
                let steps = w.iter().map(|w| w.steps as f64).sum::<f64>() / w.len().max(1) as f64;
                (muls, steps, w.iter().all(|w| w.correct))
            };
            let (m0, s0, ok0) = walk(RhoFold::None);
            let (m1, s1, ok1) = walk(RhoFold::Negation);
            let (m2, s2, ok2) = walk(RhoFold::Psi);
            let matched = if family == Family4::Psi { m2 } else { m1 };
            let ratio = |la: u64, m: f64| la as f64 * 16.0 / m;
            println!(
                "e18 {} p={} seed={} n=2^{:.1} h={} |F|={} cols {}/{} ({:.2}/col) | residuals {} rate {:.4} | control: rel {} core {} φ {:.3} LA {:.3e} | folded: rel {} core {} φ {:.3} LA {:.3e} | LA ratio {:.2} | rho steps none/neg/psi {:.0}/{:.0}/{:.0} | r control/folded vs matched {:.3}/{:.3}, vs unfolded {:.3}/{:.3} | correct {}/{} rho {}/{}/{} | {:.0} s",
                family.name(), r.p, r.seed, r.bits, r.cofactor, r.base_points, r.control.columns, r.folded.columns, r.points_per_column,
                r.residuals, r.decomposition_rate,
                r.control.relations, r.control.unknowns, r.control.phi, r.control.la_ops as f64,
                r.folded.relations, r.folded.unknowns, r.folded.phi, r.folded.la_ops as f64,
                r.control.la_ops as f64 / r.folded.la_ops.max(1) as f64,
                s0, s1, s2,
                ratio(r.control.la_ops, matched), ratio(r.folded.la_ops, matched),
                ratio(r.control.la_ops, m0), ratio(r.folded.la_ops, m0),
                r.control.correct, r.folded.correct, ok0, ok1, if family == Family4::Psi { ok2.to_string() } else { "—".into() },
                r.wall_ms / 1e3
            );
            let mut v = serde_json::to_value(&r).unwrap();
            v["experiment"] = Value::from("e18");
            rows.push(v);
            if let Some(path) = &json {
                fs::write(path, serde_json::to_string_pretty(&rows).unwrap()).expect("write json");
            }
        }
    }
}
