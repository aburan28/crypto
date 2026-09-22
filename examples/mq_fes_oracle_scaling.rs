//! Stage diagnostic: per-call cost of the quadratic Semaev oracle against
//! ℓ = dim V, previous oracle against the packed linear split.
//!
//! The linear split's packed template is built (and timed) once per factor
//! base before its calls, as the pipeline pays it; the previous oracle has
//! no template and rebuilds its system every call.
//!
//! Only curves with `n ≥ 2ℓ` are measured: there the system is
//! overdetermined, almost every target is refuted, and a call pays its full
//! enumeration, so the exponent is not confounded by early exits.  Both
//! oracles see the same targets and must agree on which ones decompose.
//!
//!     cargo run --release --example mq_fes_oracle_scaling > out.jsonl
use crypto_lib::cryptanalysis::koblitz_groebner::FieldStructure;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_frobenius_factor_base, point_key, KoblitzCurve,
};
use crypto_lib::cryptanalysis::mq_fes::{
    mq_fes_decompose, mq_fes_decompose_reference, mq_fes_profile, mq_fes_profile_reset,
};
use crypto_lib::cryptanalysis::mq_fes_semaev::packed_template;
use num_bigint::BigUint;
use std::collections::HashMap;
use std::time::Instant;

const TARGETS: u64 = 12;
const MAX_ELL: usize = 14;

fn main() {
    let mut seen = std::collections::BTreeSet::new();
    for n in (7u32..=63).step_by(2) {
        for a in [0u8, 1] {
            let Some(kc) = KoblitzCurve::new(a, n) else { continue };
            for index in 0..4 {
                let Some(fb) = build_frobenius_factor_base(&kc, index) else { break };
                let ell = fb.subspace_basis.len();
                if ell < 3 || ell > MAX_ELL || (n as usize) < 2 * ell || !seen.insert((ell, n)) {
                    continue;
                }
                let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
                let index_of: HashMap<_, _> = fb
                    .points
                    .iter()
                    .enumerate()
                    .map(|(i, p)| (point_key(p), i))
                    .collect();
                let targets: Vec<_> = (0..TARGETS)
                    .map(|t| kc.mul(kc.generator(), &BigUint::from(1_000_003u64 * (t + 1) + 17)))
                    .collect();
                // One-off, target-independent: paid once per factor base.
                let template_start = Instant::now();
                packed_template(&fb.subspace_basis, &kc.curve.b, &st).expect("template");
                let template_us = template_start.elapsed().as_nanos() as f64 / 1e3;
                let mut arms = Vec::new();
                let mut answers = Vec::new();
                for arm in ["baseline", "mqfes-linear"] {
                    mq_fes_profile_reset();
                    let start = Instant::now();
                    let found: Vec<bool> = targets
                        .iter()
                        .map(|t| {
                            if arm == "baseline" {
                                mq_fes_decompose_reference(&kc, &fb, &index_of, &st, t, 2).0.is_some()
                            } else {
                                mq_fes_decompose(&kc, &fb, &index_of, &st, t, 2).0.is_some()
                            }
                        })
                        .collect();
                    let wall = start.elapsed().as_nanos() as f64;
                    let p = mq_fes_profile();
                    arms.push(serde_json::json!({
                        "arm": arm,
                        "calls": p.calls,
                        "decomposed": found.iter().filter(|&&f| f).count(),
                        "word_ops_per_call": p.word_ops as f64 / p.calls.max(1) as f64,
                        "wall_us_per_call": wall / 1e3 / TARGETS as f64,
                        "build_us_per_call": p.build_ns as f64 / 1e3 / p.calls.max(1) as f64,
                        "walk_us_per_call": p.walk_ns as f64 / 1e3 / p.calls.max(1) as f64,
                        "lift_us_per_call": p.lift_ns as f64 / 1e3 / p.calls.max(1) as f64,
                    }));
                    answers.push(found);
                }
                println!(
                    "{}",
                    serde_json::json!({
                        "degree": n, "curve_a": a, "factor_index": index, "ell": ell,
                        "template_build_us": template_us,
                        "targets": TARGETS, "agree": answers[0] == answers[1], "arms": arms,
                    })
                );
            }
        }
    }
}
