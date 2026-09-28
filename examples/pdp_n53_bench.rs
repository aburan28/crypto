//! **Algebraic PDP at `n = 53`** over a non-invariant factor base.
//!
//! At `n = 53` the order of 2 modulo 53 is 52, so `x^53 − 1` is `(x + 1)`
//! times one irreducible of degree 52 and the only Frobenius-stable
//! subspaces of `F_{2^53}` have dimension 0, 1, 52 or 53: the invariant
//! ladder of `examples/pdp_bench.rs` has no rung here.  This bench uses the
//! classical Faugère–Perret–Petit–Renault / Petit–Quisquater base instead,
//! `F_V = { P : x(P) ∈ V }` for a frozen random `ℓ`-dimensional `V ∋ 1`
//! ([`random_subspace_containing_one`], [`plain_subspace_factor_base`]).
//!
//! ```bash
//! cargo run --release --example pdp_n53_bench -- M ELL ENGINE [PLANTED RANDOM] [BUDGET_S] [NODES]
//! ```
//!
//! `ENGINE` is `groebner` (64-bit monomials), `wide` (`u128`), or
//! `enumerate` (the exact reference: `|F|^{m−1}` group steps).  Planted
//! targets are sums of `m` base points, so a `None` is a miss; random
//! targets are `[k]G`.  Every decomposition returned is re-checked in the
//! group.  A budget stop is recorded as unreached, never as a "no".

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_groebner::{FieldStructure, SolverEngine};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    enumerate_decompose, groebner_decompose, KoblitzCurve,
};
use crypto_lib::cryptanalysis::koblitz_symmetrised::{
    plain_subspace_factor_base, random_subspace_containing_one,
};
use crypto_lib::cryptanalysis::wide_groebner::wide_groebner_decompose;
use num_bigint::BigUint;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let arg = |i: usize, d: &str| a.get(i).cloned().unwrap_or_else(|| d.to_string());
    let m: usize = arg(1, "2").parse().expect("M");
    let ell: usize = arg(2, "10").parse().expect("ELL");
    let engine = arg(3, "wide");
    let planted: usize = arg(4, "2").parse().expect("PLANTED");
    let random: usize = arg(5, "2").parse().expect("RANDOM");
    let budget: f64 = arg(6, "300").parse().expect("BUDGET");
    let nodes: usize = arg(7, "20000").parse().expect("NODES");

    let n = 53u32;
    let kc = KoblitzCurve::new(0, n).expect("K_0/2^53");
    let mut rng = StdRng::seed_from_u64(0x5335_0000 ^ ell as u64);
    let basis = random_subspace_containing_one(&kc, ell, &mut rng);
    let t0 = Instant::now();
    let fb = plain_subspace_factor_base(&kc, &basis).expect("basis spans 2^ell");
    let setup = t0.elapsed().as_secs_f64();
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let index_of = fb.index_map();

    let mut trng = StdRng::seed_from_u64(0x5044_5335 ^ (ell as u64) << 8 ^ m as u64);
    let mut targets: Vec<(bool, BinaryPoint)> = Vec::new();
    for _ in 0..planted {
        let t = (0..m).fold(BinaryPoint::Infinity, |acc, _| {
            kc.add(&acc, &fb.points[trng.gen_range(0..fb.points.len())])
        });
        targets.push((true, t));
    }
    for _ in 0..random {
        let k = BigUint::from(trng.gen_range(1u64..u64::MAX));
        targets.push((false, kc.mul(kc.generator(), &k)));
    }

    println!(
        "n=53 m={m} ell={ell} |F|={} vars={} engine={engine} setup={setup:.2}s",
        fb.points.len(),
        m * ell + m.saturating_sub(2) * n as usize
    );
    let started = Instant::now();
    for (i, (is_planted, target)) in targets.iter().enumerate() {
        if started.elapsed().as_secs_f64() > budget {
            println!("UNREACHED budget {budget}s after {i} targets");
            break;
        }
        if *target == BinaryPoint::Infinity {
            continue;
        }
        let t = Instant::now();
        let (found, note) = match engine.as_str() {
            "groebner" => {
                let (f, s) = groebner_decompose(
                    &kc,
                    &fb,
                    &index_of,
                    &st,
                    target,
                    m,
                    SolverEngine::default(),
                    nodes,
                );
                (f, format!("{s:?}"))
            }
            "wide" => {
                let (f, s) = wide_groebner_decompose(
                    &kc,
                    &fb,
                    &index_of,
                    &st,
                    target,
                    m,
                    nodes.max(1 << 24),
                );
                (f, format!("{s:?}"))
            }
            "enumerate" => (
                enumerate_decompose(&kc, &fb, &index_of, target, m),
                String::new(),
            ),
            other => panic!("unknown engine {other}"),
        };
        if let Some(idxs) = &found {
            let sum = idxs
                .iter()
                .fold(BinaryPoint::Infinity, |s, &k| kc.add(&s, &fb.points[k]));
            assert_eq!(&sum, target, "WRONG decomposition");
        }
        println!(
            "target {i} {} found={} {:.3}s {}",
            if *is_planted { "planted" } else { "random" },
            found.is_some(),
            t.elapsed().as_secs_f64(),
            &note[..note.len().min(220)]
        );
    }
}
