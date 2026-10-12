//! Round eight: a reduced factor base with single large primes, against
//! the full base and against Pollard rho, on the same Jacobians.
//!
//! Protocol (registered before this file existed):
//! `research/notes/index-calculus/RESEARCH_HYPERELLIPTIC_IC_RHO.md`,
//! "Round eight".  Curves and instance selection are round seven's
//! (`examples/hyperelliptic_ic_vs_rho.rs`), copied rather than shared so
//! that the frozen round-seven command stays byte-for-byte what it was.
//!
//! Each row prints a human-readable line and a `JSON {...}` line; the
//! analysis reads the latter.
//!
//! ```sh
//! cargo run --release --example hyperelliptic_ic_large_prime
//! ```

use std::io::{self, Write};

use num_bigint::BigUint;
use num_traits::{One, Zero};

use crypto_lib::cryptanalysis::hyperelliptic_ic_bench::{
    head_to_head, HeadToHeadTrials, RhoVariant,
};
use crypto_lib::cryptanalysis::hyperelliptic_index_calculus::{
    build_factor_base, divisor_order_bsgs, hasse_interval, largest_prime_factor,
    HecIndexCalculusParams, LinearAlgebra, RelationSearch, SmoothnessTest,
};
use crypto_lib::prime_hyperelliptic::{FpPoly, HyperellipticCurveP, MumfordDivisorP};

/// Genus 2: `y² = x⁵ + 3x³ + 2x² + x + c`.  Genus 3: `y² = x⁷ + x³ + cx + 1`.
/// Genus 4: `y² = x⁹ + x³ + cx + 1`.  As in round seven.
fn curve_over(p: u64, c: u64, genus: u32) -> HyperellipticCurveP {
    let p = BigUint::from(p);
    let coeffs = match genus {
        2 => vec![
            BigUint::from(c) % &p,
            BigUint::from(1u32),
            BigUint::from(2u32),
            BigUint::from(3u32),
            BigUint::zero(),
            BigUint::from(1u32),
        ],
        3 => vec![
            BigUint::one(),
            BigUint::from(c) % &p,
            BigUint::zero(),
            BigUint::one(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::one(),
        ],
        4 => vec![
            BigUint::one(),
            BigUint::from(c) % &p,
            BigUint::zero(),
            BigUint::one(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::one(),
        ],
        g => panic!("no model for genus {g}; degree must be 2g+1"),
    };
    let f = FpPoly::from_coeffs(coeffs, p.clone());
    HyperellipticCurveP::new(p, f, genus)
}

fn is_squarefree(curve: &HyperellipticCurveP) -> bool {
    let p = &curve.p;
    let d: Vec<BigUint> = curve
        .f
        .coeffs
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| (c * BigUint::from(i as u64)) % p)
        .collect();
    let deriv = FpPoly::from_coeffs(d, p.clone());
    if deriv.is_zero() {
        return false;
    }
    curve.f.gcd(&deriv).degree() == Some(0)
}

/// Round seven's instance rule: the first `c` whose Jacobian has a prime
/// subgroup above 1000 carrying at least a quarter of the Hasse-Weil
/// lower bound.  Chosen before either algorithm runs.
fn pick_instance(
    p: u64,
    genus: u32,
) -> Option<(HyperellipticCurveP, MumfordDivisorP, BigUint, u64)> {
    for c in 1..p.min(60) {
        let curve = curve_over(p, c, genus);
        if !is_squarefree(&curve) {
            continue;
        }
        let fb = build_factor_base(&curve, usize::MAX);
        for e in fb.entries.iter().take(8) {
            let order = match divisor_order_bsgs(&curve, &e.divisor, 200_000) {
                Some(o) => o,
                None => continue,
            };
            let l = largest_prime_factor(&order);
            let (lo, _hi) = hasse_interval(&curve)?;
            if l < BigUint::from(1000u32) || &l * BigUint::from(4u32) < BigUint::from(lo) {
                continue;
            }
            let d1 = e.divisor.scalar_mul(&(&order / &l), &curve);
            if divisor_order_bsgs(&curve, &d1, 200_000).as_ref() != Some(&l) {
                continue;
            }
            return Some((curve, d1, l, c));
        }
    }
    None
}

fn main() {
    let sweep: [(u32, &[u64]); 3] = [
        (2, &[101, 251, 503, 1009, 2003, 4001]),
        (3, &[61, 211, 401]),
        (4, &[31, 61]),
    ];
    // (θ, large primes).  θ = 1 is round seven's configuration; large
    // primes are inert there.  θ = 0.5 without them is the control.
    let configs: [(f64, bool); 5] = [
        (1.0, true),
        (0.85, true),
        (0.7, true),
        (0.5, true),
        (0.5, false),
    ];
    // Post-hoc diagnostics only (round eight's registered sweep uses the
    // defaults): `LP_ONLY=g:p` restricts the sweep to one instance and
    // `LP_IC_TRIALS=n` raises the index-calculus seed count, so seed
    // variance can be separated from a real effect on one row.
    let only: Option<(u32, u64)> = std::env::var("LP_ONLY").ok().and_then(|v| {
        let (g, p) = v.split_once(':')?;
        Some((g.parse().ok()?, p.parse().ok()?))
    });
    let mut trials = HeadToHeadTrials::default();
    if let Some(n) = std::env::var("LP_IC_TRIALS")
        .ok()
        .and_then(|v| v.parse().ok())
    {
        trials.ic = n;
    }
    println!(
        "Round eight. Unit S = group operations / sqrt(N), in-situ conversion; \
         {} IC and {} DP-rho runs per row.",
        trials.ic, trials.rho
    );
    println!(
        "{:>2} {:>5} {:>10} {:>5} {:>5} {:>3} {:>5} {:>8} {:>8} {:>7} {:>7} {:>6} {:>6}",
        "g", "p", "N", "m", "theta", "LP", "fb", "S_ic", "S_rho", "ratio", "solve", "part", "comb"
    );

    for (genus, primes) in sweep {
        for &p in primes {
            if only.is_some_and(|o| o != (genus, p)) {
                continue;
            }
            io::stdout().flush().ok();
            let (curve, d1, n, c) = match pick_instance(p, genus) {
                Some(v) => v,
                None => {
                    println!("{genus:>2} {p:>5}   skipped: no usable prime-order subgroup");
                    continue;
                }
            };
            let m_full = build_factor_base(&curve, usize::MAX).len();
            let k = &n / BigUint::from(3u32) + BigUint::from(7u32);
            let d2 = d1.scalar_mul(&k, &curve);

            for &(theta, lp) in &configs {
                let fb_size = if theta >= 1.0 {
                    usize::MAX
                } else {
                    ((theta * m_full as f64).ceil() as usize).max(genus as usize + 1)
                };
                let params = HecIndexCalculusParams {
                    fb_size,
                    extra_relations: 8,
                    max_trials: 5_000_000,
                    seed: 20260916,
                    search: RelationSearch::factor_base_walk(),
                    linear_algebra: LinearAlgebra::Sparse,
                    smoothness: SmoothnessTest::Gcd,
                    large_primes: lp,
                };
                let row = head_to_head(
                    &curve,
                    &d1,
                    &d2,
                    &n,
                    &params,
                    0xC0FFEE,
                    200_000_000,
                    &k,
                    &trials,
                    &RhoVariant::DistinguishedPoints,
                );
                let mark = if row.ic_correct && row.rho_correct {
                    ""
                } else {
                    "  [UNSOLVED - not a result]"
                };
                println!(
                    "{:>2} {:>5} {:>10} {:>5} {:>5.2} {:>3} {:>5} {:>8.3} {:>8.3} {:>7.3} {:>7.0} {:>6.0} {:>6.0}{}",
                    genus,
                    p,
                    row.n,
                    m_full,
                    theta,
                    if lp { "on" } else { "off" },
                    row.factor_base_size,
                    row.ic_s,
                    row.rho_s,
                    row.ratio_to_reference(),
                    row.ic_la_group_equiv,
                    row.ic_partial_relations,
                    row.ic_combined_relations,
                    mark
                );
                println!(
                    "JSON {{\"genus\":{},\"p\":{},\"c\":{},\"N\":{},\"m_full\":{},\"theta\":{},\"large_primes\":{},\
                     \"fb\":{},\"S_ic\":{},\"S_rho\":{},\"S_walk\":{},\"ratio\":{},\"relation_ops\":{},\
                     \"precompute_ops\":{},\"oracle_equiv\":{},\"oracle_modmuls\":{},\"la_equiv\":{},\"la_modmuls\":{},\
                     \"partial_equiv\":{},\"partials\":{},\"combined\":{},\"smoothness\":{},\"floor_S\":{},\
                     \"wall_ic_ms\":{},\"wall_rho_ms\":{},\"rho_ops\":{},\"ic_correct\":{},\"rho_correct\":{}}}",
                    genus,
                    p,
                    c,
                    row.n,
                    m_full,
                    theta,
                    lp,
                    row.factor_base_size,
                    row.ic_s,
                    row.rho_s,
                    row.rho_walk_s,
                    row.ratio_to_reference(),
                    row.ic_relation_ops,
                    row.ic_precompute_ops,
                    row.ic_oracle_group_equiv,
                    row.ic_oracle_modmuls,
                    row.ic_la_group_equiv,
                    row.ic_la_modmuls,
                    row.ic_partial_group_equiv,
                    row.ic_partial_relations,
                    row.ic_combined_relations,
                    row.ic_smoothness_rate,
                    row.ic_floor_s,
                    row.ic_wall_ms,
                    row.rho_wall_ms,
                    row.rho_group_ops,
                    row.ic_correct,
                    row.rho_correct,
                );
                io::stdout().flush().ok();
            }
        }
    }
}
