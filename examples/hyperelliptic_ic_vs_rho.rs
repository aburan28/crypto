//! Index calculus vs Pollard rho on the same hyperelliptic Jacobian.
//!
//! One table, one unit: `S = total group operations / sqrt(N)`, with
//! every cost each side incurs inside it — precomputation, relation
//! search, and the linear algebra converted to group-operation
//! equivalents by a factor measured at run time.
//!
//! ```sh
//! cargo run --release --example hyperelliptic_ic_vs_rho
//! ```

use num_bigint::BigUint;
use num_traits::{One, Zero};

use crypto_lib::cryptanalysis::hyperelliptic_ic_bench::{head_to_head, HeadToHeadTrials};
use crypto_lib::cryptanalysis::hyperelliptic_index_calculus::{
    build_factor_base, prime_order_of, subgroup_generator, HecIndexCalculusParams, LinearAlgebra,
    RelationSearch,
};
use crypto_lib::prime_hyperelliptic::{
    brute_force_jac_order_via_lpoly, FpPoly, HyperellipticCurveP, MumfordDivisorP,
};

/// `C : y² = x⁵ + 3x³ + 2x² + x + c` over `F_p`, genus 2.
fn curve_over(p: u64, c: u64) -> HyperellipticCurveP {
    let p = BigUint::from(p);
    let f = FpPoly::from_coeffs(
        vec![
            BigUint::from(c) % &p,
            BigUint::from(1u32),
            BigUint::from(2u32),
            BigUint::from(3u32),
            BigUint::zero(),
            BigUint::from(1u32),
        ],
        p.clone(),
    );
    HyperellipticCurveP::new(p, f, 2)
}

/// `f` squarefree ⟺ `gcd(f, f') = 1`.  A singular curve would run
/// Cantor's algorithm outside its hypotheses and quietly mis-measure
/// both algorithms, so such a `c` is skipped rather than used.
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

/// The comparison is only meaningful in a subgroup large enough that
/// rho's fixed costs do not dominate its walk, so pick, for each `p`,
/// the first `c` whose Jacobian has a prime factor carrying most of the
/// group.  Scanning `c` rather than fixing it is a choice about the
/// *instance*, made before either algorithm runs and reported in the
/// table.
fn pick_instance(p: u64) -> Option<(HyperellipticCurveP, BigUint, BigUint, u64)> {
    for c in 1..p.min(60) {
        let curve = curve_over(p, c);
        if !is_squarefree(&curve) {
            continue;
        }
        let jac = brute_force_jac_order_via_lpoly(&curve);
        let n = largest_prime_factor(&jac);
        if &n * BigUint::from(4u32) >= jac && n > BigUint::from(1000u32) {
            return Some((curve, jac, n, c));
        }
    }
    None
}

fn largest_prime_factor(n: &BigUint) -> BigUint {
    let mut rest = n.clone();
    let mut best = BigUint::one();
    let mut d = BigUint::from(2u32);
    while &d * &d <= rest {
        while (&rest % &d).is_zero() {
            rest /= &d;
            if d > best {
                best = d.clone();
            }
        }
        d += 1u32;
    }
    if rest > best {
        best = rest;
    }
    best
}

fn main() {
    // Squarefree f is required; a p where f has a repeated root is
    // skipped rather than silently mis-measured.
    let primes = [41u64, 61, 101, 151, 211, 251];

    println!(
        "C : y^2 = x^5 + 3x^3 + 2x^2 + x + c over F_p, genus 2; c is the first\n\
         value giving a near-prime Jacobian order (N >= #Jac/4), chosen before\n\
         either algorithm runs.\n\
         Unit: S = total group operations / sqrt(N), N = prime order of D1.\n\
         Index calculus carries its linear algebra, converted at the measured\n\
         mul-mods-per-group-op factor shown in the last column.\n"
    );
    let trials = HeadToHeadTrials::default();
    println!(
        "Each row averages {} index-calculus runs and {} rho runs on the same instance.\n",
        trials.ic, trials.rho
    );

    println!(
        "{:>5} {:>7} {:>8} {:>6} {:>10} {:>10} {:>9} {:>9} {:>8} {:>8} {:>7}",
        "p",
        "search",
        "N",
        "m",
        "IC ops",
        "rho ops",
        "S_ic",
        "S_rho",
        "S_ic/S_rho",
        "S_walk",
        "IC/flr"
    );

    for &p in &primes {
        let (curve, jac, n, c) = match pick_instance(p) {
            Some(v) => v,
            None => {
                println!("{p:>5}   skipped: no c gives a near-prime Jacobian order");
                continue;
            }
        };
        let cofactor = &jac / &n;
        let fb = build_factor_base(&curve, usize::MAX);
        let d1: MumfordDivisorP = match subgroup_generator(&curve, &fb, &jac, &n) {
            Some(d) => d,
            None => {
                println!("{p:>5}   skipped: no generator of the order-{n} subgroup");
                continue;
            }
        };
        if prime_order_of(&curve, &d1, &jac).as_ref() != Some(&n) {
            println!("{p:>5}   skipped: D1 does not have prime order {n}");
            continue;
        }

        let k = &n / BigUint::from(3u32) + BigUint::from(7u32);
        let d2 = d1.scalar_mul(&k, &curve);

        for (label, search, la) in [
            ("rnd+dns", RelationSearch::Random, LinearAlgebra::Dense),
            ("wlk+dns", RelationSearch::walk(), LinearAlgebra::Dense),
            ("rnd+spr", RelationSearch::Random, LinearAlgebra::Sparse),
            ("wlk+spr", RelationSearch::walk(), LinearAlgebra::Sparse),
        ] {
            let params = HecIndexCalculusParams {
                fb_size: usize::MAX,
                extra_relations: 8,
                max_trials: 2_000_000,
                seed: 20260916,
                search,
                linear_algebra: la,
            };
            let row = head_to_head(
                &curve, &d1, &d2, &n, &params, 0xC0FFEE, 50_000_000, &k, &trials,
            );

            // A row without a verified answer on both sides is not a result.
            let mark = match (row.ic_correct, row.rho_correct) {
                (true, true) => "",
                (false, true) => "  [IC UNSOLVED - not a result]",
                (true, false) => "  [rho UNSOLVED - not a result]",
                (false, false) => "  [both UNSOLVED - not a result]",
            };

            println!(
                "{:>5} {:>7} {:>8} {:>6} {:>10.0} {:>10.0} {:>9.2} {:>9.2} {:>8.2} {:>8.2} {:>7.2}{}",
                p,
                label,
                row.n,
                row.factor_base_size,
                row.ic_total_group_ops,
                row.rho_group_ops,
                row.ic_s,
                row.rho_s,
                row.ratio_to_reference(),
                row.rho_walk_s,
                row.ratio_to_floor(),
                mark
            );
            println!(
                "        c = {}, #Jac = {} = {} x {}; relation stage {:.0} ops \
                 ({:.0} precompute, {:.2} ops/trial) + oracle {:.0} mul-mods \
                 ({:.0} equiv) + linear algebra {:.0} mul-mods ({:.0} equiv), conv {:.0}; \
                 smoothness {:.3}; IC floor S = {:.2}; wall {:.0} ms IC vs {:.0} ms rho",
                c,
                jac,
                cofactor,
                n,
                row.ic_relation_ops,
                row.ic_precompute_ops,
                row.ic_ops_per_trial,
                row.ic_oracle_modmuls,
                row.ic_oracle_group_equiv,
                row.ic_la_modmuls,
                row.ic_la_group_equiv,
                row.modmuls_per_group_op,
                row.ic_smoothness_rate,
                row.ic_floor_s,
                row.ic_wall_ms,
                row.rho_wall_ms,
            );
        }
    }

    println!(
        "\nS_walk is rho's walk alone (S_rho also carries the 16-branch\n\
         precomputation, a fixed cost that does not scale with sqrt(N) and so\n\
         inflates S_rho at these toy sizes). The collision itself arrives after\n\
         about sqrt(pi*N/2) = 1.25*sqrt(N) distinct points, but Floyd's cycle\n\
         finding spends 3 group operations per iteration, so S_walk ~ 3 is this\n\
         reference at its own optimum. A distinguished-point rho would cut that\n\
         by roughly 3x -- which widens the gap below, it does not close it.\n\
         S_ic/S_rho below 1 would mean index calculus won that instance."
    );
}
