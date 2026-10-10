//! Parallel Pollard rho with distinguished points (van Oorschot–Wiener),
//! the negation map, and — on Koblitz curves — Frobenius folding.
//!
//! Each worker runs `r`-adding walks `P ← rep(P + R_j)`, `j = key(P) mod r`,
//! where `rep` picks one representative of the class `{±σ^i(P)}` (σ the
//! Frobenius, `i < m`, on a Koblitz curve; the class is just `{±P}`
//! elsewhere).  Coefficients `(a, b)` with `P = [a]G + [b]Q` are carried
//! exactly: the representative equals `[(−1)^s λ^i] P`, so both coefficients
//! are multiplied by that unit mod `n`.  A walk stores its point when the
//! key has `dp_bits` low zero bits, then restarts from a fresh random
//! combination; two stored points that agree with different coefficients
//! give `k`.
//!
//! Fruitless 2-cycles from the negation map are escaped by doubling; walks
//! that run `max_walk_factor · 2^dp_bits` steps without a distinguished
//! point (a longer fruitless cycle) are abandoned.
//!
//! Expected cost is `√(π n / (4·aut))` group additions where `aut` is the
//! folded automorphism order (`1`, or `m` on a Koblitz curve); the
//! negation map alone gives `√(π n / 4)`.

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::SeedableRng;

use super::group::{mix64, EcdlpGroup};

/// Knobs for [`pollard_rho`].
#[derive(Clone, Debug)]
pub struct RhoOptions {
    /// Worker threads (each runs its own walks); `0` means one.
    pub threads: usize,
    /// Seed for the jump table and the starting coefficients.
    pub seed: u64,
    /// Distinguished-point density: a point is stored when `dp_bits` low
    /// bits of a hash of its key are zero.  `None` picks
    /// `⌊log₂(expected)/2⌋ − 1`, clamped to `[0, 40]`.
    pub dp_bits: Option<u32>,
    /// Abort after this many group additions in total (all threads).
    pub max_iterations: Option<u64>,
    /// Number of jumps in the `r`-adding walk.
    pub jumps: usize,
    /// Fold `±σ^i` (negation and Frobenius where available).  `false`
    /// folds nothing — the control arm.
    pub fold_automorphisms: bool,
    /// A walk is abandoned after `max_walk_factor · 2^dp_bits` steps
    /// without reaching a distinguished point.
    pub max_walk_factor: u32,
}

impl Default for RhoOptions {
    fn default() -> Self {
        Self {
            threads: 1,
            seed: 0x5EED_u64,
            dp_bits: None,
            max_iterations: None,
            jumps: 32,
            fold_automorphisms: true,
            max_walk_factor: 20,
        }
    }
}

/// What a rho run did.
#[derive(Clone, Debug)]
pub struct RhoReport {
    /// `k` with `Q = [k]G`, verified by a scalar multiplication.
    pub scalar: Option<BigUint>,
    /// Group additions performed (doublings included), all threads.
    pub iterations: u64,
    /// Distinguished points stored.
    pub distinguished_points: u64,
    /// Walks started (initial plus every restart after a DP or abandonment).
    pub walks: u64,
    /// Walks abandoned for exceeding the length cap.
    pub abandoned_walks: u64,
    /// Fruitless 2-cycles escaped by doubling.
    pub cycle_escapes: u64,
    /// Distinguished-point collisions that carried identical coefficients
    /// and so gave nothing.
    pub useless_collisions: u64,
    pub elapsed_ms: u128,
    /// `√(π n / (4·aut))` for the fold in force (`√(π n / 2)` unfolded).
    pub expected_iterations: f64,
    pub dp_bits: u32,
    /// Size of the folded class: `2·aut` with folding, `1` without.
    pub class_size: u32,
    pub threads: usize,
}

/// Expected additions of a rho with classes of size `class_size`.
pub fn expected_rho_iterations(n: &BigUint, class_size: u32) -> f64 {
    let n_f = biguint_to_f64(n);
    (std::f64::consts::PI * n_f / (2.0 * class_size as f64)).sqrt()
}

pub(crate) fn biguint_to_f64(n: &BigUint) -> f64 {
    // Enough precision for a cost estimate at any size.
    let bits = n.bits();
    if bits <= 53 {
        return n.iter_u64_digits().next().unwrap_or(0) as f64;
    }
    let shift = bits - 53;
    let top = (n >> shift).iter_u64_digits().next().unwrap_or(0) as f64;
    top * 2f64.powi(shift as i32)
}

fn default_dp_bits(expected: f64) -> u32 {
    if expected < 4.0 {
        return 0;
    }
    let half = (expected.log2() / 2.0).floor() as i64 - 1;
    half.clamp(0, 40) as u32
}

struct Entry<E> {
    point: E,
    a: BigUint,
    b: BigUint,
}

struct Shared<E> {
    table: Mutex<HashMap<u64, Entry<E>>>,
    done: AtomicBool,
    iterations: AtomicU64,
    dps: AtomicU64,
    walks: AtomicU64,
    abandoned: AtomicU64,
    escapes: AtomicU64,
    useless: AtomicU64,
    result: Mutex<Option<BigUint>>,
}

/// Solve `Q = [k]G` in the whole group by parallel rho.
pub fn pollard_rho<G: EcdlpGroup>(g: &G, target: &G::Elt, opts: &RhoOptions) -> RhoReport {
    let t0 = Instant::now();
    let n = g.order().clone();
    let gen = g.generator();
    let class_size = if opts.fold_automorphisms {
        2 * g.automorphism_order()
    } else {
        1
    };
    let expected = expected_rho_iterations(&n, class_size);
    let dp_bits = opts.dp_bits.unwrap_or_else(|| default_dp_bits(expected));
    let threads = opts.threads.max(1);
    let mut report = RhoReport {
        scalar: None,
        iterations: 0,
        distinguished_points: 0,
        walks: 0,
        abandoned_walks: 0,
        cycle_escapes: 0,
        useless_collisions: 0,
        elapsed_ms: 0,
        expected_iterations: expected,
        dp_bits,
        class_size,
        threads,
    };

    if g.is_identity(target) {
        report.scalar = Some(BigUint::zero());
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }
    if n.is_one() {
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }

    // Jump table R_j = [a_j]G + [b_j]Q, shared by every walk.
    let r = opts.jumps.max(2);
    let mut rng = StdRng::seed_from_u64(opts.seed);
    let jumps: Vec<(G::Elt, BigUint, BigUint)> = (0..r)
        .map(|_| {
            let a = rng.gen_biguint_below(&n);
            let b = rng.gen_biguint_below(&n);
            let p = g.add(&g.mul(&gen, &a), &g.mul(target, &b));
            (p, a, b)
        })
        .collect();
    let multipliers = if opts.fold_automorphisms {
        g.orbit_multipliers()
    } else {
        vec![BigUint::one()]
    };
    let shared = Shared::<G::Elt> {
        table: Mutex::new(HashMap::new()),
        done: AtomicBool::new(false),
        iterations: AtomicU64::new(0),
        dps: AtomicU64::new(0),
        walks: AtomicU64::new(0),
        abandoned: AtomicU64::new(0),
        escapes: AtomicU64::new(0),
        useless: AtomicU64::new(0),
        result: Mutex::new(None),
    };
    let max_walk: u64 = (opts.max_walk_factor as u64).saturating_mul(1u64 << dp_bits.min(62));
    let dp_mask: u64 = if dp_bits >= 64 {
        u64::MAX
    } else {
        (1u64 << dp_bits) - 1
    };
    let budget = opts.max_iterations.unwrap_or(u64::MAX);

    std::thread::scope(|scope| {
        for t in 0..threads {
            let shared = &shared;
            let jumps = &jumps;
            let multipliers = &multipliers;
            let n = &n;
            let gen = &gen;
            let seed = opts.seed ^ mix64(0xA11C_E000 + t as u64);
            let fold = opts.fold_automorphisms;
            scope.spawn(move || {
                worker(
                    g,
                    target,
                    gen,
                    n,
                    jumps,
                    multipliers,
                    fold,
                    dp_mask,
                    max_walk,
                    budget,
                    seed,
                    shared,
                )
            });
        }
    });

    report.scalar = shared.result.into_inner().unwrap();
    report.iterations = shared.iterations.load(Ordering::Relaxed);
    report.distinguished_points = shared.dps.load(Ordering::Relaxed);
    report.walks = shared.walks.load(Ordering::Relaxed);
    report.abandoned_walks = shared.abandoned.load(Ordering::Relaxed);
    report.cycle_escapes = shared.escapes.load(Ordering::Relaxed);
    report.useless_collisions = shared.useless.load(Ordering::Relaxed);
    report.elapsed_ms = t0.elapsed().as_millis();
    report
}

/// Apply the class multiplier `(−1)^neg · λ^i` to both coefficients.
#[inline]
fn fold_coeffs(a: &mut BigUint, b: &mut BigUint, i: u32, neg: bool, mult: &[BigUint], n: &BigUint) {
    if i != 0 {
        let m = &mult[i as usize];
        *a = &*a * m % n;
        *b = &*b * m % n;
    }
    if neg {
        if !a.is_zero() {
            *a = n - &*a;
        }
        if !b.is_zero() {
            *b = n - &*b;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn worker<G: EcdlpGroup>(
    g: &G,
    target: &G::Elt,
    gen: &G::Elt,
    n: &BigUint,
    jumps: &[(G::Elt, BigUint, BigUint)],
    multipliers: &[BigUint],
    fold: bool,
    dp_mask: u64,
    max_walk: u64,
    budget: u64,
    seed: u64,
    shared: &Shared<G::Elt>,
) {
    let mut rng = StdRng::seed_from_u64(seed);
    let r = jumps.len() as u64;
    let mut local_iters: u64 = 0;
    const FLUSH: u64 = 256;

    'walks: while !shared.done.load(Ordering::Relaxed) {
        // Fresh start: P = [a]G + [b]Q.
        let mut a = rng.gen_biguint_below(n);
        let mut b = rng.gen_biguint_below(n);
        let mut p = g.add(&g.mul(gen, &a), &g.mul(target, &b));
        if fold {
            let (rep, i, neg) = g.canonical_orbit(&p);
            fold_coeffs(&mut a, &mut b, i, neg, multipliers, n);
            p = rep;
        }
        shared.walks.fetch_add(1, Ordering::Relaxed);
        let mut prev: Option<G::Elt> = None;
        let mut steps: u64 = 0;

        loop {
            if g.is_identity(&p) {
                // [a]G + [b]Q = O ⇒ k = −a/b when b ≠ 0.
                if let Some(k) = solve_from(&a, &b, &BigUint::zero(), &BigUint::zero(), n) {
                    if g.mul(gen, &k) == *target {
                        finish(shared, k);
                        break 'walks;
                    }
                }
                continue 'walks;
            }
            let key = g.key(&p);
            if mix64(key ^ 0xD15C_7A9B_E55E_D000) & dp_mask == 0 {
                // Distinguished: store, check, restart.
                shared.dps.fetch_add(1, Ordering::Relaxed);
                let mut table = shared.table.lock().unwrap();
                if let Some(e) = table.get(&key) {
                    if e.point == p {
                        if let Some(k) = solve_from(&a, &b, &e.a, &e.b, n) {
                            if g.mul(gen, &k) == *target {
                                drop(table);
                                finish(shared, k);
                                break 'walks;
                            }
                        } else {
                            shared.useless.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                } else {
                    table.insert(
                        key,
                        Entry {
                            point: p.clone(),
                            a: a.clone(),
                            b: b.clone(),
                        },
                    );
                }
                drop(table);
                continue 'walks;
            }
            if steps >= max_walk {
                shared.abandoned.fetch_add(1, Ordering::Relaxed);
                continue 'walks;
            }

            // One r-adding step.
            let j = (key % r) as usize;
            let (rj, aj, bj) = &jumps[j];
            let sum = g.add(&p, rj);
            let mut a2 = (&a + aj) % n;
            let mut b2 = (&b + bj) % n;
            let mut next = if fold {
                let (rep, i, neg) = g.canonical_orbit(&sum);
                fold_coeffs(&mut a2, &mut b2, i, neg, multipliers, n);
                rep
            } else {
                sum
            };
            local_iters += 1;

            if fold && prev.as_ref() == Some(&next) {
                // Fruitless 2-cycle: escape by doubling the current point.
                shared.escapes.fetch_add(1, Ordering::Relaxed);
                let dbl = g.double(&p);
                a2 = (&a << 1u32) % n;
                b2 = (&b << 1u32) % n;
                let (rep, i, neg) = g.canonical_orbit(&dbl);
                fold_coeffs(&mut a2, &mut b2, i, neg, multipliers, n);
                next = rep;
                local_iters += 1;
            }

            prev = Some(std::mem::replace(&mut p, next));
            a = a2;
            b = b2;
            steps += 1;

            if local_iters >= FLUSH {
                let total =
                    shared.iterations.fetch_add(local_iters, Ordering::Relaxed) + local_iters;
                local_iters = 0;
                if total >= budget {
                    shared.done.store(true, Ordering::Relaxed);
                    break 'walks;
                }
                if shared.done.load(Ordering::Relaxed) {
                    break 'walks;
                }
            }
        }
    }
    shared.iterations.fetch_add(local_iters, Ordering::Relaxed);
}

fn finish<E>(shared: &Shared<E>, k: BigUint) {
    let mut res = shared.result.lock().unwrap();
    if res.is_none() {
        *res = Some(k);
    }
    shared.done.store(true, Ordering::Relaxed);
}

/// From `[a1]G + [b1]Q = [a2]G + [b2]Q` recover `k = (a1 − a2)/(b2 − b1)`.
pub(crate) fn solve_from(
    a1: &BigUint,
    b1: &BigUint,
    a2: &BigUint,
    b2: &BigUint,
    n: &BigUint,
) -> Option<BigUint> {
    let db = (n + b2 - b1) % n;
    if db.is_zero() {
        return None;
    }
    let da = (n + a1 - a2) % n;
    let inv = db.modpow(&(n - BigUint::from(2u32)), n);
    Some(da * inv % n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ecdlp_nist::toy;

    fn check<G: EcdlpGroup>(g: &G, k: u64, opts: &RhoOptions) -> RhoReport {
        let q = g.mul(&g.generator(), &BigUint::from(k));
        let rep = pollard_rho(g, &q, opts);
        assert_eq!(
            rep.scalar.as_ref(),
            Some(&BigUint::from(k)),
            "{}: rho failed after {} iterations ({rep:?})",
            g.name(),
            rep.iterations
        );
        rep
    }

    #[test]
    fn rho_recovers_on_toy_prime_curves() {
        let opts = RhoOptions::default();
        check(&toy::prime_small(), 7777, &opts);
        let rep = check(&toy::prime_mid(), 65_432, &opts);
        assert!(rep.iterations < 60 * rep.expected_iterations as u64 + 2000);
    }

    #[test]
    fn rho_recovers_on_toy_koblitz_with_and_without_frobenius() {
        let g = toy::koblitz(23, 1).expect("toy");
        let opts = RhoOptions {
            threads: 2,
            ..RhoOptions::default()
        };
        let folded = check(&g, 1_234_567, &opts);
        assert_eq!(folded.class_size, 2 * 23);
        let plain = RhoOptions {
            fold_automorphisms: false,
            ..opts.clone()
        };
        let unfolded = check(&g, 1_234_567, &plain);
        assert_eq!(unfolded.class_size, 1);
        assert!(unfolded.expected_iterations > folded.expected_iterations * 6.0);
    }

    #[test]
    fn rho_handles_trivial_targets_and_budgets() {
        let g = toy::prime_small();
        let rep = pollard_rho(&g, &g.identity(), &RhoOptions::default());
        assert_eq!(rep.scalar, Some(BigUint::zero()));
        let q = g.mul(&g.generator(), &BigUint::from(5u32));
        let rep = pollard_rho(
            &g,
            &q,
            &RhoOptions {
                max_iterations: Some(1),
                dp_bits: Some(30),
                ..RhoOptions::default()
            },
        );
        assert!(rep.scalar.is_none() || rep.scalar == Some(BigUint::from(5u32)));
        assert!(rep.iterations <= 300);
    }
}
