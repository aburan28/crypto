//! # Index calculus on hyperelliptic Jacobians over `F_p`.
//!
//! The Adleman–DeMarrais–Huang / Gaudry attack on the discrete
//! logarithm problem in `Jac(C)(F_p)` for `C : y² = f(x)` of genus
//! `g`, built on the Mumford / Cantor substrate in
//! [`crate::prime_hyperelliptic`].
//!
//! ## The algorithm
//!
//! Write `∞` for the point at infinity (`deg f = 2g+1`, so there is
//! exactly one).  Every degree-1 place `P = (x_0, y_0) ∈ C(F_p)`
//! gives a divisor class `[P] − [∞]`, whose Mumford representation is
//! `(x − x_0, y_0)`.  Those classes are the **factor base**
//! `F = {F_1, …, F_m}`; `m ≤ #C(F_p) ≈ p` after the `±` identification
//! `[(x, −y)] = −[(x, y)]`.
//!
//! 1. **Relation search.**  Draw `a, b` uniformly mod `N` and reduce
//!    `R = a·D_1 + b·D_2` with Cantor.  The reduced representative is
//!    `(u, v)` with `deg u ≤ g`, and the class is
//!    `Σ_i m_i ([P_i] − [∞])` exactly when `u` splits into linear
//!    factors over `F_p`, `u = Π (x − x_i)^{m_i}`, with
//!    `P_i = (x_i, v(x_i))`.  That event is *smoothness*; its
//!    probability over a random class is `≈ 1/g!` for a factor base
//!    holding every degree-1 place, which is why the attack is
//!    interesting for fixed small `g` and large `p` and is *not* a
//!    threat to genus 1.
//! 2. **Linear algebra.**  Each smooth `R` yields
//!    `Σ_i c_i · log_{D_1} F_i − b · log_{D_1} D_2 ≡ a (mod N)`.
//!    Collect `m + extra` such rows and solve mod the prime order `N`
//!    of `D_1` ([`gaussian_eliminate_mod_n`]).
//! 3. **Read off** `log_{D_1} D_2` and verify `D_2 = k·D_1` before
//!    returning it.  An unverified `k` is never returned.
//!
//! ## Cost, and what this module is for
//!
//! With the full degree-1 factor base the expected work is
//! `O(p² · polylog)` group operations for the relation stage plus a
//! dense `O(m³)` solve here, against `O(p^{g/2})` for Pollard rho on a
//! group of size `≈ p^g`.  So the crossover is at `g = 4`
//! asymptotically, and Gaudry's variant — a *reduced* factor base of
//! size `p^{2/(g+1)}`, which this module supports through
//! `fb_size` — moves it to `g = 3` with an `O(p^{2 − 2/g})` relation
//! stage.  Nothing here is competitive at cryptographic sizes: the
//! implementation is dense, the root finding is `O(p)` per candidate
//! by evaluation, and the whole pipeline is meant to be *read* and to
//! be *checked against rho on the same instances*, in the unit
//! `S = total group operations / sqrt(N)` that `AGENTS.md` fixes.
//! [`HecIndexCalculusReport`] carries the operation counts needed to
//! fill that column; it measures, it does not claim.
//!
//! ## Scope
//!
//! - Odd characteristic, `deg f = 2g+1`, `f` squarefree — the
//!   "imaginary" model, one place at infinity.
//! - `p` small enough to enumerate (`p` must fit in a `u64`; the
//!   factor-base build and the root finding are both `O(p)`).
//! - `N` = order of `D_1` must be **prime**: the solver inverts
//!   pivots mod `N`.  [`prime_order_of`] computes it for toy
//!   Jacobians and refuses a composite answer.

use std::collections::HashMap;

use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::{RngCore, SeedableRng};

use crate::cryptanalysis::ec_index_calculus::{gaussian_eliminate_mod_n, sqrt_mod_p};
use crate::prime_hyperelliptic::{FpPoly, HyperellipticCurveP, MumfordDivisorP};

// ── Factor base ────────────────────────────────────────────────────────

/// One degree-1 place `P = (x, y)` of `C`, held as the divisor class
/// `[P] − [∞]`.
///
/// Only one of `±P` is stored: the representative has
/// `y ≤ (p−1)/2`, and a decomposition that meets `(x, p−y)` records
/// the coefficient `−1` against this entry instead.  Weierstrass
/// points (`y = 0`) are their own negative and always carry `+1`.
#[derive(Clone, Debug)]
pub struct HecFactorBaseEntry {
    pub x: BigUint,
    pub y: BigUint,
    pub divisor: MumfordDivisorP,
}

/// The factor base, plus the `x ↦ index` map the decomposition needs.
#[derive(Clone, Debug)]
pub struct HecFactorBase {
    pub entries: Vec<HecFactorBaseEntry>,
    by_x: HashMap<Vec<u32>, usize>,
}

impl HecFactorBase {
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Index of the entry with this `x`-coordinate, if present.
    pub fn index_of_x(&self, x: &BigUint) -> Option<usize> {
        self.by_x.get(&x.to_u32_digits()).copied()
    }
}

/// Build the factor base: the places `(x, y)` with `y ≤ (p−1)/2`,
/// in increasing `x`, truncated to `max_size` entries.
///
/// Truncating is Gaudry's reduced factor base: relations whose
/// decomposition leaves the truncated set are simply discarded, so a
/// smaller base costs more trials per relation but a smaller solve.
/// `max_size = usize::MAX` keeps every place.
///
/// At most one of `±P` appears, because both would make the relation
/// matrix carry a column and its negative.
pub fn build_factor_base(curve: &HyperellipticCurveP, max_size: usize) -> HecFactorBase {
    let p = &curve.p;
    let p_u = p
        .to_u64_digits()
        .first()
        .copied()
        .filter(|_| p.to_u64_digits().len() <= 1)
        .expect("build_factor_base: p must fit in u64");
    let half = (p - BigUint::one()) >> 1;

    let mut entries = Vec::new();
    let mut by_x = HashMap::new();
    for xi in 0..p_u {
        if entries.len() >= max_size {
            break;
        }
        let x = BigUint::from(xi);
        let y_sq = curve.f.eval(&x);
        let y = match sqrt_mod_p(&y_sq, p) {
            Some(y) => y,
            None => continue,
        };
        // Canonical representative of the ± pair.
        let y = if y > half { p - &y } else { y };
        let divisor = mumford_of_place(curve, &x, &y);
        by_x.insert(x.to_u32_digits(), entries.len());
        entries.push(HecFactorBaseEntry { x, y, divisor });
    }
    HecFactorBase { entries, by_x }
}

/// `[P] − [∞]` in Mumford form: `(x − x_0, y_0)`.
///
/// Unlike [`MumfordDivisorP::from_point`] this also accepts the
/// Weierstrass points `y_0 = 0`, whose class `(x − x_0, 0)` is the
/// 2-torsion element the factor base needs a column for.
fn mumford_of_place(curve: &HyperellipticCurveP, x0: &BigUint, y0: &BigUint) -> MumfordDivisorP {
    debug_assert!(curve.is_on_curve(x0, y0));
    let p = &curve.p;
    let u = FpPoly::from_coeffs(vec![(p - x0) % p, BigUint::one()], p.clone());
    let v = FpPoly::constant(y0.clone(), p.clone());
    MumfordDivisorP { u, v }
}

// ── Smoothness test / decomposition ────────────────────────────────────

/// Split `u` into linear factors over `F_p` by evaluating at every
/// `x ∈ F_p`, returning `(root, multiplicity)` pairs.
///
/// `None` when `u` does **not** split completely — the multiplicities
/// found do not account for `deg u`, so some irreducible factor has
/// degree `≥ 2` and the divisor is not smooth over degree-1 places.
///
/// `O(p · deg u)`; the point of the module is the relation structure,
/// not a Cantor–Zassenhaus.
pub fn split_into_linear_factors(u: &FpPoly, p: &BigUint) -> Option<Vec<(BigUint, usize)>> {
    let deg = u.degree()?;
    if deg == 0 {
        return Some(Vec::new());
    }
    let p_u = *p.to_u64_digits().first()? as u128;
    if p.to_u64_digits().len() > 1 {
        return None;
    }

    let mut roots = Vec::new();
    let mut total = 0usize;
    let mut rest = u.clone();
    for xi in 0..p_u as u64 {
        if total == deg {
            break;
        }
        let x = BigUint::from(xi);
        if !rest.eval(&x).is_zero() {
            continue;
        }
        // Divide out (x − xi) as often as it goes.
        let lin = FpPoly::from_coeffs(vec![(p - &x) % p, BigUint::one()], p.clone());
        let mut mult = 0usize;
        loop {
            let (q, r) = rest.divrem(&lin);
            if !r.is_zero() {
                break;
            }
            rest = q;
            mult += 1;
        }
        total += mult;
        roots.push((x, mult));
    }
    if total == deg {
        Some(roots)
    } else {
        None
    }
}

/// Decompose a reduced divisor over the factor base.
///
/// Returns the coefficient vector `[(index, coefficient)]` with
/// `D ≡ Σ coefficient · F_index` in `Jac(C)(F_p)`, or `None` when `D`
/// is not smooth, or is smooth but meets a place outside a truncated
/// factor base.
pub fn decompose_over_factor_base(
    curve: &HyperellipticCurveP,
    d: &MumfordDivisorP,
    fb: &HecFactorBase,
) -> Option<Vec<(usize, i64)>> {
    let p = &curve.p;
    let half = (p - BigUint::one()) >> 1;
    let roots = split_into_linear_factors(&d.u, p)?;

    let mut acc: HashMap<usize, i64> = HashMap::new();
    for (x, mult) in roots {
        let y = d.v.eval(&x);
        // `u | v² − f` guarantees this, but the attack is only as
        // sound as the representation it reads.
        debug_assert!(curve.is_on_curve(&x, &y));
        let idx = fb.index_of_x(&x)?;
        let sign: i64 = if y.is_zero() || y <= half { 1 } else { -1 };
        debug_assert_eq!(fb.entries[idx].y, if sign > 0 { y.clone() } else { p - &y });
        *acc.entry(idx).or_insert(0) += sign * mult as i64;
    }

    let mut out: Vec<(usize, i64)> = acc.into_iter().filter(|&(_, c)| c != 0).collect();
    out.sort_unstable_by_key(|&(i, _)| i);
    Some(out)
}

// ── Relations ──────────────────────────────────────────────────────────

/// One relation `a·D_1 + b·D_2 = Σ c_i F_i`.
#[derive(Clone, Debug)]
pub struct HecRelation {
    pub coef_a: BigUint,
    pub coef_b: BigUint,
    /// `(factor-base index, coefficient)`, coefficients non-zero.
    pub entries: Vec<(usize, i64)>,
}

/// Counters for the `S = operations / sqrt(N)` accounting `AGENTS.md`
/// requires of any comparison against Pollard rho.
///
/// `jacobian_ops` counts Cantor compositions and doublings charged to
/// the relation stage; the linear-algebra stage is reported as its
/// own row-operation count, because the two are not the same unit and
/// conflating them is exactly how a cost gets relabelled rather than
/// removed.
#[derive(Clone, Debug, Default)]
pub struct HecIndexCalculusReport {
    pub factor_base_size: usize,
    pub trials: usize,
    pub smooth_trials: usize,
    /// Smooth, but met a place outside a truncated factor base.
    pub discarded_off_base: usize,
    pub relations: usize,
    pub jacobian_ops: usize,
    pub solve_row_ops: usize,
    pub solved: bool,
}

impl HecIndexCalculusReport {
    /// Observed smoothness probability — compare against `1/g!` for
    /// the full degree-1 base.
    pub fn smoothness_rate(&self) -> f64 {
        if self.trials == 0 {
            0.0
        } else {
            self.smooth_trials as f64 / self.trials as f64
        }
    }
}

/// Parameters of the relation search and solve.
#[derive(Clone, Debug)]
pub struct HecIndexCalculusParams {
    /// Factor base cap; `usize::MAX` for every degree-1 place.
    pub fb_size: usize,
    /// Rows collected beyond `fb_size + 1` unknowns.
    pub extra_relations: usize,
    /// Give up after this many `(a, b)` draws in total.
    pub max_trials: usize,
    pub seed: u64,
}

impl Default for HecIndexCalculusParams {
    fn default() -> Self {
        Self {
            fb_size: usize::MAX,
            extra_relations: 8,
            max_trials: 200_000,
            seed: 0,
        }
    }
}

/// Uniform-ish `BigUint` in `[0, n)`.
///
/// Samples 64 bits beyond `n` and reduces, so the modulo bias is
/// below `2^-64` — irrelevant to a relation search, and stated rather
/// than hidden because this is not a sampler to reuse for key
/// material.
fn rand_below(rng: &mut StdRng, n: &BigUint) -> BigUint {
    let bytes = (n.bits() as usize).div_ceil(8) + 8;
    let mut buf = vec![0u8; bytes];
    rng.fill_bytes(&mut buf);
    BigUint::from_bytes_be(&buf) % n
}

/// Collect relations by walking `R = a·D_1 + b·D_2` for random
/// `(a, b)` and keeping the smooth reduced representatives.
///
/// Stops at `wanted` relations or `params.max_trials` draws,
/// whichever comes first; the report says which.
#[allow(clippy::too_many_arguments)]
pub fn collect_relations(
    curve: &HyperellipticCurveP,
    d1: &MumfordDivisorP,
    d2: &MumfordDivisorP,
    n: &BigUint,
    fb: &HecFactorBase,
    wanted: usize,
    params: &HecIndexCalculusParams,
    report: &mut HecIndexCalculusReport,
) -> Vec<HecRelation> {
    let mut rng = StdRng::seed_from_u64(params.seed);
    let mut relations = Vec::with_capacity(wanted);

    while relations.len() < wanted && report.trials < params.max_trials {
        report.trials += 1;
        let a = rand_below(&mut rng, n);
        let b = rand_below(&mut rng, n);
        // Charged as ~2·log₂ N doublings/additions: the cost the
        // rho comparison has to see.
        report.jacobian_ops += 2 * (n.bits() as usize + 1);
        let r = d1
            .scalar_mul(&a, curve)
            .add(&d2.scalar_mul(&b, curve), curve);
        if r.is_identity() {
            // a·D₁ + b·D₂ = 0 is a relation with no factor-base part;
            // it is still a valid row and pins the logarithm directly.
            relations.push(HecRelation {
                coef_a: a,
                coef_b: b,
                entries: Vec::new(),
            });
            report.smooth_trials += 1;
            continue;
        }
        match decompose_over_factor_base(curve, &r, fb) {
            Some(entries) => {
                report.smooth_trials += 1;
                relations.push(HecRelation {
                    coef_a: a,
                    coef_b: b,
                    entries,
                });
            }
            None => {
                // Distinguish "not smooth" from "smooth off base" so a
                // truncated base can be tuned on evidence.
                if split_into_linear_factors(&r.u, &curve.p).is_some() {
                    report.smooth_trials += 1;
                    report.discarded_off_base += 1;
                }
            }
        }
    }
    report.relations = relations.len();
    relations
}

// ── End-to-end DLP ─────────────────────────────────────────────────────

/// **Solve `D_2 = k · D_1` in `Jac(C)(F_p)`** by index calculus.
///
/// `n` must be the order of `D_1` and must be **prime** — the solve
/// inverts pivots modulo it.  Use [`prime_order_of`] to obtain one.
///
/// Returns `Some(k)` only after verifying `k·D_1 = D_2`: an
/// under-determined system yields `None`, never an unchecked answer.
pub fn hec_index_calculus_dlp(
    curve: &HyperellipticCurveP,
    d1: &MumfordDivisorP,
    d2: &MumfordDivisorP,
    n: &BigUint,
    params: &HecIndexCalculusParams,
) -> (Option<BigUint>, HecIndexCalculusReport) {
    let mut report = HecIndexCalculusReport::default();
    let fb = build_factor_base(curve, params.fb_size);
    report.factor_base_size = fb.len();
    if fb.is_empty() {
        return (None, report);
    }

    let m = fb.len();
    let wanted = m + 1 + params.extra_relations;
    let relations = collect_relations(curve, d1, d2, n, &fb, wanted, params, &mut report);
    if relations.len() < m + 1 {
        return (None, report);
    }

    // Unknowns: (y_1, …, y_m, k), where y_i = log_{D₁} F_i.
    // Row: Σ c_i y_i − b k ≡ a  (mod n).
    let mut matrix: Vec<Vec<BigUint>> = Vec::with_capacity(relations.len());
    let mut rhs: Vec<BigUint> = Vec::with_capacity(relations.len());
    for rel in &relations {
        let mut row = vec![BigUint::zero(); m + 1];
        for &(j, c) in &rel.entries {
            row[j] = (&row[j] + signed_mod(c, n)) % n;
        }
        row[m] = (n - &(&rel.coef_b % n)) % n;
        matrix.push(row);
        rhs.push(&rel.coef_a % n);
    }
    report.solve_row_ops = relations.len() * (m + 1) * (m + 1);

    let solution = match gaussian_eliminate_mod_n(&mut matrix, &mut rhs, n) {
        Some(s) => s,
        None => return (None, report),
    };
    let k = solution[m].clone();
    if &d1.scalar_mul(&k, curve) == d2 {
        report.solved = true;
        (Some(k), report)
    } else {
        (None, report)
    }
}

fn signed_mod(v: i64, n: &BigUint) -> BigUint {
    if v >= 0 {
        BigUint::from(v as u64) % n
    } else {
        (n - BigUint::from(v.unsigned_abs()) % n) % n
    }
}

// ── Order of a divisor class ───────────────────────────────────────────

/// Trial-division factorisation of a toy Jacobian order.
fn factorise_small(mut n: BigUint) -> Vec<(BigUint, u32)> {
    let mut out = Vec::new();
    let mut d = BigUint::from(2u32);
    while &d * &d <= n {
        let mut e = 0;
        while (&n % &d).is_zero() {
            n /= &d;
            e += 1;
        }
        if e > 0 {
            out.push((d.clone(), e));
        }
        d += 1u32;
    }
    if n > BigUint::one() {
        out.push((n, 1));
    }
    out
}

/// Exact order of `d` in `Jac(C)(F_p)`, given the group order
/// `jac_order`, and **only if that order is prime**.
///
/// Returns `None` when the order is composite: the index-calculus
/// solve needs a field, and a composite modulus belongs to
/// Pohlig–Hellman over the prime factors, not here.
pub fn prime_order_of(
    curve: &HyperellipticCurveP,
    d: &MumfordDivisorP,
    jac_order: &BigUint,
) -> Option<BigUint> {
    let mut order = jac_order.clone();
    for (q, e) in factorise_small(jac_order.clone()) {
        for _ in 0..e {
            let candidate = &order / &q;
            if d.scalar_mul(&candidate, curve).is_identity() {
                order = candidate;
            } else {
                break;
            }
        }
    }
    let factors = factorise_small(order.clone());
    if factors.len() == 1 && factors[0].1 == 1 && order > BigUint::one() {
        Some(order)
    } else {
        None
    }
}

/// Pick a generator of the (unique) subgroup of prime order `l`:
/// multiply a place by the cofactor until the result is non-trivial.
///
/// Returns `None` if every factor-base place lands in the identity,
/// which for `l | #Jac` cannot happen unless `l ∤ #Jac`.
pub fn subgroup_generator(
    curve: &HyperellipticCurveP,
    fb: &HecFactorBase,
    jac_order: &BigUint,
    l: &BigUint,
) -> Option<MumfordDivisorP> {
    let cofactor = jac_order / l;
    for e in &fb.entries {
        let d = e.divisor.scalar_mul(&cofactor, curve);
        if !d.is_identity() {
            return Some(d);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prime_hyperelliptic::brute_force_jac_order_via_lpoly;

    /// `y² = x⁵ + 3x³ + 2x² + x + 1` over `F_p`, genus 2.
    fn toy_curve(p: u64) -> HyperellipticCurveP {
        let p = BigUint::from(p);
        let f = FpPoly::from_coeffs(
            vec![
                BigUint::from(1u32),
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

    #[test]
    fn factor_base_holds_one_of_each_pm_pair() {
        let curve = toy_curve(41);
        let fb = build_factor_base(&curve, usize::MAX);
        assert!(!fb.is_empty());
        let half = (&curve.p - BigUint::one()) >> 1;
        for e in &fb.entries {
            assert!(curve.is_on_curve(&e.x, &e.y));
            assert!(e.y <= half, "non-canonical representative in factor base");
            assert_eq!(fb.index_of_x(&e.x).map(|i| &fb.entries[i].x), Some(&e.x));
        }
    }

    #[test]
    fn splitting_agrees_with_reconstruction() {
        let p = BigUint::from(41u32);
        // (x − 3)²(x − 7)
        let lin = |c: u32| {
            FpPoly::from_coeffs(vec![(&p - BigUint::from(c)) % &p, BigUint::one()], p.clone())
        };
        let u = lin(3).mul(&lin(3)).mul(&lin(7));
        let roots = split_into_linear_factors(&u, &p).expect("splits");
        assert_eq!(
            roots,
            vec![(BigUint::from(3u32), 2), (BigUint::from(7u32), 1)]
        );
        // x² + 1 is irreducible mod 41? 41 ≡ 1 (mod 4) so it splits;
        // use a genuinely irreducible quadratic instead: x² − r for a
        // non-residue r.
        let mut non_residue = None;
        for r in 2u32..41 {
            if sqrt_mod_p(&BigUint::from(r), &p).is_none() {
                non_residue = Some(r);
                break;
            }
        }
        let r = non_residue.expect("a non-residue exists");
        let irred = FpPoly::from_coeffs(
            vec![(&p - BigUint::from(r)) % &p, BigUint::zero(), BigUint::one()],
            p.clone(),
        );
        assert!(split_into_linear_factors(&irred, &p).is_none());
    }

    #[test]
    fn decomposition_reproduces_the_divisor() {
        let curve = toy_curve(41);
        let fb = build_factor_base(&curve, usize::MAX);
        // Two distinct places: `deg u = 2 = g`, so Cantor composes
        // without reducing and the class is smooth by construction.
        // Three would reduce, and a reduced representative of a sum of
        // smooth classes is not itself smooth in general — that is the
        // whole reason the relation search is a search.
        let d = fb.entries[0].divisor.add(&fb.entries[1].divisor, &curve);
        let entries = decompose_over_factor_base(&curve, &d, &fb).expect("smooth by construction");
        let mut rebuilt = MumfordDivisorP::identity(curve.p.clone());
        for (i, c) in entries {
            let base = &fb.entries[i].divisor;
            let term = if c >= 0 {
                base.scalar_mul(&BigUint::from(c as u64), &curve)
            } else {
                base.neg(&curve)
                    .scalar_mul(&BigUint::from(c.unsigned_abs()), &curve)
            };
            rebuilt = rebuilt.add(&term, &curve);
        }
        assert_eq!(rebuilt, d);
    }

    #[test]
    fn solves_a_toy_jacobian_dlp() {
        let curve = toy_curve(41);
        let jac = brute_force_jac_order_via_lpoly(&curve);
        let fb = build_factor_base(&curve, usize::MAX);

        // Largest prime factor of #Jac, and a generator of that
        // subgroup.
        let l = factorise_small(jac.clone())
            .into_iter()
            .map(|(q, _)| q)
            .max()
            .expect("non-trivial Jacobian");
        assert!(l > BigUint::from(50u32), "toy subgroup too small: {l}");
        let d1 = subgroup_generator(&curve, &fb, &jac, &l).expect("generator exists");
        assert_eq!(prime_order_of(&curve, &d1, &jac).as_ref(), Some(&l));

        let k = &l / BigUint::from(3u32) + BigUint::from(7u32);
        let d2 = d1.scalar_mul(&k, &curve);

        let params = HecIndexCalculusParams {
            fb_size: usize::MAX,
            extra_relations: 5,
            max_trials: 50_000,
            seed: 20260915,
        };
        let (found, report) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &params);
        assert_eq!(found, Some(k), "report: {report:?}");
        assert!(report.solved);
        assert!(report.smoothness_rate() > 0.0);
    }

    #[test]
    fn truncated_factor_base_still_solves_or_reports_why() {
        let curve = toy_curve(41);
        let jac = brute_force_jac_order_via_lpoly(&curve);
        let full = build_factor_base(&curve, usize::MAX);
        let l = factorise_small(jac.clone())
            .into_iter()
            .map(|(q, _)| q)
            .max()
            .unwrap();
        let d1 = subgroup_generator(&curve, &full, &jac, &l).unwrap();
        let k = BigUint::from(11u32);
        let d2 = d1.scalar_mul(&k, &curve);

        let params = HecIndexCalculusParams {
            fb_size: full.len() / 2,
            extra_relations: 5,
            max_trials: 200_000,
            seed: 7,
        };
        let (found, report) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &params);
        // A reduced base trades trials for solve size; either it got
        // there, or the report says how many smooth divisors it threw
        // away off-base.  What it must never do is answer wrongly.
        match found {
            Some(x) => assert_eq!(x, k),
            None => assert!(report.discarded_off_base > 0 || report.relations < full.len() / 2 + 1),
        }
    }
}
