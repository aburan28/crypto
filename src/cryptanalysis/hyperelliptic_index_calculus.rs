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
//! 1. **Relation search.**  Produce candidates `R = a·D_1 + b·D_2` and
//!    reduce with Cantor — either by drawing `a, b` afresh
//!    (`RelationSearch::Random`, `2⌈log₂ N⌉` group operations each) or
//!    by stepping an `r`-adding walk (`RelationSearch::Walk`, one group
//!    operation each; the default).  The reduced representative is
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
//!    of `D_1`.  Each row has at most `g + 1` non-zeros whatever `m`
//!    is, so the default solve is sparse (`LinearAlgebra::Sparse`);
//!    dense elimination ([`gaussian_eliminate_mod_n`]) is kept as the
//!    reference it has to beat.
//! 3. **Read off** `log_{D_1} D_2` and verify `D_2 = k·D_1` before
//!    returning it.  An unverified `k` is never returned.
//!
//! ## Cost, and what this module is for
//!
//! With the full degree-1 factor base and the walk, the relation stage
//! costs `≈ (m+1)·g!` group operations — `O(p)` at fixed `g` — against
//! `O(p^{g/2})` for Pollard rho on a group of size `≈ p^g`.  So the
//! asymptotic crossover is at `g = 4`, and Gaudry's variant — a
//! *reduced* factor base of size `p^{2/(g+1)}`, which this module
//! supports through `fb_size` — moves it to `g = 3`.
//!
//! The smoothness oracle runs on every candidate, so its per-call cost
//! sits directly on that stage; it is counted
//! ([`HecIndexCalculusReport::smoothness_field_ops`]) rather than
//! assumed away, since `AGENTS.md` puts an oracle's per-call work
//! inside the cost unit.  At `deg u ≤ 2` the default oracle decides
//! smoothness with one Legendre symbol instead of `O(p)` evaluations.
//!
//! What remains largest is the walk's own precomputation: the branch
//! divisors, roughly half the relation stage at `p = 251` even with
//! short step coefficients.
//!
//! The measured comparison against rho on the same instances, in the
//! unit `S = total group operations / sqrt(N)`, lives in
//! [`crate::cryptanalysis::hyperelliptic_ic_bench`] and
//! `research/notes/index-calculus/RESEARCH_HYPERELLIPTIC_IC_RHO.md`.  This module *measures*; it does
//! not claim.
//!
//! ## Scope
//!
//! - Odd characteristic, `deg f = 2g+1`, `f` squarefree — the
//!   "imaginary" model, one place at infinity.
//! - `p` small enough to enumerate: `p` must fit in a `u64`, and the
//!   factor-base build is `O(p)`.  The root finding no longer is —
//!   `SmoothnessTest::Gcd` is `O(log p)` per candidate at genus 2 —
//!   so the factor-base build is now the binding `O(p)` step.
//! - `N` = order of `D_1` must be **prime**: the solver inverts
//!   pivots mod `N`.  [`prime_order_of`] computes it for toy
//!   Jacobians and refuses a composite answer.

use std::collections::HashMap;
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::{RngCore, SeedableRng};

use crate::cryptanalysis::ec_index_calculus::{gaussian_eliminate_mod_n, sqrt_mod_p};
use crate::prime_hyperelliptic::{FpPoly, HyperellipticCurveP, MumfordDivisorP};
use crate::utils::mod_inverse;

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
    let mut ops = 0usize;
    split_into_linear_factors_counted(u, p, &mut ops)
}

/// As [`split_into_linear_factors`], adding to `ops` the `F_p`
/// multiplications it performs.
///
/// This is the smoothness oracle, and `AGENTS.md` puts "any work an
/// oracle does per call" inside the cost unit.  It is not a detail: at
/// `p = 251` the root finding costs more field multiplications per trial
/// than the walk step it follows, so a count that omitted it would
/// report a relation stage roughly half its true price — the exact shape
/// of a relabelling.
pub fn split_into_linear_factors_counted(
    u: &FpPoly,
    p: &BigUint,
    ops: &mut usize,
) -> Option<Vec<(BigUint, usize)>> {
    split_counted_with(u, p, ops, &SmoothnessTest::default())
}

/// How smoothness is decided.
///
/// The oracle runs on **every** candidate, so its per-call cost sits
/// directly on the relation stage.  Once the walk brought a candidate
/// down to one group operation, this became the dominant field cost —
/// at `p = 251` the scan spent ~600 multiplications per trial against
/// the walk's single Jacobian operation.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub enum SmoothnessTest {
    /// Evaluate `u` at every `x ∈ F_p`.  `O(p)` multiplications per
    /// call, and the reference the fast test has to match answer-for-
    /// answer.
    Scan,
    /// `u` splits into linear factors iff its radical equals
    /// `gcd(u, x^p − x)`, and `x^p mod u` costs `O(log p)`
    /// multiplications of degree-`< g` polynomials by square-and-
    /// multiply.  Roots then come from the quadratic formula at
    /// `deg ≤ 2`; higher degrees fall back to the scan (see below).
    #[default]
    Gcd,
}

/// As [`split_into_linear_factors_counted`], choosing the oracle.
pub fn split_counted_with(
    u: &FpPoly,
    p: &BigUint,
    ops: &mut usize,
    test: &SmoothnessTest,
) -> Option<Vec<(BigUint, usize)>> {
    match test {
        SmoothnessTest::Scan => split_by_scan(u, p, ops),
        SmoothnessTest::Gcd => split_by_gcd(u, p, ops),
    }
}

/// Multiply two polynomials mod `m`, charging the coefficient
/// multiplications both the product and the reduction perform.
fn poly_mul_mod(a: &FpPoly, b: &FpPoly, m: &FpPoly, ops: &mut usize) -> FpPoly {
    let (la, lb) = (a.coeffs.len().max(1), b.coeffs.len().max(1));
    *ops += la * lb;
    let prod = a.mul(b);
    let (dp, dm) = (prod.degree().unwrap_or(0), m.degree().unwrap_or(1));
    if dp >= dm {
        *ops += (dp - dm + 1) * (dm + 1);
    }
    prod.rem(m)
}

/// `x^p mod u`, by square-and-multiply on the exponent.
fn x_pow_p_mod(u: &FpPoly, p: &BigUint, ops: &mut usize) -> FpPoly {
    let x = FpPoly::x(p.clone());
    let mut result = FpPoly::one(p.clone()).rem(u);
    let mut base = x.rem(u);
    for i in 0..p.bits() {
        if p.bit(i) {
            result = poly_mul_mod(&result, &base, u, ops);
        }
        if i + 1 < p.bits() {
            base = poly_mul_mod(&base, &base, u, ops);
        }
    }
    result
}

/// Derivative of `f` over `F_p`.
fn derivative(f: &FpPoly, p: &BigUint) -> FpPoly {
    let coeffs: Vec<BigUint> = f
        .coeffs
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| (c * BigUint::from(i as u64)) % p)
        .collect();
    FpPoly::from_coeffs(coeffs, p.clone())
}

/// Roots of a **squarefree, completely split** monic `d`.
///
/// Degrees 1 and 2 close in form.  Above that it is Cantor–Zassenhaus
/// equal-degree splitting: for a random `b`, `gcd(d, (x+b)^((p−1)/2) −
/// 1)` is a non-trivial factor with probability about `1/2`, because
/// the shifted roots split into residues and non-residues.  Genus 3
/// reaches this branch — `deg u ≤ g` — and without it the oracle would
/// fall back to the `O(p)` scan exactly where the interesting
/// measurement is.
fn roots_of_split_poly(d: &FpPoly, p: &BigUint, ops: &mut usize) -> Option<Vec<BigUint>> {
    let deg = d.degree()?;
    match deg {
        0 => Some(Vec::new()),
        1 => {
            // d = x + c  ⟹  root −c (d is monic here).
            Some(vec![(p - &d.coeff(0)) % p])
        }
        2 => {
            // Monic x² + a₁x + a₀: roots (−a₁ ± sqrt(a₁² − 4a₀))/2.
            let a1 = d.coeff(1);
            let a0 = d.coeff(0);
            let disc = (&a1 * &a1 + p * p - (BigUint::from(4u32) * &a0) % p) % p;
            *ops += 2;
            // Tonelli–Shanks is a modular exponentiation: ~1.5·log₂ p
            // multiplications, the same order as the x^p step above.
            *ops += (p.bits() as usize * 3) / 2;
            let root = sqrt_mod_p(&disc, p)?;
            let two_inv = mod_inverse(&BigUint::from(2u32), p)?;
            let neg_a1 = (p - &a1) % p;
            let r1 = ((&neg_a1 + &root) % p * &two_inv) % p;
            let r2 = ((&neg_a1 + p - &root) % p * &two_inv) % p;
            *ops += 2;
            Some(vec![r1, r2])
        }
        _ => {
            // Cantor–Zassenhaus.  Deterministic seed: the caller is a
            // measurement harness, and a root finder that returns
            // different work for the same input would make the field-op
            // column unreproducible.
            let mut rng = StdRng::seed_from_u64(0x5ca1_ab1e ^ deg as u64);
            let exp = (p - BigUint::one()) >> 1;
            for _ in 0..64 {
                let b = rand_below(&mut rng, p);
                let shifted = FpPoly::from_coeffs(vec![b, BigUint::one()], p.clone());
                let powered = poly_powmod(&shifted, &exp, d, ops);
                let h = powered.sub(&FpPoly::one(p.clone()));
                if h.is_zero() {
                    continue;
                }
                let factor = d.gcd(&h).monic();
                *ops += (deg + 1) * (deg + 1);
                let fdeg = factor.degree().unwrap_or(0);
                if fdeg == 0 || fdeg == deg {
                    continue; // trivial split; try another b
                }
                let other = d.divrem(&factor).0.monic();
                let mut roots = roots_of_split_poly(&factor, p, ops)?;
                roots.extend(roots_of_split_poly(&other, p, ops)?);
                roots.sort_unstable();
                return Some(roots);
            }
            None
        }
    }
}

/// `base^e mod m`, by square-and-multiply, charging its coefficient
/// multiplications.
fn poly_powmod(base: &FpPoly, e: &BigUint, m: &FpPoly, ops: &mut usize) -> FpPoly {
    let p = &m.p;
    let mut result = FpPoly::one(p.clone()).rem(m);
    let mut acc = base.rem(m);
    for i in 0..e.bits() {
        if e.bit(i) {
            result = poly_mul_mod(&result, &acc, m, ops);
        }
        if i + 1 < e.bits() {
            acc = poly_mul_mod(&acc, &acc, m, ops);
        }
    }
    result
}

/// Smoothness by `gcd(u, x^p − x)`, with the roots read off in closed
/// form.
///
/// `u` splits into linear factors over `F_p` exactly when its radical
/// `u/gcd(u, u')` — the product of its distinct irreducible factors —
/// equals `gcd(u, x^p − x)`, which is the product of its distinct
/// *linear* factors.  Comparing degrees is enough, both being monic
/// divisors of `u` and one dividing the other.
/// `Res(a, b)` over `F_p`, by the Euclidean formula
/// `Res(a,b) = (−1)^{deg a · deg b} · lc(b)^{deg a − deg r} · Res(b, r)`
/// with `r = a mod b`.
fn resultant(a: &FpPoly, b: &FpPoly, ops: &mut usize) -> BigUint {
    let p = &a.p;
    if b.is_zero() {
        return BigUint::zero();
    }
    let (da, db) = (a.degree().unwrap_or(0), b.degree().unwrap_or(0));
    if db == 0 {
        *ops += da;
        return b.lead().modpow(&BigUint::from(da as u64), p);
    }
    let r = a.rem(b);
    *ops += (da.saturating_sub(db) + 1) * (db + 1);
    let dr = if r.is_zero() {
        0
    } else {
        r.degree().unwrap_or(0)
    };
    let mut value = resultant(b, &r, ops);
    if r.is_zero() {
        return BigUint::zero();
    }
    let factor = b.lead().modpow(&BigUint::from((da - dr) as u64), p);
    *ops += da - dr + 1;
    value = (value * factor) % p;
    if (da * db) % 2 == 1 {
        value = (p - value) % p;
    }
    value
}

/// Discriminant of a monic `u`: `(−1)^{n(n−1)/2} · Res(u, u')`.
fn discriminant(u: &FpPoly, p: &BigUint, ops: &mut usize) -> BigUint {
    let n = u.degree().unwrap_or(0);
    let du = derivative(u, p);
    if du.is_zero() {
        return BigUint::zero();
    }
    let res = resultant(u, &du, ops);
    if (n * (n - 1) / 2) % 2 == 1 {
        (p - res) % p
    } else {
        res
    }
}

/// Is `a` a square in `F_p`?  One Euler exponentiation.
fn is_square_mod(a: &BigUint, p: &BigUint, ops: &mut usize) -> bool {
    if a.is_zero() {
        return true;
    }
    *ops += (p.bits() as usize * 3) / 2;
    a.modpow(&((p - BigUint::one()) >> 1), p) == BigUint::one()
}

fn split_by_gcd(u: &FpPoly, p: &BigUint, ops: &mut usize) -> Option<Vec<(BigUint, usize)>> {
    let deg = u.degree()?;
    if deg == 0 {
        return Some(Vec::new());
    }

    // Degrees 1 and 2 — the only ones genus 2 ever produces — need no
    // `x^p` at all: a monic quadratic splits over `F_p` exactly when its
    // discriminant is a square, which is one Legendre symbol.  Going
    // through the general machinery here would cost `log₂ p` polynomial
    // squarings to learn what one exponentiation already says.
    if deg <= 2 {
        return split_low_degree(u, p, ops);
    }

    // Discriminant pre-filter.  Frobenius acts on the roots of a
    // squarefree `u` as a permutation whose cycle type is the
    // factorisation type, and `disc(u)` is a square exactly when that
    // permutation is even.  Splitting completely is the identity
    // permutation, which is even — so a **non-square discriminant
    // proves `u` does not split**, for one resultant on a degree-`≤ g`
    // polynomial and one Euler exponentiation.
    //
    // It rejects exactly half of all candidates: the odd types are
    // `(2,1)` at degree 3 (density 1/2) and `(2,1,1) + (4)` at degree 4
    // (1/4 + 1/4).  That is half the candidates never reaching the
    // `x^p mod u` that dominates this oracle.
    //
    // `disc = 0` means repeated roots, which says nothing either way —
    // the radical path below handles it.
    let disc = discriminant(u, p, ops);
    if !disc.is_zero() && !is_square_mod(&disc, p, ops) {
        return None;
    }

    let du = derivative(u, p);
    let radical = if du.is_zero() {
        // Only possible when p | every exponent; fall back rather than
        // reason about it, since it cannot occur for deg u ≤ g < p.
        return split_by_scan(u, p, ops);
    } else {
        let g = u.gcd(&du);
        *ops += (deg + 1) * (deg + 1);
        u.divrem(&g).0.monic()
    };

    let xp = x_pow_p_mod(u, p, ops);
    let xp_minus_x = xp.sub(&FpPoly::x(p.clone()));
    let linear_part = u.gcd(&xp_minus_x).monic();
    *ops += (deg + 1) * (deg + 1);

    if linear_part.degree() != radical.degree() {
        // Some irreducible factor has degree ≥ 2: not smooth.
        return None;
    }

    let roots = match roots_of_split_poly(&linear_part, p, ops) {
        Some(r) => r,
        // Degree ≥ 3 radical: no closed form here, so use the scan
        // rather than return a wrong answer.
        None => return split_by_scan(u, p, ops),
    };

    // Multiplicities, by dividing each root out of `u` while it goes.
    let mut out = Vec::with_capacity(roots.len());
    let mut rest = u.clone();
    let mut total = 0usize;
    for x in roots {
        let lin = FpPoly::from_coeffs(vec![(p - &x) % p, BigUint::one()], p.clone());
        let mut mult = 0usize;
        loop {
            *ops += rest.degree().unwrap_or(0);
            let (q, r) = rest.divrem(&lin);
            if !r.is_zero() {
                break;
            }
            rest = q;
            mult += 1;
        }
        if mult > 0 {
            total += mult;
            out.push((x, mult));
        }
    }
    // The degree check above says this holds; assert it rather than
    // trust it, since a wrong "smooth" would poison the relation matrix
    // with a row that is not a relation.
    if total != deg {
        return None;
    }
    out.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    Some(out)
}

/// Closed-form split for `deg u ∈ {1, 2}`.
///
/// `deg 1` always splits.  For monic `x² + a₁x + a₀` the discriminant
/// `a₁² − 4a₀` decides it: a non-residue means an irreducible quadratic
/// (not smooth), zero means a double root, and a square gives the two
/// roots by the usual formula.
fn split_low_degree(u: &FpPoly, p: &BigUint, ops: &mut usize) -> Option<Vec<(BigUint, usize)>> {
    match u.degree()? {
        1 => {
            // Monic x + c has the single root −c.
            let inv_lead = mod_inverse(&u.coeff(1), p)?;
            *ops += 1;
            let c = (&u.coeff(0) * &inv_lead) % p;
            Some(vec![((p - c) % p, 1)])
        }
        2 => {
            let lead_inv = mod_inverse(&u.coeff(2), p)?;
            *ops += 2;
            let a1 = (&u.coeff(1) * &lead_inv) % p;
            let a0 = (&u.coeff(0) * &lead_inv) % p;
            let four_a0 = (BigUint::from(4u32) * &a0) % p;
            let disc = (&a1 * &a1 % p + p - four_a0) % p;
            *ops += 2;
            let neg_a1 = (p - &a1) % p;
            let two_inv = mod_inverse(&BigUint::from(2u32), p)?;

            if disc.is_zero() {
                // Double root −a₁/2.
                *ops += 1;
                let r = (&neg_a1 * &two_inv) % p;
                return Some(vec![(r, 2)]);
            }
            // One Euler/Tonelli exponentiation: ~1.5·log₂ p mul-mods.
            *ops += (p.bits() as usize * 3) / 2;
            let root = sqrt_mod_p(&disc, p)?; // None ⟹ irreducible ⟹ not smooth
            let r1 = ((&neg_a1 + &root) % p * &two_inv) % p;
            let r2 = ((&neg_a1 + p - &root) % p * &two_inv) % p;
            *ops += 2;
            let (lo, hi) = if r1 <= r2 { (r1, r2) } else { (r2, r1) };
            Some(vec![(lo, 1), (hi, 1)])
        }
        _ => None,
    }
}

fn split_by_scan(u: &FpPoly, p: &BigUint, ops: &mut usize) -> Option<Vec<(BigUint, usize)>> {
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
        // Horner: one multiplication per coefficient.
        *ops += rest.degree().unwrap_or(0);
        if !rest.eval(&x).is_zero() {
            continue;
        }
        // Divide out (x − xi) as often as it goes.
        let lin = FpPoly::from_coeffs(vec![(p - &x) % p, BigUint::one()], p.clone());
        let mut mult = 0usize;
        loop {
            // Dividing by a monic linear factor is synthetic division:
            // one multiplication per remaining coefficient.
            *ops += rest.degree().unwrap_or(0);
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
    let mut ops = 0usize;
    decompose_over_factor_base_counted(curve, d, fb, &mut ops)
}

/// As [`decompose_over_factor_base`], adding its `F_p` multiplications
/// to `ops`.
pub fn decompose_over_factor_base_counted(
    curve: &HyperellipticCurveP,
    d: &MumfordDivisorP,
    fb: &HecFactorBase,
    ops: &mut usize,
) -> Option<Vec<(usize, i64)>> {
    decompose_counted_with(curve, d, fb, ops, &SmoothnessTest::default())
}

/// Why a candidate divisor did or did not yield a relation.
///
/// The three cases used to be squeezed into an `Option`, which meant the
/// caller re-ran the whole oracle on every failure just to learn which
/// kind of failure it was.  At genus 4 more than nine candidates in ten
/// fail, so that probe was about half of all the oracle's work —
/// counted, and wasted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Decomposition {
    /// Splits over degree-1 places, all of them in the factor base.
    Smooth(Vec<(usize, i64)>),
    /// Splits over degree-1 places, all but **one** of them in the
    /// factor base — a partial relation for the large-prime variation.
    /// `coef` is the off-base place's coefficient against its canonical
    /// representative `(x, y ≤ (p−1)/2)`, as for factor-base entries.
    OneLargePrime {
        entries: Vec<(usize, i64)>,
        x: BigUint,
        coef: i64,
    },
    /// Splits over degree-1 places, but meets two or more places the
    /// (truncated) factor base does not hold.
    SmoothOffBase,
    /// Does not split into degree-1 places.
    NotSmooth,
}

/// As [`decompose_over_factor_base_counted`], choosing the oracle and
/// reporting which of the three outcomes occurred.
pub fn classify_counted_with(
    curve: &HyperellipticCurveP,
    d: &MumfordDivisorP,
    fb: &HecFactorBase,
    ops: &mut usize,
    test: &SmoothnessTest,
) -> Decomposition {
    let p = &curve.p;
    let half = (p - BigUint::one()) >> 1;
    let roots = match split_counted_with(&d.u, p, ops, test) {
        Some(r) => r,
        None => return Decomposition::NotSmooth,
    };

    let mut acc: HashMap<usize, i64> = HashMap::new();
    // At most one place outside a truncated base is kept; `u` has
    // distinct roots per place (a reduced divisor never holds both `P`
    // and `−P`), so one `x` is one place.
    let mut off: Option<(BigUint, i64)> = None;
    for (x, mult) in roots {
        *ops += d.v.degree().unwrap_or(0);
        let y = d.v.eval(&x);
        debug_assert!(curve.is_on_curve(&x, &y));
        let sign: i64 = if y.is_zero() || y <= half { 1 } else { -1 };
        let idx = match fb.index_of_x(&x) {
            Some(i) => i,
            // Smooth, but this place is outside a truncated base.
            None => {
                if off.is_some() {
                    return Decomposition::SmoothOffBase;
                }
                off = Some((x, sign * mult as i64));
                continue;
            }
        };
        debug_assert_eq!(fb.entries[idx].y, if sign > 0 { y.clone() } else { p - &y });
        *acc.entry(idx).or_insert(0) += sign * mult as i64;
    }

    let mut out: Vec<(usize, i64)> = acc.into_iter().filter(|&(_, c)| c != 0).collect();
    out.sort_unstable_by_key(|&(i, _)| i);
    match off {
        None => Decomposition::Smooth(out),
        Some((x, coef)) => Decomposition::OneLargePrime {
            entries: out,
            x,
            coef,
        },
    }
}

/// As [`decompose_over_factor_base_counted`], choosing the oracle.
pub fn decompose_counted_with(
    curve: &HyperellipticCurveP,
    d: &MumfordDivisorP,
    fb: &HecFactorBase,
    ops: &mut usize,
    test: &SmoothnessTest,
) -> Option<Vec<(usize, i64)>> {
    let p = &curve.p;
    let half = (p - BigUint::one()) >> 1;
    let roots = split_counted_with(&d.u, p, ops, test)?;

    let mut acc: HashMap<usize, i64> = HashMap::new();
    for (x, mult) in roots {
        *ops += d.v.degree().unwrap_or(0);
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

impl HecIndexCalculusReport {
    /// Group operations per trial, excluding precomputation — the
    /// quantity the relation-search mode actually changes.  `Random`
    /// sits at `2⌈log₂ N⌉`; `Walk` sits at 1.
    pub fn ops_per_trial(&self) -> f64 {
        if self.trials == 0 {
            0.0
        } else {
            (self.jacobian_ops - self.precompute_ops) as f64 / self.trials as f64
        }
    }
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
    /// Smooth, but met a place outside a truncated factor base and was
    /// not kept — two or more such places, or large primes off.
    pub discarded_off_base: usize,
    /// Partial relations (one off-base place) kept for combination.
    pub partial_relations: usize,
    /// Full relations made by combining two partials.
    pub combined_relations: usize,
    pub relations: usize,
    pub jacobian_ops: usize,
    /// Of those, the walk's step precomputation — a fixed cost that
    /// amortises over trials and so dominates at toy `N`, exactly as
    /// rho's branch precomputation does.  Split out for the same
    /// reason: the per-trial price is the thing the walk changed.
    pub precompute_ops: usize,
    /// Wall-clock nanoseconds spent on the relation search's group
    /// operations, excluding the oracle.
    ///
    /// Divided by `jacobian_ops` this gives the cost of one group
    /// operation **in situ** — in the run itself, with its allocation
    /// and cache behaviour — which is what the oracle and the solve must
    /// be converted by.  A tight calibration loop under-measures it: it
    /// repeats one addition on cache-resident operands, so converting a
    /// measured oracle time by it overstates the oracle instead.  Both
    /// errors were live in this file at different times, in opposite
    /// directions, and only wall clock caught either.
    pub walk_wall_ns: u64,
    /// Wall-clock nanoseconds spent in the linear algebra.
    pub solve_wall_ns: u64,
    /// Wall-clock nanoseconds spent storing and combining large-prime
    /// partials.  Bookkeeping multiplies little mod `N` but is not free,
    /// and uncharged bookkeeping is how round six under-priced the solve
    /// by an order of magnitude, so it is timed and converted like the
    /// oracle.
    pub partial_wall_ns: u64,
    /// Wall-clock nanoseconds spent inside the smoothness oracle.
    ///
    /// The `smoothness_field_ops` count below is a hand-derived charge —
    /// coefficient multiplications — and it turned out to understate the
    /// oracle's real cost by roughly an order of magnitude: cutting that
    /// count by 60% cut wall clock by 1.5x, where a faithful charge
    /// would have predicted a few percent.  The charge misses what the
    /// implementation actually does per multiplication (allocation, the
    /// division loop inside `rem`, clones).  This field is measured
    /// instead, so a caller can convert the oracle at the same measured
    /// rate it converts everything else.
    pub smoothness_wall_ns: u64,
    /// `F_p` multiplications spent in the smoothness oracle (root
    /// finding and the decomposition's evaluations).  Reported in field
    /// operations, not group operations — the caller converts, because
    /// only a measurement can say what the ratio is on a given machine.
    pub smoothness_field_ops: usize,
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

/// How relations are drawn.
///
/// The cost difference is the whole point: `Random` pays a fresh pair of
/// scalar multiplications, `2⌈log₂ N⌉` group operations, for every
/// candidate divisor; `Walk` pays **one** group operation per candidate
/// by stepping an existing `R = a·D₁ + b·D₂` to `R + S_j` and adding the
/// step's known `(a_j, b_j)` to its coefficients — the same trick that
/// makes an `r`-adding rho walk cheap.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RelationSearch {
    /// Independent `(a, b)` per trial.  Simple, and the reference the
    /// walk has to beat.
    Random,
    /// `r`-adding walk, perturbed every `restart_interval` trials.
    ///
    /// The perturbation is not optional bookkeeping.  A deterministic
    /// walk enters a cycle after `O(sqrt(N))` steps, and this attack
    /// needs `m + 1 ≈ p/2` relations out of a group of size `N ≈ p²` —
    /// the same order — so an undisturbed walk would start re-deriving
    /// relations it already has, exactly when it still needs new ones.
    /// It escapes by adding a **randomly chosen** step rather than by
    /// building a fresh point: the choice comes from the RNG instead of
    /// the position hash, so the trajectory leaves its orbit for one
    /// group operation instead of `2⌈log₂ N⌉`.
    Walk {
        branches: usize,
        restart_interval: usize,
        /// Bit length of the step coefficients `a_j, b_j`.
        ///
        /// Each branch costs `2·step_bits` group operations to build,
        /// and precomputation is the relation stage's largest term once
        /// the per-trial cost is one operation.  Short coefficients do
        /// not weaken anything: a relation is an algebraic identity
        /// verified as it is recorded, so the distribution of the steps
        /// affects only how fast smooth divisors turn up, which the
        /// trial count measures directly.
        step_bits: u32,
    },
    /// An adding walk whose steps are **factor-base places**, not
    /// combinations of `D₁` and `D₂`.
    ///
    /// The walk carries `R = a·D₁ + b·D₂ + Σ n_j F_j` and counts the
    /// steps it took.  When `R` decomposes as `Σ c_i F_i`, the relation
    /// is `Σ (c_i − n_i) y_i − b·k ≡ a` — the same shape as before,
    /// because `log F_j` was already one of the unknowns.  Not knowing
    /// the steps' discrete logarithms costs nothing here, and buys the
    /// whole precomputation: the factor base is already built, so the
    /// steps are free.
    ///
    /// **Pollard rho cannot do this.**  Its steps must have known
    /// `(a_j, b_j)` or a collision says nothing, so it has to pay for
    /// every branch divisor it uses.  This is a structural asymmetry
    /// between the two methods, not an optimisation withheld from the
    /// reference — which is why it is worth the extra code.
    ///
    /// Cost: rows gain up to one entry per distinct step used, so they
    /// carry `≤ g + branches + 1` non-zeros instead of `≤ g + 1`.  The
    /// sparse solve absorbs that; the measurement says by how much.
    ///
    /// **`D₁` and `D₂` must stay in the step set.**  A step by a place
    /// changes neither `a` nor `b`, so a walk that only ever steps by
    /// places produces rows that all share one `(a, b)` — every row then
    /// reads `a + b·k = (something in the y's)`, subtracting any two
    /// eliminates `k`, and the system pins the `y`'s while leaving `k`
    /// free.  The rows are true identities and still say nothing about
    /// the logarithm.  `d_step_in` of the steps are `D₁` or `D₂`, which
    /// cost nothing to build either, and move `a` and `b`.
    ///
    /// Steps are drawn from the RNG rather than from a hash of the
    /// position: relation collection has no collision to detect, so
    /// there is nothing determinism buys here, and a random sequence
    /// cannot cycle — which removes the restart machinery the
    /// `(a_j, b_j)` walk needs.
    FactorBaseWalk {
        branches: usize,
        /// One step in `d_step_in` is by `D₁` or `D₂` rather than a
        /// place.  4 means a quarter of steps move `(a, b)`.
        d_step_in: u64,
    },
}

impl RelationSearch {
    /// The factor-base walk this module uses unless told otherwise.
    ///
    /// 64 branches rather than 16: they cost nothing to build, and more
    /// of them makes the walk closer to a random map.
    pub fn factor_base_walk() -> Self {
        Self::FactorBaseWalk {
            branches: 64,
            d_step_in: 4,
        }
    }

    /// The `(a_j, b_j)`-step walk, kept as the reference the
    /// factor-base walk has to beat.
    pub fn walk() -> Self {
        // 16 branches, matching the rho reference: Teske's analysis says
        // an r-adding walk is within a few percent of the random-map
        // ideal from r = 16, and every extra branch is a scalar-mult
        // pair of precomputation that has to amortise over the trials.
        Self::Walk {
            branches: 16,
            restart_interval: 64,
            step_bits: 8,
        }
    }
}

/// How the relation system is solved.
///
/// The matrix is extremely sparse by construction: a genus-`g` relation
/// touches at most `g` factor-base columns, plus the `k` column, so each
/// row carries `≤ g + 1` non-zeros regardless of `m`.  A dense `O(m³)`
/// elimination ignores that and, once the relation stage is walked
/// rather than re-drawn, becomes the dominant cost.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LinearAlgebra {
    /// Dense Gaussian elimination mod `N`.  Simple, and the reference
    /// the sparse solve has to beat.
    Dense,
    /// Sparse elimination mod `N` with Markowitz-style pivoting: at each
    /// step take the factor-base column with fewest remaining rows, and
    /// within it the shortest row.  Only `k` is read off, so columns
    /// that never become pivots are simply never formed.
    Sparse,
}

/// Parameters of the relation search and solve.
#[derive(Clone, Debug)]
pub struct HecIndexCalculusParams {
    /// Factor base cap; `usize::MAX` for every degree-1 place.
    pub fb_size: usize,
    /// Rows collected beyond `fb_size + 1` unknowns.
    pub extra_relations: usize,
    /// Give up after this many candidate divisors in total.
    pub max_trials: usize,
    pub seed: u64,
    /// How candidates are produced; see [`RelationSearch`].
    pub search: RelationSearch,
    /// How the relation system is solved; see [`LinearAlgebra`].
    pub linear_algebra: LinearAlgebra,
    /// How smoothness is decided; see [`SmoothnessTest`].
    pub smoothness: SmoothnessTest,
    /// Single large primes (Thériault): with a truncated factor base,
    /// keep a candidate that has exactly one place outside it, and
    /// combine two such partials that share the place into one full
    /// relation.  Off by default; it does nothing with the full base,
    /// which has no outside places.
    pub large_primes: bool,
}

impl Default for HecIndexCalculusParams {
    fn default() -> Self {
        Self {
            fb_size: usize::MAX,
            extra_relations: 8,
            max_trials: 200_000,
            seed: 0,
            search: RelationSearch::factor_base_walk(),
            linear_algebra: LinearAlgebra::Sparse,
            smoothness: SmoothnessTest::Gcd,
            large_primes: false,
        }
    }
}

/// `R = a·D₁ + b·D₂ + Σ n_j F_j` decomposed as `Σ c_i F_i` gives the row
/// `Σ (c_i − n_i) y_i − b·k ≡ a`, so the steps the walk took are
/// subtracted from the decomposition it found.
fn subtract_steps(
    decomposition: Vec<(usize, i64)>,
    taken: &HashMap<usize, i64>,
) -> Vec<(usize, i64)> {
    if taken.is_empty() {
        return decomposition;
    }
    let mut acc: HashMap<usize, i64> = decomposition.into_iter().collect();
    for (&idx, &count) in taken {
        *acc.entry(idx).or_insert(0) -= count;
    }
    let mut out: Vec<(usize, i64)> = acc.into_iter().filter(|&(_, c)| c != 0).collect();
    out.sort_unstable_by_key(|&(i, _)| i);
    out
}

/// Branch selector for the adding walk: a hash of the Mumford
/// representation, so equal classes take equal steps.
fn walk_branch(d: &MumfordDivisorP, branches: usize) -> usize {
    let mut acc: u64 = 0;
    for poly in [&d.u, &d.v] {
        for limb in poly.coeffs.iter().flat_map(|c| c.to_u64_digits()) {
            acc = acc.wrapping_mul(0x9e37_79b9_7f4a_7c15).wrapping_add(limb);
        }
        acc = acc.wrapping_mul(0x9e37_79b9_7f4a_7c15).wrapping_add(1);
    }
    (acc >> 32) as usize % branches.max(1)
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

    // One scalar multiplication pair, in group operations.  Both modes
    // are charged this whenever they build a divisor from scratch.
    let scalar_pair_ops = 2 * (n.bits() as usize + 1);

    // Precomputed steps S_j = a_j·D₁ + b_j·D₂ for the walk.  Charged up
    // front, like rho's branches: a precomputation left out of the
    // accounting is a cost moved, not removed.
    let (branches, restart_interval, step_bits) = match params.search {
        RelationSearch::Random => (0usize, usize::MAX, 0),
        RelationSearch::Walk {
            branches,
            restart_interval,
            step_bits,
        } => (branches.max(1), restart_interval.max(1), step_bits.max(1)),
        RelationSearch::FactorBaseWalk { branches, .. } => {
            (branches.min(fb.len()).max(1), usize::MAX, 0)
        }
    };
    let walking = !matches!(params.search, RelationSearch::Random);

    // Steps with `step_bits`-bit coefficients cost `2·step_bits`
    // operations each instead of `2⌈log₂ N⌉`.  Factor-base steps cost
    // nothing at all: they are places the factor base already holds.
    let step_bound = BigUint::one() << step_bits.max(1);
    let step_build_ops = 2 * step_bits as usize;
    let mut steps: Vec<(BigUint, BigUint, MumfordDivisorP)> = Vec::new();
    let mut fb_steps: Vec<usize> = Vec::new();
    match params.search {
        RelationSearch::FactorBaseWalk { .. } => {
            // Distinct places, so a step is never silently doubled.
            let mut chosen: Vec<usize> = (0..fb.len()).collect();
            for i in (1..chosen.len()).rev() {
                let j = (rng.next_u64() as usize) % (i + 1);
                chosen.swap(i, j);
            }
            fb_steps = chosen.into_iter().take(branches).collect();
        }
        _ => {
            for _ in 0..branches {
                let a = rand_below(&mut rng, &step_bound) + BigUint::one();
                let b = rand_below(&mut rng, &step_bound) + BigUint::one();
                let s = d1
                    .scalar_mul(&a, curve)
                    .add(&d2.scalar_mul(&b, curve), curve);
                report.jacobian_ops += step_build_ops;
                report.precompute_ops += step_build_ops;
                steps.push((a, b, s));
            }
        }
    }

    /// Walk position: coefficients of `D₁` and `D₂`, the divisor, and
    /// how many times each factor-base step has been taken.
    struct Position {
        a: BigUint,
        b: BigUint,
        r: MumfordDivisorP,
        taken: HashMap<usize, i64>,
    }

    let mut current: Option<Position> = None;
    let mut since_restart = 0usize;
    // Large-prime partials, keyed by the off-base place's `x`: the first
    // partial seen at each place, which every later one combines with.
    let mut partials: HashMap<Vec<u32>, (HecRelation, i64)> = HashMap::new();

    while relations.len() < wanted && report.trials < params.max_trials {
        report.trials += 1;

        let step_started = Instant::now();
        let mut pos = match current.take() {
            Some(pos) if walking => {
                // One step of the walk: one group operation.
                match params.search {
                    RelationSearch::FactorBaseWalk { d_step_in, .. } => {
                        // A quarter of the steps move `(a, b)`; the rest
                        // move through the factor base.  Both cost one
                        // group operation and nothing to prepare.
                        let roll = rng.next_u64() % d_step_in.max(2);
                        let mut taken = pos.taken;
                        let (mut a, mut b) = (pos.a, pos.b);
                        let next = if roll == 0 {
                            a = (a + BigUint::one()) % n;
                            report.jacobian_ops += 1;
                            pos.r.add(d1, curve)
                        } else if roll == 1 {
                            b = (b + BigUint::one()) % n;
                            report.jacobian_ops += 1;
                            pos.r.add(d2, curve)
                        } else {
                            let idx = fb_steps[(rng.next_u64() as usize) % fb_steps.len()];
                            *taken.entry(idx).or_insert(0) += 1;
                            report.jacobian_ops += 1;
                            pos.r.add(&fb.entries[idx].divisor, curve)
                        };
                        Position {
                            a,
                            b,
                            r: next,
                            taken,
                        }
                    }
                    _ => {
                        let (aj, bj, sj) = &steps[walk_branch(&pos.r, steps.len())];
                        let next = pos.r.add(sj, curve);
                        report.jacobian_ops += 1;
                        Position {
                            a: (&pos.a + aj) % n,
                            b: (&pos.b + bj) % n,
                            r: next,
                            taken: pos.taken,
                        }
                    }
                }
            }
            // Fresh start.  The factor-base walk needs no random
            // starting point — its own steps supply the randomness — so
            // it starts at `D₁ + D₂` for one operation instead of two
            // scalar multiplications.
            _ if matches!(params.search, RelationSearch::FactorBaseWalk { .. }) => {
                report.jacobian_ops += 1;
                report.precompute_ops += 1;
                Position {
                    a: BigUint::one(),
                    b: BigUint::one(),
                    r: d1.add(d2, curve),
                    taken: HashMap::new(),
                }
            }
            _ => {
                let a = rand_below(&mut rng, n);
                let b = rand_below(&mut rng, n);
                let r = d1
                    .scalar_mul(&a, curve)
                    .add(&d2.scalar_mul(&b, curve), curve);
                report.jacobian_ops += scalar_pair_ops;
                if walking {
                    // A start is precomputation: it buys position, not a
                    // candidate the cheap step could not reach.
                    report.precompute_ops += scalar_pair_ops;
                }
                since_restart = 0;
                Position {
                    a,
                    b,
                    r,
                    taken: HashMap::new(),
                }
            }
        };

        if walking && !matches!(params.search, RelationSearch::FactorBaseWalk { .. }) {
            since_restart += 1;
            if since_restart >= restart_interval {
                // Escape the orbit with one randomly chosen step: the
                // index comes from the RNG rather than the position, so
                // the next point is off the deterministic trajectory.
                match params.search {
                    RelationSearch::FactorBaseWalk { .. } => {
                        let idx = fb_steps[(rng.next_u64() as usize) % fb_steps.len()];
                        pos.r = pos.r.add(&fb.entries[idx].divisor, curve);
                        *pos.taken.entry(idx).or_insert(0) += 1;
                    }
                    _ => {
                        let j = (rng.next_u64() as usize) % steps.len();
                        let (aj, bj, sj) = &steps[j];
                        pos.r = pos.r.add(sj, curve);
                        pos.a = (&pos.a + aj) % n;
                        pos.b = (&pos.b + bj) % n;
                    }
                }
                report.jacobian_ops += 1;
                report.precompute_ops += 1;
                since_restart = 0;
            }
        }

        report.walk_wall_ns += step_started.elapsed().as_nanos() as u64;

        let (a, b, r) = (pos.a.clone(), pos.b.clone(), pos.r.clone());
        let taken = pos.taken.clone();
        if walking {
            current = Some(pos);
        }

        if r.is_identity() {
            // a·D₁ + b·D₂ + Σ n_j F_j = 0 is a relation whose only
            // factor-base part is the steps the walk took.
            relations.push(HecRelation {
                coef_a: a,
                coef_b: b,
                entries: subtract_steps(Vec::new(), &taken),
            });
            report.smooth_trials += 1;
            continue;
        }
        let mut field_ops = 0usize;
        // One oracle call per candidate.  The three outcomes come back
        // from that call rather than from re-running it on failure.
        let oracle_started = Instant::now();
        let classified = classify_counted_with(curve, &r, fb, &mut field_ops, &params.smoothness);
        report.smoothness_wall_ns += oracle_started.elapsed().as_nanos() as u64;
        report.smoothness_field_ops += field_ops;
        match classified {
            Decomposition::Smooth(entries) => {
                report.smooth_trials += 1;
                relations.push(HecRelation {
                    coef_a: a,
                    coef_b: b,
                    entries: subtract_steps(entries, &taken),
                });
            }
            Decomposition::OneLargePrime { entries, x, coef } if params.large_primes => {
                report.smooth_trials += 1;
                let started = Instant::now();
                let rel = HecRelation {
                    coef_a: a,
                    coef_b: b,
                    entries: subtract_steps(entries, &taken),
                };
                match partials.get(&x.to_u32_digits()) {
                    None => {
                        report.partial_relations += 1;
                        partials.insert(x.to_u32_digits(), (rel, coef));
                    }
                    Some((first, first_coef)) => {
                        report.partial_relations += 1;
                        if let Some(full) = combine_partials(first, *first_coef, &rel, coef, n) {
                            report.combined_relations += 1;
                            relations.push(full);
                        }
                    }
                }
                report.partial_wall_ns += started.elapsed().as_nanos() as u64;
            }
            Decomposition::OneLargePrime { .. } | Decomposition::SmoothOffBase => {
                report.smooth_trials += 1;
                report.discarded_off_base += 1;
            }
            Decomposition::NotSmooth => {}
        }
    }
    report.relations = relations.len();
    relations
}

/// `e₂·r₁ − e₁·r₂` for two partials `r₁ = … + e₁·P` and `r₂ = … + e₂·P`
/// sharing the off-base place `P`: the combination eliminates `log P`
/// and leaves a relation over the factor base alone.
///
/// Returns `None` when the combination is the zero relation — the same
/// partial reached twice — which says nothing and would only cost a row.
fn combine_partials(
    r1: &HecRelation,
    e1: i64,
    r2: &HecRelation,
    e2: i64,
    n: &BigUint,
) -> Option<HecRelation> {
    let scale = |v: &BigUint, k: i64| -> BigUint {
        let t = (v * BigUint::from(k.unsigned_abs())) % n;
        if k < 0 {
            (n - t) % n
        } else {
            t
        }
    };
    let coef_a = (scale(&r1.coef_a, e2) + scale(&r2.coef_a, -e1)) % n;
    let coef_b = (scale(&r1.coef_b, e2) + scale(&r2.coef_b, -e1)) % n;
    let mut acc: HashMap<usize, i64> = HashMap::new();
    for &(j, c) in &r1.entries {
        *acc.entry(j).or_insert(0) += e2 * c;
    }
    for &(j, c) in &r2.entries {
        *acc.entry(j).or_insert(0) -= e1 * c;
    }
    let mut entries: Vec<(usize, i64)> = acc.into_iter().filter(|&(_, c)| c != 0).collect();
    entries.sort_unstable_by_key(|&(j, _)| j);
    if entries.is_empty() && coef_a.is_zero() && coef_b.is_zero() {
        return None;
    }
    Some(HecRelation {
        coef_a,
        coef_b,
        entries,
    })
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
    // Rows from a walk are correlated — consecutive relations differ by
    // a step or two — so a fixed margin of spare rows is sometimes not
    // enough and the system pins the `y`'s while leaving `k` free.  That
    // showed up as one seed in ten failing at `p = 41`.  Collect more
    // rows and solve again rather than hope the margin was right; the
    // extra rows are counted like any others.
    let mut extra = params.extra_relations;
    let mut k = None;
    for round in 0..4 {
        let wanted = m + 1 + extra;
        let collected = collect_relations(curve, d1, d2, n, &fb, wanted, params, &mut report);
        if collected.len() < m + 1 {
            return (None, report);
        }
        let solve_started = Instant::now();
        let solved = solve_for_logarithm(curve, &collected, m, n, params, &mut report);
        report.solve_wall_ns += solve_started.elapsed().as_nanos() as u64;
        if let Some(candidate) = solved {
            if &d1.scalar_mul(&candidate, curve) == d2 {
                k = Some(candidate);
                break;
            }
        }
        // Grow the margin geometrically; `round` bounds the total.
        extra += m / 4 + 8;
        let _ = round;
    }
    let k = match k {
        Some(k) => k,
        None => return (None, report),
    };

    if &d1.scalar_mul(&k, curve) == d2 {
        report.solved = true;
        (Some(k), report)
    } else {
        (None, report)
    }
}

/// Eliminate every factor-base unknown from the sparse system and read
/// off `k` (column index `m`).
///
/// Returns `(k, mul_mods)` — the multiplication count is **measured**,
/// not estimated, because the whole point of the sparse path is that its
/// cost no longer follows `rows·m²`, so an estimate would not track it.
/// `None` when the system does not determine `k`.
fn sparse_solve_for_k(
    mut rows: Vec<Vec<(usize, BigUint)>>,
    mut rhs: Vec<BigUint>,
    m: usize,
    n: &BigUint,
) -> Option<(BigUint, usize)> {
    let mut ops = 0usize;
    // Rows stay sorted by column throughout: `solve_for_logarithm` builds
    // them sorted and the merge below preserves the order, so membership
    // is a binary search rather than a scan.
    //
    // `col_count[j]` is exact — the number of live rows carrying column
    // `j` — so choosing a pivot column is a scan of `m` integers.  The
    // previous version recomputed every column's live rows on every
    // pivot, re-testing membership by a linear search of each row; that
    // work multiplies nothing mod `N`, so it never appeared in the counted
    // mul-mods, and it was most of the solve's wall time.
    //
    // `holders[j]` may still hold stale rows (eliminated, or cancelled out
    // of `j`), filtered when `j` is pivoted on, but a row is only pushed
    // when it *gains* `j`, so the lists no longer fill with duplicates.
    let mut col_count = vec![0usize; m + 1];
    let mut holders: Vec<Vec<usize>> = vec![Vec::new(); m + 1];
    for (i, row) in rows.iter().enumerate() {
        for &(j, _) in row {
            col_count[j] += 1;
            holders[j].push(i);
        }
    }
    let mut eliminated = vec![false; rows.len()];
    let has = |row: &[(usize, BigUint)], col: usize| row.binary_search_by_key(&col, |&(j, _)| j);

    // Markowitz-lite: the factor-base column with the fewest live rows,
    // and within it the shortest row — this is what keeps fill-in from
    // turning the sparse solve back into a dense one.  Ties go to the
    // lowest column, so the elimination order — and with it the mul-mod
    // count — is deterministic; the `HashMap` iteration it replaces was
    // not.  The loop ends when nothing is left but the k column.
    while let Some(col) = (0..m)
        .filter(|&j| col_count[j] > 0)
        .min_by_key(|&j| col_count[j])
    {
        let mut live = std::mem::take(&mut holders[col]);
        live.retain(|&i| !eliminated[i] && has(&rows[i], col).is_ok());
        // A row that cancels out of `col` and later fills back into it is
        // pushed again while its stale entry is still here, so the list
        // can name a row twice.  Visiting it twice would find `col`
        // already gone on the second visit.
        live.sort_unstable();
        live.dedup();
        debug_assert_eq!(live.len(), col_count[col]);
        let pivot = *live.iter().min_by_key(|&&i| rows[i].len())?;

        // Normalise the pivot row.
        let at = has(&rows[pivot], col).ok()?;
        let inv = mod_inverse(&rows[pivot][at].1, n)?;
        for (_, v) in rows[pivot].iter_mut() {
            *v = (&*v * &inv) % n;
            ops += 1;
        }
        rhs[pivot] = (&rhs[pivot] * &inv) % n;
        ops += 1;
        eliminated[pivot] = true;
        let pivot_row = std::mem::take(&mut rows[pivot]);
        let pivot_rhs = rhs[pivot].clone();
        for &(j, _) in &pivot_row {
            col_count[j] -= 1;
        }

        // Eliminate `col` from every other live row carrying it:
        // row_i ← row_i − factor·pivot_row, as one merge of two sorted
        // rows.
        for i in live {
            if i == pivot {
                continue;
            }
            let old = std::mem::take(&mut rows[i]);
            let factor = old[has(&old, col).ok()?].1.clone();
            let mut merged: Vec<(usize, BigUint)> = Vec::with_capacity(old.len() + pivot_row.len());
            let (mut a, mut b) = (0, 0);
            while a < old.len() || b < pivot_row.len() {
                let ja = old.get(a).map_or(usize::MAX, |e| e.0);
                let jb = pivot_row.get(b).map_or(usize::MAX, |e| e.0);
                if ja < jb {
                    merged.push(old[a].clone());
                    a += 1;
                    continue;
                }
                let term = (&factor * &pivot_row[b].1) % n;
                ops += 1;
                let v = if ja == jb {
                    let v = (&old[a].1 + n - &term) % n;
                    a += 1;
                    v
                } else {
                    // Fill-in: a column the row did not carry before.
                    (n - &term) % n
                };
                if !v.is_zero() {
                    if ja != jb {
                        holders[jb].push(i);
                    }
                    merged.push((jb, v));
                }
                b += 1;
            }
            for &(j, _) in &old {
                col_count[j] -= 1;
            }
            for &(j, _) in &merged {
                col_count[j] += 1;
            }
            rhs[i] = (&rhs[i] + n - &((&factor * &pivot_rhs) % n)) % n;
            ops += 1;
            rows[i] = merged;
        }
    }

    // A surviving row is `coef·k ≡ rhs`; any of them gives `k`.
    for (i, row) in rows.iter().enumerate() {
        if eliminated[i] {
            continue;
        }
        if row.len() == 1 && row[0].0 == m {
            let inv = mod_inverse(&row[0].1, n)?;
            ops += 1;
            return Some(((&rhs[i] * &inv) % n, ops));
        }
    }
    None
}

/// Build the relation system and solve it for `k = log_{D₁} D₂`.
///
/// Returns `None` when the rows do not determine `k` — which is a
/// statement about the rows, not about the instance, so the caller may
/// collect more and try again.
fn solve_for_logarithm(
    curve: &HyperellipticCurveP,
    relations: &[HecRelation],
    m: usize,
    n: &BigUint,
    params: &HecIndexCalculusParams,
    report: &mut HecIndexCalculusReport,
) -> Option<BigUint> {
    let _ = curve;
    // Unknowns: (y_1, …, y_m, k), where y_i = log_{D₁} F_i.
    // Row: Σ c_i y_i − b k ≡ a  (mod n).
    let mut rows: Vec<Vec<(usize, BigUint)>> = Vec::with_capacity(relations.len());
    let mut rhs: Vec<BigUint> = Vec::with_capacity(relations.len());
    for rel in relations {
        let mut row: HashMap<usize, BigUint> = HashMap::new();
        for &(j, c) in &rel.entries {
            let e = row.entry(j).or_insert_with(BigUint::zero);
            *e = (&*e + signed_mod(c, n)) % n;
        }
        row.insert(m, (n - &(&rel.coef_b % n)) % n);
        let mut sparse: Vec<(usize, BigUint)> =
            row.into_iter().filter(|(_, v)| !v.is_zero()).collect();
        sparse.sort_unstable_by_key(|&(j, _)| j);
        rows.push(sparse);
        rhs.push(&rel.coef_a % n);
    }

    let k = match params.linear_algebra {
        LinearAlgebra::Sparse => {
            let (k, ops) = sparse_solve_for_k(rows, rhs, m, n)?;
            report.solve_row_ops = ops;
            k
        }
        LinearAlgebra::Dense => {
            let mut matrix: Vec<Vec<BigUint>> = Vec::with_capacity(rows.len());
            for row in &rows {
                let mut dense = vec![BigUint::zero(); m + 1];
                for (j, v) in row {
                    dense[*j] = v.clone();
                }
                matrix.push(dense);
            }
            report.solve_row_ops = rows.len() * (m + 1) * (m + 1);
            let solution = gaussian_eliminate_mod_n(&mut matrix, &mut rhs, n)?;
            solution[m].clone()
        }
    };

    Some(k)
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
pub fn factorise_small(mut n: BigUint) -> Vec<(BigUint, u32)> {
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

/// Order of `d` in `Jac(C)(F_p)` **without** knowing the group order,
/// by Baby-step Giant-step over the Hasse–Weil interval.
///
/// `#Jac(C)(F_p) ∈ [(sqrt(p) − 1)^{2g}, (sqrt(p) + 1)^{2g}]`, so some
/// multiple of `ord(d)` lies there; BSGS finds one in `O(sqrt(W))`
/// group operations, `W` being the interval width (`≈ 4g·p^{g−1/2}`).
/// Factoring that multiple and dividing out gives the exact order.
///
/// This exists because the genus-2 route to the group order goes
/// through an `L`-polynomial that is genus-2 only.  Counting the
/// Jacobian is not the point of this module; having a prime-order
/// subgroup to work in is, and that does not require the group order.
/// The Hasse–Weil interval `[(sqrt(p) − 1)^{2g}, (sqrt(p) + 1)^{2g}]`
/// containing `#Jac(C)(F_p)`.
///
/// Useful without any point counting: it says how large a subgroup has
/// to be to carry most of the group, which is the condition that makes
/// an index-calculus-vs-rho comparison meaningful — rho searches the
/// subgroup, index calculus pays for a factor base sized by the whole
/// curve, so a tiny subgroup in a large Jacobian flatters rho for a
/// reason that has nothing to do with either algorithm.
pub fn hasse_interval(curve: &HyperellipticCurveP) -> Option<(u128, u128)> {
    let g = curve.genus as i32;
    let p_f = curve.p.to_u64_digits().first().copied()? as f64;
    let root = p_f.sqrt();
    let lo = (root - 1.0).powi(2 * g).floor().max(1.0) as u128;
    let hi = (root + 1.0).powi(2 * g).ceil() as u128;
    Some((lo, hi))
}

pub fn divisor_order_bsgs(
    curve: &HyperellipticCurveP,
    d: &MumfordDivisorP,
    max_baby_steps: usize,
) -> Option<BigUint> {
    let (lo, hi) = hasse_interval(curve)?;
    let width = hi.saturating_sub(lo).max(1);

    let steps = (width as f64).sqrt().ceil() as usize + 1;
    if steps > max_baby_steps {
        return None; // would cost more than the caller allowed
    }

    // Baby steps: −j·D for j ∈ [0, steps), keyed by Mumford rep.
    let neg_d = d.neg(curve);
    let mut table: HashMap<DivisorKey, usize> = HashMap::new();
    let mut acc = MumfordDivisorP::identity(curve.p.clone());
    for j in 0..steps {
        table.entry(divisor_key(&acc)).or_insert(j);
        acc = acc.add(&neg_d, curve);
    }

    // Giant steps: (lo + i·steps)·D.
    let stride = d.scalar_mul(&BigUint::from(steps as u64), curve);
    let mut cur = d.scalar_mul(&BigUint::from(lo), curve);
    let giant_max = (width as usize) / steps + 2;
    for i in 0..=giant_max {
        if let Some(&j) = table.get(&divisor_key(&cur)) {
            // (lo + i·steps)·D = −j·D  ⟹  (lo + i·steps + j)·D = 0.
            let m = BigUint::from(lo) + BigUint::from((i * steps + j) as u64);
            if m.is_zero() {
                continue;
            }
            if !d.scalar_mul(&m, curve).is_identity() {
                continue; // f64 interval arithmetic slipped; keep going
            }
            return Some(exact_order_from_multiple(curve, d, &m));
        }
        cur = cur.add(&stride, curve);
    }
    None
}

/// A hashable stand-in for a reduced divisor: its Mumford coefficients.
type DivisorKey = (Vec<Vec<u32>>, Vec<Vec<u32>>);

fn divisor_key(d: &MumfordDivisorP) -> DivisorKey {
    (
        d.u.coeffs.iter().map(|c| c.to_u32_digits()).collect(),
        d.v.coeffs.iter().map(|c| c.to_u32_digits()).collect(),
    )
}

/// Given `m` with `m·d = 0`, divide out prime factors while the result
/// still annihilates `d`.
fn exact_order_from_multiple(
    curve: &HyperellipticCurveP,
    d: &MumfordDivisorP,
    m: &BigUint,
) -> BigUint {
    let mut order = m.clone();
    for (q, e) in factorise_small(m.clone()) {
        for _ in 0..e {
            let candidate = &order / &q;
            if candidate > BigUint::zero() && d.scalar_mul(&candidate, curve).is_identity() {
                order = candidate;
            } else {
                break;
            }
        }
    }
    order
}

/// Largest prime factor of `n`, by trial division.  Toy sizes only.
pub fn largest_prime_factor(n: &BigUint) -> BigUint {
    factorise_small(n.clone())
        .into_iter()
        .map(|(q, _)| q)
        .max()
        .unwrap_or_else(BigUint::one)
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
            FpPoly::from_coeffs(
                vec![(&p - BigUint::from(c)) % &p, BigUint::one()],
                p.clone(),
            )
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
            vec![
                (&p - BigUint::from(r)) % &p,
                BigUint::zero(),
                BigUint::one(),
            ],
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
            search: RelationSearch::Random,
            linear_algebra: LinearAlgebra::Dense,
            smoothness: SmoothnessTest::Scan,
            large_primes: false,
        };
        let (found, report) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &params);
        assert_eq!(found, Some(k), "report: {report:?}");
        assert!(report.solved);
        assert!(report.smoothness_rate() > 0.0);
    }

    /// Build a solvable instance on the toy curve: `(D₁, D₂, N, k)`.
    fn toy_instance(
        p: u64,
    ) -> (
        HyperellipticCurveP,
        MumfordDivisorP,
        MumfordDivisorP,
        BigUint,
        BigUint,
    ) {
        let curve = toy_curve(p);
        let jac = brute_force_jac_order_via_lpoly(&curve);
        let fb = build_factor_base(&curve, usize::MAX);
        let l = factorise_small(jac.clone())
            .into_iter()
            .map(|(q, _)| q)
            .max()
            .expect("non-trivial Jacobian");
        let d1 = subgroup_generator(&curve, &fb, &jac, &l).expect("generator exists");
        let k = &l / BigUint::from(3u32) + BigUint::from(7u32);
        let d2 = d1.scalar_mul(&k, &curve);
        (curve, d1, d2, l, k)
    }

    /// `C : y² = x⁷ + x³ + c x + 1` over `F_p`, genus 3.
    fn genus3_curve(p: u64, c: u64) -> HyperellipticCurveP {
        let p = BigUint::from(p);
        let f = FpPoly::from_coeffs(
            vec![
                BigUint::one(),
                BigUint::from(c) % &p,
                BigUint::zero(),
                BigUint::one(),
                BigUint::zero(),
                BigUint::zero(),
                BigUint::zero(),
                BigUint::one(),
            ],
            p.clone(),
        );
        HyperellipticCurveP::new(p, f, 3)
    }

    #[test]
    fn bsgs_order_agrees_with_the_l_polynomial_at_genus_two() {
        // Two independent routes to the same number: BSGS over the
        // Hasse-Weil interval, which knows nothing about the curve, and
        // the genus-2 L-polynomial, which counts points.  They are
        // wrong in different ways, so agreement is worth more than
        // either alone — and BSGS is the only one available at genus 3,
        // where the measurement actually goes.
        for p_u in [23u64, 41, 61] {
            let curve = toy_curve(p_u);
            let jac = brute_force_jac_order_via_lpoly(&curve);
            let fb = build_factor_base(&curve, usize::MAX);
            for e in fb.entries.iter().take(3) {
                let via_bsgs = divisor_order_bsgs(&curve, &e.divisor, 50_000)
                    .expect("BSGS should find an order at this size");
                // The L-polynomial route: the order divides #Jac, and
                // dividing out prime factors gives the same answer.
                let via_lpoly = exact_order_from_multiple(&curve, &e.divisor, &jac);
                assert_eq!(
                    via_bsgs, via_lpoly,
                    "p={p_u}: BSGS {via_bsgs} vs L-polynomial {via_lpoly}"
                );
                assert!((&jac % &via_bsgs).is_zero(), "order must divide #Jac");
            }
        }
    }

    /// A solvable genus-3 instance on `y² = x⁷ + x³ + c x + 1`.
    fn genus3_instance(
        p: u64,
    ) -> (
        HyperellipticCurveP,
        MumfordDivisorP,
        MumfordDivisorP,
        BigUint,
        BigUint,
    ) {
        for c in 1..12u64 {
            let curve = genus3_curve(p, c);
            let fb = build_factor_base(&curve, usize::MAX);
            for e in fb.entries.iter().take(4) {
                let order = match divisor_order_bsgs(&curve, &e.divisor, 50_000) {
                    Some(o) => o,
                    None => continue,
                };
                let l = largest_prime_factor(&order);
                if l < BigUint::from(500u32) {
                    continue;
                }
                let d1 = e.divisor.scalar_mul(&(&order / &l), &curve);
                if divisor_order_bsgs(&curve, &d1, 50_000).as_ref() != Some(&l) {
                    continue;
                }
                let k = &l / BigUint::from(3u32) + BigUint::from(4u32);
                let d2 = d1.scalar_mul(&k, &curve);
                return (curve, d1, d2, l, k);
            }
        }
        panic!("no usable genus-3 instance at p={p}");
    }

    #[test]
    fn solves_a_genus_three_dlp() {
        // Genus 3 is where the asymptotic crossover against rho is
        // supposed to live, so the pipeline has to work there before
        // any measurement of it means anything.  No L-polynomial is
        // available at this genus, so the order comes from BSGS over
        // the Hasse–Weil interval.
        let mut solved = false;
        // Scan `c` for a curve with a usable prime-order subgroup:
        // most genus-3 Jacobians here factor into small primes, and a
        // 37-element subgroup measures nothing.
        'outer: for c in 1..6u64 {
            let curve = genus3_curve(41, c);
            let fb = build_factor_base(&curve, usize::MAX);
            for e in fb.entries.iter().take(4) {
                let order = match divisor_order_bsgs(&curve, &e.divisor, 20_000) {
                    Some(o) => o,
                    None => continue,
                };
                assert!(e.divisor.scalar_mul(&order, &curve).is_identity());
                let l = largest_prime_factor(&order);
                if l < BigUint::from(500u32) {
                    continue;
                }
                let d1 = e.divisor.scalar_mul(&(&order / &l), &curve);
                assert_eq!(
                    divisor_order_bsgs(&curve, &d1, 20_000).as_ref(),
                    Some(&l),
                    "cofactor multiple should have exactly order l"
                );
                let k = &l / BigUint::from(3u32) + BigUint::from(4u32);
                let d2 = d1.scalar_mul(&k, &curve);
                let params = HecIndexCalculusParams {
                    fb_size: usize::MAX,
                    extra_relations: 8,
                    max_trials: 500_000,
                    seed: 4,
                    search: RelationSearch::walk(),
                    linear_algebra: LinearAlgebra::Sparse,
                    smoothness: SmoothnessTest::Gcd,
                    large_primes: false,
                };
                let (found, report) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &params);
                assert_eq!(found, Some(k), "genus-3 report: {report:?}");
                // Smoothness at genus 3 should be near 1/3! = 0.167, well
                // below genus 2's 0.5 — the cost the extra genus buys.
                assert!(
                    report.smoothness_rate() < 0.45,
                    "genus-3 smoothness {} looks like genus 2",
                    report.smoothness_rate()
                );
                solved = true;
                break 'outer;
            }
        }
        assert!(solved, "no usable genus-3 subgroup found");
    }

    #[test]
    fn cantor_zassenhaus_roots_agree_with_the_scan_at_degree_three() {
        // deg u = 3 is what genus 3 produces and what the closed forms
        // do not cover, so the CZ path is the one carrying the fast
        // oracle there.  Exhaustive over monic cubics for two primes.
        // Both residue classes mod 4.  The discriminant pre-filter
        // carries a sign `(−1)^{n(n−1)/2}`, and whether −1 is itself a
        // square depends on `p mod 4` — so a wrong sign would reject
        // split polynomials at `p ≡ 3 (mod 4)` and pass at
        // `p ≡ 1 (mod 4)`, or the reverse.  Testing one class only would
        // miss it.
        for p_u in [11u64, 13, 23, 29] {
            let p = BigUint::from(p_u);
            for a2 in 0..p_u {
                for a1 in 0..p_u {
                    for a0 in 0..p_u {
                        let u = FpPoly::from_coeffs(
                            vec![
                                BigUint::from(a0),
                                BigUint::from(a1),
                                BigUint::from(a2),
                                BigUint::one(),
                            ],
                            p.clone(),
                        );
                        let (mut o1, mut o2) = (0usize, 0usize);
                        assert_eq!(
                            split_by_scan(&u, &p, &mut o1),
                            split_by_gcd(&u, &p, &mut o2),
                            "p={p_u} u=x^3+{a2}x^2+{a1}x+{a0}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn gcd_oracle_agrees_with_the_scan_on_every_monic_quartic() {
        // Degree 4 is what genus 4 produces, and it is where the
        // discriminant pre-filter has to be right about a different set
        // of cycle types than the cubic case: the odd ones are (2,1,1)
        // and (4).  Exhaustive over monic quartics at two primes, one in
        // each class mod 4.
        for p_u in [11u64, 13] {
            let p = BigUint::from(p_u);
            let mut split_count = 0usize;
            let mut rejected = 0usize;
            for a3 in 0..p_u {
                for a2 in 0..p_u {
                    for a1 in 0..p_u {
                        for a0 in 0..p_u {
                            let u = FpPoly::from_coeffs(
                                vec![
                                    BigUint::from(a0),
                                    BigUint::from(a1),
                                    BigUint::from(a2),
                                    BigUint::from(a3),
                                    BigUint::one(),
                                ],
                                p.clone(),
                            );
                            let (mut o1, mut o2) = (0usize, 0usize);
                            let scan = split_by_scan(&u, &p, &mut o1);
                            let gcd = split_by_gcd(&u, &p, &mut o2);
                            assert_eq!(scan, gcd, "p={p_u} u=x^4+{a3}x^3+{a2}x^2+{a1}x+{a0}");
                            if scan.is_some() {
                                split_count += 1;
                            } else {
                                rejected += 1;
                            }
                        }
                    }
                }
            }
            // Sanity on the population itself: completely split monic
            // quartics should be about 1/4! of them.
            let total = (split_count + rejected) as f64;
            let rate = split_count as f64 / total;
            assert!(
                rate > 0.02 && rate < 0.08,
                "p={p_u}: split rate {rate:.3} is nowhere near 1/24"
            );
        }
    }

    #[test]
    fn the_discriminant_filter_never_rejects_a_split_polynomial() {
        // The filter is the one part of the oracle that returns "not
        // smooth" without looking for roots at all, so it is the one
        // that can silently throw away relations.  Build polynomials
        // that split by construction and require every one to survive
        // it.
        for p_u in [11u64, 13, 41, 61] {
            let p = BigUint::from(p_u);
            for r0 in 0..p_u.min(9) {
                for r1 in 0..p_u.min(9) {
                    for r2 in 0..p_u.min(9) {
                        for r3 in 0..p_u.min(5) {
                            let lin = |r: u64| {
                                FpPoly::from_coeffs(
                                    vec![(&p - BigUint::from(r)) % &p, BigUint::one()],
                                    p.clone(),
                                )
                            };
                            for u in [
                                lin(r0).mul(&lin(r1)).mul(&lin(r2)),
                                lin(r0).mul(&lin(r1)).mul(&lin(r2)).mul(&lin(r3)),
                            ] {
                                let mut ops = 0usize;
                                let disc = discriminant(&u, &p, &mut ops);
                                assert!(
                                    disc.is_zero() || is_square_mod(&disc, &p, &mut ops),
                                    "p={p_u}: split poly {:?} has non-square discriminant {disc}",
                                    u.coeffs
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn gcd_oracle_agrees_with_the_scan_on_every_monic_quadratic() {
        // Exhaustive over all monic u of degree ≤ 2 for several p: the
        // fast oracle is only worth having if it is the *same* oracle.
        // A disagreement in the "smooth" direction would put a row in
        // the relation matrix that is not a relation.
        for p_u in [11u64, 13, 41, 61] {
            let p = BigUint::from(p_u);
            for a1 in 0..p_u {
                for a0 in 0..p_u {
                    for deg in [1usize, 2] {
                        let u = if deg == 1 {
                            FpPoly::from_coeffs(vec![BigUint::from(a0), BigUint::one()], p.clone())
                        } else {
                            FpPoly::from_coeffs(
                                vec![BigUint::from(a0), BigUint::from(a1), BigUint::one()],
                                p.clone(),
                            )
                        };
                        let (mut o1, mut o2) = (0usize, 0usize);
                        let scan = split_by_scan(&u, &p, &mut o1);
                        let gcd = split_by_gcd(&u, &p, &mut o2);
                        assert_eq!(
                            scan, gcd,
                            "p={p_u} u={:?}: scan {scan:?} vs gcd {gcd:?}",
                            u.coeffs
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn gcd_oracle_is_cheaper_than_the_scan_and_gets_the_same_answer() {
        let (curve, d1, d2, l, k) = toy_instance(251);
        let scan = HecIndexCalculusParams {
            fb_size: usize::MAX,
            extra_relations: 8,
            max_trials: 200_000,
            seed: 11,
            search: RelationSearch::factor_base_walk(),
            linear_algebra: LinearAlgebra::Sparse,
            smoothness: SmoothnessTest::Scan,
            large_primes: false,
        };
        let fast = HecIndexCalculusParams {
            smoothness: SmoothnessTest::Gcd,
            ..scan.clone()
        };
        let (scan_k, scan_rep) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &scan);
        let (fast_k, fast_rep) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &fast);

        assert_eq!(scan_k.as_ref(), Some(&k));
        assert_eq!(fast_k.as_ref(), Some(&k), "fast report: {fast_rep:?}");
        // Same seed and the same oracle answers mean the same walk and
        // the same relations, so the trial counts must match exactly —
        // this is what makes the field-op comparison a like-for-like.
        assert_eq!(scan_rep.trials, fast_rep.trials);
        assert!(
            fast_rep.smoothness_field_ops * 4 < scan_rep.smoothness_field_ops,
            "gcd {} vs scan {} field ops",
            fast_rep.smoothness_field_ops,
            scan_rep.smoothness_field_ops
        );
    }

    /// Check a collected relation is a true identity in the Jacobian:
    /// `Σ c_i F_i = a·D₁ + b·D₂`.
    fn relation_holds(
        curve: &HyperellipticCurveP,
        fb: &HecFactorBase,
        d1: &MumfordDivisorP,
        d2: &MumfordDivisorP,
        rel: &HecRelation,
    ) -> bool {
        let mut lhs = MumfordDivisorP::identity(curve.p.clone());
        for &(i, c) in &rel.entries {
            let base = &fb.entries[i].divisor;
            let term = if c >= 0 {
                base.scalar_mul(&BigUint::from(c as u64), curve)
            } else {
                base.neg(curve)
                    .scalar_mul(&BigUint::from(c.unsigned_abs()), curve)
            };
            lhs = lhs.add(&term, curve);
        }
        let rhs = d1
            .scalar_mul(&rel.coef_a, curve)
            .add(&d2.scalar_mul(&rel.coef_b, curve), curve);
        lhs == rhs
    }

    #[test]
    fn every_factor_base_walk_relation_is_a_true_identity() {
        // The factor-base walk subtracts the steps it took from the
        // decomposition it found.  Get that bookkeeping wrong by one and
        // the rows are not relations at all — the solve would then
        // return a `k` that fails verification, which reads as "no
        // solution" rather than as a bug.  So check the rows directly,
        // at both genera.
        for (curve, d1, d2, _l, _k) in [toy_instance(41), genus3_instance(41)] {
            let fb = build_factor_base(&curve, usize::MAX);
            let params = HecIndexCalculusParams {
                fb_size: usize::MAX,
                extra_relations: 4,
                max_trials: 200_000,
                seed: 99,
                search: RelationSearch::factor_base_walk(),
                linear_algebra: LinearAlgebra::Sparse,
                smoothness: SmoothnessTest::Gcd,
                large_primes: false,
            };
            let mut report = HecIndexCalculusReport::default();
            let n = &_l;
            let rels =
                collect_relations(&curve, &d1, &d2, n, &fb, fb.len() + 5, &params, &mut report);
            assert!(!rels.is_empty());
            // A walked row carries step counts, so it must be wider than
            // a bare decomposition at least once, or the test is not
            // exercising the bookkeeping it claims to.
            assert!(
                rels.iter().any(|r| r.entries.len() > curve.genus as usize),
                "no row carried step counts"
            );
            for (i, rel) in rels.iter().enumerate() {
                assert!(
                    relation_holds(&curve, &fb, &d1, &d2, rel),
                    "relation {i} is not an identity: {rel:?}"
                );
            }
        }
    }

    #[test]
    fn rank_deficient_rows_are_recovered_by_collecting_more() {
        // Seed 20260919 at p = 41 produced 24 factor-base-walk rows that
        // pinned every `y_i` and left `k` free: walked rows are
        // correlated, so a fixed margin of spare rows is sometimes not
        // enough.  The driver now collects more and solves again, so
        // this seed must succeed.
        let (curve, d1, d2, l, k) = toy_instance(41);
        for seed in [20260919u64, 20260916, 20260917, 20260918, 20260920] {
            let params = HecIndexCalculusParams {
                fb_size: usize::MAX,
                extra_relations: 8,
                max_trials: 500_000,
                seed,
                search: RelationSearch::factor_base_walk(),
                linear_algebra: LinearAlgebra::Sparse,
                smoothness: SmoothnessTest::Gcd,
                large_primes: false,
            };
            let (got, report) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &params);
            assert_eq!(got.as_ref(), Some(&k), "seed {seed}: {report:?}");
        }
    }

    #[test]
    fn factor_base_walk_pays_no_step_precomputation() {
        let (curve, d1, d2, l, k) = genus3_instance(41);
        let base = HecIndexCalculusParams {
            fb_size: usize::MAX,
            extra_relations: 8,
            max_trials: 500_000,
            seed: 5,
            search: RelationSearch::walk(),
            linear_algebra: LinearAlgebra::Sparse,
            smoothness: SmoothnessTest::Gcd,
            large_primes: false,
        };
        let fbw = HecIndexCalculusParams {
            search: RelationSearch::factor_base_walk(),
            ..base.clone()
        };
        let (k_old, old) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &base);
        let (k_new, new) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &fbw);
        assert_eq!(k_old.as_ref(), Some(&k));
        assert_eq!(k_new.as_ref(), Some(&k), "factor-base walk: {new:?}");
        // The only precomputation left is the starting point and the
        // restart jumps, so it cannot come near 64 branch divisors.
        assert!(
            new.precompute_ops * 3 < old.precompute_ops,
            "factor-base walk precompute {} vs {} for the (a,b)-step walk",
            new.precompute_ops,
            old.precompute_ops
        );
        assert!(new.ops_per_trial() < 1.5, "{new:?}");
    }

    #[test]
    fn walk_solves_and_costs_far_less_than_redrawing() {
        let (curve, d1, d2, l, k) = toy_instance(41);
        let base = HecIndexCalculusParams {
            fb_size: usize::MAX,
            extra_relations: 5,
            max_trials: 50_000,
            seed: 20260916,
            search: RelationSearch::Random,
            linear_algebra: LinearAlgebra::Sparse,
            smoothness: SmoothnessTest::Gcd,
            large_primes: false,
        };
        let (rnd_k, rnd) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &base);
        let walked = HecIndexCalculusParams {
            search: RelationSearch::walk(),
            ..base
        };
        let (walk_k, walk) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &walked);

        assert_eq!(rnd_k.as_ref(), Some(&k));
        assert_eq!(walk_k.as_ref(), Some(&k), "walk report: {walk:?}");
        // The walk pays one group operation per trial where the random
        // draw pays 2⌈log₂ N⌉, so on comparable trial counts it must
        // come out well ahead.  This is the whole claim of the walk; if
        // it ever stops holding, the walk has stopped being a walk.
        // Per trial, excluding precomputation: the walk pays one group
        // operation, the random draw pays 2⌈log₂ N⌉.  Total ops are the
        // wrong comparison at this toy `N`, where the walk's 32-branch
        // precomputation has barely amortised — which is itself the
        // reason the report separates the two.
        assert!(walk.ops_per_trial() < 1.5, "walk {:?}", walk);
        assert!(
            rnd.ops_per_trial() > 4.0 * walk.ops_per_trial(),
            "random {:.1} vs walk {:.1} ops/trial",
            rnd.ops_per_trial(),
            walk.ops_per_trial()
        );
    }

    #[test]
    fn sparse_solve_agrees_with_dense_and_is_cheaper() {
        let (curve, d1, d2, l, k) = toy_instance(41);
        let dense = HecIndexCalculusParams {
            fb_size: usize::MAX,
            extra_relations: 5,
            max_trials: 50_000,
            seed: 7,
            search: RelationSearch::walk(),
            linear_algebra: LinearAlgebra::Dense,
            smoothness: SmoothnessTest::Scan,
            large_primes: false,
        };
        let sparse = HecIndexCalculusParams {
            linear_algebra: LinearAlgebra::Sparse,
            ..dense.clone()
        };
        let (dense_k, dense_rep) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &dense);
        let (sparse_k, sparse_rep) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &sparse);

        // Same seed, same relations, so the two solvers see the same
        // system and must agree on its answer.
        assert_eq!(dense_k.as_ref(), Some(&k));
        assert_eq!(sparse_k, dense_k, "sparse report: {sparse_rep:?}");
        assert!(
            sparse_rep.solve_row_ops < dense_rep.solve_row_ops,
            "sparse {} vs dense {} mul-mods",
            sparse_rep.solve_row_ops,
            dense_rep.solve_row_ops
        );
    }

    /// Rank of `rows` (dense, over `F_n`, `n < 2^32`), by plain Gaussian
    /// elimination — an oracle written independently of both solvers.
    fn rank_mod(mut rows: Vec<Vec<u64>>, n: u64) -> usize {
        let cols = rows.first().map_or(0, |r| r.len());
        let inv = |a: u64| {
            // Fermat: n is prime.
            let (mut base, mut e, mut acc) = (a % n, n - 2, 1u64);
            while e > 0 {
                if e & 1 == 1 {
                    acc = ((acc as u128 * base as u128) % n as u128) as u64;
                }
                base = ((base as u128 * base as u128) % n as u128) as u64;
                e >>= 1;
            }
            acc
        };
        let mut rank = 0;
        for c in 0..cols {
            let Some(p) = (rank..rows.len()).find(|&r| rows[r][c] != 0) else {
                continue;
            };
            rows.swap(rank, p);
            let iv = inv(rows[rank][c]);
            for r in 0..rows.len() {
                if r != rank && rows[r][c] != 0 {
                    let f = ((rows[r][c] as u128 * iv as u128) % n as u128) as u64;
                    for cc in 0..cols {
                        let t = ((f as u128 * rows[rank][cc] as u128) % n as u128) as u64;
                        rows[r][cc] = (rows[r][cc] + n - t) % n;
                    }
                }
            }
            rank += 1;
        }
        rank
    }

    /// The sparse solver on planted random systems, many more of them
    /// than a DLP run produces, against an independent rank oracle: `k`
    /// is determined exactly when `e_k` lies in the row space, i.e. when
    /// appending `e_k` does not raise the rank.  The solver must answer
    /// exactly then, and with the planted `k`.  Small `N` is included
    /// because that is where eliminations cancel to zero and a merge that
    /// mishandles a cancelled or filled-in column would show; `N = 3, 5, 7`
    /// make a row cancel out of a column and fill back into it often
    /// enough to catch a row visited twice, which larger `N` never did.
    ///
    /// Dense elimination is *not* the oracle here: on an under-determined
    /// system `gaussian_eliminate_mod_n` returns a particular solution
    /// with the free unknowns at zero rather than `None`, so "dense
    /// answered" does not mean "`k` is determined".  The DLP driver is
    /// protected by its `[k]D₁ = D₂` check; this test would not be.
    #[test]
    fn sparse_solve_matches_dense_on_planted_systems() {
        let mut rng = StdRng::seed_from_u64(0x5a_0e5e);
        let mut answered = 0;
        for &n_u64 in &[3u64, 5, 7, 101, 1009, 1_000_003, 16_790_591] {
            let n = BigUint::from(n_u64);
            for trial in 0..200 {
                let m = 5 + (rng.next_u64() % 40) as usize;
                // Half the systems are short of rows, so `k` is often
                // undetermined and the solver must decline.
                let rows_wanted = if trial % 2 == 0 {
                    m + 1 + (rng.next_u64() % 6) as usize
                } else {
                    (rng.next_u64() as usize % (m + 1)) + 1
                };
                let y: Vec<BigUint> = (0..=m).map(|_| rand_below(&mut rng, &n)).collect();
                let mut rows = Vec::new();
                let mut rhs = Vec::new();
                for _ in 0..rows_wanted {
                    let weight = 1 + (rng.next_u64() % 4) as usize;
                    let mut cols: Vec<usize> = (0..weight)
                        .map(|_| (rng.next_u64() % m as u64) as usize)
                        .collect();
                    cols.push(m);
                    cols.sort_unstable();
                    cols.dedup();
                    let row: Vec<(usize, BigUint)> = cols
                        .into_iter()
                        .map(|j| {
                            let v = BigUint::from(1 + rng.next_u64() % (n_u64 - 1));
                            (j, v)
                        })
                        .collect();
                    let b = row
                        .iter()
                        .fold(BigUint::zero(), |acc, (j, v)| (acc + v * &y[*j]) % &n);
                    rows.push(row);
                    rhs.push(b);
                }

                let dense: Vec<Vec<u64>> = rows
                    .iter()
                    .map(|row| {
                        let mut d = vec![0u64; m + 1];
                        for (j, v) in row {
                            d[*j] = v.to_u64_digits().first().copied().unwrap_or(0);
                        }
                        d
                    })
                    .collect();
                let mut with_ek = dense.clone();
                let mut ek = vec![0u64; m + 1];
                ek[m] = 1;
                with_ek.push(ek);
                let determined = rank_mod(with_ek, n_u64) == rank_mod(dense, n_u64);

                let sparse_k = sparse_solve_for_k(rows, rhs, m, &n).map(|(k, _)| k);
                let at = format!("n={n_u64} m={m} trial={trial}");
                if determined {
                    assert_eq!(sparse_k.as_ref(), Some(&y[m]), "declined or wrong: {at}");
                    answered += 1;
                } else {
                    assert_eq!(sparse_k, None, "answered an undetermined k: {at}");
                }
            }
        }
        // Guard against a test that passes by never answering.
        assert!(answered >= 400, "only {answered} systems answered");
    }

    /// Every relation collected from a halved factor base with large
    /// primes on — the combined ones included — is a true identity in
    /// the Jacobian.  A wrong sign or scale in the combination would
    /// still often solve (the retry loop and the final `[k]D₁ = D₂`
    /// check hide a few bad rows), so the rows are checked directly.
    #[test]
    fn large_prime_relations_are_jacobian_identities() {
        for (g, (curve, d1, d2, l, _k)) in [(2, toy_instance(61)), (3, genus3_instance(31))] {
            let full = build_factor_base(&curve, usize::MAX).len();
            let fb = build_factor_base(&curve, full / 2);
            let params = HecIndexCalculusParams {
                fb_size: full / 2,
                large_primes: true,
                seed: 11,
                ..HecIndexCalculusParams::default()
            };
            let mut report = HecIndexCalculusReport::default();
            let rels = collect_relations(
                &curve,
                &d1,
                &d2,
                &l,
                &fb,
                fb.len() + 9,
                &params,
                &mut report,
            );
            assert!(
                report.combined_relations > 0,
                "genus {g}: nothing combined: {report:?}"
            );
            for (i, rel) in rels.iter().enumerate() {
                assert!(
                    relation_holds(&curve, &fb, &d1, &d2, rel),
                    "genus {g}: relation {i} of {} is not an identity: {rel:?}",
                    rels.len()
                );
            }
        }
    }

    /// A halved factor base with large primes still recovers `k`, and
    /// relies on combined relations to do it.
    #[test]
    fn large_primes_solve_with_a_reduced_base() {
        let (curve, d1, d2, l, k) = toy_instance(61);
        let full = build_factor_base(&curve, usize::MAX).len();
        for seed in 0..5 {
            let params = HecIndexCalculusParams {
                fb_size: full / 2,
                large_primes: true,
                seed,
                ..HecIndexCalculusParams::default()
            };
            let (got, rep) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &params);
            assert_eq!(got.as_ref(), Some(&k), "seed {seed}: {rep:?}");
            assert!(rep.combined_relations > 0, "seed {seed}: {rep:?}");
        }
    }

    /// With the full factor base there are no off-base places, so large
    /// primes must change nothing at all: same relations, same counted
    /// work, same `k`.
    #[test]
    fn large_primes_are_inert_on_the_full_base() {
        let (curve, d1, d2, l, k) = toy_instance(61);
        let off = HecIndexCalculusParams {
            seed: 3,
            ..HecIndexCalculusParams::default()
        };
        let on = HecIndexCalculusParams {
            large_primes: true,
            ..off.clone()
        };
        let (k_off, r_off) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &off);
        let (k_on, r_on) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &on);
        assert_eq!(k_off.as_ref(), Some(&k));
        assert_eq!(k_on, k_off);
        assert_eq!(
            (
                r_on.trials,
                r_on.relations,
                r_on.jacobian_ops,
                r_on.smoothness_field_ops,
                r_on.solve_row_ops
            ),
            (
                r_off.trials,
                r_off.relations,
                r_off.jacobian_ops,
                r_off.smoothness_field_ops,
                r_off.solve_row_ops
            )
        );
        assert_eq!((r_on.partial_relations, r_on.combined_relations), (0, 0));
    }

    #[test]
    fn the_smoothness_oracle_is_charged() {
        let (curve, d1, d2, l, _k) = toy_instance(41);
        let params = HecIndexCalculusParams {
            fb_size: usize::MAX,
            extra_relations: 5,
            max_trials: 50_000,
            seed: 3,
            search: RelationSearch::factor_base_walk(),
            linear_algebra: LinearAlgebra::Sparse,
            smoothness: SmoothnessTest::Gcd,
            large_primes: false,
        };
        let (_k, report) = hec_index_calculus_dlp(&curve, &d1, &d2, &l, &params);
        // Every trial runs the oracle, so its field-operation count can
        // never be zero while trials were taken — a zero here would mean
        // the cost had gone missing from the accounting, not that it
        // was not paid.
        assert!(report.trials > 0);
        assert!(
            report.smoothness_field_ops >= report.trials,
            "oracle charged {} mul-mods over {} trials",
            report.smoothness_field_ops,
            report.trials
        );
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
            search: RelationSearch::Random,
            linear_algebra: LinearAlgebra::Dense,
            smoothness: SmoothnessTest::Scan,
            large_primes: false,
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
