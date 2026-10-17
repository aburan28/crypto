//! The volcano level at depth 1 from points over `GF(q^k)`, for primes `ℓ`
//! whose division polynomial is too large for [`super::torsion`].
//!
//! At depth 1 Frobenius has the double eigenvalue `λ ≡ t/2` on `E[ℓ]` and
//! is either the scalar `λ` (the crater) or a Jordan block (the floor).
//! With `k` the order of `λ` mod `ℓ`, the scalar satisfies `π^k = 1` on
//! `E[ℓ]` and a Jordan block does not (`ℓ ∤ k`), so
//!
//! > Frobenius is a scalar on `E[ℓ]` ⟺ `E[ℓ] ⊆ E(GF(q^k))` ⟺ the
//! > `ℓ`-part of `E(GF(q^k))` has rank 2.
//!
//! `N = #E(GF(q^k)) = q^k + 1 − V_k(t, q)` is exact from the trace.  With
//! `ℓ^v ∥ N` the `ℓ`-part is `Z/ℓ^a × Z/ℓ^b`, `a + b = v`, and either
//! outcome has a certificate found from seeded random points:
//!
//! - **rank 1** (the floor): a point whose `ℓ`-part has order `ℓ^v`, which
//!   only a cyclic `ℓ`-part has;
//! - **rank 2** (the crater): two points of order `ℓ`, the second outside
//!   the group the first generates.
//!
//! The search is randomised; whichever certificate it finds is a proof.
//! When it finds neither within [`TRIES`] points, the level stays untested.

use num_bigint::{BigInt, BigUint};
use num_traits::Zero;

use super::arith::lucas_v;
use super::torsion::{
    compose, eigen_order, jac_add, jac_double, padd, pdivrem, pinv_mod, pmul, ppow, prem, psub,
    trim, Field, Fp, FpK, Jac, Poly,
};
use super::Model;

/// Random points tried before giving up.
pub const TRIES: u64 = 24;

/// The largest extension degree `k = ord(λ)` tried.
pub const MAX_DEGREE: u64 = 64;

/// Work allowed for one test, in base-field word products: a guard on
/// build time.  The search stops once the scalar multiplications it has
/// made would exceed it, so the bound holds whatever the points do.
pub const POINT_TEST_COST: u64 = 6_000_000_000;

/// A test is not started unless the budget covers this many scalar
/// multiplications (a typical certificate needs two or three points).
const MIN_MULTIPLICATIONS: u64 = 8;

/// `GF(q^k) = F[w]/(m)` for a monic irreducible `m` of degree `k` over `F`.
pub struct Ext<'a, F: Field> {
    pub base: &'a F,
    pub modulus: Poly<F::E>,
}

impl<F: Field> Ext<'_, F> {
    fn k(&self) -> usize {
        self.modulus.len() - 1
    }

    fn constant(&self, c: &F::E) -> Poly<F::E> {
        trim(self.base, vec![c.clone()])
    }
}

impl<F: Field> Field for Ext<'_, F> {
    type E = Poly<F::E>;
    fn zero(&self) -> Poly<F::E> {
        Vec::new()
    }
    fn one(&self) -> Poly<F::E> {
        vec![self.base.one()]
    }
    fn int(&self, v: u64) -> Poly<F::E> {
        self.constant(&self.base.int(v))
    }
    fn add(&self, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
        padd(self.base, a, b)
    }
    fn neg(&self, a: &Poly<F::E>) -> Poly<F::E> {
        a.iter().map(|c| self.base.neg(c)).collect()
    }
    fn mul(&self, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
        prem(self.base, &pmul(self.base, a, b), &self.modulus)
    }
    fn inv(&self, a: &Poly<F::E>) -> Poly<F::E> {
        pinv_mod(self.base, a, &self.modulus).expect("a non-zero element of a field")
    }
    fn sample(&self, i: u64) -> Poly<F::E> {
        let k = self.k() as u64;
        trim(
            self.base,
            (0..k)
                .map(|j| self.base.sample(i.wrapping_mul(k + 1).wrapping_add(j)))
                .collect(),
        )
    }
}

fn pgcd<F: Field>(f: &F, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
    let (mut a, mut b) = (a.clone(), b.clone());
    while !b.is_empty() {
        let (_, r) = pdivrem(f, &a, &b);
        (a, b) = (b, r);
    }
    a
}

fn prime_factors(mut n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    let mut d = 2;
    while d * d <= n {
        if n.is_multiple_of(d) {
            out.push(d);
            while n.is_multiple_of(d) {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        out.push(n);
    }
    out
}

/// A monic irreducible polynomial of degree `k ≥ 2` over `F` (`|F| = q`),
/// tried from seeded random coefficients and certified by Rabin's test:
/// `x^{q^k} ≡ x (mod m)` and `gcd(x^{q^{k/r}} − x, m) = 1` for each prime
/// `r | k`.
pub fn irreducible<F: Field>(f: &F, q: &BigUint, k: usize, seed: u64) -> Option<Poly<F::E>> {
    let x = vec![f.zero(), f.one()];
    let primes = prime_factors(k);
    for attempt in 0..(64 * k as u64) {
        let mut m: Poly<F::E> = (0..k as u64)
            .map(|j| f.sample(seed.wrapping_add(attempt * (k as u64 + 1) + j)))
            .collect();
        m.push(f.one());
        let m = trim(f, m);
        // x^{q^j} for j = 1 … k, by composing with x^q.
        let x_q = ppow(f, &x, q, &m);
        let mut powers = vec![x_q.clone()];
        for _ in 1..k {
            let last = powers.last().expect("non-empty");
            powers.push(compose(f, &x_q, last, &m));
        }
        if powers[k - 1] != prem(f, &x, &m) {
            continue;
        }
        let coprime = primes.iter().all(|r| {
            let d = psub(f, &powers[k / r - 1], &x);
            pgcd(f, &m, &d).len() == 1
        });
        if coprime {
            return Some(m);
        }
    }
    None
}

/// A square root of `a` in a field of odd order `order` (Tonelli–Shanks),
/// or `None` for a non-square.
pub fn sqrt<F: Field>(f: &F, a: &F::E, order: &BigUint) -> Option<F::E> {
    if f.is_zero(a) {
        return Some(a.clone());
    }
    let pow = |b: &F::E, e: &BigUint| {
        let mut acc = f.one();
        for i in (0..e.bits()).rev() {
            acc = f.mul(&acc, &acc);
            if e.bit(i) {
                acc = f.mul(&acc, b);
            }
        }
        acc
    };
    let qm1 = order - 1u8;
    let half = &qm1 >> 1;
    let minus_one = f.neg(&f.one());
    if pow(a, &half) != f.one() {
        return None;
    }
    let s = qm1.trailing_zeros().expect("order > 1");
    let t = &qm1 >> s;
    let z = (0u64..)
        .map(|i| f.sample(i ^ 0x5eed))
        .find(|z| !f.is_zero(z) && pow(z, &half) == minus_one)?;
    let mut m = s;
    let mut c = pow(&z, &t);
    let mut tt = pow(a, &t);
    let mut r = pow(a, &((&t + 1u8) >> 1));
    while tt != f.one() {
        let mut i = 0;
        let mut sq = tt.clone();
        while sq != f.one() {
            sq = f.mul(&sq, &sq);
            i += 1;
        }
        let mut b = c.clone();
        for _ in 0..(m - i - 1) {
            b = f.mul(&b, &b);
        }
        m = i;
        c = f.mul(&b, &b);
        tt = f.mul(&tt, &c);
        r = f.mul(&r, &b);
    }
    Some(r)
}

fn infinity<F: Field>(f: &F) -> Jac<F::E> {
    Jac {
        x: f.one(),
        y: f.one(),
        z: f.zero(),
    }
}

fn mul<F: Field>(f: &F, a: &F::E, p: &Jac<F::E>, n: &BigUint) -> Jac<F::E> {
    let mut acc = infinity(f);
    for i in (0..n.bits()).rev() {
        acc = jac_double(f, a, &acc);
        if n.bit(i) {
            acc = jac_add(f, a, &acc, p);
        }
    }
    acc
}

/// Affine coordinates, `None` for the point at infinity.
fn affine<F: Field>(f: &F, p: &Jac<F::E>) -> Option<(F::E, F::E)> {
    if f.is_zero(&p.z) {
        return None;
    }
    let zi = f.inv(&p.z);
    let zi2 = f.mul(&zi, &zi);
    Some((f.mul(&p.x, &zi2), f.mul(&p.y, &f.mul(&zi2, &zi))))
}

/// Decide the `ℓ`-rank of `E(G)`, `|G| = order`, `#E(G) = n = ℓ^v·m`:
/// `Some(true)` for rank 2, `Some(false)` for rank 1.
fn rank_two<F: Field>(
    g: &F,
    a: &F::E,
    b: &F::E,
    order: &BigUint,
    n: &BigUint,
    ell: u64,
    budget: u64,
) -> Option<bool> {
    let l = BigUint::from(ell);
    // Scalar multiplications made so far, a square root counting as one.
    let mut spent = 0u64;
    let mut v = 0u32;
    let mut cof = n.clone();
    while (&cof % &l).is_zero() {
        cof /= &l;
        v += 1;
    }
    if v == 0 {
        return None;
    }
    let mut order_ell: Vec<(F::E, F::E)> = Vec::new();
    for i in 0..TRIES {
        if spent >= budget {
            return None;
        }
        let x = g.sample(i);
        let rhs = g.add(&g.mul(&g.add(&g.mul(&x, &x), a), &x), b);
        spent += 1;
        let Some(y) = sqrt(g, &rhs, order) else {
            continue;
        };
        let p = Jac { x, y, z: g.one() };
        spent += 1;
        let mut s = mul(g, a, &p, &cof);
        if g.is_zero(&s.z) {
            continue;
        }
        // s has order ℓ^e; step to the point of order ℓ below it.
        let mut e = 0u32;
        loop {
            // A multiplication by ℓ is a small fraction of one by N.
            let next = mul(g, a, &s, &l);
            e += 1;
            if g.is_zero(&next.z) {
                break;
            }
            if e == v {
                // [ℓ^v]s ≠ O: N is not the group order; refuse.
                return None;
            }
            s = next;
        }
        if e == v {
            // An ℓ-part of order ℓ^v: cyclic, rank 1.
            return Some(false);
        }
        let t = affine(g, &s).expect("order ℓ, not O");
        for (x1, y1) in &order_ell {
            let p1 = Jac {
                x: x1.clone(),
                y: y1.clone(),
                z: g.one(),
            };
            let mut multiple = p1.clone();
            let mut inside = false;
            for _ in 1..ell {
                if affine(g, &multiple).as_ref() == Some(&t) {
                    inside = true;
                    break;
                }
                multiple = jac_add(g, a, &multiple, &p1);
            }
            if !inside {
                return Some(true);
            }
        }
        order_ell.push(t);
    }
    None
}

/// Whether Frobenius acts on `E[ℓ]` as a scalar, at depth 1 (`ℓ ∥ v`),
/// decided from points over `GF(q^k)` (see the module docs).  `None` when
/// the model is not one this handles, the work would exceed
/// [`POINT_TEST_COST`], or the search finds no certificate.
pub fn scalar_by_points(model: &Model, trace: &BigInt, ell: u64) -> Option<bool> {
    // ord(λ) ≤ ℓ − 1 is found by stepping, and the extension has that
    // degree: both bound ℓ.
    if ell > super::torsion::SCALAR_TEST_MAX_ELL {
        return None;
    }
    let k = eigen_order(ell, trace);
    if k == 0 || u64::from(k) > MAX_DEGREE {
        return None;
    }
    let q = model.q();
    let qi = BigInt::from(q.clone());
    let n = (qi.pow(k) + 1u8 - lucas_v(trace, &qi, k))
        .to_biguint()
        .filter(|n| !n.is_zero())?;
    let order = q.pow(k);
    // A scalar multiple by ~N over GF(q^k) costs about 12·bits(N)
    // extension products of k² base products, each limbs² words, and the
    // search makes a few of them per point.
    let (limbs, base_k) = match model {
        Model::Prime { p, .. } => (p.bits().div_ceil(64), 1),
        Model::Extension { k: ext, .. } => (1, u64::from(*ext)),
        _ => return None,
    };
    let kk = u64::from(k) * base_k;
    let per_multiplication = 12u64
        .saturating_mul(n.bits())
        .saturating_mul(kk.saturating_mul(kk))
        .saturating_mul(limbs.saturating_mul(limbs))
        .max(1);
    let budget = POINT_TEST_COST / per_multiplication;
    if budget < MIN_MULTIPLICATIONS {
        return None;
    }
    #[allow(clippy::too_many_arguments)]
    fn over<F: Field>(
        f: &F,
        a: &F::E,
        b: &F::E,
        q: &BigUint,
        k: u32,
        order: &BigUint,
        n: &BigUint,
        ell: u64,
        budget: u64,
    ) -> Option<bool> {
        if k == 1 {
            return rank_two(f, a, b, order, n, ell, budget);
        }
        let modulus = irreducible(f, q, k as usize, 0x1e1d)?;
        let g = Ext { base: f, modulus };
        let (a, b) = (g.constant(a), g.constant(b));
        rank_two(&g, &a, &b, order, n, ell, budget)
    }
    match model {
        Model::Prime { p, a, b } => {
            if (p % ell).is_zero() || p <= &BigUint::from(3u8) {
                return None;
            }
            over(&Fp(p.clone()), a, b, &q, k, &order, &n, ell, budget)
        }
        Model::Extension {
            p, modulus, a, b, ..
        } => {
            let f = FpK::new(p, modulus)?;
            if (BigUint::from(ell) % p).is_zero() {
                return None;
            }
            let (a, b) = (f.element(a)?, f.element(b)?);
            over(&f, &a, &b, &q, k, &order, &n, ell, budget)
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::curve_traits::torsion::frobenius_is_scalar;

    /// Wherever the division-polynomial test runs, the point test must
    /// agree with it: over GF(p) for many small curves at depth 1.
    #[test]
    fn the_point_test_agrees_with_the_division_polynomial_test() {
        let mut seen = [0usize; 2];
        for p in [101u64, 103, 107, 109, 113, 127, 131, 137, 139, 149] {
            for a in 0..10u64 {
                for b in 1..10u64 {
                    if (4 * a * a * a + 27 * b * b) % p == 0 {
                        continue;
                    }
                    let n = super::super::certify::count_prime(p, a, b);
                    let t = p as i64 + 1 - n as i64;
                    let delta = t * t - 4 * p as i64;
                    for ell in [3u64, 5, 7] {
                        let l = ell as i64;
                        // Depth exactly 1 at ℓ: ℓ² | Δ, and Δ/ℓ² not again
                        // a square multiple at ℓ.
                        if delta % (l * l) != 0 || (delta / (l * l)) % (l * l) == 0 {
                            continue;
                        }
                        let model = Model::Prime {
                            p: BigUint::from(p),
                            a: BigUint::from(a),
                            b: BigUint::from(b),
                        };
                        let tt = BigInt::from(t);
                        let poly = frobenius_is_scalar(&model, &tt, ell).unwrap();
                        let Some(points) = scalar_by_points(&model, &tt, ell) else {
                            continue;
                        };
                        assert_eq!(points, poly, "p={p} a={a} b={b} ℓ={ell}");
                        seen[usize::from(points)] += 1;
                    }
                }
            }
        }
        assert!(seen[0] > 0 && seen[1] > 0, "{seen:?}");
    }

    #[test]
    fn extension_fields_are_fields() {
        let f = Fp(BigUint::from(101u8));
        let q = BigUint::from(101u8);
        for k in [2usize, 3, 4, 6] {
            let m = irreducible(&f, &q, k, 7).expect("an irreducible polynomial");
            let g = Ext {
                base: &f,
                modulus: m,
            };
            let order = q.pow(k as u32);
            for i in 0..20 {
                let e = g.sample(i);
                if g.is_zero(&e) {
                    continue;
                }
                assert_eq!(g.mul(&e, &g.inv(&e)), g.one());
                let sq = g.mul(&e, &e);
                let r = sqrt(&g, &sq, &order).expect("a square has a root");
                assert_eq!(g.mul(&r, &r), sq);
            }
        }
    }
}
