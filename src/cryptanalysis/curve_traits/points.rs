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
    trim, Big2, Field, Fp, FpK, Jac, Poly, Small2,
};
use super::Model;
use crate::cryptanalysis::semaev_decomp::Gf2;

/// Random points tried before giving up.
pub const TRIES: u64 = 24;

/// The largest extension degree `k = ord(λ)` tried.
pub const MAX_DEGREE: u64 = 64;

/// Work allowed for one test, in base-field word products: a guard on
/// build time.  The search stops once the scalar multiplications it has
/// made would exceed it, so the bound holds whatever the points do.
pub const POINT_TEST_COST: u64 = 20_000_000_000;

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

/// The group of an elliptic curve over a finite field, as the rank test
/// needs it: seeded points, multiples, sums, and a canonical form for
/// comparing points.
pub(crate) trait CurveGroup {
    type P: Clone;
    type A: PartialEq;
    /// A point at the `i`-th seeded x-coordinate, or `None` when that
    /// x-coordinate has none.
    fn sample_point(&self, i: u64) -> Option<Self::P>;
    fn mul(&self, p: &Self::P, n: &BigUint) -> Self::P;
    fn add(&self, p: &Self::P, q: &Self::P) -> Self::P;
    fn is_zero(&self, p: &Self::P) -> bool;
    /// Affine coordinates; `None` for the point at infinity.
    fn affine(&self, p: &Self::P) -> Option<Self::A>;
}

/// `y² = x³ + ax + b`, `|G| = order` odd, in Jacobian coordinates.
pub(crate) struct ShortW<'a, F: Field> {
    pub f: &'a F,
    pub a: F::E,
    pub b: F::E,
    pub order: BigUint,
}

impl<F: Field> CurveGroup for ShortW<'_, F> {
    type P = Jac<F::E>;
    type A = (F::E, F::E);
    fn sample_point(&self, i: u64) -> Option<Jac<F::E>> {
        let g = self.f;
        let x = g.sample(i);
        let rhs = g.add(&g.mul(&g.add(&g.mul(&x, &x), &self.a), &x), &self.b);
        let y = sqrt(g, &rhs, &self.order)?;
        Some(Jac { x, y, z: g.one() })
    }
    fn mul(&self, p: &Jac<F::E>, n: &BigUint) -> Jac<F::E> {
        mul(self.f, &self.a, p, n)
    }
    fn add(&self, p: &Jac<F::E>, q: &Jac<F::E>) -> Jac<F::E> {
        jac_add(self.f, &self.a, p, q)
    }
    fn is_zero(&self, p: &Jac<F::E>) -> bool {
        self.f.is_zero(&p.z)
    }
    fn affine(&self, p: &Jac<F::E>) -> Option<(F::E, F::E)> {
        affine(self.f, p)
    }
}

/// `y² + xy = x³ + ax² + b` over a field of characteristic 2 and degree
/// `m` over GF(2), in affine coordinates (`None` is the point at
/// infinity).
pub(crate) struct BinaryW<'a, F: Field> {
    pub f: &'a F,
    pub a: F::E,
    pub b: F::E,
    pub m: u64,
}

impl<F: Field> BinaryW<'_, F> {
    fn square(&self, x: &F::E) -> F::E {
        self.f.mul(x, x)
    }

    fn trace(&self, c: &F::E) -> F::E {
        let mut acc = c.clone();
        let mut pow = c.clone();
        for _ in 1..self.m {
            pow = self.square(&pow);
            acc = self.f.add(&acc, &pow);
        }
        acc
    }

    /// A root of `z² + z = c`, or `None` when `Tr(c) = 1`.  Odd degree:
    /// the half-trace.  Any degree: with `Tr(δ) = 1`,
    /// `z = Σ_{i<m−1} (Σ_{j>i} δ^{2^j}) c^{2^i}`.  The root is checked.
    fn solve_artin_schreier(&self, c: &F::E) -> Option<F::E> {
        let f = self.f;
        if !f.is_zero(&self.trace(c)) {
            return None;
        }
        let z = if self.m % 2 == 1 {
            let mut acc = c.clone();
            let mut pow = c.clone();
            for _ in 0..(self.m - 1) / 2 {
                pow = self.square(&self.square(&pow));
                acc = f.add(&acc, &pow);
            }
            acc
        } else {
            let delta = (0u64..)
                .map(|i| f.sample(i ^ 0xa57e))
                .take(256)
                .find(|d| self.trace(d) == f.one())?;
            // δ^{2^j} for j < m, then the suffix sums D_i = Σ_{j>i} δ^{2^j}.
            let mut dpow = Vec::with_capacity(self.m as usize);
            let mut d = delta;
            for _ in 0..self.m {
                dpow.push(d.clone());
                d = self.square(&d);
            }
            let mut suffix = f.zero();
            let mut acc = f.zero();
            // c^{2^i} for i from m − 2 down needs them in order: build up.
            let mut cpow = Vec::with_capacity(self.m as usize);
            let mut cc = c.clone();
            for _ in 0..self.m {
                cpow.push(cc.clone());
                cc = self.square(&cc);
            }
            for i in (0..self.m as usize - 1).rev() {
                suffix = f.add(&suffix, &dpow[i + 1]);
                acc = f.add(&acc, &f.mul(&suffix, &cpow[i]));
            }
            acc
        };
        (f.add(&self.square(&z), &z) == *c).then_some(z)
    }

    fn double(&self, p: &Option<(F::E, F::E)>) -> Option<(F::E, F::E)> {
        let f = self.f;
        let (x, y) = p.as_ref()?;
        if f.is_zero(x) {
            return None;
        }
        let lambda = f.add(x, &f.mul(y, &f.inv(x)));
        let x3 = f.add(&f.add(&self.square(&lambda), &lambda), &self.a);
        let y3 = f.add(&self.square(x), &f.mul(&f.add(&lambda, &f.one()), &x3));
        Some((x3, y3))
    }
}

impl<F: Field> CurveGroup for BinaryW<'_, F> {
    type P = Option<(F::E, F::E)>;
    type A = (F::E, F::E);
    fn sample_point(&self, i: u64) -> Option<Self::P> {
        let f = self.f;
        let x = f.sample(i);
        if f.is_zero(&x) {
            return None;
        }
        // y = xz: z² + z = x + a + b/x².
        let xi = f.inv(&x);
        let c = f.add(&f.add(&x, &self.a), &f.mul(&self.b, &self.square(&xi)));
        let z = self.solve_artin_schreier(&c)?;
        Some(Some((x.clone(), f.mul(&x, &z))))
    }
    fn mul(&self, p: &Self::P, n: &BigUint) -> Self::P {
        let mut acc: Self::P = None;
        for i in (0..n.bits()).rev() {
            acc = self.double(&acc);
            if n.bit(i) {
                acc = self.add(&acc, p);
            }
        }
        acc
    }
    fn add(&self, p: &Self::P, q: &Self::P) -> Self::P {
        let f = self.f;
        let (Some((x1, y1)), Some((x2, y2))) = (p, q) else {
            return if p.is_none() { q.clone() } else { p.clone() };
        };
        if x1 == x2 {
            // Q = P doubles; Q = −P = (x, x + y) cancels.
            return if y1 == y2 { self.double(p) } else { None };
        }
        let lambda = f.mul(&f.add(y1, y2), &f.inv(&f.add(x1, x2)));
        let x3 = f.add(
            &f.add(&f.add(&self.square(&lambda), &lambda), &f.add(x1, x2)),
            &self.a,
        );
        let y3 = f.add(&f.add(&f.mul(&lambda, &f.add(x1, &x3)), &x3), y1);
        Some((x3, y3))
    }
    fn is_zero(&self, p: &Self::P) -> bool {
        p.is_none()
    }
    fn affine(&self, p: &Self::P) -> Option<(F::E, F::E)> {
        p.clone()
    }
}

/// Decide the `ℓ`-rank of a curve group with `n = ℓ^v·m` points:
/// `Some(true)` for rank 2, `Some(false)` for rank 1.
fn rank_two<G: CurveGroup>(g: &G, n: &BigUint, ell: u64, budget: u64) -> Option<bool> {
    let l = BigUint::from(ell);
    // Scalar multiplications made so far, a point sample counting as one.
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
    let mut order_ell: Vec<G::P> = Vec::new();
    for i in 0..TRIES {
        if spent >= budget {
            return None;
        }
        spent += 1;
        let Some(p) = g.sample_point(i) else {
            continue;
        };
        spent += 1;
        let mut s = g.mul(&p, &cof);
        if g.is_zero(&s) {
            continue;
        }
        // s has order ℓ^e; step to the point of order ℓ below it.
        let mut e = 0u32;
        loop {
            // A multiplication by ℓ is a small fraction of one by N.
            let next = g.mul(&s, &l);
            e += 1;
            if g.is_zero(&next) {
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
        let t = g.affine(&s).expect("order ℓ, not O");
        for p1 in &order_ell {
            let mut multiple = p1.clone();
            let mut inside = false;
            for _ in 1..ell {
                if g.affine(&multiple).as_ref() == Some(&t) {
                    inside = true;
                    break;
                }
                multiple = g.add(&multiple, p1);
            }
            if !inside {
                return Some(true);
            }
        }
        order_ell.push(s);
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
        // GF(2^n), n ≤ 63, is one word; above that, a word per 64 bits.
        Model::Binary { n, .. } => (u64::from(*n).div_ceil(64), 1),
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
            let g = ShortW {
                f,
                a: a.clone(),
                b: b.clone(),
                order: order.clone(),
            };
            return rank_two(&g, n, ell, budget);
        }
        let modulus = irreducible(f, q, k as usize, 0x1e1d)?;
        let e = Ext { base: f, modulus };
        let g = ShortW {
            a: e.constant(a),
            b: e.constant(b),
            f: &e,
            order: order.clone(),
        };
        rank_two(&g, n, ell, budget)
    }
    /// The binary analogue, `m` the degree of the base field over GF(2).
    #[allow(clippy::too_many_arguments)]
    fn over_binary<F: Field>(
        f: &F,
        a: &F::E,
        b: &F::E,
        q: &BigUint,
        m: u64,
        k: u32,
        n: &BigUint,
        ell: u64,
        budget: u64,
    ) -> Option<bool> {
        if k == 1 {
            let g = BinaryW {
                f,
                a: a.clone(),
                b: b.clone(),
                m,
            };
            return rank_two(&g, n, ell, budget);
        }
        let modulus = irreducible(f, q, k as usize, 0x1e1d)?;
        let e = Ext { base: f, modulus };
        let g = BinaryW {
            a: e.constant(a),
            b: e.constant(b),
            f: &e,
            m: m * u64::from(k),
        };
        rank_two(&g, n, ell, budget)
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
        Model::Binary {
            n: deg, irr, a, b, ..
        } => {
            if ell == 2 {
                return None;
            }
            let m = u64::from(*deg);
            if *deg <= 63 {
                let f = Small2(Gf2::new(irr));
                let (a, b) = (f.0.from_element(a), f.0.from_element(b));
                over_binary(&f, &a, &b, &q, m, k, &n, ell, budget)
            } else {
                let f = Big2 {
                    n: *deg,
                    irr: irr.clone(),
                };
                over_binary(&f, a, b, &q, m, k, &n, ell, budget)
            }
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

    /// The same agreement over GF(2^n), through the binary group law and
    /// Artin–Schreier sampling, odd and even degree.
    #[test]
    fn the_binary_point_test_agrees_with_the_division_polynomial_test() {
        use crate::binary_ecc::{F2mElement, IrreduciblePoly};
        use crate::cryptanalysis::semaev_decomp::Gf2;
        let mut seen = [0usize; 2];
        for (n, low) in [
            (7u32, vec![0u32, 1]),
            (8, vec![0, 2, 3, 4]),
            (9, vec![0, 4]),
            (10, vec![0, 3]),
        ] {
            let irr = IrreduciblePoly {
                degree: n,
                low_terms: low,
            };
            let gf = Gf2::new(&irr);
            let q = 1i64 << n;
            for a in 0..2u64 {
                for b in 1..(q as u64).min(120) {
                    let count = super::super::certify::count_binary(&gf, a, b) as i64;
                    let t = q + 1 - count;
                    let delta = t * t - 4 * q;
                    for ell in [3u64, 5, 7] {
                        let l = ell as i64;
                        if delta % (l * l) != 0 {
                            continue;
                        }
                        let el = |v: u64| F2mElement::from_biguint(&BigUint::from(v), n);
                        let model = Model::Binary {
                            n,
                            modulus: BigUint::zero(),
                            irr: irr.clone(),
                            a: el(a),
                            b: el(b),
                        };
                        let tt = BigInt::from(t);
                        let Some(poly) = frobenius_is_scalar(&model, &tt, ell) else {
                            continue;
                        };
                        let Some(points) = scalar_by_points(&model, &tt, ell) else {
                            continue;
                        };
                        assert_eq!(points, poly, "GF(2^{n}) a={a} b={b} ℓ={ell}");
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
