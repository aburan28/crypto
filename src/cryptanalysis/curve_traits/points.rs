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
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};

use super::arith::lucas_v;
use super::torsion::{
    eigen_order, jac_add, jac_double, padd, pdivrem, pinv_mod, pmul, ppow, prem, psub,
    trim, Big2, Field, Fp, FpK, Jac, Poly, Small2,
};
use super::Model;
use crate::cryptanalysis::semaev_decomp::Gf2;

/// Random points tried before giving up.
pub const TRIES: u64 = 24;

/// The largest degree `k/2` of the twist [`scalar_by_points`] tries for
/// even `k = ord(λ)`.
pub const MAX_DEGREE: u64 = 160;

/// Work allowed for one test, in base-field word products: a guard on
/// build time.  The search stops once the scalar multiplications it has
/// made would exceed it, so the bound holds whatever the points do.
pub const POINT_TEST_COST: u64 = 20_000_000_000;

/// A test is not started unless the budget covers this many scalar
/// multiplications (a typical certificate needs two or three points).
const MIN_MULTIPLICATIONS: u64 = 8;

/// `GF(q^k) = F[w]/(m)` for a monic irreducible `m` of degree `k` over `F`
/// (`|F| = q`), with the table [`frobenius_table`] of `w^{qi}` for the
/// `q`-power map.
pub struct Ext<'a, F: Field> {
    pub base: &'a F,
    pub modulus: Poly<F::E>,
    frob: Vec<Poly<F::E>>,
}

impl<'a, F: Field> Ext<'a, F> {
    pub fn new(base: &'a F, modulus: Poly<F::E>, q: &BigUint) -> Self {
        let frob = frobenius_table(base, q, &modulus);
        Ext {
            base,
            modulus,
            frob,
        }
    }

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
    fn frobenius(&self, a: &Poly<F::E>) -> Poly<F::E> {
        apply_frobenius(self.base, &self.frob, a)
    }
}

/// `w^{qi} mod m` for `i < deg m`: the `q`-power map on `F[w]/(m)` is
/// `F`-linear (`c^q = c` on `F`), so this table is its matrix.
fn frobenius_table<F: Field>(f: &F, q: &BigUint, m: &Poly<F::E>) -> Vec<Poly<F::E>> {
    let k = m.len() - 1;
    let w_q = ppow(f, &vec![f.zero(), f.one()], q, m);
    let mut table = vec![prem(f, &vec![f.one()], m)];
    for i in 1..k {
        table.push(prem(f, &pmul(f, &table[i - 1], &w_q), m));
    }
    table
}

/// `a(w)^q = Σ a_i w^{qi}` from the table, reduced since every entry is.
fn apply_frobenius<F: Field>(f: &F, table: &[Poly<F::E>], a: &Poly<F::E>) -> Poly<F::E> {
    let mut acc = vec![f.zero(); table.len()];
    for (c, row) in a.iter().zip(table) {
        if f.is_zero(c) {
            continue;
        }
        for (slot, r) in acc.iter_mut().zip(row) {
            *slot = f.add(slot, &f.mul(c, r));
        }
    }
    trim(f, acc)
}

pub(crate) fn pgcd<F: Field>(f: &F, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
    let (mut a, mut b) = (a.clone(), b.clone());
    while !b.is_empty() {
        let (_, r) = pdivrem(f, &a, &b);
        (a, b) = (b, r);
    }
    a
}

/// A monic irreducible polynomial of degree `k ≥ 2` over `F` (`|F| = q`),
/// tried from seeded random coefficients and certified by Ben-Or's test:
/// `gcd(x^{q^i} − x, m) = 1` for every `i ≤ k/2`, since a reducible `m`
/// has a factor of degree `i ≤ k/2`, which divides `x^{q^i} − x`.  The
/// powers `x^{q^i}` come from the matrix of the `q`-power map.
pub fn irreducible<F: Field>(f: &F, q: &BigUint, k: usize, seed: u64) -> Option<Poly<F::E>> {
    let x = vec![f.zero(), f.one()];
    for attempt in 0..(64 * k as u64) {
        let mut m: Poly<F::E> = (0..k as u64)
            .map(|j| f.sample(seed.wrapping_add(attempt * (k as u64 + 1) + j)))
            .collect();
        m.push(f.one());
        let m = trim(f, m);
        let table = frobenius_table(f, q, &m);
        let mut power = prem(f, &x, &m);
        let mut coprime = true;
        for _ in 0..k / 2 {
            power = apply_frobenius(f, &table, &power);
            if pgcd(f, &m, &psub(f, &power, &x)).len() != 1 {
                coprime = false;
                break;
            }
        }
        if coprime {
            return Some(m);
        }
    }
    None
}

/// A square root of `a` in a field of odd order `order` (Tonelli–Shanks),
/// or `None` for a non-square.
pub fn sqrt<F: Field>(f: &F, a: &F::E, order: &BigUint) -> Option<F::E> {
    sqrt_with(f, a, order, None)
}

/// A non-square of the field of odd order `order`, found from seeds.
pub fn nonresidue<F: Field>(f: &F, order: &BigUint) -> Option<F::E> {
    let half = (order - 1u8) >> 1;
    let minus_one = f.neg(&f.one());
    (0u64..)
        .map(|i| f.sample(i ^ 0x5eed))
        .take(1024)
        .find(|z| !f.is_zero(z) && field_pow(f, z, &half) == minus_one)
}

fn field_pow<F: Field>(f: &F, b: &F::E, e: &BigUint) -> F::E {
    let mut acc = f.one();
    for i in (0..e.bits()).rev() {
        acc = f.mul(&acc, &acc);
        if e.bit(i) {
            acc = f.mul(&acc, b);
        }
    }
    acc
}

/// [`sqrt`] with the non-square `z` it needs supplied, when known.
pub fn sqrt_with<F: Field>(f: &F, a: &F::E, order: &BigUint, z: Option<&F::E>) -> Option<F::E> {
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
    let z = match z {
        Some(z) => z.clone(),
        None => (0u64..)
            .map(|i| f.sample(i ^ 0x5eed))
            .find(|z| !f.is_zero(z) && pow(z, &half) == minus_one)?,
    };
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
    fn neg(&self, p: &Self::P) -> Self::P;
    fn is_zero(&self, p: &Self::P) -> bool;
    /// Affine coordinates; `None` for the point at infinity.
    fn affine(&self, p: &Self::P) -> Option<Self::A>;
    /// Whether two points are equal.
    fn same(&self, p: &Self::P, q: &Self::P) -> bool {
        self.affine(p) == self.affine(q)
    }
    /// The `q`-power Frobenius `π`, `q` the order of the curve's field of
    /// definition.
    fn frobenius(&self, p: &Self::P) -> Self::P;
}

/// `y² = x³ + ax + b`, `|G| = order` odd, in Jacobian coordinates.
pub(crate) struct ShortW<'a, F: Field> {
    pub f: &'a F,
    pub a: F::E,
    pub b: F::E,
    pub order: BigUint,
    /// A non-square, found once for every square root.
    pub nonresidue: Option<F::E>,
}

impl<F: Field> CurveGroup for ShortW<'_, F> {
    type P = Jac<F::E>;
    type A = (F::E, F::E);
    fn sample_point(&self, i: u64) -> Option<Jac<F::E>> {
        let g = self.f;
        let x = g.sample(i);
        let rhs = g.add(&g.mul(&g.add(&g.mul(&x, &x), &self.a), &x), &self.b);
        let y = sqrt_with(g, &rhs, &self.order, self.nonresidue.as_ref())?;
        Some(Jac { x, y, z: g.one() })
    }
    fn mul(&self, p: &Jac<F::E>, n: &BigUint) -> Jac<F::E> {
        mul(self.f, &self.a, p, n)
    }
    fn add(&self, p: &Jac<F::E>, q: &Jac<F::E>) -> Jac<F::E> {
        jac_add(self.f, &self.a, p, q)
    }
    fn neg(&self, p: &Jac<F::E>) -> Jac<F::E> {
        Jac {
            x: p.x.clone(),
            y: self.f.neg(&p.y),
            z: p.z.clone(),
        }
    }
    fn is_zero(&self, p: &Jac<F::E>) -> bool {
        self.f.is_zero(&p.z)
    }
    fn affine(&self, p: &Jac<F::E>) -> Option<(F::E, F::E)> {
        affine(self.f, p)
    }
    /// `X₁Z₂² = X₂Z₁²` and `Y₁Z₂³ = Y₂Z₁³`, without inversions.
    fn same(&self, p: &Jac<F::E>, q: &Jac<F::E>) -> bool {
        let f = self.f;
        match (f.is_zero(&p.z), f.is_zero(&q.z)) {
            (true, true) => return true,
            (true, false) | (false, true) => return false,
            _ => {}
        }
        let (pz2, qz2) = (f.mul(&p.z, &p.z), f.mul(&q.z, &q.z));
        f.mul(&p.x, &qz2) == f.mul(&q.x, &pz2)
            && f.mul(&p.y, &f.mul(&qz2, &q.z)) == f.mul(&q.y, &f.mul(&pz2, &p.z))
    }
    fn frobenius(&self, p: &Jac<F::E>) -> Jac<F::E> {
        let f = self.f;
        Jac {
            x: f.frobenius(&p.x),
            y: f.frobenius(&p.y),
            z: f.frobenius(&p.z),
        }
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
    pub(crate) fn solve_artin_schreier(&self, c: &F::E) -> Option<F::E> {
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
    fn neg(&self, p: &Self::P) -> Self::P {
        p.as_ref().map(|(x, y)| (x.clone(), self.f.add(x, y)))
    }
    fn is_zero(&self, p: &Self::P) -> bool {
        p.is_none()
    }
    fn affine(&self, p: &Self::P) -> Option<(F::E, F::E)> {
        p.clone()
    }
    fn frobenius(&self, p: &Self::P) -> Self::P {
        p.as_ref()
            .map(|(x, y)| (self.f.frobenius(x), self.f.frobenius(y)))
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
        for p1 in &order_ell {
            let mut multiple = p1.clone();
            let mut inside = false;
            for _ in 1..ell {
                if g.same(&multiple, &s) {
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

/// `GF(p^m)` for a prime `p < 2^63` in word arithmetic ([`FpK`] on a
/// certified irreducible modulus), with the `p`-power map from the table
/// of `z^{pi}`: the extensions of prime-field curves.
struct PrimeExt {
    f: FpK,
    table: Vec<Vec<u64>>,
}

impl PrimeExt {
    /// `GF(p)` itself when `m = 1`.  The modulus is tried from seeded
    /// coefficients and certified by Ben-Or's test, as in [`irreducible`],
    /// in the word arithmetic: `gcd(z^{p^i} − z, m) = 1` for `i ≤ m/2`.
    fn new(p: &BigUint, m: usize) -> Option<Self> {
        let p64 = p.to_u64()?;
        let f = if m == 1 {
            FpK::new(p, &[BigUint::zero()])?
        } else {
            let mut found = None;
            for attempt in 0..(64 * m as u64) {
                let low: Vec<u64> = super::torsion::splitmix(0x1e1d ^ attempt.wrapping_mul(0x5851))
                    .take(m)
                    .map(|c| c % p64)
                    .collect();
                let f = FpK {
                    p: p64,
                    modulus: low,
                };
                let mut z = f.zero();
                z[1] = 1;
                let mut power = z.clone();
                let mut coprime = true;
                for _ in 0..m / 2 {
                    power = f.pow(&power, p);
                    if !f.coprime(&f.sub(&power, &z)) {
                        coprime = false;
                        break;
                    }
                }
                if coprime {
                    found = Some(f);
                    break;
                }
            }
            found?
        };
        let mut z = f.zero();
        if m > 1 {
            z[1] = 1;
        }
        let z_p = f.pow(&z, p);
        let mut table = vec![f.one()];
        for i in 1..m {
            table.push(f.mul(&table[i - 1], &z_p));
        }
        Some(PrimeExt { f, table })
    }
}

impl Field for PrimeExt {
    type E = Vec<u64>;
    fn zero(&self) -> Vec<u64> {
        self.f.zero()
    }
    fn one(&self) -> Vec<u64> {
        self.f.one()
    }
    fn int(&self, v: u64) -> Vec<u64> {
        self.f.int(v)
    }
    fn add(&self, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
        self.f.add(a, b)
    }
    fn neg(&self, a: &Vec<u64>) -> Vec<u64> {
        self.f.neg(a)
    }
    fn mul(&self, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
        self.f.mul(a, b)
    }
    fn inv(&self, a: &Vec<u64>) -> Vec<u64> {
        self.f.inv(a)
    }
    fn sample(&self, i: u64) -> Vec<u64> {
        self.f.sample(i)
    }
    fn frobenius(&self, a: &Vec<u64>) -> Vec<u64> {
        let p = u128::from(self.f.p);
        let mut acc = vec![0u128; a.len()];
        for (&c, row) in a.iter().zip(&self.table) {
            if c == 0 {
                continue;
            }
            for (slot, &r) in acc.iter_mut().zip(row) {
                *slot = (*slot + u128::from(c) * u128::from(r)) % p;
            }
        }
        acc.into_iter().map(|v| v as u64).collect()
    }
}

/// A test run on the group `E(GF(q^k))`, whatever its field arithmetic.
trait GroupTest {
    fn run<G: CurveGroup>(&self, g: &G) -> Option<bool>;
}

/// The rank test of [`scalar_by_points`].
struct RankTwo<'a> {
    n: &'a BigUint,
    ell: u64,
    budget: u64,
}

impl GroupTest for RankTwo<'_> {
    fn run<G: CurveGroup>(&self, g: &G) -> Option<bool> {
        rank_two(g, self.n, self.ell, self.budget)
    }
}

/// Work for one scalar multiplication by a number of `bits` bits over
/// `GF(q^k)`, in base-field word products: about 12 extension products
/// per bit, of `k²` base products, each `limbs²` words.  `None` for
/// models without field arithmetic.
pub(crate) fn multiplication_cost(model: &Model, k: u32, bits: u64) -> Option<u64> {
    let (limbs, base_k) = match model {
        Model::Prime { p, .. } => (p.bits().div_ceil(64), 1),
        Model::Extension { k: ext, .. } => (1, u64::from(*ext)),
        // GF(2^n), n ≤ 63, is one word; above that, a word per 64 bits.
        Model::Binary { n, .. } => (u64::from(*n).div_ceil(64), 1),
        Model::Other { .. } => return None,
    };
    let kk = u64::from(k) * base_k;
    Some(
        12u64
            .saturating_mul(bits)
            .saturating_mul(kk.saturating_mul(kk))
            .saturating_mul(limbs.saturating_mul(limbs))
            .max(1),
    )
}

/// What a word product costs in the extension [`on_extension`] builds,
/// relative to [`PrimeExt`]'s word arithmetic: the generic `F[w]/(m)` of
/// [`Ext`] allocates for every base operation (measured near 100×).
const GENERIC_OVERHEAD: u64 = 64;

/// [`multiplication_cost`] weighted by [`GENERIC_OVERHEAD`] where the
/// extension is the generic one.
fn weighted_cost(model: &Model, k: u32, bits: u64) -> Option<u64> {
    let words = matches!(model, Model::Prime { p, .. } if p.bits() < 63);
    let factor = if k == 1 || words { 1 } else { GENERIC_OVERHEAD };
    Some(multiplication_cost(model, k, bits)?.saturating_mul(factor))
}

/// `y² = x³ + ax + b` over the field `f` of order `order`, or its
/// quadratic twist `y² = x³ + ad²x + bd³` by a non-square `d`.
fn short_curve<'a, F: Field>(
    f: &'a F,
    a: F::E,
    b: F::E,
    order: BigUint,
    twist: bool,
) -> Option<ShortW<'a, F>> {
    let z = nonresidue(f, &order)?;
    let (a, b) = if twist {
        let z2 = f.mul(&z, &z);
        (f.mul(&a, &z2), f.mul(&b, &f.mul(&z2, &z)))
    } else {
        (a, b)
    };
    Some(ShortW {
        f,
        a,
        b,
        order,
        nonresidue: Some(z),
    })
}

/// `y² + xy = x³ + ax² + b` over the binary field `f` of degree `m`, or its
/// quadratic twist `y² + xy = x³ + (a + δ)x² + b`, `Tr(δ) = 1`.
fn binary_curve<'a, F: Field>(f: &'a F, a: F::E, b: F::E, m: u64, twist: bool) -> Option<BinaryW<'a, F>> {
    let mut g = BinaryW { f, a, b, m };
    if twist {
        let delta = (0u64..512)
            .map(|i| f.sample(i ^ 0x7a157))
            .find(|d| g.trace(d) == f.one())?;
        g.a = f.add(&g.a, &delta);
    }
    Some(g)
}

/// Run `test` on `E(GF(q^k))` for the curve of `model`, or on its
/// quadratic twist when `twist`, the extension built from a certified
/// irreducible modulus.  `None` for models without the arithmetic, or `ℓ`
/// the characteristic.
fn on_extension<T: GroupTest>(model: &Model, ell: u64, k: u32, twist: bool, test: &T) -> Option<bool> {
    fn over<F: Field, T: GroupTest>(
        f: &F,
        a: &F::E,
        b: &F::E,
        q: &BigUint,
        k: u32,
        twist: bool,
        test: &T,
    ) -> Option<bool> {
        let order = q.pow(k);
        if k == 1 {
            return test.run(&short_curve(f, a.clone(), b.clone(), order, twist)?);
        }
        let modulus = irreducible(f, q, k as usize, 0x1e1d)?;
        let e = Ext::new(f, modulus, q);
        test.run(&short_curve(&e, e.constant(a), e.constant(b), order, twist)?)
    }
    /// The binary analogue, `m` the degree of the base field over GF(2).
    #[allow(clippy::too_many_arguments)]
    fn over_binary<F: Field, T: GroupTest>(
        f: &F,
        a: &F::E,
        b: &F::E,
        q: &BigUint,
        m: u64,
        k: u32,
        twist: bool,
        test: &T,
    ) -> Option<bool> {
        if k == 1 {
            return test.run(&binary_curve(f, a.clone(), b.clone(), m, twist)?);
        }
        let modulus = irreducible(f, q, k as usize, 0x1e1d)?;
        let e = Ext::new(f, modulus, q);
        test.run(&binary_curve(
            &e,
            e.constant(a),
            e.constant(b),
            m * u64::from(k),
            twist,
        )?)
    }
    let q = model.q();
    match model {
        Model::Prime { p, a, b } => {
            if (p % ell).is_zero() || p <= &BigUint::from(3u8) {
                return None;
            }
            if let Some(f) = PrimeExt::new(p, k as usize) {
                // Word arithmetic, the extension built directly.
                let constant = |c: &BigUint| f.int(c.to_u64().expect("below p < 2^63"));
                return test.run(&short_curve(&f, constant(a), constant(b), q.pow(k), twist)?);
            }
            over(&Fp(p.clone()), a, b, &q, k, twist, test)
        }
        Model::Extension {
            p, modulus, a, b, ..
        } => {
            let f = FpK::new(p, modulus)?;
            if (BigUint::from(ell) % p).is_zero() {
                return None;
            }
            let (a, b) = (f.element(a)?, f.element(b)?);
            over(&f, &a, &b, &q, k, twist, test)
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
                over_binary(&f, &a, &b, &q, m, k, twist, test)
            } else {
                let f = Big2 {
                    n: *deg,
                    irr: irr.clone(),
                };
                over_binary(&f, a, b, &q, m, k, twist, test)
            }
        }
        Model::Other { .. } => None,
    }
}

/// `#E(GF(q^k)) = q^k + 1 − V_k(t, q)`, exact from the trace.
fn extension_order(q: &BigUint, trace: &BigInt, k: u32) -> Option<BigUint> {
    let qi = BigInt::from(q.clone());
    (qi.pow(k) + 1u8 - lucas_v(trace, &qi, k))
        .to_biguint()
        .filter(|n| !n.is_zero())
}

/// `v_ℓ(n)` and `n / ℓ^{v_ℓ(n)}`.
fn split_power(n: &BigUint, ell: u64) -> (u32, BigUint) {
    let l = BigUint::from(ell);
    let (mut v, mut cof) = (0u32, n.clone());
    while !cof.is_zero() && (&cof % &l).is_zero() {
        cof /= &l;
        v += 1;
    }
    (v, cof)
}

/// Whether Frobenius acts on `E[ℓ]` as a scalar, at depth 1 (`ℓ ∥ v`),
/// decided from points over `GF(q^k)` (see the module docs).  For even
/// `k` the group is the quadratic twist over `GF(q^{k/2})` instead: there
/// `λ^{k/2} = −1`, so the scalar is `π^{k/2} = −1` on `E[ℓ]`, all of `E[ℓ]`
/// in the twist's group (order `q^{k/2} + 1 + V_{k/2}(t, q)`), while a
/// Jordan block `λ(1 + N)` has `π^{k/2} + 1 = −(k/2)N` of rank 1 there
/// (`ℓ ∤ k/2`).  `None` when the model is not one this handles, the
/// degree passes [`MAX_DEGREE`], the work would exceed
/// [`POINT_TEST_COST`], or the search finds no certificate.
pub fn scalar_by_points(model: &Model, trace: &BigInt, ell: u64) -> Option<bool> {
    // ord(λ) ≤ ℓ − 1 is found by stepping, and the extension has that
    // degree: both bound ℓ.
    if ell > super::torsion::SCALAR_TEST_MAX_ELL {
        return None;
    }
    let k = eigen_order(ell, trace);
    if k == 0 {
        return None;
    }
    let q = model.q();
    let rank = |degree: u32, twist: bool, n: &BigUint, cost: u64, total: u64| {
        let budget = total / cost;
        if budget < MIN_MULTIPLICATIONS {
            return None;
        }
        on_extension(model, ell, degree, twist, &RankTwo { n, ell, budget })
    };
    // E(GF(q^k)) under the cost model it has always had, up to degree
    // [`EARLY_DEGREE`].
    if u64::from(k) <= EARLY_DEGREE {
        let n = extension_order(&q, trace, k)?;
        if let Some(r) = rank(k, false, &n, multiplication_cost(model, k, n.bits())?, POINT_TEST_COST) {
            return Some(r);
        }
    }
    // Then, for even k, the twist over GF(q^{k/2}), costed for the
    // arithmetic it gets.
    if k % 2 == 1 || u64::from(k / 2) > MAX_DEGREE {
        return None;
    }
    let degree = k / 2;
    let qi = BigInt::from(q.clone());
    let n = (qi.pow(degree) + 1u8 + lucas_v(trace, &qi, degree))
        .to_biguint()
        .filter(|n| !n.is_zero())?;
    // Each seeded point costs a square root, about three multiplications
    // the rank test does not count.
    let cost = weighted_cost(model, degree, n.bits())?.saturating_mul(4);
    rank(degree, true, &n, cost, NEW_POINT_COST)
}

/// The work allowed for one test of the twist route of
/// [`scalar_by_points`] and of [`scalar_on_power_by_points`].
pub const NEW_POINT_COST: u64 = 5_000_000_000;

/// The degree up to which [`scalar_by_points`] first tries `E(GF(q^k))`
/// itself.
const EARLY_DEGREE: u64 = 64;

/// The largest extension degree `m = ord(λ mod ℓ^j)` that
/// [`scalar_on_power_by_points`] builds.
pub const MAX_POWER_DEGREE: u64 = 512;

/// Whether Frobenius acts on `E[ℓ^j]` as the scalar `λ`, from points of
/// `E(GF(q^m))`, `m = ord(λ mod ℓ^j)`.
struct Eigen<'a> {
    n: &'a BigUint,
    ell: u64,
    j: u32,
    lambda: BigUint,
    budget: u64,
}

impl Eigen<'_> {
    /// `x` with `y = [x]b`, `b` of order `ℓ^n` (Pohlig–Hellman), or `None`
    /// when `y ∉ ⟨b⟩`.
    fn log<G: CurveGroup>(&self, g: &G, y: &G::P, b: &G::P, n: u32) -> Option<BigUint> {
        let l = BigUint::from(self.ell);
        if n == 0 {
            return g.is_zero(y).then(BigUint::zero);
        }
        let gamma = g.mul(b, &l.pow(n - 1));
        let mut x = BigUint::zero();
        let mut step = BigUint::one();
        for i in 0..n {
            let rest = g.add(y, &g.neg(&g.mul(b, &x)));
            let h = g.mul(&rest, &l.pow(n - 1 - i));
            let mut multiple = g.mul(&gamma, &BigUint::zero());
            let d = (0..self.ell).find(|_| {
                let hit = g.same(&multiple, &h);
                multiple = g.add(&multiple, &gamma);
                hit
            })?;
            x += &step * d;
            step *= &l;
        }
        g.same(&g.mul(b, &x), y).then_some(x)
    }

    /// From the point `top` of the largest order `ℓ^e1` seen and another
    /// `s`, two points of `E[ℓ^j]` that should generate it: `[ℓ^{e1−j}]top`
    /// and `[ℓ^{a−j}](s − [x]top)`, `a = v − e1`, `x` matching `s` on
    /// `[ℓ^a]`.  The caller checks whatever comes back.
    fn pair<G: CurveGroup>(
        &self,
        g: &G,
        top: &G::P,
        e1: u32,
        s: &G::P,
        v: u32,
    ) -> Option<(G::P, G::P)> {
        let l = BigUint::from(self.ell);
        let a = v.checked_sub(e1)?;
        if a < self.j {
            return None;
        }
        let y = g.mul(s, &l.pow(a));
        let b = g.mul(top, &l.pow(a));
        let x = self.log(g, &y, &b, e1.saturating_sub(a))?;
        let c = g.add(s, &g.neg(&g.mul(top, &x)));
        let t2 = g.mul(&c, &l.pow(a - self.j));
        let t1 = g.mul(top, &l.pow(e1 - self.j));
        Some((t1, t2))
    }
}

impl GroupTest for Eigen<'_> {
    fn run<G: CurveGroup>(&self, g: &G) -> Option<bool> {
        let l = BigUint::from(self.ell);
        let (v, cof) = split_power(self.n, self.ell);
        let j = self.j;
        let lj = l.pow(j);
        let lj1 = l.pow(j - 1);
        // Sylow projections seen, with their order exponents.
        let mut seen: Vec<(G::P, u32)> = Vec::new();
        let mut top: Option<usize> = None;
        let mut spent = 0u64;
        // TRIES points with a non-trivial ℓ-part, from at most 8·TRIES
        // seeds: over small fields the seeded x-coordinates repeat.
        for i in 0..8 * TRIES {
            if spent >= self.budget || seen.len() as u64 >= TRIES {
                return None;
            }
            spent += 1;
            // A square root, about three multiplications.
            spent += 2;
            let Some(p) = g.sample_point(i) else {
                continue;
            };
            spent += 1;
            let s = g.mul(&p, &cof);
            if g.is_zero(&s) {
                continue;
            }
            // The order ℓ^e of s.
            let mut e = 0u32;
            let mut cur = s.clone();
            loop {
                cur = g.mul(&cur, &l);
                e += 1;
                if g.is_zero(&cur) {
                    break;
                }
                if e == v {
                    // [ℓ^v]s ≠ O: N is not the group order; refuse.
                    return None;
                }
            }
            if e + j > v {
                // The ℓ-part is Z/ℓ^a × Z/ℓ^b with b ≥ e > v − j, so
                // a < j: E[ℓ^j] is not in E(GF(q^m)), which a scalar
                // would put there.
                return Some(false);
            }
            if e >= j {
                // A point of order ℓ^j: it must be an eigenvector.
                let t = g.mul(&s, &l.pow(e - j));
                spent += 1;
                if !g.same(&g.frobenius(&t), &g.mul(&t, &self.lambda)) {
                    return Some(false);
                }
            }
            seen.push((s, e));
            let k = seen.len() - 1;
            let fresh = match top {
                Some(t) if seen[t].1 >= e => vec![k],
                _ => {
                    top = Some(k);
                    (0..k).collect()
                }
            };
            let t = top.expect("set");
            if seen[t].1 < j {
                continue;
            }
            for other in fresh {
                spent += 2;
                let Some((t1, t2)) = self.pair(g, &seen[t].0, seen[t].1, &seen[other].0, v) else {
                    continue;
                };
                // Both in E[ℓ^j], both eigenvectors, and their order-ℓ
                // multiples independent: then ⟨t1, t2⟩ = E[ℓ^j] and π = λ
                // on it.
                if !g.is_zero(&g.mul(&t2, &lj)) {
                    continue;
                }
                if !g.same(&g.frobenius(&t2), &g.mul(&t2, &self.lambda)) {
                    return Some(false);
                }
                let (u1, u2) = (g.mul(&t1, &lj1), g.mul(&t2, &lj1));
                if g.is_zero(&u2) {
                    continue;
                }
                let mut multiple = u1.clone();
                let mut inside = false;
                for _ in 1..self.ell {
                    if g.same(&multiple, &u2) {
                        inside = true;
                        break;
                    }
                    multiple = g.add(&multiple, &u1);
                }
                if !inside {
                    return Some(true);
                }
            }
        }
        None
    }
}

/// Whether Frobenius acts on `E[ℓ^j]` as a scalar, `j ≥ 1`, decided from
/// points over `GF(q^m)` with a Frobenius check, for prime powers beyond
/// the division polynomials of [`super::torsion::scalar_on_power`].
///
/// A scalar `λ` has `2λ ≡ t (mod ℓ^j)`: one candidate for odd `ℓ`, two
/// for `ℓ = 2`.  If `π = λ` on `E[ℓ^j]` then `π^m = 1` there, `m` the
/// order of `λ` mod `ℓ^j`, so `E[ℓ^j] ⊆ E(GF(q^m))`.  With
/// `#E(GF(q^m)) = ℓ^v·c` exact from the trace, each candidate is refuted by
///
/// - `v < 2j`, or a point whose `ℓ`-part has order above `ℓ^{v−j}`: then
///   `E[ℓ^j] ⊄ E(GF(q^m))`;
/// - a point `T` of order `ℓ^j` with `πT ≠ [λ]T`, `π` applied to its
///   coordinates;
///
/// and confirmed by two points `T₁, T₂` of order `ℓ^j` with `πTᵢ = [λ]Tᵢ`
/// and `[ℓ^{j−1}]T₁, [ℓ^{j−1}]T₂` independent, which generate `E[ℓ^j]`.
/// Frobenius is a scalar when some candidate is confirmed and is not when
/// every one is refuted.  `None` otherwise: the search found neither, the
/// degree passed [`MAX_POWER_DEGREE`], or the work [`POINT_TEST_COST`].
pub fn scalar_on_power_by_points(
    model: &Model,
    trace: &BigInt,
    ell: u64,
    j: u32,
) -> Option<bool> {
    if j == 0 || ell > super::torsion::SCALAR_TEST_MAX_ELL {
        return None;
    }
    let n = BigUint::from(ell).pow(j);
    let tn = trace.mod_floor(&BigInt::from(n.clone())).to_biguint()?;
    let candidates: Vec<BigUint> = if ell == 2 {
        if tn.bit(0) {
            return None;
        }
        let half = &tn >> 1u8;
        vec![half.clone(), (half + (&n >> 1u8)) % &n]
    } else {
        // 2 is invertible mod n: λ = t·(n + 1)/2.
        vec![tn * ((&n + 1u8) >> 1u8) % &n]
    };
    let q = model.q();
    let mut verdict = Some(false);
    for lambda in candidates {
        if (&lambda % ell).is_zero() {
            // π is invertible on E[ℓ] (ℓ ∤ q), so not ≡ 0 mod ℓ.
            continue;
        }
        // m = ord(λ mod ℓ^j), by stepping.
        let mut m = 1u64;
        let mut x = lambda.clone();
        while !x.is_one() && m <= MAX_POWER_DEGREE {
            x = x * &lambda % &n;
            m += 1;
        }
        if m > MAX_POWER_DEGREE {
            verdict = None;
            continue;
        }
        let m = m as u32;
        let Some(order) = extension_order(&q, trace, m) else {
            verdict = None;
            continue;
        };
        if split_power(&order, ell).0 < 2 * j {
            continue;
        }
        let Some(cost) = weighted_cost(model, m, order.bits()) else {
            return None;
        };
        let budget = NEW_POINT_COST / cost;
        if budget < MIN_MULTIPLICATIONS {
            verdict = None;
            continue;
        }
        let test = Eigen {
            n: &order,
            ell,
            j,
            lambda,
            budget,
        };
        match on_extension(model, ell, m, false, &test) {
            Some(true) => return Some(true),
            Some(false) => {}
            None => verdict = None,
        }
    }
    verdict
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

    /// The prime-power point test agrees with the division polynomials of
    /// `torsion.rs` (themselves checked against brute-force torsion
    /// counts) on `E[ℓ^j]`, `ℓ^j ≤ 64`, including `ℓ = 2` with its two
    /// candidate scalars, and at `j = 1` with the scalar test on `E[ℓ]`.
    #[test]
    fn the_power_point_test_agrees_with_the_division_polynomials() {
        use crate::cryptanalysis::curve_traits::torsion::scalar_on_power;
        let mut seen = [0usize; 2];
        let primes = (101u64..700).filter(|&p| (2..p).take_while(|d| d * d <= p).all(|d| p % d != 0));
        for p in primes {
            for a in 0..6u64 {
                for b in 1..6u64 {
                    if (4 * a * a * a + 27 * b * b) % p == 0 {
                        continue;
                    }
                    let n = super::super::certify::count_prime(p, a, b);
                    let t = p as i64 + 1 - n as i64;
                    let delta = t * t - 4 * p as i64;
                    for (ell, j) in [(2u64, 1u32), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 2)] {
                        // Interesting only where Frobenius is a scalar on
                        // E[ℓ^{j−1}]: ℓ^{2(j−1)} | Δ at least.
                        let below = (ell as i64).pow(2 * (j - 1));
                        if delta % below != 0 || delta % (ell as i64).pow(2 * j) != 0 {
                            continue;
                        }
                        let model = Model::Prime {
                            p: BigUint::from(p),
                            a: BigUint::from(a),
                            b: BigUint::from(b),
                        };
                        let tt = BigInt::from(t);
                        let poly = if j == 1 {
                            frobenius_is_scalar(&model, &tt, ell)
                        } else {
                            scalar_on_power(&model, &tt, ell, j)
                        };
                        let Some(poly) = poly else {
                            continue;
                        };
                        let points = scalar_on_power_by_points(&model, &tt, ell, j)
                            .unwrap_or_else(|| panic!("undecided: p={p} a={a} b={b} ℓ^{j}={ell}^{j}"));
                        assert_eq!(points, poly, "p={p} a={a} b={b} ℓ^j={ell}^{j}");
                        seen[usize::from(points)] += 1;
                    }
                }
            }
        }
        assert!(seen[0] > 20 && seen[1] > 20, "{seen:?}");
    }

    /// Over GF(2^n) at `j = 1`, against the binary division polynomials.
    #[test]
    fn the_power_point_test_agrees_over_binary_fields() {
        use crate::binary_ecc::{F2mElement, IrreduciblePoly};
        use crate::cryptanalysis::semaev_decomp::Gf2;
        let mut seen = [0usize; 2];
        for (n, low) in [(7u32, vec![0u32, 1]), (8, vec![0, 2, 3, 4]), (9, vec![0, 4])] {
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
                        if delta % (ell * ell) as i64 != 0 {
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
                        let points = scalar_on_power_by_points(&model, &tt, ell, 1)
                            .unwrap_or_else(|| panic!("undecided: GF(2^{n}) a={a} b={b} ℓ={ell}"));
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
            let g = Ext::new(&f, m, &q);
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
