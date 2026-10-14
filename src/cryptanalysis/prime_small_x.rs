//! # Small-x factor bases over prime fields, and their exact pair-sum yield.
//!
//! The binary-field work in this repository chooses a factor base by the
//! number of subgroup targets its sums cover.  This module asks the same
//! question of a prime-field curve `y² = x³ − 3x + b` with the classical
//! **small-x factor base**: the points whose abscissa is one of the smallest
//! integers `0, 1, 2, …` that lie on the curve, with both signs.
//!
//! For a base of `B` abscissae the sums of two base points are, up to sign,
//! `x(a+b)` and `x(a−b)` for the `B(B−1)/2` pairs of distinct abscissae and
//! `x(2a)` for each of the `B`, which is `B²` abscissae if no two coincide.
//! A target `T` is a signed sum of two base points exactly when `x(T)` is
//! one of them, and then so is `−T`.  The number of covered targets is
//! therefore twice the number of **distinct** sum abscissae `D`, and the
//! probability that a uniformly random group element decomposes in one trial
//! is exactly
//!
//! ```text
//!     yield = 2·D / (N − 1),        N = the group order,
//! ```
//!
//! with no sampling.  `D ≤ B²` with equality when no two sums collide, so
//! `2B²/(N − 1)` is the collision-free model and the gap between the two
//! counts the collisions.  [`measure`] computes `D` from every pair sum
//! (one batched inversion per row, about seven field multiplications a
//! pair), draws random group elements to check the count by Monte Carlo
//! where the yield is large enough to see, and plants known decompositions
//! to show the detector fires.
//!
//! The field is a trait so the same code runs on the native P-192 field
//! ([`super::p192_native`]) and on small toy primes, where collisions and
//! hits are common and the method can be tested against brute force.

use std::collections::HashMap;
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use rand::{rngs::StdRng, Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;

use super::p192_native::{
    p192_order, p192_prime, AffinePoint, Fp, JacobianPoint, P192Arithmetic, B_HEX, GX_HEX, GY_HEX,
};
use super::secg_standards::{SecgCurve, SecgField};
use crate::ct_bignum::{MontgomeryContext, Uint};

/// An exact, canonical key for a field element.
pub type Key = [u64; 3];

// ── Prime fields ───────────────────────────────────────────────────

/// What the yield measurement needs of a prime field.
pub trait PrimeField: Clone + Send + Sync {
    type E: Copy + Eq + Send + Sync;
    /// An exact, canonical, orderable key for an element.
    type K: Copy + Ord + std::hash::Hash + Send + Sync + std::fmt::Debug;
    fn modulus(&self) -> BigUint;
    fn zero(&self) -> Self::E;
    fn one(&self) -> Self::E;
    fn embed_u64(&self, v: u64) -> Self::E;
    fn embed_big(&self, v: &BigUint) -> Result<Self::E, String>;
    fn add(&mut self, a: Self::E, b: Self::E) -> Self::E;
    fn sub(&mut self, a: Self::E, b: Self::E) -> Self::E;
    fn mul(&mut self, a: Self::E, b: Self::E) -> Self::E;
    fn sqr(&mut self, a: Self::E) -> Self::E;
    fn inv(&mut self, a: Self::E) -> Option<Self::E>;
    /// Canonical value as three little-endian limbs.
    fn key(&self, a: Self::E) -> Self::K;

    fn neg(&mut self, a: Self::E) -> Self::E {
        let z = self.zero();
        self.sub(z, a)
    }

    fn is_zero(&self, a: Self::E) -> bool {
        a == self.zero()
    }

    /// `base^exp` by square and multiply.
    fn pow(&mut self, base: Self::E, exp: &BigUint) -> Self::E {
        let mut acc = self.one();
        for i in (0..exp.bits()).rev() {
            acc = self.sqr(acc);
            if exp.bit(i) {
                acc = self.mul(acc, base);
            }
        }
        acc
    }

    /// A uniform element, from a generator the caller seeds.
    fn random(&self, rng: &mut StdRng) -> Self::E {
        let m = self.modulus();
        let limbs = (m.bits() as usize).div_ceil(64) + 1;
        let bytes: Vec<u8> = (0..limbs * 8).map(|_| rng.gen()).collect();
        let v = BigUint::from_bytes_le(&bytes) % &m;
        self.embed_big(&v).expect("reduced value is canonical")
    }
}

/// The native P-192 field, three-limb Montgomery arithmetic.
#[derive(Clone)]
pub struct P192Field {
    ar: P192Arithmetic,
    p: BigUint,
}

impl P192Field {
    pub fn new() -> Self {
        Self {
            ar: P192Arithmetic::new(),
            p: p192_prime(),
        }
    }
}

impl Default for P192Field {
    fn default() -> Self {
        Self::new()
    }
}

impl PrimeField for P192Field {
    type E = Fp;
    type K = Key;
    fn modulus(&self) -> BigUint {
        self.p.clone()
    }
    fn zero(&self) -> Fp {
        self.ar.zero()
    }
    fn one(&self) -> Fp {
        self.ar.one()
    }
    fn embed_u64(&self, v: u64) -> Fp {
        self.ar.from_u64(v)
    }
    fn embed_big(&self, v: &BigUint) -> Result<Fp, String> {
        self.ar.from_biguint(v)
    }
    fn add(&mut self, a: Fp, b: Fp) -> Fp {
        self.ar.add(a, b)
    }
    fn sub(&mut self, a: Fp, b: Fp) -> Fp {
        self.ar.sub(a, b)
    }
    fn mul(&mut self, a: Fp, b: Fp) -> Fp {
        self.ar.mul(a, b)
    }
    fn sqr(&mut self, a: Fp) -> Fp {
        self.ar.square(a)
    }
    fn inv(&mut self, a: Fp) -> Option<Fp> {
        // Fermat with the cached exponent, not the module's per-call parse.
        if self.is_zero(a) {
            return None;
        }
        let e = &self.p - BigUint::from(2u8);
        Some(self.pow(a, &e))
    }
    fn key(&self, a: Fp) -> Key {
        self.ar.canonical(a).0
    }
}

/// A prime field of fewer than 63 bits, for toy curves.
#[derive(Clone)]
pub struct SmallPrime {
    p: u64,
}

impl SmallPrime {
    pub fn new(p: u64) -> Self {
        assert!(p > 3 && p < (1 << 62));
        Self { p }
    }
}

impl PrimeField for SmallPrime {
    type E = u64;
    type K = Key;
    fn modulus(&self) -> BigUint {
        BigUint::from(self.p)
    }
    fn zero(&self) -> u64 {
        0
    }
    fn one(&self) -> u64 {
        1
    }
    fn embed_u64(&self, v: u64) -> u64 {
        v % self.p
    }
    fn embed_big(&self, v: &BigUint) -> Result<u64, String> {
        v.to_u64()
            .filter(|&x| x < self.p)
            .ok_or_else(|| "element is not canonical".into())
    }
    fn add(&mut self, a: u64, b: u64) -> u64 {
        (a + b) % self.p
    }
    fn sub(&mut self, a: u64, b: u64) -> u64 {
        (a + self.p - b) % self.p
    }
    fn mul(&mut self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    fn sqr(&mut self, a: u64) -> u64 {
        self.mul(a, a)
    }
    fn inv(&mut self, a: u64) -> Option<u64> {
        if a == 0 {
            return None;
        }
        let (mut t, mut new_t) = (0i128, 1i128);
        let (mut r, mut new_r) = (self.p as i128, a as i128);
        while new_r != 0 {
            let q = r / new_r;
            (t, new_t) = (new_t, t - q * new_t);
            (r, new_r) = (new_r, r - q * new_r);
        }
        if r != 1 {
            return None;
        }
        Some(t.rem_euclid(self.p as i128) as u64)
    }
    fn key(&self, a: u64) -> Key {
        [a, 0, 0]
    }
}

/// An element of [`MontField`]: `L` little-endian limbs in Montgomery form.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MElem<const L: usize>(pub [u64; L]);

/// Any odd prime field of at most `64·L` bits, `L`-limb Montgomery
/// arithmetic: `L = 3` for the SECG `secp112…secp192` curves, 4 for
/// `secp224…secp256`, 6 for `secp384r1`, 9 for `secp521r1`.
#[derive(Clone)]
pub struct MontField<const L: usize> {
    ctx: MontgomeryContext<L>,
    p: BigUint,
    p_minus_2: BigUint,
}

impl<const L: usize> MontField<L> {
    pub fn new(p: &BigUint) -> Result<Self, String> {
        if p.bits() > 64 * L as u64 {
            return Err(format!("{} bits do not fit {L} limbs", p.bits()));
        }
        let ctx = MontgomeryContext::<L>::new(Uint::<L>::from_biguint(p))
            .ok_or("the modulus must be odd and nonzero")?;
        Ok(Self {
            ctx,
            p: p.clone(),
            p_minus_2: p - BigUint::from(2u8),
        })
    }

    fn n(&self) -> Uint<L> {
        self.ctx.n
    }
}

impl<const L: usize> PrimeField for MontField<L> {
    type E = MElem<L>;
    type K = [u64; L];
    fn modulus(&self) -> BigUint {
        self.p.clone()
    }
    fn zero(&self) -> MElem<L> {
        MElem([0; L])
    }
    fn one(&self) -> MElem<L> {
        MElem(self.ctx.r_mod_n.0)
    }
    fn embed_u64(&self, v: u64) -> MElem<L> {
        let mut a = [0u64; L];
        a[0] = v;
        MElem(self.ctx.to_montgomery(&Uint(a)).0)
    }
    fn embed_big(&self, v: &BigUint) -> Result<MElem<L>, String> {
        if v >= &self.p {
            return Err("element is not canonical".into());
        }
        Ok(MElem(self.ctx.to_montgomery(&Uint::<L>::from_biguint(v)).0))
    }
    fn add(&mut self, a: MElem<L>, b: MElem<L>) -> MElem<L> {
        MElem(Uint(a.0).add_mod(&Uint(b.0), &self.n()).0)
    }
    fn sub(&mut self, a: MElem<L>, b: MElem<L>) -> MElem<L> {
        MElem(Uint(a.0).sub_mod(&Uint(b.0), &self.n()).0)
    }
    fn mul(&mut self, a: MElem<L>, b: MElem<L>) -> MElem<L> {
        MElem(self.ctx.mont_mul(&Uint(a.0), &Uint(b.0)).0)
    }
    fn sqr(&mut self, a: MElem<L>) -> MElem<L> {
        MElem(self.ctx.mont_sqr(&Uint(a.0)).0)
    }
    fn inv(&mut self, a: MElem<L>) -> Option<MElem<L>> {
        if self.is_zero(a) {
            return None;
        }
        let e = self.p_minus_2.clone();
        Some(self.pow(a, &e))
    }
    fn key(&self, a: MElem<L>) -> [u64; L] {
        self.ctx.from_montgomery(&Uint(a.0)).0
    }
}

// ── Curves and points ──────────────────────────────────────────────

/// An affine point `(x, y)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Pt<E> {
    pub x: E,
    pub y: E,
}

/// `y² = x³ + a·x + b` over a prime field, with its group order `N`.
#[derive(Clone)]
pub struct Curve<F: PrimeField> {
    pub name: String,
    pub field: F,
    pub a: F::E,
    pub b: F::E,
    /// Order of the subgroup the base and targets live in: the whole group
    /// when the cofactor is 1, else the prime-order subgroup `[h]E`.
    pub order: BigUint,
    /// Whether that order is known to be prime.
    pub order_is_prime: bool,
    /// Cofactor `h`: the group has `h·order` points.  Base points enter as
    /// `[h]P`, targets as `[h]` of a uniform point, so both lie in the
    /// subgroup.
    pub cofactor: u64,
}

impl Curve<P192Field> {
    /// NIST P-192 (secp192r1): `a = −3`, prime order, cofactor 1.
    pub fn p192() -> Self {
        let mut field = P192Field::new();
        let a = field.neg(field.embed_u64(3));
        let b = field.ar.from_hex(B_HEX);
        Self {
            name: "secp192r1 (NIST P-192)".into(),
            field,
            a,
            b,
            order: p192_order(),
            order_is_prime: true,
            cofactor: 1,
        }
    }

    /// The curve's base point.
    pub fn p192_generator(&self) -> Pt<Fp> {
        Pt {
            x: self.field.ar.from_hex(GX_HEX),
            y: self.field.ar.from_hex(GY_HEX),
        }
    }
}

impl<const L: usize> Curve<MontField<L>> {
    /// A SECG prime-field curve of at most 192 bits from the standards
    /// database.  Needs `p ≡ 3 (mod 4)` (the square-root formula used).
    pub fn from_secg(c: &SecgCurve) -> Result<Self, String> {
        let SecgField::Prime { p } = &c.field else {
            return Err(format!("{} is not a prime-field curve", c.name));
        };
        let field = MontField::<L>::new(p)?;
        let curve = Self {
            name: c.name.clone(),
            a: field.embed_big(&c.a)?,
            b: field.embed_big(&c.b)?,
            order: c.order.clone(),
            order_is_prime: true,
            cofactor: c.cofactor,
            field,
        };
        let mut f = curve.field.clone();
        let g = curve.secg_generator(c)?;
        if !curve.on_curve(&mut f, &g) {
            return Err(format!("{}: the generator is not on the curve", c.name));
        }
        Ok(curve)
    }

    /// The published base point of a curve built by [`Curve::from_secg`].
    pub fn secg_generator(&self, c: &SecgCurve) -> Result<Pt<MElem<L>>, String> {
        Ok(Pt {
            x: self.field.embed_big(&c.gx)?,
            y: self.field.embed_big(&c.gy)?,
        })
    }
}

impl<F: PrimeField> Curve<F> {
    /// `[k]P` for `k ≥ 1`; `None` for the point at infinity.
    pub fn mul_big(&self, f: &mut F, p: &Pt<F::E>, k: &BigUint) -> Option<Pt<F::E>> {
        let mut acc: Option<Pt<F::E>> = None;
        for i in (0..k.bits()).rev() {
            acc = match acc {
                Some(a) => self.add_points(f, &a, &a),
                None => None,
            };
            if k.bit(i) {
                acc = match acc {
                    Some(a) => self.add_points(f, &a, p),
                    None => Some(*p),
                };
            }
        }
        acc
    }

    /// A uniform element of the subgroup of order [`Curve::order`] other
    /// than `O`: a uniform point pushed through `[h]`.  With cofactor 1 it
    /// is [`Curve::random_point`] exactly.
    pub fn random_subgroup_point(&self, f: &mut F, rng: &mut StdRng, exp: &BigUint) -> Pt<F::E> {
        if self.cofactor == 1 {
            return self.random_point(f, rng, exp);
        }
        let h = BigUint::from(self.cofactor);
        loop {
            let p = self.random_point(f, rng, exp);
            let p = if rng.gen::<bool>() {
                self.negate(f, &p)
            } else {
                p
            };
            if let Some(q) = self.mul_big(f, &p, &h) {
                return q;
            }
        }
    }

    /// `x³ + a·x + b`.
    pub fn rhs(&self, f: &mut F, x: F::E) -> F::E {
        let x2 = f.sqr(x);
        let t = f.add(x2, self.a);
        let t = f.mul(t, x);
        f.add(t, self.b)
    }

    fn euler_exponent(&self, f: &F) -> BigUint {
        (f.modulus() - BigUint::one()) >> 1
    }

    /// Whether `v` is a nonzero square.
    pub fn is_square(&self, f: &mut F, v: F::E, exp: &BigUint) -> bool {
        !f.is_zero(v) && f.pow(v, exp) == f.one()
    }

    /// A square root of `v`, or `None` when `v` is not a square: the
    /// `(p + 1)/4` power for `p ≡ 3 (mod 4)`, Tonelli–Shanks otherwise.
    pub fn sqrt(&self, f: &mut F, v: F::E) -> Option<F::E> {
        let p = f.modulus();
        if f.is_zero(v) {
            return Some(v);
        }
        if (&p % 4u8) == BigUint::from(3u8) {
            let e = (p + BigUint::one()) >> 2;
            let r = f.pow(v, &e);
            return (f.sqr(r) == v).then_some(r);
        }
        let one = f.one();
        let pm1 = &p - BigUint::one();
        if f.pow(v, &(&pm1 >> 1)) != one {
            return None;
        }
        // p − 1 = q·2^s with q odd.
        let s = pm1.trailing_zeros().unwrap_or(0) as usize;
        let q = &pm1 >> s;
        // A quadratic non-residue.
        let mut z = 2u64;
        let nonres = loop {
            let zz = f.embed_u64(z);
            if f.pow(zz, &(&pm1 >> 1)) != one {
                break zz;
            }
            z += 1;
        };
        let mut m = s;
        let mut c = f.pow(nonres, &q);
        let mut t = f.pow(v, &q);
        let mut r = f.pow(v, &((&q + BigUint::one()) >> 1));
        while t != one {
            let mut i = 0usize;
            let mut t2 = t;
            while t2 != one {
                t2 = f.sqr(t2);
                i += 1;
                if i == m {
                    return None;
                }
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = f.sqr(b);
            }
            m = i;
            c = f.sqr(b);
            t = f.mul(t, c);
            r = f.mul(r, b);
        }
        (f.sqr(r) == v).then_some(r)
    }

    pub fn on_curve(&self, f: &mut F, p: &Pt<F::E>) -> bool {
        let lhs = f.sqr(p.y);
        lhs == self.rhs(f, p.x)
    }

    /// `P + Q`, or `None` for the point at infinity.
    pub fn add_points(&self, f: &mut F, p: &Pt<F::E>, q: &Pt<F::E>) -> Option<Pt<F::E>> {
        let lambda = if p.x == q.x {
            let sum_y = f.add(p.y, q.y);
            if f.is_zero(sum_y) {
                return None;
            }
            let x2 = f.sqr(p.x);
            let three = f.embed_u64(3);
            let num = f.mul(x2, three);
            let num = f.add(num, self.a);
            let two = f.embed_u64(2);
            let den = f.mul(two, p.y);
            let inv = f.inv(den)?;
            f.mul(num, inv)
        } else {
            let num = f.sub(q.y, p.y);
            let den = f.sub(q.x, p.x);
            let inv = f.inv(den)?;
            f.mul(num, inv)
        };
        let l2 = f.sqr(lambda);
        let x3 = f.sub(l2, p.x);
        let x3 = f.sub(x3, q.x);
        let dx = f.sub(p.x, x3);
        let y3 = f.mul(lambda, dx);
        let y3 = f.sub(y3, p.y);
        Some(Pt { x: x3, y: y3 })
    }

    pub fn negate(&self, f: &mut F, p: &Pt<F::E>) -> Pt<F::E> {
        Pt {
            x: p.x,
            y: f.neg(p.y),
        }
    }

    /// A uniform group element other than `O`, from `rng`.
    pub fn random_point(&self, f: &mut F, rng: &mut StdRng, exp: &BigUint) -> Pt<F::E> {
        loop {
            let x = f.random(rng);
            let r = self.rhs(f, x);
            if self.is_square(f, r, exp) {
                let y = self.sqrt(f, r).expect("square has a root");
                return Pt { x, y };
            }
        }
    }
}

/// The first `size` points with the smallest abscissae `x = 0, 1, 2, …`
/// that lie on the curve (one of the two signs of each).  With a cofactor
/// `h > 1` each enters as `[h]P`, so the base lies in the prime-order
/// subgroup; infinity and repeated abscissae are skipped.
pub fn small_x_base<F: PrimeField>(c: &Curve<F>, size: usize) -> Result<Vec<Pt<F::E>>, String> {
    small_x_base_with_extent(c, size).map(|(base, _)| base)
}

/// [`small_x_base`] and the largest `x` examined.
pub fn small_x_base_with_extent<F: PrimeField>(
    c: &Curve<F>,
    size: usize,
) -> Result<(Vec<Pt<F::E>>, u64), String> {
    let mut f = c.field.clone();
    let exp = c.euler_exponent(&f);
    let h = BigUint::from(c.cofactor);
    let mut out: Vec<Pt<F::E>> = Vec::with_capacity(size);
    let mut seen: std::collections::HashSet<F::K> = std::collections::HashSet::new();
    let mut x = 0u64;
    let mut largest = 0u64;
    while out.len() < size {
        if x > 64 * size as u64 + 1024 {
            return Err("ran out of small abscissae on this curve".into());
        }
        let xe = f.embed_u64(x);
        let r = c.rhs(&mut f, xe);
        if c.is_square(&mut f, r, &exp) {
            let y = c.sqrt(&mut f, r).ok_or("square without a root")?;
            let p = Pt { x: xe, y };
            largest = x;
            let q = if c.cofactor == 1 {
                Some(p)
            } else {
                c.mul_big(&mut f, &p, &h)
            };
            if let Some(q) = q {
                if seen.insert(f.key(q.x)) {
                    out.push(q);
                }
            }
        }
        x += 1;
    }
    Ok((out, largest))
}

/// The abscissae of `a ± b` for every pair of distinct base abscissae and of
/// `2a` for every `a`: `B²` keys with repeats where sums collide.
pub fn pair_sum_keys<F: PrimeField>(c: &Curve<F>, base: &[Pt<F::E>]) -> Vec<F::K> {
    let n = base.len();
    let mut keys: Vec<F::K> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut f = c.field.clone();
            let mut out: Vec<F::K> = Vec::with_capacity(2 * (n - i) + 1);
            let a = base[i];
            // Doubling: x(2a) = λ² − 2x with λ = (3x² + a)/(2y).
            let x2 = f.sqr(a.x);
            let three = f.embed_u64(3);
            let num = f.mul(x2, three);
            let num = f.add(num, c.a);
            let two = f.embed_u64(2);
            let den = f.mul(two, a.y);
            let inv = f.inv(den).expect("an on-curve point has y != 0 here");
            let l = f.mul(num, inv);
            let l2 = f.sqr(l);
            let t = f.sub(l2, a.x);
            let x3 = f.sub(t, a.x);
            out.push(f.key(x3));
            // a ± b for j > i, with one inversion for the whole row.
            let rest = &base[i + 1..];
            if !rest.is_empty() {
                let dens: Vec<F::E> = rest.iter().map(|b| f.sub(a.x, b.x)).collect();
                let mut prefix = Vec::with_capacity(dens.len());
                let mut acc = f.one();
                for d in &dens {
                    acc = f.mul(acc, *d);
                    prefix.push(acc);
                }
                let mut inv_acc = f.inv(acc).expect("abscissae are distinct");
                let mut invs = vec![f.zero(); dens.len()];
                for k in (0..dens.len()).rev() {
                    let before = if k == 0 { f.one() } else { prefix[k - 1] };
                    invs[k] = f.mul(inv_acc, before);
                    inv_acc = f.mul(inv_acc, dens[k]);
                }
                for (k, b) in rest.iter().enumerate() {
                    let sx = f.add(a.x, b.x);
                    for sign_plus in [false, true] {
                        let ny = if sign_plus {
                            f.add(a.y, b.y)
                        } else {
                            f.sub(a.y, b.y)
                        };
                        let lam = f.mul(ny, invs[k]);
                        let l2 = f.sqr(lam);
                        let x3 = f.sub(l2, sx);
                        out.push(f.key(x3));
                    }
                }
            }
            out
        })
        .flatten()
        .collect();
    keys.par_sort_unstable();
    // Keep the repeats: the caller counts them, then dedups.
    keys
}

// ── Measurement ────────────────────────────────────────────────────

/// What to measure.
#[derive(Clone, Debug, Serialize)]
pub struct Options {
    /// Abscissae in the base; the base has twice as many points.
    pub size: usize,
    /// Random group elements tested by Monte Carlo.
    pub trials: u64,
    /// Planted decompositions that must be detected and recovered.
    pub controls: usize,
    pub seed: u64,
}

/// How much of the raw small-x base's pair-sum yield survives in the
/// prime-order subgroup, by sampling (see [`raw_subgroup_fraction`]).
#[derive(Clone, Debug, Serialize)]
pub struct RawBase {
    pub samples: u64,
    pub in_subgroup: u64,
    /// `in_subgroup / samples`: the share of signed pair sums of raw small-x
    /// points that lie in the subgroup.
    pub fraction: f64,
    pub ci95: [f64; 2],
    /// `1/h`, the share if the cofactor classes were uniform.
    pub one_over_h: f64,
    /// `log₂` of `fraction · 2B²/(n − 1)`, the raw-base yield when no sums collide.
    pub raw_yield_log2: f64,
}

/// The measured cost of deciding whether a target is a signed sum of two
/// base points by solving Semaev's summation polynomial `S₃(x₁, x₂, X) = 0`
/// for `x₂` at every base abscissa `x₁`, against the pair-sum table lookup,
/// the unit of work that the yield counts. See `semaev_s3_cost`.
#[derive(Clone, Debug, Serialize)]
pub struct S3Cost {
    pub base_abscissae: usize,
    /// Random subgroup targets decided by both methods.
    pub trials: u64,
    /// Wall microseconds per target for the S₃ solve, one thread.
    pub s3_us_per_trial: f64,
    /// Wall nanoseconds per target for the sorted pair-sum lookup.
    pub table_ns_per_trial: f64,
    /// Nanoseconds per field multiplication, one thread.
    pub field_mul_ns: f64,
    /// S₃ cost in field multiplications per target.
    pub s3_field_muls_per_trial: f64,
    /// Microseconds per affine group addition, one thread.
    pub group_add_us: f64,
    /// Targets on which the two methods disagree (must be 0).
    pub disagreements: u64,
    /// Planted decompositions the S₃ solve found, of `planted`.
    pub planted_found: usize,
    pub planted: usize,
    /// `log₂` seconds of S₃ solving to collect `B + 1` relations at the exact
    /// yield: trials for a full system times the cost of one solve.
    pub s3_full_system_log2_s: f64,
    /// `log₂` seconds of Pollard rho at `√n` affine additions.
    pub rho_log2_s: f64,
    /// `s3_full_system_log2_s − rho_log2_s`: the factor by which this
    /// method is slower than rho, in bits, on the same machine.
    pub bits_slower_than_rho: f64,
}

/// The exact yield of a small-x base and its checks.
#[derive(Clone, Debug, Serialize)]
pub struct YieldReport {
    pub curve: String,
    /// `"prime"` or `"binary"`.
    pub field_kind: &'static str,
    /// Cofactor `h`; the base and targets lie in the subgroup of order
    /// [`YieldReport::order`], reached through `[h]`.
    pub cofactor: u64,
    pub field_bits: u64,
    pub order_bits: f64,
    pub order: String,
    pub base_abscissae: usize,
    pub base_points: usize,
    /// The largest `x` used by the base.
    pub largest_small_x: u64,
    /// `B²`: pair-sum abscissae if no two coincide.
    pub sums: u64,
    pub distinct_sums: u64,
    pub collisions: u64,
    /// `2·D`.
    pub covered_targets: String,
    /// `log₂` of the exact yield per trial, `2D/(N − 1)`.
    pub exact_yield_log2: f64,
    pub exact_yield: f64,
    /// `log₂` of the collision-free model, `2B²/(N − 1)`.
    pub model_yield_log2: f64,
    pub model_yield: f64,
    /// `log₂` of the expected trials per relation, `(N − 1)/(2D)`.
    pub expected_trials_log2: f64,
    /// `log₂` of the trials to collect `B + 1` relations.
    pub trials_for_full_system_log2: f64,
    /// `log₂ √N`, the cost of Pollard rho in group operations to leading order.
    pub rho_log2: f64,
    pub mc_trials: u64,
    pub mc_hits: u64,
    pub mc_hits_verified: u64,
    pub mc_rate: f64,
    pub mc_ci95: [f64; 2],
    /// Monte Carlo rate over the exact yield, when there are hits.
    pub mc_over_exact: Option<f64>,
    pub controls: usize,
    pub controls_detected: usize,
    pub controls_recovered: usize,
    pub timing_ms: Timing,
    pub threads: usize,
    /// Present when the raw small-x variant was sampled.
    pub raw_base: Option<RawBase>,
    /// Present when the cost of deciding a target by a Semaev S₃ solve was measured.
    pub s3: Option<S3Cost>,
}

#[derive(Clone, Debug, Serialize)]
pub struct Timing {
    pub base: f64,
    pub pair_sums: f64,
    pub monte_carlo: f64,
    /// Wall microseconds per Monte Carlo trial, all threads together.
    pub us_per_trial_wall: f64,
    pub controls: f64,
}

pub(crate) fn log2_big(x: &BigUint) -> f64 {
    if x.is_zero() {
        return f64::NEG_INFINITY;
    }
    let b = x.bits();
    if b <= 64 {
        return x.to_f64().unwrap_or(f64::MAX).log2();
    }
    let top = (x >> (b - 64) as usize).to_f64().unwrap_or(f64::MAX);
    top.log2() + (b - 64) as f64
}

/// 95% interval for a binomial rate; exact one-sided bound when no hit.
pub(crate) fn binomial_ci(hits: u64, trials: u64) -> [f64; 2] {
    if trials == 0 {
        return [0.0, 1.0];
    }
    let n = trials as f64;
    if hits == 0 {
        return [0.0, 1.0 - 0.025f64.powf(1.0 / n)];
    }
    let p = hits as f64 / n;
    let z = 1.959_964;
    let denom = 1.0 + z * z / n;
    let centre = (p + z * z / (2.0 * n)) / denom;
    let half = z * (p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt() / denom;
    [(centre - half).max(0.0), (centre + half).min(1.0)]
}

/// Whether `target` is `±a ± b` for base points; returns the signs and
/// indices when it is.
fn recover<F: PrimeField>(
    c: &Curve<F>,
    f: &mut F,
    base: &[Pt<F::E>],
    index: &HashMap<F::K, usize>,
    target: &Pt<F::E>,
) -> Option<(bool, usize, usize)> {
    for (i, a) in base.iter().enumerate() {
        for negate in [false, true] {
            let a_signed = if negate { c.negate(f, a) } else { *a };
            // Q = target − a_signed, so target = a_signed + Q.
            let minus_a = c.negate(f, &a_signed);
            let Some(q) = c.add_points(f, target, &minus_a) else {
                continue;
            };
            if let Some(&j) = index.get(&f.key(q.x)) {
                return Some((negate, i, j));
            }
        }
    }
    None
}

/// Measure the exact yield of the small-x base of `c`.
pub fn measure<F: PrimeField>(c: &Curve<F>, opts: &Options) -> Result<YieldReport, String> {
    if opts.size < 2 {
        return Err("a base needs at least two abscissae".into());
    }
    let mut f = c.field.clone();
    let exp = c.euler_exponent(&f);

    let t0 = Instant::now();
    let (base, largest_x) = small_x_base_with_extent(c, opts.size)?;
    for p in &base {
        if !c.on_curve(&mut f, p) {
            return Err("a base point is not on the curve".into());
        }
    }
    let base_ms = t0.elapsed().as_secs_f64() * 1000.0;

    check_memory(opts.size, std::mem::size_of::<F::K>())?;
    let t1 = Instant::now();
    let mut keys = pair_sum_keys(c, &base);
    let total = keys.len() as u64;
    keys.dedup();
    let distinct = keys.len() as u64;
    let pair_ms = t1.elapsed().as_secs_f64() * 1000.0;

    let order_minus_one = &c.order - BigUint::one();
    let covered = BigUint::from(2 * distinct);
    let exact_log2 = log2_big(&covered) - log2_big(&order_minus_one);
    let model_log2 =
        log2_big(&BigUint::from(2 * (opts.size as u64).pow(2))) - log2_big(&order_minus_one);
    let exact_yield = 2.0f64.powf(exact_log2);
    let model_yield = 2.0f64.powf(model_log2);

    // Index of base abscissae, for recovering decompositions.
    let index: HashMap<F::K, usize> = base
        .iter()
        .enumerate()
        .map(|(i, p)| (f.key(p.x), i))
        .collect();

    // Monte Carlo on uniform group elements.
    let t2 = Instant::now();
    let chunk = 2048u64;
    let chunks = opts.trials.div_ceil(chunk);
    let (hits, samples) = (0..chunks)
        .into_par_iter()
        .map(|ci| {
            let mut f = c.field.clone();
            let mut rng =
                StdRng::seed_from_u64(opts.seed ^ 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(ci + 1));
            let here = chunk.min(opts.trials - ci * chunk);
            let mut hits = 0u64;
            let mut samples: Vec<Pt<F::E>> = Vec::new();
            for _ in 0..here {
                let p = c.random_subgroup_point(&mut f, &mut rng, &exp);
                if keys.binary_search(&f.key(p.x)).is_ok() {
                    hits += 1;
                    if samples.len() < 4 {
                        samples.push(p);
                    }
                }
            }
            (hits, samples)
        })
        .reduce(
            || (0u64, Vec::new()),
            |(h1, mut s1), (h2, s2)| {
                s1.extend(s2);
                s1.truncate(16);
                (h1 + h2, s1)
            },
        );
    let mc_ms = t2.elapsed().as_secs_f64() * 1000.0;

    // A hit is a real decomposition: recover it.
    let mut verified = 0u64;
    for p in samples.iter().take(8) {
        if recover(c, &mut f, &base, &index, p).is_some() {
            verified += 1;
        }
    }

    // Positive control: plant a decomposition and require the detector fire.
    let t3 = Instant::now();
    let mut rng = StdRng::seed_from_u64(opts.seed ^ 0x434f_4e54_524f_4c53);
    let (mut detected, mut recovered) = (0usize, 0usize);
    let mut planted = 0usize;
    while planted < opts.controls {
        let i = rng.gen_range(0..base.len());
        let j = rng.gen_range(0..base.len());
        let a = base[i];
        let mut b = base[j];
        if rng.gen::<bool>() {
            b = c.negate(&mut f, &b);
        }
        let Some(target) = c.add_points(&mut f, &a, &b) else {
            continue;
        };
        planted += 1;
        if keys.binary_search(&f.key(target.x)).is_ok() {
            detected += 1;
        }
        if recover(c, &mut f, &base, &index, &target).is_some() {
            recovered += 1;
        }
    }
    let ctl_ms = t3.elapsed().as_secs_f64() * 1000.0;

    let trials = opts.trials;
    let order_bits = log2_big(&c.order);
    let trials_full = (opts.size as f64 + 1.0).log2() - exact_log2;
    Ok(YieldReport {
        curve: c.name.clone(),
        field_kind: "prime",
        cofactor: c.cofactor,
        field_bits: f.modulus().bits(),
        order_bits,
        order: c.order.to_string(),
        base_abscissae: opts.size,
        base_points: 2 * opts.size,
        largest_small_x: largest_x,
        sums: total,
        distinct_sums: distinct,
        collisions: total - distinct,
        covered_targets: covered.to_string(),
        exact_yield_log2: exact_log2,
        exact_yield,
        model_yield_log2: model_log2,
        model_yield,
        expected_trials_log2: -exact_log2,
        trials_for_full_system_log2: trials_full,
        rho_log2: order_bits / 2.0,
        mc_trials: trials,
        mc_hits: hits,
        mc_hits_verified: verified,
        mc_rate: if trials == 0 {
            0.0
        } else {
            hits as f64 / trials as f64
        },
        mc_ci95: binomial_ci(hits, trials),
        mc_over_exact: (hits > 0 && exact_yield > 0.0)
            .then(|| (hits as f64 / trials as f64) / exact_yield),
        controls: opts.controls,
        controls_detected: detected,
        controls_recovered: recovered,
        timing_ms: Timing {
            base: base_ms,
            pair_sums: pair_ms,
            monte_carlo: mc_ms,
            us_per_trial_wall: if trials == 0 {
                0.0
            } else {
                mc_ms * 1000.0 / trials as f64
            },
            controls: ctl_ms,
        },
        threads: rayon::current_num_threads(),
        raw_base: None,
        s3: None,
    })
}

/// Refuse a base whose pair-sum keys would not fit in memory: the keys are
/// held twice at the peak (per-row vectors, then the flattened set).
pub(crate) fn check_memory(size: usize, key_bytes: usize) -> Result<(), String> {
    let bytes = 2.0 * (size as f64).powi(2) * key_bytes as f64;
    if bytes > 9.0e9 {
        return Err(format!(
            "a base of {size} abscissae needs about {:.1} GB of pair-sum keys at this field size; use fewer",
            bytes / 1e9
        ));
    }
    Ok(())
}

/// The share of signed pair sums `±a ± b` of *raw* small-x points (no
/// cofactor multiplication) that lie in the prime-order subgroup, from
/// `samples` random pairs: `[n]S = O` is the membership test. With the
/// raw base the yield over the subgroup is this fraction times the
/// subgroup-base yield, since a target of the subgroup can only be a sum
/// that lands in it.
pub fn raw_subgroup_fraction<F: PrimeField>(
    c: &Curve<F>,
    size: usize,
    samples: u64,
    seed: u64,
) -> Result<RawBase, String> {
    let mut raw_curve = c.clone();
    raw_curve.cofactor = 1;
    let base = small_x_base(&raw_curve, size)?;
    let hits: u64 = (0..samples)
        .into_par_iter()
        .map(|i| {
            let mut f = c.field.clone();
            let mut rng = StdRng::seed_from_u64(seed ^ 0x5241_5742_4153_45u64.wrapping_mul(i + 1));
            let a = base[rng.gen_range(0..base.len())];
            let mut b = base[rng.gen_range(0..base.len())];
            if rng.gen::<bool>() {
                b = c.negate(&mut f, &b);
            }
            match c.add_points(&mut f, &a, &b) {
                None => 1, // a + b = O is in the subgroup
                Some(sum) => u64::from(c.mul_big(&mut f, &sum, &c.order).is_none()),
            }
        })
        .sum();
    Ok(raw_base_summary(c.cofactor, &c.order, size, hits, samples))
}

pub(crate) fn raw_base_summary(
    cofactor: u64,
    order: &BigUint,
    size: usize,
    hits: u64,
    samples: u64,
) -> RawBase {
    let fraction = hits as f64 / samples.max(1) as f64;
    let order_minus_one = order - BigUint::one();
    let sub_log2 = log2_big(&BigUint::from(2 * (size as u64).pow(2))) - log2_big(&order_minus_one);
    RawBase {
        samples,
        in_subgroup: hits,
        fraction,
        ci95: binomial_ci(hits, samples),
        one_over_h: 1.0 / cofactor as f64,
        raw_yield_log2: fraction.log2() + sub_log2,
    }
}

// ── Toy curves ─────────────────────────────────────────────────────

fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    let mulmod = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    let powmod = |mut b: u64, mut e: u64| {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = mulmod(r, b);
            }
            b = mulmod(b, b);
            e >>= 1;
        }
        r
    };
    'witness: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..s - 1 {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// A toy curve `y² = x³ − 3x + b` over a prime `p ≡ 3 (mod 4)` just below
/// `2^bits`, with `b` chosen from `seed` so the cubic has no root (no
/// rational 2-torsion, so the group order is odd and every abscissa on the
/// curve carries two points).  The order is counted exactly.  Needs
/// `bits ≤ 28`: the count walks the field.
pub fn toy_curve(bits: u32, seed: u64) -> Result<Curve<SmallPrime>, String> {
    if !(8..=28).contains(&bits) {
        return Err("toy curves use 8 to 28 bit primes".into());
    }
    let mut p = (1u64 << bits) - 1;
    while !(p % 4 == 3 && is_prime_u64(p)) {
        p -= 2;
    }
    let f = SmallPrime::new(p);
    let mut squares = vec![0u64; (p as usize) / 64 + 1];
    for y in 0..p {
        let s = ((y as u128 * y as u128) % p as u128) as usize;
        squares[s / 64] |= 1 << (s % 64);
    }
    let is_sq = |v: u64| (squares[v as usize / 64] >> (v % 64)) & 1 == 1;
    let mut rng = StdRng::seed_from_u64(seed ^ u64::from(bits));
    for _ in 0..200 {
        let b = rng.gen_range(1..p);
        let a = p - 3;
        let (mut points, mut roots) = (1u64, 0u64);
        for x in 0..p {
            let x2 = ((x as u128 * x as u128) % p as u128) as u64;
            let t = (x2 + a) % p;
            let r = (((t as u128 * x as u128) % p as u128) as u64 + b) % p;
            if r == 0 {
                roots += 1;
                points += 1;
            } else if is_sq(r) {
                points += 2;
            }
        }
        if roots == 0 {
            return Ok(Curve {
                name: format!("toy y²=x³−3x+{b} over F_{p} ({bits} bits)"),
                field: f,
                a,
                b,
                order: BigUint::from(points),
                order_is_prime: is_prime_u64(points),
                cofactor: 1,
            });
        }
    }
    Err("no root-free cubic found".into())
}

/// A toy curve `y² = x³ − 3x + b` over a prime `p ≡ 3 (mod 4)` just below
/// `2^bits` whose group order is `h·n` with `n` prime, `h ∈ {2, 4}`, counted
/// exactly by walking the field. The subgroup of order `n` is `[h]E`, so
/// this exercises the cofactor path where hits are common.
pub fn toy_curve_with_cofactor(bits: u32, seed: u64, h: u64) -> Result<Curve<SmallPrime>, String> {
    if !(8..=28).contains(&bits) {
        return Err("toy curves use 8 to 28 bit primes".into());
    }
    if h < 2 || !h.is_power_of_two() || h > 8 {
        return Err("toy cofactors are 2, 4 or 8".into());
    }
    let mut p = (1u64 << bits) - 1;
    while !(p % 4 == 3 && is_prime_u64(p)) {
        p -= 2;
    }
    let f = SmallPrime::new(p);
    let mut squares = vec![0u64; (p as usize) / 64 + 1];
    for y in 0..p {
        let s = ((y as u128 * y as u128) % p as u128) as usize;
        squares[s / 64] |= 1 << (s % 64);
    }
    let is_sq = |v: u64| (squares[v as usize / 64] >> (v % 64)) & 1 == 1;
    let mut rng = StdRng::seed_from_u64(seed ^ u64::from(bits) ^ (h << 40));
    for _ in 0..4000 {
        let b = rng.gen_range(1..p);
        let a = p - 3;
        let mut points = 1u64;
        for x in 0..p {
            let x2 = ((x as u128 * x as u128) % p as u128) as u64;
            let t = (x2 + a) % p;
            let r = (((t as u128 * x as u128) % p as u128) as u64 + b) % p;
            if r == 0 {
                points += 1;
            } else if is_sq(r) {
                points += 2;
            }
        }
        if points.is_multiple_of(h) && is_prime_u64(points / h) && points / h > 7 {
            let n = points / h;
            return Ok(Curve {
                name: format!("toy y²=x³−3x+{b} over F_{p} ({bits} bits), #E = {h}·{n}"),
                field: f,
                a,
                b,
                order: BigUint::from(n),
                order_is_prime: true,
                cofactor: h,
            });
        }
    }
    Err("no toy curve of that cofactor found".into())
}

/// `a + b` by the native Jacobian arithmetic, for cross-checks.
pub fn p192_native_sum(a: &Pt<Fp>, b: &Pt<Fp>, ar: &mut P192Arithmetic) -> Option<AffinePoint> {
    let j = JacobianPoint::from_affine(AffinePoint { x: a.x, y: a.y }, ar);
    j.add_affine(AffinePoint { x: b.x, y: b.y }, ar)
        .to_affine(ar)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    /// Every point of a toy curve, found by walking the field.
    fn all_points(c: &Curve<SmallPrime>) -> Vec<Pt<u64>> {
        let mut f = c.field.clone();
        let exp = c.euler_exponent(&f);
        let mut out = Vec::new();
        for x in 0..f.p {
            let r = c.rhs(&mut f, x);
            if c.is_square(&mut f, r, &exp) {
                let y = c.sqrt(&mut f, r).unwrap();
                out.push(Pt { x, y });
                out.push(Pt { x, y: f.neg(y) });
            }
        }
        out
    }

    #[test]
    fn toy_order_is_the_number_of_points_plus_infinity() {
        let c = toy_curve(12, 1).expect("toy");
        let n = all_points(&c).len() as u64 + 1;
        assert_eq!(c.order, BigUint::from(n));
        assert!(n % 2 == 1, "no 2-torsion");
    }

    #[test]
    fn exact_yield_equals_brute_force_coverage_on_a_toy_curve() {
        // Brute force: every sum of two signed base points, by the group
        // law, over the whole curve; the covered targets are the distinct
        // points reached, and must be twice the distinct abscissae.
        let c = toy_curve(14, 3).expect("toy");
        let mut f = c.field.clone();
        let base = small_x_base(&c, 40).expect("base");
        let mut signed: Vec<Pt<u64>> = Vec::new();
        for p in &base {
            signed.push(*p);
            signed.push(c.negate(&mut f, p));
        }
        let mut covered: HashSet<(u64, u64)> = HashSet::new();
        for a in &signed {
            for b in &signed {
                if let Some(s) = c.add_points(&mut f, a, b) {
                    covered.insert((s.x, s.y));
                }
            }
        }
        let mut keys = pair_sum_keys(&c, &base);
        keys.dedup();
        assert_eq!(covered.len(), 2 * keys.len());
        // And the report agrees with it.
        let rep = measure(
            &c,
            &Options {
                size: 40,
                trials: 2000,
                controls: 16,
                seed: 5,
            },
        )
        .expect("measure");
        assert_eq!(rep.distinct_sums as usize, keys.len());
        assert_eq!(rep.covered_targets, (2 * keys.len()).to_string());
        assert_eq!(rep.controls_detected, rep.controls);
        assert_eq!(rep.controls_recovered, rep.controls);
    }

    #[test]
    fn monte_carlo_rate_matches_the_exact_yield() {
        let c = toy_curve(16, 7).expect("toy");
        let rep = measure(
            &c,
            &Options {
                size: 300,
                trials: 400_000,
                controls: 8,
                seed: 11,
            },
        )
        .expect("measure");
        assert!(rep.collisions > 0, "this regime has collisions");
        assert!(rep.mc_hits > 1000, "hits {}", rep.mc_hits);
        let [lo, hi] = rep.mc_ci95;
        assert!(
            lo <= rep.exact_yield && rep.exact_yield <= hi,
            "exact {} outside [{lo}, {hi}]",
            rep.exact_yield
        );
        // Every sampled hit was recovered as a real decomposition.
        assert_eq!(rep.mc_hits_verified, 8);
    }

    #[test]
    fn p192_base_is_on_the_curve_and_has_no_collisions() {
        let c = Curve::p192();
        let mut f = c.field.clone();
        let base = small_x_base(&c, 64).expect("base");
        assert!(base.iter().all(|p| c.on_curve(&mut f, p)));
        let rep = measure(
            &c,
            &Options {
                size: 64,
                trials: 2000,
                controls: 8,
                seed: 1,
            },
        )
        .expect("measure");
        assert_eq!(rep.collisions, 0);
        assert_eq!(rep.distinct_sums, 64 * 64);
        assert_eq!(rep.mc_hits, 0);
        assert_eq!(rep.controls_detected, 8);
        assert_eq!(rep.controls_recovered, 8);
        // 2·64² / (N − 1), N ≈ 2^192.
        assert!(
            (rep.exact_yield_log2 + 179.0).abs() < 1.0,
            "{}",
            rep.exact_yield_log2
        );
    }

    #[test]
    fn p192_sum_abscissae_agree_with_the_native_jacobian_arithmetic() {
        let c = Curve::p192();
        let mut ar = P192Arithmetic::new();
        let base = small_x_base(&c, 6).expect("base");
        let mut f = c.field.clone();
        // Keys from the batched formula, against x(a ± b) from Jacobian.
        let keys: HashSet<Key> = pair_sum_keys(&c, &base).into_iter().collect();
        for i in 0..base.len() {
            for j in (i + 1)..base.len() {
                for neg in [false, true] {
                    let b = if neg {
                        c.negate(&mut f, &base[j])
                    } else {
                        base[j]
                    };
                    let s = p192_native_sum(&base[i], &b, &mut ar).expect("sum");
                    assert!(keys.contains(&ar.canonical(s.x).0));
                }
            }
        }
    }

    #[test]
    fn small_prime_inverse_and_p192_inverse_are_inverses() {
        let mut f = SmallPrime::new(1_000_003);
        for a in [1u64, 2, 12345, 999_999] {
            let i = f.inv(a).unwrap();
            assert_eq!(f.mul(a, i), 1);
        }
        let mut g = P192Field::new();
        let a = g.embed_u64(123_456_789);
        let i = g.inv(a).unwrap();
        assert!(g.mul(a, i) == g.one());
    }

    fn secg_prime(name: &str) -> Curve<MontField<3>> {
        use crate::cryptanalysis::secg_standards::secg;
        Curve::from_secg(secg(name).expect(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
    }

    #[test]
    fn mont_field_agrees_with_biguint_arithmetic_at_160_and_192_bits() {
        use crate::cryptanalysis::secg_standards::{secg, SecgField};
        for name in [
            "secp160k1",
            "secp160r1",
            "secp160r2",
            "secp192k1",
            "secp128r1",
        ] {
            let SecgField::Prime { p } = secg(name).unwrap().field.clone() else {
                panic!()
            };
            let mut f = MontField::<3>::new(&p).unwrap();
            let mut rng = StdRng::seed_from_u64(4);
            for _ in 0..50 {
                let (a, b) = (f.random(&mut rng), f.random(&mut rng));
                let (ba, bb) = (
                    BigUint::from_slice(
                        &f.key(a)
                            .iter()
                            .flat_map(|w| [*w as u32, (*w >> 32) as u32])
                            .collect::<Vec<_>>(),
                    ),
                    BigUint::from_slice(
                        &f.key(b)
                            .iter()
                            .flat_map(|w| [*w as u32, (*w >> 32) as u32])
                            .collect::<Vec<_>>(),
                    ),
                );
                let big = |k: Key| {
                    BigUint::from_slice(
                        &k.iter()
                            .flat_map(|w| [*w as u32, (*w >> 32) as u32])
                            .collect::<Vec<_>>(),
                    )
                };
                let r = f.mul(a, b);
                assert_eq!(big(f.key(r)), (&ba * &bb) % &p, "{name} mul");
                let r = f.add(a, b);
                assert_eq!(big(f.key(r)), (&ba + &bb) % &p, "{name} add");
                let r = f.sub(a, b);
                assert_eq!(big(f.key(r)), (&ba + &p - &bb) % &p, "{name} sub");
                if !f.is_zero(a) {
                    let i = f.inv(a).unwrap();
                    assert!(f.mul(a, i) == f.one(), "{name} inverse");
                }
            }
        }
    }

    #[test]
    fn secg_prime_generators_have_the_published_order() {
        for name in [
            "secp112r1",
            "secp128r1",
            "secp160k1",
            "secp160r1",
            "secp160r2",
            "secp192k1",
            "secp192r1",
        ] {
            let c = secg_prime(name);
            let mut f = c.field.clone();
            let g = c
                .secg_generator(crate::cryptanalysis::secg_standards::secg(name).unwrap())
                .unwrap();
            assert!(c.on_curve(&mut f, &g), "{name}");
            assert!(
                c.mul_big(&mut f, &g, &c.order).is_none(),
                "{name}: n·G != O"
            );
            assert!(c.mul_big(&mut f, &g, &(&c.order - 1u8)).is_some(), "{name}");
        }
    }

    #[test]
    fn montgomery_p192_gives_the_same_yield_as_the_native_p192_field() {
        let opts = Options {
            size: 48,
            trials: 500,
            controls: 8,
            seed: 3,
        };
        let native = measure(&Curve::p192(), &opts).unwrap();
        let generic = measure(&secg_prime("secp192r1"), &opts).unwrap();
        assert_eq!(native.distinct_sums, generic.distinct_sums);
        assert_eq!(native.largest_small_x, generic.largest_small_x);
        assert_eq!(native.exact_yield_log2, generic.exact_yield_log2);
        assert_eq!(generic.controls_recovered, generic.controls);
    }

    #[test]
    fn cofactor_four_base_lies_in_the_subgroup() {
        let c = secg_prime("secp112r2");
        assert_eq!(c.cofactor, 4);
        let mut f = c.field.clone();
        let base = small_x_base(&c, 5).unwrap();
        assert_eq!(base.len(), 5);
        for p in &base {
            assert!(c.on_curve(&mut f, p));
            assert!(c.mul_big(&mut f, p, &c.order).is_none(), "n·P != O");
        }
        let rep = measure(
            &c,
            &Options {
                size: 32,
                trials: 400,
                controls: 8,
                seed: 2,
            },
        )
        .unwrap();
        assert_eq!(rep.cofactor, 4);
        assert_eq!(rep.controls_recovered, rep.controls);
    }

    fn mont_matches_biguint<const L: usize>(name: &str) {
        use crate::cryptanalysis::secg_standards::{secg, SecgField};
        let SecgField::Prime { p } = secg(name).unwrap().field.clone() else {
            panic!()
        };
        let mut f = MontField::<L>::new(&p).unwrap();
        let big = |k: [u64; L]| {
            k.iter()
                .rev()
                .fold(BigUint::from(0u8), |a, &w| (a << 64u32) | BigUint::from(w))
        };
        let mut rng = StdRng::seed_from_u64(6);
        for _ in 0..30 {
            let (a, b) = (f.random(&mut rng), f.random(&mut rng));
            let (ba, bb) = (big(f.key(a)), big(f.key(b)));
            let r = f.mul(a, b);
            assert_eq!(big(f.key(r)), (&ba * &bb) % &p, "{name} mul");
            let r = f.add(a, b);
            assert_eq!(big(f.key(r)), (&ba + &bb) % &p, "{name} add");
            let r = f.sub(a, b);
            assert_eq!(big(f.key(r)), (&ba + &p - &bb) % &p, "{name} sub");
            if !f.is_zero(a) {
                let i = f.inv(a).unwrap();
                assert!(f.mul(a, i) == f.one(), "{name} inverse");
            }
        }
    }

    #[test]
    fn mont_field_agrees_with_biguint_at_every_limb_count() {
        mont_matches_biguint::<4>("secp224r1");
        mont_matches_biguint::<4>("secp256k1");
        mont_matches_biguint::<6>("secp384r1");
        mont_matches_biguint::<9>("secp521r1");
    }

    fn generator_has_the_published_order<const L: usize>(name: &str) {
        use crate::cryptanalysis::secg_standards::secg;
        let c = Curve::<MontField<L>>::from_secg(secg(name).unwrap()).unwrap();
        let mut f = c.field.clone();
        let g = c.secg_generator(secg(name).unwrap()).unwrap();
        assert!(c.on_curve(&mut f, &g), "{name}");
        assert!(
            c.mul_big(&mut f, &g, &c.order).is_none(),
            "{name}: n·G != O"
        );
    }

    #[test]
    fn wide_prime_generators_have_the_published_order() {
        generator_has_the_published_order::<4>("secp224k1");
        generator_has_the_published_order::<4>("secp224r1");
        generator_has_the_published_order::<4>("secp256k1");
        generator_has_the_published_order::<4>("secp256r1");
        generator_has_the_published_order::<6>("secp384r1");
        generator_has_the_published_order::<9>("secp521r1");
    }

    #[test]
    fn tonelli_shanks_takes_square_roots_when_p_is_1_mod_4() {
        use crate::cryptanalysis::secg_standards::secg;
        // secp224r1 has p ≡ 1 (mod 4) with 2^96 | p − 1.
        let c = Curve::<MontField<4>>::from_secg(secg("secp224r1").unwrap()).unwrap();
        let mut f = c.field.clone();
        assert_eq!(f.modulus() % 4u8, BigUint::from(1u8));
        let mut rng = StdRng::seed_from_u64(8);
        for _ in 0..20 {
            let x = f.random(&mut rng);
            let sq = f.sqr(x);
            let r = c.sqrt(&mut f, sq).expect("a square has a root");
            assert!(f.sqr(r) == sq);
        }
        // A non-square has none.
        let exp = c.euler_exponent(&f);
        let mut z = 2u64;
        let nonsquare = loop {
            let e = f.embed_u64(z);
            if !c.is_square(&mut f, e, &exp) {
                break e;
            }
            z += 1;
        };
        assert!(c.sqrt(&mut f, nonsquare).is_none());
    }

    #[test]
    fn toy_prime_curves_with_cofactor_match_monte_carlo() {
        for (bits, seed, h, size) in [(18u32, 3u64, 2u64, 90usize), (20, 5, 4, 200)] {
            let c = toy_curve_with_cofactor(bits, seed, h).expect("toy");
            assert_eq!(c.cofactor, h);
            let rep = measure(
                &c,
                &Options {
                    size,
                    trials: 400_000,
                    controls: 8,
                    seed: 11,
                },
            )
            .expect("measure");
            assert!(rep.mc_hits > 300, "h={h}: hits {}", rep.mc_hits);
            let [lo, hi] = rep.mc_ci95;
            assert!(
                lo <= rep.exact_yield && rep.exact_yield <= hi,
                "h={h}: exact {} outside [{lo}, {hi}]",
                rep.exact_yield
            );
            assert_eq!(rep.controls_recovered, rep.controls);
            assert_eq!(rep.mc_hits_verified, 8.min(rep.mc_hits));
        }
    }

    #[test]
    fn raw_base_fraction_is_one_over_h_on_a_toy_prime_curve() {
        let c = toy_curve_with_cofactor(18, 3, 2).expect("toy");
        let r = raw_subgroup_fraction(&c, 60, 3000, 4).expect("raw");
        let [lo, hi] = r.ci95;
        assert!(lo <= 0.5 && 0.5 <= hi, "{r:?}");
    }
}
