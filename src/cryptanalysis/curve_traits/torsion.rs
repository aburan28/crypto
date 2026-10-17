//! Frobenius on `E[ℓ]`: whether it acts as a scalar, decided from the
//! `ℓ`-division polynomial.
//!
//! When `ℓ | v = [O_K : Z[π]]`, `t² − 4q ≡ 0 (mod ℓ)` and Frobenius has
//! the double eigenvalue `λ ≡ t/2` on `E[ℓ]` (`λ = 1` for `ℓ = 2`).  It
//! is then either the scalar `λ` or conjugate to a Jordan block, and it
//! is the scalar exactly when `(π − λ)/ℓ ∈ End(E)`, that is when
//! `ℓ | [End(E) : Z[π]]`.  In volcano terms: the scalar case is every
//! level above the floor, the Jordan case is the floor.
//!
//! The test: with `k` the multiplicative order of `λ` mod `ℓ`, the scalar
//! satisfies `π^k = 1` on `E[ℓ]`, while a Jordan block `J` has `J^k = λ^k
//! [[1, kλ⁻¹], [0, 1]] ≠ ±1` (`ℓ ∤ k`).  And `π^k = ±1` on `E[ℓ]` exactly
//! when every root of the (squarefree) division polynomial `ψ_ℓ` lies in
//! `GF(q^k)`, i.e. `x^{q^k} ≡ x (mod ψ_ℓ)`.  So the test is exact.

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};

use super::Model;
use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::semaev_decomp::Gf2;

/// Field operations the division-polynomial arithmetic needs.
pub trait Field {
    type E: Clone + PartialEq;
    fn zero(&self) -> Self::E;
    fn one(&self) -> Self::E;
    fn int(&self, v: u64) -> Self::E;
    fn add(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn neg(&self, a: &Self::E) -> Self::E;
    fn mul(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn inv(&self, a: &Self::E) -> Self::E;
    /// A pseudo-random element fixed by `i` (splitmix64), so every use
    /// is deterministic.
    fn sample(&self, i: u64) -> Self::E;
    /// `a^q`, `q` the order of the curve's field of definition: the
    /// identity on that field itself, which every implementation here is
    /// except the extensions of `points.rs`.
    fn frobenius(&self, a: &Self::E) -> Self::E {
        a.clone()
    }
    fn is_zero(&self, a: &Self::E) -> bool {
        *a == self.zero()
    }
    fn sub(&self, a: &Self::E, b: &Self::E) -> Self::E {
        self.add(a, &self.neg(b))
    }
}

/// The splitmix64 sequence seeded by `i`.
pub(crate) fn splitmix(i: u64) -> impl Iterator<Item = u64> {
    let mut state = i.wrapping_mul(0x9e37_79b9_7f4a_7c15);
    std::iter::repeat_with(move || {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    })
}

/// `GF(p)`, `p` odd.
pub struct Fp(pub BigUint);

impl Field for Fp {
    type E = BigUint;
    fn zero(&self) -> BigUint {
        BigUint::zero()
    }
    fn one(&self) -> BigUint {
        BigUint::one()
    }
    fn int(&self, v: u64) -> BigUint {
        BigUint::from(v) % &self.0
    }
    fn add(&self, a: &BigUint, b: &BigUint) -> BigUint {
        let s = a + b;
        if s >= self.0 {
            s - &self.0
        } else {
            s
        }
    }
    fn neg(&self, a: &BigUint) -> BigUint {
        if a.is_zero() {
            BigUint::zero()
        } else {
            &self.0 - a
        }
    }
    fn mul(&self, a: &BigUint, b: &BigUint) -> BigUint {
        a * b % &self.0
    }
    fn inv(&self, a: &BigUint) -> BigUint {
        a.modpow(&(&self.0 - 2u8), &self.0)
    }
    fn sample(&self, i: u64) -> BigUint {
        let words = (self.0.bits() / 64 + 2) as usize;
        let v = splitmix(i)
            .take(words)
            .fold(BigUint::zero(), |acc, w| (acc << 64) + w);
        v % &self.0
    }
}

/// `GF(2^n)`, `n ≤ 63`, on the table-free kernel of `semaev_decomp`.
pub struct Small2(pub Gf2);

impl Field for Small2 {
    type E = u64;
    fn zero(&self) -> u64 {
        0
    }
    fn one(&self) -> u64 {
        1
    }
    fn int(&self, v: u64) -> u64 {
        v & 1
    }
    fn add(&self, a: &u64, b: &u64) -> u64 {
        a ^ b
    }
    fn neg(&self, a: &u64) -> u64 {
        *a
    }
    fn mul(&self, a: &u64, b: &u64) -> u64 {
        self.0.mul(*a, *b)
    }
    fn inv(&self, a: &u64) -> u64 {
        self.0.inv(*a)
    }
    fn sample(&self, i: u64) -> u64 {
        let n = self.0.n;
        splitmix(i).next().expect("infinite") & ((1u64 << n) - 1)
    }
}

/// `GF(2^n)` for any `n`.
pub struct Big2 {
    pub n: u32,
    pub irr: IrreduciblePoly,
}

impl Field for Big2 {
    type E = F2mElement;
    fn zero(&self) -> F2mElement {
        F2mElement::zero(self.n)
    }
    fn one(&self) -> F2mElement {
        F2mElement::one(self.n)
    }
    fn int(&self, v: u64) -> F2mElement {
        if v & 1 == 1 {
            self.one()
        } else {
            self.zero()
        }
    }
    fn add(&self, a: &F2mElement, b: &F2mElement) -> F2mElement {
        a.add(b)
    }
    fn neg(&self, a: &F2mElement) -> F2mElement {
        a.clone()
    }
    fn mul(&self, a: &F2mElement, b: &F2mElement) -> F2mElement {
        a.mul(b, &self.irr)
    }
    fn inv(&self, a: &F2mElement) -> F2mElement {
        a.flt_inverse(&self.irr).expect("non-zero")
    }
    fn sample(&self, i: u64) -> F2mElement {
        let words: Vec<u64> = splitmix(i).take(self.n.div_ceil(64) as usize).collect();
        let v = words
            .iter()
            .rev()
            .fold(BigUint::zero(), |acc, w| (acc << 64) + w);
        F2mElement::from_biguint(&(v % (BigUint::one() << self.n)), self.n)
    }
}

/// `GF(p^k) = GF(p)[z]/(m)` for `p < 2^63`, `m` monic of degree `k`
/// given by its `k` low coefficients (the registry's representation).
/// An element is its `k` coefficients, low to high.
pub struct FpK {
    pub p: u64,
    pub modulus: Vec<u64>,
}

impl FpK {
    /// `None` when `p` does not fit or the data is not canonical.
    pub fn new(p: &BigUint, modulus: &[BigUint]) -> Option<Self> {
        let p64 = p.to_u64().filter(|&v| v < 1 << 63 && v > 3)?;
        let modulus = modulus
            .iter()
            .map(|c| c.to_u64().filter(|&v| v < p64))
            .collect::<Option<Vec<_>>>()?;
        Some(FpK { p: p64, modulus })
    }

    pub fn k(&self) -> usize {
        self.modulus.len()
    }

    /// An element from the registry's coefficient vector.
    pub fn element(&self, v: &[BigUint]) -> Option<Vec<u64>> {
        (v.len() == self.k())
            .then(|| {
                v.iter()
                    .map(|c| c.to_u64().filter(|&x| x < self.p))
                    .collect()
            })
            .flatten()
    }

    fn mulmod(&self, a: u64, b: u64) -> u64 {
        (u128::from(a) * u128::from(b) % u128::from(self.p)) as u64
    }

    pub fn pow(&self, a: &[u64], e: &BigUint) -> Vec<u64> {
        let mut acc = self.one();
        for i in (0..e.bits()).rev() {
            acc = self.mul(&acc, &acc);
            if e.bit(i) {
                acc = self.mul(&acc, &a.to_vec());
            }
        }
        acc
    }
}

impl Field for FpK {
    type E = Vec<u64>;
    fn zero(&self) -> Vec<u64> {
        vec![0; self.k()]
    }
    fn one(&self) -> Vec<u64> {
        let mut v = self.zero();
        v[0] = 1;
        v
    }
    fn int(&self, v: u64) -> Vec<u64> {
        let mut e = self.zero();
        e[0] = v % self.p;
        e
    }
    fn add(&self, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
        a.iter()
            .zip(b)
            .map(|(x, y)| {
                let s = u128::from(*x) + u128::from(*y);
                (s % u128::from(self.p)) as u64
            })
            .collect()
    }
    fn neg(&self, a: &Vec<u64>) -> Vec<u64> {
        a.iter()
            .map(|&x| if x == 0 { 0 } else { self.p - x })
            .collect()
    }
    fn mul(&self, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
        let k = self.k();
        let p = u128::from(self.p);
        if 2 * u64::from(64 - self.p.leading_zeros()) + u64::from(64 - (4 * k as u64).leading_zeros())
            <= 126
        {
            // Products are below p² and no slot sums more than 2k of
            // them (and k more in the reduction), so u128 holds everything
            // unreduced: one reduction per coefficient instead of one per
            // product.
            let mut t = vec![0u128; 2 * k - 1];
            for (i, &x) in a.iter().enumerate() {
                if x == 0 {
                    continue;
                }
                for (j, &y) in b.iter().enumerate() {
                    t[i + j] += u128::from(x) * u128::from(y);
                }
            }
            for i in (k..2 * k - 1).rev() {
                let c = t[i] % p;
                if c == 0 {
                    continue;
                }
                for (j, &m) in self.modulus.iter().enumerate() {
                    // z^k = −Σ m_j z^j.
                    t[i - k + j] += c * ((p - u128::from(m)) % p);
                }
            }
            return t[..k].iter().map(|&v| (v % p) as u64).collect();
        }
        // Schoolbook product, each coefficient reduced mod p.
        let mut t = vec![0u64; 2 * k - 1];
        for (i, &x) in a.iter().enumerate() {
            if x == 0 {
                continue;
            }
            for (j, &y) in b.iter().enumerate() {
                let s = u128::from(t[i + j]) + u128::from(x) * u128::from(y) % p;
                t[i + j] = (s % p) as u64;
            }
        }
        // z^k = −Σ m_j z^j.
        for i in (k..2 * k - 1).rev() {
            let c = t[i];
            if c == 0 {
                continue;
            }
            for (j, &m) in self.modulus.iter().enumerate() {
                let sub = self.mulmod(c, m);
                let cur = t[i - k + j];
                t[i - k + j] = if cur >= sub {
                    cur - sub
                } else {
                    cur + self.p - sub
                };
            }
        }
        t.truncate(k);
        t
    }
    fn sample(&self, i: u64) -> Vec<u64> {
        // Every coefficient from the sequence, so samples leave every subfield.
        splitmix(i).take(self.k()).map(|z| z % self.p).collect()
    }
    fn inv(&self, a: &Vec<u64>) -> Vec<u64> {
        if a.iter().all(|&c| c == 0) {
            // As a^(q − 2) has it.
            return a.clone();
        }
        if self.k() == 1 {
            return vec![self.inv_mod(a[0])];
        }
        // Extended Euclid in GF(p)[z]: s·a ≡ gcd (mod m), the gcd a unit.
        let p = self.p;
        let strip = |mut v: Vec<u64>| {
            while v.last() == Some(&0) {
                v.pop();
            }
            v
        };
        let mut m = self.modulus.clone();
        m.push(1);
        let (mut r0, mut r1) = (m, strip(a.clone()));
        let (mut s0, mut s1): (Vec<u64>, Vec<u64>) = (Vec::new(), vec![1]);
        while r1.len() > 1 {
            // r0 = qt·r1 + rem.
            let lead_inv = self.inv_mod(*r1.last().expect("non-zero"));
            let mut rem = r0.clone();
            let mut qt = vec![0u64; rem.len() + 1 - r1.len()];
            while rem.len() >= r1.len() {
                let c = self.mulmod(*rem.last().expect("non-empty"), lead_inv);
                let shift = rem.len() - r1.len();
                qt[shift] = c;
                for (i, &x) in r1.iter().enumerate() {
                    let sub = self.mulmod(c, x);
                    rem[shift + i] = (rem[shift + i] + p - sub) % p;
                }
                rem = strip(rem);
                if rem.is_empty() {
                    break;
                }
            }
            // s2 = s0 − qt·s1.
            let mut s2 = vec![0u64; (qt.len() + s1.len()).max(s0.len())];
            for (i, &x) in s0.iter().enumerate() {
                s2[i] = x;
            }
            for (i, &x) in qt.iter().enumerate() {
                if x == 0 {
                    continue;
                }
                for (j, &y) in s1.iter().enumerate() {
                    let sub = self.mulmod(x, y);
                    s2[i + j] = (s2[i + j] + p - sub) % p;
                }
            }
            (r0, r1) = (r1, rem);
            (s0, s1) = (s1, strip(s2));
            if r1.is_empty() {
                // A common factor with m: no inverse in GF(p)[z]/(m), so
                // keep the old a^(q − 2) rather than fail.
                let q = BigUint::from(p).pow(self.k() as u32);
                return self.pow(a, &(q - 2u8));
            }
        }
        let c = self.inv_mod(r1[0]);
        let mut out = vec![0u64; self.k()];
        for (i, &x) in s1.iter().enumerate() {
            out[i] = self.mulmod(x, c);
        }
        out
    }
}

impl FpK {
    /// Whether `gcd(a, m) = 1` in `GF(p)[z]`, `m` the modulus (which need
    /// not be irreducible here): Euclid on the remainders.
    pub fn coprime(&self, a: &[u64]) -> bool {
        let p = self.p;
        let strip = |mut v: Vec<u64>| {
            while v.last() == Some(&0) {
                v.pop();
            }
            v
        };
        let mut m = self.modulus.clone();
        m.push(1);
        let (mut r0, mut r1) = (m, strip(a.to_vec()));
        while r1.len() > 1 {
            let lead_inv = self.inv_mod(*r1.last().expect("non-zero"));
            let mut rem = r0;
            while rem.len() >= r1.len() {
                let c = self.mulmod(*rem.last().expect("non-empty"), lead_inv);
                let shift = rem.len() - r1.len();
                for (i, &x) in r1.iter().enumerate() {
                    let sub = self.mulmod(c, x);
                    rem[shift + i] = (rem[shift + i] + p - sub) % p;
                }
                rem = strip(rem);
                if rem.is_empty() {
                    break;
                }
            }
            (r0, r1) = (r1, rem);
        }
        r1.len() == 1
    }

    fn inv_mod(&self, a: u64) -> u64 {
        let (mut r0, mut r1) = (i128::from(self.p), i128::from(a % self.p));
        let (mut s0, mut s1) = (0i128, 1i128);
        while r1 != 0 {
            let q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (s0, s1) = (s1, s0 - q * s1);
        }
        s0.rem_euclid(i128::from(self.p)) as u64
    }
}

/// A polynomial, coefficients low to high, no trailing zeros.
pub(crate) type Poly<E> = Vec<E>;

pub(crate) fn trim<F: Field>(f: &F, mut p: Poly<F::E>) -> Poly<F::E> {
    while p.last().is_some_and(|c| f.is_zero(c)) {
        p.pop();
    }
    p
}

pub(crate) fn padd<F: Field>(f: &F, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
    let n = a.len().max(b.len());
    let z = f.zero();
    let out = (0..n)
        .map(|i| f.add(a.get(i).unwrap_or(&z), b.get(i).unwrap_or(&z)))
        .collect();
    trim(f, out)
}

pub(crate) fn psub<F: Field>(f: &F, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
    let nb: Poly<F::E> = b.iter().map(|c| f.neg(c)).collect();
    padd(f, a, &nb)
}

pub(crate) fn pmul<F: Field>(f: &F, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![f.zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if f.is_zero(x) {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] = f.add(&out[i + j], &f.mul(x, y));
        }
    }
    trim(f, out)
}

pub(crate) fn pscale<F: Field>(f: &F, a: &Poly<F::E>, c: &F::E) -> Poly<F::E> {
    trim(f, a.iter().map(|x| f.mul(x, c)).collect())
}

/// `a mod m` for monic `m`.
pub(crate) fn prem<F: Field>(f: &F, a: &Poly<F::E>, m: &Poly<F::E>) -> Poly<F::E> {
    let d = m.len() - 1;
    let mut r = a.clone();
    while r.len() > d {
        let top = r.pop().expect("non-empty");
        if f.is_zero(&top) {
            continue;
        }
        let shift = r.len() - d;
        for (i, c) in m[..d].iter().enumerate() {
            r[shift + i] = f.sub(&r[shift + i], &f.mul(&top, c));
        }
    }
    trim(f, r)
}

/// Exact division by `x`.
fn div_x<F: Field>(f: &F, a: &Poly<F::E>) -> Poly<F::E> {
    assert!(a.first().is_none_or(|c| f.is_zero(c)), "not divisible by x");
    a.iter().skip(1).cloned().collect()
}

pub(crate) fn monic<F: Field>(f: &F, a: &Poly<F::E>) -> Poly<F::E> {
    let lead = f.inv(a.last().expect("non-zero polynomial"));
    pscale(f, a, &lead)
}

/// `base^e mod m`.
pub(crate) fn ppow<F: Field>(f: &F, base: &Poly<F::E>, e: &BigUint, m: &Poly<F::E>) -> Poly<F::E> {
    let mut acc = vec![f.one()];
    for i in (0..e.bits()).rev() {
        acc = prem(f, &pmul(f, &acc, &acc), m);
        if e.bit(i) {
            acc = prem(f, &pmul(f, &acc, base), m);
        }
    }
    acc
}

/// `g(h) mod m` by Horner.
pub(crate) fn compose<F: Field>(
    f: &F,
    g: &Poly<F::E>,
    h: &Poly<F::E>,
    m: &Poly<F::E>,
) -> Poly<F::E> {
    let mut acc: Poly<F::E> = Vec::new();
    for c in g.iter().rev() {
        acc = padd(f, &prem(f, &pmul(f, &acc, h), m), &trim(f, vec![c.clone()]));
    }
    acc
}

impl FpK {
    /// The `i`-th element in base-`p` digit order (`0, 1, …, p − 1, z, …`).
    pub fn nth(&self, mut i: u64) -> Vec<u64> {
        (0..self.k())
            .map(|_| {
                let d = i % self.p;
                i /= self.p;
                d
            })
            .collect()
    }

    /// A square root of `a`, or `None` when `a` is a non-square
    /// (Tonelli–Shanks; `q = p^k` is odd).
    pub fn sqrt(&self, a: &[u64]) -> Option<Vec<u64>> {
        let a = a.to_vec();
        if self.is_zero(&a) {
            return Some(a);
        }
        let q = BigUint::from(self.p).pow(self.k() as u32);
        let qm1 = &q - 1u8;
        let half = &qm1 >> 1;
        let minus_one = self.neg(&self.one());
        if self.pow(&a, &half) != self.one() {
            return None;
        }
        let s = qm1.trailing_zeros().expect("q − 1 > 0");
        let t = &qm1 >> s;
        // A non-square, found deterministically.  For even k every
        // element of GF(p) is a square, so search i + z (or i when k = 1).
        let z = (0..)
            .map(|i| {
                let mut z = self.int(i);
                if self.k() > 1 {
                    z[1] = 1;
                }
                z
            })
            .find(|z| self.pow(z, &half) == minus_one)
            .expect("half of GF(q)* is non-square");
        let mut m = s;
        let mut c = self.pow(&z, &t);
        let mut tt = self.pow(&a, &t);
        let mut r = self.pow(&a, &((&t + 1u8) >> 1));
        while tt != self.one() {
            let mut i = 0;
            let mut sq = tt.clone();
            while sq != self.one() {
                sq = self.mul(&sq, &sq);
                i += 1;
            }
            let mut b = c.clone();
            for _ in 0..(m - i - 1) {
                b = self.mul(&b, &b);
            }
            m = i;
            c = self.mul(&b, &b);
            tt = self.mul(&tt, &c);
            r = self.mul(&r, &b);
        }
        Some(r)
    }
}

/// `(x, y)` lies on `y² = x³ + ax + b`.
pub fn on_curve<F: Field>(f: &F, a: &F::E, b: &F::E, x: &F::E, y: &F::E) -> bool {
    let rhs = f.add(&f.mul(&f.add(&f.mul(x, x), a), x), b);
    f.mul(y, y) == rhs
}

/// A point of `y² = x³ + ax + b` in Jacobian coordinates
/// `(X : Y : Z) ↦ (X/Z², Y/Z³)`; `Z = 0` is the point at infinity.
#[derive(Clone)]
pub(crate) struct Jac<E> {
    pub(crate) x: E,
    pub(crate) y: E,
    pub(crate) z: E,
}

pub(crate) fn jac_double<F: Field>(f: &F, a: &F::E, p: &Jac<F::E>) -> Jac<F::E> {
    if f.is_zero(&p.z) || f.is_zero(&p.y) {
        return Jac {
            x: f.one(),
            y: f.one(),
            z: f.zero(),
        };
    }
    let xx = f.mul(&p.x, &p.x);
    let yy = f.mul(&p.y, &p.y);
    let zz = f.mul(&p.z, &p.z);
    let s = f.mul(&f.int(4), &f.mul(&p.x, &yy));
    let m = f.add(&f.mul(&f.int(3), &xx), &f.mul(a, &f.mul(&zz, &zz)));
    let x3 = f.sub(&f.mul(&m, &m), &f.add(&s, &s));
    let y3 = f.sub(
        &f.mul(&m, &f.sub(&s, &x3)),
        &f.mul(&f.int(8), &f.mul(&yy, &yy)),
    );
    let z3 = f.mul(&f.int(2), &f.mul(&p.y, &p.z));
    Jac {
        x: x3,
        y: y3,
        z: z3,
    }
}

pub(crate) fn jac_add<F: Field>(f: &F, a: &F::E, p: &Jac<F::E>, q: &Jac<F::E>) -> Jac<F::E> {
    if f.is_zero(&p.z) {
        return q.clone();
    }
    if f.is_zero(&q.z) {
        return p.clone();
    }
    let z1z1 = f.mul(&p.z, &p.z);
    let z2z2 = f.mul(&q.z, &q.z);
    let u1 = f.mul(&p.x, &z2z2);
    let u2 = f.mul(&q.x, &z1z1);
    let s1 = f.mul(&p.y, &f.mul(&q.z, &z2z2));
    let s2 = f.mul(&q.y, &f.mul(&p.z, &z1z1));
    if u1 == u2 {
        return if s1 == s2 {
            jac_double(f, a, p)
        } else {
            Jac {
                x: f.one(),
                y: f.one(),
                z: f.zero(),
            }
        };
    }
    let h = f.sub(&u2, &u1);
    let r = f.sub(&s2, &s1);
    let hh = f.mul(&h, &h);
    let hhh = f.mul(&h, &hh);
    let v = f.mul(&u1, &hh);
    let x3 = f.sub(&f.sub(&f.mul(&r, &r), &hhh), &f.add(&v, &v));
    let y3 = f.sub(&f.mul(&r, &f.sub(&v, &x3)), &f.mul(&s1, &hhh));
    let z3 = f.mul(&h, &f.mul(&p.z, &q.z));
    Jac {
        x: x3,
        y: y3,
        z: z3,
    }
}

/// Whether `[n](x, y) = O` on `y² = x³ + ax + b`.
pub fn kills<F: Field>(f: &F, a: &F::E, x: &F::E, y: &F::E, n: &BigUint) -> bool {
    let g = Jac {
        x: x.clone(),
        y: y.clone(),
        z: f.one(),
    };
    let mut acc = Jac {
        x: f.one(),
        y: f.one(),
        z: f.zero(),
    };
    for i in (0..n.bits()).rev() {
        acc = jac_double(f, a, &acc);
        if n.bit(i) {
            acc = jac_add(f, a, &acc, &g);
        }
    }
    f.is_zero(&acc.z)
}

/// The division polynomials `f_0 … f_top` of `y² = x³ + ax + b` (odd
/// characteristic), with `ψ_n = f_n` for odd `n` and `ψ_n = 2y·f_n` for
/// even `n`.
pub fn division_prime<F: Field>(f: &F, a: &F::E, b: &F::E, top: usize) -> Vec<Poly<F::E>> {
    let c = |v: u64| f.int(v);
    let (a2, ab, b2) = (f.mul(a, a), f.mul(a, b), f.mul(b, b));
    let a3 = f.mul(&a2, a);
    let mut fs: Vec<Poly<F::E>> = vec![
        Vec::new(),
        vec![f.one()],
        vec![f.one()],
        trim(
            f,
            vec![
                f.neg(&a2),
                f.mul(&c(12), b),
                f.mul(&c(6), a),
                f.zero(),
                c(3),
            ],
        ),
        pscale(
            f,
            &trim(
                f,
                vec![
                    f.neg(&f.add(&f.mul(&c(8), &b2), &a3)),
                    f.neg(&f.mul(&c(4), &ab)),
                    f.neg(&f.mul(&c(5), &a2)),
                    f.mul(&c(20), b),
                    f.mul(&c(5), a),
                    f.zero(),
                    f.one(),
                ],
            ),
            &c(2),
        ),
    ];
    // R² = (2y)⁴ = 16 (x³ + ax + b)².
    let cubic = trim(f, vec![b.clone(), a.clone(), f.zero(), f.one()]);
    let r2 = pscale(f, &pmul(f, &cubic, &cubic), &c(16));
    for n in 5..=top {
        let m = n / 2;
        let p = if n % 2 == 1 {
            let t1 = pmul(f, &fs[m + 2], &pmul(f, &fs[m], &pmul(f, &fs[m], &fs[m])));
            let t2 = pmul(
                f,
                &fs[m - 1],
                &pmul(f, &fs[m + 1], &pmul(f, &fs[m + 1], &fs[m + 1])),
            );
            if m % 2 == 0 {
                psub(f, &pmul(f, &r2, &t1), &t2)
            } else {
                psub(f, &t1, &pmul(f, &r2, &t2))
            }
        } else {
            let t1 = pmul(f, &fs[m + 2], &pmul(f, &fs[m - 1], &fs[m - 1]));
            let t2 = pmul(f, &fs[m - 2], &pmul(f, &fs[m + 1], &fs[m + 1]));
            pmul(f, &fs[m], &psub(f, &t1, &t2))
        };
        fs.push(p);
    }
    fs.truncate(top + 1);
    fs
}

/// The division polynomials `f_0 … f_top` of `y² + xy = x³ + ax² + b`
/// over `GF(2^n)`; they do not depend on `a`.
pub fn division_binary<F: Field>(f: &F, b: &F::E, top: usize) -> Vec<Poly<F::E>> {
    let mut fs: Vec<Poly<F::E>> = vec![
        Vec::new(),
        vec![f.one()],
        vec![f.zero(), f.one()],
        trim(f, vec![b.clone(), f.zero(), f.zero(), f.one(), f.one()]),
        trim(
            f,
            vec![
                f.zero(),
                f.zero(),
                b.clone(),
                f.zero(),
                f.zero(),
                f.zero(),
                f.one(),
            ],
        ),
    ];
    for n in 5..=top {
        let m = n / 2;
        let p = if n % 2 == 1 {
            let t1 = pmul(f, &fs[m + 2], &pmul(f, &fs[m], &pmul(f, &fs[m], &fs[m])));
            let t2 = pmul(
                f,
                &fs[m - 1],
                &pmul(f, &fs[m + 1], &pmul(f, &fs[m + 1], &fs[m + 1])),
            );
            padd(f, &t1, &t2)
        } else {
            let t1 = pmul(f, &fs[m + 2], &pmul(f, &fs[m - 1], &fs[m - 1]));
            let t2 = pmul(f, &fs[m - 2], &pmul(f, &fs[m + 1], &fs[m + 1]));
            div_x(f, &pmul(f, &fs[m], &padd(f, &t1, &t2)))
        };
        fs.push(p);
    }
    fs.truncate(top + 1);
    fs
}

/// Whether `x^{q^k} ≡ x (mod m)`, `m` monic.
fn frobenius_fixes<F: Field>(f: &F, q: &BigUint, k: u32, m: &Poly<F::E>) -> bool {
    let x = vec![f.zero(), f.one()];
    let deg = m.len() - 1;
    let x_q = ppow(f, &x, q, m);
    let mut acc = x_q.clone();
    for _ in 1..k {
        // x^{q^{j+1}} = (x^q)(x^{q^j}): composition costs `deg` products,
        // powering costs `log q`; take the cheaper.
        acc = if (deg as u64) < q.bits() {
            compose(f, &x_q, &acc, m)
        } else {
            ppow(f, &acc, q, m)
        };
    }
    acc == trim(f, x)
}

/// `(a div b, a mod b)` for any non-zero `b`.
pub(crate) fn pdivrem<F: Field>(f: &F, a: &Poly<F::E>, b: &Poly<F::E>) -> (Poly<F::E>, Poly<F::E>) {
    let db = b.len() - 1;
    let lead_inv = f.inv(b.last().expect("non-zero divisor"));
    let mut r = a.clone();
    let mut quo = vec![f.zero(); a.len().saturating_sub(db).max(1)];
    while r.len() > db && !r.is_empty() {
        let shift = r.len() - 1 - db;
        let c = f.mul(r.last().expect("non-empty"), &lead_inv);
        for (i, bc) in b.iter().enumerate() {
            r[shift + i] = f.sub(&r[shift + i], &f.mul(&c, bc));
        }
        quo[shift] = c;
        r = trim(f, r);
    }
    (trim(f, quo), r)
}

/// `a⁻¹ mod m` by the extended Euclidean algorithm; `None` when
/// `gcd(a, m) ≠ 1`.
pub(crate) fn pinv_mod<F: Field>(f: &F, a: &Poly<F::E>, m: &Poly<F::E>) -> Option<Poly<F::E>> {
    let (mut r0, mut r1) = (m.clone(), prem(f, a, m));
    let (mut s0, mut s1): (Poly<F::E>, Poly<F::E>) = (Vec::new(), vec![f.one()]);
    while !r1.is_empty() {
        let (quo, rem) = pdivrem(f, &r0, &r1);
        let s2 = psub(f, &s0, &pmul(f, &quo, &s1));
        (r0, r1) = (r1, rem);
        (s0, s1) = (s1, s2);
    }
    if r0.len() != 1 {
        return None;
    }
    let c = f.inv(&r0[0]);
    Some(prem(f, &pscale(f, &s0, &c), m))
}

/// `x([m]P) mod M` on `y² = x³ + ax + b` from the division polynomials
/// `fs` (`fs.len() > m + 1`): `x − ψ_{m−1}ψ_{m+1}/ψ_m²`, with `ψ` even
/// index `= 2y·f` and `4y² = r = 4(x³ + ax + b)`.  `None` when `ψ_m²` is
/// not invertible mod `M`.
fn x_multiple<F: Field>(
    f: &F,
    fs: &[Poly<F::E>],
    r: &Poly<F::E>,
    m: usize,
    modulus: &Poly<F::E>,
) -> Option<Poly<F::E>> {
    let x = vec![f.zero(), f.one()];
    if m == 1 {
        return Some(prem(f, &x, modulus));
    }
    let around = prem(f, &pmul(f, &fs[m - 1], &fs[m + 1]), modulus);
    let square = prem(f, &pmul(f, &fs[m], &fs[m]), modulus);
    let (num, den) = if m % 2 == 1 {
        (prem(f, &pmul(f, r, &around), modulus), square)
    } else {
        (around, prem(f, &pmul(f, r, &square), modulus))
    };
    let inv = pinv_mod(f, &den, modulus)?;
    Some(psub(
        f,
        &prem(f, &x, modulus),
        &prem(f, &pmul(f, &num, &inv), modulus),
    ))
}

/// The largest `ℓ^j` [`scalar_on_power`] tries.
pub const POWER_TEST_MAX: u64 = 64;

/// Whether Frobenius acts on `E[ℓ^j]` as a scalar, `j ≥ 2`, which holds
/// exactly when `j ≤ v_ℓ([End(E) : Z[π]])`; so with the volcano depth `d`
/// the level is `d` minus the largest such `j`.
///
/// Write `n = ℓ^j`.  A scalar `λ` has `2λ ≡ t (mod n)` (the trace of an
/// element of `n·End(E)` is divisible by `n`), which fixes `λ` for odd `ℓ`
/// and leaves two candidates mod `2^j` for `ℓ = 2`.  The test is
/// `x^q ≡ x([λ]P) (mod f_n)`, the roots of `f_n` being the x-coordinates
/// of `E[n]` outside `E[2]`: it says `πP = ±λP` for every such `P`, and
/// applied to a basis `e₁, e₂` and to `e₁ + e₂` (all of order `n`, and
/// `2e ≠ O`) that forces `π = ±λ`, a scalar.  So the test is exact.
/// `None` for models without the arithmetic, `ℓ` the characteristic, or a
/// test beyond [`SCALAR_TEST_COST`].
pub fn scalar_on_power(model: &Model, trace: &BigInt, ell: u64, j: u32) -> Option<bool> {
    let n = ell
        .checked_pow(j)
        .filter(|&n| n <= POWER_TEST_MAX && j >= 2)?;
    let tn = trace.mod_floor(&BigInt::from(n)).to_u64().expect("small");
    let candidates: Vec<u64> = if ell == 2 {
        if tn % 2 == 1 {
            return None;
        }
        vec![tn / 2, (tn / 2 + n / 2) % n]
    } else {
        // 2 is invertible mod n: λ = t·(n + 1)/2.
        vec![(u128::from(tn) * u128::from(n.div_ceil(2)) % u128::from(n)) as u64]
    };
    let q = model.q();
    let deg = (n * n) / 2;
    let limbs = match model {
        Model::Extension { k: ext, .. } => u64::from(*ext),
        _ => q.bits().div_ceil(64),
    };
    if cost(deg, q.bits(), 1, limbs).saturating_mul(candidates.len() as u64) > SCALAR_TEST_COST {
        return None;
    }
    fn run<F: Field>(
        f: &F,
        a: &F::E,
        b: &F::E,
        q: &BigUint,
        n: u64,
        lambdas: &[u64],
    ) -> Option<bool> {
        let top = lambdas.iter().copied().max().unwrap_or(1).max(n) as usize + 1;
        let fs = division_prime(f, a, b, top);
        let modulus = monic(f, &fs[n as usize]);
        let r = pscale(
            f,
            &trim(f, vec![b.clone(), a.clone(), f.zero(), f.one()]),
            &f.int(4),
        );
        let x_q = ppow(f, &vec![f.zero(), f.one()], q, &modulus);
        for &lambda in lambdas {
            if lambda == 0 {
                continue;
            }
            if x_multiple(f, &fs, &r, lambda as usize, &modulus)? == x_q {
                return Some(true);
            }
        }
        Some(false)
    }
    match model {
        Model::Prime { p, a, b } => {
            if (p % ell).is_zero() || p <= &BigUint::from(3u8) {
                return None;
            }
            run(&Fp(p.clone()), a, b, &q, n, &candidates)
        }
        Model::Extension {
            p, modulus, a, b, ..
        } => {
            let f = FpK::new(p, modulus)?;
            if (BigUint::from(ell) % p).is_zero() {
                return None;
            }
            let (a, b) = (f.element(a)?, f.element(b)?);
            run(&f, &a, &b, &q, n, &candidates)
        }
        Model::Binary { .. } | Model::Other { .. } => None,
    }
}

/// The order `k` of the double eigenvalue `λ` of Frobenius mod `ℓ`; 0
/// when `λ ≡ 0`, which only happens for `ℓ` the characteristic.
pub(crate) fn eigen_order(ell: u64, trace: &BigInt) -> u32 {
    let lambda = if ell == 2 {
        1
    } else {
        let t = trace.mod_floor(&BigInt::from(ell)).to_u64().expect("small");
        t * ell.div_ceil(2) % ell
    };
    if lambda == 0 {
        return 0;
    }
    let mut k = 1;
    let mut x = lambda;
    while x != 1 {
        x = x * lambda % ell;
        k += 1;
    }
    k
}

/// Polynomial products the test may spend, in field multiplications
/// weighted by limb count: a guard on build time, not on correctness.
pub const SCALAR_TEST_COST: u64 = 40_000_000;

/// Primes `ℓ` above this are never tested (`ψ_ℓ` has degree `(ℓ² − 1)/2`).
pub const SCALAR_TEST_MAX_ELL: u64 = 1000;

fn cost(deg: u64, q_bits: u64, k: u32, limbs: u64) -> u64 {
    let step = q_bits.min(deg).max(1);
    deg.saturating_mul(deg)
        .saturating_mul(q_bits.saturating_add(u64::from(k).saturating_mul(step)))
        .saturating_mul(limbs.saturating_mul(limbs))
}

/// Whether Frobenius acts on `E[ℓ]` as a scalar.  Valid only when
/// `ℓ | t² − 4q` (the caller checks `ℓ` divides `v`).  `None` when `ℓ` is
/// the characteristic, the model is not one this module handles, or the
/// test would exceed [`SCALAR_TEST_COST`].
pub fn frobenius_is_scalar(model: &Model, trace: &BigInt, ell: u64) -> Option<bool> {
    if ell > SCALAR_TEST_MAX_ELL {
        return None;
    }
    let q = model.q();
    let k = eigen_order(ell, trace);
    if k == 0 {
        return None;
    }
    let deg = if ell == 2 { 3 } else { (ell * ell - 1) / 2 };
    // A field multiplication costs about limbs² word products; over
    // GF(p^k) with one-word coefficients, k².
    let limbs = match model {
        Model::Extension { k: ext, .. } => u64::from(*ext),
        _ => q.bits().div_ceil(64),
    };
    if cost(deg, q.bits(), k, limbs) > SCALAR_TEST_COST {
        return None;
    }
    match model {
        Model::Prime { p, a, b } => {
            if ell == 2 {
                let f = Fp(p.clone());
                let m = vec![b.clone(), a.clone(), BigUint::zero(), BigUint::one()];
                return Some(frobenius_fixes(&f, &q, k, &trim(&f, m)));
            }
            if (p % ell).is_zero() || p <= &BigUint::from(3u8) {
                return None;
            }
            let f = Fp(p.clone());
            let psi = division_prime(&f, a, b, ell as usize).pop().expect("ψ_ℓ");
            Some(frobenius_fixes(&f, &q, k, &monic(&f, &psi)))
        }
        Model::Binary { n, irr, b, .. } => {
            if ell == 2 {
                return None;
            }
            if *n <= 63 {
                let f = Small2(Gf2::new(irr));
                let bb = f.0.from_element(b);
                let psi = division_binary(&f, &bb, ell as usize).pop().expect("ψ_ℓ");
                Some(frobenius_fixes(&f, &q, k, &monic(&f, &psi)))
            } else {
                let f = Big2 {
                    n: *n,
                    irr: irr.clone(),
                };
                let psi = division_binary(&f, b, ell as usize).pop().expect("ψ_ℓ");
                Some(frobenius_fixes(&f, &q, k, &monic(&f, &psi)))
            }
        }
        Model::Extension {
            p, modulus, a, b, ..
        } => {
            let f = FpK::new(p, modulus)?;
            if (BigUint::from(ell) % p).is_zero() {
                return None;
            }
            let (a, b) = (f.element(a)?, f.element(b)?);
            let psi = if ell == 2 {
                trim(&f, vec![b, a, f.zero(), f.one()])
            } else {
                division_prime(&f, &a, &b, ell as usize).pop().expect("ψ_ℓ")
            };
            Some(frobenius_fixes(&f, &q, k, &monic(&f, &psi)))
        }
        Model::Other { .. } => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Points `P` of `y² = x³ + ax + b` over `GF(p)` with `[ℓ]P = O`.
    fn ell_torsion_prime(p: u64, a: u64, b: u64, ell: u64) -> usize {
        use crate::ecc::field::FieldElement;
        use crate::ecc::point::Point;
        let pb = BigUint::from(p);
        let fa = FieldElement::new(BigUint::from(a), pb.clone());
        let mut count = 1;
        for x in 0..p {
            for y in 0..p {
                if (y * y) % p == (x * x % p * x + a * x + b) % p {
                    let pt = Point::Affine {
                        x: FieldElement::new(BigUint::from(x), pb.clone()),
                        y: FieldElement::new(BigUint::from(y), pb.clone()),
                    };
                    if pt.scalar_mul_vartime(&BigUint::from(ell), &fa) == Point::Infinity {
                        count += 1;
                    }
                }
            }
        }
        count
    }

    #[test]
    fn division_polynomials_have_the_expected_degrees() {
        let f = Fp(BigUint::from(1_000_003u64));
        let fs = division_prime(&f, &BigUint::from(5u8), &BigUint::from(7u8), 13);
        for n in [3usize, 5, 7, 11, 13] {
            assert_eq!(fs[n].len() - 1, (n * n - 1) / 2, "deg f_{n}");
        }
        let g = Small2(Gf2::new(&IrreduciblePoly {
            degree: 7,
            low_terms: vec![0, 1],
        }));
        let fb = division_binary(&g, &1, 13);
        for n in [3usize, 5, 7, 11, 13] {
            assert_eq!(fb[n].len() - 1, (n * n - 1) / 2, "binary deg f_{n}");
        }
    }

    /// When `λ ≡ 1 (mod ℓ)` the test asks whether all of `E[ℓ]` is
    /// rational; compare with a brute-force count of `E(GF(p))[ℓ]`.
    #[test]
    fn the_scalar_test_agrees_with_counting_rational_torsion() {
        let mut checked = [0usize; 2];
        for p in [101u64, 103, 107, 109, 113, 127, 131, 137] {
            for a in 0..12u64 {
                for b in 1..12u64 {
                    let disc = (4 * a * a * a + 27 * b * b) % p;
                    if disc == 0 {
                        continue;
                    }
                    let n = super::super::certify::count_prime(p, a, b);
                    let t = p as i64 + 1 - n as i64;
                    for ell in [2u64, 3, 5] {
                        let delta = t * t - 4 * p as i64;
                        let l = ell as i64;
                        if delta % (l * l) != 0 || eigen_order(ell, &BigInt::from(t)) != 1 {
                            continue;
                        }
                        let model = Model::Prime {
                            p: BigUint::from(p),
                            a: BigUint::from(a),
                            b: BigUint::from(b),
                        };
                        let scalar = frobenius_is_scalar(&model, &BigInt::from(t), ell).unwrap();
                        let full = ell_torsion_prime(p, a, b, ell) == (ell * ell) as usize;
                        assert_eq!(scalar, full, "p={p} a={a} b={b} ℓ={ell}");
                        checked[usize::from(scalar)] += 1;
                    }
                }
            }
        }
        assert!(checked[0] > 0 && checked[1] > 0, "{checked:?}");
    }

    /// The same comparison over `GF(p^k)`, `k > 1`, which also checks the
    /// field and the Jacobian arithmetic: `[#E]P = O` for every point.
    #[test]
    fn the_extension_field_scalar_test_agrees_with_counting_rational_torsion() {
        let mut checked = [0usize; 2];
        let mut certified = [0usize; 2];
        // (p, low coefficients of an irreducible monic modulus)
        for (p, low) in [
            (5u64, vec![3u64, 0]),
            (7, vec![1, 0]),
            (13, vec![11, 0]),
            (5, vec![1, 1, 0]),
        ] {
            let f = FpK { p, modulus: low };
            let k = f.k() as u32;
            let q = p.pow(k);
            let els: Vec<Vec<u64>> = (0..q)
                .map(|mut i| {
                    (0..k)
                        .map(|_| {
                            let d = i % p;
                            i /= p;
                            d
                        })
                        .collect()
                })
                .collect();
            for e in els.iter().skip(1) {
                assert_eq!(f.mul(e, &f.inv(e)), f.one(), "inverse in GF({p}^{k})");
            }
            let mut roots: std::collections::HashMap<Vec<u64>, Vec<Vec<u64>>> = Default::default();
            for y in &els {
                roots.entry(f.mul(y, y)).or_default().push(y.clone());
            }
            for e in &els {
                match f.sqrt(e) {
                    Some(r) => assert_eq!(&f.mul(&r, &r), e, "sqrt in GF({p}^{k})"),
                    None => assert!(!roots.contains_key(e), "missed a square root"),
                }
            }
            for a in els.iter().take(12) {
                for b in els.iter().skip(1).take(12) {
                    let disc = f.add(
                        &f.mul(&f.int(4), &f.mul(a, &f.mul(a, a))),
                        &f.mul(&f.int(27), &f.mul(b, b)),
                    );
                    if f.is_zero(&disc) {
                        continue;
                    }
                    let mut points = Vec::new();
                    for x in &els {
                        let rhs = f.add(&f.mul(&f.add(&f.mul(x, x), a), x), b);
                        for y in roots.get(&rhs).into_iter().flatten() {
                            assert!(on_curve(&f, a, b, x, y));
                            points.push((x.clone(), y.clone()));
                        }
                    }
                    let n = points.len() as u64 + 1;
                    for (x, y) in &points {
                        assert!(kills(&f, a, x, y, &BigUint::from(n)), "[#E]P ≠ O");
                    }
                    // The point-order certificate accepts #E, and never
                    // certifies a wrong order.
                    let model_of = |a: &Vec<u64>, b: &Vec<u64>| Model::Extension {
                        p: BigUint::from(p),
                        k,
                        q: BigUint::from(q),
                        field_id: String::new(),
                        modulus: f.modulus.iter().map(|&c| BigUint::from(c)).collect(),
                        a: a.iter().map(|&c| BigUint::from(c)).collect(),
                        b: b.iter().map(|&c| BigUint::from(c)).collect(),
                    };
                    let cert = |m: u64| {
                        let m = BigUint::from(m);
                        let fac = crate::cryptanalysis::curve_traits::arith::factor(&m, 1 << 12);
                        crate::cryptanalysis::curve_traits::certify::point_order_certificate(
                            &model_of(a, b),
                            &m,
                            &fac,
                        )
                    };
                    let (status, _) = cert(n).unwrap().unwrap();
                    assert!(status <= crate::cryptanalysis::curve_traits::Status::Bounded);
                    certified[usize::from(
                        status == crate::cryptanalysis::curve_traits::Status::Proved,
                    )] += 1;
                    for wrong in [n + 1, n + 2, n - 1] {
                        if let Ok(Some((status, _))) = cert(wrong) {
                            assert_eq!(
                                status,
                                crate::cryptanalysis::curve_traits::Status::Bounded,
                                "certified a wrong order {wrong} ≠ {n}"
                            );
                        }
                    }
                    let t = q as i64 + 1 - n as i64;
                    if t % p as i64 == 0 {
                        continue;
                    }
                    let model = Model::Extension {
                        p: BigUint::from(p),
                        k,
                        q: BigUint::from(q),
                        field_id: String::new(),
                        modulus: f.modulus.iter().map(|&c| BigUint::from(c)).collect(),
                        a: a.iter().map(|&c| BigUint::from(c)).collect(),
                        b: b.iter().map(|&c| BigUint::from(c)).collect(),
                    };
                    for ell in [2u64, 3, 5, 7] {
                        let l = ell as i64;
                        let delta = t * t - 4 * q as i64;
                        if ell == p
                            || delta % (l * l) != 0
                            || eigen_order(ell, &BigInt::from(t)) != 1
                        {
                            continue;
                        }
                        let torsion = 1 + points
                            .iter()
                            .filter(|(x, y)| kills(&f, a, x, y, &BigUint::from(ell)))
                            .count() as u64;
                        let scalar = frobenius_is_scalar(&model, &BigInt::from(t), ell).unwrap();
                        assert_eq!(
                            scalar,
                            torsion == ell * ell,
                            "GF({p}^{k}) a={a:?} b={b:?} ℓ={ell}"
                        );
                        checked[usize::from(scalar)] += 1;
                    }
                }
            }
        }
        assert!(checked[0] > 0 && checked[1] > 0, "{checked:?}");
        assert!(certified[1] > 0, "{certified:?}");
    }

    /// `E[n] ⊆ E(GF(q))` by brute force, over `f` (`GF(p)` as `FpK` with
    /// `k = 1`, or `GF(p^k)`).
    fn full_torsion(f: &FpK, els: &[Vec<u64>], a: &Vec<u64>, b: &Vec<u64>, n: u64) -> bool {
        let mut count = 1u64;
        for x in els {
            let rhs = f.add(&f.mul(&f.add(&f.mul(x, x), a), x), b);
            for y in els {
                if f.mul(y, y) == rhs && kills(f, a, x, y, &BigUint::from(n)) {
                    count += 1;
                }
            }
        }
        count == n * n
    }

    /// Frobenius is a scalar on `E[4]` exactly when `E[4]` is rational on
    /// `E` or on its quadratic twist (the units mod 4 are `±1`); for
    /// `E[8]` and `E[9]` a rational full torsion on either side forces
    /// "scalar".  Checked over GF(p) and GF(p²).
    #[test]
    fn the_prime_power_scalar_test_agrees_with_rational_torsion() {
        let mut seen = [0usize; 2];
        let fields: Vec<(u64, Vec<u64>)> = vec![
            (13, vec![0]),
            (17, vec![0]),
            (29, vec![0]),
            (37, vec![0]),
            (41, vec![0]),
            (5, vec![3, 0]),
            (7, vec![1, 0]),
        ];
        for (p, low) in fields {
            let f = FpK { p, modulus: low };
            let k = f.k() as u32;
            let q = p.pow(k);
            let els: Vec<Vec<u64>> = (0..q).map(|i| f.nth(i)).collect();
            // A non-square c: the twist is y² = x³ + c²a·x + c³b.
            let half = BigUint::from((q - 1) / 2);
            let c = els
                .iter()
                .find(|e| f.pow(e, &half) == f.neg(&f.one()))
                .unwrap()
                .clone();
            let (c2, c3) = (f.mul(&c, &c), f.mul(&c, &f.mul(&c, &c)));
            let wide = if k == 1 { 12 } else { 25 };
            for a in els.iter().take(wide) {
                for b in els.iter().skip(1).take(wide) {
                    let disc = f.add(
                        &f.mul(&f.int(4), &f.mul(a, &f.mul(a, a))),
                        &f.mul(&f.int(27), &f.mul(b, b)),
                    );
                    if f.is_zero(&disc) {
                        continue;
                    }
                    let mut count = 1i64;
                    for x in &els {
                        let rhs = f.add(&f.mul(&f.add(&f.mul(x, x), a), x), b);
                        count += els.iter().filter(|y| f.mul(y, y) == rhs).count() as i64;
                    }
                    let t = q as i64 + 1 - count;
                    if t % p as i64 == 0 {
                        continue;
                    }
                    let (ta, tb) = (f.mul(&c2, a), f.mul(&c3, b));
                    let model = if k == 1 {
                        Model::Prime {
                            p: BigUint::from(p),
                            a: BigUint::from(a[0]),
                            b: BigUint::from(b[0]),
                        }
                    } else {
                        Model::Extension {
                            p: BigUint::from(p),
                            k,
                            q: BigUint::from(q),
                            field_id: String::new(),
                            modulus: f.modulus.iter().map(|&c| BigUint::from(c)).collect(),
                            a: a.iter().map(|&c| BigUint::from(c)).collect(),
                            b: b.iter().map(|&c| BigUint::from(c)).collect(),
                        }
                    };
                    let tt = BigInt::from(t);
                    for (ell, j) in [(2u64, 2u32), (2, 3), (3, 2)] {
                        let n = ell.pow(j);
                        if (q as i64 + 1 - t) % (n * n) as i64 != 0
                            && (q as i64 + 1 + t) % (n * n) as i64 != 0
                        {
                            // Neither side can hold all of E[n]; test anyway
                            // for n = 4, where the equivalence is exact.
                            if n != 4 {
                                continue;
                            }
                        }
                        let Some(scalar) = scalar_on_power(&model, &tt, ell, j) else {
                            continue;
                        };
                        let rational =
                            full_torsion(&f, &els, a, b, n) || full_torsion(&f, &els, &ta, &tb, n);
                        if n == 4 {
                            assert_eq!(scalar, rational, "GF({p}^{k}) a={a:?} b={b:?} n={n}");
                        } else if rational {
                            assert!(scalar, "GF({p}^{k}) a={a:?} b={b:?} n={n}");
                        }
                        seen[usize::from(scalar)] += 1;
                    }
                }
            }
        }
        assert!(seen[0] > 0 && seen[1] > 0, "{seen:?}");
    }

    /// The same comparison over `GF(2^n)`, through the binary division
    /// polynomials.
    #[test]
    fn the_binary_scalar_test_agrees_with_counting_rational_torsion() {
        use crate::binary_ecc::curve::{scalar_mul, BinaryCurve, BinaryPoint};
        let mut checked = [0usize; 2];
        for (n, low) in [(7u32, vec![0u32, 1]), (8, vec![0, 2, 3, 4])] {
            let irr = IrreduciblePoly {
                degree: n,
                low_terms: low,
            };
            let gf = Gf2::new(&irr);
            let q = 1u64 << n;
            for a in 0..2u64 {
                for b in 1..q {
                    let count = super::super::certify::count_binary(&gf, a, b);
                    let t = q as i64 + 1 - count as i64;
                    let el = |v: u64| F2mElement::from_biguint(&BigUint::from(v), n);
                    let curve = BinaryCurve {
                        m: n,
                        irreducible: irr.clone(),
                        a: el(a),
                        b: el(b),
                        generator: BinaryPoint::Infinity,
                        order: BigUint::from(count),
                        cofactor: BigUint::one(),
                    };
                    for ell in [3u64, 5, 7] {
                        let l = ell as i64;
                        let delta = t * t - 4 * q as i64;
                        if delta % (l * l) != 0 || eigen_order(ell, &BigInt::from(t)) != 1 {
                            continue;
                        }
                        let mut torsion = 1;
                        for x in 0..q {
                            for y in 0..q {
                                let pt = BinaryPoint::Affine { x: el(x), y: el(y) };
                                if curve.is_on_curve(&pt)
                                    && scalar_mul(&curve, &pt, &BigUint::from(ell))
                                        == BinaryPoint::Infinity
                                {
                                    torsion += 1;
                                }
                            }
                        }
                        let model = Model::Binary {
                            n,
                            modulus: BigUint::zero(),
                            irr: irr.clone(),
                            a: el(a),
                            b: el(b),
                        };
                        let scalar = frobenius_is_scalar(&model, &BigInt::from(t), ell).unwrap();
                        assert_eq!(scalar, torsion == ell * ell, "n={n} a={a} b={b} ℓ={ell}");
                        checked[usize::from(scalar)] += 1;
                    }
                }
            }
        }
        assert!(checked[0] > 0 && checked[1] > 0, "{checked:?}");
    }
}
