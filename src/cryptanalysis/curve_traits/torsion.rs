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
    fn is_zero(&self, a: &Self::E) -> bool {
        *a == self.zero()
    }
    fn sub(&self, a: &Self::E, b: &Self::E) -> Self::E {
        self.add(a, &self.neg(b))
    }
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
}

/// A polynomial, coefficients low to high, no trailing zeros.
type Poly<E> = Vec<E>;

fn trim<F: Field>(f: &F, mut p: Poly<F::E>) -> Poly<F::E> {
    while p.last().is_some_and(|c| f.is_zero(c)) {
        p.pop();
    }
    p
}

fn padd<F: Field>(f: &F, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
    let n = a.len().max(b.len());
    let z = f.zero();
    let out = (0..n)
        .map(|i| f.add(a.get(i).unwrap_or(&z), b.get(i).unwrap_or(&z)))
        .collect();
    trim(f, out)
}

fn psub<F: Field>(f: &F, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
    let nb: Poly<F::E> = b.iter().map(|c| f.neg(c)).collect();
    padd(f, a, &nb)
}

fn pmul<F: Field>(f: &F, a: &Poly<F::E>, b: &Poly<F::E>) -> Poly<F::E> {
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

fn pscale<F: Field>(f: &F, a: &Poly<F::E>, c: &F::E) -> Poly<F::E> {
    trim(f, a.iter().map(|x| f.mul(x, c)).collect())
}

/// `a mod m` for monic `m`.
fn prem<F: Field>(f: &F, a: &Poly<F::E>, m: &Poly<F::E>) -> Poly<F::E> {
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

fn monic<F: Field>(f: &F, a: &Poly<F::E>) -> Poly<F::E> {
    let lead = f.inv(a.last().expect("non-zero polynomial"));
    pscale(f, a, &lead)
}

/// `base^e mod m`.
fn ppow<F: Field>(f: &F, base: &Poly<F::E>, e: &BigUint, m: &Poly<F::E>) -> Poly<F::E> {
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
fn compose<F: Field>(f: &F, g: &Poly<F::E>, h: &Poly<F::E>, m: &Poly<F::E>) -> Poly<F::E> {
    let mut acc: Poly<F::E> = Vec::new();
    for c in g.iter().rev() {
        acc = padd(f, &prem(f, &pmul(f, &acc, h), m), &trim(f, vec![c.clone()]));
    }
    acc
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

/// The order `k` of the double eigenvalue `λ` of Frobenius mod `ℓ`; 0
/// when `λ ≡ 0`, which only happens for `ℓ` the characteristic.
fn eigen_order(ell: u64, trace: &BigInt) -> u32 {
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
    let limbs = q.bits().div_ceil(64);
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
        Model::Other { .. } | Model::Extension { .. } => None,
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
