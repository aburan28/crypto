//! Variable-time Jacobian-coordinate arithmetic on GOST R 34.10-2001
//! **CryptoPro-B** (RFC 4357 §11.4.4) and the two scalar-multiplication
//! arms of the GLV measurement in
//! `research/cryptopro_b_glv_chain_20261005/`.
//!
//! **This is a measurement harness, not a production implementation.**
//! Every routine here branches on point and scalar values (identity
//! checks, wNAF digit scans, the `H = 0` test in the addition law, the
//! `num-bigint` decomposition), so its running time depends on its
//! inputs.  Never use it with a secret scalar; the constant-time code for
//! this crate's curves lives in [`crate::ecc::p256_point`] and
//! [`crate::ecc::secp256k1_point`].
//!
//! Coordinates are Jacobian, `x = X/Z²`, `y = Y/Z³`, with the `a = -3`
//! formulas of the Explicit-Formulas Database: `dbl-2001-b` (3M + 5S),
//! `add-2007-bl` (11M + 5S) and the mixed addition `madd-2007-bl`
//! (7M + 4S) for an affine second operand.  (The P-256 module of this
//! crate uses the complete projective Renes–Costello–Batina formulas
//! instead; the Jacobian set here is what both arms of the measurement
//! share, so the comparison is between scalar-multiplication strategies,
//! not between point-arithmetic styles.)
//!
//! The two arms:
//!
//! * **baseline** — [`scalar_mul_wnaf`]: width-`w` NAF with an affine
//!   table of the odd multiples `P, 3P, …, (2^(w-1) - 1)·P` built with one
//!   batched inversion, and mixed additions in the main loop.
//! * **GLV-2** — [`scalar_mul_glv`]: `k = k1 + k2·λ (mod n)` by Babai
//!   rounding against the frozen LLL-reduced basis, `φ(P)` evaluated as
//!   the degree-175 isogeny chain `5·5·7` of
//!   [`crate::ecc::cryptopro_b_chain_consts::CHAIN_5_5_7`] (two 5-isogenies
//!   and a 7-isogeny through neighbouring curves, then the isomorphism
//!   back to `E`), and an interleaved two-scalar width-`w` NAF with shared
//!   doublings and one batched inversion for both tables.
//!
//! The chain evaluation is generic over [`Chain`]: the constants module
//! carries both the degree-175 chain `5·5·7` (the element `4 + ω`) and the
//! degree-155 chain `5·31` (the element `ω` itself, a 5-isogeny then a
//! 31-isogeny), and each is a GLV arm of the measurement with its own
//! `λ`, reduced basis and test vectors.

use super::cryptopro_b_chain_consts::{Chain, ChainStep, A, B, N};
use super::cryptopro_b_field::{CryptoProBFieldElement, P};
use super::field::FieldElement;
use super::point::Point;
use crate::ct_bignum::{Uint, U256};
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::Zero;

type Fe = CryptoProBFieldElement;

/// An affine point `(x, y)` on CryptoPro-B.  Never the identity.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct CryptoProBAffine {
    pub x: Fe,
    pub y: Fe,
}

/// A Jacobian point `(X : Y : Z)`, `x = X/Z²`, `y = Y/Z³`.  The identity
/// is any point with `Z = 0`.
#[derive(Copy, Clone, Debug)]
pub struct CryptoProBJacobian {
    pub x: Fe,
    pub y: Fe,
    pub z: Fe,
}

type Aff = CryptoProBAffine;
type Jac = CryptoProBJacobian;

/// The curve coefficient `a = p - 3` in Montgomery form.
fn a_fe() -> Fe {
    Fe::from_limbs(&A)
}

/// The curve coefficient `b` in Montgomery form.
fn b_fe() -> Fe {
    Fe::from_limbs(&B)
}

impl CryptoProBAffine {
    /// From canonical little-endian limbs (the constants-module encoding).
    pub fn from_limbs(x: &[u64; 4], y: &[u64; 4]) -> Self {
        CryptoProBAffine {
            x: Fe::from_limbs(x),
            y: Fe::from_limbs(y),
        }
    }

    pub fn neg(&self) -> Self {
        CryptoProBAffine {
            x: self.x,
            y: self.y.neg(),
        }
    }

    /// `y² == x³ - 3x + b` on CryptoPro-B.
    pub fn is_on_curve(&self) -> bool {
        self.is_on_curve_with(&a_fe(), &b_fe())
    }

    /// `y² == x³ + a·x + b` for arbitrary Weierstrass coefficients (the
    /// intermediate curves of an isogeny chain).
    pub fn is_on_curve_with(&self, a: &Fe, b: &Fe) -> bool {
        let lhs = self.y.sqr();
        let rhs = self.x.sqr().mul(&self.x).add(&a.mul(&self.x)).add(b);
        lhs == rhs
    }

    /// `None` for the point at infinity.
    pub fn from_textbook(p: &Point) -> Option<Self> {
        match p {
            Point::Infinity => None,
            Point::Affine { x, y } => Some(CryptoProBAffine {
                x: Fe::from_biguint(&x.value),
                y: Fe::from_biguint(&y.value),
            }),
        }
    }

    pub fn to_textbook(&self) -> Point {
        let p = P.to_biguint();
        Point::Affine {
            x: FieldElement::new(self.x.to_biguint(), p.clone()),
            y: FieldElement::new(self.y.to_biguint(), p),
        }
    }
}

impl CryptoProBJacobian {
    pub const IDENTITY: CryptoProBJacobian = CryptoProBJacobian {
        x: Fe::ONE,
        y: Fe::ONE,
        z: Fe::ZERO,
    };

    pub fn from_affine(p: &Aff) -> Self {
        CryptoProBJacobian {
            x: p.x,
            y: p.y,
            z: Fe::ONE,
        }
    }

    pub fn is_identity(&self) -> bool {
        self.z.is_zero()
    }

    /// One inversion.  `None` for the identity.
    pub fn to_affine(&self) -> Option<Aff> {
        if self.is_identity() {
            return None;
        }
        let zi = self.z.inv();
        let zi2 = zi.sqr();
        let zi3 = zi2.mul(&zi);
        Some(CryptoProBAffine {
            x: self.x.mul(&zi2),
            y: self.y.mul(&zi3),
        })
    }

    pub fn to_textbook(&self) -> Point {
        match self.to_affine() {
            None => Point::Infinity,
            Some(a) => a.to_textbook(),
        }
    }

    pub fn neg(&self) -> Self {
        CryptoProBJacobian {
            x: self.x,
            y: self.y.neg(),
            z: self.z,
        }
    }

    /// Equality as curve points, by cross-multiplication (no inversion):
    /// `X1·Z2² == X2·Z1²` and `Y1·Z2³ == Y2·Z1³`.
    pub fn eq_point(&self, other: &Self) -> bool {
        match (self.is_identity(), other.is_identity()) {
            (true, true) => true,
            (true, false) | (false, true) => false,
            (false, false) => {
                let z1z1 = self.z.sqr();
                let z2z2 = other.z.sqr();
                let u1 = self.x.mul(&z2z2);
                let u2 = other.x.mul(&z1z1);
                let s1 = self.y.mul(&other.z).mul(&z2z2);
                let s2 = other.y.mul(&self.z).mul(&z1z1);
                u1 == u2 && s1 == s2
            }
        }
    }

    /// EFD `dbl-2001-b` for `a = -3`: 3M + 5S.
    ///
    /// `Y = 0` (a point of order two) gives `Z3 = 0`; CryptoPro-B has
    /// prime odd order so no such point exists on it.
    pub fn double(&self) -> Self {
        if self.is_identity() {
            return *self;
        }
        let delta = self.z.sqr(); // S
        let gamma = self.y.sqr(); // S
        let beta = self.x.mul(&gamma); // M
        let t = self.x.sub(&delta).mul(&self.x.add(&delta)); // M
        let alpha = t.add(&t).add(&t); // 3(X1-δ)(X1+δ)
        let beta2 = beta.add(&beta);
        let beta4 = beta2.add(&beta2);
        let beta8 = beta4.add(&beta4);
        let x3 = alpha.sqr().sub(&beta8); // S
        let z3 = self.y.add(&self.z).sqr().sub(&gamma).sub(&delta); // S
        let gamma2 = gamma.sqr(); // S
        let g2 = gamma2.add(&gamma2);
        let g4 = g2.add(&g2);
        let g8 = g4.add(&g4);
        let y3 = alpha.mul(&beta4.sub(&x3)).sub(&g8); // M
        CryptoProBJacobian {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// EFD `add-2007-bl`: 11M + 5S, with the identity and the `P = ±Q`
    /// cases handled by branches.
    pub fn add(&self, other: &Self) -> Self {
        if self.is_identity() {
            return *other;
        }
        if other.is_identity() {
            return *self;
        }
        let z1z1 = self.z.sqr(); // S
        let z2z2 = other.z.sqr(); // S
        let u1 = self.x.mul(&z2z2); // M
        let u2 = other.x.mul(&z1z1); // M
        let s1 = self.y.mul(&other.z).mul(&z2z2); // 2M
        let s2 = other.y.mul(&self.z).mul(&z1z1); // 2M
        let h = u2.sub(&u1);
        let r = s2.sub(&s1);
        let r = r.add(&r); // 2(S2-S1)
        if h.is_zero() {
            return if r.is_zero() {
                self.double()
            } else {
                CryptoProBJacobian::IDENTITY
            };
        }
        let h2 = h.add(&h);
        let i = h2.sqr(); // S  ((2H)²)
        let j = h.mul(&i); // M
        let v = u1.mul(&i); // M
        let x3 = r.sqr().sub(&j).sub(&v).sub(&v); // S
        let s1j = s1.mul(&j); // M
        let y3 = r.mul(&v.sub(&x3)).sub(&s1j).sub(&s1j); // M
        let z3 = self
            .z
            .add(&other.z)
            .sqr() // S
            .sub(&z1z1)
            .sub(&z2z2)
            .mul(&h); // M
        CryptoProBJacobian {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// EFD `madd-2007-bl`, mixed addition with an affine second operand
    /// (`Z2 = 1`): 7M + 4S, same special cases as [`Self::add`].
    pub fn add_affine(&self, other: &Aff) -> Self {
        if self.is_identity() {
            return CryptoProBJacobian::from_affine(other);
        }
        let z1z1 = self.z.sqr(); // S
        let u2 = other.x.mul(&z1z1); // M
        let s2 = other.y.mul(&self.z).mul(&z1z1); // 2M
        let h = u2.sub(&self.x);
        let r = s2.sub(&self.y);
        let r = r.add(&r); // 2(S2-Y1)
        if h.is_zero() {
            return if r.is_zero() {
                self.double()
            } else {
                CryptoProBJacobian::IDENTITY
            };
        }
        let hh = h.sqr(); // S
        let i = hh.add(&hh);
        let i = i.add(&i); // 4HH
        let j = h.mul(&i); // M
        let v = self.x.mul(&i); // M
        let x3 = r.sqr().sub(&j).sub(&v).sub(&v); // S
        let y1j = self.y.mul(&j); // M
        let y3 = r.mul(&v.sub(&x3)).sub(&y1j).sub(&y1j); // M
        let z3 = self.z.add(&h).sqr().sub(&z1z1).sub(&hh); // S
        CryptoProBJacobian {
            x: x3,
            y: y3,
            z: z3,
        }
    }
}

/// Affine forms of several Jacobian points with one inversion
/// (Montgomery's trick): `3(n-1)` multiplications plus one inversion for
/// the `Z` inverses, then `3M + 1S` per point.
///
/// Panics on the identity, which has no affine form; the callers here
/// only pass odd multiples of a point of prime order, and `φ(P)`.
pub fn batch_to_affine(points: &[Jac]) -> Vec<Aff> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    let mut prefix = Vec::with_capacity(n); // prefix[i] = z_0 · … · z_i
    let mut acc = Fe::ONE;
    for p in points {
        assert!(
            !p.is_identity(),
            "batch_to_affine: the identity has no affine form"
        );
        acc = acc.mul(&p.z);
        prefix.push(acc);
    }
    let mut inv = acc.inv(); // 1 / (z_0 · … · z_{n-1})
    let mut out = vec![
        CryptoProBAffine {
            x: Fe::ZERO,
            y: Fe::ZERO,
        };
        n
    ];
    for i in (0..n).rev() {
        let zi = if i == 0 { inv } else { inv.mul(&prefix[i - 1]) };
        if i > 0 {
            inv = inv.mul(&points[i].z);
        }
        let zi2 = zi.sqr();
        let zi3 = zi2.mul(&zi);
        out[i] = CryptoProBAffine {
            x: points[i].x.mul(&zi2),
            y: points[i].y.mul(&zi3),
        };
    }
    out
}

// ── width-w NAF ────────────────────────────────────────────────────────────

/// Width-`w` non-adjacent form of `k`, least-significant digit first.
/// Digits are zero or odd with `|d| < 2^(w-1)`, and any `w` consecutive
/// digits hold at most one non-zero.  `2 <= w <= 7`.
pub fn wnaf_digits(k: &U256, w: u32) -> Vec<i8> {
    assert!((2..=7).contains(&w), "wNAF width must be in 2..=7");
    // Five limbs: `k + (2^w - m)` can carry past 2^256.
    let mut k = [k.0[0], k.0[1], k.0[2], k.0[3], 0u64];
    let modulus = 1u64 << w;
    let half = 1u64 << (w - 1);
    let mut digits = Vec::with_capacity(258);
    while k.iter().any(|&limb| limb != 0) {
        let digit = if k[0] & 1 == 1 {
            let m = k[0] & (modulus - 1);
            if m >= half {
                // digit = m - 2^w < 0, so k -= digit adds 2^w - m
                add_small(&mut k, modulus - m);
                (m as i64 - modulus as i64) as i8
            } else {
                sub_small(&mut k, m);
                m as i8
            }
        } else {
            0
        };
        digits.push(digit);
        shr1(&mut k);
    }
    digits
}

fn add_small(k: &mut [u64; 5], v: u64) {
    let mut carry = v;
    for limb in k.iter_mut() {
        let (sum, c) = limb.overflowing_add(carry);
        *limb = sum;
        carry = u64::from(c);
        if carry == 0 {
            break;
        }
    }
}

fn sub_small(k: &mut [u64; 5], v: u64) {
    let mut borrow = v;
    for limb in k.iter_mut() {
        let (diff, b) = limb.overflowing_sub(borrow);
        *limb = diff;
        borrow = u64::from(b);
        if borrow == 0 {
            break;
        }
    }
}

fn shr1(k: &mut [u64; 5]) {
    for i in 0..4 {
        k[i] = (k[i] >> 1) | (k[i + 1] << 63);
    }
    k[4] >>= 1;
}

/// `P, 3P, 5P, …, (2^(w-1) - 1)·P` in Jacobian coordinates: one doubling
/// and `2^(w-2) - 1` full additions.
fn odd_multiples(p: &Jac, w: u32) -> Vec<Jac> {
    let count = 1usize << (w - 2);
    let mut table = Vec::with_capacity(count);
    table.push(*p);
    if count > 1 {
        let p2 = p.double();
        for i in 1..count {
            let next = table[i - 1].add(&p2);
            table.push(next);
        }
    }
    table
}

/// The table point for a non-zero wNAF digit, negated when the digit is
/// negative, and negated again when the scalar itself is negative.
#[inline]
fn table_entry(table: &[Aff], digit: i8, negate: bool) -> Aff {
    let entry = table[(digit.unsigned_abs() >> 1) as usize];
    if (digit < 0) != negate {
        entry.neg()
    } else {
        entry
    }
}

/// **Baseline arm**: `k·P` by width-`w` NAF.  Builds the affine table of
/// odd multiples with one batched inversion, then one doubling per digit
/// and one mixed addition per non-zero digit.  `k < 2^256`; the result is
/// Jacobian (no final inversion).
pub fn scalar_mul_wnaf(p: &Aff, k: &BigUint, w: u32) -> Jac {
    assert!(k.bits() <= 256, "scalar must fit in 256 bits");
    let digits = wnaf_digits(&U256::from_biguint(k), w);
    let table = batch_to_affine(&odd_multiples(&Jac::from_affine(p), w));
    let mut q = Jac::IDENTITY;
    for &d in digits.iter().rev() {
        if !q.is_identity() {
            q = q.double();
        }
        if d != 0 {
            q = q.add_affine(&table_entry(&table, d, false));
        }
    }
    q
}

// ── the chain endomorphism ─────────────────────────────────────────────────

/// One isogeny step with its polynomials in Montgomery form.
struct PreparedStep {
    ell: u32,
    /// Monic kernel polynomial, degree `(ell-1)/2`, low first.
    psi: Vec<Fe>,
    /// Numerator of the x-map, degree `ell`, low first.
    n: Vec<Fe>,
    /// Numerator of the y-map, degree `(3·ell-3)/2`, low first.
    m: Vec<Fe>,
    codomain_a: Fe,
    codomain_b: Fe,
}

/// A [`Chain`] with every field constant converted to Montgomery form
/// once, so that evaluating `φ` costs no conversions.  Construction
/// validates the frozen input: every coefficient is a canonical value
/// below `p`, each polynomial has the degree the step's `ell` dictates,
/// the kernel polynomial is monic, and the steps compose (the codomain
/// of step `i` is the domain of step `i+1`, the first domain is `E`).
pub struct PreparedChain {
    pub name: &'static str,
    pub degree: u64,
    steps: Vec<PreparedStep>,
    /// `u²` and `u³` of the final isomorphism `(x, y) -> (u²x, u³y)`.
    u2: Fe,
    u3: Fe,
}

fn checked_fe(limbs: &[u64; 4], what: &str) -> Fe {
    assert!(
        bool::from(Uint(*limbs).ct_lt(&P)),
        "chain constant {what} is not a canonical value below p"
    );
    Fe::from_limbs(limbs)
}

fn prepared_poly(coeffs: &[[u64; 4]], want_degree: usize, what: &str) -> Vec<Fe> {
    assert_eq!(
        coeffs.len(),
        want_degree + 1,
        "chain polynomial {what} does not have degree {want_degree}"
    );
    coeffs
        .iter()
        .map(|c| checked_fe(c, what))
        .collect::<Vec<Fe>>()
}

impl PreparedChain {
    pub fn new(chain: &Chain) -> Self {
        assert!(!chain.steps.is_empty(), "a chain has at least one step");
        let mut steps = Vec::with_capacity(chain.steps.len());
        let mut expect_a = A;
        let mut expect_b = B;
        let mut degree = 1u64;
        for (idx, step) in chain.steps.iter().enumerate() {
            let ChainStep { ell, .. } = *step;
            assert!(
                ell >= 3 && ell % 2 == 1,
                "step {idx}: ell must be an odd prime"
            );
            assert_eq!(
                step.domain_a, expect_a,
                "step {idx}: domain a does not compose"
            );
            assert_eq!(
                step.domain_b, expect_b,
                "step {idx}: domain b does not compose"
            );
            let ell_us = ell as usize;
            let psi = prepared_poly(step.psi, (ell_us - 1) / 2, "psi");
            let n = prepared_poly(step.n, ell_us, "N");
            let m = prepared_poly(step.m, (3 * ell_us - 3) / 2, "M");
            assert_eq!(
                psi.last().copied(),
                Some(Fe::ONE),
                "step {idx}: kernel polynomial is not monic"
            );
            steps.push(PreparedStep {
                ell,
                psi,
                n,
                m,
                codomain_a: checked_fe(&step.codomain_a, "codomain a"),
                codomain_b: checked_fe(&step.codomain_b, "codomain b"),
            });
            expect_a = step.codomain_a;
            expect_b = step.codomain_b;
            degree *= u64::from(ell);
        }
        assert_eq!(
            degree, chain.degree,
            "chain degree is not the product of the steps"
        );
        let u = checked_fe(&chain.iso_u, "iso_u");
        assert!(!u.is_zero(), "iso_u must be a unit");
        let u2 = u.sqr();
        let u3 = u2.mul(&u);
        // The isomorphism must land back on E: a = u^4·a_last, b = u^6·b_last.
        let last = steps.last().expect("non-empty");
        let u4 = u2.sqr();
        let u6 = u4.mul(&u2);
        assert_eq!(
            u4.mul(&last.codomain_a),
            a_fe(),
            "iso_u does not map a back to E"
        );
        assert_eq!(
            u6.mul(&last.codomain_b),
            b_fe(),
            "iso_u does not map b back to E"
        );
        PreparedChain {
            name: chain.name,
            degree: chain.degree,
            steps,
            u2,
            u3,
        }
    }

    pub fn step_count(&self) -> usize {
        self.steps.len()
    }

    /// `φ(P)` in Jacobian coordinates.
    ///
    /// Each step runs in projective `P²` coordinates `(X : Y : Z)`,
    /// `x = X/Z`, `y = Y/Z`, with the homogenised polynomials evaluated by
    /// Horner over precomputed powers of `Z`:
    /// `X' = N_h·ψ_h`, `Y' = Y·M_h`, `Z' = Z·ψ_h³`.  After the last step the
    /// isomorphism `X'' = u²X`, `Y'' = u³Y` lands on `E`, and the result
    /// is converted to Jacobian without an inversion as
    /// `(X·Z, Y·Z², Z)`.
    pub fn apply(&self, p: &Aff) -> Jac {
        let (mut x, mut y, mut z) = (p.x, p.y, Fe::ONE);
        for step in &self.steps {
            (x, y, z) = step_projective(step, &x, &y, &z);
        }
        x = x.mul(&self.u2);
        y = y.mul(&self.u3);
        let z2 = z.sqr();
        CryptoProBJacobian {
            x: x.mul(&z),
            y: y.mul(&z2),
            z,
        }
    }

    /// The affine image of `P` after each step, with that step's codomain
    /// coefficients: for checking that the chain walks through the curves
    /// the constants declare.  One inversion per step; tests only.
    #[cfg(test)]
    fn step_images(&self, p: &Aff) -> Vec<(Aff, Fe, Fe)> {
        let (mut x, mut y, mut z) = (p.x, p.y, Fe::ONE);
        let mut out = Vec::with_capacity(self.steps.len());
        for step in &self.steps {
            (x, y, z) = step_projective(step, &x, &y, &z);
            let zi = z.inv();
            out.push((
                CryptoProBAffine {
                    x: x.mul(&zi),
                    y: y.mul(&zi),
                },
                step.codomain_a,
                step.codomain_b,
            ));
        }
        out
    }
}

/// `Σ c_i · X^i · Z^(d-i)` for `d = coeffs.len() - 1`, Horner in `X` with
/// `zpow[j] = Z^j`.
#[inline]
fn horner_homogeneous(coeffs: &[Fe], x: &Fe, zpow: &[Fe]) -> Fe {
    let d = coeffs.len() - 1;
    let mut acc = coeffs[d];
    for i in (0..d).rev() {
        acc = acc.mul(x).add(&coeffs[i].mul(&zpow[d - i]));
    }
    acc
}

fn step_projective(step: &PreparedStep, x: &Fe, y: &Fe, z: &Fe) -> (Fe, Fe, Fe) {
    let max_deg = step.m.len() - 1; // (3·ell-3)/2 >= ell >= (ell-1)/2
    let mut zpow = Vec::with_capacity(max_deg + 1);
    zpow.push(Fe::ONE);
    zpow.push(*z);
    for i in 2..=max_deg {
        let next = zpow[i - 1].mul(z);
        zpow.push(next);
    }
    let psi = horner_homogeneous(&step.psi, x, &zpow);
    let num = horner_homogeneous(&step.n, x, &zpow);
    let mnum = horner_homogeneous(&step.m, x, &zpow);
    let psi3 = psi.sqr().mul(&psi);
    debug_assert!(step.ell >= 3);
    (num.mul(&psi), y.mul(&mnum), z.mul(&psi3))
}

// ── GLV decomposition and the GLV arm ──────────────────────────────────────

/// Everything the GLV arm needs that does not depend on the input: the
/// prepared chain, `λ`, `n`, the reduced lattice basis and its
/// determinant.  Build it once; it is curve-constant setup, like the
/// Montgomery-form curve coefficients of the other point modules.
pub struct GlvContext {
    pub chain: PreparedChain,
    /// Eigenvalue of `φ` on the prime-order group.
    pub lambda: BigUint,
    /// The group order `n`.
    pub order: BigUint,
    /// Rows `b1 = (b11, b12)`, `b2 = (b21, b22)` of the reduced basis of
    /// `{(v1, v2) : v1 + v2·λ ≡ 0 (mod n)}`.
    basis: [[BigInt; 2]; 2],
    det: BigInt,
    pub babai_bound_bits: u32,
}

/// Nearest integer to `num / den` (`den != 0`); exact halves round up.
fn round_div(num: &BigInt, den: &BigInt) -> BigInt {
    let (num, den) = if den.sign() == Sign::Minus {
        (-num, -den)
    } else {
        (num.clone(), den.clone())
    };
    (num * 2u8 + &den).div_floor(&(den * 2u8))
}

impl GlvContext {
    pub fn new(chain: &Chain) -> Self {
        let parse = |s: &str| -> BigInt {
            BigInt::parse_bytes(s.as_bytes(), 10)
                .unwrap_or_else(|| panic!("glv_basis entry {s:?} is not a signed decimal"))
        };
        let basis = [
            [parse(chain.glv_basis[0][0]), parse(chain.glv_basis[0][1])],
            [parse(chain.glv_basis[1][0]), parse(chain.glv_basis[1][1])],
        ];
        let det = &basis[0][0] * &basis[1][1] - &basis[0][1] * &basis[1][0];
        assert!(!det.is_zero(), "glv_basis is singular");
        GlvContext {
            chain: PreparedChain::new(chain),
            lambda: Uint(chain.lambda).to_biguint(),
            order: Uint(N).to_biguint(),
            basis,
            det,
            babai_bound_bits: chain.babai_bound_bits,
        }
    }

    pub fn basis(&self) -> &[[BigInt; 2]; 2] {
        &self.basis
    }

    /// Babai rounding: write `(k, 0) = c1·b1 + c2·b2` over the rationals
    /// (`c1 = k·b22/det`, `c2 = -k·b12/det`), round to `r1, r2`, and return
    /// `(k1, k2) = (k, 0) - r1·b1 - r2·b2`, so that `k1 + k2·λ ≡ k (mod n)`.
    pub fn decompose(&self, k: &BigUint) -> (BigInt, BigInt) {
        let k = BigInt::from(k.clone());
        let r1 = round_div(&(&k * &self.basis[1][1]), &self.det);
        let r2 = round_div(&(-&k * &self.basis[0][1]), &self.det);
        let k1 = &k - &r1 * &self.basis[0][0] - &r2 * &self.basis[1][0];
        let k2 = -(&r1 * &self.basis[0][1]) - &r2 * &self.basis[1][1];
        (k1, k2)
    }

    /// `k1 + k2·λ ≡ k (mod n)` and both parts within `babai_bound_bits`.
    pub fn decomposition_is_valid(&self, k: &BigUint, k1: &BigInt, k2: &BigInt) -> bool {
        let n = BigInt::from(self.order.clone());
        let lambda = BigInt::from(self.lambda.clone());
        let lhs = (k1 + k2 * lambda).mod_floor(&n);
        let rhs = BigInt::from(k.clone()).mod_floor(&n);
        let bound = u64::from(self.babai_bound_bits);
        lhs == rhs && k1.bits() <= bound && k2.bits() <= bound
    }
}

fn split_sign(k: &BigInt) -> (bool, U256) {
    assert!(k.bits() <= 256, "GLV part must fit in 256 bits");
    (k.sign() == Sign::Minus, U256::from_biguint(k.magnitude()))
}

/// The interleaved two-scalar width-`w` NAF: `k1·P + k2·φ(P)` with shared
/// doublings, one table of odd multiples per point (both converted to
/// affine in one batched inversion), and a negated table entry for a
/// negative digit or a negative scalar.
pub fn glv_combine(p: &Aff, phi_p: &Jac, k1: &BigInt, k2: &BigInt, w: u32) -> Jac {
    let (neg1, m1) = split_sign(k1);
    let (neg2, m2) = split_sign(k2);
    let d1 = wnaf_digits(&m1, w);
    let d2 = wnaf_digits(&m2, w);
    let count = 1usize << (w - 2);
    let mut both = odd_multiples(&Jac::from_affine(p), w);
    both.extend(odd_multiples(phi_p, w));
    let tables = batch_to_affine(&both);
    let (t1, t2) = tables.split_at(count);
    let len = d1.len().max(d2.len());
    let mut q = Jac::IDENTITY;
    for i in (0..len).rev() {
        if !q.is_identity() {
            q = q.double();
        }
        let a = d1.get(i).copied().unwrap_or(0);
        if a != 0 {
            q = q.add_affine(&table_entry(t1, a, neg1));
        }
        let b = d2.get(i).copied().unwrap_or(0);
        if b != 0 {
            q = q.add_affine(&table_entry(t2, b, neg2));
        }
    }
    q
}

/// **GLV arm**: `k·P` as decomposition, `φ(P)` through the chain, then
/// [`glv_combine`].  Everything input-dependent is inside; only `ctx` is
/// prepared ahead of time.
pub fn scalar_mul_glv(ctx: &GlvContext, p: &Aff, k: &BigUint, w: u32) -> Jac {
    let (k1, k2) = ctx.decompose(k);
    let phi_p = ctx.chain.apply(p);
    glv_combine(p, &phi_p, &k1, &k2, w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::cryptopro_b_chain_consts::{CHAINS, CHAIN_5_31, CHAIN_5_5_7, N_HEX};
    use crate::ecc::curve::CurveParams;
    use num_bigint::RandBigInt;
    use num_traits::One;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn curve() -> CurveParams {
        CurveParams::gost_cryptopro_b()
    }

    fn generator() -> Aff {
        CryptoProBAffine::from_textbook(&curve().generator()).unwrap()
    }

    fn order() -> BigUint {
        curve().n
    }

    fn random_scalar(rng: &mut StdRng) -> BigUint {
        loop {
            let k = rng.gen_biguint_below(&order());
            if !k.is_zero() {
                return k;
            }
        }
    }

    /// A random point of the group, through the textbook implementation
    /// so that the Jacobian code under test plays no part in making it.
    fn random_point_textbook(rng: &mut StdRng) -> Aff {
        let c = curve();
        let r = random_scalar(rng);
        let p = c.generator().scalar_mul_vartime(&r, &c.a_fe());
        CryptoProBAffine::from_textbook(&p).unwrap()
    }

    /// A random point of the group through the baseline arm (fast; the
    /// baseline is itself checked against the textbook code elsewhere).
    fn random_point(rng: &mut StdRng) -> Aff {
        let r = random_scalar(rng);
        scalar_mul_wnaf(&generator(), &r, 5).to_affine().unwrap()
    }

    fn jac_eq_textbook(j: &Jac, p: &Point) -> bool {
        match (j.to_affine(), p) {
            (None, Point::Infinity) => true,
            (Some(a), Point::Affine { x, y }) => {
                a.x.to_biguint() == x.value && a.y.to_biguint() == y.value
            }
            _ => false,
        }
    }

    #[test]
    fn constants_match_curve_zoo() {
        let c = curve();
        assert_eq!(Uint(A).to_biguint(), c.a);
        assert_eq!(Uint(B).to_biguint(), c.b);
        assert_eq!(Uint(N).to_biguint(), c.n);
        assert_eq!(c.a, &c.p - 3u8);
        let n_hex = N_HEX.trim_start_matches("0x");
        assert_eq!(BigUint::parse_bytes(n_hex.as_bytes(), 16).unwrap(), c.n);
        assert_eq!(c.h, 1);
    }

    #[test]
    fn generator_is_on_curve_and_has_order_n() {
        let g = generator();
        assert!(g.is_on_curve());
        assert!(scalar_mul_wnaf(&g, &order(), 5).is_identity());
    }

    #[test]
    fn test_vector_points_are_on_curve() {
        for chain in CHAINS {
            assert!(!chain.vectors.is_empty(), "{}: no test vectors", chain.name);
            for v in chain.vectors {
                assert!(CryptoProBAffine::from_limbs(&v.p[0], &v.p[1]).is_on_curve());
                assert!(CryptoProBAffine::from_limbs(&v.phi_p[0], &v.phi_p[1]).is_on_curve());
                assert!(CryptoProBAffine::from_limbs(&v.k_p[0], &v.k_p[1]).is_on_curve());
            }
        }
    }

    #[test]
    fn double_and_add_match_textbook_for_random_points() {
        let c = curve();
        let a = c.a_fe();
        let mut rng = StdRng::seed_from_u64(0x6f5_7a1);
        for _ in 0..20 {
            let p = random_point_textbook(&mut rng);
            let q = random_point_textbook(&mut rng);
            let (tp, tq) = (p.to_textbook(), q.to_textbook());
            let (jp, jq) = (Jac::from_affine(&p), Jac::from_affine(&q));

            assert!(jac_eq_textbook(&jp.double(), &tp.double(&a)));
            assert!(jac_eq_textbook(&jp.add(&jq), &tp.add(&tq, &a)));
            assert!(jac_eq_textbook(&jp.add_affine(&q), &tp.add(&tq, &a)));
            // Mixed and full addition agree with Z != 1 on the left.
            let jp2 = jp.double();
            assert!(jp2.add_affine(&q).eq_point(&jp2.add(&jq)));
            assert!(jac_eq_textbook(&jp2.add(&jq), &tp.double(&a).add(&tq, &a)));
            // Adding through the formulas versus doubling.
            assert!(jp.add(&jp).eq_point(&jp.double()));
            assert!(jp.add_affine(&p).eq_point(&jp.double()));
            assert!(jp2.add(&jp2).eq_point(&jp2.double()));
            // Identity and inverses.
            assert!(jp.add(&Jac::IDENTITY).eq_point(&jp));
            assert!(Jac::IDENTITY.add(&jp).eq_point(&jp));
            assert!(Jac::IDENTITY.add_affine(&p).eq_point(&jp));
            assert!(jp.add(&jp.neg()).is_identity());
            assert!(jp.add_affine(&p.neg()).is_identity());
            assert!(jp2
                .add_affine(&jp2.to_affine().unwrap().neg())
                .is_identity());
            assert!(Jac::IDENTITY.double().is_identity());
            assert!(!jp.eq_point(&jq));
        }
    }

    #[test]
    fn scalar_mul_wnaf_matches_textbook() {
        let c = curve();
        let a = c.a_fe();
        let mut rng = StdRng::seed_from_u64(0x5ca1a);
        for w in [4u32, 5] {
            for _ in 0..10 {
                let p = random_point_textbook(&mut rng);
                let k = random_scalar(&mut rng);
                let want = p.to_textbook().scalar_mul(&k, &a);
                let got = scalar_mul_wnaf(&p, &k, w);
                assert!(jac_eq_textbook(&got, &want), "w={w} k={k}");
            }
        }
    }

    #[test]
    fn scalar_mul_wnaf_edge_cases() {
        let g = generator();
        let n = order();
        for w in 2..=7 {
            assert!(scalar_mul_wnaf(&g, &BigUint::zero(), w).is_identity());
            assert!(scalar_mul_wnaf(&g, &BigUint::one(), w).eq_point(&Jac::from_affine(&g)));
            assert!(scalar_mul_wnaf(&g, &BigUint::from(2u8), w)
                .eq_point(&Jac::from_affine(&g).double()));
            assert!(scalar_mul_wnaf(&g, &n, w).is_identity());
            assert!(scalar_mul_wnaf(&g, &(&n - 1u8), w).eq_point(&Jac::from_affine(&g).neg()));
            assert!(scalar_mul_wnaf(&g, &(&n + 1u8), w).eq_point(&Jac::from_affine(&g)));
            // Scalars that exercise every digit value of the width.
            for small in 1u32..70 {
                let want = scalar_mul_wnaf(&g, &BigUint::from(small), 2);
                assert!(scalar_mul_wnaf(&g, &BigUint::from(small), w).eq_point(&want));
            }
        }
    }

    #[test]
    fn wnaf_digits_reconstruct_the_scalar() {
        let mut rng = StdRng::seed_from_u64(0x3af);
        let mut samples: Vec<BigUint> = (0..50).map(|_| rng.gen_biguint(256)).collect();
        samples.push((BigUint::one() << 256) - 1u8);
        samples.push(BigUint::one() << 255);
        samples.push(BigUint::from(15u8));
        samples.push(BigUint::zero());
        for w in 2..=7u32 {
            let half = 1i64 << (w - 1);
            for k in &samples {
                let digits = wnaf_digits(&U256::from_biguint(k), w);
                let mut acc = BigInt::zero();
                for (i, &d) in digits.iter().enumerate() {
                    assert!(d == 0 || (d % 2 != 0 && (i64::from(d)).abs() < half));
                    acc += BigInt::from(d) << i;
                }
                assert_eq!(acc, BigInt::from(k.clone()), "w={w}");
                for window in digits.windows(w as usize) {
                    assert!(window.iter().filter(|&&d| d != 0).count() <= 1);
                }
                if let Some(&top) = digits.last() {
                    assert!(top > 0);
                }
            }
        }
    }

    #[test]
    fn batch_to_affine_matches_single_inversions() {
        let mut rng = StdRng::seed_from_u64(0xba7c4);
        let g = generator();
        let points: Vec<Jac> = (0..9)
            .map(|_| scalar_mul_wnaf(&g, &random_scalar(&mut rng), 5))
            .collect();
        let batch = batch_to_affine(&points);
        for (j, a) in points.iter().zip(&batch) {
            assert_eq!(j.to_affine().unwrap(), *a);
        }
        assert!(batch_to_affine(&[]).is_empty());
        assert_eq!(
            batch_to_affine(&points[..1])[0],
            points[0].to_affine().unwrap()
        );
    }

    #[test]
    fn chain_constants_prepare_and_walk_through_declared_curves() {
        let mut rng = StdRng::seed_from_u64(0xc4a1);
        for chain_consts in CHAINS {
            let chain = PreparedChain::new(chain_consts);
            assert_eq!(chain.step_count(), chain_consts.steps.len());
            assert_eq!(chain.degree, chain_consts.degree);
            for _ in 0..5 {
                let p = random_point(&mut rng);
                for (image, a, b) in chain.step_images(&p) {
                    assert!(image.is_on_curve_with(&a, &b));
                }
                assert!(chain.apply(&p).to_affine().unwrap().is_on_curve());
            }
        }
        assert_eq!(CHAIN_5_5_7.steps.len(), 3);
        assert_eq!(CHAIN_5_5_7.degree, 175);
        assert_eq!(CHAIN_5_31.steps.len(), 2);
        assert_eq!(CHAIN_5_31.degree, 155);
        assert_eq!(CHAIN_5_31.steps[1].ell, 31);
    }

    #[test]
    fn chain_matches_test_vectors() {
        for chain_consts in CHAINS {
            let chain = PreparedChain::new(chain_consts);
            for v in chain_consts.vectors {
                let p = CryptoProBAffine::from_limbs(&v.p[0], &v.p[1]);
                let want = CryptoProBAffine::from_limbs(&v.phi_p[0], &v.phi_p[1]);
                assert_eq!(
                    chain.apply(&p).to_affine().unwrap(),
                    want,
                    "{}: phi(P) differs from the test vector",
                    chain_consts.name
                );
            }
        }
    }

    #[test]
    fn lambda_is_a_root_of_the_characteristic_polynomial() {
        // ω = (1 + √-619)/2 has trace 1 and norm 155, so a + b·ω has trace
        // 2a + b and norm a² + ab + 155b².  The frozen JSON names each chain
        // by an element but records the element the composite was matched
        // to on points as `matched_element`: (-5, 1) = -5 + ω for the 5·5·7
        // chain (the named 4 + ω up to sign and conjugation, same norm 175)
        // and (0, -1) = -ω for the 5·31 chain.  λ is the eigenvalue of the
        // matched element, so λ² - trace·λ + norm ≡ 0 (mod n) with the
        // matched element's trace, and the norm is the chain degree.
        let n = BigInt::from(order());
        let mut eigen = Vec::new();
        for (chain, (a, b)) in [(&CHAIN_5_5_7, (-5i64, 1i64)), (&CHAIN_5_31, (0, -1))] {
            let trace = 2 * a + b;
            let norm = a * a + a * b + 155 * b * b;
            assert_eq!(chain.degree, u64::try_from(norm).unwrap());
            let lambda = BigInt::from(Uint(chain.lambda).to_biguint());
            let value = (&lambda * &lambda - trace * &lambda + norm).mod_floor(&n);
            assert!(value.is_zero(), "{}", chain.name);
            eigen.push(lambda);
        }
        // Both eigenvalues come from one embedding of ω: with λ(-5 + ω) and
        // λ(-ω) the sum is -5 (mod n).
        let sum = (&eigen[0] + &eigen[1] + BigInt::from(5)).mod_floor(&n);
        assert!(sum.is_zero());
    }

    #[test]
    fn glv_basis_rows_lie_in_the_lattice() {
        for chain in CHAINS {
            let ctx = GlvContext::new(chain);
            let n = BigInt::from(ctx.order.clone());
            let lambda = BigInt::from(ctx.lambda.clone());
            for row in ctx.basis() {
                assert!((&row[0] + &row[1] * &lambda).mod_floor(&n).is_zero());
                assert!(row[0].bits() <= 129 && row[1].bits() <= 129);
            }
        }
    }

    #[test]
    fn chain_is_multiplication_by_lambda() {
        let mut rng = StdRng::seed_from_u64(0xe16e);
        for chain in CHAINS {
            let ctx = GlvContext::new(chain);
            let mut points: Vec<Aff> = chain
                .vectors
                .iter()
                .map(|v| CryptoProBAffine::from_limbs(&v.p[0], &v.p[1]))
                .collect();
            points.extend((0..50).map(|_| random_point(&mut rng)));
            for p in &points {
                let phi = ctx.chain.apply(p);
                let lambda_p = scalar_mul_wnaf(p, &ctx.lambda, 5);
                assert!(phi.eq_point(&lambda_p), "{}", chain.name);
            }
        }
    }

    #[test]
    fn decomposition_matches_test_vectors() {
        for chain in CHAINS {
            let ctx = GlvContext::new(chain);
            for v in chain.vectors {
                let k = Uint(v.k).to_biguint();
                let (k1, k2) = ctx.decompose(&k);
                assert!(ctx.decomposition_is_valid(&k, &k1, &k2));
                let want1 = signed(v.k1_neg, &v.k1);
                let want2 = signed(v.k2_neg, &v.k2);
                assert_eq!(k1, want1, "{}", chain.name);
                assert_eq!(k2, want2, "{}", chain.name);
            }
        }
    }

    fn signed(neg: bool, mag: &[u64; 4]) -> BigInt {
        let m = BigInt::from(Uint(*mag).to_biguint());
        if neg {
            -m
        } else {
            m
        }
    }

    #[test]
    fn glv_matches_baseline_and_k_p_on_test_vectors() {
        for chain in CHAINS {
            let ctx = GlvContext::new(chain);
            for v in chain.vectors {
                let p = CryptoProBAffine::from_limbs(&v.p[0], &v.p[1]);
                let k = Uint(v.k).to_biguint();
                let want = CryptoProBAffine::from_limbs(&v.k_p[0], &v.k_p[1]);
                for w in [4u32, 5] {
                    let base = scalar_mul_wnaf(&p, &k, w);
                    let glv = scalar_mul_glv(&ctx, &p, &k, w);
                    assert_eq!(base.to_affine().unwrap(), want);
                    assert_eq!(glv.to_affine().unwrap(), want, "{} w={w}", chain.name);
                    assert!(glv.eq_point(&base));
                }
            }
        }
    }

    #[test]
    fn glv_matches_baseline_for_random_scalars() {
        let mut rng = StdRng::seed_from_u64(0x91_4e5);
        let g = generator();
        for chain in CHAINS {
            let ctx = GlvContext::new(chain);
            for i in 0..200 {
                // Half on the generator, half on random points.
                let p = if i % 2 == 0 {
                    g
                } else {
                    random_point(&mut rng)
                };
                let k = random_scalar(&mut rng);
                let (k1, k2) = ctx.decompose(&k);
                assert!(ctx.decomposition_is_valid(&k, &k1, &k2));
                assert!(k1.bits() <= u64::from(chain.babai_bound_bits));
                assert!(k2.bits() <= u64::from(chain.babai_bound_bits));
                let base = scalar_mul_wnaf(&p, &k, 5);
                let glv = scalar_mul_glv(&ctx, &p, &k, 5);
                assert!(glv.eq_point(&base), "{} k={k}", chain.name);
            }
        }
    }

    #[test]
    fn glv_handles_degenerate_scalars() {
        let g = generator();
        let n = order();
        for chain in CHAINS {
            let ctx = GlvContext::new(chain);
            for k in [
                BigUint::zero(),
                BigUint::one(),
                BigUint::from(2u8),
                ctx.lambda.clone(),
                &n - 1u8,
                &n - &ctx.lambda,
                (&ctx.lambda * 2u8) % &n,
            ] {
                let base = scalar_mul_wnaf(&g, &k, 5);
                let glv = scalar_mul_glv(&ctx, &g, &k, 5);
                assert!(glv.eq_point(&base), "{} k={k}", chain.name);
            }
            assert!(scalar_mul_glv(&ctx, &g, &BigUint::zero(), 5).is_identity());
        }
    }
}
