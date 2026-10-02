//! # Curves over small extension fields `F_{p^k}` as counted groups.
//!
//! The GLS twist lives over `F_{p²}` and a subfield curve's interesting
//! subgroup over `F_{p³}`; both need the same thing from the field —
//! arithmetic, the `p`-power Frobenius, square roots, and `F_p`
//! coordinates so that a line `c·F_p` can be walked and a summation
//! polynomial can be Weil-descended — and the same thing from the
//! curve: a [`CountedGroup`] so the framework's oracles, relation loop
//! and rho reference run on it unchanged.
//!
//! [`ExtField`] is that contract, [`Fp2`] (`s² = ν`) and [`Fp3`]
//! (`s³ = ν`, `p ≡ 1 (mod 3)`) implement it, and [`ExtCurve`] is the
//! short Weierstrass curve over any of them.  `p^k < 2^62` so that a
//! packed element fits the framework's `u64` keys.
//!
//! Companion to `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`.

use std::cell::Cell;

use serde::Serialize;

use crate::cryptanalysis::glv_invariant_base::{addm, mulm as mulm_raw, subm};
use crate::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crate::cryptanalysis::residual_walk::{
    inv_mod as inv_mod_raw, is_prime_u64, pow_mod, sqrt_mod,
};

thread_local! {
    /// `F_p` multiplications and inversions performed by this module's
    /// field arithmetic on the current thread: the unit every phase of
    /// an end-to-end `S` on an extension-field group is priced in
    /// (note §8, E8).  Square roots taken through [`sqrt_mod`] and
    /// [`pow_mod`] are not counted; the generic Tonelli–Shanks of
    /// [`ExtField::sqrt`] is, since it runs on the field's own `mul`.
    static FIELD_COUNTERS: Cell<(u64, u64)> = const { Cell::new((0, 0)) };
}

/// Return and reset the thread's `(multiplications, inversions)` in
/// `F_p` performed by [`ExtField`] arithmetic since the last call.
pub fn take_field_counters() -> (u64, u64) {
    FIELD_COUNTERS.with(|c| c.replace((0, 0)))
}

#[inline]
fn mulm(a: u64, b: u64, p: u64) -> u64 {
    FIELD_COUNTERS.with(|c| {
        let (m, i) = c.get();
        c.set((m + 1, i));
    });
    mulm_raw(a, b, p)
}

#[inline]
fn inv_mod(a: u64, p: u64) -> u64 {
    FIELD_COUNTERS.with(|c| {
        let (m, i) = c.get();
        c.set((m, i + 1));
    });
    inv_mod_raw(a, p)
}

/// Jacobi symbol `(a / n)` for odd `n`.
pub fn jacobi(mut a: u64, mut n: u64) -> i32 {
    let mut result = 1i32;
    a %= n;
    while a != 0 {
        while a.is_multiple_of(2) {
            a /= 2;
            if n % 8 == 3 || n % 8 == 5 {
                result = -result;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 {
            result = -result;
        }
        a %= n;
    }
    if n == 1 {
        result
    } else {
        0
    }
}

/// A finite field `F_{p^k}` presented over `F_p`.
// `from_fp` and `from_coords` take `&self` because the field carries the
// modulus and the presentation; they are constructors of elements, not
// of the field.
#[allow(clippy::wrong_self_convention)]
pub trait ExtField: Copy + std::fmt::Debug + Serialize {
    type El: Copy + PartialEq + Eq + std::fmt::Debug + Serialize;
    fn p(&self) -> u64;
    fn degree(&self) -> u32;
    /// `p^k`, which must fit a word.
    fn order(&self) -> u64 {
        (0..self.degree()).fold(1u64, |q, _| q * self.p())
    }
    fn zero(&self) -> Self::El;
    fn one(&self) -> Self::El;
    /// The generator `s` of the presentation, `F_{p^k} = F_p(s)`.
    fn s(&self) -> Self::El;
    fn from_fp(&self, a: u64) -> Self::El;
    /// The `F_p` coordinates in the power basis `1, s, …, s^{k−1}`.
    fn coords(&self, a: Self::El) -> Vec<u64>;
    fn from_coords(&self, c: &[u64]) -> Self::El;
    fn is_zero(&self, a: Self::El) -> bool;
    fn add(&self, a: Self::El, b: Self::El) -> Self::El;
    fn sub(&self, a: Self::El, b: Self::El) -> Self::El;
    fn neg(&self, a: Self::El) -> Self::El;
    fn mul(&self, a: Self::El, b: Self::El) -> Self::El;
    /// `k·a` for `k ∈ F_p`.
    fn scale(&self, a: Self::El, k: u64) -> Self::El;
    fn sqr(&self, a: Self::El) -> Self::El {
        self.mul(a, a)
    }
    fn inv(&self, a: Self::El) -> Option<Self::El>;
    /// `a^p`.
    fn frob(&self, a: Self::El) -> Self::El;
    /// A quadratic non-residue of the field.
    fn nonresidue(&self) -> Self::El;
    fn pow(&self, a: Self::El, e: u64) -> Self::El {
        let mut acc = self.one();
        let mut b = a;
        let mut e = e;
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul(acc, b);
            }
            b = self.sqr(b);
            e >>= 1;
        }
        acc
    }
    fn is_square(&self, a: Self::El) -> bool {
        self.is_zero(a) || self.pow(a, (self.order() - 1) / 2) == self.one()
    }
    /// A square root, by Tonelli–Shanks over the field.
    fn sqrt(&self, a: Self::El) -> Option<Self::El> {
        if self.is_zero(a) {
            return Some(self.zero());
        }
        let q = self.order();
        if !self.is_square(a) {
            return None;
        }
        let mut m = q - 1;
        let mut e = 0u32;
        while m.is_multiple_of(2) {
            m /= 2;
            e += 1;
        }
        let z = self.nonresidue();
        let mut c = self.pow(z, m);
        let mut t = self.pow(a, m);
        let mut r = self.pow(a, m.div_ceil(2));
        let mut mm = e;
        loop {
            if t == self.one() {
                return Some(r);
            }
            let mut i = 0u32;
            let mut tt = t;
            while tt != self.one() {
                tt = self.sqr(tt);
                i += 1;
                if i == mm {
                    return None;
                }
            }
            let mut b = c;
            for _ in 0..(mm - i - 1) {
                b = self.sqr(b);
            }
            mm = i;
            c = self.sqr(b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
    }
    /// `Σ cᵢ pⁱ`: an injective packing into a word.
    fn pack(&self, a: Self::El) -> u64 {
        let p = self.p();
        self.coords(a)
            .iter()
            .rev()
            .fold(0u64, |acc, &c| acc * p + c)
    }
}

// ── F_{p²} ─────────────────────────────────────────────────────────

/// `a + b·s` with `s² = ν`, as `[a, b]`.
pub type Fp2El = [u64; 2];

/// `F_{p²} = F_p[s] / (s² − ν)`, `ν` the smallest non-residue.
#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
pub struct Fp2 {
    pub p: u64,
    pub nu: u64,
}

impl Fp2 {
    /// The field for an odd prime `p < 2^31`.
    pub fn new(p: u64) -> Result<Self, String> {
        if !(3..(1 << 31)).contains(&p) || !is_prime_u64(p) {
            return Err(format!("p = {p} must be an odd prime below 2^31"));
        }
        let nu = (2..p)
            .find(|&n| jacobi(n, p) == -1)
            .ok_or("no quadratic non-residue")?;
        Ok(Self { p, nu })
    }
    /// `a₀² − ν a₁²`, the norm to `F_p`.
    pub fn norm(&self, a: Fp2El) -> u64 {
        subm(
            mulm(a[0], a[0], self.p),
            mulm(self.nu, mulm(a[1], a[1], self.p), self.p),
            self.p,
        )
    }
}

impl ExtField for Fp2 {
    type El = Fp2El;
    fn p(&self) -> u64 {
        self.p
    }
    fn degree(&self) -> u32 {
        2
    }
    fn zero(&self) -> Fp2El {
        [0, 0]
    }
    fn one(&self) -> Fp2El {
        [1, 0]
    }
    fn s(&self) -> Fp2El {
        [0, 1]
    }
    fn from_fp(&self, a: u64) -> Fp2El {
        [a % self.p, 0]
    }
    fn coords(&self, a: Fp2El) -> Vec<u64> {
        a.to_vec()
    }
    fn from_coords(&self, c: &[u64]) -> Fp2El {
        [c[0] % self.p, c.get(1).copied().unwrap_or(0) % self.p]
    }
    fn is_zero(&self, a: Fp2El) -> bool {
        a == [0, 0]
    }
    fn add(&self, a: Fp2El, b: Fp2El) -> Fp2El {
        [addm(a[0], b[0], self.p), addm(a[1], b[1], self.p)]
    }
    fn sub(&self, a: Fp2El, b: Fp2El) -> Fp2El {
        [subm(a[0], b[0], self.p), subm(a[1], b[1], self.p)]
    }
    fn neg(&self, a: Fp2El) -> Fp2El {
        [subm(0, a[0], self.p), subm(0, a[1], self.p)]
    }
    fn mul(&self, a: Fp2El, b: Fp2El) -> Fp2El {
        let p = self.p;
        let a0b0 = mulm(a[0], b[0], p);
        let a1b1 = mulm(a[1], b[1], p);
        let cross = addm(mulm(a[0], b[1], p), mulm(a[1], b[0], p), p);
        [addm(a0b0, mulm(self.nu, a1b1, p), p), cross]
    }
    fn scale(&self, a: Fp2El, k: u64) -> Fp2El {
        [mulm(a[0], k, self.p), mulm(a[1], k, self.p)]
    }
    fn inv(&self, a: Fp2El) -> Option<Fp2El> {
        let n = self.norm(a);
        if n == 0 {
            return None;
        }
        let ni = inv_mod(n, self.p);
        Some(self.scale(self.frob(a), ni))
    }
    /// `a^p = a₀ − a₁ s`: the conjugate.
    fn frob(&self, a: Fp2El) -> Fp2El {
        [a[0], subm(0, a[1], self.p)]
    }
    fn nonresidue(&self) -> Fp2El {
        // s is a non-square in F_{p²}: N(s) = −ν is a non-square in F_p
        // exactly when −1 is a square, so pick by the criterion.
        let s = self.s();
        if !self.is_square(s) {
            return s;
        }
        let one_plus_s = [1, 1];
        if !self.is_square(one_plus_s) {
            return one_plus_s;
        }
        (2..self.p)
            .map(|k| [k, 1])
            .find(|&e| !self.is_square(e))
            .expect("a non-residue exists")
    }
    /// The closed form: `(x₀ + x₁ s)² = a` with `x₀² = (a₀ ± √N)/2`.
    fn sqrt(&self, a: Fp2El) -> Option<Fp2El> {
        let p = self.p;
        if self.is_zero(a) {
            return Some(self.zero());
        }
        if a[1] == 0 {
            return match sqrt_mod(a[0], p) {
                Some(r) => Some([r, 0]),
                None => sqrt_mod(mulm(a[0], inv_mod(self.nu, p), p), p).map(|r| [0, r]),
            };
        }
        let n = sqrt_mod(self.norm(a), p)?;
        let half = inv_mod(2, p);
        for sign in [n, subm(0, n, p)] {
            let x0_sq = mulm(addm(a[0], sign, p), half, p);
            if let Some(x0) = sqrt_mod(x0_sq, p) {
                if x0 == 0 {
                    continue;
                }
                let x1 = mulm(a[1], inv_mod(mulm(2, x0, p), p), p);
                let root = [x0, x1];
                if self.sqr(root) == a {
                    return Some(root);
                }
            }
        }
        None
    }
}

// ── F_{p³} ─────────────────────────────────────────────────────────

/// `a + b·s + c·s²` with `s³ = ν`, as `[a, b, c]`.
pub type Fp3El = [u64; 3];

/// `F_{p³} = F_p[s] / (s³ − ν)`, `ν` a non-cube; needs `p ≡ 1 (mod 3)`
/// so that a non-cube exists.  The Frobenius is `s ↦ ω s` with
/// `ω = ν^{(p−1)/3}` a primitive cube root of unity.
#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq)]
pub struct Fp3 {
    pub p: u64,
    pub nu: u64,
    pub omega: u64,
    nonresidue: u64,
}

impl Fp3 {
    /// The field for a prime `p ≡ 1 (mod 3)`, `p < 2^20`.
    pub fn new(p: u64) -> Result<Self, String> {
        if !(7..(1 << 20)).contains(&p) || !is_prime_u64(p) {
            return Err(format!("p = {p} must be a prime below 2^20"));
        }
        if p % 3 != 1 {
            return Err(format!("p = {p} ≢ 1 (mod 3): every element is a cube"));
        }
        let nu = (2..p)
            .find(|&n| pow_mod(n, (p - 1) / 3, p) != 1)
            .ok_or("no non-cube")?;
        let omega = pow_mod(nu, (p - 1) / 3, p);
        let nonresidue = (2..p)
            .find(|&n| jacobi(n, p) == -1)
            .ok_or("no quadratic non-residue")?;
        Ok(Self {
            p,
            nu,
            omega,
            nonresidue,
        })
    }
    /// `a · a^p · a^{p²} ∈ F_p`.
    pub fn norm(&self, a: Fp3El) -> u64 {
        let c = self.frob(a);
        let cc = self.frob(c);
        let n = self.mul(self.mul(a, c), cc);
        debug_assert!(n[1] == 0 && n[2] == 0, "the norm is in F_p");
        n[0]
    }
}

impl ExtField for Fp3 {
    type El = Fp3El;
    fn p(&self) -> u64 {
        self.p
    }
    fn degree(&self) -> u32 {
        3
    }
    fn zero(&self) -> Fp3El {
        [0, 0, 0]
    }
    fn one(&self) -> Fp3El {
        [1, 0, 0]
    }
    fn s(&self) -> Fp3El {
        [0, 1, 0]
    }
    fn from_fp(&self, a: u64) -> Fp3El {
        [a % self.p, 0, 0]
    }
    fn coords(&self, a: Fp3El) -> Vec<u64> {
        a.to_vec()
    }
    fn from_coords(&self, c: &[u64]) -> Fp3El {
        let g = |i: usize| c.get(i).copied().unwrap_or(0) % self.p;
        [g(0), g(1), g(2)]
    }
    fn is_zero(&self, a: Fp3El) -> bool {
        a == [0, 0, 0]
    }
    fn add(&self, a: Fp3El, b: Fp3El) -> Fp3El {
        let p = self.p;
        [
            addm(a[0], b[0], p),
            addm(a[1], b[1], p),
            addm(a[2], b[2], p),
        ]
    }
    fn sub(&self, a: Fp3El, b: Fp3El) -> Fp3El {
        let p = self.p;
        [
            subm(a[0], b[0], p),
            subm(a[1], b[1], p),
            subm(a[2], b[2], p),
        ]
    }
    fn neg(&self, a: Fp3El) -> Fp3El {
        let p = self.p;
        [subm(0, a[0], p), subm(0, a[1], p), subm(0, a[2], p)]
    }
    fn mul(&self, a: Fp3El, b: Fp3El) -> Fp3El {
        let p = self.p;
        let mut c = [0u64; 5];
        for i in 0..3 {
            for j in 0..3 {
                c[i + j] = addm(c[i + j], mulm(a[i], b[j], p), p);
            }
        }
        [
            addm(c[0], mulm(self.nu, c[3], p), p),
            addm(c[1], mulm(self.nu, c[4], p), p),
            c[2],
        ]
    }
    fn scale(&self, a: Fp3El, k: u64) -> Fp3El {
        let p = self.p;
        [mulm(a[0], k, p), mulm(a[1], k, p), mulm(a[2], k, p)]
    }
    fn inv(&self, a: Fp3El) -> Option<Fp3El> {
        let n = self.norm(a);
        if n == 0 {
            return None;
        }
        let c = self.frob(a);
        let cc = self.frob(c);
        Some(self.scale(self.mul(c, cc), inv_mod(n, self.p)))
    }
    /// `a^p = a₀ + a₁ ω s + a₂ ω² s²`.
    fn frob(&self, a: Fp3El) -> Fp3El {
        let p = self.p;
        [
            a[0],
            mulm(a[1], self.omega, p),
            mulm(a[2], mulm(self.omega, self.omega, p), p),
        ]
    }
    fn nonresidue(&self) -> Fp3El {
        [self.nonresidue, 0, 0]
    }
}

// ── The curve ──────────────────────────────────────────────────────

/// A point of `E(F_{p^k})`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct ExtPoint<El> {
    pub x: El,
    pub y: El,
    pub infinity: bool,
}

impl<El: Copy> ExtPoint<El> {
    pub fn affine(x: El, y: El) -> Self {
        Self {
            x,
            y,
            infinity: false,
        }
    }
    pub fn infinity(zero: El) -> Self {
        Self {
            x: zero,
            y: zero,
            infinity: true,
        }
    }
}

/// `y² = x³ + Ax + B` over `F_{p^k}`, every operation charged.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct ExtCurve<F: ExtField> {
    pub f: F,
    pub a: F::El,
    pub b: F::El,
}

impl<F: ExtField> ExtCurve<F> {
    /// `x³ + Ax + B`.
    pub fn rhs(&self, x: F::El) -> F::El {
        let f = &self.f;
        f.add(f.add(f.mul(f.sqr(x), x), f.mul(self.a, x)), self.b)
    }
    pub fn is_on_curve(&self, pt: ExtPoint<F::El>) -> bool {
        pt.infinity || self.f.sqr(pt.y) == self.rhs(pt.x)
    }
    /// The two points above `x`, if any (one if `y = 0`).
    pub fn lift_x(&self, x: F::El) -> Vec<ExtPoint<F::El>> {
        match self.f.sqrt(self.rhs(x)) {
            None => Vec::new(),
            Some(y) if self.f.is_zero(y) => vec![ExtPoint::affine(x, y)],
            Some(y) => vec![ExtPoint::affine(x, y), ExtPoint::affine(x, self.f.neg(y))],
        }
    }
    fn add_raw(&self, a: ExtPoint<F::El>, b: ExtPoint<F::El>) -> ExtPoint<F::El> {
        if a.infinity {
            return b;
        }
        if b.infinity {
            return a;
        }
        let f = &self.f;
        if a.x == b.x {
            if f.is_zero(f.add(a.y, b.y)) {
                return ExtPoint::infinity(f.zero());
            }
            return self.double_raw(a);
        }
        let lambda = f.mul(f.sub(b.y, a.y), f.inv(f.sub(b.x, a.x)).expect("distinct x"));
        let x3 = f.sub(f.sub(f.sqr(lambda), a.x), b.x);
        let y3 = f.sub(f.mul(lambda, f.sub(a.x, x3)), a.y);
        ExtPoint::affine(x3, y3)
    }
    fn double_raw(&self, a: ExtPoint<F::El>) -> ExtPoint<F::El> {
        if a.infinity || self.f.is_zero(a.y) {
            return ExtPoint::infinity(self.f.zero());
        }
        let f = &self.f;
        let num = f.add(f.scale(f.sqr(a.x), 3), self.a);
        let lambda = f.mul(num, f.inv(f.scale(a.y, 2)).expect("y ≠ 0"));
        let x3 = f.sub(f.sqr(lambda), f.scale(a.x, 2));
        let y3 = f.sub(f.mul(lambda, f.sub(a.x, x3)), a.y);
        ExtPoint::affine(x3, y3)
    }
}

impl<F: ExtField> CountedGroup for ExtCurve<F> {
    type Elt = ExtPoint<F::El>;
    fn identity(&self) -> Self::Elt {
        ExtPoint::infinity(self.f.zero())
    }
    fn is_identity(&self, p: &Self::Elt) -> bool {
        p.infinity
    }
    fn add(&self, ops: &mut GroupOps, p: Self::Elt, q: Self::Elt) -> Self::Elt {
        ops.adds += 1;
        self.add_raw(p, q)
    }
    fn double(&self, ops: &mut GroupOps, p: Self::Elt) -> Self::Elt {
        ops.doubles += 1;
        self.double_raw(p)
    }
    fn neg(&self, p: Self::Elt) -> Self::Elt {
        if p.infinity {
            p
        } else {
            ExtPoint::affine(p.x, self.f.neg(p.y))
        }
    }
    fn key(&self, p: &Self::Elt) -> u64 {
        if p.infinity {
            return 0;
        }
        let sign = u64::from(self.f.pack(p.y) > self.f.pack(self.f.neg(p.y)));
        ((self.f.pack(p.x) + 1) << 1) | sign
    }
}

/// `(x, y) ↦ (cx·x, cy·y)` on a curve over an extension field: the
/// lift of a `j = 0` or `j = 1728` automorphism to a twist, or any
/// other diagonal automorphism.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct ExtDiagonalAutomorphism<El> {
    pub cx: El,
    pub cy: El,
    pub eigenvalue: u64,
    pub order: u32,
    pub label: &'static str,
}

impl<F: ExtField> crate::cryptanalysis::glv_invariant_base::Endomorphism<ExtCurve<F>>
    for ExtDiagonalAutomorphism<F::El>
{
    fn name(&self) -> String {
        self.label.to_string()
    }
    fn degree(&self) -> u64 {
        1
    }
    fn eigenvalue(&self) -> u64 {
        self.eigenvalue
    }
    fn apply(&self, g: &ExtCurve<F>, p: ExtPoint<F::El>) -> ExtPoint<F::El> {
        if p.infinity {
            return p;
        }
        let f = &g.f;
        ExtPoint::affine(scale_by(f, self.cx, p.x), scale_by(f, self.cy, p.y))
    }
}

/// `c·a`, free when `c = 1`: the plain Frobenius and the `y`-coordinate
/// of `ζ` carry a unit constant, and a walk canonicalises on every
/// step, so the multiplication by one is not paid.
#[inline]
fn scale_by<F: ExtField>(f: &F, c: F::El, a: F::El) -> F::El {
    if c == f.one() {
        a
    } else {
        f.mul(c, a)
    }
}

/// `(x, y) ↦ (cx·x^p, cy·y^p)`: a Frobenius-type map — the GLS `ψ`
/// with its twist constants, or the plain Frobenius with `cx = cy = 1`.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct FrobeniusType<El> {
    pub cx: El,
    pub cy: El,
    pub eigenvalue: u64,
    pub p: u64,
    pub label: &'static str,
}

impl<F: ExtField> crate::cryptanalysis::glv_invariant_base::Endomorphism<ExtCurve<F>>
    for FrobeniusType<F::El>
{
    fn name(&self) -> String {
        self.label.to_string()
    }
    fn degree(&self) -> u64 {
        self.p
    }
    fn eigenvalue(&self) -> u64 {
        self.eigenvalue
    }
    fn apply(&self, g: &ExtCurve<F>, pt: ExtPoint<F::El>) -> ExtPoint<F::El> {
        if pt.infinity {
            return pt;
        }
        let f = &g.f;
        ExtPoint::affine(
            scale_by(f, self.cx, f.frob(pt.x)),
            scale_by(f, self.cy, f.frob(pt.y)),
        )
    }
}

/// A random affine point with `y ≠ 0`.
pub fn random_point<F: ExtField>(
    curve: &ExtCurve<F>,
    rng: &mut rand::rngs::StdRng,
) -> ExtPoint<F::El> {
    use rand::Rng;
    loop {
        let coords: Vec<u64> = (0..curve.f.degree())
            .map(|_| rng.gen_range(0..curve.f.p()))
            .collect();
        let x = curve.f.from_coords(&coords);
        let pts = curve.lift_x(x);
        if let Some(&pt) = pts.first() {
            if !curve.f.is_zero(pt.y) {
                return pt;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    #[test]
    fn fp3_is_a_field_with_a_frobenius_of_order_three() {
        let f = Fp3::new(1009).unwrap(); // 1009 ≡ 1 (mod 3)
        let mut rng = StdRng::seed_from_u64(5);
        for _ in 0..200 {
            let a = f.from_coords(&[
                rng.gen_range(0..1009),
                rng.gen_range(0..1009),
                rng.gen_range(0..1009),
            ]);
            let b = f.from_coords(&[
                rng.gen_range(0..1009),
                rng.gen_range(0..1009),
                rng.gen_range(0..1009),
            ]);
            assert_eq!(f.mul(a, b), f.mul(b, a));
            assert_eq!(f.mul(f.add(a, b), b), f.add(f.mul(a, b), f.mul(b, b)));
            if !f.is_zero(a) {
                assert_eq!(f.mul(a, f.inv(a).unwrap()), f.one());
            }
            assert_eq!(f.pow(a, 1009), f.frob(a), "a^p is the Frobenius");
            assert_eq!(f.frob(f.frob(f.frob(a))), a, "order 3");
            assert_eq!(f.frob(f.mul(a, b)), f.mul(f.frob(a), f.frob(b)));
            let sq = f.sqr(a);
            let root = f.sqrt(sq).expect("a square has a root");
            assert!(root == a || root == f.neg(a));
            assert_eq!(f.is_square(a), f.sqrt(a).is_some());
        }
        assert_eq!(f.order(), 1009 * 1009 * 1009);
        assert!(!f.is_square(f.nonresidue()));
    }

    #[test]
    fn fp2_generic_and_closed_form_roots_agree() {
        let f = Fp2::new(1013).unwrap();
        let mut rng = StdRng::seed_from_u64(6);
        for _ in 0..200 {
            let a = [rng.gen_range(0..1013), rng.gen_range(0..1013)];
            let closed = f.sqrt(a);
            // The trait's generic Tonelli–Shanks, forced by calling it on a
            // wrapper that does not override sqrt.
            #[derive(Clone, Copy, Debug, Serialize)]
            struct Generic(Fp2);
            impl ExtField for Generic {
                type El = Fp2El;
                fn p(&self) -> u64 {
                    self.0.p
                }
                fn degree(&self) -> u32 {
                    2
                }
                fn zero(&self) -> Fp2El {
                    self.0.zero()
                }
                fn one(&self) -> Fp2El {
                    self.0.one()
                }
                fn s(&self) -> Fp2El {
                    self.0.s()
                }
                fn from_fp(&self, a: u64) -> Fp2El {
                    self.0.from_fp(a)
                }
                fn coords(&self, a: Fp2El) -> Vec<u64> {
                    self.0.coords(a)
                }
                fn from_coords(&self, c: &[u64]) -> Fp2El {
                    self.0.from_coords(c)
                }
                fn is_zero(&self, a: Fp2El) -> bool {
                    self.0.is_zero(a)
                }
                fn add(&self, a: Fp2El, b: Fp2El) -> Fp2El {
                    self.0.add(a, b)
                }
                fn sub(&self, a: Fp2El, b: Fp2El) -> Fp2El {
                    self.0.sub(a, b)
                }
                fn neg(&self, a: Fp2El) -> Fp2El {
                    self.0.neg(a)
                }
                fn mul(&self, a: Fp2El, b: Fp2El) -> Fp2El {
                    self.0.mul(a, b)
                }
                fn scale(&self, a: Fp2El, k: u64) -> Fp2El {
                    self.0.scale(a, k)
                }
                fn inv(&self, a: Fp2El) -> Option<Fp2El> {
                    self.0.inv(a)
                }
                fn frob(&self, a: Fp2El) -> Fp2El {
                    self.0.frob(a)
                }
                fn nonresidue(&self) -> Fp2El {
                    self.0.nonresidue()
                }
            }
            let generic = Generic(f).sqrt(a);
            assert_eq!(closed.is_some(), generic.is_some());
            if let (Some(c), Some(g)) = (closed, generic) {
                assert!(c == g || c == f.neg(g));
            }
        }
    }

    #[test]
    fn a_curve_over_fp3_is_a_group() {
        let f = Fp3::new(1009).unwrap();
        let curve = ExtCurve {
            f,
            a: f.from_fp(3),
            b: f.from_fp(7),
        };
        let mut rng = StdRng::seed_from_u64(7);
        let mut ops = GroupOps::default();
        for _ in 0..20 {
            let p = random_point(&curve, &mut rng);
            let q = random_point(&curve, &mut rng);
            let s = curve.add(&mut ops, p, q);
            assert!(curve.is_on_curve(s));
            assert_eq!(curve.add(&mut ops, s, curve.neg(q)), p);
            assert_eq!(curve.add(&mut ops, p, p), curve.double(&mut ops, p));
            let r = random_point(&curve, &mut rng);
            let pq = curve.add(&mut ops, p, q);
            let lhs = curve.add(&mut ops, pq, r);
            let qr = curve.add(&mut ops, q, r);
            let rhs = curve.add(&mut ops, p, qr);
            assert_eq!(lhs, rhs, "associative");
        }
    }
}
