//! # Genus-2 Jacobians over `F_p` with `u64` arithmetic.
//!
//! A self-contained toolkit for the Jacobian of an imaginary genus-2 curve
//!
//! ```text
//!     C : y² = f(x),   deg f = 5,   f squarefree,   p odd prime < 2^31
//! ```
//!
//! over a small prime field. It was moved out of the round-4 research example
//! `examples/r4_16_6_genus2_census.rs` (lead 16.6 of
//! `RESEARCH_ISOGENY_FIBERS_ECDLP.md`, §17.6 and §17.8), where it ran the
//! affine-saturation census on `y² = x⁵ + b`.
//!
//! What it provides:
//!
//! - [`Genus2Curve`]: the curve, validated at construction (prime `p`, exact
//!   degree 5, squarefree `f`).
//! - [`Divisor`]: a reduced Mumford pair `(u, v)` with `u` monic, `deg u ≤ 2`,
//!   `deg v < deg u` and `u | f − v²`, stored as fixed `u64` arrays so it is
//!   `Copy`, `Hash` and `Ord`.
//! - Cantor addition with reduction, doubling, negation, subtraction and
//!   double-and-add scalar multiplication.
//! - The automorphism `ζ : (x, y) ↦ (ζ₅x, y)` of `y² = c·x⁵ + b` when
//!   `p ≡ 1 (mod 5)`, and with the hyperelliptic involution the cyclic
//!   automorphism group of order 10 ([`Genus2Curve::automorphism`]).
//! - `#C(F_p)`, `#C(F_{p²})`, the `L`-polynomial and `#J(F_p) = L(1)`
//!   ([`Genus2Curve::l_polynomial`]).
//! - Exhaustive enumeration of `J(F_p)` for small `p`: every monic `u` of
//!   degree at most 2 is visited and all its `v` are solved for and checked.
//! - Uniform random sampling of `J(F_p)` ([`Genus2Curve::random`]).
//!
//! ## Relation to the other genus-2 code
//!
//! [`crate::prime_hyperelliptic`] has the genus-generic, `BigUint`-backed
//! Mumford/Cantor implementation, and [`crate::cryptanalysis::jv_cover::Hyp`]
//! is generic over a field trait (used there over `F_{p²}` and `F_{p³}`).
//! This module is the specialised fast path: genus 2, one prime below `2^31`,
//! machine-word coefficients, plus the `ζ`-action and the exhaustive
//! enumerator that neither of those has. The unit tests cross-check its
//! Cantor addition against both and its `#J` against
//! [`crate::prime_hyperelliptic::brute_force_jac_order_via_lpoly`].
//!
//! ## Bounds
//!
//! Every coefficient is kept below `p < 2^31` ([`MAX_PRIME`]), so each
//! product of two coefficients, and each sum of two such products, fits in a
//! `u64`. The routines whose cost grows with `p` refuse larger primes with
//! [`Genus2Error::PrimeTooLargeFor`] instead of running for hours:
//! point counting and `L(1)` need `p ≤` [`MAX_POINT_COUNT_PRIME`]
//! (`O(p²)` work), exhaustive counting needs `p ≤` [`MAX_ENUMERATION_PRIME`]
//! (`O(p²)` square roots), and listing every element needs
//! `p ≤` [`MAX_LIST_PRIME`] (memory `O(p²)`).
//!
//! ## Example
//!
//! ```
//! use crypto_lib::cryptanalysis::genus2_fp::Genus2Curve;
//!
//! // y² = x⁵ + 3 over F_31 (31 ≡ 1 mod 5, so ζ₅ ∈ F_31).
//! let c = Genus2Curve::x5_plus_b(31, 3).unwrap();
//! let order = c.jacobian_order().unwrap();
//! assert_eq!(order, c.count_exhaustive().unwrap());
//!
//! let mut rng = rand::thread_rng();
//! let d = c.random(&mut rng);
//! assert!(c.mul(&d, order).is_identity());
//! let z = c.zeta(&d).unwrap();
//! assert_eq!(c.zeta_pow(&z, 4).unwrap(), d); // ζ⁵ = id
//! ```

use std::fmt;

use rand::Rng;

use super::residual_walk::is_prime_u64;

/// Exclusive upper bound on the prime: `p < 2^31`. Below it every product of
/// two reduced coefficients plus another such product fits in a `u64`.
pub const MAX_PRIME: u64 = 1 << 31;

/// Largest prime accepted by the point counters and by `L(1)`; counting over
/// `F_{p²}` visits all `p²` abscissae.
pub const MAX_POINT_COUNT_PRIME: u64 = 1 << 16;

/// Largest prime accepted by [`Genus2Curve::for_each_element`] and
/// [`Genus2Curve::count_exhaustive`]; they visit all `p² + p + 1` monic `u`.
pub const MAX_ENUMERATION_PRIME: u64 = 1 << 13;

/// Largest prime accepted by [`Genus2Curve::enumerate`], which stores all
/// `#J(F_p) ≈ p²` elements.
pub const MAX_LIST_PRIME: u64 = 1 << 10;

/// Why a curve, divisor or request was refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Genus2Error {
    /// `p` is not an odd prime below [`MAX_PRIME`].
    PrimeOutOfRange {
        /// The rejected modulus.
        p: u64,
    },
    /// A coefficient of `f` is not reduced modulo `p`.
    CoefficientNotReduced {
        /// The offending coefficient.
        value: u64,
    },
    /// The coefficient of `x⁵` is zero, so the model is not imaginary of
    /// genus 2.
    NotQuintic,
    /// `f` has a repeated root over the algebraic closure, so `C` is
    /// singular.
    NotSquarefree,
    /// A polynomial pair is not a reduced Mumford representation on this
    /// curve.
    InvalidDivisor,
    /// The `ζ₅`-action needs `f = c·x⁵ + b` and `p ≡ 1 (mod 5)`.
    NoZetaAction,
    /// The requested routine is too costly at this `p`.
    PrimeTooLargeFor {
        /// The curve's prime.
        p: u64,
        /// The largest prime the routine accepts.
        limit: u64,
        /// Which routine refused.
        what: &'static str,
    },
}

impl fmt::Display for Genus2Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Genus2Error::PrimeOutOfRange { p } => {
                write!(f, "p = {p} is not an odd prime below 2^31")
            }
            Genus2Error::CoefficientNotReduced { value } => {
                write!(f, "coefficient {value} is not reduced modulo p")
            }
            Genus2Error::NotQuintic => write!(f, "f must have exact degree 5"),
            Genus2Error::NotSquarefree => write!(f, "f is not squarefree"),
            Genus2Error::InvalidDivisor => write!(f, "not a reduced Mumford pair on this curve"),
            Genus2Error::NoZetaAction => {
                write!(f, "the ζ₅ action needs f = c·x⁵ + b and p ≡ 1 (mod 5)")
            }
            Genus2Error::PrimeTooLargeFor { p, limit, what } => {
                write!(f, "{what} supports p ≤ {limit}, got p = {p}")
            }
        }
    }
}

impl std::error::Error for Genus2Error {}

// ── F_p and F_{p²} ──────────────────────────────────────────────────────────

/// `F_p` with a fixed quadratic non-residue `w`, `F_{p²} = F_p[t]/(t² − w)`.
#[derive(Clone, Debug)]
struct Fp {
    p: u64,
    w: u64,
}

type Fp2 = (u64, u64);

impl Fp {
    fn new(p: u64) -> Fp {
        let mut k = Fp { p, w: 0 };
        // Half of F_p^* are non-residues, so the search ends quickly.
        k.w = (2..p).find(|&w| k.legendre(w) == -1).unwrap_or(p - 1);
        k
    }
    #[inline]
    fn add(&self, a: u64, b: u64) -> u64 {
        (a + b) % self.p
    }
    #[inline]
    fn sub(&self, a: u64, b: u64) -> u64 {
        (a + self.p - b) % self.p
    }
    #[inline]
    fn neg(&self, a: u64) -> u64 {
        (self.p - a) % self.p
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        a * b % self.p
    }
    fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1 % self.p;
        a %= self.p;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
    /// Inverse of a nonzero element (Fermat); `0` maps to `0`.
    fn inv(&self, a: u64) -> u64 {
        self.pow(a, self.p - 2)
    }
    /// Legendre symbol by Euler's criterion: `0`, `1` or `-1`.
    fn legendre(&self, a: u64) -> i8 {
        let a = a % self.p;
        if a == 0 {
            return 0;
        }
        if self.pow(a, (self.p - 1) / 2) == 1 {
            1
        } else {
            -1
        }
    }
    /// Tonelli–Shanks; `None` for a non-residue.
    fn sqrt(&self, a: u64) -> Option<u64> {
        let p = self.p;
        let a = a % p;
        if a == 0 {
            return Some(0);
        }
        if self.legendre(a) != 1 {
            return None;
        }
        let mut q = p - 1;
        let mut s = 0u32;
        while q.is_multiple_of(2) {
            q /= 2;
            s += 1;
        }
        let mut m = s;
        let mut c = self.pow(self.w, q);
        let mut t = self.pow(a, q);
        let mut r = self.pow(a, q.div_ceil(2));
        while t != 1 {
            let mut i = 0u32;
            let mut tt = t;
            while tt != 1 && i < m {
                tt = self.mul(tt, tt);
                i += 1;
            }
            if i >= m {
                return None; // unreachable for a residue
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = self.mul(b, b);
            }
            m = i;
            c = self.mul(b, b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
        Some(r)
    }
    fn m2(&self, a: Fp2, b: Fp2) -> Fp2 {
        let p = self.p;
        (
            (a.0 * b.0 % p + self.w * (a.1 * b.1 % p)) % p,
            (a.0 * b.1 % p + a.1 * b.0 % p) % p,
        )
    }
    fn norm2(&self, a: Fp2) -> u64 {
        self.sub(self.mul(a.0, a.0), self.mul(self.w, self.mul(a.1, a.1)))
    }
    /// Square root in `F_{p²}` by the norm method, verified before return.
    fn sqrt2(&self, a: Fp2) -> Option<Fp2> {
        if a.1 == 0 {
            if let Some(s) = self.sqrt(a.0) {
                return Some((s, 0));
            }
            let s = self.sqrt(self.mul(a.0, self.inv(self.w)))?;
            return Some((0, s));
        }
        let s = self.sqrt(self.norm2(a))?;
        let half = self.inv(2);
        // The two candidates for x0² multiply to w·a1²/4, a non-residue, so
        // exactly one is a residue.
        let mut x0sq = self.mul(self.add(a.0, s), half);
        if self.legendre(x0sq) != 1 {
            x0sq = self.mul(self.sub(a.0, s), half);
        }
        let x0 = self.sqrt(x0sq)?;
        if x0 == 0 {
            return None;
        }
        let x1 = self.mul(a.1, self.inv(self.mul(2, x0)));
        let r = (x0, x1);
        (self.m2(r, r) == a).then_some(r)
    }
}

// ── Polynomials over F_p (trimmed Vec, lowest degree first) ────────────────

type Poly = Vec<u64>;

fn trim(mut a: Poly) -> Poly {
    while a.last() == Some(&0) {
        a.pop();
    }
    a
}

fn padd(k: &Fp, a: &[u64], b: &[u64]) -> Poly {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|i| k.add(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0)))
            .collect(),
    )
}

fn psub(k: &Fp, a: &[u64], b: &[u64]) -> Poly {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|i| k.sub(*a.get(i).unwrap_or(&0), *b.get(i).unwrap_or(&0)))
            .collect(),
    )
}

fn pmul(k: &Fp, a: &[u64], b: &[u64]) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = (r[i + j] + x * y) % k.p;
        }
    }
    trim(r)
}

fn pscale(k: &Fp, a: &[u64], c: u64) -> Poly {
    trim(a.iter().map(|&x| k.mul(x, c)).collect())
}

/// Quotient and remainder; division by the zero polynomial returns `(0, a)`.
fn pdivrem(k: &Fp, a: &[u64], b: &[u64]) -> (Poly, Poly) {
    let mut r: Poly = trim(a.to_vec());
    let Some(&lb) = b.last() else {
        return (Vec::new(), r);
    };
    if lb == 0 || r.len() < b.len() {
        return (Vec::new(), r);
    }
    let db = b.len() - 1;
    let li = k.inv(lb);
    let mut q = vec![0u64; r.len() - db];
    while r.len() >= b.len() {
        let d = r.len() - 1;
        let c = k.mul(r[d], li);
        q[d - db] = c;
        for (j, &bj) in b.iter().enumerate() {
            r[d - db + j] = k.sub(r[d - db + j], k.mul(c, bj));
        }
        r = trim(r);
    }
    (trim(q), r)
}

fn pmonic(k: &Fp, a: &[u64]) -> Poly {
    match a.last() {
        None => Vec::new(),
        Some(&lc) => pscale(k, a, k.inv(lc)),
    }
}

/// `(d, s, t)` with `d = s·a + t·b` and `d` monic (or zero).
fn pxgcd(k: &Fp, a: &[u64], b: &[u64]) -> (Poly, Poly, Poly) {
    let (mut r0, mut r1) = (trim(a.to_vec()), trim(b.to_vec()));
    let (mut s0, mut s1): (Poly, Poly) = (vec![1], Vec::new());
    let (mut t0, mut t1): (Poly, Poly) = (Vec::new(), vec![1]);
    while !r1.is_empty() {
        let (q, r) = pdivrem(k, &r0, &r1);
        let s = psub(k, &s0, &pmul(k, &q, &s1));
        let t = psub(k, &t0, &pmul(k, &q, &t1));
        r0 = std::mem::replace(&mut r1, r);
        s0 = std::mem::replace(&mut s1, s);
        t0 = std::mem::replace(&mut t1, t);
    }
    match r0.last() {
        None => (r0, s0, t0),
        Some(&lc) => {
            let i = k.inv(lc);
            (pscale(k, &r0, i), pscale(k, &s0, i), pscale(k, &t0, i))
        }
    }
}

fn peval(k: &Fp, a: &[u64], x: u64) -> u64 {
    a.iter().rev().fold(0, |acc, &c| k.add(k.mul(acc, x), c))
}

fn peval2(k: &Fp, a: &[u64], x: Fp2) -> Fp2 {
    a.iter().rev().fold((0, 0), |acc, &c| {
        let m = k.m2(acc, x);
        (k.add(m.0, c), m.1)
    })
}

fn pderiv(k: &Fp, a: &[u64]) -> Poly {
    trim(
        a.iter()
            .enumerate()
            .skip(1)
            .map(|(i, &c)| k.mul(i as u64 % k.p, c))
            .collect(),
    )
}

// ── Divisors ────────────────────────────────────────────────────────────────

/// A reduced Mumford pair `(u, v)` on a genus-2 curve: `u` monic of degree at
/// most 2, `deg v < deg u`, `u | f − v²`. The identity is `(1, 0)`.
///
/// Stored with fixed arrays (`u = x^deg + u[1]·x + u[0]` truncated to its
/// degree, `v = v[1]·x + v[0]`; unused slots are zero), so equality is
/// equality of divisor classes. Values are obtained from a [`Genus2Curve`],
/// which validates them; a `Divisor` does not record its curve.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct Divisor {
    deg: u8,
    u: [u64; 2],
    v: [u64; 2],
}

impl Divisor {
    /// The identity class `(1, 0)`.
    pub const IDENTITY: Divisor = Divisor {
        deg: 0,
        u: [0, 0],
        v: [0, 0],
    };

    /// `deg u` (0, 1 or 2).
    pub fn degree(&self) -> usize {
        self.deg as usize
    }

    /// Whether this is the identity class.
    pub fn is_identity(&self) -> bool {
        self.deg == 0
    }

    /// The coefficients of `u`, lowest degree first, including the leading 1
    /// (length `deg u + 1`).
    pub fn u_coeffs(&self) -> Vec<u64> {
        let d = self.degree();
        let mut out = self.u[..d].to_vec();
        out.push(1);
        out
    }

    /// The coefficients of `v`, lowest degree first, without trailing zeros
    /// (empty for `v = 0`).
    pub fn v_coeffs(&self) -> Vec<u64> {
        trim(self.v[..self.degree()].to_vec())
    }

    /// Builds the fixed representation from trimmed polynomials with
    /// `deg u ≤ 2` and `deg v < deg u`; other shapes give `None`.
    fn from_polys(u: &[u64], v: &[u64]) -> Option<Divisor> {
        let deg = u.len().checked_sub(1)?;
        if deg > 2 || u[deg] != 1 || v.len() > deg {
            return None;
        }
        let mut d = Divisor {
            deg: deg as u8,
            ..Divisor::IDENTITY
        };
        d.u[..deg].copy_from_slice(&u[..deg]);
        d.v[..v.len()].copy_from_slice(v);
        Some(d)
    }
}

/// The `L`-polynomial `L(T) = 1 + a1·T + a2·T² + p·a1·T³ + p²·T⁴` of a
/// genus-2 curve over `F_p`, from its point counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LPolynomial {
    /// The field size `p`.
    pub p: u64,
    /// `#C(F_p)`, counting the single point at infinity.
    pub n1: u64,
    /// `#C(F_{p²})`, counting the single point at infinity.
    pub n2: u64,
    /// `a1 = n1 − p − 1`.
    pub a1: i64,
    /// `a2 = (n2 − p² − 1 + a1²) / 2`.
    pub a2: i64,
}

impl LPolynomial {
    /// `L(1) = #J(F_p) = 1 + a1 + a2 + p·a1 + p²`.
    pub fn jacobian_order(&self) -> u64 {
        let p = self.p as i128;
        let (a1, a2) = (self.a1 as i128, self.a2 as i128);
        // Weil: (√p − 1)⁴ ≤ L(1) ≤ (√p + 1)⁴, positive and below 2^64 here.
        (1 + a1 + a2 + p * a1 + p * p).clamp(0, u64::MAX as i128) as u64
    }

    /// `L(T)` evaluated at an integer `t`.
    pub fn eval(&self, t: i128) -> i128 {
        let p = self.p as i128;
        let (a1, a2) = (self.a1 as i128, self.a2 as i128);
        1 + a1 * t + a2 * t * t + p * a1 * t * t * t + p * p * t * t * t * t
    }
}

// ── The curve ───────────────────────────────────────────────────────────────

/// A genus-2 curve `y² = f(x)` over `F_p` in the imaginary model: `deg f = 5`,
/// `f` squarefree, `p` an odd prime below [`MAX_PRIME`], one point at
/// infinity. All Jacobian operations live here and take [`Divisor`] values.
///
/// Operations on divisors that are not valid for this curve (see
/// [`Genus2Curve::is_valid`]) never panic, but their results are unspecified.
#[derive(Clone, Debug)]
pub struct Genus2Curve {
    k: Fp,
    f: [u64; 6],
    /// A primitive fifth root of unity when the `ζ`-action exists.
    zeta5: Option<u64>,
}

impl Genus2Curve {
    /// The curve `y² = f[0] + f[1]·x + … + f[5]·x⁵` over `F_p`.
    ///
    /// Requires an odd prime `p < 2^31`, coefficients reduced below `p`,
    /// `f[5] ≠ 0` and `f` squarefree.
    pub fn new(p: u64, f: [u64; 6]) -> Result<Genus2Curve, Genus2Error> {
        if !(3..MAX_PRIME).contains(&p) || !is_prime_u64(p) {
            return Err(Genus2Error::PrimeOutOfRange { p });
        }
        if let Some(&value) = f.iter().find(|&&c| c >= p) {
            return Err(Genus2Error::CoefficientNotReduced { value });
        }
        if f[5] == 0 {
            return Err(Genus2Error::NotQuintic);
        }
        let k = Fp::new(p);
        let (g, _, _) = pxgcd(&k, &f, &pderiv(&k, &f));
        if g.len() != 1 {
            return Err(Genus2Error::NotSquarefree);
        }
        let zeta5 = if f[1..5].iter().all(|&c| c == 0) && (p - 1).is_multiple_of(5) {
            (2..p).map(|g| k.pow(g, (p - 1) / 5)).find(|&z| z != 1)
        } else {
            None
        };
        Ok(Genus2Curve { k, f, zeta5 })
    }

    /// The curve `y² = x⁵ + b` (`b` is reduced modulo `p`). It is singular
    /// for `b ≡ 0` or `p = 5`, and has the `ζ`-action when `p ≡ 1 (mod 5)`.
    pub fn x5_plus_b(p: u64, b: u64) -> Result<Genus2Curve, Genus2Error> {
        if p == 0 {
            return Err(Genus2Error::PrimeOutOfRange { p });
        }
        Genus2Curve::new(p, [b % p, 0, 0, 0, 0, 1])
    }

    /// The field size `p`.
    pub fn p(&self) -> u64 {
        self.k.p
    }

    /// The coefficients of `f`, lowest degree first.
    pub fn f(&self) -> [u64; 6] {
        self.f
    }

    /// `f(x)` for `x ∈ F_p`.
    pub fn eval_f(&self, x: u64) -> u64 {
        peval(&self.k, &self.f, x % self.k.p)
    }

    /// Whether the affine point `(x, y)` (coordinates reduced mod `p`) lies
    /// on the curve.
    pub fn is_on_curve(&self, x: u64, y: u64) -> bool {
        let y = y % self.k.p;
        self.k.mul(y, y) == self.eval_f(x)
    }

    /// The identity class.
    pub fn identity(&self) -> Divisor {
        Divisor::IDENTITY
    }

    /// Whether `d` is a reduced Mumford pair on this curve: coefficients
    /// below `p` and `u | f − v²`.
    pub fn is_valid(&self, d: &Divisor) -> bool {
        let p = self.k.p;
        if d.deg > 2 || d.u.iter().chain(d.v.iter()).any(|&c| c >= p) {
            return false;
        }
        let deg = d.degree();
        if d.u[deg..].iter().chain(d.v[deg..].iter()).any(|&c| c != 0) {
            return false;
        }
        let (u, v) = (d.u_coeffs(), d.v_coeffs());
        let k = &self.k;
        pdivrem(k, &psub(k, &self.f, &pmul(k, &v, &v)), &u)
            .1
            .is_empty()
    }

    /// The divisor with Mumford pair `(u, v)`, coefficients lowest degree
    /// first. `u` must be monic of degree at most 2 (trailing zeros are
    /// ignored), `v` must reduce to degree below `deg u`, all coefficients
    /// below `p`, and `u | f − v²`.
    pub fn divisor(&self, u: &[u64], v: &[u64]) -> Result<Divisor, Genus2Error> {
        let p = self.k.p;
        if u.iter().chain(v.iter()).any(|&c| c >= p) {
            return Err(Genus2Error::InvalidDivisor);
        }
        let d = Divisor::from_polys(&trim(u.to_vec()), &trim(v.to_vec()))
            .ok_or(Genus2Error::InvalidDivisor)?;
        if self.is_valid(&d) {
            Ok(d)
        } else {
            Err(Genus2Error::InvalidDivisor)
        }
    }

    /// The class `[(x, y) − ∞]` of an affine point on the curve.
    pub fn point(&self, x: u64, y: u64) -> Result<Divisor, Genus2Error> {
        let p = self.k.p;
        if x >= p || y >= p || !self.is_on_curve(x, y) {
            return Err(Genus2Error::InvalidDivisor);
        }
        self.divisor(&[self.k.neg(x), 1], &[y])
    }

    /// Mumford polynomials of `d` with coefficients reduced mod `p`, so that
    /// arithmetic never overflows even on foreign input.
    fn polys(&self, d: &Divisor) -> (Poly, Poly) {
        let p = self.k.p;
        let deg = d.degree().min(2);
        let mut u: Poly = d.u[..deg].iter().map(|&c| c % p).collect();
        u.push(1);
        let v = trim(d.v[..deg].iter().map(|&c| c % p).collect());
        (u, v)
    }

    /// `−d = (u, −v)`, the image under the hyperelliptic involution.
    pub fn neg(&self, d: &Divisor) -> Divisor {
        let mut out = *d;
        for c in out.v.iter_mut() {
            *c = self.k.neg(*c % self.k.p);
        }
        out
    }

    /// Cantor addition: composition followed by reduction to `deg u ≤ 2`.
    pub fn add(&self, a: &Divisor, b: &Divisor) -> Divisor {
        if a.is_identity() {
            return *b;
        }
        if b.is_identity() {
            return *a;
        }
        let k = &self.k;
        let (au, av) = self.polys(a);
        let (bu, bv) = self.polys(b);
        // Composition.
        let (d1, e1, e2) = pxgcd(k, &au, &bu);
        let vs = padd(k, &av, &bv);
        let (d, c1, c2) = pxgcd(k, &d1, &vs);
        let s1 = pmul(k, &c1, &e1);
        let s2 = pmul(k, &c1, &e2);
        let (u, _) = pdivrem(k, &pmul(k, &au, &bu), &pmul(k, &d, &d));
        let t1 = pmul(k, &pmul(k, &s1, &au), &bv);
        let t2 = pmul(k, &pmul(k, &s2, &bu), &av);
        let t3 = pmul(k, &c2, &padd(k, &pmul(k, &av, &bv), &self.f));
        let (vq, _) = pdivrem(k, &padd(k, &padd(k, &t1, &t2), &t3), &d);
        let mut u = pmonic(k, &u);
        if u.is_empty() {
            return Divisor::IDENTITY; // only for invalid input
        }
        let mut v = pdivrem(k, &vq, &u).1;
        // Reduction: (u, v) ↦ ((f − v²)/u, −v mod ·) until deg u ≤ 2.
        while u.len() > 3 {
            let (u2, _) = pdivrem(k, &psub(k, &self.f, &pmul(k, &v, &v)), &u);
            u = pmonic(k, &u2);
            if u.is_empty() {
                return Divisor::IDENTITY; // only for invalid input
            }
            let nv: Poly = trim(v.iter().map(|&c| k.neg(c)).collect());
            v = pdivrem(k, &nv, &u).1;
        }
        Divisor::from_polys(&u, &v).unwrap_or(Divisor::IDENTITY)
    }

    /// `2·d`.
    pub fn double(&self, d: &Divisor) -> Divisor {
        self.add(d, d)
    }

    /// `a − b`.
    pub fn sub(&self, a: &Divisor, b: &Divisor) -> Divisor {
        self.add(a, &self.neg(b))
    }

    /// `[n]·d` by double-and-add.
    pub fn mul(&self, d: &Divisor, mut n: u64) -> Divisor {
        let mut acc = Divisor::IDENTITY;
        let mut base = *d;
        while n > 0 {
            if n & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            n >>= 1;
            if n > 0 {
                base = self.double(&base);
            }
        }
        acc
    }

    /// `[n]·d` for a signed `n`.
    pub fn mul_signed(&self, d: &Divisor, n: i64) -> Divisor {
        let m = self.mul(d, n.unsigned_abs());
        if n < 0 {
            self.neg(&m)
        } else {
            m
        }
    }

    // ── ζ-action ───────────────────────────────────────────────────────────

    /// The primitive fifth root of unity `ζ₅ ∈ F_p` used by [`Self::zeta`]:
    /// `g^((p−1)/5)` for the least `g ≥ 2` giving a value other than 1.
    /// `None` unless `f = c·x⁵ + b` and `p ≡ 1 (mod 5)`.
    pub fn zeta5(&self) -> Option<u64> {
        self.zeta5
    }

    /// `ζ : (x, y) ↦ (ζ₅x, y)` on divisor classes, an automorphism of order 5
    /// of `y² = c·x⁵ + b`. On Mumford pairs it is `u(x) ↦ ζ^{deg u}·u(x/ζ)`,
    /// `v(x) ↦ v(x/ζ)`.
    pub fn zeta(&self, d: &Divisor) -> Result<Divisor, Genus2Error> {
        self.zeta_pow(d, 1)
    }

    /// `ζ^e` applied to `d` (`e` is taken modulo 5).
    pub fn zeta_pow(&self, d: &Divisor, e: u32) -> Result<Divisor, Genus2Error> {
        let z = self.zeta5.ok_or(Genus2Error::NoZetaAction)?;
        let k = &self.k;
        let z = k.pow(z, u64::from(e % 5));
        let zi = k.inv(z);
        let deg = d.degree().min(2);
        let mut out = Divisor {
            deg: deg as u8,
            ..Divisor::IDENTITY
        };
        for i in 0..deg {
            out.u[i] = k.mul(d.u[i] % k.p, k.pow(z, (deg - i) as u64));
            out.v[i] = k.mul(d.v[i] % k.p, k.pow(zi, i as u64));
        }
        Ok(out)
    }

    /// The order-10 automorphism `α = ι∘ζ` (with `ι` the hyperelliptic
    /// involution, i.e. negation) raised to the power `e`: `α^e = ι^e ζ^e`.
    /// `α⁵ = ι` and `α` generates the cyclic group `⟨ζ, ι⟩ ≅ Z/10`.
    pub fn automorphism(&self, d: &Divisor, e: u32) -> Result<Divisor, Genus2Error> {
        let z = self.zeta_pow(d, e % 5)?;
        Ok(if e % 2 == 1 { self.neg(&z) } else { z })
    }

    // ── Solving for v ──────────────────────────────────────────────────────

    /// Every reduced class whose first Mumford polynomial is the monic `u`
    /// (coefficients lowest degree first, `deg u ≤ 2`): at most 1, 2 or 4
    /// classes for `deg u = 0, 1, 2`. Each returned pair is checked against
    /// `u | f − v²`.
    pub fn divisors_with_u(&self, u: &[u64]) -> Result<Vec<Divisor>, Genus2Error> {
        let u = trim(u.to_vec());
        if u.is_empty() || u.len() > 3 || u.last() != Some(&1) || u.iter().any(|&c| c >= self.k.p) {
            return Err(Genus2Error::InvalidDivisor);
        }
        Ok(self.solve(&u))
    }

    /// Solutions for a monic, reduced, trimmed `u` of degree ≤ 2.
    fn solve(&self, u: &[u64]) -> Vec<Divisor> {
        let k = &self.k;
        let mut cand: Vec<Poly> = Vec::new();
        match u.len() {
            1 => return vec![Divisor::IDENTITY],
            2 => {
                let a = k.neg(u[0]);
                if let Some(s) = k.sqrt(self.eval_f(a)) {
                    cand.push(trim(vec![s]));
                    cand.push(trim(vec![k.neg(s)]));
                }
            }
            _ => {
                let (a0, a1) = (u[0], u[1]);
                let half = k.inv(2);
                let disc = k.sub(k.mul(a1, a1), k.mul(4, a0));
                if disc == 0 {
                    // u = (x − e)²: v(e) = s, v'(e) = f'(e) / (2s), s ≠ 0.
                    let e = k.mul(k.neg(a1), half);
                    let fe = self.eval_f(e);
                    if fe != 0 {
                        if let Some(s) = k.sqrt(fe) {
                            let fd = peval(k, &pderiv(k, &self.f), e);
                            for s in [s, k.neg(s)] {
                                let v1 = k.mul(fd, k.inv(k.mul(2, s)));
                                let v0 = k.sub(s, k.mul(v1, e));
                                cand.push(trim(vec![v0, v1]));
                            }
                        }
                    }
                } else if let Some(sq) = k.sqrt(disc) {
                    // Two rational roots: interpolate v through (e_i, ±s_i).
                    let e1 = k.mul(k.add(k.neg(a1), sq), half);
                    let e2 = k.mul(k.sub(k.neg(a1), sq), half);
                    if let (Some(s1), Some(s2)) = (k.sqrt(self.eval_f(e1)), k.sqrt(self.eval_f(e2)))
                    {
                        let di = k.inv(k.sub(e1, e2));
                        for y1 in [s1, k.neg(s1)] {
                            for y2 in [s2, k.neg(s2)] {
                                let v1 = k.mul(k.sub(y1, y2), di);
                                let v0 = k.sub(y1, k.mul(v1, e1));
                                cand.push(trim(vec![v0, v1]));
                            }
                        }
                    }
                } else {
                    // Conjugate roots θ, θ̄ ∈ F_{p²}: θ = (−a1 + κ·t)/2 with
                    // κ² = disc / w; v(θ) = s with s² = f(θ).
                    if let Some(kap) = k.sqrt(k.mul(disc, k.inv(k.w))) {
                        let th = (k.mul(k.neg(a1), half), k.mul(kap, half));
                        let fth = peval2(k, &self.f, th);
                        let roots: Vec<Fp2> = if fth == (0, 0) {
                            vec![(0, 0)]
                        } else {
                            match k.sqrt2(fth) {
                                Some(s) => vec![s, (k.neg(s.0), k.neg(s.1))],
                                None => Vec::new(),
                            }
                        };
                        let ti = k.inv(th.1);
                        for s in roots {
                            let v1 = k.mul(s.1, ti);
                            let v0 = k.sub(s.0, k.mul(v1, th.0));
                            cand.push(trim(vec![v0, v1]));
                        }
                    }
                }
            }
        }
        let mut out: Vec<Divisor> = Vec::with_capacity(cand.len());
        for v in cand {
            if let Some(d) = Divisor::from_polys(u, &v) {
                if self.is_valid(&d) && !out.contains(&d) {
                    out.push(d);
                }
            }
        }
        out
    }

    /// A uniformly random element of `J(F_p)`: draw a monic `u` of degree
    /// ≤ 2 uniformly from the `p² + p + 1` candidates and a slot in `0..4`,
    /// and accept the class in that slot if `u` has that many. Every class
    /// is drawn with probability `1 / (4(p² + p + 1))` per trial; the
    /// expected number of trials is `4(p² + p + 1) / #J ≈ 4`.
    pub fn random<R: Rng + ?Sized>(&self, rng: &mut R) -> Divisor {
        let p = self.k.p;
        let total = p * p + p + 1;
        loop {
            let code = rng.gen_range(0..total);
            let u: Poly = if code == 0 {
                vec![1]
            } else if code <= p {
                vec![code - 1, 1]
            } else {
                let c = code - p - 1;
                vec![c % p, c / p, 1]
            };
            let slot = rng.gen_range(0..4usize);
            let sols = self.solve(&u);
            if let Some(d) = sols.get(slot) {
                return *d;
            }
        }
    }

    // ── Counting ───────────────────────────────────────────────────────────

    fn guard(&self, limit: u64, what: &'static str) -> Result<(), Genus2Error> {
        if self.k.p > limit {
            Err(Genus2Error::PrimeTooLargeFor {
                p: self.k.p,
                limit,
                what,
            })
        } else {
            Ok(())
        }
    }

    /// Quadratic character table of `F_p`: `chi[a] ∈ {0, 1, −1}`.
    fn chi_table(&self) -> Vec<i8> {
        let p = self.k.p;
        let mut chi = vec![-1i8; p as usize];
        chi[0] = 0;
        for x in 1..=(p - 1) / 2 {
            chi[(x * x % p) as usize] = 1;
        }
        chi
    }

    /// `#C(F_p)`, including the single point at infinity. Needs
    /// `p ≤` [`MAX_POINT_COUNT_PRIME`].
    pub fn count_points_fp(&self) -> Result<u64, Genus2Error> {
        self.guard(MAX_POINT_COUNT_PRIME, "count_points_fp")?;
        let chi = self.chi_table();
        Ok(1 + (0..self.k.p)
            .map(|x| (1 + i64::from(chi[self.eval_f(x) as usize])) as u64)
            .sum::<u64>())
    }

    /// `#C(F_{p²})`, including the single point at infinity, by evaluating
    /// `f` at all `p²` abscissae. Needs `p ≤` [`MAX_POINT_COUNT_PRIME`].
    pub fn count_points_fp2(&self) -> Result<u64, Genus2Error> {
        self.guard(MAX_POINT_COUNT_PRIME, "count_points_fp2")?;
        let k = &self.k;
        let chi = self.chi_table();
        let mut n = 1u64;
        for x0 in 0..k.p {
            for x1 in 0..k.p {
                let fx = peval2(k, &self.f, (x0, x1));
                // An element of F_{p²} is a square iff its norm is a square
                // in F_p.
                n += if fx == (0, 0) {
                    1
                } else {
                    (1 + i64::from(chi[k.norm2(fx) as usize])) as u64
                };
            }
        }
        Ok(n)
    }

    /// The `L`-polynomial from `#C(F_p)` and `#C(F_{p²})`. Needs
    /// `p ≤` [`MAX_POINT_COUNT_PRIME`].
    pub fn l_polynomial(&self) -> Result<LPolynomial, Genus2Error> {
        let n1 = self.count_points_fp()?;
        let n2 = self.count_points_fp2()?;
        let p = self.k.p as i64;
        let a1 = n1 as i64 - p - 1;
        let a2 = (n2 as i64 - p * p - 1 + a1 * a1) / 2;
        Ok(LPolynomial {
            p: self.k.p,
            n1,
            n2,
            a1,
            a2,
        })
    }

    /// `#J(F_p) = L(1)`. Needs `p ≤` [`MAX_POINT_COUNT_PRIME`].
    pub fn jacobian_order(&self) -> Result<u64, Genus2Error> {
        Ok(self.l_polynomial()?.jacobian_order())
    }

    /// Calls `visit` once on every element of `J(F_p)`: the identity, then
    /// every class with `deg u = 1`, then every class with `deg u = 2`, each
    /// group in increasing order of the coefficients of `u`. Needs
    /// `p ≤` [`MAX_ENUMERATION_PRIME`].
    pub fn for_each_element<F: FnMut(&Divisor)>(&self, mut visit: F) -> Result<(), Genus2Error> {
        self.guard(MAX_ENUMERATION_PRIME, "for_each_element")?;
        let p = self.k.p;
        visit(&Divisor::IDENTITY);
        for c in 0..p {
            for d in self.solve(&[c, 1]) {
                visit(&d);
            }
        }
        for a1 in 0..p {
            for a0 in 0..p {
                for d in self.solve(&[a0, a1, 1]) {
                    visit(&d);
                }
            }
        }
        Ok(())
    }

    /// `#J(F_p)` by exhaustive enumeration of reduced Mumford pairs,
    /// independent of the point counts. Needs `p ≤` [`MAX_ENUMERATION_PRIME`].
    pub fn count_exhaustive(&self) -> Result<u64, Genus2Error> {
        let mut n = 0u64;
        self.for_each_element(|_| n += 1)?;
        Ok(n)
    }

    /// Every element of `J(F_p)`, in the order of
    /// [`Self::for_each_element`]. Needs `p ≤` [`MAX_LIST_PRIME`].
    pub fn enumerate(&self) -> Result<Vec<Divisor>, Genus2Error> {
        self.guard(MAX_LIST_PRIME, "enumerate")?;
        let mut out = Vec::new();
        self.for_each_element(|d| out.push(*d))?;
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::jv_cover::{Div as JvDiv, Fld, Fq, Hyp, E2};
    use crate::prime_hyperelliptic::{
        brute_force_jac_order_via_lpoly, FpPoly, HyperellipticCurveP, MumfordDivisorP,
    };
    use num_bigint::BigUint;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use std::collections::HashSet;

    /// The x⁵ + b curves, and generic quintics, used across the tests.
    fn curves() -> Vec<Genus2Curve> {
        let mut out = Vec::new();
        for (p, b) in [(31, 3), (61, 2), (101, 7), (151, 5), (41, 1)] {
            out.push(Genus2Curve::x5_plus_b(p, b).unwrap());
        }
        // Generic quintics, including p ≢ 1 (mod 5) and a non-monic one.
        for (p, f) in [
            (37u64, [1u64, 2, 0, 3, 0, 1]),
            (43, [5, 0, 1, 0, 7, 1]),
            (101, [17, 3, 9, 0, 2, 4]),
        ] {
            out.push(Genus2Curve::new(p, f).unwrap());
        }
        out
    }

    #[test]
    fn genus2_construction_errors() {
        assert_eq!(
            Genus2Curve::new(35, [1, 0, 0, 0, 0, 1]).unwrap_err(),
            Genus2Error::PrimeOutOfRange { p: 35 }
        );
        assert_eq!(
            Genus2Curve::new(2, [1, 0, 0, 0, 0, 1]).unwrap_err(),
            Genus2Error::PrimeOutOfRange { p: 2 }
        );
        assert!(matches!(
            Genus2Curve::new((1 << 31) + 11, [1, 0, 0, 0, 0, 1]),
            Err(Genus2Error::PrimeOutOfRange { .. })
        ));
        assert_eq!(
            Genus2Curve::new(31, [1, 0, 0, 0, 1, 0]).unwrap_err(),
            Genus2Error::NotQuintic
        );
        assert_eq!(
            Genus2Curve::new(31, [40, 0, 0, 0, 0, 1]).unwrap_err(),
            Genus2Error::CoefficientNotReduced { value: 40 }
        );
        // x⁵ + b in characteristic 5, and x²·(x³ + 1), are singular.
        assert_eq!(
            Genus2Curve::x5_plus_b(5, 1).unwrap_err(),
            Genus2Error::NotSquarefree
        );
        assert_eq!(
            Genus2Curve::new(31, [0, 0, 1, 0, 0, 1]).unwrap_err(),
            Genus2Error::NotSquarefree
        );
        assert_eq!(
            Genus2Curve::x5_plus_b(31, 0).unwrap_err(),
            Genus2Error::NotSquarefree
        );
        assert_eq!(
            Genus2Curve::x5_plus_b(0, 1).unwrap_err(),
            Genus2Error::PrimeOutOfRange { p: 0 }
        );
        // A large prime is accepted, but the costly routines refuse it.
        let big = Genus2Curve::x5_plus_b(2_147_483_647, 3).unwrap();
        assert!(matches!(
            big.jacobian_order(),
            Err(Genus2Error::PrimeTooLargeFor { .. })
        ));
        assert!(matches!(
            big.count_exhaustive(),
            Err(Genus2Error::PrimeTooLargeFor { .. })
        ));
        let mut rng = StdRng::seed_from_u64(9);
        let d = big.random(&mut rng);
        assert!(big.is_valid(&d));
        let e = big.random(&mut rng);
        let s = big.add(&d, &e);
        assert!(big.is_valid(&s));
        assert_eq!(big.sub(&s, &e), d);
    }

    #[test]
    fn genus2_divisor_validation() {
        let c = Genus2Curve::x5_plus_b(31, 3).unwrap();
        assert_eq!(c.divisor(&[1], &[]).unwrap(), Divisor::IDENTITY);
        assert_eq!(c.divisor(&[2, 1], &[]), Err(Genus2Error::InvalidDivisor));
        assert_eq!(c.divisor(&[0, 2], &[]), Err(Genus2Error::InvalidDivisor));
        assert_eq!(
            c.divisor(&[0, 1], &[1, 1]),
            Err(Genus2Error::InvalidDivisor)
        );
        assert_eq!(
            c.divisor(&[0, 0, 0, 1], &[]),
            Err(Genus2Error::InvalidDivisor)
        );
        assert_eq!(c.divisor(&[99, 1], &[]), Err(Genus2Error::InvalidDivisor));
        assert!(c.divisors_with_u(&[0, 0, 0, 1]).is_err());
        // Points: f(x) = x⁵ + 3, find one and round-trip it.
        let (x, y) = (0..31)
            .flat_map(|x| (0..31).map(move |y| (x, y)))
            .find(|&(x, y)| c.is_on_curve(x, y))
            .unwrap();
        let d = c.point(x, y).unwrap();
        assert_eq!(d.degree(), 1);
        assert_eq!(d.u_coeffs(), vec![(31 - x) % 31, 1]);
        assert_eq!(c.divisor(&d.u_coeffs(), &d.v_coeffs()).unwrap(), d);
        assert!(c.point(x, (y + 1) % 31).is_err() || c.is_on_curve(x, (y + 1) % 31));
        // Operations on a foreign divisor do not panic.
        let other = Genus2Curve::x5_plus_b(151, 5).unwrap();
        let mut rng = StdRng::seed_from_u64(2);
        let foreign = other.random(&mut rng);
        let _ = c.add(&foreign, &d);
        let _ = c.mul(&foreign, 1000);
        let _ = c.zeta(&foreign);
    }

    #[test]
    fn genus2_enumeration_equals_l1() {
        for c in curves() {
            let lp = c.l_polynomial().unwrap();
            let all = c.enumerate().unwrap();
            assert_eq!(all.len() as u64, lp.jacobian_order(), "p = {}", c.p());
            assert_eq!(c.count_exhaustive().unwrap(), lp.jacobian_order());
            assert_eq!(lp.eval(1), lp.jacobian_order() as i128);
            let set: HashSet<Divisor> = all.iter().copied().collect();
            assert_eq!(set.len(), all.len(), "duplicates at p = {}", c.p());
            assert!(all.iter().all(|d| c.is_valid(d)));
            // Weil bounds on the trace.
            assert!((lp.a1.unsigned_abs() as f64) <= 4.0 * (c.p() as f64).sqrt());
        }
    }

    #[test]
    fn genus2_l1_matches_prime_hyperelliptic() {
        for c in curves() {
            let p = BigUint::from(c.p());
            let f =
                FpPoly::from_coeffs(c.f().iter().map(|&x| BigUint::from(x)).collect(), p.clone());
            // prime_hyperelliptic wants f monic; skip the non-monic curve.
            if c.f()[5] != 1 {
                continue;
            }
            let lib = HyperellipticCurveP::try_new(p, f, 2).unwrap();
            let order = brute_force_jac_order_via_lpoly(&lib);
            assert_eq!(order, BigUint::from(c.jacobian_order().unwrap()));
        }
    }

    /// The solver-based enumeration agrees with naive search over every
    /// `(u, v)` pair at a tiny prime.
    #[test]
    fn genus2_enumeration_matches_naive_search() {
        for (p, f) in [(11u64, [1u64, 0, 0, 0, 0, 1]), (13, [2, 1, 0, 0, 3, 1])] {
            let c = Genus2Curve::new(p, f).unwrap();
            let mut naive: HashSet<Divisor> = HashSet::new();
            naive.insert(Divisor::IDENTITY);
            for a in 0..p {
                for v0 in 0..p {
                    if let Ok(d) = c.divisor(&[a, 1], &[v0]) {
                        naive.insert(d);
                    }
                }
            }
            for a0 in 0..p {
                for a1 in 0..p {
                    for v0 in 0..p {
                        for v1 in 0..p {
                            if let Ok(d) = c.divisor(&[a0, a1, 1], &[v0, v1]) {
                                naive.insert(d);
                            }
                        }
                    }
                }
            }
            let ours: HashSet<Divisor> = c.enumerate().unwrap().into_iter().collect();
            assert_eq!(ours, naive, "p = {p}");
        }
    }

    #[test]
    fn genus2_group_laws() {
        let mut rng = StdRng::seed_from_u64(1);
        for c in curves() {
            let order = c.jacobian_order().unwrap();
            for _ in 0..40 {
                let a = c.random(&mut rng);
                let b = c.random(&mut rng);
                let d = c.random(&mut rng);
                assert!(c.is_valid(&a));
                assert_eq!(c.add(&a, &b), c.add(&b, &a));
                assert_eq!(
                    c.add(&c.add(&a, &b), &d),
                    c.add(&a, &c.add(&b, &d)),
                    "associativity at p = {}",
                    c.p()
                );
                assert!(c.add(&a, &c.neg(&a)).is_identity());
                assert_eq!(c.double(&a), c.add(&a, &a));
                assert_eq!(c.add(&a, &Divisor::IDENTITY), a);
                assert!(c.mul(&a, order).is_identity());
                assert_eq!(c.mul(&a, order + 1), a);
                assert_eq!(c.mul_signed(&a, -3), c.neg(&c.mul(&a, 3)));
                assert_eq!(c.mul(&a, 5), c.add(&c.double(&c.double(&a)), &a));
            }
        }
    }

    /// Cantor here agrees with `jv_cover::Hyp` (over `F_{p²} ⊃ F_p`) and with
    /// `prime_hyperelliptic::MumfordDivisorP`.
    #[test]
    fn genus2_add_matches_existing_implementations() {
        let mut rng = StdRng::seed_from_u64(7);
        for c in curves().into_iter().filter(|c| c.f()[5] == 1) {
            let fq = Fq::new(c.p());
            let lift = |a: &[u64]| -> Vec<E2> { a.iter().map(|&x| fq.from_fp(x)).collect() };
            let hyp = Hyp {
                f: &fq,
                h: lift(&c.f()),
                g: 2,
            };
            let norm = |d: &JvDiv<E2>| -> (Vec<E2>, Vec<E2>) {
                let t = |a: &Vec<E2>| {
                    let mut a = a.clone();
                    while a.last() == Some(&fq.zero()) {
                        a.pop();
                    }
                    a
                };
                (t(&d.u), t(&d.v))
            };
            let pb = BigUint::from(c.p());
            let big = HyperellipticCurveP::try_new(
                pb.clone(),
                FpPoly::from_coeffs(
                    c.f().iter().map(|&x| BigUint::from(x)).collect(),
                    pb.clone(),
                ),
                2,
            )
            .unwrap();
            let to_big = |d: &Divisor| MumfordDivisorP {
                u: FpPoly::from_coeffs(
                    d.u_coeffs().into_iter().map(BigUint::from).collect(),
                    pb.clone(),
                ),
                v: FpPoly::from_coeffs(
                    d.v_coeffs().into_iter().map(BigUint::from).collect(),
                    pb.clone(),
                ),
            };
            for i in 0..60 {
                let a = c.random(&mut rng);
                // Exercise doubling and the inverse pair, too.
                let b = match i % 6 {
                    0 => a,
                    1 => c.neg(&a),
                    _ => c.random(&mut rng),
                };
                let s = c.add(&a, &b);
                let ja = JvDiv {
                    u: lift(&a.u_coeffs()),
                    v: lift(&a.v_coeffs()),
                };
                let jb = JvDiv {
                    u: lift(&b.u_coeffs()),
                    v: lift(&b.v_coeffs()),
                };
                let js = hyp.add(&ja, &jb);
                let mine = JvDiv {
                    u: lift(&s.u_coeffs()),
                    v: lift(&s.v_coeffs()),
                };
                assert_eq!(norm(&js), norm(&mine), "jv_cover at p = {}", c.p());
                let bs = to_big(&a).checked_add(&to_big(&b), &big).unwrap();
                assert_eq!(bs, to_big(&s), "prime_hyperelliptic at p = {}", c.p());
            }
        }
    }

    #[test]
    fn genus2_zeta_action() {
        let mut rng = StdRng::seed_from_u64(3);
        for c in curves() {
            let Some(z) = c.zeta5() else {
                assert!(matches!(
                    c.zeta(&Divisor::IDENTITY),
                    Err(Genus2Error::NoZetaAction)
                ));
                continue;
            };
            let p = c.p();
            assert_eq!(p % 5, 1);
            assert_ne!(z, 1);
            assert_eq!(c.k.pow(z, 5), 1);
            let mut moved = false;
            for _ in 0..60 {
                let a = c.random(&mut rng);
                let b = c.random(&mut rng);
                let za = c.zeta(&a).unwrap();
                assert!(c.is_valid(&za));
                assert_eq!(
                    c.zeta(&c.add(&a, &b)).unwrap(),
                    c.add(&za, &c.zeta(&b).unwrap()),
                    "ζ is not a homomorphism at p = {p}"
                );
                assert_eq!(c.zeta_pow(&a, 5).unwrap(), a);
                assert_eq!(c.zeta_pow(&za, 4).unwrap(), a);
                moved |= za != a;
                // α = ι∘ζ: α⁵ = −1, α^10 = id, α² = ζ².
                assert_eq!(c.automorphism(&a, 5).unwrap(), c.neg(&a));
                assert_eq!(c.automorphism(&a, 10).unwrap(), a);
                assert_eq!(c.automorphism(&a, 2).unwrap(), c.zeta_pow(&a, 2).unwrap());
                let mut it = a;
                for _ in 0..3 {
                    it = c.automorphism(&it, 1).unwrap();
                }
                assert_eq!(it, c.automorphism(&a, 3).unwrap());
            }
            assert!(moved, "ζ acts trivially at p = {p}");
            // ζ permutes J(F_p).
            if p <= 61 {
                let all: HashSet<Divisor> = c.enumerate().unwrap().into_iter().collect();
                let image: HashSet<Divisor> = all.iter().map(|d| c.zeta(d).unwrap()).collect();
                assert_eq!(all, image);
            }
        }
    }

    #[test]
    fn genus2_point_counts_by_naive_search() {
        // #C(F_p) by trying every (x, y), against the character sum.
        for c in curves() {
            let p = c.p();
            let naive = 1
                + (0..p)
                    .flat_map(|x| (0..p).map(move |y| (x, y)))
                    .filter(|&(x, y)| c.is_on_curve(x, y))
                    .count() as u64;
            assert_eq!(c.count_points_fp().unwrap(), naive);
        }
    }

    #[test]
    fn genus2_random_is_roughly_uniform() {
        // At p = 11 every element should appear in 40·#J draws.
        let c = Genus2Curve::x5_plus_b(11, 1).unwrap();
        let n = c.jacobian_order().unwrap();
        let mut rng = StdRng::seed_from_u64(5);
        let seen: HashSet<Divisor> = (0..40 * n).map(|_| c.random(&mut rng)).collect();
        assert_eq!(seen.len() as u64, n);
    }
}
