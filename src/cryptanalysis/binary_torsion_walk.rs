//! # Walking a whole Koblitz isogeny class by building `E[ℓ]` where it lives
//!
//! [`crate::cryptanalysis::binary_velu`] computes an isogeny once its kernel is
//! known, and found kernels only at `ℓ = 3`, where the kernel polynomial is
//! linear.  For larger `ℓ` it named the obstacle as factoring `ψ_ℓ` — degree
//! 36720 at `n = 17` and 104424 at `n = 19`.  This module does not factor
//! anything.  It builds the `ℓ`-torsion directly, in the field it is defined
//! over, and reads the kernels off points.
//!
//! ## Why that field is small
//!
//! The Koblitz curve `K_a` is defined over `F_2`, so the small Frobenius `τ`
//! (`τ² ∓ τ + 2 = 0`) is an endomorphism and `End(K_a) = Z[τ] = O_K` is the
//! **maximal** order.  `K_a` is the crater of its volcano, not the floor.  The
//! big Frobenius is `π = τ^n` and `c = [O_K : Z[π]]` is the conductor, so for
//! every `ℓ | c`
//!
//! ```text
//!   π ∈ Z + ℓ·O_K   ⇒   π acts on E[ℓ] as a scalar λ,
//!   2λ ≡ t (mod ℓ),   λ² ≡ q (mod ℓ).
//! ```
//!
//! A scalar Frobenius fixes every line of `E[ℓ]`, so **all `ℓ + 1` subgroups
//! are rational**, and `π^m = λ^m = 1` on `E[ℓ]` with `m = ord_ℓ(λ)` puts the
//! whole of `E[ℓ]` inside `E(F_{q^m})`.  For the classes the cost sweep
//! measured, `m` is 135 at `(n, ℓ) = (17, 271)` and 8 at `(19, 457)` —
//! `E[457]` of a curve over `F_{2^19}` lives in `F_{2^152}`.
//!
//! The same argument applies one level down.  Below the crater `End(E)` has
//! conductor `f | c`, and `π` is scalar on `E[ℓ]` exactly when `ℓ | c/f`, so
//! every **descending** edge of the volcano comes from a vertex where this
//! construction applies.  A breadth-first walk of descending edges from
//! `K_a` therefore reaches the whole class, and the walk tests that claim by
//! comparing what it reached against the exhaustive trace-scan census.
//!
//! ## What is computed, and what is checked
//!
//! For each vertex and each `ℓ | c`:
//!
//! 1. Find two independent points of order `ℓ` in `E(F_{q^m})` by clearing
//!    the cofactor off random points — an `x`-only Montgomery ladder
//!    (López–Dahab) with the group order from the trace recurrence.  Finding
//!    only one line is recorded as that vertex being on the `ℓ`-floor, which
//!    is the volcano prediction for vertices reached by an `ℓ`-edge.
//! 2. List the `ℓ + 1` kernels `⟨P₁⟩, ⟨P₂ + kP₁⟩` and, for each, the
//!    abscissae `x(jG)`, `j = 1..(ℓ−1)/2`, by `x`-only differential addition.
//! 3. Apply Vélu from [`crate::cryptanalysis::binary_velu`]:
//!    `t = Σ x(jG)`, `a₆' = a₆ + t + t²`, `X = x + A + A²`.
//!
//! Every one of those is computed in `F_{q^m}`, and every one **must land in
//! `F_q`** — `t` because the kernel is rational, `X` because the isogeny is.
//! Nothing forces that except the mathematics, so it is checked on every
//! kernel and every transported point, and a failure is counted rather than
//! skipped.  On top of that: the codomain's order must equal `#E`, the
//! transported `φ(P)` must have order `r`, and `[d]φ(P)` must sit at
//! `x(φ(Q))` on the codomain.

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use rayon::prelude::*;

use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::binary_ecc::IrreduciblePoly;
use crate::cryptanalysis::binary_velu::Curve;
use crate::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crate::cryptanalysis::koblitz_fast::FastPoint;
use crate::cryptanalysis::semaev_decomp::Gf2;

// ── polynomials over F_q, coefficients as words ─────────────────────

fn trim(v: &mut Vec<u64>) {
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    if v.is_empty() {
        v.push(0);
    }
}

fn deg(v: &[u64]) -> Option<usize> {
    (0..v.len()).rev().find(|&i| v[i] != 0)
}

fn padd(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = vec![0u64; a.len().max(b.len())];
    for (i, x) in a.iter().enumerate() {
        out[i] ^= x;
    }
    for (i, x) in b.iter().enumerate() {
        out[i] ^= x;
    }
    trim(&mut out);
    out
}

fn pmul(gf: &Gf2, a: &[u64], b: &[u64]) -> Vec<u64> {
    let (Some(da), Some(db)) = (deg(a), deg(b)) else {
        return vec![0];
    };
    let mut out = vec![0u64; da + db + 1];
    for i in 0..=da {
        if a[i] == 0 {
            continue;
        }
        for j in 0..=db {
            if b[j] != 0 {
                out[i + j] ^= gf.mul(a[i], b[j]);
            }
        }
    }
    out
}

/// `a mod f`, for any nonzero `f`.
fn pmod(gf: &Gf2, a: &[u64], f: &[u64]) -> Vec<u64> {
    let df = deg(f).expect("nonzero modulus");
    let lead_inv = gf.inv(f[df]);
    let mut r = a.to_vec();
    trim(&mut r);
    while let Some(dr) = deg(&r) {
        if dr < df {
            break;
        }
        let c = gf.mul(r[dr], lead_inv);
        let shift = dr - df;
        for j in 0..=df {
            if f[j] != 0 {
                r[shift + j] ^= gf.mul(c, f[j]);
            }
        }
    }
    trim(&mut r);
    r
}

fn pgcd(gf: &Gf2, a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    trim(&mut a);
    trim(&mut b);
    while deg(&b).is_some() {
        let r = pmod(gf, &a, &b);
        a = b;
        b = r;
    }
    a
}

/// `a⁻¹ mod f` by the extended Euclidean algorithm, or `None` if not a unit.
fn pinv(gf: &Gf2, a: &[u64], f: &[u64]) -> Option<Vec<u64>> {
    // Invariants: s_i · a ≡ r_i (mod f).
    let mut r0 = f.to_vec();
    let mut r1 = pmod(gf, a, f);
    let mut s0: Vec<u64> = vec![0];
    let mut s1: Vec<u64> = vec![1];
    while let Some(d1) = deg(&r1) {
        // one division step r0 = q·r1 + r
        let lead_inv = gf.inv(r1[d1]);
        let mut rem = r0.clone();
        trim(&mut rem);
        let mut quo = vec![0u64; rem.len().max(1)];
        while let Some(dr) = deg(&rem) {
            if dr < d1 {
                break;
            }
            let c = gf.mul(rem[dr], lead_inv);
            let sh = dr - d1;
            quo[sh] ^= c;
            for j in 0..=d1 {
                if r1[j] != 0 {
                    rem[sh + j] ^= gf.mul(c, r1[j]);
                }
            }
        }
        trim(&mut rem);
        trim(&mut quo);
        let s2 = padd(&s0, &pmul(gf, &quo, &s1));
        r0 = r1;
        r1 = rem;
        s0 = s1;
        s1 = s2;
    }
    // r0 = gcd; unit iff degree 0
    if deg(&r0) != Some(0) {
        return None;
    }
    let c = gf.inv(r0[0]);
    let mut out: Vec<u64> = s0.iter().map(|&x| gf.mul(x, c)).collect();
    let mut out = pmod(gf, &std::mem::take(&mut out), f);
    trim(&mut out);
    Some(out)
}

// ── the extension F_{q^m} = F_q[z]/(f) ──────────────────────────────

/// `F_{q^m}` as `F_q[z]/(f)` with `f` monic irreducible of degree `m`.
/// Elements are length-`m` coefficient vectors over `F_q`.
pub struct Ext {
    pub gf: Gf2,
    pub n: u32,
    pub m: usize,
    f: Vec<u64>,
}

pub type El = Vec<u64>;

impl Ext {
    /// Find a monic irreducible of degree `m` over `F_q` by Ben-Or's test and
    /// build the field.  Deterministic in `seed`.
    pub fn new(irr: &IrreduciblePoly, m: usize, seed: u64) -> Self {
        let gf = Gf2::new(irr);
        let n = irr.degree;
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        if m == 1 {
            return Ext {
                gf,
                n,
                m,
                f: vec![0, 1],
            };
        }
        let mut rng = StdRng::seed_from_u64(seed ^ ((m as u64) << 32) ^ n as u64);
        loop {
            let mut f: Vec<u64> = (0..m).map(|_| rng.gen::<u64>() & mask).collect();
            f.push(1);
            if f[0] == 0 {
                continue;
            }
            if Self::ben_or_irreducible(&gf, n, &f) {
                return Ext { gf, n, m, f };
            }
        }
    }

    fn ben_or_irreducible(gf: &Gf2, n: u32, f: &[u64]) -> bool {
        let m = f.len() - 1;
        let z: Vec<u64> = vec![0, 1];
        let mut h = z.clone();
        for _ in 1..=m / 2 {
            // h ← h^q mod f, q = 2^n: n squarings.
            for _ in 0..n {
                let mut s = vec![0u64; 2 * h.len()];
                for (i, &c) in h.iter().enumerate() {
                    s[2 * i] = gf.sqr(c);
                }
                h = pmod(gf, &s, f);
            }
            let g = pgcd(gf, f, &padd(&h, &z));
            if deg(&g).unwrap_or(0) > 0 {
                return false;
            }
        }
        true
    }

    pub fn zero(&self) -> El {
        vec![0; self.m]
    }
    pub fn one(&self) -> El {
        self.base(1)
    }
    pub fn base(&self, c: u64) -> El {
        let mut v = self.zero();
        v[0] = c;
        v
    }
    pub fn is_zero(&self, a: &El) -> bool {
        a.iter().all(|&x| x == 0)
    }
    /// Whether `a` lies in the base field `F_q`.
    pub fn is_base(&self, a: &El) -> bool {
        a[1..].iter().all(|&x| x == 0)
    }
    pub fn add(&self, a: &El, b: &El) -> El {
        a.iter().zip(b).map(|(x, y)| x ^ y).collect()
    }
    fn reduce(&self, mut v: Vec<u64>) -> El {
        let m = self.m;
        for i in (m..v.len()).rev() {
            let c = v[i];
            if c == 0 {
                continue;
            }
            v[i] = 0;
            let sh = i - m;
            for j in 0..m {
                if self.f[j] != 0 {
                    v[sh + j] ^= self.gf.mul(c, self.f[j]);
                }
            }
        }
        v.truncate(m);
        v.resize(m, 0);
        v
    }
    pub fn mul(&self, a: &El, b: &El) -> El {
        let m = self.m;
        let mut out = vec![0u64; 2 * m - 1];
        for i in 0..m {
            if a[i] == 0 {
                continue;
            }
            for j in 0..m {
                if b[j] != 0 {
                    out[i + j] ^= self.gf.mul(a[i], b[j]);
                }
            }
        }
        self.reduce(out)
    }
    pub fn sqr(&self, a: &El) -> El {
        let m = self.m;
        let mut out = vec![0u64; 2 * m - 1];
        for i in 0..m {
            out[2 * i] = self.gf.sqr(a[i]);
        }
        self.reduce(out)
    }
    pub fn inv(&self, a: &El) -> Option<El> {
        if self.is_zero(a) {
            return None;
        }
        if self.m == 1 {
            return Some(vec![self.gf.inv(a[0])]);
        }
        let r = pinv(&self.gf, a, &self.f)?;
        let mut out = r;
        out.resize(self.m, 0);
        Some(out)
    }
    pub fn random(&self, rng: &mut StdRng) -> El {
        let mask = (1u64 << self.n) - 1;
        (0..self.m).map(|_| rng.gen::<u64>() & mask).collect()
    }
    /// Absolute trace `F_{q^m} → F_2`, as 0 or 1.
    pub fn trace2(&self, a: &El) -> u64 {
        let total = self.n as usize * self.m;
        let mut acc = a.clone();
        let mut cur = a.clone();
        for _ in 1..total {
            cur = self.sqr(&cur);
            acc = self.add(&acc, &cur);
        }
        debug_assert!(self.is_base(&acc) && acc[0] <= 1, "trace must land in F_2");
        acc[0] & 1
    }
    /// Solve `w² + w = c` (requires `Tr(c) = 0`).
    pub fn solve_as(&self, c: &El, rng: &mut StdRng) -> Option<El> {
        if self.trace2(c) != 0 {
            return None;
        }
        let total = self.n as usize * self.m;
        if total % 2 == 1 {
            // half-trace
            let mut acc = c.clone();
            let mut cur = c.clone();
            for _ in 0..(total - 1) / 2 {
                cur = self.sqr(&self.sqr(&cur));
                acc = self.add(&acc, &cur);
            }
            return Some(acc);
        }
        // Even degree: w = Σ_{i=0}^{N−2} (Σ_{j=i+1}^{N−1} ω^{2^j}) c^{2^i}, Tr(ω) = 1.
        let omega = loop {
            let w = self.random(rng);
            if self.trace2(&w) == 1 {
                break w;
            }
        };
        let mut wpow = Vec::with_capacity(total);
        let mut cur = omega;
        for _ in 0..total {
            wpow.push(cur.clone());
            cur = self.sqr(&cur);
        }
        // suffix[i] = Σ_{j>i} ω^{2^j}
        let mut suffix = vec![self.zero(); total];
        let mut run = self.zero();
        for i in (0..total).rev() {
            suffix[i] = run.clone();
            run = self.add(&run, &wpow[i]);
        }
        let mut acc = self.zero();
        let mut cp = c.clone();
        for s in suffix.iter().take(total - 1) {
            acc = self.add(&acc, &self.mul(s, &cp));
            cp = self.sqr(&cp);
        }
        Some(acc)
    }
}

// ── the curve y² + xy = x³ + a₂x² + a₆ over F_{q^m} ─────────────────

/// An affine point, or the identity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Pt {
    O,
    A(El, El),
}

pub struct ExtCurve<'a> {
    pub k: &'a Ext,
    pub a2: El,
    pub a6: El,
}

impl<'a> ExtCurve<'a> {
    pub fn new(k: &'a Ext, a2: u64, a6: u64) -> Self {
        ExtCurve {
            k,
            a2: k.base(a2),
            a6: k.base(a6),
        }
    }

    /// `x(2P) = x² + a₆/x²`.
    pub fn x_double(&self, x: &El) -> Option<El> {
        let k = self.k;
        let x2 = k.sqr(x);
        let inv = k.inv(&x2)?;
        Some(k.add(&x2, &k.mul(&self.a6, &inv)))
    }

    /// `x(P+Q) = x(P−Q) + x_P x_Q / (x_P + x_Q)²`.
    pub fn x_diff_add(&self, xp: &El, xq: &El, xdiff: &El) -> Option<El> {
        let k = self.k;
        let s = k.add(xp, xq);
        let inv = k.inv(&k.sqr(&s))?;
        Some(k.add(xdiff, &k.mul(&k.mul(xp, xq), &inv)))
    }

    /// `x([s]P)` by the López–Dahab Montgomery ladder in projective `(X:Z)`;
    /// `None` when the result is the identity.
    pub fn ladder_x(&self, x: &El, s: &BigUint) -> Option<El> {
        let k = self.k;
        if s.is_zero() || k.is_zero(x) && !s.bit(0) {
            return None;
        }
        let dbl = |xx: &El, zz: &El| -> (El, El) {
            let x2 = k.sqr(xx);
            let z2 = k.sqr(zz);
            (
                k.add(&k.sqr(&x2), &k.mul(&self.a6, &k.sqr(&z2))),
                k.mul(&x2, &z2),
            )
        };
        let madd = |x0: &El, z0: &El, x1: &El, z1: &El| -> (El, El) {
            let a = k.mul(x0, z1);
            let b = k.mul(x1, z0);
            let z3 = k.sqr(&k.add(&a, &b));
            (k.add(&k.mul(x, &z3), &k.mul(&a, &b)), z3)
        };
        let (mut x0, mut z0) = (x.clone(), k.one());
        let (mut x1, mut z1) = dbl(x, &k.one());
        let bits = s.bits();
        for i in (0..bits - 1).rev() {
            if s.bit(i) {
                let (ax, az) = madd(&x0, &z0, &x1, &z1);
                x0 = ax;
                z0 = az;
                let (dx, dz) = dbl(&x1, &z1);
                x1 = dx;
                z1 = dz;
            } else {
                let (ax, az) = madd(&x0, &z0, &x1, &z1);
                x1 = ax;
                z1 = az;
                let (dx, dz) = dbl(&x0, &z0);
                x0 = dx;
                z0 = dz;
            }
        }
        if k.is_zero(&z0) {
            return None;
        }
        Some(k.mul(&x0, &k.inv(&z0)?))
    }

    /// A `y` with `(x, y)` on the curve, if one exists: `y = x·w`,
    /// `w² + w = x + a₂ + a₆/x²`.
    pub fn lift_x(&self, x: &El, rng: &mut StdRng) -> Option<El> {
        let k = self.k;
        if k.is_zero(x) {
            return None;
        }
        let inv2 = k.inv(&k.sqr(x))?;
        let c = k.add(&k.add(x, &self.a2), &k.mul(&self.a6, &inv2));
        let w = k.solve_as(&c, rng)?;
        Some(k.mul(x, &w))
    }

    /// A random abscissa that carries a point.
    pub fn random_x(&self, rng: &mut StdRng) -> El {
        let k = self.k;
        loop {
            let x = k.random(rng);
            if k.is_zero(&x) {
                continue;
            }
            let Some(inv2) = k.inv(&k.sqr(&x)) else {
                continue;
            };
            let c = k.add(&k.add(&x, &self.a2), &k.mul(&self.a6, &inv2));
            if k.trace2(&c) == 0 {
                return x;
            }
        }
    }

    pub fn neg(&self, p: &Pt) -> Pt {
        match p {
            Pt::O => Pt::O,
            Pt::A(x, y) => Pt::A(x.clone(), self.k.add(x, y)),
        }
    }

    pub fn double(&self, p: &Pt) -> Pt {
        let k = self.k;
        let Pt::A(x1, y1) = p else { return Pt::O };
        if k.is_zero(x1) {
            return Pt::O;
        }
        let lam = k.add(x1, &k.mul(y1, &k.inv(x1).expect("x ≠ 0")));
        let x3 = k.add(&k.add(&k.sqr(&lam), &lam), &self.a2);
        let y3 = k.add(&k.sqr(x1), &k.mul(&k.add(&lam, &k.one()), &x3));
        Pt::A(x3, y3)
    }

    pub fn add(&self, p: &Pt, q: &Pt) -> Pt {
        let k = self.k;
        let (Pt::A(x1, y1), Pt::A(x2, y2)) = (p, q) else {
            return if *p == Pt::O { q.clone() } else { p.clone() };
        };
        if x1 == x2 {
            if y1 == y2 {
                return self.double(p);
            }
            return Pt::O;
        }
        let lam = k.mul(&k.add(y1, y2), &k.inv(&k.add(x1, x2)).expect("x1 ≠ x2"));
        let x3 = k.add(&k.add(&k.add(&k.sqr(&lam), &lam), &k.add(x1, x2)), &self.a2);
        let y3 = k.add(&k.add(&k.mul(&lam, &k.add(x1, &x3)), &x3), y1);
        Pt::A(x3, y3)
    }
}

// ── group orders over extensions ────────────────────────────────────

/// `#E(F_{q^m})` from the trace `t` of `π`: `s_k = t·s_{k−1} − q·s_{k−2}`.
pub fn order_over_extension(q: u64, t: i64, m: usize) -> BigUint {
    let q = BigInt::from(q);
    let t = BigInt::from(t);
    let (mut s0, mut s1) = (BigInt::from(2), t.clone());
    for _ in 1..m {
        let s2 = &t * &s1 - &q * &s0;
        s0 = s1;
        s1 = s2;
    }
    let qm = q.pow(m as u32);
    (qm + BigInt::one() - s1)
        .to_biguint()
        .expect("a group order is positive")
}

/// `ord_ℓ(λ)` for `λ ≠ 0 mod ℓ`.
pub fn mult_order(lam: u64, ell: u64) -> u64 {
    let mut k = 1u64;
    let mut x = lam % ell;
    while x != 1 {
        x = x * lam % ell;
        k += 1;
    }
    k
}

/// The scalar Frobenius eigenvalue mod `ℓ`, when the arithmetic allows one:
/// `2λ ≡ t` and `λ² ≡ q`.  Necessary for `π` to act as a scalar on `E[ℓ]`;
/// whether it actually does is a property of `End(E)`, tested by
/// [`torsion_basis`] finding two independent lines.
pub fn scalar_eigenvalue(q: u64, t: i64, ell: u64) -> Option<u64> {
    let l = ell as i64;
    let inv2 = (l + 1) / 2; // 2⁻¹ mod odd ℓ
    let lam = ((t % l + l) % l * inv2 % l) as u64;
    ((lam * lam) % ell == q % ell && lam != 0).then_some(lam)
}

// ── the ℓ-torsion and its kernels ───────────────────────────────────

/// What the `ℓ`-Sylow of `E(F_{q^m})` turned out to be, each outcome
/// **proved** rather than guessed.
pub enum Torsion {
    /// Two independent points of order `ℓ`: all `ℓ + 1` lines are here, the
    /// vertex has descending `ℓ`-edges.
    Rank2(Pt, Pt),
    /// The Sylow is cyclic — proved by exhibiting an element whose order is
    /// the whole Sylow's.  Only one line of `E[ℓ]` is rational: the
    /// `ℓ`-floor.
    Rank1,
    /// `ℓ ∤ #E(F_{q^m})`.
    None,
    /// Neither proof completed within the sampling budget.  Reported, never
    /// folded into either count.
    Undecided,
}

impl<'a> ExtCurve<'a> {
    /// `[s]P` by affine double-and-add.
    pub fn mul_pt(&self, p: &Pt, s: &BigUint) -> Pt {
        if s.is_zero() || *p == Pt::O {
            return Pt::O;
        }
        let mut acc = Pt::O;
        for i in (0..s.bits()).rev() {
            acc = self.double(&acc);
            if s.bit(i) {
                acc = self.add(&acc, p);
            }
        }
        acc
    }

    fn random_point(&self, rng: &mut StdRng) -> Pt {
        let x = self.random_x(rng);
        let y = self.lift_x(&x, rng).expect("random_x carries a point");
        Pt::A(x, y)
    }

    /// Exponent `k` of the order `ℓ^k` of a point known to lie in the
    /// `ℓ`-Sylow, by full arithmetic.
    fn ell_exponent(&self, p: &Pt, l: &BigUint) -> u32 {
        let mut k = 0;
        let mut cur = p.clone();
        while cur != Pt::O {
            cur = self.mul_pt(&cur, l);
            k += 1;
        }
        k
    }
}

/// Search `E(F_{q^m})` for an `ℓ`-torsion basis, with both outcomes proved.
///
/// **Rank 1** is proved by an element of order `ℓ^v`, `v = v_ℓ(#E(F_{q^m}))`:
/// then the Sylow is cyclic and has one line of order `ℓ`.  That needs only
/// the `x`-only ladder, so it is tried first and costs one sample in the
/// common case.
///
/// **Rank 2** is proved by two points of order `ℓ` on different lines.  They
/// cannot be found by clearing the cofactor and multiplying down by `ℓ`:
/// when the Sylow is `Z/ℓ^a × Z/ℓ^b` with `b < a`, `ℓ^{a−1}·S` is a single
/// line and almost every sample lands on it.  An earlier revision did exactly
/// that and undercounted the rank-2 vertices at `ℓ = 3`, `n = 16` (31 and 32
/// against the 33 the volcano predicts).  Instead, a second sample is reduced
/// against the first — `s₂ ← s₂ − [c·ℓ^{a−k}]s₁` whenever `s₂`'s bottom
/// element lies on `s₁`'s line — which strictly lowers its order until it
/// either vanishes (resample) or exposes a second line.
pub fn torsion_basis(
    curve: &ExtCurve<'_>,
    group_order: &BigUint,
    ell: u64,
    attempts: usize,
    rng: &mut StdRng,
) -> Torsion {
    let l = BigUint::from(ell);
    let mut u = group_order.clone();
    let mut v = 0u32;
    while (&u % &l).is_zero() {
        u /= &l;
        v += 1;
    }
    if v == 0 {
        return Torsion::None;
    }
    let d = ((ell - 1) / 2) as usize;
    let attempts = attempts.max(3);

    // ── rank 1: an element of full Sylow order, x-only ──
    for _ in 0..attempts {
        let x = curve.random_x(rng);
        let Some(mut s) = curve.ladder_x(&x, &u) else {
            continue;
        };
        let mut k = 1u32;
        while let Some(next) = curve.ladder_x(&s, &l) {
            s = next;
            k += 1;
        }
        if k == v {
            return Torsion::Rank1;
        }
    }

    // ── rank 2: reduce samples against a maximal-order one ──
    let sylow_sample = |rng: &mut StdRng| -> Pt { curve.mul_pt(&curve.random_point(rng), &u) };
    let pow = |e: u32| -> BigUint { num_traits::pow(l.clone(), e as usize) };

    let mut s1 = sylow_sample(rng);
    let mut a = curve.ell_exponent(&s1, &l);
    for _ in 0..attempts {
        if a == v {
            return Torsion::Rank1;
        }
        if a > 0 {
            break;
        }
        s1 = sylow_sample(rng);
        a = curve.ell_exponent(&s1, &l);
    }
    if a == 0 {
        return Torsion::Undecided;
    }
    let mut t1 = curve.mul_pt(&s1, &pow(a - 1));
    let mut line: BTreeSet<El> = match &t1 {
        Pt::A(x, _) => subgroup_abscissae(curve, x, d).into_iter().collect(),
        Pt::O => unreachable!("order ℓ^a with a ≥ 1"),
    };

    for _ in 0..attempts * 4 {
        let mut s2 = sylow_sample(rng);
        loop {
            if s2 == Pt::O {
                break; // s2 ∈ ⟨s1⟩: resample
            }
            let k = curve.ell_exponent(&s2, &l);
            if k == v {
                return Torsion::Rank1;
            }
            if k > a {
                // a larger-order element: make it the reference
                s1 = s2;
                a = k;
                t1 = curve.mul_pt(&s1, &pow(a - 1));
                line = match &t1 {
                    Pt::A(x, _) => subgroup_abscissae(curve, x, d).into_iter().collect(),
                    Pt::O => unreachable!(),
                };
                break;
            }
            let t = curve.mul_pt(&s2, &pow(k - 1));
            let Pt::A(tx, _) = &t else {
                unreachable!("order ℓ^k, k ≥ 1")
            };
            if !line.contains(tx) {
                return Torsion::Rank2(t1, t);
            }
            // t = [c]t1 for some c in 1..ℓ; subtract its lift from s2.
            let mut c = 0u64;
            let mut acc = Pt::O;
            for j in 1..ell {
                acc = curve.add(&acc, &t1);
                if acc == t {
                    c = j;
                    break;
                }
            }
            assert!(c != 0, "t lies on t1's line, so it is a multiple of t1");
            let lift = curve.mul_pt(&s1, &(BigUint::from(c) * pow(a - k)));
            s2 = curve.add(&s2, &curve.neg(&lift));
        }
    }
    Torsion::Undecided
}

/// `x(jG)` for `j = 1..d`, by `x`-only differential addition.
pub fn subgroup_abscissae(curve: &ExtCurve<'_>, xg: &El, d: usize) -> Vec<El> {
    let mut xs = Vec::with_capacity(d);
    xs.push(xg.clone());
    if d == 1 {
        return xs;
    }
    xs.push(curve.x_double(xg).expect("order > 2, so x ≠ 0"));
    for j in 2..d {
        // x((j+1)G) = x((j−1)G) + x(jG)x(G)/(x(jG)+x(G))²
        let next = curve
            .x_diff_add(&xs[j - 1], xg, &xs[j - 2])
            .expect("jG ≠ ±G for j < (ℓ−1)/2");
        xs.push(next);
    }
    xs
}

/// One kernel, ready for Vélu.
pub struct Kernel {
    /// `x(jG)`, `j = 1..(ℓ−1)/2`, in `F_{q^m}`.
    pub abscissae: Vec<El>,
    /// `t = Σ x(jG)`, which must lie in `F_q`.
    pub t: El,
    pub t_is_rational: bool,
}

/// The `ℓ + 1` kernel generators of a basis: `P₁` and `P₂ + kP₁`,
/// `k = 0..ℓ−1`.  Only generators are held; each kernel's abscissae are
/// built on demand by [`kernel_of`], because all of them at once is
/// `(ℓ+1)(ℓ−1)/2` extension elements — 48 GB at `(n, ℓ) = (31, 7193)`.
pub fn kernel_generators(curve: &ExtCurve<'_>, p1: &Pt, p2: &Pt, ell: u64) -> Vec<El> {
    let mut gens = Vec::with_capacity(ell as usize + 1);
    let x_of = |g: &Pt| match g {
        Pt::A(x, _) => x.clone(),
        Pt::O => unreachable!("generators have order ℓ"),
    };
    gens.push(x_of(p1));
    let mut cur = p2.clone();
    for _ in 0..ell {
        gens.push(x_of(&cur));
        cur = curve.add(&cur, p1);
    }
    gens
}

/// The kernel generated by the point with abscissa `xg`.
pub fn kernel_of(curve: &ExtCurve<'_>, xg: &El, ell: u64) -> Kernel {
    let k = curve.k;
    let abscissae = subgroup_abscissae(curve, xg, ((ell - 1) / 2) as usize);
    let t = abscissae.iter().fold(k.zero(), |acc, x| k.add(&acc, x));
    let t_is_rational = k.is_base(&t);
    Kernel {
        abscissae,
        t,
        t_is_rational,
    }
}

/// All `ℓ + 1` kernels from a basis.  Materialises every abscissa; for
/// small `ℓ` only — the walk streams [`kernel_of`] instead.
pub fn all_kernels(curve: &ExtCurve<'_>, p1: &Pt, p2: &Pt, ell: u64) -> Vec<Kernel> {
    kernel_generators(curve, p1, p2, ell)
        .iter()
        .map(|xg| kernel_of(curve, xg, ell))
        .collect()
}

/// Vélu's abscissa map at a base-field point, `X = x + A + A²`,
/// `A = (d mod 2) + x·Σ 1/(x + x_j)`.  Returns the value and whether it lies
/// in `F_q`, which it must for a rational isogeny.
pub fn velu_x_from_abscissae(k: &Ext, abscissae: &[El], x0: u64) -> Option<(u64, bool)> {
    let x = k.base(x0);
    // batch inversion of (x + x_j)
    let dens: Vec<El> = abscissae.iter().map(|xj| k.add(&x, xj)).collect();
    if dens.iter().any(|v| k.is_zero(v)) {
        return None; // x0 is a kernel abscissa: maps to O
    }
    let mut prefix = Vec::with_capacity(dens.len());
    let mut acc = k.one();
    for v in &dens {
        prefix.push(acc.clone());
        acc = k.mul(&acc, v);
    }
    let mut inv_all = k.inv(&acc)?;
    let mut sum = k.zero();
    for i in (0..dens.len()).rev() {
        let inv_i = k.mul(&inv_all, &prefix[i]);
        sum = k.add(&sum, &inv_i);
        inv_all = k.mul(&inv_all, &dens[i]);
    }
    let mut a = k.mul(&x, &sum);
    if abscissae.len() % 2 == 1 {
        a = k.add(&a, &k.one());
    }
    let xx = k.add(&k.add(&x, &a), &k.sqr(&a));
    Some((xx[0], k.is_base(&xx)))
}

// ── the walk ────────────────────────────────────────────────────────

/// Kernels built in parallel per batch; also the checkpoint granularity.
const KERNEL_CHUNK: usize = 64;
/// Degrees large enough that per-batch progress is worth printing.
const LARGE_ELL: u64 = 1000;

/// What one kernel contributes, reduced to base-field data so it can be
/// checkpointed: the extension-field abscissae are dropped once used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum KernelOutcome {
    /// `t ∉ F_q` — a failure of the construction, counted.
    IrrationalT,
    /// `x(P)` or `x(Q)` is a kernel abscissa; no image to carry.
    Singular,
    Edge {
        a6p: u64,
        xpp: u64,
        xqq: u64,
        rp: bool,
        rq: bool,
    },
}

impl KernelOutcome {
    fn of(ext: &Ext, ker: &Kernel, a6: u64, xp: u64, xq: u64) -> Self {
        if !ker.t_is_rational {
            return KernelOutcome::IrrationalT;
        }
        let tv = ker.t[0];
        let a6p = a6 ^ tv ^ ext.gf.mul(tv, tv);
        let img_p = velu_x_from_abscissae(ext, &ker.abscissae, xp);
        let img_q = velu_x_from_abscissae(ext, &ker.abscissae, xq);
        match (img_p, img_q) {
            (Some((xpp, rp)), Some((xqq, rq))) => KernelOutcome::Edge {
                a6p,
                xpp,
                xqq,
                rp,
                rq,
            },
            _ => KernelOutcome::Singular,
        }
    }

    fn encode(&self) -> String {
        match self {
            KernelOutcome::IrrationalT => "I".into(),
            KernelOutcome::Singular => "S".into(),
            KernelOutcome::Edge {
                a6p,
                xpp,
                xqq,
                rp,
                rq,
            } => {
                format!("E {a6p} {xpp} {xqq} {} {}", *rp as u8, *rq as u8)
            }
        }
    }

    fn decode(f: &[&str]) -> Option<Self> {
        match *f.first()? {
            "I" => Some(KernelOutcome::IrrationalT),
            "S" => Some(KernelOutcome::Singular),
            "E" => Some(KernelOutcome::Edge {
                a6p: f.get(1)?.parse().ok()?,
                xpp: f.get(2)?.parse().ok()?,
                xqq: f.get(3)?.parse().ok()?,
                rp: *f.get(4)? == "1",
                rq: *f.get(5)? == "1",
            }),
            _ => None,
        }
    }
}

/// Decide `E[ℓ]`'s rank at each vertex in parallel and record it; a rank-2
/// vertex also gets its kernel generators stored.  One RNG per
/// `(vertex, ℓ)`, so a vertex's basis — and the kernel indices the
/// checkpoint is keyed by — never depends on what a resumed walk skipped.
#[allow(clippy::too_many_arguments)]
fn classify(
    ckpt: &mut Checkpoint,
    a2: u8,
    vertices: &[u64],
    ell: u64,
    ext: &Ext,
    order_m: &BigUint,
    attempts: usize,
    seed: u64,
) {
    let results: Vec<(u64, u8, Option<Vec<El>>)> = vertices
        .par_iter()
        .map(|&a6| {
            let curve = ExtCurve::new(ext, a2 as u64, a6);
            let mut rng = StdRng::seed_from_u64(
                seed ^ a6.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ ell.rotate_left(32),
            );
            match torsion_basis(&curve, order_m, ell, attempts, &mut rng) {
                Torsion::Rank2(p1, p2) => (a6, 2, Some(kernel_generators(&curve, &p1, &p2, ell))),
                Torsion::Rank1 => (a6, 1, None),
                Torsion::None => (a6, 0, None),
                Torsion::Undecided => (a6, 3, None),
            }
        })
        .collect();
    for (a6, r, gens) in results {
        if let Some(g) = gens {
            ckpt.store_generators(a6, ell, &g);
        }
        ckpt.record_rank(a6, ell, r);
    }
}

/// Append-only record of finished kernels, keyed by `(a₆, ℓ, index)`.
/// The walk is deterministic in its seed — same BFS order, same bases, same
/// kernel indices — so a resumed walk reads back exactly what it would have
/// recomputed.  A torn last line is ignored.
struct Checkpoint {
    done: BTreeMap<(u64, u64, usize), KernelOutcome>,
    /// `(a₆, ℓ) ↦` rank of `E[ℓ]`: 0 none, 1, 2, 3 undecided.
    rank: BTreeMap<(u64, u64), u8>,
    /// Generators held in memory when there is no checkpoint directory.
    gens: BTreeMap<(u64, u64), Vec<El>>,
    file: Option<std::fs::File>,
    path: Option<std::path::PathBuf>,
}

impl Checkpoint {
    fn open(path: Option<&Path>) -> Self {
        let mut done = BTreeMap::new();
        let mut rank = BTreeMap::new();
        let Some(path) = path else {
            return Checkpoint {
                done,
                rank,
                gens: BTreeMap::new(),
                file: None,
                path: None,
            };
        };
        if let Ok(text) = std::fs::read_to_string(path) {
            for line in text.lines() {
                let f: Vec<&str> = line.split_whitespace().collect();
                if f.first() == Some(&"T") && f.len() == 4 {
                    if let (Ok(a6), Ok(ell), Ok(r)) = (f[1].parse(), f[2].parse(), f[3].parse()) {
                        rank.insert((a6, ell), r);
                    }
                    continue;
                }
                if f.len() < 4 {
                    continue;
                }
                let (Ok(a6), Ok(ell), Ok(i)) = (f[0].parse(), f[1].parse(), f[2].parse()) else {
                    continue;
                };
                if let Some(o) = KernelOutcome::decode(&f[3..]) {
                    done.insert((a6, ell, i), o);
                }
            }
        }
        let file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .ok();
        Checkpoint {
            done,
            rank,
            gens: BTreeMap::new(),
            file,
            path: Some(path.to_path_buf()),
        }
    }

    /// Kernel generators of a rank-2 vertex are the expensive setup — a
    /// torsion basis and `ℓ` additions in `F_{q^m}`, minutes at `ℓ = 7193` —
    /// so they are stored beside the checkpoint and a resumed walk starts
    /// computing kernels at once.  Written to a temporary name and renamed,
    /// so a file that exists is complete.
    fn generators_path(&self, a6: u64, ell: u64) -> Option<std::path::PathBuf> {
        let p = self.path.as_ref()?;
        Some(p.with_extension(format!("gens_{a6}_{ell}")))
    }

    fn record_rank(&mut self, a6: u64, ell: u64, r: u8) {
        // A rank-2 record is only trusted once its generators are on disk,
        // so it is written after them.
        if let Some(f) = self.file.as_mut() {
            let _ = writeln!(f, "T {a6} {ell} {r}");
        }
        self.rank.insert((a6, ell), r);
    }

    fn load_generators(&self, a6: u64, ell: u64, m: usize) -> Option<Vec<El>> {
        if let Some(g) = self.gens.get(&(a6, ell)) {
            return Some(g.clone());
        }
        let bytes = std::fs::read(self.generators_path(a6, ell)?).ok()?;
        let words: Vec<u64> = bytes
            .chunks_exact(8)
            .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        if words.len() != (ell as usize + 1) * m {
            return None;
        }
        Some(words.chunks(m).map(|c| c.to_vec()).collect())
    }

    fn store_generators(&mut self, a6: u64, ell: u64, gens: &[El]) {
        let Some(p) = self.generators_path(a6, ell) else {
            self.gens.insert((a6, ell), gens.to_vec());
            return;
        };
        let bytes: Vec<u8> = gens
            .iter()
            .flat_map(|g| g.iter().flat_map(|w| w.to_le_bytes()))
            .collect();
        let tmp = p.with_extension("tmp");
        if std::fs::write(&tmp, bytes).is_ok() {
            let _ = std::fs::rename(tmp, p);
        }
    }

    fn record(&mut self, a6: u64, ell: u64, i: usize, o: KernelOutcome) {
        if let Some(f) = self.file.as_mut() {
            let _ = writeln!(f, "{a6} {ell} {i} {}", o.encode());
        }
        self.done.insert((a6, ell, i), o);
    }
}

/// One edge the walk took, with everything checked on it.
#[derive(Clone, Debug)]
pub struct Edge {
    pub ell: u64,
    pub from: u64,
    pub to: u64,
    pub t_rational: bool,
    pub x_rational: bool,
    pub order_preserved: bool,
    pub image_order_ok: bool,
    pub transported: bool,
}

/// The result of walking a class.
#[derive(Clone, Debug, Default)]
pub struct WalkReport {
    pub n: u32,
    pub a2: u8,
    pub ells: Vec<u64>,
    /// `ℓ ↦ m = ord_ℓ(λ)`, the degree of the field `E[ℓ]` was built in.
    pub ext_degree: BTreeMap<u64, usize>,
    /// `a₆` of every vertex reached, with its depth.
    pub reached: BTreeMap<u64, u32>,
    pub edges: Vec<Edge>,
    /// Vertices where `E[ℓ]` had rank 2 (descending edges exist) and rank 1
    /// (the `ℓ`-floor), per `ℓ`.
    pub rank2: BTreeMap<u64, usize>,
    pub rank1: BTreeMap<u64, usize>,
    /// Vertices where neither rank could be proved within the budget.
    pub undecided: BTreeMap<u64, usize>,
    /// Kernels whose `t` or transported `X` left `F_q` — each one is a
    /// failure of the construction and is counted, never dropped.
    pub irrational_kernels: usize,
    pub irrational_images: usize,
    pub seconds: f64,
}

impl WalkReport {
    pub fn transport_failures(&self) -> usize {
        self.edges
            .iter()
            .filter(|e| !(e.order_preserved && e.image_order_ok && e.transported))
            .count()
    }
}

/// **Walk the isogeny class of `K_a` by descending edges**, carrying the DLP
/// instance `(P, Q = [d]P)` along every edge.
///
/// `instance` is `(x(P), x(Q), d, r)` on the Koblitz curve (`a₆ = 1`).  Only
/// abscissae are carried: `φ(±P)` share one, and the checks below are all
/// sign-insensitive.
#[allow(clippy::too_many_arguments)]
pub fn walk_class(
    n: u32,
    irr: &IrreduciblePoly,
    a2: u8,
    ells: &[u64],
    instance: (u64, u64, u64, u64),
    attempts: usize,
    max_vertices: usize,
    seed: u64,
    checkpoint: Option<&Path>,
    mut progress: impl FnMut(&WalkReport),
) -> WalkReport {
    let started = Instant::now();
    let mut ckpt = Checkpoint::open(checkpoint);
    let q = 1u64 << n;
    let base = Curve::new(n, irr, a2, 1).expect("the Koblitz curve");
    let t = (q + 1) as i64 - base.order as i64;
    let (xp0, xq0, dlog, r) = instance;

    let mut report = WalkReport {
        n,
        a2,
        ells: ells.to_vec(),
        ..Default::default()
    };

    // One extension per ℓ, shared by the whole class: t, hence λ and m, is a
    // class invariant.
    let mut fields: Vec<(u64, Ext, BigUint)> = Vec::new();
    for &ell in ells {
        let Some(lam) = scalar_eigenvalue(q, t, ell) else {
            continue;
        };
        let m = mult_order(lam, ell) as usize;
        let ext = Ext::new(irr, m, seed ^ ell);
        let order_m = order_over_extension(q, t, m);
        report.ext_degree.insert(ell, m);
        fields.push((ell, ext, order_m));
    }

    let mut queue: VecDeque<(u64, u64, u64, u32)> = VecDeque::new();
    report.reached.insert(1, 0);
    queue.push_back((1, xp0, xq0, 0));

    while let Some((a6, xp, xq, depth)) = queue.pop_front() {
        if report.reached.len() >= max_vertices {
            break;
        }
        // Classify every queued vertex at once, in parallel: at a large
        // class nearly all of them are the floor, where proving rank 1 is the
        // whole cost and no edge depends on another vertex's result.
        if !queue.is_empty() {
            for (ell, ext, order_m) in &fields {
                let pending: Vec<u64> = std::iter::once(a6)
                    .chain(queue.iter().map(|v| v.0))
                    .filter(|v| {
                        !ckpt.rank.contains_key(&(*v, *ell))
                            && ckpt.load_generators(*v, *ell, ext.m).is_none()
                    })
                    .collect();
                for (c, chunk) in pending.chunks(KERNEL_CHUNK).enumerate() {
                    classify(&mut ckpt, a2, chunk, *ell, ext, order_m, attempts, seed);
                    if *ell >= LARGE_ELL {
                        eprintln!(
                            "     … ℓ={ell}: classified {} of {} queued vertices, {:.0}s",
                            (c * KERNEL_CHUNK + chunk.len()),
                            pending.len(),
                            started.elapsed().as_secs_f64()
                        );
                    }
                }
            }
        }
        for (ell, ext, order_m) in &fields {
            let curve = ExtCurve::new(ext, a2 as u64, a6);
            if !ckpt.rank.contains_key(&(a6, *ell))
                && ckpt.load_generators(a6, *ell, ext.m).is_some()
            {
                ckpt.record_rank(a6, *ell, 2);
            }
            if !ckpt.rank.contains_key(&(a6, *ell)) {
                classify(&mut ckpt, a2, &[a6], *ell, ext, order_m, attempts, seed);
            }
            let gens = match ckpt.rank[&(a6, *ell)] {
                2 => Ok(ckpt
                    .load_generators(a6, *ell, ext.m)
                    .expect("a rank-2 vertex has stored generators")),
                1 => Err(Torsion::Rank1),
                0 => Err(Torsion::None),
                _ => Err(Torsion::Undecided),
            };
            match gens {
                Ok(gens) => {
                    *report.rank2.entry(*ell).or_insert(0) += 1;
                    let mut outcomes = Vec::with_capacity(gens.len());
                    for (c, chunk) in gens.chunks(KERNEL_CHUNK).enumerate() {
                        let first = c * KERNEL_CHUNK;
                        let todo: Vec<usize> = (first..first + chunk.len())
                            .filter(|i| !ckpt.done.contains_key(&(a6, *ell, *i)))
                            .collect();
                        let fresh: Vec<(usize, KernelOutcome)> = todo
                            .par_iter()
                            .map(|&i| {
                                let ker = kernel_of(&curve, &gens[i], *ell);
                                (i, KernelOutcome::of(ext, &ker, a6, xp, xq))
                            })
                            .collect();
                        for (i, o) in fresh {
                            ckpt.record(a6, *ell, i, o);
                        }
                        for i in first..first + chunk.len() {
                            outcomes.push(ckpt.done[&(a6, *ell, i)]);
                        }
                        if *ell >= LARGE_ELL {
                            eprintln!(
                                "     … a6={a6}, ℓ={ell}: {} of {} kernels, {:.0}s",
                                outcomes.len(),
                                gens.len(),
                                started.elapsed().as_secs_f64()
                            );
                        }
                    }
                    for o in outcomes {
                        let (a6p, xpp, xqq, rp, rq) = match o {
                            KernelOutcome::IrrationalT => {
                                report.irrational_kernels += 1;
                                continue;
                            }
                            KernelOutcome::Singular => continue,
                            KernelOutcome::Edge {
                                a6p,
                                xpp,
                                xqq,
                                rp,
                                rq,
                            } => (a6p, xpp, xqq, rp, rq),
                        };
                        if a6p == 0 || report.reached.contains_key(&a6p) {
                            continue;
                        }
                        if !(rp && rq) {
                            report.irrational_images += 1;
                        }
                        // `#E' = #E` is certified without counting when r
                        // exceeds the Hasse width 4√q: a point of order r on
                        // E' puts r | #E', and #E is the only multiple of r
                        // in the interval.  Otherwise count, as before.
                        let certify_by_r = (r as f64) > 4.0 * (q as f64).sqrt() + 2.0;
                        let codomain = if certify_by_r {
                            Curve::with_order(n, irr, a2, a6p, base.order)
                        } else {
                            Curve::new(n, irr, a2, a6p)
                        }
                        .expect("codomain");
                        let (mut image_order_ok, mut transported) = (false, false);
                        if let Some(&pp) = codomain.points_with_x(xpp).first() {
                            let g = codomain.group();
                            let mut ops = GroupOps::default();
                            image_order_ok = !pp.infinity && g.mul(&mut ops, pp, r).infinity;
                            let dp: FastPoint = g.mul(&mut ops, pp, dlog % r);
                            transported = !dp.infinity && dp.x == xqq;
                        }
                        let order_preserved = if certify_by_r {
                            image_order_ok
                        } else {
                            codomain.order == base.order
                        };
                        report.edges.push(Edge {
                            ell: *ell,
                            from: a6,
                            to: a6p,
                            t_rational: true,
                            x_rational: rp && rq,
                            order_preserved,
                            image_order_ok,
                            transported,
                        });
                        report.reached.insert(a6p, depth + 1);
                        queue.push_back((a6p, xpp, xqq, depth + 1));
                    }
                }
                Err(Torsion::Rank1) => {
                    *report.rank1.entry(*ell).or_insert(0) += 1;
                }
                Err(Torsion::Undecided) => {
                    *report.undecided.entry(*ell).or_insert(0) += 1;
                }
                Err(_) => {}
            }
        }
        report.seconds = started.elapsed().as_secs_f64();
        progress(&report);
    }
    report.seconds = started.elapsed().as_secs_f64();
    report
}

/// The Koblitz curve's DLP instance with the class-wide secret the cost
/// sweep used, so the walk carries *that* instance: `(x(P), x(Q), d, r)`.
pub fn koblitz_instance(
    n: u32,
    irr: &IrreduciblePoly,
    a2: u8,
    seed: u64,
) -> Option<(u64, u64, u64, u64)> {
    let curve = Curve::new(n, irr, a2, 1)?;
    let mut r = curve.order;
    let mut f = 2u64;
    let mut largest = 1u64;
    while f * f <= r {
        while r % f == 0 {
            largest = largest.max(f);
            r /= f;
        }
        f += 1;
    }
    if r > 1 {
        largest = largest.max(r);
    }
    let r = largest;
    let cof = curve.order / r;
    let g = curve.group();
    let mut ops = GroupOps::default();
    for x in 1..(1u64 << n) {
        let Some(&p) = curve.points_with_x(x).first() else {
            continue;
        };
        let gen = g.mul(&mut ops, p, cof);
        if gen.infinity || !g.mul(&mut ops, gen, r).infinity {
            continue;
        }
        let d = crate::cryptanalysis::koblitz_isogeny_cost::class_secret(seed, n, r);
        let q = g.mul(&mut ops, gen, d);
        return Some((gen.x, q.x, d, r));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::binary_isogeny::find_roots_in_f2m;
    use crate::cryptanalysis::binary_velu::{division_polynomial, elt, to_u64};
    use crate::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
    use num_traits::{Signed, ToPrimitive};

    fn field(n: u32) -> IrreduciblePoly {
        find_irreducible_sparse(n).expect("irreducible")
    }

    /// The `x`-only formulas against the repository's own point arithmetic,
    /// over the base field (`m = 1`), where both are available.
    #[test]
    fn x_only_formulas_match_point_arithmetic() {
        let n = 10;
        let irr = field(n);
        let ext = Ext::new(&irr, 1, 1);
        for a2 in [0u8, 1] {
            let base = Curve::new(n, &irr, a2, 7).expect("curve");
            let ec = ExtCurve::new(&ext, a2 as u64, 7);
            let g = base.group();
            let mut ops = GroupOps::default();
            let mut checked = 0;
            for x in 1..(1u64 << n) {
                let Some(&p) = base.points_with_x(x).first() else {
                    continue;
                };
                for s in [2u64, 3, 5, 7, 12, 97] {
                    let want = g.mul(&mut ops, p, s);
                    let got = ec.ladder_x(&ext.base(x), &BigUint::from(s));
                    match (want.infinity, got) {
                        (true, None) => {}
                        (false, Some(v)) => assert_eq!(v[0], want.x, "ladder [{s}]·x={x}"),
                        (w, g2) => panic!("ladder disagrees at x={x} s={s}: O={w} got={g2:?}"),
                    }
                }
                // differential chain against repeated addition
                let xs = subgroup_abscissae(&ec, &ext.base(x), 6);
                let mut acc = p;
                for (j, xj) in xs.iter().enumerate() {
                    if acc.infinity {
                        break;
                    }
                    assert_eq!(xj[0], acc.x, "x(({})P) at x={x}", j + 1);
                    acc = g.add(&mut ops, acc, p);
                }
                checked += 1;
                if checked > 40 {
                    break;
                }
            }
            assert!(checked > 0);
        }
    }

    #[test]
    fn extension_arithmetic_is_a_field() {
        let n = 7;
        let irr = field(n);
        let ext = Ext::new(&irr, 5, 3);
        let mut rng = StdRng::seed_from_u64(9);
        for _ in 0..50 {
            let a = ext.random(&mut rng);
            if ext.is_zero(&a) {
                continue;
            }
            let ai = ext.inv(&a).expect("nonzero is a unit in a field");
            assert_eq!(ext.mul(&a, &ai), ext.one());
            let b = ext.random(&mut rng);
            assert_eq!(
                ext.sqr(&ext.add(&a, &b)),
                ext.add(&ext.sqr(&a), &ext.sqr(&b))
            );
        }
    }

    #[test]
    fn artin_schreier_solves_in_odd_and_even_degree() {
        for (n, m) in [(7u32, 3usize), (6, 4)] {
            let irr = field(n);
            let ext = Ext::new(&irr, m, 5);
            let mut rng = StdRng::seed_from_u64(11);
            let mut solved = 0;
            for _ in 0..40 {
                let c = ext.random(&mut rng);
                if let Some(w) = ext.solve_as(&c, &mut rng) {
                    assert_eq!(ext.add(&ext.sqr(&w), &w), c, "w² + w = c");
                    solved += 1;
                }
            }
            assert!(solved > 0, "about half of all c have trace 0");
        }
    }

    #[test]
    fn affine_addition_agrees_with_the_ladder() {
        let n = 9;
        let irr = field(n);
        let ext = Ext::new(&irr, 2, 4);
        let ec = ExtCurve::new(&ext, 1, 5);
        let mut rng = StdRng::seed_from_u64(3);
        for _ in 0..10 {
            let x = ec.random_x(&mut rng);
            let y = ec.lift_x(&x, &mut rng).expect("random_x carries a point");
            let p = Pt::A(x.clone(), y);
            let mut acc = p.clone();
            for s in 2..8u64 {
                acc = ec.add(&acc, &p);
                let via_ladder = ec.ladder_x(&x, &BigUint::from(s));
                match (&acc, via_ladder) {
                    (Pt::O, None) => {}
                    (Pt::A(ax, _), Some(lx)) => assert_eq!(ax, &lx, "[{s}]P"),
                    _ => panic!("affine and ladder disagree at s = {s}"),
                }
            }
        }
    }

    /// The load-bearing test: at `ℓ = 3`, the kernels this module builds
    /// from points must be exactly the kernels `ψ₃`'s roots name — and those
    /// were checked against the tabulated `Φ₃` in `binary_velu`.  So this
    /// chains the new construction to an independent oracle.
    #[test]
    fn torsion_kernels_equal_psi3_kernels_at_the_crater() {
        let n = 16;
        let irr = field(n);
        let q = 1u64 << n;
        for a2 in [0u8, 1] {
            let base = Curve::new(n, &irr, a2, 1).expect("K_a");
            let t = (q + 1) as i64 - base.order as i64;
            let lam = scalar_eigenvalue(q, t, 3).expect("3 | c at n = 16");
            let m = mult_order(lam, 3) as usize;
            let ext = Ext::new(&irr, m, 17);
            let om = order_over_extension(q, t, m);
            let ec = ExtCurve::new(&ext, a2 as u64, 1);
            let mut rng = StdRng::seed_from_u64(23);
            let Torsion::Rank2(p1, p2) = torsion_basis(&ec, &om, 3, 8, &mut rng) else {
                panic!("K_a is the crater: E[3] must have rank 2 over F_(q^{m})");
            };
            let kers = all_kernels(&ec, &p1, &p2, 3);
            assert_eq!(kers.len(), 4);
            let from_points: BTreeSet<u64> = kers
                .iter()
                .map(|k| {
                    assert!(k.t_is_rational, "every kernel is rational at the crater");
                    k.t[0]
                })
                .collect();
            let psi3 = division_polynomial(&elt(1, n), 3, n, &irr);
            let from_psi: BTreeSet<u64> = find_roots_in_f2m(&psi3, n, &irr)
                .iter()
                .map(to_u64)
                .collect();
            assert_eq!(from_points, from_psi, "a₂ = {a2}");
        }
    }

    #[test]
    fn group_order_over_the_base_field_matches_point_counting() {
        let n = 12;
        let irr = field(n);
        let q = 1u64 << n;
        let c = Curve::new(n, &irr, 1, 1).expect("curve");
        let t = (q + 1) as i64 - c.order as i64;
        assert_eq!(order_over_extension(q, t, 1), BigUint::from(c.order));
        let _ = (BigInt::from(1).abs(), BigUint::from(1u8).to_u64());
    }
}
