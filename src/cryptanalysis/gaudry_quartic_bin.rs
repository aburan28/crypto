//! **Gaudry's `k = 4` decomposition in characteristic 2, with FGHR's
//! 2-torsion symmetry: `E/F_{q⁴}`, `q = 2^m`, the factor base on a line
//! in a coordinate where the 2-torsion translation is `t ↦ t + 1`.**
//!
//! The characteristic-2 counterpart of [`super::gaudry_quartic_d4`].  An
//! ordinary curve `E: y² + xy = x³ + ax² + b` over `F_{q⁴}` always has the
//! rational 2-torsion point `T = (0, √b)`, and `x(P + T) = √b / x(P)`.  With
//! `β = b^{1/4}` and
//!
//! ```text
//!     t = β / (x + β),        x = β (t + 1) / t,
//! ```
//!
//! the translation is `t(P + T) = t(P) + 1` for every `a, b` (no subfield
//! condition), so the base `{P : t(P) ∈ F_q}` is `T`-stable.  Since
//! `2T = 0`, translating an even number of the four summands keeps their
//! sum, and `H(t₁, …, t₄; t_R) = S₅(x(t₁), …, x(t_R)) · Π t_i⁸ · t_R⁸` is
//! invariant under permutations of the `t_i` and even numbers of
//! translations `t_i ↦ t_i + 1`.  The translations' invariants are
//! `w_i = t_i² + t_i` and `π = t₁ + t₂ + t₃ + t₄` with `π² + π = Σ w_i`;
//! with the symmetric group, the invariant ring is polynomial in
//! `s₂, s₃, s₄ = e_k(w)` and `π`, where `s₁ = π² + π`.  Each `s_k` has
//! degree 2 in every `t_i` and `π` degree 1: the weights `(2, 2, 2, 1)` of
//! the odd-characteristic system, the same 105 invariant monomials, the same
//! Hilbert series and the same weighted Macaulay degree 26.
//!
//! Everything else is over `F_q = F_{2^m}` (table arithmetic, every
//! multiplication counted): the interpolation of `H` on the 105 monomials
//! (checked on fresh points), its Weil restriction to four equations over
//! `F_q`, their weighted Macaulay solve with the eigen step along the chains
//! of `π` ([`solve_bin`], the binary form of `gaudry_quartic`'s
//! `solve_along_chains`), the split of each solution back to four points,
//! and an independent meet-in-the-middle enumeration over the same base.

use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet, HashMap};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;

use super::ic_boundary::{CountedGroup, FactorBase, GroupOps, OracleCounters, Relation};
use super::ic_framework::stages::{
    DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, SolverTotals, SystemShape,
};

// ── F_{2^m} ───────────────────────────────────────────────────────────────

/// An element of `F_q`, `q = 2^m ≤ 2^20`, as its bits.  It was a `u16`
/// while `m ≤ 15`; the values, and so every count, are the same.
pub type Fe = u32;

/// Primitive polynomials of degree `m` (the `x^m` bit included).
const PRIMITIVE: [u32; 21] = [
    0, 0, 0x7, 0xB, 0x13, 0x25, 0x43, 0x89, 0x11D, 0x211, 0x409, 0x805, 0x1053, 0x201B, 0x4443,
    0x8003, 0x1100B, 0x20009, 0x40081, 0x80027, 0x100009,
];

/// The logarithm tables of `F_{2^m}`: shareable across threads, no
/// counter (the solver's loops count for themselves).
pub struct GfTab {
    pub q: u32,
    exp: Vec<Fe>,
    log: Vec<u32>,
}

impl GfTab {
    #[inline(always)]
    pub fn mul(&self, a: Fe, b: Fe) -> Fe {
        if a == 0 || b == 0 {
            0
        } else {
            self.exp[(self.log[a as usize] + self.log[b as usize]) as usize]
        }
    }
    #[inline(always)]
    pub fn inv(&self, a: Fe) -> Fe {
        assert!(a != 0, "inverse of 0");
        self.exp[((self.q - 1 - self.log[a as usize]) % (self.q - 1)) as usize]
    }
    /// `row ← row + f · piv` elementwise (`f ≠ 0`).
    #[inline(always)]
    pub fn axpy(&self, row: &mut [Fe], piv: &[Fe], f: Fe) {
        let lf = self.log[f as usize];
        for (r, &v) in row.iter_mut().zip(piv) {
            if v != 0 {
                *r ^= self.exp[(lf + self.log[v as usize]) as usize];
            }
        }
    }
}

/// `F_{2^m}`, `2 ≤ m ≤ 20`, by logarithm tables.  `mul` and `inv` count
/// one multiplication each on [`Gf2q::muls`]; the `_nc` forms are for
/// loops that count for themselves.
pub struct Gf2q {
    pub m: u32,
    pub q: u32,
    pub modulus: u32,
    pub tab: GfTab,
    /// `t` with `t² + t = w`, or `Fe::MAX` when `Tr(w) ≠ 0`.
    artin: Vec<Fe>,
    muls: Cell<u64>,
}

impl Gf2q {
    pub fn new(m: u32) -> Gf2q {
        assert!((2..=20).contains(&m), "2 ≤ m ≤ 20");
        let q = 1u32 << m;
        let modulus = PRIMITIVE[m as usize];
        let mut exp = vec![0u32; 2 * q as usize];
        let mut log = vec![0u32; q as usize];
        let mut v = 1u32;
        for i in 0..(q - 1) {
            exp[i as usize] = v as Fe;
            assert!(i == 0 || v != 1, "x is primitive modulo {modulus:#x}");
            log[v as usize] = i;
            v <<= 1;
            if v & q != 0 {
                v ^= modulus;
            }
        }
        assert_eq!(v, 1, "x has order q − 1");
        for i in (q - 1)..(2 * q) {
            exp[i as usize] = exp[(i - (q - 1)) as usize];
        }
        let mut f = Gf2q {
            m,
            q,
            modulus,
            tab: GfTab { q, exp, log },
            artin: vec![Fe::MAX; q as usize],
            muls: Cell::new(0),
        };
        for t in 0..q as Fe {
            let w = f.mul_nc(t, t) ^ t;
            f.artin[w as usize] = t.min(t ^ 1);
        }
        f
    }
    pub fn muls(&self) -> u64 {
        self.muls.get()
    }
    pub fn count(&self, k: u64) {
        self.muls.set(self.muls.get() + k);
    }
    #[inline(always)]
    pub fn mul_nc(&self, a: Fe, b: Fe) -> Fe {
        self.tab.mul(a, b)
    }
    #[inline(always)]
    pub fn mul(&self, a: Fe, b: Fe) -> Fe {
        self.count(1);
        self.tab.mul(a, b)
    }
    pub fn inv_nc(&self, a: Fe) -> Fe {
        self.tab.inv(a)
    }
    pub fn inv(&self, a: Fe) -> Fe {
        self.count(1);
        self.tab.inv(a)
    }
    /// The smaller `t` with `t² + t = w`, if any (a table: no arithmetic).
    pub fn artin_schreier(&self, w: Fe) -> Option<Fe> {
        let t = self.artin[w as usize];
        (t != Fe::MAX).then_some(t)
    }
    /// `√a = a^{2^{m−1}}`: `m − 1` squarings, counted.
    pub fn sqrt(&self, a: Fe) -> Fe {
        let mut r = a;
        for _ in 1..self.m {
            r = self.mul(r, r);
        }
        r
    }
}

// ── F_{q⁴} = F_q[z] / (z⁴ − c₃z³ − c₂z² − c₁z − c₀) ───────────────────────

/// An element of `F_{q⁴}`: coefficients of `1, z, z², z³` over `F_q`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Eb(pub [Fe; 4]);

impl Eb {
    pub const ZERO: Eb = Eb([0; 4]);
    pub const ONE: Eb = Eb([1, 0, 0, 0]);
    pub fn is_zero(&self) -> bool {
        self.0 == [0; 4]
    }
    pub fn in_fq(&self) -> bool {
        self.0[1] == 0 && self.0[2] == 0 && self.0[3] == 0
    }
    pub fn add(&self, o: &Eb) -> Eb {
        Eb(core::array::from_fn(|i| self.0[i] ^ o.0[i]))
    }
}

pub struct Fq4 {
    pub f: Gf2q,
    /// `z⁴ = c₀ + c₁z + c₂z² + c₃z³`.
    pub c: [Fe; 4],
    /// `(zⁱ)^q`, for the Frobenius.
    zq: [Eb; 4],
}

impl Fq4 {
    /// The first `z⁴ + c₃z³ + c₁z + c₀`, `c₃ ∈ {0, 1}` (by `(c₃, c₁, c₀)` in
    /// order) irreducible over `F_q`: no factor of degree 1 or 2, i.e.
    /// `z^{q²} − z` a unit modulo it.
    pub fn new(m: u32) -> Fq4 {
        assert!((4..=20).contains(&m), "4 ≤ m ≤ 20");
        let mut fq = Fq4 {
            f: Gf2q::new(m),
            c: [0; 4],
            zq: [Eb::ONE; 4],
        };
        // z⁴ + c₁z + c₀ is an affine linearised polynomial in characteristic
        // 2, often reducible for every (c₁, c₀); the z³ term avoids that.
        for (c3, c1) in (0..2u32).flat_map(|c3| (0..fq.f.q as Fe).map(move |c1| (c3, c1))) {
            for c0 in 1..fq.f.q as Fe {
                fq.c = [c0, c1, 0, c3];
                let z = Eb([0, 1, 0, 0]);
                let d = fq.pow2k(&z, 2 * m).add(&z);
                if d.is_zero() || fq.mat_rank(&d) != 4 {
                    continue;
                }
                let zq = fq.pow2k(&z, m);
                let mut acc = Eb::ONE;
                for i in 0..4 {
                    fq.zq[i] = acc;
                    acc = fq.mul_nc(&acc, &zq);
                }
                fq.f.muls.set(0);
                return fq;
            }
        }
        unreachable!("an irreducible quartic exists")
    }
    /// Rank over `F_q` of multiplication by `d` (4 means a unit).
    fn mat_rank(&self, d: &Eb) -> usize {
        let f = &self.f;
        let mut m: Vec<[Fe; 4]> = (0..4)
            .map(|i| {
                let mut e = Eb::ZERO;
                e.0[i] = 1;
                self.mul_nc(d, &e).0
            })
            .collect();
        let mut rank = 0;
        for c in 0..4 {
            let Some(r) = (rank..4).find(|&r| m[r][c] != 0) else {
                continue;
            };
            m.swap(r, rank);
            let inv = f.inv_nc(m[rank][c]);
            let piv = m[rank];
            for rr in 0..4 {
                if rr != rank && m[rr][c] != 0 {
                    let k = f.mul_nc(m[rr][c], inv);
                    for j in 0..4 {
                        m[rr][j] ^= f.mul_nc(k, piv[j]);
                    }
                }
            }
            rank += 1;
        }
        rank
    }
    pub fn from_fq(&self, v: Fe) -> Eb {
        Eb([v, 0, 0, 0])
    }
    fn reduce(&self, p: &mut [Fe; 7], counted: bool) -> Eb {
        let f = &self.f;
        for k in (4..7).rev() {
            let v = p[k];
            if v == 0 {
                continue;
            }
            for i in 0..4 {
                if self.c[i] != 0 {
                    p[k - 4 + i] ^= f.mul_nc(v, self.c[i]);
                    if counted {
                        f.count(1);
                    }
                }
            }
        }
        Eb([p[0], p[1], p[2], p[3]])
    }
    fn mul_impl(&self, a: &Eb, b: &Eb, counted: bool) -> Eb {
        let f = &self.f;
        let mut p = [0u32; 7];
        for i in 0..4 {
            if a.0[i] == 0 {
                continue;
            }
            for j in 0..4 {
                if b.0[j] != 0 {
                    p[i + j] ^= f.mul_nc(a.0[i], b.0[j]);
                    if counted {
                        f.count(1);
                    }
                }
            }
        }
        self.reduce(&mut p, counted)
    }
    pub fn mul(&self, a: &Eb, b: &Eb) -> Eb {
        self.mul_impl(a, b, true)
    }
    fn mul_nc(&self, a: &Eb, b: &Eb) -> Eb {
        self.mul_impl(a, b, false)
    }
    /// `a²`: the coefficients squared into even positions, then reduced.
    pub fn sq(&self, a: &Eb) -> Eb {
        let f = &self.f;
        let mut p = [0u32; 7];
        for i in 0..4 {
            if a.0[i] != 0 {
                p[2 * i] = f.mul(a.0[i], a.0[i]);
            }
        }
        self.reduce(&mut p, true)
    }
    fn sq_nc(&self, a: &Eb) -> Eb {
        let f = &self.f;
        let mut p = [0u32; 7];
        for i in 0..4 {
            p[2 * i] = f.mul_nc(a.0[i], a.0[i]);
        }
        self.reduce(&mut p, false)
    }
    /// `a^{2^k}` (uncounted: setup only).
    fn pow2k(&self, a: &Eb, k: u32) -> Eb {
        let mut r = *a;
        for _ in 0..k {
            r = self.sq_nc(&r);
        }
        r
    }
    pub fn scale(&self, a: &Eb, k: Fe) -> Eb {
        Eb(core::array::from_fn(|i| self.f.mul(a.0[i], k)))
    }
    /// `a^q`: `F_q`-linear, `Σ aᵢ (zⁱ)^q`.
    pub fn frob(&self, a: &Eb) -> Eb {
        let mut r = Eb([a.0[0], 0, 0, 0]);
        for i in 1..4 {
            if a.0[i] != 0 {
                r = r.add(&self.scale(&self.zq[i], a.0[i]));
            }
        }
        r
    }
    /// `a⁻¹ = a^{q + q² + q³} / N(a)`, `N(a) = a · a^{q + q² + q³} ∈ F_q`.
    pub fn inv(&self, a: &Eb) -> Eb {
        assert!(!a.is_zero(), "inverse of 0");
        let a1 = self.frob(a);
        let a2 = self.frob(&a1);
        let a3 = self.frob(&a2);
        let c = self.mul(&self.mul(&a1, &a2), &a3);
        let n = self.mul(a, &c);
        debug_assert!(n.in_fq(), "the norm lies in F_q");
        self.scale(&c, self.f.inv(n.0[0]))
    }
    /// Whether `a ∈ F_{q²}` (`a^{q²} = a`; uncounted, a diagnostic).
    pub fn in_fq2(&self, a: &Eb) -> bool {
        self.pow2k(a, 2 * self.f.m) == *a
    }
    /// `a^{2^k}`, counted.
    pub fn pow2k_counted(&self, a: &Eb, k: u32) -> Eb {
        let mut r = *a;
        for _ in 0..k {
            r = self.sq(&r);
        }
        r
    }
    /// `√a = a^{2^{4m − 1}}`.
    pub fn sqrt(&self, a: &Eb) -> Eb {
        self.pow2k_counted(a, 4 * self.f.m - 1)
    }
    pub fn random(&self, rng: &mut StdRng) -> Eb {
        Eb(core::array::from_fn(|_| rng.gen_range(0..self.f.q) as Fe))
    }
    /// Pack into one word: coefficient `i` in bits `[im, (i + 1)m)`
    /// (`4m ≤ 80` bits; below `2^64` while `m ≤ 16`).
    pub fn pack(&self, a: &Eb) -> u128 {
        (0..4).fold(0u128, |acc, i| {
            acc | (a.0[i] as u128) << (i as u32 * self.f.m)
        })
    }
    pub fn unpack(&self, v: u128) -> Eb {
        let mask = (1u128 << self.f.m) - 1;
        Eb(core::array::from_fn(|i| {
            ((v >> (i as u32 * self.f.m)) & mask) as Fe
        }))
    }
}

// ── The curve y² + xy = x³ + ax² + b over F_{q⁴} ─────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PtB {
    pub x: Eb,
    pub y: Eb,
    pub inf: bool,
}

impl PtB {
    pub const INF: PtB = PtB {
        x: Eb::ZERO,
        y: Eb::ZERO,
        inf: true,
    };
}

pub struct CurveB {
    pub f: Fq4,
    pub a: Eb,
    pub b: Eb,
    /// `β = b^{1/4}`.
    pub beta: Eb,
    /// `u ↦ u² + u` over `F_2` on the packed bits, row-reduced with its
    /// transformation, for the lift.
    as_solver: Vec<(u128, u128)>,
    ops: Cell<u64>,
}

impl CurveB {
    pub fn new(f: Fq4, a: Eb, b: Eb) -> CurveB {
        assert!(!b.is_zero(), "b ≠ 0: the curve is ordinary and smooth");
        let beta = f.pow2k(&b, 4 * f.f.m - 2);
        assert_eq!(f.sq_nc(&f.sq_nc(&beta)), b, "β⁴ = b");
        let n = 4 * f.f.m;
        // Rows: image of each basis bit, with the bit itself as the
        // transformation; reduce to echelon form by leading bit.
        let mut rows: Vec<(u128, u128)> = (0..n)
            .map(|i| {
                let e = f.unpack(1u128 << i);
                let img = f.sq_nc(&e).add(&e);
                (f.pack(&img), 1u128 << i)
            })
            .collect();
        let mut basis: Vec<(u128, u128)> = Vec::new();
        for (mut v, mut tr) in rows.drain(..) {
            // Basis kept by leading bit, descending; each lead unique.
            for &(bv, btr) in &basis {
                let lead = 127 - bv.leading_zeros();
                if v >> lead & 1 == 1 {
                    v ^= bv;
                    tr ^= btr;
                }
            }
            if v != 0 {
                basis.push((v, tr));
                basis.sort_by_key(|x| std::cmp::Reverse(127 - x.0.leading_zeros()));
            }
        }
        f.f.muls.set(0);
        CurveB {
            f,
            a,
            b,
            beta,
            as_solver: basis,
            ops: Cell::new(0),
        }
    }
    pub fn ops(&self) -> u64 {
        self.ops.get()
    }
    /// A `u` with `u² + u = c`, if `Tr(c) = 0` (linear algebra over `F_2`,
    /// no field multiplication).
    fn artin_schreier(&self, c: &Eb) -> Option<Eb> {
        let mut v = self.f.pack(c);
        let mut u = 0u128;
        for &(bv, btr) in &self.as_solver {
            let lead = 127 - bv.leading_zeros();
            if v >> lead & 1 == 1 {
                v ^= bv;
                u ^= btr;
            }
        }
        (v == 0).then(|| self.f.unpack(u))
    }
    /// A point with abscissa `x`, if one exists; the ordinate's choice of
    /// the two is the sampler's.
    pub fn lift_x(&self, x: &Eb, rng: &mut StdRng) -> Option<PtB> {
        let f = &self.f;
        if x.is_zero() {
            return Some(PtB {
                x: *x,
                y: f.sqrt(&self.b),
                inf: false,
            });
        }
        let x2 = f.sq(x);
        let c = x.add(&self.a).add(&f.mul(&self.b, &f.inv(&x2)));
        let u = self.artin_schreier(&c)?;
        let u = if rng.gen_bool(0.5) {
            u
        } else {
            u.add(&Eb::ONE)
        };
        Some(PtB {
            x: *x,
            y: f.mul(x, &u),
            inf: false,
        })
    }
    pub fn on_curve(&self, p: &PtB) -> bool {
        if p.inf {
            return true;
        }
        let f = &self.f;
        let lhs = f.sq(&p.y).add(&f.mul(&p.x, &p.y));
        let x2 = f.sq(&p.x);
        let rhs = f.mul(&x2, &p.x).add(&f.mul(&self.a, &x2)).add(&self.b);
        lhs == rhs
    }
    pub fn neg(&self, p: &PtB) -> PtB {
        if p.inf {
            return *p;
        }
        PtB {
            x: p.x,
            y: p.x.add(&p.y),
            inf: false,
        }
    }
    pub fn add(&self, p: &PtB, q: &PtB) -> PtB {
        self.ops.set(self.ops.get() + 1);
        if p.inf {
            return *q;
        }
        if q.inf {
            return *p;
        }
        let f = &self.f;
        if p.x == q.x {
            if p.y != q.y || p.x.is_zero() {
                return PtB::INF;
            }
            // Doubling: λ = x + y/x, x₃ = λ² + λ + a, y₃ = x² + (λ + 1)x₃.
            let l = p.x.add(&f.mul(&p.y, &f.inv(&p.x)));
            let x3 = f.sq(&l).add(&l).add(&self.a);
            let y3 = f.sq(&p.x).add(&f.mul(&l.add(&Eb::ONE), &x3));
            return PtB {
                x: x3,
                y: y3,
                inf: false,
            };
        }
        let l = f.mul(&p.y.add(&q.y), &f.inv(&p.x.add(&q.x)));
        let x3 = f.sq(&l).add(&l).add(&p.x).add(&q.x).add(&self.a);
        let y3 = f.mul(&l, &p.x.add(&x3)).add(&x3).add(&p.y);
        PtB {
            x: x3,
            y: y3,
            inf: false,
        }
    }
    pub fn sub(&self, p: &PtB, q: &PtB) -> PtB {
        self.add(p, &self.neg(q))
    }
    pub fn mul(&self, p: &PtB, mut k: u128) -> PtB {
        let mut acc = PtB::INF;
        let mut base = *p;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            base = self.add(&base, &base);
            k >>= 1;
        }
        acc
    }
    /// `T = (0, √b)`.
    pub fn torsion(&self) -> PtB {
        PtB {
            x: Eb::ZERO,
            y: self.f.sqrt(&self.b),
            inf: false,
        }
    }
    /// `t = β / (x + β)` (`x = β` is the order-4 point, `t = ∞`).
    pub fn t_of_x(&self, x: &Eb) -> Eb {
        let f = &self.f;
        f.mul(&self.beta, &f.inv(&x.add(&self.beta)))
    }
    /// `x = β (t + 1) / t`.
    pub fn x_of_t(&self, t: &Eb) -> Eb {
        let f = &self.f;
        f.mul(&f.mul(&self.beta, &t.add(&Eb::ONE)), &f.inv(t))
    }

    /// `S₃` as a quadratic in its last argument `U`: `(x₁² + x₂²)U² +
    /// x₁x₂U + x₁²x₂² + b`.
    fn s3_quadratic(&self, x1: &Eb, x2: &Eb) -> [Eb; 3] {
        let f = &self.f;
        let s1 = f.sq(x1);
        let s2 = f.sq(x2);
        [f.mul(&s1, &s2).add(&self.b), f.mul(x1, x2), s1.add(&s2)]
    }
    /// `S₃(x₁, x₂, x₃) = x₁²x₂² + x₁²x₃² + x₂²x₃² + x₁x₂x₃ + b`.
    pub fn s3(&self, x1: &Eb, x2: &Eb, x3: &Eb) -> Eb {
        let f = &self.f;
        let [c, bq, a] = self.s3_quadratic(x1, x2);
        f.mul(&a, &f.sq(x3)).add(&f.mul(&bq, x3)).add(&c)
    }
    /// `S₄ = Res_U(S₃(x₁, x₂, U), S₃(x₃, x₄, U))`: for quadratics
    /// `(A, B, C)`, `(A₁C₂ + A₂C₁)² + (A₁B₂ + A₂B₁)(B₁C₂ + B₂C₁)`.
    pub fn s4(&self, x1: &Eb, x2: &Eb, x3: &Eb, x4: &Eb) -> Eb {
        let f = &self.f;
        let [c1, b1, a1] = self.s3_quadratic(x1, x2);
        let [c2, b2, a2] = self.s3_quadratic(x3, x4);
        let ac = f.mul(&a1, &c2).add(&f.mul(&a2, &c1));
        let ab = f.mul(&a1, &b2).add(&f.mul(&a2, &b1));
        let bc = f.mul(&b1, &c2).add(&f.mul(&b2, &c1));
        f.sq(&ac).add(&f.mul(&ab, &bc))
    }
    /// `S₅(x₁, …, x₅) = Res_U(S₃(x₁, x₂, U), S₄(x₃, x₄, x₅, U))`, `S₄` of
    /// formal degree 4 in `U` interpolated on five nodes, the resultant as
    /// the 6 × 6 Sylvester determinant.
    pub fn s5(&self, x: &[Eb; 5]) -> Eb {
        let f = &self.f;
        let nodes: [Eb; 5] = core::array::from_fn(|i| Eb([i as Fe, 1, 0, 0]));
        let vals: [Eb; 5] = core::array::from_fn(|i| self.s4(&x[2], &x[3], &x[4], &nodes[i]));
        let g = interpolate(f, &nodes, &vals);
        let [c, bq, a] = self.s3_quadratic(&x[0], &x[1]);
        // Sylvester: 4 rows of (A, B, C) shifted, 2 rows of g (degree 4).
        let mut m = vec![vec![Eb::ZERO; 6]; 6];
        for r in 0..4 {
            m[r][r] = a;
            m[r][r + 1] = bq;
            m[r][r + 2] = c;
        }
        for r in 0..2 {
            for k in 0..5 {
                m[4 + r][r + k] = g[4 - k];
            }
        }
        det_e(f, m)
    }
}

/// The polynomial of degree `< n` through `(nodes[i], vals[i])`, ascending.
fn interpolate(f: &Fq4, nodes: &[Eb], vals: &[Eb]) -> Vec<Eb> {
    let n = nodes.len();
    let mut out = vec![Eb::ZERO; n];
    for i in 0..n {
        // basis polynomial Π_{j≠i} (U − u_j) / (u_i − u_j)
        let mut poly = vec![Eb::ONE];
        let mut den = Eb::ONE;
        for j in 0..n {
            if j == i {
                continue;
            }
            let mut next = vec![Eb::ZERO; poly.len() + 1];
            for (k, c) in poly.iter().enumerate() {
                next[k + 1] = next[k + 1].add(c);
                next[k] = next[k].add(&f.mul(c, &nodes[j]));
            }
            poly = next;
            den = f.mul(&den, &nodes[i].add(&nodes[j]));
        }
        let s = f.mul(&vals[i], &f.inv(&den));
        for (k, c) in poly.iter().enumerate() {
            out[k] = out[k].add(&f.mul(c, &s));
        }
    }
    out
}

fn det_e(f: &Fq4, mut m: Vec<Vec<Eb>>) -> Eb {
    let n = m.len();
    let mut det = Eb::ONE;
    for c in 0..n {
        let Some(r) = (c..n).find(|&r| !m[r][c].is_zero()) else {
            return Eb::ZERO;
        };
        m.swap(r, c);
        det = f.mul(&det, &m[c][c]);
        let inv = f.inv(&m[c][c]);
        for r in (c + 1)..n {
            if m[r][c].is_zero() {
                continue;
            }
            let k = f.mul(&m[r][c], &inv);
            for j in c..n {
                let v = f.mul(&k, &m[c][j]);
                m[r][j] = m[r][j].add(&v);
            }
        }
    }
    det
}

// ── The D₄-symmetrised S₅ and its Weil restriction ───────────────────────

/// `s₂^a s₃^b s₄^c π^d` with `2(a + b + c) + d ≤ 8`: 105 monomials, as
/// `[a, b, c, d]`.
pub fn invariant_monomials() -> Vec<[u8; 4]> {
    let mut v = Vec::new();
    for d in 0..=8u8 {
        let k = (8 - d) / 2;
        for a in 0..=k {
            for b in 0..=(k - a) {
                for c in 0..=(k - a - b) {
                    v.push([a, b, c, d]);
                }
            }
        }
    }
    v
}

/// `(s₂, s₃, s₄, π)` of `(t₁, …, t₄)`.
fn invariants(f: &Fq4, t: &[Eb; 4]) -> [Eb; 4] {
    let w: Vec<Eb> = t.iter().map(|v| f.sq(v).add(v)).collect();
    let mut e = [Eb::ONE, Eb::ZERO, Eb::ZERO, Eb::ZERO, Eb::ZERO];
    for wi in &w {
        for k in (1..5).rev() {
            e[k] = e[k].add(&f.mul(&e[k - 1], wi));
        }
    }
    let pi = t[0].add(&t[1]).add(&t[2]).add(&t[3]);
    [e[2], e[3], e[4], pi]
}

fn mono_values(f: &Fq4, inv: &[Eb; 4], monos: &[[u8; 4]]) -> Vec<Eb> {
    let mut pows = [[Eb::ONE; 9]; 4];
    for i in 0..4 {
        for k in 1..9 {
            pows[i][k] = f.mul(&pows[i][k - 1], &inv[i]);
        }
    }
    monos
        .iter()
        .map(|m| {
            f.mul(
                &f.mul(&pows[0][m[0] as usize], &pows[1][m[1] as usize]),
                &f.mul(&pows[2][m[2] as usize], &pows[3][m[3] as usize]),
            )
        })
        .collect()
}

/// `H = Σ_m coef[m](t_R) · m(s₂, s₃, s₄, π)`, nine `F_{q⁴}` coefficients in
/// `t_R` per monomial.
pub struct SymmetrisedBin {
    pub monos: Vec<[u8; 4]>,
    pub coef: Vec<[Eb; 9]>,
    pub checked_points: usize,
    pub precompute_muls: u64,
}

impl SymmetrisedBin {
    fn h(curve: &CurveB, t: &[Eb; 4], tr: &Eb) -> Eb {
        let f = &curve.f;
        let x: [Eb; 5] = [
            curve.x_of_t(&t[0]),
            curve.x_of_t(&t[1]),
            curve.x_of_t(&t[2]),
            curve.x_of_t(&t[3]),
            curve.x_of_t(tr),
        ];
        let mut v = curve.s5(&x);
        for ti in t.iter().chain(std::iter::once(tr)) {
            v = f.mul(&v, &f.pow2k_counted(ti, 3));
        }
        v
    }

    /// Nine distinct nonzero nodes for `t_R`, outside `F_q`.
    fn nodes() -> [Eb; 9] {
        core::array::from_fn(|i| Eb([i as Fe, 1, 0, 0]))
    }

    fn coefficients_in_last(curve: &CurveB, t: &[Eb; 4], vinv: &[[Eb; 9]; 9]) -> [Eb; 9] {
        let f = &curve.f;
        let vals: Vec<Eb> = Self::nodes().iter().map(|n| Self::h(curve, t, n)).collect();
        core::array::from_fn(|j| {
            vals.iter()
                .enumerate()
                .fold(Eb::ZERO, |acc, (k, v)| acc.add(&f.mul(v, &vinv[j][k])))
        })
    }

    fn random_t(f: &Fq4, rng: &mut StdRng) -> [Eb; 4] {
        core::array::from_fn(|_| loop {
            let v = f.random(rng);
            if !v.is_zero() && v != Eb::ONE {
                break v;
            }
        })
    }

    /// Interpolate on 105 random points and check on `checks` fresh ones;
    /// `Err` if `H` is not in the span of the invariant monomials.
    pub fn precompute(curve: &CurveB, rng: &mut StdRng, checks: usize) -> Result<Self, String> {
        let f = &curve.f;
        let before = f.f.muls();
        let nodes = Self::nodes();
        // Inverse Vandermonde: coefficient j of the polynomial through the
        // nodes is Σ_k vinv[j][k] · value_k.
        let vinv: [[Eb; 9]; 9] = {
            let mut out = [[Eb::ZERO; 9]; 9];
            for k in 0..9 {
                let vals: Vec<Eb> = (0..9)
                    .map(|i| if i == k { Eb::ONE } else { Eb::ZERO })
                    .collect();
                let poly = interpolate(f, &nodes, &vals);
                for j in 0..9 {
                    out[j][k] = poly[j];
                }
            }
            out
        };
        let monos = invariant_monomials();
        let n = monos.len();
        for _attempt in 0..8 {
            let mut sys: Vec<Vec<Eb>> = Vec::with_capacity(n);
            for _ in 0..n {
                let t = Self::random_t(f, rng);
                let mut row = mono_values(f, &invariants(f, &t), &monos);
                row.extend(Self::coefficients_in_last(curve, &t, &vinv));
                sys.push(row);
            }
            let mut ok = true;
            for c in 0..n {
                let Some(r) = (c..n).find(|&r| !sys[r][c].is_zero()) else {
                    ok = false;
                    break;
                };
                sys.swap(r, c);
                let inv = f.inv(&sys[c][c]);
                let pivot: Vec<Eb> = sys[c].iter().map(|v| f.mul(v, &inv)).collect();
                sys[c] = pivot.clone();
                for r in 0..n {
                    if r == c || sys[r][c].is_zero() {
                        continue;
                    }
                    let k = sys[r][c];
                    for j in c..pivot.len() {
                        if !pivot[j].is_zero() {
                            let v = f.mul(&k, &pivot[j]);
                            sys[r][j] = sys[r][j].add(&v);
                        }
                    }
                }
            }
            if !ok {
                continue;
            }
            let coef: Vec<[Eb; 9]> = (0..n)
                .map(|r| core::array::from_fn(|j| sys[r][n + j]))
                .collect();
            for t_ in 0..checks {
                let t = Self::random_t(f, rng);
                let want = Self::coefficients_in_last(curve, &t, &vinv);
                let mv = mono_values(f, &invariants(f, &t), &monos);
                for j in 0..9 {
                    let got = mv
                        .iter()
                        .zip(&coef)
                        .fold(Eb::ZERO, |acc, (m, c)| acc.add(&f.mul(m, &c[j])));
                    if got != want[j] {
                        return Err(format!(
                            "H is not in the span of the invariants (check point {t_}, coefficient t_R^{j})"
                        ));
                    }
                }
            }
            return Ok(SymmetrisedBin {
                monos,
                coef,
                checked_points: checks,
                precompute_muls: f.f.muls() - before,
            });
        }
        Err("the interpolation system was singular eight times".into())
    }

    /// `H(·; t_R)` at every monomial, in `F_{q⁴}`.
    fn at(&self, f: &Fq4, t_r: &Eb) -> Vec<Eb> {
        let mut pows = [Eb::ONE; 9];
        for k in 1..9 {
            pows[k] = f.mul(&pows[k - 1], t_r);
        }
        self.coef
            .iter()
            .map(|c| (0..9).fold(Eb::ZERO, |acc, j| acc.add(&f.mul(&c[j], &pows[j]))))
            .collect()
    }

    /// The four Weil components in `(s₂, s₃, s₄, π)`: weighted form,
    /// weights [`WEIGHTS`].
    pub fn system_weighted(&self, f: &Fq4, t_r: &Eb) -> Vec<HashMap<[u8; 4], Fe>> {
        let vals = self.at(f, t_r);
        let mut out: Vec<HashMap<[u8; 4], Fe>> = vec![HashMap::new(); 4];
        for (m, v) in self.monos.iter().zip(&vals) {
            for (i, comp) in out.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    *comp.entry(*m).or_insert(0) ^= v.0[i];
                }
            }
        }
        for c in &mut out {
            c.retain(|_, v| *v != 0);
        }
        out
    }

    /// The unweighted form in `(s₁, s₂, s₃, s₄, π)`: `π^d` reduced modulo
    /// `π² + π + s₁` to `A_d(s₁) + B_d(s₁) π` (so every term has total
    /// degree ≤ 4), and the fifth equation `π² + π + s₁`.
    pub fn system_plain(&self, f: &Fq4, t_r: &Eb) -> Vec<HashMap<[u8; 5], Fe>> {
        // π^d = Σ_k A[d][k] s₁^k + Σ_k B[d][k] s₁^k π, coefficients in F_2.
        let mut ab: Vec<(Vec<u8>, Vec<u8>)> = vec![(vec![1], vec![0])];
        for d in 1..=8usize {
            let (a, b) = &ab[d - 1];
            // π^{d} = A π + B π² = A π + B (π + s₁) = B s₁ + (A + B) π
            let mut na = vec![0u8; b.len() + 1];
            for (k, &c) in b.iter().enumerate() {
                na[k + 1] ^= c;
            }
            let len = a.len().max(b.len());
            let nb: Vec<u8> = (0..len)
                .map(|k| a.get(k).copied().unwrap_or(0) ^ b.get(k).copied().unwrap_or(0))
                .collect();
            ab.push((na, nb));
        }
        let vals = self.at(f, t_r);
        let mut out: Vec<HashMap<[u8; 5], Fe>> = vec![HashMap::new(); 4];
        for (m, v) in self.monos.iter().zip(&vals) {
            let (a, b) = &ab[m[3] as usize];
            for (part, pe) in [(a, 0u8), (b, 1u8)] {
                for (k, &c) in part.iter().enumerate() {
                    if c == 0 {
                        continue;
                    }
                    let key = [k as u8, m[0], m[1], m[2], pe];
                    for (i, comp) in out.iter_mut().enumerate() {
                        if v.0[i] != 0 {
                            *comp.entry(key).or_insert(0) ^= v.0[i];
                        }
                    }
                }
            }
        }
        for c in &mut out {
            c.retain(|_, v| *v != 0);
        }
        let mut rel = HashMap::new();
        rel.insert([0, 0, 0, 0, 2], 1);
        rel.insert([0, 0, 0, 0, 1], 1);
        rel.insert([1, 0, 0, 0, 0], 1);
        out.push(rel);
        out
    }
}

/// The weights of `(s₂, s₃, s₄, π)`.
pub const WEIGHTS: [u8; 4] = [2, 2, 2, 1];

/// The weights of `(S₂, S₃, s₄, P) = (s₂², s₃², s₄, π²)` in the
/// Frobenius-reduced form.
pub const FROBENIUS_WEIGHTS: [u8; 4] = [2, 2, 1, 1];

impl SymmetrisedBin {
    /// Whether every monomial carrying `H` has even exponents in `s₂`, `s₃`
    /// and `π` (measured on every instance from the interpolated form: true
    /// on every curve tried, not proved).  Then, since squaring is a
    /// bijection of `F_q`, `H` is a polynomial in `S₂ = s₂²`, `S₃ = s₃²`,
    /// `s₄` and `P = π²` with the same `F_q`-rational solutions: four
    /// equations of weighted degree 4 in the weights `(2, 2, 1, 1)`, Hilbert
    /// series `(1 + t²)²(1 + t + t² + t³)²`, 64 solutions over the closure
    /// against 512, weighted Macaulay degree 11.
    pub fn frobenius_reducible(&self) -> bool {
        self.monos
            .iter()
            .zip(&self.coef)
            .filter(|(_, c)| c.iter().any(|v| !v.is_zero()))
            .all(|(m, _)| m[0] % 2 == 0 && m[1] % 2 == 0 && m[3] % 2 == 0)
    }

    /// The four Weil components in `(S₂, S₃, s₄, P)` (requires
    /// [`Self::frobenius_reducible`]).
    pub fn system_frobenius(&self, f: &Fq4, t_r: &Eb) -> Vec<HashMap<[u8; 4], Fe>> {
        assert!(
            self.frobenius_reducible(),
            "H has an odd power of s₂, s₃ or π"
        );
        let vals = self.at(f, t_r);
        let mut out: Vec<HashMap<[u8; 4], Fe>> = vec![HashMap::new(); 4];
        for (m, v) in self.monos.iter().zip(&vals) {
            if v.is_zero() {
                continue;
            }
            let key = [m[0] / 2, m[1] / 2, m[2], m[3] / 2];
            for (i, comp) in out.iter_mut().enumerate() {
                if v.0[i] != 0 {
                    *comp.entry(key).or_insert(0) ^= v.0[i];
                }
            }
        }
        for c in &mut out {
            c.retain(|_, v| *v != 0);
        }
        out
    }
}

// ── Instances, the base, and the independent enumeration ─────────────────

/// A random ordinary curve over `F_{q⁴}` with `j = 1/b` outside `F_{q²}`
/// (no subfield structure).  Deterministic in `(m, seed)`.
pub fn generate_instance(m: u32, seed: u64) -> CurveB {
    let mut rng = StdRng::seed_from_u64(seed ^ 0xb1d4);
    let f = Fq4::new(m);
    loop {
        let b = f.random(&mut rng);
        // b ∉ F_{q²}: b^{q²} ≠ b.
        if b.is_zero() || f.pow2k(&b, 2 * m) == b {
            continue;
        }
        let a = f.random(&mut rng);
        return CurveB::new(f, a, b);
    }
}

/// The base: one point per `t ∈ F_q \ {0, 1}` that lifts, with its `t`.
pub fn t_base(curve: &CurveB, rng: &mut StdRng) -> Vec<(Fe, PtB)> {
    let f = &curve.f;
    (2..f.f.q as Fe)
        .filter_map(|t| {
            let x = curve.x_of_t(&f.from_fq(t));
            curve.lift_x(&x, rng).map(|p| (t, p))
        })
        .collect()
}

/// Every set of four distinct base points with `±P_i ± P_j ± P_k ± P_l =
/// R`, as sorted `t`-quadruples, by meet in the middle over pairs.
pub fn mitm_t_decompositions(curve: &CurveB, base: &[(Fe, PtB)], r: &PtB) -> BTreeSet<[Fe; 4]> {
    let n = base.len();
    let mut table: HashMap<Eb, Vec<(usize, usize)>> = HashMap::new();
    for i in 0..n {
        for j in (i + 1)..n {
            for q in [
                curve.add(&base[i].1, &base[j].1),
                curve.sub(&base[i].1, &base[j].1),
            ] {
                if !q.inf {
                    table.entry(q.x).or_default().push((i, j));
                }
            }
        }
    }
    let mut out = BTreeSet::new();
    for k in 0..n {
        for l in (k + 1)..n {
            for (sk, sl) in [(false, false), (false, true), (true, false), (true, true)] {
                let pk = if sk { curve.neg(&base[k].1) } else { base[k].1 };
                let pl = if sl { curve.neg(&base[l].1) } else { base[l].1 };
                let v = curve.sub(&curve.sub(r, &pk), &pl);
                if v.inf {
                    continue;
                }
                if let Some(pairs) = table.get(&v.x) {
                    for &(i, j) in pairs {
                        if i == k || i == l || j == k || j == l {
                            continue;
                        }
                        let mut ts = [base[i].0, base[j].0, base[k].0, base[l].0];
                        ts.sort_unstable();
                        out.insert(ts);
                    }
                }
            }
        }
    }
    out
}

/// The signs that make the four points sum to `r`, applied, if any do.
fn signed_sum(curve: &CurveB, pts: &[PtB; 4], r: &PtB) -> Option<[PtB; 4]> {
    for s in 0..16u32 {
        let signed: [PtB; 4] = core::array::from_fn(|i| {
            if s >> i & 1 == 1 {
                curve.neg(&pts[i])
            } else {
                pts[i]
            }
        });
        let mut acc = PtB::INF;
        for q in &signed {
            acc = curve.add(&acc, q);
        }
        if !acc.inf && acc.x == r.x && acc.y == r.y {
            return Some(signed);
        }
    }
    None
}

// ── Weighted monomials ───────────────────────────────────────────────────

fn wdeg<const N: usize>(m: &[u8; N], w: &[u8; N]) -> Fe {
    m.iter().zip(w).map(|(&x, &k)| x as Fe * k as Fe).sum()
}

/// Weighted degree, then reverse lexicographic (grevlex for weights 1).
fn wcmp<const N: usize>(a: &[u8; N], b: &[u8; N], w: &[u8; N]) -> std::cmp::Ordering {
    let (da, db) = (wdeg(a, w), wdeg(b, w));
    if da != db {
        return da.cmp(&db);
    }
    for i in (0..N).rev() {
        if a[i] != b[i] {
            return b[i].cmp(&a[i]);
        }
    }
    std::cmp::Ordering::Equal
}

/// Monomials of weighted degree `≤ d`, descending.
fn monomials_desc<const N: usize>(d: u8, w: &[u8; N]) -> Vec<[u8; N]> {
    fn rec<const N: usize>(
        i: usize,
        left: u8,
        w: &[u8; N],
        cur: &mut [u8; N],
        out: &mut Vec<[u8; N]>,
    ) {
        if i == N {
            out.push(*cur);
            return;
        }
        for k in 0..=left / w[i] {
            cur[i] = k;
            rec(i + 1, left - k * w[i], w, cur, out);
        }
        cur[i] = 0;
    }
    let mut out = Vec::new();
    rec(0, d, w, &mut [0u8; N], &mut out);
    out.sort_by(|x, y| wcmp(y, x, w));
    out
}

fn lead<const N: usize>(f: &HashMap<[u8; N], Fe>, w: &[u8; N]) -> Option<[u8; N]> {
    f.iter()
        .filter(|(_, &c)| c != 0)
        .map(|(m, _)| *m)
        .max_by(|a, b| wcmp(a, b, w))
}

/// Make the leading monomials distinct by eliminating equal leads.
fn triangularise<const N: usize>(
    g: &GfTab,
    comps: &[HashMap<[u8; N], Fe>],
    w: &[u8; N],
    muls: &mut u64,
) -> (Vec<HashMap<[u8; N], Fe>>, Vec<Option<[u8; N]>>) {
    let mut fs = comps.to_vec();
    for _ in 0..16 {
        fs.sort_by(|a, b| match (lead(a, w), lead(b, w)) {
            (Some(x), Some(y)) => wcmp(&y, &x, w),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        });
        let mut changed = false;
        for i in 1..fs.len() {
            let (Some(li), Some(lj)) = (lead(&fs[i], w), lead(&fs[i - 1], w)) else {
                continue;
            };
            if li == lj {
                let k = g.mul(fs[i][&li], g.inv(fs[i - 1][&lj]));
                *muls += 1;
                let prev = fs[i - 1].clone();
                for (&m, &c) in &prev {
                    *fs[i].entry(m).or_insert(0) ^= g.mul(k, c);
                    *muls += 1;
                }
                fs[i].retain(|_, c| *c != 0);
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let leads = fs.iter().map(|f| lead(f, w)).collect();
    (fs, leads)
}

// ── Elimination ──────────────────────────────────────────────────────────

struct Echelon {
    rows: Vec<Vec<Fe>>,
    pivot_row: Vec<u32>,
    zero_rows: usize,
    zero_ids: Vec<u32>,
}

/// Forward elimination by leading column (as `gaudry_quartic`'s
/// `echelon_by_lead_with`): one multiplication per nonzero entry of a pivot
/// row per row it is subtracted from, and per entry normalised.  With
/// `sparsest`, each column's pivot is the sparsest candidate row.
fn echelon(g: &GfTab, mut rows: Vec<Vec<Fe>>, muls: &mut u64, sparsest: bool) -> Echelon {
    let ncols = rows.first().map_or(0, |r| r.len());
    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); ncols];
    let mut zero_rows = 0usize;
    let mut zero_ids = Vec::new();
    let mut nnz: Vec<u32> = rows
        .par_iter()
        .map(|r| r.iter().filter(|&&x| x != 0).count() as u32)
        .collect();
    for (r, row) in rows.iter().enumerate() {
        match row.iter().position(|&x| x != 0) {
            Some(c) => buckets[c].push(r as u32),
            None => {
                zero_rows += 1;
                zero_ids.push(r as u32);
            }
        }
    }
    let mut pivot_row = vec![u32::MAX; ncols];
    for c in 0..ncols {
        let mut b = std::mem::take(&mut buckets[c]);
        if b.is_empty() {
            continue;
        }
        b.sort_unstable();
        if sparsest && b.len() > 1 {
            let best = (0..b.len())
                .min_by_key(|&i| (nnz[b[i] as usize], i))
                .unwrap_or(0);
            b.swap(0, best);
        }
        let pr = b[0] as usize;
        let mut piv = std::mem::take(&mut rows[pr]);
        let inv = g.inv(piv[c]);
        let mut nz: Vec<u32> = Vec::new();
        for j in c..ncols {
            if piv[j] != 0 {
                piv[j] = g.mul(piv[j], inv);
                nz.push(j as u32);
            }
        }
        *muls += nz.len() as u64;
        let dense = nz.len() * 8 > ncols - c;
        let mut others: Vec<(u32, Vec<Fe>, Option<usize>, u32)> = b[1..]
            .iter()
            .map(|&r| (r, std::mem::take(&mut rows[r as usize]), None, 0))
            .collect();
        *muls += (others.len() * nz.len()) as u64;
        others.par_iter_mut().for_each(|(_, row, ld, cnt)| {
            let f = row[c];
            if dense {
                g.axpy(&mut row[c..], &piv[c..], f);
            } else {
                for &j in &nz {
                    let j = j as usize;
                    row[j] ^= g.mul(f, piv[j]);
                }
            }
            *ld = row[c + 1..].iter().position(|&x| x != 0).map(|k| k + c + 1);
            if let Some(l) = *ld {
                *cnt = row[l..].iter().filter(|&&x| x != 0).count() as u32;
            }
        });
        for (r, row, ld, cnt) in others {
            nnz[r as usize] = cnt;
            rows[r as usize] = row;
            match ld {
                Some(l) => buckets[l].push(r),
                None => {
                    zero_rows += 1;
                    zero_ids.push(r);
                }
            }
        }
        rows[pr] = piv;
        pivot_row[c] = pr as u32;
    }
    zero_ids.sort_unstable();
    Echelon {
        rows,
        pivot_row,
        zero_rows,
        zero_ids,
    }
}

/// Reduce held rows into the echelon until `need` new pivots are found.
fn top_up(
    g: &GfTab,
    ech: &mut Echelon,
    held: Vec<Vec<Fe>>,
    mut need: usize,
    muls: &mut u64,
) -> u64 {
    let ncols = ech.pivot_row.len();
    let mut used = 0;
    for mut v in held {
        if need == 0 {
            break;
        }
        used += 1;
        let mut ld = None;
        for c in 0..ncols {
            if v[c] == 0 {
                continue;
            }
            let pr = ech.pivot_row[c];
            if pr == u32::MAX {
                ld = Some(c);
                break;
            }
            let f = v[c];
            let piv = &ech.rows[pr as usize];
            *muls += piv[c..].iter().filter(|&&x| x != 0).count() as u64;
            g.axpy(&mut v[c..], &piv[c..], f);
        }
        let Some(c) = ld else {
            ech.zero_rows += 1;
            continue;
        };
        let inv = g.inv(v[c]);
        for x in &mut v[c..] {
            if *x != 0 {
                *x = g.mul(*x, inv);
                *muls += 1;
            }
        }
        ech.pivot_row[c] = ech.rows.len() as u32;
        ech.rows.push(v);
        need -= 1;
    }
    used
}

// ── Small linear algebra and polynomials over F_q (counted) ──────────────

/// Reduced row echelon form; returns pivot columns.
fn rref(g: &GfTab, m: &mut [Vec<Fe>], muls: &mut u64) -> Vec<usize> {
    let rows = m.len();
    let cols = if rows > 0 { m[0].len() } else { 0 };
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r >= rows {
            break;
        }
        let Some(pr) = (r..rows).find(|&i| m[i][c] != 0) else {
            continue;
        };
        m.swap(r, pr);
        let inv = g.inv(m[r][c]);
        for v in &mut m[r][c..] {
            *v = g.mul(*v, inv);
        }
        *muls += (cols - c) as u64;
        let piv = m[r].clone();
        for i in 0..rows {
            if i != r && m[i][c] != 0 {
                let f = m[i][c];
                g.axpy(&mut m[i][c..], &piv[c..], f);
                *muls += (cols - c) as u64;
            }
        }
        pivots.push(c);
        r += 1;
    }
    pivots
}

/// A basis of the kernel of `m − λI` (`m` square).
fn kernel_shifted(g: &GfTab, m: &[Vec<Fe>], lambda: Fe, muls: &mut u64) -> Vec<Vec<Fe>> {
    let n = m.len();
    let mut a = m.to_vec();
    for (i, row) in a.iter_mut().enumerate() {
        row[i] ^= lambda;
    }
    let pivots = rref(g, &mut a, muls);
    let mut is_piv = vec![false; n];
    for &c in &pivots {
        is_piv[c] = true;
    }
    (0..n)
        .filter(|&c| !is_piv[c])
        .map(|fc| {
            let mut v = vec![0u32; n];
            v[fc] = 1;
            for (r, &pc) in pivots.iter().enumerate() {
                v[pc] = a[r][fc];
            }
            v
        })
        .collect()
}

fn det(g: &GfTab, mut a: Vec<Vec<Fe>>, muls: &mut u64) -> Fe {
    let k = a.len();
    let mut d = 1u32;
    for c in 0..k {
        let Some(r) = (c..k).find(|&r| a[r][c] != 0) else {
            return 0;
        };
        a.swap(r, c);
        d = g.mul(d, a[c][c]);
        let inv = g.inv(a[c][c]);
        *muls += 1;
        for r in (c + 1)..k {
            if a[r][c] == 0 {
                continue;
            }
            let f = g.mul(a[r][c], inv);
            let piv = a[c].clone();
            g.axpy(&mut a[r][c + 1..], &piv[c + 1..], f);
            *muls += (k - c) as u64;
        }
    }
    d
}

fn ptrim(a: &mut Vec<Fe>) {
    while a.last() == Some(&0) {
        a.pop();
    }
}

/// `a mod h` (`h` monic).
fn pmod(g: &GfTab, a: &[Fe], h: &[Fe], muls: &mut u64) -> Vec<Fe> {
    let mut r = a.to_vec();
    ptrim(&mut r);
    let dh = h.len() - 1;
    while r.len() > dh {
        let c = *r.last().unwrap();
        let shift = r.len() - 1 - dh;
        for (i, &hi) in h.iter().enumerate() {
            if hi != 0 {
                r[shift + i] ^= g.mul(c, hi);
                *muls += 1;
            }
        }
        ptrim(&mut r);
    }
    r
}

fn pmul(g: &GfTab, a: &[Fe], b: &[Fe], muls: &mut u64) -> Vec<Fe> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0u32; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            if y != 0 {
                out[i + j] ^= g.mul(x, y);
                *muls += 1;
            }
        }
    }
    out
}

fn pmonic(g: &GfTab, mut a: Vec<Fe>, muls: &mut u64) -> Vec<Fe> {
    ptrim(&mut a);
    if let Some(&l) = a.last() {
        if l != 1 {
            let li = g.inv(l);
            for c in a.iter_mut() {
                *c = g.mul(*c, li);
            }
            *muls += a.len() as u64;
        }
    }
    a
}

fn pgcd(g: &GfTab, a: &[Fe], b: &[Fe], muls: &mut u64) -> Vec<Fe> {
    let mut a = pmonic(g, a.to_vec(), muls);
    let mut b = pmonic(g, b.to_vec(), muls);
    while !b.is_empty() {
        let r = pmod(g, &a, &b, muls);
        a = b;
        b = pmonic(g, r, muls);
    }
    a
}

fn pdiv(g: &GfTab, a: &[Fe], h: &[Fe], muls: &mut u64) -> Vec<Fe> {
    let mut r = a.to_vec();
    ptrim(&mut r);
    let dh = h.len() - 1;
    let mut qt = vec![0u32; r.len().saturating_sub(dh)];
    while r.len() > dh {
        let c = *r.last().unwrap();
        let shift = r.len() - 1 - dh;
        qt[shift] = c;
        for (i, &hi) in h.iter().enumerate() {
            if hi != 0 {
                r[shift + i] ^= g.mul(c, hi);
                *muls += 1;
            }
        }
        ptrim(&mut r);
    }
    qt
}

/// Distinct roots in `F_q` of `poly` (ascending coefficients): `gcd` with
/// `X^q − X` (by `m` squarings modulo the polynomial), then splitting by
/// `Tr(δX)`; every multiplication counted.
pub fn roots(f: &Gf2q, poly: &[Fe], muls: &mut u64) -> Vec<Fe> {
    let g = &f.tab;
    let p = pmonic(g, poly.to_vec(), muls);
    if p.len() <= 1 {
        return Vec::new();
    }
    let mut r = pmod(g, &[0, 1], &p, muls);
    for _ in 0..f.m {
        let sq = pmul(g, &r, &r, muls);
        r = pmod(g, &sq, &p, muls);
    }
    if r.len() < 2 {
        r.resize(2, 0);
    }
    r[1] ^= 1;
    ptrim(&mut r);
    let h = if r.is_empty() {
        p.clone()
    } else {
        pgcd(g, &p, &r, muls)
    };
    let mut out = Vec::new();
    split_roots(f, h, &mut out, 1, muls);
    out.sort_unstable();
    out
}

fn split_roots(f: &Gf2q, h: Vec<Fe>, out: &mut Vec<Fe>, mut seed: u32, muls: &mut u64) {
    let g = &f.tab;
    match h.len() {
        0 | 1 => {}
        2 => out.push(h[0]),
        _ => loop {
            let delta = g.exp[(seed % (f.q - 1)) as usize];
            seed += 1;
            let mut x = pmod(g, &[0, delta], &h, muls);
            let mut acc = x.clone();
            for _ in 1..f.m {
                let sq = pmul(g, &x, &x, muls);
                x = pmod(g, &sq, &h, muls);
                let len = acc.len().max(x.len());
                acc.resize(len, 0);
                for (i, &c) in x.iter().enumerate() {
                    acc[i] ^= c;
                }
                ptrim(&mut acc);
            }
            if acc.is_empty() {
                continue;
            }
            let d = pgcd(g, &h, &acc, muls);
            if d.len() > 1 && d.len() < h.len() {
                let e = pdiv(g, &h, &d, muls);
                split_roots(f, d, out, seed, muls);
                split_roots(f, e, out, seed + 7, muls);
                return;
            }
        },
    }
}

// ── The solve ────────────────────────────────────────────────────────────

/// What one solve cost, by phase, in `F_q` multiplications.
#[derive(Clone, Debug, Default, Serialize)]
pub struct BinSolve {
    pub macaulay_degree: u8,
    pub rows: usize,
    pub cols: usize,
    pub pivots: usize,
    pub zero_rows: usize,
    pub dim: usize,
    pub chains: usize,
    pub weil_muls: u64,
    pub triangularise_muls: u64,
    pub echelon_muls: u64,
    pub nf_muls: u64,
    pub charpoly_muls: u64,
    pub roots_muls: u64,
    pub eigenvector_muls: u64,
    pub verify_muls: u64,
    pub split_muls: u64,
    pub total_muls: u64,
    pub rational_eigenvalues: usize,
    pub degenerate_eigenspaces: usize,
    pub e_solutions: usize,
    pub topped_up: u64,
    pub replay: String,
    pub outcome: String,
    pub total_ms: f64,
}

impl BinSolve {
    fn close(&mut self) {
        self.total_muls = self.weil_muls
            + self.triangularise_muls
            + self.echelon_muls
            + self.nf_muls
            + self.charpoly_muls
            + self.roots_muls
            + self.eigenvector_muls
            + self.verify_muls
            + self.split_muls;
    }
}

/// The recorded shape of a solved system (F4Remake's idea).
#[derive(Clone, Debug, Default)]
pub struct ShapeB {
    leads: Vec<Option<Vec<u8>>>,
    rows_built: usize,
    zero_ids: Vec<u32>,
    pivots: usize,
}

/// A shape carried from solve to solve, with its outcomes.
#[derive(Clone, Debug, Default)]
pub struct BinTrace {
    pub shape: Option<ShapeB>,
    pub replays: u64,
    pub fallbacks: u64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct BinOptions {
    /// `(s₂, s₃, s₄, π)`, weights `(2, 2, 2, 1)`, degree 26; otherwise the
    /// unweighted `(s₁, s₂, s₃, s₄, π)` with `π² + π + s₁`, degree 14.
    pub weighted: bool,
    /// The Frobenius-reduced form `(S₂, S₃, s₄, P)`, weights
    /// `(2, 2, 1, 1)`, degree 11 (implies the weighted invariants).
    #[serde(default)]
    pub frobenius: bool,
    pub sparsest: bool,
    pub replay: bool,
}

struct Chains {
    pos: Vec<(u32, u32)>,
    len: Vec<u32>,
    nf: Vec<Vec<(u32, u32, Fe)>>,
}

impl Chains {
    fn powers(g: &GfTab, lam: Fe, top: usize, muls: &mut u64) -> Vec<Fe> {
        let mut pw = vec![1u32; top + 1];
        for t in 1..=top {
            pw[t] = g.mul(pw[t - 1], lam);
        }
        *muls += top as u64;
        pw
    }
    fn top(&self) -> usize {
        *self.len.iter().max().unwrap_or(&0) as usize
    }
    fn matrix(&self, g: &GfTab, lam: Fe, muls: &mut u64) -> Vec<Vec<Fe>> {
        let pw = Self::powers(g, lam, self.top(), muls);
        let k = self.len.len();
        let mut a = vec![vec![0u32; k]; k];
        for (j, row) in self.nf.iter().enumerate() {
            for &(i, t, c) in row {
                a[j][i as usize] ^= g.mul(c, pw[t as usize]);
            }
            *muls += row.len() as u64;
            a[j][j] ^= pw[self.len[j] as usize];
        }
        a
    }
    fn eigenvectors(&self, g: &GfTab, lam: Fe, muls: &mut u64) -> Vec<Vec<Fe>> {
        let a = self.matrix(g, lam, muls);
        let ker = kernel_shifted(g, &a, 0, muls);
        let pw = Self::powers(g, lam, self.top(), muls);
        ker.iter()
            .map(|u| {
                *muls += self.pos.len() as u64;
                self.pos
                    .iter()
                    .map(|&(j, t)| g.mul(u[j as usize], pw[t as usize]))
                    .collect()
            })
            .collect()
    }
}

/// The characteristic polynomial of a small square matrix, by its
/// determinant at `k + 1` points interpolated (needs `q > k`).
fn small_charpoly(f: &Gf2q, a: &[Vec<Fe>], muls: &mut u64) -> Vec<Fe> {
    let g = &f.tab;
    let k = a.len();
    let xs: Vec<Fe> = (0..=k as Fe).collect();
    let vals: Vec<Fe> = xs
        .iter()
        .map(|&x| {
            let mut m = a.to_vec();
            for (i, row) in m.iter_mut().enumerate() {
                row[i] ^= x;
            }
            det(g, m, muls)
        })
        .collect();
    newton_interpolate(g, &xs, &vals, muls)
}

/// The polynomial through `(xs[i], vals[i])`, ascending, by Newton's
/// divided differences (one multiplication per division).
fn newton_interpolate(g: &GfTab, xs: &[Fe], vals: &[Fe], muls: &mut u64) -> Vec<Fe> {
    let n = xs.len();
    let mut dd = vals.to_vec();
    for lvl in 1..n {
        for i in (lvl..n).rev() {
            let num = dd[i] ^ dd[i - 1];
            dd[i] = g.mul(num, g.inv(xs[i] ^ xs[i - lvl]));
        }
        *muls += (n - lvl) as u64;
    }
    let mut poly = vec![0u32; n];
    for i in (0..n).rev() {
        let mut nxt = vec![0u32; n];
        for d in 0..n - 1 {
            if poly[d] != 0 {
                nxt[d + 1] ^= poly[d];
                nxt[d] ^= g.mul(poly[d], xs[i]);
            }
        }
        *muls += (n - 1 - i) as u64;
        nxt[0] ^= dd[i];
        poly = nxt;
    }
    poly
}

/// Split an eigenspace by another unknown `v`, known only on the standard
/// monomials `b` with `v·b` standard (`known`: `(b, v·b)` indices).
fn split_partial(
    f: &Gf2q,
    w: &[Vec<Fe>],
    known: &[(usize, usize)],
    muls: &mut u64,
) -> Vec<Vec<Fe>> {
    let g = &f.tab;
    let k = w.len();
    let mut a = vec![vec![0u32; k]; k];
    for i in 0..k {
        let mut sys: Vec<Vec<Fe>> = known
            .iter()
            .map(|&(l, t)| {
                let mut row: Vec<Fe> = (0..k).map(|j| w[j][l]).collect();
                row.push(w[i][t]);
                row
            })
            .collect();
        let pivots = rref(g, &mut sys, muls);
        if pivots.len() != k || pivots.iter().any(|&c| c >= k) {
            return w.to_vec();
        }
        for (r, &pc) in pivots.iter().enumerate() {
            a[i][pc] = sys[r][k];
        }
    }
    let at: Vec<Vec<Fe>> = (0..k).map(|j| (0..k).map(|i| a[i][j]).collect()).collect();
    let cp = small_charpoly(f, &a, muls);
    let mus = roots(f, &cp, muls);
    let dim = w[0].len();
    let mut out = Vec::new();
    for mu in mus {
        for v in kernel_shifted(g, &at, mu, muls) {
            let mut u = vec![0u32; dim];
            for j in 0..k {
                if v[j] == 0 {
                    continue;
                }
                for r in 0..dim {
                    if w[j][r] != 0 {
                        u[r] ^= g.mul(v[j], w[j][r]);
                        *muls += 1;
                    }
                }
            }
            out.push(u);
        }
    }
    out
}

fn eval<const N: usize>(
    g: &GfTab,
    comp: &HashMap<[u8; N], Fe>,
    e: &[Fe; N],
    muls: &mut u64,
) -> Fe {
    let mut acc = 0u32;
    for (m, &c) in comp {
        let mut t = c;
        for k in 0..N {
            for _ in 0..m[k] {
                t = g.mul(t, e[k]);
                *muls += 1;
            }
        }
        acc ^= t;
    }
    acc
}

/// All `F_q`-rational solutions of `comps` (weights `w`), through a
/// weighted Macaulay matrix at `macaulay`, the eigen step along the
/// chains of unknown `mvar` (weight 1).  `Err(())` only from a replay
/// whose leads or row count differ.
#[allow(clippy::too_many_arguments)]
fn solve_once<const N: usize>(
    f: &Gf2q,
    comps: &[HashMap<[u8; N], Fe>],
    w: &[u8; N],
    mvar: usize,
    macaulay: u8,
    sparsest: bool,
    st: &mut BinSolve,
    replay: Option<&ShapeB>,
    record: Option<&mut ShapeB>,
) -> Result<Option<Vec<[Fe; N]>>, ()> {
    let g = &f.tab;
    let q = f.q as u64;
    let t0 = std::time::Instant::now();
    st.macaulay_degree = macaulay;
    let cols = monomials_desc(macaulay, w);
    let col_index: HashMap<[u8; N], usize> =
        cols.iter().enumerate().map(|(i, m)| (*m, i)).collect();
    let mut tri = 0;
    let (comps, leads) = triangularise(g, comps, w, &mut tri);
    st.triangularise_muls += tri;
    let leads_v: Vec<Option<Vec<u8>>> = leads.iter().map(|l| l.map(|m| m.to_vec())).collect();
    if let Some(t) = replay {
        if t.leads != leads_v {
            return Err(());
        }
    }
    let divides = |lm: &[u8; N], m: &[u8; N]| (0..N).all(|k| lm[k] <= m[k]);
    let mut rows: Vec<Vec<Fe>> = Vec::new();
    let mut held: Vec<Vec<Fe>> = Vec::new();
    let (mut built, mut skip_at) = (0usize, 0usize);
    for (i, comp) in comps.iter().enumerate() {
        let dg = comp.keys().map(|m| wdeg(m, w)).max().unwrap_or(0) as u8;
        if dg > macaulay {
            st.outcome = "an equation above the Macaulay degree".into();
            return Ok(None);
        }
        for sh in monomials_desc::<N>(macaulay - dg, w) {
            if leads
                .iter()
                .take(i)
                .any(|lm| lm.is_some_and(|lm| divides(&lm, &sh)))
            {
                continue;
            }
            let index = built;
            built += 1;
            let mut row = vec![0u32; cols.len()];
            for (&e, &c) in comp {
                let m: [u8; N] = core::array::from_fn(|k| e[k] + sh[k]);
                row[col_index[&m]] ^= c;
            }
            if let Some(t) = replay {
                if skip_at < t.zero_ids.len() && t.zero_ids[skip_at] as usize == index {
                    skip_at += 1;
                    held.push(row);
                    continue;
                }
            }
            rows.push(row);
        }
    }
    if let Some(t) = replay {
        if t.rows_built != built {
            return Err(());
        }
    }
    st.rows = rows.len();
    st.cols = cols.len();
    let mut em = 0;
    let mut ech = echelon(g, rows, &mut em, sparsest);
    if let Some(t) = replay {
        let have = ech.pivot_row.iter().filter(|&&r| r != u32::MAX).count();
        if have < t.pivots {
            st.topped_up += top_up(g, &mut ech, held, t.pivots - have, &mut em);
        }
    }
    st.echelon_muls += em;
    st.zero_rows = ech.zero_rows;
    st.pivots = ech.pivot_row.iter().filter(|&&r| r != u32::MAX).count();
    let is_std: Vec<bool> = (0..cols.len())
        .map(|c| ech.pivot_row[c] == u32::MAX && wdeg(&cols[c], w) < macaulay as Fe)
        .collect();
    let standard: Vec<usize> = (0..cols.len()).filter(|&c| is_std[c]).collect();
    let dim = standard.len();
    st.dim = dim;
    if let Some(out) = record {
        *out = ShapeB {
            leads: leads_v,
            rows_built: built,
            zero_ids: ech.zero_ids.clone(),
            pivots: st.pivots,
        };
    }
    let mut std_index = vec![u32::MAX; cols.len()];
    for (i, &c) in standard.iter().enumerate() {
        std_index[c] = i as u32;
    }
    // Columns of 1 and of each unknown.
    let mut idx = Vec::with_capacity(N + 1);
    for v in 0..=N {
        let m: [u8; N] = core::array::from_fn(|k| u8::from(v > 0 && k == v - 1));
        match col_index.get(&m) {
            Some(&c) if is_std[c] => idx.push(std_index[c] as usize),
            _ => {
                st.outcome = "1 or a variable is not standard".into();
                return Ok(None);
            }
        }
    }
    // Chains of the multiplication unknown.
    let mut next: Vec<Result<usize, usize>> = Vec::with_capacity(dim);
    for &bc in &standard {
        let mut m = cols[bc];
        m[mvar] += 1;
        let Some(&t) = col_index.get(&m) else {
            st.outcome = "border unreachable at this degree".into();
            return Ok(None);
        };
        next.push(if is_std[t] {
            Ok(std_index[t] as usize)
        } else {
            Err(t)
        });
    }
    let mut pos = vec![(u32::MAX, u32::MAX); dim];
    let (mut len, mut ends) = (Vec::new(), Vec::new());
    for s in 0..dim {
        if cols[standard[s]][mvar] != 0 {
            continue;
        }
        let j = len.len() as u32;
        let (mut cur, mut t) = (s, 0u32);
        loop {
            if pos[cur].0 != u32::MAX {
                st.outcome = "chain structure: a monomial on two chains".into();
                return Ok(None);
            }
            pos[cur] = (j, t);
            t += 1;
            match next[cur] {
                Ok(n) => cur = n,
                Err(b) => {
                    ends.push(b);
                    break;
                }
            }
        }
        len.push(t);
    }
    if pos.iter().any(|&(j, _)| j == u32::MAX) {
        st.outcome = "chain structure: standard monomials are not an order ideal".into();
        return Ok(None);
    }
    st.chains = len.len();
    let ncols = cols.len();
    let nnz: Vec<u64> = (0..ncols)
        .map(|c| match ech.pivot_row[c] {
            u32::MAX => 0,
            r => ech.rows[r as usize][c..]
                .iter()
                .filter(|&&v| v != 0)
                .count() as u64,
        })
        .collect();
    let reduced: Vec<Result<(Vec<Fe>, u64), ()>> = ends
        .par_iter()
        .map(|&b| {
            let mut v = vec![0u32; ncols];
            v[b] = 1;
            let mut m = 0u64;
            for c in b..ncols {
                if v[c] == 0 || is_std[c] {
                    continue;
                }
                let pr = ech.pivot_row[c];
                if pr == u32::MAX {
                    return Err(());
                }
                let fct = v[c];
                g.axpy(&mut v[c..], &ech.rows[pr as usize][c..], fct);
                m += nnz[c];
            }
            Ok((v, m))
        })
        .collect();
    let mut nf = Vec::with_capacity(ends.len());
    for r in reduced {
        let Ok((v, m)) = r else {
            st.outcome = "border unreachable at this degree".into();
            return Ok(None);
        };
        st.nf_muls += m;
        nf.push(
            standard
                .iter()
                .enumerate()
                .filter(|&(_, &c)| v[c] != 0)
                .map(|(s, &c)| (pos[s].0, pos[s].1, v[c]))
                .collect(),
        );
    }
    let ch = Chains { pos, len, nf };
    // Rational eigenvalues: zeros of det A on F_q, or its roots.
    let npts = if q <= dim as u64 + 1 {
        q
    } else {
        dim as u64 + 1
    };
    let vals: Vec<(Fe, u64)> = (0..npts as Fe)
        .into_par_iter()
        .map(|x| {
            let mut m = 0;
            (det(g, ch.matrix(g, x, &mut m), &mut m), m)
        })
        .collect();
    st.charpoly_muls += vals.iter().map(|v| v.1).sum::<u64>();
    let lambdas: Vec<Fe> = if npts == q {
        (0..q as Fe).filter(|&x| vals[x as usize].0 == 0).collect()
    } else {
        let xs: Vec<Fe> = (0..npts as Fe).collect();
        let ys: Vec<Fe> = vals.iter().map(|v| v.0).collect();
        let mut cm = 0;
        let poly = newton_interpolate(g, &xs, &ys, &mut cm);
        st.charpoly_muls += cm;
        if poly.len() != dim + 1 || poly[dim] != 1 {
            st.outcome = "det A is not monic of degree dim".into();
            return Ok(None);
        }
        let mut rm = 0;
        let r = roots(f, &poly, &mut rm);
        st.roots_muls += rm;
        r
    };
    st.rational_eigenvalues = lambdas.len();
    let known: Vec<Vec<(usize, usize)>> = (0..N)
        .filter(|&v| v != mvar)
        .map(|v| {
            (0..dim)
                .filter_map(|s| {
                    let mut m = cols[standard[s]];
                    m[v] += 1;
                    let &t = col_index.get(&m)?;
                    is_std[t].then(|| (s, std_index[t] as usize))
                })
                .collect()
        })
        .collect();
    let mut out = Vec::new();
    for lam in lambdas {
        let mut ev = 0;
        let mut ker = ch.eigenvectors(g, lam, &mut ev);
        if ker.len() > 1 {
            st.degenerate_eigenspaces += 1;
            for kn in &known {
                if ker.len() <= 1 {
                    break;
                }
                ker = split_partial(f, &ker, kn, &mut ev);
            }
        }
        st.eigenvector_muls += ev;
        for v in &ker {
            if v[idx[0]] == 0 {
                continue;
            }
            let inv = g.inv(v[idx[0]]);
            let e: [Fe; N] = core::array::from_fn(|k| g.mul(v[idx[k + 1]], inv));
            if e[mvar] != lam {
                continue;
            }
            let mut vm = 0;
            let ok = comps.iter().all(|c| eval(g, c, &e, &mut vm) == 0);
            st.verify_muls += vm;
            if ok {
                out.push(e);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    st.e_solutions = out.len();
    st.outcome = "solved".into();
    st.total_ms += t0.elapsed().as_secs_f64() * 1e3;
    Ok(Some(out))
}

/// [`solve_once`] with a shape trace: replay when one is recorded (a
/// mismatch of leads or rows runs the full solve too, both charged),
/// record on the first solved system.
#[allow(clippy::too_many_arguments)]
pub fn solve_bin<const N: usize>(
    f: &Gf2q,
    comps: &[HashMap<[u8; N], Fe>],
    w: &[u8; N],
    mvar: usize,
    macaulay: u8,
    sparsest: bool,
    st: &mut BinSolve,
    trace: Option<&mut BinTrace>,
) -> Option<Vec<[Fe; N]>> {
    let Some(tr) = trace else {
        st.replay = "none".into();
        return solve_once(f, comps, w, mvar, macaulay, sparsest, st, None, None)
            .ok()
            .flatten();
    };
    if let Some(shape) = &tr.shape {
        if let Ok(s) = solve_once(f, comps, w, mvar, macaulay, sparsest, st, Some(shape), None) {
            tr.replays += 1;
            st.replay = "replayed".into();
            return s;
        }
        tr.fallbacks += 1;
        st.replay = "fell back".into();
        return solve_once(f, comps, w, mvar, macaulay, sparsest, st, None, None)
            .ok()
            .flatten();
    }
    let mut shape = ShapeB::default();
    let s = solve_once(
        f,
        comps,
        w,
        mvar,
        macaulay,
        sparsest,
        st,
        None,
        Some(&mut shape),
    )
    .ok()
    .flatten();
    st.replay = "recorded".into();
    if s.is_some() {
        tr.shape = Some(shape);
    }
    s
}

// ── One decomposition ────────────────────────────────────────────────────

#[derive(Clone, Debug, Default, Serialize)]
pub struct BinD4Solve {
    pub solve: BinSolve,
    pub rational_solutions: usize,
    /// The decompositions as sorted `t`-quadruples.
    pub decompositions: Vec<[Fe; 4]>,
    #[serde(skip)]
    pub witness: Option<[PtB; 4]>,
    /// Every decomposition's signed points, in the order found (the
    /// large-prime collector picks among them).
    #[serde(skip)]
    pub witnesses: Vec<[PtB; 4]>,
}

/// Decompose `r` as `±P₁ ± P₂ ± P₃ ± P₄` with `t(P_i) ∈ F_q`: the Weil
/// restriction of `H(·; t(r))`, its solve, and for each rational solution
/// the split `w_i` (roots of `X⁴ + s₁X³ + s₂X² + s₃X + s₄`), `t_i² + t_i =
/// w_i` with `Σ t_i = π`, and the check that the points sum to `r`.
/// `Err` carries what an unsolved attempt cost.
pub fn solve_d4_bin(
    curve: &CurveB,
    pre: &SymmetrisedBin,
    r: &PtB,
    opts: BinOptions,
    macaulay: u8,
    rng: &mut StdRng,
    trace: Option<&mut BinTrace>,
) -> Result<BinD4Solve, BinSolve> {
    let f = &curve.f;
    let g = &f.f.tab;
    let mut out = BinD4Solve::default();
    let before = f.f.muls();
    let t_r = curve.t_of_x(&r.x);
    let sols: Option<Vec<[Fe; 5]>> = if opts.frobenius {
        let comps = pre.system_frobenius(f, &t_r);
        out.solve.weil_muls = f.f.muls() - before;
        let fq = &f.f;
        solve_bin::<4>(
            fq,
            &comps,
            &FROBENIUS_WEIGHTS,
            3,
            macaulay,
            opts.sparsest,
            &mut out.solve,
            trace,
        )
        .map(|v| {
            v.into_iter()
                .map(|s| {
                    // s₂ = √S₂, s₃ = √S₃, π = √P, s₁ = P + π.
                    let before = fq.muls();
                    let (s2, s3, pi) = (fq.sqrt(s[0]), fq.sqrt(s[1]), fq.sqrt(s[3]));
                    out.solve.split_muls += fq.muls() - before;
                    [s[3] ^ pi, s2, s3, s[2], pi]
                })
                .collect()
        })
    } else if opts.weighted {
        let comps = pre.system_weighted(f, &t_r);
        out.solve.weil_muls = f.f.muls() - before;
        solve_bin::<4>(
            &f.f,
            &comps,
            &WEIGHTS,
            3,
            macaulay,
            opts.sparsest,
            &mut out.solve,
            trace,
        )
        .map(|v| {
            v.into_iter()
                .map(|s| {
                    out.solve.split_muls += 1;
                    [g.mul(s[3], s[3]) ^ s[3], s[0], s[1], s[2], s[3]]
                })
                .collect()
        })
    } else {
        let comps = pre.system_plain(f, &t_r);
        out.solve.weil_muls = f.f.muls() - before;
        solve_bin::<5>(
            &f.f,
            &comps,
            &[1; 5],
            4,
            macaulay,
            opts.sparsest,
            &mut out.solve,
            trace,
        )
    };
    let Some(sols) = sols else {
        out.solve.close();
        return Err(out.solve);
    };
    out.rational_solutions = sols.len();
    let mut found = BTreeSet::new();
    for s in sols {
        // s = (s₁, s₂, s₃, s₄, π).
        let mut sm = 0;
        let quartic = [s[3], s[2], s[1], s[0], 1];
        let rts = roots(&f.f, &quartic, &mut sm);
        out.solve.split_muls += sm;
        // The multisets of roots with these symmetric functions.
        let mut multisets = Vec::new();
        for (i1, &a) in rts.iter().enumerate() {
            for (i2, &b) in rts.iter().enumerate().skip(i1) {
                for (i3, &c) in rts.iter().enumerate().skip(i2) {
                    for &d in rts.iter().skip(i3) {
                        let x = [a, b, c, d];
                        let mut e = [1u32, 0, 0, 0, 0];
                        for &v in &x {
                            for k in (1..5).rev() {
                                e[k] ^= g.mul(e[k - 1], v);
                            }
                        }
                        out.solve.split_muls += 16;
                        if e[1] == s[0] && e[2] == s[1] && e[3] == s[2] && e[4] == s[3] {
                            multisets.push(x);
                        }
                    }
                }
            }
        }
        for ws in multisets {
            let base_t: Option<Vec<Fe>> = ws.iter().map(|&w| f.f.artin_schreier(w)).collect();
            let Some(base_t) = base_t else {
                continue;
            };
            for flips in 0..16u32 {
                let ts: [Fe; 4] = core::array::from_fn(|i| base_t[i] ^ ((flips >> i & 1) as Fe));
                if ts.iter().fold(0u32, |a, &v| a ^ v) != s[4] {
                    continue;
                }
                if ts.iter().any(|&v| v < 2) {
                    continue;
                }
                let mut sorted = ts;
                sorted.sort_unstable();
                if sorted.windows(2).any(|w| w[0] == w[1]) {
                    continue;
                }
                // The lifts' field work is the split's (the sign checks
                // below are charged as additions).
                let lift_before = f.f.muls();
                let pts: Option<Vec<PtB>> = ts
                    .iter()
                    .map(|&v| curve.lift_x(&curve.x_of_t(&f.from_fq(v)), rng))
                    .collect();
                out.solve.split_muls += f.f.muls() - lift_before;
                let Some(pts) = pts else {
                    continue;
                };
                let pts = [pts[0], pts[1], pts[2], pts[3]];
                if let Some(signed) = signed_sum(curve, &pts, r) {
                    if found.insert(sorted) {
                        out.witnesses.push(signed);
                    }
                    if out.witness.is_none() {
                        out.witness = Some(signed);
                    }
                }
            }
        }
    }
    out.decompositions = found.into_iter().collect();
    out.solve.close();
    Ok(out)
}

// ── Curves with a prime-order subgroup, for ecbench ──────────────────────

/// `#E(F_{q⁴})`: the one multiple in the Hasse interval
/// `[q⁴ + 1 − 2q², q⁴ + 1 + 2q²]` that kills a few random points, each
/// point's killing multiples found by baby steps and giant steps
/// (`≈ 2q` of each over the interval's `4q² + 1` values).
pub fn group_order_bin(c: &CurveB, rng: &mut StdRng) -> Option<u128> {
    let q = c.f.f.q as u128;
    let q4 = q.pow(4);
    let (lo, hi) = (q4 + 1 - 2 * q * q, q4 + 1 + 2 * q * q);
    let width = hi - lo + 1;
    let b = (width as f64).sqrt().ceil() as u128;
    let key = |p: &PtB| -> u128 {
        if p.inf {
            u128::MAX
        } else {
            c.f.pack(&p.x) * 2 + u128::from(c.f.pack(&p.y) > c.f.pack(&p.x.add(&p.y)))
        }
    };
    let mut cands: Option<Vec<u128>> = None;
    for _ in 0..24 {
        let pt = loop {
            let x = c.f.random(rng);
            if let Some(p) = c.lift_x(&x, rng) {
                break p;
            }
        };
        // Baby steps: j·P for 0 ≤ j < b.
        let mut baby: HashMap<u128, Vec<u128>> = HashMap::new();
        let mut acc = PtB::INF;
        for j in 0..b {
            baby.entry(key(&acc)).or_default().push(j);
            acc = c.add(&acc, &pt);
        }
        // Giant steps: (lo + i·b)·P; mP = O for m = lo + i·b + j iff
        // (lo + i·b)·P = −j·P.
        let step = c.mul(&pt, b);
        let mut g = c.mul(&pt, lo);
        let mut here = Vec::new();
        let mut i = 0u128;
        while lo + i * b <= hi {
            let target = key(&c.neg(&g));
            if let Some(js) = baby.get(&target) {
                for &j in js {
                    let m = lo + i * b + j;
                    if m <= hi && c.mul(&pt, m).inf {
                        here.push(m);
                    }
                }
            }
            g = c.add(&g, &step);
            i += 1;
        }
        here.sort_unstable();
        here.dedup();
        cands = Some(match cands {
            None => here,
            Some(v) => v.into_iter().filter(|m| here.contains(m)).collect(),
        });
        if cands.as_ref().is_some_and(|v| v.len() == 1) {
            return cands.map(|v| v[0]);
        }
    }
    None
}

/// Deterministic Miller–Rabin: the bases 2, 3, …, 37 decide every
/// `n < 3.3·10²⁴` (above `2^81`), which covers `#E(F_{q⁴})` for `m ≤ 20`.
/// Below `2^64` it is the one-word test every committed curve used.
fn is_prime(n: u128) -> bool {
    if let Ok(n) = u64::try_from(n) {
        return super::residual_walk::is_prime_u64(n);
    }
    use super::koblitz_strong_rho::RhoScalar;
    let mulm = |a: u128, b: u128| <u128 as RhoScalar>::mul_mod(a, b, n);
    let powm = |mut b: u128, mut e: u128| {
        let mut acc = 1u128;
        while e > 0 {
            if e & 1 == 1 {
                acc = mulm(acc, b);
            }
            b = mulm(b, b);
            e >>= 1;
        }
        acc
    };
    if n % 2 == 0 {
        return false;
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'bases: for a in [2u128, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powm(a % n, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulm(x, x);
            if x == n - 1 {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

/// A curve of [`generate_instance`]'s family with `#E = h·r`, `r` prime.
pub struct BinCurve {
    pub curve: CurveB,
    pub instance_seed: u64,
    /// Orders as `u128`: `#E(F_{q⁴}) ≈ 2^{4m}` leaves one word at `m = 16`.
    pub group_order: u128,
    pub r: u128,
    pub cofactor: u128,
    pub generator: PtB,
}

/// The first instance seed from `seed` with `#E = h·r`, `r` prime,
/// `h ≤ max_cofactor`.  Deterministic in `(m, seed, max_cofactor)`.
pub fn find_bin_curve(m: u32, seed: u64, max_cofactor: u64) -> Result<BinCurve, String> {
    let max_cofactor = u128::from(max_cofactor);
    for s in seed..seed + 2_000 {
        let curve = generate_instance(m, s);
        let mut rng = StdRng::seed_from_u64(s ^ 0x0c0f);
        let Some(order) = group_order_bin(&curve, &mut rng) else {
            continue;
        };
        let Some(h) = (2..=max_cofactor).find(|&h| order % h == 0 && is_prime(order / h)) else {
            continue;
        };
        let r = order / h;
        let generator = loop {
            let x = curve.f.random(&mut rng);
            let Some(p) = curve.lift_x(&x, &mut rng) else {
                continue;
            };
            let g = curve.mul(&p, h);
            if !g.inf {
                break g;
            }
        };
        if !curve.mul(&generator, r).inf {
            continue;
        }
        curve.f.f.muls.set(0);
        return Ok(BinCurve {
            curve,
            instance_seed: s,
            group_order: order,
            r,
            cofactor: h,
            generator,
        });
    }
    Err(format!(
        "no curve with cofactor ≤ {max_cofactor} over F_(2^{m})^4 from seed {seed}"
    ))
}

/// `F_q` multiplications in one affine addition of two distinct points.
pub fn fq_muls_per_add(c: &CurveB) -> u64 {
    let mut rng = StdRng::seed_from_u64(0xadd);
    let p = loop {
        let x = c.f.random(&mut rng);
        if let Some(p) = c.lift_x(&x, &mut rng) {
            break p;
        }
    };
    let q = c.add(&p, &p);
    let before = c.f.f.muls();
    let _ = c.add(&p, &q);
    c.f.f.muls() - before
}

// ── In ecbench: the counted group, the identity's basis, the base, the oracle

impl CountedGroup for CurveB {
    type Elt = PtB;
    fn identity(&self) -> PtB {
        PtB::INF
    }
    fn is_identity(&self, p: &PtB) -> bool {
        p.inf
    }
    fn add(&self, ops: &mut GroupOps, p: PtB, q: PtB) -> PtB {
        if p.inf || q.inf || p.x != q.x {
            ops.adds += 1;
        } else {
            ops.doubles += 1;
        }
        CurveB::add(self, &p, &q)
    }
    fn double(&self, ops: &mut GroupOps, p: PtB) -> PtB {
        ops.doubles += 1;
        CurveB::add(self, &p, &p)
    }
    fn neg(&self, p: PtB) -> PtB {
        CurveB::neg(self, &p)
    }
    /// The packed abscissa and one bit for the ordinate (`y` or `x + y`, the
    /// larger packed value): exact for `4m ≤ 62`, the key every committed
    /// session used.  Past that ([`wide_key`]) it is a 64-bit mix of the
    /// exact `4m + 1`-bit key.
    fn key(&self, p: &PtB) -> u64 {
        if p.inf {
            return u64::MAX;
        }
        let y = self.f.pack(&p.y);
        let ny = self.f.pack(&p.x.add(&p.y));
        let exact = self.f.pack(&p.x) * 2 + u128::from(y > ny);
        if 4 * self.f.f.m < 63 {
            exact as u64
        } else {
            wide_key(exact)
        }
    }
}

/// A 64-bit key for an exact key wider than a word (`m ≥ 16`): two rounds
/// of the splitmix64 finaliser over the two halves.  Distinct points can
/// share a key with probability about `2^{-64}` per pair; every use is
/// followed by an exact check (a factor-base index is checked by the
/// group law, a rho collision by `[d]G = Q`), so a shared key costs a
/// wasted step, never a wrong answer.
pub fn wide_key(exact: u128) -> u64 {
    let mix = |mut z: u64| {
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    };
    let k = mix(mix(exact as u64) ^ (exact >> 64) as u64);
    // u64::MAX is the identity's key.
    if k == u64::MAX {
        k - 1
    } else {
        k
    }
}

/// The tower as `GF(2)[X]/(g)`: `θ`, the first of `z, z + 1, z + x, …`
/// (`x` the generator of `F_q`) whose minimal polynomial over `F_2` has
/// degree `4m`, gives `g` (bit `i` for `X^i`, the `X^{4m}` bit included)
/// and the coordinates of an element in `1, θ, …, θ^{4m−1}`.  The curve's
/// identity is `curve_id::binary` of `(g, a, b)` in that basis.
pub struct TowerBasis {
    pub modulus: u128,
    /// Rows of the inverse of the matrix whose column `i` is `θ^i` packed.
    inv: Vec<u128>,
}

impl TowerBasis {
    pub fn new(f: &Fq4) -> TowerBasis {
        let n = 4 * f.f.m as usize;
        for c0 in 0..f.f.q as Fe {
            for c1 in 1..f.f.q as Fe {
                let theta = Eb([c0, c1, 0, 0]);
                // Powers θ⁰ … θⁿ packed; reduce θⁿ by the first n.
                let mut pw = Vec::with_capacity(n + 1);
                let mut acc = Eb::ONE;
                for _ in 0..=n {
                    pw.push(f.pack(&acc));
                    acc = f.mul_nc(&acc, &theta);
                }
                // Gaussian elimination over F_2 on columns: express each
                // packed bit vector in the θ-power basis.
                let mut rows: Vec<(u128, u128)> = Vec::new(); // (vector, combination)
                let mut ok = true;
                for (i, &v) in pw.iter().take(n).enumerate() {
                    let (mut v, mut comb) = (v, 1u128 << i);
                    for &(bv, bc) in &rows {
                        let lead = 127 - bv.leading_zeros();
                        if v >> lead & 1 == 1 {
                            v ^= bv;
                            comb ^= bc;
                        }
                    }
                    if v == 0 {
                        ok = false;
                        break;
                    }
                    rows.push((v, comb));
                    rows.sort_by_key(|x| std::cmp::Reverse(127 - x.0.leading_zeros()));
                }
                if !ok {
                    continue;
                }
                let solve = |mut v: u128| -> u128 {
                    let mut comb = 0u128;
                    for &(bv, bc) in &rows {
                        let lead = 127 - bv.leading_zeros();
                        if v >> lead & 1 == 1 {
                            v ^= bv;
                            comb ^= bc;
                        }
                    }
                    debug_assert_eq!(v, 0);
                    comb
                };
                let modulus = (1u128 << n) | solve(pw[n]);
                // The inverse map: packed bit j ↦ its θ-coordinates.
                let inv: Vec<u128> = (0..n).map(|j| solve(1u128 << j)).collect();
                return TowerBasis { modulus, inv };
            }
        }
        unreachable!("the tower has a generator over F_2")
    }
    /// `e` in the basis `1, θ, …, θ^{4m−1}`.
    pub fn coords(&self, f: &Fq4, e: &Eb) -> u128 {
        let v = f.pack(e);
        (0..self.inv.len())
            .filter(|&j| v >> j & 1 == 1)
            .fold(0u128, |acc, j| acc ^ self.inv[j])
    }
}

/// The coefficient of base point `i` of [`TLineBase`] modulo `r`: the
/// points of a column come as `P, −P, P + T, −(P + T)`.
pub fn t_line_coef(i: usize, r: u128) -> u128 {
    if i % 2 == 0 {
        1
    } else {
        r - 1
    }
}

/// **`bin-t-line`**: `{P : t(P) ∈ F_q \ {0, 1}}`, a column per pair
/// `{t, t + 1}`: the points `±P, ±(P + T)`, `P + T` computed (one
/// addition, charged).  `[h](P + T) = [h]P` since `h` is even.
pub struct TLineBase<'a> {
    pub curve: &'a BinCurve,
}

impl FactorBaseBuilder<CurveB> for TLineBase<'_> {
    fn name(&self) -> &str {
        "bin-t-line"
    }
    fn describe(&self, _params: &Params) -> String {
        "FGHR's t-line in characteristic 2: t(P) ∈ F_q, columns {±P, ±(P + T)}".into()
    }
    fn build(
        &self,
        _ctx: &InstanceCtx<CurveB>,
        _params: &Params,
        ops: &mut GroupOps,
    ) -> Result<FactorBase<PtB>, String> {
        let d = self.curve;
        let c = &d.curve;
        let t = c.torsion();
        let mut rng = StdRng::seed_from_u64(0x9b);
        let (mut points, mut col_of, mut coef_of) = (Vec::new(), Vec::new(), Vec::new());
        let (mut columns, mut lifts) = (0usize, 0u64);
        for tv in (2..c.f.f.q as Fe).step_by(2) {
            let x = c.x_of_t(&c.f.from_fq(tv));
            lifts += 1;
            let Some(pt) = c.lift_x(&x, &mut rng) else {
                continue;
            };
            let pt = if CountedGroup::key(c, &pt) & 1 == 1 {
                pt
            } else {
                CurveB::neg(c, &pt)
            };
            let pt_t = CountedGroup::add(c, ops, pt, t);
            // −1 mod r; a wide `r` (`m ≥ 16`) has no one-word image, and
            // its pipeline reads the sign from the index ([`t_line_coef`]).
            let minus = u64::try_from(d.r - 1).unwrap_or(0);
            for (q, coef) in [
                (pt, 1),
                (CurveB::neg(c, &pt), minus),
                (pt_t, 1),
                (CurveB::neg(c, &pt_t), minus),
            ] {
                points.push(q);
                col_of.push(columns);
                coef_of.push(coef);
            }
            columns += 1;
        }
        if columns == 0 {
            return Err("no t in F_q carries a point".into());
        }
        let m = c.f.f.m;
        let mut fb = FactorBase::from_column_map(
            format!("t-line over F_(2^{m})^4, {columns} columns {{±P, ±(P + T)}}"),
            points,
            col_of,
            coef_of,
            columns,
            |q| CountedGroup::key(c, q),
            |q| CountedGroup::key(c, &CurveB::neg(c, q)),
            |q| {
                let v = c.f.pack(&q.x);
                if 4 * c.f.f.m < 63 {
                    v as u64
                } else {
                    wide_key(v)
                }
            },
        )?;
        fb.cost.count("sqrt_solves", lifts);
        Ok(fb)
    }
}

/// **`gaudry-d4-bin`**: the four-summand decomposition by [`solve_d4_bin`].
/// Its `F_q` multiplications are solver work in the unit `F_q
/// multiplications`, with the curve's own `fq_muls_per_add` in `extra` so
/// the runner prices them by count; the additions that check a
/// decomposition are charged as additions.  Unsolved systems are charged
/// when `charge_unsolved`.
pub struct BinD4Oracle<'a> {
    pub curve: &'a BinCurve,
    pub macaulay: u8,
    pub seed: u64,
    pub opts: BinOptions,
    pub charge_unsolved: bool,
    trace: BinTrace,
    pre: Option<SymmetrisedBin>,
    totals: SolverTotals,
    per_add: u64,
    /// Large primes (`large_primes`, `small_permille`); off by default.
    pub lp: LargePrimes,
    last_witnesses: Vec<[PtB; 4]>,
}

/// **One and two large primes** for `gaudry-d4-bin` (Gaudry–Thomé–
/// Thériault–Diem's variant, with Thériault's single large prime as its
/// `max_large = 1` case).
///
/// The base is cut to its first `small_permille`‰ columns `F′` (the
/// t-line's columns are in increasing `t`, so this is an arbitrary,
/// fixed subset); the other columns are *large primes*. A decomposition
/// with every column in `F′` is a full relation. One with one or two
/// large columns (up to `max_large`) is a partial and is eliminated
/// against stored partials modulo `r`: each stored partial is the pivot
/// of its largest large column, and reducing a new partial by the pivot
/// of its largest large column replaces that column by smaller ones, so
/// the loop ends, either with no large column left (a full relation,
/// returned as a [`Relation::Row`] over `F′`) or as a new pivot. This is
/// the graph method's cycle-closing step written as elimination; a
/// two-large-prime partial is an edge, a one-large-prime partial an edge
/// to the root. Every combination is a linear combination of group
/// equations, so a returned row holds as one.
///
/// Of a target's decompositions the one with fewest large columns is
/// used; the solver enumerates them all anyway, so the choice is free.
/// A merge costs one multiplication modulo `r` per entry of the pivot
/// it subtracts, plus an inversion; those are counted
/// (`lp_merge_mulmods`) and charged like solver multiplications.
#[derive(Clone, Debug, Default)]
pub struct LargePrimes {
    pub max_large: u8,
    pub small_permille: u64,
    small_cols: usize,
    pivots: HashMap<usize, LpRow>,
    pub counts: BTreeMap<String, u64>,
}

#[derive(Clone, Debug)]
struct LpRow {
    small: Vec<(usize, u128)>,
    large: Vec<(usize, u128)>,
    a: u128,
    b: u128,
}

/// A full relation over `F′`: `h·a + h·b·d = Σ coef · x_col`.
#[derive(Clone, Debug)]
pub struct LpFull {
    pub columns: Vec<(usize, u128)>,
    pub a: u128,
    pub b: u128,
}

/// Arithmetic modulo `r < 2^127`; exact at every width, so a one-word
/// `r` gives the values a one-word implementation would.
fn lp_mul(a: u128, b: u128, r: u128) -> u128 {
    use super::koblitz_strong_rho::RhoScalar;
    <u128 as RhoScalar>::mul_mod(a, b, r)
}

fn lp_inv(a: u128, r: u128) -> u128 {
    // r is prime: a^(r−2).
    let (mut base, mut e, mut acc) = (a % r, r - 2, 1u128);
    while e > 0 {
        if e & 1 == 1 {
            acc = lp_mul(acc, base, r);
        }
        base = lp_mul(base, base, r);
        e >>= 1;
    }
    acc
}

fn lp_sub(a: u128, b: u128, r: u128) -> u128 {
    if a >= b {
        a - b
    } else {
        r - (b - a)
    }
}

fn lp_add(a: u128, b: u128, r: u128) -> u128 {
    lp_sub(a, r - b, r)
}

/// `x − f·y` on sorted sparse vectors modulo `r`, dropping zeros.
fn lp_axpy(x: &[(usize, u128)], f: u128, y: &[(usize, u128)], r: u128) -> Vec<(usize, u128)> {
    let mut out = Vec::with_capacity(x.len() + y.len());
    let (mut i, mut j) = (0, 0);
    while i < x.len() || j < y.len() {
        let take_x = j >= y.len() || (i < x.len() && x[i].0 < y[j].0);
        let take_y = i >= x.len() || (j < y.len() && y[j].0 < x[i].0);
        if take_x {
            out.push(x[i]);
            i += 1;
        } else if take_y {
            let v = lp_sub(0, lp_mul(f, y[j].1, r), r);
            if v != 0 {
                out.push((y[j].0, v));
            }
            j += 1;
        } else {
            let fy = lp_mul(f, y[j].1, r);
            let v = if x[i].1 >= fy { x[i].1 - fy } else { r - (fy - x[i].1) };
            if v != 0 {
                out.push((x[i].0, v));
            }
            i += 1;
            j += 1;
        }
    }
    out
}

impl LargePrimes {
    pub fn enabled(&self) -> bool {
        self.max_large > 0 || (self.small_permille > 0 && self.small_permille < 1000)
    }

    fn bump(&mut self, k: &str, v: u64) {
        *self.counts.entry(k.into()).or_insert(0) += v;
    }

    /// Feed one decomposition (columns and coefficients); a full
    /// relation comes back when one closes. Returns the multiplications
    /// modulo `r` it spent.
    pub fn feed(
        &mut self,
        cols: &[(usize, u128)],
        a: u128,
        b: u128,
        r: u128,
    ) -> (Option<LpFull>, u64) {
        let mut merged: BTreeMap<usize, u128> = BTreeMap::new();
        for &(c, v) in cols {
            let e = merged.entry(c).or_insert(0);
            *e = lp_add(*e, v % r, r);
        }
        let (mut small, mut large) = (Vec::new(), Vec::new());
        for (c, v) in merged {
            if v == 0 {
                continue;
            }
            if c < self.small_cols {
                small.push((c, v));
            } else {
                large.push((c, v));
            }
        }
        if large.len() > self.max_large as usize {
            self.bump("lp_discarded", 1);
            return (None, 0);
        }
        self.bump(&format!("lp_decompositions_{}_large", large.len()), 1);
        let mut row = LpRow { small, large, a, b };
        let mut ops = 0u64;
        let mut merges = 0u64;
        loop {
            let Some(&(lc, lv)) = row.large.last() else {
                self.bump(
                    if merges == 0 {
                        "lp_full_direct"
                    } else {
                        "lp_full_merged"
                    },
                    1,
                );
                self.bump("lp_full_weight", row.small.len() as u64);
                self.bump("lp_merges", merges);
                return (
                    Some(LpFull {
                        columns: row.small,
                        a: row.a,
                        b: row.b,
                    }),
                    ops,
                );
            };
            let Some(p) = self.pivots.get(&lc) else {
                self.bump(&format!("lp_stored_{}_large", row.large.len()), 1);
                self.pivots.insert(lc, row);
                self.bump("lp_merges", merges);
                return (None, ops);
            };
            let pv = p.large.last().expect("a pivot has its column").1;
            let f = lp_mul(lv, lp_inv(pv, r), r);
            ops += 64 + p.small.len() as u64 + p.large.len() as u64 + 2;
            row.small = lp_axpy(&row.small, f, &p.small, r);
            row.large = lp_axpy(&row.large, f, &p.large, r);
            row.a = lp_sub(row.a, lp_mul(f, p.a, r), r);
            row.b = lp_sub(row.b, lp_mul(f, p.b, r), r);
            merges += 1;
        }
    }
}

/// A relation from [`BinD4Oracle::relations_for`], at any width.
pub enum BinRelation {
    Summands(Vec<usize>),
    Row(LpFull),
}

impl<'a> BinD4Oracle<'a> {
    /// The large-prime seam at any width (`r < 2^127`): [`decompose`]
    /// as one set of summands when large primes are off, else the
    /// decomposition with fewest large columns fed to [`LargePrimes`].
    /// The one-word pipeline and `ic_wide` both come through here.
    ///
    /// [`decompose`]: DecompositionOracle::decompose
    #[allow(clippy::too_many_arguments)]
    pub fn relations_for(
        &mut self,
        ctx: &InstanceCtx<CurveB>,
        fb: &FactorBase<PtB>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: PtB,
        a: u128,
        b: u128,
        coef: impl Fn(usize) -> u128,
    ) -> Vec<BinRelation> {
        if !self.lp.enabled() {
            return match self.decompose(ctx, fb, ops, counters, point) {
                Some(summands) => vec![BinRelation::Summands(summands)],
                None => Vec::new(),
            };
        }
        if self.lp.small_cols == 0 {
            let permille = if self.lp.small_permille == 0 {
                1000
            } else {
                self.lp.small_permille
            };
            self.lp.small_cols = (fb.columns as u64 * permille / 1000).max(1) as usize;
        }
        self.last_witnesses.clear();
        if self.decompose(ctx, fb, ops, counters, point).is_none() && self.last_witnesses.is_empty()
        {
            return Vec::new();
        }
        let c = &self.curve.curve;
        let r = self.curve.r;
        let small_cols = self.lp.small_cols;
        // The decomposition with fewest large columns.
        let best = self
            .last_witnesses
            .iter()
            .filter_map(|w| {
                w.iter()
                    .map(|q| {
                        fb.index_of_key(CountedGroup::key(c, q))
                            .map(|i| (fb.col_of[i], coef(i)))
                    })
                    .collect::<Option<Vec<_>>>()
            })
            .min_by_key(|cols| cols.iter().filter(|&&(col, _)| col >= small_cols).count());
        let Some(cols) = best else {
            counters.lift_failures += 1;
            return Vec::new();
        };
        let (rel, mulmods) = self.lp.feed(&cols, a % r, b % r, r);
        self.totals.ops += mulmods;
        *self
            .totals
            .extra
            .entry("lp_merge_mulmods".into())
            .or_insert(0) += mulmods;
        rel.into_iter().map(BinRelation::Row).collect()
    }

    pub fn new(curve: &'a BinCurve, macaulay: u8, seed: u64, opts: BinOptions) -> Self {
        BinD4Oracle {
            curve,
            macaulay,
            seed,
            opts,
            charge_unsolved: false,
            trace: BinTrace::default(),
            pre: None,
            totals: SolverTotals {
                solver: "d4-bin-macaulay".into(),
                op_unit: "F_q multiplications".into(),
                ..Default::default()
            },
            per_add: fq_muls_per_add(&curve.curve),
            lp: LargePrimes::default(),
            last_witnesses: Vec::new(),
        }
    }
    fn shape(&self) -> SystemShape {
        if self.opts.frobenius {
            SystemShape {
                n_vars: 4,
                n_equations: 4,
                degrees: vec![4, 4, 4, 4],
                semi_regular_degree: Some(11),
            }
        } else if self.opts.weighted {
            SystemShape {
                n_vars: 4,
                n_equations: 4,
                degrees: vec![8, 8, 8, 8],
                semi_regular_degree: Some(26),
            }
        } else {
            SystemShape {
                n_vars: 5,
                n_equations: 5,
                degrees: vec![4, 4, 4, 4, 2],
                semi_regular_degree: Some(14),
            }
        }
    }
}

impl DecompositionOracle<CurveB> for BinD4Oracle<'_> {
    fn name(&self) -> &str {
        "gaudry-d4-bin"
    }
    fn summands(&self) -> u32 {
        4
    }
    fn describe(&self, _params: &Params) -> String {
        let form = if self.opts.frobenius {
            "(s₂², s₃², s₄, π²), weights (2, 2, 1, 1): 64 solutions"
        } else if self.opts.weighted {
            "(s₂, s₃, s₄, π), s₁ = π² + π, weights (2, 2, 2, 1): 512 solutions"
        } else {
            "(s₁, s₂, s₃, s₄, π) and π² + π + s₁, unweighted: 512 solutions"
        };
        format!(
            "S₅ in the t-line's D₄ invariants over F_q, {form}, Macaulay degree {}, eigen step along the chains{}{}",
            self.macaulay,
            if self.opts.sparsest { ", sparsest pivots" } else { "" },
            if self.opts.replay { ", shape replay" } else { "" }
        )
    }
    fn prepare(
        &mut self,
        _ctx: &InstanceCtx<CurveB>,
        _fb: &FactorBase<PtB>,
        _params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<(), String> {
        let mut rng = StdRng::seed_from_u64(self.seed ^ 0x9e);
        let pre = SymmetrisedBin::precompute(&self.curve.curve, &mut rng, 8)?;
        if self.opts.frobenius && !pre.frobenius_reducible() {
            return Err(
                "gaudry-d4-bin: frobenius=1, but H has an odd power of s₂, s₃ or π on this curve"
                    .into(),
            );
        }
        self.totals.ops += pre.precompute_muls;
        self.totals
            .extra
            .insert("precompute_fq_muls".into(), pre.precompute_muls);
        self.totals
            .extra
            .insert("fp_muls_per_add".into(), self.per_add);
        self.pre = Some(pre);
        Ok(())
    }
    fn decompose(
        &mut self,
        _ctx: &InstanceCtx<CurveB>,
        fb: &FactorBase<PtB>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: PtB,
    ) -> Option<Vec<usize>> {
        let pre = self.pre.as_ref()?;
        let c = &self.curve.curve;
        let mut rng =
            StdRng::seed_from_u64(self.seed ^ self.totals.calls.wrapping_mul(0x9e37_79b9));
        let adds_before = c.ops();
        let t = std::time::Instant::now();
        let sol = solve_d4_bin(
            c,
            pre,
            &point,
            self.opts,
            self.macaulay,
            &mut rng,
            self.opts.replay.then_some(&mut self.trace),
        );
        self.totals.wall_ns += t.elapsed().as_nanos() as u64;
        // The additions that test a candidate's signs (their field
        // multiplications are inside the additions' price, not the solve's).
        ops.adds += c.ops() - adds_before;
        self.totals.calls += 1;
        self.totals.shape = Some(self.shape());
        if self.opts.replay {
            self.totals
                .extra
                .insert("shape_replays".into(), self.trace.replays);
            self.totals
                .extra
                .insert("shape_fallbacks".into(), self.trace.fallbacks);
        }
        let sol = match sol {
            Ok(sol) => sol,
            Err(cost) => {
                *self
                    .totals
                    .extra
                    .entry("degenerate_systems".into())
                    .or_insert(0) += 1;
                if self.charge_unsolved {
                    self.totals.ops += cost.total_muls;
                    *self
                        .totals
                        .extra
                        .entry("degenerate_fq_muls".into())
                        .or_insert(0) += cost.total_muls;
                }
                return None;
            }
        };
        self.totals.ops += sol.solve.total_muls;
        for (k, v) in [
            ("echelon_fq_muls", sol.solve.echelon_muls),
            ("normal_form_fq_muls", sol.solve.nf_muls),
            ("charpoly_fq_muls", sol.solve.charpoly_muls),
            ("rational_solutions", sol.rational_solutions as u64),
        ] {
            *self.totals.extra.entry(k.into()).or_insert(0) += v;
        }
        self.totals.solving_degree_sum += u64::from(self.macaulay);
        self.totals.solving_degree_max =
            self.totals.solving_degree_max.max(u32::from(self.macaulay));
        self.totals.calls_with_a_degree += 1;
        self.last_witnesses = sol.witnesses.clone();
        let w = sol.witness?;
        let idx: Option<Vec<usize>> = w
            .iter()
            .map(|q| fb.index_of_key(CountedGroup::key(c, q)))
            .collect();
        if idx.is_none() {
            counters.lift_failures += 1;
        }
        idx
    }
    fn decompose_relations(
        &mut self,
        ctx: &InstanceCtx<CurveB>,
        fb: &FactorBase<PtB>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: PtB,
        a: u64,
        b: u64,
    ) -> Vec<Relation> {
        let r = self.curve.r;
        let coef = |i: usize| u128::from(fb.coef_of[i]);
        // This pipeline is the one-word one: r < 2^64 here.
        self.relations_for(ctx, fb, ops, counters, point, u128::from(a), u128::from(b), coef)
            .into_iter()
            .map(|rel| match rel {
                BinRelation::Summands(s) => Relation::Summands(s),
                BinRelation::Row(f) => Relation::Row {
                    columns: f.columns.iter().map(|&(c, v)| (c, (v % r) as u64)).collect(),
                    a: f.a as u64,
                    b: f.b as u64,
                },
            })
            .collect()
    }
    fn relation_counters(&self) -> BTreeMap<String, u64> {
        let mut out = self.lp.counts.clone();
        if self.lp.enabled() {
            out.insert("lp_small_columns".into(), self.lp.small_cols as u64);
            out.insert("lp_stored_pivots".into(), self.lp.pivots.len() as u64);
        }
        out
    }
    fn last_system(&self) -> Option<SystemShape> {
        Some(self.shape())
    }
    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_point(c: &CurveB, rng: &mut StdRng) -> PtB {
        loop {
            let x = c.f.random(rng);
            if let Some(p) = c.lift_x(&x, rng) {
                return p;
            }
        }
    }

    #[test]
    fn fields_and_the_group_law() {
        let c = generate_instance(5, 1);
        let f = &c.f;
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..50 {
            let a = f.random(&mut rng);
            if a.is_zero() {
                continue;
            }
            assert_eq!(f.mul(&a, &f.inv(&a)), Eb::ONE);
            let s = f.sqrt(&a);
            assert_eq!(f.sq(&s), a);
            let p = random_point(&c, &mut rng);
            let q = random_point(&c, &mut rng);
            assert!(c.on_curve(&p) && c.on_curve(&c.add(&p, &q)));
            assert_eq!(c.add(&c.add(&p, &q), &c.neg(&q)), p);
            assert!(c.add(&p, &c.neg(&p)).inf);
        }
    }

    #[test]
    fn the_two_torsion_translation_is_t_plus_one() {
        let c = generate_instance(5, 2);
        let mut rng = StdRng::seed_from_u64(2);
        let t = c.torsion();
        assert!(c.on_curve(&t) && c.add(&t, &t).inf);
        for _ in 0..30 {
            let p = random_point(&c, &mut rng);
            let pt = c.add(&p, &t);
            if pt.inf || p.x == c.beta {
                continue;
            }
            assert_eq!(c.t_of_x(&pt.x), c.t_of_x(&p.x).add(&Eb::ONE));
            assert_eq!(c.x_of_t(&c.t_of_x(&p.x)), p.x);
        }
    }

    #[test]
    fn summation_polynomials_vanish_on_sums() {
        let c = generate_instance(5, 3);
        let mut rng = StdRng::seed_from_u64(3);
        for _ in 0..10 {
            let p: Vec<PtB> = (0..4).map(|_| random_point(&c, &mut rng)).collect();
            let s3 = c.neg(&c.add(&p[0], &p[1]));
            assert!(c.s3(&p[0].x, &p[1].x, &s3.x).is_zero());
            let s5 = c.neg(&c.add(&c.add(&p[0], &p[1]), &c.add(&p[2], &p[3])));
            assert!(c.s5(&[p[0].x, p[1].x, p[2].x, p[3].x, s5.x]).is_zero());
            assert!(!c.s5(&[p[0].x, p[1].x, p[2].x, p[3].x, p[0].x]).is_zero());
        }
    }

    #[test]
    fn s5_in_t_lies_in_the_invariants_and_solves_exactly() {
        let c = generate_instance(5, 4);
        let mut rng = StdRng::seed_from_u64(4);
        let pre = SymmetrisedBin::precompute(&c, &mut rng, 6).expect("invariant");
        assert_eq!(pre.monos.len(), 105);
        let base = t_base(&c, &mut rng);
        let mut trace = BinTrace::default();
        let mut checked = 0;
        for k in 0..6u64 {
            // A planted target.
            let idx: Vec<usize> = {
                let mut v = Vec::new();
                while v.len() < 4 {
                    let i = rng.gen_range(0..base.len());
                    if !v.contains(&i) {
                        v.push(i);
                    }
                }
                v
            };
            let mut r = PtB::INF;
            for &i in &idx {
                r = c.add(&r, &base[i].1);
            }
            if r.inf {
                continue;
            }
            let mitm = mitm_t_decompositions(&c, &base, &r);
            for opts in [
                BinOptions {
                    weighted: true,
                    frobenius: false,
                    sparsest: true,
                    replay: true,
                },
                BinOptions {
                    weighted: false,
                    frobenius: false,
                    sparsest: true,
                    replay: false,
                },
                BinOptions {
                    weighted: true,
                    frobenius: true,
                    sparsest: true,
                    replay: false,
                },
            ] {
                let mut s = StdRng::seed_from_u64(k);
                let d = if opts.frobenius {
                    11
                } else if opts.weighted {
                    26
                } else {
                    14
                };
                let tr = if opts.replay { Some(&mut trace) } else { None };
                if let Ok(sol) = solve_d4_bin(&c, &pre, &r, opts, d, &mut s, tr) {
                    assert_eq!(
                        sol.solve.dim,
                        if opts.frobenius { 64 } else { 512 },
                        "{opts:?}"
                    );
                    let got: BTreeSet<[Fe; 4]> = sol.decompositions.into_iter().collect();
                    assert_eq!(got, mitm, "{opts:?}");
                    assert!(!got.is_empty());
                    checked += 1;
                }
            }
        }
        assert!(checked >= 12, "{checked} solves checked");
        assert!(pre.frobenius_reducible());
    }
}

#[cfg(test)]
mod ecbench_tests {
    use super::*;

    #[test]
    fn group_order_and_tower_identity() {
        let c = find_bin_curve(5, 1, 16).unwrap();
        assert_eq!(c.cofactor * c.r, c.group_order);
        assert!(c.curve.mul(&c.generator, c.r).inf);
        // Every random point is killed by the order.
        let mut rng = StdRng::seed_from_u64(5);
        for _ in 0..5 {
            let x = c.curve.f.random(&mut rng);
            if let Some(p) = c.curve.lift_x(&x, &mut rng) {
                assert!(c.curve.mul(&p, c.group_order).inf);
            }
        }
        // The tower's polynomial basis: coordinates are a ring isomorphism
        // onto GF(2)[X]/(modulus).
        let f = &c.curve.f;
        let tb = TowerBasis::new(f);
        let n = 4 * f.f.m;
        let mulmod = |mut a: u128, mut b: u128| -> u128 {
            let mut r = 0u128;
            while b != 0 {
                if b & 1 == 1 {
                    r ^= a;
                }
                b >>= 1;
                a <<= 1;
                if a >> n & 1 == 1 {
                    a ^= tb.modulus;
                }
            }
            r
        };
        for _ in 0..20 {
            let (a, b) = (f.random(&mut rng), f.random(&mut rng));
            assert_eq!(
                tb.coords(f, &f.mul(&a, &b)),
                mulmod(tb.coords(f, &a), tb.coords(f, &b))
            );
        }
    }
}
