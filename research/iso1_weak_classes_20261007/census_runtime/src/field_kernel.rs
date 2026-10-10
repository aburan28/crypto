//! Field and affine group arithmetic frozen from 62becf9572fe74cbe8b3d8cebee3bf8a240708a1.
use crate::scalar_kernel::{am, inv_mod, mm, pow_mod, sm};
use rand::rngs::StdRng;
use rand::Rng;
use std::cell::Cell;

pub trait Fld {
    type E: Copy + Eq + std::fmt::Debug;
    fn zero(&self) -> Self::E;
    fn one(&self) -> Self::E;
    fn add(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn sub(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn neg(&self, a: &Self::E) -> Self::E;
    fn mul(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn inv(&self, a: &Self::E) -> Self::E;
    fn is_zero(&self, a: &Self::E) -> bool {
        *a == self.zero()
    }
}

// ── F_{p²} = F_p[t] / (t² − ω) ───────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
pub struct E2(pub [u64; 2]);

impl E2 {
    pub const ZERO: E2 = E2([0, 0]);
    pub const ONE: E2 = E2([1, 0]);
    pub fn in_fp(&self) -> bool {
        self.0[1] == 0
    }
}

/// `F_{p²}` with a counter of `F_p` multiplications: a product is `5` (four
/// products and the scaling by `ω`), an inverse `21` (the norm, one `F_p`
/// inversion at the harness's `16`, two scalings).
#[derive(Clone)]
pub struct Fq {
    pub p: u64,
    pub w: u64,
    muls: Cell<u64>,
}

impl Fq {
    pub fn new(p: u64) -> Fq {
        assert!(p % 2 == 1 && p > 3);
        let w = (2..p)
            .find(|&w| pow_mod(w, (p - 1) / 2, p) != 1)
            .expect("a non-residue");
        Fq {
            p,
            w,
            muls: Cell::new(0),
        }
    }
    fn count(&self, k: u64) {
        self.muls.set(self.muls.get() + k);
    }
    pub fn count_public(&self, k: u64) {
        self.count(k);
    }
    pub fn muls(&self) -> u64 {
        self.muls.get()
    }
    pub fn reset_muls(&self) {
        self.muls.set(0);
    }
    pub fn from_fp(&self, x: u64) -> E2 {
        E2([x % self.p, 0])
    }
    /// Multiplication by an element of `F_p`: two products.
    pub fn scale(&self, a: &E2, k: u64) -> E2 {
        self.count(2);
        E2([mm(a.0[0], k, self.p), mm(a.0[1], k, self.p)])
    }
    pub fn sq(&self, a: &E2) -> E2 {
        self.mul(a, a)
    }
    pub fn conj(&self, a: &E2) -> E2 {
        E2([a.0[0], (self.p - a.0[1]) % self.p])
    }
    pub fn pow(&self, a: &E2, mut e: u128) -> E2 {
        let mut base = *a;
        let mut acc = E2::ONE;
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul(&acc, &base);
            }
            base = self.sq(&base);
            e >>= 1;
        }
        acc
    }
    fn order(&self) -> u128 {
        (self.p as u128) * (self.p as u128) - 1
    }
    pub fn is_square(&self, a: &E2) -> bool {
        *a == E2::ZERO || self.pow(a, self.order() / 2) == E2::ONE
    }
    pub fn random(&self, rng: &mut StdRng) -> E2 {
        E2([rng.gen_range(0..self.p), rng.gen_range(0..self.p)])
    }
    /// Tonelli–Shanks in `F_{p²}^×`.
    pub fn sqrt(&self, a: &E2, rng: &mut StdRng) -> Option<E2> {
        if *a == E2::ZERO {
            return Some(E2::ZERO);
        }
        if !self.is_square(a) {
            return None;
        }
        let q1 = self.order();
        let s = q1.trailing_zeros();
        let odd = q1 >> s;
        let z = loop {
            let z = self.random(rng);
            if z != E2::ZERO && !self.is_square(&z) {
                break z;
            }
        };
        let mut m = s;
        let mut cc = self.pow(&z, odd);
        let mut t = self.pow(a, odd);
        let mut r = self.pow(a, odd.div_ceil(2));
        while t != E2::ONE {
            let mut i = 0;
            let mut tt = t;
            while tt != E2::ONE {
                tt = self.sq(&tt);
                i += 1;
            }
            let mut b = cc;
            for _ in 0..(m - i - 1) {
                b = self.sq(&b);
            }
            m = i;
            cc = self.sq(&b);
            t = self.mul(&t, &cc);
            r = self.mul(&r, &b);
        }
        Some(r)
    }
}

impl Fld for Fq {
    type E = E2;
    fn zero(&self) -> E2 {
        E2::ZERO
    }
    fn one(&self) -> E2 {
        E2::ONE
    }
    fn add(&self, a: &E2, b: &E2) -> E2 {
        E2([am(a.0[0], b.0[0], self.p), am(a.0[1], b.0[1], self.p)])
    }
    fn sub(&self, a: &E2, b: &E2) -> E2 {
        E2([sm(a.0[0], b.0[0], self.p), sm(a.0[1], b.0[1], self.p)])
    }
    fn neg(&self, a: &E2) -> E2 {
        E2([(self.p - a.0[0]) % self.p, (self.p - a.0[1]) % self.p])
    }
    fn mul(&self, a: &E2, b: &E2) -> E2 {
        let p = self.p;
        self.count(5);
        let m00 = mm(a.0[0], b.0[0], p);
        let m11 = mm(a.0[1], b.0[1], p);
        E2([
            am(m00, mm(self.w, m11, p), p),
            am(mm(a.0[0], b.0[1], p), mm(a.0[1], b.0[0], p), p),
        ])
    }
    fn inv(&self, a: &E2) -> E2 {
        assert!(*a != E2::ZERO, "inverse of zero in F_p²");
        let p = self.p;
        self.count(3 + 16 + 2);
        let norm = sm(
            mm(a.0[0], a.0[0], p),
            mm(self.w, mm(a.0[1], a.0[1], p), p),
            p,
        );
        let ni = inv_mod(norm, p);
        E2([mm(a.0[0], ni, p), mm((p - a.0[1]) % p, ni, p)])
    }
}

// ── F_{p⁶} = F_{p²}[θ] / (θ³ − s) ────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
pub struct E6(pub [E2; 3]);

impl E6 {
    pub const ZERO: E6 = E6([E2::ZERO; 3]);
    pub const ONE: E6 = E6([E2::ONE, E2::ZERO, E2::ZERO]);
    pub fn in_fq(&self) -> bool {
        self.0[1] == E2::ZERO && self.0[2] == E2::ZERO
    }
}

/// `F_{q³}`, `q = p²`, over `F_q`: `s ∈ F_q` a non-cube (so `θ³ − s` is
/// irreducible; it needs `3 | q − 1`, i.e. `p ≢ 0 (mod 3)`), `ζ = s^{(q−1)/3}`
/// a primitive cube root of unity with `θ^q = ζθ`, so the `q`-Frobenius is
/// `σ(a₀ + a₁θ + a₂θ²) = a₀ + ζa₁θ + ζ²a₂θ²`.
#[derive(Clone)]
pub struct Fq3 {
    pub f: Fq,
    pub s: E2,
    pub zeta: E2,
    pub zeta2: E2,
}

impl Fq3 {
    pub fn new(p: u64) -> Fq3 {
        assert!(!p.is_multiple_of(3));
        let f = Fq::new(p);
        let q1 = f.order();
        let s = (1..p * p)
            .map(|k| E2([k % p, (k / p) % p]))
            .find(|s| *s != E2::ZERO && f.pow(s, q1 / 3) != E2::ONE)
            .expect("a non-cube");
        let zeta = f.pow(&s, q1 / 3);
        let zeta2 = f.mul(&zeta, &zeta);
        f.reset_muls();
        Fq3 { f, s, zeta, zeta2 }
    }
    pub fn muls(&self) -> u64 {
        self.f.muls()
    }
    pub fn reset_muls(&self) {
        self.f.reset_muls()
    }
    pub fn from_fq(&self, x: E2) -> E6 {
        E6([x, E2::ZERO, E2::ZERO])
    }
    /// Scale by an element of `F_q`: three products.
    pub fn scale(&self, a: &E6, k: &E2) -> E6 {
        E6(core::array::from_fn(|i| self.f.mul(&a.0[i], k)))
    }
    pub fn sq(&self, a: &E6) -> E6 {
        self.mul(a, a)
    }
    /// The `q`-Frobenius, two products.
    pub fn sigma(&self, a: &E6) -> E6 {
        E6([
            a.0[0],
            self.f.mul(&a.0[1], &self.zeta),
            self.f.mul(&a.0[2], &self.zeta2),
        ])
    }
    pub fn pow(&self, a: &E6, mut e: u128) -> E6 {
        let mut base = *a;
        let mut acc = E6::ONE;
        while e > 0 {
            if e & 1 == 1 {
                acc = self.mul(&acc, &base);
            }
            base = self.sq(&base);
            e >>= 1;
        }
        acc
    }
    fn order(&self) -> u128 {
        (self.f.p as u128).pow(6) - 1
    }
    pub fn is_square(&self, a: &E6) -> bool {
        *a == E6::ZERO || self.pow(a, self.order() / 2) == E6::ONE
    }
    pub fn random(&self, rng: &mut StdRng) -> E6 {
        E6(core::array::from_fn(|_| self.f.random(rng)))
    }
    pub fn sqrt(&self, a: &E6, rng: &mut StdRng) -> Option<E6> {
        if *a == E6::ZERO {
            return Some(E6::ZERO);
        }
        if !self.is_square(a) {
            return None;
        }
        let q1 = self.order();
        let s = q1.trailing_zeros();
        let odd = q1 >> s;
        let z = loop {
            let z = self.random(rng);
            if z != E6::ZERO && !self.is_square(&z) {
                break z;
            }
        };
        let mut m = s;
        let mut cc = self.pow(&z, odd);
        let mut t = self.pow(a, odd);
        let mut r = self.pow(a, odd.div_ceil(2));
        while t != E6::ONE {
            let mut i = 0;
            let mut tt = t;
            while tt != E6::ONE {
                tt = self.sq(&tt);
                i += 1;
            }
            let mut b = cc;
            for _ in 0..(m - i - 1) {
                b = self.sq(&b);
            }
            m = i;
            cc = self.sq(&b);
            t = self.mul(&t, &cc);
            r = self.mul(&r, &b);
        }
        Some(r)
    }
}

impl Fld for Fq3 {
    type E = E6;
    fn zero(&self) -> E6 {
        E6::ZERO
    }
    fn one(&self) -> E6 {
        E6::ONE
    }
    fn add(&self, a: &E6, b: &E6) -> E6 {
        E6(core::array::from_fn(|i| self.f.add(&a.0[i], &b.0[i])))
    }
    fn sub(&self, a: &E6, b: &E6) -> E6 {
        E6(core::array::from_fn(|i| self.f.sub(&a.0[i], &b.0[i])))
    }
    fn neg(&self, a: &E6) -> E6 {
        E6(core::array::from_fn(|i| self.f.neg(&a.0[i])))
    }
    /// Schoolbook with `θ³ = s`: nine products and two reductions.
    fn mul(&self, a: &E6, b: &E6) -> E6 {
        let f = &self.f;
        let mut d = [E2::ZERO; 5];
        for i in 0..3 {
            for j in 0..3 {
                d[i + j] = f.add(&d[i + j], &f.mul(&a.0[i], &b.0[j]));
            }
        }
        E6([
            f.add(&d[0], &f.mul(&self.s, &d[3])),
            f.add(&d[1], &f.mul(&self.s, &d[4])),
            d[2],
        ])
    }
    /// Through the norm to `F_q`: `a⁻¹ = σ(a)σ²(a) / N(a)`.
    fn inv(&self, a: &E6) -> E6 {
        assert!(*a != E6::ZERO, "inverse of zero in F_p⁶");
        let s1 = self.sigma(a);
        let s2 = self.sigma(&s1);
        let b = self.mul(&s1, &s2);
        let n = self.mul(a, &b);
        debug_assert!(n.in_fq());
        let ni = self.f.inv(&n.0[0]);
        self.scale(&b, &ni)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PtE6 {
    pub x: E6,
    pub y: E6,
    pub inf: bool,
}

impl PtE6 {
    pub const INF: PtE6 = PtE6 {
        x: E6::ZERO,
        y: E6::ZERO,
        inf: true,
    };
}

pub struct EllE<'a> {
    pub f: &'a Fq3,
    pub a2: E6,
    pub a4: E6,
    ops: Cell<u64>,
}

impl<'a> EllE<'a> {
    pub fn new(f: &'a Fq3, alpha: &E6) -> EllE<'a> {
        let sa = f.sigma(alpha);
        EllE {
            f,
            a2: f.neg(&f.add(alpha, &sa)),
            a4: f.mul(alpha, &sa),
            ops: Cell::new(0),
        }
    }
    /// `y² = x³ + a₂x² + a₄x` from its coefficients (the model of a curve
    /// with full 2-torsion after one root is moved to `0`; §17's walk).
    pub fn from_a2_a4(f: &'a Fq3, a2: E6, a4: E6) -> EllE<'a> {
        EllE {
            f,
            a2,
            a4,
            ops: Cell::new(0),
        }
    }
    pub fn ops(&self) -> u64 {
        self.ops.get()
    }
    pub fn reset_ops(&self) {
        self.ops.set(0);
    }
    pub fn rhs(&self, x: &E6) -> E6 {
        let f = self.f;
        let x2 = f.sq(x);
        f.add(
            &f.add(&f.mul(&x2, x), &f.mul(&self.a2, &x2)),
            &f.mul(&self.a4, x),
        )
    }
    pub fn on_curve(&self, p: &PtE6) -> bool {
        p.inf || self.f.sq(&p.y) == self.rhs(&p.x)
    }
    pub fn neg(&self, p: &PtE6) -> PtE6 {
        if p.inf {
            *p
        } else {
            PtE6 {
                x: p.x,
                y: self.f.neg(&p.y),
                inf: false,
            }
        }
    }
    /// Affine addition: one inversion, three products.
    pub fn add(&self, p: &PtE6, q: &PtE6) -> PtE6 {
        let f = self.f;
        if p.inf {
            return *q;
        }
        if q.inf {
            return *p;
        }
        self.ops.set(self.ops.get() + 1);
        let lam = if p.x == q.x {
            if p.y == f.neg(&q.y) || p.y == E6::ZERO {
                return PtE6::INF;
            }
            // (3x² + 2a₂x + a₄) / 2y
            let x2 = f.sq(&p.x);
            let num = f.add(
                &f.add(
                    &f.add(&x2, &f.add(&x2, &x2)),
                    &f.mul(&f.add(&self.a2, &self.a2), &p.x),
                ),
                &self.a4,
            );
            f.mul(&num, &f.inv(&f.add(&p.y, &p.y)))
        } else {
            f.mul(&f.sub(&q.y, &p.y), &f.inv(&f.sub(&q.x, &p.x)))
        };
        let x3 = f.sub(&f.sub(&f.sub(&f.sq(&lam), &self.a2), &p.x), &q.x);
        let y3 = f.sub(&f.mul(&lam, &f.sub(&p.x, &x3)), &p.y);
        PtE6 {
            x: x3,
            y: y3,
            inf: false,
        }
    }
    pub fn sub(&self, p: &PtE6, q: &PtE6) -> PtE6 {
        self.add(p, &self.neg(q))
    }
    /// [`EllE::mul`] for a scalar above 2⁶⁴ (the group order itself).
    pub fn mul_u128(&self, pt: &PtE6, mut k: u128) -> PtE6 {
        let mut acc = PtE6::INF;
        let mut base = *pt;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            k >>= 1;
            if k > 0 {
                base = self.add(&base, &base);
            }
        }
        acc
    }
    pub fn mul(&self, pt: &PtE6, mut k: u64) -> PtE6 {
        let mut acc = PtE6::INF;
        let mut base = *pt;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(&acc, &base);
            }
            k >>= 1;
            if k > 0 {
                base = self.add(&base, &base);
            }
        }
        acc
    }
    pub fn random_point(&self, rng: &mut StdRng) -> PtE6 {
        loop {
            let x = self.f.random(rng);
            if let Some(y) = self.f.sqrt(&self.rhs(&x), rng) {
                if y != E6::ZERO {
                    return PtE6 { x, y, inf: false };
                }
            }
        }
    }
}
