//! Round 4, lead 16.2: a closed-form `4E` membership test on the a = 0
//! Koblitz curve `y^2 + xy = x^3 + 1` over `F_{2^n}`, n odd, checked
//! against the `[r]P == O` subgroup check (defender side).
//!
//! ```sh
//! cargo run --release --example r4_16_2_fast_4e -- exhaustive --nmax 21 --out ex.json
//! cargo run --release --example r4_16_2_fast_4e -- sample --n 233 --uniform 200000 \
//!     --balanced 200000 --lib 2000 --seed 1 --out s233.json
//! ```
//!
//! # The test (Prop. 14.2 of `RESEARCH_ISOGENY_FIBERS_ECDLP.md`, re-derived)
//!
//! Doubling `Q = (u, v)`, `u != 0`, on `y^2 + xy = x^3 + 1` gives
//! `lambda = u + v/u`, `x = lambda^2 + lambda`, `y = u^2 + (lambda + 1) x`.
//! So `P = (x, y)` has a rational half iff `lambda^2 + lambda = x` is
//! solvable, iff `Tr(x) = 0`; for odd n the half-trace `H(x)` is a root.
//! A half then has `u^2 = y + (lambda + 1) x`, and since
//! `Tr(u) = Tr(u^2)` and `Tr(x) = 0`,
//! `Tr(u) = Tr(y) + Tr(lambda x)`.  The two roots `lambda`, `lambda + 1`
//! give the two halves `Q`, `Q + T` (`T = (0, 1)`), and their `Tr(u)`
//! differ by `Tr(x) = 0`.  Applying the halving criterion to the half:
//!
//! `P in 4E  <=>  Tr(x) = 0  and  Tr(y) + Tr(x H(x)) = 0`.
//!
//! `O` is in `4E`.  `T = (0, 1)` gives `Tr(0) = 0`, `Tr(1) = n mod 2 = 1`:
//! rejected, correctly.  The 2-part of `E(F_{2^n})` is exactly `Z/4` for
//! odd n (`P4 = (1, 0)` has order 4 and `Tr(1) = 1` means it has no half),
//! so with `#E = 4r`, `r` prime, `4E` is the order-`r` subgroup.
//!
//! Everything here is a native computation; the field and curve arithmetic
//! is a fixed four-word (n <= 256) carry-less implementation local to this
//! example, cross-checked against `crypto_lib::binary_ecc`.  Field
//! multiplications, squarings, inversions, traces and half-trace table
//! lookups are counted.

use std::cell::Cell;
use std::collections::BTreeMap;
use std::hint::black_box;
use std::time::Instant;

use crypto_lib::binary_ecc::curve::{point_add, point_double, scalar_mul};
use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::curve_construction::weil::{curve_order, is_probable_prime};
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

const W: usize = 4;
type Fe = [u64; W];
const ZERO: Fe = [0; W];
const ONE: Fe = [1, 0, 0, 0];

// ───────────────────────────── counters ─────────────────────────────

#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
struct Count {
    mul: u64,
    sqr: u64,
    inv: u64,
    tr: u64,
    htl: u64,
}

impl Count {
    fn sub(self, o: Count) -> Count {
        Count {
            mul: self.mul - o.mul,
            sqr: self.sqr - o.sqr,
            inv: self.inv - o.inv,
            tr: self.tr - o.tr,
            htl: self.htl - o.htl,
        }
    }
}

#[derive(Default)]
struct Ops {
    mul: Cell<u64>,
    sqr: Cell<u64>,
    inv: Cell<u64>,
    tr: Cell<u64>,
    htl: Cell<u64>,
}

fn bump(c: &Cell<u64>, k: u64) {
    c.set(c.get() + k);
}

// ───────────────────────────── the field ─────────────────────────────

fn clmul_sw(a: u64, b: u64) -> (u64, u64) {
    let mut acc = 0u128;
    let a = a as u128;
    let mut b = b;
    let mut i = 0;
    while b != 0 {
        if b & 1 == 1 {
            acc ^= a << i;
        }
        b >>= 1;
        i += 1;
    }
    (acc as u64, (acc >> 64) as u64)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "pclmulqdq")]
unsafe fn clmul_hw(a: u64, b: u64) -> (u64, u64) {
    use std::arch::x86_64::*;
    let x = _mm_set_epi64x(0, a as i64);
    let y = _mm_set_epi64x(0, b as i64);
    let z = _mm_clmulepi64_si128::<0x00>(x, y);
    let lo = _mm_cvtsi128_si64(z) as u64;
    let hi = _mm_cvtsi128_si64(_mm_srli_si128::<8>(z)) as u64;
    (lo, hi)
}

/// `GF(2)[z]/(z^n + sum z^t)`, n <= 256, polynomial basis.
struct Field {
    n: u32,
    nw: usize,
    low: Vec<u32>,
    mask: Fe,
    /// Bit i set iff `Tr(z^i) = 1`.
    tmask: Fe,
    /// `ht_tab[k][v] = H(v * z^{8k})` for the byte value v.
    ht_tab: Vec<[Fe; 256]>,
    hw: bool,
    ops: Ops,
}

fn shr8(c: &[u64; 2 * W], s: u32) -> [u64; 2 * W] {
    let ws = (s / 64) as usize;
    let bs = s % 64;
    let mut out = [0u64; 2 * W];
    for (i, o) in out.iter_mut().enumerate() {
        let j = i + ws;
        if j < 2 * W {
            let mut w = c[j] >> bs;
            if bs != 0 && j + 1 < 2 * W {
                w |= c[j + 1] << (64 - bs);
            }
            *o = w;
        }
    }
    out
}

fn xor_shl8(c: &mut [u64; 2 * W], h: &[u64; 2 * W], t: u32) {
    let ws = (t / 64) as usize;
    let bs = t % 64;
    for i in 0..2 * W {
        if h[i] == 0 {
            continue;
        }
        let j = i + ws;
        if j < 2 * W {
            c[j] ^= h[i] << bs;
        }
        if bs != 0 && j + 1 < 2 * W {
            c[j + 1] ^= h[i] >> (64 - bs);
        }
    }
}

fn add(a: Fe, b: Fe) -> Fe {
    [a[0] ^ b[0], a[1] ^ b[1], a[2] ^ b[2], a[3] ^ b[3]]
}

fn is_zero(a: Fe) -> bool {
    a == ZERO
}

impl Field {
    fn new(n: u32, low: &[u32]) -> Field {
        assert!(n <= 64 * W as u32 && n >= 3);
        let nw = n.div_ceil(64) as usize;
        let mut mask = ZERO;
        for i in 0..n {
            mask[(i / 64) as usize] |= 1u64 << (i % 64);
        }
        #[cfg(target_arch = "x86_64")]
        let hw = std::arch::is_x86_feature_detected!("pclmulqdq");
        #[cfg(not(target_arch = "x86_64"))]
        let hw = false;
        let mut f = Field {
            n,
            nw,
            low: low.to_vec(),
            mask,
            tmask: ZERO,
            ht_tab: Vec::new(),
            hw,
            ops: Ops::default(),
        };
        // Trace mask: Tr(z^i) = sum_{j<n} (z^i)^{2^j}, which is 0 or 1.
        let mut basis_ht = Vec::with_capacity(n as usize);
        for i in 0..n {
            let mut e = ZERO;
            e[(i / 64) as usize] = 1u64 << (i % 64);
            let mut s = e;
            let mut acc = e;
            for _ in 1..n {
                s = f.sqr_raw(s);
                acc = add(acc, s);
            }
            assert!(acc == ZERO || acc == ONE, "trace not in F_2");
            if acc == ONE {
                f.tmask[(i / 64) as usize] |= 1u64 << (i % 64);
            }
            if n % 2 == 1 {
                basis_ht.push(f.half_trace_raw(e));
            }
        }
        if n % 2 == 1 {
            let nbytes = n.div_ceil(8) as usize;
            let mut tab = vec![[ZERO; 256]; nbytes];
            for (k, t) in tab.iter_mut().enumerate() {
                for v in 1..256usize {
                    let mut acc = ZERO;
                    for j in 0..8 {
                        let bit = 8 * k + j;
                        if v >> j & 1 == 1 && bit < n as usize {
                            acc = add(acc, basis_ht[bit]);
                        }
                    }
                    t[v] = acc;
                }
            }
            f.ht_tab = tab;
        }
        f
    }

    fn clmul(&self, a: u64, b: u64) -> (u64, u64) {
        #[cfg(target_arch = "x86_64")]
        if self.hw {
            // SAFETY: pclmulqdq availability was checked at construction.
            return unsafe { clmul_hw(a, b) };
        }
        clmul_sw(a, b)
    }

    fn reduce(&self, mut c: [u64; 2 * W]) -> Fe {
        loop {
            let hi = shr8(&c, self.n);
            if hi.iter().all(|&w| w == 0) {
                break;
            }
            for (i, w) in c.iter_mut().enumerate() {
                if i >= W {
                    *w = 0;
                } else {
                    *w &= self.mask[i];
                }
            }
            for &t in &self.low {
                xor_shl8(&mut c, &hi, t);
            }
        }
        [c[0], c[1], c[2], c[3]]
    }

    fn mul_raw(&self, a: Fe, b: Fe) -> Fe {
        let mut c = [0u64; 2 * W];
        for i in 0..self.nw {
            if a[i] == 0 {
                continue;
            }
            for j in 0..self.nw {
                let (lo, hi) = self.clmul(a[i], b[j]);
                c[i + j] ^= lo;
                c[i + j + 1] ^= hi;
            }
        }
        self.reduce(c)
    }

    fn sqr_raw(&self, a: Fe) -> Fe {
        let mut c = [0u64; 2 * W];
        for i in 0..self.nw {
            let (lo, hi) = self.clmul(a[i], a[i]);
            c[2 * i] = lo;
            c[2 * i + 1] = hi;
        }
        self.reduce(c)
    }

    fn mul(&self, a: Fe, b: Fe) -> Fe {
        bump(&self.ops.mul, 1);
        self.mul_raw(a, b)
    }

    fn sqr(&self, a: Fe) -> Fe {
        bump(&self.ops.sqr, 1);
        self.sqr_raw(a)
    }

    /// Itoh–Tsujii, `a^{2^n - 2}`; counted as one inversion (its inner
    /// products are not added to the multiplication counter).
    fn inv(&self, a: Fe) -> Fe {
        assert!(!is_zero(a));
        bump(&self.ops.inv, 1);
        let m = self.n - 1;
        let bits = 32 - m.leading_zeros();
        let mut beta = a;
        let mut len = 1u32;
        for i in (0..bits - 1).rev() {
            let mut s = beta;
            for _ in 0..len {
                s = self.sqr_raw(s);
            }
            beta = self.mul_raw(beta, s);
            len *= 2;
            if m >> i & 1 == 1 {
                beta = self.mul_raw(self.sqr_raw(beta), a);
                len += 1;
            }
        }
        debug_assert_eq!(len, m);
        self.sqr_raw(beta)
    }

    fn trace(&self, a: Fe) -> u8 {
        bump(&self.ops.tr, 1);
        let mut p = 0u32;
        for i in 0..W {
            p ^= (a[i] & self.tmask[i]).count_ones();
        }
        (p & 1) as u8
    }

    fn half_trace_raw(&self, a: Fe) -> Fe {
        let mut h = a;
        let mut t = a;
        for _ in 0..(self.n - 1) / 2 {
            t = self.sqr_raw(self.sqr_raw(t));
            h = add(h, t);
        }
        h
    }

    /// `H(a) = sum_{i=0}^{(n-1)/2} a^{2^{2i}}`: n - 1 counted squarings.
    fn half_trace_chain(&self, a: Fe) -> Fe {
        let mut h = a;
        let mut t = a;
        for _ in 0..(self.n - 1) / 2 {
            t = self.sqr(self.sqr(t));
            h = add(h, t);
        }
        h
    }

    /// `H` as a linear map: ceil(n/8) byte-table lookups, no field ops.
    fn half_trace_table(&self, a: Fe) -> Fe {
        let nbytes = self.ht_tab.len();
        bump(&self.ops.htl, nbytes as u64);
        let mut acc = ZERO;
        for (k, t) in self.ht_tab.iter().enumerate() {
            let byte = (a[k / 8] >> (8 * (k % 8))) & 0xff;
            acc = add(acc, t[byte as usize]);
        }
        acc
    }

    fn snapshot(&self) -> Count {
        Count {
            mul: self.ops.mul.get(),
            sqr: self.ops.sqr.get(),
            inv: self.ops.inv.get(),
            tr: self.ops.tr.get(),
            htl: self.ops.htl.get(),
        }
    }

    fn random(&self, rng: &mut StdRng) -> Fe {
        let mut a = ZERO;
        for (i, w) in a.iter_mut().enumerate() {
            *w = rng.gen::<u64>() & self.mask[i];
        }
        a
    }

    fn fe_to_lib(&self, a: Fe) -> F2mElement {
        F2mElement::from_words(&a[..self.nw], self.n)
    }

    fn lib_to_fe(&self, e: &F2mElement) -> Fe {
        let mut a = ZERO;
        for (i, w) in e.raw_bits().iter().enumerate().take(W) {
            a[i] = *w;
        }
        a
    }
}

// ───────────────────────────── the curve ─────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Aff {
    Inf,
    Pt(Fe, Fe),
}

/// López–Dahab projective `(X, Y, Z) <-> (X/Z, Y/Z^2)`; `Z = 0` is `O`.
#[derive(Clone, Copy, Debug)]
struct Ld {
    x: Fe,
    y: Fe,
    z: Fe,
}

const LD_INF: Ld = Ld {
    x: ONE,
    y: ZERO,
    z: ZERO,
};

const P4: Aff = Aff::Pt(ONE, ZERO);
const P4_NEG: Aff = Aff::Pt(ONE, ONE);
const T2: Aff = Aff::Pt(ZERO, ONE);

/// The a = 0 Koblitz curve `y^2 + xy = x^3 + 1` over the given field.
struct Curve<'a> {
    f: &'a Field,
}

impl<'a> Curve<'a> {
    fn on_curve(&self, p: Aff) -> bool {
        match p {
            Aff::Inf => true,
            Aff::Pt(x, y) => {
                let f = self.f;
                let lhs = add(f.sqr_raw(y), f.mul_raw(x, y));
                let rhs = add(f.mul_raw(f.sqr_raw(x), x), ONE);
                lhs == rhs
            }
        }
    }

    fn neg(p: Aff) -> Aff {
        match p {
            Aff::Inf => Aff::Inf,
            Aff::Pt(x, y) => Aff::Pt(x, add(x, y)),
        }
    }

    fn aff_double(&self, p: Aff) -> Aff {
        let f = self.f;
        match p {
            Aff::Inf => Aff::Inf,
            Aff::Pt(x, y) => {
                if is_zero(x) {
                    return Aff::Inf;
                }
                let lam = add(x, f.mul(y, f.inv(x)));
                let x3 = add(f.sqr(lam), lam);
                let y3 = add(f.sqr(x), f.mul(add(lam, ONE), x3));
                Aff::Pt(x3, y3)
            }
        }
    }

    fn aff_add(&self, p: Aff, q: Aff) -> Aff {
        let f = self.f;
        match (p, q) {
            (Aff::Inf, r) | (r, Aff::Inf) => r,
            (Aff::Pt(x1, y1), Aff::Pt(x2, y2)) => {
                if x1 == x2 {
                    if add(y1, y2) == x1 {
                        return Aff::Inf;
                    }
                    return self.aff_double(p);
                }
                let lam = f.mul(add(y1, y2), f.inv(add(x1, x2)));
                let x3 = add(add(add(f.sqr(lam), lam), x1), x2);
                let y3 = add(add(f.mul(lam, add(x1, x3)), x3), y1);
                Aff::Pt(x3, y3)
            }
        }
    }

    /// LD doubling, a = 0, b = 1: 3M + 5S.
    fn ld_double(&self, p: Ld) -> Ld {
        if is_zero(p.z) {
            return p;
        }
        let f = self.f;
        let x2 = f.sqr(p.x);
        let z2 = f.sqr(p.z);
        let z3 = f.mul(x2, z2);
        let x4 = f.sqr(x2);
        let z4 = f.sqr(z2);
        let x3 = add(x4, z4);
        let y3 = add(f.mul(z4, z3), f.mul(x3, add(f.sqr(p.y), z4)));
        Ld {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// LD + affine mixed addition, a = 0: 8M + 5S (HMV Alg. 3.25 shape).
    fn ld_madd(&self, p: Ld, q: Aff) -> Ld {
        let (x2, y2) = match q {
            Aff::Inf => return p,
            Aff::Pt(x, y) => (x, y),
        };
        if is_zero(p.z) {
            return Ld {
                x: x2,
                y: y2,
                z: ONE,
            };
        }
        let f = self.f;
        let z1s = f.sqr(p.z);
        let a = add(f.mul(y2, z1s), p.y);
        let b = add(f.mul(x2, p.z), p.x);
        if is_zero(b) {
            if is_zero(a) {
                return self.ld_double(Ld {
                    x: x2,
                    y: y2,
                    z: ONE,
                });
            }
            return LD_INF;
        }
        let c = f.mul(p.z, b);
        let d = f.mul(f.sqr(b), c);
        let z3 = f.sqr(c);
        let e = f.mul(a, c);
        let x3 = add(add(f.sqr(a), d), e);
        let ff = add(x3, f.mul(x2, z3));
        let g = f.mul(add(x2, y2), f.sqr(z3));
        let y3 = add(f.mul(add(e, z3), ff), g);
        Ld {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Frobenius on LD coordinates: 3S.
    fn ld_frob(&self, p: Ld) -> Ld {
        if is_zero(p.z) {
            return p;
        }
        let f = self.f;
        Ld {
            x: f.sqr(p.x),
            y: f.sqr(p.y),
            z: f.sqr(p.z),
        }
    }

    fn ld_to_aff(&self, p: Ld) -> Aff {
        if is_zero(p.z) {
            return Aff::Inf;
        }
        let f = self.f;
        let zi = f.inv(p.z);
        Aff::Pt(f.mul_raw(p.x, zi), f.mul_raw(p.y, f.sqr_raw(zi)))
    }

    /// Uncounted equality of an LD point with an affine point.
    fn ld_eq(&self, p: Ld, q: Aff) -> bool {
        let f = self.f;
        match q {
            Aff::Inf => is_zero(p.z),
            Aff::Pt(x, y) => {
                !is_zero(p.z) && p.x == f.mul_raw(x, p.z) && p.y == f.mul_raw(y, f.sqr_raw(p.z))
            }
        }
    }

    /// Baseline 1: `[k]P` by left-to-right binary double-and-add in LD.
    fn mul_binary(&self, p: Aff, k: &BigUint) -> Ld {
        let mut q = LD_INF;
        for i in (0..k.bits()).rev() {
            q = self.ld_double(q);
            if k.bit(i) {
                q = self.ld_madd(q, p);
            }
        }
        q
    }

    /// Baseline 2: `[k]P` from a precomputed tau-NAF of `k mod (tau^n - 1)`
    /// (digits least significant first), Frobenius 3S, mixed add 8M + 5S.
    fn mul_tnaf(&self, p: Aff, digits: &[i8]) -> Ld {
        let neg = Self::neg(p);
        let mut q = LD_INF;
        for &d in digits.iter().rev() {
            q = self.ld_frob(q);
            match d {
                1 => q = self.ld_madd(q, p),
                -1 => q = self.ld_madd(q, neg),
                _ => {}
            }
        }
        q
    }

    /// Both trace bits of the closed form, always computed (diagnostic):
    /// `(Tr(x), Tr(y) + Tr(x H(x)))`.
    fn test_bits(&self, p: Aff) -> (u8, u8) {
        let f = self.f;
        match p {
            Aff::Inf => (0, 0),
            Aff::Pt(x, y) => {
                let h = f.half_trace_table(x);
                (f.trace(x), f.trace(y) ^ f.trace(f.mul(x, h)))
            }
        }
    }

    /// The closed-form test, half-trace by the squaring chain.
    fn in_4e_chain(&self, p: Aff) -> bool {
        let f = self.f;
        match p {
            Aff::Inf => true,
            Aff::Pt(x, y) => {
                if f.trace(x) != 0 {
                    return false;
                }
                let h = f.half_trace_chain(x);
                f.trace(y) ^ f.trace(f.mul(x, h)) == 0
            }
        }
    }

    /// The closed-form test, half-trace by byte tables.
    fn in_4e_table(&self, p: Aff) -> bool {
        let f = self.f;
        match p {
            Aff::Inf => true,
            Aff::Pt(x, y) => {
                if f.trace(x) != 0 {
                    return false;
                }
                let h = f.half_trace_table(x);
                f.trace(y) ^ f.trace(f.mul(x, h)) == 0
            }
        }
    }

    /// A point with the given x, if any (n odd): w^2 + w = x + 1/x^2,
    /// y = w x.
    fn lift_x(&self, x: Fe, flip: bool) -> Option<Aff> {
        let f = self.f;
        if is_zero(x) {
            return Some(T2);
        }
        let xi = f.inv(x);
        let c = add(x, f.sqr(xi));
        if f.trace(c) != 0 {
            return None;
        }
        let w = f.half_trace_table(c);
        let mut y = f.mul(w, x);
        if flip {
            y = add(y, x);
        }
        Some(Aff::Pt(x, y))
    }

    fn random_point(&self, rng: &mut StdRng) -> Aff {
        loop {
            let x = self.f.random(rng);
            if is_zero(x) {
                continue;
            }
            if let Some(p) = self.lift_x(x, rng.gen::<bool>()) {
                return p;
            }
        }
    }
}

// ───────────────────────────── tau-NAF ─────────────────────────────

fn round_div(a: &BigInt, n: &BigInt) -> BigInt {
    let mut q = a / n;
    let rem = a - &q * n;
    let two_rem: BigInt = &rem * 2;
    if &two_rem > n {
        q += 1;
    } else if two_rem < -n.clone() {
        q -= 1;
    }
    q
}

/// `Z[tau]`, `tau^2 = mu tau - 2`, mu = -1 for a = 0.
fn zt_mul(a: &(BigInt, BigInt), b: &(BigInt, BigInt)) -> (BigInt, BigInt) {
    let mu = BigInt::from(-1);
    let bd = &a.1 * &b.1;
    (&a.0 * &b.0 - &bd * 2, &a.0 * &b.1 + &a.1 * &b.0 + &mu * &bd)
}

fn zt_norm(a: &(BigInt, BigInt)) -> BigInt {
    &a.0 * &a.0 - &a.0 * &a.1 + &a.1 * &a.1 * 2
}

/// tau-NAF digits of `k mod (tau^n - 1)`; returns digits and `N(tau^n-1)`.
fn tnaf_digits(k: &BigUint, n: u32) -> (Vec<i8>, BigInt, BigInt) {
    // tau^n = a + b tau.
    let (mut a, mut b) = (BigInt::one(), BigInt::zero());
    for _ in 0..n {
        let na: BigInt = -(&b * BigInt::from(2));
        let nb: BigInt = &a - &b;
        a = na;
        b = nb;
    }
    let delta = (a - 1, b);
    let norm = zt_norm(&delta);
    let kk = BigInt::from_biguint(Sign::Plus, k.clone());
    // conj(d0 + d1 tau) = (d0 + mu d1) - d1 tau.
    let l0 = round_div(&(&kk * (&delta.0 - &delta.1)), &norm);
    let l1 = round_div(&(-(&kk * &delta.1)), &norm);
    let kd = zt_mul(&(l0, l1), &delta);
    let mut r0 = &kk - kd.0;
    let mut r1 = -kd.1;
    let rho_norm = zt_norm(&(r0.clone(), r1.clone()));
    let mut digits = Vec::new();
    let four = BigInt::from(4);
    while !r0.is_zero() || !r1.is_zero() {
        let u: i8 = if (&r0 % 2u32) != BigInt::zero() {
            let m: BigInt = (((&r0 - &r1 * BigInt::from(2)) % &four) + &four) % &four;
            let u = 2 - m.to_i64().expect("small") as i8;
            r0 -= u;
            u
        } else {
            0
        };
        digits.push(u);
        let half: BigInt = &r0 / 2; // exact
        let nr0 = &r1 - &half; // r1 + mu * r0/2
        r1 = -half;
        r0 = nr0;
    }
    (digits, norm, rho_norm)
}

// ───────────────────────────── helpers ─────────────────────────────

fn small_irreducible(n: u32) -> Vec<u32> {
    // Ben-Or test over u64 polynomials (n <= 62).
    fn mulmod(a: u64, b: u64, f: u64, n: u32) -> u64 {
        let mut r = 0u64;
        for i in (0..n).rev() {
            r <<= 1;
            if r >> n & 1 == 1 {
                r ^= f;
            }
            if b >> i & 1 == 1 {
                r ^= a;
            }
        }
        r
    }
    fn deg(a: u64) -> i32 {
        63 - a.leading_zeros() as i32
    }
    fn gcd(mut a: u64, mut b: u64) -> u64 {
        while b != 0 {
            while a != 0 && deg(a) >= deg(b) {
                a ^= b << (deg(a) - deg(b));
            }
            std::mem::swap(&mut a, &mut b);
        }
        a
    }
    let irreducible = |f: u64| {
        let mut g = 2u64;
        for _ in 1..=n / 2 {
            g = mulmod(g, g, f, n);
            if gcd(f, g ^ 2) != 1 {
                return false;
            }
        }
        true
    };
    for k in 1..n {
        let f = (1u64 << n) | (1u64 << k) | 1;
        if irreducible(f) {
            return vec![0, k];
        }
    }
    for k1 in 3..n {
        for k2 in 2..k1 {
            for k3 in 1..k2 {
                let f = (1u64 << n) | (1u64 << k1) | (1u64 << k2) | (1u64 << k3) | 1;
                if irreducible(f) {
                    return vec![0, k3, k2, k1];
                }
            }
        }
    }
    panic!("no irreducible found for n = {n}");
}

fn cnt_json(c: Count, checks: u64) -> Value {
    let d = checks.max(1) as f64;
    json!({
        "checks": checks,
        "mul_per_check": c.mul as f64 / d,
        "sqr_per_check": c.sqr as f64 / d,
        "inv_per_check": c.inv as f64 / d,
        "trace_per_check": c.tr as f64 / d,
        "halftrace_table_lookups_per_check": c.htl as f64 / d,
    })
}

fn hex(a: Fe, nw: usize) -> String {
    a[..nw]
        .iter()
        .rev()
        .map(|w| format!("{w:016x}"))
        .collect::<String>()
}

fn arg<'a>(args: &'a [String], key: &str) -> Option<&'a str> {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1))
        .map(|s| s.as_str())
}

// ───────────────────────────── modes ─────────────────────────────

/// Exhaustive check over every point for small odd n.
fn exhaustive(nmax: u32) -> Value {
    let mut rows = Vec::new();
    let mut n = 5;
    while n <= nmax {
        let low = small_irreducible(n);
        let f = Field::new(n, &low);
        let c = Curve { f: &f };
        let order = curve_order(-1, n).to_biguint().expect("positive");
        let q = &order / 4u32;
        let two_part_is_4 = (&order % 8u32) == BigUint::from(4u32);
        let mut points = 1u64; // O
        let mut agree = 0u64;
        let mut fp = 0u64;
        let mut fneg = 0u64;
        let mut in4e = 1u64; // O
        let mut test_in = 1u64;
        let t0 = Instant::now();
        let mut handle = |p: Aff| {
            let truth = is_zero(c.mul_binary(p, &q).z);
            let t = c.in_4e_chain(p);
            let t2 = c.in_4e_table(p);
            assert_eq!(t, t2);
            in4e += truth as u64;
            test_in += t as u64;
            if t == truth {
                agree += 1;
            } else if t {
                fp += 1;
            } else {
                fneg += 1;
            }
        };
        for xv in 0u64..(1u64 << n) {
            let x: Fe = [xv, 0, 0, 0];
            if xv == 0 {
                points += 1;
                handle(T2);
                continue;
            }
            if let Some(p) = c.lift_x(x, false) {
                assert!(c.on_curve(p));
                points += 2;
                handle(p);
                handle(Curve::neg(p));
            }
        }
        let ok = BigUint::from(points) == order && fp == 0 && fneg == 0 && BigUint::from(in4e) == q;
        println!(
            "exhaustive n={n} modulus low terms {low:?}: #E={points} (expected {order}), \
             |4E|={in4e}, test-accepted={test_in}, disagreements fp={fp} fn={fneg} [{:.1}s] {}",
            t0.elapsed().as_secs_f64(),
            if ok { "OK" } else { "FAIL" }
        );
        rows.push(json!({
            "n": n,
            "modulus_low_terms": low,
            "points_enumerated": points,
            "curve_order_formula": order.to_string(),
            "order_over_4": q.to_string(),
            "order_over_4_is_prime": is_probable_prime(&q),
            "two_part_exactly_4": two_part_is_4,
            "members_4E_by_scalar": in4e,
            "accepted_by_test": test_in,
            "false_accept": fp,
            "false_reject": fneg,
            "agree": agree,
            "ok": ok,
        }));
        n += 2;
    }
    json!({"mode": "exhaustive", "rows": rows})
}

#[allow(clippy::too_many_arguments)]
fn sample(
    n: u32,
    n_uniform: usize,
    n_bal: usize,
    n_lib: usize,
    seed: u64,
    batch: usize,
    reps: usize,
) -> Value {
    let (irr, lib_curve_opt) = match n {
        131 => (IrreduciblePoly::deg_131(), None),
        233 => (IrreduciblePoly::deg_233(), Some(BinaryCurve::sect233k1())),
        _ => panic!("n must be 131 or 233"),
    };
    let f = Field::new(n, &irr.low_terms);
    let c = Curve { f: &f };
    let mut rng = StdRng::seed_from_u64(seed);
    let mut checks = BTreeMap::<String, Value>::new();

    // ── field cross-checks against crypto_lib ──
    let mut field_ok = true;
    for _ in 0..2000 {
        let a = f.random(&mut rng);
        let b = f.random(&mut rng);
        let la = f.fe_to_lib(a);
        let lb = f.fe_to_lib(b);
        field_ok &= f.lib_to_fe(&la.mul(&lb, &irr)) == f.mul_raw(a, b);
        field_ok &= f.lib_to_fe(&la.square(&irr)) == f.sqr_raw(a);
        if !is_zero(a) {
            field_ok &= f.mul_raw(f.inv(a), a) == ONE;
        }
        let h = f.half_trace_chain(a);
        field_ok &= h == f.half_trace_table(a);
        let tr = f.trace(a);
        let rhs = if tr == 1 { add(a, ONE) } else { a };
        field_ok &= add(f.sqr_raw(h), h) == rhs;
    }
    checks.insert(
        "field_vs_library_and_halftrace_identity".into(),
        json!(field_ok),
    );

    // ── order, r, tau-NAF ──
    let order = curve_order(-1, n).to_biguint().expect("positive");
    let r = &order / 4u32;
    let r_prime = is_probable_prime(&r);
    let mut r_checks = json!({
        "order_formula": order.to_string(),
        "order_mod_8": (&order % 8u32).to_string(),
        "r": r.to_string(),
        "r_bits": r.bits(),
        "r_popcount": r.count_ones(),
        "r_probable_prime": r_prime,
    });
    if let Some(lc) = &lib_curve_opt {
        r_checks["r_equals_library_sect233k1_order"] = json!(lc.order == r);
        r_checks["library_cofactor"] = json!(lc.cofactor.to_string());
    } else {
        let lit = BigUint::parse_bytes(b"680564733841876926932320129493409985129", 10).unwrap();
        r_checks["r_equals_ECC2K130_literature_value"] = json!(r == lit);
    }
    let (digits, dnorm, rho_norm) = tnaf_digits(&r, n);
    let nz = digits.iter().filter(|&&d| d != 0).count();
    r_checks["norm_tau_n_minus_1_equals_order"] = json!(dnorm == BigInt::from(order.clone()));
    r_checks["tnaf_length"] = json!(digits.len());
    r_checks["tnaf_nonzero"] = json!(nz);
    r_checks["log2_norm_rho"] = json!(rho_norm.bits());
    println!("n={n}: {r_checks}");
    checks.insert("order".into(), r_checks);

    // ── generator of G ──
    let g = if let Some(lc) = &lib_curve_opt {
        match &lc.generator {
            BinaryPoint::Affine { x, y } => Aff::Pt(f.lib_to_fe(x), f.lib_to_fe(y)),
            BinaryPoint::Infinity => unreachable!(),
        }
    } else {
        let mut grng = StdRng::seed_from_u64(0x4e);
        loop {
            let p = c.random_point(&mut grng);
            let q = c.aff_double(c.aff_double(p));
            if q != Aff::Inf {
                break q;
            }
        }
    };
    assert!(c.on_curve(g));
    assert!(is_zero(c.mul_binary(g, &r).z), "[r]G != O");

    // ── LD / tau-NAF arithmetic cross-checks ──
    let lib_curve = lib_curve_opt.clone().unwrap_or(BinaryCurve {
        m: n,
        irreducible: irr.clone(),
        a: F2mElement::zero(n),
        b: F2mElement::one(n),
        generator: match g {
            Aff::Pt(x, y) => BinaryPoint::Affine {
                x: f.fe_to_lib(x),
                y: f.fe_to_lib(y),
            },
            Aff::Inf => unreachable!(),
        },
        order: r.clone(),
        cofactor: BigUint::from(4u32),
    });
    let to_lib_pt = |p: Aff| match p {
        Aff::Inf => BinaryPoint::Infinity,
        Aff::Pt(x, y) => BinaryPoint::Affine {
            x: f.fe_to_lib(x),
            y: f.fe_to_lib(y),
        },
    };
    let mut arith_ok = true;
    for _ in 0..20 {
        let p = c.random_point(&mut rng);
        let q = c.random_point(&mut rng);
        arith_ok &=
            to_lib_pt(c.aff_add(p, q)) == point_add(&lib_curve, &to_lib_pt(p), &to_lib_pt(q));
        arith_ok &= to_lib_pt(c.aff_double(p)) == point_double(&lib_curve, &to_lib_pt(p));
        let ld = c.ld_madd(c.ld_madd(LD_INF, p), q);
        arith_ok &= c.ld_to_aff(ld) == c.aff_add(p, q);
        arith_ok &= c.ld_to_aff(c.ld_double(c.ld_madd(LD_INF, p))) == c.aff_double(p);
        let k: BigUint = BigUint::from(rng.gen::<u64>()) << 40u32 | BigUint::from(rng.gen::<u64>());
        arith_ok &= to_lib_pt(c.ld_to_aff(c.mul_binary(p, &k)))
            == scalar_mul(&lib_curve, &to_lib_pt(p), &k);
        arith_ok &= c.ld_to_aff(c.mul_binary(p, &r)) == c.ld_to_aff(c.mul_tnaf(p, &digits));
    }
    checks.insert("ld_tnaf_vs_library_affine".into(), json!(arith_ok));
    assert!(field_ok && arith_ok);

    // coset label of P in E/4E = Z/4 relative to P4: [r]P = [k r] P4.
    let r_mod4 = (&r % 4u32).to_u64().unwrap();
    let label = |rp: Ld| -> usize {
        let j = if is_zero(rp.z) {
            0
        } else if c.ld_eq(rp, P4) {
            1
        } else if c.ld_eq(rp, T2) {
            2
        } else if c.ld_eq(rp, P4_NEG) {
            3
        } else {
            panic!("[r]P has order not dividing 4: #E != 4r");
        };
        (j as u64 * r_mod4 % 4) as usize
    };

    // ── build the sample ──
    let mut pts: Vec<(Aff, &'static str)> = Vec::new();
    let edge = [
        (Aff::Inf, "edge:O"),
        (T2, "edge:T"),
        (P4, "edge:P4"),
        (P4_NEG, "edge:-P4"),
        (g, "edge:G"),
        (Curve::neg(g), "edge:-G"),
        (c.aff_add(g, T2), "edge:G+T"),
        (c.aff_add(g, P4), "edge:G+P4"),
        (c.aff_add(g, P4_NEG), "edge:G-P4"),
    ];
    pts.extend(edge.iter().copied());
    for _ in 0..n_uniform {
        pts.push((c.random_point(&mut rng), "uniform"));
    }
    let p4m = [Aff::Inf, P4, T2, P4_NEG];
    for _ in 0..n_bal {
        let p = c.random_point(&mut rng);
        let q = c.aff_double(c.aff_double(p));
        for (k, m) in p4m.iter().enumerate() {
            let tag = ["balanced:k0", "balanced:k1", "balanced:k2", "balanced:k3"][k];
            pts.push((c.aff_add(q, *m), tag));
        }
    }
    for (p, _) in &pts {
        assert!(c.on_curve(*p));
    }
    println!("n={n}: {} points sampled", pts.len());

    // ── classify ──
    let mut confusion = [[0u64; 2]; 4]; // [true coset][test says in 4E]
    let mut bits_by_coset = [[0u64; 4]; 4]; // [coset][2*b1 + b2]
    let mut by_stream: BTreeMap<&str, [u64; 4]> = BTreeMap::new();
    let mut tnaf_disagree = 0u64;
    let mut chain_table_disagree = 0u64;
    let (mut c_bin_g, mut c_tnaf_g, mut c_chain_g, mut c_table_g) = (
        Count::default(),
        Count::default(),
        Count::default(),
        Count::default(),
    );
    let (mut c_bin_all, mut c_tnaf_all, mut c_chain_all, mut c_table_all) = (
        Count::default(),
        Count::default(),
        Count::default(),
        Count::default(),
    );
    let mut n_g = 0u64;
    let mut chain_minmax = (u64::MAX, 0u64);
    let mut bin_minmax_g = (u64::MAX, 0u64);
    let t_class = Instant::now();
    for (p, stream) in &pts {
        let s0 = f.snapshot();
        let rb = c.mul_binary(*p, &r);
        let s1 = f.snapshot();
        let rt = c.mul_tnaf(*p, &digits);
        let s2 = f.snapshot();
        let ic = c.in_4e_chain(*p);
        let s3 = f.snapshot();
        let it = c.in_4e_table(*p);
        let s4 = f.snapshot();
        let (db, dt, dc, dtab) = (s1.sub(s0), s2.sub(s1), s3.sub(s2), s4.sub(s3));
        let k = label(rb);
        if is_zero(rb.z) != is_zero(rt.z) || (!is_zero(rb.z) && label(rt) != k) {
            tnaf_disagree += 1;
        }
        if ic != it {
            chain_table_disagree += 1;
        }
        confusion[k][ic as usize] += 1;
        let (b1, b2) = c.test_bits(*p);
        bits_by_coset[k][(2 * b1 + b2) as usize] += 1;
        by_stream.entry(stream).or_insert([0; 4])[k] += 1;
        let acc = |a: &mut Count, d: Count| {
            a.mul += d.mul;
            a.sqr += d.sqr;
            a.inv += d.inv;
            a.tr += d.tr;
            a.htl += d.htl;
        };
        acc(&mut c_bin_all, db);
        acc(&mut c_tnaf_all, dt);
        acc(&mut c_chain_all, dc);
        acc(&mut c_table_all, dtab);
        chain_minmax = (chain_minmax.0.min(dc.sqr), chain_minmax.1.max(dc.sqr));
        if k == 0 && *p != Aff::Inf {
            n_g += 1;
            acc(&mut c_bin_g, db);
            acc(&mut c_tnaf_g, dt);
            acc(&mut c_chain_g, dc);
            acc(&mut c_table_g, dtab);
            bin_minmax_g = (bin_minmax_g.0.min(db.mul), bin_minmax_g.1.max(db.mul));
        }
    }
    let class_secs = t_class.elapsed().as_secs_f64();
    let errors = confusion[0][0] + confusion[1][1] + confusion[2][1] + confusion[3][1];
    println!("n={n}: confusion [coset][test_in] = {confusion:?}; errors = {errors}; tnaf disagreements = {tnaf_disagree} [{class_secs:.1}s]");

    // ── library oracle on a subsample (is_valid_public_point rejects O by design) ──
    let mut lib_agree = 0u64;
    let mut lib_disagree = 0u64;
    let mut lib_checked = 0u64;
    let t_lib = Instant::now();
    for (p, _) in pts.iter().filter(|(p, _)| *p != Aff::Inf).take(n_lib) {
        let lv = lib_curve.is_valid_public_point(&to_lib_pt(*p));
        if lv == c.in_4e_table(*p) {
            lib_agree += 1;
        } else {
            lib_disagree += 1;
        }
        lib_checked += 1;
    }
    let lib_ns = t_lib.elapsed().as_nanos() as f64 / lib_checked.max(1) as f64;
    println!("n={n}: library is_valid_public_point agreement {lib_agree}/{lib_checked}");

    // ── analytic expectations for members of G ──
    let l = r.bits();
    let w = r.count_ones();
    let analytic = json!({
        "binary_ld": {"mul": 3 * (l - 1) + 8 * (w - 1), "sqr": 5 * (l - 1) + 5 * (w - 1),
                       "formula": "3(L-1)+8(w-1) M, 5(L-1)+5(w-1) S, L=bitlen(r), w=popcount(r)"},
        "tnaf_ld": {"mul": 8 * (nz as u64 - 1), "sqr": 3 * (digits.len() as u64 - 1) + 5 * (nz as u64 - 1),
                     "formula": "8(z-1) M, 3(len-1)+5(z-1) S, len/z = tau-NAF length/nonzeros"},
        "closed_form_chain": {"mul": 1, "sqr": n - 1, "trace": 3},
        "closed_form_table": {"mul": 1, "sqr": 0, "trace": 3, "halftrace_table_lookups": n.div_ceil(8)},
    });

    // ── contended wall time ──
    let bench_uniform: Vec<Aff> = pts
        .iter()
        .filter(|(_, s)| *s == "uniform")
        .map(|(p, _)| *p)
        .take(batch)
        .collect();
    let bench_g: Vec<Aff> = pts
        .iter()
        .filter(|(_, s)| *s == "balanced:k0")
        .map(|(p, _)| *p)
        .take(batch)
        .collect();
    let time = |name: &str, set: &[Aff], fun: &dyn Fn(Aff) -> bool| -> Value {
        let mut per = Vec::new();
        for _ in 0..reps {
            let t = Instant::now();
            let mut acc = 0u64;
            for p in set {
                acc += fun(black_box(*p)) as u64;
            }
            black_box(acc);
            per.push(t.elapsed().as_nanos() as f64 / set.len() as f64);
        }
        per.sort_by(|a, b| a.partial_cmp(b).unwrap());
        json!({"method": name, "points": set.len(), "reps": reps,
               "contended_wall_ns_per_check_median": per[per.len() / 2],
               "contended_wall_ns_per_check_min": per[0],
               "contended_wall_ns_per_check_max": per[per.len() - 1]})
    };
    let mut bench = Vec::new();
    for (setname, set) in [("uniform", &bench_uniform), ("members_of_G", &bench_g)] {
        let rows = vec![
            time("closed_form_chain", set, &|p| c.in_4e_chain(p)),
            time("closed_form_table", set, &|p| c.in_4e_table(p)),
            time("binary_ld_rP", set, &|p| is_zero(c.mul_binary(p, &r).z)),
            time("tnaf_ld_rP", set, &|p| is_zero(c.mul_tnaf(p, &digits).z)),
        ];
        for row in &rows {
            println!("n={n} bench[{setname}]: {row}");
        }
        bench.push(json!({"set": setname, "rows": rows}));
    }

    let confusion_json: Vec<Value> = (0..4)
        .map(|k| {
            json!({"true_coset": k, "in_4E": k == 0,
                   "test_reject": confusion[k][0], "test_accept": confusion[k][1],
                   "bits_00": bits_by_coset[k][0], "bits_01": bits_by_coset[k][1],
                   "bits_10": bits_by_coset[k][2], "bits_11": bits_by_coset[k][3]})
        })
        .collect();
    let streams: BTreeMap<String, Value> = by_stream
        .iter()
        .map(|(s, v)| ((*s).to_string(), json!(v)))
        .collect();
    json!({
        "mode": "sample",
        "curve": format!("y^2 + xy = x^3 + 1 over F_2^{n} (a = 0 Koblitz){}",
                         if n == 233 { ", sect233k1 parameters from crypto_lib" } else { ", the ECC2K-130 curve" }),
        "modulus_low_terms": irr.low_terms,
        "seed": seed,
        "clmul_kernel": if f.hw { "pclmulqdq" } else { "portable" },
        "generator_x": match g { Aff::Pt(x, _) => hex(x, f.nw), Aff::Inf => String::new() },
        "checks": checks,
        "points_total": pts.len(),
        "points_by_stream_and_true_coset": streams,
        "confusion": confusion_json,
        "errors_total": errors,
        "tnaf_vs_binary_disagreements": tnaf_disagree,
        "chain_vs_table_disagreements": chain_table_disagree,
        "library_oracle": {"checked": lib_checked, "agree": lib_agree, "disagree": lib_disagree,
                           "contended_wall_ns_per_check": lib_ns,
                           "note": "BinaryCurve::is_valid_public_point (on-curve + affine [r]P); O excluded since it returns false for O by design"},
        "counted_ops_members_of_G": {
            "points": n_g,
            "binary_ld_rP": cnt_json(c_bin_g, n_g),
            "tnaf_ld_rP": cnt_json(c_tnaf_g, n_g),
            "closed_form_chain": cnt_json(c_chain_g, n_g),
            "closed_form_table": cnt_json(c_table_g, n_g),
            "binary_mul_min_max": [bin_minmax_g.0, bin_minmax_g.1],
        },
        "counted_ops_all_points": {
            "points": pts.len(),
            "binary_ld_rP": cnt_json(c_bin_all, pts.len() as u64),
            "tnaf_ld_rP": cnt_json(c_tnaf_all, pts.len() as u64),
            "closed_form_chain": cnt_json(c_chain_all, pts.len() as u64),
            "closed_form_table": cnt_json(c_table_all, pts.len() as u64),
            "closed_form_chain_sqr_min_max": [chain_minmax.0, chain_minmax.1],
        },
        "analytic_members_of_G": analytic,
        "bench": bench,
        "classification_wall_seconds": class_secs,
    })
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mode = args.first().map(|s| s.as_str()).unwrap_or("sample");
    let out = arg(&args, "--out");
    let v = match mode {
        "exhaustive" => exhaustive(arg(&args, "--nmax").map_or(19, |s| s.parse().unwrap())),
        "sample" => {
            let p = |k: &str, d: usize| arg(&args, k).map_or(d, |s| s.parse().unwrap());
            sample(
                p("--n", 233) as u32,
                p("--uniform", 1000),
                p("--balanced", 1000),
                p("--lib", 200),
                p("--seed", 1) as u64,
                p("--bench-batch", 2048),
                p("--bench-reps", 5),
            )
        }
        _ => panic!("mode must be exhaustive or sample"),
    };
    let s = serde_json::to_string_pretty(&v).unwrap();
    if let Some(path) = out {
        std::fs::write(path, &s).unwrap();
        println!("wrote {path}");
    } else {
        println!("{s}");
    }
}
