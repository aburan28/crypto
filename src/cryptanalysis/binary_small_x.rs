//! # Small-x factor bases over binary fields, and their exact pair-sum yield.
//!
//! The binary counterpart of [`super::prime_small_x`], for the SECG curves
//! `y² + x·y = x³ + a·x² + b` over `F_{2^m}` with `m` odd and at most 256
//! (`sect113r1`, `sect113r2`, `sect131r1/r2`, `sect163k1/r1/r2`,
//! `sect193r1/r2`, `sect233k1/r1`, `sect239k1`).
//!
//! The base is the points with the smallest abscissae `x = 1, 2, 3, …`
//! (the polynomial whose coefficient bits are the binary digits of `x`)
//! that lie on the curve.  Every SECG binary curve has cofactor `h ≥ 2`, and
//! a decomposition `T = ±a ± b` of an element `T` of the prime-order
//! subgroup `⟨G⟩ = [h]E` is a relation in that subgroup only when the base
//! points are in it too.  So each small-x point `P` enters the base as
//! `[h]P`; the abscissa is then no longer small, but the count below is the
//! count of a base that lies in `⟨G⟩`, which is what a relation needs.
//! (Keeping the raw small-x points and accepting only the pairs that sum into
//! `⟨G⟩` gives about `1/h` of this yield.)
//!
//! A pair of base points `a, b` has the sums `a + b` and `a − b`, with
//! abscissae from one shared inverse: with `λ = (y_a + y_b)/(x_a + x_b)`,
//! `x(a + b) = λ² + λ + x_a + x_b + a`, and `−b = (x_b, x_b + y_b)` gives
//! `x(a − b)` from `λ' = (y_a + y_b + x_b)/(x_a + x_b)`.  Doubling is
//! `x(2a) = x² + b/x²`.  The covered targets of a negation-closed base are
//! twice the number of distinct sum abscissae `D`, and the yield per trial is
//! `2D/(n − 1)` for `n` the order of `⟨G⟩`, exact as in the prime case.

use std::collections::HashMap;
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::{rngs::StdRng, Rng, SeedableRng};
use rayon::prelude::*;

use super::prime_small_x::{
    binomial_ci, check_memory, log2_big, raw_base_summary, Options, RawBase, Timing, YieldReport,
};
use super::secg_standards::{SecgCurve, SecgField};

/// A field element: coefficient bit `i` is bit `i % 64` of limb `i / 64`.
pub type Elem<const L: usize> = [u64; L];

/// Largest extension degree held (four limbs).
pub const MAX_DEGREE: u32 = 576;

/// Limbs of the double-width product buffer, enough for nine-limb fields.
const BUF: usize = 18;

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

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "aes")]
unsafe fn clmul_hw(a: u64, b: u64) -> (u64, u64) {
    let z = std::arch::aarch64::vmull_p64(a, b);
    (z as u64, (z >> 64) as u64)
}

fn clmul_portable(a: u64, b: u64) -> (u64, u64) {
    let (mut lo, mut hi) = (0u64, 0u64);
    let mut b = b;
    let mut i = 0;
    while b != 0 {
        if b & 1 == 1 {
            lo ^= a << i;
            if i > 0 {
                hi ^= a >> (64 - i);
            }
        }
        b >>= 1;
        i += 1;
    }
    (lo, hi)
}

/// `dst ^= src << s` on the low `n` limbs of a double-width value.
fn xor_shl(dst: &mut [u64; BUF], src: &[u64; BUF], s: u32, n: usize) {
    let (w, b) = ((s / 64) as usize, s % 64);
    for i in (w..n).rev() {
        let mut v = src[i - w] << b;
        if b > 0 && i > w {
            v |= src[i - w - 1] >> (64 - b);
        }
        dst[i] ^= v;
    }
}

/// `src >> s` on the low `n` limbs of a double-width value.
fn shr(src: &[u64; BUF], s: u32, n: usize) -> [u64; BUF] {
    let (w, b) = ((s / 64) as usize, s % 64);
    let mut out = [0u64; BUF];
    for i in 0..n.saturating_sub(w) {
        let mut v = src[i + w] >> b;
        if b > 0 && i + w + 1 < n {
            v |= src[i + w + 1] << (64 - b);
        }
        out[i] = v;
    }
    out
}

/// The field `F_2[z]/(z^m + Σ z^t)` with `3 ≤ m ≤ 64·L`.
#[derive(Clone, Debug)]
pub struct Gf2m<const L: usize> {
    pub m: u32,
    limbs: usize,
    low: Vec<u32>,
    has_clmul: bool,
}

impl<const L: usize> Gf2m<L> {
    pub fn new(m: u32, low_terms: &[u32]) -> Result<Self, String> {
        if !(3..=MAX_DEGREE).contains(&m) || m as usize > 64 * L {
            return Err(format!(
                "degree {m} is outside 3..={} for {L} limbs",
                64 * L
            ));
        }
        let mut low: Vec<u32> = low_terms.to_vec();
        low.sort_unstable();
        low.dedup();
        if low.is_empty() || low.iter().any(|&t| t >= m) || low[0] != 0 {
            return Err("the modulus needs a constant term and terms below z^m".into());
        }
        #[cfg(target_arch = "x86_64")]
        let has_clmul = std::arch::is_x86_feature_detected!("pclmulqdq");
        #[cfg(target_arch = "aarch64")]
        let has_clmul = std::arch::is_aarch64_feature_detected!("aes");
        #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
        let has_clmul = false;
        Ok(Self {
            m,
            limbs: (m as usize).div_ceil(64),
            low,
            has_clmul,
        })
    }

    /// The same field with the carry-less kernel chosen by the caller, so a
    /// test can compare both kernels on one machine.
    pub fn with_portable_kernel(&self) -> Self {
        Self {
            has_clmul: false,
            ..self.clone()
        }
    }

    pub fn zero(&self) -> Elem<L> {
        [0; L]
    }

    pub fn one(&self) -> Elem<L> {
        let mut e = [0u64; L];
        e[0] = 1;
        e
    }

    pub fn is_zero(&self, a: &Elem<L>) -> bool {
        a.iter().all(|&w| w == 0)
    }

    pub fn embed_u64(&self, v: u64) -> Elem<L> {
        let mut e = [0u64; L];
        e[0] = v;
        if self.m < 64 {
            e[0] &= (1u64 << self.m) - 1;
        }
        e
    }

    pub fn embed_big(&self, v: &BigUint) -> Result<Elem<L>, String> {
        if v.bits() > u64::from(self.m) {
            return Err("value has more bits than the field degree".into());
        }
        let mut e = [0u64; L];
        for (i, d) in v.to_u64_digits().into_iter().enumerate() {
            e[i] = d;
        }
        Ok(e)
    }

    pub fn to_big(&self, a: &Elem<L>) -> BigUint {
        a.iter()
            .rev()
            .fold(BigUint::zero(), |acc, &w| (acc << 64u32) | BigUint::from(w))
    }

    #[inline(always)]
    pub fn add(&self, a: &Elem<L>, b: &Elem<L>) -> Elem<L> {
        let mut r = *a;
        for i in 0..L {
            r[i] ^= b[i];
        }
        r
    }

    #[inline(always)]
    fn clmul(&self, a: u64, b: u64) -> (u64, u64) {
        #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
        if self.has_clmul {
            // SAFETY: `has_clmul` is set only when the CPU reports the
            // instruction.
            return unsafe { clmul_hw(a, b) };
        }
        clmul_portable(a, b)
    }

    fn reduce(&self, p: &mut [u64; BUF]) -> Elem<L> {
        let m = self.m;
        let n2 = 2 * self.limbs;
        loop {
            let top = shr(p, m, n2);
            if top.iter().all(|&w| w == 0) {
                break;
            }
            // Clear every bit at or above z^m.
            let (w, b) = ((m / 64) as usize, m % 64);
            if b == 0 {
                for x in p.iter_mut().skip(w) {
                    *x = 0;
                }
            } else {
                p[w] &= (1u64 << b) - 1;
                for x in p.iter_mut().skip(w + 1) {
                    *x = 0;
                }
            }
            for &t in &self.low {
                xor_shl(p, &top, t, n2);
            }
        }
        let mut out = [0u64; L];
        out.copy_from_slice(&p[..L]);
        out
    }

    pub fn mul(&self, a: &Elem<L>, b: &Elem<L>) -> Elem<L> {
        let mut p = [0u64; BUF];
        for i in 0..self.limbs {
            if a[i] == 0 {
                continue;
            }
            for j in 0..self.limbs {
                let (lo, hi) = self.clmul(a[i], b[j]);
                p[i + j] ^= lo;
                p[i + j + 1] ^= hi;
            }
        }
        self.reduce(&mut p)
    }

    pub fn sqr(&self, a: &Elem<L>) -> Elem<L> {
        self.mul(a, a)
    }

    fn sqr_k(&self, a: &Elem<L>, k: u32) -> Elem<L> {
        let mut r = *a;
        for _ in 0..k {
            r = self.sqr(&r);
        }
        r
    }

    /// `a⁻¹` by Itoh–Tsujii; `0` maps to `0`.
    pub fn inv(&self, a: &Elem<L>) -> Elem<L> {
        if self.is_zero(a) {
            return self.zero();
        }
        let e = self.m - 1;
        let bits = 32 - e.leading_zeros();
        let mut c = *a;
        let mut k = 1u32;
        for i in (0..bits - 1).rev() {
            let s = self.sqr_k(&c, k);
            c = self.mul(&c, &s);
            k *= 2;
            if (e >> i) & 1 == 1 {
                let s = self.sqr(&c);
                c = self.mul(&s, a);
                k += 1;
            }
        }
        debug_assert_eq!(k, e);
        self.sqr(&c)
    }

    /// The absolute trace, `0` or `1`.
    pub fn trace(&self, a: &Elem<L>) -> u64 {
        let mut t = *a;
        let mut acc = *a;
        for _ in 1..self.m {
            t = self.sqr(&t);
            acc = self.add(&acc, &t);
        }
        debug_assert!(acc[1..].iter().all(|&w| w == 0) && acc[0] <= 1);
        acc[0] & 1
    }

    /// A root `w` of `w² + w = c` for odd `m` and `Tr(c) = 0`, by the
    /// half-trace; `None` otherwise.
    pub fn solve_quadratic(&self, c: &Elem<L>) -> Option<Elem<L>> {
        if self.m.is_multiple_of(2) || self.trace(c) != 0 {
            return None;
        }
        let mut t = *c;
        let mut acc = *c;
        for _ in 0..(self.m - 1) / 2 {
            t = self.sqr(&self.sqr(&t));
            acc = self.add(&acc, &t);
        }
        debug_assert_eq!(self.add(&self.sqr(&acc), &acc), *c);
        Some(acc)
    }

    /// A uniform element.
    pub fn random(&self, rng: &mut StdRng) -> Elem<L> {
        let mut e = [0u64; L];
        for w in e.iter_mut().take(self.limbs) {
            *w = rng.gen();
        }
        let r = self.m % 64;
        if r != 0 {
            e[self.limbs - 1] &= (1u64 << r) - 1;
        }
        e
    }
}

/// An affine point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Pt<const L: usize> {
    pub x: Elem<L>,
    pub y: Elem<L>,
}

/// `y² + x·y = x³ + a·x² + b` with the order `n` of a prime-order subgroup
/// and the cofactor `h`, so the group has `h·n` points.
#[derive(Clone, Debug)]
pub struct BinCurve<const L: usize> {
    pub name: String,
    pub f: Gf2m<L>,
    pub a: Elem<L>,
    pub b: Elem<L>,
    pub order: BigUint,
    pub cofactor: u64,
    pub generator: Pt<L>,
}

impl<const L: usize> BinCurve<L> {
    /// A SECG binary curve from the standards database.
    pub fn from_secg(c: &SecgCurve) -> Result<Self, String> {
        let SecgField::Binary { degree, low_terms } = &c.field else {
            return Err(format!("{} is not a binary-field curve", c.name));
        };
        if degree % 2 == 0 {
            return Err("the half-trace needs an odd extension degree".into());
        }
        let f = Gf2m::new(*degree, low_terms)?;
        let curve = Self {
            name: c.name.clone(),
            a: f.embed_big(&c.a)?,
            b: f.embed_big(&c.b)?,
            order: c.order.clone(),
            cofactor: c.cofactor,
            generator: Pt {
                x: f.embed_big(&c.gx)?,
                y: f.embed_big(&c.gy)?,
            },
            f,
        };
        if !curve.on_curve(&curve.generator) {
            return Err(format!("{}: the generator is not on the curve", c.name));
        }
        Ok(curve)
    }

    pub fn on_curve(&self, p: &Pt<L>) -> bool {
        let f = &self.f;
        let lhs = f.add(&f.sqr(&p.y), &f.mul(&p.x, &p.y));
        let x2 = f.sqr(&p.x);
        let rhs = f.add(&f.add(&f.mul(&x2, &p.x), &f.mul(&self.a, &x2)), &self.b);
        lhs == rhs
    }

    pub fn negate(&self, p: &Pt<L>) -> Pt<L> {
        Pt {
            x: p.x,
            y: self.f.add(&p.x, &p.y),
        }
    }

    /// `P + Q`, or `None` for the point at infinity.
    pub fn add_points(&self, p: &Pt<L>, q: &Pt<L>) -> Option<Pt<L>> {
        let f = &self.f;
        if p.x == q.x {
            if p.y != q.y || f.is_zero(&p.x) {
                // Q = −P (or P of order two).
                return None;
            }
            // Doubling: λ = x + y/x.
            let lam = f.add(&p.x, &f.mul(&p.y, &f.inv(&p.x)));
            let x3 = f.add(&f.add(&f.sqr(&lam), &lam), &self.a);
            let y3 = f.add(&f.sqr(&p.x), &f.mul(&f.add(&lam, &f.one()), &x3));
            return Some(Pt { x: x3, y: y3 });
        }
        let lam = f.mul(&f.add(&p.y, &q.y), &f.inv(&f.add(&p.x, &q.x)));
        let x3 = f.add(
            &f.add(&f.add(&f.sqr(&lam), &lam), &f.add(&p.x, &q.x)),
            &self.a,
        );
        let y3 = f.add(&f.add(&f.mul(&lam, &f.add(&p.x, &x3)), &x3), &p.y);
        Some(Pt { x: x3, y: y3 })
    }

    /// `[k]P` for a small `k ≥ 1`; `None` for infinity.
    pub fn mul_small(&self, p: &Pt<L>, k: u64) -> Option<Pt<L>> {
        let mut acc: Option<Pt<L>> = None;
        for i in (0..64 - k.leading_zeros()).rev() {
            acc = match acc {
                Some(a) => self.add_points(&a, &a),
                None => None,
            };
            if (k >> i) & 1 == 1 {
                acc = match acc {
                    Some(a) => self.add_points(&a, p),
                    None => Some(*p),
                };
            }
        }
        acc
    }

    /// `[k]P` for any `k`; `None` for infinity.
    pub fn mul_big(&self, p: &Pt<L>, k: &BigUint) -> Option<Pt<L>> {
        let mut acc: Option<Pt<L>> = None;
        for i in (0..k.bits()).rev() {
            acc = match acc {
                Some(a) => self.add_points(&a, &a),
                None => None,
            };
            if k.bit(i) {
                acc = match acc {
                    Some(a) => self.add_points(&a, p),
                    None => Some(*p),
                };
            }
        }
        acc
    }

    /// The point with abscissa `x ≠ 0` and the chosen root, if any.
    /// `y = x·w` with `w² + w = x + a + b/x²`.
    fn lift_x(&self, x: &Elem<L>, second_root: bool) -> Option<Pt<L>> {
        let f = &self.f;
        if f.is_zero(x) {
            return None;
        }
        let x2inv = f.inv(&f.sqr(x));
        let c = f.add(&f.add(x, &self.a), &f.mul(&self.b, &x2inv));
        let mut w = f.solve_quadratic(&c)?;
        if second_root {
            w = f.add(&w, &f.one());
        }
        Some(Pt {
            x: *x,
            y: f.mul(x, &w),
        })
    }

    /// A uniform element of `⟨G⟩` other than infinity: a uniform curve
    /// point pushed through `[h]`.
    pub fn random_subgroup_point(&self, rng: &mut StdRng) -> Pt<L> {
        loop {
            let x = self.f.random(rng);
            let second = rng.gen::<bool>();
            if let Some(p) = self.lift_x(&x, second) {
                if let Some(q) = self.mul_small(&p, self.cofactor) {
                    return q;
                }
            }
        }
    }
}

/// The base: for each of the smallest abscissae `x = 1, 2, …` on the curve,
/// the point `[h]P` (infinity and repeated abscissae skipped), until `size`
/// are held.  Returns the base and the largest `x` examined.
pub fn small_x_base<const L: usize>(
    c: &BinCurve<L>,
    size: usize,
) -> Result<(Vec<Pt<L>>, u64), String> {
    let mut out: Vec<Pt<L>> = Vec::with_capacity(size);
    let mut seen: std::collections::HashSet<Elem<L>> = std::collections::HashSet::new();
    let limit = 64 * size as u64 + 1024;
    let mut x = 1u64;
    let mut largest = 0u64;
    while out.len() < size {
        if x > limit {
            return Err("ran out of small abscissae on this curve".into());
        }
        let xe = c.f.embed_u64(x);
        if let Some(p) = c.lift_x(&xe, false) {
            largest = x;
            if let Some(q) = c.mul_small(&p, c.cofactor) {
                if seen.insert(q.x) {
                    out.push(q);
                }
            }
        }
        x += 1;
    }
    Ok((out, largest))
}

/// The abscissae of `a ± b` for every pair of distinct base abscissae and of
/// `2a` for every `a`, sorted: `B²` keys with repeats where sums collide.
pub fn pair_sum_keys<const L: usize>(c: &BinCurve<L>, base: &[Pt<L>]) -> Vec<Elem<L>> {
    let n = base.len();
    let mut keys: Vec<Elem<L>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let f = &c.f;
            let mut out: Vec<Elem<L>> = Vec::with_capacity(2 * (n - i) + 1);
            let a = base[i];
            // Doubling: x(2a) = x² + b/x².
            let x2 = f.sqr(&a.x);
            out.push(f.add(&x2, &f.mul(&c.b, &f.inv(&x2))));
            let rest = &base[i + 1..];
            if !rest.is_empty() {
                let dens: Vec<Elem<L>> = rest.iter().map(|b| f.add(&a.x, &b.x)).collect();
                let mut prefix: Vec<Elem<L>> = Vec::with_capacity(dens.len());
                let mut acc = f.one();
                for d in &dens {
                    acc = f.mul(&acc, d);
                    prefix.push(acc);
                }
                let mut inv_acc = f.inv(&acc);
                let mut invs = vec![f.zero(); dens.len()];
                for k in (0..dens.len()).rev() {
                    let before = if k == 0 { f.one() } else { prefix[k - 1] };
                    invs[k] = f.mul(&inv_acc, &before);
                    inv_acc = f.mul(&inv_acc, &dens[k]);
                }
                for (k, b) in rest.iter().enumerate() {
                    let sx = f.add(&a.x, &b.x);
                    let base_y = f.add(&a.y, &b.y);
                    // a + b uses y_a + y_b; a − b uses y_a + y_b + x_b.
                    for ny in [base_y, f.add(&base_y, &b.x)] {
                        let lam = f.mul(&ny, &invs[k]);
                        let x3 = f.add(&f.add(&f.add(&f.sqr(&lam), &lam), &sx), &c.a);
                        out.push(x3);
                    }
                }
            }
            out
        })
        .flatten()
        .collect();
    keys.par_sort_unstable();
    keys
}

/// Whether `target` is `±a ± b` for base points.
fn recover<const L: usize>(
    c: &BinCurve<L>,
    base: &[Pt<L>],
    index: &HashMap<Elem<L>, usize>,
    target: &Pt<L>,
) -> Option<(bool, usize, usize)> {
    for (i, a) in base.iter().enumerate() {
        for negate in [false, true] {
            let a_signed = if negate { c.negate(a) } else { *a };
            let minus_a = c.negate(&a_signed);
            // Q = target − a_signed, so target = a_signed + Q.
            let q = match c.add_points(target, &minus_a) {
                Some(q) => q,
                None => continue,
            };
            if let Some(&j) = index.get(&q.x) {
                return Some((negate, i, j));
            }
        }
    }
    None
}

/// Measure the exact yield of the small-x base of `c`.
pub fn measure<const L: usize>(c: &BinCurve<L>, opts: &Options) -> Result<YieldReport, String> {
    if opts.size < 2 {
        return Err("a base needs at least two abscissae".into());
    }
    let f = &c.f;

    let t0 = Instant::now();
    let (base, largest_x) = small_x_base(c, opts.size)?;
    for p in &base {
        if !c.on_curve(p) {
            return Err("a base point is not on the curve".into());
        }
    }
    let base_ms = t0.elapsed().as_secs_f64() * 1000.0;

    let t1 = Instant::now();
    check_memory(opts.size, 8 * L)?;
    let keys = pair_sum_keys(c, &base);
    let total = keys.len() as u64;
    let distinct = 1 + keys.windows(2).filter(|w| w[0] != w[1]).count() as u64;
    let pair_ms = t1.elapsed().as_secs_f64() * 1000.0;

    let order_minus_one = &c.order - BigUint::one();
    let covered = BigUint::from(2 * distinct);
    let exact_log2 = log2_big(&covered) - log2_big(&order_minus_one);
    let model_log2 =
        log2_big(&BigUint::from(2 * (opts.size as u64).pow(2))) - log2_big(&order_minus_one);
    let exact_yield = 2.0f64.powf(exact_log2);
    let model_yield = 2.0f64.powf(model_log2);

    let index: HashMap<Elem<L>, usize> = base.iter().enumerate().map(|(i, p)| (p.x, i)).collect();

    let t2 = Instant::now();
    let chunk = 256u64;
    let chunks = opts.trials.div_ceil(chunk);
    let (hits, samples) = (0..chunks)
        .into_par_iter()
        .map(|ci| {
            let mut rng =
                StdRng::seed_from_u64(opts.seed ^ 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(ci + 1));
            let here = chunk.min(opts.trials - ci * chunk);
            let mut hits = 0u64;
            let mut samples: Vec<Pt<L>> = Vec::new();
            for _ in 0..here {
                let p = c.random_subgroup_point(&mut rng);
                if keys.binary_search(&p.x).is_ok() {
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

    let mut verified = 0u64;
    for p in samples.iter().take(8) {
        if recover(c, &base, &index, p).is_some() {
            verified += 1;
        }
    }

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
            b = c.negate(&b);
        }
        let Some(target) = c.add_points(&a, &b) else {
            continue;
        };
        planted += 1;
        if keys.binary_search(&target.x).is_ok() {
            detected += 1;
        }
        if recover(c, &base, &index, &target).is_some() {
            recovered += 1;
        }
    }
    let ctl_ms = t3.elapsed().as_secs_f64() * 1000.0;

    let trials = opts.trials;
    let order_bits = log2_big(&c.order);
    let trials_full = (opts.size as f64 + 1.0).log2() - exact_log2;
    Ok(YieldReport {
        curve: c.name.clone(),
        field_kind: "binary",
        cofactor: c.cofactor,
        field_bits: u64::from(f.m),
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
    })
}

/// The share of signed pair sums of raw small-x points (no `[h]`) that lie
/// in the prime-order subgroup, from `samples` random pairs; see
/// [`super::prime_small_x::raw_subgroup_fraction`].
pub fn raw_subgroup_fraction<const L: usize>(
    c: &BinCurve<L>,
    size: usize,
    samples: u64,
    seed: u64,
) -> Result<RawBase, String> {
    let mut raw_curve = c.clone();
    raw_curve.cofactor = 1;
    let (base, _) = small_x_base(&raw_curve, size)?;
    let hits: u64 = (0..samples)
        .into_par_iter()
        .map(|i| {
            let mut rng = StdRng::seed_from_u64(seed ^ 0x5241_5742_4153_45u64.wrapping_mul(i + 1));
            let a = base[rng.gen_range(0..base.len())];
            let mut b = base[rng.gen_range(0..base.len())];
            if rng.gen::<bool>() {
                b = c.negate(&b);
            }
            match c.add_points(&a, &b) {
                None => 1,
                Some(sum) => u64::from(c.mul_big(&sum, &c.order).is_none()),
            }
        })
        .sum();
    Ok(raw_base_summary(c.cofactor, &c.order, size, hits, samples))
}

// ── Toy curves ─────────────────────────────────────────────────────

/// Irreducible polynomials of prime degree with an odd number of terms, for
/// toy curves: `(degree, low terms)`.  A prime degree is what makes the
/// check `z^(2^m) = z` prove irreducibility (a factor has degree 1 or `m`,
/// and `z`, `z + 1` do not divide a polynomial with a constant term and an
/// odd number of terms).
const TOY_FIELDS: &[(u32, &[u32])] = &[
    (11, &[0, 2]),
    (13, &[0, 1, 3, 4]),
    (17, &[0, 3]),
    (19, &[0, 1, 2, 5]),
];

/// A toy binary curve of cofactor 2 (`a = 1`): see [`toy_curve_with_cofactor`].
pub fn toy_curve<const L: usize>(degree: u32, seed: u64) -> Result<BinCurve<L>, String> {
    toy_curve_with_cofactor(degree, seed, 2)
}

/// A toy binary curve `y² + xy = x³ + a·x² + b` over a small odd-degree field
/// whose group order is `h·n` with `n` prime, counted exactly by walking the
/// field, so cofactor handling and the yield can be checked against brute
/// force.  `h = 2` uses `a = 1` and `h = 4` uses `a = 0` (the Koblitz
/// shape: the order of such a curve is a multiple of 4).  `degree` is one of
/// 11, 13, 17, 19.
pub fn toy_curve_with_cofactor<const L: usize>(
    degree: u32,
    seed: u64,
    h: u64,
) -> Result<BinCurve<L>, String> {
    if h != 2 && h != 4 {
        return Err("toy binary cofactors are 2 and 4".into());
    }
    let &(_, low) = TOY_FIELDS
        .iter()
        .find(|(d, _)| *d == degree)
        .ok_or("toy binary degrees are 11, 13, 17 or 19")?;
    let f = Gf2m::new(degree, low)?;
    // Irreducible: z^(2^m) = z, and see the table note on prime degrees.
    let z = f.embed_u64(2);
    let mut t = z;
    for _ in 0..degree {
        t = f.sqr(&t);
    }
    if t != z {
        return Err("the toy modulus is not irreducible".into());
    }
    let mut rng = StdRng::seed_from_u64(seed ^ u64::from(degree) ^ (h << 40));
    for _ in 0..4000 {
        let b = f.random(&mut rng);
        if f.is_zero(&b) {
            continue;
        }
        let a = if h == 2 { f.one() } else { f.zero() };
        let curve = BinCurve {
            name: String::new(),
            f: f.clone(),
            a,
            b,
            order: BigUint::zero(),
            cofactor: h,
            generator: Pt {
                x: f.zero(),
                y: f.zero(),
            },
        };
        // Count points: x = 0 gives one point (y = √b); other x give 2 or 0.
        let mut points = 1u64 + 1;
        for x in 1..(1u64 << degree) {
            let xe = f.embed_u64(x);
            if curve.lift_x(&xe, false).is_some() {
                points += 2;
            }
        }
        if points.is_multiple_of(h) && is_prime(points / h) && points / h > 7 {
            let n = points / h;
            return Ok(BinCurve {
                name: format!(
                    "toy y²+xy=x³+{}x²+b over F_2^{degree}, #E = {h}·{n}",
                    u8::from(h == 2)
                ),
                order: BigUint::from(n),
                ..curve
            });
        }
    }
    Err("no toy curve of that cofactor found".into())
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    let mut d = 2;
    while d * d <= n {
        if n.is_multiple_of(d) {
            return false;
        }
        d += 1;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement};
    use crate::cryptanalysis::secg_standards::secg;
    use std::collections::HashSet;

    fn secg_curve<const L: usize>(name: &str) -> BinCurve<L> {
        BinCurve::from_secg(secg(name).expect(name)).expect(name)
    }

    fn field_vs_repository<const L: usize>(name: &str, irr: crate::binary_ecc::IrreduciblePoly) {
        let c = secg_curve::<L>(name);
        let m = irr.degree;
        let mut rng = StdRng::seed_from_u64(9);
        for _ in 0..20 {
            let (a, b) = (c.f.random(&mut rng), c.f.random(&mut rng));
            let (fa, fb) = (
                F2mElement::from_biguint(&c.f.to_big(&a), m),
                F2mElement::from_biguint(&c.f.to_big(&b), m),
            );
            let want = fa.mul(&fb, &irr).to_biguint();
            assert_eq!(c.f.to_big(&c.f.mul(&a, &b)), want, "{name} mul");
            assert_eq!(
                c.f.to_big(&c.f.mul(&a, &b)),
                c.f.to_big(&c.f.with_portable_kernel().mul(&a, &b)),
                "{name} kernels"
            );
            assert_eq!(c.f.to_big(&c.f.sqr(&a)), fa.square(&irr).to_biguint());
            if !c.f.is_zero(&a) {
                assert_eq!(c.f.mul(&a, &c.f.inv(&a)), c.f.one(), "{name} inverse");
            }
        }
    }

    #[test]
    fn field_matches_the_repository_f2m_at_every_limb_count() {
        use crate::binary_ecc::IrreduciblePoly as P;
        field_vs_repository::<2>("sect113r1", P::deg_113());
        field_vs_repository::<3>("sect163r2", P::deg_163());
        field_vs_repository::<4>("sect233r1", P::deg_233());
        field_vs_repository::<5>("sect283r1", P::deg_283());
        field_vs_repository::<7>("sect409r1", P::deg_409());
        field_vs_repository::<9>("sect571r1", P::deg_571());
    }

    fn generator_has_the_published_order<const L: usize>(name: &str) {
        let c = secg_curve::<L>(name);
        assert!(c.on_curve(&c.generator), "{name}");
        // n·G = infinity, (n − 1)·G = −G.
        assert!(c.mul_big(&c.generator, &c.order).is_none(), "{name} n·G");
        let nm1 = c.mul_big(&c.generator, &(&c.order - 1u8)).expect(name);
        assert_eq!(nm1, c.negate(&c.generator), "{name} (n−1)·G");
    }

    #[test]
    fn the_published_generators_are_on_the_curve_with_the_published_order() {
        generator_has_the_published_order::<2>("sect113r1");
        generator_has_the_published_order::<2>("sect113r2");
        generator_has_the_published_order::<3>("sect163r1");
        generator_has_the_published_order::<3>("sect163r2");
        generator_has_the_published_order::<4>("sect233k1");
        generator_has_the_published_order::<5>("sect283k1");
        generator_has_the_published_order::<7>("sect409k1");
        generator_has_the_published_order::<9>("sect571k1");
    }

    #[test]
    fn group_law_agrees_with_the_repository_curve_on_sect113r1() {
        let c = secg_curve::<2>("sect113r1");
        let hand = BinaryCurve::sect113r1();
        let g2 = c.add_points(&c.generator, &c.generator).unwrap();
        let g3 = c.add_points(&g2, &c.generator).unwrap();
        let want =
            crate::binary_ecc::curve::scalar_mul(&hand, &hand.generator, &BigUint::from(3u8));
        let BinaryPoint::Affine { x, y } = want else {
            panic!("3G is finite")
        };
        assert_eq!(c.f.to_big(&g3.x), x.to_biguint());
        assert_eq!(c.f.to_big(&g3.y), y.to_biguint());
    }

    #[test]
    fn small_x_base_lies_in_the_prime_order_subgroup() {
        for name in ["sect113r1", "sect163r2"] {
            let c = secg_curve::<3>(name);
            let (base, _) = small_x_base(&c, 6).unwrap();
            assert_eq!(base.len(), 6);
            for p in &base {
                assert!(c.on_curve(p), "{name}");
                assert!(c.mul_big(p, &c.order).is_none(), "{name}: n·P != O");
            }
            let xs: HashSet<Elem<3>> = base.iter().map(|p| p.x).collect();
            assert_eq!(xs.len(), 6);
        }
    }

    #[test]
    fn toy_order_is_the_counted_one() {
        let c = toy_curve::<1>(11, 1).expect("toy");
        // Every subgroup element has order dividing n.
        let (base, _) = small_x_base(&c, 8).unwrap();
        for p in &base {
            assert!(c.mul_big(p, &c.order).is_none());
        }
    }

    #[test]
    fn exact_yield_equals_brute_force_coverage_on_a_toy_curve() {
        let c = toy_curve::<1>(13, 3).expect("toy");
        let (base, _) = small_x_base(&c, 40).unwrap();
        let mut signed: Vec<Pt<1>> = Vec::new();
        for p in &base {
            signed.push(*p);
            signed.push(c.negate(p));
        }
        let mut covered: HashSet<(Elem<1>, Elem<1>)> = HashSet::new();
        for a in &signed {
            for b in &signed {
                if let Some(s) = c.add_points(a, b) {
                    covered.insert((s.x, s.y));
                }
            }
        }
        let mut keys = pair_sum_keys(&c, &base);
        keys.dedup();
        assert_eq!(covered.len(), 2 * keys.len());
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
    fn monte_carlo_rate_matches_the_exact_yield_on_a_toy_curve() {
        let c = toy_curve::<1>(17, 7).expect("toy");
        let rep = measure(
            &c,
            &Options {
                size: 200,
                trials: 300_000,
                controls: 8,
                seed: 11,
            },
        )
        .expect("measure");
        assert!(rep.collisions > 0 || rep.mc_hits > 200, "{rep:?}");
        assert!(rep.mc_hits > 200, "hits {}", rep.mc_hits);
        let [lo, hi] = rep.mc_ci95;
        assert!(
            lo <= rep.exact_yield && rep.exact_yield <= hi,
            "exact {} outside [{lo}, {hi}]",
            rep.exact_yield
        );
        assert_eq!(rep.mc_hits_verified, 8);
    }

    #[test]
    fn sect113r1_has_a_null_monte_carlo_and_working_controls() {
        let c = secg_curve::<2>("sect113r1");
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
        // 2·64² / (n − 1), n ≈ 2^112.
        assert!(
            (rep.exact_yield_log2 + 99.0).abs() < 1.0,
            "{}",
            rep.exact_yield_log2
        );
    }

    #[test]
    fn toy_cofactor_four_binary_curve_matches_monte_carlo() {
        let c = toy_curve_with_cofactor::<1>(17, 5, 4).expect("toy");
        assert_eq!(c.cofactor, 4);
        let rep = measure(
            &c,
            &Options {
                size: 150,
                trials: 300_000,
                controls: 8,
                seed: 13,
            },
        )
        .expect("measure");
        assert!(rep.mc_hits > 200, "hits {}", rep.mc_hits);
        let [lo, hi] = rep.mc_ci95;
        assert!(
            lo <= rep.exact_yield && rep.exact_yield <= hi,
            "exact {} outside [{lo}, {hi}]",
            rep.exact_yield
        );
        assert_eq!(rep.controls_recovered, rep.controls);
        assert_eq!(rep.mc_hits_verified, 8);
    }

    #[test]
    fn raw_base_fraction_on_toy_binary_curves() {
        for h in [2u64, 4] {
            let c = toy_curve_with_cofactor::<1>(17, 5, h).expect("toy");
            let r = raw_subgroup_fraction(&c, 60, 20000, 3).expect("raw");
            // Whatever the class structure, the fraction is a probability and
            // not above the subgroup base's share.
            assert!(r.fraction > 0.0 && r.fraction < 1.0, "h={h}: {r:?}");
            eprintln!(
                "h={h}: raw fraction {:.4} (1/h = {:.4})",
                r.fraction, r.one_over_h
            );
        }
    }
}
