//! **Binary fields of degree up to 127 in two words** (`ic` tool
//! programme, plan §9, B3; brought into the tree from
//! `research/ic_tool_program/track-b`, branch `b3-int` at `2b7306bb`,
//! for `ecbench`'s m = 83 gate group, `koblitz_wide`).
//!
//! [`crate::cryptanalysis::semaev_decomp::Gf2`] holds a field element in
//! one `u64` and stops at degree 63.  This holds it in a `u128` and
//! reaches 127: AGENTS.md §8a's `m = 83` gate and everything between.
//!
//! - **Multiplication** is four carry-less 64-bit products, PCLMULQDQ
//!   (or PMULL) where the CPU has it and a portable loop otherwise.
//! - **Reduction** folds when the modulus's low terms `f − zⁿ` fit one
//!   word, as they do for every trinomial and pentanomial the standards
//!   use: `H·zⁿ ≡ H·(f − zⁿ)`, one carry-less product per word of `H`,
//!   repeated while anything is left above degree `n` (three times on
//!   the gate's modulus, once on a trinomial with a low middle term).
//!   Any other modulus takes polynomial Barrett: with `μ = ⌊z^{2n}/f⌋`
//!   the quotient of a product `P` by `f` is `⌊⌊P/zⁿ⌋·μ/zⁿ⌋` exactly,
//!   because GF(2)[z] has no carries.  Two more products and no table.
//! - **Inversion** is Itoh–Tsujii: `a⁻¹ = (a^{2^{n−1}−1})²`, `n − 1`
//!   squarings and about `log₂ n` multiplications.
//! - **The trace** is one AND and a parity, against a mask of the basis's
//!   traces.
//!
//! **Kernels.**  [`Kernel`] is the arithmetic a hot loop is compiled
//! against.  On x86-64 [`Xmm`] keeps elements in vector registers from
//! product to reduced result — the carry-less products, the shifts by
//! `n` and the folds never pass through the general registers — and a
//! loop generic over the kernel, called from a function compiled with
//! the features on, inlines every operation.  [`Generic`] is the same
//! arithmetic on `u128` with any [`Clmul`]: the portable loop anywhere,
//! PMULL on Arm64.  The scalar methods ([`GfWide::mul`] and the rest)
//! dispatch to the best kernel one operation at a time.
//!
//! The tests check every operation against the general big-integer field
//! (`binary_ecc::F2mElement`) at degrees 31 to 127 and on the gate's
//! modulus, and the kernels against each other.
use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use num_bigint::BigUint;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// The largest degree two words hold with room for the product's
/// reduction.
pub const MAX_WIDE_DEGREE: u32 = 127;

/// `GF(2)[z] / (f)` for `deg f ≤ 127`, elements below `2^n` in a `u128`.
#[derive(Clone, Debug)]
pub struct GfWide {
    pub n: u32,
    /// `f − zⁿ`: the modulus without its leading term.
    low: u128,
    /// `⌊z^{2n} / f⌋`, of degree `n`.
    mu: u128,
    mask: u128,
    /// Bit `i` set when `Tr(z^i) = 1`.
    trace_mask: u128,
    /// How many folds clear a product when `f − zⁿ` fits one word and
    /// that is cheaper than Barrett; `0` selects Barrett.
    folds: u32,
    #[cfg(target_arch = "x86_64")]
    xmm: Option<Xmm>,
    #[cfg(target_arch = "aarch64")]
    pmull: Option<Pmull>,
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "aes")]
#[inline]
unsafe fn clmul64_hw(a: u64, b: u64) -> u128 {
    std::arch::aarch64::vmull_p64(a, b)
}

/// A 64-bit carry-less product, for the [`Generic`] kernel.
pub trait Clmul: Copy {
    fn clmul(self, a: u64, b: u64) -> u128;
}

/// The portable product, on any CPU.
#[derive(Clone, Copy, Debug)]
pub struct Portable;

impl Clmul for Portable {
    #[inline(always)]
    fn clmul(self, a: u64, b: u64) -> u128 {
        clmul64_portable(a, b)
    }
}

/// PMULL on Arm64.  Only [`Pmull::detect`] makes one, after checking the
/// CPU, so holding one is the proof the product is safe to issue.
#[cfg(target_arch = "aarch64")]
#[derive(Clone, Copy, Debug)]
pub struct Pmull(());

#[cfg(target_arch = "aarch64")]
impl Pmull {
    pub fn detect() -> Option<Self> {
        std::arch::is_aarch64_feature_detected!("aes").then_some(Self(()))
    }
}

#[cfg(target_arch = "aarch64")]
impl Clmul for Pmull {
    #[inline(always)]
    fn clmul(self, a: u64, b: u64) -> u128 {
        // SAFETY: a `Pmull` exists only where `detect` found the feature.
        unsafe { clmul64_hw(a, b) }
    }
}

#[inline]
fn clmul64_portable(a: u64, b: u64) -> u128 {
    let mut out = 0u128;
    let mut b = b;
    let a = a as u128;
    let mut i = 0;
    while b != 0 {
        if b & 1 == 1 {
            out ^= a << i;
        }
        b >>= 1;
        i += 1;
    }
    out
}

/// The full 256-bit carry-less product of two 128-bit polynomials, as
/// `(low, high)` halves, by `clmul64` four times.
#[inline(always)]
fn clmul128(a: u128, b: u128, clmul64: impl Fn(u64, u64) -> u128) -> (u128, u128) {
    let (a0, a1) = (a as u64, (a >> 64) as u64);
    let (b0, b1) = (b as u64, (b >> 64) as u64);
    let lo = clmul64(a0, b0);
    let hi = clmul64(a1, b1);
    let mid = clmul64(a0, b1) ^ clmul64(a1, b0);
    (lo ^ (mid << 64), hi ^ (mid >> 64))
}

/// `(lo, hi) >> s` for a 256-bit value, `0 < s < 128`.
#[inline(always)]
fn shr256(lo: u128, hi: u128, s: u32) -> u128 {
    (lo >> s) | (hi << (128 - s))
}

/// How many folds reduce a product modulo `zⁿ + (f − zⁿ)` with
/// `d = deg(f − zⁿ)`, or `0` when `d ≥ 64` or more than five would be
/// needed and Barrett is cheaper.  The first fold leaves at most degree
/// `d − 2` above `zⁿ`, and each later one lowers that by `n − d`.
fn fold_count(n: u32, d: u32) -> u32 {
    if d >= 64 || d >= n {
        return 0;
    }
    let mut folds = 1;
    let mut above = i64::from(d) - 2;
    while above >= 0 {
        folds += 1;
        above += i64::from(d) - i64::from(n);
        if folds > 5 {
            return 0;
        }
    }
    folds
}

/// Polynomial division over GF(2) of a 256-bit dividend, `n + 1`-bit
/// divisor: the quotient's low 128 bits.  Construction only.
fn poly_div_z2n(n: u32, f: u128) -> u128 {
    // z^{2n} as (lo, hi), reduced by long division.
    let (mut lo, mut hi) = if 2 * n >= 128 {
        (0u128, 1u128 << (2 * n - 128))
    } else {
        (1u128 << (2 * n), 0)
    };
    let mut q = 0u128;
    for shift in (0..=n).rev() {
        let bit = n + shift;
        let set = if bit >= 128 {
            (hi >> (bit - 128)) & 1
        } else {
            (lo >> bit) & 1
        };
        if set == 1 {
            q |= 1u128 << shift;
            // Subtract f · z^shift.
            let (flo, fhi) = if shift == 0 {
                (f, 0)
            } else {
                (f << shift, f >> (128 - shift))
            };
            let fhi = if shift == 0 { 0 } else { fhi };
            lo ^= flo;
            hi ^= fhi;
        }
    }
    q
}

// ── Kernels ─────────────────────────────────────────────────────────

/// Field arithmetic as a hot loop is compiled against it: elements in
/// the kernel's own form `E`, converted where the loop reads and writes
/// memory.  Every method inlines, so a loop generic over a kernel and
/// called from a function compiled with the kernel's CPU features
/// compiles to straight-line arithmetic.
pub trait Kernel: Copy {
    type E: Copy;
    fn degree(self) -> u32;
    fn load(self, v: u128) -> Self::E;
    fn store(self, e: Self::E) -> u128;
    fn xor(self, a: Self::E, b: Self::E) -> Self::E;
    fn mul(self, a: Self::E, b: Self::E) -> Self::E;
    fn sqr(self, a: Self::E) -> Self::E;

    /// A linear map applied by byte tables: the XOR of `tables[i][b_i]`
    /// over the bytes `b_i` of `x`.  Indexed by the bytes themselves, so
    /// nothing is shifted.
    #[inline(always)]
    fn apply(self, tables: &[[u128; 256]], x: u128) -> u128 {
        let bytes = x.to_le_bytes();
        let mut acc = 0u128;
        for (table, &b) in tables.iter().zip(bytes.iter()) {
            acc ^= table[usize::from(b)];
        }
        acc
    }

    /// `a⁻¹` by Itoh–Tsujii, for `a ≠ 0`: `β_k = a^{2^k − 1}` built along
    /// the binary expansion of `n − 1`, then `a⁻¹ = β_{n−1}²`.
    #[inline(always)]
    fn inv(self, a: Self::E) -> Self::E {
        let m = self.degree() - 1;
        let mut beta = a;
        let mut k = 1u32;
        for bit in (0..31 - m.leading_zeros()).rev() {
            // β_{2k} = β_k^{2^k} · β_k.
            let mut s = beta;
            for _ in 0..k {
                s = self.sqr(s);
            }
            beta = self.mul(s, beta);
            k *= 2;
            if (m >> bit) & 1 == 1 {
                // β_{k+1} = β_k² · a.
                beta = self.mul(self.sqr(beta), a);
                k += 1;
            }
        }
        self.sqr(beta)
    }
}

/// The arithmetic on `u128` with a given carry-less product.
#[derive(Clone, Copy, Debug)]
pub struct Generic<'a, C: Clmul> {
    field: &'a GfWide,
    clmul: C,
}

impl<C: Clmul> Kernel for Generic<'_, C> {
    type E = u128;
    #[inline(always)]
    fn degree(self) -> u32 {
        self.field.n
    }
    #[inline(always)]
    fn load(self, v: u128) -> u128 {
        v
    }
    #[inline(always)]
    fn store(self, e: u128) -> u128 {
        e
    }
    #[inline(always)]
    fn xor(self, a: u128, b: u128) -> u128 {
        a ^ b
    }
    #[inline(always)]
    fn mul(self, a: u128, b: u128) -> u128 {
        self.field.mul_by(self.clmul, a, b)
    }
    #[inline(always)]
    fn sqr(self, a: u128) -> u128 {
        self.field.sqr_by(self.clmul, a)
    }
}

/// The arithmetic in SSE registers with PCLMULQDQ.  Only
/// [`GfWide::xmm`] makes one, on a CPU it has checked for PCLMULQDQ,
/// SSSE3 and POPCNT, so holding one is the proof they are safe to issue;
/// a function that instantiates a loop with it enables those features.
#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy, Debug)]
pub struct Xmm {
    n: u32,
    folds: u32,
    /// `n ≥ 64`: a shift by `n` crosses a whole word first.
    wide: bool,
    low: __m128i,
    mu: __m128i,
    mask: __m128i,
    /// The shift counts within a word: `n mod 64`-ish and its complement,
    /// so that `>> n` is two word shifts and an OR.
    s1: __m128i,
    s2: __m128i,
}

#[cfg(target_arch = "x86_64")]
impl Xmm {
    #[target_feature(enable = "pclmulqdq,ssse3")]
    unsafe fn new(field: &GfWide) -> Self {
        let n = field.n;
        let wide = n >= 64;
        let s1 = if wide { n - 64 } else { n };
        Self {
            n,
            folds: field.folds,
            wide,
            low: std::mem::transmute::<u128, __m128i>(field.low),
            mu: std::mem::transmute::<u128, __m128i>(field.mu),
            mask: std::mem::transmute::<u128, __m128i>(field.mask),
            s1: _mm_cvtsi64_si128(i64::from(s1)),
            s2: _mm_cvtsi64_si128(i64::from(64 - s1)),
        }
    }

    /// Whether a shift by `n` crosses a whole word, and the fold count:
    /// the two things [`XmmShaped`] fixes at compile time.
    pub fn shape(&self) -> (bool, u32) {
        (self.wide, self.folds)
    }

    /// `(lo, hi) >> n`, the low 128 bits.  `W` fixes the shape: 1 for
    /// `n < 64`, 2 for `n ≥ 64`, 0 to read it at run time.
    #[target_feature(enable = "pclmulqdq,ssse3")]
    #[inline]
    unsafe fn shr_n<const W: u8>(self, lo: __m128i, hi: __m128i) -> __m128i {
        // [lo₁, hi₀]: the words a shift by n ≥ 64 starts from.
        let mid = _mm_alignr_epi8::<8>(hi, lo);
        let wide = match W {
            1 => false,
            2 => true,
            _ => self.wide,
        };
        let (x, y) = if wide { (mid, hi) } else { (lo, mid) };
        _mm_or_si128(_mm_srl_epi64(x, self.s1), _mm_sll_epi64(y, self.s2))
    }

    #[target_feature(enable = "pclmulqdq,ssse3")]
    #[inline]
    unsafe fn product(a: __m128i, b: __m128i) -> (__m128i, __m128i) {
        let lo = _mm_clmulepi64_si128::<0x00>(a, b);
        let hi = _mm_clmulepi64_si128::<0x11>(a, b);
        let mid = _mm_xor_si128(
            _mm_clmulepi64_si128::<0x01>(a, b),
            _mm_clmulepi64_si128::<0x10>(a, b),
        );
        (
            _mm_xor_si128(lo, _mm_slli_si128::<8>(mid)),
            _mm_xor_si128(hi, _mm_srli_si128::<8>(mid)),
        )
    }

    /// `P mod f`.  `F` fixes the fold count (1 to 5), 0 to read it —
    /// and Barrett when it is 0 — at run time.
    #[target_feature(enable = "pclmulqdq,ssse3")]
    #[inline]
    unsafe fn reduce<const W: u8, const F: u32>(self, lo: __m128i, hi: __m128i) -> __m128i {
        let h = self.shr_n::<W>(lo, hi);
        let folds = if F == 0 { self.folds } else { F };
        if folds != 0 {
            // H·low over both words of H, then single words.
            let p0 = _mm_clmulepi64_si128::<0x00>(h, self.low);
            let p1 = _mm_clmulepi64_si128::<0x01>(h, self.low);
            let wlo = _mm_xor_si128(p0, _mm_slli_si128::<8>(p1));
            let whi = _mm_srli_si128::<8>(p1);
            let mut acc = _mm_and_si128(_mm_xor_si128(lo, wlo), self.mask);
            let mut h = self.shr_n::<W>(wlo, whi);
            for _ in 1..folds {
                let w = _mm_clmulepi64_si128::<0x00>(h, self.low);
                acc = _mm_xor_si128(acc, _mm_and_si128(w, self.mask));
                h = self.shr_n::<W>(w, _mm_setzero_si128());
            }
            acc
        } else {
            let (ql, qh) = Self::product(h, self.mu);
            let q = self.shr_n::<W>(ql, qh);
            let (rl, _) = Self::product(q, self.low);
            _mm_and_si128(_mm_xor_si128(lo, rl), self.mask)
        }
    }

    #[target_feature(enable = "pclmulqdq,ssse3")]
    #[inline]
    unsafe fn mul_raw<const W: u8, const F: u32>(self, a: __m128i, b: __m128i) -> __m128i {
        let (lo, hi) = Self::product(a, b);
        self.reduce::<W, F>(lo, hi)
    }

    #[target_feature(enable = "pclmulqdq,ssse3")]
    #[inline]
    unsafe fn sqr_raw<const W: u8, const F: u32>(self, a: __m128i) -> __m128i {
        let lo = _mm_clmulepi64_si128::<0x00>(a, a);
        let hi = _mm_clmulepi64_si128::<0x11>(a, a);
        self.reduce::<W, F>(lo, hi)
    }

    #[target_feature(enable = "pclmulqdq,ssse3")]
    #[inline]
    unsafe fn apply_raw(tables: &[[u128; 256]], x: u128) -> u128 {
        let bytes = x.to_le_bytes();
        let mut acc = _mm_setzero_si128();
        for (table, &b) in tables.iter().zip(bytes.iter()) {
            // `b < 256`, so the entry is inside the table.
            let entry = _mm_loadu_si128(table.as_ptr().add(usize::from(b)).cast());
            acc = _mm_xor_si128(acc, entry);
        }
        std::mem::transmute::<__m128i, u128>(acc)
    }

    #[target_feature(enable = "pclmulqdq,ssse3")]
    #[inline]
    unsafe fn xor_raw(a: __m128i, b: __m128i) -> __m128i {
        _mm_xor_si128(a, b)
    }

    #[target_feature(enable = "pclmulqdq,ssse3,popcnt")]
    unsafe fn mul_u128(self, a: u128, b: u128) -> u128 {
        self.store(self.mul(self.load(a), self.load(b)))
    }

    #[target_feature(enable = "pclmulqdq,ssse3,popcnt")]
    unsafe fn sqr_u128(self, a: u128) -> u128 {
        self.store(self.sqr(self.load(a)))
    }

    #[target_feature(enable = "pclmulqdq,ssse3,popcnt")]
    unsafe fn inv_u128(self, a: u128) -> u128 {
        self.store(self.inv(self.load(a)))
    }
}

/// [`Xmm`] with the modulus's shape fixed at compile time: `W` whether a
/// shift by `n` crosses a whole word (1 no, 2 yes) and `F` the fold
/// count (1 to 5), so that neither costs a branch or a loop per
/// operation.  `XmmShaped<0, 0>` reads both at run time and takes any
/// modulus, Barrett's included; [`Xmm`] itself is that kernel.
#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy, Debug)]
pub struct XmmShaped<const W: u8, const F: u32>(pub Xmm);

#[cfg(target_arch = "x86_64")]
impl<const W: u8, const F: u32> Kernel for XmmShaped<W, F> {
    type E = __m128i;
    #[inline(always)]
    fn degree(self) -> u32 {
        self.0.n
    }
    #[inline(always)]
    fn load(self, v: u128) -> __m128i {
        // SAFETY: both are sixteen plain bytes, the low word first.
        unsafe { std::mem::transmute::<u128, __m128i>(v) }
    }
    #[inline(always)]
    fn store(self, e: __m128i) -> u128 {
        // SAFETY: as in `load`.
        unsafe { std::mem::transmute::<__m128i, u128>(e) }
    }
    #[inline(always)]
    fn xor(self, a: __m128i, b: __m128i) -> __m128i {
        // SAFETY: an `Xmm` exists only where the features were found.
        unsafe { Xmm::xor_raw(a, b) }
    }
    #[inline(always)]
    fn mul(self, a: __m128i, b: __m128i) -> __m128i {
        // SAFETY: as in `xor`.
        unsafe { self.0.mul_raw::<W, F>(a, b) }
    }
    #[inline(always)]
    fn sqr(self, a: __m128i) -> __m128i {
        // SAFETY: as in `xor`.
        unsafe { self.0.sqr_raw::<W, F>(a) }
    }
    #[inline(always)]
    fn apply(self, tables: &[[u128; 256]], x: u128) -> u128 {
        // SAFETY: as in `xor`.
        unsafe { Xmm::apply_raw(tables, x) }
    }
}

#[cfg(target_arch = "x86_64")]
impl Kernel for Xmm {
    type E = __m128i;
    #[inline(always)]
    fn degree(self) -> u32 {
        self.n
    }
    #[inline(always)]
    fn load(self, v: u128) -> __m128i {
        XmmShaped::<0, 0>(self).load(v)
    }
    #[inline(always)]
    fn store(self, e: __m128i) -> u128 {
        XmmShaped::<0, 0>(self).store(e)
    }
    #[inline(always)]
    fn xor(self, a: __m128i, b: __m128i) -> __m128i {
        XmmShaped::<0, 0>(self).xor(a, b)
    }
    #[inline(always)]
    fn mul(self, a: __m128i, b: __m128i) -> __m128i {
        XmmShaped::<0, 0>(self).mul(a, b)
    }
    #[inline(always)]
    fn sqr(self, a: __m128i) -> __m128i {
        XmmShaped::<0, 0>(self).sqr(a)
    }
    #[inline(always)]
    fn apply(self, tables: &[[u128; 256]], x: u128) -> u128 {
        XmmShaped::<0, 0>(self).apply(tables, x)
    }
}

impl GfWide {
    /// The field `GF(2)[z]/(f)` for the given modulus, which the caller
    /// has found irreducible; `None` past degree 127 or below 2.
    pub fn new(irr: &IrreduciblePoly) -> Option<Self> {
        let n = irr.degree;
        if !(2..=MAX_WIDE_DEGREE).contains(&n) {
            return None;
        }
        let low = irr
            .low_terms
            .iter()
            .fold(0u128, |acc, &t| acc | (1u128 << t));
        let f = low | (1u128 << n);
        let mut field = Self {
            n,
            low,
            folds: fold_count(n, 127 - low.leading_zeros()),
            mu: poly_div_z2n(n, f),
            mask: (1u128 << n) - 1,
            trace_mask: 0,
            #[cfg(target_arch = "x86_64")]
            xmm: None,
            #[cfg(target_arch = "aarch64")]
            pmull: Pmull::detect(),
        };
        #[cfg(target_arch = "x86_64")]
        if std::arch::is_x86_feature_detected!("pclmulqdq")
            && std::arch::is_x86_feature_detected!("ssse3")
            && std::arch::is_x86_feature_detected!("popcnt")
        {
            // SAFETY: the features were found just above.
            field.xmm = Some(unsafe { Xmm::new(&field) });
        }
        field.trace_mask = (0..n).fold(0u128, |acc, i| {
            acc | (u128::from(field.trace_slow(1u128 << i)) << i)
        });
        Some(field)
    }

    /// Which carry-less product the scalar methods use.
    pub fn kernel_name(&self) -> &'static str {
        #[cfg(target_arch = "x86_64")]
        if self.xmm.is_some() {
            return "pclmulqdq";
        }
        #[cfg(target_arch = "aarch64")]
        if self.pmull.is_some() {
            return "pmull";
        }
        "portable"
    }

    /// The vector kernel, on a CPU that has its features.
    #[cfg(target_arch = "x86_64")]
    pub fn xmm(&self) -> Option<Xmm> {
        self.xmm
    }

    /// The `u128` kernel with PMULL, on a CPU that has it.
    #[cfg(target_arch = "aarch64")]
    pub fn pmull(&self) -> Option<Generic<'_, Pmull>> {
        self.pmull.map(|clmul| Generic { field: self, clmul })
    }

    /// The `u128` kernel with the portable product.
    pub fn portable(&self) -> Generic<'_, Portable> {
        Generic {
            field: self,
            clmul: Portable,
        }
    }

    /// `P mod f` by folding, `f − zⁿ` in one word: `H·zⁿ ≡ H·(f − zⁿ)`.
    /// The first fold multiplies both words of `H`; every later `H` is
    /// below `z^{d−1}` for `d = deg(f − zⁿ) < 64`, one product each.
    #[inline(always)]
    fn reduce_fold(&self, lo: u128, hi: u128, clmul: impl Fn(u64, u64) -> u128) -> u128 {
        let n = self.n;
        let low = self.low as u64;
        let h = shr256(lo, hi, n);
        let p0 = clmul(h as u64, low);
        let p1 = clmul((h >> 64) as u64, low);
        let (wlo, whi) = (p0 ^ (p1 << 64), p1 >> 64);
        let mut acc = (lo ^ wlo) & self.mask;
        let mut h = shr256(wlo, whi, n) as u64;
        for _ in 1..self.folds {
            let w = clmul(h, low);
            acc ^= w & self.mask;
            h = (w >> n) as u64;
        }
        acc
    }

    /// `P mod f` by Barrett, for a modulus folding would not clear
    /// quickly.
    #[inline(always)]
    fn reduce_barrett(&self, lo: u128, hi: u128, clmul: impl Fn(u64, u64) -> u128 + Copy) -> u128 {
        let n = self.n;
        // H = ⌊P / zⁿ⌋, of degree at most n − 2; Q = ⌊H·μ / zⁿ⌋.
        let h = shr256(lo, hi, n);
        let (ql, qh) = clmul128(h, self.mu, clmul);
        let q = shr256(ql, qh, n);
        // R = (P − Q·f) mod zⁿ = (lo − Q·low) mod zⁿ.
        let (rl, _) = clmul128(q, self.low, clmul);
        (lo ^ rl) & self.mask
    }

    #[inline(always)]
    fn reduce(&self, lo: u128, hi: u128, clmul: impl Fn(u64, u64) -> u128 + Copy) -> u128 {
        if self.folds != 0 {
            self.reduce_fold(lo, hi, clmul)
        } else {
            self.reduce_barrett(lo, hi, clmul)
        }
    }

    /// `a·b` with the given product, inlined into the caller.
    #[inline(always)]
    pub fn mul_by<C: Clmul>(&self, c: C, a: u128, b: u128) -> u128 {
        let (lo, hi) = clmul128(a, b, |x, y| c.clmul(x, y));
        self.reduce(lo, hi, |x, y| c.clmul(x, y))
    }

    /// `a²`: a square has no cross terms, two products, not four.
    #[inline(always)]
    pub fn sqr_by<C: Clmul>(&self, c: C, a: u128) -> u128 {
        let (a0, a1) = (a as u64, (a >> 64) as u64);
        self.reduce(c.clmul(a0, a0), c.clmul(a1, a1), |x, y| c.clmul(x, y))
    }

    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "aes")]
    unsafe fn mul_arm(&self, c: Pmull, a: u128, b: u128) -> u128 {
        self.mul_by(c, a, b)
    }

    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "aes")]
    unsafe fn sqr_arm(&self, c: Pmull, a: u128) -> u128 {
        self.sqr_by(c, a)
    }

    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "aes")]
    unsafe fn inv_arm(&self, c: Pmull, a: u128) -> u128 {
        Generic {
            field: self,
            clmul: c,
        }
        .inv(a)
    }

    /// `a·b`, through the CPU's product where it has one.
    #[inline]
    pub fn mul(&self, a: u128, b: u128) -> u128 {
        #[cfg(target_arch = "x86_64")]
        if let Some(x) = self.xmm {
            // SAFETY: the `Xmm` is the feature check.
            return unsafe { x.mul_u128(a, b) };
        }
        #[cfg(target_arch = "aarch64")]
        if let Some(c) = self.pmull {
            // SAFETY: the `Pmull` token is the feature check.
            return unsafe { self.mul_arm(c, a, b) };
        }
        self.mul_by(Portable, a, b)
    }

    #[inline]
    pub fn sqr(&self, a: u128) -> u128 {
        #[cfg(target_arch = "x86_64")]
        if let Some(x) = self.xmm {
            // SAFETY: the `Xmm` is the feature check.
            return unsafe { x.sqr_u128(a) };
        }
        #[cfg(target_arch = "aarch64")]
        if let Some(c) = self.pmull {
            // SAFETY: the `Pmull` token is the feature check.
            return unsafe { self.sqr_arm(c, a) };
        }
        self.sqr_by(Portable, a)
    }

    /// `a^{2^k}`.
    pub fn sqr_k(&self, mut a: u128, k: u32) -> u128 {
        for _ in 0..k {
            a = self.sqr(a);
        }
        a
    }

    /// `a⁻¹` by Itoh–Tsujii; `0` for `0`.
    pub fn inv(&self, a: u128) -> u128 {
        if a == 0 {
            return 0;
        }
        #[cfg(target_arch = "x86_64")]
        if let Some(x) = self.xmm {
            // SAFETY: the `Xmm` is the feature check.
            return unsafe { x.inv_u128(a) };
        }
        #[cfg(target_arch = "aarch64")]
        if let Some(c) = self.pmull {
            // SAFETY: the `Pmull` token is the feature check.
            return unsafe { self.inv_arm(c, a) };
        }
        self.portable().inv(a)
    }

    fn trace_slow(&self, a: u128) -> bool {
        let mut t = a;
        let mut s = a;
        for _ in 1..self.n {
            s = self.sqr(s);
            t ^= s;
        }
        debug_assert!(t <= 1, "a trace lies in GF(2)");
        t == 1
    }

    /// `Tr(a) = Σ a^{2^i}`, as one AND and a parity.
    #[inline]
    pub fn trace(&self, a: u128) -> bool {
        (a & self.trace_mask).count_ones() & 1 == 1
    }

    /// `√a = a^{2^{n−1}}`.
    pub fn sqrt(&self, a: u128) -> u128 {
        self.sqr_k(a, self.n - 1)
    }

    /// A root `z` of `z² + z = c` with bit 0 clear, or `None` when
    /// `Tr(c) = 1`.  Odd degrees only (the half-trace); construction and
    /// point lifting use it, never the walk.
    pub fn solve_quadratic(&self, c: u128) -> Option<u128> {
        if self.n.is_multiple_of(2) || self.trace(c) {
            return None;
        }
        let mut z = c;
        let mut s = c;
        for _ in 0..(self.n - 1) / 2 {
            s = self.sqr_k(s, 2);
            z ^= s;
        }
        debug_assert_eq!(self.sqr(z) ^ z, c);
        Some(if z & 1 == 1 { z ^ 1 } else { z })
    }

    pub fn from_element(&self, e: &F2mElement) -> u128 {
        let v = e.to_biguint();
        let digits = v.to_u64_digits();
        let lo = digits.first().copied().unwrap_or(0) as u128;
        let hi = digits.get(1).copied().unwrap_or(0) as u128;
        (lo | (hi << 64)) & self.mask
    }

    pub fn to_element(&self, v: u128) -> F2mElement {
        F2mElement::from_biguint(&BigUint::from(v), self.n)
    }

    pub fn mask(&self) -> u128 {
        self.mask
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The least sparse irreducible of each degree, by the general
    /// field's own irreducibility test through its arithmetic.
    fn modulus(n: u32, low: &[u32]) -> IrreduciblePoly {
        IrreduciblePoly {
            degree: n,
            low_terms: low.to_vec(),
        }
    }

    fn rng(seed: &mut u64) -> u64 {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        *seed
    }

    #[test]
    fn wide_field_agrees_with_the_general_field() {
        // The gate's modulus (AGENTS.md §8a), ECC2K-95's, and trinomials
        // and pentanomials at the ends of the range.
        let cases = [
            modulus(83, &[0, 1, 2, 45]),
            modulus(97, &[0, 6]),
            modulus(127, &[0, 1]),
            modulus(64, &[0, 1, 3, 4]),
            modulus(65, &[0, 18]),
            modulus(113, &[0, 9]),
            modulus(31, &[0, 3]),
        ];
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        for irr in cases {
            let f = GfWide::new(&irr).expect("a wide field");
            let n = irr.degree;
            for _ in 0..200 {
                let a = (((rng(&mut seed) as u128) << 64) | rng(&mut seed) as u128) & f.mask();
                let b = (((rng(&mut seed) as u128) << 64) | rng(&mut seed) as u128) & f.mask();
                let (ea, eb) = (f.to_element(a), f.to_element(b));
                assert_eq!(f.mul(a, b), f.from_element(&ea.mul(&eb, &irr)), "n = {n}");
                assert_eq!(f.sqr(a), f.from_element(&ea.square(&irr)), "n = {n}");
                if a != 0 {
                    assert_eq!(f.mul(a, f.inv(a)), 1, "n = {n}");
                }
                assert_eq!(f.sqr(f.sqrt(a)), a, "n = {n}");
                assert_eq!(f.trace(a), f.trace_slow(a), "n = {n}");
                if let Some(z) = f.solve_quadratic(a) {
                    assert_eq!(f.sqr(z) ^ z, a, "n = {n}");
                    assert_eq!(z & 1, 0);
                } else if n % 2 == 1 {
                    assert!(f.trace(a), "n = {n}");
                }
            }
        }
    }

    #[test]
    fn the_kernels_agree() {
        // Every kernel this CPU has, against the portable one, on the
        // moduli the agreement test uses and one dense enough for Barrett.
        let mut seed = 0x1234_5678_9abc_def1u64;
        for (n, low) in [
            (83u32, vec![0u32, 1, 2, 45]),
            (97, vec![0, 6]),
            (127, vec![0, 1]),
            (64, vec![0, 1, 3, 4]),
            (65, vec![0, 18]),
            (31, vec![0, 3]),
            (2, vec![0, 1]),
        ] {
            let f = GfWide::new(&modulus(n, &low)).unwrap();
            let portable = f.portable();
            for _ in 0..500 {
                let a = (((rng(&mut seed) as u128) << 64) | rng(&mut seed) as u128) & f.mask();
                let b = (((rng(&mut seed) as u128) << 64) | rng(&mut seed) as u128) & f.mask();
                let want = (portable.mul(a, b), portable.sqr(a));
                assert_eq!((f.mul(a, b), f.sqr(a)), want, "n = {n}");
                if a != 0 {
                    assert_eq!(f.inv(a), portable.inv(a), "n = {n}");
                }
                #[cfg(target_arch = "x86_64")]
                if let Some(x) = f.xmm() {
                    let got = (
                        x.store(x.mul(x.load(a), x.load(b))),
                        x.store(x.sqr(x.load(a))),
                    );
                    assert_eq!(got, want, "n = {n}");
                    // The shaped kernel for this modulus, as the walk picks it.
                    fn shaped<K: Kernel>(k: K, a: u128, b: u128) -> (u128, u128) {
                        (
                            k.store(k.mul(k.load(a), k.load(b))),
                            k.store(k.sqr(k.load(a))),
                        )
                    }
                    let got = match x.shape() {
                        (true, 1) => shaped(XmmShaped::<2, 1>(x), a, b),
                        (true, 2) => shaped(XmmShaped::<2, 2>(x), a, b),
                        (true, 3) => shaped(XmmShaped::<2, 3>(x), a, b),
                        (false, 1) => shaped(XmmShaped::<1, 1>(x), a, b),
                        (false, 2) => shaped(XmmShaped::<1, 2>(x), a, b),
                        (false, 3) => shaped(XmmShaped::<1, 3>(x), a, b),
                        _ => shaped(XmmShaped::<0, 0>(x), a, b),
                    };
                    assert_eq!(got, want, "n = {n}, shaped");
                }
            }
        }
    }

    #[test]
    fn folding_and_barrett_reduce_alike() {
        // Every modulus the agreement test uses, reduced both ways; the
        // fold counts are the ones the degree bound predicts.
        assert_eq!(fold_count(83, 45), 3);
        assert_eq!(fold_count(127, 1), 1);
        assert_eq!(fold_count(97, 6), 2);
        assert_eq!(fold_count(64, 63), 0);
        let mut seed = 0x0bad_5eed_1234_5678u64;
        for (n, low) in [
            (83u32, vec![0u32, 1, 2, 45]),
            (97, vec![0, 6]),
            (127, vec![0, 1]),
            (65, vec![0, 18]),
            (113, vec![0, 9]),
            (31, vec![0, 3]),
        ] {
            let f = GfWide::new(&modulus(n, &low)).unwrap();
            assert_ne!(f.folds, 0, "n = {n}");
            for _ in 0..500 {
                let a = (((rng(&mut seed) as u128) << 64) | rng(&mut seed) as u128) & f.mask();
                let b = (((rng(&mut seed) as u128) << 64) | rng(&mut seed) as u128) & f.mask();
                let (lo, hi) = clmul128(a, b, clmul64_portable);
                assert_eq!(
                    f.reduce_fold(lo, hi, clmul64_portable),
                    f.reduce_barrett(lo, hi, clmul64_portable),
                    "n = {n}"
                );
            }
        }
    }
}
