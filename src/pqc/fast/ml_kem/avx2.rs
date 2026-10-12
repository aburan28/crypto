//! AVX2 ring arithmetic for [`super`]: NTT, inverse NTT, the pointwise
//! product and Barrett reduction, sixteen 16-bit coefficients per vector.
//!
//! # Bit-identical, not merely congruent
//!
//! Every function here returns *exactly* the array its scalar counterpart
//! returns, not just the same values mod q. Two facts make that possible:
//!
//! * **Montgomery.** The scalar `montgomery_reduce(a·b)` is
//!   `(a·b − t·q) >> 16` with `t = lo16(a·b)·q⁻¹`. The low halves of `a·b`
//!   and `t·q` agree by construction, so that is exactly
//!   `hi16(a·b) − hi16(t·q)`, which is what `vpmulhw`/`vpmullw` compute.
//! * **Barrett.** The scalar `⌊(v·a + 2²⁵) / 2²⁶⌋` equals
//!   `⌊(⌊v·a / 2¹⁶⌋ + 2⁹) / 2¹⁰⌋` by the nested-floor identity, which is
//!   `vpmulhw` then an add and an arithmetic shift.
//!
//! So the coefficient order is the scalar code's too. The pq-crystals AVX2
//! Kyber permutes coefficients inside its NTT to save shuffles and repacks at
//! the byte boundary; this keeps the standard order and pays three shuffles
//! per vector pair on the three short layers instead, which is what lets the
//! scalar and vector paths be mixed freely and tested for equality.
//!
//! # The short layers
//!
//! Layers with butterfly distance 16 or more pair whole vectors. Distances 8,
//! 4 and 2 pair halves, quarters and dword-pairs of the *same* vector, so two
//! adjacent vectors are first regrouped so that the two operands of each
//! butterfly sit in the same lane of two different registers. Each regrouping
//! (`shuffle8`, `shuffle4`, `shuffle2`) is its own inverse, so applying it
//! again restores the order.

#![allow(clippy::missing_safety_doc)]

use super::{Poly, INVNTT_F, MONT2, N, Q, QINV, ZETAS};

// The parent's twiddle tables are `static`, which a `const fn` cannot read,
// so the tables below are built from the same generators.
const ZT: [i16; 128] = super::zetas_table();
const BZT: [i16; 128] = super::basemul_zetas_table();
use core::arch::x86_64::*;

// ── constant tables ──────────────────────────────────────────────────────────

/// `lo16(z · q⁻¹)`: the second operand Montgomery multiplication by the
/// constant `z` needs.
const fn qinv_of(z: i16) -> i16 {
    (z as i32).wrapping_mul(QINV as i32) as i16
}

const fn with_qinv<const M: usize>(t: [[i16; 16]; M]) -> [[i16; 16]; M] {
    let mut r = [[0i16; 16]; M];
    let mut m = 0;
    while m < M {
        let mut l = 0;
        while l < 16 {
            r[m][l] = qinv_of(t[m][l]);
            l += 1;
        }
        m += 1;
    }
    r
}

/// Distance-8 layer: vector pair `(2m, 2m+1)` regrouped as `[A.lo, B.lo]`,
/// `[A.hi, B.hi]`. Each vector is one block, so the left half takes A's
/// twiddle and the right half B's. `first` is the twiddle index of block 0
/// and `step` is +1 for the forward transform, −1 for the inverse.
const fn tab8(first: i32, step: i32) -> [[i16; 16]; 8] {
    let mut r = [[0i16; 16]; 8];
    let mut m = 0;
    while m < 8 {
        let a = (first + step * (2 * m as i32)) as usize;
        let b = (first + step * (2 * m as i32 + 1)) as usize;
        let mut l = 0;
        while l < 16 {
            r[m][l] = ZT[if l < 8 { a } else { b }];
            l += 1;
        }
        m += 1;
    }
    r
}

/// Distance-4 layer: `[A.q0, B.q0, A.q2, B.q2]` against
/// `[A.q1, B.q1, A.q3, B.q3]`. A holds blocks 4m and 4m+1, B holds 4m+2 and
/// 4m+3, so the quarter-groups take blocks 4m, 4m+2, 4m+1, 4m+3.
const fn tab4(first: i32, step: i32) -> [[i16; 16]; 8] {
    let mut r = [[0i16; 16]; 8];
    let order = [0, 2, 1, 3];
    let mut m = 0;
    while m < 8 {
        let mut l = 0;
        while l < 16 {
            let blk = 4 * m as i32 + order[l / 4];
            r[m][l] = ZT[(first + step * blk) as usize];
            l += 1;
        }
        m += 1;
    }
    r
}

/// Distance-2 layer: `[A.d0, B.d0, A.d2, B.d2 | A.d4, B.d4, A.d6, B.d6]`
/// against the odd dwords. A holds blocks 8m..8m+4, B holds 8m+4..8m+8, so
/// the coefficient pairs take blocks a0, b0, a1, b1, a2, b2, a3, b3.
const fn tab2(first: i32, step: i32) -> [[i16; 16]; 8] {
    let mut r = [[0i16; 16]; 8];
    let order = [0, 4, 1, 5, 2, 6, 3, 7];
    let mut m = 0;
    while m < 8 {
        let mut l = 0;
        while l < 16 {
            let blk = 8 * m as i32 + order[l / 2];
            r[m][l] = ZT[(first + step * blk) as usize];
            l += 1;
        }
        m += 1;
    }
    r
}

/// `BASEMUL_ZETAS[i]` on both coefficients of pair `i`, in vector rows.
const fn gamma_table() -> [[i16; 16]; 16] {
    let mut r = [[0i16; 16]; 16];
    let mut i = 0;
    while i < 128 {
        r[i / 8][2 * (i % 8)] = BZT[i];
        r[i / 8][2 * (i % 8) + 1] = BZT[i];
        i += 1;
    }
    r
}

// Forward twiddles run upward from 16, 32 and 64; inverse ones downward from
// 31, 63 and 127, exactly as the scalar loops consume them.
static FWD8: [[i16; 16]; 8] = tab8(16, 1);
static FWD8Q: [[i16; 16]; 8] = with_qinv(tab8(16, 1));
static FWD4: [[i16; 16]; 8] = tab4(32, 1);
static FWD4Q: [[i16; 16]; 8] = with_qinv(tab4(32, 1));
static FWD2: [[i16; 16]; 8] = tab2(64, 1);
static FWD2Q: [[i16; 16]; 8] = with_qinv(tab2(64, 1));
static INV2: [[i16; 16]; 8] = tab2(127, -1);
static INV2Q: [[i16; 16]; 8] = with_qinv(tab2(127, -1));
static INV4: [[i16; 16]; 8] = tab4(63, -1);
static INV4Q: [[i16; 16]; 8] = with_qinv(tab4(63, -1));
static INV8: [[i16; 16]; 8] = tab8(31, -1);
static INV8Q: [[i16; 16]; 8] = with_qinv(tab8(31, -1));
static GAMMA: [[i16; 16]; 16] = gamma_table();
static GAMMAQ: [[i16; 16]; 16] = with_qinv(gamma_table());

// ── lane arithmetic ──────────────────────────────────────────────────────────

#[inline(always)]
unsafe fn splat(x: i16) -> __m256i {
    _mm256_set1_epi16(x)
}

#[inline(always)]
unsafe fn row(t: &[i16; 16]) -> __m256i {
    _mm256_loadu_si256(t.as_ptr() as *const __m256i)
}

/// `fqmul(x, z)` for a constant `z` with `zq = lo16(z·q⁻¹)` precomputed.
#[inline(always)]
unsafe fn fqmul_c(x: __m256i, z: __m256i, zq: __m256i) -> __m256i {
    let t = _mm256_mullo_epi16(x, zq);
    let hi = _mm256_mulhi_epi16(x, z);
    let t = _mm256_mulhi_epi16(t, splat(Q));
    _mm256_sub_epi16(hi, t)
}

/// `fqmul(x, y)` for two variables.
#[inline(always)]
unsafe fn fqmul_v(x: __m256i, y: __m256i) -> __m256i {
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let t = _mm256_mullo_epi16(lo, splat(QINV));
    let t = _mm256_mulhi_epi16(t, splat(Q));
    _mm256_sub_epi16(hi, t)
}

/// The scalar `barrett_reduce`, lane for lane.
#[inline(always)]
unsafe fn barrett(a: __m256i) -> __m256i {
    const V: i16 = (((1i32 << 26) / Q as i32) + 1) as i16;
    let h = _mm256_mulhi_epi16(a, splat(V));
    let t = _mm256_srai_epi16::<10>(_mm256_add_epi16(h, splat(1 << 9)));
    _mm256_sub_epi16(a, _mm256_mullo_epi16(t, splat(Q)))
}

/// Cooley–Tukey butterfly: `(a, b) ← (a + ζb, a − ζb)`.
#[inline(always)]
unsafe fn ct(a: &mut __m256i, b: &mut __m256i, z: __m256i, zq: __m256i) {
    let t = fqmul_c(*b, z, zq);
    *b = _mm256_sub_epi16(*a, t);
    *a = _mm256_add_epi16(*a, t);
}

/// Gentleman–Sande butterfly: `(a, b) ← (barrett(a + b), ζ(b − a))`.
#[inline(always)]
unsafe fn gs(a: &mut __m256i, b: &mut __m256i, z: __m256i, zq: __m256i) {
    let t = *a;
    *a = barrett(_mm256_add_epi16(t, *b));
    *b = fqmul_c(_mm256_sub_epi16(*b, t), z, zq);
}

#[inline(always)]
unsafe fn shuffle8(a: __m256i, b: __m256i) -> (__m256i, __m256i) {
    (
        _mm256_permute2x128_si256::<0x20>(a, b),
        _mm256_permute2x128_si256::<0x31>(a, b),
    )
}

#[inline(always)]
unsafe fn shuffle4(a: __m256i, b: __m256i) -> (__m256i, __m256i) {
    (_mm256_unpacklo_epi64(a, b), _mm256_unpackhi_epi64(a, b))
}

#[inline(always)]
unsafe fn shuffle2(a: __m256i, b: __m256i) -> (__m256i, __m256i) {
    (
        _mm256_blend_epi32::<0xaa>(a, _mm256_slli_epi64::<32>(b)),
        _mm256_blend_epi32::<0xaa>(_mm256_srli_epi64::<32>(a), b),
    )
}

#[inline(always)]
unsafe fn load(r: &Poly) -> [__m256i; 16] {
    let p = r.0.as_ptr() as *const __m256i;
    let mut v = [_mm256_setzero_si256(); 16];
    for (i, x) in v.iter_mut().enumerate() {
        *x = _mm256_loadu_si256(p.add(i));
    }
    v
}

#[inline(always)]
unsafe fn store(r: &mut Poly, v: &[__m256i; 16]) {
    let p = r.0.as_mut_ptr() as *mut __m256i;
    for (i, x) in v.iter().enumerate() {
        _mm256_storeu_si256(p.add(i), *x);
    }
}

// ── the transforms ───────────────────────────────────────────────────────────

/// Forward NTT, then `poly_reduce`: the scalar `ntt`, bit for bit.
#[target_feature(enable = "avx2")]
pub unsafe fn ntt(r: &mut Poly) {
    let mut v = load(r);
    // Distances 128, 64, 32, 16: whole vectors, twiddles 1..16 in order.
    let mut k = 1;
    let mut len = 8; // in vectors
    while len >= 1 {
        let mut start = 0;
        while start < 16 {
            let z = splat(ZETAS[k]);
            let zq = splat(qinv_of(ZETAS[k]));
            k += 1;
            for j in start..start + len {
                let (lo, hi) = v.split_at_mut(j + len);
                ct(&mut lo[j], &mut hi[0], z, zq);
            }
            start += 2 * len;
        }
        len >>= 1;
    }
    // Distances 8, 4, 2: regroup each vector pair, butterfly, restore.
    for m in 0..8 {
        let (mut a, mut b) = shuffle8(v[2 * m], v[2 * m + 1]);
        ct(&mut a, &mut b, row(&FWD8[m]), row(&FWD8Q[m]));
        let (x, y) = shuffle8(a, b);

        let (mut a, mut b) = shuffle4(x, y);
        ct(&mut a, &mut b, row(&FWD4[m]), row(&FWD4Q[m]));
        let (x, y) = shuffle4(a, b);

        let (mut a, mut b) = shuffle2(x, y);
        ct(&mut a, &mut b, row(&FWD2[m]), row(&FWD2Q[m]));
        let (x, y) = shuffle2(a, b);

        v[2 * m] = barrett(x);
        v[2 * m + 1] = barrett(y);
    }
    store(r, &v);
}

/// Inverse NTT including the final `2^16/128` scaling: the scalar `invntt`,
/// bit for bit.
#[target_feature(enable = "avx2")]
pub unsafe fn invntt(r: &mut Poly) {
    let mut v = load(r);
    // Distances 2, 4, 8 first (Gentleman–Sande runs the layers backwards).
    for m in 0..8 {
        let (mut a, mut b) = shuffle2(v[2 * m], v[2 * m + 1]);
        gs(&mut a, &mut b, row(&INV2[m]), row(&INV2Q[m]));
        let (x, y) = shuffle2(a, b);

        let (mut a, mut b) = shuffle4(x, y);
        gs(&mut a, &mut b, row(&INV4[m]), row(&INV4Q[m]));
        let (x, y) = shuffle4(a, b);

        let (mut a, mut b) = shuffle8(x, y);
        gs(&mut a, &mut b, row(&INV8[m]), row(&INV8Q[m]));
        let (x, y) = shuffle8(a, b);

        v[2 * m] = x;
        v[2 * m + 1] = y;
    }
    // Distances 16, 32, 64, 128: twiddles 15 down to 1.
    let mut k = 15usize;
    let mut len = 1; // in vectors
    while len <= 8 {
        let mut start = 0;
        while start < 16 {
            let z = splat(ZETAS[k]);
            let zq = splat(qinv_of(ZETAS[k]));
            k -= 1;
            for j in start..start + len {
                let (lo, hi) = v.split_at_mut(j + len);
                gs(&mut lo[j], &mut hi[0], z, zq);
            }
            start += 2 * len;
        }
        len <<= 1;
    }
    let f = splat(INVNTT_F);
    let fq = splat(qinv_of(INVNTT_F));
    for x in v.iter_mut() {
        *x = fqmul_c(*x, f, fq);
    }
    store(r, &v);
}

/// `r += a ∘ b`: the scalar `poly_basemul_acc`, bit for bit.
///
/// Coefficients stay interleaved as pairs `(c₀, c₁)`. With `P = a·b` and
/// `S = a·swap(b)` lane-wise, the even output `a₀b₀ + γ·a₁b₁` is
/// `P + swap(γP)` and the odd output `a₀b₁ + a₁b₀` is `S + swap(S)`, each
/// read in its own parity of lanes; one blend takes both.
#[target_feature(enable = "avx2")]
pub unsafe fn basemul_acc(r: &mut Poly, a: &Poly, b: &Poly) {
    let swap = _mm256_setr_epi8(
        2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13, 2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9,
        14, 15, 12, 13,
    );
    let pa = a.0.as_ptr() as *const __m256i;
    let pb = b.0.as_ptr() as *const __m256i;
    let pr = r.0.as_mut_ptr() as *mut __m256i;
    for i in 0..N / 16 {
        let x = _mm256_loadu_si256(pa.add(i));
        let y = _mm256_loadu_si256(pb.add(i));
        let p = fqmul_v(x, y);
        let s = fqmul_v(x, _mm256_shuffle_epi8(y, swap));
        let g = fqmul_c(p, row(&GAMMA[i]), row(&GAMMAQ[i]));
        let even = _mm256_add_epi16(p, _mm256_shuffle_epi8(g, swap));
        let odd = _mm256_add_epi16(s, _mm256_shuffle_epi8(s, swap));
        let t = _mm256_blend_epi16::<0xaa>(even, odd);
        let acc = _mm256_loadu_si256(pr.add(i));
        _mm256_storeu_si256(pr.add(i), _mm256_add_epi16(acc, t));
    }
}

/// The scalar `poly_reduce`.
#[target_feature(enable = "avx2")]
pub unsafe fn reduce(r: &mut Poly) {
    let p = r.0.as_mut_ptr() as *mut __m256i;
    for i in 0..N / 16 {
        _mm256_storeu_si256(p.add(i), barrett(_mm256_loadu_si256(p.add(i))));
    }
}

/// The scalar `poly_tomont`.
#[target_feature(enable = "avx2")]
pub unsafe fn tomont(r: &mut Poly) {
    let f = splat(MONT2);
    let fq = splat(qinv_of(MONT2));
    let p = r.0.as_mut_ptr() as *mut __m256i;
    for i in 0..N / 16 {
        _mm256_storeu_si256(p.add(i), fqmul_c(_mm256_loadu_si256(p.add(i)), f, fq));
    }
}

/// `Compress_d(to_canonical(barrett_reduce(f)))` for every coefficient: what
/// the scalar `poly_compress_encode` feeds its packer, value for value.
///
/// The scalar code divides by q, which the compiler turns into a 64-bit
/// multiply per coefficient. Here the division is done four lanes per
/// `vpmuludq` with the reciprocal `M = ⌈2³⁵/q⌉ = 10321340`. Writing
/// `M·q = 2³⁵ + e` with `e = 2492`, `⌊n·M/2³⁵⌋ = ⌊n/q⌋` whenever `n·e < 2³⁵`,
/// that is for `n` below about 13.7 million; the largest `n` here is
/// `(q−1)·2¹¹ + 1664 < 2²³`. The tests also check every input exhaustively.
#[target_feature(enable = "avx2")]
pub unsafe fn compress(f: &Poly, d: u32, out: &mut [u16; N]) {
    const M: i32 = 10_321_340;
    debug_assert!((1..=11).contains(&d));
    let m = _mm256_set1_epi32(M);
    let half = _mm256_set1_epi32(Q as i32 / 2);
    let dv = _mm256_set1_epi32(d as i32);
    let mask = _mm256_set1_epi32((1 << d) - 1);
    let pf = f.0.as_ptr() as *const __m256i;
    let po = out.as_mut_ptr() as *mut __m256i;

    for i in 0..N / 16 {
        let a = barrett(_mm256_loadu_si256(pf.add(i)));
        // to_canonical: add q where negative.
        let a = _mm256_add_epi16(a, _mm256_and_si256(_mm256_srai_epi16::<15>(a), splat(Q)));
        let lo = compress8(
            _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a)),
            dv,
            half,
            m,
            mask,
        );
        let hi = compress8(
            _mm256_cvtepu16_epi32(_mm256_extracti128_si256::<1>(a)),
            dv,
            half,
            m,
            mask,
        );
        // packus interleaves 128-bit halves; the permute restores order.
        let r = _mm256_permute4x64_epi64::<0xd8>(_mm256_packus_epi32(lo, hi));
        _mm256_storeu_si256(po.add(i), r);
    }
}

/// `⌊(x·2^d + q/2) / q⌋ mod 2^d` for eight 32-bit lanes; see [`compress`].
#[inline(always)]
unsafe fn compress8(x: __m256i, dv: __m256i, half: __m256i, m: __m256i, mask: __m256i) -> __m256i {
    let n = _mm256_add_epi32(_mm256_sllv_epi32(x, dv), half);
    let even = _mm256_srli_epi64::<35>(_mm256_mul_epu32(n, m));
    let odd = _mm256_srli_epi64::<35>(_mm256_mul_epu32(_mm256_srli_epi64::<32>(n), m));
    let quot = _mm256_blend_epi32::<0xaa>(even, _mm256_slli_epi64::<32>(odd));
    _mm256_and_si256(quot, mask)
}

/// For each 8-bit accept mask, the `pshufb` control that moves the accepted
/// 16-bit lanes of a 128-bit half to its front, in order.
const fn compact_table() -> [[u8; 16]; 256] {
    let mut t = [[0x80u8; 16]; 256];
    let mut m = 0;
    while m < 256 {
        let mut k = 0;
        let mut j = 0;
        while j < 8 {
            if (m >> j) & 1 == 1 {
                t[m][2 * k] = 2 * j as u8;
                t[m][2 * k + 1] = 2 * j as u8 + 1;
                k += 1;
            }
            j += 1;
        }
        m += 1;
    }
    t
}

static COMPACT: [[u8; 16]; 256] = compact_table();

/// The vector half of `RejBuf::parse`: consume 24-byte groups of `buf` (16
/// candidates each) while at least 16 output slots remain below N, appending
/// the accepted values to `out` from position `n`, in order. Returns the new
/// count and the number of bytes consumed; the scalar loop finishes the rest.
///
/// Each group is loaded as two overlapping 16-byte halves, bytes 0..16 and
/// 8..24, so nothing past the 24 bytes is read. A byte shuffle puts the two
/// bytes each 12-bit candidate straddles into one 16-bit lane; even lanes are
/// masked to 12 bits, odd lanes shifted down by 4. Accepted lanes are then
/// compacted per half with a table-driven `pshufb` and stored unaligned, so
/// `out` needs 8 spare slots past N.
#[target_feature(enable = "avx2")]
pub unsafe fn rej_uniform(out: &mut [i16; N + 16], mut n: usize, buf: &[u8]) -> (usize, usize) {
    let spread = _mm256_setr_epi8(
        0, 1, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 10, 10, 11, 4, 5, 5, 6, 7, 8, 8, 9, 10, 11, 11, 12,
        13, 14, 14, 15,
    );
    let q = splat(Q);
    let low12 = splat(0x0fff);
    let mut pos = 0;
    while n + 16 <= N && pos + 24 <= buf.len() {
        let p = buf.as_ptr().add(pos);
        let lo = _mm_loadu_si128(p as *const __m128i);
        let hi = _mm_loadu_si128(p.add(8) as *const __m128i);
        let v = _mm256_shuffle_epi8(_mm256_set_m128i(hi, lo), spread);
        let v = _mm256_blend_epi16::<0xaa>(_mm256_and_si256(v, low12), _mm256_srli_epi16::<4>(v));
        let good = _mm256_cmpgt_epi16(q, v);
        // packs puts each lane's accept bit in one byte, per 128-bit half.
        let bits = _mm256_movemask_epi8(_mm256_packs_epi16(good, _mm256_setzero_si256())) as u32;
        let (m0, m1) = ((bits & 0xff) as usize, ((bits >> 16) & 0xff) as usize);
        let o = out.as_mut_ptr();
        let v0 = _mm256_castsi256_si128(v);
        let v1 = _mm256_extracti128_si256::<1>(v);
        let c0 = _mm_loadu_si128(COMPACT[m0].as_ptr() as *const __m128i);
        _mm_storeu_si128(o.add(n) as *mut __m128i, _mm_shuffle_epi8(v0, c0));
        n += m0.count_ones() as usize;
        let c1 = _mm_loadu_si128(COMPACT[m1].as_ptr() as *const __m128i);
        _mm_storeu_si128(o.add(n) as *mut __m128i, _mm_shuffle_epi8(v1, c1));
        n += m1.count_ones() as usize;
        pos += 24;
    }
    (n, pos)
}
