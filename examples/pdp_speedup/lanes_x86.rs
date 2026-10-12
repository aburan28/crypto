//! AVX2 (256 lanes) and AVX-512F (512 lanes) vector lanes for the bit-sliced
//! kernel.  The vector types are only ever *run* inside functions compiled
//! with the matching `target_feature` (see `bskernel`), where every method
//! below inlines; `mul` is a separate, feature-compiled function per width so
//! that one copy of the schoolbook product serves every call.

use crate::bitslice::*;
use std::arch::x86_64::*;

#[derive(Clone, Copy)]
pub struct W256(pub __m256i);

impl Lane for W256 {
    const LANES: usize = 256;
    const WORDS: usize = 4;
    #[inline(always)]
    fn from_words(w: &[u64]) -> Self {
        assert!(w.len() >= 4);
        // SAFETY: four words are readable; the load is unaligned.
        unsafe { W256(_mm256_loadu_si256(w.as_ptr() as *const __m256i)) }
    }
    #[inline(always)]
    fn to_words(self, w: &mut [u64]) {
        assert!(w.len() >= 4);
        // SAFETY: four words are writable; the store is unaligned.
        unsafe { _mm256_storeu_si256(w.as_mut_ptr() as *mut __m256i, self.0) }
    }
    #[inline(always)]
    fn zero() -> Self {
        // SAFETY: only reached in code compiled with avx2 (see module comment).
        unsafe { W256(_mm256_setzero_si256()) }
    }
    #[inline(always)]
    fn ones() -> Self {
        // SAFETY: as above.
        unsafe { W256(_mm256_set1_epi64x(-1)) }
    }
    #[inline(always)]
    fn xor(self, o: Self) -> Self {
        // SAFETY: as above.
        unsafe { W256(_mm256_xor_si256(self.0, o.0)) }
    }
    #[inline(always)]
    fn and(self, o: Self) -> Self {
        // SAFETY: as above.
        unsafe { W256(_mm256_and_si256(self.0, o.0)) }
    }
    #[inline(always)]
    fn or(self, o: Self) -> Self {
        // SAFETY: as above.
        unsafe { W256(_mm256_or_si256(self.0, o.0)) }
    }
    #[inline(always)]
    fn not(self) -> Self {
        // SAFETY: as above.
        unsafe { W256(_mm256_xor_si256(self.0, _mm256_set1_epi64x(-1))) }
    }
    #[inline(always)]
    fn xand(acc: Self, a: Self, b: Self) -> Self {
        // SAFETY: as above.
        unsafe { W256(_mm256_xor_si256(acc.0, _mm256_and_si256(a.0, b.0))) }
    }
    #[inline(always)]
    fn is_zero(self) -> bool {
        // SAFETY: as above.
        unsafe { _mm256_testz_si256(self.0, self.0) == 1 }
    }
    #[inline(always)]
    fn mul<const N: usize>(f: &BsField<N>, a: &E<Self>, b: &E<Self>) -> E<Self> {
        // SAFETY: callers run inside an avx2 function.
        unsafe { mul_avx2::<N>(f, a, b) }
    }
}

#[target_feature(enable = "avx2")]
#[inline(never)]
unsafe fn mul_avx2<const N: usize>(f: &BsField<N>, a: &E<W256>, b: &E<W256>) -> E<W256> {
    tiled_mul::<W256, N, 4>(f, a, b)
}

#[derive(Clone, Copy)]
pub struct W512(pub __m512i);

impl Lane for W512 {
    const LANES: usize = 512;
    const WORDS: usize = 8;
    #[inline(always)]
    fn from_words(w: &[u64]) -> Self {
        assert!(w.len() >= 8);
        // SAFETY: eight words are readable; the load is unaligned.
        unsafe { W512(_mm512_loadu_si512(w.as_ptr() as *const _)) }
    }
    #[inline(always)]
    fn to_words(self, w: &mut [u64]) {
        assert!(w.len() >= 8);
        // SAFETY: eight words are writable; the store is unaligned.
        unsafe { _mm512_storeu_si512(w.as_mut_ptr() as *mut _, self.0) }
    }
    #[inline(always)]
    fn zero() -> Self {
        // SAFETY: only reached in code compiled with avx512f.
        unsafe { W512(_mm512_setzero_si512()) }
    }
    #[inline(always)]
    fn ones() -> Self {
        // SAFETY: as above.
        unsafe { W512(_mm512_set1_epi64(-1)) }
    }
    #[inline(always)]
    fn xor(self, o: Self) -> Self {
        // SAFETY: as above.
        unsafe { W512(_mm512_xor_si512(self.0, o.0)) }
    }
    #[inline(always)]
    fn and(self, o: Self) -> Self {
        // SAFETY: as above.
        unsafe { W512(_mm512_and_si512(self.0, o.0)) }
    }
    #[inline(always)]
    fn or(self, o: Self) -> Self {
        // SAFETY: as above.
        unsafe { W512(_mm512_or_si512(self.0, o.0)) }
    }
    #[inline(always)]
    fn not(self) -> Self {
        // SAFETY: as above.  imm 0x55 = NOT of the third operand.
        unsafe { W512(_mm512_ternarylogic_epi64::<0x55>(self.0, self.0, self.0)) }
    }
    /// `acc ^ (a & b)` in one `vpternlogq`: with operands `(acc, a, b)` the truth table
    /// `A ^ (B & C)` is `0x78`.
    #[inline(always)]
    fn xand(acc: Self, a: Self, b: Self) -> Self {
        // SAFETY: as above.
        unsafe { W512(_mm512_ternarylogic_epi64::<0x78>(acc.0, a.0, b.0)) }
    }
    #[inline(always)]
    fn is_zero(self) -> bool {
        // SAFETY: as above.
        unsafe { _mm512_test_epi64_mask(self.0, self.0) == 0 }
    }
    #[inline(always)]
    fn mul<const N: usize>(f: &BsField<N>, a: &E<Self>, b: &E<Self>) -> E<Self> {
        // SAFETY: callers run inside an avx512f function.
        unsafe { mul_avx512::<N>(f, a, b) }
    }
}

#[target_feature(enable = "avx512f")]
#[inline(never)]
unsafe fn mul_avx512<const N: usize>(f: &BsField<N>, a: &E<W512>, b: &E<W512>) -> E<W512> {
    tiled_mul::<W512, N, 8>(f, a, b)
}
