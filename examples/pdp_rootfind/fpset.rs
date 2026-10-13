//! Fingerprint buckets for the 2-sum tables: one 64-byte line per probe.
//!
//! Each key hashes to a bucket of 32 sixteen-bit fingerprints and to a
//! fingerprint; a probe reads the one line and compares all 32 lanes
//! (SWAR).  No false negatives: a bucket that overflowed is flagged and
//! answers "maybe" to everything.  The caller confirms a "maybe" in the
//! sorted key array, so membership stays exact.  Cost: 4 bytes a key at
//! the chosen load (about 16 fingerprints a bucket), false-positive rate
//! about `16 / 2^16` a probe.

#[repr(C, align(64))]
#[derive(Clone, Copy)]
struct Bucket([u64; 8]);

pub struct FpSet {
    buckets: Vec<Bucket>,
    overflow: Vec<u64>,
    shift: u32,
}

const C1: u64 = 0x9E37_79B9_7F4A_7C15;
const C2: u64 = 0xC2B2_AE3D_27D4_EB4F;
const LANE_LO: u64 = 0x0001_0001_0001_0001;
const LANE_HI: u64 = 0x8000_8000_8000_8000;

#[inline(always)]
fn split(x: u64, shift: u32) -> (usize, u64) {
    let b = (x.wrapping_mul(C1) >> shift) as usize;
    let f = x.wrapping_mul(C2) >> 48;
    (b, if f == 0 { 1 } else { f })
}

impl FpSet {
    /// A set for `entries` keys: `entries / 16` buckets, rounded up to a
    /// power of two (at least two).
    pub fn new(entries: usize) -> Self {
        let nb = (entries / 16).max(2).next_power_of_two();
        FpSet {
            buckets: vec![Bucket([0; 8]); nb],
            overflow: vec![0u64; nb.div_ceil(64)],
            shift: 64 - nb.trailing_zeros(),
        }
    }

    pub fn insert(&mut self, x: u64) {
        let (b, f) = split(x, self.shift);
        let words = &mut self.buckets[b].0;
        for w in words.iter_mut() {
            for lane in 0..4 {
                let v = (*w >> (16 * lane)) & 0xFFFF;
                if v == f {
                    return;
                }
                if v == 0 {
                    *w |= f << (16 * lane);
                    return;
                }
            }
        }
        self.overflow[b >> 6] |= 1u64 << (b & 63);
    }

    /// Prefetch the line a later probe of `x` will read.
    #[inline(always)]
    pub fn prefetch(&self, x: u64) {
        #[cfg(target_arch = "x86_64")]
        {
            let b = (x.wrapping_mul(C1) >> self.shift) as usize;
            // SAFETY: `b < buckets.len()` by construction of `shift`; a
            // prefetch has no effect on program state.
            unsafe {
                use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
                _mm_prefetch(self.buckets.as_ptr().add(b) as *const i8, _MM_HINT_T0);
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = x;
        }
    }

    /// No false negatives.
    #[inline(always)]
    pub fn maybe(&self, x: u64) -> bool {
        let (b, f) = split(x, self.shift);
        let words = &self.buckets[b].0;
        let pat = f * LANE_LO;
        let mut acc = 0u64;
        for &w in words {
            let v = w ^ pat;
            acc |= v.wrapping_sub(LANE_LO) & !v & LANE_HI;
        }
        acc != 0 || (self.overflow[b >> 6] >> (b & 63)) & 1 == 1
    }

    pub fn bytes(&self) -> usize {
        self.buckets.len() * 64 + self.overflow.len() * 8
    }
}

/// Exact membership of every candidate in a sorted key array through the
/// fingerprint set, lines prefetched `D` probes ahead.
pub fn contains_batch(fp: &FpSet, keys: &[u64], cands: &[u64], hits: &mut Vec<bool>) {
    const D: usize = 16;
    hits.clear();
    hits.reserve(cands.len());
    for &c in cands.iter().take(D) {
        fp.prefetch(c);
    }
    for i in 0..cands.len() {
        if i + D < cands.len() {
            fp.prefetch(cands[i + D]);
        }
        let c = cands[i];
        hits.push(fp.maybe(c) && keys.binary_search(&c).is_ok());
    }
}
