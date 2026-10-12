//! Small constant-time helpers on top of `subtle`.

pub use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

/// All-ones mask if `c` is true, zero otherwise.
#[inline(always)]
pub fn mask(c: Choice) -> u64 {
    0u64.wrapping_sub(c.unwrap_u8() as u64)
}

/// Choice from a word that is zero or one.  The value passes through
/// `black_box` so the optimiser cannot turn later selections into branches.
#[inline(always)]
pub fn choice_from_bit(bit: u64) -> Choice {
    Choice::from(core::hint::black_box(bit as u8 & 1))
}

/// `a == 0` as a Choice, branch-free.
#[inline(always)]
pub fn is_zero_u64(a: u64) -> Choice {
    // (a | -a) has its top bit set iff a != 0.
    choice_from_bit(((a | a.wrapping_neg()) >> 63) ^ 1)
}

/// Constant-time select of u64 words: `b` if `c`, else `a`.
#[inline(always)]
pub fn sel(a: u64, b: u64, c: Choice) -> u64 {
    u64::conditional_select(&a, &b, c)
}

/// Constant-time comparison of byte strings of equal length.
pub fn bytes_eq(a: &[u8], b: &[u8]) -> Choice {
    if a.len() != b.len() {
        return Choice::from(0);
    }
    a.ct_eq(b)
}
