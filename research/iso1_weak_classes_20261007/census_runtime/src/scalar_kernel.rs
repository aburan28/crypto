//! Scalar helpers frozen from 62becf9572fe74cbe8b3d8cebee3bf8a240708a1; mul_mod uses exact u128 reduction.
#[inline]
pub(crate) fn mm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
#[inline]
pub(crate) fn am(a: u64, b: u64, p: u64) -> u64 {
    let s = a + b;
    if s >= p {
        s - p
    } else {
        s
    }
}
#[inline]
pub(crate) fn sm(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + p - b
    }
}
fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    mm(a, b, p)
}
pub fn pow_mod(mut base: u64, mut e: u64, p: u64) -> u64 {
    let mut acc = 1u64 % p;
    base %= p;
    while e > 0 {
        if e & 1 == 1 {
            acc = mul_mod(acc, base, p);
        }
        base = mul_mod(base, base, p);
        e >>= 1;
    }
    acc
}
pub fn inv_mod(a: u64, p: u64) -> u64 {
    let (mut old_r, mut r) = (a % p, p);
    // |s_{i−1}| and |s_i|, with s_i = (−1)^i·|s_i| (s_1 = 0).
    let (mut old_s, mut s) = (1u64, 0u64);
    let mut odd = false;
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_s, s) = (s, old_s + q * s);
        odd = !odd;
    }
    debug_assert_eq!(old_r, 1, "inv_mod of a non-unit");
    if odd && old_s != 0 {
        p - old_s
    } else {
        old_s
    }
}
