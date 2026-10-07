//! Exact integer arithmetic for Koblitz-curve Weil polynomials.
//!
//! `E_a : y² + xy = x³ + a·x² + 1` over `F_2` has Frobenius trace
//! `t = (−1)^{a+1}`; `s_k` is the trace of the `2^k`-power Frobenius.  The Weil
//! restriction from `F_{2^n}` has characteristic polynomial
//! `T^{2n} − s_n·T^n + 2^n`, which splits as `(T² − tT + 2) · P_A(T)` with
//! `P_A` the trace-zero part.

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// Frobenius trace of `E_a` over `F_2`: −1 for `a = 0`, +1 for `a = 1`.
pub fn koblitz_trace(a: u8) -> i64 {
    if a == 0 {
        -1
    } else {
        1
    }
}

/// `s_k = α^k + ᾱ^k` for the Frobenius of a curve over `F_2` with trace `t`.
pub fn frobenius_trace(t: i64, k: u32) -> BigInt {
    let (mut s0, mut s1) = (BigInt::from(2), BigInt::from(t));
    if k == 0 {
        return s0;
    }
    for _ in 1..k {
        let next = &s1 * t - &s0 * 2;
        s0 = s1;
        s1 = next;
    }
    s1
}

/// `#E(F_{2^k}) = 2^k + 1 − s_k`.
pub fn curve_order(t: i64, k: u32) -> BigInt {
    (BigInt::one() << k) + 1 - frobenius_trace(t, k)
}

/// `T² − tT + 2`, little-endian.
pub fn curve_poly(t: i64) -> Vec<BigInt> {
    vec![BigInt::from(2), BigInt::from(-t), BigInt::one()]
}

/// `T^{2n} − s_n T^n + 2^n`, little-endian.
pub fn weil_restriction_poly(t: i64, n: u32) -> Vec<BigInt> {
    let mut p = vec![BigInt::zero(); 2 * n as usize + 1];
    p[0] = BigInt::one() << n;
    p[n as usize] = -frobenius_trace(t, n);
    p[2 * n as usize] = BigInt::one();
    p
}

/// Division by a monic polynomial, returning quotient and remainder.
pub fn divide_monic(num: &[BigInt], den: &[BigInt]) -> (Vec<BigInt>, Vec<BigInt>) {
    assert!(
        den.last().is_some_and(|c| c.is_one()),
        "denominator must be monic"
    );
    let mut rem = num.to_vec();
    if num.len() < den.len() {
        return (vec![BigInt::zero()], rem);
    }
    let mut quo = vec![BigInt::zero(); num.len() - den.len() + 1];
    for i in (0..quo.len()).rev() {
        let c = rem[i + den.len() - 1].clone();
        if !c.is_zero() {
            for (j, d) in den.iter().enumerate() {
                rem[i + j] -= &c * d;
            }
        }
        quo[i] = c;
    }
    while rem.len() > 1 && rem.last().is_some_and(|c| c.is_zero()) {
        rem.pop();
    }
    (quo, rem)
}

pub fn eval(p: &[BigInt], x: &BigInt) -> BigInt {
    p.iter().rev().fold(BigInt::zero(), |acc, c| acc * x + c)
}

/// The trace-zero part `P_A = P_W / P_E` of the Weil restriction from `F_{2^n}`.
pub fn trace_zero_poly(t: i64, n: u32) -> Vec<BigInt> {
    let (quo, rem) = divide_monic(&weil_restriction_poly(t, n), &curve_poly(t));
    assert!(rem.iter().all(|c| c.is_zero()), "P_E must divide P_W");
    quo
}

pub fn to_i128(p: &[BigInt]) -> Vec<i128> {
    p.iter()
        .map(|c| c.to_i128().expect("coefficient fits i128"))
        .collect()
}

/// `log₂ x` for a positive big integer, accurate to ~1e-15 relative.
pub fn log2_big(x: &BigUint) -> f64 {
    assert!(!x.is_zero(), "log of zero");
    let bits = x.bits();
    if bits <= 64 {
        return (x.to_u64().expect("fits") as f64).log2();
    }
    let shift = bits - 64;
    let top = (x >> shift).to_u64().expect("fits");
    (top as f64).log2() + shift as f64
}

pub fn log2_bigint(x: &BigInt) -> f64 {
    log2_big(&x.abs().to_biguint().expect("abs is non-negative"))
}

/// Deterministic Miller–Rabin over the first twelve primes (ample below 2^128,
/// and a strong probable-prime test above it).
pub fn is_probable_prime(n: &BigUint) -> bool {
    let two = BigUint::from(2u32);
    if n < &two {
        return false;
    }
    let bases = [2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for &p in &bases {
        let p = BigUint::from(p);
        if n == &p {
            return true;
        }
        if (n % &p).is_zero() {
            return false;
        }
    }
    let n1 = n - 1u32;
    let s = n1.trailing_zeros().expect("n > 1");
    let d = &n1 >> s;
    'outer: for &a in &bases {
        let mut x = BigUint::from(a).modpow(&d, n);
        if x.is_one() || x == n1 {
            continue;
        }
        for _ in 1..s {
            x = x.modpow(&two, n);
            if x == n1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

/// Tonelli–Shanks square root modulo an odd prime.
pub fn sqrt_mod(a: &BigUint, p: &BigUint) -> Option<BigUint> {
    let a = a % p;
    if a.is_zero() {
        return Some(a);
    }
    let one = BigUint::one();
    let p1 = p - 1u32;
    let half = &p1 >> 1;
    if a.modpow(&half, p) != one {
        return None;
    }
    let s = p1.trailing_zeros().expect("p > 1");
    let q = &p1 >> s;
    let mut z = BigUint::from(2u32);
    while z.modpow(&half, p) != p1 {
        z += 1u32;
    }
    let mut m = s;
    let mut c = z.modpow(&q, p);
    let mut t = a.modpow(&q, p);
    let mut r = a.modpow(&((&q + 1u32) >> 1), p);
    while !t.is_one() {
        let mut i = 0u64;
        let mut t2 = t.clone();
        while !t2.is_one() {
            t2 = &t2 * &t2 % p;
            i += 1;
        }
        let b = c.modpow(&(BigUint::one() << (m - i - 1)), p);
        m = i;
        c = &b * &b % p;
        t = &t * &c % p;
        r = &r * &b % p;
    }
    Some(r)
}

/// Multiplicative order of `a` modulo `m` (`gcd(a, m) = 1`).
pub fn mult_order(a: u64, m: u64) -> u64 {
    assert_eq!(a.gcd(&m), 1, "a must be a unit");
    let (mut x, mut k) = (a % m, 1u64);
    while x != 1 % m {
        x = x * a % m;
        k += 1;
    }
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_koblitz_orders() {
        // E_0 over F_2, F_4, F_8: 4, 8, 4 points; E_1: 2, 8, 14.
        let e0: Vec<i64> = (1..=3)
            .map(|k| curve_order(-1, k).to_i64().unwrap())
            .collect();
        let e1: Vec<i64> = (1..=3)
            .map(|k| curve_order(1, k).to_i64().unwrap())
            .collect();
        assert_eq!(e0, vec![4, 8, 4]);
        assert_eq!(e1, vec![2, 8, 14]);
    }

    #[test]
    fn trace_zero_part_has_r_points() {
        let pa = trace_zero_poly(-1, 131);
        assert_eq!(pa.len(), 261);
        let r = curve_order(-1, 131) / 4;
        assert_eq!(eval(&pa, &BigInt::one()), r);
        assert!(is_probable_prime(&r.to_biguint().unwrap()));
    }

    #[test]
    fn tonelli_shanks_round_trip() {
        let p = BigUint::from(1_000_000_007u64);
        let a = BigUint::from(5u32);
        let r = sqrt_mod(&(&a * &a % &p), &p).unwrap();
        assert!(r == a || r == &p - &a);
    }
}
