//! Variable-time number theory for key generation and tests only
//! (num-bigint is not constant time; nothing here touches decapsulation).

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use rand_core::RngCore;

pub fn random_below<R: RngCore>(rng: &mut R, n: &BigUint) -> BigUint {
    let bytes = (n.bits() as usize).div_ceil(8) + 16;
    let mut buf = vec![0u8; bytes];
    rng.fill_bytes(&mut buf);
    BigUint::from_bytes_le(&buf) % n
}

const SMALL_PRIMES: [u32; 54] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
    113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239,
    241, 251,
];

/// Miller–Rabin with deterministic small bases plus `extra` random bases.
pub fn is_probable_prime<R: RngCore>(n: &BigUint, rng: &mut R, extra: usize) -> bool {
    let two = BigUint::from(2u32);
    if *n < two {
        return false;
    }
    for &q in SMALL_PRIMES.iter() {
        let q = BigUint::from(q);
        if *n == q {
            return true;
        }
        if (n % &q).is_zero() {
            return false;
        }
    }
    let nm1 = n - 1u32;
    let s = nm1.trailing_zeros().unwrap_or(0);
    let d = &nm1 >> s;
    let mut bases: Vec<BigUint> = [2u32, 3, 5, 7, 11, 13, 17].iter().map(|&b| BigUint::from(b)).collect();
    for _ in 0..extra {
        bases.push(random_below(rng, &(n - 3u32)) + 2u32);
    }
    'outer: for a in bases {
        let mut x = a.modpow(&d, n);
        if x.is_one() || x == nm1 {
            continue;
        }
        for _ in 1..s {
            x = x.modpow(&two, n);
            if x == nm1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

/// sqrt(a) mod prime q (Tonelli–Shanks); None if a is a non-residue.
pub fn sqrt_mod_prime<R: RngCore>(a: &BigUint, q: &BigUint, rng: &mut R) -> Option<BigUint> {
    let a = a % q;
    if a.is_zero() {
        return Some(BigUint::zero());
    }
    let one = BigUint::one();
    let qm1 = q - 1u32;
    let legendre = |x: &BigUint| x.modpow(&(&qm1 >> 1), q);
    if legendre(&a) != one {
        return None;
    }
    if (q % 4u32) == BigUint::from(3u32) {
        return Some(a.modpow(&((q + 1u32) >> 2), q));
    }
    let s = qm1.trailing_zeros().unwrap();
    let t = &qm1 >> s;
    // z must be a quadratic non-residue (Legendre symbol -1).  Testing only
    // "!= 1" accepted z = 0, after which the loop below never terminated
    // (keygen hung for some RNG streams on small q = 1 mod 8, e.g. q = 17).
    let mut z = BigUint::from(2u32);
    while legendre(&z) != qm1 {
        z = random_below(rng, q);
    }
    let mut m = s;
    let mut c = z.modpow(&t, q);
    let mut x = a.modpow(&((&t + 1u32) >> 1), q);
    let mut b = a.modpow(&t, q);
    while !b.is_one() {
        let mut i = 0u64;
        let mut bb = b.clone();
        while !bb.is_one() {
            bb = (&bb * &bb) % q;
            i += 1;
            if i >= m {
                return None; // q not prime: no element of 2-power order here
            }
        }
        let e = BigUint::one() << (m - i - 1);
        let f = c.modpow(&e, q);
        x = (x * &f) % q;
        c = (&f * &f) % q;
        b = (b * &c) % q;
        m = i;
    }
    Some(x)
}

/// Solve x^2 + D y^2 = m for prime m (Cornacchia), D > 0.
pub fn cornacchia<R: RngCore>(d: u64, m: &BigUint, rng: &mut R) -> Option<(BigUint, BigUint)> {
    let dm = BigUint::from(d) % m;
    let r = sqrt_mod_prime(&((m - &dm) % m), m, rng)?;
    let mut r0 = m.clone();
    let mut r1 = if r.clone() * 2u32 > *m { m - &r } else { r };
    let bound = m.sqrt();
    while r1 > bound {
        let t = &r0 % &r1;
        r0 = r1;
        r1 = t;
    }
    let rest = m - &r1 * &r1;
    if (&rest % d) != BigUint::zero() {
        return None;
    }
    let y2 = rest / d;
    let y = y2.sqrt();
    if &y * &y != y2 {
        return None;
    }
    Some((r1, y))
}

/// Write m = x^2 + y^2 with m = (product of small primes) * (prime q = 1 mod 4).
/// Returns None if m has a prime factor 3 mod 4 to an odd power among the
/// small primes, or the cofactor is not a prime 1 mod 4.
pub fn two_squares<R: RngCore>(m: &BigUint, rng: &mut R, trial_bound: u32) -> Option<(BigUint, BigUint)> {
    let mut rest = m.clone();
    // Gaussian integer accumulator (a + b i)
    let mut acc = (BigInt::one(), BigInt::zero());
    let gmul = |x: &(BigInt, BigInt), y: &(BigInt, BigInt)| (&x.0 * &y.0 - &x.1 * &y.1, &x.0 * &y.1 + &x.1 * &y.0);
    let mut q = 2u32;
    while q <= trial_bound {
        let qb = BigUint::from(q);
        let mut e = 0;
        while (&rest % &qb).is_zero() {
            rest /= &qb;
            e += 1;
        }
        if e > 0 {
            if q == 2 {
                for _ in 0..e {
                    acc = gmul(&acc, &(BigInt::one(), BigInt::one()));
                }
            } else if q % 4 == 3 {
                if e % 2 == 1 {
                    return None;
                }
                let s = BigInt::from(q).pow(e / 2);
                acc = (&acc.0 * &s, &acc.1 * &s);
            } else {
                let (x, y) = cornacchia(1, &qb, rng)?;
                let g = (BigInt::from_biguint(Sign::Plus, x), BigInt::from_biguint(Sign::Plus, y));
                for _ in 0..e {
                    acc = gmul(&acc, &g);
                }
            }
        }
        q += if q == 2 { 1 } else { 2 };
    }
    if !rest.is_one() {
        if (&rest % 4u32) != BigUint::one() || !is_probable_prime(&rest, rng, 4) {
            return None;
        }
        let (x, y) = cornacchia(1, &rest, rng)?;
        acc = gmul(&acc, &(BigInt::from_biguint(Sign::Plus, x), BigInt::from_biguint(Sign::Plus, y)));
    }
    Some((acc.0.abs().to_biguint().unwrap(), acc.1.abs().to_biguint().unwrap()))
}

/// s with s^2 = a mod 2^k for a = 1 mod 8.
pub fn sqrt_mod_2power(a: &BigUint, k: u64) -> BigUint {
    let m = BigUint::one() << k;
    let a = a % &m;
    assert_eq!(&a % 8u32, BigUint::one());
    let mut s = BigUint::one();
    for j in 3..k {
        let mj1 = BigUint::one() << (j + 1);
        let t = (&s * &s) % &mj1;
        if t != (&a % &mj1) {
            s += BigUint::one() << (j - 1);
        }
    }
    s
}

pub fn v2(x: &BigInt) -> u64 {
    if x.is_zero() {
        return 999;
    }
    x.abs().to_biguint().unwrap().trailing_zeros().unwrap()
}

pub fn modpos(x: &BigInt, m: &BigUint) -> BigUint {
    let mm = BigInt::from_biguint(Sign::Plus, m.clone());
    x.mod_floor(&mm).to_biguint().unwrap()
}
