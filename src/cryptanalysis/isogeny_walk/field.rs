//! `GF(p)` for any odd prime `p < 2^256`, in four-limb Montgomery form.
//!
//! Variable time: the walker handles only public curve data.  One
//! implementation serves P-256, P-224 and the small primes the tests
//! cross-check against, so the code a test exercises is the code a walk
//! runs.

use num_bigint::BigUint;
use num_traits::{One, Zero};

/// An element in Montgomery form (`x·R mod p`, `R = 2^256`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Fe(pub [u64; 4]);

/// The field `GF(p)`.
#[derive(Clone, Debug)]
pub struct Field {
    p: [u64; 4],
    /// `-p^{-1} mod 2^64`.
    n0: u64,
    /// `R^2 mod p`.
    r2: Fe,
    one: Fe,
    p_big: BigUint,
    /// `(p - 1) / 2`, for the Legendre symbol.
    half: BigUint,
    /// Tonelli–Shanks: `p - 1 = 2^s · q`, `q` odd, and a non-residue's
    /// `q`-th power.
    ts_s: u32,
    ts_q: BigUint,
    ts_z: Fe,
    /// The least quadratic non-residue, as an integer.
    nonresidue: u64,
}

fn to_limbs(v: &BigUint) -> [u64; 4] {
    let mut out = [0u64; 4];
    for (i, d) in v.to_u64_digits().iter().enumerate().take(4) {
        out[i] = *d;
    }
    out
}

fn from_limbs(v: &[u64; 4]) -> BigUint {
    let mut bytes = Vec::with_capacity(32);
    for limb in v {
        bytes.extend_from_slice(&limb.to_le_bytes());
    }
    BigUint::from_bytes_le(&bytes)
}

fn geq(a: &[u64; 4], b: &[u64; 4]) -> bool {
    for i in (0..4).rev() {
        if a[i] != b[i] {
            return a[i] > b[i];
        }
    }
    true
}

fn sub_limbs(a: &[u64; 4], b: &[u64; 4]) -> ([u64; 4], bool) {
    let mut out = [0u64; 4];
    let mut borrow = false;
    for i in 0..4 {
        let (d, b1) = a[i].overflowing_sub(b[i]);
        let (d, b2) = d.overflowing_sub(borrow as u64);
        out[i] = d;
        borrow = b1 || b2;
    }
    (out, borrow)
}

fn add_limbs(a: &[u64; 4], b: &[u64; 4]) -> ([u64; 4], bool) {
    let mut out = [0u64; 4];
    let mut carry = false;
    for i in 0..4 {
        let (s, c1) = a[i].overflowing_add(b[i]);
        let (s, c2) = s.overflowing_add(carry as u64);
        out[i] = s;
        carry = c1 || c2;
    }
    (out, carry)
}

impl Field {
    /// The field of the odd prime `p < 2^256`.  Primality is the caller's
    /// claim; [`Field::new`] checks it with Miller–Rabin and refuses a
    /// composite.
    pub fn new(p: &BigUint) -> Option<Self> {
        if p.bits() > 256 || p < &BigUint::from(5u8) || !p.bit(0) || !is_probable_prime(p) {
            return None;
        }
        let limbs = to_limbs(p);
        // Newton iteration for p^{-1} mod 2^64.
        let mut inv: u64 = 1;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(limbs[0].wrapping_mul(inv)));
        }
        let r = BigUint::one() << 256;
        let mut f = Field {
            p: limbs,
            n0: inv.wrapping_neg(),
            r2: Fe(to_limbs(&((&r * &r) % p))),
            one: Fe(to_limbs(&(&r % p))),
            p_big: p.clone(),
            half: (p - 1u8) >> 1,
            ts_s: 0,
            ts_q: BigUint::zero(),
            ts_z: Fe([0; 4]),
            nonresidue: 0,
        };
        let pm1 = p - 1u8;
        let s = pm1.trailing_zeros().unwrap_or(0) as u32;
        f.ts_s = s;
        f.ts_q = &pm1 >> s;
        let mut z = 2u64;
        while f.legendre(&f.from_u64(z)) != -1 {
            z += 1;
        }
        f.nonresidue = z;
        f.ts_z = f.pow(&f.from_u64(z), &f.ts_q.clone());
        Some(f)
    }

    pub fn modulus(&self) -> &BigUint {
        &self.p_big
    }

    /// The least quadratic non-residue, `≥ 2`.
    pub fn least_nonresidue(&self) -> u64 {
        self.nonresidue
    }

    pub fn zero(&self) -> Fe {
        Fe([0; 4])
    }

    pub fn one(&self) -> Fe {
        self.one
    }

    pub fn is_zero(&self, a: &Fe) -> bool {
        a.0 == [0; 4]
    }

    pub fn mul(&self, a: &Fe, b: &Fe) -> Fe {
        // CIOS Montgomery multiplication; the result is < 2p before the
        // final subtraction for any odd p < 2^256.
        let p = &self.p;
        let mut t = [0u64; 6];
        for i in 0..4 {
            let mut c: u128 = 0;
            for j in 0..4 {
                let s = t[j] as u128 + (a.0[j] as u128) * (b.0[i] as u128) + c;
                t[j] = s as u64;
                c = s >> 64;
            }
            let s = t[4] as u128 + c;
            t[4] = s as u64;
            t[5] = (s >> 64) as u64;
            let m = t[0].wrapping_mul(self.n0);
            let s = t[0] as u128 + (m as u128) * (p[0] as u128);
            let mut c = s >> 64;
            for j in 1..4 {
                let s = t[j] as u128 + (m as u128) * (p[j] as u128) + c;
                t[j - 1] = s as u64;
                c = s >> 64;
            }
            let s = t[4] as u128 + c;
            t[3] = s as u64;
            t[4] = t[5] + (s >> 64) as u64;
            t[5] = 0;
        }
        let r = [t[0], t[1], t[2], t[3]];
        if t[4] != 0 || geq(&r, p) {
            Fe(sub_limbs(&r, p).0)
        } else {
            Fe(r)
        }
    }

    pub fn sqr(&self, a: &Fe) -> Fe {
        self.mul(a, a)
    }

    pub fn add(&self, a: &Fe, b: &Fe) -> Fe {
        let (s, carry) = add_limbs(&a.0, &b.0);
        if carry || geq(&s, &self.p) {
            Fe(sub_limbs(&s, &self.p).0)
        } else {
            Fe(s)
        }
    }

    pub fn sub(&self, a: &Fe, b: &Fe) -> Fe {
        let (d, borrow) = sub_limbs(&a.0, &b.0);
        if borrow {
            Fe(add_limbs(&d, &self.p).0)
        } else {
            Fe(d)
        }
    }

    pub fn neg(&self, a: &Fe) -> Fe {
        self.sub(&self.zero(), a)
    }

    pub fn from_u64(&self, v: u64) -> Fe {
        self.from_big(&BigUint::from(v))
    }

    pub fn from_i64(&self, v: i64) -> Fe {
        let m = self.from_u64(v.unsigned_abs());
        if v < 0 {
            self.neg(&m)
        } else {
            m
        }
    }

    pub fn from_u128(&self, v: u128) -> Fe {
        self.from_big(&BigUint::from(v))
    }

    pub fn from_big(&self, v: &BigUint) -> Fe {
        let v = v % &self.p_big;
        self.mul(&Fe(to_limbs(&v)), &self.r2)
    }

    pub fn to_big(&self, a: &Fe) -> BigUint {
        from_limbs(&self.mul(a, &Fe([1, 0, 0, 0])).0)
    }

    pub fn pow(&self, a: &Fe, e: &BigUint) -> Fe {
        let mut acc = self.one;
        for i in (0..e.bits()).rev() {
            acc = self.sqr(&acc);
            if e.bit(i) {
                acc = self.mul(&acc, a);
            }
        }
        acc
    }

    /// `a^{-1}`; `None` for zero.
    pub fn inv(&self, a: &Fe) -> Option<Fe> {
        if self.is_zero(a) {
            return None;
        }
        Some(self.pow(a, &(&self.p_big - 2u8)))
    }

    pub fn div(&self, a: &Fe, b: &Fe) -> Option<Fe> {
        Some(self.mul(a, &self.inv(b)?))
    }

    /// `1`, `-1`, or `0` for zero.
    pub fn legendre(&self, a: &Fe) -> i32 {
        if self.is_zero(a) {
            return 0;
        }
        if self.pow(a, &self.half) == self.one {
            1
        } else {
            -1
        }
    }

    /// A square root by Tonelli–Shanks, or `None` for a non-residue.
    pub fn sqrt(&self, a: &Fe) -> Option<Fe> {
        match self.legendre(a) {
            0 => return Some(self.zero()),
            -1 => return None,
            _ => {}
        }
        let mut m = self.ts_s;
        let mut c = self.ts_z;
        let mut t = self.pow(a, &self.ts_q);
        let mut r = self.pow(a, &((&self.ts_q + 1u8) >> 1));
        while t != self.one {
            let mut i = 0;
            let mut t2 = t;
            while t2 != self.one {
                t2 = self.sqr(&t2);
                i += 1;
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = self.sqr(&b);
            }
            m = i;
            c = self.sqr(&b);
            t = self.mul(&t, &c);
            r = self.mul(&r, &b);
        }
        Some(r)
    }

    /// The canonical integer of `a` as lower-case hex, no prefix.
    pub fn hex(&self, a: &Fe) -> String {
        format!("{:x}", self.to_big(a))
    }
}

/// Miller–Rabin with the first twenty prime bases: deterministic far
/// beyond 64 bits and an error below `4^-20` for the moduli used here.
pub fn is_probable_prime(n: &BigUint) -> bool {
    let small = [
        2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    ];
    if n < &BigUint::from(2u8) {
        return false;
    }
    for &q in &small {
        let q = BigUint::from(q);
        if n == &q {
            return true;
        }
        if (n % &q).is_zero() {
            return false;
        }
    }
    let nm1 = n - 1u8;
    let s = nm1.trailing_zeros().unwrap_or(0);
    let d = &nm1 >> s;
    'base: for &a in &small {
        let mut x = BigUint::from(a).modpow(&d, n);
        if x.is_one() || x == nm1 {
            continue;
        }
        for _ in 1..s {
            x = x.modpow(&BigUint::from(2u8), n);
            if x == nm1 {
                continue 'base;
            }
        }
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_matches_biguint_on_p224_and_small_primes() {
        let p224 = BigUint::parse_bytes(
            b"ffffffffffffffffffffffffffffffff000000000000000000000001",
            16,
        )
        .unwrap();
        for p in [p224, BigUint::from(1_000_003u32), BigUint::from(101u32)] {
            let f = Field::new(&p).unwrap();
            let a = BigUint::from(123_456_789_123u64) % &p;
            let b = BigUint::from(987_654_321_987u64) % &p;
            let (fa, fb) = (f.from_big(&a), f.from_big(&b));
            assert_eq!(f.to_big(&f.mul(&fa, &fb)), (&a * &b) % &p);
            assert_eq!(f.to_big(&f.add(&fa, &fb)), (&a + &b) % &p);
            assert_eq!(f.to_big(&f.sub(&fa, &fb)), (&a + &p - &b) % &p);
            let inv = f.inv(&fa).unwrap();
            assert_eq!(f.mul(&inv, &fa), f.one());
            let sq = f.sqr(&fb);
            let r = f.sqrt(&sq).unwrap();
            assert_eq!(f.sqr(&r), sq);
            let nr = f.from_u64(f.least_nonresidue());
            assert_eq!(f.legendre(&nr), -1);
            assert!(f.sqrt(&nr).is_none());
        }
    }
}
