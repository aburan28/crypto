use serde::{Deserialize, Serialize};
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

// Prime modulus. Must satisfy P^2 < 2^64 for the naive u64 mul to be safe.
// 1_000_003 ≈ 10^6 (15× bigger than the original toy 65521), so the curve has
// ~10^6 points — large enough that group_order takes ~1 s instead of ms, and
// brute-force sqrt becomes infeasible (hence Tonelli-Shanks below).
pub const P: u64 = 1_000_003;

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct Fp(pub u64);

impl Fp {
    pub fn new(x: i64) -> Self {
        Fp(x.rem_euclid(P as i64) as u64)
    }
    pub fn zero() -> Self { Fp(0) }
    pub fn one() -> Self { Fp(1) }
    pub fn is_zero(&self) -> bool { self.0 == 0 }

    // Fermat inverse: a^(p-2) mod p
    pub fn inv(&self) -> Option<Self> {
        if self.is_zero() { return None; }
        Some(pow(*self, P - 2))
    }

    // Square root via Tonelli-Shanks. Returns None if `self` is not a QR.
    // Used by curve::point_at_x once P is bigger than ~10^4 (brute-force y
    // scan becomes prohibitively slow above that).
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() { return Some(Fp::zero()); }
        // Fast path for p ≡ 3 (mod 4): sqrt(n) = n^((p+1)/4) if it's a QR.
        if P % 4 == 3 {
            let r = pow(*self, (P + 1) / 4);
            return if r * r == *self { Some(r) } else { None };
        }
        // General Tonelli-Shanks (works for any odd prime p).
        // 1. Factor p-1 = q * 2^s with q odd.
        let mut q = P - 1;
        let mut s: u32 = 0;
        while q & 1 == 0 { q >>= 1; s += 1; }
        // 2. Find a quadratic non-residue z.
        let mut z = Fp::new(2);
        while is_qr(z) {
            z = z + Fp::one();
        }
        // 3. Initialize.
        let mut m = s;
        let mut c = pow(z, q);
        let mut t = pow(*self, q);
        let mut r = pow(*self, (q + 1) / 2);
        // 4. Loop.
        loop {
            if t == Fp::one() { return Some(r); }
            // Smallest i, 0 < i < m, with t^(2^i) == 1.
            let mut t_pow = t;
            let mut i: u32 = 0;
            while t_pow != Fp::one() {
                t_pow = t_pow * t_pow;
                i += 1;
                if i == m { return None; } // self was not a QR after all.
            }
            // b = c^(2^(m-i-1))
            let mut b = c;
            for _ in 0..(m - i - 1) { b = b * b; }
            m = i;
            c = b * b;
            t = t * c;
            r = r * b;
        }
    }
}

// Modular exponentiation.
pub fn pow(base: Fp, exp: u64) -> Fp {
    let mut result = Fp::one();
    let mut b = base;
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 { result = result * b; }
        b = b * b;
        e >>= 1;
    }
    result
}

// Euler-criterion quadratic-residue test. True iff `a` is 0 or a nonzero square.
pub fn is_qr(a: Fp) -> bool {
    if a.is_zero() { return true; }
    pow(a, (P - 1) / 2) == Fp::one()
}

impl Add for Fp { type Output = Self; fn add(self, o: Self) -> Self { Fp((self.0 + o.0) % P) } }
impl Sub for Fp { type Output = Self; fn sub(self, o: Self) -> Self { Fp((self.0 + P - o.0) % P) } }
impl Mul for Fp { type Output = Self; fn mul(self, o: Self) -> Self { Fp(self.0 * o.0 % P) } }
impl Neg for Fp { type Output = Self; fn neg(self) -> Self { Fp((P - self.0) % P) } }

impl fmt::Display for Fp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { write!(f, "{}", self.0) }
}
