//! Signed big integers on top of `bigint::Big` (sign + magnitude), with the operations the
//! quaternion code needs: ring operations, floor division, gcd / extended gcd, integer square
//! roots, modular powers, inverses and square roots. Values compare and hash by value.
use crate::bigint::{is_probable_prime, Big};
use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Int {
    neg: bool,
    mag: Big,
}

impl fmt::Debug for Int {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self)
    }
}
impl fmt::Display for Int {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.neg {
            write!(f, "-")?;
        }
        write!(f, "{}", self.mag.to_dec())
    }
}

impl Int {
    pub fn zero() -> Int {
        Int {
            neg: false,
            mag: Big::zero(),
        }
    }
    pub fn one() -> Int {
        Int::from(1i64)
    }
    fn make(neg: bool, mag: Big) -> Int {
        let neg = neg && !mag.is_zero();
        Int { neg, mag }
    }
    pub fn from_big(b: &Big) -> Int {
        Int::make(false, b.clone())
    }
    pub fn mag(&self) -> &Big {
        &self.mag
    }
    pub fn is_zero(&self) -> bool {
        self.mag.is_zero()
    }
    pub fn is_neg(&self) -> bool {
        self.neg
    }
    pub fn signum(&self) -> i32 {
        if self.is_zero() {
            0
        } else if self.neg {
            -1
        } else {
            1
        }
    }
    pub fn abs(&self) -> Int {
        Int::make(false, self.mag.clone())
    }
    pub fn is_even(&self) -> bool {
        !self.mag.bit(0)
    }
    pub fn bits(&self) -> usize {
        self.mag.bits()
    }
    pub fn to_i128(&self) -> Option<i128> {
        let m = self.mag.to_u128()?;
        if m > i128::MAX as u128 {
            return None;
        }
        Some(if self.neg { -(m as i128) } else { m as i128 })
    }
    pub fn to_f64(&self) -> f64 {
        let mut v = 0f64;
        for &l in self.mag.0.iter().rev() {
            v = v * 18446744073709551616.0 + l as f64;
        }
        if self.neg {
            -v
        } else {
            v
        }
    }
    /// Truncated division (quotient rounds towards zero).
    pub fn divrem_trunc(&self, d: &Int) -> (Int, Int) {
        let (q, r) = self.mag.divrem(&d.mag);
        (Int::make(self.neg != d.neg, q), Int::make(self.neg, r))
    }
    /// Floor division and the matching remainder (sign of the divisor).
    pub fn div_floor(&self, d: &Int) -> Int {
        let (q, r) = self.divrem_trunc(d);
        if !r.is_zero() && (r.neg != d.neg) {
            q - Int::one()
        } else {
            q
        }
    }
    /// Remainder in [0, |m|).
    pub fn modulo(&self, m: &Int) -> Int {
        let (_, r) = self.divrem_trunc(m);
        if r.neg {
            r + m.abs()
        } else {
            r
        }
    }
    /// Nearest integer to self/d (ties away from zero is fine for our use).
    pub fn div_round(&self, d: &Int) -> Int {
        let two = Int::from(2i64);
        (&(self * &two) + d).div_floor(&(d * &two))
    }
    pub fn mod_u64(&self, m: u64) -> u64 {
        let r = self.mag.rem_small(m);
        if self.neg && r != 0 {
            m - r
        } else {
            r
        }
    }
    pub fn gcd(&self, o: &Int) -> Int {
        let (mut a, mut b) = (self.abs(), o.abs());
        while !b.is_zero() {
            let r = a.modulo(&b);
            a = b;
            b = r;
        }
        a
    }
    /// (g, x, y) with a x + b y = g = gcd(a, b) >= 0.
    pub fn xgcd(a: &Int, b: &Int) -> (Int, Int, Int) {
        let (mut r0, mut r1) = (a.clone(), b.clone());
        let (mut s0, mut s1) = (Int::one(), Int::zero());
        let (mut t0, mut t1) = (Int::zero(), Int::one());
        while !r1.is_zero() {
            let q = r0.div_floor(&r1);
            let r2 = &r0 - &(&q * &r1);
            r0 = std::mem::replace(&mut r1, r2);
            let s2 = &s0 - &(&q * &s1);
            s0 = std::mem::replace(&mut s1, s2);
            let t2 = &t0 - &(&q * &t1);
            t0 = std::mem::replace(&mut t1, t2);
        }
        if r0.neg {
            (-r0, -s0, -t0)
        } else {
            (r0, s0, t0)
        }
    }
    pub fn inv_mod(&self, m: &Int) -> Option<Int> {
        let (g, x, _) = Int::xgcd(&self.modulo(m), m);
        if g == Int::one() {
            Some(x.modulo(m))
        } else {
            None
        }
    }
    pub fn pow(&self, mut e: u32) -> Int {
        let mut r = Int::one();
        let mut b = self.clone();
        while e > 0 {
            if e & 1 == 1 {
                r = &r * &b;
            }
            b = &b * &b;
            e >>= 1;
        }
        r
    }
    pub fn pow_mod(&self, e: &Int, m: &Int) -> Int {
        let mut r = Int::one().modulo(m);
        let b = self.modulo(m);
        for i in (0..e.bits()).rev() {
            r = (&r * &r).modulo(m);
            if e.mag.bit(i) {
                r = (&r * &b).modulo(m);
            }
        }
        r
    }
    /// floor(sqrt(self)) for self >= 0.
    pub fn isqrt(&self) -> Int {
        assert!(!self.neg, "isqrt of a negative number");
        Int::from_big(&self.mag.isqrt())
    }
    pub fn is_square(&self) -> bool {
        !self.neg && {
            let r = self.isqrt();
            &r * &r == *self
        }
    }
    pub fn is_probable_prime(&self) -> bool {
        if self.neg || self.mag.bits() < 2 {
            return false;
        }
        if let Some(v) = self.mag.to_u128() {
            if v < 4 {
                return v == 2 || v == 3;
            }
            if v % 2 == 0 {
                return false;
            }
        }
        is_probable_prime(&self.mag)
    }
    /// A square root of a modulo an odd prime m (Tonelli-Shanks), None for non-residues.
    pub fn sqrt_mod_prime(a: &Int, m: &Int) -> Option<Int> {
        let a = a.modulo(m);
        if a.is_zero() {
            return Some(a);
        }
        let one = Int::one();
        let m1 = m - &one;
        let half = m1.div_floor(&Int::from(2i64));
        if a.pow_mod(&half, m) != one {
            return None;
        }
        let mut s = 0usize;
        let mut q = m1.clone();
        while q.is_even() {
            q = q.div_floor(&Int::from(2i64));
            s += 1;
        }
        let mut z = Int::from(2i64);
        while z.pow_mod(&half, m) == one {
            z = &z + &one;
        }
        let mut c = z.pow_mod(&q, m);
        let mut t = a.pow_mod(&q, m);
        let mut r = a.pow_mod(&(&q + &one).div_floor(&Int::from(2i64)), m);
        let mut mm = s;
        while t != one {
            let (mut i, mut tt) = (0usize, t.clone());
            while tt != one {
                tt = (&tt * &tt).modulo(m);
                i += 1;
            }
            let mut b = c.clone();
            for _ in 0..mm - i - 1 {
                b = (&b * &b).modulo(m);
            }
            mm = i;
            c = (&b * &b).modulo(m);
            t = (&t * &c).modulo(m);
            r = (&r * &b).modulo(m);
        }
        Some(r)
    }
}

impl From<i64> for Int {
    fn from(v: i64) -> Int {
        Int::make(v < 0, Big::from_u64(v.unsigned_abs()))
    }
}
impl From<i128> for Int {
    fn from(v: i128) -> Int {
        Int::make(v < 0, Big::from_u128(v.unsigned_abs()))
    }
}
impl From<u64> for Int {
    fn from(v: u64) -> Int {
        Int::make(false, Big::from_u64(v))
    }
}

impl PartialOrd for Int {
    fn partial_cmp(&self, o: &Int) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Int {
    fn cmp(&self, o: &Int) -> Ordering {
        match (self.neg, o.neg) {
            (false, true) => Ordering::Greater,
            (true, false) => Ordering::Less,
            (false, false) => self.mag.cmp(&o.mag),
            (true, true) => o.mag.cmp(&self.mag),
        }
    }
}

fn add_signed(a: &Int, b: &Int, negate_b: bool) -> Int {
    let bneg = b.neg != negate_b && !b.is_zero();
    if a.neg == bneg {
        Int::make(a.neg, a.mag.add(&b.mag))
    } else if a.mag >= b.mag {
        Int::make(a.neg, a.mag.sub(&b.mag))
    } else {
        Int::make(bneg, b.mag.sub(&a.mag))
    }
}

impl Add<&Int> for &Int {
    type Output = Int;
    fn add(self, o: &Int) -> Int {
        add_signed(self, o, false)
    }
}
impl Sub<&Int> for &Int {
    type Output = Int;
    fn sub(self, o: &Int) -> Int {
        add_signed(self, o, true)
    }
}
impl Mul<&Int> for &Int {
    type Output = Int;
    fn mul(self, o: &Int) -> Int {
        Int::make(self.neg != o.neg, self.mag.mul(&o.mag))
    }
}
impl Neg for &Int {
    type Output = Int;
    fn neg(self) -> Int {
        Int::make(!self.neg, self.mag.clone())
    }
}
impl Neg for Int {
    type Output = Int;
    fn neg(self) -> Int {
        Int::make(!self.neg, self.mag)
    }
}
macro_rules! owned_ops {
    ($tr:ident, $m:ident) => {
        impl $tr<Int> for Int {
            type Output = Int;
            fn $m(self, o: Int) -> Int {
                (&self).$m(&o)
            }
        }
        impl $tr<&Int> for Int {
            type Output = Int;
            fn $m(self, o: &Int) -> Int {
                (&self).$m(o)
            }
        }
        impl $tr<Int> for &Int {
            type Output = Int;
            fn $m(self, o: Int) -> Int {
                self.$m(&o)
            }
        }
    };
}
owned_ops!(Add, add);
owned_ops!(Sub, sub);
owned_ops!(Mul, mul);
