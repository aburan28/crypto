//! Minimal unsigned big integers (little-endian u64 limbs) for exponents, field sizes and moduli.
//! Only what the field code needs; not performance critical.
use std::cmp::Ordering;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Big(pub Vec<u64>);

impl Big {
    pub fn zero() -> Big {
        Big(vec![])
    }
    pub fn from_u64(v: u64) -> Big {
        Big::norm(vec![v])
    }
    pub fn from_u128(v: u128) -> Big {
        Big::norm(vec![v as u64, (v >> 64) as u64])
    }
    pub fn from_limbs(l: &[u64]) -> Big {
        Big::norm(l.to_vec())
    }
    fn norm(mut v: Vec<u64>) -> Big {
        while v.last() == Some(&0) {
            v.pop();
        }
        Big(v)
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    pub fn bits(&self) -> usize {
        match self.0.last() {
            None => 0,
            Some(&t) => 64 * (self.0.len() - 1) + (64 - t.leading_zeros() as usize),
        }
    }
    pub fn bit(&self, i: usize) -> bool {
        let (w, b) = (i / 64, i % 64);
        w < self.0.len() && (self.0[w] >> b) & 1 == 1
    }
    pub fn limbs(&self, n: usize) -> Vec<u64> {
        let mut v = self.0.clone();
        v.resize(n, 0);
        v
    }
    pub fn to_u128(&self) -> Option<u128> {
        if self.0.len() > 2 {
            return None;
        }
        let l = self.limbs(2);
        Some(l[0] as u128 | (l[1] as u128) << 64)
    }
    pub fn add(&self, o: &Big) -> Big {
        let n = self.0.len().max(o.0.len());
        let (a, b) = (self.limbs(n), o.limbs(n));
        let mut r = Vec::with_capacity(n + 1);
        let mut c = 0u128;
        for i in 0..n {
            let s = a[i] as u128 + b[i] as u128 + c;
            r.push(s as u64);
            c = s >> 64;
        }
        r.push(c as u64);
        Big::norm(r)
    }
    /// self - o, requires self >= o.
    pub fn sub(&self, o: &Big) -> Big {
        assert!(self.cmp(o) != Ordering::Less, "Big::sub underflow");
        let n = self.0.len();
        let b = o.limbs(n);
        let mut r = Vec::with_capacity(n);
        let mut borrow = 0u64;
        for i in 0..n {
            let (d1, o1) = self.0[i].overflowing_sub(b[i]);
            let (d2, o2) = d1.overflowing_sub(borrow);
            r.push(d2);
            borrow = (o1 || o2) as u64;
        }
        Big::norm(r)
    }
    pub fn mul(&self, o: &Big) -> Big {
        if self.is_zero() || o.is_zero() {
            return Big::zero();
        }
        let mut r = vec![0u64; self.0.len() + o.0.len()];
        for (i, &a) in self.0.iter().enumerate() {
            let mut c = 0u128;
            for (j, &b) in o.0.iter().enumerate() {
                let s = r[i + j] as u128 + a as u128 * b as u128 + c;
                r[i + j] = s as u64;
                c = s >> 64;
            }
            r[i + o.0.len()] = c as u64;
        }
        Big::norm(r)
    }
    pub fn mul_small(&self, m: u64) -> Big {
        self.mul(&Big::from_u64(m))
    }
    pub fn add_small(&self, m: u64) -> Big {
        self.add(&Big::from_u64(m))
    }
    pub fn sub_small(&self, m: u64) -> Big {
        self.sub(&Big::from_u64(m))
    }
    pub fn shr(&self, k: usize) -> Big {
        let (w, b) = (k / 64, k % 64);
        if w >= self.0.len() {
            return Big::zero();
        }
        let src = &self.0[w..];
        let mut r = vec![0u64; src.len()];
        for i in 0..src.len() {
            r[i] = src[i] >> b;
            if b > 0 && i + 1 < src.len() {
                r[i] |= src[i + 1] << (64 - b);
            }
        }
        Big::norm(r)
    }
    /// (quotient, remainder) by a small divisor.
    pub fn divrem_small(&self, d: u64) -> (Big, u64) {
        let mut q = vec![0u64; self.0.len()];
        let mut r = 0u128;
        for i in (0..self.0.len()).rev() {
            let cur = (r << 64) | self.0[i] as u128;
            q[i] = (cur / d as u128) as u64;
            r = cur % d as u128;
        }
        (Big::norm(q), r as u64)
    }
    pub fn rem_small(&self, d: u64) -> u64 {
        self.divrem_small(d).1
    }
    pub fn from_dec(s: &str) -> Big {
        let mut r = Big::zero();
        for ch in s.chars().filter(|c| !c.is_whitespace() && *c != '_') {
            let d = ch.to_digit(10).expect("decimal digit") as u64;
            r = r.mul_small(10).add_small(d);
        }
        r
    }
    pub fn to_dec(&self) -> String {
        if self.is_zero() {
            return "0".into();
        }
        let mut digits = vec![];
        let mut cur = self.clone();
        while !cur.is_zero() {
            let (q, r) = cur.divrem_small(10_000_000_000_000_000_000);
            digits.push(r);
            cur = q;
        }
        let mut s = digits.last().unwrap().to_string();
        for d in digits.iter().rev().skip(1) {
            s.push_str(&format!("{:019}", d));
        }
        s
    }
}

impl PartialOrd for Big {
    fn partial_cmp(&self, o: &Big) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Big {
    fn cmp(&self, o: &Big) -> Ordering {
        if self.0.len() != o.0.len() {
            return self.0.len().cmp(&o.0.len());
        }
        for i in (0..self.0.len()).rev() {
            if self.0[i] != o.0[i] {
                return self.0[i].cmp(&o.0[i]);
            }
        }
        Ordering::Equal
    }
}

/// Miller-Rabin for big integers (bases 2..=37 plus `extra` random bases); `n` odd > 3.
pub fn is_probable_prime(n: &Big) -> bool {
    if let Some(v) = n.to_u128() {
        if v < (1u128 << 63) {
            return crate::field::is_prime(v as u64);
        }
    }
    for &q in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.rem_small(q) == 0 {
            return false;
        }
    }
    // n - 1 = d 2^s
    let nm1 = n.sub_small(1);
    let mut s = 0;
    while !nm1.bit(s) {
        s += 1;
    }
    let d = nm1.shr(s);
    let modmul = |a: &Big, b: &Big| -> Big { big_mod(&a.mul(b), n) };
    'outer: for &base in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = Big::from_u64(1);
        let b = Big::from_u64(base);
        for i in (0..d.bits()).rev() {
            x = modmul(&x, &x);
            if d.bit(i) {
                x = modmul(&x, &b);
            }
        }
        if x == Big::from_u64(1) || x == nm1 {
            continue;
        }
        for _ in 0..s - 1 {
            x = modmul(&x, &x);
            if x == nm1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

/// a mod n by shift-and-subtract (slow; for setup only).
pub fn big_mod(a: &Big, n: &Big) -> Big {
    if a < n {
        return a.clone();
    }
    let mut r = Big::zero();
    for i in (0..a.bits()).rev() {
        r = r.add(&r);
        if a.bit(i) {
            r = r.add_small(1);
        }
        if r >= *n {
            r = r.sub(n);
        }
    }
    r
}
