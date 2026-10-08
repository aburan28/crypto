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
    pub fn shl(&self, k: usize) -> Big {
        if self.is_zero() {
            return Big::zero();
        }
        let (w, b) = (k / 64, k % 64);
        let mut r = vec![0u64; self.0.len() + w + 1];
        for (i, &x) in self.0.iter().enumerate() {
            r[i + w] |= x << b;
            if b > 0 {
                r[i + w + 1] |= x >> (64 - b);
            }
        }
        Big::norm(r)
    }
    /// (quotient, remainder) by Knuth's Algorithm D (64-bit limbs, normalised divisor).
    pub fn divrem(&self, d: &Big) -> (Big, Big) {
        assert!(!d.is_zero(), "division by zero");
        if self < d {
            return (Big::zero(), self.clone());
        }
        if d.0.len() == 1 {
            let (q, r) = self.divrem_small(d.0[0]);
            return (q, Big::from_u64(r));
        }
        let n = d.0.len();
        let m = self.0.len() - n;
        let s = d.0[n - 1].leading_zeros();
        let shift = |v: &[u64], len: usize| -> Vec<u64> {
            let mut r = vec![0u64; len];
            for i in 0..v.len() {
                r[i] |= v[i] << s;
                if s > 0 && i + 1 < len {
                    r[i + 1] |= v[i] >> (64 - s);
                }
            }
            r
        };
        let vn = shift(&d.0, n);
        let mut un = shift(&self.0, m + n + 1);
        let mut q = vec![0u64; m + 1];
        let b: u128 = 1 << 64;
        for j in (0..=m).rev() {
            let num = ((un[j + n] as u128) << 64) | un[j + n - 1] as u128;
            let mut qhat = num / vn[n - 1] as u128;
            let mut rhat = num % vn[n - 1] as u128;
            while qhat >= b || qhat * vn[n - 2] as u128 > ((rhat << 64) | un[j + n - 2] as u128) {
                qhat -= 1;
                rhat += vn[n - 1] as u128;
                if rhat >= b {
                    break;
                }
            }
            // un[j..j+n+1] -= qhat * vn
            let mut k: i128 = 0;
            for i in 0..n {
                let p = qhat * vn[i] as u128;
                let t = un[i + j] as i128 - k - (p as u64) as i128;
                un[i + j] = t as u64;
                k = (p >> 64) as i128 - (t >> 64);
            }
            let t = un[j + n] as i128 - k;
            un[j + n] = t as u64;
            if t < 0 {
                qhat -= 1;
                let mut c = 0u128;
                for i in 0..n {
                    let s2 = un[i + j] as u128 + vn[i] as u128 + c;
                    un[i + j] = s2 as u64;
                    c = s2 >> 64;
                }
                un[j + n] = un[j + n].wrapping_add(c as u64);
            }
            q[j] = qhat as u64;
        }
        let mut r = vec![0u64; n];
        for i in 0..n {
            r[i] = un[i] >> s;
            if s > 0 {
                r[i] |= un[i + 1] << (64 - s);
            }
        }
        (Big::norm(q), Big::norm(r))
    }
    /// floor(sqrt(self)) by Newton's iteration.
    pub fn isqrt(&self) -> Big {
        if self.is_zero() {
            return Big::zero();
        }
        let mut x = Big::from_u64(1).shl(self.bits().div_ceil(2));
        loop {
            let y = x.add(&self.divrem(&x).0).shr(1);
            if y >= x {
                return x;
            }
            x = y;
        }
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

/// Montgomery arithmetic with a runtime limb count (for primality tests on candidates of any size).
struct MontRt {
    n: Vec<u64>,
    ninv: u64,
}
impl MontRt {
    fn new(n: &Big) -> Self {
        let mut inv = 1u64;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(n.0[0].wrapping_mul(inv)));
        }
        MontRt { n: n.0.clone(), ninv: inv.wrapping_neg() }
    }
    fn geq(a: &[u64], b: &[u64]) -> bool {
        for i in (0..b.len()).rev() {
            if a[i] != b[i] {
                return a[i] > b[i];
            }
        }
        true
    }
    fn sub_in(a: &mut [u64], b: &[u64]) {
        let mut br = false;
        for i in 0..b.len() {
            let (d1, o1) = a[i].overflowing_sub(b[i]);
            let (d2, o2) = d1.overflowing_sub(br as u64);
            a[i] = d2;
            br = o1 || o2;
        }
    }
    /// 2a mod n for a < n.
    fn dbl(&self, a: &mut [u64]) {
        let mut c = 0u64;
        for x in a.iter_mut() {
            let nc = *x >> 63;
            *x = (*x << 1) | c;
            c = nc;
        }
        if c == 1 || Self::geq(a, &self.n) {
            Self::sub_in(a, &self.n);
        }
    }
    fn mul(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        let k = self.n.len();
        let mut t = vec![0u64; k + 2];
        for i in 0..k {
            let mut c: u128 = 0;
            for j in 0..k {
                let s = t[j] as u128 + a[j] as u128 * b[i] as u128 + c;
                t[j] = s as u64;
                c = s >> 64;
            }
            let s = t[k] as u128 + c;
            t[k] = s as u64;
            t[k + 1] = (s >> 64) as u64;
            let m = t[0].wrapping_mul(self.ninv);
            let mut c = (t[0] as u128 + m as u128 * self.n[0] as u128) >> 64;
            for j in 1..k {
                let s = t[j] as u128 + m as u128 * self.n[j] as u128 + c;
                t[j - 1] = s as u64;
                c = s >> 64;
            }
            let s = t[k] as u128 + c;
            t[k - 1] = s as u64;
            t[k] = t[k + 1] + (s >> 64) as u64;
        }
        let top = t[k];
        t.truncate(k);
        if top != 0 || Self::geq(&t, &self.n) {
            Self::sub_in(&mut t, &self.n);
        }
        t
    }
}

/// Small odd primes for trial division before Miller-Rabin.
fn small_odd_primes() -> &'static [u64] {
    use std::sync::OnceLock;
    static P: OnceLock<Vec<u64>> = OnceLock::new();
    P.get_or_init(|| (3u64..2000).step_by(2).filter(|&q| crate::field::is_prime(q)).collect())
}

/// Miller-Rabin for big integers (bases 2..=37) after trial division by the primes below 2000;
/// Montgomery arithmetic with the candidate's own limb count. `n` odd > 3.
pub fn is_probable_prime(n: &Big) -> bool {
    if let Some(v) = n.to_u128() {
        if v < (1u128 << 63) {
            return crate::field::is_prime(v as u64);
        }
    }
    if !n.bit(0) {
        return false;
    }
    for &q in small_odd_primes() {
        if n.rem_small(q) == 0 {
            return false;
        }
    }
    let mt = MontRt::new(n);
    let k = n.0.len();
    // R mod n (Montgomery one) by doubling 1 64k times
    let mut one = vec![0u64; k];
    one[0] = 1;
    for _ in 0..64 * k {
        mt.dbl(&mut one);
    }
    // n - 1 = d 2^s; in Montgomery form, -1 is n - R
    let nm1 = n.sub_small(1);
    let mut s = 0;
    while !nm1.bit(s) {
        s += 1;
    }
    let d = nm1.shr(s);
    let mut minus_one = n.0.clone();
    MontRt::sub_in(&mut minus_one, &one);
    'outer: for &base in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        // base R mod n = base * (R mod n) by additions
        let mut b = vec![0u64; k];
        for _ in 0..base {
            let mut c = false;
            for i in 0..k {
                let (s1, o1) = b[i].overflowing_add(one[i]);
                let (s2, o2) = s1.overflowing_add(c as u64);
                b[i] = s2;
                c = o1 || o2;
            }
            if c || MontRt::geq(&b, &mt.n) {
                MontRt::sub_in(&mut b, &mt.n);
            }
        }
        let mut x = one.clone();
        for i in (0..d.bits()).rev() {
            x = mt.mul(&x, &x);
            if d.bit(i) {
                x = mt.mul(&x, &b);
            }
        }
        if x == one || x == minus_one {
            continue;
        }
        for _ in 0..s - 1 {
            x = mt.mul(&x, &x);
            if x == minus_one {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

/// a mod n by shift-and-subtract on preallocated limbs (setup only).
pub fn big_mod(a: &Big, n: &Big) -> Big {
    if a < n {
        return a.clone();
    }
    let k = n.0.len();
    let mut r = vec![0u64; k + 1];
    let mut nn = n.0.clone();
    nn.push(0);
    for i in (0..a.bits()).rev() {
        let mut c = a.bit(i) as u64;
        for x in r.iter_mut() {
            let nc = *x >> 63;
            *x = (*x << 1) | c;
            c = nc;
        }
        if MontRt::geq(&r, &nn) {
            MontRt::sub_in(&mut r, &nn);
        }
    }
    Big::from_limbs(&r)
}
