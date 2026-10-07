//! Field abstraction. Baseline arithmetic is deliberately naive
//! (u128 `%` reduction, extended-Euclid inversion) so later iterations have
//! something to beat.
use crate::bigint::Big;
use std::fmt::Debug;
use std::hash::Hash;

#[derive(Clone)]
pub struct Rng(u64);
impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) ^ 0xD1B54A32D192ED03)
    }
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    pub fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

pub trait Field: Clone + Send + Sync + 'static {
    type E: Copy + Eq + Hash + Debug + Send + Sync;
    fn zero(&self) -> Self::E;
    fn one(&self) -> Self::E;
    fn add(&self, a: Self::E, b: Self::E) -> Self::E;
    fn sub(&self, a: Self::E, b: Self::E) -> Self::E;
    fn neg(&self, a: Self::E) -> Self::E;
    fn mul(&self, a: Self::E, b: Self::E) -> Self::E;
    /// Inverse; panics on zero.
    fn inv(&self, a: Self::E) -> Self::E;
    fn from_u64(&self, n: u64) -> Self::E;
    fn char(&self) -> u64;
    /// Field size q (p or p^2).
    fn size(&self) -> u128;
    fn random(&self, rng: &mut Rng) -> Self::E;
    fn is_zero(&self, a: Self::E) -> bool {
        a == self.zero()
    }
    fn sq(&self, a: Self::E) -> Self::E {
        self.mul(a, a)
    }
    /// Operand length below which polynomial products are schoolbook (`conv`).
    fn karatsuba_threshold(&self) -> usize {
        32
    }
    fn pow(&self, a: Self::E, mut e: u128) -> Self::E {
        let mut r = self.one();
        let mut b = a;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, b);
            }
            b = self.mul(b, b);
            e >>= 1;
        }
        r
    }
    /// Field size as a big integer (defaults to `size()`; big fields override).
    fn q(&self) -> Big {
        Big::from_u128(self.size())
    }
    /// a^e by left-to-right sliding windows over the odd powers a, a^3, ..., a^(2^w - 1)
    /// (w = 5 above 256 bits): ~bits/(w+1) multiplications instead of ~bits/2.
    fn pow_big(&self, a: Self::E, e: &Big) -> Self::E {
        let nb = e.bits();
        if nb <= 64 {
            return self.pow(a, e.to_u128().unwrap());
        }
        let w = if nb > 256 { 5 } else if nb > 128 { 4 } else { 3 };
        let a2 = self.sq(a);
        let mut tbl = Vec::with_capacity(1 << (w - 1));
        tbl.push(a);
        for k in 1..(1usize << (w - 1)) {
            tbl.push(self.mul(tbl[k - 1], a2));
        }
        let mut r = self.one();
        let mut started = false;
        let mut i = nb as isize - 1;
        while i >= 0 {
            if !e.bit(i as usize) {
                if started {
                    r = self.sq(r);
                }
                i -= 1;
                continue;
            }
            let mut j = (i - w as isize + 1).max(0);
            while !e.bit(j as usize) {
                j += 1;
            }
            let mut val = 0usize;
            for k in (j..=i).rev() {
                val = 2 * val + e.bit(k as usize) as usize;
            }
            if started {
                for _ in j..=i {
                    r = self.sq(r);
                }
                r = self.mul(r, tbl[(val - 1) / 2]);
            } else {
                r = tbl[(val - 1) / 2];
                started = true;
            }
            i = j - 1;
        }
        r
    }
    /// Full product of two coefficient vectors (length a.len()+b.len()-1). Default: schoolbook;
    /// single-word fields override it with lazy (accumulate-then-reduce) arithmetic.
    fn conv(&self, a: &[Self::E], b: &[Self::E]) -> Vec<Self::E> {
        let mut r = vec![self.zero(); a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            if self.is_zero(x) {
                continue;
            }
            for (j, &y) in b.iter().enumerate() {
                r[i + j] = self.add(r[i + j], self.mul(x, y));
            }
        }
        r
    }
    /// First `n` coefficients of the product (truncated power-series multiplication).
    fn conv_trunc(&self, a: &[Self::E], b: &[Self::E], n: usize) -> Vec<Self::E> {
        let mut r = vec![self.zero(); n];
        for i in 0..a.len().min(n) {
            if self.is_zero(a[i]) {
                continue;
            }
            for j in 0..b.len().min(n - i) {
                r[i + j] = self.add(r[i + j], self.mul(a[i], b[j]));
            }
        }
        r
    }
    /// sum_i a[i] * b[len-1-i] over the common length (reversed dot product).
    fn dot_rev(&self, a: &[Self::E], b: &[Self::E]) -> Self::E {
        let n = a.len().min(b.len());
        let mut s = self.zero();
        for i in 0..n {
            s = self.add(s, self.mul(a[i], b[b.len() - 1 - i]));
        }
        s
    }
    /// Inverses of all non-zero entries with one field inversion (Montgomery's trick).
    fn batch_inv(&self, xs: &[Self::E]) -> Vec<Self::E> {
        let n = xs.len();
        let mut pre = Vec::with_capacity(n);
        let mut acc = self.one();
        for &x in xs {
            pre.push(acc);
            acc = self.mul(acc, x);
        }
        let mut inv = self.inv(acc);
        let mut out = vec![self.zero(); n];
        for i in (0..n).rev() {
            out[i] = self.mul(inv, pre[i]);
            inv = self.mul(inv, xs[i]);
        }
        out
    }
    fn from_i64(&self, n: i64) -> Self::E {
        if n >= 0 {
            self.from_u64(n as u64)
        } else {
            self.neg(self.from_u64((-n) as u64))
        }
    }
    fn div(&self, a: Self::E, b: Self::E) -> Self::E {
        self.mul(a, self.inv(b))
    }
    /// Square root in a field of odd size (generic Tonelli-Shanks); None for non-residues.
    fn sqrt(&self, a: Self::E) -> Option<Self::E> {
        if self.is_zero(a) {
            return Some(a);
        }
        let q = self.q();
        if q.bit(0) && q.bit(1) {
            // q = 3 mod 4: one exponentiation, then check
            let r = self.pow_big(a, &q.add_small(1).shr(2));
            return if self.sq(r) == a { Some(r) } else { None };
        }
        let qm1 = q.sub_small(1);
        let half = qm1.shr(1);
        if self.pow_big(a, &half) != self.one() {
            return None;
        }
        let mut rng = Rng::new(0xC0FFEE);
        let z = loop {
            let z = self.random(&mut rng);
            if !self.is_zero(z) && self.pow_big(z, &half) != self.one() {
                break z;
            }
        };
        let mut s = 0u32;
        while !qm1.bit(s as usize) {
            s += 1;
        }
        let m = qm1.shr(s as usize);
        let mut c = self.pow_big(z, &m);
        let mut t = self.pow_big(a, &m);
        let mut r = self.pow_big(a, &m.add_small(1).shr(1));
        let mut mm = s;
        while t != self.one() {
            let (mut i, mut tt) = (0u32, t);
            while tt != self.one() {
                tt = self.sq(tt);
                i += 1;
            }
            let mut b = c;
            for _ in 0..(mm - i - 1) {
                b = self.sq(b);
            }
            mm = i;
            c = self.sq(b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
        Some(r)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Zp {
    pub p: u64,
}

impl Zp {
    pub fn new(p: u64) -> Self {
        assert!(p > 3 && p < (1u64 << 62));
        Zp { p }
    }
    pub fn legendre(&self, a: u64) -> i32 {
        if a == 0 {
            return 0;
        }
        if self.pow(a, ((self.p - 1) / 2) as u128) == 1 {
            1
        } else {
            -1
        }
    }
    /// Tonelli-Shanks.
    pub fn sqrt(&self, a: u64) -> Option<u64> {
        if a == 0 {
            return Some(0);
        }
        if self.legendre(a) != 1 {
            return None;
        }
        let p = self.p;
        if p % 4 == 3 {
            return Some(self.pow(a, ((p + 1) / 4) as u128));
        }
        let mut q = p - 1;
        let mut s = 0;
        while q % 2 == 0 {
            q /= 2;
            s += 1;
        }
        let mut z = 2;
        while self.legendre(z) != -1 {
            z += 1;
        }
        let mut m = s;
        let mut c = self.pow(z, q as u128);
        let mut t = self.pow(a, q as u128);
        let mut r = self.pow(a, ((q + 1) / 2) as u128);
        while t != 1 {
            let mut i = 0;
            let mut tt = t;
            while tt != 1 {
                tt = self.mul(tt, tt);
                i += 1;
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = self.mul(b, b);
            }
            m = i;
            c = self.mul(b, b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
        Some(r)
    }
}

impl Field for Zp {
    type E = u64;
    fn zero(&self) -> u64 {
        0
    }
    fn one(&self) -> u64 {
        1
    }
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    fn mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % (self.p as u128)) as u64
    }
    fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "inverse of zero");
        crate::fpm::inv_u64(a, self.p)
    }
    fn from_u64(&self, n: u64) -> u64 {
        n % self.p
    }
    fn char(&self) -> u64 {
        self.p
    }
    fn size(&self) -> u128 {
        self.p as u128
    }
    fn random(&self, rng: &mut Rng) -> u64 {
        rng.below(self.p)
    }
    fn sqrt(&self, a: u64) -> Option<u64> {
        Zp::sqrt(self, a)
    }
    fn conv(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        lazy_conv(a, b, self.p, a.len() + b.len() - 1, |acc| acc as u64)
    }
    fn conv_trunc(&self, a: &[u64], b: &[u64], n: usize) -> Vec<u64> {
        lazy_conv(a, b, self.p, n, |acc| acc as u64)
    }
    fn dot_rev(&self, a: &[u64], b: &[u64]) -> u64 {
        lazy_dot_rev(a, b, self.p) as u64
    }
}

/// sum a[i] b[len-1-i] mod p with lazy u128 accumulation (p < 2^62).
pub fn lazy_dot_rev(a: &[u64], b: &[u64], p: u64) -> u128 {
    let n = a.len().min(b.len());
    let pp = p as u128;
    let mut acc: u128 = 0;
    let mut cnt = 0;
    let bl = b.len();
    for i in 0..n {
        acc += a[i] as u128 * b[bl - 1 - i] as u128;
        cnt += 1;
        if cnt == 15 {
            acc %= pp;
            cnt = 1;
        }
    }
    acc % pp
}

/// Convolution with u128 accumulation: products are < p^2 < 2^124 for p < 2^62, so 15 of them can be
/// summed before reducing; `fin` maps the reduced accumulator (< p) to the output representation.
pub fn lazy_conv(a: &[u64], b: &[u64], p: u64, nout: usize, fin: impl Fn(u128) -> u64) -> Vec<u64> {
    let (n, m) = (a.len(), b.len());
    let mut r = Vec::with_capacity(nout);
    let pp = p as u128;
    if n == 0 || m == 0 {
        return vec![fin(0); nout];
    }
    for k in 0..nout {
        if k > n + m - 2 {
            r.push(fin(0));
            continue;
        }
        let lo = k.saturating_sub(m - 1);
        let hi = k.min(n - 1);
        let mut acc: u128 = 0;
        let mut cnt = 0;
        for i in lo..=hi {
            acc += a[i] as u128 * b[k - i] as u128;
            cnt += 1;
            if cnt == 15 {
                acc %= pp;
                cnt = 1;
            }
        }
        r.push(fin(acc % pp));
    }
    r
}

/// F_{p^2} = F_p[i]/(i^2 - nr), nr a non-residue.
#[derive(Clone, Copy, Debug)]
pub struct Zp2 {
    pub base: Zp,
    pub nr: u64,
}

impl Zp2 {
    pub fn new(p: u64) -> Self {
        let base = Zp::new(p);
        let mut nr = p - 1;
        if base.legendre(nr) != -1 {
            nr = 2;
            while base.legendre(nr) != -1 {
                nr += 1;
            }
        }
        Zp2 { base, nr }
    }
    pub fn embed(&self, a: u64) -> (u64, u64) {
        (a, 0)
    }
}

impl Field for Zp2 {
    type E = (u64, u64);
    fn zero(&self) -> Self::E {
        (0, 0)
    }
    fn one(&self) -> Self::E {
        (1, 0)
    }
    fn add(&self, a: Self::E, b: Self::E) -> Self::E {
        (self.base.add(a.0, b.0), self.base.add(a.1, b.1))
    }
    fn sub(&self, a: Self::E, b: Self::E) -> Self::E {
        (self.base.sub(a.0, b.0), self.base.sub(a.1, b.1))
    }
    fn neg(&self, a: Self::E) -> Self::E {
        (self.base.neg(a.0), self.base.neg(a.1))
    }
    fn mul(&self, a: Self::E, b: Self::E) -> Self::E {
        let f = &self.base;
        let t0 = f.mul(a.0, b.0);
        let t1 = f.mul(a.1, b.1);
        let c0 = f.add(t0, f.mul(self.nr, t1));
        let c1 = f.add(f.mul(a.0, b.1), f.mul(a.1, b.0));
        (c0, c1)
    }
    fn inv(&self, a: Self::E) -> Self::E {
        let f = &self.base;
        let n = f.sub(f.mul(a.0, a.0), f.mul(self.nr, f.mul(a.1, a.1)));
        let ni = f.inv(n);
        (f.mul(a.0, ni), f.neg(f.mul(a.1, ni)))
    }
    fn from_u64(&self, n: u64) -> Self::E {
        (n % self.base.p, 0)
    }
    fn char(&self) -> u64 {
        self.base.p
    }
    fn size(&self) -> u128 {
        (self.base.p as u128) * (self.base.p as u128)
    }
    fn random(&self, rng: &mut Rng) -> Self::E {
        (rng.below(self.base.p), rng.below(self.base.p))
    }
}

/// Deterministic Miller-Rabin for u64.
pub fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &q in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 {
            return n == q;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    let mulm = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    'w: for &a in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = 1u64;
        let mut b = a % n;
        let mut e = d;
        while e > 0 {
            if e & 1 == 1 {
                x = mulm(x, b);
            }
            b = mulm(b, b);
            e >>= 1;
        }
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..s - 1 {
            x = mulm(x, x);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}

pub fn next_prime(mut n: u64) -> u64 {
    if n % 2 == 0 {
        n += 1;
    }
    while !is_prime(n) {
        n += 2;
    }
    n
}

/// Prime factorisation of n < 2^64 (trial division to 1000, then Pollard-Brent rho), as
/// (prime, exponent) pairs in increasing order.
pub fn factor_u64(n: u64) -> Vec<(u64, u32)> {
    fn mulmod(a: u64, b: u64, m: u64) -> u64 {
        (a as u128 * b as u128 % m as u128) as u64
    }
    fn is_prime64(n: u64) -> bool {
        if n < 2 {
            return false;
        }
        for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
            if n % p == 0 {
                return n == p;
            }
        }
        let (mut d, mut s) = (n - 1, 0);
        while d % 2 == 0 {
            d /= 2;
            s += 1;
        }
        'outer: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
            let mut x = 1u64;
            let (mut b, mut e) = (a % n, d);
            while e > 0 {
                if e & 1 == 1 {
                    x = mulmod(x, b, n);
                }
                b = mulmod(b, b, n);
                e >>= 1;
            }
            if x == 1 || x == n - 1 {
                continue;
            }
            for _ in 1..s {
                x = mulmod(x, x, n);
                if x == n - 1 {
                    continue 'outer;
                }
            }
            return false;
        }
        true
    }
    fn gcd(mut a: u64, mut b: u64) -> u64 {
        while b != 0 {
            let t = a % b;
            a = b;
            b = t;
        }
        a
    }
    fn rho(n: u64) -> u64 {
        if n % 2 == 0 {
            return 2;
        }
        let mut c = 1u64;
        loop {
            let f = |x: u64| (mulmod(x, x, n) + c) % n;
            let (mut y, mut r, mut q, mut g) = (2u64, 1u64, 1u64, 1u64);
            let (mut x, mut ys) = (0u64, 0u64);
            while g == 1 {
                x = y;
                for _ in 0..r {
                    y = f(y);
                }
                let mut k = 0;
                while k < r && g == 1 {
                    ys = y;
                    for _ in 0..(128.min(r - k)) {
                        y = f(y);
                        q = mulmod(q, x.abs_diff(y), n);
                    }
                    g = gcd(q, n);
                    k += 128;
                }
                r *= 2;
            }
            if g == n {
                loop {
                    ys = f(ys);
                    g = gcd(x.abs_diff(ys), n);
                    if g > 1 {
                        break;
                    }
                }
            }
            if g != n {
                return g;
            }
            c += 1;
        }
    }
    fn rec(n: u64, out: &mut Vec<u64>) {
        if n == 1 {
            return;
        }
        if is_prime64(n) {
            out.push(n);
            return;
        }
        let d = rho(n);
        rec(d, out);
        rec(n / d, out);
    }
    let mut n = n;
    let mut ps = vec![];
    for p in 2..1000u64 {
        while n % p == 0 {
            ps.push(p);
            n /= p;
        }
    }
    rec(n, &mut ps);
    ps.sort();
    let mut out: Vec<(u64, u32)> = vec![];
    for p in ps {
        match out.last_mut() {
            Some((q, e)) if *q == p => *e += 1,
            _ => out.push((p, 1)),
        }
    }
    out
}
