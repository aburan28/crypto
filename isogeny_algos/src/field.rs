//! Field abstraction. Baseline arithmetic is deliberately naive
//! (u128 `%` reduction, extended-Euclid inversion) so later iterations have
//! something to beat.
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
        let q = self.size();
        let half = (q - 1) / 2;
        if self.pow(a, half) != self.one() {
            return None;
        }
        let mut rng = Rng::new(0xC0FFEE);
        let z = loop {
            let z = self.random(&mut rng);
            if !self.is_zero(z) && self.pow(z, half) != self.one() {
                break z;
            }
        };
        let (mut s, mut m) = (0u32, q - 1);
        while m % 2 == 0 {
            m /= 2;
            s += 1;
        }
        let mut c = self.pow(z, m);
        let mut t = self.pow(a, m);
        let mut r = self.pow(a, m.div_ceil(2));
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
        let (mut t, mut nt) = (0i128, 1i128);
        let (mut r, mut nr) = (self.p as i128, a as i128);
        while nr != 0 {
            let q = r / nr;
            (t, nt) = (nt, t - q * nt);
            (r, nr) = (nr, r - q * nr);
        }
        if t < 0 {
            t += self.p as i128;
        }
        t as u64
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
