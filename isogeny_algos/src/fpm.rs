//! Prime fields in Montgomery representation with N 64-bit limbs (`FpM<1>` for primes below 2^63,
//! `FpM<4>` for 256-bit, `FpM<8>` for 512-bit). Elements are `[u64; N]` in Montgomery form
//! a R mod p (R = 2^{64N}); equality and hashing work on that form. The simple reference field
//! remains `field::Zp` (canonical u64 residues).
use crate::bigint::{big_mod, Big};
use crate::field::{Field, Rng};

const MAXN: usize = 16;

#[derive(Clone, Debug)]
pub struct FpM<const N: usize> {
    pub p: [u64; N],
    pinv: u64,
    one: [u64; N],
    r2: [u64; N],
    r3: [u64; N],
    pbig: Big,
    pbits: u32,
    fast1: bool,
}

#[inline(always)]
fn add_raw<const N: usize>(a: &[u64; N], b: &[u64; N]) -> ([u64; N], bool) {
    let mut r = [0u64; N];
    let mut c = false;
    for i in 0..N {
        let (s1, o1) = a[i].overflowing_add(b[i]);
        let (s2, o2) = s1.overflowing_add(c as u64);
        r[i] = s2;
        c = o1 || o2;
    }
    (r, c)
}
#[inline(always)]
fn sub_raw<const N: usize>(a: &[u64; N], b: &[u64; N]) -> ([u64; N], bool) {
    let mut r = [0u64; N];
    let mut br = false;
    for i in 0..N {
        let (d1, o1) = a[i].overflowing_sub(b[i]);
        let (d2, o2) = d1.overflowing_sub(br as u64);
        r[i] = d2;
        br = o1 || o2;
    }
    (r, br)
}
#[inline(always)]
fn geq<const N: usize>(a: &[u64; N], b: &[u64; N]) -> bool {
    for i in (0..N).rev() {
        if a[i] != b[i] {
            return a[i] > b[i];
        }
    }
    true
}

/// Inverse of a in [1, p) modulo p < 2^63: extended Euclid on u64 with i64 Bezout coefficients
/// (every true coefficient is bounded by p, so wrapping two's-complement arithmetic is exact).
/// Measured 2x faster than the i128 version and 2.7x faster than binary gcd on this machine.
pub fn inv_u64(a: u64, p: u64) -> u64 {
    let (mut t, mut nt) = (0i64, 1i64);
    let (mut r, mut nr) = (p, a);
    while nr != 0 {
        let q = r / nr;
        let tmp = t.wrapping_sub((q as i64).wrapping_mul(nt));
        t = nt;
        nt = tmp;
        let tmp2 = r - q * nr;
        r = nr;
        nr = tmp2;
    }
    if t < 0 {
        (t + p as i64) as u64
    } else {
        t as u64
    }
}

impl<const N: usize> FpM<N> {
    pub fn new(p: &Big) -> Self {
        assert!(N >= 1 && N <= MAXN);
        assert!(p.bit(0) && p.bits() >= 3, "p must be odd and > 3");
        assert!(p.bits() <= 64 * N, "p does not fit in {N} limbs");
        if N > 1 {
            assert!(p.bits() > 64 * (N - 1), "use fewer limbs for this p");
        }
        let pl: [u64; N] = p.limbs(N).try_into().unwrap();
        let mut inv = 1u64;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(pl[0].wrapping_mul(inv)));
        }
        let mut f = FpM {
            p: pl,
            pinv: inv.wrapping_neg(),
            one: [0; N],
            r2: [0; N],
            r3: [0; N],
            pbig: p.clone(),
            pbits: p.bits() as u32,
            fast1: N == 1 && p.bits() <= 63,
        };
        // R mod p, R^2 mod p, R^3 mod p by repeated doubling of 1
        let mut x = [0u64; N];
        x[0] = 1;
        let mut out = vec![];
        for _ in 0..3 {
            for _ in 0..64 * N {
                x = f.add_mod(&x, &x);
            }
            out.push(x);
        }
        f.one = out[0];
        f.r2 = out[1];
        f.r3 = out[2];
        f
    }
    pub fn from_u64_modulus(p: u64) -> Self {
        FpM::new(&Big::from_u64(p))
    }
    pub fn from_dec(s: &str) -> Self {
        FpM::new(&Big::from_dec(s))
    }
    pub fn modulus(&self) -> &Big {
        &self.pbig
    }
    #[inline(always)]
    fn add_mod(&self, a: &[u64; N], b: &[u64; N]) -> [u64; N] {
        let (s, c) = add_raw(a, b);
        if c || geq(&s, &self.p) {
            sub_raw(&s, &self.p).0
        } else {
            s
        }
    }
    #[inline(always)]
    fn sub_mod(&self, a: &[u64; N], b: &[u64; N]) -> [u64; N] {
        let (d, br) = sub_raw(a, b);
        if br {
            add_raw(&d, &self.p).0
        } else {
            d
        }
    }
    #[inline(always)]
    fn mont_mul(&self, a: &[u64; N], b: &[u64; N]) -> [u64; N] {
        if self.fast1 {
            let p0 = self.p[0];
            let x = a[0] as u128 * b[0] as u128;
            let m = (x as u64).wrapping_mul(self.pinv);
            let t = ((x + m as u128 * p0 as u128) >> 64) as u64;
            let mut r = [0u64; N];
            r[0] = if t >= p0 { t - p0 } else { t };
            return r;
        }
        let mut t = [0u64; MAXN + 2];
        for i in 0..N {
            let mut c: u128 = 0;
            for j in 0..N {
                let s = t[j] as u128 + a[j] as u128 * b[i] as u128 + c;
                t[j] = s as u64;
                c = s >> 64;
            }
            let s = t[N] as u128 + c;
            t[N] = s as u64;
            t[N + 1] = (s >> 64) as u64;
            let m = t[0].wrapping_mul(self.pinv);
            let mut c = (t[0] as u128 + m as u128 * self.p[0] as u128) >> 64;
            for j in 1..N {
                let s = t[j] as u128 + m as u128 * self.p[j] as u128 + c;
                t[j - 1] = s as u64;
                c = s >> 64;
            }
            let s = t[N] as u128 + c;
            t[N - 1] = s as u64;
            t[N] = t[N + 1] + (s >> 64) as u64;
        }
        let mut r = [0u64; N];
        r.copy_from_slice(&t[..N]);
        if t[N] != 0 || geq(&r, &self.p) {
            r = sub_raw(&r, &self.p).0;
        }
        r
    }
    /// Canonical residue (out of Montgomery form).
    pub fn to_canonical(&self, a: &[u64; N]) -> [u64; N] {
        let mut one = [0u64; N];
        one[0] = 1;
        self.mont_mul(a, &one)
    }
    pub fn to_big(&self, a: &[u64; N]) -> Big {
        Big::from_limbs(&self.to_canonical(a))
    }
    pub fn from_big(&self, v: &Big) -> [u64; N] {
        let c: [u64; N] = big_mod(v, &self.pbig).limbs(N).try_into().unwrap();
        self.mont_mul(&c, &self.r2)
    }
}

impl<const N: usize> Field for FpM<N> {
    type E = [u64; N];
    fn zero(&self) -> Self::E {
        [0u64; N]
    }
    fn one(&self) -> Self::E {
        self.one
    }
    #[inline(always)]
    fn add(&self, a: Self::E, b: Self::E) -> Self::E {
        self.add_mod(&a, &b)
    }
    #[inline(always)]
    fn sub(&self, a: Self::E, b: Self::E) -> Self::E {
        self.sub_mod(&a, &b)
    }
    #[inline(always)]
    fn neg(&self, a: Self::E) -> Self::E {
        if a == [0u64; N] {
            a
        } else {
            sub_raw(&self.p, &a).0
        }
    }
    #[inline(always)]
    fn mul(&self, a: Self::E, b: Self::E) -> Self::E {
        self.mont_mul(&a, &b)
    }
    fn inv(&self, a: Self::E) -> Self::E {
        assert!(a != [0u64; N], "inverse of zero");
        if N == 1 && self.pbits <= 63 {
            // (aR)^{-1} computed on integers, then * R^3 R^{-1} = a^{-1} R
            let c = inv_u64(a[0], self.p[0]);
            let mut x = [0u64; N];
            x[0] = c;
            return self.mont_mul(&x, &self.r3);
        }
        self.pow_big(a, &self.pbig.sub_small(2))
    }
    fn from_u64(&self, n: u64) -> Self::E {
        let mut c = [0u64; N];
        c[0] = if N == 1 { n % self.p[0] } else { n };
        self.mont_mul(&c, &self.r2)
    }
    fn char(&self) -> u64 {
        self.p[0]
    }
    fn conv(&self, a: &[Self::E], b: &[Self::E]) -> Vec<Self::E> {
        self.conv_trunc(a, b, a.len() + b.len() - 1)
    }
    fn dot_rev(&self, a: &[Self::E], b: &[Self::E]) -> Self::E {
        if N == 1 && self.pbits <= 62 {
            let a1: Vec<u64> = a.iter().map(|x| x[0]).collect();
            let b1: Vec<u64> = b.iter().map(|x| x[0]).collect();
            let v = crate::field::lazy_dot_rev(&a1, &b1, self.p[0]) as u64;
            let mut t = [0u64; N];
            t[0] = v;
            let mut one = [0u64; N];
            one[0] = 1;
            return self.mont_mul(&t, &one);
        }
        let n = a.len().min(b.len());
        let mut s = [0u64; N];
        for i in 0..n {
            s = self.add_mod(&s, &self.mont_mul(&a[i], &b[b.len() - 1 - i]));
        }
        s
    }
    fn conv_trunc(&self, a: &[Self::E], b: &[Self::E], nout: usize) -> Vec<Self::E> {
        if N == 1 && self.pbits <= 62 {
            // sum of (aR)(bR) reduced mod p, then one REDC: (sum ab) R
            let a1: Vec<u64> = a.iter().map(|x| x[0]).collect();
            let b1: Vec<u64> = b.iter().map(|x| x[0]).collect();
            let p0 = self.p[0];
            let v = crate::field::lazy_conv(&a1, &b1, p0, nout, |acc| acc as u64);
            return v
                .into_iter()
                .map(|x| {
                    let mut t = [0u64; N];
                    t[0] = x;
                    let mut one = [0u64; N];
                    one[0] = 1;
                    self.mont_mul(&t, &one)
                })
                .collect();
        }
        let mut r = vec![[0u64; N]; nout];
        for (i, x) in a.iter().enumerate().take(nout) {
            for (j, y) in b.iter().enumerate().take(nout - i) {
                r[i + j] = self.add_mod(&r[i + j], &self.mont_mul(x, y));
            }
        }
        r
    }
    fn size(&self) -> u128 {
        self.pbig
            .to_u128()
            .expect("FpM::size: modulus exceeds 128 bits; use q()")
    }
    fn q(&self) -> Big {
        self.pbig.clone()
    }
    fn random(&self, rng: &mut Rng) -> Self::E {
        let top_bits = self.pbits as usize - 64 * (N - 1);
        loop {
            let mut r = [0u64; N];
            for x in r.iter_mut() {
                *x = rng.next();
            }
            if top_bits < 64 {
                r[N - 1] &= (1u64 << top_bits) - 1;
            }
            if !geq(&r, &self.p) {
                return r;
            }
        }
    }
}
