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
    /// p[N-1] < 2^63 - 2: the no-carry CIOS variant applies (top bit of the top limb free).
    nocarry: bool,
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
            nocarry: N > 1 && pl[N - 1] < (1u64 << 63) - 2,
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
        if N == 1 && self.fast1 {
            let p0 = self.p[0];
            let x = a[0] as u128 * b[0] as u128;
            let m = (x as u64).wrapping_mul(self.pinv);
            let t = ((x + m as u128 * p0 as u128) >> 64) as u64;
            let mut r = [0u64; N];
            r[0] = if t >= p0 { t - p0 } else { t };
            return r;
        }
        if self.nocarry {
            // CIOS without the overflow word (Gautam-Botrel-Piellard, gnark): when the top limb of
            // p is below 2^63 - 2 the running value always fits in N words. Measured ~14% faster
            // than the general loop at N = 8.
            let mut t = [0u64; N];
            for i in 0..N {
                let bi = b[i] as u128;
                let s = t[0] as u128 + a[0] as u128 * bi;
                let mut ca = (s >> 64) as u64;
                let t0 = s as u64;
                let m = t0.wrapping_mul(self.pinv) as u128;
                let mut cb = ((t0 as u128 + m * self.p[0] as u128) >> 64) as u64;
                for j in 1..N {
                    let s = t[j] as u128 + a[j] as u128 * bi + ca as u128;
                    ca = (s >> 64) as u64;
                    let s2 = (s as u64) as u128 + m * self.p[j] as u128 + cb as u128;
                    cb = (s2 >> 64) as u64;
                    t[j - 1] = s2 as u64;
                }
                t[N - 1] = ca.wrapping_add(cb);
            }
            if geq(&t, &self.p) {
                t = sub_raw(&t, &self.p).0;
            }
            return t;
        }
        // CIOS with the running value in t[0..N] plus one overflow word (t_n); every inner loop
        // has constant trip count N, so it is fully unrolled per instantiation.
        let mut t = [0u64; N];
        let mut t_n = 0u64;
        for i in 0..N {
            let bi = b[i] as u128;
            let mut c = 0u64;
            for j in 0..N {
                let s = t[j] as u128 + a[j] as u128 * bi + c as u128;
                t[j] = s as u64;
                c = (s >> 64) as u64;
            }
            let s = t_n as u128 + c as u128;
            let (t_n0, t_n1) = (s as u64, (s >> 64) as u64);
            let m = t[0].wrapping_mul(self.pinv) as u128;
            let mut c = ((t[0] as u128 + m * self.p[0] as u128) >> 64) as u64;
            for j in 1..N {
                let s = t[j] as u128 + m * self.p[j] as u128 + c as u128;
                t[j - 1] = s as u64;
                c = (s >> 64) as u64;
            }
            let s = t_n0 as u128 + c as u128;
            t[N - 1] = s as u64;
            t_n = t_n1 + (s >> 64) as u64;
        }
        if t_n != 0 || geq(&t, &self.p) {
            t = sub_raw(&t, &self.p).0;
        }
        t
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

/// Signed (two's complement) value on N+1 limbs, the scratch type of the binary GCD.
type Wide = [u64; MAXN + 1];

#[inline(always)]
fn lincomb<const N: usize>(x: &[u64; N], fx: i64, y: &[u64; N], fy: i64) -> Wide {
    let mut r = [0u64; MAXN + 1];
    let mut carry: i128 = 0;
    for i in 0..N {
        let s = x[i] as i128 * fx as i128 + y[i] as i128 * fy as i128 + carry;
        r[i] = s as u64;
        carry = s >> 64;
    }
    r[N] = carry as u64;
    r
}
#[inline(always)]
fn wide_neg(r: &mut Wide, n: usize) {
    let mut c = true;
    for x in r.iter_mut().take(n) {
        let (v, o) = (!*x).overflowing_add(c as u64);
        *x = v;
        c = o;
    }
}
/// Arithmetic shift right by k < 64 of an (n)-limb signed value.
#[inline(always)]
fn wide_shr(r: &mut Wide, n: usize, k: u32) {
    for i in 0..n - 1 {
        r[i] = (r[i] >> k) | (r[i + 1] << (64 - k));
    }
    r[n - 1] = ((r[n - 1] as i64) >> k) as u64;
}
#[inline(always)]
fn wide_add_mul_small<const N: usize>(r: &mut Wide, m: &[u64; N], t: u64) {
    let mut c: u128 = 0;
    for i in 0..N {
        let s = r[i] as u128 + m[i] as u128 * t as u128 + c;
        r[i] = s as u64;
        c = s >> 64;
    }
    r[N] = r[N].wrapping_add(c as u64);
}
#[inline(always)]
fn bitlen<const N: usize>(a: &[u64; N]) -> usize {
    for i in (0..N).rev() {
        if a[i] != 0 {
            return 64 * i + 64 - a[i].leading_zeros() as usize;
        }
    }
    0
}
/// Bits [s, s+32) of a (s may exceed the length; missing bits are zero).
#[inline(always)]
fn bits32<const N: usize>(a: &[u64; N], s: usize) -> u64 {
    let (w, o) = (s / 64, s % 64);
    let lo = if w < N { a[w] >> o } else { 0 };
    let hi = if o > 32 && w + 1 < N { a[w + 1] << (64 - o) } else { 0 };
    (lo | hi) & 0xFFFF_FFFF
}

impl<const N: usize> FpM<N> {
    /// Plain-integer inverse of y in [1, p) by Pornin's optimized binary GCD (ePrint 2020/972):
    /// 31 binary-GCD steps at a time on 62-bit approximations of (a, b) (top 32 and low 30 bits),
    /// the accumulated 2x2 update applied once to the full-length values. Invariants a = u y,
    /// b = v y (mod p); ends with b = gcd = 1, v = 1/y. Not constant time.
    fn bingcd_inv(&self, y: &[u64; N]) -> [u64; N] {
        const K: u32 = 31;
        let p = &self.p;
        // -p^{-1} mod 2^(K-1)
        let minv = self.pinv & ((1u64 << (K - 1)) - 1);
        let (mut a, mut b) = (*y, *p);
        let mut u = [0u64; N];
        u[0] = 1;
        let mut v = [0u64; N];
        let n1 = N + 1;
        // (u f + v g) / 2^(K-1) mod p, reduced into [0, p)
        let upd = |u: &[u64; N], f: i64, v: &[u64; N], g: i64| -> [u64; N] {
            let mut w = lincomb(u, f, v, g);
            let t = (w[0].wrapping_mul(minv)) & ((1u64 << (K - 1)) - 1);
            wide_add_mul_small(&mut w, p, t);
            wide_shr(&mut w, n1, K - 1);
            // |w| < 2p: fold into [0, p)
            let mut pw = [0u64; MAXN + 1];
            pw[..N].copy_from_slice(p);
            while (w[N] as i64) < 0 {
                let mut c = false;
                for i in 0..n1 {
                    let (s1, o1) = w[i].overflowing_add(pw[i]);
                    let (s2, o2) = s1.overflowing_add(c as u64);
                    w[i] = s2;
                    c = o1 || o2;
                }
            }
            loop {
                let mut lo = [0u64; N];
                lo.copy_from_slice(&w[..N]);
                if w[N] == 0 && !geq(&lo, p) {
                    return lo;
                }
                let mut br = false;
                for i in 0..n1 {
                    let (d1, o1) = w[i].overflowing_sub(pw[i]);
                    let (d2, o2) = d1.overflowing_sub(br as u64);
                    w[i] = d2;
                    br = o1 || o2;
                }
            }
        };
        let cap = 4 * (64 * N) / (K as usize - 1) + 8;
        for _ in 0..cap {
            if a.iter().all(|&x| x == 0) {
                break;
            }
            let n = bitlen(&a).max(bitlen(&b)).max(2 * K as usize);
            let low = (1u64 << (K - 1)) - 1;
            let mut aa = (a[0] & low) | (bits32(&a, n - K as usize - 1) << (K - 1));
            let mut bb = (b[0] & low) | (bits32(&b, n - K as usize - 1) << (K - 1));
            let (mut f0, mut g0, mut f1, mut g1) = (1i64, 0i64, 0i64, 1i64);
            for _ in 0..K - 1 {
                if aa & 1 == 1 {
                    if aa < bb {
                        std::mem::swap(&mut aa, &mut bb);
                        std::mem::swap(&mut f0, &mut f1);
                        std::mem::swap(&mut g0, &mut g1);
                    }
                    aa -= bb;
                    f0 -= f1;
                    g0 -= g1;
                }
                aa >>= 1;
                f1 <<= 1;
                g1 <<= 1;
            }
            let mut na = lincomb(&a, f0, &b, g0);
            let mut nb = lincomb(&a, f1, &b, g1);
            wide_shr(&mut na, n1, K - 1);
            wide_shr(&mut nb, n1, K - 1);
            if (na[N] as i64) < 0 {
                wide_neg(&mut na, n1);
                f0 = -f0;
                g0 = -g0;
            }
            if (nb[N] as i64) < 0 {
                wide_neg(&mut nb, n1);
                f1 = -f1;
                g1 = -g1;
            }
            a.copy_from_slice(&na[..N]);
            b.copy_from_slice(&nb[..N]);
            let (nu, nv) = (upd(&u, f0, &v, g0), upd(&u, f1, &v, g1));
            u = nu;
            v = nv;
        }
        let mut one = [0u64; N];
        one[0] = 1;
        if a.iter().any(|&x| x != 0) || b != one {
            // no convergence within the cap (not expected): fall back to Fermat on R y
            let yr = self.mont_mul(y, &self.r2);
            let inv_m = self.pow_big(yr, &self.pbig.sub_small(2));
            // inv_m = y^{-1} R; return the plain integer y^{-1}
            return self.mont_mul(&inv_m, &one);
        }
        v
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
        let c = self.bingcd_inv(&a);
        self.mont_mul(&c, &self.r3)
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
