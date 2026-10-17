//! GF(2^m) for m <= 64 (elements are u64 bit vectors, constant term in bit 0),
//! the X9.62 / FIPS 186 lowest-weight reduction polynomial, and the ordinary
//! binary Weierstrass curve y^2 + xy = x^3 + a2 x^2 + a6.

#[derive(Clone, Debug)]
pub struct F2m {
    pub m: u32,
    /// Exponents of the reduction polynomial, descending, e.g. [53, 6, 2, 1, 0].
    pub modulus: Vec<u32>,
    low: u64, // f - x^m
    mask: u64,
    trace_mask: u64,
    clmul_hw: bool,
}

#[inline]
fn clmul_sw(a: u64, b: u64) -> u128 {
    let mut r = 0u128;
    let a = a as u128;
    for i in 0..64 {
        let bit = ((b >> i) & 1) as u128;
        r ^= (a << i) & (0u128.wrapping_sub(bit));
    }
    r
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "pclmulqdq,sse2")]
unsafe fn clmul_hw(a: u64, b: u64) -> u128 {
    use std::arch::x86_64::*;
    let va = _mm_set_epi64x(0, a as i64);
    let vb = _mm_set_epi64x(0, b as i64);
    let r = _mm_clmulepi64_si128(va, vb, 0);
    let lo = _mm_cvtsi128_si64(r) as u64 as u128;
    let hi = _mm_cvtsi128_si64(_mm_unpackhi_epi64(r, r)) as u64 as u128;
    lo | (hi << 64)
}

/// Polynomial (over F_2) helpers on u128 used for irreducibility testing.
fn poly_mulmod_f2(a: u128, b: u128, f: u128, deg: u32) -> u128 {
    let mut r = 0u128;
    let mut a = a;
    let mut b = b;
    while b != 0 {
        if b & 1 == 1 {
            r ^= a;
        }
        b >>= 1;
        a <<= 1;
        if (a >> deg) & 1 == 1 {
            a ^= f;
        }
    }
    r
}

/// Rabin's test for prime-degree or general m (only m prime is used here).
pub fn is_irreducible(exps: &[u32]) -> bool {
    let m = exps[0];
    let mut f = 0u128;
    for &e in exps {
        f ^= 1u128 << e;
    }
    if f & 1 == 0 {
        return false;
    }
    // x^(2^k) mod f for k = 1..m
    let x = 2u128;
    let mut pows = vec![x];
    let mut cur = x;
    for _ in 0..m {
        cur = poly_mulmod_f2(cur, cur, f, m);
        pows.push(cur);
    }
    if pows[m as usize] != x {
        return false;
    }
    let mut primes = vec![];
    let mut n = m;
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            primes.push(d);
            while n % d == 0 {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        primes.push(n);
    }
    for r in primes {
        let g = pows[(m / r) as usize] ^ x;
        if poly_gcd_f2(g, f) != 1 {
            return false;
        }
    }
    true
}

fn poly_gcd_f2(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        while a != 0 && 127 - a.leading_zeros() >= 127 - b.leading_zeros() {
            let s = (127 - a.leading_zeros()) - (127 - b.leading_zeros());
            a ^= b << s;
        }
        std::mem::swap(&mut a, &mut b);
    }
    a
}

/// FIPS 186 / X9.62 rule: the trinomial x^m + x^k + 1 with least k, else the
/// pentanomial x^m + x^k3 + x^k2 + x^k1 + 1 with least k3, then k2, then k1.
pub fn standard_modulus(m: u32) -> Vec<u32> {
    for k in 1..m {
        let e = vec![m, k, 0];
        if is_irreducible(&e) {
            return e;
        }
    }
    for k3 in 3..m {
        for k2 in 2..k3 {
            for k1 in 1..k2 {
                let e = vec![m, k3, k2, k1, 0];
                if is_irreducible(&e) {
                    return e;
                }
            }
        }
    }
    panic!("no pentanomial for m={m}");
}

impl F2m {
    pub fn new(modulus: Vec<u32>) -> Self {
        let m = modulus[0];
        assert!((2..=64).contains(&m));
        assert!(is_irreducible(&modulus), "reducible modulus");
        let mut low = 0u64;
        for &e in &modulus[1..] {
            low ^= 1u64 << e;
        }
        let mask = if m == 64 { u64::MAX } else { (1u64 << m) - 1 };
        #[cfg(target_arch = "x86_64")]
        let clmul_hw = std::arch::is_x86_feature_detected!("pclmulqdq");
        #[cfg(not(target_arch = "x86_64"))]
        let clmul_hw = false;
        let mut k = F2m { m, modulus, low, mask, trace_mask: 0, clmul_hw };
        let mut tm = 0u64;
        for i in 0..m {
            if k.trace_slow(1u64 << i) == 1 {
                tm |= 1u64 << i;
            }
        }
        k.trace_mask = tm;
        k
    }

    #[inline]
    fn clmul(&self, a: u64, b: u64) -> u128 {
        #[cfg(target_arch = "x86_64")]
        if self.clmul_hw {
            return unsafe { clmul_hw(a, b) };
        }
        clmul_sw(a, b)
    }

    #[inline]
    fn reduce(&self, mut x: u128) -> u64 {
        let m = self.m;
        while (x >> m) != 0 {
            let h = (x >> m) as u64;
            x = (x & self.mask as u128) ^ self.clmul(h, self.low);
        }
        x as u64
    }

    #[inline]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.reduce(self.clmul(a, b))
    }

    #[inline]
    pub fn sqr(&self, a: u64) -> u64 {
        self.reduce(self.clmul(a, a))
    }

    pub fn pow2k(&self, mut a: u64, k: u32) -> u64 {
        for _ in 0..k {
            a = self.sqr(a);
        }
        a
    }

    /// Inverse by the binary extended Euclidean algorithm (Hankerson Alg. 2.48).
    pub fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "inverse of zero");
        let mut f = 0u128;
        for &e in &self.modulus {
            f ^= 1u128 << e;
        }
        let (mut u, mut v) = (a as u128, f);
        let (mut g1, mut g2) = (1u128, 0u128);
        while u != 1 {
            let du = 127 - u.leading_zeros() as i32;
            let dv = 127 - v.leading_zeros() as i32;
            let mut j = du - dv;
            if j < 0 {
                std::mem::swap(&mut u, &mut v);
                std::mem::swap(&mut g1, &mut g2);
                j = -j;
            }
            u ^= v << j;
            g1 ^= g2 << j;
        }
        // g1 may exceed degree m-1 by construction only transiently; reduce.
        self.reduce(g1)
    }

    fn trace_slow(&self, a: u64) -> u64 {
        let mut t = a;
        let mut s = a;
        for _ in 1..self.m {
            s = self.sqr(s);
            t ^= s;
        }
        assert!(t <= 1, "trace not in F_2");
        t
    }

    #[inline]
    pub fn trace(&self, a: u64) -> u64 {
        ((a & self.trace_mask).count_ones() & 1) as u64
    }

    /// A root z of z^2 + z = c for odd m when Tr(c) = 0.
    pub fn half_trace(&self, c: u64) -> u64 {
        assert!(self.m % 2 == 1);
        let mut h = c;
        let mut s = c;
        for _ in 0..(self.m - 1) / 2 {
            s = self.sqr(self.sqr(s));
            h ^= s;
        }
        h
    }

    pub fn modulus_hex(&self) -> String {
        let mut f = 0u128;
        for &e in &self.modulus {
            f ^= 1u128 << e;
        }
        format!("0x{f:x}")
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BPt {
    Inf,
    A(u64, u64),
}

/// y^2 + xy = x^3 + a2 x^2 + a6 over GF(2^m).
#[derive(Clone, Debug)]
pub struct BinCurve {
    pub k: F2m,
    pub a2: u64,
    pub a6: u64,
}

impl BinCurve {
    pub fn on_curve(&self, p: BPt) -> bool {
        match p {
            BPt::Inf => true,
            BPt::A(x, y) => {
                let k = &self.k;
                let lhs = k.sqr(y) ^ k.mul(x, y);
                let x2 = k.sqr(x);
                let rhs = k.mul(x2, x) ^ k.mul(self.a2, x2) ^ self.a6;
                lhs == rhs
            }
        }
    }

    pub fn neg(&self, p: BPt) -> BPt {
        match p {
            BPt::Inf => BPt::Inf,
            BPt::A(x, y) => BPt::A(x, x ^ y),
        }
    }

    pub fn dbl(&self, p: BPt) -> BPt {
        match p {
            BPt::Inf => BPt::Inf,
            BPt::A(x, y) => {
                if x == 0 {
                    return BPt::Inf;
                }
                let k = &self.k;
                let l = x ^ k.mul(y, k.inv(x));
                let x3 = k.sqr(l) ^ l ^ self.a2;
                let y3 = k.sqr(x) ^ k.mul(l ^ 1, x3);
                BPt::A(x3, y3)
            }
        }
    }

    pub fn add(&self, p: BPt, q: BPt) -> BPt {
        match (p, q) {
            (BPt::Inf, _) => q,
            (_, BPt::Inf) => p,
            (BPt::A(x1, y1), BPt::A(x2, y2)) => {
                if x1 == x2 {
                    if y1 == y2 {
                        return self.dbl(p);
                    }
                    return BPt::Inf;
                }
                let k = &self.k;
                let l = k.mul(y1 ^ y2, k.inv(x1 ^ x2));
                let x3 = k.sqr(l) ^ l ^ x1 ^ x2 ^ self.a2;
                let y3 = k.mul(l, x1 ^ x3) ^ x3 ^ y1;
                BPt::A(x3, y3)
            }
        }
    }

    pub fn mul(&self, p: BPt, n: u128) -> BPt {
        let mut r = BPt::Inf;
        if n == 0 {
            return r;
        }
        for i in (0..128 - n.leading_zeros()).rev() {
            r = self.dbl(r);
            if (n >> i) & 1 == 1 {
                r = self.add(r, p);
            }
        }
        r
    }

    /// Point with the given x (if any), using z^2 + z = x + a2 + a6/x^2.
    pub fn lift_x(&self, x: u64) -> Option<BPt> {
        let k = &self.k;
        if x == 0 {
            return None;
        }
        let c = x ^ self.a2 ^ k.mul(self.a6, k.inv(k.sqr(x)));
        if k.trace(c) != 0 {
            return None;
        }
        let z = k.half_trace(c);
        let p = BPt::A(x, k.mul(x, z));
        debug_assert!(self.on_curve(p));
        Some(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_moduli_match_tables() {
        // SEC 2 sect113 and the k0n7 registry modulus 0x83; Seroussi's table for 47 and 53.
        assert_eq!(standard_modulus(113), vec![113, 9, 0]);
        assert_eq!(standard_modulus(7), vec![7, 1, 0]);
        assert_eq!(standard_modulus(47), vec![47, 5, 0]);
        assert_eq!(standard_modulus(53), vec![53, 6, 2, 1, 0]);
    }

    #[test]
    fn field_inverse_and_half_trace() {
        let k = F2m::new(standard_modulus(53));
        let mut x = 0x1234_5678_9abc_defu64 & ((1 << 53) - 1);
        for _ in 0..200 {
            x = k.mul(x, x ^ 0x5a5a5) ^ 1;
            let i = k.inv(x);
            assert_eq!(k.mul(x, i), 1);
            if k.trace(x) == 0 {
                let z = k.half_trace(x);
                assert_eq!(k.sqr(z) ^ z, x);
            }
        }
    }
}
