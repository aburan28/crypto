//! Binary fields GF(2^n) = GF(2)[x]/(f), n <= 63, elements as u64 bit vectors (bit i = coefficient
//! of x^i). Multiplication by carry-less multiply (PCLMULQDQ when the CPU has it, else a
//! shift-and-xor loop) and reduction by the sparse tail of f; inversion by the binary-polynomial
//! extended Euclid algorithm (Hankerson-Menezes-Vanstone, Alg. 2.48); square roots by the
//! inverse Frobenius a^(2^(n-1)); z^2 + z = c through a precomputed echelon form of the GF(2)-linear
//! map z -> z^2 + z; polynomial products accumulate unreduced 126-bit products and reduce once
//! per output coefficient.
use crate::bigint::Big;
use crate::field::{Field, Rng};

#[derive(Clone, Debug)]
pub struct GF2n {
    pub n: u32,
    /// f(x) - x^n (bits below n)
    pub red: u64,
    mask: u64,
    #[cfg(target_arch = "x86_64")]
    hw: bool,
    /// trace(x^i) as bit i: Tr(a) = parity(a & tmask)
    tmask: u64,
    /// echelon form of z -> z^2 + z: (image with pivot bit, preimage)
    qsolve: Vec<(u64, u64)>,
    /// the echelon solve is GF(2)-linear in c: byte tables of its values on each byte
    qtab: Vec<[u64; 256]>,
    /// sqrt(x), for sqrt(a) = even-part + sqrt(x) * odd-part
    sqrt_x: u64,
}

#[inline(always)]
fn clmul_soft(a: u64, b: u64) -> u128 {
    let mut r = 0u128;
    let a = a as u128;
    let mut b = b;
    let mut i = 0;
    while b != 0 {
        if b & 1 == 1 {
            r ^= a << i;
        }
        b >>= 1;
        i += 1;
    }
    r
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "pclmulqdq")]
unsafe fn clmul_hw(a: u64, b: u64) -> u128 {
    use std::arch::x86_64::*;
    let r = _mm_clmulepi64_si128(
        _mm_set_epi64x(0, a as i64),
        _mm_set_epi64x(0, b as i64),
        0x00,
    );
    std::mem::transmute::<__m128i, u128>(r)
}

/// Product and reduction in one call (the intrinsics inline here, not across the call).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "pclmulqdq")]
unsafe fn reduce_hw(mut t: u128, n: u32, mask: u64, red: u64) -> u64 {
    loop {
        let hi = (t >> n) as u64;
        if hi == 0 {
            return t as u64;
        }
        t = (t & mask as u128) ^ clmul_hw(hi, red);
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "pclmulqdq")]
unsafe fn mul_hw(a: u64, b: u64, n: u32, mask: u64, red: u64) -> u64 {
    reduce_hw(clmul_hw(a, b), n, mask, red)
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "pclmulqdq")]
unsafe fn dot_hw(
    a: &[u64],
    b: &[u64],
    lo: usize,
    hi: usize,
    k: usize,
    n: u32,
    mask: u64,
    red: u64,
) -> u64 {
    let mut acc = 0u128;
    for i in lo..=hi {
        acc ^= clmul_hw(a[i], b[k - i]);
    }
    reduce_hw(acc, n, mask, red)
}

/// Low-weight irreducible polynomials (x^n + tail): trinomials where one exists, else
/// pentanomials; `None` above the table (use `GF2n::with_tail`).
pub fn default_tail(n: u32) -> Option<u64> {
    // (n, [k]) for x^n + x^k + 1 or (n, [k1, k2, k3]) for x^n + x^k1 + x^k2 + x^k3 + 1
    let t: &[(u32, &[u32])] = &[
        (2, &[1]),
        (3, &[1]),
        (4, &[1]),
        (5, &[2]),
        (6, &[1]),
        (7, &[1]),
        (8, &[4, 3, 1]),
        (9, &[1]),
        (10, &[3]),
        (11, &[2]),
        (12, &[3]),
        (13, &[4, 3, 1]),
        (14, &[5]),
        (15, &[1]),
        (16, &[5, 3, 1]),
        (17, &[3]),
        (18, &[3]),
        (19, &[5, 2, 1]),
        (20, &[3]),
        (21, &[2]),
        (22, &[1]),
        (23, &[5]),
        (24, &[4, 3, 1]),
        (25, &[3]),
        (26, &[4, 3, 1]),
        (27, &[5, 2, 1]),
        (28, &[1]),
        (29, &[2]),
        (30, &[1]),
        (31, &[3]),
        (32, &[7, 3, 2]),
        (33, &[10]),
        (34, &[7]),
        (35, &[2]),
        (36, &[9]),
        (37, &[6, 4, 1]),
        (38, &[6, 5, 1]),
        (39, &[4]),
        (40, &[5, 4, 3]),
        (41, &[3]),
        (42, &[7]),
        (43, &[6, 4, 3]),
        (44, &[5]),
        (45, &[4, 3, 1]),
        (46, &[1]),
        (47, &[5]),
        (48, &[5, 3, 2]),
        (49, &[9]),
        (50, &[4, 3, 2]),
        (51, &[6, 3, 1]),
        (52, &[3]),
        (53, &[6, 2, 1]),
        (54, &[9]),
        (55, &[7]),
        (56, &[7, 4, 2]),
        (57, &[4]),
        (58, &[19]),
        (59, &[7, 4, 2]),
        (60, &[1]),
        (61, &[5, 2, 1]),
        (62, &[29]),
        (63, &[1]),
    ];
    t.iter()
        .find(|e| e.0 == n)
        .map(|e| e.1.iter().fold(1u64, |acc, &k| acc | (1u64 << k)))
}

fn deg(a: u64) -> i32 {
    63 - a.leading_zeros() as i32
}

impl GF2n {
    pub fn new(n: u32) -> Self {
        GF2n::with_tail(
            n,
            default_tail(n).expect("no tabulated irreducible for this n"),
        )
    }

    /// GF(2)[x]/(x^n + tail); panics if the polynomial is reducible.
    pub fn with_tail(n: u32, tail: u64) -> Self {
        assert!((2..=63).contains(&n) && tail < (1u64 << n) && tail & 1 == 1);
        #[cfg(target_arch = "x86_64")]
        let hw = std::arch::is_x86_feature_detected!("pclmulqdq");
        let mut f = GF2n {
            n,
            red: tail,
            mask: (1u64 << n) - 1,
            #[cfg(target_arch = "x86_64")]
            hw,
            tmask: 0,
            qsolve: vec![],
            qtab: vec![],
            sqrt_x: 0,
        };
        assert!(f.is_irreducible(), "x^{n} + {tail:#x} is reducible");
        // trace of each basis element
        for i in 0..n {
            let mut a = 1u64 << i;
            let mut s = 0u64;
            for _ in 0..n {
                s ^= a;
                a = f.sq(a);
            }
            debug_assert!(s <= 1);
            f.tmask |= s << i;
        }
        // echelon form of L(z) = z^2 + z over the basis x^i
        let mut rows: Vec<(u64, u64)> = (0..n)
            .map(|i| (f.sq(1u64 << i) ^ (1u64 << i), 1u64 << i))
            .collect();
        let mut piv = vec![];
        for bit in (0..n).rev() {
            if let Some(k) = rows.iter().position(|r| (r.0 >> bit) & 1 == 1) {
                let r = rows.swap_remove(k);
                for o in rows.iter_mut() {
                    if (o.0 >> bit) & 1 == 1 {
                        o.0 ^= r.0;
                        o.1 ^= r.1;
                    }
                }
                piv.push(r);
            }
        }
        f.qsolve = piv;
        let nbytes = (n as usize).div_ceil(8);
        f.qtab = (0..nbytes)
            .map(|k| {
                let mut t = [0u64; 256];
                for (v, e) in t.iter_mut().enumerate() {
                    *e = f.solve_echelon((v as u64) << (8 * k)).0;
                }
                t
            })
            .collect();
        let mut r = 2u64;
        for _ in 0..n - 1 {
            r = f.sq(r);
        }
        f.sqrt_x = r;
        f
    }

    /// Run the echelon elimination: (z, residual). Linear in c.
    fn solve_echelon(&self, c: u64) -> (u64, u64) {
        let (mut c, mut z) = (c, 0u64);
        for &(img, pre) in &self.qsolve {
            let bit = deg(img);
            if (c >> bit) & 1 == 1 {
                c ^= img;
                z ^= pre;
            }
        }
        (z, c)
    }

    #[inline(always)]
    pub fn clmul(&self, a: u64, b: u64) -> u128 {
        #[cfg(target_arch = "x86_64")]
        if self.hw {
            return unsafe { clmul_hw(a, b) };
        }
        clmul_soft(a, b)
    }

    #[inline(always)]
    pub fn reduce(&self, mut t: u128) -> u64 {
        let n = self.n;
        loop {
            let hi = (t >> n) as u64;
            if hi == 0 {
                return t as u64;
            }
            t = (t & self.mask as u128) ^ self.clmul(hi, self.red);
        }
    }

    /// Absolute trace to GF(2).
    pub fn trace(&self, a: u64) -> u64 {
        ((a & self.tmask).count_ones() & 1) as u64
    }

    /// A root z of z^2 + z = c (the other is z + 1), or None when Tr(c) = 1.
    pub fn solve_quadratic(&self, c: u64) -> Option<u64> {
        if self.trace(c) == 1 {
            return None;
        }
        let mut z = 0u64;
        for (k, t) in self.qtab.iter().enumerate() {
            z ^= t[((c >> (8 * k)) & 0xff) as usize];
        }
        Some(z)
    }

    /// Rabin's test: x^(2^n) = x mod f and gcd(x^(2^(n/r)) - x, f) = 1 for prime r | n.
    fn is_irreducible(&self) -> bool {
        let n = self.n;
        let pow2k = |k: u32| {
            let mut a = 2u64; // x
            for _ in 0..k {
                a = self.sq(a);
            }
            a
        };
        if pow2k(n) != 2 {
            return false;
        }
        let full = (1u64 << n) | self.red; // fits: n <= 63
        let gcd = |mut a: u64, mut b: u64| {
            while b != 0 {
                while a != 0 && deg(a) >= deg(b) {
                    a ^= b << (deg(a) - deg(b));
                }
                std::mem::swap(&mut a, &mut b);
            }
            a
        };
        let mut m = n;
        let mut r = 2;
        while m > 1 {
            if m.is_multiple_of(r) {
                if gcd(full, pow2k(n / r) ^ 2) != 1 {
                    return false;
                }
                while m.is_multiple_of(r) {
                    m /= r;
                }
            }
            r += 1;
        }
        true
    }
}

impl Field for GF2n {
    type E = u64;
    fn zero(&self) -> u64 {
        0
    }
    fn one(&self) -> u64 {
        1
    }
    #[inline(always)]
    fn add(&self, a: u64, b: u64) -> u64 {
        a ^ b
    }
    #[inline(always)]
    fn sub(&self, a: u64, b: u64) -> u64 {
        a ^ b
    }
    #[inline(always)]
    fn neg(&self, a: u64) -> u64 {
        a
    }
    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        #[cfg(target_arch = "x86_64")]
        if self.hw {
            return unsafe { mul_hw(a, b, self.n, self.mask, self.red) };
        }
        self.reduce(self.clmul(a, b))
    }
    fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "inverse of zero");
        let (mut u, mut v) = (a, (1u64 << self.n) | self.red);
        let (mut g1, mut g2) = (1u64, 0u64);
        while u != 1 {
            let mut j = deg(u) - deg(v);
            if j < 0 {
                std::mem::swap(&mut u, &mut v);
                std::mem::swap(&mut g1, &mut g2);
                j = -j;
            }
            u ^= v << j;
            g1 ^= g2 << j;
        }
        g1 & self.mask
    }
    fn from_u64(&self, n: u64) -> u64 {
        n & 1
    }
    fn char(&self) -> u64 {
        2
    }
    fn size(&self) -> u128 {
        1u128 << self.n
    }
    fn q(&self) -> Big {
        Big::from_u128(1u128 << self.n)
    }
    fn random(&self, rng: &mut Rng) -> u64 {
        rng.next() & self.mask
    }
    /// sqrt(sum a_i x^i) = sum a_{2i} x^i + sqrt(x) sum a_{2i+1} x^i (squaring is additive).
    fn sqrt(&self, a: u64) -> Option<u64> {
        fn compress(mut v: u64) -> u64 {
            // gather bits 0, 2, 4, ... into the low half
            v &= 0x5555_5555_5555_5555;
            v = (v | (v >> 1)) & 0x3333_3333_3333_3333;
            v = (v | (v >> 2)) & 0x0f0f_0f0f_0f0f_0f0f;
            v = (v | (v >> 4)) & 0x00ff_00ff_00ff_00ff;
            v = (v | (v >> 8)) & 0x0000_ffff_0000_ffff;
            (v | (v >> 16)) & 0x0000_0000_ffff_ffff
        }
        Some(compress(a) ^ self.mul(self.sqrt_x, compress(a >> 1)))
    }
    fn conv_trunc(&self, a: &[u64], b: &[u64], nout: usize) -> Vec<u64> {
        (0..nout)
            .map(|k| {
                let lo = k.saturating_sub(b.len() - 1);
                let hi = k.min(a.len() - 1);
                if lo > hi {
                    return 0;
                }
                #[cfg(target_arch = "x86_64")]
                if self.hw {
                    return unsafe { dot_hw(a, b, lo, hi, k, self.n, self.mask, self.red) };
                }
                let mut acc = 0u128;
                for i in lo..=hi {
                    acc ^= self.clmul(a[i], b[k - i]);
                }
                self.reduce(acc)
            })
            .collect()
    }
    fn conv(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        self.conv_trunc(a, b, a.len() + b.len() - 1)
    }
    fn dot_rev(&self, a: &[u64], b: &[u64]) -> u64 {
        let n = a.len().min(b.len());
        let mut acc = 0u128;
        for i in 0..n {
            acc ^= self.clmul(a[i], b[b.len() - 1 - i]);
        }
        self.reduce(acc)
    }
}
