//! # `Fp52`: `F_p` for `p < 2⁵²` with an AVX-512 IFMA row kernel
//!
//! A [`PrimeField`] backend for the F4 engines' dense row updates at primes
//! between `2³²` and `2⁵²` — the Kummer primes `c·2^k + 1` with
//! `k ∈ {40, 48}` of the PKM tower study among them. Scalar arithmetic is
//! [`Fp64`]'s, element for element: an element is the same `u64`
//! Montgomery representative (`a·2⁶⁴ mod p`), so `Fp52` and `Fp64` agree
//! on every value and every report. Only the lazy accumulator differs.
//!
//! **The accumulator.** `Wide = [u64; 2] = (lo, hi)` holds the integer
//! `V = lo + hi·2⁵²`, and the element whose representative is
//! `V·2⁻⁶⁴ mod p` (so a product of two representatives, `a·b·2¹²⁸`, adds
//! the element `a·b`, as in `Fp64`). A product `x·y < 2¹⁰⁴` of two
//! representatives is added as its low 52 bits to `lo` and its high 52
//! bits to `hi` — exactly what one `vpmadd52luq` and one `vpmadd52huq`
//! do, eight 64-bit lanes at a time, so the row kernel
//! `acc[i] += c·row[i]` processes four entries per pair of instructions.
//! A load stores `lo = a·2⁶⁴ mod p < 2⁵²`, `hi = 0`.
//!
//! **Capacity.** Each product adds less than `2⁵²` to each word, so after
//! a load and `k` products `lo < (k + 1)·2⁵²` and `hi < k·2⁵²`: both stay
//! below `2⁶⁴` for `k ≤ 4095` ([`PrimeField::lazy_capacity`]). The
//! engines count products per row and reload before that
//! ([`f4_tower_generic`](super::f4_tower_generic)).
//!
//! **Reduction.** With `hi = h₁·2¹² + h₀` (`h₀ < 2¹²`),
//! `V = L + h₁·2⁶⁴` where `L = lo + h₀·2⁵² < 2⁶⁵`, so
//! `V·2⁻⁶⁴ ≡ redc(L) + h₁ (mod p)`: one Montgomery reduction (valid for
//! `L < p·2⁶⁴`) and `h₁ mod p`.
//!
//! Without AVX-512 IFMA the same accumulator runs as a scalar loop
//! (correct, slower than `Fp64`; prefer `Fp64` there).

use num_bigint::BigUint;

use super::fp_backend::{Fp64, PrimeField, SqrtPre};

/// `F_p` for an odd prime `2 < p < 2⁵²`; see the module documentation.
#[derive(Clone, Debug)]
pub struct Fp52 {
    inner: Fp64,
    p: u64,
    /// `⌊2⁶⁴ / p⌋`, for `h₁ mod p` (Barrett, `h₁ < 2⁵²`).
    barrett: u64,
    ifma: bool,
}

const MASK52: u64 = (1 << 52) - 1;

impl Fp52 {
    /// The field for an odd prime `3 ≤ p < 2⁵²` (primality is the caller's
    /// promise). Uses the IFMA kernel when this CPU has AVX-512 IFMA.
    pub fn new(p: u64) -> Self {
        assert!(
            p & 1 == 1 && (3..1u64 << 52).contains(&p),
            "Fp52 needs an odd prime below 2^52"
        );
        Fp52 {
            inner: Fp64::new(p),
            p,
            barrett: ((1u128 << 64) / u128::from(p)) as u64,
            ifma: has_ifma(),
        }
    }

    /// The scalar accumulator loop only (for tests and comparisons).
    pub fn new_scalar(p: u64) -> Self {
        let mut f = Self::new(p);
        f.ifma = false;
        f
    }

    /// Whether the row kernel uses AVX-512 IFMA.
    pub fn uses_ifma(&self) -> bool {
        self.ifma
    }

    pub fn p(&self) -> u64 {
        self.p
    }

    /// `x mod p` for `x < 2⁶⁴`.
    #[inline(always)]
    fn reduce_u64(&self, x: u64) -> u64 {
        let q = ((u128::from(x) * u128::from(self.barrett)) >> 64) as u64;
        let r = x - q * self.p;
        if r >= self.p {
            r - self.p
        } else {
            r
        }
    }
}

fn has_ifma() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("avx512f") && std::is_x86_feature_detected!("avx512ifma")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

#[inline(always)]
fn axpy_scalar(acc: &mut [[u64; 2]], c: u64, row: &[u64]) {
    for (w, &y) in acc.iter_mut().zip(row) {
        let prod = u128::from(c) * u128::from(y);
        w[0] += (prod as u64) & MASK52;
        w[1] += (prod >> 52) as u64;
    }
}

/// `acc[i] += c·row[i]` with IFMA, four entries per 512-bit vector: the
/// row's entries are duplicated into lane pairs, the low half of each
/// product goes to the even (`lo`) lanes and the high half to the odd
/// (`hi`) lanes through merge masks.
///
/// # Safety
/// The CPU must support AVX-512F and AVX-512 IFMA.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512ifma")]
unsafe fn axpy_ifma(acc: &mut [[u64; 2]], c: u64, row: &[u64]) {
    use std::arch::x86_64::*;
    let n = acc.len().min(row.len());
    let cv = _mm512_set1_epi64(c as i64);
    let dup = _mm512_set_epi64(3, 3, 2, 2, 1, 1, 0, 0);
    let ap = acc.as_mut_ptr() as *mut u64;
    let rp = row.as_ptr();
    let mut i = 0usize;
    while i + 8 <= n {
        // SAFETY: i + 8 ≤ n bounds both slices; unaligned loads.
        let y0 = _mm256_loadu_si256(rp.add(i) as *const __m256i);
        let y1 = _mm256_loadu_si256(rp.add(i + 4) as *const __m256i);
        let y0 = _mm512_permutexvar_epi64(dup, _mm512_castsi256_si512(y0));
        let y1 = _mm512_permutexvar_epi64(dup, _mm512_castsi256_si512(y1));
        let a0 = _mm512_loadu_si512(ap.add(2 * i) as *const __m512i);
        let a1 = _mm512_loadu_si512(ap.add(2 * i + 8) as *const __m512i);
        let a0 = _mm512_mask_madd52lo_epu64(a0, 0x55, cv, y0);
        let a1 = _mm512_mask_madd52lo_epu64(a1, 0x55, cv, y1);
        let a0 = _mm512_mask_madd52hi_epu64(a0, 0xAA, cv, y0);
        let a1 = _mm512_mask_madd52hi_epu64(a1, 0xAA, cv, y1);
        _mm512_storeu_si512(ap.add(2 * i) as *mut __m512i, a0);
        _mm512_storeu_si512(ap.add(2 * i + 8) as *mut __m512i, a1);
        i += 8;
    }
    while i + 4 <= n {
        let y0 = _mm256_loadu_si256(rp.add(i) as *const __m256i);
        let y0 = _mm512_permutexvar_epi64(dup, _mm512_castsi256_si512(y0));
        let a0 = _mm512_loadu_si512(ap.add(2 * i) as *const __m512i);
        let a0 = _mm512_mask_madd52lo_epu64(a0, 0x55, cv, y0);
        let a0 = _mm512_mask_madd52hi_epu64(a0, 0xAA, cv, y0);
        _mm512_storeu_si512(ap.add(2 * i) as *mut __m512i, a0);
        i += 4;
    }
    axpy_scalar(&mut acc[i..n], c, &row[i..n]);
}

impl PrimeField for Fp52 {
    type Elem = u64;
    type Wide = [u64; 2];

    fn name(&self) -> &'static str {
        if self.ifma {
            "Fp52-ifma"
        } else {
            "Fp52-scalar"
        }
    }
    fn modulus(&self) -> &BigUint {
        self.inner.modulus()
    }
    #[inline(always)]
    fn zero(&self) -> u64 {
        0
    }
    #[inline(always)]
    fn one(&self) -> u64 {
        self.inner.one()
    }
    #[inline(always)]
    fn from_u64(&self, x: u64) -> u64 {
        self.inner.from_u64(x)
    }
    fn from_biguint(&self, x: &BigUint) -> u64 {
        self.inner.from_biguint(x)
    }
    fn to_biguint(&self, a: u64) -> BigUint {
        self.inner.to_biguint(a)
    }
    #[inline(always)]
    fn to_u64(&self, a: u64) -> Option<u64> {
        self.inner.to_u64(a)
    }
    #[inline(always)]
    fn add(&self, a: u64, b: u64) -> u64 {
        self.inner.add(a, b)
    }
    #[inline(always)]
    fn sub(&self, a: u64, b: u64) -> u64 {
        self.inner.sub(a, b)
    }
    #[inline(always)]
    fn neg(&self, a: u64) -> u64 {
        self.inner.neg(a)
    }
    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.inner.mul(a, b)
    }
    fn inv(&self, a: u64) -> Option<u64> {
        self.inner.inv(a)
    }
    fn sqrt_pre(&self) -> &SqrtPre<u64> {
        self.inner.sqrt_pre()
    }
    fn legendre(&self, a: u64) -> i32 {
        self.inner.legendre(a)
    }
    fn sqrt(&self, a: u64) -> Option<u64> {
        self.inner.sqrt(a)
    }

    fn lazy_capacity(&self) -> u64 {
        4095
    }
    #[inline(always)]
    fn wide_zero(&self) -> [u64; 2] {
        [0, 0]
    }
    #[inline(always)]
    fn wide_from(&self, a: u64) -> [u64; 2] {
        // V = a·2⁶⁴ mod p, so V·2⁻⁶⁴ ≡ a.
        [self.inner.to_mont(a), 0]
    }
    #[inline(always)]
    fn wide_mul_add(&self, w: &mut [u64; 2], a: u64, b: u64) {
        let prod = u128::from(a) * u128::from(b);
        debug_assert!(w[0] < u64::MAX - MASK52 && w[1] < u64::MAX - MASK52);
        w[0] += (prod as u64) & MASK52;
        w[1] += (prod >> 52) as u64;
    }
    #[inline(always)]
    fn wide_reduce(&self, w: &[u64; 2]) -> u64 {
        let (lo, hi) = (w[0], w[1]);
        let l = u128::from(lo) + (u128::from(hi & 0xfff) << 52);
        let h1 = hi >> 12;
        self.inner.add(self.inner.redc(l), self.reduce_u64(h1))
    }
    #[inline]
    fn row_axpy(&self, acc: &mut [[u64; 2]], c: u64, row: &[u64]) {
        #[cfg(target_arch = "x86_64")]
        if self.ifma {
            // SAFETY: `has_ifma` found AVX-512F and AVX-512 IFMA.
            unsafe { axpy_ifma(acc, c, row) };
            return;
        }
        axpy_scalar(acc, c, row);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    /// Deterministic Miller–Rabin below `2⁶⁴`.
    fn is_prime(n: u64) -> bool {
        const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
        if n < 2 {
            return false;
        }
        for q in BASES {
            if n % q == 0 {
                return n == q;
            }
        }
        let m = |a: u64, b: u64| ((u128::from(a) * u128::from(b)) % u128::from(n)) as u64;
        let pw = |mut a: u64, mut e: u64| {
            let mut r = 1;
            while e > 0 {
                if e & 1 == 1 {
                    r = m(r, a);
                }
                a = m(a, a);
                e >>= 1;
            }
            r
        };
        let s = (n - 1).trailing_zeros();
        let d = (n - 1) >> s;
        'b: for a in BASES {
            let mut x = pw(a, d);
            if x == 1 || x == n - 1 {
                continue;
            }
            for _ in 1..s {
                x = m(x, x);
                if x == n - 1 {
                    continue 'b;
                }
            }
            return false;
        }
        true
    }

    /// Small, 31-bit, 32-bit, the `k = 40, 48` Kummer primes and the
    /// largest prime below `2⁵²`.
    fn primes() -> Vec<u64> {
        let proth = |k: u32| {
            (1u64..)
                .step_by(2)
                .map(|c| (c << k) + 1)
                .find(|&p| is_prime(p))
                .unwrap()
        };
        let top = (1u64 << 52) - 1;
        let below = (0..1000u64)
            .map(|i| top - 2 * i)
            .find(|&p| is_prime(p))
            .unwrap();
        vec![
            3,
            786_433,
            2_013_265_921,
            3_221_225_473,
            proth(40),
            proth(48),
            below,
        ]
    }

    /// Every operation against `BigUint`, and the row kernels (both
    /// dispatches) against exact arithmetic up to the full capacity.
    #[test]
    fn fp52_matches_biguint_and_its_capacity_holds() {
        let mut rng = StdRng::seed_from_u64(52);
        for p in primes() {
            let pb = BigUint::from(p);
            for f in [Fp52::new(p), Fp52::new_scalar(p)] {
                let edge = [0, 1, 2, p - 1, p - 2, p / 2, (p + 1) / 2];
                let vals: Vec<u64> = edge
                    .iter()
                    .copied()
                    .chain((0..200).map(|_| rng.gen_range(0..p)))
                    .collect();
                for &a in &vals {
                    for &b in vals.iter().take(20) {
                        let (ea, eb) = (f.from_u64(a), f.from_u64(b));
                        let want = (BigUint::from(a) * BigUint::from(b)) % &pb;
                        assert_eq!(f.to_biguint(f.mul(ea, eb)), want);
                        let mut w = f.wide_from(ea);
                        f.wide_mul_add(&mut w, ea, eb);
                        let want = (BigUint::from(a) + BigUint::from(a) * BigUint::from(b)) % &pb;
                        assert_eq!(f.to_biguint(f.wide_reduce(&w)), want, "p {p}");
                    }
                }
                // Rows: a load and `cap` products, every entry worst-case
                // (p − 1) and random, against exact sums.
                for trial in 0..4 {
                    let n = rng.gen_range(1..40usize);
                    let cap = f.lazy_capacity() as usize;
                    let k = if trial == 0 {
                        cap
                    } else {
                        rng.gen_range(1..300)
                    };
                    let worst = trial % 2 == 0;
                    let draw = |rng: &mut StdRng| if worst { p - 1 } else { rng.gen_range(0..p) };
                    let init: Vec<u64> = (0..n).map(|_| draw(&mut rng)).collect();
                    let mut exact: Vec<BigUint> = init.iter().map(|&x| BigUint::from(x)).collect();
                    let mut acc: Vec<[u64; 2]> =
                        init.iter().map(|&x| f.wide_from(f.from_u64(x))).collect();
                    for _ in 0..k {
                        let c = draw(&mut rng);
                        let row: Vec<u64> = (0..n).map(|_| draw(&mut rng)).collect();
                        let er: Vec<u64> = row.iter().map(|&x| f.from_u64(x)).collect();
                        f.row_sub(&mut acc, f.from_u64(c), &er);
                        for (e, &y) in exact.iter_mut().zip(&row) {
                            *e = (&*e + &pb * &pb - BigUint::from(c) * BigUint::from(y)) % &pb;
                        }
                    }
                    let mut out = vec![0u64; n];
                    f.row_reduce(&acc, &mut out);
                    let got: Vec<BigUint> = out.iter().map(|&x| f.to_biguint(x)).collect();
                    assert_eq!(got, exact, "p {p} {}", f.name());
                }
            }
        }
    }
}
