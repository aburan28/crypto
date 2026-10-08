//! Ternary fields GF(3^n) = F_3[t]/(f), elements as a packed base-3 digit vector in a u64
//! (coefficient of t^i is the i-th ternary digit; n <= 40 so 3^n fits in u64 and the packed form
//! fits too). The packed form is the external encoding only: arithmetic runs bitsliced, a digit
//! vector held as two bit masks (digit 1 sets a bit of `lo`, digit 2 a bit of `hi`), converted
//! eight digits at a time through tables. Addition is a few bitwise operations, multiplication is
//! shift-and-add over the set bits of one operand followed by reduction with the sparse
//! t^n = red(t), cubing is the Frobenius (spread digit i to 3i), and inversion is Itoh–Tsujii
//! (a^(3^n - 2) by an addition chain in Frobenius powers). Square roots when they exist.
//! Characteristic 3: the elliptic-curve work on top uses the general Weierstrass model.
use crate::bigint::Big;
use crate::field::{Field, Rng};
use std::ops::{BitAnd, BitOr, BitXor, Shl, Shr};

/// A digit vector c[0..n] (each 0, 1, 2) packed as sum c[i] 3^i.
fn to_digits(mut v: u64, n: u32) -> Vec<u8> {
    let mut d = vec![0u8; n as usize];
    for x in d.iter_mut() {
        *x = (v % 3) as u8;
        v /= 3;
    }
    d
}

/// 3^8: eight packed digits per table step.
const CHUNK: u64 = 6561;

/// SPLIT[v] for v < 3^8: the eight base-3 digits of v as masks, `lo` in the low byte, `hi` in
/// the high byte.
const fn split_table() -> [u16; CHUNK as usize] {
    let mut t = [0u16; CHUNK as usize];
    let mut v = 0;
    while v < CHUNK as usize {
        let (mut x, mut i, mut lo, mut hi) = (v, 0, 0u16, 0u16);
        while i < 8 {
            match x % 3 {
                1 => lo |= 1 << i,
                2 => hi |= 1 << i,
                _ => {}
            }
            x /= 3;
            i += 1;
        }
        t[v] = lo | (hi << 8);
        v += 1;
    }
    t
}
/// JOIN[m] = sum of 3^i over the set bits i of the byte m.
const fn join_table() -> [u16; 256] {
    let mut t = [0u16; 256];
    let mut m = 0;
    while m < 256 {
        let (mut i, mut p, mut s) = (0, 1u16, 0u16);
        while i < 8 {
            if m >> i & 1 == 1 {
                s += p;
            }
            p *= 3;
            i += 1;
        }
        t[m] = s;
        m += 1;
    }
    t
}
/// SPREAD3[m]: bit i of the byte m moved to bit 3i (the Frobenius on digit positions).
const fn spread3_table() -> [u32; 256] {
    let mut t = [0u32; 256];
    let mut m = 0;
    while m < 256 {
        let (mut i, mut s) = (0, 0u32);
        while i < 8 {
            if m >> i & 1 == 1 {
                s |= 1 << (3 * i);
            }
            i += 1;
        }
        t[m] = s;
        m += 1;
    }
    t
}
static SPLIT: [u16; CHUNK as usize] = split_table();
static JOIN: [u16; 256] = join_table();
static SPREAD3: [u32; 256] = spread3_table();

/// A word holding a bitsliced digit vector (u64 when the product fits, else u128).
trait Lane:
    Copy
    + Eq
    + BitOr<Output = Self>
    + BitXor<Output = Self>
    + BitAnd<Output = Self>
    + Shl<u32, Output = Self>
    + Shr<u32, Output = Self>
{
    const ZERO: Self;
    fn from64(x: u64) -> Self;
    fn low64(self) -> u64;
}
impl Lane for u64 {
    const ZERO: Self = 0;
    #[inline(always)]
    fn from64(x: u64) -> Self {
        x
    }
    #[inline(always)]
    fn low64(self) -> u64 {
        self
    }
}
impl Lane for u128 {
    const ZERO: Self = 0;
    #[inline(always)]
    fn from64(x: u64) -> Self {
        x as u128
    }
    #[inline(always)]
    fn low64(self) -> u64 {
        self as u64
    }
}

/// Digitwise sum mod 3 of two bitsliced digit vectors (lo = digit 1, hi = digit 2).
#[inline(always)]
fn tadd<T: Lane>(a: (T, T), b: (T, T)) -> (T, T) {
    let t = (a.0 | b.1) ^ (a.1 | b.0);
    ((a.1 | b.1) ^ t, (a.0 | b.0) ^ t)
}

/// Product of two bitsliced polynomials of degree < 64 (no reduction): shift-and-add over the
/// set bits of `a`, a digit 2 adding the negated (mask-swapped) shift of `b`.
#[inline(always)]
fn shift_add<T: Lane>(a: (u64, u64), b: (u64, u64)) -> (T, T) {
    let (b1, b2) = (T::from64(b.0), T::from64(b.1));
    let mut acc = (T::ZERO, T::ZERO);
    let mut m = a.0;
    while m != 0 {
        let i = m.trailing_zeros();
        acc = tadd(acc, (b1 << i, b2 << i));
        m &= m - 1;
    }
    let mut m = a.1;
    while m != 0 {
        let i = m.trailing_zeros();
        acc = tadd(acc, (b2 << i, b1 << i));
        m &= m - 1;
    }
    acc
}

/// Bit i -> bit 3i for the low 40 bits of `m` (result below 2^120).
#[inline(always)]
fn spread3(m: u64) -> u128 {
    let mut s = 0u128;
    let mut j = 0;
    let mut m = m;
    while m != 0 {
        s |= (SPREAD3[(m & 0xff) as usize] as u128) << (24 * j);
        m >>= 8;
        j += 1;
    }
    s
}

#[derive(Clone, Debug)]
pub struct GF3n {
    pub n: u32,
    /// reduction polynomial f = t^n - red(t); red as digit vector of degree < n
    red: Vec<u8>,
    /// the non-zero terms (k, digit) of red, for reduction
    terms: Vec<(u32, u8)>,
    pow3n: u64,
    mask: u64,
    /// number of eight-digit chunks in the packed form
    chunks: u32,
}

impl GF3n {
    /// GF(3^n) with a low-weight monic irreducible t^n - red found by search (Rabin's test):
    /// the fewest non-zero tail terms, and among those the smallest packed tail.
    pub fn new(n: u32) -> Self {
        assert!((1..=40).contains(&n));
        let pow3n = 3u64.pow(n);
        if n == 1 {
            return Self::with_red(n, vec![0], pow3n);
        }
        for weight in 1..=n as usize {
            let mut tails = vec![];
            exact_weight_tails(n as usize, weight, &mut vec![], &mut tails);
            tails.sort_unstable();
            for packed in tails {
                let cand = Self::with_red(n, to_digits(packed, n), pow3n);
                if cand.is_irreducible() {
                    return cand;
                }
            }
        }
        unreachable!("an irreducible of degree n over F_3 exists")
    }

    fn with_red(n: u32, red: Vec<u8>, pow3n: u64) -> Self {
        let terms = red
            .iter()
            .enumerate()
            .filter(|(_, &r)| r != 0)
            .map(|(k, &r)| (k as u32, r))
            .collect();
        GF3n {
            n,
            red,
            terms,
            pow3n,
            mask: (1u64 << n) - 1,
            chunks: n.div_ceil(8),
        }
    }

    /// packed base-3 -> bitsliced (lo, hi)
    #[inline(always)]
    fn unpack(&self, mut v: u64) -> (u64, u64) {
        let (mut lo, mut hi, mut s) = (0u64, 0u64, 0);
        while v != 0 {
            let c = SPLIT[(v % CHUNK) as usize] as u64;
            lo |= (c & 0xff) << s;
            hi |= (c >> 8) << s;
            v /= CHUNK;
            s += 8;
        }
        (lo, hi)
    }
    /// bitsliced (lo, hi) -> packed base-3
    #[inline(always)]
    fn pack(&self, a: (u64, u64)) -> u64 {
        let mut v = 0u64;
        for c in (0..self.chunks).rev() {
            let s = 8 * c;
            v = v * CHUNK
                + JOIN[(a.0 >> s & 0xff) as usize] as u64
                + 2 * JOIN[(a.1 >> s & 0xff) as usize] as u64;
        }
        v
    }

    /// Reduce a bitsliced polynomial (of any degree the word holds) modulo t^n - red.
    #[inline(always)]
    fn reduce<T: Lane>(&self, mut c: (T, T)) -> (u64, u64) {
        let n = self.n;
        let mask = T::from64(self.mask);
        loop {
            let h = (c.0 >> n, c.1 >> n);
            if (h.0 | h.1) == T::ZERO {
                return (c.0.low64(), c.1.low64());
            }
            c = (c.0 & mask, c.1 & mask);
            for &(k, r) in &self.terms {
                // h t^n = h red: add r h t^k (r = 2 swaps the masks)
                let s = if r == 1 {
                    (h.0 << k, h.1 << k)
                } else {
                    (h.1 << k, h.0 << k)
                };
                c = tadd(c, s);
            }
        }
    }

    #[inline(always)]
    fn mul_bits(&self, a: (u64, u64), b: (u64, u64)) -> (u64, u64) {
        if self.n <= 32 {
            self.reduce(shift_add::<u64>(a, b)) // degree <= 62
        } else {
            self.reduce(shift_add::<u128>(a, b))
        }
    }

    /// a^3, the Frobenius.
    #[inline(always)]
    fn cube_bits(&self, a: (u64, u64)) -> (u64, u64) {
        let c = (spread3(a.0), spread3(a.1));
        if self.n <= 22 {
            self.reduce((c.0 as u64, c.1 as u64)) // degree <= 63
        } else {
            self.reduce(c)
        }
    }

    fn pow_bits(&self, a: (u64, u64), mut e: u64) -> (u64, u64) {
        let mut r = (1u64, 0u64);
        let mut b = a;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul_bits(r, b);
            }
            e >>= 1;
            if e > 0 {
                b = self.mul_bits(b, b);
            }
        }
        r
    }

    fn pow(&self, a: u64, e: u64) -> u64 {
        self.pack(self.pow_bits(self.unpack(a), e))
    }

    /// Rabin irreducibility test for the monic degree-n polynomial t^n - red over F_3.
    fn is_irreducible(&self) -> bool {
        // x^(3^n) = x (mod f), and for each maximal proper divisor n/q: gcd(x^(3^(n/q)) - x, f) = 1
        let n = self.n;
        let x = self.unpack(self.coeff_x());
        let xq = |e: u32| (0..e).fold(x, |v, _| self.cube_bits(v)); // x^(3^e) mod f
        if xq(n) != x {
            return false;
        }
        for p in distinct_primes(n) {
            let d = self.pack(tadd(xq(n / p), (x.1, x.0)));
            if !self.poly_gcd_is_unit(d) {
                return false;
            }
        }
        true
    }

    fn coeff_x(&self) -> u64 {
        if self.n == 1 {
            0 // in F_3 = F_3[t]/(t), t is 0
        } else {
            3 // digit vector [0,1] = t
        }
    }

    /// Whether gcd(g, f) is a unit, g a reduced element (degree < n) given as a true polynomial
    /// over F_3; f may be reducible here, so this is a real gcd.
    fn poly_gcd_is_unit(&self, g: u64) -> bool {
        if g == 0 {
            return false;
        }
        let mut a = self.full_f();
        let mut b = to_digits(g, self.n);
        trim3(&mut b);
        while !b.is_empty() {
            let r = polyrem3(&a, &b);
            a = b;
            b = r;
        }
        a.len() == 1 // gcd degree 0
    }
    fn full_f(&self) -> Vec<u8> {
        let mut f = vec![0u8; self.n as usize + 1];
        f[self.n as usize] = 1;
        for (k, &r) in self.red.iter().enumerate() {
            f[k] = (3 - r) % 3; // f = t^n - red
        }
        trim3(&mut f);
        f
    }

    /// a^-1 = a^(3^n - 2) = a (b^3)^2 with b = a^(1 + 3 + ... + 3^(n-2)) (Itoh–Tsujii): the
    /// exponents beta_k = a^((3^k - 1)/2) satisfy beta_(j+k) = beta_j^(3^k) beta_k, so beta_(n-1)
    /// costs O(log n) multiplications and about n Frobenius cubings.
    fn inv_bits(&self, a: (u64, u64)) -> (u64, u64) {
        let m = self.n - 1;
        let mut beta = (1u64, 0u64); // beta_0
        if m > 0 {
            beta = a;
            let mut k = 1u32;
            for bit in (0..31 - m.leading_zeros()).rev() {
                let mut f = beta;
                for _ in 0..k {
                    f = self.cube_bits(f);
                }
                beta = self.mul_bits(f, beta);
                k *= 2;
                if m >> bit & 1 == 1 {
                    beta = self.mul_bits(self.cube_bits(beta), a);
                    k += 1;
                }
            }
        }
        let s = self.cube_bits(beta);
        self.mul_bits(a, self.mul_bits(s, s))
    }
}

/// All tails of exactly `w` non-zero digits among positions 0..n (position 0 non-zero, or t
/// divides f), as packed values.
fn exact_weight_tails(n: usize, w: usize, cur: &mut Vec<usize>, out: &mut Vec<u64>) {
    if cur.len() == w {
        if cur[0] != 0 {
            return;
        }
        // every 1/2 assignment of the chosen positions
        for signs in 0..1u64 << w {
            let mut v = 0u64;
            for (i, &pos) in cur.iter().enumerate() {
                v += (1 + (signs >> i & 1)) * 3u64.pow(pos as u32);
            }
            out.push(v);
        }
        return;
    }
    let start = cur.last().map_or(0, |&x| x + 1);
    for pos in start..n {
        if n - pos < w - cur.len() {
            break;
        }
        cur.push(pos);
        exact_weight_tails(n, w, cur, out);
        cur.pop();
    }
}

fn trim3(v: &mut Vec<u8>) {
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    if v.len() == 1 && v[0] == 0 {
        v.clear();
    }
}

/// remainder of a mod b over F_3 (b monic-izable).
fn polyrem3(a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut r = a.to_vec();
    let bl = b.len();
    let binv = inv3(b[bl - 1]);
    while r.len() >= bl && !(r.len() == 1 && r[0] == 0) {
        trim3(&mut r);
        if r.len() < bl {
            break;
        }
        let shift = r.len() - bl;
        let c = (r[r.len() - 1] * binv) % 3;
        if c == 0 {
            r.pop();
            continue;
        }
        for i in 0..bl {
            r[shift + i] = (r[shift + i] + 3 - (c * b[i]) % 3) % 3;
        }
        trim3(&mut r);
    }
    trim3(&mut r);
    r
}
fn inv3(x: u8) -> u8 {
    match x {
        1 => 1,
        2 => 2,
        _ => 0,
    }
}
fn distinct_primes(mut n: u32) -> Vec<u32> {
    let mut ps = vec![];
    let mut d = 2;
    while d * d <= n {
        if n.is_multiple_of(d) {
            ps.push(d);
            while n.is_multiple_of(d) {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        ps.push(n);
    }
    ps
}

impl Field for GF3n {
    type E = u64;
    fn zero(&self) -> u64 {
        0
    }
    fn one(&self) -> u64 {
        1
    }
    fn add(&self, a: u64, b: u64) -> u64 {
        self.pack(tadd(self.unpack(a), self.unpack(b)))
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        let b = self.unpack(b);
        self.pack(tadd(self.unpack(a), (b.1, b.0)))
    }
    fn neg(&self, a: u64) -> u64 {
        let a = self.unpack(a);
        self.pack((a.1, a.0))
    }
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.pack(self.mul_bits(self.unpack(a), self.unpack(b)))
    }
    fn sq(&self, a: u64) -> u64 {
        let a = self.unpack(a);
        self.pack(self.mul_bits(a, a))
    }
    fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "inverse of zero in GF(3^n)");
        self.pack(self.inv_bits(self.unpack(a)))
    }
    fn from_u64(&self, n: u64) -> u64 {
        n % 3
    }
    fn char(&self) -> u64 {
        3
    }
    fn size(&self) -> u128 {
        self.pow3n as u128
    }
    fn q(&self) -> Big {
        Big::from_u64(self.pow3n)
    }
    fn random(&self, rng: &mut Rng) -> u64 {
        rng.below(self.pow3n)
    }
    fn sqrt(&self, a: u64) -> Option<u64> {
        if a == 0 {
            return Some(0);
        }
        // q = 3^n is odd; a is a square iff a^((q-1)/2) = 1, then r = a^((q+1)/4) if q = 3 mod 4,
        // else Tonelli-Shanks. q = 3^n: 3^n mod 4 is 3 if n odd, 1 if n even.
        if self.pow(a, (self.pow3n - 1) / 2) != 1 {
            return None;
        }
        if self.pow3n % 4 == 3 {
            let r = self.pow(a, (self.pow3n + 1) / 4);
            return Some(r);
        }
        // Tonelli-Shanks over F_q
        let qm1 = self.pow3n - 1;
        let mut s = 0u32;
        let mut m = qm1;
        while m.is_multiple_of(2) {
            m /= 2;
            s += 1;
        }
        let mut rng = Rng::new(0x3333);
        let z = loop {
            let z = 1 + rng.below(self.pow3n - 1);
            if self.pow(z, qm1 / 2) != 1 {
                break z;
            }
        };
        let sq = |x: (u64, u64)| self.mul_bits(x, x);
        let mut c = self.unpack(self.pow(z, m));
        let mut t = self.unpack(self.pow(a, m));
        let mut r = self.unpack(self.pow(a, m.div_ceil(2)));
        let one = (1u64, 0u64);
        let mut mm = s;
        while t != one {
            let mut i = 0u32;
            let mut tt = t;
            while tt != one {
                tt = sq(tt);
                i += 1;
            }
            let mut b = c;
            for _ in 0..(mm - i - 1) {
                b = sq(b);
            }
            mm = i;
            c = sq(b);
            t = self.mul_bits(t, c);
            r = self.mul_bits(r, b);
        }
        Some(self.pack(r))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The bitsliced digit addition and the table conversions agree with digit arithmetic.
    #[test]
    fn bitsliced_digits() {
        let enc = |d: u64| ((d == 1) as u64, (d == 2) as u64);
        for x in 0..3 {
            for y in 0..3 {
                assert_eq!(tadd(enc(x), enc(y)), enc((x + y) % 3));
            }
        }
        let f = GF3n::with_red(40, vec![0; 40], 3u64.pow(40));
        let mut rng = Rng::new(9);
        for _ in 0..1000 {
            let v = rng.below(f.pow3n);
            let (lo, hi) = f.unpack(v);
            let d = to_digits(v, 40);
            for (i, &di) in d.iter().enumerate() {
                assert_eq!(
                    (lo >> i & 1, hi >> i & 1),
                    ((di == 1) as u64, (di == 2) as u64)
                );
            }
            assert_eq!(f.pack((lo, hi)), v);
        }
    }

    /// The constructor picks the same reduction polynomials as the earlier exhaustive scan
    /// (which visited every packed tail in increasing order, by weight; recorded for n <= 17,
    /// where that scan still finished).
    #[test]
    fn moduli_match_exhaustive_scan() {
        let expect: [&[u8]; 17] = [
            &[0],
            &[2, 0],
            &[1, 1, 0],
            &[1, 1, 0, 0],
            &[1, 1, 0, 0, 0],
            &[1, 1, 0, 0, 0, 0],
            &[2, 0, 1, 0, 0, 0, 0],
            &[1, 0, 1, 0, 0, 0, 0, 0],
            &[2, 0, 0, 0, 1, 0, 0, 0, 0],
            &[2, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            &[2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            &[1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            &[1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            &[1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            &[2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            &[1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            &[1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ];
        for (i, red) in expect.iter().enumerate() {
            assert_eq!(GF3n::new(i as u32 + 1).red, *red, "n = {}", i + 1);
        }
    }

    /// Bitsliced multiplication and Itoh–Tsujii inversion agree with schoolbook digit
    /// arithmetic (convolution mod 3, then reduction by t^n = red, digit by digit).
    #[test]
    fn mul_inv_match_digit_schoolbook() {
        for n in [2u32, 7, 16, 25, 33, 40] {
            let f = GF3n::new(n);
            let nn = n as usize;
            let school = |a: u64, b: u64| -> u64 {
                let (da, db) = (to_digits(a, n), to_digits(b, n));
                let mut c = vec![0u8; 2 * nn];
                for i in 0..nn {
                    for j in 0..nn {
                        c[i + j] = (c[i + j] + da[i] * db[j]) % 3;
                    }
                }
                for i in (nn..2 * nn).rev() {
                    let ci = std::mem::take(&mut c[i]);
                    for (k, &r) in f.red.iter().enumerate() {
                        c[i - nn + k] = (c[i - nn + k] + ci * r) % 3;
                    }
                }
                c[..nn].iter().rev().fold(0u64, |v, &x| v * 3 + x as u64)
            };
            let mut rng = Rng::new(40 + n as u64);
            for _ in 0..300 {
                let (a, b) = (f.random(&mut rng), f.random(&mut rng));
                assert_eq!(f.mul(a, b), school(a, b), "n = {n}");
                if a != 0 {
                    assert_eq!(school(a, f.inv(a)), 1, "n = {n}");
                }
            }
        }
    }
}
