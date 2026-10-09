//! Independent exact-size Montgomery field and binary extended-GCD inversion.
//! Arithmetic for public curve data; based on the separately validated walker field.
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// An element in Montgomery form (`x·R mod p`, `R = 2^256`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Element<const N: usize>(pub [u64; N]);

/// The field `GF(p)`.
#[derive(Clone, Debug)]
pub struct PrimeField<const N: usize> {
    p: [u64; N],
    /// `-p^{-1} mod 2^64`.
    n0: u64,
    /// `R^2 mod p`.
    r2: Element<N>,
    one: Element<N>,
    p_big: BigUint,
    /// `(p - 1) / 2`, for the Legendre symbol.
    half: BigUint,
    /// Tonelli–Shanks: `p - 1 = 2^s · q`, `q` odd, and a non-residue's
    /// `q`-th power.
    ts_s: u32,
    ts_q: BigUint,
    ts_z: Element<N>,
    /// The least quadratic non-residue, as an integer.
    nonresidue: u64,
}

fn to_limbs<const N: usize>(v: &BigUint) -> [u64; N] {
    let mut out = [0u64; N];
    for (i, d) in v.to_u64_digits().iter().enumerate().take(N) {
        out[i] = *d;
    }
    out
}

fn from_limbs<const N: usize>(v: &[u64; N]) -> BigUint {
    let mut bytes = Vec::with_capacity(8 * N);
    for limb in v {
        bytes.extend_from_slice(&limb.to_le_bytes());
    }
    BigUint::from_bytes_le(&bytes)
}

fn geq<const N: usize>(a: &[u64; N], b: &[u64; N]) -> bool {
    for i in (0..N).rev() {
        if a[i] != b[i] {
            return a[i] > b[i];
        }
    }
    true
}

fn sub_limbs<const N: usize>(a: &[u64; N], b: &[u64; N]) -> ([u64; N], bool) {
    let mut out = [0u64; N];
    let mut borrow = false;
    for i in 0..N {
        let (d, b1) = a[i].overflowing_sub(b[i]);
        let (d, b2) = d.overflowing_sub(borrow as u64);
        out[i] = d;
        borrow = b1 || b2;
    }
    (out, borrow)
}

fn add_limbs<const N: usize>(a: &[u64; N], b: &[u64; N]) -> ([u64; N], bool) {
    let mut out = [0u64; N];
    let mut carry = false;
    for i in 0..N {
        let (s, c1) = a[i].overflowing_add(b[i]);
        let (s, c2) = s.overflowing_add(carry as u64);
        out[i] = s;
        carry = c1 || c2;
    }
    (out, carry)
}

fn shr1<const N: usize>(a: &mut [u64; N], carry: bool) {
    let mut high = carry as u64;
    for i in (0..N).rev() {
        let next = a[i] & 1;
        a[i] = (a[i] >> 1) | (high << 63);
        high = next;
    }
}

impl<const N: usize> PrimeField<N> {
    /// The field of the odd prime `p < 2^256`.  Primality is the caller's
    /// claim; [`PrimeField<N>::new`] checks it with Miller–Rabin and refuses a
    /// composite.
    pub fn new(p: &BigUint) -> Option<Self> {
        if N == 0
            || N > 16
            || p.bits() > (64 * N) as u64
            || p < &BigUint::from(5u8)
            || !p.bit(0)
            || !is_probable_prime(p)
        {
            return None;
        }
        let limbs = to_limbs(p);
        // Newton iteration for p^{-1} mod 2^64.
        let mut inv: u64 = 1;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(limbs[0].wrapping_mul(inv)));
        }
        let r = BigUint::one() << (64 * N);
        let mut f = PrimeField::<N> {
            p: limbs,
            n0: inv.wrapping_neg(),
            r2: Element::<N>(to_limbs(&((&r * &r) % p))),
            one: Element::<N>(to_limbs(&(&r % p))),
            p_big: p.clone(),
            half: (p - 1u8) >> 1,
            ts_s: 0,
            ts_q: BigUint::zero(),
            ts_z: Element::<N>([0; N]),
            nonresidue: 0,
        };
        let pm1 = p - 1u8;
        let s = pm1.trailing_zeros().unwrap_or(0) as u32;
        f.ts_s = s;
        f.ts_q = &pm1 >> s;
        let mut z = 2u64;
        while f.legendre(&f.from_u64(z)) != -1 {
            z += 1;
        }
        f.nonresidue = z;
        f.ts_z = f.pow(&f.from_u64(z), &f.ts_q.clone());
        Some(f)
    }

    pub fn modulus(&self) -> &BigUint {
        &self.p_big
    }

    /// The least quadratic non-residue, `≥ 2`.
    pub fn least_nonresidue(&self) -> u64 {
        self.nonresidue
    }

    pub fn zero(&self) -> Element<N> {
        Element::<N>([0; N])
    }

    pub fn one(&self) -> Element<N> {
        self.one
    }

    pub fn is_zero(&self, a: &Element<N>) -> bool {
        a.0 == [0; N]
    }

    pub fn mul(&self, a: &Element<N>, b: &Element<N>) -> Element<N> {
        // CIOS Montgomery multiplication; the result is < 2p before the
        // final subtraction for any odd p < 2^256.
        let p = &self.p;
        let mut t = [0u64; 18];
        for i in 0..N {
            let mut c: u128 = 0;
            for j in 0..N {
                let s = t[j] as u128 + (a.0[j] as u128) * (b.0[i] as u128) + c;
                t[j] = s as u64;
                c = s >> 64;
            }
            let s = t[N] as u128 + c;
            t[N] = s as u64;
            t[N + 1] = (s >> 64) as u64;
            let m = t[0].wrapping_mul(self.n0);
            let s = t[0] as u128 + (m as u128) * (p[0] as u128);
            let mut c = s >> 64;
            for j in 1..N {
                let s = t[j] as u128 + (m as u128) * (p[j] as u128) + c;
                t[j - 1] = s as u64;
                c = s >> 64;
            }
            let s = t[N] as u128 + c;
            t[N - 1] = s as u64;
            t[N] = t[N + 1] + (s >> 64) as u64;
            t[N + 1] = 0;
        }
        let mut r = [0u64; N];
        r.copy_from_slice(&t[..N]);
        if t[N] != 0 || geq(&r, p) {
            Element::<N>(sub_limbs(&r, p).0)
        } else {
            Element::<N>(r)
        }
    }

    pub fn sqr(&self, a: &Element<N>) -> Element<N> {
        self.mul(a, a)
    }

    pub fn add(&self, a: &Element<N>, b: &Element<N>) -> Element<N> {
        let (s, carry) = add_limbs(&a.0, &b.0);
        if carry || geq(&s, &self.p) {
            Element::<N>(sub_limbs(&s, &self.p).0)
        } else {
            Element::<N>(s)
        }
    }

    pub fn sub(&self, a: &Element<N>, b: &Element<N>) -> Element<N> {
        let (d, borrow) = sub_limbs(&a.0, &b.0);
        if borrow {
            Element::<N>(add_limbs(&d, &self.p).0)
        } else {
            Element::<N>(d)
        }
    }

    pub fn neg(&self, a: &Element<N>) -> Element<N> {
        self.sub(&self.zero(), a)
    }

    pub fn from_u64(&self, v: u64) -> Element<N> {
        let mut raw = [0u64; N];
        raw[0] = if self.p_big.bits() <= 64 {
            v % self.p[0]
        } else {
            v
        };
        self.mul(&Element::<N>(raw), &self.r2)
    }

    pub fn from_i64(&self, v: i64) -> Element<N> {
        let m = self.from_u64(v.unsigned_abs());
        if v < 0 {
            self.neg(&m)
        } else {
            m
        }
    }

    pub fn from_u128(&self, v: u128) -> Element<N> {
        self.from_big(&BigUint::from(v))
    }

    pub fn from_big(&self, v: &BigUint) -> Element<N> {
        let v = v % &self.p_big;
        self.mul(&Element::<N>(to_limbs(&v)), &self.r2)
    }

    pub fn to_big(&self, a: &Element<N>) -> BigUint {
        {
            let mut unit = [0u64; N];
            unit[0] = 1;
            from_limbs(&self.mul(a, &Element::<N>(unit)).0)
        }
    }

    pub fn pow(&self, a: &Element<N>, e: &BigUint) -> Element<N> {
        let mut acc = self.one;
        for i in (0..e.bits()).rev() {
            acc = self.sqr(&acc);
            if e.bit(i) {
                acc = self.mul(&acc, a);
            }
        }
        acc
    }

    /// `a^{-1}`; `None` for zero.
    pub fn inv(&self, a: &Element<N>) -> Option<Element<N>> {
        if self.is_zero(a) {
            return None;
        }
        if *a == self.one {
            return Some(self.one);
        }
        let mut u = a.0;
        let mut v = self.p;
        let mut x = self.r2;
        let mut y = self.zero();
        let mut unit = [0u64; N];
        unit[0] = 1;
        while u != unit && v != unit {
            if u == [0; N] || v == [0; N] {
                return None;
            }
            while u[0] & 1 == 0 {
                shr1(&mut u, false);
                x = self.half_coefficient(x);
            }
            while v[0] & 1 == 0 {
                shr1(&mut v, false);
                y = self.half_coefficient(y);
            }
            if geq(&u, &v) {
                u = sub_limbs(&u, &v).0;
                x = self.sub(&x, &y);
            } else {
                v = sub_limbs(&v, &u).0;
                y = self.sub(&y, &x);
            }
        }
        Some(if u == unit { x } else { y })
    }

    fn half_coefficient(&self, mut a: Element<N>) -> Element<N> {
        let mut carry = false;
        if a.0[0] & 1 != 0 {
            let sum = add_limbs(&a.0, &self.p);
            a.0 = sum.0;
            carry = sum.1;
        }
        shr1(&mut a.0, carry);
        a
    }

    pub fn div(&self, a: &Element<N>, b: &Element<N>) -> Option<Element<N>> {
        Some(self.mul(a, &self.inv(b)?))
    }

    /// `1`, `-1`, or `0` for zero.
    pub fn legendre(&self, a: &Element<N>) -> i32 {
        if self.is_zero(a) {
            return 0;
        }
        if self.pow(a, &self.half) == self.one {
            1
        } else {
            -1
        }
    }

    /// A square root by Tonelli–Shanks, or `None` for a non-residue.
    pub fn sqrt(&self, a: &Element<N>) -> Option<Element<N>> {
        match self.legendre(a) {
            0 => return Some(self.zero()),
            -1 => return None,
            _ => {}
        }
        let mut m = self.ts_s;
        let mut c = self.ts_z;
        let mut t = self.pow(a, &self.ts_q);
        let mut r = self.pow(a, &((&self.ts_q + 1u8) >> 1));
        while t != self.one {
            let mut i = 0;
            let mut t2 = t;
            while t2 != self.one {
                t2 = self.sqr(&t2);
                i += 1;
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = self.sqr(&b);
            }
            m = i;
            c = self.sqr(&b);
            t = self.mul(&t, &c);
            r = self.mul(&r, &b);
        }
        Some(r)
    }

    /// The canonical integer of `a` as lower-case hex, no prefix.
    pub fn hex(&self, a: &Element<N>) -> String {
        format!("{:x}", self.to_big(a))
    }
}

/// Probable-prime screening with twenty fixed Miller–Rabin bases.
/// This is not a universal primality proof for arbitrary large inputs.
pub fn is_probable_prime(n: &BigUint) -> bool {
    let small = [
        2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    ];
    if n < &BigUint::from(2u8) {
        return false;
    }
    for &q in &small {
        let q = BigUint::from(q);
        if n == &q {
            return true;
        }
        if (n % &q).is_zero() {
            return false;
        }
    }
    let nm1 = n - 1u8;
    let s = nm1.trailing_zeros().unwrap_or(0);
    let d = &nm1 >> s;
    'base: for &a in &small {
        let mut x = BigUint::from(a).modpow(&d, n);
        if x.is_one() || x == nm1 {
            continue;
        }
        for _ in 1..s {
            x = x.modpow(&BigUint::from(2u8), n);
            if x == nm1 {
                continue 'base;
            }
        }
        return false;
    }
    true
}
