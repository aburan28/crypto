//! Curve-generic Montgomery-form prime-field arithmetic for the native
//! GLV-chain benchmark (`research/endo_native_bench_20261016`).
//!
//! [`GenericField<L>`] holds an odd prime `p < 2^(64·L)` and its Montgomery
//! constants, all computed at run time from `p` (the modulus comes from a
//! frozen chain-set file, not from a compile-time constant), and [`Fe<L>`] is
//! an element in Montgomery form `a·R mod p`, `R = 2^(64·L)`.  Arithmetic
//! delegates to the crate's [`crate::ct_bignum::Uint`] primitives
//! (`mont_mul`, `mont_sqr`, `add_mod`, `sub_mod`), exactly as
//! [`crate::ecc::cryptopro_b_field`] does for its fixed prime; the only
//! difference is that `p` is a run-time value, so the compiler cannot fold
//! it.  Four limbs serve 255/256-bit fields, six serve 384-bit fields and
//! thirteen serve the 782-bit field of CP6-782.
//!
//! The field operations are the crate's constant-time ones, but
//! [`GenericField::inv`] scans the fixed public exponent `p - 2` with a
//! plain branch (public constant, never the input), and the point code
//! built on this module is variable time throughout: a measurement harness,
//! never for secrets.
//!
//! With the `endo-chain-opcount` cargo feature, every multiplication,
//! squaring, inversion and addition on the current thread is counted
//! ([`opcount`]), the same accounting as the `cryptopro-b-opcount` counter:
//! an inversion is also counted through the squarings and multiplications it
//! performs, and conversions into and out of Montgomery form are not counted.

use crate::ct_bignum::{compute_minv64, Uint};
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// Field-operation counters for the counting run of `endo-chain-bench`,
/// compiled in only with the `endo-chain-opcount` feature.  Without it none
/// of this exists and the arithmetic is unchanged.
#[cfg(feature = "endo-chain-opcount")]
pub mod opcount {
    use std::cell::Cell;

    /// Field operations on this thread since the last [`reset`].
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
    pub struct OpCounts {
        /// Montgomery multiplications, including those inside inversions.
        pub mul: u64,
        /// Montgomery squarings, including those inside inversions.
        pub sqr: u64,
        /// Inversions.
        pub inv: u64,
        /// Additions, subtractions and negations.
        pub add: u64,
    }

    thread_local! {
        static COUNTS: Cell<OpCounts> = const {
            Cell::new(OpCounts { mul: 0, sqr: 0, inv: 0, add: 0 })
        };
    }

    #[inline]
    pub(crate) fn bump(f: impl FnOnce(&mut OpCounts)) {
        COUNTS.with(|c| {
            let mut v = c.get();
            f(&mut v);
            c.set(v);
        });
    }

    /// Zero this thread's counters.
    pub fn reset() {
        COUNTS.with(|c| c.set(OpCounts::default()));
    }

    /// This thread's counters.
    pub fn read() -> OpCounts {
        COUNTS.with(|c| c.get())
    }
}

/// A field element in Montgomery form, canonical (`< p`).  Only meaningful
/// together with the [`GenericField`] that made it.
#[derive(Copy, Clone, Debug)]
pub struct Fe<const L: usize>(Uint<L>);

impl<const L: usize> PartialEq for Fe<L> {
    fn eq(&self, other: &Self) -> bool {
        self.0 .0 == other.0 .0
    }
}

impl<const L: usize> Eq for Fe<L> {}

impl<const L: usize> Fe<L> {
    /// Zero (the same in every representation).
    pub const ZERO: Self = Fe(Uint::<L>::ZERO);

    pub fn is_zero(&self) -> bool {
        self.0 .0.iter().all(|&l| l == 0)
    }

    /// The Montgomery-form limbs.
    pub fn montgomery_limbs(&self) -> &[u64; L] {
        &self.0 .0
    }
}

/// An odd prime field `F_p` with `p < 2^(64·L)`, with its Montgomery
/// constants.
#[derive(Clone, Debug)]
pub struct GenericField<const L: usize> {
    p: Uint<L>,
    p_inv_low: u64,
    /// `R mod p`: one in Montgomery form.
    one: Uint<L>,
    /// `R² mod p`.
    r2: Uint<L>,
    /// `p - 2`, the Fermat exponent, and its bit length.
    exp: Uint<L>,
    exp_bits: usize,
    p_big: BigUint,
}

impl<const L: usize> GenericField<L> {
    /// The field of the odd modulus `p`, which must satisfy `3 < p < 2^(64·L)`.
    /// Primality is the caller's input (the chain-set loader runs a
    /// probable-prime test); nothing here depends on it except `inv`.
    pub fn new(p: &BigUint) -> Result<Self, String> {
        if p.bits() as usize > 64 * L {
            return Err(format!(
                "p has {} bits, more than the {} of {} limbs",
                p.bits(),
                64 * L,
                L
            ));
        }
        if p <= &BigUint::from(3u8) || (p % 2u8).is_zero() {
            return Err("p must be an odd prime above 3".into());
        }
        let r = BigUint::one() << (64 * L);
        let pu = Uint::<L>::from_biguint(p);
        let one = Uint::<L>::from_biguint(&(&r % p));
        let r2 = Uint::<L>::from_biguint(&((&r * &r) % p));
        let e = p - 2u8;
        Ok(GenericField {
            p: pu,
            p_inv_low: compute_minv64(pu.0[0]),
            one,
            r2,
            exp: Uint::<L>::from_biguint(&e),
            exp_bits: e.bits() as usize,
            p_big: p.clone(),
        })
    }

    pub fn modulus(&self) -> &BigUint {
        &self.p_big
    }

    pub fn one(&self) -> Fe<L> {
        Fe(self.one)
    }

    pub fn zero(&self) -> Fe<L> {
        Fe(Uint::<L>::ZERO)
    }

    /// Montgomery form of `v mod p`.
    pub fn from_biguint(&self, v: &BigUint) -> Fe<L> {
        let c = Uint::<L>::from_biguint(&(v % &self.p_big));
        Fe(Uint::mont_mul(&c, &self.r2, &self.p, self.p_inv_low))
    }

    /// Montgomery form of a canonical value; `Err` unless `v < p`.
    pub fn from_canonical(&self, v: &BigUint) -> Result<Fe<L>, String> {
        if v >= &self.p_big {
            return Err("value is not below the modulus".into());
        }
        Ok(self.from_biguint(v))
    }

    /// The signed integer `v` (e.g. a small curve coefficient) mod p.
    pub fn from_i64(&self, v: i64) -> Fe<L> {
        let m = self.from_biguint(&BigUint::from(v.unsigned_abs()));
        if v < 0 {
            Fe(Uint::<L>::ZERO.sub_mod(&m.0, &self.p))
        } else {
            m
        }
    }

    pub fn to_biguint(&self, a: &Fe<L>) -> BigUint {
        Uint::mont_mul(&a.0, &Uint::<L>::ONE, &self.p, self.p_inv_low).to_biguint()
    }

    #[inline]
    pub fn add(&self, a: &Fe<L>, b: &Fe<L>) -> Fe<L> {
        #[cfg(feature = "endo-chain-opcount")]
        opcount::bump(|c| c.add += 1);
        Fe(a.0.add_mod(&b.0, &self.p))
    }

    #[inline]
    pub fn sub(&self, a: &Fe<L>, b: &Fe<L>) -> Fe<L> {
        #[cfg(feature = "endo-chain-opcount")]
        opcount::bump(|c| c.add += 1);
        Fe(a.0.sub_mod(&b.0, &self.p))
    }

    #[inline]
    pub fn neg(&self, a: &Fe<L>) -> Fe<L> {
        #[cfg(feature = "endo-chain-opcount")]
        opcount::bump(|c| c.add += 1);
        Fe(Uint::<L>::ZERO.sub_mod(&a.0, &self.p))
    }

    #[inline]
    pub fn mul(&self, a: &Fe<L>, b: &Fe<L>) -> Fe<L> {
        #[cfg(feature = "endo-chain-opcount")]
        opcount::bump(|c| c.mul += 1);
        Fe(Uint::mont_mul(&a.0, &b.0, &self.p, self.p_inv_low))
    }

    #[inline]
    pub fn sqr(&self, a: &Fe<L>) -> Fe<L> {
        #[cfg(feature = "endo-chain-opcount")]
        opcount::bump(|c| c.sqr += 1);
        Fe(Uint::mont_sqr(&a.0, &self.p, self.p_inv_low))
    }

    /// `a^(p-2)` by left-to-right square-and-multiply over every bit of
    /// `p - 2` from its top bit, starting from one: `bits(p-2)` squarings and
    /// `popcount(p-2)` multiplications, for every `a` (the branch is on the
    /// public exponent).  For CryptoPro-B that is 256 S + 8 M, the count of
    /// [`crate::ecc::cryptopro_b_field::CryptoProBFieldElement::inv`].
    pub fn inv(&self, a: &Fe<L>) -> Fe<L> {
        #[cfg(feature = "endo-chain-opcount")]
        opcount::bump(|c| c.inv += 1);
        let mut acc = self.one();
        for i in (0..self.exp_bits).rev() {
            acc = self.sqr(&acc);
            if (self.exp.0[i / 64] >> (i % 64)) & 1 == 1 {
                acc = self.mul(&acc, a);
            }
        }
        acc
    }

    /// `(multiplications, squarings)` of one [`Self::inv`].
    pub fn inv_cost(&self) -> (u64, u64) {
        let pop: u32 = self.exp.0.iter().map(|l| l.count_ones()).sum();
        (u64::from(pop), self.exp_bits as u64)
    }

    /// `a · c` for a small signed integer `c`, by additions (double and
    /// add on `|c|`, then a negation): no multiplication is counted.
    pub fn mul_small(&self, a: &Fe<L>, c: i64) -> Fe<L> {
        let mut m = c.unsigned_abs();
        let mut acc: Option<Fe<L>> = None;
        let mut base = *a;
        while m > 0 {
            if m & 1 == 1 {
                acc = Some(match acc {
                    None => base,
                    Some(x) => self.add(&x, &base),
                });
            }
            m >>= 1;
            if m > 0 {
                base = self.add(&base, &base);
            }
        }
        let r = acc.unwrap_or_else(|| self.zero());
        if c < 0 {
            self.neg(&r)
        } else {
            r
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn primes() -> Vec<BigUint> {
        let h = |s: &str| BigUint::parse_bytes(s.as_bytes(), 16).unwrap();
        vec![
            // GOST CryptoPro-B, 2^255 + 3225
            (BigUint::one() << 255u32) + 3225u32,
            // the BLS12-381 scalar field (Bandersnatch's base field), 255 bits
            h("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"),
            // P-384's prime
            h("fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff"),
        ]
    }

    fn samples(p: &BigUint) -> Vec<BigUint> {
        vec![
            BigUint::zero(),
            BigUint::one(),
            BigUint::from(2u8),
            BigUint::from(0xdead_beef_u64),
            BigUint::one() << 64u32,
            BigUint::one() << 200u32,
            p - 1u8,
            p - 2u8,
            p / 2u8,
            p / 3u8 + 12345u32,
        ]
    }

    fn check<const L: usize>(p: &BigUint) {
        let f = GenericField::<L>::new(p).unwrap();
        for a in samples(p) {
            let fa = f.from_biguint(&a);
            assert_eq!(f.to_biguint(&fa), &a % p);
            assert_eq!(f.to_biguint(&f.neg(&fa)), (p - (&a % p)) % p);
            assert_eq!(
                f.to_biguint(&f.mul_small(&fa, -5)),
                (p * 5u8 - (&a % p) * 5u8) % p
            );
            assert_eq!(f.to_biguint(&f.mul_small(&fa, 7)), (&a * 7u8) % p);
            for b in samples(p) {
                let fb = f.from_biguint(&b);
                assert_eq!(f.to_biguint(&f.add(&fa, &fb)), (&a + &b) % p);
                assert_eq!(
                    f.to_biguint(&f.sub(&fa, &fb)),
                    (p + (&a % p) - (&b % p)) % p
                );
                assert_eq!(f.to_biguint(&f.mul(&fa, &fb)), (&a * &b) % p);
            }
            assert_eq!(f.sqr(&fa), f.mul(&fa, &fa));
            if !(&a % p).is_zero() {
                assert_eq!(f.to_biguint(&f.inv(&fa)), a.modpow(&(p - 2u8), p));
                assert_eq!(f.mul(&fa, &f.inv(&fa)), f.one());
            }
        }
        assert_eq!(f.from_i64(-3), f.from_biguint(&(p - 3u8)));
        assert!(f.from_canonical(p).is_err());
    }

    #[test]
    fn arithmetic_matches_biguint_at_four_and_six_limbs() {
        let ps = primes();
        check::<4>(&ps[0]);
        check::<4>(&ps[1]);
        check::<6>(&ps[2]);
        // a 255-bit prime in six limbs too (p far below R)
        check::<6>(&ps[1]);
    }

    #[test]
    fn arithmetic_matches_biguint_at_thirteen_limbs() {
        // A 782-bit odd modulus (not necessarily prime): the ring identities
        // hold for any odd modulus; inversion needs a prime and is tested at
        // four and six limbs.
        let m = (BigUint::one() << 781u32) + 0x1_2345u32;
        let f = GenericField::<13>::new(&m).unwrap();
        for a in samples(&m) {
            let fa = f.from_biguint(&a);
            assert_eq!(f.to_biguint(&fa), &a % &m);
            for b in samples(&m) {
                let fb = f.from_biguint(&b);
                assert_eq!(f.to_biguint(&f.mul(&fa, &fb)), (&a * &b) % &m);
                assert_eq!(f.to_biguint(&f.add(&fa, &fb)), (&a + &b) % &m);
                assert_eq!(
                    f.to_biguint(&f.sub(&fa, &fb)),
                    (&m + (&a % &m) - (&b % &m)) % &m
                );
            }
        }
    }

    #[test]
    fn cryptopro_b_inversion_cost_is_256_squarings_and_8_multiplications() {
        let f = GenericField::<4>::new(&primes()[0]).unwrap();
        assert_eq!(f.inv_cost(), (8, 256));
    }

    #[test]
    fn oversized_or_even_moduli_are_refused() {
        assert!(GenericField::<4>::new(&(BigUint::one() << 256u32)).is_err());
        assert!(GenericField::<4>::new(&BigUint::from(1000u32)).is_err());
        assert!(GenericField::<4>::new(&((BigUint::one() << 256u32) + 1u8)).is_err());
    }

    #[cfg(feature = "endo-chain-opcount")]
    #[test]
    fn counters_count_inversions_through_their_operations() {
        let f = GenericField::<4>::new(&primes()[0]).unwrap();
        let a = f.from_biguint(&BigUint::from(12345u32));
        opcount::reset();
        let _ = f.inv(&a);
        let c = opcount::read();
        assert_eq!((c.mul, c.sqr, c.inv), (8, 256, 1));
    }
}
