//! Montgomery-form field arithmetic over the prime field of the GOST R
//! 34.10-2001 **CryptoPro-B** parameter set (RFC 4357 §11.4.4).
//!
//! `p = 2^255 + 3225`
//!
//! Mirrors [`crate::ecc::p256_field`]: every value is held in Montgomery
//! form `a·R mod p` where `R = 2^256`, and arithmetic delegates to the
//! [`crate::ct_bignum::U256`] primitives.  `m' = -p^(-1) mod 2^64` is
//! computed at compile time via [`crate::ct_bignum::compute_minv64`];
//! `R mod p` and `R² mod p` were computed once with `num-bigint` and are
//! hardcoded below, and this module's tests recompute them and fail if
//! the constants ever drift.
//!
//! The module exists for the CryptoPro-B GLV measurement in
//! [`crate::ecc::cryptopro_b_point`] (a variable-time benchmark harness).
//! The add/sub/mul/sqr primitives are the crate's constant-time ones, but
//! [`CryptoProBFieldElement::inv`] scans the fixed public exponent `p - 2`
//! with a plain branch, which is fine for a public constant and is not
//! the cmov ladder `p256_field` uses.

use crate::ct_bignum::{compute_minv64, Uint, U256};
use num_bigint::BigUint;
use subtle::{Choice, ConstantTimeEq};

/// CryptoPro-B prime modulus, as little-endian limbs.
///
/// `p = 0x8000000000000000000000000000000000000000000000000000000000000C99`
/// `  = 2^255 + 3225`
pub const P: U256 = Uint([
    0x0000_0000_0000_0C99,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x8000_0000_0000_0000,
]);

/// Montgomery constant `m' = -p^(-1) mod 2^64`.
pub const P_INV_LOW: u64 = compute_minv64(P.0[0]);

/// `R mod p`, where `R = 2^256`.  Equal to "1" in Montgomery form.
///
/// `2^256 - p = 2^255 - 3225`.
pub const R_MOD_P: U256 = Uint([
    0xFFFF_FFFF_FFFF_F367,
    0xFFFF_FFFF_FFFF_FFFF,
    0xFFFF_FFFF_FFFF_FFFF,
    0x7FFF_FFFF_FFFF_FFFF,
]);

/// `R^2 mod p`.  Used to convert canonical values into Montgomery form
/// via `mont_mul(a, R^2 mod p)`.
///
/// `2^255 ≡ -3225`, so `2^512 ≡ 4·3225² = 41602500 = 0x27ACDC4`.
pub const R_SQUARED_MOD_P: U256 = Uint([
    0x0000_0000_027A_CDC4,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
]);

/// A field element over the CryptoPro-B prime, stored in Montgomery form.
#[derive(Copy, Clone, Debug)]
pub struct CryptoProBFieldElement(U256);

impl CryptoProBFieldElement {
    pub const ZERO: CryptoProBFieldElement = CryptoProBFieldElement(U256::ZERO);
    pub const ONE: CryptoProBFieldElement = CryptoProBFieldElement(R_MOD_P);

    /// Montgomery form of a canonical value `a < p`.
    pub fn from_canonical(a: &U256) -> Self {
        CryptoProBFieldElement(U256::mont_mul(a, &R_SQUARED_MOD_P, &P, P_INV_LOW))
    }

    /// Montgomery form of a canonical little-endian limb array `a < p`
    /// (the encoding of every field element in
    /// [`crate::ecc::cryptopro_b_chain_consts`]).
    pub fn from_limbs(a: &[u64; 4]) -> Self {
        Self::from_canonical(&Uint(*a))
    }

    pub fn to_canonical(&self) -> U256 {
        U256::mont_mul(&self.0, &U256::ONE, &P, P_INV_LOW)
    }

    pub fn from_biguint(v: &BigUint) -> Self {
        let p_bu = P.to_biguint();
        let reduced = U256::from_biguint(&(v % &p_bu));
        Self::from_canonical(&reduced)
    }

    pub fn to_biguint(&self) -> BigUint {
        self.to_canonical().to_biguint()
    }

    pub fn to_bytes_be(&self) -> [u8; 32] {
        self.to_canonical().to_bytes_be()
    }

    pub fn from_bytes_be(bytes: &[u8; 32]) -> Self {
        Self::from_canonical(&U256::from_bytes_be(bytes))
    }

    pub fn add(&self, other: &Self) -> Self {
        CryptoProBFieldElement(self.0.add_mod(&other.0, &P))
    }

    pub fn sub(&self, other: &Self) -> Self {
        CryptoProBFieldElement(self.0.sub_mod(&other.0, &P))
    }

    pub fn mul(&self, other: &Self) -> Self {
        CryptoProBFieldElement(U256::mont_mul(&self.0, &other.0, &P, P_INV_LOW))
    }

    pub fn sqr(&self) -> Self {
        CryptoProBFieldElement(U256::mont_sqr(&self.0, &P, P_INV_LOW))
    }

    pub fn neg(&self) -> Self {
        CryptoProBFieldElement(U256::ZERO.sub_mod(&self.0, &P))
    }

    /// Multiplicative inverse via Fermat: `a^(p-2) mod p`.
    ///
    /// Square-and-multiply over the fixed exponent `p - 2 = 2^255 + 3223`,
    /// which has only eight set bits, so an inversion costs 256 squarings
    /// and 8 multiplications.  The branch is on the public constant `p`,
    /// never on `a`; the per-input work is the same for every `a`.  (The
    /// cmov ladder of [`crate::ecc::p256_field::P256FieldElement::inv`]
    /// would spend 256 multiplications on the same exponent.)
    pub fn inv(&self) -> Self {
        // exp = p - 2.  P[0] = 0xC99, so subtract 2 from the low limb
        // without affecting the higher ones.
        let exp = Uint([P.0[0] - 2, P.0[1], P.0[2], P.0[3]]);
        let mut acc = CryptoProBFieldElement::ONE;
        for i in (0..256).rev() {
            acc = acc.sqr();
            if (exp.0[i / 64] >> (i % 64)) & 1 == 1 {
                acc = acc.mul(self);
            }
        }
        acc
    }

    pub fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq_full(&other.0)
    }

    pub fn ct_is_zero(&self) -> Choice {
        self.0.ct_is_zero()
    }

    pub fn is_zero(&self) -> bool {
        bool::from(self.ct_is_zero())
    }

    pub fn cmov(a: &Self, b: &Self, choice: Choice) -> Self {
        CryptoProBFieldElement(U256::cmov(&a.0, &b.0, choice))
    }

    #[inline]
    pub fn as_montgomery(&self) -> &U256 {
        &self.0
    }

    pub const fn from_montgomery_unchecked(limbs: U256) -> Self {
        CryptoProBFieldElement(limbs)
    }
}

impl PartialEq for CryptoProBFieldElement {
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl Eq for CryptoProBFieldElement {}

impl ConstantTimeEq for CryptoProBFieldElement {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.ct_eq(other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use num_traits::{One, Zero};

    fn p_bu() -> BigUint {
        P.to_biguint()
    }

    fn samples() -> Vec<BigUint> {
        let p = p_bu();
        vec![
            BigUint::zero(),
            BigUint::one(),
            BigUint::from(2u8),
            BigUint::from(3225u32),
            BigUint::from(0xdeadbeefu64),
            BigUint::from(1u8) << 32,
            BigUint::from(1u8) << 64,
            BigUint::from(1u8) << 128,
            BigUint::from(1u8) << 192,
            BigUint::from(1u8) << 255, // the top bit of p is set, so this is < p
            (BigUint::from(1u8) << 255) - 1u8,
            &p - 1u8,
            &p - 2u8,
            &p / 2u8,
            (&p / 2u8) + 1u8,
            BigUint::parse_bytes(
                b"3E1AF419A269A5F866A7D3C25C3DF80AE979259373FF2B182F49D4CE7E1BBC8B",
                16,
            )
            .unwrap(),
            BigUint::parse_bytes(
                b"7EDCBA9876543210FEDCBA9876543210FEDCBA9876543210FEDCBA9876543210",
                16,
            )
            .unwrap(),
        ]
    }

    #[test]
    fn p_constant_matches_curve_params() {
        let want = U256::from_biguint(&crate::ecc::CurveParams::gost_cryptopro_b().p);
        assert!(bool::from(P.ct_eq_full(&want)));
    }

    #[test]
    fn p_constant_matches_chain_constants() {
        let want = Uint(crate::ecc::cryptopro_b_chain_consts::P);
        assert!(bool::from(P.ct_eq_full(&want)));
        let p_hex = crate::ecc::cryptopro_b_chain_consts::P_HEX.trim_start_matches("0x");
        assert_eq!(BigUint::parse_bytes(p_hex.as_bytes(), 16).unwrap(), p_bu());
    }

    #[test]
    fn p_is_two_to_the_255_plus_3225() {
        assert_eq!(p_bu(), (BigUint::one() << 255) + BigUint::from(3225u32));
    }

    #[test]
    fn r_mod_p_constant_is_correct() {
        let r = BigUint::one() << 256;
        let want = U256::from_biguint(&(&r % p_bu()));
        assert!(bool::from(R_MOD_P.ct_eq_full(&want)));
    }

    #[test]
    fn r_squared_constant_is_correct() {
        let r = BigUint::one() << 256;
        let want = U256::from_biguint(&((&r * &r) % p_bu()));
        assert!(bool::from(R_SQUARED_MOD_P.ct_eq_full(&want)));
    }

    #[test]
    fn p_inv_low_consistent_with_compute_minv64() {
        let prod = P.0[0].wrapping_mul(P_INV_LOW);
        assert_eq!(prod.wrapping_add(1), 0);
    }

    #[test]
    fn from_to_canonical_roundtrips() {
        for v in samples() {
            let fe = CryptoProBFieldElement::from_biguint(&v);
            assert_eq!(fe.to_biguint(), v);
            let limbs = U256::from_biguint(&v);
            assert_eq!(CryptoProBFieldElement::from_limbs(&limbs.0).to_biguint(), v);
        }
    }

    #[test]
    fn bytes_roundtrip() {
        for v in samples() {
            let fe = CryptoProBFieldElement::from_biguint(&v);
            let back = CryptoProBFieldElement::from_bytes_be(&fe.to_bytes_be());
            assert_eq!(back, fe);
        }
    }

    #[test]
    fn add_matches_biguint() {
        let p = p_bu();
        for a in samples() {
            for b in samples() {
                let fa = CryptoProBFieldElement::from_biguint(&a);
                let fb = CryptoProBFieldElement::from_biguint(&b);
                assert_eq!(fa.add(&fb).to_biguint(), (&a + &b) % &p);
            }
        }
    }

    #[test]
    fn sub_matches_biguint() {
        let p = p_bu();
        for a in samples() {
            for b in samples() {
                let fa = CryptoProBFieldElement::from_biguint(&a);
                let fb = CryptoProBFieldElement::from_biguint(&b);
                let want = (&p + (&a % &p) - (&b % &p)) % &p;
                assert_eq!(fa.sub(&fb).to_biguint(), want);
            }
        }
    }

    #[test]
    fn mul_matches_biguint() {
        let p = p_bu();
        for a in samples() {
            for b in samples() {
                let fa = CryptoProBFieldElement::from_biguint(&a);
                let fb = CryptoProBFieldElement::from_biguint(&b);
                assert_eq!(fa.mul(&fb).to_biguint(), (&a * &b) % &p);
            }
        }
    }

    #[test]
    fn sqr_matches_mul_self() {
        let p = p_bu();
        for a in samples() {
            let fa = CryptoProBFieldElement::from_biguint(&a);
            assert_eq!(fa.sqr().to_biguint(), (&a * &a) % &p);
            assert_eq!(fa.sqr(), fa.mul(&fa));
        }
    }

    #[test]
    fn neg_matches_biguint() {
        let p = p_bu();
        for a in samples() {
            let fa = CryptoProBFieldElement::from_biguint(&a);
            let want = (&p - (&a % &p)) % &p;
            assert_eq!(fa.neg().to_biguint(), want);
            assert_eq!(fa.add(&fa.neg()), CryptoProBFieldElement::ZERO);
        }
    }

    #[test]
    fn inv_times_self_is_one() {
        for a in samples() {
            if a.is_zero() {
                continue;
            }
            let fa = CryptoProBFieldElement::from_biguint(&a);
            assert_eq!(fa.mul(&fa.inv()).to_biguint(), BigUint::one());
            // Round trip through the inverse and back.
            assert_eq!(fa.inv().inv(), fa);
        }
    }

    #[test]
    fn inv_matches_biguint_modpow() {
        let p = p_bu();
        let exp = &p - 2u8;
        for a in samples() {
            if a.is_zero() {
                continue;
            }
            let fa = CryptoProBFieldElement::from_biguint(&a);
            assert_eq!(fa.inv().to_biguint(), a.modpow(&exp, &p));
        }
    }
}
