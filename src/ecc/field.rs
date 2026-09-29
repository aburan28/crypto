//! Finite field arithmetic over F_p (prime field).
//!
//! Every `FieldElement` stores its value reduced mod p and carries p along.
//! The modulus p must be prime for `inv()` to work correctly (Fermat's
//! little theorem: a^(p-2) ≡ a^(-1) mod p).

use crate::utils::mod_pow_ct;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::fmt;

#[derive(Clone, Debug)]
pub struct FieldElement {
    /// Value in [0, p)
    pub value: BigUint,
    /// The field prime p
    pub modulus: BigUint,
}

impl FieldElement {
    /// Construct a new field element, reducing `value` mod `modulus`.
    pub fn new(value: BigUint, modulus: BigUint) -> Self {
        let value = value % &modulus;
        FieldElement { value, modulus }
    }

    pub fn zero(modulus: BigUint) -> Self {
        FieldElement {
            value: BigUint::zero(),
            modulus,
        }
    }

    pub fn one(modulus: BigUint) -> Self {
        FieldElement {
            value: BigUint::one(),
            modulus,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.value.is_zero()
    }

    // ── Arithmetic ────────────────────────────────────────────────────────────

    pub fn add(&self, rhs: &Self) -> Self {
        Self::new(&self.value + &rhs.value, self.modulus.clone())
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        // Always compute (a + p - b) mod p instead of branching on
        // (a >= b).  Algebraically the same as (a - b) mod p, but
        // without the data-dependent comparison that the branch-and-
        // negate variant uses.  The underlying `BigUint` add/sub still
        // leak via limb count — see SECURITY.md.
        let v = (&self.value + &self.modulus - &rhs.value) % &self.modulus;
        FieldElement {
            value: v,
            modulus: self.modulus.clone(),
        }
    }

    /// [`sub`](Self::sub) in **variable time**, for public values only:
    /// `a − b` when `a ≥ b`, else `a + p − b`.  For reduced operands that
    /// is exactly `sub`'s `(a + p − b) mod p`, without its division; the
    /// comparison it branches on is the one `sub` avoids, so protocol code
    /// keeps `sub`.  Used by the `*_vartime` point operations.
    ///
    /// An unreduced operand (no constructor makes one, but the fields are
    /// public) can leave that difference at or above `p`; it is reduced
    /// then, so the value is `sub`'s for every input.
    pub fn sub_vartime(&self, rhs: &Self) -> Self {
        let mut value = if self.value >= rhs.value {
            &self.value - &rhs.value
        } else {
            &self.value + &self.modulus - &rhs.value
        };
        if value >= self.modulus {
            value %= &self.modulus;
        }
        FieldElement {
            value,
            modulus: self.modulus.clone(),
        }
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        Self::new(&self.value * &rhs.value, self.modulus.clone())
    }

    /// Additive inverse: -a mod p.  Always computes (p - a) mod p
    /// instead of branching on `a == 0`.  For `a = 0` the result is
    /// `p mod p = 0`; for `a != 0` it is `p - a` (unreduced, since
    /// `0 < p - a < p`).
    pub fn neg(&self) -> Self {
        let v = (&self.modulus - &self.value) % &self.modulus;
        FieldElement {
            value: v,
            modulus: self.modulus.clone(),
        }
    }

    /// [`neg`](Self::neg) in **variable time**, for public values only:
    /// `0` for zero, else `p − a`, the value `neg` computes without its
    /// division, and on words when `p` fits one.  Branches on `a = 0`,
    /// which `neg` avoids, so protocol code keeps `neg`.
    pub fn neg_vartime(&self) -> Self {
        let value = match (single_word(&self.value), single_word(&self.modulus)) {
            (Some(a), Some(p)) if a <= p => BigUint::from(if a == 0 { 0 } else { p - a }),
            _ if self.value.is_zero() => BigUint::zero(),
            _ => &self.modulus - &self.value,
        };
        FieldElement {
            value,
            modulus: self.modulus.clone(),
        }
    }

    /// Modular exponentiation via the Montgomery-ladder
    /// [`mod_pow_ct`].  Iteration count is fixed at `modulus.bits()`,
    /// and every step performs one multiplication and one squaring
    /// regardless of the bit value.  Caveat: the underlying
    /// `BigUint` arithmetic is still not constant-time at the limb
    /// level — see SECURITY.md.
    pub fn pow(&self, exp: &BigUint) -> Self {
        let bits = self.modulus.bits() as usize;
        let value = mod_pow_ct(&self.value, exp, &self.modulus, bits);
        FieldElement {
            value,
            modulus: self.modulus.clone(),
        }
    }

    /// Multiplicative inverse via Fermat's little theorem: a^(p-2) mod p.
    /// Returns `None` if `self` is zero.  Uses the laddered [`pow`].
    pub fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        let exp = &self.modulus - BigUint::from(2u32);
        Some(self.pow(&exp))
    }

    /// Multiplicative inverse in **variable time**; `None` if `self` is
    /// zero.  For public values only.
    ///
    /// [`inv`](Self::inv) exponentiates through the constant-time ladder,
    /// three multiplications per bit of `p` whatever the value, because
    /// protocol code (ECDSA, Schnorr, ElGamal, the proofs) inverts values
    /// that depend on secrets.  The cryptanalysis walks invert public
    /// coordinates one step at a time, and there that ladder was 80–90%
    /// of the run.  This runs the extended Euclidean algorithm instead —
    /// on machine words when `p` fits one, through
    /// [`BigUint::modinv`] otherwise — whose running time depends on the
    /// value.  Never call it on anything derived from a secret; the
    /// callers are the `*_vartime` point operations, which carry the same
    /// rule, and the cryptanalysis group's batch inversion.
    ///
    /// Over a prime `p`, the modulus this type is for (see the module
    /// note), a nonzero element has exactly one inverse, so this returns
    /// exactly `self.inv()`.  Where Euclid finds none — a nonzero multiple
    /// of `p` left unreduced, or a non-unit when `p` is not prime — it
    /// returns `inv`'s value instead, so it is `None` exactly when `inv`
    /// is and the operations built on it panic nowhere `inv`'s callers
    /// did not.  The one input on which the two differ is a unit modulo
    /// a composite `p`: this returns its inverse there, `inv` returns
    /// `a^(p−2)`, which is not one.
    pub fn inv_vartime(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        let value = match (single_word(&self.modulus), single_word(&self.value)) {
            (Some(p), Some(a)) if a < p => inv_mod_u64(a, p).map(BigUint::from),
            _ => self.value.modinv(&self.modulus),
        };
        match value {
            Some(value) => Some(FieldElement {
                value,
                modulus: self.modulus.clone(),
            }),
            None => self.inv(),
        }
    }
}

/// The value of `n` when it fits one 64-bit word.
pub(crate) fn single_word(n: &BigUint) -> Option<u64> {
    let mut digits = n.iter_u64_digits();
    let lo = digits.next().unwrap_or(0);
    match digits.next() {
        None => Some(lo),
        Some(_) => None,
    }
}

/// `a⁻¹ mod m` for `0 < a < m`, or `None` when `gcd(a, m) ≠ 1`.  Variable
/// time (see [`FieldElement::inv_vartime`]).
///
/// The extended Euclidean algorithm on `(m, a)`, keeping only `a`'s
/// coefficient `t`.  Those coefficients alternate in sign and never
/// exceed `m` in magnitude, so the loop carries their magnitudes in a
/// `u64` and recovers the sign of the last one from the parity of the
/// step count: no signed or double-width arithmetic.
pub(crate) fn inv_mod_u64(a: u64, m: u64) -> Option<u64> {
    let (mut r0, mut r1) = (m, a);
    let (mut t0, mut t1) = (0u64, 1u64);
    let mut odd = false;
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 + q * t1);
        odd = !odd;
    }
    if r0 != 1 {
        return None;
    }
    // `t0` is the coefficient from before the last step; it is positive
    // exactly when the step count is odd.
    Some(if odd { t0 } else { m - t0 })
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fe(v: u64, p: u64) -> FieldElement {
        FieldElement::new(BigUint::from(v), BigUint::from(p))
    }

    #[test]
    fn add_mod() {
        // In F_7: 5 + 4 = 9 ≡ 2
        assert_eq!(fe(5, 7).add(&fe(4, 7)), fe(2, 7));
    }

    #[test]
    fn sub_wrap() {
        // In F_7: 2 - 5 = -3 ≡ 4
        assert_eq!(fe(2, 7).sub(&fe(5, 7)), fe(4, 7));
    }

    #[test]
    fn mul_mod() {
        // In F_7: 3 * 5 = 15 ≡ 1
        assert_eq!(fe(3, 7).mul(&fe(5, 7)), fe(1, 7));
    }

    #[test]
    fn inv_fermat() {
        // In F_7: 3^(-1) = 5 because 3*5=15≡1 mod 7
        let inv = fe(3, 7).inv().unwrap();
        assert_eq!(inv, fe(5, 7));
    }

    /// `inv_vartime` is `inv` — on the word path, across word sizes up
    /// to the largest 64-bit prime, and on the `BigUint` path — and both
    /// refuse zero.
    #[test]
    fn inv_vartime_matches_inv() {
        let primes: [u64; 6] = [
            2,
            7,
            10_007,
            4_294_967_291,
            140_737_529_461_027,
            u64::MAX - 58, // 2^64 - 59, the largest 64-bit prime
        ];
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        for &p in &primes {
            assert!(fe(0, p).inv_vartime().is_none());
            for i in 0..200u64 {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                let v = match i {
                    0 => 1,
                    1 => p - 1,
                    _ => s % p,
                };
                let x = fe(v, p);
                assert_eq!(x.inv_vartime(), x.inv(), "p={p} v={v}");
            }
        }
        // A word above p (no constructor makes one) takes the `modinv`
        // path, which reduces it first.
        let unreduced = FieldElement {
            value: BigUint::from(10_010u32),
            modulus: BigUint::from(10_007u32),
        };
        assert_eq!(unreduced.inv_vartime(), unreduced.inv());
        // An unreduced nonzero multiple of p has no inverse; both return
        // the Fermat value 0.
        let multiple = FieldElement {
            value: BigUint::from(20_014u32),
            modulus: BigUint::from(10_007u32),
        };
        assert_eq!(multiple.inv_vartime(), multiple.inv());
        assert_eq!(
            multiple.inv_vartime().map(|v| v.value),
            Some(BigUint::zero())
        );
        let p256 = crate::ecc::curve::CurveParams::p256().p;
        let x = FieldElement::new(BigUint::from(0xdead_beefu32).pow(9), p256.clone());
        assert_eq!(x.inv_vartime(), x.inv());
        assert!(FieldElement::zero(p256).inv_vartime().is_none());
    }

    /// `sub_vartime` is `sub` on reduced operands, on both sides of the
    /// wrap and on a multi-word modulus, and on unreduced ones.
    #[test]
    fn sub_vartime_matches_sub() {
        for p in [2u64, 7, 10_007, u64::MAX - 58] {
            for a in [0, 1, p / 2, p - 1] {
                for b in [0, 1, p / 3, p - 1] {
                    let (x, y) = (fe(a, p), fe(b, p));
                    let (s, t) = (x.sub_vartime(&y), x.sub(&y));
                    assert_eq!((&s.value, &s.modulus), (&t.value, &t.modulus));
                }
            }
        }
        let p256 = crate::ecc::curve::CurveParams::p256().p;
        let a = FieldElement::new(BigUint::from(3u32).pow(150), p256.clone());
        let b = FieldElement::new(BigUint::from(5u32).pow(100), p256);
        assert_eq!(a.sub_vartime(&b).value, a.sub(&b).value);
        assert_eq!(b.sub_vartime(&a).value, b.sub(&a).value);
        // Unreduced operands (no constructor makes them) come out reduced,
        // as from `sub`: differences at and above p on both branches.
        let p = BigUint::from(10_007u32);
        for (a, b) in [
            (20_020u32, 3u32),
            (10_007, 0),
            (30_000, 29_999),
            (3, 5),
            (10_010, 10_012),
        ] {
            let x = FieldElement {
                value: BigUint::from(a),
                modulus: p.clone(),
            };
            let y = FieldElement {
                value: BigUint::from(b),
                modulus: p.clone(),
            };
            assert_eq!(x.sub_vartime(&y).value, x.sub(&y).value, "{a} - {b}");
        }
    }

    /// `neg_vartime` is `neg` on reduced values, zero included, on one
    /// word and on a multi-word modulus.
    #[test]
    fn neg_vartime_matches_neg() {
        for p in [2u64, 7, 10_007, u64::MAX - 58] {
            for a in [0, 1, p / 2, p - 1] {
                let (s, t) = (fe(a, p).neg_vartime(), fe(a, p).neg());
                assert_eq!((&s.value, &s.modulus), (&t.value, &t.modulus));
            }
        }
        let p256 = crate::ecc::curve::CurveParams::p256().p;
        for v in [BigUint::zero(), BigUint::from(3u32).pow(150)] {
            let x = FieldElement::new(v, p256.clone());
            assert_eq!(x.neg_vartime().value, x.neg().value);
        }
    }

    #[test]
    fn pow_identity() {
        // In F_7: 2^6 ≡ 1 (Fermat's little theorem)
        let r = fe(2, 7).pow(&BigUint::from(6u32));
        assert_eq!(r, fe(1, 7));
    }
}
