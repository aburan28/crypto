//! Fixed-width arithmetic modulo the NIST P-521 prime, 2^521 - 1.
//!
//! Nine little-endian u64 limbs hold a canonical value. Multiplication uses
//! 81 word products and the Mersenne identity 2^521 = 1 (mod p).

const TOP: u64 = (1 << 9) - 1;
const P: [u64; 9] = [
    u64::MAX,
    u64::MAX,
    u64::MAX,
    u64::MAX,
    u64::MAX,
    u64::MAX,
    u64::MAX,
    u64::MAX,
    TOP,
];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct P521Fe([u64; 9]);

impl P521Fe {
    pub const ZERO: Self = Self([0; 9]);
    pub const ONE: Self = Self([1, 0, 0, 0, 0, 0, 0, 0, 0]);

    pub fn from_le_bytes(bytes: &[u8; 66]) -> Self {
        let mut h = [0u64; 9];
        for i in 0..8 {
            h[i] = u64::from_le_bytes(bytes[i * 8..i * 8 + 8].try_into().unwrap());
        }
        h[8] = u16::from_le_bytes([bytes[64], bytes[65]]) as u64;
        Self::normalize(h)
    }

    pub fn to_le_bytes(&self) -> [u8; 66] {
        let mut out = [0u8; 66];
        for i in 0..8 {
            out[i * 8..i * 8 + 8].copy_from_slice(&self.0[i].to_le_bytes());
        }
        out[64..66].copy_from_slice(&(self.0[8] as u16).to_le_bytes());
        out
    }

    #[inline]
    fn normalize(mut h: [u64; 9]) -> Self {
        // The caller's top limb is at most 16 bits. Two fixed end-around
        // passes absorb its excess and any carry through the low 512 bits.
        for _ in 0..2 {
            let wrap = h[8] >> 9;
            h[8] &= TOP;
            let (sum, mut carry) = h[0].overflowing_add(wrap);
            h[0] = sum;
            for limb in h.iter_mut().skip(1) {
                let (sum, next) = limb.overflowing_add(carry as u64);
                *limb = sum;
                carry = next;
            }
        }
        // The all-ones representation equals zero modulo p. Keep the
        // selection mask opaque to avoid a compiler branch on field data.
        let mut different = 0u64;
        for i in 0..9 {
            different |= h[i] ^ P[i];
        }
        let equal = core::hint::black_box((different == 0) as u64);
        let keep = 0u64.wrapping_sub(1 - equal);
        for limb in &mut h {
            *limb &= keep;
        }
        Self(h)
    }

    pub fn add(&self, rhs: &Self) -> Self {
        let mut h = [0u64; 9];
        let mut carry = 0u64;
        for i in 0..9 {
            let sum = self.0[i] as u128 + rhs.0[i] as u128 + carry as u128;
            h[i] = sum as u64;
            carry = (sum >> 64) as u64;
        }
        debug_assert_eq!(carry, 0);
        Self::normalize(h)
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        let mut h = [0u64; 9];
        let mut borrow = 0u64;
        for i in 0..9 {
            let (part, first) = self.0[i].overflowing_sub(rhs.0[i]);
            let (part, second) = part.overflowing_sub(borrow);
            h[i] = part;
            borrow = (first | second) as u64;
        }
        let mask = 0u64.wrapping_sub(borrow);
        let mut carry = 0u128;
        for i in 0..9 {
            let sum = h[i] as u128 + (P[i] & mask) as u128 + carry;
            h[i] = sum as u64;
            carry = sum >> 64;
        }
        Self::normalize(h)
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        let mut product = [0u64; 18];
        for i in 0..9 {
            let mut carry = 0u128;
            for j in 0..9 {
                let sum = self.0[i] as u128 * rhs.0[j] as u128 + product[i + j] as u128 + carry;
                product[i + j] = sum as u64;
                carry = sum >> 64;
            }
            product[i + 9] = carry as u64;
        }
        let mut h = [0u64; 9];
        let mut carry = 0u128;
        for i in 0..9 {
            let low = if i == 8 { product[i] & TOP } else { product[i] };
            let high = (product[i + 8] >> 9) | (product[i + 9] << 55);
            let sum = low as u128 + high as u128 + carry;
            h[i] = sum as u64;
            carry = sum >> 64;
        }
        debug_assert_eq!(carry, 0);
        Self::normalize(h)
    }

    pub fn square(&self) -> Self {
        self.mul(self)
    }

    #[inline]
    fn square_n(&self, count: usize) -> Self {
        let mut value = *self;
        for _ in 0..count {
            value = value.square();
        }
        value
    }

    /// Fermat inversion using a fixed chain for 2^521 - 3.
    /// Zero maps to zero; callers that require an inverse must reject zero.
    pub fn invert(&self) -> Self {
        // f[k] = a^(2^k - 1). Doubling: f[2k] = f[k]^(2^k) f[k].
        let f1 = *self;
        let f2 = f1.square().mul(&f1);
        let f4 = f2.square_n(2).mul(&f2);
        let f8 = f4.square_n(4).mul(&f4);
        let f16 = f8.square_n(8).mul(&f8);
        let f32 = f16.square_n(16).mul(&f16);
        let f64 = f32.square_n(32).mul(&f32);
        let f128 = f64.square_n(64).mul(&f64);
        let f256 = f128.square_n(128).mul(&f128);
        let f512 = f256.square_n(256).mul(&f256);
        let f3 = f2.square().mul(&f1);
        let f7 = f4.square_n(3).mul(&f3);
        let f519 = f512.square_n(7).mul(&f7);
        f519.square_n(2).mul(&f1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use num_traits::{One, Zero};

    #[test]
    fn arithmetic_matches_integer_reference() {
        let p: BigUint = (BigUint::one() << 521usize) - BigUint::one();
        let mut samples = vec![[0u8; 66], [0xff; 66], [0x55; 66], [0xaa; 66]];
        let mut one = [0u8; 66];
        one[0] = 1;
        samples.push(one);
        let mut state = 0x41b7_6d8e_579d_1531u64;
        for _ in 0..32 {
            let mut bytes = [0u8; 66];
            for chunk in bytes.chunks_mut(8) {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                chunk.copy_from_slice(&state.to_le_bytes()[..chunk.len()]);
            }
            samples.push(bytes);
        }
        for a_bytes in &samples {
            let a = P521Fe::from_le_bytes(a_bytes);
            let a_ref = BigUint::from_bytes_le(a_bytes) % &p;
            assert_eq!(BigUint::from_bytes_le(&a.to_le_bytes()), a_ref);
            assert_eq!(
                BigUint::from_bytes_le(&a.square().to_le_bytes()),
                (&a_ref * &a_ref) % &p
            );
            for b_bytes in samples.iter().take(9) {
                let b = P521Fe::from_le_bytes(b_bytes);
                let b_ref = BigUint::from_bytes_le(b_bytes) % &p;
                assert_eq!(
                    BigUint::from_bytes_le(&a.add(&b).to_le_bytes()),
                    (&a_ref + &b_ref) % &p
                );
                assert_eq!(
                    BigUint::from_bytes_le(&a.sub(&b).to_le_bytes()),
                    (&a_ref + &p - &b_ref) % &p
                );
                assert_eq!(
                    BigUint::from_bytes_le(&a.mul(&b).to_le_bytes()),
                    (&a_ref * &b_ref) % &p
                );
            }
        }
        for a_bytes in samples.iter().take(8) {
            let a = P521Fe::from_le_bytes(a_bytes);
            let a_ref = BigUint::from_bytes_le(a_bytes) % &p;
            let expected = if a_ref.is_zero() {
                BigUint::zero()
            } else {
                a_ref.modpow(&(&p - 2u8), &p)
            };
            assert_eq!(BigUint::from_bytes_le(&a.invert().to_le_bytes()), expected);
        }
    }
}
