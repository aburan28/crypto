//! Five-limb radix-2^51 arithmetic for the X25519 Montgomery ladder.
//!
//! The limbs are always reduced below 2^51 after each operation. This module
//! intentionally keeps the public RFC 7748 input decoder separate from the
//! strict Ed25519 point decoder.

const MASK: u64 = (1u64 << 51) - 1;
const P: [u64; 5] = [MASK - 18, MASK, MASK, MASK, MASK];

#[derive(Clone, Copy, Debug)]
pub(super) struct Fe(pub(super) [u64; 5]);

impl Fe {
    pub const ZERO: Self = Self([0; 5]);
    pub const ONE: Self = Self([1, 0, 0, 0, 0]);

    /// RFC 7748 masks bit 255 and accepts noncanonical coordinates.
    pub fn from_bytes(input: &[u8; 32]) -> Self {
        let mut bytes = *input;
        bytes[31] &= 0x7f;
        let load = |offset: usize| -> u64 {
            u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
        };
        Self([
            load(0) & MASK,
            (load(6) >> 3) & MASK,
            (load(12) >> 6) & MASK,
            (load(19) >> 1) & MASK,
            (load(24) >> 12) & MASK,
        ])
    }

    pub fn to_bytes(self) -> [u8; 32] {
        let h = self.0;
        // Select h-p iff h >= p. All limbs are below 2^51 here.
        let mut difference = [0u64; 5];
        let mut borrow = 0u64;
        for i in 0..5 {
            let (part, first) = h[i].overflowing_sub(P[i]);
            let (part, second) = part.overflowing_sub(borrow);
            difference[i] = part;
            borrow = (first as u64) | (second as u64);
        }
        // Keep the all-zero/all-one mask opaque to LLVM. Without this
        // barrier it turned the canonical selection into data-dependent
        // branches on ARM64 even though the Rust expression used bitmasks.
        let mask = core::hint::black_box(0u64.wrapping_sub(1 - borrow));
        let mut canonical = [0u64; 5];
        for i in 0..5 {
            canonical[i] = (h[i] & !mask) | (difference[i] & mask);
        }
        let words = [
            canonical[0] | (canonical[1] << 51),
            (canonical[1] >> 13) | (canonical[2] << 38),
            (canonical[2] >> 26) | (canonical[3] << 25),
            (canonical[3] >> 39) | (canonical[4] << 12),
        ];
        let mut result = [0u8; 32];
        for (i, word) in words.iter().enumerate() {
            result[i * 8..i * 8 + 8].copy_from_slice(&word.to_le_bytes());
        }
        result
    }

    #[inline]
    fn reduce(mut h: [u128; 5]) -> Self {
        // Three fixed carry passes cover both a full multiplication and the
        // final wrap from bit 255 back into limb zero.
        for _ in 0..3 {
            for i in 0..4 {
                let carry = h[i] >> 51;
                h[i] &= MASK as u128;
                h[i + 1] += carry;
            }
            let carry = h[4] >> 51;
            h[4] &= MASK as u128;
            h[0] += 19 * carry;
        }
        debug_assert!(h.iter().all(|&limb| limb <= MASK as u128));
        Self(h.map(|limb| limb as u64))
    }

    #[inline]
    pub fn add(&self, b: &Self) -> Self {
        Self::reduce(core::array::from_fn(|i| (self.0[i] + b.0[i]) as u128))
    }

    #[inline]
    pub fn sub(&self, b: &Self) -> Self {
        Self::reduce(core::array::from_fn(|i| {
            (self.0[i] as u128) + 2 * (P[i] as u128) - (b.0[i] as u128)
        }))
    }

    #[inline]
    pub fn mul(&self, b: &Self) -> Self {
        let a = self.0.map(u128::from);
        let b = b.0.map(u128::from);
        let h0 = a[0] * b[0] + 19 * (a[1] * b[4] + a[2] * b[3] + a[3] * b[2] + a[4] * b[1]);
        let h1 = a[0] * b[1] + a[1] * b[0] + 19 * (a[2] * b[4] + a[3] * b[3] + a[4] * b[2]);
        let h2 = a[0] * b[2] + a[1] * b[1] + a[2] * b[0] + 19 * (a[3] * b[4] + a[4] * b[3]);
        let h3 = a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0] + 19 * a[4] * b[4];
        let h4 = a[0] * b[4] + a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + a[4] * b[0];
        Self::reduce([h0, h1, h2, h3, h4])
    }

    #[inline]
    pub fn square(&self) -> Self {
        let a = self.0.map(u128::from);
        let h0 = a[0] * a[0] + 38 * (a[1] * a[4] + a[2] * a[3]);
        let h1 = 2 * a[0] * a[1] + 38 * a[2] * a[4] + 19 * a[3] * a[3];
        let h2 = 2 * a[0] * a[2] + a[1] * a[1] + 38 * a[3] * a[4];
        let h3 = 2 * (a[0] * a[3] + a[1] * a[2]) + 19 * a[4] * a[4];
        let h4 = 2 * (a[0] * a[4] + a[1] * a[3]) + a[2] * a[2];
        Self::reduce([h0, h1, h2, h3, h4])
    }

    #[inline]
    pub fn mul_a24(&self) -> Self {
        Self::reduce(self.0.map(|limb| (limb as u128) * 121665))
    }

    #[inline]
    fn square_n(&self, count: usize) -> Self {
        let mut result = *self;
        for _ in 0..count {
            result = result.square();
        }
        result
    }

    pub fn invert(&self) -> Self {
        // Fixed addition chain for p-2 = 2^255-21: 254 squares, 11
        // multiplies, and no scalar-dependent decisions.
        let z2 = self.square();
        let z9 = z2.square_n(2).mul(self);
        let z11 = z9.mul(&z2);
        let z2_5_0 = z11.square().mul(&z9);
        let z2_10_0 = z2_5_0.square_n(5).mul(&z2_5_0);
        let z2_20_0 = z2_10_0.square_n(10).mul(&z2_10_0);
        let z2_40_0 = z2_20_0.square_n(20).mul(&z2_20_0);
        let z2_50_0 = z2_40_0.square_n(10).mul(&z2_10_0);
        let z2_100_0 = z2_50_0.square_n(50).mul(&z2_50_0);
        let z2_200_0 = z2_100_0.square_n(100).mul(&z2_100_0);
        let z2_250_0 = z2_200_0.square_n(50).mul(&z2_50_0);
        z2_250_0.square_n(5).mul(&z11)
    }

    /// Conditional swap without a scalar-dependent branch.
    pub fn cswap(a: &mut Self, b: &mut Self, bit: u8) {
        let mask = 0u64.wrapping_sub((bit & 1) as u64);
        for i in 0..5 {
            let selected = mask & (a.0[i] ^ b.0[i]);
            a.0[i] ^= selected;
            b.0[i] ^= selected;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use num_traits::One;

    fn reference(fe: Fe) -> BigUint {
        BigUint::from_bytes_le(&fe.to_bytes())
    }

    #[test]
    fn arithmetic_matches_biguint_at_boundaries_and_seeded_inputs() {
        let p = (BigUint::one() << 255) - BigUint::from(19u8);
        let mut samples = vec![[0u8; 32], [0xff; 32], [0x55; 32], [0xaa; 32]];
        let mut one = [0u8; 32];
        one[0] = 1;
        samples.push(one);
        let mut state = 0x0123_4567_89ab_cdefu64;
        for _ in 0..64 {
            let mut bytes = [0u8; 32];
            for chunk in bytes.chunks_exact_mut(8) {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                chunk.copy_from_slice(&state.to_le_bytes());
            }
            samples.push(bytes);
        }
        for a_bytes in &samples {
            let a = Fe::from_bytes(a_bytes);
            let mut masked_a = *a_bytes;
            masked_a[31] &= 0x7f;
            let a_ref = BigUint::from_bytes_le(&masked_a) % &p;
            assert_eq!(reference(a), a_ref);
            assert_eq!(reference(a.square()), (&a_ref * &a_ref) % &p);
            assert_eq!(
                reference(a.mul_a24()),
                (&a_ref * BigUint::from(121665u32)) % &p
            );
            for b_bytes in samples.iter().take(9) {
                let b = Fe::from_bytes(b_bytes);
                let mut masked_b = *b_bytes;
                masked_b[31] &= 0x7f;
                let b_ref = BigUint::from_bytes_le(&masked_b) % &p;
                assert_eq!(reference(a.add(&b)), (&a_ref + &b_ref) % &p);
                assert_eq!(reference(a.sub(&b)), (&a_ref + &p - &b_ref) % &p);
                assert_eq!(reference(a.mul(&b)), (&a_ref * &b_ref) % &p);
            }
        }
        for a_bytes in samples.iter().take(12).skip(1) {
            let a = Fe::from_bytes(a_bytes);
            let expected = reference(a).modpow(&(&p - BigUint::from(2u8)), &p);
            assert_eq!(reference(a.invert()), expected);
        }
    }
}
