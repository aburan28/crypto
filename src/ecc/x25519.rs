//! X25519 — Diffie-Hellman key agreement on Curve25519 (RFC 7748).
//!
//! The Montgomery-form curve is `B y² = x³ + A x² + x` with
//! `A = 486662, B = 1` over `F_{2²⁵⁵ − 19}`.  X25519 operates on
//! x-coordinates only via the **Montgomery ladder**, which is
//! implemented with fixed-limb field arithmetic and branchless swaps.
//!
//! # Wire format
//!
//! - **Public key**: 32 bytes, little-endian x-coordinate
//! - **Private key (scalar)**: 32 bytes, "clamped":
//!   - `byte[0] &= 248` (clear bottom 3 bits)
//!   - `byte[31] &= 127` (clear top bit)
//!   - `byte[31] |= 64` (set bit 254)
//!     This forces the scalar into the prime-order subgroup and
//!     prevents small-subgroup attacks.
//! - **Shared secret**: 32 bytes, little-endian x-coordinate of
//!   `[scalar] · public`
//!
//! # Test vectors
//!
//! Verified against RFC 7748 §5.2 vectors.

use super::f25519::Fe;

/// "Clamp" a 32-byte private scalar per RFC 7748 §5.
pub fn clamp(scalar: &mut [u8; 32]) {
    scalar[0] &= 248;
    scalar[31] &= 127;
    scalar[31] |= 64;
}

/// X25519 scalar multiplication: `[k] · u`.  Both inputs are
/// 32-byte little-endian; output is 32 bytes.
pub fn x25519(scalar: &[u8; 32], u_coord: &[u8; 32]) -> [u8; 32] {
    let mut k = *scalar;
    clamp(&mut k);
    // RFC 7748 §5: implementations MUST mask the most-significant bit
    // of the u-coordinate before reducing.
    let mut u_masked = *u_coord;
    u_masked[31] &= 0x7F;
    ladder(&k, &Fe::from_bytes(&u_masked)).to_bytes()
}

/// X25519 key agreement with RFC 7748's all-zero shared-secret check.
///
/// Low-order public inputs can force the Montgomery ladder to return the
/// all-zero string. Callers doing ECDH should reject that result rather than
/// treating it as a successful shared secret.
pub fn x25519_checked(scalar: &[u8; 32], u_coord: &[u8; 32]) -> Option<[u8; 32]> {
    let shared = x25519(scalar, u_coord);
    if shared.iter().all(|&b| b == 0) {
        None
    } else {
        Some(shared)
    }
}

/// X25519 base-point scalar multiplication: derive a public key
/// from a clamped 32-byte scalar.  Base point `u = 9`.
pub fn x25519_base(scalar: &[u8; 32]) -> [u8; 32] {
    let mut base = [0u8; 32];
    base[0] = 9;
    x25519(scalar, &base)
}

/// RFC 7748 Montgomery ladder on five radix-2^51 limbs.
fn ladder(scalar: &[u8; 32], u: &Fe) -> Fe {
    let x1 = *u;
    let mut x2 = Fe::ONE;
    let mut z2 = Fe::ZERO;
    let mut x3 = *u;
    let mut z3 = Fe::ONE;
    let mut swap = 0u8;

    // Iterate from MSB to LSB (bit 254 down to 0).
    for t in (0..255).rev() {
        let k_t = (scalar[t / 8] >> (t % 8)) & 1;
        let s = swap ^ k_t;
        Fe::cswap(&mut x2, &mut x3, s);
        Fe::cswap(&mut z2, &mut z3, s);
        swap = k_t;

        let a = x2.add(&z2);
        let aa = a.square();
        let b = x2.sub(&z2);
        let bb = b.square();
        let e = aa.sub(&bb);
        let c = x3.add(&z3);
        let d = x3.sub(&z3);
        let da = d.mul(&a);
        let cb = c.mul(&b);

        x3 = da.add(&cb).square();
        z3 = x1.mul(&da.sub(&cb).square());

        x2 = aa.mul(&bb);
        // RFC 7748: z_2 = E · (AA + a24 · E)
        z2 = e.mul(&aa.add(&e.mul_a24()));
    }
    Fe::cswap(&mut x2, &mut x3, swap);
    Fe::cswap(&mut z2, &mut z3, swap);
    x2.mul(&z2.invert())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::curve25519::{fe_add, fe_from_bytes, fe_inv, fe_mul, fe_sub, fe_to_bytes};
    use num_bigint::BigUint;

    fn hex_decode(s: &str) -> Vec<u8> {
        let s: String = s.chars().filter(|c| !c.is_whitespace()).collect();
        (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
            .collect()
    }

    fn hex32(s: &str) -> [u8; 32] {
        let v = hex_decode(s);
        let mut out = [0u8; 32];
        out.copy_from_slice(&v);
        out
    }

    fn hex(bytes: &[u8]) -> String {
        bytes.iter().map(|b| format!("{:02x}", b)).collect()
    }

    /// Preserve the original BigUint ladder as an independent output oracle
    /// for the replacement field and scalar path.
    fn reference_x25519(scalar: &[u8; 32], u_coord: &[u8; 32]) -> [u8; 32] {
        let mut k = *scalar;
        clamp(&mut k);
        let scalar = BigUint::from_bytes_le(&k);
        let mut u_masked = *u_coord;
        u_masked[31] &= 0x7f;
        let u = fe_from_bytes(&u_masked);
        let mut x2 = BigUint::from(1u32);
        let mut z2 = BigUint::from(0u32);
        let mut x3 = u.clone();
        let mut z3 = BigUint::from(1u32);
        let mut swap = 0u32;
        let a24 = BigUint::from(121665u32);
        for t in (0..255).rev() {
            let bit = scalar.bit(t) as u32;
            if swap ^ bit == 1 {
                std::mem::swap(&mut x2, &mut x3);
                std::mem::swap(&mut z2, &mut z3);
            }
            swap = bit;
            let a = fe_add(&x2, &z2);
            let aa = fe_mul(&a, &a);
            let b = fe_sub(&x2, &z2);
            let bb = fe_mul(&b, &b);
            let e = fe_sub(&aa, &bb);
            let c = fe_add(&x3, &z3);
            let d = fe_sub(&x3, &z3);
            let da = fe_mul(&d, &a);
            let cb = fe_mul(&c, &b);
            x3 = fe_mul(&fe_add(&da, &cb), &fe_add(&da, &cb));
            z3 = fe_mul(&u, &fe_mul(&fe_sub(&da, &cb), &fe_sub(&da, &cb)));
            x2 = fe_mul(&aa, &bb);
            z2 = fe_mul(&e, &fe_add(&aa, &fe_mul(&a24, &e)));
        }
        if swap == 1 {
            std::mem::swap(&mut x2, &mut x3);
            std::mem::swap(&mut z2, &mut z3);
        }
        fe_to_bytes(&fe_mul(&x2, &fe_inv(&z2)))
    }

    #[test]
    fn fixed_limb_ladder_matches_original_reference() {
        let mut cases = vec![[0u8; 32], [0xff; 32], [0x55; 32], [0xaa; 32]];
        let mut modulus = [0xffu8; 32];
        modulus[0] = 0xed;
        modulus[31] = 0x7f;
        cases.push(modulus);
        let mut state = 0x91e1_0da5_c79e_7b1du64;
        for _ in 0..16 {
            let mut bytes = [0u8; 32];
            for chunk in bytes.chunks_exact_mut(8) {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                chunk.copy_from_slice(&state.to_le_bytes());
            }
            cases.push(bytes);
        }
        for (i, scalar) in cases.iter().enumerate() {
            for u in [&cases[i], &cases[(i + 7) % cases.len()]] {
                assert_eq!(
                    x25519(scalar, u),
                    reference_x25519(scalar, u),
                    "scalar case {i}, u case mismatch"
                );
            }
        }
    }

    /// **RFC 7748 §5.2 test vector 1**.  Scalar and u given, expected output.
    #[test]
    fn rfc7748_vector_1() {
        let scalar = hex32("a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4");
        let u = hex32("e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c");
        let expected = "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552";
        let result = x25519(&scalar, &u);
        assert_eq!(hex(&result), expected, "RFC 7748 vector 1 mismatch");
    }

    /// **RFC 7748 §5.2 test vector 2**.
    #[test]
    fn rfc7748_vector_2() {
        let scalar = hex32("4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d");
        let u = hex32("e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493");
        let expected = "95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957";
        let result = x25519(&scalar, &u);
        assert_eq!(hex(&result), expected, "RFC 7748 vector 2 mismatch");
    }

    /// X25519 ECDH round-trip: `[a]([b]G) = [b]([a]G)`.
    #[test]
    fn ecdh_roundtrip() {
        let a_priv = [0x01u8; 32];
        let b_priv = [0x02u8; 32];
        let a_pub = x25519_base(&a_priv);
        let b_pub = x25519_base(&b_priv);
        let s_ab = x25519(&a_priv, &b_pub);
        let s_ba = x25519(&b_priv, &a_pub);
        assert_eq!(s_ab, s_ba, "ECDH shared-secret asymmetry — bug");
    }

    #[test]
    fn checked_x25519_rejects_low_order_input() {
        let scalar = [0x11u8; 32];
        let low_order = [0u8; 32];
        assert!(
            x25519_checked(&scalar, &low_order).is_none(),
            "low-order public inputs must not be accepted as ECDH secrets"
        );
        assert_eq!(
            x25519(&scalar, &low_order),
            [0u8; 32],
            "regression guard: raw primitive still exposes the RFC all-zero output"
        );
    }

    /// Clamping zeros the right bits.
    #[test]
    fn clamp_zeros_correct_bits() {
        let mut s = [0xFFu8; 32];
        clamp(&mut s);
        assert_eq!(s[0] & 0x07, 0, "low 3 bits not cleared");
        assert_eq!(s[31] & 0x80, 0, "top bit not cleared");
        assert_eq!(s[31] & 0x40, 0x40, "bit 254 not set");
    }
}
