//! ECDSA — Elliptic Curve Digital Signature Algorithm (FIPS 186-4).
//!
//! # Sign
//! Given a private key d, curve order n, and message hash z:
//!   1. Choose a random nonce k ∈ [1, n-1].
//!   2. Compute (x₁, _) = k·G.
//!   3. r = x₁ mod n  (restart if r = 0).
//!   4. s = k⁻¹(z + r·d) mod n  (restart if s = 0).
//!   Output: (r, s).
//!
//! # Verify
//! Given public key Q, signature (r, s), and message hash z:
//!   1. Compute w = s⁻¹ mod n.
//!   2. u₁ = z·w mod n,  u₂ = r·w mod n.
//!   3. (x₁, _) = u₁·G + u₂·Q.
//!   4. Valid iff x₁ mod n = r.

use super::curve::CurveParams;
use super::keys::{EccPrivateKey, EccPublicKey};
use super::point::Point;
use crate::hash::sha256::sha256;
use crate::utils::random::random_scalar;
use crate::utils::mod_inverse;
use num_bigint::BigUint;
use num_traits::Zero;
/// An ECDSA signature pair (r, s).
#[derive(Clone, Debug)]
pub struct EcdsaSignature {
    pub r: BigUint,
    pub s: BigUint,
}

/// Sign `message` with `private_key` on `curve`.
///
/// The message is SHA-256 hashed before signing.
pub fn sign(message: &[u8], private_key: &EccPrivateKey, curve: &CurveParams) -> EcdsaSignature {
    sign_hash(&sha256(message), private_key, curve)
}

/// Sign a pre-computed 32-byte hash with `private_key`.
pub fn sign_hash(
    hash: &[u8],
    private_key: &EccPrivateKey,
    curve: &CurveParams,
) -> EcdsaSignature {
    let z = hash_to_scalar(hash, &curve.n);
    let a = curve.a_fe();

    loop {
        let k = random_scalar(&curve.n);
        let kg = curve.generator().scalar_mul(&k, &a);

        let x1 = match &kg {
            Point::Affine { x, .. } => x.value.clone(),
            Point::Infinity => continue,
        };

        let r = &x1 % &curve.n;
        if r.is_zero() {
            continue;
        }

        // s = k⁻¹ · (z + r·d) mod n
        let rd = (&r * &private_key.scalar) % &curve.n;
        let z_plus_rd = (&z + &rd) % &curve.n;
        let k_inv = match mod_inverse(&k, &curve.n) {
            Some(v) => v,
            None => continue,
        };
        let s = (&k_inv * &z_plus_rd) % &curve.n;
        if s.is_zero() {
            continue;
        }

        return EcdsaSignature { r, s };
    }
}

/// Verify that `sig` is a valid ECDSA signature over `message` for `public_key`.
pub fn verify(
    message: &[u8],
    public_key: &EccPublicKey,
    sig: &EcdsaSignature,
    curve: &CurveParams,
) -> bool {
    verify_hash(&sha256(message), public_key, sig, curve)
}

/// Verify a signature against a pre-computed hash.
pub fn verify_hash(
    hash: &[u8],
    public_key: &EccPublicKey,
    sig: &EcdsaSignature,
    curve: &CurveParams,
) -> bool {
    // Reject out-of-range values immediately.
    if sig.r.is_zero() || sig.r >= curve.n || sig.s.is_zero() || sig.s >= curve.n {
        return false;
    }

    let z = hash_to_scalar(hash, &curve.n);
    let a = curve.a_fe();

    let w = match mod_inverse(&sig.s, &curve.n) {
        Some(v) => v,
        None => return false,
    };

    let u1 = (&z * &w) % &curve.n;
    let u2 = (&sig.r * &w) % &curve.n;

    let point = curve
        .generator()
        .scalar_mul(&u1, &a)
        .add(&public_key.point.scalar_mul(&u2, &a), &a);

    match point {
        Point::Affine { x, .. } => (&x.value % &curve.n) == sig.r,
        Point::Infinity => false,
    }
}

/// Reduce a hash to a scalar in [0, n) by interpreting it as a big-endian integer.
fn hash_to_scalar(hash: &[u8], n: &BigUint) -> BigUint {
    BigUint::from_bytes_be(hash) % n
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::{curve::CurveParams, keys::EccKeyPair};

    #[test]
    fn sign_and_verify() {
        let curve = CurveParams::secp256k1();
        let kp = EccKeyPair::generate(&curve);
        let msg = b"hello, elliptic curves!";
        let sig = sign(msg, &kp.private, &curve);
        assert!(verify(msg, &kp.public, &sig, &curve));
    }

    #[test]
    fn wrong_message_fails() {
        let curve = CurveParams::secp256k1();
        let kp = EccKeyPair::generate(&curve);
        let sig = sign(b"original", &kp.private, &curve);
        assert!(!verify(b"tampered", &kp.public, &sig, &curve));
    }

    #[test]
    fn wrong_key_fails() {
        let curve = CurveParams::secp256k1();
        let kp1 = EccKeyPair::generate(&curve);
        let kp2 = EccKeyPair::generate(&curve);
        let sig = sign(b"message", &kp1.private, &curve);
        assert!(!verify(b"message", &kp2.public, &sig, &curve));
    }
}
