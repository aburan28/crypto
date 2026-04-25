//! ECC key generation and representation.
//!
//! A private key is a random scalar d in [1, n-1].
//! The corresponding public key is the point Q = d·G.

use super::curve::CurveParams;
use super::point::Point;
use crate::utils::random::random_scalar;
use num_bigint::BigUint;
use num_traits::Zero;

/// An ECC private key: a random scalar d ∈ [1, n-1].
///
/// On drop, the scalar is best-effort overwritten with zero.  Note: the
/// underlying `BigUint` heap allocation is owned by `num-bigint` and the
/// crate does not expose its internal `Vec<u64>`; we can only re-assign
/// the high-level value to zero, which the allocator may or may not
/// scrub before reuse.  For genuinely defence-in-depth zeroization, swap
/// `num-bigint` for `crypto-bigint` (a constant-time, zeroize-aware
/// arbitrary-precision integer crate).
#[derive(Clone, Debug)]
pub struct EccPrivateKey {
    pub scalar: BigUint,
    pub curve_name: String,
}

impl Drop for EccPrivateKey {
    fn drop(&mut self) {
        self.scalar.set_zero();
        // String contents are public (curve name) — no need to scrub.
    }
}

/// An ECC public key: the point Q = d·G on the curve.
#[derive(Clone, Debug)]
pub struct EccPublicKey {
    pub point: Point,
    pub curve_name: String,
}

/// A private/public key pair.
pub struct EccKeyPair {
    pub private: EccPrivateKey,
    pub public: EccPublicKey,
}

impl EccKeyPair {
    /// Generate a fresh key pair on the given curve.
    pub fn generate(curve: &CurveParams) -> Self {
        let d = random_scalar(&curve.n);
        let q = curve.generator().scalar_mul(&d, &curve.a_fe());
        EccKeyPair {
            private: EccPrivateKey { scalar: d, curve_name: curve.name.to_string() },
            public: EccPublicKey { point: q, curve_name: curve.name.to_string() },
        }
    }

    /// Derive the public key from a known private scalar (useful for testing).
    pub fn from_private(d: BigUint, curve: &CurveParams) -> Self {
        let q = curve.generator().scalar_mul(&d, &curve.a_fe());
        EccKeyPair {
            private: EccPrivateKey { scalar: d, curve_name: curve.name.to_string() },
            public: EccPublicKey { point: q, curve_name: curve.name.to_string() },
        }
    }
}

impl EccPublicKey {
    /// Serialize the public key as uncompressed SEC1 encoding (04 || X || Y),
    /// where X and Y are 32-byte big-endian coordinates.
    pub fn to_sec1_uncompressed(&self) -> Option<Vec<u8>> {
        match &self.point {
            Point::Affine { x, y } => {
                let mut out = vec![0x04];
                let x_bytes = crate::utils::encoding::bigint_to_bytes_be(&x.value, 32);
                let y_bytes = crate::utils::encoding::bigint_to_bytes_be(&y.value, 32);
                out.extend_from_slice(&x_bytes);
                out.extend_from_slice(&y_bytes);
                Some(out)
            }
            Point::Infinity => None,
        }
    }
}
