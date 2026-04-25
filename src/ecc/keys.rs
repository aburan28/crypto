//! ECC key generation and representation.
//!
//! A private key is a random scalar d in [1, n-1].
//! The corresponding public key is the point Q = d·G.

use super::curve::CurveParams;
use super::point::Point;
use crate::utils::random::random_scalar;
use num_bigint::BigUint;
/// An ECC private key: a random scalar d ∈ [1, n-1].
#[derive(Clone, Debug)]
pub struct EccPrivateKey {
    pub scalar: BigUint,
    pub curve_name: String,
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
