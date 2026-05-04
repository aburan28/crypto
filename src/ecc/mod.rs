//! Elliptic curve cryptography: field arithmetic, point operations,
//! key generation, ECDSA, and ECDH over secp256k1 / P-256.

pub mod ct;
pub mod curve;
pub mod ecdh;
pub mod ecdsa;
pub mod field;
pub mod keys;
pub mod point;
pub mod p256_field;
pub mod p256_point;
pub mod schnorr;
pub mod secp256k1_field;
pub mod secp256k1_point;

pub use curve::CurveParams;
pub use ecdh::{ecdh, ecdh_raw};
pub use ecdsa::{sign, sign_hash, verify, verify_hash, EcdsaSignature};
pub use field::FieldElement;
pub use keys::{EccKeyPair, EccPrivateKey, EccPublicKey};
pub use point::Point;
pub use schnorr::{
    schnorr_keypair, schnorr_sign, schnorr_verify, tagged_hash, xonly_pubkey,
};
