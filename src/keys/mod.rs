//! Key management: unified key types and serialization helpers.
//!
//! This module provides `KeyBundle` — a single owner's set of keys across
//! all algorithms supported by this library.

use crate::asymmetric::rsa::RsaKeyPair;
use crate::ecc::curve::CurveParams;
use crate::ecc::keys::{EccKeyPair, EccPrivateKey, EccPublicKey};
use crate::utils::encoding::{bigint_to_bytes_be, to_hex};

/// A complete bundle of keys for one identity: ECC + RSA.
pub struct KeyBundle {
    pub ecc: EccKeyPair,
    pub rsa: RsaKeyPair,
}

impl KeyBundle {
    /// Generate a new key bundle (secp256k1 ECC + 2048-bit RSA).
    ///
    /// Note: RSA key generation is slow (~seconds for 2048 bits).
    pub fn generate() -> Self {
        let curve = CurveParams::secp256k1();
        KeyBundle {
            ecc: EccKeyPair::generate(&curve),
            rsa: RsaKeyPair::generate(2048),
        }
    }

    /// Print a human-readable key summary.
    pub fn summary(&self) {
        println!("=== ECC (secp256k1) ===");
        println!("Private: {}", to_hex(&bigint_to_bytes_be(&self.ecc.private.scalar, 32)));
        if let Some(pub_bytes) = self.ecc.public.to_sec1_uncompressed() {
            println!("Public:  {}", to_hex(&pub_bytes));
        }
        println!();
        println!("=== RSA-{} ===", self.rsa.public.bits);
        let n_hex = to_hex(&bigint_to_bytes_be(&self.rsa.public.n, (self.rsa.public.bits as usize + 7) / 8));
        println!("n (first 32 bytes): {}", &n_hex[..64]);
        println!("e: {}", self.rsa.public.e);
    }
}

/// Serialize an ECC private key as a 32-byte hex string.
pub fn ecc_private_to_hex(k: &EccPrivateKey) -> String {
    to_hex(&bigint_to_bytes_be(&k.scalar, 32))
}

/// Serialize an ECC public key as uncompressed SEC1 hex.
pub fn ecc_public_to_hex(k: &EccPublicKey) -> Option<String> {
    k.to_sec1_uncompressed().map(|b| to_hex(&b))
}
