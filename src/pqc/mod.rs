//! Post-quantum cryptography: simplified Kyber/ML-KEM key encapsulation.
//!
//! See individual modules for algorithm documentation and limitations.

pub mod kyber;

pub use kyber::{
    KyberCiphertext, KyberPrivateKey, KyberPublicKey,
    kyber_decapsulate, kyber_encapsulate, kyber_keygen,
};
