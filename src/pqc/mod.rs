//! Post-quantum cryptography: simplified Kyber/ML-KEM key encapsulation and
//! a code-based McEliece cryptosystem.
//!
//! See individual modules for algorithm documentation and limitations.

pub mod kyber;
pub mod mceliece;

pub use kyber::{
    KyberCiphertext, KyberPrivateKey, KyberPublicKey,
    kyber_decapsulate, kyber_encapsulate, kyber_keygen,
};

pub use mceliece::{
    McElieceKeyPair, McEliecePrivateKey, McElieceePublicKey,
    mceliece_decrypt, mceliece_encrypt,
};
