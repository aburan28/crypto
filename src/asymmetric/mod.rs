//! Asymmetric cryptography: RSA key generation, encryption, and signatures.

pub mod rsa;

pub use rsa::{
    RsaKeyPair, RsaPrivateKey, RsaPublicKey,
    is_prime, random_prime,
    rsa_decrypt, rsa_encrypt,
    rsa_sign, rsa_verify,
};
// Note: textbook RSA (`rsa_encrypt_raw`/`rsa_decrypt_raw`) and the PKCS#1 v1.5
// padding helpers are deliberately not re-exported.  They are `pub(crate)`
// internal helpers for `rsa_encrypt`/`rsa_decrypt`; exposing them invites
// misuse (textbook RSA is malleable and not semantically secure).
