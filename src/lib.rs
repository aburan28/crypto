//! # crypto-lib — a comprehensive educational cryptography library
//!
//! Implemented-from-scratch algorithms across all major cryptographic families:
//!
//! | Module         | Algorithms                                              |
//! |----------------|---------------------------------------------------------|
//! | `ecc`          | secp256k1/P-256 field ops, ECDSA, ECDH                  |
//! | `symmetric`    | AES-128/256 (ECB/CTR/GCM), ChaCha20-Poly1305            |
//! | `asymmetric`   | RSA (keygen, encrypt, sign/verify)                      |
//! | `hash`         | SHA-256, SHA-3 (Keccak), BLAKE3                         |
//! | `kdf`          | HKDF-SHA256, PBKDF2-HMAC-SHA256                         |
//! | `pqc`          | Simplified Kyber/ML-KEM (educational)                   |
//! | `keys`         | Unified key management                                  |
//! | `utils`        | Modular arithmetic, randomness, encoding                |

pub mod asymmetric;
pub mod ecc;
pub mod hash;
pub mod kdf;
pub mod keys;
pub mod pqc;
pub mod symmetric;
pub mod utils;
