//! # crypto-lib — a comprehensive educational cryptography library
//!
//! Implemented-from-scratch algorithms across all major cryptographic families:
//!
//! | Module         | Algorithms                                              |
//! |----------------|---------------------------------------------------------|
//! | `ecc`          | secp256k1/P-256 field ops, ECDSA, ECDH                  |
//! | `symmetric`    | AES, ChaCha20-Poly1305, Serpent, Threefish-256/512/1024 |
//! | `asymmetric`   | RSA (keygen/enc/sig, CRT-CT), Paillier + ElGamal + EC-ElGamal HE |
//! | `hash`         | SHA-256, SHA-3 (Keccak), BLAKE3                         |
//! | `kdf`          | HKDF-SHA256, PBKDF2-HMAC-SHA256                         |
//! | `pqc`          | Simplified Kyber/ML-KEM and binary-Goppa McEliece       |
//! | `keys`         | Unified key management                                  |
//! | `cryptanalysis`| LLL/BKZ + HNP + multi-key-HNP + ρ-variants + Bleichenbacher + Smart + corpus-sweep |
//! | `ecc_safety`   | ECDH/ECDSA curve-parameter safety auditor               |
//! | `utils`        | Modular arithmetic, randomness, encoding                |

pub mod asymmetric;
pub mod cryptanalysis;
pub mod ct_bignum;
pub mod ecc;
pub mod ecc_safety;
pub mod hash;
pub mod kdf;
pub mod keys;
pub mod pqc;
pub mod symmetric;
pub mod utils;
