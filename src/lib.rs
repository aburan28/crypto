//! # crypto-lib — a comprehensive educational cryptography library
//!
//! Implemented-from-scratch algorithms across all major cryptographic families:
//!
//! | Module         | Algorithms                                              |
//! |----------------|---------------------------------------------------------|
//! | `ecc`          | secp256k1/P-256/SM2/GOST-3410-2012 field ops, ECDSA, ECDH, SM2 sig/enc, GOST sig |
//! | `symmetric`    | AES, ChaCha20-Poly1305, Serpent, Threefish, SM4, Kuznyechik, GOST Magma |
//! | `asymmetric`   | RSA (keygen/enc/sig, CRT-CT), Paillier + ElGamal + EC-ElGamal HE |
//! | `hash`         | SHA-256, SHA-3 (Keccak), BLAKE3, SM3, Streebog          |
//! | `kdf`          | HKDF-SHA256, PBKDF2-HMAC-SHA256                         |
//! | `pqc`          | Simplified Kyber/ML-KEM and binary-Goppa McEliece       |
//! | `keys`         | Unified key management                                  |
//! | `cryptanalysis`| LLL/BKZ + HNP + multi-key-HNP + ρ-variants + Bleichenbacher + Smart + corpus-sweep |
//! | `ecc_safety`   | ECDH/ECDSA curve-parameter safety auditor               |
//! | `utils`        | Modular arithmetic, randomness, encoding                |
//! | `zk`           | Schnorr ZKP, Pedersen commitments, Chaum-Pedersen, Merkle tree |

pub mod asymmetric;
pub mod binary_ecc;
pub mod bls12_381;
pub mod cryptanalysis;
pub mod ct_bignum;
pub mod ecc;
pub mod ecc_safety;
pub mod encoding;
pub mod hash;
pub mod kdf;
pub mod keys;
pub mod pqc;
pub mod symmetric;
pub mod utils;
pub mod zk;
