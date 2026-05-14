//! Hash functions.
//!
//! Modern from-scratch family: SHA-256, SHA-3/Keccak, SHA-512,
//! RIPEMD-160, BLAKE3, SipHash, SM3, Streebog.
//!
//! Plus the **classic broken** Merkle–Damgård hashes (MD4, MD5)
//! shipped as study targets for the cryptanalysis side, with the
//! standard "do not use for new designs" disclaimers in their module
//! docstrings.

pub mod blake3;
pub mod md4;
pub mod md5;
pub mod ripemd160;
pub mod sha256;
pub mod sha3;
pub mod sha512;
pub mod siphash;
pub mod sm3;
pub mod streebog;

pub use blake3::{blake3_derive_key, blake3_hash, blake3_keyed};
pub use md4::{md4, md4_compress, md4_pad};
pub use md5::{md5, md5_pad};
pub use ripemd160::{hash160, ripemd160};
pub use sha256::{sha224, sha256};
pub use sha3::{sha3_256, sha3_512, shake128, shake256};
pub use sha512::{sha384, sha512};
pub use siphash::{siphash, SipKey};
pub use sm3::{sm3, Sm3};
pub use streebog::{streebog_256, streebog_512, Streebog};
