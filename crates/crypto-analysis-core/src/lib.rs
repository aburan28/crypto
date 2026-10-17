//! Small, reusable analysis surface from the crypto repository.
//!
//! The implementations live in `src/cryptanalysis/` so the main `crypto_lib`
//! crate and this path dependency compile the same source files.

#[path = "../../../src/cryptanalysis/avalanche.rs"]
pub mod avalanche;

#[path = "../../../src/cryptanalysis/statistical.rs"]
pub mod statistical;

#[cfg(test)]
mod hash {
    use sha2::{Digest, Sha256};

    pub fn sha256(data: &[u8]) -> [u8; 32] {
        Sha256::digest(data).into()
    }
}
