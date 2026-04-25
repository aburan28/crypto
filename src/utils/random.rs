//! Cryptographically secure random number generation utilities.
//!
//! Wraps `rand::thread_rng` (which uses the OS CSPRNG) and provides
//! helpers for generating scalars in a specific range.

use num_bigint::{BigUint, RandBigInt};
use num_traits::Zero;
use rand::RngCore;

/// Fill `buf` with cryptographically random bytes.
pub fn random_bytes(buf: &mut [u8]) {
    rand::thread_rng().fill_bytes(buf);
}

/// Generate a random `BigUint` in the range `[1, max)`.
///
/// Uses rejection sampling so the distribution is uniform — no modular bias.
pub fn random_scalar(max: &BigUint) -> BigUint {
    let mut rng = rand::thread_rng();
    loop {
        let candidate = rng.gen_biguint_below(max);
        if !candidate.is_zero() {
            return candidate;
        }
    }
}

/// Generate `len` random bytes as a `Vec<u8>`.
pub fn random_bytes_vec(len: usize) -> Vec<u8> {
    let mut buf = vec![0u8; len];
    rand::thread_rng().fill_bytes(&mut buf);
    buf
}
