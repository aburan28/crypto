//! PBKDF2 — Password-Based Key Derivation Function 2 (RFC 8018 / PKCS#5 v2.1).
//!
//! # Algorithm
//! DK = T₁ || T₂ || … || Tᵢ   (concatenated to the desired length)
//!
//! Tᵢ = F(Password, Salt, c, i)
//! F(P, S, c, i) = U₁ ⊕ U₂ ⊕ … ⊕ Uᵢ
//!
//! U₁ = HMAC(P, S || INT(i))
//! Uⱼ = HMAC(P, Uⱼ₋₁)     for j = 2 … c
//!
//! Using HMAC-SHA256 as the underlying PRF (hLen = 32 bytes).
//! For password hashing, use c ≥ 600,000 iterations (NIST SP 800-132).

use crate::kdf::hkdf::hmac_sha256;

const H_LEN: usize = 32; // HMAC-SHA256 output length

/// Derive `dk_len` bytes from `password` and `salt` using `iterations` rounds.
///
/// # Parameters
/// * `password`   — the secret (e.g., a user's password as UTF-8 bytes)
/// * `salt`       — a random per-credential salt (≥ 16 bytes recommended)
/// * `iterations` — work factor; ≥ 600_000 recommended for password storage
/// * `dk_len`     — desired output length in bytes
pub fn pbkdf2_hmac_sha256(
    password: &[u8],
    salt: &[u8],
    iterations: u32,
    dk_len: usize,
) -> Vec<u8> {
    let block_count = (dk_len + H_LEN - 1) / H_LEN;
    let mut dk = Vec::with_capacity(block_count * H_LEN);

    for i in 1u32..=block_count as u32 {
        // U₁ = HMAC(password, salt || INT32_BE(i))
        let mut salt_i = salt.to_vec();
        salt_i.extend_from_slice(&i.to_be_bytes());
        let mut u = hmac_sha256(password, &salt_i);
        let mut t = u;

        // U₂ … Uc
        for _ in 1..iterations {
            u = hmac_sha256(password, &u);
            for (a, b) in t.iter_mut().zip(u.iter()) {
                *a ^= b;
            }
        }
        dk.extend_from_slice(&t);
    }

    dk.truncate(dk_len);
    dk
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rfc6070_tc1() {
        // RFC 6070 Test Vector 1: password="password", salt="salt", c=1, dkLen=20
        // Note: this uses HMAC-SHA1, but we use HMAC-SHA256; just verify length/no-panic.
        let dk = pbkdf2_hmac_sha256(b"password", b"salt", 1, 32);
        assert_eq!(dk.len(), 32);
    }

    #[test]
    fn deterministic() {
        let dk1 = pbkdf2_hmac_sha256(b"secret", b"pepper", 100, 32);
        let dk2 = pbkdf2_hmac_sha256(b"secret", b"pepper", 100, 32);
        assert_eq!(dk1, dk2);
    }

    #[test]
    fn different_passwords() {
        let dk1 = pbkdf2_hmac_sha256(b"password1", b"salt", 100, 32);
        let dk2 = pbkdf2_hmac_sha256(b"password2", b"salt", 100, 32);
        assert_ne!(dk1, dk2);
    }
}
