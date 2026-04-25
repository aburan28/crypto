//! HKDF — HMAC-based Extract-and-Expand Key Derivation Function (RFC 5869).
//!
//! # Two-step process
//!
//! ## Extract
//! PRK = HMAC-SHA256(salt, IKM)
//! Converts the (possibly non-uniform) input keying material (IKM) into a
//! pseudorandom key PRK using `salt` as the HMAC key.  If no salt is provided,
//! a zero-filled string of HashLen bytes is used.
//!
//! ## Expand
//! OKM = T(1) || T(2) || … || T(ceil(L/HashLen))
//! where T(i) = HMAC-SHA256(PRK, T(i-1) || info || i)
//! Produces output keying material of the desired length from the PRK.

use crate::hash::sha256::sha256;

// ── HMAC-SHA256 ───────────────────────────────────────────────────────────────

const BLOCK_SIZE: usize = 64; // SHA-256 block size in bytes
const HASH_LEN: usize = 32;   // SHA-256 output length in bytes

/// HMAC-SHA256 as defined in RFC 2104.
pub fn hmac_sha256(key: &[u8], data: &[u8]) -> [u8; 32] {
    // Normalize the key to exactly BLOCK_SIZE bytes
    let mut k = [0u8; BLOCK_SIZE];
    if key.len() > BLOCK_SIZE {
        let hk = sha256(key);
        k[..HASH_LEN].copy_from_slice(&hk);
    } else {
        k[..key.len()].copy_from_slice(key);
    }

    // Inner hash: SHA256((k ⊕ ipad) || data)
    let mut ipad_data = [0x36u8; BLOCK_SIZE].to_vec();
    for (a, b) in ipad_data.iter_mut().zip(k.iter()) { *a ^= b; }
    ipad_data.extend_from_slice(data);
    let inner = sha256(&ipad_data);

    // Outer hash: SHA256((k ⊕ opad) || inner)
    let mut opad_data = [0x5cu8; BLOCK_SIZE].to_vec();
    for (a, b) in opad_data.iter_mut().zip(k.iter()) { *a ^= b; }
    opad_data.extend_from_slice(&inner);
    sha256(&opad_data)
}

// ── HKDF ─────────────────────────────────────────────────────────────────────

/// HKDF-Extract: PRK = HMAC-SHA256(salt, IKM).
///
/// `salt`: optional random non-secret value; pass `None` to use the default (32 zero bytes).
pub fn hkdf_extract(salt: Option<&[u8]>, ikm: &[u8]) -> [u8; 32] {
    let default_salt = [0u8; HASH_LEN];
    let s = salt.unwrap_or(&default_salt);
    hmac_sha256(s, ikm)
}

/// HKDF-Expand: derive `length` bytes of OKM from `prk` and `info`.
///
/// `length` must be ≤ 255 × HashLen = 8160 bytes.
pub fn hkdf_expand(prk: &[u8; 32], info: &[u8], length: usize) -> Vec<u8> {
    assert!(length <= 255 * HASH_LEN, "HKDF output too long");
    let n = (length + HASH_LEN - 1) / HASH_LEN;
    let mut okm = Vec::with_capacity(n * HASH_LEN);
    let mut t = Vec::new(); // T(0) = empty

    for i in 1..=n {
        let mut input = t.clone();
        input.extend_from_slice(info);
        input.push(i as u8);
        t = hmac_sha256(prk, &input).to_vec();
        okm.extend_from_slice(&t);
    }
    okm.truncate(length);
    okm
}

/// One-shot HKDF: extract then expand.
pub fn hkdf(
    salt: Option<&[u8]>,
    ikm: &[u8],
    info: &[u8],
    length: usize,
) -> Vec<u8> {
    let prk = hkdf_extract(salt, ikm);
    hkdf_expand(&prk, info, length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hmac_sha256_rfc4231_tc1() {
        // RFC 4231 Test Case 1
        let key = [0x0b; 20];
        let data = b"Hi There";
        let mac = hmac_sha256(&key, data);
        let expected = hex::decode(
            "b0344c61d8db38535ca8afceaf0bf12b881dc200c9833da726e9376c2e32cff7",
        ).unwrap();
        assert_eq!(&mac, expected.as_slice());
    }

    #[test]
    fn hkdf_rfc5869_tc1() {
        // RFC 5869 Appendix A, Test Case 1 (SHA-256)
        let ikm  = hex::decode("0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b").unwrap();
        let salt = hex::decode("000102030405060708090a0b0c").unwrap();
        let info = hex::decode("f0f1f2f3f4f5f6f7f8f9").unwrap();
        let okm  = hkdf(Some(&salt), &ikm, &info, 42);
        let expected = hex::decode(
            "3cb25f25faacd57a90434f64d0362f2a2d2d0a90cf1a5a4c5db02d56ecc4c5bf34007208d5b887185865",
        ).unwrap();
        assert_eq!(okm, expected);
    }
}
