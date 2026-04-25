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

    fn h(s: &str) -> Vec<u8> { hex::decode(s).unwrap() }

    // ── HMAC-SHA256 KATs (RFC 4231) ──────────────────────────────────────────

    #[test]
    fn hmac_sha256_rfc4231_tc1() {
        // RFC 4231 §4.2 — short key, "Hi There"
        let mac = hmac_sha256(&[0x0b; 20], b"Hi There");
        assert_eq!(
            &mac,
            h("b0344c61d8db38535ca8afceaf0bf12b881dc200c9833da726e9376c2e32cff7").as_slice(),
        );
    }

    #[test]
    fn hmac_sha256_rfc4231_tc2() {
        // RFC 4231 §4.3 — short key
        let mac = hmac_sha256(b"Jefe", b"what do ya want for nothing?");
        assert_eq!(
            &mac,
            h("5bdcc146bf60754e6a042426089575c75a003f089d2739839dec58b964ec3843").as_slice(),
        );
    }

    #[test]
    fn hmac_sha256_rfc4231_tc3() {
        // RFC 4231 §4.4 — combined key/data length spans block
        let mac = hmac_sha256(&[0xaa; 20], &[0xdd; 50]);
        assert_eq!(
            &mac,
            h("773ea91e36800e46854db8ebd09181a72959098b3ef8c122d9635514ced565fe").as_slice(),
        );
    }

    #[test]
    fn hmac_sha256_rfc4231_tc4() {
        // RFC 4231 §4.5 — key 1..25
        let key: Vec<u8> = (1u8..=25).collect();
        let mac = hmac_sha256(&key, &[0xcd; 50]);
        assert_eq!(
            &mac,
            h("82558a389a443c0ea4cc819899f2083a85f0faa3e578f8077a2e3ff46729665b").as_slice(),
        );
    }

    #[test]
    fn hmac_sha256_rfc4231_tc6() {
        // RFC 4231 §4.7 — long key (131 bytes 0xaa) gets hashed first
        let mac = hmac_sha256(
            &[0xaa; 131],
            b"Test Using Larger Than Block-Size Key - Hash Key First",
        );
        assert_eq!(
            &mac,
            h("60e431591ee0b67f0d8a26aacbf5b77f8e0bc6213728c5140546040f0ee37f54").as_slice(),
        );
    }

    // ── HKDF KATs (RFC 5869 Appendix A) ──────────────────────────────────────

    #[test]
    fn hkdf_rfc5869_tc1_basic() {
        // RFC 5869 §A.1
        let ikm  = h("0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b");
        let salt = h("000102030405060708090a0b0c");
        let info = h("f0f1f2f3f4f5f6f7f8f9");
        assert_eq!(
            hkdf_extract(Some(&salt), &ikm).as_slice(),
            h("077709362c2e32df0ddc3f0dc47bba6390b6c73bb50f9c3122ec844ad7c2b3e5").as_slice(),
        );
        assert_eq!(
            hkdf(Some(&salt), &ikm, &info, 42),
            h("3cb25f25faacd57a90434f64d0362f2a2d2d0a90cf1a5a4c5db02d56ecc4c5bf34007208d5b887185865"),
        );
    }

    #[test]
    fn hkdf_rfc5869_tc2_long() {
        // RFC 5869 §A.2 — longer inputs/outputs
        let ikm: Vec<u8> = (0u8..80).collect();
        let salt: Vec<u8> = (0x60u8..0xb0).collect();
        let info: Vec<u8> = (0xb0u8..0xff).chain(std::iter::once(0xff)).collect();
        let okm = hkdf(Some(&salt), &ikm, &info, 82);
        assert_eq!(
            okm,
            h("b11e398dc80327a1c8e7f78c596a49344f012eda2d4efad8a050cc4c19afa97c\
               59045a99cac7827271cb41c65e590e09da3275600c2f09b8367793a9aca3db71\
               cc30c58179ec3e87c14c01d5c1f3434f1d87"),
        );
    }

    #[test]
    fn hkdf_rfc5869_tc3_no_salt() {
        // RFC 5869 §A.3 — empty salt and info
        let ikm = h("0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b");
        let okm = hkdf(None, &ikm, &[], 42);
        assert_eq!(
            okm,
            h("8da4e775a563c18f715f802a063c5a31b8a11f5c5ee1879ec3454e5f3c738d2d\
               9d201395faa4b61a96c8"),
        );
    }
}
