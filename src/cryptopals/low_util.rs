//! Shared helpers for Cryptopals Sets 1-6.
//!
//! Tiny utilities (XOR, English-frequency scoring, PKCS#7) that
//! show up across many challenges.  Kept here so each challenge
//! file stays focused on the *attack* rather than boilerplate.

use base64::Engine;

/// XOR two byte slices of equal length.  Panics on mismatch.
pub fn xor_bytes(a: &[u8], b: &[u8]) -> Vec<u8> {
    assert_eq!(a.len(), b.len(), "xor_bytes: length mismatch");
    a.iter().zip(b).map(|(x, y)| x ^ y).collect()
}

/// XOR each byte of `data` with the corresponding byte of `key`,
/// cycling `key`.
pub fn xor_repeating(data: &[u8], key: &[u8]) -> Vec<u8> {
    data.iter()
        .enumerate()
        .map(|(i, b)| b ^ key[i % key.len()])
        .collect()
}

/// Decode a base64 string, ignoring whitespace.
pub fn b64_decode(s: &str) -> Vec<u8> {
    let stripped: String = s.chars().filter(|c| !c.is_whitespace()).collect();
    base64::engine::general_purpose::STANDARD
        .decode(stripped)
        .expect("invalid base64")
}

/// Encode bytes as base64.
pub fn b64_encode(b: &[u8]) -> String {
    base64::engine::general_purpose::STANDARD.encode(b)
}

/// Decode a hex string.
pub fn hex_decode(s: &str) -> Vec<u8> {
    hex::decode(s).expect("invalid hex")
}

/// Encode bytes as hex.
pub fn hex_encode(b: &[u8]) -> String {
    hex::encode(b)
}

/// English-character frequency.  Source: Wikipedia "Letter frequency"
/// — roughly representative of English prose.  Lowercase only; we
/// fold case before scoring.
const ENGLISH_FREQ: [(u8, f64); 27] = [
    (b' ', 13.0), // most common in real English text by a wide margin
    (b'e', 12.7),
    (b't', 9.1),
    (b'a', 8.2),
    (b'o', 7.5),
    (b'i', 7.0),
    (b'n', 6.7),
    (b's', 6.3),
    (b'h', 6.1),
    (b'r', 6.0),
    (b'd', 4.3),
    (b'l', 4.0),
    (b'c', 2.8),
    (b'u', 2.8),
    (b'm', 2.4),
    (b'w', 2.4),
    (b'f', 2.2),
    (b'g', 2.0),
    (b'y', 2.0),
    (b'p', 1.9),
    (b'b', 1.5),
    (b'v', 1.0),
    (b'k', 0.8),
    (b'j', 0.15),
    (b'x', 0.15),
    (b'q', 0.10),
    (b'z', 0.07),
];

/// Score `text` by how "English-like" it is.  Higher is more
/// English.  Used by Set 1 to spot the winning single-byte XOR key.
pub fn score_english(text: &[u8]) -> f64 {
    let mut s = 0.0;
    for &b in text {
        let lc = b.to_ascii_lowercase();
        if let Some((_, w)) = ENGLISH_FREQ.iter().find(|(c, _)| *c == lc) {
            s += w;
        } else if b < 0x20 && b != b'\n' && b != b'\r' && b != b'\t' {
            // Heavy penalty for non-printable bytes.
            s -= 10.0;
        }
    }
    s
}

/// Try every single-byte XOR key against `ct` and return the
/// `(key, plaintext, score)` with the best English score.
pub fn break_single_xor(ct: &[u8]) -> (u8, Vec<u8>, f64) {
    let mut best_key = 0u8;
    let mut best_pt: Vec<u8> = Vec::new();
    let mut best_score = f64::NEG_INFINITY;
    for k in 0..=255u8 {
        let pt: Vec<u8> = ct.iter().map(|b| b ^ k).collect();
        let s = score_english(&pt);
        if s > best_score {
            best_score = s;
            best_key = k;
            best_pt = pt;
        }
    }
    (best_key, best_pt, best_score)
}

/// Hamming distance: number of differing bits across two equal-length
/// byte slices.  Used by Set 1 challenge 6 to guess Vigenère key
/// length.
pub fn hamming(a: &[u8], b: &[u8]) -> u32 {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).map(|(x, y)| (x ^ y).count_ones()).sum()
}

/// PKCS#7 pad to a multiple of `block_size`.  Always adds at least
/// one byte.
pub fn pkcs7_pad(data: &[u8], block_size: usize) -> Vec<u8> {
    let pad = block_size - (data.len() % block_size);
    let mut out = data.to_vec();
    out.extend(std::iter::repeat(pad as u8).take(pad));
    out
}

/// Strip and validate PKCS#7 padding.  Returns `None` on bad
/// padding.
pub fn pkcs7_unpad(data: &[u8], block_size: usize) -> Option<Vec<u8>> {
    if data.is_empty() || data.len() % block_size != 0 {
        return None;
    }
    let pad = *data.last()? as usize;
    if pad == 0 || pad > block_size {
        return None;
    }
    let n = data.len();
    for &b in &data[n - pad..] {
        if b as usize != pad {
            return None;
        }
    }
    Some(data[..n - pad].to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xor_round_trip() {
        let a = b"YELLOW SUBMARINE";
        let key = b"ICE";
        let ct = xor_repeating(a, key);
        let pt = xor_repeating(&ct, key);
        assert_eq!(pt, a);
    }

    #[test]
    fn hamming_distance_known_pair() {
        // Cryptopals' own example.
        assert_eq!(hamming(b"this is a test", b"wokka wokka!!!"), 37);
    }

    #[test]
    fn pkcs7_round_trip() {
        let data = b"hello!";
        let padded = pkcs7_pad(data, 16);
        assert_eq!(padded.len(), 16);
        assert_eq!(pkcs7_unpad(&padded, 16).unwrap(), data.to_vec());
    }
}
