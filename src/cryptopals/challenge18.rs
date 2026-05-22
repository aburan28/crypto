//! # Challenge 18 — AES-CTR (decrypt + verify round-trip)
//!
//! CTR uses a 64-bit little-endian nonce and a 64-bit little-endian
//! block counter (cryptopals format) to derive a keystream.

use crate::cryptopals::low_util::b64_decode;
use crate::cryptopals::Report;
use crate::symmetric::aes::{encrypt_block, AesKey};

/// Apply CTR keystream to `data` in place.  `nonce` is 8 bytes
/// little-endian; counter starts at 0 and is also 8 bytes LE.
pub fn ctr_xor(data: &[u8], key: &AesKey, nonce: u64) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len());
    let mut counter: u64 = 0;
    for chunk in data.chunks(16) {
        let mut input = [0u8; 16];
        input[..8].copy_from_slice(&nonce.to_le_bytes());
        input[8..].copy_from_slice(&counter.to_le_bytes());
        let stream = encrypt_block(&input, key);
        for (i, &b) in chunk.iter().enumerate() {
            out.push(b ^ stream[i]);
        }
        counter += 1;
    }
    out
}

const CT_B64: &str =
    "L77na/nrFsKvynd6HzOoG7GHTLXsTVu9qvY/2syLXzhPweyyMTJULu/6/kXX0KSvoOLSFQ==";

pub fn run() -> Report {
    let mut r = Report::new(18, "AES-CTR mode");
    let ct = b64_decode(CT_B64);
    let key = AesKey::new(b"YELLOW SUBMARINE").unwrap();
    let pt = ctr_xor(&ct, &key, 0);
    r.line(format!("Plaintext: {}", String::from_utf8_lossy(&pt)));
    assert!(pt.starts_with(b"Yo, VIP Let's kick it"));
    // Round-trip.
    let round = ctr_xor(&pt, &key, 0);
    assert_eq!(round, ct);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn ctr_round_trips() {
        assert!(super::run().success);
    }
}
