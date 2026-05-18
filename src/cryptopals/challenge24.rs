//! # Challenge 24 — MT19937 stream cipher + token check
//!
//! Use MT19937 as a keystream generator (XOR each output byte
//! against a plaintext byte).  Recover the 16-bit seed from a
//! plaintext-prefix-padded ciphertext.

use crate::cryptopals::challenge21::Mt19937;
use crate::cryptopals::Report;

pub fn mt_xor(data: &[u8], seed: u16) -> Vec<u8> {
    let mut rng = Mt19937::new(seed as u32);
    let mut out = Vec::with_capacity(data.len());
    let mut buf = [0u8; 4];
    let mut buf_idx = 4;
    for &b in data {
        if buf_idx == 4 {
            buf = rng.next_u32().to_le_bytes();
            buf_idx = 0;
        }
        out.push(b ^ buf[buf_idx]);
        buf_idx += 1;
    }
    out
}

pub fn run() -> Report {
    let mut r = Report::new(24, "MT19937 stream cipher");
    let seed: u16 = 0xC0DE;
    // Encrypt random-prefix || 14·"A"
    let mut pt = b"prefix_xxx".to_vec();
    pt.extend_from_slice(b"AAAAAAAAAAAAAA");
    let ct = mt_xor(&pt, seed);
    // Recover seed: brute-force 2^16 values.
    let known_suffix = b"AAAAAAAAAAAAAA";
    let suffix_len = known_suffix.len();
    let mut recovered = None;
    for s in 0u32..=0xFFFFu32 {
        let cand = mt_xor(&ct, s as u16);
        if cand.ends_with(known_suffix) {
            recovered = Some(s as u16);
            break;
        }
    }
    let got = recovered.expect("seed in 2^16 range");
    r.line(format!("True seed     : 0x{:04x}", seed));
    r.line(format!("Recovered seed: 0x{:04x}", got));
    let _ = suffix_len;
    assert_eq!(got, seed);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn brute_force_16bit_seed() {
        assert!(super::run().success);
    }
}
