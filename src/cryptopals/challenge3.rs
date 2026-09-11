//! # Challenge 3 — Single-byte XOR cipher
//!
//! A ciphertext was XORed against a single byte.  Brute-force all
//! 256 keys and pick the most English-looking decryption.

use crate::cryptopals::low_util::{break_single_xor, hex_decode};
use crate::cryptopals::Report;

const CT_HEX: &str =
    "1b37373331363f78151b7f2b783431333d78397828372d363c78373e783a393b3736";

pub fn run() -> Report {
    let mut r = Report::new(3, "Single-byte XOR cipher");
    let ct = hex_decode(CT_HEX);
    let (key, pt, score) = break_single_xor(&ct);
    r.line(format!("ct hex : {}", CT_HEX));
    r.line(format!("key    : 0x{key:02x} ('{}')", key as char));
    r.line(format!("pt     : {:?}", std::str::from_utf8(&pt).unwrap_or("?")));
    r.line(format!("score  : {score:.1}"));
    let expected_pt = b"Cooking MC's like a pound of bacon";
    assert_eq!(pt, expected_pt);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_known_plaintext() {
        assert!(super::run().success);
    }
}
