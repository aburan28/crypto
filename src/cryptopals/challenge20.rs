//! # Challenge 20 — Break fixed-nonce CTR statistically
//!
//! Same setup as #19 but a larger corpus; truncate all ciphertexts
//! to the shortest length and recover with a Vigenère-style attack.

use crate::cryptopals::challenge18::ctr_xor;
use crate::cryptopals::low_util::{b64_decode, break_single_xor};
use crate::cryptopals::Report;
use crate::symmetric::aes::AesKey;

const DATA: &str = include_str!("data_20.txt");

pub fn run() -> Report {
    let mut r = Report::new(20, "Break fixed-nonce CTR statistically");
    let key = AesKey::new(b"YELLOW SUBMARINE").unwrap();
    let cts: Vec<Vec<u8>> = DATA
        .lines()
        .map(|l| ctr_xor(&b64_decode(l), &key, 0))
        .collect();
    let n = cts.iter().map(|c| c.len()).min().unwrap();
    let mut keystream = vec![0u8; n];
    for pos in 0..n {
        let col: Vec<u8> = cts.iter().map(|c| c[pos]).collect();
        keystream[pos] = break_single_xor(&col).0;
    }
    let recovered_first: String = cts[0][..n]
        .iter()
        .zip(&keystream)
        .map(|(b, k)| (b ^ k) as char)
        .collect();
    r.line(format!("Line 1 head: {:?}", recovered_first));
    assert!(recovered_first.to_ascii_lowercase().contains("i'm rated") || recovered_first.chars().filter(|c| c.is_ascii()).count() > n * 9 / 10);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn breaks_statistically() {
        assert!(super::run().success);
    }
}
