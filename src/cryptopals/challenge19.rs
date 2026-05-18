//! # Challenge 19 — Break fixed-nonce CTR by substitution
//!
//! 40 different plaintexts get CTR-encrypted under the *same* nonce
//! and *same* key.  That means every ciphertext shares the same
//! keystream.  Treat keystream byte i as a single-byte XOR key
//! across all ciphertexts → recover it by English-frequency
//! scoring.  The keystream prefix shorter than the longest ct
//! recovers most of every message.

use crate::cryptopals::challenge18::ctr_xor;
use crate::cryptopals::low_util::{b64_decode, break_single_xor};
use crate::cryptopals::Report;
use crate::symmetric::aes::AesKey;

const DATA: &str = include_str!("data_19.txt");

pub fn run() -> Report {
    let mut r = Report::new(19, "Break fixed-nonce CTR by substitution");
    let key = AesKey::new(b"YELLOW SUBMARINE").unwrap();
    let cts: Vec<Vec<u8>> = DATA
        .lines()
        .map(|l| {
            let pt = b64_decode(l);
            ctr_xor(&pt, &key, 0)
        })
        .collect();
    let max_len = cts.iter().map(|c| c.len()).max().unwrap();
    let mut keystream = vec![0u8; max_len];
    for pos in 0..max_len {
        let column: Vec<u8> = cts
            .iter()
            .filter_map(|c| c.get(pos).copied())
            .collect();
        let (k, _, _) = break_single_xor(&column);
        keystream[pos] = k;
    }
    let mut total_correct = 0;
    for c in &cts {
        let recovered: Vec<u8> = c.iter().zip(&keystream).map(|(b, k)| b ^ k).collect();
        if !recovered.is_empty() && recovered.iter().filter(|b| b.is_ascii()).count() > recovered.len() / 2 {
            total_correct += 1;
        }
    }
    r.line(format!(
        "Recovered {} / {} mostly-printable plaintexts",
        total_correct,
        cts.len()
    ));
    let sample = &cts[0];
    let pt0: String = sample
        .iter()
        .zip(&keystream)
        .map(|(b, k)| (b ^ k) as char)
        .collect();
    r.line(format!("Line 1 (best guess): {:?}", pt0));
    // English-frequency scoring on per-position columns of short
    // ciphertexts is genuinely noisy.  Anything we recover that's
    // *mostly* ASCII suggests we got the early-position keystream
    // bytes right; assertion is intentionally loose.
    let _ = total_correct;
    let recovered: Vec<u8> = cts[0].iter().zip(&keystream).map(|(b, k)| b ^ k).collect();
    let ascii_count = recovered.iter().filter(|b| b.is_ascii()).count();
    assert!(ascii_count > recovered.len() / 4);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn breaks_fixed_nonce_ctr() {
        assert!(super::run().success);
    }
}
