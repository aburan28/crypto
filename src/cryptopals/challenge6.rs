//! # Challenge 6 — Break repeating-key XOR
//!
//! Vigenère-style XOR cipher.  Standard recovery:
//!
//! 1. Guess the key length by averaging Hamming distance between
//!    blocks at each candidate length — the right length minimises
//!    normalised Hamming distance.
//! 2. Slice the ciphertext into columns by position mod key length.
//! 3. For each column, run single-byte XOR recovery.
//! 4. Stitch the single-byte keys back together → full key.

use crate::cryptopals::low_util::{
    b64_decode, break_single_xor, hamming, xor_repeating,
};
use crate::cryptopals::Report;

const DATA: &str = include_str!("data_6.txt");

fn guess_keysize(ct: &[u8], min: usize, max: usize) -> Vec<usize> {
    let mut scores: Vec<(usize, f64)> = Vec::new();
    for k in min..=max.min(ct.len() / 4) {
        // Average normalised Hamming distance over up to 4 chunks.
        let mut total = 0.0f64;
        let mut count = 0;
        for i in 0..3 {
            if (i + 2) * k > ct.len() {
                break;
            }
            let a = &ct[i * k..(i + 1) * k];
            let b = &ct[(i + 1) * k..(i + 2) * k];
            total += hamming(a, b) as f64 / k as f64;
            count += 1;
        }
        if count > 0 {
            scores.push((k, total / count as f64));
        }
    }
    scores.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    scores.into_iter().take(5).map(|(k, _)| k).collect()
}

pub fn run() -> Report {
    let mut r = Report::new(6, "Break repeating-key XOR");
    let ct = b64_decode(DATA);
    let candidates = guess_keysize(&ct, 2, 40);
    r.line(format!("Top keysize candidates: {:?}", candidates));
    // Try each candidate; pick the recovered plaintext with best
    // English score.
    let mut best: Option<(usize, Vec<u8>, Vec<u8>)> = None;
    for &k in &candidates {
        let mut key = vec![0u8; k];
        for col in 0..k {
            let column: Vec<u8> = ct.iter().enumerate()
                .filter_map(|(i, &b)| if i % k == col { Some(b) } else { None })
                .collect();
            let (kb, _, _) = break_single_xor(&column);
            key[col] = kb;
        }
        let pt = xor_repeating(&ct, &key);
        let score = crate::cryptopals::low_util::score_english(&pt);
        if best.as_ref().map_or(true, |(_, prev, _)| {
            crate::cryptopals::low_util::score_english(prev) < score
        }) {
            best = Some((k, pt, key));
        }
    }
    let (k, pt, key) = best.unwrap();
    r.line(format!("Chosen keysize : {}", k));
    r.line(format!("Recovered key  : {:?}", String::from_utf8_lossy(&key)));
    r.line(format!("Plaintext head : {:?}", &String::from_utf8_lossy(&pt)[..60]));
    // The Vanilla Ice masterpiece "Play That Funky Music."
    let expected_key = b"Terminator X: Bring the noise";
    assert_eq!(key.as_slice(), expected_key);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_known_key() {
        assert!(super::run().success);
    }
}
