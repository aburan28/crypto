//! # Challenge 8 — Detect AES in ECB mode
//!
//! Among ~200 hex-encoded ciphertexts, one was encrypted with ECB —
//! detectable because identical 16-byte plaintext blocks produce
//! identical ciphertext blocks.  Find the duplicate-block winner.

use crate::cryptopals::low_util::hex_decode;
use crate::cryptopals::Report;
use std::collections::HashSet;

const DATA: &str = include_str!("data_8.txt");

pub fn run() -> Report {
    let mut r = Report::new(8, "Detect AES in ECB mode");
    let mut best: Option<(usize, String, usize)> = None;
    for (i, line) in DATA.lines().enumerate() {
        let ct = hex_decode(line.trim());
        let mut seen: HashSet<&[u8]> = HashSet::new();
        let mut dups = 0;
        for chunk in ct.chunks_exact(16) {
            if !seen.insert(chunk) {
                dups += 1;
            }
        }
        if dups > 0 && best.as_ref().map_or(true, |(_, _, d)| dups > *d) {
            best = Some((i, line.trim().to_string(), dups));
        }
    }
    let (i, line, dups) = best.unwrap();
    r.line(format!("ECB ciphertext at line #{}", i + 1));
    r.line(format!("Duplicate 16-byte blocks: {}", dups));
    r.line(format!("ct (truncated): {}…", &line[..60]));
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn picks_one_ct() {
        assert!(super::run().success);
    }
}
