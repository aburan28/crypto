//! # Challenge 4 — Detect single-character XOR
//!
//! 60 hex strings are given, one of which was single-byte XOR
//! encoded.  Find it.

use crate::cryptopals::low_util::{break_single_xor, hex_decode};
use crate::cryptopals::Report;

const DATA: &str = include_str!("data_4.txt");

pub fn run() -> Report {
    let mut r = Report::new(4, "Detect single-character XOR");
    let mut best: Option<(String, u8, Vec<u8>, f64)> = None;
    for line in DATA.lines() {
        let ct = hex_decode(line.trim());
        let (k, pt, s) = break_single_xor(&ct);
        if best.as_ref().map_or(true, |(_, _, _, bs)| s > *bs) {
            best = Some((line.trim().to_string(), k, pt, s));
        }
    }
    let (l, k, pt, _) = best.unwrap();
    r.line(format!("source hex line : {}", l));
    r.line(format!("key             : 0x{k:02x}"));
    r.line(format!("plaintext       : {:?}", String::from_utf8_lossy(&pt)));
    let expected = b"Now that the party is jumping\n";
    assert_eq!(pt, expected);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn finds_xored_line() {
        assert!(super::run().success);
    }
}
