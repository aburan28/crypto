//! # Challenge 9 — Implement PKCS#7 padding
//!
//! `pkcs7_pad("YELLOW SUBMARINE", 20)` → `"YELLOW SUBMARINE\x04\x04\x04\x04"`.

use crate::cryptopals::low_util::pkcs7_pad;
use crate::cryptopals::Report;

pub fn run() -> Report {
    let mut r = Report::new(9, "PKCS#7 padding");
    let padded = pkcs7_pad(b"YELLOW SUBMARINE", 20);
    r.line(format!("padded: {:?}", padded));
    assert_eq!(padded, b"YELLOW SUBMARINE\x04\x04\x04\x04");
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn pads_correctly() {
        assert!(super::run().success);
    }
}
