//! # Challenge 1 — Convert hex to base64
//!
//! The first cryptopals challenge.  Decode the hex, re-encode as
//! base64.  Trivial; mostly a setup check for the rest of the set.

use crate::cryptopals::low_util::{b64_encode, hex_decode};
use crate::cryptopals::Report;

const INPUT_HEX: &str = "49276d206b696c6c696e6720796f757220627261696e206c696b65206120706f69736f6e6f7573206d757368726f6f6d";
const EXPECTED_B64: &str =
    "SSdtIGtpbGxpbmcgeW91ciBicmFpbiBsaWtlIGEgcG9pc29ub3VzIG11c2hyb29t";

pub fn run() -> Report {
    let mut r = Report::new(1, "Convert hex to base64");
    let bytes = hex_decode(INPUT_HEX);
    let got = b64_encode(&bytes);
    r.line(format!("input  : {}", INPUT_HEX));
    r.line(format!("decoded: {:?}", std::str::from_utf8(&bytes).unwrap()));
    r.line(format!("base64 : {}", got));
    r.line(format!("match  : {}", got == EXPECTED_B64));
    assert_eq!(got, EXPECTED_B64);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn matches_expected() {
        let r = super::run();
        assert!(r.success);
    }
}
