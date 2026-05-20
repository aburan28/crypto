//! # Challenge 2 — Fixed XOR
//!
//! XOR two equal-length hex strings and report the result.

use crate::cryptopals::low_util::{hex_decode, hex_encode, xor_bytes};
use crate::cryptopals::Report;

const A: &str = "1c0111001f010100061a024b53535009181c";
const B: &str = "686974207468652062756c6c277320657965";
const EXPECTED: &str = "746865206b696420646f6e277420706c6179";

pub fn run() -> Report {
    let mut r = Report::new(2, "Fixed XOR");
    let a = hex_decode(A);
    let b = hex_decode(B);
    let x = xor_bytes(&a, &b);
    let hx = hex_encode(&x);
    r.line(format!("a   : {}", A));
    r.line(format!("b   : {}", B));
    r.line(format!("a⊕b : {}", hx));
    r.line(format!("text: {:?}", std::str::from_utf8(&x).unwrap()));
    assert_eq!(hx, EXPECTED);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn xor_matches() {
        assert!(super::run().success);
    }
}
