//! # Challenge 5 — Implement repeating-key XOR
//!
//! Encrypt the given lyrics under the repeating key `ICE`; verify
//! the output matches the expected hex.

use crate::cryptopals::low_util::{hex_encode, xor_repeating};
use crate::cryptopals::Report;

const PT: &[u8] =
    b"Burning 'em, if you ain't quick and nimble\nI go crazy when I hear a cymbal";
const KEY: &[u8] = b"ICE";
const EXPECTED_HEX: &str = "0b3637272a2b2e63622c2e69692a23693a2a3c6324202d623d63343c2a26226324272765272a282b2f20430a652e2c652a3124333a653e2b2027630c692b20283165286326302e27282f";

pub fn run() -> Report {
    let mut r = Report::new(5, "Implement repeating-key XOR");
    let ct = xor_repeating(PT, KEY);
    let hx = hex_encode(&ct);
    r.line(format!("plaintext: {}", String::from_utf8_lossy(PT)));
    r.line(format!("key      : {:?}", std::str::from_utf8(KEY).unwrap()));
    r.line(format!("ct hex   : {}", hx));
    r.line(format!("match    : {}", hx == EXPECTED_HEX));
    assert_eq!(hx, EXPECTED_HEX);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn matches_known_ct() {
        assert!(super::run().success);
    }
}
