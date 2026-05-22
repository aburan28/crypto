//! # Challenge 15 — PKCS#7 padding validation
//!
//! Strict unpad: validate the trailing padding bytes and return an
//! error on malformed padding.

use crate::cryptopals::low_util::pkcs7_unpad;
use crate::cryptopals::Report;

pub fn run() -> Report {
    let mut r = Report::new(15, "PKCS#7 padding validation");

    let good = b"ICE ICE BABY\x04\x04\x04\x04";
    let bad1 = b"ICE ICE BABY\x05\x05\x05\x05";
    let bad2 = b"ICE ICE BABY\x01\x02\x03\x04";

    let pad_good = pkcs7_unpad(good, 16);
    let pad_bad1 = pkcs7_unpad(bad1, 16);
    let pad_bad2 = pkcs7_unpad(bad2, 16);
    r.line(format!("good        : {:?}", pad_good.as_ref().map(|v| std::str::from_utf8(v).unwrap())));
    r.line(format!("bad #1      : {:?}", pad_bad1.is_none()));
    r.line(format!("bad #2      : {:?}", pad_bad2.is_none()));
    assert_eq!(pad_good.unwrap(), b"ICE ICE BABY");
    assert!(pad_bad1.is_none());
    assert!(pad_bad2.is_none());
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn validates_padding() {
        assert!(super::run().success);
    }
}
