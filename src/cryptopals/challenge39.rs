//! # Challenge 39 — Implement RSA
//!
//! Use the existing RSA module to keygen / encrypt / decrypt a
//! small plaintext and verify round-trip.

use crate::cryptopals::Report;

pub fn run() -> Report {
    let mut r = Report::new(39, "Implement RSA");
    let kp = crate::asymmetric::rsa::RsaKeyPair::generate(512);
    let msg = b"hello, RSA";
    let ct = crate::asymmetric::rsa::rsa_encrypt(msg, &kp.public).unwrap();
    let pt = crate::asymmetric::rsa::rsa_decrypt(&ct, &kp.private).unwrap();
    r.line(format!("plaintext after round-trip: {:?}", std::str::from_utf8(&pt).unwrap()));
    assert_eq!(pt, msg);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn rsa_round_trip() {
        assert!(super::run().success);
    }
}
