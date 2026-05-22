//! # Challenge 7 — AES-128-ECB decryption
//!
//! Decrypt the supplied base64 ciphertext with AES-128 in ECB mode
//! and the key `YELLOW SUBMARINE`.

use crate::cryptopals::low_util::{b64_decode, pkcs7_unpad};
use crate::cryptopals::Report;
use crate::symmetric::aes::{decrypt_block, AesKey};

const DATA: &str = include_str!("data_7.txt");
const KEY: &[u8; 16] = b"YELLOW SUBMARINE";

pub fn run() -> Report {
    let mut r = Report::new(7, "AES-128-ECB decryption");
    let ct = b64_decode(DATA);
    let key = AesKey::new(KEY).unwrap();
    let mut pt = Vec::with_capacity(ct.len());
    for chunk in ct.chunks_exact(16) {
        let mut block = [0u8; 16];
        block.copy_from_slice(chunk);
        let dec = decrypt_block(&block, &key);
        pt.extend_from_slice(&dec);
    }
    let unpadded = pkcs7_unpad(&pt, 16).unwrap();
    r.line(format!("Plaintext (first line): {:?}", &String::from_utf8_lossy(&unpadded)[..40]));
    assert!(unpadded.starts_with(b"I'm back and I'm ringin' the bell"));
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn decrypts_known_ct() {
        assert!(super::run().success);
    }
}
