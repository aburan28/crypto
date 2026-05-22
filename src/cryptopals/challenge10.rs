//! # Challenge 10 — Implement CBC mode
//!
//! Decrypt the supplied base64 CBC ciphertext (key `YELLOW SUBMARINE`,
//! IV all-zero).  The repo already has a CBC implementation; we
//! reuse it to validate.

use crate::cryptopals::low_util::{b64_decode, pkcs7_unpad, xor_bytes};
use crate::cryptopals::Report;
use crate::symmetric::aes::{decrypt_block, encrypt_block, AesKey};

const DATA: &str = include_str!("data_10.txt");
const KEY: &[u8; 16] = b"YELLOW SUBMARINE";

pub fn cbc_encrypt_no_iv_prefix(pt: &[u8], key: &AesKey, iv: &[u8; 16]) -> Vec<u8> {
    let padded = crate::cryptopals::low_util::pkcs7_pad(pt, 16);
    let mut out = Vec::with_capacity(padded.len());
    let mut prev = *iv;
    for chunk in padded.chunks_exact(16) {
        let mut block = [0u8; 16];
        block.copy_from_slice(chunk);
        let xored = xor_bytes(&block, &prev);
        let mut x = [0u8; 16];
        x.copy_from_slice(&xored);
        let ct = encrypt_block(&x, key);
        out.extend_from_slice(&ct);
        prev = ct;
    }
    out
}

pub fn cbc_decrypt_no_iv_prefix(ct: &[u8], key: &AesKey, iv: &[u8; 16]) -> Option<Vec<u8>> {
    if ct.len() % 16 != 0 {
        return None;
    }
    let mut out = Vec::with_capacity(ct.len());
    let mut prev = *iv;
    for chunk in ct.chunks_exact(16) {
        let mut block = [0u8; 16];
        block.copy_from_slice(chunk);
        let dec = decrypt_block(&block, key);
        let pt = xor_bytes(&dec, &prev);
        out.extend_from_slice(&pt);
        prev = block;
    }
    pkcs7_unpad(&out, 16)
}

pub fn run() -> Report {
    let mut r = Report::new(10, "Implement CBC mode");
    let ct = b64_decode(DATA);
    let key = AesKey::new(KEY).unwrap();
    let iv = [0u8; 16];
    let pt = cbc_decrypt_no_iv_prefix(&ct, &key, &iv).expect("CBC decrypt");
    r.line(format!("plaintext head: {:?}", &String::from_utf8_lossy(&pt)[..40]));
    // Round-trip sanity.
    let re_ct = cbc_encrypt_no_iv_prefix(&pt, &key, &iv);
    assert_eq!(re_ct, ct);
    assert!(pt.starts_with(b"I'm back and I'm ringin' the bell"));
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn cbc_round_trips_and_decrypts() {
        assert!(super::run().success);
    }
}
