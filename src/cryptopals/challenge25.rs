//! # Challenge 25 — Break "random-access" CTR
//!
//! An API lets the client edit any range of plaintext given the
//! original ciphertext (it re-XORs against the same keystream).
//! That's catastrophic: send `edit(ct, 0, "\x00"·n)` and the
//! returned bytes ARE the keystream.  XOR with the original
//! ciphertext and you have the plaintext.

use crate::cryptopals::challenge18::ctr_xor;
use crate::cryptopals::low_util::{b64_decode, pkcs7_unpad};
use crate::cryptopals::Report;
use crate::symmetric::aes::{decrypt_block, AesKey};

const DATA: &str = include_str!("data_25.txt");

fn ctr_keystream_at(key: &AesKey, nonce: u64, offset: usize, len: usize) -> Vec<u8> {
    // Wrapper using ctr_xor against a zero buffer at the right offset.
    let mut zero = vec![0u8; offset + len];
    let stream = ctr_xor(&zero, key, nonce);
    // We can't seek with the simple ctr_xor; reconstruct via offset.
    zero[..].copy_from_slice(&stream[..offset + len]);
    zero[offset..offset + len].to_vec()
}

pub fn edit(ct: &[u8], key: &AesKey, nonce: u64, offset: usize, new_pt: &[u8]) -> Vec<u8> {
    let ks = ctr_keystream_at(key, nonce, offset, new_pt.len());
    let mut out = ct.to_vec();
    for i in 0..new_pt.len() {
        out[offset + i] = new_pt[i] ^ ks[i];
    }
    out
}

pub fn run() -> Report {
    let mut r = Report::new(25, "Break random-access CTR");
    // Decrypt provided ECB data first.
    let key_ecb = AesKey::new(b"YELLOW SUBMARINE").unwrap();
    let raw = b64_decode(DATA);
    let mut ecb_pt = Vec::new();
    for chunk in raw.chunks_exact(16) {
        let mut b = [0u8; 16];
        b.copy_from_slice(chunk);
        ecb_pt.extend_from_slice(&decrypt_block(&b, &key_ecb));
    }
    let pt = pkcs7_unpad(&ecb_pt, 16).unwrap();
    // Re-encrypt under CTR with a fresh key.
    let key = AesKey::new(b"ATTACKMECTR_____").unwrap();
    let ct = ctr_xor(&pt, &key, 0xCAFEBABE);
    // Attacker: query edit() with all-zero plaintext to obtain ks.
    let zero = vec![0u8; ct.len()];
    let leaked = edit(&ct, &key, 0xCAFEBABE, 0, &zero);
    // leaked = new_pt XOR ks = ks.  Recover original pt by XORing ct with ks.
    let recovered: Vec<u8> = ct.iter().zip(&leaked).map(|(c, k)| c ^ k).collect();
    r.line(format!("Recovered head: {:?}", &String::from_utf8_lossy(&recovered)[..40]));
    assert_eq!(recovered, pt);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn leaks_keystream() {
        assert!(super::run().success);
    }
}
