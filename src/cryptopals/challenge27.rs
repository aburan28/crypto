//! # Challenge 27 — Recover key from CBC where IV = key
//!
//! If a CBC implementation uses the same value for both the
//! encryption key AND the IV (looking at you, real-world TLS
//! implementations c.2010), feeding a 3-block ciphertext
//! `C₁, 0, C₁` and inspecting the decryption error message gives
//! you back `P₁ ⊕ P₃ = K`.

use crate::cryptopals::challenge10::cbc_encrypt_no_iv_prefix;
use crate::cryptopals::Report;
use crate::symmetric::aes::{decrypt_block, AesKey};

fn cbc_decrypt_raw(ct: &[u8], key: &AesKey, iv: &[u8; 16]) -> Vec<u8> {
    let mut out = Vec::with_capacity(ct.len());
    let mut prev = *iv;
    for chunk in ct.chunks_exact(16) {
        let mut b = [0u8; 16];
        b.copy_from_slice(chunk);
        let dec = decrypt_block(&b, key);
        for i in 0..16 {
            out.push(dec[i] ^ prev[i]);
        }
        prev = b;
    }
    out
}

pub fn run() -> Report {
    let mut r = Report::new(27, "CBC with IV = key — key recovery");
    let secret_key: [u8; 16] = *b"key=iv-secret-yz";
    let key = AesKey::new(&secret_key).unwrap();
    let iv = secret_key;
    // Choose a 48-byte plaintext so ct is exactly 3 blocks (no pad).
    let pt = b"ABCDEFGHIJKLMNOPABCDEFGHIJKLMNOPABCDEFGHIJKLMNOP"; // 48
    let ct = cbc_encrypt_no_iv_prefix(pt, &key, &iv);
    assert_eq!(ct.len(), 48 + 16); // 48 bytes pt + one block PKCS7 pad
    // Take only the first 48 bytes (3 blocks).
    let mut malformed = vec![0u8; 48];
    malformed[..16].copy_from_slice(&ct[..16]);
    // block 1 = all zeros
    malformed[32..48].copy_from_slice(&ct[..16]); // copy of block 1
    let decoded = cbc_decrypt_raw(&malformed, &key, &iv);
    // Plaintext leakage: P1 ⊕ P3 = K (since block 2 is 0,
    // P2 = D(C1) ⊕ C0_block_2 = D(C1) ⊕ 0 = D(C1), and P3 = D(C1) ⊕ 0...
    // Actually the standard formula: K = P1 ⊕ P3.
    let p1 = &decoded[0..16];
    let p3 = &decoded[32..48];
    let mut recovered = [0u8; 16];
    for i in 0..16 {
        recovered[i] = p1[i] ^ p3[i];
    }
    r.line(format!("True key       : {}", hex::encode(secret_key)));
    r.line(format!("Recovered key  : {}", hex::encode(recovered)));
    assert_eq!(recovered, secret_key);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_key() {
        assert!(super::run().success);
    }
}
