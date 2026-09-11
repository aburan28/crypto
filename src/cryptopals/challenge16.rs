//! # Challenge 16 — CBC bit-flipping
//!
//! Oracle wraps attacker data in `"comment1=cooking%20MCs;userdata="`
//! and `";comment2=%20like%20a%20pound%20of%20bacon"`, quotes any
//! `;` and `=` chars, then CBC-encrypts.  The decryption oracle
//! returns true iff the plaintext contains `;admin=true;`.
//!
//! Attack: pick attacker input "AAAAAAAAAAAAAAAA" (one block), then
//! flip bits in the *previous* ciphertext block so that the
//! plaintext block decrypts to `;admin=true;....`.  In CBC, flipping
//! ciphertext block N's byte affects plaintext block N+1's same-
//! position byte directly (XOR).

use crate::cryptopals::challenge10::{cbc_decrypt_no_iv_prefix, cbc_encrypt_no_iv_prefix};
use crate::cryptopals::Report;
use crate::symmetric::aes::AesKey;

const PREFIX: &[u8] = b"comment1=cooking%20MCs;userdata=";
const SUFFIX: &[u8] = b";comment2=%20like%20a%20pound%20of%20bacon";

pub fn oracle_encrypt(user_data: &[u8], key: &AesKey, iv: &[u8; 16]) -> Vec<u8> {
    let mut sanitised: Vec<u8> = user_data
        .iter()
        .filter(|b| **b != b';' && **b != b'=')
        .copied()
        .collect();
    let mut full = PREFIX.to_vec();
    full.append(&mut sanitised);
    full.extend_from_slice(SUFFIX);
    cbc_encrypt_no_iv_prefix(&full, key, iv)
}

pub fn oracle_check_admin(ct: &[u8], key: &AesKey, iv: &[u8; 16]) -> bool {
    match cbc_decrypt_no_iv_prefix(ct, key, iv) {
        None => false,
        Some(pt) => pt.windows(13).any(|w| w == b";admin=true;a="[..12].as_ref())
            || pt.windows(12).any(|w| w == b";admin=true;"),
    }
}

pub fn run() -> Report {
    let mut r = Report::new(16, "CBC bit-flipping");
    let key_bytes: [u8; 16] = *b"keyfor16-xx-yyzz";
    let key = AesKey::new(&key_bytes).unwrap();
    let iv = [0u8; 16];

    // 1. Prefix is exactly 32 bytes = 2 blocks.  Attacker block
    //    (block index 2) is fully attacker-controlled.
    //
    // 2. Strategy: put a known plaintext "XXXXXXXXXXXXXXXX" (16
    //    'X's) into block 2, then flip bits in ciphertext block 1
    //    so that block 2 decrypts to ";admin=true;XXXX".  Block 1's
    //    decryption is destroyed but the bank doesn't check that.
    let attacker = b"AAAAAAAAAAAAAAAA"; // 16 'A's
    let mut ct = oracle_encrypt(attacker, &key, &iv);
    // Block indexes in ct: block 0 = comment1=cooking, block 1 = %20MCs;userdata=,
    // block 2 = attacker bytes.  PREFIX is 32 = 2 blocks.
    // To turn "AAAAAAAAAAAAAAAA" into ";admin=true;AAAA", XOR
    // ciphertext block 1 byte k with ('A' ^ desired[k]) for k in 0..12.
    let desired = b";admin=true;AAAA";
    for i in 0..12 {
        ct[16 + i] ^= attacker[i] ^ desired[i];
    }
    let ok = oracle_check_admin(&ct, &key, &iv);
    r.line(format!("Oracle reports admin=true : {}", ok));
    assert!(ok);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn flips_in_admin_true() {
        assert!(super::run().success);
    }
}
