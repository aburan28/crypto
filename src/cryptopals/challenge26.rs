//! # Challenge 26 — CTR bit-flipping
//!
//! Same wrapper as #16 but the cipher is CTR.  In CTR, every
//! ciphertext byte XORs against an independent keystream byte, so
//! flipping ciphertext bit `b` flips plaintext bit `b` directly —
//! no block-corruption side effect.

use crate::cryptopals::challenge18::ctr_xor;
use crate::cryptopals::Report;
use crate::symmetric::aes::AesKey;

const PREFIX: &[u8] = b"comment1=cooking%20MCs;userdata=";
const SUFFIX: &[u8] = b";comment2=%20like%20a%20pound%20of%20bacon";

fn oracle_encrypt(user: &[u8], key: &AesKey, nonce: u64) -> Vec<u8> {
    let mut sanitised: Vec<u8> = user.iter().filter(|b| **b != b';' && **b != b'=').copied().collect();
    let mut full = PREFIX.to_vec();
    full.append(&mut sanitised);
    full.extend_from_slice(SUFFIX);
    ctr_xor(&full, key, nonce)
}

fn is_admin(ct: &[u8], key: &AesKey, nonce: u64) -> bool {
    let pt = ctr_xor(ct, key, nonce);
    pt.windows(12).any(|w| w == b";admin=true;")
}

pub fn run() -> Report {
    let mut r = Report::new(26, "CTR bit-flipping");
    let key = AesKey::new(b"CTR16-attack--KK").unwrap();
    let nonce = 0xC0FFEE_C0FFEEu64;
    let attacker = b"AAAAAAAAAAAAAAAA";
    let mut ct = oracle_encrypt(attacker, &key, nonce);
    // Position of attacker bytes in plaintext: starts at index 32
    // (PREFIX.len()).  Flip in place: ct[32+i] ^= attacker[i] ^ desired[i].
    let desired = b";admin=true;AAAA";
    for i in 0..12 {
        ct[32 + i] ^= attacker[i] ^ desired[i];
    }
    let ok = is_admin(&ct, &key, nonce);
    r.line(format!("Admin true: {}", ok));
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
