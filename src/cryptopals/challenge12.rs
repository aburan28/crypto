//! # Challenge 12 — Byte-at-a-time ECB decryption (simple)
//!
//! An oracle ECB-encrypts `attacker || target_secret` under a fixed
//! key.  Recover `target_secret` one byte at a time:
//!
//! 1. Find block size: feed inputs of length 1, 2, 3, … and watch
//!    when the ciphertext length jumps.
//! 2. Detect ECB.
//! 3. For each unknown byte: craft an input of length
//!    `block_size - 1 - (known_len % block_size)` so the next byte
//!    of the secret is the last byte of a block.  Try all 256
//!    candidates as the trailing byte; whichever produces a match
//!    is the right one.

use crate::cryptopals::low_util::b64_decode;
use crate::cryptopals::Report;
use crate::symmetric::aes::{encrypt_block, AesKey};

const TARGET_B64: &str = "Um9sbGluJyBpbiBteSA1LjAKV2l0aCBteSByYWctdG9wIGRvd24gc28gbXkgaGFpciBjYW4gYmxvdwpUaGUgZ2lybGllcyBvbiBzdGFuZGJ5IHdhdmluZyBqdXN0IHRvIHNheSBoaQpEaWQgeW91IHN0b3A/IE5vLCBJIGp1c3QgZHJvdmUgYnkK";

fn ecb_encrypt(pt: &[u8], key: &AesKey) -> Vec<u8> {
    let padded = crate::cryptopals::low_util::pkcs7_pad(pt, 16);
    let mut out = Vec::new();
    for chunk in padded.chunks_exact(16) {
        let mut b = [0u8; 16];
        b.copy_from_slice(chunk);
        out.extend_from_slice(&encrypt_block(&b, key));
    }
    out
}

pub fn run() -> Report {
    let mut r = Report::new(12, "Byte-at-a-time ECB decryption (simple)");
    let key_bytes: [u8; 16] = *b"YELLOW_SUBMARINE";
    let key = AesKey::new(&key_bytes).unwrap();
    let secret = b64_decode(TARGET_B64);
    let oracle = |attacker: &[u8]| {
        let mut input = attacker.to_vec();
        input.extend_from_slice(&secret);
        ecb_encrypt(&input, &key)
    };

    // 1. Block size detection.
    let baseline = oracle(&[]).len();
    let mut block_size = 0;
    for i in 1..=64 {
        let l = oracle(&vec![b'A'; i]).len();
        if l != baseline {
            block_size = l - baseline;
            break;
        }
    }
    r.line(format!("Detected block size: {}", block_size));
    assert_eq!(block_size, 16);

    // 2. ECB confirmation.
    let probe = oracle(&vec![b'A'; 32]);
    assert_eq!(&probe[0..16], &probe[16..32]);
    r.line("ECB confirmed (duplicate cipher blocks on duplicate input).");

    // 3. Byte-at-a-time recovery.
    let mut recovered: Vec<u8> = Vec::new();
    let secret_len = secret.len();
    while recovered.len() < secret_len {
        let block_idx = recovered.len() / block_size;
        let pad_len = block_size - 1 - (recovered.len() % block_size);
        let prefix = vec![b'A'; pad_len];
        let target_ct = oracle(&prefix);
        let target_block =
            &target_ct[block_idx * block_size..(block_idx + 1) * block_size];
        let mut found = false;
        for b in 0..=255u8 {
            let mut candidate = prefix.clone();
            candidate.extend_from_slice(&recovered);
            candidate.push(b);
            let cand_ct = oracle(&candidate);
            let cand_block =
                &cand_ct[block_idx * block_size..(block_idx + 1) * block_size];
            if cand_block == target_block {
                recovered.push(b);
                found = true;
                break;
            }
        }
        if !found {
            // PKCS#7 trailer pad — we hit it; stop.
            break;
        }
    }
    r.line(format!("Recovered head: {:?}", &String::from_utf8_lossy(&recovered)[..40]));
    assert!(recovered.starts_with(b"Rollin' in my 5.0"));
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_secret() {
        assert!(super::run().success);
    }
}
