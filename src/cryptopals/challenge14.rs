//! # Challenge 14 — Byte-at-a-time ECB decryption (harder)
//!
//! Same as #12, but the oracle prepends a random *fixed* prefix
//! before the attacker's bytes and the secret.  Strategy:
//!
//! 1. Find the prefix length: send increasingly long attacker
//!    inputs and watch when two adjacent ciphertext blocks become
//!    identical.  Use the alignment of that pair to deduce
//!    `prefix_len mod block_size` and `prefix_len / block_size`.
//! 2. Pad your attacker input to align the secret to a block
//!    boundary, then do the byte-at-a-time attack as in #12 with
//!    an offset.

use crate::cryptopals::low_util::b64_decode;
use crate::cryptopals::Report;
use crate::symmetric::aes::{encrypt_block, AesKey};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

const TARGET_B64: &str = "Um9sbGluJyBpbiBteSA1LjAKV2l0aCBteSByYWctdG9wIGRvd24gc28gbXkgaGFpciBjYW4gYmxvdwpUaGUgZ2lybGllcyBvbiBzdGFuZGJ5IHdhdmluZyBqdXN0IHRvIHNheSBoaQpEaWQgeW91IHN0b3A/IE5vLCBJIGp1c3QgZHJvdmUgYnkK";

fn ecb_enc(pt: &[u8], key: &AesKey) -> Vec<u8> {
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
    let mut r = Report::new(14, "Byte-at-a-time ECB decryption (harder)");
    let key_bytes: [u8; 16] = *b"keyforC14XXXyyyy";
    let key = AesKey::new(&key_bytes).unwrap();
    let mut rng = StdRng::seed_from_u64(54321);
    let prefix_len = rng.gen_range(1..=64);
    let prefix: Vec<u8> = (0..prefix_len).map(|_| rng.gen::<u8>()).collect();
    let secret = b64_decode(TARGET_B64);
    let oracle = |attacker: &[u8]| {
        let mut input = prefix.clone();
        input.extend_from_slice(attacker);
        input.extend_from_slice(&secret);
        ecb_enc(&input, &key)
    };
    let block_size = 16usize;
    r.line(format!("(hidden) prefix len = {}", prefix_len));

    // Step 1: find prefix_len.  Send N attacker bytes; pick the
    // smallest N such that some two adjacent ciphertext blocks
    // contain repetition we control.  Concretely: send
    // 2·block_size attacker bytes + extra to align.  Iterate
    // attacker length from 0..block_size and look for two
    // consecutive identical blocks in the cipher.
    let mut prefix_blocks_partial = 0;
    let mut align_pad = 0;
    'outer: for pad in 0..block_size {
        let pt = vec![b'A'; pad + 2 * block_size];
        let ct = oracle(&pt);
        // Look for two consecutive equal blocks.
        for i in 0..(ct.len() / block_size - 1) {
            if ct[i * block_size..(i + 1) * block_size]
                == ct[(i + 1) * block_size..(i + 2) * block_size]
            {
                prefix_blocks_partial = i;
                align_pad = pad;
                break 'outer;
            }
        }
    }
    // For attacker-pad = `align_pad`, the duplicate blocks start at
    // block index `prefix_blocks_partial`.  That block index times
    // block_size minus the pad gives the true prefix length:
    //   prefix_len = i · block_size − pad   (in [0, block_size·i))
    let detected_prefix_len = prefix_blocks_partial * block_size - align_pad;
    r.line(format!("Detected prefix len: {}", detected_prefix_len));
    assert_eq!(detected_prefix_len, prefix_len);

    // Step 2: recover the secret byte-by-byte, accounting for
    // the prefix.
    let prefix_pad = (block_size - prefix_len % block_size) % block_size;
    let prefix_blocks = (prefix_len + prefix_pad) / block_size;
    let pad_attacker = vec![b'A'; prefix_pad];
    let mut recovered: Vec<u8> = Vec::new();
    let secret_len = secret.len();
    while recovered.len() < secret_len {
        let block_idx = prefix_blocks + recovered.len() / block_size;
        let pad_extra = block_size - 1 - (recovered.len() % block_size);
        let mut prefix_in = pad_attacker.clone();
        prefix_in.extend(std::iter::repeat(b'A').take(pad_extra));
        let target_ct = oracle(&prefix_in);
        let tgt = &target_ct[block_idx * block_size..(block_idx + 1) * block_size];
        let mut found = false;
        for b in 0..=255u8 {
            let mut cand = prefix_in.clone();
            cand.extend_from_slice(&recovered);
            cand.push(b);
            let cand_ct = oracle(&cand);
            let cb = &cand_ct[block_idx * block_size..(block_idx + 1) * block_size];
            if cb == tgt {
                recovered.push(b);
                found = true;
                break;
            }
        }
        if !found {
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
    fn recovers_secret_with_prefix() {
        assert!(super::run().success);
    }
}
