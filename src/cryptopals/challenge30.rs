//! # Challenge 30 — MD4 length extension
//!
//! Same shape as #29 but with MD4 (little-endian state instead of
//! big-endian).  The repo already ships an `md4_compress` we can
//! reuse.

use crate::cryptopals::Report;
use crate::hash::md4::md4_compress;

pub fn md4_mac(key: &[u8], msg: &[u8]) -> [u8; 16] {
    let mut buf = key.to_vec();
    buf.extend_from_slice(msg);
    crate::hash::md4(&buf)
}

fn md4_padding(len: usize) -> Vec<u8> {
    let mut g = vec![0x80u8];
    while (len + g.len()) % 64 != 56 {
        g.push(0);
    }
    let bit_len = (len as u64).wrapping_mul(8);
    g.extend_from_slice(&bit_len.to_le_bytes());
    g
}

fn md4_from_state(state_in: &[u32; 4], msg: &[u8], prior_bytes: u64) -> [u8; 16] {
    let mut padded = msg.to_vec();
    let total = (prior_bytes + msg.len() as u64).wrapping_mul(8);
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&total.to_le_bytes());
    let mut s = *state_in;
    for chunk in padded.chunks_exact(64) {
        let mut block = [0u8; 64];
        block.copy_from_slice(chunk);
        md4_compress(&mut s, &block);
    }
    let mut out = [0u8; 16];
    for i in 0..4 {
        out[4 * i..4 * i + 4].copy_from_slice(&s[i].to_le_bytes());
    }
    out
}

pub fn run() -> Report {
    let mut r = Report::new(30, "MD4 length extension");
    let key = b"YELLOW SUBMARINE";
    let msg: &[u8] = b"comment1=cooking%20MCs;userdata=foo;comment2=%20like%20a%20pound%20of%20bacon";
    let tag = md4_mac(key, msg);
    let suffix = b";admin=true";
    let kl = key.len();
    let glue = md4_padding(kl + msg.len());
    let mut forged_msg = msg.to_vec();
    forged_msg.extend_from_slice(&glue);
    forged_msg.extend_from_slice(suffix);
    let mut state = [0u32; 4];
    for i in 0..4 {
        state[i] = u32::from_le_bytes(tag[4 * i..4 * i + 4].try_into().unwrap());
    }
    let prior = (kl + msg.len() + glue.len()) as u64;
    let forged_tag = md4_from_state(&state, suffix, prior);
    let oracle_tag = md4_mac(key, &forged_msg);
    r.line(format!("forged_tag = {}", hex::encode(forged_tag)));
    r.line(format!("oracle_tag = {}", hex::encode(oracle_tag)));
    assert_eq!(forged_tag, oracle_tag);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn extends_md4_mac() {
        assert!(super::run().success);
    }
}
