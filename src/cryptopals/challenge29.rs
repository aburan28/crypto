//! # Challenge 29 — SHA-1 keyed-MAC length extension
//!
//! Given `(message, sha1_mac(key, message))` but no `key`,
//! forge a MAC for `message || glue || suffix` where `glue` is the
//! padding the original would have appended.  Works because SHA-1's
//! Merkle-Damgård state at the end of the original message is
//! recoverable from the MAC itself.

use crate::cryptopals::challenge28::sha1_mac;
use crate::cryptopals::Report;

fn sha1_with_state(msg: &[u8], state: [u32; 5], prior_bytes: u64) -> [u8; 20] {
    // Resume SHA-1 from a known state, pretending we already
    // processed `prior_bytes` bytes.
    let mut padded = msg.to_vec();
    let total_bits = (prior_bytes + msg.len() as u64).wrapping_mul(8);
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&total_bits.to_be_bytes());
    let mut h0 = state[0];
    let mut h1 = state[1];
    let mut h2 = state[2];
    let mut h3 = state[3];
    let mut h4 = state[4];
    for chunk in padded.chunks_exact(64) {
        let mut w = [0u32; 80];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([
                chunk[4 * i],
                chunk[4 * i + 1],
                chunk[4 * i + 2],
                chunk[4 * i + 3],
            ]);
        }
        for i in 16..80 {
            w[i] = (w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16]).rotate_left(1);
        }
        let (mut a, mut b, mut c, mut d, mut e) = (h0, h1, h2, h3, h4);
        for i in 0..80 {
            let (f, k) = match i {
                0..=19 => ((b & c) | ((!b) & d), 0x5A82_7999u32),
                20..=39 => (b ^ c ^ d, 0x6ED9_EBA1),
                40..=59 => ((b & c) | (b & d) | (c & d), 0x8F1B_BCDC),
                _ => (b ^ c ^ d, 0xCA62_C1D6),
            };
            let temp = a
                .rotate_left(5)
                .wrapping_add(f)
                .wrapping_add(e)
                .wrapping_add(k)
                .wrapping_add(w[i]);
            e = d;
            d = c;
            c = b.rotate_left(30);
            b = a;
            a = temp;
        }
        h0 = h0.wrapping_add(a);
        h1 = h1.wrapping_add(b);
        h2 = h2.wrapping_add(c);
        h3 = h3.wrapping_add(d);
        h4 = h4.wrapping_add(e);
    }
    let mut out = [0u8; 20];
    out[0..4].copy_from_slice(&h0.to_be_bytes());
    out[4..8].copy_from_slice(&h1.to_be_bytes());
    out[8..12].copy_from_slice(&h2.to_be_bytes());
    out[12..16].copy_from_slice(&h3.to_be_bytes());
    out[16..20].copy_from_slice(&h4.to_be_bytes());
    out
}

fn glue_padding(len: usize) -> Vec<u8> {
    let mut g = vec![0x80u8];
    while (len + g.len()) % 64 != 56 {
        g.push(0);
    }
    let bit_len = (len as u64).wrapping_mul(8);
    g.extend_from_slice(&bit_len.to_be_bytes());
    g
}

pub fn run() -> Report {
    let mut r = Report::new(29, "SHA-1 length extension");
    let key = b"YELLOW SUBMARINE";
    let msg: &[u8] = b"comment1=cooking%20MCs;userdata=foo;comment2=%20like%20a%20pound%20of%20bacon";
    let tag = sha1_mac(key, msg);
    let suffix = b";admin=true";
    // Attack: try keys of plausible length (we don't know it).
    // Cryptopals' is 16, we'll just use that.
    let assumed_key_len = key.len();
    let glue = glue_padding(assumed_key_len + msg.len());
    let mut forged_msg = msg.to_vec();
    forged_msg.extend_from_slice(&glue);
    forged_msg.extend_from_slice(suffix);

    let mut state = [0u32; 5];
    for i in 0..5 {
        state[i] = u32::from_be_bytes(tag[4 * i..4 * i + 4].try_into().unwrap());
    }
    let prior = (assumed_key_len + msg.len() + glue.len()) as u64;
    let forged_tag = sha1_with_state(suffix, state, prior);

    let oracle_tag = sha1_mac(key, &forged_msg);
    r.line(format!("forged_tag = {}", hex::encode(forged_tag)));
    r.line(format!("oracle_tag = {}", hex::encode(oracle_tag)));
    assert_eq!(forged_tag, oracle_tag);
    assert!(forged_msg.ends_with(b";admin=true"));
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn extends_known_mac() {
        assert!(super::run().success);
    }
}
