//! # Challenge 28 — SHA-1 keyed MAC
//!
//! Compute `SHA1(key || message)` as a MAC.  Verify it detects
//! tampering.

use crate::cryptopals::Report;
use crate::hash::sha256::sha256;

/// Minimal SHA-1 implementation (just for this challenge — the
/// repo doesn't ship SHA-1).  100-line FIPS 180-1 transcription.
pub fn sha1(msg: &[u8]) -> [u8; 20] {
    let mut padded = msg.to_vec();
    let bit_len = (msg.len() as u64).wrapping_mul(8);
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());
    let mut h0: u32 = 0x6745_2301;
    let mut h1: u32 = 0xEFCD_AB89;
    let mut h2: u32 = 0x98BA_DCFE;
    let mut h3: u32 = 0x1032_5476;
    let mut h4: u32 = 0xC3D2_E1F0;
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

pub fn sha1_mac(key: &[u8], msg: &[u8]) -> [u8; 20] {
    let mut buf = key.to_vec();
    buf.extend_from_slice(msg);
    sha1(&buf)
}

pub fn run() -> Report {
    let mut r = Report::new(28, "SHA-1 keyed MAC");
    let _ = sha256;
    // RFC 3174 test vector: "abc" → a9993e36...d39
    let h = sha1(b"abc");
    r.line(format!("SHA-1(\"abc\") = {}", hex::encode(h)));
    assert_eq!(hex::encode(h), "a9993e364706816aba3e25717850c26c9cd0d89d");
    let key = b"YELLOW SUBMARINE";
    let m1 = b"hello";
    let m2 = b"hello!";
    let t1 = sha1_mac(key, m1);
    let t2 = sha1_mac(key, m2);
    assert_ne!(t1, t2);
    r.line("Different messages produce different MACs.");
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn sha1_known_vector() {
        assert!(super::run().success);
    }
}
