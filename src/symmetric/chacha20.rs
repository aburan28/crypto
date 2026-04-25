//! ChaCha20-Poly1305 AEAD — RFC 8439, implemented from scratch.
//!
//! # ChaCha20 stream cipher
//! State: 16 32-bit little-endian words.
//!   words  0- 3: "expa nd 3 2-by te k" (ASCII constants)
//!   words  4-11: 256-bit key
//!   word  12: 32-bit block counter
//!   words 13-15: 96-bit nonce
//!
//! Each block applies 20 rounds of "quarter rounds" (ARX operations),
//! then adds the original state to produce 64 bytes of keystream.
//!
//! # Poly1305 MAC
//! One-time authenticator over GF(2¹³⁰ - 5).
//! Key: (r, s) derived from the first 32 bytes of ChaCha20 output (counter=0).
//! r is clamped by zeroing specific bits; s is used as-is.
//! The tag is: ((clamp(r)·m₁ + clamp(r)·m₂ + … + s) mod 2¹³⁰-5) mod 2¹²⁸

use num_bigint::BigUint;
use num_traits::{One, Zero};

// ── ChaCha20 ─────────────────────────────────────────────────────────────────

const SIGMA: [u32; 4] = [0x61707865, 0x3320646e, 0x79622d32, 0x6b206574];

/// Apply one ChaCha20 quarter round to four words.
macro_rules! qr {
    ($a:expr, $b:expr, $c:expr, $d:expr) => {
        $a = $a.wrapping_add($b); $d ^= $a; $d = $d.rotate_left(16);
        $c = $c.wrapping_add($d); $b ^= $c; $b = $b.rotate_left(12);
        $a = $a.wrapping_add($b); $d ^= $a; $d = $d.rotate_left(8);
        $c = $c.wrapping_add($d); $b ^= $c; $b = $b.rotate_left(7);
    };
}

/// Generate one 64-byte ChaCha20 block at the given counter value.
pub fn chacha20_block(key: &[u8; 32], nonce: &[u8; 12], counter: u32) -> [u8; 64] {
    let mut state = [0u32; 16];
    state[..4].copy_from_slice(&SIGMA);

    for i in 0..8 {
        state[4 + i] = u32::from_le_bytes(key[i*4..i*4+4].try_into().unwrap());
    }
    state[12] = counter;
    for i in 0..3 {
        state[13 + i] = u32::from_le_bytes(nonce[i*4..i*4+4].try_into().unwrap());
    }

    let mut working = state;

    // 20 rounds = 10 column rounds + 10 diagonal rounds
    for _ in 0..10 {
        // Column rounds
        qr!(working[0], working[4], working[8],  working[12]);
        qr!(working[1], working[5], working[9],  working[13]);
        qr!(working[2], working[6], working[10], working[14]);
        qr!(working[3], working[7], working[11], working[15]);
        // Diagonal rounds
        qr!(working[0], working[5], working[10], working[15]);
        qr!(working[1], working[6], working[11], working[12]);
        qr!(working[2], working[7], working[8],  working[13]);
        qr!(working[3], working[4], working[9],  working[14]);
    }

    let mut out = [0u8; 64];
    for (i, w) in working.iter().enumerate() {
        let w_final = w.wrapping_add(state[i]);
        out[i*4..i*4+4].copy_from_slice(&w_final.to_le_bytes());
    }
    out
}

/// Encrypt or decrypt `data` using ChaCha20 (XOR with keystream, counter starts at 1).
pub fn chacha20_xor(data: &[u8], key: &[u8; 32], nonce: &[u8; 12]) -> Vec<u8> {
    let mut out = data.to_vec();
    let mut counter = 1u32;

    for chunk in out.chunks_mut(64) {
        let block = chacha20_block(key, nonce, counter);
        for (b, k) in chunk.iter_mut().zip(block.iter()) {
            *b ^= k;
        }
        counter = counter.wrapping_add(1);
    }
    out
}

// ── Poly1305 ─────────────────────────────────────────────────────────────────

fn p1305() -> BigUint {
    // 2^130 - 5
    (BigUint::one() << 130u32) - BigUint::from(5u32)
}

/// Clamp `r` per RFC 8439: zero out specific bits to reduce the key space.
fn clamp_r(r: &mut [u8; 16]) {
    r[3]  &= 0x0f;
    r[7]  &= 0x0f;
    r[11] &= 0x0f;
    r[15] &= 0x0f;
    r[4]  &= 0xfc;
    r[8]  &= 0xfc;
    r[12] &= 0xfc;
}

/// Compute a Poly1305 MAC over `msg` with a 32-byte one-time key (r || s).
pub fn poly1305(key: &[u8; 32], msg: &[u8]) -> [u8; 16] {
    let mut r_bytes: [u8; 16] = key[..16].try_into().unwrap();
    clamp_r(&mut r_bytes);

    let r = BigUint::from_bytes_le(&r_bytes);
    let s = BigUint::from_bytes_le(&key[16..]);
    let p = p1305();

    let mut acc = BigUint::zero();

    for block in msg.chunks(16) {
        let mut block_bytes = block.to_vec();
        block_bytes.push(0x01); // Append one bit above the block
        let n = BigUint::from_bytes_le(&block_bytes);
        acc = ((acc + n) * &r) % &p;
    }

    acc += s;
    // Take mod 2^128 (truncate to 16 bytes little-endian)
    let bytes = acc.to_bytes_le();
    let mut tag = [0u8; 16];
    let copy_len = bytes.len().min(16);
    tag[..copy_len].copy_from_slice(&bytes[..copy_len]);
    tag
}

// ── ChaCha20-Poly1305 AEAD ────────────────────────────────────────────────────

/// Derive the Poly1305 one-time key from the ChaCha20 key + nonce (counter=0).
fn poly1305_key_gen(key: &[u8; 32], nonce: &[u8; 12]) -> [u8; 32] {
    let block = chacha20_block(key, nonce, 0);
    block[..32].try_into().unwrap()
}

/// Build the Poly1305 input: PAD(AAD) || PAD(CT) || len(AAD) || len(CT).
fn poly1305_pad(aad: &[u8], ciphertext: &[u8]) -> Vec<u8> {
    let mut msg = Vec::new();
    msg.extend_from_slice(aad);
    pad16(&mut msg);
    msg.extend_from_slice(ciphertext);
    pad16(&mut msg);
    msg.extend_from_slice(&(aad.len() as u64).to_le_bytes());
    msg.extend_from_slice(&(ciphertext.len() as u64).to_le_bytes());
    msg
}

fn pad16(v: &mut Vec<u8>) {
    let rem = v.len() % 16;
    if rem != 0 { v.resize(v.len() + (16 - rem), 0); }
}

/// ChaCha20-Poly1305 encrypt. Returns `ciphertext || 16-byte tag`.
pub fn chacha20_poly1305_encrypt(
    plaintext: &[u8],
    key: &[u8; 32],
    nonce: &[u8; 12],
    aad: &[u8],
) -> Vec<u8> {
    let ciphertext = chacha20_xor(plaintext, key, nonce);
    let otk = poly1305_key_gen(key, nonce);
    let mac_data = poly1305_pad(aad, &ciphertext);
    let tag = poly1305(&otk, &mac_data);

    let mut out = ciphertext;
    out.extend_from_slice(&tag);
    out
}

/// ChaCha20-Poly1305 decrypt. Returns `Ok(plaintext)` or `Err(())` on tag failure.
pub fn chacha20_poly1305_decrypt(
    ciphertext_and_tag: &[u8],
    key: &[u8; 32],
    nonce: &[u8; 12],
    aad: &[u8],
) -> Result<Vec<u8>, ()> {
    if ciphertext_and_tag.len() < 16 { return Err(()); }
    let (ciphertext, tag_bytes) = ciphertext_and_tag.split_at(ciphertext_and_tag.len() - 16);

    let otk = poly1305_key_gen(key, nonce);
    let mac_data = poly1305_pad(aad, ciphertext);
    let expected_tag = poly1305(&otk, &mac_data);

    let ok = expected_tag.iter().zip(tag_bytes.iter()).fold(0u8, |a, (x, y)| a | (x ^ y)) == 0;
    if !ok { return Err(()); }

    Ok(chacha20_xor(ciphertext, key, nonce))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chacha20_block_rfc8439() {
        // RFC 8439 §2.1.2 test vector
        let key = [0u8; 32];
        let nonce = [0u8; 12];
        let block = chacha20_block(&key, &nonce, 0);
        // First 4 bytes of the output keystream
        assert_eq!(&block[..4], &[0x76, 0xb8, 0xe0, 0xad]);
    }

    #[test]
    fn chacha20_roundtrip() {
        let key = [0x42u8; 32];
        let nonce = [0u8; 12];
        let msg = b"Hello, ChaCha20!";
        let ct = chacha20_xor(msg, &key, &nonce);
        let pt = chacha20_xor(&ct, &key, &nonce);
        assert_eq!(&pt, msg);
    }

    #[test]
    fn chacha20_poly1305_roundtrip() {
        let key = [0x42u8; 32];
        let nonce = [0u8; 12];
        let msg = b"authenticated message";
        let aad = b"header data";
        let ct = chacha20_poly1305_encrypt(msg, &key, &nonce, aad);
        let pt = chacha20_poly1305_decrypt(&ct, &key, &nonce, aad).unwrap();
        assert_eq!(&pt, msg);
    }

    #[test]
    fn chacha20_poly1305_tamper_fails() {
        let key = [0x42u8; 32];
        let nonce = [0u8; 12];
        let mut ct = chacha20_poly1305_encrypt(b"secret", &key, &nonce, b"");
        ct[0] ^= 0xff;
        assert!(chacha20_poly1305_decrypt(&ct, &key, &nonce, b"").is_err());
    }
}
