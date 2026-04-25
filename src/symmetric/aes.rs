//! AES (Advanced Encryption Standard) — FIPS 197, implemented from scratch.
//!
//! Supports AES-128 and AES-256 block encryption/decryption.
//! CTR streaming mode is built on top of the block cipher.
//! GCM authenticated encryption is provided via GHASH over GF(2¹²⁸).
//!
//! # Round structure
//! Each round applies four transformations to a 4×4 byte state:
//!   SubBytes  — non-linear byte substitution via S-box (built from GF(2⁸) inverse + affine map).
//!   ShiftRows — cyclic rotation of each row by its index.
//!   MixColumns — GF(2⁸) matrix multiply of each column (skipped in the final round).
//!   AddRoundKey — XOR with the round key derived by the key schedule.

// ── S-box / inverse S-box ────────────────────────────────────────────────────

#[rustfmt::skip]
const SBOX: [u8; 256] = [
    0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
    0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
    0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
    0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
    0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
    0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
    0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
    0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
    0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
    0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
    0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
    0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
    0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
    0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
    0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
    0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16,
];

#[rustfmt::skip]
const INV_SBOX: [u8; 256] = [
    0x52,0x09,0x6a,0xd5,0x30,0x36,0xa5,0x38,0xbf,0x40,0xa3,0x9e,0x81,0xf3,0xd7,0xfb,
    0x7c,0xe3,0x39,0x82,0x9b,0x2f,0xff,0x87,0x34,0x8e,0x43,0x44,0xc4,0xde,0xe9,0xcb,
    0x54,0x7b,0x94,0x32,0xa6,0xc2,0x23,0x3d,0xee,0x4c,0x95,0x0b,0x42,0xfa,0xc3,0x4e,
    0x08,0x2e,0xa1,0x66,0x28,0xd9,0x24,0xb2,0x76,0x5b,0xa2,0x49,0x6d,0x8b,0xd1,0x25,
    0x72,0xf8,0xf6,0x64,0x86,0x68,0x98,0x16,0xd4,0xa4,0x5c,0xcc,0x5d,0x65,0xb6,0x92,
    0x6c,0x70,0x48,0x50,0xfd,0xed,0xb9,0xda,0x5e,0x15,0x46,0x57,0xa7,0x8d,0x9d,0x84,
    0x90,0xd8,0xab,0x00,0x8c,0xbc,0xd3,0x0a,0xf7,0xe4,0x58,0x05,0xb8,0xb3,0x45,0x06,
    0xd0,0x2c,0x1e,0x8f,0xca,0x3f,0x0f,0x02,0xc1,0xaf,0xbd,0x03,0x01,0x13,0x8a,0x6b,
    0x3a,0x91,0x11,0x41,0x4f,0x67,0xdc,0xea,0x97,0xf2,0xcf,0xce,0xf0,0xb4,0xe6,0x73,
    0x96,0xac,0x74,0x22,0xe7,0xad,0x35,0x85,0xe2,0xf9,0x37,0xe8,0x1c,0x75,0xdf,0x6e,
    0x47,0xf1,0x1a,0x71,0x1d,0x29,0xc5,0x89,0x6f,0xb7,0x62,0x0e,0xaa,0x18,0xbe,0x1b,
    0xfc,0x56,0x3e,0x4b,0xc6,0xd2,0x79,0x20,0x9a,0xdb,0xc0,0xfe,0x78,0xcd,0x5a,0xf4,
    0x1f,0xdd,0xa8,0x33,0x88,0x07,0xc7,0x31,0xb1,0x12,0x10,0x59,0x27,0x80,0xec,0x5f,
    0x60,0x51,0x7f,0xa9,0x19,0xb5,0x4a,0x0d,0x2d,0xe5,0x7a,0x9f,0x93,0xc9,0x9c,0xef,
    0xa0,0xe0,0x3b,0x4d,0xae,0x2a,0xf5,0xb0,0xc8,0xeb,0xbb,0x3c,0x83,0x53,0x99,0x61,
    0x17,0x2b,0x04,0x7e,0xba,0x77,0xd6,0x26,0xe1,0x69,0x14,0x63,0x55,0x21,0x0c,0x7d,
];

// Round constant (Rcon) for key schedule — powers of 2 in GF(2⁸).
const RCON: [u8; 11] = [0x00, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36];

// ── GF(2⁸) arithmetic ────────────────────────────────────────────────────────

/// Multiply two bytes in GF(2⁸) using the peasant's algorithm.
fn gmul(mut a: u8, mut b: u8) -> u8 {
    let mut p = 0u8;
    for _ in 0..8 {
        if b & 1 != 0 { p ^= a; }
        let hi = a & 0x80;
        a <<= 1;
        if hi != 0 { a ^= 0x1b; }
        b >>= 1;
    }
    p
}

// ── AES key types ─────────────────────────────────────────────────────────────

/// AES key variants.
pub enum AesKey {
    Aes128([u8; 16]),
    Aes256([u8; 32]),
}

impl AesKey {
    pub fn key_len(&self) -> usize {
        match self { AesKey::Aes128(_) => 16, AesKey::Aes256(_) => 32 }
    }
    pub fn nk(&self) -> usize { self.key_len() / 4 }
    pub fn nr(&self) -> usize { self.nk() + 6 }
    pub fn as_bytes(&self) -> &[u8] {
        match self { AesKey::Aes128(k) => k, AesKey::Aes256(k) => k }
    }
}

// ── Key expansion ─────────────────────────────────────────────────────────────

/// Expand a key into the full round key schedule (Nb*(Nr+1) words = 4*(Nr+1) 32-bit words).
pub fn key_expansion(key: &AesKey) -> Vec<[u8; 4]> {
    let nk = key.nk();
    let nr = key.nr();
    let total = 4 * (nr + 1);
    let mut w: Vec<[u8; 4]> = Vec::with_capacity(total);
    let k = key.as_bytes();

    for i in 0..nk {
        w.push([k[4*i], k[4*i+1], k[4*i+2], k[4*i+3]]);
    }

    for i in nk..total {
        let mut temp = w[i - 1];
        if i % nk == 0 {
            temp = sub_word(rot_word(temp));
            temp[0] ^= RCON[i / nk];
        } else if nk > 6 && i % nk == 4 {
            temp = sub_word(temp);
        }
        let prev = w[i - nk];
        w.push([prev[0]^temp[0], prev[1]^temp[1], prev[2]^temp[2], prev[3]^temp[3]]);
    }
    w
}

fn rot_word(w: [u8; 4]) -> [u8; 4] { [w[1], w[2], w[3], w[0]] }
fn sub_word(w: [u8; 4]) -> [u8; 4] {
    [SBOX[w[0] as usize], SBOX[w[1] as usize], SBOX[w[2] as usize], SBOX[w[3] as usize]]
}

// ── State helpers ─────────────────────────────────────────────────────────────

/// AES state: 4×4 bytes, stored column-major.
type State = [[u8; 4]; 4];

fn bytes_to_state(block: &[u8]) -> State {
    let mut s = [[0u8; 4]; 4];
    for c in 0..4 { for r in 0..4 { s[c][r] = block[c * 4 + r]; } }
    s
}

fn state_to_bytes(s: &State) -> [u8; 16] {
    let mut out = [0u8; 16];
    for c in 0..4 { for r in 0..4 { out[c * 4 + r] = s[c][r]; } }
    out
}

fn add_round_key(s: &mut State, round_key: &[[u8; 4]]) {
    for c in 0..4 {
        for r in 0..4 {
            s[c][r] ^= round_key[c][r];
        }
    }
}

// ── Forward transformations ───────────────────────────────────────────────────

fn sub_bytes(s: &mut State) {
    for col in s.iter_mut() { for b in col.iter_mut() { *b = SBOX[*b as usize]; } }
}

fn shift_rows(s: &mut State) {
    // Row 0: no shift. Row 1: shift left 1. Row 2: shift left 2. Row 3: shift left 3.
    let r1 = [s[1][1], s[2][1], s[3][1], s[0][1]];
    let r2 = [s[2][2], s[3][2], s[0][2], s[1][2]];
    let r3 = [s[3][3], s[0][3], s[1][3], s[2][3]];
    for c in 0..4 { s[c][1] = r1[c]; s[c][2] = r2[c]; s[c][3] = r3[c]; }
}

fn mix_columns(s: &mut State) {
    for c in 0..4 {
        let a = s[c];
        s[c][0] = gmul(0x02, a[0]) ^ gmul(0x03, a[1]) ^ a[2]             ^ a[3];
        s[c][1] = a[0]             ^ gmul(0x02, a[1]) ^ gmul(0x03, a[2]) ^ a[3];
        s[c][2] = a[0]             ^ a[1]             ^ gmul(0x02, a[2]) ^ gmul(0x03, a[3]);
        s[c][3] = gmul(0x03, a[0]) ^ a[1]             ^ a[2]             ^ gmul(0x02, a[3]);
    }
}

// ── Inverse transformations ───────────────────────────────────────────────────

fn inv_sub_bytes(s: &mut State) {
    for col in s.iter_mut() { for b in col.iter_mut() { *b = INV_SBOX[*b as usize]; } }
}

fn inv_shift_rows(s: &mut State) {
    let r1 = [s[3][1], s[0][1], s[1][1], s[2][1]];
    let r2 = [s[2][2], s[3][2], s[0][2], s[1][2]];
    let r3 = [s[1][3], s[2][3], s[3][3], s[0][3]];
    for c in 0..4 { s[c][1] = r1[c]; s[c][2] = r2[c]; s[c][3] = r3[c]; }
}

fn inv_mix_columns(s: &mut State) {
    for c in 0..4 {
        let a = s[c];
        s[c][0] = gmul(0x0e,a[0])^gmul(0x0b,a[1])^gmul(0x0d,a[2])^gmul(0x09,a[3]);
        s[c][1] = gmul(0x09,a[0])^gmul(0x0e,a[1])^gmul(0x0b,a[2])^gmul(0x0d,a[3]);
        s[c][2] = gmul(0x0d,a[0])^gmul(0x09,a[1])^gmul(0x0e,a[2])^gmul(0x0b,a[3]);
        s[c][3] = gmul(0x0b,a[0])^gmul(0x0d,a[1])^gmul(0x09,a[2])^gmul(0x0e,a[3]);
    }
}

// ── Block encrypt / decrypt ───────────────────────────────────────────────────

/// Encrypt a single 16-byte block with `key`.
pub fn encrypt_block(block: &[u8; 16], key: &AesKey) -> [u8; 16] {
    let nr = key.nr();
    let w = key_expansion(key);
    let mut s = bytes_to_state(block);

    add_round_key(&mut s, &w[..4]);

    for round in 1..nr {
        sub_bytes(&mut s);
        shift_rows(&mut s);
        mix_columns(&mut s);
        add_round_key(&mut s, &w[round * 4..(round + 1) * 4]);
    }

    // Final round (no MixColumns)
    sub_bytes(&mut s);
    shift_rows(&mut s);
    add_round_key(&mut s, &w[nr * 4..(nr + 1) * 4]);

    state_to_bytes(&s)
}

/// Decrypt a single 16-byte block.
pub fn decrypt_block(block: &[u8; 16], key: &AesKey) -> [u8; 16] {
    let nr = key.nr();
    let w = key_expansion(key);
    let mut s = bytes_to_state(block);

    add_round_key(&mut s, &w[nr * 4..(nr + 1) * 4]);

    for round in (1..nr).rev() {
        inv_shift_rows(&mut s);
        inv_sub_bytes(&mut s);
        add_round_key(&mut s, &w[round * 4..(round + 1) * 4]);
        inv_mix_columns(&mut s);
    }

    inv_shift_rows(&mut s);
    inv_sub_bytes(&mut s);
    add_round_key(&mut s, &w[..4]);

    state_to_bytes(&s)
}

// ── CTR mode ─────────────────────────────────────────────────────────────────

/// AES-CTR encryption/decryption (symmetric operation).
///
/// Counter block layout: `nonce` (12 bytes) || counter (4 bytes, big-endian).
/// The counter starts at 1 (matching RFC 3686 / GCM conventions).
pub fn aes_ctr(data: &[u8], key: &AesKey, nonce: &[u8; 12]) -> Vec<u8> {
    let mut out = data.to_vec();
    let mut counter = 1u32;

    for chunk in out.chunks_mut(16) {
        let mut ctr_block = [0u8; 16];
        ctr_block[..12].copy_from_slice(nonce);
        ctr_block[12..].copy_from_slice(&counter.to_be_bytes());

        let keystream = encrypt_block(&ctr_block, key);
        for (b, k) in chunk.iter_mut().zip(keystream.iter()) {
            *b ^= k;
        }
        counter = counter.wrapping_add(1);
    }
    out
}

// ── GCM (GHASH + CTR) ─────────────────────────────────────────────────────────

/// Multiply two 128-bit values in GF(2¹²⁸) with the GCM modulus.
/// GCM uses the "bit-reflected" representation; the modulus is x¹²⁸+x⁷+x²+x+1.
fn gcm_mult(x: &[u8; 16], y: &[u8; 16]) -> [u8; 16] {
    let mut z = [0u8; 16];
    let mut v = *y;

    for i in 0..16 {
        for bit in (0..8).rev() {
            if (x[i] >> bit) & 1 == 1 {
                for j in 0..16 { z[j] ^= v[j]; }
            }
            let lsb = v[15] & 1;
            // Right-shift v by 1 bit
            for j in (1..16).rev() { v[j] = (v[j] >> 1) | (v[j-1] << 7); }
            v[0] >>= 1;
            // Reduce mod the GCM polynomial if the shifted-out bit was 1
            if lsb == 1 { v[0] ^= 0xe1; }
        }
    }
    z
}

/// GHASH: authenticate `data` using hash key `h`, returning a 16-byte tag.
fn ghash(h: &[u8; 16], aad: &[u8], ciphertext: &[u8]) -> [u8; 16] {
    let mut y = [0u8; 16];

    // Process AAD padded to 16-byte boundary
    let padded_aad = pad_to_block(aad);
    for block in padded_aad.chunks(16) {
        xor_into(&mut y, block);
        y = gcm_mult(&y, h);
    }

    // Process ciphertext
    let padded_ct = pad_to_block(ciphertext);
    for block in padded_ct.chunks(16) {
        xor_into(&mut y, block);
        y = gcm_mult(&y, h);
    }

    // Length block: len(AAD)||len(ciphertext) in bits, each 64-bit big-endian
    let mut len_block = [0u8; 16];
    len_block[..8].copy_from_slice(&((aad.len() as u64 * 8).to_be_bytes()));
    len_block[8..].copy_from_slice(&((ciphertext.len() as u64 * 8).to_be_bytes()));
    xor_into(&mut y, &len_block);
    y = gcm_mult(&y, h);

    y
}

fn pad_to_block(data: &[u8]) -> Vec<u8> {
    let mut v = data.to_vec();
    let rem = v.len() % 16;
    if rem != 0 { v.resize(v.len() + (16 - rem), 0); }
    if v.is_empty() { v.resize(16, 0); }
    v
}

fn xor_into(dst: &mut [u8; 16], src: &[u8]) {
    for (a, b) in dst.iter_mut().zip(src.iter()) { *a ^= b; }
}

/// AES-GCM encrypt. Returns `ciphertext || 16-byte tag`.
///
/// `key`:   AES-128 or AES-256 key.
/// `nonce`: 12-byte IV (use a unique value for each message).
/// `aad`:   Additional Authenticated Data (not encrypted, but authenticated).
pub fn aes_gcm_encrypt(plaintext: &[u8], key: &AesKey, nonce: &[u8; 12], aad: &[u8]) -> Vec<u8> {
    // Derive GHASH key H = E(K, 0¹²⁸)
    let h: [u8; 16] = encrypt_block(&[0u8; 16], key);

    // Encrypt
    let ciphertext = aes_ctr(plaintext, key, nonce);

    // Compute tag: GHASH(H, AAD, CT) XOR E(K, J₀)
    let mut tag = ghash(&h, aad, &ciphertext);
    let j0 = encrypt_j0(key, nonce);
    for (t, j) in tag.iter_mut().zip(j0.iter()) { *t ^= j; }

    let mut out = ciphertext;
    out.extend_from_slice(&tag);
    out
}

/// AES-GCM decrypt. Returns `Ok(plaintext)` if the tag verifies, `Err(())` otherwise.
pub fn aes_gcm_decrypt(
    ciphertext_and_tag: &[u8],
    key: &AesKey,
    nonce: &[u8; 12],
    aad: &[u8],
) -> Result<Vec<u8>, ()> {
    if ciphertext_and_tag.len() < 16 { return Err(()); }
    let (ciphertext, tag_bytes) = ciphertext_and_tag.split_at(ciphertext_and_tag.len() - 16);

    let h: [u8; 16] = encrypt_block(&[0u8; 16], key);
    let mut expected_tag = ghash(&h, aad, ciphertext);
    let j0 = encrypt_j0(key, nonce);
    for (t, j) in expected_tag.iter_mut().zip(j0.iter()) { *t ^= j; }

    // Constant-time tag comparison to prevent timing attacks
    let tag_ok = expected_tag.iter().zip(tag_bytes.iter()).fold(0u8, |acc, (a, b)| acc | (a ^ b)) == 0;
    if !tag_ok { return Err(()); }

    Ok(aes_ctr(ciphertext, key, nonce))
}

/// Encrypt the J₀ counter block (nonce || 0x00000001) for the authentication tag.
fn encrypt_j0(key: &AesKey, nonce: &[u8; 12]) -> [u8; 16] {
    let mut j0 = [0u8; 16];
    j0[..12].copy_from_slice(nonce);
    j0[15] = 1;
    encrypt_block(&j0, key)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aes128_roundtrip() {
        let key = AesKey::Aes128([0u8; 16]);
        let plaintext = [0u8; 16];
        let ct = encrypt_block(&plaintext, &key);
        let pt = decrypt_block(&ct, &key);
        assert_eq!(pt, plaintext);
    }

    #[test]
    fn aes128_known_answer() {
        // FIPS 197 Appendix B: key=2b7e..., pt=3243...
        let key = AesKey::Aes128([
            0x2b,0x7e,0x15,0x16,0x28,0xae,0xd2,0xa6,
            0xab,0xf7,0x15,0x88,0x09,0xcf,0x4f,0x3c,
        ]);
        let pt = [
            0x32,0x43,0xf6,0xa8,0x88,0x5a,0x30,0x8d,
            0x31,0x31,0x98,0xa2,0xe0,0x37,0x07,0x34,
        ];
        let ct = encrypt_block(&pt, &key);
        let expected = [
            0x39,0x25,0x84,0x1d,0x02,0xdc,0x09,0xfb,
            0xdc,0x11,0x85,0x97,0x19,0x6a,0x0b,0x32,
        ];
        assert_eq!(ct, expected);
    }

    #[test]
    fn aes256_roundtrip() {
        let key = AesKey::Aes256([0u8; 32]);
        let plaintext = [0x42u8; 16];
        let ct = encrypt_block(&plaintext, &key);
        let pt = decrypt_block(&ct, &key);
        assert_eq!(pt, plaintext);
    }

    #[test]
    fn ctr_roundtrip() {
        let key = AesKey::Aes128([0x2b,0x7e,0x15,0x16,0x28,0xae,0xd2,0xa6,
                                   0xab,0xf7,0x15,0x88,0x09,0xcf,0x4f,0x3c]);
        let nonce = [0u8; 12];
        let msg = b"Hello, AES-CTR world!";
        let ct = aes_ctr(msg, &key, &nonce);
        let pt = aes_ctr(&ct, &key, &nonce);
        assert_eq!(&pt, msg);
    }

    #[test]
    fn gcm_roundtrip() {
        let key = AesKey::Aes128([0u8; 16]);
        let nonce = [0u8; 12];
        let pt = b"secret message!";
        let aad = b"header";
        let ct = aes_gcm_encrypt(pt, &key, &nonce, aad);
        let recovered = aes_gcm_decrypt(&ct, &key, &nonce, aad).unwrap();
        assert_eq!(recovered, pt);
    }

    #[test]
    fn gcm_tag_failure() {
        let key = AesKey::Aes128([0u8; 16]);
        let nonce = [0u8; 12];
        let mut ct = aes_gcm_encrypt(b"msg", &key, &nonce, b"aad");
        // Tamper with the ciphertext
        ct[0] ^= 0xff;
        assert!(aes_gcm_decrypt(&ct, &key, &nonce, b"aad").is_err());
    }
}
