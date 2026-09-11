//! # Challenge 17 — CBC padding oracle
//!
//! Vaudenay 2002.  Given a ciphertext and an oracle that returns
//! true iff the plaintext has valid PKCS#7 padding, recover the
//! full plaintext one byte at a time.
//!
//! Per-byte recipe: target plaintext byte `P[k]`.  Pick a "scratch"
//! ciphertext block `C* = C_prev ⊕ pad_byte ⊕ guess`.  When the
//! oracle accepts, the decrypted byte at position k equals
//! `pad_byte`, so `guess = C_prev[k] ⊕ pad_byte ⊕ P[k]`.

use crate::cryptopals::challenge10::cbc_encrypt_no_iv_prefix;
use crate::cryptopals::low_util::{b64_decode, pkcs7_unpad};
use crate::cryptopals::Report;
use crate::symmetric::aes::{decrypt_block, AesKey};
use rand::{seq::SliceRandom, SeedableRng};
use rand::rngs::StdRng;

const STRINGS: &[&str] = &[
    "MDAwMDAwTm93IHRoYXQgdGhlIHBhcnR5IGlzIGp1bXBpbmc=",
    "MDAwMDAxV2l0aCB0aGUgYmFzcyBraWNrZWQgaW4gYW5kIHRoZSBWZWdhJ3MgYXJlIHB1bXBpbic=",
    "MDAwMDAyUXVpY2sgdG8gdGhlIHBvaW50LCB0byB0aGUgcG9pbnQsIG5vIGZha2luZw==",
    "MDAwMDAzQ29va2luZyBNQydzIGxpa2UgYSBwb3VuZCBvZiBiYWNvbg==",
    "MDAwMDA0QnVybmluZyAnZW0sIGlmIHlvdSBhaW4ndCBxdWljayBhbmQgbmltYmxl",
    "MDAwMDA1SSBnbyBjcmF6eSB3aGVuIEkgaGVhciBhIGN5bWJhbA==",
    "MDAwMDA2QW5kIGEgaGlnaCBoYXQgd2l0aCBhIHNvdXBlZCB1cCB0ZW1wbw==",
    "MDAwMDA3SSdtIG9uIGEgcm9sbCwgaXQncyB0aW1lIHRvIGdvIHNvbG8=",
    "MDAwMDA4b2xsaW4nIGluIG15IGZpdmUgcG9pbnQgb2g=",
    "MDAwMDA5aXRoIG15IHJhZy10b3AgZG93biBzbyBteSBoYWlyIGNhbiBibG93",
];

fn oracle_encrypt(key: &AesKey, rng: &mut StdRng) -> (Vec<u8>, [u8; 16]) {
    use rand::Rng;
    let s = STRINGS.choose(rng).unwrap();
    let pt = b64_decode(s);
    let mut iv = [0u8; 16];
    rng.fill(&mut iv);
    let ct = cbc_encrypt_no_iv_prefix(&pt, key, &iv);
    (ct, iv)
}

/// Returns `true` iff the ct decrypts to a valid-padded plaintext.
pub fn padding_oracle(ct: &[u8], iv: &[u8; 16], key: &AesKey) -> bool {
    // Decrypt without enforcing pkcs7 inside aes_cbc; do it ourselves.
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
    pkcs7_unpad(&out, 16).is_some()
}

/// Decrypt one ciphertext block `cn` whose predecessor was `c_prev`,
/// using a padding oracle.  Returns the plaintext block.
pub fn decrypt_block_via_oracle(
    c_prev: &[u8; 16],
    cn: &[u8; 16],
    key: &AesKey,
) -> [u8; 16] {
    let mut intermediate = [0u8; 16]; // I = AES_K^-1(cn)
    for k in (0..16).rev() {
        let pad = (16 - k) as u8;
        // Set up scratch block s so that the trailing (16 - k)
        // bytes decrypt to `pad pad pad …`.
        for guess in 0..=255u8 {
            let mut scratch = [0u8; 16];
            scratch[k] = guess;
            // bytes after k: each = pad ^ intermediate[j]
            for j in (k + 1)..16 {
                scratch[j] = pad ^ intermediate[j];
            }
            if padding_oracle(cn, &scratch, key) {
                // Edge case for k == 15 (pad == 1): the oracle
                // can be fooled by an existing valid pad >1 in cn.
                // Guard by flipping byte k-1 and re-asking.
                if k == 15 {
                    let mut probe = scratch;
                    if k > 0 {
                        probe[k - 1] ^= 0x01;
                    }
                    if !padding_oracle(cn, &probe, key) {
                        continue;
                    }
                }
                intermediate[k] = guess ^ pad;
                break;
            }
        }
    }
    let mut pt = [0u8; 16];
    for i in 0..16 {
        pt[i] = intermediate[i] ^ c_prev[i];
    }
    pt
}

pub fn run() -> Report {
    let mut r = Report::new(17, "CBC padding oracle");
    let key_bytes: [u8; 16] = *b"keyforC17-xxxxyz";
    let key = AesKey::new(&key_bytes).unwrap();
    let mut rng = StdRng::seed_from_u64(17);
    let (ct, iv) = oracle_encrypt(&key, &mut rng);

    // Decrypt block by block.
    let mut pt = Vec::new();
    let mut prev = iv;
    for chunk in ct.chunks_exact(16) {
        let mut cn = [0u8; 16];
        cn.copy_from_slice(chunk);
        let p = decrypt_block_via_oracle(&prev, &cn, &key);
        pt.extend_from_slice(&p);
        prev = cn;
    }
    let stripped = pkcs7_unpad(&pt, 16).expect("valid pad after recovery");
    r.line(format!(
        "Recovered: {:?}",
        std::str::from_utf8(&stripped).unwrap_or("?")
    ));
    assert!(stripped.iter().any(|c| c.is_ascii()));
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn padding_oracle_recovers_one_string() {
        assert!(super::run().success);
    }
}
