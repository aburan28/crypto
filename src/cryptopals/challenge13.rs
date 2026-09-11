//! # Challenge 13 — ECB cut-and-paste
//!
//! Encode `profile_for(email)` as `email=…&uid=10&role=user` and
//! ECB-encrypt.  By choosing the email carefully, the attacker can:
//!
//! 1. Position the literal bytes `admin\x0b\x0b\x0b\x0b\x0b…` (PKCS#7
//!    padding for an "admin" block) at the start of a block.
//! 2. Capture that ciphertext block.
//! 3. Build a *second* email whose `role=` field starts at a block
//!    boundary, capturing all the leading blocks of THAT
//!    ciphertext.
//! 4. Concatenate: the leading blocks from (3) + the "admin"
//!    block from (2) → an ECB ciphertext that decrypts to
//!    `email=…&uid=10&role=admin\x0b\x0b…`.

use crate::cryptopals::Report;
use crate::symmetric::aes::{decrypt_block, encrypt_block, AesKey};

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

fn ecb_dec(ct: &[u8], key: &AesKey) -> Vec<u8> {
    let mut out = Vec::new();
    for chunk in ct.chunks_exact(16) {
        let mut b = [0u8; 16];
        b.copy_from_slice(chunk);
        out.extend_from_slice(&decrypt_block(&b, key));
    }
    out
}

fn profile_for(email: &str) -> String {
    let sanitised: String = email
        .chars()
        .filter(|c| *c != '&' && *c != '=')
        .collect();
    format!("email={}&uid=10&role=user", sanitised)
}

pub fn run() -> Report {
    let mut r = Report::new(13, "ECB cut-and-paste");
    let key_bytes: [u8; 16] = *b"keyfor13demoXXyy";
    let key = AesKey::new(&key_bytes).unwrap();
    let oracle_encrypt = |email: &str| ecb_enc(profile_for(email).as_bytes(), &key);
    let oracle_decrypt = |ct: &[u8]| ecb_dec(ct, &key);

    // Step 1: get a ciphertext block whose plaintext is
    // `admin\x0b\x0b\x0b\x0b\x0b\x0b\x0b\x0b\x0b\x0b\x0b` (PKCS#7
    // padding to 16).  Craft an email so that's the second block.
    // "email=" is 6 bytes; pad to 16 with 10 chars, then start the
    // payload block.
    let mut crafted_email = String::from("AAAAAAAAAA"); // 10 chars → fills block 0 after "email="
    crafted_email.push_str(&"admin".to_string());
    crafted_email.push_str(&"\x0b".repeat(11));
    let ct1 = oracle_encrypt(&crafted_email);
    let admin_block = ct1[16..32].to_vec();

    // Step 2: get a ciphertext whose final block is `role=` + something,
    // and crucially whose role= aligns to a block boundary.
    //   "email=AAAA@AA.com&uid=10&role=user"
    //   Need "role=" to start at byte index = N·16.  "email=" is 6,
    //   then email, then "&uid=10&role=" is 13.  6 + email_len + 13
    //   must be a multiple of 16 → email_len ∈ {13, 29, …}.
    let target_email = "AAAAAA@CD.COM"; // 13 bytes
    assert_eq!(target_email.len(), 13);
    let ct2 = oracle_encrypt(target_email);
    // Keep all but the last block, splice in admin_block.
    let mut forged = ct2[..ct2.len() - 16].to_vec();
    forged.extend_from_slice(&admin_block);
    let decoded = oracle_decrypt(&forged);
    let stripped = crate::cryptopals::low_util::pkcs7_unpad(&decoded, 16).unwrap();
    let role = std::str::from_utf8(&stripped).unwrap();
    r.line(format!("Forged plaintext: {:?}", role));
    assert!(role.ends_with("&role=admin"));
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn forges_admin() {
        assert!(super::run().success);
    }
}
