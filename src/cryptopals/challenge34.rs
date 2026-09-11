//! # Challenge 34 — MITM key-fixing attack on DH
//!
//! Alice → Bob: (p, g, A). Bob → Alice: B.  Both derive `K = pub^priv`.
//!
//! MITM Mallory replaces both `A` and `B` with `p`.  Then both
//! Alice and Bob compute `K = p^priv mod p = 0`, and Mallory knows
//! the key in advance: `K = SHA-256(0)`.

use crate::cryptopals::challenge33::nist_p_g;
use crate::cryptopals::challenge10::cbc_decrypt_no_iv_prefix;
use crate::cryptopals::Report;
use crate::hash::sha256::sha256;
use crate::symmetric::aes::AesKey;
use num_bigint::{BigUint, RandBigInt};
use num_traits::Zero;
use rand::SeedableRng;
use rand::rngs::StdRng;

pub fn run() -> Report {
    let mut r = Report::new(34, "MITM key-fixing on DH");
    let (p, g) = nist_p_g();
    let mut rng = StdRng::seed_from_u64(34);
    let a = rng.gen_biguint_range(&BigUint::from(2u32), &p);
    let b = rng.gen_biguint_range(&BigUint::from(2u32), &p);
    let _pub_a = g.modpow(&a, &p);
    let _pub_b = g.modpow(&b, &p);
    // Mallory replaces both with p.  Alice computes p^a mod p = 0.
    // Bob computes p^b mod p = 0.
    let alice_shared = p.modpow(&a, &p);
    let bob_shared = p.modpow(&b, &p);
    assert!(alice_shared.is_zero() && bob_shared.is_zero());
    // Both derive same AES key from SHA-256(0).
    let key_alice = AesKey::new(&sha256(b"\x00")[..16]).unwrap();
    let key_bob = AesKey::new(&sha256(b"\x00")[..16]).unwrap();
    let _ = (key_alice, key_bob);
    // Alice sends Bob a message; Mallory intercepts and decrypts.
    let msg = b"the eagle has landed";
    let iv = [0u8; 16];
    let ct = crate::cryptopals::challenge10::cbc_encrypt_no_iv_prefix(
        msg,
        &AesKey::new(&sha256(b"\x00")[..16]).unwrap(),
        &iv,
    );
    let mallory_key = AesKey::new(&sha256(b"\x00")[..16]).unwrap();
    let pt = cbc_decrypt_no_iv_prefix(&ct, &mallory_key, &iv).unwrap();
    r.line(format!("Mallory decrypts: {:?}", std::str::from_utf8(&pt).unwrap()));
    assert_eq!(pt, msg);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn mitm_recovers_plaintext() {
        assert!(super::run().success);
    }
}
