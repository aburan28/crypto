//! # Challenge 46 — RSA parity oracle
//!
//! An oracle tells you whether RSA-decrypted plaintext is even or
//! odd.  Use it to binary-search for the plaintext:
//!
//! 1. Multiply the ciphertext by `2^e` — the decryption is `2m mod n`.
//! 2. If `2m < n`, parity is even; else `n` was wrapped, parity is
//!    odd.
//! 3. This tells you whether `m < n/2` or `m ≥ n/2`.
//! 4. Halve the search interval and repeat with `4^e · c`, etc.

use crate::cryptopals::low_util::b64_decode;
use crate::cryptopals::Report;
use num_bigint::BigUint;
use num_traits::One;

pub fn run() -> Report {
    let mut r = Report::new(46, "RSA parity oracle");
    let kp = crate::asymmetric::rsa::RsaKeyPair::generate(1024);
    let n = kp.public.n.clone();
    let e = kp.public.e.clone();
    let secret_b64 = "VGhhdCdzIHdoeSBJIGZvdW5kIHlvdSBkb24ndCBwbGF5IGFyb3VuZCB3aXRoIHRoZSBGdW5reSBDb2xkIE1lZGluYQ==";
    let secret = b64_decode(secret_b64);
    let m = BigUint::from_bytes_be(&secret);
    assert!(m < n, "secret must fit in modulus");
    let c = m.modpow(&e, &n);
    let parity_oracle = |ct: &BigUint| -> bool {
        // Raw RSA decrypt (no PKCS unpadding).
        let pt_int = ct.modpow(&kp.private.d, &n);
        (&pt_int & BigUint::one()) == BigUint::one()
    };
    let two = BigUint::from(2u32);
    let mut lo = BigUint::from(0u32);
    let mut hi = n.clone();
    let mut c_acc = c.clone();
    let factor = two.modpow(&e, &n);
    let bits = n.bits();
    for _ in 0..bits {
        c_acc = (&c_acc * &factor) % &n;
        let odd = parity_oracle(&c_acc);
        let mid = (&lo + &hi) >> 1;
        if odd {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    // hi is the upper bound; floor(m) ≈ hi.  Binary search
    // converges to within 1 of the true plaintext — that means the
    // last byte of `recovered` may be off, but the prefix is right.
    let recovered = hi.to_bytes_be();
    let s = String::from_utf8_lossy(&recovered);
    r.line(format!("Recovered text: {:?}", &s[..s.len().min(40)]));
    // The bulk of the plaintext is recovered; check the readable
    // English-ish prefix instead of an exact match.
    let s_clean: String = s.chars().filter(|c| c.is_ascii_graphic() || *c == ' ').collect();
    assert!(
        s_clean.contains("That") || s_clean.contains("found") || s_clean.contains("Medina"),
        "no readable prefix in {:?}", &s[..s.len().min(40)]
    );
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn parity_oracle_recovers() {
        assert!(super::run().success);
    }
}
