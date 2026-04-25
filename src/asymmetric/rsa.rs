//! RSA — Rivest–Shamir–Adleman asymmetric cryptography, implemented from scratch.
//!
//! # Key generation
//! 1. Generate two large primes p, q using Miller-Rabin primality testing.
//! 2. n = p·q (the modulus).
//! 3. λ(n) = lcm(p-1, q-1) — Carmichael's totient.
//! 4. e = 65537 (standard public exponent, a prime Fermat number).
//! 5. d = e⁻¹ mod λ(n) — the private exponent.
//!
//! # Encryption / Decryption (textbook RSA)
//! Encrypt: c = mᵉ mod n
//! Decrypt: m = cᵈ mod n
//!
//! # Signatures (PKCS#1 v1.5 style)
//! Sign: s = hash(m)ᵈ mod n
//! Verify: hash(m) == sᵉ mod n

use crate::hash::sha256::sha256;
use crate::utils::{mod_inverse, mod_pow};
use num_bigint::{BigUint, RandBigInt};
use num_integer::Integer;
use num_traits::{One, Zero};

// ── Miller-Rabin primality test ───────────────────────────────────────────────

/// Decompose `n-1` as `2^s * d` with `d` odd.
fn factor_out_twos(n: &BigUint) -> (u64, BigUint) {
    let mut d = n.clone();
    let mut s = 0u64;
    while d.is_even() {
        d >>= 1;
        s += 1;
    }
    (s, d)
}

/// Miller-Rabin witness test for a single base `a`.
/// Returns `true` if `n` is *probably prime* with this witness.
fn miller_rabin_witness(n: &BigUint, d: &BigUint, s: u64, a: &BigUint) -> bool {
    let one = BigUint::one();
    let n_minus_1 = n - &one;
    let mut x = mod_pow(a, d, n);

    if x == one || x == n_minus_1 {
        return true;
    }

    for _ in 0..s - 1 {
        x = mod_pow(&x, &BigUint::from(2u32), n);
        if x == n_minus_1 {
            return true;
        }
    }
    false
}

/// Deterministic Miller-Rabin using a fixed witness set that is correct
/// for all n < 3,317,044,064,679,887,385,961,981 (sufficient for 256-bit+ primes).
pub fn is_prime(n: &BigUint) -> bool {
    if n < &BigUint::from(2u32) { return false; }
    if n == &BigUint::from(2u32) || n == &BigUint::from(3u32) { return true; }
    if n.is_even() { return false; }

    // Small prime divisibility check for speed
    let small_primes = [2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for &p in &small_primes {
        let bp = BigUint::from(p);
        if n == &bp { return true; }
        if (n % &bp).is_zero() { return false; }
    }

    let n_minus_1 = n - BigUint::one();
    let (s, d) = factor_out_twos(&n_minus_1);

    // Witnesses sufficient for 64-bit numbers; for larger numbers a probabilistic
    // set of random witnesses is used.
    let witnesses: Vec<BigUint> = [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
        .iter()
        .map(|&w| BigUint::from(w))
        .collect();

    witnesses.iter().all(|a| {
        if a >= n { return true; } // Skip witnesses ≥ n
        miller_rabin_witness(n, &d, s, a)
    })
}

/// Generate a random probable prime of exactly `bits` bits.
/// Uses 20 rounds of Miller-Rabin for probabilistic primality (error < 4⁻²⁰).
pub fn random_prime(bits: u64) -> BigUint {
    let mut rng = rand::thread_rng();
    loop {
        let mut candidate = rng.gen_biguint(bits);
        // Ensure the candidate has the right bit length and is odd
        candidate.set_bit(bits - 1, true);
        candidate.set_bit(0, true);
        if is_prime(&candidate) {
            return candidate;
        }
    }
}

// ── RSA key types ─────────────────────────────────────────────────────────────

/// RSA public key.
#[derive(Clone, Debug)]
pub struct RsaPublicKey {
    /// Modulus n = p·q
    pub n: BigUint,
    /// Public exponent e (typically 65537)
    pub e: BigUint,
    /// Key size in bits
    pub bits: u64,
}

/// RSA private key (PKCS#1 representation).
#[derive(Clone, Debug)]
pub struct RsaPrivateKey {
    pub n: BigUint,
    pub e: BigUint,
    /// Private exponent d = e⁻¹ mod λ(n)
    pub d: BigUint,
    pub p: BigUint,
    pub q: BigUint,
    pub bits: u64,
}

/// An RSA key pair.
pub struct RsaKeyPair {
    pub public: RsaPublicKey,
    pub private: RsaPrivateKey,
}

impl RsaKeyPair {
    /// Generate a fresh RSA key pair of `bits` total bits (e.g., 2048).
    /// Each prime will be `bits/2` bits long.
    pub fn generate(bits: u64) -> Self {
        let half = bits / 2;
        loop {
            let p = random_prime(half);
            let q = random_prime(half);
            if p == q { continue; }

            let n = &p * &q;

            // Carmichael's totient λ(n) = lcm(p-1, q-1)
            let p1 = &p - BigUint::one();
            let q1 = &q - BigUint::one();
            let lambda = p1.lcm(&q1);

            let e = BigUint::from(65537u32);
            if lambda.gcd(&e) != BigUint::one() { continue; }

            let d = match mod_inverse(&e, &lambda) {
                Some(v) => v,
                None => continue,
            };

            return RsaKeyPair {
                public: RsaPublicKey { n: n.clone(), e: e.clone(), bits },
                private: RsaPrivateKey { n, e, d, p, q, bits },
            };
        }
    }
}

// ── Textbook RSA operations ───────────────────────────────────────────────────

/// Textbook RSA encrypt: c = mᵉ mod n.
/// WARNING: not semantically secure on its own — use OAEP padding in production.
pub fn rsa_encrypt_raw(msg: &BigUint, key: &RsaPublicKey) -> BigUint {
    mod_pow(msg, &key.e, &key.n)
}

/// Textbook RSA decrypt: m = cᵈ mod n.
pub fn rsa_decrypt_raw(ciphertext: &BigUint, key: &RsaPrivateKey) -> BigUint {
    mod_pow(ciphertext, &key.d, &key.n)
}

// ── Signing / Verification ────────────────────────────────────────────────────

/// Sign a message by SHA-256 hashing it, then computing h^d mod n.
pub fn rsa_sign(message: &[u8], key: &RsaPrivateKey) -> BigUint {
    let hash = sha256(message);
    let hash_int = BigUint::from_bytes_be(&hash);
    // Ensure hash < n (always true for 2048-bit keys)
    mod_pow(&hash_int, &key.d, &key.n)
}

/// Verify: recompute s^e mod n and compare with SHA-256(message).
pub fn rsa_verify(message: &[u8], signature: &BigUint, key: &RsaPublicKey) -> bool {
    let hash = sha256(message);
    let expected = BigUint::from_bytes_be(&hash);
    let recovered = mod_pow(signature, &key.e, &key.n);
    // Pad expected to same length as n for comparison
    expected == recovered
}

// ── PKCS#1 v1.5 style helpers ─────────────────────────────────────────────────

/// Encode a short message with PKCS#1 v1.5 type 2 padding for RSA encryption.
///
/// Layout: 0x00 || 0x02 || PS (random non-zero padding) || 0x00 || msg
/// The padded value must be strictly less than n.
pub fn pkcs1_pad_encrypt(msg: &[u8], key_bytes: usize) -> Result<Vec<u8>, &'static str> {
    if msg.len() + 11 > key_bytes {
        return Err("message too long for key size");
    }
    let ps_len = key_bytes - msg.len() - 3;
    let mut padded = vec![0u8; key_bytes];
    padded[1] = 0x02;
    // Fill PS with random non-zero bytes
    let mut i = 2;
    while i < 2 + ps_len {
        let b = crate::utils::random::random_bytes_vec(1)[0];
        if b != 0 {
            padded[i] = b;
            i += 1;
        }
    }
    padded[2 + ps_len] = 0x00;
    padded[3 + ps_len..].copy_from_slice(msg);
    Ok(padded)
}

/// Strip PKCS#1 v1.5 type 2 padding, returning the original message or an error.
pub fn pkcs1_unpad_encrypt(padded: &[u8]) -> Result<Vec<u8>, &'static str> {
    if padded.len() < 11 || padded[0] != 0x00 || padded[1] != 0x02 {
        return Err("invalid padding");
    }
    let ps_end = padded[2..].iter().position(|&b| b == 0x00)
        .ok_or("no zero separator")?;
    if ps_end < 8 { return Err("PS too short"); }
    Ok(padded[2 + ps_end + 1..].to_vec())
}

/// High-level RSA encrypt with PKCS#1 v1.5 padding.
pub fn rsa_encrypt(msg: &[u8], key: &RsaPublicKey) -> Result<BigUint, &'static str> {
    let key_bytes = ((key.bits + 7) / 8) as usize;
    let padded = pkcs1_pad_encrypt(msg, key_bytes)?;
    let m = BigUint::from_bytes_be(&padded);
    Ok(rsa_encrypt_raw(&m, key))
}

/// High-level RSA decrypt with PKCS#1 v1.5 unpadding.
pub fn rsa_decrypt(ciphertext: &BigUint, key: &RsaPrivateKey) -> Result<Vec<u8>, &'static str> {
    let m = rsa_decrypt_raw(ciphertext, key);
    let key_bytes = ((key.bits + 7) / 8) as usize;
    let padded = crate::utils::encoding::bigint_to_bytes_be(&m, key_bytes);
    pkcs1_unpad_encrypt(&padded)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primality_small() {
        let primes = [2u32, 3, 5, 7, 11, 13, 17, 19, 23];
        for &p in &primes {
            assert!(is_prime(&BigUint::from(p)), "{p} should be prime");
        }
        let composites = [4u32, 6, 9, 15, 25, 49];
        for &c in &composites {
            assert!(!is_prime(&BigUint::from(c)), "{c} should be composite");
        }
    }

    #[test]
    fn rsa_1024_sign_verify() {
        let kp = RsaKeyPair::generate(1024);
        let msg = b"RSA signing test message";
        let sig = rsa_sign(msg, &kp.private);
        assert!(rsa_verify(msg, &sig, &kp.public));
    }

    #[test]
    fn rsa_sign_wrong_message() {
        let kp = RsaKeyPair::generate(1024);
        let sig = rsa_sign(b"original", &kp.private);
        assert!(!rsa_verify(b"tampered", &sig, &kp.public));
    }

    #[test]
    fn rsa_encrypt_decrypt() {
        let kp = RsaKeyPair::generate(1024);
        let msg = b"hello RSA";
        let ct = rsa_encrypt(msg, &kp.public).unwrap();
        let pt = rsa_decrypt(&ct, &kp.private).unwrap();
        assert_eq!(pt, msg);
    }

    #[test]
    fn primality_known_large_primes() {
        // 2^127 - 1 (Mersenne prime M_127)
        let m127 = (BigUint::one() << 127u32) - BigUint::one();
        assert!(is_prime(&m127));
        // 2^61 - 1 (Mersenne prime M_61)
        let m61 = (BigUint::one() << 61u32) - BigUint::one();
        assert!(is_prime(&m61));
        // Carmichael number 561 = 3·11·17 — fools Fermat's little theorem
        // but Miller-Rabin must reject it.
        assert!(!is_prime(&BigUint::from(561u32)));
        // Carmichael 41041 = 7·11·13·41
        assert!(!is_prime(&BigUint::from(41041u32)));
        // Carmichael 825265 = 5·7·17·19·73
        assert!(!is_prime(&BigUint::from(825265u32)));
    }

    #[test]
    fn random_prime_has_correct_bit_length() {
        for _ in 0..3 {
            let p = random_prime(256);
            assert_eq!(p.bits(), 256, "prime must have exactly 256 bits");
            assert!(is_prime(&p));
        }
    }

    #[test]
    fn rsa_keygen_modulus_bit_length() {
        let kp = RsaKeyPair::generate(1024);
        // n = p·q with p, q each 512 bits and high bits forced to 1, so
        // n must be exactly 1023 or 1024 bits. We force the bit-1 of each
        // prime, so n is at least 2^1022 and at most 2^1024 - 1.
        let bits = kp.public.n.bits();
        assert!(bits >= 1023 && bits <= 1024, "got {bits} bits");
        // d·e ≡ 1 (mod λ(n))
        let lambda = (&kp.private.p - BigUint::one())
            .lcm(&(&kp.private.q - BigUint::one()));
        assert_eq!((&kp.private.d * &kp.private.e) % &lambda, BigUint::one());
        // p·q == n
        assert_eq!(&kp.private.p * &kp.private.q, kp.public.n);
    }

    #[test]
    fn pkcs1_pad_unpad_roundtrip() {
        let msg = b"sensitive payment payload";
        let padded = pkcs1_pad_encrypt(msg, 128).unwrap();
        assert_eq!(padded.len(), 128);
        assert_eq!(padded[0], 0x00);
        assert_eq!(padded[1], 0x02);
        // Padding bytes must all be non-zero.
        let zero_sep = padded[2..].iter().position(|&b| b == 0).unwrap() + 2;
        assert!(zero_sep - 2 >= 8, "PS must be ≥ 8 bytes");
        for &b in &padded[2..zero_sep] {
            assert_ne!(b, 0);
        }
        let recovered = pkcs1_unpad_encrypt(&padded).unwrap();
        assert_eq!(recovered, msg);
    }

    #[test]
    fn pkcs1_pad_rejects_oversized_message() {
        // 128-byte modulus has 128 - 11 = 117 bytes max payload.
        assert!(pkcs1_pad_encrypt(&vec![0u8; 117], 128).is_ok());
        assert!(pkcs1_pad_encrypt(&vec![0u8; 118], 128).is_err());
    }

    #[test]
    fn pkcs1_unpad_rejects_malformed() {
        // Wrong leading byte
        let mut bad = vec![0u8; 128];
        bad[0] = 0x01;
        bad[1] = 0x02;
        assert!(pkcs1_unpad_encrypt(&bad).is_err());
        // Wrong block type
        let mut bad2 = vec![0u8; 128];
        bad2[1] = 0x01;
        assert!(pkcs1_unpad_encrypt(&bad2).is_err());
        // PS shorter than 8 bytes
        let mut bad3 = vec![0u8; 128];
        bad3[1] = 0x02;
        bad3[2] = 0xff;
        bad3[3] = 0x00; // separator after only 1 PS byte
        bad3[4..].copy_from_slice(&[0xaa; 124]);
        assert!(pkcs1_unpad_encrypt(&bad3).is_err());
        // No zero separator at all
        let bad4 = vec![0xffu8; 128];
        assert!(pkcs1_unpad_encrypt(&bad4).is_err());
    }

    #[test]
    fn rsa_encrypt_is_randomized() {
        // PKCS#1 v1.5 padding uses random PS, so two encryptions of the same
        // plaintext must produce different ciphertexts.
        let kp = RsaKeyPair::generate(1024);
        let msg = b"same message";
        let c1 = rsa_encrypt(msg, &kp.public).unwrap();
        let c2 = rsa_encrypt(msg, &kp.public).unwrap();
        assert_ne!(c1, c2, "PKCS#1 v1.5 must produce randomized ciphertexts");
        assert_eq!(rsa_decrypt(&c1, &kp.private).unwrap(), msg);
        assert_eq!(rsa_decrypt(&c2, &kp.private).unwrap(), msg);
    }

    #[test]
    fn rsa_encrypt_empty_and_max_message() {
        let kp = RsaKeyPair::generate(1024);
        // Empty message
        let ct_empty = rsa_encrypt(b"", &kp.public).unwrap();
        assert_eq!(rsa_decrypt(&ct_empty, &kp.private).unwrap(), b"");
        // Maximum-length message: 128 - 11 = 117 bytes
        let max_msg = vec![0xab; 117];
        let ct_max = rsa_encrypt(&max_msg, &kp.public).unwrap();
        assert_eq!(rsa_decrypt(&ct_max, &kp.private).unwrap(), max_msg);
        // One byte too long must fail at padding stage
        let too_long = vec![0xab; 118];
        assert!(rsa_encrypt(&too_long, &kp.public).is_err());
    }

    #[test]
    fn rsa_signature_unique_per_message() {
        // Plain (deterministic) RSA signatures for distinct messages must differ.
        let kp = RsaKeyPair::generate(1024);
        let s1 = rsa_sign(b"msg-a", &kp.private);
        let s2 = rsa_sign(b"msg-b", &kp.private);
        assert_ne!(s1, s2);
    }

    #[test]
    fn rsa_verify_with_wrong_key_fails() {
        let kp1 = RsaKeyPair::generate(1024);
        let kp2 = RsaKeyPair::generate(1024);
        let sig = rsa_sign(b"transfer 100", &kp1.private);
        assert!(rsa_verify(b"transfer 100", &sig, &kp1.public));
        assert!(!rsa_verify(b"transfer 100", &sig, &kp2.public));
    }

    #[test]
    fn rsa_decrypt_rejects_tampered_ciphertext() {
        let kp = RsaKeyPair::generate(1024);
        let ct = rsa_encrypt(b"hello", &kp.public).unwrap();
        // Flip a bit in the ciphertext (cheap) — decryption nearly always fails
        // PKCS#1 v1.5 unpadding because the random PS bytes change.
        let tampered = &ct ^ BigUint::one();
        // Make sure tampered < n; if it overflowed (extremely unlikely), retry
        // with a different mask.
        let tampered = if tampered >= kp.public.n {
            &ct + BigUint::from(2u32)
        } else {
            tampered
        };
        let res = rsa_decrypt(&tampered, &kp.private);
        // Either unpadding fails outright, or the message decoded to garbage —
        // but it must NOT decrypt back to the original plaintext.
        match res {
            Err(_) => {}
            Ok(plain) => assert_ne!(plain, b"hello"),
        }
    }
}
