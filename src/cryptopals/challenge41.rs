//! # Challenge 41 — Unpadded RSA message-recovery oracle
//!
//! Server decrypts any RSA ciphertext except a specific one it
//! recorded.  Attacker submits `C' = (S^e · C) mod N` for any
//! `S > 1`.  Server returns `P' = S · P mod N`; attacker divides by
//! `S` to get `P`.

use crate::cryptopals::Report;
use num_bigint::BigUint;
use num_integer::Integer;

fn mod_inv(a: &BigUint, n: &BigUint) -> BigUint {
    use num_bigint::ToBigInt;
    let a_i = a.to_bigint().unwrap();
    let n_i = n.to_bigint().unwrap();
    let g = a_i.extended_gcd(&n_i);
    let r = ((g.x % &n_i) + &n_i) % &n_i;
    r.to_biguint().unwrap()
}

pub fn run() -> Report {
    let mut r = Report::new(41, "Unpadded RSA message-recovery oracle");
    let kp = crate::asymmetric::rsa::RsaKeyPair::generate(384);
    let n = kp.public.n.clone();
    let e = kp.public.e.clone();
    let secret = b"sensitive message";
    let m = BigUint::from_bytes_be(secret);
    let c = m.modpow(&e, &n);
    // Attacker picks S = 2.
    let s = BigUint::from(2u32);
    let c_prime = (&s.modpow(&e, &n) * &c) % &n;
    // Raw decrypt — no PKCS#1 unpadding.  Cryptopals assumes a
    // textbook RSA oracle here.
    let p_prime = c_prime.modpow(&kp.private.d, &n);
    // Recover m = p_prime / s mod n
    let m_recovered = (&p_prime * mod_inv(&s, &n)) % &n;
    let bytes = m_recovered.to_bytes_be();
    r.line(format!("Recovered: {:?}", std::str::from_utf8(&bytes).unwrap_or("?")));
    assert_eq!(bytes, secret);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_unpadded() {
        assert!(super::run().success);
    }
}
