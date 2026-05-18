//! # Challenge 38 — Offline dictionary attack on simplified SRP
//!
//! Simplified SRP drops the `k·v` term.  Now `B = g^b` and
//! `S = B^(a + u·x)`.  Malicious "server" picks any `b, u, salt`
//! and learns `MAC(S)` from the client.  Then offline:
//! `S_guess = (A · g^(u·x_guess))^b` for every candidate password.
//! Pick the candidate whose `MAC(SHA256(S_guess))` matches.

use crate::cryptopals::challenge33::nist_p_g;
use crate::cryptopals::set8_util::{biguint_to_bytes_be, hmac_sha256};
use crate::cryptopals::Report;
use crate::hash::sha256::sha256;
use num_bigint::{BigUint, RandBigInt};
use rand::SeedableRng;
use rand::rngs::StdRng;

pub fn run() -> Report {
    let mut r = Report::new(38, "Offline dict attack on simplified SRP");
    let (n, g) = nist_p_g();
    let password = b"hunter2"; // weak
    let dictionary: Vec<&[u8]> = vec![b"password", b"qwerty", b"hunter2", b"abc123"];
    let mut rng = StdRng::seed_from_u64(38);
    let salt = rng.gen_biguint(64);
    let salt_bytes = biguint_to_bytes_be(&salt, 8);
    let a = rng.gen_biguint_range(&BigUint::from(2u32), &n);
    let cap_a = g.modpow(&a, &n);
    // Malicious server picks b, u, salt freely.
    let b = rng.gen_biguint_range(&BigUint::from(2u32), &n);
    let cap_b = g.modpow(&b, &n);
    let u = BigUint::from(1u32);
    // Client computes S, K, MAC.
    let mut x_in = salt_bytes.clone();
    x_in.extend_from_slice(password);
    let x = BigUint::from_bytes_be(&sha256(&x_in));
    let s_client = cap_b.modpow(&(&a + &u * &x), &n);
    let k_client = sha256(&biguint_to_bytes_be(&s_client, 256));
    let mac_client = hmac_sha256(&k_client, &salt_bytes);
    // Attacker tries every dict entry.
    let mut found = None;
    for cand in &dictionary {
        let mut x_in = salt_bytes.clone();
        x_in.extend_from_slice(cand);
        let x_cand = BigUint::from_bytes_be(&sha256(&x_in));
        let v = g.modpow(&x_cand, &n);
        // S_attack = (A · v^u)^b
        let s = (&cap_a * v.modpow(&u, &n) % &n).modpow(&b, &n);
        let k = sha256(&biguint_to_bytes_be(&s, 256));
        let mac = hmac_sha256(&k, &salt_bytes);
        if mac == mac_client {
            found = Some(*cand);
            break;
        }
    }
    let p = found.expect("dictionary contains the password");
    r.line(format!("Recovered password: {:?}", std::str::from_utf8(p).unwrap()));
    assert_eq!(p, password);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_dict_password() {
        assert!(super::run().success);
    }
}
