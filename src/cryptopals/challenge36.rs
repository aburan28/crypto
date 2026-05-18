//! # Challenge 36 — Implement SRP
//!
//! SRP-6a auth, with `K = SHA-256(S)` and a final `HMAC(K, salt) ==`
//! check on the wire.

use crate::cryptopals::challenge33::nist_p_g;
use crate::cryptopals::set8_util::{biguint_to_bytes_be, hmac_sha256};
use crate::cryptopals::Report;
use crate::hash::sha256::sha256;
use num_bigint::{BigUint, RandBigInt};
use rand::SeedableRng;
use rand::rngs::StdRng;

fn h_concat(parts: &[&BigUint]) -> BigUint {
    let mut buf = Vec::new();
    for p in parts {
        buf.extend(biguint_to_bytes_be(p, 256));
    }
    BigUint::from_bytes_be(&sha256(&buf))
}

pub fn srp_handshake() -> bool {
    let (n, g) = nist_p_g();
    let k = BigUint::from(3u32);
    let i = b"alice@example.com";
    let p = b"correct horse battery staple";
    let mut rng = StdRng::seed_from_u64(36);
    // Server setup.
    let salt = rng.gen_biguint(64);
    let salt_bytes = biguint_to_bytes_be(&salt, 8);
    let mut x_in = salt_bytes.clone();
    x_in.extend_from_slice(p);
    let x = BigUint::from_bytes_be(&sha256(&x_in));
    let v = g.modpow(&x, &n);
    // Client → server: I, A = g^a
    let a = rng.gen_biguint_range(&BigUint::from(2u32), &n);
    let cap_a = g.modpow(&a, &n);
    // Server → client: salt, B = kv + g^b
    let b = rng.gen_biguint_range(&BigUint::from(2u32), &n);
    let cap_b = (&k * &v + g.modpow(&b, &n)) % &n;
    let u = h_concat(&[&cap_a, &cap_b]);
    // Client S = (B - k g^x)^(a + ux)
    let s_client = {
        let mut term = (&cap_b + &n - (&k * g.modpow(&x, &n)) % &n) % &n;
        let exp = &a + &u * &x;
        term = term.modpow(&exp, &n);
        term
    };
    // Server S = (A v^u)^b
    let s_server = {
        let term = (&cap_a * v.modpow(&u, &n)) % &n;
        term.modpow(&b, &n)
    };
    assert_eq!(s_client, s_server);
    let key_client = sha256(&biguint_to_bytes_be(&s_client, 256));
    let key_server = sha256(&biguint_to_bytes_be(&s_server, 256));
    let mac_client = hmac_sha256(&key_client, &salt_bytes);
    let mac_server = hmac_sha256(&key_server, &salt_bytes);
    let _ = i;
    mac_client == mac_server
}

pub fn run() -> Report {
    let mut r = Report::new(36, "Implement SRP");
    let ok = srp_handshake();
    r.line(format!("SRP-6a auth verifies: {}", ok));
    assert!(ok);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn srp_auths() {
        assert!(super::run().success);
    }
}
