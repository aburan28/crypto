//! # Challenge 37 — Break SRP with A = 0 / N / kN
//!
//! Client sends `A = 0`.  Server computes
//! `S = (A · v^u)^b = (0 · …)^b = 0`.  Client knows `S = 0` too,
//! without ever knowing the password.

use crate::cryptopals::challenge33::nist_p_g;
use crate::cryptopals::set8_util::{biguint_to_bytes_be, hmac_sha256};
use crate::cryptopals::Report;
use crate::hash::sha256::sha256;
use num_bigint::{BigUint, RandBigInt};
use num_traits::Zero;
use rand::SeedableRng;
use rand::rngs::StdRng;

pub fn run() -> Report {
    let mut r = Report::new(37, "Break SRP with A = 0 / N");
    let (n, g) = nist_p_g();
    let k = BigUint::from(3u32);
    let password = b"any password the attacker pleases";
    let mut rng = StdRng::seed_from_u64(37);
    let salt = rng.gen_biguint(64);
    let salt_bytes = biguint_to_bytes_be(&salt, 8);
    let mut x_in = salt_bytes.clone();
    x_in.extend_from_slice(password);
    let x = BigUint::from_bytes_be(&sha256(&x_in));
    let v = g.modpow(&x, &n);

    for evil in [BigUint::zero(), n.clone(), &n * BigUint::from(2u32)] {
        // Server still computes B normally.
        let b_priv = rng.gen_biguint_range(&BigUint::from(2u32), &n);
        let cap_b = (&k * &v + g.modpow(&b_priv, &n)) % &n;
        let cap_a = evil.clone();
        let _u = {
            let mut buf = biguint_to_bytes_be(&cap_a, 256);
            buf.extend(biguint_to_bytes_be(&cap_b, 256));
            BigUint::from_bytes_be(&sha256(&buf))
        };
        // Server S = (A v^u)^b = (0)^b = 0
        let s_server = (&cap_a * v.modpow(&_u, &n) % &n).modpow(&b_priv, &n);
        assert!(s_server.is_zero());
        // Attacker uses S = 0 to derive K.
        let key_attacker = sha256(&biguint_to_bytes_be(&BigUint::zero(), 256));
        let key_server = sha256(&biguint_to_bytes_be(&s_server, 256));
        let mac_a = hmac_sha256(&key_attacker, &salt_bytes);
        let mac_s = hmac_sha256(&key_server, &salt_bytes);
        assert_eq!(mac_a, mac_s);
    }
    r.line("All evil A values forged auth without password.");
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn srp_zero_a_works() {
        assert!(super::run().success);
    }
}
