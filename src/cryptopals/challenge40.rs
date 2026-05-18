//! # Challenge 40 — Broadcast RSA e=3 attack (Håstad)
//!
//! Encrypt the same message under three different `e=3` RSA keys.
//! Use CRT to recover `m^3 (mod n0·n1·n2)`, then take the integer
//! cube root — that's `m`.

use crate::cryptopals::Report;
use crate::cryptopals::set8_util::crt_combine;
use num_bigint::BigUint;
use num_traits::One;

fn int_cube_root(n: &BigUint) -> BigUint {
    // Newton's method.
    let two = BigUint::from(2u32);
    let three = BigUint::from(3u32);
    let mut x = BigUint::one() << ((n.bits() / 3) as usize + 1);
    loop {
        let next = (&x * &two + n / (&x * &x)) / &three;
        if next >= x {
            return x;
        }
        x = next;
    }
}

pub fn run() -> Report {
    let mut r = Report::new(40, "Hastad broadcast e=3");
    let mut keys = Vec::new();
    for _ in 0..3 {
        keys.push(crate::asymmetric::rsa::RsaKeyPair::generate(384));
    }
    // Make sure all e = 3; current RSA module uses 65537 by default,
    // but the attack just needs c_i = m^e mod n_i for shared small e.
    // If e != 3, the demo still illustrates CRT but needs the same e.
    // Force e = 3 by replacing the public exponent in a fresh key.
    let msg = b"the cube root of unity";
    let m = BigUint::from_bytes_be(msg);

    // For demo purposes, manually construct three (N, 3) pairs.
    // We need each n_i coprime to 3 and m³ < ∏ n_i.
    let mut pairs: Vec<(BigUint, BigUint)> = Vec::new();
    for i in 0..3 {
        let kp = &keys[i];
        let n = kp.public.n.clone();
        let c = m.modpow(&BigUint::from(3u32), &n);
        pairs.push((c, n));
    }
    let (cube, _) = crt_combine(&pairs);
    let cube_root = int_cube_root(&cube);
    let recovered = cube_root.to_bytes_be();
    r.line(format!("Recovered: {:?}", std::str::from_utf8(&recovered).unwrap_or("?")));
    assert!(recovered == msg.to_vec() || recovered == msg.to_vec().split_off(0));
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_via_crt() {
        assert!(super::run().success);
    }
}
