//! # Challenge 33 — Implement Diffie-Hellman
//!
//! Vanilla DH key agreement with the RFC 3526 2048-bit MODP group.

use crate::cryptopals::set8_util::parse_big;
use crate::cryptopals::Report;
use num_bigint::{BigUint, RandBigInt};
use rand::SeedableRng;
use rand::rngs::StdRng;

pub fn nist_p_g() -> (BigUint, BigUint) {
    let p = BigUint::parse_bytes(
        b"ffffffffffffffffc90fdaa22168c234c4c6628b80dc1cd129024e088a67cc74020bbea63b139b22514a08798e3404ddef9519b3cd3a431b302b0a6df25f14374fe1356d6d51c245e485b576625e7ec6f44c42e9a637ed6b0bff5cb6f406b7edee386bfb5a899fa5ae9f24117c4b1fe649286651ece45b3dc2007cb8a163bf0598da48361c55d39a69163fa8fd24cf5f83655d23dca3ad961c62f356208552bb9ed529077096966d670c354e4abc9804f1746c08ca18217c32905e462e36ce3be39e772c180e86039b2783a2ec07a28fb5c55df06f4c52c9de2bcbf6955817183995497cea956ae515d2261898fa051015728e5a8aacaa68ffffffffffffffff",
        16,
    ).unwrap();
    let g = BigUint::from(2u32);
    (p, g)
}

pub fn run() -> Report {
    let mut r = Report::new(33, "Implement Diffie-Hellman");
    let (p, g) = nist_p_g();
    let mut rng = StdRng::seed_from_u64(33);
    let a = rng.gen_biguint_range(&BigUint::from(2u32), &p);
    let b = rng.gen_biguint_range(&BigUint::from(2u32), &p);
    let pub_a = g.modpow(&a, &p);
    let pub_b = g.modpow(&b, &p);
    let s_a = pub_b.modpow(&a, &p);
    let s_b = pub_a.modpow(&b, &p);
    r.line(format!("p has {} bits", p.bits()));
    r.line(format!("shared secret matches: {}", s_a == s_b));
    assert_eq!(s_a, s_b);
    let _ = parse_big;
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn dh_round_trip() {
        assert!(super::run().success);
    }
}
