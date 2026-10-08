//! Time `Φ_ℓ mod p` construction over the P-256 field for a few ℓ.
use crypto_lib::cryptanalysis::isogeny_walk::{field::Field, modpoly};
use num_bigint::BigUint;
use std::time::Instant;

fn main() {
    let p = BigUint::parse_bytes(
        b"115792089210356248762697446949407573530086143415290314195533631308867097853951",
        10,
    )
    .unwrap();
    let f = Field::new(&p).unwrap();
    let ells: Vec<u64> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    for l in if ells.is_empty() { vec![11, 17, 23, 31, 41, 47, 61] } else { ells } {
        let t = Instant::now();
        let phi = modpoly::modular_polynomial(&f, l).unwrap();
        println!("ell={l} ms={} checks={}", t.elapsed().as_millis(), phi.checks_passed);
    }
}
