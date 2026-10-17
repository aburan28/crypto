//! Enumerate the X9.62 candidates of one 1997–98 NIST seed-generation run,
//! for the point-count census in `research/nist_seed_provenance_20261016/`.
//!
//! ```bash
//! cargo run --release --example nist_seed_candidates -- --run 871696354 --blocks 1842
//! ```
//!
//! For each 20-byte block `i = 1..=blocks` of `(int)(drand48() * 256)` after
//! `srand48(run)`, it derives the X9.62 A.3.3.1 value `r` for the run's field
//! (`p` from [`PRIME_SEEDS`]), drops the degenerate cases `r ≡ 0` and
//! `4r + 27 ≡ 0`, and prints the blocks where `b² = −27/r` is a square, as
//! `i seed_hex r b2 b0` with `b0 = b2^((p+1)/4) mod p` when `p ≡ 3 (mod 4)`
//! and `0` otherwise (the census then takes the root in GP).  The output is
//! the input of `sea_worker.gp`.

use std::process::exit;

use crypto_lib::cryptanalysis::nist_seed_provenance::{Drand48, PRIME_SEEDS};
use crypto_lib::cryptanalysis::p256_backdoor_map::x962_candidate_r;
use num_bigint::BigUint;

fn arg(name: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1).cloned())
}

fn main() {
    let (Some(run), Some(blocks)) = (
        arg("--run").and_then(|s| s.parse::<u32>().ok()),
        arg("--blocks").and_then(|s| s.parse::<u32>().ok()),
    ) else {
        eprintln!("usage: nist_seed_candidates --run <srand48 value> --blocks <n>");
        exit(2);
    };
    let Some(curve) = PRIME_SEEDS.iter().find(|s| s.srand48 == run) else {
        eprintln!("srand48({run}) is not one of the published runs");
        exit(2);
    };
    let p = BigUint::parse_bytes(curve.p_hex.as_bytes(), 16).expect("p hex");
    let t_bits = p.bits() as u32;
    let zero = BigUint::from(0u8);
    let one = BigUint::from(1u8);
    let m27 = &p - BigUint::from(27u8);
    let half = (&p - &one) >> 1;
    let p_mod4 = (&p % BigUint::from(4u8)) == BigUint::from(3u8);
    let quarter = (&p + &one) >> 2;

    let mut g = Drand48::srand48(run);
    for i in 1..=blocks {
        let mut seed = [0u8; 20];
        for b in seed.iter_mut() {
            *b = g.byte();
        }
        let (_, _, r) = x962_candidate_r(&seed, t_bits);
        let rp = &r % &p;
        if rp == zero || (BigUint::from(4u8) * &rp + BigUint::from(27u8)) % &p == zero {
            continue;
        }
        let r_inv = rp.modpow(&(&p - BigUint::from(2u8)), &p);
        let b2 = &m27 * &r_inv % &p;
        if b2.modpow(&half, &p) != one {
            continue;
        }
        let b0 = if p_mod4 {
            b2.modpow(&quarter, &p)
        } else {
            zero.clone()
        };
        let hex: String = seed.iter().map(|b| format!("{b:02x}")).collect();
        println!("{i} {hex} {r} {b2} {b0}");
    }
}
