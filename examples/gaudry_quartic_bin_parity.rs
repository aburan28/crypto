//! How often is `H` Frobenius-reducible?  For each `m` and instance seed:
//! build `generate_instance(m, seed)`, interpolate `H` in the D₄
//! invariants (`SymmetrisedBin::precompute`, checked on 8 fresh points),
//! and report whether every monomial carrying `H` is even in `s₂, s₃, π`
//! (the parity `gaudry-d4-bin:frobenius=1` relies on), with the counts.
//!
//! `cargo run --release --example gaudry_quartic_bin_parity -- 4 18 25`
//! (m from, m to, seeds per m); one JSON line per curve.

use crypto_lib::cryptanalysis::gaudry_quartic_bin::{generate_instance, SymmetrisedBin};
use rand::rngs::StdRng;
use rand::SeedableRng;

fn main() {
    let a: Vec<u64> = std::env::args()
        .skip(1)
        .map(|s| s.parse().unwrap())
        .collect();
    let (lo, hi, n) = (a[0] as u32, a[1] as u32, a[2]);
    for m in lo..=hi {
        for seed in 1..=n {
            let curve = generate_instance(m, seed);
            let mut rng = StdRng::seed_from_u64(seed ^ 0x9e);
            let line = match SymmetrisedBin::precompute(&curve, &mut rng, 8) {
                Ok(pre) => {
                    let (carrying, odd) = pre.parity_counts();
                    format!(
                        "{{\"m\":{m},\"seed\":{seed},\"interpolated\":true,\"monomials_carrying_h\":{carrying},\"odd_in_s2_s3_pi\":{odd},\"frobenius_reducible\":{}}}",
                        pre.frobenius_reducible()
                    )
                }
                Err(e) => {
                    format!("{{\"m\":{m},\"seed\":{seed},\"interpolated\":false,\"error\":{e:?}}}")
                }
            };
            println!("{line}");
        }
    }
}
