//! Curves for the prime-field ecbench ladder of
//! `research/pkm_tower_ecbench_20261013/`: for each requested size, a prime
//! `p = c·2^s + 1` of `bits` bits and a curve `y² = x³ + ax + b` over it
//! with `#E = h·r`, `r` prime and `h ≤ max_cofactor`, found from a seed;
//! printed as `ecbench` `prime_explicit` curve specs (one JSON object a
//! line) with a generator of the order-`r` subgroup.  Deterministic in
//! the seed.  Point counting is exhaustive (`PrimeCurve::point_count`),
//! which is why the sizes stay at or below 28 bits.
//!
//! ```text
//! cargo run --release --example prime_tower_curves -- --bits 20,22,24,26,28 --two-adicity 12 --seed 20261013 --max-cofactor 4
//! ```

use std::env;

use crypto_lib::cryptanalysis::challenge_ladder::is_prime_u64;
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps, PrimeCurve, PrimePoint};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut bits: Vec<u32> = vec![20, 22, 24, 26, 28];
    let mut s = 12u32;
    let mut seed = 20261013u64;
    let mut max_h = 4u64;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--bits" => {
                i += 1;
                bits = args[i].split(',').map(|v| v.parse().unwrap()).collect();
            }
            "--two-adicity" => {
                i += 1;
                s = args[i].parse().unwrap();
            }
            "--seed" => {
                i += 1;
                seed = args[i].parse().unwrap();
            }
            "--max-cofactor" => {
                i += 1;
                max_h = args[i].parse().unwrap();
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    for b in bits {
        let mut rng = StdRng::seed_from_u64(seed ^ ((b as u64) << 40));
        // p = c·2^s + 1 with c odd and p of exactly b bits.
        let p = loop {
            let c_bits = b - s;
            let c = (rng.gen::<u64>() & ((1u64 << c_bits) - 1)) | (1u64 << (c_bits - 1)) | 1;
            let p = (c << s) + 1;
            if is_prime_u64(p) {
                break p;
            }
        };
        // Candidates in batches of 8, counted in parallel; the first that
        // qualifies in candidate order wins, so the result is deterministic.
        let found = loop {
            let cands: Vec<(u64, u64)> = (0..8)
                .map(|_| (rng.gen_range(0..p), rng.gen_range(1..p)))
                .collect();
            let orders: Vec<Option<(u64, u64, u64, u64)>> = cands
                .par_iter()
                .map(|&(a, bb)| {
                    let disc = (4 * (a as u128).pow(3) + 27 * (bb as u128).pow(2)) % p as u128;
                    if disc == 0 {
                        return None;
                    }
                    let curve = PrimeCurve { p, a, b: bb };
                    let n = curve.point_count();
                    (1..=max_h)
                        .find(|h| n.is_multiple_of(*h) && is_prime_u64(n / h) && n / h > max_h)
                        .map(|h| (a, bb, n, h))
                })
                .collect();
            if let Some(hit) = orders.into_iter().flatten().next() {
                break hit;
            }
        };
        let (a, bb, n, h) = found;
        let curve = PrimeCurve { p, a, b: bb };
        let r = n / h;
        let mut ops = GroupOps::default();
        let g = loop {
            let x = rng.gen_range(0..p);
            if let Some(y) = curve.sqrt(curve.rhs(x)) {
                let q = curve.mul(&mut ops, PrimePoint::affine(x, y), h);
                if !curve.is_identity(&q) {
                    break q;
                }
            }
        };
        assert!(
            curve.is_identity(&curve.mul(&mut ops, g, r)),
            "generator order"
        );
        println!(
            "{{\"kind\": \"prime_explicit\", \"p\": {p}, \"a\": {a}, \"b\": {bb}, \"group_order\": {n}, \"r\": {r}, \"gx\": {}, \"gy\": {}}}",
            g.x, g.y
        );
        eprintln!(
            "bits {b}: p = {p} (2-adicity {}), #E = {n} = {h}·{r}",
            (p - 1).trailing_zeros()
        );
    }
}
