//! Prime-order `a = −3` curves for the prime ladder of
//! `research/ic_frontier_campaign_20261015/`: for each requested size the
//! deterministic P-256-class rung of [`prime_fast::find_a3_curve`] (`p` the
//! largest prime `≡ 3 (mod 4)` below `2^bits`, `b` the smallest giving prime
//! order, order by Hasse-interval BSGS), printed as one `ecbench`
//! `prime_explicit` curve spec per line.
//!
//! ```text
//! cargo run --release --example ic_frontier_prime_curves -- --bits 20,24,28,32,36,40
//! ```
use std::env;

use crypto_lib::cryptanalysis::prime_fast::find_a3_curve;

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut bits: Vec<u32> = vec![20, 24, 28, 32, 36, 40];
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--bits" => {
                i += 1;
                bits = args[i].split(',').map(|v| v.parse().unwrap()).collect();
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    for b in bits {
        let (curve, g) = find_a3_curve(b, 3).unwrap_or_else(|| panic!("no curve at {b} bits"));
        let p = curve.f.p;
        let (gx, gy) = curve.canonical(g);
        println!(
            "{{\"kind\": \"prime_explicit\", \"p\": {p}, \"a\": {}, \"b\": {}, \"group_order\": {n}, \"r\": {n}, \"gx\": {gx}, \"gy\": {gy}}}",
            p - 3,
            curve.f.from_mont(curve.b),
            n = curve.n,
        );
    }
}
