//! Wall time of the single-word parallel Pollard rho
//! ([`crypto_lib::cryptanalysis::prime_fast::rho_parallel`]) on one
//! public known-answer target, for the baseline / candidate arms of
//! `research/generic_dlp_memory_parallelism_20261017`.
//!
//! Uses only API present at the study's baseline (`b1bf5b893`), so the
//! same source builds there and on the candidate branch.  The curve is
//! [`find_a3_curve`]`(bits, 3)`, the `a = −3` prime-order curve the
//! `prime_ecdlp_fast_bench` ladder uses above 28 bits; generating it and
//! the target is outside every timer.
//!
//! ```sh
//! cargo run --release --example generic_dlp_rho_time -- <bits> <threads> <seed> <label>
//! ```
//!
//! Prints one JSON line:
//!
//! - `solve_s`: the whole `rho_parallel` call, which is what a caller
//!   waits for: the `r`-adding table (1024 [`FastCurve::combine`] calls),
//!   every walk start (one `combine` each), the walk, and the collision's
//!   verification.
//! - `setup_s`: the same number of `combine` calls as the table and the
//!   initial walk starts, timed on their own after the solve (a derived
//!   estimate of the start-up share, not a phase of the solve itself).
//! - the counted outputs, which must be identical between arms at one
//!   thread: `total_steps`, `distinguished_points`, the recovered log and
//!   its check `[log]G = Q`.

use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crypto_lib::cryptanalysis::prime_fast::{find_a3_curve, rho_parallel, RhoConfig};

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    assert!(a.len() >= 4, "usage: bits threads seed label");
    let bits: u32 = a[0].parse().expect("bits");
    let threads: usize = a[1].parse().expect("threads");
    let seed: u64 = a[2].parse().expect("seed");
    let label = &a[3];

    let (fc, g) = find_a3_curve(bits, 3).expect("a = -3 prime-order curve");
    let mut rng = StdRng::seed_from_u64(0x6765_6e64_6c70 ^ seed ^ ((bits as u64) << 32));
    let k = rng.gen_range(1..fc.n);
    let q = fc.mul(Some(g), k).expect("target is finite");

    let mut cfg = RhoConfig::for_bits(bits, true);
    cfg.threads = threads.max(1);
    cfg.seed = seed;

    let t0 = Instant::now();
    let res = rho_parallel(&fc, g, q, &cfg).expect("rho solves");
    let solve_s = t0.elapsed().as_secs_f64();
    let verified = res.log == k && fc.mul(Some(g), res.log) == Some(q);
    assert!(verified, "recovered log must verify");

    // The start-up work, repeated on its own: one `combine` per table
    // entry and per initial walk start, scalars drawn like rho draws them.
    let calls = (1usize << cfg.table_bits) + cfg.threads * cfg.walkers;
    let mut srng = StdRng::seed_from_u64(seed ^ 0x5e70_0000);
    let scalars: Vec<(u64, u64)> = (0..calls)
        .map(|_| (srng.gen_range(1..fc.n), srng.gen_range(1..fc.n)))
        .collect();
    let t1 = Instant::now();
    let mut sink = 0u64;
    for &(u, v) in &scalars {
        if let Some(p) = fc.combine(g, q, u, v) {
            sink ^= p.x;
        }
    }
    let setup_s = t1.elapsed().as_secs_f64();

    let (p, a3, b3) = (fc.f.p, fc.f.from_mont(fc.a), fc.f.from_mont(fc.b));
    let (gx, gy) = fc.canonical(g);
    println!(
        "{{\"label\":\"{label}\",\"bits\":{bits},\"p\":{p},\"a\":{a3},\"b\":{b3},\"n\":{},\"gx\":{gx},\"gy\":{gy},\"seed\":{seed},\"threads\":{},\"walkers_per_thread\":{},\"dp_bits\":{},\"table_bits\":{},\"k\":{k},\"log\":{},\"verified\":{verified},\"total_steps\":{},\"distinguished_points\":{},\"solve_s\":{solve_s:.6},\"setup_calls\":{calls},\"setup_s\":{setup_s:.6},\"sink\":{}}}",
        fc.n,
        cfg.threads,
        cfg.walkers,
        cfg.dp_bits,
        cfg.table_bits,
        res.log,
        res.total_steps,
        res.distinguished_points,
        sink & 1
    );
}
