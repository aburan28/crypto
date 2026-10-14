//! lit2026 W11 (e): the Dallaire-Demers–Doyle–Foo `secpkk1` ECDLP challenge
//! ladder (arXiv 2508.14011, Algorithm 3), derived in Rust and solved by the
//! repository's parallel distinguished-point rho with the negation map
//! ([`crypto_lib::cryptanalysis::prime_fast::rho_parallel`]).
//!
//! ```sh
//! cargo run --release --example lit2026_w11_challenge_ladder -- derive 30 60 2
//! cargo run --release --example lit2026_w11_challenge_ladder -- solve 30 60 2 <reps>
//! ```
//!
//! `derive lo hi step` prints one JSON line per level (prime, order, the
//! two NUMS points in SEC1, embedding degree, the curve as BSGS-timing
//! arguments).  `solve lo hi step reps` runs rho `reps` times per level
//! (seeds 1..=reps), verifies `κQ = P`, and prints the iterations against
//! `√(πn/4)`, the negation-map expectation.

use std::time::Instant;

use crypto_lib::cryptanalysis::challenge_ladder::{self, derive};
use crypto_lib::cryptanalysis::prime_fast::{rho_parallel, FastCurve, Pt, RhoConfig};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mode = args.first().map(String::as_str).unwrap_or("derive");
    let num = |i: usize, d: u32| -> u32 {
        args.get(i)
            .map(|s| s.parse().expect("integer argument"))
            .unwrap_or(d)
    };
    let (lo, hi, step) = (num(1, 30), num(2, 60), num(3, 2).max(1));
    let reps = num(4, 4) as u64;
    let mut bits = lo;
    while bits <= hi {
        let t0 = Instant::now();
        let c = derive(bits).expect("instance");
        let derive_s = t0.elapsed().as_secs_f64();
        let r = c
            .embedding_degree_le_64
            .map(|v| v.to_string())
            .unwrap_or_else(|| "\">64\"".into());
        if mode == "derive" {
            println!(
                "{{\"bits\":{bits},\"p\":{},\"n\":{},\"P\":\"{}\",\"Q\":\"{}\",\"Px\":{},\"Py\":{},\"Qx\":{},\"Qy\":{},\"delta_P\":{},\"delta_Q\":{},\"primes_rejected\":{},\"embedding_degree\":{r},\"derive_s\":{derive_s:.3}}}",
                c.p,
                c.n,
                c.sec1(c.p_pt),
                c.sec1(c.q_pt),
                c.p_pt.0,
                c.p_pt.1,
                c.q_pt.0,
                c.q_pt.1,
                c.delta_p,
                c.delta_q,
                c.primes_rejected
            );
        } else {
            let curve = FastCurve::new(&format!("secp{bits}k1"), c.p, 0, 7, c.n);
            let f = curve.f;
            let g = Pt {
                x: f.to_mont(c.q_pt.0),
                y: f.to_mont(c.q_pt.1),
            };
            let q = Pt {
                x: f.to_mont(c.p_pt.0),
                y: f.to_mont(c.p_pt.1),
            };
            let expected = (std::f64::consts::PI * c.n as f64 / 4.0).sqrt();
            for seed in 1..=reps {
                let mut cfg = RhoConfig::for_bits(64 - c.n.leading_zeros(), true);
                cfg.seed = seed;
                let res = rho_parallel(&curve, g, q, &cfg).expect("rho solved");
                let kappa = res.log;
                let ok = challenge_ladder::mul(c.p, Some(c.q_pt), kappa) == Some(c.p_pt);
                assert!(ok, "κQ must equal P");
                let overhead =
                    (cfg.threads * cfg.walkers) as f64 * 2f64.powi(cfg.dp_bits as i32) / expected;
                println!(
                    "{{\"bits\":{bits},\"p\":{},\"n\":{},\"seed\":{seed},\"kappa\":{kappa},\"verified\":{ok},\"iterations\":{},\"expected\":{expected:.1},\"ratio\":{:.4},\"dps\":{},\"threads\":{},\"walkers_per_thread\":{},\"dp_bits\":{},\"dp_overhead_est\":{overhead:.4},\"wall_s\":{:.3},\"steps_per_s\":{:.0}}}",
                    c.p,
                    c.n,
                    res.total_steps,
                    res.total_steps as f64 / expected,
                    res.distinguished_points,
                    res.threads,
                    res.walkers_per_thread,
                    res.dp_bits,
                    res.wall_secs,
                    res.steps_per_sec
                );
            }
        }
        bits += step;
    }
}
