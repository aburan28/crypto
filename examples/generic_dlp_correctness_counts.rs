//! Before/after counts for `research/generic_dlp_correctness_20261017`.
//!
//! Uses only API present at `b1bf5b893` (`pollard_collab::{demo_curve,
//! JobSpec, run_walker, WalkerOutcome}`, `pohlig_hellman::pohlig_hellman_curve`),
//! so the same source builds on `main` (before) and on the fix branch
//! (after).  Everything printed is a deterministic count; nothing is timed.
//!
//! ```sh
//! cargo run --release --example generic_dlp_correctness_counts -- <label>
//! ```
//!
//! Prints JSON lines:
//!
//! - `collab`: per demo curve and mode, how many of `walkers` walkers reached
//!   a distinguished point and how many were abandoned at the step cap (dead
//!   trails), their steps, and whether every DP record verifies.
//! - `ph`: per target on the 47-bit smooth-order curve the pinned walk test
//!   uses (`3 · 19⁴ · 151 · 181 · 13171`) and per smoothness bound, the
//!   recovered logarithm, whether it equals the planted one, and the
//!   report's `total_steps`.

use num_bigint::BigUint;

use crypto_lib::cryptanalysis::pohlig_hellman::pohlig_hellman_curve;
use crypto_lib::cryptanalysis::pollard_collab::{demo_curve, run_walker, JobSpec, WalkerOutcome};
use crypto_lib::ecc::CurveParams;

fn ph_curve() -> CurveParams {
    CurveParams {
        name: "pinned-j0-smooth47",
        p: BigUint::from(140_737_529_461_027u64),
        a: BigUint::from(0u32),
        b: BigUint::from(4u32),
        gx: BigUint::from(16_535_943_672_630u64),
        gy: BigUint::from(38_148_592_391_465u64),
        n: BigUint::from(140_737_531_856_763u64),
        h: 1,
    }
}

fn main() {
    let label = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "unlabelled".into());

    for (curve_name, walkers) in [("demo-small", 512u64), ("demo-mid", 512), ("demo-32", 256)] {
        let curve = demo_curve(curve_name).expect("demo curve");
        let q = curve
            .generator()
            .scalar_mul(&(&curve.n * 2u32 / 3u32), &curve.a_fe());
        for negation in [false, true] {
            let mut spec =
                JobSpec::new(&curve, &q, "correctness-counts", 0xC0_11AB).expect("job spec");
            spec.negation_map = negation;
            let dp_bits = spec.dp_bits;
            let ctx = spec.build().expect("job context");
            let (mut dps, mut dead, mut dp_steps, mut dead_steps, mut verified) =
                (0u64, 0u64, 0u64, 0u64, 0u64);
            for i in 0..walkers {
                match run_walker(&ctx, i) {
                    WalkerOutcome::Dp(r) => {
                        dps += 1;
                        dp_steps += r.steps;
                        if r.verify(&ctx).is_ok() {
                            verified += 1;
                        }
                    }
                    WalkerOutcome::DeadTrail { steps } => {
                        dead += 1;
                        dead_steps += steps;
                    }
                }
            }
            println!(
                "{{\"label\":\"{label}\",\"kind\":\"collab\",\"curve\":\"{curve_name}\",\"n_bits\":{},\"negation_map\":{negation},\"dp_bits\":{dp_bits},\"walkers\":{walkers},\"dps\":{dps},\"dead_trails\":{dead},\"dp_steps\":{dp_steps},\"dead_steps\":{dead_steps},\"dps_verified\":{verified}}}",
                curve.n.bits()
            );
        }
    }

    let curve = ph_curve();
    let g = curve.generator();
    let a = curve.a_fe();
    let mut x = 0x9E37_79B9_7F4A_7C15u64;
    for t in 0..24u32 {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        let d = BigUint::from(x) % &curve.n;
        let q = g.scalar_mul(&d, &a);
        for bound in [1u64 << 14, 1 << 17] {
            let r = pohlig_hellman_curve(&curve, &g, &q, &curve.n, bound);
            let ok = r.recovered_d.as_ref() == Some(&d);
            println!(
                "{{\"label\":\"{label}\",\"kind\":\"ph\",\"curve\":\"{}\",\"target\":{t},\"bound\":{bound},\"recovered\":{},\"correct\":{ok},\"total_steps\":{}}}",
                curve.name,
                r.recovered_d.map_or("null".to_string(), |v| format!("\"{v}\"")),
                r.total_steps
            );
        }
    }
}
