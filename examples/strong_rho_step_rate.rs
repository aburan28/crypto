//! What the strong signed-Frobenius rho reference costs per step on the
//! m = 83 confidence-gate curve (AGENTS.md §8a), where a complete walk —
//! about 2^39 steps — cannot run.
//!
//! Two measurements, both on `E_0: y² + xy = x³ + 1` over
//! `GF(2)[z]/(z^83 + z^45 + z² + z + 1)`:
//!
//! - **Canonical forms.** A fixed sequence of states (a public start point
//!   stepped by a public jump) is canonicalized, and every result is
//!   hashed.  Two builds that print the same digest computed the same
//!   canonical point and the same orbit multiplier on every state.
//! - **A capped walk.** [`WideStrongRho::solve`] at the default lanes and
//!   distinguished-point bits with `step_cap_factor = 0`, which leaves the
//!   step cap at its floor of 10^6.  Its wall time is the cost of those
//!   steps; a walk from the same seeds takes the same path in any build
//!   whose canonical forms agree.
//!
//! The target is the public synthetic point `[0x5eed_0083]G`.  Nothing here
//! recovers a logarithm: the walk stops at the cap.
//!
//! `cargo run --release --example strong_rho_step_rate [states] [seed]`
//!
//! Prints one JSON object.
use std::hint::black_box;
use std::time::Instant;

use crypto_lib::binary_ecc::IrreduciblePoly;
use crypto_lib::cryptanalysis::koblitz_strong_rho::{
    RawPointG, StrongRhoCharges, StrongRhoParams, WalkStateG,
};
use crypto_lib::cryptanalysis::koblitz_wide::WideKoblitz;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::json;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let states: u64 = args.get(1).map_or(1 << 20, |s| s.parse().expect("states"));
    let seed: u64 = args
        .get(2)
        .map_or(0x0083_5eed, |s| s.parse().expect("seed"));

    let irr = IrreduciblePoly {
        degree: 83,
        low_terms: vec![0, 1, 2, 45],
    };
    let kc = WideKoblitz::new(0, &irr).expect("E_0 over GF(2^83)");
    let rho = kc.strong_rho();
    let g = rho.generator();

    // Canonical forms of a fixed walk of states.
    let mut point = rho.scalar_mul(g, 0x1234_5678_9abc_def0_u128 % kc.r);
    let jump = rho.scalar_mul(g, 0x0fed_cba9_8765_4321_u128 % kc.r);
    let mut walk = Vec::with_capacity(states as usize);
    for i in 0..states {
        point = rho.add(point, jump);
        walk.push(WalkStateG {
            point,
            a: u128::from(i) % kc.r,
            b: 1,
        });
    }
    let mut digest = blake3::Hasher::new();
    let t = Instant::now();
    for &state in &walk {
        let c = rho.canonicalize(black_box(state));
        if let RawPointG::Affine { x, y } = c.point {
            digest.update(&x.to_le_bytes());
            digest.update(&y.to_le_bytes());
        }
        digest.update(&c.a.to_le_bytes());
        digest.update(&c.b.to_le_bytes());
    }
    let canon_ns = t.elapsed().as_nanos();

    // A walk capped at 10^6 steps.
    let params = StrongRhoParams {
        step_cap_factor: 0,
        ..StrongRhoParams::default()
    };
    let mut charges = StrongRhoCharges::default();
    let jumps = rho.jumps(seed, &mut charges);
    let q = rho.scalar_mul(g, 0x5eed_0083);
    let mut start = StdRng::seed_from_u64(seed);
    let t = Instant::now();
    let outcome = rho.solve(q, &jumps, &mut start, &params, charges);
    let walk_ns = t.elapsed().as_nanos();

    println!(
        "{}",
        json!({
            "curve": kc.curve_id().slug,
            "n": kc.n,
            "states": states,
            "canonical_digest": digest.finalize().to_hex().to_string(),
            "canonicalize_ns_per_state": canon_ns as f64 / states as f64,
            "walk": {
                "lanes": params.lanes,
                "dp_bits": params.dp_bits,
                "step_cap": 1_000_000u64,
                "seed": seed,
                "reached_cap": outcome.is_none(),
                "wall_ns": walk_ns as u64,
            },
        })
    );
}
