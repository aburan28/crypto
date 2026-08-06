//! Empirical LLL recovery threshold for `hnp_recover_key` on P-256.
//!
//! Thread 23 (2026-08-06) showed that the Boneh–Venkatesan lattice built by
//! `hnp_ecdsa` encoded the nonce unknowns UNCENTERED (`k_i ∈ [0, K)`), which
//! costs a factor 2 in the target vector's length versus the centered
//! encoding (`E[x²] = K²/3` vs `K²/12`).  `hnp_ecdsa.rs` now centers.
//!
//! Measured P-256 thresholds (3 seeds/cell, smallest `m` recovering 3/3):
//!
//! | bias bits | uncentered | centered |
//! |-----------|-----------|----------|
//! | 64        | 5         | 5        |
//! | 32        | 9         | 9        |
//! | 16        | 18        | **17**   |
//!
//! The gain is smaller here than on the GLV-HNP Phase-2 lattice (which gets
//! 2–3×) because only the `m` nonce coordinates are centered: the `d`
//! coordinate and the constant Kannan anchor `n·2^l` are each the same
//! order as one nonce coordinate and are not centerable, so `‖w‖²` drops
//! from `(m/3 + 4/3)` to `(m/12 + 4/3)` rather than by a full factor 4.
//! At m=17 that is 1.60× in length — worth about one signature.  At 64 and
//! 32 bias bits one signature is a coarser step than the gain, so the
//! threshold does not move there.
//!
//! Run: `cargo test --test hnp_centering_threshold -- --ignored --nocapture`

use crypto_lib::cryptanalysis::hnp_ecdsa::{hnp_recover_key, BiasedSignature};
use crypto_lib::ecc::curve::CurveParams;
use crypto_lib::ecc::keys::EccKeyPair;
use crypto_lib::ecc::point::Point;
use crypto_lib::utils::mod_inverse;
use num_bigint::BigUint;
use num_traits::Zero;
use rand::rngs::StdRng;
use rand::{RngCore, SeedableRng};

fn sign_with_nonce(
    z: &BigUint,
    k: &BigUint,
    d: &BigUint,
    curve: &CurveParams,
) -> Option<(BigUint, BigUint)> {
    let g = curve.generator();
    let a = curve.a_fe();
    let kg = g.scalar_mul(k, &a);
    let x1 = match &kg {
        Point::Affine { x, .. } => x.value.clone(),
        Point::Infinity => return None,
    };
    let r = &x1 % &curve.n;
    if r.is_zero() {
        return None;
    }
    let rd = (&r * d) % &curve.n;
    let z_plus_rd = (z + &rd) % &curve.n;
    let k_inv = mod_inverse(k, &curve.n)?;
    let s = (&k_inv * &z_plus_rd) % &curve.n;
    if s.is_zero() {
        return None;
    }
    Some((r, s))
}

fn biased_nonce<R: RngCore>(rng: &mut R, k_bits: u32) -> BigUint {
    loop {
        let bytes = ((k_bits + 7) / 8) as usize;
        let mut buf = vec![0u8; bytes];
        rng.fill_bytes(&mut buf);
        let extra = (bytes as u32) * 8 - k_bits;
        if extra > 0 {
            buf[0] &= 0xff >> extra;
        }
        let k = BigUint::from_bytes_be(&buf);
        if !k.is_zero() {
            return k;
        }
    }
}

/// One (bias_bits, m, seed) trial.  Returns true on exact recovery.
fn trial(bias_bits: u32, m: usize, seed: u64) -> bool {
    let curve = CurveParams::p256();
    let n = curve.n.clone();
    let n_bits = n.bits() as u32;
    let k_bits = n_bits - bias_bits;
    let mut rng = StdRng::seed_from_u64(seed);
    let d = {
        let mut buf = [0u8; 32];
        rng.fill_bytes(&mut buf);
        BigUint::from_bytes_be(&buf) % (&n - 1u32) + 1u32
    };
    let kp = EccKeyPair::from_private(d.clone(), &curve);

    let mut sigs: Vec<BiasedSignature> = Vec::new();
    let mut z_seed = 0xDEAD_BEEF_u64 ^ seed;
    while sigs.len() < m {
        let z = BigUint::from(z_seed) % &n;
        z_seed = z_seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let k = biased_nonce(&mut rng, k_bits);
        if let Some((r, s)) = sign_with_nonce(&z, &k, &d, &curve) {
            sigs.push(BiasedSignature { r, s, z, k_bits });
        }
    }
    hnp_recover_key(&curve, &kp.public, &sigs).ok() == Some(d)
}

/// Print the success rate over a (bias_bits, m) grid.  Diagnostic only.
#[test]
#[ignore = "slow: dozens of LLL reductions on up to 34-dim P-256 lattices"]
fn lll_threshold_sweep() {
    const SEEDS: [u64; 3] = [1, 2, 3];
    println!("\nP-256, plain LLL (delta=0.75), 3 seeds per cell");
    println!("{:>10} {:>5} {:>7}", "bias_bits", "m", "recov");
    for &(bias, ms) in &[
        (64u32, &[4usize, 5, 6][..]),
        (32, &[8, 9, 10][..]),
        (16, &[16, 17, 18][..]),
    ] {
        for &m in ms {
            let w = SEEDS.iter().filter(|&&s| trial(bias, m, s)).count();
            println!("{:>10} {:>5} {:>7}", bias, m, format!("{}/3", w));
            use std::io::Write;
            std::io::stdout().flush().ok();
        }
    }
}

/// Regression guard for the centering change.
///
/// `(16, 17)` is the discriminating cell: measured 2026-08-06, the
/// uncentered lattice gets 1/3 there and the centered one 3/3.  `(64, 5)`
/// and `(32, 9)` are thresholds both formulations meet — they are here to
/// catch an outright break rather than a reverted centering.
#[test]
#[ignore = "slow: 9 LLL reductions on up to 19-dim P-256 lattices (~4 min)"]
fn centered_lattice_meets_threshold() {
    for &(bias, m) in &[(64u32, 5usize), (32, 9), (16, 17)] {
        let w = [1u64, 2, 3].iter().filter(|&&s| trial(bias, m, s)).count();
        assert_eq!(
            w, 3,
            "bias={bias} bits, m={m}: expected 3/3 recoveries, got {w}/3 \
             (has the nonce centering in hnp_ecdsa.rs been reverted?)"
        );
    }
}
