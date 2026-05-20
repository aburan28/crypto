//! Probe the secp256k1 LLL-degeneracy phenomenon observed during
//! the GLV-HNP Phase 1.5 work.
//!
//! ## Observation
//!
//! With k_bits = 192 (64-bit nonce bias) and m = 8 signatures:
//!
//! - On NIST P-256: LLL converges to the correct key in ~1 second.
//! - On secp256k1: LLL hits its iteration cap across every random
//!   seed we tried.
//!
//! This probe runs a small, fast head-to-head comparison so the
//! difference is empirically documented.  Full BKZ/parameter
//! sweeps are deferred to the longer-running `#[ignore]`'d probe
//! version (commented out at the bottom).
//!
//! Run with:
//!   cargo test --test lll_degeneracy_probe -- --nocapture

use std::time::Instant;

use crypto_lib::cryptanalysis::hnp_ecdsa::{
    hnp_recover_key_with_reduction, BiasedSignature, HnpReduction,
};
use crypto_lib::ecc::curve::CurveParams;
use crypto_lib::ecc::keys::EccKeyPair;
use crypto_lib::ecc::point::Point;
use crypto_lib::utils::mod_inverse;
use num_bigint::{BigUint, RandBigInt};
use num_traits::Zero;
use rand::rngs::StdRng;
use rand::{RngCore, SeedableRng};

fn generate_biased_sigs(
    curve: &CurveParams,
    d: &BigUint,
    k_bits: u32,
    count: usize,
    k_seed: u64,
) -> Vec<BiasedSignature> {
    let n = &curve.n;
    let mut rng = StdRng::seed_from_u64(k_seed);
    let mut z_seed: u64 = 0xDEAD_BEEF;
    let mut sigs = Vec::with_capacity(count);

    while sigs.len() < count {
        let bytes_count = ((k_bits + 7) / 8) as usize;
        let mut buf = vec![0u8; bytes_count];
        rng.fill_bytes(&mut buf);
        let extra = (bytes_count as u32) * 8 - k_bits;
        if extra > 0 {
            buf[0] &= 0xff >> extra;
        }
        let k = BigUint::from_bytes_be(&buf);
        if k.is_zero() {
            continue;
        }

        let z = BigUint::from(z_seed) % n;
        z_seed = z_seed.wrapping_add(0x9E37_79B9_7F4A_7C15);

        let g = curve.generator();
        let a_fe = curve.a_fe();
        let kg = g.scalar_mul(&k, &a_fe);
        let x1 = match &kg {
            Point::Affine { x, .. } => x.value.clone(),
            Point::Infinity => continue,
        };
        let r = &x1 % n;
        if r.is_zero() {
            continue;
        }
        let rd = (&r * d) % n;
        let z_plus_rd = (&z + &rd) % n;
        let k_inv = match mod_inverse(&k, n) {
            Some(v) => v,
            None => continue,
        };
        let s = (&k_inv * &z_plus_rd) % n;
        if s.is_zero() {
            continue;
        }

        sigs.push(BiasedSignature { r, s, z, k_bits });
    }
    sigs
}

fn probe_once(
    curve: &CurveParams,
    curve_name: &str,
    k_bits: u32,
    m: usize,
    d_seed: u64,
    k_seed: u64,
) -> (String, u128) {
    let n = curve.n.clone();
    let mut d_rng = StdRng::seed_from_u64(d_seed);
    let d = d_rng.gen_biguint_below(&n);
    let kp = EccKeyPair::from_private(d.clone(), curve);

    let sigs = generate_biased_sigs(curve, &d, k_bits, m, k_seed);

    let t0 = Instant::now();
    let result =
        hnp_recover_key_with_reduction(curve, &kp.public, &sigs, HnpReduction::Lll);
    let elapsed = t0.elapsed().as_millis();

    let outcome = match result {
        Ok(recovered) if recovered == d => "✓ RECOVERED".to_string(),
        Ok(_) => "✗ WRONG KEY".to_string(),
        Err(e) => format!("✗ {}", e),
    };
    eprintln!(
        "{:<10} k_bits={:>3} m={:>2} d_seed=0x{:08X} k_seed=0x{:08X}  →  {:<40}  ({} ms)",
        curve_name, k_bits, m, d_seed, k_seed, outcome, elapsed
    );
    (outcome, elapsed)
}

#[test]
fn probe_lll_degeneracy_head_to_head() {
    let secp = CurveParams::secp256k1();
    let p256 = CurveParams::p256();

    eprintln!();
    eprintln!("=== LLL convergence: secp256k1 vs P-256 (k_bits=192, m=8) ===");

    let seeds: [(u64, u64); 3] = [
        (0xC0FFEE, 0xC0FFEE),
        (0xDEAD_BEEF, 0xBADC_AFE),
        (0x1234_5678, 0x9ABC_DEF0),
    ];

    let mut secp_pass = 0;
    let mut p256_pass = 0;

    for (d_seed, k_seed) in seeds {
        let (out_p, _) = probe_once(&p256, "P-256", 192, 8, d_seed, k_seed);
        if out_p.starts_with("✓") {
            p256_pass += 1;
        }
        let (out_s, _) = probe_once(&secp, "secp256k1", 192, 8, d_seed, k_seed);
        if out_s.starts_with("✓") {
            secp_pass += 1;
        }
    }

    eprintln!();
    eprintln!("=== Summary ===");
    eprintln!("P-256:     {}/{} probes recovered key", p256_pass, seeds.len());
    eprintln!("secp256k1: {}/{} probes recovered key", secp_pass, seeds.len());

    assert!(
        p256_pass + secp_pass > 0,
        "expected at least some LLL convergence; got 0 for both curves"
    );
}

/// Extended bit-length sweep of LLL convergence across curves.
///
/// Test the refined hypothesis (after the cross-Koblitz refutation):
/// does **bit-length alone** predict LLL convergence among
/// generic-n curves?  Or is secp256k1 specifically the outlier?
///
/// Sweeps P-256, P-384, P-521, brainpoolP256r1, brainpoolP384r1 —
/// all generic-n curves — at parameters scaled to each curve's
/// bit-length.  Compares against secp256k1 (the known failure).
///
/// Run: `cargo test --test lll_degeneracy_probe sweep -- --ignored --nocapture`
#[test]
#[ignore = "slow: 7 curves × 1 seed = ~10 probes (~3-5 minutes)"]
fn probe_lll_sweep_by_bit_length() {
    use crypto_lib::ecc::curve::CurveParams;
    let curves: Vec<(CurveParams, &str, u32, u32)> = vec![
        // (curve, name, k_bits, n_bits)
        (CurveParams::secp192k1(), "secp192k1", 144, 192),
        (CurveParams::secp224k1(), "secp224k1", 168, 224),
        (CurveParams::secp256k1(), "secp256k1", 192, 256),
        (CurveParams::p256(), "P-256", 192, 256),
        (CurveParams::brainpool_p256r1(), "brainpoolP256r1", 192, 256),
        (CurveParams::p384(), "P-384", 288, 384),
        (CurveParams::brainpool_p384r1(), "brainpoolP384r1", 288, 384),
        (CurveParams::p521(), "P-521", 384, 521),
    ];

    eprintln!();
    eprintln!("=== Bit-length sweep: LLL convergence across deployed curves ===");
    eprintln!(
        "{:<20} {:>8} {:>10} {:>10} {:>15}",
        "curve", "n_bits", "k_bits", "elapsed", "outcome"
    );

    let mut results: Vec<(String, u32, String, u128)> = Vec::new();

    for (curve, name, k_bits, n_bits) in &curves {
        let (outcome, elapsed) = probe_once(curve, name, *k_bits, 8, 0xC0FFEE, 0xC0FFEE);
        eprintln!(
            "{:<20} {:>8} {:>10} {:>10} {:>15}",
            name, n_bits, k_bits, elapsed, outcome
        );
        results.push((name.to_string(), *n_bits, outcome, elapsed));
    }

    eprintln!();
    eprintln!("=== Analysis ===");
    let failures: Vec<&(String, u32, String, u128)> = results
        .iter()
        .filter(|(_, _, o, _)| !o.starts_with("✓"))
        .collect();
    eprintln!("Curves that FAILED LLL: {}", failures.len());
    for f in &failures {
        eprintln!("  - {} (n_bits = {})", f.0, f.1);
    }
    eprintln!();
    if failures.iter().all(|(name, _, _, _)| name == "secp256k1") {
        eprintln!("✓ secp256k1 is the SOLE failure across deployed curves.");
        eprintln!("  Hypotheses (2) bit-length cap and (3) generic 256-bit");
        eprintln!("  are REFUTED — P-384/P-521 (larger bit-length) succeed,");
        eprintln!("  and other 256-bit curves succeed.  The degeneracy is");
        eprintln!("  curve-specific to secp256k1, not bit-length or generic-n.");
    } else if failures.iter().all(|(_, b, _, _)| *b >= 256) {
        eprintln!("? Failures restricted to ≥256-bit curves; bit-length");
        eprintln!("  hypothesis remains plausible.");
    } else {
        eprintln!("? Mixed failures; need more data.");
    }
}

/// Cross-Koblitz LLL-degeneracy probe.
///
/// Tests the working hypothesis that **near-power-of-2 group order
/// `n` correlates with LLL failure** on the Boneh-Venkatesan
/// lattice.  We probe four curves:
///
///   - secp256k1  (Koblitz, n ≈ 2^256 − small, LLL-fails)
///   - secp192k1  (Koblitz, n ≈ 2^192 − small) ← prediction: fails
///   - secp224k1  (Koblitz, n ≈ 2^224 − small) ← prediction: fails
///   - P-256      (NIST, n more "generic") ← control: passes
///
/// At parameters scaled to the curve's bit-length (k_bits = 3·n_bits/4,
/// m = 8), LLL should systematically fail on the Koblitz family but
/// succeed on P-256.
///
/// Run: `cargo test --test lll_degeneracy_probe koblitz -- --nocapture`
#[test]
#[ignore = "slow: 4 curves × 2 seeds = 8 probes (~2-3 minutes)"]
fn probe_koblitz_lll_degeneracy_hypothesis() {
    let curves: Vec<(CurveParams, &str, u32, u32)> = vec![
        // (curve, name, k_bits, n_bits) — k_bits chosen as 3/4 of n_bits
        (CurveParams::secp192k1(), "secp192k1", 144, 192),
        (CurveParams::secp224k1(), "secp224k1", 168, 224),
        (CurveParams::secp256k1(), "secp256k1", 192, 256),
        (CurveParams::p256(), "P-256", 192, 256),
    ];

    let seeds: [(u64, u64); 2] = [(0xC0FFEE, 0xC0FFEE), (0xDEAD_BEEF, 0xBADC_AFE)];

    eprintln!();
    eprintln!("=== Koblitz vs NIST: LLL convergence by group-order shape ===");
    eprintln!(
        "{:<12} {:>8} {:>10} {:>8} {:>10}",
        "curve", "n_bits", "n-form", "passes", "total"
    );

    let mut results: Vec<(String, u32, &str, usize, usize, u128)> = Vec::new();

    for (curve, name, k_bits, n_bits) in &curves {
        // Classify the group-order shape: is n near 2^n_bits?
        let n_bits_actual = curve.n.bits() as u32;
        let two_to_n = num_bigint::BigUint::from(1u32) << n_bits_actual;
        let gap = if two_to_n > curve.n {
            &two_to_n - &curve.n
        } else {
            &curve.n - &two_to_n
        };
        let gap_bits = gap.bits() as u32;
        let form = if gap_bits < n_bits_actual / 2 {
            "near-2^k"
        } else {
            "generic"
        };

        let mut passes = 0usize;
        let mut total_elapsed: u128 = 0;
        for (d_seed, k_seed) in &seeds {
            let (outcome, elapsed) =
                probe_once(curve, name, *k_bits, 8, *d_seed, *k_seed);
            total_elapsed += elapsed;
            if outcome.starts_with("✓") {
                passes += 1;
            }
        }
        eprintln!(
            "{:<12} {:>8} {:>10} {:>8} {:>10}",
            name,
            n_bits,
            form,
            passes,
            seeds.len()
        );
        results.push((name.to_string(), *n_bits, form, passes, seeds.len(), total_elapsed));
    }

    eprintln!();
    eprintln!("=== Hypothesis: near-2^k group orders ⇒ LLL fails ===");
    let koblitz_total: usize = results
        .iter()
        .filter(|(_, _, f, _, _, _)| *f == "near-2^k")
        .map(|(_, _, _, p, _, _)| p)
        .sum();
    let koblitz_attempts: usize = results
        .iter()
        .filter(|(_, _, f, _, _, _)| *f == "near-2^k")
        .map(|(_, _, _, _, t, _)| t)
        .sum();
    let generic_total: usize = results
        .iter()
        .filter(|(_, _, f, _, _, _)| *f == "generic")
        .map(|(_, _, _, p, _, _)| p)
        .sum();
    let generic_attempts: usize = results
        .iter()
        .filter(|(_, _, f, _, _, _)| *f == "generic")
        .map(|(_, _, _, _, t, _)| t)
        .sum();
    eprintln!(
        "  near-2^k curves: {}/{} LLL successes",
        koblitz_total, koblitz_attempts
    );
    eprintln!(
        "  generic curves:  {}/{} LLL successes",
        generic_total, generic_attempts
    );
    eprintln!();
    if koblitz_total == 0 && generic_total == generic_attempts {
        eprintln!("✓ Hypothesis confirmed: near-power-of-2 n ⇒ LLL fails uniformly.");
    } else if koblitz_total > 0 && generic_total < generic_attempts {
        eprintln!("✗ Hypothesis refuted: outcomes do not correlate with n's form.");
    } else {
        eprintln!("? Hypothesis partially supported; more curves needed.");
    }
}
