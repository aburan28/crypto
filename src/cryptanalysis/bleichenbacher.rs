//! Bleichenbacher's FFT-based bias detector for ECDSA nonces.
//!
//! Daniel Bleichenbacher, "On the generation of one-time keys in
//! DL signature schemes" (2000 unpublished talk; canonical refs:
//! de Mulder–Hutter–Marson–Pearson CHES 2014, Aranha–Tibouchi–Yarom
//! "LadderLeak" CCS 2020).  The lattice-attack approach (HNP via
//! LLL) requires `m · b > n_bits` signatures with `b` bits of bias.
//! Bleichenbacher's spectral method works *below* this threshold —
//! it can detect bias of less than 1 bit per signature given
//! enough samples.
//!
//! # The algorithm
//!
//! Given many ECDSA signatures `(r_i, s_i, z_i)` with possibly-
//! biased nonces `k_i`, define the HNP form:
//!
//! ```text
//! h_i = s_i⁻¹ · r_i  (mod n)
//! t_i = s_i⁻¹ · z_i  (mod n)
//! ```
//!
//! such that `k_i = t_i + h_i · d (mod n)`.  For the **correct**
//! private key `d`, the values `(t_i + h_i · d) mod n` cluster
//! near 0 (since each `k_i < 2^k_bits ≪ n`).  For a wrong `d`,
//! they look uniform.
//!
//! Bleichenbacher's bias function:
//!
//! ```text
//! Z(d) = Σᵢ exp(2π·i · (t_i + h_i · d) / n)
//! ```
//!
//! - **Correct `d`**: angles cluster, `|Z(d)|` ≈ `m · sinc(π · 2^k_bits / n)`.
//! - **Wrong `d`**: random walk in C, `|Z(d)|` ≈ `√m`.
//!
//! For toy `n` (small enough to evaluate `Z(d)` for all candidates),
//! the peak in `|Z(d)|` reveals `d` directly.  For cryptographic
//! `n`, this requires the **lattice + FFT hybrid**: use lattice
//! reduction on `(h_i)` to find linear combinations whose
//! transformed `h` values fit in a small window `[0, L)`, then FFT
//! to find `d mod L`, repeat at multiple `L`'s, CRT-combine.  This
//! module ships the **toy-scale direct version** with documentation
//! on the lattice-reduction extension.
//!
//! # What this proves
//!
//! - At toy scale (`n ~ 2^{20}`), Bleichenbacher's spectral method
//!   detects nonce bias and identifies the correct `d` by peak-
//!   finding even when standard HNP would be borderline.
//! - The signal-to-noise ratio scales as `√m · sinc(π · 2^k_bits / n)`,
//!   so for fixed bias, more signatures push the peak above the noise
//!   floor.  This is the "fractional-bit" regime: arbitrarily small
//!   `bias_bits` can be detected with arbitrarily many signatures.
//!
//! # What this does NOT do
//!
//! Implement the lattice + FFT hybrid required for cryptographic-
//! size `n`.  That requires (a) lattice reduction to compress `h_i`
//! into a small window, (b) FFT over that window, (c) CRT
//! combination across many windows.  Each step alone is in this
//! library (LLL/BKZ in `cryptanalysis::lattice`); their composition
//! into a full Bleichenbacher attack is left as future work.

use crate::utils::mod_inverse;
use num_bigint::BigUint;
use num_traits::Zero;

/// One ECDSA signature in HNP-derived form.
#[derive(Clone, Debug)]
pub struct BleichenbacherSample {
    /// `t_i = s_i⁻¹ · z_i mod n`
    pub t: BigUint,
    /// `h_i = s_i⁻¹ · r_i mod n`
    pub h: BigUint,
}

/// Convert a raw ECDSA `(r, s, z)` triple to Bleichenbacher form.
pub fn signature_to_sample(
    r: &BigUint,
    s: &BigUint,
    z: &BigUint,
    n: &BigUint,
) -> Option<BleichenbacherSample> {
    if s.is_zero() || r.is_zero() || s >= n || r >= n {
        return None;
    }
    let s_inv = mod_inverse(s, n)?;
    Some(BleichenbacherSample {
        t: (&s_inv * z) % n,
        h: (&s_inv * r) % n,
    })
}

/// Result of a peak search.
#[derive(Clone, Debug)]
pub struct BleichenbacherPeak {
    /// Candidate `d` value.
    pub d: BigUint,
    /// `|Z(d)|` (peak magnitude).
    pub magnitude: f64,
    /// Sample noise floor (`√m`).
    pub noise_floor: f64,
    /// Peak SNR (`magnitude / noise_floor`).
    pub snr: f64,
}

/// Compute the Bleichenbacher bias function `|Z(d)|` for a given
/// candidate `d`.
///
/// `n` is the curve order.  Returns `|Z(d)|` as `f64`.  For the
/// correct `d`, this should be ≈ `m · sinc(π · 2^k_bits / n)` —
/// significantly above the `√m` noise floor.
pub fn bias_magnitude(samples: &[BleichenbacherSample], d: &BigUint, n: &BigUint) -> f64 {
    let n_f = biguint_to_f64(n);
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut re = 0.0_f64;
    let mut im = 0.0_f64;
    for s in samples {
        let arg_int = (&s.t + &s.h * d) % n;
        let arg = biguint_to_f64(&arg_int) / n_f;
        let theta = two_pi * arg;
        re += theta.cos();
        im += theta.sin();
    }
    (re * re + im * im).sqrt()
}

/// **Direct (toy-scale) Bleichenbacher attack**: scan all `d ∈ [0, n)`
/// and return the one with the highest `|Z(d)|`.
///
/// Cost: `O(n · m)` where `n` is the group order and `m` the number
/// of samples.  Suitable for `n ≤ ~2^{22}` (a few seconds on a
/// laptop).  For larger `n`, use the lattice + FFT hybrid (not
/// implemented here).
///
/// `snr_threshold` rejects "no bias detected" results: returns
/// `None` if the peak SNR (peak-magnitude / noise-floor) is below
/// the threshold.  Default heuristic: 3.0 (peak should be ≥ 3× the
/// noise floor to be a reliable signal).
pub fn bleichenbacher_direct(
    samples: &[BleichenbacherSample],
    n: &BigUint,
    snr_threshold: f64,
) -> Option<BleichenbacherPeak> {
    let m = samples.len();
    if m < 2 {
        return None;
    }
    let n_u64 = n.to_u64_digits().first().copied().unwrap_or(0);
    if n_u64 == 0 || n_u64 > (1u64 << 24) {
        return None; // beyond toy scale
    }
    let noise_floor = (m as f64).sqrt();

    let mut best_d = 0u64;
    let mut best_mag = 0.0_f64;
    for d_u in 0..n_u64 {
        let d = BigUint::from(d_u);
        let mag = bias_magnitude(samples, &d, n);
        if mag > best_mag {
            best_mag = mag;
            best_d = d_u;
        }
    }
    let snr = best_mag / noise_floor.max(1e-12);
    if snr < snr_threshold {
        return None;
    }
    Some(BleichenbacherPeak {
        d: BigUint::from(best_d),
        magnitude: best_mag,
        noise_floor,
        snr,
    })
}

/// Convert `BigUint` to `f64`, lossy for cryptographic-size values.
fn biguint_to_f64(x: &BigUint) -> f64 {
    let bits = x.bits() as i32;
    if bits < 53 {
        x.iter_u64_digits().next().unwrap_or(0) as f64
    } else {
        // Use leading 53 bits as mantissa.
        let bytes = x.to_bytes_be();
        let mut acc = 0.0_f64;
        for b in bytes.iter().take(8) {
            acc = acc * 256.0 + *b as f64;
        }
        acc * 2f64.powi((bits - 53).max(0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use num_bigint::RandBigInt;
    use rand::rngs::StdRng;
    use rand::{RngCore, SeedableRng};

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

    /// Synthesize a Bleichenbacher sample directly from
    /// `(t = k − h·d) mod n` rather than via full ECDSA — avoids
    /// the cost of point arithmetic at scale and lets us test
    /// the spectral method in isolation.
    fn synth_sample<R: RngCore>(
        rng: &mut R,
        n: &BigUint,
        d: &BigUint,
        k_bits: u32,
    ) -> BleichenbacherSample {
        let h = rng.gen_biguint_below(n);
        let k = biased_nonce(rng, k_bits);
        // t = (k − h·d) mod n
        let hd = (&h * d) % n;
        let t = if k >= hd {
            (&k - &hd) % n
        } else {
            n - ((&hd - &k) % n)
        };
        BleichenbacherSample { t, h }
    }

    /// Verify the bias-magnitude function: at the correct `d`,
    /// `|Z(d)|` should be `≈ m · sinc(π · 2^k_bits / n)`; at random
    /// `d`, `≈ √m`.
    #[test]
    fn bias_magnitude_peaks_at_correct_d() {
        let n = BigUint::from(1009u32);
        let d_true = BigUint::from(314u32);
        let mut rng = StdRng::seed_from_u64(0xC0FFEEu64);
        let samples: Vec<BleichenbacherSample> = (0..200)
            .map(|_| synth_sample(&mut rng, &n, &d_true, 7)) // 7-bit k = 3-bit bias
            .collect();
        let mag_true = bias_magnitude(&samples, &d_true, &n);
        // Sample the magnitude at 5 wrong d values; max of those should
        // be well below the true-d magnitude.
        let mut wrong_max = 0.0_f64;
        for &d_wrong in &[0u32, 100, 500, 700, 900] {
            let m = bias_magnitude(&samples, &BigUint::from(d_wrong), &n);
            if m > wrong_max {
                wrong_max = m;
            }
        }
        assert!(
            mag_true > 2.0 * wrong_max,
            "expected mag_true > 2× max-wrong, got mag_true={} wrong_max={}",
            mag_true, wrong_max
        );
    }

    /// **Headline test**: direct Bleichenbacher attack on toy scale.
    /// At `n = 1009` (10-bit) with planted `d` and 3-bit nonce
    /// bias, the peak finder recovers `d`.
    #[test]
    fn bleichenbacher_direct_recovers_d_at_toy_scale() {
        let n = BigUint::from(1009u32);
        let d_true = BigUint::from(314u32);
        let mut rng = StdRng::seed_from_u64(0xBADCAFEu64);
        // 7-bit k = 3-bit bias on a 10-bit n
        let samples: Vec<BleichenbacherSample> = (0..400)
            .map(|_| synth_sample(&mut rng, &n, &d_true, 7))
            .collect();
        let peak = bleichenbacher_direct(&samples, &n, 3.0)
            .expect("peak should clear the 3.0 SNR threshold");
        assert_eq!(
            peak.d, d_true,
            "Bleichenbacher peak found d = {}, expected {}",
            peak.d, d_true
        );
        assert!(peak.snr > 3.0, "SNR too low: {}", peak.snr);
    }

    /// **Sub-bit bias regime**: with very weak bias (1 bit on a 10-bit
    /// curve, i.e., `k_bits = 9` on `n = 1009`) and many samples,
    /// the spectral method still finds `d`.
    #[test]
    fn bleichenbacher_handles_one_bit_bias() {
        let n = BigUint::from(1009u32);
        let d_true = BigUint::from(412u32);
        let mut rng = StdRng::seed_from_u64(0xFEEDFACEu64);
        // k_bits = 9 means k < 512; n = 1009 ≈ 2^10, so 1-bit bias
        let samples: Vec<BleichenbacherSample> = (0..2000)
            .map(|_| synth_sample(&mut rng, &n, &d_true, 9))
            .collect();
        let peak = bleichenbacher_direct(&samples, &n, 2.5)
            .expect("1-bit bias with 2000 samples should clear SNR 2.5");
        assert_eq!(peak.d, d_true);
    }

    /// **Negative control**: full-entropy (uniform) nonces produce
    /// no peak above the noise floor.
    #[test]
    fn bleichenbacher_no_bias_no_peak() {
        let n = BigUint::from(1009u32);
        let d_true = BigUint::from(500u32);
        let mut rng = StdRng::seed_from_u64(0xDEADu64);
        // k_bits = 10 ≈ log₂(n) → essentially uniform
        let samples: Vec<BleichenbacherSample> = (0..400)
            .map(|_| {
                // Use a fully uniform k.
                let h = rng.gen_biguint_below(&n);
                let k = rng.gen_biguint_below(&n); // FULL uniform
                let hd = (&h * &d_true) % &n;
                let t = if k >= hd {
                    (&k - &hd) % &n
                } else {
                    &n - ((&hd - &k) % &n)
                };
                BleichenbacherSample { t, h }
            })
            .collect();
        let result = bleichenbacher_direct(&samples, &n, 5.0);
        assert!(
            result.is_none(),
            "uniform nonces should not produce a high-SNR peak; got {:?}",
            result
        );
    }

    /// `signature_to_sample` rejects malformed inputs.
    #[test]
    fn rejects_malformed_signatures() {
        let n = BigUint::from(101u32);
        assert!(signature_to_sample(
            &BigUint::zero(),
            &BigUint::from(2u32),
            &BigUint::from(3u32),
            &n,
        )
        .is_none());
        assert!(signature_to_sample(
            &BigUint::from(2u32),
            &BigUint::zero(),
            &BigUint::from(3u32),
            &n,
        )
        .is_none());
    }

    /// Reject too-small sample counts.
    #[test]
    fn rejects_short_transcripts() {
        let n = BigUint::from(1009u32);
        let samples: Vec<BleichenbacherSample> = vec![BleichenbacherSample {
            t: BigUint::from(1u32),
            h: BigUint::from(2u32),
        }];
        assert!(bleichenbacher_direct(&samples, &n, 3.0).is_none());
    }
}
