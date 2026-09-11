//! # Challenge 74 — Correlation Power Analysis on AES-128
//!
//! Simulated power-trace oracle.  Each trace has, for each AES key
//! byte position, a leakage value
//!
//! ```text
//!     leakage[i] = HW(SBox(pt[i] ⊕ K[i]))  +  Gaussian(σ)
//! ```
//!
//! Brier-Clavier-Olivier 2004's *Correlation Power Analysis*:
//! for each key-byte hypothesis `k̂`, compute the predicted
//! Hamming weight of `SBox(pt[i] ⊕ k̂)` over all traces, take
//! the Pearson correlation with the observed leakage at byte `i`,
//! and pick the `k̂` with the highest correlation.  The correct
//! key byte shows a sharp peak (≈ 0.8); wrong guesses sit near
//! the noise floor.
//!
//! This is *the* canonical embedded-systems attack and the reason
//! every modern smartcard ships masking + shuffling countermeasures.

use crate::cryptopals::Report;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

const AES_SBOX: [u8; 256] = [
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16,
];

fn hamming_weight(b: u8) -> u32 {
    b.count_ones()
}

/// Sample Gaussian `N(0, σ)` via the Box-Muller transform (uses
/// two uniform samples).  Avoids pulling in `rand_distr` for one
/// line of code.
fn gauss(rng: &mut StdRng, sigma: f64) -> f64 {
    let u1: f64 = rng.gen_range(1e-12..1.0);
    let u2: f64 = rng.gen_range(0.0..1.0);
    sigma * (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// Capture `n_traces` simulated power traces.  Each trace is a
/// 16-element `f64` vector; trace `t`'s `i`th component is the
/// noisy Hamming weight of `SBox(pt[t][i] ^ key[i])`.
pub fn collect_traces(
    key: &[u8; 16],
    n_traces: usize,
    sigma: f64,
    seed: u64,
) -> (Vec<[u8; 16]>, Vec<[f64; 16]>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut pts = Vec::with_capacity(n_traces);
    let mut traces = Vec::with_capacity(n_traces);
    for _ in 0..n_traces {
        let mut pt = [0u8; 16];
        rng.fill(&mut pt);
        let mut trace = [0f64; 16];
        for i in 0..16 {
            let sbox_out = AES_SBOX[(pt[i] ^ key[i]) as usize];
            trace[i] = hamming_weight(sbox_out) as f64 + gauss(&mut rng, sigma);
        }
        pts.push(pt);
        traces.push(trace);
    }
    (pts, traces)
}

/// Pearson correlation between two equal-length sequences.
fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mean_x: f64 = x.iter().sum::<f64>() / n;
    let mean_y: f64 = y.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut sx = 0.0;
    let mut sy = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        num += dx * dy;
        sx += dx * dx;
        sy += dy * dy;
    }
    if sx == 0.0 || sy == 0.0 {
        0.0
    } else {
        num / (sx * sy).sqrt()
    }
}

/// Recover one AES key byte at position `byte_idx` from
/// `(plaintexts, traces)` via Pearson CPA.
pub fn recover_key_byte(
    pts: &[[u8; 16]],
    traces: &[[f64; 16]],
    byte_idx: usize,
) -> (u8, f64) {
    let n = traces.len();
    let leakage: Vec<f64> = (0..n).map(|t| traces[t][byte_idx]).collect();
    let mut best = (0u8, 0.0f64);
    for k_hat in 0..=255u8 {
        let predictions: Vec<f64> = (0..n)
            .map(|t| {
                let s = AES_SBOX[(pts[t][byte_idx] ^ k_hat) as usize];
                hamming_weight(s) as f64
            })
            .collect();
        let corr = pearson(&predictions, &leakage).abs();
        if corr > best.1 {
            best = (k_hat, corr);
        }
    }
    best
}

/// Recover the full 16-byte AES key.
pub fn cpa_recover(pts: &[[u8; 16]], traces: &[[f64; 16]]) -> ([u8; 16], [f64; 16]) {
    let mut key = [0u8; 16];
    let mut corrs = [0f64; 16];
    for i in 0..16 {
        let (k, c) = recover_key_byte(pts, traces, i);
        key[i] = k;
        corrs[i] = c;
    }
    (key, corrs)
}

pub fn run() -> Report {
    let mut r = Report::new(74, "Correlation Power Analysis on AES-128");
    let key: [u8; 16] = [
        0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f,
        0x3c,
    ];
    let n_traces = 4000;
    let sigma = 1.5;
    r.line(format!("True key      : {}", hex::encode(key)));
    r.line(format!("traces        : {} (σ = {})", n_traces, sigma));
    let (pts, traces) = collect_traces(&key, n_traces, sigma, 0xC0FFEE);
    let (recovered, corrs) = cpa_recover(&pts, &traces);
    r.line(format!("Recovered key : {}", hex::encode(recovered)));
    let mean_corr = corrs.iter().sum::<f64>() / 16.0;
    r.line(format!("Mean |corr|   : {:.3}", mean_corr));
    let ok = recovered == key;
    r.line(format!("Match         : {}", ok));
    if ok {
        return r.succeed();
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpa_recovers_aes_key() {
        let key: [u8; 16] = [
            0xde, 0xad, 0xbe, 0xef, 0xca, 0xfe, 0xba, 0xbe, 0xfe, 0xed, 0xfa, 0xce, 0x12, 0x34,
            0x56, 0x78,
        ];
        let (pts, traces) = collect_traces(&key, 4000, 1.5, 42);
        let (recovered, _) = cpa_recover(&pts, &traces);
        assert_eq!(recovered, key);
    }

    #[test]
    fn pearson_basic_cases() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((pearson(&x, &y) - 1.0).abs() < 1e-9);
        let z = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        assert!((pearson(&x, &z) + 1.0).abs() < 1e-9);
    }
}
