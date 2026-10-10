//! EXP7: complete Fourier sweep of cheap point predicates on toy prime curves.
//! Protocols: research/prime_fourier_flatness_exp7_20261007/PROTOCOL.md
//! and research/prime_fourier_flatness_exp7_20261007/PROTOCOL_ALL_PATTERNS.md

use crypto_lib::cryptanalysis::ic_boundary::{
    find_prime_order_curve, CountedGroup, GroupOps, PrimeCurve, PrimeInstance, PrimePoint,
};
use crypto_lib::hash::sha256::sha256;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::f64::consts::PI;
use std::fs;
use std::ops::{Add, Mul, Sub};
use std::path::Path;

const CURVE_BITS: [u32; 4] = [8, 10, 12, 14];
const CURVE_SEED: u64 = 20261007;
const NULL_REPS: usize = 2048;
const CASES: [&str; 16] = [
    "smooth",
    "cf16",
    "cantor3",
    "hamming",
    "farey",
    "legendre+++",
    "legendre++-",
    "legendre+-+",
    "legendre+--",
    "legendre-++",
    "legendre-+-",
    "legendre--+",
    "legendre---",
    "sha-x",
    "random-pairs",
    "log-interval",
];

#[derive(Clone, Copy, Debug, Default)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    fn scale(self, c: f64) -> Self {
        Self::new(self.re * c, self.im * c)
    }
}

impl Add for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl Sub for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl Mul for Complex {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

fn fft_power_of_two(values: &mut [Complex], inverse: bool) {
    let n = values.len();
    assert!(n.is_power_of_two());
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            values.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let sign = if inverse { 1.0 } else { -1.0 };
        let (sin, cos) = (sign * 2.0 * PI / len as f64).sin_cos();
        let root = Complex::new(cos, sin);
        for chunk in values.chunks_exact_mut(len) {
            let mut w = Complex::new(1.0, 0.0);
            for k in 0..len / 2 {
                let u = chunk[k];
                let v = chunk[k + len / 2] * w;
                chunk[k] = u + v;
                chunk[k + len / 2] = u - v;
                w = w * root;
            }
        }
        len *= 2;
    }
    if inverse {
        for value in values {
            *value = value.scale(1.0 / n as f64);
        }
    }
}

/// Arbitrary-length DFT, using Bluestein's chirp convolution.
fn dft(values: &[Complex], inverse: bool) -> Vec<Complex> {
    let n = values.len();
    if n == 0 {
        return Vec::new();
    }
    if inverse {
        let conjugated: Vec<_> = values.iter().map(|v| v.conj()).collect();
        return dft(&conjugated, false)
            .into_iter()
            .map(|v| v.conj().scale(1.0 / n as f64))
            .collect();
    }
    let size = (2 * n - 1).next_power_of_two();
    let mut a = vec![Complex::default(); size];
    let mut b = vec![Complex::default(); size];
    for k in 0..n {
        // k^2 modulo 2n gives the same complex phase with smaller arguments.
        let phase = PI * ((k * k) % (2 * n)) as f64 / n as f64;
        let (sin, cos) = phase.sin_cos();
        let chirp = Complex::new(cos, -sin);
        a[k] = values[k] * chirp;
        b[k] = chirp.conj();
        if k > 0 {
            b[size - k] = b[k];
        }
    }
    fft_power_of_two(&mut a, false);
    fft_power_of_two(&mut b, false);
    for k in 0..size {
        a[k] = a[k] * b[k];
    }
    fft_power_of_two(&mut a, true);
    (0..n)
        .map(|k| {
            let phase = PI * ((k * k) % (2 * n)) as f64 / n as f64;
            let (sin, cos) = phase.sin_cos();
            a[k] * Complex::new(cos, -sin)
        })
        .collect()
}

fn spectrum(members: &[bool]) -> Vec<Complex> {
    dft(
        &members
            .iter()
            .map(|&member| Complex::new(f64::from(member), 0.0))
            .collect::<Vec<_>>(),
        false,
    )
}

fn peak(fft: &[Complex], members: usize) -> Option<(usize, f64, f64)> {
    let n = fft.len();
    if members == 0 || members == n {
        return None;
    }
    let (j, value) = (1..=n / 2)
        .map(|j| (j, fft[j].abs()))
        .max_by(|a, b| a.1.total_cmp(&b.1))?;
    let norm = (members as f64 * (1.0 - members as f64 / n as f64)).sqrt();
    Some((j, value, value / norm))
}

fn pair_sum_max_deviation(fft: &[Complex], members: usize) -> (f64, f64) {
    let square: Vec<_> = fft.iter().map(|v| *v * *v).collect();
    let convolution = dft(&square, true);
    let expected = (members * members) as f64 / fft.len() as f64;
    let mut deviation = 0.0_f64;
    let mut numeric_error = 0.0_f64;
    for value in convolution {
        deviation = deviation.max((value.re - expected).abs());
        numeric_error = numeric_error.max((value.re - value.re.round()).abs().max(value.im.abs()));
    }
    (deviation, numeric_error)
}

fn isqrt(n: u64) -> u64 {
    let mut r = (n as f64).sqrt() as u64;
    while (r + 1) * (r + 1) <= n {
        r += 1;
    }
    while r * r > n {
        r -= 1;
    }
    r
}

fn is_smooth(mut x: u64, bound: u64) -> bool {
    if x == 0 {
        return false;
    }
    for factor in 2..=bound {
        while x.is_multiple_of(factor) {
            x /= factor;
        }
    }
    x == 1
}

fn cf_bounded(x: u64, p: u64, bound: u64) -> bool {
    if x == 0 {
        return false;
    }
    let (mut a, mut b) = (p, x);
    while b > 0 {
        if a / b > bound {
            return false;
        }
        (a, b) = (b, a % b);
    }
    true
}

fn cantor_digits(mut x: u64) -> bool {
    while x > 0 {
        if x % 3 == 1 {
            return false;
        }
        x /= 3;
    }
    true
}

fn farey_height(x: u64, p: u64, height: u64) -> bool {
    (1..=height).any(|b| {
        let residue = (x * b) % p;
        residue <= height || p - residue <= height
    })
}

fn coordinate_member(case: &str, x: u64, curve: &PrimeCurve) -> bool {
    let p = curve.p;
    match case {
        "smooth" => is_smooth(x, isqrt(isqrt(p)).max(2)),
        "cf16" => cf_bounded(x, p, 16),
        "cantor3" => cantor_digits(x),
        "hamming" => x.count_ones() <= (63 - p.leading_zeros()) / 2,
        "farey" => farey_height(x, p, (isqrt(p) / 3).max(2)),
        name if name.starts_with("legendre") => {
            let signs = name.trim_start_matches("legendre").as_bytes();
            assert_eq!(signs.len(), 3);
            signs.iter().enumerate().all(|(i, &sign)| {
                let desired = match sign {
                    b'+' => 1,
                    b'-' => -1,
                    _ => panic!("bad Legendre sign"),
                };
                curve.legendre((x + i as u64) % p) == desired
            })
        }
        "sha-x" => {
            let mut input = b"exp7-sha-x-v1".to_vec();
            input.extend_from_slice(&p.to_le_bytes());
            input.extend_from_slice(&x.to_le_bytes());
            sha256(&input)[0] & 0xc0 == 0
        }
        _ => panic!("unknown coordinate case {case}"),
    }
}

fn null_seed(p: u64, case: &str, replicate: usize) -> [u8; 32] {
    let mut input = b"exp7-null-v1".to_vec();
    input.extend_from_slice(&p.to_le_bytes());
    input.extend_from_slice(case.as_bytes());
    input.extend_from_slice(&(replicate as u64).to_le_bytes());
    sha256(&input)
}

fn random_set(n: usize, m: usize, paired: bool, seed: [u8; 32]) -> Vec<bool> {
    let mut rng = StdRng::from_seed(seed);
    let mut selected = vec![false; n];
    if paired {
        assert!(n % 2 == 1 && m.is_multiple_of(2) && m < n);
        let mut pair_indices: Vec<_> = (1..=n / 2).collect();
        pair_indices.shuffle(&mut rng);
        for &k in pair_indices.iter().take(m / 2) {
            selected[k] = true;
            selected[n - k] = true;
        }
    } else {
        let mut indices: Vec<_> = (1..n).collect();
        indices.shuffle(&mut rng);
        for &k in indices.iter().take(m) {
            selected[k] = true;
        }
    }
    selected
}

fn enumerate_group(instance: &PrimeInstance) -> (Vec<PrimePoint>, GroupOps) {
    assert_eq!(instance.cofactor, 1);
    assert_eq!(instance.group_order, instance.r);
    let n = instance.r as usize;
    let generator = instance.generator_point();
    assert!(instance.curve.is_on_curve(generator));
    let mut points = Vec::with_capacity(n);
    let mut seen = HashSet::with_capacity(n);
    let mut point = PrimePoint::INFINITY;
    let mut ops = GroupOps::default();
    for k in 0..n {
        assert!(instance.curve.is_on_curve(point), "off curve at log {k}");
        assert!(
            seen.insert((point.infinity, point.x, point.y)),
            "duplicate at log {k}"
        );
        points.push(point);
        point = instance.curve.add(&mut ops, point, generator);
    }
    assert_eq!(point, PrimePoint::INFINITY, "[n]G must be identity");
    assert_eq!(seen.len(), n);
    (points, ops)
}

fn bitset_hex(bits: &[bool]) -> String {
    let mut bytes = vec![0u8; bits.len().div_ceil(8)];
    for (k, &member) in bits.iter().enumerate() {
        if member {
            bytes[k / 8] |= 1 << (k % 8);
        }
    }
    hex::encode(bytes)
}

fn quantile(sorted: &[f64], fraction: f64) -> f64 {
    let index = ((fraction * sorted.len() as f64).ceil() as usize)
        .saturating_sub(1)
        .min(sorted.len() - 1);
    sorted[index]
}

fn case_result(case: &str, bits: &[bool], p: u64, reps: usize) -> Value {
    let n = bits.len();
    let members = bits.iter().filter(|&&member| member).count();
    let transform = spectrum(bits);
    let (pair_dev, pair_numeric_error) = pair_sum_max_deviation(&transform, members);
    assert!(
        pair_numeric_error < 1e-6,
        "pair convolution lost integer precision"
    );
    let energy: f64 = transform.iter().map(|v| v.abs().powi(2)).sum();
    let parseval_relerr = (energy - (n * members) as f64).abs() / (n * members).max(1) as f64;
    assert!(parseval_relerr < 1e-8, "Parseval check failed");
    let Some((frequency, amplitude, z)) = peak(&transform, members) else {
        return json!({
            "name": case, "status": "degenerate", "members": members,
            "members_bitset_hex_lsb": bitset_hex(bits), "pair_sum_max_abs_deviation": pair_dev,
            "pair_numeric_error": pair_numeric_error, "parseval_relerr": parseval_relerr
        });
    };
    let mut top: Vec<_> = (1..=n / 2).map(|j| (j, transform[j])).collect();
    top.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
    let top_eight: Vec<_> = top
        .iter()
        .take(8)
        .map(|(j, v)| json!({"frequency": j, "real": v.re, "imag": v.im, "absolute": v.abs()}))
        .collect();
    let paired = case != "log-interval";
    let mut null_peaks = Vec::with_capacity(reps);
    for replicate in 0..reps {
        let null_bits = random_set(n, members, paired, null_seed(p, case, replicate));
        let null_fft = spectrum(&null_bits);
        null_peaks.push(peak(&null_fft, members).unwrap().2);
    }
    let exceedances = null_peaks.iter().filter(|&&v| v >= z).count();
    let empirical_p = (exceedances + 1) as f64 / (reps + 1) as f64;
    let mut sorted = null_peaks.clone();
    sorted.sort_by(f64::total_cmp);
    let median = (sorted[(reps - 1) / 2] + sorted[reps / 2]) / 2.0;
    json!({
        "name": case, "status": "measured", "members": members,
        "density": members as f64 / n as f64,
        "members_bitset_hex_lsb": bitset_hex(bits),
        "peak_frequency": frequency, "peak_absolute": amplitude, "peak_z": z,
        "top_eight": top_eight,
        "pair_sum_max_abs_deviation": pair_dev,
        "pair_numeric_error": pair_numeric_error,
        "parseval_relerr": parseval_relerr,
        "null_paired": paired, "null_replicates": reps,
        "null_peaks_z": null_peaks,
        "null_median_z": median, "null_q95_z": quantile(&sorted, 0.95),
        "null_max_z": sorted[reps - 1], "null_exceedances": exceedances,
        "empirical_p": empirical_p
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let output = args
        .next()
        .ok_or("usage: exp7_fourier_flatness RESULTS.json")?;
    if args.next().is_some() {
        return Err("usage: exp7_fourier_flatness RESULTS.json".into());
    }
    let source_revision = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()?;
    if !source_revision.status.success() {
        return Err("git rev-parse HEAD failed".into());
    }
    let source_revision = String::from_utf8(source_revision.stdout)?.trim().to_owned();
    let mut curves = Vec::new();
    for bits in CURVE_BITS {
        let instance = find_prime_order_curve(bits, CURVE_SEED);
        let (points, enumeration_ops) = enumerate_group(&instance);
        let n = points.len();
        let mut rows = Vec::new();
        for case in CASES {
            let members: Vec<bool> = match case {
                "random-pairs" => {
                    let count = 2 * ((n - 1) / 8);
                    random_set(n, count, true, null_seed(instance.curve.p, case, NULL_REPS))
                }
                "log-interval" => (0..n).map(|k| (1..=(n - 1) / 4).contains(&k)).collect(),
                _ => points
                    .iter()
                    .map(|pt| !pt.infinity && coordinate_member(case, pt.x, &instance.curve))
                    .collect(),
            };
            let row = case_result(case, &members, instance.curve.p, NULL_REPS);
            eprintln!(
                "bits={bits} p={} n={} {:<14} m={:<6} z={:.4} null_p={:.5}",
                instance.curve.p,
                n,
                case,
                row["members"].as_u64().unwrap_or(0),
                row["peak_z"].as_f64().unwrap_or(f64::NAN),
                row["empirical_p"].as_f64().unwrap_or(f64::NAN)
            );
            rows.push(row);
        }
        let point_rows: Vec<_> = points
            .iter()
            .map(|pt| {
                if pt.infinity {
                    Value::Null
                } else {
                    json!([pt.x, pt.y])
                }
            })
            .collect();
        curves.push(json!({
            "requested_bits": bits, "icv1": instance.curve_id().icv1,
            "slug": instance.curve_id().slug,
            "p": instance.curve.p, "a": instance.curve.a, "b": instance.curve.b,
            "n": instance.r, "cofactor": instance.cofactor,
            "generator": [instance.generator.0, instance.generator.1],
            "enumeration_group_ops": enumeration_ops,
            "points_by_log": point_rows, "cases": rows
        }));
    }
    let threshold = 0.05 / 52.0;
    let mut leads = Vec::new();
    let mut controls_pass = true;
    for curve in &curves {
        for case in curve["cases"].as_array().unwrap() {
            let name = case["name"].as_str().unwrap();
            if name == "log-interval" {
                controls_pass &= case["peak_z"].as_f64().unwrap_or(0.0)
                    > case["null_q95_z"].as_f64().unwrap_or(f64::INFINITY);
            } else if name != "sha-x" && name != "random-pairs" {
                if let Some(p) = case["empirical_p"].as_f64() {
                    if p <= threshold {
                        leads.push(json!({
                            "curve": curve["icv1"], "case": name,
                            "empirical_p": p, "peak_z": case["peak_z"]
                        }));
                    }
                }
            }
        }
    }
    let document = json!({
        "schema": "exp7-fourier-flatness-v2", "source_revision": source_revision,
        "curve_seed": CURVE_SEED, "null_replicates": NULL_REPS,
        "candidate_family_cells": 52, "bonferroni_threshold": threshold,
        "positive_controls_pass": controls_pass, "leads": leads,
        "bitset_order": "byte k/8, bit k%8; index 0 is identity",
        "null_quantile": "nearest-rank ceil(q*R)-1; median averages middle two",
        "curves": curves
    });
    let bytes = serde_json::to_vec_pretty(&document)?;
    let path = Path::new(&output);
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, bytes)?;
    println!(
        "wrote {output}: {} curves, {} adjusted leads, controls_pass={controls_pass}",
        CURVE_BITS.len(),
        document["leads"].as_array().unwrap().len()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn direct_dft(input: &[Complex]) -> Vec<Complex> {
        (0..input.len())
            .map(|j| {
                input
                    .iter()
                    .enumerate()
                    .fold(Complex::default(), |acc, (k, value)| {
                        let phase = -2.0 * PI * (j * k) as f64 / input.len() as f64;
                        let (sin, cos) = phase.sin_cos();
                        acc + *value * Complex::new(cos, sin)
                    })
            })
            .collect()
    }

    #[test]
    fn bluestein_matches_direct_dft_and_inverse_for_nonpowers_of_two() {
        for n in 2..38 {
            let input: Vec<_> = (0..n)
                .map(|k| Complex::new((k * 7 % 13) as f64 / 13.0, (k * 5 % 11) as f64 / 11.0))
                .collect();
            let actual = dft(&input, false);
            let expected = direct_dft(&input);
            for (a, e) in actual.iter().zip(&expected) {
                assert!((*a - *e).abs() < 1e-8, "n={n}");
            }
            for (a, e) in dft(&actual, true).iter().zip(&input) {
                assert!((*a - *e).abs() < 1e-8, "inverse n={n}");
            }
        }
    }

    #[test]
    fn pair_sum_inverse_matches_direct_ordered_counts() {
        for n in 3..25 {
            let members: Vec<bool> = (0..n).map(|k| k % 3 == 0 || k % 5 == 1).collect();
            let transform = spectrum(&members);
            let squared: Vec<_> = transform.iter().map(|v| *v * *v).collect();
            let counts = dft(&squared, true);
            for t in 0..n {
                let direct = (0..n)
                    .filter(|&k| members[k] && members[(t + n - k) % n])
                    .count();
                assert!((counts[t].re - direct as f64).abs() < 1e-8);
                assert!(counts[t].im.abs() < 1e-8);
            }
        }
    }

    #[test]
    fn predicates_and_matched_nulls_obey_the_frozen_rules() {
        assert!(is_smooth(72, 3));
        assert!(!is_smooth(77, 3));
        assert!(cf_bounded(2, 17, 16));
        assert!(!cf_bounded(1, 17, 16));
        assert!(cantor_digits(20)); // 202 base 3
        assert!(!cantor_digits(21)); // 210 base 3
        assert!(farey_height(50, 101, 2)); // 2*x = -1 mod 101
        let set = random_set(101, 24, true, null_seed(101, "smooth", 0));
        assert!(!set[0]);
        assert_eq!(set.iter().filter(|&&v| v).count(), 24);
        for k in 1..101 {
            assert_eq!(set[k], set[101 - k]);
        }
        assert_eq!(set, random_set(101, 24, true, null_seed(101, "smooth", 0)));
        let curve = find_prime_order_curve(8, CURVE_SEED).curve;
        for x in 0..curve.p {
            let matching_patterns = CASES
                .iter()
                .filter(|name| name.starts_with("legendre"))
                .filter(|name| coordinate_member(name, x, &curve))
                .count();
            let expected = usize::from((0..3).all(|i| curve.legendre((x + i) % curve.p) != 0));
            assert_eq!(matching_patterns, expected, "x={x}");
        }
    }

    #[test]
    fn repository_prime_curve_enumerates_once_and_returns_to_identity() {
        let curve = find_prime_order_curve(8, CURVE_SEED);
        let (points, ops) = enumerate_group(&curve);
        assert_eq!(points.len(), curve.r as usize);
        assert_eq!(points[0], PrimePoint::INFINITY);
        // The counted API charges the O + G and final (n-1)G + G calls.
        assert_eq!(ops.adds + ops.doubles, curve.r);
    }
}
