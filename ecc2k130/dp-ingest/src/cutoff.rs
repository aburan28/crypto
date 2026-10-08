//! The campaign's distinguished-point weight and what a slot's ratio says.

pub const CAMPAIGN_DP_WEIGHT: u32 = 32;
pub const CAMPAIGN_ITER_PER_DP_LOG2: f64 = 28.41;
pub const DP_RATIO_TOLERANCE_LOG2: f64 = 1.0;
pub const DP_RATIO_MIN_RECORDS: u64 = 50_000;
pub const FIELD_BITS: u32 = 131;

/// `log2(C(n, k))` via summing log differences.
fn log2_binomial(n: u32, k: u32) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    let k = k.min(n - k);
    let mut acc = 0.0;
    for i in 0..k {
        acc += f64::from(n - i).log2() - f64::from(i + 1).log2();
    }
    acc
}

/// Iterations per point for `HW(x) ≤ weight` over uniform `m`-bit strings.
pub fn theoretical_iter_per_dp_log2(weight: u32, m: u32) -> f64 {
    let weight = weight.min(m);
    let logs: Vec<f64> = (0..=weight).map(|k| log2_binomial(m, k)).collect();
    let max = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = logs.iter().map(|l| 2f64.powf(l - max)).sum();
    f64::from(m) - (max + sum.log2())
}

/// The cutoff that gives about this interval.
///
/// The curve's x-coordinates run ~0.6 shorter in the exponent than uniform
/// strings, so the binomial tail is shifted by the live-fleet measurement
/// before the nearest weight is read.
pub fn estimate_dp_weight(log2_iter_per_dp: f64, m: u32) -> u32 {
    let shift = CAMPAIGN_ITER_PER_DP_LOG2 - theoretical_iter_per_dp_log2(CAMPAIGN_DP_WEIGHT, m);
    let mut best = 1u32;
    let mut best_gap = f64::INFINITY;
    for k in 1..m {
        let gap = (theoretical_iter_per_dp_log2(k, m) + shift - log2_iter_per_dp).abs();
        if gap < best_gap {
            best = k;
            best_gap = gap;
        }
    }
    best
}

/// `(log2 iterations per point, estimated weight, at campaign weight?)`.
///
/// The last is `Some(true/false)`, or `None` when there is too little to judge.
pub fn dp_weight_verdict(
    iterations: u128,
    records: u64,
) -> (Option<f64>, Option<u32>, Option<bool>) {
    if records < DP_RATIO_MIN_RECORDS || iterations == 0 {
        return (None, None, None);
    }
    let ratio = (iterations as f64 / records as f64).log2();
    let ok = (ratio - CAMPAIGN_ITER_PER_DP_LOG2).abs() <= DP_RATIO_TOLERANCE_LOG2;
    let rounded = (ratio * 1000.0).round() / 1000.0;
    (
        Some(rounded),
        Some(estimate_dp_weight(ratio, FIELD_BITS)),
        Some(ok),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_theoretical_interval_is_the_binomial_tail() {
        assert!((theoretical_iter_per_dp_log2(32, FIELD_BITS) - 29.01).abs() < 0.005);
        assert!((theoretical_iter_per_dp_log2(34, FIELD_BITS) - 25.84).abs() < 0.005);
    }

    #[test]
    fn the_estimate_reads_the_live_ratios() {
        assert_eq!(estimate_dp_weight(28.41, FIELD_BITS), 32);
        assert_eq!(estimate_dp_weight(25.27, FIELD_BITS), 34);
        assert_eq!(estimate_dp_weight(25.30, FIELD_BITS), 34);
        assert_eq!(estimate_dp_weight(23.58, FIELD_BITS), 35);
    }

    #[test]
    fn a_campaign_slot_passes_and_a_loose_one_does_not() {
        let (ratio, weight, ok) = dp_weight_verdict(76_797_696u128 * 6_160_384, 1_324_201);
        assert_eq!(ok, Some(true));
        assert_eq!(weight, Some(32));
        assert!((ratio.unwrap() - 28.41).abs() < 0.05);

        let (ratio, weight, ok) = dp_weight_verdict(8_533_504u128 * 4_000_000, 1_699_148);
        assert_eq!(ok, Some(false));
        assert_eq!(weight, Some(35));
        let _ = ratio;

        let (_, weight, ok) = dp_weight_verdict(19_576_320u128 * 4_000_000, 1_900_000);
        assert_eq!(ok, Some(false));
        assert_eq!(weight, Some(34));
    }

    #[test]
    fn too_few_records_is_no_verdict() {
        assert_eq!(dp_weight_verdict(10u128.pow(15), 100), (None, None, None));
        assert_eq!(dp_weight_verdict(0, 1_000_000), (None, None, None));
    }
}
