//! Conservative Bernoulli Chernoff bounds, with finite-query PRP coupling slack.
//! No normal approximation or chi-square tail is used for sparse events.
pub fn null_probability(rank: u32, zero_target: bool) -> f64 {
    let excluded_zero = 2.0_f64.powi(-128);
    (2.0_f64.powi(-(rank as i32)) - if zero_target { excluded_zero } else { 0.0 })
        / (1.0 - excluded_zero)
}

pub fn kl(x: f64, p: f64) -> f64 {
    if x == p {
        return 0.0;
    }
    if p == 0.0 || p == 1.0 {
        return f64::INFINITY;
    }
    let a = if x == 0.0 { 0.0 } else { x * (x / p).ln() };
    let b = if x == 1.0 {
        0.0
    } else {
        (1.0 - x) * ((1.0 - x) / (1.0 - p)).ln()
    };
    (a + b).max(0.0)
}

/// Two-sided union of Chernoff tails. It is an upper bound, not an exact p-value.
pub fn tail_bound(hits: u64, samples: u64, p: f64) -> f64 {
    if samples == 0 {
        return 1.0;
    }
    (2.0 * (-(samples as f64) * kl(hits as f64 / samples as f64, p)).exp()).min(1.0)
}

/// IID independent output pairs coupled to a lazy uniform permutation on at
/// most q distinct inputs. Union bound on cross-pair output collisions; the
/// denominator is conservatively reduced for within-pair exclusions.
pub fn coupling_slack(q: u64) -> f64 {
    let q = q as f64;
    (q * (q - 1.0).max(0.0) / (2.0 * (2.0_f64.powi(128) - 1.0))).min(1.0)
}

pub fn corrected_bound(hits: u64, n: u64, p: f64, tests: usize, q: u64) -> f64 {
    ((tail_bound(hits, n, p) + coupling_slack(q)) * tests as f64).min(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_and_boundary_cases() {
        assert_eq!(tail_bound(0, 1, 0.0), 1.0);
        assert_eq!(tail_bound(1, 1, 0.0), 0.0);
        assert_eq!(tail_bound(0, 0, 0.5), 1.0);
        assert_eq!(tail_bound(500, 1000, 0.5), 1.0);
        assert!(tail_bound(0, 8192, 1.0 / 256.0) < 1e-12);
        assert_eq!(null_probability(128, true), 0.0);
        assert!(coupling_slack(1_000_000) < 1e-25);
        assert!(
            corrected_bound(0, 8192, 1.0 / 256.0, 100, 16384) >= tail_bound(0, 8192, 1.0 / 256.0)
        );
    }
}
