//! The statistics a comparison reports, deterministic in their seed.
//!
//! Rho, BSGS on a random target, and the kangaroo are Las Vegas
//! algorithms: their cost is a random variable over the target and the
//! walk's seed, with a spread of the same order as its mean.  A single
//! run says little; a ratio of two single runs says less.  Every ratio
//! here therefore comes with a bootstrap interval, resampled within
//! workloads (stratified) and over the pairs a round produced, so the
//! interval reflects the variance the design actually has.

/// Arithmetic mean; `None` when empty.
pub fn mean(xs: &[f64]) -> Option<f64> {
    (!xs.is_empty()).then(|| xs.iter().sum::<f64>() / xs.len() as f64)
}

/// The `q`-quantile (`0 ≤ q ≤ 1`) by linear interpolation of the sorted
/// sample; `None` when empty.
pub fn quantile(xs: &[f64], q: f64) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }
    let mut v = xs.to_vec();
    v.sort_by(f64::total_cmp);
    let pos = q.clamp(0.0, 1.0) * (v.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    Some(v[lo] + (v[hi] - v[lo]) * (pos - lo as f64))
}

pub fn median(xs: &[f64]) -> Option<f64> {
    quantile(xs, 0.5)
}

/// Sample standard deviation; `None` below two points.
pub fn stdev(xs: &[f64]) -> Option<f64> {
    if xs.len() < 2 {
        return None;
    }
    let m = mean(xs)?;
    Some((xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (xs.len() - 1) as f64).sqrt())
}

/// Geometric mean of positive values; `None` when empty or any is not
/// positive.
pub fn geomean(xs: &[f64]) -> Option<f64> {
    if xs.is_empty() || xs.iter().any(|x| *x <= 0.0 || !x.is_finite()) {
        return None;
    }
    Some((xs.iter().map(|x| x.ln()).sum::<f64>() / xs.len() as f64).exp())
}

/// A small deterministic generator for resampling.
pub struct Resampler(u64);

impl Resampler {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }
    pub fn below(&mut self, n: usize) -> usize {
        self.0 = crate::cryptanalysis::ecbench::canonical::splitmix64(self.0);
        (self.0 % n as u64) as usize
    }
}

/// Percentile bootstrap of `stat` over strata: each resample draws, in
/// every stratum, as many items as it has, with replacement, and
/// evaluates `stat` on the resampled strata.  Returns the `(2.5 %, 97.5
/// %)` interval, or `None` when a stratum is empty or the statistic is
/// undefined on too many resamples.
///
/// This treats the strata as fixed: it answers "how sure are we about
/// *these* workloads".  [`cluster_bootstrap_ci`] also resamples the
/// strata, which is the question a claim about a curve asks.
pub fn bootstrap_ci<T: Clone>(
    strata: &[Vec<T>],
    resamples: usize,
    seed: u64,
    stat: impl Fn(&[Vec<T>]) -> Option<f64>,
) -> Option<(f64, f64)> {
    resample(strata, resamples, seed, false, stat)
}

/// Two-stage (cluster) bootstrap: each resample first draws as many
/// strata as there are, with replacement, then items within each drawn
/// stratum.  With workloads as strata, the interval covers both the
/// variation between targets (BSGS's cost is fixed by its target, so all
/// its variance is here) and the variation between seeds within a
/// target (rho's).  Undefined below two strata, where it would be the
/// one-stage interval in disguise.
pub fn cluster_bootstrap_ci<T: Clone>(
    strata: &[Vec<T>],
    resamples: usize,
    seed: u64,
    stat: impl Fn(&[Vec<T>]) -> Option<f64>,
) -> Option<(f64, f64)> {
    if strata.len() < 2 {
        return None;
    }
    resample(strata, resamples, seed, true, stat)
}

fn resample<T: Clone>(
    strata: &[Vec<T>],
    resamples: usize,
    seed: u64,
    clusters: bool,
    stat: impl Fn(&[Vec<T>]) -> Option<f64>,
) -> Option<(f64, f64)> {
    if strata.is_empty() || strata.iter().any(|s| s.is_empty()) {
        return None;
    }
    let mut rng = Resampler::new(seed);
    let mut values = Vec::with_capacity(resamples);
    for _ in 0..resamples {
        let picked: Vec<&Vec<T>> = if clusters {
            (0..strata.len())
                .map(|_| &strata[rng.below(strata.len())])
                .collect()
        } else {
            strata.iter().collect()
        };
        let draw: Vec<Vec<T>> = picked
            .iter()
            .map(|s| {
                (0..s.len())
                    .map(|_| s[rng.below(s.len())].clone())
                    .collect()
            })
            .collect();
        if let Some(v) = stat(&draw) {
            if v.is_finite() {
                values.push(v);
            }
        }
    }
    if values.len() < resamples * 9 / 10 {
        return None;
    }
    Some((quantile(&values, 0.025)?, quantile(&values, 0.975)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantiles_interpolate() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), Some(2.0));
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), Some(2.5));
        assert_eq!(quantile(&[], 0.5), None);
        assert!((geomean(&[1.0, 4.0]).unwrap() - 2.0).abs() < 1e-12);
        assert_eq!(geomean(&[1.0, 0.0]), None);
    }

    #[test]
    fn bootstrap_is_deterministic_and_brackets_the_mean() {
        let strata = vec![vec![1.0, 2.0, 3.0, 4.0, 5.0], vec![10.0, 11.0, 12.0]];
        let f = |s: &[Vec<f64>]| mean(&s.concat());
        let a = bootstrap_ci(&strata, 2000, 7, f).unwrap();
        let b = bootstrap_ci(&strata, 2000, 7, f).unwrap();
        assert_eq!(a, b);
        let m = mean(&strata.concat()).unwrap();
        assert!(a.0 < m && m < a.1);
    }

    #[test]
    fn cluster_bootstrap_sees_between_stratum_spread() {
        // No spread inside strata, all of it between them: the one-stage
        // interval collapses to a point, the two-stage one does not.
        let strata = vec![vec![1.0; 4], vec![2.0; 4], vec![3.0; 4]];
        let f = |s: &[Vec<f64>]| mean(&s.concat());
        let one = bootstrap_ci(&strata, 2000, 1, f).unwrap();
        let two = cluster_bootstrap_ci(&strata, 2000, 1, f).unwrap();
        assert_eq!(one.0, one.1);
        assert!(two.1 - two.0 > 0.5);
        assert!(cluster_bootstrap_ci(&strata[..1], 100, 1, f).is_none());
    }
}
