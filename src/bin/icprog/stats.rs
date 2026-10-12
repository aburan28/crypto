//! The statistics the programme's rounds share (IC_TOOL_PROGRAM.md §5).
//!
//! The primary interval is ledger §22's: the geometric mean of paired
//! ratios, with a Student `t` interval on their logarithms.  The `t`
//! quantile comes from the regularised incomplete beta function, so no
//! round depends on a table.
//!
//! This is the native form of the retired `harness/stats.py`, with the
//! same arithmetic step for step: `fmean` is an exactly rounded sum
//! divided by `n` (Shewchuk's algorithm, as `math.fsum`), `stdev` is the
//! correctly rounded square root of the exact sample variance (as
//! Python 3.11's `statistics.stdev`), and `lgamma` is Lanczos' formula
//! with the coefficients `math.lgamma` uses.  The earlier rounds'
//! analyses are therefore reproduced to the last bit
//! (`icprog analyse r03`).

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::json::J;

// ── lgamma, as CPython's m_lgamma ───────────────────────────────────

// CPython's literals, digit for digit, so each rounds to the same double.
#[allow(clippy::excessive_precision)]
const LANCZOS_G: f64 = 6.024_680_040_776_729_583_740_234_375;
#[allow(clippy::excessive_precision)]
const LANCZOS_NUM: [f64; 13] = [
    23531376880.410759688572007674451636754734846804940,
    42919803642.649098768957899047001988850926355848959,
    35711959237.355668049440185451547166705960488635843,
    17921034426.037209699919755754458931112671403265390,
    6039542586.3520280050642916443072979210699388420708,
    1439720407.3117216736632230727949123939715485786772,
    248874557.86205415651146038641322942321632125127801,
    31426415.585400194380614231628318205362874684987640,
    2876370.6289353724412254090516208496135991145378768,
    186056.26539522349504029498971604569928220784236328,
    8071.6720023658162106380029022722506138218516325024,
    210.82427775157934587250973392071336271166969580291,
    2.5066282746310002701649081771338373386264310793408,
];
const LANCZOS_DEN: [f64; 13] = [
    0.0,
    39916800.0,
    120543840.0,
    150917976.0,
    105258076.0,
    45995730.0,
    13339535.0,
    2637558.0,
    357423.0,
    32670.0,
    1925.0,
    66.0,
    1.0,
];

fn lanczos_sum(x: f64) -> f64 {
    let (mut num, mut den) = (0.0f64, 0.0f64);
    if x < 5.0 {
        for i in (0..13).rev() {
            num = num * x + LANCZOS_NUM[i];
            den = den * x + LANCZOS_DEN[i];
        }
    } else {
        for i in 0..13 {
            num = num / x + LANCZOS_NUM[i];
            den = den / x + LANCZOS_DEN[i];
        }
    }
    num / den
}

/// `ln Γ(x)` for `x > 0`, the only arguments the beta function asks for.
pub fn lgamma(x: f64) -> f64 {
    assert!(x > 0.0 && x.is_finite(), "lgamma({x})");
    if x == x.floor() && x <= 2.0 {
        return 0.0;
    }
    if x < 1e-20 {
        return -x.ln();
    }
    let mut r = lanczos_sum(x).ln() - LANCZOS_G;
    r += (x - 0.5) * ((x + LANCZOS_G - 0.5).ln() - 1.0);
    r
}

// ── the incomplete beta and Student's t ─────────────────────────────

/// Continued fraction for the incomplete beta (Numerical Recipes, betacf).
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let (tiny, eps) = (1e-300, 3e-16);
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let guard = |v: f64| if v.abs() > tiny { v } else { tiny };
    let mut c = 1.0;
    let mut d = 1.0 / guard(1.0 - qab * x / qap);
    let mut h = d;
    for m in 1..400 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 / guard(1.0 + aa * d);
        c = guard(1.0 + aa / c);
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 / guard(1.0 + aa * d);
        c = guard(1.0 + aa / c);
        let delta = d * c;
        h *= delta;
        if (delta - 1.0).abs() < eps {
            break;
        }
    }
    h
}

pub fn betai(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front = (lgamma(a + b) - lgamma(a) - lgamma(b) + a * x.ln() + b * (-x).ln_1p()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * betacf(a, b, x) / a
    } else {
        1.0 - front * betacf(b, a, 1.0 - x) / b
    }
}

pub fn t_cdf(t: f64, dof: usize) -> f64 {
    let dof = dof as f64;
    let x = dof / (dof + t * t);
    let tail = 0.5 * betai(dof / 2.0, 0.5, x);
    if t >= 0.0 {
        1.0 - tail
    } else {
        tail
    }
}

/// The `p` quantile of Student's t with `dof` degrees of freedom, by bisection.
pub fn t_quantile(p: f64, dof: usize) -> f64 {
    let (mut lo, mut hi) = (-1e3f64, 1e3f64);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if t_cdf(mid, dof) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

// ── sums, means, medians ────────────────────────────────────────────

/// The exactly rounded sum (Shewchuk; CPython's `math.fsum`).
pub fn fsum(xs: &[f64]) -> f64 {
    let mut p: Vec<f64> = Vec::new();
    for &x0 in xs {
        assert!(x0.is_finite(), "fsum of {x0}");
        let mut x = x0;
        let mut i = 0;
        for j in 0..p.len() {
            let mut y = p[j];
            if x.abs() < y.abs() {
                std::mem::swap(&mut x, &mut y);
            }
            let hi = x + y;
            let yr = hi - x;
            let lo = y - yr;
            if lo != 0.0 {
                p[i] = lo;
                i += 1;
            }
            x = hi;
        }
        p.truncate(i);
        if x != 0.0 {
            assert!(x.is_finite(), "intermediate overflow in fsum");
            p.push(x);
        }
    }
    let mut n = p.len();
    let mut hi = 0.0;
    if n > 0 {
        n -= 1;
        hi = p[n];
        let mut lo = 0.0;
        while n > 0 {
            let x = hi;
            n -= 1;
            let y = p[n];
            hi = x + y;
            let yr = hi - x;
            lo = y - yr;
            if lo != 0.0 {
                break;
            }
        }
        if n > 0 && ((lo < 0.0 && p[n - 1] < 0.0) || (lo > 0.0 && p[n - 1] > 0.0)) {
            let y = lo * 2.0;
            let x = hi + y;
            let yr = x - hi;
            if y == yr {
                hi = x;
            }
        }
    }
    hi
}

pub fn fmean(xs: &[f64]) -> f64 {
    assert!(!xs.is_empty(), "fmean of nothing");
    fsum(xs) / xs.len() as f64
}

/// The median: the middle value, or the mean of the two middle values.
pub fn median(xs: &[f64]) -> f64 {
    assert!(!xs.is_empty(), "median of nothing");
    let mut v = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

// ── the exact sample deviation ──────────────────────────────────────

/// `x` as the exact fraction `num / 2^k`.
fn exact(x: f64) -> (BigInt, u32) {
    assert!(x.is_finite());
    if x == 0.0 {
        return (BigInt::zero(), 0);
    }
    let bits = x.to_bits();
    let sign = if bits >> 63 == 1 { -1i64 } else { 1 };
    let exp = ((bits >> 52) & 0x7ff) as i64;
    let frac = bits & ((1u64 << 52) - 1);
    let (mant, e) = if exp == 0 {
        (frac, -1074i64)
    } else {
        (frac | (1u64 << 52), exp - 1075)
    };
    let m = BigInt::from(sign) * BigInt::from(mant);
    if e >= 0 {
        (m << e as usize, 0)
    } else {
        (m, (-e) as u32)
    }
}

/// `n / d` correctly rounded to a float (Python's `int / int`).
fn div_to_f64(n: &BigUint, d: &BigUint) -> f64 {
    if n.is_zero() {
        return 0.0;
    }
    // Scale so that the quotient has 55 significant bits, then round
    // half to even with the remainder as the sticky bit.
    let shift = n.bits() as i64 - d.bits() as i64 - 55;
    let (nn, dd) = if shift >= 0 {
        (n.clone(), d << shift as usize)
    } else {
        (n << (-shift) as usize, d.clone())
    };
    let (mut q, r) = nn.div_rem(&dd);
    let mut e = shift;
    // q has 55 or 56 bits; bring it to exactly 53 with round-half-even.
    let extra = q.bits() as i64 - 53;
    let sticky = !r.is_zero();
    let low_mask = (BigUint::one() << extra as usize) - BigUint::one();
    let low = &q & &low_mask;
    q >>= extra as usize;
    e += extra;
    let half = BigUint::one() << (extra as usize - 1);
    if low > half || (low == half && (sticky || q.bit(0))) {
        q += BigUint::one();
        if q.bits() > 53 {
            q >>= 1;
            e += 1;
        }
    }
    let m = q.to_u64().expect("53 bits");
    (m as f64) * 2f64.powi(e as i32)
}

/// The integer square root of `n / m`, rounded to odd.
fn isqrt_frac_rto(n: &BigUint, m: &BigUint) -> BigUint {
    let a = (n / m).sqrt();
    if &(&a * &a) * m != *n {
        a | BigUint::one()
    } else {
        a
    }
}

/// `√(n / m)` correctly rounded (Python's `_float_sqrt_of_frac`).
fn float_sqrt_of_frac(n: &BigUint, m: &BigUint) -> f64 {
    const WIDTH: i64 = 2 * 53 + 3;
    let q = (n.bits() as i64 - m.bits() as i64 - WIDTH).div_euclid(2);
    let (num, den) = if q >= 0 {
        (
            isqrt_frac_rto(n, &(m << (2 * q) as usize)) << q as usize,
            BigUint::one(),
        )
    } else {
        (
            isqrt_frac_rto(&(n << (-2 * q) as usize), m),
            BigUint::one() << (-q) as usize,
        )
    };
    div_to_f64(&num, &den)
}

/// The sample standard deviation, as Python 3.11's `statistics.stdev`:
/// the exact variance of the exact data, then a correctly rounded root.
pub fn stdev(xs: &[f64]) -> f64 {
    let n = xs.len();
    assert!(n >= 2, "stdev needs two points");
    let parts: Vec<(BigInt, u32)> = xs.iter().map(|&x| exact(x)).collect();
    let k = parts.iter().map(|(_, k)| *k).max().unwrap_or(0);
    // Every x as X / 2^k; then ss = (n ΣX² − (ΣX)²) / (n 2^{2k}).
    let big: Vec<BigInt> = parts.iter().map(|(m, kk)| m << (k - kk) as usize).collect();
    let sx: BigInt = big.iter().sum();
    let sxx: BigInt = big.iter().map(|x| x * x).sum();
    let nn = BigInt::from(n);
    let ss_num = &nn * &sxx - &sx * &sx; // over n · 2^{2k}
    debug_assert!(!ss_num.is_negative());
    // mss = ss / (n − 1) = ss_num / (n (n−1) 2^{2k}), in lowest terms.
    let num = ss_num.to_biguint().expect("non-negative");
    let den = (BigUint::from(n) * BigUint::from(n - 1)) << (2 * k) as usize;
    let g = num.gcd(&den);
    if num.is_zero() {
        return 0.0;
    }
    float_sqrt_of_frac(&(&num / &g), &(&den / &g))
}

// ── the round's figures ─────────────────────────────────────────────

/// The geometric mean of `ratios` with a `t` interval on their logarithms.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GeoCi {
    pub n: usize,
    pub geomean: Option<f64>,
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub lo: Option<f64>,
    pub hi: Option<f64>,
}

pub fn geo_ci(ratios: &[f64], level: f64) -> GeoCi {
    let mut out = GeoCi {
        n: ratios.len(),
        ..GeoCi::default()
    };
    if ratios.is_empty() {
        return out;
    }
    let logs: Vec<f64> = ratios.iter().map(|x| x.ln()).collect();
    let m = fmean(&logs);
    out.geomean = Some(m.exp());
    out.min = ratios.iter().copied().reduce(f64::min);
    out.max = ratios.iter().copied().reduce(f64::max);
    if logs.len() > 1 {
        let h = t_quantile(0.5 + level / 2.0, logs.len() - 1) * stdev(&logs)
            / (logs.len() as f64).sqrt();
        out.lo = Some((m - h).exp());
        out.hi = Some((m + h).exp());
    }
    out
}

impl GeoCi {
    /// The record's keys in the order the rounds have always written them.
    pub fn fields(&self) -> Vec<(String, J)> {
        let mut kv = vec![("n".to_string(), J::Int(self.n as i128))];
        for (k, v) in [
            ("geomean", self.geomean),
            ("min", self.min),
            ("max", self.max),
            ("lo", self.lo),
            ("hi", self.hi),
        ] {
            if let Some(v) = v {
                kv.push((k.to_string(), J::Float(v)));
            }
        }
        kv
    }
}

/// Least squares of `ys` on `xs` with a `t` interval on the slope, as the
/// harness's `stats.fit`; `None` below three points.
pub fn fit(xs: &[f64], ys: &[f64], level: f64) -> Option<J> {
    let n = xs.len();
    if n < 3 || ys.len() != n {
        return None;
    }
    let (xm, ym) = (fmean(xs), fmean(ys));
    let sxx: f64 = xs.iter().map(|x| (x - xm).powi(2)).sum();
    let beta = xs
        .iter()
        .zip(ys)
        .map(|(x, y)| (x - xm) * (y - ym))
        .sum::<f64>()
        / sxx;
    let alpha = ym - beta * xm;
    let dof = n - 2;
    let s2 = xs
        .iter()
        .zip(ys)
        .map(|(x, y)| (y - alpha - beta * x).powi(2))
        .sum::<f64>()
        / dof as f64;
    let h = t_quantile(0.5 + level / 2.0, dof) * (s2 / sxx).sqrt();
    Some(J::Obj(vec![
        ("beta".into(), J::Float(beta)),
        ("lo".into(), J::Float(beta - h)),
        ("hi".into(), J::Float(beta + h)),
        ("alpha".into(), J::Float(alpha)),
        ("points".into(), J::Int(n as i128)),
        ("dof".into(), J::Int(dof as i128)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn t_quantiles_match_the_table() {
        // Two-sided 95%: t(0.975; dof), to the table's six decimals.
        for (dof, want) in [(1, 12.706205), (4, 2.776445), (9, 2.262157), (39, 2.022691)] {
            let t = t_quantile(0.975, dof);
            assert!((t - want).abs() < 1e-6, "dof {dof}: {t}");
        }
        // The closed forms at one and two degrees of freedom.
        let p: f64 = 0.975;
        let t1 = (std::f64::consts::PI * (p - 0.5)).tan();
        let t2 = (2.0 * p - 1.0) / (2.0 * p * (1.0 - p)).sqrt();
        assert!((t_quantile(p, 1) - t1).abs() < 1e-9 * t1);
        assert!((t_quantile(p, 2) - t2).abs() < 1e-12 * t2);
    }

    #[test]
    fn lgamma_agrees_with_the_factorials() {
        for (x, fact) in [(3.0, 2.0f64), (5.0, 24.0), (10.0, 362_880.0)] {
            assert!((lgamma(x) - fact.ln()).abs() < 1e-13, "{x}");
        }
        assert!((lgamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-14);
        assert_eq!(lgamma(1.0), 0.0);
        assert_eq!(lgamma(2.0), 0.0);
    }

    #[test]
    fn fsum_is_exactly_rounded() {
        assert_eq!(fsum(&[0.1; 10]), 1.0);
        assert_eq!(
            fsum(&[1e100, 1.0, -1e100, 1e-100, 1e50, -1.0, -1e50]),
            1e-100
        );
        assert_eq!(fsum(&[1e-16, 1.0, 1e16]), 10000000000000002.0);
    }

    #[test]
    fn stdev_is_the_correctly_rounded_root_of_the_exact_variance() {
        // Python: statistics.stdev([1.5, 2.5, 2.5, 2.75, 3.25, 4.75])
        assert_eq!(
            stdev(&[1.5, 2.5, 2.5, 2.75, 3.25, 4.75]),
            1.0810874155219827
        );
        assert_eq!(stdev(&[2.0, 2.0]), 0.0);
        // The exact variance of [0, 1] is 1/2, whose root IEEE sqrt rounds correctly.
        assert_eq!(stdev(&[0.0, 1.0]), 0.5f64.sqrt());
        assert_eq!(stdev(&[1.0, 3.0]), 2f64.sqrt());
    }

    #[test]
    fn division_rounds_half_to_even() {
        let b = |x: u64| BigUint::from(x);
        assert_eq!(div_to_f64(&b(1), &b(3)), 1.0 / 3.0);
        assert_eq!(div_to_f64(&b(2), &b(3)), 2.0 / 3.0);
        assert_eq!(div_to_f64(&b((1 << 53) + 1), &b(1)), 9007199254740992.0);
        assert_eq!(div_to_f64(&b((1 << 53) + 3), &b(1)), 9007199254740996.0);
        assert_eq!(div_to_f64(&b(7), &b(1 << 40)), 7.0 / (1u64 << 40) as f64);
    }

    #[test]
    fn the_geometric_interval_brackets_its_mean() {
        let ci = geo_ci(&[1.1, 1.2, 1.15, 1.3, 1.05], 0.95);
        let (g, lo, hi) = (ci.geomean.unwrap(), ci.lo.unwrap(), ci.hi.unwrap());
        assert!(lo < g && g < hi);
        assert_eq!(ci.min, Some(1.05));
        assert_eq!(ci.max, Some(1.3));
        assert_eq!(geo_ci(&[], 0.95).geomean, None);
        assert_eq!(geo_ci(&[2.0], 0.95).lo, None);
    }
}
