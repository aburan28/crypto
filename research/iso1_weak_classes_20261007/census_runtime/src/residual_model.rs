//! Native evaluation of the registered ISO-1 residual-zero model V2
//! (`../residual_model_20261016/PROTOCOL_V2.md`), added 2026-10-16.
//!
//! N_8(t) = H_odd(t) * L_2, with k = v2(f_pi) >= 2, s = (D_K/2) and
//!   L_2 = (2 - s) 2^(k-2)                    if t != +-(Q+1) (mod 32),
//!   L_2 = 5 (2 - s) 2^(k-2) - 6 (1 - s)      if t == +-(Q+1) (mod 32).
//! Weak orbits ~ Binomial(N_8/12, r), r = 8/(p^2 - 1); P(zero) = (1-r)^(N_8/12).
//!
//! usage:
//!   iso1_residual_model exact OMEGA8.csv LEVELS.csv
//!   iso1_residual_model holdout P CENSUS.csv LEVELS.csv
//!   iso1_residual_model predict P LEVELS.csv
//! LEVELS.csv is the output of `class_levels.gp`; OMEGA8.csv of
//! `iso1_lambda_census` with 8th powers; CENSUS.csv of `iso1_class_census`.
use std::collections::BTreeMap;
use std::fs;

/// A rational with denominator dividing 6, stored as 6 * value.
type Sixths = i128;

fn sixths(s: &str) -> Sixths {
    match s.split_once('/') {
        Some((n, d)) => {
            let (n, d): (i128, i128) = (n.parse().unwrap(), d.parse().unwrap());
            assert_eq!(6 % d, 0, "denominator {d} does not divide 6");
            n * (6 / d)
        }
        None => s.parse::<i128>().unwrap() * 6,
    }
}

fn v2(n: i128) -> u32 {
    n.unsigned_abs().trailing_zeros()
}

struct Level {
    t: i128,
    dk: i128,
    f_pi: i128,
    h_odd: Sixths,
}

fn levels(path: &str) -> Vec<Level> {
    let text = fs::read_to_string(path).expect("levels csv");
    let mut lines = text.lines();
    let head: Vec<_> = lines.next().unwrap().split(',').collect();
    let col = |name: &str| head.iter().position(|h| *h == name).unwrap();
    let (ct, cd, cf, ch) = (col("t"), col("D_K"), col("f_pi"), col("H_odd"));
    lines
        .map(|l| {
            let c: Vec<_> = l.split(',').collect();
            Level {
                t: c[ct].parse().unwrap(),
                dk: c[cd].parse().unwrap(),
                f_pi: c[cf].parse().unwrap(),
                h_odd: sixths(c[ch]),
            }
        })
        .collect()
}

/// Registered N_8(t), times 6; `None` when v2(f_pi) < 2.
fn n8_sixths(p: i128, row: &Level) -> Option<Sixths> {
    let q6 = p.pow(6);
    let k = v2(row.f_pi);
    if k < 2 {
        return None;
    }
    let s: i128 = if row.dk.rem_euclid(2) == 0 {
        0
    } else if row.dk.rem_euclid(8) == 1 {
        1
    } else {
        -1
    };
    let deep = (row.t - (q6 + 1)).rem_euclid(32) == 0 || (row.t + q6 + 1).rem_euclid(32) == 0;
    let base = (2 - s) << (k - 2);
    let l2 = if deep { 5 * base - 6 * (1 - s) } else { base };
    Some(row.h_odd * l2)
}

fn exact(omega: &str, lev: &str) {
    let text = fs::read_to_string(omega).expect("omega csv");
    let mut lines = text.lines();
    let head: Vec<_> = lines.next().unwrap().split(',').collect();
    let col = |name: &str| head.iter().position(|h| *h == name).unwrap();
    let (cp, ce, ct, co, cc) = (
        col("p"),
        col("e"),
        col("trace"),
        col("in_omega"),
        col("count"),
    );
    let mut p = 0i128;
    let mut n: BTreeMap<i128, i128> = BTreeMap::new();
    for l in lines {
        let c: Vec<_> = l.split(',').collect();
        p = c[cp].parse().unwrap();
        assert_eq!(c[ce], "3", "exactness test needs 8th powers (e = 3)");
        assert_eq!(c[co], "1");
        *n.entry(c[ct].parse::<i128>().unwrap().abs()).or_default() +=
            c[cc].parse::<i128>().unwrap();
    }
    let rows = levels(lev);
    let (mut bad, mut count) = (0, 0);
    for r in rows.iter().filter(|r| r.t % p != 0) {
        count += 1;
        let pred = n8_sixths(p, r).unwrap_or(0);
        let obs = 6 * n.get(&r.t).copied().unwrap_or(0);
        if pred != obs {
            bad += 1;
            if bad <= 5 {
                println!(
                    "  mismatch t={} predicted={} counted={}",
                    r.t,
                    pred as f64 / 6.0,
                    obs / 6
                );
            }
        }
    }
    let listed: std::collections::BTreeSet<i128> = rows.iter().map(|r| r.t).collect();
    let nonord: Vec<_> = n.iter().filter(|(t, _)| *t % p == 0).collect();
    let extra = n
        .keys()
        .filter(|t| *t % p != 0 && !listed.contains(t))
        .count();
    println!(
        "p={p}: {count} ordinary rows, {bad} mismatches, {extra} unexpected traces, nonordinary traces with 8th-power lambda {nonord:?} -> {}",
        if bad == 0 && extra == 0 { "PASS" } else { "FAIL" }
    );
}

/// Rank-based ROC AUC with average ranks for ties.
fn auc(scores: &[f64], labels: &[bool]) -> f64 {
    let mut pairs: Vec<(f64, bool)> = scores.iter().copied().zip(labels.iter().copied()).collect();
    pairs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let pos = labels.iter().filter(|&&l| l).count() as f64;
    let neg = labels.len() as f64 - pos;
    let (mut rank_sum, mut i) = (0.0, 0);
    while i < pairs.len() {
        let mut j = i;
        while j < pairs.len() && pairs[j].0 == pairs[i].0 {
            j += 1;
        }
        let avg = (i + j + 1) as f64 / 2.0;
        rank_sum += avg * pairs[i..j].iter().filter(|x| x.1).count() as f64;
        i = j;
    }
    (rank_sum - pos * (pos + 1.0) / 2.0) / (pos * neg)
}

fn census(path: &str) -> BTreeMap<i128, i128> {
    let text = fs::read_to_string(path).expect("census csv");
    let mut lines = text.lines();
    let head: Vec<_> = lines.next().unwrap().split(',').collect();
    let col = |name: &str| head.iter().position(|h| *h == name).unwrap();
    let (ct, cw, cs) = (
        col("trace"),
        col("weak_representatives"),
        col("trace_status"),
    );
    lines
        .map(|l| l.split(',').collect::<Vec<_>>())
        .filter(|c| c[cs] == "ordinary")
        .map(|c| (c[ct].parse().unwrap(), c[cw].parse().unwrap()))
        .collect()
}

struct Prediction {
    p_zero: Vec<f64>,
    rows: usize,
    expected: f64,
    sd: f64,
    weak_total: f64,
}

fn predict_rows(p: i128, rows: &[Level], mut each: impl FnMut(&Level, f64)) -> Prediction {
    let q = p * p;
    let r = 8.0 / (q - 1) as f64;
    let mut out = Prediction {
        p_zero: Vec::new(),
        rows: 0,
        expected: 0.0,
        sd: 0.0,
        weak_total: 0.0,
    };
    let mut var = 0.0;
    for row in rows {
        let Some(n6) = n8_sixths(p, row) else {
            continue;
        };
        let n8 = n6 as f64 / 6.0;
        let orbits = n8 / 12.0;
        let pz = (1.0 - r).powf(orbits);
        out.rows += 1;
        out.weak_total += r * n8;
        out.expected += pz;
        var += pz * (1.0 - pz);
        out.p_zero.push(pz);
        each(row, orbits);
    }
    out.sd = var.sqrt();
    out
}

fn holdout(p: i128, census_path: &str, lev: &str) {
    let q = p * p;
    let r = 8.0 / (q - 1) as f64;
    let w = census(census_path);
    let rows = levels(lev);
    let mut zero = Vec::new();
    let mut chi = Vec::new();
    let pred = predict_rows(p, &rows, |row, orbits| {
        let weak = w[&row.t];
        zero.push(weak == 0);
        let mean = orbits * r;
        if mean >= 5.0 {
            let x = weak as f64 / 12.0 - mean;
            chi.push(x * x / (orbits * r * (1.0 - r)));
        }
    });
    let obs = zero.iter().filter(|&&z| z).count();
    let (lo, hi) = (
        pred.expected - 1.96 * pred.sd,
        pred.expected + 1.96 * pred.sd,
    );
    let total_obs: i128 = w.iter().filter(|(t, _)| **t > 0).map(|(_, c)| c).sum();
    let pass = lo <= obs as f64 && obs as f64 <= hi;
    println!(
        "p={p} rows={} zeros obs={obs} exp={:.1} [{lo:.1},{hi:.1}] {} dispersion={:.3} (n={}) auc={:.4} weak_total obs={total_obs} pred={:.0}",
        pred.rows,
        pred.expected,
        if pass { "PASS" } else { "FAIL" },
        chi.iter().sum::<f64>() / chi.len() as f64,
        chi.len(),
        auc(&pred.p_zero, &zero),
        pred.weak_total
    );
    let mut pts: Vec<(f64, bool)> = pred
        .p_zero
        .iter()
        .copied()
        .zip(zero.iter().copied())
        .filter(|x| x.0 > 1e-6)
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = pts.len();
    for k in 0..5 {
        let ch = &pts[k * m / 5..(k + 1) * m / 5];
        println!(
            "   P(zero) in [{:.2e},{:.2e}]: exp={:.1} obs={}",
            ch[0].0,
            ch[ch.len() - 1].0,
            ch.iter().map(|x| x.0).sum::<f64>(),
            ch.iter().filter(|x| x.1).count()
        );
    }
    let tail = pred
        .p_zero
        .iter()
        .zip(&zero)
        .filter(|(a, z)| **a <= 1e-6 && **z)
        .count();
    println!("   P(zero) <= 1e-6: obs={tail}");
}

fn predict(p: i128, lev: &str) {
    let q = p * p;
    let pred = predict_rows(p, &levels(lev), |_, _| {});
    println!(
        "p={p} high_depth_rows={} predicted_zero_rows={:.1} interval95=[{:.1},{:.1}] predicted_weak_total={:.0} exact_weak_total={}",
        pred.rows,
        pred.expected,
        pred.expected - 1.96 * pred.sd,
        pred.expected + 1.96 * pred.sd,
        pred.weak_total,
        q * q + q
    );
}

/// Pooled observed-minus-predicted zero rows over scored primes (optionally only
/// those listed after the path, e.g. the registered follow-up primes), overall and
/// grouped by e = v2(p^2 - 1). Reads figures/zero_rows.csv; a row is scored
/// when its observed column is non-empty and its role is pristine or forward.
/// Each prime's variance is ((hi - lo) / 3.92)^2, its 95% interval half-width.
fn pool(path: &str, only: &[i128]) {
    let text = fs::read_to_string(path).expect("zero_rows csv");
    let mut groups: BTreeMap<u32, Vec<(i128, f64, f64)>> = BTreeMap::new();
    for l in text.lines().skip(1) {
        let c: Vec<&str> = l.split(',').collect();
        if c[2].is_empty() || !(c[1] == "pristine" || c[1] == "forward") {
            continue;
        }
        let p: i128 = c[0].parse().unwrap();
        if !only.is_empty() && !only.contains(&p) {
            continue;
        }
        let (obs, exp): (f64, f64) = (c[2].parse().unwrap(), c[3].parse().unwrap());
        let (lo, hi): (f64, f64) = (c[4].parse().unwrap(), c[5].parse().unwrap());
        let sd = (hi - lo) / 3.92;
        groups
            .entry(v2(p * p - 1))
            .or_default()
            .push((p, obs - exp, sd));
    }
    let summarize = |label: String, rows: Vec<&(i128, f64, f64)>| {
        let d: f64 = rows.iter().map(|r| r.1).sum();
        let sd = rows.iter().map(|r| r.2 * r.2).sum::<f64>().sqrt();
        let ps: Vec<String> = rows.iter().map(|r| r.0.to_string()).collect();
        println!(
            "{label}: primes [{}] obs-exp = {d:+.1}, sd = {sd:.1}, z = {:.2}",
            ps.join(","),
            d / sd
        );
    };
    summarize("all".into(), groups.values().flatten().collect());
    for (e, rows) in &groups {
        summarize(format!("e = {e}"), rows.iter().collect());
    }
    let hi: Vec<_> = groups
        .iter()
        .filter(|(e, _)| **e >= 4)
        .flat_map(|(_, r)| r)
        .collect();
    summarize("e >= 4".into(), hi);
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    match a.first().map(String::as_str) {
        Some("exact") if a.len() == 3 => exact(&a[1], &a[2]),
        Some("holdout") if a.len() == 4 => holdout(a[1].parse().unwrap(), &a[2], &a[3]),
        Some("predict") if a.len() == 3 => predict(a[1].parse().unwrap(), &a[2]),
        Some("pool") if a.len() >= 2 => {
            let only: Vec<i128> = a[2..].iter().map(|x| x.parse().unwrap()).collect();
            pool(&a[1], &only)
        }
        _ => panic!("usage: iso1_residual_model exact OMEGA8 LEVELS | holdout P CENSUS LEVELS | predict P LEVELS"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_factor_matches_exhaustive_table() {
        // Cells read from the exhaustive p = 7, 11, 13, 17 censuses (H_odd = 1).
        let row = |t: i128, dk: i128, f: i128| Level {
            t,
            dk,
            f_pi: f,
            h_odd: 6,
        };
        let p: i128 = 11;
        let q6 = p.pow(6);
        let deep_t = q6 + 1 - 32 * 7; // t == Q+1 (mod 32)
        assert_eq!(n8_sixths(p, &row(6, -40, 4)), Some(6 * 2)); // k=2, s=0, shallow
        assert_eq!(n8_sixths(p, &row(6, -3, 8)), Some(6 * 6)); // k=3, s=-1
        assert_eq!(n8_sixths(p, &row(deep_t, -40, 8)), Some(6 * 14)); // k=3, s=0, deep
        assert_eq!(n8_sixths(p, &row(deep_t, -7, 16)), Some(6 * 20)); // k=4, s=1, deep
        assert_eq!(n8_sixths(p, &row(deep_t, -3, 32)), Some(6 * 108)); // k=5, s=-1, deep
        assert_eq!(n8_sixths(p, &row(6, -40, 2)), None); // k=1: proved empty
    }
}
