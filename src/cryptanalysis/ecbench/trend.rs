//! How a paired ratio moves across curves: `ln ratio` against `ln r` and
//! the field degree `n`, from `ecbench compare --save` files.
//!
//! Each curve contributes its paired ratio `B / A` and the cluster
//! bootstrap interval `compare` recorded; the standard error of
//! `ln ratio` is read from that interval (`(ln hi − ln lo) / 3.92`).
//! Three weighted least-squares fits:
//!
//! - `free_r`: `ln ratio = c + β ln r`;
//! - `free_r_n`: `ln ratio = c + β ln r + γ ln n` (needs four curves whose
//!   `n` and `r` do not move together);
//! - `model`: `β`, `γ` fixed (for Frobenius-orbit index calculus against a
//!   Frobenius rho, `β = 1/6`, `γ = −1/2`), only `c` fitted, with the
//!   `χ²` per degree of freedom saying whether the curves agree with it.
//!
//! Coefficient intervals are `±1.96` WLS standard errors, inflated by
//! `√(χ²/dof)` when that exceeds one.  The fit is descriptive: a handful
//! of toy curves, one host.

use serde::Serialize;
use serde_json::Value;

#[derive(Clone, Debug, Serialize)]
pub struct Point {
    pub slug: String,
    pub n: u32,
    pub log2_r: f64,
    pub ratio: f64,
    pub ci95: [f64; 2],
    pub se_ln: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Fit {
    pub name: String,
    pub formula: String,
    /// Fitted coefficients by name, with their 95 % intervals.
    pub coefficients: Vec<(String, f64, [f64; 2])>,
    pub chi2_per_dof: Option<f64>,
    pub dof: i64,
    /// Fitted ratio per curve, in input order.
    pub fitted: Vec<f64>,
}

fn n_of(slug: &str) -> Result<u32, String> {
    let rest = slug
        .strip_prefix("icv1-f2m")
        .ok_or_else(|| format!("{slug}: not a binary-field ICV1 slug"))?;
    let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
    digits
        .parse()
        .map_err(|_| format!("{slug}: no field degree"))
}

/// The per-curve points of comparison files (`ecbench.comparison/v1`).
pub fn points(docs: &[Value]) -> Result<Vec<Point>, String> {
    let mut out = Vec::new();
    for d in docs {
        let curves = d["curves"].as_array().ok_or("comparison without curves")?;
        for c in curves {
            let slug = c["slug"].as_str().ok_or("curve without slug")?.to_string();
            let ratio = c["ratio_b_over_a"].as_f64().ok_or("curve without ratio")?;
            let lo = c["ci95"][0].as_f64().ok_or("curve without ci95")?;
            let hi = c["ci95"][1].as_f64().ok_or("curve without ci95")?;
            if !(ratio > 0.0 && lo > 0.0 && hi >= lo) {
                return Err(format!("{slug}: ratio or interval not positive"));
            }
            // A degenerate interval would get infinite weight.
            let se_ln = ((hi.ln() - lo.ln()) / 3.92).max(1e-3);
            out.push(Point {
                n: n_of(&slug)?,
                slug,
                log2_r: c["log2_r"].as_f64().ok_or("curve without log2_r")?,
                ratio,
                ci95: [lo, hi],
                se_ln,
            });
        }
    }
    Ok(out)
}

/// Solve the symmetric positive system `a x = b` (Gauss–Jordan, small).
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let k = b.len();
    let mut inv: Vec<Vec<f64>> = (0..k)
        .map(|i| (0..k).map(|j| f64::from(u8::from(i == j))).collect())
        .collect();
    for col in 0..k {
        let piv = (col..k).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[piv][col].abs() < 1e-12 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        inv.swap(col, piv);
        let p = a[col][col];
        for j in 0..k {
            a[col][j] /= p;
            inv[col][j] /= p;
        }
        b[col] /= p;
        for i in 0..k {
            if i != col {
                let f = a[i][col];
                for j in 0..k {
                    a[i][j] -= f * a[col][j];
                    inv[i][j] -= f * inv[col][j];
                }
                b[i] -= f * b[col];
            }
        }
    }
    Some((b, inv))
}

/// Weighted least squares of `y − offset` on the given regressors.
fn wls(
    name: &str,
    formula: &str,
    pts: &[Point],
    regressors: &[(&str, &dyn Fn(&Point) -> f64)],
    offset: &dyn Fn(&Point) -> f64,
) -> Option<Fit> {
    let k = regressors.len();
    let m = pts.len();
    if m < k {
        return None;
    }
    let mut xtwx = vec![vec![0.0; k]; k];
    let mut xtwy = vec![0.0; k];
    for p in pts {
        let w = 1.0 / (p.se_ln * p.se_ln);
        let y = p.ratio.ln() - offset(p);
        let x: Vec<f64> = regressors.iter().map(|(_, f)| f(p)).collect();
        for i in 0..k {
            xtwy[i] += w * x[i] * y;
            for j in 0..k {
                xtwx[i][j] += w * x[i] * x[j];
            }
        }
    }
    let (beta, cov) = solve(xtwx, xtwy)?;
    let mut chi2 = 0.0;
    let mut fitted = Vec::new();
    for p in pts {
        let yhat: f64 = regressors
            .iter()
            .zip(&beta)
            .map(|((_, f), b)| b * f(p))
            .sum::<f64>()
            + offset(p);
        chi2 += ((p.ratio.ln() - yhat) / p.se_ln).powi(2);
        fitted.push(yhat.exp());
    }
    let dof = m as i64 - k as i64;
    let chi2_per_dof = (dof > 0).then(|| chi2 / dof as f64);
    let inflate = chi2_per_dof.map_or(1.0, |c| c.max(1.0).sqrt());
    let coefficients = regressors
        .iter()
        .enumerate()
        .map(|(i, (label, _))| {
            let half = 1.96 * cov[i][i].sqrt() * inflate;
            (label.to_string(), beta[i], [beta[i] - half, beta[i] + half])
        })
        .collect();
    Some(Fit {
        name: name.into(),
        formula: formula.into(),
        coefficients,
        chi2_per_dof,
        dof,
        fitted,
    })
}

/// The three fits; `model` carries the fixed `β` and `γ`.
pub fn fits(pts: &[Point], beta: f64, gamma: f64) -> Vec<Fit> {
    let one = |_: &Point| 1.0;
    let ln_r = |p: &Point| p.log2_r * std::f64::consts::LN_2;
    let ln_n = |p: &Point| f64::from(p.n).ln();
    let zero = |_: &Point| 0.0;
    let fixed = move |p: &Point| beta * ln_r(p) + gamma * ln_n(p);
    let mut out = Vec::new();
    if let Some(f) = wls(
        "free_r",
        "ln ratio = c + beta ln r",
        pts,
        &[("c", &one), ("beta", &ln_r)],
        &zero,
    ) {
        out.push(f);
    }
    if pts.len() >= 4 {
        if let Some(f) = wls(
            "free_r_n",
            "ln ratio = c + beta ln r + gamma ln n",
            pts,
            &[("c", &one), ("beta", &ln_r), ("gamma", &ln_n)],
            &zero,
        ) {
            out.push(f);
        }
    }
    if let Some(f) = wls(
        "model",
        &format!("ln ratio = c + {beta} ln r + {gamma} ln n"),
        pts,
        &[("c", &one)],
        &fixed,
    ) {
        out.push(f);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_power_laws_are_recovered() {
        let pts: Vec<Point> = [
            (13u32, 11.0),
            (17, 16.0),
            (19, 17.0),
            (23, 21.0),
            (29, 15.4),
        ]
        .iter()
        .map(|&(n, l2r)| {
            let r: f64 = 2f64.powf(l2r);
            let ratio = 0.3 * r.powf(1.0 / 6.0) / f64::from(n).sqrt();
            Point {
                slug: format!("icv1-f2m{n}-x"),
                n,
                log2_r: l2r,
                ratio,
                ci95: [ratio * 0.9, ratio * 1.1],
                se_ln: 0.05,
            }
        })
        .collect();
        let f = fits(&pts, 1.0 / 6.0, -0.5);
        let rn = f.iter().find(|f| f.name == "free_r_n").unwrap();
        assert!((rn.coefficients[1].1 - 1.0 / 6.0).abs() < 1e-9);
        assert!((rn.coefficients[2].1 + 0.5).abs() < 1e-9);
        let model = f.iter().find(|f| f.name == "model").unwrap();
        assert!((model.coefficients[0].1 - 0.3f64.ln()).abs() < 1e-9);
        assert!(model.chi2_per_dof.unwrap() < 1e-12);
        assert_eq!(n_of("icv1-f2m41-tm2308219-7f48b14a").unwrap(), 41);
    }
}
