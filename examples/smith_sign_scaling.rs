//! Does the Prop. 4 sign bias of Smith's construction shrink like 1/√p?
//! (`research/smith_sign_scaling_20261016/PROTOCOL.md`).
//!
//! ```text
//! cargo run --release --example smith_sign_scaling -- run --p 1009 --n 1000000 --out DIR
//! cargo run --release --example smith_sign_scaling -- verdict --dir DIR 101 1009 10007
//! ```
//!
//! `run` draws random octics (uniform coefficients, non-zero leading
//! coefficient), keeps the squarefree ones, and regresses the Prop. 4 outcome
//! of every non-degenerate Prop. 2 pass on `x = a1 / (2p + 2)`, `a1 = #H(F_p) −
//! (p + 1)`, with curve-clustered (CR1) standard errors. It writes
//! `scaling_p<p>.json`. `verdict` fits `ln β` on `ln p` and writes
//! `verdict.json`.

use crypto_lib::cryptanalysis::smith_prop24 as sp;
use rayon::prelude::*;
use serde_json::{json, Value};
use std::path::PathBuf;
use std::{fs, time::Instant};

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter()
        .position(|a| a == k)
        .and_then(|i| args.get(i + 1).cloned())
}

const SEED: u64 = 20261016;

/// `#H(F_p)` for `y^2 = f`, `deg f = 8`, using a table of squares.
fn point_count(f: &[u64], p: u64, sq: &[i8]) -> i64 {
    let mut n = 0i64;
    for x in 0..p {
        let mut v = 0u64;
        for &c in f.iter().rev() {
            v = (v * x + c) % p;
        }
        n += 1 + sq[v as usize] as i64;
    }
    n + 1 + sq[f[8] as usize] as i64
}

struct Obs {
    x: f64,
    k1: u32,
    m: u32,
}

fn run(p: u64, n: u64) -> Value {
    let t0 = Instant::now();
    let mut sq = vec![-1i8; p as usize];
    sq[0] = 0;
    for a in 1..p {
        sq[(a * a % p) as usize] = 1;
    }
    let seed = SEED + p;
    let per: Vec<Option<(i64, u32, u32, u32, bool)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut rng = sp::SplitMix64(seed ^ i.wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let mut f: Vec<u64> = (0..8).map(|_| rng.below(p)).collect();
            f.push(1 + rng.below(p - 1));
            let mut arng = sp::SplitMix64(seed.wrapping_add(i));
            let o = sp::analyze_curve(&f, p, &mut arng)?;
            let n1 = point_count(&f, p, &sq);
            let ctl = i < 4000 && n1 as u64 != sp::h_point_count(&f, p);
            Some((
                n1,
                o.n_trig() as u32,
                o.n_iso(false) as u32,
                o.n_iso(true) as u32,
                ctl,
            ))
        })
        .collect();
    let mut obs = vec![];
    let (mut sqf, mut trig, mut degen, mut over, mut ctl_bad, mut ctl_checked) =
        (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    // bins of t = a1 / (2 sqrt p), width 0.25 over [-3, 3]
    let mut bins = vec![[0u64; 2]; 24];
    let s = (p as f64).sqrt();
    for (i, r) in per.iter().enumerate() {
        let Some((n1, nt, k1, k0, bad)) = *r else {
            continue;
        };
        sqf += 1;
        if sqf <= 1000 && (i as u64) < 4000 {
            ctl_checked += 1;
            ctl_bad += bad as u64;
        }
        trig += nt as u64;
        if k1 + k0 > nt {
            over += 1;
        }
        degen += (nt.saturating_sub(k1 + k0)) as u64;
        let m = k1 + k0;
        if m == 0 {
            continue;
        }
        let a1 = n1 - (p as i64 + 1);
        obs.push(Obs {
            x: a1 as f64 / (2 * p + 2) as f64,
            k1,
            m,
        });
        let t = a1 as f64 / (2.0 * s);
        let b = (((t + 3.0) / 0.25).floor() as i64).clamp(0, 23) as usize;
        bins[b][0] += m as u64;
        bins[b][1] += k1 as u64;
    }
    // OLS y = a + b x with curve-clustered CR1 errors.
    let (mut sw, mut sx, mut sxx, mut sy, mut sxy) = (0f64, 0f64, 0f64, 0f64, 0f64);
    for o in &obs {
        let m = o.m as f64;
        sw += m;
        sx += m * o.x;
        sxx += m * o.x * o.x;
        sy += o.k1 as f64;
        sxy += o.k1 as f64 * o.x;
    }
    let det = sw * sxx - sx * sx;
    let beta = (sw * sxy - sx * sy) / det;
    let alpha = (sy - beta * sx) / sw;
    let inv = [[sxx / det, -sx / det], [-sx / det, sw / det]];
    let mut meat = [[0f64; 2]; 2];
    for o in &obs {
        let e = o.k1 as f64 - o.m as f64 * (alpha + beta * o.x);
        let g = [e, e * o.x];
        for a in 0..2 {
            for b in 0..2 {
                meat[a][b] += g[a] * g[b];
            }
        }
    }
    let gcl = obs.len() as f64;
    let adj = gcl / (gcl - 1.0) * (sw - 1.0) / (sw - 2.0);
    let mut v = [[0f64; 2]; 2];
    for a in 0..2 {
        for b in 0..2 {
            for c in 0..2 {
                for d in 0..2 {
                    v[a][b] += inv[a][c] * meat[c][d] * inv[d][b];
                }
            }
        }
    }
    let (se_a, se_b) = ((v[0][0] * adj).sqrt(), (v[1][1] * adj).sqrt());
    let alpha_ok = (alpha - 0.5).abs() <= 0.02;
    let valid = alpha_ok && ctl_bad == 0 && over == 0;
    let bins_json: Vec<Value> = bins
        .iter()
        .enumerate()
        .filter(|(_, b)| b[0] > 0)
        .map(|(i, b)| {
            json!({"t_lo": -3.0 + 0.25 * i as f64, "t_hi": -2.75 + 0.25 * i as f64,
                   "n": b[0], "pass": b[1], "rate": b[1] as f64 / b[0] as f64})
        })
        .collect();
    let git = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default();
    json!({
        "p": p, "n_drawn": n, "seed": seed, "git_head": git,
        "squarefree": sqf, "prop2_passes": trig, "degenerate": degen,
        "nondegenerate": sw as u64, "curves_with_nondegenerate": obs.len(),
        "alpha": alpha, "alpha_se": se_a, "beta": beta, "beta_se": se_b,
        "beta_ci": [beta - 1.96 * se_b, beta + 1.96 * se_b],
        "beta_ci_contains_1": (beta - 1.96 * se_b) <= 1.0 && 1.0 <= (beta + 1.96 * se_b),
        "sd_x": ((sxx / sw) - (sx / sw).powi(2)).sqrt(),
        "controls": {
            "alpha_within_0.02": alpha_ok,
            "point_count_checked": ctl_checked, "point_count_mismatches": ctl_bad,
            "curves_k0_plus_k1_over_trig": over,
        },
        "valid": valid,
        "bins_t": bins_json,
        "seconds": t0.elapsed().as_secs_f64(),
        "threads": rayon::current_num_threads(),
    })
}

fn verdict(dir: &std::path::Path, ps: &[u64]) -> Value {
    let rows: Vec<Value> = ps
        .iter()
        .map(|p| {
            serde_json::from_str(
                &fs::read_to_string(dir.join(format!("scaling_p{p}.json"))).unwrap(),
            )
            .unwrap()
        })
        .collect();
    let valid = rows.iter().all(|r| r["valid"].as_bool() == Some(true));
    let b: Vec<(f64, f64, f64)> = rows
        .iter()
        .map(|r| {
            (
                (r["p"].as_u64().unwrap() as f64).ln(),
                r["beta"].as_f64().unwrap(),
                r["beta_se"].as_f64().unwrap(),
            )
        })
        .collect();
    if b.iter().any(|t| t.1 <= 0.0) {
        return json!({"valid": valid, "h3": "inconclusive: non-positive β"});
    }
    // WLS of ln beta on ln p.
    let (mut sw, mut sx, mut sy, mut sxx, mut sxy) = (0f64, 0f64, 0f64, 0f64, 0f64);
    for &(lp, be, se) in &b {
        let w = (be / se).powi(2);
        let y = be.ln();
        sw += w;
        sx += w * lp;
        sy += w * y;
        sxx += w * lp * lp;
        sxy += w * lp * y;
    }
    let det = sw * sxx - sx * sx;
    let gamma = (sw * sxy - sx * sy) / det;
    let gamma_se = (sw / det).sqrt();
    let (lo, hi) = (gamma - 1.96 * gamma_se, gamma + 1.96 * gamma_se);
    let h3 = if !valid {
        "control failed - invalid"
    } else if hi < 0.25 {
        "bias shrinks like 1/√p (β bounded)"
    } else if lo > 0.25 {
        "bias persists beyond 1/√p"
    } else {
        "inconclusive"
    };
    json!({
        "primes": ps, "valid": valid,
        "gamma": gamma, "gamma_se": gamma_se, "gamma_ci": [lo, hi],
        "h3": h3,
        "h1prime_beta_ci_contains_1": rows.iter().map(|r| (r["p"].clone(), r["beta_ci_contains_1"].clone())).collect::<Vec<_>>(),
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("run") => {
            let p: u64 = arg(&args, "--p").expect("--p").parse().unwrap();
            let n: u64 = arg(&args, "--n").expect("--n").parse().unwrap();
            let out = PathBuf::from(arg(&args, "--out").expect("--out DIR"));
            fs::create_dir_all(&out).unwrap();
            let r = run(p, n);
            fs::write(
                out.join(format!("scaling_p{p}.json")),
                serde_json::to_string_pretty(&r).unwrap() + "\n",
            )
            .unwrap();
            eprintln!(
                "p={p}: beta {:.3} ± {:.3}, alpha {:.4}, valid {}, {:.0}s",
                r["beta"].as_f64().unwrap(),
                1.96 * r["beta_se"].as_f64().unwrap(),
                r["alpha"].as_f64().unwrap(),
                r["valid"],
                r["seconds"].as_f64().unwrap()
            );
        }
        Some("verdict") => {
            let dir = PathBuf::from(arg(&args, "--dir").expect("--dir DIR"));
            let ps: Vec<u64> = args[4..].iter().map(|s| s.parse().unwrap()).collect();
            let v = verdict(&dir, &ps);
            fs::write(
                dir.join("verdict.json"),
                serde_json::to_string_pretty(&v).unwrap() + "\n",
            )
            .unwrap();
            println!("{}", serde_json::to_string_pretty(&v).unwrap());
        }
        _ => eprintln!(
            "usage: smith_sign_scaling run --p P --n N --out DIR | verdict --dir DIR P..."
        ),
    }
}
