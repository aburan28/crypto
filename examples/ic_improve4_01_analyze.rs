//! Analysis for IC improvement round 3, item 01 (the relation phase of genus-3 IC after the cubic
//! GHS transfer), adapted from `ic_improve2_06_analyze.rs`.
//!
//! - `gate <run_dir> <item06_run_dir>`: every `sched:lpsgebw2` record (b <= 21) must equal item 06's.
//! - `gatetrain <train_dir> <item06_train_dir>`: relation-only incumbent reruns must equal item 06's
//!   training relation costs on the same (b, curve, kappa).
//! - `diag <train_dir> <out.json>`: step-1 decomposition of the relation phase (trials x cost per
//!   trial; primitives; smoothness and 2LP yield; q^a (log q)^beta against q^c fits).
//! - `r2rule <train_dir> <item06_train_dir> <out.json>`: the frozen R2 kappa rule (PROTOCOL.md).
//! - `acct <run_dir>`: `+R1s` (R1 driving the walk) against the shadow-accounted `+R1` arm.
//! - `logfit <run_dir> <arm> <metric> <lo> <hi>`: q^a (log q)^beta against q^c on test records.
//! - `analyze <run_dir> <lo> <hi> <upper_lo> <out.json> <arm>...`: exponents (whole, relation, LA)
//!   with paired Deltas against the first arm (the incumbent), local exponents, their trend and the
//!   paired trend difference, paired ratios per size. 2000 bootstrap resamples of curves within each
//!   size, paired across arms; a size uses the curves on which every listed arm is correct.
//!
//! ```sh
//! cargo run --release --example ic_improve3_01_analyze -- diag research/ic_improve3_20261015/01-genus3-relations/results/train d.json
//! ```
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;

fn read_jsonl(dir: &str) -> Vec<Value> {
    let mut out = vec![];
    let mut names: Vec<_> = fs::read_dir(dir)
        .unwrap()
        .flatten()
        .map(|e| e.path())
        .collect();
    names.sort();
    for p in names {
        if p.extension().map(|e| e == "jsonl").unwrap_or(false) {
            for line in fs::read_to_string(&p).unwrap().lines() {
                if !line.trim().is_empty() {
                    out.push(serde_json::from_str(line).unwrap());
                }
            }
        }
    }
    out
}

/// (arm, bits, curve) -> result record (with "q" copied in).
type Data = BTreeMap<(String, u64, u64), Value>;

fn load(dir: &str) -> Data {
    let mut d = Data::new();
    for rec in read_jsonl(dir) {
        let b = rec["bits"].as_u64().unwrap();
        let c = rec["curve"].as_u64().unwrap();
        for r in rec["results"].as_array().unwrap() {
            let mut r = r.clone();
            r["q"] = rec["q"].clone();
            r["ell"] = rec["ell"].clone();
            let mut arm = r["arm"].as_str().unwrap().to_string();
            // the cx of every arm here is a step function of b: one arm name over all sizes
            // (`sched:lpsgebw2@16+R1` -> `sched:lpsgebw2+R1`)
            // round 4: only the schedule's cx(b) is dropped from the name; any other cx stays
            if let Some((base, rest)) = arm.split_once('@') {
                let suffix = rest.find('+').map(|i| &rest[i..]).unwrap_or("");
                let cxs = &rest[..rest.len() - suffix.len()];
                let sched = if b <= 17 { 16.0 } else { 8.0 };
                if cxs.parse::<f64>().map(|x| x == sched).unwrap_or(false) {
                    arm = format!("{base}{suffix}");
                }
            }
            d.insert((arm, b, c), r);
        }
    }
    d
}

fn f(v: &Value, k: &str) -> f64 {
    v[k].as_f64()
        .unwrap_or_else(|| v[k].as_u64().map(|x| x as f64).unwrap_or(f64::NAN))
}

fn gm(xs: &[f64]) -> f64 {
    (xs.iter().map(|x| x.ln()).sum::<f64>() / xs.len() as f64).exp()
}

fn ols(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    for (a, b) in x.iter().zip(y) {
        sxy += (a - mx) * (b - my);
        sxx += (a - mx) * (a - mx);
    }
    sxy / sxx
}

fn ci(mut v: Vec<f64>) -> [f64; 2] {
    v.retain(|x| x.is_finite());
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if v.is_empty() {
        return [f64::NAN, f64::NAN];
    }
    let p = |q: f64| v[((q * (v.len() - 1) as f64).round() as usize).min(v.len() - 1)];
    [p(0.025), p(0.975)]
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Counted unit: M_q + ceil(log2 q) I_q + c_ell M_ell (`cost_units` of the instrument, c_ell = 1).
fn units(v: &Value, q: f64) -> f64 {
    let a = v.as_array().unwrap();
    let c = |i: usize| a[i].as_u64().unwrap() as f64;
    let bl = 64.0 - ((q as u64) - 1).leading_zeros() as f64;
    c(0) + bl * c(1) + c(2)
}

/// Gate 0: every R1-driven incumbent record `r2:lpsgebw2` must equal round 3's shadow-accounted
/// `r2:lpsgebw2+R1` record on the same (b, curve): k, cost, relation cost, LA counters, steps.
fn gate(run: &str, old: &str) {
    let d = load(run);
    let o = load(old);
    let (mut ok, mut bad) = (0, 0);
    println!("| log2 q | curve | k equal | cost (round 3 R2) | cost (this item) | relation cost equal | LA counters equal | steps equal | identical |");
    println!("|---|---|---|---|---|---|---|---|---|");
    for ((arm, b, c), r) in &d {
        if arm != "r2:lpsgebw2" {
            continue;
        }
        let Some(p) = o.get(&("r2:lpsgebw2+R1".to_string(), *b, *c)) else {
            continue;
        };
        let eq = r["cost"] == p["cost"]
            && r["k"] == p["k"]
            && r["rel_cost"] == p["rel_cost"]
            && r["la_counters"] == p["la_counters"]
            && r["steps"] == p["steps"];
        if eq {
            ok += 1
        } else {
            bad += 1
        }
        println!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            b,
            c,
            r["k"] == p["k"],
            p["cost"],
            r["cost"],
            r["rel_cost"] == p["rel_cost"],
            r["la_counters"] == p["la_counters"],
            r["steps"] == p["steps"],
            if eq { "yes" } else { "NO" }
        );
    }
    println!("\ngate 0: {ok} identical, {bad} different");
}

/// Per-trial cost by primitive (R1, C1, C2), from every record with os = 0 whose arm starts
/// with `prefix` (the base arm; its diag carries all variants' counters). Fits A + B log2 q
/// (OLS, every record a point) and the q^c exponent of the per-trial cost.
fn pertrial(dir: &str, prefix: &str, out: &str) {
    let d = load(dir);
    let mut by_b: BTreeMap<u64, Vec<Value>> = BTreeMap::new();
    for ((arm, b, _c), r) in &d {
        if arm != prefix {
            continue;
        }
        by_b.entry(*b).or_default().push(r.clone());
    }
    let q_of = |r: &Value| f(r, "q");
    let pr = |r: &Value, k: &str| units(&r["diag"][k], q_of(r)) / f(&r["diag"], "trials");
    // variants: (name, walk key, split key, roots key)
    let vars = [
        ("R1", "walk_add_r1", "split_r1", "roots_r1"),
        ("C1", "walk_add_r1", "split_c1", "roots_c1"),
        ("C2", "walk_add_c2", "split_r1", "roots_r1"),
        ("C1+C2", "walk_add_c2", "split_c1", "roots_c1"),
    ];
    println!("| log2 q | n | trials | R1 walk | R1 split | R1 roots | R1 total | C1 split | C1 roots | C1 total | C2 walk | C2 total | C1+C2 total | C1/R1 | C2/R1 | C1+C2/R1 | batched adds / batch |");
    println!("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|");
    let mut rows = vec![];
    let mut pts: BTreeMap<&str, Vec<(f64, f64)>> = BTreeMap::new();
    for (b, rs) in &by_b {
        let mean = |fun: &dyn Fn(&Value) -> f64| rs.iter().map(fun).sum::<f64>() / rs.len() as f64;
        let tot = |r: &Value, v: &(&str, &str, &str, &str)| pr(r, v.1) + pr(r, v.2) + pr(r, v.3);
        for v in &vars {
            for r in rs {
                pts.entry(v.0)
                    .or_default()
                    .push(((q_of(r)).log2(), tot(r, v)));
            }
        }
        let t = |i: usize| mean(&|r| tot(r, &vars[i]));
        let bpb = mean(&|r| f(&r["diag"], "batched_adds") / f(&r["diag"], "batches").max(1.0));
        println!(
            "| {} | {} | {:.3e} | {:.1} | {:.1} | {:.1} | {:.1} | {:.1} | {:.1} | {:.1} | {:.1} | {:.1} | {:.1} | {:.4} | {:.4} | {:.4} | {:.1} |",
            b,
            rs.len(),
            mean(&|r| f(&r["diag"], "trials")),
            mean(&|r| pr(r, "walk_add_r1")),
            mean(&|r| pr(r, "split_r1")),
            mean(&|r| pr(r, "roots_r1")),
            t(0),
            mean(&|r| pr(r, "split_c1")),
            mean(&|r| pr(r, "roots_c1")),
            t(1),
            mean(&|r| pr(r, "walk_add_c2")),
            t(2),
            t(3),
            t(1) / t(0),
            t(2) / t(0),
            t(3) / t(0),
            bpb
        );
        rows.push(json!({"bits": b, "n": rs.len(), "R1": t(0), "C1": t(1), "C2": t(2), "C1C2": t(3),
            "r1_walk": mean(&|r| pr(r, "walk_add_r1")), "r1_split": mean(&|r| pr(r, "split_r1")),
            "r1_roots": mean(&|r| pr(r, "roots_r1")), "c1_split": mean(&|r| pr(r, "split_c1")),
            "c1_roots": mean(&|r| pr(r, "roots_c1")), "c2_walk": mean(&|r| pr(r, "walk_add_c2"))}));
    }
    println!("\nPer-trial cost fitted as A + B log2 q (OLS, every record a point), and the share B / ((A + B b) ln 2) it adds to a local exponent at b = 21:\n");
    println!("| variant | A | B per bit | log share at 2^21 |");
    println!("|---|---|---|---|");
    let mut fits = json!({});
    for v in &vars {
        let p = &pts[v.0];
        let x: Vec<f64> = p.iter().map(|t| t.0).collect();
        let y: Vec<f64> = p.iter().map(|t| t.1).collect();
        let bb = ols(&x, &y);
        let a = y.iter().sum::<f64>() / y.len() as f64 - bb * x.iter().sum::<f64>() / x.len() as f64;
        let share = bb / ((a + 21.0 * bb) * std::f64::consts::LN_2);
        println!("| {} | {:.1} | {:.2} | {:.4} |", v.0, a, bb, share);
        fits[v.0] = json!({"A": a, "B": bb, "share21": share});
    }
    fs::write(
        out,
        serde_json::to_string_pretty(&json!({"rows": rows, "fits": fits})).unwrap(),
    )
    .unwrap();
}

/// Training rules (PROTOCOL.md): W for C2, epsilon for C3, epsilon and kappa for C4.
/// Records: `w` stream (W walks) at kappa = kappa_R2(b) x {grid}, curves 5-7.
fn tune(dir: &str, out: &str) {
    let d = load(dir);
    // (arm without stream, b) -> [(kappa, curve, cost)]
    let mut m: BTreeMap<(String, u64), Vec<(f64, u64, f64)>> = BTreeMap::new();
    for ((arm, b, c), r) in &d {
        if r["correct"] != json!(true) || !arm.starts_with("w:") {
            continue;
        }
        let name = arm[2..].to_string();
        m.entry((name, *b))
            .or_default()
            .push((f(r, "kappa"), *c, f(r, "cost")));
    }
    let bits: BTreeSet<u64> = m.keys().map(|k| k.1).collect();
    let kr2 = |b: u64| (0.2669 - 0.0941 * (b as f64 - 20.0)).exp();
    // grid points: kappa / kappa_R2 rounded to 0.05
    let mult = |k: f64, b: u64| ((k / kr2(b)) * 20.0).round() / 20.0;
    let names: BTreeSet<String> = m.keys().map(|k| k.0.clone()).collect();
    // geometric mean over curves at (name, b, grid multiplier); all 3 curves required
    let gmc = |name: &str, b: u64, mu: f64| -> Option<f64> {
        let v = m.get(&(name.to_string(), b))?;
        let xs: Vec<f64> = v
            .iter()
            .filter(|t| (mult(t.0, b) - mu).abs() < 1e-9)
            .map(|t| t.2)
            .collect();
        if xs.len() >= 3 {
            Some(gm(&xs))
        } else {
            None
        }
    };
    let mut res = json!({});
    println!("## C3: epsilon at kappa_R2 (R1 primitives, W-walk collection), ratio to os = 0\n");
    println!("| arm | sizes | gm ratio over sizes |");
    println!("|---|---|---|");
    let mut best = (f64::INFINITY, String::new());
    for name in &names {
        if name.contains('+') {
            continue;
        }
        let mut lr = vec![];
        for &b in &bits {
            if let (Some(x), Some(y)) = (gmc(name, b, 1.0), gmc("lpsgebw2", b, 1.0)) {
                lr.push(x.ln() - y.ln());
            }
        }
        if lr.is_empty() {
            continue;
        }
        let r = (lr.iter().sum::<f64>() / lr.len() as f64).exp();
        println!("| {name} | {} | {r:.4} |", lr.len());
        if r < best.0 - 1e-12 {
            best = (r, name.clone());
        }
    }
    println!("\nC3 rule: epsilon* = `{}` (ratio {:.4})\n", best.1, best.0);
    res["c3"] = json!({"arm": best.1, "ratio": best.0});

    println!("## C4: C1 + C2 + C3, kappa re-tuned: quadratic in ln kappa per size, per epsilon\n");
    println!("| arm | log2 q | grid points | kappa* | predicted min cost | at kappa_R2 |");
    println!("|---|---|---|---|---|---|");
    let mut best4 = (f64::INFINITY, String::new(), vec![]);
    for name in &names {
        if !name.ends_with("+C1C2") {
            continue;
        }
        let mut tot = 0.0;
        let mut pts = vec![];
        let mut nb = 0;
        for &b in &bits {
            let Some(v) = m.get(&(name.clone(), b)) else {
                continue;
            };
            let mus: BTreeSet<i64> = v.iter().map(|t| (mult(t.0, b) * 20.0) as i64).collect();
            let (mut x, mut y) = (vec![], vec![]);
            for mu in &mus {
                let mu = *mu as f64 / 20.0;
                if let Some(g) = gmc(name, b, mu) {
                    x.push((mu * kr2(b)).ln());
                    y.push(g.ln());
                }
            }
            if x.len() < 3 {
                continue;
            }
            let (c0, c1, c2) = quad(&x, &y);
            let (xlo, xhi) = (
                x.iter().cloned().fold(f64::INFINITY, f64::min),
                x.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            );
            let mut xs = if c2 > 0.0 { -c1 / (2.0 * c2) } else { f64::NAN };
            if !xs.is_finite() {
                // not convex: grid argmin
                let i = (0..y.len())
                    .min_by(|&a, &b| y[a].partial_cmp(&y[b]).unwrap())
                    .unwrap();
                xs = x[i];
            }
            let xs = xs.clamp(xlo, xhi);
            let ymin = c0 + c1 * xs + c2 * xs * xs;
            tot += ymin;
            nb += 1;
            pts.push((b, xs));
            println!(
                "| {name} | {b} | {} | {:.4} | {:.4e} | {:.4e} |",
                x.len(),
                xs.exp(),
                ymin.exp(),
                gmc(name, b, 1.0).unwrap_or(f64::NAN)
            );
        }
        if nb == bits.len() && tot < best4.0 {
            best4 = (tot, name.clone(), pts.clone());
        }
    }
    let (bx, ly): (Vec<f64>, Vec<f64>) = best4.2.iter().map(|&(b, x)| (b as f64 - 20.0, x)).unzip();
    let beta = ols(&bx, &ly);
    let alpha = ly.iter().sum::<f64>() / ly.len() as f64 - beta * bx.iter().sum::<f64>() / bx.len() as f64;
    println!(
        "\nC4 rule: arm `{}`; ln kappa_C4(b) = {:.4} + ({:.4}) (b - 20)\n",
        best4.1, alpha, beta
    );
    println!("| b | 12 | 14 | 16 | 18 | 20 | 21 | 22 |");
    println!("|---|---|---|---|---|---|---|---|");
    print!("| kappa_C4(b) |");
    for b in [12, 14, 16, 18, 20, 21, 22] {
        print!(" {:.3} |", (alpha + beta * (b as f64 - 20.0)).exp());
    }
    println!();
    print!("| kappa_R2(b) |");
    for b in [12, 14, 16, 18, 20, 21, 22] {
        print!(" {:.3} |", kr2(b));
    }
    println!();
    res["c4"] = json!({"arm": best4.1, "alpha": alpha, "beta": beta, "points": best4.2});
    fs::write(out, serde_json::to_string_pretty(&res).unwrap()).unwrap();
}

/// Relation-phase cost of the W-walk collections (relation-only training jobs at one size):
/// trials to the first check and relation cost (R1 unbatched and C2 batched), per W.
fn wtune(dir: &str) {
    let d = load(dir);
    // (b, W) -> curve -> (trials, R1 relation cost, C2 relation cost)
    let mut m: BTreeMap<(u64, u64), BTreeMap<u64, (f64, f64, f64)>> = BTreeMap::new();
    for ((arm, b, c), r) in &d {
        let w = f(r, "walks") as u64;
        let e = m.entry((*b, w)).or_default().entry(*c).or_insert((0.0, 0.0, 0.0));
        if arm.ends_with(":none") {
            e.0 = f(r, "steps");
            e.1 = f(r, "rel_cost");
            if w == 1 {
                e.2 = e.1;
            }
        } else if arm.ends_with(":none+C2") {
            e.2 = f(r, "rel_cost");
        }
    }
    println!("| log2 q | W | n | gm trials | gm relation cost, R1 (unbatched) | gm relation cost, batched (C2) | C2 / R1 |");
    println!("|---|---|---|---|---|---|---|");
    for ((b, w), v) in &m {
        let t = gm(&v.values().map(|x| x.0).collect::<Vec<_>>());
        let r1 = gm(&v.values().map(|x| x.1).collect::<Vec<_>>());
        let c2 = gm(&v.values().map(|x| x.2).collect::<Vec<_>>());
        println!(
            "| {b} | {w} | {} | {t:.4e} | {r1:.4e} | {c2:.4e} | {:.4} |",
            v.len(),
            c2 / r1
        );
    }
}

/// least-squares quadratic y = a + b x + c x^2; returns (a, b, c)
fn quad(x: &[f64], y: &[f64]) -> (f64, f64, f64) {
    let mut s = [[0.0f64; 4]; 3];
    for (&xi, &yi) in x.iter().zip(y) {
        let p = [1.0, xi, xi * xi];
        for r in 0..3 {
            for c in 0..3 {
                s[r][c] += p[r] * p[c];
            }
            s[r][3] += p[r] * yi;
        }
    }
    solve3(s)
}

fn solve3(mut s: [[f64; 4]; 3]) -> (f64, f64, f64) {
    for c in 0..3 {
        let piv = (c..3)
            .max_by(|&a, &b| s[a][c].abs().partial_cmp(&s[b][c].abs()).unwrap())
            .unwrap();
        s.swap(c, piv);
        for r in 0..3 {
            if r != c {
                let fac = s[r][c] / s[c][c];
                for k in c..4 {
                    s[r][k] -= fac * s[c][k];
                }
            }
        }
    }
    (s[0][3] / s[0][0], s[1][3] / s[1][1], s[2][3] / s[2][2])
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a[1].as_str() {
        "gate" => gate(&a[2], &a[3]),
        "pertrial" => pertrial(&a[2], &a[3], &a[4]),
        "tune" => tune(&a[2], &a[3]),
        "wtune" => wtune(&a[2]),
        "analyze" => analyze(
            &a[2],
            a[3].parse().unwrap(),
            a[4].parse().unwrap(),
            a[5].parse().unwrap(),
            &a[6],
            &a[7..],
        ),
        _ => panic!("usage"),
    }
}

fn analyze(dir: &str, lo: u64, hi: u64, up: u64, out: &str, arms_in: &[String]) {
    let d = load(dir);
    let arms: Vec<String> = arms_in.to_vec();
    let inc = arms[0].clone();
    println!("## Correctness (every record in the run directory)\n");
    println!("| arm | correct / solves | sizes |");
    println!("|---|---|---|");
    let all_arms: BTreeSet<String> = d.keys().map(|k| k.0.clone()).collect();
    for arm in &all_arms {
        let rs: Vec<_> = d.iter().filter(|(k, _)| &k.0 == arm).collect();
        let ok = rs
            .iter()
            .filter(|(_, r)| r["correct"] == json!(true))
            .count();
        let bs: BTreeSet<u64> = rs.iter().map(|(k, _)| k.1).collect();
        println!(
            "| {} | {}/{} | {}..{} |",
            arm,
            ok,
            rs.len(),
            bs.iter().next().unwrap(),
            bs.iter().last().unwrap()
        );
    }
    let good = |arm: &str, b: u64, c: u64| {
        d.get(&(arm.to_string(), b, c))
            .map(|r| r["correct"] == json!(true))
            .unwrap_or(false)
    };
    // curves per size on which every listed arm is correct
    let common = |b: u64| -> Vec<u64> {
        (0..10u64)
            .filter(|&c| arms.iter().all(|a| good(a, b, c)))
            .collect()
    };
    let all_bits: Vec<u64> = (lo..=hi).filter(|&b| common(b).len() >= 3).collect();
    println!(
        "\nsizes with >= 3 common curves: {:?}; curves per size: {:?}",
        all_bits,
        all_bits
            .iter()
            .map(|&b| common(b).len())
            .collect::<Vec<_>>()
    );
    let lnq = |b: u64| -> f64 { f(&d[&(inc.clone(), b, common(b)[0])], "q").ln() };
    let nboot = 2000;
    let mut summary = json!({"arms": arms, "bits": all_bits, "incumbent": inc, "per_arm": {}});
    let mut rng = Rng(0x9e3779b97f4a7c15);
    let base: BTreeMap<u64, Vec<u64>> = all_bits.iter().map(|&b| (b, common(b))).collect();
    let mut plans: Vec<BTreeMap<u64, Vec<u64>>> = vec![base.clone()];
    for _ in 0..nboot {
        plans.push(
            base.iter()
                .map(|(&b, cs)| (b, (0..cs.len()).map(|_| cs[rng.below(cs.len())]).collect()))
                .collect(),
        );
    }
    let val = |arm: &str, b: u64, c: u64, m: &str| -> Option<f64> {
        d.get(&(arm.to_string(), b, c))
            .filter(|r| r["correct"] == json!(true))
            .map(|r| f(r, m))
    };
    let slope = |arm: &str, m: &str, l: u64, h: u64, plan: &BTreeMap<u64, Vec<u64>>| -> f64 {
        let (mut x, mut y) = (vec![], vec![]);
        for (&b, cs) in plan {
            if b < l || b > h {
                continue;
            }
            for &c in cs {
                if let Some(v) = val(arm, b, c, m) {
                    x.push(lnq(b));
                    y.push(v.ln());
                }
            }
        }
        if x.len() < 2 {
            return f64::NAN;
        }
        ols(&x, &y)
    };
    let gmp = |arm: &str, m: &str, b: u64, cs: &[u64]| -> f64 {
        let v: Vec<f64> = cs.iter().filter_map(|&c| val(arm, b, c, m)).collect();
        if v.len() < cs.len() || v.is_empty() {
            f64::NAN
        } else {
            gm(&v)
        }
    };
    let est_ci = |fun: &dyn Fn(&BTreeMap<u64, Vec<u64>>) -> f64| -> (f64, [f64; 2]) {
        let e = fun(&plans[0]);
        let v: Vec<f64> = plans[1..].iter().map(fun).collect();
        (e, ci(v))
    };
    let fmt = |e: (f64, [f64; 2])| format!("{:.3} [{:.3}, {:.3}]", e.0, e.1[0], e.1[1]);
    let fmts = |e: (f64, [f64; 2])| format!("{:+.4} [{:+.4}, {:+.4}]", e.0, e.1[0], e.1[1]);
    let h = *all_bits.last().unwrap();
    let l = *all_bits.first().unwrap();
    let up = up.max(l);

    println!("\n## Exponents, b = {l}..{h} and {up}..{h} (Delta paired against {inc})\n");
    println!("| arm | metric | slope {l}-{h} | Delta | slope {up}-{h} | Delta |");
    println!("|---|---|---|---|---|---|");
    for arm in &arms {
        let mut ent = json!({});
        for m in ["cost", "rel_cost", "la_cost", "cost_cl9", "la_cost_cl9"] {
            let mut cells = vec![arm.clone(), m.to_string()];
            for ll in [l, up] {
                let e = est_ci(&|p| slope(arm, m, ll, h, p));
                let dl = est_ci(&|p| slope(arm, m, ll, h, p) - slope(&inc, m, ll, h, p));
                ent[format!("slope_{m}_{ll}_{h}")] = json!([e.0, e.1[0], e.1[1]]);
                ent[format!("delta_{m}_{ll}_{h}")] = json!([dl.0, dl.1[0], dl.1[1]]);
                cells.push(fmt(e));
                cells.push(if arm == &inc { "-".into() } else { fmts(dl) });
            }
            println!("| {} |", cells.join(" | "));
        }
        summary["per_arm"][arm] = ent;
    }

    println!("\n## Per size (geometric means over the common curves; ratios paired against {inc}, 95% bootstrap CI)\n");
    for arm in &arms {
        let mut rows = vec![];
        for &b in &all_bits {
            let cs = base[&b].clone();
            let g = gmp(arm, "cost", b, &cs);
            let gla = gmp(arm, "la_cost", b, &cs);
            let grel = gmp(arm, "rel_cost", b, &cs);
            let ratio = |m: &str| {
                est_ci(&|p| {
                    let cs = &p[&b];
                    let v: Vec<f64> = cs
                        .iter()
                        .filter_map(|&c| Some(val(arm, b, c, m)? / val(&inc, b, c, m)?))
                        .collect();
                    if v.len() < cs.len() {
                        f64::NAN
                    } else {
                        gm(&v)
                    }
                })
            };
            let rw = ratio("cost");
            let rla = ratio("la_cost");
            let rrel = ratio("rel_cost");
            let rw9 = ratio("cost_cl9");
            let mean = |k: &str| {
                cs.iter().filter_map(|&c| val(arm, b, c, k)).sum::<f64>() / cs.len() as f64
            };
            // product split: matrix-vector products (last check) vs the rest of the LA products
            let mv_prod = cs
                .iter()
                .map(|&c| val(arm, b, c, "matvecs").unwrap() * val(arm, b, c, "la_nnz").unwrap())
                .sum::<f64>()
                / cs.len() as f64;
            rows.push(json!({"bits": b, "curves": cs, "gm_cost": g, "gm_la": gla, "gm_rel": grel,
                "ratio": [rw.0, rw.1[0], rw.1[1]], "ratio_la": [rla.0, rla.1[0], rla.1[1]],
                "ratio_rel": [rrel.0, rrel.1[0], rrel.1[1]], "ratio_cl9": [rw9.0, rw9.1[0], rw9.1[1]],
                "la_share": gla / g, "la_cols": mean("la_cols"), "la_nnz": mean("la_nnz"), "nfb": mean("nfb"),
                "checks": mean("checks"), "la_wall": mean("la_wall"), "matvecs": mean("matvecs"),
                "mv_products": mv_prod, "la_products": mean("la_cost"), "kappa": mean("kappa"), "cx": mean("cx"),
                "costs": cs.iter().map(|&c| val(arm, b, c, "cost").unwrap()).collect::<Vec<_>>(),
                "q": lnq(b).exp()}));
        }
        summary["per_arm"][arm]["sizes"] = json!(rows);
        println!("### {arm}\n");
        println!("| log2 q | n | kappa | cx | cost | LA share | ratio whole [95%] | ratio relation | ratio LA [95%] | ratio c_ell=9 | solver cols | nnz/col | matvecs / cols | mat-vec share of LA products | checks |");
        println!("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|");
        for r in &rows {
            let tr = |k: &str| {
                fmt((
                    r[k][0].as_f64().unwrap(),
                    [r[k][1].as_f64().unwrap(), r[k][2].as_f64().unwrap()],
                ))
            };
            println!(
                "| {} | {} | {:.3} | {} | {:.3e} | {:.2} | {} | {:.3} | {} | {:.3} | {:.0} | {:.1} | {:.3} | {:.3} | {:.2} |",
                r["bits"], r["curves"].as_array().unwrap().len(), f(r, "kappa"), f(r, "cx"), f(r, "gm_cost"), f(r, "la_share"),
                tr("ratio"), r["ratio_rel"][0].as_f64().unwrap(), tr("ratio_la"), r["ratio_cl9"][0].as_f64().unwrap(),
                f(r, "la_cols"), f(r, "la_nnz") / f(r, "la_cols"), f(r, "matvecs") / f(r, "la_cols"),
                f(r, "mv_products") / f(r, "la_products"), f(r, "checks")
            );
        }
        println!();
    }

    println!("## Local exponents ({l}..{h}), their trend per bit, and the paired trend difference against {inc}\n");
    for m in ["cost", "rel_cost", "la_cost"] {
        println!("### {m}\n");
        print!("| interval |");
        for a in &arms {
            print!(" {a} |");
        }
        println!();
        print!("|---|");
        for _ in &arms {
            print!("---|");
        }
        println!();
        let loc_of = |a: &str, b0: u64, b1: u64, p: &BTreeMap<u64, Vec<u64>>| {
            (gmp(a, m, b1, &p[&b1]).ln() - gmp(a, m, b0, &p[&b0]).ln()) / (lnq(b1) - lnq(b0))
        };
        for w in all_bits.windows(2) {
            let (b0, b1) = (w[0], w[1]);
            print!("| {b0}-{b1} |");
            for a in &arms {
                let e = est_ci(&|p| loc_of(a, b0, b1, p));
                if let Some(v) = summary["per_arm"][a][format!("local_{m}")].as_array_mut() {
                    v.push(json!([b0, e.0, e.1[0], e.1[1]]))
                }
                if summary["per_arm"][a][format!("local_{m}")].is_null() {
                    summary["per_arm"][a][format!("local_{m}")] =
                        json!([[b0, e.0, e.1[0], e.1[1]]]);
                }
                print!(" {} |", fmt(e));
            }
            println!();
        }
        let trend = |a: &str, p: &BTreeMap<u64, Vec<u64>>| -> f64 {
            let w: Vec<(f64, f64)> = all_bits
                .windows(2)
                .map(|w| ((w[0] + w[1]) as f64 / 2.0, loc_of(a, w[0], w[1], p)))
                .collect();
            ols(
                &w.iter().map(|x| x.0).collect::<Vec<_>>(),
                &w.iter().map(|x| x.1).collect::<Vec<_>>(),
            )
        };
        print!("| trend per bit |");
        for a in &arms {
            let e = est_ci(&|p| trend(a, p));
            summary["per_arm"][a][format!("trend_{m}")] = json!([e.0, e.1[0], e.1[1]]);
            print!(" {} |", fmts(e));
        }
        println!();
        print!("| trend difference |");
        for a in &arms {
            if a == &inc {
                print!(" - |");
                continue;
            }
            let e = est_ci(&|p| trend(a, p) - trend(&inc, p));
            summary["per_arm"][a][format!("trend_diff_{m}")] = json!([e.0, e.1[0], e.1[1]]);
            print!(" {} |", fmts(e));
        }
        println!("\n");
    }
    fs::write(out, serde_json::to_string_pretty(&summary).unwrap()).unwrap();
}
