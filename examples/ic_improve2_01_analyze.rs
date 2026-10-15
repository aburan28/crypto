//! Analysis for IC improvement round 2, item 01 (linear algebra of genus-3 IC-2LP).
//!
//! `gate <run_dir> <item04_run_dir>`: `main:inc` must equal item 04's `ic2lp` records (k, cost).
//! `tune1 <dir>`: geometric-mean whole cost per (arm, cx, size) on the training curves.
//! `tune2 <dir>`: geometric-mean whole cost per (arm, kappa, size) on the training curves.
//! `analyze <run_dir> <lo> <hi> <upper_lo> <out.json>`: per-arm exponents of the whole cost and of
//! the linear-algebra cost (least squares of ln cost on ln q, every solve a point, 95% percentile
//! bootstrap over curves within each size, 2000 resamples, paired across arms), exponent
//! differences against `main:inc`, local exponents and their trend, paired ratios per size.
//! Prints Markdown tables and writes a JSON summary for the figure script.
//!
//! ```sh
//! cargo run --release --example ic_improve2_01_analyze -- analyze \
//!     research/ic_improve_20261015/01-genus3-linalg/results/run 10 20 14 summary.json
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
            let arm = r["arm"].as_str().unwrap().to_string();
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

fn gate(run: &str, old: &str) {
    let d = load(run);
    let o = load(old);
    let (mut ok, mut bad) = (0, 0);
    println!("| log2 q | curve | k equal | cost (item 04) | cost (main:inc) | equal |");
    println!("|---|---|---|---|---|---|");
    for ((arm, b, c), r) in &d {
        if arm != "main:inc" {
            continue;
        }
        let Some(p) = o.get(&("ic2lp".to_string(), *b, *c)) else {
            continue;
        };
        let eq = r["cost"] == p["cost"] && r["k"] == p["k"] && r["la_counters"] == p["la_counters"];
        if eq {
            ok += 1
        } else {
            bad += 1
        }
        println!(
            "| {} | {} | {} | {} | {} | {} |",
            b,
            c,
            r["k"] == p["k"],
            p["cost"],
            r["cost"],
            if eq { "yes" } else { "NO" }
        );
    }
    println!("\ngate 0: {ok} identical, {bad} different");
}

fn tune(dir: &str, key: &str) {
    let recs = read_jsonl(dir);
    // (arm, key value, bits) -> costs
    let mut m: BTreeMap<(String, String, u64), Vec<f64>> = BTreeMap::new();
    let mut fails = 0;
    for rec in &recs {
        let b = rec["bits"].as_u64().unwrap();
        for r in rec["results"].as_array().unwrap() {
            if r["correct"] != json!(true) {
                fails += 1;
            }
            let kv = format!("{}", r[key]);
            m.entry((r["la"].as_str().unwrap().to_string(), kv, b))
                .or_default()
                .push(f(r, "cost"));
        }
    }
    let bits: BTreeSet<u64> = m.keys().map(|k| k.2).collect();
    print!("| arm | {key} |");
    for b in &bits {
        print!(" b={b} |");
    }
    println!();
    print!("|---|---|");
    for _ in &bits {
        print!("---|");
    }
    println!();
    let rows: BTreeSet<(String, String)> = m.keys().map(|k| (k.0.clone(), k.1.clone())).collect();
    for (a, kv) in rows {
        print!("| {a} | {kv} |");
        for b in &bits {
            match m.get(&(a.clone(), kv.clone(), *b)) {
                Some(v) => print!(" {:.4e} (n={}) |", gm(v), v.len()),
                None => print!(" - |"),
            }
        }
        println!();
    }
    println!("\nincorrect results: {fails}");
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a[1].as_str() {
        "gate" => gate(&a[2], &a[3]),
        "tune1" => tune(&a[2], "cx"),
        "tune2" => tune(&a[2], "kappa"),
        "analyze" => analyze(
            &a[2],
            a[3].parse().unwrap(),
            a[4].parse().unwrap(),
            a[5].parse().unwrap(),
            &a[6],
        ),
        _ => panic!("usage"),
    }
}

fn analyze(dir: &str, lo: u64, hi: u64, up: u64, out: &str) {
    let d = load(dir);
    let inc = "main:inc".to_string();
    let arms: Vec<String> = {
        let s: BTreeSet<String> = d.keys().map(|k| k.0.clone()).collect();
        let mut v: Vec<String> = s.into_iter().collect();
        v.sort_by_key(|x| (x != &inc, x.clone()));
        v
    };
    // correctness
    println!("## Correctness\n");
    println!("| arm | correct / solves | sizes |");
    println!("|---|---|---|");
    for arm in &arms {
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
    // curves per size, common to inc and the arm
    let curves = |arm: &str, b: u64| -> Vec<u64> {
        (0..10u64)
            .filter(|c| {
                d.get(&(arm.to_string(), b, *c))
                    .map(|r| r["correct"] == json!(true))
                    .unwrap_or(false)
                    && d.get(&(inc.clone(), b, *c))
                        .map(|r| r["correct"] == json!(true))
                        .unwrap_or(false)
            })
            .collect()
    };
    let all_bits: Vec<u64> = (lo..=hi).filter(|b| !curves(&inc, *b).is_empty()).collect();
    let lnq = |b: u64| -> f64 { f(&d[&(inc.clone(), b, curves(&inc, b)[0])], "q").ln() };
    let metrics = ["cost", "la_cost", "cost_cl9", "la_cost_cl9", "rel_cost"];
    let nboot = 2000;
    let mut summary = json!({"arms": arms, "bits": all_bits, "per_arm": {}});
    // resample plan: per size, a vector of curve indices (into the inc curve list)
    let mut rng = Rng(0x9e3779b97f4a7c15);
    let base: BTreeMap<u64, Vec<u64>> = all_bits.iter().map(|&b| (b, curves(&inc, b))).collect();
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
    // slope over [l, h] for arm/metric under plan
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
        let nb: BTreeSet<u64> = plan
            .keys()
            .filter(|&&b| b >= l && b <= h)
            .copied()
            .collect();
        if nb.len() < 2 || x.len() < 2 {
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
        let v: Vec<f64> = plans[1..].iter().map(|p| fun(p)).collect();
        (e, ci(v))
    };
    let fmt = |e: (f64, [f64; 2])| format!("{:.3} [{:.3}, {:.3}]", e.0, e.1[0], e.1[1]);
    let fmtr = |e: (f64, [f64; 2])| format!("{:.3} [{:.3}, {:.3}]", e.0, e.1[0], e.1[1]);
    let arm_bits = |arm: &str| -> Vec<u64> {
        all_bits
            .iter()
            .copied()
            .filter(|&b| !curves(arm, b).is_empty() && curves(arm, b).len() == base[&b].len())
            .collect()
    };

    println!("\n## Exponents (whole cost and linear algebra), b = {lo}..{hi} and {up}..{hi}\n");
    println!("| arm | sizes | whole {lo}-{hi} | Delta vs inc | whole {up}-{hi} | Delta vs inc | LA {lo}-{hi} | Delta vs inc | LA {up}-{hi} | whole c_ell=9 {lo}-{hi} | Delta vs inc | LA c_ell=9 {lo}-{hi} | relation phase {lo}-{hi} |");
    println!("|---|---|---|---|---|---|---|---|---|---|---|---|---|");
    for arm in &arms {
        let bs = arm_bits(arm);
        if bs.len() < 4 {
            continue;
        }
        let (l, h) = (*bs.first().unwrap(), *bs.last().unwrap());
        let h = h.min(hi);
        let mut ent = json!({});
        let mut cells = vec![arm.clone(), format!("{l}-{h}")];
        for (m, ll) in [
            ("cost", l),
            ("cost", up.max(l)),
            ("la_cost", l),
            ("la_cost", up.max(l)),
            ("cost_cl9", l),
            ("la_cost_cl9", l),
            ("rel_cost", l),
        ] {
            let e = est_ci(&|p| slope(arm, m, ll, h, p));
            let dl = est_ci(&|p| slope(arm, m, ll, h, p) - slope(&inc, m, ll, h, p));
            ent[format!("slope_{m}_{ll}_{h}")] = json!([e.0, e.1[0], e.1[1]]);
            ent[format!("delta_{m}_{ll}_{h}")] = json!([dl.0, dl.1[0], dl.1[1]]);
            cells.push(fmt(e));
            if !(m == "la_cost" && ll != l) && m != "la_cost_cl9" && m != "rel_cost" {
                cells.push(if arm == &inc { "-".into() } else { fmt(dl) });
            }
        }
        println!("| {} |", cells.join(" | "));
        summary["per_arm"][arm] = ent;
    }

    // per-size table
    println!("\n## Per size (geometric means over curves; ratios paired against main:inc, 95% bootstrap CI)\n");
    for arm in &arms {
        let mut rows = vec![];
        for &b in &all_bits {
            let cs = curves(arm, b);
            if cs.len() != base[&b].len() {
                continue;
            }
            let g = gmp(arm, "cost", b, &cs);
            let gla = gmp(arm, "la_cost", b, &cs);
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
            let rw9 = ratio("cost_cl9");
            let mean = |k: &str| {
                cs.iter().filter_map(|&c| val(arm, b, c, k)).sum::<f64>() / cs.len() as f64
            };
            let share = gla / g;
            rows.push(json!({"bits": b, "gm_cost": g, "gm_la": gla, "ratio": [rw.0, rw.1[0], rw.1[1]],
                "ratio_la": [rla.0, rla.1[0], rla.1[1]], "ratio_cl9": [rw9.0, rw9.1[0], rw9.1[1]],
                "la_share": share, "la_cols": mean("la_cols"), "la_nnz": mean("la_nnz"), "la_in_cols": mean("la_in_cols"),
                "la_in_nnz": mean("la_in_nnz"), "checks": mean("checks"), "la_wall": mean("la_wall"),
                "costs": cs.iter().map(|&c| val(arm, b, c, "cost").unwrap()).collect::<Vec<_>>(),
                "la_costs": cs.iter().map(|&c| val(arm, b, c, "la_cost").unwrap()).collect::<Vec<_>>(),
                "q": lnq(b).exp()}));
        }
        summary["per_arm"][arm]["sizes"] = json!(rows);
        println!("### {arm}\n");
        println!("| log2 q | cost | LA cost | LA share | ratio whole [95%] | ratio LA [95%] | ratio whole c_ell=9 | solver cols | solver nnz/col | input cols | input nnz/col | checks | LA wall s |");
        println!("|---|---|---|---|---|---|---|---|---|---|---|---|---|");
        for r in &rows {
            let tr = |k: &str| {
                fmtr((
                    r[k][0].as_f64().unwrap(),
                    [r[k][1].as_f64().unwrap(), r[k][2].as_f64().unwrap()],
                ))
            };
            println!(
                "| {} | {:.3e} | {:.3e} | {:.2} | {} | {} | {} | {:.0} | {:.1} | {:.0} | {:.1} | {:.2} | {:.1} |",
                r["bits"], f(r, "gm_cost"), f(r, "gm_la"), f(r, "la_share"), tr("ratio"), tr("ratio_la"), tr("ratio_cl9"),
                f(r, "la_cols"), f(r, "la_nnz") / f(r, "la_cols"), f(r, "la_in_cols"), f(r, "la_in_nnz") / f(r, "la_in_cols"),
                f(r, "checks"), f(r, "la_wall")
            );
        }
        println!();
    }

    // local exponents
    println!("## Local exponents (whole cost; LA cost in the second table)\n");
    for m in ["cost", "la_cost"] {
        print!("| interval |");
        let shown: Vec<&String> = arms.iter().filter(|a| arm_bits(a).len() >= 4).collect();
        for a in &shown {
            print!(" {a} |");
        }
        println!();
        print!("|---|");
        for _ in &shown {
            print!("---|");
        }
        println!();
        let mut loc: BTreeMap<String, Vec<(u64, (f64, [f64; 2]))>> = BTreeMap::new();
        for w in all_bits.windows(2) {
            let (b0, b1) = (w[0], w[1]);
            print!("| {b0}-{b1} |");
            for a in &shown {
                if !arm_bits(a).contains(&b0) || !arm_bits(a).contains(&b1) {
                    print!(" - |");
                    continue;
                }
                let e = est_ci(&|p| {
                    (gmp(a, m, b1, &p[&b1]).ln() - gmp(a, m, b0, &p[&b0]).ln())
                        / (lnq(b1) - lnq(b0))
                });
                loc.entry(a.to_string()).or_default().push((b0, e));
                print!(" {} |", fmt(e));
            }
            println!();
        }
        // trend
        print!("| trend per bit |");
        for a in &shown {
            let bs: Vec<u64> = loc
                .get(*a)
                .map(|v| v.iter().map(|x| x.0).collect())
                .unwrap_or_default();
            if bs.len() < 3 {
                print!(" - |");
                continue;
            }
            let e = est_ci(&|p| {
                let ys: Vec<f64> = bs
                    .iter()
                    .map(|&b0| {
                        (gmp(a, m, b0 + 1, &p[&(b0 + 1)]).ln() - gmp(a, m, b0, &p[&b0]).ln())
                            / (lnq(b0 + 1) - lnq(b0))
                    })
                    .collect();
                let xs: Vec<f64> = bs.iter().map(|&b| b as f64).collect();
                ols(&xs, &ys)
            });
            summary["per_arm"][*a][format!("trend_{m}")] = json!([e.0, e.1[0], e.1[1]]);
            print!(" {:+.4} [{:+.4}, {:+.4}] |", e.0, e.1[0], e.1[1]);
        }
        println!();
        for (a, v) in loc {
            summary["per_arm"][&a][format!("local_{m}")] = json!(v
                .iter()
                .map(|(b, e)| json!([b, e.0, e.1[0], e.1[1]]))
                .collect::<Vec<_>>());
        }
        println!();
    }
    fs::write(out, serde_json::to_string_pretty(&summary).unwrap()).unwrap();
}
