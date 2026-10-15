//! Analysis for IC improvement round 2, item 06 (kappa schedule and block Wiedemann for genus-3 IC),
//! adapted from `ic_improve2_01_analyze.rs`.
//!
//! `gate <run_dir> <item01_run_dir>`: `main:lpsgewd` must equal item 01's `main:lpsgewd` records.
//! `gatetrain <train_dir> <item01_tune2_dir>`: training `lpsgewd@8` at b = 20 must equal item 01's `tune2`.
//! `train <dir>`: the frozen training rules of PROTOCOL.md (m*, scalar or block, cx(b), kappa(b) schedule).
//! `analyze <run_dir> <lo> <hi> <upper_lo> <out.json> <arm>...`: exponents (whole, relation, LA) with
//! paired Deltas against the first arm (the incumbent), local exponents, their trend and the paired trend
//! difference, paired ratios per size, block-Wiedemann product split. 2000 bootstrap resamples of curves
//! within each size, paired across arms; a size uses the curves on which every listed arm is correct.
//!
//! ```sh
//! cargo run --release --example ic_improve2_06_analyze -- train research/ic_improve_20261015/06-genus3-kappa-bw/results/train
//! ```

use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;

fn read_jsonl(dir: &str) -> Vec<Value> {
    let mut out = vec![];
    let mut names: Vec<_> = fs::read_dir(dir).unwrap().flatten().map(|e| e.path()).collect();
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
    v[k].as_f64().unwrap_or_else(|| v[k].as_u64().map(|x| x as f64).unwrap_or(f64::NAN))
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
    println!("| log2 q | curve | k equal | cost (item 01) | cost (item 06) | equal |");
    println!("|---|---|---|---|---|---|");
    for ((arm, b, c), r) in &d {
        if arm != "main:lpsgewd" {
            continue;
        }
        let Some(p) = o.get(&(arm.clone(), *b, *c)) else { continue };
        let eq = r["cost"] == p["cost"] && r["k"] == p["k"] && r["la_counters"] == p["la_counters"];
        if eq {
            ok += 1
        } else {
            bad += 1
        }
        println!("| {} | {} | {} | {} | {} | {} |", b, c, r["k"] == p["k"], p["cost"], r["cost"], if eq { "yes" } else { "NO" });
    }
    println!("\ngate 0: {ok} identical, {bad} different");
}


/// (bits, kappa, la-name) -> (curve, record)
fn train_table(dir: &str) -> (BTreeMap<(u64, String, String), Vec<(u64, Value)>>, usize) {
    let mut m: BTreeMap<(u64, String, String), Vec<(u64, Value)>> = BTreeMap::new();
    let mut fails = 0;
    for rec in read_jsonl(dir) {
        let b = rec["bits"].as_u64().unwrap();
        let c = rec["curve"].as_u64().unwrap();
        for r in rec["results"].as_array().unwrap() {
            if r["correct"] != json!(true) {
                fails += 1;
                continue;
            }
            let kv = format!("{:.2}", f(r, "kappa"));
            m.entry((b, kv, r["la"].as_str().unwrap().to_string())).or_default().push((c, r.clone()));
        }
    }
    (m, fails)
}

fn gm_of(v: &[(u64, Value)], k: &str) -> f64 {
    gm(&v.iter().map(|(_, r)| f(r, k)).collect::<Vec<_>>())
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
    for c in 0..3 {
        let piv = (c..3).max_by(|&a, &b| s[a][c].abs().partial_cmp(&s[b][c].abs()).unwrap()).unwrap();
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

fn train(dir: &str, out: Option<&String>) {
    let (m, fails) = train_table(dir);
    let sizes: BTreeSet<u64> = m.keys().map(|k| k.0).collect();
    let las: BTreeSet<String> = m.keys().map(|k| k.2.clone()).collect();
    println!("## Training table (geometric-mean whole cost; LA cost in brackets; curves 5-7)\n");
    println!("incorrect or missing k: {fails}\n");
    for &b in &sizes {
        let ks: BTreeSet<String> = m.keys().filter(|k| k.0 == b).map(|k| k.1.clone()).collect();
        print!("### b = {b}\n\n| kappa |");
        for l in &las {
            print!(" {l} |");
        }
        println!();
        print!("|---|");
        for _ in &las {
            print!("---|");
        }
        println!();
        for k in &ks {
            print!("| {k} |");
            for l in &las {
                match m.get(&(b, k.clone(), l.clone())) {
                    Some(v) => print!(" {:.4e} [{:.3e}] (n={}) |", gm_of(v, "cost"), gm_of(v, "la_cost"), v.len()),
                    None => print!(" - |"),
                }
            }
            println!();
        }
        println!();
    }
    let complete = |b: u64, k: &str, l: &str| m.get(&(b, k.to_string(), l.to_string())).map(|v| v.len() >= 3).unwrap_or(false);
    // rule 1: m*
    let top = *sizes.iter().filter(|&&b| complete(b, "2.00", "lpsgebw2@8") && complete(b, "2.00", "lpsgebw3@8")).max().unwrap();
    let la_gm = |b: u64, k: &str, l: &str| gm_of(&m[&(b, k.to_string(), l.to_string())], "la_cost");
    let w_gm = |b: u64, k: &str, l: &str| gm_of(&m[&(b, k.to_string(), l.to_string())], "cost");
    let mut mstar = if la_gm(top, "2.00", "lpsgebw2@8") <= la_gm(top, "2.00", "lpsgebw3@8") { 2 } else { 3 };
    println!("## Frozen rules\n");
    println!("1. m*: at b = {top}, kappa 2: LA gm `lpsgebw2@8` {:.4e}, `lpsgebw3@8` {:.4e}", la_gm(top, "2.00", "lpsgebw2@8"), la_gm(top, "2.00", "lpsgebw3@8"));
    if complete(18, "2.00", "lpsgebw4@8") {
        let l4 = la_gm(18, "2.00", "lpsgebw4@8");
        let (l2, l3) = (la_gm(18, "2.00", "lpsgebw2@8"), la_gm(18, "2.00", "lpsgebw3@8"));
        println!("   b = 18: m=2 {l2:.4e}, m=3 {l3:.4e}, m=4 {l4:.4e}");
        if l4 < l2 && l4 < l3 {
            mstar = 4;
        }
    }
    println!("   **m* = {mstar}**");
    let bw = format!("lpsgebw{mstar}@8");
    let ws = w_gm(top, "2.00", "lpsgewd@8");
    let wb = w_gm(top, "2.00", &bw);
    println!("2. scalar or block at b = {top}, kappa 2: whole gm `lpsgewd@8` {ws:.4e}, `{bw}` {wb:.4e} -> **{}**", if wb < ws { "block" } else { "scalar" });
    // rules 3, 4
    let mut pts = vec![];
    let mut curves_json = vec![];
    println!("3./4. per size: best kappa per cx (grid), chosen cx, quadratic kappa*:\n");
    println!("| b | cx | grid-best kappa | gm cost there | quadratic minimiser |");
    println!("|---|---|---|---|---|");
    for &b in &sizes {
        let ks: Vec<String> = m.keys().filter(|k| k.0 == b && k.2 == "lpsgewd@8").map(|k| k.1.clone()).collect();
        let mut best: Option<(f64, String, String)> = None;
        for cxs in ["8", "16", "32"] {
            let l = format!("lpsgewd@{cxs}");
            let kk: Vec<&String> = ks.iter().filter(|k| complete(b, k, &l)).collect();
            if kk.is_empty() {
                continue;
            }
            let (bk, bc) = kk.iter().map(|k| (k.to_string(), w_gm(b, k, &l))).min_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();
            if best.as_ref().map(|x| bc < x.0).unwrap_or(true) {
                best = Some((bc, cxs.to_string(), bk));
            }
        }
        let (bc, cxs, bk) = best.unwrap();
        let l = format!("lpsgewd@{cxs}");
        let kk: Vec<&String> = ks.iter().filter(|k| complete(b, k, &l)).collect();
        let x: Vec<f64> = kk.iter().map(|k| k.parse::<f64>().unwrap().ln()).collect();
        let y: Vec<f64> = kk.iter().map(|k| w_gm(b, k, &l).ln()).collect();
        let (lo, hi) = (x.iter().cloned().fold(f64::MAX, f64::min), x.iter().cloned().fold(f64::MIN, f64::max));
        let (_, qb, qc) = quad(&x, &y);
        let kstar = if qc > 0.0 { (-qb / (2.0 * qc)).clamp(lo, hi) } else { bk.parse::<f64>().unwrap().ln() };
        println!("| {b} | {cxs} | {bk} | {bc:.4e} | {:.3} |", kstar.exp());
        pts.push((b as f64, kstar));
        curves_json.push(json!({"bits": b, "cx": cxs, "kappa_star": kstar.exp(),
            "grid": kk.iter().map(|k| json!([k.parse::<f64>().unwrap(), w_gm(b, k, &l)])).collect::<Vec<_>>()}));
    }
    let xs: Vec<f64> = pts.iter().map(|p| p.0 - 20.0).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.1).collect();
    let beta = ols(&xs, &ys);
    let alpha = ys.iter().sum::<f64>() / ys.len() as f64 - beta * xs.iter().sum::<f64>() / xs.len() as f64;
    println!("\n**schedule: ln kappa(b) = {alpha:.4} + ({beta:.4}) (b - 20)**, i.e. kappa(b) = {:.4} * {:.4}^(b - 20)\n", alpha.exp(), beta.exp());
    print!("| b |");
    for b in 12..=22 {
        print!(" {b} |");
    }
    println!();
    print!("|---|");
    for _ in 12..=22 {
        print!("---|");
    }
    println!();
    print!("| kappa(b) |");
    for b in 12..=22 {
        print!(" {:.3} |", (alpha + beta * (b as f64 - 20.0)).exp());
    }
    println!();
    if let Some(o) = out {
        let j = json!({"alpha": alpha, "beta": beta, "mstar": mstar, "block": wb < ws, "sizes": curves_json});
        fs::write(o, serde_json::to_string_pretty(&j).unwrap()).unwrap();
    }
}

fn gatetrain(tr: &str, old: &str) {
    let (m, _) = train_table(tr);
    let (o, _) = train_table(old);
    let (mut ok, mut bad) = (0, 0);
    for ((b, k, l), v) in &m {
        if l != "lpsgewd@8" {
            continue;
        }
        let Some(ov) = o.get(&(*b, k.clone(), "lpsgewd".to_string())) else { continue };
        for (c, r) in v {
            if let Some((_, p)) = ov.iter().find(|(cc, _)| cc == c) {
                let eq = r["cost"] == p["cost"] && r["k"] == p["k"] && r["la_counters"] == p["la_counters"];
                println!("| {b} | {k} | {c} | {} | {} | {} |", p["cost"], r["cost"], if eq { "yes" } else { "NO" });
                if eq { ok += 1 } else { bad += 1 }
            }
        }
    }
    println!("\ntraining gate: {ok} identical, {bad} different");
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a[1].as_str() {
        "gate" => gate(&a[2], &a[3]),
        "gatetrain" => gatetrain(&a[2], &a[3]),
        "train" => train(&a[2], a.get(3)),
        "analyze" => analyze(&a[2], a[3].parse().unwrap(), a[4].parse().unwrap(), a[5].parse().unwrap(), &a[6], &a[7..]),
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
        let ok = rs.iter().filter(|(_, r)| r["correct"] == json!(true)).count();
        let bs: BTreeSet<u64> = rs.iter().map(|(k, _)| k.1).collect();
        println!("| {} | {}/{} | {}..{} |", arm, ok, rs.len(), bs.iter().next().unwrap(), bs.iter().last().unwrap());
    }
    let good = |arm: &str, b: u64, c: u64| d.get(&(arm.to_string(), b, c)).map(|r| r["correct"] == json!(true)).unwrap_or(false);
    // curves per size on which every listed arm is correct
    let common = |b: u64| -> Vec<u64> { (0..10u64).filter(|&c| arms.iter().all(|a| good(a, b, c))).collect() };
    let all_bits: Vec<u64> = (lo..=hi).filter(|&b| common(b).len() >= 3).collect();
    println!("\nsizes with >= 3 common curves: {:?}; curves per size: {:?}", all_bits, all_bits.iter().map(|&b| common(b).len()).collect::<Vec<_>>());
    let lnq = |b: u64| -> f64 { f(&d[&(inc.clone(), b, common(b)[0])], "q").ln() };
    let nboot = 2000;
    let mut summary = json!({"arms": arms, "bits": all_bits, "incumbent": inc, "per_arm": {}});
    let mut rng = Rng(0x9e3779b97f4a7c15);
    let base: BTreeMap<u64, Vec<u64>> = all_bits.iter().map(|&b| (b, common(b))).collect();
    let mut plans: Vec<BTreeMap<u64, Vec<u64>>> = vec![base.clone()];
    for _ in 0..nboot {
        plans.push(base.iter().map(|(&b, cs)| (b, (0..cs.len()).map(|_| cs[rng.below(cs.len())]).collect())).collect());
    }
    let val = |arm: &str, b: u64, c: u64, m: &str| -> Option<f64> {
        d.get(&(arm.to_string(), b, c)).filter(|r| r["correct"] == json!(true)).map(|r| f(r, m))
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
        if v.len() < cs.len() || v.is_empty() { f64::NAN } else { gm(&v) }
    };
    let est_ci = |fun: &dyn Fn(&BTreeMap<u64, Vec<u64>>) -> f64| -> (f64, [f64; 2]) {
        let e = fun(&plans[0]);
        let v: Vec<f64> = plans[1..].iter().map(|p| fun(p)).collect();
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
            let ratio = |m: &str| est_ci(&|p| {
                let cs = &p[&b];
                let v: Vec<f64> = cs.iter().filter_map(|&c| Some(val(arm, b, c, m)? / val(&inc, b, c, m)?)).collect();
                if v.len() < cs.len() { f64::NAN } else { gm(&v) }
            });
            let rw = ratio("cost");
            let rla = ratio("la_cost");
            let rrel = ratio("rel_cost");
            let rw9 = ratio("cost_cl9");
            let mean = |k: &str| cs.iter().filter_map(|&c| val(arm, b, c, k)).sum::<f64>() / cs.len() as f64;
            // product split: matrix-vector products (last check) vs the rest of the LA products
            let mv_prod = cs.iter().map(|&c| val(arm, b, c, "matvecs").unwrap() * val(arm, b, c, "la_nnz").unwrap()).sum::<f64>() / cs.len() as f64;
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
            let tr = |k: &str| fmt((r[k][0].as_f64().unwrap(), [r[k][1].as_f64().unwrap(), r[k][2].as_f64().unwrap()]));
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
        let loc_of = |a: &str, b0: u64, b1: u64, p: &BTreeMap<u64, Vec<u64>>| (gmp(a, m, b1, &p[&b1]).ln() - gmp(a, m, b0, &p[&b0]).ln()) / (lnq(b1) - lnq(b0));
        for w in all_bits.windows(2) {
            let (b0, b1) = (w[0], w[1]);
            print!("| {b0}-{b1} |");
            for a in &arms {
                let e = est_ci(&|p| loc_of(a, b0, b1, p));
                summary["per_arm"][a][format!("local_{m}")].as_array_mut().map(|v| v.push(json!([b0, e.0, e.1[0], e.1[1]])));
                if summary["per_arm"][a][format!("local_{m}")].is_null() {
                    summary["per_arm"][a][format!("local_{m}")] = json!([[b0, e.0, e.1[0], e.1[1]]]);
                }
                print!(" {} |", fmt(e));
            }
            println!();
        }
        let trend = |a: &str, p: &BTreeMap<u64, Vec<u64>>| -> f64 {
            let w: Vec<(f64, f64)> = all_bits.windows(2).map(|w| ((w[0] + w[1]) as f64 / 2.0, loc_of(a, w[0], w[1], p))).collect();
            ols(&w.iter().map(|x| x.0).collect::<Vec<_>>(), &w.iter().map(|x| x.1).collect::<Vec<_>>())
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
