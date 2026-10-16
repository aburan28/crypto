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
            if let Some((base, rest)) = arm.split_once('@') {
                let suffix = rest.find('+').map(|i| &rest[i..]).unwrap_or("");
                arm = format!("{base}{suffix}");
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

fn gate(run: &str, old: &str) {
    let d = load(run);
    let o = load(old);
    let (mut ok, mut bad) = (0, 0);
    println!("| log2 q | curve | k equal | cost (item 06) | cost (this item) | relation cost equal | LA counters equal | identical |");
    println!("|---|---|---|---|---|---|---|---|");
    for ((arm, b, c), r) in &d {
        if arm != "sched:lpsgebw2" {
            continue;
        }
        let Some(p) = o.get(&(arm.clone(), *b, *c)) else {
            continue;
        };
        let eq = r["cost"] == p["cost"]
            && r["k"] == p["k"]
            && r["rel_cost"] == p["rel_cost"]
            && r["la_counters"] == p["la_counters"];
        if eq {
            ok += 1
        } else {
            bad += 1
        }
        println!(
            "| {} | {} | {} | {} | {} | {} | {} | {} |",
            b,
            c,
            r["k"] == p["k"],
            p["cost"],
            r["cost"],
            r["rel_cost"] == p["rel_cost"],
            r["la_counters"] == p["la_counters"],
            if eq { "yes" } else { "NO" }
        );
    }
    println!("\ngate 0: {ok} identical, {bad} different");
}

/// (bits, kappa "x.xx", la-name) -> (curve, record)
type Table = BTreeMap<(u64, String, String), Vec<(u64, Value)>>;

fn train_table(dir: &str) -> Table {
    let mut m: Table = BTreeMap::new();
    for rec in read_jsonl(dir) {
        let b = rec["bits"].as_u64().unwrap();
        let c = rec["curve"].as_u64().unwrap();
        for r in rec["results"].as_array().unwrap() {
            let mut r = r.clone();
            r["q"] = rec["q"].clone();
            let kv = format!("{:.2}", f(&r, "kappa"));
            let name = r["arm"]
                .as_str()
                .unwrap()
                .split(':')
                .nth(1)
                .unwrap()
                .to_string();
            m.entry((b, kv, name)).or_default().push((c, r));
        }
    }
    m
}

fn rec_of<'a>(m: &'a Table, b: u64, k: &str, la: &str, c: u64) -> Option<&'a Value> {
    m.get(&(b, k.to_string(), la.to_string()))
        .and_then(|v| v.iter().find(|(cc, _)| *cc == c).map(|(_, r)| r))
}

fn gatetrain(tr: &str, old: &str) {
    let m = train_table(tr);
    let o = train_table(old);
    let (mut ok, mut bad, mut na) = (0, 0, 0);
    println!("| log2 q | kappa | curve | relation cost (item 06, `lpsgebw2@8`) | relation-only rerun | equal |");
    println!("|---|---|---|---|---|---|");
    for ((b, k, la), v) in &m {
        if la != "none" {
            continue;
        }
        for (c, r) in v {
            let Some(p) = rec_of(&o, *b, k, "lpsgebw2@8", *c) else {
                continue;
            };
            if p["checks"] != json!(1) {
                na += 1;
                println!(
                    "| {b} | {k} | {c} | {} | {} | n/a (checks > 1) |",
                    p["rel_cost"], r["rel_cost"]
                );
                continue;
            }
            let eq = p["rel_cost"] == r["rel_cost"] && p["steps"] == r["steps"];
            println!(
                "| {b} | {k} | {c} | {} | {} | {} |",
                p["rel_cost"],
                r["rel_cost"],
                if eq { "yes" } else { "NO" }
            );
            if eq {
                ok += 1
            } else {
                bad += 1
            }
        }
    }
    println!("\ntraining gate: {ok} identical, {bad} different, {na} not comparable");
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

/// OLS y = c0 + a x1 + beta x2; returns (a, beta, rss)
fn ols2(x1: &[f64], x2: &[f64], y: &[f64]) -> (f64, f64, f64) {
    let mut s = [[0.0f64; 4]; 3];
    for i in 0..y.len() {
        let p = [1.0, x1[i], x2[i]];
        for r in 0..3 {
            for c in 0..3 {
                s[r][c] += p[r] * p[c];
            }
            s[r][3] += p[r] * y[i];
        }
    }
    let (c0, a, b) = solve3(s);
    let rss = (0..y.len())
        .map(|i| (y[i] - c0 - a * x1[i] - b * x2[i]).powi(2))
        .sum();
    (a, b, rss)
}

fn rss1(x: &[f64], y: &[f64]) -> f64 {
    let a = ols(x, y);
    let n = x.len() as f64;
    let c = y.iter().sum::<f64>() / n - a * x.iter().sum::<f64>() / n;
    x.iter()
        .zip(y)
        .map(|(xi, yi)| (yi - c - a * xi).powi(2))
        .sum()
}

/// Per-record decomposition of the relation phase. Units per trial unless stated.
fn decomp(r: &Value) -> BTreeMap<&'static str, f64> {
    let q = f(r, "q");
    let d = &r["diag"];
    let t = f(d, "trials");
    let mut m = BTreeMap::new();
    let rel = f(r, "rel_cost");
    let tr = units(&r["transfer_counters"], q);
    let u = |k: &str| units(&d[k], q);
    m.insert("q", q);
    m.insert("kappa", f(r, "kappa"));
    m.insert("nfb", f(r, "nfb"));
    m.insert("trials", t);
    m.insert("t_kappa_q43", t * f(r, "kappa") / q.powf(4.0 / 3.0));
    m.insert("rel", rel);
    m.insert("rel_per_trial", rel / t);
    m.insert("walk", u("walk_add") / t);
    m.insert("split", u("split") / t);
    m.insert("roots", u("roots") / t);
    m.insert("cand_walk", u("cand_walk_add") / t);
    m.insert("cand_split", u("cand_split") / t);
    m.insert("cand_roots", u("cand_roots") / t);
    m.insert("fb_share", u("fb") / rel);
    m.insert("setup_share", u("setup") / rel);
    // the remainder: the y-coordinate evaluations at the roots (bookkeeping), per trial
    m.insert(
        "other",
        (rel - tr - u("fb") - u("setup") - u("walk_add") - u("split") - u("roots")) / t,
    );
    let walk_i = d["walk_add"][1].as_u64().unwrap() as f64;
    let walk_ci = d["cand_walk_add"][1].as_u64().unwrap() as f64;
    let adds = d["walk_add"][4].as_u64().unwrap() as f64;
    m.insert("inv_per_add", walk_i / adds);
    m.insert("cand_inv_per_add", walk_ci / adds);
    m.insert(
        "mq_per_add",
        d["walk_add"][0].as_u64().unwrap() as f64 / adds,
    );
    m.insert(
        "cand_mq_per_add",
        d["cand_walk_add"][0].as_u64().unwrap() as f64 / adds,
    );
    let sp = f(d, "splits_true");
    m.insert("p_split", sp / t);
    let x = (q - 1.0) / 2.0;
    let fr = f(r, "nfb") / x;
    m.insert("f", fr);
    let lp2 = f(d, "lp2");
    m.insert("p2lp_given_split", lp2 / f(d, "split_new"));
    m.insert("p2lp_theory", 3.0 * fr * (1.0 - fr) * (1.0 - fr));
    m.insert("p1lp_given_split", f(d, "lp1") / f(d, "split_new"));
    m.insert("p1lp_theory", 3.0 * fr * fr * (1.0 - fr));
    m.insert("pfull_given_split", f(d, "full") / f(d, "split_new"));
    m.insert("p3lp_given_split", f(d, "lp3") / f(d, "split_new"));
    let v = x - f(r, "nfb");
    m.insert("mean_degree", 2.0 * lp2 / v);
    m.insert("cycles_over_nfb", f(d, "cycles") / f(r, "nfb"));
    // trials an independent-trial model needs to reach the measured 2LP edge count
    m.insert(
        "t_model_over_t",
        lp2 / ((1.0 / 6.0) * 3.0 * fr * (1.0 - fr) * (1.0 - fr)) / t,
    );
    m
}

fn diag(dir: &str, out: &str) {
    let m = train_table(dir);
    // (b, curve) -> (incumbent record, R1 record) at kappa 2
    let mut rows: BTreeMap<u64, Vec<(u64, Value, Value)>> = BTreeMap::new();
    for ((b, k, la), v) in &m {
        if k != "2.00" || la != "none" {
            continue;
        }
        for (c, r) in v {
            let r1 = rec_of(&m, *b, k, "none+R1", *c).unwrap().clone();
            rows.entry(*b).or_default().push((*c, r.clone(), r1));
        }
    }
    let keys = [
        "trials",
        "t_kappa_q43",
        "rel_per_trial",
        "walk",
        "split",
        "roots",
        "cand_walk",
        "cand_split",
        "cand_roots",
        "inv_per_add",
        "cand_inv_per_add",
        "mq_per_add",
        "cand_mq_per_add",
        "fb_share",
        "setup_share",
        "other",
        "p_split",
        "f",
        "p2lp_given_split",
        "p2lp_theory",
        "p1lp_given_split",
        "p1lp_theory",
        "pfull_given_split",
        "p3lp_given_split",
        "mean_degree",
        "cycles_over_nfb",
        "t_model_over_t",
    ];
    println!("## Step-1 decomposition (training curves 5-7, kappa = 2, relation phase only; geometric means over curves)\n");
    println!("Per-trial columns are counted units per trial (`M_q + ceil(log2 q) I_q`).\n");
    let mut js = vec![];
    print!("| log2 q | n |");
    for k in keys {
        print!(" {k} |");
    }
    println!(" R1/inc relation |");
    print!("|---|---|");
    for _ in keys {
        print!("---|");
    }
    println!("---|");
    for (b, v) in &rows {
        let ds: Vec<_> = v.iter().map(|(_, r, _)| decomp(r)).collect();
        let ratio = gm(&v
            .iter()
            .map(|(_, r, r1)| f(r1, "rel_cost") / f(r, "rel_cost"))
            .collect::<Vec<_>>());
        print!("| {b} | {} |", v.len());
        let mut e = json!({"bits": b, "rel_r1_over_inc": ratio});
        for k in keys {
            let xs: Vec<f64> = ds.iter().map(|d| d[k]).collect();
            let val = gm(&xs);
            e[k] = json!(val);
            if val.abs() >= 100.0 {
                print!(" {val:.4e} |");
            } else {
                print!(" {val:.4} |");
            }
        }
        println!(" {ratio:.4} |");
        js.push(e);
    }
    // exponent fits, every solve a point, bootstrap curves within size
    let bits: Vec<u64> = rows.keys().cloned().collect();
    let mut rng = Rng(0x9e3779b97f4a7c15);
    let get = |b: u64, c: u64, key: &str| -> f64 {
        let (_, r, r1) = rows[&b].iter().find(|x| x.0 == c).unwrap();
        match key {
            "rel" => f(r, "rel_cost"),
            "rel_r1" => f(r1, "rel_cost"),
            "trials" => f(&r["diag"], "trials"),
            "per_trial" => f(r, "rel_cost") / f(&r["diag"], "trials"),
            "per_trial_r1" => f(r1, "rel_cost") / f(&r["diag"], "trials"),
            _ => unreachable!(),
        }
    };
    let curves: BTreeMap<u64, Vec<u64>> = rows
        .iter()
        .map(|(b, v)| (*b, v.iter().map(|x| x.0).collect()))
        .collect();
    let qof = |b: u64| f(&rows[&b][0].1, "q");
    let mut plans = vec![curves.clone()];
    for _ in 0..2000 {
        plans.push(
            curves
                .iter()
                .map(|(&b, cs)| (b, (0..cs.len()).map(|_| cs[rng.below(cs.len())]).collect()))
                .collect(),
        );
    }
    println!("\n## Exponent fits over {}..{} (OLS of ln y on ln q, every solve a point; 95% bootstrap CI over curves within size)\n", bits[0], bits[bits.len() - 1]);
    println!(
        "| quantity | q^c: c | q^a (log q)^beta: a | beta | RSS ratio (two-term / one-term) |"
    );
    println!("|---|---|---|---|---|");
    let mut fits = json!({});
    for key in ["trials", "per_trial", "per_trial_r1", "rel", "rel_r1"] {
        let fit = |p: &BTreeMap<u64, Vec<u64>>| {
            let (mut x, mut x2, mut y) = (vec![], vec![], vec![]);
            for (&b, cs) in p {
                for &c in cs {
                    let lq = qof(b).ln();
                    x.push(lq);
                    x2.push(lq.ln());
                    y.push(get(b, c, key).ln());
                }
            }
            let c1 = ols(&x, &y);
            let (a, be, rss2) = ols2(&x, &x2, &y);
            (c1, a, be, rss2 / rss1(&x, &y))
        };
        let e = fit(&plans[0]);
        let bs: Vec<_> = plans[1..].iter().map(fit).collect();
        let cc = ci(bs.iter().map(|x| x.0).collect());
        let ca = ci(bs.iter().map(|x| x.1).collect());
        let cb = ci(bs.iter().map(|x| x.2).collect());
        println!(
            "| {key} | {:.4} [{:.4}, {:.4}] | {:.3} [{:.3}, {:.3}] | {:.2} [{:.2}, {:.2}] | {:.3} |",
            e.0, cc[0], cc[1], e.1, ca[0], ca[1], e.2, cb[0], cb[1], e.3
        );
        fits[key] = json!({"c": [e.0, cc[0], cc[1]], "a": [e.1, ca[0], ca[1]], "beta": [e.2, cb[0], cb[1]], "rss_ratio": e.3});
    }
    // per-trial cost = A + B log2 q (geometric means per size)
    let (mut xb, mut ya, mut yr) = (vec![], vec![], vec![]);
    for r in &js {
        xb.push(r["bits"].as_u64().unwrap() as f64);
        ya.push(f(r, "rel_per_trial"));
        yr.push(f(r, "rel_per_trial") * f(r, "rel_r1_over_inc"));
    }
    let lin = |y: &[f64]| {
        let bb = ols(&xb, y);
        let n = xb.len() as f64;
        (
            y.iter().sum::<f64>() / n - bb * xb.iter().sum::<f64>() / n,
            bb,
        )
    };
    let (a0, b0) = lin(&ya);
    let (a1, b1) = lin(&yr);
    println!("\nPer-trial cost (geometric mean per size) fitted as A + B log2 q: incumbent A = {a0:.1}, B = {b0:.2}; R1 A = {a1:.1}, B = {b1:.2}.");
    println!("Implied log-factor contribution to the local exponent at log2 q = 21, B / ((A + 21 B) ln 2): incumbent {:.4}, R1 {:.4}.", b0 / ((a0 + 21.0 * b0) * 2f64.ln()), b1 / ((a1 + 21.0 * b1) * 2f64.ln()));
    // local exponents of trials and of the relation cost (gm per size)
    println!("\n| interval | trials local exponent | relation local exponent (incumbent) | R1 | per-trial (incumbent) | per-trial (R1) |");
    println!("|---|---|---|---|---|---|");
    let gmk = |b: u64, key: &str| {
        gm(&curves[&b]
            .iter()
            .map(|&c| get(b, c, key))
            .collect::<Vec<_>>())
    };
    let mut locs = vec![];
    for w in bits.windows(2) {
        let dl = qof(w[1]).ln() - qof(w[0]).ln();
        let le = |k: &str| (gmk(w[1], k).ln() - gmk(w[0], k).ln()) / dl;
        println!(
            "| {}-{} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} |",
            w[0],
            w[1],
            le("trials"),
            le("rel"),
            le("rel_r1"),
            le("per_trial"),
            le("per_trial_r1")
        );
        locs.push(json!([
            w[0],
            le("trials"),
            le("rel"),
            le("rel_r1"),
            le("per_trial"),
            le("per_trial_r1")
        ]));
    }
    // kappa dependence of rho at b = 16, 18, 20
    println!("\n## R1 / incumbent relation-cost ratio rho against kappa (geometric mean over curves 5-7)\n");
    println!("| log2 q | kappa | rho | trials kappa / q^(4/3) | n |");
    println!("|---|---|---|---|---|");
    let mut rhos = vec![];
    for ((b, k, la), v) in &m {
        if la != "none" {
            continue;
        }
        let rr: Vec<f64> = v
            .iter()
            .map(|(c, r)| {
                f(rec_of(&m, *b, k, "none+R1", *c).unwrap(), "rel_cost") / f(r, "rel_cost")
            })
            .collect();
        let tk: Vec<f64> = v.iter().map(|(_, r)| decomp(r)["t_kappa_q43"]).collect();
        println!(
            "| {b} | {k} | {:.4} | {:.4} | {} |",
            gm(&rr),
            gm(&tk),
            v.len()
        );
        rhos.push(json!([b, k, gm(&rr), gm(&tk)]));
    }
    let j = json!({"sizes": js, "fits": fits, "per_trial_lin": {"inc": [a0, b0], "r1": [a1, b1]}, "local": locs, "rho": rhos});
    fs::write(out, serde_json::to_string_pretty(&j).unwrap()).unwrap();
}

fn r2rule(tr: &str, old: &str, out: &str) {
    let m = train_table(tr);
    let o = train_table(old);
    println!(
        "## Frozen R2 rule (PROTOCOL.md): predicted R1 whole cost on item 06's training grid\n"
    );
    let mut pts = vec![];
    let mut sizes_json = vec![];
    for b in [16u64, 18, 20, 21] {
        // measured kappas (relation-only) at this size
        let meas: Vec<String> = m
            .keys()
            .filter(|k| k.0 == b && k.2 == "none")
            .map(|k| k.1.clone())
            .collect();
        let grid: Vec<String> = o
            .keys()
            .filter(|k| k.0 == b && k.2 == "lpsgebw2@8")
            .map(|k| k.1.clone())
            .collect();
        println!("### b = {b} (rho measured at kappa {:?})\n", meas);
        println!(
            "| kappa | rho used (kappa') | predicted R1 whole (gm) | incumbent whole (gm) | n |"
        );
        println!("|---|---|---|---|---|");
        let (mut x, mut y) = (vec![], vec![]);
        let mut gridj = vec![];
        for k in &grid {
            let kf: f64 = k.parse().unwrap();
            let kp = meas
                .iter()
                .min_by(|a, bb| {
                    let da = (a.parse::<f64>().unwrap().ln() - kf.ln()).abs();
                    let db = (bb.parse::<f64>().unwrap().ln() - kf.ln()).abs();
                    da.partial_cmp(&db).unwrap()
                })
                .unwrap();
            let mut pred = vec![];
            let mut incw = vec![];
            let mut rhos = vec![];
            for c in [5u64, 6, 7] {
                let Some(p) = rec_of(&o, b, k, "lpsgebw2@8", c) else {
                    continue;
                };
                let ri = rec_of(&m, b, kp, "none", c).unwrap();
                let r1 = rec_of(&m, b, kp, "none+R1", c).unwrap();
                let rho = f(r1, "rel_cost") / f(ri, "rel_cost");
                // gate: the item-06 relation cost is used only if the rerun at this kappa equals it
                if kp == k && ri["rel_cost"] != p["rel_cost"] {
                    panic!(
                        "gate: b={b} kappa={k} curve={c}: rerun relation cost differs from item 06"
                    );
                }
                rhos.push(rho);
                pred.push(rho * f(p, "rel_cost") + f(p, "la_cost"));
                incw.push(f(p, "cost"));
            }
            if pred.len() < 3 {
                continue;
            }
            println!(
                "| {k} | {:.4} ({kp}) | {:.4e} | {:.4e} | {} |",
                gm(&rhos),
                gm(&pred),
                gm(&incw),
                pred.len()
            );
            x.push(kf.ln());
            y.push(gm(&pred).ln());
            gridj.push(json!([kf, gm(&pred), gm(&incw), gm(&rhos)]));
        }
        let (lo, hi) = (
            x.iter().cloned().fold(f64::MAX, f64::min),
            x.iter().cloned().fold(f64::MIN, f64::max),
        );
        let (_, qb, qc) = quad(&x, &y);
        let argmin = x[y
            .iter()
            .enumerate()
            .min_by(|a, bb| a.1.partial_cmp(bb.1).unwrap())
            .unwrap()
            .0];
        let kstar = if qc > 0.0 {
            (-qb / (2.0 * qc)).clamp(lo, hi)
        } else {
            argmin
        };
        println!(
            "\nkappa*_R1({b}) = {:.4} (grid argmin {:.2}{})\n",
            kstar.exp(),
            argmin.exp(),
            if qc > 0.0 && ((-qb / (2.0 * qc)) < lo || (-qb / (2.0 * qc)) > hi) {
                ", clamped"
            } else {
                ""
            }
        );
        pts.push((b as f64, kstar));
        sizes_json.push(json!({"bits": b, "kappa_star": kstar.exp(), "grid": gridj}));
    }
    let xs: Vec<f64> = pts.iter().map(|p| p.0 - 20.0).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.1).collect();
    let beta = ols(&xs, &ys);
    let alpha =
        ys.iter().sum::<f64>() / ys.len() as f64 - beta * xs.iter().sum::<f64>() / xs.len() as f64;
    println!("**R2 schedule: ln kappa(b) = {alpha:.4} + ({beta:.4}) (b - 20)**, i.e. kappa(b) = {:.4} * {:.4}^(b - 20)\n", alpha.exp(), beta.exp());
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
    for (name, a, be) in [
        ("kappa_R2(b)", alpha, beta),
        ("item 06 kappa(b)", 0.5533, -0.0806),
    ] {
        print!("| {name} |");
        for b in 12..=22 {
            print!(" {:.3} |", (a + be * (b as f64 - 20.0)).exp());
        }
        println!();
    }
    let j = json!({"alpha": alpha, "beta": beta, "sizes": sizes_json});
    fs::write(out, serde_json::to_string_pretty(&j).unwrap()).unwrap();
}

fn acct(dir: &str) {
    let d = load(dir);
    println!("| arm | log2 q | curve | relation cost, R1 driving the walk (`+R1s`) | shadow-accounted `+R1` | equal | k correct |");
    println!("|---|---|---|---|---|---|---|");
    for ((arm, b, c), r) in &d {
        let Some(base) = arm.strip_suffix("+R1s") else {
            continue;
        };
        let Some(p) = d.get(&(format!("{base}+R1"), *b, *c)) else {
            continue;
        };
        println!(
            "| {arm} | {b} | {c} | {} | {} | {} | {} |",
            r["rel_cost"],
            p["rel_cost"],
            if r["rel_cost"] == p["rel_cost"] {
                "yes"
            } else {
                "NO"
            },
            r["correct"]
        );
    }
}

fn logfit(dir: &str, arm: &str, metric: &str, lo: u64, hi: u64) {
    let d = load(dir);
    let mut by: BTreeMap<u64, Vec<(f64, f64)>> = BTreeMap::new();
    for ((a, b, _), r) in &d {
        if a == arm && *b >= lo && *b <= hi && r["correct"] == json!(true) {
            by.entry(*b)
                .or_default()
                .push((f(r, "q").ln(), f(r, metric).ln()));
        }
    }
    let mut rng = Rng(0x9e3779b97f4a7c15);
    let fit = |p: &BTreeMap<u64, Vec<(f64, f64)>>| {
        let (mut x, mut x2, mut y) = (vec![], vec![], vec![]);
        for v in p.values() {
            for &(lq, ly) in v {
                x.push(lq);
                x2.push(lq.ln());
                y.push(ly);
            }
        }
        let c = ols(&x, &y);
        let (a, be, r2) = ols2(&x, &x2, &y);
        (c, a, be, r2 / rss1(&x, &y))
    };
    let e = fit(&by);
    let bs: Vec<_> = (0..2000)
        .map(|_| {
            let p: BTreeMap<u64, Vec<(f64, f64)>> = by
                .iter()
                .map(|(b, v)| (*b, (0..v.len()).map(|_| v[rng.below(v.len())]).collect()))
                .collect();
            fit(&p)
        })
        .collect();
    let cc = ci(bs.iter().map(|x| x.0).collect());
    let ca = ci(bs.iter().map(|x| x.1).collect());
    let cb = ci(bs.iter().map(|x| x.2).collect());
    println!(
        "| {arm} | {metric} | {lo}-{hi} | {:.4} [{:.4}, {:.4}] | {:.3} [{:.3}, {:.3}] | {:.2} [{:.2}, {:.2}] | {:.3} |",
        e.0, cc[0], cc[1], e.1, ca[0], ca[1], e.2, cb[0], cb[1], e.3
    );
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a[1].as_str() {
        "gate" => gate(&a[2], &a[3]),
        "gatetrain" => gatetrain(&a[2], &a[3]),
        "diag" => diag(&a[2], &a[3]),
        "r2rule" => r2rule(&a[2], &a[3], &a[4]),
        "acct" => acct(&a[2]),
        "logfit" => logfit(
            &a[2],
            &a[3],
            &a[4],
            a[5].parse().unwrap(),
            a[6].parse().unwrap(),
        ),
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
