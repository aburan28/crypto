//! Analysis for IC improvement item 04 (genus-3 IC after the cubic GHS transfer).
//!
//! `xcheck <python_run_dir> <rust_dir>`: compares the Rust port's `--mode python`
//! records with item 7's `ghs3.py` records field by field.
//!
//! `analyze <run_dir> <common_lo> <common_hi> <upper_lo> <out.json>`: per-arm
//! exponents (least squares of ln cost on ln q, 95% percentile bootstrap over
//! curves within each size, 2000 resamples), paired exponent differences, local
//! exponents, paired cost ratios per size, the local-exponent trend, and rho
//! reference values. Writes a JSON summary and prints Markdown tables.
//!
//! ```sh
//! cargo run --release --example ic_improve_04_analyze -- analyze \
//!     research/ic_improve_20261014/04-genus3-2lp/results/run-20261014 10 18 14 summary.json
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

fn xcheck(pydir: &str, rsdir: &str) {
    let py = read_jsonl(pydir);
    let rs = read_jsonl(rsdir);
    let key = |v: &Value| (v["bits"].as_u64().unwrap(), v["curve"].as_u64().unwrap());
    let pym: BTreeMap<_, _> = py
        .iter()
        .filter(|v| v["results"][0]["arm"] != "rho")
        .map(|v| (key(v), v))
        .collect();
    let inst_fields = [
        "q",
        "N",
        "ell",
        "c",
        "h",
        "a",
        "b",
        "r0",
        "Kmod",
        "build_tries",
        "checks",
    ];
    let res_fields = [
        "arm",
        "k",
        "correct",
        "nfb",
        "rows",
        "lp_vertices",
        "checks",
        "counters",
        "setup_counters",
        "steps",
        "cand_split",
        "full",
        "lp1",
        "lp2",
        "lp3",
        "restarts",
        "cycles",
        "weier",
        "transfer_counters",
        "total_counters",
        "cost",
        "cost_cl9",
    ];
    let (mut n_ok, mut n_bad) = (0, 0);
    println!(
        "| log2 q | curve | arm | k (both) | steps | cost (Python) | cost (Rust) | fields equal |"
    );
    println!("|---|---|---|---|---|---|---|---|");
    for r in &rs {
        let kk = key(r);
        let Some(p) = pym.get(&kk) else {
            println!("missing python record for {:?}", kk);
            n_bad += 1;
            continue;
        };
        let mut bad = vec![];
        for f in inst_fields {
            if r[f] != p[f] {
                bad.push(f.to_string());
            }
        }
        let rr = r["results"].as_array().unwrap();
        let pr = p["results"].as_array().unwrap();
        if rr.len() != pr.len() {
            bad.push("n_results".into());
        }
        for (a, b) in rr.iter().zip(pr.iter()) {
            let mut badr = vec![];
            for f in res_fields {
                if a[f] != b[f] {
                    badr.push(f);
                }
            }
            let ok = bad.is_empty() && badr.is_empty();
            if ok {
                n_ok += 1;
            } else {
                n_bad += 1;
            }
            println!(
                "| {} | {} | {} | {} | {} | {} | {} | {} |",
                kk.0,
                kk.1,
                a["arm"].as_str().unwrap(),
                a["k"],
                a["steps"],
                b["cost"],
                a["cost"],
                if ok {
                    "all".to_string()
                } else {
                    format!("DIFF {:?} {:?}", bad, badr)
                }
            );
        }
    }
    println!("\nxcheck: {} solves identical, {} differing", n_ok, n_bad);
}

// ------------------------------------------------------------- statistics
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

fn ols(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    for (a, b) in x.iter().zip(y.iter()) {
        sxy += (a - mx) * (b - my);
        sxx += (a - mx) * (a - mx);
    }
    sxy / sxx
}

fn pct(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i = ((v.len() - 1) as f64 * q).round() as usize;
    v[i]
}

fn ci(mut v: Vec<f64>) -> [f64; 2] {
    let lo = pct(&mut v, 0.025);
    let hi = pct(&mut v, 0.975);
    [lo, hi]
}

/// data: arm -> bits -> curve -> (ln q, ln cost)
type Data = BTreeMap<String, BTreeMap<u64, BTreeMap<u64, (f64, f64)>>>;

fn curves_common(d: &Data, arms: &[&str], b: u64) -> Vec<u64> {
    let mut s: Option<BTreeSet<u64>> = None;
    for a in arms {
        let cs: BTreeSet<u64> = d
            .get(*a)
            .and_then(|m| m.get(&b))
            .map(|m| m.keys().cloned().collect())
            .unwrap_or_default();
        s = Some(match s {
            None => cs,
            Some(x) => x.intersection(&cs).cloned().collect(),
        });
    }
    s.unwrap_or_default().into_iter().collect()
}

/// slopes for `arms` over `bits`, on a per-size curve resample `sel` (indices into common curves).
fn slopes(d: &Data, arms: &[&str], bits: &[u64], sel: &BTreeMap<u64, Vec<u64>>) -> Vec<f64> {
    arms.iter()
        .map(|a| {
            let (mut x, mut y) = (vec![], vec![]);
            for b in bits {
                for c in &sel[b] {
                    let (lx, ly) = d[*a][b][c];
                    x.push(lx);
                    y.push(ly);
                }
            }
            ols(&x, &y)
        })
        .collect()
}

fn resample(cs: &BTreeMap<u64, Vec<u64>>, rng: &mut Rng) -> BTreeMap<u64, Vec<u64>> {
    cs.iter()
        .map(|(b, v)| (*b, (0..v.len()).map(|_| v[rng.below(v.len())]).collect()))
        .collect()
}

const NB: usize = 2000;

/// exponent of each arm and pairwise differences arms[i] - arms[0], over bits, paired.
fn exp_block(d: &Data, arms: &[&str], bits: &[u64], seed: u64) -> Value {
    let mut cs = BTreeMap::new();
    for b in bits {
        cs.insert(*b, curves_common(d, arms, *b));
    }
    if cs.values().any(|v| v.is_empty()) {
        return json!(null);
    }
    let point = slopes(d, arms, bits, &cs);
    let mut rng = Rng(seed);
    let mut boots: Vec<Vec<f64>> = vec![vec![]; arms.len()];
    let mut dboots: Vec<Vec<f64>> = vec![vec![]; arms.len()];
    for _ in 0..NB {
        let s = resample(&cs, &mut rng);
        let sl = slopes(d, arms, bits, &s);
        for i in 0..arms.len() {
            boots[i].push(sl[i]);
            dboots[i].push(sl[i] - sl[0]);
        }
    }
    let mut out = serde_json::Map::new();
    for (i, a) in arms.iter().enumerate() {
        out.insert(
            a.to_string(),
            json!({"slope": point[i], "ci": ci(boots[i].clone()),
                   "diff_vs_first": point[i] - point[0], "diff_ci": ci(dboots[i].clone())}),
        );
    }
    json!({"bits": bits, "n_curves": cs.values().map(|v| v.len()).collect::<Vec<_>>(), "arms": out})
}

fn gmean(d: &Data, a: &str, b: u64, cs: &[u64]) -> f64 {
    (cs.iter().map(|c| d[a][&b][c].1).sum::<f64>() / cs.len() as f64).exp()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args[1] == "tune" {
        // geometric-mean cost per (arm, kappa, bits) on the training curves
        let mut m: BTreeMap<(String, String, u64), Vec<f64>> = BTreeMap::new();
        for r in read_jsonl(&args[2]) {
            for x in r["results"].as_array().unwrap() {
                assert!(x["correct"].as_bool().unwrap(), "incorrect training solve");
                m.entry((
                    x["arm"].as_str().unwrap().into(),
                    format!("{}", x["kappa"]),
                    r["bits"].as_u64().unwrap(),
                ))
                .or_default()
                .push(x["cost"].as_f64().unwrap().ln());
            }
        }
        println!("| arm | kappa | log2 q | n | geometric-mean cost |");
        println!("|---|---|---|---|---|");
        for ((a, k, b), v) in &m {
            println!(
                "| {} | {} | {} | {} | {:.4e} |",
                a,
                k,
                b,
                v.len(),
                (v.iter().sum::<f64>() / v.len() as f64).exp()
            );
        }
        return;
    }
    if args[1] == "xcheck" {
        xcheck(&args[2], &args[3]);
        return;
    }
    let dir = &args[2];
    let clo: u64 = args[3].parse().unwrap();
    let chi: u64 = args[4].parse().unwrap();
    let ulo: u64 = args[5].parse().unwrap();
    let outp = &args[6];
    let recs = read_jsonl(dir);
    let mut d: Data = BTreeMap::new();
    let mut d9: Data = BTreeMap::new();
    let mut correct: BTreeMap<(String, u64), (u64, u64)> = BTreeMap::new();
    let mut rho: BTreeMap<u64, Vec<(f64, f64, f64, u64)>> = BTreeMap::new(); // cost, model, steps/sqrt(ell), curve
    let mut extra: BTreeMap<(String, u64), Vec<Value>> = BTreeMap::new();
    let mut inst: BTreeMap<(u64, u64), Value> = BTreeMap::new();
    for r in &recs {
        let b = r["bits"].as_u64().unwrap();
        let c = r["curve"].as_u64().unwrap();
        let q = r["q"].as_f64().unwrap();
        inst.insert(
            (b, c),
            json!({"N": r["N"], "ell": r["ell"], "c": r["c"], "checks": r["checks"]}),
        );
        for x in r["results"].as_array().unwrap() {
            let arm = x["arm"].as_str().unwrap().to_string();
            let e = correct.entry((arm.clone(), b)).or_insert((0, 0));
            e.1 += 1;
            if x["correct"].as_bool().unwrap() {
                e.0 += 1;
            }
            if arm == "rho" {
                rho.entry(b).or_default().push((
                    x["cost"].as_f64().unwrap(),
                    x["rho_model_neg"].as_f64().unwrap(),
                    x["steps_over_sqrt_ell"].as_f64().unwrap(),
                    c,
                ));
                continue;
            }
            let cost = x["cost"].as_f64().unwrap();
            if x["correct"].as_bool().unwrap()
                && (arm == "ic2lp" || arm == "basic" || arm == "ic1lp")
            {
                // phase split (derived): linear algebra = counted mults mod ell; the rest = relation phase
                let la = x["counters"][2].as_f64().unwrap();
                d.entry(format!("{}_la", arm))
                    .or_default()
                    .entry(b)
                    .or_default()
                    .insert(c, (q.ln(), la.ln()));
                d.entry(format!("{}_rel", arm))
                    .or_default()
                    .entry(b)
                    .or_default()
                    .insert(c, (q.ln(), (cost - la).ln()));
            }
            if x["correct"].as_bool().unwrap() {
                d.entry(arm.clone())
                    .or_default()
                    .entry(b)
                    .or_default()
                    .insert(c, (q.ln(), cost.ln()));
                d9.entry(arm.clone())
                    .or_default()
                    .entry(b)
                    .or_default()
                    .insert(c, (q.ln(), x["cost_cl9"].as_f64().unwrap().ln()));
            }
            extra.entry((arm.clone(), b)).or_default().push(json!({
                "curve": c, "cost": cost, "steps": x["steps"], "nfb": x["nfb"],
                "la_share": x["counters"][2].as_f64().unwrap() / cost,
                "lp_vertices": x["lp_vertices"], "rows": x["rows"], "la_nnz": x["la_nnz"],
                "la_ncols": x["la_ncols"], "wall": x["wall_total"], "checks": x["checks"],
            }));
        }
    }
    let arms_all: Vec<String> = d.keys().cloned().collect();
    let common: Vec<u64> = (clo..=chi).collect();
    let mut summary = serde_json::Map::new();

    // 1. exponents over the common range: candidate arms vs basic (paired)
    let mut cand: Vec<&str> = vec!["basic"];
    for a in ["ic2lp", "ic1lp", "ic2lpt", "aa"] {
        if arms_all.iter().any(|x| x == a) {
            cand.push(a);
        }
    }
    summary.insert("common_range".into(), exp_block(&d, &cand, &common, 1));
    summary.insert("common_range_cl9".into(), exp_block(&d9, &cand, &common, 2));

    // 2. each arm over its own full range, and over the upper range
    let mut own = serde_json::Map::new();
    let mut upper = serde_json::Map::new();
    for a in &arms_all {
        let bits: Vec<u64> = d[a]
            .keys()
            .cloned()
            .filter(|b| *b >= clo && d[a][b].len() >= 3)
            .collect();
        if bits.len() >= 3 {
            own.insert(a.clone(), exp_block(&d, &[a.as_str()], &bits, 3));
        }
        let ub: Vec<u64> = bits.iter().cloned().filter(|b| *b >= ulo).collect();
        if ub.len() >= 3 {
            upper.insert(a.clone(), exp_block(&d, &[a.as_str()], &ub, 4));
        }
    }
    summary.insert("own_range".into(), Value::Object(own));
    summary.insert("upper_range".into(), Value::Object(upper));
    // A/A and 2LP vs 1LP over the 2LP range
    let lp_bits: Vec<u64> = d
        .get("ic2lp")
        .map(|m| m.keys().cloned().filter(|b| *b >= clo).collect())
        .unwrap_or_default();
    let mut lparms = vec!["ic2lp"];
    for a in ["aa", "ic1lp", "ic2lpt"] {
        if arms_all.iter().any(|x| x == a) {
            lparms.push(a);
        }
    }
    summary.insert(
        "lp_range_vs_ic2lp".into(),
        exp_block(&d, &lparms, &lp_bits, 5),
    );

    // 3. local exponents per arm, with bootstrap
    let mut local = serde_json::Map::new();
    let mut rng = Rng(7);
    for a in &arms_all {
        let bits: Vec<u64> = d[a].keys().cloned().filter(|b| *b >= clo).collect();
        let mut rows = vec![];
        let mut pts_for_trend: Vec<(f64, Vec<f64>, f64)> = vec![];
        for w in bits.windows(2) {
            let (b0, b1) = (w[0], w[1]);
            let c0: Vec<u64> = d[a][&b0].keys().cloned().collect();
            let c1: Vec<u64> = d[a][&b1].keys().cloned().collect();
            let lq0 = d[a][&b0][&c0[0]].0;
            let lq1 = d[a][&b1][&c1[0]].0;
            let est = |s0: &[u64], s1: &[u64]| {
                (gmean(&d, a, b1, s1).ln() - gmean(&d, a, b0, s0).ln()) / (lq1 - lq0)
            };
            let pt = est(&c0, &c1);
            let mut bs = vec![];
            for _ in 0..NB {
                let s0: Vec<u64> = (0..c0.len()).map(|_| c0[rng.below(c0.len())]).collect();
                let s1: Vec<u64> = (0..c1.len()).map(|_| c1[rng.below(c1.len())]).collect();
                bs.push(est(&s0, &s1));
            }
            rows.push(json!({"from": b0, "to": b1, "local": pt, "ci": ci(bs.clone())}));
            pts_for_trend.push(((b0 as f64 + b1 as f64) / 2.0, bs, pt));
        }
        // trend: OLS slope of local exponent on mid-bits, over intervals with mid >= clo, bootstrap
        let mut trend = Value::Null;
        if pts_for_trend.len() >= 3 {
            let xs: Vec<f64> = pts_for_trend.iter().map(|t| t.0).collect();
            let ys: Vec<f64> = pts_for_trend.iter().map(|t| t.2).collect();
            let pt = ols(&xs, &ys);
            let mut tb = vec![];
            for i in 0..NB {
                let ys: Vec<f64> = pts_for_trend.iter().map(|t| t.1[i]).collect();
                tb.push(ols(&xs, &ys));
            }
            trend = json!({"slope_per_bit": pt, "ci": ci(tb)});
        }
        local.insert(a.clone(), json!({"intervals": rows, "trend": trend}));
    }
    summary.insert("local_exponents".into(), Value::Object(local));

    // 4. paired ratios per size vs basic, and aa vs ic2lp, ic2lp vs ic1lp
    let pairs = [
        ("ic2lp", "basic"),
        ("ic1lp", "basic"),
        ("ic2lpt", "basic"),
        ("aa", "ic2lp"),
        ("ic2lp", "ic1lp"),
        ("ic2lpt", "ic2lp"),
    ];
    let mut ratios = serde_json::Map::new();
    let mut rng = Rng(11);
    for (num, den) in pairs {
        if !d.contains_key(num) || !d.contains_key(den) {
            continue;
        }
        let mut rows = vec![];
        let allb: BTreeSet<u64> = d[num].keys().chain(d[den].keys()).cloned().collect();
        for b in allb {
            let cs = curves_common(&d, &[num, den], b);
            if cs.is_empty() {
                continue;
            }
            let lr: Vec<f64> = cs
                .iter()
                .map(|c| d[num][&b][c].1 - d[den][&b][c].1)
                .collect();
            let pt = (lr.iter().sum::<f64>() / lr.len() as f64).exp();
            let mut bs = vec![];
            for _ in 0..NB {
                let s: f64 =
                    (0..lr.len()).map(|_| lr[rng.below(lr.len())]).sum::<f64>() / lr.len() as f64;
                bs.push(s.exp());
            }
            rows.push(json!({"bits": b, "n": cs.len(), "ratio": pt, "ci": ci(bs), "per_curve": lr.iter().map(|x| x.exp()).collect::<Vec<_>>()}));
        }
        ratios.insert(format!("{}/{}", num, den), Value::Array(rows));
    }
    summary.insert("paired_ratios".into(), Value::Object(ratios));

    // 5. per-size table
    let mut table = vec![];
    let allb: BTreeSet<u64> = d
        .values()
        .flat_map(|m| m.keys().cloned())
        .chain(rho.keys().cloned())
        .collect();
    for b in &allb {
        let mut row = serde_json::Map::new();
        row.insert("bits".into(), json!(b));
        for a in &arms_all {
            if let Some(m) = d[a].get(b) {
                let cs: Vec<u64> = m.keys().cloned().collect();
                row.insert(format!("{}_gmean", a), json!(gmean(&d, a, *b, &cs)));
                row.insert(format!("{}_n", a), json!(cs.len()));
            }
            if let Some(e) = extra.get(&(a.clone(), *b)) {
                row.insert(format!("{}_rows", a), json!(e));
            }
        }
        if let Some(v) = rho.get(b) {
            let n = v.len() as f64;
            row.insert(
                "rho_cost_mean".into(),
                json!(v.iter().map(|t| t.0).sum::<f64>() / n),
            );
            row.insert(
                "rho_model_neg_mean".into(),
                json!(v.iter().map(|t| t.1).sum::<f64>() / n),
            );
            row.insert(
                "rho_steps_over_sqrt_ell_mean".into(),
                json!(v.iter().map(|t| t.2).sum::<f64>() / n),
            );
            row.insert("rho_n".into(), json!(v.len()));
            // ic2lp / rho model, per curve
            if let Some(m) = d.get("ic2lp").and_then(|m| m.get(b)) {
                let mut lr = vec![];
                for t in v {
                    if let Some(&(_, lc)) = m.get(&t.3) {
                        lr.push(lc - t.1.ln());
                    }
                }
                if !lr.is_empty() {
                    row.insert(
                        "ic2lp_over_rho_model".into(),
                        json!((lr.iter().sum::<f64>() / lr.len() as f64).exp()),
                    );
                }
            }
        }
        let mut corr = serde_json::Map::new();
        for ((a, bb), (ok, tot)) in &correct {
            if bb == b {
                corr.insert(a.clone(), json!(format!("{}/{}", ok, tot)));
            }
        }
        row.insert("correct".into(), Value::Object(corr));
        table.push(Value::Object(row));
    }
    summary.insert("per_size".into(), Value::Array(table));
    summary.insert(
        "instances".into(),
        json!(inst
            .iter()
            .map(|((b, c), v)| json!({"bits": b, "curve": c, "inst": v}))
            .collect::<Vec<_>>()),
    );
    let s = Value::Object(summary);
    fs::write(outp, serde_json::to_string_pretty(&s).unwrap()).unwrap();
    println!(
        "{}",
        serde_json::to_string_pretty(&s["common_range"]).unwrap()
    );
    println!(
        "{}",
        serde_json::to_string_pretty(&s["upper_range"]).unwrap()
    );
}
