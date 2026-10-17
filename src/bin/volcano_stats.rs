//! volcano-stats: the statistics of an isogeny-volcano cost experiment
//! (`research/isogeny-volcano-cost_20261017/PLAN.md`).
//!
//! Reads ecbench `records.jsonl` files and the `curves.jsonl` sidecars that
//! `volcano-rho --ecbench` writes, and decomposes each arm's cost into:
//!
//! - **vertical**: the floor/crater ratio of means, pooled over classes
//!   after normalising each class by its own mean, with a two-stage
//!   bootstrap interval (curves within level, runs within curve) and a
//!   permutation p (level labels shuffled among curves within each class);
//! - **horizontal**: heterogeneity between curves of the same level (one
//!   endomorphism ring, so horizontally isogenous), as a Cochran Q with a
//!   permutation null (runs reassigned among the level's curves) and the
//!   DerSimonian–Laird spread τ beside the null's 95th percentile;
//! - **model**: the same test between isomorphic models of one curve, the
//!   coordinate effect an isogeny effect has to beat;
//! - **exactness**: for BSGS, the number of distinct costs per paired target
//!   across every curve and model (1 when transport and pairing hold);
//! - **pairs**: the yield or cost ratio of two arms (the planted control).
//!
//! A run's value is a ratio `num/den`: S is `(s, 1)`, relation yield is
//! `(relations, trials)`, so means of means are never taken. All numbers
//! are computed here; Python only draws `analysis.json`.
//!
//!   volcano-stats --records R1.jsonl[,R2…] --curves C1.jsonl[,C2…] [--margin 0.05]
//!                 [--boot 2000] [--perm 10000] [--pair ic-mitm-plus1:ic-mitm] [--out analysis.json]
//!                 [--plant-level L:FACTOR] [--plant-curve V:FACTOR]   (positive controls of the analysis)

use std::collections::{BTreeMap, HashMap};
use std::fs;

use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

#[derive(Clone)]
struct Run {
    class: u64,
    vertex: u64,
    model: u64,
    level: u64,
    /// Paired target index: the same k on every curve of the class.
    target: u64,
    num: f64,
    den: f64,
}

/// A records or sidecar file, gzipped (`.gz`, as committed) or plain.
fn read_text(path: &str) -> String {
    use std::io::Read;
    let bytes = fs::read(path).unwrap_or_else(|e| die(&format!("{path}: {e}")));
    if !path.ends_with(".gz") {
        return String::from_utf8(bytes).unwrap_or_else(|e| die(&format!("{path}: {e}")));
    }
    let mut out = String::new();
    flate2::read::GzDecoder::new(&bytes[..])
        .read_to_string(&mut out)
        .unwrap_or_else(|e| die(&format!("{path}: {e}")));
    out
}

fn die(msg: &str) -> ! {
    eprintln!("volcano-stats: {msg}");
    std::process::exit(2);
}

fn ratio(rs: &[&Run]) -> f64 {
    let (n, d) = rs
        .iter()
        .fold((0.0, 0.0), |(n, d), r| (n + r.num, d + r.den));
    n / d
}

/// Ratio estimator and its variance (runs are the sampling units).
fn ratio_se2(rs: &[&Run]) -> (f64, f64) {
    let m = ratio(rs);
    let d: f64 = rs.iter().map(|r| r.den).sum();
    let k = rs.len() as f64;
    let ss: f64 = rs.iter().map(|r| (r.num - m * r.den).powi(2)).sum();
    (m, ss / (k - 1.0).max(1.0) * k / (d * d))
}

fn quantile(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((q * (v.len() - 1) as f64).round() as usize).min(v.len() - 1)]
}

/// Runs normalised by their class mean, so classes of different size pool.
fn normalise(runs: &[Run]) -> Vec<Run> {
    let mut by: HashMap<u64, Vec<&Run>> = HashMap::new();
    for r in runs {
        by.entry(r.class).or_default().push(r);
    }
    let mean: HashMap<u64, f64> = by.iter().map(|(c, rs)| (*c, ratio(rs))).collect();
    runs.iter()
        .map(|r| Run {
            num: r.num / mean[&r.class],
            ..r.clone()
        })
        .collect()
}

fn group<K: Ord + Clone>(runs: &[Run], key: impl Fn(&Run) -> K) -> BTreeMap<K, Vec<&Run>> {
    let mut m: BTreeMap<K, Vec<&Run>> = BTreeMap::new();
    for r in runs {
        m.entry(key(r)).or_default().push(r);
    }
    m
}

/// Floor/crater ratio of the pooled, class-normalised means.
/// `--vertical stratified` (PROTOCOL v1.1): the floor/crater ratio within
/// each class, combined as the equal-weight mean of log ratios. The v1
/// pooled ratio mixes classes in different floor and crater proportions, so
/// a shift in one class's mean (a resampled target) moves it even when no
/// curve differs from another: BSGS, identical on every curve, gets a
/// spurious interval.
static STRATIFIED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// `--alternative`: 0 two-sided (|ln V|), 1 less (V < 1), 2 greater (V > 1)
/// for the vertical permutation p (PROTOCOL v2's one-sided hypothesis).
static ALTERNATIVE: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(0);

/// Floor/crater ratio: pooled over classes (v1) or stratified by class (v1.1).
fn vertical(runs: &[Run], floor: &HashMap<u64, u64>) -> Option<f64> {
    if STRATIFIED.load(std::sync::atomic::Ordering::Relaxed) {
        let (mut acc, mut k) = (0.0, 0usize);
        for (c, rs) in group(runs, |r| r.class) {
            let crater: Vec<&Run> = rs.iter().copied().filter(|r| r.level == 0).collect();
            let bottom: Vec<&Run> = rs
                .iter()
                .copied()
                .filter(|r| r.level == floor[&c] && r.level > 0)
                .collect();
            if !crater.is_empty() && !bottom.is_empty() {
                acc += (ratio(&bottom) / ratio(&crater)).ln();
                k += 1;
            }
        }
        return (k > 0).then(|| (acc / k as f64).exp());
    }
    let crater: Vec<&Run> = runs.iter().filter(|r| r.level == 0).collect();
    let bottom: Vec<&Run> = runs
        .iter()
        .filter(|r| r.level == floor[&r.class] && r.level > 0)
        .collect();
    (!crater.is_empty() && !bottom.is_empty()).then(|| ratio(&bottom) / ratio(&crater))
}

/// Cochran Q of cell means around their stratum means, the DL spread τ, and
/// the largest standardised deviation of one cell (the test for a single
/// deviant curve, which the omnibus Q dilutes over its degrees of freedom).
fn cochran(
    runs: &[Run],
    stratum: impl Fn(&Run) -> (u64, u64, u64),
    cell: impl Fn(&Run) -> (u64, u64, u64, u64),
) -> (f64, f64, usize, f64) {
    let mut q = 0.0;
    let mut maxz = 0f64;
    let (mut sw, mut sw2, mut df) = (0.0, 0.0, 0usize);
    for (_, rs) in group(runs, &stratum) {
        let owned: Vec<Run> = rs.iter().map(|r| (*r).clone()).collect();
        let cells = group(&owned, &cell);
        if cells.len() < 2 {
            continue;
        }
        let stats: Vec<(f64, f64)> = cells
            .values()
            .map(|c| ratio_se2(c))
            .filter(|s| s.1 > 0.0)
            .collect();
        let w: Vec<f64> = stats.iter().map(|s| 1.0 / s.1).collect();
        let tw: f64 = w.iter().sum();
        let mbar = stats.iter().zip(&w).map(|(s, w)| s.0 * w).sum::<f64>() / tw;
        q += stats
            .iter()
            .zip(&w)
            .map(|(s, w)| w * (s.0 - mbar).powi(2))
            .sum::<f64>();
        maxz = stats
            .iter()
            .zip(&w)
            .map(|(s, w)| (s.0 - mbar).abs() * w.sqrt())
            .fold(maxz, f64::max);
        sw += tw;
        sw2 += w.iter().map(|w| w * w).sum::<f64>();
        df += stats.len() - 1;
    }
    let tau2 = ((q - df as f64) / (sw - sw2 / sw.max(1e-300))).max(0.0);
    (q, tau2.sqrt(), df, maxz)
}

/// Q's permutation null: runs reassigned among the cells of each stratum.
fn cochran_null(
    runs: &[Run],
    stratum: impl Fn(&Run) -> (u64, u64, u64) + Copy,
    cell: impl Fn(&Run) -> (u64, u64, u64, u64) + Copy,
    perms: usize,
    rng: &mut SmallRng,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let strata = group(runs, stratum);
    let mut qs = Vec::with_capacity(perms);
    let mut taus = Vec::with_capacity(perms);
    let mut zs = Vec::with_capacity(perms);
    for _ in 0..perms {
        let mut shuffled = Vec::with_capacity(runs.len());
        for rs in strata.values() {
            let mut vals: Vec<(f64, f64)> = rs.iter().map(|r| (r.num, r.den)).collect();
            vals.shuffle(rng);
            for (r, (n, d)) in rs.iter().zip(vals) {
                shuffled.push(Run {
                    num: n,
                    den: d,
                    ..(*r).clone()
                });
            }
        }
        let (q, t, _, z) = cochran(&shuffled, stratum, cell);
        qs.push(q);
        taus.push(t);
        zs.push(z);
    }
    (qs, taus, zs)
}

fn analyse(
    arm: &str,
    raw: &[Run],
    floor: &HashMap<u64, u64>,
    margin: f64,
    boot: usize,
    perms: usize,
    rng: &mut SmallRng,
) -> Value {
    let runs = normalise(raw);
    // Level profile.
    let profile: Vec<Value> = group(&runs, |r| r.level)
        .iter()
        .map(|(l, rs)| {
            let (m, se2) = ratio_se2(rs);
            let curves = group(&rs.iter().map(|r| (*r).clone()).collect::<Vec<_>>(), |r| (r.class, r.vertex)).len();
            json!({"level": l, "runs": rs.len(), "curves": curves, "rel_mean": m, "ci95_half": 1.96 * se2.sqrt()})
        })
        .collect();

    // Vertical: bootstrap (curves within class×level, then runs within curve) and permutation.
    let v_obs = vertical(&runs, floor);
    let mut vert = json!(null);
    if let Some(v) = v_obs {
        let cells = group(&runs, |r| (r.class, r.level, r.vertex));
        let mut by_stratum: BTreeMap<(u64, u64), Vec<&Vec<&Run>>> = BTreeMap::new();
        for ((c, l, _), rs) in &cells {
            by_stratum.entry((*c, *l)).or_default().push(rs);
        }
        // Paired targets: one draw of target indices per class, used by every
        // curve, so a cost fixed by its target (BSGS) has no spurious spread.
        let targets: BTreeMap<u64, Vec<u64>> = group(&runs, |r| r.class)
            .iter()
            .map(|(c, rs)| {
                let mut t: Vec<u64> = rs.iter().map(|r| r.target).collect();
                t.sort();
                t.dedup();
                (*c, t)
            })
            .collect();
        let mut bs = Vec::with_capacity(boot);
        for _ in 0..boot {
            let draw: HashMap<u64, Vec<u64>> = targets
                .iter()
                .map(|(c, t)| {
                    (
                        *c,
                        (0..t.len()).map(|_| t[rng.gen_range(0..t.len())]).collect(),
                    )
                })
                .collect();
            let mut sample = Vec::new();
            for ((class, _), curves) in &by_stratum {
                for _ in 0..curves.len() {
                    let c = curves[rng.gen_range(0..curves.len())];
                    let owned: Vec<Run> = c.iter().map(|r| (*r).clone()).collect();
                    let by_t = group(&owned, |r| r.target);
                    for t in &draw[class] {
                        if let Some(rs) = by_t.get(t) {
                            sample.extend(rs.iter().map(|r| (*r).clone()));
                        }
                    }
                }
            }
            if let Some(b) = vertical(&sample, floor) {
                bs.push(b);
            }
        }
        // Permutation: each class's level labels shuffled among its curves.
        let by_class = group(&runs, |r| r.class);
        let mut ge = 0usize;
        for _ in 0..perms {
            let mut shuffled = Vec::with_capacity(runs.len());
            for rs in by_class.values() {
                let owned: Vec<Run> = rs.iter().map(|r| (*r).clone()).collect();
                let curves: Vec<(u64, u64)> = group(&owned, |r| (r.vertex, r.level))
                    .keys()
                    .copied()
                    .collect();
                let mut labels: Vec<u64> = curves.iter().map(|c| c.1).collect();
                labels.shuffle(rng);
                let relabel: HashMap<u64, u64> =
                    curves.iter().zip(labels).map(|(c, l)| (c.0, l)).collect();
                for r in rs {
                    shuffled.push(Run {
                        level: relabel[&r.vertex],
                        ..(*r).clone()
                    });
                }
            }
            if let Some(p) = vertical(&shuffled, floor) {
                let extreme = match ALTERNATIVE.load(std::sync::atomic::Ordering::Relaxed) {
                    1 => p.ln() <= v.ln() + 1e-12,
                    2 => p.ln() >= v.ln() - 1e-12,
                    _ => p.ln().abs() >= v.ln().abs() - 1e-12,
                };
                if extreme {
                    ge += 1;
                }
            }
        }
        let (lo90, hi90) = (
            quantile(&mut bs.clone(), 0.05),
            quantile(&mut bs.clone(), 0.95),
        );
        let (lo95, hi95) = (
            quantile(&mut bs.clone(), 0.025),
            quantile(&mut bs.clone(), 0.975),
        );
        vert = json!({
            "floor_over_crater": v, "ci90": [lo90, hi90], "ci95": [lo95, hi95],
            "perm_p": (ge + 1) as f64 / (perms + 1) as f64, "perms": perms, "boot": bs.len(),
            "equivalent_within_margin": lo90 >= 1.0 - margin && hi90 <= 1.0 + margin, "margin": margin,
        });
    }

    // Horizontal: curves within (class, level); model: models within curve.
    let hs = |r: &Run| (r.class, r.level, 0u64);
    let hc = |r: &Run| (r.class, r.level, r.vertex, 0u64);
    let (qh, tauh, dfh, zh) = cochran(&runs, hs, hc);
    let (mut qhn, mut tauhn, zhn) = cochran_null(&runs, hs, hc, perms.min(2000), rng);
    let ms = |r: &Run| (r.class, r.vertex, 0u64);
    let mc = |r: &Run| (r.class, r.vertex, r.model, 0u64);
    let (qm, taum, dfm, zm) = cochran(&runs, ms, mc);
    let (mut qmn, mut taumn, zmn) = cochran_null(&runs, ms, mc, perms.min(2000), rng);
    let pval = |obs: f64, null: &[f64]| {
        (null.iter().filter(|&&q| q >= obs - 1e-12).count() + 1) as f64 / (null.len() + 1) as f64
    };
    let horizontal = json!({
        "q": qh, "df": dfh, "perm_p": pval(qh, &qhn), "tau_rel": tauh,
        "max_abs_z": zh, "max_z_perm_p": pval(zh, &zhn),
        "tau_null_p95": if tauhn.is_empty() { 0.0 } else { quantile(&mut tauhn, 0.95) },
        "q_null_p95": if qhn.is_empty() { 0.0 } else { quantile(&mut qhn, 0.95) },
    });
    let model = json!({
        "q": qm, "df": dfm, "perm_p": pval(qm, &qmn), "tau_rel": taum,
        "max_abs_z": zm, "max_z_perm_p": pval(zm, &zmn),
        "tau_null_p95": if taumn.is_empty() { 0.0 } else { quantile(&mut taumn, 0.95) },
        "q_null_p95": if qmn.is_empty() { 0.0 } else { quantile(&mut qmn, 0.95) },
    });

    // Adjacent same-level pairs (direct horizontal ℓ-edges), descriptive.
    let curve_cells = group(&runs, |r| (r.class, r.vertex));
    let curve_means: HashMap<(u64, u64), (f64, f64)> = curve_cells
        .iter()
        .map(|(k, rs)| (*k, ratio_se2(rs)))
        .collect();
    let curves: Vec<Value> = curve_cells
        .iter()
        .map(|((c, v), rs)| {
            let (m, se2) = curve_means[&(*c, *v)];
            json!({"class": c, "vertex": v, "level": rs[0].level, "runs": rs.len(), "rel_mean": m, "ci95_half": 1.96 * se2.sqrt()})
        })
        .collect();
    json!({
        "arm": arm, "runs": raw.len(), "levels": profile, "vertical": vert,
        "horizontal": horizontal, "model": model, "curves": curves,
        "curve_means": curve_means.len(),
    })
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| {
        argv.iter()
            .position(|a| a == k)
            .and_then(|i| argv.get(i + 1))
            .cloned()
    };
    if argv.iter().any(|a| a == "-h" || a == "--help") || get("--records").is_none() {
        println!("usage: volcano-stats --records R1.jsonl[,R2...] --curves C1.jsonl[,C2...] [--margin 0.05] [--boot 2000] [--perm 10000] [--pair A:B]... [--out analysis.json]");
        return;
    }
    let margin: f64 = get("--margin").map_or(0.05, |v| v.parse().unwrap());
    let boot: usize = get("--boot").map_or(2000, |v| v.parse().unwrap());
    let perms: usize = get("--perm").map_or(10000, |v| v.parse().unwrap());
    let seed: u64 = get("--seed").map_or(20261017, |v| v.parse().unwrap());
    let mut rng = SmallRng::seed_from_u64(seed);
    let alternative = get("--alternative").unwrap_or_else(|| "two-sided".into());
    ALTERNATIVE.store(
        match alternative.as_str() {
            "two-sided" => 0,
            "less" => 1,
            "greater" => 2,
            _ => die("--alternative takes two-sided, less or greater"),
        },
        std::sync::atomic::Ordering::Relaxed,
    );
    let vmode = get("--vertical").unwrap_or_else(|| "pooled".into());
    match vmode.as_str() {
        "pooled" => {}
        "stratified" => STRATIFIED.store(true, std::sync::atomic::Ordering::Relaxed),
        _ => die("--vertical takes pooled (PROTOCOL v1) or stratified (v1.1)"),
    }
    let plant = |k: &str| {
        get(k).map(|v| {
            let (a, b) = v
                .split_once(':')
                .unwrap_or_else(|| die(&format!("{k} takes INDEX:FACTOR")));
            (a.parse::<u64>().unwrap(), b.parse::<f64>().unwrap())
        })
    };
    let (plant_level, plant_curve) = (plant("--plant-level"), plant("--plant-curve"));

    // Sidecars: (p, a, b, gx, gy) → (class, vertex, model, level). Each
    // sidecar file is one class (two classes may share p).
    let mut side: HashMap<(u64, u64, u64, u64, u64), (u64, u64, u64, u64)> = HashMap::new();
    let mut floor: HashMap<u64, u64> = HashMap::new();
    let mut class_names = Vec::new();
    for (ci, path) in get("--curves")
        .unwrap_or_default()
        .split(',')
        .filter(|s| !s.is_empty())
        .enumerate()
    {
        let ci = ci as u64;
        class_names.push(json!({"class": ci, "curves": path}));
        for line in read_text(path).lines() {
            let c: Value = serde_json::from_str(line).unwrap();
            let u = |k: &str| c[k].as_u64().unwrap();
            side.insert(
                (u("p"), u("a"), u("b"), u("gx"), u("gy")),
                (ci, u("vertex"), u("model"), u("level")),
            );
            let f = floor.entry(ci).or_insert(0);
            *f = (*f).max(u("level"));
        }
    }

    let mut by_arm: BTreeMap<String, Vec<Run>> = BTreeMap::new();
    let mut ic_yield: BTreeMap<String, Vec<Run>> = BTreeMap::new();
    let mut bsgs: BTreeMap<String, BTreeMap<(u64, u64), Vec<f64>>> = BTreeMap::new();
    let (mut skipped, mut unverified) = (0usize, 0usize);
    for path in get("--records").unwrap().split(',') {
        for line in read_text(path).lines() {
            let r: Value = serde_json::from_str(line).unwrap();
            if r["warmup"].as_bool() == Some(true) {
                continue;
            }
            if r["outcome"]["status"] != "verified" {
                unverified += 1;
                continue;
            }
            let cs = &r["workload"]["curve_spec"];
            let key = |k: &str| cs[k].as_u64().unwrap_or(u64::MAX);
            let Some(&(class, vertex, model, level)) =
                side.get(&(key("p"), key("a"), key("b"), key("gx"), key("gy")))
            else {
                skipped += 1;
                continue;
            };
            let arm = r["arm"].as_str().unwrap().to_string();
            let target = r["workload"]["target_index"].as_u64().unwrap();
            let s = r["cost"]["s"].as_f64().unwrap();
            let mut base = Run {
                class,
                vertex,
                model,
                level,
                target,
                num: s,
                den: 1.0,
            };
            // Planted effects (--plant-level L:x, --plant-curve V:x): the
            // analysis's own positive control, applied to S and to yield.
            if plant_level.is_some_and(|(l, _)| l == level) {
                base.num *= plant_level.unwrap().1;
            }
            if plant_curve.is_some_and(|(v, _)| v == vertex) {
                base.num *= plant_curve.unwrap().1;
            }
            if r["method"]["family"] == "bsgs" {
                bsgs.entry(arm.clone())
                    .or_default()
                    .entry((base.class, target))
                    .or_default()
                    .push(s);
            }
            if let (Some(rel), Some(tr)) = (
                r["counters"]["relations_found"].as_f64(),
                r["counters"]["targets_tried"].as_f64(),
            ) {
                // A planted yield effect scales relations per trial.
                let k = base.num / s;
                ic_yield.entry(arm.clone()).or_default().push(Run {
                    num: rel * k,
                    den: tr,
                    ..base.clone()
                });
            }
            by_arm.entry(arm).or_default().push(base);
        }
    }
    if unverified > 0 {
        eprintln!("volcano-stats: {unverified} unverified runs: every comparison involving them is incomplete");
    }

    let mut metrics = Vec::new();
    for (arm, runs) in &by_arm {
        eprintln!("{arm}: S over {} runs", runs.len());
        let mut v = analyse(arm, runs, &floor, margin, boot, perms, &mut rng);
        v["metric"] = json!("S");
        metrics.push(v);
    }
    for (arm, runs) in &ic_yield {
        eprintln!("{arm}: relation yield over {} runs", runs.len());
        let mut v = analyse(arm, runs, &floor, margin, boot, perms, &mut rng);
        v["metric"] = json!("relation_yield");
        metrics.push(v);
    }
    let exact: Vec<Value> = bsgs
        .iter()
        .map(|(arm, t)| {
            let worst = t.values().map(|v| {
                let mut d: Vec<u64> = v.iter().map(|x| x.to_bits()).collect();
                d.sort();
                d.dedup();
                d.len()
            });
            let worst = worst.max().unwrap_or(0);
            json!({"arm": arm, "paired_targets": t.len(), "max_distinct_cost_per_target": worst, "exact": worst == 1})
        })
        .collect();

    // Planted / named arm pairs: ratio of pooled class-normalised yield (or S), bootstrap over runs.
    let mut pairs = Vec::new();
    for spec in argv
        .windows(2)
        .filter(|w| w[0] == "--pair")
        .map(|w| w[1].clone())
    {
        // Matched by workload (class, curve, model, target), so classes with
        // different group orders never mix in the ratio.
        let (a, b) = spec
            .split_once(':')
            .unwrap_or_else(|| die("--pair takes ARM_A:ARM_B"));
        let src = if ic_yield.contains_key(a) {
            &ic_yield
        } else {
            &by_arm
        };
        let (Some(ra), Some(rb)) = (src.get(a), src.get(b)) else {
            die(&format!("--pair {spec}: unknown arm"))
        };
        let key = |r: &Run| (r.class, r.vertex, r.model, r.target);
        let bmap: HashMap<_, &Run> = rb.iter().map(|r| (key(r), r)).collect();
        let matched: Vec<(&Run, &Run)> = ra
            .iter()
            .filter_map(|r| bmap.get(&key(r)).map(|o| (r, *o)))
            .collect();
        if matched.is_empty() {
            die(&format!("--pair {spec}: no matched workloads"));
        }
        // Per class, then the class-weighted mean of log ratios.
        let classes: Vec<u64> = {
            let mut c: Vec<u64> = matched.iter().map(|m| m.0.class).collect();
            c.sort();
            c.dedup();
            c
        };
        let stat = |ms: &[(&Run, &Run)]| -> f64 {
            let mut acc = 0.0;
            for c in &classes {
                let xs: Vec<&(&Run, &Run)> = ms.iter().filter(|m| m.0.class == *c).collect();
                let (na, da, nb, db) = xs.iter().fold((0.0, 0.0, 0.0, 0.0), |s, m| {
                    (s.0 + m.0.num, s.1 + m.0.den, s.2 + m.1.num, s.3 + m.1.den)
                });
                acc += ((na / da) / (nb / db)).ln();
            }
            (acc / classes.len() as f64).exp()
        };
        let obs = stat(&matched);
        let mut bs: Vec<f64> = (0..boot)
            .map(|_| {
                let s: Vec<(&Run, &Run)> = (0..matched.len())
                    .map(|_| matched[rng.gen_range(0..matched.len())])
                    .collect();
                stat(&s)
            })
            .collect();
        let ci = [quantile(&mut bs, 0.025), quantile(&mut bs, 0.975)];
        pairs.push(json!({"a": a, "b": b, "metric": if std::ptr::eq(src, &ic_yield) { "relation_yield" } else { "S" },
                          "matched_workloads": matched.len(), "classes": classes.len(),
                          "ratio": obs, "ci95": ci, "excludes_1": ci[0] > 1.0 || ci[1] < 1.0}));
    }

    let out = json!({
        "schema": "volcano-stats/v1", "vertical_statistic": vmode, "alternative": alternative, "margin": margin, "boot": boot, "perm": perms, "seed": seed,
        "classes": class_names, "skipped_unmatched_records": skipped, "unverified": unverified,
        "metrics": metrics, "bsgs_exactness": exact, "pairs": pairs,
        "argv": std::env::args().collect::<Vec<_>>(),
    });
    let text = serde_json::to_string_pretty(&out).unwrap();
    match get("--out") {
        Some(p) => fs::write(&p, text).unwrap(),
        None => println!("{text}"),
    }
    // A one-screen summary.
    for m in out["metrics"].as_array().unwrap() {
        let v = &m["vertical"];
        eprintln!(
            "{:>16} {:>15}: floor/crater {:.4} 90% [{:.4}, {:.4}] p={:.4} equiv={} | horizontal tau {:.4} (null p95 {:.4}) p={:.4} maxz p={:.4} | model tau {:.4} (null p95 {:.4}) p={:.4} maxz p={:.4}",
            m["arm"].as_str().unwrap(), m["metric"].as_str().unwrap(),
            v["floor_over_crater"].as_f64().unwrap_or(f64::NAN), v["ci90"][0].as_f64().unwrap_or(f64::NAN), v["ci90"][1].as_f64().unwrap_or(f64::NAN),
            v["perm_p"].as_f64().unwrap_or(f64::NAN), v["equivalent_within_margin"],
            m["horizontal"]["tau_rel"].as_f64().unwrap_or(f64::NAN), m["horizontal"]["tau_null_p95"].as_f64().unwrap_or(f64::NAN), m["horizontal"]["perm_p"].as_f64().unwrap_or(f64::NAN),
            m["horizontal"]["max_z_perm_p"].as_f64().unwrap_or(f64::NAN),
            m["model"]["tau_rel"].as_f64().unwrap_or(f64::NAN), m["model"]["tau_null_p95"].as_f64().unwrap_or(f64::NAN), m["model"]["perm_p"].as_f64().unwrap_or(f64::NAN),
            m["model"]["max_z_perm_p"].as_f64().unwrap_or(f64::NAN),
        );
    }
    for e in &exact {
        eprintln!("bsgs exactness: {e}");
    }
    for p in &pairs {
        eprintln!("pair: {p}");
    }
}
