//! Lead D4 / A2 falsifier: exact relation yield of the cyclotomic γ-weight
//! base `C_w` against a null population of random Frobenius-stable bases
//! with the same number of projected signed Frobenius columns, the
//! compact-orbit-scan base, the matched ONB Hamming-weight base, and a
//! planted positive control (the λ-progression base).
//!
//! Every base is moved into discrete-log space: each projected point `P`
//! of the order-`r` subgroup is replaced by `log_G P ∈ Z/r` (full table of
//! `[k]G`, exact).  Then, over **all** `r` targets:
//!
//! * `cov2` = |B + B| / r   (targets with at least one ordered 2-sum),
//! * `cov3` = |B + B + B| / r,
//! * `energy` = #{(a,b,c,d) ∈ B⁴ : a + b = c + d}  (additive energy),
//!
//! with `B` the signed, Frobenius-closed point set.  The mean number of
//! ordered `m`-decompositions per target is `|B|^m / r` for *every* base
//! of `|B|` points (an identity), so yield differences at matched size can
//! only come from how sums collide; that is what these three numbers see.
//!
//! The cyclotomic and compact-orbit-scan bases are built by the
//! `ic_framework` plug-ins themselves (`ecbench fb`'s
//! `methods::dump_factor_base`), so they are the bases `ic.pipeline` runs.
//! Protocol: `research/cyclotomic_falsifier_20261015/PROTOCOL.md`.
//!
//! ```bash
//! cargo run --release --example cyclotomic_falsifier_yield -- \
//!     --a 1 --n 23 --nulls 1999 --seed 20261015 --out OUT.json
//! ```
//!
//! Public synthetic Koblitz curves only; no target logarithm is sought.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::time::Instant;

use crypto_lib::cryptanalysis::cyclotomic_factor_base::CyclotomicField;
use crypto_lib::cryptanalysis::ecbench::methods::dump_factor_base;
use crypto_lib::cryptanalysis::ecbench::workload::{CurveSpec, Instance};
use crypto_lib::cryptanalysis::ic_boundary::BinaryInstance;
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};

struct Args {
    a: u8,
    n: u32,
    nulls: usize,
    seed: u64,
    out: String,
    /// Extra column counts to evaluate besides the pre-registered ones.
    extra_k: Vec<usize>,
    alpha: f64,
    family: usize,
}

fn parse_args() -> Args {
    let mut args = Args {
        a: 1,
        n: 23,
        nulls: 1999,
        seed: 20261015,
        out: "yield.json".into(),
        extra_k: Vec::new(),
        alpha: 0.05,
        family: 15,
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut value = || it.next().expect("flag value");
        match flag.as_str() {
            "--a" => args.a = value().parse().expect("a"),
            "--n" => args.n = value().parse().expect("n"),
            "--nulls" => args.nulls = value().parse().expect("nulls"),
            "--seed" => args.seed = value().parse().expect("seed"),
            "--out" => args.out = value(),
            "--extra-k" => {
                args.extra_k = value()
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(|s| s.parse().expect("k"))
                    .collect()
            }
            "--alpha" => args.alpha = value().parse().expect("alpha"),
            "--family" => args.family = value().parse().expect("family size"),
            other => panic!("unknown flag {other}"),
        }
    }
    args
}

/// splitmix64: the seeded stream for the null draws.
struct SplitMix(u64);
impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// Exact discrete logs in the order-`r` subgroup from a full table.
struct DlTable {
    r: u64,
    map: HashMap<u64, (u32, u64)>,
}

impl DlTable {
    fn new(inst: &BinaryInstance) -> Self {
        let r = inst.r;
        let half = (r - 1) / 2;
        let mut map = HashMap::with_capacity(half as usize + 1);
        let g = inst.generator;
        let mut p = g;
        for k in 1..=half {
            assert!(!p.infinity, "[k]G = O before r");
            let prev = map.insert(p.x, (k as u32, p.y));
            assert!(prev.is_none(), "x-collision in the subgroup table");
            p = inst.fast.add(p, g);
        }
        Self { r, map }
    }

    fn log(&self, p: &FastPoint) -> u64 {
        if p.infinity {
            return 0;
        }
        let &(k, y) = self.map.get(&p.x).expect("point outside the subgroup");
        if p.y == y {
            k as u64
        } else {
            self.r - k as u64
        }
    }
}

/// Logs of one signed Frobenius orbit (`±λ^j l`).
fn signed_orbit(l: u64, lambda: u64, n: u32, r: u64) -> Vec<u64> {
    let mut out = Vec::with_capacity(2 * n as usize);
    let mut v = l;
    for _ in 0..n {
        out.push(v);
        out.push((r - v) % r);
        v = mulmod(v, lambda, r);
    }
    assert_eq!(v, l, "orbit not closed");
    out
}

fn orbit_key(l: u64, lambda: u64, n: u32, r: u64) -> u64 {
    *signed_orbit(l, lambda, n, r).iter().min().unwrap()
}

struct Ctx<'a> {
    inst: &'a BinaryInstance,
    dl: DlTable,
    lambda: u64,
    n: u32,
    r: u64,
    h: u64,
}

impl Ctx<'_> {
    /// Projected column representatives (logs) from abscissae, in the
    /// order given, deduplicated by signed orbit — the plug-ins' rule.
    fn columns_from_abscissae(
        &self,
        xs: impl Iterator<Item = u64>,
        k: usize,
        seen: &mut HashSet<u64>,
        reps: &mut Vec<u64>,
    ) {
        for x in xs {
            if reps.len() >= k {
                return;
            }
            for p in self.inst.points_with_x(x) {
                let q = self.inst.fast.mul_u64(p, self.h);
                if q.infinity {
                    continue;
                }
                let l = self.dl.log(&q);
                if seen.insert(orbit_key(l, self.lambda, self.n, self.r)) {
                    reps.push(l);
                    if reps.len() >= k {
                        return;
                    }
                }
            }
        }
    }

    fn points_of(&self, reps: &[u64]) -> Vec<u64> {
        let mut pts: Vec<u64> = reps
            .iter()
            .flat_map(|&l| signed_orbit(l, self.lambda, self.n, self.r))
            .collect();
        pts.sort_unstable();
        pts.dedup();
        pts
    }
}

/// Scratch space for the exact sumset metrics.
struct Metrics {
    r: usize,
    c2: Vec<u32>,
    touched: Vec<u32>,
    s2: Vec<u64>,
    d: Vec<u64>,
    s3: Vec<u64>,
}

#[derive(Clone, Copy, Debug)]
struct Yield {
    points: usize,
    cov2: u64,
    cov3: u64,
    energy: u128,
}

impl Metrics {
    fn new(r: u64) -> Self {
        let r = r as usize;
        let words = r.div_ceil(64);
        Self {
            r,
            c2: vec![0; r],
            touched: Vec::new(),
            s2: vec![0; words],
            d: vec![0; 2 * words + 2],
            s3: vec![0; words],
        }
    }

    fn get_bit(v: &[u64], i: usize) -> u64 {
        (v[i >> 6] >> (i & 63)) & 1
    }

    fn measure(&mut self, pts: &[u64]) -> Yield {
        self.measure_opt(pts, true)
    }

    /// `want3 = false` skips the 3-sum pass (`cov3` reported as 0) on
    /// column counts where no pre-registered test reads it.
    fn measure_opt(&mut self, pts: &[u64], want3: bool) -> Yield {
        let r = self.r;
        // Ordered 2-sums.
        self.touched.clear();
        for &a in pts {
            for &b in pts {
                let mut t = a as usize + b as usize;
                if t >= r {
                    t -= r;
                }
                if self.c2[t] == 0 {
                    self.touched.push(t as u32);
                }
                self.c2[t] += 1;
            }
        }
        let cov2 = self.touched.len() as u64;
        let mut energy: u128 = 0;
        for w in self.s2.iter_mut() {
            *w = 0;
        }
        for &t in &self.touched {
            let c = self.c2[t as usize] as u128;
            energy += c * c;
            self.c2[t as usize] = 0;
            self.s2[t as usize >> 6] |= 1u64 << (t & 63);
        }
        if !want3 {
            return Yield {
                points: pts.len(),
                cov2,
                cov3: 0,
                energy,
            };
        }
        // 3-sums: S3 = ∪_b (S2 + b).
        for w in self.s3.iter_mut() {
            *w = 0;
        }
        let words = self.s2.len();
        if (cov2 as usize) * 64 < r {
            for &t in &self.touched {
                for &b in pts {
                    let mut u = t as usize + b as usize;
                    if u >= r {
                        u -= r;
                    }
                    self.s3[u >> 6] |= 1u64 << (u & 63);
                }
            }
        } else {
            // D = S2 ‖ S2 as a 2r-bit string; (S2 + b)[t] = D[t + r − b].
            for w in self.d.iter_mut() {
                *w = 0;
            }
            for i in 0..r {
                if Self::get_bit(&self.s2, i) == 1 {
                    self.d[i >> 6] |= 1 << (i & 63);
                    let j = i + r;
                    self.d[j >> 6] |= 1 << (j & 63);
                }
            }
            for &b in pts {
                let o = r - b as usize; // b ∈ [0, r)
                let o = if o == r { 0 } else { o };
                let (w0, s) = (o >> 6, o & 63);
                for i in 0..words {
                    let lo = self.d[w0 + i] >> s;
                    let hi = if s == 0 { 0 } else { self.d[w0 + i + 1] << (64 - s) };
                    self.s3[i] |= lo | hi;
                }
            }
            let tail = r & 63;
            if tail != 0 {
                self.s3[words - 1] &= (1u64 << tail) - 1;
            }
        }
        let cov3 = self.s3.iter().map(|w| w.count_ones() as u64).sum();
        Yield {
            points: pts.len(),
            cov2,
            cov3,
            energy,
        }
    }
}

/// Transition column counts: the K at which a random base of `2nK`
/// points would cover half the targets with `m`-sums,
/// `1 − exp(−(2nK)^m / (m!·r)) = 1/2`.
fn transition_k(r: u64, n: u32, m: u32) -> usize {
    let fact = (1..=m).product::<u32>() as f64;
    let pts = (fact * r as f64 * std::f64::consts::LN_2).powf(1.0 / m as f64);
    ((pts / (2.0 * n as f64)).round() as usize).max(1)
}

fn rank_stats(value: f64, null: &[f64]) -> Value {
    let n = null.len();
    let le = null.iter().filter(|&&v| v <= value).count();
    let ge = null.iter().filter(|&&v| v >= value).count();
    let p_lo = (le + 1) as f64 / (n + 1) as f64;
    let p_hi = (ge + 1) as f64 / (n + 1) as f64;
    let p_two = (2.0 * p_lo.min(p_hi)).min(1.0);
    let mut sorted = null.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |f: f64| sorted[((f * (n - 1) as f64).round() as usize).min(n - 1)];
    let (lo, hi) = (sorted[0], sorted[n - 1]);
    let mean = null.iter().sum::<f64>() / n as f64;
    let (med, q025, q975) = (q(0.5), q(0.025), q(0.975));
    json!({
        "value": value,
        "null_median": med, "null_q025": q025, "null_q975": q975,
        "null_min": lo, "null_max": hi,
        "null_mean": mean,
        "rank_le": le, "rank_ge": ge, "n_null": n,
        "p_lower": p_lo, "p_upper": p_hi, "p_two_sided": p_two,
        "all_null_equal": lo == hi,
    })
}

fn main() {
    let args = parse_args();
    let started = Instant::now();
    let spec = CurveSpec::Koblitz {
        a: args.a,
        n: args.n,
    };
    let spec_json = serde_json::to_string(&spec).unwrap();
    let inst = match spec.build().expect("curve") {
        Instance::Binary(i) => i,
        _ => panic!("single-word Koblitz curve expected"),
    };
    let kc = inst.koblitz.as_ref().expect("Koblitz");
    let lambda = kc.lambda.to_u64().expect("λ fits");
    let (n, r, h) = (inst.n, inst.r, inst.cofactor);
    eprintln!(
        "K_{}/F_2^{n}: #E = {}, r = {r}, h = {h}, λ = {lambda}",
        args.a, inst.group_order
    );
    let t = Instant::now();
    let dl = DlTable::new(&inst);
    eprintln!("dl table: {} entries [{:.1}s]", dl.map.len(), t.elapsed().as_secs_f64());
    let ctx = Ctx {
        inst: &inst,
        dl,
        lambda,
        n,
        r,
        h,
    };
    // Self-check: Frobenius acts as λ on logs.
    {
        let p = inst.fast.mul_u64(inst.generator, 12345 % r);
        let fp = inst.fast.frobenius_k(p, 1);
        assert_eq!(ctx.dl.log(&fp), mulmod(ctx.dl.log(&p), lambda, r), "Frobenius ≠ λ");
    }

    // ── Plug-in bases (what ic.pipeline runs) ──
    let mut checks = Vec::new();
    let plugin_reps = |fb_spec: &str, checks: &mut Vec<Value>| -> Option<(Vec<u64>, Value)> {
        let dump = match dump_factor_base(&spec, fb_spec) {
            Ok(d) => d,
            Err(e) => {
                eprintln!("  {fb_spec}: {e}");
                checks.push(json!({"spec": fb_spec, "error": e}));
                return None;
            }
        };
        let mut rep_of_col: HashMap<u64, u64> = HashMap::new();
        let mut logs = Vec::new();
        for p in &dump.points {
            let x = u64::from_str_radix(p.x.trim_start_matches("0x"), 16).unwrap();
            let y = u64::from_str_radix(p.y.trim_start_matches("0x"), 16).unwrap();
            let fp = FastPoint::affine(x, y);
            let l = ctx.dl.log(&fp);
            logs.push((p.col, p.coef, l));
            if p.coef == 1 {
                rep_of_col.insert(p.col, l);
            }
        }
        // Self-check: log(point) = coef · log(column representative).
        let consistent = logs
            .iter()
            .all(|&(c, coef, l)| rep_of_col.get(&c).map(|&rl| mulmod(rl, coef, r)) == Some(l));
        assert!(consistent, "{fb_spec}: dump coefficients disagree with logs");
        let mut cols: Vec<u64> = rep_of_col.keys().copied().collect();
        cols.sort_unstable();
        let reps: Vec<u64> = cols.iter().map(|c| rep_of_col[c]).collect();
        let facts = json!({
            "spec": fb_spec, "fb_id": dump.factor_base.fb_id,
            "columns": dump.factor_base.columns, "signed_points": dump.factor_base.signed_points,
            "coef_log_consistent": consistent,
        });
        checks.push(facts.clone());
        Some((reps, facts))
    };

    // Full cyclotomic bases (census) for w = 1, 2, 3.
    let mut cyc_full: Vec<(usize, Vec<u64>)> = Vec::new();
    let mut census = Vec::new();
    for w in 1..=3usize {
        if let Some((reps, facts)) = plugin_reps(&format!("cyclotomic:w={w}"), &mut checks) {
            census.push(json!({"w": w, "full_columns": reps.len(), "facts": facts}));
            cyc_full.push((w, reps));
        } else {
            census.push(json!({"w": w, "full_columns": 0}));
        }
    }
    let k2 = transition_k(r, n, 2);
    let k3 = transition_k(r, n, 3);
    eprintln!("transition K: m=2 → {k2}, m=3 → {k3}");

    // Pre-registered cells: (w, K, metrics tested at that K).
    // K ∈ {min(K2, Kfull) for cov2/energy, min(K3, Kfull) for cov3}.
    let mut cells: Vec<(usize, usize, Vec<&str>)> = Vec::new();
    for (w, reps) in &cyc_full {
        let full = reps.len();
        let ka = k2.min(full);
        let kb = k3.min(full);
        if ka == kb {
            cells.push((*w, ka, vec!["cov2", "energy", "cov3"]));
        } else {
            cells.push((*w, ka, vec!["cov2", "energy"]));
            cells.push((*w, kb, vec!["cov3"]));
        }
    }
    let mut ks: BTreeSet<usize> = cells.iter().map(|c| c.1).collect();
    ks.extend(args.extra_k.iter().copied());
    let kmax = *ks.iter().max().unwrap_or(&1);
    let ks3: BTreeSet<usize> = cells
        .iter()
        .filter(|c| c.2.contains(&"cov3"))
        .map(|c| c.1)
        .chain(args.extra_k.iter().copied())
        .collect();

    let mut m = Metrics::new(r);

    // ── Comparator bases at every K ──
    let onb_field = CyclotomicField::new(&inst.irreducible, n).expect("cyclotomic field");
    let mut comparators: Vec<Value> = Vec::new();
    let mut comparator_values: HashMap<(String, usize), Yield> = HashMap::new();
    for &k in &ks {
        let mut push = |kind: String, reps: &[u64], facts: Value, m: &mut Metrics| {
            let pts = ctx.points_of(reps);
            let y = m.measure(&pts);
            comparators.push(json!({
                "kind": kind, "k_requested": k, "columns": reps.len(), "points": y.points,
                "cov2": y.cov2, "cov3": y.cov3, "energy": y.energy.to_string(),
                "facts": facts,
            }));
            comparator_values.insert((kind, k), y);
        };
        for (w, full) in &cyc_full {
            if k <= full.len() {
                let spec_k = format!("cyclotomic:w={w},columns={k}");
                if let Some((reps, facts)) = plugin_reps(&spec_k, &mut checks) {
                    push(format!("cyclotomic_w{w}"), &reps, facts, &mut m);
                }
            }
        }
        let cos = format!("compact-orbit-scan:columns={k},raw_x_cap={}", (1u64 << n) - 1);
        if let Some((reps, facts)) = plugin_reps(&cos, &mut checks) {
            push("compact_orbit_scan".into(), &reps, facts, &mut m);
        }
        // Matched ONB Hamming-weight base: least weight v reaching K
        // projected columns, first K columns in the set's order.
        for v in 1..=(n as usize / 2) {
            let set = onb_field.onb_weight_set(v).expect("onb set");
            let mut seen = HashSet::new();
            let mut reps = Vec::new();
            ctx.columns_from_abscissae(
                set.iter().map(|e| e.raw_bits().first().copied().unwrap_or(0)),
                k,
                &mut seen,
                &mut reps,
            );
            if reps.len() == k {
                push("onb_weight".into(), &reps, json!({"onb_weight": v}), &mut m);
                break;
            }
        }
        // Positive control: the λ-progression base, column reps 1, 2, 3, …
        {
            let mut seen = HashSet::new();
            let mut reps = Vec::new();
            let mut l = 1u64;
            while reps.len() < k {
                if seen.insert(orbit_key(l, lambda, n, r)) {
                    reps.push(l);
                }
                l += 1;
            }
            push(
                "positive_lambda_progression".into(),
                &reps,
                json!({"reps": "logs 1,2,3,… (deduplicated by signed orbit)"}),
                &mut m,
            );
        }
    }
    eprintln!("comparators done [{:.1}s]", started.elapsed().as_secs_f64());

    // ── Null population: random abscissae, same projection and fold ──
    let t = Instant::now();
    let mut null: HashMap<usize, Vec<Yield>> = HashMap::new();
    let mask = (1u64 << n) - 1;
    for draw in 0..args.nulls {
        let mut rng = SplitMix(args.seed ^ (draw as u64).wrapping_mul(0x2545_f491_4f6c_dd1d));
        let mut seen = HashSet::new();
        let mut reps = Vec::new();
        while reps.len() < kmax {
            let x = loop {
                let x = rng.next() & mask;
                if x != 0 {
                    break x;
                }
            };
            ctx.columns_from_abscissae(std::iter::once(x), kmax, &mut seen, &mut reps);
        }
        for &k in &ks {
            let pts = ctx.points_of(&reps[..k]);
            null.entry(k).or_default().push(m.measure_opt(&pts, ks3.contains(&k)));
        }
        if draw % 200 == 0 {
            eprintln!("  null draw {draw} [{:.1}s]", t.elapsed().as_secs_f64());
        }
    }
    let null_secs = t.elapsed().as_secs_f64();

    // ── Tests ──
    let crit = args.alpha / args.family as f64;
    let mut tests = Vec::new();
    for (w, k, metrics) in &cells {
        let nulls = &null[k];
        for kind in [
            format!("cyclotomic_w{w}"),
            "compact_orbit_scan".into(),
            "onb_weight".into(),
            "positive_lambda_progression".into(),
        ] {
            let Some(y) = comparator_values.get(&(kind.clone(), *k)) else {
                continue;
            };
            for &metric in metrics {
                let pick = |y: &Yield| match metric {
                    "cov2" => y.cov2 as f64,
                    "cov3" => y.cov3 as f64,
                    _ => y.energy as f64,
                };
                let nv: Vec<f64> = nulls.iter().map(pick).collect();
                let mut st = rank_stats(pick(y), &nv);
                let p = st["p_two_sided"].as_f64().unwrap();
                let med = st["null_median"].as_f64().unwrap();
                let v = pick(y);
                let direction = if v < med {
                    "below"
                } else if v > med {
                    "above"
                } else {
                    "at"
                };
                st["metric"] = json!(metric);
                st["base"] = json!(kind);
                st["w_cell"] = json!(w);
                st["k"] = json!(k);
                st["direction"] = json!(direction);
                st["significant"] = json!(p < crit);
                st["effect_vs_median"] = json!(if med != 0.0 { v / med - 1.0 } else { 0.0 });
                tests.push(st);
            }
        }
    }

    let null_dump: Value = ks
        .iter()
        .map(|k| {
            let v = &null[k];
            (
                k.to_string(),
                json!({
                    "cov2": v.iter().map(|y| y.cov2).collect::<Vec<_>>(),
                    "cov3": v.iter().map(|y| y.cov3).collect::<Vec<_>>(),
                    "energy": v.iter().map(|y| y.energy.to_string()).collect::<Vec<_>>(),
                    "points": v.iter().map(|y| y.points).collect::<Vec<_>>(),
                }),
            )
        })
        .collect::<serde_json::Map<String, Value>>()
        .into();

    let report = json!({
        "schema": "cyclotomic_falsifier_yield/v1",
        "curve": spec_json, "a": args.a, "n": n, "r": r, "cofactor": h,
        "group_order": inst.group_order, "lambda": lambda,
        "seed": args.seed, "nulls": args.nulls,
        "alpha": args.alpha, "family": args.family, "critical_p": crit,
        "transition_k": {"m2": k2, "m3": k3},
        "census": census, "cells": cells.iter().map(|(w,k,ms)| json!({"w": w, "k": k, "metrics": ms})).collect::<Vec<_>>(),
        "plugin_checks": checks,
        "comparators": comparators,
        "tests": tests,
        "null": null_dump,
        "null_seconds": null_secs,
        "total_seconds": started.elapsed().as_secs_f64(),
    });
    std::fs::write(&args.out, serde_json::to_string(&report).unwrap() + "\n").expect("write");
    eprintln!("wrote {} [{:.1}s]", args.out, started.elapsed().as_secs_f64());
}
