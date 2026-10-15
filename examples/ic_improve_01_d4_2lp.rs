//! IC improvement round, item 01: the `D₄` Y-line arm stacked with double
//! large primes (2LP) for `k = 4` index calculus on `E(F_{p⁴})`, against the
//! `D₄` incumbent, on one shared trial stream per curve.  Protocol:
//! `research/ic_improve_20261014/01-d4-2lp/PROTOCOL.md`.
//!
//! ```text
//! ic_improve_01_d4_2lp train   --p 521 --seed 4 --json raw/train.jsonl
//! ic_improve_01_d4_2lp tune    --in raw/train.jsonl --out raw/frozen.json
//! ic_improve_01_d4_2lp test    --p 521 --seed 1 --frozen raw/frozen.json --json raw/test.jsonl
//! ic_improve_01_d4_2lp analyze --in raw/test.jsonl --prior <item-1 ladder.jsonl> --out raw/summary.json
//! ```
//!
//! The oracle, the Y-line base, the trial stream and the charge are item 1's
//! (`examples/ic_next10_01_fp4_ladder.rs`), copied unchanged so the
//! incumbent reproduces item 1's trial counts.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;

use crypto_lib::cryptanalysis::gaudry_cubic::{
    square_core, wiedemann_u64, LargePrimeEliminator, SparseRel,
};
use crypto_lib::cryptanalysis::gaudry_quartic::{Curve4, Pt4, E4};
use crypto_lib::cryptanalysis::gaudry_quartic_d4::{y_base, D4Instance};
use crypto_lib::cryptanalysis::quartic_folds::{
    generate_translation_instance, rho4_folded, FoldInstance4, RhoFold,
};
use crypto_lib::cryptanalysis::residual_walk::mix64;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

/// `F_p` multiplications per `E(F_{p⁴})` addition (§11.16, measured §11.17).
const C_ADD: f64 = 97.0;
/// `F_p` multiplications charged per multiplication mod `r` (§11.16).
const C_MODN: f64 = 16.0;
/// `gaudry_quartic_d4_20261013`'s measured `C₄` for the `D₄` solve.
const C4_D4: f64 = 2.54e10;

fn fp(x: &E4) -> u64 {
    mix64(x.0[0] ^ mix64(x.0[1] ^ mix64(x.0[2] ^ mix64(x.0[3] ^ 0x5eed))))
}

fn signed(curve: &Curve4, p: &Pt4, s: i64) -> Pt4 {
    if s == 1 {
        *p
    } else {
        curve.neg(p)
    }
}

/// The exact 4-term decomposition oracle over one base.
struct Oracle {
    pts: Vec<Pt4>,
    /// `x(P_i + s P_j) ↦ (i, j, s, y)`.
    pairs: HashMap<E4, Vec<(u32, u32, i8, E4)>>,
    /// Sorted fingerprints of `x(±P_i ± P_j ± P_k)`, `i < j < k`.
    triples: Vec<u64>,
    setup_ops: u64,
}

impl Oracle {
    fn build(curve: &Curve4, pts: Vec<Pt4>) -> Oracle {
        let before = curve.ops();
        let m = pts.len();
        let mut pairs: HashMap<E4, Vec<(u32, u32, i8, E4)>> = HashMap::new();
        let cap = 4 * m * (m.saturating_sub(1)) * (m.saturating_sub(2)) / 6;
        let mut triples: Vec<u64> = Vec::with_capacity(cap);
        for i in 0..m {
            for j in (i + 1)..m {
                for s in [1i8, -1] {
                    let q = curve.add(&pts[i], &signed(curve, &pts[j], s as i64));
                    if q.inf {
                        continue;
                    }
                    pairs
                        .entry(q.x)
                        .or_default()
                        .push((i as u32, j as u32, s, q.y));
                    for pk in &pts[(j + 1)..] {
                        for t in [curve.add(&q, pk), curve.sub(&q, pk)] {
                            if !t.inf {
                                triples.push(fp(&t.x));
                            }
                        }
                    }
                }
            }
        }
        triples.sort_unstable();
        Oracle {
            pts,
            pairs,
            triples,
            setup_ops: curve.ops() - before,
        }
    }

    /// Every `R = Σ s_t P_{i_t}` over four distinct base points, signs
    /// included, sorted by index.
    fn decompose(&self, curve: &Curve4, r: &Pt4) -> Vec<Vec<(usize, i64)>> {
        let m = self.pts.len();
        let mut out: BTreeSet<Vec<(usize, i64)>> = BTreeSet::new();
        for l in 0..m {
            for sl in [1i64, -1] {
                let v = curve.sub(r, &signed(curve, &self.pts[l], sl));
                if v.inf || self.triples.binary_search(&fp(&v.x)).is_err() {
                    continue;
                }
                for k in 0..m {
                    if k == l {
                        continue;
                    }
                    for sk in [1i64, -1] {
                        let w = curve.sub(&v, &signed(curve, &self.pts[k], sk));
                        if w.inf {
                            continue;
                        }
                        let Some(c) = self.pairs.get(&w.x) else {
                            continue;
                        };
                        for &(i, j, s, y) in c {
                            let (i, j) = (i as usize, j as usize);
                            if i == k || i == l || j == k || j == l {
                                continue;
                            }
                            let t: i64 = if w.y == y { 1 } else { -1 };
                            let mut terms = vec![(l, sl), (k, sk), (i, t), (j, t * s as i64)];
                            terms.sort_unstable();
                            out.insert(terms);
                        }
                    }
                }
            }
        }
        out.into_iter().collect()
    }
}

/// One arm's configuration.
#[derive(Clone, Debug)]
struct Cfg {
    name: String,
    /// `None`: no large primes (the incumbent).  `Some(α)`: `s = ⌈ncols^α⌉`;
    /// `Some(f64::INFINITY)` is the null (`s = ncols`).
    alpha: Option<f64>,
    cap: u8,
    minlp: bool,
    omega: f64,
    /// Use only trials with `t % stride == 0` (positive control: 2).
    stride: u64,
}

impl Cfg {
    fn incumbent() -> Cfg {
        Cfg {
            name: "d4".into(),
            alpha: None,
            cap: 0,
            minlp: false,
            omega: 1.0,
            stride: 1,
        }
    }
    fn lp(name: String, alpha: f64, cap: u8, minlp: bool, omega: f64) -> Cfg {
        Cfg {
            name,
            alpha: Some(alpha),
            cap,
            minlp,
            omega,
            stride: 1,
        }
    }
    fn key(alpha: f64, cap: u8, minlp: bool, omega: f64) -> String {
        format!(
            "a{alpha:.2}-c{cap}-{}-w{omega:.1}",
            if minlp { "minlp" } else { "uniform" }
        )
    }
    fn json(&self) -> Value {
        json!({"alpha": self.alpha.map(|a| if a.is_finite() { json!(a) } else { json!("ncols") }),
               "merge_cap": self.cap, "pick": if self.minlp {"minlp"} else {"uniform"},
               "omega": self.omega, "stride": self.stride})
    }
}

/// One IC arm's state on the shared trial stream.
struct Arm {
    cfg: Cfg,
    col: Vec<(usize, i64)>,
    small: usize,
    max_lp: usize,
    elim: Option<LargePrimeEliminator>,
    full: Vec<SparseRel>,
    next_attempt: usize,
    la_ops: u64,
    merge_ops: u64,
    la_attempts: u64,
    targets_used: u64,
    lp_hist: [u64; 5],
    solved: bool,
    correct: bool,
    recovered: Option<u64>,
    unknowns: usize,
    la_rng: StdRng,
    snapshot: BTreeMap<&'static str, u64>,
}

impl Arm {
    fn new(cfg: Cfg, ycol: &[(usize, i64)], ncols: usize, r: u64, seed: u64) -> Arm {
        let (small, col, max_lp) = match cfg.alpha {
            None => (ncols, ycol.to_vec(), 0),
            Some(a) => {
                let s = if a.is_finite() {
                    ((ncols as f64).powf(a).ceil() as usize).min(ncols)
                } else {
                    ncols
                };
                let col = ycol
                    .iter()
                    .map(|&(c, k)| (if c < s { c } else { c + 1 }, k))
                    .collect();
                (s, col, 2)
            }
        };
        let elim = (max_lp > 0).then(|| LargePrimeEliminator::with_cap(small + 1, r, cfg.cap));
        let first = ((cfg.omega * (small + 1) as f64).ceil() as usize).max(small + 1);
        Arm {
            cfg,
            col,
            small,
            max_lp,
            elim,
            full: Vec::new(),
            next_attempt: first,
            la_ops: 0,
            merge_ops: 0,
            la_attempts: 0,
            targets_used: 0,
            lp_hist: [0; 5],
            solved: false,
            correct: false,
            recovered: None,
            unknowns: 0,
            la_rng: StdRng::seed_from_u64(seed ^ 0x1A1A),
            snapshot: BTreeMap::new(),
        }
    }

    fn row(&self, terms: &[(usize, i64)], a: u64, b: u64, r: u64) -> (SparseRel, usize) {
        let mut acc: BTreeMap<usize, u64> = BTreeMap::new();
        for &(i, s) in terms {
            let (c, k) = self.col[i];
            let v = if s * k == 1 { 1 } else { r - 1 };
            let e = acc.entry(c).or_insert(0);
            *e = (*e + v) % r;
        }
        let e = acc.entry(self.small).or_insert(0);
        *e = (*e + (r - b) % r) % r;
        let cols: Vec<(usize, u64)> = acc.into_iter().filter(|&(_, v)| v != 0).collect();
        let lp = cols.iter().filter(|&&(c, _)| c > self.small).count();
        (SparseRel { cols, rhs: a }, lp)
    }

    /// Feed trial `t`'s verified decompositions; true once solved.
    #[allow(clippy::too_many_arguments)]
    fn feed(
        &mut self,
        t: u64,
        decs: &[Vec<(usize, i64)>],
        a: u64,
        b: u64,
        inst: &FoldInstance4,
        curve: &Curve4,
        pick: u64,
    ) -> bool {
        let r = inst.n;
        if t % self.cfg.stride != 0 {
            return false;
        }
        let mut rows: Vec<(SparseRel, usize)> = decs
            .iter()
            .map(|d| self.row(d, a, b, r))
            .filter(|(_, lp)| *lp <= self.max_lp)
            .collect();
        if self.cfg.minlp {
            if let Some(m) = rows.iter().map(|x| x.1).min() {
                rows.retain(|x| x.1 == m);
            }
        }
        if rows.is_empty() {
            return false;
        }
        let (rel, lp) = rows[(pick % rows.len() as u64) as usize].clone();
        self.targets_used += 1;
        self.lp_hist[lp.min(4)] += 1;
        let full = if lp == 0 {
            Some(rel)
        } else {
            self.elim
                .as_mut()
                .expect("large primes need an eliminator")
                .feed(rel, &mut self.merge_ops)
        };
        let Some(full) = full else { return false };
        self.full.push(full);
        let unknowns = self.small + 1;
        if self.full.len() < self.next_attempt {
            return false;
        }
        let Some(core) = square_core(&self.full, self.small) else {
            self.next_attempt = self.full.len() + (unknowns / 20).max(1);
            return false;
        };
        self.la_attempts += 1;
        let mut map = vec![usize::MAX; unknowns];
        let mut k = 0;
        for (c, slot) in map.iter_mut().enumerate() {
            if core.columns.get(c).copied().unwrap_or(false) {
                *slot = k;
                k += 1;
            }
        }
        let sel: Vec<SparseRel> = core
            .rows
            .iter()
            .map(|&i| SparseRel {
                cols: self.full[i]
                    .cols
                    .iter()
                    .map(|&(c, v)| (map[c], v))
                    .collect(),
                rhs: self.full[i].rhs,
            })
            .collect();
        let x = wiedemann_u64(&sel, sel.len(), r, &mut self.la_rng, &mut self.la_ops);
        let dd = x.map(|x| x[map[self.small]]);
        if dd.is_some_and(|dd| curve.mul(&inst.g, dd) == inst.q) {
            self.solved = true;
            self.recovered = dd;
            self.correct = dd == Some(inst.d);
            self.unknowns = sel.len();
            return true;
        }
        self.next_attempt = self.full.len() + (sel.len() / 20).max(1);
        false
    }
}

fn trial_ab(seed: u64, t: u64, r: u64) -> (u64, u64, u64) {
    let mut rng = StdRng::seed_from_u64(mix64(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ t));
    (rng.gen_range(0..r), rng.gen_range(1..r), rng.gen())
}

/// Item 1's x-line lifts, run only to consume the RNG in the same order.
fn x_line_base(curve: &Curve4, rng: &mut StdRng) -> Vec<Pt4> {
    (0..curve.f.p)
        .filter_map(|x| curve.lift_x(&curve.f.from_fp(x), rng))
        .collect()
}

/// Y-line base and its fold: `P_{p−y} = ±(P_y + T)`.
fn y_line(
    inst: &FoldInstance4,
    d4: &D4Instance,
    rng: &mut StdRng,
) -> (Vec<(u64, Pt4)>, Vec<(usize, i64)>, usize) {
    let curve = &d4.curve;
    let p = inst.p;
    let base = y_base(d4, rng);
    let t = Pt4 {
        x: inst.x0,
        y: E4::ZERO,
        inf: false,
    };
    let idx: HashMap<u64, usize> = base.iter().enumerate().map(|(i, (y, _))| (*y, i)).collect();
    let mut col = vec![(usize::MAX, 1i64); base.len()];
    let mut ncols = 0;
    for (i, (y, pt)) in base.iter().enumerate() {
        if col[i].0 != usize::MAX {
            continue;
        }
        col[i] = (ncols, 1);
        if let Some(&j) = idx.get(&((p - y) % p)) {
            if j != i {
                let q = curve.add(pt, &t);
                let s = if q == base[j].1 {
                    1
                } else if q == curve.neg(&base[j].1) {
                    -1
                } else {
                    panic!("P_y + T is not ±P_(p−y)")
                };
                col[j] = (ncols, s);
            }
        }
        ncols += 1;
    }
    (base, col, ncols)
}

fn d4_of(inst: &FoldInstance4) -> D4Instance {
    D4Instance {
        curve: Curve4::with_params(inst.p, inst.a, inst.b),
        x0: inst.x0,
        c: inst.c,
    }
}

fn arg<T: std::str::FromStr>(args: &[String], k: &str, d: T) -> T {
    args.iter()
        .position(|a| a == k)
        .and_then(|i| args.get(i + 1))
        .map(|v| v.parse().ok().expect(k))
        .unwrap_or(d)
}

fn emit(path: &str, v: &Value) {
    let mut f = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .expect("json");
    writeln!(f, "{v}").unwrap();
    println!("{v}");
}

const ALPHAS: [f64; 5] = [0.75, 0.80, 0.85, 0.90, 0.95];
const CAPS: [u8; 3] = [0, 2, 4];
const OMEGAS: [f64; 2] = [1.0, 1.1];

fn grid() -> Vec<Cfg> {
    let mut v = Vec::new();
    for &a in &ALPHAS {
        for &c in &CAPS {
            for m in [false, true] {
                for &w in &OMEGAS {
                    v.push(Cfg::lp(Cfg::key(a, c, m, w), a, c, m, w));
                }
            }
        }
    }
    v
}

fn null_cfg() -> Cfg {
    Cfg::lp("d4-2lp-null".into(), f64::INFINITY, 0, false, 1.0)
}

/// One curve: build the Y-line oracle, run every arm on the shared stream.
fn run_curve(args: &[String], cfgs: Vec<Cfg>, stop_mult: Option<f64>, split: &str) {
    let p: u64 = arg(args, "--p", 269);
    let seed: u64 = arg(args, "--seed", 1);
    let rho_runs: u64 = arg(args, "--rho-runs", 0);
    let json_path: String = arg(args, "--json", "out.jsonl".to_string());
    let t0 = Instant::now();
    let inst = generate_translation_instance(p, seed, 16).expect("instance");
    let r = inst.n;
    let sqrt_r = (r as f64).sqrt();
    let xcurve = inst.curve();
    let mut rng = StdRng::seed_from_u64(seed ^ 0x01F4_1AD0);
    let _ = x_line_base(&xcurve, &mut rng);
    let d4 = d4_of(&inst);
    let curve = &d4.curve;
    curve.f.reset_muls();
    let (yb, ycol, ncols) = y_line(&inst, &d4, &mut rng);
    let fb_gae = curve.f.muls() as f64 / C_ADD;
    let pts: Vec<Pt4> = yb.iter().map(|(_, q)| *q).collect();
    let tb = Instant::now();
    let oracle = Oracle::build(curve, pts);
    let build_s = tb.elapsed().as_secs_f64();
    let mut arms: Vec<Arm> = cfgs
        .into_iter()
        .map(|c| Arm::new(c, &ycol, ncols, r, seed))
        .collect();
    let cap = 2000 * ncols as u64;
    let (mut trials, mut scalar_ops, mut verify_ops, mut decs_total) = (0u64, 0u64, 0u64, 0u64);
    while trials < cap && arms.iter().any(|a| !a.solved) {
        if let Some(m) = stop_mult {
            if arms[0].solved && trials as f64 >= m * arms[0].snapshot["trials"] as f64 {
                break;
            }
        }
        let t = trials;
        let (a, b, pick) = trial_ab(seed, t, r);
        trials += 1;
        let o0 = curve.ops();
        let rr = curve.add(&curve.mul(&inst.g, a), &curve.mul(&inst.q, b));
        let o1 = curve.ops();
        let decs = oracle.decompose(curve, &rr);
        let o2 = curve.ops();
        let mut good = Vec::new();
        for terms in decs {
            let mut acc = rr;
            for &(i, s) in &terms {
                acc = curve.sub(&acc, &signed(curve, &oracle.pts[i], s));
            }
            if acc.inf {
                good.push(terms);
            }
        }
        let o3 = curve.ops();
        scalar_ops += o1 - o0;
        verify_ops += o3 - o2;
        decs_total += good.len() as u64;
        for arm in arms.iter_mut().filter(|a| !a.solved) {
            if arm.feed(t, &good, a, b, &inst, curve, pick) {
                arm.snapshot.insert("trials", trials);
                arm.snapshot.insert("scalar_ops", scalar_ops);
                arm.snapshot.insert("verify_ops", verify_ops);
            }
        }
    }
    let check = 1.5 * (r as f64).log2();
    let out: Vec<Value> = arms
        .iter()
        .map(|arm| {
            let g = |k: &str, d: u64| *arm.snapshot.get(k).unwrap_or(&d) as f64;
            let tr = g("trials", trials);
            let scalar = g("scalar_ops", scalar_ops);
            let ver = g("verify_ops", verify_ops);
            let solver = tr * C4_D4 / C_ADD;
            let la = (arm.la_ops + arm.merge_ops) as f64 * C_MODN / C_ADD;
            let rel = scalar + solver + ver;
            let whole = fb_gae + rel + la + check;
            json!({
                "arm": arm.cfg.name, "cfg": arm.cfg.json(),
                "solved": arm.solved, "correct": arm.correct, "recovered": arm.recovered,
                "trials": tr, "small_columns": arm.small, "unknowns_core": arm.unknowns,
                "full_relations": arm.full.len(), "targets_used": arm.targets_used,
                "lp_hist": arm.lp_hist, "la_attempts": arm.la_attempts,
                "la_modr_mults": arm.la_ops, "merge_modr_mults": arm.merge_ops,
                "pivots": arm.elim.as_ref().map(|e| e.pivot_count()),
                "abandoned": arm.elim.as_ref().map(|e| e.abandoned),
                "mean_merge_depth": arm.elim.as_ref().map(|e| e.mean_merge_depth()),
                "gae": {"factor_base": fb_gae, "trial_scalar": scalar, "solver_charged": solver,
                        "verify": ver, "relation": rel, "la": la, "final_check": check, "whole": whole},
                "S_charged": whole / sqrt_r,
            })
        })
        .collect();
    let mut row = json!({
        "item": "ic_improve_01_d4_2lp", "split": split,
        "p": p, "seed": seed, "r": r, "log2_r": (r as f64).log2(),
        "group_order": inst.group_order, "log2_E": (inst.group_order as f64).log2(),
        "cofactor": inst.cofactor, "planted_d": inst.d, "c_add": C_ADD, "c_modn": C_MODN, "c4": C4_D4,
        "y_base": yb.len(), "y_columns": ncols, "oracle_setup_ops": oracle.setup_ops,
        "oracle_build_s": build_s, "stream_trials": trials, "stream_decompositions": decs_total,
        "stop_mult": stop_mult, "arms": out,
    });
    drop(oracle);
    if rho_runs > 0 {
        let runs: Vec<Value> = (0..rho_runs)
            .map(|k| {
                let w = rho4_folded(&inst, RhoFold::Negation, seed * 1000 + k);
                json!({"steps": w.steps, "group_ops": w.group_ops, "correct": w.correct,
                       "S": w.group_ops as f64 / sqrt_r})
            })
            .collect();
        let m = runs.iter().map(|v| v["S"].as_f64().unwrap()).sum::<f64>() / runs.len() as f64;
        row["rho_neg"] = json!({"mean_S": m, "runs": runs});
    }
    row["wall_s"] = json!(t0.elapsed().as_secs_f64());
    emit(&json_path, &row);
}

fn train(args: &[String]) {
    let mut cfgs = vec![Cfg::incumbent(), null_cfg()];
    cfgs.extend(grid());
    run_curve(args, cfgs, Some(8.0), "train");
}

fn frozen_cfg(path: &str) -> Cfg {
    let v: Value = serde_json::from_str(&std::fs::read_to_string(path).expect("frozen")).unwrap();
    let f = &v["frozen"];
    Cfg::lp(
        "d4-2lp".into(),
        f["alpha"].as_f64().unwrap(),
        f["merge_cap"].as_u64().unwrap() as u8,
        f["pick"] == "minlp",
        f["omega"].as_f64().unwrap(),
    )
}

fn test(args: &[String]) {
    let frozen: String = arg(args, "--frozen", "raw/frozen.json".to_string());
    let mut half = Cfg::incumbent();
    half.name = "d4-half".into();
    half.stride = 2;
    let cfgs = vec![
        Cfg::incumbent(),
        frozen_cfg(&frozen),
        Cfg::lp("d4-2lp-gaudry".into(), 0.75, 0, false, 1.0),
        null_cfg(),
        half,
    ];
    run_curve(args, cfgs, None, "test");
}

// ── Analysis (Rust; Python only plots).

fn load(path: &str) -> Vec<Value> {
    std::fs::read_to_string(path)
        .expect("jsonl")
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}

fn arm_of<'a>(row: &'a Value, name: &str) -> &'a Value {
    row["arms"]
        .as_array()
        .unwrap()
        .iter()
        .find(|a| a["arm"] == name)
        .unwrap_or_else(|| panic!("arm {name}"))
}

fn null_equal(row: &Value) -> bool {
    let (a, b) = (arm_of(row, "d4"), arm_of(row, "d4-2lp-null"));
    a["trials"] == b["trials"]
        && a["full_relations"] == b["full_relations"]
        && a["la_modr_mults"] == b["la_modr_mults"]
        && a["recovered"] == b["recovered"]
}

fn tune(args: &[String]) {
    let input: String = arg(args, "--in", "raw/train.jsonl".to_string());
    let out_path: String = arg(args, "--out", "raw/frozen.json".to_string());
    let min_p: u64 = arg(args, "--min-p", 521);
    let rows = load(&input);
    let mut nulls = true;
    let mut inc_ok = true;
    let mut table: Vec<Value> = Vec::new();
    let mut best: Option<(f64, Cfg)> = None;
    for c in grid() {
        let (mut lw, mut lw_all, mut lrel, mut lla) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut unsolved = 0;
        for row in &rows {
            let i = arm_of(row, "d4");
            let a = arm_of(row, &c.name);
            let rw =
                (a["gae"]["whole"].as_f64().unwrap() / i["gae"]["whole"].as_f64().unwrap()).log2();
            lw_all.push(rw);
            if row["p"].as_u64().unwrap() >= min_p {
                lw.push(rw);
                lrel.push(
                    (a["gae"]["relation"].as_f64().unwrap()
                        / i["gae"]["relation"].as_f64().unwrap())
                    .log2(),
                );
                lla.push(
                    (a["gae"]["la"].as_f64().unwrap() / i["gae"]["la"].as_f64().unwrap()).log2(),
                );
            }
            if a["solved"] != true {
                unsolved += 1;
            } else if a["correct"] != true {
                unsolved += 1000;
            }
        }
        let gm = |v: &[f64]| 2f64.powf(v.iter().sum::<f64>() / v.len() as f64);
        let score = gm(&lw);
        table.push(json!({"config": c.name, "cfg": c.json(), "score_gm_ratio_p_ge_min": score,
            "gm_ratio_all_train": gm(&lw_all), "gm_relation_ratio": gm(&lrel), "gm_la_ratio": gm(&lla),
            "unsolved_or_wrong": unsolved, "n": lw.len()}));
        if best.as_ref().is_none_or(|(s, _)| score < *s) {
            best = Some((score, c.clone()));
        }
    }
    for row in &rows {
        nulls &= null_equal(row);
        let i = arm_of(row, "d4");
        inc_ok &= i["solved"] == true && i["correct"] == true;
    }
    table.sort_by(|a, b| {
        a["score_gm_ratio_p_ge_min"]
            .as_f64()
            .partial_cmp(&b["score_gm_ratio_p_ge_min"].as_f64())
            .unwrap()
    });
    let (score, c) = best.unwrap();
    let out = json!({
        "rule": format!("argmin over the 60-config grid of the geometric-mean paired whole-cost ratio vs d4 on training curves with p >= {min_p}"),
        "training_curves": rows.len(), "null_identical_all_train": nulls, "incumbent_all_correct": inc_ok,
        "frozen": {"name": c.name, "alpha": c.alpha.unwrap(), "merge_cap": c.cap,
                   "pick": if c.minlp {"minlp"} else {"uniform"}, "omega": c.omega, "score": score},
        "ranking": table,
    });
    std::fs::write(&out_path, serde_json::to_string_pretty(&out).unwrap()).unwrap();
    println!("{}", serde_json::to_string_pretty(&json!({"frozen": out["frozen"], "null": nulls, "top5": &out["ranking"].as_array().unwrap()[..5], "gaudry": out["ranking"].as_array().unwrap().iter().find(|t| t["config"] == "a0.75-c0-uniform-w1.0")})).unwrap());
}

fn ols(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
    let b = sxy / sxx;
    (b, my - b * mx)
}

fn quant(mut v: Vec<f64>, f: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() as f64 - 1.0) * f).round() as usize]
}

/// OLS of `ly` (already log2) on `xs` with a 95 % bootstrap over rows.
fn fit_log(xs: &[f64], ly: &[f64], rng: &mut StdRng) -> Value {
    let (b, a) = ols(xs, ly);
    let n = xs.len();
    let mut bs = Vec::new();
    let mut is = Vec::new();
    for _ in 0..2000 {
        let idx: Vec<usize> = (0..n).map(|_| rng.gen_range(0..n)).collect();
        let x: Vec<f64> = idx.iter().map(|&i| xs[i]).collect();
        if x.iter().all(|v| (v - x[0]).abs() < 1e-9) {
            continue;
        }
        let y: Vec<f64> = idx.iter().map(|&i| ly[i]).collect();
        let (bb, aa) = ols(&x, &y);
        bs.push(bb);
        is.push(aa);
    }
    json!({"slope": b, "intercept_log2": a, "ci95": [quant(bs.clone(), 0.025), quant(bs, 0.975)], "n": n})
}

/// Geometric mean of `2^ly` with a 95 % bootstrap over entries.
fn gm_ci(ly: &[f64], rng: &mut StdRng) -> Value {
    let n = ly.len();
    let m = ly.iter().sum::<f64>() / n as f64;
    let bs: Vec<f64> = (0..2000)
        .map(|_| (0..n).map(|_| ly[rng.gen_range(0..n)]).sum::<f64>() / n as f64)
        .collect();
    json!({"ratio": 2f64.powf(m), "ci95": [2f64.powf(quant(bs.clone(), 0.025)), 2f64.powf(quant(bs, 0.975))], "n": n})
}

fn verdict_ratio(v: &Value) -> &'static str {
    let lo = v["ci95"][0].as_f64().unwrap();
    let hi = v["ci95"][1].as_f64().unwrap();
    if hi < 1.0 {
        "ADVANCES"
    } else if lo > 1.0 {
        "REGRESSES"
    } else {
        "MATCHES"
    }
}

fn verdict_slope(v: &Value) -> &'static str {
    let lo = v["ci95"][0].as_f64().unwrap();
    let hi = v["ci95"][1].as_f64().unwrap();
    if hi < 0.0 {
        "ADVANCES"
    } else if lo > 0.0 {
        "REGRESSES"
    } else {
        "MATCHES"
    }
}

fn analyze(args: &[String]) {
    let input: String = arg(args, "--in", "raw/test.jsonl".to_string());
    let prior: String = arg(args, "--prior", String::new());
    let out_path: String = arg(args, "--out", "raw/summary.json".to_string());
    let rows = load(&input);
    let mut rng = StdRng::seed_from_u64(20261014);
    // Prior item-1 rows: d4 trials and rho-negation S, keyed by (p, seed).
    let mut prior_map: HashMap<(u64, u64), (f64, f64)> = HashMap::new();
    if !prior.is_empty() {
        for r in load(&prior) {
            let d4 = r["arms"]
                .as_array()
                .unwrap()
                .iter()
                .find(|a| a["arm"] == "d4")
                .unwrap();
            prior_map.insert(
                (r["p"].as_u64().unwrap(), r["seed"].as_u64().unwrap()),
                (
                    d4["trials"].as_f64().unwrap(),
                    r["rho"]["rho-neg"]["mean_S"].as_f64().unwrap(),
                ),
            );
        }
    }
    let arms = ["d4", "d4-2lp", "d4-2lp-gaudry", "d4-2lp-null", "d4-half"];
    let phases = ["whole", "relation", "la"];
    let mut verified = true;
    let mut nulls = true;
    let mut repro = (0, 0);
    let mut xs = Vec::new();
    let mut costs: BTreeMap<(&str, &str), Vec<f64>> = BTreeMap::new();
    let mut table = Vec::new();
    for row in &rows {
        let l2e = row["log2_E"].as_f64().unwrap();
        xs.push(l2e);
        nulls &= null_equal(row);
        let key = (row["p"].as_u64().unwrap(), row["seed"].as_u64().unwrap());
        let inc = arm_of(row, "d4");
        let pr = prior_map.get(&key);
        if let Some(&(t, _)) = pr {
            repro.1 += 1;
            if (t - inc["trials"].as_f64().unwrap()).abs() < 0.5 {
                repro.0 += 1;
            }
        }
        let rho_neg = row["rho_neg"]["mean_S"].as_f64().or(pr.map(|x| x.1));
        let mut cells = serde_json::Map::new();
        for a in arms {
            let v = arm_of(row, a);
            verified &= v["solved"] == true && v["correct"] == true;
            for ph in phases {
                costs
                    .entry((a, ph))
                    .or_default()
                    .push(v["gae"][ph].as_f64().unwrap());
            }
            cells.insert(a.into(), json!({"S_charged": v["S_charged"], "trials": v["trials"],
                "small_columns": v["small_columns"], "unknowns_core": v["unknowns_core"],
                "lp_hist": v["lp_hist"], "relation_gae": v["gae"]["relation"], "la_gae": v["gae"]["la"],
                "whole_gae": v["gae"]["whole"], "solved": v["solved"], "correct": v["correct"],
                "ratio_vs_d4": v["gae"]["whole"].as_f64().unwrap() / inc["gae"]["whole"].as_f64().unwrap(),
                "S_over_rho_neg_reference": rho_neg.map(|x| v["S_charged"].as_f64().unwrap() / x)}));
        }
        table.push(json!({"p": key.0, "seed": key.1, "log2_E": l2e, "log2_r": row["log2_r"],
            "cofactor": row["cofactor"], "y_columns": row["y_columns"], "rho_neg_S_reference": rho_neg,
            "prior_d4_trials": pr.map(|x| x.0), "arms": cells}));
    }
    // Per-arm, per-phase exponents.
    let mut fits = serde_json::Map::new();
    for a in arms {
        let mut m = serde_json::Map::new();
        for ph in phases {
            let ly: Vec<f64> = costs[&(a, ph)].iter().map(|c| c.log2()).collect();
            m.insert(ph.into(), fit_log(&xs, &ly, &mut rng));
        }
        fits.insert(a.into(), Value::Object(m));
    }
    // Paired ratios against the incumbent.
    let ps: Vec<u64> = {
        let mut v: Vec<u64> = rows.iter().map(|r| r["p"].as_u64().unwrap()).collect();
        v.dedup();
        v
    };
    let mut paired = serde_json::Map::new();
    for a in ["d4-2lp", "d4-2lp-gaudry", "d4-half", "d4-2lp-null"] {
        let mut m = serde_json::Map::new();
        for ph in phases {
            let lr: Vec<f64> = costs[&(a, ph)]
                .iter()
                .zip(&costs[&("d4", ph)])
                .map(|(c, i)| (c / i).log2())
                .collect();
            let pooled = gm_ci(&lr, &mut rng);
            let trend = fit_log(&xs, &lr, &mut rng);
            let per_size: Vec<Value> = ps
                .iter()
                .map(|&p| {
                    let v: Vec<f64> = rows
                        .iter()
                        .zip(&lr)
                        .filter(|(r, _)| r["p"].as_u64() == Some(p))
                        .map(|(_, x)| *x)
                        .collect();
                    let l2e = rows
                        .iter()
                        .filter(|r| r["p"].as_u64() == Some(p))
                        .map(|r| r["log2_E"].as_f64().unwrap())
                        .sum::<f64>()
                        / v.len() as f64;
                    let mut c = gm_ci(&v, &mut rng);
                    c["p"] = json!(p);
                    c["log2_E"] = json!(l2e);
                    c["per_curve"] = json!(v.iter().map(|x| 2f64.powf(*x)).collect::<Vec<_>>());
                    c
                })
                .collect();
            // Fitted handover: where the trend line of log2 ratio crosses 0.
            let b = trend["slope"].as_f64().unwrap();
            let i0 = trend["intercept_log2"].as_f64().unwrap();
            let cross = if b.abs() > 1e-12 { Some(-i0 / b) } else { None };
            m.insert(ph.into(), json!({
                "pooled": pooled, "constant_verdict": verdict_ratio(&pooled),
                "trend_slope_log2_ratio_per_log2_E": trend, "exponent_verdict": verdict_slope(&trend),
                "fitted_ratio_one_at_log2_E": cross, "per_size": per_size}));
        }
        paired.insert(a.into(), Value::Object(m));
    }
    // Phase-model handover: each arm's relation and LA phases on their own
    // fitted power laws; first log2 #E in [20, 600] where the candidate's
    // modelled whole cost is at or below the incumbent's.
    let model = |a: &str, e: f64| -> f64 {
        ["relation", "la"]
            .iter()
            .map(|ph| {
                let f = &fits[a][*ph];
                2f64.powf(f["intercept_log2"].as_f64().unwrap() + f["slope"].as_f64().unwrap() * e)
            })
            .sum()
    };
    let mut handover = serde_json::Map::new();
    for a in ["d4-2lp", "d4-2lp-gaudry"] {
        let mut flips: Vec<f64> = Vec::new();
        let mut e = 20.0;
        let below_now = model(a, 20.0) <= model("d4", 20.0);
        let mut state = below_now;
        while e <= 600.0 {
            if (model(a, e) <= model("d4", e)) != state {
                state = !state;
                flips.push((e * 10.0).round() / 10.0);
            }
            e += 0.1;
        }
        let first = json!(flips);
        let la_share = {
            let f = &fits["d4"]["la"];
            2f64.powf(f["intercept_log2"].as_f64().unwrap() + f["slope"].as_f64().unwrap() * 128.0)
                / model("d4", 128.0)
        };
        handover.insert(a.into(), json!({"modelled_below_d4_at_2^20": below_now,
            "log2_E_where_order_flips": first, "model_ratio_at_2^256": model(a, 256.0) / model("d4", 256.0),
            "model_ratio_at_2^40": model(a, 40.0) / model("d4", 40.0),
            "model_ratio_at_2^64": model(a, 64.0) / model("d4", 64.0),
            "model_ratio_at_2^128": model(a, 128.0) / model("d4", 128.0),
            "la_share_d4_at_2^128": la_share}));
    }
    let summary = json!({"rows": table, "fits": fits, "paired_vs_d4": paired, "phase_model_handover": handover,
        "every_logarithm_verified": verified, "null_identical": nulls,
        "incumbent_trials_reproduce_item1": {"equal": repro.0, "of": repro.1},
        "note": "fits: OLS of log2(cost) on log2 #E over test curves, 95% bootstrap over curves (2000). Paired ratios: geometric means, bootstrap over curves. Handover figures are extrapolations."});
    std::fs::write(&out_path, serde_json::to_string_pretty(&summary).unwrap()).unwrap();
    let brief = json!({"verified": verified, "null": nulls, "repro": summary["incumbent_trials_reproduce_item1"],
        "fits": summary["fits"], "handover": summary["phase_model_handover"],
        "pooled": paired.iter().map(|(k, v)| (k.clone(), json!({
            "whole": [v["whole"]["pooled"].clone(), v["whole"]["constant_verdict"].clone(), v["whole"]["trend_slope_log2_ratio_per_log2_E"]["slope"].clone(), v["whole"]["trend_slope_log2_ratio_per_log2_E"]["ci95"].clone(), v["whole"]["exponent_verdict"].clone()],
            "relation": [v["relation"]["pooled"].clone(), v["relation"]["constant_verdict"].clone(), v["relation"]["exponent_verdict"].clone()],
            "la": [v["la"]["pooled"].clone(), v["la"]["constant_verdict"].clone(), v["la"]["exponent_verdict"].clone()]}))).collect::<serde_json::Map<_, _>>()});
    println!("{}", serde_json::to_string_pretty(&brief).unwrap());
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("train") => train(&args),
        Some("test") => test(&args),
        Some("tune") => tune(&args),
        Some("analyze") => analyze(&args),
        _ => eprintln!("usage: train|test|tune|analyze …"),
    }
}
