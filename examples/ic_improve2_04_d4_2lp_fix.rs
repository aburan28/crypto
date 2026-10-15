//! IC improvement round 2, item 04: fixing the `D₄` + double-large-prime
//! regression for `k = 4` index calculus on `E(F_{p⁴})`.  Protocol:
//! `research/ic_improve_20261015/04-d4-2lp-fix/PROTOCOL.md`.
//!
//! ```text
//! ic_improve2_04_d4_2lp_fix train   --p 521 --seed 4 --json raw/train.jsonl
//! ic_improve2_04_d4_2lp_fix fit     --train raw/train.jsonl --prior-train <item-01 train.jsonl> --out raw/frozen.json
//! ic_improve2_04_d4_2lp_fix test    --p 521 --seed 1 --frozen raw/frozen.json --json raw/test.jsonl
//! ic_improve2_04_d4_2lp_fix analyze --in raw/test.jsonl --prior <item-01 test.jsonl> --prior-summary <item-01 summary.json> --out raw/summary.json
//! ```
//!
//! Item 01's oracle, Y-line base, trial stream, charge and arms are copied
//! unchanged (`examples/ic_improve_01_d4_2lp.rs`).  Two changes:
//! the oracle keeps 48-bit bucketed fingerprints (output-identical, half the
//! memory, so p = 1601 fits in 3 GB), and arms gain the eager
//! first-dependency stop (`d4-ed`), all-rows feeding and the adaptive chooser.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;

use crypto_lib::cryptanalysis::gaudry_cubic::{wiedemann_u64, LargePrimeEliminator, SparseRel};
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
/// Margin before the adaptive arm leaves plain `D₄` (PROTOCOL.md).
const ADAPT_MARGIN: f64 = 0.95;

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

const NB_BITS: u32 = 16;
const NB: usize = 1 << NB_BITS;

/// 48 bits of the fingerprint: bucket (top 16) and value (next 32).
fn fkey(x: &E4) -> (usize, u32) {
    let f = fp(x);
    ((f >> (64 - NB_BITS)) as usize, (f >> 16) as u32)
}

/// The exact 4-term decomposition oracle over one base (item 01's), with the
/// triple table stored as 48-bit bucketed fingerprints.  A fingerprint hit
/// only triggers the exact pair lookup, so outputs are identical.
struct Oracle {
    pts: Vec<Pt4>,
    pairs: HashMap<E4, Vec<(u32, u32, i8, E4)>>,
    tri: Vec<u32>,
    cnt: Vec<u32>,
    bcap: usize,
    setup_ops: u64,
}

impl Oracle {
    fn build(curve: &Curve4, pts: Vec<Pt4>) -> Oracle {
        let before = curve.ops();
        let m = pts.len();
        let mut pairs: HashMap<E4, Vec<(u32, u32, i8, E4)>> = HashMap::new();
        let n = 4.0 * m as f64 * (m.saturating_sub(1)) as f64 * (m.saturating_sub(2)) as f64 / 6.0;
        let mean = n / NB as f64;
        // D₄ makes x-values repeat up to 4× (even sign flips by T), so the
        // bucket counts are clumped: allow 16 clumped standard deviations.
        let bcap = (mean + 16.0 * (4.0 * mean).sqrt() + 256.0).ceil() as usize;
        let mut tri = vec![0u32; NB * bcap];
        let mut cnt = vec![0u32; NB];
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
                                let (b, v) = fkey(&t.x);
                                let c = cnt[b] as usize;
                                assert!(c < bcap, "oracle bucket overflow");
                                tri[b * bcap + c] = v;
                                cnt[b] += 1;
                            }
                        }
                    }
                }
            }
        }
        for b in 0..NB {
            tri[b * bcap..b * bcap + cnt[b] as usize].sort_unstable();
        }
        Oracle {
            pts,
            pairs,
            tri,
            cnt,
            bcap,
            setup_ops: curve.ops() - before,
        }
    }

    fn has_triple(&self, x: &E4) -> bool {
        let (b, v) = fkey(x);
        self.tri[b * self.bcap..b * self.bcap + self.cnt[b] as usize]
            .binary_search(&v)
            .is_ok()
    }

    /// Every `R = Σ s_t P_{i_t}` over four distinct base points, signs
    /// included, sorted by index.
    fn decompose(&self, curve: &Curve4, r: &Pt4) -> Vec<Vec<(usize, i64)>> {
        let m = self.pts.len();
        let mut out: BTreeSet<Vec<(usize, i64)>> = BTreeSet::new();
        for l in 0..m {
            for sl in [1i64, -1] {
                let v = curve.sub(r, &signed(curve, &self.pts[l], sl));
                if v.inf || !self.has_triple(&v.x) {
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

/// `gaudry_cubic::square_core`, copied line for line, counting the index
/// operations of the filter (uncharged bookkeeping in every arm).
fn square_core_counted(
    rels: &[SparseRel],
    keep: usize,
    ops: &mut u64,
) -> Option<(Vec<usize>, Vec<bool>)> {
    let ncols = rels
        .iter()
        .flat_map(|r| r.cols.iter().map(|&(c, _)| c + 1))
        .max()
        .unwrap_or(0)
        .max(keep + 1);
    let mut alive: Vec<bool> = vec![true; rels.len()];
    let weights = |alive: &[bool], ops: &mut u64| -> Vec<usize> {
        let mut w = vec![0usize; ncols];
        for (i, r) in rels.iter().enumerate() {
            if alive[i] {
                *ops += r.cols.len() as u64;
                for &(c, _) in &r.cols {
                    w[c] += 1;
                }
            }
        }
        w
    };
    loop {
        loop {
            let w = weights(&alive, ops);
            let mut dropped = false;
            for (i, r) in rels.iter().enumerate() {
                if alive[i] {
                    *ops += r.cols.len() as u64;
                    if r.cols.iter().any(|&(c, _)| c != keep && w[c] == 1) {
                        alive[i] = false;
                        dropped = true;
                    }
                }
            }
            if !dropped {
                break;
            }
        }
        let w = weights(&alive, ops);
        let n_rows = alive.iter().filter(|&&a| a).count();
        let n_cols = w.iter().filter(|&&x| x > 0).count();
        if n_rows == 0 || w[keep] == 0 {
            return None;
        }
        if n_rows == n_cols {
            return Some((
                (0..rels.len()).filter(|&i| alive[i]).collect(),
                (0..ncols).map(|c| w[c] > 0).collect(),
            ));
        }
        if n_rows < n_cols {
            return None;
        }
        let victim = (0..rels.len()).filter(|&i| alive[i]).min_by_key(|&i| {
            rels[i]
                .cols
                .iter()
                .filter(|&&(c, _)| c != keep && w[c] == 2)
                .count()
        })?;
        alive[victim] = false;
    }
}

fn mulm(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 * b as u128) % r as u128) as u64
}

fn powm(mut a: u64, mut e: u64, r: u64) -> u64 {
    let mut x = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            x = mulm(x, a, r);
        }
        a = mulm(a, a, r);
        e >>= 1;
    }
    x
}

/// Hindsight diagnostic: the number of rows after which the base parts
/// (columns `< small`) first become linearly dependent mod `r`.
fn first_dependency(rows: &[SparseRel], small: usize, r: u64) -> Option<usize> {
    let mut piv: Vec<Option<Vec<u64>>> = vec![None; small];
    for (k, row) in rows.iter().enumerate() {
        let mut v = vec![0u64; small];
        for &(c, x) in &row.cols {
            if c < small {
                v[c] = x;
            }
        }
        let mut lead = None;
        for c in 0..small {
            if v[c] == 0 {
                continue;
            }
            match &piv[c] {
                Some(p) => {
                    let f = v[c];
                    for j in c..small {
                        if p[j] != 0 {
                            v[j] = (v[j] + r - mulm(f, p[j], r)) % r;
                        }
                    }
                }
                None => {
                    lead = Some(c);
                    break;
                }
            }
        }
        match lead {
            None => return Some(k + 1),
            Some(c) => {
                // Finish reducing the tail is unnecessary for rank; normalise.
                let inv = powm(v[c], r - 2, r);
                for x in v.iter_mut() {
                    *x = mulm(*x, inv, r);
                }
                piv[c] = Some(v);
            }
        }
    }
    None
}

/// One arm's configuration.
#[derive(Clone, Debug)]
struct Cfg {
    name: String,
    /// `None`: no large primes.  `Some(α)`: `s = ⌈ncols^α⌉`;
    /// `Some(f64::INFINITY)` is the null (`s = ncols`).
    alpha: Option<f64>,
    cap: u8,
    minlp: bool,
    omega: f64,
    stride: u64,
    /// Attempt the filter after every new full relation (first-dependency stop).
    eager: bool,
    /// Feed every distinct row of a decomposed trial.
    all_rows: bool,
    note: Value,
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
            eager: false,
            all_rows: false,
            note: Value::Null,
        }
    }
    fn lp(name: String, alpha: f64, cap: u8, minlp: bool, omega: f64) -> Cfg {
        Cfg {
            name,
            alpha: Some(alpha),
            cap,
            minlp,
            omega,
            ..Cfg::incumbent()
        }
    }
    fn named(mut self, n: &str) -> Cfg {
        self.name = n.into();
        self
    }
    fn json(&self) -> Value {
        json!({"alpha": self.alpha.map(|a| if a.is_finite() { json!(a) } else { json!("ncols") }),
               "merge_cap": self.cap, "pick": if self.minlp {"minlp"} else {"uniform"},
               "omega": self.omega, "stride": self.stride, "eager": self.eager,
               "all_rows": self.all_rows, "note": self.note})
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
    filter_ops: u64,
    filter_calls: u64,
    la_attempts: u64,
    targets_used: u64,
    rows_fed: u64,
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
        let first = if cfg.eager {
            1
        } else {
            ((cfg.omega * (small + 1) as f64).ceil() as usize).max(small + 1)
        };
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
            filter_ops: 0,
            filter_calls: 0,
            la_attempts: 0,
            targets_used: 0,
            rows_fed: 0,
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

    fn step(&self, base: usize) -> usize {
        if self.cfg.eager {
            1
        } else {
            (base / 20).max(1)
        }
    }

    /// Filter and solve if the schedule says so; true once solved.
    fn try_solve(&mut self, inst: &FoldInstance4, curve: &Curve4) -> bool {
        let r = inst.n;
        if self.full.len() < self.next_attempt {
            return false;
        }
        let unknowns = self.small + 1;
        self.filter_calls += 1;
        let Some((rows, columns)) =
            square_core_counted(&self.full, self.small, &mut self.filter_ops)
        else {
            self.next_attempt = self.full.len() + self.step(unknowns);
            return false;
        };
        self.la_attempts += 1;
        let mut map = vec![usize::MAX; unknowns];
        let mut k = 0;
        for (c, slot) in map.iter_mut().enumerate() {
            if columns.get(c).copied().unwrap_or(false) {
                *slot = k;
                k += 1;
            }
        }
        let sel: Vec<SparseRel> = rows
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
        self.next_attempt = self.full.len() + self.step(sel.len());
        false
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
        let idx = (pick % rows.len() as u64) as usize;
        let mut chosen = vec![rows[idx].clone()];
        if self.cfg.all_rows {
            for (i, x) in rows.iter().enumerate() {
                if i != idx && !chosen.iter().any(|c| c.0.cols == x.0.cols) {
                    chosen.push(x.clone());
                }
            }
        }
        self.targets_used += 1;
        for (rel, lp) in chosen {
            self.rows_fed += 1;
            self.lp_hist[lp.min(4)] += 1;
            let full = if lp == 0 {
                Some(rel)
            } else {
                self.elim
                    .as_mut()
                    .expect("large primes need an eliminator")
                    .feed(rel, &mut self.merge_ops)
            };
            let Some(full) = full else { continue };
            self.full.push(full);
            if self.try_solve(inst, curve) {
                return true;
            }
        }
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
    // One write per row, so two lanes appending to one file never interleave.
    f.write_all(format!("{v}\n").as_bytes()).unwrap();
    println!("{v}");
}

const ALPHAS: [f64; 5] = [0.75, 0.80, 0.85, 0.90, 0.95];

fn alpha_key(a: f64) -> String {
    format!("a{a:.2}-c0-uniform-w1.0")
}

fn null_cfg() -> Cfg {
    Cfg::lp("d4-2lp-null".into(), f64::INFINITY, 0, false, 1.0)
}

fn ed_cfg() -> Cfg {
    let mut c = Cfg::incumbent().named("d4-ed");
    c.eager = true;
    c
}

fn ed_all_cfg() -> Cfg {
    let mut c = ed_cfg().named("d4-ed-all");
    c.all_rows = true;
    c
}

fn late_cfg() -> Cfg {
    let mut c = Cfg::incumbent().named("d4-late110");
    c.omega = 1.1;
    c
}

/// One curve: build the Y-line oracle, run every arm on the shared stream.
fn run_curve(args: &[String], cfgs: Vec<Cfg>, split: &str, adapt: Option<&Value>) {
    let p: u64 = arg(args, "--p", 269);
    let seed: u64 = arg(args, "--seed", 1);
    let rho_runs: u64 = arg(args, "--rho-runs", 0);
    let json_path: String = arg(args, "--json", "out.jsonl".to_string());
    let t0 = Instant::now();
    let inst = generate_translation_instance(p, seed, 16).expect("instance");
    let r = inst.n;
    let sqrt_r = (r as f64).sqrt();
    let log2_e = (inst.group_order as f64).log2();
    let mut cfgs = cfgs;
    if let Some(model) = adapt {
        let (choice, ratio) = adapt_choice(model, log2_e);
        let mut c = match choice {
            None => Cfg::incumbent(),
            Some(a) => Cfg::lp(String::new(), a, 0, false, 1.0),
        };
        c.name = "d4-adapt".into();
        c.note =
            json!({"chosen_alpha": choice, "modelled_ratio_vs_plain": ratio, "log2_E": log2_e});
        cfgs.push(c);
    }
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
    let mut decomposed_trials = 0u64;
    let mut multi_row_trials = 0u64;
    while trials < cap && arms.iter().any(|a| !a.solved) {
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
        if !good.is_empty() {
            decomposed_trials += 1;
            let distinct: BTreeSet<Vec<(usize, u64)>> = good
                .iter()
                .map(|d| arms[0].row(d, a, b, r).0.cols)
                .collect();
            if distinct.len() > 1 {
                multi_row_trials += 1;
            }
        }
        for arm in arms.iter_mut().filter(|a| !a.solved) {
            if arm.feed(t, &good, a, b, &inst, curve, pick) {
                arm.snapshot.insert("trials", trials);
                arm.snapshot.insert("scalar_ops", scalar_ops);
                arm.snapshot.insert("verify_ops", verify_ops);
            }
        }
    }
    let check = 1.5 * (r as f64).log2();
    // Hindsight first dependency of the incumbent's row sequence.
    let fd = arms
        .iter()
        .find(|a| a.cfg.name == "d4")
        .and_then(|a| first_dependency(&a.full, a.small, r));
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
                "rows_fed": arm.rows_fed,
                "lp_hist": arm.lp_hist, "la_attempts": arm.la_attempts,
                "filter_calls": arm.filter_calls, "filter_index_ops_uncharged": arm.filter_ops,
                "la_modr_mults": arm.la_ops, "merge_modr_mults": arm.merge_ops,
                "pivots": arm.elim.as_ref().map(|e| e.pivot_count()),
                "abandoned": arm.elim.as_ref().map(|e| e.abandoned),
                "gae": {"factor_base": fb_gae, "trial_scalar": scalar, "solver_charged": solver,
                        "verify": ver, "relation": rel, "la": la, "final_check": check, "whole": whole},
                "S_charged": whole / sqrt_r,
            })
        })
        .collect();
    let mut row = json!({
        "item": "ic_improve2_04_d4_2lp_fix", "split": split,
        "p": p, "seed": seed, "r": r, "log2_r": (r as f64).log2(),
        "group_order": inst.group_order, "log2_E": log2_e,
        "cofactor": inst.cofactor, "planted_d": inst.d, "c_add": C_ADD, "c_modn": C_MODN, "c4": C4_D4,
        "y_base": yb.len(), "y_columns": ncols, "oracle_setup_ops": oracle.setup_ops,
        "oracle_build_s": build_s, "stream_trials": trials, "stream_decompositions": decs_total,
        "stream_decomposed_trials": decomposed_trials, "stream_multi_row_trials": multi_row_trials,
        "d4_first_dependency_rows_hindsight": fd, "arms": out,
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
    let cfgs = vec![
        Cfg::incumbent(),
        ed_cfg(),
        ed_all_cfg(),
        null_cfg(),
        late_cfg(),
        Cfg::lp(alpha_key(0.95), 0.95, 0, false, 1.0),
    ];
    run_curve(args, cfgs, "train", None);
}

fn test(args: &[String]) {
    let frozen: String = arg(args, "--frozen", "raw/frozen.json".to_string());
    let model: Value =
        serde_json::from_str(&std::fs::read_to_string(&frozen).expect("frozen")).unwrap();
    let cfgs = vec![
        Cfg::incumbent(),
        ed_cfg(),
        ed_all_cfg(),
        Cfg::lp("d4-2lp".into(), 0.95, 0, false, 1.0),
        null_cfg(),
        late_cfg(),
    ];
    run_curve(args, cfgs, "test", Some(&model["adaptive"]));
}

// ── Fit (Rust; Python only plots).

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

fn ols(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
    let b = sxy / sxx;
    (b, my - b * mx)
}

/// Modelled whole cost of one configuration at `x = log₂ #E`.
fn model_cost(m: &Value, x: f64) -> f64 {
    ["relation", "la"]
        .iter()
        .map(|ph| 2f64.powf(m[*ph][1].as_f64().unwrap() + m[*ph][0].as_f64().unwrap() * x))
        .sum()
}

/// The adaptive rule: `(None | Some(α), modelled ratio vs plain)`.
fn adapt_choice(model: &Value, x: f64) -> (Option<f64>, f64) {
    let plain = model_cost(&model["plain"], x);
    let mut best: (Option<f64>, f64) = (None, 1.0);
    for a in ALPHAS {
        let c = model_cost(&model[alpha_key(a)], x) / plain;
        if c < best.1 {
            best = (Some(a), c);
        }
    }
    if best.1 < ADAPT_MARGIN {
        best
    } else {
        (None, best.1)
    }
}

fn fit(args: &[String]) {
    let train_path: String = arg(args, "--train", "raw/train.jsonl".to_string());
    let prior_path: String = arg(args, "--prior-train", String::new());
    let out_path: String = arg(args, "--out", "raw/frozen.json".to_string());
    let mine = load(&train_path);
    let prior = load(&prior_path);
    let key = |r: &Value| (r["p"].as_u64().unwrap(), r["seed"].as_u64().unwrap());
    // Reproduction: d4 and α = 0.95 trial counts on the same curves.
    let mut repro = Vec::new();
    let mut ok = true;
    for r in &mine {
        let pr = prior
            .iter()
            .find(|q| key(q) == key(r))
            .expect("prior training curve");
        for (a, b) in [
            ("d4", "d4"),
            (alpha_key(0.95).as_str(), alpha_key(0.95).as_str()),
        ] {
            let x = arm_of(r, a)["trials"].as_f64().unwrap();
            let y = arm_of(pr, b)["trials"].as_f64().unwrap();
            ok &= x == y;
            repro.push(
                json!({"p": key(r).0, "seed": key(r).1, "arm": a, "this_run": x, "item01": y}),
            );
        }
    }
    let src = if ok { &prior } else { &mine };
    let xs: Vec<f64> = src.iter().map(|r| r["log2_E"].as_f64().unwrap()).collect();
    let mut model = serde_json::Map::new();
    let mut names: Vec<(String, String)> = vec![("plain".into(), "d4".into())];
    for a in ALPHAS {
        names.push((alpha_key(a), alpha_key(a)));
    }
    for (k, arm) in &names {
        let mut m = serde_json::Map::new();
        for ph in ["relation", "la"] {
            let ys: Vec<f64> = src
                .iter()
                .map(|r| arm_of(r, arm)["gae"][ph].as_f64().unwrap().log2())
                .collect();
            let (b, a) = ols(&xs, &ys);
            m.insert(ph.into(), json!([b, a]));
        }
        model.insert(k.clone(), Value::Object(m));
    }
    let model = Value::Object(model);
    // Where the frozen rule first leaves plain D₄ (extrapolation).
    // Every flip of the rule's choice from 20 to 600 (below the smallest
    // training size, 26.3, the fits are extrapolated too).
    let mut leave = Vec::new();
    let mut prev: Option<f64> = None;
    let mut x = 20.0;
    while x <= 600.0 {
        let (ch, c) = adapt_choice(&model, x);
        if ch != prev {
            leave.push(
                json!({"log2_E": (x * 10.0f64).round() / 10.0, "alpha": ch, "modelled_ratio": c}),
            );
            prev = ch;
        }
        x += 0.1;
    }
    let at: Vec<Value> = [26.0, 30.0, 35.0, 40.0, 41.4, 42.6, 64.0, 128.0, 256.0]
        .iter()
        .map(|&x| {
            let (c, r) = adapt_choice(&model, x);
            json!({"log2_E": x, "choice": c, "modelled_best_ratio_vs_plain": r})
        })
        .collect();
    let out = json!({
        "rule": "fit OLS of log2(relation GAE) and log2(LA GAE) on log2 #E per configuration over the 14 training curves; at a curve's log2 #E choose the alpha with the least modelled whole cost; leave plain D4 only if that is < 0.95 x plain",
        "source": if ok {"item-01 raw/train.jsonl (reproduced exactly by this run)"} else {"this run's train.jsonl (item-01 not reproduced)"},
        "training_reproduces_item01": ok, "reproduction": repro,
        "margin": ADAPT_MARGIN,
        "adaptive": model, "first_leaves_plain": leave, "choice_at": at,
    });
    std::fs::write(&out_path, serde_json::to_string_pretty(&out).unwrap()).unwrap();
    println!(
        "{}",
        serde_json::to_string_pretty(
            &json!({"repro_ok": ok, "leave": out["first_leaves_plain"], "at": out["choice_at"]})
        )
        .unwrap()
    );
    // Training sanity for the candidate (not a test result).
    for r in &mine {
        let i = arm_of(r, "d4");
        let e = arm_of(r, "d4-ed");
        let ea = arm_of(r, "d4-ed-all");
        println!(
            "train p={} s={} cols={} fd={} d4 rel={} tr={} | ed rel={} tr={} ok={} | ed-all tr={} ok={} | null={} late={}",
            r["p"], r["seed"], r["y_columns"], r["d4_first_dependency_rows_hindsight"],
            i["full_relations"], i["trials"], e["full_relations"], e["trials"], e["correct"],
            ea["trials"], ea["correct"],
            arm_of(r, "d4-2lp-null")["trials"] == i["trials"],
            arm_of(r, "d4-late110")["trials"].as_f64().unwrap() / i["trials"].as_f64().unwrap()
        );
    }
}

// ── Analysis.

fn quant(mut v: Vec<f64>, f: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() as f64 - 1.0) * f).round() as usize]
}

/// OLS of `ly` (already log2) on `xs` with a 95 % bootstrap over rows.
fn fit_log(xs: &[f64], ly: &[f64], rng: &mut StdRng) -> Value {
    let (b, a) = ols(xs, ly);
    let n = xs.len();
    let mut bs = Vec::new();
    for _ in 0..2000 {
        let idx: Vec<usize> = (0..n).map(|_| rng.gen_range(0..n)).collect();
        let x: Vec<f64> = idx.iter().map(|&i| xs[i]).collect();
        if x.iter().all(|v| (v - x[0]).abs() < 1e-9) {
            continue;
        }
        let y: Vec<f64> = idx.iter().map(|&i| ly[i]).collect();
        bs.push(ols(&x, &y).0);
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
    // Exact ties (ratio 1 on every curve) are MATCHES; 1e-12 absorbs the
    // rounding of exp2(mean log2).
    if hi < 1.0 - 1e-12 {
        "ADVANCES"
    } else if lo > 1.0 + 1e-12 {
        "REGRESSES"
    } else {
        "MATCHES"
    }
}

fn verdict_slope(v: &Value) -> &'static str {
    let lo = v["ci95"][0].as_f64().unwrap();
    let hi = v["ci95"][1].as_f64().unwrap();
    if hi < -1e-12 {
        "ADVANCES"
    } else if lo > 1e-12 {
        "REGRESSES"
    } else {
        "MATCHES"
    }
}

fn null_equal(row: &Value) -> bool {
    let (a, b) = (arm_of(row, "d4"), arm_of(row, "d4-2lp-null"));
    a["trials"] == b["trials"]
        && a["full_relations"] == b["full_relations"]
        && a["la_modr_mults"] == b["la_modr_mults"]
        && a["recovered"] == b["recovered"]
}

const ARMS: [&str; 7] = [
    "d4",
    "d4-ed",
    "d4-ed-all",
    "d4-adapt",
    "d4-2lp",
    "d4-2lp-null",
    "d4-late110",
];

fn analyze(args: &[String]) {
    let input: String = arg(args, "--in", "raw/test.jsonl".to_string());
    let prior: String = arg(args, "--prior", String::new());
    let prior_summary: String = arg(args, "--prior-summary", String::new());
    let out_path: String = arg(args, "--out", "raw/summary.json".to_string());
    let mut rows = load(&input);
    rows.sort_by_key(|r| (r["p"].as_u64().unwrap(), r["seed"].as_u64().unwrap()));
    let mut rng = StdRng::seed_from_u64(20261015);
    let key = |r: &Value| (r["p"].as_u64().unwrap(), r["seed"].as_u64().unwrap());
    let prior_rows = if prior.is_empty() {
        vec![]
    } else {
        load(&prior)
    };
    let prior_rho: HashMap<(u64, u64), f64> = if prior_summary.is_empty() {
        HashMap::new()
    } else {
        let s: Value =
            serde_json::from_str(&std::fs::read_to_string(&prior_summary).unwrap()).unwrap();
        s["rows"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|r| Some((key(r), r["rho_neg_S_reference"].as_f64()?)))
            .collect()
    };
    let phases = ["whole", "relation", "la"];
    let mut verified = true;
    let mut nulls = true;
    let mut repro = (0, 0);
    let mut dominance = Vec::new();
    let mut xs = Vec::new();
    let mut costs: BTreeMap<(&str, &str), Vec<f64>> = BTreeMap::new();
    let mut sens: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
    let mut table = Vec::new();
    for row in &rows {
        let l2e = row["log2_E"].as_f64().unwrap();
        xs.push(l2e);
        nulls &= null_equal(row);
        let k = key(row);
        let inc = arm_of(row, "d4");
        if let Some(pr) = prior_rows.iter().find(|q| key(q) == k) {
            repro.1 += 1;
            if arm_of(pr, "d4")["trials"] == inc["trials"] {
                repro.0 += 1;
            }
        }
        let lp_tr = arm_of(row, "d4-2lp")["trials"].as_f64().unwrap();
        for a in ["d4-ed", "d4-ed-all"] {
            let t = arm_of(row, a)["trials"].as_f64().unwrap();
            if t > lp_tr {
                dominance.push(
                    json!({"p": k.0, "seed": k.1, "arm": a, "trials": t, "d4_2lp_trials": lp_tr}),
                );
            }
        }
        let rho = row["rho_neg"]["mean_S"]
            .as_f64()
            .or(prior_rho.get(&k).copied());
        let mut cells = serde_json::Map::new();
        let inc_whole = inc["gae"]["whole"].as_f64().unwrap();
        let inc_sens = inc_whole + inc["filter_index_ops_uncharged"].as_f64().unwrap() / C_ADD;
        for a in ARMS {
            let v = arm_of(row, a);
            verified &= v["solved"] == true && v["correct"] == true;
            for ph in phases {
                costs
                    .entry((a, ph))
                    .or_default()
                    .push(v["gae"][ph].as_f64().unwrap());
            }
            let s = v["gae"]["whole"].as_f64().unwrap()
                + v["filter_index_ops_uncharged"].as_f64().unwrap() / C_ADD;
            sens.entry(a).or_default().push((s / inc_sens).log2());
            cells.insert(a.into(), json!({"S_charged": v["S_charged"], "trials": v["trials"],
                "full_relations": v["full_relations"], "rows_fed": v["rows_fed"],
                "small_columns": v["small_columns"], "unknowns_core": v["unknowns_core"],
                "la_attempts": v["la_attempts"], "filter_calls": v["filter_calls"],
                "filter_index_ops_uncharged": v["filter_index_ops_uncharged"],
                "lp_hist": v["lp_hist"], "relation_gae": v["gae"]["relation"], "la_gae": v["gae"]["la"],
                "whole_gae": v["gae"]["whole"], "solved": v["solved"], "correct": v["correct"],
                "cfg": v["cfg"],
                "ratio_vs_d4": v["gae"]["whole"].as_f64().unwrap() / inc_whole,
                "S_over_rho_neg_reference": rho.map(|x| v["S_charged"].as_f64().unwrap() / x)}));
        }
        let cols = row["y_columns"].as_f64().unwrap();
        table.push(json!({"p": k.0, "seed": k.1, "log2_E": l2e, "log2_r": row["log2_r"],
            "y_columns": row["y_columns"], "rho_neg_S_reference": rho,
            "d4_first_dependency_rows_hindsight": row["d4_first_dependency_rows_hindsight"],
            "fd_over_cols_plus_1": row["d4_first_dependency_rows_hindsight"].as_f64().map(|f| f / (cols + 1.0)),
            "stream_decomposed_trials": row["stream_decomposed_trials"],
            "stream_multi_row_trials": row["stream_multi_row_trials"],
            "wall_s": row["wall_s"], "oracle_build_s": row["oracle_build_s"], "arms": cells}));
    }
    let mut fits = serde_json::Map::new();
    for a in ARMS {
        let mut m = serde_json::Map::new();
        for ph in phases {
            let ly: Vec<f64> = costs[&(a, ph)].iter().map(|c| c.log2()).collect();
            m.insert(ph.into(), fit_log(&xs, &ly, &mut rng));
        }
        fits.insert(a.into(), Value::Object(m));
    }
    let mut ps: Vec<u64> = rows.iter().map(|r| r["p"].as_u64().unwrap()).collect();
    ps.dedup();
    let mut paired = serde_json::Map::new();
    for a in &ARMS[1..] {
        let mut m = serde_json::Map::new();
        for ph in phases {
            let lr: Vec<f64> = costs[&(*a, ph)]
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
                    c["verdict"] = json!(verdict_ratio(&c));
                    c["p"] = json!(p);
                    c["log2_E"] = json!(l2e);
                    c["per_curve"] = json!(v.iter().map(|x| 2f64.powf(*x)).collect::<Vec<_>>());
                    c
                })
                .collect();
            m.insert(ph.into(), json!({
                "pooled": pooled, "constant_verdict": verdict_ratio(&pooled),
                "trend_slope_log2_ratio_per_log2_E": trend, "exponent_verdict": verdict_slope(&trend),
                "per_size": per_size}));
        }
        let sp = gm_ci(&sens[a], &mut rng);
        m.insert(
            "whole_with_filter_ops_charged".into(),
            json!({"pooled": sp, "constant_verdict": verdict_ratio(&sp)}),
        );
        paired.insert((*a).into(), Value::Object(m));
    }
    // Relation-count accounting.
    let rel_ratio = |a: &str| -> Vec<f64> {
        rows.iter()
            .map(|r| {
                arm_of(r, a)["full_relations"].as_f64().unwrap()
                    / (r["y_columns"].as_f64().unwrap() + 1.0)
            })
            .collect()
    };
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let fdr: Vec<f64> = rows
        .iter()
        .filter_map(|r| {
            Some(
                r["d4_first_dependency_rows_hindsight"].as_f64()?
                    / (r["y_columns"].as_f64().unwrap() + 1.0),
            )
        })
        .collect();
    let ed_gap: Vec<f64> = rows
        .iter()
        .filter_map(|r| {
            Some(
                arm_of(r, "d4-ed")["full_relations"].as_f64().unwrap()
                    - r["d4_first_dependency_rows_hindsight"].as_f64()?,
            )
        })
        .collect();
    let accounting = json!({
        "mean_full_relations_over_cols_plus_1": {"d4": mean(&rel_ratio("d4")), "d4-ed": mean(&rel_ratio("d4-ed")),
            "d4-ed-all": mean(&rel_ratio("d4-ed-all")), "d4-2lp": mean(&rel_ratio("d4-2lp")), "d4-late110": mean(&rel_ratio("d4-late110"))},
        "mean_first_dependency_over_cols_plus_1": mean(&fdr),
        "d4_ed_relations_minus_first_dependency": {"mean": mean(&ed_gap), "max": ed_gap.iter().cloned().fold(f64::MIN, f64::max), "min": ed_gap.iter().cloned().fold(f64::MAX, f64::min)},
        "multi_row_trials_total": rows.iter().map(|r| r["stream_multi_row_trials"].as_u64().unwrap()).sum::<u64>(),
        "decomposed_trials_total": rows.iter().map(|r| r["stream_decomposed_trials"].as_u64().unwrap()).sum::<u64>(),
    });
    let summary = json!({"rows": table, "fits": fits, "paired_vs_d4": paired, "accounting": accounting,
        "every_logarithm_verified": verified, "null_identical": nulls,
        "incumbent_trials_reproduce_item01": {"equal": repro.0, "of": repro.1},
        "dominance_violations_ed_vs_2lp": dominance,
        "note": "fits: OLS of log2(cost) on log2 #E over test curves, 95% bootstrap over curves (2000). Paired ratios: geometric means, bootstrap over curves."});
    std::fs::write(&out_path, serde_json::to_string_pretty(&summary).unwrap()).unwrap();
    let brief = json!({"verified": verified, "null": nulls, "repro": summary["incumbent_trials_reproduce_item01"],
        "dominance_violations": summary["dominance_violations_ed_vs_2lp"], "accounting": accounting,
        "pooled": paired.iter().map(|(k, v)| (k.clone(), json!({
            "whole": [v["whole"]["pooled"].clone(), v["whole"]["constant_verdict"].clone(), v["whole"]["trend_slope_log2_ratio_per_log2_E"]["slope"].clone(), v["whole"]["trend_slope_log2_ratio_per_log2_E"]["ci95"].clone(), v["whole"]["exponent_verdict"].clone()],
            "relation": [v["relation"]["pooled"].clone(), v["relation"]["constant_verdict"].clone(), v["relation"]["exponent_verdict"].clone()],
            "la": [v["la"]["pooled"].clone(), v["la"]["constant_verdict"].clone(), v["la"]["trend_slope_log2_ratio_per_log2_E"]["slope"].clone(), v["la"]["exponent_verdict"].clone()],
            "per_size_whole": v["whole"]["per_size"].as_array().unwrap().iter().map(|c| json!([c["p"], c["ratio"], c["ci95"], c["verdict"]])).collect::<Vec<_>>(),
            "sens": v["whole_with_filter_ops_charged"].clone()}))).collect::<serde_json::Map<_, _>>()});
    println!("{}", serde_json::to_string_pretty(&brief).unwrap());
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("train") => train(&args),
        Some("test") => test(&args),
        Some("fit") => fit(&args),
        Some("analyze") => analyze(&args),
        _ => eprintln!("usage: train|fit|test|analyze …"),
    }
}
