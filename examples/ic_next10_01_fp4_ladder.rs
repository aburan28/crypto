//! IC next-10, item 1: the `E(F_{p⁴})` ladder end to end — plain `k = 4`,
//! `D₄`-folded (Y-line), double-large-prime, and the matched rho walks, on
//! one curve and one planted target per row, all in ecbench's unit
//! (group-addition equivalents, GAE; one `E(F_{p⁴})` addition = 97 `F_p`
//! multiplications).  Protocol:
//! `research/ic_next10_20261014/01-fp4-ladder/PROTOCOL.md`.
//!
//! ```text
//! ic_next10_01_fp4_ladder ladder  --p 269 --seed 1 --rho-runs 16 --json out.jsonl
//! ic_next10_01_fp4_ladder control --p 269 --seed 1 --planted 3 --uniform 3 --json ctl.jsonl
//! ic_next10_01_fp4_ladder census  --p 101 --curves 200 --json census.jsonl
//! ```
//!
//! Relations come from an exact enumeration oracle (a table of
//! `x(±P_i ± P_j ± P_k)` probed with `R − s·P_l`, hits recovered through the
//! pair table); the relation phase is reported both as executed and charged
//! at the algebraic solver's per-trial cost `C₄` (PROTOCOL.md, "Oracle and
//! charge").  `control` runs the real `D₄` solver against the oracle.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;

use crypto_lib::cryptanalysis::gaudry_cubic::{
    square_core, wiedemann_u64, LargePrimeEliminator, SparseRel,
};
use crypto_lib::cryptanalysis::gaudry_quartic::{Curve4, Pt4, E4};
use crypto_lib::cryptanalysis::gaudry_quartic_d4::{solve_d4, y_base, D4Instance, SymmetrisedD4};
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
/// §11.17's measured `C₄` for the `S₄` (x-line) solve.
const C4_S4: f64 = 1.213e12;
/// `gaudry_quartic_d4_20261013`'s measured `C₄` for the `D₄` solve
/// (mid of 2.535–2.552·10¹⁰); replaced by the in-situ control when given.
const C4_D4_COMMITTED: f64 = 2.54e10;

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

/// One IC arm's state on the shared trial stream.
struct Arm {
    name: &'static str,
    /// Base index → (column, coefficient ±1 as ±1).
    col: Vec<(usize, i64)>,
    /// Columns `< small` are the small base; `d` is column `small`; large
    /// primes are `> small`.
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
    trials: u64,
    unknowns: usize,
    filtered_out: usize,
    snapshot: BTreeMap<&'static str, u64>,
}

impl Arm {
    fn new(name: &'static str, col: Vec<(usize, i64)>, small: usize, max_lp: usize, r: u64) -> Arm {
        Arm {
            name,
            col,
            small,
            max_lp,
            elim: (max_lp > 0).then(|| LargePrimeEliminator::new(small + 1, r)),
            full: Vec::new(),
            next_attempt: small + 1,
            la_ops: 0,
            merge_ops: 0,
            la_attempts: 0,
            targets_used: 0,
            lp_hist: [0; 5],
            solved: false,
            correct: false,
            trials: 0,
            unknowns: 0,
            filtered_out: 0,
            snapshot: BTreeMap::new(),
        }
    }

    /// The relation row of a verified decomposition of `aG + bQ`.
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

    /// Feed one target's verified decompositions; true once solved.
    #[allow(clippy::too_many_arguments)]
    fn feed(
        &mut self,
        decs: &[Vec<(usize, i64)>],
        a: u64,
        b: u64,
        inst: &FoldInstance4,
        curve: &Curve4,
        pick: u64,
        la_rng: &mut StdRng,
    ) -> bool {
        let r = inst.n;
        self.trials += 1;
        let rows: Vec<(SparseRel, usize)> = decs
            .iter()
            .map(|t| self.row(t, a, b, r))
            .filter(|(_, lp)| *lp <= self.max_lp)
            .collect();
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
        let x = wiedemann_u64(&sel, sel.len(), r, la_rng, &mut self.la_ops);
        let dd = x.map(|x| x[map[self.small]]);
        if dd.is_some_and(|dd| curve.mul(&inst.g, dd) == inst.q) {
            self.solved = true;
            self.correct = dd == Some(inst.d);
            self.unknowns = sel.len();
            self.filtered_out = self.full.len() - sel.len();
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

struct Stream {
    scalar_ops: u64,
    oracle_ops: u64,
    verify_ops: u64,
    trials: u64,
    decs: u64,
    decomposed: u64,
}

/// Run arms on one base until all are solved or the cap is hit.
#[allow(clippy::too_many_arguments)]
fn run_stream(
    inst: &FoldInstance4,
    curve: &Curve4,
    oracle: &Oracle,
    arms: &mut [Arm],
    seed: u64,
    cap: u64,
    c4: f64,
    la_rng: &mut StdRng,
) -> (Stream, Vec<Value>) {
    let r = inst.n;
    let mut st = Stream {
        scalar_ops: 0,
        oracle_ops: 0,
        verify_ops: 0,
        trials: 0,
        decs: 0,
        decomposed: 0,
    };
    while st.trials < cap && arms.iter().any(|a| !a.solved) {
        let (a, b, pick) = trial_ab(seed, st.trials, r);
        st.trials += 1;
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
        st.scalar_ops += o1 - o0;
        st.oracle_ops += o2 - o1;
        st.verify_ops += o3 - o2;
        st.decs += good.len() as u64;
        st.decomposed += u64::from(!good.is_empty());
        for arm in arms.iter_mut().filter(|a| !a.solved) {
            if arm.feed(&good, a, b, inst, curve, pick, la_rng) {
                arm.snapshot.insert("trials", st.trials);
                arm.snapshot.insert("scalar_ops", st.scalar_ops);
                arm.snapshot.insert("oracle_ops", st.oracle_ops);
                arm.snapshot.insert("verify_ops", st.verify_ops);
                arm.snapshot.insert("decs", st.decs);
            }
        }
        // Unsolved arms see every trial; trials counted inside feed.
    }
    let sqrt_r = (r as f64).sqrt();
    let out = arms
        .iter()
        .map(|arm| {
            let g = |k: &str| *arm.snapshot.get(k).unwrap_or(&0) as f64;
            let trials = if arm.solved {
                g("trials")
            } else {
                st.trials as f64
            };
            let scalar = if arm.solved {
                g("scalar_ops")
            } else {
                st.scalar_ops as f64
            };
            let orc = if arm.solved {
                g("oracle_ops")
            } else {
                st.oracle_ops as f64
            };
            let ver = if arm.solved {
                g("verify_ops")
            } else {
                st.verify_ops as f64
            };
            let la_gae = (arm.la_ops + arm.merge_ops) as f64 * C_MODN / C_ADD;
            let solver_gae = trials * c4 / C_ADD;
            let check = 1.5 * (r as f64).log2();
            let charged = scalar + solver_gae + ver + la_gae + check;
            let executed = scalar + orc + ver + la_gae + check;
            json!({
                "arm": arm.name,
                "solved": arm.solved,
                "correct": arm.correct,
                "trials": trials,
                "targets_used": arm.targets_used,
                "small_columns": arm.small,
                "unknowns_core": arm.unknowns,
                "filtered_out": arm.filtered_out,
                "full_relations": arm.full.len(),
                "la_attempts": arm.la_attempts,
                "la_modr_mults": arm.la_ops,
                "merge_modr_mults": arm.merge_ops,
                "lp_hist": arm.lp_hist,
                "pivots": arm.elim.as_ref().map(|e| e.pivot_count()),
                "c4_fp_mults": c4,
                "gae": {
                    "trial_scalar": scalar,
                    "oracle_executed": orc,
                    "solver_charged": solver_gae,
                    "verify": ver,
                    "la": la_gae,
                    "final_check": check,
                },
                "S_charged_wo_fb": charged / sqrt_r,
                "S_executed_wo_fb": executed / sqrt_r,
            })
        })
        .collect();
    (st, out)
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

fn d4_of(inst: &FoldInstance4) -> D4Instance {
    D4Instance {
        curve: Curve4::with_params(inst.p, inst.a, inst.b),
        x0: inst.x0,
        c: inst.c,
    }
}

fn ladder(args: &[String]) {
    let p: u64 = arg(args, "--p", 269);
    let seed: u64 = arg(args, "--seed", 1);
    let rho_runs: u64 = arg(args, "--rho-runs", 16);
    let c4_d4: f64 = arg(args, "--c4-d4", C4_D4_COMMITTED);
    let json_path: String = arg(args, "--json", "ladder.jsonl".to_string());
    let only: String = arg(args, "--arms", "x,d4,rho".to_string());
    let t0 = Instant::now();
    let inst = generate_translation_instance(p, seed, 16).expect("instance");
    let r = inst.n;
    let sqrt_r = (r as f64).sqrt();
    let curve = inst.curve();
    let mut rng = StdRng::seed_from_u64(seed ^ 0x01F4_1AD0);
    let mut la_rng = StdRng::seed_from_u64(seed ^ 0x1A1A);
    let mut row = json!({
        "item": "ic_next10_01_fp4_ladder",
        "p": p, "seed": seed, "r": r, "log2_r": (r as f64).log2(),
        "group_order": inst.group_order, "log2_E": (inst.group_order as f64).log2(),
        "cofactor": inst.cofactor, "c_add": C_ADD, "c_modn": C_MODN,
    });
    // Measured F_p multiplications per addition on this field.
    {
        curve.f.reset_muls();
        let mut acc = inst.g;
        for _ in 0..64 {
            acc = curve.add(&acc, &inst.q);
        }
        row["fp_muls_per_add_measured"] = json!(curve.f.muls() as f64 / 64.0);
    }
    let mut arms_out: Vec<Value> = Vec::new();
    if only.contains('x') {
        curve.f.reset_muls();
        let base = x_line_base(&curve, &mut rng);
        let fb_gae = curve.f.muls() as f64 / C_ADD;
        let m = base.len();
        let s = ((m as f64).powf(0.75)).ceil() as usize;
        let tb = Instant::now();
        let oracle = Oracle::build(&curve, base);
        row["x_base"] = json!(m);
        row["x_oracle_setup_ops"] = json!(oracle.setup_ops);
        row["x_oracle_build_s"] = json!(tb.elapsed().as_secs_f64());
        row["x_fb_gae"] = json!(fb_gae);
        let ident: Vec<(usize, i64)> = (0..m).map(|i| (i, 1)).collect();
        let dlp_col: Vec<(usize, i64)> =
            (0..m).map(|i| (if i < s { i } else { i + 1 }, 1)).collect();
        let mut arms = vec![
            Arm::new("plain", ident.clone(), m, 0, r),
            Arm::new("dlp", dlp_col, s, 2, r),
            Arm::new("dlp-null", ident, m, 2, r),
        ];
        let cap = 600 * m as u64;
        let (st, mut out) = run_stream(
            &inst,
            &curve,
            &oracle,
            &mut arms,
            seed,
            cap,
            C4_S4,
            &mut la_rng,
        );
        for o in out.iter_mut() {
            o["fb_gae"] = json!(fb_gae);
            o["oracle_setup_gae"] = json!(oracle.setup_ops);
            let c = o["S_charged_wo_fb"].as_f64().unwrap() + fb_gae / sqrt_r;
            let e = o["S_executed_wo_fb"].as_f64().unwrap()
                + (fb_gae + oracle.setup_ops as f64) / sqrt_r;
            o["S_charged"] = json!(c);
            o["S_executed"] = json!(e);
        }
        row["x_stream"] = json!({"trials": st.trials, "decompositions": st.decs, "decomposed_targets": st.decomposed});
        arms_out.extend(out);
        drop(oracle);
    }
    if only.contains("d4") {
        let d4 = d4_of(&inst);
        let dcurve = &d4.curve;
        dcurve.f.reset_muls();
        let (yb, col, ncols) = y_line(&inst, &d4, &mut rng);
        let fb_gae = dcurve.f.muls() as f64 / C_ADD;
        let pts: Vec<Pt4> = yb.iter().map(|(_, q)| *q).collect();
        let tb = Instant::now();
        let oracle = Oracle::build(dcurve, pts);
        row["y_base"] = json!(yb.len());
        row["y_columns"] = json!(ncols);
        row["y_oracle_setup_ops"] = json!(oracle.setup_ops);
        row["y_oracle_build_s"] = json!(tb.elapsed().as_secs_f64());
        let mut arms = vec![Arm::new("d4", col, ncols, 0, r)];
        let cap = 2000 * ncols as u64;
        let (st, mut out) = run_stream(
            &inst,
            dcurve,
            &oracle,
            &mut arms,
            seed,
            cap,
            c4_d4,
            &mut la_rng,
        );
        for o in out.iter_mut() {
            o["fb_gae"] = json!(fb_gae);
            o["oracle_setup_gae"] = json!(oracle.setup_ops);
            let c = o["S_charged_wo_fb"].as_f64().unwrap() + fb_gae / sqrt_r;
            let e = o["S_executed_wo_fb"].as_f64().unwrap()
                + (fb_gae + oracle.setup_ops as f64) / sqrt_r;
            o["S_charged"] = json!(c);
            o["S_executed"] = json!(e);
        }
        row["y_stream"] = json!({"trials": st.trials, "decompositions": st.decs, "decomposed_targets": st.decomposed});
        arms_out.extend(out);
    }
    row["arms"] = json!(arms_out);
    if only.contains("rho") {
        let mut rho = serde_json::Map::new();
        for (name, fold) in [("rho", RhoFold::None), ("rho-neg", RhoFold::Negation)] {
            let runs: Vec<Value> = (0..rho_runs)
                .map(|k| {
                    let w = rho4_folded(&inst, fold, seed * 1000 + k);
                    json!({"seed": w.seed, "steps": w.steps, "group_ops": w.group_ops,
                           "fp_muls": w.fp_muls, "cycle_escapes": w.cycle_escapes,
                           "correct": w.correct, "S": w.group_ops as f64 / sqrt_r,
                           "S_fp": w.fp_muls as f64 / C_ADD / sqrt_r})
                })
                .collect();
            let ss: Vec<f64> = runs.iter().map(|v| v["S"].as_f64().unwrap()).collect();
            let mean = ss.iter().sum::<f64>() / ss.len() as f64;
            rho.insert(name.into(), json!({"mean_S": mean, "runs": runs}));
        }
        row["rho"] = Value::Object(rho);
    }
    row["wall_s"] = json!(t0.elapsed().as_secs_f64());
    emit(&json_path, &row);
}

/// The real `D₄` solver against the oracle on planted and uniform targets.
fn control(args: &[String]) {
    let p: u64 = arg(args, "--p", 269);
    let seed: u64 = arg(args, "--seed", 1);
    let planted: usize = arg(args, "--planted", 3);
    let uniform: usize = arg(args, "--uniform", 3);
    let json_path: String = arg(args, "--json", "control.jsonl".to_string());
    let inst = generate_translation_instance(p, seed, 16).expect("instance");
    let d4 = d4_of(&inst);
    let curve = &d4.curve;
    let mut rng = StdRng::seed_from_u64(seed ^ 0x01F4_1AD0);
    let (yb, _, _) = y_line(&inst, &d4, &mut rng);
    let pts: Vec<Pt4> = yb.iter().map(|(_, q)| *q).collect();
    let oracle = Oracle::build(curve, pts.clone());
    let tp = Instant::now();
    let pre = SymmetrisedD4::precompute(&d4, &mut rng, 8).expect("D4 precompute");
    let pre_s = tp.elapsed().as_secs_f64();
    let mut trng = StdRng::seed_from_u64(seed ^ 0xC0_7701);
    for k in 0..(planted + uniform) {
        let target = if k < planted {
            let mut idx: Vec<usize> = Vec::new();
            while idx.len() < 4 {
                let i = trng.gen_range(0..pts.len());
                if !idx.contains(&i) {
                    idx.push(i);
                }
            }
            let mut acc = Pt4::INF;
            for (n, &i) in idx.iter().enumerate() {
                let q = if n == 2 { curve.neg(&pts[i]) } else { pts[i] };
                acc = curve.add(&acc, &q);
            }
            acc
        } else {
            let (a, b, _) = trial_ab(seed ^ 0xFACE, k as u64, inst.n);
            curve.add(&curve.mul(&inst.g, a), &curve.mul(&inst.q, b))
        };
        let ts = Instant::now();
        let sol = solve_d4(&d4, &pre, &target, 14, &mut rng, false);
        let wall = ts.elapsed().as_secs_f64();
        let orc: BTreeSet<[u64; 4]> = oracle
            .decompose(curve, &target)
            .into_iter()
            .map(|t| {
                let mut ys: [u64; 4] = core::array::from_fn(|n| yb[t[n].0].0);
                ys.sort_unstable();
                ys
            })
            .collect();
        let (solver_set, c4, dim, closed) = match &sol {
            Some(s) => (
                s.decompositions.iter().copied().collect::<BTreeSet<_>>(),
                s.solve.total_muls,
                s.solve.dim,
                true,
            ),
            None => (BTreeSet::new(), 0, 0, false),
        };
        emit(
            &json_path,
            &json!({"p": p, "seed": seed, "kind": if k < planted {"planted"} else {"uniform"},
                    "closed": closed, "dim": dim, "c4_fp_mults": c4, "wall_s": wall,
                    "precompute_s": pre_s,
                    "solver_decompositions": solver_set.len(), "oracle_decompositions": orc.len(),
                    "equal": closed && solver_set == orc}),
        );
    }
}

/// `#E mod 2^k` census of the construction: does `f'(x₀) = c²` force
/// `4 | #E`?
fn census(args: &[String]) {
    let p: u64 = arg(args, "--p", 101);
    let curves: u64 = arg(args, "--curves", 50);
    let json_path: String = arg(args, "--json", "census.jsonl".to_string());
    let mut hist: BTreeMap<u64, u64> = BTreeMap::new();
    for seed in 0..curves {
        if let Ok(inst) = generate_translation_instance(p, 10_000 + seed, 1 << 20) {
            *hist.entry(inst.cofactor).or_default() += 1;
        }
    }
    emit(
        &json_path,
        &json!({"p": p, "curves": curves, "cofactor_hist": hist}),
    );
}

// ── Analysis: fits and pooled constants from the recorded rows (Rust, so
// Python only plots; AGENTS.md "Implementation language").

fn ols(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
    let b = sxy / sxx;
    (b, my - b * mx)
}

/// Slope and intercept of log2(y) on x with a 95% bootstrap interval of
/// the slope over rows.
fn fit(xs: &[f64], ys: &[f64], rng: &mut StdRng) -> Value {
    let ly: Vec<f64> = ys.iter().map(|y| y.log2()).collect();
    let (b, a) = ols(xs, &ly);
    let n = xs.len();
    let mut bs: Vec<f64> = Vec::new();
    for _ in 0..2000 {
        let idx: Vec<usize> = (0..n).map(|_| rng.gen_range(0..n)).collect();
        let x: Vec<f64> = idx.iter().map(|&i| xs[i]).collect();
        if x.iter().all(|v| (v - x[0]).abs() < 1e-9) {
            continue;
        }
        let y: Vec<f64> = idx.iter().map(|&i| ly[i]).collect();
        bs.push(ols(&x, &y).0);
    }
    bs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |f: f64| bs[((bs.len() as f64 - 1.0) * f).round() as usize];
    json!({"slope": b, "intercept_log2": a, "ci95": [q(0.025), q(0.975)], "n": n})
}

fn mean_sd(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let sd = (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0).max(1.0)).sqrt();
    (m, sd)
}

fn analyze(args: &[String]) {
    let input: String = arg(args, "--in", "raw/ladder.jsonl".to_string());
    let out_path: String = arg(args, "--out", "raw/summary.json".to_string());
    let text = std::fs::read_to_string(&input).expect("ladder jsonl");
    let rows: Vec<Value> = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let mut rng = StdRng::seed_from_u64(20261014);
    let mut table: Vec<Value> = Vec::new();
    let mut null_ok = true;
    let mut all_verified = true;
    // arm → (log2E, charged cost, rel gae, la gae, trials, ratio vs rho-neg, la ratio)
    let mut per: BTreeMap<String, Vec<[f64; 7]>> = BTreeMap::new();
    let mut rho_cost: Vec<(f64, f64, f64)> = Vec::new();
    for row in &rows {
        let l2e = row["log2_E"].as_f64().unwrap();
        let r = row["r"].as_f64().unwrap();
        let sr = r.sqrt();
        let rho = row["rho"]["rho"]["mean_S"].as_f64().unwrap();
        let rhoneg = row["rho"]["rho-neg"]["mean_S"].as_f64().unwrap();
        for k in ["rho", "rho-neg"] {
            all_verified &= row["rho"][k]["runs"]
                .as_array()
                .unwrap()
                .iter()
                .all(|w| w["correct"].as_bool() == Some(true));
        }
        rho_cost.push((l2e, rho * sr, rhoneg * sr));
        let arms = row["arms"].as_array().unwrap();
        let get = |n: &str| arms.iter().find(|a| a["arm"] == n).cloned();
        if let (Some(a), Some(b)) = (get("plain"), get("dlp-null")) {
            null_ok &= a["trials"] == b["trials"]
                && a["la_modr_mults"] == b["la_modr_mults"]
                && a["full_relations"] == b["full_relations"];
        }
        let mut cells = serde_json::Map::new();
        for a in arms {
            let name = a["arm"].as_str().unwrap().to_string();
            all_verified &= a["correct"].as_bool() == Some(true);
            let g = &a["gae"];
            let rel = g["trial_scalar"].as_f64().unwrap()
                + g["solver_charged"].as_f64().unwrap()
                + g["verify"].as_f64().unwrap();
            let la = g["la"].as_f64().unwrap();
            let sc = a["S_charged"].as_f64().unwrap();
            let ratio = sc / rhoneg;
            per.entry(name.clone()).or_default().push([
                l2e,
                sc * sr,
                rel,
                la,
                a["trials"].as_f64().unwrap(),
                ratio,
                la / sr / rhoneg,
            ]);
            cells.insert(name, json!({"S_charged": sc, "S_executed": a["S_executed"], "ratio_vs_rho_neg": ratio,
                "ratio_vs_rho": sc / rho, "la_ratio_vs_rho_neg": la / sr / rhoneg, "trials": a["trials"],
                "unknowns_core": a["unknowns_core"], "small_columns": a["small_columns"], "lp_hist": a["lp_hist"],
                "rel_gae": rel, "la_gae": la}));
        }
        table.push(json!({"p": row["p"], "seed": row["seed"], "log2_E": l2e, "log2_r": row["log2_r"], "cofactor": row["cofactor"],
            "rho_S": rho, "rho_neg_S": rhoneg, "arms": cells}));
    }
    let mut fits = serde_json::Map::new();
    for (name, v) in &per {
        let x: Vec<f64> = v.iter().map(|c| c[0]).collect();
        let col = |k: usize| -> Vec<f64> { v.iter().map(|c| c[k]).collect() };
        let la_r = col(6);
        let (m, sd) = mean_sd(&la_r);
        fits.insert(name.clone(), json!({
            "total_charged_cost": fit(&x, &col(1), &mut rng),
            "relation_phase": fit(&x, &col(2), &mut rng),
            "linear_algebra": fit(&x, &col(3), &mut rng),
            "trials": fit(&x, &col(4), &mut rng),
            "ratio_vs_rho_neg": fit(&x, &col(5), &mut rng),
            "la_ratio_vs_rho_neg": fit(&x, &la_r, &mut rng),
            "la_ratio_vs_rho_neg_pooled": {"mean": m, "sd": sd, "se": sd / (la_r.len() as f64).sqrt()},
        }));
    }
    let x: Vec<f64> = rho_cost.iter().map(|c| c.0).collect();
    fits.insert(
        "rho".into(),
        json!({"cost": fit(&x, &rho_cost.iter().map(|c| c.1).collect::<Vec<_>>(), &mut rng)}),
    );
    fits.insert(
        "rho-neg".into(),
        json!({"cost": fit(&x, &rho_cost.iter().map(|c| c.2).collect::<Vec<_>>(), &mut rng)}),
    );
    // Extrapolated handover (labelled as such).  Rho is taken at its
    // large-size S (rows with log2 r ≥ 34, where the walk's set-up is
    // under 10 %), times √r = √(#E / h); each IC phase is its own fitted
    // power law in #E.  First log2 #E in [41, 600] with IC ≤ rho.
    let big: Vec<f64> = table
        .iter()
        .filter(|t| t["log2_r"].as_f64().unwrap() >= 34.0)
        .map(|t| t["rho_neg_S"].as_f64().unwrap())
        .collect();
    let (rho_inf, rho_inf_sd) = mean_sd(&big);
    let h = rows
        .first()
        .map(|r| r["cofactor"].as_f64().unwrap())
        .unwrap_or(4.0);
    let rho_c = |e: f64| rho_inf * 2f64.powf((e - h.log2()) / 2.0);
    let term = |name: &str, k: &str, e: f64| -> f64 {
        let f = &fits[name][k];
        2f64.powf(f["intercept_log2"].as_f64().unwrap() + f["slope"].as_f64().unwrap() * e)
    };
    // Plain and D₄: the linear algebra is N² ∝ √#E, so it enters as its
    // large-size ratio to rho (r_LA); the DLP arm's has its own exponent.
    let la_large = |name: &str| -> f64 {
        let v: Vec<f64> = table
            .iter()
            .filter(|t| t["log2_r"].as_f64().unwrap() >= 34.0)
            .filter_map(|t| t["arms"][name]["la_ratio_vs_rho_neg"].as_f64())
            .collect();
        mean_sd(&v).0
    };
    let total = |name: &str, e: f64| {
        let la = if name == "dlp" {
            term(name, "linear_algebra", e)
        } else {
            la_large(name) * rho_c(e)
        };
        term(name, "relation_phase", e) + la
    };
    let first = |pred: &dyn Fn(f64) -> bool| -> Value {
        let mut e = 41.0;
        while e <= 600.0 {
            if pred(e) {
                return json!((e * 10.0).round() / 10.0);
            }
            e += 0.1;
        }
        Value::Null
    };
    let mut handover = serde_json::Map::new();
    for name in ["plain", "d4", "dlp"] {
        if !fits.contains_key(name) {
            continue;
        }
        handover.insert(
            name.into(),
            json!({
                "log2_E_ic_below_rho_neg": first(&|e| total(name, e) <= rho_c(e)),
                "ratio_at_2^128": total(name, 128.0) / rho_c(128.0),
                "ratio_at_2^160": total(name, 160.0) / rho_c(160.0),
                "ratio_at_2^256": total(name, 256.0) / rho_c(256.0),
                "la_share_at_2^160": 1.0 - term(name, "relation_phase", 160.0) / total(name, 160.0),
                "la_ratio_large_sizes": la_large(name),
            }),
        );
    }
    // Sensitivity: the relation phase anchored at the largest size's mean
    // and carried at the theoretical exponent (1/4; 3/8 for DLP) instead
    // of the fitted one.
    let emax = table
        .iter()
        .map(|t| t["log2_E"].as_f64().unwrap())
        .fold(0.0, f64::max);
    let anchored = |name: &str, e: f64| -> f64 {
        let v: Vec<f64> = table
            .iter()
            .filter(|t| (t["log2_E"].as_f64().unwrap() - emax).abs() < 0.5)
            .filter_map(|t| t["arms"][name]["rel_gae"].as_f64())
            .collect();
        let a = if name == "dlp" { 0.375 } else { 0.25 };
        let rel = mean_sd(&v).0 * 2f64.powf(a * (e - emax));
        let la = if name == "dlp" {
            term(name, "linear_algebra", e)
        } else {
            la_large(name) * rho_c(e)
        };
        rel + la
    };
    for name in ["plain", "d4", "dlp"] {
        if let Some(Value::Object(m)) = handover.get_mut(name) {
            m.insert(
                "anchored_theory_exponent_log2_E_below_rho_neg".into(),
                first(&|e| anchored(name, e) <= rho_c(e)),
            );
        }
    }
    if fits.contains_key("plain") && fits.contains_key("dlp") {
        handover.insert(
            "dlp_below_plain".into(),
            json!({
            "log2_E": first(&|e| total("dlp", e) <= total("plain", e))}),
        );
    }
    handover.insert(
        "rho_neg_S_large".into(),
        json!({"mean": rho_inf, "sd": rho_inf_sd, "n": big.len(), "cofactor": h}),
    );
    let summary = json!({"rows": table, "fits": fits, "null_control_identical": null_ok,
        "every_logarithm_verified": all_verified, "extrapolated_handover": handover,
        "note": "fits: OLS of log2(quantity) on log2 #E over curves, 95% bootstrap over curves (2000); handover is an extrapolation of the fitted power laws, not a measurement"});
    std::fs::write(&out_path, serde_json::to_string_pretty(&summary).unwrap()).unwrap();
    println!("{}", serde_json::to_string_pretty(&json!({"fits": summary["fits"], "null": null_ok, "verified": all_verified, "handover": summary["extrapolated_handover"]})).unwrap());
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("ladder") => ladder(&args),
        Some("control") => control(&args),
        Some("census") => census(&args),
        Some("analyze") => analyze(&args),
        _ => eprintln!("usage: ladder|control|census …"),
    }
}
