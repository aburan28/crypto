//! # `E(F_{p⁴})` in ecbench: the curve, its counted group, and Gaudry's
//! `k = 4` index calculus
//!
//! Item 1 of `research/ic_next10_20261014` (`01-fp4-ladder/PROTOCOL.md`).
//!
//! The curve is `quartic_folds::generate_translation_instance(p, seed, 16)`.
//! It is the `D₄` module's construction: a rational 2-torsion point
//! `T = (x₀, 0)` with `f'(x₀) = c²`, `c ∈ F_p^*`, coefficients outside
//! `F_{p²}`, and `#E = h·r` with `r` prime.  The planted logarithm of that
//! constructor is dropped: methods see `G`, the target and a seed only.
//! The constructor needs `p ≡ 1 (mod 4)` and `p < 2^15`.
//!
//! A point is written as one 64-bit word, its four `F_p` limbs at 16 bits
//! each.  So the harness's hex targets, verification and replays carry an
//! `E(F_{p⁴})` point unchanged.
//!
//! [`Fp4Group`] implements [`CountedGroup`]: one addition or doubling is
//! one group operation (97 `F_p` multiplications, §11.17), so every
//! generic method — `rho.plain`, `rho.negation`, BSGS, the kangaroo — runs
//! on `E(F_{p⁴})` under its own registered accounting.
//!
//! [`solve_k4`] is `ic.k4_gaudry`, Gaudry's `k = 4` pipeline:
//!
//! - **Base.** Either the `x`-line `{P : x(P) ∈ F_p}` or the `D₄` `Y`-line
//!   `{P : Y(P) ∈ F_p}`, `Y = (x − x₀ − c)/(x − x₀ + c)`.  The `Y`-line
//!   folds `{±P, ±(P + T)}` into one column.
//! - **Relations.** Four-term decompositions of the trials
//!   `R = [a]G + [b]Q`.  Each target contributes one decomposition, drawn
//!   uniformly from its verified ones.
//! - **Large primes (optional).** Up to two points outside the small base
//!   `F'`, `|F'| = ⌈|F|^{lp_exponent}⌉`, cancelled by
//!   `gaudry_cubic::LargePrimeEliminator`.
//! - **Linear algebra.** Filtering to a square core, then Wiedemann
//!   (§11.7).  The answer is accepted only when `[d]G = Q`.
//!
//! **The charge.**  Decompositions are found by an exact enumeration
//! oracle: a sorted table of `x(±P_i ± P_j ± P_k)`, probed with
//! `R − s·P_l`, with hits recovered through the pair table.  It returns
//! the decomposition set the algebraic solvers return; the `D₄` solve
//! matched it on 19 of 19 targets in `gaudry_quartic_d4_20261013`, and on
//! this item's controls.  Each trial is charged the algebraic solve's
//! measured per-trial cost `c4` (parameter, `F_p` multiplications,
//! ÷ 97).  The oracle's own table and probe work is counted as
//! `*_uncharged`, so the record is a lower bound and can promote nothing.
//! Linear algebra is charged 16 `F_p` multiplications per multiplication
//! mod `r` (§11.16).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

use crate::cryptanalysis::gaudry_cubic::{
    square_core, wiedemann_u64, LargePrimeEliminator, SparseRel,
};
use crate::cryptanalysis::gaudry_quartic::{Curve4, Pt4, E4};
use crate::cryptanalysis::gaudry_quartic_d4::{y_base, D4Instance};
use crate::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crate::cryptanalysis::quartic_folds::generate_translation_instance;
use crate::cryptanalysis::residual_walk::mix64;
use crate::hash::sha256::sha256;

/// `F_p` multiplications per `E(F_{p⁴})` affine addition (§11.16–11.17).
pub const C_ADD: f64 = 97.0;
/// `F_p` multiplications charged per multiplication mod `r` (§11.16).
pub const C_MODN: f64 = 16.0;

/// A curve over `F_{p⁴}` with rational 2-torsion and a prime subgroup.
#[derive(Clone, Debug)]
pub struct Fp4Instance {
    pub p: u64,
    pub seed: u64,
    pub a: E4,
    pub b: E4,
    pub x0: E4,
    pub c: u64,
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    pub g: Pt4,
}

fn pack(e: &E4) -> u64 {
    e.0[0] | e.0[1] << 16 | e.0[2] << 32 | e.0[3] << 48
}

fn unpack(v: u64) -> E4 {
    E4(core::array::from_fn(|i| (v >> (16 * i)) & 0xffff))
}

impl Fp4Instance {
    pub fn build(p: u64, seed: u64) -> Result<Fp4Instance, String> {
        if p >= 1 << 15 {
            return Err(format!("fp4_translation takes p < 2^15, not {p}"));
        }
        let fi = generate_translation_instance(p, seed, 16)?;
        Ok(Fp4Instance {
            p,
            seed,
            a: fi.a,
            b: fi.b,
            x0: fi.x0,
            c: fi.c,
            group_order: fi.group_order,
            r: fi.n,
            cofactor: fi.cofactor,
            g: fi.g,
        })
    }

    pub fn curve(&self) -> Curve4 {
        Curve4::with_params(self.p, self.a, self.b)
    }

    pub fn encode(pt: &Pt4) -> Option<[String; 2]> {
        (!pt.inf).then(|| {
            [
                format!("0x{:x}", pack(&pt.x)),
                format!("0x{:x}", pack(&pt.y)),
            ]
        })
    }

    /// The point a hex pair names, checked to lie on the curve.
    pub fn decode(&self, x: u64, y: u64) -> Result<Pt4, String> {
        let (x, y) = (unpack(x), unpack(y));
        if x.0.iter().chain(y.0.iter()).any(|&v| v >= self.p) {
            return Err("a limb is not reduced mod p".into());
        }
        let pt = Pt4 { x, y, inf: false };
        if !self.curve().on_curve(&pt) {
            return Err("the target is not on the curve".into());
        }
        Ok(pt)
    }

    pub fn mul_generator_hex(&self, k: u64) -> Option<[String; 2]> {
        Self::encode(&self.curve().mul(&self.g, k % self.r))
    }

    /// The canonical model JSON and its SHA-256, the last field of the name.
    pub fn model(&self) -> (String, String) {
        let limbs = |e: &E4| format!("[{},{},{},{}]", e.0[0], e.0[1], e.0[2], e.0[3]);
        let f = self.curve().f;
        let json = format!(
            "{{\"field\":\"F_p[t]/(t^4-{})\",\"p\":{},\"a\":{},\"b\":{},\"order\":{}}}",
            f.c,
            self.p,
            limbs(&self.a),
            limbs(&self.b),
            self.group_order
        );
        let h: String = sha256(json.as_bytes())
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        (json, h)
    }

    /// `(icv1, slug)`.  ICV1 has no extension-field form, so the name keeps
    /// ICV1's shape with a `fp4` field tag; it is never registered.
    pub fn names(&self) -> (String, String) {
        let (_, h) = self.model();
        let q = (self.p as u128).pow(4);
        let t = q as i128 + 1 - self.group_order as i128;
        let tt = if t < 0 {
            format!("tm{}", -t)
        } else {
            format!("t{t}")
        };
        (
            format!(
                "ICV1:fp4-p{}:{t}:{}:-:-:-:-:{}",
                self.p,
                self.group_order,
                &h[..12]
            ),
            format!("icv1-fp4p{}-{tt}-{}", self.p, &h[..8]),
        )
    }
}

/// `E(F_{p⁴})` as a [`CountedGroup`]: every addition and doubling counted.
pub struct Fp4Group(pub Curve4);

fn fp(e: &E4) -> u64 {
    mix64(e.0[0] ^ mix64(e.0[1] ^ mix64(e.0[2] ^ mix64(e.0[3] ^ 0x5eed))))
}

impl CountedGroup for Fp4Group {
    type Elt = Pt4;
    fn identity(&self) -> Pt4 {
        Pt4::INF
    }
    fn is_identity(&self, p: &Pt4) -> bool {
        p.inf
    }
    fn add(&self, ops: &mut GroupOps, p: Pt4, q: Pt4) -> Pt4 {
        if p == q {
            ops.doubles += 1;
        } else {
            ops.adds += 1;
        }
        self.0.add(&p, &q)
    }
    fn double(&self, ops: &mut GroupOps, p: Pt4) -> Pt4 {
        ops.doubles += 1;
        self.0.add(&p, &p)
    }
    fn neg(&self, p: Pt4) -> Pt4 {
        self.0.neg(&p)
    }
    fn key(&self, p: &Pt4) -> u64 {
        if p.inf {
            return 0;
        }
        mix64(fp(&p.x) ^ fp(&p.y).rotate_left(23))
    }
}

fn signed(curve: &Curve4, p: &Pt4, s: i64) -> Pt4 {
    if s == 1 {
        *p
    } else {
        curve.neg(p)
    }
}

/// The exact four-term decomposition oracle over one base (the stand-in).
pub struct Oracle {
    pub pts: Vec<Pt4>,
    pairs: HashMap<E4, Vec<(u32, u32, i8, E4)>>,
    triples: Vec<u64>,
    pub setup_adds: u64,
}

impl Oracle {
    pub fn build(curve: &Curve4, pts: Vec<Pt4>) -> Oracle {
        let before = curve.ops();
        let m = pts.len();
        let mut pairs: HashMap<E4, Vec<(u32, u32, i8, E4)>> = HashMap::new();
        let cap = 4 * m * m.saturating_sub(1) * m.saturating_sub(2) / 6;
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
            setup_adds: curve.ops() - before,
        }
    }

    /// Every `R = Σ s_t P_{i_t}` over four distinct base points.
    pub fn decompose(&self, curve: &Curve4, r: &Pt4) -> Vec<Vec<(usize, i64)>> {
        let m = self.pts.len();
        let mut out: BTreeSet<Vec<(usize, i64)>> = BTreeSet::new();
        for l in 0..m {
            for sl in [1i64, -1] {
                let v = curve.sub(r, &signed(curve, &self.pts[l], sl));
                if v.inf || self.triples.binary_search(&fp(&v.x)).is_err() {
                    continue;
                }
                for k in (0..m).filter(|&k| k != l) {
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

/// Which base `ic.k4_gaudry` uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum K4Base {
    XLine,
    YLineD4,
}

#[derive(Clone, Copy, Debug)]
pub struct K4Params {
    pub base: K4Base,
    /// `0` (full relations only) or `2` (double large primes).
    pub large_primes: usize,
    /// `|F'| = ⌈|F|^{lp_exponent}⌉` when `large_primes > 0`.
    pub lp_exponent: f64,
    /// The algebraic solve's cost per trial, `F_p` multiplications.
    pub c4: f64,
    /// Trials before the run counts as exhausted, per column.
    pub trial_cap_per_column: u64,
    /// Item ic_improve_20261015/04: attempt filtering after every new full
    /// relation and solve as soon as the core is square (first-dependency
    /// stop), instead of first at `small + 1` and then every +5 %.
    pub eager_la: bool,
}

/// What [`solve_k4`] did, phase by phase, in GAE.
#[derive(Clone, Debug)]
pub struct K4Outcome {
    pub recovered: Option<u64>,
    pub exhausted: bool,
    /// `(phase, adds, gae)`.
    pub phases: Vec<(&'static str, u64, f64)>,
    pub counters: BTreeMap<String, u64>,
    pub detail: Value,
}

/// `(a, b, pick)` of trial `t`.
pub fn trial_ab(seed: u64, t: u64, r: u64) -> (u64, u64, u64) {
    let mut rng = StdRng::seed_from_u64(mix64(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ t));
    (rng.gen_range(0..r), rng.gen_range(1..r), rng.gen())
}

/// `ic.k4_gaudry` on `inst` for `target`.  Never sees the logarithm.
pub fn solve_k4(inst: &Fp4Instance, target: Pt4, seed: u64, prm: &K4Params) -> K4Outcome {
    let r = inst.r;
    let curve = inst.curve();
    let mut rng = StdRng::seed_from_u64(seed ^ 0x01F4_1AD0);
    let mut la_rng = StdRng::seed_from_u64(seed ^ 0x1A1A);
    // Factor base: native square roots, at 97 F_p multiplications a GAE.
    curve.f.reset_muls();
    #[allow(clippy::type_complexity)]
    let (pts, col, ncols, fb_muls): (Vec<Pt4>, Vec<(usize, i64)>, usize, u64) = match prm.base {
        K4Base::XLine => {
            let pts: Vec<Pt4> = (0..inst.p)
                .filter_map(|x| curve.lift_x(&curve.f.from_fp(x), &mut rng))
                .collect();
            let m = pts.len();
            (pts, (0..m).map(|i| (i, 1)).collect(), m, curve.f.muls())
        }
        K4Base::YLineD4 => {
            let d4 = D4Instance {
                curve: inst.curve(),
                x0: inst.x0,
                c: inst.c,
            };
            let yb = y_base(&d4, &mut rng);
            let t = Pt4 {
                x: inst.x0,
                y: E4::ZERO,
                inf: false,
            };
            let idx: HashMap<u64, usize> =
                yb.iter().enumerate().map(|(i, (y, _))| (*y, i)).collect();
            let mut col = vec![(usize::MAX, 1i64); yb.len()];
            let mut n = 0;
            for (i, (y, pt)) in yb.iter().enumerate() {
                if col[i].0 != usize::MAX {
                    continue;
                }
                col[i] = (n, 1);
                if let Some(&j) = idx.get(&((inst.p - y) % inst.p)) {
                    if j != i {
                        let q = curve.add(pt, &t);
                        col[j] = (n, if q == yb[j].1 { 1 } else { -1 });
                    }
                }
                n += 1;
            }
            let pts = yb.iter().map(|(_, q)| *q).collect();
            (pts, col, n, d4.curve.f.muls() + curve.f.muls())
        }
    };
    let fb_adds = curve.ops();
    curve.reset_ops();
    let m = pts.len();
    // Columns: small base `< small`, `d` at `small`, large primes above.
    let (small, colmap): (usize, Vec<(usize, i64)>) = if prm.large_primes > 0 {
        let s = ((ncols as f64).powf(prm.lp_exponent)).ceil() as usize;
        (
            s,
            col.iter()
                .map(|&(c, k)| (if c < s { c } else { c + 1 }, k))
                .collect(),
        )
    } else {
        (ncols, col)
    };
    let oracle = Oracle::build(&curve, pts);
    let setup_adds = oracle.setup_adds;
    curve.reset_ops();
    let mut elim = (prm.large_primes > 0).then(|| LargePrimeEliminator::new(small + 1, r));
    let mut full: Vec<SparseRel> = Vec::new();
    let mut next_attempt = if prm.eager_la { 1 } else { small + 1 };
    let (mut la_ops, mut merge_ops, mut la_attempts) = (0u64, 0u64, 0u64);
    let (mut scalar_adds, mut probe_adds, mut verify_adds, mut check_adds) = (0u64, 0, 0, 0);
    let (mut decs_total, mut used, mut decomposed) = (0u64, 0u64, 0u64);
    let mut lp_hist = [0u64; 5];
    let cap = prm.trial_cap_per_column * ncols as u64;
    let mut trials = 0u64;
    let mut recovered = None;
    let mut core_unknowns = 0usize;
    while trials < cap && recovered.is_none() {
        let (a, b, pick) = trial_ab(seed, trials, r);
        trials += 1;
        let o0 = curve.ops();
        let rr = curve.add(&curve.mul(&inst.g, a), &curve.mul(&target, b));
        let o1 = curve.ops();
        let decs = oracle.decompose(&curve, &rr);
        let o2 = curve.ops();
        let mut rows: Vec<(SparseRel, usize)> = Vec::new();
        for terms in decs {
            let mut acc = rr;
            for &(i, s) in &terms {
                acc = curve.sub(&acc, &signed(&curve, &oracle.pts[i], s));
            }
            if !acc.inf {
                continue;
            }
            decs_total += 1;
            let mut row: BTreeMap<usize, u64> = BTreeMap::new();
            for &(i, s) in &terms {
                let (c, k) = colmap[i];
                let e = row.entry(c).or_insert(0);
                *e = (*e + if s * k == 1 { 1 } else { r - 1 }) % r;
            }
            let e = row.entry(small).or_insert(0);
            *e = (*e + (r - b) % r) % r;
            let cols: Vec<(usize, u64)> = row.into_iter().filter(|&(_, v)| v != 0).collect();
            let lp = cols.iter().filter(|&&(c, _)| c > small).count();
            if lp <= prm.large_primes {
                rows.push((SparseRel { cols, rhs: a }, lp));
            }
        }
        let o3 = curve.ops();
        scalar_adds += o1 - o0;
        probe_adds += o2 - o1;
        verify_adds += o3 - o2;
        if rows.is_empty() {
            continue;
        }
        decomposed += 1;
        let (rel, lp) = rows.swap_remove((pick % rows.len() as u64) as usize);
        used += 1;
        lp_hist[lp.min(4)] += 1;
        let f = if lp == 0 {
            Some(rel)
        } else {
            elim.as_mut()
                .expect("an eliminator when large primes are on")
                .feed(rel, &mut merge_ops)
        };
        let Some(f) = f else { continue };
        full.push(f);
        if full.len() < next_attempt {
            continue;
        }
        let unknowns = small + 1;
        let Some(core) = square_core(&full, small) else {
            next_attempt = full.len()
                + if prm.eager_la {
                    1
                } else {
                    (unknowns / 20).max(1)
                };
            continue;
        };
        la_attempts += 1;
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
                cols: full[i].cols.iter().map(|&(c, v)| (map[c], v)).collect(),
                rhs: full[i].rhs,
            })
            .collect();
        let x = wiedemann_u64(&sel, sel.len(), r, &mut la_rng, &mut la_ops);
        let o4 = curve.ops();
        let ok = x
            .map(|x| x[map[small]])
            .filter(|&dd| curve.mul(&inst.g, dd) == target);
        check_adds += curve.ops() - o4;
        if ok.is_some() {
            recovered = ok;
            core_unknowns = sel.len();
        } else {
            next_attempt = full.len()
                + if prm.eager_la {
                    1
                } else {
                    (sel.len() / 20).max(1)
                };
        }
    }
    let solver_gae = trials as f64 * prm.c4 / C_ADD;
    let la_gae = (la_ops + merge_ops) as f64 * C_MODN / C_ADD;
    let phases = vec![
        (
            "factor_base",
            fb_adds,
            fb_muls as f64 / C_ADD + fb_adds as f64,
        ),
        (
            "relations",
            scalar_adds + verify_adds,
            (scalar_adds + verify_adds) as f64 + solver_gae,
        ),
        ("linear_algebra", check_adds, la_gae + check_adds as f64),
    ];
    let mut counters = BTreeMap::new();
    counters.insert("base_points".into(), m as u64);
    counters.insert("columns".into(), ncols as u64);
    counters.insert("small_columns".into(), small as u64);
    counters.insert("trials".into(), trials);
    counters.insert("decompositions".into(), decs_total);
    counters.insert("targets_decomposed".into(), decomposed);
    counters.insert("relations_used".into(), used);
    counters.insert("full_relations".into(), full.len() as u64);
    counters.insert("core_unknowns".into(), core_unknowns as u64);
    counters.insert("la_attempts".into(), la_attempts);
    counters.insert("wiedemann_modr_mults".into(), la_ops);
    counters.insert("merge_modr_mults".into(), merge_ops);
    counters.insert("solver_charged_trials".into(), trials);
    counters.insert("oracle_table_adds_uncharged".into(), setup_adds);
    counters.insert("oracle_probe_adds_uncharged".into(), probe_adds);
    for (i, v) in lp_hist.iter().enumerate() {
        counters.insert(format!("relations_with_{i}_large_primes"), *v);
    }
    K4Outcome {
        recovered,
        exhausted: recovered.is_none(),
        phases,
        counters,
        detail: json!({
            "c4_fp_mults": prm.c4,
            "c_add": C_ADD,
            "c_modn": C_MODN,
            "base": format!("{:?}", prm.base),
            "large_primes": prm.large_primes,
            "lp_exponent": prm.lp_exponent,
            "solver_charged_gae": solver_gae,
            "la_gae": la_gae,
        }),
    }
}
