//! IC improvement round 2, item 03: fall-row augmentation in the production
//! Koblitz decomposition oracle, and a cheap per-system gate.
//!
//! The incumbent is the pipeline's oracle solve: the equations
//! `decomposition_solver_equations` hands the solver (the target-fixed
//! decomposition system in the interleaved chain order), solved by
//! `solve_boolean_system_filtered` with `SolverEngine::default()` and
//! `split_rule_default()`.  Arm F appends the degree-fall rows of one
//! degree-3 Macaulay pass (`matrix_f4_f2_counted_with`) to the generators.
//! The gate decides per system, before any elimination, whether to do so.
//!
//! The barrel helpers (gauge-fixed G2 barrel systems) are copied from
//! `examples/ic_improve_02_equivariant_oracle.rs`.
//!
//! See `research/ic_improve_20261015/03-fallrow-gate/PROTOCOL.md`.
//!
//! ```bash
//! cargo run --release --example ic_improve2_03_fallrow_gate -- --family pipe --n 17 --ell 6 --m 3 --targets 4 --mode diag
//! cargo run --release --example ic_improve2_03_fallrow_gate -- --family barrel --n 11 --targets 4 --mode measure
//! ```
#![allow(clippy::too_many_arguments, clippy::type_complexity, dead_code)]

use crypto_lib::binary_ecc::f2m::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, f4_profile, f4_profile_reset, f4_word_ops_thread, matrix_f4_f2_counted_with,
    solve_boolean_system_filtered, split_rule_default, FieldStructure, RowCriterion,
    SolveOptions, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, decomposition_solver_equations, find_irreducible,
    KoblitzCurve,
};
use crypto_lib::cryptanalysis::inherited_f4::{InheritCost, ReducedBasis};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn tcpu() -> f64 {
    let mut ts = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    unsafe {
        libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut ts);
    }
    ts.tv_sec as f64 + ts.tv_nsec as f64 * 1e-9
}

// ---------------------------------------------------------------------------
// Solver arms
// ---------------------------------------------------------------------------

#[derive(Default, Clone)]
struct Arm {
    roots: Vec<u64>,
    cpu: f64,
    splits: usize,
    reductions: usize,
    eliminated: usize,
    word_ops: u64,
    calls: u64,
    rank_sum: u64,
    cols_sum: u64,
    spec_ops: u64,
    exhausted: bool,
}

fn solve_full(eqs: &[F2BoolPoly], nv: usize) -> Arm {
    let opts = SolveOptions {
        engine: SolverEngine::default(),
        max_solutions: usize::MAX,
        node_budget: usize::MAX,
        split_rule: split_rule_default(),
    };
    f4_profile_reset();
    let w0 = f4_word_ops_thread();
    let c0 = tcpu();
    let (mut roots, st) = solve_boolean_system_filtered(eqs, nv, &opts, |_| false);
    let cpu = tcpu() - c0;
    let p = f4_profile();
    roots.sort_unstable();
    roots.dedup();
    Arm {
        roots,
        cpu,
        splits: st.splits,
        reductions: st.reductions,
        eliminated: st.eliminated,
        word_ops: f4_word_ops_thread() - w0,
        calls: p.calls,
        rank_sum: p.rows,
        cols_sum: p.cols,
        spec_ops: p.specialise_word_ops,
        exhausted: st.exhausted,
    }
}

fn poly_degree(p: &F2BoolPoly) -> u32 {
    p.terms.iter().map(|t| t.mask.count_ones()).max().unwrap_or(0)
}

fn support(p: &F2BoolPoly) -> u64 {
    p.terms.iter().fold(0u64, |a, t| a | t.mask)
}

/// Degree-fall rows of the degree-3 Macaulay pass (the reduced rows whose
/// every term has degree < 3), with the pass's CPU seconds and word XORs.
fn fall_rows(eqs: &[F2BoolPoly], nv: usize) -> (Vec<F2BoolPoly>, f64, u64) {
    let c0 = tcpu();
    let (rows, ops) =
        matrix_f4_f2_counted_with(eqs, nv, 3, RowCriterion::None).expect("matrix within caps");
    let falls: Vec<F2BoolPoly> = rows.into_iter().filter(|p| poly_degree(p) < 3).collect();
    (falls, tcpu() - c0, ops)
}

/// The incumbent root's own linear tail: the support-local degree-3 basis
/// the solver builds at its root, read with `linear_tail`, keeping the
/// linear rows the solver discards (two or more variables).  CPU seconds of
/// the root build and the read are returned.
fn root_linear_rows(eqs: &[F2BoolPoly], nv: usize) -> (Vec<F2BoolPoly>, f64) {
    let c0 = tcpu();
    let mut cost = InheritCost::default();
    let rows = match ReducedBasis::from_system_with(eqs, nv, 3, true) {
        Some((mut basis, _)) => basis
            .linear_tail(&mut cost)
            .into_iter()
            .filter(|p| p.terms.iter().filter(|t| t.mask != 0).count() >= 2)
            .collect(),
        None => Vec::new(),
    };
    (rows, tcpu() - c0)
}

fn binom(n: u64, k: u64) -> u64 {
    if k > n {
        return 0;
    }
    let mut r = 1u64;
    for i in 0..k {
        r = r * (n - i) / (i + 1);
    }
    r
}

/// The gate's feature, computed by counting alone (no elimination): the
/// expansion ratio `E` of the incumbent's degree-3 support-local Macaulay
/// matrix, rows per generator.  A generator of degree `d` over `s`
/// variables contributes `Σ_{k ≤ 3−d} C(s, k)` rows.  `E ≈ 1` means the
/// incumbent's root basis is just its generators (a cubic system); a large
/// `E` means it already carries the products `x·f`.
fn expansion_ratio(eqs: &[F2BoolPoly]) -> f64 {
    let mut rows = 0u64;
    let mut gens = 0u64;
    for p in eqs {
        if p.is_zero() {
            continue;
        }
        let d = poly_degree(p).min(3);
        let s = support(p).count_ones() as u64;
        gens += 1;
        rows += (0..=(3 - d) as u64).map(|k| binom(s, k)).sum::<u64>();
    }
    rows as f64 / gens.max(1) as f64
}

fn with_falls(eqs: &[F2BoolPoly], falls: &[F2BoolPoly]) -> Vec<F2BoolPoly> {
    let mut v = eqs.to_vec();
    v.extend(falls.iter().cloned());
    v
}

/// Deferred variant: split by hand on the `k` highest occurring variables
/// (the solver's own first `HighestFree` choices), and in each of the `2^k`
/// leaves append the leaf's own fall rows before solving it.  The leaf
/// carries `x_v + c` generators so its roots are roots of the full system.
fn deferred(eqs: &[F2BoolPoly], nv: usize, k: usize) -> (Arm, usize) {
    let occ = eqs.iter().fold(0u64, |a, p| a | support(p));
    let vars: Vec<u32> = (0..nv as u32).rev().filter(|v| occ >> v & 1 == 1).take(k).collect();
    let mut total = Arm::default();
    let mut extra = 0usize;
    let c0 = tcpu();
    for bits in 0..(1u64 << vars.len()) {
        let mut sub: Vec<F2BoolPoly> = eqs.to_vec();
        for (i, &v) in vars.iter().enumerate() {
            let val = bits >> i & 1 == 1;
            sub = sub.iter().map(|p| p.substitute(v, val)).filter(|p| !p.is_zero()).collect();
        }
        if sub.iter().any(|p| p.terms.len() == 1 && p.terms[0].mask == 0) {
            continue;
        }
        let (falls, _, _) = fall_rows(&sub, nv);
        extra += falls.len();
        let mut leaf = with_falls(&sub, &falls);
        for (i, &v) in vars.iter().enumerate() {
            let mut monos = vec![F2BoolMono::from_mask(1u64 << v)];
            if bits >> i & 1 == 1 {
                monos.push(F2BoolMono::from_mask(0));
            }
            leaf.push(F2BoolPoly::from_monos(monos, nv));
        }
        let a = solve_full(&leaf, nv);
        total.roots.extend(a.roots);
        total.splits += a.splits;
        total.reductions += a.reductions;
        total.word_ops += a.word_ops;
        total.calls += a.calls;
        total.rank_sum += a.rank_sum;
        total.cols_sum += a.cols_sum;
        total.exhausted |= a.exhausted;
    }
    total.cpu = tcpu() - c0;
    total.roots.sort_unstable();
    total.roots.dedup();
    (total, extra)
}

// ---------------------------------------------------------------------------
// Systems
// ---------------------------------------------------------------------------

struct PipeCell {
    kc: KoblitzCurve,
    fb: crypto_lib::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    st: FieldStructure,
    m: usize,
}

impl PipeCell {
    fn new(a: u8, n: u32, ell: u32, m: usize) -> Self {
        // a = 0 first; a = 1 where a = 0 has no usable subgroup at this n.
        let kc = KoblitzCurve::new(a, n)
            .or_else(|| KoblitzCurve::new(1 - a, n))
            .expect("curve");
        let fb = build_standard_subspace_factor_base(&kc, ell).expect("factor base");
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        PipeCell { kc, fb, st, m }
    }
    /// Target `t`: a uniform multiple of the generator, seeded per (n, ℓ, m, t).
    fn system(&self, t: u64, ell: u32, plain: bool) -> (Vec<F2BoolPoly>, usize, u64) {
        let ord = self.kc.subgroup_order.to_u64_digits()[0];
        let mut rng =
            StdRng::seed_from_u64(20261015 ^ (self.kc.n as u64 * 1_000_000 + ell as u64 * 1000 + self.m as u64 * 100 + t));
        loop {
            let k = rng.gen_range(1..ord);
            let target = self.kc.mul(self.kc.generator(), &BigUint::from(k));
            if plain {
                // Item 02's lead-probe layout: the system in its own
                // (non-interleaved) variable order.
                let x_r = match &target {
                    crypto_lib::binary_ecc::BinaryPoint::Affine { x, .. } => x.clone(),
                    _ => continue,
                };
                if let Some(sys) = build_decomposition_system(
                    &self.fb.subspace_basis,
                    &x_r,
                    &self.kc.curve.b,
                    self.m,
                    &self.st,
                ) {
                    return (sys.equations, sys.n_vars, k);
                }
                continue;
            }
            if let Some((eqs, nv)) = decomposition_solver_equations(
                &self.kc,
                &self.fb,
                &self.st,
                &target,
                self.m,
                SolverEngine::default(),
            ) {
                return (eqs, nv, k);
            }
        }
    }
}

struct BarrelCtx {
    n: u32,
    nn: usize,
    l: u32,
    irr: IrreduciblePoly,
    f: Field,
    normal: Vec<u64>,
    pinv: Vec<u64>,
}

impl BarrelCtx {
    fn new(n: u32, l: u32) -> Self {
        let irr = find_irreducible(n).expect("irreducible");
        let f = Field::new(n, &irr);
        let mut rng = StdRng::seed_from_u64(0x62 + n as u64);
        let (normal, pinv) = normal_basis(n, &irr, &mut rng);
        BarrelCtx {
            n,
            nn: n as usize,
            l,
            irr,
            f,
            normal,
            pinv,
        }
    }
    fn target(&self, t: u64) -> u64 {
        let mut rng = StdRng::seed_from_u64(20261015 ^ (self.n as u64 * 1000 + t) ^ 0xba44e1);
        rng.gen_range(1..(1u64 << self.n))
    }
    fn system(&self, r: u64) -> (Vec<F2BoolPoly>, usize) {
        let rel = elem_from_bits(r, self.n);
        let (lay, eqs, _) =
            build_symmetric(self.n, self.l, &self.irr, &self.f, 1, &rel, &self.normal, &self.pinv);
        let gnv = lay.n_vars - self.nn;
        assert!(gnv <= 64);
        let geqs = gauge_system(&eqs, self.nn);
        (geqs.iter().map(|p| to_u64(p, gnv)).collect(), gnv)
    }
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

fn arm_json(name: &str, a: &Arm) -> String {
    let per_node = if a.calls > 0 { a.rank_sum as f64 / a.calls as f64 } else { 0.0 };
    format!(
        r#""{name}_cpu":{:.6},"{name}_splits":{},"{name}_reductions":{},"{name}_elim":{},"{name}_word_ops":{},"{name}_calls":{},"{name}_rank_per_call":{:.1},"{name}_cols_per_call":{:.1},"{name}_spec_ops":{},"{name}_roots":{}"#,
        a.cpu,
        a.splits,
        a.reductions,
        a.eliminated,
        a.word_ops,
        a.calls,
        per_node,
        if a.calls > 0 { a.cols_sum as f64 / a.calls as f64 } else { 0.0 },
        a.spec_ops,
        a.roots.len()
    )
}

/// One target, every arm interleaved in this process, in a fixed order.
///
/// * `diag`: I, F, LR, FLR, FL (linear falls only), FQ (degree-2 falls
///   only), D2/D4 (deferred), I2 (A/A), plus the fall-row anatomy.
/// * `train`: I, F, LR, I2 — what the gate is fitted on.
/// * `measure`: I, G (the gated candidate), F, LR, I2.
///
/// The gate (`gens_per_var ≥ theta`, counting only) sends a system to F;
/// otherwise G takes the `nonfire` path (`lr` or `i`).
fn run_target(
    family: &str,
    label: &str,
    t: u64,
    eqs: &[F2BoolPoly],
    nv: usize,
    mode: &str,
    theta: f64,
    nonfire: &str,
    truth: Option<&Vec<u64>>,
) {
    // The gate's feature, by counting alone; its CPU is charged to G.
    let c0 = tcpu();
    let ngens = eqs.iter().filter(|p| !p.is_zero()).count();
    let occ_vars = eqs.iter().fold(0u64, |a, p| a | support(p)).count_ones();
    let gpv = ngens as f64 / occ_vars.max(1) as f64;
    let fires = gpv >= theta;
    let gate_cpu = tcpu() - c0;
    let e = expansion_ratio(eqs);
    let sys_deg = eqs.iter().map(poly_degree).max().unwrap_or(0);

    let i = solve_full(eqs, nv);
    let mut ok = !i.exhausted;
    let theta_s = if theta.is_finite() { format!("{theta}") } else { "null".to_string() };
    let mut out = format!(
        r#"{{"family":"{family}","cell":"{label}","mode":"{mode}","target":{t},"n_vars":{nv},"gens":{ngens},"sys_degree":{sys_deg},"occ_vars":{occ_vars},"gens_per_var":{gpv:.4},"expansion":{e:.3},"gate_cpu":{gate_cpu:.7},"theta":{theta_s},"nonfire":"{nonfire}","fires":{fires},"#
    );
    out += &arm_json("i", &i);

    // G, the gated candidate (measure mode), right after I.
    if mode == "measure" {
        let mut g = if fires {
            let (gf, gp, _) = fall_rows(eqs, nv);
            let mut a = solve_full(&with_falls(eqs, &gf), nv);
            a.cpu += gp;
            a
        } else if nonfire == "lr" {
            let (rows, rp) = root_linear_rows(eqs, nv);
            let mut a = solve_full(&with_falls(eqs, &rows), nv);
            a.cpu += rp;
            a
        } else {
            solve_full(eqs, nv)
        };
        g.cpu += gate_cpu;
        out += ",";
        out += &arm_json("g", &g);
        ok &= !g.exhausted && g.roots == i.roots;
    }

    // F: every degree-fall row of the degree-3 pass appended.
    let (falls, pass_cpu, pass_ops) = fall_rows(eqs, nv);
    let fsys = with_falls(eqs, &falls);
    let f_solve = solve_full(&fsys, nv);
    let mut f = f_solve.clone();
    f.cpu += pass_cpu;
    let mut deg_hist = [0usize; 3];
    let mut supp_sum = 0u64;
    for p in &falls {
        deg_hist[poly_degree(p) as usize] += 1;
        supp_sum += support(p).count_ones() as u64;
    }
    out += &format!(
        r#","falls":{},"falls_deg0":{},"falls_deg1":{},"falls_deg2":{},"fall_mean_support":{:.2},"pass_cpu":{pass_cpu:.6},"pass_word_ops":{pass_ops},"f_solve_cpu":{:.6},"#,
        falls.len(),
        deg_hist[0],
        deg_hist[1],
        deg_hist[2],
        supp_sum as f64 / falls.len().max(1) as f64,
        f_solve.cpu,
    );
    out += &arm_json("f", &f);
    ok &= !f.exhausted && f.roots == i.roots;

    // LR: the incumbent root's own multi-variable linear tail appended.
    let (lrows, lr_pass) = root_linear_rows(eqs, nv);
    let mut lr = solve_full(&with_falls(eqs, &lrows), nv);
    lr.cpu += lr_pass;
    out += &format!(r#","lr_rows":{},"lr_pass_cpu":{lr_pass:.6},"#, lrows.len());
    out += &arm_json("lr", &lr);
    ok &= !lr.exhausted && lr.roots == i.roots;

    if mode == "diag" {
        // FLR: the linear tail of the root basis of (generators + falls) —
        // one MutantXL round at the root — appended to the bare generators.
        let c0 = tcpu();
        let (frows, _) = root_linear_rows(&fsys, nv);
        let flr_pass = pass_cpu + (tcpu() - c0);
        let mut flr = solve_full(&with_falls(eqs, &frows), nv);
        flr.cpu += flr_pass;
        out += &format!(r#","flr_rows":{},"flr_pass_cpu":{flr_pass:.6},"#, frows.len());
        out += &arm_json("flr", &flr);
        let lin: Vec<F2BoolPoly> = falls.iter().filter(|p| poly_degree(p) <= 1).cloned().collect();
        let mut fl = solve_full(&with_falls(eqs, &lin), nv);
        fl.cpu += pass_cpu;
        let quad: Vec<F2BoolPoly> = falls.iter().filter(|p| poly_degree(p) == 2).cloned().collect();
        let mut fq = solve_full(&with_falls(eqs, &quad), nv);
        fq.cpu += pass_cpu;
        let (d2, d2x) = deferred(eqs, nv, 2);
        let (d4, d4x) = deferred(eqs, nv, 4);
        out += ",";
        out += &arm_json("fl", &fl);
        out += ",";
        out += &arm_json("fq", &fq);
        out += ",";
        out += &arm_json("d2", &d2);
        out += &format!(r#","d2_extra":{d2x},"#);
        out += &arm_json("d4", &d4);
        out += &format!(r#","d4_extra":{d4x}"#);
        for a in [&flr, &fl, &fq, &d2, &d4] {
            ok &= !a.exhausted && a.roots == i.roots;
        }
    }
    let i2 = solve_full(eqs, nv);
    out += ",";
    out += &arm_json("i2", &i2);
    ok &= i2.roots == i.roots;
    if let Some(tr) = truth {
        ok &= &i.roots == tr;
    }
    out += &format!(r#","gate_ok":{ok}}}"#);
    println!("{out}");
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let flag = |k: &str| argv.iter().position(|x| x == k).map(|i| argv[i + 1].clone());
    let family = flag("--family").unwrap_or_else(|| "pipe".into());
    let mode = flag("--mode").unwrap_or_else(|| "diag".into());
    let n: u32 = flag("--n").map(|v| v.parse().unwrap()).unwrap_or(13);
    let ell: u32 = flag("--ell").map(|v| v.parse().unwrap()).unwrap_or(4);
    let m: usize = flag("--m").map(|v| v.parse().unwrap()).unwrap_or(3);
    let a: u8 = flag("--a").map(|v| v.parse().unwrap()).unwrap_or(0);
    let targets: u64 = flag("--targets").map(|v| v.parse().unwrap()).unwrap_or(4);
    let first: u64 = flag("--first").map(|v| v.parse().unwrap()).unwrap_or(0);
    let theta: f64 = flag("--theta").map(|v| v.parse().unwrap()).unwrap_or(f64::INFINITY);
    let nonfire = flag("--nonfire").unwrap_or_else(|| "i".into());
    match family.as_str() {
        "pipe" => {
            let cell = PipeCell::new(a, n, ell, m);
            let plain = flag("--layout").as_deref() == Some("plain");
            let label = format!(
                "pipe{}-n{n}-l{ell}-m{m}-a{}",
                if plain { "plain" } else { "" },
                cell.kc.curve.a.raw_bits().first().copied().unwrap_or(0)
            );
            for t in first..first + targets {
                let (eqs, nv, _k) = cell.system(t, ell, plain);
                run_target("pipe", &label, t, &eqs, nv, &mode, theta, &nonfire, None);
            }
        }
        "barrel" => {
            let cx = BarrelCtx::new(n, 3);
            let label = format!("barrel-n{n}-l3");
            for t in first..first + targets {
                let r = cx.target(t);
                let (eqs, nv) = cx.system(r);
                let truth = cx.brute_roots(r);
                run_target("barrel", &label, t, &eqs, nv, &mode, theta, &nonfire, Some(&truth));
            }
        }
        _ => panic!("--family pipe|barrel"),
    }
}

// ---------------------------------------------------------------------------
// Barrel helpers, copied from examples/ic_improve_02_equivariant_oracle.rs
// ---------------------------------------------------------------------------

const WORDS: usize = 3;
type Mask = [u64; WORDS];

fn mask_var(v: usize) -> Mask {
    let mut m = [0u64; WORDS];
    m[v / 64] |= 1u64 << (v % 64);
    m
}
fn mask_or(a: &Mask, b: &Mask) -> Mask {
    let mut m = [0u64; WORDS];
    for i in 0..WORDS {
        m[i] = a[i] | b[i];
    }
    m
}
fn mask_deg(a: &Mask) -> u32 {
    a.iter().map(|w| w.count_ones()).sum()
}

/// Boolean polynomial = sorted set of monomials (F_2 coefficients).
type Poly = Vec<Mask>;

fn poly_add(a: &Poly, b: &Poly) -> Poly {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}
fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    let mut all: Vec<Mask> = Vec::with_capacity(a.len() * b.len());
    for x in a {
        for y in b {
            all.push(mask_or(x, y));
        }
    }
    normalize(all)
}
fn normalize(mut all: Vec<Mask>) -> Poly {
    all.sort_unstable();
    let mut out = Vec::with_capacity(all.len());
    let mut i = 0;
    while i < all.len() {
        let mut j = i;
        while j < all.len() && all[j] == all[i] {
            j += 1;
        }
        if (j - i) % 2 == 1 {
            out.push(all[i]);
        }
        i = j;
    }
    out
}
fn poly_mul_mono(p: &Poly, m: &Mask) -> Poly {
    normalize(p.iter().map(|t| mask_or(t, m)).collect())
}

/// Field structure constants for F_2[z]/(f): reduced[i][j] = z^{i+j} mod f, squares[k] = z^{2k} mod f.
struct Field {
    n: usize,
    reduced: Vec<Vec<u64>>,
    squares: Vec<u64>,
}
fn bits_of(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}
fn elem_from_bits(bits: u64, n: u32) -> F2mElement {
    let positions: Vec<u32> = (0..n).filter(|k| (bits >> k) & 1 == 1).collect();
    F2mElement::from_bit_positions(&positions, n)
}
impl Field {
    fn new(n: u32, irr: &IrreduciblePoly) -> Self {
        let mono = |k: u32| F2mElement::from_bit_positions(&[k], n);
        let nn = n as usize;
        let mut reduced = vec![vec![0u64; nn]; nn];
        for i in 0..nn {
            for j in 0..nn {
                reduced[i][j] = bits_of(&mono(i as u32).mul(&mono(j as u32), irr));
            }
        }
        let squares = (0..nn).map(|k| reduced[k][k]).collect();
        Field {
            n: nn,
            reduced,
            squares,
        }
    }
}

/// Symbolic field element: n coordinate polynomials in the poly basis.
#[derive(Clone)]
struct Sym {
    coords: Vec<Poly>,
}
impl Sym {
    fn zero(n: usize) -> Self {
        Sym {
            coords: vec![Vec::new(); n],
        }
    }
    fn constant(v: u64, n: usize) -> Self {
        Sym {
            coords: (0..n)
                .map(|k| {
                    if (v >> k) & 1 == 1 {
                        vec![[0u64; WORDS]]
                    } else {
                        Vec::new()
                    }
                })
                .collect(),
        }
    }
    fn from_subspace_vars(basis: &[u64], offset: usize, n: usize) -> Self {
        let mut coords: Vec<Vec<Mask>> = vec![Vec::new(); n];
        for (t, &b) in basis.iter().enumerate() {
            for k in 0..n {
                if (b >> k) & 1 == 1 {
                    coords[k].push(mask_var(offset + t));
                }
            }
        }
        Sym {
            coords: coords.into_iter().map(normalize).collect(),
        }
    }
    fn add(&self, o: &Self) -> Self {
        Sym {
            coords: self
                .coords
                .iter()
                .zip(o.coords.iter())
                .map(|(a, b)| poly_add(a, b))
                .collect(),
        }
    }
    fn mul(&self, o: &Self, f: &Field) -> Self {
        let n = f.n;
        let mut acc: Vec<Vec<Mask>> = vec![Vec::new(); n];
        for i in 0..n {
            if self.coords[i].is_empty() {
                continue;
            }
            for j in 0..n {
                if o.coords[j].is_empty() {
                    continue;
                }
                let prod = poly_mul(&self.coords[i], &o.coords[j]);
                let red = f.reduced[i][j];
                for k in 0..n {
                    if (red >> k) & 1 == 1 {
                        acc[k].extend_from_slice(&prod);
                    }
                }
            }
        }
        Sym {
            coords: acc.into_iter().map(normalize).collect(),
        }
    }
    fn square(&self, f: &Field) -> Self {
        let n = f.n;
        let mut acc: Vec<Vec<Mask>> = vec![Vec::new(); n];
        for i in 0..n {
            let sq = f.squares[i];
            for k in 0..n {
                if (sq >> k) & 1 == 1 {
                    acc[k].extend_from_slice(&self.coords[i]);
                }
            }
        }
        Sym {
            coords: acc.into_iter().map(normalize).collect(),
        }
    }
    /// Multiply every coordinate by a Boolean variable.
    fn scale_var(&self, v: usize) -> Self {
        let m = mask_var(v);
        Sym {
            coords: self.coords.iter().map(|p| poly_mul_mono(p, &m)).collect(),
        }
    }
}

/// Binary S3 for y² + xy = x³ + ax² + b: (x1 + x2)² x3² + x1 x2 x3 + (x1 x2)² + b.
fn sym_s3(x1: &Sym, x2: &Sym, x3: &Sym, b: u64, f: &Field) -> Vec<Poly> {
    let s = x1.add(x2).square(f);
    let t1 = s.mul(&x3.square(f), f);
    let p12 = x1.mul(x2, f);
    let t2 = p12.mul(x3, f);
    let t3 = p12.square(f);
    let total = t1.add(&t2).add(&t3).add(&Sym::constant(b, f.n));
    total.coords
}

fn rank_u64(vectors: &[u64]) -> usize {
    let mut rows = vectors.to_vec();
    let mut rank = 0;
    for bit in (0..64).rev() {
        if let Some(p) = (rank..rows.len()).find(|&r| (rows[r] >> bit) & 1 == 1) {
            rows.swap(rank, p);
            let pv = rows[rank];
            for r in 0..rows.len() {
                if r != rank && (rows[r] >> bit) & 1 == 1 {
                    rows[r] ^= pv;
                }
            }
            rank += 1;
        }
    }
    rank
}

/// Normal basis (conjugates of a normal element) and P^{-1} rows.
fn normal_basis(n: u32, irr: &IrreduciblePoly, rng: &mut StdRng) -> (Vec<u64>, Vec<u64>) {
    loop {
        let alpha = elem_from_bits(rng.gen_range(1..(1u64 << n)), n);
        let conj: Vec<u64> = (0..n)
            .map(|j| bits_of(&alpha.square_k_times(j, irr)))
            .collect();
        if rank_u64(&conj) as u32 != n {
            continue;
        }
        let nn = n as usize;
        let mut rows: Vec<(u64, u64)> = (0..nn)
            .map(|i| {
                let mut p_row = 0u64;
                for (j, &c) in conj.iter().enumerate() {
                    if (c >> i) & 1 == 1 {
                        p_row |= 1 << j;
                    }
                }
                (p_row, 1u64 << i)
            })
            .collect();
        for col in 0..nn {
            let piv = (col..nn).find(|&r| (rows[r].0 >> col) & 1 == 1).unwrap();
            rows.swap(col, piv);
            let pv = rows[col];
            for r in 0..nn {
                if r != col && (rows[r].0 >> col) & 1 == 1 {
                    rows[r].0 ^= pv.0;
                    rows[r].1 ^= pv.1;
                }
            }
        }
        return (conj, rows.iter().map(|r| r.1).collect());
    }
}
fn to_normal(coords: &[Poly], pinv: &[u64]) -> Vec<Poly> {
    pinv.iter()
        .map(|row| {
            let mut acc: Poly = Vec::new();
            for (i, c) in coords.iter().enumerate() {
                if (row >> i) & 1 == 1 {
                    acc = poly_add(&acc, c);
                }
            }
            acc
        })
        .collect()
}

struct Layout {
    n: usize,
    s0: usize,
    s1: usize,
    s2: usize,
    c1: usize,
    c2: usize,
    u1: usize,
    u2: usize,
    n_vars: usize,
}
impl Layout {
    fn new(n: usize, l: usize) -> Self {
        let (s0, s1, s2) = (0, n, 2 * n);
        let (c1, c2) = (3 * n, 3 * n + l);
        let (u1, u2) = (3 * n + 2 * l, 4 * n + 2 * l);
        Layout {
            n,
            s0,
            s1,
            s2,
            c1,
            c2,
            u1,
            u2,
            n_vars: 5 * n + 2 * l,
        }
    }
    fn shift_mask(&self, m: &Mask) -> Mask {
        let mut out = *m;
        for base in [self.s0, self.s1, self.s2, self.u1, self.u2] {
            // rotate the n-bit block starting at `base` by one position upward
            let mut bits = Vec::with_capacity(self.n);
            for k in 0..self.n {
                let v = base + k;
                bits.push((out[v / 64] >> (v % 64)) & 1);
            }
            for k in 0..self.n {
                let v = base + ((k + 1) % self.n);
                if bits[k] == 1 {
                    out[v / 64] |= 1 << (v % 64);
                } else {
                    out[v / 64] &= !(1 << (v % 64));
                }
            }
        }
        out
    }
}

#[derive(Clone, Copy)]
enum Fam {
    S3(usize),
    U(u8, usize),
    OneLin(u8),
    OneQuad(u8, usize, usize),
}

fn build_symmetric(
    n: u32,
    l: u32,
    irr: &IrreduciblePoly,
    f: &Field,
    b: u64,
    r: &F2mElement,
    normal: &[u64],
    pinv: &[u64],
) -> (Layout, Vec<Poly>, Vec<Fam>) {
    let nn = n as usize;
    let lay = Layout::new(nn, l as usize);
    let window: Vec<F2mElement> = (0..l)
        .map(|k| F2mElement::from_bit_positions(&[k], n))
        .collect();
    let u1 = Sym::from_subspace_vars(normal, lay.u1, nn);
    let u2 = Sym::from_subspace_vars(normal, lay.u2, nn);
    let barrel = |s_off: usize, c_off: usize| -> Sym {
        let mut acc = Sym::zero(nn);
        for k in 0..nn {
            let shifted: Vec<u64> = window
                .iter()
                .map(|w| bits_of(&w.square_k_times(k as u32, irr)))
                .collect();
            acc = acc.add(&Sym::from_subspace_vars(&shifted, c_off, nn).scale_var(s_off + k));
        }
        acc
    };
    let mut xr = Sym::zero(nn);
    for k in 0..nn {
        xr = xr.add(
            &Sym::constant(bits_of(&r.square_k_times(k as u32, irr)), nn).scale_var(lay.s0 + k),
        );
    }
    let mut eqs = Vec::new();
    let mut fam = Vec::new();
    for (j, e) in to_normal(&sym_s3(&u1, &u2, &xr, b, f), pinv)
        .into_iter()
        .enumerate()
    {
        eqs.push(e);
        fam.push(Fam::S3(j));
    }
    for (i, (u, s_off, c_off)) in [(&u1, lay.s1, lay.c1), (&u2, lay.s2, lay.c2)]
        .into_iter()
        .enumerate()
    {
        let diff = u.add(&barrel(s_off, c_off));
        for (j, e) in to_normal(&diff.coords, pinv).into_iter().enumerate() {
            eqs.push(e);
            fam.push(Fam::U(i as u8, j));
        }
    }
    for (blk, off) in [(0u8, lay.s0), (1u8, lay.s1), (2u8, lay.s2)] {
        let mut monos: Vec<Mask> = (0..nn).map(|k| mask_var(off + k)).collect();
        monos.push([0u64; WORDS]);
        eqs.push(normalize(monos));
        fam.push(Fam::OneLin(blk));
        for k in 0..nn {
            for k2 in (k + 1)..nn {
                eqs.push(vec![mask_or(&mask_var(off + k), &mask_var(off + k2))]);
                fam.push(Fam::OneQuad(blk, k, k2));
            }
        }
    }
    (lay, eqs, fam)
}
fn spec_mask(m: &Mask, n: usize) -> Option<Mask> {
    let low = m[0] & ((1u64 << n) - 1);
    if low & !1 != 0 {
        return None;
    }
    let mut out = [0u64; WORDS];
    for i in 0..WORDS {
        out[i] = m[i] >> n;
        if i + 1 < WORDS {
            out[i] |= m[i + 1] << (64 - n);
        }
    }
    Some(out)
}

fn spec_poly(p: &Poly, n: usize) -> Poly {
    normalize(p.iter().filter_map(|m| spec_mask(m, n)).collect())
}

/// Degree-`degree` Macaulay rows `t·f` (as monomial lists), as `build_macaulay`
/// enumerates them: every multiplier of degree `≤ degree − deg f`.
fn gauge_system(eqs: &[Poly], n: usize) -> Vec<Poly> {
    eqs.iter()
        .map(|p| spec_poly(p, n))
        .filter(|p| !p.is_empty())
        .collect()
}

fn to_u64(p: &Poly, nv: usize) -> F2BoolPoly {
    let monos = p
        .iter()
        .map(|m| {
            assert!(m[1] == 0 && m[2] == 0, "gauge system exceeds 64 unknowns");
            F2BoolMono::from_mask(m[0])
        })
        .collect();
    F2BoolPoly::from_monos(monos, nv)
}
impl BarrelCtx {
    fn fmul(&self, a: u64, b: u64) -> u64 {
        let mut acc = 0u64;
        for i in 0..self.nn {
            if (a >> i) & 1 == 0 {
                continue;
            }
            for j in 0..self.nn {
                if (b >> j) & 1 == 1 {
                    acc ^= self.f.reduced[i][j];
                }
            }
        }
        acc
    }
    fn normal_coords(&self, x: u64) -> u64 {
        let mut u = 0u64;
        for (j, row) in self.pinv.iter().enumerate() {
            if (row & x).count_ones() & 1 == 1 {
                u |= 1 << j;
            }
        }
        u
    }
    /// Every root of the gauge-fixed system, by enumerating
    /// `(k₁, c₁, k₂, c₂)` and testing `S₃(x₁, x₂, x_R) = 0` directly.
    fn brute_roots(&self, r: u64) -> Vec<u64> {
        let (nn, l) = (self.nn, self.l as usize);
        let xs: Vec<Vec<u64>> = (0..nn)
            .map(|k| {
                (0..(1u64 << l))
                    .map(|c| bits_of(&elem_from_bits(c, self.n).square_k_times(k as u32, &self.irr)))
                    .collect()
            })
            .collect();
        let r2 = self.fmul(r, r);
        let mut roots = Vec::new();
        for k1 in 0..nn {
            for c1 in 0..(1u64 << l) {
                let x1 = xs[k1][c1 as usize];
                for k2 in 0..nn {
                    for c2 in 0..(1u64 << l) {
                        let x2 = xs[k2][c2 as usize];
                        let s = x1 ^ x2;
                        let p12 = self.fmul(x1, x2);
                        let v = self.fmul(self.fmul(s, s), r2)
                            ^ self.fmul(p12, r)
                            ^ self.fmul(p12, p12)
                            ^ 1;
                        if v == 0 {
                            let mut mask = 0u64;
                            mask |= 1 << k1;
                            mask |= 1 << (nn + k2);
                            mask |= c1 << (2 * nn);
                            mask |= c2 << (2 * nn + l);
                            mask |= self.normal_coords(x1) << (2 * nn + 2 * l);
                            mask |= self.normal_coords(x2) << (3 * nn + 2 * l);
                            roots.push(mask);
                        }
                    }
                }
            }
        }
        roots.sort_unstable();
        roots
    }
}
