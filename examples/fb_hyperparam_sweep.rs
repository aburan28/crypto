//! Factor-base hyperparameter sweep: relation yield vs. first-fall degree
//! `d_ff` vs. refutation degree `D*` vs. measured decomposition cost, on one
//! grid and one metric.
//!
//! Study: `research/fb_hyperparam_sweep_20261014/README.md`.
//!
//! Three arms, each joined into one JSON table:
//!
//! - **binary** — `F_{2ⁿ}`, 2-decompositions (`S₃`), factor base = an
//!   `F_2`-subspace `V` of dimension `n'` from
//!   [`descent_lowgamma::FactorSubspace`]: `Coordinate` (Gaudry/Diem
//!   `span{1, z, …, z^{n'−1}}`), `Random` (full-rank basis, control) and
//!   `Subfield` (`V = F_{2^{n'}}`, `n' | n`, the Frobenius-stable factor
//!   base). Two curve classes: generic `b` and Koblitz `b = 1`. Per cell:
//!   measured yield (exact root-finding oracle, cross-checked against
//!   brute force over `V × V`), `d_ff` and `D*` on non-decomposable targets
//!   (`pc_degree_harness::refutation_scan`, unchanged), the calibration that
//!   decomposable targets never refute, and the measured per-attempt time of
//!   each oracle (root finding, brute force, Macaulay up to `D*`).
//! - **prime** — the negative control: complete DLP solves on `a = −3`
//!   prime-order curves with the small-`x` factor base
//!   ([`prime_fast::ic_solve`], `m ∈ {2, 3}`) across factor-base sizes,
//!   against parallel rho on the same curve and the same secret.
//!
//! After the grid, [`analyze`] applies the protocol's decision rules
//! (`research/fb_hyperparam_sweep_20261014/PROTOCOL.md`) and the derived
//! cost model to the recorded values and writes them under `analysis` and
//! `derived_cost_model`; the plotting script only lays them out.
//!
//! ```bash
//! cargo run --release --example fb_hyperparam_sweep -- <out_dir> [quick]
//! ```

use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use crypto_lib::binary_ecc::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::descent_lowgamma::{
    descend_on_subspace, is_decomposable_on_subspace, BasisFamily, FactorSubspace,
};
use crypto_lib::cryptanalysis::ffd_harness::{choose_irreducible, random_nonzero_f2m};
use crypto_lib::cryptanalysis::pc_degree_harness::refutation_scan;
use crypto_lib::cryptanalysis::prime_fast::{find_a3_curve, ic_solve, rho_parallel, RhoConfig};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

const SEED: u64 = 0x_FB_5EE9_2026_1014;

// ── F_2-linear helpers (n ≤ 64) ─────────────────────────────────────

fn word(e: &F2mElement) -> u64 {
    e.raw_bits()[0]
}

fn elem(w: u64, n: u32) -> F2mElement {
    F2mElement::from_words(&[w], n)
}

/// Echelon basis keyed by leading bit; each pivot carries a tag recording
/// which inputs were combined, so `reduce` also returns a preimage.
struct Echelon {
    piv: Vec<Option<(u64, u64)>>,
}

impl Echelon {
    fn new(n: u32) -> Self {
        Echelon {
            piv: vec![None; n as usize],
        }
    }

    fn insert(&mut self, mut v: u64, mut tag: u64) {
        while v != 0 {
            let lead = 63 - v.leading_zeros() as usize;
            match self.piv[lead] {
                Some((pv, pt)) => {
                    v ^= pv;
                    tag ^= pt;
                }
                None => {
                    self.piv[lead] = Some((v, tag));
                    return;
                }
            }
        }
    }

    /// Returns `(residual, tag)`: `residual == 0` iff `v` is in the span,
    /// and then `v = Σ_{i ∈ tag} input_i`.
    fn reduce(&self, mut v: u64) -> (u64, u64) {
        let mut tag = 0u64;
        while v != 0 {
            let lead = 63 - v.leading_zeros() as usize;
            match self.piv[lead] {
                Some((pv, pt)) => {
                    v ^= pv;
                    tag ^= pt;
                }
                None => return (v, tag),
            }
        }
        (0, tag)
    }
}

/// Exact `S₃` root-finding decomposition oracle over an arbitrary
/// subspace `V`: for each `x₁ ∈ V`, `S₃(x₁, X, x₃)` is the quadratic
/// `A·X² + C·X + D` with `A = x₃² + x₁²`, `C = x₁x₃`, `D = x₁²x₃² + b`;
/// solve it in `F_{2ⁿ}` (linear algebra on `Y² + Y`) and test both roots
/// for membership in `V`. `|V|` field-equation solves instead of the
/// `|V|²` evaluations of brute force.
struct RootFinder {
    n: u32,
    irr: IrreduciblePoly,
    /// `Y ↦ Y² + Y` as an echelon basis over the monomial images.
    as_map: Echelon,
    /// The factor-base subspace `V`.
    v_span: Echelon,
    v_elems: Vec<F2mElement>,
}

impl RootFinder {
    fn new(n: u32, irr: &IrreduciblePoly, v: &FactorSubspace) -> Self {
        let mut as_map = Echelon::new(n);
        for i in 0..n {
            let zi = F2mElement::from_bit_positions(&[i], n);
            as_map.insert(word(&zi.square(irr).add(&zi)), 1u64 << i);
        }
        let mut v_span = Echelon::new(n);
        for c in &v.cols {
            v_span.insert(word(c), 0);
        }
        let v_elems = (0..1u32 << v.n_sub).map(|m| v.element(m)).collect();
        RootFinder {
            n,
            irr: irr.clone(),
            as_map,
            v_span,
            v_elems,
        }
    }

    fn in_v(&self, x: &F2mElement) -> bool {
        self.v_span.reduce(word(x)).0 == 0
    }

    fn sqrt(&self, x: &F2mElement) -> F2mElement {
        x.square_k_times(self.n - 1, &self.irr)
    }

    /// Roots in `F_{2ⁿ}` of `A·X² + C·X + D` (`None` = every `X`, which
    /// only happens for the zero polynomial).
    fn quadratic_roots(
        &self,
        a: &F2mElement,
        c: &F2mElement,
        d: &F2mElement,
    ) -> Option<Vec<F2mElement>> {
        let irr = &self.irr;
        Some(if a.is_zero() {
            // Linear: C·X = D.
            if c.is_zero() {
                if d.is_zero() {
                    return None;
                }
                vec![]
            } else {
                vec![d.mul(&c.flt_inverse(irr).unwrap(), irr)]
            }
        } else if c.is_zero() {
            // A·X² = D: unique root √(D/A).
            vec![self.sqrt(&d.mul(&a.flt_inverse(irr).unwrap(), irr))]
        } else {
            // X = (C/A)·Y, Y² + Y = A·D/C².
            let cinv = c.flt_inverse(irr).unwrap();
            let rhs = a.mul(d, irr).mul(&cinv.square(irr), irr);
            let (res, tag) = self.as_map.reduce(word(&rhs));
            if res != 0 {
                vec![]
            } else {
                let scale = c.mul(&a.flt_inverse(irr).unwrap(), irr);
                let r0 = scale.mul(&elem(tag, self.n), irr);
                let r1 = r0.add(&scale);
                vec![r0, r1]
            }
        })
    }

    fn decomposes(&self, b: &F2mElement, x3: &F2mElement) -> bool {
        let irr = &self.irr;
        let x3sq = x3.square(irr);
        for x1 in &self.v_elems {
            let x1sq = x1.square(irr);
            let a = x3sq.add(&x1sq);
            let c = x1.mul(x3, irr);
            let d = x1sq.mul(&x3sq, irr).add(b);
            let hit = match self.quadratic_roots(&a, &c, &d) {
                None => true,
                Some(rs) => rs.iter().any(|r| self.in_v(r)),
            };
            if hit {
                return true;
            }
        }
        false
    }

    /// Plant a decomposable target: pick `x₁, x₂ ∈ V` and return a nonzero
    /// root `x₃` of `S₃(x₁, x₂, X₃)`, when one exists in `F_{2ⁿ}`.
    fn plant<R: Rng>(&self, rng: &mut R, b: &F2mElement) -> Option<F2mElement> {
        let irr = &self.irr;
        let k = self.v_elems.len();
        let x1 = &self.v_elems[rng.gen_range(0..k)];
        let x2 = &self.v_elems[rng.gen_range(0..k)];
        let a = x1.add(x2).square(irr);
        let p = x1.mul(x2, irr);
        let d = p.square(irr).add(b);
        self.quadratic_roots(&a, &p, &d)?
            .into_iter()
            .find(|r| !r.is_zero())
    }
}

// ── Binary arm ──────────────────────────────────────────────────────

#[derive(Clone, Copy)]
enum CurveClass {
    Generic,
    Koblitz,
}

impl CurveClass {
    fn name(self) -> &'static str {
        match self {
            CurveClass::Generic => "generic_b",
            CurveClass::Koblitz => "koblitz_b1",
        }
    }
}

fn family_name(f: BasisFamily) -> &'static str {
    match f {
        BasisFamily::Coordinate => "coordinate",
        BasisFamily::Subfield => "subfield",
        BasisFamily::Random => "random",
    }
}

struct BinaryBudget {
    yield_targets: u32,
    crosscheck_targets: u32,
    dstar_nondecomp: u32,
    dstar_decomp: u32,
    d_max: u32,
    curves: u32,
    planted: u32,
}

fn mean(v: &[f64]) -> Option<f64> {
    if v.is_empty() {
        None
    } else {
        Some(v.iter().sum::<f64>() / v.len() as f64)
    }
}

fn sd(v: &[f64]) -> Option<f64> {
    let m = mean(v)?;
    if v.len() < 2 {
        return None;
    }
    Some((v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (v.len() - 1) as f64).sqrt())
}

/// Wilson 95% interval for a binomial proportion.
fn wilson(k: u32, n: u32) -> (f64, f64) {
    if n == 0 {
        return (0.0, 1.0);
    }
    let z = 1.959964f64;
    let nf = n as f64;
    let p = k as f64 / nf;
    let den = 1.0 + z * z / nf;
    let centre = (p + z * z / (2.0 * nf)) / den;
    let half = z * ((p * (1.0 - p) / nf + z * z / (4.0 * nf * nf)).sqrt()) / den;
    ((centre - half).max(0.0), (centre + half).min(1.0))
}

fn binary_cell(
    class: CurveClass,
    family: BasisFamily,
    n: u32,
    n_sub: u32,
    budget: &BinaryBudget,
    seed: u64,
) -> Option<Value> {
    let irr = choose_irreducible(n);
    let v = FactorSubspace::build(family, n, n_sub, &irr, seed)?;
    let rf = RootFinder::new(n, &irr, &v);
    let mut rng = StdRng::seed_from_u64(seed);

    let bs: Vec<F2mElement> = match class {
        CurveClass::Koblitz => vec![F2mElement::one(n)],
        CurveClass::Generic => (0..budget.curves)
            .map(|_| random_nonzero_f2m(&mut rng, n))
            .collect(),
    };

    // 1. Yield by root finding; brute-force cross-check on a prefix.
    let mut decomposed = 0u32;
    let mut disagreements = 0u32;
    let mut crosschecked = 0u32;
    let mut t_root = 0f64;
    let mut t_brute = 0f64;
    let mut brute_n = 0u32;
    let mut pool_nd: Vec<(F2mElement, F2mElement)> = Vec::new();
    let mut pool_d: Vec<(F2mElement, F2mElement)> = Vec::new();
    for t in 0..budget.yield_targets {
        let b = &bs[t as usize % bs.len()];
        let x3 = random_nonzero_f2m(&mut rng, n);
        let s = Instant::now();
        let dec = rf.decomposes(b, &x3);
        t_root += s.elapsed().as_secs_f64();
        if dec {
            decomposed += 1;
        }
        if t < budget.crosscheck_targets {
            let s = Instant::now();
            let brute = is_decomposable_on_subspace(&v, &irr, b, &x3);
            t_brute += s.elapsed().as_secs_f64();
            brute_n += 1;
            crosschecked += 1;
            if brute != dec {
                disagreements += 1;
            }
        }
        if dec && pool_d.len() < budget.dstar_decomp as usize {
            pool_d.push((b.clone(), x3));
        } else if !dec && pool_nd.len() < budget.dstar_nondecomp as usize {
            pool_nd.push((b.clone(), x3));
        }
    }

    // Positive control: planted decomposable targets must all be found.
    let mut planted = 0u32;
    let mut planted_found = 0u32;
    let mut planted_brute_found = 0u32;
    let mut attempts = 0u32;
    while planted < budget.planted && attempts < 16 * budget.planted {
        attempts += 1;
        let b = &bs[attempts as usize % bs.len()];
        if let Some(x3) = rf.plant(&mut rng, b) {
            planted += 1;
            if rf.decomposes(b, &x3) {
                planted_found += 1;
            }
            if planted <= 16 && is_decomposable_on_subspace(&v, &irr, b, &x3) {
                planted_brute_found += 1;
            }
        }
    }

    // 2. d_ff and D* on non-decomposable targets; calibration on
    //    decomposable ones (they must never refute).
    let nvars = 2 * n_sub;
    let mut ff: Vec<f64> = Vec::new();
    let mut ds: Vec<f64> = Vec::new();
    let mut ds_hist = std::collections::BTreeMap::<String, u32>::new();
    let mut censored = 0u32;
    let mut t_mac: Vec<f64> = Vec::new();
    let mut mac_cols_at_dstar: Vec<f64> = Vec::new();
    for (b, x3) in &pool_nd {
        let s = Instant::now();
        let eqs = descend_on_subspace(n, &v, &irr, b, x3);
        let (f, d, per) = refutation_scan(&eqs, nvars, budget.d_max);
        t_mac.push(s.elapsed().as_secs_f64());
        if let Some(f) = f {
            ff.push(f as f64);
        }
        match d {
            Some(d) => {
                ds.push(d as f64);
                *ds_hist.entry(d.to_string()).or_default() += 1;
                if let Some(last) = per.last() {
                    mac_cols_at_dstar.push(last.cols as f64);
                }
            }
            None => {
                censored += 1;
                *ds_hist.entry(format!(">{}", budget.d_max)).or_default() += 1;
            }
        }
    }
    // Censored targets enter the comparisons as d_max + 1 (lower bound).
    let ds_with_censored: Vec<f64> = ds
        .iter()
        .cloned()
        .chain(std::iter::repeat_n(
            (budget.d_max + 1) as f64,
            censored as usize,
        ))
        .collect();
    let mut calib_refuted = 0u32;
    for (b, x3) in &pool_d {
        let eqs = descend_on_subspace(n, &v, &irr, b, x3);
        let (_, d, _) = refutation_scan(&eqs, nvars, budget.d_max);
        if d.is_some() {
            calib_refuted += 1;
        }
    }

    let (lo, hi) = wilson(decomposed, budget.yield_targets);
    let v_size = 1u64 << n_sub;
    Some(json!({
        "arm": "binary",
        "curve_class": class.name(),
        "family": family_name(family),
        "n": n,
        "n_sub": n_sub,
        "irreducible_low_terms": irr.low_terms,
        "curves": bs.len(),
        "fb_x_values": v_size,
        "yield": {
            "targets": budget.yield_targets,
            "decomposed": decomposed,
            "rate": decomposed as f64 / budget.yield_targets as f64,
            "wilson95": [lo, hi],
            "heuristic_V2_over_2n": (v_size * v_size) as f64 / 2f64.powi(n as i32),
        },
        "crosscheck": {"targets": crosschecked, "disagreements": disagreements},
        "positive_control": {
            "planted": planted,
            "found_root_find": planted_found,
            "brute_checked": planted.min(16),
            "found_brute": planted_brute_found,
        },
        "degrees": {
            "d_max": budget.d_max,
            "nondecomp_targets": pool_nd.len(),
            "first_fall_mean": mean(&ff),
            "first_fall_defined": ff.len(),
            "dstar_mean": mean(&ds),
            "dstar_sd": sd(&ds),
            "dstar_min": ds.iter().cloned().fold(None, |a: Option<f64>, x| Some(a.map_or(x, |a| a.min(x)))),
            "dstar_max": ds.iter().cloned().fold(None, |a: Option<f64>, x| Some(a.map_or(x, |a| a.max(x)))),
            "dstar_censored": censored,
            "dstar_hist": ds_hist,
            "dstar_values_censored_as_dmax_plus1": ds_with_censored,
            "first_fall_values": ff,
            "macaulay_cols_at_dstar_mean": mean(&mac_cols_at_dstar),
            "calibration_decomp_targets": pool_d.len(),
            "calibration_decomp_refuted": calib_refuted,
        },
        "seconds_per_attempt": {
            "root_find": t_root / budget.yield_targets as f64,
            "brute_force": if brute_n > 0 { Some(t_brute / brute_n as f64) } else { None },
            "macaulay_to_dstar_mean": mean(&t_mac),
            "macaulay_to_dstar_sd": sd(&t_mac),
            "macaulay_samples": t_mac.len(),
        },
    }))
}

/// Per-operation timings used by the cost model (field mul/inverse in
/// each `F_{2ⁿ}`, measured on this host).
fn field_op_timings(n: u32) -> Value {
    let irr = choose_irreducible(n);
    let mut rng = StdRng::seed_from_u64(SEED ^ n as u64);
    let xs: Vec<F2mElement> = (0..256).map(|_| random_nonzero_f2m(&mut rng, n)).collect();
    let reps = 200_000usize;
    let s = Instant::now();
    let mut acc = xs[0].clone();
    for i in 0..reps {
        acc = acc.mul(&xs[i & 255], &irr);
        if acc.is_zero() {
            acc = xs[1].clone();
        }
    }
    let t_mul = s.elapsed().as_secs_f64() / reps as f64;
    let reps_inv = 50_000usize;
    let s = Instant::now();
    let mut sink = 0u64;
    for i in 0..reps_inv {
        sink ^= word(&xs[i & 255].flt_inverse(&irr).unwrap());
    }
    let t_inv = s.elapsed().as_secs_f64() / reps_inv as f64;
    json!({"n": n, "mul_s": t_mul, "inv_s": t_inv, "sink": sink ^ word(&acc)})
}

// ── Prime arm (negative control) ────────────────────────────────────

fn prime_arm(bits_list: &[u32], fb_sizes: &[usize], reps: u64) -> Vec<Value> {
    let mut out = Vec::new();
    for &bits in bits_list {
        let (curve, g) = find_a3_curve(bits, 3).expect("a=-3 bench curve");
        for rep in 0..reps {
            let mut rng = StdRng::seed_from_u64(SEED ^ ((bits as u64) << 32) ^ rep);
            let secret = rng.gen_range(2..curve.n - 1);
            let q = curve.mul(Some(g), secret).unwrap();

            let mut cfg = RhoConfig::for_bits(bits, true);
            cfg.threads = 1;
            cfg.seed = SEED ^ rep;
            let s = Instant::now();
            let rho = rho_parallel(&curve, g, q, &cfg);
            let rho_s = s.elapsed().as_secs_f64();
            out.push(json!({
                "arm": "prime", "method": "rho_negation_1thread", "bits": bits,
                "p": curve.f.p, "order": curve.n, "rep": rep,
                "solved": rho.as_ref().map(|r| r.log == secret).unwrap_or(false),
                "steps": rho.as_ref().map(|r| r.total_steps),
                "seconds": rho_s,
            }));

            for &m in &[2u32, 3] {
                for &fb in fb_sizes {
                    if m == 3 && fb > 4096 {
                        continue;
                    }
                    let max_targets = 1u64 << (bits + 2).min(36);
                    let s = Instant::now();
                    let res = ic_solve(&curve, g, q, fb, 8, m, max_targets, SEED ^ rep);
                    let secs = s.elapsed().as_secs_f64();
                    let row = match res {
                        Some((x, r)) => json!({
                            "arm": "prime", "method": format!("ic_m{m}"), "bits": bits,
                            "p": curve.f.p, "order": curve.n, "rep": rep, "fb_size": r.factor_base_size,
                            "solved": x == secret, "targets": r.targets, "relations": r.relations,
                            "yield": r.relations as f64 / r.targets as f64,
                            "relations_ms": r.relations_ms, "linear_algebra_ms": r.linear_algebra_ms,
                            "pair_table_ms": r.pair_table_ms, "seconds": secs,
                        }),
                        None => json!({
                            "arm": "prime", "method": format!("ic_m{m}"), "bits": bits,
                            "p": curve.f.p, "order": curve.n, "rep": rep, "fb_size": fb,
                            "solved": false, "failure": "relation budget exhausted or LA failed",
                            "max_targets": max_targets, "seconds": secs,
                        }),
                    };
                    eprintln!("prime bits={bits} rep={rep} m={m} fb={fb}: {secs:.3}s");
                    out.push(row);
                }
            }
        }
    }
    out
}

// ── Analysis (protocol decision rules and the derived cost model) ────

/// Inverse standard-normal CDF (Acklam), for Bonferroni-adjusted z.
#[allow(clippy::excessive_precision)]
fn norm_ppf(p: f64) -> f64 {
    let a = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    let b = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    let c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    let d = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    let pl = 0.02425;
    if p < pl {
        let q = (-2.0 * p.ln()).sqrt();
        (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    } else if p <= 1.0 - pl {
        let q = p - 0.5;
        let r = q * q;
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    } else {
        -norm_ppf(1.0 - p)
    }
}

fn wilson_z(k: f64, n: f64, z: f64) -> (f64, f64) {
    let p = k / n;
    let den = 1.0 + z * z / n;
    let centre = (p + z * z / (2.0 * n)) / den;
    let half = z * ((p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt()) / den;
    ((centre - half).max(0.0), (centre + half).min(1.0))
}

/// Percentile bootstrap of mean(x) − mean(y); returns (diff, lo, hi) at
/// two-sided level `alpha`.
fn bootstrap_diff(
    x: &[f64],
    y: &[f64],
    alpha: f64,
    reps: usize,
    seed: u64,
) -> Option<(f64, f64, f64)> {
    if x.is_empty() || y.is_empty() {
        return None;
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let diff = mean(x)? - mean(y)?;
    let mut ds: Vec<f64> = (0..reps)
        .map(|_| {
            let mx: f64 = (0..x.len())
                .map(|_| x[rng.gen_range(0..x.len())])
                .sum::<f64>()
                / x.len() as f64;
            let my: f64 = (0..y.len())
                .map(|_| y[rng.gen_range(0..y.len())])
                .sum::<f64>()
                / y.len() as f64;
            mx - my
        })
        .collect();
    ds.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let lo = ds[((alpha / 2.0) * reps as f64).floor() as usize];
    let hi = ds[(((1.0 - alpha / 2.0) * reps as f64).ceil() as usize).min(reps - 1)];
    Some((diff, lo, hi))
}

/// Newcombe hybrid-score interval for p1 − p2 at normal quantile `z`.
fn newcombe(k1: f64, n1: f64, k2: f64, n2: f64, z: f64) -> (f64, f64, f64) {
    let (p1, p2) = (k1 / n1, k2 / n2);
    let (l1, u1) = wilson_z(k1, n1, z);
    let (l2, u2) = wilson_z(k2, n2, z);
    let d = p1 - p2;
    let lo = d - ((p1 - l1).powi(2) + (u2 - p2).powi(2)).sqrt();
    let hi = d + ((u1 - p1).powi(2) + (p2 - l2).powi(2)).sqrt();
    (d, lo, hi)
}

fn controls_ok(r: &Value) -> bool {
    let pc = &r["positive_control"];
    r["crosscheck"]["disagreements"].as_u64() == Some(0)
        && pc["planted"].as_u64() == pc["found_root_find"].as_u64()
        && pc["brute_checked"].as_u64() == pc["found_brute"].as_u64()
        && r["degrees"]["calibration_decomp_refuted"].as_u64() == Some(0)
}

fn fvec(v: &Value) -> Vec<f64> {
    v.as_array()
        .map(|a| a.iter().filter_map(Value::as_f64).collect())
        .unwrap_or_default()
}

fn analyze(rows: &mut [Value], prime: &[Value], field_ops: &[Value], la_op_s: f64) -> Value {
    // Derived cost model per cell.
    for r in rows.iter_mut() {
        let n = r["n"].as_u64().unwrap() as u32;
        let n_sub = r["n_sub"].as_u64().unwrap() as f64;
        let v = r["fb_x_values"].as_f64().unwrap();
        let frob = r["curve_class"] == "koblitz_b1" && r["family"] == "subfield";
        let f = if frob { n_sub } else { 1.0 };
        let rel = v / f + 10.0;
        let y = r["yield"]["rate"].as_f64().unwrap();
        let fo = field_ops
            .iter()
            .find(|x| x["n"].as_u64() == Some(n as u64))
            .unwrap();
        let t_add = fo["inv_s"].as_f64().unwrap() + 2.0 * fo["mul_s"].as_f64().unwrap();
        let t_root = r["seconds_per_attempt"]["root_find"].as_f64().unwrap();
        let t_mac = r["seconds_per_attempt"]["macaulay_to_dstar_mean"].as_f64();
        let attempts = if y > 0.0 { Some(rel / y) } else { None };
        let la = 2.0 * rel * rel * la_op_s;
        let dec_root = attempts.map(|a| a * t_root);
        let dec_mac = attempts.and_then(|a| t_mac.map(|t| a * t));
        r["derived_cost_model"] = json!({
            "frobenius_orbit_factor": f,
            "relations_needed": rel,
            "attempts": attempts,
            "decomposition_seconds_root_find": dec_root,
            "decomposition_seconds_macaulay": dec_mac,
            "linear_algebra_seconds": la,
            "total_seconds_root_find": dec_root.map(|d| d + la),
            "total_seconds_macaulay": dec_mac.map(|d| d + la),
            "rho_reference_seconds": (std::f64::consts::PI * 2f64.powi(n as i32) / 4.0).sqrt() * t_add,
            "macaulay_censored_lower_bound": r["degrees"]["dstar_censored"].as_u64().unwrap_or(0) > 0,
        });
    }

    let failed: Vec<Value> = rows
        .iter()
        .filter(|r| !controls_ok(r))
        .map(|r| json!({"curve_class": r["curve_class"], "family": r["family"], "n": r["n"], "n_sub": r["n_sub"]}))
        .collect();
    let ok: Vec<&Value> = rows.iter().filter(|r| controls_ok(r)).collect();

    // H1: heuristic inside Bonferroni-adjusted Wilson interval.
    let z1 = norm_ppf(1.0 - 0.05 / (2.0 * ok.len().max(1) as f64));
    let mut h1_consistent = 0;
    let mut h1 = Vec::new();
    for r in &ok {
        let k = r["yield"]["decomposed"].as_f64().unwrap();
        let nn = r["yield"]["targets"].as_f64().unwrap();
        let h = r["yield"]["heuristic_V2_over_2n"].as_f64().unwrap();
        let (lo, hi) = wilson_z(k, nn, z1);
        let cons = h >= lo && h <= hi;
        if cons {
            h1_consistent += 1;
        }
        h1.push(json!({"curve_class": r["curve_class"], "family": r["family"], "n": r["n"], "n_sub": r["n_sub"],
            "measured": k / nn, "heuristic": h, "ratio": if h > 0.0 { Some(k / nn / h) } else { None },
            "interval_bonferroni": [lo, hi], "consistent": cons}));
    }

    let find = |class: &str, fam: &str, n: &Value, ns: &Value| -> Option<&Value> {
        ok.iter().copied().find(|r| {
            r["curve_class"] == class && r["family"] == fam && &r["n"] == n && &r["n_sub"] == ns
        })
    };

    // H2: family vs coordinate at matched (class, n, n').
    let mut h2_pairs: Vec<(&Value, &Value)> = Vec::new();
    for r in &ok {
        if r["family"] != "coordinate" {
            if let Some(c) = find(
                r["curve_class"].as_str().unwrap(),
                "coordinate",
                &r["n"],
                &r["n_sub"],
            ) {
                h2_pairs.push((r, c));
            }
        }
    }
    // H3: koblitz vs generic at matched (family, n, n').
    let mut h3_pairs: Vec<(&Value, &Value)> = Vec::new();
    for r in &ok {
        if r["curve_class"] == "koblitz_b1" {
            if let Some(g) = find(
                "generic_b",
                r["family"].as_str().unwrap(),
                &r["n"],
                &r["n_sub"],
            ) {
                h3_pairs.push((r, g));
            }
        }
    }
    let compare = |pairs: &[(&Value, &Value)], tag: u64| -> Vec<Value> {
        let alpha = 0.05 / pairs.len().max(1) as f64;
        let z = norm_ppf(1.0 - alpha / 2.0);
        pairs
            .iter()
            .enumerate()
            .map(|(i, (a, b))| {
                let xa = fvec(&a["degrees"]["dstar_values_censored_as_dmax_plus1"]);
                let xb = fvec(&b["degrees"]["dstar_values_censored_as_dmax_plus1"]);
                let ds = bootstrap_diff(&xa, &xb, alpha, 10_000, SEED ^ tag ^ i as u64);
                let (yd, ylo, yhi) = newcombe(
                    a["yield"]["decomposed"].as_f64().unwrap(), a["yield"]["targets"].as_f64().unwrap(),
                    b["yield"]["decomposed"].as_f64().unwrap(), b["yield"]["targets"].as_f64().unwrap(), z);
                json!({
                    "a": {"curve_class": a["curve_class"], "family": a["family"]},
                    "b": {"curve_class": b["curve_class"], "family": b["family"]},
                    "n": a["n"], "n_sub": a["n_sub"],
                    "dstar_diff": ds.map(|d| d.0), "dstar_interval": ds.map(|d| [d.1, d.2]),
                    "dstar_different": ds.map(|d| d.1 > 0.0 || d.2 < 0.0),
                    "yield_diff": yd, "yield_interval": [ylo, yhi], "yield_different": ylo > 0.0 || yhi < 0.0,
                    "bonferroni_alpha": alpha,
                })
            })
            .collect()
    };
    let h2 = compare(&h2_pairs, 0x22);
    let h3 = compare(&h3_pairs, 0x33);

    // H4: spread of d_ff vs spread of mean D* across n'.
    let mut h4 = Vec::new();
    let mut keys: Vec<(String, String, u64)> = ok
        .iter()
        .map(|r| {
            (
                r["curve_class"].as_str().unwrap().to_string(),
                r["family"].as_str().unwrap().to_string(),
                r["n"].as_u64().unwrap(),
            )
        })
        .collect();
    keys.sort();
    keys.dedup();
    for (c, f, n) in keys {
        let cells: Vec<&&Value> = ok
            .iter()
            .filter(|r| {
                r["curve_class"] == c.as_str()
                    && r["family"] == f.as_str()
                    && r["n"].as_u64() == Some(n)
            })
            .collect();
        if cells.len() < 2 {
            continue;
        }
        let ffm: Vec<f64> = cells
            .iter()
            .filter_map(|r| r["degrees"]["first_fall_mean"].as_f64())
            .collect();
        let dsm: Vec<f64> = cells
            .iter()
            .filter_map(|r| mean(&fvec(&r["degrees"]["dstar_values_censored_as_dmax_plus1"])))
            .collect();
        let range = |v: &[f64]| {
            v.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                - v.iter().cloned().fold(f64::INFINITY, f64::min)
        };
        h4.push(
            json!({"curve_class": c, "family": f, "n": n, "cells": cells.len(),
            "first_fall_mean_range": if ffm.len() > 1 { Some(range(&ffm)) } else { None },
            "dstar_mean_range": if dsm.len() > 1 { Some(range(&dsm)) } else { None }}),
        );
    }

    // H5: prime control.
    let mut h5 = Vec::new();
    let mut bits_list: Vec<u64> = prime.iter().filter_map(|r| r["bits"].as_u64()).collect();
    bits_list.sort();
    bits_list.dedup();
    for bits in bits_list {
        let rho: Vec<f64> = prime
            .iter()
            .filter(|r| {
                r["bits"].as_u64() == Some(bits)
                    && r["method"] == "rho_negation_1thread"
                    && r["solved"] == true
            })
            .filter_map(|r| r["seconds"].as_f64())
            .collect();
        let rho_max = rho.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mut cfgs: Vec<(String, u64)> = prime
            .iter()
            .filter(|r| r["bits"].as_u64() == Some(bits) && r["method"] != "rho_negation_1thread")
            .map(|r| {
                (
                    r["method"].as_str().unwrap().to_string(),
                    r["fb_size"].as_u64().unwrap(),
                )
            })
            .collect();
        cfgs.sort();
        cfgs.dedup();
        for (m, fb) in cfgs {
            let runs: Vec<&Value> = prime
                .iter()
                .filter(|r| {
                    r["bits"].as_u64() == Some(bits)
                        && r["method"] == m.as_str()
                        && r["fb_size"].as_u64() == Some(fb)
                })
                .collect();
            let solved: Vec<f64> = runs
                .iter()
                .filter(|r| r["solved"] == true)
                .filter_map(|r| r["seconds"].as_f64())
                .collect();
            let ic_min = solved.iter().cloned().fold(f64::INFINITY, f64::min);
            h5.push(json!({"bits": bits, "method": m, "fb_size": fb, "runs": runs.len(), "solved": solved.len(),
                "ic_fastest_s": if solved.is_empty() { None } else { Some(ic_min) },
                "rho_slowest_s": if rho.is_empty() { None } else { Some(rho_max) },
                "ic_slower_than_rho": if solved.is_empty() { None } else { Some(rho_max < ic_min) }}));
        }
    }

    json!({
        "cells": rows.len(),
        "cells_controls_ok": ok.len(),
        "cells_failed_controls": failed,
        "h1_yield_model": {"bonferroni_z": z1, "consistent_cells": h1_consistent, "cells": h1},
        "h2_family_vs_coordinate": h2,
        "h3_koblitz_vs_generic": h3,
        "h4_dff_vs_dstar_spread": h4,
        "h5_prime_control": h5,
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out_dir = PathBuf::from(args.get(1).map(String::as_str).unwrap_or("fb_sweep_out"));
    let quick = args.get(2).map(|s| s == "quick").unwrap_or(false);
    fs::create_dir_all(&out_dir).unwrap();

    let budget = if quick {
        BinaryBudget {
            yield_targets: 512,
            crosscheck_targets: 64,
            dstar_nondecomp: 4,
            dstar_decomp: 2,
            d_max: 6,
            curves: 2,
            planted: 16,
        }
    } else {
        BinaryBudget {
            yield_targets: 4096,
            crosscheck_targets: 256,
            dstar_nondecomp: 16,
            dstar_decomp: 4,
            d_max: 6,
            curves: 4,
            planted: 64,
        }
    };
    let (n_range, nsub_max): (Vec<u32>, u32) = if quick {
        (vec![8, 9], 4)
    } else {
        ((8..=16).collect(), 6)
    };

    let started = Instant::now();
    let mut rows = Vec::new();
    for &n in &n_range {
        for class in [CurveClass::Generic, CurveClass::Koblitz] {
            for family in [
                BasisFamily::Coordinate,
                BasisFamily::Random,
                BasisFamily::Subfield,
            ] {
                for n_sub in 2..=nsub_max.min(n - 1) {
                    if family == BasisFamily::Subfield && n % n_sub != 0 {
                        continue;
                    }
                    let seed = SEED
                        ^ ((n as u64) << 40)
                        ^ ((n_sub as u64) << 32)
                        ^ (family as u64) << 8
                        ^ class as u64;
                    let s = Instant::now();
                    if let Some(row) = binary_cell(class, family, n, n_sub, &budget, seed) {
                        eprintln!(
                            "binary {} {} n={n} n'={n_sub}: yield={:.4} D*={} ({:.1}s)",
                            class.name(),
                            family_name(family),
                            row["yield"]["rate"].as_f64().unwrap(),
                            row["degrees"]["dstar_mean"],
                            s.elapsed().as_secs_f64()
                        );
                        rows.push(row);
                    }
                }
            }
        }
    }
    let field_ops: Vec<Value> = n_range.iter().map(|&n| field_op_timings(n)).collect();

    let prime = if quick {
        prime_arm(&[20], &[64, 256], 1)
    } else {
        prime_arm(&[20, 24, 28], &[64, 256, 1024, 4096, 16384], 3)
    };

    let la_reps = 20_000_000u64;
    let s = Instant::now();
    let mut acc = 1u64;
    let p = (1u64 << 61) - 1;
    for i in 0..la_reps {
        acc = ((acc as u128 * 0x9E37_79B9_7F4A_7C15u128 + i as u128) % p as u128) as u64;
    }
    let la_op_s = s.elapsed().as_secs_f64() / la_reps as f64;

    let analysis = analyze(&mut rows, &prime, &field_ops, la_op_s);

    let doc = json!({
        "study": "fb_hyperparam_sweep_20261014",
        "seed": SEED,
        "quick": quick,
        "host_threads": std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
        "wall_seconds": started.elapsed().as_secs_f64(),
        "budget": {
            "yield_targets": budget.yield_targets,
            "crosscheck_targets": budget.crosscheck_targets,
            "dstar_nondecomp": budget.dstar_nondecomp,
            "dstar_decomp": budget.dstar_decomp,
            "d_max": budget.d_max,
            "generic_curves_per_cell": budget.curves,
            "planted_per_cell": budget.planted,
        },
        "field_ops": field_ops,
        "la_op_seconds": la_op_s,
        "la_op_sink": acc,
        "analysis": analysis,
        "binary": rows,
        "prime": prime,
    });
    let path = out_dir.join(if quick {
        "sweep_quick.json"
    } else {
        "sweep.json"
    });
    fs::write(&path, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
    eprintln!(
        "wrote {} in {:.1}s",
        path.display(),
        started.elapsed().as_secs_f64()
    );
}
