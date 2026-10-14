//! lit2026 W8: summation-polynomial model checks.
//!
//! Four diagnostics from the 2025–26 literature search
//! (`research/literature_search_20261011/`), each run on the repository's own
//! Weil-descended decomposition systems and solved by the same matrix-F4 with
//! splitting (`koblitz_groebner::solve_boolean_system`, all roots enumerated):
//!
//! * `eg` — Dryło, ePrint 2026/1299: a `G`-invariant compression is `f∘Φ` with
//!   `Φ: E → E/G`.  On a Koblitz curve the FGHR invariant `w = u² + u`,
//!   `u = 1/(x + 1)`, is `1/X∘Φ` with `X` the Vélu abscissa of `E/⟨T⟩`.  Arms on
//!   the same factor base `F_u`: `u-plain` (Semaev in `u`, no quotient), `eg`
//!   (Semaev of `E/⟨T⟩` in `W = 1/X`), `sym` (the repo's symmetrised system:
//!   `eg` plus the parity cover `s = Σu`), and `x-direct` (Semaev in `x` on the
//!   Frobenius base of the same divisor, the repo's context arm).
//! * `wu` — Wu, arXiv 2608.01675: `S₃` in the native Kummer coordinates of the
//!   binary state model `(u² + u)(v² + v) = d` (`u` and the Kummer ratio
//!   `k = (u + 1)/u`), against Weierstrass `x` and binary Edwards `w = x + y`.
//! * `frob` — Bariant, ePrint 2026/1927: the Frobenius auxiliary variable
//!   `y = x^{2^k}` added to the descended `S_{m+1}` on a Frobenius-stable base.
//! * `lucas` — La Scala–Marchesin–Tiwari, arXiv 2609.18866: convolution
//!   splitting of `e_i(X₁, …, X_m)` with auxiliaries, against the repo's
//!   `binary_semaev_s4` correspondence.
//! * `cert` — emit the `identity.certificate/v1` record for
//!   `Res_s(S^sym_{m+1}, s² + s + Σw) = S^{E/⟨T⟩}_{m+1}(1/w)·Πw^d` over `F₂`.
//!
//! ```bash
//! cargo run --release --example lit2026_w8_summation_models -- eg    --rounds 5 --out eg.json
//! cargo run --release --example lit2026_w8_summation_models -- wu    --rounds 5 --out wu.json
//! cargo run --release --example lit2026_w8_summation_models -- frob  --rounds 5 --out frob.json
//! cargo run --release --example lit2026_w8_summation_models -- lucas --rounds 5 --out lucas.json
//! cargo run --release --example lit2026_w8_summation_models -- cert  --out cert.json
//! ```
//!
//! Every arm's root set is checked: `eg` against `sym` (projection dropping
//! the parity bit) and `u-plain` (image under `AS`), `wu` against brute-force
//! enumeration of the polynomial over `V × V`, `frob` against the plain
//! system's roots, `lucas` against brute force over `V³`.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::time::Instant;

use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

use crypto_lib::binary_ecc::curve::{point_add, point_neg};
use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_semaev_s4::{
    elementary_symmetric_3, symmetrised_s4_eval, weil_descend_s4, AnfF2m, AnfPoly,
};
use crypto_lib::cryptanalysis::identity_certificate::{
    check, issue, syntactic_degree, Expr, IdentityCertificate, IdentityStatement, Verdict,
};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    f4_word_ops_thread, first_fall_degree, solve_boolean_system, solving_degree,
    split_rule_default, system_degree, FieldStructure, SolveOptions, SolveStats, SolverEngine,
    SymElement,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_frobenius_factor_base_from_divisor, points_with_x, KoblitzCurve,
};
use crypto_lib::cryptanalysis::koblitz_symmetrised::{
    artin_schreier, build_direct_x_system, build_symmetrised_factor_base, build_symmetrised_system,
    divisor_for_dimension, enumerate_symmetrised, lift_symmetrised_root, plain_terms,
    symmetrised_terms, u_of_x, SymmetrisedFactorBase,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};

// ════════════════════════════════════════════════════════════════════
// Field helpers
// ════════════════════════════════════════════════════════════════════

fn bits(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

fn fe(v: u64, n: u32) -> F2mElement {
    F2mElement::from_biguint(&BigUint::from(v), n)
}

fn fpow(a: &F2mElement, e: u32, irr: &IrreduciblePoly) -> F2mElement {
    let n = a.m_value();
    let mut r = F2mElement::one(n);
    for _ in 0..e {
        r = r.mul(a, irr);
    }
    r
}

/// Absolute trace `Σ_{i<n} z^{2^i}` as a bit.
#[allow(dead_code)]
fn trace(z: &F2mElement, irr: &IrreduciblePoly) -> bool {
    let n = z.m_value();
    let mut acc = F2mElement::zero(n);
    let mut t = z.clone();
    for _ in 0..n {
        acc = acc.add(&t);
        t = t.square(irr);
    }
    !acc.is_zero()
}

/// A root of `s² + s = c` for odd `n` (half-trace), if `Tr(c) = 0`.
fn solve_as(c: &F2mElement, irr: &IrreduciblePoly) -> Option<F2mElement> {
    let n = c.m_value();
    assert!(n % 2 == 1, "half-trace needs odd n");
    let mut acc = F2mElement::zero(n);
    let mut t = c.clone();
    for _ in 0..=(n - 1) / 2 {
        acc = acc.add(&t);
        t = t.square(irr).square(irr);
    }
    (acc.square(irr).add(&acc) == *c).then_some(acc)
}

fn fsqrt(c: &F2mElement, irr: &IrreduciblePoly) -> F2mElement {
    c.square_k_times(c.m_value() - 1, irr)
}

// ════════════════════════════════════════════════════════════════════
// Polynomials over F_{2^n} with exponent vectors, and Möbius substitution
// ════════════════════════════════════════════════════════════════════

/// `Σ coef · Π v_i^{e_i}`, keyed by exponent vector.
type FPoly = BTreeMap<Vec<u32>, F2mElement>;

fn fpoly_add_term(p: &mut FPoly, e: Vec<u32>, c: F2mElement) {
    if c.is_zero() {
        return;
    }
    let n = c.m_value();
    let s = p.remove(&e).unwrap_or_else(|| F2mElement::zero(n)).add(&c);
    if !s.is_zero() {
        p.insert(e, s);
    }
}

/// Coefficient-one polynomial from exponent vectors (a constant term of `b`
/// when `const_b` is given).
fn fpoly_from_terms(terms: &[Vec<u32>], n: u32, const_b: Option<&F2mElement>) -> FPoly {
    let mut p = FPoly::new();
    for t in terms {
        let c = if t.iter().all(|&e| e == 0) {
            const_b.cloned().unwrap_or_else(|| F2mElement::one(n))
        } else {
            F2mElement::one(n)
        };
        fpoly_add_term(&mut p, t.clone(), c);
    }
    p
}

/// `(a·c + b)^e` as coefficients of `c^0 … c^e`.
fn upow(a: &F2mElement, b: &F2mElement, e: u32, irr: &IrreduciblePoly) -> Vec<F2mElement> {
    let n = a.m_value();
    let mut r = vec![F2mElement::one(n)];
    for _ in 0..e {
        let mut next = vec![F2mElement::zero(n); r.len() + 1];
        for (k, c) in r.iter().enumerate() {
            next[k] = next[k].add(&c.mul(b, irr));
            next[k + 1] = next[k + 1].add(&c.mul(a, irr));
        }
        r = next;
    }
    r
}

fn umul(p: &[F2mElement], q: &[F2mElement], irr: &IrreduciblePoly) -> Vec<F2mElement> {
    let n = p[0].m_value();
    let mut r = vec![F2mElement::zero(n); p.len() + q.len() - 1];
    for (i, a) in p.iter().enumerate() {
        for (j, b) in q.iter().enumerate() {
            r[i + j] = r[i + j].add(&a.mul(b, irr));
        }
    }
    r
}

/// Substitute `x_i = (A c_i + B)/(C c_i + D)` (`maps[i] = [A, B, C, D]`) and
/// clear denominators by `Π (C c_i + D)^{deg_i}`.
fn mobius_substitute(
    p: &FPoly,
    degs: &[u32],
    maps: &[[F2mElement; 4]],
    irr: &IrreduciblePoly,
) -> FPoly {
    let mut out = FPoly::new();
    for (exps, coef) in p {
        let factors: Vec<Vec<F2mElement>> = exps
            .iter()
            .enumerate()
            .map(|(i, &e)| {
                let [a, b, c, d] = &maps[i];
                umul(&upow(a, b, e, irr), &upow(c, d, degs[i] - e, irr), irr)
            })
            .collect();
        // Cartesian product.
        let mut acc: Vec<(Vec<u32>, F2mElement)> = vec![(Vec::new(), coef.clone())];
        for f in &factors {
            let mut next = Vec::new();
            for (ev, c) in &acc {
                for (k, fc) in f.iter().enumerate() {
                    if fc.is_zero() {
                        continue;
                    }
                    let mut e2 = ev.clone();
                    e2.push(k as u32);
                    next.push((e2, c.mul(fc, irr)));
                }
            }
            acc = next;
        }
        for (e, c) in acc {
            fpoly_add_term(&mut out, e, c);
        }
    }
    out
}

fn fpoly_eval(p: &FPoly, vals: &[F2mElement], irr: &IrreduciblePoly) -> F2mElement {
    let n = vals[0].m_value();
    let mut acc = F2mElement::zero(n);
    for (e, c) in p {
        let mut t = c.clone();
        for (i, &k) in e.iter().enumerate() {
            if k > 0 {
                t = t.mul(&fpow(&vals[i], k, irr), irr);
            }
        }
        acc = acc.add(&t);
    }
    acc
}

#[derive(Clone, Debug, serde::Serialize)]
struct PolyShape {
    terms: usize,
    total_degree: u32,
    per_var_degree: Vec<u32>,
    f2_coefficients: bool,
}

fn poly_shape(p: &FPoly) -> PolyShape {
    let nv = p.keys().next().map(|e| e.len()).unwrap_or(0);
    let mut per = vec![0u32; nv];
    let mut total = 0;
    let mut f2 = true;
    for (e, c) in p {
        total = total.max(e.iter().sum());
        for (i, &k) in e.iter().enumerate() {
            per[i] = per[i].max(k);
        }
        if bits(c) != 1 {
            f2 = false;
        }
    }
    PolyShape {
        terms: p.len(),
        total_degree: total,
        per_var_degree: per,
        f2_coefficients: f2,
    }
}

// ════════════════════════════════════════════════════════════════════
// Weil descent of an FPoly whose last variable is the known target
// ════════════════════════════════════════════════════════════════════

fn sym_pow_cached(
    cache: &mut HashMap<(usize, u32), SymElement>,
    vars: &[SymElement],
    i: usize,
    e: u32,
    st: &FieldStructure,
) -> SymElement {
    if let Some(v) = cache.get(&(i, e)) {
        return v.clone();
    }
    let v = if e == 1 {
        vars[i].clone()
    } else if e.is_multiple_of(2) {
        sym_pow_cached(cache, vars, i, e / 2, st).square(st)
    } else {
        sym_pow_cached(cache, vars, i, e - 1, st).mul(&vars[i], st)
    };
    cache.insert((i, e), v.clone());
    v
}

/// Descend `p(v_0, …, v_{k−1}, t)` with `t = target` known: `n` Boolean
/// equations over the variables of `vars` (one symbolic element per unknown).
fn descend(
    p: &FPoly,
    vars: &[SymElement],
    target: &F2mElement,
    st: &FieldStructure,
    irr: &IrreduciblePoly,
    n_vars: usize,
) -> Vec<F2BoolPoly> {
    let n = st.n;
    let k = vars.len();
    // Fold the target into the coefficients.
    let mut folded = FPoly::new();
    for (e, c) in p {
        let c2 = c.mul(&fpow(target, e[k], irr), irr);
        fpoly_add_term(&mut folded, e[..k].to_vec(), c2);
    }
    let mut cache = HashMap::new();
    let mut acc = SymElement::zero(n, n_vars);
    for (e, c) in &folded {
        let mut prod: Option<SymElement> = None;
        for (i, &ei) in e.iter().enumerate() {
            if ei == 0 {
                continue;
            }
            let pw = sym_pow_cached(&mut cache, vars, i, ei, st);
            prod = Some(match prod {
                None => pw,
                Some(q) => q.mul(&pw, st),
            });
        }
        let term = match prod {
            None => SymElement::constant(c, n, n_vars),
            Some(q) if bits(c) == 1 => q,
            Some(q) => q.mul_constant(c, st),
        };
        acc = acc.add(&term);
    }
    acc.coords.into_iter().filter(|p| !p.is_zero()).collect()
}

#[derive(Clone, Debug, serde::Serialize)]
struct SystemShape {
    n_vars: usize,
    n_equations: usize,
    boolean_degree: u32,
    distinct_monomials: usize,
    total_terms: usize,
}

fn system_shape(eqs: &[F2BoolPoly], n_vars: usize) -> SystemShape {
    let mut monos = HashSet::new();
    let mut total = 0;
    for p in eqs {
        for t in &p.terms {
            monos.insert(t.mask);
            total += 1;
        }
    }
    SystemShape {
        n_vars,
        n_equations: eqs.iter().filter(|p| !p.is_zero()).count(),
        boolean_degree: system_degree(eqs),
        distinct_monomials: monos.len(),
        total_terms: total,
    }
}

#[derive(Clone, Debug, serde::Serialize)]
struct Degrees {
    /// Lowest Macaulay degree at which bounded multiples refute the system
    /// or pin every occurring variable (`koblitz_groebner::solving_degree`);
    /// `None` = not resolved up to `scan_max` (expected for systems with
    /// several roots, see that function's docs).
    resolution_degree: Option<u32>,
    first_fall_degree: Option<u32>,
    scan_max: u32,
    /// Degrees at which the resolution scan actually built a Macaulay
    /// matrix; a scan that stops short of `scan_max` without resolving hit
    /// the size caps.
    resolution_degrees_built: Vec<u32>,
}

fn degrees(eqs: &[F2BoolPoly], n_vars: usize, dmax: u32) -> Degrees {
    let (sd, profiles) = solving_degree(eqs, n_vars, dmax);
    let (ffd, _) = first_fall_degree(eqs, n_vars, dmax);
    Degrees {
        resolution_degree: sd,
        first_fall_degree: ffd,
        scan_max: dmax,
        resolution_degrees_built: profiles.iter().map(|p| p.degree).collect(),
    }
}

/// `--engine inherited` selects `InheritedF4` (the repository default, which
/// also eliminates linear generators); the default here is `MatrixF4`.
static INHERITED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

/// `--arms a,b` restricts the `eg` part to the named arms (for running the
/// `original` arm from a binary built against the program-baseline library);
/// `--no-degrees` skips the untimed degree scans.
static ARMS: std::sync::OnceLock<Option<Vec<String>>> = std::sync::OnceLock::new();
static NO_DEGREES: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

fn arm_selected(name: &str) -> bool {
    match ARMS.get().and_then(|a| a.as_ref()) {
        None => true,
        Some(list) => list.iter().any(|a| a == name),
    }
}

fn degrees_enabled() -> bool {
    !*NO_DEGREES.get().unwrap_or(&false)
}

/// FNV-1a over the sorted key set: lets two processes compare root sets.
fn key_digest(keys: &BTreeSet<Vec<u64>>) -> String {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for k in keys {
        for &v in k {
            for b in v.to_le_bytes() {
                h ^= b as u64;
                h = h.wrapping_mul(0x0100_0000_01b3);
            }
        }
        h ^= 0xff;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    format!("{h:016x}:{}", keys.len())
}

fn engine_for(max_degree: u32) -> SolverEngine {
    if *INHERITED.get().unwrap_or(&false) {
        SolverEngine::InheritedF4 { max_degree }
    } else {
        SolverEngine::MatrixF4 { max_degree }
    }
}

fn engine_label(max_degree: u32) -> String {
    format!("{:?}", engine_for(max_degree))
}

#[derive(Clone, Debug, serde::Serialize)]
struct SolveSample {
    ms: f64,
    word_ops: u64,
    splits: usize,
    reductions: usize,
    max_degree_built: u32,
    oversize: usize,
    exhausted: bool,
    roots: usize,
}

fn timed_solve(
    eqs: &[F2BoolPoly],
    n_vars: usize,
    max_degree: u32,
    node_budget: usize,
) -> (Vec<u64>, SolveSample) {
    let opts = SolveOptions {
        engine: engine_for(max_degree),
        max_solutions: usize::MAX,
        node_budget,
        split_rule: split_rule_default(),
    };
    let w0 = f4_word_ops_thread();
    let t0 = Instant::now();
    let (mut roots, stats): (Vec<u64>, SolveStats) = solve_boolean_system(eqs, n_vars, &opts);
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    let w = f4_word_ops_thread() - w0;
    roots.sort_unstable();
    roots.dedup();
    let s = SolveSample {
        ms,
        word_ops: w,
        splits: stats.splits,
        reductions: stats.reductions,
        max_degree_built: stats.max_degree_built,
        oversize: stats.oversize,
        exhausted: stats.exhausted,
        roots: roots.len(),
    };
    (roots, s)
}

fn uptime() -> String {
    std::process::Command::new("uptime")
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default()
}

fn span_elem(basis: &[F2mElement], root: u64, offset: usize, n: u32) -> F2mElement {
    let mut acc = F2mElement::zero(n);
    for (t, b) in basis.iter().enumerate() {
        if (root >> (offset + t)) & 1 == 1 {
            acc = acc.add(b);
        }
    }
    acc
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if s.is_empty() {
        f64::NAN
    } else {
        s[s.len() / 2]
    }
}

// ════════════════════════════════════════════════════════════════════
// Part A — E/G check of the FGHR 2-torsion symmetrisation
// ════════════════════════════════════════════════════════════════════

/// `X = x + √b/x`, Vélu's abscissa of `E/⟨T⟩`, `T = (0, √b)`.
fn velu_x(x: &F2mElement, b: &F2mElement, irr: &IrreduciblePoly) -> Option<F2mElement> {
    Some(x.add(&fsqrt(b, irr).mul(&x.flt_inverse(irr)?, irr)))
}

/// `b''` of `E/⟨T⟩` read off one relation `P₁ + P₂ + P₃ = O` through
/// `S₃ = e₂(X)² + e₃(X) + b''`, then checked on `checks` further relations.
/// `None` if any check fails.
fn quotient_b(curve: &BinaryCurve, seed: u64, checks: usize) -> Option<F2mElement> {
    let irr = &curve.irreducible;
    let n = curve.m;
    let mut rng = StdRng::seed_from_u64(seed);
    let rand_point = |rng: &mut StdRng| loop {
        let x = fe(rng.gen_range(2..(1u64 << n)), n);
        if let Some(p) = points_with_x(curve, &x).into_iter().next() {
            return p;
        }
    };
    let mut found: Option<F2mElement> = None;
    let mut done = 0;
    while done <= checks {
        let p1 = rand_point(&mut rng);
        let p2 = rand_point(&mut rng);
        let p3 = point_neg(&point_add(curve, &p1, &p2));
        let xs: Vec<F2mElement> = [p1, p2, p3]
            .iter()
            .filter_map(|p| match p {
                BinaryPoint::Affine { x, .. } => velu_x(x, &curve.b, irr),
                _ => None,
            })
            .collect();
        if xs.len() != 3 {
            continue;
        }
        let e2 = xs[0]
            .mul(&xs[1], irr)
            .add(&xs[0].mul(&xs[2], irr))
            .add(&xs[1].mul(&xs[2], irr));
        let e3 = xs[0].mul(&xs[1], irr).mul(&xs[2], irr);
        let bpp = e2.square(irr).add(&e3);
        match &found {
            None => found = Some(bpp),
            Some(f) if *f != bpp => return None,
            _ => {}
        }
        done += 1;
    }
    found
}

/// `S^{E/⟨T⟩}_{m+1}` in `W = 1/X` (Koblitz, `b'' = 1`): Möbius `X = 1/W`.
fn eg_poly(m: usize, n: u32, irr: &IrreduciblePoly) -> FPoly {
    let d = 1u32 << (m - 1);
    let p = fpoly_from_terms(&plain_terms(m).unwrap(), n, None);
    let one = F2mElement::one(n);
    let zero = F2mElement::zero(n);
    let map = [zero.clone(), one.clone(), one.clone(), zero.clone()];
    mobius_substitute(&p, &vec![d; m + 1], &vec![map; m + 1], irr)
}

/// Semaev `S_{m+1}` (`b = 1`) in `u = 1/(x + 1)`: `x = (u + 1)/u`.
fn u_plain_poly(m: usize, n: u32, irr: &IrreduciblePoly) -> FPoly {
    let d = 1u32 << (m - 1);
    let p = fpoly_from_terms(&plain_terms(m).unwrap(), n, None);
    let one = F2mElement::one(n);
    let zero = F2mElement::zero(n);
    let map = [one.clone(), one.clone(), one.clone(), zero.clone()];
    mobius_substitute(&p, &vec![d; m + 1], &vec![map; m + 1], irr)
}

/// The symmetrised polynomial as an [`FPoly`] over `(w₁..w_m, w_R, s)`.
fn sym_fpoly(m: usize, n: u32) -> FPoly {
    fpoly_from_terms(&symmetrised_terms(m).unwrap(), n, None)
}

struct ArmSys {
    name: &'static str,
    eqs: Vec<F2BoolPoly>,
    n_vars: usize,
    /// Root → comparison key (ordered tuple of field-element bit patterns).
    keys: Box<dyn Fn(u64) -> Vec<u64>>,
    /// Whether the keys are comparable across arms of this cell.
    comparable: bool,
}

#[derive(Clone, Debug, serde::Serialize)]
struct ArmReport {
    arm: String,
    poly: Option<PolyShape>,
    system: SystemShape,
    degrees: Vec<Degrees>,
    samples: Vec<Vec<SolveSample>>, // [target][round]
    median_ms: f64,
    median_word_ops: f64,
    key_sets_identical_to_reference: Vec<bool>,
    /// Per target: digest of the arm's `w`-tuple key set (comparable across
    /// processes and binaries for the `sym`, `eg` and `u-plain` arms).
    key_set_digests: Vec<String>,
}

fn run_eg_cell(
    a: u8,
    n: u32,
    m: usize,
    targets: usize,
    rounds: usize,
    seed: u64,
    max_degree: u32,
    dmax: u32,
    node_budget: usize,
    with_x_direct: bool,
) -> Option<Value> {
    let kc = KoblitzCurve::new(a, n)?;
    let irr = kc.curve.irreducible.clone();
    let st = FieldStructure::new(n, &irr);
    let target_dim = (n + 1).div_ceil(m as u32);
    let divisor = divisor_for_dimension(n, target_dim)?;
    let fb: SymmetrisedFactorBase = build_symmetrised_factor_base(&kc, &divisor)?;
    let fb_x = build_frobenius_factor_base_from_divisor(&kc, &divisor)?;
    let ell = fb.ell;
    let bpp = quotient_b(&kc.curve, seed ^ 0xE6, 32)?;

    let p_eg = eg_poly(m, n, &irr);
    let p_u = u_plain_poly(m, n, &irr);
    let p_sym = sym_fpoly(m, n);
    let p_x = fpoly_from_terms(&plain_terms(m).unwrap(), n, None);

    let mut rng = StdRng::seed_from_u64(seed);
    let g = kc.generator().clone();
    let r = kc
        .subgroup_order
        .to_u64_digits()
        .first()
        .copied()
        .unwrap_or(2)
        .max(2);
    let mut tgts = Vec::new();
    while tgts.len() < targets {
        let p = kc.mul(&g, &BigUint::from(rng.gen_range(1..r)));
        if let BinaryPoint::Affine { x, .. } = &p {
            if *x == F2mElement::one(n) || x.is_zero() {
                continue;
            }
            tgts.push(p);
        }
    }

    let mut reports: Vec<ArmReport> = Vec::new();
    let mut truths = Vec::new();
    let mut verdict_ok = Vec::new();
    let mut uptimes = Vec::new();

    // Build every arm for every target first (not timed).
    let mut cells: Vec<Vec<ArmSys>> = Vec::new();
    for t in &tgts {
        let x_r = match t {
            BinaryPoint::Affine { x, .. } => x.clone(),
            _ => unreachable!(),
        };
        let u_r = u_of_x(&x_r, &kc)?;
        let w_r = artin_schreier(&u_r, &kc);
        let mut arms = Vec::new();
        // sym
        let sys = build_symmetrised_system(&kc, &fb, t, m, &st)?;
        let wb = fb.w_basis.clone();
        let lw = ell - 1;
        arms.push(ArmSys {
            name: "sym",
            eqs: sys.equations.clone(),
            n_vars: sys.n_vars,
            keys: Box::new(move |root| {
                (0..m)
                    .map(|i| bits(&span_elem(&wb, root, i * lw, n)))
                    .collect()
            }),
            comparable: true,
        });
        // eg
        let nv = m * lw;
        let vars: Vec<SymElement> = (0..m)
            .map(|i| SymElement::from_subspace_vars(&fb.w_basis, i * lw, n, nv))
            .collect();
        let eqs = descend(&p_eg, &vars, &w_r, &st, &irr, nv);
        let wb = fb.w_basis.clone();
        arms.push(ArmSys {
            name: "eg",
            eqs,
            n_vars: nv,
            keys: Box::new(move |root| {
                (0..m)
                    .map(|i| bits(&span_elem(&wb, root, i * lw, n)))
                    .collect()
            }),
            comparable: true,
        });
        // u-plain
        let nv = m * ell;
        let vars: Vec<SymElement> = (0..m)
            .map(|i| SymElement::from_subspace_vars(&fb.v_basis, i * ell, n, nv))
            .collect();
        let eqs = descend(&p_u, &vars, &u_r, &st, &irr, nv);
        let vb = fb.v_basis.clone();
        let irr2 = irr.clone();
        arms.push(ArmSys {
            name: "u-plain",
            eqs,
            n_vars: nv,
            keys: Box::new(move |root| {
                (0..m)
                    .map(|i| {
                        let u = span_elem(&vb, root, i * ell, n);
                        bits(&u.square(&irr2).add(&u))
                    })
                    .collect()
            }),
            comparable: true,
        });
        if with_x_direct {
            let (eqs, nv) = build_direct_x_system(&fb_x.subspace_basis, &x_r, m, &st)?;
            let xb = fb_x.subspace_basis.clone();
            let lx = xb.len();
            arms.push(ArmSys {
                name: "x-direct",
                eqs,
                n_vars: nv,
                keys: Box::new(move |root| {
                    (0..m)
                        .map(|i| bits(&span_elem(&xb, root, i * lx, n)))
                        .collect()
                }),
                comparable: false,
            });
        }
        // Truth on F_u.
        truths.push(enumerate_symmetrised(&kc, &fb, t, m).is_some());
        arms.retain(|a| arm_selected(a.name));
        cells.push(arms);
    }
    let n_arms = cells[0].len();
    for ai in 0..n_arms {
        let arm = &cells[0][ai];
        let poly = match arm.name {
            "sym" => Some(poly_shape(&p_sym)),
            "eg" => Some(poly_shape(&p_eg)),
            "u-plain" => Some(poly_shape(&p_u)),
            "x-direct" => Some(poly_shape(&p_x)),
            _ => None,
        };
        reports.push(ArmReport {
            arm: arm.name.to_string(),
            poly,
            system: system_shape(&arm.eqs, arm.n_vars),
            degrees: if degrees_enabled() {
                cells
                    .iter()
                    .map(|c| degrees(&c[ai].eqs, c[ai].n_vars, dmax))
                    .collect()
            } else {
                Vec::new()
            },
            samples: vec![Vec::new(); tgts.len()],
            median_ms: 0.0,
            median_word_ops: 0.0,
            key_sets_identical_to_reference: Vec::new(),
            key_set_digests: Vec::new(),
        });
    }

    // Interleaved rounds: arm order rotates with the round.
    let mut ref_keys: Vec<Option<BTreeSet<Vec<u64>>>> = vec![None; tgts.len()];
    let mut arm_keys: Vec<Vec<Option<BTreeSet<Vec<u64>>>>> = vec![vec![None; tgts.len()]; n_arms];
    let mut root_stable = true;
    for round in 0..rounds {
        uptimes.push(uptime());
        for (ti, arms) in cells.iter().enumerate() {
            for k in 0..n_arms {
                let ai = (k + round) % n_arms;
                let arm = &arms[ai];
                let (roots, s) = timed_solve(&arm.eqs, arm.n_vars, max_degree, node_budget);
                let ks: BTreeSet<Vec<u64>> = roots.iter().map(|&r| (arm.keys)(r)).collect();
                if arm_keys[ai][ti].is_none() {
                    arm_keys[ai][ti] = Some(ks.clone());
                } else if arm_keys[ai][ti].as_ref() != Some(&ks) {
                    root_stable = false;
                }
                if arm.name == "sym" && ref_keys[ti].is_none() {
                    ref_keys[ti] = Some(ks);
                    // Verdict gate: does some sym root lift?
                    let t = &tgts[ti];
                    let sys = build_symmetrised_system(&kc, &fb, t, m, &st).unwrap();
                    let lifts = roots
                        .iter()
                        .any(|&rt| lift_symmetrised_root(&kc, &fb, &sys, rt, t).is_some());
                    verdict_ok.push(lifts == truths[ti]);
                }
                reports[ai].samples[ti].push(s);
            }
        }
    }
    for (ai, rep) in reports.iter_mut().enumerate() {
        let comparable = cells[0][ai].comparable;
        rep.key_sets_identical_to_reference = (0..tgts.len())
            .map(|ti| !comparable || ref_keys[ti].is_none() || arm_keys[ai][ti] == ref_keys[ti])
            .collect();
        rep.key_set_digests = (0..tgts.len())
            .map(|ti| {
                arm_keys[ai][ti]
                    .as_ref()
                    .map(key_digest)
                    .unwrap_or_default()
            })
            .collect();
        let ms: Vec<f64> = rep.samples.iter().flatten().map(|s| s.ms).collect();
        let wo: Vec<f64> = rep
            .samples
            .iter()
            .flatten()
            .map(|s| s.word_ops as f64)
            .collect();
        rep.median_ms = median(&ms);
        rep.median_word_ops = median(&wo);
    }
    Some(json!({
        "cell": format!("K{}/F2^{} m={}", a, n, m),
        "a": a, "n": n, "m": m, "ell": ell,
        "subgroup_order": kc.subgroup_order.to_string(),
        "cofactor": kc.cofactor.to_string(),
        "quotient_b_double_prime": bits(&bpp),
        "targets": tgts.len(),
        "rounds": rounds,
        "engine": engine_label(max_degree),
        "node_budget": node_budget,
        "truth_decomposable": truths,
        "sym_verdict_matches_enumeration": verdict_ok,
        "roots_stable_across_rounds": root_stable,
        "uptime_per_round": uptimes,
        "arms": reports,
    }))
}

// ════════════════════════════════════════════════════════════════════
// Part B — Wu's binary state model and binary Edwards, S₃ only
// ════════════════════════════════════════════════════════════════════

/// Solve the homogeneous `rows × 4` system over `F_{2^n}` for a nonzero kernel
/// vector (Gaussian elimination, last free column set to 1).
/// Reduced row echelon over `F_{2^n}`; returns the nullity and one kernel
/// vector (last free column set to 1) when the nullity is positive.
fn kernel_vec(
    mut rows: Vec<Vec<F2mElement>>,
    ncols: usize,
    irr: &IrreduciblePoly,
) -> (usize, Option<Vec<F2mElement>>) {
    let n = rows[0][0].m_value();
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..ncols {
        let Some(p) = (r..rows.len()).find(|&i| !rows[i][c].is_zero()) else {
            continue;
        };
        rows.swap(r, p);
        let inv = rows[r][c].flt_inverse(irr).unwrap();
        for k in 0..ncols {
            rows[r][k] = rows[r][k].mul(&inv, irr);
        }
        for i in 0..rows.len() {
            if i != r && !rows[i][c].is_zero() {
                let f = rows[i][c].clone();
                for k in 0..ncols {
                    let t = rows[r][k].mul(&f, irr);
                    rows[i][k] = rows[i][k].add(&t);
                }
            }
        }
        pivots.push(c);
        r += 1;
    }
    let nullity = ncols - pivots.len();
    let Some(free) = (0..ncols).rev().find(|c| !pivots.contains(c)) else {
        return (0, None);
    };
    let mut v = vec![F2mElement::zero(n); ncols];
    v[free] = F2mElement::one(n);
    for (ri, &pc) in pivots.iter().enumerate() {
        v[pc] = rows[ri][free].clone();
    }
    (nullity, Some(v))
}

/// A binary Edwards model `d₁(x+y) + d₂(x²+y²) = xy + xy(x+y) + x²y²` of
/// `y² + xy = x³ + a x² + b` (Bernstein–Lange–Rezaeian Farashahi 2008):
/// returns `(d₁, d₂, s)` with the Weierstrass curve of `(d₁, d₂)` equal to
/// `y² + xy = x³ + (d₁² + d₂)x² + b` and `s² + s = a + d₁² + d₂`.
fn binary_edwards_params(
    curve: &BinaryCurve,
    seed: u64,
) -> Option<(F2mElement, F2mElement, F2mElement)> {
    let irr = &curve.irreducible;
    let n = curve.m;
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..10_000 {
        let d1 = fe(rng.gen_range(2..(1u64 << n)), n);
        let d1_2 = d1.square(irr);
        let d1_4 = d1_2.square(irr);
        let inv = d1_4.flt_inverse(irr)?;
        let d2_sq = curve.b.mul(&inv, irr).add(&d1_4).add(&d1_2);
        let d2 = fsqrt(&d2_sq, irr);
        if d2 == d1_2.add(&d1) {
            continue;
        }
        let a2 = d1_2.add(&d2);
        if let Some(s) = solve_as(&curve.a.add(&a2), irr) {
            return Some((d1, d2, s));
        }
    }
    None
}

/// Map `(x, y)` on `y² + xy = x³ + ax² + b` to the binary Edwards model and
/// return `(x_e, y_e)`, checking the Edwards equation.
fn to_edwards(
    p: &BinaryPoint,
    d1: &F2mElement,
    d2: &F2mElement,
    s: &F2mElement,
    irr: &IrreduciblePoly,
) -> Option<(F2mElement, F2mElement)> {
    let BinaryPoint::Affine { x, y } = p else {
        return None;
    };
    let u = x.clone();
    let v = y.add(&s.mul(x, irr));
    let d1_2 = d1.square(irr);
    let k = d1_2.add(d1).add(d2); // d1² + d1 + d2
    let c = d1_2.add(d1).mul(&k, irr); // (d1² + d1)(d1² + d1 + d2)
    let num = d1.mul(&u.add(&k), irr);
    let xe = num.mul(&u.add(&v).add(&c).flt_inverse(irr)?, irr);
    let ye = num.mul(&v.add(&c).flt_inverse(irr)?, irr);
    // Check the Edwards equation.
    let lhs = d1
        .mul(&xe.add(&ye), irr)
        .add(&d2.mul(&xe.square(irr).add(&ye.square(irr)), irr));
    let xy = xe.mul(&ye, irr);
    let rhs = xy.add(&xy.mul(&xe.add(&ye), irr)).add(&xy.square(irr));
    (lhs == rhs).then_some((xe, ye))
}

struct Model {
    name: &'static str,
    /// The coordinate is a Möbius function of the Vélu abscissa
    /// `X = x + √b/x` of `E/⟨T⟩` (a Dryło `G`-compression) rather than of
    /// `x`; its `S₃` is then the quotient curve's, with `b''`.
    quotient: bool,
    /// `base = (A c + B)/(C c + D)`, `base = x` or `X`.
    to_x: [F2mElement; 4],
    /// `c = (A' base + B')/(C' base + D')`.
    from_x: [F2mElement; 4],
}

/// The base abscissa of a model at `x`: `x` itself, or `X = x + √b/x`.
fn model_base(
    md: &Model,
    x: &F2mElement,
    b: &F2mElement,
    irr: &IrreduciblePoly,
) -> Option<F2mElement> {
    if md.quotient {
        velu_x(x, b, irr)
    } else {
        Some(x.clone())
    }
}

/// Abscissae `x` over a base value: `[x]`, or the roots of `x² + Xx + √b`.
fn base_to_xs(
    md: &Model,
    v: &F2mElement,
    b: &F2mElement,
    irr: &IrreduciblePoly,
) -> Vec<F2mElement> {
    if !md.quotient {
        return vec![v.clone()];
    }
    let Some(inv) = v.flt_inverse(irr) else {
        return Vec::new();
    };
    let c = fsqrt(b, irr).mul(&inv.square(irr), irr);
    match solve_as(&c, irr) {
        Some(t) => {
            let x1 = v.mul(&t, irr);
            let x2 = x1.add(v);
            vec![x1, x2]
        }
        None => Vec::new(),
    }
}

fn mobius_apply(m: &[F2mElement; 4], v: &F2mElement, irr: &IrreduciblePoly) -> Option<F2mElement> {
    let num = m[0].mul(v, irr).add(&m[1]);
    let den = m[2].mul(v, irr).add(&m[3]);
    Some(num.mul(&den.flt_inverse(irr)?, irr))
}

fn wu_models(curve: &BinaryCurve, seed: u64) -> Result<(Vec<Model>, Value), String> {
    let irr = &curve.irreducible;
    let n = curve.m;
    let one = F2mElement::one(n);
    let zero = F2mElement::zero(n);
    // Wu: d = b^{1/4}; X = d(u+1)/u is the Weierstrass x (shown below).
    let d = curve.b.square_k_times(n - 2, irr);
    assert_eq!(d.square(irr).square(irr), curve.b);
    let d_inv = d.flt_inverse(irr).unwrap();
    // Check the state-model equation on points: u = d/(x + d), v = (y + d²)/x
    // when Tr(a) = 0 (s² + s = a), the model's Weierstrass form being
    // Y² + XY = X³ + d²X.
    let s_wu = solve_as(&curve.a, irr);
    let mut rng = StdRng::seed_from_u64(seed ^ 0x57);
    let mut sample = Vec::new();
    while sample.len() < 64 {
        let x = fe(rng.gen_range(2..(1u64 << n)), n);
        sample.extend(points_with_x(curve, &x));
    }
    let mut wu_eq_ok = 0usize;
    if let Some(s) = &s_wu {
        for p in &sample {
            if let BinaryPoint::Affine { x, y } = p {
                if *x == d {
                    continue;
                }
                let u = d.mul(&x.add(&d).flt_inverse(irr).unwrap(), irr);
                // Y = y + s x + d² on Y² + XY = X³ + d²X; v = Y/X.
                let yy = y.add(&s.mul(x, irr)).add(&d.square(irr));
                let v = yy.mul(&x.flt_inverse(irr).unwrap(), irr);
                let lhs = u.square(irr).add(&u).mul(&v.square(irr).add(&v), irr);
                if lhs == d {
                    wu_eq_ok += 1;
                }
            }
        }
    }
    // Binary Edwards: fit w = x_e + y_e as a Möbius function of x, and
    // failing that of the Vélu abscissa X = x + √b/x of E/⟨T⟩.
    let (d1, d2, s) = binary_edwards_params(curve, seed ^ 0xBE)
        .ok_or_else(|| "no binary Edwards parameters".to_string())?;
    let mut fit_pts = Vec::new();
    for p in &sample {
        if let Some((xe, ye)) = to_edwards(p, &d1, &d2, &s, irr) {
            if let BinaryPoint::Affine { x, .. } = p {
                fit_pts.push((x.clone(), xe.add(&ye)));
            }
        }
    }
    if fit_pts.len() < 8 {
        return Err("binary Edwards map failed its curve-equation check".into());
    }
    // Is w invariant under translation by T = (0, √b)?
    let t_pt = points_with_x(curve, &zero).into_iter().next().unwrap();
    let mut t_invariant = 0usize;
    let mut t_tested = 0usize;
    for p in sample.iter().take(16) {
        let q = point_add(curve, p, &t_pt);
        if let (Some((a1, b1)), Some((a2, b2))) = (
            to_edwards(p, &d1, &d2, &s, irr),
            to_edwards(&q, &d1, &d2, &s, irr),
        ) {
            t_tested += 1;
            if a1.add(&b1) == a2.add(&b2) {
                t_invariant += 1;
            }
        }
    }
    let try_fit = |base: &dyn Fn(&F2mElement) -> Option<F2mElement>| -> Option<[F2mElement; 4]> {
        let pts: Vec<(F2mElement, F2mElement)> = fit_pts
            .iter()
            .filter_map(|(x, w)| Some((base(x)?, w.clone())))
            .collect();
        // w (γv + δ) = αv + β  ⇔  α v + β + γ v w + δ w = 0.
        // All points (the sample holds ±P pairs, so any three may repeat).
        let rows: Vec<Vec<F2mElement>> = pts
            .iter()
            .map(|(v, w)| vec![v.clone(), one.clone(), v.mul(w, irr), w.clone()])
            .collect();
        let (nullity, kv) = kernel_vec(rows, 4, irr);
        if nullity != 1 {
            return None;
        }
        let kv = kv?;
        let k: [F2mElement; 4] = [kv[0].clone(), kv[1].clone(), kv[2].clone(), kv[3].clone()];
        pts.iter()
            .all(|(v, w)| mobius_apply(&k, v, irr).as_ref() == Some(w))
            .then_some(k)
    };
    let (be_quotient, be_from_x) = if let Some(k) = try_fit(&|x| Some(x.clone())) {
        (false, k)
    } else if let Some(k) = try_fit(&|x| velu_x(x, &curve.b, irr)) {
        (true, k)
    } else {
        return Err("binary Edwards w is Möbius in neither x nor X = x + sqrt(b)/x".into());
    };
    // Inverse of [α β γ δ] in char 2 is [δ β γ α].
    let be_to_x = [
        be_from_x[3].clone(),
        be_from_x[1].clone(),
        be_from_x[2].clone(),
        be_from_x[0].clone(),
    ];
    let models = vec![
        Model {
            name: "weierstrass-x",
            quotient: false,
            to_x: [one.clone(), zero.clone(), zero.clone(), one.clone()],
            from_x: [one.clone(), zero.clone(), zero.clone(), one.clone()],
        },
        Model {
            // k = (u + 1)/u = x/d.
            name: "wu-kummer-ratio",
            quotient: false,
            to_x: [d.clone(), zero.clone(), zero.clone(), one.clone()],
            from_x: [d_inv.clone(), zero.clone(), zero.clone(), one.clone()],
        },
        Model {
            // u = d/(x + d), x = d(u + 1)/u.
            name: "wu-u",
            quotient: false,
            to_x: [d.clone(), d.clone(), one.clone(), zero.clone()],
            from_x: [zero.clone(), d.clone(), one.clone(), d.clone()],
        },
        Model {
            name: "binary-edwards-w",
            quotient: be_quotient,
            to_x: be_to_x,
            from_x: be_from_x,
        },
    ];
    let info = json!({
        "wu_d": bits(&d),
        "wu_model_applies": s_wu.is_some(),
        "wu_equation_checked_points": if s_wu.is_some() { wu_eq_ok } else { 0 },
        "wu_sample_points": sample.len(),
        "be_d1": bits(&d1), "be_d2": bits(&d2),
        "be_points_checked": fit_pts.len(),
        "be_w_invariant_under_T": [t_invariant, t_tested],
        "be_w_is_mobius_of": if be_quotient { "X = x + sqrt(b)/x (E/<T>)" } else { "x" },
        "be_w_mobius_in_x": [bits(&models[3].from_x[0]), bits(&models[3].from_x[1]), bits(&models[3].from_x[2]), bits(&models[3].from_x[3])],
    });
    Ok((models, info))
}

fn rand_point_any(curve: &BinaryCurve, rng: &mut StdRng) -> BinaryPoint {
    let n = curve.m;
    loop {
        let x = fe(rng.gen_range(2..(1u64 << n)), n);
        let ps = points_with_x(curve, &x);
        if !ps.is_empty() {
            return ps[rng.gen_range(0..ps.len())].clone();
        }
    }
}

/// Does some sign choice of points with these abscissae sum to `target`?
fn lifts(curve: &BinaryCurve, xs: &[Vec<F2mElement>], target: &BinaryPoint) -> bool {
    fn walk(
        c: &BinaryCurve,
        xs: &[Vec<F2mElement>],
        d: usize,
        acc: &BinaryPoint,
        t: &BinaryPoint,
    ) -> bool {
        if d == xs.len() {
            return acc == t;
        }
        xs[d].iter().any(|x| {
            points_with_x(c, x)
                .iter()
                .any(|p| walk(c, xs, d + 1, &point_add(c, acc, p), t))
        })
    }
    walk(curve, xs, 0, &BinaryPoint::Infinity, target)
}

#[allow(clippy::too_many_arguments)]
fn run_wu_cell(
    slug: &str,
    curve: BinaryCurve,
    ell: usize,
    targets: usize,
    rounds: usize,
    seed: u64,
    max_degree: u32,
    dmax: u32,
    node_budget: usize,
) -> Value {
    let irr = curve.irreducible.clone();
    let n = curve.m;
    let st = FieldStructure::new(n, &irr);
    let (models, info) = match wu_models(&curve, seed) {
        Ok(v) => v,
        Err(e) => return json!({"curve": slug, "error": e}),
    };
    let basis: Vec<F2mElement> = (0..ell as u32)
        .map(|k| F2mElement::from_bit_positions(&[k], n))
        .collect();
    let s3 = fpoly_from_terms(&plain_terms(2).unwrap(), n, Some(&curve.b));
    let bpp = quotient_b(&curve, seed ^ 0xE6, 32);
    let s3q = bpp
        .as_ref()
        .map(|bq| fpoly_from_terms(&plain_terms(2).unwrap(), n, Some(bq)));
    let polys: Vec<FPoly> = models
        .iter()
        .map(|md| {
            let base = if md.quotient {
                s3q.as_ref().expect("quotient b''")
            } else {
                &s3
            };
            mobius_substitute(
                base,
                &[2, 2, 2],
                &[md.to_x.clone(), md.to_x.clone(), md.to_x.clone()],
                &irr,
            )
        })
        .collect();
    // Sanity: each model polynomial vanishes on random relations.
    let mut rng = StdRng::seed_from_u64(seed);
    let mut poly_checks = vec![0usize; models.len()];
    for _ in 0..32 {
        let p1 = rand_point_any(&curve, &mut rng);
        let p2 = rand_point_any(&curve, &mut rng);
        let p3 = point_neg(&point_add(&curve, &p1, &p2));
        let xs: Vec<F2mElement> = [&p1, &p2, &p3]
            .iter()
            .filter_map(|p| match p {
                BinaryPoint::Affine { x, .. } => Some(x.clone()),
                _ => None,
            })
            .collect();
        if xs.len() < 3 {
            continue;
        }
        for (mi, md) in models.iter().enumerate() {
            let cs: Option<Vec<F2mElement>> = xs
                .iter()
                .map(|x| mobius_apply(&md.from_x, &model_base(md, x, &curve.b, &irr)?, &irr))
                .collect();
            if let Some(cs) = cs {
                if fpoly_eval(&polys[mi], &cs, &irr).is_zero() {
                    poly_checks[mi] += 1;
                }
            }
        }
    }
    // Targets: uniform affine points (no subgroup restriction: these models
    // are compared as polynomials and systems, not as an index calculus).
    let mut tgts = Vec::new();
    while tgts.len() < targets {
        let p = rand_point_any(&curve, &mut rng);
        if let BinaryPoint::Affine { x, .. } = &p {
            if models.iter().all(|md| {
                model_base(md, x, &curve.b, &irr)
                    .and_then(|v| mobius_apply(&md.from_x, &v, &irr))
                    .is_some()
            }) {
                tgts.push(p);
            }
        }
    }
    let span: Vec<F2mElement> = (0..(1u64 << ell)).map(|v| fe(v, n)).collect();
    let mut arms_out = Vec::new();
    let mut cells: Vec<Vec<(Vec<F2BoolPoly>, usize, BTreeSet<u64>, usize)>> = Vec::new();
    for t in &tgts {
        let BinaryPoint::Affine { x: x_r, .. } = t else {
            unreachable!()
        };
        let mut row = Vec::new();
        for (mi, md) in models.iter().enumerate() {
            let c_r = mobius_apply(
                &md.from_x,
                &model_base(md, x_r, &curve.b, &irr).unwrap(),
                &irr,
            )
            .unwrap();
            let nv = 2 * ell;
            let vars: Vec<SymElement> = (0..2)
                .map(|i| SymElement::from_subspace_vars(&basis, i * ell, n, nv))
                .collect();
            let eqs = descend(&polys[mi], &vars, &c_r, &st, &irr, nv);
            // Brute force over V × V.
            let mut truth = BTreeSet::new();
            let mut lifting = 0usize;
            for (i1, c1) in span.iter().enumerate() {
                for (i2, c2) in span.iter().enumerate() {
                    if fpoly_eval(&polys[mi], &[c1.clone(), c2.clone(), c_r.clone()], &irr)
                        .is_zero()
                    {
                        truth.insert(i1 as u64 | ((i2 as u64) << ell));
                        let xs: Option<Vec<Vec<F2mElement>>> = [c1, c2]
                            .iter()
                            .map(|c| {
                                Some(base_to_xs(
                                    md,
                                    &mobius_apply(&md.to_x, c, &irr)?,
                                    &curve.b,
                                    &irr,
                                ))
                            })
                            .collect();
                        if let Some(xs) = xs {
                            if lifts(&curve, &xs, t) {
                                lifting += 1;
                            }
                        }
                    }
                }
            }
            row.push((eqs, nv, truth, lifting));
        }
        cells.push(row);
    }
    let mut uptimes = Vec::new();
    let mut samples: Vec<Vec<Vec<SolveSample>>> = vec![vec![Vec::new(); tgts.len()]; models.len()];
    let mut identical: Vec<Vec<bool>> = vec![vec![true; tgts.len()]; models.len()];
    for round in 0..rounds {
        uptimes.push(uptime());
        for (ti, row) in cells.iter().enumerate() {
            for k in 0..models.len() {
                let mi = (k + round) % models.len();
                let (eqs, nv, truth, _) = &row[mi];
                let (roots, s) = timed_solve(eqs, *nv, max_degree, node_budget);
                let set: BTreeSet<u64> = roots.into_iter().collect();
                if set != *truth {
                    identical[mi][ti] = false;
                }
                samples[mi][ti].push(s);
            }
        }
    }
    for (mi, md) in models.iter().enumerate() {
        let ms: Vec<f64> = samples[mi].iter().flatten().map(|s| s.ms).collect();
        let wo: Vec<f64> = samples[mi]
            .iter()
            .flatten()
            .map(|s| s.word_ops as f64)
            .collect();
        arms_out.push(json!({
            "arm": md.name,
            "poly": poly_shape(&polys[mi]),
            "poly_vanishes_on_random_relations": poly_checks[mi],
            "system": system_shape(&cells[0][mi].0, cells[0][mi].1),
            "degrees": cells.iter().map(|r| degrees(&r[mi].0, r[mi].1, dmax)).collect::<Vec<_>>(),
            "roots_per_target": cells.iter().map(|r| r[mi].2.len()).collect::<Vec<_>>(),
            "lifting_roots_per_target": cells.iter().map(|r| r[mi].3).collect::<Vec<_>>(),
            "root_sets_equal_brute_force": identical[mi],
            "samples": samples[mi],
            "median_ms": median(&ms),
            "median_word_ops": median(&wo),
        }));
    }
    json!({
        "curve": slug,
        "n": n, "a": bits(&curve.a), "b": bits(&curve.b), "ell": ell,
        "quotient_b_double_prime": bpp.as_ref().map(bits),
        "subspace": "polynomial basis <1, z, ..., z^(ell-1)>",
        "targets": tgts.len(), "rounds": rounds,
        "engine": engine_label(max_degree),
        "model_info": info,
        "uptime_per_round": uptimes,
        "arms": arms_out,
    })
}

// ════════════════════════════════════════════════════════════════════
// Part C — Frobenius auxiliary variable y = x^{2^k}
// ════════════════════════════════════════════════════════════════════

/// The direct `S_{m+1}` (`b = 1`) with every power `x_i^{2^k}` (`k ≥ 1`) of a
/// summand replaced by a fresh block `y_{i,k}` ranging over the same
/// Frobenius-stable `V`, tied to `x_i` by the `n` linear equations
/// `y_{i,k} = x_i^{2^k}`.  Variables: `x` blocks first (`m·ℓ`), then the
/// `y_{i,k}` blocks.
fn build_frob_aux_system(
    basis: &[F2mElement],
    x_r: &F2mElement,
    m: usize,
    st: &FieldStructure,
    irr: &IrreduciblePoly,
) -> (Vec<F2BoolPoly>, usize) {
    let n = st.n;
    let ell = basis.len();
    let terms = plain_terms(m).unwrap();
    let max_e = terms
        .iter()
        .flat_map(|t| t[..m].iter().copied())
        .max()
        .unwrap();
    let ks: Vec<u32> = (1..).take_while(|k| (1u32 << k) <= max_e).collect();
    let n_vars = m * ell + m * ks.len() * ell;
    let xs: Vec<SymElement> = (0..m)
        .map(|i| SymElement::from_subspace_vars(basis, i * ell, n, n_vars))
        .collect();
    let y_off = |i: usize, kk: usize| m * ell + (i * ks.len() + kk) * ell;
    let ys: Vec<Vec<SymElement>> = (0..m)
        .map(|i| {
            (0..ks.len())
                .map(|kk| SymElement::from_subspace_vars(basis, y_off(i, kk), n, n_vars))
                .collect()
        })
        .collect();
    let mut eqs = Vec::new();
    // Couplings.
    for i in 0..m {
        let mut pw = xs[i].clone();
        for (kk, _) in ks.iter().enumerate() {
            pw = pw.square(st);
            let c = ys[i][kk].add(&pw);
            eqs.extend(c.coords.into_iter().filter(|p| !p.is_zero()));
        }
    }
    // S_{m+1} with x_i^e = Π_{bit j of e} (j = 0 ? x_i : y_{i,j}).
    let mut acc = SymElement::zero(n, n_vars);
    for t in &terms {
        let mut coef = F2mElement::one(n);
        coef = coef.mul(&fpow(x_r, t[m], irr), irr);
        let mut prod: Option<SymElement> = None;
        for i in 0..m {
            let e = t[i];
            for j in 0..32u32 {
                if (e >> j) & 1 == 0 {
                    continue;
                }
                let f = if j == 0 {
                    xs[i].clone()
                } else {
                    ys[i][(j - 1) as usize].clone()
                };
                prod = Some(match prod {
                    None => f,
                    Some(q) => q.mul(&f, st),
                });
            }
        }
        let term = match prod {
            None => SymElement::constant(&coef, n, n_vars),
            Some(q) if bits(&coef) == 1 => q,
            Some(q) => q.mul_constant(&coef, st),
        };
        acc = acc.add(&term);
    }
    eqs.extend(acc.coords.into_iter().filter(|p| !p.is_zero()));
    (eqs, n_vars)
}

#[allow(clippy::too_many_arguments)]
fn run_frob_cell(
    a: u8,
    n: u32,
    m: usize,
    targets: usize,
    rounds: usize,
    seed: u64,
    max_degree: u32,
    dmax: u32,
    node_budget: usize,
) -> Option<Value> {
    let kc = KoblitzCurve::new(a, n)?;
    let irr = kc.curve.irreducible.clone();
    let st = FieldStructure::new(n, &irr);
    let divisor = divisor_for_dimension(n, (n + 1).div_ceil(m as u32))?;
    let fb_x = build_frobenius_factor_base_from_divisor(&kc, &divisor)?;
    let basis = fb_x.subspace_basis.clone();
    let ell = basis.len();
    let mut rng = StdRng::seed_from_u64(seed);
    let g = kc.generator().clone();
    let r = kc
        .subgroup_order
        .to_u64_digits()
        .first()
        .copied()
        .unwrap_or(2)
        .max(2);
    let mut cells = Vec::new();
    while cells.len() < targets {
        let p = kc.mul(&g, &BigUint::from(rng.gen_range(1..r)));
        let BinaryPoint::Affine { x, .. } = &p else {
            continue;
        };
        let (e0, n0) = build_direct_x_system(&basis, x, m, &st)?;
        let (e1, n1) = build_frob_aux_system(&basis, x, m, &st, &irr);
        cells.push([(e0, n0), (e1, n1)]);
    }
    let names = ["x-direct", "x-frob-aux"];
    let mask = if m * ell >= 64 {
        u64::MAX
    } else {
        (1u64 << (m * ell)) - 1
    };
    let mut uptimes = Vec::new();
    let mut samples: Vec<Vec<Vec<SolveSample>>> = vec![vec![Vec::new(); targets]; 2];
    let mut sets: Vec<Vec<Option<BTreeSet<u64>>>> = vec![vec![None; targets]; 2];
    let mut stable = true;
    for round in 0..rounds {
        uptimes.push(uptime());
        for (ti, cell) in cells.iter().enumerate() {
            for k in 0..2 {
                let ai = (k + round) % 2;
                let (eqs, nv) = &cell[ai];
                let (roots, s) = timed_solve(eqs, *nv, max_degree, node_budget);
                let set: BTreeSet<u64> = roots.iter().map(|r| r & mask).collect();
                if sets[ai][ti].is_none() {
                    sets[ai][ti] = Some(set);
                } else if sets[ai][ti].as_ref() != Some(&set) {
                    stable = false;
                }
                samples[ai][ti].push(s);
            }
        }
    }
    let identical: Vec<bool> = (0..targets).map(|ti| sets[0][ti] == sets[1][ti]).collect();
    let arms: Vec<Value> = (0..2)
        .map(|ai| {
            let ms: Vec<f64> = samples[ai].iter().flatten().map(|s| s.ms).collect();
            let wo: Vec<f64> = samples[ai].iter().flatten().map(|s| s.word_ops as f64).collect();
            json!({
                "arm": names[ai],
                "system": system_shape(&cells[0][ai].0, cells[0][ai].1),
                "degrees": cells.iter().map(|c| degrees(&c[ai].0, c[ai].1, dmax)).collect::<Vec<_>>(),
                "samples": samples[ai],
                "median_ms": median(&ms),
                "median_word_ops": median(&wo),
            })
        })
        .collect();
    Some(json!({
        "cell": format!("K{}/F2^{} m={}", a, n, m),
        "a": a, "n": n, "m": m, "ell": ell,
        "subgroup_order": kc.subgroup_order.to_string(),
        "targets": targets, "rounds": rounds,
        "engine": engine_label(max_degree),
        "x_root_sets_identical": identical,
        "roots_stable_across_rounds": stable,
        "uptime_per_round": uptimes,
        "arms": arms,
    }))
}

// ════════════════════════════════════════════════════════════════════
// Part D — convolution (La Scala–Marchesin–Tiwari) splitting of e_i
// ════════════════════════════════════════════════════════════════════

fn anf_to_bool(p: &AnfPoly, offset: u32, n_vars: usize) -> F2BoolPoly {
    let monos: Vec<F2BoolMono> = p
        .monomials()
        .map(|vs| F2BoolMono::from_mask(vs.iter().fold(0u64, |m, &v| m | (1u64 << (v + offset)))))
        .collect();
    F2BoolPoly::from_monos(monos, n_vars)
}

/// Same, for an ANF whose variables live in two spaces (`< split` → `off_lo`,
/// `≥ split` → shifted by `off_hi − split`).
fn anf_to_bool_two(p: &AnfPoly, split: u32, off_lo: u32, off_hi: u32, n_vars: usize) -> F2BoolPoly {
    let monos: Vec<F2BoolMono> = p
        .monomials()
        .map(|vs| {
            F2BoolMono::from_mask(vs.iter().fold(0u64, |m, &v| {
                let idx = if v < split {
                    v + off_lo
                } else {
                    v - split + off_hi
                };
                m | (1u64 << idx)
            }))
        })
        .collect();
    F2BoolPoly::from_monos(monos, n_vars)
}

fn e_len(i: usize, l: usize) -> usize {
    i * l - (i - 1)
}

/// The repo's split: `x` (3l) then `e` (6l − 3); correspondence
/// `e_{i,d} = σ_{i,d}(X)` plus the `n` descended `S₄` equations in `e`.
fn lucas_repo_system(
    n: u32,
    l: u32,
    irr: &IrreduciblePoly,
    x_r: &F2mElement,
) -> (Vec<F2BoolPoly>, usize) {
    let sys = weil_descend_s4(n, l, irr, &F2mElement::one(n), x_r);
    let nx = sys.n_x_vars();
    let ne = sys.n_e_vars();
    let nv = (nx + ne) as usize;
    let mut eqs = Vec::new();
    for (i, row) in sys.correspondence.iter().enumerate() {
        for (d, sigma) in row.iter().enumerate() {
            let mut p = anf_to_bool(sigma, 0, nv);
            p = p.add(&F2BoolPoly::from_monos(
                vec![F2BoolMono::var(nx + sys.e_var(i, d))],
                nv,
            ));
            eqs.push(p);
        }
    }
    for p in &sys.semaev {
        let q = anf_to_bool(p, nx, nv);
        if !q.is_zero() {
            eqs.push(q);
        }
    }
    (eqs, nv)
}

/// The convolution split: `x` (3l), `e` (6l − 3), then `f = X₁X₂` (2l − 1):
/// `f = X₁X₂`, `e₁ = X₁ + X₂ + X₃`, `e₂ = f + (X₁ + X₂)X₃`, `e₃ = f·X₃`, all
/// of degree ≤ 2, plus the same descended `S₄` in `e`.
fn lucas_conv_system(
    n: u32,
    l: u32,
    irr: &IrreduciblePoly,
    x_r: &F2mElement,
) -> (Vec<F2BoolPoly>, usize) {
    let sys = weil_descend_s4(n, l, irr, &F2mElement::one(n), x_r);
    let lu = l as usize;
    let nx = 3 * l;
    let ne = sys.n_e_vars();
    let nf = 2 * l - 1;
    let nv = (nx + ne + nf) as usize;
    let f_off = nx + ne;
    // Work in one ANF variable space: x at 0..3l, f at 3l..3l+nf.
    let xs: Vec<AnfF2m> = (0..3)
        .map(|i| AnfF2m::from_vars(i as u32 * l, lu))
        .collect();
    let f = AnfF2m::from_vars(nx, nf as usize);
    let x12 = xs[0].mul(&xs[1]);
    let s12 = xs[0].xor(&xs[1]);
    let e1 = s12.xor(&xs[2]);
    let e2 = f.xor(&s12.mul(&xs[2]));
    let e3 = f.mul(&xs[2]);
    let to_bool = |p: &AnfPoly| anf_to_bool_two(p, nx, 0, f_off, nv);
    let mut eqs = Vec::new();
    for d in 0..nf as usize {
        let rhs = x12.coeffs.get(d).cloned().unwrap_or_default();
        let p = to_bool(&rhs).add(&F2BoolPoly::from_monos(
            vec![F2BoolMono::var(f_off + d as u32)],
            nv,
        ));
        eqs.push(p);
    }
    for (i, ei) in [e1, e2, e3].iter().enumerate() {
        for d in 0..e_len(i + 1, lu) {
            let rhs = ei.coeffs.get(d).cloned().unwrap_or_default();
            let p = to_bool(&rhs).add(&F2BoolPoly::from_monos(
                vec![F2BoolMono::var(nx + sys.e_var(i, d))],
                nv,
            ));
            eqs.push(p);
        }
    }
    for p in &sys.semaev {
        let q = anf_to_bool(p, nx, nv);
        if !q.is_zero() {
            eqs.push(q);
        }
    }
    (eqs, nv)
}

/// Correspondence-only shapes at `m` summands: direct `σ_i(X₁..X_m)` against
/// the convolution split (`m = 3`: `f = X₁X₂`; `m = 4`: `f = X₁X₂`,
/// `g = X₃X₄`).  Returns `(equations, max degree, distinct monomials)`.
fn correspondence_shapes(m: usize, l: u32) -> Value {
    let lu = l as usize;
    let shape = |polys: &[AnfPoly]| {
        let mut monos = BTreeSet::new();
        let mut deg = 0;
        let mut terms = 0;
        for p in polys {
            deg = deg.max(p.degree());
            for mo in p.monomials() {
                monos.insert(mo.clone());
                terms += 1;
            }
        }
        json!({"equations": polys.len(), "max_degree": deg, "distinct_monomials": monos.len(), "total_terms": terms})
    };
    let xs: Vec<AnfF2m> = (0..m)
        .map(|i| AnfF2m::from_vars(i as u32 * l, lu))
        .collect();
    // Direct elementary symmetric functions by the recurrence e_k ← e_k + e_{k−1}·X.
    let mut es: Vec<AnfF2m> = vec![AnfF2m::from_const(&F2mElement::one(1), 1)];
    for x in &xs {
        let mut next = es.clone();
        next.push(AnfF2m::zero(0));
        for k in (1..next.len()).rev() {
            let add = es[k - 1].mul(x);
            next[k] = if k < es.len() { es[k].xor(&add) } else { add };
        }
        es = next;
    }
    let direct: Vec<AnfPoly> = es[1..].iter().flat_map(|e| e.coeffs.clone()).collect();
    let base = (m as u32) * l;
    let conv: Vec<AnfPoly> = match m {
        2 => direct.clone(),
        3 => {
            let f = AnfF2m::from_vars(base, 2 * lu - 1);
            let x12 = xs[0].mul(&xs[1]);
            let s12 = xs[0].xor(&xs[1]);
            let mut v: Vec<AnfPoly> = x12.xor(&f).coeffs;
            v.extend(s12.xor(&xs[2]).coeffs);
            v.extend(f.xor(&s12.mul(&xs[2])).coeffs);
            v.extend(f.mul(&xs[2]).coeffs);
            v
        }
        4 => {
            let f = AnfF2m::from_vars(base, 2 * lu - 1);
            let g = AnfF2m::from_vars(base + 2 * l - 1, 2 * lu - 1);
            let s12 = xs[0].xor(&xs[1]);
            let s34 = xs[2].xor(&xs[3]);
            let mut v: Vec<AnfPoly> = xs[0].mul(&xs[1]).xor(&f).coeffs;
            v.extend(xs[2].mul(&xs[3]).xor(&g).coeffs);
            v.extend(s12.xor(&s34).coeffs);
            v.extend(f.xor(&g).xor(&s12.mul(&s34)).coeffs);
            v.extend(f.mul(&s34).xor(&g.mul(&s12)).coeffs);
            v.extend(f.mul(&g).coeffs);
            v
        }
        _ => unreachable!(),
    };
    json!({"m": m, "l": l, "direct": shape(&direct), "convolution": shape(&conv)})
}

#[allow(clippy::too_many_arguments)]
fn run_lucas_cell(
    n: u32,
    l: u32,
    targets: usize,
    rounds: usize,
    seed: u64,
    max_degree: u32,
    dmax: u32,
    node_budget: usize,
) -> Option<Value> {
    let irr = crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse(n)?;
    let mut rng = StdRng::seed_from_u64(seed);
    let lu = l as usize;
    let span: Vec<F2mElement> = (0..(1u64 << l)).map(|v| fe(v, n)).collect();
    let mut cells = Vec::new();
    let mut truths = Vec::new();
    while cells.len() < targets {
        let x_r = fe(rng.gen_range(2..(1u64 << n)), n);
        let rep = lucas_repo_system(n, l, &irr, &x_r);
        let conv = lucas_conv_system(n, l, &irr, &x_r);
        let mut truth = BTreeSet::new();
        for (a, x1) in span.iter().enumerate() {
            for (b, x2) in span.iter().enumerate() {
                for (c, x3) in span.iter().enumerate() {
                    let (e1, e2, e3) = elementary_symmetric_3(x1, x2, x3, &irr);
                    if symmetrised_s4_eval(&e1, &e2, &e3, &x_r, &irr).is_zero() {
                        truth.insert(a as u64 | ((b as u64) << l) | ((c as u64) << (2 * l)));
                    }
                }
            }
        }
        truths.push(truth);
        cells.push([rep, conv]);
    }
    let names = ["e-split-repo", "e-split-convolution"];
    let mask = (1u64 << (3 * lu)) - 1;
    let mut uptimes = Vec::new();
    let mut samples: Vec<Vec<Vec<SolveSample>>> = vec![vec![Vec::new(); targets]; 2];
    let mut ok: Vec<Vec<bool>> = vec![vec![true; targets]; 2];
    for round in 0..rounds {
        uptimes.push(uptime());
        for (ti, cell) in cells.iter().enumerate() {
            for k in 0..2 {
                let ai = (k + round) % 2;
                let (eqs, nv) = &cell[ai];
                let (roots, s) = timed_solve(eqs, *nv, max_degree, node_budget);
                let set: BTreeSet<u64> = roots.iter().map(|r| r & mask).collect();
                if set != truths[ti] || set.len() != roots.len() {
                    ok[ai][ti] = false;
                }
                samples[ai][ti].push(s);
            }
        }
    }
    let arms: Vec<Value> = (0..2)
        .map(|ai| {
            let ms: Vec<f64> = samples[ai].iter().flatten().map(|s| s.ms).collect();
            let wo: Vec<f64> = samples[ai].iter().flatten().map(|s| s.word_ops as f64).collect();
            json!({
                "arm": names[ai],
                "system": system_shape(&cells[0][ai].0, cells[0][ai].1),
                "degrees": cells.iter().map(|c| degrees(&c[ai].0, c[ai].1, dmax)).collect::<Vec<_>>(),
                "x_root_sets_equal_brute_force": ok[ai],
                "samples": samples[ai],
                "median_ms": median(&ms),
                "median_word_ops": median(&wo),
            })
        })
        .collect();
    Some(json!({
        "cell": format!("S4 symmetrised, b=1, n={n}, l={l}"),
        "n": n, "l": l, "modulus_low_terms": irr.low_terms,
        "targets": targets, "rounds": rounds,
        "engine": engine_label(max_degree),
        "roots_per_target": truths.iter().map(|t| t.len()).collect::<Vec<_>>(),
        "uptime_per_round": uptimes,
        "arms": arms,
    }))
}

// ════════════════════════════════════════════════════════════════════
// Identity certificate: Res_s(S^sym, s² + s + Σw) = S^{E/⟨T⟩}(1/w)·Πw^d over F₂
// ════════════════════════════════════════════════════════════════════

type ZPoly = BTreeMap<Vec<u32>, i64>;

fn zmul(a: &ZPoly, b: &ZPoly) -> ZPoly {
    let mut r = ZPoly::new();
    for (ea, ca) in a {
        for (eb, cb) in b {
            let e: Vec<u32> = ea.iter().zip(eb).map(|(x, y)| x + y).collect();
            *r.entry(e).or_insert(0) += ca * cb;
        }
    }
    r.retain(|_, c| *c != 0);
    r
}

fn zadd(a: &ZPoly, b: &ZPoly, sb: i64) -> ZPoly {
    let mut r = a.clone();
    for (e, c) in b {
        *r.entry(e.clone()).or_insert(0) += sb * c;
    }
    r.retain(|_, c| *c != 0);
    r
}

fn zexpr(p: &ZPoly, names: &[String]) -> Expr {
    Expr::sum(
        p.iter()
            .map(|(e, &c)| {
                let mut f = vec![Expr::Const(c)];
                for (i, &k) in e.iter().enumerate() {
                    for _ in 0..k {
                        f.push(Expr::var(&names[i]));
                    }
                }
                Expr::product(f)
            })
            .collect(),
    )
}

/// Build the statement for arity `m` (2 or 3): variables `w1..wm, wR`.
/// `G + s·H` is the symmetrised polynomial (`s` occurs linearly), and the
/// integer identity certified is
/// `G² − G·H + H²·Σw = E + 2·Q`, with `E` the coefficient-one lift of
/// `S^{E/⟨T⟩}_{m+1}` in `w` and `Q` the integer quotient computed here —
/// which, reduced mod 2, is the `F₂` identity
/// `Res_s(G + sH, s² + s + Σw) = S^{E/⟨T⟩}_{m+1}(1/w)·Π w^{2^{m−1}}`.
fn eg_identity_statement(m: usize) -> IdentityStatement {
    let k = m + 1;
    let sym = symmetrised_terms(m).unwrap();
    let mut g = ZPoly::new();
    let mut h = ZPoly::new();
    for t in &sym {
        let (w, s) = (t[..k].to_vec(), t[k]);
        match s {
            0 => *g.entry(w).or_insert(0) += 1,
            1 => *h.entry(w).or_insert(0) += 1,
            _ => panic!("s occurs linearly"),
        }
    }
    let mut sigma = ZPoly::new();
    for i in 0..k {
        let mut e = vec![0; k];
        e[i] = 1;
        sigma.insert(e, 1);
    }
    let lhs = zadd(
        &zadd(&zmul(&g, &g), &zmul(&g, &h), -1),
        &zmul(&zmul(&h, &h), &sigma),
        1,
    );
    let d = 1u32 << (m - 1);
    let mut eg = ZPoly::new();
    for t in plain_terms(m).unwrap() {
        *eg.entry(t.iter().map(|&e| d - e).collect()).or_insert(0) += 1;
    }
    let diff = zadd(&lhs, &eg, -1);
    assert!(diff.values().all(|c| c % 2 == 0), "not an identity mod 2");
    let q: ZPoly = diff.iter().map(|(e, c)| (e.clone(), c / 2)).collect();
    let names: Vec<String> = (0..m)
        .map(|i| format!("w{}", i + 1))
        .chain(["wR".to_string()])
        .collect();
    let lhs_expr = {
        let ge = zexpr(&g, &names);
        let he = zexpr(&h, &names);
        let se = zexpr(&sigma, &names);
        Expr::sum(vec![
            Expr::times(ge.clone(), ge.clone()),
            Expr::minus(Expr::times(ge, he.clone())),
            Expr::product(vec![he.clone(), he, se]),
        ])
    };
    let rhs_expr = Expr::plus(
        zexpr(&eg, &names),
        Expr::times(Expr::Const(2), zexpr(&q, &names)),
    );
    let deg = syntactic_degree(&lhs_expr).max(syntactic_degree(&rhs_expr));
    IdentityStatement {
        name: format!(
            "Resultant in s of the FGHR symmetrised S_{} (Koblitz, b = 1) equals the Semaev S_{} of E/<T> in w = 1/X, mod 2",
            m + 1,
            m + 1
        ),
        source: "lit2026 W8 (Drylo, ePrint 2026/1299): examples/lit2026_w8_summation_models.rs eg_identity_statement; \
                 symmetrised terms from koblitz_symmetrised::symmetrised_terms, Semaev terms from plain_terms. \
                 Certified over Z as LHS = E + 2Q (the F_2 claim is its reduction mod 2)."
            .to_string(),
        variables: names,
        degree_bound: deg,
        lhs: lhs_expr,
        rhs: rhs_expr,
    }
}

fn eg_certificate(m: usize) -> IdentityCertificate {
    let cert = issue(
        eg_identity_statement(m),
        0x1299_2026_0000_0000 + m as u64,
        32,
        "issued by examples/lit2026_w8_summation_models.rs (cert); checked by \
         cryptanalysis::identity_certificate::check replaying from the sealed id and seed",
        None,
    )
    .expect("identity holds");
    match check(&cert).expect("well-formed") {
        Verdict::Accept { .. } => {}
        other => panic!("self-check failed: {other:?}"),
    }
    cert
}

// ════════════════════════════════════════════════════════════════════
// CLI
// ════════════════════════════════════════════════════════════════════

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).cloned().unwrap_or_default();
    let val = |f: &str| {
        args.iter()
            .position(|a| a == f)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let num = |f: &str, d: usize| val(f).and_then(|v| v.parse::<usize>().ok()).unwrap_or(d);
    let rounds = num("--rounds", 5);
    let targets = num("--targets", 8);
    let seed = num("--seed", 0x5EED) as u64;
    let budget = num("--node-budget", 200_000);
    let _ = INHERITED.set(val("--engine").as_deref() == Some("inherited"));
    let _ = ARMS.set(val("--arms").map(|a| a.split(',').map(str::to_string).collect()));
    let _ = NO_DEGREES.set(args.iter().any(|a| a == "--no-degrees"));
    let only = val("--only");
    let out_path = val("--out");
    let started = chrono_like_now();
    let result = match cmd.as_str() {
        "eg" => {
            // (a, n, m, max_degree, dmax, x-direct)
            let cells: Vec<(u8, u32, usize, u32, u32, bool)> = vec![
                (0, 9, 2, 4, 5, true),
                (1, 15, 2, 4, 5, true),
                (0, 15, 2, 4, 5, true),
                (1, 17, 2, 4, 4, true),
                (0, 23, 2, 4, 4, true),
                (1, 23, 2, 4, 4, true),
                (0, 9, 3, 4, 6, true),
                (0, 15, 3, 4, 6, true),
                (1, 15, 3, 4, 6, true),
                (1, 17, 3, 4, 5, false),
            ];
            let mut out = Vec::new();
            for (a, n, m, md, dmax, xd) in cells {
                let tag = format!("{a}-{n}-{m}");
                if only
                    .as_ref()
                    .is_some_and(|o| !o.split(',').any(|x| x == tag))
                {
                    continue;
                }
                eprintln!("eg cell K{a}/F2^{n} m={m}");
                match run_eg_cell(a, n, m, targets, rounds, seed, md, dmax, budget, xd) {
                    Some(v) => out.push(v),
                    None => out.push(json!({"cell": tag, "error": "cell could not be built"})),
                }
                checkpoint(&out_path, "eg", &out);
            }
            json!({"part": "eg", "cells": out})
        }
        "wu" => {
            // (slug, a, b, modulus low terms, n, ell)
            let curves: Vec<(&str, u64, u64, Vec<u32>, u32, usize)> = vec![
                ("icv1-f2m15-tm275-2d22ff5d", 0, 0x1, vec![0, 1], 15, 8),
                ("icv1-f2m15-t5-d9fb183e", 0, 0x2493, vec![0, 1], 15, 8),
                ("icv1-f2m17-t101-e6c4b64d", 0, 0x1, vec![0, 3], 17, 9),
                ("icv1-f2m17-t209-0eba1a1e", 0, 0x11a8e, vec![0, 3], 17, 9),
            ];
            let mut out = Vec::new();
            for (slug, a, b, low, n, ell) in curves {
                if only
                    .as_ref()
                    .is_some_and(|o| !o.split(',').any(|x| x == slug))
                {
                    continue;
                }
                eprintln!("wu curve {slug}");
                let irr = IrreduciblePoly {
                    degree: n,
                    low_terms: low,
                };
                let curve = BinaryCurve {
                    m: n,
                    irreducible: irr,
                    a: fe(a, n),
                    b: fe(b, n),
                    generator: BinaryPoint::Infinity,
                    order: BigUint::from(0u32),
                    cofactor: BigUint::from(1u32),
                };
                out.push(run_wu_cell(
                    slug, curve, ell, targets, rounds, seed, 4, 5, budget,
                ));
                checkpoint(&out_path, "wu", &out);
            }
            json!({"part": "wu", "cells": out})
        }
        "frob" => {
            // The aux arm has (1 + K)·m·ℓ variables, K = 1 (m = 2) or 2
            // (m = 3); cells are kept where its Macaulay matrices fit.
            let cells: Vec<(u8, u32, usize, u32, u32)> = vec![
                (0, 9, 2, 3, 4),
                (1, 15, 2, 3, 4),
                (0, 15, 2, 3, 4),
                (0, 9, 3, 6, 6),
            ];
            let mut out = Vec::new();
            for (a, n, m, md, dmax) in cells {
                let tag = format!("{a}-{n}-{m}");
                if only
                    .as_ref()
                    .is_some_and(|o| !o.split(',').any(|x| x == tag))
                {
                    continue;
                }
                eprintln!("frob cell K{a}/F2^{n} m={m}");
                match run_frob_cell(a, n, m, targets, rounds, seed, md, dmax, budget) {
                    Some(v) => out.push(v),
                    None => out.push(json!({"cell": tag, "error": "cell could not be built"})),
                }
                checkpoint(&out_path, "frob", &out);
            }
            json!({"part": "frob", "cells": out})
        }
        "lucas" => {
            let mut out = Vec::new();
            for (n, l) in [(9u32, 3u32), (12, 4), (15, 5)] {
                let tag = format!("{n}-{l}");
                if only
                    .as_ref()
                    .is_some_and(|o| !o.split(',').any(|x| x == tag))
                {
                    continue;
                }
                eprintln!("lucas cell n={n} l={l}");
                match run_lucas_cell(n, l, targets, rounds, seed, 3, 4, budget) {
                    Some(v) => out.push(v),
                    None => out.push(json!({"cell": tag, "error": "no modulus"})),
                }
                checkpoint(&out_path, "lucas", &out);
            }
            let shapes: Vec<Value> = [(2usize, 4u32), (3, 4), (3, 5), (4, 4), (4, 5)]
                .iter()
                .map(|&(m, l)| correspondence_shapes(m, l))
                .collect();
            json!({"part": "lucas", "cells": out, "correspondence_shapes": shapes})
        }
        "compare" => compare(&val("--inputs").expect("--inputs label=path,...")),
        "cert" => {
            let c2 = eg_certificate(2);
            let c3 = eg_certificate(3);
            json!({"part": "cert", "certificates": [c2, c3]})
        }
        _ => {
            eprintln!("usage: lit2026_w8_summation_models (eg|wu|frob|lucas|cert) [--rounds R] [--targets T] [--seed S] [--node-budget B] [--only tags] [--out path]");
            std::process::exit(2);
        }
    };
    let wrapped = json!({
        "schema": "lit2026-w8-run/v1",
        "command": args[1..].join(" "),
        "started_unix": started,
        "env": {
            "KIC_F4_MAX_DEGREE_ONLY": std::env::var("KIC_F4_MAX_DEGREE_ONLY").ok(),
            "KIC_F4_INHERIT": std::env::var("KIC_F4_INHERIT").ok(),
        },
        "git_head": std::process::Command::new("git").args(["rev-parse", "HEAD"]).output()
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string()).unwrap_or_default(),
        "result": result,
    });
    let s = serde_json::to_string_pretty(&wrapped).unwrap();
    match out_path {
        Some(p) => std::fs::write(&p, s).expect("write output"),
        None => println!("{s}"),
    }
}

/// `compare --inputs label=path,...`: read `eg` outputs of single-round runs
/// (one process per binary per round, interleaved by the driver) and report,
/// per cell and per `label/arm`, the per-round total milliseconds over all
/// targets, the deterministic word-op totals, and whether the root-set
/// digests agree with every other run of a comparable arm.  Ratios are
/// computed here, in the harness, and only plotted downstream.
fn compare(inputs: &str) -> Value {
    // (cell, label/arm) -> per-round (total ms, total word ops, digests)
    let mut table: BTreeMap<(String, String), Vec<(f64, u64, Vec<String>, String)>> =
        BTreeMap::new();
    let mut uptimes: Vec<(String, String)> = Vec::new();
    for item in inputs.split(',') {
        let (label, path) = item.split_once('=').expect("label=path");
        let v: Value = serde_json::from_str(&std::fs::read_to_string(path).expect("read input"))
            .expect("json");
        for cell in v["result"]["cells"].as_array().into_iter().flatten() {
            let name = cell["cell"].as_str().unwrap_or("?").to_string();
            for u in cell["uptime_per_round"].as_array().into_iter().flatten() {
                uptimes.push((path.to_string(), u.as_str().unwrap_or("").to_string()));
            }
            for arm in cell["arms"].as_array().into_iter().flatten() {
                let aname = arm["arm"].as_str().unwrap_or("?");
                let samples = arm["samples"].as_array().cloned().unwrap_or_default();
                let rounds = samples
                    .first()
                    .and_then(|t| t.as_array())
                    .map(|t| t.len())
                    .unwrap_or(0);
                for r in 0..rounds {
                    let ms: f64 = samples
                        .iter()
                        .map(|t| t[r]["ms"].as_f64().unwrap_or(0.0))
                        .sum();
                    let wo: u64 = samples
                        .iter()
                        .map(|t| t[r]["word_ops"].as_u64().unwrap_or(0))
                        .sum();
                    let dig: Vec<String> = arm["key_set_digests"]
                        .as_array()
                        .into_iter()
                        .flatten()
                        .map(|d| d.as_str().unwrap_or("").to_string())
                        .collect();
                    table
                        .entry((name.clone(), format!("{label}/{aname}")))
                        .or_default()
                        .push((ms, wo, dig, path.to_string()));
                }
            }
        }
    }
    let mut cells: BTreeMap<String, Vec<Value>> = BTreeMap::new();
    let mut by_cell: BTreeMap<String, BTreeMap<String, Vec<(f64, u64, Vec<String>, String)>>> =
        BTreeMap::new();
    for ((cell, key), rows) in table {
        by_cell.entry(cell).or_default().insert(key, rows);
    }
    for (cell, arms) in &by_cell {
        let keys: Vec<&String> = arms.keys().collect();
        let mut out = Vec::new();
        for k in &keys {
            let rows = &arms[*k];
            let ms: Vec<f64> = rows.iter().map(|r| r.0).collect();
            let wo: BTreeSet<u64> = rows.iter().map(|r| r.1).collect();
            out.push(json!({
                "arm": k,
                "rounds": rows.len(),
                "total_ms_per_round": ms,
                "median_total_ms": median(&ms),
                "word_ops_totals_seen": wo,
                "digests": rows.first().map(|r| r.2.clone()).unwrap_or_default(),
                "digests_stable_across_rounds": rows.iter().all(|r| r.2 == rows[0].2),
            }));
        }
        // Pairwise ratios a/b of per-round totals (round i of a against round i of b).
        let mut ratios = Vec::new();
        for a in &keys {
            for b in &keys {
                if a == b {
                    continue;
                }
                let ra = &arms[*a];
                let rb = &arms[*b];
                let n = ra.len().min(rb.len());
                let per: Vec<f64> = (0..n).map(|i| ra[i].0 / rb[i].0).collect();
                let wo_a = ra.first().map(|r| r.1).unwrap_or(0) as f64;
                let wo_b = rb.first().map(|r| r.1).unwrap_or(0) as f64;
                let same_digests = ra.first().map(|r| &r.2) == rb.first().map(|r| &r.2);
                ratios.push(json!({
                    "numerator": a, "denominator": b,
                    "median_ratio_of_median_ms": median(&ra.iter().map(|r| r.0).collect::<Vec<_>>()) / median(&rb.iter().map(|r| r.0).collect::<Vec<_>>()),
                    "per_round_ratio": per,
                    "per_round_ratio_min": per.iter().cloned().fold(f64::INFINITY, f64::min),
                    "per_round_ratio_max": per.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                    "word_ops_ratio": if wo_b > 0.0 { wo_a / wo_b } else { f64::NAN },
                    "root_digests_identical": same_digests,
                }));
            }
        }
        cells.insert(cell.clone(), vec![json!({"arms": out, "ratios": ratios})]);
    }
    json!({"part": "compare", "inputs": inputs, "cells": cells, "uptimes": uptimes.iter().map(|(p, u)| json!([p, u])).collect::<Vec<_>>()})
}

/// Write the cells finished so far, so a long part leaves evidence if stopped.
fn checkpoint(out_path: &Option<String>, part: &str, cells: &[Value]) {
    if let Some(p) = out_path {
        let v =
            json!({"schema": "lit2026-w8-run/v1", "partial": true, "part": part, "cells": cells});
        let _ = std::fs::write(
            format!("{p}.partial"),
            serde_json::to_string_pretty(&v).unwrap(),
        );
    }
}

fn chrono_like_now() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

// ════════════════════════════════════════════════════════════════════
// Tests
// ════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    /// The quotient by `T` of a Koblitz curve is again `b'' = 1` in the Vélu
    /// abscissa `X = x + 1/x`, and `X = 1/w` with `w = u² + u`.
    #[test]
    fn koblitz_quotient_is_b_one_and_w_is_inverse_velu_x() {
        for (a, n) in [(0u8, 9u32), (1, 9), (1, 15)] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let bpp = quotient_b(&kc.curve, 7, 40).expect("consistent b''");
            assert_eq!(bits(&bpp), 1, "K{a}/F2^{n}");
            let irr = &kc.curve.irreducible;
            for xv in 2..200u64 {
                let x = fe(xv, n);
                let u = u_of_x(&x, &kc).unwrap();
                let w = artin_schreier(&u, &kc);
                let big_x = velu_x(&x, &kc.curve.b, irr).unwrap();
                if big_x.is_zero() {
                    continue;
                }
                assert_eq!(w.mul(&big_x, irr), F2mElement::one(n));
            }
        }
    }

    /// `Res_s(sym, s² + s + Σw)` equals the `E/⟨T⟩` Semaev polynomial in `w`
    /// mod 2, for `m = 2, 3`, and the certificate checks.
    #[test]
    fn eg_resultant_identity_certificates_check() {
        for m in [2usize, 3] {
            let c = eg_certificate(m);
            assert!(matches!(check(&c).unwrap(), Verdict::Accept { .. }));
        }
    }

    /// The committed certificates in the study directory still check and
    /// are what this harness issues.
    #[test]
    fn committed_eg_certificates_still_check() {
        let dir = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/research/lit2026_w8_summation_models_20261011/identity-certificates"
        );
        for (m, file) in [
            (2usize, "eg-resultant-s3.json"),
            (3, "eg-resultant-s4.json"),
        ] {
            let text =
                std::fs::read_to_string(format!("{dir}/{file}")).expect("committed certificate");
            let c: IdentityCertificate = serde_json::from_str(&text).unwrap();
            assert!(matches!(check(&c).unwrap(), Verdict::Accept { .. }));
            assert_eq!(c, eg_certificate(m));
        }
    }

    /// The Möbius-substituted `E/⟨T⟩` polynomial in `w` is the
    /// exponent-reversed Semaev term list.
    #[test]
    fn eg_poly_is_reversed_semaev() {
        let kc = KoblitzCurve::new(1, 9).unwrap();
        for m in [2usize, 3] {
            let p = eg_poly(m, 9, &kc.curve.irreducible);
            let d = 1u32 << (m - 1);
            let want: BTreeSet<Vec<u32>> = plain_terms(m)
                .unwrap()
                .iter()
                .map(|t| t.iter().map(|&e| d - e).collect())
                .collect();
            let got: BTreeSet<Vec<u32>> = p.keys().cloned().collect();
            assert_eq!(got, want);
            assert!(p.values().all(|c| bits(c) == 1));
        }
    }

    /// On a small cell the `eg`, `sym` and `u-plain` arms have identical
    /// `w`-tuple root sets and `sym`'s verdicts match enumeration.
    #[test]
    fn eg_cell_root_sets_agree_small() {
        let v = run_eg_cell(0, 9, 2, 4, 1, 3, 3, 4, 100_000, false).unwrap();
        for arm in v["arms"].as_array().unwrap() {
            for ok in arm["key_sets_identical_to_reference"].as_array().unwrap() {
                assert!(ok.as_bool().unwrap(), "{}", arm["arm"]);
            }
        }
        for ok in v["sym_verdict_matches_enumeration"].as_array().unwrap() {
            assert!(ok.as_bool().unwrap());
        }
    }

    /// Frobenius auxiliary blocks leave the projected root set unchanged.
    #[test]
    fn frob_aux_root_sets_agree_small() {
        let v = run_frob_cell(0, 9, 2, 3, 1, 5, 4, 4, 100_000).unwrap();
        for ok in v["x_root_sets_identical"].as_array().unwrap() {
            assert!(ok.as_bool().unwrap());
        }
    }

    /// Both `e`-splittings reproduce the brute-force root set at `l = 3`.
    #[test]
    fn lucas_splits_match_brute_force_small() {
        let v = run_lucas_cell(9, 3, 2, 1, 11, 3, 3, 100_000).unwrap();
        for arm in v["arms"].as_array().unwrap() {
            for ok in arm["x_root_sets_equal_brute_force"].as_array().unwrap() {
                assert!(ok.as_bool().unwrap(), "{}", arm["arm"]);
            }
        }
    }

    /// Every model's `S₃` vanishes on random relations and its descended
    /// system has exactly the brute-force root set.
    #[test]
    fn wu_models_small() {
        let n = 9;
        let kc = KoblitzCurve::new(0, n).unwrap();
        let v = run_wu_cell(
            "icv1-f2m9-t5-81e744be",
            kc.curve.clone(),
            5,
            2,
            1,
            9,
            4,
            3,
            100_000,
        );
        assert!(v.get("error").is_none(), "{v}");
        let info = &v["model_info"];
        // Wu's state-model equation holds at every sampled point (Tr(a) = 0).
        assert_eq!(info["wu_equation_checked_points"], info["wu_sample_points"]);
        // Binary Edwards w is a Dryło compression: Möbius in E/<T>'s abscissa,
        // invariant under T.
        assert_eq!(info["be_w_is_mobius_of"], "X = x + sqrt(b)/x (E/<T>)");
        assert_eq!(
            info["be_w_invariant_under_T"][0],
            info["be_w_invariant_under_T"][1]
        );
        for arm in v["arms"].as_array().unwrap() {
            assert_eq!(
                arm["poly_vanishes_on_random_relations"].as_u64().unwrap() > 0,
                true
            );
            for ok in arm["root_sets_equal_brute_force"].as_array().unwrap() {
                assert!(ok.as_bool().unwrap(), "{}", arm["arm"]);
            }
        }
    }
}
