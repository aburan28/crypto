//! Factor-base closure sweep — EXP-H of the FFD programme.
//!
//! **Question.** Among `F₂`-subspaces `V ⊂ F_{2ⁿ}` of one fixed dimension
//! `l`, which *shape* gives the lowest refutation degree `D*` for the `m = 2`
//! Semaev system, and is the ordering explained by how fast `V` grows under
//! multiplication (`dim V·V`, `dim V³`) and under Frobenius (`dim V ∩ σV`)?
//!
//! The FFD programme measured three shapes (EXP-E/G): the subfield
//! (`V·V = V = σV`, mean `D* ≈ 2.0`), the coordinate subspace
//! `⟨1, z, …, z^{l−1}⟩` (`≈ 2.5`) and random (`≈ 4.0`), and found the
//! early Macaulay defect predicts `D*` (P3-alg). This sweep adds the shapes
//! that sit between those corners and that still exist at prime `n`, where
//! no subfield and no Frobenius-stable subspace exists:
//!
//! - `gp`: geometric progression `a·⟨1, g, …, g^{l−1}⟩`, random ratio `g`.
//!   By the linear Kneser theorem (Hou–Leung–Xiang) every subspace at prime
//!   `n` has `dim V·V ≥ 2l − 1`, and by the linear Vosper theorem
//!   (Bachoc–Serra–Zémor) the geometric progressions are the only ones that
//!   attain it. The coordinate subspace is `g = z`.
//! - `gpsym`: symmetric progression `⟨g^{−k}, …, g^{k}⟩` (inversion-closed).
//! - `gpord`: progression whose ratio has small multiplicative order.
//! - `fp`: Frobenius progression `⟨α, α², α⁴, …, α^{2^{l−1}}⟩`, the subspace
//!   with `dim V ∩ σV = l − 1` (almost Frobenius-stable); `fp2` is the
//!   stride-2 variant `⟨α^{4^i}⟩`.
//! - `mix`: product shape `⟨g^i · α^{2^j}⟩`, `a·b = l`, between `gp` and `fp`.
//!
//! Every cell also records the decomposability rate (the index-calculus
//! yield) of the shape, for a random curve constant and for the Koblitz
//! constant `b = 1`, because a shape that lowers `D*` by becoming a subgroup
//! has bought nothing.
//!
//! **Trace class.** The identity `Tr(S₃(X₁,X₂,x₃)/x₃²) = Tr(X₁) + Tr(X₂) +
//! Tr(√b/x₃)` holds for all `X₁, X₂` (the Kosters–Yeo trace equation, with
//! the certificate `c = 1/x₃²`). So when `V ⊂ ker Tr` every target with
//! `Tr(√b/x₃) = 1` is refuted at degree 2 for free: it lies in the wrong
//! `2E`-coset. Those targets say nothing about the algebra of `V`. Each
//! target is therefore tagged with its class `t = Tr(√b/x₃)`, every cell
//! records whether `V ⊂ ker Tr`, and `D*` is reported per class; the
//! shape comparison uses class-0 targets. The `tz*` shapes are trace-zero
//! versions of the plain shapes, built to the same dimension.
//!
//! Stage diagnostic: `D*` at `m = 2` on toy fields. No `S`, no speedup.
//!
//! Run: `cargo run --release --example factor_base_closure_sweep [seed] [scale]`

use crypto_lib::binary_ecc::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::descent_algebraic::{early_defect, rank_profile};
use crypto_lib::cryptanalysis::descent_expansion::enumerate_irreducibles;
use crypto_lib::cryptanalysis::descent_lowgamma::{
    descend_on_subspace, measure_on_subspace, BasisFamily, FactorSubspace,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

// ── F₂-linear algebra on field elements ─────────────────────────────────

fn to_u128(e: &F2mElement, n: u32) -> u128 {
    let raw = e.raw_bits();
    let mut v = 0u128;
    for j in 0..n as usize {
        if (raw.get(j / 64).copied().unwrap_or(0) >> (j % 64)) & 1 == 1 {
            v |= 1u128 << j;
        }
    }
    v
}

fn rank_u128(mut rows: Vec<u128>) -> usize {
    let mut r = 0usize;
    for bit in 0..128u32 {
        if r == rows.len() {
            break;
        }
        let Some(p) = (r..rows.len()).find(|&i| (rows[i] >> bit) & 1 == 1) else {
            continue;
        };
        rows.swap(r, p);
        let pivot = rows[r];
        for (i, row) in rows.iter_mut().enumerate() {
            if i != r && (*row >> bit) & 1 == 1 {
                *row ^= pivot;
            }
        }
        r += 1;
    }
    r
}

fn span_dim(elems: &[F2mElement], n: u32) -> usize {
    rank_u128(elems.iter().map(|e| to_u128(e, n)).collect())
}

fn products(a: &[F2mElement], b: &[F2mElement], irr: &IrreduciblePoly) -> Vec<F2mElement> {
    let mut out = Vec::with_capacity(a.len() * b.len());
    for x in a {
        for y in b {
            out.push(x.mul(y, irr));
        }
    }
    out
}

fn frobenius(v: &[F2mElement], irr: &IrreduciblePoly) -> Vec<F2mElement> {
    v.iter().map(|e| e.square(irr)).collect()
}

fn pow(g: &F2mElement, mut e: u64, irr: &IrreduciblePoly) -> F2mElement {
    let mut acc = F2mElement::one(g.m_value());
    let mut base = g.clone();
    while e > 0 {
        if e & 1 == 1 {
            acc = acc.mul(&base, irr);
        }
        base = base.square(irr);
        e >>= 1;
    }
    acc
}

/// Absolute trace `Tr(e) = Σ_i e^{2^i} ∈ F₂`.
fn trace(e: &F2mElement, irr: &IrreduciblePoly) -> bool {
    let n = e.m_value();
    let mut acc = F2mElement::zero(n);
    let mut t = e.clone();
    for _ in 0..n {
        acc = acc.add(&t);
        t = t.square(irr);
    }
    acc == F2mElement::one(n)
}

/// A random nonzero element `a` with `Tr(a·e) = 0` for every `e` in
/// `elems` (the trace-annihilator of their span), or `None` if it is `{0}`.
fn random_trace_annihilator(
    elems: &[F2mElement],
    n: u32,
    irr: &IrreduciblePoly,
    rng: &mut StdRng,
) -> Option<F2mElement> {
    // Row i: the linear functional a ↦ Tr(a·e_i) in the coordinates of a.
    let rows: Vec<u128> = elems
        .iter()
        .map(|e| {
            let mut r = 0u128;
            for j in 0..n {
                let zj = F2mElement::from_bit_positions(&[j], n);
                if trace(&zj.mul(e, irr), irr) {
                    r |= 1u128 << j;
                }
            }
            r
        })
        .collect();
    // Kernel of the row space by elimination: reduce rows to echelon form,
    // then sample a random vector orthogonal (over F₂, dot product) to all.
    let mut basis: Vec<u128> = Vec::new();
    for r in rows {
        let mut v = r;
        for b in &basis {
            let pivot = 1u128 << (127 - b.leading_zeros());
            if v & pivot != 0 {
                v ^= b;
            }
        }
        if v != 0 {
            basis.push(v);
            basis.sort_by(|a, b| b.cmp(a));
        }
    }
    for _ in 0..64 {
        let mut a: u128 = 0;
        for j in 0..n {
            if rng.gen::<bool>() {
                a |= 1u128 << j;
            }
        }
        // Project onto the orthogonal complement: for each basis row with
        // pivot p, if <a, row> = 1 flip the pivot bit (rows are echelon, so
        // flipping a pivot bit only affects that row's dot product once
        // processed from highest pivot down with fully reduced rows).
        let mut rows_sorted = basis.clone();
        rows_sorted.sort_by(|x, y| y.cmp(x));
        // Full reduction so each pivot occurs in exactly one row.
        for i in 0..rows_sorted.len() {
            let pivot = 1u128 << (127 - rows_sorted[i].leading_zeros());
            for k in 0..rows_sorted.len() {
                if k != i && rows_sorted[k] & pivot != 0 {
                    rows_sorted[k] ^= rows_sorted[i];
                }
            }
        }
        for r in &rows_sorted {
            if (a & r).count_ones() % 2 == 1 {
                let pivot = 1u128 << (127 - r.leading_zeros());
                a ^= pivot;
            }
        }
        if a == 0 {
            continue;
        }
        let bits: Vec<u32> = (0..n).filter(|&j| (a >> j) & 1 == 1).collect();
        let e = F2mElement::from_bit_positions(&bits, n);
        if elems.iter().all(|x| !trace(&x.mul(&e, irr), irr)) {
            return Some(e);
        }
    }
    None
}

fn rand_nz(rng: &mut StdRng, n: u32) -> F2mElement {
    loop {
        let bits: Vec<u32> = (0..n).filter(|_| rng.gen::<bool>()).collect();
        let e = F2mElement::from_bit_positions(&bits, n);
        if !e.is_zero() {
            return e;
        }
    }
}

fn rand_not_01(rng: &mut StdRng, n: u32) -> F2mElement {
    loop {
        let e = rand_nz(rng, n);
        if e != F2mElement::one(n) {
            return e;
        }
    }
}

fn small_prime_factors(mut x: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut p = 2u64;
    while p * p <= x {
        if x.is_multiple_of(p) {
            out.push(p);
            while x.is_multiple_of(p) {
                x /= p;
            }
        }
        p += 1;
    }
    if x > 1 {
        out.push(x);
    }
    out
}

fn irr_string(irr: &IrreduciblePoly) -> String {
    let mut terms: Vec<String> = vec![format!("z^{}", irr.degree)];
    for &k in irr.low_terms.iter().rev() {
        terms.push(match k {
            0 => "1".to_string(),
            1 => "z".to_string(),
            k => format!("z^{k}"),
        });
    }
    terms.join("+")
}

// ── Shapes ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Shape {
    Coord,
    Subfield,
    Random,
    Gp,
    GpScaled,
    GpSym,
    GpOrd,
    Fp,
    Fp2,
    Mix,
    TzCoord,
    TzGp,
    TzFp,
    TzRandom,
}

impl Shape {
    fn name(self) -> &'static str {
        match self {
            Shape::Coord => "coord",
            Shape::Subfield => "subfield",
            Shape::Random => "random",
            Shape::Gp => "gp",
            Shape::GpScaled => "gpscaled",
            Shape::GpSym => "gpsym",
            Shape::GpOrd => "gpord",
            Shape::Fp => "fp",
            Shape::Fp2 => "fp2",
            Shape::Mix => "mix",
            Shape::TzCoord => "tzcoord",
            Shape::TzGp => "tzgp",
            Shape::TzFp => "tzfp",
            Shape::TzRandom => "tzrandom",
        }
    }
    /// How many independently drawn instances of the shape to run per point.
    fn reseeds(self) -> u32 {
        match self {
            Shape::Coord | Shape::Subfield | Shape::TzCoord => 1,
            Shape::GpOrd => 1,
            _ => 3,
        }
    }
}

/// One draw of a shape as an explicit basis. `None` if the shape does not
/// exist at `(n, l)` or 64 draws failed to reach full rank.
fn draw_shape(
    shape: Shape,
    n: u32,
    l: u32,
    irr: &IrreduciblePoly,
    rng: &mut StdRng,
) -> Option<(Vec<F2mElement>, String)> {
    let lu = l as usize;
    for _attempt in 0..64 {
        let (cols, detail): (Vec<F2mElement>, String) = match shape {
            Shape::Coord => (
                (0..l)
                    .map(|p| F2mElement::from_bit_positions(&[p], n))
                    .collect(),
                "g=z".to_string(),
            ),
            Shape::Subfield => {
                let v = FactorSubspace::build(BasisFamily::Subfield, n, l, irr, 0)?;
                (v.cols, format!("F_{{2^{l}}}"))
            }
            Shape::Random => {
                let seed: u64 = rng.gen();
                let v = FactorSubspace::build(BasisFamily::Random, n, l, irr, seed)?;
                (v.cols, format!("seed={seed:#x}"))
            }
            Shape::Gp => {
                let g = rand_not_01(rng, n);
                let mut cols = Vec::with_capacity(lu);
                let mut acc = F2mElement::one(n);
                for _ in 0..l {
                    cols.push(acc.clone());
                    acc = acc.mul(&g, irr);
                }
                (cols, format!("g={:#x}", to_u128(&g, n)))
            }
            Shape::GpScaled => {
                let a = rand_not_01(rng, n);
                let z = F2mElement::z(n);
                let mut cols = Vec::with_capacity(lu);
                let mut acc = a.clone();
                for _ in 0..l {
                    cols.push(acc.clone());
                    acc = acc.mul(&z, irr);
                }
                (cols, format!("a={:#x} g=z", to_u128(&a, n)))
            }
            Shape::GpSym => {
                let g = rand_not_01(rng, n);
                let ginv = g.flt_inverse(irr)?;
                // exponents −k, …, k (l odd) or −(k−1), …, k (l even)
                let k = l as i64 / 2;
                let lo = if l % 2 == 1 { -k } else { -(k - 1) };
                let mut cols = Vec::with_capacity(lu);
                for e in lo..=k {
                    let v = if e >= 0 {
                        pow(&g, e as u64, irr)
                    } else {
                        pow(&ginv, (-e) as u64, irr)
                    };
                    cols.push(v);
                }
                (cols, format!("g={:#x} sym", to_u128(&g, n)))
            }
            Shape::GpOrd => {
                // Ratio of the smallest prime order d | 2ⁿ−1 with d ≥ 2l−1,
                // so the progression and its square do not wrap.
                let order = (1u64 << n) - 1;
                let d = *small_prime_factors(order)
                    .iter()
                    .find(|&&d| d >= (2 * l as u64).saturating_sub(1).max(2))?;
                let g = rand_not_01(rng, n);
                let zeta = pow(&g, order / d, irr);
                if zeta == F2mElement::one(n) {
                    continue;
                }
                let mut cols = Vec::with_capacity(lu);
                let mut acc = F2mElement::one(n);
                for _ in 0..l {
                    cols.push(acc.clone());
                    acc = acc.mul(&zeta, irr);
                }
                (cols, format!("ord(g)={d}"))
            }
            Shape::Fp | Shape::Fp2 => {
                let stride = if shape == Shape::Fp { 1 } else { 2 };
                let alpha = rand_not_01(rng, n);
                let cols = (0..l)
                    .map(|i| alpha.square_k_times(stride * i, irr))
                    .collect();
                (
                    cols,
                    format!("alpha={:#x} stride={stride}", to_u128(&alpha, n)),
                )
            }
            Shape::Mix => {
                // a·b = l with 1 < a ≤ b; pick the most balanced split.
                let a = (2..=l)
                    .rev()
                    .filter(|a| l.is_multiple_of(*a) && a * a <= l)
                    .max()?;
                let b = l / a;
                let g = rand_not_01(rng, n);
                let alpha = rand_not_01(rng, n);
                let mut cols = Vec::with_capacity(lu);
                for i in 0..a {
                    let gi = pow(&g, i as u64, irr);
                    for j in 0..b {
                        cols.push(gi.mul(&alpha.square_k_times(j, irr), irr));
                    }
                }
                (cols, format!("gp{a}×fp{b}"))
            }
            Shape::TzCoord => {
                // Smallest window ⟨z^s, …, z^{s+l−1}⟩ inside ker Tr.
                let zpow: Vec<F2mElement> = (0..n)
                    .map(|j| F2mElement::from_bit_positions(&[j], n))
                    .collect();
                let tr: Vec<bool> = zpow.iter().map(|e| trace(e, irr)).collect();
                let s = (0..=(n - l)).find(|&s| (0..l).all(|i| !tr[(s + i) as usize]))?;
                (
                    (0..l).map(|i| zpow[(s + i) as usize].clone()).collect(),
                    format!("window z^{s}..z^{}", s + l - 1),
                )
            }
            Shape::TzGp => {
                let g = rand_not_01(rng, n);
                let mut powers = Vec::with_capacity(lu);
                let mut acc = F2mElement::one(n);
                for _ in 0..l {
                    powers.push(acc.clone());
                    acc = acc.mul(&g, irr);
                }
                let Some(a) = random_trace_annihilator(&powers, n, irr, rng) else {
                    continue;
                };
                (
                    powers.iter().map(|p| p.mul(&a, irr)).collect(),
                    format!("a·gp g={:#x}", to_u128(&g, n)),
                )
            }
            Shape::TzFp => {
                let alpha = loop {
                    let a = rand_not_01(rng, n);
                    if !trace(&a, irr) {
                        break a;
                    }
                };
                (
                    (0..l).map(|i| alpha.square_k_times(i, irr)).collect(),
                    format!("alpha={:#x} Tr=0", to_u128(&alpha, n)),
                )
            }
            Shape::TzRandom => {
                let mut cols = Vec::with_capacity(lu);
                while cols.len() < lu {
                    let e = rand_nz(rng, n);
                    if !trace(&e, irr) {
                        cols.push(e);
                    }
                }
                (cols, "random ⊂ ker Tr".to_string())
            }
        };
        if cols.len() == lu && span_dim(&cols, n) == lu {
            return Some((cols, detail));
        }
        if matches!(shape, Shape::Coord | Shape::Subfield | Shape::TzCoord) {
            return None;
        }
    }
    None
}

// ── Decomposability (yield) ─────────────────────────────────────────────

fn subspace_elements(v: &FactorSubspace) -> Vec<F2mElement> {
    (0..(1u32 << v.n_sub)).map(|mask| v.element(mask)).collect()
}

/// `∃ (x₁, x₂) ∈ V²` with `S₃(x₁, x₂, x₃) = 0`, by enumeration. This is
/// exactly "the descended system is satisfiable", including `x₁ = x₂`.
fn decomposable(
    elems: &[F2mElement],
    x3: &F2mElement,
    b: &F2mElement,
    irr: &IrreduciblePoly,
) -> bool {
    let x3_sq = x3.square(irr);
    let sq_x3sq: Vec<F2mElement> = elems
        .iter()
        .map(|e| e.square(irr).mul(&x3_sq, irr))
        .collect();
    for (i, x1) in elems.iter().enumerate() {
        for (j, x2) in elems.iter().enumerate().skip(i) {
            let p = x1.mul(x2, irr);
            let s3 = sq_x3sq[i]
                .add(&sq_x3sq[j])
                .add(&p.mul(x3, irr))
                .add(&p.square(irr))
                .add(b);
            if s3.is_zero() {
                return true;
            }
        }
    }
    false
}

// ── One cell ────────────────────────────────────────────────────────────

struct Cell {
    n: u32,
    l: u32,
    shape: Shape,
    tag: String,
    detail: String,
    irr: String,
    irr_weight: usize,
    dim_vv: usize,
    dim_vvv: usize,
    dim_v_sigv: usize,
    dim_v_cap_sigv: usize,
    yield_rand_b: f64,
    yield_b1: f64,
    yield_draws: u32,
    trace_zero: bool,
    early_defect: f64,
    /// Class-0 targets (`Tr(√b/x₃) = 0`): the ones a trace-zero base cannot
    /// refute for free.
    dstar0_mean: Option<f64>,
    dstar0_hist: Vec<(u32, u32)>,
    n0_measured: u32,
    censored0: u32,
    /// Class-1 targets (`Tr(√b/x₃) = 1`).
    dstar1_mean: Option<f64>,
    dstar1_hist: Vec<(u32, u32)>,
    n1_measured: u32,
    censored1: u32,
    skipped_decomposable: u32,
    seconds: f64,
}

fn hist_of(ds: &[u32]) -> Vec<(u32, u32)> {
    let mut hist: Vec<(u32, u32)> = Vec::new();
    for &d in ds {
        match hist.iter_mut().find(|(k, _)| *k == d) {
            Some((_, c)) => *c += 1,
            None => hist.push((d, 1)),
        }
    }
    hist.sort();
    hist
}

fn mean_of(ds: &[u32]) -> Option<f64> {
    if ds.is_empty() {
        None
    } else {
        Some(ds.iter().map(|&d| d as f64).sum::<f64>() / ds.len() as f64)
    }
}

#[allow(clippy::too_many_arguments)]
fn run_cell(
    n: u32,
    l: u32,
    shape: Shape,
    tag: String,
    cols: Vec<F2mElement>,
    detail: String,
    irr: &IrreduciblePoly,
    targets: u32,
    d_cap: u32,
    yield_draws: u32,
    rng: &mut StdRng,
) -> Cell {
    let t0 = Instant::now();
    let irr_name = irr_string(irr);
    let irr_weight = irr.low_terms.len() + 1;
    let v = FactorSubspace {
        n,
        n_sub: l,
        cols: cols.clone(),
        family: BasisFamily::Random,
    };
    // Growth profile.
    let vv = products(&cols, &cols, irr);
    let dim_vv = span_dim(&vv, n);
    let vvv = products(&vv, &cols, irr);
    let dim_vvv = span_dim(&vvv, n);
    let sig = frobenius(&cols, irr);
    let dim_v_sigv = span_dim(&products(&cols, &sig, irr), n);
    let mut union = cols.clone();
    union.extend(sig.iter().cloned());
    let dim_v_plus_sigv = span_dim(&union, n);
    let dim_v_cap_sigv = 2 * l as usize - dim_v_plus_sigv;

    // Yield, random b and Koblitz b = 1.
    let elems = subspace_elements(&v);
    let one = F2mElement::one(n);
    let mut hits_rand = 0u32;
    let mut hits_b1 = 0u32;
    for _ in 0..yield_draws {
        let b = rand_nz(rng, n);
        let x3 = rand_nz(rng, n);
        if decomposable(&elems, &x3, &b, irr) {
            hits_rand += 1;
        }
        let x3 = rand_nz(rng, n);
        if decomposable(&elems, &x3, &one, irr) {
            hits_b1 += 1;
        }
    }

    // D* over non-decomposable targets only (as pc_degree_avg does), tagged
    // by trace class t = Tr(√b/x₃) = Tr(b/x₃²). Collect `targets` class-0
    // instances; class-1 instances are measured as they come, up to `targets`.
    let trace_zero = cols.iter().all(|c| !trace(c, irr));
    let prof_dmax = 3.min(2 * l);
    let scan_dmax = d_cap.min(2 * l + 2);
    let mut de = Vec::new();
    let mut ds0: Vec<u32> = Vec::new();
    let mut ds1: Vec<u32> = Vec::new();
    let (mut censored0, mut censored1) = (0u32, 0u32);
    let (mut n0, mut n1) = (0u32, 0u32);
    let mut skipped = 0u32;
    let max_attempts = targets.saturating_mul(40).max(64);
    let mut attempts = 0u32;
    while n0 < targets && attempts < max_attempts {
        attempts += 1;
        let b = rand_nz(rng, n);
        let x3 = rand_nz(rng, n);
        let inv = x3.flt_inverse(irr).expect("x3 ≠ 0");
        let class1 = trace(&b.mul(&inv.square(irr), irr), irr);
        if class1 && n1 >= targets {
            continue;
        }
        if decomposable(&elems, &x3, &b, irr) {
            skipped += 1;
            continue;
        }
        let eqs = descend_on_subspace(n, &v, irr, &b, &x3);
        let prof = rank_profile(&eqs, 2 * l, n, prof_dmax);
        let pt = measure_on_subspace(n, &v, irr, &b, &x3, scan_dmax);
        if class1 {
            n1 += 1;
            match pt.refutation_degree {
                Some(d) => ds1.push(d),
                None => censored1 += 1,
            }
        } else {
            n0 += 1;
            de.push(early_defect(&prof, 3));
            match pt.refutation_degree {
                Some(d) => ds0.push(d),
                None => censored0 += 1,
            }
        }
    }
    Cell {
        n,
        l,
        shape,
        tag,
        detail,
        irr: irr_name,
        irr_weight,
        dim_vv,
        dim_vvv,
        dim_v_sigv,
        dim_v_cap_sigv,
        yield_rand_b: hits_rand as f64 / yield_draws as f64,
        yield_b1: hits_b1 as f64 / yield_draws as f64,
        yield_draws,
        trace_zero,
        early_defect: if de.is_empty() {
            f64::NAN
        } else {
            de.iter().sum::<f64>() / de.len() as f64
        },
        dstar0_mean: mean_of(&ds0),
        dstar0_hist: hist_of(&ds0),
        n0_measured: n0,
        censored0,
        dstar1_mean: mean_of(&ds1),
        dstar1_hist: hist_of(&ds1),
        n1_measured: n1,
        censored1,
        skipped_decomposable: skipped,
        seconds: t0.elapsed().as_secs_f64(),
    }
}

fn spearman(xs: &[f64], ys: &[f64]) -> Option<f64> {
    if xs.len() < 3 || xs.len() != ys.len() {
        return None;
    }
    let rank = |v: &[f64]| -> Vec<f64> {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
        let mut r = vec![0.0; v.len()];
        let mut i = 0;
        while i < idx.len() {
            let mut j = i;
            while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
                j += 1;
            }
            let avg = (i + j) as f64 / 2.0 + 1.0;
            for &k in &idx[i..=j] {
                r[k] = avg;
            }
            i = j + 1;
        }
        r
    };
    let rx = rank(xs);
    let ry = rank(ys);
    let m = rx.len() as f64;
    let mx = rx.iter().sum::<f64>() / m;
    let my = ry.iter().sum::<f64>() / m;
    let cov: f64 = rx.iter().zip(&ry).map(|(a, b)| (a - mx) * (b - my)).sum();
    let vx: f64 = rx.iter().map(|a| (a - mx).powi(2)).sum();
    let vy: f64 = ry.iter().map(|b| (b - my).powi(2)).sum();
    if vx == 0.0 || vy == 0.0 {
        None
    } else {
        Some(cov / (vx * vy).sqrt())
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let seed: u64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(20261001);
    let scale: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let out_path = args
        .get(3)
        .cloned()
        .unwrap_or_else(|| "experiments/factor_base_closure_sweep.json".to_string());
    let moduli_limit: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(400);

    // (n, l, targets, d_cap, yield_draws). Critical operating points
    // 2l ∈ {n, n−1}, including prime n where only coord/random/gp/fp exist.
    let points: &[(u32, u32, u32, u32, u32)] = &[
        (7, 3, 20, 8, 200),
        (8, 4, 20, 8, 200),
        (9, 4, 20, 8, 200),
        (10, 5, 20, 7, 200),
        (11, 5, 20, 7, 200),
        (12, 6, 12, 6, 100),
        (13, 6, 12, 6, 100),
        (14, 7, 8, 6, 60),
        (15, 7, 8, 6, 60),
        (16, 8, 6, 5, 40),
    ];
    let shapes = [
        Shape::Subfield,
        Shape::Coord,
        Shape::Gp,
        Shape::GpScaled,
        Shape::GpSym,
        Shape::GpOrd,
        Shape::Mix,
        Shape::Fp,
        Shape::Fp2,
        Shape::Random,
        Shape::TzCoord,
        Shape::TzGp,
        Shape::TzFp,
        Shape::TzRandom,
    ];

    let mut rng = StdRng::seed_from_u64(seed);
    let mut cells: Vec<Cell> = Vec::new();
    let t_all = Instant::now();

    println!("factor_base_closure_sweep  seed={seed}  scale={scale}  moduli_limit={moduli_limit}");
    println!(
        "{:>3} {:>2} {:<10} {:>4} {:>4} {:>5} {:>5} {:>6} {:>6} {:>7} {:>6} {:>4} {:>6} {:>4} {:>2} {:>6}  detail",
        "n", "l", "tag", "VV", "VVV", "VσV", "V∩σV", "yld_b", "yld_1", "defect", "D*0", "n0", "D*1", "n1", "tz", "sec"
    );
    let print_cell = |cell: &Cell| {
        let fmt = |d: Option<f64>| d.map(|d| format!("{d:.3}")).unwrap_or_else(|| "-".into());
        println!(
            "{:>3} {:>2} {:<10} {:>4} {:>4} {:>5} {:>5} {:>6.3} {:>6.3} {:>7.4} {:>6} {:>4} {:>6} {:>4} {:>2} {:>6.1}  {} [{}]",
            cell.n,
            cell.l,
            cell.tag,
            cell.dim_vv,
            cell.dim_vvv,
            cell.dim_v_sigv,
            cell.dim_v_cap_sigv,
            cell.yield_rand_b,
            cell.yield_b1,
            cell.early_defect,
            fmt(cell.dstar0_mean),
            cell.n0_measured,
            fmt(cell.dstar1_mean),
            cell.n1_measured,
            if cell.trace_zero { "y" } else { "n" },
            cell.seconds,
            cell.detail,
            cell.irr
        );
    };
    for &(n, l, targets, d_cap, yield_draws) in points {
        let targets = ((targets as f64 * scale).round() as u32).max(2);
        let yield_draws = ((yield_draws as f64 * scale).round() as u32).max(10);
        let irrs = enumerate_irreducibles(n, moduli_limit);
        let Some(default_irr) = irrs.first().cloned() else {
            continue;
        };
        // Modulus panel for the coordinate shape: every trinomial, the
        // lightest pentanomial, and the two heaviest moduli enumerated.
        // Because D* is invariant under F₂-recombination of the equations,
        // "coord under modulus f" is "the progression whose ratio has
        // minimal polynomial f": the modulus panel is a ratio panel.
        let mut panel: Vec<IrreduciblePoly> = Vec::new();
        let push = |f: &IrreduciblePoly, panel: &mut Vec<IrreduciblePoly>| {
            if !panel.iter().any(|g| g.low_terms == f.low_terms) {
                panel.push(f.clone());
            }
        };
        for f in irrs.iter().filter(|f| f.low_terms.len() == 2) {
            push(f, &mut panel);
        }
        if let Some(f) = irrs.iter().find(|f| f.low_terms.len() == 4) {
            push(f, &mut panel);
        }
        let mut by_weight: Vec<&IrreduciblePoly> = irrs.iter().collect();
        by_weight.sort_by_key(|f| std::cmp::Reverse(f.low_terms.len()));
        for f in by_weight.iter().take(2) {
            push(f, &mut panel);
        }
        push(&default_irr, &mut panel);
        for (mi, irr) in panel.iter().enumerate() {
            if irr.low_terms == default_irr.low_terms {
                continue; // run below with the full shape set
            }
            let Some((cols, detail)) = draw_shape(Shape::Coord, n, l, irr, &mut rng) else {
                continue;
            };
            let cell = run_cell(
                n,
                l,
                Shape::Coord,
                format!("coord_m{mi}"),
                cols,
                detail,
                irr,
                targets,
                d_cap,
                yield_draws,
                &mut rng,
            );
            print_cell(&cell);
            cells.push(cell);
        }
        let irr = &default_irr;
        for &shape in &shapes {
            for k in 0..shape.reseeds() {
                let Some((cols, detail)) = draw_shape(shape, n, l, irr, &mut rng) else {
                    continue;
                };
                let tag = if shape.reseeds() == 1 {
                    shape.name().to_string()
                } else {
                    format!("{}{}", shape.name(), k)
                };
                let cell = run_cell(
                    n,
                    l,
                    shape,
                    tag,
                    cols,
                    detail,
                    irr,
                    targets,
                    d_cap,
                    yield_draws,
                    &mut rng,
                );
                print_cell(&cell);
                cells.push(cell);
            }
        }
    }

    // Pooled analysis over non-subfield cells with a D*.
    let pool: Vec<&Cell> = cells
        .iter()
        .filter(|c| c.shape != Shape::Subfield && c.dstar0_mean.is_some())
        .collect();
    let d: Vec<f64> = pool.iter().map(|c| c.dstar0_mean.unwrap()).collect();
    let vv: Vec<f64> = pool
        .iter()
        .map(|c| c.dim_vv as f64 - (2 * c.l) as f64 + 1.0)
        .collect();
    let cap: Vec<f64> = pool.iter().map(|c| c.dim_v_cap_sigv as f64).collect();
    let de: Vec<f64> = pool.iter().map(|c| c.early_defect).collect();
    let rho_vv = spearman(&vv, &d);
    let rho_cap = spearman(&cap, &d);
    let rho_de = spearman(&de, &d);
    println!();
    println!(
        "pooled non-subfield cells (class 0): {}  spearman(D*0, dimVV−(2l−1)) = {:?}  spearman(D*0, dim V∩σV) = {:?}  spearman(D*0, early_defect) = {:?}",
        pool.len(),
        rho_vv,
        rho_cap,
        rho_de
    );

    // Per-point ranking by mean D* (shape means over reseeds).
    let mut ranking = Vec::new();
    for &(n, l, _, _, _) in points {
        let mut by_shape: Vec<(String, f64, usize)> = Vec::new();
        for &shape in &shapes {
            let vals: Vec<f64> = cells
                .iter()
                .filter(|c| {
                    c.n == n && c.l == l && c.shape == shape && !c.tag.starts_with("coord_m")
                })
                .filter_map(|c| c.dstar0_mean)
                .collect();
            if !vals.is_empty() {
                by_shape.push((
                    shape.name().to_string(),
                    vals.iter().sum::<f64>() / vals.len() as f64,
                    vals.len(),
                ));
            }
        }
        by_shape.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let line: Vec<String> = by_shape
            .iter()
            .map(|(s, d, k)| format!("{s}={d:.2}({k})"))
            .collect();
        println!("rank(class 0) n={n} l={l}: {}", line.join("  "));
        ranking.push(json!({ "n": n, "l": l, "by_mean_dstar0": by_shape.iter().map(|(s, d, k)| json!({"shape": s, "mean_dstar0": d, "cells": k})).collect::<Vec<_>>() }));
    }

    let generated_at = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let doc = json!({
        "schema": "factor_base_closure_sweep/v2",
        "seed": seed,
        "scale": scale,
        "moduli_limit": moduli_limit,
        "generated_at": generated_at,
        "total_seconds": t_all.elapsed().as_secs_f64(),
        "protocol": {
            "system": "binary Semaev S3 Weil-descended, X1,X2 in V, random b and x3 (ffd_harness::weil_descend_s3_subspace via descent_lowgamma)",
            "dstar": "refutation degree over non-decomposable targets only, split by trace class t = Tr(sqrt(b)/x3); censored = no refutation by d_cap",
            "trace_identity": "Tr(S3(X1,X2,x3)/x3^2) = Tr(X1) + Tr(X2) + Tr(sqrt(b)/x3) for all X1, X2",
            "yield": "fraction of random x3 with some (x1,x2) in V^2 solving S3; random b and b=1",
            "growth": "dim V·V, dim V³, dim V·σV, dim V∩σV over F2",
            "points": points.iter().map(|&(n,l,t,d,y)| json!({"n":n,"l":l,"targets":t,"d_cap":d,"yield_draws":y})).collect::<Vec<_>>(),
        },
        "pooled": {
            "cells": pool.len(),
            "spearman_dstar_vs_vv_excess": rho_vv,
            "spearman_dstar_vs_v_cap_sigv": rho_cap,
            "spearman_dstar_vs_early_defect": rho_de,
        },
        "ranking": ranking,
        "cells": cells.iter().map(|c| json!({
            "n": c.n, "l": c.l, "shape": c.shape.name(), "tag": c.tag, "detail": c.detail,
            "irr": c.irr, "irr_weight": c.irr_weight,
            "dim_vv": c.dim_vv, "dim_vvv": c.dim_vvv, "dim_v_sigv": c.dim_v_sigv, "dim_v_cap_sigv": c.dim_v_cap_sigv,
            "yield_rand_b": c.yield_rand_b, "yield_b1": c.yield_b1, "yield_draws": c.yield_draws,
            "trace_zero": c.trace_zero, "early_defect": c.early_defect,
            "dstar0_mean": c.dstar0_mean, "dstar0_hist": c.dstar0_hist.iter().map(|(d,k)| json!([d,k])).collect::<Vec<_>>(),
            "n0_measured": c.n0_measured, "censored0": c.censored0,
            "dstar1_mean": c.dstar1_mean, "dstar1_hist": c.dstar1_hist.iter().map(|(d,k)| json!([d,k])).collect::<Vec<_>>(),
            "n1_measured": c.n1_measured, "censored1": c.censored1,
            "skipped_decomposable": c.skipped_decomposable,
            "seconds": c.seconds,
        })).collect::<Vec<_>>(),
    });
    if let Some(dir) = std::path::Path::new(&out_path).parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    std::fs::write(&out_path, serde_json::to_string_pretty(&doc).unwrap())
        .unwrap_or_else(|e| panic!("write {out_path}: {e}"));
    println!("wrote {out_path}  ({:.1}s)", t_all.elapsed().as_secs_f64());
}
