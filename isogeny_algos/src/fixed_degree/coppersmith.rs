//! Coron (bivariate) and Bauer–Joux-setting (trivariate) integer Coppersmith.
//!
//! Step 4 of the fixed-degree isogeny algorithm of Benčina–Kutas–Merz–Petit–Stopar–Weitkämper
//! (2023): after guessing some coordinates of a solution of `Q(x) = d`, the remaining small
//! coordinates are a small *integer* root (not a modular one) of a degree-2 polynomial over Z.
//!
//! # The lattice (Coron, CRYPTO 2007, "Finding small roots of bivariate integer polynomial
//! equations: a direct approach")
//!
//! Let `p ∈ Z[x_1..x_v]` with `p(0) = p_00 ≠ 0` and bounds `X_i` with `gcd(X_i, p_00) = 1` (the
//! solvers first apply a shift `x ↦ x + a`, `a_i ∈ {0, ±1, ±2, ±3}`, chosen so that `p(a) ≠ 0`
//! and the bounds, enlarged by `|a_i|` and bumped up to be coprime to `p(a)`, are smallest).  Let
//! `W = ||p(x_1 X_1, .., x_v X_v)||_∞`, `u = W + ((1 - W) mod |p_00|)` (so `u ≡ 1 mod p_00`),
//! `n = u · Π X_i^k` and `q = p_00^{-1} p mod n` (so `q(0) = 1`).  With a set `S` of shift
//! monomials and `M = S + supp(p)` the lattice is spanned by the coefficient vectors of
//!
//! * `g_s(x) = x^s · Π X_i^{k - s_i} · q(x)` for `s ∈ S`,
//! * `n · x^m` for `m ∈ M \ S`,
//!
//! each evaluated at `x_i ↦ x_i X_i`.  Every lattice polynomial `h` satisfies `h(r) ≡ 0 mod n`
//! at every integer root `r` of `p`.  Coron's choice is the box `S = [0, k]^2` for a
//! polynomial of degree `δ` in each variable; this module also offers the simplex
//! `S = {s : Σ s_i <= k}`, which is the natural shape for the total-degree-2 polynomials of
//! the norm equation (Coron's own extension to other Newton polygons follows the same pattern;
//! cf. Jochemsz–May 2006 for the general integer strategy).
//!
//! **Certificate.**  For each LLL vector we test exactly `||h(xX)||_1 < n`.  Then
//! `|h(r)| <= ||h(xX)||_1 < n` for every `r` with `|r_i| <= X_i`, so `h(r) = 0` for *every*
//! root of `p` in the box, not just a planted one.  Only certified vectors are used.
//!
//! **Bound.**  Coron proves that for `p` irreducible with degree `δ` in each variable, the
//! box lattice yields such an `h` that is not a multiple of `p` as soon as
//! `XY < W^{2/(3δ) - ε}` with `ε → 0` as `k → ∞` (lattice dimension `(δ + k + 1)^2`), giving
//! all small roots in time polynomial in `log W` and `2^{1/ε}`.  For `δ = 2` this is
//! `XY < W^{1/3}`.  For the simplex shape and total degree `δ` the reference point is Coppersmith's
//! classical total-degree bound `XY < W^{1/δ}` (`W^{1/2}` here); for three variables the
//! tests compare against `XYZ < W^{1/2}`.  **These two simplex figures are quoted, not derived
//! for this implementation**, and they are far from reached at the small `k` that are
//! practical here (see the measured table: at `k = 4` the bivariate simplex lattice reaches
//! `X ≈ 2^9` where `W^{1/2}` would allow `2^15`).  All bounds are only approached as `k → ∞`; [`lattice_margin_bits`] gives the
//! finite-`k` determinant heuristic and the tests print the measured limits.
//!
//! # Root extraction
//!
//! * Bivariate ([`coron_bivariate`]): for the first certified `h` with `Res_y(h, p) ≠ 0`
//!   (exact, Bareiss over Z[x]) all integer roots `x0` in range are found exactly
//!   ([`upoly::integer_roots`]), then `y0` from `p(x0, y)`.  When `p` is irreducible over Q
//!   and `h` is not a multiple of `p` the resultant is nonzero (Coron's theorem guarantees such
//!   an `h` inside the bound).  Given a certified `h` with nonzero resultant the output is
//!   *complete*: every root in the box is returned (the report says whether this held).
//! * Trivariate ([`trivariate_small_roots`]): two certified `h1, h2`; `R_i = Res_z(p, h_i)`,
//!   `V = Res_y(R_1, R_2) ∈ Z[x]`.  `V` is computed modulo one or more 61-bit primes by
//!   evaluation/interpolation (formal Sylvester degrees, Bézout degree bounds
//!   `tdeg Res_z(f, g) <= tdeg f · tdeg g`), its roots modulo the primes are found and
//!   combined by CRT into candidates `|x0| < X`, and each candidate is then solved exactly
//!   (`Res_z(p(x0,·,·), h(x0,·,·))` over Z, integer roots, back-substitution).  `V mod P` lies in
//!   the ideal `(p, h_1, h_2) mod P`, so every root in the box appears among the candidates
//!   whenever `V mod P ≠ 0`; spurious candidates are eliminated by the exact solve, and every
//!   returned root is verified by exact evaluation.
//!
//! **Heuristic part (trivariate).**  Whether `p, h1, h2` are algebraically independent (i.e.
//! `V ≠ 0`) is *not* guaranteed; it is checked at run time (`V mod P ≠ 0` for some prime; up
//! to `max_pairs` pairs of certified vectors are tried).  Bauer–Joux (EUROCRYPT 2007) make
//! this rigorous by building the second polynomial from a Gröbner basis of the ideal
//! `(p, h1)` so that independence is guaranteed; **that construction is not implemented
//! here**.  If no independent pair is found the report says so and no roots are claimed.
//!
//! # Norm equations
//!
//! [`solve_norm_eq_coppersmith`] uses the convention `Q(x) = Σ_{i,j} G[i][j] x_i x_j` (so a
//! symmetric `G` has cross coefficient `2 G[i][j]`; `G` need not be symmetric), guesses the
//! last one (trivariate) or two (bivariate) coordinates in `[-bound, bound]` (up to the sign
//! symmetry `x ↦ -x`), skips guesses for which the real conic/quadric has no real point
//! (exact Schur-complement test), and runs the Coppersmith solver for the rest with
//! `X_i = bound + 1`.
use crate::int::Int;
use crate::lattice::{lll_with, LllParams};
use crate::mpoly::{upoly, MPoly};
use std::collections::{BTreeSet, HashMap};
use std::time::{Duration, Instant};

fn int(v: i64) -> Int {
    Int::from(v)
}

/// Shape of the shift-monomial set `S`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    /// `S = [0, k]^v` (Coron 2007).
    Box,
    /// `S = {s : Σ s_i <= k}` (total-degree shape).
    Simplex,
}

/// Parameters of the Coron-style lattice.
#[derive(Clone, Debug)]
pub struct CoppersmithParams {
    /// Shift parameter `k`.
    pub k: u32,
    pub shape: Shape,
    pub lll: LllParams,
    /// Trivariate: number of pairs of certified vectors to try.
    pub max_pairs: usize,
    /// [`solve_norm_eq_coppersmith_with`] only: if too few vectors pass the certificate at
    /// `k`, retry with `k + 1, .., k_max`.
    pub k_max: u32,
}

impl CoppersmithParams {
    pub fn new(k: u32, shape: Shape) -> Self {
        CoppersmithParams {
            k,
            shape,
            lll: LllParams::default(),
            max_pairs: 6,
            k_max: k,
        }
    }
}

/// Diagnostics of one solver call.
#[derive(Clone, Debug, Default)]
pub struct SolveReport {
    pub dim: usize,
    /// Bounds actually used (after making them coprime to `p(0)` and after a shift).
    pub bounds: Vec<Int>,
    pub w_bits: usize,
    pub n_bits: usize,
    /// Number of LLL vectors with `||h(xX)||_1 < n`.
    pub certified: usize,
    /// All roots in the box are guaranteed to have been returned (certificate + nonzero
    /// elimination polynomial + every candidate solved completely).
    pub complete: bool,
    pub note: String,
    pub lll_time: Duration,
    pub total_time: Duration,
}

/// Certified short polynomials from the lattice.
#[derive(Clone, Debug)]
pub struct SmallRootPolys {
    pub n: Int,
    pub bounds: Vec<Int>,
    /// Certified polynomials (in the original, unscaled variables), LLL order.
    pub polys: Vec<MPoly>,
    pub dim: usize,
    pub w_bits: usize,
    pub lll_time: Duration,
}

fn log2(x: &Int) -> f64 {
    let b = x.bits();
    if b == 0 {
        return f64::NEG_INFINITY;
    }
    if b <= 60 {
        return x.abs().to_f64().log2();
    }
    x.abs().shr(b - 60).to_f64().log2() + (b - 60) as f64
}

fn shift_set(nv: usize, k: u32, shape: Shape) -> Vec<Vec<u32>> {
    let mut out = vec![];
    let mut e = vec![0u32; nv];
    loop {
        let ok = match shape {
            Shape::Box => true,
            Shape::Simplex => e.iter().sum::<u32>() <= k,
        };
        if ok {
            out.push(e.clone());
        }
        let mut i = 0;
        loop {
            if i == nv {
                return out;
            }
            e[i] += 1;
            if e[i] <= k {
                break;
            }
            e[i] = 0;
            i += 1;
        }
    }
}

fn pow_table(b: &[Int], max: u32) -> Vec<Vec<Int>> {
    b.iter()
        .map(|x| {
            let mut v = vec![Int::one()];
            for _ in 0..max {
                let t = v.last().unwrap() * x;
                v.push(t);
            }
            v
        })
        .collect()
}

fn mono_val(pw: &[Vec<Int>], e: &[u32]) -> Int {
    let mut t = Int::one();
    for (i, &k) in e.iter().enumerate() {
        if k > 0 {
            t = &t * &pw[i][k as usize];
        }
    }
    t
}

/// Bounds made coprime to `c` (each bumped up as needed).
fn coprime_bounds(bounds: &[Int], c: &Int) -> Vec<Int> {
    bounds
        .iter()
        .map(|x| {
            let mut x = x.clone();
            if x < Int::one() {
                x = Int::one();
            }
            while x.gcd(c) != Int::one() {
                x = &x + &Int::one();
            }
            x
        })
        .collect()
}

/// `(W, u, n)` for the lattice (see the module docs); `bounds` must be coprime to `p(0) ≠ 0`.
fn modulus(p: &MPoly, bounds: &[Int], k: u32) -> (Int, Int) {
    let p00 = p.constant_term();
    let w = p.scale_vars(bounds).norm_inf();
    let u = &w + &(&Int::one() - &w).modulo(&p00.abs());
    let mut n = u;
    for x in bounds {
        n = &n * &x.pow(k);
    }
    (w, n)
}

/// Build the Coron-style lattice for `p` (with `p(0) ≠ 0`), LLL-reduce it, and return the
/// certified polynomials (see the module docs).
pub fn coron_lattice_polys(
    p: &MPoly,
    bounds: &[Int],
    params: &CoppersmithParams,
) -> SmallRootPolys {
    let nv = p.nvars();
    assert_eq!(bounds.len(), nv);
    let p00 = p.constant_term();
    assert!(!p00.is_zero(), "coron_lattice_polys needs p(0) != 0");
    let k = params.k;
    let bounds = coprime_bounds(bounds, &p00);
    let (w, n) = modulus(p, &bounds, k);
    let inv = p00.inv_mod(&n).expect("gcd(n, p00) = 1 by construction");
    let half = n.div_floor(&int(2));
    let mut q = MPoly::zero(nv);
    for (e, c) in p.terms() {
        let mut v = (c * &inv).modulo(&n);
        if v > half {
            v = &v - &n;
        }
        q.add_term(e, &v);
    }
    let s_set = shift_set(nv, k, params.shape);
    let supp: Vec<Vec<u32>> = p.terms().map(|(e, _)| e.clone()).collect();
    let mut m_set: BTreeSet<Vec<u32>> = BTreeSet::new();
    for s in &s_set {
        for t in &supp {
            m_set.insert(s.iter().zip(t).map(|(a, b)| a + b).collect());
        }
    }
    for s in &s_set {
        m_set.insert(s.clone());
    }
    let cols: Vec<Vec<u32>> = m_set.into_iter().collect();
    let col_of: HashMap<Vec<u32>, usize> = cols
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, e)| (e, i))
        .collect();
    let maxe = cols
        .iter()
        .flat_map(|e| e.iter().copied())
        .max()
        .unwrap_or(0)
        + k;
    let pw = pow_table(&bounds, maxe);
    let s_lookup: BTreeSet<Vec<u32>> = s_set.iter().cloned().collect();
    let mut rows: Vec<Vec<Int>> = vec![];
    let lvec = vec![k; nv];
    for s in &s_set {
        let mut row = vec![Int::zero(); cols.len()];
        for (t, qt) in q.terms() {
            let m: Vec<u32> = s.iter().zip(t).map(|(a, b)| a + b).collect();
            let lt: Vec<u32> = lvec.iter().zip(t).map(|(a, b)| a + b).collect();
            row[col_of[&m]] = qt * &mono_val(&pw, &lt);
        }
        rows.push(row);
    }
    for m in &cols {
        if !s_lookup.contains(m) {
            let mut row = vec![Int::zero(); cols.len()];
            row[col_of[m]] = &n * &mono_val(&pw, m);
            rows.push(row);
        }
    }
    let dim = rows.len();
    let t0 = Instant::now();
    let red = lll_with(&rows, &params.lll).basis;
    let lll_time = t0.elapsed();
    let mut polys = vec![];
    for v in &red {
        let mut l1 = Int::zero();
        for x in v {
            l1 = &l1 + &x.abs();
        }
        if l1 >= n {
            continue;
        }
        let mut h = MPoly::zero(nv);
        for (j, x) in v.iter().enumerate() {
            if !x.is_zero() {
                let (c, r) = x.divrem_trunc(&mono_val(&pw, &cols[j]));
                debug_assert!(r.is_zero());
                h.add_term(&cols[j], &c);
            }
        }
        polys.push(h);
    }
    SmallRootPolys {
        n,
        bounds,
        polys,
        dim,
        w_bits: w.bits(),
        lll_time,
    }
}

/// Finite-`k` determinant heuristic for the lattice of [`coron_lattice_polys`]:
/// `ω log2 n − log2 det L − ω log2 √ω − ω^2 log2(1.02)` (practical LLL factor `1.02^ω`).
/// Positive means the first LLL vector is expected to be certified.
pub fn lattice_margin_bits(p: &MPoly, bounds: &[Int], params: &CoppersmithParams) -> f64 {
    let nv = p.nvars();
    let (p, _, bounds) = choose_shift(p, bounds);
    let k = params.k;
    let (_, n) = modulus(&p, &bounds, k);
    let s_set = shift_set(nv, k, params.shape);
    let mut m_set: BTreeSet<Vec<u32>> = BTreeSet::new();
    for s in &s_set {
        for (t, _) in p.terms() {
            m_set.insert(s.iter().zip(t).map(|(a, b)| a + b).collect());
        }
        m_set.insert(s.clone());
    }
    let lb: Vec<f64> = bounds.iter().map(log2).collect();
    let ln = log2(&n);
    let xl: f64 = lb.iter().map(|b| b * k as f64).sum();
    let s_lookup: BTreeSet<Vec<u32>> = s_set.iter().cloned().collect();
    let mut ldet = 0.0;
    for m in &m_set {
        if s_lookup.contains(m) {
            ldet += xl;
        } else {
            ldet += ln + m.iter().zip(&lb).map(|(&e, b)| e as f64 * b).sum::<f64>();
        }
    }
    let w = m_set.len() as f64;
    w * ln - ldet - w * w.sqrt().log2() - w * w * 1.02f64.log2()
}

/// Choose a small shift `a` (entries in `{0, ±1, ±2, ±3}`) with `p(a) ≠ 0` that minimises the
/// product of the lattice bounds after they are enlarged by `|a_i|` and made coprime to `p(a)`
/// (a constant term divisible by many small primes would otherwise inflate the bounds).
/// Returns `(p(x + a), a, bounds')`.
fn choose_shift(p: &MPoly, bounds: &[Int]) -> (MPoly, Vec<Int>, Vec<Int>) {
    let nv = p.nvars();
    let cand = [0i64, 1, -1, 2, -2, 3, -3];
    let total = cand.len().pow(nv as u32);
    let mut best: Option<(f64, Vec<Int>, Vec<Int>)> = None;
    for idx in 0..total {
        let mut a = vec![];
        let mut t = idx;
        for _ in 0..nv {
            a.push(int(cand[t % cand.len()]));
            t /= cand.len();
        }
        let c = p.eval(&a);
        if c.is_zero() {
            continue;
        }
        let b: Vec<Int> = bounds.iter().zip(&a).map(|(x, s)| x + &s.abs()).collect();
        let b = coprime_bounds(&b, &c);
        let cost: f64 = b.iter().map(log2).sum::<f64>() + 1e-6 * log2(&c);
        if best.as_ref().is_none_or(|(bc, _, _)| cost < *bc) {
            best = Some((cost, a, b));
        }
        // the unshifted polynomial with already-coprime bounds cannot be beaten
        if idx == 0 && best.as_ref().unwrap().2 == bounds {
            break;
        }
    }
    let (_, a, b) = best.expect("polynomial vanishes on a small grid: zero polynomial?");
    let mut q = p.clone();
    for (i, ai) in a.iter().enumerate() {
        if !ai.is_zero() {
            q = q.shift_var(i, ai);
        }
    }
    (q, a, b)
}

fn in_box(r: &[Int], bounds: &[Int]) -> bool {
    r.iter().zip(bounds).all(|(x, b)| x.abs() < *b)
}

/// Coron's direct method (box lattice, parameter `k`): all integer roots `(x0, y0)` of `f` with
/// `|x0| < X`, `|y0| < Y` that the lattice finds (all of them when the report says
/// `complete`; see the module docs).  Every returned root satisfies `f(x0, y0) = 0` exactly.
pub fn coron_bivariate(f: &MPoly, x_bound: &Int, y_bound: &Int, k: u32) -> Vec<(Int, Int)> {
    coron_bivariate_with(f, x_bound, y_bound, &CoppersmithParams::new(k, Shape::Box)).0
}

/// [`coron_bivariate`] with explicit parameters (shape, LLL) and a report.
pub fn coron_bivariate_with(
    f: &MPoly,
    x_bound: &Int,
    y_bound: &Int,
    params: &CoppersmithParams,
) -> (Vec<(Int, Int)>, SolveReport) {
    assert_eq!(f.nvars(), 2);
    let t0 = Instant::now();
    let user = [x_bound.clone(), y_bound.clone()];
    let (p, a, b2) = choose_shift(f, &user);
    let lat = coron_lattice_polys(&p, &b2, params);
    let mut rep = SolveReport {
        dim: lat.dim,
        bounds: lat.bounds.clone(),
        w_bits: lat.w_bits,
        n_bits: lat.n.bits(),
        certified: lat.polys.len(),
        lll_time: lat.lll_time,
        ..Default::default()
    };
    let mut out: Vec<(Int, Int)> = vec![];
    let xb = &lat.bounds[0] - &Int::one();
    let yb = &lat.bounds[1] - &Int::one();
    let mut found_res = false;
    for h in &lat.polys {
        let r = h.resultant_univariate(&p, 1, 0).unwrap();
        if upoly::trim(r.clone()).is_empty() {
            continue;
        }
        found_res = true;
        let mut complete = true;
        for x0 in upoly::integer_roots(&r, Some(&xb)) {
            let py = p.subs_value(0, &x0).to_univariate(1).unwrap();
            let py = upoly::trim(py);
            if py.is_empty() {
                complete = false;
                rep.note = format!("p(x0, y) == 0 for x0 = {x0}: y not determined");
                continue;
            }
            for y0 in upoly::integer_roots(&py, Some(&yb)) {
                let r0 = [&x0 + &a[0], &y0 + &a[1]];
                if f.eval(&r0).is_zero() && in_box(&r0, &user) {
                    out.push((r0[0].clone(), r0[1].clone()));
                }
            }
        }
        rep.complete = complete;
        break;
    }
    if !found_res {
        rep.note = if lat.polys.is_empty() {
            "no LLL vector passed the certificate ||h(xX)||_1 < n (bound too large for k)".into()
        } else {
            "every certified h shares a factor with p (p reducible?)".into()
        };
    }
    out.sort();
    out.dedup();
    rep.total_time = t0.elapsed();
    (out, rep)
}

// ------------------------------------------------------------------------------------------
// modular helpers for the trivariate elimination

fn mm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
fn pw_mod(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1;
    while e > 0 {
        if e & 1 == 1 {
            r = mm(r, a, p);
        }
        a = mm(a, a, p);
        e >>= 1;
    }
    r
}
fn inv_mod(a: u64, p: u64) -> u64 {
    pw_mod(a, p - 2, p)
}

/// Sylvester resultant mod `p` with formal degrees (`f.len() - 1`, `g.len() - 1`).
fn res_mod(f: &[u64], g: &[u64], p: u64) -> u64 {
    let (a, b) = (f.len() - 1, g.len() - 1);
    if a == 0 {
        return pw_mod(f[0], b as u64, p);
    }
    if b == 0 {
        return pw_mod(g[0], a as u64, p);
    }
    let s = a + b;
    let mut m = vec![vec![0u64; s]; s];
    for r in 0..b {
        for (t, &c) in f.iter().rev().enumerate() {
            m[r][r + t] = c;
        }
    }
    for r in 0..a {
        for (t, &c) in g.iter().rev().enumerate() {
            m[b + r][r + t] = c;
        }
    }
    let mut det = 1u64;
    for c in 0..s {
        let Some(piv) = (c..s).find(|&r| m[r][c] != 0) else {
            return 0;
        };
        if piv != c {
            m.swap(piv, c);
            det = (p - det) % p;
        }
        det = mm(det, m[c][c], p);
        let iv = inv_mod(m[c][c], p);
        for r in c + 1..s {
            if m[r][c] != 0 {
                let fct = mm(m[r][c], iv, p);
                for j in c..s {
                    let t = mm(fct, m[c][j], p);
                    m[r][j] = (m[r][j] + p - t) % p;
                }
            }
        }
    }
    det
}

/// Coefficients (monomial basis) of the polynomial of degree `< vals.len()` taking value
/// `vals[i]` at `i` (mod `p`, `p > vals.len()`).
fn interpolate(vals: &[u64], p: u64) -> Vec<u64> {
    let n = vals.len();
    let mut dd = vals.to_vec();
    for j in 1..n {
        let iv = inv_mod(j as u64 % p, p);
        for i in (j..n).rev() {
            dd[i] = mm((dd[i] + p - dd[i - 1]) % p, iv, p);
        }
    }
    let mut poly = vec![dd[n - 1]];
    for j in (0..n - 1).rev() {
        // poly = poly * (x - j) + dd[j]
        let mut np = vec![0u64; poly.len() + 1];
        for (i, &c) in poly.iter().enumerate() {
            np[i + 1] = (np[i + 1] + c) % p;
            np[i] = (np[i] + p - mm(c, j as u64 % p, p)) % p;
        }
        np[0] = (np[0] + dd[j]) % p;
        poly = np;
    }
    while poly.last() == Some(&0) {
        poly.pop();
    }
    poly
}

struct ModPoly {
    terms: Vec<([u32; 3], u64)>,
    dz: usize,
}
impl ModPoly {
    fn new(f: &MPoly, p: u64) -> ModPoly {
        ModPoly {
            terms: f
                .terms()
                .map(|(e, c)| ([e[0], e[1], e[2]], c.mod_u64(p)))
                .collect(),
            dz: f.deg_in(2).unwrap_or(0) as usize,
        }
    }
    /// As a polynomial in z (formal degree dz) at x = a, y = b.
    fn at(&self, a: u64, b: u64, p: u64) -> Vec<u64> {
        let mut v = vec![0u64; self.dz + 1];
        for (e, c) in &self.terms {
            let t = mm(
                mm(*c, pw_mod(a, e[0] as u64, p), p),
                pw_mod(b, e[1] as u64, p),
                p,
            );
            v[e[2] as usize] = (v[e[2] as usize] + t) % p;
        }
        v
    }
}

/// `V = Res_y(Res_z(f, g1), Res_z(f, g2)) mod p` (a polynomial in x lying in the ideal
/// `(f, g1, g2) mod p`), by evaluation/interpolation.  `None` if `Res_z(f, g_i) ≡ 0 mod p`.
fn elim_mod(f: &MPoly, g1: &MPoly, g2: &MPoly, p: u64) -> Option<Vec<u64>> {
    let tf = f.total_degree().unwrap_or(0) as usize;
    let t1 = tf * g1.total_degree().unwrap_or(0) as usize;
    let t2 = tf * g2.total_degree().unwrap_or(0) as usize;
    let dmax = t1 * t2;
    let (mf, m1, m2) = (ModPoly::new(f, p), ModPoly::new(g1, p), ModPoly::new(g2, p));
    let mut r1s = vec![];
    let mut r2s = vec![];
    let (mut d1, mut d2) = (None::<usize>, None::<usize>);
    for a in 0..=dmax as u64 {
        let ty = t1.max(t2);
        let mut v1 = vec![];
        let mut v2 = vec![];
        for b in 0..=ty as u64 {
            let fz = mf.at(a, b, p);
            v1.push(res_mod(&fz, &m1.at(a, b, p), p));
            v2.push(res_mod(&fz, &m2.at(a, b, p), p));
        }
        let r1 = interpolate(&v1[..=t1], p);
        let r2 = interpolate(&v2[..=t2], p);
        if !r1.is_empty() {
            d1 = Some(d1.map_or(r1.len() - 1, |d: usize| d.max(r1.len() - 1)));
        }
        if !r2.is_empty() {
            d2 = Some(d2.map_or(r2.len() - 1, |d: usize| d.max(r2.len() - 1)));
        }
        r1s.push(r1);
        r2s.push(r2);
    }
    let (d1, d2) = (d1?, d2?);
    let vals: Vec<u64> = (0..=dmax)
        .map(|i| {
            let mut a = r1s[i].clone();
            a.resize(d1 + 1, 0);
            let mut b = r2s[i].clone();
            b.resize(d2 + 1, 0);
            res_mod(&a, &b, p)
        })
        .collect();
    Some(interpolate(&vals, p))
}

fn primes_below_2_61() -> impl Iterator<Item = u64> {
    let mut c = (1u64 << 61) - 1;
    std::iter::from_fn(move || {
        while !upoly::is_prime_u64(c) {
            c -= 2;
        }
        let r = c;
        c -= 2;
        Some(r)
    })
}

/// Candidates `|x0| < xb` for the first coordinate of the common roots of `f, g1, g2`
/// (complete, possibly with spurious extra values), or `None` if `V ≡ 0` modulo the primes
/// tried (dependent triple).
fn x_candidates(f: &MPoly, g1: &MPoly, g2: &MPoly, xb: &Int) -> Option<Vec<Int>> {
    // need prod(primes) > 2 xb * 2^20
    let need = xb.bits() + 22;
    let mut cands: Vec<Int> = vec![Int::zero()];
    let mut modulus = Int::one();
    let mut zero_hits = 0;
    let mut primes = primes_below_2_61();
    while modulus.bits() < need {
        let p = primes.next().unwrap();
        let v = match elim_mod(f, g1, g2, p) {
            Some(v) if !v.is_empty() => v,
            _ => {
                zero_hits += 1;
                if zero_hits >= 2 {
                    return None;
                }
                continue;
            }
        };
        let roots = upoly::roots_mod_p(&v, p);
        let pi = Int::from(p);
        let minv = modulus.modulo(&pi).inv_mod(&pi).unwrap();
        let mut next = vec![];
        for c in &cands {
            for &r in &roots {
                let t = (&(&Int::from(r) - c) * &minv).modulo(&pi);
                next.push(c + &(&modulus * &t));
            }
        }
        modulus = &modulus * &pi;
        cands = next;
        if cands.len() > 100_000 {
            return None;
        }
    }
    let half = modulus.div_floor(&int(2));
    let mut out: Vec<Int> = cands
        .into_iter()
        .map(|c| if c > half { &c - &modulus } else { c })
        .filter(|c| c.abs() <= *xb)
        .collect();
    out.sort();
    out.dedup();
    Some(out)
}

/// Small integer roots of a trivariate integer polynomial (any degree; degree 2 in the use
/// case): all `[x0, y0, z0]` with `|x0| < X`, `|y0| < Y`, `|z0| < Z` and `f = 0` found by the
/// Coron-style lattice of parameter `params.k`/`params.shape` plus resultant elimination
/// with a run-time independence check (heuristic; the rigorous Bauer–Joux Gröbner step is not
/// implemented — see the module docs).
pub fn trivariate_small_roots(
    f: &MPoly,
    bounds: &[Int; 3],
    params: &CoppersmithParams,
) -> Vec<[Int; 3]> {
    trivariate_small_roots_report(f, bounds, params).0
}

/// [`trivariate_small_roots`] with a report.
pub fn trivariate_small_roots_report(
    f: &MPoly,
    bounds: &[Int; 3],
    params: &CoppersmithParams,
) -> (Vec<[Int; 3]>, SolveReport) {
    assert_eq!(f.nvars(), 3);
    let t0 = Instant::now();
    let (p, a, b2) = choose_shift(f, bounds);
    let lat = coron_lattice_polys(&p, &b2, params);
    let mut rep = SolveReport {
        dim: lat.dim,
        bounds: lat.bounds.clone(),
        w_bits: lat.w_bits,
        n_bits: lat.n.bits(),
        certified: lat.polys.len(),
        lll_time: lat.lll_time,
        ..Default::default()
    };
    let lim: Vec<Int> = lat.bounds.iter().map(|x| x - &Int::one()).collect();
    let mut out: Vec<[Int; 3]> = vec![];
    let hs = &lat.polys;
    let mut pairs = vec![];
    for j in 1..hs.len() {
        for i in 0..j {
            pairs.push((i, j));
        }
    }
    pairs.truncate(params.max_pairs);
    let mut solved = false;
    for (i, j) in pairs {
        let Some(xs) = x_candidates(&p, &hs[i], &hs[j], &lim[0]) else {
            continue;
        };
        solved = true;
        let mut complete = true;
        let polys = [&p, &hs[i], &hs[j]];
        for x0 in xs {
            let sub: Vec<MPoly> = polys.iter().map(|q| q.subs_value(0, &x0)).collect();
            // y: first nonzero Res_z among the pairs
            let mut ry = None;
            for (u, v) in [(0, 1), (0, 2), (1, 2)] {
                let r = upoly::trim(sub[u].resultant_univariate(&sub[v], 2, 1).unwrap());
                if !r.is_empty() {
                    ry = Some(r);
                    break;
                }
            }
            let Some(ry) = ry else {
                complete = false;
                rep.note = format!("y not determined for x0 = {x0}");
                continue;
            };
            for y0 in upoly::integer_roots(&ry, Some(&lim[1])) {
                let mut zp = None;
                for s in &sub {
                    let z = upoly::trim(s.subs_value(1, &y0).to_univariate(2).unwrap());
                    if !z.is_empty() {
                        zp = Some(z);
                        break;
                    }
                }
                let Some(zp) = zp else {
                    complete = false;
                    rep.note = format!("z not determined for (x0, y0) = ({x0}, {y0})");
                    continue;
                };
                for z0 in upoly::integer_roots(&zp, Some(&lim[2])) {
                    let r0 = [&x0 + &a[0], &y0 + &a[1], &z0 + &a[2]];
                    if f.eval(&r0).is_zero() && in_box(&r0, bounds) {
                        out.push(r0);
                    }
                }
            }
        }
        rep.complete = complete;
        break;
    }
    if !solved {
        rep.note = if hs.len() < 2 {
            format!(
                "only {} LLL vector(s) passed the certificate (need 2)",
                hs.len()
            )
        } else {
            "no algebraically independent pair (p, h_i, h_j) among the certified vectors".into()
        };
    }
    out.sort();
    out.dedup();
    rep.total_time = t0.elapsed();
    (out, rep)
}

// ------------------------------------------------------------------------------------------
// norm equations

/// Which coordinates are guessed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NormEqMode {
    /// Guess `x3, x4`; bivariate Coppersmith in `x1, x2`.
    Bivariate,
    /// Guess `x4`; trivariate Coppersmith in `x1, x2, x3`.
    Trivariate,
}

/// `Q(x) = Σ_{i,j} G[i][j] x_i x_j` evaluated exactly.
pub fn quad_form(gram: &[[Int; 4]; 4], x: &[Int; 4]) -> Int {
    let mut s = Int::zero();
    for i in 0..4 {
        for j in 0..4 {
            s = &s + &(&(&gram[i][j] * &x[i]) * &x[j]);
        }
    }
    s
}

/// The polynomial `Q(x_free, guessed) − d` in the free variables (`nfree` of them, the first
/// coordinates).
pub fn norm_eq_poly(gram: &[[Int; 4]; 4], d: &Int, guessed: &[Int]) -> MPoly {
    let nfree = 4 - guessed.len();
    let val = |i: usize| -> Option<&Int> {
        if i < nfree {
            None
        } else {
            Some(&guessed[i - nfree])
        }
    };
    let mut f = MPoly::constant(nfree, &-d);
    for i in 0..4 {
        for j in 0..4 {
            let g = &gram[i][j];
            let mut e = vec![0u32; nfree];
            let mut c = g.clone();
            match val(i) {
                None => e[i] += 1,
                Some(v) => c = &c * v,
            }
            match val(j) {
                None => e[j] += 1,
                Some(v) => c = &c * v,
            }
            f.add_term(&e, &c);
        }
    }
    f
}

fn det_int(mut m: Vec<Vec<Int>>) -> Int {
    let n = m.len();
    let mut sign = false;
    let mut prev = Int::one();
    for k in 0..n {
        if m[k][k].is_zero() {
            match (k + 1..n).find(|&r| !m[r][k].is_zero()) {
                Some(r) => {
                    m.swap(k, r);
                    sign = !sign;
                }
                None => return Int::zero(),
            }
        }
        for i in k + 1..n {
            for j in k + 1..n {
                let t = &(&m[i][j] * &m[k][k]) - &(&m[i][k] * &m[k][j]);
                m[i][j] = t.divrem_trunc(&prev).0;
            }
        }
        prev = m[k][k].clone();
    }
    if sign {
        -prev
    } else {
        prev
    }
}

/// Whether `Q(u, guessed) = d` has a real solution `u` (exact; needs `Q` positive definite on
/// the free coordinates): `min_u 2(Q − d) = det(H) / det(A2)` with `A2 = (G + G^T)_free`.
fn real_solvable(gram: &[[Int; 4]; 4], d: &Int, guessed: &[Int]) -> bool {
    let nfree = 4 - guessed.len();
    let sym = |i: usize, j: usize| &gram[i][j] + &gram[j][i];
    let mut h = vec![vec![Int::zero(); nfree + 1]; nfree + 1];
    for a in 0..nfree {
        for b in 0..nfree {
            h[a][b] = sym(a, b);
        }
        let mut s = Int::zero();
        for (c, g) in guessed.iter().enumerate() {
            s = &s + &(&sym(a, nfree + c) * g);
        }
        h[a][nfree] = s.clone();
        h[nfree][a] = s;
    }
    let mut full = [Int::zero(), Int::zero(), Int::zero(), Int::zero()];
    for (c, g) in guessed.iter().enumerate() {
        full[nfree + c] = g.clone();
    }
    h[nfree][nfree] = &(&quad_form(gram, &full) - d) * &int(2);
    det_int(h) <= Int::zero()
}

/// Default lattice parameters per mode (simplex shape, which matches total degree 2).
pub fn default_norm_eq_params(mode: NormEqMode) -> CoppersmithParams {
    match mode {
        NormEqMode::Bivariate => CoppersmithParams {
            k_max: 4,
            ..CoppersmithParams::new(2, Shape::Simplex)
        },
        NormEqMode::Trivariate => CoppersmithParams {
            k_max: 3,
            ..CoppersmithParams::new(1, Shape::Simplex)
        },
    }
}

/// Find `x ∈ Z^4` with `Q(x) = Σ G[i][j] x_i x_j = d` and `|x_i| <= bound`, by guessing the last
/// one or two coordinates and running the Coppersmith solver on the rest (default
/// parameters: [`default_norm_eq_params`]).  The result is verified exactly.  `gram` must be
/// positive definite (used for the real-solvability pruning).
pub fn solve_norm_eq_coppersmith(
    gram: &[[Int; 4]; 4],
    d: &Int,
    bound: &Int,
    mode: NormEqMode,
) -> Option<[Int; 4]> {
    solve_norm_eq_coppersmith_with(gram, d, bound, mode, &default_norm_eq_params(mode))
}

/// [`solve_norm_eq_coppersmith`] with explicit lattice parameters.
pub fn solve_norm_eq_coppersmith_with(
    gram: &[[Int; 4]; 4],
    d: &Int,
    bound: &Int,
    mode: NormEqMode,
    params: &CoppersmithParams,
) -> Option<[Int; 4]> {
    let b = bound
        .to_i128()
        .expect("bound must fit in i64 for the guessing loop") as i64;
    let xb = bound + &Int::one();
    let mut guesses: Vec<Vec<Int>> = vec![];
    match mode {
        NormEqMode::Trivariate => {
            for x4 in 0..=b {
                guesses.push(vec![int(x4)]);
            }
        }
        NormEqMode::Bivariate => {
            // (x3, x4) up to sign: x4 > 0, or x4 = 0 and x3 >= 0
            for x4 in 0..=b {
                for x3 in -b..=b {
                    if x4 == 0 && x3 < 0 {
                        continue;
                    }
                    guesses.push(vec![int(x3), int(x4)]);
                }
            }
        }
    }
    for g in guesses {
        if !real_solvable(gram, d, &g) {
            continue;
        }
        let f = norm_eq_poly(gram, d, &g);
        let mut sols: Vec<Vec<Int>> = vec![];
        let need = match mode {
            NormEqMode::Bivariate => 1,
            NormEqMode::Trivariate => 2,
        };
        for k in params.k..=params.k_max.max(params.k) {
            let pk = CoppersmithParams {
                k,
                ..params.clone()
            };
            let (s, rep) = match mode {
                NormEqMode::Bivariate => {
                    let (r, rep) = coron_bivariate_with(&f, &xb, &xb, &pk);
                    (
                        r.into_iter().map(|(x, y)| vec![x, y]).collect::<Vec<_>>(),
                        rep,
                    )
                }
                NormEqMode::Trivariate => {
                    let (r, rep) = trivariate_small_roots_report(
                        &f,
                        &[xb.clone(), xb.clone(), xb.clone()],
                        &pk,
                    );
                    (r.into_iter().map(|r| r.to_vec()).collect(), rep)
                }
            };
            sols = s;
            if rep.certified >= need {
                break;
            }
        }
        for s in sols {
            let mut x: Vec<Int> = s;
            x.extend(g.iter().cloned());
            let x: [Int; 4] = std::array::from_fn(|i| x[i].clone());
            if quad_form(gram, &x) == *d && x.iter().all(|v| v.abs() <= *bound) {
                return Some(x);
            }
        }
    }
    None
}
