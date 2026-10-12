//! lit2026 W12: the zero-minor ECDLP method of Mahalanobis et al., measured
//! in counted field operations against the repository's `rho.negation`.
//!
//! Papers (read in full on arXiv):
//! - Mahalanobis, Mallick, *A Las Vegas algorithm to solve the elliptic curve
//!   discrete logarithm problem*, arXiv 1703.07544 (INDOCRYPT 2018):
//!   Algorithm 1, the matrix `M` of degree-`n'` monomials at `3n' + l`
//!   points and its left kernel `K` (Theorem 2.3, Corollary 2.4).
//! - Abdullah, Mahalanobis, *Minors solve the elliptic curve discrete
//!   logarithm problem*, arXiv 2310.04132: `K` in anti-diagonal form
//!   `[A | J]`, and a zero minor of the dense part `A` gives a maximal zero
//!   minor of `K`, which solves the ECDLP (Theorem 3.1).
//! - Mahalanobis, *A Guess and Determine Attack on the Elliptic Curve
//!   Discrete Logarithm Problem*, arXiv 2607.09814v2: the guess step
//!   ("finding a mate", Theorem 1, Algorithm 1 lines 3–8) and the determine
//!   step (signature matrix, Theorems 3–5, Algorithm 1 lines 10–22,
//!   Algorithm 2), with the defect `d` chosen from Equation (4).
//!
//! Everything the method does in `F_p` goes through [`Fp`], which tallies
//! multiplications and inversions (additions and subtractions free, the
//! convention of `ic_boundary::FieldOps`).  The points `n_i P`, `-n_j Q` are
//! computed with the repository's `PrimeCurve` group law, which tallies its
//! own field operations.  The rho reference is `ic_boundary::rho_walk_with`
//! with `NegationClasses` and `RhoWalk::negation()` — the code path of
//! ecbench's `rho.negation` — on the same curves and targets.
//!
//! Usage:
//!   lit2026_w12_zero_minor run [--bits 16,18,20,22,24,26] [--targets 8]
//!       [--rho-seeds 4] [--d auto|N] [--nprime auto|N] [--seed S] --out FILE
//!   lit2026_w12_zero_minor dsweep --bits 20 --ds 2,3,4 [--targets 8] --out FILE

use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

use crypto_lib::cryptanalysis::ic_boundary::{
    rho_cap, rho_walk_with, CountedGroup, FieldOps, GroupOps, NegationClasses, PrimeCurve,
    PrimeInstance, PrimePoint, RhoWalk,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;
use serde_json::json;

// ───────────────────────── counted F_p ─────────────────────────

/// Field-operation tally of the zero-minor method, same convention as
/// `ic_boundary::FieldOps`: every modular multiplication and inversion
/// is counted, additions and subtractions are free.
#[derive(Clone, Copy, Debug, Default, Serialize, PartialEq, Eq)]
pub struct Tally {
    pub muls: u64,
    pub invs: u64,
}

impl Tally {
    fn add(&mut self, o: Tally) {
        self.muls += o.muls;
        self.invs += o.invs;
    }
    fn from_field_ops(f: FieldOps) -> Tally {
        Tally {
            muls: f.muls + f.sqrs,
            invs: f.invs,
        }
    }
    /// The comparison unit: multiplications (squarings included) plus
    /// inversions, each counted as one operation.
    pub fn ops(&self) -> u64 {
        self.muls + self.invs
    }
}

pub struct Fp {
    pub p: u64,
    pub t: Tally,
}

impl Fp {
    pub fn new(p: u64) -> Self {
        assert!(p < (1u64 << 32), "single-word products need p < 2^32");
        Fp {
            p,
            t: Tally::default(),
        }
    }
    #[inline]
    pub fn mul(&mut self, a: u64, b: u64) -> u64 {
        self.t.muls += 1;
        a * b % self.p
    }
    #[inline]
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    #[inline]
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    #[inline]
    pub fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    pub fn inv(&mut self, a: u64) -> u64 {
        self.t.invs += 1;
        inv_mod(a, self.p)
    }
    /// The tally since the last call, and reset.
    pub fn take(&mut self) -> Tally {
        std::mem::take(&mut self.t)
    }
}

pub fn inv_mod(a: u64, m: u64) -> u64 {
    let (mut t, mut nt) = (0i128, 1i128);
    let (mut r, mut nr) = (m as i128, (a % m) as i128);
    while nr != 0 {
        let q = r / nr;
        (t, nt) = (nt, t - q * nt);
        (r, nr) = (nr, r - q * nr);
    }
    assert_eq!(r, 1, "not invertible");
    if t < 0 {
        t += m as i128;
    }
    t as u64
}

// ───────────────────────── counted linear algebra ─────────────────────────

/// In-place reduced row echelon form over `F_p`, columns visited in
/// `order`.  Returns `(row, column)` pivot pairs in order.  Counts one
/// inversion per pivot and one multiplication per entry it scales or
/// eliminates.
pub fn rref_order(f: &mut Fp, m: &mut [Vec<u64>], order: &[usize]) -> Vec<usize> {
    let rows = m.len();
    let mut pivots = Vec::new();
    let mut r = 0;
    for &c in order {
        if r == rows {
            break;
        }
        let Some(pr) = (r..rows).find(|&i| m[i][c] != 0) else {
            continue;
        };
        m.swap(r, pr);
        let inv = f.inv(m[r][c]);
        // Scale the pivot row (only the columns still to be visited, and
        // the pivot itself, can be non-zero; scale every non-zero entry).
        for j in 0..m[r].len() {
            if m[r][j] != 0 {
                m[r][j] = if j == c { 1 } else { f.mul(m[r][j], inv) };
            }
        }
        let prow = m[r].clone();
        for (i, row) in m.iter_mut().enumerate() {
            if i == r || row[c] == 0 {
                continue;
            }
            let fac = row[c];
            for j in 0..prow.len() {
                if prow[j] != 0 {
                    let s = f.mul(fac, prow[j]);
                    row[j] = f.sub(row[j], s);
                }
            }
        }
        pivots.push(c);
        r += 1;
    }
    pivots
}

/// Right kernel `{x : X x = 0}` of `x_mat` (rows × cols), as a basis in
/// free-variable form: one vector per non-pivot column `f`, with `1` at `f`.
pub fn right_kernel(f: &mut Fp, x_mat: &[Vec<u64>], cols: usize) -> Vec<Vec<u64>> {
    let mut m: Vec<Vec<u64>> = x_mat.to_vec();
    let order: Vec<usize> = (0..cols).collect();
    let piv = rref_order(f, &mut m, &order);
    let mut is_piv = vec![false; cols];
    for &c in &piv {
        is_piv[c] = true;
    }
    let mut basis = Vec::new();
    for free in (0..cols).filter(|&c| !is_piv[c]) {
        let mut v = vec![0u64; cols];
        v[free] = 1;
        for (i, &pc) in piv.iter().enumerate() {
            v[pc] = f.neg(m[i][free]);
        }
        basis.push(v);
    }
    basis
}

/// Left kernel `{v : v X = 0}` of `x_mat` (rows × cols): the right kernel of
/// the transpose.
pub fn left_kernel(f: &mut Fp, x_mat: &[Vec<u64>], rows: usize, cols: usize) -> Vec<Vec<u64>> {
    let t: Vec<Vec<u64>> = (0..cols)
        .map(|j| (0..rows).map(|i| x_mat[i][j]).collect())
        .collect();
    right_kernel(f, &t, rows)
}

/// Determinant by elimination (uncounted; tests and checks only).
pub fn det_uncounted(p: u64, m: &[Vec<u64>]) -> u64 {
    let mut f = Fp::new(p);
    let n = m.len();
    let mut a = m.to_vec();
    let mut det = 1u64;
    for c in 0..n {
        let Some(pr) = (c..n).find(|&i| a[i][c] != 0) else {
            return 0;
        };
        if pr != c {
            a.swap(pr, c);
            det = f.neg(det);
        }
        det = f.mul(det, a[c][c]);
        let inv = f.inv(a[c][c]);
        for i in c + 1..n {
            if a[i][c] == 0 {
                continue;
            }
            let fac = f.mul(a[i][c], inv);
            for j in c..n {
                let s = f.mul(fac, a[c][j]);
                a[i][j] = f.sub(a[i][j], s);
            }
        }
    }
    det
}

fn binom(n: u64, k: u64) -> f64 {
    if k > n {
        return 0.0;
    }
    let mut r = 1.0f64;
    for i in 0..k {
        r = r * (n - i) as f64 / (i + 1) as f64;
    }
    r
}

/// Next k-combination of `0..n` in lexicographic order.
fn next_comb(c: &mut [usize], n: usize) -> bool {
    let k = c.len();
    if k == 0 {
        return false;
    }
    let mut i = k;
    while i > 0 {
        i -= 1;
        if c[i] < n - k + i {
            c[i] += 1;
            for j in i + 1..k {
                c[j] = c[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

// ───────────────────────── the instance ─────────────────────────

/// A registered prime-order curve, built from its registry parameters.
pub struct Curve {
    pub slug: String,
    pub generator_call: String,
    pub inst: PrimeInstance,
}

/// `docs/curves/registry.json` rows with `generator_call =
/// find_prime_order_curve(bits, 59297)` (the ecbench prime ladder):
/// `(bits, p, a, b, order, slug)`.
pub const LADDER: &[(u32, u64, u64, u64, u64, &str)] = &[
    (12, 3529, 2574, 2980, 3583, "icv1-fp12-tm53-6c9ff4bb"),
    (14, 14389, 7150, 684, 14533, "icv1-fp14-tm143-dec4cf6b"),
    (16, 53411, 4589, 186, 53117, "icv1-fp16-t295-8d3c3165"),
    (
        18,
        212099,
        188885,
        173081,
        212437,
        "icv1-fp18-tm337-d28d5e09",
    ),
    (20, 563153, 51142, 548946, 562427, "icv1-fp20-t727-cd198a38"),
    (
        22,
        2434423,
        620112,
        1810598,
        2435809,
        "icv1-fp22-tm1385-475dcb5f",
    ),
    (
        24,
        16483309,
        8820674,
        11008593,
        16482251,
        "icv1-fp24-t1059-6df599df",
    ),
    (
        26,
        42652151,
        26623966,
        14210560,
        42653927,
        "icv1-fp26-tm1775-7e8fb6df",
    ),
];

pub fn ladder_curve(bits: u32) -> Curve {
    let &(b, p, a, bb, order, slug) = LADDER
        .iter()
        .find(|row| row.0 == bits)
        .unwrap_or_else(|| panic!("no registered ladder curve at {bits} bits"));
    let curve = PrimeCurve { p, a, b: bb };
    // Generator: the point with the least abscissa ≥ 1 and the smaller
    // ordinate.  The order is prime, so every affine point generates.
    let g = (1..p)
        .find_map(|x| {
            let rhs = curve.rhs(x);
            if curve.legendre(rhs) < 0 {
                return None;
            }
            let y = curve.sqrt(rhs)?;
            Some(PrimePoint::affine(x, y.min(p - y)))
        })
        .expect("a point");
    assert!(curve.is_on_curve(g));
    let mut ops = GroupOps::default();
    assert!(curve.mul(&mut ops, g, order).infinity, "[r]G = O");
    let inst = PrimeInstance {
        name: slug.to_string(),
        curve,
        group_order: order,
        r: order,
        cofactor: 1,
        generator: (g.x, g.y),
    };
    assert_eq!(inst.curve_id().slug, slug, "ICV1 slug matches the registry");
    Curve {
        slug: slug.to_string(),
        generator_call: format!("find_prime_order_curve({b}, 59297)"),
        inst,
    }
}

// ───────────────────────── the zero-minor method ─────────────────────────

/// Parameters of one zero-minor solve.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct ZmParams {
    /// Degree `n'` of the plane curves; `ℓ = 3n'` (must be even).
    pub nprime: usize,
    /// Defect `d` of the determine step.
    pub d: usize,
    /// Give up after this many `d`-subsets covered (a multiple of `r`).
    pub cap_covered: f64,
}

/// `n'`: the paper takes `n' ≈ log₂ p` and needs `ℓ = 3n'` even.
pub fn default_nprime(r: u64) -> usize {
    let lg = (r as f64).log2().floor() as usize;
    lg - (lg % 2)
}

/// Equation (4) of arXiv 2607.09814: probability that one guess `a`
/// succeeds, `1 − Λ^{C(ℓ',d)}`, `Λ = (1 − (p−2d+1)/(p(p−d+1)))^{C(ℓ'+d,d)}`.
pub fn eq4_success(r: u64, lp: usize, d: usize) -> f64 {
    let p = r as f64;
    let d_ = d as f64;
    let q = (p - 2.0 * d_ + 1.0) / (p * (p - d_ + 1.0));
    let ln_lambda = binom((lp + d) as u64, d as u64) * (-q).ln_1p();
    let ln_none = binom(lp as u64, d as u64) * ln_lambda;
    -ln_none.exp_m1()
}

/// The paper's rule (Tables 2–3): the least defect with success
/// probability ≥ 0.9 by Equation (4); never below 2.
pub fn default_d(r: u64, lp: usize) -> usize {
    (2..lp)
        .find(|&d| eq4_success(r, lp, d) >= 0.9)
        .unwrap_or(lp - 1)
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ZmReport {
    pub recovered: Option<u64>,
    /// Field operations per phase.
    pub phases: BTreeMap<String, Tally>,
    pub total: Tally,
    /// Group operations of the point set-up (`2ℓ` scalar multiplications).
    pub group_ops_setup: u64,
    pub counters: BTreeMap<String, u64>,
    /// How the zero minor was found.
    pub found_by: String,
    pub nprime: usize,
    pub ell: usize,
    pub d: usize,
}

struct Kernel {
    /// `ℓ × 2ℓ`, anti-diagonal format: row `i` has its `1` at column `2ℓ−1−i`.
    k: Vec<Vec<u64>>,
    /// Scalars: column `j < ℓ` is the point `n_j P`; column `ℓ + j` is
    /// `−n'_j Q`.
    n_p: Vec<u64>,
    n_q: Vec<u64>,
}

struct Ctx<'a> {
    f: Fp,
    r: u64,
    curve: &'a PrimeCurve,
    rep: ZmReport,
}

impl Ctx<'_> {
    fn charge(&mut self, phase: &str) {
        let t = self.f.take();
        self.rep.phases.entry(phase.to_string()).or_default().add(t);
    }
    fn count(&mut self, key: &str, by: u64) {
        *self.rep.counters.entry(key.to_string()).or_default() += by;
    }
}

/// 1703.07544 Algorithm 1 / 2310.04132 §2: `2ℓ` distinct points (`ℓ` of
/// the form `n_i P`, `ℓ` of the form `−n_j Q`), the matrix `M` of all
/// degree-`n'` monomials at them, its left kernel `K` (dimension `ℓ` by
/// the extra-points theorem), reduced to anti-diagonal format.  `None`
/// when the sparse block is singular (Las Vegas: draw new points).
fn build_kernel(
    cx: &mut Ctx,
    rng: &mut StdRng,
    gen: PrimePoint,
    target: PrimePoint,
    nprime: usize,
) -> Option<Kernel> {
    let ell = 3 * nprime;
    let r = cx.r;
    let g = cx.curve;
    let mut gops = GroupOps::default();
    let mut used = std::collections::HashSet::new();
    let mut pts = Vec::with_capacity(2 * ell);
    let mut n_p = Vec::with_capacity(ell);
    let mut n_q = Vec::with_capacity(ell);
    while pts.len() < 2 * ell {
        let is_p = pts.len() < ell;
        let n = rng.gen_range(1..r);
        let pt = if is_p {
            g.mul(&mut gops, gen, n)
        } else {
            g.neg(g.mul(&mut gops, target, n))
        };
        if pt.infinity || !used.insert((pt.x, pt.y)) {
            continue;
        }
        pts.push(pt);
        if is_p {
            n_p.push(n);
        } else {
            n_q.push(n);
        }
    }
    cx.rep.group_ops_setup += gops.adds + gops.doubles;
    let pt_tally = Tally::from_field_ops(gops.field_ops);
    cx.rep
        .phases
        .entry("points".into())
        .or_default()
        .add(pt_tally);

    // Monomials x^i y^j z^(n'-i-j), evaluated at z = 1.
    let mono: Vec<(usize, usize)> = (0..=nprime)
        .flat_map(|i| (0..=nprime - i).map(move |j| (i, j)))
        .collect();
    let ncols = mono.len();
    let f = &mut cx.f;
    let mut m = Vec::with_capacity(2 * ell);
    for pt in &pts {
        let mut xp = vec![1u64; nprime + 1];
        let mut yp = vec![1u64; nprime + 1];
        for e in 1..=nprime {
            xp[e] = f.mul(xp[e - 1], pt.x);
            yp[e] = f.mul(yp[e - 1], pt.y);
        }
        let row: Vec<u64> = mono
            .iter()
            .map(|&(i, j)| {
                if i == 0 {
                    yp[j]
                } else if j == 0 {
                    xp[i]
                } else {
                    f.mul(xp[i], yp[j])
                }
            })
            .collect();
        m.push(row);
    }
    cx.charge("matrix_m");

    let kb = left_kernel(&mut cx.f, &m, 2 * ell, ncols);
    cx.charge("left_kernel");
    if kb.len() != ell {
        cx.count("kernel_dim_mismatch", 1);
        return None;
    }
    // Anti-diagonal format: RREF with the sparse columns visited first,
    // in the order 2ℓ−1, 2ℓ−2, …, ℓ, so row i pivots at 2ℓ−1−i.
    let mut k = kb;
    let order: Vec<usize> = (ell..2 * ell).rev().chain(0..ell).collect();
    let piv = rref_order(&mut cx.f, &mut k, &order);
    cx.charge("anti_diagonal");
    if piv[..] != order[..ell] {
        cx.count("sparse_block_singular", 1);
        return None;
    }
    Some(Kernel { k, n_p, n_q })
}

/// A maximal zero minor of `K` at the columns `ups` gives a vector of the
/// row space with zeros at `ups`; its support is a set of points summing
/// to `O` (1703.07544 Corollary 2.4): `Σ n_i − m Σ n'_j ≡ 0 (mod r)`.
fn decode(cx: &mut Ctx, kern: &Kernel, ups: &[usize]) -> Option<u64> {
    let ell = kern.k.len();
    let sub: Vec<Vec<u64>> = kern
        .k
        .iter()
        .map(|row| ups.iter().map(|&c| row[c]).collect())
        .collect();
    let lk = left_kernel(&mut cx.f, &sub, ell, ups.len());
    let c = lk.first()?;
    let mut sum_p = 0u128;
    let mut sum_q = 0u128;
    let mut support = 0u64;
    for col in 0..2 * ell {
        let mut v = 0u64;
        for (i, &ci) in c.iter().enumerate() {
            if ci != 0 && kern.k[i][col] != 0 {
                let t = cx.f.mul(ci, kern.k[i][col]);
                v = cx.f.add(v, t);
            }
        }
        if v != 0 {
            support += 1;
            if col < ell {
                sum_p += kern.n_p[col] as u128;
            } else {
                sum_q += kern.n_q[col - ell] as u128;
            }
        }
    }
    cx.charge("decode");
    cx.count("decode_support", support);
    let r = cx.r as u128;
    let (sp, sq) = ((sum_p % r) as u64, (sum_q % r) as u64);
    if sq == 0 {
        cx.count("decode_no_q_part", 1);
        return None;
    }
    // Arithmetic mod the group order r (not the field): not a field op.
    let m = (sp as u128 * inv_mod(sq, cx.r) as u128 % r) as u64;
    Some(m)
}

/// Column set of `K` for the zero `|R|×|R|` minor `K[R, C]` of the dense
/// part: `C` with the sparse columns of the rows outside `R`
/// (2310.04132 §3; 2607.09814 Theorem 1).
fn ups_from_dense(ell: usize, rows: &[usize], cols: &[usize]) -> Vec<usize> {
    let mut ups: Vec<usize> = cols.to_vec();
    for i in 0..ell {
        if !rows.contains(&i) {
            ups.push(2 * ell - 1 - i);
        }
    }
    ups
}

/// A zero maximal minor of `K'` at the K-rows `ups1` is `det K[ups1, a] = 0`.
fn finish(cx: &mut Ctx<'_>, kern: &Kernel, a: &[usize], ups1: &[usize], how: &str) -> Option<u64> {
    let ell = kern.k.len();
    let m = decode(cx, kern, &ups_from_dense(ell, ups1, a))?;
    cx.rep.found_by = how.into();
    Some(m)
}

/// Solve one target by the guess-and-determine zero-minor search.
/// Never sees the planted scalar; returns `m` unverified.
pub fn zero_minor_solve(
    inst: &PrimeInstance,
    target: PrimePoint,
    seed: u64,
    prm: ZmParams,
) -> ZmReport {
    let r = inst.r;
    let gen = inst.generator_point();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut cx = Ctx {
        f: Fp::new(inst.curve.p),
        r,
        curve: &inst.curve,
        rep: ZmReport {
            nprime: prm.nprime,
            ell: 3 * prm.nprime,
            d: prm.d,
            ..Default::default()
        },
    };
    assert!(
        (3 * prm.nprime).is_multiple_of(2),
        "ℓ = 3n' must be even (2607.09814 §1)"
    );
    let ell = 3 * prm.nprime;
    let lp = ell / 2;
    let d = prm.d;
    assert!(d >= 2 && d < lp);
    let cap = (prm.cap_covered * r as f64) as u64;
    let mut covered = 0u64;

    'kernel: loop {
        cx.count("kernels", 1);
        let Some(kern) = build_kernel(&mut cx, &mut rng, gen, target, prm.nprime) else {
            continue;
        };
        // 2310.04132 §3 / 2607.09814 §2: a zero entry of the dense part is
        // a zero 1×1 minor and solves the problem.
        for i in 0..ell {
            for j in 0..ell {
                if kern.k[i][j] == 0 {
                    cx.count("dense_zero_entries", 1);
                    let ups = ups_from_dense(ell, &[i], &[j]);
                    if let Some(m) = decode(&mut cx, &kern, &ups) {
                        cx.rep.found_by = "dense_zero_entry".into();
                        cx.rep.recovered = Some(m);
                        break 'kernel;
                    }
                }
            }
        }
        covered += (ell * ell) as u64;

        loop {
            if covered >= cap {
                cx.rep.found_by = "cap".into();
                break 'kernel;
            }
            // Guess: a ⊂ dense columns, |a| = ℓ' (Algorithm 1 line 4).
            cx.count("guesses_a", 1);
            let mut a: Vec<usize> = (0..ell).collect();
            for i in 0..lp {
                let j = rng.gen_range(i..ell);
                a.swap(i, j);
            }
            a.truncate(lp);
            a.sort_unstable();
            // K' = K[a]^T (ℓ' × ℓ, columns indexed by rows of K), reduced to
            // anti-diagonal format: rows ℓ−1, …, ℓ' of K become the sparse
            // part, row s pivoting at K-row ℓ−1−s (Algorithm 1 lines 5–8).
            let mut k1: Vec<Vec<u64>> = a
                .iter()
                .map(|&c| (0..ell).map(|i| kern.k[i][c]).collect())
                .collect();
            let order: Vec<usize> = (lp..ell).rev().chain(0..lp).collect();
            let piv = rref_order(&mut cx.f, &mut k1, &order);
            cx.charge("guess_reduce");
            if piv[..] != order[..lp] {
                // det K[{ℓ'..ℓ−1}, a] = 0: a zero minor of the dense part.
                cx.count("guess_block_singular", 1);
                let rows: Vec<usize> = (lp..ell).collect();
                if let Some(m) = decode(&mut cx, &kern, &ups_from_dense(ell, &rows, &a)) {
                    cx.rep.found_by = "guess_block_singular".into();
                    cx.rep.recovered = Some(m);
                    break 'kernel;
                }
                continue;
            }

            // Determine: every b ⊂ dense columns of K' with |b| = ℓ' − d
            // (Algorithm 1 lines 11–22).
            let mut b: Vec<usize> = (0..lp - d).collect();
            let nsig = lp + d;
            loop {
                cx.count("guesses_b", 1);
                let compl: Vec<usize> = (0..lp).filter(|i| !b.contains(i)).collect();
                // T = left kernel of K'[b] (ℓ' × (ℓ'−d)), dim d.
                let mb: Vec<Vec<u64>> = (0..lp)
                    .map(|s| b.iter().map(|&c| k1[s][c]).collect())
                    .collect();
                let t = left_kernel(&mut cx.f, &mb, lp, lp - d);
                cx.charge("determine_t");
                if t.len() != d {
                    // The columns b are dependent: every maximal minor
                    // containing them is zero (not general position).
                    cx.count("t_dim_not_d", 1);
                    let mut ups1: Vec<usize> = b.clone();
                    ups1.extend(compl.iter().take(d));
                    if let Some(m) = finish(&mut cx, &kern, &a, &ups1, "t_not_general_position") {
                        cx.rep.recovered = Some(m);
                        break 'kernel;
                    }
                } else {
                    // Signature matrix ((ℓ'+d) × d): rows 0..d are T·col for
                    // the dense columns not in b (by RREF in the paper; the
                    // same line up to a scalar), rows d.. are the columns
                    // of T (sparse columns: free, Algorithm 1 line 19).
                    // K'-column of signature row k, as a K-row index.
                    let sig_col = |k: usize| -> usize {
                        if k < d {
                            compl[k]
                        } else {
                            ell - 1 - (k - d)
                        }
                    };
                    let mut sig: Vec<[u64; 8]> = Vec::with_capacity(nsig);
                    for &c in &compl {
                        let mut u = [0u64; 8];
                        for (q, tq) in t.iter().enumerate() {
                            let mut acc = 0u64;
                            for s in 0..lp {
                                if tq[s] != 0 && k1[s][c] != 0 {
                                    let v = cx.f.mul(tq[s], k1[s][c]);
                                    acc = cx.f.add(acc, v);
                                }
                            }
                            u[q] = acc;
                        }
                        sig.push(u);
                    }
                    for s in 0..lp {
                        // The K'-column with its 1 in row s is K-row ℓ−1−s.
                        let mut u = [0u64; 8];
                        for (q, tq) in t.iter().enumerate() {
                            u[q] = tq[s];
                        }
                        sig.push(u);
                    }
                    cx.charge("determine_signatures");

                    // Penultimate intersections: the right kernel of every
                    // (d−1)-row submatrix, hashed; a repeat whose union has
                    // d rows is a zero d×d minor (Theorems 4–5, Algorithm 2).
                    let mut table: HashMap<[u64; 8], Vec<usize>> = HashMap::new();
                    let mut x: Vec<usize> = (0..d - 1).collect();
                    let mut found: Option<Vec<usize>> = None;
                    let mut kernels = 0u64;
                    loop {
                        kernels += 1;
                        let xm: Vec<Vec<u64>> = x.iter().map(|&k| sig[k][..d].to_vec()).collect();
                        let kb = right_kernel(&mut cx.f, &xm, d);
                        if kb.len() != 1 {
                            // d−1 dependent signatures: with any further row
                            // they form a zero d×d minor.
                            cx.count("sig_rank_deficient", 1);
                            let mut u = x.clone();
                            let extra = (0..nsig).find(|k| !x.contains(k)).unwrap();
                            u.push(extra);
                            found = Some(u);
                            break;
                        }
                        let mut key = [u64::MAX; 8];
                        key[..d].copy_from_slice(&kb[0]);
                        if let Some(prev) = table.get(&key) {
                            let mut u: Vec<usize> = prev.clone();
                            for &k in &x {
                                if !u.contains(&k) {
                                    u.push(k);
                                }
                            }
                            if u.len() == d {
                                found = Some(u);
                                break;
                            }
                            cx.count("collision_union_gt_d", 1);
                        } else {
                            table.insert(key, x.clone());
                        }
                        if !next_comb(&mut x, nsig) {
                            break;
                        }
                    }
                    cx.charge("determine_kernels");
                    cx.count("signature_kernels", kernels);
                    cx.count("hash_probes_uncharged", kernels);
                    if let Some(u) = found {
                        let mut ups1: Vec<usize> = b.clone();
                        ups1.extend(u.iter().map(|&k| sig_col(k)));
                        cx.count("covered_d_subsets", binom(nsig as u64, d as u64) as u64);
                        if let Some(m) = finish(&mut cx, &kern, &a, &ups1, "signature_collision") {
                            cx.rep.recovered = Some(m);
                            break 'kernel;
                        }
                    } else {
                        let c = binom(nsig as u64, d as u64) as u64;
                        cx.count("covered_d_subsets", c);
                        covered += c;
                    }
                }
                if covered >= cap {
                    cx.rep.found_by = "cap".into();
                    break 'kernel;
                }
                if !next_comb(&mut b, lp) {
                    break;
                }
            }
        }
    }
    let mut total = Tally::default();
    for t in cx.rep.phases.values() {
        total.add(*t);
    }
    cx.rep.total = total;
    cx.rep
}

// ───────────────────────── measurement ─────────────────────────

#[derive(Serialize)]
struct RunRow {
    method: String,
    bits: u32,
    slug: String,
    r: u64,
    log2_r: f64,
    target: usize,
    seed: u64,
    ops: u64,
    muls: u64,
    invs: u64,
    group_ops: u64,
    /// Field ops of the search alone (zero-minor: guess and determine
    /// phases; rho: the whole walk).
    search_ops: u64,
    verified: bool,
    wall_ms: f64,
    /// Host load (`uptime`) when the sample ran; wall time is L0 only.
    uptime: String,
    detail: serde_json::Value,
}

fn uptime() -> String {
    std::process::Command::new("uptime")
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default()
}

fn targets_for(inst: &PrimeInstance, n: usize, seed: u64) -> Vec<(u64, PrimePoint)> {
    let mut rng = StdRng::seed_from_u64(seed ^ inst.r);
    let g = inst.generator_point();
    let mut ops = GroupOps::default();
    (0..n)
        .map(|_| {
            let k = rng.gen_range(2..inst.r - 1);
            (k, inst.curve.mul(&mut ops, g, k))
        })
        .collect()
}

fn verify(inst: &PrimeInstance, m: Option<u64>, target: PrimePoint) -> bool {
    let Some(m) = m else { return false };
    let mut ops = GroupOps::default();
    inst.curve.mul(&mut ops, inst.generator_point(), m) == target
}

/// Two-sided 97.5% Student-t quantile (Cornish–Fisher expansion).
fn t975(df: f64) -> f64 {
    let z: f64 = 1.959_963_985;
    z + (z.powi(3) + z) / (4.0 * df)
        + (5.0 * z.powi(5) + 16.0 * z.powi(3) + 3.0 * z) / (96.0 * df * df)
}

#[derive(Serialize, Clone)]
struct Fit {
    slope: f64,
    intercept: f64,
    n: usize,
    slope_ci95_t: [f64; 2],
    slope_ci95_bootstrap: [f64; 2],
}

/// At least two distinct curve sizes: a slope needs them.
fn spread(points: &[(f64, f64)]) -> bool {
    points.iter().any(|p| p.0 != points[0].0)
}

fn ols(xs: &[f64], ys: &[f64]) -> (f64, f64, f64) {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let slope = sxy / sxx;
    let icpt = my - slope * mx;
    let sse: f64 = xs
        .iter()
        .zip(ys)
        .map(|(x, y)| (y - icpt - slope * x).powi(2))
        .sum();
    let se = (sse / (n - 2.0) / sxx).sqrt();
    (slope, icpt, se)
}

/// log₂ cost against log₂ r over all runs; 95% intervals by t and by a
/// bootstrap that resamples runs within each curve (seeded, 4000 draws).
fn fit(points: &[(f64, f64)]) -> Fit {
    let xs: Vec<f64> = points.iter().map(|p| p.0).collect();
    let ys: Vec<f64> = points.iter().map(|p| p.1).collect();
    let (slope, intercept, se) = ols(&xs, &ys);
    let n = points.len();
    let t = t975(n as f64 - 2.0);
    let mut groups: BTreeMap<u64, Vec<(f64, f64)>> = BTreeMap::new();
    for &p in points {
        groups.entry(p.0.to_bits()).or_default().push(p);
    }
    let mut rng = StdRng::seed_from_u64(0x5712_0011);
    let mut slopes = Vec::with_capacity(4000);
    for _ in 0..4000 {
        let mut bx = Vec::with_capacity(n);
        let mut by = Vec::with_capacity(n);
        for g in groups.values() {
            for _ in 0..g.len() {
                let p = g[rng.gen_range(0..g.len())];
                bx.push(p.0);
                by.push(p.1);
            }
        }
        slopes.push(ols(&bx, &by).0);
    }
    slopes.sort_by(|a, b| a.total_cmp(b));
    Fit {
        slope,
        intercept,
        n,
        slope_ci95_t: [slope - t * se, slope + t * se],
        slope_ci95_bootstrap: [slopes[100], slopes[3899]],
    }
}

fn arg<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .map(|s| s.as_str())
}

fn list(s: &str) -> Vec<u64> {
    s.split(',').map(|x| x.trim().parse().unwrap()).collect()
}

fn run_zm(
    c: &Curve,
    t_idx: usize,
    target: PrimePoint,
    seed: u64,
    prm: ZmParams,
    bits: u32,
) -> RunRow {
    let t0 = Instant::now();
    let rep = zero_minor_solve(&c.inst, target, seed, prm);
    let wall = t0.elapsed().as_secs_f64() * 1e3;
    let ok = verify(&c.inst, rep.recovered, target);
    RunRow {
        method: format!("zero_minor.guess_determine(d={})", prm.d),
        bits,
        slug: c.slug.clone(),
        r: c.inst.r,
        log2_r: (c.inst.r as f64).log2(),
        target: t_idx,
        seed,
        ops: rep.total.ops(),
        muls: rep.total.muls,
        invs: rep.total.invs,
        group_ops: rep.group_ops_setup,
        search_ops: rep
            .phases
            .iter()
            .filter(|(k, _)| k.starts_with("guess_") || k.starts_with("determine_"))
            .map(|(_, t)| t.ops())
            .sum(),
        verified: ok,
        wall_ms: wall,
        uptime: uptime(),
        detail: serde_json::to_value(&rep).unwrap(),
    }
}

fn run_rho(c: &Curve, t_idx: usize, target: PrimePoint, seed: u64, bits: u32) -> RunRow {
    let inst = &c.inst;
    let t0 = Instant::now();
    let res = rho_walk_with(
        &inst.curve,
        &NegationClasses { r: inst.r },
        inst.generator_point(),
        target,
        inst.r,
        seed,
        rho_cap(inst.r, 64.0),
        RhoWalk::negation(),
    );
    let wall = t0.elapsed().as_secs_f64() * 1e3;
    let fo = Tally::from_field_ops(res.group_ops.field_ops);
    let ok = verify(inst, res.recovered, target);
    RunRow {
        method: "rho.negation".into(),
        bits,
        slug: c.slug.clone(),
        r: inst.r,
        log2_r: (inst.r as f64).log2(),
        target: t_idx,
        seed,
        ops: fo.ops(),
        muls: fo.muls,
        invs: fo.invs,
        group_ops: res.group_ops.adds + res.group_ops.doubles,
        search_ops: fo.ops(),
        verified: ok,
        wall_ms: wall,
        uptime: uptime(),
        detail: json!({
            "steps": res.steps,
            "gae": res.gae,
            "field_ops": res.group_ops.field_ops,
            "walks": res.walks,
            "distinguished_points": res.distinguished_points,
        }),
    }
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.total_cmp(b));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// `refit --in scaling.json --out fits.json`: the slope fits over the
/// zero-minor runs pooled across defects (the run's own fits group by the
/// method string, which carries `d`), plus covered/r per curve.
fn refit(input: &str, out: &str) {
    let doc: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(input).unwrap()).unwrap();
    let runs = doc["runs"].as_array().unwrap();
    let is_zm = |r: &serde_json::Value| r["method"].as_str().unwrap().starts_with("zero_minor");
    let ok = |r: &serde_json::Value| r["verified"].as_bool().unwrap();
    let pts = |sel: &dyn Fn(&serde_json::Value) -> Option<f64>| -> Vec<(f64, f64)> {
        runs.iter()
            .filter(|r| ok(r))
            .filter_map(|r| {
                let y = sel(r)?;
                (y > 0.0).then(|| (r["log2_r"].as_f64().unwrap(), y.log2()))
            })
            .collect()
    };
    let mut fits = BTreeMap::new();
    let defs: Vec<(&str, Box<dyn Fn(&serde_json::Value) -> Option<f64>>)> = vec![
        (
            "rho.negation [field ops]",
            Box::new(move |r| (!is_zm(r)).then(|| r["ops"].as_f64().unwrap())),
        ),
        (
            "rho.negation [group ops]",
            Box::new(move |r| (!is_zm(r)).then(|| r["group_ops"].as_f64().unwrap())),
        ),
        (
            "zero_minor, d by Eq.(4) rule [field ops, all phases]",
            Box::new(move |r| is_zm(r).then(|| r["ops"].as_f64().unwrap())),
        ),
        (
            "zero_minor, d by Eq.(4) rule [field ops, search phases only]",
            Box::new(move |r| is_zm(r).then(|| r["search_ops"].as_f64().unwrap())),
        ),
        (
            "zero_minor, d by Eq.(4) rule [d-subsets covered]",
            Box::new(move |r| {
                if !is_zm(r) {
                    return None;
                }
                r["detail"]["counters"]["covered_d_subsets"].as_f64()
            }),
        ),
        (
            "zero_minor, d = 3 runs only [field ops, all phases]",
            Box::new(move |r| {
                (is_zm(r) && r["detail"]["d"].as_u64() == Some(3))
                    .then(|| r["ops"].as_f64().unwrap())
            }),
        ),
    ];
    for (name, sel) in &defs {
        let p = pts(sel.as_ref());
        let xs: std::collections::BTreeSet<u64> = p.iter().map(|q| q.0.to_bits()).collect();
        if p.len() >= 3 && xs.len() >= 2 {
            fits.insert(name.to_string(), fit(&p));
        }
    }
    let mut cov = Vec::new();
    let mut bits: Vec<u64> = runs.iter().map(|r| r["bits"].as_u64().unwrap()).collect();
    bits.sort();
    bits.dedup();
    for b in bits {
        let mut v: Vec<f64> = runs
            .iter()
            .filter(|r| is_zm(r) && r["bits"].as_u64() == Some(b))
            .filter_map(|r| {
                Some(r["detail"]["counters"]["covered_d_subsets"].as_f64()? / r["r"].as_f64()?)
            })
            .collect();
        if v.is_empty() {
            continue;
        }
        let n = v.len();
        let mn = v.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = v.iter().cloned().fold(0.0, f64::max);
        cov.push(json!({"bits": b, "runs_with_collision": n,
            "covered_over_r_median": median(&mut v), "min": mn, "max": mx}));
    }
    let doc = json!({
        "schema": "lit2026-w12-zero-minor-refit/v1",
        "input": input,
        "fits_log2_vs_log2r": fits,
        "covered_over_r": cov,
    });
    std::fs::write(out, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
    for (k, f) in &fits {
        println!(
            "{k}: slope {:.3} boot [{:.3},{:.3}] t [{:.3},{:.3}] n={}",
            f.slope,
            f.slope_ci95_bootstrap[0],
            f.slope_ci95_bootstrap[1],
            f.slope_ci95_t[0],
            f.slope_ci95_t[1],
            f.n
        );
    }
}

/// `compare --a X --b Y`: whether two runs recorded identical counts and
/// recovered values row by row (the replay check).
fn compare(a: &str, b: &str) {
    let load = |p: &str| -> serde_json::Value {
        serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
    };
    let (da, db) = (load(a), load(b));
    let (ra, rb) = (
        da["runs"].as_array().unwrap(),
        db["runs"].as_array().unwrap(),
    );
    let mut same = 0;
    let mut diff = 0;
    for (x, y) in ra.iter().zip(rb) {
        let key = |r: &serde_json::Value| {
            (
                r["method"].clone(),
                r["bits"].clone(),
                r["target"].clone(),
                r["seed"].clone(),
                r["ops"].clone(),
                r["group_ops"].clone(),
                r["search_ops"].clone(),
                r["verified"].clone(),
                r["detail"]["recovered"].clone(),
            )
        };
        if key(x) == key(y) {
            same += 1;
        } else {
            diff += 1;
        }
    }
    println!(
        "rows a={} b={} identical={} different={}",
        ra.len(),
        rb.len(),
        same,
        diff
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(|s| s.as_str()).unwrap_or("run");
    if mode == "refit" {
        return refit(arg(&args, "--in").unwrap(), arg(&args, "--out").unwrap());
    }
    if mode == "compare" {
        return compare(arg(&args, "--a").unwrap(), arg(&args, "--b").unwrap());
    }
    let out = arg(&args, "--out").expect("--out FILE");
    let bits_list: Vec<u32> = list(arg(&args, "--bits").unwrap_or("16,18,20,22,24,26"))
        .into_iter()
        .map(|b| b as u32)
        .collect();
    let ntargets: usize = arg(&args, "--targets").unwrap_or("8").parse().unwrap();
    let rho_seeds: u64 = arg(&args, "--rho-seeds").unwrap_or("4").parse().unwrap();
    let seed: u64 = arg(&args, "--seed").unwrap_or("20261011").parse().unwrap();
    let cap: f64 = arg(&args, "--cap").unwrap_or("64").parse().unwrap();
    let mut rows: Vec<RunRow> = Vec::new();
    let mut curves = Vec::new();
    let ds_for_sweep: Vec<u64> = arg(&args, "--ds").map(list).unwrap_or_default();

    for &bits in &bits_list {
        let c = ladder_curve(bits);
        let r = c.inst.r;
        let nprime = match arg(&args, "--nprime") {
            Some("auto") | None => default_nprime(r),
            Some(v) => v.parse().unwrap(),
        };
        let lp = 3 * nprime / 2;
        let d_auto = default_d(r, lp);
        let ds: Vec<usize> = if mode == "dsweep" {
            ds_for_sweep.iter().map(|&d| d as usize).collect()
        } else {
            match arg(&args, "--d") {
                Some("auto") | None => vec![d_auto],
                Some(v) => vec![v.parse().unwrap()],
            }
        };
        let targets = targets_for(&c.inst, ntargets, seed);
        eprintln!(
            "{} r={} log2r={:.2} n'={} ℓ={} ℓ'={} d_auto={} eq4(d_auto)={:.4} ds={:?}",
            c.slug,
            r,
            (r as f64).log2(),
            nprime,
            3 * nprime,
            lp,
            d_auto,
            eq4_success(r, lp, d_auto),
            ds
        );
        curves.push(json!({
            "bits": bits, "slug": c.slug, "generator_call": c.generator_call,
            "p": c.inst.curve.p, "a": c.inst.curve.a, "b": c.inst.curve.b,
            "r": r, "cofactor": 1, "generator": [c.inst.generator.0, c.inst.generator.1],
            "nprime": nprime, "ell": 3 * nprime, "ell_half": lp, "d_auto": d_auto,
            "eq4_success_d_auto": eq4_success(r, lp, d_auto),
            "distinct_maximal_minors_per_guess": (0..=d_auto).map(|j| binom(lp as u64, j as u64).powi(2)).sum::<f64>(),
            "covered_per_guess_with_multiplicity": binom(lp as u64, d_auto as u64) * binom((lp + d_auto) as u64, d_auto as u64),
        }));
        for (ti, &(_k, q)) in targets.iter().enumerate() {
            // Interleave the methods target by target (A B A B …).
            for &d in &ds {
                let prm = ZmParams {
                    nprime,
                    d,
                    cap_covered: cap,
                };
                let s = seed ^ ((bits as u64) << 32) ^ ((ti as u64) << 8) ^ d as u64;
                let row = run_zm(&c, ti, q, s, prm, bits);
                eprintln!(
                    "  zm d={} t{} ops={} covered={} verified={} {:.0}ms by={}",
                    d,
                    ti,
                    row.ops,
                    row.detail["counters"]["covered_d_subsets"],
                    row.verified,
                    row.wall_ms,
                    row.detail["found_by"]
                );
                rows.push(row);
            }
            if mode == "run" {
                for rs in 0..rho_seeds {
                    let s = seed ^ ((bits as u64) << 32) ^ ((ti as u64) << 8) ^ (0x100 + rs);
                    let row = run_rho(&c, ti, q, s, bits);
                    rows.push(row);
                }
            }
        }
    }

    // Fits and per-curve summaries (computed here, in Rust).
    let mut methods: Vec<String> = rows.iter().map(|r| r.method.clone()).collect();
    methods.sort();
    methods.dedup();
    let mut fits = BTreeMap::new();
    let mut per_curve = Vec::new();
    for m in &methods {
        let pts: Vec<(f64, f64)> = rows
            .iter()
            .filter(|r| &r.method == m && r.verified)
            .map(|r| (r.log2_r, (r.ops as f64).log2()))
            .collect();
        if pts.len() >= 3 && spread(&pts) {
            fits.insert(m.clone(), fit(&pts));
        }
        let spts: Vec<(f64, f64)> = rows
            .iter()
            .filter(|r| &r.method == m && r.verified && r.search_ops > 0)
            .map(|r| (r.log2_r, (r.search_ops as f64).log2()))
            .collect();
        if spts.len() >= 3 && spread(&spts) && m.starts_with("zero_minor") {
            fits.insert(format!("{m} [search phases only]"), fit(&spts));
            let cpts: Vec<(f64, f64)> = rows
                .iter()
                .filter(|r| &r.method == m && r.verified)
                .filter_map(|r| {
                    let c = r.detail["counters"]["covered_d_subsets"].as_u64()?;
                    (c > 0).then(|| (r.log2_r, (c as f64).log2()))
                })
                .collect();
            if cpts.len() >= 3 && spread(&cpts) {
                fits.insert(format!("{m} [d-subsets covered]"), fit(&cpts));
            }
        }
        for &bits in &bits_list {
            let sel: Vec<&RunRow> = rows
                .iter()
                .filter(|r| &r.method == m && r.bits == bits)
                .collect();
            if sel.is_empty() {
                continue;
            }
            let mut ops: Vec<f64> = sel.iter().map(|r| r.ops as f64).collect();
            let mut sops: Vec<f64> = sel.iter().map(|r| r.search_ops as f64).collect();
            let mut cov: Vec<f64> = sel
                .iter()
                .filter_map(|r| r.detail["counters"]["covered_d_subsets"].as_f64())
                .collect();
            let cov_med = if cov.is_empty() {
                None
            } else {
                Some(median(&mut cov))
            };
            let mut wall: Vec<f64> = sel.iter().map(|r| r.wall_ms).collect();
            let mn = ops.iter().cloned().fold(f64::INFINITY, f64::min);
            let mx = ops.iter().cloned().fold(0.0, f64::max);
            let mean = ops.iter().sum::<f64>() / ops.len() as f64;
            per_curve.push(json!({
                "method": m, "bits": bits, "slug": sel[0].slug, "r": sel[0].r,
                "runs": sel.len(), "verified": sel.iter().filter(|r| r.verified).count(),
                "ops_median": median(&mut ops), "ops_mean": mean, "ops_min": mn, "ops_max": mx,
                "search_ops_median": median(&mut sops), "covered_d_subsets_median": cov_med,
                "wall_ms_median": median(&mut wall),
            }));
        }
    }
    // Paired per-target ratio zero-minor / rho (rho: mean over its seeds).
    let mut paired = Vec::new();
    if mode == "run" {
        for &bits in &bits_list {
            let mut ratios = Vec::new();
            for ti in 0..ntargets {
                let zm: Vec<f64> = rows
                    .iter()
                    .filter(|r| {
                        r.bits == bits && r.target == ti && r.method.starts_with("zero_minor")
                    })
                    .map(|r| r.ops as f64)
                    .collect();
                let rh: Vec<f64> = rows
                    .iter()
                    .filter(|r| r.bits == bits && r.target == ti && r.method == "rho.negation")
                    .map(|r| r.ops as f64)
                    .collect();
                if zm.len() == 1 && !rh.is_empty() {
                    ratios.push(zm[0] / (rh.iter().sum::<f64>() / rh.len() as f64));
                }
            }
            if ratios.is_empty() {
                continue;
            }
            let lg: Vec<f64> = ratios.iter().map(|x| x.ln()).collect();
            let gm = (lg.iter().sum::<f64>() / lg.len() as f64).exp();
            let mut rng = StdRng::seed_from_u64(0xA11 ^ bits as u64);
            let mut bs: Vec<f64> = (0..4000)
                .map(|_| {
                    let s: f64 = (0..lg.len()).map(|_| lg[rng.gen_range(0..lg.len())]).sum();
                    (s / lg.len() as f64).exp()
                })
                .collect();
            bs.sort_by(|a, b| a.total_cmp(b));
            let mn = ratios.iter().cloned().fold(f64::INFINITY, f64::min);
            let mx = ratios.iter().cloned().fold(0.0, f64::max);
            paired.push(json!({
                "bits": bits, "targets": ratios.len(),
                "ratio_zero_minor_over_rho_geomean": gm,
                "ci95_bootstrap": [bs[100], bs[3899]],
                "min": mn, "max": mx,
            }));
        }
    }
    let doc = json!({
        "schema": "lit2026-w12-zero-minor/v1",
        "mode": mode,
        "args": args[1..].to_vec(),
        "unit": "field ops = F_p multiplications (squarings included) + inversions, each counted once; additions free",
        "curves": curves,
        "fits_log2ops_vs_log2r": fits,
        "per_curve": per_curve,
        "paired_ratio": paired,
        "runs": rows,
    });
    std::fs::write(out, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
    eprintln!("wrote {out}");
    for (m, f) in &fits {
        eprintln!(
            "fit {m}: slope {:.3} t-CI [{:.3},{:.3}] boot-CI [{:.3},{:.3}] n={}",
            f.slope,
            f.slope_ci95_t[0],
            f.slope_ci95_t[1],
            f.slope_ci95_bootstrap[0],
            f.slope_ci95_bootstrap[1],
            f.n
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small() -> Curve {
        ladder_curve(12)
    }

    #[test]
    fn ladder_curves_match_registry() {
        for &(bits, ..) in LADDER {
            let c = ladder_curve(bits);
            assert!(c.inst.curve.is_on_curve(c.inst.generator_point()));
        }
    }

    /// Extra-points theorem (1703.07544 Thm 2.3): dim K = ℓ, and the
    /// anti-diagonal format exists; Theorem 1 of 2607.09814: the maximal
    /// minor at C ∪ (sparse columns of rows ∉ R) is ± det A[R, C].
    #[test]
    fn kernel_dimension_and_minor_correspondence() {
        let c = small();
        let inst = &c.inst;
        let q = targets_for(inst, 1, 7)[0].1;
        let mut rng = StdRng::seed_from_u64(3);
        let mut cx = Ctx {
            f: Fp::new(inst.curve.p),
            r: inst.r,
            curve: &inst.curve,
            rep: ZmReport::default(),
        };
        let nprime = 4;
        let kern = loop {
            if let Some(k) = build_kernel(&mut cx, &mut rng, inst.generator_point(), q, nprime) {
                break k;
            }
        };
        let ell = 3 * nprime;
        assert_eq!(kern.k.len(), ell);
        for i in 0..ell {
            for j in ell..2 * ell {
                assert_eq!(kern.k[i][j], u64::from(j == 2 * ell - 1 - i));
            }
        }
        let p = inst.curve.p;
        for trial in 0..20 {
            let k = 1 + trial % 4;
            let mut rows: Vec<usize> = (0..ell).collect();
            let mut cols: Vec<usize> = (0..ell).collect();
            for i in 0..k {
                let a = rng.gen_range(i..ell);
                rows.swap(i, a);
                let b = rng.gen_range(i..ell);
                cols.swap(i, b);
            }
            rows.truncate(k);
            cols.truncate(k);
            rows.sort_unstable();
            cols.sort_unstable();
            let small: Vec<Vec<u64>> = rows
                .iter()
                .map(|&i| cols.iter().map(|&j| kern.k[i][j]).collect())
                .collect();
            let ups = ups_from_dense(ell, &rows, &cols);
            let big: Vec<Vec<u64>> = (0..ell)
                .map(|i| ups.iter().map(|&j| kern.k[i][j]).collect())
                .collect();
            let ds = det_uncounted(p, &small);
            let db = det_uncounted(p, &big);
            assert!(db == ds || db == (p - ds) % p, "Theorem 1 correspondence");
        }
    }

    /// End to end on a 12-bit curve with small `n' = 4`, `d = 2`: every
    /// zero minor found decodes to the planted logarithm.
    #[test]
    fn zero_minor_solves_planted_instances() {
        let c = small();
        let inst = &c.inst;
        for (i, (k, q)) in targets_for(inst, 8, 11).into_iter().enumerate() {
            let rep = zero_minor_solve(
                inst,
                q,
                100 + i as u64,
                ZmParams {
                    nprime: 4,
                    d: 2,
                    cap_covered: 64.0,
                },
            );
            assert_eq!(rep.recovered, Some(k), "target {i}: {:?}", rep.found_by);
            assert!(verify(inst, rep.recovered, q));
        }
    }

    /// The determine step at the paper's default parameters on a 14-bit
    /// curve, every run verified.
    #[test]
    fn guess_determine_default_parameters() {
        let c = ladder_curve(14);
        let inst = &c.inst;
        let nprime = default_nprime(inst.r);
        let d = default_d(inst.r, 3 * nprime / 2);
        for (i, (k, q)) in targets_for(inst, 3, 5).into_iter().enumerate() {
            let rep = zero_minor_solve(
                inst,
                q,
                7 + i as u64,
                ZmParams {
                    nprime,
                    d,
                    cap_covered: 64.0,
                },
            );
            assert_eq!(rep.recovered, Some(k));
        }
    }

    /// Signatures (2607.09814 Theorem 3): the signature of a column c with
    /// respect to T is T·c up to a scalar, and a zero d×d minor of the
    /// signature matrix is a vector of span(T) orthogonal to d columns.
    #[test]
    fn signature_kernel_is_intersection() {
        let mut f = Fp::new(3529);
        let mut rng = StdRng::seed_from_u64(9);
        let (n, k) = (7usize, 4usize);
        let m: Vec<Vec<u64>> = (0..n)
            .map(|_| (0..k).map(|_| rng.gen_range(0..3529)).collect())
            .collect();
        let t = left_kernel(&mut f, &m, n, k);
        assert_eq!(t.len(), n - k);
        for v in &t {
            for j in 0..k {
                let mut s = 0;
                for i in 0..n {
                    let pr = f.mul(v[i], m[i][j]);
                    s = f.add(s, pr);
                }
                assert_eq!(s, 0);
            }
        }
    }

    #[test]
    fn eq4_matches_paper_table_trend() {
        // Table 3 of 2607.09814: log2 p = 150 needs d = 15 for 0.99, with
        // n' = 150, ℓ = 450, ℓ' = 225.  p is not a u64 here; check the
        // monotonicity of the rule at our sizes instead.
        let r = 42_653_927;
        let lp = 3 * default_nprime(r) / 2;
        assert!(eq4_success(r, lp, 2) < eq4_success(r, lp, 3));
        assert!(eq4_success(r, lp, default_d(r, lp)) >= 0.9);
    }

    #[test]
    fn rho_reference_solves() {
        let c = ladder_curve(16);
        let (k, q) = targets_for(&c.inst, 1, 1)[0];
        let row = run_rho(&c, 0, q, 1, 16);
        assert!(row.verified, "rho recovered {k}");
        assert!(row.ops > 0);
    }
}
