//! **The `D₄`-symmetrised `S₅` as four weighted equations in four unknowns
//! (lead 5b of `docs/ic/factor-base-encyclopedia` §8).**
//!
//! Companion to `research/d4_symmetrised_s5_20261015/`.  The curve, the
//! base and the invariant form of `S₅` are those of
//! [`super::gaudry_quartic_d4`] (candidate C8, 2026-10-13): a curve over
//! `F_{p⁴}` with a rational 2-torsion point `T = (x₀, 0)`, `f'(x₀) = c²`,
//! FGHR's coordinate `Y = (x − x₀ − c)/(x − x₀ + c)` with
//! `Y(P + T) = −Y(P)`, and `H(Y₁, …, Y₄; Y_R) = S₅ · Π(1 − Y)⁸` written on
//! the 105 monomials `σ₁^a σ₂^b σ₃^c π^d` of the Weyl group
//! `D₄ = (Z/2)³ ⋊ S₄` (`σ_k = e_k(Y²)`, `π = Y₁Y₂Y₃Y₄`, `2(a+b+c) + d ≤ 8`).
//!
//! C8 solved that system as five equations in `(σ₁, σ₂, σ₃, σ₄, π)`
//! (degrees `4, 4, 4, 4` and `π² − σ₄`), unweighted, Macaulay degree 14.
//! Here the same four Weil components are kept in four unknowns
//! `(σ₁, σ₂, σ₃, π)` with **weights `(2, 2, 2, 1)`**: every monomial of `H`
//! has weighted degree `≤ 8`, so the weighted Bézout number is
//! `8⁴ / (2·2·2·1) = 512` and the weighted Macaulay bound is
//! `Σ dᵢ − Σ wⱼ + 1 = 32 − 7 + 1 = 26`.  The Macaulay matrix is graded by
//! weighted degree and ordered by weighted grevlex; everything after the
//! echelon (normal forms, multiplication matrix of `π`, characteristic
//! polynomial, roots, eigenvectors) is [`super::gaudry_quartic`]'s counted
//! code.  `π` has weight 1, so its multiplication matrix closes at the
//! bound; a degenerate eigenspace needs a weight-2 unknown, whose
//! multiplication matrix needs degree `bound + 1`: the solve then reruns
//! there, both runs charged.
//!
//! [`lift_back`] recovers the decompositions from the invariants and
//! charges every multiplication it makes (the split into `Y_i²`, square
//! roots in `F_p`, the points in `F_{p⁴}`, and the group check of all eight
//! even `T`-translates).  [`brute_force_decompositions`] enumerates every
//! four-subset of the base with every sign: the ground truth on tiny fields.

use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use super::gaudry_cubic::{am, mm, sm, split_eigenspace, PolyRing};
use super::gaudry_quartic::{
    charpoly_u16, echelon_by_lead, eigenvectors_u16, split_quartic, Curve4, Fp4, NormalForms4, Pt4,
    QuarticSolve, E4,
};
use super::gaudry_quartic_d4::{D4Instance, SymmetrisedD4};
use super::residual_walk::{inv_mod, pow_mod};

/// Weights of `(σ₁, σ₂, σ₃, π)`.
pub const WEIGHTS: [u8; 4] = [2, 2, 2, 1];
/// Weighted degree of every Weil component of `H`.
pub const EQUATION_WDEG: u16 = 8;
/// `Σ dᵢ − Σ wⱼ + 1` for four equations of weighted degree 8.
pub const MACAULAY_BOUND: u16 = 4 * 8 - (2 + 2 + 2 + 1) + 1;
/// The unknown whose multiplication matrix is diagonalised (`π`).
const EIG_VAR: usize = 3;

pub fn wdeg(m: &[u8; 4]) -> u16 {
    m.iter()
        .zip(WEIGHTS)
        .map(|(&e, w)| e as u16 * w as u16)
        .sum()
}

/// Weighted grevlex: weighted degree first, then the reverse-lexicographic
/// tie-break of `gaudry_quartic`'s grevlex.
pub fn wcmp(a: &[u8; 4], b: &[u8; 4]) -> Ordering {
    let (da, db) = (wdeg(a), wdeg(b));
    if da != db {
        return da.cmp(&db);
    }
    for i in (0..4).rev() {
        if a[i] != b[i] {
            return b[i].cmp(&a[i]);
        }
    }
    Ordering::Equal
}

/// Monomials of weighted degree `≤ d`, descending in weighted grevlex.
pub fn weighted_monomials_desc(d: u16) -> Vec<[u8; 4]> {
    let mut v = Vec::new();
    for a in 0..=d / 2 {
        for b in 0..=(d - 2 * a) / 2 {
            for c in 0..=(d - 2 * a - 2 * b) / 2 {
                for e in 0..=(d - 2 * a - 2 * b - 2 * c) {
                    v.push([a as u8, b as u8, c as u8, e as u8]);
                }
            }
        }
    }
    v.sort_by(|x, y| wcmp(y, x));
    v
}

fn leading(f: &HashMap<[u8; 4], u64>) -> Option<[u8; 4]> {
    f.iter()
        .filter(|(_, &c)| c != 0)
        .map(|(m, _)| *m)
        .max_by(wcmp)
}

/// Make the leading monomials pairwise distinct (as `gaudry_quartic`'s
/// `triangularise_leads4`, in the weighted order).
fn triangularise(
    comps: &[HashMap<[u8; 4], u64>],
    p: u64,
    muls: &mut u64,
) -> (Vec<HashMap<[u8; 4], u64>>, Vec<Option<[u8; 4]>>) {
    let mut fs: Vec<HashMap<[u8; 4], u64>> = comps.to_vec();
    for _round in 0..16 {
        fs.sort_by(|a, b| match (leading(a), leading(b)) {
            (Some(x), Some(y)) => wcmp(&y, &x),
            (Some(_), None) => Ordering::Less,
            (None, Some(_)) => Ordering::Greater,
            (None, None) => Ordering::Equal,
        });
        let mut changed = false;
        for i in 1..fs.len() {
            let (Some(li), Some(lj)) = (leading(&fs[i]), leading(&fs[i - 1])) else {
                continue;
            };
            if li == lj {
                let k = mm(fs[i][&li], inv_mod(fs[i - 1][&lj], p), p);
                let prev = fs[i - 1].clone();
                for (&m, &c) in &prev {
                    let e = fs[i].entry(m).or_insert(0);
                    *e = sm(*e, mm(k, c, p), p);
                    *muls += 1;
                }
                fs[i].retain(|_, c| *c != 0);
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let leads = fs.iter().map(leading).collect();
    (fs, leads)
}

fn eval(comp: &HashMap<[u8; 4], u64>, e: &[u64; 4], p: u64, muls: &mut u64) -> u64 {
    let mut acc = 0u64;
    for (m, &c) in comp {
        let mut t = c;
        for k in 0..4 {
            for _ in 0..m[k] {
                t = mm(t, e[k], p);
                *muls += 1;
            }
        }
        acc = am(acc, t, p);
    }
    acc
}

/// The four Weil components of `H(·; Y_R)` over `F_p` in `(σ₁, σ₂, σ₃, π)`
/// (the multiplications are charged to `f`).
pub fn weighted_system(pre: &SymmetrisedD4, f: &Fp4, y_r: &E4) -> Vec<HashMap<[u8; 4], u64>> {
    let mut pows = [E4::ONE; 9];
    for k in 1..9 {
        pows[k] = f.mul(&pows[k - 1], y_r);
    }
    let mut out: Vec<HashMap<[u8; 4], u64>> = vec![HashMap::new(); 4];
    for (m, c) in pre.monos.iter().zip(&pre.coef) {
        let v = (0..9).fold(E4::ZERO, |acc, j| f.add(&acc, &f.mul(&c[j], &pows[j])));
        for (i, comp) in out.iter_mut().enumerate() {
            if v.0[i] != 0 {
                comp.insert(*m, v.0[i]);
            }
        }
    }
    out
}

/// Matrix shape at weighted Macaulay degree `d` for four equations of
/// weighted degree 8, with Macaulay's row selection on generic leads
/// (derivation aid; the solve records its own shape).
pub fn columns_at(d: u16) -> usize {
    weighted_monomials_desc(d).len()
}

enum Inner {
    Solved(Vec<[u64; 4]>),
    Unsolved(String),
    NeedHigher,
}

#[allow(clippy::too_many_arguments)]
fn solve_at(
    comps: &[HashMap<[u8; 4], u64>],
    d: u16,
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    row_selection: bool,
    debug: bool,
) -> Inner {
    assert!(p < 1 << 15, "rows are u16 and products must fit u32");
    let t0 = Instant::now();
    st.equation_degree = EQUATION_WDEG as u8;
    st.macaulay_degree = d as u8;
    let cols = weighted_monomials_desc(d);
    let col_index: HashMap<[u8; 4], usize> =
        cols.iter().enumerate().map(|(i, m)| (*m, i)).collect();
    let mut tri = 0u64;
    let (comps, leads) = triangularise(comps, p, &mut tri);
    st.triangularise_muls += tri;
    let divides = |lm: &[u8; 4], m: &[u8; 4]| (0..4).all(|k| lm[k] <= m[k]);
    let mut rows: Vec<Vec<u16>> = Vec::new();
    for (i, comp) in comps.iter().enumerate() {
        let Some(dg) = comp
            .iter()
            .filter(|(_, &v)| v != 0)
            .map(|(m, _)| wdeg(m))
            .max()
        else {
            continue;
        };
        if dg > d {
            return Inner::Unsolved("an equation above the Macaulay degree".into());
        }
        for sh in weighted_monomials_desc(d - dg) {
            if row_selection
                && leads
                    .iter()
                    .take(i)
                    .any(|lm| lm.is_some_and(|lm| divides(&lm, &sh)))
            {
                continue;
            }
            let mut row = vec![0u16; cols.len()];
            for (&e, &c) in comp {
                let m: [u8; 4] = core::array::from_fn(|k| e[k] + sh[k]);
                let ci = col_index[&m];
                row[ci] = am(row[ci] as u64, c, p) as u16;
            }
            rows.push(row);
        }
    }
    st.rows = rows.len();
    st.cols = cols.len();
    st.build_ms += t0.elapsed().as_secs_f64() * 1e3;
    if debug {
        eprintln!(
            "  weighted degree {d}: {} rows × {} columns ({:.3} GB)",
            st.rows,
            st.cols,
            (st.rows * st.cols * 2) as f64 / 1e9
        );
    }
    let t1 = Instant::now();
    let mut ech_muls = 0u64;
    let ech = echelon_by_lead(rows, p as u32, &mut ech_muls, false);
    st.echelon_muls += ech_muls;
    st.echelon_ms += t1.elapsed().as_secs_f64() * 1e3;
    st.zero_rows = ech.zero_rows;
    st.pivots = ech.pivot_row.iter().filter(|&&r| r != u32::MAX).count();
    let is_std: Vec<bool> = (0..cols.len())
        .map(|c| ech.pivot_row[c] == u32::MAX && wdeg(&cols[c]) < d)
        .collect();
    let standard: Vec<usize> = (0..cols.len()).filter(|&c| is_std[c]).collect();
    let dim = standard.len();
    st.dim = dim;
    let mut std_index = vec![u32::MAX; cols.len()];
    for (i, &c) in standard.iter().enumerate() {
        std_index[c] = i as u32;
    }
    let mut low: Vec<usize> = Vec::with_capacity(5);
    for v in 0..=4usize {
        let m: [u8; 4] = core::array::from_fn(|k| u8::from(v > 0 && k == v - 1));
        low.push(col_index[&m]);
    }
    if low.iter().any(|&c| !is_std[c]) {
        return Inner::Unsolved("1 or an unknown is not standard".into());
    }
    let idx: Vec<usize> = low.iter().map(|&c| std_index[c] as usize).collect();
    let t2 = Instant::now();
    let mut nfs = NormalForms4 {
        ech: &ech,
        is_std,
        std_index,
        dim,
        p: p as u32,
        memo: vec![None; cols.len()],
        reached: 0,
    };
    let mult_matrix =
        |nfs: &mut NormalForms4, var: usize, muls: &mut u64| -> Option<Vec<Vec<u64>>> {
            let mut targets = Vec::with_capacity(dim);
            for &bc in &standard {
                let mut m = cols[bc];
                m[var] += 1;
                targets.push(*col_index.get(&m)?);
            }
            nfs.ensure(&targets, muls).ok()?;
            let mut out = vec![vec![0u64; dim]; dim];
            for (bi, &tc) in targets.iter().enumerate() {
                for (ri, v) in nfs.vector(tc).into_iter().enumerate() {
                    out[ri][bi] = v;
                }
            }
            Some(out)
        };
    let mut nf_muls = 0u64;
    let m_pi = mult_matrix(&mut nfs, EIG_VAR, &mut nf_muls);
    st.nf_muls += nf_muls;
    st.nf_ms += t2.elapsed().as_secs_f64() * 1e3;
    st.reached_pivots = nfs.reached;
    let Some(m_pi) = m_pi else {
        return Inner::Unsolved("border unreachable at this degree".into());
    };
    let t3 = Instant::now();
    let mut cp_muls = 0u64;
    let cp = charpoly_u16(&m_pi, p as u32, &mut cp_muls);
    st.charpoly_muls += cp_muls;
    st.charpoly_ms += t3.elapsed().as_secs_f64() * 1e3;
    let t4 = Instant::now();
    let ring = PolyRing::new(p);
    let lambdas = ring.roots(&cp, rng);
    st.roots_muls += ring.muls.get();
    st.rational_eigenvalues = lambdas.len();
    let mut mt = vec![vec![0u16; dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            mt[j][i] = m_pi[i][j] as u16;
        }
    }
    drop(m_pi);
    let mut others: Vec<Option<Vec<Vec<u64>>>> = vec![None; 3];
    let mut out = Vec::new();
    for lam in lambdas {
        let mut ev = 0u64;
        let mut ker = eigenvectors_u16(&mt, lam, p as u32, &mut ev);
        if ker.len() > 1 {
            st.degenerate_eigenspaces += 1;
            for (k, slot) in others.iter_mut().enumerate() {
                if ker.len() <= 1 {
                    break;
                }
                if slot.is_none() {
                    let mut nm = 0u64;
                    *slot = mult_matrix(&mut nfs, k, &mut nm);
                    st.nf_muls += nm;
                }
                let Some(mv) = slot.as_ref() else {
                    st.eigenvector_muls += ev;
                    return Inner::NeedHigher;
                };
                ker = split_eigenspace(&ker, mv, p, rng, &mut ev);
            }
        }
        st.eigenvector_muls += ev;
        for w in &ker {
            if w[idx[0]] == 0 {
                continue;
            }
            let inv = inv_mod(w[idx[0]], p);
            let e: [u64; 4] = core::array::from_fn(|k| mm(w[idx[k + 1]], inv, p));
            if e[EIG_VAR] != lam {
                continue;
            }
            let mut vm = 0u64;
            let ok = comps.iter().all(|c| eval(c, &e, p, &mut vm) == 0);
            st.verify_muls += vm;
            if ok {
                out.push(e);
            }
        }
    }
    st.eigen_ms += t4.elapsed().as_secs_f64() * 1e3;
    out.sort_unstable();
    out.dedup();
    st.e_solutions = out.len();
    st.outcome = "solved".into();
    Inner::Solved(out)
}

/// The weighted solve: at `d` (normally [`MACAULAY_BOUND`]), and at `d + 1`
/// if a degenerate eigenspace of `M_π` needs a weight-2 unknown to split.
/// Both runs are charged in `st`.  Returns the `F_p`-rational solutions
/// `(σ₁, σ₂, σ₃, π)` and the degree that produced them.
pub fn solve_weighted(
    comps: &[HashMap<[u8; 4], u64>],
    d: u16,
    p: u64,
    rng: &mut StdRng,
    st: &mut QuarticSolve,
    row_selection: bool,
    debug: bool,
) -> (Option<Vec<[u64; 4]>>, u16) {
    let t0 = Instant::now();
    let mut deg = d;
    let res = loop {
        match solve_at(comps, deg, p, rng, st, row_selection, debug) {
            Inner::Solved(s) => break Some(s),
            Inner::Unsolved(why) => {
                st.outcome = why;
                break None;
            }
            Inner::NeedHigher if deg < d + 2 => {
                if debug {
                    eprintln!(
                        "  degenerate eigenspace: rerun at weighted degree {}",
                        deg + 1
                    );
                }
                deg += 1;
            }
            Inner::NeedHigher => {
                st.outcome = "degenerate eigenspace unresolved".into();
                break None;
            }
        }
    };
    st.total_ms = t0.elapsed().as_secs_f64() * 1e3;
    st.close_public();
    (res, deg)
}

/// What recovering the decompositions from the invariant solutions cost.
#[derive(Clone, Debug, Default, Serialize)]
pub struct LiftBack {
    /// `F_p` multiplications: split of `T⁴ − σ₁T³ + σ₂T² − σ₃T + π²`,
    /// square roots in `F_p`, `x(Y)` and the lift to `E(F_{p⁴})`, and the
    /// signed-sum checks (group additions at their counted cost).
    pub muls: u64,
    /// Rational invariant solutions whose `Y_i²` split in `F_p`.
    pub split: usize,
    /// ... and whose `Y_i` are in `F_p` (all four squares are squares).
    pub rooted: usize,
    /// Even-sign representatives tried (each an unordered `Y`-quadruple).
    pub candidates: usize,
    /// Group additions spent on the signed-sum checks.
    pub group_additions: u64,
    /// The verified decompositions, as sorted `Y`-quadruples.
    pub decompositions: Vec<[u64; 4]>,
}

/// Tonelli–Shanks in `F_p` with every multiplication counted.
fn sqrt_counted(n: u64, p: u64, muls: &mut u64) -> Option<u64> {
    let n = n % p;
    if n == 0 {
        return Some(0);
    }
    let pw = |b: u64, e: u64, muls: &mut u64| -> u64 {
        *muls += 2 * (64 - e.leading_zeros() as u64);
        pow_mod(b, e, p)
    };
    if pw(n, (p - 1) / 2, muls) != 1 {
        return None;
    }
    let mut q = p - 1;
    let mut s = 0u32;
    while q.is_multiple_of(2) {
        q /= 2;
        s += 1;
    }
    let mut z = 2u64;
    while pw(z, (p - 1) / 2, muls) != p - 1 {
        z += 1;
    }
    let mut m = s;
    let mut c = pw(z, q, muls);
    let mut t = pw(n, q, muls);
    let mut r = pw(n, q.div_ceil(2), muls);
    while t != 1 {
        let mut i = 0u32;
        let mut tt = t;
        while tt != 1 {
            tt = mm(tt, tt, p);
            *muls += 1;
            i += 1;
            if i == m {
                return None;
            }
        }
        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = mm(b, b, p);
            *muls += 1;
        }
        m = i;
        c = mm(b, b, p);
        t = mm(t, c, p);
        r = mm(r, b, p);
        *muls += 3;
    }
    Some(r)
}

/// `±P₁ ± P₂ ± P₃ ± P₄ = R` for some signs (the first fixed to `+`, the
/// match taken on `x`, so `−R` is the all-flipped pattern).
fn signed_sum_hits(curve: &Curve4, pts: &[Pt4; 4], r: &Pt4, adds: &mut u64) -> bool {
    for s in 0..8u32 {
        let mut acc = pts[0];
        for (i, q) in pts.iter().enumerate().skip(1) {
            let q = if s >> (i - 1) & 1 == 1 {
                curve.neg(q)
            } else {
                *q
            };
            acc = curve.add(&acc, &q);
            *adds += 1;
        }
        if !acc.inf && acc.x == r.x {
            return true;
        }
    }
    false
}

/// Recover the decompositions of `R` from the rational invariant solutions:
/// split for `Y_i²`, square roots, the eight sign vectors with product `π`
/// (the even `T`-translates of one decomposition), the points, and the
/// group check of every translate.  Every multiplication is charged.
pub fn lift_back(inst: &D4Instance, sols: &[[u64; 4]], r: &Pt4, rng: &mut StdRng) -> LiftBack {
    let f = &inst.curve.f;
    let p = f.p;
    let mut lb = LiftBack::default();
    let before = f.muls();
    let mut found = BTreeSet::new();
    let mut fp_muls = 0u64;
    for s in sols {
        let pi2 = mm(s[3], s[3], p);
        fp_muls += 1;
        let Some(sq) = split_quartic(&[s[0], s[1], s[2], pi2], p, rng, &mut fp_muls) else {
            continue;
        };
        lb.split += 1;
        let roots: Option<Vec<u64>> = sq
            .iter()
            .map(|&v| sqrt_counted(v, p, &mut fp_muls))
            .collect();
        let Some(roots) = roots else {
            continue;
        };
        lb.rooted += 1;
        for signs in 0..16u32 {
            let ys: [u64; 4] = core::array::from_fn(|i| {
                if signs >> i & 1 == 1 {
                    (p - roots[i]) % p
                } else {
                    roots[i]
                }
            });
            let prod = ys.iter().fold(1u64, |a, &v| a * v % p);
            fp_muls += 3;
            if prod != s[3] {
                continue;
            }
            if ys.iter().any(|&v| v == 1 || v == p - 1) {
                continue;
            }
            let mut sorted = ys;
            sorted.sort_unstable();
            if sorted.windows(2).any(|w| w[0] == w[1]) {
                continue;
            }
            lb.candidates += 1;
            let pts: Option<Vec<Pt4>> = ys
                .iter()
                .map(|&v| inst.curve.lift_x(&inst.x_of_y(&f.from_fp(v)), rng))
                .collect();
            let Some(pts) = pts else {
                continue;
            };
            let pts: [Pt4; 4] = [pts[0], pts[1], pts[2], pts[3]];
            if signed_sum_hits(&inst.curve, &pts, r, &mut lb.group_additions) {
                found.insert(sorted);
            }
        }
    }
    lb.muls = fp_muls + (f.muls() - before);
    lb.decompositions = found.into_iter().collect();
    lb
}

/// Every four-subset of distinct base points with a signed sum equal to
/// `R`, as sorted `Y`-quadruples: `C(|F|, 4) · 8` sums, the ground truth
/// on tiny fields.
pub fn brute_force_decompositions(
    curve: &Curve4,
    base: &[(u64, Pt4)],
    r: &Pt4,
) -> BTreeSet<[u64; 4]> {
    let n = base.len();
    let mut out = BTreeSet::new();
    let pm = |q: &Pt4, neg: bool| if neg { curve.neg(q) } else { *q };
    for a in 0..n {
        for b in (a + 1)..n {
            for sb in [false, true] {
                let ab = curve.add(&base[a].1, &pm(&base[b].1, sb));
                for c in (b + 1)..n {
                    for sc in [false, true] {
                        let abc = curve.add(&ab, &pm(&base[c].1, sc));
                        for (dd, item) in base.iter().enumerate().skip(c + 1) {
                            for sd in [false, true] {
                                let s = curve.add(&abc, &pm(&item.1, sd));
                                if !s.inf && s.x == r.x {
                                    let mut ys = [base[a].0, base[b].0, base[c].0, base[dd].0];
                                    ys.sort_unstable();
                                    out.insert(ys);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

/// A planted generic system for the solver's positive control: four
/// equations with random coefficients on every monomial of weighted degree
/// `≤ 8` and a common root at a random point of `F_p⁴`.
pub fn planted_generic(p: u64, seed: u64) -> (Vec<HashMap<[u8; 4], u64>>, [u64; 4]) {
    let mut rng = StdRng::seed_from_u64(seed ^ 0x6e6e);
    let root: [u64; 4] = core::array::from_fn(|_| rng.gen_range(1..p));
    let monos = weighted_monomials_desc(EQUATION_WDEG);
    let comps = (0..4)
        .map(|_| {
            let mut c: HashMap<[u8; 4], u64> = HashMap::new();
            for m in &monos {
                if *m != [0, 0, 0, 0] {
                    c.insert(*m, rng.gen_range(0..p));
                }
            }
            let mut sc = 0u64;
            let v = eval(&c, &root, p, &mut sc);
            c.insert([0, 0, 0, 0], (p - v) % p);
            c.retain(|_, v| *v != 0);
            c
        })
        .collect();
    (comps, root)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::gaudry_quartic_d4::{generate_d4_instance, y_base};

    #[test]
    fn weighted_bound_and_monomial_count() {
        assert_eq!(MACAULAY_BOUND, 26);
        // The 105 invariant monomials of H are exactly weighted degree ≤ 8.
        assert_eq!(weighted_monomials_desc(8).len(), 105);
        assert_eq!(columns_at(26), 4200);
    }

    #[test]
    fn planted_generic_system_has_quotient_512_and_its_root() {
        let p = 269;
        let (comps, root) = planted_generic(p, 3);
        let mut rng = StdRng::seed_from_u64(1);
        let mut st = QuarticSolve::default();
        let (sols, _) = solve_weighted(&comps, MACAULAY_BOUND, p, &mut rng, &mut st, true, false);
        assert_eq!(st.dim, 512);
        assert!(sols.expect("solved").contains(&root));
    }

    #[test]
    fn weighted_solve_agrees_with_brute_force_on_a_tiny_field() {
        let p = 13;
        let inst = generate_d4_instance(p, 2);
        let mut rng = StdRng::seed_from_u64(5);
        let pre = SymmetrisedD4::precompute(&inst, &mut rng, 4).expect("invariant");
        let base = y_base(&inst, &mut rng);
        // A planted target: a signed sum of four distinct base points.
        let r = {
            let c = &inst.curve;
            let s = c.add(
                &c.add(&base[0].1, &base[1].1),
                &c.sub(&base[2].1, &base[3].1),
            );
            assert!(!s.inf);
            s
        };
        let y_r = inst.y_of_x(&r.x);
        let comps = weighted_system(&pre, &inst.curve.f, &y_r);
        let mut st = QuarticSolve::default();
        let (sols, _) = solve_weighted(&comps, MACAULAY_BOUND, p, &mut rng, &mut st, true, false);
        let sols = sols.expect("solved");
        let lb = lift_back(&inst, &sols, &r, &mut rng);
        let truth = brute_force_decompositions(&inst.curve, &base, &r);
        let mine: BTreeSet<[u64; 4]> = lb.decompositions.iter().copied().collect();
        assert!(!truth.is_empty());
        assert_eq!(mine, truth);
    }
}
