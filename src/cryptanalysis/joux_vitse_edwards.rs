//! **The `k = 5` Joux–Vitse decomposition with the 2-torsion symmetry, end
//! to end: F4 trace replay on the Edwards four-point test, complete DLPs.**
//!
//! Companion to `research/joux_vitse_n1_20261015/` (the plain Weierstrass
//! route, measured end to end) and to [`super::jv_quintic_edwards`], whose
//! construction this module uses unchanged: an Edwards curve
//! `x² + y² = 1 + d x² y²` over `F_{p⁵}` with `#E = 4n`, `n` prime; the
//! factor base `{P : y(P) ∈ F_p}`, one point per `y²` class; the symmetrised
//! `S₅` in the elementary symmetric functions of the `y_i²` and the product
//! variable `π` (`π² = e₄`); the meet-in-the-middle oracle
//! ([`mitm4_signed`]).  `RESEARCH_K5_TORSION_JOUX_VITSE.md` measured one
//! test at `C″ = 7.3·10⁷` `F_p` multiplications and named two things not
//! done (§5–§6): the trace-driven elimination (Joux–Vitse's "F4 remake")
//! and an end-to-end run.  This module adds both:
//!
//! * [`trial`]: the descended system (five Weil components and
//!   `π² − e₄`), consistency by replaying an F4 trace recorded on the
//!   curve's first inconsistent system ([`f4_fp::f4_replay`]); a replay
//!   that diverges, or a consistent system, goes to the committed test's
//!   solve (`f4_fp::solve`, both charged), and the decompositions are
//!   split, resolved and re-added in the group as
//!   [`super::jv_quintic_edwards::jv4_decompose`] does;
//! * [`dlp_rows`] and [`super::joux_vitse::solve_logs`]: a relation
//!   `R = Σ s_i P_i + εT` with `R = aG + bQ` in the subgroup of order `n`
//!   becomes, through `[4]` (which kills `T`), `Σ s_i L_i − 4 b k ≡ 4 a
//!   (mod n)` with `L_i = log_G(4 P_i)`;
//! * [`rho`]: Pollard rho on the same subgroup, the paired reference.

use std::collections::{BTreeSet, HashMap};

use rand::rngs::StdRng;
use rand::Rng;
use serde::Serialize;

use super::f4_fp::{self, F4Options, F4Trace, Ordering as F4Ordering, Verdict};
use super::gaudry_cubic::{PolyRing, UPoly};
use super::jv_quintic_edwards::{
    mitm4_signed, polys_of, resolve, verify_quad, Ed5, PairTableE, PtE, QuadE, SymmetrisedS5Ed,
};

/// The committed test's F4 options without the wall-clock budget (a budget
/// would make a verdict depend on the host): grevlex, pairs bounded at
/// degree 24, as `experiments/27_jv_quintic_edwards_csecond.json` ran.
pub fn trial_options() -> F4Options {
    F4Options::new(F4Ordering::Grevlex, 24)
}

/// The descended system of a residual: five `F_p`-components of
/// `A(e; y_R) + π B(e; y_R)` and `π² − e₄`, in `(e₁, e₂, e₃, e₄, π)`.
pub fn system(curve: &Ed5, pre: &SymmetrisedS5Ed, r: &PtE) -> Vec<f4_fp::Poly> {
    polys_of(&pre.weil_restrict(&curve.f, &r.y), curve.f.p)
}

/// Record a trace on the first inconsistent system among residuals drawn
/// by `draw`; returns the trace, the row-reduction multiplications spent and
/// the number of systems tried.
pub fn record_trace(
    curve: &Ed5,
    pre: &SymmetrisedS5Ed,
    opts: &F4Options,
    mut draw: impl FnMut() -> PtE,
) -> (F4Trace, u64, u64) {
    let (mut muls, mut tries) = (0u64, 0u64);
    loop {
        tries += 1;
        let r = draw();
        let polys = system(curve, pre, &r);
        let (rep, tr) = f4_fp::f4_record(&polys, 5, curve.f.p, opts);
        muls += rep.field_ops;
        if rep.inconsistent {
            return (tr, muls, tries);
        }
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct TrialE {
    pub consistent: bool,
    pub replayed: bool,
    pub replay_mismatch: bool,
    pub undetermined: bool,
    /// `F_p` multiplications building the system (Weil restriction).
    pub weil_muls: u64,
    /// Row reductions of the replay (zero without a trace).
    pub replay_muls: u64,
    /// Row reductions of the committed solve (after a mismatch, on a
    /// consistent system, or every time without a trace).
    pub solve_muls: u64,
    /// Root finding of the quartics in `y²`.
    pub roots_muls: u64,
    /// Group operations resolving signs and `T`.
    pub group_ops: u64,
    pub e_solutions: usize,
    pub not_split: usize,
    pub decompositions: usize,
}

/// One trial on `r`: replay `trace` when given; otherwise, or on a
/// mismatch or a consistent system, the committed test's solve.  The
/// per-report row-reduction counts are exact when trials run one at a time;
/// under concurrency read [`f4_fp::field_ops_total`] around the batch.
pub fn trial(
    curve: &Ed5,
    base: &[PtE],
    by_y2: &HashMap<u64, usize>,
    pre: &SymmetrisedS5Ed,
    r: &PtE,
    opts: &F4Options,
    trace: Option<&F4Trace>,
    rng: &mut StdRng,
) -> (Vec<QuadE>, TrialE) {
    let f = &curve.f;
    let p = f.p;
    let mut out = TrialE::default();
    let m0 = f.muls();
    let polys = system(curve, pre, r);
    out.weil_muls = f.muls() - m0;
    if let Some(tr) = trace {
        let rr = f4_fp::f4_replay(&polys, 5, p, opts, tr);
        out.replay_muls = rr.field_ops;
        if rr.trace_mismatch {
            out.replay_mismatch = true;
        } else if rr.inconsistent {
            out.replayed = true;
            return (Vec::new(), out);
        } else {
            out.replayed = true;
        }
    }
    let rep = f4_fp::solve(&polys, 5, p, opts);
    out.solve_muls = rep.field_ops;
    let es = match rep.verdict {
        Verdict::Inconsistent => Vec::new(),
        Verdict::Undetermined => {
            out.undetermined = true;
            Vec::new()
        }
        Verdict::Solutions(s) => s,
    };
    out.consistent = !es.is_empty();
    out.e_solutions = es.len();
    let ring = PolyRing::new(p);
    let g0 = curve.ops();
    let mut found: BTreeSet<QuadE> = BTreeSet::new();
    for e in &es {
        let quartic = UPoly(vec![e[3], (p - e[2]) % p, e[1], (p - e[0]) % p, 1]);
        let roots = ring.roots(&quartic, rng);
        if roots.len() != 4 {
            out.not_split += 1;
            continue;
        }
        let Some(idx) = roots
            .iter()
            .map(|x| by_y2.get(x).copied())
            .collect::<Option<Vec<usize>>>()
        else {
            continue;
        };
        let mut sorted = idx.clone();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() != 4 {
            continue;
        }
        if let Some(q) = resolve(curve, base, r, [idx[0], idx[1], idx[2], idx[3]]) {
            found.insert(q);
        }
    }
    out.roots_muls = ring.muls.get();
    out.group_ops = curve.ops() - g0;
    out.decompositions = found.len();
    (found.into_iter().collect(), out)
}

/// The trial and the oracle agree, and every decomposition re-adds to `r`.
pub fn check(
    curve: &Ed5,
    base: &[PtE],
    table: &PairTableE,
    r: &PtE,
    quads: &[QuadE],
) -> (bool, bool) {
    let oracle = mitm4_signed(curve, base, table, r);
    let verified = quads.iter().all(|q| verify_quad(curve, base, r, q));
    (oracle == quads, verified)
}

/// The logarithm row of a relation `aG + bQ = Σ s_i P_i + εT`, for
/// [`super::joux_vitse::solve_logs`] with cofactor `h = 4`.
pub fn dlp_row(q: &QuadE, a: u64, b: u64) -> (Vec<(usize, i8)>, u128, u128) {
    (
        q.terms.iter().map(|&(i, s)| (i, s as i8)).collect(),
        a as u128,
        b as u128,
    )
}

/// Pollard rho on `⟨G⟩` of prime order `n`: an `r`-adding walk (`r = 20`)
/// with Brent's cycle detection and no negation map, as
/// [`super::joux_vitse::rho`] runs on the Weierstrass curves.  Every group
/// operation is counted, the multipliers' set-up included.
pub fn rho(curve: &Ed5, g: &PtE, q: &PtE, n: u64, rng: &mut StdRng) -> (Option<u64>, u64) {
    const R: usize = 20;
    let ops0 = curve.ops();
    let mulm = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    let inv = |a: u64| -> u64 {
        let (mut t, mut nt) = (0i128, 1i128);
        let (mut r, mut nr) = (n as i128, a as i128);
        while nr != 0 {
            let k = r / nr;
            (t, nt) = (nt, t - k * nt);
            (r, nr) = (nr, r - k * nr);
        }
        t.rem_euclid(n as i128) as u64
    };
    loop {
        let mults: Vec<(u64, u64, PtE)> = (0..R)
            .map(|_| {
                let a = rng.gen_range(0..n);
                let b = rng.gen_range(0..n);
                (a, b, curve.add(&curve.mul(g, a), &curve.mul(q, b)))
            })
            .collect();
        let idx = |pt: &PtE| -> usize {
            ((pt.y.0[0] * 7 + pt.y.0[1] * 3 + pt.x.0[0]) % R as u64) as usize
        };
        let (mut a, mut b) = (rng.gen_range(0..n), rng.gen_range(0..n));
        let mut x = curve.add(&curve.mul(g, a), &curve.mul(q, b));
        let (mut sa, mut sb, mut sx) = (a, b, x);
        let (mut power, mut lam) = (1u64, 0u64);
        loop {
            let j = idx(&x);
            x = curve.add(&x, &mults[j].2);
            a = (a + mults[j].0) % n;
            b = (b + mults[j].1) % n;
            lam += 1;
            if x == sx {
                let db = (sb + n - b) % n;
                if db == 0 {
                    break;
                }
                let k = mulm((a + n - sa) % n, inv(db));
                return (Some(k), curve.ops() - ops0);
            }
            if lam == power {
                (sa, sb, sx) = (a, b, x);
                power *= 2;
                lam = 0;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::jv_quintic_edwards::{
        factor_base, generate_instance_edwards, pair_table, y2_key,
    };
    use rand::SeedableRng;

    /// Replay agrees with the committed test and the oracle on random and
    /// planted residuals, and a planted relation re-adds and enters a row.
    #[test]
    fn replayed_trial_matches_oracle_and_committed_test() {
        let inst = generate_instance_edwards(31, 1);
        let curve = &inst.curve;
        let mut rng = StdRng::seed_from_u64(7);
        let base = factor_base(curve, &mut rng);
        assert!(base.len() >= 4);
        let by_y2: HashMap<u64, usize> = base
            .iter()
            .enumerate()
            .map(|(i, q)| (y2_key(curve, q), i))
            .collect();
        let pre = SymmetrisedS5Ed::precompute(curve, &mut rng);
        let table = pair_table(curve, &base);
        let opts = trial_options();
        let n = inst.n;
        let (tr, _, _) = record_trace(curve, &pre, &opts, || {
            let a = rng.gen_range(0..n);
            curve.mul(&inst.g, a)
        });
        let planted = curve.sub(
            &curve.add(&base[0], &base[1]),
            &curve.add(&base[2], &base[3]),
        );
        let mut rng2 = StdRng::seed_from_u64(9);
        for r in [planted, curve.mul(&inst.g, 12345), curve.mul(&inst.q, 777)] {
            let (q1, t1) = trial(curve, &base, &by_y2, &pre, &r, &opts, Some(&tr), &mut rng2);
            let (q0, _) = trial(curve, &base, &by_y2, &pre, &r, &opts, None, &mut rng2);
            assert_eq!(q1, q0);
            let (agree, verified) = check(curve, &base, &table, &r, &q1);
            assert!(agree && verified, "{t1:?}");
        }
        let (q, _) = trial(
            curve,
            &base,
            &by_y2,
            &pre,
            &planted,
            &opts,
            Some(&tr),
            &mut rng2,
        );
        assert!(!q.is_empty());
    }
}
