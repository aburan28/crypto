//! Bit-sliced pairs-and-solve by the group law (PS-2 gate count, PS-3
//! execution).  One *task* is one (row, block) cell: a row `i1` (two row
//! points `U± = P₁ ± R`, broadcast over lanes) against a block of base
//! points `P₂` laid out across the lanes.
//!
//! Per lane the stage computes, for each of `U₊`, `U₋`,
//!
//! ```text
//!   d = x_U + x₂,   j = 1/d,   λ = (y_U + y₂)·j,   λ̃ = λ + x₂·j
//!   t  = λ² + λ + d + 1,       t̃ = λ̃² + λ̃ + d + 1
//! ```
//!
//! and the lane hits when any of the four `t` lies in `V`.  The inverses
//! are taken by Montgomery's trick across the tasks of a batch (one
//! bit-sliced inversion per batch and chain).  A lane with `d = 0` is
//! *exceptional*: its `d` is replaced by 1 so the chain stays invertible,
//! and the caller recomputes that pair exactly.

use crate::bitslice::*;

/// A scalar field element broadcast to every lane.
pub fn broadcast<L: Lane>(v: u64, n: usize) -> E<L> {
    let mut o = zero_elt::<L>();
    for (i, plane) in o.iter_mut().enumerate().take(n) {
        *plane = if (v >> i) & 1 == 1 {
            L::ones()
        } else {
            L::zero()
        };
    }
    o
}

/// The row's two points, broadcast.
pub struct RowPlanes<L: Lane> {
    pub xu: E<L>,
    pub yu: E<L>,
    pub xv: E<L>,
    pub yv: E<L>,
}

impl<L: Lane> RowPlanes<L> {
    pub fn new(n: usize, xu: u64, yu: u64, xv: u64, yv: u64) -> Self {
        Self {
            xu: broadcast(xu, n),
            yu: broadcast(yu, n),
            xv: broadcast(xv, n),
            yv: broadcast(yv, n),
        }
    }
}

/// Stage A output for one task.
pub struct StageA<L: Lane> {
    pub du: E<L>,
    pub dv: E<L>,
    /// Lanes whose `d_U` or `d_V` is zero (before the fix-up).
    pub exceptional: L,
}

pub fn stage_a<L: Lane>(f: &BsField, row: &RowPlanes<L>, x2: &E<L>) -> StageA<L> {
    let du = f.add(&row.xu, x2);
    let dv = f.add(&row.xv, x2);
    let zu = f.zero_mask(&du);
    let zv = f.zero_mask(&dv);
    let mut du = du;
    let mut dv = dv;
    // d = 0 becomes d = 1 (plane 0 ^= mask), so the product chain survives.
    du[0] = du[0].xor(zu);
    dv[0] = dv[0].xor(zv);
    StageA {
        du,
        dv,
        exceptional: zu.or(zv),
    }
}

/// Montgomery's trick over `ds` (per-lane inverses of every element; the
/// elements must be non-zero in every lane).  One bit-sliced inversion.
pub fn montgomery<L: Lane>(f: &BsField, ds: &[&E<L>]) -> Vec<E<L>> {
    let t = ds.len();
    assert!(t >= 1);
    // pre[k] = d_0 · … · d_{k−1} (pre[0] = 1, represented implicitly).
    let mut pre: Vec<E<L>> = Vec::with_capacity(t);
    let mut acc = *ds[0];
    for k in 1..t {
        pre.push(acc); // pre[k] for k >= 1
        acc = f.mul(&acc, ds[k]);
    }
    let mut inv = f.inv(&acc);
    let mut out: Vec<E<L>> = vec![zero_elt::<L>(); t];
    for k in (1..t).rev() {
        out[k] = f.mul(&inv, &pre[k - 1]);
        inv = f.mul(&inv, ds[k]);
    }
    out[0] = inv;
    out
}

/// Stage C: the four candidates and the hit mask.
pub fn stage_c<L: Lane>(
    f: &BsField,
    row: &RowPlanes<L>,
    x2: &E<L>,
    y2: &E<L>,
    a: &StageA<L>,
    ju: &E<L>,
    jv: &E<L>,
    l: usize,
) -> L {
    let one_cand = |yu: &E<L>, d: &E<L>, j: &E<L>| -> (E<L>, E<L>) {
        let lam = f.mul(&f.add(yu, y2), j);
        let sq = f.sqr(&lam);
        let t = f.add_one(&f.add(&f.add(&sq, &lam), d));
        let lam2 = f.add(&lam, &f.mul(x2, j));
        let sq2 = f.sqr(&lam2);
        let t2 = f.add_one(&f.add(&f.add(&sq2, &lam2), d));
        (t, t2)
    };
    let (t1, t2) = one_cand(&row.yu, &a.du, ju);
    let (t3, t4) = one_cand(&row.yv, &a.dv, jv);
    let m1 = f.in_subspace_mask(&t1, l);
    let m2 = f.in_subspace_mask(&t2, l);
    let m3 = f.in_subspace_mask(&t3, l);
    let m4 = f.in_subspace_mask(&t4, l);
    m1.or(m2).or(m3.or(m4))
}

/// A batch of `T` tasks: stage A for each, one Montgomery inversion per chain
/// across the batch, stage C for each.  Returns per task the hit mask and the
/// mask of exceptional lanes (zero denominator, to be recomputed exactly).
pub fn batch<L: Lane>(
    f: &BsField,
    rows: &[RowPlanes<L>],
    x2s: &[E<L>],
    y2s: &[E<L>],
    l: usize,
) -> Vec<(L, L)> {
    let t = rows.len();
    let stages: Vec<StageA<L>> = (0..t).map(|k| stage_a(f, &rows[k], &x2s[k])).collect();
    let dus: Vec<&E<L>> = stages.iter().map(|s| &s.du).collect();
    let dvs: Vec<&E<L>> = stages.iter().map(|s| &s.dv).collect();
    let ju = montgomery(f, &dus);
    let jv = montgomery(f, &dvs);
    (0..t)
        .map(|k| {
            let hit = stage_c(f, &rows[k], &x2s[k], &y2s[k], &stages[k], &ju[k], &jv[k], l);
            (hit, stages[k].exceptional)
        })
        .collect()
}
