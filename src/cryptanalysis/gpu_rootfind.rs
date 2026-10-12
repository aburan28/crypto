//! CPU half of the GPU root-finding 2-decomposition oracle
//! (`gpu/rootfind/`): the [`Relation`] type the oracle reports, the CPU
//! reference it is validated against, and the fail-closed group check.
//!
//! The oracle answers, for a batch of targets `R_k` and a factor base
//! `F = {F_0 … F_{N−1}}` on `y² = x³ + ax + b`, every
//! `R_k = ε·F_i + τ·F_j` (`j ≥ i`) by the root-finding identity of
//! [`crate::cryptanalysis::prime_fast::collect_relations`]: the roots of
//! Semaev's `S₃(x_R, x_i, X)` are `x(R ± F_i)`, so each `(k, i)` pair is
//! two slope computations and two table probes.  The C/CUDA side
//! implements that per `(target, slice)` thread.
//!
//! The FFI bindings (`decompose64`, `decompose_p256`) and the validation
//! harness live in the separate crate `gpu/rootfind/rs/` (`gpu-rootfind`),
//! which depends on this crate and compiles the native library in its own
//! `build.rs`:
//!
//! ```bash
//! cargo test --release --manifest-path gpu/rootfind/rs/Cargo.toml
//! cargo run --release --manifest-path gpu/rootfind/rs/Cargo.toml --bin gpu_rootfind_validate
//! ```
//!
//! That split is deliberate: a build script in the root package makes
//! cargo fingerprint the whole repository (which fails on the CI runners),
//! and an optional path dependency has to exist even when its feature is
//! off (which breaks the sparse-checkout workflows).  Here the root crate
//! carries no feature, dependency or build script for the oracle.
//!
//! This module provides [`cpu_reference64`], the identical sweep written
//! directly on [`prime_fast`] primitives (`Fp64::batch_inv`,
//! `FastFactorBase::lookup`), so a harness can demand set equality
//! between the two on every target, and [`verify_relation64`], the
//! re-add-in-the-group check every reported relation must pass.
//!
//! Claim hygiene: this is the decomposition *subroutine* of index
//! calculus on public synthetic factor bases.  Nothing here changes the
//! asymptotics of the ECDLP or says anything about deployed curves.
//!
//! [`prime_fast`]: crate::cryptanalysis::prime_fast

use std::collections::BTreeSet;

use crate::cryptanalysis::prime_fast::{FastCurve, FastFactorBase, Pt};

/// One decomposition `R_target = eps·F_i + tau·F_j`.  `tau == 0` marks
/// the degenerate `R = eps·F_i` (then `j == i`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Relation {
    pub target: u32,
    pub i: u32,
    pub j: u32,
    pub eps: i8,
    pub tau: i8,
}

// ── CPU reference (prime_fast primitives) ─────────────────────────────

/// The same sweep as the kernel, written on [`prime_fast`] primitives:
/// one `Fp64::batch_inv` per target over `x_R − x_i`, both slopes, and a
/// `FastFactorBase::lookup` per `(i, sign)`.  This is the inner loop of
/// `collect_relations` with the output turned into index-level
/// relations, so agreement with the GPU is agreement with the production
/// CPU oracle.
///
/// [`prime_fast`]: crate::cryptanalysis::prime_fast
pub fn cpu_reference64(
    curve: &FastCurve,
    fb: &FastFactorBase,
    targets: &[Pt],
) -> BTreeSet<Relation> {
    let f = curve.f;
    let n = fb.len();
    let mut d = vec![0u64; n];
    let mut degenerate = vec![false; n];
    let mut scratch = Vec::with_capacity(n);
    let mut out = BTreeSet::new();
    for (k, &r) in targets.iter().enumerate() {
        for i in 0..n {
            let dd = f.sub(r.x, fb.pts[i].x);
            degenerate[i] = dd == 0;
            d[i] = if dd == 0 { f.one } else { dd };
        }
        f.batch_inv(&mut d, &mut scratch);
        for i in 0..n {
            let fi = fb.pts[i];
            if degenerate[i] {
                let eps = if r.y == fi.y { 1i8 } else { -1 };
                out.insert(Relation {
                    target: k as u32,
                    i: i as u32,
                    j: i as u32,
                    eps,
                    tau: 0,
                });
                continue;
            }
            let inv = d[i];
            let lam_minus = f.mul(f.add(r.y, fi.y), inv);
            let lam_plus = f.mul(f.sub(r.y, fi.y), inv);
            for (eps, lam) in [(1i8, lam_minus), (-1i8, lam_plus)] {
                let x3 = f.sub(f.sub(f.sqr(lam), r.x), fi.x);
                if let Some(j) = fb.lookup(f.from_mont(x3)) {
                    if j >= i {
                        let y3 = f.sub(f.mul(lam, f.sub(r.x, x3)), r.y);
                        let tau = if y3 == fb.pts[j].y { 1i8 } else { -1 };
                        out.insert(Relation {
                            target: k as u32,
                            i: i as u32,
                            j: j as u32,
                            eps,
                            tau,
                        });
                    }
                }
            }
        }
    }
    out
}

/// Fail-closed check: re-add `eps·F_i + tau·F_j` in the group and compare
/// with the target.
pub fn verify_relation64(
    curve: &FastCurve,
    fb: &FastFactorBase,
    target: Pt,
    rel: &Relation,
) -> bool {
    let (i, j) = (rel.i as usize, rel.j as usize);
    if i >= fb.len() || j >= fb.len() || j < i {
        return false;
    }
    let sign = |p: Pt, s: i8| if s < 0 { curve.neg(p) } else { p };
    let pi = sign(fb.pts[i], rel.eps);
    let sum = if rel.tau == 0 {
        if i != j {
            return false;
        }
        Some(pi)
    } else {
        curve.add(Some(pi), Some(sign(fb.pts[j], rel.tau)))
    };
    sum == Some(target)
}

/// Canonical coordinate arrays for a factor base, in the layout the FFI
/// takes.
pub fn factor_base_canonical(curve: &FastCurve, fb: &FastFactorBase) -> (Vec<u64>, Vec<u64>) {
    let mut xs = Vec::with_capacity(fb.len());
    let mut ys = Vec::with_capacity(fb.len());
    for &p in &fb.pts {
        let (x, y) = curve.canonical(p);
        xs.push(x);
        ys.push(y);
    }
    (xs, ys)
}

/// Canonical coordinate arrays for a target batch.
pub fn targets_canonical(curve: &FastCurve, targets: &[Pt]) -> (Vec<u64>, Vec<u64>) {
    let mut xs = Vec::with_capacity(targets.len());
    let mut ys = Vec::with_capacity(targets.len());
    for &p in targets {
        let (x, y) = curve.canonical(p);
        xs.push(x);
        ys.push(y);
    }
    (xs, ys)
}

/// Public synthetic known-answer instance for harnesses and tests: the
/// `a = −3` prime-order curve at `bits` bits, a small-x factor base of
/// `fb_size` points, and a target batch of up to 40 positive controls
/// `±F_i ± F_j`, `n_random` random on-curve points, and the degenerate
/// `R = F_d` last, with `d = min(3, fb_size − 1)` (so `R = F_3` for every
/// factor base of at least four points).  Returns
/// `(curve, fb, targets, n_controls)`.
///
/// # Panics
///
/// If `fb_size == 0` or the curve yields fewer points than requested.
pub fn synthetic_instance(
    bits: u32,
    fb_size: usize,
    n_random: usize,
    seed: u64,
) -> (FastCurve, FastFactorBase, Vec<Pt>, usize) {
    use crate::cryptanalysis::prime_fast::find_a3_curve;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    let (curve, _g) = find_a3_curve(bits, 3).expect("curve");
    let fb = FastFactorBase::build(&curve, fb_size);
    assert!(
        fb_size > 0 && fb.len() == fb_size,
        "synthetic_instance: factor base of {fb_size} points requested, {} built",
        fb.len()
    );
    let mut rng = StdRng::seed_from_u64(seed);
    let mut targets = Vec::new();
    for _ in 0..40 {
        let i = rng.gen_range(0..fb.len());
        let j = rng.gen_range(0..fb.len());
        let pi = if rng.gen_bool(0.5) {
            curve.neg(fb.pts[i])
        } else {
            fb.pts[i]
        };
        let pj = if rng.gen_bool(0.5) {
            curve.neg(fb.pts[j])
        } else {
            fb.pts[j]
        };
        if let Some(r) = curve.add(Some(pi), Some(pj)) {
            targets.push(r);
        }
    }
    let n_ctl = targets.len();
    while targets.len() < n_ctl + n_random {
        let x = rng.gen_range(1..curve.f.p);
        if let Some(p) = curve.lift_x(x) {
            targets.push(p);
        }
    }
    targets.push(fb.pts[3.min(fb.len() - 1)]);
    (curve, fb, targets, n_ctl)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpu_reference_relations_re_add_and_cover_controls() {
        let (curve, fb, targets, n_ctl) = synthetic_instance(40, 200, 100, 7);
        let rels = cpu_reference64(&curve, &fb, &targets);
        for r in &rels {
            assert!(
                verify_relation64(&curve, &fb, targets[r.target as usize], r),
                "{r:?}"
            );
        }
        for k in 0..n_ctl {
            assert!(
                rels.iter().any(|r| r.target as usize == k),
                "control {k} not decomposed"
            );
        }
        let last = (targets.len() - 1) as u32;
        assert!(rels.contains(&Relation {
            target: last,
            i: 3,
            j: 3,
            eps: 1,
            tau: 0
        }));
    }
}
