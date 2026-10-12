//! The arms of the PDP-oracle speedup comparison: every arm decides one
//! target (a rational point `R`) and reports whether it decomposes and how
//! many pairs it visited.

use crate::common::*;
use crate::group::*;
use crate::kernel::ScalarKernel;
use crypto_lib::cryptanalysis::semaev_decomp::{decompose, subspace_poly, SubspaceOracle};

pub enum Arm<'a> {
    /// `semaev_decomp::decompose`: every abscissa of `V`.
    Baseline(&'a Rung),
    /// `SubspaceOracle` over all of `V` (the baseline algorithm, general-`b` code path).
    SubspaceFull(&'a Rung, SubspaceOracle),
    /// `SubspaceOracle` over the liftable abscissae of `V` only: the baseline's
    /// root test on the candidate set the group-law kernel uses (matched control).
    SubspaceLiftable(&'a Rung, SubspaceOracle),
    /// The group-law row kernel, scalar.
    GroupScalar(&'a Rung, ScalarKernel<'a>),
}

impl<'a> Arm<'a> {
    pub fn build(name: &str, rung: &'a Rung) -> Self {
        match name {
            "baseline" => Arm::Baseline(rung),
            "so-full" => Arm::SubspaceFull(rung, SubspaceOracle::low_order(rung.l, 1, &rung.gf)),
            "so-liftable" => {
                let solver = ArtinSchreier::new(&rung.gf);
                let span: Vec<u64> = base_points(&rung.gf, &solver, rung.l)
                    .into_iter()
                    .map(|p| p.0)
                    .collect();
                let lv = subspace_poly(rung.l, &rung.gf);
                Arm::SubspaceLiftable(
                    rung,
                    SubspaceOracle {
                        l: rung.l,
                        b: 1,
                        span,
                        lv,
                    },
                )
            }
            "group-scalar" => {
                let solver = ArtinSchreier::new(&rung.gf);
                Arm::GroupScalar(rung, ScalarKernel::new(&rung.gf, &solver, rung.l))
            }
            other => panic!("unknown arm {other}"),
        }
    }

    /// Pairs the arm visits when it refutes.
    pub fn pairs_total(&self) -> u64 {
        match self {
            Arm::Baseline(r) => r.pairs_total(),
            Arm::SubspaceFull(r, _) => r.pairs_total(),
            Arm::SubspaceLiftable(_, o) => {
                let m = o.span.len() as u64;
                m * (m + 1) / 2
            }
            Arm::GroupScalar(_, k) => k.pairs_total(),
        }
    }

    /// Decide `R = (xr, yr)`.  A target for which the group-law formulation
    /// may disagree with the baseline (`degenerate_target`) is decided by
    /// the baseline, inside this call.
    pub fn decide(&self, xr: u64, yr: u64) -> bool {
        match self {
            Arm::Baseline(r) => decompose(xr, r.l, &r.gf).is_some(),
            Arm::SubspaceFull(r, o) => o.decompose(xr, &r.gf).0.is_some(),
            Arm::SubspaceLiftable(r, o) => {
                if degenerate_target(&r.gf, r.l, xr) {
                    decompose(xr, r.l, &r.gf).is_some()
                } else {
                    o.decompose(xr, &r.gf).0.is_some()
                }
            }
            Arm::GroupScalar(r, k) => {
                if degenerate_target(&r.gf, r.l, xr) {
                    decompose(xr, r.l, &r.gf).is_some()
                } else {
                    let mut out = Vec::new();
                    k.run((xr, yr), false, &mut out).0
                }
            }
        }
    }

    /// Every sorted root, for the arms that can list them (group kernel).
    pub fn all_roots(&self, xr: u64, yr: u64) -> Option<Vec<[u64; 3]>> {
        match self {
            Arm::GroupScalar(r, k) if !degenerate_target(&r.gf, r.l, xr) => {
                let mut out = Vec::new();
                k.run((xr, yr), true, &mut out);
                out.sort_unstable();
                out.dedup();
                Some(out)
            }
            _ => None,
        }
    }
}
