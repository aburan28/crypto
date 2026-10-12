//! The arms of the PDP-oracle speedup comparison: every arm decides one
//! target (a rational point `R`) and reports whether it decomposes and how
//! many pairs it visited.

use crate::bskernel::{BsKernel, Dispatch};
use crate::common::*;
use crate::group::*;
use crate::kernel::ScalarKernel;
use crate::mitm::MitmTable;
use crypto_lib::cryptanalysis::semaev_decomp::{decompose, subspace_poly, SubspaceOracle};

/// A bit-sliced kernel of one lane width and field degree, behind a trait object.
pub trait BsRunner: Sync {
    fn pairs_total(&self) -> u64;
    fn run(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, u64);
}

struct Runner<L: Dispatch, const N: usize> {
    k: BsKernel<L, N>,
    threads: usize,
}

impl<L: Dispatch, const N: usize> BsRunner for Runner<L, N> {
    fn pairs_total(&self) -> u64 {
        self.k.pairs_total()
    }
    fn run(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, u64) {
        self.k.run(r, all, self.threads, out)
    }
}

macro_rules! runner_for_n {
    ($lane:ty, $rung:expr, $solver:expr, $threads:expr) => {{
        let (gf, l) = (&$rung.gf, $rung.l);
        match $rung.n {
            24 => Box::new(Runner::<$lane, 24> {
                k: BsKernel::new(gf, &$solver, l),
                threads: $threads,
            }) as Box<dyn BsRunner>,
            27 => Box::new(Runner::<$lane, 27> {
                k: BsKernel::new(gf, &$solver, l),
                threads: $threads,
            }),
            30 => Box::new(Runner::<$lane, 30> {
                k: BsKernel::new(gf, &$solver, l),
                threads: $threads,
            }),
            33 => Box::new(Runner::<$lane, 33> {
                k: BsKernel::new(gf, &$solver, l),
                threads: $threads,
            }),
            36 => Box::new(Runner::<$lane, 36> {
                k: BsKernel::new(gf, &$solver, l),
                threads: $threads,
            }),
            39 => Box::new(Runner::<$lane, 39> {
                k: BsKernel::new(gf, &$solver, l),
                threads: $threads,
            }),
            n => panic!("no bit-sliced instantiation for n = {n}"),
        }
    }};
}

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
    /// The bit-sliced group-law kernel: `bs64`, `bs256` or `bs512`, `-tK` threads.
    BitSliced(&'a Rung, Box<dyn BsRunner>),
    /// PS-8: the meet-in-the-middle pair table, probed once per signed base point.
    Mitm(&'a Rung, MitmTable<'a>),
}

impl<'a> Arm<'a> {
    pub fn build(name: &str, rung: &'a Rung) -> Self {
        if let Some(rest) = name.strip_prefix("bs") {
            let (width, threads) = match rest.split_once("-t") {
                Some((w, t)) => (w, t.parse::<usize>().expect("threads")),
                None => (rest, 1),
            };
            let solver = ArtinSchreier::new(&rung.gf);
            let runner: Box<dyn BsRunner> = match width {
                "64" => runner_for_n!(u64, rung, solver, threads),
                #[cfg(target_arch = "x86_64")]
                "256" => runner_for_n!(crate::lanes_x86::W256, rung, solver, threads),
                #[cfg(target_arch = "x86_64")]
                "512" => runner_for_n!(crate::lanes_x86::W512, rung, solver, threads),
                other => panic!("unknown bit-sliced width {other}"),
            };
            return Arm::BitSliced(rung, runner);
        }
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
            "mitm" => {
                let solver = ArtinSchreier::new(&rung.gf);
                Arm::Mitm(rung, MitmTable::new(&rung.gf, &solver, rung.l))
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
            Arm::BitSliced(_, k) => k.pairs_total(),
            // The table arm's unit is a probe, not a pair.
            Arm::Mitm(_, t) => t.probes_per_target(),
        }
    }

    /// Keys and bytes of a per-run table the arm built (zero for the arms without one).
    pub fn aux(&self) -> (u64, u64) {
        match self {
            Arm::Mitm(_, t) => (t.entries, t.bytes()),
            _ => (0, 0),
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
            Arm::BitSliced(r, k) => {
                if degenerate_target(&r.gf, r.l, xr) {
                    decompose(xr, r.l, &r.gf).is_some()
                } else {
                    let mut out = Vec::new();
                    k.run((xr, yr), false, &mut out).0
                }
            }
            Arm::Mitm(r, t) => {
                if degenerate_target(&r.gf, r.l, xr) {
                    decompose(xr, r.l, &r.gf).is_some()
                } else {
                    let mut out = Vec::new();
                    t.run((xr, yr), false, &mut out).0
                }
            }
        }
    }

    /// Every sorted root, for the arms that can list them (the group kernels).
    pub fn all_roots(&self, xr: u64, yr: u64) -> Option<Vec<[u64; 3]>> {
        let list = |k: &dyn Fn(&mut Vec<[u64; 3]>)| {
            let mut out = Vec::new();
            k(&mut out);
            out.sort_unstable();
            out.dedup();
            Some(out)
        };
        match self {
            Arm::GroupScalar(r, k) if !degenerate_target(&r.gf, r.l, xr) => list(&|o| {
                k.run((xr, yr), true, o);
            }),
            Arm::BitSliced(r, k) if !degenerate_target(&r.gf, r.l, xr) => list(&|o| {
                k.run((xr, yr), true, o);
            }),
            Arm::Mitm(r, t) if !degenerate_target(&r.gf, r.l, xr) => list(&|o| {
                t.run((xr, yr), true, o);
            }),
            _ => None,
        }
    }
}
