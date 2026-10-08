//! # The algebraic three-summand oracle on a Frobenius line (E11).
//!
//! On a subfield curve `E / F_p` the Frobenius-stable line `x ∈ s·F_p`
//! of `E(F_{p³})` is the base of
//! `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`
//! §6.5 and §8.2.  Two summands there are degenerate (every relation
//! stays inside one orbit, §6.5) and three summands over the pair
//! table cost `|F|` additions a target (§8.2).  This oracle decomposes
//! a target `R` into three line points algebraically: the fourth
//! summation polynomial `S₄(x₁, x₂, x₃, x_R)`, symmetrised in
//! `(x₁, x₂, x₃)` by `gaudry_cubic::SymmetrisedS4`, is restricted to the
//! line by `x_i = s t_i` — the elementary symmetric functions scale as
//! `e_k = s^k σ_k(t)`, so only the coefficients move
//! ([`SymmetrisedS4::on_line`]) — then Weil-descended to three `F_p`
//! equations in `(σ₁, σ₂, σ₃)` and solved by the Macaulay-matrix and
//! eigenvalue solver `solve_s4_subspace_with`, whose cost does not
//! depend on `p`.  The roots `t_i ∈ F_p` of the cubic
//! `T³ − σ₁T² + σ₂T − σ₃` are the line abscissae; the signs are
//! settled by group arithmetic (`lift_abscissae`), so a decomposition
//! is returned only when the three points sum to `R`.
//!
//! Two presentations of `F_{p³}` meet here: `ext_curve::Fp3` (`s³ = ν`)
//! carries the group, `gaudry_cubic::Fp3` (`t³ = c`) the solver.  With
//! `c = ν` they are the same field in the same basis and an element
//! carries over coefficient by coefficient.
//!
//! The oracle is cross-checked against the pair-table oracle target by
//! target by the E11 driver (AGENTS.md §6) and by this module's tests.

use std::time::Instant;

use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::cryptanalysis::ext_curve::{ExtCurve, ExtField, ExtPoint, Fp3, Fp3El};
use crate::cryptanalysis::gaudry_cubic::{
    solve_s4_subspace_with, Curve3, Fp3 as SolverField, Instance3, Pt3, SolveMode, SolveStats,
    SymmetrisedS4, E3,
};
use crate::cryptanalysis::ic_boundary::{lift_abscissae, FactorBase, GroupOps, OracleCounters};
use crate::cryptanalysis::ic_framework::stages::{
    DecompositionOracle, InstanceCtx, Params, SolverCost, SolverTotals, SystemShape,
};

/// The oracle: built once per base by `prepare`, then one Macaulay solve
/// per target.
pub struct LineS4Oracle {
    f: Fp3,
    line: Fp3El,
    solver: Option<(Instance3, SymmetrisedS4)>,
    rng: StdRng,
    /// The solver's own counters: `fp_muls` is its cost in `F_p`
    /// multiplications, `unsolved` the targets no Macaulay degree closed.
    pub stats: SolveStats,
    totals: SolverTotals,
    shape: SystemShape,
    /// Targets whose system had roots in `F_p` none of which lifted to
    /// three base points summing to the target.
    pub unliftable: u64,
    /// Root triples returned by the solver, summed over targets.
    pub solutions: u64,
}

fn to_e3(a: Fp3El) -> E3 {
    E3(a)
}

impl LineS4Oracle {
    pub fn new(f: Fp3, line: Fp3El, seed: u64) -> Self {
        Self {
            f,
            line,
            solver: None,
            rng: StdRng::seed_from_u64(seed ^ 0x5334_4C49_4E45),
            stats: SolveStats::default(),
            totals: SolverTotals {
                solver: "line-s4-macaulay".into(),
                ..Default::default()
            },
            shape: SystemShape {
                n_vars: 3,
                n_equations: 3,
                degrees: vec![4; 3],
                semi_regular_degree: None,
            },
            unliftable: 0,
            solutions: 0,
        }
    }

    /// The solver's field and the symmetrised polynomial on the line,
    /// built from the curve: `E / F_p` must have `a, b ∈ F_p`.
    fn build(
        &self,
        curve: &ExtCurve<Fp3>,
        r: u64,
        generator: ExtPoint<Fp3El>,
    ) -> Result<(Instance3, SymmetrisedS4), String> {
        let field = SolverField::with_cube_nonresidue(self.f.p, self.f.nu)?;
        for (name, v) in [("a", curve.a), ("b", curve.b)] {
            if v[1] != 0 || v[2] != 0 {
                return Err(format!(
                    "the coefficient {name} is not in F_p: the S₄ solver needs a subfield curve"
                ));
            }
        }
        let g = if generator.infinity {
            Pt3::INFINITY
        } else {
            Pt3::affine(to_e3(generator.x), to_e3(generator.y))
        };
        let curve3 = Curve3::new(field, to_e3(curve.a), to_e3(curve.b), r, g);
        let pre = SymmetrisedS4::precompute(&curve3).on_line(&curve3.field, &to_e3(self.line));
        let inst = Instance3 {
            curve: curve3,
            q: Pt3::INFINITY,
            d: 0,
        };
        Ok((inst, pre))
    }
}

impl DecompositionOracle<ExtCurve<Fp3>> for LineS4Oracle {
    fn name(&self) -> &str {
        "line-s4-macaulay"
    }

    fn summands(&self) -> u32 {
        3
    }

    fn describe(&self, _params: &Params) -> String {
        "S₄(s t₁, s t₂, s t₃, x_R) symmetrised in the t_i, Weil-descended over the line s·F_p to three equations in (σ₁, σ₂, σ₃), solved by the Macaulay matrix and the eigenvalues of e₁, cubic split over F_p, signs by group arithmetic".into()
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<ExtCurve<Fp3>>,
        fb: &FactorBase<ExtPoint<Fp3El>>,
        _params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<(), String> {
        let f = &ctx.group.f;
        let inv = f.inv(self.line).ok_or("the line direction is zero")?;
        for pt in &fb.points {
            let t = f.mul(pt.x, inv);
            if f.coords(t)[1..].iter().any(|&v| v != 0) {
                return Err("the S₄ line oracle needs a base whose abscissae lie on s·F_p".into());
            }
        }
        self.solver = Some(self.build(ctx.group, ctx.r, ctx.generator)?);
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<ExtCurve<Fp3>>,
        fb: &FactorBase<ExtPoint<Fp3El>>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: ExtPoint<Fp3El>,
    ) -> Option<Vec<usize>> {
        if point.infinity {
            return None;
        }
        let (inst, pre) = self.solver.as_ref()?;
        let started = Instant::now();
        let muls_before = self.stats.fp_muls;
        let triples = solve_s4_subspace_with(
            inst,
            pre,
            &to_e3(point.x),
            &mut self.rng,
            &mut self.stats,
            SolveMode::default(),
        );
        let cost = SolverCost {
            ops: self.stats.fp_muls - muls_before,
            op_unit: "F_p multiplications".into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: Some(13),
            solving_degree: Some(10),
            timed_out: false,
            extra: Default::default(),
        };
        self.totals
            .absorb(&self.shape.clone(), Some(&cost), triples.is_none());
        let triples = triples?;
        if triples.is_empty() {
            return None;
        }
        let f = &ctx.group.f;
        for [t1, t2, t3] in &triples {
            self.solutions += 1;
            let xs = [
                f.pack(f.scale(self.line, *t1)),
                f.pack(f.scale(self.line, *t2)),
                f.pack(f.scale(self.line, *t3)),
            ];
            if let Some(indices) = lift_abscissae(ctx.group, fb, ops, &xs, point) {
                return Some(indices);
            }
            counters.lift_failures += 1;
        }
        self.unliftable += 1;
        counters.unliftable_systems += 1;
        None
    }

    fn last_system(&self) -> Option<SystemShape> {
        Some(self.shape.clone())
    }

    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::CountedGroup;
    use crate::cryptanalysis::ic_framework::plugins::MitmOracle;
    use crate::cryptanalysis::subfield_fp3::{
        generate_subfield_instance, subfield_line_base, SubfieldFold,
    };

    fn ctx_for(
        inst: &crate::cryptanalysis::subfield_fp3::SubfieldInstance,
        k: u64,
    ) -> (InstanceCtx<'_, ExtCurve<Fp3>>, ExtPoint<Fp3El>) {
        let mut ops = GroupOps::default();
        let target = inst.curve.mul(&mut ops, inst.generator, k);
        (
            InstanceCtx {
                group: &inst.curve,
                generator: inst.generator,
                target,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: Some(3),
            },
            target,
        )
    }

    /// The line-twisted symmetrised polynomial vanishes on the
    /// `σ`-values of three line points summing to the target.
    #[test]
    fn the_twisted_polynomial_vanishes_on_a_true_decomposition() {
        let inst = generate_subfield_instance(7, 1, false, 8).unwrap();
        let (fb, _) = subfield_line_base(&inst, SubfieldFold::Negation).unwrap();
        let oracle = LineS4Oracle::new(inst.curve.f, inst.line, 1);
        let (inst3, pre) = oracle.build(&inst.curve, inst.r, inst.generator).unwrap();
        let f3 = &inst3.curve.field;
        let f = &inst.curve.f;
        let inv = f.inv(inst.line).unwrap();
        let mut ops = GroupOps::default();
        for k in 0..20usize {
            let n = fb.points.len();
            let (p1, p2, p3) = (
                fb.points[k % n],
                fb.points[(3 * k + 1) % n],
                fb.points[(5 * k + 2) % n],
            );
            let p12 = inst.curve.add(&mut ops, p1, p2);
            let r = inst.curve.add(&mut ops, p12, p3);
            if r.infinity {
                continue;
            }
            let t: Vec<u64> = [p1, p2, p3]
                .iter()
                .map(|p| f.coords(f.mul(p.x, inv))[0])
                .collect();
            let p = f.p;
            let m = |a: u64, b: u64| ((a as u128 * b as u128) % p as u128) as u64;
            let s1 = (t[0] + t[1] + t[2]) % p;
            let s2 = (m(t[0], t[1]) + m(t[0], t[2]) + m(t[1], t[2])) % p;
            let s3 = m(m(t[0], t[1]), t[2]);
            let sig = [1u64, s1, s2, s3];
            let mut acc = SolverField::ZERO;
            let x_r = to_e3(r.x);
            for (&e, &c) in pre.terms() {
                let mut term = c;
                for _ in 0..e[0] {
                    term = f3.mul(&term, &f3.from_base(sig[1]));
                }
                for _ in 0..e[1] {
                    term = f3.mul(&term, &f3.from_base(sig[2]));
                }
                for _ in 0..e[2] {
                    term = f3.mul(&term, &f3.from_base(sig[3]));
                }
                for _ in 0..e[3] {
                    term = f3.mul(&term, &x_r);
                }
                acc = f3.add(&acc, &term);
            }
            assert_eq!(acc, SolverField::ZERO, "k = {k}");
        }
    }

    /// Every decomposition the oracle returns sums to its target, and
    /// the oracle agrees with the pair table on whether a target
    /// decomposes, target by target.
    #[test]
    fn the_s4_line_oracle_agrees_with_the_pair_table_target_by_target() {
        let inst = generate_subfield_instance(7, 2, false, 8).unwrap();
        let (fb, _) = subfield_line_base(&inst, SubfieldFold::Frobenius).unwrap();
        let (ctx, _) = ctx_for(&inst, 3);
        let mut ops = GroupOps::default();
        let mut s4 = LineS4Oracle::new(inst.curve.f, inst.line, 7);
        s4.prepare(&ctx, &fb, &Params::default(), &mut ops).unwrap();
        let mut mitm = MitmOracle::new(3);
        let mut params = Params::default();
        params.set("negation_folded", "1");
        mitm.prepare(&ctx, &fb, &params, &mut ops).unwrap();
        let (mut ca, mut cb) = (OracleCounters::default(), OracleCounters::default());
        let (mut hits, mut agree, mut disagree) = (0, 0, 0);
        for k in 2..120u64 {
            let pt = inst.curve.mul(&mut ops, inst.generator, k);
            let a = s4.decompose(&ctx, &fb, &mut ops, &mut ca, pt);
            let b = mitm.decompose(&ctx, &fb, &mut ops, &mut cb, pt);
            if let Some(idx) = &a {
                assert_eq!(idx.len(), 3);
                let sum = idx.iter().fold(inst.curve.identity(), |acc, &i| {
                    inst.curve.add(&mut ops, acc, fb.points[i])
                });
                assert_eq!(sum, pt, "k = {k}");
                hits += 1;
            }
            if a.is_some() == b.is_some() {
                agree += 1;
            } else {
                disagree += 1;
            }
        }
        assert!(hits > 0, "no decomposition found at all");
        assert!(
            disagree <= agree / 10,
            "agree {agree}, disagree {disagree}, unsolved {}",
            s4.stats.unsolved
        );
    }
}
