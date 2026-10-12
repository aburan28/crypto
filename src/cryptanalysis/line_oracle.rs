//! # An `O(log p)` decomposition oracle for a line base over `F_{p^k}`.
//!
//! A base `{P : x(P) ∈ c·F_p}` on `E(F_{p^k})` has one `F_p`-unknown
//! per point, so the two-summand point-decomposition problem for a
//! target `R` is the summation polynomial
//!
//! ```text
//!     S₃(c t₁, c t₂, x_R) = 0,        t₁, t₂ ∈ F_p,
//! ```
//!
//! a polynomial of bidegree `(2, 2)` over `F_{p^k}` whose `k` `F_p`-
//! coordinates are `k` equations over `F_p` in two unknowns (the Weil
//! descent).  Two of them have a resultant in `t₂` of degree at most
//! `8` in `t₁`; its roots in `F_p` come from `gcd(t^p − t, ·)`, each
//! root's `t₂` from the common roots of the descended quadratics, and
//! the pair lifts to signed base points through the framework's own
//! `lift_abscissae`.  The whole thing costs `O(log p)` field operations
//! a target where the `subtract` oracle costs `|F| ≈ p` group
//! additions, and it is **complete**: `S₃` vanishes exactly when the
//! three abscissae carry points summing to zero over the algebraic
//! closure, so every decomposition over the base is found, and a root
//! whose points live on the twist (`y ∉ F_{p^k}`) fails to lift and is
//! counted as such, never turned into a relation.
//!
//! The coefficients are extracted by evaluation and interpolation at
//! `t ∈ {0, 1, 2}` (nine evaluations of `S₃`), and the resultant by
//! evaluating the `4 × 4` Sylvester determinant at nine values of `t₁`
//! and interpolating — no symbolic polynomial arithmetic over the
//! extension field is needed, which keeps the oracle generic in the
//! field.  Its cost is reported through `solver_totals` in exact `F_p`
//! multiplications of the descent, resultant and root-finding stages.
//!
//! Used by E2b (the GLS line over `F_{p²}`) and E5 (the Frobenius line
//! over `F_{p³}`) of `RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`.

use std::time::Instant;

use crate::cryptanalysis::ext_curve::{ExtCurve, ExtField, ExtPoint};
use crate::cryptanalysis::glv_invariant_base::{addm, mulm, poly_roots_fp_counted, subm};
use crate::cryptanalysis::ic_boundary::{
    lift_abscissae, CountedGroup, FactorBase, GroupOps, OracleCounters,
};
use crate::cryptanalysis::ic_framework::stages::{
    DecompositionOracle, InstanceCtx, Params, SolverCost, SolverTotals, SystemShape,
};
use crate::cryptanalysis::residual_walk::inv_mod;

/// Semaev's third summation polynomial for `y² = x³ + ax + b`:
/// `(x₁ − x₂)² x₃² − 2((x₁ + x₂)(x₁x₂ + a) + 2b) x₃ + (x₁x₂ − a)² − 4b(x₁ + x₂)`.
pub fn semaev_s3<F: ExtField>(f: &F, a: F::El, b: F::El, x1: F::El, x2: F::El, x3: F::El) -> F::El {
    let d = f.sub(x1, x2);
    let sum = f.add(x1, x2);
    let prod = f.mul(x1, x2);
    let term1 = f.mul(f.sqr(d), f.sqr(x3));
    let inner = f.add(f.mul(sum, f.add(prod, a)), f.scale(b, 2));
    let term2 = f.mul(f.scale(inner, 2), x3);
    let term3 = f.sqr(f.sub(prod, a));
    let term4 = f.mul(f.scale(b, 4), sum);
    f.sub(f.add(f.sub(term1, term2), term3), term4)
}

/// Coefficients of the quadratic through `(0, v₀), (1, v₁), (2, v₂)`,
/// low to high, over the field.
fn quadratic_from_values<F: ExtField>(f: &F, v: [F::El; 3]) -> [F::El; 3] {
    let p = f.p();
    let half = inv_mod(2, p);
    // c₂ = (v₀ − 2v₁ + v₂)/2, c₁ = (−3v₀ + 4v₁ − v₂)/2, c₀ = v₀.
    let c2 = f.scale(f.add(f.sub(v[0], f.scale(v[1], 2)), v[2]), half);
    let c1 = f.scale(
        f.sub(f.add(f.scale(v[1], 4), f.neg(f.scale(v[0], 3))), v[2]),
        half,
    );
    [v[0], c1, c2]
}

/// Determinant of a `4 × 4` matrix over `F_p`, with a multiplication
/// counter.
fn det4(m: &[[u64; 4]; 4], p: u64, muls: &mut u64) -> u64 {
    let mut a = *m;
    let mut det = 1u64;
    for col in 0..4 {
        let Some(piv) = (col..4).find(|&r| a[r][col] != 0) else {
            return 0;
        };
        if piv != col {
            a.swap(piv, col);
            det = subm(0, det, p);
        }
        det = mulm(det, a[col][col], p);
        *muls += 1;
        let inv = inv_mod(a[col][col], p);
        for r in (col + 1)..4 {
            if a[r][col] == 0 {
                continue;
            }
            let factor = mulm(a[r][col], inv, p);
            *muls += 1;
            for c in col..4 {
                a[r][c] = subm(a[r][c], mulm(factor, a[col][c], p), p);
                *muls += 1;
            }
        }
    }
    det
}

/// Newton interpolation through `(0, v₀), …, (n − 1, v_{n−1})`, low to
/// high.
fn interpolate(values: &[u64], p: u64, muls: &mut u64) -> Vec<u64> {
    let n = values.len();
    // Divided differences.
    let mut dd: Vec<u64> = values.to_vec();
    let mut coef = vec![dd[0]];
    for level in 1..n {
        let inv = inv_mod(level as u64, p);
        for i in 0..(n - level) {
            dd[i] = mulm(subm(dd[i + 1], dd[i], p), inv, p);
            *muls += 1;
        }
        coef.push(dd[0]);
    }
    // Expand Σ coef_k Π_{j<k} (t − j).
    let mut out = vec![0u64; n];
    let mut basis = vec![1u64]; // Π (t − j)
    for (k, &c) in coef.iter().enumerate() {
        for (i, &bi) in basis.iter().enumerate() {
            out[i] = addm(out[i], mulm(c, bi, p), p);
            *muls += 1;
        }
        // basis *= (t − k)
        let mut next = vec![0u64; basis.len() + 1];
        for (i, &bi) in basis.iter().enumerate() {
            next[i + 1] = addm(next[i + 1], bi, p);
            next[i] = subm(next[i], mulm(bi, k as u64, p), p);
            *muls += 1;
        }
        basis = next;
    }
    while out.last() == Some(&0) {
        out.pop();
    }
    out
}

/// Evaluate a polynomial (low to high) at `t`.
fn eval_poly(c: &[u64], t: u64, p: u64, muls: &mut u64) -> u64 {
    let mut acc = 0u64;
    for &ci in c.iter().rev() {
        acc = addm(mulm(acc, t, p), ci, p);
        *muls += 1;
    }
    acc
}

/// The oracle: `c` is the line's direction.
pub struct LineOracle<F: ExtField> {
    f: F,
    c: F::El,
    totals: SolverTotals,
    shape: SystemShape,
    /// Systems with roots in `F_p` none of whose pairs lifted to base
    /// points summing to the target (twist abscissae, or `t = 0`).
    pub unliftable: u64,
    /// Targets whose resultant vanished identically, so the descent
    /// could not isolate `t₁`; never turned into a relation.
    pub degenerate: u64,
}

impl<F: ExtField> LineOracle<F> {
    pub fn new(f: F, c: F::El) -> Self {
        let k = f.degree() as usize;
        Self {
            f,
            c,
            totals: SolverTotals {
                solver: "line-resultant".into(),
                ..Default::default()
            },
            shape: SystemShape {
                n_vars: 2,
                n_equations: k,
                degrees: vec![4; k],
                semi_regular_degree: None,
            },
            unliftable: 0,
            degenerate: 0,
        }
    }

    /// The `k` descended polynomials `f_j(t₁, t₂)`, each as `[i][j]`
    /// coefficients of `t₁ⁱ t₂ʲ`.
    fn descend(&self, curve: &ExtCurve<F>, x_r: F::El) -> Vec<[[u64; 3]; 3]> {
        let f = &self.f;
        // g(t₁, t₂) = S₃(c t₁, c t₂, x_R); evaluate on the 3 × 3 grid.
        let mut grid = [[f.zero(); 3]; 3];
        for (i, row) in grid.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                let x1 = f.scale(self.c, i as u64);
                let x2 = f.scale(self.c, j as u64);
                *cell = semaev_s3(f, curve.a, curve.b, x1, x2, x_r);
            }
        }
        // Interpolate in t₁ for each t₂ = j, then in t₂.
        let mut by_j = [[f.zero(); 3]; 3]; // by_j[j][i]: coefficient of t₁ⁱ at t₂ = j
        for j in 0..3 {
            by_j[j] = quadratic_from_values(f, [grid[0][j], grid[1][j], grid[2][j]]);
        }
        let mut m = [[f.zero(); 3]; 3]; // m[i][j]: coefficient of t₁ⁱ t₂ʲ
        for i in 0..3 {
            let q = quadratic_from_values(f, [by_j[0][i], by_j[1][i], by_j[2][i]]);
            m[i] = q;
        }
        let k = f.degree() as usize;
        let mut out = vec![[[0u64; 3]; 3]; k];
        for i in 0..3 {
            for j in 0..3 {
                let coords = f.coords(m[i][j]);
                for (comp, &v) in coords.iter().enumerate() {
                    out[comp][i][j] = v;
                }
            }
        }
        out
    }

    /// The pairs `(t₁, t₂)` in `F_p²` at which every descended
    /// polynomial vanishes.
    fn solve(&mut self, polys: &[[[u64; 3]; 3]], seed: u64) -> (Vec<(u64, u64)>, u64, bool) {
        let p = self.f.p();
        let mut muls = 0u64;
        // Coefficients of t₂ʲ as polynomials in t₁, for the first two components.
        let coef_in_t1 =
            |q: &[[u64; 3]; 3], j: usize| -> Vec<u64> { vec![q[0][j], q[1][j], q[2][j]] };
        // The first pair of components with a non-zero resultant; a
        // component that vanishes identically (a `j = 0` curve on a
        // line with `x³ ∈ F_p` loses one) is skipped, and a system with
        // no such pair is degenerate.
        let nonzero: Vec<&[[u64; 3]; 3]> = polys
            .iter()
            .filter(|q| q.iter().any(|row| row.iter().any(|&v| v != 0)))
            .collect();
        let mut res: Vec<u64> = Vec::new();
        'pairs: for (ia, f0) in nonzero.iter().enumerate() {
            for f1 in nonzero.iter().skip(ia + 1) {
                // Resultant in t₂ by evaluation at t₁ = 0..8 and interpolation.
                let mut values = Vec::with_capacity(9);
                for t1 in 0..9u64 {
                    let a: Vec<u64> = (0..3)
                        .map(|j| eval_poly(&coef_in_t1(f0, j), t1, p, &mut muls))
                        .collect();
                    let b: Vec<u64> = (0..3)
                        .map(|j| eval_poly(&coef_in_t1(f1, j), t1, p, &mut muls))
                        .collect();
                    // Sylvester matrix of a₂t² + a₁t + a₀ and b₂t² + b₁t + b₀.
                    let m = [
                        [a[2], a[1], a[0], 0],
                        [0, a[2], a[1], a[0]],
                        [b[2], b[1], b[0], 0],
                        [0, b[2], b[1], b[0]],
                    ];
                    values.push(det4(&m, p, &mut muls));
                }
                let candidate = interpolate(&values, p, &mut muls);
                if !candidate.is_empty() {
                    res = candidate;
                    break 'pairs;
                }
            }
        }
        if res.is_empty() {
            return (Vec::new(), muls, true);
        }
        let t1s = poly_roots_fp_counted(&res, p, seed, &mut muls);
        let mut pairs = Vec::new();
        for t1 in t1s {
            // Candidate t₂: roots of the first component that is not
            // identically zero at this t₁.
            let mut cands: Option<Vec<u64>> = None;
            for q in polys {
                let a: Vec<u64> = (0..3)
                    .map(|j| eval_poly(&coef_in_t1(q, j), t1, p, &mut muls))
                    .collect();
                if a.iter().all(|&v| v == 0) {
                    continue;
                }
                cands = Some(poly_roots_fp_counted(&a, p, seed ^ t1, &mut muls));
                break;
            }
            let Some(cands) = cands else {
                continue;
            };
            for t2 in cands {
                let all = polys.iter().all(|q| {
                    let a: Vec<u64> = (0..3)
                        .map(|j| eval_poly(&coef_in_t1(q, j), t1, p, &mut muls))
                        .collect();
                    eval_poly(&a, t2, p, &mut muls) == 0
                });
                if all {
                    pairs.push((t1, t2));
                }
            }
        }
        (pairs, muls, false)
    }
}

impl<F: ExtField> DecompositionOracle<ExtCurve<F>> for LineOracle<F> {
    fn name(&self) -> &str {
        "line-resultant"
    }

    fn summands(&self) -> u32 {
        2
    }

    fn describe(&self, _params: &Params) -> String {
        format!(
            "Weil-descend S₃(c t₁, c t₂, x_R) over the line c·F_p to {} equations in (t₁, t₂), resultant in t₂, roots by gcd(t^p − t, ·)",
            self.f.degree()
        )
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<ExtCurve<F>>,
        fb: &FactorBase<ExtPoint<F::El>>,
        _params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<(), String> {
        // Every base abscissa must lie on the line.
        let f = &ctx.group.f;
        let c_inv = f.inv(self.c).ok_or("the line direction is zero")?;
        for pt in &fb.points {
            let t = f.mul(pt.x, c_inv);
            let coords = f.coords(t);
            if coords[1..].iter().any(|&v| v != 0) {
                return Err("the line oracle needs a base whose abscissae lie on c·F_p".into());
            }
        }
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<ExtCurve<F>>,
        fb: &FactorBase<ExtPoint<F::El>>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: ExtPoint<F::El>,
    ) -> Option<Vec<usize>> {
        let started = Instant::now();
        let polys = self.descend(ctx.group, point.x);
        let seed = ctx.group.key(&point);
        let (pairs, muls, degenerate) = self.solve(&polys, seed);
        let cost = SolverCost {
            ops: muls,
            op_unit: "F_p multiplications".into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: Some(8),
            solving_degree: Some(8),
            timed_out: false,
            extra: Default::default(),
        };
        self.totals.absorb(&self.shape.clone(), Some(&cost), false);
        if degenerate {
            self.degenerate += 1;
            return None;
        }
        if pairs.is_empty() {
            return None;
        }
        let f = &ctx.group.f;
        for (t1, t2) in pairs {
            let xs = [f.pack(f.scale(self.c, t1)), f.pack(f.scale(self.c, t2))];
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
    use crate::cryptanalysis::ext_curve::{random_point, Fp2};
    use crate::cryptanalysis::gls_fp2::{generate_gls_instance, gls_line_base};
    use crate::cryptanalysis::ic_boundary::CountedGroup;
    use crate::cryptanalysis::ic_framework::plugins::SubtractOracle;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn s3_vanishes_on_the_abscissae_of_three_points_summing_to_zero() {
        let inst = generate_gls_instance(9, 4, 16).unwrap();
        let curve = &inst.curve;
        let mut rng = StdRng::seed_from_u64(1);
        let mut ops = GroupOps::default();
        for _ in 0..20 {
            let p = random_point(curve, &mut rng);
            let q = random_point(curve, &mut rng);
            let s = curve.add(&mut ops, p, q);
            assert_eq!(
                semaev_s3(&curve.f, curve.a, curve.b, p.x, q.x, s.x),
                curve.f.zero()
            );
            let t = random_point(curve, &mut rng);
            if t.x != s.x {
                assert_ne!(
                    semaev_s3(&curve.f, curve.a, curve.b, p.x, q.x, t.x),
                    curve.f.zero()
                );
            }
        }
    }

    #[test]
    fn interpolation_and_determinant_round_trip() {
        let p = 1009;
        let mut muls = 0;
        let poly = vec![5u64, 0, 7, 1, 0, 3, 0, 0, 2];
        let values: Vec<u64> = (0..9).map(|t| eval_poly(&poly, t, p, &mut muls)).collect();
        assert_eq!(interpolate(&values, p, &mut muls), poly);
        let m = [[2, 0, 0, 0], [0, 3, 0, 0], [0, 0, 5, 0], [0, 0, 0, 7]];
        assert_eq!(det4(&m, p, &mut muls), 210);
        let m = [[1, 2, 3, 4], [2, 4, 6, 8], [1, 0, 1, 0], [0, 1, 0, 1]];
        assert_eq!(det4(&m, p, &mut muls), 0);
    }

    /// **The oracle agrees with `subtract` target by target** (AGENTS.md
    /// §6): both are complete over the base, so they must find a
    /// decomposition for exactly the same targets.
    #[test]
    fn the_line_oracle_agrees_with_subtract_on_every_target() {
        let inst = generate_gls_instance(9, 5, 16).unwrap();
        let (fb, _) = gls_line_base(&inst, true).unwrap();
        let mut ops = GroupOps::default();
        let target = inst.curve.mul(&mut ops, inst.generator, 777 % inst.r);
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: inst.generator,
            target,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(2),
        };
        let mut line: LineOracle<Fp2> = LineOracle::new(inst.curve.f, inst.line);
        let mut subtract = SubtractOracle;
        line.prepare(&ctx, &fb, &Params::default(), &mut ops)
            .unwrap();
        let (mut ca, mut cb) = (OracleCounters::default(), OracleCounters::default());
        let mut hits = 0;
        for k in 1..=300u64 {
            let pt = inst.curve.mul(&mut ops, inst.generator, k);
            let a = line.decompose(&ctx, &fb, &mut ops, &mut ca, pt);
            let b = subtract.decompose(&ctx, &fb, &mut ops, &mut cb, pt);
            assert_eq!(
                a.is_some(),
                b.is_some(),
                "target [{k}]G: line {a:?}, subtract {b:?}"
            );
            if let Some(idx) = a {
                let sum = inst
                    .curve
                    .add(&mut ops, fb.points[idx[0]], fb.points[idx[1]]);
                assert_eq!(sum, pt, "a decomposition that does not sum to its target");
                hits += 1;
            }
        }
        assert!(hits > 20, "only {hits} hits in 300 targets");
        assert_eq!(line.degenerate, 0, "no degenerate resultant");
        let totals = line.solver_totals().unwrap();
        assert_eq!(totals.calls, 300);
        assert!(totals.ops > 0);
    }
}
