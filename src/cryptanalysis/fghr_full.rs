//! # FGHR's 2-torsion symmetry on the full group `E(F_{p³})` (E17–E18).
//!
//! E13 (`fghr_line`) measured the cheapest decomposition oracle of
//! `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`
//! — `S₄` in a coordinate `Y` with `Y(P + T) = −Y(P)`, symmetrised by
//! the even sign changes and permutations of `(Y₁, Y₂, Y₃)` and solved
//! as three conics over `F_p[q₃]` — on a *subfield* curve, where §8.10
//! shows the large prime lives in a group of order `≈ p²` and no line
//! base can beat rho by more than a constant.  This module moves the
//! same oracle to the setting where index calculus has a sub-rho
//! exponent: a curve over `F_{p³}` that is **not** a subfield curve,
//! with a large prime `r ≈ p³ / h` for a small cofactor `h`, and
//! Gaudry's base of abscissae in an affine `F_p`-line (`n^{1/3}`
//! relations against rho's `n^{1/2}`;
//! `RESEARCH_EXTENSION_FIELD_BOUNDARIES.md` Theorems 2–3).
//!
//! * **The curve.**  `x₀ ∈ F_{p³} \ F_p` and `c ∈ F_p^*` at random,
//!   `a = c² − 3x₀²`, `b = −x₀³ − a x₀`: then `x₀` is a root of
//!   `x³ + ax + b` and `f'(x₀) = 3x₀² + a = c²`, so `T = (x₀, 0)` is a
//!   rational 2-torsion point and `x(P + T) − x₀ = c² / (x − x₀)`.
//!   Instances whose `j`-invariant lies in `F_p` (twists of subfield
//!   curves) are refused.  `#E(F_{p³})` is the unique multiple of a
//!   random point's order in the Hasse interval (baby-step giant-step);
//!   it is even, and the instance keeps it when `#E = h·r` with `r`
//!   prime and `h ≤ max_cofactor`.
//! * **The coordinate.**  `Y = (x − x₀ − c)/(x − x₀ + c)`, so
//!   `Y(P + T) = −Y(P)` and `Y(−P) = Y(P)`, and
//!   `x = ((c − x₀) Y + (x₀ + c)) / (1 − Y)`: a Möbius map with
//!   coefficients in `F_{p³}` (in `F_p` on E13's subfield curves).
//! * **The base.**  `{P : Y(P) ∈ F_p \ {0, ±1}}`, i.e. `x(P) ∈ x₀ +
//!   F_p` up to the Möbius map's two excluded values: Gaudry's affine
//!   `F_p`-line of abscissae, stable under negation and under
//!   `τ_T: P ↦ P + T` (`Y ↦ −Y`).  `T` lies in the cofactor, so `P` and
//!   `P + T` share their `⟨G⟩`-component and `τ_T` folds with
//!   coefficient `1`: `4` points a column under `⟨−1, τ_T⟩` against `2`
//!   under `⟨−1⟩`.
//! * **The oracles.**  `fghr_line`'s `D₃` conic-resultant oracle and its
//!   `S₃` Macaulay oracle, unchanged, on the line `Y ∈ F_p` (`L = 1`);
//!   the instance enters them only through [`YLine`].

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use crate::cryptanalysis::ext_curve::{random_point, ExtCurve, ExtField, ExtPoint, Fp3, Fp3El};
use crate::cryptanalysis::fghr_line::{Translation, YLine};
use crate::cryptanalysis::gaudry_cubic::{
    unique_hasse_multiple3, Curve3, Fp3 as SolverField, Pt3, E3,
};
use crate::cryptanalysis::glv_invariant_base::{
    factor_u64, fold_by_endomorphisms, Closure, Endomorphism, FoldReport, Negation,
};
use crate::cryptanalysis::ic_boundary::{CountedGroup, FactorBase, GroupOps};
use crate::cryptanalysis::residual_walk::is_prime_u64;
use crate::cryptanalysis::subfield_fp3::{Fp3Curve, Fp3Point};

/// A non-subfield curve over `F_{p³}` with a rational 2-torsion point
/// `T = (x₀, 0)`, `f'(x₀) = c²`, and a large prime `r` of small cofactor.
#[derive(Clone, Debug, Serialize)]
pub struct FullFghrInstance {
    pub name: String,
    pub p: u64,
    pub curve: Fp3Curve,
    /// `#E(F_{p³}) = h·r`.
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    pub generator: Fp3Point,
    pub x0: Fp3El,
    pub c: u64,
    pub t: Fp3Point,
    /// `x = (A Y + B)/(C Y + D)`: `A = c − x₀`, `B = x₀ + c`, `C = −1`,
    /// `D = 1`.
    pub mobius: [Fp3El; 4],
    /// `j(E)`, checked not to lie in `F_p`.
    pub j: Fp3El,
}

impl FullFghrInstance {
    /// `x` from `Y`.
    pub fn x_of_y(&self, y: Fp3El) -> Fp3El {
        let f = &self.curve.f;
        let [a, b, c, d] = self.mobius;
        let num = f.add(f.mul(a, y), b);
        let den = f.add(f.mul(c, y), d);
        f.mul(num, f.inv(den).expect("Y = 1 is not on the line"))
    }

    /// `Y` from `x`.
    pub fn y_of_x(&self, x: Fp3El) -> Fp3El {
        let f = &self.curve.f;
        let u = f.sub(x, self.x0);
        let c = f.from_fp(self.c);
        let num = f.sub(u, c);
        let den = f.add(u, c);
        f.mul(num, f.inv(den).expect("x = x₀ − c has no Y"))
    }
}

impl YLine for FullFghrInstance {
    fn x_at(&self, t: u64) -> Fp3El {
        self.x_of_y(self.curve.f.from_fp(t))
    }
}

fn in_fp(a: Fp3El) -> bool {
    a[1] == 0 && a[2] == 0
}

/// A full-group instance at `p_bits`: a prime `p ≡ 1 (mod 3)` of that
/// size, `#E(F_{p³}) = h·r` with `r` prime and `h ≤ max_cofactor`
/// (`h` is even: `T` is rational).  Deterministic in `seed`.
pub fn generate_full_fghr_instance(
    p_bits: u32,
    seed: u64,
    max_cofactor: u64,
) -> Result<FullFghrInstance, String> {
    if !(5..=13).contains(&p_bits) {
        return Err(format!("p_bits = {p_bits} outside 5..=13"));
    }
    let mut rng = StdRng::seed_from_u64(seed ^ 0x4647_4852_4655_4C4C ^ (p_bits as u64));
    for _ in 0..100_000u32 {
        let p = loop {
            let c = rng.gen_range((1u64 << (p_bits - 1))..(1u64 << p_bits)) | 1;
            if is_prime_u64(c) && c % 3 == 1 {
                break c;
            }
        };
        let f = Fp3::new(p)?;
        let x0 = f.from_coords(&[
            rng.gen_range(0..p),
            rng.gen_range(0..p),
            rng.gen_range(0..p),
        ]);
        if in_fp(x0) {
            continue;
        }
        let c = rng.gen_range(1..p);
        let x0sq = f.sqr(x0);
        let a = f.sub(f.from_fp(c * c % p), f.scale(x0sq, 3));
        let b = f.neg(f.add(f.mul(x0sq, x0), f.mul(a, x0)));
        // Δ = −16(4a³ + 27b²) ≠ 0 and j = 1728·4a³/(4a³ + 27b²) ∉ F_p.
        let a3x4 = f.scale(f.mul(f.sqr(a), a), 4);
        let disc = f.add(a3x4, f.scale(f.sqr(b), 27));
        let Some(dinv) = f.inv(disc) else {
            continue;
        };
        let j = f.scale(f.mul(a3x4, dinv), 1728 % p);
        if in_fp(j) {
            continue;
        }
        let curve: Fp3Curve = ExtCurve { f, a, b };
        let field = SolverField::with_cube_nonresidue(p, f.nu)?;
        let c3 = Curve3::new(field, E3(a), E3(b), 0, Pt3::INFINITY);
        let pt = random_point(&curve, &mut rng);
        let Some(order) = unique_hasse_multiple3(&c3, &Pt3::affine(E3(pt.x), E3(pt.y))) else {
            continue;
        };
        if order % 2 != 0 {
            return Err(format!("#E = {order} is odd although T = (x₀, 0) is rational"));
        }
        let Some(&(r, _)) = factor_u64(order).last() else {
            continue;
        };
        let h = order / r;
        if h > max_cofactor || r < 1 << 10 || order % (r * r) == 0 {
            continue;
        }
        let mut ops = GroupOps::default();
        let generator = loop {
            let g = curve.mul(&mut ops, random_point(&curve, &mut rng), h);
            if !g.infinity {
                break g;
            }
        };
        if !curve.mul(&mut ops, generator, r).infinity {
            return Err("r·G ≠ O: the order is wrong".into());
        }
        let t = ExtPoint::affine(x0, f.zero());
        if !curve.is_on_curve(t) {
            return Err("T is not on the curve".into());
        }
        let mobius = [
            f.sub(f.from_fp(c), x0),
            f.add(x0, f.from_fp(c)),
            f.from_fp(p - 1),
            f.one(),
        ];
        let tag = rng.gen::<u32>();
        let inst = FullFghrInstance {
            name: format!("fghrfull-p{p_bits}bit-p{p}-h{h}-{tag:x}"),
            p,
            curve,
            group_order: order,
            r,
            cofactor: h,
            generator,
            x0,
            c,
            t,
            mobius,
            j,
        };
        // The translation acts as Y ↦ −Y: checked on points, not assumed.
        for _ in 0..8 {
            let q = random_point(&inst.curve, &mut rng);
            let s = inst.curve.add(&mut ops, q, t);
            let den = f.add(f.sub(s.x, x0), f.from_fp(c));
            if s.infinity || f.is_zero(den) || f.is_zero(f.add(f.sub(q.x, x0), f.from_fp(c))) {
                continue;
            }
            if f.add(inst.y_of_x(q.x), inst.y_of_x(s.x)) != f.zero() {
                return Err("Y(P + T) ≠ −Y(P): the coordinate is wrong".into());
            }
            if inst.x_of_y(inst.y_of_x(q.x)) != q.x {
                return Err("x(Y(x)) ≠ x: the Möbius map is wrong".into());
            }
        }
        return Ok(inst);
    }
    Err(format!(
        "no full-group instance with cofactor ≤ {max_cofactor} at p_bits = {p_bits}"
    ))
}

/// Which group the full-group `Y`-line base is folded by.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum FullFold {
    /// `⟨−1⟩`: `2` points a column, `P` and `P + T` in two columns.
    Negation,
    /// `⟨−1, τ_T⟩`: `4` points a column.
    Translation,
}

impl FullFold {
    pub fn name(self) -> &'static str {
        match self {
            Self::Negation => "negation",
            Self::Translation => "negation+translation",
        }
    }
}

/// The base `{P : Y(P) ∈ F_p \ {0, ±1}}`, folded by `fold`; the set is
/// checked to be invariant under the group, not assumed.  The same
/// points in the same order for every fold.
pub fn full_line_base(
    inst: &FullFghrInstance,
    fold: FullFold,
) -> Result<(FactorBase<Fp3Point>, FoldReport), String> {
    let start = std::time::Instant::now();
    let curve = &inst.curve;
    let f = &curve.f;
    let mut seed = Vec::with_capacity(2 * inst.p as usize);
    for t in 2..inst.p - 1 {
        seed.extend(curve.lift_x(inst.x_at(t)));
    }
    let neg = Negation { r: inst.r };
    let tau = Translation { t: inst.t };
    let mut gens: Vec<&dyn Endomorphism<Fp3Curve>> = vec![&neg];
    if fold == FullFold::Translation {
        gens.push(&tau);
    }
    let (mut fb, report) = fold_by_endomorphisms(
        curve,
        inst.r,
        inst.cofactor,
        seed,
        &gens,
        Closure::Strict,
        32,
        |p| curve.key(p),
        |p| f.pack(p.x),
        format!(
            "Gaudry's affine F_p-line x ∈ x0 + F_p on a non-subfield curve, as Y ∈ F_p \\ {{0, ±1}} ({} values), Y = (x − x0 − c)/(x − x0 + c), folded by {}",
            inst.p - 3,
            fold.name()
        ),
    )?;
    fb.cost.wall_ns = start.elapsed().as_nanos() as u64;
    fb.cost.count("sqrt_solves", inst.p - 3);
    Ok((fb, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::fghr_line::{
        fghr_polynomials_on, FghrOracle, Mobius, YLineS4Oracle,
    };
    use crate::cryptanalysis::ic_boundary::OracleCounters;
    use crate::cryptanalysis::ic_framework::plugins::MitmOracle;
    use crate::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx, Params};

    #[test]
    fn the_instance_is_a_non_subfield_curve_with_a_folding_translation() {
        let inst = generate_full_fghr_instance(7, 3, 8).unwrap();
        assert!(!in_fp(inst.j));
        assert_eq!(inst.group_order, inst.cofactor * inst.r);
        assert!(inst.r > inst.p * inst.p, "r ≈ p³/h, not p²");
        let (b4, r4) = full_line_base(&inst, FullFold::Translation).unwrap();
        let (b2, _) = full_line_base(&inst, FullFold::Negation).unwrap();
        assert_eq!(b4.points, b2.points);
        assert_eq!(r4.points_per_orbit, 4.0);
        assert!(b4.columns * 2 <= b2.columns + 2);
    }

    #[test]
    fn the_d3_oracle_agrees_with_the_s3_oracle_and_the_pair_table() {
        let inst = generate_full_fghr_instance(6, 1, 8).unwrap();
        let (fb, _) = full_line_base(&inst, FullFold::Translation).unwrap();
        let polys = fghr_polynomials_on(
            &inst.curve,
            inst.r,
            inst.generator,
            inst.curve.f.one(),
            Mobius::Fp3(inst.mobius),
        )
        .unwrap();
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: inst.generator,
            target: inst.generator,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(3),
        };
        let mut d3 = FghrOracle::new(&inst, &polys, 7);
        let mut s3 = YLineS4Oracle::new(&inst, &polys, 7);
        let mut mitm = MitmOracle::new(3);
        let mut params = Params::default();
        params.set("negation_folded", "1");
        let mut ops = GroupOps::default();
        mitm.prepare(&ctx, &fb, &params, &mut ops).unwrap();
        let mut c = [
            OracleCounters::default(),
            OracleCounters::default(),
            OracleCounters::default(),
        ];
        let (mut hits, mut found) = (0u32, 0u32);
        for k in 2..120u64 {
            let pt = inst.curve.mul(&mut ops, inst.generator, k);
            let a = d3.decompose(&ctx, &fb, &mut ops, &mut c[0], pt);
            let b = s3.decompose(&ctx, &fb, &mut ops, &mut c[1], pt);
            let m = mitm.decompose(&ctx, &fb, &mut ops, &mut c[2], pt);
            for idx in [&a, &b].into_iter().flatten() {
                let sum = idx.iter().fold(inst.curve.identity(), |acc, &i| {
                    inst.curve.add(&mut ops, acc, fb.points[i])
                });
                assert_eq!(sum, pt);
            }
            assert_eq!(a.is_some(), m.is_some(), "D₃ against the pair table at k = {k}");
            assert_eq!(b.is_some(), m.is_some(), "S₃ against the pair table at k = {k}");
            hits += m.is_some() as u32;
            found += a.is_some() as u32;
        }
        assert!(hits > 0 && found == hits);
        assert!(d3.stats.fp_muls < s3.stats.fp_muls);
    }
}
