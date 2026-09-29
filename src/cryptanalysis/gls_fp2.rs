//! # GLS curves over `F_{p²}` and their `ψ`-invariant line base (type B).
//!
//! Galbraith–Lin–Scott take a curve `E / F_p` and its quadratic twist
//! `E' / F_{p²}`, `E' : y² = x³ + a u² x + b u³` for a non-square
//! `u ∈ F_{p²}`, and conjugate the `p`-power Frobenius through the twist
//! isomorphism `τ_u : (x, y) ↦ (u x, u^{3/2} y)`:
//!
//! ```text
//!     ψ = τ_u ∘ π_p ∘ τ_u⁻¹ :  (x, y) ↦ (u^{1−p} x^p,  u^{3(1−p)/2} y^p),
//!     ψ² = −1 on E'(F_{p²}),   #E'(F_{p²}) = (p − 1)² + t².
//! ```
//!
//! `ψ` is a Frobenius-type endomorphism of degree `p` whose eigenvalue
//! has order **4** modulo `r` (`λ² ≡ −1`), so it folds a factor base
//! exactly as the Koblitz Frobenius does, by `⟨−1, ψ⟩` of order 4.  The
//! base that is invariant by construction is a **line**: writing
//! `F_{p²} = F_p(s)`, `s² = ν`, the abscissa map `x ↦ u^{1−p} x^p` sends
//! `u·s·t ↦ −u·s·t` for `t ∈ F_p` (because `s^p = −s`), so
//!
//! ```text
//!     F = { P ∈ E'(F_{p²}) : x(P) ∈ u·s·F_p }
//! ```
//!
//! is `ψ`-stable with `4` signed points per column.  (The other
//! `ψ`-stable line, `u·F_p`, carries no points: `x = u t` gives
//! `y² = u³(t³ + at + b)`, a non-square times a square.)  Both facts are
//! checked at run time, not assumed: [`gls_line_base`] folds under
//! [`Closure::Strict`], which refuses a base that is not invariant, and
//! [`generate_gls_instance`] verifies `ψ` on random points.
//!
//! This is the same construction as the Koblitz fold in a different
//! field: there the `π`-stable sets are `F_2`-subspaces of `F_{2^n}`
//! that are `π`-stable, here the `ψ`-stable sets are the `F_p`-lines of
//! `F_{p²}` that `x ↦ u^{1−p}x^p` preserves.  A `4`-dimensional GLV+GLS
//! curve (`j = 1728` twisted over `F_{p²}`, FourQ-like) would fold by
//! `⟨ι, ψ⟩` of order 8 through the same [`fold_by_endomorphisms`] call
//! with both generators; that is the next family the plan lists.
//!
//! The group is a [`CountedGroup`], so the framework's oracles, relation
//! loop and reports run on it unchanged.  `p < 2^31` so that a packed
//! `F_{p²}` element fits the framework's `u64` keys.
//!
//! Companion to `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use crate::cryptanalysis::ext_curve::random_point;
use crate::cryptanalysis::glv_invariant_base::{
    addm, factor_u64, fold_by_endomorphisms, mulm, verify_endomorphism, Closure, Endomorphism,
    FoldReport, Negation,
};
use crate::cryptanalysis::ic_boundary::{CountedGroup, FactorBase, GroupOps};
use crate::cryptanalysis::residual_walk::{inv_mod, is_prime_u64, pow_mod, sqrt_mod};

use crate::cryptanalysis::ext_curve::ExtDiagonalAutomorphism;
pub use crate::cryptanalysis::ext_curve::{jacobi, ExtCurve, ExtField, ExtPoint, Fp2, Fp2El};

/// The GLS twist as a counted group: `ExtCurve` over `F_{p²}`.
pub type Fp2Curve = ExtCurve<Fp2>;
/// A point of the twist.
pub type Fp2Point = ExtPoint<Fp2El>;

/// `ψ(x, y) = (cx · x^p, cy · y^p)`.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct GlsEndomorphism {
    pub cx: Fp2El,
    pub cy: Fp2El,
    pub eigenvalue: u64,
    pub p: u64,
}

impl Endomorphism<Fp2Curve> for GlsEndomorphism {
    fn name(&self) -> String {
        "gls-psi".into()
    }
    fn degree(&self) -> u64 {
        self.p
    }
    fn eigenvalue(&self) -> u64 {
        self.eigenvalue
    }
    fn apply(&self, g: &Fp2Curve, pt: Fp2Point) -> Fp2Point {
        if pt.infinity {
            return pt;
        }
        let f = &g.f;
        Fp2Point::affine(f.mul(self.cx, f.frob(pt.x)), f.mul(self.cy, f.frob(pt.y)))
    }
}

/// A GLS instance: the twist, its certified order, the prime-order
/// subgroup, and the verified `ψ`.
#[derive(Clone, Debug, Serialize)]
pub struct GlsInstance {
    pub name: String,
    pub p: u64,
    /// The base curve `y² = x³ + ax + b` over `F_p` and its trace.
    pub base_a: u64,
    pub base_b: u64,
    pub trace: i64,
    pub u: Fp2El,
    pub curve: Fp2Curve,
    /// `(p − 1)² + t²`.
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    pub generator: Fp2Point,
    pub psi: GlsEndomorphism,
    /// The `ψ`-stable line's direction, `u·s`.
    pub line: Fp2El,
    /// `generic`, `j0` or `j1728`: the base curve's family.
    pub family: &'static str,
    /// The base curve's automorphism lifted to the twist (`ζ` on
    /// `j = 0`, `ι` on `j = 1728`), with its eigenvalue decided on `G`.
    pub automorphism: Option<ExtDiagonalAutomorphism<Fp2El>>,
}

/// The base curve's family for [`generate_gls_instance_of`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum GlsFamily {
    Generic,
    /// `y² = x³ + b` over `F_p`, `p ≡ 1 (mod 3)`, twisted: carries `ζ`.
    J0,
    /// `y² = x³ + ax` over `F_p`, `p ≡ 1 (mod 4)`, twisted: carries `ι`.
    J1728,
}

impl GlsFamily {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "generic" => Ok(Self::Generic),
            "j0" => Ok(Self::J0),
            "j1728" => Ok(Self::J1728),
            other => Err(format!(
                "unknown GLS family `{other}`; try generic, j0 or j1728"
            )),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::Generic => "generic",
            Self::J0 => "j0",
            Self::J1728 => "j1728",
        }
    }
}

/// Which generators [`gls_line_base_by`] folds the line by (always
/// with negation).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum GlsFold {
    Negation,
    Psi,
    Aut,
    PsiAut,
}

impl GlsFold {
    pub fn name(self) -> &'static str {
        match self {
            Self::Negation => "negation",
            Self::Psi => "psi",
            Self::Aut => "aut",
            Self::PsiAut => "psi+aut",
        }
    }
}

/// **A GLS instance with `p` of about `p_bits` bits**, deterministic in
/// `seed`.  The base curve's order is counted in `O(p)`, the twist's
/// order `(p − 1)² + t²` is certified as `[N]P = O` on random points,
/// `r` is its largest prime factor with cofactor at most
/// `max_cofactor`, and `ψ` is verified on random points before it is
/// returned.  `p_bits ≤ 24`.
pub fn generate_gls_instance(
    p_bits: u32,
    seed: u64,
    max_cofactor: u64,
) -> Result<GlsInstance, String> {
    generate_gls_instance_of(GlsFamily::Generic, p_bits, seed, max_cofactor)
}

/// [`generate_gls_instance`] with the base curve drawn from `family`:
/// `j = 0` (`p ≡ 1 (mod 3)`) or `j = 1728` (`p ≡ 1 (mod 4)`) base
/// curves carry their automorphism, lifted to the twist and verified,
/// so that the composite fold `⟨−1, ψ, aut⟩` can be measured (E3).
pub fn generate_gls_instance_of(
    family: GlsFamily,
    p_bits: u32,
    seed: u64,
    max_cofactor: u64,
) -> Result<GlsInstance, String> {
    if !(6..=24).contains(&p_bits) {
        return Err(format!("p_bits = {p_bits} outside 6..=24"));
    }
    let mut rng = StdRng::seed_from_u64(seed ^ 0x474C_535F_5053_4900 ^ (p_bits as u64));
    for _ in 0..20_000u32 {
        let p = loop {
            let c = rng.gen_range((1u64 << (p_bits - 1))..(1u64 << p_bits)) | 1;
            let admissible = match family {
                GlsFamily::Generic => true,
                GlsFamily::J0 => c % 3 == 1,
                GlsFamily::J1728 => c % 4 == 1,
            };
            if admissible && is_prime_u64(c) {
                break c;
            }
        };
        let f = Fp2::new(p)?;
        let (a, b) = match family {
            GlsFamily::Generic => (rng.gen_range(1..p), rng.gen_range(1..p)),
            GlsFamily::J0 => (0, rng.gen_range(1..p)),
            GlsFamily::J1728 => (rng.gen_range(1..p), 0),
        };
        let disc = addm(
            mulm(4, mulm(mulm(a, a, p), a, p), p),
            mulm(27, mulm(b, b, p), p),
            p,
        );
        if disc == 0 {
            continue;
        }
        // #E(F_p) by the Legendre symbol at every abscissa.
        let mut count = 1u64;
        for x in 0..p {
            let rhs = addm(addm(mulm(mulm(x, x, p), x, p), mulm(a, x, p), p), b, p);
            count += (1 + jacobi(rhs, p)) as u64;
        }
        let t = p as i64 + 1 - count as i64;
        let order = ((p - 1) as u128 * (p - 1) as u128 + (t as i128 * t as i128) as u128) as u64;
        let factors = factor_u64(order);
        let Some(&(r, _)) = factors.last() else {
            continue;
        };
        let h = order / r;
        if h > max_cofactor || r < 64 || r % 4 != 1 {
            continue;
        }
        // A non-square u ∈ F_{p²} and the twist.
        let u = loop {
            let c = [rng.gen_range(1..p), rng.gen_range(1..p)];
            if !f.is_square(c) {
                break c;
            }
        };
        let u2 = f.sqr(u);
        let u3 = f.mul(u2, u);
        let curve = Fp2Curve {
            f,
            a: f.mul(f.from_fp(a), u2),
            b: f.mul(f.from_fp(b), u3),
        };
        let mut ops = GroupOps::default();
        let pts: Vec<Fp2Point> = (0..3).map(|_| random_point(&curve, &mut rng)).collect();
        if !pts
            .iter()
            .all(|&pt| curve.mul(&mut ops, pt, order).infinity)
        {
            return Err(format!(
                "the twist's order (p − 1)² + t² = {order} is wrong at p = {p}: [N]P ≠ O"
            ));
        }
        let generator = loop {
            let g = curve.mul(&mut ops, random_point(&curve, &mut rng), h);
            if !g.infinity {
                break g;
            }
        };
        // ψ: cx = u^{1−p} = u / u^p, cy = √(cx³).
        let cx = f.mul(u, f.inv(f.frob(u)).expect("u ≠ 0"));
        let Some(cy0) = f.sqrt(f.mul(f.sqr(cx), cx)) else {
            return Err("cx³ is not a square in F_{p²}; the GLS derivation is wrong".into());
        };
        let s = sqrt_mod(r - 1, r).ok_or("r ≢ 1 (mod 4)")?;
        let mut psi = None;
        'outer: for cy in [cy0, f.neg(cy0)] {
            for lambda in [s, r - s] {
                let cand = GlsEndomorphism {
                    cx,
                    cy,
                    eigenvalue: lambda,
                    p,
                };
                if cand.apply(&curve, generator) == curve.mul(&mut ops, generator, lambda) {
                    psi = Some(cand);
                    break 'outer;
                }
            }
        }
        let Some(psi) = psi else {
            return Err(format!("no sign of ψ acts as √−1 on G at p = {p}"));
        };
        verify_endomorphism(&curve, generator, r, &psi, 20, seed)?;
        // The base curve's automorphism, lifted: it commutes with the
        // twist isomorphism since its constants lie in F_p.
        let automorphism = match family {
            GlsFamily::Generic => None,
            GlsFamily::J0 => {
                if r % 3 != 1 {
                    continue;
                }
                let zeta = (2..p)
                    .map(|g| pow_mod(g, (p - 1) / 3, p))
                    .find(|&z| z != 1)
                    .ok_or("no cube root of unity")?;
                let sq3 = sqrt_mod(r - 3, r).ok_or("−3 is not a square mod r")?;
                let half = inv_mod(2, r);
                let cands = [
                    mulm(addm(r - 1, sq3, r), half, r),
                    mulm(addm(r - 1, r - sq3, r), half, r),
                ];
                let img = ExtPoint::affine(f.scale(generator.x, zeta), generator.y);
                let eig = cands
                    .into_iter()
                    .find(|&l| curve.mul(&mut ops, generator, l) == img)
                    .ok_or("neither root of λ² + λ + 1 is the eigenvalue of ζ on G")?;
                Some(ExtDiagonalAutomorphism {
                    cx: f.from_fp(zeta),
                    cy: f.one(),
                    eigenvalue: eig,
                    order: 3,
                    label: "zeta3",
                })
            }
            GlsFamily::J1728 => {
                let i = (2..p)
                    .map(|g| pow_mod(g, (p - 1) / 4, p))
                    .find(|&z| mulm(z, z, p) == p - 1)
                    .ok_or("no fourth root of unity")?;
                let img = ExtPoint::affine(f.neg(generator.x), f.scale(generator.y, i));
                let eig = [s, r - s]
                    .into_iter()
                    .find(|&l| curve.mul(&mut ops, generator, l) == img)
                    .ok_or("neither square root of −1 is the eigenvalue of ι on G")?;
                Some(ExtDiagonalAutomorphism {
                    cx: f.from_fp(p - 1),
                    cy: f.from_fp(i),
                    eigenvalue: eig,
                    order: 4,
                    label: "iota4",
                })
            }
        };
        if let Some(aut) = &automorphism {
            verify_endomorphism(&curve, generator, r, aut, 20, seed)?;
        }
        return Ok(GlsInstance {
            name: format!("gls-{}-p{p_bits}bit-p{p}", family.name()),
            p,
            base_a: a,
            base_b: b,
            trace: t,
            u,
            curve,
            group_order: order,
            r,
            cofactor: h,
            generator,
            psi,
            line: f.mul(u, f.s()),
            family: family.name(),
            automorphism,
        });
    }
    Err(format!(
        "no GLS instance found at p_bits = {p_bits} with cofactor ≤ {max_cofactor}"
    ))
}

/// **The `ψ`-stable line base** `{P : x(P) ∈ u·s·F_p^*}`, folded by
/// `⟨−1, ψ⟩` (`fold = true`, four signed points per column) or by
/// negation alone (the control on the same points).  The set is
/// checked to be invariant, not assumed.
pub fn gls_line_base(
    inst: &GlsInstance,
    fold: bool,
) -> Result<(FactorBase<Fp2Point>, FoldReport), String> {
    gls_line_base_by(
        inst,
        if fold {
            GlsFold::Psi
        } else {
            GlsFold::Negation
        },
    )
}

/// The line base folded by the chosen generators.  Under `Aut` and
/// `PsiAut` the seed is closed under the maps (`Closure::Close`): `ζ`
/// with `ζ ∈ F_p` and `ι` preserve the line, so nothing is added and
/// the closure is a check.
pub fn gls_line_base_by(
    inst: &GlsInstance,
    fold: GlsFold,
) -> Result<(FactorBase<Fp2Point>, FoldReport), String> {
    let start = std::time::Instant::now();
    let curve = &inst.curve;
    let f = &curve.f;
    let mut seed = Vec::with_capacity(inst.p as usize);
    let mut sqrts = 0u64;
    for t in 1..inst.p {
        let x = f.scale(inst.line, t);
        sqrts += 1;
        seed.extend(curve.lift_x(x));
    }
    let neg = Negation { r: inst.r };
    let mut gens: Vec<&dyn Endomorphism<Fp2Curve>> = vec![&neg];
    match fold {
        GlsFold::Negation => {}
        GlsFold::Psi => gens.push(&inst.psi),
        GlsFold::Aut => gens.push(
            inst.automorphism
                .as_ref()
                .ok_or("no automorphism on this family")?,
        ),
        GlsFold::PsiAut => {
            gens.push(&inst.psi);
            gens.push(
                inst.automorphism
                    .as_ref()
                    .ok_or("no automorphism on this family")?,
            );
        }
    }
    let description = format!(
        "the ψ-stable line x ∈ u·s·F_p ({} abscissae), folded by {}",
        inst.p - 1,
        fold.name()
    );
    let (mut fb, report) = fold_by_endomorphisms(
        curve,
        inst.r,
        inst.cofactor,
        seed,
        &gens,
        Closure::Close,
        16,
        |p| curve.key(p),
        |p| f.pack(p.x),
        description,
    )?;
    fb.cost.wall_ns = start.elapsed().as_nanos() as u64;
    fb.cost.count("abscissae_scanned", inst.p - 1);
    fb.cost.count("sqrt_solves", sqrts);
    Ok((fb, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::glv_invariant_base::eigenvalue_order;
    use crate::cryptanalysis::ic_boundary::{collect_and_solve, RestartPool, TargetSource};
    use crate::cryptanalysis::ic_framework::plugins::SubtractOracle;
    use crate::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx};

    #[test]
    fn fp2_arithmetic_is_a_field() {
        let f = Fp2::new(1009).unwrap();
        let mut rng = StdRng::seed_from_u64(3);
        for _ in 0..200 {
            let a = [rng.gen_range(0..1009), rng.gen_range(0..1009)];
            let b = [rng.gen_range(0..1009), rng.gen_range(0..1009)];
            assert_eq!(f.mul(a, b), f.mul(b, a));
            if !f.is_zero(a) {
                assert_eq!(f.mul(a, f.inv(a).unwrap()), f.one());
            }
            // Frobenius is a field automorphism of order 2.
            assert_eq!(f.frob(f.mul(a, b)), f.mul(f.frob(a), f.frob(b)));
            assert_eq!(f.frob(f.frob(a)), a);
            assert_eq!(f.pow(a, 1009), f.frob(a), "a^p is the conjugate");
            // Every square has a root, and non-squares have none.
            let sq = f.sqr(a);
            let root = f.sqrt(sq).expect("a square has a root");
            assert!(root == a || root == f.neg(a));
            assert_eq!(f.is_square(a), f.sqrt(a).is_some());
        }
        assert_eq!(f.pow(f.s(), 1009), f.neg(f.s()), "s^p = −s");
    }

    #[test]
    fn a_gls_instance_has_a_verified_psi_of_order_four() {
        let inst = generate_gls_instance(10, 1, 16).unwrap();
        assert_eq!(
            inst.group_order as i128,
            (inst.p as i128 - 1).pow(2) + (inst.trace as i128).pow(2)
        );
        let check =
            verify_endomorphism(&inst.curve, inst.generator, inst.r, &inst.psi, 30, 1).unwrap();
        assert_eq!(check.eigenvalue_order, 4);
        assert_eq!(eigenvalue_order(inst.psi.eigenvalue, inst.r), 4);
        // ψ² = −1 on random points of the whole group, not only ⟨G⟩.
        let mut rng = StdRng::seed_from_u64(9);
        for _ in 0..20 {
            let pt = random_point(&inst.curve, &mut rng);
            let psi2 = inst.psi.apply(&inst.curve, inst.psi.apply(&inst.curve, pt));
            assert_eq!(psi2, inst.curve.neg(pt), "ψ² = −1");
        }
    }

    /// The line `u·s·F_p` is `ψ`-stable and `u·F_p` carries no points:
    /// the two facts the module header derives.
    #[test]
    fn the_line_is_psi_stable_and_the_other_line_is_empty() {
        let inst = generate_gls_instance(9, 2, 16).unwrap();
        let f = &inst.curve.f;
        for t in 1..inst.p {
            let x = f.scale(inst.line, t);
            let image = f.mul(inst.psi.cx, f.frob(x));
            assert_eq!(image, f.neg(x), "ψ_x(u s t) = −u s t");
        }
        let mut on_other_line = 0;
        for t in 1..inst.p {
            on_other_line += inst.curve.lift_x(f.scale(inst.u, t)).len();
        }
        assert!(
            on_other_line <= 3,
            "u·F_p carries only 2-torsion: {on_other_line}"
        );
    }

    /// **E3, the composite fold.**  On a `j = 0` twist the eigenvalues of
    /// `ζ` (order 3) and `ψ` (order 4) generate the twelfth roots of
    /// unity, so `⟨−1, ψ, ζ⟩` folds twelve points to a column; on a
    /// `j = 1728` twist `ψ` and `ι` both square to `−1` and coincide on
    /// the subgroup up to sign, so `⟨−1, ψ, ι⟩` folds four, as `⟨−1, ψ⟩`
    /// alone does.
    #[test]
    fn composite_folds_are_the_distinct_roots_of_unity_the_generators_realise() {
        let inst = generate_gls_instance_of(GlsFamily::J0, 8, 7, 64).unwrap();
        let aut = inst.automorphism.as_ref().unwrap();
        assert_eq!(eigenvalue_order(aut.eigenvalue, inst.r), 3);
        let (_, psi) = gls_line_base_by(&inst, GlsFold::Psi).unwrap();
        let (_, both) = gls_line_base_by(&inst, GlsFold::PsiAut).unwrap();
        assert!((psi.points_per_orbit - 4.0).abs() < 1e-9);
        assert!(
            (both.points_per_orbit - 12.0).abs() < 1e-9,
            "{}",
            both.points_per_orbit
        );

        // A j = 1728 twist's order always carries a cofactor of about p
        // (its `ψ` and `ι` coincide up to sign as maps of the subgroup,
        // and the twist is isogenous to a curve over F_p), so r ≈ p.
        let inst = generate_gls_instance_of(GlsFamily::J1728, 8, 8, 4096).unwrap();
        let aut = inst.automorphism.as_ref().unwrap();
        assert_eq!(eigenvalue_order(aut.eigenvalue, inst.r), 4);
        assert!(
            aut.eigenvalue == inst.psi.eigenvalue || aut.eigenvalue == inst.r - inst.psi.eigenvalue,
            "ι = ±ψ on the subgroup"
        );
        let (_, psi) = gls_line_base_by(&inst, GlsFold::Psi).unwrap();
        let (_, both) = gls_line_base_by(&inst, GlsFold::PsiAut).unwrap();
        assert!((psi.points_per_orbit - 4.0).abs() < 1e-9);
        assert!(
            (both.points_per_orbit - 4.0).abs() < 1e-9,
            "{}",
            both.points_per_orbit
        );
    }

    #[test]
    fn the_line_base_folds_four_to_one_and_recovers_the_logarithm() {
        let inst = generate_gls_instance(9, 3, 16).unwrap();
        let (folded, rep) = gls_line_base(&inst, true).unwrap();
        let (control, crep) = gls_line_base(&inst, false).unwrap();
        assert_eq!(folded.points.len(), control.points.len());
        assert!(
            (rep.points_per_orbit - 4.0).abs() < 1e-9,
            "{}",
            rep.points_per_orbit
        );
        assert!((crep.points_per_orbit - 2.0).abs() < 1e-9);
        assert_eq!(folded.columns * 2, control.columns);
        assert!(
            folded.points.len() > inst.p as usize / 2,
            "about p signed points"
        );

        let planted = 4321 % inst.r;
        for fb in [&folded, &control] {
            let mut ops = GroupOps::default();
            let target = inst.curve.mul(&mut ops, inst.generator, planted);
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
            let mut oracle = SubtractOracle;
            let out = collect_and_solve(
                &inst.curve,
                inst.generator,
                target,
                inst.r,
                inst.cofactor,
                fb,
                5,
                2_000_000,
                TargetSource::Walk,
                RestartPool::Lazy,
                |ops, ctr, pt| oracle.decompose(&ctx, fb, ops, ctr, pt),
            );
            assert_eq!(out.recovered, Some(planted), "{}", fb.description);
            assert!(out.verified);
        }
    }
}
