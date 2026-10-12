//! # Subfield curves on `E(F_{p³})` and the Frobenius fold (E5).
//!
//! A curve `E / F_p` viewed over `F_{p³}` carries the `p`-power
//! Frobenius `π` as an endomorphism of order 3 on the points, and on the
//! prime-order subgroup outside `E(F_p)` it acts as a primitive cube
//! root of unity `λ_π` modulo `r` (`π³ = 1`, `π ≠ 1` there).  This is
//! the odd-characteristic twin of the Koblitz fold: a `π`-stable base
//! folds by `⟨−1, π⟩` of order 6.
//!
//! The `π`-stable `F_p`-subspaces of `F_{p³} = F_p(s)`, `s³ = ν`, are
//! spanned by the eigenlines `F_p`, `s·F_p` (`s^p = ωs`) and `s²·F_p`.
//! `F_p` is useless — its points are `E(F_p)`, which `π` fixes and which
//! meets the target subgroup only at `O` — so the base here is the line
//! `x ∈ s·F_p^*`, on which `π` acts as `t ↦ ωt`.
//!
//! On a `j = 0` curve (`a = 0`) the automorphism `ζ(x, y) = (ζx, y)`
//! with `ζ ∈ F_p` (here `p ≡ 1 (mod 3)` always) also preserves the
//! line.  Its eigenvalue is a primitive cube root of unity modulo `r`
//! too, so `λ_ζ ∈ {λ_π, λ_π²}`: on `⟨G⟩`, `ζ` **is** `π` or `π²`, and
//! `⟨−1, π, ζ⟩` folds no further than `⟨−1, π⟩`.  That is the general
//! rule the plan's E3 and E5 predictions missed and this module
//! measures: the fold is the number of distinct roots of unity the
//! generators realise modulo `r`, not the product of their orders.
//!
//! Companion to `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md` §E5.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use crate::cryptanalysis::ext_curve::{
    jacobi, random_point, ExtCurve, ExtDiagonalAutomorphism, ExtField, ExtPoint, Fp3, Fp3El,
    FrobeniusType,
};
use crate::cryptanalysis::glv_invariant_base::{
    addm, factor_u64, fold_by_endomorphisms, mulm, subm, verify_endomorphism, Closure,
    Endomorphism, FoldReport, Negation,
};
use crate::cryptanalysis::ic_boundary::{CountedGroup, FactorBase, GroupOps};
use crate::cryptanalysis::residual_walk::{inv_mod, is_prime_u64, sqrt_mod};

pub type Fp3Curve = ExtCurve<Fp3>;
pub type Fp3Point = ExtPoint<Fp3El>;

/// A subfield instance: `E / F_p` on `E(F_{p³})`, its certified order,
/// the prime-order subgroup outside `E(F_p)`, and the verified
/// Frobenius (and `ζ` on `j = 0`).
#[derive(Clone, Debug, Serialize)]
pub struct SubfieldInstance {
    pub name: String,
    pub p: u64,
    pub base_a: u64,
    pub base_b: u64,
    /// `#E(F_p) = p + 1 − t`.
    pub base_order: u64,
    pub trace: i64,
    pub curve: Fp3Curve,
    /// `#E(F_{p³}) = p³ + 1 − (t³ − 3pt)`.
    pub group_order: u64,
    pub r: u64,
    pub cofactor: u64,
    pub generator: Fp3Point,
    pub frobenius: FrobeniusType<Fp3El>,
    /// `ζ` when the base curve has `j = 0`.
    pub zeta: Option<ExtDiagonalAutomorphism<Fp3El>>,
    /// The line's direction, `s`.
    pub line: Fp3El,
}

/// Which generators the line base is folded by.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum SubfieldFold {
    Negation,
    Frobenius,
    Zeta,
    FrobeniusZeta,
}

impl SubfieldFold {
    pub fn name(self) -> &'static str {
        match self {
            Self::Negation => "negation",
            Self::Frobenius => "frobenius",
            Self::Zeta => "zeta",
            Self::FrobeniusZeta => "frobenius+zeta",
        }
    }
}

/// **A subfield instance with `p` of about `p_bits` bits**, `j = 0`
/// when `j0`, deterministic in `seed`.  `#E(F_p)` is counted in
/// `O(p)`, `#E(F_{p³})` follows from the trace and is certified as
/// `[N]P = O` on random points, `r` is the largest prime factor of
/// `#E(F_{p³}) / #E(F_p)` with `r ∤ #E(F_p)` and cofactor at most
/// `max_cofactor_ratio · #E(F_p)`, and the endomorphisms are verified
/// before they are returned.  `p_bits ≤ 20`.
pub fn generate_subfield_instance(
    p_bits: u32,
    seed: u64,
    j0: bool,
    max_cofactor_ratio: u64,
) -> Result<SubfieldInstance, String> {
    if !(5..=20).contains(&p_bits) {
        return Err(format!("p_bits = {p_bits} outside 5..=20"));
    }
    let mut rng = StdRng::seed_from_u64(seed ^ 0x5355_4246_4945_4C44 ^ (p_bits as u64));
    for _ in 0..20_000u32 {
        let p = loop {
            let c = rng.gen_range((1u64 << (p_bits - 1))..(1u64 << p_bits)) | 1;
            if is_prime_u64(c) && c % 3 == 1 {
                break c;
            }
        };
        let f = Fp3::new(p)?;
        let (a, b) = if j0 {
            (0, rng.gen_range(1..p))
        } else {
            (rng.gen_range(1..p), rng.gen_range(1..p))
        };
        let disc = addm(
            mulm(4, mulm(mulm(a, a, p), a, p), p),
            mulm(27, mulm(b, b, p), p),
            p,
        );
        if disc == 0 {
            continue;
        }
        let mut count = 1u64;
        for x in 0..p {
            let rhs = addm(addm(mulm(mulm(x, x, p), x, p), mulm(a, x, p), p), b, p);
            count += (1 + jacobi(rhs, p)) as u64;
        }
        let t = p as i128 + 1 - count as i128;
        let t3 = t * t * t - 3 * (p as i128) * t;
        let n3 = ((p as i128).pow(3) + 1 - t3) as u64;
        if !n3.is_multiple_of(count) {
            return Err(format!("#E(F_p) = {count} does not divide #E(F_p³) = {n3}"));
        }
        let quotient = n3 / count;
        let Some(&(r, _)) = factor_u64(quotient).last() else {
            continue;
        };
        if r < 64 || count.is_multiple_of(r) || r % 3 != 1 {
            continue;
        }
        let h = n3 / r;
        if h > max_cofactor_ratio * count {
            continue;
        }
        let curve = ExtCurve {
            f,
            a: f.from_fp(a),
            b: f.from_fp(b),
        };
        let mut ops = GroupOps::default();
        let pts: Vec<Fp3Point> = (0..3).map(|_| random_point(&curve, &mut rng)).collect();
        if !pts.iter().all(|&pt| curve.mul(&mut ops, pt, n3).infinity) {
            return Err(format!("#E(F_p³) = {n3} is wrong at p = {p}: [N]P ≠ O"));
        }
        let generator = loop {
            let g = curve.mul(&mut ops, random_point(&curve, &mut rng), h);
            if !g.infinity {
                break g;
            }
        };
        // λ_π: a primitive cube root of unity modulo r, decided on G.
        let sq = sqrt_mod(r - 3, r).ok_or("−3 is not a square mod r")?;
        let half = inv_mod(2, r);
        let cands = [
            mulm(subm(sq, 1, r), half, r),
            mulm(subm(r - sq, 1, r), half, r),
        ];
        let frob_g = ExtPoint::affine(f.frob(generator.x), f.frob(generator.y));
        let lambda = cands
            .into_iter()
            .find(|&l| curve.mul(&mut ops, generator, l) == frob_g)
            .ok_or("neither cube root of unity is the eigenvalue of π on G")?;
        let frobenius = FrobeniusType {
            cx: f.one(),
            cy: f.one(),
            eigenvalue: lambda,
            p,
            label: "frobenius",
        };
        verify_endomorphism(&curve, generator, r, &frobenius, 20, seed)?;
        // The two π-stable lines s·F_p and s²·F_p carry the points of
        // the two cubic twists of E over F_p (`x = s t` gives
        // `y² = ν t³ + a s t + b`); on j = 0 each line lies in one of the
        // two eigen-subgroups N(π − ω), N(π − ω²), and only one of them
        // meets ⟨G⟩.  Take the line whose points survive [h].
        let line = {
            let s1 = f.s();
            let s2 = f.mul(s1, s1);
            let survivors = |line: Fp3El| -> usize {
                let mut n = 0usize;
                let mut o = GroupOps::default();
                for t in 1..p.min(64) {
                    for pt in curve.lift_x(f.scale(line, t)) {
                        if !curve.mul(&mut o, pt, h).infinity {
                            n += 1;
                        }
                    }
                }
                n
            };
            if survivors(s2) > survivors(s1) {
                s2
            } else {
                s1
            }
        };
        let zeta = if j0 {
            let zeta_fp = f.omega;
            let zeta_g = ExtPoint::affine(f.scale(generator.x, zeta_fp), generator.y);
            let eig = cands
                .into_iter()
                .find(|&l| curve.mul(&mut ops, generator, l) == zeta_g)
                .ok_or("neither cube root of unity is the eigenvalue of ζ on G")?;
            let z = ExtDiagonalAutomorphism {
                cx: f.from_fp(zeta_fp),
                cy: f.one(),
                eigenvalue: eig,
                order: 3,
                label: "zeta3",
            };
            verify_endomorphism(&curve, generator, r, &z, 20, seed)?;
            Some(z)
        } else {
            None
        };
        return Ok(SubfieldInstance {
            name: format!("subfield{}-p{p_bits}bit-p{p}", if j0 { "-j0" } else { "" }),
            p,
            base_a: a,
            base_b: b,
            base_order: count,
            trace: t as i64,
            curve,
            group_order: n3,
            r,
            cofactor: h,
            generator,
            frobenius,
            zeta,
            line,
        });
    }
    Err(format!(
        "no subfield instance at p_bits = {p_bits} with cofactor ratio ≤ {max_cofactor_ratio}"
    ))
}

/// **The `π`-stable line base** `{P : x(P) ∈ s·F_p^*}` folded by the
/// chosen generators (always with negation).  The set is checked to be
/// invariant, not assumed.
pub fn subfield_line_base(
    inst: &SubfieldInstance,
    fold: SubfieldFold,
) -> Result<(FactorBase<Fp3Point>, FoldReport), String> {
    let start = std::time::Instant::now();
    let curve = &inst.curve;
    let f = &curve.f;
    let mut seed = Vec::with_capacity(inst.p as usize);
    for t in 1..inst.p {
        seed.extend(curve.lift_x(f.scale(inst.line, t)));
    }
    let neg = Negation { r: inst.r };
    let mut gens: Vec<&dyn Endomorphism<Fp3Curve>> = vec![&neg];
    match fold {
        SubfieldFold::Negation => {}
        SubfieldFold::Frobenius => gens.push(&inst.frobenius),
        SubfieldFold::Zeta => gens.push(inst.zeta.as_ref().ok_or("no ζ: not a j = 0 instance")?),
        SubfieldFold::FrobeniusZeta => {
            gens.push(&inst.frobenius);
            gens.push(inst.zeta.as_ref().ok_or("no ζ: not a j = 0 instance")?);
        }
    }
    let description = format!(
        "the π-stable line x ∈ s·F_p ({} abscissae), folded by {}",
        inst.p - 1,
        fold.name()
    );
    let (mut fb, report) = fold_by_endomorphisms(
        curve,
        inst.r,
        inst.cofactor,
        seed,
        &gens,
        Closure::Strict,
        16,
        |p| curve.key(p),
        |p| f.pack(p.x),
        description,
    )?;
    fb.cost.wall_ns = start.elapsed().as_nanos() as u64;
    fb.cost.count("abscissae_scanned", inst.p - 1);
    fb.cost.count("sqrt_solves", inst.p - 1);
    Ok((fb, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::glv_invariant_base::eigenvalue_order;

    #[test]
    fn a_subfield_instance_has_a_frobenius_of_order_three_on_the_subgroup() {
        let inst = generate_subfield_instance(8, 1, false, 8).unwrap();
        assert_eq!(eigenvalue_order(inst.frobenius.eigenvalue, inst.r), 3);
        assert!(inst.group_order.is_multiple_of(inst.base_order));
        assert!(!inst.base_order.is_multiple_of(inst.r));
        // π³ = 1 on every point of E(F_{p³}).
        let mut rng = StdRng::seed_from_u64(2);
        for _ in 0..10 {
            let pt = random_point(&inst.curve, &mut rng);
            let mut q = pt;
            for _ in 0..3 {
                q = inst.frobenius.apply(&inst.curve, q);
            }
            assert_eq!(q, pt);
        }
    }

    #[test]
    fn the_frobenius_line_folds_six_to_one() {
        let inst = generate_subfield_instance(7, 3, false, 8).unwrap();
        let (folded, rep) = subfield_line_base(&inst, SubfieldFold::Frobenius).unwrap();
        let (control, crep) = subfield_line_base(&inst, SubfieldFold::Negation).unwrap();
        assert_eq!(folded.points.len(), control.points.len());
        assert!(
            (rep.points_per_orbit - 6.0).abs() < 1e-9,
            "{}",
            rep.points_per_orbit
        );
        assert!((crep.points_per_orbit - 2.0).abs() < 1e-9);
        assert_eq!(folded.columns * 3, control.columns);
    }

    /// **`ζ` is `π` or `π²` on the subgroup**: adding it to the
    /// Frobenius fold changes nothing, and its eigenvalue is a power of
    /// the Frobenius's.
    #[test]
    fn zeta_coincides_with_a_power_of_frobenius_on_the_subgroup() {
        // On j = 0 the new part N(π² + π + 1) splits as N(π − ω)·N(π − ω²)
        // in Z[ω], two pieces of about p each, so r ≈ p and h ≈ p·#E(F_p).
        let inst = generate_subfield_instance(7, 5, true, 1024).unwrap();
        let zeta = inst.zeta.as_ref().unwrap();
        let lp = inst.frobenius.eigenvalue;
        let lz = zeta.eigenvalue;
        assert!(
            lz == lp || lz == mulm(lp, lp, inst.r),
            "λ_ζ ∈ {{λ_π, λ_π²}}"
        );
        let (both, rep_both) = subfield_line_base(&inst, SubfieldFold::FrobeniusZeta).unwrap();
        let (frob, rep_frob) = subfield_line_base(&inst, SubfieldFold::Frobenius).unwrap();
        let (zeta_only, rep_zeta) = subfield_line_base(&inst, SubfieldFold::Zeta).unwrap();
        assert_eq!(
            both.columns, frob.columns,
            "⟨π, ζ⟩ folds no further than ⟨π⟩"
        );
        assert_eq!(zeta_only.columns, frob.columns);
        assert!(
            (rep_both.points_per_orbit - 6.0).abs() < 1e-9,
            "{}",
            rep_both.points_per_orbit
        );
        assert!((rep_frob.points_per_orbit - 6.0).abs() < 1e-9);
        assert!((rep_zeta.points_per_orbit - 6.0).abs() < 1e-9);
    }
}
