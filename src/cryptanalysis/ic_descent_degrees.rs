//! # The degree a Weil-descent system actually reaches.
//!
//! Petit and Quisquater's Table 2 (*On Polynomial Systems Arising from
//! a Weil Descent*, ASIACRYPT 2012, p. 461) reports, for a handful of
//! `(curve family, n, n', m)` cells, the **average maximal degree
//! reached** in a Gröbner-basis computation, the average time and the
//! peak memory.  Its point is not the timings: it is that in every cell
//! the degree reached came out *below* the first-fall-degree bound,
//! because Semaev's polynomials are sparse and the bound is derived for
//! a generic system.  The bound used here is the semi-regular degree
//! of a boolean system with the same equation degrees, which plays the
//! same role: it is what a system with no exploitable structure would
//! reach.
//!
//! This module measures the same three quantities on the descent
//! systems this repository builds, so that the phenomenon can be
//! checked here rather than cited.
//!
//! ## What is and is not being reproduced
//!
//! **Not reproduced**: Petit–Quisquater's *numbers*.  Their Table 2
//! solves a symmetrised system in `mt + 1 = m² + 1` variables — five at
//! `m = 2`, ten at `m = 3` — obtained from the block structure of
//! Section 4 of that paper.  [`weil_descend_s3`] and [`weil_descend_s4`]
//! build the plain descent, `m·n'` boolean variables with no
//! symmetrisation, so at `(n, n', m) = (11, 6, 2)` this module solves a
//! 12-variable system where they solve a 5-variable one.  The degrees
//! are therefore not comparable row by row and are not presented as
//! though they were.
//!
//! **Reproduced**: the *shape* of the measurement and the comparison
//! that gives it meaning — a derived degree bound in one column, the
//! degree actually reached beside it, and their ratio.
//!
//! ## What the columns mean
//!
//! | column | meaning |
//! |:--|:--|
//! | `family` | `K` for the Koblitz curve `y² + xy = x³ + a x² + 1`, `R` for a random binary curve of the same degree |
//! | `n` | the field degree: the system is over `F_{2^n}` and descends to `n` equations |
//! | `n'` | the `F_2`-dimension of the subspace `V ⊂ F_{2^n}` the summands are drawn from |
//! | `m` | summands: `m = 2` descends `S₃`, `m = 3` descends `S₄` |
//! | `vars` | `m · n'`, the boolean unknowns |
//! | `D_av` | mean over targets of the **solving degree**: the highest degree at which the run produced a new basis element |
//! | `D_max` | the largest solving degree over the targets |
//! | `D_pair` | mean highest degree of a pair *processed*, which Buchberger's strategy pushes above the solving degree |
//! | `D_sr` | the semi-regular degree of a system with the same equation degrees — the derived bound `D_av` is measured against |
//! | `ops` | monomial operations, the metric |
//! | `enumerate` | the reference in the same unit: evaluating every equation at every point of the subspace |
//! | `ops/enum` | above one means the Gröbner basis costs more than enumerating |
//! | `ms` | wall time, the practicality note |
//! | `KiB` | the engine's own peak footprint: monomials held × 8 bytes |
//!
//! ## This is a stage diagnostic
//!
//! `AGENTS.md` §2 is explicit that a number pricing one phase is never a
//! speed.  Everything here prices **one decomposition oracle call** on
//! one target.  It says nothing on its own about whether index calculus
//! beats rho; the whole-pipeline unit `S` and its boundaries live in
//! `ic_boundary` and the ledger note.  No row of this table may be
//! quoted as a speedup.

use serde::Serialize;

use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement};
use crate::cryptanalysis::ic_boundary::{koblitz_instance, random_binary_instance, BinaryInstance};
use crate::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use crate::cryptanalysis::pq_descent::{weil_descend_s3, weil_descend_s4};
use crate::cryptanalysis::pq_groebner_f2::{groebner_basis_f2_stats, GbStats};
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// The caps the descent's truth-table construction imposes: it
/// enumerates `2^{m·n'}` inputs to build the algebraic normal form.
pub const MAX_N_PRIME_M2: u32 = 8;
pub const MAX_N_PRIME_M3: u32 = 5;

/// The largest subspace dimension this module will descend at `m`
/// summands.
pub fn max_n_prime(m: u32) -> u32 {
    match m {
        2 => MAX_N_PRIME_M2,
        3 => MAX_N_PRIME_M3,
        _ => 0,
    }
}

/// One Gröbner run on one target.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct TargetRun {
    pub stats: GbStats,
    /// Whether the basis is `{1}`: the system is inconsistent, so this
    /// target has no decomposition over `V`.
    pub inconsistent: bool,
}

/// One `(family, n, n', m)` cell, over every target drawn for it.
#[derive(Clone, Debug, Serialize)]
pub struct DescentCell {
    pub family: String,
    pub curve: String,
    pub n: u32,
    pub n_prime: u32,
    pub summands: u32,
    pub n_vars: usize,
    pub equations: usize,
    pub targets: usize,
    /// Targets whose system was inconsistent — no decomposition in `V`.
    pub inconsistent: usize,
    /// Mean and max over targets of the **solving degree**: the
    /// highest degree at which the run produced a new basis element.
    /// This is the column the bound is for.
    pub d_av: f64,
    pub d_max: u32,
    /// The same for the highest degree of a pair *processed*, which
    /// Buchberger's strategy pushes above the solving degree; it is a
    /// property of the pair selection, not of the system.
    pub d_pair_av: f64,
    pub d_pair_max: u32,
    /// The same for the highest degree of any intermediate polynomial,
    /// which a reduction can push above the pair degree.
    pub d_poly_av: f64,
    pub d_poly_max: u32,
    /// The derived bound: the semi-regular degree of a boolean system
    /// with the same equation degrees, over the targets drawn.
    pub d_semireg_min: Option<u32>,
    pub d_semireg_max: Option<u32>,
    pub d_semireg_unbounded: usize,
    /// `d_av / d_semireg_min`: below one is Petit–Quisquater's phenomenon.
    pub d_av_over_semireg: Option<f64>,
    /// The metric: monomial operations per target.
    pub mono_ops_mean: f64,
    pub mono_ops_max: u64,
    /// **The reference**, in the same unit: what it costs to solve the
    /// same problem by evaluating every equation at every point of the
    /// subspace, `2^vars · Σ_i |terms_i|` monomial tests.  `AGENTS.md`
    /// §1 asks a method to be priced against the best algorithm that
    /// already solves the problem, and at these sizes that is
    /// exhaustive search.
    pub brute_force_ops: f64,
    /// `mono_ops_mean / brute_force_ops`.  Above one means the Gröbner
    /// basis costs more than enumerating the subspace.
    pub ops_over_brute_force: f64,
    pub spolys_mean: f64,
    pub basis_len_mean: f64,
    /// The practicality note, never the metric.
    pub ms_mean: f64,
    pub ms_max: f64,
    /// The engine's own peak footprint, monomials held × 8 bytes.
    pub peak_kib_max: f64,
    pub per_target: Vec<TargetRun>,
}

fn mean(xs: impl Iterator<Item = f64>) -> f64 {
    let v: Vec<f64> = xs.collect();
    if v.is_empty() {
        f64::NAN
    } else {
        v.iter().sum::<f64>() / v.len() as f64
    }
}

/// Build the `BinaryCurve` an instance describes.
fn curve_of(inst: &BinaryInstance) -> Option<BinaryCurve> {
    let irreducible = find_irreducible_sparse(inst.n)?;
    Some(BinaryCurve {
        m: inst.n,
        irreducible,
        a: inst.gf.to_element(inst.a),
        b: inst.gf.to_element(inst.b),
        generator: BinaryPoint::Infinity,
        order: BigUint::from(inst.r),
        cofactor: BigUint::from(inst.cofactor),
    })
}

/// `{1, z, …, z^{n'-1}}`, the standard low-degree subspace of
/// `F_{2^n}`.  It is the subspace Gaudry's and Diem's factor bases use
/// and the one the descent is cheapest on; a random subspace would
/// change the constants but not the shape.
fn standard_basis(n: u32, n_prime: u32) -> Vec<F2mElement> {
    (0..n_prime)
        .map(|k| F2mElement::from_bit_positions(&[k], n))
        .collect()
}

/// The `K` or `R` instance at degree `n`.
pub fn instance_for(family: &str, n: u32, seed: u64) -> Option<BinaryInstance> {
    match family {
        // Petit–Quisquater use `y² + xy = x³ + x² + 1`, which is `K_1`.
        "K" => koblitz_instance(1, n).or_else(|| koblitz_instance(0, n)),
        "R" => random_binary_instance(n, seed, 1 << 20),
        _ => None,
    }
}

/// Measure one cell: descend `targets` random targets and run the
/// boolean Gröbner basis on each.
pub fn price_descent_cell(
    family: &str,
    n: u32,
    n_prime: u32,
    summands: u32,
    targets: usize,
    seed: u64,
) -> Option<DescentCell> {
    if n_prime == 0 || n_prime > max_n_prime(summands) {
        return None;
    }
    let inst = instance_for(family, n, seed)?;
    let curve = curve_of(&inst)?;
    let v_basis = standard_basis(n, n_prime);
    let mut rng = StdRng::seed_from_u64(seed ^ ((n as u64) << 32) ^ ((summands as u64) << 16));

    let mut runs: Vec<TargetRun> = Vec::with_capacity(targets);
    let mut bounds: Vec<Option<u32>> = Vec::new();
    let mut brute: Vec<f64> = Vec::new();
    let mut n_vars = 0usize;
    let mut equations = 0usize;

    for _ in 0..targets {
        // A target abscissa drawn uniformly from the field.
        let x_r = inst.gf.to_element(rng.gen::<u64>() & inst.gf.mask);
        let (eqs, vars) = match summands {
            2 => {
                let sys = weil_descend_s3(&curve, &x_r, &v_basis);
                (sys.equations, sys.n_vars)
            }
            3 => {
                let sys = weil_descend_s4(&curve, &x_r, &v_basis);
                (sys.equations, sys.n_vars)
            }
            _ => return None,
        };
        n_vars = vars;
        equations = eqs.len();
        // The derived bound, on the same equations the run will solve.
        bounds.push(semi_regular_degree_of(&eqs, vars));
        // The reference: evaluate every equation at every point of the
        // subspace.  One monomial test per term per point.
        let terms: usize = eqs.iter().map(|p| p.terms.len()).sum();
        brute.push((terms as f64) * 2f64.powi(vars as i32));
        let (gb, stats) = groebner_basis_f2_stats(eqs, vars);
        let inconsistent = gb.len() == 1 && gb[0].terms.len() == 1 && gb[0].terms[0].degree() == 0;
        runs.push(TargetRun { stats, inconsistent });
    }

    let seen: Vec<u32> = bounds.iter().filter_map(|f| *f).collect();
    let d_semireg_min = seen.iter().copied().min();
    let d_av = mean(runs.iter().map(|r| r.stats.solving_degree as f64));
    let ops_mean = mean(runs.iter().map(|r| r.stats.mono_ops as f64));
    let brute_mean = mean(brute.iter().copied());
    Some(DescentCell {
        family: family.into(),
        curve: inst.name.clone(),
        n,
        n_prime,
        summands,
        n_vars,
        equations,
        targets,
        inconsistent: runs.iter().filter(|r| r.inconsistent).count(),
        d_av,
        d_max: runs.iter().map(|r| r.stats.solving_degree).max().unwrap_or(0),
        d_pair_av: mean(runs.iter().map(|r| r.stats.max_pair_degree as f64)),
        d_pair_max: runs.iter().map(|r| r.stats.max_pair_degree).max().unwrap_or(0),
        d_poly_av: mean(runs.iter().map(|r| r.stats.max_poly_degree as f64)),
        d_poly_max: runs.iter().map(|r| r.stats.max_poly_degree).max().unwrap_or(0),
        d_semireg_min,
        d_semireg_max: seen.iter().copied().max(),
        d_semireg_unbounded: bounds.iter().filter(|f| f.is_none()).count(),
        d_av_over_semireg: d_semireg_min.map(|f| d_av / f as f64),
        mono_ops_mean: ops_mean,
        brute_force_ops: brute_mean,
        ops_over_brute_force: ops_mean / brute_mean,
        mono_ops_max: runs.iter().map(|r| r.stats.mono_ops).max().unwrap_or(0),
        spolys_mean: mean(runs.iter().map(|r| r.stats.spolys as f64)),
        basis_len_mean: mean(runs.iter().map(|r| r.stats.basis_len as f64)),
        ms_mean: mean(runs.iter().map(|r| r.stats.wall_ns as f64 / 1e6)),
        ms_max: runs
            .iter()
            .map(|r| r.stats.wall_ns as f64 / 1e6)
            .fold(0.0, f64::max),
        peak_kib_max: runs
            .iter()
            .map(|r| r.stats.peak_bytes() as f64 / 1024.0)
            .fold(0.0, f64::max),
        per_target: runs,
    })
}

/// **The derived boundary: the semi-regular degree of a boolean system.**
///
/// For a semi-regular sequence of equations of degrees `d_1, …, d_k` in
/// `v` boolean variables — that is, one behaving as generically as the
/// field equations `v_i² = v_i` allow — the Hilbert series of the
/// quotient is
///
/// ```text
///     H(t) = (1 + t)^v / Π_i (1 + t^{d_i})
/// ```
///
/// and the degree of regularity is the index of its first non-positive
/// coefficient (Bardet–Faugère–Salvy).  A Gröbner computation on a
/// semi-regular system reaches that degree; one on a system with extra
/// structure falls earlier.  That gap is the whole subject of
/// Petit–Quisquater's Table 2, so this is the column their `D_av` is
/// worth reporting against.
///
/// Returns `None` when the series stays positive out to degree `v`,
/// which is the boolean ring's own ceiling.
pub fn semi_regular_degree(n_vars: usize, degrees: &[u32]) -> Option<u32> {
    if n_vars == 0 || n_vars > 512 {
        return None;
    }
    let top = n_vars + 1;
    // Numerator (1 + t)^v by Pascal's triangle, truncated at t^top.
    let mut num = vec![0i128; top + 1];
    num[0] = 1;
    for _ in 0..n_vars {
        for k in (1..=top).rev() {
            num[k] += num[k - 1];
        }
    }
    // Denominator Π (1 + t^{d_i}), truncated likewise.
    let mut den = vec![0i128; top + 1];
    den[0] = 1;
    for &d in degrees {
        if d == 0 || d as usize > top {
            continue;
        }
        for k in (d as usize..=top).rev() {
            den[k] += den[k - d as usize];
        }
    }
    // Long division: num = den · quot, with den[0] = 1.
    let mut quot = vec![0i128; top + 1];
    for k in 0..=top {
        let mut acc = num[k];
        for j in 1..=k {
            acc -= den[j] * quot[k - j];
        }
        quot[k] = acc;
        if k >= 1 && quot[k] <= 0 {
            return Some(k as u32);
        }
    }
    None
}

/// The derived boundary for one descent system, from the degrees its
/// equations actually have.
fn semi_regular_degree_of(
    eqs: &[crate::cryptanalysis::pq_groebner_f2::F2BoolPoly],
    n_vars: usize,
) -> Option<u32> {
    let degrees: Vec<u32> = eqs
        .iter()
        .filter_map(|p| p.terms.iter().map(|t| t.degree()).max())
        .filter(|d| *d > 0)
        .collect();
    semi_regular_degree(n_vars, &degrees)
}

/// The table, in the shape of Petit–Quisquater's Table 2.
pub fn format_markdown(cells: &[DescentCell]) -> String {
    let mut out = String::new();
    out.push_str("| E | n | n' | m | vars | eqs | D_av | D_pair | D_sr | D_av/D_sr | ops | enumerate | ops/enum | ms | KiB | no decomp |\n");
    out.push_str("|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|\n");
    for c in cells {
        let bound = match (c.d_semireg_min, c.d_semireg_max) {
            (Some(a), Some(b)) if a == b => format!("{a}"),
            (Some(a), Some(b)) => format!("{a}–{b}"),
            _ => "—".into(),
        };
        let ratio = match c.d_av_over_semireg {
            Some(r) => format!("{r:.2}"),
            None => "—".into(),
        };
        out.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {:.1} | {:.1} | {} | {} | {:.3e} | {:.3e} | {:.1}x | {:.1} | {:.0} | {}/{} |\n",
            c.family,
            c.n,
            c.n_prime,
            c.summands,
            c.n_vars,
            c.equations,
            c.d_av,
            c.d_pair_av,
            bound,
            ratio,
            c.mono_ops_mean,
            c.brute_force_ops,
            c.ops_over_brute_force,
            c.ms_mean,
            c.peak_kib_max,
            c.inconsistent,
            c.targets,
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The instrumentation must report a degree that the run could
    /// actually have reached, a basis it actually returned, and a cost
    /// that grows with the work done.
    #[test]
    fn a_descent_cell_reports_a_degree_and_a_cost() {
        let cell = price_descent_cell("K", 11, 4, 2, 3, 7).expect("K_1 over GF(2^11) at n'=4");
        assert_eq!(cell.n_vars, 8, "two summands over a 4-dimensional subspace");
        assert_eq!(cell.equations, 11, "one equation per field coordinate");
        assert_eq!(cell.per_target.len(), 3);
        assert!(cell.d_av >= 1.0, "a processed pair has degree at least one");
        assert!(
            cell.d_max as f64 >= cell.d_av,
            "the max cannot be below the mean"
        );
        assert!(
            cell.d_poly_max >= cell.d_max,
            "an intermediate polynomial is at least as high as the pair degree"
        );
        assert!(cell.mono_ops_mean > 0.0, "the run touched monomials");
        assert!(cell.peak_kib_max > 0.0, "the run held a basis");
    }

    /// Both families must be constructible at the same degree, so the
    /// `K` against `R` comparison the table draws is on matched sizes.
    #[test]
    fn both_curve_families_descend_at_the_same_degree() {
        for family in ["K", "R"] {
            let cell = price_descent_cell(family, 11, 3, 2, 2, 11)
                .unwrap_or_else(|| panic!("{family} at n = 11"));
            assert_eq!(cell.equations, 11);
            assert_eq!(cell.n_vars, 6);
            assert_eq!(cell.family, family);
        }
    }

    /// **The pruning must not lose a solution.**  The chain criterion
    /// skips S-polynomials on the argument that they cannot contribute;
    /// that argument is standard, but it is being applied here to the
    /// systems this table measures, so it is checked on them: the
    /// variety of the Gröbner basis must equal the variety of the
    /// original equations, computed by enumerating every point of the
    /// subspace.  A criterion that dropped a needed pair would show up
    /// as a basis with solutions the system does not have, or the
    /// reverse.
    #[test]
    fn the_pruned_basis_has_exactly_the_solutions_the_system_has() {
        use crate::cryptanalysis::pq_descent::weil_descend_s3;
        use crate::cryptanalysis::pq_groebner_f2::{groebner_basis_f2_stats, solve_system_f2};

        let inst = instance_for("K", 11, 5).expect("K_1 over GF(2^11)");
        let curve = curve_of(&inst).unwrap();
        let v_basis = standard_basis(11, 5);
        let mut rng = StdRng::seed_from_u64(20260922);
        let mut consistent = 0;
        for _ in 0..6 {
            let x_r = inst.gf.to_element(rng.gen::<u64>() & inst.gf.mask);
            let sys = weil_descend_s3(&curve, &x_r, &v_basis);
            let (gb, _) = groebner_basis_f2_stats(sys.equations.clone(), sys.n_vars);
            let from_basis: std::collections::HashSet<u64> =
                solve_system_f2(&gb, sys.n_vars).into_iter().collect();
            // Every point of the subspace, checked against the original
            // equations rather than the basis.
            let direct: std::collections::HashSet<u64> = (0..1u64 << sys.n_vars)
                .filter(|v| sys.equations.iter().all(|e| e.eval(*v) == 0))
                .collect();
            assert_eq!(from_basis, direct, "the pruned basis changed the variety");
            if !direct.is_empty() {
                consistent += 1;
            }
        }
        assert!(
            consistent > 0,
            "every target was inconsistent; the fixture proves only that {{1}} has no roots"
        );
    }

    /// The series is short enough to divide by hand, so it is:
    /// `(1+t)^4 / (1+t²)² = 1 + 4t + 4t² − 4t³ + …`, whose first
    /// non-positive coefficient sits at degree three.
    #[test]
    fn the_semi_regular_degree_is_the_first_non_positive_coefficient() {
        assert_eq!(semi_regular_degree(4, &[2, 2]), Some(3));
        // More equations at the same variable count fall earlier, and
        // the degree never exceeds the boolean ring's own ceiling.
        let square = semi_regular_degree(16, &[2; 16]).unwrap();
        let over = semi_regular_degree(16, &[2; 32]).unwrap();
        assert!(over < square, "{over} should fall before {square}");
        assert!(square <= 16, "a boolean monomial cannot exceed the variable count");
        // Degree-1 equations are linear: they cut the ring down at once.
        assert_eq!(semi_regular_degree(8, &[1; 8]), Some(1));
    }

    /// The truth-table descent caps the subspace dimension; asking past
    /// the cap must decline rather than attempt a `2^{3·6}`-entry table.
    #[test]
    fn the_subspace_dimension_cap_is_enforced() {
        assert_eq!(max_n_prime(2), 8);
        assert_eq!(max_n_prime(3), 5);
        assert!(price_descent_cell("K", 17, 9, 2, 1, 3).is_none());
        assert!(price_descent_cell("K", 17, 6, 3, 1, 3).is_none());
    }
}
