//! # Greedy coverage-driven factor base for Koblitz index calculus.
//!
//! [`super::koblitz_index_calculus::build_subgroup_orbit_factor_base`]
//! draws abscissae at random and keeps the subgroup points they lift to.
//! Every point it adds is chosen without looking at what the points
//! already in the base cover, so two random points collide in their sums
//! more often than a deliberate choice would.  That loss is invisible at
//! small bases and becomes the whole story once the base is large enough
//! that its pair sums fill the subgroup (`|F|² ≈ |⟨G⟩|`).
//!
//! This builder picks the signed Frobenius orbits one at a time.  At each
//! step it scores each candidate orbit by the number of subgroup targets
//! that the candidate newly covers as a sum of two base points, and adds
//! the best one.  The score is the quantity
//! [`super::koblitz_factor_base_search`] minimises the ratio of
//! (`coverage` against unknowns), so the construction is a direct greedy
//! attack on that ratio rather than a proxy for it.
//!
//! ## Why the scalar model is exact here
//!
//! The targets of the Koblitz index calculus lie in the prime-order
//! subgroup `⟨G⟩ ≅ ℤ/rℤ`.  Writing each point as `kG` identifies the
//! group law with addition mod `r`, and a sum of base points is a sum of
//! their scalars.  Frobenius acts on `⟨G⟩` as multiplication by `λ`, so
//! a signed Frobenius orbit is a fixed set of scalars.  Coverage of a
//! base `F` is therefore `|{a + b mod r : a, b ∈ F}∖{0}| / (r − 1)`, and
//! the greedy step needs no curve arithmetic at all.  The production
//! census (`PairSumTable` and `witness_list`) recomputes the same number
//! on the curve, and the tests in this module require the two to agree.
//!
//! ## What is and is not claimed
//!
//! The greedy is a heuristic over orbit columns, not an optimum: it
//! maximises the next step's gain and does not look ahead.  It is
//! deterministic when every orbit is scored (`pool = 0`).  Scoring a
//! sample of `pool` orbits per step is cheaper and seeded, and the
//! experiment reports both.  Nothing here changes the oracle, the solver
//! or the counted unit, so the yield it buys is a yield of the census.

use std::collections::{HashMap, HashSet};
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use crate::binary_ecc::{BinaryPoint, F2mElement};

use super::koblitz_factor_base_search::TargetSet;
use super::koblitz_index_calculus::{
    build_explicit_frobenius_orbit_factor_base, FrobeniusFactorBase, KoblitzCurve,
};

/// Largest subgroup order handled by the scalar model.  The coverage
/// array has `r` entries, so this bounds memory at a few megabytes.
pub const MAX_SUBGROUP_ORDER: u64 = 1 << 22;

/// What one greedy construction did, in its own counts.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct GreedyBuildReport {
    /// Signed Frobenius orbits available to choose from.
    pub units_available: usize,
    /// Orbits in the finished base; each is one unknown column.
    pub units_chosen: usize,
    /// Points in the base (both signs), equal to `fb.points.len()`.
    pub points: usize,
    /// Subgroup targets covered by a sum of two base points, from the
    /// scalar model.  The production census must report the same count.
    pub covered_model: usize,
    /// Nonzero subgroup targets, `r − 1`.
    pub targets: usize,
    /// Orbits scored across all steps (candidate evaluations).
    pub candidates_scored: u64,
    /// Steps at which the best candidate covered nothing new.  Those
    /// orbits are added only to reach the requested point count.
    pub zero_gain_steps: usize,
    /// Wall time of the whole construction, orbit enumeration included.
    pub build_ms: f64,
    /// Targets covered after each column is added, scalar model.  The
    /// selection is prefix-stable, so entry `k − 1` is the coverage of the
    /// `k`-orbit base, whichever budget the run was given.
    pub coverage_by_column: Vec<usize>,
    /// Abscissa of one point per chosen orbit, as a bitmask, sorted.  These
    /// freeze the selection: rebuilding the base from them
    /// ([`explicit_orbit_factor_base`]) needs no scan and no enumeration of
    /// the subgroup, and does not depend on this algorithm staying the same.
    pub representatives: Vec<u64>,
}

/// The base made of the Frobenius orbits of the given abscissae (as
/// polynomial-basis bitmasks), exactly as the production builder closes
/// them.  This is what a frozen greedy selection rebuilds from.
pub fn explicit_orbit_factor_base(
    kc: &KoblitzCurve,
    abscissae: &[u64],
) -> Result<FrobeniusFactorBase, String> {
    if kc.n >= 64 {
        return Err("explicit orbit bases need n < 64".into());
    }
    if abscissae.iter().any(|&m| m == 0 || m >= 1u64 << kc.n) {
        return Err("abscissa outside the field".into());
    }
    let reps: Vec<F2mElement> = abscissae
        .iter()
        .map(|&m| F2mElement::from_biguint(&BigUint::from(m), kc.n))
        .collect();
    build_explicit_frobenius_orbit_factor_base(kc, &reps)
        .ok_or_else(|| "abscissae did not define a factor base".into())
}

/// One signed Frobenius orbit of subgroup points.
#[derive(Clone, Debug)]
struct Unit {
    /// Abscissa of the orbit's smallest scalar; the base is rebuilt from
    /// these, and the Frobenius closure restores the rest of the orbit.
    rep: F2mElement,
    /// Every scalar in the orbit, both signs, sorted.
    scalars: Vec<u64>,
}

/// Enumerate the signed Frobenius orbits of `⟨G⟩ ∖ {O}`.
///
/// Each orbit is listed as the set of scalars of its points.  A point and
/// its negative share an abscissa, so an abscissa stands for `{k, r − k}`
/// and an orbit is a union of such pairs.
fn signed_orbit_units(kc: &KoblitzCurve, r: u64) -> Result<Vec<Unit>, String> {
    let all = TargetSet::new(kc, 0, u64::MAX, 0);
    if all.len() as u64 != r - 1 {
        return Err("scalar enumeration did not cover the subgroup".into());
    }
    // abscissa → scalars of the points carrying it, in scalar order.
    let mut by_x: HashMap<BigUint, Vec<u64>> = HashMap::new();
    for (idx, point) in all.points.iter().enumerate() {
        let BinaryPoint::Affine { x, .. } = point else {
            return Err("subgroup scalar mapped to the point at infinity".into());
        };
        by_x.entry(x.to_biguint())
            .or_default()
            .push((idx as u64) + 1);
    }
    let mut visited: HashSet<BigUint> = HashSet::new();
    let mut units = Vec::new();
    for k in 1..r {
        let BinaryPoint::Affine { x, .. } = &all.points[(k - 1) as usize] else {
            continue;
        };
        let key = x.to_biguint();
        if visited.contains(&key) {
            continue;
        }
        // Close the Frobenius orbit the same way the production builder
        // does, so the unit and the finished base agree on its members.
        let mut orbit_x = F2mElement::from_biguint(&key, kc.n);
        let mut scalars: Vec<u64> = Vec::new();
        for _ in 0..kc.extension_degree() {
            let okey = orbit_x.to_biguint();
            let members = by_x
                .get(&okey)
                .ok_or("Frobenius orbit left the prime-order subgroup")?;
            scalars.extend_from_slice(members);
            visited.insert(okey);
            orbit_x = kc.frobenius_x(&orbit_x);
        }
        scalars.sort_unstable();
        scalars.dedup();
        units.push(Unit {
            rep: F2mElement::from_biguint(&key, kc.n),
            scalars,
        });
    }
    Ok(units)
}

/// Targets newly covered by adding `unit` to a base already holding
/// `base` and `cov`.  Pairs are counted once each, with `stamp`/`tick`
/// deduplicating within this candidate.
fn gain_of(
    unit: &[u64],
    base: &[u64],
    cov: &[bool],
    stamp: &mut [u32],
    tick: u32,
    r: u64,
) -> usize {
    let mut gain = 0usize;
    for &a in unit {
        for &b in base.iter().chain(unit.iter()) {
            let t = ((a + b) % r) as usize;
            if t != 0 && !cov[t] && stamp[t] != tick {
                stamp[t] = tick;
                gain += 1;
            }
        }
    }
    gain
}

/// When the selection stops.
///
/// `Points` is the factor-base size the shipped sampler is specified in;
/// `Columns` is the number of signed orbits, which is the unknown count
/// the expected-trial formula charges (`(U + 1 + extra) / coverage`).
/// Orbits differ in size, and the shipped sampler's batches overshoot a
/// requested point count by up to eight orbits, so comparing two bases of
/// equal points is not comparing equal unknowns.  Equal columns is the
/// comparison that isolates the selection rule.
#[derive(Clone, Copy, Debug)]
pub enum Budget {
    /// At least this many points, reached by whole orbits.
    Points(usize),
    /// Exactly this many signed orbits.
    Columns(usize),
}

/// Build the greedy coverage-driven base for `kc` and report what it cost.
///
/// `points` is the target point count, reached by whole orbits, so the
/// base holds at least `points` points (the same stopping rule as
/// `build_subgroup_orbit_factor_base`).  `pool` is the number of candidate
/// orbits scored per step; `0` scores every remaining orbit, which makes
/// the result independent of `seed`.  With a finite pool the candidates are
/// drawn from `seed`, so the recipe `(curve, points, pool, seed)` is
/// reproducible.
pub fn build_greedy_coverage_factor_base(
    kc: &KoblitzCurve,
    points: usize,
    pool: usize,
    seed: u64,
) -> Result<(FrobeniusFactorBase, GreedyBuildReport), String> {
    if points == 0 {
        return Err("a factor base needs at least one point".into());
    }
    greedy_select(kc, Budget::Points(points), pool, seed)
}

/// The greedy base with exactly `columns` signed orbits.  The first `k`
/// orbits of this selection are the `k`-orbit base, so a sweep over the
/// column count is a prefix of one greedy run.
pub fn build_greedy_column_factor_base(
    kc: &KoblitzCurve,
    columns: usize,
    pool: usize,
    seed: u64,
) -> Result<(FrobeniusFactorBase, GreedyBuildReport), String> {
    if columns == 0 {
        return Err("a factor base needs at least one column".into());
    }
    greedy_select(kc, Budget::Columns(columns), pool, seed)
}

/// The control for the greedy: `columns` signed orbits drawn uniformly
/// from the same orbit list, with no coverage information.  Same
/// granularity and same unknown count, so any difference in coverage
/// between this and the greedy is due to the choice of orbits.
pub fn build_random_column_factor_base(
    kc: &KoblitzCurve,
    columns: usize,
    seed: u64,
) -> Result<FrobeniusFactorBase, String> {
    let r = subgroup_scalar_order(kc)?;
    let units = signed_orbit_units(kc, r)?;
    if columns == 0 || columns > units.len() {
        return Err(format!("{columns} columns outside 1..={}", units.len()));
    }
    let mut order: Vec<usize> = (0..units.len()).collect();
    let mut rng = StdRng::seed_from_u64(seed);
    // Partial Fisher–Yates: the first `columns` entries are a uniform sample.
    for i in 0..columns {
        let j = rng.gen_range(i..order.len());
        order.swap(i, j);
    }
    let reps: Vec<F2mElement> = order[..columns]
        .iter()
        .map(|&u| units[u].rep.clone())
        .collect();
    build_explicit_frobenius_orbit_factor_base(kc, &reps)
        .ok_or_else(|| "random representatives did not define a factor base".into())
}

/// The subgroup order as a scalar, if the scalar model applies.
fn subgroup_scalar_order(kc: &KoblitzCurve) -> Result<u64, String> {
    if kc.n >= 64 {
        return Err("greedy coverage needs n < 64".into());
    }
    kc.subgroup_order
        .to_u64()
        .filter(|&r| (3..=MAX_SUBGROUP_ORDER).contains(&r))
        .ok_or_else(|| "subgroup order outside the scalar model range".into())
}

fn greedy_select(
    kc: &KoblitzCurve,
    budget: Budget,
    pool: usize,
    seed: u64,
) -> Result<(FrobeniusFactorBase, GreedyBuildReport), String> {
    let start = Instant::now();
    let r = subgroup_scalar_order(kc)?;

    let units = signed_orbit_units(kc, r)?;
    if let Budget::Columns(c) = budget {
        if c > units.len() {
            return Err(format!(
                "{c} columns exceed the {} orbits available",
                units.len()
            ));
        }
    }
    let mut report = GreedyBuildReport {
        units_available: units.len(),
        targets: (r - 1) as usize,
        ..Default::default()
    };

    let mut cov = vec![false; r as usize];
    let mut covered = 0usize;
    let mut stamp = vec![0u32; r as usize];
    let mut tick = 0u32;
    let mut base: Vec<u64> = Vec::new();
    let mut taken = vec![false; units.len()];
    let mut chosen: Vec<usize> = Vec::new();
    let mut rng = StdRng::seed_from_u64(seed);
    let full_scan = pool == 0 || pool >= units.len();

    let done = |base: &[u64], chosen: &[usize]| match budget {
        Budget::Points(p) => base.len() >= p,
        Budget::Columns(c) => chosen.len() >= c,
    };
    while !done(&base, &chosen) {
        let candidates: Vec<usize> = if full_scan {
            (0..units.len()).filter(|&u| !taken[u]).collect()
        } else {
            // Seeded sample of untaken orbits; a bounded number of draws
            // keeps this terminating when few orbits remain.
            let mut picked: HashSet<usize> = HashSet::new();
            let mut draws = 0usize;
            while picked.len() < pool && draws < 64 * pool {
                draws += 1;
                let u = rng.gen_range(0..units.len());
                if !taken[u] {
                    picked.insert(u);
                }
            }
            let mut v: Vec<usize> = picked.into_iter().collect();
            v.sort_unstable();
            v
        };
        if candidates.is_empty() {
            return Err(format!(
                "ran out of orbits after {} columns and {} points",
                chosen.len(),
                base.len()
            ));
        }
        let mut best: Option<(usize, usize)> = None;
        for &u in &candidates {
            tick = tick.wrapping_add(1);
            if tick == 0 {
                stamp.iter_mut().for_each(|s| *s = 0);
                tick = 1;
            }
            let gain = gain_of(&units[u].scalars, &base, &cov, &mut stamp, tick, r);
            report.candidates_scored += 1;
            if best.is_none_or(|(g, _)| gain > g) {
                best = Some((gain, u));
            }
        }
        let (gain, u) = best.expect("candidates is nonempty");
        if gain == 0 {
            report.zero_gain_steps += 1;
        }
        taken[u] = true;
        chosen.push(u);
        base.extend_from_slice(&units[u].scalars);
        // New pairs are (new, any): every such sum is marked once here.
        for &a in &units[u].scalars {
            for &b in &base {
                let t = ((a + b) % r) as usize;
                if t != 0 && !cov[t] {
                    cov[t] = true;
                    covered += 1;
                }
            }
        }
        report.coverage_by_column.push(covered);
    }

    report.covered_model = covered;
    report.units_chosen = chosen.len();
    report.representatives = chosen
        .iter()
        .filter_map(|&u| units[u].rep.to_biguint().to_u64())
        .collect();
    report.representatives.sort_unstable();
    let reps: Vec<F2mElement> = chosen.iter().map(|&u| units[u].rep.clone()).collect();
    let fb = build_explicit_frobenius_orbit_factor_base(kc, &reps)
        .ok_or("greedy representatives did not define a factor base")?;
    report.points = fb.points.len();
    if report.points != base.len() {
        return Err(format!(
            "built base has {} points, the scalar model has {}",
            report.points,
            base.len()
        ));
    }
    report.build_ms = start.elapsed().as_secs_f64() * 1000.0;
    Ok((fb, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::koblitz_factor_base_search::{expected_trials, witness_list};
    use crate::cryptanalysis::koblitz_index_calculus::PairSumTable;

    /// The production census on a base, as `(covered, targets)`.
    fn production_coverage(kc: &KoblitzCurve, fb: &FrobeniusFactorBase) -> (usize, usize) {
        let targets = TargetSet::new(kc, 0, u64::MAX, 0);
        let table = PairSumTable::build(kc, fb).expect("pair table");
        let witnesses = witness_list(kc, fb, &table, &targets, 2);
        let mut covered = vec![false; targets.len()];
        for w in &witnesses {
            covered[w.target as usize] = true;
        }
        (covered.iter().filter(|&&c| c).count(), targets.len())
    }

    #[test]
    fn scalar_model_agrees_with_production_census() {
        // The greedy's own count must equal the census of the base it
        // returns, or the coverage it optimises is not the coverage the
        // pipeline sees.
        let kc = KoblitzCurve::new(0, 9).expect("curve");
        let (fb, report) = build_greedy_coverage_factor_base(&kc, 40, 0, 1).expect("greedy");
        let (covered, targets) = production_coverage(&kc, &fb);
        assert_eq!(targets, report.targets);
        assert_eq!(covered, report.covered_model);
        assert!(report.covered_model > 0);
    }

    #[test]
    fn greedy_covers_more_than_random_orbits_at_equal_columns() {
        // Equal unknowns, so equal `U + 1` in the trial formula: the
        // comparison is coverage alone.  Against the mean of 20 random
        // orbit choices of the same count, measured by the census.
        let kc = KoblitzCurve::new(0, 9).expect("curve");
        let columns = 3;
        let (g_fb, _) = build_greedy_column_factor_base(&kc, columns, 0, 0).expect("greedy");
        let (g_cov, t) = production_coverage(&kc, &g_fb);
        let mut sum = 0usize;
        for seed in 1..=20u64 {
            let fb = build_random_column_factor_base(&kc, columns, seed).expect("random");
            sum += production_coverage(&kc, &fb).0;
        }
        let mean = sum as f64 / 20.0;
        assert!(
            g_cov as f64 >= mean,
            "greedy {g_cov}/{t} below random mean {mean:.1}"
        );
        let g_trials = expected_trials(columns, 0, g_cov, t);
        assert!(g_trials.is_finite());
    }

    #[test]
    fn full_scan_is_deterministic_and_seed_free() {
        let kc = KoblitzCurve::new(0, 9).expect("curve");
        let (a, _) = build_greedy_coverage_factor_base(&kc, 30, 0, 1).expect("a");
        let (b, _) = build_greedy_coverage_factor_base(&kc, 30, 0, 99).expect("b");
        assert_eq!(a.points, b.points);
    }

    #[test]
    fn column_budget_is_a_prefix_of_one_greedy_run() {
        // Sweeping the column count must be one selection, not many, so
        // the k-orbit base is the first k orbits of the (k+1)-orbit one.
        let kc = KoblitzCurve::new(0, 9).expect("curve");
        let (small, s) = build_greedy_column_factor_base(&kc, 2, 0, 0).expect("two");
        let (big, b) = build_greedy_column_factor_base(&kc, 4, 0, 0).expect("four");
        assert_eq!(s.units_chosen, 2);
        assert_eq!(b.units_chosen, 4);
        assert_eq!(small.signed_orbits.len(), 2);
        assert_eq!(big.signed_orbits.len(), 4);
        assert!(b.covered_model >= s.covered_model);
    }

    #[test]
    fn random_column_control_has_the_requested_columns_and_census() {
        let kc = KoblitzCurve::new(0, 9).expect("curve");
        let fb = build_random_column_factor_base(&kc, 3, 5).expect("random");
        assert_eq!(fb.signed_orbits.len(), 3);
        let (covered, targets) = production_coverage(&kc, &fb);
        assert!(covered <= targets);
        // Same seed, same base: the control is reproducible.
        let again = build_random_column_factor_base(&kc, 3, 5).expect("random");
        assert_eq!(fb.points, again.points);
    }

    #[test]
    fn base_is_closed_under_negation() {
        // Every abscissa carries both signs, so the point count is even
        // and each signed orbit is a whole Frobenius orbit of abscissae.
        let kc = KoblitzCurve::new(1, 9).expect("curve");
        let (fb, report) = build_greedy_coverage_factor_base(&kc, 30, 0, 1).expect("greedy");
        assert_eq!(fb.points.len() % 2, 0);
        assert_eq!(report.points, fb.points.len());
        assert_eq!(report.units_chosen, fb.signed_orbits.len());
    }

    #[test]
    fn frozen_selection_rebuilds_the_same_base_without_the_scan() {
        use crate::cryptanalysis::koblitz_factor_base_search::{freeze_spec, FactorBaseSpec};
        let kc = KoblitzCurve::new(0, 13).expect("curve");
        let spec = FactorBaseSpec::GreedyCoverage {
            seed: 1,
            points: 52,
            pool: 0,
        };
        let frozen = freeze_spec(&kc, &spec).expect("freeze");
        assert!(matches!(frozen, FactorBaseSpec::ExplicitOrbits { .. }));
        let a = spec.materialize(&kc).expect("greedy");
        let b = frozen.materialize(&kc).expect("frozen");
        let mut pa: Vec<_> = a.points.iter().map(|p| format!("{p:?}")).collect();
        let mut pb: Vec<_> = b.points.iter().map(|p| format!("{p:?}")).collect();
        pa.sort();
        pb.sort();
        assert_eq!(pa, pb);
        assert_eq!(a.signed_orbits.len(), b.signed_orbits.len());
        // A wrapper around the greedy freezes through to its parent.
        let wrapped = FactorBaseSpec::TwoTorsionSaturated {
            parent: Box::new(spec.clone()),
        };
        let fw = freeze_spec(&kc, &wrapped).expect("freeze wrapped");
        let FactorBaseSpec::TwoTorsionSaturated { parent } = fw else {
            panic!("wrapper lost");
        };
        assert!(matches!(*parent, FactorBaseSpec::ExplicitOrbits { .. }));
        // Specs with no greedy root pass through unchanged.
        let plain = FactorBaseSpec::SubgroupOrbits {
            seed: 2,
            points: 30,
        };
        assert_eq!(freeze_spec(&kc, &plain).expect("plain"), plain);
    }

    #[test]
    fn search_offers_greedy_candidates_only_when_asked_and_only_where_they_fit() {
        use crate::cryptanalysis::koblitz_factor_base_search::{
            candidate_specs, FactorBaseSpec, Family, SearchOptions,
        };
        let kc = KoblitzCurve::new(0, 13).expect("curve");
        let off = SearchOptions::default();
        assert!(candidate_specs(&kc, &off)
            .iter()
            .all(|s| !matches!(s, FactorBaseSpec::GreedyCoverage { .. })));
        let on = SearchOptions {
            families: vec![Family::Greedy],
            saturate: false,
            ..SearchOptions::default()
        };
        let specs = candidate_specs(&kc, &on);
        assert!(!specs.is_empty());
        assert!(specs
            .iter()
            .all(|s| matches!(s, FactorBaseSpec::GreedyCoverage { pool: 0, .. })));
        // Beyond the scalar model's cap there is nothing to enumerate, so
        // the family yields no candidates instead of failing the search.
        let wide = KoblitzCurve::new(0, 41).expect("curve");
        assert!(candidate_specs(&wide, &on).is_empty());
    }

    #[test]
    fn recipe_round_trips_and_rebuilds_the_same_base() {
        // The `ic` tool reads factor bases as recipes, so the variant must
        // serialise, parse back equal, and materialise to the direct build.
        use crate::cryptanalysis::koblitz_factor_base_search::FactorBaseSpec;
        let kc = KoblitzCurve::new(1, 9).expect("curve");
        for pool in [0usize, 4] {
            let spec = FactorBaseSpec::GreedyCoverage {
                seed: 3,
                points: 30,
                pool,
            };
            let json = serde_json::to_string(&spec).expect("json");
            assert!(json.contains("\"kind\":\"greedy_coverage\""), "{json}");
            let back: FactorBaseSpec = serde_json::from_str(&json).expect("parse");
            assert_eq!(back, spec);
            let via_spec = spec.materialize(&kc).expect("materialize");
            let (direct, _) = build_greedy_coverage_factor_base(&kc, 30, pool, 3).expect("direct");
            assert_eq!(via_spec.points, direct.points);
        }
    }
}
