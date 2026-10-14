//! # The greedy coverage selector without enumerating the subgroup.
//!
//! [`super::koblitz_greedy_base`] works in the scalar model `ℤ/rℤ`: it lists
//! every orbit of the prime-order subgroup and keeps a coverage array of
//! length `r`, so it stops at `r = 2²²`.  This module makes the same choice
//! with curve arithmetic and a hashed set, and has no bound on `r` beyond
//! the single-word field (`n ≤ 62`) and the memory of the set.
//!
//! * **Candidates are sampled.**  Each step draws `pool` random subgroup
//!   points (a random abscissa, lifted, multiplied by the cofactor — the
//!   shipped sampler's draw, which is uniform on `⟨G⟩`) and takes the signed
//!   Frobenius orbit of each.
//! * **Coverage is a set of sum abscissae.**  A base `F` closed under
//!   negation covers a target `T` exactly when `T = a + b` for `a, b ∈ F`,
//!   and then `−T = (−a) + (−b)` is covered as well.  Sums therefore come in
//!   `±` pairs sharing an abscissa, so the number of covered targets is
//!   exactly twice the number of distinct abscissae of the sums (a sum equal
//!   to `O` is not a target).  No target set is enumerated or sampled, and
//!   the count is exact.
//! * **A candidate's gain** is the number of its sums `a + b`, over `a` in
//!   the orbit and `b` in the base or the orbit, whose abscissa is not yet
//!   covered.  One batched addition (a single field inversion for the whole
//!   slice) computes `a + b` for every `b` at once.
//!
//! The cost of a step is `pool · |orbit| · |F|` point additions, and the set
//! holds about `|F|²/2` abscissae, so the base is capped at
//! [`MAX_SAMPLED_POINTS`].  That is the limit that matters: the greedy gains
//! on a random base only when collisions among its sums are common, which
//! needs `|F|² ≳ r`, and a base of a few thousand points reaches that only
//! for `r` of a few million to a few hundred million.  Past that the sums
//! cover a vanishing fraction of the subgroup and the two choices are
//! indistinguishable.  `examples/koblitz_greedy_scaling.rs` measures this.

use std::time::Instant;

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::{rngs::StdRng, Rng, SeedableRng};

use crate::binary_ecc::F2mElement;

use super::fx_hash::FxSet;
use super::koblitz_fast::{BatchScratch, FastCurve, FastPoint};
use super::koblitz_greedy_base::{Budget, GreedyBuildReport};
use super::koblitz_index_calculus::{
    build_explicit_frobenius_orbit_factor_base, points_with_x_fast, FrobeniusFactorBase,
    KoblitzCurve,
};

/// Largest base, in points, the sampled selector builds.  The covered-sum
/// set holds about half the square of this many abscissae.
pub const MAX_SAMPLED_POINTS: usize = 8192;

/// Candidate orbits scored per step when the caller asks for the full scan,
/// which does not exist here.
pub const DEFAULT_SAMPLED_POOL: usize = 64;

/// A signed Frobenius orbit of subgroup points.
struct Orbit {
    /// The least abscissa in the orbit; names the orbit.
    key: u64,
    points: Vec<FastPoint>,
}

/// Draws uniform signed orbits of the prime-order subgroup.
struct Sampler {
    curve: FastCurve,
    rng: StdRng,
    cofactor: BigUint,
    field_size: u64,
    n: u32,
    k: u32,
    extension_degree: u32,
}

impl Sampler {
    fn new(kc: &KoblitzCurve, seed: u64) -> Result<Self, String> {
        if kc.n % 2 == 0 {
            return Err("the sampled selector needs odd n".into());
        }
        if !kc.frobenius_is_endomorphism {
            return Err("the sampled selector needs Frobenius to be an endomorphism".into());
        }
        let curve = FastCurve::new(&kc.curve).ok_or("field too wide for single-word arithmetic")?;
        Ok(Self {
            curve,
            rng: StdRng::seed_from_u64(seed),
            cofactor: kc.cofactor.clone(),
            field_size: 1u64 << kc.n,
            n: kc.n,
            k: kc.k,
            extension_degree: kc.extension_degree(),
        })
    }

    /// One random orbit, or `None` after a long run of draws that did not
    /// lift (a degenerate curve).
    fn draw(&mut self) -> Option<Orbit> {
        for _ in 0..4096 {
            let raw = self.rng.gen_range(1..self.field_size);
            let x = F2mElement::from_biguint(&BigUint::from(raw), self.n);
            let Some(lift) = points_with_x_fast(&self.curve, &x).into_iter().next() else {
                continue;
            };
            // [h]P lies in the prime-order subgroup for every P, so no
            // draw is wasted on the wrong coset.
            let projected = self.curve.mul(self.curve.lift(&lift), &self.cofactor);
            if projected.infinity {
                continue;
            }
            let mut points: Vec<FastPoint> = Vec::new();
            let mut key = u64::MAX;
            let mut q = projected;
            for _ in 0..self.extension_degree {
                key = key.min(q.x);
                points.push(q);
                points.push(self.curve.neg(q));
                q = self.curve.frobenius_k(q, self.k);
            }
            points.sort_unstable_by_key(|p| p.pack());
            points.dedup_by_key(|p| p.pack());
            return Some(Orbit { key, points });
        }
        None
    }
}

/// The subgroup order as a target count, saturating at `u64::MAX`.
fn target_count(kc: &KoblitzCurve) -> u64 {
    kc.subgroup_order
        .to_u64()
        .map_or(u64::MAX, |r| r.saturating_sub(1))
}

fn finish(
    kc: &KoblitzCurve,
    keys: &[u64],
    point_count: usize,
    mut report: GreedyBuildReport,
    start: Instant,
) -> Result<(FrobeniusFactorBase, GreedyBuildReport), String> {
    let reps: Vec<F2mElement> = keys
        .iter()
        .map(|&m| F2mElement::from_biguint(&BigUint::from(m), kc.n))
        .collect();
    let fb = build_explicit_frobenius_orbit_factor_base(kc, &reps)
        .ok_or("sampled representatives did not define a factor base")?;
    if fb.points.len() != point_count {
        return Err(format!(
            "built base has {} points, the sampled model has {point_count}",
            fb.points.len()
        ));
    }
    report.points = fb.points.len();
    report.units_chosen = keys.len();
    report.representatives = keys.to_vec();
    report.representatives.sort_unstable();
    report.build_ms = start.elapsed().as_secs_f64() * 1000.0;
    Ok((fb, report))
}

/// Select by the greedy rule over sampled candidate orbits.
///
/// `pool` is the number of orbits drawn and scored per step (`0` means
/// [`DEFAULT_SAMPLED_POOL`]).  The result depends on `seed`.
pub fn sampled_select(
    kc: &KoblitzCurve,
    budget: Budget,
    pool: usize,
    seed: u64,
) -> Result<(FrobeniusFactorBase, GreedyBuildReport), String> {
    let start = Instant::now();
    if let Budget::Points(p) = budget {
        if p > MAX_SAMPLED_POINTS {
            return Err(format!(
                "{p} points exceed the sampled selector's cap of {MAX_SAMPLED_POINTS}"
            ));
        }
    }
    let pool = if pool == 0 {
        DEFAULT_SAMPLED_POOL
    } else {
        pool
    };
    let mut sampler = Sampler::new(kc, seed)?;
    let mut report = GreedyBuildReport {
        targets: target_count(kc) as usize,
        ..Default::default()
    };
    let mut covered: FxSet<u64> = FxSet::default();
    let mut base: Vec<FastPoint> = Vec::new();
    let mut keys: Vec<u64> = Vec::new();
    let mut scratch = BatchScratch::default();
    let mut out: Vec<FastPoint> = Vec::new();

    let done = |base: &[FastPoint], keys: &[u64]| match budget {
        Budget::Points(p) => base.len() >= p,
        Budget::Columns(c) => keys.len() >= c,
    };
    while !done(&base, &keys) {
        let mut best: Option<(usize, Orbit, FxSet<u64>)> = None;
        let mut drawn: FxSet<u64> = FxSet::default();
        for _ in 0..pool {
            let Some(orbit) = sampler.draw() else {
                continue;
            };
            if keys.contains(&orbit.key) || !drawn.insert(orbit.key) {
                continue;
            }
            // Addends: the base and the candidate, so sums of two new
            // points are scored with the sums of a new and an old one.
            let mut addends = base.clone();
            addends.extend_from_slice(&orbit.points);
            let mut fresh: FxSet<u64> = FxSet::default();
            for &a in &orbit.points {
                out.clear();
                sampler.curve.add_many(a, &addends, &mut out, &mut scratch);
                for p in &out {
                    if !p.infinity && !covered.contains(&p.x) {
                        fresh.insert(p.x);
                    }
                }
            }
            report.candidates_scored += 1;
            if best.as_ref().is_none_or(|(g, _, _)| fresh.len() > *g) {
                best = Some((fresh.len(), orbit, fresh));
            }
        }
        let (gain, orbit, fresh) = best.ok_or("no candidate orbit could be drawn")?;
        if gain == 0 {
            report.zero_gain_steps += 1;
        }
        if base.len() + orbit.points.len() > MAX_SAMPLED_POINTS {
            return Err(format!(
                "the next orbit would pass the sampled selector's cap of {MAX_SAMPLED_POINTS} points"
            ));
        }
        covered.extend(fresh);
        base.extend_from_slice(&orbit.points);
        keys.push(orbit.key);
        report.coverage_by_column.push(2 * covered.len());
    }
    report.covered_model = 2 * covered.len();
    finish(kc, &keys, base.len(), report, start)
}

/// The control for the sampled greedy: `columns` orbits drawn uniformly,
/// with no coverage information.
pub fn sampled_random(
    kc: &KoblitzCurve,
    columns: usize,
    seed: u64,
) -> Result<FrobeniusFactorBase, String> {
    if columns == 0 {
        return Err("a factor base needs at least one column".into());
    }
    let mut sampler = Sampler::new(kc, seed)?;
    let mut keys: Vec<u64> = Vec::new();
    let mut points = 0usize;
    while keys.len() < columns {
        let orbit = sampler.draw().ok_or("no orbit could be drawn")?;
        if !keys.contains(&orbit.key) {
            points += orbit.points.len();
            keys.push(orbit.key);
        }
    }
    let report = GreedyBuildReport::default();
    finish(kc, &keys, points, report, Instant::now()).map(|(fb, _)| fb)
}

/// Targets covered by sums of two points of `fb`, and the number of
/// targets, computed from the base alone.
///
/// This does not use the selector's bookkeeping, so it is an independent
/// check of its count, and it scores a base from any source.  It needs
/// `|F|²/2` point additions and a set of about that many abscissae.
pub fn pair_sum_coverage(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
) -> Result<(u64, u64), String> {
    if fb.points.len() > MAX_SAMPLED_POINTS {
        return Err(format!(
            "{} points exceed the pair-sum cap of {MAX_SAMPLED_POINTS}",
            fb.points.len()
        ));
    }
    let curve = FastCurve::new(&kc.curve).ok_or("field too wide for single-word arithmetic")?;
    let points: Vec<FastPoint> = fb.points.iter().map(|p| curve.lift(p)).collect();
    let mut keys: FxSet<u64> = FxSet::default();
    let mut scratch = BatchScratch::default();
    let mut out: Vec<FastPoint> = Vec::new();
    for i in 0..points.len() {
        out.clear();
        curve.add_many(points[i], &points[i..], &mut out, &mut scratch);
        keys.extend(out.iter().filter(|p| !p.infinity).map(|p| p.x));
    }
    Ok((2 * keys.len() as u64, target_count(kc)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::koblitz_factor_base_search::{witness_list, TargetSet};
    use crate::cryptanalysis::koblitz_index_calculus::PairSumTable;

    /// The production census coverage on exhaustive targets.
    fn census_covered(kc: &KoblitzCurve, fb: &FrobeniusFactorBase) -> (usize, usize) {
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
    fn sampled_count_matches_the_census_and_the_independent_pair_sums() {
        // Three counts of the same quantity: the selector's running set,
        // the stand-alone pair-sum count, and the production census over
        // every target of the subgroup.
        let kc = KoblitzCurve::new(0, 13).expect("curve");
        let (fb, report) = sampled_select(&kc, Budget::Columns(3), 16, 1).expect("select");
        let (pairs, targets) = pair_sum_coverage(&kc, &fb).expect("pairs");
        let (census, census_targets) = census_covered(&kc, &fb);
        assert_eq!(targets as usize, census_targets);
        assert_eq!(report.covered_model as u64, pairs);
        assert_eq!(pairs as usize, census);
        assert_eq!(report.coverage_by_column.len(), 3);
        assert_eq!(report.units_chosen, fb.signed_orbits.len());
    }

    #[test]
    fn sampled_greedy_beats_sampled_random_at_equal_columns() {
        let kc = KoblitzCurve::new(1, 17).expect("curve");
        let columns = 8;
        let (g_fb, _) = sampled_select(&kc, Budget::Columns(columns), 32, 1).expect("select");
        let (g, _) = pair_sum_coverage(&kc, &g_fb).expect("greedy");
        let mut sum = 0u64;
        for seed in 1..=10u64 {
            let fb = sampled_random(&kc, columns, seed).expect("random");
            assert_eq!(fb.points.len(), g_fb.points.len());
            sum += pair_sum_coverage(&kc, &fb).expect("random").0;
        }
        assert!(
            g as f64 >= sum as f64 / 10.0,
            "greedy {g} below the random mean {}",
            sum as f64 / 10.0
        );
    }

    #[test]
    fn sampled_selection_is_reproducible_from_its_seed() {
        let kc = KoblitzCurve::new(0, 13).expect("curve");
        let (a, ra) = sampled_select(&kc, Budget::Columns(3), 8, 7).expect("a");
        let (b, rb) = sampled_select(&kc, Budget::Columns(3), 8, 7).expect("b");
        assert_eq!(ra.representatives, rb.representatives);
        assert_eq!(a.points, b.points);
    }

    #[test]
    fn the_point_cap_is_enforced() {
        let kc = KoblitzCurve::new(0, 13).expect("curve");
        assert!(sampled_select(&kc, Budget::Points(MAX_SAMPLED_POINTS + 1), 4, 1).is_err());
    }
}
