//! # Large-prime collector for the Frobenius pair-table oracle.
//!
//! The classic three-summand oracle (`mitm-frobenius:m=3`) walks the
//! base `F`, forms `T = R − P_i` for every signed base point and probes
//! the Frobenius-folded pair table for `T = P_j + P_k`.  This collector
//! keeps that loop and that table unchanged and adds a **large-prime
//! domain** `D`: a set of signed Frobenius orbits of the prime-order
//! subgroup, disjoint from `F`, of about `lp_domain_factor · |F|` signed
//! points — concretely, the next orbits the archived compact-orbit raw-x
//! scan would have taken after the base's.  Three partial shapes are
//! accepted on the same `T`, all with the trial's own `(a, b)`:
//!
//! | shape | found by | small | large |
//! |:--|:--|--:|--:|
//! | `R = P_i + L` | one orbit-key lookup of `T` in `D` | 1 | 1 |
//! | `R = P_i + P_j + L` | the large table, pair `F × D` | 2 | 1 |
//! | `R = P_i + L_1 + L_2` | the large table, pair `D × D` (`large_primes = 2` only) | 1 | 2 |
//!
//! The large table is a second [`FrobeniusPairTable`] over the envelope
//! `F ∪ D`, built with the `F × F` pairs left out (those are the classic
//! table's) and, with one large prime, the `D × D` pairs left out too.
//! Its build is charged to the oracle's setup phase like the classic
//! table's; a probe is one canonicalisation (shared by all three
//! lookups of a `T`), one lookup per table and one addition on a hit.
//!
//! Partials go to [`WideLargePrimeEliminator`], which proves each one as
//! a point equation in the subgroup, names each residual by its signed
//! Frobenius orbit with the exact `±λ^t` coefficient, and reduces
//! against stored pivots modulo `r`; a row with no residual left is a
//! **complete relation** with its own merged `(a, b)`.  Before such a
//! row reaches the matrix the collector checks it again, independently,
//! in the fast counted group: `[a]G + [b]Q = Σ_col [c_col]·U_col` with
//! `U_col` the column's unit point (`coef_of = 1`).  Those scalar
//! multiplications are charged to the relation phase.  A trial yields at
//! most one complete relation, direct or merged, as the classic oracle
//! does; the partials it stored before that stay in the graph.
//!
//! What this is for is a measurement (`docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md`
//! §7 item 4): whether bounded one- and two-large-prime partials change
//! the trials and operations per complete relation on the compact-orbit
//! pipeline.  It is a constant-factor device inside the pair-table
//! family and makes no claim beyond that.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::ToPrimitive;

use super::stages::{InstanceCtx, Params};
use crate::cryptanalysis::ic_boundary::{
    BinaryGroup, BinaryInstance, CountedGroup, FactorBase, FrobeniusPairTable, GroupOps,
    OracleCounters, Relation,
};
use crate::cryptanalysis::koblitz_fast::FastPoint;
use crate::cryptanalysis::koblitz_index_calculus::build_explicit_frobenius_point_factor_base;
use crate::cryptanalysis::koblitz_large_prime::{
    FeedOutcome, LargePrimeLimits, PartialCandidate, WideFullRelation, WideLargePrimeEliminator,
};

/// What the oracle parameters ask for.
#[derive(Clone, Copy, Debug)]
pub struct LargePrimeConfig {
    /// `1` or `2`.
    pub large_primes: u32,
    /// `|D| / |F|` in signed points; `D` is rounded to whole orbits.
    pub domain_factor: f64,
    /// Inclusive raw-abscissa cap of the scan that materialises `D`.
    pub raw_x_cap: u64,
    pub limits: LargePrimeLimits,
}

impl LargePrimeConfig {
    /// Read `large_primes`, `lp_domain_factor` and `lp_raw_x_cap` from
    /// the oracle's parameters; `None` when `large_primes` is absent or
    /// zero, which is the classic oracle.
    pub fn from_params(params: &Params) -> Result<Option<Self>, String> {
        let large_primes = params.u64_or("large_primes", 0)?;
        if large_primes == 0 {
            return Ok(None);
        }
        if large_primes > 2 {
            return Err("large_primes must be 0, 1 or 2".into());
        }
        let domain_factor = match params.get("lp_domain_factor") {
            None => 10.0,
            Some(raw) => raw
                .parse::<f64>()
                .map_err(|_| format!("parameter `lp_domain_factor` is not a number: `{raw}`"))?,
        };
        // NaN fails every comparison, so it is named.
        if domain_factor.is_nan() || domain_factor <= 0.0 || !domain_factor.is_finite() {
            return Err("lp_domain_factor must be a positive number".into());
        }
        let raw_x_cap = params.u64_or("lp_raw_x_cap", 1 << 40)?;
        Ok(Some(Self {
            large_primes: large_primes as u32,
            domain_factor,
            raw_x_cap,
            limits: LargePrimeLimits {
                max_merge_depth: 4096,
                max_row_terms: 1 << 16,
                max_pivots: 1 << 24,
            },
        }))
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Counters {
    partials_one_large: u64,
    partials_mixed: u64,
    partials_two_large: u64,
    partials_fed: u64,
    merged_relations: u64,
    direct_two_sum: u64,
    direct_three_sum: u64,
    capped: u64,
    rejected: u64,
    trivial_rows: u64,
    target_only_rows: u64,
    verification_scalar_mults: u64,
    verification_failures: u64,
    merge_depth_total: u64,
    merge_sources_total: u64,
    filter_wall_ns: u64,
}

/// What an orbit key of `T` names.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Membership {
    Base,
    Domain,
}

/// The collector: the classic table, the large table over `F ∪ D`, the
/// exact eliminator and the column map back to the framework base.
pub struct LargePrimeCollector<'i> {
    instance: &'i BinaryInstance,
    config: LargePrimeConfig,
    /// `|F|`: envelope indices below this are base points.
    small_points: usize,
    /// `F ∪ D` as one base; the `D` columns are placeholders the loop
    /// never writes a row over.
    envelope: FactorBase<FastPoint>,
    /// The classic table over `F`, identical to the plain oracle's.
    table_small: FrobeniusPairTable,
    /// The envelope table without `F × F` (and, with one large prime,
    /// without `D × D`).
    table_large: FrobeniusPairTable,
    membership: HashMap<u64, Membership>,
    domain_orbits: usize,
    eliminator: WideLargePrimeEliminator<'i>,
    /// Signed orbit of the eliminator's base → (framework column, the
    /// representative's `coef_of`).
    orbit_column: Vec<(usize, u64)>,
    /// Per framework column, a point with `coef_of = 1`.
    column_unit: Vec<usize>,
    next_source: u64,
    counters: Counters,
    setup: BTreeMap<String, u64>,
    last_error: Option<String>,
}

impl<'i> LargePrimeCollector<'i> {
    /// Build the domain, both tables and the eliminator for one run.
    /// All group work is charged to `ops` (the oracle's setup phase).
    pub fn new(
        instance: &'i BinaryInstance,
        ctx: &InstanceCtx<BinaryGroup<'_>>,
        fb: &FactorBase<FastPoint>,
        config: LargePrimeConfig,
        ops: &mut GroupOps,
    ) -> Result<Self, String> {
        let kc = instance
            .koblitz
            .as_ref()
            .ok_or("large primes need a Koblitz instance")?;
        if ctx.r != instance.r || ctx.cofactor != instance.cofactor {
            return Err("large-prime context does not match its instance".into());
        }
        if instance.fast.lower(ctx.generator) != *kc.generator() {
            return Err("large primes need the instance's own generator".into());
        }
        if !(1..=2).contains(&config.large_primes) {
            return Err("large_primes must be 1 or 2".into());
        }
        if fb.points.is_empty() || fb.columns == 0 {
            return Err("large primes need a non-empty factor base".into());
        }
        let g = ctx.group;
        let n = instance.n as usize;
        let r = instance.r;

        // One unit point per column, and the subgroup check on it: the
        // base is closed under the signed Frobenius (and the multiplier
        // closure), so a column's points all have the unit's order.
        let mut column_unit = vec![usize::MAX; fb.columns];
        for (i, (&col, &coef)) in fb.col_of.iter().zip(&fb.coef_of).enumerate() {
            if coef == 1 && column_unit[col] == usize::MAX {
                column_unit[col] = i;
            }
        }
        if let Some(col) = column_unit.iter().position(|&u| u == usize::MAX) {
            return Err(format!(
                "large primes need a unit point (coef_of = 1) in every column; column {col} has none"
            ));
        }
        for &u in &column_unit {
            if !g.is_identity(&g.mul(ops, fb.points[u], r)) {
                return Err(
                    "large primes need a factor base inside the prime-order subgroup".into(),
                );
            }
        }

        // The classic table first: with it the direct path is the plain
        // oracle's, probe for probe.
        let table_small = FrobeniusPairTable::build(instance, fb)
            .ok_or("this instance has no Koblitz structure to fold by")?;
        ops.merge(table_small.build_ops);

        // Orbit keys of the base (one per signed orbit): the scan below
        // skips them and the probe recognises `T ∈ F`.
        let mut membership: HashMap<u64, Membership> = HashMap::new();
        let mut base_keys = 0u64;
        for p in &fb.points {
            if membership
                .insert(table_small.orbit_key(p.x), Membership::Base)
                .is_none()
            {
                base_keys += 1;
            }
        }

        // The domain: the next projected signed Frobenius orbits of the
        // raw-x scan that are not the base's, up to the requested size.
        let wanted_orbits = ((config.domain_factor * fb.points.len() as f64) / (2 * n) as f64)
            .round()
            .max(1.0) as usize;
        let mut domain_points: Vec<FastPoint> = Vec::with_capacity(wanted_orbits * 2 * n);
        let mut domain_col: Vec<usize> = Vec::with_capacity(wanted_orbits * 2 * n);
        let mut domain_coef: Vec<u64> = Vec::with_capacity(wanted_orbits * 2 * n);
        let lambda = kc
            .lambda
            .to_u64()
            .ok_or("Frobenius eigenvalue does not fit in u64")?;
        let mut domain_orbits = 0usize;
        let (mut scanned, mut lifts, mut projected, mut frobenius_maps) = (0u64, 0u64, 0u64, 0u64);
        let max_x = config.raw_x_cap.min((1u64 << instance.n) - 1);
        let mut raw_x = 1u64;
        while domain_orbits < wanted_orbits && raw_x <= max_x {
            scanned += 1;
            let candidates = instance.points_with_x(raw_x);
            raw_x += 1;
            let Some(&lift) = candidates.first() else {
                continue;
            };
            lifts += 1;
            let seed = g.mul(ops, lift, instance.cofactor);
            if g.is_identity(&seed) || seed.x <= 1 {
                continue;
            }
            projected += 1;
            let key = table_small.orbit_key(seed.x);
            if membership.contains_key(&key) {
                continue;
            }
            if !g.is_identity(&g.mul(ops, seed, r)) {
                return Err(format!(
                    "projected point from raw x {} has wrong order",
                    raw_x - 1
                ));
            }
            membership.insert(key, Membership::Domain);
            let column = fb.columns + domain_orbits;
            let mut member = seed;
            let mut coefficient = 1u64;
            for _ in 0..n {
                domain_points.push(member);
                domain_col.push(column);
                domain_coef.push(coefficient);
                domain_points.push(g.neg(member));
                domain_col.push(column);
                domain_coef.push((r - coefficient) % r);
                member = instance.fast.frobenius_k(member, 1);
                frobenius_maps += 1;
                coefficient = ((coefficient as u128 * lambda as u128) % r as u128) as u64;
            }
            if member != seed {
                return Err(format!(
                    "Frobenius orbit from raw x {} is not closed",
                    raw_x - 1
                ));
            }
            domain_orbits += 1;
        }
        if domain_orbits < wanted_orbits {
            return Err(format!(
                "lp_raw_x_cap {} exhausted after {scanned} abscissae with {domain_orbits} of {wanted_orbits} domain orbits",
                config.raw_x_cap
            ));
        }
        let domain_point_count = domain_points.len();

        // The envelope and the large table over it.
        let small_points = fb.points.len();
        let mut points = fb.points.clone();
        points.extend(domain_points);
        let mut col_of = fb.col_of.clone();
        col_of.extend(domain_col);
        let mut coef_of = fb.coef_of.clone();
        coef_of.extend(domain_coef);
        let fast = &instance.fast;
        let envelope = FactorBase::from_column_map(
            format!(
                "{} plus {domain_orbits} large-prime orbits ({domain_point_count} signed points)",
                fb.description
            ),
            points,
            col_of,
            coef_of,
            fb.columns + domain_orbits,
            |p| p.pack(),
            |p| fast.neg(*p).pack(),
            |p| p.x,
        )?;
        let two = config.large_primes == 2;
        let table_large = FrobeniusPairTable::build_filtered(instance, &envelope, |rep, j| {
            let rep_large = rep >= small_points;
            let j_large = j >= small_points;
            (rep_large || j_large) && (two || !(rep_large && j_large))
        })
        .ok_or("the envelope has no Koblitz structure to fold by")?;
        ops.merge(table_large.build_ops);
        // Both canonicalisations must name orbits the same way, or the
        // shared key would probe the wrong entry.
        for p in fb.points.iter().take(8) {
            if table_small.orbit_key(p.x) != table_large.orbit_key(p.x) {
                return Err("the two pair tables canonicalise differently".into());
            }
        }

        // The eliminator's signed-orbit base over the same points, in the
        // same order, so a framework index is an eliminator index.
        let lowered: Vec<_> = fb.points.iter().map(|&p| fast.lower(p)).collect();
        let ffb = build_explicit_frobenius_point_factor_base(kc, lowered.clone())
            .ok_or("the factor base is not a complete signed-Frobenius base")?;
        if ffb.points != lowered {
            return Err("the signed-orbit base reordered the factor base".into());
        }
        let orbit_column: Vec<(usize, u64)> = ffb
            .signed_orbits
            .iter()
            .map(|orbit| {
                let rep = orbit[0];
                (fb.col_of[rep], fb.coef_of[rep])
            })
            .collect();
        let eliminator =
            WideLargePrimeEliminator::new_owned(kc, ffb, fast.lower(ctx.target), config.limits)
                .map_err(|e| format!("large-prime eliminator: {e}"))?;

        let mut setup = BTreeMap::new();
        setup.insert("lp_domain_orbits".into(), domain_orbits as u64);
        setup.insert("lp_domain_points".into(), domain_point_count as u64);
        setup.insert("lp_base_orbit_keys".into(), base_keys);
        setup.insert("lp_domain_scan_abscissae_uncharged".into(), scanned);
        setup.insert("lp_domain_scan_lifts_uncharged".into(), lifts);
        setup.insert("lp_domain_scan_projected".into(), projected);
        setup.insert("lp_domain_frobenius_maps".into(), frobenius_maps);
        setup.insert("lp_table_entries".into(), table_large.entries);
        setup.insert("lp_table_adds".into(), table_large.build_ops.adds);
        setup.insert(
            "lp_table_canonicalisations".into(),
            table_large.build_canonicalisations,
        );
        setup.insert("lp_table_lookups".into(), table_large.build_lookups);
        setup.insert(
            "lp_table_frobenius_maps".into(),
            table_large.build_frobenius_maps,
        );
        setup.insert("pair_table_entries".into(), table_small.entries);
        setup.insert(
            "canonicalisations".into(),
            table_small.build_canonicalisations,
        );
        setup.insert("frobenius_maps".into(), table_small.build_frobenius_maps);
        setup.insert("lookups".into(), table_small.build_lookups);
        Ok(Self {
            instance,
            config,
            small_points,
            envelope,
            table_small,
            table_large,
            membership,
            domain_orbits,
            eliminator,
            orbit_column,
            column_unit,
            next_source: 1,
            counters: Counters::default(),
            setup,
            last_error: None,
        })
    }

    pub fn config(&self) -> LargePrimeConfig {
        self.config
    }

    pub fn domain_orbits(&self) -> usize {
        self.domain_orbits
    }

    pub fn domain_points(&self) -> usize {
        self.envelope.points.len() - self.small_points
    }

    /// Native setup counters, for the oracle's `setup_native`.
    pub fn setup_counters(&self) -> BTreeMap<String, u64> {
        self.setup.clone()
    }

    /// The first error the eliminator returned, if any; a non-empty value
    /// means a partial the collector formed did not pass the exact check.
    pub fn last_error(&self) -> Option<&str> {
        self.last_error.as_deref()
    }

    /// Native relation-phase counters, for the oracle's `relation_counters`.
    pub fn relation_counters(&self) -> BTreeMap<String, u64> {
        let c = &self.counters;
        let stats = self.eliminator.stats();
        let mut out = BTreeMap::new();
        out.insert("lp_partials_one_large".into(), c.partials_one_large);
        out.insert("lp_partials_mixed".into(), c.partials_mixed);
        out.insert("lp_partials_two_large".into(), c.partials_two_large);
        out.insert("lp_partials_fed".into(), c.partials_fed);
        out.insert("lp_partials_stored".into(), stats.pivot_rows as u64);
        out.insert("lp_pivot_reductions".into(), stats.pivot_reductions);
        out.insert("lp_merged_relations".into(), c.merged_relations);
        out.insert("lp_merge_depth_total".into(), c.merge_depth_total);
        out.insert("lp_merge_sources_total".into(), c.merge_sources_total);
        out.insert("lp_direct_two_sum".into(), c.direct_two_sum);
        out.insert("lp_direct_three_sum".into(), c.direct_three_sum);
        out.insert("lp_capped".into(), c.capped);
        out.insert("lp_rejected".into(), c.rejected);
        out.insert("lp_trivial_rows".into(), c.trivial_rows);
        out.insert("lp_target_only_rows".into(), c.target_only_rows);
        out.insert(
            "lp_verification_scalar_mults".into(),
            c.verification_scalar_mults,
        );
        out.insert("lp_verification_failures".into(), c.verification_failures);
        out.insert("lp_filter_wall_ns".into(), c.filter_wall_ns);
        out
    }

    /// Feed one partial; `Some(row)` when the graph closed on it and the
    /// row passed both checks.
    fn feed(
        &mut self,
        g: &BinaryGroup<'_>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        a: u64,
        b: u64,
        small: Vec<usize>,
        residual: Vec<FastPoint>,
    ) -> Option<Relation> {
        let instance = self.instance;
        let candidate = PartialCandidate {
            source_id: self.next_source,
            coef_a: BigUint::from(a),
            coef_b: BigUint::from(b),
            small_point_indices: small,
            residual_points: residual.iter().map(|&p| instance.fast.lower(p)).collect(),
        };
        self.next_source += 1;
        self.counters.partials_fed += 1;
        let started = Instant::now();
        let outcome = self.eliminator.feed(candidate);
        self.counters.filter_wall_ns += started.elapsed().as_nanos() as u64;
        match outcome {
            Err(e) => {
                self.counters.rejected += 1;
                if self.last_error.is_none() {
                    self.last_error = Some(e.to_string());
                }
                None
            }
            Ok(FeedOutcome::PivotStored { .. }) => None,
            Ok(FeedOutcome::UnknownCap { .. }) => {
                self.counters.capped += 1;
                None
            }
            Ok(FeedOutcome::Full(full)) => self.complete(g, fb, ops, full),
        }
    }

    /// Convert a completed eliminator row to a framework row and check it
    /// in the counted group.
    fn complete(
        &mut self,
        g: &BinaryGroup<'_>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        full: WideFullRelation,
    ) -> Option<Relation> {
        let instance = self.instance;
        let r = instance.r;
        let a = full.coef_a.to_u64()?;
        let b = full.coef_b.to_u64()?;
        let mut columns: BTreeMap<usize, u64> = BTreeMap::new();
        for (&orbit, coefficient) in &full.small {
            let (col, rep_coef) = *self.orbit_column.get(orbit)?;
            let c = coefficient.to_u64()?;
            let weighted = ((c as u128 * rep_coef as u128) % r as u128) as u64;
            let slot = columns.entry(col).or_insert(0);
            *slot = ((*slot as u128 + weighted as u128) % r as u128) as u64;
        }
        columns.retain(|_, c| *c != 0);
        if columns.is_empty() && a == 0 && b == 0 {
            self.counters.trivial_rows += 1;
            return None;
        }
        // Independent replay in the counted group: [a]G + [b]Q = Σ [c]·U_col.
        let generator = instance.generator;
        let target = instance.fast.lift(self.eliminator.target());
        let lhs_a = g.mul(ops, generator, a);
        let lhs_b = g.mul(ops, target, b);
        let lhs = g.add(ops, lhs_a, lhs_b);
        let mut rhs = g.identity();
        for (&col, &c) in &columns {
            let term = g.mul(ops, fb.points[self.column_unit[col]], c);
            rhs = g.add(ops, rhs, term);
        }
        self.counters.verification_scalar_mults += 2 + columns.len() as u64;
        if lhs != rhs {
            self.counters.verification_failures += 1;
            if self.last_error.is_none() {
                self.last_error = Some("merged row failed counted-group replay".into());
            }
            return None;
        }
        if columns.is_empty() {
            self.counters.target_only_rows += 1;
        }
        self.counters.merged_relations += 1;
        self.counters.merge_depth_total += full.merge_depth as u64;
        self.counters.merge_sources_total += full.source_coeffs.len() as u64;
        Some(Relation::Row {
            columns: columns.into_iter().collect(),
            a,
            b,
        })
    }

    /// The large-prime form of the three-summand probe loop.
    #[allow(clippy::too_many_arguments)]
    pub fn decompose(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'_>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
        a: u64,
        b: u64,
    ) -> Vec<Relation> {
        let g = ctx.group;
        let mut out = Vec::new();
        for i in 0..self.small_points {
            let t = g.add(ops, point, g.neg(fb.points[i]));
            if g.is_identity(&t) {
                continue;
            }
            // The classic probe, unchanged.
            if let Some((j, k)) = self.table_small.probe(g, fb, ops, counters, t) {
                self.counters.direct_three_sum += 1;
                out.push(Relation::Summands(vec![i, j, k]));
                return out;
            }
            // The large table: F × D and, with two large primes, D × D.
            if let Some((j, k)) = self.table_large.probe(g, &self.envelope, ops, counters, t) {
                let (j_large, k_large) = (j >= self.small_points, k >= self.small_points);
                let (small, residual) = match (j_large, k_large) {
                    (false, false) => {
                        unreachable!("F × F pairs are filtered out of the large table")
                    }
                    (true, true) => {
                        self.counters.partials_two_large += 1;
                        (
                            vec![i],
                            vec![self.envelope.points[j], self.envelope.points[k]],
                        )
                    }
                    (false, true) => {
                        self.counters.partials_mixed += 1;
                        (vec![i, j], vec![self.envelope.points[k]])
                    }
                    (true, false) => {
                        self.counters.partials_mixed += 1;
                        (vec![i, k], vec![self.envelope.points[j]])
                    }
                };
                if let Some(row) = self.feed(g, fb, ops, a, b, small, residual) {
                    out.push(row);
                    return out;
                }
            }
            // Membership of T itself: in D, a one-large-prime partial; in
            // F, a two-summand relation the classic m = 3 loop never looks
            // for.  One canonicalisation and one lookup.
            counters.canonicalisations += 1;
            let key = self.table_small.orbit_key(t.x);
            counters.lookups += 1;
            match self.membership.get(&key) {
                Some(Membership::Domain) => {
                    self.counters.partials_one_large += 1;
                    if let Some(row) = self.feed(g, fb, ops, a, b, vec![i], vec![t]) {
                        out.push(row);
                        return out;
                    }
                }
                Some(Membership::Base) => {
                    counters.lookups += 1;
                    if let Some(&j) = fb.point_index.get(&t.pack()) {
                        self.counters.direct_two_sum += 1;
                        out.push(Relation::Summands(vec![i, j]));
                        return out;
                    }
                }
                None => {}
            }
        }
        out
    }
}

/// Keys of the orbits a base occupies, for tests and diagnostics.
pub fn orbit_key_set(table: &FrobeniusPairTable, fb: &FactorBase<FastPoint>) -> HashSet<u64> {
    fb.points.iter().map(|p| table.orbit_key(p.x)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::{
        collect_and_solve_relations, collect_and_solve_with_completion, koblitz_instance,
        Calibration, CompletionRule, IncrementalGauss, RestartPool, TargetSource,
    };
    use crate::cryptanalysis::ic_framework::plugins::{CompactOrbitScanBase, FrobeniusMitmOracle};
    use crate::cryptanalysis::ic_framework::stages::{
        DecompositionOracle, FactorBaseBuilder, Targets,
    };
    use crate::cryptanalysis::ic_framework::{run_pipeline, PipelineSpec};

    fn ctx<'a>(
        inst: &'a BinaryInstance,
        g: &'a BinaryGroup<'a>,
        q: FastPoint,
    ) -> InstanceCtx<'a, BinaryGroup<'a>> {
        InstanceCtx {
            group: g,
            generator: inst.generator,
            target: q,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(inst.n),
        }
    }

    fn multiple(inst: &BinaryInstance, k: u64) -> FastPoint {
        let g = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        g.mul(&mut ops, inst.generator, k % inst.r)
    }

    fn scan_base(
        ctx: &InstanceCtx<BinaryGroup<'_>>,
        inst: &BinaryInstance,
        columns: u32,
    ) -> FactorBase<FastPoint> {
        let (_, params) = Params::parse_spec(&format!(
            "compact-orbit-scan:columns={columns},raw_x_cap=100000000"
        ))
        .unwrap();
        let mut ops = GroupOps::default();
        CompactOrbitScanBase { instance: inst }
            .build(ctx, &params, &mut ops)
            .unwrap()
    }

    /// A base of known logarithms: the signed Frobenius orbits of `[k]G`
    /// for the listed `k`, one column each, unit `[k]G` (so column `c`'s
    /// unknown is `k_c`).
    fn known_log_base(inst: &BinaryInstance, seeds: &[u64]) -> FactorBase<FastPoint> {
        let g = BinaryGroup(&inst.fast);
        let r = inst.r;
        let lambda = inst.koblitz.as_ref().unwrap().lambda.to_u64().unwrap();
        let (mut points, mut col_of, mut coef_of) = (Vec::new(), Vec::new(), Vec::new());
        for (column, &k) in seeds.iter().enumerate() {
            let seed = multiple(inst, k);
            let mut member = seed;
            let mut coefficient = 1u64;
            for _ in 0..inst.n {
                points.push(member);
                col_of.push(column);
                coef_of.push(coefficient);
                points.push(g.neg(member));
                col_of.push(column);
                coef_of.push((r - coefficient) % r);
                member = inst.fast.frobenius_k(member, 1);
                coefficient = ((coefficient as u128 * lambda as u128) % r as u128) as u64;
            }
            assert_eq!(member, seed);
        }
        let fast = &inst.fast;
        FactorBase::from_column_map(
            format!("orbits of [k]G for k in {seeds:?}"),
            points,
            col_of,
            coef_of,
            seeds.len(),
            |p| p.pack(),
            |p| fast.neg(*p).pack(),
            |p| p.x,
        )
        .unwrap()
    }

    /// Planted partials over a base of known logarithms (columns `G` and
    /// `[3]G`), targets `[k]G` with `b = 0`:
    ///
    /// * one large prime: `[8]G = G + L` and `[3 + 7λ]G = [3]G + σ(L)`
    ///   with `L = [7]G` share the residual orbit and merge to a row
    ///   relating the two columns, which the counted group confirms;
    /// * two large primes: `[20]G = G + L₁ + L₂`, `[6]G = G + L₁`,
    ///   `[3 + 14λ]G = [3]G + σ(L₂)` with `L₁ = [5]G`, `L₂ = [14]G` form
    ///   a cycle that closes at depth two;
    /// * `[9]G = G + L` with `L = [7]G` is false and is rejected without
    ///   storing anything.
    #[test]
    fn planted_one_and_two_large_prime_partials_merge_and_verify() {
        let inst = koblitz_instance(0, 13).expect("K_0 at n = 13");
        let g = BinaryGroup(&inst.fast);
        let planted = 777 % inst.r;
        let q = multiple(&inst, planted);
        let ctx = ctx(&inst, &g, q);
        let fb = known_log_base(&inst, &[1, 3]);
        let r = inst.r;
        let lambda = inst.koblitz.as_ref().unwrap().lambda.to_u64().unwrap();
        let lam = |k: u64| ((k as u128 * lambda as u128) % r as u128) as u64;
        let config = LargePrimeConfig {
            large_primes: 2,
            domain_factor: 3.0,
            raw_x_cap: 1 << 13,
            limits: LargePrimeLimits::default(),
        };
        let mut ops = GroupOps::default();
        let mut lp = LargePrimeCollector::new(&inst, &ctx, &fb, config, &mut ops).unwrap();
        assert_eq!(lp.domain_orbits(), 6, "2 columns × factor 3");
        assert_eq!(lp.domain_points(), 6 * 2 * 13);
        for k in [5u64, 7, 14] {
            assert!(
                !fb.point_index.contains_key(&multiple(&inst, k).pack()),
                "[{k}]G must not lie in the base"
            );
        }
        let i_g = fb.point_index[&inst.generator.pack()];
        let i_3g = fb.point_index[&multiple(&inst, 3).pack()];

        // One large prime.
        let l = multiple(&inst, 7);
        let sigma_l = inst.fast.frobenius_k(l, 1);
        assert!(lp
            .feed(&g, &fb, &mut ops, 8, 0, vec![i_g], vec![l])
            .is_none());
        let merged = lp
            .feed(
                &g,
                &fb,
                &mut ops,
                (3 + lam(7)) % r,
                0,
                vec![i_3g],
                vec![sigma_l],
            )
            .expect("two partials on one residual orbit close");
        let Relation::Row { columns, a, b } = &merged else {
            panic!("expected a merged row, got {merged:?}");
        };
        assert_eq!(*b, 0);
        // The merged row is  (3 + 7λ) − λ·8 = x₁ − λ·x₀  with x₀ = 1, x₁ = 3.
        let cols: BTreeMap<usize, u64> = columns.iter().copied().collect();
        assert_eq!(cols.get(&1), Some(&1));
        assert_eq!(cols.get(&0), Some(&((r - lambda) % r)));
        assert_eq!(*a, (3 + lam(7) + r - lam(8)) % r);
        assert_eq!(lp.counters.merged_relations, 1);
        assert_eq!(lp.counters.verification_failures, 0);

        // Two large primes.
        let l1 = multiple(&inst, 5);
        let l2 = multiple(&inst, 14);
        let sigma_l2 = inst.fast.frobenius_k(l2, 1);
        assert!(lp
            .feed(&g, &fb, &mut ops, 20, 0, vec![i_g], vec![l1, l2])
            .is_none());
        assert!(lp
            .feed(&g, &fb, &mut ops, 6, 0, vec![i_g], vec![l1])
            .is_none());
        let third = lp
            .feed(
                &g,
                &fb,
                &mut ops,
                (3 + lam(14)) % r,
                0,
                vec![i_3g],
                vec![sigma_l2],
            )
            .expect("the two-large-prime cycle closes");
        let Relation::Row { columns, a, b } = &third else {
            panic!("expected a merged row, got {third:?}");
        };
        assert_eq!(*b, 0);
        // (3 + 14λ) − λ·(20 − 6) = x₁ − λ·(x₀ − x₀) = x₁  ⇒ a = 3, row = {1: 1}.
        let cols: BTreeMap<usize, u64> = columns.iter().copied().collect();
        assert_eq!(cols, BTreeMap::from([(1usize, 1u64)]));
        assert_eq!(*a, 3);
        assert_eq!(lp.counters.merged_relations, 2);
        // Two sources closed the first row, three the cycle; each closing
        // row was reduced against one stored pivot (the eliminator's
        // `merge_depth` counts the reductions of the incoming row only).
        assert_eq!(lp.counters.merge_sources_total, 2 + 3);
        assert_eq!(lp.counters.merge_depth_total, 2);
        assert_eq!(lp.counters.verification_failures, 0);
        assert_eq!(lp.counters.rejected, 0);

        // A wrong coefficient is rejected before anything is stored.
        let stored_before = lp.eliminator.stats().pivot_rows;
        assert!(lp
            .feed(&g, &fb, &mut ops, 9, 0, vec![i_g], vec![l])
            .is_none());
        assert_eq!(lp.counters.rejected, 1);
        assert_eq!(lp.eliminator.stats().pivot_rows, stored_before);
        assert!(lp.last_error().is_some());
        assert!(
            ops.scalar_mults > 0,
            "merged rows were replayed in the counted group"
        );
    }

    /// With `large_primes` unset the framework loop runs the classic
    /// oracle through `decompose_relations`, and its counters equal the
    /// classic `collect_and_solve_with_completion` run of the same
    /// oracle on the same seed, counter for counter.
    #[test]
    fn large_primes_zero_reproduces_the_classic_path_exactly() {
        let inst = koblitz_instance(0, 19).expect("K_0 at n = 19");
        let g = BinaryGroup(&inst.fast);
        let planted = 4242 % inst.r;
        let q = multiple(&inst, planted);
        let ctx = ctx(&inst, &g, q);
        let fb = scan_base(&ctx, &inst, 4);
        let (_, params) = Params::parse_spec("mitm-frobenius:m=3").unwrap();
        let run = |seam: bool| {
            let mut oracle = FrobeniusMitmOracle::new(3, &inst);
            let mut ops = GroupOps::default();
            oracle.prepare(&ctx, &fb, &params, &mut ops).unwrap();
            let mut matrix = IncrementalGauss::new(fb.columns + 1, inst.r);
            if seam {
                collect_and_solve_relations(
                    &g,
                    inst.generator,
                    q,
                    inst.r,
                    inst.cofactor,
                    &fb,
                    99,
                    200_000,
                    TargetSource::Walk,
                    RestartPool::Lazy,
                    &mut matrix,
                    CompletionRule::TargetPinned,
                    |ops, ctr, point, a, b| {
                        oracle.decompose_relations(&ctx, &fb, ops, ctr, point, a, b)
                    },
                )
            } else {
                collect_and_solve_with_completion(
                    &g,
                    inst.generator,
                    q,
                    inst.r,
                    inst.cofactor,
                    &fb,
                    99,
                    200_000,
                    TargetSource::Walk,
                    RestartPool::Lazy,
                    &mut matrix,
                    CompletionRule::TargetPinned,
                    |ops, ctr, point| oracle.decompose(&ctx, &fb, ops, ctr, point),
                )
            }
        };
        let classic = run(false);
        let seam = run(true);
        assert!(classic.verified && classic.recovered == Some(planted));
        assert_eq!(classic.trials, seam.trials);
        assert_eq!(classic.relations_found, seam.relations_found);
        assert_eq!(classic.rank, seam.rank);
        assert_eq!(classic.recovered, seam.recovered);
        assert_eq!(
            classic.relations.group_ops.adds,
            seam.relations.group_ops.adds
        );
        assert_eq!(
            classic.relations.group_ops.scalar_mults,
            seam.relations.group_ops.scalar_mults
        );
        assert_eq!(classic.relations.native, seam.relations.native);
        assert_eq!(classic.linear_algebra.native, seam.linear_algebra.native);
        assert_eq!(classic.counters.lookups, seam.counters.lookups);
        assert_eq!(
            classic.counters.canonicalisations,
            seam.counters.canonicalisations
        );
    }

    /// End to end through `run_pipeline`: one and two large primes recover
    /// and verify the planted logarithm, every merged row passed both
    /// checks, and no partial the collector formed was rejected.
    #[test]
    fn end_to_end_large_prime_runs_recover_the_planted_logarithm() {
        let inst = koblitz_instance(0, 19).expect("K_0 at n = 19");
        let g = BinaryGroup(&inst.fast);
        for large_primes in [1u64, 2] {
            let planted = 31_337 % inst.r;
            let q = multiple(&inst, planted);
            let ctx = ctx(&inst, &g, q);
            let spec = PipelineSpec {
                factor_base: "compact-orbit-scan".into(),
                factor_base_params: Params::parse_spec(
                    "compact-orbit-scan:columns=4,raw_x_cap=100000000",
                )
                .unwrap()
                .1,
                oracle: "mitm-frobenius".into(),
                oracle_params: Params::parse_spec(&format!(
                    "mitm-frobenius:m=3,large_primes={large_primes},lp_domain_factor=3"
                ))
                .unwrap()
                .1,
                targets: Targets::Walk,
                max_trials: 400_000,
                seed: 5,
                ..Default::default()
            };
            let base = CompactOrbitScanBase { instance: &inst };
            let mut oracle = FrobeniusMitmOracle::new(3, &inst);
            let calib = Calibration::default();
            let report =
                run_pipeline(&ctx, &spec, &base, &mut oracle, planted, &calib, None).unwrap();
            assert!(report.verified, "large_primes={large_primes}: {report:?}");
            let rel = &report.decomposition.cost;
            assert_eq!(rel.get("lp_verification_failures"), 0);
            assert_eq!(
                rel.get("lp_rejected"),
                0,
                "{:?}",
                oracle.large_prime_error()
            );
            assert!(rel.get("lp_partials_fed") > 0);
            if large_primes == 1 {
                assert_eq!(rel.get("lp_partials_two_large"), 0);
            }
            assert!(
                report
                    .decomposition
                    .setup
                    .native
                    .get("lp_table_entries")
                    .copied()
                    .unwrap_or(0)
                    > 0
            );
        }
    }
}
