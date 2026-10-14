//! `rootfind-mitm` and `rootfind-pairs`: the root-finding decomposition
//! oracles of `research/pdp_root_finder_20261012` as `ic_framework`
//! plug-ins, so that `ic.pipeline` can run them on any counted group
//! and charge them in the harness's unit.
//!
//! Both decide a target `R` by the group law alone, never by an
//! algebraic system.  For a base with abscissae `x_1 < … < x_A` (two
//! signed points above each, one above `x = 0`):
//!
//! - **`rootfind-mitm`** (`m = 3`): a pair table of `P_i + P_j`
//!   (negation folded, the harness's [`PairTable`]) built once per run;
//!   per target the `2A` points `R ∓ P_k`, one representative `P_k` an
//!   abscissa, computed in blocks of 32 abscissae and probed in the table,
//!   stopping at the first hit.  `m = 2` probes `R` itself.
//! - **`rootfind-pairs`** (`m = 3`): no table; the row points `U = R ∓ P_i`,
//!   and for each row the candidates `U ∓ P_j` over the abscissae
//!   `x_j ≥ x_i`, tested for membership in the base.  The enumeration
//!   floor, kept as the memoryless control.
//!
//! **Charged:** every point addition, one each, through the ledger: the
//! table build, the additions a target of `rootfind-mitm` (`2A` for a
//! refutation, the count `mitm:m=3` charges over the `2A` signed points;
//! at most a block more than `mitm` on a hit), the probe's confirming
//! addition on a key hit, and every candidate of `rootfind-pairs`.
//! **Counted, not charged:** table lookups and membership tests
//! (`lookups`).  Where the group batches ([`PairedAdd`]: the binary
//! group, through Montgomery's trick over a block) the additions are
//! still charged one each: the shared inversion is how they are done,
//! not work beside them, so the recorded cost is complete.  The witness
//! indices returned are exact: a probe re-adds the pair it found and a
//! membership test compares full point keys.
//!
//! The binary group batches through [`FastCurve::add_many`]; the prime
//! group uses the default one-inversion-an-addition path, so on a prime
//! curve the two ids measure the same additions as `mitm` and `subtract`
//! with the memoryless three-summand arm added.
//!
//! [`FastCurve::add_many`]: crate::cryptanalysis::koblitz_fast::FastCurve::add_many

use std::collections::BTreeMap;

use super::stages::{DecompositionOracle, InstanceCtx, Params};
use crate::cryptanalysis::ic_boundary::{
    BinaryGroup, CountedGroup, FactorBase, GroupOps, OracleCounters, PairTable, PrimeCurve,
};
use crate::cryptanalysis::koblitz_fast::BatchScratch;

/// `R + Q_j` for every `j`, appended to `out`, charged one addition
/// each.  A group that can share one field inversion over the slice
/// overrides the default and counts that inversion.
pub trait PairedAdd: CountedGroup {
    fn add_many(
        &self,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        r: Self::Elt,
        qs: &[Self::Elt],
        out: &mut Vec<Self::Elt>,
    ) {
        for &q in qs {
            out.push(self.add(ops, r, q));
        }
        let _ = counters;
    }
}

impl PairedAdd for PrimeCurve {}

impl PairedAdd for BinaryGroup<'_> {
    fn add_many(
        &self,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        r: Self::Elt,
        qs: &[Self::Elt],
        out: &mut Vec<Self::Elt>,
    ) {
        if qs.is_empty() {
            return;
        }
        ops.adds += qs.len() as u64;
        let _ = counters;
        let mut scratch = BatchScratch::default();
        self.0.add_many(r, qs, out, &mut scratch);
    }
}

/// Abscissae a batch of `rootfind-mitm`: `R ∓ P_k` for a block are
/// computed together, then probed, so a hit stops the scan within its
/// block and a refutation costs the same `2A` additions however it is
/// batched.
const BLOCK: usize = 32;

/// One signed representative an abscissa: the base point's index, the
/// point, its negative, and whether the point is its own negative (the
/// 2-torsion point above `x = 0`, which gets one addition, not two).
struct Reps<E> {
    index: Vec<usize>,
    point: Vec<E>,
    negated: Vec<E>,
    self_negated: Vec<bool>,
}

impl<E: Copy> Reps<E> {
    fn of<G: CountedGroup<Elt = E>>(g: &G, fb: &FactorBase<E>) -> Self {
        let mut reps = Reps {
            index: Vec::with_capacity(fb.abscissa_list.len()),
            point: Vec::with_capacity(fb.abscissa_list.len()),
            negated: Vec::with_capacity(fb.abscissa_list.len()),
            self_negated: Vec::with_capacity(fb.abscissa_list.len()),
        };
        for x in &fb.abscissa_list {
            let above = fb.x_index.get(x).expect("every listed abscissa is indexed");
            let i = above[0];
            reps.index.push(i);
            reps.point.push(fb.points[i]);
            reps.negated.push(g.neg(fb.points[i]));
            reps.self_negated.push(fb.neg_index[i] == i);
        }
        reps
    }
}

/// `R ∓ P_k` for every representative, as `(signed index of P_k, R − P_k)`
/// pairs: `R − P_k` carries `k`'s own index, `R + P_k` the index of `−P_k`,
/// so that `R = P_{index} + (the returned point)` in both cases.
fn row_points<G: PairedAdd>(
    g: &G,
    reps: &Reps<G::Elt>,
    fb: &FactorBase<G::Elt>,
    ops: &mut GroupOps,
    counters: &mut OracleCounters,
    r: G::Elt,
    from: usize,
    to: usize,
) -> Vec<(usize, G::Elt)> {
    let mut minus = Vec::with_capacity(to - from);
    g.add_many(ops, counters, r, &reps.negated[from..to], &mut minus);
    let plus_src: Vec<G::Elt> = (from..to)
        .filter(|&k| !reps.self_negated[k])
        .map(|k| reps.point[k])
        .collect();
    let mut plus = Vec::with_capacity(plus_src.len());
    g.add_many(ops, counters, r, &plus_src, &mut plus);
    let mut out = Vec::with_capacity(minus.len() + plus.len());
    for (k, s) in (from..).zip(minus) {
        out.push((reps.index[k], s));
    }
    let mut it = plus.into_iter();
    for k in from..to {
        if !reps.self_negated[k] {
            out.push((
                fb.neg_index[reps.index[k]],
                it.next().expect("one sum a point"),
            ));
        }
    }
    out
}

/// The `rootfind-mitm` oracle.
pub struct RootfindMitmOracle<E> {
    m: u32,
    table: Option<PairTable>,
    reps: Option<Reps<E>>,
}

impl<E> RootfindMitmOracle<E> {
    pub fn new(m: u32) -> Self {
        Self {
            m,
            table: None,
            reps: None,
        }
    }
}

impl<G: PairedAdd> DecompositionOracle<G> for RootfindMitmOracle<G::Elt> {
    fn name(&self) -> &str {
        "rootfind-mitm"
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "root-finding by the group law over a pair table ({}): the 2A points R ∓ P_k, one \
             representative an abscissa, computed in one batch and probed; m = {}",
            if params.flag("full_table") {
                "built in full"
            } else {
                "negation folded"
            },
            self.m
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("m", "summands a relation: 2 (probe R) or 3 (probe R ∓ P_k)"),
            (
                "full_table",
                "1 to build every i ≤ j pair instead of one pair per orbit under negation",
            ),
        ]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        if !matches!(self.m, 2 | 3) {
            return Err(format!("rootfind-mitm runs at m = 2 or 3, not {}", self.m));
        }
        let table = if params.flag("full_table") {
            PairTable::build(ctx.group, fb)
        } else {
            PairTable::build_negation_folded(ctx.group, fb)
        };
        ops.merge(table.build_ops);
        self.table = Some(table);
        self.reps = Some(Reps::of(ctx.group, fb));
        Ok(())
    }

    fn setup_native(&self) -> BTreeMap<String, u64> {
        let mut counts = BTreeMap::new();
        if let Some(table) = &self.table {
            counts.insert("pair_table_entries".into(), table.entries);
        }
        if let Some(reps) = &self.reps {
            counts.insert("abscissae".into(), reps.index.len() as u64);
        }
        counts
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
    ) -> Option<Vec<usize>> {
        let g = ctx.group;
        let table = self.table.as_ref()?;
        let reps = self.reps.as_ref()?;
        if self.m == 2 {
            return table
                .probe(g, fb, ops, &mut counters.lookups, point)
                .map(|(j, l)| vec![j, l]);
        }
        let a = reps.index.len();
        let mut from = 0;
        while from < a {
            let to = (from + BLOCK).min(a);
            for (k, s) in row_points(g, reps, fb, ops, counters, point, from, to) {
                if g.is_identity(&s) {
                    continue;
                }
                if let Some((j, l)) = table.probe(g, fb, ops, &mut counters.lookups, s) {
                    return Some(vec![k, j, l]);
                }
            }
            from = to;
        }
        None
    }
}

/// The `rootfind-pairs` oracle: the memoryless three-summand search.
pub struct RootfindPairsOracle<E> {
    m: u32,
    reps: Option<Reps<E>>,
}

impl<E> RootfindPairsOracle<E> {
    pub fn new(m: u32) -> Self {
        Self { m, reps: None }
    }
}

impl<G: PairedAdd> DecompositionOracle<G> for RootfindPairsOracle<G::Elt> {
    fn name(&self) -> &str {
        "rootfind-pairs"
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, _params: &Params) -> String {
        format!(
            "root-finding by the group law with no table: row points U = R ∓ P_i, candidates \
             U ∓ P_j over x_j ≥ x_i tested in the base; the enumeration floor; m = {}",
            self.m
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[("m", "summands a relation: 3")]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        _params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<(), String> {
        if self.m != 3 {
            return Err(format!("rootfind-pairs runs at m = 3, not {}", self.m));
        }
        self.reps = Some(Reps::of(ctx.group, fb));
        Ok(())
    }

    fn setup_native(&self) -> BTreeMap<String, u64> {
        let mut counts = BTreeMap::new();
        if let Some(reps) = &self.reps {
            counts.insert("abscissae".into(), reps.index.len() as u64);
        }
        counts
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
    ) -> Option<Vec<usize>> {
        let g = ctx.group;
        let reps = self.reps.as_ref()?;
        let a = reps.index.len();
        let rows = row_points(g, reps, fb, ops, counters, point, 0, a);
        // Row `i` pairs with the abscissae `x_j ≥ x_i`; the position of a
        // row point's abscissa among the representatives is recovered
        // from its index.
        let position: BTreeMap<usize, usize> = reps
            .index
            .iter()
            .enumerate()
            .flat_map(|(pos, &i)| [(i, pos), (fb.neg_index[i], pos)])
            .collect();
        for (si, u) in rows {
            if g.is_identity(&u) {
                continue;
            }
            let from = position[&si];
            for (sj, c) in row_points(g, reps, fb, ops, counters, u, from, a) {
                counters.lookups += 1;
                if g.is_identity(&c) {
                    continue;
                }
                if let Some(l) = fb.index_of_key(g.key(&c)) {
                    return Some(vec![si, sj, l]);
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::super::plugins::{BinarySubspaceBase, MitmOracle, PrimeAbscissaBase};
    use super::super::stages::FactorBaseBuilder;
    use super::*;
    use crate::cryptanalysis::ic_boundary::{koblitz_instance, prime_instance_for};

    fn check<G: PairedAdd>(ctx: &InstanceCtx<G>, fb: &FactorBase<G::Elt>, targets: &[G::Elt]) {
        let g = ctx.group;
        let params = Params::default();
        let mut ops = GroupOps::default();
        let mut mitm = MitmOracle::new(3);
        let mut rf_mitm = RootfindMitmOracle::new(3);
        let mut rf_pairs = RootfindPairsOracle::new(3);
        let (_, folded) = Params::parse_spec("mitm:negation_folded=1").unwrap();
        mitm.prepare(ctx, fb, &folded, &mut ops).unwrap();
        rf_mitm.prepare(ctx, fb, &params, &mut ops).unwrap();
        rf_pairs.prepare(ctx, fb, &params, &mut ops).unwrap();
        let (mut found, mut refuted) = (0, 0);
        for &t in targets {
            let (mut o1, mut o2, mut o3) = (
                GroupOps::default(),
                GroupOps::default(),
                GroupOps::default(),
            );
            let (mut c1, mut c2, mut c3) = (
                OracleCounters::default(),
                OracleCounters::default(),
                OracleCounters::default(),
            );
            let a = mitm.decompose(ctx, fb, &mut o1, &mut c1, t);
            let b = rf_mitm.decompose(ctx, fb, &mut o2, &mut c2, t);
            let c = rf_pairs.decompose(ctx, fb, &mut o3, &mut c3, t);
            assert_eq!(a.is_some(), b.is_some(), "verdicts differ on a target");
            assert_eq!(a.is_some(), c.is_some(), "the memoryless verdict differs");
            for w in [a.as_ref(), b.as_ref(), c.as_ref()].into_iter().flatten() {
                let mut scratch = GroupOps::default();
                let mut sum = g.identity();
                for &i in w {
                    sum = g.add(&mut scratch, sum, fb.points[i]);
                }
                assert_eq!(sum, t, "a witness does not sum to its target");
            }
            if a.is_none() {
                // A refutation visits every signed point in both table
                // oracles: the same additions, however they are batched.
                assert_eq!(o1.adds, o2.adds, "refutation charge differs");
                refuted += 1;
            } else {
                found += 1;
            }
        }
        assert!(found > 0 && refuted > 0, "found {found}, refuted {refuted}");
    }

    #[test]
    fn binary_rootfind_oracles_agree_with_mitm_on_a_koblitz_subspace_base() {
        let inst = koblitz_instance(1, 17).expect("registered n17 Koblitz instance");
        let g = BinaryGroup(&inst.fast);
        let ctx = InstanceCtx {
            group: &g,
            generator: inst.generator,
            target: inst.generator,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(inst.n),
        };
        let (_, params) = Params::parse_spec("binary-subspace:dimension=6").unwrap();
        let mut ops = GroupOps::default();
        let fb = BinarySubspaceBase { instance: &inst }
            .build(&ctx, &params, &mut ops)
            .unwrap();
        let mut scratch = GroupOps::default();
        let targets: Vec<_> = (1..=96u64)
            .map(|k| {
                g.mul(
                    &mut scratch,
                    inst.generator,
                    k.wrapping_mul(0x9E37_79B9) % inst.r,
                )
            })
            .collect();
        check(&ctx, &fb, &targets);
    }

    #[test]
    fn prime_rootfind_oracles_agree_with_mitm_on_an_abscissa_base() {
        let inst = prime_instance_for(18, 7);
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: inst.generator_point(),
            target: inst.generator_point(),
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        };
        let (_, params) = Params::parse_spec("prime-abscissa:size=24").unwrap();
        let mut ops = GroupOps::default();
        let fb = PrimeAbscissaBase { instance: &inst }
            .build(&ctx, &params, &mut ops)
            .unwrap();
        let mut scratch = GroupOps::default();
        let targets: Vec<_> = (1..=96u64)
            .map(|k| {
                inst.curve.mul(
                    &mut scratch,
                    inst.generator_point(),
                    k.wrapping_mul(0x9E37_79B9) % inst.r,
                )
            })
            .collect();
        check(&ctx, &fb, &targets);
    }
}
