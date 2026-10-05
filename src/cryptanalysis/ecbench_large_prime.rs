//! Exact binary-curve import and a bounded `m`-summand index-calculus
//! collector with the single- and double-large-prime variations.
//!
//! This is deliberately a correctness engine, not a scaling claim.  It
//! enumerates the meet-in-the-middle state space exactly and refuses before
//! allocation when `max_states` is too small.  In particular, `m = n - 1`
//! means exactly that; it is never replaced by a cheaper decomposition.

use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::Path;
use std::time::Instant;

use num_bigint::BigUint;
use serde::{Deserialize, Serialize};

use crate::binary_ecc::curve::scalar_mul;
use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crate::cryptanalysis::gaudry_cubic::{LargePrimeEliminator, SparseRel};
use crate::cryptanalysis::ic_boundary::{
    ArtinSchreier, BinaryGroup, BinaryInstance, CountedGroup, GroupOps, IncrementalGauss, RowStatus,
};
use crate::cryptanalysis::ic_measurement as measurement;
use crate::cryptanalysis::koblitz_fast::{FastCurve, FastPoint};
use crate::cryptanalysis::koblitz_index_calculus::{
    frobenius_eigenvalue, is_irreducible_f2, koblitz_point_count, KoblitzCurve,
};
use crate::cryptanalysis::residual_walk::is_prime_u64;
use crate::cryptanalysis::semaev_decomp::Gf2;

pub const MANIFEST_SCHEMA: &str = "ecbench.curve-corpus/v1";

#[derive(Clone, Debug, Deserialize)]
pub struct CorpusManifest {
    #[serde(default)]
    pub schema: Option<String>,
    pub families: Vec<CorpusFamily>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct CorpusFamily {
    pub id: String,
    pub field: CorpusField,
    pub curve: CorpusCurve,
    pub factor_base: CorpusFactorBase,
    pub dlog_challenge: CorpusChallenge,
}

#[derive(Clone, Debug, Deserialize)]
pub struct CorpusField {
    pub characteristic: u32,
    pub degree: u32,
    pub modulus: u64,
}

#[derive(Clone, Debug, Deserialize)]
pub struct CorpusCurve {
    pub equation: String,
    pub order: u64,
}

#[derive(Clone, Debug, Deserialize)]
pub struct CorpusFactorBase {
    pub dimension: u32,
    pub small_dimension: u32,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct CorpusPoint {
    pub x: u64,
    pub y: u64,
}

#[derive(Clone, Debug, Deserialize)]
pub struct CorpusChallenge {
    pub base_point: CorpusPoint,
    pub target_point: CorpusPoint,
    pub group_order: u64,
    pub subgroup_order: u64,
    pub cofactor: u64,
    pub planted_log: u64,
}

pub struct ImportedFamily {
    pub id: String,
    pub instance: BinaryInstance,
    pub target: FastPoint,
    pub small_dimension: u32,
    pub envelope_dimension: u32,
    /// Used by the CLI only after the solver returns; never passed to
    /// [`solve`].
    pub expected: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ValidationFacts {
    pub id: String,
    pub n: u32,
    pub modulus: u64,
    pub group_order: u64,
    pub subgroup_order: u64,
    pub cofactor: u64,
    pub small_dimension: u32,
    pub envelope_dimension: u32,
    pub base_point: CorpusPoint,
    pub target_point: CorpusPoint,
    pub checks: Vec<String>,
}

/// Parse the external corpus without accepting a silently different shape.
pub fn load_manifest(path: &Path) -> Result<CorpusManifest, String> {
    let text = fs::read_to_string(path)
        .map_err(|e| format!("read corpus manifest {}: {e}", path.display()))?;
    let manifest: CorpusManifest =
        serde_json::from_str(&text).map_err(|e| format!("parse corpus manifest: {e}"))?;
    if manifest.families.is_empty() {
        return Err("corpus manifest contains no families".into());
    }
    if let Some(schema) = &manifest.schema {
        if schema != MANIFEST_SCHEMA {
            return Err(format!("unsupported corpus schema `{schema}`"));
        }
    }
    Ok(manifest)
}

/// Construct an exact explicit binary instance and check every fact that can
/// be checked without re-counting all `2^n` abscissae.
pub fn binary_explicit_instance(
    n: u32,
    modulus: u64,
    a: u64,
    b: u64,
    group_order: u64,
    r: u64,
    gx: u64,
    gy: u64,
) -> Result<BinaryInstance, String> {
    if !(3..=FastCurve::MAX_DEGREE).contains(&n) {
        return Err(format!("binary_explicit takes 3 <= n <= 62, not {n}"));
    }
    let top = 1u64 << n;
    if modulus & top == 0 || modulus >> (n + 1) != 0 {
        return Err(format!("modulus 0x{modulus:x} does not have degree {n}"));
    }
    if !is_irreducible_f2(modulus) {
        return Err(format!("modulus 0x{modulus:x} is reducible over F_2"));
    }
    if a >= top || b == 0 || b >= top {
        return Err("binary_explicit requires a,b in the field and b != 0".into());
    }
    if r < 3 || !is_prime_u64(r) || !group_order.is_multiple_of(r) {
        return Err(format!(
            "subgroup order {r} must be prime and divide group order {group_order}"
        ));
    }
    let low_terms = (0..n).filter(|i| (modulus >> i) & 1 == 1).collect();
    let irreducible = IrreduciblePoly {
        degree: n,
        low_terms,
    };
    let gf = Gf2::new(&irreducible);
    let artin_schreier = ArtinSchreier::new(&gf);
    let shell = BinaryCurve {
        m: n,
        irreducible: irreducible.clone(),
        a: gf.to_element(a),
        b: gf.to_element(b),
        generator: BinaryPoint::Infinity,
        order: BigUint::from(r),
        cofactor: BigUint::from(group_order / r),
    };
    let fast = FastCurve::new(&shell).ok_or("explicit binary curve is too wide")?;
    let generator = FastPoint::affine(gx, gy);
    if !fast.is_on_curve(generator) {
        return Err("explicit binary generator is not on the curve".into());
    }
    if !fast.mul_u64(generator, r).infinity {
        return Err("[r]G is not the identity".into());
    }
    let general = BinaryCurve {
        generator: fast.lower(generator),
        ..shell
    };
    // A curve with a,b in F_2 and b=1 is genuinely Koblitz regardless of
    // which irreducible polynomial represents F_{2^n}. Certify its point
    // count and Frobenius eigenvalue so ecbench can expose the matched 2n
    // automorphisms while this module still uses its polynomial-basis base.
    let koblitz = if a <= 1 && b == 1 {
        let expected = koblitz_point_count(a as u8, n);
        if expected != BigUint::from(group_order) {
            return Err(format!(
                "Koblitz recurrence gives #E={expected}, not the supplied {group_order}"
            ));
        }
        let trace = if a == 0 { -1 } else { 1 };
        let order = BigUint::from(r);
        let lambda = frobenius_eigenvalue(&general, trace, &order)
            .ok_or("could not certify the Frobenius eigenvalue")?;
        Some(KoblitzCurve {
            a: a as u8,
            n,
            curve: general,
            trace,
            group_order: BigUint::from(group_order),
            subgroup_order: order,
            cofactor: BigUint::from(group_order / r),
            lambda,
            k: 1,
            q: 2,
            a_index: a,
            b_index: 1,
            subfield_basis: vec![F2mElement::one(n)],
        })
    } else {
        None
    };
    let mut instance = BinaryInstance {
        name: String::new(),
        n,
        irreducible,
        fast,
        gf,
        artin_schreier,
        a,
        b,
        group_order,
        r,
        cofactor: group_order / r,
        generator,
        koblitz,
    };
    instance.name = instance.curve_id().slug;
    Ok(instance)
}

pub fn import_family(family: &CorpusFamily) -> Result<(ImportedFamily, ValidationFacts), String> {
    if family.field.characteristic != 2 {
        return Err(format!("{} is not a binary-field family", family.id));
    }
    if family.curve.equation.replace(' ', "") != "y^2+x*y=x^3+x^2+1" {
        return Err(format!(
            "{} has unsupported equation `{}`",
            family.id, family.curve.equation
        ));
    }
    let c = &family.dlog_challenge;
    if family.curve.order != c.group_order {
        return Err(format!("{} disagrees on the curve group order", family.id));
    }
    if c.cofactor.checked_mul(c.subgroup_order) != Some(c.group_order) {
        return Err(format!("{} does not satisfy #E = h*r", family.id));
    }
    if family.factor_base.small_dimension > family.factor_base.dimension
        || family.factor_base.dimension > family.field.degree
    {
        return Err(format!(
            "{} has invalid nested factor-base dimensions",
            family.id
        ));
    }
    let instance = binary_explicit_instance(
        family.field.degree,
        family.field.modulus,
        1,
        1,
        c.group_order,
        c.subgroup_order,
        c.base_point.x,
        c.base_point.y,
    )?;
    let target = FastPoint::affine(c.target_point.x, c.target_point.y);
    if !instance.fast.is_on_curve(target) {
        return Err(format!("{} target is not on the curve", family.id));
    }
    if !instance.fast.mul_u64(target, instance.r).infinity {
        return Err(format!(
            "{} target is not in the order-r subgroup",
            family.id
        ));
    }
    if instance
        .fast
        .mul_u64(instance.generator, c.planted_log % instance.r)
        != target
    {
        return Err(format!(
            "{} planted public equality does not hold",
            family.id
        ));
    }
    independent_verify(&instance, target, c.planted_log % instance.r)?;
    let facts = ValidationFacts {
        id: family.id.clone(),
        n: family.field.degree,
        modulus: family.field.modulus,
        group_order: c.group_order,
        subgroup_order: c.subgroup_order,
        cofactor: c.cofactor,
        small_dimension: family.factor_base.small_dimension,
        envelope_dimension: family.factor_base.dimension,
        base_point: c.base_point,
        target_point: c.target_point,
        checks: vec![
            "irreducible_modulus".into(),
            "supported_curve_equation".into(),
            "group_order_equals_cofactor_times_prime_subgroup_order".into(),
            "base_and_target_on_curve".into(),
            "base_and_target_in_subgroup".into(),
            "planted_public_equality_fast_and_independent".into(),
        ],
    };
    Ok((
        ImportedFamily {
            id: family.id.clone(),
            instance,
            target,
            small_dimension: family.factor_base.small_dimension,
            envelope_dimension: family.factor_base.dimension,
            expected: c.planted_log % c.subgroup_order,
        },
        facts,
    ))
}

fn general_curve(instance: &BinaryInstance) -> BinaryCurve {
    BinaryCurve {
        m: instance.n,
        irreducible: instance.irreducible.clone(),
        a: instance.gf.to_element(instance.a),
        b: instance.gf.to_element(instance.b),
        generator: instance.fast.lower(instance.generator),
        order: BigUint::from(instance.r),
        cofactor: BigUint::from(instance.cofactor),
    }
}

pub fn independent_verify(
    instance: &BinaryInstance,
    target: FastPoint,
    scalar: u64,
) -> Result<(), String> {
    let curve = general_curve(instance);
    let got = scalar_mul(&curve, &curve.generator, &BigUint::from(scalar));
    if got != instance.fast.lower(target) {
        return Err("independent general-arithmetic verification failed".into());
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum RelationClass {
    Full,
    SingleLargePrime,
    DoubleLargePrime,
}

#[derive(Clone, Debug, Serialize)]
pub struct PartialRelation {
    pub trial: u64,
    pub a: u64,
    pub b: u64,
    pub class: RelationClass,
    pub raw_terms: Vec<usize>,
    pub columns: Vec<(usize, u64)>,
}

#[derive(Clone, Debug, Serialize)]
pub struct SolveConfig {
    pub small_dimension: u32,
    pub envelope_dimension: u32,
    pub summands: u32,
    pub max_large_primes: u8,
    pub max_trials: u64,
    pub max_states: u64,
    pub seed: u64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SolveCounters {
    pub trials: u64,
    pub decompositions: u64,
    pub rejected_large_prime_count: u64,
    pub relation_histogram: [u64; 4],
    pub accepted_relations: u64,
    pub eliminated_full_relations: u64,
    pub independent_rows: u64,
    pub matrix_rank: u64,
    pub mitm_lookups_uncharged: u64,
    pub combination_states_uncharged: u64,
    pub lp_merge_ops_uncharged: u64,
    pub row_ops_uncharged: u64,
    pub lp_pivots: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct FactorBaseReport {
    pub raw_points: usize,
    pub small_raw_points: usize,
    pub columns: usize,
    pub small_columns: usize,
    pub large_columns: usize,
    /// Canonical projected point keys in column order. The native ecbench
    /// adapter hashes these into its factor-base identity; reports omit them.
    #[serde(skip)]
    pub column_keys: Vec<u64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct MitmReport {
    pub left_terms: u32,
    pub right_terms: u32,
    pub left_states: u64,
    pub right_states: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct LargePrimeReport {
    pub recovered: Option<u64>,
    pub exhausted: bool,
    pub verified_fast: bool,
    pub verified_independent: bool,
    pub verification_performed: bool,
    pub config: SolveConfig,
    pub factor_base: FactorBaseReport,
    pub mitm: MitmReport,
    pub counters: SolveCounters,
    pub phase_group_ops: BTreeMap<String, GroupOps>,
    pub solve_wall_ns: u64,
    pub relation_samples: Vec<PartialRelation>,
}

#[derive(Clone)]
struct Atom {
    raw: FastPoint,
    column: Option<usize>,
    coefficient: u64,
    small: bool,
}

struct FactorBase {
    atoms: Vec<Atom>,
    small_columns: usize,
    report: FactorBaseReport,
}

fn canonical(curve: &FastCurve, p: FastPoint) -> (FastPoint, bool) {
    let neg = curve.neg(p);
    if (p.x, p.y) <= (neg.x, neg.y) {
        (p, false)
    } else {
        (neg, true)
    }
}

fn build_factor_base(
    instance: &BinaryInstance,
    small_dimension: u32,
    envelope_dimension: u32,
    ops: &mut GroupOps,
) -> Result<FactorBase, String> {
    if small_dimension > envelope_dimension || envelope_dimension > instance.n {
        return Err("factor-base dimensions must satisfy small <= envelope <= n".into());
    }
    if envelope_dimension >= 63 {
        return Err("factor-base envelope does not fit a word".into());
    }
    let small_x = 1u64 << small_dimension;
    let envelope_x = 1u64 << envelope_dimension;
    let group = BinaryGroup(&instance.fast);
    let mut raw = Vec::new();
    let mut keys: BTreeMap<u64, bool> = BTreeMap::new();
    for x in 0..envelope_x {
        for p in instance.points_with_x(x) {
            let projected = group.mul(ops, p, instance.cofactor);
            if !projected.infinity {
                let (c, _) = canonical(&instance.fast, projected);
                keys.entry(c.pack())
                    .and_modify(|s| *s |= x < small_x)
                    .or_insert(x < small_x);
            }
            raw.push((p, projected, x < small_x));
        }
    }
    let mut col_of = HashMap::new();
    let mut next = 0usize;
    for (&key, &small) in &keys {
        if small {
            col_of.insert(key, next);
            next += 1;
        }
    }
    let small_columns = next;
    for (&key, &small) in &keys {
        if !small {
            col_of.insert(key, next);
            next += 1;
        }
    }
    let atoms = raw
        .into_iter()
        .map(|(raw, projected, small)| {
            if projected.infinity {
                Atom {
                    raw,
                    column: None,
                    coefficient: 0,
                    small,
                }
            } else {
                let (c, negated) = canonical(&instance.fast, projected);
                Atom {
                    raw,
                    column: Some(col_of[&c.pack()]),
                    coefficient: if negated { instance.r - 1 } else { 1 },
                    small,
                }
            }
        })
        .collect::<Vec<_>>();
    let mut column_keys = vec![0u64; next];
    for (key, column) in &col_of {
        column_keys[*column] = *key;
    }
    let small_raw_points = atoms.iter().filter(|a| a.small).count();
    Ok(FactorBase {
        report: FactorBaseReport {
            raw_points: atoms.len(),
            small_raw_points,
            columns: next,
            small_columns,
            large_columns: next - small_columns,
            column_keys,
        },
        atoms,
        small_columns,
    })
}

#[derive(Clone)]
struct State {
    sum: FastPoint,
    terms: Vec<usize>,
}

struct Decomposer {
    left: HashMap<FastPoint, Vec<Vec<usize>>>,
    right: Vec<State>,
    report: MitmReport,
}

fn combinations(n: usize, k: usize) -> u128 {
    if k == 0 {
        return 1;
    }
    let mut v = 1u128;
    for i in 1..=k {
        v = v.saturating_mul((n + i - 1) as u128) / i as u128;
    }
    v
}

fn enumerate_states(
    group: &BinaryGroup<'_>,
    atoms: &[Atom],
    terms: usize,
    ops: &mut GroupOps,
) -> Vec<State> {
    fn visit(
        group: &BinaryGroup<'_>,
        atoms: &[Atom],
        left: usize,
        start: usize,
        acc: FastPoint,
        current: &mut Vec<usize>,
        out: &mut Vec<State>,
        ops: &mut GroupOps,
    ) {
        if left == 0 {
            out.push(State {
                sum: acc,
                terms: current.clone(),
            });
            return;
        }
        for i in start..atoms.len() {
            let sum = group.add(ops, acc, atoms[i].raw);
            current.push(i);
            visit(group, atoms, left - 1, i, sum, current, out, ops);
            current.pop();
        }
    }
    let cap = combinations(atoms.len(), terms).min(usize::MAX as u128) as usize;
    let mut out = Vec::with_capacity(cap);
    visit(
        group,
        atoms,
        terms,
        0,
        FastPoint::INFINITY,
        &mut Vec::with_capacity(terms),
        &mut out,
        ops,
    );
    out
}

impl Decomposer {
    fn new(
        group: &BinaryGroup<'_>,
        base: &FactorBase,
        summands: u32,
        max_states: u64,
        ops: &mut GroupOps,
    ) -> Result<Self, String> {
        if summands < 2 {
            return Err("summands must be at least 2".into());
        }
        let left_terms = summands.div_ceil(2);
        let right_terms = summands - left_terms;
        let lc = combinations(base.atoms.len(), left_terms as usize);
        let rc = combinations(base.atoms.len(), right_terms as usize);
        let total = lc.saturating_add(rc);
        if total > max_states as u128 {
            return Err(format!(
                "exact {summands}-summand MITM needs {lc} + {rc} = {total} states, exceeding max_states={max_states}"
            ));
        }
        let left_states = enumerate_states(group, &base.atoms, left_terms as usize, ops);
        let right = enumerate_states(group, &base.atoms, right_terms as usize, ops);
        let mut left: HashMap<FastPoint, Vec<Vec<usize>>> = HashMap::new();
        for state in left_states {
            left.entry(state.sum).or_default().push(state.terms);
        }
        Ok(Self {
            left,
            report: MitmReport {
                left_terms,
                right_terms,
                left_states: lc as u64,
                right_states: rc as u64,
            },
            right,
        })
    }

    fn find(
        &self,
        group: &BinaryGroup<'_>,
        base: &FactorBase,
        target: FastPoint,
        max_large_primes: u8,
        r: u64,
        ops: &mut GroupOps,
        lookups: &mut u64,
    ) -> Option<(Vec<usize>, Vec<(usize, u64)>, usize)> {
        for right in &self.right {
            let need = group.add(ops, target, group.neg(right.sum));
            *lookups += 1;
            let Some(lefts) = self.left.get(&need) else {
                continue;
            };
            for left in lefts {
                let mut terms = left.clone();
                terms.extend_from_slice(&right.terms);
                let mut coeffs: BTreeMap<usize, u64> = BTreeMap::new();
                for &i in &terms {
                    let atom = &base.atoms[i];
                    if let Some(c) = atom.column {
                        let v = coeffs.entry(c).or_insert(0);
                        *v = ((*v as u128 + atom.coefficient as u128) % r as u128) as u64;
                    }
                }
                let projected: Vec<(usize, u64)> =
                    coeffs.into_iter().filter(|(_, v)| *v != 0).collect();
                let lp = projected
                    .iter()
                    .filter(|(c, _)| *c >= base.small_columns)
                    .count();
                if lp <= max_large_primes as usize {
                    return Some((terms, projected, lp));
                }
            }
        }
        None
    }
}

struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}

fn mod_mul(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// Solve one target.  The planted logarithm is intentionally absent from the
/// signature.
pub fn solve(
    instance: &BinaryInstance,
    target: FastPoint,
    config: SolveConfig,
) -> Result<LargePrimeReport, String> {
    solve_with_verification(instance, target, config, true)
}

/// Return the algebraic candidate without checking it against the target.
/// Native ecbench uses this entry point and verifies in its runner process.
pub fn solve_candidate(
    instance: &BinaryInstance,
    target: FastPoint,
    config: SolveConfig,
) -> Result<LargePrimeReport, String> {
    solve_with_verification(instance, target, config, false)
}

fn solve_with_verification(
    instance: &BinaryInstance,
    target: FastPoint,
    config: SolveConfig,
    verify_candidate: bool,
) -> Result<LargePrimeReport, String> {
    if config.max_large_primes > 2 {
        return Err("max_large_primes must be 0, 1, or 2".into());
    }
    if !instance.fast.is_on_curve(target) || !instance.fast.mul_u64(target, instance.r).infinity {
        return Err("target is not a point in the configured prime-order subgroup".into());
    }
    let start = Instant::now();
    let group = BinaryGroup(&instance.fast);
    let mut phase_group_ops = BTreeMap::new();

    let mut fb_ops = GroupOps::default();
    let base = build_factor_base(
        instance,
        config.small_dimension,
        config.envelope_dimension,
        &mut fb_ops,
    )?;
    if base.small_columns == 0 {
        return Err("projected small factor base has no columns".into());
    }
    phase_group_ops.insert("factor_base".into(), fb_ops);

    let mut setup_ops = GroupOps::default();
    let decomposer = Decomposer::new(
        &group,
        &base,
        config.summands,
        config.max_states,
        &mut setup_ops,
    )?;
    phase_group_ops.insert("oracle_setup".into(), setup_ops);

    let d_col = base.small_columns;
    let mut matrix = IncrementalGauss::new(base.small_columns + 1, instance.r);
    let mut eliminator = LargePrimeEliminator::new(base.small_columns + 1, instance.r);
    let mut rng = SplitMix64(config.seed ^ 0x4a6f_7578_4c50_7631);
    let mut rel_ops = GroupOps::default();
    let mut verify_ops = GroupOps::default();
    let mut counters = SolveCounters {
        combination_states_uncharged: decomposer.report.left_states
            + decomposer.report.right_states,
        ..Default::default()
    };
    let mut samples = Vec::new();
    let mut recovered = None;
    measurement::begin_online(measurement::Phase::TargetQuery);
    for trial in 1..=config.max_trials {
        counters.trials = trial;
        let a = rng.next() % instance.r;
        let b = 1 + rng.next() % (instance.r - 1);
        let ap = group.mul(&mut rel_ops, instance.generator, a);
        let bq = group.mul(&mut rel_ops, target, b);
        let residual = group.add(&mut rel_ops, ap, bq);
        let decomposition = {
            let _pdp = measurement::scope(measurement::Phase::TargetPdp);
            decomposer.find(
                &group,
                &base,
                residual,
                config.max_large_primes,
                instance.r,
                &mut rel_ops,
                &mut counters.mitm_lookups_uncharged,
            )
        };
        let Some((terms, projected, lp_count)) = decomposition else {
            continue;
        };
        counters.decompositions += 1;
        {
            let _check = measurement::scope(measurement::Phase::TargetRelationCheck);
            let mut check = FastPoint::INFINITY;
            for &i in &terms {
                check = group.add(&mut verify_ops, check, base.atoms[i].raw);
            }
            if check != residual {
                return Err("MITM returned a non-relation".into());
            }
        }
        counters.relation_histogram[lp_count.min(3)] += 1;
        if lp_count > config.max_large_primes as usize {
            counters.rejected_large_prime_count += 1;
            continue;
        }
        counters.accepted_relations += 1;
        let class = match lp_count {
            0 => RelationClass::Full,
            1 => RelationClass::SingleLargePrime,
            _ => RelationClass::DoubleLargePrime,
        };
        let mut cols = Vec::with_capacity(projected.len() + 1);
        for &(c, v) in &projected {
            cols.push((if c < base.small_columns { c } else { c + 1 }, v));
        }
        let hb = mod_mul(instance.cofactor % instance.r, b, instance.r);
        if hb != 0 {
            cols.push((d_col, instance.r - hb));
        }
        cols.sort_unstable();
        let rhs = mod_mul(instance.cofactor % instance.r, a, instance.r);
        if samples.len() < 12 {
            samples.push(PartialRelation {
                trial,
                a,
                b,
                class,
                raw_terms: terms,
                columns: cols.clone(),
            });
        }
        let _descent = measurement::scope(measurement::Phase::TargetDescent);
        let sparse = SparseRel { cols, rhs };
        let full = if lp_count == 0 {
            Some(sparse)
        } else {
            let before = counters.lp_merge_ops_uncharged;
            let out = eliminator.feed(sparse, &mut counters.lp_merge_ops_uncharged);
            if out.is_some() {
                counters.eliminated_full_relations += 1;
            }
            debug_assert!(counters.lp_merge_ops_uncharged >= before);
            out
        };
        if let Some(full) = full {
            let mut row = vec![0u64; base.small_columns + 1];
            for (c, v) in full.cols {
                if c >= row.len() {
                    return Err("large-prime elimination left a large column".into());
                }
                row[c] = v;
            }
            match matrix.add_row(row, full.rhs) {
                RowStatus::Independent => counters.independent_rows += 1,
                RowStatus::Dependent => {}
                RowStatus::Inconsistent => return Err("relation system became inconsistent".into()),
            }
            if let Some(d) = matrix.pinned(d_col) {
                recovered = Some(d);
                break;
            }
        }
    }
    measurement::end_online();
    counters.matrix_rank = matrix.rank() as u64;
    counters.row_ops_uncharged = matrix.row_ops;
    counters.lp_pivots = eliminator.pivot_count() as u64;
    phase_group_ops.insert("relations".into(), rel_ops);
    phase_group_ops.insert("verification".into(), verify_ops);
    let verified_fast = verify_candidate
        && recovered.is_some_and(|d| instance.fast.mul_u64(instance.generator, d) == target);
    let verified_independent = verify_candidate
        && recovered.is_some_and(|d| independent_verify(instance, target, d).is_ok());
    Ok(LargePrimeReport {
        recovered,
        exhausted: recovered.is_none(),
        verified_fast,
        verified_independent,
        verification_performed: verify_candidate,
        config,
        factor_base: base.report,
        mitm: decomposer.report,
        counters,
        phase_group_ops,
        solve_wall_ns: start.elapsed().as_nanos() as u64,
        relation_samples: samples,
    })
}
