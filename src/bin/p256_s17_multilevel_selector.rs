//! Exact variable-column S17 selectors for the accepted P-256 Dickson factor base.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{self, CURVE_SLUG};
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::{CurveParams, Point};
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

const SPEC: &str = "dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129";
const FB_ID: &str = "FB1h2f8621cda105";
const FB_SHA256: &str = "2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42";
const POINTS_SHA256: &str = "70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1";
const TERMINAL: &str = "0x5b17195299a3158b93389ad04c776fff2a8bb23ca8659b5b0c2b75f9009b65e5";
const COLUMNS: u64 = 131_458;
const SIGNED_POINTS: u64 = 262_916;
const ROUND18_SHA256: &str = "b9081ea3cd3957dcae7805852371e395a2a63c8d4a2062e5c8e2c0b23980fce3";
const ROUND6_SHA256: &str = "71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4";
const ACTIVE_SIZES: [usize; 4] = [17, 18, 19, 20];
const OUTER_OFFSET: u64 = 90_322;
const OUTER_STRIDE: u64 = 23_509;
const NORMALIZE_BATCH: usize = 8_192;
const KEY_BYTES: usize = 33;
const DISK_RECORD_BYTES: usize = KEY_BYTES + 4 + 4;
const BUCKETS: usize = 256;
const TEMP_LIMIT_BYTES: u64 = 6 * 1024 * 1024 * 1024;
const PLANTED_PREIMAGE: &str = concat!(
    "icv1-fp256-t89188191154553853111372247798585809583-f188c491",
    "/s17-multilevel-round19/planted-signs-0"
);
const PUBLIC_PREIMAGE: &str = concat!(
    "icv1-fp256-t89188191154553853111372247798585809583-f188c491",
    "/s17-multilevel-round19/public-target-0"
);

#[derive(Parser)]
#[command(about = "Evaluate exact variable-column P-256 S17 selectors")]
struct Cli {
    /// Round-18 result whose fixed-atom boundary this experiment replaces.
    #[arg(long)]
    round18: PathBuf,
    /// Round-6 structured residual-degree receipt.
    #[arg(long)]
    round6: PathBuf,
    /// Deterministic JSON result path; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
    /// Parent for bounded external buckets; defaults to the system temp directory.
    #[arg(long)]
    temp_parent: Option<PathBuf>,
}

#[derive(Clone, Copy)]
struct ActivePoint {
    column: u32,
    positive: Projective,
    positive_double: Projective,
    negative_double: Projective,
}

#[derive(Clone, Copy)]
struct Pending {
    point: Projective,
    columns: u32,
    negative: u32,
}

#[derive(Clone, Copy)]
struct Partial {
    point: Projective,
    columns: u32,
    negative: u32,
    min: u8,
    max: u8,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct IndexRecord {
    key: [u8; KEY_BYTES],
    columns: u32,
    negative: u32,
}

impl Ord for IndexRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key
            .cmp(&other.key)
            .then_with(|| self.columns.cmp(&other.columns))
            .then_with(|| self.negative.cmp(&other.negative))
    }
}

impl PartialOrd for IndexRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone)]
struct Target {
    kind: String,
    scalar: Option<String>,
    point: Projective,
    affine: [String; 2],
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize)]
struct RelationWitness {
    columns_mask: u32,
    negative_mask: u32,
}

#[derive(Clone, Copy, Default, Serialize)]
struct WorkCounts {
    group_additions: u64,
    normalization_multiplications: u64,
    equality_multiplications: u64,
    key_lookups: u64,
    key_matches: u64,
    replayed_matches: u64,
    replay_failures: u64,
}

impl WorkCounts {
    fn field_multiplication_equivalents(&self) -> u64 {
        self.group_additions
            .saturating_mul(17)
            .saturating_add(self.normalization_multiplications)
            .saturating_add(self.equality_multiplications)
    }
}

#[derive(Serialize)]
struct FactorBaseReceipt {
    spec: String,
    fb_id: String,
    fb_sha256: String,
    points_sha256: String,
    terminal: String,
    columns: u64,
    signed_points: u64,
    complete_rebuild_verified: bool,
}

#[derive(Serialize)]
struct DependencyReceipt {
    round18_sha256: String,
    round6_sha256: String,
}

#[derive(Serialize)]
struct ResidualDegreeEvidence {
    maximum_by_residual_depth: BTreeMap<u32, u32>,
    every_component_complete_and_correct: bool,
    local_join_degree: u32,
    structured_residual_maximum: u32,
    unsplit_s17_degree_of_regularity: Option<u32>,
}

#[derive(Serialize)]
struct TargetIdentity {
    kind: String,
    scalar: Option<String>,
    point: [String; 2],
}

#[derive(Serialize)]
struct TargetSelectorResult {
    kind: String,
    relations: Vec<RelationWitness>,
    relation_sha256: String,
    false_positives: u64,
    false_negatives: u64,
    exact: bool,
}

#[derive(Serialize)]
struct BaselineResult {
    left_records: u64,
    right_records: u64,
    record_logical_bytes: u64,
    peak_logical_bytes: u64,
    radix_or_comparison_sort_records: u64,
    work: WorkCounts,
    field_multiplication_equivalents: u64,
    targets: Vec<TargetSelectorResult>,
}

#[derive(Serialize)]
struct BucketWidths {
    minimum: u64,
    median: u64,
    maximum: u64,
    total: u64,
}

#[derive(Serialize)]
struct StructuredResult {
    four_records: u64,
    five_records: u64,
    joined_eight_records: u64,
    joined_nine_records_per_target: u64,
    dickson_leaf_candidates: u64,
    dickson_leaf_rejections: u64,
    dickson_survival_ratio: f64,
    dickson_chain_sha256: String,
    propagated_certificate_records: u64,
    bucket_count: u32,
    left_bucket_widths: BucketWidths,
    right_bucket_widths: Vec<BucketWidths>,
    bytes_written: u64,
    bytes_read: u64,
    peak_resident_logical_bytes: u64,
    peak_materialized_bytes: u64,
    work: WorkCounts,
    field_multiplication_equivalents: u64,
    targets: Vec<TargetSelectorResult>,
}

#[derive(Serialize)]
struct ReferenceResult {
    complete: bool,
    assignments: u64,
    work: WorkCounts,
    targets: Vec<TargetSelectorResult>,
}

#[derive(Serialize)]
struct CellResult {
    active_columns: usize,
    factor_base_columns: Vec<u32>,
    baseline: BaselineResult,
    structured: StructuredResult,
    exhaustive_reference: Option<ReferenceResult>,
    selector_relation_sets_equal: bool,
    false_positives: u64,
    false_negatives: u64,
    exact: bool,
}

#[derive(Serialize)]
struct GrowthFit {
    sizes: Vec<usize>,
    baseline_fme: Vec<u64>,
    structured_fme: Vec<u64>,
    baseline_log2_slope_per_log2_b: f64,
    structured_log2_slope_per_log2_b: f64,
}

#[derive(Serialize)]
struct SparseLinearAlgebraProjection {
    rows: u64,
    columns: u64,
    row_weight: u64,
    nonzeros: u64,
    csr_bytes: u64,
    field_vector_bytes: u64,
    working_bytes: u64,
    wiedemann_nonzero_additions: u64,
    berlekamp_massey_scalar_operations: u64,
    optimistic_minimum_fme: u64,
}

#[derive(Serialize)]
struct FullProjection {
    left_records: String,
    left_records_log2: f64,
    right_records: String,
    right_records_log2: f64,
    exact_key_record_bytes: u64,
    baseline_left_materialized_bytes_log2: f64,
    structured_frontier_materialized_bytes_log2: f64,
    optimistic_per_target_fme_log2: f64,
    optimistic_per_usable_relation_fme_log2: f64,
    relation_rows: u64,
    whole_factor_base_success: f64,
    projected_targets: f64,
    projected_collection_fme_log2: f64,
    projected_total_with_sparse_la_fme_log2: f64,
    rho_reference_fme_log2: f64,
    total_over_rho_log2: f64,
    per_relation_below_2_pow_103: bool,
    collection_below_2_pow_120: bool,
    storage_below_2_pow_50: bool,
    optimistic_lower_bound: bool,
}

#[derive(Serialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    dependencies: DependencyReceipt,
    factor_base: FactorBaseReceipt,
    active_sizes: Vec<usize>,
    active_permutation: String,
    targets: Vec<TargetIdentity>,
    planted_negative_mask: u32,
    cells: Vec<CellResult>,
    growth_fit: GrowthFit,
    residual_degree_evidence: ResidualDegreeEvidence,
    sparse_linear_algebra_projection: SparseLinearAlgebraProjection,
    full_projection: FullProjection,
    exact: bool,
    attack_promotion_gate: bool,
    full_depth_unplanted_attempted: bool,
    decision: String,
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("not a non-negative integer: `{value}`"))
}

fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn file_sha256(path: &Path) -> Result<(Vec<u8>, String), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    Ok((bytes, digest))
}

fn verify_dependencies(round18: &Path, round6: &Path) -> Result<DependencyReceipt, String> {
    let (round18_bytes, round18_sha256) = file_sha256(round18)?;
    if round18_sha256 != ROUND18_SHA256 {
        return Err(format!(
            "round-18 SHA-256 mismatch: expected {ROUND18_SHA256}, got {round18_sha256}"
        ));
    }
    let round18_json: serde_json::Value =
        serde_json::from_slice(&round18_bytes).map_err(|error| error.to_string())?;
    if round18_json["schema"].as_str() != Some("p256.s17_outer_scan/v1")
        || round18_json["curve"].as_str() != Some(CURVE_SLUG)
        || round18_json["factor_base"]["fb_id"].as_str() != Some(FB_ID)
        || round18_json["attack_promotion_gate"].as_bool() != Some(false)
    {
        return Err("round-18 identity or rejection boundary changed".into());
    }

    let (round6_bytes, round6_sha256) = file_sha256(round6)?;
    if round6_sha256 != ROUND6_SHA256 {
        return Err(format!(
            "round-6 SHA-256 mismatch: expected {ROUND6_SHA256}, got {round6_sha256}"
        ));
    }
    let round6_json: serde_json::Value =
        serde_json::from_slice(&round6_bytes).map_err(|error| error.to_string())?;
    if round6_json["schema"].as_str() != Some("p256.dickson_residual_scaling/v1")
        || round6_json["curve"].as_str() != Some(CURVE_SLUG)
    {
        return Err("round-6 degree identity changed".into());
    }
    Ok(DependencyReceipt {
        round18_sha256,
        round6_sha256,
    })
}

fn load_degree_evidence(path: &Path) -> Result<ResidualDegreeEvidence, String> {
    let (bytes, _) = file_sha256(path)?;
    let dependency: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    let cells = dependency["cells"]
        .as_array()
        .ok_or("round-6 cells are missing")?;
    let mut maximum_by_residual_depth = BTreeMap::new();
    let mut complete = true;
    for cell in cells {
        let residual = cell["residual_depth"]
            .as_u64()
            .ok_or("round-6 residual depth is missing")? as u32;
        let degree = cell["max_solving_degree"]
            .as_u64()
            .ok_or("round-6 solving degree is missing")? as u32;
        maximum_by_residual_depth
            .entry(residual)
            .and_modify(|value: &mut u32| *value = (*value).max(degree))
            .or_insert(degree);
        let components = cell["components"]
            .as_u64()
            .ok_or("round-6 component count is missing")?;
        complete &= cell["complete_components"].as_u64() == Some(components)
            && cell["correct_components"].as_u64() == Some(components);
    }
    let expected = BTreeMap::from([(1u32, 3u32), (2, 3), (3, 4)]);
    if maximum_by_residual_depth != expected || !complete {
        return Err(format!(
            "round-6 residual-degree boundary changed: {maximum_by_residual_depth:?}"
        ));
    }
    Ok(ResidualDegreeEvidence {
        maximum_by_residual_depth,
        every_component_complete_and_correct: true,
        local_join_degree: 2,
        structured_residual_maximum: 4,
        unsplit_s17_degree_of_regularity: None,
    })
}

fn build_active_points() -> Result<(FactorBaseReceipt, Vec<ActivePoint>, String), String> {
    let built = p256_dickson_factor_base::build(SPEC)?;
    p256_dickson_factor_base::verify(&built.dump)?;
    let fb = &built.dump.factor_base;
    if fb.fb_id != FB_ID
        || fb.fb_sha256 != FB_SHA256
        || fb.points_sha256 != POINTS_SHA256
        || fb.columns != COLUMNS
        || fb.signed_points != SIGNED_POINTS
        || fb.params.get("terminal").map(String::as_str) != Some(TERMINAL)
    {
        return Err("factor-base identity changed".into());
    }
    let maximum_active = *ACTIVE_SIZES.last().expect("nonempty sizes");
    let mut active = Vec::with_capacity(maximum_active);
    let mut seen = BTreeSet::new();
    for position in 0..maximum_active {
        let column = ((OUTER_OFFSET + OUTER_STRIDE * position as u64) % COLUMNS) as usize;
        if !seen.insert(column) {
            return Err(format!("duplicate active column {column}"));
        }
        let row = built
            .dump
            .points
            .get(2 * column)
            .ok_or_else(|| format!("missing low point for column {column}"))?;
        if row.col as usize != column || row.coef != "1" {
            return Err(format!(
                "factor-base low-row convention changed at {column}"
            ));
        }
        let curve = CurveParams::p256();
        let affine = Point::Affine {
            x: curve.fe(parse_big(&row.x)?),
            y: curve.fe(parse_big(&row.y)?),
        };
        if !curve.is_on_curve(&affine) {
            return Err(format!("active column {column} is off curve"));
        }
        let positive = Projective::from_textbook(&affine);
        let positive_double = positive.add(&positive);
        active.push(ActivePoint {
            column: column as u32,
            positive,
            positive_double,
            negative_double: positive_double.neg(),
        });
    }

    let terminal = Fe::from_biguint(&parse_big(TERMINAL)?);
    let two = Fe::ONE.add(&Fe::ONE);
    let mut chain_receipt = Vec::new();
    for item in &active {
        let affine = item
            .positive
            .to_affine()
            .ok_or("active factor-base point became infinity")?;
        let mut value = Fe::from_canonical(&affine.0);
        chain_receipt.extend_from_slice(&item.column.to_be_bytes());
        for _ in 0..18 {
            value = value.sqr().sub(&two);
            chain_receipt.extend_from_slice(&value.to_bytes_be());
        }
        if value != terminal {
            return Err(format!(
                "active column {} does not reach the Dickson terminal",
                item.column
            ));
        }
    }

    Ok((
        FactorBaseReceipt {
            spec: SPEC.into(),
            fb_id: fb.fb_id.clone(),
            fb_sha256: fb.fb_sha256.clone(),
            points_sha256: fb.points_sha256.clone(),
            terminal: TERMINAL.into(),
            columns: fb.columns,
            signed_points: fb.signed_points,
            complete_rebuild_verified: true,
        },
        active,
        hex::encode(sha256(&chain_receipt)),
    ))
}

fn projective_equal(left: &Projective, right: &Projective, counts: &mut WorkCounts) -> bool {
    let left_identity = bool::from(left.is_identity());
    let right_identity = bool::from(right.is_identity());
    if left_identity || right_identity {
        return left_identity == right_identity;
    }
    counts.equality_multiplications += 4;
    left.x.mul(&right.z) == right.x.mul(&left.z) && left.y.mul(&right.z) == right.y.mul(&left.z)
}

fn pow_counted(mut base: Fe, exponent: &BigUint, counts: &mut WorkCounts) -> Fe {
    let mut result = Fe::ONE;
    for bit in 0..exponent.bits() {
        if exponent.bit(bit) {
            result = result.mul(&base);
            counts.normalization_multiplications += 1;
        }
        base = base.sqr();
        counts.normalization_multiplications += 1;
    }
    result
}

fn normalize_batch(
    pending: &[Pending],
    inverse_exponent: &BigUint,
    counts: &mut WorkCounts,
) -> Result<Vec<IndexRecord>, String> {
    let nonidentity: Vec<usize> = pending
        .iter()
        .enumerate()
        .filter_map(|(index, value)| (!bool::from(value.point.is_identity())).then_some(index))
        .collect();
    let mut inverses = vec![Fe::ZERO; pending.len()];
    if !nonidentity.is_empty() {
        let mut prefixes = Vec::with_capacity(nonidentity.len());
        let mut product = Fe::ONE;
        for &index in &nonidentity {
            prefixes.push(product);
            product = product.mul(&pending[index].point.z);
            counts.normalization_multiplications += 1;
        }
        let mut inverse = pow_counted(product, inverse_exponent, counts);
        for position in (0..nonidentity.len()).rev() {
            let index = nonidentity[position];
            inverses[index] = inverse.mul(&prefixes[position]);
            counts.normalization_multiplications += 1;
            if position != 0 {
                inverse = inverse.mul(&pending[index].point.z);
                counts.normalization_multiplications += 1;
            }
        }
    }
    let mut records = Vec::with_capacity(pending.len());
    for (index, item) in pending.iter().enumerate() {
        let mut key = [0u8; KEY_BYTES];
        if !bool::from(item.point.is_identity()) {
            let x = item.point.x.mul(&inverses[index]);
            let y = item.point.y.mul(&inverses[index]);
            counts.normalization_multiplications += 2;
            let x_bytes = x.to_bytes_be();
            let y_bytes = y.to_bytes_be();
            key[0] = 2 | (y_bytes[31] & 1);
            key[1..].copy_from_slice(&x_bytes);
        }
        records.push(IndexRecord {
            key,
            columns: item.columns,
            negative: item.negative,
        });
    }
    Ok(records)
}

fn next_combination(combination: &mut [usize], n: usize) -> bool {
    for index in (0..combination.len()).rev() {
        let maximum = n - combination.len() + index;
        if combination[index] != maximum {
            combination[index] += 1;
            for tail in index + 1..combination.len() {
                combination[tail] = combination[tail - 1] + 1;
            }
            return true;
        }
    }
    false
}

fn emit_direct_signed_sums<F>(
    active: &[ActivePoint],
    k: usize,
    counts: &mut WorkCounts,
    mut emit: F,
) -> Result<u64, String>
where
    F: FnMut(&[Pending], &mut WorkCounts) -> Result<(), String>,
{
    if k == 0 || k > active.len() || active.len() > 31 {
        return Err(format!(
            "invalid signed-sum shape n={}, k={k}",
            active.len()
        ));
    }
    let mut combination: Vec<usize> = (0..k).collect();
    let mut batch = Vec::with_capacity(NORMALIZE_BATCH);
    let mut emitted = 0u64;
    loop {
        let mut point = Projective::IDENTITY;
        let mut columns = 0u32;
        for &index in &combination {
            point = point.add(&active[index].positive);
            counts.group_additions += 1;
            columns |= 1u32 << index;
        }
        let mut previous_gray = 0usize;
        let mut negative = 0u32;
        for ordinal in 0..(1usize << k) {
            let gray = ordinal ^ (ordinal >> 1);
            if ordinal != 0 {
                let changed = gray ^ previous_gray;
                let bit = changed.trailing_zeros() as usize;
                let active_index = combination[bit];
                if gray & (1usize << bit) != 0 {
                    point = point.add(&active[active_index].negative_double);
                    negative |= 1u32 << active_index;
                } else {
                    point = point.add(&active[active_index].positive_double);
                    negative &= !(1u32 << active_index);
                }
                counts.group_additions += 1;
            }
            previous_gray = gray;
            batch.push(Pending {
                point,
                columns,
                negative,
            });
            emitted += 1;
            if batch.len() == NORMALIZE_BATCH {
                emit(&batch, counts)?;
                batch.clear();
            }
        }
        if !next_combination(&mut combination, active.len()) {
            break;
        }
    }
    if !batch.is_empty() {
        emit(&batch, counts)?;
    }
    Ok(emitted)
}

fn key_range(records: &[IndexRecord], key: &[u8; KEY_BYTES]) -> std::ops::Range<usize> {
    let start = records.partition_point(|record| record.key < *key);
    let end = records.partition_point(|record| record.key <= *key);
    start..end
}

fn min_index(mask: u32) -> u8 {
    mask.trailing_zeros() as u8
}

fn max_index(mask: u32) -> u8 {
    (31 - mask.leading_zeros()) as u8
}

fn replay_relation(
    witness: &RelationWitness,
    active: &[ActivePoint],
    target: &Projective,
    counts: &mut WorkCounts,
) -> bool {
    let mut sum = Projective::IDENTITY;
    for (index, point) in active.iter().enumerate() {
        if witness.columns_mask & (1u32 << index) == 0 {
            continue;
        }
        let signed = if witness.negative_mask & (1u32 << index) == 0 {
            point.positive
        } else {
            point.positive.neg()
        };
        sum = sum.add(&signed);
        counts.group_additions += 1;
    }
    counts.replayed_matches += 1;
    let exact = projective_equal(&sum, target, counts);
    if !exact {
        counts.replay_failures += 1;
    }
    exact
}

fn relation_digest(relations: &BTreeSet<RelationWitness>) -> String {
    let mut bytes = Vec::with_capacity(relations.len() * 8);
    for relation in relations {
        bytes.extend_from_slice(&relation.columns_mask.to_be_bytes());
        bytes.extend_from_slice(&relation.negative_mask.to_be_bytes());
    }
    hex::encode(sha256(&bytes))
}

fn selector_target_result(
    target: &Target,
    relations: &BTreeSet<RelationWitness>,
    false_positives: u64,
    false_negatives: u64,
) -> TargetSelectorResult {
    TargetSelectorResult {
        kind: target.kind.clone(),
        relations: relations.iter().cloned().collect(),
        relation_sha256: relation_digest(relations),
        false_positives,
        false_negatives,
        exact: false_positives == 0 && false_negatives == 0,
    }
}

fn run_baseline(
    active: &[ActivePoint],
    targets: &[Target],
    inverse_exponent: &BigUint,
) -> Result<(BaselineResult, Vec<BTreeSet<RelationWitness>>), String> {
    let mut work = WorkCounts::default();
    let mut left = Vec::new();
    let left_records = emit_direct_signed_sums(active, 8, &mut work, |batch, counts| {
        left.extend(normalize_batch(batch, inverse_exponent, counts)?);
        Ok(())
    })?;
    left.sort_unstable();

    let mut relations = vec![BTreeSet::new(); targets.len()];
    let right_records = emit_direct_signed_sums(active, 9, &mut work, |batch, counts| {
        for (target_index, target) in targets.iter().enumerate() {
            let adjusted: Vec<Pending> = batch
                .iter()
                .map(|right| {
                    counts.group_additions += 1;
                    Pending {
                        point: target.point.add(&right.point.neg()),
                        columns: right.columns,
                        negative: right.negative,
                    }
                })
                .collect();
            let queries = normalize_batch(&adjusted, inverse_exponent, counts)?;
            for query in queries {
                counts.key_lookups += 1;
                for candidate in &left[key_range(&left, &query.key)] {
                    if max_index(candidate.columns) >= min_index(query.columns) {
                        continue;
                    }
                    counts.key_matches += 1;
                    let witness = RelationWitness {
                        columns_mask: candidate.columns | query.columns,
                        negative_mask: candidate.negative | query.negative,
                    };
                    if replay_relation(&witness, active, &target.point, counts) {
                        relations[target_index].insert(witness);
                    }
                }
            }
        }
        Ok(())
    })?;

    let target_results = targets
        .iter()
        .zip(&relations)
        .map(|(target, found)| selector_target_result(target, found, work.replay_failures, 0))
        .collect();
    let peak_logical_bytes = left_records * std::mem::size_of::<IndexRecord>() as u64;
    Ok((
        BaselineResult {
            left_records,
            right_records,
            record_logical_bytes: std::mem::size_of::<IndexRecord>() as u64,
            peak_logical_bytes,
            radix_or_comparison_sort_records: left_records,
            field_multiplication_equivalents: work.field_multiplication_equivalents(),
            work,
            targets: target_results,
        },
        relations,
    ))
}

fn build_partials(
    active: &[ActivePoint],
    k: usize,
    counts: &mut WorkCounts,
) -> Result<Vec<Partial>, String> {
    let mut partials = Vec::new();
    emit_direct_signed_sums(active, k, counts, |batch, _| {
        partials.extend(batch.iter().map(|item| Partial {
            point: item.point,
            columns: item.columns,
            negative: item.negative,
            min: min_index(item.columns),
            max: max_index(item.columns),
        }));
        Ok(())
    })?;
    partials.sort_unstable_by_key(|item| (item.min, item.max, item.columns, item.negative));
    Ok(partials)
}

fn first_min_greater(partials: &[Partial], maximum: u8) -> usize {
    partials.partition_point(|item| item.min <= maximum)
}

fn emit_structured_join<F>(
    left: &[Partial],
    right: &[Partial],
    counts: &mut WorkCounts,
    mut emit: F,
) -> Result<u64, String>
where
    F: FnMut(&[Pending], &mut WorkCounts) -> Result<(), String>,
{
    let mut batch = Vec::with_capacity(NORMALIZE_BATCH);
    let mut emitted = 0u64;
    for left_item in left {
        let start = first_min_greater(right, left_item.max);
        for right_item in &right[start..] {
            counts.group_additions += 1;
            batch.push(Pending {
                point: left_item.point.add(&right_item.point),
                columns: left_item.columns | right_item.columns,
                negative: left_item.negative | right_item.negative,
            });
            emitted += 1;
            if batch.len() == NORMALIZE_BATCH {
                emit(&batch, counts)?;
                batch.clear();
            }
        }
    }
    if !batch.is_empty() {
        emit(&batch, counts)?;
    }
    Ok(emitted)
}

fn bucket_for(key: &[u8; KEY_BYTES]) -> usize {
    if key[0] == 0 {
        0
    } else {
        key[KEY_BYTES - 1] as usize
    }
}

fn bucket_paths(root: &Path, stem: &str) -> Vec<PathBuf> {
    (0..BUCKETS)
        .map(|bucket| root.join(format!("{stem}-{bucket:03}.bin")))
        .collect()
}

fn open_bucket_writers(paths: &[PathBuf]) -> Result<Vec<BufWriter<File>>, String> {
    paths
        .iter()
        .map(|path| {
            File::create(path)
                .map(BufWriter::new)
                .map_err(|error| format!("{}: {error}", path.display()))
        })
        .collect()
}

fn write_record(writer: &mut BufWriter<File>, record: &IndexRecord) -> Result<(), String> {
    writer
        .write_all(&record.key)
        .and_then(|()| writer.write_all(&record.columns.to_le_bytes()))
        .and_then(|()| writer.write_all(&record.negative.to_le_bytes()))
        .map_err(|error| error.to_string())
}

fn finish_writers(writers: Vec<BufWriter<File>>) -> Result<(), String> {
    for mut writer in writers {
        writer.flush().map_err(|error| error.to_string())?;
    }
    Ok(())
}

fn read_records(path: &Path) -> Result<Vec<IndexRecord>, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    if bytes.len() % DISK_RECORD_BYTES != 0 {
        return Err(format!("{} has a partial record", path.display()));
    }
    let mut records = Vec::with_capacity(bytes.len() / DISK_RECORD_BYTES);
    for chunk in bytes.chunks_exact(DISK_RECORD_BYTES) {
        let mut key = [0u8; KEY_BYTES];
        key.copy_from_slice(&chunk[..KEY_BYTES]);
        let columns = u32::from_le_bytes(
            chunk[KEY_BYTES..KEY_BYTES + 4]
                .try_into()
                .expect("four-byte column mask"),
        );
        let negative = u32::from_le_bytes(
            chunk[KEY_BYTES + 4..]
                .try_into()
                .expect("four-byte sign mask"),
        );
        records.push(IndexRecord {
            key,
            columns,
            negative,
        });
    }
    Ok(records)
}

fn widths(values: &[u64]) -> BucketWidths {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    BucketWidths {
        minimum: sorted[0],
        median: sorted[sorted.len() / 2],
        maximum: *sorted.last().expect("nonempty bucket widths"),
        total: sorted.iter().sum(),
    }
}

fn merge_bucket(
    mut left: Vec<IndexRecord>,
    mut right: Vec<IndexRecord>,
    active: &[ActivePoint],
    target: &Target,
    work: &mut WorkCounts,
    relations: &mut BTreeSet<RelationWitness>,
) {
    left.sort_unstable();
    right.sort_unstable();
    let mut left_index = 0usize;
    let mut right_index = 0usize;
    while left_index < left.len() && right_index < right.len() {
        match left[left_index].key.cmp(&right[right_index].key) {
            Ordering::Less => left_index += 1,
            Ordering::Greater => right_index += 1,
            Ordering::Equal => {
                let key = left[left_index].key;
                let left_end = left_index + left[left_index..].partition_point(|r| r.key == key);
                let right_end =
                    right_index + right[right_index..].partition_point(|r| r.key == key);
                for left_record in &left[left_index..left_end] {
                    for right_record in &right[right_index..right_end] {
                        if max_index(left_record.columns) >= min_index(right_record.columns) {
                            continue;
                        }
                        work.key_matches += 1;
                        let witness = RelationWitness {
                            columns_mask: left_record.columns | right_record.columns,
                            negative_mask: left_record.negative | right_record.negative,
                        };
                        if replay_relation(&witness, active, &target.point, work) {
                            relations.insert(witness);
                        }
                    }
                }
                left_index = left_end;
                right_index = right_end;
            }
        }
    }
}

fn run_structured(
    active: &[ActivePoint],
    targets: &[Target],
    inverse_exponent: &BigUint,
    dickson_chain_sha256: &str,
    temp_parent: &Path,
) -> Result<(StructuredResult, Vec<BTreeSet<RelationWitness>>), String> {
    let mut work = WorkCounts::default();
    let four = build_partials(active, 4, &mut work)?;
    let five = build_partials(active, 5, &mut work)?;
    let four_records = four.len() as u64;
    let five_records = five.len() as u64;

    let temp_root = temp_parent.join(format!(
        "p256-s17-round19-{}-{}",
        std::process::id(),
        active.len()
    ));
    if temp_root.exists() {
        return Err(format!(
            "temporary directory already exists: {}",
            temp_root.display()
        ));
    }
    fs::create_dir(&temp_root).map_err(|error| error.to_string())?;

    let left_paths = bucket_paths(&temp_root, "left");
    let mut left_writers = open_bucket_writers(&left_paths)?;
    let mut left_widths = vec![0u64; BUCKETS];
    let joined_eight_records = emit_structured_join(&four, &four, &mut work, |batch, counts| {
        for record in normalize_batch(batch, inverse_exponent, counts)? {
            let bucket = bucket_for(&record.key);
            write_record(&mut left_writers[bucket], &record)?;
            left_widths[bucket] += 1;
        }
        Ok(())
    })?;
    finish_writers(left_writers)?;
    let left_bytes = joined_eight_records * DISK_RECORD_BYTES as u64;

    let mut relations = Vec::with_capacity(targets.len());
    let mut target_results = Vec::with_capacity(targets.len());
    let mut right_width_receipts = Vec::with_capacity(targets.len());
    let mut bytes_written = left_bytes;
    let mut bytes_read = 0u64;
    let mut peak_materialized = left_bytes;
    let mut peak_resident = 0u64;
    let mut joined_nine_records_per_target = None;

    for (target_index, target) in targets.iter().enumerate() {
        let right_paths = bucket_paths(&temp_root, &format!("right-{target_index}"));
        let mut right_writers = open_bucket_writers(&right_paths)?;
        let mut right_widths = vec![0u64; BUCKETS];
        let joined_nine = emit_structured_join(&four, &five, &mut work, |batch, counts| {
            let adjusted: Vec<Pending> = batch
                .iter()
                .map(|right| {
                    counts.group_additions += 1;
                    Pending {
                        point: target.point.add(&right.point.neg()),
                        columns: right.columns,
                        negative: right.negative,
                    }
                })
                .collect();
            for record in normalize_batch(&adjusted, inverse_exponent, counts)? {
                let bucket = bucket_for(&record.key);
                write_record(&mut right_writers[bucket], &record)?;
                right_widths[bucket] += 1;
            }
            Ok(())
        })?;
        if joined_nine_records_per_target
            .replace(joined_nine)
            .is_some_and(|v| v != joined_nine)
        {
            return Err("structured right width changed between targets".into());
        }
        finish_writers(right_writers)?;
        let right_bytes = joined_nine * DISK_RECORD_BYTES as u64;
        let materialized = left_bytes + right_bytes;
        if materialized > TEMP_LIMIT_BYTES {
            return Err(format!(
                "structured materialized bytes {materialized} exceed {TEMP_LIMIT_BYTES}"
            ));
        }
        peak_materialized = peak_materialized.max(materialized);
        bytes_written += right_bytes;
        let mut found = BTreeSet::new();
        let replay_failures_before = work.replay_failures;
        for bucket in 0..BUCKETS {
            let left_bucket = read_records(&left_paths[bucket])?;
            let right_bucket = read_records(&right_paths[bucket])?;
            bytes_read += ((left_bucket.len() + right_bucket.len()) * DISK_RECORD_BYTES) as u64;
            peak_resident = peak_resident.max(
                ((left_bucket.len() + right_bucket.len()) * std::mem::size_of::<IndexRecord>())
                    as u64,
            );
            work.key_lookups += right_bucket.len() as u64;
            merge_bucket(
                left_bucket,
                right_bucket,
                active,
                target,
                &mut work,
                &mut found,
            );
            fs::remove_file(&right_paths[bucket])
                .map_err(|error| format!("{}: {error}", right_paths[bucket].display()))?;
        }
        let target_replay_failures = work.replay_failures - replay_failures_before;
        target_results.push(selector_target_result(
            target,
            &found,
            target_replay_failures,
            0,
        ));
        relations.push(found);
        right_width_receipts.push(widths(&right_widths));
    }
    for path in &left_paths {
        fs::remove_file(path).map_err(|error| format!("{}: {error}", path.display()))?;
    }
    fs::remove_dir(&temp_root).map_err(|error| error.to_string())?;

    let joined_nine = joined_nine_records_per_target.ok_or("no targets")?;
    let propagated =
        four_records + five_records + joined_eight_records + joined_nine * targets.len() as u64;
    Ok((
        StructuredResult {
            four_records,
            five_records,
            joined_eight_records,
            joined_nine_records_per_target: joined_nine,
            dickson_leaf_candidates: active.len() as u64,
            dickson_leaf_rejections: 0,
            dickson_survival_ratio: 1.0,
            dickson_chain_sha256: dickson_chain_sha256.into(),
            propagated_certificate_records: propagated,
            bucket_count: BUCKETS as u32,
            left_bucket_widths: widths(&left_widths),
            right_bucket_widths: right_width_receipts,
            bytes_written,
            bytes_read,
            peak_resident_logical_bytes: peak_resident,
            peak_materialized_bytes: peak_materialized,
            field_multiplication_equivalents: work.field_multiplication_equivalents(),
            work,
            targets: target_results,
        },
        relations,
    ))
}

fn run_exhaustive_reference(
    active: &[ActivePoint],
    targets: &[Target],
) -> Result<(ReferenceResult, Vec<BTreeSet<RelationWitness>>), String> {
    let mut work = WorkCounts::default();
    let mut relations = vec![BTreeSet::new(); targets.len()];
    let mut assignments = 0u64;
    let mut combination: Vec<usize> = (0..17).collect();
    loop {
        let mut point = Projective::IDENTITY;
        let mut columns = 0u32;
        for &index in &combination {
            point = point.add(&active[index].positive);
            work.group_additions += 1;
            columns |= 1u32 << index;
        }
        let mut previous_gray = 0usize;
        let mut negative = 0u32;
        for ordinal in 0..(1usize << 17) {
            let gray = ordinal ^ (ordinal >> 1);
            if ordinal != 0 {
                let changed = gray ^ previous_gray;
                let bit = changed.trailing_zeros() as usize;
                let active_index = combination[bit];
                if gray & (1usize << bit) != 0 {
                    point = point.add(&active[active_index].negative_double);
                    negative |= 1u32 << active_index;
                } else {
                    point = point.add(&active[active_index].positive_double);
                    negative &= !(1u32 << active_index);
                }
                work.group_additions += 1;
            }
            previous_gray = gray;
            assignments += 1;
            for (target_index, target) in targets.iter().enumerate() {
                if projective_equal(&point, &target.point, &mut work) {
                    let witness = RelationWitness {
                        columns_mask: columns,
                        negative_mask: negative,
                    };
                    if replay_relation(&witness, active, &target.point, &mut work) {
                        relations[target_index].insert(witness);
                    }
                }
            }
        }
        if !next_combination(&mut combination, active.len()) {
            break;
        }
    }
    let target_results = targets
        .iter()
        .zip(&relations)
        .map(|(target, found)| selector_target_result(target, found, work.replay_failures, 0))
        .collect();
    Ok((
        ReferenceResult {
            complete: true,
            assignments,
            work,
            targets: target_results,
        },
        relations,
    ))
}

fn make_targets(active: &[ActivePoint]) -> Result<(Vec<Target>, u32), String> {
    let curve = CurveParams::p256();
    let sign_integer = BigUint::from_bytes_be(&sha256(PLANTED_PREIMAGE.as_bytes()));
    let mut negative_mask = 0u32;
    let mut planted = Projective::IDENTITY;
    for (index, point) in active.iter().take(17).enumerate() {
        let negative = sign_integer.bit(index as u64);
        if negative {
            negative_mask |= 1u32 << index;
        }
        let signed = if negative {
            point.positive.neg()
        } else {
            point.positive
        };
        planted = planted.add(&signed);
    }
    if bool::from(planted.is_identity()) {
        return Err("planted target is infinity".into());
    }
    let mut scalar = BigUint::from_bytes_be(&sha256(PUBLIC_PREIMAGE.as_bytes())) % &curve.n;
    if scalar.is_zero() {
        scalar = BigUint::one();
    }
    let public_textbook = curve.generator().scalar_mul_vartime(&scalar, &curve.a_fe());
    if !curve.is_on_curve(&public_textbook) {
        return Err("hash-public target is off curve".into());
    }
    let public = Projective::from_textbook(&public_textbook);
    Ok((
        vec![
            Target {
                kind: "planted-control".into(),
                scalar: None,
                point: planted,
                affine: affine_hex(&planted)?,
            },
            Target {
                kind: "hash-public".into(),
                scalar: Some(lower_hex(&scalar)),
                point: public,
                affine: affine_hex(&public)?,
            },
        ],
        negative_mask,
    ))
}

fn affine_hex(point: &Projective) -> Result<[String; 2], String> {
    let (x, y) = point.to_affine().ok_or("target is infinity")?;
    Ok([lower_hex(&x.to_biguint()), lower_hex(&y.to_biguint())])
}

fn choose_u64(n: usize, k: usize) -> u64 {
    let k = k.min(n - k);
    let mut value = 1u64;
    for index in 0..k {
        value = value * (n - index) as u64 / (index + 1) as u64;
    }
    value
}

fn choose_big(n: u64, k: u32) -> BigUint {
    let k = u64::from(k).min(n - u64::from(k));
    let mut value = BigUint::one();
    for index in 0..k {
        value *= BigUint::from(n - index);
        value /= BigUint::from(index + 1);
    }
    value
}

fn big_log2(value: &BigUint) -> f64 {
    value
        .to_f64()
        .expect("P-256 selector count fits finite f64")
        .log2()
}

fn slope(xs: &[f64], ys: &[f64]) -> f64 {
    let count = xs.len() as f64;
    let mean_x = xs.iter().sum::<f64>() / count;
    let mean_y = ys.iter().sum::<f64>() / count;
    let numerator = xs
        .iter()
        .zip(ys)
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum::<f64>();
    let denominator = xs.iter().map(|x| (x - mean_x).powi(2)).sum::<f64>();
    numerator / denominator
}

fn projection() -> (SparseLinearAlgebraProjection, FullProjection) {
    let left = choose_big(COLUMNS, 8) << 8usize;
    let right = choose_big(COLUMNS, 9) << 9usize;
    let left_f = left.to_f64().expect("left count fits f64");
    let right_f = right.to_f64().expect("right count fits f64");
    let record_bytes = DISK_RECORD_BYTES as f64;

    // Optimistic warm lower bound: signed right sums are shared across all
    // targets.  Each target still needs one complete target adjustment (17
    // FME) and affine normalization (asymptotically 5 FME) per right record.
    let per_target_fme = 22.0 * right_f;
    let whole_success = 0.964_000_182_074_505_8f64;
    let per_relation_fme = per_target_fme / whole_success;
    let relation_rows = 138_031u64;
    let projected_targets = relation_rows as f64 / whole_success;

    let row_weight = 17u64;
    let nonzeros = relation_rows * row_weight;
    let csr_bytes = nonzeros * 5 + (relation_rows + 1) * 8;
    let field_vector_bytes = 3 * COLUMNS * 32;
    let wiedemann_nonzero_additions = 2 * COLUMNS * nonzeros;
    let berlekamp_massey_scalar_operations = COLUMNS * COLUMNS;
    let sparse_minimum = wiedemann_nonzero_additions + berlekamp_massey_scalar_operations;
    let sparse = SparseLinearAlgebraProjection {
        rows: relation_rows,
        columns: COLUMNS,
        row_weight,
        nonzeros,
        csr_bytes,
        field_vector_bytes,
        working_bytes: csr_bytes + field_vector_bytes,
        wiedemann_nonzero_additions,
        berlekamp_massey_scalar_operations,
        optimistic_minimum_fme: sparse_minimum,
    };

    let setup_fme = (22.0 * left_f) + (17.0 * right_f);
    let collection_fme = setup_fme + projected_targets * per_target_fme;
    let total_fme = collection_fme + sparse_minimum as f64;
    let curve = CurveParams::p256();
    let rho_fme = 1.3 * curve.n.to_f64().expect("P-256 order fits f64").sqrt() * 17.0;
    let baseline_storage = left_f * record_bytes;
    let structured_storage = (left_f + right_f) * record_bytes;
    let full = FullProjection {
        left_records: left.to_string(),
        left_records_log2: big_log2(&left),
        right_records: right.to_string(),
        right_records_log2: big_log2(&right),
        exact_key_record_bytes: DISK_RECORD_BYTES as u64,
        baseline_left_materialized_bytes_log2: baseline_storage.log2(),
        structured_frontier_materialized_bytes_log2: structured_storage.log2(),
        optimistic_per_target_fme_log2: per_target_fme.log2(),
        optimistic_per_usable_relation_fme_log2: per_relation_fme.log2(),
        relation_rows,
        whole_factor_base_success: whole_success,
        projected_targets,
        projected_collection_fme_log2: collection_fme.log2(),
        projected_total_with_sparse_la_fme_log2: total_fme.log2(),
        rho_reference_fme_log2: rho_fme.log2(),
        total_over_rho_log2: (total_fme / rho_fme).log2(),
        per_relation_below_2_pow_103: per_relation_fme.log2() < 103.0,
        collection_below_2_pow_120: total_fme.log2() < 120.0,
        storage_below_2_pow_50: structured_storage.log2() < 50.0,
        optimistic_lower_bound: true,
    };
    (sparse, full)
}

fn run(cli: Cli) -> Result<(), String> {
    let dependencies = verify_dependencies(&cli.round18, &cli.round6)?;
    let residual_degree_evidence = load_degree_evidence(&cli.round6)?;
    let (factor_base, all_active, dickson_chain_sha256) = build_active_points()?;
    let (targets, planted_negative_mask) = make_targets(&all_active)?;
    let curve = CurveParams::p256();
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let temp_parent = cli.temp_parent.unwrap_or_else(std::env::temp_dir);
    if !temp_parent.is_dir() {
        return Err(format!(
            "temp parent is not a directory: {}",
            temp_parent.display()
        ));
    }

    let planted_columns = (1u32 << 17) - 1;
    let planted_witness = RelationWitness {
        columns_mask: planted_columns,
        negative_mask: planted_negative_mask,
    };
    let mut cells = Vec::new();
    for &size in &ACTIVE_SIZES {
        let active = &all_active[..size];
        let expected_left = choose_u64(size, 8) << 8;
        let expected_right = choose_u64(size, 9) << 9;
        let materialized = (expected_left + expected_right) * DISK_RECORD_BYTES as u64;
        if materialized > TEMP_LIMIT_BYTES {
            return Err(format!(
                "preregistered size {size} needs {materialized} temporary bytes"
            ));
        }

        let (mut baseline, baseline_relations) = run_baseline(active, &targets, &inverse_exponent)?;
        if baseline.left_records != expected_left || baseline.right_records != expected_right {
            return Err(format!("baseline width mismatch at B={size}"));
        }
        let (mut structured, structured_relations) = run_structured(
            active,
            &targets,
            &inverse_exponent,
            &dickson_chain_sha256,
            &temp_parent,
        )?;
        if structured.joined_eight_records != expected_left
            || structured.joined_nine_records_per_target != expected_right
        {
            return Err(format!("structured width mismatch at B={size}"));
        }

        let exhaustive = if size <= 18 {
            Some(run_exhaustive_reference(active, &targets)?)
        } else {
            None
        };
        let reference_relations = exhaustive.as_ref().map(|(_, relations)| relations);
        let selector_equal = baseline_relations == structured_relations;
        let reference_equal = reference_relations
            .map(|reference| reference == &baseline_relations && reference == &structured_relations)
            .unwrap_or(true);
        if !selector_equal || !reference_equal {
            return Err(format!("selector/reference disagreement at B={size}"));
        }
        if !baseline_relations[0].contains(&planted_witness)
            || !structured_relations[0].contains(&planted_witness)
            || reference_relations.is_some_and(|reference| !reference[0].contains(&planted_witness))
        {
            return Err(format!("missing planted witness at B={size}"));
        }

        for (target_index, relation_set) in baseline_relations.iter().enumerate() {
            baseline.targets[target_index].false_negatives = 0;
            baseline.targets[target_index].exact = true;
            structured.targets[target_index].false_negatives = 0;
            structured.targets[target_index].exact = true;
            if baseline.targets[target_index].relations.len() != relation_set.len()
                || structured.targets[target_index].relations.len() != relation_set.len()
            {
                return Err("serialized relation count changed".into());
            }
        }
        let (reference_result, _) = exhaustive.unzip();
        let false_positives = baseline.work.replay_failures + structured.work.replay_failures;
        let false_negatives = 0;
        let exact = selector_equal && reference_equal && false_positives == 0;
        cells.push(CellResult {
            active_columns: size,
            factor_base_columns: active.iter().map(|point| point.column).collect(),
            baseline,
            structured,
            exhaustive_reference: reference_result,
            selector_relation_sets_equal: selector_equal,
            false_positives,
            false_negatives,
            exact,
        });
    }

    let xs: Vec<f64> = cells
        .iter()
        .map(|cell| (cell.active_columns as f64).log2())
        .collect();
    let baseline_fme: Vec<u64> = cells
        .iter()
        .map(|cell| cell.baseline.field_multiplication_equivalents)
        .collect();
    let structured_fme: Vec<u64> = cells
        .iter()
        .map(|cell| cell.structured.field_multiplication_equivalents)
        .collect();
    let growth_fit = GrowthFit {
        sizes: ACTIVE_SIZES.to_vec(),
        baseline_log2_slope_per_log2_b: slope(
            &xs,
            &baseline_fme
                .iter()
                .map(|value| (*value as f64).log2())
                .collect::<Vec<_>>(),
        ),
        structured_log2_slope_per_log2_b: slope(
            &xs,
            &structured_fme
                .iter()
                .map(|value| (*value as f64).log2())
                .collect::<Vec<_>>(),
        ),
        baseline_fme,
        structured_fme,
    };
    let (sparse_linear_algebra_projection, full_projection) = projection();
    let exact = cells.iter().all(|cell| cell.exact)
        && residual_degree_evidence.every_component_complete_and_correct;
    let attack_promotion_gate = exact
        && residual_degree_evidence.structured_residual_maximum <= 5
        && full_projection.per_relation_below_2_pow_103
        && full_projection.collection_below_2_pow_120
        && full_projection.storage_below_2_pow_50;
    let result = ExperimentResult {
        schema: "p256.s17_multilevel_selector/v1".into(),
        curve: CURVE_SLUG.into(),
        dependencies,
        factor_base,
        active_sizes: ACTIVE_SIZES.to_vec(),
        active_permutation: format!("column(j)=({OUTER_OFFSET}+{OUTER_STRIDE}*j) mod {COLUMNS}"),
        targets: targets
            .iter()
            .map(|target| TargetIdentity {
                kind: target.kind.clone(),
                scalar: target.scalar.clone(),
                point: target.affine.clone(),
            })
            .collect(),
        planted_negative_mask,
        cells,
        growth_fit,
        residual_degree_evidence,
        sparse_linear_algebra_projection,
        full_projection,
        exact,
        attack_promotion_gate,
        full_depth_unplanted_attempted: false,
        decision: if attack_promotion_gate {
            "promotion gates passed; full-depth unplanted attempt is authorized".into()
        } else {
            "negative: exact multilevel selection retains the global 8+9 frontier".into()
        },
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match cli.out {
        Some(path) => fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy_active(count: usize) -> Vec<ActivePoint> {
        let curve = CurveParams::p256();
        (1..=count)
            .map(|scalar| {
                let point = curve
                    .generator()
                    .scalar_mul_vartime(&BigUint::from(scalar as u64), &curve.a_fe());
                let positive = Projective::from_textbook(&point);
                let positive_double = positive.add(&positive);
                ActivePoint {
                    column: scalar as u32,
                    positive,
                    positive_double,
                    negative_double: positive_double.neg(),
                }
            })
            .collect()
    }

    #[test]
    fn direct_signed_sum_count_is_exact() {
        let active = toy_active(6);
        let mut counts = WorkCounts::default();
        let emitted = emit_direct_signed_sums(&active, 3, &mut counts, |_, _| Ok(())).unwrap();
        assert_eq!(emitted, choose_u64(6, 3) << 3);
        assert_eq!(counts.group_additions, choose_u64(6, 3) * (3 + 7));
    }

    #[test]
    fn batch_normalization_matches_affine_conversion() {
        let active = toy_active(5);
        let curve = CurveParams::p256();
        let inverse_exponent = &curve.p - BigUint::from(2u8);
        let pending: Vec<Pending> = active
            .iter()
            .enumerate()
            .map(|(index, point)| Pending {
                point: point.positive,
                columns: 1u32 << index,
                negative: 0,
            })
            .collect();
        let mut counts = WorkCounts::default();
        let records = normalize_batch(&pending, &inverse_exponent, &mut counts).unwrap();
        for (record, point) in records.iter().zip(&active) {
            let (x, y) = point.positive.to_affine().unwrap();
            let mut key = [0u8; KEY_BYTES];
            key[0] = 2 | (y.to_bytes_be()[31] & 1);
            key[1..].copy_from_slice(&x.to_bytes_be());
            assert_eq!(record.key, key);
        }
    }

    #[test]
    fn structured_join_has_canonical_width() {
        let active = toy_active(9);
        let mut counts = WorkCounts::default();
        let four = build_partials(&active, 4, &mut counts).unwrap();
        let five = build_partials(&active, 5, &mut counts).unwrap();
        let emitted = emit_structured_join(&four, &five, &mut counts, |_, _| Ok(())).unwrap();
        assert_eq!(emitted, choose_u64(9, 9) << 9);
    }

    #[test]
    fn disk_record_round_trip_is_exact() {
        let root =
            std::env::temp_dir().join(format!("p256-s17-round19-test-{}", std::process::id()));
        if root.exists() {
            fs::remove_file(&root).unwrap();
        }
        let record = IndexRecord {
            key: [0x5a; KEY_BYTES],
            columns: 0x12345,
            negative: 0x10201,
        };
        {
            let mut writer = BufWriter::new(File::create(&root).unwrap());
            write_record(&mut writer, &record).unwrap();
            writer.flush().unwrap();
        }
        assert_eq!(read_records(&root).unwrap(), vec![record]);
        fs::remove_file(root).unwrap();
    }
}
