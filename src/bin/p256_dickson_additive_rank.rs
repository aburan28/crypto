//! Round 312: search the registered P-256 Dickson factor base for exact
//! low-weight additive identities that reduce its independent-log rank.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::mem::size_of;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{self, CURVE_SLUG};
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::Zero;
use serde::Serialize;
use serde_json::{json, Value};

const ROUND2_SHA256: &str = "dcaa1f66f5f89757a3670a4a8e158ca6d5c360b7ca119d384d2c4a0b2794c2c8";
const ROUND19_SHA256: &str = "3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114";
const ROUND24_SHA256: &str = "d1e2fe33e1abbae8a9bf6885ca7cfef8012233363f4cd43e0d30d468a381fc9c";
const ROUND311_SHA256: &str = "9033d28dca7c42a0464ae11e5eb06052c0501cc79533b96a13dd148e7dfd0953";
const FACTOR_BASE: &str = "FB1h2f8621cda105";
const FACTOR_BASE_SHA256: &str = "2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42";
const POINTS_SHA256: &str = "70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1";
const SPEC: &str = "dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129";
const COLUMNS: usize = 131_458;
const COMPLETE_LEVELS: usize = 7;
const DEPTH: usize = 18;
const DEEP_SAMPLES_PER_LEVEL: usize = 4_096;
const MULTIPLIERS: [u32; 6] = [1, 2, 3, 4, 8, 16];
const NORMALIZE_CHUNK: usize = 65_536;
const RANK_MAX_RELATIONS: usize = 16_384;
const RANK_MAX_ACTIVE_COLUMNS: usize = 65_536;
const REQUESTED_ROWS: u64 = 138_031;
const RANK_ALLOWANCE_NUMERATOR: u64 = 105;
const RANK_ALLOWANCE_DENOMINATOR: u64 = 100;
const RHO_S: f64 = 1.3;
const RANK_MODULUS: u64 = (1u64 << 61) - 1;

#[derive(Parser)]
#[command(about = "Screen Dickson-local additive rank on the registered P-256 base")]
struct Cli {
    #[arg(long)]
    round2: PathBuf,
    #[arg(long)]
    round19: PathBuf,
    #[arg(long)]
    round24: PathBuf,
    #[arg(long)]
    round311: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    assessment: PathBuf,
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    bytes: u64,
    sha256: String,
    schema: String,
}

#[derive(Clone, Copy)]
struct Column {
    point: Projective,
}

#[derive(Clone, Copy)]
struct TempAtom {
    point: Projective,
    left: u32,
    right: u32,
    right_sign: i8,
    level: u8,
    complete: bool,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct AtomRecord {
    x: [u8; 32],
    left: u32,
    right: u32,
    left_sign: i8,
    right_sign: i8,
    level: u8,
    complete: bool,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct MultipleRecord {
    x: [u8; 32],
    column: u32,
    coefficient: i32,
    multiplier: u32,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct RelationKey(Vec<(u32, i32)>);

#[derive(Clone)]
struct RelationEvidence {
    source: String,
    complete: bool,
}

#[derive(Default)]
struct RelationSet {
    rows: BTreeMap<RelationKey, RelationEvidence>,
    trivial_rows: u64,
}

#[derive(Clone, Default, Serialize)]
struct OperationCounts {
    dickson_squares: u64,
    projective_additions: u64,
    projective_doublings: u64,
    batch_inversions: u64,
    normalization_multiplications: u64,
    relation_replay_additions: u64,
    relation_replay_doublings: u64,
    hierarchy_pairs_examined: u64,
    atom_records_materialized: u64,
    small_multiple_records_materialized: u64,
    hash_pair_draws: u64,
}

#[derive(Clone, Serialize)]
struct LevelCensus {
    level: usize,
    complete: bool,
    parent_buckets: u64,
    nonempty_two_child_buckets: u64,
    maximum_left_width: usize,
    maximum_right_width: usize,
    unordered_pairs: u64,
    atom_records: u64,
    identity_atoms: u64,
    sample_draws: u64,
}

#[derive(Clone, Serialize)]
struct ToyControl {
    modulus: i64,
    columns: usize,
    atom_records: usize,
    multiple_records: usize,
    optimized_relations: usize,
    reference_relations: usize,
    independent_rank: usize,
    false_positives: usize,
    false_negatives: usize,
}

#[derive(Clone, Serialize)]
struct P256DetectorControl {
    scalar_columns: usize,
    atom_records: usize,
    multiple_records: usize,
    detected_relations: usize,
    replayed_relations: usize,
    replay_failures: usize,
    rank_mod_2_61_minus_1: usize,
    rank_mod_p256_order: usize,
    group_operations: OperationCounts,
}

#[derive(Clone, Serialize)]
struct RelationReceipt {
    terms: Vec<(u32, i32)>,
    source: String,
    complete_scope: bool,
    replayed: bool,
}

#[derive(Clone, Serialize)]
struct RankReceipt {
    complete_relation_rows: usize,
    sampled_relation_rows: usize,
    active_columns: usize,
    rank_mod_2_61_minus_1: Option<usize>,
    reverse_rank_mod_2_61_minus_1: Option<usize>,
    rank_mod_p256_order: Option<usize>,
    reverse_rank_mod_p256_order: Option<usize>,
    promotion_credited_rank: usize,
    exact_rank_materialization_limit_relations: usize,
    exact_rank_materialization_limit_active_columns: usize,
    rank_unset_due_to_limit: bool,
}

#[derive(Clone, Serialize)]
struct CostProjection {
    original_columns: usize,
    certified_static_rank: usize,
    effective_log_classes: usize,
    original_requested_dynamic_rows: u64,
    adjusted_dynamic_rows_with_same_five_percent_allowance: u64,
    free_oracle_s: f64,
    free_oracle_ratio_to_rho: f64,
    imported_unmodified_ratio_to_rho: f64,
    imported_s17_success_probability: f64,
    imported_per_usable_relation_fme_log2: f64,
    projected_collection_fme_log2: f64,
    collection_below_2_pow_120: bool,
    per_usable_relation_below_2_pow_103: bool,
    complete_cold_ratio_to_rho: Option<f64>,
    projection_is_measurement: bool,
    sampled_branches_receive_projection_credit: bool,
}

#[derive(Clone, Serialize)]
struct Gates {
    dependencies_and_factor_base_verified: bool,
    complete_toy_zero_false_positives_and_negatives: bool,
    p256_detector_controls_replayed: bool,
    every_reported_p256_relation_replayed: bool,
    structured_residual_degree_at_most_five: bool,
    exact_effective_log_rank_established: bool,
    projected_collection_below_2_pow_120: bool,
    projected_per_usable_relation_below_2_pow_103: bool,
    projected_storage_below_2_pow_50: bool,
    complete_cold_cost_at_or_below_rho: bool,
    sampled_branches_excluded_from_projection: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultDocument {
    schema: String,
    round: u32,
    curve: String,
    factor_base: Value,
    dependencies: Vec<Dependency>,
    screen: Value,
    toy_control: ToyControl,
    p256_detector_control: P256DetectorControl,
    levels: Vec<LevelCensus>,
    operation_counts: OperationCounts,
    atom_stream_blake3: String,
    complete_atom_stream_blake3: String,
    sampled_atom_stream_blake3: String,
    relation_stream_sha256: String,
    complete_relations: Vec<RelationReceipt>,
    sampled_relations: Vec<RelationReceipt>,
    rank: RankReceipt,
    structured_degree: Value,
    cost_projection: CostProjection,
    storage: Value,
    gates: Gates,
    decision: String,
    dominant_obstruction: String,
    full_depth_unplanted_p256_relation_attempted: bool,
    semantic_evidence_sha256: String,
}

#[derive(Serialize)]
struct TransferAssessment {
    schema: String,
    curve: String,
    screening_round: u32,
    skill_profile: String,
    required_companion_resources_available: bool,
    typed_correspondence_graph: Vec<String>,
    obligations: Vec<Value>,
    controls: Vec<String>,
    cost_accounting: Vec<String>,
    exploration_boundary: String,
    weakest_open_obligation: String,
    narrowest_supported_finding: String,
    semantic_evidence_sha256: String,
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("invalid non-negative integer `{value}`"))
}

fn checked_json(
    path: &Path,
    expected_sha256: &str,
    expected_schema: &str,
) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != expected_sha256 {
        return Err(format!(
            "{} hash {digest}, expected {expected_sha256}",
            path.display()
        ));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if value.get("schema").and_then(Value::as_str) != Some(expected_schema)
        || value.get("curve").and_then(Value::as_str) != Some(CURVE_SLUG)
    {
        return Err(format!("{} schema or curve mismatch", path.display()));
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            sha256: digest,
            schema: expected_schema.into(),
        },
        value,
    ))
}

fn load_dependencies(cli: &Cli) -> Result<(Vec<Dependency>, Value), String> {
    let (dep2, round2) = checked_json(&cli.round2, ROUND2_SHA256, "p256.dickson_coset_search/v1")?;
    let (dep19, round19) = checked_json(
        &cli.round19,
        ROUND19_SHA256,
        "p256.s17_multilevel_selector/v1",
    )?;
    let (dep24, round24) = checked_json(
        &cli.round24,
        ROUND24_SHA256,
        "p256.negation_cross_colour/v1",
    )?;
    let (dep311, round311) = checked_json(
        &cli.round311,
        ROUND311_SHA256,
        "p256-involution-atom-screen/v1",
    )?;
    let candidate = &round2["candidates"][4];
    if round2["selected_index"].as_u64() != Some(4)
        || candidate["fb_id"].as_str() != Some(FACTOR_BASE)
        || candidate["fb_sha256"].as_str() != Some(FACTOR_BASE_SHA256)
        || candidate["points_sha256"].as_str() != Some(POINTS_SHA256)
        || candidate["columns"].as_u64() != Some(COLUMNS as u64)
        || round19
            .pointer("/factor_base/fb_id")
            .and_then(Value::as_str)
            != Some(FACTOR_BASE)
        || round19
            .pointer("/residual_degree_evidence/structured_residual_maximum")
            .and_then(Value::as_u64)
            != Some(4)
        || round24["all_reported_relations_replayed"].as_bool() != Some(true)
        || round311
            .pointer("/gates/actual_nonlabelled_dimension_reduction")
            .and_then(Value::as_bool)
            != Some(false)
    {
        return Err("imported conclusion or identity mismatch".into());
    }
    Ok((vec![dep2, dep19, dep24, dep311], round19))
}

fn canonical_relation<I>(terms: I) -> Option<RelationKey>
where
    I: IntoIterator<Item = (u32, i32)>,
{
    let mut combined = BTreeMap::<u32, i64>::new();
    for (column, coefficient) in terms {
        *combined.entry(column).or_default() += i64::from(coefficient);
    }
    combined.retain(|_, coefficient| *coefficient != 0);
    if combined.is_empty() {
        return None;
    }
    let mut gcd = 0i64;
    for coefficient in combined.values() {
        gcd = gcd_i64(gcd, coefficient.abs());
    }
    let sign = if *combined.values().next().expect("nonempty") < 0 {
        -1i64
    } else {
        1i64
    };
    let row = combined
        .into_iter()
        .map(|(column, coefficient)| {
            let reduced = sign * coefficient / gcd;
            (
                column,
                i32::try_from(reduced).expect("small relation coefficient"),
            )
        })
        .collect();
    Some(RelationKey(row))
}

fn gcd_i64(mut left: i64, mut right: i64) -> i64 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left.abs()
}

impl RelationSet {
    fn insert<I>(&mut self, terms: I, source: String, complete: bool)
    where
        I: IntoIterator<Item = (u32, i32)>,
    {
        let Some(key) = canonical_relation(terms) else {
            self.trivial_rows += 1;
            return;
        };
        self.rows
            .entry(key)
            .and_modify(|existing| {
                if complete && !existing.complete {
                    existing.source = source.clone();
                }
                existing.complete |= complete;
            })
            .or_insert(RelationEvidence { source, complete });
    }
}

fn atom_terms(atom: &AtomRecord, scale: i32) -> [(u32, i32); 2] {
    [
        (atom.left, scale * i32::from(atom.left_sign)),
        (atom.right, scale * i32::from(atom.right_sign)),
    ]
}

fn multiple_terms(record: &MultipleRecord, scale: i32) -> [(u32, i32); 1] {
    [(record.column, scale * record.coefficient)]
}

fn low_orientation(y: &Fe) -> (i8, [u8; 32]) {
    let positive = y.to_bytes_be();
    let negative = y.neg().to_bytes_be();
    if positive <= negative {
        (1, positive)
    } else {
        (-1, negative)
    }
}

fn normalize_points(points: &[Projective], operations: &mut OperationCounts) -> Vec<(Fe, Fe)> {
    if points.is_empty() {
        return Vec::new();
    }
    assert!(points.iter().all(|point| !bool::from(point.is_identity())));
    let mut prefixes = Vec::with_capacity(points.len());
    let mut product = Fe::ONE;
    for point in points {
        prefixes.push(product);
        product = product.mul(&point.z);
    }
    let mut inverse = product.inv();
    operations.batch_inversions += 1;
    operations.normalization_multiplications += 5 * points.len() as u64;
    let mut output = vec![(Fe::ZERO, Fe::ZERO); points.len()];
    for index in (0..points.len()).rev() {
        let z_inverse = inverse.mul(&prefixes[index]);
        inverse = inverse.mul(&points[index].z);
        output[index] = (
            points[index].x.mul(&z_inverse),
            points[index].y.mul(&z_inverse),
        );
    }
    output
}

fn flush_atoms(
    pending: &mut Vec<TempAtom>,
    records: &mut Vec<AtomRecord>,
    all_hasher: &mut blake3::Hasher,
    complete_hasher: &mut blake3::Hasher,
    sampled_hasher: &mut blake3::Hasher,
    operations: &mut OperationCounts,
) {
    if pending.is_empty() {
        return;
    }
    let points: Vec<Projective> = pending.iter().map(|atom| atom.point).collect();
    let affine = normalize_points(&points, operations);
    for (atom, (x, y)) in pending.drain(..).zip(affine) {
        let (orientation, _) = low_orientation(&y);
        let record = AtomRecord {
            x: x.to_bytes_be(),
            left: atom.left,
            right: atom.right,
            left_sign: orientation,
            right_sign: orientation * atom.right_sign,
            level: atom.level,
            complete: atom.complete,
        };
        let encoded = encode_atom(&record);
        all_hasher.update(&encoded);
        if record.complete {
            complete_hasher.update(&encoded);
        } else {
            sampled_hasher.update(&encoded);
        }
        records.push(record);
    }
}

fn encode_atom(record: &AtomRecord) -> [u8; 44] {
    let mut encoded = [0u8; 44];
    encoded[..32].copy_from_slice(&record.x);
    encoded[32..36].copy_from_slice(&record.left.to_be_bytes());
    encoded[36..40].copy_from_slice(&record.right.to_be_bytes());
    encoded[40] = record.left_sign as u8;
    encoded[41] = record.right_sign as u8;
    encoded[42] = record.level;
    encoded[43] = u8::from(record.complete);
    encoded
}

fn load_factor_base(
    operations: &mut OperationCounts,
) -> Result<(Vec<Column>, Vec<Vec<[u8; 32]>>), String> {
    let built = p256_dickson_factor_base::build(SPEC)?;
    if built.dump.factor_base.fb_id != FACTOR_BASE
        || built.dump.factor_base.fb_sha256 != FACTOR_BASE_SHA256
        || built.dump.factor_base.points_sha256 != POINTS_SHA256
        || built.dump.factor_base.columns != COLUMNS as u64
        || built.dump.points.len() != 2 * COLUMNS
    {
        return Err("native factor-base rebuild identity mismatch".into());
    }
    let verified = p256_dickson_factor_base::verify(&built.dump)?;
    if verified.dump != built.dump {
        return Err("independent factor-base rebuild mismatch".into());
    }

    let curve = CurveParams::p256();
    let mut columns = Vec::with_capacity(COLUMNS);
    let mut keys = (0..=DEPTH)
        .map(|_| Vec::with_capacity(COLUMNS))
        .collect::<Vec<_>>();
    let two = Fe::from_canonical(&U256::from_biguint(&BigUint::from(2u8)));
    for pair in built.dump.points.chunks_exact(2) {
        let x = parse_big(&pair[0].x)?;
        let y = parse_big(&pair[0].y)?;
        let neg_x = parse_big(&pair[1].x)?;
        let neg_y = parse_big(&pair[1].y)?;
        if x != neg_x || (&y + &neg_y) % &curve.p != BigUint::zero() {
            return Err("factor-base signed pair mismatch".into());
        }
        let x_u = U256::from_biguint(&x);
        let y_u = U256::from_biguint(&y);
        let point = Projective::from_affine(&x_u, &y_u);
        let Some((roundtrip_x, roundtrip_y)) = point.to_affine() else {
            return Err("factor-base point converted to identity".into());
        };
        if roundtrip_x.0 != x_u.0 || roundtrip_y.0 != y_u.0 {
            return Err("factor-base point conversion mismatch".into());
        }
        columns.push(Column { point });
        let mut value = Fe::from_canonical(&x_u);
        keys[0].push(value.to_bytes_be());
        for level_keys in keys.iter_mut().take(DEPTH + 1).skip(1) {
            value = value.sqr().sub(&two);
            operations.dickson_squares += 1;
            level_keys.push(value.to_bytes_be());
        }
    }
    if columns.len() != COLUMNS
        || keys.iter().any(|level| level.len() != COLUMNS)
        || keys[DEPTH].windows(2).any(|pair| pair[0] != pair[1])
    {
        return Err("Dickson hierarchy reconstruction mismatch".into());
    }
    Ok((columns, keys))
}

fn push_atom(
    point: Projective,
    left: usize,
    right: usize,
    right_sign: i8,
    level: usize,
    complete: bool,
    pending: &mut Vec<TempAtom>,
    records: &mut Vec<AtomRecord>,
    relations: &mut RelationSet,
    all_hasher: &mut blake3::Hasher,
    complete_hasher: &mut blake3::Hasher,
    sampled_hasher: &mut blake3::Hasher,
    operations: &mut OperationCounts,
) -> u64 {
    if bool::from(point.is_identity()) {
        relations.insert(
            [(left as u32, 1), (right as u32, i32::from(right_sign))],
            format!("level-{level} identity atom"),
            complete,
        );
        return 1;
    }
    pending.push(TempAtom {
        point,
        left: left as u32,
        right: right as u32,
        right_sign,
        level: level as u8,
        complete,
    });
    if pending.len() >= NORMALIZE_CHUNK {
        flush_atoms(
            pending,
            records,
            all_hasher,
            complete_hasher,
            sampled_hasher,
            operations,
        );
    }
    0
}

fn enumerate_complete_levels(
    columns: &[Column],
    keys: &[Vec<[u8; 32]>],
    records: &mut Vec<AtomRecord>,
    relations: &mut RelationSet,
    all_hasher: &mut blake3::Hasher,
    complete_hasher: &mut blake3::Hasher,
    sampled_hasher: &mut blake3::Hasher,
    operations: &mut OperationCounts,
) -> Result<Vec<LevelCensus>, String> {
    let mut census = Vec::new();
    let mut pending = Vec::with_capacity(NORMALIZE_CHUNK);
    for level in 1..=COMPLETE_LEVELS {
        let mut entries = (0..columns.len())
            .map(|column| (keys[level][column], keys[level - 1][column], column))
            .collect::<Vec<_>>();
        entries.sort_unstable();
        let mut parent_buckets = 0u64;
        let mut two_child_buckets = 0u64;
        let mut maximum_left = 0usize;
        let mut maximum_right = 0usize;
        let mut unordered_pairs = 0u64;
        let mut identity_atoms = 0u64;
        let before = records.len();
        let mut start = 0usize;
        while start < entries.len() {
            let mut end = start + 1;
            while end < entries.len() && entries[end].0 == entries[start].0 {
                end += 1;
            }
            parent_buckets += 1;
            let mut child_ranges = Vec::new();
            let mut child_start = start;
            while child_start < end {
                let mut child_end = child_start + 1;
                while child_end < end && entries[child_end].1 == entries[child_start].1 {
                    child_end += 1;
                }
                child_ranges.push((child_start, child_end));
                child_start = child_end;
            }
            if child_ranges.len() > 2 {
                return Err(format!(
                    "level {level} parent has {} quadratic children",
                    child_ranges.len()
                ));
            }
            if child_ranges.len() == 2 {
                two_child_buckets += 1;
                let (left_start, left_end) = child_ranges[0];
                let (right_start, right_end) = child_ranges[1];
                maximum_left = maximum_left.max(left_end - left_start);
                maximum_right = maximum_right.max(right_end - right_start);
                for left_entry in &entries[left_start..left_end] {
                    for right_entry in &entries[right_start..right_end] {
                        let left = left_entry.2;
                        let right = right_entry.2;
                        unordered_pairs += 1;
                        operations.hierarchy_pairs_examined += 1;
                        operations.projective_additions += 2;
                        identity_atoms += push_atom(
                            columns[left].point.add(&columns[right].point),
                            left,
                            right,
                            1,
                            level,
                            true,
                            &mut pending,
                            records,
                            relations,
                            all_hasher,
                            complete_hasher,
                            sampled_hasher,
                            operations,
                        );
                        identity_atoms += push_atom(
                            columns[left].point.add(&columns[right].point.neg()),
                            left,
                            right,
                            -1,
                            level,
                            true,
                            &mut pending,
                            records,
                            relations,
                            all_hasher,
                            complete_hasher,
                            sampled_hasher,
                            operations,
                        );
                    }
                }
            }
            start = end;
        }
        flush_atoms(
            &mut pending,
            records,
            all_hasher,
            complete_hasher,
            sampled_hasher,
            operations,
        );
        census.push(LevelCensus {
            level,
            complete: true,
            parent_buckets,
            nonempty_two_child_buckets: two_child_buckets,
            maximum_left_width: maximum_left,
            maximum_right_width: maximum_right,
            unordered_pairs,
            atom_records: (records.len() - before) as u64,
            identity_atoms,
            sample_draws: 0,
        });
    }
    Ok(census)
}

fn minimum_coalescence(keys: &[Vec<[u8; 32]>], left: usize, right: usize) -> Option<usize> {
    (1..=DEPTH).find(|level| keys[*level][left] == keys[*level][right])
}

fn enumerate_deep_samples(
    columns: &[Column],
    keys: &[Vec<[u8; 32]>],
    records: &mut Vec<AtomRecord>,
    relations: &mut RelationSet,
    all_hasher: &mut blake3::Hasher,
    complete_hasher: &mut blake3::Hasher,
    sampled_hasher: &mut blake3::Hasher,
    operations: &mut OperationCounts,
) -> Result<Vec<LevelCensus>, String> {
    let mut census = Vec::new();
    let mut pending = Vec::with_capacity(NORMALIZE_CHUNK);
    for level in COMPLETE_LEVELS + 1..=DEPTH {
        let before = records.len();
        let mut selected = BTreeSet::<(usize, usize)>::new();
        let mut counter = 0u64;
        while selected.len() < DEEP_SAMPLES_PER_LEVEL {
            let mut message = b"round312/deep-pair/v1/".to_vec();
            message.extend_from_slice(&(level as u64).to_be_bytes());
            message.extend_from_slice(&counter.to_be_bytes());
            let digest = sha256(&message);
            counter += 1;
            let left_raw = u64::from_be_bytes(digest[..8].try_into().expect("eight bytes"));
            let right_raw = u64::from_be_bytes(digest[8..16].try_into().expect("eight bytes"));
            let mut left = left_raw as usize % columns.len();
            let mut right = right_raw as usize % columns.len();
            if left == right {
                continue;
            }
            if left > right {
                std::mem::swap(&mut left, &mut right);
            }
            if minimum_coalescence(keys, left, right) == Some(level) {
                selected.insert((left, right));
            }
            if counter > 200_000_000 {
                return Err(format!("level {level} deep sampler exceeded draw bound"));
            }
        }
        operations.hash_pair_draws += counter;
        let mut identity_atoms = 0u64;
        for (left, right) in selected {
            operations.hierarchy_pairs_examined += 1;
            operations.projective_additions += 2;
            identity_atoms += push_atom(
                columns[left].point.add(&columns[right].point),
                left,
                right,
                1,
                level,
                false,
                &mut pending,
                records,
                relations,
                all_hasher,
                complete_hasher,
                sampled_hasher,
                operations,
            );
            identity_atoms += push_atom(
                columns[left].point.add(&columns[right].point.neg()),
                left,
                right,
                -1,
                level,
                false,
                &mut pending,
                records,
                relations,
                all_hasher,
                complete_hasher,
                sampled_hasher,
                operations,
            );
        }
        flush_atoms(
            &mut pending,
            records,
            all_hasher,
            complete_hasher,
            sampled_hasher,
            operations,
        );
        census.push(LevelCensus {
            level,
            complete: false,
            parent_buckets: 0,
            nonempty_two_child_buckets: 0,
            maximum_left_width: 0,
            maximum_right_width: 0,
            unordered_pairs: DEEP_SAMPLES_PER_LEVEL as u64,
            atom_records: (records.len() - before) as u64,
            identity_atoms,
            sample_draws: counter,
        });
    }
    Ok(census)
}

fn flush_multiples(
    pending: &mut Vec<(Projective, u32, u32)>,
    records: &mut Vec<MultipleRecord>,
    hasher: &mut blake3::Hasher,
    operations: &mut OperationCounts,
) {
    if pending.is_empty() {
        return;
    }
    let points: Vec<Projective> = pending.iter().map(|entry| entry.0).collect();
    let affine = normalize_points(&points, operations);
    for ((_, column, multiplier), (x, y)) in pending.drain(..).zip(affine) {
        let (orientation, _) = low_orientation(&y);
        let record = MultipleRecord {
            x: x.to_bytes_be(),
            column,
            coefficient: i32::from(orientation) * multiplier as i32,
            multiplier,
        };
        let mut encoded = [0u8; 44];
        encoded[..32].copy_from_slice(&record.x);
        encoded[32..36].copy_from_slice(&record.column.to_be_bytes());
        encoded[36..40].copy_from_slice(&record.coefficient.to_be_bytes());
        encoded[40..44].copy_from_slice(&record.multiplier.to_be_bytes());
        hasher.update(&encoded);
        records.push(record);
    }
}

fn build_small_multiples(
    columns: &[Column],
    operations: &mut OperationCounts,
) -> (Vec<MultipleRecord>, String) {
    let mut records = Vec::with_capacity(columns.len() * MULTIPLIERS.len());
    let mut pending = Vec::with_capacity(NORMALIZE_CHUNK);
    let mut hasher = blake3::Hasher::new();
    for (column, base) in columns.iter().enumerate() {
        let p1 = base.point;
        let p2 = p1.double();
        let p3 = p2.add(&p1);
        let p4 = p2.double();
        let p8 = p4.double();
        let p16 = p8.double();
        operations.projective_doublings += 4;
        operations.projective_additions += 1;
        for (multiplier, point) in [(1, p1), (2, p2), (3, p3), (4, p4), (8, p8), (16, p16)] {
            pending.push((point, column as u32, multiplier));
            if pending.len() >= NORMALIZE_CHUNK {
                flush_multiples(&mut pending, &mut records, &mut hasher, operations);
            }
        }
    }
    flush_multiples(&mut pending, &mut records, &mut hasher, operations);
    (records, hasher.finalize().to_hex().to_string())
}

#[derive(Default, Serialize)]
struct CollisionCensus {
    atom_output_groups_with_multiplicity: u64,
    atom_collision_pairs: u64,
    complete_atom_relations_before_deduplication: u64,
    sampled_atom_relations_before_deduplication: u64,
    atom_to_multiple_matches: u64,
    complete_atom_to_multiple_matches: u64,
    sampled_atom_to_multiple_matches: u64,
    multiple_output_groups_with_multiplicity: u64,
    multiple_collision_pairs: u64,
}

fn discover_relations(
    atoms: &mut [AtomRecord],
    multiples: &mut [MultipleRecord],
    relations: &mut RelationSet,
) -> CollisionCensus {
    atoms.sort_unstable();
    multiples.sort_unstable();
    let mut census = CollisionCensus::default();

    let mut start = 0usize;
    while start < atoms.len() {
        let mut end = start + 1;
        while end < atoms.len() && atoms[end].x == atoms[start].x {
            end += 1;
        }
        if end - start > 1 {
            census.atom_output_groups_with_multiplicity += 1;
            census.atom_collision_pairs += ((end - start) * (end - start - 1) / 2) as u64;
            for left in start..end {
                for right in left + 1..end {
                    let mut terms = atom_terms(&atoms[left], 1).to_vec();
                    terms.extend(atom_terms(&atoms[right], -1));
                    let complete = atoms[left].complete && atoms[right].complete;
                    relations.insert(
                        terms,
                        format!(
                            "{} atom collision levels {} and {}",
                            if complete { "complete" } else { "sampled" },
                            atoms[left].level,
                            atoms[right].level
                        ),
                        complete,
                    );
                    if complete {
                        census.complete_atom_relations_before_deduplication += 1;
                    } else {
                        census.sampled_atom_relations_before_deduplication += 1;
                    }
                }
            }
        }
        start = end;
    }

    start = 0;
    while start < multiples.len() {
        let mut end = start + 1;
        while end < multiples.len() && multiples[end].x == multiples[start].x {
            end += 1;
        }
        if end - start > 1 {
            census.multiple_output_groups_with_multiplicity += 1;
            census.multiple_collision_pairs += ((end - start) * (end - start - 1) / 2) as u64;
            for left in start..end {
                for right in left + 1..end {
                    let mut terms = multiple_terms(&multiples[left], 1).to_vec();
                    terms.extend(multiple_terms(&multiples[right], -1));
                    relations.insert(
                        terms,
                        format!(
                            "small-multiple collision [{}]P_{} and [{}]P_{}",
                            multiples[left].multiplier,
                            multiples[left].column,
                            multiples[right].multiplier,
                            multiples[right].column
                        ),
                        true,
                    );
                }
            }
        }
        start = end;
    }

    let mut atom_index = 0usize;
    let mut multiple_index = 0usize;
    while atom_index < atoms.len() && multiple_index < multiples.len() {
        match atoms[atom_index].x.cmp(&multiples[multiple_index].x) {
            std::cmp::Ordering::Less => atom_index += 1,
            std::cmp::Ordering::Greater => multiple_index += 1,
            std::cmp::Ordering::Equal => {
                let atom_start = atom_index;
                let multiple_start = multiple_index;
                while atom_index < atoms.len() && atoms[atom_index].x == atoms[atom_start].x {
                    atom_index += 1;
                }
                while multiple_index < multiples.len()
                    && multiples[multiple_index].x == multiples[multiple_start].x
                {
                    multiple_index += 1;
                }
                for atom in &atoms[atom_start..atom_index] {
                    for multiple in &multiples[multiple_start..multiple_index] {
                        let mut terms = atom_terms(atom, 1).to_vec();
                        terms.extend(multiple_terms(multiple, -1));
                        relations.insert(
                            terms,
                            format!(
                                "{} atom/multiple match at level {} against [{}]P_{}",
                                if atom.complete { "complete" } else { "sampled" },
                                atom.level,
                                multiple.multiplier,
                                multiple.column
                            ),
                            atom.complete,
                        );
                        census.atom_to_multiple_matches += 1;
                        if atom.complete {
                            census.complete_atom_to_multiple_matches += 1;
                        } else {
                            census.sampled_atom_to_multiple_matches += 1;
                        }
                    }
                }
            }
        }
    }
    census
}

fn coefficient_mod(coefficient: i32, modulus: &BigUint) -> BigUint {
    if coefficient >= 0 {
        BigUint::from(coefficient as u32) % modulus
    } else {
        let magnitude = BigUint::from(coefficient.unsigned_abs()) % modulus;
        if magnitude.is_zero() {
            BigUint::zero()
        } else {
            modulus - magnitude
        }
    }
}

fn sparse_rank(rows: &[RelationKey], modulus: &BigUint, reverse: bool) -> usize {
    let mut pivots = BTreeMap::<u32, BTreeMap<u32, BigUint>>::new();
    let order: Box<dyn Iterator<Item = &RelationKey>> = if reverse {
        Box::new(rows.iter().rev())
    } else {
        Box::new(rows.iter())
    };
    for relation in order {
        let mut row = relation
            .0
            .iter()
            .map(|(column, coefficient)| (*column, coefficient_mod(*coefficient, modulus)))
            .filter(|(_, coefficient)| !coefficient.is_zero())
            .collect::<BTreeMap<_, _>>();
        while let Some((&pivot, coefficient)) = row.first_key_value() {
            let coefficient = coefficient.clone();
            if let Some(base) = pivots.get(&pivot).cloned() {
                for (column, base_coefficient) in base {
                    let product = (&coefficient * base_coefficient) % modulus;
                    let current = row.remove(&column).unwrap_or_else(BigUint::zero);
                    let reduced = if current >= product {
                        current - product
                    } else {
                        current + modulus - product
                    };
                    if !reduced.is_zero() {
                        row.insert(column, reduced);
                    }
                }
            } else {
                let inverse = coefficient.modpow(&(modulus - BigUint::from(2u8)), modulus);
                for value in row.values_mut() {
                    *value = (&*value * &inverse) % modulus;
                }
                pivots.insert(pivot, row);
                break;
            }
        }
    }
    pivots.len()
}

fn rank_rows(rows: &[RelationKey], modulus: &BigUint) -> (usize, usize) {
    (
        sparse_rank(rows, modulus, false),
        sparse_rank(rows, modulus, true),
    )
}

fn small_scalar_mul(
    point: &Projective,
    mut scalar: u32,
    operations: &mut OperationCounts,
) -> Projective {
    let mut result = Projective::IDENTITY;
    let mut base = *point;
    while scalar != 0 {
        if scalar & 1 == 1 {
            result = result.add(&base);
            operations.relation_replay_additions += 1;
        }
        scalar >>= 1;
        if scalar != 0 {
            base = base.double();
            operations.relation_replay_doublings += 1;
        }
    }
    result
}

fn replay_relation(
    relation: &RelationKey,
    columns: &[Column],
    operations: &mut OperationCounts,
) -> bool {
    let mut sum = Projective::IDENTITY;
    for (column, coefficient) in &relation.0 {
        let mut term = small_scalar_mul(
            &columns[*column as usize].point,
            coefficient.unsigned_abs(),
            operations,
        );
        if *coefficient < 0 {
            term = term.neg();
        }
        sum = sum.add(&term);
        operations.relation_replay_additions += 1;
    }
    bool::from(sum.is_identity())
}

fn relation_receipts(
    relations: &RelationSet,
    columns: &[Column],
    operations: &mut OperationCounts,
) -> Result<(Vec<RelationReceipt>, Vec<RelationReceipt>, String), String> {
    let mut complete = Vec::new();
    let mut sampled = Vec::new();
    let mut bytes = Vec::new();
    for (relation, evidence) in &relations.rows {
        let replayed = replay_relation(relation, columns, operations);
        if !replayed {
            return Err(format!("relation replay failed: {:?}", relation.0));
        }
        bytes.push(u8::from(evidence.complete));
        for (column, coefficient) in &relation.0 {
            bytes.extend_from_slice(&column.to_be_bytes());
            bytes.extend_from_slice(&coefficient.to_be_bytes());
        }
        bytes.extend_from_slice(&u32::MAX.to_be_bytes());
        let receipt = RelationReceipt {
            terms: relation.0.clone(),
            source: evidence.source.clone(),
            complete_scope: evidence.complete,
            replayed,
        };
        if evidence.complete {
            complete.push(receipt);
        } else {
            sampled.push(receipt);
        }
    }
    Ok((complete, sampled, sha256_hex(&bytes)))
}

fn rank_receipt(relations: &RelationSet) -> RankReceipt {
    let complete_rows = relations
        .rows
        .iter()
        .filter(|(_, evidence)| evidence.complete)
        .map(|(row, _)| row.clone())
        .collect::<Vec<_>>();
    let sampled_rows = relations.rows.len() - complete_rows.len();
    let active = complete_rows
        .iter()
        .flat_map(|row| row.0.iter().map(|term| term.0))
        .collect::<BTreeSet<_>>()
        .len();
    let limited = complete_rows.len() > RANK_MAX_RELATIONS || active > RANK_MAX_ACTIVE_COLUMNS;
    let (rank61, reverse61, rankn, reversen) = if limited {
        (None, None, None, None)
    } else {
        let modulus61 = BigUint::from(RANK_MODULUS);
        let modulusn = CurveParams::p256().n;
        let (rank61, reverse61) = rank_rows(&complete_rows, &modulus61);
        let (rankn, reversen) = rank_rows(&complete_rows, &modulusn);
        (Some(rank61), Some(reverse61), Some(rankn), Some(reversen))
    };
    let credit = match (rank61, reverse61, rankn, reversen) {
        (Some(a), Some(b), Some(c), Some(d)) if a == b && c == d => a.min(c),
        _ => 0,
    };
    RankReceipt {
        complete_relation_rows: complete_rows.len(),
        sampled_relation_rows: sampled_rows,
        active_columns: active,
        rank_mod_2_61_minus_1: rank61,
        reverse_rank_mod_2_61_minus_1: reverse61,
        rank_mod_p256_order: rankn,
        reverse_rank_mod_p256_order: reversen,
        promotion_credited_rank: credit,
        exact_rank_materialization_limit_relations: RANK_MAX_RELATIONS,
        exact_rank_materialization_limit_active_columns: RANK_MAX_ACTIVE_COLUMNS,
        rank_unset_due_to_limit: limited,
    }
}

fn toy_key(value: i64, modulus: i64) -> (i8, [u8; 32]) {
    let value = value.rem_euclid(modulus);
    let negative = (-value).rem_euclid(modulus);
    let (orientation, canonical) = if value <= negative {
        (1, value)
    } else {
        (-1, negative)
    };
    let mut key = [0u8; 32];
    key[24..].copy_from_slice(&(canonical as u64).to_be_bytes());
    (orientation, key)
}

fn toy_control() -> ToyControl {
    const MODULUS: i64 = 101;
    let values = [1i64, 2, 4, 6, 7, 9, 12, 20];
    let mut atoms = Vec::new();
    for left in 0..values.len() {
        for right in left + 1..values.len() {
            for right_sign in [1i8, -1i8] {
                let value = values[left] + i64::from(right_sign) * values[right];
                let (orientation, x) = toy_key(value, MODULUS);
                atoms.push(AtomRecord {
                    x,
                    left: left as u32,
                    right: right as u32,
                    left_sign: orientation,
                    right_sign: orientation * right_sign,
                    level: 3,
                    complete: true,
                });
            }
        }
    }
    let mut multiples = Vec::new();
    for (column, value) in values.iter().enumerate() {
        for multiplier in MULTIPLIERS {
            let (orientation, x) = toy_key(i64::from(multiplier) * value, MODULUS);
            multiples.push(MultipleRecord {
                x,
                column: column as u32,
                coefficient: i32::from(orientation) * multiplier as i32,
                multiplier,
            });
        }
    }

    let mut optimized = RelationSet::default();
    discover_relations(&mut atoms, &mut multiples, &mut optimized);

    enum Rep<'a> {
        Atom(&'a AtomRecord),
        Multiple(&'a MultipleRecord),
    }
    let representatives = atoms
        .iter()
        .map(Rep::Atom)
        .chain(multiples.iter().map(Rep::Multiple))
        .collect::<Vec<_>>();
    let mut reference = BTreeSet::new();
    for left in 0..representatives.len() {
        for right in left + 1..representatives.len() {
            let left_key = match representatives[left] {
                Rep::Atom(record) => record.x,
                Rep::Multiple(record) => record.x,
            };
            let right_key = match representatives[right] {
                Rep::Atom(record) => record.x,
                Rep::Multiple(record) => record.x,
            };
            if left_key != right_key {
                continue;
            }
            let mut terms = match representatives[left] {
                Rep::Atom(record) => atom_terms(record, 1).to_vec(),
                Rep::Multiple(record) => multiple_terms(record, 1).to_vec(),
            };
            match representatives[right] {
                Rep::Atom(record) => terms.extend(atom_terms(record, -1)),
                Rep::Multiple(record) => terms.extend(multiple_terms(record, -1)),
            }
            if let Some(row) = canonical_relation(terms) {
                reference.insert(row);
            }
        }
    }
    let optimized_keys = optimized.rows.keys().cloned().collect::<BTreeSet<_>>();
    let false_positives = optimized_keys.difference(&reference).count();
    let false_negatives = reference.difference(&optimized_keys).count();
    let rows = optimized_keys.iter().cloned().collect::<Vec<_>>();
    let (rank, reverse) = rank_rows(&rows, &BigUint::from(MODULUS as u64));
    assert_eq!(rank, reverse);
    ToyControl {
        modulus: MODULUS,
        columns: values.len(),
        atom_records: atoms.len(),
        multiple_records: multiples.len(),
        optimized_relations: optimized_keys.len(),
        reference_relations: reference.len(),
        independent_rank: rank,
        false_positives,
        false_negatives,
    }
}

fn p256_detector_control(operations: &mut OperationCounts) -> Result<P256DetectorControl, String> {
    let scalars = [1u32, 2, 4, 6, 7, 9, 12, 20];
    let curve = CurveParams::p256();
    let generator = Projective::from_affine(
        &U256::from_biguint(&curve.gx),
        &U256::from_biguint(&curve.gy),
    );
    operations.projective_additions += scalars.len() as u64 * curve.order_bits() as u64;
    operations.projective_doublings += scalars.len() as u64 * 2 * curve.order_bits() as u64;
    let columns = scalars
        .iter()
        .map(|scalar| Column {
            point: generator.scalar_mul_ct(&BigUint::from(*scalar), curve.order_bits()),
        })
        .collect::<Vec<_>>();
    let mut atoms = Vec::new();
    let mut pending = Vec::new();
    let mut relations = RelationSet::default();
    let mut all_hasher = blake3::Hasher::new();
    let mut complete_hasher = blake3::Hasher::new();
    let mut sampled_hasher = blake3::Hasher::new();
    for left in 0..columns.len() {
        for right in left + 1..columns.len() {
            operations.projective_additions += 2;
            push_atom(
                columns[left].point.add(&columns[right].point),
                left,
                right,
                1,
                3,
                true,
                &mut pending,
                &mut atoms,
                &mut relations,
                &mut all_hasher,
                &mut complete_hasher,
                &mut sampled_hasher,
                operations,
            );
            push_atom(
                columns[left].point.add(&columns[right].point.neg()),
                left,
                right,
                -1,
                3,
                true,
                &mut pending,
                &mut atoms,
                &mut relations,
                &mut all_hasher,
                &mut complete_hasher,
                &mut sampled_hasher,
                operations,
            );
        }
    }
    flush_atoms(
        &mut pending,
        &mut atoms,
        &mut all_hasher,
        &mut complete_hasher,
        &mut sampled_hasher,
        operations,
    );
    let (mut multiples, _) = build_small_multiples(&columns, operations);
    discover_relations(&mut atoms, &mut multiples, &mut relations);
    let rows = relations.rows.keys().cloned().collect::<Vec<_>>();
    let modulus61 = BigUint::from(RANK_MODULUS);
    let (rank61, reverse61) = rank_rows(&rows, &modulus61);
    let (rankn, reversen) = rank_rows(&rows, &curve.n);
    if rank61 != reverse61 || rankn != reversen {
        return Err("P-256 detector control rank order mismatch".into());
    }
    let mut replayed = 0usize;
    let mut failures = 0usize;
    for row in &rows {
        replayed += 1;
        failures += usize::from(!replay_relation(row, &columns, operations));
    }
    Ok(P256DetectorControl {
        scalar_columns: columns.len(),
        atom_records: atoms.len(),
        multiple_records: multiples.len(),
        detected_relations: rows.len(),
        replayed_relations: replayed,
        replay_failures: failures,
        rank_mod_2_61_minus_1: rank61,
        rank_mod_p256_order: rankn,
        group_operations: operations.clone(),
    })
}

fn gamma_ratio(collisions: u64) -> f64 {
    assert!(collisions >= 1);
    let mut ratio = std::f64::consts::PI.sqrt() / 2.0;
    for k in 1..collisions {
        ratio *= (k as f64 + 0.5) / k as f64;
    }
    ratio
}

fn cost_projection(rank: usize, round19: &Value) -> Result<CostProjection, String> {
    if rank >= COLUMNS {
        return Err("homogeneous static rank cannot consume every log column".into());
    }
    let effective = COLUMNS - rank;
    let adjusted_rows =
        (RANK_ALLOWANCE_NUMERATOR * effective as u64).div_ceil(RANK_ALLOWANCE_DENOMINATOR);
    let free_s = 2.0f64.sqrt() * gamma_ratio(effective as u64);
    let imported_ratio = 394.425_280_345_252_2;
    let success = round19
        .pointer("/full_projection/whole_factor_base_success")
        .and_then(Value::as_f64)
        .ok_or("Round 19 relation success missing")?;
    let per_relation = round19
        .pointer("/full_projection/optimistic_per_usable_relation_fme_log2")
        .and_then(Value::as_f64)
        .ok_or("Round 19 per-relation cost missing")?;
    let collection_log2 = per_relation + (adjusted_rows as f64).log2();
    Ok(CostProjection {
        original_columns: COLUMNS,
        certified_static_rank: rank,
        effective_log_classes: effective,
        original_requested_dynamic_rows: REQUESTED_ROWS,
        adjusted_dynamic_rows_with_same_five_percent_allowance: adjusted_rows,
        free_oracle_s: free_s,
        free_oracle_ratio_to_rho: free_s / RHO_S,
        imported_unmodified_ratio_to_rho: imported_ratio,
        imported_s17_success_probability: success,
        imported_per_usable_relation_fme_log2: per_relation,
        projected_collection_fme_log2: collection_log2,
        collection_below_2_pow_120: collection_log2 < 120.0,
        per_usable_relation_below_2_pow_103: per_relation < 103.0,
        complete_cold_ratio_to_rho: None,
        projection_is_measurement: false,
        sampled_branches_receive_projection_credit: false,
    })
}

fn obligation(name: &str, status: &str, evidence: String, scope: &str) -> Value {
    json!({
        "name": name,
        "status": status,
        "evidence": evidence,
        "scope": scope,
    })
}

fn build_assessment(
    rank: &RankReceipt,
    complete_relations: &[RelationReceipt],
    sampled_relations: &[RelationReceipt],
    cost: &CostProjection,
    replay_failures: usize,
) -> TransferAssessment {
    let rank_status = if rank.rank_unset_due_to_limit {
        "unknown"
    } else if rank.promotion_credited_rank == 0 {
        "refuted in complete scope"
    } else {
        "supported in complete scope"
    };
    let narrowest = if rank.promotion_credited_rank == 0 {
        format!(
            "The complete Dickson-local levels-1--{COMPLETE_LEVELS} and small-multiple screen found no nontrivial static log relation; the registered base retains all {COLUMNS} independent classes within this scope."
        )
    } else {
        format!(
            "The complete screen certifies static rank {}, leaving {} log classes, but the executable S17 relation cost remains above the registered gates.",
            rank.promotion_credited_rank, cost.effective_log_classes
        )
    };
    let semantic = json!({
        "complete_rows": complete_relations.len(),
        "sampled_rows": sampled_relations.len(),
        "rank": rank.promotion_credited_rank,
        "effective": cost.effective_log_classes,
        "replay_failures": replay_failures,
    });
    TransferAssessment {
        schema: "p256-dickson-additive-rank-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 312,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        typed_correspondence_graph: vec![
            "FB1h2f8621cda105 columns --[Dickson D_2 hierarchy]--> exact local sibling domains".into(),
            "local column pairs --[P_i +/- P_j]--> canonical P-256 atom outputs".into(),
            "equal atom or small-multiple outputs --[subtract sparse rows]--> homogeneous factor-base log relation".into(),
            "static relation matrix --[dual-modulus exact rank]--> reduced independent-log quotient".into(),
            "reduced quotient --[S17 relation collection and recovery]--> one-target DLP only if every complete cost is charged".into(),
        ],
        obligations: vec![
            obligation(
                "source and subgroup identity",
                "supported",
                format!("Native rebuild fixed {COLUMNS} columns on the registered prime-order P-256 subgroup."),
                "Exact registered factor base and ICV1 curve.",
            ),
            obligation(
                "Dickson-local non-scalar log transport",
                rank_status,
                format!(
                    "{} complete-scope relation rows, certified rank {}, {} replay failures.",
                    complete_relations.len(), rank.promotion_credited_rank, replay_failures
                ),
                &format!("All minimal coalescence pairs through level {COMPLETE_LEVELS}, with six complete small-multiple tables."),
            ),
            obligation(
                "deeper hierarchy",
                "unknown",
                format!("{} sampled relation rows receive no rank or projection credit.", sampled_relations.len()),
                "Exactly 4,096 hash-selected pairs per level 8--18; not exhaustive.",
            ),
            obligation(
                "end-to-end advantage",
                "refuted for the measured candidate",
                format!(
                    "Executable per-usable-relation projection is 2^{:.6} FME; complete cold S is unset.",
                    cost.imported_per_usable_relation_fme_log2
                ),
                "Unchanged Round-19 S17 selector with certified static-rank credit only.",
            ),
        ],
        controls: vec![
            "Complete order-101 scalar reference compares every atom and multiplier representation pair.".into(),
            "Deterministic scalar-labelled P-256 control plants and replays low-weight identities.".into(),
            "Every actual reported relation is independently summed to the P-256 identity.".into(),
            "Ranks are recomputed in both row orders over 2^61-1 and the P-256 subgroup order.".into(),
        ],
        cost_accounting: vec![
            "Measured: hierarchy candidates, group operations, normalization field operations, materialized records, isolation telemetry, and artifact bytes.".into(),
            "Imported projection: Round-19 S17 success, per-relation FME, sparse algebra, and storage; not reclassified as measurement.".into(),
            "Deep sampled pairs are excluded from rank and all projected improvements.".into(),
            "Complete cold DLP cost remains unset because no qualifying relation/recovery pipeline completed.".into(),
        ],
        exploration_boundary: format!(
            "Exact for all pair atoms with minimal Dickson coalescence at levels 1--{COMPLETE_LEVELS}, all atom collisions in that set, atom matches to [1,2,3,4,8,16] times every base column, and all collisions within those multiple tables. Levels 8--18 are sampled only."
        ),
        weakest_open_obligation: "Supply a complete non-scalar static relation family outside the screened local domain, or a non-generic collector that lowers the unchanged 2^148.080980-FME usable-relation cost while retaining exact recovery.".into(),
        narrowest_supported_finding: narrowest,
        semantic_evidence_sha256: sha256_hex(
            &serde_json::to_vec(&semantic).expect("assessment semantic JSON"),
        ),
    }
}

fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<u64, String> {
    let text = serde_json::to_string_pretty(value).map_err(|error| error.to_string())? + "\n";
    fs::write(path, text.as_bytes()).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok(text.len() as u64)
}

fn run(cli: Cli) -> Result<(), String> {
    let (dependencies, round19) = load_dependencies(&cli)?;
    let toy = toy_control();
    if toy.false_positives != 0 || toy.false_negatives != 0 || toy.independent_rank == 0 {
        return Err("complete toy detector control failed".into());
    }

    let mut control_operations = OperationCounts::default();
    let p256_control = p256_detector_control(&mut control_operations)?;
    if p256_control.detected_relations == 0 || p256_control.replay_failures != 0 {
        return Err("P-256 detector control failed".into());
    }

    let mut operations = OperationCounts::default();
    let (columns, keys) = load_factor_base(&mut operations)?;
    let mut atoms = Vec::<AtomRecord>::new();
    let mut relations = RelationSet::default();
    let mut all_hasher = blake3::Hasher::new();
    let mut complete_hasher = blake3::Hasher::new();
    let mut sampled_hasher = blake3::Hasher::new();
    let mut levels = enumerate_complete_levels(
        &columns,
        &keys,
        &mut atoms,
        &mut relations,
        &mut all_hasher,
        &mut complete_hasher,
        &mut sampled_hasher,
        &mut operations,
    )?;
    levels.extend(enumerate_deep_samples(
        &columns,
        &keys,
        &mut atoms,
        &mut relations,
        &mut all_hasher,
        &mut complete_hasher,
        &mut sampled_hasher,
        &mut operations,
    )?);
    operations.atom_records_materialized = atoms.len() as u64;
    let atom_stream_blake3 = all_hasher.finalize().to_hex().to_string();
    let complete_atom_stream_blake3 = complete_hasher.finalize().to_hex().to_string();
    let sampled_atom_stream_blake3 = sampled_hasher.finalize().to_hex().to_string();

    let (mut multiples, multiple_stream_blake3) = build_small_multiples(&columns, &mut operations);
    operations.small_multiple_records_materialized = multiples.len() as u64;
    let collision_census = discover_relations(&mut atoms, &mut multiples, &mut relations);
    let rank = rank_receipt(&relations);
    let (complete_relations, sampled_relations, relation_stream_sha256) =
        relation_receipts(&relations, &columns, &mut operations)?;
    if complete_relations.len() != rank.complete_relation_rows
        || sampled_relations.len() != rank.sampled_relation_rows
    {
        return Err("relation partition mismatch".into());
    }
    let replay_failures = complete_relations
        .iter()
        .chain(&sampled_relations)
        .filter(|receipt| !receipt.replayed)
        .count();
    let cost = cost_projection(rank.promotion_credited_rank, &round19)?;

    let atom_bytes = atoms.len() as u64 * size_of::<AtomRecord>() as u64;
    let multiple_bytes = multiples.len() as u64 * size_of::<MultipleRecord>() as u64;
    let hierarchy_key_bytes = (DEPTH as u64 + 1) * COLUMNS as u64 * 32;
    let input_point_bytes = COLUMNS as u64 * size_of::<Column>() as u64;
    let screen_peak_logical_bytes =
        atom_bytes + multiple_bytes + hierarchy_key_bytes + input_point_bytes;
    let imported_attack_storage_log2 = round19
        .pointer("/full_projection/structured_frontier_materialized_bytes_log2")
        .and_then(Value::as_f64)
        .ok_or("Round 19 storage projection missing")?;

    let exact_rank = !rank.rank_unset_due_to_limit
        && rank.rank_mod_2_61_minus_1 == rank.reverse_rank_mod_2_61_minus_1
        && rank.rank_mod_p256_order == rank.reverse_rank_mod_p256_order;
    let projected_storage_below = imported_attack_storage_log2 < 50.0;
    let gates = Gates {
        dependencies_and_factor_base_verified: true,
        complete_toy_zero_false_positives_and_negatives: true,
        p256_detector_controls_replayed: true,
        every_reported_p256_relation_replayed: replay_failures == 0,
        structured_residual_degree_at_most_five: true,
        exact_effective_log_rank_established: exact_rank,
        projected_collection_below_2_pow_120: cost.collection_below_2_pow_120,
        projected_per_usable_relation_below_2_pow_103: cost.per_usable_relation_below_2_pow_103,
        projected_storage_below_2_pow_50: projected_storage_below,
        complete_cold_cost_at_or_below_rho: cost
            .complete_cold_ratio_to_rho
            .is_some_and(|ratio| ratio <= 1.0),
        sampled_branches_excluded_from_projection: !cost.sampled_branches_receive_projection_credit,
        promoted: false,
    };
    let promoted = gates.dependencies_and_factor_base_verified
        && gates.complete_toy_zero_false_positives_and_negatives
        && gates.p256_detector_controls_replayed
        && gates.every_reported_p256_relation_replayed
        && gates.structured_residual_degree_at_most_five
        && gates.exact_effective_log_rank_established
        && gates.projected_collection_below_2_pow_120
        && gates.projected_per_usable_relation_below_2_pow_103
        && gates.projected_storage_below_2_pow_50
        && gates.complete_cold_cost_at_or_below_rho
        && gates.sampled_branches_excluded_from_projection;
    let mut gates = gates;
    gates.promoted = promoted;

    let decision = if rank.promotion_credited_rank == 0 {
        format!(
            "Negative in the preregistered complete domain: no Dickson-local or small-multiple identity reduces the {COLUMNS}-column log rank. The free-oracle boundary remains {:.6} times rho, while the executable S17 per-relation projection remains 2^{:.6} FME.",
            cost.free_oracle_ratio_to_rho, cost.imported_per_usable_relation_fme_log2
        )
    } else {
        format!(
            "Static rank {} is certified, leaving {} log classes, but the executable per-relation, collection, storage, and complete-cold gates do not all pass.",
            rank.promotion_credited_rank, cost.effective_log_classes
        )
    };
    let obstruction = if rank.promotion_credited_rank == 0 {
        "The explicit Dickson hierarchy does not induce a low-weight P-256 group-law identity in the complete screened domain; coordinate coalescence is not logarithm transport. The unchanged S17 selector still materializes a 2^143.568654 right frontier.".to_string()
    } else {
        "Certified static identities do not lower the unchanged S17 decomposition cost enough to meet the per-row or complete collection gates.".to_string()
    };
    let semantic = json!({
        "complete_levels": COMPLETE_LEVELS,
        "deep_samples_per_level": DEEP_SAMPLES_PER_LEVEL,
        "levels": &levels,
        "atom_stream_blake3": atom_stream_blake3,
        "multiple_stream_blake3": multiple_stream_blake3,
        "relations_sha256": relation_stream_sha256,
        "complete_relations": complete_relations.len(),
        "sampled_relations": sampled_relations.len(),
        "rank": rank.promotion_credited_rank,
        "effective_log_classes": cost.effective_log_classes,
        "free_oracle_ratio_to_rho": cost.free_oracle_ratio_to_rho,
        "promoted": promoted,
    });
    let semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).map_err(|error| error.to_string())?);
    let storage = json!({
        "atom_record_size_bytes": size_of::<AtomRecord>(),
        "multiple_record_size_bytes": size_of::<MultipleRecord>(),
        "atom_records_logical_bytes": atom_bytes,
        "small_multiple_records_logical_bytes": multiple_bytes,
        "hierarchy_keys_logical_bytes": hierarchy_key_bytes,
        "input_points_logical_bytes": input_point_bytes,
        "screen_peak_logical_bytes_upper_sum": screen_peak_logical_bytes,
        "screen_peak_logical_bytes_log2": (screen_peak_logical_bytes as f64).log2(),
        "imported_attack_materialized_bytes_log2": imported_attack_storage_log2,
        "imported_attack_storage_below_2_pow_50": projected_storage_below,
        "result_json_bytes": 0,
        "assessment_json_bytes": 0,
        "disk_traffic_is_artifact_writes_only": true,
    });
    let mut document = ResultDocument {
        schema: "p256-dickson-additive-rank-screen/v1".into(),
        round: 312,
        curve: CURVE_SLUG.into(),
        factor_base: json!({
            "fb_id": FACTOR_BASE,
            "fb_sha256": FACTOR_BASE_SHA256,
            "points_sha256": POINTS_SHA256,
            "spec": SPEC,
            "columns": COLUMNS,
            "signed_points": 2 * COLUMNS,
            "complete_native_rebuild_verified": true,
        }),
        dependencies,
        screen: json!({
            "complete_levels": format!("1..={COMPLETE_LEVELS}"),
            "complete_domain": "all unordered pairs whose minimal Dickson coalescence is at most 7; both P_i+P_j and P_i-P_j",
            "deep_sample_levels": format!("{}..={DEPTH}", COMPLETE_LEVELS + 1),
            "deep_samples_per_level": DEEP_SAMPLES_PER_LEVEL,
            "deep_sample_domain": "round312/deep-pair/v1/<level>/<counter>",
            "small_multipliers": MULTIPLIERS,
            "collision_census": collision_census,
            "trivial_rows_removed": relations.trivial_rows,
            "multiple_stream_blake3": multiple_stream_blake3,
        }),
        toy_control: toy,
        p256_detector_control: p256_control,
        levels,
        operation_counts: operations,
        atom_stream_blake3,
        complete_atom_stream_blake3,
        sampled_atom_stream_blake3,
        relation_stream_sha256,
        complete_relations,
        sampled_relations,
        rank,
        structured_degree: json!({
            "factor_base_unchanged": true,
            "round19_local_join_degree": 2,
            "round19_structured_residual_maximum": 4,
            "structured_residual_degree_gate_at_most_five": true,
            "unsplit_s17_degree_of_regularity": null,
        }),
        cost_projection: cost,
        storage,
        gates,
        decision,
        dominant_obstruction: obstruction,
        full_depth_unplanted_p256_relation_attempted: false,
        semantic_evidence_sha256,
    };
    let assessment = build_assessment(
        &document.rank,
        &document.complete_relations,
        &document.sampled_relations,
        &document.cost_projection,
        replay_failures,
    );

    let assessment_bytes = write_json(&cli.assessment, &assessment)?;
    if let Some(storage) = document.storage.as_object_mut() {
        storage.insert("assessment_json_bytes".into(), assessment_bytes.into());
    }
    let mut result_bytes = 0u64;
    for _ in 0..4 {
        if let Some(storage) = document.storage.as_object_mut() {
            storage.insert("result_json_bytes".into(), result_bytes.into());
        }
        let text =
            serde_json::to_string_pretty(&document).map_err(|error| error.to_string())? + "\n";
        let next = text.len() as u64;
        fs::write(&cli.out, text.as_bytes())
            .map_err(|error| format!("{}: {error}", cli.out.display()))?;
        if next == result_bytes {
            break;
        }
        result_bytes = next;
    }
    eprintln!(
        "round 312: atoms={}, complete_relations={}, sampled_relations={}, rank={}, effective={}, promoted={}",
        document.operation_counts.atom_records_materialized,
        document.complete_relations.len(),
        document.sampled_relations.len(),
        document.rank.promotion_credited_rank,
        document.cost_projection.effective_log_classes,
        document.gates.promoted,
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_dickson_additive_rank: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relation_canonicalization_combines_gcd_and_sign() {
        assert_eq!(
            canonical_relation([(3, -2), (1, 4), (3, -2)]),
            Some(RelationKey(vec![(1, 1), (3, -1)]))
        );
        assert_eq!(canonical_relation([(1, 1), (1, -1)]), None);
    }

    #[test]
    fn complete_toy_detector_matches_direct_reference() {
        let control = toy_control();
        assert_eq!(control.false_positives, 0);
        assert_eq!(control.false_negatives, 0);
        assert!(control.optimized_relations > 0);
        assert!(control.independent_rank > 0);
    }

    #[test]
    fn batch_normalization_matches_individual_affine_conversion() {
        let curve = CurveParams::p256();
        let generator = Projective::from_affine(
            &U256::from_biguint(&curve.gx),
            &U256::from_biguint(&curve.gy),
        );
        let points = vec![
            generator,
            generator.double(),
            generator.double().add(&generator),
        ];
        let mut operations = OperationCounts::default();
        let normalized = normalize_points(&points, &mut operations);
        for (point, (x, y)) in points.iter().zip(normalized) {
            let (want_x, want_y) = point.to_affine().unwrap();
            assert_eq!(x.to_canonical().0, want_x.0);
            assert_eq!(y.to_canonical().0, want_y.0);
        }
        assert_eq!(operations.batch_inversions, 1);
    }

    #[test]
    fn collision_boundary_matches_registered_values() {
        let ratio1 = 2.0f64.sqrt() * gamma_ratio(1) / RHO_S;
        let ratio2 = 2.0f64.sqrt() * gamma_ratio(2) / RHO_S;
        assert!((ratio1 - 0.964_087_797_935_000_2).abs() < 1e-12);
        assert!((ratio2 - 1.446_131_696_902_500_4).abs() < 1e-12);
    }
}
