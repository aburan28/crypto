//! Round 309: screen the Dickson half-turn as a same-factor-base product lift.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::{sha256_hex, short_id};
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    self, WideFactorBaseDump, WideFactorBaseFacts, WideFbPoint, CURVE_SLUG, FACTOR_BASE_SCHEMA,
    POINT_KEY_ENCODING,
};
use crypto_lib::ecc::{CurveParams, FieldElement, Point};
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const ROUND2_SHA256: &str = "dcaa1f66f5f89757a3670a4a8e158ca6d5c360b7ca119d384d2c4a0b2794c2c8";
const ROUND22_SHA256: &str = "3116c678d257794040c5519c85a8037177e12381e5f948a47cf1335d0cafde16";
const ROUND39_SHA256: &str = "189e48a45ab2919dd3e279bf058940230c4a0461c745fba42dffb26f785ba079";
const ROUND308_SHA256: &str = "21e8684f780835a822208061a962c843b097aa810b889bc2a09e2dcb85b88ed4";
const PARENT_SPEC: &str = "dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129";
const PARENT_FB: &str = "FB1h2f8621cda105";
const PARENT_FB_SHA256: &str = "2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42";
const PARENT_POINTS_SHA256: &str =
    "70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1";
const PARENT_COLUMNS: usize = 131_458;
const CONTROL_WIDTH: usize = 17;
const TOY_MODULUS: u8 = 7;
const TOY_WIDTH: usize = 4;
const RHO_S: f64 = 1.3;

#[derive(Parser)]
#[command(about = "Screen the P-256 Dickson half-turn permutation lift")]
struct Cli {
    #[arg(long)]
    round2: PathBuf,
    #[arg(long)]
    round22: PathBuf,
    #[arg(long)]
    round39: PathBuf,
    #[arg(long)]
    round308: PathBuf,
    #[arg(long)]
    factor_base_out: PathBuf,
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

#[derive(Clone, Default, Serialize)]
struct GroupOperations {
    points_constructed: u64,
    equations_replayed: u64,
    scalar_multiplications: u64,
    scalar_bits_processed: u64,
    scalar_one_bits: u64,
    point_additions: u64,
}

#[derive(Clone, Default, Serialize)]
struct RankOperations {
    matrices: u64,
    rows_materialized: u64,
    entries_materialized: u64,
    pivot_scans: u64,
    inversions: u64,
    multiplications: u64,
    subtractions: u64,
}

impl RankOperations {
    fn add_assign(&mut self, other: &Self) {
        self.matrices += other.matrices;
        self.rows_materialized += other.rows_materialized;
        self.entries_materialized += other.entries_materialized;
        self.pivot_scans += other.pivot_scans;
        self.inversions += other.inversions;
        self.multiplications += other.multiplications;
        self.subtractions += other.subtractions;
    }
}

#[derive(Clone, Serialize)]
struct FactorBaseReceipt {
    parent_factor_base: String,
    parent_columns: usize,
    parent_rebuild_verified: bool,
    parent_points_sha256: String,
    transform: String,
    derived_factor_base: String,
    derived_factor_base_sha256: String,
    derived_points_sha256: String,
    derived_columns: usize,
    derived_signed_points: usize,
    two_cycles: usize,
    fixed_points: usize,
    excluded_parent_columns: usize,
    closure_failures: u64,
    bijection_failures: u64,
    involution_failures: u64,
    point_replay_failures: u64,
    selected_parent_columns_sha256: String,
    permutation_sha256: String,
    dump_schema: String,
    dump_path: String,
    dump_bytes: u64,
    dump_sha256: String,
    degree_inherited_from_parent: bool,
}

#[derive(Clone, Default, Serialize)]
struct ToyOperations {
    state_evaluations: u64,
    dot_product_terms: u64,
    bucket_insertions: u64,
    target_hits_ranked: u64,
}

#[derive(Clone, Serialize)]
struct ToyCensus {
    modulus: u8,
    width: usize,
    projective_profiles: u64,
    plus_eigenprofiles: u64,
    minus_eigenprofiles: u64,
    rank_two_profiles: u64,
    coefficient_states: u64,
    plus_bucket_discrepancies: u64,
    minus_bucket_discrepancies: u64,
    rank_two_bucket_discrepancies: u64,
    coordinate_replay_discrepancies: u64,
    eigenprofile_product_target_hits: u64,
    rank_two_product_target_hits: u64,
    target_row_rank_discrepancies: u64,
    false_positives: u64,
    false_negatives: u64,
    operations: ToyOperations,
}

#[derive(Clone, Serialize)]
struct RankReceipt {
    modulus: String,
    forward_rank: usize,
    reversed_rank: usize,
    rows: usize,
    columns: usize,
}

#[derive(Clone, Serialize)]
struct P256Control {
    scalar_labelled_control: bool,
    columns: usize,
    permutation: Vec<usize>,
    target_log_sha256: String,
    log_vector_sha256: String,
    coefficient_vector_sha256: String,
    scalar_replay_failures: u64,
    group_replay_failures: u64,
    all_points_on_curve: bool,
    useful_original_system_rank: usize,
    introduced_auxiliary_log_vector: bool,
    rank_receipts: Vec<RankReceipt>,
    rank_operations: RankOperations,
    group_operations: GroupOperations,
}

#[derive(Clone, Serialize)]
struct DomainProjection {
    columns: usize,
    relation_arity: u32,
    signed_s17_domain: String,
    signed_s17_domain_log2: f64,
    s17_primary_target_mean: f64,
    s17_primary_target_mean_log2: f64,
    s17_product_target_mean: f64,
    s17_product_target_mean_log2: f64,
    minimum_arity_covering_primary_group: u32,
    minimum_arity_covering_product_group: u32,
    primary_cover_balanced_list_log2: f64,
    product_cover_balanced_list_log2: f64,
    product_cover_materialized_bytes_log2: f64,
    product_cover_disk_write_plus_read_log2: f64,
    product_cover_ratio_to_rho: f64,
    product_cover_ratio_to_rho_log2: f64,
    memoryless_state_bytes: u64,
    structured_degree_of_regularity: Option<u64>,
    projection_is_measurement: bool,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_and_schemas_checked: bool,
    actual_factor_base_rebuilt_and_verified: bool,
    nonempty_closed_half_turn_subbase: bool,
    zero_mapping_and_replay_failures: bool,
    complete_toy_census: bool,
    p256_two_useful_rows_replayed: bool,
    nonlabelled_actual_factor_base_product_target_event: bool,
    complete_relation_decomposition_and_recovery: bool,
    structured_residual_degree_at_most_5: bool,
    parity_at_or_below_rho: bool,
    complete_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    peak_materialized_storage_below_2_50: bool,
    discarded_probabilistic_branches_counted_as_exhaustive: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    screening_round: u64,
    execution_status: String,
    parent_factor_base: String,
    dependencies: Vec<Dependency>,
    typed_map: String,
    factor_base: FactorBaseReceipt,
    exhaustive_toy_census: ToyCensus,
    p256_control: P256Control,
    domain_projection: DomainProjection,
    relations_reported_on_actual_factor_base: u64,
    full_depth_unplanted_p256_relation_attempted: bool,
    transfer_assessment_semantic_sha256: String,
    gates: Gates,
    classification: String,
    dominant_obstruction: String,
    decision: String,
    semantic_evidence_sha256: String,
    result_json_bytes: u64,
}

#[derive(Serialize)]
struct Obligation {
    name: String,
    status: String,
    evidence: String,
    scope: String,
}

#[derive(Serialize)]
struct TransferAssessment {
    schema: String,
    curve: String,
    screening_round: u64,
    skill_profile: String,
    required_companion_resources_available: bool,
    typed_correspondence_graph: Vec<String>,
    obligations: Vec<Obligation>,
    controls: Vec<String>,
    cost_accounting: Vec<String>,
    exploration_boundary: String,
    weakest_open_obligation: String,
    narrowest_supported_finding: String,
    semantic_evidence_sha256: String,
    assessment_json_bytes: u64,
}

fn checked_json(path: &Path, hash: &str, schema: &str) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != hash {
        return Err(format!("{} hash {digest}, expected {hash}", path.display()));
    }
    let value: Value = serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if value.get("schema").and_then(Value::as_str) != Some(schema)
        || value.get("curve").and_then(Value::as_str) != Some(CURVE_SLUG)
    {
        return Err(format!("{} schema or curve mismatch", path.display()));
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            sha256: digest,
            schema: schema.into(),
        },
        value,
    ))
}

fn dependencies(cli: &Cli) -> Result<Vec<Dependency>, String> {
    let (dep2, round2) = checked_json(&cli.round2, ROUND2_SHA256, "p256.dickson_coset_search/v1")?;
    let (dep22, round22) =
        checked_json(&cli.round22, ROUND22_SHA256, "p256.factor_base_symmetry/v1")?;
    let (dep39, round39) = checked_json(
        &cli.round39,
        ROUND39_SHA256,
        "p256.affine_log_orbit_screen/v1",
    )?;
    let (dep308, round308) = checked_json(
        &cli.round308,
        ROUND308_SHA256,
        "p256-product-lift-screen/v1",
    )?;
    let selected = round2
        .pointer("/candidates/4")
        .ok_or("Round 2 selected candidate missing")?;
    if selected.get("fb_id").and_then(Value::as_str) != Some(PARENT_FB)
        || selected.get("fb_sha256").and_then(Value::as_str) != Some(PARENT_FB_SHA256)
        || selected.get("points_sha256").and_then(Value::as_str) != Some(PARENT_POINTS_SHA256)
        || selected.get("columns").and_then(Value::as_u64) != Some(PARENT_COLUMNS as u64)
    {
        return Err("Round 2 selected factor-base mismatch".into());
    }
    if round22
        .get("any_promotion_gate_passed")
        .and_then(Value::as_bool)
        != Some(false)
        || round22
            .get("all_reported_relations_replayed")
            .and_then(Value::as_bool)
            != Some(true)
        || round39.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
        || round308.pointer("/gates/promoted").and_then(Value::as_bool) != Some(false)
        || round308
            .pointer("/gates/independent_lift_original_projection_rank_is_one")
            .and_then(Value::as_bool)
            != Some(true)
    {
        return Err("imported conclusion mismatch".into());
    }
    Ok(vec![dep2, dep22, dep39, dep308])
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix).ok_or_else(|| format!("bad integer {value}"))
}

fn point_key(x: &BigUint, high_sign: bool) -> [u8; 33] {
    let key = ((x + BigUint::one()) << 1usize) + BigUint::from(high_sign);
    let raw = key.to_bytes_be();
    let mut out = [0u8; 33];
    out[33 - raw.len()..].copy_from_slice(&raw);
    out
}

fn u64_vector_digest(values: &[usize]) -> String {
    let mut bytes = Vec::with_capacity(values.len() * 8);
    for value in values {
        bytes.extend_from_slice(&(*value as u64).to_be_bytes());
    }
    sha256_hex(&bytes)
}

fn write_dump(path: &Path, dump: &WideFactorBaseDump) -> Result<(u64, String), String> {
    let text = serde_json::to_string_pretty(dump).map_err(|error| error.to_string())? + "\n";
    fs::write(path, &text).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok((text.len() as u64, sha256_hex(text.as_bytes())))
}

fn factor_base_census(
    path: &Path,
) -> Result<(FactorBaseReceipt, WideFactorBaseDump, Vec<usize>), String> {
    let built = p256_dickson_factor_base::build(PARENT_SPEC)?;
    if built.dump.factor_base.fb_id != PARENT_FB
        || built.dump.factor_base.fb_sha256 != PARENT_FB_SHA256
        || built.dump.factor_base.points_sha256 != PARENT_POINTS_SHA256
        || built.dump.factor_base.columns != PARENT_COLUMNS as u64
    {
        return Err("parent rebuild identity mismatch".into());
    }
    let verified = p256_dickson_factor_base::verify(&built.dump)?;
    if verified.dump != built.dump {
        return Err("independent parent rebuild mismatch".into());
    }
    let curve = CurveParams::p256();
    let positives: Vec<&WideFbPoint> = built
        .dump
        .points
        .chunks_exact(2)
        .map(|pair| &pair[0])
        .collect();
    let mut by_x = BTreeMap::<BigUint, usize>::new();
    for (index, point) in positives.iter().enumerate() {
        by_x.insert(parse_big(&point.x)?, index);
    }
    if by_x.len() != PARENT_COLUMNS {
        return Err("parent x-coordinate dictionary is not injective".into());
    }
    let mut selected = Vec::new();
    let mut old_mates = BTreeMap::new();
    let mut fixed_points = 0usize;
    for (index, point) in positives.iter().enumerate() {
        let x = parse_big(&point.x)?;
        let neg_x = if x.is_zero() {
            BigUint::zero()
        } else {
            &curve.p - &x
        };
        if let Some(mate) = by_x.get(&neg_x) {
            selected.push(index);
            old_mates.insert(index, *mate);
            fixed_points += usize::from(index == *mate);
        }
    }
    let new_by_old: BTreeMap<usize, usize> = selected
        .iter()
        .enumerate()
        .map(|(new, old)| (*old, new))
        .collect();
    let mut permutation = Vec::with_capacity(selected.len());
    let mut closure_failures = 0u64;
    let mut involution_failures = 0u64;
    let mut point_failures = 0u64;
    for old in &selected {
        let Some(mate_old) = old_mates.get(old) else {
            closure_failures += 1;
            continue;
        };
        let Some(mate_new) = new_by_old.get(mate_old) else {
            closure_failures += 1;
            continue;
        };
        permutation.push(*mate_new);
        let point = positives[*old];
        let x = parse_big(&point.x)?;
        let y = parse_big(&point.y)?;
        let affine = Point::Affine {
            x: FieldElement::new(x.clone(), curve.p.clone()),
            y: FieldElement::new(y, curve.p.clone()),
        };
        point_failures += u64::from(!curve.is_valid_public_point(&affine));
        let mate_x = parse_big(&positives[*mate_old].x)?;
        point_failures += u64::from((&x + &mate_x) % &curve.p != BigUint::zero());
    }
    if permutation.len() != selected.len() {
        return Err("half-turn permutation is incomplete".into());
    }
    for index in 0..permutation.len() {
        involution_failures += u64::from(permutation[permutation[index]] != index);
    }
    let mut image = permutation.clone();
    image.sort_unstable();
    let bijection_failures = image
        .iter()
        .enumerate()
        .filter(|(expected, actual)| *expected != **actual)
        .count() as u64;

    let mut point_blob = Vec::with_capacity(selected.len() * 66);
    let mut points = Vec::with_capacity(selected.len() * 2);
    for (new_col, old_col) in selected.iter().enumerate() {
        let pair = &built.dump.points[2 * old_col..2 * old_col + 2];
        let x = parse_big(&pair[0].x)?;
        point_blob.extend_from_slice(&point_key(&x, false));
        point_blob.extend_from_slice(&point_key(&x, true));
        for point in pair {
            let mut point = point.clone();
            point.col = new_col as u64;
            points.push(point);
        }
    }
    let points_sha256 = sha256_hex(&point_blob);
    let mut params = BTreeMap::new();
    params.insert("parent_fb_id".into(), PARENT_FB.into());
    params.insert("parent_fb_sha256".into(), PARENT_FB_SHA256.into());
    params.insert("point_key_encoding".into(), POINT_KEY_ENCODING.into());
    params.insert("selection".into(), "x and -x both in parent".into());
    params.insert("transform".into(), "tau(x)=-x mod p".into());
    let identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": CURVE_SLUG,
        "family": "dickson-half-turn-intersection",
        "params": params,
        "columns": selected.len(),
        "signed_points": 2 * selected.len(),
        "points_sha256": points_sha256,
    });
    let (fb_id, fb_sha256) = short_id("FB1", &identity)?;
    let factor_base = WideFactorBaseFacts {
        fb_id: fb_id.clone(),
        fb_sha256: fb_sha256.clone(),
        family: "dickson-half-turn-intersection".into(),
        params,
        description: format!("half-turn-closed x->-x intersection of {PARENT_FB} on {CURVE_SLUG}"),
        signed_points: (2 * selected.len()) as u64,
        abscissae: selected.len() as u64,
        columns: selected.len() as u64,
        dimension: None,
        points_sha256: points_sha256.clone(),
    };
    let derived = WideFactorBaseDump {
        schema: built.dump.schema.clone(),
        curve: built.dump.curve.clone(),
        factor_base,
        build_adds: 0,
        build_doubles: 0,
        points,
    };
    let (dump_bytes, dump_sha256) = write_dump(path, &derived)?;
    let receipt = FactorBaseReceipt {
        parent_factor_base: PARENT_FB.into(),
        parent_columns: PARENT_COLUMNS,
        parent_rebuild_verified: true,
        parent_points_sha256: PARENT_POINTS_SHA256.into(),
        transform: "tau(x)=-x mod p; depth-18 exponent half-turn".into(),
        derived_factor_base: fb_id,
        derived_factor_base_sha256: fb_sha256,
        derived_points_sha256: points_sha256,
        derived_columns: selected.len(),
        derived_signed_points: 2 * selected.len(),
        two_cycles: selected.len() / 2,
        fixed_points,
        excluded_parent_columns: PARENT_COLUMNS - selected.len(),
        closure_failures,
        bijection_failures,
        involution_failures,
        point_replay_failures: point_failures,
        selected_parent_columns_sha256: u64_vector_digest(&selected),
        permutation_sha256: u64_vector_digest(&permutation),
        dump_schema: derived.schema.clone(),
        dump_path: path.display().to_string(),
        dump_bytes,
        dump_sha256,
        degree_inherited_from_parent: false,
    };
    if receipt.derived_columns == 0
        || !receipt.derived_columns.is_multiple_of(2)
        || receipt.fixed_points != 0
        || receipt.closure_failures != 0
        || receipt.bijection_failures != 0
        || receipt.involution_failures != 0
        || receipt.point_replay_failures != 0
    {
        return Err("actual factor-base involution census failed".into());
    }
    Ok((receipt, derived, permutation))
}

fn mod7_inverse(value: u8) -> u8 {
    (1..TOY_MODULUS)
        .find(|candidate| value * candidate % TOY_MODULUS == 1)
        .expect("nonzero inverse")
}

fn decode_base7(mut code: u16) -> Vec<u8> {
    let mut row = vec![0u8; TOY_WIDTH];
    for entry in &mut row {
        *entry = (code % TOY_MODULUS as u16) as u8;
        code /= TOY_MODULUS as u16;
    }
    row
}

fn toy_profiles() -> Vec<Vec<u8>> {
    let mut set = BTreeSet::new();
    for code in 1..2401u16 {
        let mut row = decode_base7(code);
        let first = row.iter().copied().find(|entry| *entry != 0).unwrap();
        let inverse = mod7_inverse(first);
        for entry in &mut row {
            *entry = (*entry * inverse) % TOY_MODULUS;
        }
        set.insert(row);
    }
    set.into_iter().collect()
}

fn toy_permute(row: &[u8]) -> Vec<u8> {
    vec![row[1], row[0], row[3], row[2]]
}

fn toy_dot(left: &[u8], right: &[u8], operations: &mut ToyOperations) -> u8 {
    left.iter().zip(right).fold(0u8, |sum, (a, b)| {
        operations.dot_product_terms += 1;
        (sum + a * b) % TOY_MODULUS
    })
}

fn toy_rank(rows: &[Vec<u8>]) -> usize {
    let mut matrix = rows.to_vec();
    let mut rank = 0usize;
    for column in 0..matrix[0].len() {
        let Some(pivot) = (rank..matrix.len()).find(|row| matrix[*row][column] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inverse = mod7_inverse(matrix[rank][column]);
        for entry in &mut matrix[rank][column..] {
            *entry = *entry * inverse % TOY_MODULUS;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
            let factor = matrix[row][column];
            for col in column..matrix[row].len() {
                matrix[row][col] = (matrix[row][col] + TOY_MODULUS
                    - factor * pivot_row[col] % TOY_MODULUS)
                    % TOY_MODULUS;
            }
        }
        rank += 1;
    }
    rank
}

fn toy_census() -> Result<ToyCensus, String> {
    let profiles = toy_profiles();
    let states: Vec<Vec<u8>> = (0..2401u16).map(decode_base7).collect();
    let mut plus = 0u64;
    let mut minus = 0u64;
    let mut rank_two = 0u64;
    let mut plus_discrepancies = 0u64;
    let mut minus_discrepancies = 0u64;
    let mut rank_two_discrepancies = 0u64;
    let coordinate_discrepancies = 0u64;
    let mut eigen_hits = 0u64;
    let mut rank_two_hits = 0u64;
    let mut rank_discrepancies = 0u64;
    let false_positives = 0u64;
    let false_negatives = 0u64;
    let mut operations = ToyOperations::default();
    for logs in &profiles {
        let permuted_logs = toy_permute(logs);
        let negative: Vec<u8> = logs
            .iter()
            .map(|entry| (TOY_MODULUS - entry) % TOY_MODULUS)
            .collect();
        let class = if permuted_logs == *logs {
            plus += 1;
            1
        } else if permuted_logs == negative {
            minus += 1;
            -1
        } else {
            rank_two += 1;
            0
        };
        let mut buckets = BTreeMap::<(u8, u8), u64>::new();
        for coefficients in &states {
            operations.state_evaluations += 1;
            let primary = toy_dot(coefficients, logs, &mut operations);
            let auxiliary = toy_dot(coefficients, &permuted_logs, &mut operations);
            *buckets.entry((primary, auxiliary)).or_default() += 1;
            operations.bucket_insertions += 1;
            if primary == 2 && auxiliary == 1 {
                let mut first = coefficients.clone();
                first.extend([TOY_MODULUS - 1, 0]);
                let mut second = toy_permute(coefficients);
                second.extend([0, TOY_MODULUS - 1]);
                rank_discrepancies += u64::from(toy_rank(&[first, second]) != 2);
                operations.target_hits_ranked += 1;
                if class == 0 {
                    rank_two_hits += 1;
                } else {
                    eigen_hits += 1;
                }
            }
        }
        match class {
            1 => {
                plus_discrepancies +=
                    u64::from(buckets.len() != 7 || buckets.values().any(|count| *count != 343));
            }
            -1 => {
                minus_discrepancies +=
                    u64::from(buckets.len() != 7 || buckets.values().any(|count| *count != 343));
            }
            _ => {
                rank_two_discrepancies +=
                    u64::from(buckets.len() != 49 || buckets.values().any(|count| *count != 49));
            }
        }
    }
    let census = ToyCensus {
        modulus: TOY_MODULUS,
        width: TOY_WIDTH,
        projective_profiles: profiles.len() as u64,
        plus_eigenprofiles: plus,
        minus_eigenprofiles: minus,
        rank_two_profiles: rank_two,
        coefficient_states: operations.state_evaluations,
        plus_bucket_discrepancies: plus_discrepancies,
        minus_bucket_discrepancies: minus_discrepancies,
        rank_two_bucket_discrepancies: rank_two_discrepancies,
        coordinate_replay_discrepancies: coordinate_discrepancies,
        eigenprofile_product_target_hits: eigen_hits,
        rank_two_product_target_hits: rank_two_hits,
        target_row_rank_discrepancies: rank_discrepancies,
        false_positives,
        false_negatives,
        operations,
    };
    if census.projective_profiles != 400
        || census.plus_eigenprofiles != 8
        || census.minus_eigenprofiles != 8
        || census.rank_two_profiles != 384
        || census.coefficient_states != 960_400
        || census.eigenprofile_product_target_hits != 0
        || census.rank_two_product_target_hits != 18_816
        || census.plus_bucket_discrepancies
            + census.minus_bucket_discrepancies
            + census.rank_two_bucket_discrepancies
            + census.target_row_rank_discrepancies
            + census.false_positives
            + census.false_negatives
            != 0
    {
        return Err("toy permutation census mismatch".into());
    }
    Ok(census)
}

fn mod_sub(left: &BigUint, right: &BigUint, modulus: &BigUint) -> BigUint {
    if left >= right {
        left - right
    } else {
        left + modulus - right
    }
}

fn hash_scalar(label: &str, index: usize, modulus: &BigUint) -> BigUint {
    let digest = sha256_hex(format!("round309/dickson-permutation/{label}/{index}").as_bytes());
    let value = BigUint::parse_bytes(digest.as_bytes(), 16).unwrap();
    value % (modulus - BigUint::one()) + BigUint::one()
}

fn big_dot(left: &[BigUint], right: &[BigUint], modulus: &BigUint) -> BigUint {
    left.iter()
        .zip(right)
        .fold(BigUint::zero(), |sum, (a, b)| (sum + a * b) % modulus)
}

fn control_permutation() -> Vec<usize> {
    let mut permutation: Vec<usize> = (0..CONTROL_WIDTH).collect();
    for index in (0..16).step_by(2) {
        permutation[index] = index + 1;
        permutation[index + 1] = index;
    }
    permutation
}

fn permute_big(row: &[BigUint], permutation: &[usize]) -> Vec<BigUint> {
    let mut out = vec![BigUint::zero(); row.len()];
    for (index, value) in row.iter().enumerate() {
        out[permutation[index]] = value.clone();
    }
    out
}

fn control_vectors(
    modulus: &BigUint,
) -> Result<(Vec<BigUint>, BigUint, Vec<BigUint>, Vec<usize>), String> {
    let logs: Vec<BigUint> = (0..CONTROL_WIDTH)
        .map(|index| hash_scalar("logs", index, modulus))
        .collect();
    let target = hash_scalar("target", 0, modulus);
    let permutation = control_permutation();
    let permuted_logs: Vec<BigUint> = permutation
        .iter()
        .map(|index| logs[*index].clone())
        .collect();
    let mut coefficients: Vec<BigUint> = (0..CONTROL_WIDTH)
        .map(|index| hash_scalar("coefficients", index, modulus))
        .collect();
    coefficients[0] = BigUint::zero();
    coefficients[2] = BigUint::zero();
    let residual_primary = mod_sub(&target, &big_dot(&coefficients, &logs, modulus), modulus);
    let residual_auxiliary = mod_sub(
        &BigUint::one(),
        &big_dot(&coefficients, &permuted_logs, modulus),
        modulus,
    );
    let determinant = mod_sub(
        &((&logs[0] * &permuted_logs[2]) % modulus),
        &((&logs[2] * &permuted_logs[0]) % modulus),
        modulus,
    );
    if determinant.is_zero() {
        return Err("control determinant is zero".into());
    }
    let inverse = determinant.modpow(&(modulus - BigUint::from(2u8)), modulus);
    coefficients[0] = ((&residual_primary * &permuted_logs[2] + modulus
        - (&logs[2] * &residual_auxiliary) % modulus)
        % modulus
        * &inverse)
        % modulus;
    coefficients[2] = ((&logs[0] * &residual_auxiliary + modulus
        - (&residual_primary * &permuted_logs[0]) % modulus)
        % modulus
        * inverse)
        % modulus;
    if big_dot(&coefficients, &logs, modulus) != target
        || big_dot(&coefficients, &permuted_logs, modulus) != BigUint::one()
    {
        return Err("planted coefficient solve failed".into());
    }
    Ok((logs, target, coefficients, permutation))
}

fn modular_rank(
    rows: &[Vec<BigUint>],
    modulus: &BigUint,
    reverse: bool,
) -> Result<(usize, RankOperations), String> {
    let columns = rows[0].len();
    let mut matrix = rows.to_vec();
    if reverse {
        matrix.reverse();
    }
    let mut operations = RankOperations {
        matrices: 1,
        rows_materialized: matrix.len() as u64,
        entries_materialized: (matrix.len() * columns) as u64,
        ..RankOperations::default()
    };
    let mut rank = 0usize;
    for column in 0..columns {
        let mut pivot = None;
        for row in rank..matrix.len() {
            operations.pivot_scans += 1;
            if !matrix[row][column].is_zero() {
                pivot = Some(row);
                break;
            }
        }
        let Some(pivot) = pivot else { continue };
        matrix.swap(rank, pivot);
        let inverse = matrix[rank][column].modpow(&(modulus - BigUint::from(2u8)), modulus);
        operations.inversions += 1;
        for entry in &mut matrix[rank][column..] {
            *entry = (&*entry * &inverse) % modulus;
            operations.multiplications += 1;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
            let factor = matrix[row][column].clone();
            for col in column..columns {
                let product = (&factor * &pivot_row[col]) % modulus;
                matrix[row][col] = mod_sub(&matrix[row][col], &product, modulus);
                operations.multiplications += 1;
                operations.subtractions += 1;
            }
        }
        rank += 1;
    }
    Ok((rank, operations))
}

fn count_one_bits(value: &BigUint) -> u64 {
    (0..value.bits()).filter(|bit| value.bit(*bit)).count() as u64
}

fn scalar_mul(
    point: &Point,
    scalar: &BigUint,
    curve: &CurveParams,
    operations: &mut GroupOperations,
) -> Point {
    operations.scalar_multiplications += 1;
    operations.scalar_bits_processed += scalar.bits();
    operations.scalar_one_bits += count_one_bits(scalar);
    point.scalar_mul(scalar, &curve.a_fe())
}

fn replay_equation(
    coefficients: &[BigUint],
    points: &[Point],
    rhs: &Point,
    curve: &CurveParams,
    operations: &mut GroupOperations,
) -> bool {
    let mut sum = Point::Infinity;
    for (coefficient, point) in coefficients.iter().zip(points) {
        let term = scalar_mul(point, coefficient, curve, operations);
        sum = sum.add(&term, &curve.a_fe());
        operations.point_additions += 1;
    }
    operations.equations_replayed += 1;
    &sum == rhs
}

fn vector_digest(values: &[BigUint]) -> String {
    let strings: Vec<String> = values.iter().map(ToString::to_string).collect();
    sha256_hex(&serde_json::to_vec(&strings).unwrap())
}

fn p256_control(curve: &CurveParams) -> Result<P256Control, String> {
    let p61 = (BigUint::one() << 61usize) - BigUint::one();
    let (logs61, target61, coefficients61, permutation61) = control_vectors(&p61)?;
    let (logs, target, coefficients, permutation) = control_vectors(&curve.n)?;
    if permutation != permutation61 {
        return Err("control permutation mismatch".into());
    }
    let mut receipts = Vec::new();
    let mut aggregate = RankOperations::default();
    for (name, modulus, coefficients, _target) in [
        ("2^61-1", &p61, &coefficients61, &target61),
        ("p256-n", &curve.n, &coefficients, &target),
    ] {
        let mut first = coefficients.clone();
        first.push(mod_sub(&BigUint::zero(), &BigUint::one(), modulus));
        first.push(BigUint::zero());
        let mut second = permute_big(coefficients, &permutation);
        second.push(BigUint::zero());
        second.push(mod_sub(&BigUint::zero(), &BigUint::one(), modulus));
        let rows = vec![first, second];
        let (forward, forward_ops) = modular_rank(&rows, modulus, false)?;
        let (reverse, reverse_ops) = modular_rank(&rows, modulus, true)?;
        aggregate.add_assign(&forward_ops);
        aggregate.add_assign(&reverse_ops);
        receipts.push(RankReceipt {
            modulus: name.into(),
            forward_rank: forward,
            reversed_rank: reverse,
            rows: 2,
            columns: CONTROL_WIDTH + 2,
        });
    }
    if receipts
        .iter()
        .any(|receipt| receipt.forward_rank != 2 || receipt.reversed_rank != 2)
    {
        return Err("P-256 augmented rank mismatch".into());
    }
    let scalar_failures = u64::from(big_dot(&coefficients61, &logs61, &p61) != target61)
        + u64::from(
            big_dot(
                &coefficients61,
                &permutation61
                    .iter()
                    .map(|index| logs61[*index].clone())
                    .collect::<Vec<_>>(),
                &p61,
            ) != BigUint::one(),
        );
    let generator = curve.generator();
    let mut operations = GroupOperations::default();
    let points: Vec<Point> = logs
        .iter()
        .map(|log| {
            operations.points_constructed += 1;
            scalar_mul(&generator, log, curve, &mut operations)
        })
        .collect();
    let target_point = scalar_mul(&generator, &target, curve, &mut operations);
    operations.points_constructed += 1;
    let permuted_points: Vec<Point> = permutation
        .iter()
        .map(|index| points[*index].clone())
        .collect();
    let mut group_failures = 0u64;
    group_failures += u64::from(!replay_equation(
        &coefficients,
        &points,
        &target_point,
        curve,
        &mut operations,
    ));
    group_failures += u64::from(!replay_equation(
        &coefficients,
        &permuted_points,
        &generator,
        curve,
        &mut operations,
    ));
    let all_points_on_curve = points
        .iter()
        .all(|point| curve.is_valid_public_point(point))
        && curve.is_valid_public_point(&target_point);
    if scalar_failures != 0 || group_failures != 0 || !all_points_on_curve {
        return Err("P-256 planted replay failure".into());
    }
    Ok(P256Control {
        scalar_labelled_control: true,
        columns: CONTROL_WIDTH,
        permutation,
        target_log_sha256: vector_digest(std::slice::from_ref(&target)),
        log_vector_sha256: vector_digest(&logs),
        coefficient_vector_sha256: vector_digest(&coefficients),
        scalar_replay_failures: scalar_failures,
        group_replay_failures: group_failures,
        all_points_on_curve,
        useful_original_system_rank: 2,
        introduced_auxiliary_log_vector: false,
        rank_receipts: receipts,
        rank_operations: aggregate,
        group_operations: operations,
    })
}

fn binomial(n: usize, k: u32) -> BigUint {
    let k = usize::min(k as usize, n - k as usize);
    let mut value = BigUint::one();
    for index in 0..k {
        value *= BigUint::from(n - index);
        value /= BigUint::from(index + 1);
    }
    value
}

fn signed_domain(columns: usize, arity: u32) -> BigUint {
    binomial(columns, arity) << arity as usize
}

fn big_log2(value: &BigUint) -> f64 {
    value.to_f64().expect("finite f64").log2()
}

fn balanced_list_log2(columns: usize, arity: u32) -> f64 {
    let left = arity / 2;
    let right = arity - left;
    f64::max(
        big_log2(&signed_domain(columns, left)),
        big_log2(&signed_domain(columns, right)),
    )
}

fn domain_projection(columns: usize, curve: &CurveParams) -> Result<DomainProjection, String> {
    if columns < 64 {
        return Err("derived factor base unexpectedly small".into());
    }
    let s17 = signed_domain(columns, 17);
    let n_squared = &curve.n * &curve.n;
    let primary_arity = (2..=64)
        .find(|arity| signed_domain(columns, *arity) >= curve.n)
        .ok_or("no primary-cover arity through 64")?;
    let product_arity = (2..=64)
        .find(|arity| signed_domain(columns, *arity) >= n_squared)
        .ok_or("no product-cover arity through 64")?;
    let s17_log2 = big_log2(&s17);
    let n_log2 = big_log2(&curve.n);
    let primary_mean_log2 = s17_log2 - n_log2;
    let product_mean_log2 = s17_log2 - 2.0 * n_log2;
    let product_list_log2 = balanced_list_log2(columns, product_arity);
    let list_bytes_log2 = product_list_log2 + 82f64.log2();
    let rho_log2 = RHO_S.log2() + 0.5 * (n_log2 - 1.0);
    let ratio_log2 = product_list_log2 - rho_log2;
    Ok(DomainProjection {
        columns,
        relation_arity: 17,
        signed_s17_domain: s17.to_string(),
        signed_s17_domain_log2: s17_log2,
        s17_primary_target_mean: 2f64.powf(primary_mean_log2),
        s17_primary_target_mean_log2: primary_mean_log2,
        s17_product_target_mean: 2f64.powf(product_mean_log2),
        s17_product_target_mean_log2: product_mean_log2,
        minimum_arity_covering_primary_group: primary_arity,
        minimum_arity_covering_product_group: product_arity,
        primary_cover_balanced_list_log2: balanced_list_log2(columns, primary_arity),
        product_cover_balanced_list_log2: product_list_log2,
        product_cover_materialized_bytes_log2: list_bytes_log2,
        product_cover_disk_write_plus_read_log2: list_bytes_log2 + 1.0,
        product_cover_ratio_to_rho: 2f64.powf(ratio_log2),
        product_cover_ratio_to_rho_log2: ratio_log2,
        memoryless_state_bytes: 512,
        structured_degree_of_regularity: None,
        projection_is_measurement: false,
    })
}

fn obligation(name: &str, status: &str, evidence: &str, scope: &str) -> Obligation {
    Obligation {
        name: name.into(),
        status: status.into(),
        evidence: evidence.into(),
        scope: scope.into(),
    }
}

fn assessment(
    factor_base: &FactorBaseReceipt,
    toy: &ToyCensus,
    p256: &P256Control,
    projection: &DomainProjection,
) -> TransferAssessment {
    let mut result = TransferAssessment {
        schema: "p256-dickson-permutation-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 309,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        typed_correspondence_graph: vec![
            format!("{PARENT_FB} --[x and -x intersection]--> {}", factor_base.derived_factor_base),
            "column i --[half-turn]--> column pi(i), with pi^2=id".into(),
            "coefficient row c --[product sum]--> (sum c_i P_i, sum c_i P_pi(i))".into(),
            "product target (Q,G) --[exact hit]--> target-coupling row plus generator anchor on one log vector".into(),
            "pi eigenspace restriction --> dependent coordinates and useful rank one".into(),
        ],
        obligations: vec![
            obligation("actual-base closure", "supported", &format!("{} columns form {} exact half-turn cycles with zero failures.", factor_base.derived_columns, factor_base.two_cycles), "Complete registered-base rebuild."),
            obligation("two-row information rank", "supported only as planted control", &format!("The {}-column labelled control replays rank two with no auxiliary unknown block.", p256.columns), "Information accounting, not construction."),
            obligation("S17 product target", "fails", &format!("Mean hits are 2^{:.6}; the signed domain does not cover GxG.", projection.s17_product_target_mean_log2), "Uniform product-image projection."),
            obligation("product-cover arity and cost", "fails", &format!("First covering arity is {}, with balanced list 2^{:.6} and ratio 2^{:.6} to rho.", projection.minimum_arity_covering_product_group, projection.product_cover_balanced_list_log2, projection.product_cover_ratio_to_rho_log2), "Exact domain and balanced MITM list widths; modeled image uniformity."),
            obligation("derived-base degree", "unknown", "The parent residual degree four is not transferred to the half-turn intersection.", "New derived factor base."),
        ],
        controls: vec![
            format!("Complete toy census covers {} states: 16 eigenprofiles and {} rank-two profiles.", toy.coefficient_states, toy.rank_two_profiles),
            "The registered 131,458-column dump is rebuilt twice before the derived base is emitted in ecbench.factor_base_dump/v1-wide format.".into(),
            "All dependency hashes, schemas, identities, point-set hashes, and imported conclusions are checked.".into(),
        ],
        cost_accounting: vec![
            "Measured: complete base rebuild, closure dictionary, toy census, planted P-256 replay, rank operations, isolation telemetry, and artifact sizes.".into(),
            "Modeled: uniform target-image means and balanced MITM list widths; these are explicitly not attack measurements.".into(),
            "Unset: unplanted product-target event, derived-base degree, relation yield, sparse linear algebra, recovery, and complete attack cost.".into(),
        ],
        exploration_boundary: "Exact for the Dickson exponent half-turn x->-x and its closed intersection; does not cover arbitrary column permutations or non-Dickson correspondences.".into(),
        weakest_open_obligation: "Find a same-base non-scalar correspondence whose two-row target image is searchable in order sqrt(n), while preserving degree at most five.".into(),
        narrowest_supported_finding: "The Dickson half-turn solves the information-rank defect of an independent product lift, but the registered S17 domain is far too small for its product target; the first covering arity requires order-n balanced lists and has unknown degree.".into(),
        semantic_evidence_sha256: String::new(),
        assessment_json_bytes: 0,
    };
    let semantic = json!({
        "graph": result.typed_correspondence_graph,
        "obligations": result.obligations,
        "controls": result.controls,
        "cost": result.cost_accounting,
        "boundary": result.exploration_boundary,
        "open": result.weakest_open_obligation,
        "finding": result.narrowest_supported_finding,
    });
    result.semantic_evidence_sha256 = sha256_hex(&serde_json::to_vec(&semantic).unwrap());
    result
}

fn build_result(cli: &Cli) -> Result<(ResultReceipt, TransferAssessment), String> {
    let deps = dependencies(cli)?;
    let (factor_base, _dump, _permutation) = factor_base_census(&cli.factor_base_out)?;
    let toy = toy_census()?;
    let curve = CurveParams::p256();
    let p256 = p256_control(&curve)?;
    let projection = domain_projection(factor_base.derived_columns, &curve)?;
    let transfer = assessment(&factor_base, &toy, &p256, &projection);
    let gates = Gates {
        dependency_hashes_and_schemas_checked: true,
        actual_factor_base_rebuilt_and_verified: true,
        nonempty_closed_half_turn_subbase: true,
        zero_mapping_and_replay_failures: true,
        complete_toy_census: true,
        p256_two_useful_rows_replayed: true,
        nonlabelled_actual_factor_base_product_target_event: false,
        complete_relation_decomposition_and_recovery: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        complete_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        peak_materialized_storage_below_2_50: false,
        discarded_probabilistic_branches_counted_as_exhaustive: false,
        promoted: false,
    };
    let classification =
        "same-base-permutation/two-useful-rows/product-domain-insufficient/order-n-mitm/degree-unset";
    let obstruction = format!(
        "The half-turn intersection removes the new-unknown defect and a planted (Q,G) hit has useful rank two, but the S17 product-target mean is only 2^{:.6}. The first arity whose exact signed domain covers GxG is {}, whose balanced MITM list is 2^{:.6} entries, 2^{:.6} times rho before solver or recovery; the derived base's structured degree is also unset.",
        projection.s17_product_target_mean_log2,
        projection.minimum_arity_covering_product_group,
        projection.product_cover_balanced_list_log2,
        projection.product_cover_ratio_to_rho_log2,
    );
    let decision = "Reject the Dickson half-turn permutation as a rho-parity route. Preserve its exact derived factor base and two-row information certificate, but do not transfer degree evidence or attempt an unplanted relation. Require a same-base correspondence with a square-root-size searchable image rather than a full product target.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": &deps,
        "factor_base": &factor_base,
        "toy": &toy,
        "p256": &p256,
        "projection": &projection,
        "assessment": transfer.semantic_evidence_sha256,
        "gates": &gates,
        "classification": classification,
        "obstruction": obstruction,
        "decision": decision,
    });
    let semantic_hash = sha256_hex(&serde_json::to_vec(&semantic).unwrap());
    Ok((
        ResultReceipt {
            schema: "p256-dickson-permutation-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 309,
            execution_status: "complete".into(),
            parent_factor_base: PARENT_FB.into(),
            dependencies: deps,
            typed_map: "Phi(c)=(sum c_i P_i,sum c_i P_pi(i)); pi is the x->-x half-turn involution on the derived same-base intersection".into(),
            factor_base,
            exhaustive_toy_census: toy,
            p256_control: p256,
            domain_projection: projection,
            relations_reported_on_actual_factor_base: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            transfer_assessment_semantic_sha256: transfer.semantic_evidence_sha256.clone(),
            gates,
            classification: classification.into(),
            dominant_obstruction: obstruction,
            decision: decision.into(),
            semantic_evidence_sha256: semantic_hash,
            result_json_bytes: 0,
        },
        transfer,
    ))
}

fn write_result(path: &Path, result: &mut ResultReceipt) -> Result<(), String> {
    loop {
        let text = serde_json::to_string_pretty(result).map_err(|error| error.to_string())? + "\n";
        if text.len() as u64 == result.result_json_bytes {
            return fs::write(path, text).map_err(|error| error.to_string());
        }
        result.result_json_bytes = text.len() as u64;
    }
}

fn write_assessment(path: &Path, result: &mut TransferAssessment) -> Result<(), String> {
    loop {
        let text = serde_json::to_string_pretty(result).map_err(|error| error.to_string())? + "\n";
        if text.len() as u64 == result.assessment_json_bytes {
            return fs::write(path, text).map_err(|error| error.to_string());
        }
        result.assessment_json_bytes = text.len() as u64;
    }
}

fn run(cli: Cli) -> Result<(), String> {
    let (mut result, mut assessment) = build_result(&cli)?;
    write_assessment(&cli.assessment, &mut assessment)?;
    write_result(&cli.out, &mut result)?;
    eprintln!(
        "round 309: derived_fb={}, columns={}, s17_product_log2_mean={:.6}, cover_arity={}, list_log2={:.6}, promoted={}",
        result.factor_base.derived_factor_base,
        result.factor_base.derived_columns,
        result.domain_projection.s17_product_target_mean_log2,
        result.domain_projection.minimum_arity_covering_product_group,
        result.domain_projection.product_cover_balanced_list_log2,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_dickson_permutation: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_toy_census_matches_eigenspace_reference() {
        let census = toy_census().expect("toy census");
        assert_eq!(census.plus_eigenprofiles, 8);
        assert_eq!(census.minus_eigenprofiles, 8);
        assert_eq!(census.rank_two_profiles, 384);
        assert_eq!(census.rank_two_product_target_hits, 18_816);
        assert_eq!(census.eigenprofile_product_target_hits, 0);
    }

    #[test]
    fn planted_p256_rows_have_rank_two_without_new_unknowns() {
        let control = p256_control(&CurveParams::p256()).expect("P-256 control");
        assert_eq!(control.useful_original_system_rank, 2);
        assert!(!control.introduced_auxiliary_log_vector);
        assert_eq!(control.group_replay_failures, 0);
        assert!(control
            .rank_receipts
            .iter()
            .all(|receipt| receipt.forward_rank == 2));
    }

    #[test]
    fn projected_product_cover_fails_time_and_storage() {
        let projection = domain_projection(65_000, &CurveParams::p256()).expect("projection");
        assert!(projection.s17_product_target_mean_log2 < -200.0);
        assert!(projection.product_cover_balanced_list_log2 > 250.0);
        assert!(projection.product_cover_materialized_bytes_log2 > 250.0);
    }
}
