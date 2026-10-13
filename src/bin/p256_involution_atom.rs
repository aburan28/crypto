//! Round 311: exact eigenspace-atom audit for the P-256 Dickson half-turn.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::{sha256_hex, short_id};
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    WideFactorBaseDump, WideFactorBaseFacts, WideFbPoint, CURVE_SLUG, FACTOR_BASE_SCHEMA,
    POINT_KEY_ENCODING,
};
use crypto_lib::ecc::{CurveParams, FieldElement, Point};
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const ROUND294_SHA256: &str = "9716fc1ba7d24a023627e8701094ca870b2121a84f7a0148784f70565640815a";
const ROUND309_SHA256: &str = "1fa17da4d11345f2dce55d6704355d193e081eae6912d3b45a9605e2efa32bc1";
const ROUND309_DUMP_SHA256: &str =
    "677f5ddbe0f2ea1e0c2bd29272d76461acc0d0197d11843d3ebf77d91eaf15ac";
const ROUND310_SHA256: &str = "f8c39500b97394e633dcbcd36a92bf2cb896176aff7e4ffcb357e201b4dbfdc8";
const REGISTERED_FB: &str = "FB1h2f8621cda105";
const HALF_TURN_FB: &str = "FB1h6255ce9746fe";
const HALF_TURN_COLUMNS: usize = 65_724;
const EXPECTED_CYCLES: usize = 32_862;
const REQUESTED_ROWS: u64 = 138_031;
const ATOM_ARITY: u32 = 20;
const HALF_ARITY: u32 = 10;
const ENTRY_BYTES: u64 = 82;
const RHO_S: f64 = 1.3;
const TOY_MODULUS: u8 = 7;
const TOY_WIDTH: usize = 4;

#[derive(Parser)]
#[command(about = "Audit eigenspace atoms of the P-256 Dickson half-turn")]
struct Cli {
    #[arg(long)]
    round294: PathBuf,
    #[arg(long)]
    round309: PathBuf,
    #[arg(long)]
    factor_base: PathBuf,
    #[arg(long)]
    round310: PathBuf,
    #[arg(long)]
    atom_base_out: PathBuf,
    #[arg(long)]
    out: PathBuf,
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
    atom_additions: u64,
    identity_replay_additions: u64,
    target_split_additions: u64,
    target_split_scalar_multiplications: u64,
    small_scalar_additions: u64,
}

#[derive(Clone, Serialize)]
struct ScalarScreen {
    scalar: u32,
    plus_to_plus_matches: usize,
    plus_to_minus_matches: usize,
    minus_to_plus_matches: usize,
    minus_to_minus_matches: usize,
    complete_cross_transport: bool,
}

#[derive(Clone, Serialize)]
struct AtomCensus {
    parent_columns: usize,
    two_cycles: usize,
    plus_atoms: usize,
    minus_atoms: usize,
    plus_unique_up_to_negation: usize,
    minus_unique_up_to_negation: usize,
    union_unique_up_to_negation: usize,
    plus_internal_duplicate_columns: usize,
    minus_internal_duplicate_columns: usize,
    cross_channel_overlaps: usize,
    plus_parent_overlaps: usize,
    minus_parent_overlaps: usize,
    infinity_atoms: u64,
    invalid_parent_points: u64,
    invalid_atom_points: u64,
    cycle_rebuild_failures: u64,
    group_identity_replay_failures: u64,
    coordinate_identity_failures: u64,
    plus_points_sha256: String,
    minus_points_sha256: String,
    union_points_sha256: String,
    local_transform_determinant: String,
    aggregate_transform_rank: usize,
    independent_log_dimension_reduction: usize,
    complete_small_scalar_cross_transport: bool,
    scalar_screens: Vec<ScalarScreen>,
    group_operations: GroupOperations,
}

#[derive(Clone, Serialize)]
struct AtomFactorBaseReceipt {
    factor_base: String,
    factor_base_sha256: String,
    family: String,
    columns: usize,
    signed_points: usize,
    points_sha256: String,
    dump_schema: String,
    dump_path: String,
    dump_bytes: u64,
    dump_sha256: String,
}

#[derive(Clone, Default, Serialize)]
struct ToyOperations {
    profiles: u64,
    coefficient_states: u64,
    direct_dot_terms: u64,
    atom_dot_terms: u64,
    rank_matrices: u64,
}

#[derive(Clone, Serialize)]
struct ToyCensus {
    modulus: u8,
    parent_width: usize,
    cycles: usize,
    projective_profiles: u64,
    coefficient_states_per_profile: u64,
    total_states: u64,
    factorization_discrepancies: u64,
    inverse_transform_discrepancies: u64,
    forward_reverse_rank_discrepancies: u64,
    augmented_rank_discrepancies: u64,
    false_positives: u64,
    false_negatives: u64,
    operations: ToyOperations,
}

#[derive(Clone, Serialize)]
struct P256TargetControl {
    scalar_labelled_control: bool,
    target_log_sha256: String,
    inverse_two: String,
    split_replay_failures: u64,
    all_points_on_curve: bool,
    group_operations: GroupOperations,
}

#[derive(Clone, Serialize)]
struct TargetDiagnostic {
    channel_columns: usize,
    atom_arity_per_channel: u32,
    original_point_terms_upper_bound: u32,
    deterministic_partition: String,
    left_columns: usize,
    right_columns: usize,
    half_arity: u32,
    left_list_entries: String,
    right_list_entries: String,
    left_list_log2: f64,
    right_list_log2: f64,
    one_target_domain: String,
    one_target_domain_log2: f64,
    one_target_mean_log2: f64,
    two_target_mean_log2: f64,
    two_target_total_entries: String,
    two_target_total_entries_log2: f64,
    optimistic_stage_s: f64,
    optimistic_stage_ratio_to_rho: f64,
    optimistic_stage_ratio_to_rho_log2: f64,
    materialized_bytes: String,
    materialized_bytes_log2: f64,
    disk_write_plus_read_bytes: String,
    disk_write_plus_read_log2: f64,
    assumes_atom_logs_available: bool,
    is_s17_relation: bool,
    is_end_to_end_measurement: bool,
}

#[derive(Clone, Serialize)]
struct CollectionProjection {
    requested_usable_rows: u64,
    larger_channel_rows: u64,
    smaller_channel_rows: u64,
    atom_log_variables: usize,
    nonzeros_per_atom_row: u32,
    exact_kth_collision_formula: String,
    optimistic_free_joint_sample_s: f64,
    optimistic_free_joint_sample_ratio_to_rho: f64,
    charged_two_group_sum_s: f64,
    charged_two_group_sum_ratio_to_rho: f64,
    charged_collection_operations_log2: f64,
    charged_collection_per_usable_row_log2: f64,
    sparse_linear_algebra_operations: String,
    sparse_linear_algebra_operations_log2: f64,
    target_decomposition_operations_log2: f64,
    optimistic_complete_operations_log2: f64,
    optimistic_complete_s: f64,
    optimistic_complete_ratio_to_rho: f64,
    memoryless_state_bytes: u64,
    materialized_target_storage_log2: f64,
    projection_is_measurement: bool,
    projection_is_exhaustive_search: bool,
    omitted_favourable_costs: Vec<String>,
}

#[derive(Clone, Serialize)]
struct BoundaryRow {
    variant: String,
    s_group_addition_equivalents: f64,
    ratio_to_rho: f64,
    correctness: String,
    evidence_kind: String,
    class: String,
}

#[derive(Clone, Serialize)]
struct ImportedBoundary {
    round294_registered_base_exact_kth_ratio_to_rho: f64,
    round309_s17_product_target_mean_log2: f64,
    round309_product_cover_ratio_to_rho_log2: f64,
    round309_structured_degree_of_regularity: Option<u64>,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_schemas_and_curve_checked: bool,
    complete_actual_half_turn_census: bool,
    zero_actual_replay_and_coordinate_failures: bool,
    complete_toy_census: bool,
    zero_toy_false_positives_and_false_negatives: bool,
    actual_nonlabelled_dimension_reduction: bool,
    complete_known_scalar_cross_transport: bool,
    target_stage_at_or_below_rho: bool,
    structured_residual_degree_at_most_5: bool,
    complete_relation_decomposition_and_recovery: bool,
    complete_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    peak_materialized_storage_below_2_50: bool,
    parity_at_or_below_rho: bool,
    discarded_probabilistic_branches_counted_as_exhaustive: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    screening_round: u64,
    execution_status: String,
    registered_factor_base: String,
    half_turn_factor_base: String,
    dependencies: Vec<Dependency>,
    construction: String,
    correction_scope: String,
    atom_census: AtomCensus,
    atom_factor_base: AtomFactorBaseReceipt,
    exhaustive_toy_census: ToyCensus,
    p256_target_control: P256TargetControl,
    target_diagnostic: TargetDiagnostic,
    relation_collection_projection: CollectionProjection,
    imported_boundary: ImportedBoundary,
    one_unit_boundary_table: Vec<BoundaryRow>,
    structured_residual_degree_of_regularity: Option<u64>,
    relations_reported_on_actual_factor_base: u64,
    full_depth_unplanted_p256_relation_attempted: bool,
    gates: Gates,
    classification: String,
    dominant_obstruction: String,
    decision: String,
    semantic_evidence_sha256: String,
    result_json_bytes: u64,
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("not a non-negative integer: {value}"))
}

fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn checked_json(path: &Path, expected: &str, schema: &str) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != expected {
        return Err(format!(
            "{} hash {digest}, expected {expected}",
            path.display()
        ));
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

fn checked_dump(path: &Path) -> Result<(Dependency, WideFactorBaseDump), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != ROUND309_DUMP_SHA256 {
        return Err(format!(
            "{} hash {digest}, expected {ROUND309_DUMP_SHA256}",
            path.display()
        ));
    }
    let dump: WideFactorBaseDump =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if dump.schema != "ecbench.factor_base_dump/v1-wide"
        || dump.curve.slug != CURVE_SLUG
        || dump.factor_base.fb_id != HALF_TURN_FB
        || dump.factor_base.columns != HALF_TURN_COLUMNS as u64
        || dump.points.len() != 2 * HALF_TURN_COLUMNS
    {
        return Err("Round 309 dump identity mismatch".into());
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            sha256: digest,
            schema: dump.schema.clone(),
        },
        dump,
    ))
}

fn dependencies(
    cli: &Cli,
) -> Result<(Vec<Dependency>, ImportedBoundary, WideFactorBaseDump), String> {
    let (dep294, round294) = checked_json(
        &cli.round294,
        ROUND294_SHA256,
        "p256-variable-arity-closure/v2",
    )?;
    let (dep309, round309) = checked_json(
        &cli.round309,
        ROUND309_SHA256,
        "p256-dickson-permutation-screen/v1",
    )?;
    let (dep_dump, dump) = checked_dump(&cli.factor_base)?;
    let (dep310, round310) = checked_json(
        &cli.round310,
        ROUND310_SHA256,
        "p256-permutation-dichotomy-screen/v1",
    )?;
    if round294
        .get("comparison_factor_base")
        .and_then(Value::as_str)
        != Some(REGISTERED_FB)
        || round309
            .pointer("/factor_base/derived_factor_base")
            .and_then(Value::as_str)
            != Some(HALF_TURN_FB)
        || round309
            .pointer("/factor_base/derived_columns")
            .and_then(Value::as_u64)
            != Some(HALF_TURN_COLUMNS as u64)
        || round309
            .pointer("/gates/p256_two_useful_rows_replayed")
            .and_then(Value::as_bool)
            != Some(true)
        || round310
            .pointer("/gates/affine_transport_iff_rank_one_verified")
            .and_then(Value::as_bool)
            != Some(true)
        || round310
            .pointer("/p256_control/half_turn_difference_image_rank")
            .and_then(Value::as_u64)
            != Some(2)
    {
        return Err("imported conclusion mismatch".into());
    }
    let registered_ratio = round294
        .get("dickson_inventory")
        .and_then(Value::as_array)
        .and_then(|rows| {
            rows.iter()
                .find(|row| row.get("fb_id").and_then(Value::as_str) == Some(REGISTERED_FB))
        })
        .and_then(|row| row.get("exact_kth_collision_ratio_to_rho"))
        .and_then(Value::as_f64)
        .ok_or("Round 294 registered boundary missing")?;
    let product_mean = round309
        .pointer("/domain_projection/s17_product_target_mean_log2")
        .and_then(Value::as_f64)
        .ok_or("Round 309 S17 product mean missing")?;
    let product_ratio = round309
        .pointer("/domain_projection/product_cover_ratio_to_rho_log2")
        .and_then(Value::as_f64)
        .ok_or("Round 309 product ratio missing")?;
    if (registered_ratio - 394.425_280_345_252_2).abs() > 1e-12
        || (product_mean + 271.270_333_552_737_44).abs() > 1e-12
        || (product_ratio - 151.122_586_704_556_74).abs() > 1e-12
    {
        return Err("imported numeric boundary mismatch".into());
    }
    Ok((
        vec![dep294, dep309, dep_dump, dep310],
        ImportedBoundary {
            round294_registered_base_exact_kth_ratio_to_rho: registered_ratio,
            round309_s17_product_target_mean_log2: product_mean,
            round309_product_cover_ratio_to_rho_log2: product_ratio,
            round309_structured_degree_of_regularity: round309
                .pointer("/domain_projection/structured_degree_of_regularity")
                .and_then(Value::as_u64),
        },
        dump,
    ))
}

fn point_from_row(row: &WideFbPoint, curve: &CurveParams) -> Result<Point, String> {
    Ok(Point::Affine {
        x: FieldElement::new(parse_big(&row.x)?, curve.p.clone()),
        y: FieldElement::new(parse_big(&row.y)?, curve.p.clone()),
    })
}

fn point_xy(point: &Point) -> Result<(&BigUint, &BigUint), String> {
    match point {
        Point::Affine { x, y } => Ok((&x.value, &y.value)),
        Point::Infinity => Err("point is infinity".into()),
    }
}

fn canonical_y(y: &BigUint, p: &BigUint) -> BigUint {
    let neg = p - y;
    if y <= &neg {
        y.clone()
    } else {
        neg
    }
}

fn point_key(x: &BigUint, high_sign: bool) -> [u8; 33] {
    let key = ((x + BigUint::one()) << 1usize) + BigUint::from(high_sign);
    let raw = key.to_bytes_be();
    let mut out = [0u8; 33];
    out[33 - raw.len()..].copy_from_slice(&raw);
    out
}

fn x_digest(values: &BTreeSet<BigUint>) -> String {
    let mut bytes = Vec::with_capacity(values.len() * 32);
    for value in values {
        let raw = value.to_bytes_be();
        bytes.extend(std::iter::repeat_n(0, 32 - raw.len()));
        bytes.extend(raw);
    }
    sha256_hex(&bytes)
}

fn selected_multiples(
    point: &Point,
    a: &FieldElement,
    operations: &mut GroupOperations,
) -> [Point; 5] {
    let two = point.double_vartime(a);
    let three = two.add_vartime(point, a);
    let four = two.double_vartime(a);
    let eight = four.double_vartime(a);
    let sixteen = eight.double_vartime(a);
    operations.small_scalar_additions += 5;
    [two, three, four, eight, sixteen]
}

fn write_atom_dump(
    path: &Path,
    parent: &WideFactorBaseDump,
    union: &BTreeMap<BigUint, BigUint>,
    curve: &CurveParams,
) -> Result<AtomFactorBaseReceipt, String> {
    let mut point_blob = Vec::with_capacity(union.len() * 66);
    let mut points = Vec::with_capacity(union.len() * 2);
    for (column, (x, low_y)) in union.iter().enumerate() {
        let high_y = &curve.p - low_y;
        point_blob.extend_from_slice(&point_key(x, false));
        point_blob.extend_from_slice(&point_key(x, true));
        points.push(WideFbPoint {
            x: lower_hex(x),
            y: lower_hex(low_y),
            col: column as u64,
            coef: "1".into(),
        });
        points.push(WideFbPoint {
            x: lower_hex(x),
            y: lower_hex(&high_y),
            col: column as u64,
            coef: (&curve.n - BigUint::one()).to_string(),
        });
    }
    let points_sha256 = sha256_hex(&point_blob);
    let mut params = BTreeMap::new();
    params.insert("parent_fb_id".into(), HALF_TURN_FB.into());
    params.insert("point_key_encoding".into(), POINT_KEY_ENCODING.into());
    params.insert("transform".into(), "A=P+P_pi;B=P-P_pi".into());
    params.insert(
        "deduplication".into(),
        "global negation folded by affine x".into(),
    );
    let identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": CURVE_SLUG,
        "family": "dickson-half-turn-eigenspace-atoms",
        "params": params,
        "columns": union.len(),
        "signed_points": 2 * union.len(),
        "points_sha256": points_sha256,
    });
    let (fb_id, fb_sha256) = short_id("FB1", &identity)?;
    let dump = WideFactorBaseDump {
        schema: parent.schema.clone(),
        curve: parent.curve.clone(),
        factor_base: WideFactorBaseFacts {
            fb_id: fb_id.clone(),
            fb_sha256: fb_sha256.clone(),
            family: "dickson-half-turn-eigenspace-atoms".into(),
            params,
            description: format!(
                "negation-folded union of P+P_pi and P-P_pi atoms from {HALF_TURN_FB} on {CURVE_SLUG}"
            ),
            signed_points: (2 * union.len()) as u64,
            abscissae: union.len() as u64,
            columns: union.len() as u64,
            dimension: None,
            points_sha256: points_sha256.clone(),
        },
        build_adds: (2 * EXPECTED_CYCLES) as u64,
        build_doubles: 0,
        points,
    };
    let text = serde_json::to_string_pretty(&dump).map_err(|error| error.to_string())? + "\n";
    fs::write(path, &text).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok(AtomFactorBaseReceipt {
        factor_base: fb_id,
        factor_base_sha256: fb_sha256,
        family: dump.factor_base.family,
        columns: union.len(),
        signed_points: 2 * union.len(),
        points_sha256,
        dump_schema: dump.schema,
        dump_path: path.display().to_string(),
        dump_bytes: text.len() as u64,
        dump_sha256: sha256_hex(text.as_bytes()),
    })
}

fn atom_census(
    dump: &WideFactorBaseDump,
    path: &Path,
    curve: &CurveParams,
) -> Result<(AtomCensus, AtomFactorBaseReceipt), String> {
    let positives: Vec<&WideFbPoint> = dump.points.chunks_exact(2).map(|pair| &pair[0]).collect();
    let mut by_x = BTreeMap::<BigUint, usize>::new();
    let mut parent_x = BTreeSet::new();
    for (index, row) in positives.iter().enumerate() {
        let x = parse_big(&row.x)?;
        if by_x.insert(x.clone(), index).is_some() {
            return Err("parent factor base has duplicate x".into());
        }
        parent_x.insert(x);
    }
    let a_fe = curve.a_fe();
    let mut plus_points = Vec::with_capacity(EXPECTED_CYCLES);
    let mut minus_points = Vec::with_capacity(EXPECTED_CYCLES);
    let mut plus_x = BTreeSet::new();
    let mut minus_x = BTreeSet::new();
    let mut union = BTreeMap::<BigUint, BigUint>::new();
    let mut cycles = 0usize;
    let mut infinity = 0u64;
    let mut invalid_parent = 0u64;
    let mut invalid_atom = 0u64;
    let mut cycle_failures = 0u64;
    let mut group_failures = 0u64;
    let mut coordinate_failures = 0u64;
    let mut operations = GroupOperations::default();
    for (index, row) in positives.iter().enumerate() {
        let x = parse_big(&row.x)?;
        let neg_x = if x.is_zero() {
            BigUint::zero()
        } else {
            &curve.p - &x
        };
        let Some(mate) = by_x.get(&neg_x).copied() else {
            cycle_failures += 1;
            continue;
        };
        if index >= mate {
            continue;
        }
        cycles += 1;
        let p = point_from_row(row, curve)?;
        let r = point_from_row(positives[mate], curve)?;
        invalid_parent += u64::from(!curve.is_valid_public_point(&p));
        invalid_parent += u64::from(!curve.is_valid_public_point(&r));
        let plus = p.add_vartime(&r, &a_fe);
        let minus = p.add_vartime(&r.neg_vartime(), &a_fe);
        operations.atom_additions += 2;
        if matches!(plus, Point::Infinity) || matches!(minus, Point::Infinity) {
            infinity += 1;
            continue;
        }
        invalid_atom += u64::from(!curve.is_valid_public_point(&plus));
        invalid_atom += u64::from(!curve.is_valid_public_point(&minus));
        let p2 = p.double_vartime(&a_fe);
        let r2 = r.double_vartime(&a_fe);
        let plus_plus_minus = plus.add_vartime(&minus, &a_fe);
        let plus_minus_minus = plus.add_vartime(&minus.neg_vartime(), &a_fe);
        operations.identity_replay_additions += 4;
        group_failures += u64::from(plus_plus_minus != p2);
        group_failures += u64::from(plus_minus_minus != r2);

        let (xp, yp) = point_xy(&plus)?;
        let (xm, ym) = point_xy(&minus)?;
        let x2 = (&x * &x) % &curve.p;
        let lhs_sum = (((xp + xm) % &curve.p) * &x2) % &curve.p;
        let lhs_product = (((BigUint::from(4u8) * &x2) % &curve.p * xp) % &curve.p * xm) % &curve.p;
        let x2_plus_a = (&x2 + &curve.a) % &curve.p;
        let rhs_product = (&x2_plus_a * &x2_plus_a) % &curve.p;
        coordinate_failures += u64::from(lhs_sum != curve.b);
        coordinate_failures += u64::from(lhs_product != rhs_product);

        plus_x.insert(xp.clone());
        minus_x.insert(xm.clone());
        union
            .entry(xp.clone())
            .or_insert_with(|| canonical_y(yp, &curve.p));
        union
            .entry(xm.clone())
            .or_insert_with(|| canonical_y(ym, &curve.p));
        plus_points.push(plus);
        minus_points.push(minus);
    }
    if cycles != EXPECTED_CYCLES {
        return Err(format!(
            "rebuilt {cycles} cycles, expected {EXPECTED_CYCLES}"
        ));
    }

    let scalars = [2u32, 3, 4, 8, 16];
    let mut counts = vec![[0usize; 4]; scalars.len()];
    for point in &plus_points {
        for (slot, multiple) in selected_multiples(point, &a_fe, &mut operations)
            .into_iter()
            .enumerate()
        {
            let (x, _) = point_xy(&multiple)?;
            counts[slot][0] += usize::from(plus_x.contains(x));
            counts[slot][1] += usize::from(minus_x.contains(x));
        }
    }
    for point in &minus_points {
        for (slot, multiple) in selected_multiples(point, &a_fe, &mut operations)
            .into_iter()
            .enumerate()
        {
            let (x, _) = point_xy(&multiple)?;
            counts[slot][2] += usize::from(plus_x.contains(x));
            counts[slot][3] += usize::from(minus_x.contains(x));
        }
    }
    let scalar_screens: Vec<ScalarScreen> = scalars
        .iter()
        .enumerate()
        .map(|(slot, scalar)| ScalarScreen {
            scalar: *scalar,
            plus_to_plus_matches: counts[slot][0],
            plus_to_minus_matches: counts[slot][1],
            minus_to_plus_matches: counts[slot][2],
            minus_to_minus_matches: counts[slot][3],
            complete_cross_transport: counts[slot][1] == plus_points.len()
                || counts[slot][2] == minus_points.len(),
        })
        .collect();
    let complete_transport = scalar_screens
        .iter()
        .any(|screen| screen.complete_cross_transport);
    let cross = plus_x.intersection(&minus_x).count();
    let plus_parent = plus_x.intersection(&parent_x).count();
    let minus_parent = minus_x.intersection(&parent_x).count();
    let union_x: BTreeSet<BigUint> = union.keys().cloned().collect();
    let atom_dump = write_atom_dump(path, dump, &union, curve)?;
    let reduction = 2 * cycles - union.len();
    let census = AtomCensus {
        parent_columns: positives.len(),
        two_cycles: cycles,
        plus_atoms: plus_points.len(),
        minus_atoms: minus_points.len(),
        plus_unique_up_to_negation: plus_x.len(),
        minus_unique_up_to_negation: minus_x.len(),
        union_unique_up_to_negation: union.len(),
        plus_internal_duplicate_columns: plus_points.len() - plus_x.len(),
        minus_internal_duplicate_columns: minus_points.len() - minus_x.len(),
        cross_channel_overlaps: cross,
        plus_parent_overlaps: plus_parent,
        minus_parent_overlaps: minus_parent,
        infinity_atoms: infinity,
        invalid_parent_points: invalid_parent,
        invalid_atom_points: invalid_atom,
        cycle_rebuild_failures: cycle_failures,
        group_identity_replay_failures: group_failures,
        coordinate_identity_failures: coordinate_failures,
        plus_points_sha256: x_digest(&plus_x),
        minus_points_sha256: x_digest(&minus_x),
        union_points_sha256: x_digest(&union_x),
        local_transform_determinant: (&curve.n - BigUint::from(2u8)).to_string(),
        aggregate_transform_rank: 2 * cycles,
        independent_log_dimension_reduction: reduction,
        complete_small_scalar_cross_transport: complete_transport,
        scalar_screens,
        group_operations: operations,
    };
    if census.infinity_atoms
        + census.invalid_parent_points
        + census.invalid_atom_points
        + census.cycle_rebuild_failures
        + census.group_identity_replay_failures
        + census.coordinate_identity_failures
        != 0
        || atom_dump.columns != census.union_unique_up_to_negation
    {
        return Err("actual atom census correctness failure".into());
    }
    Ok((census, atom_dump))
}

fn toy_inverse(value: u8) -> u8 {
    (1..TOY_MODULUS)
        .find(|candidate| value * candidate % TOY_MODULUS == 1)
        .expect("nonzero inverse")
}

fn toy_sub(left: u8, right: u8) -> u8 {
    (left + TOY_MODULUS - right) % TOY_MODULUS
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
        let inverse = toy_inverse(first);
        for entry in &mut row {
            *entry = *entry * inverse % TOY_MODULUS;
        }
        set.insert(row);
    }
    set.into_iter().collect()
}

fn toy_dot(left: &[u8], right: &[u8], terms: &mut u64) -> u8 {
    left.iter().zip(right).fold(0u8, |sum, (a, b)| {
        *terms += 1;
        (sum + a * b) % TOY_MODULUS
    })
}

fn toy_rank(rows: &[Vec<u8>], reverse: bool) -> usize {
    let mut matrix = rows.to_vec();
    if reverse {
        matrix.reverse();
    }
    let mut rank = 0usize;
    for column in 0..matrix[0].len() {
        let Some(pivot) = (rank..matrix.len()).find(|row| matrix[*row][column] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inverse = toy_inverse(matrix[rank][column]);
        for entry in &mut matrix[rank][column..] {
            *entry = *entry * inverse % TOY_MODULUS;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..matrix.len() {
            let factor = matrix[row][column];
            for col in column..matrix[row].len() {
                matrix[row][col] = toy_sub(matrix[row][col], factor * pivot_row[col] % TOY_MODULUS);
            }
        }
        rank += 1;
    }
    rank
}

fn toy_census() -> Result<ToyCensus, String> {
    let profiles = toy_profiles();
    let states: Vec<Vec<u8>> = (0..2401u16).map(decode_base7).collect();
    let inv2 = toy_inverse(2);
    let mut factorization = 0u64;
    let mut inverse_discrepancies = 0u64;
    let mut order_discrepancies = 0u64;
    let mut rank_discrepancies = 0u64;
    let mut false_positives = 0u64;
    let mut false_negatives = 0u64;
    let mut operations = ToyOperations {
        profiles: profiles.len() as u64,
        ..ToyOperations::default()
    };
    for logs in &profiles {
        let permuted = vec![logs[1], logs[0], logs[3], logs[2]];
        let atom_plus = vec![
            (logs[0] + logs[1]) % TOY_MODULUS,
            (logs[2] + logs[3]) % TOY_MODULUS,
        ];
        let atom_minus = vec![toy_sub(logs[0], logs[1]), toy_sub(logs[2], logs[3])];
        for state in &states {
            operations.coefficient_states += 1;
            let (a0, a1, b0, b1) = (state[0], state[1], state[2], state[3]);
            let plus = [a0, a0, a1, a1];
            let minus = [b0, toy_sub(0, b0), b1, toy_sub(0, b1)];
            let coefficients: Vec<u8> = plus
                .iter()
                .zip(&minus)
                .map(|(left, right)| (left + right) % TOY_MODULUS)
                .collect();
            let primary = toy_dot(&coefficients, logs, &mut operations.direct_dot_terms);
            let auxiliary = toy_dot(&coefficients, &permuted, &mut operations.direct_dot_terms);
            let plus_sum = toy_dot(&[a0, a1], &atom_plus, &mut operations.atom_dot_terms);
            let minus_sum = toy_dot(&[b0, b1], &atom_minus, &mut operations.atom_dot_terms);
            factorization += u64::from(
                primary != (plus_sum + minus_sum) % TOY_MODULUS
                    || auxiliary != toy_sub(plus_sum, minus_sum),
            );
            let recovered_plus = (primary + auxiliary) % TOY_MODULUS * inv2 % TOY_MODULUS;
            let recovered_minus = toy_sub(primary, auxiliary) * inv2 % TOY_MODULUS;
            inverse_discrepancies +=
                u64::from(recovered_plus != plus_sum || recovered_minus != minus_sum);
            let direct_target_hit = primary == 2 && auxiliary == 1;
            let factorized_target_hit = plus_sum == 3 * inv2 % TOY_MODULUS && minus_sum == inv2;
            false_positives += u64::from(factorized_target_hit && !direct_target_hit);
            false_negatives += u64::from(direct_target_hit && !factorized_target_hit);
            let mut row_plus: Vec<u8> = plus.iter().map(|value| 2 * value % TOY_MODULUS).collect();
            row_plus.extend([TOY_MODULUS - 1, TOY_MODULUS - 1]);
            let mut row_minus: Vec<u8> =
                minus.iter().map(|value| 2 * value % TOY_MODULUS).collect();
            row_minus.extend([TOY_MODULUS - 1, 1]);
            let rows = vec![row_plus, row_minus];
            let rank = toy_rank(&rows, false);
            let reverse = toy_rank(&rows, true);
            operations.rank_matrices += 2;
            order_discrepancies += u64::from(rank != reverse);
            rank_discrepancies += u64::from(rank != 2);
        }
    }
    let census = ToyCensus {
        modulus: TOY_MODULUS,
        parent_width: TOY_WIDTH,
        cycles: TOY_WIDTH / 2,
        projective_profiles: profiles.len() as u64,
        coefficient_states_per_profile: states.len() as u64,
        total_states: operations.coefficient_states,
        factorization_discrepancies: factorization,
        inverse_transform_discrepancies: inverse_discrepancies,
        forward_reverse_rank_discrepancies: order_discrepancies,
        augmented_rank_discrepancies: rank_discrepancies,
        false_positives,
        false_negatives,
        operations,
    };
    if census.projective_profiles != 400
        || census.coefficient_states_per_profile != 2_401
        || census.total_states != 960_400
        || census.factorization_discrepancies
            + census.inverse_transform_discrepancies
            + census.forward_reverse_rank_discrepancies
            + census.augmented_rank_discrepancies
            + census.false_positives
            + census.false_negatives
            != 0
    {
        return Err("complete toy atom factorization mismatch".into());
    }
    Ok(census)
}

fn hash_scalar(label: &str, modulus: &BigUint) -> BigUint {
    let digest = sha256_hex(format!("round311/involution-atom/{label}").as_bytes());
    let value = BigUint::parse_bytes(digest.as_bytes(), 16).unwrap();
    value % (modulus - BigUint::one()) + BigUint::one()
}

fn p256_target_control(curve: &CurveParams) -> P256TargetControl {
    let d = hash_scalar("target", &curve.n);
    let inv2 = (&curve.n + BigUint::one()) >> 1usize;
    let x_log = ((&d + BigUint::one()) * &inv2) % &curve.n;
    let y_log = ((&d + &curve.n - BigUint::one()) * &inv2) % &curve.n;
    let g = curve.generator();
    let a = curve.a_fe();
    let q = g.scalar_mul(&d, &a);
    let x = g.scalar_mul(&x_log, &a);
    let y = g.scalar_mul(&y_log, &a);
    let mut failures = 0u64;
    failures += u64::from(x.add(&y, &a) != q);
    failures += u64::from(x.add(&y.neg(), &a) != g);
    P256TargetControl {
        scalar_labelled_control: true,
        target_log_sha256: sha256_hex(d.to_string().as_bytes()),
        inverse_two: inv2.to_string(),
        split_replay_failures: failures,
        all_points_on_curve: curve.is_valid_public_point(&q)
            && curve.is_valid_public_point(&x)
            && curve.is_valid_public_point(&y),
        group_operations: GroupOperations {
            target_split_additions: 2,
            target_split_scalar_multiplications: 3,
            ..GroupOperations::default()
        },
    }
}

fn choose(n: usize, k: u32) -> BigUint {
    let mut value = BigUint::one();
    for index in 0..k as usize {
        value *= BigUint::from(n - index);
        value /= BigUint::from(index + 1);
    }
    value
}

fn log2_big(value: &BigUint) -> f64 {
    let bits = value.bits();
    let shift = bits.saturating_sub(53);
    let top = (value >> shift).to_u64().unwrap() as f64;
    top.log2() + shift as f64
}

fn target_diagnostic(channel_columns: usize, order: &BigUint) -> TargetDiagnostic {
    let left_columns = channel_columns / 2;
    let right_columns = channel_columns - left_columns;
    let left = choose(left_columns, HALF_ARITY) << HALF_ARITY as usize;
    let right = choose(right_columns, HALF_ARITY) << HALF_ARITY as usize;
    let domain = &left * &right;
    let entries = BigUint::from(2u8) * (&left + &right);
    let storage = &entries * BigUint::from(ENTRY_BYTES);
    let disk = &storage << 1usize;
    let log_order = log2_big(order);
    let log_entries = log2_big(&entries);
    let stage_s = 2.0f64.powf(log_entries - log_order / 2.0);
    TargetDiagnostic {
        channel_columns,
        atom_arity_per_channel: ATOM_ARITY,
        original_point_terms_upper_bound: 4 * ATOM_ARITY,
        deterministic_partition: "sorted channel indices split floor(N/2)|ceil(N/2); exactly 10 signed atoms from each half".into(),
        left_columns,
        right_columns,
        half_arity: HALF_ARITY,
        left_list_entries: left.to_string(),
        right_list_entries: right.to_string(),
        left_list_log2: log2_big(&left),
        right_list_log2: log2_big(&right),
        one_target_domain: domain.to_string(),
        one_target_domain_log2: log2_big(&domain),
        one_target_mean_log2: log2_big(&domain) - log_order,
        two_target_mean_log2: 2.0 * (log2_big(&domain) - log_order),
        two_target_total_entries: entries.to_string(),
        two_target_total_entries_log2: log_entries,
        optimistic_stage_s: stage_s,
        optimistic_stage_ratio_to_rho: stage_s / RHO_S,
        optimistic_stage_ratio_to_rho_log2: (stage_s / RHO_S).log2(),
        materialized_bytes: storage.to_string(),
        materialized_bytes_log2: log2_big(&storage),
        disk_write_plus_read_bytes: disk.to_string(),
        disk_write_plus_read_log2: log2_big(&disk),
        assumes_atom_logs_available: true,
        is_s17_relation: false,
        is_end_to_end_measurement: false,
    }
}

fn gamma_ratio(collisions: u64) -> f64 {
    let mut ratio = std::f64::consts::PI.sqrt() / 2.0;
    for k in 1..collisions {
        ratio *= (k as f64 + 0.5) / k as f64;
    }
    ratio
}

fn log2_sum(values: &[f64]) -> f64 {
    let maximum = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    maximum
        + values
            .iter()
            .map(|value| 2.0f64.powf(value - maximum))
            .sum::<f64>()
            .log2()
}

fn collection_projection(
    variables: usize,
    target: &TargetDiagnostic,
    order: &BigUint,
) -> CollectionProjection {
    let larger = REQUESTED_ROWS.div_ceil(2);
    let smaller = REQUESTED_ROWS / 2;
    let free_s = 2.0f64.sqrt() * gamma_ratio(larger);
    let charged_s = 2.0 * free_s;
    let log_order = log2_big(order);
    let collection_log = log_order / 2.0 + charged_s.log2();
    let nnz = BigUint::from(REQUESTED_ROWS) * BigUint::from(ATOM_ARITY);
    let k = BigUint::from(variables);
    let sparse = BigUint::from(2u8) * &k * &nnz + &k * &k;
    let sparse_log = log2_big(&sparse);
    let complete_log = log2_sum(&[
        collection_log,
        sparse_log,
        target.two_target_total_entries_log2,
    ]);
    let complete_s = 2.0f64.powf(complete_log - log_order / 2.0);
    CollectionProjection {
        requested_usable_rows: REQUESTED_ROWS,
        larger_channel_rows: larger,
        smaller_channel_rows: smaller,
        atom_log_variables: variables,
        nonzeros_per_atom_row: ATOM_ARITY,
        exact_kth_collision_formula: "T_K=sqrt(2*n)*Gamma(K+1/2)/Gamma(K); K=69016; disjointness and simultaneous-channel stopping penalties omitted".into(),
        optimistic_free_joint_sample_s: free_s,
        optimistic_free_joint_sample_ratio_to_rho: free_s / RHO_S,
        charged_two_group_sum_s: charged_s,
        charged_two_group_sum_ratio_to_rho: charged_s / RHO_S,
        charged_collection_operations_log2: collection_log,
        charged_collection_per_usable_row_log2: collection_log
            - (REQUESTED_ROWS as f64).log2(),
        sparse_linear_algebra_operations: sparse.to_string(),
        sparse_linear_algebra_operations_log2: sparse_log,
        target_decomposition_operations_log2: target.two_target_total_entries_log2,
        optimistic_complete_operations_log2: complete_log,
        optimistic_complete_s: complete_s,
        optimistic_complete_ratio_to_rho: complete_s / RHO_S,
        memoryless_state_bytes: 512,
        materialized_target_storage_log2: target.materialized_bytes_log2,
        projection_is_measurement: false,
        projection_is_exhaustive_search: false,
        omitted_favourable_costs: vec![
            "extra stopping time for both collision channels to reach quota".into(),
            "coefficient generation and distinct-support rejection".into(),
            "relation verification and duplicate/rank loss beyond the requested row allowance".into(),
            "structured solver, setup, recovery, and exact sparse-field constants".into(),
        ],
    }
}

fn boundary_table(
    target: &TargetDiagnostic,
    collection: &CollectionProjection,
    imported: &ImportedBoundary,
) -> Vec<BoundaryRow> {
    vec![
        BoundaryRow {
            variant: "Pollard rho".into(),
            s_group_addition_equivalents: RHO_S,
            ratio_to_rho: 1.0,
            correctness: "repository reference".into(),
            evidence_kind: "fixed reference".into(),
            class: "reference".into(),
        },
        BoundaryRow {
            variant: "Round 294 registered-base free collision oracle".into(),
            s_group_addition_equivalents: RHO_S
                * imported.round294_registered_base_exact_kth_ratio_to_rho,
            ratio_to_rho: imported.round294_registered_base_exact_kth_ratio_to_rho,
            correctness: "hash-pinned imported boundary".into(),
            evidence_kind: "optimistic lower projection".into(),
            class: "boundary".into(),
        },
        BoundaryRow {
            variant: "atom split: two target decompositions only".into(),
            s_group_addition_equivalents: target.optimistic_stage_s,
            ratio_to_rho: target.optimistic_stage_ratio_to_rho,
            correctness: "factorization replayed; decompositions projected".into(),
            evidence_kind: "target-only stage projection".into(),
            class: "accounting".into(),
        },
        BoundaryRow {
            variant: "atom collection: free paired-output samples".into(),
            s_group_addition_equivalents: collection.optimistic_free_joint_sample_s,
            ratio_to_rho: collection.optimistic_free_joint_sample_ratio_to_rho,
            correctness: "favourable lower boundary".into(),
            evidence_kind: "non-exhaustive projection".into(),
            class: "boundary".into(),
        },
        BoundaryRow {
            variant: "atom collection: two charged group sums per sample".into(),
            s_group_addition_equivalents: collection.charged_two_group_sum_s,
            ratio_to_rho: collection.charged_two_group_sum_ratio_to_rho,
            correctness: "favourable lower boundary".into(),
            evidence_kind: "non-exhaustive projection".into(),
            class: "accounting".into(),
        },
        BoundaryRow {
            variant: "atom optimistic cold total".into(),
            s_group_addition_equivalents: collection.optimistic_complete_s,
            ratio_to_rho: collection.optimistic_complete_ratio_to_rho,
            correctness: "no recovered DLP; omitted costs listed".into(),
            evidence_kind: "lower projection, not measurement".into(),
            class: "negative".into(),
        },
    ]
}

fn build_result(cli: &Cli) -> Result<ResultReceipt, String> {
    let (dependencies, imported, dump) = dependencies(cli)?;
    let curve = CurveParams::p256();
    let (atoms, atom_base) = atom_census(&dump, &cli.atom_base_out, &curve)?;
    let toy = toy_census()?;
    let target_control = p256_target_control(&curve);
    if target_control.split_replay_failures != 0 || !target_control.all_points_on_curve {
        return Err("P-256 target split control failed".into());
    }
    let channel_columns = usize::min(
        atoms.plus_unique_up_to_negation,
        atoms.minus_unique_up_to_negation,
    );
    let effective_variables = if atoms.complete_small_scalar_cross_transport {
        usize::max(
            atoms.plus_unique_up_to_negation,
            atoms.minus_unique_up_to_negation,
        )
    } else {
        atom_base.columns
    };
    let target = target_diagnostic(channel_columns, &curve.n);
    let collection = collection_projection(effective_variables, &target, &curve.n);
    let table = boundary_table(&target, &collection, &imported);
    let dimension_reduction = atoms.independent_log_dimension_reduction > 0;
    let target_pass = target.optimistic_stage_ratio_to_rho <= 1.0;
    let collection_pass = collection.optimistic_complete_operations_log2 < 120.0;
    let per_row_pass = collection.charged_collection_per_usable_row_log2 < 103.0;
    let storage_pass = target.materialized_bytes_log2 < 50.0;
    let parity_pass = collection.optimistic_complete_ratio_to_rho <= 1.0;
    let gates = Gates {
        dependency_hashes_schemas_and_curve_checked: true,
        complete_actual_half_turn_census: true,
        zero_actual_replay_and_coordinate_failures: true,
        complete_toy_census: true,
        zero_toy_false_positives_and_false_negatives: true,
        actual_nonlabelled_dimension_reduction: dimension_reduction,
        complete_known_scalar_cross_transport: atoms.complete_small_scalar_cross_transport,
        target_stage_at_or_below_rho: target_pass,
        structured_residual_degree_at_most_5: false,
        complete_relation_decomposition_and_recovery: false,
        complete_collection_below_2_120: collection_pass,
        per_usable_relation_below_2_103: per_row_pass,
        peak_materialized_storage_below_2_50: storage_pass,
        parity_at_or_below_rho: parity_pass,
        discarded_probabilistic_branches_counted_as_exhaustive: false,
        promoted: false,
    };
    let classification = "involution-eigenspace/target-factorization/invertible-variable-relabeling/relation-boundary-blocked";
    let obstruction = format!(
        "The half-turn product target separates into two atom targets, but the A/B transform has aggregate rank {} and reduces {} original variables by only {}. The target-only stage is {:.6} times rho and is not S17; collecting {} rows has a favourable charged lower projection of 2^{:.6} operations ({:.6} times rho), before omitted solver and recovery costs. Structured degree remains unset.",
        atoms.aggregate_transform_rank,
        atoms.parent_columns,
        atoms.independent_log_dimension_reduction,
        target.optimistic_stage_ratio_to_rho,
        REQUESTED_ROWS,
        collection.charged_collection_operations_log2,
        collection.charged_two_group_sum_ratio_to_rho,
    );
    let decision = "Reject involution atomization as an end-to-end rho-parity construction unless an independently verified atom collision or non-generic collector reduces the unknown-log dimension and full collection cost. Preserve the exact target factorization as an accounting correction: it lowers a separable target stage for a different arity-20-plus-20 family, not the S17 source or complete DLP cost.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": &dependencies,
        "atoms": &atoms,
        "atom_base": &atom_base,
        "toy": &toy,
        "target_control": &target_control,
        "target": &target,
        "collection": &collection,
        "imported": &imported,
        "table": &table,
        "gates": &gates,
        "classification": classification,
        "obstruction": obstruction,
        "decision": decision,
    });
    Ok(ResultReceipt {
        schema: "p256-involution-atom-screen/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 311,
        execution_status: "complete".into(),
        registered_factor_base: REGISTERED_FB.into(),
        half_turn_factor_base: HALF_TURN_FB.into(),
        dependencies,
        construction: "For every half-turn cycle, A=P+P_pi and B=P-P_pi; paired target (Q,G) maps to ((Q+G)/2,(Q-G)/2).".into(),
        correction_scope: "Round 309's fixed S17 domain remains exact. Its generic order-n product-search interpretation does not apply to a direct-product eigenspace coefficient family, whose targets separate; that family uses two arity-20 atom equations and is not S17.".into(),
        atom_census: atoms,
        atom_factor_base: atom_base,
        exhaustive_toy_census: toy,
        p256_target_control: target_control,
        target_diagnostic: target,
        relation_collection_projection: collection,
        imported_boundary: imported,
        one_unit_boundary_table: table,
        structured_residual_degree_of_regularity: None,
        relations_reported_on_actual_factor_base: 0,
        full_depth_unplanted_p256_relation_attempted: false,
        gates,
        classification: classification.into(),
        dominant_obstruction: obstruction,
        decision: decision.into(),
        semantic_evidence_sha256: sha256_hex(&serde_json::to_vec(&semantic).unwrap()),
        result_json_bytes: 0,
    })
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

fn run(cli: Cli) -> Result<(), String> {
    let mut result = build_result(&cli)?;
    write_result(&cli.out, &mut result)?;
    eprintln!(
        "round 311: cycles={}, atom_columns={}, reduction={}, target_ratio={:.6}, complete_ratio={:.6}, promoted={}",
        result.atom_census.two_cycles,
        result.atom_factor_base.columns,
        result.atom_census.independent_log_dimension_reduction,
        result.target_diagnostic.optimistic_stage_ratio_to_rho,
        result.relation_collection_projection.optimistic_complete_ratio_to_rho,
        result.gates.promoted,
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_involution_atom: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_toy_factorization_passes() {
        let census = toy_census().expect("toy census");
        assert_eq!(census.total_states, 960_400);
        assert_eq!(census.factorization_discrepancies, 0);
        assert_eq!(census.augmented_rank_discrepancies, 0);
    }

    #[test]
    fn target_split_replays_on_p256() {
        let control = p256_target_control(&CurveParams::p256());
        assert_eq!(control.split_replay_failures, 0);
        assert!(control.all_points_on_curve);
    }

    #[test]
    fn fixed_partition_domain_is_deterministic() {
        let curve = CurveParams::p256();
        let diagnostic = target_diagnostic(EXPECTED_CYCLES, &curve.n);
        assert_eq!(diagnostic.left_columns, 16_431);
        assert_eq!(diagnostic.right_columns, 16_431);
        assert_eq!(diagnostic.original_point_terms_upper_bound, 80);
        assert!(!diagnostic.is_s17_relation);
    }
}
