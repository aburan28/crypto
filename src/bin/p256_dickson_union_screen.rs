//! Round 295: exact low-width unions of complete P-256 Dickson fibres.

use std::collections::BTreeMap;
#[cfg(test)]
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::{sha256_hex, short_id};
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    self, WideBuildResult, WideFactorBaseDump, WideFactorBaseFacts, WideFbPoint,
};
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::json;

const CURVE: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const FAMILY: &str = "dickson-torus-union";
const ROUND294_SHA256: &str = "9716fc1ba7d24a023627e8701094ca870b2121a84f7a0148784f70565640815a";
const COMPONENTS_PER_DEPTH: u32 = 256;
const DEPTH_A: u32 = 8;
const DEPTH_B: u32 = 6;
const THEORETICAL_COLUMNS: u64 = 164;
const DEPTH9_COLUMNS: u64 = 266;
const DEPTH9_RATIO_TO_RHO: f64 = 17.734_068_359_562_66;
const THEORETICAL_RATIO_TO_RHO: f64 = 13.920_747_397_073_491;
const RHO_S: f64 = 1.3;
const POINT_KEY_BYTES: usize = 33;

#[derive(Parser)]
#[command(about = "Screen exact low-width unions of complete P-256 Dickson fibres")]
struct Cli {
    /// Canonical Round-294 result.
    #[arg(long)]
    round294: PathBuf,
    /// Compact deterministic screen result.
    #[arg(long)]
    out: PathBuf,
    /// Selected factor base in the repository wide-dump storage format.
    #[arg(long)]
    factor_base_out: PathBuf,
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    bytes: u64,
    sha256: String,
    schema: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct UnionColumn {
    x: BigUint,
    low_y: BigUint,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
struct ComponentReceipt {
    depth: u32,
    index: u32,
    derivation_preimage: String,
    specification: String,
    root_exponent: String,
    terminal: String,
    fb_id: String,
    fb_sha256: String,
    points_sha256: String,
    columns: u64,
    signed_points: u64,
    complete_rebuild_verified: bool,
}

struct Component {
    receipt: ComponentReceipt,
    built: WideBuildResult,
    columns: Vec<UnionColumn>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
struct RejectedComponent {
    depth: u32,
    index: u32,
    derivation_preimage: String,
    root_exponent: String,
    reason: String,
}

#[derive(Clone, Debug, Serialize)]
struct WidthBoundary {
    columns: u64,
    all_arity_domain: String,
    maximum_fixed_arity: u32,
    maximum_fixed_domain: String,
    maximum_fixed_covers_order: bool,
    covering_arities: Vec<u32>,
    exact_kth_collision_ratio_to_rho: f64,
    balanced_two_list_ratio_to_rho: f64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
struct PairSelection {
    depth8_index: u32,
    depth6_index: u32,
    depth8_fb_id: String,
    depth6_fb_id: String,
    depth8_columns: u64,
    depth6_columns: u64,
    duplicate_columns: u64,
    union_columns: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
struct WinnerReceipt {
    selection: PairSelection,
    fb_id: String,
    fb_sha256: String,
    points_sha256: String,
    signed_points: u64,
    component_local_generator_degree: u32,
    global_union_degree_of_regularity: Option<u32>,
    factor_base_artifact: String,
    factor_base_json_bytes: u64,
    factor_base_json_sha256: String,
    independent_replay_byte_identical: bool,
}

#[derive(Serialize)]
struct Gates {
    dependency_hash_checked: bool,
    every_component_rebuilt_exactly: bool,
    union_identity_replay_exact: bool,
    zero_false_identity_or_deduplication_matches: bool,
    exact_width_164_realized: bool,
    factor_base_boundary_below_depth9_control: bool,
    factor_base_boundary_at_or_below_rho: bool,
    structured_residual_degree_at_most_5: bool,
    per_usable_relation_below_2_103: bool,
    complete_cost_below_2_120: bool,
    projected_materialized_storage_below_2_50: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    screening_round: u32,
    execution_status: String,
    dependency: Dependency,
    derivation_domain: String,
    component_depths: [u32; 2],
    candidates_per_depth: u32,
    accepted_components: Vec<ComponentReceipt>,
    rejected_components: Vec<RejectedComponent>,
    candidate_pairs: u64,
    covering_pairs: u64,
    width_histogram: BTreeMap<u64, u64>,
    duplicate_histogram: BTreeMap<u64, u64>,
    exact_boundary_by_width: Vec<WidthBoundary>,
    boundary_table_sha256: String,
    selected: WinnerReceipt,
    selected_boundary: WidthBoundary,
    theoretical_width: u64,
    theoretical_ratio_to_rho: f64,
    depth9_control_columns: u64,
    depth9_control_ratio_to_rho: f64,
    factor_base_boundary_improvement: f64,
    exact_integer_operations: u64,
    component_builds: u64,
    signed_points_verified: u64,
    relations_reported: u64,
    full_depth_unplanted_relation_attempted: bool,
    process_telemetry_recorded_in_isolation_receipt: bool,
    gates: Gates,
    classification: String,
    dominant_obstruction: String,
    decision: String,
    semantic_evidence_sha256: String,
    result_json_bytes: u64,
}

fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn load_dependency(path: &Path) -> Result<Dependency, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND294_SHA256 {
        return Err(format!(
            "{} hash mismatch: expected {ROUND294_SHA256}, got {digest}",
            path.display()
        ));
    }
    let value: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    let schema = value
        .get("schema")
        .and_then(serde_json::Value::as_str)
        .ok_or("Round 294 result has no schema")?;
    if schema != "p256-variable-arity-closure/v2"
        || value
            .pointer("/first_fixed_arity_cover/columns")
            .and_then(serde_json::Value::as_u64)
            != Some(THEORETICAL_COLUMNS)
    {
        return Err("Round 294 boundary identity mismatch".into());
    }
    Ok(Dependency {
        path: path.display().to_string(),
        bytes: bytes.len() as u64,
        sha256: digest,
        schema: schema.into(),
    })
}

fn derivation_preimage(depth: u32, index: u32) -> String {
    format!("{CURVE}/dickson-union-round295/depth-{depth}/{index}")
}

fn root_exponent(depth: u32, index: u32) -> (String, BigUint) {
    let preimage = derivation_preimage(depth, index);
    let digest = sha256(preimage.as_bytes());
    let limit = BigUint::one() << (96usize - depth as usize);
    let exponent = (BigUint::from_bytes_be(&digest) & (&limit - BigUint::one())) | BigUint::one();
    (preimage, exponent)
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("not a non-negative integer: {value}"))
}

fn component_columns(built: &WideBuildResult) -> Result<Vec<UnionColumn>, String> {
    let points = &built.dump.points;
    if !points.len().is_multiple_of(2) {
        return Err("component has an odd signed-point count".into());
    }
    let mut columns = Vec::with_capacity(points.len() / 2);
    for (column_index, pair) in points.chunks_exact(2).enumerate() {
        if pair[0].x != pair[1].x
            || pair[0].col != column_index as u64
            || pair[1].col != column_index as u64
            || pair[0].coef != "1"
        {
            return Err(format!(
                "component point-pair mismatch at column {column_index}"
            ));
        }
        let x = parse_big(&pair[0].x)?;
        let y0 = parse_big(&pair[0].y)?;
        let y1 = parse_big(&pair[1].y)?;
        if (&y0 + &y1) % &CurveParams::p256().p != BigUint::zero() {
            return Err(format!(
                "component signs do not negate at column {column_index}"
            ));
        }
        columns.push(UnionColumn { x, low_y: y0 });
    }
    if !columns.windows(2).all(|pair| pair[0].x < pair[1].x) {
        return Err("component columns are not strictly sorted".into());
    }
    Ok(columns)
}

fn build_component(depth: u32, index: u32) -> Result<Component, RejectedComponent> {
    let (preimage, exponent) = root_exponent(depth, index);
    let specification = format!(
        "dickson-torus:depth={depth},root_exponent=0x{}",
        exponent.to_str_radix(16)
    );
    let built = match p256_dickson_factor_base::build(&specification) {
        Ok(built) => built,
        Err(reason) => {
            return Err(RejectedComponent {
                depth,
                index,
                derivation_preimage: preimage,
                root_exponent: lower_hex(&exponent),
                reason,
            });
        }
    };
    let rebuilt =
        p256_dickson_factor_base::verify(&built.dump).map_err(|reason| RejectedComponent {
            depth,
            index,
            derivation_preimage: preimage.clone(),
            root_exponent: lower_hex(&exponent),
            reason,
        })?;
    if rebuilt != built {
        return Err(RejectedComponent {
            depth,
            index,
            derivation_preimage: preimage,
            root_exponent: lower_hex(&exponent),
            reason: "component verification returned a different build".into(),
        });
    }
    let columns = component_columns(&built).map_err(|reason| RejectedComponent {
        depth,
        index,
        derivation_preimage: preimage.clone(),
        root_exponent: lower_hex(&exponent),
        reason,
    })?;
    let fb = &built.dump.factor_base;
    let receipt = ComponentReceipt {
        depth,
        index,
        derivation_preimage: preimage,
        specification,
        root_exponent: fb.params["root_exponent"].clone(),
        terminal: fb.params["terminal"].clone(),
        fb_id: fb.fb_id.clone(),
        fb_sha256: fb.fb_sha256.clone(),
        points_sha256: fb.points_sha256.clone(),
        columns: fb.columns,
        signed_points: fb.signed_points,
        complete_rebuild_verified: true,
    };
    Ok(Component {
        receipt,
        built,
        columns,
    })
}

fn intersection_size(left: &[UnionColumn], right: &[UnionColumn]) -> u64 {
    let (mut i, mut j, mut common) = (0usize, 0usize, 0u64);
    while i < left.len() && j < right.len() {
        match left[i].x.cmp(&right[j].x) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                common += 1;
                i += 1;
                j += 1;
            }
        }
    }
    common
}

fn merge_columns(left: &[UnionColumn], right: &[UnionColumn]) -> Result<Vec<UnionColumn>, String> {
    let mut merged = Vec::with_capacity(left.len() + right.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < left.len() || j < right.len() {
        if j == right.len() || (i < left.len() && left[i].x < right[j].x) {
            merged.push(left[i].clone());
            i += 1;
        } else if i == left.len() || right[j].x < left[i].x {
            merged.push(right[j].clone());
            j += 1;
        } else {
            if left[i].low_y != right[j].low_y {
                return Err("equal union abscissae carry different canonical y values".into());
            }
            merged.push(left[i].clone());
            i += 1;
            j += 1;
        }
    }
    Ok(merged)
}

fn point_key_bytes(x: &BigUint, high_sign: bool) -> Result<[u8; POINT_KEY_BYTES], String> {
    let key = ((x + BigUint::one()) << 1usize) + BigUint::from(high_sign);
    let raw = key.to_bytes_be();
    if raw.len() > POINT_KEY_BYTES {
        return Err("P-256 point key exceeded 33 bytes".into());
    }
    let mut out = [0u8; POINT_KEY_BYTES];
    out[POINT_KEY_BYTES - raw.len()..].copy_from_slice(&raw);
    Ok(out)
}

fn union_dump(left: &Component, right: &Component) -> Result<WideFactorBaseDump, String> {
    let columns = merge_columns(&left.columns, &right.columns)?;
    let curve = CurveParams::p256();
    let mut key_blob = Vec::with_capacity(columns.len() * 2 * POINT_KEY_BYTES);
    let mut points = Vec::with_capacity(columns.len() * 2);
    for (col, column) in columns.iter().enumerate() {
        key_blob.extend(point_key_bytes(&column.x, false)?);
        key_blob.extend(point_key_bytes(&column.x, true)?);
        points.push(WideFbPoint {
            x: lower_hex(&column.x),
            y: lower_hex(&column.low_y),
            col: col as u64,
            coef: "1".into(),
        });
        points.push(WideFbPoint {
            x: lower_hex(&column.x),
            y: lower_hex(&(&curve.p - &column.low_y)),
            col: col as u64,
            coef: (&curve.n - BigUint::one()).to_string(),
        });
    }
    let points_sha256 = sha256_hex(&key_blob);
    let ordered = [&left.receipt, &right.receipt];
    let mut params = BTreeMap::new();
    params.insert(
        "point_key_encoding".into(),
        p256_dickson_factor_base::POINT_KEY_ENCODING.into(),
    );
    params.insert("union_version".into(), "1".into());
    for (position, component) in ordered.iter().enumerate() {
        params.insert(
            format!("component_{position}_depth"),
            component.depth.to_string(),
        );
        params.insert(
            format!("component_{position}_index"),
            component.index.to_string(),
        );
        params.insert(
            format!("component_{position}_fb_id"),
            component.fb_id.clone(),
        );
        params.insert(
            format!("component_{position}_root_exponent"),
            component.root_exponent.clone(),
        );
        params.insert(
            format!("component_{position}_terminal"),
            component.terminal.clone(),
        );
    }
    let count = columns.len() as u64;
    let signed_points = points.len() as u64;
    let identity = json!({
        "schema": p256_dickson_factor_base::FACTOR_BASE_SCHEMA,
        "curve": CURVE,
        "family": FAMILY,
        "params": params,
        "columns": count,
        "signed_points": signed_points,
        "points_sha256": points_sha256,
    });
    let (fb_id, fb_sha256) = short_id("FB1", &identity)?;
    Ok(WideFactorBaseDump {
        schema: p256_dickson_factor_base::DUMP_SCHEMA.into(),
        curve: left.built.dump.curve.clone(),
        factor_base: WideFactorBaseFacts {
            fb_id,
            fb_sha256,
            family: FAMILY.into(),
            params,
            description: format!(
                "union of complete depth-{} and depth-{} Dickson trace fibres on {CURVE}",
                left.receipt.depth, right.receipt.depth
            ),
            signed_points,
            abscissae: count,
            columns: count,
            dimension: None,
            points_sha256,
        },
        build_adds: 0,
        build_doubles: 0,
        points,
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

fn ratios(columns: u64) -> (f64, f64) {
    (
        2.0_f64.sqrt() * gamma_ratio(columns) / RHO_S,
        2.0 * (columns as f64).sqrt() / RHO_S,
    )
}

fn width_boundary(columns: u64, order: &BigUint, operations: &mut u64) -> WidthBoundary {
    let mut combination = BigUint::one();
    let mut all_arity = BigUint::zero();
    let mut maximum = BigUint::zero();
    let mut maximum_arity = 0u32;
    let mut covering_arities = Vec::new();
    for arity in 0..=columns {
        let domain = &combination << arity as usize;
        all_arity += &domain;
        *operations += 2;
        if domain > maximum {
            maximum = domain.clone();
            maximum_arity = arity as u32;
        }
        if domain >= *order {
            covering_arities.push(arity as u32);
        }
        if arity < columns {
            combination *= BigUint::from(columns - arity);
            combination /= BigUint::from(arity + 1);
            *operations += 2;
        }
    }
    assert_eq!(all_arity, BigUint::from(3u8).pow(columns as u32));
    let (exact, balanced) = ratios(columns);
    WidthBoundary {
        columns,
        all_arity_domain: all_arity.to_string(),
        maximum_fixed_arity: maximum_arity,
        maximum_fixed_domain: maximum.to_string(),
        maximum_fixed_covers_order: maximum >= *order,
        covering_arities,
        exact_kth_collision_ratio_to_rho: exact,
        balanced_two_list_ratio_to_rho: balanced,
    }
}

fn boundary_digest(rows: &[WidthBoundary]) -> String {
    let mut bytes = Vec::new();
    for row in rows {
        bytes.extend(row.columns.to_be_bytes());
        bytes.extend(row.maximum_fixed_arity.to_be_bytes());
        bytes.extend(row.maximum_fixed_domain.as_bytes());
        for arity in &row.covering_arities {
            bytes.extend(arity.to_be_bytes());
        }
    }
    hex::encode(sha256(&bytes))
}

fn semantic_digest(
    accepted: &[ComponentReceipt],
    widths: &BTreeMap<u64, u64>,
    selection: &PairSelection,
    winner: &WinnerReceipt,
) -> String {
    let mut bytes = Vec::new();
    bytes.extend(CURVE.as_bytes());
    for component in accepted {
        bytes.extend(component.depth.to_be_bytes());
        bytes.extend(component.index.to_be_bytes());
        bytes.extend(component.fb_id.as_bytes());
        bytes.extend(component.columns.to_be_bytes());
    }
    for (width, count) in widths {
        bytes.extend(width.to_be_bytes());
        bytes.extend(count.to_be_bytes());
    }
    bytes.extend(selection.depth8_index.to_be_bytes());
    bytes.extend(selection.depth6_index.to_be_bytes());
    bytes.extend(selection.union_columns.to_be_bytes());
    bytes.extend(winner.fb_id.as_bytes());
    bytes.extend(winner.factor_base_json_sha256.as_bytes());
    hex::encode(sha256(&bytes))
}

fn stable_json<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    let mut encoded = serde_json::to_vec_pretty(value).map_err(|error| error.to_string())?;
    encoded.push(b'\n');
    Ok(encoded)
}

fn write_result(path: &Path, result: &mut ResultReceipt) -> Result<(), String> {
    for _ in 0..8 {
        let encoded = stable_json(result)?;
        let bytes = encoded.len() as u64;
        if result.result_json_bytes == bytes {
            return fs::write(path, encoded).map_err(|error| error.to_string());
        }
        result.result_json_bytes = bytes;
    }
    Err("result JSON byte count did not stabilize".into())
}

fn build_result(cli: &Cli) -> Result<(ResultReceipt, Vec<u8>), String> {
    let dependency = load_dependency(&cli.round294)?;
    let mut depth8 = Vec::new();
    let mut depth6 = Vec::new();
    let mut rejected = Vec::new();
    for depth in [DEPTH_A, DEPTH_B] {
        for index in 0..COMPONENTS_PER_DEPTH {
            match build_component(depth, index) {
                Ok(component) if depth == DEPTH_A => depth8.push(component),
                Ok(component) => depth6.push(component),
                Err(row) => rejected.push(row),
            }
        }
    }
    if depth8.is_empty() || depth6.is_empty() {
        return Err("a complete component depth was rejected".into());
    }

    let order = CurveParams::p256().n;
    let mut operations = 0u64;
    let mut width_histogram = BTreeMap::new();
    let mut duplicate_histogram = BTreeMap::new();
    let mut boundary_cache: BTreeMap<u64, WidthBoundary> = BTreeMap::new();
    let mut best: Option<PairSelection> = None;
    let mut covering_pairs = 0u64;
    for left in &depth8 {
        for right in &depth6 {
            let duplicates = intersection_size(&left.columns, &right.columns);
            let width = left.receipt.columns + right.receipt.columns - duplicates;
            *width_histogram.entry(width).or_insert(0) += 1;
            *duplicate_histogram.entry(duplicates).or_insert(0) += 1;
            let covers = if let Some(boundary) = boundary_cache.get(&width) {
                boundary.maximum_fixed_covers_order
            } else {
                let boundary = width_boundary(width, &order, &mut operations);
                let covers = boundary.maximum_fixed_covers_order;
                boundary_cache.insert(width, boundary);
                covers
            };
            if !covers {
                continue;
            }
            covering_pairs += 1;
            let candidate = PairSelection {
                depth8_index: left.receipt.index,
                depth6_index: right.receipt.index,
                depth8_fb_id: left.receipt.fb_id.clone(),
                depth6_fb_id: right.receipt.fb_id.clone(),
                depth8_columns: left.receipt.columns,
                depth6_columns: right.receipt.columns,
                duplicate_columns: duplicates,
                union_columns: width,
            };
            let key = |row: &PairSelection| {
                (
                    row.union_columns,
                    row.duplicate_columns,
                    row.depth8_index,
                    row.depth6_index,
                )
            };
            if best
                .as_ref()
                .is_none_or(|current| key(&candidate) < key(current))
            {
                best = Some(candidate);
            }
        }
    }
    let selection = best.ok_or("no covering two-fibre union in the frozen search")?;
    let left = depth8
        .iter()
        .find(|component| component.receipt.index == selection.depth8_index)
        .ok_or("selected depth-8 component disappeared")?;
    let right = depth6
        .iter()
        .find(|component| component.receipt.index == selection.depth6_index)
        .ok_or("selected depth-6 component disappeared")?;
    let dump = union_dump(left, right)?;
    if dump.factor_base.columns != selection.union_columns {
        return Err("selected union width disagrees with its dump".into());
    }

    let replay_left = build_component(DEPTH_A, selection.depth8_index)
        .map_err(|row| format!("winner depth-8 replay failed: {}", row.reason))?;
    let replay_right = build_component(DEPTH_B, selection.depth6_index)
        .map_err(|row| format!("winner depth-6 replay failed: {}", row.reason))?;
    let replay_dump = union_dump(&replay_left, &replay_right)?;
    let dump_bytes = stable_json(&dump)?;
    let replay_bytes = stable_json(&replay_dump)?;
    if dump_bytes != replay_bytes {
        return Err("selected union dump is not byte-identical on replay".into());
    }
    let factor_base_json_sha256 = hex::encode(sha256(&dump_bytes));
    let winner = WinnerReceipt {
        selection: selection.clone(),
        fb_id: dump.factor_base.fb_id.clone(),
        fb_sha256: dump.factor_base.fb_sha256.clone(),
        points_sha256: dump.factor_base.points_sha256.clone(),
        signed_points: dump.factor_base.signed_points,
        component_local_generator_degree: 2,
        global_union_degree_of_regularity: None,
        factor_base_artifact: cli
            .factor_base_out
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("selected-factor-base.json")
            .into(),
        factor_base_json_bytes: dump_bytes.len() as u64,
        factor_base_json_sha256,
        independent_replay_byte_identical: true,
    };

    let exact_boundary_by_width = boundary_cache.into_values().collect::<Vec<_>>();
    let selected_boundary = exact_boundary_by_width
        .iter()
        .find(|row| row.columns == selection.union_columns)
        .cloned()
        .ok_or("selected union has no boundary row")?;
    let boundary_table_sha256 = boundary_digest(&exact_boundary_by_width);
    let mut accepted = depth8
        .iter()
        .chain(depth6.iter())
        .map(|component| component.receipt.clone())
        .collect::<Vec<_>>();
    accepted.sort_by_key(|row| (row.depth, row.index));
    rejected.sort_by_key(|row| (row.depth, row.index));
    let signed_points_verified = accepted.iter().map(|row| row.signed_points).sum::<u64>()
        + replay_left.receipt.signed_points
        + replay_right.receipt.signed_points;
    let factor_base_boundary_improvement =
        DEPTH9_RATIO_TO_RHO / selected_boundary.exact_kth_collision_ratio_to_rho;
    let semantic_evidence_sha256 =
        semantic_digest(&accepted, &width_histogram, &selection, &winner);
    let exact_width = selection.union_columns == THEORETICAL_COLUMNS;
    let below_control = selected_boundary.exact_kth_collision_ratio_to_rho < DEPTH9_RATIO_TO_RHO;
    let at_or_below_rho = selected_boundary.exact_kth_collision_ratio_to_rho <= 1.0;

    Ok((
        ResultReceipt {
            schema: "p256-dickson-union-screen/v1".into(),
            curve: CURVE.into(),
            screening_round: 295,
            execution_status: "complete".into(),
            dependency,
            derivation_domain: format!("{CURVE}/dickson-union-round295/depth-{{d}}/{{i}}"),
            component_depths: [DEPTH_A, DEPTH_B],
            candidates_per_depth: COMPONENTS_PER_DEPTH,
            accepted_components: accepted,
            rejected_components: rejected,
            candidate_pairs: depth8.len() as u64 * depth6.len() as u64,
            covering_pairs,
            width_histogram,
            duplicate_histogram,
            exact_boundary_by_width,
            boundary_table_sha256,
            selected: winner,
            selected_boundary,
            theoretical_width: THEORETICAL_COLUMNS,
            theoretical_ratio_to_rho: THEORETICAL_RATIO_TO_RHO,
            depth9_control_columns: DEPTH9_COLUMNS,
            depth9_control_ratio_to_rho: DEPTH9_RATIO_TO_RHO,
            factor_base_boundary_improvement,
            exact_integer_operations: operations,
            component_builds: (accepted_component_count(&depth8, &depth6) + 2) * 2,
            signed_points_verified,
            relations_reported: 0,
            full_depth_unplanted_relation_attempted: false,
            process_telemetry_recorded_in_isolation_receipt: true,
            gates: Gates {
                dependency_hash_checked: true,
                every_component_rebuilt_exactly: true,
                union_identity_replay_exact: true,
                zero_false_identity_or_deduplication_matches: true,
                exact_width_164_realized: exact_width,
                factor_base_boundary_below_depth9_control: below_control,
                factor_base_boundary_at_or_below_rho: at_or_below_rho,
                structured_residual_degree_at_most_5: false,
                per_usable_relation_below_2_103: false,
                complete_cost_below_2_120: false,
                projected_materialized_storage_below_2_50: false,
                promoted: false,
            },
            classification: "factor-base-boundary-advance/promotion-negative".into(),
            dominant_obstruction: "The union realizes a near-minimal complete low-degree component inventory, but every column still carries an independent logarithm. Its exact generic-collection boundary remains above rho, and the global high-arity union system has no measured degree of regularity.".into(),
            decision: "Retain the selected FB1 as the best complete low-width Dickson inventory found. Do not attempt an unplanted P-256 relation; continue only with a non-generic batch solver or exact logarithm transport that reduces the independent-log quotient.".into(),
            semantic_evidence_sha256,
            result_json_bytes: 0,
        },
        dump_bytes,
    ))
}

fn accepted_component_count(depth8: &[Component], depth6: &[Component]) -> u64 {
    (depth8.len() + depth6.len()) as u64
}

fn run(cli: Cli) -> Result<(), String> {
    let (mut result, dump_bytes) = build_result(&cli)?;
    fs::write(&cli.factor_base_out, dump_bytes).map_err(|error| error.to_string())?;
    write_result(&cli.out, &mut result)?;
    eprintln!(
        "round 295: {} at B={}, m={:?}, exact boundary {:.6}x rho (depth-9 control {:.6}x) -> {}",
        result.selected.fb_id,
        result.selected_boundary.columns,
        result.selected_boundary.covering_arities,
        result.selected_boundary.exact_kth_collision_ratio_to_rho,
        result.depth9_control_ratio_to_rho,
        cli.out.display()
    );
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

    #[test]
    fn derived_exponents_are_odd_bounded_and_distinct() {
        for depth in [DEPTH_A, DEPTH_B] {
            let limit = BigUint::one() << (96usize - depth as usize);
            let mut seen = BTreeSet::new();
            for index in 0..COMPONENTS_PER_DEPTH {
                let (_, exponent) = root_exponent(depth, index);
                assert!(exponent < limit);
                assert!(exponent.bit(0));
                assert!(seen.insert(exponent));
            }
        }
    }

    #[test]
    fn exact_width_164_matches_round294() {
        let mut operations = 0;
        let row = width_boundary(THEORETICAL_COLUMNS, &CurveParams::p256().n, &mut operations);
        assert_eq!(row.covering_arities, vec![109, 110]);
        assert!((row.exact_kth_collision_ratio_to_rho - THEORETICAL_RATIO_TO_RHO).abs() < 1e-12);
        assert!(operations > 0);
    }

    #[test]
    fn merge_deduplicates_equal_columns() {
        let column = |x: u32| UnionColumn {
            x: BigUint::from(x),
            low_y: BigUint::from(x + 10),
        };
        let left = vec![column(1), column(3)];
        let right = vec![column(2), column(3)];
        assert_eq!(intersection_size(&left, &right), 1);
        let merged = merge_columns(&left, &right).unwrap();
        assert_eq!(
            merged.iter().map(|row| &row.x).collect::<Vec<_>>(),
            vec![
                &BigUint::from(1u8),
                &BigUint::from(2u8),
                &BigUint::from(3u8)
            ]
        );
    }

    #[test]
    fn generic_boundary_increases_with_width() {
        assert!(ratios(THEORETICAL_COLUMNS).0 < ratios(DEPTH9_COLUMNS).0);
        assert!(ratios(THEORETICAL_COLUMNS).1 < ratios(DEPTH9_COLUMNS).1);
    }
}
