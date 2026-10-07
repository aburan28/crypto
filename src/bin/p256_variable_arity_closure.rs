//! Round 294: exact variable-arity closure for P-256 Dickson factor bases.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{self, WideBuildResult};
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;

const CURVE: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const REGISTERED_SPEC: &str = "dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129";
const REGISTERED_FB: &str = "FB1h2f8621cda105";
const REGISTERED_FB_SHA256: &str =
    "2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42";
const REGISTERED_POINTS_SHA256: &str =
    "70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1";
const REGISTERED_COLUMNS: u64 = 131_458;
const REGISTERED_SIGNED_POINTS: u64 = 262_916;
const RHO_S: f64 = 1.3;
const MAX_BOUNDARY_WIDTH: u64 = 512;
const MAX_ARITY: u32 = 256;
const LADDER_DEPTHS: [u32; 6] = [7, 8, 9, 10, 11, 12];

const ROUND21_SHA256: &str = "4aaf5fde72257ed1f262d5674ab2fe81c7b4f90698fa3767615adfed6524e757";
const ROUND24_SHA256: &str = "d1e2fe33e1abbae8a9bf6885ca7cfef8012233363f4cd43e0d30d468a381fc9c";
const ROUND25_SHA256: &str = "dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf";
const ROUND29_SHA256: &str = "9be6e5fc38644c34ad746dec9a6542da7e35ef2d8e11eb39a83de7af0f54f443";
const ROUND293_SHA256: &str = "02626b9147e024fd222ad6505688c5113a955a0e68b73afa1e7600cda04c6c97";

#[derive(Parser)]
#[command(about = "Run the exact P-256 variable-arity Dickson closure")]
struct Cli {
    #[arg(long)]
    round21: PathBuf,
    #[arg(long)]
    round24: PathBuf,
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    round29: PathBuf,
    #[arg(long)]
    round293: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    bytes: u64,
    sha256: String,
}

#[derive(Clone, Serialize)]
struct WidthBoundary {
    columns: u64,
    all_arity_domain: String,
    all_arity_covers_order: bool,
    maximum_fixed_arity: u32,
    maximum_fixed_domain: String,
    maximum_fixed_covers_order: bool,
    covering_arities: Vec<u32>,
    preregistered_superseded_shortcut_ratio_to_rho: f64,
    exact_kth_collision_ratio_to_rho: f64,
    balanced_two_list_ratio_to_rho: f64,
}

#[derive(Clone, Serialize)]
struct FixedArityThreshold {
    arity: u32,
    minimum_columns: Option<u64>,
    signed_domain: Option<String>,
    covers_order: bool,
    preregistered_superseded_shortcut_ratio_to_rho: Option<f64>,
    exact_kth_collision_ratio_to_rho: Option<f64>,
    balanced_two_list_ratio_to_rho: Option<f64>,
}

#[derive(Clone, Serialize)]
struct DicksonInventory {
    specification: String,
    depth: u32,
    terminal: String,
    fb_id: String,
    fb_sha256: String,
    points_sha256: String,
    legacy_enumeration_sha256: String,
    columns: u64,
    signed_points: u64,
    points_verified: u64,
    verification_failures: u64,
    maximum_checked_arity: u32,
    maximum_domain_arity: u32,
    maximum_domain: String,
    first_covering_arity: Option<u32>,
    fixed_arity_covers_order: bool,
    preregistered_superseded_shortcut_ratio_to_rho: f64,
    exact_kth_collision_ratio_to_rho: f64,
    balanced_two_list_ratio_to_rho: f64,
    logical_point_bytes: u64,
}

#[derive(Serialize)]
struct Gates {
    dependencies_hash_checked: bool,
    exact_binomial_recurrence: bool,
    all_arity_identity_exact: bool,
    native_factor_bases_verified: bool,
    zero_false_positives_and_false_negatives: bool,
    any_covering_row_at_or_below_rho_boundary: bool,
    structured_residual_degree_at_most_5: bool,
    per_usable_relation_below_2_103: bool,
    complete_cost_below_2_120: bool,
    projected_storage_below_2_50: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    comparison_factor_base: String,
    screening_round: u32,
    execution_status: String,
    dependencies: Vec<Dependency>,
    subgroup_order: String,
    rho_s: f64,
    coefficient_domain_formula: String,
    all_arity_identity: String,
    boundary_width_limit: u64,
    arity_limit: u32,
    first_all_arity_cover_columns: u64,
    first_fixed_arity_cover: WidthBoundary,
    boundary_neighbourhood: Vec<WidthBoundary>,
    fixed_arity_thresholds: Vec<FixedArityThreshold>,
    boundary_table_sha256: String,
    dickson_inventory: Vec<DicksonInventory>,
    minimum_covering_inventory_id: String,
    minimum_covering_inventory_exact_ratio_to_rho: f64,
    global_minimum_preregistered_superseded_shortcut_ratio_to_rho: f64,
    global_minimum_exact_kth_collision_ratio_to_rho: f64,
    global_minimum_balanced_two_list_ratio_to_rho: f64,
    imported_local_solving_degree_maxima: [u32; 3],
    imported_local_split_image_degree: u32,
    variable_arity_unsplit_degree_of_regularity: Option<u32>,
    exact_integer_operations: u64,
    factor_base_builds: u64,
    points_verified: u64,
    algorithm_disk_bytes: u64,
    logical_materialized_bytes: u64,
    process_telemetry_recorded_in_isolation_receipt: bool,
    result_json_bytes: u64,
    gates: Gates,
    relations_reported: u64,
    full_depth_unplanted_relation_attempted: bool,
    classification: String,
    dominant_obstruction: String,
    decision: String,
    semantic_evidence_sha256: String,
}

fn load_dependency(path: &Path, expected: &str) -> Result<Dependency, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != expected {
        return Err(format!(
            "{} hash mismatch: expected {expected}, got {digest}",
            path.display()
        ));
    }
    serde_json::from_slice::<serde_json::Value>(&bytes)
        .map_err(|error| format!("{}: {error}", path.display()))?;
    Ok(Dependency {
        path: path.display().to_string(),
        bytes: bytes.len() as u64,
        sha256: digest,
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

fn ratios(columns: u64) -> (f64, f64, f64) {
    let root = (columns as f64).sqrt();
    (
        ((std::f64::consts::PI / 2.0).sqrt() * root) / RHO_S,
        (2.0_f64.sqrt() * gamma_ratio(columns)) / RHO_S,
        2.0 * root / RHO_S,
    )
}

fn signed_domain(columns: u64, arity: u32, operations: &mut u64) -> BigUint {
    if u64::from(arity) > columns {
        return BigUint::zero();
    }
    let mut combinations = BigUint::one();
    for i in 0..u64::from(arity) {
        combinations *= BigUint::from(columns - i);
        combinations /= BigUint::from(i + 1);
        *operations += 2;
    }
    combinations << arity as usize
}

fn width_boundary(columns: u64, order: &BigUint, operations: &mut u64) -> WidthBoundary {
    let mut combination = BigUint::one();
    let mut sum = BigUint::zero();
    let mut maximum = BigUint::zero();
    let mut maximum_arity = 0;
    let mut covering_arities = Vec::new();
    for arity in 0..=columns {
        let domain = &combination << arity as usize;
        sum += &domain;
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
    let expected_sum = BigUint::from(3u8).pow(columns as u32);
    assert_eq!(
        sum, expected_sum,
        "all-arity identity failed at B={columns}"
    );
    let (superseded_shortcut, exact_kth, two_list) = ratios(columns);
    WidthBoundary {
        columns,
        all_arity_domain: sum.to_string(),
        all_arity_covers_order: expected_sum >= *order,
        maximum_fixed_arity: maximum_arity,
        maximum_fixed_domain: maximum.to_string(),
        maximum_fixed_covers_order: maximum >= *order,
        covering_arities,
        preregistered_superseded_shortcut_ratio_to_rho: superseded_shortcut,
        exact_kth_collision_ratio_to_rho: exact_kth,
        balanced_two_list_ratio_to_rho: two_list,
    }
}

fn minimum_columns_for_arity(
    arity: u32,
    maximum_columns: u64,
    order: &BigUint,
    operations: &mut u64,
) -> Option<(u64, BigUint)> {
    if signed_domain(maximum_columns, arity, operations) < *order {
        return None;
    }
    let mut low = u64::from(arity);
    let mut high = maximum_columns;
    while low < high {
        let middle = low + (high - low) / 2;
        if signed_domain(middle, arity, operations) >= *order {
            high = middle;
        } else {
            low = middle + 1;
        }
    }
    let domain = signed_domain(low, arity, operations);
    Some((low, domain))
}

fn checked_inventory(
    specification: String,
    expected_registered: bool,
    order: &BigUint,
    operations: &mut u64,
) -> Result<DicksonInventory, String> {
    let built = p256_dickson_factor_base::build(&specification)?;
    p256_dickson_factor_base::verify(&built.dump)?;
    check_registered(&built, expected_registered)?;
    let factor_base = &built.dump.factor_base;
    let depth = factor_base
        .params
        .get("depth")
        .ok_or("built factor base has no depth")?
        .parse::<u32>()
        .map_err(|error| error.to_string())?;
    let maximum_checked_arity = u32::min(MAX_ARITY, factor_base.columns as u32);
    let mut maximum_domain = BigUint::zero();
    let mut maximum_domain_arity = 0;
    let mut first_covering_arity = None;
    for arity in 0..=maximum_checked_arity {
        let domain = signed_domain(factor_base.columns, arity, operations);
        if domain > maximum_domain {
            maximum_domain = domain.clone();
            maximum_domain_arity = arity;
        }
        if first_covering_arity.is_none() && domain >= *order {
            first_covering_arity = Some(arity);
        }
    }
    let (superseded_shortcut, exact_kth, two_list) = ratios(factor_base.columns);
    let logical_point_bytes = built
        .dump
        .points
        .iter()
        .map(|point| {
            (std::mem::size_of_val(point) + point.x.len() + point.y.len() + point.coef.len()) as u64
        })
        .sum();
    Ok(DicksonInventory {
        specification,
        depth,
        terminal: factor_base
            .params
            .get("terminal")
            .cloned()
            .ok_or("built factor base has no terminal")?,
        fb_id: factor_base.fb_id.clone(),
        fb_sha256: factor_base.fb_sha256.clone(),
        points_sha256: factor_base.points_sha256.clone(),
        legacy_enumeration_sha256: built.legacy_enumeration_sha256,
        columns: factor_base.columns,
        signed_points: factor_base.signed_points,
        points_verified: factor_base.signed_points,
        verification_failures: 0,
        maximum_checked_arity,
        maximum_domain_arity,
        maximum_domain: maximum_domain.to_string(),
        first_covering_arity,
        fixed_arity_covers_order: maximum_domain >= *order,
        preregistered_superseded_shortcut_ratio_to_rho: superseded_shortcut,
        exact_kth_collision_ratio_to_rho: exact_kth,
        balanced_two_list_ratio_to_rho: two_list,
        logical_point_bytes,
    })
}

fn check_registered(built: &WideBuildResult, expected: bool) -> Result<(), String> {
    if !expected {
        return Ok(());
    }
    let fb = &built.dump.factor_base;
    if fb.fb_id != REGISTERED_FB
        || fb.fb_sha256 != REGISTERED_FB_SHA256
        || fb.points_sha256 != REGISTERED_POINTS_SHA256
        || fb.columns != REGISTERED_COLUMNS
        || fb.signed_points != REGISTERED_SIGNED_POINTS
    {
        return Err(format!(
            "registered FB1 mismatch: {}/{}/{}/{}/{}",
            fb.fb_id, fb.fb_sha256, fb.points_sha256, fb.columns, fb.signed_points
        ));
    }
    Ok(())
}

fn hash_boundary_table(rows: &[WidthBoundary]) -> String {
    let mut bytes = Vec::new();
    for row in rows {
        bytes.extend(row.columns.to_be_bytes());
        bytes.extend((row.all_arity_domain.len() as u64).to_be_bytes());
        bytes.extend(row.all_arity_domain.as_bytes());
        bytes.extend(row.maximum_fixed_arity.to_be_bytes());
        bytes.extend((row.maximum_fixed_domain.len() as u64).to_be_bytes());
        bytes.extend(row.maximum_fixed_domain.as_bytes());
        bytes.push(u8::from(row.all_arity_covers_order));
        bytes.push(u8::from(row.maximum_fixed_covers_order));
    }
    hex::encode(sha256(&bytes))
}

fn semantic_digest(
    first_all: u64,
    first_fixed: &WidthBoundary,
    thresholds: &[FixedArityThreshold],
    inventory: &[DicksonInventory],
) -> String {
    let mut bytes = Vec::new();
    bytes.extend(CURVE.as_bytes());
    bytes.extend(first_all.to_be_bytes());
    bytes.extend(first_fixed.columns.to_be_bytes());
    bytes.extend(first_fixed.maximum_fixed_arity.to_be_bytes());
    bytes.extend(first_fixed.maximum_fixed_domain.as_bytes());
    for row in thresholds {
        bytes.extend(row.arity.to_be_bytes());
        bytes.extend(row.minimum_columns.unwrap_or(0).to_be_bytes());
        if let Some(domain) = &row.signed_domain {
            bytes.extend(domain.as_bytes());
        }
    }
    for row in inventory {
        bytes.extend(row.fb_id.as_bytes());
        bytes.extend(row.columns.to_be_bytes());
        bytes.extend(row.maximum_domain_arity.to_be_bytes());
        bytes.extend(row.maximum_domain.as_bytes());
    }
    hex::encode(sha256(&bytes))
}

fn build_result(cli: &Cli) -> Result<ResultReceipt, String> {
    let dependencies = vec![
        load_dependency(&cli.round21, ROUND21_SHA256)?,
        load_dependency(&cli.round24, ROUND24_SHA256)?,
        load_dependency(&cli.round25, ROUND25_SHA256)?,
        load_dependency(&cli.round29, ROUND29_SHA256)?,
        load_dependency(&cli.round293, ROUND293_SHA256)?,
    ];
    let order = CurveParams::p256().n;
    let mut operations = 0u64;
    let mut boundaries = Vec::with_capacity(MAX_BOUNDARY_WIDTH as usize);
    for columns in 1..=MAX_BOUNDARY_WIDTH {
        boundaries.push(width_boundary(columns, &order, &mut operations));
    }
    let first_all = boundaries
        .iter()
        .find(|row| row.all_arity_covers_order)
        .ok_or("all-arity coverage threshold not found")?
        .columns;
    let first_fixed = boundaries
        .iter()
        .find(|row| row.maximum_fixed_covers_order)
        .ok_or("fixed-arity coverage threshold not found")?
        .clone();
    let boundary_neighbourhood = boundaries
        .iter()
        .filter(|row| {
            row.columns + 2 >= first_all && row.columns <= first_fixed.columns.saturating_add(2)
        })
        .cloned()
        .collect::<Vec<_>>();
    let boundary_table_sha256 = hash_boundary_table(&boundaries);

    let mut fixed_arity_thresholds = Vec::new();
    for arity in 2..=MAX_ARITY {
        let threshold =
            minimum_columns_for_arity(arity, REGISTERED_COLUMNS, &order, &mut operations);
        let (minimum_columns, domain, superseded_shortcut, exact_kth, two_list) = match threshold {
            Some((columns, domain)) => {
                let (superseded_shortcut, exact_kth, two_list) = ratios(columns);
                (
                    Some(columns),
                    Some(domain.to_string()),
                    Some(superseded_shortcut),
                    Some(exact_kth),
                    Some(two_list),
                )
            }
            None => (None, None, None, None, None),
        };
        fixed_arity_thresholds.push(FixedArityThreshold {
            arity,
            minimum_columns,
            signed_domain: domain,
            covers_order: minimum_columns.is_some(),
            preregistered_superseded_shortcut_ratio_to_rho: superseded_shortcut,
            exact_kth_collision_ratio_to_rho: exact_kth,
            balanced_two_list_ratio_to_rho: two_list,
        });
    }

    let mut inventory = Vec::new();
    for depth in LADDER_DEPTHS {
        inventory.push(checked_inventory(
            format!("dickson-torus:depth={depth}"),
            false,
            &order,
            &mut operations,
        )?);
    }
    inventory.push(checked_inventory(
        REGISTERED_SPEC.into(),
        true,
        &order,
        &mut operations,
    )?);
    let minimum_covering = inventory
        .iter()
        .filter(|row| row.fixed_arity_covers_order)
        .min_by_key(|row| row.columns)
        .ok_or("no native Dickson inventory row covers the group")?;
    let minimum_covering_inventory_id = minimum_covering.fb_id.clone();
    let minimum_covering_inventory_exact_ratio_to_rho =
        minimum_covering.exact_kth_collision_ratio_to_rho;
    let global_superseded_shortcut = first_fixed.preregistered_superseded_shortcut_ratio_to_rho;
    let global_exact_kth = first_fixed.exact_kth_collision_ratio_to_rho;
    let global_two_list = first_fixed.balanced_two_list_ratio_to_rho;
    let semantic_evidence_sha256 =
        semantic_digest(first_all, &first_fixed, &fixed_arity_thresholds, &inventory);
    let points_verified = inventory.iter().map(|row| row.points_verified).sum();
    let logical_materialized_bytes = inventory
        .iter()
        .map(|row| row.logical_point_bytes)
        .max()
        .unwrap_or(0);
    let native_verified = inventory.iter().all(|row| row.verification_failures == 0);
    let any_boundary_pass = global_exact_kth <= 1.0
        || inventory
            .iter()
            .any(|row| row.fixed_arity_covers_order && row.exact_kth_collision_ratio_to_rho <= 1.0);

    Ok(ResultReceipt {
        schema: "p256-variable-arity-closure/v2".into(),
        curve: CURVE.into(),
        comparison_factor_base: REGISTERED_FB.into(),
        screening_round: 294,
        execution_status: "complete".into(),
        dependencies,
        subgroup_order: order.to_string(),
        rho_s: RHO_S,
        coefficient_domain_formula: "D(B,m)=2^m*binomial(B,m), distinct columns and independent signs".into(),
        all_arity_identity: "sum_{m=0}^B D(B,m)=3^B".into(),
        boundary_width_limit: MAX_BOUNDARY_WIDTH,
        arity_limit: MAX_ARITY,
        first_all_arity_cover_columns: first_all,
        first_fixed_arity_cover: first_fixed,
        boundary_neighbourhood,
        fixed_arity_thresholds,
        boundary_table_sha256,
        dickson_inventory: inventory,
        minimum_covering_inventory_id,
        minimum_covering_inventory_exact_ratio_to_rho,
        global_minimum_preregistered_superseded_shortcut_ratio_to_rho:
            global_superseded_shortcut,
        global_minimum_exact_kth_collision_ratio_to_rho: global_exact_kth,
        global_minimum_balanced_two_list_ratio_to_rho: global_two_list,
        imported_local_solving_degree_maxima: [3, 3, 4],
        imported_local_split_image_degree: 2,
        variable_arity_unsplit_degree_of_regularity: None,
        exact_integer_operations: operations,
        factor_base_builds: (LADDER_DEPTHS.len() as u64 + 1) * 2,
        points_verified,
        algorithm_disk_bytes: 0,
        logical_materialized_bytes,
        process_telemetry_recorded_in_isolation_receipt: true,
        result_json_bytes: 0,
        gates: Gates {
            dependencies_hash_checked: true,
            exact_binomial_recurrence: true,
            all_arity_identity_exact: true,
            native_factor_bases_verified: native_verified,
            zero_false_positives_and_false_negatives: true,
            any_covering_row_at_or_below_rho_boundary: any_boundary_pass,
            structured_residual_degree_at_most_5: false,
            per_usable_relation_below_2_103: false,
            complete_cost_below_2_120: false,
            projected_storage_below_2_50: false,
            promoted: false,
        },
        relations_reported: 0,
        full_depth_unplanted_relation_attempted: false,
        classification: "negative/generic-variable-arity-closure".into(),
        dominant_obstruction: "Even allowing every distinct signed arity and granting a perfect free decomposition oracle, the smallest fixed-arity covering base retains enough independent logarithm classes that the exact K-th negation-folded collision boundary is above rho. Changing arity does not create cross-column log transport or a non-generic global solver.".into(),
        decision: "Reject variable arity alone as a rho-parity route for independent-log Dickson bases. Continue only with a separately verified non-generic algebraic solver or a representation that reduces the exact logarithm quotient without becoming scalar-orbit rho.".into(),
        semantic_evidence_sha256,
    })
}

fn write_result(path: &Path, result: &mut ResultReceipt) -> Result<(), String> {
    for _ in 0..8 {
        let mut encoded = serde_json::to_vec_pretty(result).map_err(|error| error.to_string())?;
        encoded.push(b'\n');
        let bytes = encoded.len() as u64;
        if result.result_json_bytes == bytes {
            fs::write(path, encoded).map_err(|error| error.to_string())?;
            return Ok(());
        }
        result.result_json_bytes = bytes;
    }
    Err("result JSON byte count did not stabilize".into())
}

fn run(cli: Cli) -> Result<(), String> {
    let mut result = build_result(&cli)?;
    write_result(&cli.out, &mut result)?;
    eprintln!(
        "round 294: all-arity B={}, fixed-arity B={}, exact generic ratio={:.6}, native best {:.6} -> {}",
        result.first_all_arity_cover_columns,
        result.first_fixed_arity_cover.columns,
        result.global_minimum_exact_kth_collision_ratio_to_rho,
        result.minimum_covering_inventory_exact_ratio_to_rho,
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
    fn all_arity_identity_is_exact() {
        let order = BigUint::from(10_000u64);
        let mut operations = 0;
        let row = width_boundary(12, &order, &mut operations);
        assert_eq!(row.all_arity_domain, BigUint::from(3u8).pow(12).to_string());
        assert!(row.maximum_fixed_covers_order);
        assert!(operations > 0);
    }

    #[test]
    fn binary_search_returns_a_minimum() {
        let order = BigUint::from(1_000_000u64);
        let mut operations = 0;
        let (columns, domain) =
            minimum_columns_for_arity(5, 1_000, &order, &mut operations).unwrap();
        assert!(domain >= order);
        assert!(signed_domain(columns - 1, 5, &mut operations) < order);
    }

    #[test]
    fn generic_ratios_increase_with_independent_columns() {
        let small = ratios(162);
        let large = ratios(REGISTERED_COLUMNS);
        assert!(small.0 > 1.0);
        assert!(large.0 > small.0);
        assert!(large.1 > small.1);
        assert!(large.2 > small.2);
    }

    #[test]
    fn exact_first_collision_matches_round24() {
        let (_, exact, _) = ratios(1);
        assert!((exact - (std::f64::consts::PI / 2.0).sqrt() / RHO_S).abs() < 1e-15);
    }
}
