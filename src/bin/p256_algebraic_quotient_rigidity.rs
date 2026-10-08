//! Round 304: algebraic self-map and affine factor-base quotient rigidity on P-256.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::CurveParams;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const ROUND25_SHA256: &str = "dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf";
const ROUND297_SHA256: &str = "53949881daf4a48769c9ef8f869fad633c3dc69d11b0d341a5b7416de2e43982";
const ROUND298_SHA256: &str = "8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965";
const ROUND303_SHA256: &str = "6b023f47fbb2c57eb79f870f1c25be48e0b4230d974c3823caf76784a095caa0";
const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const REGISTERED_S17_FACTOR_BASE_ID: &str = "FB1h2f8621cda105";
const COLUMNS: usize = 164;
const TOY_MODULUS: u8 = 7;
const TOY_FUNCTIONS: u64 = 823_543;
const UNQUOTIENTED_RATIO: f64 = 13.920_747_397_073_491;
const REGISTERED_S17_RATIO: f64 = 394.425_280;
const KNOWN_ANCHOR_ORACLE_RATIO: f64 = 0.964_336_477_130_181;
const KNOWN_ANCHOR_DIRECT_RATIO: f64 = 7.714_691_817_041_447;
const UNKNOWN_ANCHOR_ORACLE_RATIO: f64 = 1.446_504_715_695_271_7;
const UNKNOWN_ANCHOR_DIRECT_RATIO: f64 = 11.572_037_725_562_172;

#[derive(Parser)]
#[command(about = "Certify algebraic factor-base quotient rigidity on P-256")]
struct Cli {
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    round297: PathBuf,
    #[arg(long)]
    round298: PathBuf,
    #[arg(long)]
    round303: PathBuf,
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

#[derive(Clone, Serialize)]
struct FunctionExample {
    values: Vec<u8>,
    multiplier: u8,
    offset: u8,
}

#[derive(Clone, Serialize)]
struct FunctionCensus {
    modulus: u8,
    functions_enumerated: u64,
    additive_law_evaluations: u64,
    affine_law_evaluations: u64,
    additive_functions: u64,
    affine_functions: u64,
    bijective_affine_functions: u64,
    non_scalar_additive_functions: u64,
    non_affine_affine_law_functions: u64,
    admitted_table_replays: u64,
    replay_failures: u64,
    first_nonzero_additive: FunctionExample,
    first_nontrivial_affine: FunctionExample,
}

#[derive(Clone, Default, Serialize)]
struct RankOperations {
    matrices: u64,
    rows_materialized: u64,
    entries_materialized: u64,
    pivot_row_scans: u64,
    pivot_swaps: u64,
    modular_inversions: u64,
    modular_multiplications: u64,
    modular_subtractions: u64,
}

impl RankOperations {
    fn add_assign(&mut self, other: &Self) {
        self.matrices += other.matrices;
        self.rows_materialized += other.rows_materialized;
        self.entries_materialized += other.entries_materialized;
        self.pivot_row_scans += other.pivot_row_scans;
        self.pivot_swaps += other.pivot_swaps;
        self.modular_inversions += other.modular_inversions;
        self.modular_multiplications += other.modular_multiplications;
        self.modular_subtractions += other.modular_subtractions;
    }
}

#[derive(Clone)]
struct Equation {
    coefficients: Vec<BigUint>,
    rhs: BigUint,
}

#[derive(Clone, Serialize)]
struct CycleReceipt {
    modulus: String,
    chain_multiplier: String,
    chain_offset: String,
    closing_multiplier: String,
    closing_offset: String,
    composed_multiplier: String,
    composed_offset: String,
    composition_is_identity: bool,
    composition_fixes_witness_anchor: bool,
}

#[derive(Clone, Serialize)]
struct GraphProfile {
    variant: String,
    nodes: usize,
    edges: usize,
    known_affine_edge_labels: bool,
    rank_mod_2_61_minus_1: usize,
    rank_mod_p256_n: usize,
    reversed_rank_mod_2_61_minus_1: usize,
    reversed_rank_mod_p256_n: usize,
    independent_log_quotient_dimension: usize,
    anchor_status: String,
    equation_replays: u64,
    replay_failures: u64,
    rank_replay_discrepancies: u64,
    materialized_bytes_mod_p256_n: u64,
    p61_cycle: Option<CycleReceipt>,
    p256_cycle: Option<CycleReceipt>,
    operations: RankOperations,
}

#[derive(Clone, Serialize)]
struct BoundaryRow {
    variant: String,
    independent_log_classes: usize,
    free_oracle_ratio_to_rho: f64,
    direct_measured_ratio_to_rho: Option<f64>,
    complete_attack_measured: bool,
    status: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_and_schemas_checked: bool,
    exhaustive_function_census_passed: bool,
    zero_non_affine_exact_transports: bool,
    dual_modulus_graph_ranks_and_replays_passed: bool,
    non_scalar_non_affine_p256_log_quotient: bool,
    complete_relation_decomposition_and_recovery: bool,
    structured_residual_degree_at_most_5: bool,
    parity_at_or_below_rho: bool,
    complete_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    projected_storage_below_2_50: bool,
    discarded_probabilistic_branches_counted_as_exhaustive: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    screening_round: u64,
    execution_status: String,
    factor_base: String,
    registered_s17_factor_base: String,
    columns: usize,
    dependencies: Vec<Dependency>,
    rigidity_statement: String,
    induced_log_action: String,
    exhaustive_function_control: FunctionCensus,
    graph_profiles: Vec<GraphProfile>,
    aggregate_rank_operations: RankOperations,
    boundary_table: Vec<BoundaryRow>,
    structured_residual_degree_of_regularity: Option<u64>,
    relations_reported_on_p256: u64,
    full_depth_unplanted_p256_relation_attempted: bool,
    exploration_boundary: String,
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
    methodology_resource: String,
    template_resource: String,
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

fn checked_json(
    path: &Path,
    expected_hash: &str,
    schema: &str,
) -> Result<(Dependency, Value), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != expected_hash {
        return Err(format!(
            "{} SHA-256 is {digest}, expected {expected_hash}",
            path.display()
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("{}: {error}", path.display()))?;
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

fn require_ratio(value: &Value, pointer: &str, expected: f64, label: &str) -> Result<(), String> {
    let actual = value
        .pointer(pointer)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("{label} missing at {pointer}"))?;
    if (actual - expected).abs() > 1e-14 {
        return Err(format!("{label} is {actual}, expected {expected}"));
    }
    Ok(())
}

fn dependency_checks(cli: &Cli) -> Result<Vec<Dependency>, String> {
    let (round25_dep, round25) = checked_json(
        &cli.round25,
        ROUND25_SHA256,
        "p256.scalar_orbit_factor_base_screen/v1",
    )?;
    let (round297_dep, round297) = checked_json(
        &cli.round297,
        ROUND297_SHA256,
        "p256-union-scalar-stabilizer/v1",
    )?;
    let (round298_dep, round298) = checked_json(
        &cli.round298,
        ROUND298_SHA256,
        "p256-parity-escape-screen/v2",
    )?;
    let (round303_dep, round303) = checked_json(
        &cli.round303,
        ROUND303_SHA256,
        "p256-jacobian-row-capacity-screen/v1",
    )?;
    for (label, value, schema) in [
        (
            "Round 25",
            &round25,
            "p256.scalar_orbit_factor_base_screen/v1",
        ),
        ("Round 297", &round297, "p256-union-scalar-stabilizer/v1"),
        ("Round 298", &round298, "p256-parity-escape-screen/v2"),
        (
            "Round 303",
            &round303,
            "p256-jacobian-row-capacity-screen/v1",
        ),
    ] {
        if value.pointer("/schema").and_then(Value::as_str) != Some(schema)
            || value.pointer("/curve").and_then(Value::as_str) != Some(CURVE_SLUG)
        {
            return Err(format!("{label} schema or curve mismatch"));
        }
    }
    let action = round25
        .pointer("/cm_screen/rational_point_action")
        .and_then(Value::as_str)
        .ok_or("Round 25 rational-point action missing")?;
    if !action.contains("known scalar")
        || round25
            .pointer("/factor_base/independent_log_quotient_dimension")
            .and_then(Value::as_u64)
            != Some(1)
    {
        return Err("Round 25 scalar quotient mismatch".into());
    }
    require_ratio(
        &round25,
        "/projections/0/oracle_ratio_to_rho",
        KNOWN_ANCHOR_ORACLE_RATIO,
        "known-anchor oracle ratio",
    )?;
    require_ratio(
        &round25,
        "/projections/0/direct_ratio_to_rho",
        KNOWN_ANCHOR_DIRECT_RATIO,
        "known-anchor direct ratio",
    )?;
    require_ratio(
        &round25,
        "/projections/1/oracle_ratio_to_rho",
        UNKNOWN_ANCHOR_ORACLE_RATIO,
        "unknown-anchor oracle ratio",
    )?;
    require_ratio(
        &round25,
        "/projections/1/direct_ratio_to_rho",
        UNKNOWN_ANCHOR_DIRECT_RATIO,
        "unknown-anchor direct ratio",
    )?;
    if round297
        .pointer("/stabilizer/effective_folded_action_order")
        .and_then(Value::as_u64)
        != Some(1)
        || round297
            .pointer("/stabilizer/independent_log_quotient_dimension")
            .and_then(Value::as_u64)
            != Some(COLUMNS as u64)
    {
        return Err("Round 297 scalar stabilizer mismatch".into());
    }
    require_ratio(
        &round297,
        "/boundary_table/2/ratio_to_rho",
        UNQUOTIENTED_RATIO,
        "Round 297 boundary",
    )?;
    if round298
        .pointer("/p256_boundary/minimum_independent_rows_per_event")
        .and_then(Value::as_u64)
        != Some(COLUMNS as u64)
        || round303
            .pointer("/minimum_capacity/minimum_genus")
            .and_then(Value::as_u64)
            != Some(165)
        || round303
            .pointer("/gates/parity_at_or_below_rho")
            .and_then(Value::as_bool)
            != Some(false)
    {
        return Err("Round 298 or Round 303 boundary mismatch".into());
    }
    Ok(vec![round25_dep, round297_dep, round298_dep, round303_dep])
}

fn decode_function(mut code: u64) -> [u8; TOY_MODULUS as usize] {
    let mut values = [0u8; TOY_MODULUS as usize];
    for value in &mut values {
        *value = (code % u64::from(TOY_MODULUS)) as u8;
        code /= u64::from(TOY_MODULUS);
    }
    values
}

fn mod7_sub(left: u8, right: u8) -> u8 {
    (left + TOY_MODULUS - right) % TOY_MODULUS
}

fn affine_coefficients(values: &[u8; TOY_MODULUS as usize]) -> (u8, u8) {
    (mod7_sub(values[1], values[0]), values[0])
}

fn matches_affine(values: &[u8; TOY_MODULUS as usize], a: u8, b: u8) -> bool {
    values
        .iter()
        .enumerate()
        .all(|(x, value)| *value == ((a as usize * x + b as usize) % TOY_MODULUS as usize) as u8)
}

fn is_bijective(values: &[u8; TOY_MODULUS as usize]) -> bool {
    let mut seen = [false; TOY_MODULUS as usize];
    for value in values {
        seen[*value as usize] = true;
    }
    seen.into_iter().all(|value| value)
}

fn function_census() -> Result<FunctionCensus, String> {
    let mut additive_functions = 0u64;
    let mut affine_functions = 0u64;
    let mut bijective_affine_functions = 0u64;
    let mut non_scalar_additive = 0u64;
    let mut non_affine_affine_law = 0u64;
    let mut replays = 0u64;
    let mut replay_failures = 0u64;
    let mut first_nonzero_additive = None;
    let mut first_nontrivial_affine = None;
    let mut additive_evaluations = 0u64;
    let mut affine_evaluations = 0u64;

    for code in 0..TOY_FUNCTIONS {
        let values = decode_function(code);
        let mut additive = true;
        let mut affine_law = true;
        for x in 0..TOY_MODULUS as usize {
            for y in 0..TOY_MODULUS as usize {
                let sum = (x + y) % TOY_MODULUS as usize;
                additive_evaluations += 1;
                affine_evaluations += 1;
                if values[sum]
                    != ((u16::from(values[x]) + u16::from(values[y])) % u16::from(TOY_MODULUS))
                        as u8
                {
                    additive = false;
                }
                let left = mod7_sub(values[sum], values[0]);
                let right = mod7_sub(
                    ((u16::from(values[x]) + u16::from(values[y])) % u16::from(TOY_MODULUS)) as u8,
                    ((2 * u16::from(values[0])) % u16::from(TOY_MODULUS)) as u8,
                );
                if left != right {
                    affine_law = false;
                }
            }
        }
        let (a, b) = affine_coefficients(&values);
        let affine_replay = matches_affine(&values, a, b);
        if additive {
            additive_functions += 1;
            replays += 1;
            if b != 0 || !affine_replay {
                non_scalar_additive += 1;
                replay_failures += 1;
            }
            if a != 0 && first_nonzero_additive.is_none() {
                first_nonzero_additive = Some(FunctionExample {
                    values: values.to_vec(),
                    multiplier: a,
                    offset: b,
                });
            }
        }
        if affine_law {
            affine_functions += 1;
            replays += 1;
            if !affine_replay {
                non_affine_affine_law += 1;
                replay_failures += 1;
            }
            if is_bijective(&values) {
                bijective_affine_functions += 1;
            }
            if b != 0 && first_nontrivial_affine.is_none() {
                first_nontrivial_affine = Some(FunctionExample {
                    values: values.to_vec(),
                    multiplier: a,
                    offset: b,
                });
            }
        }
    }
    if additive_functions != 7
        || affine_functions != 49
        || bijective_affine_functions != 42
        || non_scalar_additive != 0
        || non_affine_affine_law != 0
        || replay_failures != 0
    {
        return Err("exhaustive affine-function census failed".into());
    }
    Ok(FunctionCensus {
        modulus: TOY_MODULUS,
        functions_enumerated: TOY_FUNCTIONS,
        additive_law_evaluations: additive_evaluations,
        affine_law_evaluations: affine_evaluations,
        additive_functions,
        affine_functions,
        bijective_affine_functions,
        non_scalar_additive_functions: non_scalar_additive,
        non_affine_affine_law_functions: non_affine_affine_law,
        admitted_table_replays: replays,
        replay_failures,
        first_nonzero_additive: first_nonzero_additive.ok_or("no additive example")?,
        first_nontrivial_affine: first_nontrivial_affine.ok_or("no affine example")?,
    })
}

fn mod_sub(left: &BigUint, right: &BigUint, modulus: &BigUint) -> BigUint {
    if left >= right {
        left - right
    } else {
        left + modulus - right
    }
}

fn witness_vector(modulus: &BigUint) -> Vec<BigUint> {
    (0..COLUMNS)
        .map(|index| {
            let index = BigUint::from(index as u64);
            ((&index * &index * &index) + BigUint::from(17u8) * &index + BigUint::from(42u8))
                % modulus
        })
        .collect()
}

fn affine_multiplier(index: usize, modulus: &BigUint) -> BigUint {
    BigUint::from((index % 13 + 1) as u64) % modulus
}

fn chain_system(modulus: &BigUint) -> (Vec<Equation>, BigUint, BigUint, Vec<BigUint>) {
    let witness = witness_vector(modulus);
    let mut equations = Vec::with_capacity(COLUMNS - 1);
    let mut composed_a = BigUint::one();
    let mut composed_b = BigUint::zero();
    for index in 0..(COLUMNS - 1) {
        let a = affine_multiplier(index, modulus);
        let b = mod_sub(
            &witness[index + 1],
            &((&a * &witness[index]) % modulus),
            modulus,
        );
        let mut coefficients = vec![BigUint::zero(); COLUMNS];
        coefficients[index] = mod_sub(&BigUint::zero(), &a, modulus);
        coefficients[index + 1] = BigUint::one();
        equations.push(Equation {
            coefficients,
            rhs: b.clone(),
        });
        composed_b = ((&a * &composed_b) + &b) % modulus;
        composed_a = (&a * &composed_a) % modulus;
    }
    (equations, composed_a, composed_b, witness)
}

fn cycle_equation(
    modulus: &BigUint,
    witness: &[BigUint],
    closing_multiplier: &BigUint,
) -> Equation {
    let mut coefficients = vec![BigUint::zero(); COLUMNS];
    coefficients[0] = BigUint::one();
    coefficients[COLUMNS - 1] = mod_sub(&BigUint::zero(), closing_multiplier, modulus);
    let rhs = mod_sub(
        &witness[0],
        &((closing_multiplier * &witness[COLUMNS - 1]) % modulus),
        modulus,
    );
    Equation { coefficients, rhs }
}

fn dot(coefficients: &[BigUint], values: &[BigUint], modulus: &BigUint) -> BigUint {
    coefficients
        .iter()
        .zip(values)
        .fold(BigUint::zero(), |accumulator, (coefficient, value)| {
            (accumulator + coefficient * value) % modulus
        })
}

fn replay_equations(equations: &[Equation], witness: &[BigUint], modulus: &BigUint) -> u64 {
    equations
        .iter()
        .filter(|equation| dot(&equation.coefficients, witness, modulus) != equation.rhs)
        .count() as u64
}

fn modular_rank(
    equations: &[Equation],
    modulus: &BigUint,
    reverse: bool,
) -> Result<(usize, RankOperations), String> {
    let mut matrix: Vec<Vec<BigUint>> = equations
        .iter()
        .map(|equation| equation.coefficients.clone())
        .collect();
    if reverse {
        matrix.reverse();
    }
    let rows = matrix.len();
    let mut ops = RankOperations {
        matrices: 1,
        rows_materialized: rows as u64,
        entries_materialized: (rows * COLUMNS) as u64,
        ..RankOperations::default()
    };
    let exponent = modulus - BigUint::from(2u8);
    let mut rank = 0usize;
    for column in 0..COLUMNS {
        let mut pivot = None;
        for row in rank..rows {
            ops.pivot_row_scans += 1;
            if !matrix[row][column].is_zero() {
                pivot = Some(row);
                break;
            }
        }
        let Some(pivot) = pivot else {
            continue;
        };
        if pivot != rank {
            matrix.swap(pivot, rank);
            ops.pivot_swaps += 1;
        }
        let inverse = matrix[rank][column].modpow(&exponent, modulus);
        ops.modular_inversions += 1;
        if (&matrix[rank][column] * &inverse) % modulus != BigUint::one() {
            return Err("rank pivot inversion failed".into());
        }
        for entry in &mut matrix[rank][column..] {
            *entry = (&*entry * &inverse) % modulus;
            ops.modular_multiplications += 1;
        }
        let pivot_row = matrix[rank].clone();
        for row in (rank + 1)..rows {
            let factor = matrix[row][column].clone();
            if factor.is_zero() {
                continue;
            }
            for entry_column in column..COLUMNS {
                let product = (&factor * &pivot_row[entry_column]) % modulus;
                matrix[row][entry_column] = mod_sub(&matrix[row][entry_column], &product, modulus);
                ops.modular_multiplications += 1;
                ops.modular_subtractions += 1;
            }
        }
        rank += 1;
        if rank == rows {
            break;
        }
    }
    Ok((rank, ops))
}

fn cycle_receipt(
    modulus: &BigUint,
    chain_a: &BigUint,
    chain_b: &BigUint,
    closing: &Equation,
    closing_a: &BigUint,
    witness: &[BigUint],
) -> CycleReceipt {
    let composed_a = (closing_a * chain_a) % modulus;
    let composed_b = ((closing_a * chain_b) + &closing.rhs) % modulus;
    let fixed = ((&composed_a * &witness[0]) + &composed_b) % modulus == witness[0];
    CycleReceipt {
        modulus: modulus.to_string(),
        chain_multiplier: chain_a.to_string(),
        chain_offset: chain_b.to_string(),
        closing_multiplier: closing_a.to_string(),
        closing_offset: closing.rhs.to_string(),
        composed_multiplier: composed_a.to_string(),
        composed_offset: composed_b.to_string(),
        composition_is_identity: composed_a.is_one() && composed_b.is_zero(),
        composition_fixes_witness_anchor: fixed,
    }
}

fn graph_equations(
    variant: &str,
    modulus: &BigUint,
) -> Result<(Vec<Equation>, Option<CycleReceipt>, Vec<BigUint>), String> {
    let (mut equations, chain_a, chain_b, witness) = chain_system(modulus);
    match variant {
        "independent" => Ok((Vec::new(), None, witness)),
        "connected-chain" => Ok((equations, None, witness)),
        "dependent-cycle" | "anchor-solving-cycle" => {
            let inverse_a = chain_a.modpow(&(modulus - BigUint::from(2u8)), modulus);
            if (&chain_a * &inverse_a) % modulus != BigUint::one() {
                return Err("chain multiplier inversion failed".into());
            }
            let closing_a = if variant == "dependent-cycle" {
                inverse_a
            } else {
                (BigUint::from(2u8) * inverse_a) % modulus
            };
            let closing = cycle_equation(modulus, &witness, &closing_a);
            let receipt =
                cycle_receipt(modulus, &chain_a, &chain_b, &closing, &closing_a, &witness);
            equations.push(closing);
            Ok((equations, Some(receipt), witness))
        }
        _ => Err(format!("unknown graph variant {variant}")),
    }
}

fn graph_profile(variant: &str, p256_n: &BigUint) -> Result<GraphProfile, String> {
    let p61 = (BigUint::one() << 61usize) - BigUint::one();
    let (eq61, cycle61, witness61) = graph_equations(variant, &p61)?;
    let (eqn, cyclen, witnessn) = graph_equations(variant, p256_n)?;
    let failures =
        replay_equations(&eq61, &witness61, &p61) + replay_equations(&eqn, &witnessn, p256_n);
    let (rank61, ops61) = modular_rank(&eq61, &p61, false)?;
    let (reverse61, reverse_ops61) = modular_rank(&eq61, &p61, true)?;
    let (rankn, opsn) = modular_rank(&eqn, p256_n, false)?;
    let (reversen, reverse_opsn) = modular_rank(&eqn, p256_n, true)?;
    let expected = match variant {
        "independent" => 0,
        "connected-chain" | "dependent-cycle" => COLUMNS - 1,
        "anchor-solving-cycle" => COLUMNS,
        _ => return Err("unknown graph variant".into()),
    };
    let discrepancies = [rank61, reverse61, rankn, reversen]
        .into_iter()
        .filter(|rank| *rank != expected)
        .count() as u64;
    if failures != 0 || discrepancies != 0 {
        return Err(format!("{variant} graph replay or rank failed"));
    }
    if variant == "dependent-cycle"
        && (!cycle61
            .as_ref()
            .is_some_and(|receipt| receipt.composition_is_identity)
            || !cyclen
                .as_ref()
                .is_some_and(|receipt| receipt.composition_is_identity))
    {
        return Err("dependent cycle does not compose to identity".into());
    }
    if variant == "anchor-solving-cycle"
        && (cycle61
            .as_ref()
            .is_some_and(|receipt| receipt.composition_is_identity)
            || cyclen
                .as_ref()
                .is_some_and(|receipt| receipt.composition_is_identity)
            || !cycle61
                .as_ref()
                .is_some_and(|receipt| receipt.composition_fixes_witness_anchor)
            || !cyclen
                .as_ref()
                .is_some_and(|receipt| receipt.composition_fixes_witness_anchor))
    {
        return Err("anchor-solving cycle certificate failed".into());
    }
    let mut operations = RankOperations::default();
    for counts in [&ops61, &reverse_ops61, &opsn, &reverse_opsn] {
        operations.add_assign(counts);
    }
    Ok(GraphProfile {
        variant: variant.into(),
        nodes: COLUMNS,
        edges: eqn.len(),
        known_affine_edge_labels: variant != "independent",
        rank_mod_2_61_minus_1: rank61,
        rank_mod_p256_n: rankn,
        reversed_rank_mod_2_61_minus_1: reverse61,
        reversed_rank_mod_p256_n: reversen,
        independent_log_quotient_dimension: COLUMNS - expected,
        anchor_status: match variant {
            "independent" => "164 independent logs".into(),
            "connected-chain" | "dependent-cycle" => "one unknown anchor".into(),
            "anchor-solving-cycle" => "anchor solved by an independent affine label".into(),
            _ => unreachable!(),
        },
        equation_replays: (eq61.len() + eqn.len()) as u64,
        replay_failures: failures,
        rank_replay_discrepancies: discrepancies,
        materialized_bytes_mod_p256_n: (eqn.len() * COLUMNS * 32) as u64,
        p61_cycle: cycle61,
        p256_cycle: cyclen,
        operations,
    })
}

fn boundary_table() -> Vec<BoundaryRow> {
    vec![
        BoundaryRow {
            variant: "Pollard rho".into(),
            independent_log_classes: 0,
            free_oracle_ratio_to_rho: 1.0,
            direct_measured_ratio_to_rho: Some(1.0),
            complete_attack_measured: true,
            status: "reference".into(),
        },
        BoundaryRow {
            variant: format!("unquotiented geometry-defined {FACTOR_BASE_ID}"),
            independent_log_classes: COLUMNS,
            free_oracle_ratio_to_rho: UNQUOTIENTED_RATIO,
            direct_measured_ratio_to_rho: None,
            complete_attack_measured: false,
            status: "no replayable log transport".into(),
        },
        BoundaryRow {
            variant: "connected affine orbit, unknown anchor".into(),
            independent_log_classes: 1,
            free_oracle_ratio_to_rho: UNKNOWN_ANCHOR_ORACLE_RATIO,
            direct_measured_ratio_to_rho: Some(UNKNOWN_ANCHOR_DIRECT_RATIO),
            complete_attack_measured: false,
            status: "free oracle already above rho; direct Round-25 control higher".into(),
        },
        BoundaryRow {
            variant: "affine orbit, known or cycle-solved anchor".into(),
            independent_log_classes: 0,
            free_oracle_ratio_to_rho: KNOWN_ANCHOR_ORACLE_RATIO,
            direct_measured_ratio_to_rho: Some(KNOWN_ANCHOR_DIRECT_RATIO),
            complete_attack_measured: false,
            status: "scalar-labelled generic control; Round-25 direct implementation above rho"
                .into(),
        },
        BoundaryRow {
            variant: format!("registered 17-term {REGISTERED_S17_FACTOR_BASE_ID}"),
            independent_log_classes: COLUMNS,
            free_oracle_ratio_to_rho: REGISTERED_S17_RATIO,
            direct_measured_ratio_to_rho: None,
            complete_attack_measured: false,
            status: "unchanged free-perfect-oracle comparison".into(),
        },
    ]
}

fn obligation(name: &str, status: &str, evidence: &str, scope: &str) -> Obligation {
    Obligation {
        name: name.into(),
        status: status.into(),
        evidence: evidence.into(),
        scope: scope.into(),
    }
}

fn build_assessment(census: &FunctionCensus, graphs: &[GraphProfile]) -> TransferAssessment {
    let mut assessment = TransferAssessment {
        schema: "p256-algebraic-quotient-rigidity-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 304,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            format!("{CURVE_SLUG}/E --[rational self-map f]--> E"),
            "f(P)=phi(P)+T, with phi in End(E) and T=f(O)".into(),
            "E(F_p)[n] --[Round-25 action]--> log(f(P))=a*log(P)+b".into(),
            "factor-base affine transport graph --[connected components]--> one anchor variable per component".into(),
            "independent affine cycle --[known label]--> solved anchor / scalar-labelled generic control".into(),
        ],
        obligations: vec![
            obligation("algebraic self-map classification", "supported", "Elliptic-curve morphism rigidity gives translation plus endomorphism; Round 25 gives scalar rational-subgroup action.", "Rational algebraic self-maps of P-256."),
            obligation("finite affine-law control", "supported", &format!("All {} functions Z/7Z->Z/7Z were checked; exactly 7 additive and 49 affine functions occurred.", census.functions_enumerated), "Finite exact control for the induced log algebra."),
            obligation("connected quotient dimension", "supported", "The 164-node connected chain and dependent cycle have dual-modulus rank 163; one anchor remains.", "Known affine transport labels."),
            obligation("anchor-solving cycle", "supported as algebra", "An independent labelled cycle raises dual-modulus rank to 164 and replays the deterministic solution.", "A scalar-labelled control, not evidence that geometry reveals the labels."),
            obligation("non-affine global P-256 quotient", "unknown", "No replayable non-affine log transport is constructed; coordinate permutations receive no log quotient credit.", "Future factor-base mechanisms."),
            obligation("complete below-rho recovery", "unknown", "The known-anchor Round-25 implementation is above rho and no new relation/decomposition pipeline exists.", "End-to-end P-256 attack."),
        ],
        controls: vec![
            format!("The exhaustive function census performed {} additive and {} affine law evaluations with zero replay failures.", census.additive_law_evaluations, census.affine_law_evaluations),
            format!("{} graph variants were ranked over two fields in forward and reversed row order.", graphs.len()),
            "Round 25, 297, 298, and 303 were imported only after exact byte-hash, schema, curve, and conclusion checks.".into(),
        ],
        cost_accounting: vec![
            "Measured: exhaustive function tables and law evaluations; graph rows, matrix entries, pivot scans, swaps, inversions, multiplications, subtractions, equation replays, process telemetry, and artifact sizes.".into(),
            "Imported without extrapolation: Round-25 known/unknown-anchor oracle and direct ratios, and the Round-297 unquotiented boundary.".into(),
            "Unset: any non-affine factor base, structured relation solver, yield, sparse linear algebra, target recovery, and complete new attack cost.".into(),
        ],
        exploration_boundary: "Exact for algebraic P-256 self-map transports and finite affine-log quotient graphs; not a universal impossibility theorem for non-homomorphic relation mechanisms.".into(),
        weakest_open_obligation: "Produce a replayable P-256 factor-base quotient that is not affine scalar transport, or a non-homomorphic relation mechanism whose rank and complete recovery cost beat rho.".into(),
        narrowest_supported_finding: "Algebraic P-256 self-map symmetry induces only affine scalar log transport: a connected base retains one unknown anchor unless known scalar labels solve it, and neither registered anchor case beats rho end to end.".into(),
        semantic_evidence_sha256: String::new(),
        assessment_json_bytes: 0,
    };
    let semantic = json!({
        "curve": assessment.curve,
        "typed_correspondence_graph": assessment.typed_correspondence_graph,
        "obligations": assessment.obligations,
        "controls": assessment.controls,
        "cost_accounting": assessment.cost_accounting,
        "exploration_boundary": assessment.exploration_boundary,
        "weakest_open_obligation": assessment.weakest_open_obligation,
        "narrowest_supported_finding": assessment.narrowest_supported_finding,
    });
    assessment.semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).expect("assessment semantic JSON"));
    assessment
}

fn build_result(cli: &Cli) -> Result<(ResultReceipt, TransferAssessment), String> {
    let dependencies = dependency_checks(cli)?;
    let curve = CurveParams::p256();
    if curve.h != 1 || curve.n.bits() != 256 {
        return Err("P-256 subgroup prerequisites failed".into());
    }
    let census = function_census()?;
    let mut graphs = Vec::new();
    let mut aggregate = RankOperations::default();
    for variant in [
        "independent",
        "connected-chain",
        "dependent-cycle",
        "anchor-solving-cycle",
    ] {
        let graph = graph_profile(variant, &curve.n)?;
        aggregate.add_assign(&graph.operations);
        graphs.push(graph);
    }
    let gates = Gates {
        dependency_hashes_and_schemas_checked: true,
        exhaustive_function_census_passed: true,
        zero_non_affine_exact_transports: true,
        dual_modulus_graph_ranks_and_replays_passed: true,
        non_scalar_non_affine_p256_log_quotient: false,
        complete_relation_decomposition_and_recovery: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        complete_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        projected_storage_below_2_50: false,
        discarded_probabilistic_branches_counted_as_exhaustive: false,
        promoted: false,
    };
    let assessment = build_assessment(&census, &graphs);
    let boundaries = boundary_table();
    let classification =
        "algebraic-factor-base-quotient/affine-scalar-rigidity/anchor-conservation/parity-blocked";
    let obstruction = "Every algebraic P-256 self-map induces log transport x->a*x+b. A connected affine transport graph reduces 164 logs only to one unknown anchor, whose free-oracle boundary is already 1.446505 times rho. Solving the anchor requires an independent known affine label and yields a scalar-labelled generic control; the registered direct instance costs 7.714692 times rho. No replayable non-affine quotient or complete below-rho pipeline is constructed.";
    let decision = "Reject algebraic coordinate symmetry as an unpriced logarithm quotient. Credit only replayable affine scalar edges, retain one variable per unknown-anchor component, and classify solved-anchor orbits as scalar-labelled generic controls. Keep genuinely non-homomorphic relation mechanisms open; do not attempt an unplanted full-depth relation.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": &dependencies,
        "rigidity_statement": "every algebraic self-map is translation plus endomorphism",
        "exhaustive_function_control": &census,
        "graph_profiles": &graphs,
        "boundary_table": &boundaries,
        "assessment": assessment.semantic_evidence_sha256,
        "gates": &gates,
        "classification": classification,
        "dominant_obstruction": obstruction,
        "decision": decision,
    });
    let semantic_hash =
        sha256_hex(&serde_json::to_vec(&semantic).map_err(|error| error.to_string())?);
    Ok((
        ResultReceipt {
            schema: "p256-algebraic-quotient-rigidity-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 304,
            execution_status: "complete".into(),
            factor_base: FACTOR_BASE_ID.into(),
            registered_s17_factor_base: REGISTERED_S17_FACTOR_BASE_ID.into(),
            columns: COLUMNS,
            dependencies,
            rigidity_statement: "Every rational algebraic self-map f:E->E has f(P)=phi(P)+T; on E(F_p)[n], Round 25 gives phi(P)=[a]P, so the induced log action is affine x->a*x+b.".into(),
            induced_log_action: "log_G(f(P))=a*log_G(P)+b mod n; unknown b is an anchor variable, known b is a scalar label.".into(),
            exhaustive_function_control: census,
            graph_profiles: graphs,
            aggregate_rank_operations: aggregate,
            boundary_table: boundaries,
            structured_residual_degree_of_regularity: None,
            relations_reported_on_p256: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            exploration_boundary: "Exact for algebraic self-map transports and affine factor-base quotient graphs; not universal over non-homomorphic relation mechanisms.".into(),
            transfer_assessment_semantic_sha256: assessment.semantic_evidence_sha256.clone(),
            gates,
            classification: classification.into(),
            dominant_obstruction: obstruction.into(),
            decision: decision.into(),
            semantic_evidence_sha256: semantic_hash,
            result_json_bytes: 0,
        },
        assessment,
    ))
}

fn write_result(path: &Path, result: &mut ResultReceipt) -> Result<(), String> {
    loop {
        let text = serde_json::to_string_pretty(result).map_err(|error| error.to_string())? + "\n";
        let bytes = text.len() as u64;
        if bytes == result.result_json_bytes {
            fs::write(path, text).map_err(|error| format!("{}: {error}", path.display()))?;
            return Ok(());
        }
        result.result_json_bytes = bytes;
    }
}

fn write_assessment(path: &Path, assessment: &mut TransferAssessment) -> Result<(), String> {
    loop {
        let text =
            serde_json::to_string_pretty(assessment).map_err(|error| error.to_string())? + "\n";
        let bytes = text.len() as u64;
        if bytes == assessment.assessment_json_bytes {
            fs::write(path, text).map_err(|error| format!("{}: {error}", path.display()))?;
            return Ok(());
        }
        assessment.assessment_json_bytes = bytes;
    }
}

fn run(cli: Cli) -> Result<(), String> {
    let (mut result, mut assessment) = build_result(&cli)?;
    write_assessment(&cli.assessment, &mut assessment)?;
    write_result(&cli.out, &mut result)?;
    eprintln!(
        "round 304: functions={}, affine={}, graphs={}, promoted={}",
        result.exhaustive_function_control.functions_enumerated,
        result.exhaustive_function_control.affine_functions,
        result.graph_profiles.len(),
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_algebraic_quotient_rigidity: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exhaustive_affine_function_counts() {
        let census = function_census().expect("function census");
        assert_eq!(census.functions_enumerated, TOY_FUNCTIONS);
        assert_eq!(census.additive_functions, 7);
        assert_eq!(census.affine_functions, 49);
        assert_eq!(census.bijective_affine_functions, 42);
        assert_eq!(census.replay_failures, 0);
    }

    #[test]
    fn graph_ranks_have_expected_anchor_dimensions() {
        let modulus = (BigUint::one() << 61usize) - BigUint::one();
        for (variant, expected_rank) in [
            ("independent", 0usize),
            ("connected-chain", 163),
            ("dependent-cycle", 163),
            ("anchor-solving-cycle", 164),
        ] {
            let (equations, _, witness) = graph_equations(variant, &modulus).expect("graph");
            assert_eq!(replay_equations(&equations, &witness, &modulus), 0);
            assert_eq!(
                modular_rank(&equations, &modulus, false).expect("rank").0,
                expected_rank
            );
        }
    }

    #[test]
    fn cycle_compositions_distinguish_anchor() {
        let modulus = (BigUint::one() << 61usize) - BigUint::one();
        let (_, dependent, _) = graph_equations("dependent-cycle", &modulus).expect("dependent");
        let (_, solved, _) = graph_equations("anchor-solving-cycle", &modulus).expect("solved");
        assert!(dependent.expect("cycle").composition_is_identity);
        let solved = solved.expect("cycle");
        assert!(!solved.composition_is_identity);
        assert!(solved.composition_fixes_witness_anchor);
    }
}
