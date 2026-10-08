//! Round 303: exact homomorphic cover/Jacobian row-capacity certificate for P-256.

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
const ROUND298_SHA256: &str = "8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965";
const ROUND300_SHA256: &str = "7d65f3e015c6d6c6d24b64ea2e59a6e9e8653c88a1900835838b1bc0a1ec36e5";
const ROUND302_SHA256: &str = "a9a12d823cbca6df965d28f8372217b08e65c04326e13e408bba9dac9f1819c1";
const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const REGISTERED_S17_FACTOR_BASE_ID: &str = "FB1h2f8621cda105";
const COLUMNS: usize = 164;
const ARITY: usize = 109;
const RHO_S: f64 = 1.3;
const REGISTERED_S17_RATIO_TO_RHO: f64 = 394.425_280;
const GENUS_PROFILES: [usize; 11] = [2, 3, 5, 9, 17, 33, 65, 83, 129, 164, 165];

#[derive(Parser)]
#[command(about = "Certify the P-256 homomorphic Jacobian row-capacity bound")]
struct Cli {
    #[arg(long)]
    round25: PathBuf,
    #[arg(long)]
    round298: PathBuf,
    #[arg(long)]
    round300: PathBuf,
    #[arg(long)]
    round302: PathBuf,
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

#[derive(Clone, Serialize)]
struct CapacityProfile {
    genus: usize,
    prym_dimension: usize,
    maximum_p256_prym_multiplicity: usize,
    maximum_same_target_homogeneous_rank: usize,
    collision_events_required: usize,
    collision_s: f64,
    optimistic_ratio_to_rho: f64,
    parity_capacity_possible: bool,
    character_bits: usize,
    complete_character_sheets: usize,
    character_columns: usize,
    rank_mod_2_61_minus_1: usize,
    rank_mod_p256_n: usize,
    replay_rank_mod_2_61_minus_1: usize,
    replay_rank_mod_p256_n: usize,
    rank_replay_discrepancies: u64,
    materialized_bytes_mod_2_61_minus_1: u64,
    materialized_bytes_mod_p256_n: u64,
    operations: RankOperations,
}

#[derive(Serialize)]
struct MinimumCapacity {
    required_homogeneous_rank: usize,
    minimum_genus: usize,
    minimum_prym_dimension: usize,
    minimum_p256_prym_multiplicity: usize,
    minimum_distinct_normalized_fibre_preimages: usize,
    minimum_geometric_cover_degree: usize,
    minimum_total_ramification_degree_at_minimum_genus: usize,
    minimum_raw_folded_preimages_under_round298_sign_model: usize,
    explicit_cover_constructed: bool,
    capacity_is_attack: bool,
}

#[derive(Serialize)]
struct NegativeControl {
    construction: String,
    cover_genus: usize,
    prym_dimension: usize,
    rational_p256_order_n_prym_factors: usize,
    maximum_same_target_sheet_rank_gain: usize,
    imported_round300_maximum_quotient_rank_gain: usize,
    imported_round302_order_n_transport_exists: bool,
    status: String,
}

#[derive(Serialize)]
struct BoundaryRow {
    variant: String,
    maximum_homogeneous_rank_per_event: Option<usize>,
    collision_events_required: Option<usize>,
    optimistic_ratio_to_rho: f64,
    complete_attack_measured: bool,
    status: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_and_schemas_checked: bool,
    exact_capacity_bound_derived: bool,
    dual_modulus_character_ranks_and_replays_pass: bool,
    explicit_genus_165_degree_165_p256_cover: bool,
    p256_prym_multiplicity_at_least_164: bool,
    rank_164_same_target_fibre_replayed: bool,
    explicit_forward_inverse_kernel_and_recovery: bool,
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
    arity: usize,
    dependencies: Vec<Dependency>,
    conditional_capacity_statement: String,
    scalar_action_import: String,
    capacity_profiles: Vec<CapacityProfile>,
    minimum_capacity: MinimumCapacity,
    registered_split_cover_negative_control: NegativeControl,
    aggregate_rank_operations: RankOperations,
    peak_projected_character_matrix_bytes: u64,
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

fn dependency_checks(cli: &Cli) -> Result<Vec<Dependency>, String> {
    let (round25_dep, round25) = checked_json(
        &cli.round25,
        ROUND25_SHA256,
        "p256.scalar_orbit_factor_base_screen/v1",
    )?;
    let (round298_dep, round298) = checked_json(
        &cli.round298,
        ROUND298_SHA256,
        "p256-parity-escape-screen/v2",
    )?;
    let (round300_dep, round300) =
        checked_json(&cli.round300, ROUND300_SHA256, "p256-cover-fiber-screen/v1")?;
    let (round302_dep, round302) = checked_json(
        &cli.round302,
        ROUND302_SHA256,
        "p256-complementary-quotient-screen/v1",
    )?;

    for (label, value, schema) in [
        (
            "Round 25",
            &round25,
            "p256.scalar_orbit_factor_base_screen/v1",
        ),
        ("Round 298", &round298, "p256-parity-escape-screen/v2"),
        ("Round 300", &round300, "p256-cover-fiber-screen/v1"),
        (
            "Round 302",
            &round302,
            "p256-complementary-quotient-screen/v1",
        ),
    ] {
        if value.pointer("/schema").and_then(Value::as_str) != Some(schema)
            || value.pointer("/curve").and_then(Value::as_str) != Some(CURVE_SLUG)
        {
            return Err(format!("{label} schema or curve mismatch"));
        }
    }

    let scalar_action = round25
        .pointer("/cm_screen/rational_point_action")
        .and_then(Value::as_str)
        .ok_or("Round 25 scalar action missing")?;
    if !scalar_action.contains("known scalar")
        || round25
            .pointer("/cm_screen/low_degree_noninteger_transport_exists")
            .and_then(Value::as_bool)
            != Some(false)
    {
        return Err("Round 25 scalar-action conclusion mismatch".into());
    }
    let boundary = round298
        .pointer("/p256_boundary")
        .ok_or("Round 298 P-256 boundary missing")?;
    let expected_n = CurveParams::p256().n.to_string();
    if boundary.pointer("/columns").and_then(Value::as_u64) != Some(COLUMNS as u64)
        || boundary.pointer("/arity").and_then(Value::as_u64) != Some(ARITY as u64)
        || boundary.pointer("/subgroup_order").and_then(Value::as_str) != Some(expected_n.as_str())
        || boundary
            .pointer("/parity_collision_events")
            .and_then(Value::as_u64)
            != Some(1)
        || boundary
            .pointer("/minimum_independent_rows_per_event")
            .and_then(Value::as_u64)
            != Some(COLUMNS as u64)
        || boundary
            .pointer("/minimum_distinct_bucket_occupancy_for_collision_rows")
            .and_then(Value::as_u64)
            != Some((COLUMNS + 1) as u64)
    {
        return Err("Round 298 parity boundary mismatch".into());
    }
    let one = boundary
        .pointer("/parity_event_ratio_to_rho")
        .and_then(Value::as_f64)
        .ok_or("Round 298 one-event ratio missing")?;
    let two = boundary
        .pointer("/two_event_ratio_to_rho")
        .and_then(Value::as_f64)
        .ok_or("Round 298 two-event ratio missing")?;
    if (one - collision_ratio(1)).abs() > 1e-14 || (two - collision_ratio(2)).abs() > 1e-14 {
        return Err("Round 298 collision recurrence mismatch".into());
    }
    if round300
        .pointer("/maximum_observed_quotient_rank_gain")
        .and_then(Value::as_u64)
        != Some(0)
        || round302
            .pointer("/subgroup_certificate/rational_order_n_transport_exists")
            .and_then(Value::as_bool)
            != Some(false)
        || round302
            .pointer("/distinct_p256_order_n_rows_per_cover_event")
            .and_then(Value::as_u64)
            != Some(1)
    {
        return Err("Round 300 or Round 302 negative control mismatch".into());
    }
    Ok(vec![round25_dep, round298_dep, round300_dep, round302_dep])
}

fn collision_s(events: usize) -> f64 {
    assert!(events >= 1);
    let mut s = (std::f64::consts::PI / 2.0).sqrt();
    for q in 1..events {
        s *= (q as f64 + 0.5) / q as f64;
    }
    s
}

fn collision_ratio(events: usize) -> f64 {
    collision_s(events) / RHO_S
}

fn character_bits(columns: usize) -> usize {
    let required_sheets = columns + 1;
    let mut bits = 0usize;
    let mut sheets = 1usize;
    while sheets < required_sheets {
        sheets <<= 1;
        bits += 1;
    }
    bits
}

fn character_matrix(columns: usize, modulus: &BigUint, reverse: bool) -> Vec<Vec<BigUint>> {
    let bits = character_bits(columns);
    let sheets = 1usize << bits;
    let minus_two = modulus - BigUint::from(2u8);
    let mut order: Vec<usize> = (1..sheets).collect();
    if reverse {
        order.reverse();
    }
    order
        .into_iter()
        .map(|sheet| {
            (1..=columns)
                .map(|mask| {
                    if (mask & sheet).count_ones() & 1 == 1 {
                        minus_two.clone()
                    } else {
                        BigUint::zero()
                    }
                })
                .collect()
        })
        .collect()
}

fn modular_rank(
    columns: usize,
    modulus: &BigUint,
    reverse: bool,
) -> Result<(usize, RankOperations), String> {
    let mut matrix = character_matrix(columns, modulus, reverse);
    let rows = matrix.len();
    let mut ops = RankOperations {
        matrices: 1,
        rows_materialized: rows as u64,
        entries_materialized: (rows * columns) as u64,
        ..RankOperations::default()
    };
    let mut rank = 0usize;
    let exponent = modulus - BigUint::from(2u8);
    for column in 0..columns {
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
            return Err("pivot inversion replay failed".into());
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
            for entry_column in column..columns {
                let product = (&factor * &pivot_row[entry_column]) % modulus;
                ops.modular_multiplications += 1;
                matrix[row][entry_column] = if matrix[row][entry_column] >= product {
                    &matrix[row][entry_column] - &product
                } else {
                    &matrix[row][entry_column] + modulus - &product
                };
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

fn capacity_profile(genus: usize, p256_n: &BigUint) -> Result<CapacityProfile, String> {
    let cap = genus.checked_sub(1).ok_or("genus must be positive")?;
    let events = COLUMNS.div_ceil(cap);
    let bits = character_bits(cap);
    let sheets = 1usize << bits;
    let p61 = (BigUint::one() << 61usize) - BigUint::one();
    let (rank61, ops61) = modular_rank(cap, &p61, false)?;
    let (rank_n, ops_n) = modular_rank(cap, p256_n, false)?;
    let (replay61, replay_ops61) = modular_rank(cap, &p61, true)?;
    let (replay_n, replay_ops_n) = modular_rank(cap, p256_n, true)?;
    let discrepancies = [rank61, rank_n, replay61, replay_n]
        .into_iter()
        .filter(|rank| *rank != cap)
        .count() as u64;
    if discrepancies != 0 {
        return Err(format!("genus {genus} character-rank replay failed"));
    }
    let mut operations = RankOperations::default();
    for counts in [&ops61, &ops_n, &replay_ops61, &replay_ops_n] {
        operations.add_assign(counts);
    }
    let entries = (sheets - 1) * cap;
    Ok(CapacityProfile {
        genus,
        prym_dimension: cap,
        maximum_p256_prym_multiplicity: cap,
        maximum_same_target_homogeneous_rank: cap,
        collision_events_required: events,
        collision_s: collision_s(events),
        optimistic_ratio_to_rho: collision_ratio(events),
        parity_capacity_possible: events == 1,
        character_bits: bits,
        complete_character_sheets: sheets,
        character_columns: cap,
        rank_mod_2_61_minus_1: rank61,
        rank_mod_p256_n: rank_n,
        replay_rank_mod_2_61_minus_1: replay61,
        replay_rank_mod_p256_n: replay_n,
        rank_replay_discrepancies: discrepancies,
        materialized_bytes_mod_2_61_minus_1: (entries * 8) as u64,
        materialized_bytes_mod_p256_n: (entries * 32) as u64,
        operations,
    })
}

fn boundary_table(profiles: &[CapacityProfile]) -> Result<Vec<BoundaryRow>, String> {
    let row = |genus: usize| -> Result<&CapacityProfile, String> {
        profiles
            .iter()
            .find(|profile| profile.genus == genus)
            .ok_or_else(|| format!("missing genus {genus} profile"))
    };
    let profile_row = |genus: usize, status: &str| -> Result<BoundaryRow, String> {
        let profile = row(genus)?;
        Ok(BoundaryRow {
            variant: format!("homomorphic cover capacity at genus {genus}"),
            maximum_homogeneous_rank_per_event: Some(profile.maximum_same_target_homogeneous_rank),
            collision_events_required: Some(profile.collision_events_required),
            optimistic_ratio_to_rho: profile.optimistic_ratio_to_rho,
            complete_attack_measured: false,
            status: status.into(),
        })
    };
    Ok(vec![
        BoundaryRow {
            variant: "Pollard rho".into(),
            maximum_homogeneous_rank_per_event: None,
            collision_events_required: None,
            optimistic_ratio_to_rho: 1.0,
            complete_attack_measured: true,
            status: "reference".into(),
        },
        profile_row(
            2,
            "degree-two/self-gluing ceiling; one Prym direction does not reduce events",
        )?,
        profile_row(83, "82-row event cap still needs two events")?,
        profile_row(164, "163-row event cap still needs two events")?,
        profile_row(
            165,
            "first capacity below rho; abstract character control only, no cover or attack",
        )?,
        BoundaryRow {
            variant: format!("registered 17-term {REGISTERED_S17_FACTOR_BASE_ID}"),
            maximum_homogeneous_rank_per_event: Some(1),
            collision_events_required: Some(COLUMNS),
            optimistic_ratio_to_rho: REGISTERED_S17_RATIO_TO_RHO,
            complete_attack_measured: false,
            status: "unchanged free-perfect-oracle comparison".into(),
        },
    ])
}

fn obligation(name: &str, status: &str, evidence: &str, scope: &str) -> Obligation {
    Obligation {
        name: name.into(),
        status: status.into(),
        evidence: evidence.into(),
        scope: scope.into(),
    }
}

fn build_assessment(profiles: &[CapacityProfile]) -> TransferAssessment {
    let full = profiles
        .iter()
        .find(|profile| profile.genus == 165)
        .expect("genus-165 profile");
    let mut assessment = TransferAssessment {
        schema: "p256-jacobian-row-capacity-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 303,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            format!("C/F_p --[separable pi]--> {CURVE_SLUG}/E"),
            "J(C) --[pi_*]--> E; Prym(pi)=ker(pi_*)^0 with dim g-1".into(),
            "Prym(pi) --[E-isotypic quotient]--> E^m, with m<=g-1".into(),
            "E^m(F_p)[n] --[Round-25 scalar quotient]--> at most m independent log coordinates".into(),
            "one pi-fibre --[Abel-Jacobi differences]--> Prym(pi)(F_p)[n]".into(),
        ],
        obligations: vec![
            obligation("Prym dimension ceiling", "supported", "For a separable cover C->E, dim Prym(pi)=g-1 and an E-isogeny factor has dimension one.", "Homomorphic quotient maps from a smooth curve Jacobian."),
            obligation("scalar-action quotient", "supported", "Round 25 certifies every rational P-256 endomorphism acts as a known scalar on E(F_p)[n].", "The P-256 rational prime-order subgroup."),
            obligation("rank-164 capacity threshold", "supported", "Dual-modulus complete character controls reach rank 164 only after supplying 164 abstract Prym coordinates; the first permitted genus is 165.", "Capacity, not existence of a curve."),
            obligation("explicit genus-165 cover", "unknown", "No genus-165, degree-at-least-165 cover with 164 P-256 Prym factors is constructed.", "Future homomorphic cover work."),
            obligation("non-homomorphic relation mechanism", "unknown", "The dimension theorem does not apply to a mechanism that does not factor through rational homomorphisms J(C)->E.", "Outside this screen."),
            obligation("inverse recovery and complete cost", "unknown", "No qualifying cover, fibre, relation solver, or recovery pipeline exists in the evidence.", "End-to-end P-256 attack."),
        ],
        controls: vec![
            format!("{} genus profiles were completely materialized and ranked over two coefficient fields in forward and reversed sheet order.", profiles.len()),
            format!("The full control uses {} sheets, {} character columns, and rank {} under both moduli.", full.complete_character_sheets, full.character_columns, full.rank_mod_p256_n),
            "Round 300 rank-gain zero and Round 302 complementary transport failure are exact negative controls.".into(),
            "Round 25, 298, 300, and 302 were imported only after exact byte-hash, curve, schema, and conclusion checks.".into(),
        ],
        cost_accounting: vec![
            "Measured: character-matrix rows, entries, pivot scans, swaps, inversions, multiplications, subtractions, and deterministic artifact/process telemetry.".into(),
            "Derived exactly: g>=165, Prym dimension and P-256 multiplicity >=164, degree>=165, and minimum total ramification 328 at minimum genus.".into(),
            "Unset: cover construction, rational fibre yield, relation-system degree, unsuccessful attempts, sparse linear algebra, inverse recovery, and complete attack cost.".into(),
        ],
        exploration_boundary: "Exact for relation rows obtained from same-target fibre differences through rational homomorphisms from a smooth cover Jacobian to P-256 factors; not a universal bound on non-homomorphic encodings or algorithms.".into(),
        weakest_open_obligation: "Construct and replay a genus-at-least-165, degree-at-least-165 P-256 cover whose Prym contains at least 164 rational P-256 order-n factors, or supply a demonstrably non-homomorphic mechanism outside the dimension bound.".into(),
        narrowest_supported_finding: "A homomorphic same-target cover/Jacobian route cannot meet the one-event rank-164 parity gate below genus 165; the registered genus-two split cover has zero P-256 Prym directions.".into(),
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
    let mut profiles = Vec::new();
    let mut aggregate = RankOperations::default();
    let mut peak_bytes = 0u64;
    for genus in GENUS_PROFILES {
        let profile = capacity_profile(genus, &curve.n)?;
        aggregate.add_assign(&profile.operations);
        peak_bytes = peak_bytes.max(profile.materialized_bytes_mod_p256_n);
        profiles.push(profile);
    }
    let first_parity = profiles
        .iter()
        .find(|profile| profile.parity_capacity_possible)
        .ok_or("no parity-capable profile")?;
    if first_parity.genus != 165
        || profiles
            .iter()
            .any(|profile| profile.genus < 165 && profile.parity_capacity_possible)
    {
        return Err("minimum parity genus is not 165".into());
    }
    let minimum = MinimumCapacity {
        required_homogeneous_rank: COLUMNS,
        minimum_genus: 165,
        minimum_prym_dimension: 164,
        minimum_p256_prym_multiplicity: 164,
        minimum_distinct_normalized_fibre_preimages: 165,
        minimum_geometric_cover_degree: 165,
        minimum_total_ramification_degree_at_minimum_genus: 328,
        minimum_raw_folded_preimages_under_round298_sign_model: 330,
        explicit_cover_constructed: false,
        capacity_is_attack: false,
    };
    let negative = NegativeControl {
        construction: "registered split degree-two cover H:v^2=u^6+a*u^2+b".into(),
        cover_genus: 2,
        prym_dimension: 1,
        rational_p256_order_n_prym_factors: 0,
        maximum_same_target_sheet_rank_gain: 0,
        imported_round300_maximum_quotient_rank_gain: 0,
        imported_round302_order_n_transport_exists: false,
        status: "exact negative control".into(),
    };
    let gates = Gates {
        dependency_hashes_and_schemas_checked: true,
        exact_capacity_bound_derived: true,
        dual_modulus_character_ranks_and_replays_pass: true,
        explicit_genus_165_degree_165_p256_cover: false,
        p256_prym_multiplicity_at_least_164: false,
        rank_164_same_target_fibre_replayed: false,
        explicit_forward_inverse_kernel_and_recovery: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        complete_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        projected_storage_below_2_50: false,
        discarded_probabilistic_branches_counted_as_exhaustive: false,
        promoted: false,
    };
    let assessment = build_assessment(&profiles);
    let boundaries = boundary_table(&profiles)?;
    let classification =
        "homomorphic-cover/jacobian-prym-dimension-bound/genus-165-minimum/parity-blocked";
    let obstruction = "Same-target sheet differences lie in the Prym of C->E. After Round 25 quotients scalar endomorphism transports, each rational P-256 elliptic factor supplies at most one independent order-n coordinate, so rank<=dim Prym=g-1. Rank 164 therefore requires genus at least 165, cover degree at least 165, and at least 164 P-256 Prym factors; no such cover, fibre, solver, or recovery pipeline is constructed.";
    let decision = "Reject every homomorphic cover/Jacobian capacity below genus 165 as a rho-parity route. Treat the abstract genus-165 character rank as a necessary capacity threshold, not an attack. Keep an explicit qualifying high-genus cover and genuinely non-homomorphic relation mechanisms open; do not attempt an unplanted full-depth relation.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": &dependencies,
        "conditional_capacity_statement": "same-target homomorphic sheet-row rank <= P-256 Prym multiplicity <= g-1",
        "capacity_profiles": &profiles,
        "minimum_capacity": &minimum,
        "registered_split_cover_negative_control": &negative,
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
            schema: "p256-jacobian-row-capacity-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 303,
            execution_status: "complete".into(),
            factor_base: FACTOR_BASE_ID.into(),
            registered_s17_factor_base: REGISTERED_S17_FACTOR_BASE_ID.into(),
            columns: COLUMNS,
            arity: ARITY,
            dependencies,
            conditional_capacity_statement: "For a separable cover pi:C->E whose same-target rows factor through rational homomorphisms from J(C), independent P-256 order-n sheet-row rank is at most the P-256 isogeny multiplicity in Prym(pi), at most dim Prym(pi)=genus(C)-1.".into(),
            scalar_action_import: "Round 25: rational endomorphisms of every P-256 factor act as known scalars on E(F_p)[n].".into(),
            capacity_profiles: profiles,
            minimum_capacity: minimum,
            registered_split_cover_negative_control: negative,
            aggregate_rank_operations: aggregate,
            peak_projected_character_matrix_bytes: peak_bytes,
            boundary_table: boundaries,
            structured_residual_degree_of_regularity: None,
            relations_reported_on_p256: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            exploration_boundary: "Exact for homomorphic same-target fibre differences from smooth cover Jacobians; not universal over non-homomorphic relation mechanisms.".into(),
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
        "round 303: minimum_genus={}, profiles={}, peak_matrix_bytes={}, promoted={}",
        result.minimum_capacity.minimum_genus,
        result.capacity_profiles.len(),
        result.peak_projected_character_matrix_bytes,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_jacobian_row_capacity: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collision_recurrence_matches_frozen_boundaries() {
        assert!((collision_ratio(1) - 0.964_087_797_935_000_2).abs() < 1e-14);
        assert!((collision_ratio(2) - 1.446_131_696_902_500_4).abs() < 1e-14);
        assert!((collision_ratio(164) - 13.920_747_397_073_491).abs() < 1e-13);
    }

    #[test]
    fn complete_character_controls_have_declared_rank() {
        let p61 = (BigUint::one() << 61usize) - BigUint::one();
        for columns in [1usize, 2, 4, 8, 16] {
            let (rank, _) = modular_rank(columns, &p61, false).expect("forward rank");
            let (replay, _) = modular_rank(columns, &p61, true).expect("replay rank");
            assert_eq!(rank, columns);
            assert_eq!(replay, columns);
        }
    }

    #[test]
    fn genus_165_is_first_one_event_capacity() {
        for genus in 2usize..165 {
            let cap = genus - 1;
            assert!(COLUMNS.div_ceil(cap) >= 2);
        }
        assert_eq!(COLUMNS.div_ceil(164), 1);
        assert_eq!(2 * 165 - 2, 328);
    }
}
