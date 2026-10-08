//! Exhaustive scalar-stabilizer screen for the 164-column P-256 Dickson union.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    CURVE_SLUG, DUMP_SCHEMA, WideFactorBaseDump,
};
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::CurveParams;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::json;

const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const FACTOR_BASE_SHA256: &str = "d27516ca40a612ecf3ebabfa8ae04776084c7948110e20f221da438e1e64d8f8";
const ROUND296_SHA256: &str = "52cbea224df1e1d4d6e8406f368ba413fd4e508f5309fcde65f67a98bc91c5ba";
const COLUMNS: usize = 164;
const SIGNED_POINTS: u64 = 328;
const ROOT_ORDER: u64 = 8;
const RHO_S: f64 = 1.3;

#[derive(Parser)]
#[command(about = "Exhaustively screen scalar stabilizers of the P-256 Dickson union")]
struct Cli {
    #[arg(long)]
    factor_base: PathBuf,
    #[arg(long)]
    round296: PathBuf,
    #[arg(long)]
    out: PathBuf,
}

#[derive(Clone, Serialize)]
struct Dependency {
    path: String,
    bytes: u64,
    sha256: String,
}

#[derive(Serialize)]
struct TheoremCertificate {
    subgroup_order: String,
    subgroup_order_is_registered_prime: bool,
    n_minus_one: String,
    signed_set_size: u64,
    signed_set_size_divisors: Vec<u64>,
    gcd_n_minus_one_signed_set_size: String,
    admissible_scalar_orders: Vec<u64>,
    candidate_subgroup_order: u64,
    free_action_argument: String,
    candidate_enumeration_exhaustive: bool,
}

#[derive(Serialize)]
struct InputCertificate {
    dump_schema: String,
    factor_base_id: String,
    factor_base_sha256: String,
    columns: u64,
    signed_points: u64,
    all_rows_in_canonical_pairs: bool,
    all_points_nonidentity: bool,
    all_points_on_curve: bool,
    all_sign_pairs_exact: bool,
    unique_folded_columns: bool,
}

#[derive(Clone, Serialize)]
struct MissingImage {
    source_column: u64,
    x: String,
    y: String,
}

#[derive(Clone, Serialize)]
struct RootReceipt {
    exponent: u64,
    scalar: String,
    exact_order: u64,
    images_checked: u64,
    image_hits: u64,
    unique_destinations: u64,
    all_images_in_signed_set: bool,
    destination_permutation: bool,
    stabilizes_signed_set: bool,
    first_missing_image: Option<MissingImage>,
    screen_transcript_sha256: String,
    accepted_action_sha256: Option<String>,
    accepted_action_replays: u64,
    accepted_action_replay_failures: u64,
}

#[derive(Default, Serialize)]
struct OperationCounts {
    scalar_candidates: u64,
    scalar_multiplications: u64,
    scalar_bit_rounds: u64,
    group_additions: u64,
    group_doublings: u64,
    affine_normalizations: u64,
    factor_base_lookups: u64,
    action_replay_scalar_multiplications: u64,
    action_replay_bit_rounds: u64,
    action_replay_group_additions: u64,
    action_replay_group_doublings: u64,
    action_replay_affine_normalizations: u64,
    subgroup_composition_checks: u64,
}

#[derive(Serialize)]
struct StabilizerCertificate {
    least_generator_seed: u64,
    order_eight_generator: String,
    distinct_roots: u64,
    accepted_scalars: Vec<String>,
    signed_stabilizer_order: u64,
    negation_present: bool,
    effective_folded_action_order: u64,
    accepted_scalars_form_subgroup: bool,
    signed_actions_compose_exactly: bool,
    folded_orbit_sizes: Vec<u64>,
    folded_orbits_cover_columns_once: bool,
    independent_log_quotient_dimension: u64,
}

#[derive(Clone, Serialize)]
struct BoundaryRow {
    variant: String,
    independent_log_classes: Option<u64>,
    s: f64,
    ratio_to_rho: f64,
    correctness: String,
    classification: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_checked: bool,
    theorem_enumeration_exhaustive: bool,
    every_candidate_image_checked: bool,
    accepted_actions_replay_exactly: bool,
    accepted_actions_form_exact_subgroup: bool,
    scalar_quotient_reduces_boundary: bool,
    scalar_quotient_reaches_rho: bool,
    zero_false_positives_and_false_negatives_on_complete_relation_instances: bool,
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
    screening_round: u64,
    execution_status: String,
    factor_base_dependency: Dependency,
    round296_dependency: Dependency,
    theorem: TheoremCertificate,
    input: InputCertificate,
    roots: Vec<RootReceipt>,
    stabilizer: StabilizerCertificate,
    operations: OperationCounts,
    boundary_unit: String,
    boundary_table: Vec<BoundaryRow>,
    boundary_table_sha256: String,
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

#[derive(Clone)]
struct Column {
    x: BigUint,
    low_y: BigUint,
    positive: Projective,
}

#[derive(Clone)]
struct Action {
    scalar: BigUint,
    mapping: Vec<(usize, bool)>,
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("not a non-negative integer: `{value}`"))
}

fn fixed_hex(value: &BigUint) -> String {
    format!("0x{:0>64}", value.to_str_radix(16))
}

fn fixed_32(value: &BigUint) -> Result<[u8; 32], String> {
    let bytes = value.to_bytes_be();
    if bytes.len() > 32 {
        return Err("integer exceeds 32 bytes".into());
    }
    let mut out = [0u8; 32];
    out[32 - bytes.len()..].copy_from_slice(&bytes);
    Ok(out)
}

fn checked_dependency(path: &Path, expected: &str) -> Result<(Dependency, Vec<u8>), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != expected {
        return Err(format!(
            "{} SHA-256 is {digest}, expected {expected}",
            path.display()
        ));
    }
    serde_json::from_slice::<serde_json::Value>(&bytes)
        .map_err(|error| format!("{}: {error}", path.display()))?;
    Ok((
        Dependency {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            sha256: digest,
        },
        bytes,
    ))
}

fn divisors(value: u64) -> Vec<u64> {
    let mut out = Vec::new();
    for candidate in 1..=value {
        if value.is_multiple_of(candidate) {
            out.push(candidate);
        }
    }
    out
}

fn exact_order(value: &BigUint, modulus: &BigUint, maximum: u64) -> Result<u64, String> {
    if value.is_zero() || value >= modulus {
        return Err("scalar is outside F_n^*".into());
    }
    for order in divisors(maximum) {
        if value.modpow(&BigUint::from(order), modulus).is_one() {
            return Ok(order);
        }
    }
    Err(format!("scalar does not have order dividing {maximum}"))
}

fn root_candidates(n: &BigUint) -> Result<(u64, BigUint, Vec<BigUint>), String> {
    let exponent = (n - BigUint::one()) / BigUint::from(ROOT_ORDER);
    let mut seed = 2u64;
    let generator = loop {
        let candidate = BigUint::from(seed).modpow(&exponent, n);
        if exact_order(&candidate, n, ROOT_ORDER)? == ROOT_ORDER {
            break candidate;
        }
        seed = seed.checked_add(1).ok_or("root seed overflow")?;
    };
    let mut roots = Vec::new();
    let mut value = BigUint::one();
    for _ in 0..ROOT_ORDER {
        roots.push(value.clone());
        value = (value * &generator) % n;
    }
    let distinct = roots.iter().cloned().collect::<BTreeSet<_>>();
    if distinct.len() != ROOT_ORDER as usize || !value.is_one() {
        return Err("eighth-root enumeration is not an exact cycle".into());
    }
    Ok((seed, generator, roots))
}

fn load_columns(dump: &WideFactorBaseDump) -> Result<Vec<Column>, String> {
    let curve = CurveParams::p256();
    if dump.schema != DUMP_SCHEMA || dump.curve.slug != CURVE_SLUG {
        return Err("factor-base dump schema or curve mismatch".into());
    }
    if dump.factor_base.fb_id != FACTOR_BASE_ID
        || dump.factor_base.columns != COLUMNS as u64
        || dump.factor_base.signed_points != SIGNED_POINTS
        || dump.points.len() != SIGNED_POINTS as usize
    {
        return Err("factor-base identity or cardinality mismatch".into());
    }
    let mut rows = BTreeMap::<u64, Vec<_>>::new();
    for row in &dump.points {
        rows.entry(row.col).or_default().push(row);
    }
    let mut columns = Vec::with_capacity(COLUMNS);
    let negative_one = (&curve.n - BigUint::one()).to_string();
    let mut seen_x = BTreeSet::new();
    for col in 0..COLUMNS {
        let pair = rows
            .get(&(col as u64))
            .ok_or_else(|| format!("missing column {col}"))?;
        if pair.len() != 2 {
            return Err(format!("column {col} has {} signed rows", pair.len()));
        }
        let positive = pair
            .iter()
            .find(|row| row.coef == "1")
            .ok_or_else(|| format!("column {col} has no +1 row"))?;
        let negative = pair
            .iter()
            .find(|row| row.coef == negative_one)
            .ok_or_else(|| format!("column {col} has no -1 row"))?;
        let x = parse_big(&positive.x)?;
        let low_y = parse_big(&positive.y)?;
        let negative_x = parse_big(&negative.x)?;
        let high_y = parse_big(&negative.y)?;
        if x != negative_x || (&low_y + &high_y) % &curve.p != BigUint::zero() {
            return Err(format!("column {col} is not an exact sign pair"));
        }
        if x >= curve.p || low_y >= curve.p || high_y >= curve.p || low_y >= high_y {
            return Err(format!("column {col} is not in canonical low-y form"));
        }
        if !seen_x.insert(x.clone()) {
            return Err(format!("duplicate folded abscissa in column {col}"));
        }
        let textbook = crypto_lib::ecc::Point::Affine {
            x: crypto_lib::ecc::FieldElement::new(x.clone(), curve.p.clone()),
            y: crypto_lib::ecc::FieldElement::new(low_y.clone(), curve.p.clone()),
        };
        if !curve.is_on_curve(&textbook) {
            return Err(format!("column {col} is off curve"));
        }
        columns.push(Column {
            x: x.clone(),
            low_y: low_y.clone(),
            positive: Projective::from_affine(&U256::from_biguint(&x), &U256::from_biguint(&low_y)),
        });
    }
    if rows.len() != COLUMNS {
        return Err("factor base contains an out-of-range column".into());
    }
    Ok(columns)
}

fn action_hash(scalar: &BigUint, mapping: &[(usize, bool)]) -> Result<String, String> {
    let mut stream = Vec::with_capacity(32 + mapping.len() * 9);
    stream.extend(fixed_32(scalar)?);
    for (source, (destination, high_sign)) in mapping.iter().enumerate() {
        stream.extend((source as u32).to_be_bytes());
        stream.extend((*destination as u32).to_be_bytes());
        stream.push(u8::from(*high_sign));
    }
    Ok(sha256_hex(&stream))
}

fn screen_root(
    exponent: u64,
    scalar: &BigUint,
    order: u64,
    columns: &[Column],
    curve: &CurveParams,
    operations: &mut OperationCounts,
) -> Result<(RootReceipt, Option<Action>), String> {
    let lookup = columns
        .iter()
        .enumerate()
        .map(|(index, column)| (column.x.clone(), index))
        .collect::<BTreeMap<_, _>>();
    let mut mapping = Vec::with_capacity(columns.len());
    let mut first_missing = None;
    let mut hits = 0u64;
    let mut destinations = BTreeSet::new();
    let mut transcript = Vec::new();
    transcript.extend(fixed_32(scalar)?);
    for (source, column) in columns.iter().enumerate() {
        let image = column.positive.scalar_mul_ct(scalar, curve.order_bits());
        operations.scalar_multiplications += 1;
        operations.scalar_bit_rounds += curve.order_bits() as u64;
        operations.group_additions += curve.order_bits() as u64;
        operations.group_doublings += 2 * curve.order_bits() as u64;
        operations.affine_normalizations += 1;
        operations.factor_base_lookups += 1;
        let (x_raw, y_raw) = image
            .to_affine()
            .ok_or_else(|| format!("nonzero scalar mapped column {source} to identity"))?;
        let x = x_raw.to_biguint();
        let y = y_raw.to_biguint();
        transcript.extend((source as u32).to_be_bytes());
        transcript.extend(fixed_32(&x)?);
        transcript.extend(fixed_32(&y)?);
        if let Some(&destination) = lookup.get(&x) {
            let target = &columns[destination];
            let high_sign = if y == target.low_y {
                false
            } else if y == (&curve.p - &target.low_y) {
                true
            } else {
                return Err(format!(
                    "column {source} image has a known x but neither y sign"
                ));
            };
            transcript.push(1);
            transcript.extend((destination as u32).to_be_bytes());
            transcript.push(u8::from(high_sign));
            mapping.push((destination, high_sign));
            destinations.insert(destination);
            hits += 1;
        } else {
            transcript.push(0);
            if first_missing.is_none() {
                first_missing = Some(MissingImage {
                    source_column: source as u64,
                    x: fixed_hex(&x),
                    y: fixed_hex(&y),
                });
            }
        }
    }
    let all_images = hits == columns.len() as u64;
    let permutation = all_images && destinations.len() == columns.len();
    let accepted = all_images && permutation;
    let accepted_hash = accepted
        .then(|| action_hash(scalar, &mapping))
        .transpose()?;
    let action = accepted.then(|| Action {
        scalar: scalar.clone(),
        mapping,
    });
    Ok((
        RootReceipt {
            exponent,
            scalar: fixed_hex(scalar),
            exact_order: order,
            images_checked: columns.len() as u64,
            image_hits: hits,
            unique_destinations: destinations.len() as u64,
            all_images_in_signed_set: all_images,
            destination_permutation: permutation,
            stabilizes_signed_set: accepted,
            first_missing_image: first_missing,
            screen_transcript_sha256: sha256_hex(&transcript),
            accepted_action_sha256: accepted_hash,
            accepted_action_replays: 0,
            accepted_action_replay_failures: 0,
        },
        action,
    ))
}

fn replay_action(
    action: &Action,
    columns: &[Column],
    curve: &CurveParams,
    operations: &mut OperationCounts,
) -> Result<(u64, u64), String> {
    let mut replays = 0u64;
    let mut failures = 0u64;
    for (source, &(destination, high_sign)) in action.mapping.iter().enumerate() {
        let image = columns[source]
            .positive
            .scalar_mul_ct(&action.scalar, curve.order_bits());
        operations.action_replay_scalar_multiplications += 1;
        operations.action_replay_bit_rounds += curve.order_bits() as u64;
        operations.action_replay_group_additions += curve.order_bits() as u64;
        operations.action_replay_group_doublings += 2 * curve.order_bits() as u64;
        operations.action_replay_affine_normalizations += 1;
        let expected = if high_sign {
            columns[destination].positive.neg()
        } else {
            columns[destination].positive
        };
        let image_affine = image
            .to_affine()
            .map(|(x, y)| (x.to_biguint(), y.to_biguint()));
        let expected_affine = expected
            .to_affine()
            .map(|(x, y)| (x.to_biguint(), y.to_biguint()));
        if image_affine != expected_affine {
            failures += 1;
        }
        replays += 1;
    }
    Ok((replays, failures))
}

fn certify_actions(
    actions: &[Action],
    n: &BigUint,
    operations: &mut OperationCounts,
) -> Result<(bool, bool, Vec<u64>), String> {
    let index = actions
        .iter()
        .enumerate()
        .map(|(position, action)| (action.scalar.clone(), position))
        .collect::<BTreeMap<_, _>>();
    let mut subgroup = true;
    let mut composition = true;
    for left in actions {
        for right in actions {
            operations.subgroup_composition_checks += 1;
            let product = (&left.scalar * &right.scalar) % n;
            let Some(&product_index) = index.get(&product) else {
                subgroup = false;
                composition = false;
                continue;
            };
            let product_action = &actions[product_index];
            for source in 0..left.mapping.len() {
                let (middle, first_sign) = left.mapping[source];
                let (destination, second_sign) = right.mapping[middle];
                if product_action.mapping[source] != (destination, first_sign ^ second_sign) {
                    composition = false;
                }
            }
        }
    }

    let mut seen = [false; COLUMNS];
    let mut orbit_sizes = Vec::new();
    for start in 0..COLUMNS {
        if seen[start] {
            continue;
        }
        let orbit = actions
            .iter()
            .map(|action| action.mapping[start].0)
            .collect::<BTreeSet<_>>();
        if orbit.is_empty() {
            return Err("empty folded orbit".into());
        }
        for &column in &orbit {
            if seen[column] {
                return Err("folded scalar orbits overlap".into());
            }
            seen[column] = true;
        }
        orbit_sizes.push(orbit.len() as u64);
    }
    orbit_sizes.sort_unstable();
    Ok((subgroup, composition, orbit_sizes))
}

fn gamma_ratio(collisions: u64) -> f64 {
    assert!(collisions >= 1);
    let mut ratio = std::f64::consts::PI.sqrt() / 2.0;
    for k in 1..collisions {
        ratio *= (k as f64 + 0.5) / k as f64;
    }
    ratio
}

fn collision_s(classes: u64) -> f64 {
    2.0_f64.sqrt() * gamma_ratio(classes)
}

fn boundary_row(variant: &str, classes: u64, classification: &str) -> BoundaryRow {
    let s = collision_s(classes);
    BoundaryRow {
        variant: variant.into(),
        independent_log_classes: Some(classes),
        s,
        ratio_to_rho: s / RHO_S,
        correctness: "exact favourable lower boundary; no relation solver".into(),
        classification: classification.into(),
    }
}

fn boundary_digest(rows: &[BoundaryRow]) -> Result<String, String> {
    serde_json::to_vec(rows)
        .map(|bytes| sha256_hex(&bytes))
        .map_err(|error| error.to_string())
}

fn build_result(cli: &Cli) -> Result<ResultReceipt, String> {
    let (factor_base_dependency, factor_base_bytes) =
        checked_dependency(&cli.factor_base, FACTOR_BASE_SHA256)?;
    let (round296_dependency, _) = checked_dependency(&cli.round296, ROUND296_SHA256)?;
    let dump: WideFactorBaseDump =
        serde_json::from_slice(&factor_base_bytes).map_err(|error| error.to_string())?;
    let columns = load_columns(&dump)?;
    let curve = CurveParams::p256();
    let n_minus_one = &curve.n - BigUint::one();
    let set_divisors = divisors(SIGNED_POINTS);
    let admissible_orders = set_divisors
        .iter()
        .copied()
        .filter(|order| (&n_minus_one % BigUint::from(*order)).is_zero())
        .collect::<Vec<_>>();
    let gcd = n_minus_one.gcd(&BigUint::from(SIGNED_POINTS));
    if gcd != BigUint::from(ROOT_ORDER) || admissible_orders != vec![1, 2, 4, 8] {
        return Err("unexpected scalar-order divisibility certificate".into());
    }
    let (least_seed, generator, root_values) = root_candidates(&curve.n)?;
    let root_orders = root_values
        .iter()
        .map(|root| exact_order(root, &curve.n, ROOT_ORDER))
        .collect::<Result<Vec<_>, _>>()?;
    if root_orders.iter().copied().collect::<BTreeSet<_>>()
        != admissible_orders.iter().copied().collect()
    {
        return Err("root enumeration does not cover every admissible order".into());
    }

    let mut operations = OperationCounts {
        scalar_candidates: ROOT_ORDER,
        ..OperationCounts::default()
    };
    let mut roots = Vec::new();
    let mut actions = Vec::new();
    for (exponent, (scalar, order)) in root_values.iter().zip(root_orders).enumerate() {
        let (receipt, action) = screen_root(
            exponent as u64,
            scalar,
            order,
            &columns,
            &curve,
            &mut operations,
        )?;
        roots.push(receipt);
        if let Some(action) = action {
            actions.push(action);
        }
    }
    for action in &actions {
        let position = root_values
            .iter()
            .position(|scalar| scalar == &action.scalar)
            .ok_or("accepted scalar missing from root enumeration")?;
        let (replays, failures) = replay_action(action, &columns, &curve, &mut operations)?;
        roots[position].accepted_action_replays = replays;
        roots[position].accepted_action_replay_failures = failures;
    }
    let replay_failures = roots
        .iter()
        .map(|root| root.accepted_action_replay_failures)
        .sum::<u64>();
    let (subgroup, composition, orbit_sizes) =
        certify_actions(&actions, &curve.n, &mut operations)?;
    let orbits_cover = orbit_sizes.iter().sum::<u64>() == COLUMNS as u64;
    let quotient_dimension = orbit_sizes.len() as u64;
    let negative_one = &curve.n - BigUint::one();
    let negation_present = actions.iter().any(|action| action.scalar == negative_one);
    if !negation_present || actions.is_empty() || actions.len() % 2 != 0 {
        return Err("accepted scalar set lacks the required negation subgroup".into());
    }
    let effective_order = actions.len() as u64 / 2;
    let accepted_scalars = actions
        .iter()
        .map(|action| fixed_hex(&action.scalar))
        .collect::<Vec<_>>();
    let original = boundary_row("Round-295 unquotiented base", COLUMNS as u64, "control");
    let quotient = boundary_row(
        "Round-297 exact scalar quotient",
        quotient_dimension,
        if quotient_dimension < COLUMNS as u64 {
            "advance"
        } else {
            "negative"
        },
    );
    let boundary_reduced = quotient.ratio_to_rho < original.ratio_to_rho;
    let at_rho = quotient.ratio_to_rho <= 1.0;
    let boundary_table = vec![
        BoundaryRow {
            variant: "Pollard rho".into(),
            independent_log_classes: None,
            s: RHO_S,
            ratio_to_rho: 1.0,
            correctness: "frozen reference".into(),
            classification: "reference".into(),
        },
        original,
        quotient,
    ];
    let boundary_table_sha256 = boundary_digest(&boundary_table)?;
    let exact = replay_failures == 0 && subgroup && composition && orbits_cover;
    let theorem = TheoremCertificate {
        subgroup_order: curve.n.to_string(),
        subgroup_order_is_registered_prime: true,
        n_minus_one: n_minus_one.to_string(),
        signed_set_size: SIGNED_POINTS,
        signed_set_size_divisors: set_divisors,
        gcd_n_minus_one_signed_set_size: gcd.to_string(),
        admissible_scalar_orders: admissible_orders,
        candidate_subgroup_order: ROOT_ORDER,
        free_action_argument: "For nonzero P in the prime-order subgroup, [a]^j P=P iff a^j=1 mod n; therefore a scalar stabilizing a finite nonidentity point set acts in free cycles of length ord_n(a), which must divide 328.".into(),
        candidate_enumeration_exhaustive: true,
    };
    let input = InputCertificate {
        dump_schema: dump.schema.clone(),
        factor_base_id: dump.factor_base.fb_id.clone(),
        factor_base_sha256: dump.factor_base.fb_sha256.clone(),
        columns: columns.len() as u64,
        signed_points: dump.points.len() as u64,
        all_rows_in_canonical_pairs: true,
        all_points_nonidentity: true,
        all_points_on_curve: true,
        all_sign_pairs_exact: true,
        unique_folded_columns: true,
    };
    let stabilizer = StabilizerCertificate {
        least_generator_seed: least_seed,
        order_eight_generator: fixed_hex(&generator),
        distinct_roots: root_values.len() as u64,
        accepted_scalars,
        signed_stabilizer_order: actions.len() as u64,
        negation_present,
        effective_folded_action_order: effective_order,
        accepted_scalars_form_subgroup: subgroup,
        signed_actions_compose_exactly: composition,
        folded_orbit_sizes: orbit_sizes,
        folded_orbits_cover_columns_once: orbits_cover,
        independent_log_quotient_dimension: quotient_dimension,
    };
    let gates = Gates {
        dependency_hashes_checked: true,
        theorem_enumeration_exhaustive: true,
        every_candidate_image_checked: operations.scalar_multiplications
            == ROOT_ORDER * COLUMNS as u64,
        accepted_actions_replay_exactly: replay_failures == 0,
        accepted_actions_form_exact_subgroup: subgroup && composition,
        scalar_quotient_reduces_boundary: boundary_reduced,
        scalar_quotient_reaches_rho: at_rho,
        zero_false_positives_and_false_negatives_on_complete_relation_instances: false,
        structured_residual_degree_at_most_5: false,
        per_usable_relation_below_2_103: false,
        complete_cost_below_2_120: false,
        projected_materialized_storage_below_2_50: false,
        promoted: false,
    };
    if !exact {
        return Err("scalar action certificate failed exact replay or composition".into());
    }
    let classification = if boundary_reduced {
        "exact-scalar-quotient-advance/promotion-negative"
    } else {
        "exhaustive-negative/promotion-negative"
    };
    let dominant_obstruction = if quotient_dimension == COLUMNS as u64 {
        "The exhaustive scalar stabilizer is only {+1,-1}; negation was already folded, so all 164 factor-base logarithms remain independent and the exact favourable collection boundary stays 13.920747 times rho."
    } else {
        "A nontrivial folded scalar quotient exists, but its exact favourable collection boundary remains above rho and no complete relation solver or end-to-end cost is established."
    };
    let decision = if boundary_reduced {
        "Retain the exact quotient as structural evidence, but do not attempt an unplanted P-256 relation because the rho and solver promotion gates fail."
    } else {
        "Close scalar set-stabilizer transport for this factor base. Do not attempt an unplanted P-256 relation; a successor needs a non-scalar correspondence or a different globally symmetric base."
    };

    let semantic = json!({
        "curve": CURVE_SLUG,
        "factor_base_sha256": factor_base_dependency.sha256,
        "round296_sha256": round296_dependency.sha256,
        "theorem": theorem,
        "input": input,
        "roots": roots,
        "stabilizer": stabilizer,
        "operations": operations,
        "boundary_table": boundary_table,
        "gates": gates,
        "classification": classification,
        "dominant_obstruction": dominant_obstruction,
        "decision": decision,
    });
    let semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).map_err(|error| error.to_string())?);
    Ok(ResultReceipt {
        schema: "p256-union-scalar-stabilizer/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 297,
        execution_status: "complete".into(),
        factor_base_dependency,
        round296_dependency,
        theorem,
        input,
        roots,
        stabilizer,
        operations,
        boundary_unit: "favourable P-256 group-addition equivalents / sqrt(n)".into(),
        boundary_table,
        boundary_table_sha256,
        relations_reported: 0,
        full_depth_unplanted_relation_attempted: false,
        process_telemetry_recorded_in_isolation_receipt: true,
        gates,
        classification: classification.into(),
        dominant_obstruction: dominant_obstruction.into(),
        decision: decision.into(),
        semantic_evidence_sha256,
        result_json_bytes: 0,
    })
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

fn run(cli: Cli) -> Result<(), String> {
    let mut result = build_result(&cli)?;
    write_result(&cli.out, &mut result)?;
    eprintln!(
        "round 297: signed stabilizer {}, folded K={}, boundary {:.6}x rho -> {}",
        result.stabilizer.signed_stabilizer_order,
        result.stabilizer.independent_log_quotient_dimension,
        result
            .boundary_table
            .last()
            .expect("quotient row")
            .ratio_to_rho,
        result.classification
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_union_scalar_stabilizer: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_set_divisor_intersection_is_exactly_eight() {
        let n_minus_one = CurveParams::p256().n - BigUint::one();
        assert_eq!(
            n_minus_one.gcd(&BigUint::from(SIGNED_POINTS)),
            BigUint::from(8u8)
        );
        let admissible = divisors(SIGNED_POINTS)
            .into_iter()
            .filter(|order| (&n_minus_one % BigUint::from(*order)).is_zero())
            .collect::<Vec<_>>();
        assert_eq!(admissible, vec![1, 2, 4, 8]);
    }

    #[test]
    fn eighth_roots_are_distinct_and_cover_admissible_orders() {
        let curve = CurveParams::p256();
        let (_, generator, roots) = root_candidates(&curve.n).expect("roots");
        assert_eq!(
            exact_order(&generator, &curve.n, ROOT_ORDER).expect("order"),
            8
        );
        assert_eq!(roots.iter().cloned().collect::<BTreeSet<_>>().len(), 8);
        let orders = roots
            .iter()
            .map(|root| exact_order(root, &curve.n, ROOT_ORDER).expect("order"))
            .collect::<BTreeSet<_>>();
        assert_eq!(orders, BTreeSet::from([1, 2, 4, 8]));
    }

    #[test]
    fn frozen_k164_boundary_matches_round295() {
        let ratio = collision_s(164) / RHO_S;
        assert!((ratio - 13.920747397073491).abs() < 1e-12);
    }
}
