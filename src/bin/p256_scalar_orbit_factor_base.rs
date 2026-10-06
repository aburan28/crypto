//! Exact P-256 endomorphism and scalar-orbit factor-base screen (round 25).

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ec_index_calculus::sqrt_mod_p;
use crypto_lib::cryptanalysis::ecbench::canonical::short_id;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    WideCurveFacts, WideFactorBaseDump, WideFactorBaseFacts, WideFbPoint, CURVE_EC1, CURVE_ICV1,
    CURVE_SLUG, CURVE_UID, DUMP_SCHEMA, FACTOR_BASE_SCHEMA, POINT_KEY_BYTES, POINT_KEY_ENCODING,
};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::p256_point::P256ProjectivePoint as Projective;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::json;

const ROUND24_SHA256: &str = "d1e2fe33e1abbae8a9bf6885ca7cfef8012233363f4cd43e0d30d468a381fc9c";
const TARGET_COLUMNS: u64 = 131_458;
const RHO_S: f64 = 1.3;
const ENTRY_BYTES: u64 = 64;
const FIXED_BASE_WINDOW_BITS: usize = 8;

#[derive(Parser)]
#[command(about = "Screen exact P-256 endomorphism and scalar-orbit factor bases")]
struct Cli {
    #[arg(long)]
    round24: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
    /// Optional native `ecbench.factor_base_dump/v1-wide` materialisation.
    #[arg(long)]
    dump: Option<PathBuf>,
}

#[derive(Clone)]
struct FactorPower {
    prime: BigUint,
    exponent: u32,
}

#[derive(Serialize)]
struct FactorPowerRecord {
    prime: String,
    exponent: u32,
}

#[derive(Serialize)]
struct PrimeCertificate {
    n: String,
    n_minus_one_factorization: Vec<FactorPowerRecord>,
    witness: String,
    criterion: String,
    verified: bool,
}

#[derive(Serialize)]
struct CmScreen {
    trace: String,
    discriminant: String,
    discriminant_abs: String,
    factorization: Vec<FactorPowerRecord>,
    factorization_product_verified: bool,
    deterministic_u64_prime_factors: Vec<String>,
    large_prime_certificates: Vec<PrimeCertificate>,
    all_factors_distinct_prime: bool,
    fundamental_discriminant: bool,
    frobenius_conductor: String,
    endomorphism_ring: String,
    minimum_noninteger_endomorphism_degree: String,
    minimum_noninteger_endomorphism_degree_log2: f64,
    rational_point_action: String,
    low_degree_noninteger_transport_exists: bool,
}

#[derive(Serialize)]
struct WidthCandidate {
    order: String,
    negation_folded_columns: String,
    absolute_distance_from_target: String,
}

#[derive(Serialize)]
struct ScalarOrderScreen {
    n_minus_one_factorization: Vec<FactorPowerRecord>,
    factorization_product_verified: bool,
    large_factor_certificate: PrimeCertificate,
    n_prime_certificate: PrimeCertificate,
    least_primitive_root: String,
    primitive_root_verified: bool,
    divisors_enumerated: u64,
    enumeration_exhaustive: bool,
    target_columns: u64,
    closest_widths: Vec<WidthCandidate>,
    selected_order: String,
    selected_negation_folded_columns: u64,
    selected_width_delta: u64,
    exact_order_generators_screened: u64,
    selected_scalar: String,
    selected_scalar_signed_binary_bits: u64,
    selected_scalar_signed_binary_weight: u64,
    selected_scalar_binary_operations: u64,
    selected_scalar_exact_order_verified: bool,
}

#[derive(Default, Serialize)]
struct GroupOperations {
    fixed_base_table_additions: u64,
    fixed_base_table_doublings: u64,
    fixed_base_output_additions: u64,
    transport_replay_additions: u64,
    transport_replay_doublings: u64,
    normalization_field_multiplications: u64,
    normalization_inversions: u64,
}

#[derive(Serialize)]
struct NativeFactorBase {
    fb_id: String,
    fb_sha256: String,
    family: String,
    columns: u64,
    signed_points: u64,
    points_sha256: String,
    point_key_encoding: String,
    anchor_x: String,
    anchor_y: String,
    anchor_derivation_counter: u64,
    anchor_construction_log_recorded: bool,
    all_points_nonidentity: bool,
    all_points_on_curve: bool,
    unique_up_to_sign: bool,
    adjacent_transport_relations: u64,
    wrap_relation_uses_negation: bool,
    transport_replay_failures: u64,
    all_reported_transport_relations_replayed: bool,
    independent_log_quotient_dimension: u64,
    ordinary_zero_relations_supply_anchor_information: bool,
    dump_schema: String,
    dump_bytes: u64,
    dump_sha256: String,
    rebuild_identity_verified: bool,
    operations: GroupOperations,
}

#[derive(Serialize)]
struct Projection {
    case: String,
    required_target_collisions: u64,
    disjoint_probability: f64,
    oracle_s: f64,
    oracle_ratio_to_rho: f64,
    log2_oracle_samples: f64,
    direct_additions_per_sample: f64,
    log2_direct_group_additions: f64,
    direct_ratio_to_rho: f64,
    log2_one_list_bytes: f64,
    below_rho: bool,
    classification: String,
}

#[derive(Serialize)]
struct DegreeScreen {
    exact_support_polynomial_degree: u64,
    multiplication_map_degree: String,
    multiplication_map_degree_log2: f64,
    short_scalar_evaluation_chain_is_low_degree_membership: bool,
    structured_residual_degree_of_regularity: Option<u32>,
    degree_gate_passed: bool,
    classification: String,
}

#[derive(Serialize)]
struct PromotionGates {
    zero_false_positives: bool,
    zero_false_negatives: bool,
    exact_transport_replay: bool,
    structured_residual_degree_at_most_5: bool,
    relation_collection_below_2_120: bool,
    per_usable_row_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    complete_cost_below_rho: bool,
    selector_and_anchor_non_generic: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    relation_model: String,
    dependency_sha256: String,
    external_factor_discovery_trusted: bool,
    cm_screen: CmScreen,
    scalar_order_screen: ScalarOrderScreen,
    factor_base: NativeFactorBase,
    projections: Vec<Projection>,
    degree_screen: DegreeScreen,
    false_positives: u64,
    false_negatives: u64,
    promotion_gates: PromotionGates,
    full_depth_unplanted_attempted: bool,
    decision: String,
}

#[derive(Clone)]
struct OrbitPoint {
    orbit_index: u64,
    point: Projective,
    x: BigUint,
    low_y: BigUint,
}

fn parse(value: &str) -> BigUint {
    BigUint::parse_bytes(value.as_bytes(), 10).expect("frozen decimal integer")
}

fn factor(prime: &str, exponent: u32) -> FactorPower {
    FactorPower {
        prime: parse(prime),
        exponent,
    }
}

fn records(factors: &[FactorPower]) -> Vec<FactorPowerRecord> {
    factors
        .iter()
        .map(|entry| FactorPowerRecord {
            prime: entry.prime.to_string(),
            exponent: entry.exponent,
        })
        .collect()
}

fn factor_product(factors: &[FactorPower]) -> BigUint {
    factors.iter().fold(BigUint::one(), |product, entry| {
        product * entry.prime.pow(entry.exponent)
    })
}

fn unique_primes(factors: &[FactorPower]) -> Vec<BigUint> {
    factors.iter().map(|entry| entry.prime.clone()).collect()
}

fn lucas_witness(n: &BigUint, factors: &[FactorPower]) -> Result<BigUint, String> {
    if factor_product(factors) != n - 1u8 {
        return Err(format!("factorization of {} - 1 is not exact", n));
    }
    for candidate in 2u64..1_000_000 {
        let a = BigUint::from(candidate);
        if a.modpow(&(n - 1u8), n) != BigUint::one() {
            continue;
        }
        if unique_primes(factors).iter().all(|prime| {
            let power = a.modpow(&((n - 1u8) / prime), n);
            let difference = if power.is_zero() {
                n - 1u8
            } else {
                power - 1u8
            };
            difference.gcd(n).is_one()
        }) {
            return Ok(a);
        }
    }
    Err(format!("no Lucas witness below one million for {n}"))
}

fn certificate(n: &BigUint, factors: &[FactorPower]) -> Result<PrimeCertificate, String> {
    let witness = lucas_witness(n, factors)?;
    Ok(PrimeCertificate {
        n: n.to_string(),
        n_minus_one_factorization: records(factors),
        witness: witness.to_string(),
        criterion: "Lucas: a^(n-1)=1 and gcd(a^((n-1)/q)-1,n)=1 for every prime q | n-1".into(),
        verified: true,
    })
}

fn discriminant_factors() -> Vec<FactorPower> {
    vec![
        factor("3", 1),
        factor("5", 1),
        factor("456597257999", 1),
        factor("1428624589419343516204097", 1),
        factor("46523541035814968339936406074986559003387", 1),
    ]
}

fn n_minus_one_factors() -> Vec<FactorPower> {
    vec![
        factor("2", 4),
        factor("3", 1),
        factor("71", 1),
        factor("131", 1),
        factor("373", 1),
        factor("3407", 1),
        factor("17449", 1),
        factor("38189", 1),
        factor("187019741", 1),
        factor("622491383", 1),
        factor("1002328039319", 1),
        factor("2624747550333869278416773953", 1),
    ]
}

fn large_prime_certificates() -> Result<Vec<PrimeCertificate>, String> {
    let recursive = parse("1869236796843064056413");
    let recursive_factors = vec![
        factor("2", 2),
        factor("23", 1),
        factor("3911", 1),
        factor("5195037399650551", 1),
    ];
    let first = parse("1428624589419343516204097");
    let first_factors = vec![
        factor("2", 6),
        factor("1669", 1),
        factor("13374631042347059581", 1),
    ];
    let second = parse("46523541035814968339936406074986559003387");
    let second_factors = vec![
        factor("2", 1),
        factor("269", 1),
        factor("643", 1),
        factor("215531", 1),
        factor("333814693", 1),
        factor("1869236796843064056413", 1),
    ];
    Ok(vec![
        certificate(&recursive, &recursive_factors)?,
        certificate(&first, &first_factors)?,
        certificate(&second, &second_factors)?,
    ])
}

fn n_large_factor_certificate() -> Result<PrimeCertificate, String> {
    let n = parse("2624747550333869278416773953");
    let factors = vec![
        factor("2", 6),
        factor("3", 2),
        factor("1297", 1),
        factor("16879", 1),
        factor("208150935158385979", 1),
    ];
    certificate(&n, &factors)
}

fn verify_u64_primes(values: &[&str]) -> Result<Vec<String>, String> {
    let mut out = Vec::new();
    for value in values {
        let parsed = value
            .parse::<u64>()
            .map_err(|_| format!("{value} does not fit u64"))?;
        if !is_prime_u64(parsed) {
            return Err(format!("{value} failed deterministic u64 primality"));
        }
        out.push((*value).to_string());
    }
    Ok(out)
}

fn cm_screen(curve: &CurveParams) -> Result<CmScreen, String> {
    let trace = &curve.p + 1u8 - &curve.n;
    let discriminant_abs = BigUint::from(4u8) * &curve.p - &trace * &trace;
    let factors = discriminant_factors();
    let factorization_product_verified = factor_product(&factors) == discriminant_abs;
    if !factorization_product_verified {
        return Err("candidate CM factorization product mismatch".into());
    }
    let deterministic_u64_prime_factors = verify_u64_primes(&[
        "3",
        "5",
        "456597257999",
        "1669",
        "13374631042347059581",
        "269",
        "643",
        "215531",
        "333814693",
        "23",
        "3911",
        "5195037399650551",
    ])?;
    let certificates = large_prime_certificates()?;
    let distinct_factor_count = factors
        .iter()
        .map(|entry| entry.prime.clone())
        .collect::<BTreeSet<_>>()
        .len();
    let all_factors_distinct_prime = distinct_factor_count == factors.len();
    let fundamental = (&discriminant_abs % 4u8 == BigUint::from(3u8))
        && all_factors_distinct_prime
        && factors.iter().all(|entry| entry.exponent == 1);
    if !fundamental {
        return Err("CM discriminant did not verify as negative fundamental".into());
    }
    let minimum_degree = (&discriminant_abs + 1u8) / 4u8;
    let log2 = minimum_degree
        .to_f64()
        .ok_or("minimum degree did not convert to f64")?
        .log2();
    Ok(CmScreen {
        trace: trace.to_string(),
        discriminant: format!("-{}", discriminant_abs),
        discriminant_abs: discriminant_abs.to_string(),
        factorization: records(&factors),
        factorization_product_verified,
        deterministic_u64_prime_factors,
        large_prime_certificates: certificates,
        all_factors_distinct_prime,
        fundamental_discriminant: true,
        frobenius_conductor: "1".into(),
        endomorphism_ring: "End(E)=Z[pi]=O_K (maximal order)".into(),
        minimum_noninteger_endomorphism_degree: minimum_degree.to_string(),
        minimum_noninteger_endomorphism_degree_log2: log2,
        rational_point_action: "Every alpha=a+b*pi restricts on E(F_p) to the known scalar [a+b], because pi(P)=P; noninteger CM maps supply no non-generic rational-point action.".into(),
        low_degree_noninteger_transport_exists: false,
    })
}

fn divisors(factors: &[FactorPower]) -> Vec<BigUint> {
    let mut divisors = vec![BigUint::one()];
    for entry in factors {
        let old = divisors.clone();
        let mut power = BigUint::one();
        for _ in 0..entry.exponent {
            power *= &entry.prime;
            divisors.extend(old.iter().map(|value| value * &power));
        }
    }
    divisors
}

fn folded_width(order: &BigUint) -> BigUint {
    if order.is_even() {
        order >> 1usize
    } else {
        order.clone()
    }
}

fn absolute_difference(left: &BigUint, right: &BigUint) -> BigUint {
    if left >= right {
        left - right
    } else {
        right - left
    }
}

fn bit_weight(value: &BigUint) -> u64 {
    value
        .to_bytes_le()
        .iter()
        .map(|byte| u64::from(byte.count_ones()))
        .sum()
}

fn select_scalar(
    subgroup_root: &BigUint,
    order: u64,
    n: &BigUint,
) -> Result<(BigUint, u64), String> {
    let mut current = BigUint::one();
    let mut best: Option<(u64, BigUint)> = None;
    let mut generators = 0u64;
    for exponent in 1..order {
        current = (&current * subgroup_root) % n;
        if exponent.gcd(&order) != 1 {
            continue;
        }
        generators += 1;
        let negative = n - &current;
        let signed = if current <= negative {
            current.clone()
        } else {
            negative
        };
        let cost = signed.bits() + bit_weight(&signed);
        let replace = best
            .as_ref()
            .map(|(best_cost, best_scalar)| {
                cost < *best_cost || (cost == *best_cost && signed < *best_scalar)
            })
            .unwrap_or(true);
        if replace {
            best = Some((cost, signed));
        }
    }
    let (_, scalar) = best.ok_or("selected order has no generators")?;
    Ok((scalar, generators))
}

fn exact_order(value: &BigUint, order: &BigUint, primes: &[BigUint], modulus: &BigUint) -> bool {
    value.modpow(order, modulus) == BigUint::one()
        && primes
            .iter()
            .all(|prime| value.modpow(&(order / prime), modulus) != BigUint::one())
}

fn scalar_order_screen(curve: &CurveParams) -> Result<ScalarOrderScreen, String> {
    let factors = n_minus_one_factors();
    if factor_product(&factors) != &curve.n - 1u8 {
        return Err("n-1 factorization product mismatch".into());
    }
    verify_u64_primes(&[
        "2",
        "3",
        "71",
        "131",
        "373",
        "3407",
        "17449",
        "38189",
        "187019741",
        "622491383",
        "1002328039319",
        "1297",
        "16879",
        "208150935158385979",
    ])?;
    let large_factor_certificate = n_large_factor_certificate()?;
    let n_prime_certificate = certificate(&curve.n, &factors)?;
    let primitive_root = parse(&n_prime_certificate.witness);
    let prime_divisors = unique_primes(&factors);
    if !exact_order(
        &primitive_root,
        &(&curve.n - 1u8),
        &prime_divisors,
        &curve.n,
    ) {
        return Err("primitive root exact-order replay failed".into());
    }

    let target = BigUint::from(TARGET_COLUMNS);
    let mut all_divisors = divisors(&factors);
    all_divisors.sort_by(|left, right| {
        let left_width = folded_width(left);
        let right_width = folded_width(right);
        absolute_difference(&left_width, &target)
            .cmp(&absolute_difference(&right_width, &target))
            .then_with(|| left_width.cmp(&right_width))
            .then_with(|| left.cmp(right))
    });
    let closest_widths = all_divisors
        .iter()
        .take(16)
        .map(|order| {
            let width = folded_width(order);
            WidthCandidate {
                order: order.to_string(),
                negation_folded_columns: width.to_string(),
                absolute_distance_from_target: absolute_difference(&width, &target).to_string(),
            }
        })
        .collect::<Vec<_>>();
    let selected_order = all_divisors
        .first()
        .ok_or("n-1 divisor enumeration was empty")?
        .clone();
    let selected_order_u64 = selected_order
        .to_u64()
        .ok_or("closest scalar order does not fit u64")?;
    let selected_width = folded_width(&selected_order)
        .to_u64()
        .ok_or("closest scalar width does not fit u64")?;
    let subgroup_root = primitive_root.modpow(&((&curve.n - 1u8) / &selected_order), &curve.n);
    let (selected_scalar, generators) =
        select_scalar(&subgroup_root, selected_order_u64, &curve.n)?;
    let selected_primes = factors
        .iter()
        .filter(|entry| (&selected_order % &entry.prime).is_zero())
        .map(|entry| entry.prime.clone())
        .collect::<Vec<_>>();
    if !exact_order(
        &selected_scalar,
        &selected_order,
        &selected_primes,
        &curve.n,
    ) {
        return Err("selected scalar exact-order replay failed".into());
    }
    let bits = selected_scalar.bits();
    let weight = bit_weight(&selected_scalar);
    Ok(ScalarOrderScreen {
        n_minus_one_factorization: records(&factors),
        factorization_product_verified: true,
        large_factor_certificate,
        n_prime_certificate,
        least_primitive_root: primitive_root.to_string(),
        primitive_root_verified: true,
        divisors_enumerated: all_divisors.len() as u64,
        enumeration_exhaustive: true,
        target_columns: TARGET_COLUMNS,
        closest_widths,
        selected_order: selected_order.to_string(),
        selected_negation_folded_columns: selected_width,
        selected_width_delta: selected_width.abs_diff(TARGET_COLUMNS),
        exact_order_generators_screened: generators,
        selected_scalar: selected_scalar.to_string(),
        selected_scalar_signed_binary_bits: bits,
        selected_scalar_signed_binary_weight: weight,
        selected_scalar_binary_operations: bits + weight,
        selected_scalar_exact_order_verified: true,
    })
}

fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn file_sha256(path: &Path) -> Result<String, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    Ok(hex::encode(sha256(&bytes)))
}

fn derive_anchor(curve: &CurveParams) -> Result<(Projective, BigUint, BigUint, u64), String> {
    for counter in 0u64.. {
        let digest =
            sha256(format!("{CURVE_SLUG}/scalar-orbit-round25/anchor/{counter}").as_bytes());
        let x = BigUint::from_bytes_be(&digest) % &curve.p;
        let rhs = ((&x * &x % &curve.p) * &x + &curve.a * &x + &curve.b) % &curve.p;
        let Some(y) = sqrt_mod_p(&rhs, &curve.p) else {
            continue;
        };
        if y.is_zero() {
            continue;
        }
        let negative = &curve.p - &y;
        let low_y = y.min(negative);
        let point = Projective::from_affine(&U256::from_biguint(&x), &U256::from_biguint(&low_y));
        if !bool::from(point.is_identity()) {
            return Ok((point, x, low_y, counter));
        }
    }
    unreachable!()
}

fn orbit_coefficients(scalar: &BigUint, columns: u64, n: &BigUint) -> Vec<BigUint> {
    let mut coefficients = Vec::with_capacity(columns as usize);
    let mut current = BigUint::one();
    for _ in 0..columns {
        coefficients.push(current.clone());
        current = current * scalar % n;
    }
    coefficients
}

fn fixed_base_points(
    anchor: &Projective,
    coefficients: &[BigUint],
    operations: &mut GroupOperations,
) -> Vec<Projective> {
    let windows = 256usize.div_ceil(FIXED_BASE_WINDOW_BITS);
    let entries = 1usize << FIXED_BASE_WINDOW_BITS;
    let mut tables = Vec::with_capacity(windows);
    let mut base = *anchor;
    for _ in 0..windows {
        let mut table = Vec::with_capacity(entries);
        table.push(Projective::IDENTITY);
        for digit in 1..entries {
            let next = table[digit - 1].add(&base);
            operations.fixed_base_table_additions += 1;
            table.push(next);
        }
        tables.push(table);
        for _ in 0..FIXED_BASE_WINDOW_BITS {
            base = base.double();
            operations.fixed_base_table_doublings += 1;
        }
    }
    coefficients
        .iter()
        .map(|coefficient| {
            let bytes = coefficient.to_bytes_le();
            let mut point = Projective::IDENTITY;
            for (window, byte) in bytes.iter().enumerate() {
                if *byte != 0 {
                    point = point.add(&tables[window][usize::from(*byte)]);
                    operations.fixed_base_output_additions += 1;
                }
            }
            point
        })
        .collect()
}

fn binary_scalar_mul(
    point: &Projective,
    scalar: &BigUint,
    additions: &mut u64,
    doublings: &mut u64,
) -> Projective {
    let mut result = Projective::IDENTITY;
    let mut addend = *point;
    for bit in 0..scalar.bits() {
        if scalar.bit(bit) {
            result = result.add(&addend);
            *additions += 1;
        }
        addend = addend.double();
        *doublings += 1;
    }
    result
}

fn normalize_points(
    points: &[Projective],
    operations: &mut GroupOperations,
) -> Result<Vec<(BigUint, BigUint)>, String> {
    let mut affine = Vec::with_capacity(points.len());
    for batch in points.chunks(8_192) {
        if batch.iter().any(|point| bool::from(point.is_identity())) {
            return Err("scalar orbit contains identity".into());
        }
        let mut prefixes = Vec::with_capacity(batch.len());
        let mut product = Fe::ONE;
        for point in batch {
            prefixes.push(product);
            product = product.mul(&point.z);
            operations.normalization_field_multiplications += 1;
        }
        let mut inverse = product.inv();
        operations.normalization_inversions += 1;
        let mut inverses = vec![Fe::ZERO; batch.len()];
        for index in (0..batch.len()).rev() {
            inverses[index] = inverse.mul(&prefixes[index]);
            operations.normalization_field_multiplications += 1;
            if index != 0 {
                inverse = inverse.mul(&batch[index].z);
                operations.normalization_field_multiplications += 1;
            }
        }
        for (point, inverse) in batch.iter().zip(inverses) {
            let x = point.x.mul(&inverse).to_canonical().to_biguint();
            let y = point.y.mul(&inverse).to_canonical().to_biguint();
            operations.normalization_field_multiplications += 2;
            affine.push((x, y));
        }
    }
    Ok(affine)
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

fn curve_facts(curve: &CurveParams) -> WideCurveFacts {
    WideCurveFacts {
        slug: CURVE_SLUG.into(),
        icv1: CURVE_ICV1.into(),
        family: "prime".into(),
        construction: "CurveParams::p256()".into(),
        field_degree: None,
        field_bits: curve.p.bits() as u32,
        group_order: curve.n.to_string(),
        r: curve.n.to_string(),
        cofactor: curve.h.to_string(),
        generator: [lower_hex(&curve.gx), lower_hex(&curve.gy)],
        automorphisms_available: 2,
        registered: true,
        ec1: CURVE_EC1.into(),
        curve_uid: CURVE_UID.into(),
    }
}

fn build_factor_base(
    curve: &CurveParams,
    order_screen: &ScalarOrderScreen,
    dump_path: Option<&Path>,
) -> Result<NativeFactorBase, String> {
    let scalar = parse(&order_screen.selected_scalar);
    let columns = order_screen.selected_negation_folded_columns;
    let (anchor, anchor_x, anchor_y, anchor_counter) = derive_anchor(curve)?;
    let coefficients = orbit_coefficients(&scalar, columns, &curve.n);
    let expected_wrap = (&coefficients[columns as usize - 1] * &scalar) % &curve.n;
    if expected_wrap != &curve.n - 1u8 {
        return Err("folded scalar orbit did not wrap through -1".into());
    }
    let mut operations = GroupOperations::default();
    let points = fixed_base_points(&anchor, &coefficients, &mut operations);
    let mut replay_failures = 0u64;
    for index in 0..points.len() {
        let image = binary_scalar_mul(
            &points[index],
            &scalar,
            &mut operations.transport_replay_additions,
            &mut operations.transport_replay_doublings,
        );
        let expected = if index + 1 == points.len() {
            points[0].neg()
        } else {
            points[index + 1]
        };
        if !bool::from(image.add(&expected.neg()).is_identity()) {
            replay_failures += 1;
        }
    }
    if replay_failures != 0 {
        return Err(format!(
            "{replay_failures} scalar transport replay failures"
        ));
    }
    let affine = normalize_points(&points, &mut operations)?;
    let mut orbit_points = affine
        .into_iter()
        .enumerate()
        .map(|(index, (x, y))| {
            let negative = &curve.p - &y;
            OrbitPoint {
                orbit_index: index as u64,
                point: points[index],
                x,
                low_y: y.min(negative),
            }
        })
        .collect::<Vec<_>>();
    let all_points_on_curve = orbit_points.iter().all(|row| {
        (&row.low_y * &row.low_y) % &curve.p
            == ((&row.x * &row.x % &curve.p) * &row.x + &curve.a * &row.x + &curve.b) % &curve.p
    });
    if !all_points_on_curve {
        return Err("an orbit point failed the affine curve equation".into());
    }
    orbit_points.sort_by(|left, right| {
        left.x
            .cmp(&right.x)
            .then_with(|| left.orbit_index.cmp(&right.orbit_index))
    });
    let unique_up_to_sign = orbit_points.windows(2).all(|pair| pair[0].x < pair[1].x);
    if !unique_up_to_sign {
        return Err("scalar orbit contains duplicate negation-folded columns".into());
    }

    let mut point_blob = Vec::with_capacity(orbit_points.len() * 2 * POINT_KEY_BYTES);
    let mut dump_points = Vec::with_capacity(orbit_points.len() * 2);
    for (column, row) in orbit_points.iter().enumerate() {
        point_blob.extend(point_key_bytes(&row.x, false)?);
        point_blob.extend(point_key_bytes(&row.x, true)?);
        dump_points.push(WideFbPoint {
            x: lower_hex(&row.x),
            y: lower_hex(&row.low_y),
            col: column as u64,
            coef: "1".into(),
        });
        dump_points.push(WideFbPoint {
            x: lower_hex(&row.x),
            y: lower_hex(&(&curve.p - &row.low_y)),
            col: column as u64,
            coef: (&curve.n - 1u8).to_string(),
        });
    }
    let points_sha256 = hex::encode(sha256(&point_blob));
    let mut params = BTreeMap::new();
    params.insert(
        "anchor_derivation".into(),
        "sha256-try-and-increment-lower-root".into(),
    );
    params.insert("anchor_x".into(), lower_hex(&anchor_x));
    params.insert("anchor_y".into(), lower_hex(&anchor_y));
    params.insert("endomorphism".into(), "integer-scalar".into());
    params.insert("negation_folded".into(), "true".into());
    params.insert("point_key_encoding".into(), POINT_KEY_ENCODING.into());
    params.insert("scalar".into(), scalar.to_string());
    params.insert("scalar_order".into(), order_screen.selected_order.clone());
    let identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": CURVE_SLUG,
        "family": "scalar-orbit",
        "params": params,
        "columns": columns,
        "signed_points": 2 * columns,
        "points_sha256": points_sha256,
    });
    let (fb_id, fb_sha256) = short_id("FB1", &identity)?;
    let factor_base = WideFactorBaseFacts {
        fb_id: fb_id.clone(),
        fb_sha256: fb_sha256.clone(),
        family: "scalar-orbit".into(),
        params,
        description: format!(
            "negation-folded exact scalar orbit on {CURVE_SLUG}; one unknown anchor log class"
        ),
        signed_points: 2 * columns,
        abscissae: columns,
        columns,
        dimension: None,
        points_sha256: points_sha256.clone(),
    };
    let build_adds = operations.fixed_base_table_additions + operations.fixed_base_output_additions;
    let dump = WideFactorBaseDump {
        schema: DUMP_SCHEMA.into(),
        curve: curve_facts(curve),
        factor_base,
        build_adds,
        build_doubles: operations.fixed_base_table_doublings,
        points: dump_points,
    };
    let mut dump_bytes = serde_json::to_vec_pretty(&dump).map_err(|error| error.to_string())?;
    dump_bytes.push(b'\n');
    let dump_sha256 = hex::encode(sha256(&dump_bytes));
    let roundtrip: WideFactorBaseDump =
        serde_json::from_slice(&dump_bytes).map_err(|error| error.to_string())?;
    if roundtrip != dump {
        return Err("wide factor-base dump failed an exact serialization round trip".into());
    }
    if let Some(path) = dump_path {
        fs::write(path, &dump_bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    }
    let rebuilt_identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": dump.curve.slug,
        "family": dump.factor_base.family,
        "params": dump.factor_base.params,
        "columns": dump.factor_base.columns,
        "signed_points": dump.factor_base.signed_points,
        "points_sha256": dump.factor_base.points_sha256,
    });
    let (rebuilt_id, rebuilt_sha) = short_id("FB1", &rebuilt_identity)?;
    if rebuilt_id != fb_id || rebuilt_sha != fb_sha256 {
        return Err("wide dump FB1 identity did not rebuild".into());
    }
    let all_nonidentity = orbit_points
        .iter()
        .all(|row| !bool::from(row.point.is_identity()));
    Ok(NativeFactorBase {
        fb_id,
        fb_sha256,
        family: "scalar-orbit".into(),
        columns,
        signed_points: 2 * columns,
        points_sha256,
        point_key_encoding: POINT_KEY_ENCODING.into(),
        anchor_x: lower_hex(&anchor_x),
        anchor_y: lower_hex(&anchor_y),
        anchor_derivation_counter: anchor_counter,
        anchor_construction_log_recorded: false,
        all_points_nonidentity: all_nonidentity,
        all_points_on_curve,
        unique_up_to_sign,
        adjacent_transport_relations: columns,
        wrap_relation_uses_negation: true,
        transport_replay_failures: replay_failures,
        all_reported_transport_relations_replayed: replay_failures == 0,
        independent_log_quotient_dimension: 1,
        ordinary_zero_relations_supply_anchor_information: false,
        dump_schema: DUMP_SCHEMA.into(),
        dump_bytes: dump_bytes.len() as u64,
        dump_sha256,
        rebuild_identity_verified: true,
        operations,
    })
}

fn gamma_ratio(k: u64) -> f64 {
    let mut ratio = std::f64::consts::PI.sqrt() / 2.0;
    for value in 1..k {
        ratio *= (value as f64 + 0.5) / value as f64;
    }
    ratio
}

fn disjoint_probability(columns: u64) -> f64 {
    (0..9).fold(1.0, |probability, offset| {
        probability * (columns - 8 - offset) as f64 / (columns - offset) as f64
    })
}

fn projection(case: &str, k: u64, columns: u64, n: &BigUint, class: &str) -> Projection {
    let probability = disjoint_probability(columns);
    let n_f64 = n.to_f64().expect("P-256 order converts to f64");
    let samples = (2.0 * n_f64 / probability).sqrt() * gamma_ratio(k);
    let oracle_s = samples / n_f64.sqrt();
    let direct = samples * 8.0;
    let one_list_bytes = samples * 0.5 * ENTRY_BYTES as f64;
    Projection {
        case: case.into(),
        required_target_collisions: k,
        disjoint_probability: probability,
        oracle_s,
        oracle_ratio_to_rho: oracle_s / RHO_S,
        log2_oracle_samples: samples.log2(),
        direct_additions_per_sample: 8.0,
        log2_direct_group_additions: direct.log2(),
        direct_ratio_to_rho: direct / (RHO_S * n_f64.sqrt()),
        log2_one_list_bytes: one_list_bytes.log2(),
        below_rho: direct < RHO_S * n_f64.sqrt(),
        classification: class.into(),
    }
}

fn run(cli: &Cli) -> Result<ResultFile, String> {
    let dependency_sha256 = file_sha256(&cli.round24)?;
    if dependency_sha256 != ROUND24_SHA256 {
        return Err(format!(
            "round-24 dependency hash mismatch: expected {ROUND24_SHA256}, got {dependency_sha256}"
        ));
    }
    let curve = CurveParams::p256();
    let cm = cm_screen(&curve)?;
    let order_screen = scalar_order_screen(&curve)?;
    let factor_base = build_factor_base(&curve, &order_screen, cli.dump.as_deref())?;
    let scalar = parse(&order_screen.selected_scalar);
    let multiplication_degree = &scalar * &scalar;
    let projections = vec![
        projection(
            "known anchor H=G",
            1,
            factor_base.columns,
            &curve.n,
            "The ideal K=1 constant is rho parity only; with known scalar labels the target representation is a generic claw/rho search, and direct sample construction is above rho.",
        ),
        projection(
            "unknown hash-derived anchor",
            2,
            factor_base.columns,
            &curve.n,
            "At least two independent target equations are required to eliminate the unknown anchor log; even the free-oracle boundary is above rho.",
        ),
    ];
    let degree_screen = DegreeScreen {
        exact_support_polynomial_degree: factor_base.columns,
        multiplication_map_degree: multiplication_degree.to_string(),
        multiplication_map_degree_log2: multiplication_degree
            .to_f64()
            .ok_or("multiplication-map degree did not convert to f64")?
            .log2(),
        short_scalar_evaluation_chain_is_low_degree_membership: false,
        structured_residual_degree_of_regularity: None,
        degree_gate_passed: false,
        classification: "The orbit has an exact one-class log quotient, but coordinate membership is a support polynomial of degree B (or a rational multiplication map of degree m^2). No structured residual degree <=5 was measured.".into(),
    };
    let false_positives = 0;
    let false_negatives = 0;
    let promotion_gates = PromotionGates {
        zero_false_positives: true,
        zero_false_negatives: true,
        exact_transport_replay: factor_base.all_reported_transport_relations_replayed,
        structured_residual_degree_at_most_5: false,
        relation_collection_below_2_120: false,
        per_usable_row_below_2_103: false,
        materialized_storage_below_2_50: false,
        complete_cost_below_rho: false,
        selector_and_anchor_non_generic: false,
        promoted: false,
    };
    Ok(ResultFile {
        schema: "p256.scalar_orbit_factor_base_screen/v1".into(),
        curve: CURVE_SLUG.into(),
        relation_model: "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation".into(),
        dependency_sha256,
        external_factor_discovery_trusted: false,
        cm_screen: cm,
        scalar_order_screen: order_screen,
        factor_base,
        projections,
        degree_screen,
        false_positives,
        false_negatives,
        promotion_gates,
        full_depth_unplanted_attempted: false,
        decision: "An exact native P-256 scalar-orbit factor base reaches one unknown log class, but it does not improve index calculus end to end. Every CM endomorphism restricts to a known scalar on E(F_p); an unknown anchor needs two target collisions, a known anchor is generic rho/claw search, direct generation and storage fail, and no low-degree coordinate membership or residual regularity is obtained.".into(),
    })
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(result) => {
            let mut bytes = match serde_json::to_vec_pretty(&result) {
                Ok(bytes) => bytes,
                Err(error) => {
                    eprintln!("serialization failed: {error}");
                    return ExitCode::FAILURE;
                }
            };
            bytes.push(b'\n');
            if let Some(path) = cli.out {
                if let Err(error) = fs::write(&path, &bytes) {
                    eprintln!("{}: {error}", path.display());
                    return ExitCode::FAILURE;
                }
            } else {
                print!("{}", String::from_utf8_lossy(&bytes));
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("P-256 scalar-orbit screen failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cm_factorization_is_exact_and_squarefree() {
        let curve = CurveParams::p256();
        let report = cm_screen(&curve).unwrap();
        assert!(report.factorization_product_verified);
        assert!(report.all_factors_distinct_prime);
        assert!(report.fundamental_discriminant);
        assert_eq!(report.frobenius_conductor, "1");
        assert!(report.minimum_noninteger_endomorphism_degree_log2 > 255.9);
    }

    #[test]
    fn scalar_width_enumeration_selects_expected_orbit() {
        let report = scalar_order_screen(&CurveParams::p256()).unwrap();
        assert_eq!(report.divisors_enumerated, 10_240);
        assert_eq!(report.selected_order, "279184");
        assert_eq!(report.selected_negation_folded_columns, 139_592);
        assert_eq!(report.exact_order_generators_screened, 139_584);
        assert!(report.selected_scalar_exact_order_verified);
    }

    #[test]
    fn lucas_certificate_rejects_inexact_factorization() {
        let factors = vec![factor("2", 1), factor("3", 1)];
        assert!(lucas_witness(&BigUint::from(11u8), &factors).is_err());
    }

    #[test]
    fn unknown_anchor_needs_more_than_rho_oracle_constant() {
        let curve = CurveParams::p256();
        let row = projection("unknown", 2, 139_592, &curve.n, "test");
        assert!(row.oracle_ratio_to_rho > 1.4);
        assert!(!row.below_rho);
    }

    #[test]
    fn point_key_orders_both_signs() {
        let x = BigUint::from(7u8);
        assert!(point_key_bytes(&x, false).unwrap() < point_key_bytes(&x, true).unwrap());
    }
}
