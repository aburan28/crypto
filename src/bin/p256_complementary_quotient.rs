//! Round 302: exact complementary elliptic quotient of the P-256 quadratic cover.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::CurveParams;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;
use serde_json::{json, Value};

const REGISTRY_SHA256: &str = "b1b8439baf4d18307bbc81626f28fb840bd48ba69eda585fde4a9090fcc10fdf";
const COVERS_SHA256: &str = "39811d799f2b12982f2b9e4cd0a4408c21fb1b533f117181e98f4fe506afeef5";
const ROUND300_SHA256: &str = "7d65f3e015c6d6c6d24b64ea2e59a6e9e8653c88a1900835838b1bc0a1ec36e5";
const ROUND301_SHA256: &str = "ae4b44315181c18a40a790ae1bc784ffdb05e3061d442b8f1be78a1f25dc20e9";
const ISOGENOUS_SLUG: &str = "icv1-fp256-t89188191154553853111372247798585809583-8ae9296a";
const P256_EC1: &str = "EC1P256Cp256h0523b774e066";
const ISOGENOUS_EC1: &str = "EC1P256Cfph72dde3c076e4";
const FACTOR_BASE_ID: &str = "FB1hc72514a2a8d3";
const REGISTERED_S17_FACTOR_BASE_ID: &str = "FB1h2f8621cda105";
const MINIMUM_WIDTH_RATIO_TO_RHO: f64 = 13.920_747_397_073_491;
const REGISTERED_S17_RATIO_TO_RHO: f64 = 394.425_280;
const P256_MAP_SAMPLES: u64 = 64;

#[derive(Parser)]
#[command(about = "Certify the complementary elliptic quotient of the P-256 split cover")]
struct Cli {
    #[arg(long)]
    registry: PathBuf,
    #[arg(long)]
    covers: PathBuf,
    #[arg(long)]
    round300: PathBuf,
    #[arg(long)]
    round301: PathBuf,
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
struct OperationCounts {
    field_additions: u64,
    field_subtractions: u64,
    field_multiplications: u64,
    field_squarings: u64,
    field_inversions: u64,
    modular_exponentiations: u64,
    modular_exponentiation_multiplications: u64,
    modular_exponentiation_squarings: u64,
    legendre_tests: u64,
    square_root_exponentiations: u64,
    group_additions: u64,
    group_doublings: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum Point {
    Infinity,
    Affine { x: BigUint, y: BigUint },
}

#[derive(Clone, Serialize)]
struct PointReceipt {
    infinity: bool,
    x: Option<String>,
    y: Option<String>,
}

#[derive(Clone, Serialize)]
struct MapSample {
    u: String,
    v: String,
    primary_x: String,
    primary_y: String,
    quartic_x: String,
    quartic_y: String,
    complementary_u: String,
    complementary_v: String,
}

#[derive(Serialize)]
struct FibreSummary {
    distinct_images: u64,
    minimum_affine_fibre: u64,
    maximum_affine_fibre: u64,
    histogram: BTreeMap<u64, u64>,
}

#[derive(Serialize)]
struct ToyControl {
    p: u64,
    a: u64,
    b: u64,
    selection_rule: String,
    nonsingular: bool,
    elliptic_primary_order: u64,
    elliptic_complementary_order: u64,
    affine_cover_points: u64,
    primary_map_replays: u64,
    quartic_map_replays: u64,
    complementary_map_replays: u64,
    complementary_exceptional_points: u64,
    replay_failures: u64,
    primary_fibres: FibreSummary,
    complementary_fibres: FibreSummary,
    operations: OperationCounts,
}

#[derive(Serialize)]
struct P256MapControl {
    requested_points: u64,
    scanned_u_candidates: u64,
    accepted_cover_points: u64,
    primary_map_replays: u64,
    quartic_map_replays: u64,
    complementary_map_replays: u64,
    exceptional_points: u64,
    replay_failures: u64,
    sample_digest_encoding: String,
    sample_sha256: String,
    first_sample: MapSample,
    last_sample: MapSample,
    operations: OperationCounts,
}

#[derive(Serialize)]
struct ComplementaryCurve {
    model: String,
    p: String,
    a2: String,
    a4: String,
    a6: String,
    short_a: String,
    short_b: String,
    nonsingular: bool,
}

#[derive(Serialize)]
struct SubgroupCertificate {
    witness_selection: String,
    scanned_u_candidates: u64,
    witness: PointReceipt,
    witness_on_curve: bool,
    scalar: String,
    primary_n_times_witness: PointReceipt,
    replay_n_times_witness: PointReceipt,
    scalar_algorithms_match: bool,
    n_times_witness_is_infinity: bool,
    hasse_upper_bound_used: String,
    hasse_upper_bound_below_two_n: bool,
    only_positive_multiple_of_n_in_hasse_interval_is_n: bool,
    n_divides_complementary_group_order: bool,
    rational_order_n_transport_exists: bool,
    search_operations: OperationCounts,
    primary_scalar_operations: OperationCounts,
    replay_scalar_operations: OperationCounts,
}

#[derive(Serialize)]
struct BoundaryRow {
    variant: String,
    independent_rows_per_cover_event: Option<u64>,
    optimistic_ratio_to_rho: Option<f64>,
    complete_cost_measured: bool,
    status: String,
}

#[derive(Serialize)]
struct Gates {
    dependency_hashes_and_schemas_checked: bool,
    zero_map_and_replay_failures: bool,
    two_distinct_order_n_rows_per_cover_event: bool,
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
    dependencies: Vec<Dependency>,
    exact_maps: Vec<String>,
    complementary_curve: ComplementaryCurve,
    toy_control: ToyControl,
    p256_map_control: P256MapControl,
    subgroup_certificate: SubgroupCertificate,
    boundary_table: Vec<BoundaryRow>,
    distinct_p256_order_n_rows_per_cover_event: u64,
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

fn addm(left: &BigUint, right: &BigUint, p: &BigUint, ops: &mut OperationCounts) -> BigUint {
    ops.field_additions += 1;
    (left + right) % p
}

fn subm(left: &BigUint, right: &BigUint, p: &BigUint, ops: &mut OperationCounts) -> BigUint {
    ops.field_subtractions += 1;
    if left >= right {
        left - right
    } else {
        left + p - right
    }
}

fn mulm(left: &BigUint, right: &BigUint, p: &BigUint, ops: &mut OperationCounts) -> BigUint {
    ops.field_multiplications += 1;
    (left * right) % p
}

fn sqrm(value: &BigUint, p: &BigUint, ops: &mut OperationCounts) -> BigUint {
    ops.field_squarings += 1;
    (value * value) % p
}

fn powm(base: &BigUint, exponent: &BigUint, p: &BigUint, ops: &mut OperationCounts) -> BigUint {
    ops.modular_exponentiations += 1;
    let mut value = exponent.clone();
    let mut power = base % p;
    let mut out = BigUint::one();
    while !value.is_zero() {
        if value.bit(0) {
            out = (&out * &power) % p;
            ops.modular_exponentiation_multiplications += 1;
        }
        power = (&power * &power) % p;
        ops.modular_exponentiation_squarings += 1;
        value >>= 1usize;
    }
    out
}

fn invm(value: &BigUint, p: &BigUint, ops: &mut OperationCounts) -> Result<BigUint, String> {
    if value.is_zero() {
        return Err("attempted inversion of zero".into());
    }
    ops.field_inversions += 1;
    Ok(powm(value, &(p - BigUint::from(2u8)), p, ops))
}

fn legendre(value: &BigUint, p: &BigUint, ops: &mut OperationCounts) -> BigUint {
    ops.legendre_tests += 1;
    powm(value, &((p - BigUint::one()) >> 1usize), p, ops)
}

fn canonical_sqrt(
    value: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> Result<Option<BigUint>, String> {
    if value.is_zero() {
        return Ok(Some(BigUint::zero()));
    }
    if p % BigUint::from(4u8) != BigUint::from(3u8) {
        return Err("canonical_sqrt requires p = 3 mod 4".into());
    }
    if legendre(value, p, ops) != BigUint::one() {
        return Ok(None);
    }
    ops.square_root_exponentiations += 1;
    let root = powm(value, &((p + BigUint::one()) >> 2usize), p, ops);
    if sqrm(&root, p, ops) != value % p {
        return Err("square-root replay failed".into());
    }
    let negative = p - &root;
    Ok(Some(root.min(negative)))
}

fn primary_rhs(
    x: &BigUint,
    a: &BigUint,
    b: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> BigUint {
    let x2 = sqrm(x, p, ops);
    let x3 = mulm(&x2, x, p, ops);
    let ax = mulm(a, x, p, ops);
    addm(&addm(&x3, &ax, p, ops), b, p, ops)
}

fn cover_rhs(
    u: &BigUint,
    a: &BigUint,
    b: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> BigUint {
    let u2 = sqrm(u, p, ops);
    let u4 = sqrm(&u2, p, ops);
    let u6 = mulm(&u4, &u2, p, ops);
    let au2 = mulm(a, &u2, p, ops);
    addm(&addm(&u6, &au2, p, ops), b, p, ops)
}

fn quartic_rhs(
    x: &BigUint,
    a: &BigUint,
    b: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> BigUint {
    let x2 = sqrm(x, p, ops);
    let x4 = sqrm(&x2, p, ops);
    let ax2 = mulm(a, &x2, p, ops);
    let bx = mulm(b, x, p, ops);
    addm(&addm(&x4, &ax2, p, ops), &bx, p, ops)
}

fn complementary_rhs(
    u: &BigUint,
    a2: &BigUint,
    b2: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> BigUint {
    let u2 = sqrm(u, p, ops);
    let u3 = mulm(&u2, u, p, ops);
    let a2u2 = mulm(a2, &u2, p, ops);
    addm(&addm(&u3, &a2u2, p, ops), b2, p, ops)
}

fn replay_maps(
    u: &BigUint,
    v: &BigUint,
    a: &BigUint,
    b: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> Result<Option<MapSample>, String> {
    let v2 = sqrm(v, p, ops);
    if v2 != cover_rhs(u, a, b, p, ops) {
        return Err("point is not on the cover".into());
    }
    let u2 = sqrm(u, p, ops);
    if v2 != primary_rhs(&u2, a, b, p, ops) {
        return Err("primary quotient equation failed".into());
    }
    let quartic_y = mulm(u, v, p, ops);
    if sqrm(&quartic_y, p, ops) != quartic_rhs(&u2, a, b, p, ops) {
        return Err("quartic quotient equation failed".into());
    }
    if u.is_zero() {
        return Ok(None);
    }
    let inv_u2 = invm(&u2, p, ops)?;
    let complement_u = mulm(b, &inv_u2, p, ops);
    let u3 = mulm(&u2, u, p, ops);
    let inv_u3 = invm(&u3, p, ops)?;
    let bv = mulm(b, v, p, ops);
    let complement_v = mulm(&bv, &inv_u3, p, ops);
    let b2 = sqrm(b, p, ops);
    if sqrm(&complement_v, p, ops) != complementary_rhs(&complement_u, a, &b2, p, ops) {
        return Err("complementary quotient equation failed".into());
    }
    Ok(Some(MapSample {
        u: u.to_string(),
        v: v.to_string(),
        primary_x: u2.to_string(),
        primary_y: v.to_string(),
        quartic_x: u2.to_string(),
        quartic_y: quartic_y.to_string(),
        complementary_u: complement_u.to_string(),
        complementary_v: complement_v.to_string(),
    }))
}

fn point_receipt(point: &Point) -> PointReceipt {
    match point {
        Point::Infinity => PointReceipt {
            infinity: true,
            x: None,
            y: None,
        },
        Point::Affine { x, y } => PointReceipt {
            infinity: false,
            x: Some(x.to_string()),
            y: Some(y.to_string()),
        },
    }
}

fn point_on_complement(
    point: &Point,
    a2: &BigUint,
    b2: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> bool {
    match point {
        Point::Infinity => true,
        Point::Affine { x, y } => sqrm(y, p, ops) == complementary_rhs(x, a2, b2, p, ops),
    }
}

fn point_double(
    point: &Point,
    a2: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> Result<Point, String> {
    let Point::Affine { x, y } = point else {
        return Ok(Point::Infinity);
    };
    ops.group_doublings += 1;
    if y.is_zero() {
        return Ok(Point::Infinity);
    }
    let x2 = sqrm(x, p, ops);
    let three_x2 = mulm(&BigUint::from(3u8), &x2, p, ops);
    let two_a2_x = mulm(&BigUint::from(2u8), &mulm(a2, x, p, ops), p, ops);
    let numerator = addm(&three_x2, &two_a2_x, p, ops);
    let denominator = mulm(&BigUint::from(2u8), y, p, ops);
    let lambda = mulm(&numerator, &invm(&denominator, p, ops)?, p, ops);
    let lambda2 = sqrm(&lambda, p, ops);
    let two_x = addm(x, x, p, ops);
    let x3 = subm(&subm(&lambda2, a2, p, ops), &two_x, p, ops);
    let y3 = subm(&mulm(&lambda, &subm(x, &x3, p, ops), p, ops), y, p, ops);
    Ok(Point::Affine { x: x3, y: y3 })
}

fn point_add(
    left: &Point,
    right: &Point,
    a2: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> Result<Point, String> {
    match (left, right) {
        (Point::Infinity, _) => Ok(right.clone()),
        (_, Point::Infinity) => Ok(left.clone()),
        (Point::Affine { x: x1, y: y1 }, Point::Affine { x: x2, y: y2 }) => {
            ops.group_additions += 1;
            if x1 == x2 {
                if y1 == y2 {
                    return point_double(left, a2, p, ops);
                }
                return Ok(Point::Infinity);
            }
            let numerator = subm(y2, y1, p, ops);
            let denominator = subm(x2, x1, p, ops);
            let lambda = mulm(&numerator, &invm(&denominator, p, ops)?, p, ops);
            let x3 = subm(
                &subm(&subm(&sqrm(&lambda, p, ops), a2, p, ops), x1, p, ops),
                x2,
                p,
                ops,
            );
            let y3 = subm(&mulm(&lambda, &subm(x1, &x3, p, ops), p, ops), y1, p, ops);
            Ok(Point::Affine { x: x3, y: y3 })
        }
    }
}

fn scalar_mul_lsb(
    point: &Point,
    scalar: &BigUint,
    a2: &BigUint,
    b2: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> Result<Point, String> {
    let mut value = scalar.clone();
    let mut base = point.clone();
    let mut out = Point::Infinity;
    while !value.is_zero() {
        if value.bit(0) {
            out = point_add(&out, &base, a2, p, ops)?;
        }
        base = point_double(&base, a2, p, ops)?;
        value >>= 1usize;
    }
    if !point_on_complement(&out, a2, b2, p, ops) {
        return Err("LSB scalar result is off curve".into());
    }
    Ok(out)
}

fn scalar_mul_msb(
    point: &Point,
    scalar: &BigUint,
    a2: &BigUint,
    b2: &BigUint,
    p: &BigUint,
    ops: &mut OperationCounts,
) -> Result<Point, String> {
    let mut out = Point::Infinity;
    for bit in (0..scalar.bits()).rev() {
        out = point_double(&out, a2, p, ops)?;
        if scalar.bit(bit) {
            out = point_add(&out, point, a2, p, ops)?;
        }
    }
    if !point_on_complement(&out, a2, b2, p, ops) {
        return Err("MSB scalar result is off curve".into());
    }
    Ok(out)
}

fn count_curve_order(
    p: u64,
    a2: u64,
    a4: u64,
    a6: u64,
    ops: &mut OperationCounts,
) -> Result<u64, String> {
    let modulus = BigUint::from(p);
    let mut order = 1u64;
    for x in 0..p {
        let x = BigUint::from(x);
        let x2 = sqrm(&x, &modulus, ops);
        let x3 = mulm(&x2, &x, &modulus, ops);
        let rhs = addm(
            &addm(
                &x3,
                &mulm(&BigUint::from(a2), &x2, &modulus, ops),
                &modulus,
                ops,
            ),
            &addm(
                &mulm(&BigUint::from(a4), &x, &modulus, ops),
                &BigUint::from(a6),
                &modulus,
                ops,
            ),
            &modulus,
            ops,
        );
        if rhs.is_zero() {
            order += 1;
        } else if legendre(&rhs, &modulus, ops) == BigUint::one() {
            order += 2;
        }
    }
    Ok(order)
}

fn fibre_summary(images: BTreeMap<String, u64>) -> FibreSummary {
    let mut histogram = BTreeMap::new();
    let mut minimum = u64::MAX;
    let mut maximum = 0;
    for size in images.values().copied() {
        minimum = minimum.min(size);
        maximum = maximum.max(size);
        *histogram.entry(size).or_insert(0) += 1;
    }
    FibreSummary {
        distinct_images: images.len() as u64,
        minimum_affine_fibre: if images.is_empty() { 0 } else { minimum },
        maximum_affine_fibre: maximum,
        histogram,
    }
}

fn toy_control() -> Result<ToyControl, String> {
    let (p_u64, a_u64, b_u64) = (1019u64, 2u64, 3u64);
    let p = BigUint::from(p_u64);
    let a = BigUint::from(a_u64);
    let b = BigUint::from(b_u64);
    let discriminant =
        (4u128 * u128::from(a_u64).pow(3) + 27u128 * u128::from(b_u64).pow(2)) % u128::from(p_u64);
    if discriminant == 0 {
        return Err("toy curve is singular".into());
    }
    let mut ops = OperationCounts::default();
    let b2 = sqrm(&b, &p, &mut ops);
    let mut affine = 0u64;
    let mut primary_replays = 0u64;
    let mut quartic_replays = 0u64;
    let mut complement_replays = 0u64;
    let mut exceptional = 0u64;
    let mut failures = 0u64;
    let mut primary_images = BTreeMap::new();
    let mut complement_images = BTreeMap::new();
    for raw_u in 0..p_u64 {
        let u = BigUint::from(raw_u);
        let rhs = cover_rhs(&u, &a, &b, &p, &mut ops);
        let roots = if rhs.is_zero() {
            vec![BigUint::zero()]
        } else if let Some(root) = canonical_sqrt(&rhs, &p, &mut ops)? {
            let negative = &p - &root;
            if root == negative {
                vec![root]
            } else {
                vec![root, negative]
            }
        } else {
            Vec::new()
        };
        for v in roots {
            affine += 1;
            match replay_maps(&u, &v, &a, &b, &p, &mut ops) {
                Ok(sample) => {
                    primary_replays += 1;
                    quartic_replays += 1;
                    let x = (&u * &u) % &p;
                    *primary_images.entry(format!("{x}:{}", v)).or_insert(0) += 1;
                    if let Some(sample) = sample {
                        complement_replays += 1;
                        *complement_images
                            .entry(format!(
                                "{}:{}",
                                sample.complementary_u, sample.complementary_v
                            ))
                            .or_insert(0) += 1;
                    } else {
                        exceptional += 1;
                    }
                }
                Err(_) => failures += 1,
            }
        }
    }
    if failures != 0 || primary_replays != affine || complement_replays + exceptional != affine {
        return Err("complete toy cover replay failed".into());
    }
    let primary_order = count_curve_order(p_u64, 0, a_u64, b_u64, &mut ops)?;
    let complementary_order = count_curve_order(
        p_u64,
        a_u64,
        0,
        b2.to_u64().ok_or("toy b^2 does not fit u64")?,
        &mut ops,
    )?;
    Ok(ToyControl {
        p: p_u64,
        a: a_u64,
        b: b_u64,
        selection_rule: "first literal preregistered cell [(1019,2,3)]".into(),
        nonsingular: true,
        elliptic_primary_order: primary_order,
        elliptic_complementary_order: complementary_order,
        affine_cover_points: affine,
        primary_map_replays: primary_replays,
        quartic_map_replays: quartic_replays,
        complementary_map_replays: complement_replays,
        complementary_exceptional_points: exceptional,
        replay_failures: failures,
        primary_fibres: fibre_summary(primary_images),
        complementary_fibres: fibre_summary(complement_images),
        operations: ops,
    })
}

fn fixed_be(value: &BigUint, bytes: usize) -> Result<Vec<u8>, String> {
    let encoded = value.to_bytes_be();
    if encoded.len() > bytes {
        return Err("integer exceeds fixed encoding".into());
    }
    let mut out = vec![0u8; bytes - encoded.len()];
    out.extend(encoded);
    Ok(out)
}

fn p256_map_control(curve: &CurveParams) -> Result<P256MapControl, String> {
    let mut ops = OperationCounts::default();
    let mut accepted = 0u64;
    let mut scanned = 0u64;
    let mut u = BigUint::one();
    let mut samples = Vec::new();
    let mut digest_input = Vec::new();
    while accepted < P256_MAP_SAMPLES {
        scanned += 1;
        let rhs = cover_rhs(&u, &curve.a, &curve.b, &curve.p, &mut ops);
        if let Some(v) = canonical_sqrt(&rhs, &curve.p, &mut ops)? {
            let sample = replay_maps(&u, &v, &curve.a, &curve.b, &curve.p, &mut ops)?
                .ok_or("nonzero P-256 u became an exceptional quotient point")?;
            for coordinate in [
                &sample.u,
                &sample.v,
                &sample.primary_x,
                &sample.primary_y,
                &sample.quartic_x,
                &sample.quartic_y,
                &sample.complementary_u,
                &sample.complementary_v,
            ] {
                let value = BigUint::parse_bytes(coordinate.as_bytes(), 10)
                    .ok_or("sample coordinate is not decimal")?;
                digest_input.extend(fixed_be(&value, 32)?);
            }
            samples.push(sample);
            accepted += 1;
        }
        u += BigUint::one();
    }
    let first = samples.first().cloned().ok_or("no P-256 map samples")?;
    let last = samples.last().cloned().ok_or("no P-256 map samples")?;
    Ok(P256MapControl {
        requested_points: P256_MAP_SAMPLES,
        scanned_u_candidates: scanned,
        accepted_cover_points: accepted,
        primary_map_replays: accepted,
        quartic_map_replays: accepted,
        complementary_map_replays: accepted,
        exceptional_points: 0,
        replay_failures: 0,
        sample_digest_encoding: "64 records * (u,v,pi1.x,pi1.y,Q.x,Q.y,pi2.U,pi2.V) * be32".into(),
        sample_sha256: sha256_hex(&digest_input),
        first_sample: first,
        last_sample: last,
        operations: ops,
    })
}

fn complement_curve(curve: &CurveParams) -> Result<ComplementaryCurve, String> {
    let p = &curve.p;
    let mut ops = OperationCounts::default();
    let b2 = sqrm(&curve.b, p, &mut ops);
    let inv3 = invm(&BigUint::from(3u8), p, &mut ops)?;
    let inv27 = invm(&BigUint::from(27u8), p, &mut ops)?;
    let a_sq = sqrm(&curve.a, p, &mut ops);
    let short_a = subm(
        &BigUint::zero(),
        &mulm(&a_sq, &inv3, p, &mut ops),
        p,
        &mut ops,
    );
    let a_cubed = mulm(&a_sq, &curve.a, p, &mut ops);
    let two_a3_over_27 = mulm(
        &mulm(&BigUint::from(2u8), &a_cubed, p, &mut ops),
        &inv27,
        p,
        &mut ops,
    );
    let short_b = addm(&b2, &two_a3_over_27, p, &mut ops);
    let four_a3 = mulm(
        &BigUint::from(4u8),
        &mulm(&sqrm(&short_a, p, &mut ops), &short_a, p, &mut ops),
        p,
        &mut ops,
    );
    let twenty_seven_b2 = mulm(
        &BigUint::from(27u8),
        &sqrm(&short_b, p, &mut ops),
        p,
        &mut ops,
    );
    let nonsingular = !addm(&four_a3, &twenty_seven_b2, p, &mut ops).is_zero();
    if !nonsingular {
        return Err("complementary P-256 quotient is singular".into());
    }
    Ok(ComplementaryCurve {
        model: "V^2=U^3+a*U^2+b^2; short x=U+a/3".into(),
        p: p.to_string(),
        a2: curve.a.to_string(),
        a4: "0".into(),
        a6: b2.to_string(),
        short_a: short_a.to_string(),
        short_b: short_b.to_string(),
        nonsingular,
    })
}

fn subgroup_certificate(curve: &CurveParams) -> Result<SubgroupCertificate, String> {
    let p = &curve.p;
    let a2 = &curve.a;
    let mut search_ops = OperationCounts::default();
    let b2 = sqrm(&curve.b, p, &mut search_ops);
    let mut scanned = 0u64;
    let mut raw_u = BigUint::zero();
    let witness = loop {
        scanned += 1;
        let rhs = complementary_rhs(&raw_u, a2, &b2, p, &mut search_ops);
        if !rhs.is_zero() {
            if let Some(v) = canonical_sqrt(&rhs, p, &mut search_ops)? {
                if !v.is_zero() {
                    break Point::Affine {
                        x: raw_u.clone(),
                        y: v,
                    };
                }
            }
        }
        raw_u += BigUint::one();
    };
    let witness_on_curve = point_on_complement(&witness, a2, &b2, p, &mut search_ops);
    if !witness_on_curve {
        return Err("complementary witness is off curve".into());
    }
    let mut primary_ops = OperationCounts::default();
    let primary = scalar_mul_lsb(&witness, &curve.n, a2, &b2, p, &mut primary_ops)?;
    let mut replay_ops = OperationCounts::default();
    let replay = scalar_mul_msb(&witness, &curve.n, a2, &b2, p, &mut replay_ops)?;
    if primary != replay {
        return Err("independent scalar algorithms disagree".into());
    }
    let n_result_is_infinity = primary == Point::Infinity;
    let hasse_upper = p + BigUint::one() + (BigUint::one() << 129usize);
    let two_n = &curve.n << 1usize;
    let hasse_below_two_n = hasse_upper < two_n;
    if !hasse_below_two_n {
        return Err("coarse Hasse upper bound is not below 2n".into());
    }
    if n_result_is_infinity {
        return Err("preregistered witness did not refute complementary n-torsion".into());
    }
    Ok(SubgroupCertificate {
        witness_selection: "first U>=0 with nonzero quadratic-residue RHS and canonical nonzero V"
            .into(),
        scanned_u_candidates: scanned,
        witness: point_receipt(&witness),
        witness_on_curve,
        scalar: curve.n.to_string(),
        primary_n_times_witness: point_receipt(&primary),
        replay_n_times_witness: point_receipt(&replay),
        scalar_algorithms_match: true,
        n_times_witness_is_infinity: false,
        hasse_upper_bound_used: hasse_upper.to_string(),
        hasse_upper_bound_below_two_n: hasse_below_two_n,
        only_positive_multiple_of_n_in_hasse_interval_is_n: true,
        n_divides_complementary_group_order: false,
        rational_order_n_transport_exists: false,
        search_operations: search_ops,
        primary_scalar_operations: primary_ops,
        replay_scalar_operations: replay_ops,
    })
}

fn dependency_checks(cli: &Cli) -> Result<(Vec<Dependency>, Value, Value, Value, Value), String> {
    let (registry_dep, registry) =
        checked_json(&cli.registry, REGISTRY_SHA256, "ICV1 registry/v1")?;
    let (covers_dep, covers) = checked_json(&cli.covers, COVERS_SHA256, "curve-covers/v1")?;
    let (round300_dep, round300) =
        checked_json(&cli.round300, ROUND300_SHA256, "p256-cover-fiber-screen/v1")?;
    let (round301_dep, round301) = checked_json(
        &cli.round301,
        ROUND301_SHA256,
        "p256-algebraic-escape-census/v1",
    )?;
    if registry.pointer("/schema_version").and_then(Value::as_u64) != Some(1) {
        return Err("registry schema mismatch".into());
    }
    let curves = registry
        .pointer("/curves")
        .and_then(Value::as_array)
        .ok_or("registry curves missing")?;
    for (slug, ec1) in [(CURVE_SLUG, P256_EC1), (ISOGENOUS_SLUG, ISOGENOUS_EC1)] {
        let row = curves
            .iter()
            .find(|row| row.pointer("/slug").and_then(Value::as_str) == Some(slug))
            .ok_or_else(|| format!("registry lacks {slug}"))?;
        let found = row
            .pointer("/representations")
            .and_then(Value::as_array)
            .is_some_and(|rows| {
                rows.iter()
                    .any(|item| item.pointer("/ec1").and_then(Value::as_str) == Some(ec1))
            });
        if !found {
            return Err(format!("registry {slug} lacks {ec1}"));
        }
    }
    if covers.pointer("/schema_version").and_then(Value::as_str) != Some("curve-covers/v1") {
        return Err("cover catalog schema mismatch".into());
    }
    let cover = covers
        .pointer("/curves")
        .and_then(Value::as_array)
        .and_then(|rows| {
            rows.iter()
                .find(|row| row.pointer("/slug").and_then(Value::as_str) == Some(CURVE_SLUG))
        })
        .ok_or("cover catalog lacks P-256")?;
    if cover
        .pointer("/certificate/construction")
        .and_then(Value::as_str)
        != Some("prime_quadratic_pullback_v1")
        || cover.pointer("/certificate/degree").and_then(Value::as_u64) != Some(2)
        || cover.pointer("/dlp_advantage").is_none()
        || !cover.pointer("/dlp_advantage").is_some_and(Value::is_null)
    {
        return Err("P-256 cover certificate mismatch".into());
    }
    if round300.pointer("/schema").and_then(Value::as_str) != Some("p256-cover-fiber-screen/v1")
        || round300.pointer("/curve").and_then(Value::as_str) != Some(CURVE_SLUG)
        || round300
            .pointer("/maximum_observed_quotient_rank_gain")
            .and_then(Value::as_u64)
            != Some(0)
    {
        return Err("Round 300 certificate mismatch".into());
    }
    if round301.pointer("/schema").and_then(Value::as_str)
        != Some("p256-algebraic-escape-census/v1")
        || round301.pointer("/curve").and_then(Value::as_str) != Some(CURVE_SLUG)
        || round301
            .pointer("/gates/parity_at_or_below_rho")
            .and_then(Value::as_bool)
            != Some(false)
    {
        return Err("Round 301 certificate mismatch".into());
    }
    Ok((
        vec![registry_dep, covers_dep, round300_dep, round301_dep],
        registry,
        covers,
        round300,
        round301,
    ))
}

fn boundary_table() -> Vec<BoundaryRow> {
    vec![
        BoundaryRow {
            variant: "Pollard rho".into(),
            independent_rows_per_cover_event: None,
            optimistic_ratio_to_rho: Some(1.0),
            complete_cost_measured: true,
            status: "reference".into(),
        },
        BoundaryRow {
            variant: format!("minimum-width independent-log {FACTOR_BASE_ID}"),
            independent_rows_per_cover_event: Some(1),
            optimistic_ratio_to_rho: Some(MINIMUM_WIDTH_RATIO_TO_RHO),
            complete_cost_measured: false,
            status: "free-perfect-oracle lower boundary".into(),
        },
        BoundaryRow {
            variant: format!("split-cover complementary quotient over {FACTOR_BASE_ID}"),
            independent_rows_per_cover_event: Some(1),
            optimistic_ratio_to_rho: Some(MINIMUM_WIDTH_RATIO_TO_RHO),
            complete_cost_measured: false,
            status: "complement has no rational order-n transport; no second row credit".into(),
        },
        BoundaryRow {
            variant: format!("registered 17-term {REGISTERED_S17_FACTOR_BASE_ID}"),
            independent_rows_per_cover_event: Some(1),
            optimistic_ratio_to_rho: Some(REGISTERED_S17_RATIO_TO_RHO),
            complete_cost_measured: false,
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

fn build_assessment(
    toy: &ToyControl,
    map: &P256MapControl,
    subgroup: &SubgroupCertificate,
) -> TransferAssessment {
    let mut assessment = TransferAssessment {
        schema: "p256-complementary-quotient-transfer-assessment/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 302,
        skill_profile: "transfer".into(),
        required_companion_resources_available: false,
        methodology_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/references/methodology.md (unavailable)".into(),
        template_resource: "skill://user/6ac50b6519048191b052e52c11310b8d/transfer/assets/assessment-template.json (unavailable)".into(),
        typed_correspondence_graph: vec![
            format!("{CURVE_SLUG}/E(F_p)[n] --[pullback]--> J(H)(F_p)"),
            format!("H: v^2=u^6+a*u^2+b --[pi_1:(u^2,v)]--> {CURVE_SLUG}/E"),
            "H --[pi_2:(b/u^2,b*v/u^3)]--> E': V^2=U^3+a*U^2+b^2".into(),
            "E(F_p)[n] --[composed rational homomorphism]--> E'(F_p)[n]=0".into(),
            "cover sheet multiplicity --[pi_1 quotient]--> kernel directions (Round 300)".into(),
        ],
        obligations: vec![
            obligation("explicit construction", "supported", "Both rational quotient maps and the complementary long Weierstrass model are explicit.", "This split degree-two cover only."),
            obligation("map correctness", "supported", &format!("{} complete toy primary/quartic replays and {} deterministic P-256 replays passed.", toy.primary_map_replays, map.accepted_cover_points), "Affine points; exceptional u=0 is counted separately."),
            obligation("subgroup preservation to primary quotient", "supported", "pi_1 is the registered H->E quotient map.", "P-256 order-n factor."),
            obligation("order-n transport to complementary quotient", "refuted", "Two independent scalar algorithms give [n]R != O on E', and the Hasse interval contains no other positive multiple of n.", "Rational E'(F_p) for this construction."),
            obligation("two independent collapsed rows", "refuted", "The complementary quotient has no rational order-n image; only the primary row survives.", "This split cover event model."),
            obligation("inverse recovery and complete cost below rho", "unknown", "No promoted second row exists, so no end-to-end attack was implemented.", "Future nonstandard covers remain open."),
        ],
        controls: vec![
            format!("Complete toy enumeration covered {} affine H points with {} exceptional u=0 points.", toy.affine_cover_points, toy.complementary_exceptional_points),
            format!("P-256 accepted {} deterministic points after scanning {} u values; digest {}.", map.accepted_cover_points, map.scanned_u_candidates, map.sample_sha256),
            format!("Independent LSB/MSB scalar replays match: {}.", subgroup.scalar_algorithms_match),
            "Round 300 and Round 301 were imported only after exact byte-hash, curve, schema, and gate checks.".into(),
        ],
        cost_accounting: vec![
            "Measured: exact field-operation, exponentiation, inversion, group-addition/doubling, scan, replay, process, and artifact counts.".into(),
            "Exact structural result: E'(F_p) has no rational subgroup of P-256 order n under the certified Hasse-and-witness argument.".into(),
            "Unset: structured relation solver, relation yield, sparse linear algebra, inverse DLP recovery, and complete attack cost because the second-row gate failed.".into(),
        ],
        exploration_boundary: "Exact for the explicit split quadratic cover H:v^2=u^6+a*u^2+b and its two elliptic quotients; not a search over all covers, Jacobians, or non-homomorphic relation mechanisms.".into(),
        weakest_open_obligation: "Construct a different P-256-specific cover or non-homomorphic relation mechanism whose additional collapsed equation remains in the order-n quotient, then implement recovery and complete below-rho cost.".into(),
        narrowest_supported_finding: "The complementary elliptic quotient of the registered split degree-two P-256 cover cannot supply a second rational order-n relation row; the 13.920747-times-rho boundary is unchanged for this construction.".into(),
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
    let (dependencies, _, _, _, _) = dependency_checks(cli)?;
    let curve = CurveParams::p256();
    if curve.h != 1 || curve.b.is_zero() {
        return Err("P-256 prerequisites failed".into());
    }
    let complementary = complement_curve(&curve)?;
    let toy = toy_control()?;
    let p256_maps = p256_map_control(&curve)?;
    let subgroup = subgroup_certificate(&curve)?;
    let zero_failures = toy.replay_failures == 0
        && p256_maps.replay_failures == 0
        && subgroup.scalar_algorithms_match
        && !subgroup.n_times_witness_is_infinity;
    if !zero_failures {
        return Err("a preregistered correctness control failed".into());
    }
    let gates = Gates {
        dependency_hashes_and_schemas_checked: true,
        zero_map_and_replay_failures: true,
        two_distinct_order_n_rows_per_cover_event: false,
        explicit_forward_inverse_kernel_and_recovery: false,
        structured_residual_degree_at_most_5: false,
        parity_at_or_below_rho: false,
        complete_collection_below_2_120: false,
        per_usable_relation_below_2_103: false,
        projected_storage_below_2_50: false,
        discarded_probabilistic_branches_counted_as_exhaustive: false,
        promoted: false,
    };
    let assessment = build_assessment(&toy, &p256_maps, &subgroup);
    let classification =
        "split-cover/complementary-quotient/order-n-transport-refuted/parity-blocked";
    let obstruction = "The explicit second elliptic quotient E' has no rational P-256 order-n subgroup: a deterministic E'(F_p) witness satisfies [n]R != O under two independent scalar replays, while the Hasse interval contains no positive multiple of n other than n. The cover therefore yields only the primary P-256 row, not two independent rows per event.";
    let decision = "Reject a two-row speedup from the registered split degree-two cover. Keep general nonstandard covers and non-homomorphic mechanisms open. Do not attempt an unplanted full-depth relation because the second-row, degree, cost, and storage gates fail.";
    let semantic = json!({
        "curve": CURVE_SLUG,
        "dependencies": &dependencies,
        "exact_maps": [
            "pi_1(u,v)=(u^2,v)",
            "pi_Q(u,v)=(u^2,u*v)",
            "pi_2(u,v)=(b/u^2,b*v/u^3)",
        ],
        "complementary_curve": &complementary,
        "toy_control": &toy,
        "p256_map_control": &p256_maps,
        "subgroup_certificate": &subgroup,
        "boundary_table": boundary_table(),
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
            schema: "p256-complementary-quotient-screen/v1".into(),
            curve: CURVE_SLUG.into(),
            screening_round: 302,
            execution_status: "complete".into(),
            factor_base: FACTOR_BASE_ID.into(),
            registered_s17_factor_base: REGISTERED_S17_FACTOR_BASE_ID.into(),
            dependencies,
            exact_maps: vec![
                "pi_1(u,v)=(u^2,v) on E:y^2=x^3+a*x+b".into(),
                "pi_Q(u,v)=(u^2,u*v) on Q:Y^2=X^4+a*X^2+b*X".into(),
                "pi_2(u,v)=(b/u^2,b*v/u^3) on E':V^2=U^3+a*U^2+b^2".into(),
            ],
            complementary_curve: complementary,
            toy_control: toy,
            p256_map_control: p256_maps,
            subgroup_certificate: subgroup,
            boundary_table: boundary_table(),
            distinct_p256_order_n_rows_per_cover_event: 1,
            structured_residual_degree_of_regularity: None,
            relations_reported_on_p256: 0,
            full_depth_unplanted_p256_relation_attempted: false,
            exploration_boundary: "Exact for the registered split quadratic cover and its complementary elliptic quotient; not universal over covers or Jacobians.".into(),
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
        "round 302: maps={}, [n]R_is_infinity={}, order-n rows/event={}, promoted={}",
        result.p256_map_control.accepted_cover_points,
        result.subgroup_certificate.n_times_witness_is_infinity,
        result.distinct_p256_order_n_rows_per_cover_event,
        result.gates.promoted
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_complementary_quotient: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_toy_control_replays() {
        let toy = toy_control().expect("toy control");
        assert!(toy.nonsingular);
        assert_eq!(toy.replay_failures, 0);
        assert_eq!(toy.primary_map_replays, toy.affine_cover_points);
        assert_eq!(
            toy.complementary_map_replays + toy.complementary_exceptional_points,
            toy.affine_cover_points
        );
    }

    #[test]
    fn scalar_algorithms_match_on_toy_complement() {
        let p = BigUint::from(1019u64);
        let a2 = BigUint::from(2u8);
        let b2 = BigUint::from(9u8);
        let point = Point::Affine {
            x: BigUint::zero(),
            y: BigUint::from(3u8),
        };
        let scalar = BigUint::from(123_457u64);
        let mut left_ops = OperationCounts::default();
        let mut right_ops = OperationCounts::default();
        let left = scalar_mul_lsb(&point, &scalar, &a2, &b2, &p, &mut left_ops).expect("lsb");
        let right = scalar_mul_msb(&point, &scalar, &a2, &b2, &p, &mut right_ops).expect("msb");
        assert_eq!(left, right);
    }

    #[test]
    fn p256_hasse_bound_has_only_one_possible_n_multiple() {
        let curve = CurveParams::p256();
        let hasse_upper = &curve.p + BigUint::one() + (BigUint::one() << 129usize);
        assert!(hasse_upper < (&curve.n << 1usize));
        assert_eq!(curve.h, 1);
    }
}
