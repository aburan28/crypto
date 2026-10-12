use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};

use crate::{
    encoding::{
        candidate_shard_digest, candidate_shell_digest, canonical_json, disposition_shard_digest,
        disposition_shell_digest, encode_candidate_record, encode_disposition_record, encode_sint,
        parse_canonical_json, sha256, ShardDigest,
    },
    factor_base::{
        kronecker_at_prime, roots_for_prime, verify_manifest, FactorBaseEntry, VerifiedFactorBase,
    },
    ideal::{evaluate_candidate_segmented, factor_u64, is_prime_u64},
    identity::{A_DEC, B_DEC, C_DEC, GX_HEX, GY_HEX},
    params::{self, CURVE_UID, LARGE_PRIME_BOUND},
    ref0::{pair_at, RegeneratedReference},
    Result,
};

const CONTROLS_SCHEMA: &str = "p192-wcm-controls-v1";
const RESULTS_DOMAIN: &[u8] = b"P192-WCM-CONTROL-RESULTS-v1\0";

const EXPECTED: &[(&str, &[&str])] = &[
    (
        "D23-ORDER3",
        &[
            "canonical_reduced_form_norm_2",
            "ideal_matches_form",
            "alpha_certificate_valid",
            "encoding_replays",
            "relation_g2_cubed_identity",
        ],
    ),
    (
        "P192-SQRT-D",
        &[
            "alpha_coordinates_minus_t_2",
            "norm_equals_abs_D",
            "factorization_equals_5_11_31_C",
            "algebra_passes",
            "C_degree_map_unavailable",
            "C_scalar_action_rejected",
        ],
    ),
    (
        "P192-RAMIFIED-SQUARES",
        &["ell_5_scalar", "ell_11_scalar", "ell_31_scalar"],
    ),
    (
        "INERT-FACTOR-REJECTION",
        &["base_injection_rejected", "residual_injection_rejected"],
    ),
    (
        "MUTATION-REJECTION",
        &[
            "root_mutation_rejected",
            "orientation_sign_mutation_rejected",
            "exponent_mutation_rejected",
            "alpha_coordinate_mutation_rejected",
            "norm_factor_mutation_rejected",
            "lp_endpoint_mutation_rejected",
            "hnf_column_mutation_rejected",
        ],
    ),
    (
        "LP-CANCELLATION",
        &["unmatched_orientation_rejected", "exact_conjugate_cancels"],
    ),
    (
        "REF0-SCALAR-SEGMENTED",
        &[
            "candidate_bytes_equal",
            "candidate_shard_root_equal",
            "candidate_shell_root_equal",
            "disposition_bytes_equal",
            "disposition_shard_root_equal",
            "disposition_shell_root_equal",
            "retained_records_equal",
        ],
    ),
    (
        "REF0-DIRECTION-WORKERS",
        &[
            "forward_reverse_candidate_root_equal",
            "forward_reverse_disposition_root_equal",
            "one_eight_worker_candidate_root_equal",
            "one_eight_worker_disposition_root_equal",
            "independent_in_memory_ref0_regenerations_streams_roots_counters_equal",
        ],
    ),
    (
        "LP-BOUNDARY",
        &[
            "q_L_minus_1_not_prime_rejected_as_lp",
            "q_L_prime_split_accepted_as_one_lp",
            "q_L_plus_1_out_of_range_rejected_as_lp",
            "residual_L_is_one_large_prime",
            "residual_L_squared_is_two_large_prime_with_multiplicity_two",
            "residual_L_squared_plus_one_is_rejected_by_inequality",
        ],
    ),
];

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AssertionRecord {
    pub assertion_id: String,
    pub passed: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ControlRecord {
    pub control_id: String,
    pub status: String,
    pub assertions: Vec<AssertionRecord>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ControlsManifest {
    pub schema: String,
    pub experiment_id: String,
    pub protocol_version: u64,
    pub protocol_commit: String,
    pub source_commit: String,
    pub curve_uid: String,
    pub factor_base_sha256: String,
    pub controls: Vec<ControlRecord>,
    pub control_results_sha256: String,
    pub overall_status: String,
}

#[derive(Clone, Debug)]
pub struct VerifiedControls {
    pub manifest: ControlsManifest,
    pub canonical_bytes: Vec<u8>,
    pub results_sha256: [u8; 32],
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct D23Fixture {
    discriminant: i64,
    trace: i64,
    norm_pi: i64,
    roots_mod_2: [i64; 2],
    selected_root: i64,
    positive_form: [i64; 3],
    conjugate_form: [i64; 3],
    principal_form: [i64; 3],
    alpha_u: i64,
    alpha_v: i64,
    alpha_x: i64,
    norm: i64,
    rational_prime: i64,
    rational_exponent: u32,
    signed_coordinate: i64,
    principal_matrix: [[i64; 2]; 2],
    ideal_power_matrix: [[i64; 2]; 2],
    signed_coordinate_sint_hex: String,
}

fn d23_fixture() -> D23Fixture {
    D23Fixture {
        discriminant: -23,
        trace: 1,
        norm_pi: 6,
        roots_mod_2: [0, 1],
        selected_root: 1,
        positive_form: [2, -1, 3],
        conjugate_form: [2, 1, 3],
        principal_form: [1, 1, 6],
        alpha_u: 1,
        alpha_v: 1,
        alpha_x: 3,
        norm: 8,
        rational_prime: 2,
        rational_exponent: 3,
        signed_coordinate: -3,
        principal_matrix: [[1, -6], [1, 2]],
        ideal_power_matrix: [[8, -7], [0, 1]],
        signed_coordinate_sint_hex: "020000000103".to_owned(),
    }
}

fn determinant(matrix: [[i64; 2]; 2]) -> i64 {
    matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
}

fn validate_d23(fixture: &D23Fixture) -> Result<()> {
    if fixture.trace * fixture.trace - 4 * fixture.norm_pi != fixture.discriminant {
        return Err("D=-23 polynomial discriminant mismatch".to_owned());
    }
    // The reduced primitive positive forms of discriminant -23 are exactly
    // (1,1,6), (2,1,3), and its inverse (2,-1,3).  The latter pair multiply
    // to the identity and squaring either yields its inverse, proving order 3.
    let mut forms = Vec::new();
    for a in 1i64..=3 {
        for b in -a..=a {
            if (b * b + 23) % (4 * a) != 0 {
                continue;
            }
            let c = (b * b + 23) / (4 * a);
            if a > c || (a == c || b.abs() == a) && b < 0 || a.gcd(&b).gcd(&c) != 1 {
                continue;
            }
            forms.push((a, b, c));
        }
    }
    forms.sort_unstable();
    if forms != vec![(1, 1, 6), (2, -1, 3), (2, 1, 3)]
        || fixture.positive_form != [2, -1, 3]
        || fixture.conjugate_form != [2, 1, 3]
        || fixture.principal_form != [1, 1, 6]
    {
        return Err("D=-23 reduced-form enumeration failed".to_owned());
    }
    let forms_from_roots: Vec<_> = fixture
        .roots_mod_2
        .iter()
        .map(|root| {
            let b = 2 * root - fixture.trace;
            let c = (b * b - fixture.discriminant) / 8;
            [2, b, c]
        })
        .collect();
    if forms_from_roots != vec![fixture.positive_form, fixture.conjugate_form] {
        return Err("D=-23 ideals do not match their root-derived forms".to_owned());
    }
    for root in fixture.roots_mod_2 {
        if (root * root - fixture.trace * root + fixture.norm_pi).rem_euclid(2) != 0 {
            return Err("D=-23 root fixture is not a polynomial root".to_owned());
        }
    }
    if fixture.roots_mod_2 != [0, 1]
        || !fixture.roots_mod_2.contains(&fixture.selected_root)
        || (fixture.alpha_u + fixture.alpha_v * fixture.selected_root).rem_euclid(2) != 0
        || fixture
            .roots_mod_2
            .iter()
            .filter(|root| (fixture.alpha_u + fixture.alpha_v * **root).rem_euclid(2) == 0)
            .count()
            != 1
    {
        return Err("D=-23 selected ideal/root does not divide alpha uniquely".to_owned());
    }
    let computed_x = 2 * fixture.alpha_u + fixture.trace * fixture.alpha_v;
    let computed_norm = fixture.alpha_u * fixture.alpha_u
        + fixture.trace * fixture.alpha_u * fixture.alpha_v
        + fixture.norm_pi * fixture.alpha_v * fixture.alpha_v;
    if fixture.alpha_x != computed_x || fixture.norm != computed_norm {
        return Err("D=-23 principal alpha certificate failed".to_owned());
    }
    let rational_norm = fixture
        .rational_prime
        .checked_pow(fixture.rational_exponent)
        .ok_or_else(|| "D=-23 factor exponent overflow".to_owned())?;
    let expected_coordinate = if fixture.selected_root == fixture.roots_mod_2[0] {
        i64::from(fixture.rational_exponent)
    } else {
        -i64::from(fixture.rational_exponent)
    };
    if rational_norm != fixture.norm || fixture.signed_coordinate != expected_coordinate {
        return Err("D=-23 norm factor/oriented exponent mismatch".to_owned());
    }
    let mut encoded = Vec::new();
    encode_sint(&BigInt::from(fixture.signed_coordinate), &mut encoded)?;
    if hex::encode(encoded) != fixture.signed_coordinate_sint_hex {
        return Err("D=-23 encoding replay failed".to_owned());
    }
    let expected_principal = [
        [fixture.alpha_u, -fixture.norm_pi * fixture.alpha_v],
        [
            fixture.alpha_v,
            fixture.alpha_u + fixture.trace * fixture.alpha_v,
        ],
    ];
    let lifted_roots: Vec<_> = (0..fixture.norm)
        .filter(|root| {
            root.rem_euclid(2) == fixture.selected_root
                && (root * root - fixture.trace * root + fixture.norm_pi).rem_euclid(fixture.norm)
                    == 0
                && (fixture.alpha_u + fixture.alpha_v * root).rem_euclid(fixture.norm) == 0
        })
        .collect();
    if lifted_roots.len() != 1 {
        return Err("D=-23 ideal root does not lift uniquely modulo norm".to_owned());
    }
    let expected_ideal_power = [[fixture.norm, -lifted_roots[0]], [0, 1]];
    if fixture.principal_matrix != expected_principal
        || fixture.ideal_power_matrix != expected_ideal_power
        || determinant(fixture.principal_matrix) != fixture.norm
        || determinant(fixture.ideal_power_matrix) != fixture.norm
        || [
            fixture.ideal_power_matrix[0][0] + fixture.ideal_power_matrix[0][1],
            fixture.ideal_power_matrix[1][0] + fixture.ideal_power_matrix[1][1],
        ] != [fixture.alpha_u, fixture.alpha_v]
        || [
            fixture.ideal_power_matrix[0][0] + 2 * fixture.ideal_power_matrix[0][1],
            fixture.ideal_power_matrix[1][0] + 2 * fixture.ideal_power_matrix[1][1],
        ] != [expected_principal[0][1], expected_principal[1][1]]
    {
        return Err("D=-23 ideal/principal HNF equality failed".to_owned());
    }
    Ok(())
}

fn verify_d23() -> Result<()> {
    validate_d23(&d23_fixture())
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MutationAlpha {
    u: String,
    v: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MutationRationalFactor {
    ell: String,
    exponent: u32,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct D23MutationTarget {
    selected_root: i64,
    signed_coordinate: String,
    alpha: MutationAlpha,
    norm: String,
    rational_factors: Vec<MutationRationalFactor>,
    ideal_power_column_matrix_rows: [[i64; 2]; 2],
    signed_coordinate_sint_hex: String,
}

fn d23_mutation_target() -> D23MutationTarget {
    D23MutationTarget {
        selected_root: 1,
        signed_coordinate: "-3".to_owned(),
        alpha: MutationAlpha {
            u: "1".to_owned(),
            v: "1".to_owned(),
        },
        norm: "8".to_owned(),
        rational_factors: vec![MutationRationalFactor {
            ell: "2".to_owned(),
            exponent: 3,
        }],
        ideal_power_column_matrix_rows: [[8, -7], [0, 1]],
        signed_coordinate_sint_hex: "020000000103".to_owned(),
    }
}

fn d23_from_mutation_target(target: &D23MutationTarget) -> Result<D23Fixture> {
    if target.rational_factors.len() != 1 || target.rational_factors[0].ell != "2" {
        return Err("D=-23 mutation target has a noncanonical rational-factor list".to_owned());
    }
    let mut fixture = d23_fixture();
    fixture.selected_root = target.selected_root;
    fixture.signed_coordinate = target
        .signed_coordinate
        .parse()
        .map_err(|_| "D=-23 signed coordinate is not an i64".to_owned())?;
    fixture.alpha_u = target
        .alpha
        .u
        .parse()
        .map_err(|_| "D=-23 alpha.u is not an i64".to_owned())?;
    fixture.alpha_v = target
        .alpha
        .v
        .parse()
        .map_err(|_| "D=-23 alpha.v is not an i64".to_owned())?;
    fixture.norm = target
        .norm
        .parse()
        .map_err(|_| "D=-23 norm is not an i64".to_owned())?;
    fixture.rational_exponent = target.rational_factors[0].exponent;
    fixture.ideal_power_matrix = target.ideal_power_column_matrix_rows;
    fixture.signed_coordinate_sint_hex = target.signed_coordinate_sint_hex.clone();
    Ok(fixture)
}

#[derive(Clone, Debug)]
struct RamifiedSquareFixture {
    ell: u64,
    root: u64,
    alpha_u: BigInt,
    alpha_v: BigInt,
    claimed_norm: BigInt,
    claimed_non_scalar: bool,
}

fn validate_ramified_square(
    fixture: &RamifiedSquareFixture,
    factor_base: &VerifiedFactorBase,
) -> Result<()> {
    let entry = factor_base
        .manifest
        .entries
        .iter()
        .find(|entry| entry.ell == fixture.ell)
        .ok_or_else(|| format!("ramified control prime {} absent", fixture.ell))?;
    if entry.kind != 1
        || entry.kronecker != 0
        || entry.roots != vec![fixture.root]
        || entry.positive_root_index != 0
    {
        return Err(format!("ramified control prime {} malformed", fixture.ell));
    }
    let root = BigInt::from(fixture.root);
    let modulus = BigInt::from(fixture.ell);
    let polynomial = &root * &root - params::t() * &root + params::p();
    let derivative = 2u8 * &root - params::t();
    if !polynomial.mod_floor(&modulus).is_zero()
        || !derivative.mod_floor(&modulus).is_zero()
        || !params::d().mod_floor(&modulus).is_zero()
        || (params::d() / &modulus).mod_floor(&modulus).is_zero()
    {
        return Err(format!(
            "ramified ideal at {} is not squarefree/repeated-root",
            fixture.ell
        ));
    }
    let computed_norm = &fixture.alpha_u * &fixture.alpha_u
        + params::t() * &fixture.alpha_u * &fixture.alpha_v
        + params::p() * &fixture.alpha_v * &fixture.alpha_v;
    let content = fixture.alpha_u.abs().gcd(&fixture.alpha_v);
    let scalar = fixture.alpha_v.is_zero() && content == BigInt::from(fixture.ell);
    if fixture.alpha_u != BigInt::from(fixture.ell)
        || fixture.claimed_norm != computed_norm
        || fixture.claimed_norm != BigInt::from(fixture.ell).pow(2u32)
        || !scalar
        || fixture.claimed_non_scalar
    {
        return Err(format!(
            "ramified square {} was not classified as a removed scalar",
            fixture.ell
        ));
    }
    Ok(())
}

fn verify_sqrt_d(factor_base: &VerifiedFactorBase) -> Result<()> {
    let t = params::t();
    let u = -&t;
    let v = BigInt::from(2u8);
    let norm = &u * &u + &t * &u * &v + params::p() * &v * &v;
    let abs_d = -params::d();
    if norm != abs_d || u.abs().gcd(&v) != BigInt::one() {
        return Err("P-192 sqrt(D) alpha arithmetic failed".to_owned());
    }
    let c = params::bigint(C_DEC, "C")?;
    let subgroup_order = params::subgroup_order();
    if abs_d != BigInt::from(5u8) * 11u8 * 31u8 * &c {
        return Err("P-192 sqrt(D) factorization failed".to_owned());
    }
    let c_in_factor_base = factor_base
        .manifest
        .entries
        .iter()
        .any(|entry| BigInt::from(entry.ell) == c);
    let action =
        derive_sqrt_d_subgroup_action(&u, &v, &subgroup_order, &c, c_in_factor_base, false)?;
    if action.action_admitted {
        return Err("unavailable C-degree edge/scalar action was not rejected".to_owned());
    }
    for prime in [5u64, 11, 31] {
        let entry = factor_base
            .manifest
            .entries
            .iter()
            .find(|entry| entry.ell == prime)
            .ok_or_else(|| format!("ramified control prime {prime} absent"))?;
        let root = *entry
            .roots
            .first()
            .ok_or_else(|| format!("ramified control prime {prime} lacks a root"))?;
        validate_ramified_square(
            &RamifiedSquareFixture {
                ell: prime,
                root,
                alpha_u: BigInt::from(prime),
                alpha_v: BigInt::zero(),
                claimed_norm: BigInt::from(prime).pow(2u32),
                claimed_non_scalar: false,
            },
            factor_base,
        )?;
    }
    Ok(())
}

const SQRT_D_SUBGROUP_SCALAR_DEC: &str =
    "6277101735386680763835789423144451611450480845975359606884";

#[derive(Clone, Debug, PartialEq, Eq)]
struct SqrtDActionAssessment {
    scalar_action: BigInt,
    action_admitted: bool,
}

fn p192_generator_is_fixed_by_frobenius() -> Result<bool> {
    let p = params::biguint(params::P_DEC, "p")?;
    let a = params::biguint(A_DEC, "a")?;
    let b = params::biguint(B_DEC, "b")?;
    let x = BigUint::parse_bytes(GX_HEX.as_bytes(), 16)
        .ok_or_else(|| "invalid frozen P-192 generator x".to_owned())?;
    let y = BigUint::parse_bytes(GY_HEX.as_bytes(), 16)
        .ok_or_else(|| "invalid frozen P-192 generator y".to_owned())?;
    if x >= p || y >= p {
        return Ok(false);
    }
    let on_curve =
        y.modpow(&BigUint::from(2u8), &p) == (x.modpow(&BigUint::from(3u8), &p) + &a * &x + b) % &p;
    let frobenius_x = x.modpow(&p, &p);
    let frobenius_y = y.modpow(&p, &p);
    Ok(on_curve && frobenius_x == x && frobenius_y == y)
}

fn derive_sqrt_d_subgroup_action(
    u: &BigInt,
    v: &BigInt,
    subgroup_order: &BigInt,
    c: &BigInt,
    c_in_factor_base: bool,
    c_degree_map_available: bool,
) -> Result<SqrtDActionAssessment> {
    let p = params::p();
    let t = params::t();
    let d = params::d();
    if subgroup_order <= &BigInt::one() || &p + 1u8 - &t != *subgroup_order {
        return Err("sqrt(D) subgroup order does not equal p+1-t".to_owned());
    }
    if u != &-&t || v != &BigInt::from(2u8) {
        return Err("sqrt(D) coordinates are not (-t,2) in the (1,pi) basis".to_owned());
    }
    let norm = u * u + &t * u * v + &p * v * v;
    if norm != -&d {
        return Err("sqrt(D) coordinate norm is not |D|".to_owned());
    }

    // Every point in the named subgroup is F_p-rational, so arithmetic
    // Frobenius fixes it.  Check that fact on the independently compiled
    // generator, then check lambda_pi=1 in X^2-tX+p modulo n rather than
    // assuming the endomorphism's scalar label.
    if !p192_generator_is_fixed_by_frobenius()? {
        return Err("P-192 generator is not fixed by p-power Frobenius".to_owned());
    }
    let lambda_pi = BigInt::one();
    let characteristic = &lambda_pi * &lambda_pi - &t * &lambda_pi + &p;
    if !characteristic.mod_floor(subgroup_order).is_zero() {
        return Err(
            "Frobenius scalar 1 does not satisfy its characteristic polynomial mod n".to_owned(),
        );
    }

    let scalar_action = (u + v * &lambda_pi).mod_floor(subgroup_order);
    let expected = params::bigint(SQRT_D_SUBGROUP_SCALAR_DEC, "sqrt(D) subgroup scalar")?;
    if scalar_action <= BigInt::zero()
        || scalar_action >= *subgroup_order
        || scalar_action != expected
        || !(&scalar_action * &scalar_action - &d)
            .mod_floor(subgroup_order)
            .is_zero()
    {
        return Err("derived sqrt(D) subgroup scalar action is inconsistent".to_owned());
    }

    let large_prime_bound = BigInt::from(params::LARGE_PRIME_BOUND);
    let large_prime_bound_squared = &large_prime_bound * &large_prime_bound;
    if c <= &BigInt::from(params::MAPPABLE_BOUND)
        || c <= &large_prime_bound_squared
        || c_in_factor_base
        || c_degree_map_available
    {
        return Err("C-degree sqrt(D) action unexpectedly became mappable or available".to_owned());
    }
    Ok(SqrtDActionAssessment {
        scalar_action,
        action_admitted: false,
    })
}

fn verify_inert_rejection(factor_base: &VerifiedFactorBase) -> Result<()> {
    if kronecker_at_prime(&params::d(), 2) != -1
        || factor_base
            .manifest
            .entries
            .iter()
            .any(|entry| entry.ell == 2)
    {
        return Err("inert prime 2 was not rejected from the factor base".to_owned());
    }
    let mut injected = factor_base.manifest.clone();
    injected.entries.insert(
        0,
        FactorBaseEntry {
            ell: 2,
            kind: 0,
            kronecker: -1,
            roots: Vec::new(),
            positive_root_index: 0,
        },
    );
    injected.entry_count += 1;
    if verify_manifest(&injected).is_ok() {
        return Err("inert factor-base injection was accepted".to_owned());
    }

    let residual_inert = (params::ALGEBRAIC_BOUND + 1..)
        .find(|candidate| {
            is_prime_u64(*candidate) && kronecker_at_prime(&params::d(), *candidate) < 0
        })
        .ok_or_else(|| "could not construct inert residual control".to_owned())?;
    if residual_inert > LARGE_PRIME_BOUND
        || validate_retained_residual_prime(residual_inert).is_ok()
    {
        return Err("inert residual injection was not rejected inside the LP bound".to_owned());
    }
    Ok(())
}

fn validate_retained_residual_prime(prime: u64) -> Result<()> {
    if prime <= params::ALGEBRAIC_BOUND
        || prime > LARGE_PRIME_BOUND
        || !is_prime_u64(prime)
        || kronecker_at_prime(&params::d(), prime) != 1
        || roots_for_prime(prime, &params::t(), &params::p(), &params::d())?.len() != 2
    {
        return Err("residual prime is not an admitted split LP".to_owned());
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct LpFixture {
    q: u64,
    discriminant_residue: u64,
    smaller_root: u64,
    larger_root: u64,
    sign: u8,
    multiplicity: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct LpMutationTarget {
    q: u64,
    smaller_root: u64,
    larger_root: u64,
    sign: u8,
    multiplicity: u64,
}

fn lp_mutation_target() -> LpMutationTarget {
    let fixture = lp_fixture(1);
    LpMutationTarget {
        q: fixture.q,
        smaller_root: fixture.smaller_root,
        larger_root: fixture.larger_root,
        sign: fixture.sign,
        multiplicity: fixture.multiplicity,
    }
}

fn lp_from_mutation_target(target: &LpMutationTarget) -> LpFixture {
    LpFixture {
        q: target.q,
        discriminant_residue: 2_121_375_023,
        smaller_root: target.smaller_root,
        larger_root: target.larger_root,
        sign: target.sign,
        multiplicity: target.multiplicity,
    }
}

fn lp_fixture(sign: u8) -> LpFixture {
    LpFixture {
        q: LARGE_PRIME_BOUND,
        discriminant_residue: 2_121_375_023,
        smaller_root: 1_109_020_142,
        larger_root: 2_025_854_371,
        sign,
        multiplicity: 1,
    }
}

fn validate_lp(fixture: &LpFixture) -> Result<()> {
    if validate_retained_residual_prime(fixture.q).is_err()
        || !matches!(fixture.sign, 1 | 2)
        || fixture.multiplicity == 0
    {
        return Err("large-prime endpoint envelope is invalid".to_owned());
    }
    let discriminant_residue = params::d().mod_floor(&BigInt::from(fixture.q));
    let discriminant_residue = u64::try_from(discriminant_residue)
        .map_err(|_| "LP discriminant residue conversion failed".to_owned())?;
    let roots = roots_for_prime(fixture.q, &params::t(), &params::p(), &params::d())?;
    if discriminant_residue != fixture.discriminant_residue
        || roots != vec![fixture.smaller_root, fixture.larger_root]
        || fixture.smaller_root >= fixture.larger_root
    {
        return Err("large-prime endpoint roots/discriminant mismatch".to_owned());
    }
    Ok(())
}

const MUTATION_FILE_SCHEMA: &str = "p192-wcm-implementation-mutation-file-v1";
const MUTATION_RECORD_DOMAIN: &[u8] = b"P192-WCM-MUTATION-RECORD-v1\0";
const MUTATION_FILE_DOMAIN: &[u8] = b"P192-WCM-MUTATION-FILE-v1\0";

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MutationRecord {
    mutation_id: String,
    fixture_kind: String,
    target_field_path: String,
    fixture: serde_json::Value,
    record_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MutationFile {
    schema: String,
    records: Vec<MutationRecord>,
    file_sha256: String,
}

#[derive(Clone, Debug)]
struct MutationSpec {
    mutation_id: &'static str,
    fixture_kind: &'static str,
    target_field_path: &'static str,
    carrier_field_path: &'static str,
    before: serde_json::Value,
    after: serde_json::Value,
    base_fixture: serde_json::Value,
    mutated_fixture: serde_json::Value,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct MutationReplayEvidence {
    cases: u64,
    mutation_identities: Vec<(&'static str, &'static str)>,
    exact_single_scalar_changes: u64,
    stale_record_hash_rejections: u64,
    stale_file_hash_rejections: u64,
    recomputed_hash_chain_acceptances: u64,
    semantic_rejections_after_hash_acceptance: u64,
}

const EXPECTED_MUTATION_IDENTITIES: &[(&str, &str)] = &[
    (
        "root",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/selected_root",
    ),
    (
        "orientation_sign",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate",
    ),
    (
        "exponent",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate",
    ),
    (
        "alpha_u",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/alpha/u",
    ),
    (
        "norm",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/norm",
    ),
    (
        "rational_factor_exponent",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/rational_factors/0/exponent",
    ),
    (
        "lp_smaller_root",
        "/role1_interface_addendum/role1_control_fixtures/LP-BOUNDARY/endpoint_fixtures/1/smaller_root",
    ),
    (
        "hnf_matrix_entry",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/ideal_power_column_matrix_rows/0/1",
    ),
    (
        "sint_sign_byte",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate_sint_hex",
    ),
];

fn mutation_record_sha256(record: &MutationRecord) -> Result<String> {
    let projection = serde_json::json!({
        "fixture": &record.fixture,
        "fixture_kind": &record.fixture_kind,
        "mutation_id": &record.mutation_id,
        "target_field_path": &record.target_field_path,
    });
    let mut bytes = MUTATION_RECORD_DOMAIN.to_vec();
    bytes.extend_from_slice(&canonical_json(&projection)?);
    Ok(hex::encode(sha256(&bytes)))
}

fn mutation_file_sha256(file: &MutationFile) -> Result<String> {
    let projection = serde_json::json!({
        "records": &file.records,
        "schema": &file.schema,
    });
    let mut bytes = MUTATION_FILE_DOMAIN.to_vec();
    bytes.extend_from_slice(&canonical_json(&projection)?);
    Ok(hex::encode(sha256(&bytes)))
}

fn seal_mutation_file(spec: &MutationSpec, fixture: serde_json::Value) -> Result<Vec<u8>> {
    let mut record = MutationRecord {
        mutation_id: spec.mutation_id.to_owned(),
        fixture_kind: spec.fixture_kind.to_owned(),
        target_field_path: spec.target_field_path.to_owned(),
        fixture,
        record_sha256: String::new(),
    };
    record.record_sha256 = mutation_record_sha256(&record)?;
    let mut file = MutationFile {
        schema: MUTATION_FILE_SCHEMA.to_owned(),
        records: vec![record],
        file_sha256: String::new(),
    };
    file.file_sha256 = mutation_file_sha256(&file)?;
    canonical_json(&serde_json::to_value(file).map_err(|error| error.to_string())?)
}

fn verify_mutation_hash_chain(bytes: &[u8]) -> Result<MutationFile> {
    let value = parse_canonical_json(bytes)?;
    let file: MutationFile =
        serde_json::from_value(value).map_err(|error| format!("mutation file schema: {error}"))?;
    if file.schema != MUTATION_FILE_SCHEMA || file.records.len() != 1 {
        return Err("mutation file envelope is noncanonical".to_owned());
    }
    let record = &file.records[0];
    if record.record_sha256 != mutation_record_sha256(record)? {
        return Err("mutation record hash mismatch".to_owned());
    }
    if file.file_sha256 != mutation_file_sha256(&file)? {
        return Err("mutation file hash mismatch".to_owned());
    }
    Ok(file)
}

fn json_pointer_escape(key: &str) -> String {
    key.replace('~', "~0").replace('/', "~1")
}

fn scalar_differences(
    before: &serde_json::Value,
    after: &serde_json::Value,
    path: &str,
    differences: &mut Vec<(String, serde_json::Value, serde_json::Value)>,
) -> Result<()> {
    match (before, after) {
        (serde_json::Value::Object(left), serde_json::Value::Object(right)) => {
            if left.keys().collect::<Vec<_>>() != right.keys().collect::<Vec<_>>() {
                return Err("mutation changed the carrier object shape".to_owned());
            }
            for (key, left_value) in left {
                let right_value = right
                    .get(key)
                    .ok_or_else(|| "mutation removed a carrier field".to_owned())?;
                scalar_differences(
                    left_value,
                    right_value,
                    &format!("{path}/{}", json_pointer_escape(key)),
                    differences,
                )?;
            }
        }
        (serde_json::Value::Array(left), serde_json::Value::Array(right)) => {
            if left.len() != right.len() {
                return Err("mutation changed the carrier array length".to_owned());
            }
            for (index, (left_value, right_value)) in left.iter().zip(right).enumerate() {
                scalar_differences(
                    left_value,
                    right_value,
                    &format!("{path}/{index}"),
                    differences,
                )?;
            }
        }
        (left, right)
            if left.is_object() || left.is_array() || right.is_object() || right.is_array() =>
        {
            return Err("mutation changed a scalar into a container or vice versa".to_owned());
        }
        (left, right) if left != right => {
            differences.push((path.to_owned(), left.clone(), right.clone()));
        }
        _ => {}
    }
    Ok(())
}

fn exact_frozen_scalar_change(spec: &MutationSpec) -> Result<()> {
    let mut differences = Vec::new();
    scalar_differences(
        &spec.base_fixture,
        &spec.mutated_fixture,
        "",
        &mut differences,
    )?;
    if differences
        != vec![(
            spec.carrier_field_path.to_owned(),
            spec.before.clone(),
            spec.after.clone(),
        )]
    {
        return Err(format!(
            "mutation {} did not change exactly its frozen scalar",
            spec.mutation_id
        ));
    }
    Ok(())
}

fn mutation_semantically_rejected(record: &MutationRecord) -> Result<bool> {
    match record.fixture_kind.as_str() {
        "D23-ORDER3" => {
            let target: D23MutationTarget = serde_json::from_value(record.fixture.clone())
                .map_err(|error| format!("D=-23 mutation target schema: {error}"))?;
            Ok(d23_from_mutation_target(&target)
                .and_then(|fixture| validate_d23(&fixture))
                .is_err())
        }
        "LP-BOUNDARY" => {
            let target: LpMutationTarget = serde_json::from_value(record.fixture.clone())
                .map_err(|error| format!("LP mutation target schema: {error}"))?;
            Ok(validate_lp(&lp_from_mutation_target(&target)).is_err())
        }
        _ => Err("unknown mutation fixture kind".to_owned()),
    }
}

fn mutation_specs(
    d23_base: D23MutationTarget,
    lp_base: LpMutationTarget,
) -> Result<Vec<MutationSpec>> {
    let d23_value = serde_json::to_value(&d23_base).map_err(|error| error.to_string())?;
    let mut specs = Vec::new();
    let push_d23 = |specs: &mut Vec<MutationSpec>,
                    mutation_id,
                    target_field_path,
                    carrier_field_path,
                    before,
                    after,
                    mutated: D23MutationTarget|
     -> Result<()> {
        specs.push(MutationSpec {
            mutation_id,
            fixture_kind: "D23-ORDER3",
            target_field_path,
            carrier_field_path,
            before,
            after,
            base_fixture: d23_value.clone(),
            mutated_fixture: serde_json::to_value(mutated).map_err(|error| error.to_string())?,
        });
        Ok(())
    };

    let mut mutated = d23_base.clone();
    mutated.selected_root = 0;
    push_d23(
        &mut specs,
        "root",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/selected_root",
        "/selected_root",
        serde_json::json!(1),
        serde_json::json!(0),
        mutated,
    )?;
    let mut mutated = d23_base.clone();
    mutated.signed_coordinate = "3".to_owned();
    push_d23(
        &mut specs,
        "orientation_sign",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate",
        "/signed_coordinate",
        serde_json::json!("-3"),
        serde_json::json!("3"),
        mutated,
    )?;
    let mut mutated = d23_base.clone();
    mutated.signed_coordinate = "-2".to_owned();
    push_d23(
        &mut specs,
        "exponent",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate",
        "/signed_coordinate",
        serde_json::json!("-3"),
        serde_json::json!("-2"),
        mutated,
    )?;
    let mut mutated = d23_base.clone();
    mutated.alpha.u = "2".to_owned();
    push_d23(
        &mut specs,
        "alpha_u",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/alpha/u",
        "/alpha/u",
        serde_json::json!("1"),
        serde_json::json!("2"),
        mutated,
    )?;
    let mut mutated = d23_base.clone();
    mutated.norm = "4".to_owned();
    push_d23(
        &mut specs,
        "norm",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/norm",
        "/norm",
        serde_json::json!("8"),
        serde_json::json!("4"),
        mutated,
    )?;
    let mut mutated = d23_base.clone();
    mutated.rational_factors[0].exponent = 2;
    push_d23(
        &mut specs,
        "rational_factor_exponent",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/rational_factors/0/exponent",
        "/rational_factors/0/exponent",
        serde_json::json!(3),
        serde_json::json!(2),
        mutated,
    )?;

    let lp_value = serde_json::to_value(&lp_base).map_err(|error| error.to_string())?;
    let mut mutated = lp_base;
    mutated.smaller_root = 1_109_020_143;
    specs.push(MutationSpec {
        mutation_id: "lp_smaller_root",
        fixture_kind: "LP-BOUNDARY",
        target_field_path: "/role1_interface_addendum/role1_control_fixtures/LP-BOUNDARY/endpoint_fixtures/1/smaller_root",
        carrier_field_path: "/smaller_root",
        before: serde_json::json!(1_109_020_142u64),
        after: serde_json::json!(1_109_020_143u64),
        base_fixture: lp_value,
        mutated_fixture: serde_json::to_value(mutated).map_err(|error| error.to_string())?,
    });

    let mut mutated = d23_base.clone();
    mutated.ideal_power_column_matrix_rows[0][1] = -6;
    push_d23(
        &mut specs,
        "hnf_matrix_entry",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/ideal_power_column_matrix_rows/0/1",
        "/ideal_power_column_matrix_rows/0/1",
        serde_json::json!(-7),
        serde_json::json!(-6),
        mutated,
    )?;
    let mut mutated = d23_base;
    mutated.signed_coordinate_sint_hex = "010000000103".to_owned();
    push_d23(
        &mut specs,
        "sint_sign_byte",
        "/role1_interface_addendum/role1_control_fixtures/D23-ORDER3/mutation_target_record/signed_coordinate_sint_hex",
        "/signed_coordinate_sint_hex",
        serde_json::json!("020000000103"),
        serde_json::json!("010000000103"),
        mutated,
    )?;

    Ok(specs)
}

fn replay_mutation_rejections() -> Result<MutationReplayEvidence> {
    replay_mutation_rejections_with_bases(d23_mutation_target(), lp_mutation_target())
}

fn replay_mutation_rejections_with_bases(
    d23_base: D23MutationTarget,
    lp_base: LpMutationTarget,
) -> Result<MutationReplayEvidence> {
    validate_d23(&d23_fixture())?;
    validate_lp(&lp_fixture(1))?;
    validate_d23(&d23_from_mutation_target(&d23_base)?)
        .map_err(|error| format!("frozen D23 mutation base fixture failed replay: {error}"))?;
    validate_lp(&lp_from_mutation_target(&lp_base))
        .map_err(|error| format!("frozen LP mutation base fixture failed replay: {error}"))?;
    let mut evidence = MutationReplayEvidence::default();
    for spec in mutation_specs(d23_base, lp_base)? {
        evidence.cases += 1;
        exact_frozen_scalar_change(&spec)?;
        evidence.exact_single_scalar_changes += 1;

        let base_bytes = seal_mutation_file(&spec, spec.base_fixture.clone())?;
        let mut stale_record: MutationFile = serde_json::from_slice(&base_bytes)
            .map_err(|error| format!("parse sealed mutation file: {error}"))?;
        stale_record.records[0].fixture = spec.mutated_fixture.clone();
        let stale_record_bytes = canonical_json(
            &serde_json::to_value(stale_record).map_err(|error| error.to_string())?,
        )?;
        match verify_mutation_hash_chain(&stale_record_bytes) {
            Err(error) if error == "mutation record hash mismatch" => {}
            Ok(_) => {
                return Err(format!(
                    "mutation {} passed with a stale record hash",
                    spec.mutation_id
                ));
            }
            Err(error) => {
                return Err(format!(
                    "mutation {} stale record hash failed at the wrong gate: {error}",
                    spec.mutation_id
                ));
            }
        }
        evidence.stale_record_hash_rejections += 1;

        let mut stale_file: MutationFile = serde_json::from_slice(&base_bytes)
            .map_err(|error| format!("parse sealed mutation file: {error}"))?;
        stale_file.records[0].fixture = spec.mutated_fixture.clone();
        stale_file.records[0].record_sha256 = mutation_record_sha256(&stale_file.records[0])?;
        let stale_file_bytes =
            canonical_json(&serde_json::to_value(stale_file).map_err(|error| error.to_string())?)?;
        match verify_mutation_hash_chain(&stale_file_bytes) {
            Err(error) if error == "mutation file hash mismatch" => {}
            Ok(_) => {
                return Err(format!(
                    "mutation {} passed with a stale file hash",
                    spec.mutation_id
                ));
            }
            Err(error) => {
                return Err(format!(
                    "mutation {} stale file hash failed at the wrong gate: {error}",
                    spec.mutation_id
                ));
            }
        }
        evidence.stale_file_hash_rejections += 1;

        let rehashed_bytes = seal_mutation_file(&spec, spec.mutated_fixture.clone())?;
        let rehashed = verify_mutation_hash_chain(&rehashed_bytes)?;
        evidence.recomputed_hash_chain_acceptances += 1;
        let record = &rehashed.records[0];
        if record.mutation_id != spec.mutation_id
            || record.fixture_kind != spec.fixture_kind
            || record.target_field_path != spec.target_field_path
            || record.fixture != spec.mutated_fixture
        {
            return Err(format!(
                "mutation {} changed identity after hashing",
                spec.mutation_id
            ));
        }
        evidence
            .mutation_identities
            .push((spec.mutation_id, spec.target_field_path));
        if !mutation_semantically_rejected(record)? {
            return Err(format!(
                "mutation {} passed semantic replay after its hash chain was recomputed",
                spec.mutation_id
            ));
        }
        evidence.semantic_rejections_after_hash_acceptance += 1;
    }
    Ok(evidence)
}

fn verify_mutation_rejection() -> Result<()> {
    let evidence = replay_mutation_rejections()?;
    let expected = evidence.cases;
    if expected != 9
        || evidence.mutation_identities != EXPECTED_MUTATION_IDENTITIES
        || evidence.exact_single_scalar_changes != expected
        || evidence.stale_record_hash_rejections != expected
        || evidence.stale_file_hash_rejections != expected
        || evidence.recomputed_hash_chain_acceptances != expected
        || evidence.semantic_rejections_after_hash_acceptance != expected
    {
        return Err("mutation replay evidence is incomplete".to_owned());
    }
    Ok(())
}

fn validate_lp_cancellation(entries: &[LpFixture], require_zero: bool) -> Result<()> {
    let first = entries
        .first()
        .ok_or_else(|| "LP cancellation fixture is empty".to_owned())?;
    let mut coordinate = 0i64;
    for entry in entries {
        validate_lp(entry)?;
        if entry.q != first.q
            || entry.smaller_root != first.smaller_root
            || entry.larger_root != first.larger_root
        {
            return Err("LP cancellation mixed nonconjugate endpoints".to_owned());
        }
        let multiplicity = i64::try_from(entry.multiplicity)
            .map_err(|_| "LP multiplicity does not fit signed accumulator".to_owned())?;
        coordinate += if entry.sign == 1 {
            multiplicity
        } else {
            -multiplicity
        };
    }
    if (coordinate == 0) != require_zero {
        return Err("LP cancellation result differs from its required boundary state".to_owned());
    }
    Ok(())
}

fn verify_lp_cancellation() -> Result<()> {
    validate_lp_cancellation(&[lp_fixture(1)], false)?;
    validate_lp_cancellation(&[lp_fixture(1), lp_fixture(2)], true)
}

struct AlternateReference {
    candidate_records: Vec<u8>,
    disposition_records: Vec<u8>,
    certificate_store: Vec<u8>,
    status_counts: [u64; 7],
    candidate_shard: [u8; 32],
    candidate_shell: [u8; 32],
    disposition_shard: [u8; 32],
    disposition_shell: [u8; 32],
}

fn segmented_factor_marks(factor_base: &VerifiedFactorBase) -> Result<Vec<Vec<usize>>> {
    const SEGMENT_RECORDS: u64 = 128;
    let mut candidates = Vec::with_capacity(params::REF0_COUNT as usize);
    for index in 0..params::REF0_COUNT {
        let (v_box, x_box) = pair_at(index)?;
        let v = BigInt::from(v_box);
        let numerator = BigInt::from(x_box) - params::t() * &v;
        if numerator.is_odd() {
            return Err("segmented REF-0 iterator produced nonintegral u".to_owned());
        }
        let u = numerator / 2u8;
        let primitive = u.abs().gcd(&v) == BigInt::one();
        candidates.push((u, v, primitive));
    }
    let mut marks = vec![Vec::new(); params::REF0_COUNT as usize];
    let mut first = 0u64;
    while first < params::REF0_COUNT {
        let end = (first + SEGMENT_RECORDS).min(params::REF0_COUNT);
        for (entry_index, entry) in factor_base.manifest.entries.iter().enumerate() {
            let modulus = BigInt::from(entry.ell);
            for index in first..end {
                let (u, v, primitive) = &candidates[index as usize];
                if !primitive {
                    continue;
                }
                let root_hits = entry
                    .roots
                    .iter()
                    .filter(|root| (u + v * BigInt::from(**root)).mod_floor(&modulus).is_zero())
                    .count();
                if root_hits > 1 {
                    return Err(format!(
                        "segmented primitive row {index} hits both roots at {}",
                        entry.ell
                    ));
                }
                if root_hits == 1 {
                    marks[index as usize].push(entry_index);
                }
            }
        }
        first = end;
    }
    Ok(marks)
}

fn alternate_ref0(
    factor_base: &VerifiedFactorBase,
    marks: &[Vec<usize>],
    traversal: impl IntoIterator<Item = u64>,
) -> Result<AlternateReference> {
    if marks.len() as u64 != params::REF0_COUNT {
        return Err("segmented REF-0 mark table length mismatch".to_owned());
    }
    let mut rows = Vec::with_capacity(params::REF0_COUNT as usize);
    for index in traversal {
        let (v, x) = pair_at(index)?;
        let evaluation =
            evaluate_candidate_segmented(index, v, x, factor_base, &marks[index as usize])?;
        let (hashes, certificates) = match evaluation.certificates {
            Some(certificates) => {
                let hashes = Some((
                    sha256(&certificates.candidate),
                    sha256(&certificates.factorization),
                ));
                (
                    hashes,
                    Some((certificates.candidate, certificates.factorization)),
                )
            }
            None => (None, None),
        };
        rows.push((
            index,
            encode_candidate_record(v, x),
            encode_disposition_record(evaluation.status, index, hashes)?,
            certificates,
            evaluation.status,
        ));
    }
    rows.sort_by_key(|row| row.0);
    if rows.len() as u64 != params::REF0_COUNT
        || rows
            .iter()
            .enumerate()
            .any(|(index, row)| row.0 != index as u64)
    {
        return Err("alternate REF-0 traversal omitted or duplicated an index".to_owned());
    }
    let mut candidate_records = Vec::with_capacity(16_400);
    let mut disposition_records = Vec::new();
    let mut certificates = std::collections::BTreeMap::new();
    let mut status_counts = [0u64; 7];
    for (index, candidate, disposition, retained, status) in rows {
        candidate_records.extend_from_slice(&candidate);
        disposition_records.extend_from_slice(&disposition);
        status_counts[usize::from(status)] += 1;
        if let Some(retained) = retained {
            certificates.insert(index, retained);
        }
    }
    let mut certificate_store = b"P192-WCM-REF-CERT-STORE-v1\0".to_vec();
    certificate_store.extend_from_slice(
        &u64::try_from(certificates.len())
            .map_err(|_| "alternate certificate count does not fit u64".to_owned())?
            .to_be_bytes(),
    );
    for (index, (candidate, factorization)) in certificates {
        certificate_store.extend_from_slice(&index.to_be_bytes());
        certificate_store.extend_from_slice(
            &u32::try_from(candidate.len())
                .map_err(|_| "alternate CAND length does not fit u32".to_owned())?
                .to_be_bytes(),
        );
        certificate_store.extend_from_slice(&candidate);
        certificate_store.extend_from_slice(
            &u32::try_from(factorization.len())
                .map_err(|_| "alternate FACT length does not fit u32".to_owned())?
                .to_be_bytes(),
        );
        certificate_store.extend_from_slice(&factorization);
    }
    let candidate_shard = candidate_shard_digest(
        params::REF0_SHELL_ID,
        0,
        params::REF0_COUNT,
        &candidate_records,
    )?;
    let candidate_shell = candidate_shell_digest(
        params::REF0_SHELL_ID,
        &[ShardDigest {
            shard_index: 0,
            first_candidate_index: 0,
            record_count: params::REF0_COUNT,
            sha256: candidate_shard,
        }],
    )?;
    let disposition_shard = disposition_shard_digest(
        params::REF0_SHELL_ID,
        0,
        0,
        params::REF0_COUNT,
        &disposition_records,
    )?;
    let disposition_shell = disposition_shell_digest(
        params::REF0_SHELL_ID,
        &[ShardDigest {
            shard_index: 0,
            first_candidate_index: 0,
            record_count: params::REF0_COUNT,
            sha256: disposition_shard,
        }],
    )?;
    Ok(AlternateReference {
        candidate_records,
        disposition_records,
        certificate_store,
        status_counts,
        candidate_shard,
        candidate_shell,
        disposition_shard,
        disposition_shell,
    })
}

fn compare_alternate(
    reference: &RegeneratedReference,
    alternate: &AlternateReference,
) -> Result<()> {
    let reference_counts = [
        reference.status_counts.nonprimitive_duplicate,
        reference.status_counts.complete,
        reference.status_counts.one_large_prime,
        reference.status_counts.two_large_prime,
        reference.status_counts.rejected,
        reference.status_counts.invalid,
        reference.status_counts.unresolved,
    ];
    if alternate.candidate_records != reference.candidate_records
        || alternate.disposition_records != reference.disposition_records
        || alternate.certificate_store != reference.certificate_store
        || alternate.status_counts != reference_counts
        || alternate.candidate_shard != reference.candidate_shard_sha256
        || alternate.candidate_shell != reference.candidate_shell_sha256
        || alternate.disposition_shard != reference.disposition_shard_sha256
        || alternate.disposition_shell != reference.disposition_shell_sha256
    {
        return Err(
            "alternate REF-0 traversal differs in bytes, roots, certificates, or counters"
                .to_owned(),
        );
    }
    Ok(())
}

fn verify_ref0_controls(
    factor_base: &VerifiedFactorBase,
    reference: &RegeneratedReference,
) -> Result<()> {
    let marks = segmented_factor_marks(factor_base)?;
    let mut alternate = Vec::with_capacity(reference.candidate_records.len());
    for worker in 0..8u64 {
        let start = REF0_PARTITION * worker;
        let end = (REF0_PARTITION * (worker + 1)).min(params::REF0_COUNT);
        for index in start..end {
            let (v, x) = pair_at(index)?;
            alternate.extend_from_slice(&crate::encoding::encode_candidate_record(v, x));
        }
    }
    if alternate != reference.candidate_records {
        return Err("one/eight-worker candidate stream mismatch".to_owned());
    }
    let mut reversed: Vec<_> = (0..params::REF0_COUNT)
        .rev()
        .map(|index| pair_at(index).map(|pair| (index, pair)))
        .collect::<Result<_>>()?;
    reversed.sort_by_key(|entry| entry.0);
    let reverse_bytes: Vec<_> = reversed
        .into_iter()
        .flat_map(|(_, (v, x))| crate::encoding::encode_candidate_record(v, x))
        .collect();
    if reverse_bytes != reference.candidate_records {
        return Err("forward/reverse candidate stream mismatch".to_owned());
    }

    let reverse = alternate_ref0(factor_base, &marks, (0..params::REF0_COUNT).rev())?;
    compare_alternate(reference, &reverse)?;
    let worker_order = (0..8u64).flat_map(|worker| (worker..params::REF0_COUNT).step_by(8usize));
    let workers = alternate_ref0(factor_base, &marks, worker_order)?;
    compare_alternate(reference, &workers)?;

    // Role 1 has no mid-shard resume.  Two complete in-memory regenerations
    // must reproduce all three streams, four roots, and seven counters.
    let rebuild_one = alternate_ref0(factor_base, &marks, 0..params::REF0_COUNT)?;
    let rebuild_two = alternate_ref0(factor_base, &marks, 0..params::REF0_COUNT)?;
    compare_alternate(reference, &rebuild_one)?;
    compare_alternate(reference, &rebuild_two)?;
    if rebuild_one.candidate_records != rebuild_two.candidate_records
        || rebuild_one.disposition_records != rebuild_two.disposition_records
        || rebuild_one.certificate_store != rebuild_two.certificate_store
        || rebuild_one.status_counts != rebuild_two.status_counts
        || rebuild_one.candidate_shard != rebuild_two.candidate_shard
        || rebuild_one.candidate_shell != rebuild_two.candidate_shell
        || rebuild_one.disposition_shard != rebuild_two.disposition_shard
        || rebuild_one.disposition_shell != rebuild_two.disposition_shell
    {
        return Err("two fresh REF-0 rebuilds differ".to_owned());
    }
    Ok(())
}

const REF0_PARTITION: u64 = params::REF0_COUNT.div_ceil(8);

fn verify_lp_boundary() -> Result<()> {
    let limit = LARGE_PRIME_BOUND;
    if is_prime_u64(limit - 1) || limit - 1 > LARGE_PRIME_BOUND {
        return Err("q=L-1 endpoint was not rejected as composite".to_owned());
    }
    validate_lp(&lp_fixture(1))?;
    let squared = u128::from(limit) * u128::from(limit);
    let immediately_above = limit
        .checked_add(1)
        .ok_or_else(|| "LP upper endpoint overflow".to_owned())?;
    if squared != 4_611_686_014_132_420_609u128
        || squared > u128::from(u64::MAX)
        || immediately_above != 2_147_483_648
        || immediately_above <= LARGE_PRIME_BOUND
    {
        return Err("LP boundary arithmetic overflow".to_owned());
    }
    if factor_u64(limit)? != vec![limit]
        || factor_u64(u64::try_from(squared).map_err(|_| "L^2 conversion".to_owned())?)?
            != vec![limit, limit]
        || squared
            .checked_add(1)
            .ok_or_else(|| "L^2+1 overflow".to_owned())?
            != 4_611_686_014_132_420_610u128
    {
        return Err("LP residual boundary classification failed".to_owned());
    }
    Ok(())
}

pub fn verify_controls_bytes(
    bytes: &[u8],
    protocol_commit: &str,
    source_commit: &str,
    factor_base: &VerifiedFactorBase,
    reference: &RegeneratedReference,
) -> Result<VerifiedControls> {
    let value = parse_canonical_json(bytes)?;
    let manifest: ControlsManifest =
        serde_json::from_slice(bytes).map_err(|error| format!("controls schema: {error}"))?;
    if manifest.schema != CONTROLS_SCHEMA
        || manifest.experiment_id != "EXP-SCURVE-1a8daf"
        || manifest.protocol_version != 2
        || manifest.protocol_commit != protocol_commit
        || manifest.source_commit != source_commit
        || manifest.curve_uid != CURVE_UID
        || manifest.factor_base_sha256 != hex::encode(factor_base.sha256)
        || manifest.controls.len() != EXPECTED.len()
        || manifest.overall_status != "PASS"
    {
        return Err("controls envelope differs from the frozen contract".to_owned());
    }
    for (record, (control_id, assertions)) in manifest.controls.iter().zip(EXPECTED) {
        if record.control_id != *control_id
            || record.status != "PASS"
            || record.assertions.len() != assertions.len()
        {
            return Err(format!("control record {control_id} differs"));
        }
        for (actual, expected) in record.assertions.iter().zip(*assertions) {
            if actual.assertion_id != *expected || !actual.passed {
                return Err(format!("control assertion {control_id}/{expected} failed"));
            }
        }
    }

    verify_d23()?;
    verify_sqrt_d(factor_base)?;
    verify_inert_rejection(factor_base)?;
    verify_mutation_rejection()?;
    verify_lp_cancellation()?;
    verify_ref0_controls(factor_base, reference)?;
    verify_lp_boundary()?;

    let controls_value = value
        .get("controls")
        .ok_or_else(|| "controls array is absent".to_owned())?;
    let mut projection = RESULTS_DOMAIN.to_vec();
    projection.extend_from_slice(&canonical_json(controls_value)?);
    let digest = sha256(&projection);
    if manifest.control_results_sha256 != hex::encode(digest) {
        return Err("control_results_sha256 mismatch".to_owned());
    }
    Ok(VerifiedControls {
        manifest,
        canonical_bytes: bytes.to_vec(),
        results_sha256: digest,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::factor_base::{regenerate_entries, FactorBaseManifest};
    use crate::params::{ALGEBRAIC_BOUND, D_DEC, MAPPABLE_BOUND, P_DEC, T_DEC};

    #[test]
    fn d23_and_lp_boundary_replay() {
        verify_d23().unwrap();
        verify_mutation_rejection().unwrap();
        let evidence = replay_mutation_rejections().unwrap();
        assert_eq!(evidence.cases, 9);
        assert_eq!(
            evidence.mutation_identities.as_slice(),
            EXPECTED_MUTATION_IDENTITIES
        );
        assert_eq!(evidence.exact_single_scalar_changes, 9);
        assert_eq!(evidence.stale_record_hash_rejections, 9);
        assert_eq!(evidence.stale_file_hash_rejections, 9);
        assert_eq!(evidence.recomputed_hash_chain_acceptances, 9);
        assert_eq!(evidence.semantic_rejections_after_hash_acceptance, 9);
        verify_lp_cancellation().unwrap();
        verify_lp_boundary().unwrap();
        let mut wrong_conjugate = lp_fixture(2);
        wrong_conjugate.larger_root -= 1;
        assert!(validate_lp_cancellation(&[lp_fixture(1), wrong_conjugate], true).is_err());

        let mut held_fixed_d23 = d23_mutation_target();
        held_fixed_d23.alpha.v = "2".to_owned();
        assert!(
            replay_mutation_rejections_with_bases(held_fixed_d23, lp_mutation_target())
                .unwrap_err()
                .starts_with("frozen D23 mutation base fixture failed replay:")
        );
        let mut held_fixed_lp = lp_mutation_target();
        held_fixed_lp.larger_root -= 1;
        assert!(
            replay_mutation_rejections_with_bases(d23_mutation_target(), held_fixed_lp)
                .unwrap_err()
                .starts_with("frozen LP mutation base fixture failed replay:")
        );
    }

    #[test]
    fn p192_ramified_and_inert_controls_replay() {
        let entries = regenerate_entries().unwrap();
        let factor_base = VerifiedFactorBase {
            manifest: FactorBaseManifest {
                schema: crate::factor_base::FACTOR_BASE_SCHEMA.to_owned(),
                curve_uid: CURVE_UID.to_owned(),
                p: P_DEC.to_owned(),
                t: T_DEC.to_owned(),
                discriminant: D_DEC.to_owned(),
                algebraic_bound: ALGEBRAIC_BOUND,
                mappable_bound: MAPPABLE_BOUND,
                source_commit: "0".repeat(40),
                entry_count: entries.len() as u64,
                entries,
            },
            canonical_bytes: Vec::new(),
            sha256: [0u8; 32],
        };
        verify_sqrt_d(&factor_base).unwrap();
        let c = params::bigint(C_DEC, "C").unwrap();
        let u = -params::t();
        let v = BigInt::from(2u8);
        let action =
            derive_sqrt_d_subgroup_action(&u, &v, &params::subgroup_order(), &c, false, false)
                .unwrap();
        assert_eq!(action.scalar_action.to_string(), SQRT_D_SUBGROUP_SCALAR_DEC);
        assert!(!action.action_admitted);
        let mut wrong_u = u.clone();
        wrong_u += 1u8;
        assert!(derive_sqrt_d_subgroup_action(
            &wrong_u,
            &v,
            &params::subgroup_order(),
            &c,
            false,
            false,
        )
        .is_err());
        let wrong_n = params::subgroup_order() - 1u8;
        assert!(derive_sqrt_d_subgroup_action(&u, &v, &wrong_n, &c, false, false).is_err());
        assert!(
            derive_sqrt_d_subgroup_action(&u, &v, &params::subgroup_order(), &c, false, true,)
                .is_err()
        );
        verify_inert_rejection(&factor_base).unwrap();
        let marks = segmented_factor_marks(&factor_base).unwrap();
        let marked_index = marks
            .iter()
            .position(|candidate_marks| !candidate_marks.is_empty())
            .unwrap();
        let (v, x) = pair_at(marked_index as u64).unwrap();
        evaluate_candidate_segmented(
            marked_index as u64,
            v,
            x,
            &factor_base,
            &marks[marked_index],
        )
        .unwrap();
        let mut missed = marks[marked_index].clone();
        missed.remove(0);
        assert!(
            evaluate_candidate_segmented(marked_index as u64, v, x, &factor_base, &missed,)
                .is_err()
        );
        let false_mark = (0..factor_base.manifest.entries.len())
            .find(|entry| marks[marked_index].binary_search(entry).is_err())
            .unwrap();
        let mut falsely_marked = marks[marked_index].clone();
        let insertion = falsely_marked.binary_search(&false_mark).unwrap_err();
        falsely_marked.insert(insertion, false_mark);
        assert!(evaluate_candidate_segmented(
            marked_index as u64,
            v,
            x,
            &factor_base,
            &falsely_marked,
        )
        .is_err());
        let root = factor_base
            .manifest
            .entries
            .iter()
            .find(|entry| entry.ell == 5)
            .unwrap()
            .roots[0];
        let base = RamifiedSquareFixture {
            ell: 5,
            root,
            alpha_u: BigInt::from(5),
            alpha_v: BigInt::zero(),
            claimed_norm: BigInt::from(25),
            claimed_non_scalar: false,
        };
        validate_ramified_square(&base, &factor_base).unwrap();
        let mut nonscalar = base.clone();
        nonscalar.alpha_v = BigInt::one();
        assert!(validate_ramified_square(&nonscalar, &factor_base).is_err());
        let mut wrong_norm = base.clone();
        wrong_norm.claimed_norm += 1;
        assert!(validate_ramified_square(&wrong_norm, &factor_base).is_err());
        let mut mislabeled = base;
        mislabeled.claimed_non_scalar = true;
        assert!(validate_ramified_square(&mislabeled, &factor_base).is_err());
        let reference = crate::ref0::regenerate(&factor_base).unwrap();
        verify_ref0_controls(&factor_base, &reference).unwrap();
    }

    #[test]
    fn control_assertion_contract_is_exact() {
        assert_eq!(EXPECTED.len(), 9);
        assert_eq!(
            EXPECTED[7].1.last().copied(),
            Some("independent_in_memory_ref0_regenerations_streams_roots_counters_equal")
        );
        assert_eq!(
            &EXPECTED[8].1[..3],
            [
                "q_L_minus_1_not_prime_rejected_as_lp",
                "q_L_prime_split_accepted_as_one_lp",
                "q_L_plus_1_out_of_range_rejected_as_lp",
            ]
        );
    }
}
