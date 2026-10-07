//! Frozen producer-side Role-1 controls.

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde_json::{json, Value};

use super::arithmetic::{
    d23_control, hnf_from_vectors, integer, multiply_hnf, norm, p192_parameters,
    p192_sqrt_discriminant_control, pow_hnf, principal_hnf, Hnf,
};
use super::certificate::{encoded_sint, residual_prime_admissible};
use super::digest::sha256_hex;
use super::factor_base::{characteristic_roots, is_characteristic_root, Entry};
use super::primality;
use super::provenance::BuildProvenance;
use super::ref0::{
    build_reference_segmented, build_reference_with_order, ReferenceOutput, RECORD_COUNT,
};
use super::schema::{canonical_json_bytes, CURVE_UID, EXPERIMENT_ID, PROTOCOL_VERSION};
use super::Result;

fn assertion(assertion_id: &str, passed: bool) -> Value {
    json!({"assertion_id":assertion_id,"passed":passed})
}

fn control(control_id: &str, values: &[(&str, bool)]) -> Value {
    let assertions = values
        .iter()
        .map(|(assertion_id, passed)| assertion(assertion_id, *passed))
        .collect::<Vec<_>>();
    json!({
        "assertions": assertions,
        "control_id": control_id,
        "status": if values.iter().all(|(_, passed)| *passed) { "PASS" } else { "FAIL" },
    })
}

fn split_dispositions(bytes: &[u8]) -> Result<Vec<&[u8]>> {
    let mut output = Vec::new();
    let mut offset = 0usize;
    for expected_index in 0..RECORD_COUNT {
        if bytes.len().saturating_sub(offset) < 10 {
            return Err("truncated REF-0 disposition control input".to_owned());
        }
        let status = bytes[offset];
        let index = u64::from_be_bytes(
            bytes[offset + 1..offset + 9]
                .try_into()
                .map_err(|_| "disposition index width".to_owned())?,
        );
        let flag = bytes[offset + 9];
        if index != expected_index || status > 6 || flag > 1 {
            return Err("noncanonical REF-0 disposition control input".to_owned());
        }
        let retained = matches!(status, 1..=3);
        if retained != (flag == 1) {
            return Err("REF-0 disposition certificate flag mismatch".to_owned());
        }
        let length = if retained { 74 } else { 10 };
        if bytes.len().saturating_sub(offset) < length {
            return Err("truncated REF-0 disposition hashes".to_owned());
        }
        output.push(&bytes[offset..offset + length]);
        offset += length;
    }
    if offset != bytes.len() {
        return Err("trailing REF-0 disposition bytes".to_owned());
    }
    Ok(output)
}

fn certificate_store_count(bytes: &[u8]) -> Result<u64> {
    const PREFIX: &[u8] = b"P192-WCM-REF-CERT-STORE-v1\0";
    if bytes.len() < PREFIX.len() + 8 || &bytes[..PREFIX.len()] != PREFIX {
        return Err("invalid REF-0 certificate-store prefix".to_owned());
    }
    Ok(u64::from_be_bytes(
        bytes[PREFIX.len()..PREFIX.len() + 8]
            .try_into()
            .map_err(|_| "certificate-store count width".to_owned())?,
    ))
}

fn same_reference(left: &ReferenceOutput, right: &ReferenceOutput) -> bool {
    left.candidate_records == right.candidate_records
        && left.candidate_shard_sha256 == right.candidate_shard_sha256
        && left.candidate_shell_sha256 == right.candidate_shell_sha256
        && left.disposition_records == right.disposition_records
        && left.disposition_shard_sha256 == right.disposition_shard_sha256
        && left.disposition_shell_sha256 == right.disposition_shell_sha256
        && left.certificates == right.certificates
        && left.status_counts == right.status_counts
}

fn endpoint_allowed(
    candidate: u64,
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> bool {
    candidate <= bound
        && primality::is_prime(candidate)
        && residual_prime_admissible(candidate, algebraic_bound, trace, field_norm, discriminant)
}

fn residual_status(
    residual: &BigUint,
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> Result<u8> {
    let bound_squared = BigUint::from(bound) * bound;
    if residual > &bound_squared {
        return Ok(4);
    }
    let value = residual
        .to_u64()
        .ok_or_else(|| "bounded LP control residual does not fit u64".to_owned())?;
    let factors = primality::factor(value)?;
    let multiplicity = factors.iter().map(|(_, exponent)| *exponent).sum::<u32>();
    if factors.iter().any(|(prime, _)| {
        !endpoint_allowed(
            *prime,
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )
    }) {
        return Ok(4);
    }
    Ok(match multiplicity {
        1 => 2,
        2 => 3,
        _ => 4,
    })
}

#[derive(Clone)]
struct D23Fixture {
    selected_root: u64,
    signed_coordinate: BigInt,
    alpha_u: BigInt,
    alpha_v: BigInt,
    claimed_norm: BigInt,
    rational_factor_exponent: u32,
    ideal_power_top_right: BigInt,
    signed_coordinate_encoding: Vec<u8>,
}

fn d23_fixture() -> D23Fixture {
    D23Fixture {
        selected_root: 1,
        signed_coordinate: BigInt::from(-3),
        alpha_u: BigInt::one(),
        alpha_v: BigInt::one(),
        claimed_norm: BigInt::from(8u8),
        rational_factor_exponent: 3,
        ideal_power_top_right: BigInt::from(-7),
        signed_coordinate_encoding: vec![0x02, 0, 0, 0, 1, 3],
    }
}

fn validate_d23_fixture(fixture: &D23Fixture) -> Result<bool> {
    let p = BigInt::from(6u8);
    let t = BigInt::one();
    if !is_characteristic_root(fixture.selected_root, 2, &t, &p)
        || (&fixture.alpha_u + &fixture.alpha_v * BigInt::from(fixture.selected_root))
            .mod_floor(&BigInt::from(2u8))
            != BigInt::zero()
    {
        return Ok(false);
    }
    let expected_coordinate = if fixture.selected_root == 0 {
        BigInt::from(fixture.rational_factor_exponent)
    } else {
        -BigInt::from(fixture.rational_factor_exponent)
    };
    if fixture.signed_coordinate != expected_coordinate
        || encoded_sint(&fixture.signed_coordinate)? != fixture.signed_coordinate_encoding
        || norm(&fixture.alpha_u, &fixture.alpha_v, &t, &p) != fixture.claimed_norm
        || BigInt::from(2u8).pow(fixture.rational_factor_exponent) != fixture.claimed_norm
    {
        return Ok(false);
    }
    let ideal_power = pow_hnf(
        &Hnf::prime(2, fixture.selected_root)?,
        fixture.rational_factor_exponent,
        &p,
        &t,
    )?;
    let recorded_hnf = hnf_from_vectors(&[
        (BigInt::from(8u8), BigInt::zero()),
        (fixture.ideal_power_top_right.clone(), BigInt::one()),
    ])?;
    Ok(ideal_power == recorded_hnf
        && ideal_power == principal_hnf(&fixture.alpha_u, &fixture.alpha_v, &p, &t)?)
}

#[derive(Clone)]
struct LpFixture {
    q: u64,
    smaller_root: u64,
    larger_root: u64,
    sign: u8,
    multiplicity: u32,
}

fn validate_lp_fixture(
    fixture: &LpFixture,
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> bool {
    fixture.smaller_root < fixture.larger_root
        && fixture.sign == 1
        && fixture.multiplicity == 1
        && is_characteristic_root(fixture.smaller_root, fixture.q, trace, field_norm)
        && is_characteristic_root(fixture.larger_root, fixture.q, trace, field_norm)
        && endpoint_allowed(
            fixture.q,
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )
}

enum MutationFixture {
    D23(D23Fixture),
    Lp(LpFixture),
}

fn validate_mutation_fixture(
    fixture: &MutationFixture,
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> Result<bool> {
    match fixture {
        MutationFixture::D23(fixture) => validate_d23_fixture(fixture),
        MutationFixture::Lp(fixture) => Ok(validate_lp_fixture(
            fixture,
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )),
    }
}

fn mutation_rejected(
    fixture: MutationFixture,
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> Result<bool> {
    Ok(!validate_mutation_fixture(
        &fixture,
        bound,
        algebraic_bound,
        trace,
        field_norm,
        discriminant,
    )?)
}

#[derive(Debug, PartialEq, Eq)]
struct FrozenMutationResults {
    root: bool,
    orientation_sign: bool,
    exponent: bool,
    alpha_u: bool,
    norm: bool,
    rational_factor_exponent: bool,
    lp_smaller_root: bool,
    hnf_matrix_entry: bool,
    sint_sign_byte: bool,
}

impl FrozenMutationResults {
    fn all_rejected(&self) -> bool {
        self.root
            && self.orientation_sign
            && self.exponent
            && self.alpha_u
            && self.norm
            && self.rational_factor_exponent
            && self.lp_smaller_root
            && self.hnf_matrix_entry
            && self.sint_sign_byte
    }
}

fn run_frozen_mutations(
    bound: u64,
    algebraic_bound: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> Result<FrozenMutationResults> {
    let d23 = d23_fixture();
    if !validate_mutation_fixture(
        &MutationFixture::D23(d23.clone()),
        bound,
        algebraic_bound,
        trace,
        field_norm,
        discriminant,
    )? {
        return Err("frozen D23 mutation base fixture failed replay".to_owned());
    }
    let lp = LpFixture {
        q: bound,
        smaller_root: 1_109_020_142,
        larger_root: 2_025_854_371,
        sign: 1,
        multiplicity: 1,
    };
    if !validate_mutation_fixture(
        &MutationFixture::Lp(lp.clone()),
        bound,
        algebraic_bound,
        trace,
        field_norm,
        discriminant,
    )? {
        return Err("frozen LP mutation base fixture failed replay".to_owned());
    }

    let mut root = d23.clone();
    root.selected_root = 0;
    let mut orientation_sign = d23.clone();
    orientation_sign.signed_coordinate = BigInt::from(3u8);
    let mut exponent = d23.clone();
    exponent.signed_coordinate = BigInt::from(-2);
    let mut alpha_u = d23.clone();
    alpha_u.alpha_u = BigInt::from(2u8);
    let mut norm = d23.clone();
    norm.claimed_norm = BigInt::from(4u8);
    let mut rational_factor_exponent = d23.clone();
    rational_factor_exponent.rational_factor_exponent = 2;
    let mut hnf_matrix_entry = d23.clone();
    hnf_matrix_entry.ideal_power_top_right = BigInt::from(-6);
    let mut sint_sign_byte = d23;
    sint_sign_byte.signed_coordinate_encoding = vec![0x01, 0, 0, 0, 1, 3];
    let mut lp_smaller_root = lp;
    lp_smaller_root.smaller_root = 1_109_020_143;

    Ok(FrozenMutationResults {
        root: mutation_rejected(
            MutationFixture::D23(root),
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?,
        orientation_sign: mutation_rejected(
            MutationFixture::D23(orientation_sign),
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?,
        exponent: mutation_rejected(
            MutationFixture::D23(exponent),
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?,
        alpha_u: mutation_rejected(
            MutationFixture::D23(alpha_u),
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?,
        norm: mutation_rejected(
            MutationFixture::D23(norm),
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?,
        rational_factor_exponent: mutation_rejected(
            MutationFixture::D23(rational_factor_exponent),
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?,
        lp_smaller_root: mutation_rejected(
            MutationFixture::Lp(lp_smaller_root),
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?,
        hnf_matrix_entry: mutation_rejected(
            MutationFixture::D23(hnf_matrix_entry),
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?,
        sint_sign_byte: mutation_rejected(
            MutationFixture::D23(sint_sign_byte),
            bound,
            algebraic_bound,
            trace,
            field_norm,
            discriminant,
        )?,
    })
}

pub fn build(
    provenance: &BuildProvenance,
    factor_base_sha256: &str,
    factor_base_digest: &[u8; 32],
    factor_base: &[Entry],
    reference: &ReferenceOutput,
) -> Result<Value> {
    let d23 = d23_control()?;
    let d23_fixture = d23_fixture();
    let d23_fixture_valid = validate_d23_fixture(&d23_fixture)?;
    let d23_principal = d23.prime_ideal_cubed == d23.principal_alpha;
    let form_a = BigInt::from(2u8);
    let form_b = BigInt::one();
    let form_c = BigInt::from(3u8);
    let form_discriminant = &form_b * &form_b - BigInt::from(4u8) * &form_a * &form_c;
    let canonical_form = form_discriminant == BigInt::from(-23)
        && form_b.abs() <= form_a
        && form_a <= form_c
        && form_a.gcd(&form_b).gcd(&form_c).is_one();
    let d23_ideal_matches = is_characteristic_root(1, 2, &BigInt::one(), &BigInt::from(6u8))
        && Hnf::prime(2, 1)?
            == Hnf {
                a: BigInt::from(2u8),
                b: BigInt::one(),
                d: BigInt::one(),
            };
    let d23_x = BigInt::from(2u8) * &d23_fixture.alpha_u + &d23_fixture.alpha_v;
    let d23_alpha_valid =
        d23_fixture_valid && d23.norm == d23_fixture.claimed_norm && d23_x == BigInt::from(3u8);
    let d23_values = [
        ("canonical_reduced_form_norm_2", canonical_form),
        ("ideal_matches_form", d23_ideal_matches),
        ("alpha_certificate_valid", d23_alpha_valid),
        (
            "encoding_replays",
            encoded_sint(&d23_fixture.signed_coordinate)? == d23_fixture.signed_coordinate_encoding,
        ),
        ("relation_g2_cubed_identity", d23_principal),
    ];

    let sqrt_d = p192_sqrt_discriminant_control()?;
    let (p, t, d) = p192_parameters()?;
    let c = integer(super::arithmetic::C_DEC)?;
    let inverse_two_mod_c = (&c + BigInt::one()) / BigInt::from(2u8);
    let root_c = (&t * inverse_two_mod_c).mod_floor(&c);
    let prime_c = Hnf {
        a: c.clone(),
        b: (-&root_c).mod_floor(&c),
        d: BigInt::one(),
    };
    let mut sqrt_ideal = Hnf::identity();
    for (ell, root) in [(5u64, 2u64), (11, 3), (31, 7)] {
        sqrt_ideal = multiply_hnf(&sqrt_ideal, &Hnf::prime(ell, root)?, &p, &t)?;
    }
    sqrt_ideal = multiply_hnf(&sqrt_ideal, &prime_c, &p, &t)?;
    let sqrt_principal = principal_hnf(&sqrt_d.u, &sqrt_d.v, &p, &t)?;
    let p192_algebra_passes =
        sqrt_ideal == sqrt_principal && sqrt_ideal.determinant() == sqrt_d.norm;
    let c_degree_unavailable =
        c > BigInt::from(u64::MAX) && !factor_base.iter().any(|entry| BigInt::from(entry.ell) == c);
    let subgroup_order = integer(super::arithmetic::N_DEC)?;
    let scalar_action = (&sqrt_d.u + &sqrt_d.v).mod_floor(&subgroup_order);
    let scalar_action_rejected = scalar_action.to_string()
        == "6277101735386680763835789423144451611450480845975359606884"
        && !scalar_action.is_zero()
        && sqrt_d.residual_exceeds_bound_squared;
    let p192_values = [
        (
            "alpha_coordinates_minus_t_2",
            sqrt_d.u == -integer(super::arithmetic::T_DEC)? && sqrt_d.v == BigInt::from(2u8),
        ),
        (
            "norm_equals_abs_D",
            sqrt_d.norm == integer(super::arithmetic::D_DEC)?.abs(),
        ),
        (
            "factorization_equals_5_11_31_C",
            sqrt_d.factorization.factors == vec![(5, 1), (11, 1), (31, 1)]
                && sqrt_d.factorization.residual.to_string() == super::arithmetic::C_DEC,
        ),
        ("algebra_passes", p192_algebra_passes),
        ("C_degree_map_unavailable", c_degree_unavailable),
        ("C_scalar_action_rejected", scalar_action_rejected),
    ];

    let mut ramified = Vec::new();
    for (ell, root, id) in [
        (5u64, 2u64, "ell_5_scalar"),
        (11, 3, "ell_11_scalar"),
        (31, 7, "ell_31_scalar"),
    ] {
        let square = pow_hnf(&Hnf::prime(ell, root)?, 2, &p, &t)?;
        let scalar = Hnf {
            a: BigInt::from(ell),
            b: BigInt::zero(),
            d: BigInt::from(ell),
        };
        ramified.push((id, square == scalar));
    }

    let inert_absent = !factor_base.iter().any(|entry| entry.ell == 2)
        && characteristic_roots(2, &t, &p, &d).is_empty();
    let inert_residual_prime = 65_537u64;
    let inert_residual_rejected = inert_residual_prime > 65_521
        && primality::is_prime(inert_residual_prime)
        && !residual_prime_admissible(inert_residual_prime, 65_521, &t, &p, &d)
        && characteristic_roots(inert_residual_prime, &t, &p, &d).is_empty();
    let inert_values = [
        ("base_injection_rejected", inert_absent),
        ("residual_injection_rejected", inert_residual_rejected),
    ];

    let large = 2_147_483_647u64;
    let mutations = run_frozen_mutations(large, 65_521, &t, &p, &d)?;
    if !mutations.all_rejected() {
        return Err("one or more frozen Role-1 mutations were accepted".to_owned());
    }
    let mutation_values = [
        ("root_mutation_rejected", mutations.root),
        (
            "orientation_sign_mutation_rejected",
            mutations.orientation_sign && mutations.sint_sign_byte,
        ),
        ("exponent_mutation_rejected", mutations.exponent),
        ("alpha_coordinate_mutation_rejected", mutations.alpha_u),
        (
            "norm_factor_mutation_rejected",
            mutations.norm && mutations.rational_factor_exponent,
        ),
        ("lp_endpoint_mutation_rejected", mutations.lp_smaller_root),
        ("hnf_column_mutation_rejected", mutations.hnf_matrix_entry),
    ];
    let lp_roots = characteristic_roots(large, &t, &p, &d);
    if lp_roots.len() != 2 {
        return Err("LP cancellation fixture is not split".to_owned());
    }
    let lp_positive = Hnf::prime(large, lp_roots[0])?;
    let lp_negative = Hnf::prime(large, lp_roots[1])?;
    let lp_scalar = Hnf {
        a: BigInt::from(large),
        b: BigInt::zero(),
        d: BigInt::from(large),
    };
    let cancellation_values = [
        ("unmatched_orientation_rejected", lp_positive != lp_scalar),
        (
            "exact_conjugate_cancels",
            multiply_hnf(&lp_positive, &lp_negative, &p, &t)? == lp_scalar,
        ),
    ];

    let canonical_order = (0..RECORD_COUNT).collect::<Vec<_>>();
    let rebuilt = build_reference_with_order(
        factor_base,
        factor_base_digest,
        65_521,
        large,
        &canonical_order,
    )?;
    let segmented = build_reference_segmented(factor_base, factor_base_digest, 65_521, large)?;
    let disposition_records = split_dispositions(&segmented.disposition_records)?;
    let retained_record_count = disposition_records
        .iter()
        .filter(|record| matches!(record[0], 1..=3))
        .count() as u64;
    let scalar_segmented_values = [
        (
            "candidate_bytes_equal",
            segmented.candidate_records == reference.candidate_records,
        ),
        (
            "candidate_shard_root_equal",
            segmented.candidate_shard_sha256 == reference.candidate_shard_sha256,
        ),
        (
            "candidate_shell_root_equal",
            segmented.candidate_shell_sha256 == reference.candidate_shell_sha256,
        ),
        (
            "disposition_bytes_equal",
            segmented.disposition_records == reference.disposition_records,
        ),
        (
            "disposition_shard_root_equal",
            segmented.disposition_shard_sha256 == reference.disposition_shard_sha256,
        ),
        (
            "disposition_shell_root_equal",
            segmented.disposition_shell_sha256 == reference.disposition_shell_sha256,
        ),
        (
            "retained_records_equal",
            retained_record_count == reference.status_counts.retained()
                && segmented.certificates == reference.certificates
                && segmented.status_counts == reference.status_counts
                && certificate_store_count(&segmented.certificates)?
                    == reference.status_counts.retained(),
        ),
    ];
    let reverse_order = (0..RECORD_COUNT).rev().collect::<Vec<_>>();
    let reverse = build_reference_with_order(
        factor_base,
        factor_base_digest,
        65_521,
        large,
        &reverse_order,
    )?;
    let mut worker_order = Vec::with_capacity(RECORD_COUNT as usize);
    for lane in 0..8u64 {
        worker_order.extend((lane..RECORD_COUNT).step_by(8));
    }
    let workers = build_reference_with_order(
        factor_base,
        factor_base_digest,
        65_521,
        large,
        &worker_order,
    )?;
    let direction_values = [
        (
            "forward_reverse_candidate_root_equal",
            reverse.candidate_shard_sha256 == reference.candidate_shard_sha256
                && reverse.candidate_shell_sha256 == reference.candidate_shell_sha256,
        ),
        (
            "forward_reverse_disposition_root_equal",
            reverse.disposition_shard_sha256 == reference.disposition_shard_sha256
                && reverse.disposition_shell_sha256 == reference.disposition_shell_sha256,
        ),
        (
            "one_eight_worker_candidate_root_equal",
            workers.candidate_shard_sha256 == reference.candidate_shard_sha256
                && workers.candidate_shell_sha256 == reference.candidate_shell_sha256,
        ),
        (
            "one_eight_worker_disposition_root_equal",
            workers.disposition_shard_sha256 == reference.disposition_shard_sha256
                && workers.disposition_shell_sha256 == reference.disposition_shell_sha256,
        ),
        (
            "independent_in_memory_ref0_regenerations_streams_roots_counters_equal",
            same_reference(&rebuilt, reference),
        ),
    ];

    let large_roots = characteristic_roots(large, &t, &p, &d);
    let large_minus_one = large - 1;
    let large_plus_one = large + 1;
    let large_squared = BigUint::from(large) * large;
    let above_large_squared = &large_squared + BigUint::one();
    let lp_values = [
        (
            "q_L_minus_1_not_prime_rejected_as_lp",
            !endpoint_allowed(large_minus_one, large, 65_521, &t, &p, &d),
        ),
        (
            "q_L_prime_split_accepted_as_one_lp",
            endpoint_allowed(large, large, 65_521, &t, &p, &d)
                && large_roots == vec![1_109_020_142, 2_025_854_371],
        ),
        (
            "q_L_plus_1_out_of_range_rejected_as_lp",
            !endpoint_allowed(large_plus_one, large, 65_521, &t, &p, &d),
        ),
        (
            "residual_L_is_one_large_prime",
            residual_status(&BigUint::from(large), large, 65_521, &t, &p, &d)? == 2,
        ),
        (
            "residual_L_squared_is_two_large_prime_with_multiplicity_two",
            residual_status(&large_squared, large, 65_521, &t, &p, &d)? == 3
                && primality::factor(large * large)? == vec![(large, 2)],
        ),
        (
            "residual_L_squared_plus_one_is_rejected_by_inequality",
            residual_status(&above_large_squared, large, 65_521, &t, &p, &d)? == 4,
        ),
    ];

    let controls = vec![
        control("D23-ORDER3", &d23_values),
        control("P192-SQRT-D", &p192_values),
        control("P192-RAMIFIED-SQUARES", &ramified),
        control("INERT-FACTOR-REJECTION", &inert_values),
        control("MUTATION-REJECTION", &mutation_values),
        control("LP-CANCELLATION", &cancellation_values),
        control("REF0-SCALAR-SEGMENTED", &scalar_segmented_values),
        control("REF0-DIRECTION-WORKERS", &direction_values),
        control("LP-BOUNDARY", &lp_values),
    ];
    if controls.iter().any(|value| value["status"] != "PASS") {
        return Err("one or more frozen Role-1 controls failed".to_owned());
    }
    let mut digest_preimage = b"P192-WCM-CONTROL-RESULTS-v1\0".to_vec();
    digest_preimage.extend_from_slice(&canonical_json_bytes(&Value::Array(controls.clone()))?);
    let control_results_sha256 = sha256_hex(&digest_preimage);
    Ok(json!({
        "control_results_sha256": control_results_sha256,
        "controls": controls,
        "curve_uid": CURVE_UID,
        "experiment_id": EXPERIMENT_ID,
        "factor_base_sha256": factor_base_sha256,
        "overall_status": "PASS",
        "protocol_commit": provenance.protocol_commit,
        "protocol_version": PROTOCOL_VERSION,
        "schema": "p192-wcm-controls-v1",
        "source_commit": provenance.source_commit,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disposition_splitter_rejects_wrong_certificate_flag() {
        let mut bytes = Vec::new();
        for index in 0..RECORD_COUNT {
            bytes.push(4);
            bytes.extend_from_slice(&index.to_be_bytes());
            bytes.push(0);
        }
        assert_eq!(
            split_dispositions(&bytes).unwrap().len() as u64,
            RECORD_COUNT
        );
        bytes[9] = 1;
        assert!(split_dispositions(&bytes).is_err());
    }

    #[test]
    fn all_nine_frozen_mutations_replay_through_common_validator() {
        let (field_norm, trace, discriminant) = p192_parameters().unwrap();
        let results =
            run_frozen_mutations(2_147_483_647, 65_521, &trace, &field_norm, &discriminant)
                .unwrap();
        assert_eq!(
            [
                results.root,
                results.orientation_sign,
                results.exponent,
                results.alpha_u,
                results.norm,
                results.rational_factor_exponent,
                results.lp_smaller_root,
                results.hnf_matrix_entry,
                results.sint_sign_byte,
            ],
            [true; 9]
        );
        assert!(results.all_rejected());
    }

    #[test]
    fn complete_role1_control_set_recomputes() {
        let factor_base = super::super::factor_base::build(65_521).unwrap();
        let digest = [0x5au8; 32];
        let reference =
            super::super::ref0::build_reference(&factor_base, &digest, 65_521, 2_147_483_647)
                .unwrap();
        let provenance = BuildProvenance {
            protocol_commit: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            source_commit: "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
            source_dirty: false,
            git_metadata_present: true,
            build_profile: "release",
        };
        let document = build(
            &provenance,
            &super::super::ref0::hex(&digest),
            &digest,
            &factor_base,
            &reference,
        )
        .unwrap();
        assert_eq!(document["overall_status"], "PASS");
        assert_eq!(document["controls"].as_array().unwrap().len(), 9);
        assert_eq!(
            document["controls"][7]["assertions"][4]["assertion_id"],
            "independent_in_memory_ref0_regenerations_streams_roots_counters_equal"
        );
    }
}
