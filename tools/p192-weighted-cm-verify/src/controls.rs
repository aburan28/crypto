use num_bigint::BigInt;
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
    identity::C_DEC,
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

#[derive(Clone, Debug, Serialize)]
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
    signed_coordinate_sint: Vec<u8>,
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
        signed_coordinate_sint: hex::decode("020000000103").expect("frozen hex"),
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
    if encoded != fixture.signed_coordinate_sint {
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
    if abs_d != BigInt::from(5u8) * 11u8 * 31u8 * &c {
        return Err("P-192 sqrt(D) factorization failed".to_owned());
    }
    if c <= BigInt::from(params::MAPPABLE_BOUND)
        || factor_base
            .manifest
            .entries
            .iter()
            .any(|entry| BigInt::from(entry.ell) == c)
        || v.is_zero()
    {
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

#[derive(Clone, Debug, Serialize)]
struct LpFixture {
    q: u64,
    discriminant_residue: u64,
    smaller_root: u64,
    larger_root: u64,
    sign: u8,
    multiplicity: u64,
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

fn mutation_envelope_hash<T: Serialize>(fixture_kind: &str, fixture: &T) -> Result<[u8; 32]> {
    let value = serde_json::json!({
        "fixture": fixture,
        "fixture_kind": fixture_kind,
        "schema": "p192-wcm-control-mutation-envelope-v1",
    });
    Ok(sha256(&canonical_json(&value)?))
}

fn verify_mutation_rejection() -> Result<()> {
    validate_d23(&d23_fixture())?;
    let base = d23_fixture();
    let base_hash = mutation_envelope_hash("D23-ORDER3", &base)?;
    let mut mutations = Vec::new();
    let mut fixture = base.clone();
    fixture.selected_root = 0;
    mutations.push(("root", fixture));
    let mut fixture = base.clone();
    fixture.signed_coordinate = 3;
    mutations.push(("orientation_sign", fixture));
    let mut fixture = base.clone();
    fixture.signed_coordinate = -2;
    mutations.push(("exponent", fixture));
    let mut fixture = base.clone();
    fixture.alpha_u = 2;
    mutations.push(("alpha_u", fixture));
    let mut fixture = base.clone();
    fixture.norm = 4;
    mutations.push(("norm", fixture));
    let mut fixture = base.clone();
    fixture.rational_exponent = 2;
    mutations.push(("rational_factor_exponent", fixture));
    let mut fixture = base.clone();
    fixture.ideal_power_matrix[0][1] = -6;
    mutations.push(("hnf_matrix_entry", fixture));
    let mut fixture = base;
    fixture.signed_coordinate_sint = hex::decode("010000000103").expect("frozen hex");
    mutations.push(("sint_sign_byte", fixture));
    for (mutation, fixture) in mutations {
        let mutated_hash = mutation_envelope_hash("D23-ORDER3", &fixture)?;
        if mutated_hash == base_hash {
            return Err(format!(
                "D=-23 mutation {mutation} did not change its recomputed envelope hash"
            ));
        }
        if validate_d23(&fixture).is_ok() {
            return Err(format!("D=-23 mutation {mutation} was accepted"));
        }
    }

    let mut lp = lp_fixture(1);
    validate_lp(&lp)?;
    let lp_hash = mutation_envelope_hash("LP-BOUNDARY", &lp)?;
    lp.smaller_root = 1_109_020_143;
    if mutation_envelope_hash("LP-BOUNDARY", &lp)? == lp_hash {
        return Err("LP mutation did not change its recomputed envelope hash".to_owned());
    }
    if validate_lp(&lp).is_ok() {
        return Err("LP smaller-root endpoint mutation was accepted".to_owned());
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
        verify_lp_cancellation().unwrap();
        verify_lp_boundary().unwrap();
        let mut wrong_conjugate = lp_fixture(2);
        wrong_conjugate.larger_root -= 1;
        assert!(validate_lp_cancellation(&[lp_fixture(1), wrong_conjugate], true).is_err());
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
