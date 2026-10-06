//! Frozen native recovery protocol for the P-192 singular companion.
//!
//! The victim lane calls [`crate::ecc::ct::scalar_mul_secret`] on every
//! attacker-selected point and exposes only a domain-separated HMAC tag keyed
//! by the 24-byte x-coordinate.  All prediction, certification, and search
//! arithmetic is independent fixed-width code from [`super::p192_native`].

use super::p192_interval_bsgs::{
    make_target, solve_two_orientations, IntervalPlan, IntervalSolution, IntervalStats, FROZEN_M,
    FROZEN_STRIDE, FROZEN_TABLE_BITS,
};
use super::p192_native::{
    batch_torus_to_curve, is_on_p192, is_on_singular_companion, p192_order, p192_prime,
    residue_modulus, scalar_mul, singular_factor_product, singular_order, AffinePoint,
    NativeOpCounts, P192Arithmetic, TorusPoint, B_HEX, GX_HEX, GY_HEX, N_HEX, P_HEX, QMAX,
    SINGULAR_DISTINCT_FACTORS,
};
use crate::ecc::ct::scalar_mul_secret;
use crate::ecc::ecdh::ecdh_raw;
use crate::ecc::keys::{EccPrivateKey, EccPublicKey};
use crate::ecc::point::Point;
use crate::ecc::CurveParams;
use crate::hash::{hmac_sha256, sha256};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};

pub const EXPERIMENT_ID: &str = "EXP-SCURVE-29040c";
pub const REPORT_SCHEMA: &str = "crypto.p192-singular-recovery/v2";
pub const FROZEN_SEED_HEX: &str = "503139322d53494e47554c41522d5631";
pub const EXPECTED_DIGEST_HEX: &str =
    "f16fdbf87ede32e18b70f853c97659bb0cf20658693385892c918c91970c03af";
pub const EXPECTED_D_DEC: &str = "3419090462492416794041408622244871729251768983841006244224";
pub const EXPECTED_D_HEX: &str = "8b70f85429c81223e88be1a60590627a3c8e558a24dd5180";
pub const EXPECTED_RESIDUE_SMALL: &str = "15075375527598908552546067046058862913433216";
pub const EXPECTED_RESIDUE_LARGE: &str = "78222222987683236539393919745863586265518464";
pub const EXPECTED_RESIDUAL_K: u64 = 36_647_143_301_682;
pub const EXPECTED_CANDIDATE_COMPARISONS: u64 = 1_245_111;

const SEED_DOMAIN: &[u8] = b"crypto/p192-singular-recovery/seed/v1\0";
const KDF_DOMAIN: &[u8] = b"crypto/p192-singular/x-only-kdf/v1\0";
const CONFIRM_DOMAIN: &[u8] = b"crypto/p192-singular/confirm/v1\0";

fn decimal(value: &str) -> BigUint {
    BigUint::parse_bytes(value.as_bytes(), 10).expect("hard-coded decimal integer")
}

fn hexadecimal(value: &str) -> BigUint {
    BigUint::parse_bytes(value.as_bytes(), 16).expect("hard-coded hexadecimal integer")
}

fn gcd(mut left: BigUint, mut right: BigUint) -> BigUint {
    while !right.is_zero() {
        let remainder = &left % &right;
        left = right;
        right = remainder;
    }
    left
}

fn mul_mod_u64(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}

fn pow_mod_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut acc = 1u64;
    base %= modulus;
    while exponent != 0 {
        if exponent & 1 == 1 {
            acc = mul_mod_u64(acc, base, modulus);
        }
        base = mul_mod_u64(base, base, modulus);
        exponent >>= 1;
    }
    acc
}

/// Deterministic Miller-Rabin for every `u64`.
fn is_prime_u64(value: u64) -> bool {
    if value < 2 {
        return false;
    }
    for prime in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if value == prime {
            return true;
        }
        if value % prime == 0 {
            return false;
        }
    }
    let mut odd = value - 1;
    let powers = odd.trailing_zeros();
    odd >>= powers;
    for witness in [2u64, 325, 9_375, 28_178, 450_775, 9_780_504, 1_795_265_022] {
        if witness % value == 0 {
            continue;
        }
        let mut x = pow_mod_u64(witness % value, odd, value);
        if x == 1 || x == value - 1 {
            continue;
        }
        let mut passed = false;
        for _ in 1..powers {
            x = mul_mod_u64(x, x, value);
            if x == value - 1 {
                passed = true;
                break;
            }
        }
        if !passed {
            return false;
        }
    }
    true
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct PocklingtonFactor {
    pub prime: String,
    pub exponent: u32,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct PocklingtonCertificate {
    pub value: String,
    pub factorization_of_value_minus_one: Vec<PocklingtonFactor>,
    pub witnesses_by_prime: Vec<(String, String)>,
    pub factorization_exact: bool,
    pub verified_prime: bool,
}

fn prove_pocklington(
    value: &BigUint,
    factors: &[(BigUint, u32)],
    proven_large_primes: &[BigUint],
) -> Result<PocklingtonCertificate, String> {
    let mut product = BigUint::one();
    for (prime, exponent) in factors {
        let small_prime = prime
            .iter_u64_digits()
            .collect::<Vec<_>>()
            .as_slice()
            .first()
            .copied()
            .filter(|_| prime.bits() <= 64)
            .is_some_and(is_prime_u64);
        if !small_prime && !proven_large_primes.iter().any(|known| known == prime) {
            return Err(format!("unproven Pocklington factor {prime}"));
        }
        product *= prime.pow(*exponent);
    }
    if product != value - BigUint::one() {
        return Err(format!("incomplete factorization for {value}-1"));
    }

    let value_minus_one = value - BigUint::one();
    let mut witnesses = Vec::new();
    for (prime, _) in factors {
        let exponent = &value_minus_one / prime;
        let mut found = None;
        for candidate in 2u64..10_000 {
            let base = BigUint::from(candidate);
            if base.modpow(&value_minus_one, value) != BigUint::one() {
                continue;
            }
            let residue = base.modpow(&exponent, value);
            let difference = if residue.is_zero() {
                value - BigUint::one()
            } else {
                residue - BigUint::one()
            };
            if gcd(difference, value.clone()) == BigUint::one() {
                found = Some(base);
                break;
            }
        }
        let witness = found.ok_or_else(|| format!("no Pocklington witness for q={prime}"))?;
        witnesses.push((prime.to_string(), witness.to_string()));
    }
    Ok(PocklingtonCertificate {
        value: value.to_string(),
        factorization_of_value_minus_one: factors
            .iter()
            .map(|(prime, exponent)| PocklingtonFactor {
                prime: prime.to_string(),
                exponent: *exponent,
            })
            .collect(),
        witnesses_by_prime: witnesses,
        factorization_exact: true,
        verified_prime: true,
    })
}

fn pocklington_certificates() -> Result<Vec<PocklingtonCertificate>, String> {
    let p_child = decimal("288626509448065367648032903");
    let n_child_a = decimal("9564682313913860059195669");
    let n_child_b = decimal("3433859179316188682119986911");

    let p_child_factors = vec![
        (BigUint::from(2u8), 1),
        (BigUint::from(3u8), 1),
        (BigUint::from(19u8), 2),
        (BigUint::from(37u8), 1),
        (BigUint::from(607u16), 1),
        (decimal("5933177618131140283"), 1),
    ];
    let n_child_a_factors = vec![
        (BigUint::from(2u8), 2),
        (BigUint::from(3u8), 1),
        (BigUint::from(7u8), 1),
        (BigUint::from(15_716_741u32), 1),
        (decimal("7244839476697597"), 1),
    ];
    let n_child_b_factors = vec![
        (BigUint::from(2u8), 1),
        (BigUint::from(5u8), 1),
        (BigUint::from(439u16), 1),
        (BigUint::from(103_840_547u32), 1),
        (decimal("7532705587894727"), 1),
    ];
    let p_child_cert = prove_pocklington(&p_child, &p_child_factors, &[])?;
    let n_child_a_cert = prove_pocklington(&n_child_a, &n_child_a_factors, &[])?;
    let n_child_b_cert = prove_pocklington(&n_child_b, &n_child_b_factors, &[])?;

    let p_factors = vec![
        (BigUint::from(2u8), 1),
        (BigUint::from(59u8), 1),
        (BigUint::from(149_309u32), 1),
        (BigUint::from(11_393_611u32), 1),
        (decimal("108341181769254293"), 1),
        (p_child.clone(), 1),
    ];
    let n_factors = vec![
        (BigUint::from(2u8), 4),
        (BigUint::from(5u8), 1),
        (BigUint::from(2_389u16), 1),
        (n_child_a.clone(), 1),
        (n_child_b.clone(), 1),
    ];
    let p_cert = prove_pocklington(&p192_prime(), &p_factors, &[p_child])?;
    let n_cert = prove_pocklington(&p192_order(), &n_factors, &[n_child_a, n_child_b])?;
    Ok(vec![
        p_child_cert,
        n_child_a_cert,
        n_child_b_cert,
        p_cert,
        n_cert,
    ])
}

fn singular_factorization() -> Vec<(u64, u32)> {
    vec![
        (2, 64),
        (3, 1),
        (5, 1),
        (17, 1),
        (257, 1),
        (641, 1),
        (65_537, 1),
        (274_177, 1),
        (6_700_417, 1),
        (QMAX, 1),
    ]
}

/// The actual acceptance predicate used for a claimed factorization of the
/// singular smooth-locus order.  Negative controls mutate inputs to this
/// predicate; they do not merely compare two known-different constants.
fn accepts_singular_factorization(order: &BigUint, factors: &[(u64, u32)]) -> bool {
    let mut product = BigUint::one();
    for (prime, exponent) in factors {
        if !is_prime_u64(*prime) {
            return false;
        }
        product *= BigUint::from(*prime).pow(*exponent);
    }
    product == *order
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct PointHex {
    pub x: String,
    pub y: String,
}

fn point_hex(point: AffinePoint, arithmetic: &P192Arithmetic) -> PointHex {
    PointHex {
        x: hex::encode(arithmetic.to_bytes(point.x)),
        y: hex::encode(arithmetic.to_bytes(point.y)),
    }
}

fn point_from_hex(point: &PointHex, arithmetic: &P192Arithmetic) -> Result<AffinePoint, String> {
    Ok(AffinePoint {
        x: arithmetic.from_biguint(&hexadecimal(&point.x))?,
        y: arithmetic.from_biguint(&hexadecimal(&point.y))?,
    })
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CertificationReport {
    pub standardized_tuple_exact: bool,
    pub pocklington_certificates: Vec<PocklingtonCertificate>,
    pub p_and_n_proven_prime: bool,
    pub generator_on_p192: bool,
    pub generator_order_n: bool,
    pub singular_discriminant_zero: bool,
    pub singular_factorization_exact: bool,
    pub singular_factor_primes_exact: bool,
    pub r_finite_smooth_and_off_p192: bool,
    pub r_order_n_singular: bool,
    pub torus_parameterization_verified: bool,
    pub operations: NativeOpCounts,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ResidueStage {
    pub stage: u32,
    pub factor: u64,
    pub modulus_before: String,
    pub modulus_after: String,
    pub probe: PointHex,
    pub probe_exact_order: bool,
    pub primary_tag_hex: String,
    pub candidate_comparisons: u64,
    pub matched_lift: String,
    pub retained_negative: String,
    pub low_level_matches_independent_torus: bool,
    pub safe_ecdh_rejected: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ResidueReport {
    pub mode: String,
    pub primary_oracle_queries: u64,
    pub candidate_comparisons: u64,
    pub residue: String,
    pub negative_residue: String,
    pub expected_unordered_residues_match: bool,
    pub transcript_sha256: String,
    pub stages: Vec<ResidueStage>,
    pub operations: NativeOpCounts,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct IntervalCertificate {
    pub geometry: String,
    pub shifted_centers: bool,
    pub plan: IntervalPlan,
    pub plus_residue: String,
    pub minus_residue: String,
    pub plus_width: u64,
    pub minus_width: u64,
    pub solution: IntervalSolution,
    pub stats: IntervalStats,
    pub recovered_d: String,
    pub recovered_d_hex: String,
    pub final_replay_matches_public_key: bool,
    pub final_replay_operations: NativeOpCounts,
    pub frozen_tail_covered: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ControlReport {
    pub exhaustive_p59_sign_classes: bool,
    pub p59_scalars_tested: u64,
    pub p59_nested_queries: u64,
    pub p59_candidate_comparisons: u64,
    pub p59_two_sign_stage_classes: u64,
    pub p59_self_inverse_stage_classes: u64,
    pub challenge59_regression: bool,
    pub safe_ecdh_rejected_every_probe: bool,
    pub wrong_a_rejected: bool,
    pub mutated_singular_order_rejected: bool,
    pub mutated_factor_rejected: bool,
    pub mutated_target_rejected: bool,
    pub one_bit_tag_mutation_rejected: bool,
    pub raw_full_point_model_executed: bool,
    pub raw_full_point_model_recovered_mod_m: bool,
    pub raw_full_point_model_primary: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct RecoveryReport {
    pub schema: String,
    pub experiment_id: String,
    pub seed_hex: String,
    pub seed_digest_hex: String,
    pub planted_d: String,
    pub planted_d_hex: String,
    pub public_key: PointHex,
    pub certification: CertificationReport,
    pub residue: ResidueReport,
    pub interval: IntervalCertificate,
    pub controls: ControlReport,
    pub success: bool,
    pub certificate_sha256: String,
}

#[derive(Clone, Debug)]
pub struct RecoveryConfig {
    pub seed: Vec<u8>,
    pub threads: usize,
    pub batch_lanes: usize,
    pub table_bits: u32,
}

impl RecoveryConfig {
    pub fn frozen(seed: Vec<u8>, threads: usize, batch_lanes: usize, table_bits: u32) -> Self {
        Self {
            seed,
            threads,
            batch_lanes,
            table_bits,
        }
    }
}

pub fn derive_planted_key(seed: &[u8]) -> ([u8; 32], BigUint) {
    let mut input = Vec::with_capacity(SEED_DOMAIN.len() + seed.len());
    input.extend_from_slice(SEED_DOMAIN);
    input.extend_from_slice(seed);
    let digest = sha256(&input);
    let order_minus_one = p192_order() - BigUint::one();
    let scalar = BigUint::one() + BigUint::from_bytes_be(&digest) % order_minus_one;
    (digest, scalar)
}

fn zx_from_native(point: Option<AffinePoint>, arithmetic: &P192Arithmetic) -> [u8; 25] {
    let mut zx = [0u8; 25];
    if let Some(point) = point {
        zx[0] = 1;
        zx[1..].copy_from_slice(&arithmetic.to_bytes(point.x));
    }
    zx
}

fn zx_from_repo(point: &Point) -> Result<[u8; 25], String> {
    let mut zx = [0u8; 25];
    match point {
        Point::Infinity => Ok(zx),
        Point::Affine { x, .. } => {
            let bytes = x.value.to_bytes_be();
            if bytes.len() > 24 {
                return Err("victim returned a non-P-192 x-coordinate".into());
            }
            zx[0] = 1;
            zx[25 - bytes.len()..].copy_from_slice(&bytes);
            Ok(zx)
        }
    }
}

fn confirmation_message(stage: u32, probe: AffinePoint, arithmetic: &P192Arithmetic) -> Vec<u8> {
    let mut message = Vec::with_capacity(CONFIRM_DOMAIN.len() + 4 + 49);
    message.extend_from_slice(CONFIRM_DOMAIN);
    message.extend_from_slice(&stage.to_be_bytes());
    message.extend_from_slice(&probe.sec1_24(arithmetic));
    message
}

fn tag_from_zx(zx: &[u8; 25], message: &[u8]) -> [u8; 32] {
    let mut key_input = Vec::with_capacity(KDF_DOMAIN.len() + zx.len());
    key_input.extend_from_slice(KDF_DOMAIN);
    key_input.extend_from_slice(zx);
    let key = sha256(&key_input);
    hmac_sha256(&key, message)
}

fn offline_tag(
    shared: Option<AffinePoint>,
    message: &[u8],
    arithmetic: &P192Arithmetic,
) -> [u8; 32] {
    tag_from_zx(&zx_from_native(shared, arithmetic), message)
}

fn victim_tag(
    curve: &CurveParams,
    probe: AffinePoint,
    secret: &BigUint,
    stage: u32,
    arithmetic: &P192Arithmetic,
) -> Result<([u8; 32], [u8; 25]), String> {
    let repo_probe = probe.to_repo_point(arithmetic, curve);
    let shared = scalar_mul_secret(&repo_probe, secret, curve);
    let zx = zx_from_repo(&shared)?;
    let message = confirmation_message(stage, probe, arithmetic);
    Ok((tag_from_zx(&zx, &message), zx))
}

fn safe_ecdh_rejects(
    curve: &CurveParams,
    probe: AffinePoint,
    secret: &BigUint,
    arithmetic: &P192Arithmetic,
) -> bool {
    let private = EccPrivateKey {
        scalar: secret.clone(),
        curve_name: curve.name.to_string(),
    };
    let public = EccPublicKey {
        point: probe.to_repo_point(arithmetic, curve),
        curve_name: curve.name.to_string(),
    };
    ecdh_raw(&private, &public, curve).is_none()
}

fn certify() -> Result<CertificationReport, String> {
    let curve = CurveParams::p192();
    let tuple_exact = curve.p == hexadecimal(P_HEX)
        && curve.a == &curve.p - BigUint::from(3u8)
        && curve.b == hexadecimal(B_HEX)
        && curve.gx == hexadecimal(GX_HEX)
        && curve.gy == hexadecimal(GY_HEX)
        && curve.n == hexadecimal(N_HEX)
        && curve.h == 1;
    let certificates = pocklington_certificates()?;
    let p_and_n_proven_prime = certificates.len() == 5
        && certificates
            .iter()
            .all(|certificate| certificate.verified_prime);

    let mut arithmetic = P192Arithmetic::new();
    let g = AffinePoint::generator(&arithmetic);
    let generator_on_p192 = is_on_p192(g, &mut arithmetic);
    let generator_order_n = scalar_mul(g, &p192_order(), &mut arithmetic).is_infinity(&arithmetic);
    let r = AffinePoint::singular_generator(&arithmetic);
    let singular_ok = is_on_singular_companion(r, &mut arithmetic);
    let r_off_p192 = !is_on_p192(r, &mut arithmetic);
    let torus_r = TorusPoint::from_t(8, &arithmetic);
    let mapped = torus_r.to_curve(&mut arithmetic);
    let torus_parameterization_verified = mapped == Some(r);
    let r_order = torus_r
        .scalar_mul(&singular_order(), &mut arithmetic)
        .is_identity(&arithmetic)
        && SINGULAR_DISTINCT_FACTORS.iter().all(|factor| {
            !torus_r
                .scalar_mul(
                    &(singular_order() / BigUint::from(*factor)),
                    &mut arithmetic,
                )
                .is_identity(&arithmetic)
        });
    let factor_primes = SINGULAR_DISTINCT_FACTORS.iter().copied().all(is_prime_u64);
    let factorization_exact =
        accepts_singular_factorization(&singular_order(), &singular_factorization());
    let node = AffinePoint {
        x: arithmetic
            .from_biguint(&(p192_prime() - BigUint::one()))
            .expect("-1 is a canonical P-192 field element"),
        y: arithmetic.zero(),
    };
    Ok(CertificationReport {
        standardized_tuple_exact: tuple_exact,
        pocklington_certificates: certificates,
        p_and_n_proven_prime,
        generator_on_p192,
        generator_order_n,
        singular_discriminant_zero: 4 * (-3i64).pow(3) + 27 * (-2i64).pow(2) == 0,
        singular_factorization_exact: factorization_exact
            && singular_factor_product() == singular_order(),
        singular_factor_primes_exact: factor_primes,
        r_finite_smooth_and_off_p192: singular_ok && r_off_p192 && r != node,
        r_order_n_singular: r_order,
        torus_parameterization_verified,
        operations: arithmetic.counts,
    })
}

fn exact_probe_order(
    probe: TorusPoint,
    order: &BigUint,
    factors: &[u64],
    arithmetic: &mut P192Arithmetic,
) -> bool {
    probe.scalar_mul(order, arithmetic).is_identity(arithmetic)
        && factors.iter().all(|factor| {
            !probe
                .scalar_mul(&(order / BigUint::from(*factor)), arithmetic)
                .is_identity(arithmetic)
        })
}

fn factor_schedule() -> Vec<u64> {
    let mut schedule = vec![2u64; 64];
    schedule.extend([3, 5, 17, 257, 641, 65_537, 274_177, 6_700_417]);
    schedule
}

fn recover_residue_schedule(
    secret: &BigUint,
    schedule: &[u64],
    batch_lanes: usize,
) -> Result<ResidueReport, String> {
    let curve = CurveParams::p192();
    let mut arithmetic = P192Arithmetic::new();
    let r_generator = TorusPoint::from_t(8, &arithmetic);
    let n_singular = singular_order();
    let mut modulus = BigUint::one();
    let mut residue = BigUint::zero();
    let mut distinct_factors = Vec::<u64>::new();
    let mut stages = Vec::new();
    let mut total_comparisons = 0u64;
    let mut transcript = Vec::new();

    for (index, factor) in schedule.iter().copied().enumerate() {
        if n_singular.clone() % (&modulus * BigUint::from(factor)) != BigUint::zero() {
            return Err("residue schedule does not divide the singular order".into());
        }
        let before = modulus.clone();
        let after = &before * BigUint::from(factor);
        if !distinct_factors.contains(&factor) {
            distinct_factors.push(factor);
        }
        let probe_torus = r_generator.scalar_mul(&(&n_singular / &after), &mut arithmetic);
        let probe_order =
            exact_probe_order(probe_torus, &after, &distinct_factors, &mut arithmetic);
        if !probe_order {
            return Err(format!(
                "stage {} probe lacks exact order {after}",
                index + 1
            ));
        }
        let probe = probe_torus
            .to_curve(&mut arithmetic)
            .ok_or_else(|| "nontrivial residue probe mapped to Infinity".to_string())?;
        if !is_on_singular_companion(probe, &mut arithmetic) || is_on_p192(probe, &mut arithmetic) {
            return Err("residue probe failed curve-membership separation".into());
        }
        let stage_number = (index + 1) as u32;
        let (tag, victim_zx) = victim_tag(&curve, probe, secret, stage_number, &arithmetic)?;
        let independent_shared = probe_torus.scalar_mul(secret, &mut arithmetic);
        let independent_curve = independent_shared.to_curve(&mut arithmetic);
        let independent_zx = zx_from_native(independent_curve, &arithmetic);
        let low_level_match = victim_zx == independent_zx;
        if !low_level_match {
            return Err(format!("stage {stage_number} low-level/torus mismatch"));
        }
        let safe_rejected = safe_ecdh_rejects(&curve, probe, secret, &arithmetic);
        if !safe_rejected {
            return Err(format!(
                "stage {stage_number} safe ECDH accepted singular input"
            ));
        }

        let message = confirmation_message(stage_number, probe, &arithmetic);
        let step = probe_torus.scalar_mul(&before, &mut arithmetic);
        let mut current = probe_torus.scalar_mul(&residue, &mut arithmetic);
        let mut matched = None;
        let mut j = 0u64;
        while j < factor && matched.is_none() {
            let take = (factor - j).min(batch_lanes.max(1) as u64) as usize;
            let mut candidates = Vec::with_capacity(take);
            for _ in 0..take {
                candidates.push(current);
                current = current.compose(step, &mut arithmetic);
            }
            let curve_candidates = batch_torus_to_curve(&candidates, &mut arithmetic);
            for (offset, candidate) in curve_candidates.into_iter().enumerate() {
                total_comparisons += 1;
                let candidate_tag = offline_tag(candidate, &message, &arithmetic);
                if candidate_tag == tag {
                    let lift_index = j + offset as u64;
                    matched = Some(&residue + &before * BigUint::from(lift_index));
                    break;
                }
            }
            j += take as u64;
        }
        let matched = matched.ok_or_else(|| format!("stage {stage_number} found no lift"))?;
        let comparisons_this_stage = total_comparisons
            - stages
                .iter()
                .map(|stage: &ResidueStage| stage.candidate_comparisons)
                .sum::<u64>();
        residue = matched;
        modulus = after;
        let negative = if residue.is_zero() {
            BigUint::zero()
        } else {
            &modulus - &residue
        };
        transcript.extend_from_slice(&stage_number.to_be_bytes());
        transcript.extend_from_slice(&tag);
        transcript.extend_from_slice(&modulus.to_bytes_be());
        transcript.extend_from_slice(&residue.to_bytes_be());
        stages.push(ResidueStage {
            stage: stage_number,
            factor,
            modulus_before: before.to_string(),
            modulus_after: modulus.to_string(),
            probe: point_hex(probe, &arithmetic),
            probe_exact_order: probe_order,
            primary_tag_hex: hex::encode(tag),
            candidate_comparisons: comparisons_this_stage,
            matched_lift: residue.to_string(),
            retained_negative: negative.to_string(),
            low_level_matches_independent_torus: low_level_match,
            safe_ecdh_rejected: safe_rejected,
        });
    }

    let negative = if residue.is_zero() {
        BigUint::zero()
    } else {
        &modulus - &residue
    };
    let expected = [
        decimal(EXPECTED_RESIDUE_SMALL),
        decimal(EXPECTED_RESIDUE_LARGE),
    ];
    let expected_match = modulus == residue_modulus()
        && ((residue == expected[0] && negative == expected[1])
            || (residue == expected[1] && negative == expected[0]));
    Ok(ResidueReport {
        mode: "nested-composite/x-only".into(),
        primary_oracle_queries: stages.len() as u64,
        candidate_comparisons: total_comparisons,
        residue: residue.to_string(),
        negative_residue: negative.to_string(),
        expected_unordered_residues_match: expected_match,
        transcript_sha256: hex::encode(sha256(&transcript)),
        stages,
        operations: arithmetic.counts,
    })
}

pub fn recover_frozen_residue(
    secret: &BigUint,
    batch_lanes: usize,
) -> Result<ResidueReport, String> {
    recover_residue_schedule(secret, &factor_schedule(), batch_lanes)
}

/// Stronger positive-model lane: expose the complete singular-companion
/// result `(x,y)` and run the same nested smooth-factor recovery.  This is
/// deliberately separate from the primary x-only transcript and is never
/// used to seed or prune the residual BSGS.
fn recover_full_point_model_schedule(
    secret: &BigUint,
    schedule: &[u64],
    batch_lanes: usize,
) -> Result<(BigUint, BigUint, u64), String> {
    let curve = CurveParams::p192();
    let mut arithmetic = P192Arithmetic::new();
    let generator = TorusPoint::from_t(8, &arithmetic);
    let mut modulus = BigUint::one();
    let mut residue = BigUint::zero();
    let mut queries = 0u64;
    for factor in schedule {
        let before = modulus.clone();
        let after = &before * BigUint::from(*factor);
        let probe_torus = generator.scalar_mul(&(singular_order() / &after), &mut arithmetic);
        let probe = probe_torus
            .to_curve(&mut arithmetic)
            .ok_or_else(|| "full-point model probe is Infinity".to_string())?;
        let victim = scalar_mul_secret(&probe.to_repo_point(&arithmetic, &curve), secret, &curve);
        let observed = AffinePoint::from_repo_point(&victim, &arithmetic)?;
        queries += 1;

        let step = probe_torus.scalar_mul(&before, &mut arithmetic);
        let mut current = probe_torus.scalar_mul(&residue, &mut arithmetic);
        let mut matched = None;
        let mut lift_index = 0u64;
        while lift_index < *factor && matched.is_none() {
            let take = (*factor - lift_index).min(batch_lanes.max(1) as u64) as usize;
            let mut candidates = Vec::with_capacity(take);
            for _ in 0..take {
                candidates.push(current);
                current = current.compose(step, &mut arithmetic);
            }
            for (offset, candidate) in batch_torus_to_curve(&candidates, &mut arithmetic)
                .into_iter()
                .enumerate()
            {
                if candidate == observed {
                    matched = Some(&residue + &before * BigUint::from(lift_index + offset as u64));
                    break;
                }
            }
            lift_index += take as u64;
        }
        residue = matched.ok_or_else(|| "full-point model found no exact lift".to_string())?;
        modulus = after;
    }
    Ok((residue, modulus, queries))
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct P59ControlOutcome {
    passed: bool,
    scalars_tested: u64,
    nested_queries: u64,
    candidate_comparisons: u64,
    two_sign_stage_classes: u64,
    self_inverse_stage_classes: u64,
}

fn p59_control() -> P59ControlOutcome {
    const P: u64 = 59;
    const ORDER: u64 = 60;
    fn inverse(value: u64) -> u64 {
        pow_mod_u64(value, P - 2, P)
    }
    fn compose(left: (u64, u64), right: (u64, u64)) -> (u64, u64) {
        (
            (mul_mod_u64(left.0, right.0, P) + P
                - mul_mod_u64(3, mul_mod_u64(left.1, right.1, P), P))
                % P,
            (mul_mod_u64(left.0, right.1, P) + mul_mod_u64(right.0, left.1, P)) % P,
        )
    }
    fn scalar(base: (u64, u64), mut k: u64) -> (u64, u64) {
        let mut acc = (1, 0);
        let mut power = base;
        while k != 0 {
            if k & 1 == 1 {
                acc = compose(acc, power);
            }
            power = compose(power, power);
            k >>= 1;
        }
        acc
    }
    fn x_of(point: (u64, u64)) -> Option<u64> {
        if point.1 == 0 {
            return None;
        }
        let t = mul_mod_u64(point.0, inverse(point.1), P);
        Some((mul_mod_u64(t, t, P) + 2) % P)
    }
    let generator = (2u64..P).map(|t| (t, 1)).find(|point| {
        scalar(*point, ORDER).1 == 0
            && [2u64, 3, 5]
                .iter()
                .all(|q| scalar(*point, ORDER / q).1 != 0)
    });
    let Some(generator) = generator else {
        return P59ControlOutcome::default();
    };
    let mut outcome = P59ControlOutcome {
        passed: true,
        ..P59ControlOutcome::default()
    };
    for secret in 0..ORDER {
        outcome.scalars_tested += 1;
        // Independent full-order preimage census: the x-only observation has
        // exactly the expected {d,-d} class, including the self-inverse cases.
        let observed = x_of(scalar(generator, secret));
        let matches: Vec<_> = (0..ORDER)
            .filter(|candidate| x_of(scalar(generator, *candidate)) == observed)
            .collect();
        let mut expected = vec![secret, (ORDER - secret) % ORDER];
        expected.sort_unstable();
        expected.dedup();
        if matches != expected {
            outcome.passed = false;
            return outcome;
        }

        // Execute the same nested-composite protocol as the primary lane at
        // reduced width.  Every one of the 60 scalars gets four actual staged
        // x-queries, least-lift scanning, and a coherent retained sign class.
        let mut modulus = 1u64;
        let mut residue = 0u64;
        let mut queries_for_secret = 0u64;
        for factor in [2u64, 2, 3, 5] {
            let after = modulus * factor;
            let probe = scalar(generator, ORDER / after);
            if scalar(probe, after).1 != 0
                || [2u64, 3, 5]
                    .iter()
                    .filter(|prime| after % **prime == 0)
                    .any(|prime| scalar(probe, after / prime).1 == 0)
            {
                outcome.passed = false;
                return outcome;
            }
            let observed = x_of(scalar(probe, secret));
            let mut matched = None;
            for lift_index in 0..factor {
                outcome.candidate_comparisons += 1;
                let candidate = residue + lift_index * modulus;
                if x_of(scalar(probe, candidate)) == observed {
                    matched = Some(candidate);
                    break;
                }
            }
            let Some(next_residue) = matched else {
                outcome.passed = false;
                return outcome;
            };
            queries_for_secret += 1;
            outcome.nested_queries += 1;
            let d_mod = secret % after;
            let negative = (after - next_residue) % after;
            if next_residue != d_mod && negative != d_mod {
                outcome.passed = false;
                return outcome;
            }
            if next_residue == negative {
                outcome.self_inverse_stage_classes += 1;
            } else {
                outcome.two_sign_stage_classes += 1;
            }
            residue = next_residue;
            modulus = after;
        }
        if queries_for_secret != 4
            || modulus != ORDER
            || (residue != secret && (ORDER - residue) % ORDER != secret)
        {
            outcome.passed = false;
            return outcome;
        }
    }
    if outcome.nested_queries != ORDER * 4
        || outcome.two_sign_stage_classes + outcome.self_inverse_stage_classes
            != outcome.nested_queries
    {
        outcome.passed = false;
    }
    outcome
}

fn challenge59_control() -> bool {
    use crate::cryptopals::challenge59::{attack, cryptopals_curve, Alice};
    let (curve, order, _) = cryptopals_curve();
    let secret = decimal("987654321987654321987654321");
    let alice = Alice {
        curve: curve.clone(),
        d: secret.clone(),
    };
    attack(&alice, &curve, &order).is_some_and(|recovered| recovered == secret)
}

fn frozen_tail_covered() -> bool {
    let positions = QMAX.div_ceil(FROZEN_STRIDE);
    if positions != 5_800_018 {
        return false;
    }
    let k = QMAX - 1;
    let index = k / FROZEN_STRIDE;
    let center = (FROZEN_M - 1) + index * FROZEN_STRIDE;
    let baby = center.abs_diff(k);
    index < positions && baby < FROZEN_M
}

fn accepts_confirmation_tag(expected: &[u8; 32], presented: &[u8; 32]) -> bool {
    // Tags are public experiment outputs here; the production HMAC verifier
    // remains constant-time.  This predicate is the transcript acceptance
    // gate exercised by both wrong-arithmetic and bit-flip controls.
    expected == presented
}

fn accepts_recovered_target(
    candidate: &BigUint,
    target: super::p192_native::JacobianPoint,
    arithmetic: &mut P192Arithmetic,
) -> bool {
    if candidate.is_zero() || candidate >= &p192_order() {
        return false;
    }
    let replay = scalar_mul(AffinePoint::generator(arithmetic), candidate, arithmetic);
    super::p192_native::affine_equal(replay, target, arithmetic)
}

fn controls(
    secret: &BigUint,
    q: AffinePoint,
    residue: &ResidueReport,
    arithmetic: &mut P192Arithmetic,
) -> Result<ControlReport, String> {
    let curve = CurveParams::p192();
    let last = residue
        .stages
        .last()
        .ok_or_else(|| "residue report contains no stages".to_string())?;
    let probe = point_from_hex(&last.probe, arithmetic)?;
    let stage = last.stage;
    let (normal_tag, _) = victim_tag(&curve, probe, secret, stage, arithmetic)?;
    let mut wrong_a = curve.clone();
    wrong_a.a = &wrong_a.p - BigUint::from(4u8);
    let (wrong_a_tag, _) = victim_tag(&wrong_a, probe, secret, stage, arithmetic)?;

    let factors = singular_factorization();
    let mutated_order_rejected =
        !accepts_singular_factorization(&(singular_order() + BigUint::one()), &factors);
    let mut mutated_factors = factors.clone();
    mutated_factors[0].1 -= 1;
    let mutated_factor_rejected =
        !accepts_singular_factorization(&singular_order(), &mutated_factors);

    let g = AffinePoint::generator(arithmetic);
    let mutated_q =
        super::p192_native::JacobianPoint::from_affine(q, arithmetic).add_affine(g, arithmetic);
    let mutated_target_rejected = !accepts_recovered_target(secret, mutated_q, arithmetic);

    let mut flipped = normal_tag;
    flipped[0] ^= 1;
    let message = confirmation_message(stage, probe, arithmetic);
    let independent = TorusPoint::from_t(8, arithmetic)
        .scalar_mul(&(singular_order() / residue_modulus()), arithmetic)
        .scalar_mul(secret, arithmetic)
        .to_curve(arithmetic);
    let expected = offline_tag(independent, &message, arithmetic);
    let p59 = p59_control();
    let (raw_residue, raw_modulus, raw_queries) =
        recover_full_point_model_schedule(secret, &factor_schedule(), 256)?;
    let raw_model_ok = raw_modulus == residue_modulus()
        && raw_queries == 72
        && raw_residue == secret % &raw_modulus;
    Ok(ControlReport {
        exhaustive_p59_sign_classes: p59.passed,
        p59_scalars_tested: p59.scalars_tested,
        p59_nested_queries: p59.nested_queries,
        p59_candidate_comparisons: p59.candidate_comparisons,
        p59_two_sign_stage_classes: p59.two_sign_stage_classes,
        p59_self_inverse_stage_classes: p59.self_inverse_stage_classes,
        challenge59_regression: challenge59_control(),
        safe_ecdh_rejected_every_probe: residue.stages.iter().all(|row| row.safe_ecdh_rejected),
        wrong_a_rejected: accepts_confirmation_tag(&expected, &normal_tag)
            && !accepts_confirmation_tag(&expected, &wrong_a_tag),
        mutated_singular_order_rejected: mutated_order_rejected,
        mutated_factor_rejected,
        mutated_target_rejected,
        one_bit_tag_mutation_rejected: accepts_confirmation_tag(&expected, &normal_tag)
            && !accepts_confirmation_tag(&expected, &flipped),
        raw_full_point_model_executed: raw_queries == 72,
        raw_full_point_model_recovered_mod_m: raw_model_ok,
        raw_full_point_model_primary: false,
    })
}

fn report_digest(report: &RecoveryReport) -> Result<String, String> {
    let mut unsigned = report.clone();
    unsigned.certificate_sha256.clear();
    let bytes = serde_json::to_vec(&unsigned)
        .map_err(|error| format!("failed to encode recovery certificate: {error}"))?;
    Ok(hex::encode(sha256(&bytes)))
}

pub fn run_recovery(config: &RecoveryConfig) -> Result<RecoveryReport, String> {
    if hex::encode(&config.seed) != FROZEN_SEED_HEX {
        return Err("the approved protocol permits only the frozen planted-key seed".into());
    }
    if config.table_bits != FROZEN_TABLE_BITS {
        return Err("crypto-scale recovery requires --table-bits 24".into());
    }
    let plan = IntervalPlan::frozen(config.threads, config.batch_lanes)?;
    let certification = certify()?;
    let cert_ok = certification_succeeds(&certification);
    if !cert_ok {
        return Err("mandatory P-192/singular certification failed".into());
    }

    let (digest, secret) = derive_planted_key(&config.seed);
    let curve = CurveParams::p192();
    let public_repo = scalar_mul_secret(&curve.generator(), &secret, &curve);
    let mut arithmetic = P192Arithmetic::new();
    let public = AffinePoint::from_repo_point(&public_repo, &arithmetic)?
        .ok_or_else(|| "planted P-192 public key is Infinity".to_string())?;
    let independent_public = scalar_mul(
        AffinePoint::generator(&arithmetic),
        &secret,
        &mut arithmetic,
    );
    if !super::p192_native::affine_equal(
        independent_public,
        super::p192_native::JacobianPoint::from_affine(public, &arithmetic),
        &mut arithmetic,
    ) {
        return Err("actual/native planted public keys disagree".into());
    }

    let residue = recover_frozen_residue(&secret, config.batch_lanes)?;
    if residue.primary_oracle_queries != 72
        || !residue.expected_unordered_residues_match
        || residue.candidate_comparisons != EXPECTED_CANDIDATE_COMPARISONS
    {
        return Err("frozen nested-composite tail checks failed".into());
    }
    let plus_residue = decimal(&residue.residue);
    let minus_residue = decimal(&residue.negative_residue);
    let order = p192_order();
    let modulus = residue_modulus();
    let g = AffinePoint::generator(&arithmetic);
    let h = scalar_mul(g, &modulus, &mut arithmetic)
        .to_affine(&mut arithmetic)
        .ok_or_else(|| "[M]G is Infinity".to_string())?;
    let plus = make_target(
        1,
        public,
        g,
        &plus_residue,
        &modulus,
        &order,
        &mut arithmetic,
    )?;
    let minus = make_target(
        -1,
        public,
        g,
        &minus_residue,
        &modulus,
        &order,
        &mut arithmetic,
    )?;
    let plus_width = plus.max_k + 1;
    let minus_width = minus.max_k + 1;
    let replay_context = super::p192_interval_bsgs::FinalReplayContext {
        generator: g,
        public_key: public,
        plus_residue: plus_residue.clone(),
        minus_residue: minus_residue.clone(),
        modulus: modulus.clone(),
        order: order.clone(),
    };
    let (solution, interval_stats) = solve_two_orientations(h, plus, minus, plan, &replay_context)?;
    let selected_residue = if solution.orientation == 1 {
        &plus_residue
    } else {
        &minus_residue
    };
    let recovered = selected_residue + &modulus * BigUint::from(solution.k);
    arithmetic.reset_counts();
    let replay = scalar_mul(g, &recovered, &mut arithmetic);
    let replay_ok = super::p192_native::affine_equal(
        replay,
        super::p192_native::JacobianPoint::from_affine(public, &arithmetic),
        &mut arithmetic,
    );
    let final_replay_operations = arithmetic.counts;
    let interval = IntervalCertificate {
        geometry: SHIFTED_GEOMETRY.into(),
        shifted_centers: true,
        plan,
        plus_residue: plus_residue.to_string(),
        minus_residue: minus_residue.to_string(),
        plus_width,
        minus_width,
        solution,
        stats: interval_stats,
        recovered_d: recovered.to_string(),
        recovered_d_hex: format!("{recovered:048x}"),
        final_replay_matches_public_key: replay_ok,
        final_replay_operations,
        frozen_tail_covered: frozen_tail_covered(),
    };
    let controls = controls(&secret, public, &residue, &mut arithmetic)?;
    let controls_ok = controls_succeed(&controls);
    let success = cert_ok
        && controls_ok
        && interval.final_replay_matches_public_key
        && interval.frozen_tail_covered
        && recovered == secret
        && solution_k_and_orientation_ok(&interval);
    let mut report = RecoveryReport {
        schema: REPORT_SCHEMA.into(),
        experiment_id: EXPERIMENT_ID.into(),
        seed_hex: hex::encode(&config.seed),
        seed_digest_hex: hex::encode(digest),
        planted_d: secret.to_string(),
        planted_d_hex: format!("{secret:048x}"),
        public_key: point_hex(public, &arithmetic),
        certification,
        residue,
        interval,
        controls,
        success,
        certificate_sha256: String::new(),
    };
    report.certificate_sha256 = report_digest(&report)?;
    if !report.success {
        return Err("recovery completed but a mandatory success gate failed".into());
    }
    Ok(report)
}

fn solution_k_and_orientation_ok(interval: &IntervalCertificate) -> bool {
    interval.solution.orientation == -1
        && interval.solution.k == EXPECTED_RESIDUAL_K
        && interval.solution.giant_index == 3_159_226
        && interval.solution.baby_index == 989_698
        && interval.solution.baby_sign == -1
}

const SHIFTED_GEOMETRY: &str = "shifted-negation-v2";

fn certification_succeeds(certification: &CertificationReport) -> bool {
    certification.standardized_tuple_exact
        && certification.p_and_n_proven_prime
        && certification.generator_on_p192
        && certification.generator_order_n
        && certification.singular_discriminant_zero
        && certification.singular_factorization_exact
        && certification.singular_factor_primes_exact
        && certification.r_finite_smooth_and_off_p192
        && certification.r_order_n_singular
        && certification.torus_parameterization_verified
}

fn controls_succeed(controls: &ControlReport) -> bool {
    controls.exhaustive_p59_sign_classes
        && controls.p59_scalars_tested == 60
        && controls.p59_nested_queries == 240
        && controls.p59_candidate_comparisons >= controls.p59_nested_queries
        && controls.p59_two_sign_stage_classes + controls.p59_self_inverse_stage_classes == 240
        && controls.challenge59_regression
        && controls.safe_ecdh_rejected_every_probe
        && controls.wrong_a_rejected
        && controls.mutated_singular_order_rejected
        && controls.mutated_factor_rejected
        && controls.mutated_target_rejected
        && controls.one_bit_tag_mutation_rejected
        && controls.raw_full_point_model_executed
        && controls.raw_full_point_model_recovered_mod_m
        && !controls.raw_full_point_model_primary
}

fn to_u64(value: &BigUint, label: &str) -> Result<u64, String> {
    let digits: Vec<_> = value.iter_u64_digits().collect();
    match digits.as_slice() {
        [] => Ok(0),
        [value] => Ok(*value),
        _ => Err(format!("{label} does not fit u64")),
    }
}

fn interval_width(residue: &BigUint) -> Result<u64, String> {
    let max_k = (p192_order() - BigUint::one() - residue) / residue_modulus();
    to_u64(&(max_k + BigUint::one()), "residual interval width")
}

fn verify_residue_shape(residue: &ResidueReport) -> Result<(), String> {
    let schedule = factor_schedule();
    if residue.mode != "nested-composite/x-only"
        || residue.primary_oracle_queries != 72
        || residue.stages.len() != schedule.len()
        || residue.candidate_comparisons != EXPECTED_CANDIDATE_COMPARISONS
        || !residue.expected_unordered_residues_match
    {
        return Err("residue summary differs from the frozen protocol".into());
    }
    if residue.residue != EXPECTED_RESIDUE_SMALL
        || residue.negative_residue != EXPECTED_RESIDUE_LARGE
    {
        return Err("ordered least-lift residue pair differs from the frozen pair".into());
    }
    let mut before = BigUint::one();
    let mut previous_lift = BigUint::zero();
    let mut comparisons = 0u64;
    let mut transcript = Vec::new();
    let mut arithmetic = P192Arithmetic::new();
    for (index, (stage, factor)) in residue.stages.iter().zip(schedule).enumerate() {
        let after = &before * BigUint::from(factor);
        let lift = decimal(&stage.matched_lift);
        let negative = if lift.is_zero() {
            BigUint::zero()
        } else {
            &after - &lift
        };
        if stage.stage != (index + 1) as u32
            || stage.factor != factor
            || stage.modulus_before != before.to_string()
            || stage.modulus_after != after.to_string()
            || lift >= after
            || (&lift % &before) != previous_lift
            || stage.retained_negative != negative.to_string()
            || !stage.probe_exact_order
            || !stage.low_level_matches_independent_torus
            || !stage.safe_ecdh_rejected
        {
            return Err(format!(
                "residue stage {} is structurally invalid",
                index + 1
            ));
        }
        let tag = hex::decode(&stage.primary_tag_hex)
            .map_err(|error| format!("invalid stage tag hex: {error}"))?;
        if tag.len() != 32 {
            return Err("primary HMAC tag is not 32 bytes".into());
        }
        let probe = point_from_hex(&stage.probe, &arithmetic)?;
        if !is_on_singular_companion(probe, &mut arithmetic) || is_on_p192(probe, &mut arithmetic) {
            return Err("reported probe is not separated onto the singular companion".into());
        }
        comparisons = comparisons
            .checked_add(stage.candidate_comparisons)
            .ok_or_else(|| "candidate comparison count overflow".to_string())?;
        transcript.extend_from_slice(&stage.stage.to_be_bytes());
        transcript.extend_from_slice(&tag);
        transcript.extend_from_slice(&after.to_bytes_be());
        transcript.extend_from_slice(&lift.to_bytes_be());
        previous_lift = lift;
        before = after;
    }
    if comparisons != residue.candidate_comparisons
        || before != residue_modulus()
        || previous_lift.to_string() != residue.residue
        || hex::encode(sha256(&transcript)) != residue.transcript_sha256
    {
        return Err("residue transcript totals or digest are inconsistent".into());
    }
    Ok(())
}

fn verify_interval_claim(
    report: &RecoveryReport,
    secret: &BigUint,
    public: AffinePoint,
) -> Result<
    (
        AffinePoint,
        super::p192_interval_bsgs::IntervalTarget,
        super::p192_interval_bsgs::IntervalTarget,
        super::p192_interval_bsgs::FinalReplayContext,
    ),
    String,
> {
    let interval = &report.interval;
    interval.plan.validate_frozen()?;
    if interval.geometry != SHIFTED_GEOMETRY
        || !interval.shifted_centers
        || !interval.frozen_tail_covered
        || !frozen_tail_covered()
    {
        return Err("shifted-center geometry or tail certificate is invalid".into());
    }
    if interval.plus_residue != report.residue.residue
        || interval.minus_residue != report.residue.negative_residue
    {
        return Err("interval residues are not bound to the primary transcript".into());
    }
    let plus_residue = decimal(&interval.plus_residue);
    let minus_residue = decimal(&interval.minus_residue);
    if &plus_residue + &minus_residue != residue_modulus() {
        return Err("interval residues are not negatives modulo M".into());
    }
    if interval.plus_width != interval_width(&plus_residue)?
        || interval.minus_width != interval_width(&minus_residue)?
    {
        return Err("reported interval widths do not replay".into());
    }
    if !solution_k_and_orientation_ok(interval) {
        return Err("BSGS hit differs from the frozen hit certificate".into());
    }
    let solution = &interval.solution;
    if solution.baby_index >= interval.plan.m {
        return Err("BSGS baby index is outside the frozen table".into());
    }
    let center = (interval.plan.m - 1)
        .checked_add(
            solution
                .giant_index
                .checked_mul(interval.plan.stride)
                .ok_or_else(|| "BSGS center overflow".to_string())?,
        )
        .ok_or_else(|| "BSGS center overflow".to_string())?;
    let reconstructed_k = match solution.baby_sign {
        -1 => center.checked_sub(solution.baby_index),
        0 if solution.baby_index == 0 => Some(center),
        1 => center.checked_add(solution.baby_index),
        _ => None,
    }
    .ok_or_else(|| "BSGS center/index/sign certificate is invalid".to_string())?;
    if reconstructed_k != solution.k {
        return Err("BSGS center/index/sign does not reconstruct k".into());
    }
    let selected_width = if solution.orientation == 1 {
        interval.plus_width
    } else {
        interval.minus_width
    };
    if solution.k >= selected_width {
        return Err("BSGS k is outside the selected inclusive interval".into());
    }

    let stats = &interval.stats;
    if stats.table_bits != interval.plan.table_bits
        || stats.table_bytes != interval.plan.table_bytes()
        || stats.baby_entries != interval.plan.m - 1
        || stats.baby_chain_steps != interval.plan.m - 1
        || stats.giant_positions_plus != solution.giant_index + 1
        || stats.giant_positions_minus != solution.giant_index + 1
        || stats.giant_positions_plus > super::p192_interval_bsgs::FROZEN_MAX_GIANT_STEPS
        || stats.giant_positions_minus > super::p192_interval_bsgs::FROZEN_MAX_GIANT_STEPS
        || stats.fingerprint_hits == 0
        || stats.exact_replays == 0
        || stats.false_fingerprint_hits > stats.fingerprint_hits
    {
        return Err("BSGS operation/table statistics are inconsistent with the hit".into());
    }
    let mut combined_operations = stats.baby_operations;
    combined_operations.merge(stats.giant_and_replay_operations);
    if combined_operations != stats.operations {
        return Err("BSGS phase operation counts do not sum to the total".into());
    }

    let mut arithmetic = P192Arithmetic::new();
    let g = AffinePoint::generator(&arithmetic);
    let modulus = residue_modulus();
    let order = p192_order();
    let h = scalar_mul(g, &modulus, &mut arithmetic)
        .to_affine(&mut arithmetic)
        .ok_or_else(|| "[M]G is Infinity while verifying interval".to_string())?;
    let plus = make_target(
        1,
        public,
        g,
        &plus_residue,
        &modulus,
        &order,
        &mut arithmetic,
    )?;
    let minus = make_target(
        -1,
        public,
        g,
        &minus_residue,
        &modulus,
        &order,
        &mut arithmetic,
    )?;
    if plus.max_k + 1 != interval.plus_width || minus.max_k + 1 != interval.minus_width {
        return Err("make_target interval widths disagree with the report".into());
    }
    let selected_target = if solution.orientation == 1 {
        plus
    } else {
        minus
    };
    let quotient_replay = super::p192_native::scalar_mul_u64(h, solution.k, &mut arithmetic);
    let quotient_matches = match selected_target.point {
        None => quotient_replay.is_infinity(&arithmetic),
        Some(target) => super::p192_native::affine_equal(
            quotient_replay,
            super::p192_native::JacobianPoint::from_affine(target, &arithmetic),
            &mut arithmetic,
        ),
    };
    if !quotient_matches {
        return Err("BSGS quotient equation [k]H=T does not replay".into());
    }
    let selected_residue = if solution.orientation == 1 {
        &plus_residue
    } else {
        &minus_residue
    };
    let recovered = selected_residue + &modulus * BigUint::from(solution.k);
    if recovered.is_zero()
        || recovered >= order
        || recovered != *secret
        || interval.recovered_d != recovered.to_string()
        || interval.recovered_d_hex != format!("{recovered:048x}")
    {
        return Err("d=r+kM reconstruction or range check failed".into());
    }
    let target = super::p192_native::JacobianPoint::from_affine(public, &arithmetic);
    arithmetic.reset_counts();
    if !accepts_recovered_target(&recovered, target, &mut arithmetic)
        || !interval.final_replay_matches_public_key
        || arithmetic.counts != interval.final_replay_operations
    {
        return Err("final legitimate-curve replay or its operation count is invalid".into());
    }
    let replay = super::p192_interval_bsgs::FinalReplayContext {
        generator: g,
        public_key: public,
        plus_residue,
        minus_residue,
        modulus,
        order,
    };
    Ok((h, plus, minus, replay))
}

fn verify_report_preflight(report: &RecoveryReport) -> Result<(BigUint, AffinePoint), String> {
    if report.schema != REPORT_SCHEMA || report.experiment_id != EXPERIMENT_ID {
        return Err("unknown recovery report schema or experiment".into());
    }
    if report.certificate_sha256 != report_digest(report)? {
        return Err("recovery certificate digest mismatch".into());
    }
    if report.seed_hex != FROZEN_SEED_HEX {
        return Err("report does not use the single frozen seed".into());
    }
    let seed =
        hex::decode(&report.seed_hex).map_err(|error| format!("invalid seed hex: {error}"))?;
    let (digest, secret) = derive_planted_key(&seed);
    if report.seed_digest_hex != EXPECTED_DIGEST_HEX
        || report.seed_digest_hex != hex::encode(digest)
        || report.planted_d != EXPECTED_D_DEC
        || report.planted_d != secret.to_string()
        || report.planted_d_hex != EXPECTED_D_HEX
        || report.planted_d_hex != format!("{secret:048x}")
    {
        return Err("frozen planted-key derivation does not replay".into());
    }
    if !certification_succeeds(&report.certification) {
        return Err("one or more certification success gates is false".into());
    }
    verify_residue_shape(&report.residue)?;
    let mut arithmetic = P192Arithmetic::new();
    let public = point_from_hex(&report.public_key, &arithmetic)?;
    if !is_on_p192(public, &mut arithmetic) {
        return Err("reported public key is not on P-192".into());
    }
    let curve = CurveParams::p192();
    let actual_public = scalar_mul_secret(&curve.generator(), &secret, &curve);
    if AffinePoint::from_repo_point(&actual_public, &arithmetic)? != Some(public) {
        return Err("reported public key differs from the actual low-level [d]G".into());
    }
    verify_interval_claim(report, &secret, public)?;
    if !controls_succeed(&report.controls) || !report.success {
        return Err("control or overall success gate is false".into());
    }
    Ok((secret, public))
}

/// Verify a completed report.  This replays the primary 72-query residue
/// transcript through the real low-level multiplier, rebuilds the shared BSGS
/// table/search to bind every statistic, and independently re-executes every
/// negative control and final public-key equation.
pub fn verify_report(report: &RecoveryReport) -> Result<(), String> {
    let (secret, public) = verify_report_preflight(report)?;
    let certification = certify()?;
    if certification != report.certification {
        return Err("order/point certification does not replay".into());
    }
    let residue = recover_frozen_residue(&secret, report.interval.plan.batch_lanes)?;
    if residue != report.residue {
        return Err("primary x-only residue transcript does not replay".into());
    }
    let (h, plus, minus, replay_context) = verify_interval_claim(report, &secret, public)?;
    let (solution, stats) =
        solve_two_orientations(h, plus, minus, report.interval.plan, &replay_context)?;
    if solution != report.interval.solution || stats != report.interval.stats {
        return Err("full BSGS replay differs from the reported hit or statistics".into());
    }
    let mut arithmetic = P192Arithmetic::new();
    let replayed_controls = controls(&secret, public, &residue, &mut arithmetic)?;
    if replayed_controls != report.controls || !controls_succeed(&replayed_controls) {
        return Err("negative/control suite does not replay".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synthetic_preflight_report() -> RecoveryReport {
        let seed = hex::decode(FROZEN_SEED_HEX).unwrap();
        let (digest, secret) = derive_planted_key(&seed);
        let mut arithmetic = P192Arithmetic::new();
        let g = AffinePoint::generator(&arithmetic);
        let public = scalar_mul(g, &secret, &mut arithmetic)
            .to_affine(&mut arithmetic)
            .unwrap();
        let probe = AffinePoint::singular_generator(&arithmetic);

        let mut before = BigUint::one();
        let mut current = BigUint::zero();
        let mut comparisons = 0u64;
        let mut transcript = Vec::new();
        let mut stages = Vec::new();
        for (index, factor) in factor_schedule().into_iter().enumerate() {
            let after = &before * BigUint::from(factor);
            let positive = &secret % &after;
            let negative = if positive.is_zero() {
                BigUint::zero()
            } else {
                &after - &positive
            };
            let mut lifts: Vec<BigUint> = [positive, negative]
                .into_iter()
                .filter(|candidate| candidate % &before == current)
                .collect();
            lifts.sort();
            lifts.dedup();
            let lift = lifts[0].clone();
            let lift_index = (&lift - &current) / &before;
            let stage_comparisons = to_u64(&lift_index, "synthetic lift").unwrap() + 1;
            comparisons += stage_comparisons;
            let retained_negative = if lift.is_zero() {
                BigUint::zero()
            } else {
                &after - &lift
            };
            let stage_number = (index + 1) as u32;
            let tag = [stage_number as u8; 32];
            transcript.extend_from_slice(&stage_number.to_be_bytes());
            transcript.extend_from_slice(&tag);
            transcript.extend_from_slice(&after.to_bytes_be());
            transcript.extend_from_slice(&lift.to_bytes_be());
            stages.push(ResidueStage {
                stage: stage_number,
                factor,
                modulus_before: before.to_string(),
                modulus_after: after.to_string(),
                probe: point_hex(probe, &arithmetic),
                probe_exact_order: true,
                primary_tag_hex: hex::encode(tag),
                candidate_comparisons: stage_comparisons,
                matched_lift: lift.to_string(),
                retained_negative: retained_negative.to_string(),
                low_level_matches_independent_torus: true,
                safe_ecdh_rejected: true,
            });
            current = lift;
            before = after;
        }
        assert_eq!(comparisons, EXPECTED_CANDIDATE_COMPARISONS);
        assert_eq!(current.to_string(), EXPECTED_RESIDUE_SMALL);
        let residue = ResidueReport {
            mode: "nested-composite/x-only".into(),
            primary_oracle_queries: 72,
            candidate_comparisons: comparisons,
            residue: EXPECTED_RESIDUE_SMALL.into(),
            negative_residue: EXPECTED_RESIDUE_LARGE.into(),
            expected_unordered_residues_match: true,
            transcript_sha256: hex::encode(sha256(&transcript)),
            stages,
            operations: NativeOpCounts::default(),
        };

        let plus_residue = decimal(EXPECTED_RESIDUE_SMALL);
        let minus_residue = decimal(EXPECTED_RESIDUE_LARGE);
        let solution = IntervalSolution {
            orientation: -1,
            k: EXPECTED_RESIDUAL_K,
            giant_index: 3_159_226,
            baby_index: 989_698,
            baby_sign: -1,
        };
        let mut baby_operations = NativeOpCounts::default();
        baby_operations.point_additions = 1;
        let mut giant_operations = NativeOpCounts::default();
        giant_operations.point_additions = 1;
        let mut total_operations = baby_operations;
        total_operations.merge(giant_operations);
        let stats = IntervalStats {
            table_bits: FROZEN_TABLE_BITS,
            table_bytes: 8u64 << FROZEN_TABLE_BITS,
            baby_entries: FROZEN_M - 1,
            baby_chain_steps: FROZEN_M - 1,
            giant_positions_plus: solution.giant_index + 1,
            giant_positions_minus: solution.giant_index + 1,
            fingerprint_hits: 1,
            false_fingerprint_hits: 0,
            exact_replays: 1,
            baby_operations,
            giant_and_replay_operations: giant_operations,
            operations: total_operations,
        };
        let target = super::super::p192_native::JacobianPoint::from_affine(public, &arithmetic);
        arithmetic.reset_counts();
        assert!(accepts_recovered_target(&secret, target, &mut arithmetic));
        let final_replay_operations = arithmetic.counts;
        let interval = IntervalCertificate {
            geometry: SHIFTED_GEOMETRY.into(),
            shifted_centers: true,
            plan: IntervalPlan::frozen(4, 256).unwrap(),
            plus_residue: plus_residue.to_string(),
            minus_residue: minus_residue.to_string(),
            plus_width: interval_width(&plus_residue).unwrap(),
            minus_width: interval_width(&minus_residue).unwrap(),
            solution,
            stats,
            recovered_d: secret.to_string(),
            recovered_d_hex: format!("{secret:048x}"),
            final_replay_matches_public_key: true,
            final_replay_operations,
            frozen_tail_covered: true,
        };
        let p59 = p59_control();
        let controls = ControlReport {
            exhaustive_p59_sign_classes: p59.passed,
            p59_scalars_tested: p59.scalars_tested,
            p59_nested_queries: p59.nested_queries,
            p59_candidate_comparisons: p59.candidate_comparisons,
            p59_two_sign_stage_classes: p59.two_sign_stage_classes,
            p59_self_inverse_stage_classes: p59.self_inverse_stage_classes,
            challenge59_regression: true,
            safe_ecdh_rejected_every_probe: true,
            wrong_a_rejected: true,
            mutated_singular_order_rejected: true,
            mutated_factor_rejected: true,
            mutated_target_rejected: true,
            one_bit_tag_mutation_rejected: true,
            raw_full_point_model_executed: true,
            raw_full_point_model_recovered_mod_m: true,
            raw_full_point_model_primary: false,
        };
        let mut report = RecoveryReport {
            schema: REPORT_SCHEMA.into(),
            experiment_id: EXPERIMENT_ID.into(),
            seed_hex: FROZEN_SEED_HEX.into(),
            seed_digest_hex: hex::encode(digest),
            planted_d: secret.to_string(),
            planted_d_hex: format!("{secret:048x}"),
            public_key: point_hex(public, &arithmetic),
            certification: certify().unwrap(),
            residue,
            interval,
            controls,
            success: true,
            certificate_sha256: String::new(),
        };
        report.certificate_sha256 = report_digest(&report).unwrap();
        report
    }

    fn recompute_unkeyed_digest(report: &mut RecoveryReport) {
        report.certificate_sha256 = report_digest(report).unwrap();
    }

    fn assert_semantic_tamper_rejected<F>(base: &RecoveryReport, mutate: F)
    where
        F: FnOnce(&mut RecoveryReport),
    {
        let mut tampered = base.clone();
        mutate(&mut tampered);
        recompute_unkeyed_digest(&mut tampered);
        assert!(verify_report_preflight(&tampered).is_err());
    }

    #[test]
    fn frozen_seed_known_answer() {
        let seed = hex::decode(FROZEN_SEED_HEX).unwrap();
        let (digest, scalar) = derive_planted_key(&seed);
        assert_eq!(hex::encode(digest), EXPECTED_DIGEST_HEX);
        assert_eq!(scalar, decimal(EXPECTED_D_DEC));
        assert_eq!(format!("{scalar:048x}"), EXPECTED_D_HEX);
    }

    #[test]
    fn pocklington_chain_proves_p_and_n() {
        let certificates = pocklington_certificates().unwrap();
        assert_eq!(certificates.len(), 5);
        assert!(certificates
            .iter()
            .all(|certificate| { certificate.factorization_exact && certificate.verified_prime }));
        assert_eq!(certificates[3].value, p192_prime().to_string());
        assert_eq!(certificates[4].value, p192_order().to_string());
    }

    #[test]
    fn actual_low_level_x_only_oracle_matches_torus() {
        let seed = hex::decode(FROZEN_SEED_HEX).unwrap();
        let (_, secret) = derive_planted_key(&seed);
        let mut arithmetic = P192Arithmetic::new();
        let m = BigUint::from(2u8 * 2 * 2 * 3);
        let probe_torus = TorusPoint::from_t(8, &arithmetic)
            .scalar_mul(&(singular_order() / &m), &mut arithmetic);
        let probe = probe_torus.to_curve(&mut arithmetic).unwrap();
        let (tag, actual_zx) =
            victim_tag(&CurveParams::p192(), probe, &secret, 7, &arithmetic).unwrap();
        let independent = probe_torus.scalar_mul(&secret, &mut arithmetic);
        let independent_point = independent.to_curve(&mut arithmetic);
        let expected_zx = zx_from_native(independent_point, &arithmetic);
        assert_eq!(actual_zx, expected_zx);
        assert_eq!(
            tag,
            offline_tag(
                independent_point,
                &confirmation_message(7, probe, &arithmetic),
                &arithmetic
            )
        );
        assert!(safe_ecdh_rejects(
            &CurveParams::p192(),
            probe,
            &secret,
            &arithmetic
        ));
    }

    #[test]
    fn reduced_nested_schedule_recovers_sign_class() {
        let secret = BigUint::from(12_345u64);
        let schedule = [2u64, 2, 2, 2, 2, 3, 5, 17];
        let report = recover_residue_schedule(&secret, &schedule, 16).unwrap();
        let modulus = decimal(&report.stages.last().unwrap().modulus_after);
        let residue = decimal(&report.residue);
        let negative = decimal(&report.negative_residue);
        let d_mod = &secret % &modulus;
        assert!(residue == d_mod || negative == d_mod);
        assert_eq!(report.primary_oracle_queries, schedule.len() as u64);
        assert!(report.stages.iter().all(|stage| stage.probe_exact_order
            && stage.safe_ecdh_rejected
            && stage.low_level_matches_independent_torus));
    }

    #[test]
    fn exhaustive_p59_control_accounts_for_every_sign_collision() {
        let outcome = p59_control();
        assert!(outcome.passed);
        assert_eq!(outcome.scalars_tested, 60);
        assert_eq!(outcome.nested_queries, 240);
        assert_eq!(
            outcome.two_sign_stage_classes + outcome.self_inverse_stage_classes,
            240
        );
        assert!(outcome.candidate_comparisons >= outcome.nested_queries);
        assert!(outcome.self_inverse_stage_classes > 0);
        assert!(outcome.two_sign_stage_classes > 0);
    }

    #[test]
    fn reduced_full_point_positive_model_resolves_the_sign() {
        let secret = BigUint::from(12_345u64);
        let schedule = [2u64, 2, 2, 2, 2, 3, 5, 17];
        let (recovered, modulus, queries) =
            recover_full_point_model_schedule(&secret, &schedule, 16).unwrap();
        assert_eq!(queries, schedule.len() as u64);
        assert_eq!(recovered, &secret % &modulus);
    }

    #[test]
    fn shifted_centers_cover_frozen_tail_and_range_guard() {
        assert!(frozen_tail_covered());
        let positions = QMAX.div_ceil(FROZEN_STRIDE);
        assert_eq!(positions, 5_800_018);
        let last_index = positions - 1;
        let center = (FROZEN_M - 1) + last_index * FROZEN_STRIDE;
        let tail = QMAX - 1;
        assert!(center.abs_diff(tail) < FROZEN_M);
        // The last shifted block extends beyond qmax; exact-candidate range
        // checks must reject this first out-of-range scalar.
        let first_out_of_range = QMAX;
        assert!(first_out_of_range > tail);
        assert!(center.abs_diff(first_out_of_range) < FROZEN_M);
    }

    #[test]
    fn recomputed_unkeyed_digest_cannot_hide_semantic_tampering() {
        let report = synthetic_preflight_report();
        verify_report_preflight(&report).unwrap();

        assert_semantic_tamper_rejected(&report, |value| {
            value.interval.geometry.push_str("-tampered")
        });
        assert_semantic_tamper_rejected(&report, |value| value.interval.plus_residue = "1".into());
        assert_semantic_tamper_rejected(&report, |value| value.interval.plus_width -= 1);
        assert_semantic_tamper_rejected(&report, |value| value.interval.solution.baby_index += 1);
        assert_semantic_tamper_rejected(&report, |value| {
            value.interval.stats.giant_positions_minus -= 1
        });
        assert_semantic_tamper_rejected(&report, |value| {
            value.interval.stats.operations.point_additions += 1
        });
        assert_semantic_tamper_rejected(&report, |value| {
            value.interval.final_replay_matches_public_key = false
        });
        assert_semantic_tamper_rejected(&report, |value| {
            value.controls.mutated_target_rejected = false
        });
        assert_semantic_tamper_rejected(&report, |value| {
            let replacement = if value.public_key.x.starts_with('0') {
                "1"
            } else {
                "0"
            };
            value.public_key.x.replace_range(..1, replacement);
        });
    }

    #[test]
    fn negative_controls_exercise_real_acceptance_predicates() {
        let factors = singular_factorization();
        assert!(accepts_singular_factorization(&singular_order(), &factors));
        assert!(!accepts_singular_factorization(
            &(singular_order() + BigUint::one()),
            &factors
        ));
        let mut mutated_factors = factors;
        mutated_factors[0].1 -= 1;
        assert!(!accepts_singular_factorization(
            &singular_order(),
            &mutated_factors
        ));

        let expected = [0x5au8; 32];
        let mut flipped = expected;
        flipped[7] ^= 1;
        assert!(accepts_confirmation_tag(&expected, &expected));
        assert!(!accepts_confirmation_tag(&expected, &flipped));

        let seed = hex::decode(FROZEN_SEED_HEX).unwrap();
        let (_, secret) = derive_planted_key(&seed);
        let mut arithmetic = P192Arithmetic::new();
        let g = AffinePoint::generator(&arithmetic);
        let q = scalar_mul(g, &secret, &mut arithmetic);
        assert!(accepts_recovered_target(&secret, q, &mut arithmetic));
        let q_affine = q.to_affine(&mut arithmetic).unwrap();
        let mutated_q =
            super::super::p192_native::JacobianPoint::from_affine(q_affine, &arithmetic)
                .add_affine(g, &mut arithmetic);
        assert!(!accepts_recovered_target(
            &secret,
            mutated_q,
            &mut arithmetic
        ));
    }
}
