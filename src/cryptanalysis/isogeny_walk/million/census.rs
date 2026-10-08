//! Streaming, source-bound structural traits for certified P-256 grid curves.
//!
//! The source certificates already carry the expensive kernel and order
//! replay evidence. This sidecar keeps shared isogeny-class values once and
//! derives only model-, path-, and edge-specific values per curve. Grid axes
//! and isogeny-volcano directions are deliberately separate fields.

use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use flate2::read::GzDecoder;
use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use super::*;
use crate::cryptanalysis::isogeny_walk::field::is_probable_prime;
use crate::cryptanalysis::p256_isogeny_campaign::P256_ICV1;
use crate::cryptanalysis::p256_structural;

pub const TRAIT_CENSUS_SCHEMA: &str = "p256.isogeny-grid-trait-census/v1";
pub const TRAIT_CENSUS_SUMMARY_SCHEMA: &str = "p256.isogeny-grid-trait-census-summary/v1";
pub const TRAIT_CENSUS_RECEIPT_SCHEMA: &str = "p256.isogeny-grid-trait-census-receipt/v1";

const FACTOR_DISCOVERY_SOURCE: &str = "https://factordb.com/";
const FACTOR_VERIFICATION_METHOD: &str =
    "native exact product plus deterministic-u64 Miller-Rabin and recursive Pocklington n-1 certificates";

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FactorPower {
    pub prime: String,
    pub exponent: u32,
    pub prime_status: String,
    pub proof: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PocklingtonWitness {
    pub prime_divisor: String,
    pub witness: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimeCertificate {
    pub prime: String,
    pub method: String,
    pub n_minus_one_factors: Vec<FactorPower>,
    pub known_factor_product: String,
    pub witnesses: Vec<PocklingtonWitness>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Factorization {
    pub value: String,
    pub factors: Vec<FactorPower>,
    pub factor_product: String,
    pub remaining_cofactor: String,
    pub status: String,
    pub verification_method: String,
    pub discovery_source: String,
    pub largest_prime_factor: String,
    pub largest_prime_factor_bits: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EllClassTraits {
    pub ell: u64,
    pub frobenius_discriminant_mod_ell: u64,
    pub splitting: String,
    pub v_ell_frobenius_discriminant: u32,
    pub volcano_depth: u32,
    pub ell_divides_frobenius_conductor: bool,
    pub expected_fp_rational_neighbours: u32,
    pub direction_for_certified_edges: String,
    pub direction_status: String,
    pub reason: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExtensionOrder {
    pub degree: u32,
    pub order: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClassTraits {
    pub field_p: String,
    pub field_prime_status: String,
    pub field_prime_certificate: String,
    pub group_order: String,
    pub group_order_factorization: Factorization,
    pub group_structure: String,
    pub trace: String,
    pub hasse_bound_satisfied: bool,
    pub ordinary: bool,
    pub rational_torsion: String,
    pub frobenius_polynomial: String,
    pub frobenius_discriminant: String,
    pub frobenius_discriminant_abs_bits: u64,
    pub frobenius_discriminant_factorization: Factorization,
    pub cm_field_discriminant: String,
    pub cm_field_discriminant_status: String,
    pub frobenius_order_conductor: String,
    pub frobenius_order_conductor_factorization: Factorization,
    pub endomorphism_order_conductor_all_curves: String,
    pub endomorphism_order_status: String,
    pub conductor_gap: String,
    pub twist_order: String,
    pub twist_order_factorization: Factorization,
    pub extension_orders: Vec<ExtensionOrder>,
    pub embedding_degree: Option<u64>,
    pub embedding_degree_searched_to: u64,
    pub ell: Vec<EllClassTraits>,
    pub prime_certificates: Vec<PrimeCertificate>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceBinding {
    pub source_index: u32,
    pub schema: String,
    pub source_commit: String,
    pub target_curves: u64,
    pub curve_records: u64,
    pub bound_records: u64,
    pub record_chain_sha256: String,
    pub stored_sha256: String,
    pub stored_bytes: u64,
    pub content_sha256: String,
    pub content_bytes: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TraitCensusHeader {
    pub record: String,
    pub schema: String,
    pub source_commit: String,
    pub input_curves: u64,
    pub unique_input_j_invariants: u64,
    pub sources: Vec<SourceBinding>,
    pub class: ClassTraits,
    pub completeness: String,
    pub integer_factor_scope: String,
    pub path_degree_encoding: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct QuadraticCharacters {
    pub a: i8,
    pub b: i8,
    pub j: i8,
    pub model_denominator: i8,
    pub c4: i8,
    pub c6: i8,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelTraits {
    pub model_denominator: String,
    pub equation_discriminant: String,
    pub c4: String,
    pub c6: String,
    pub j_recomputed_equal: bool,
    pub quadratic_characters: QuadraticCharacters,
    pub j_zero: bool,
    pub j_1728: bool,
    pub automorphism_order_over_algebraic_closure: u32,
    pub canonical_model_branch: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PathTraits {
    pub degree_notation: String,
    pub degree_11_exponent: u32,
    pub degree_13_exponent: u32,
    pub degree_bit_length: u16,
    pub decimal_digits: u16,
    pub log10_lower_inclusive: u16,
    pub log10_upper_exclusive: u16,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EdgeTraits {
    pub source_coordinate: Coordinate,
    pub target_coordinate: Coordinate,
    pub source_curve_uid: String,
    pub target_curve_uid: String,
    pub degree: u64,
    pub construction_axis: String,
    pub predecessor_j: Option<String>,
    pub defined_over: String,
    pub separable: bool,
    pub cyclic: bool,
    pub splitting: String,
    pub ell_divides_frobenius_conductor: bool,
    pub volcano_direction: String,
    pub direction_status: String,
    pub direction_reason: String,
    pub source_endomorphism_conductor: Option<String>,
    pub target_endomorphism_conductor: Option<String>,
    pub conductor_ratio: Option<String>,
    pub conductor_status: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TraitCurveRecord {
    pub record: String,
    pub ordinal: u64,
    pub source_index: u32,
    pub source_stored_sha256: String,
    pub source_record_index: u64,
    pub coordinate: Coordinate,
    pub a: String,
    pub b: String,
    pub j: String,
    pub icv1_slug: String,
    pub icv1_identity: String,
    pub ec1_alias: String,
    pub curve_uid: String,
    pub order_status: String,
    pub detectors: DetectorRecord,
    pub model: ModelTraits,
    pub path: PathTraits,
    pub incoming_edge: Option<EdgeTraits>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CharacterCounts {
    pub non_residue: u64,
    pub zero: u64,
    pub residue: u64,
}

impl CharacterCounts {
    fn observe(&mut self, value: i8) -> Result<(), String> {
        match value {
            -1 => self.non_residue += 1,
            0 => self.zero += 1,
            1 => self.residue += 1,
            _ => return Err(format!("invalid quadratic character {value}")),
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CharacterSummary {
    pub a: CharacterCounts,
    pub b: CharacterCounts,
    pub j: CharacterCounts,
    pub model_denominator: CharacterCounts,
    pub c4: CharacterCounts,
    pub c6: CharacterCounts,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PathExtremum {
    pub coordinate: Coordinate,
    pub bit_length: u16,
    pub decimal_digits: u16,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TraitCensusSummary {
    pub record: String,
    pub schema: String,
    pub bound_records: u64,
    pub record_chain_sha256: String,
    pub curves: u64,
    pub edges: u64,
    pub row_degree_11_edges: u64,
    pub spine_degree_13_edges: u64,
    pub horizontal_edges: u64,
    pub ascending_edges: u64,
    pub descending_edges: u64,
    pub unknown_direction_edges: u64,
    pub endomorphism_conductor_one_curves: u64,
    pub unresolved_conductor_curves: u64,
    pub registered_root_models: u64,
    pub a_minus_3_models: u64,
    pub general_canonical_models: u64,
    pub j_zero_curves: u64,
    pub j_1728_curves: u64,
    pub quadratic_characters: CharacterSummary,
    pub qr_prefix_64_histogram: Vec<HistogramBin>,
    pub minimum_a_signed_bits: u64,
    pub minimum_b_signed_bits: u64,
    pub maximum_path_degree: PathExtremum,
    pub discrepancies: Vec<String>,
    pub ecdlp_speedup: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TraitCensusReceipt {
    pub schema: String,
    pub status: String,
    pub source_commit: String,
    pub elapsed_ms: u128,
    pub resources: ResourceUsage,
    pub uncompressed_bytes: u64,
    pub sources: Vec<SourceBinding>,
    pub summary: TraitCensusSummary,
}

fn number(value: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.as_bytes(), 10).ok_or_else(|| format!("invalid integer {value}"))
}

fn factor(prime: &str, exponent: u32, proof: &str) -> FactorPower {
    FactorPower {
        prime: prime.into(),
        exponent,
        prime_status: "proved_prime".into(),
        proof: proof.into(),
    }
}

fn mul_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(modulus)) as u64
}

fn pow_mod_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    base %= modulus;
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = mul_mod_u64(result, base, modulus);
        }
        base = mul_mod_u64(base, base, modulus);
        exponent >>= 1;
    }
    result
}

/// Deterministic Miller-Rabin for every `u64`.
fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for prime in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == prime {
            return true;
        }
        if n.is_multiple_of(prime) {
            return false;
        }
    }
    let mut d = n - 1;
    let s = d.trailing_zeros();
    d >>= s;
    for witness in [2u64, 325, 9_375, 28_178, 450_775, 9_780_504, 1_795_265_022] {
        let a = witness % n;
        if a == 0 {
            continue;
        }
        let mut x = pow_mod_u64(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut passed = false;
        for _ in 1..s {
            x = mul_mod_u64(x, x, n);
            if x == n - 1 {
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

fn factor_product(factors: &[FactorPower]) -> Result<BigUint, String> {
    factors.iter().try_fold(BigUint::one(), |product, factor| {
        Ok(product * number(&factor.prime)?.pow(factor.exponent))
    })
}

fn pocklington_certificate(
    prime: &str,
    raw_factors: &[(&str, u32)],
    proved_large: &HashSet<String>,
) -> Result<PrimeCertificate, String> {
    let n = number(prime)?;
    if n <= BigUint::from(u64::MAX) {
        return Err(format!(
            "{prime}: Pocklington is reserved for primes above u64"
        ));
    }
    if !is_probable_prime(&n) {
        return Err(format!(
            "{prime}: fixed-base Miller-Rabin rejected candidate"
        ));
    }
    let factors: Vec<FactorPower> = raw_factors
        .iter()
        .map(|(q, exponent)| {
            let value = number(q)?;
            let proof = if let Some(small) = value.to_u64() {
                if !is_prime_u64(small) {
                    return Err(format!("{prime}: n-1 factor {q} is composite"));
                }
                "deterministic_miller_rabin_u64".to_string()
            } else {
                if !proved_large.contains(*q) {
                    return Err(format!("{prime}: n-1 factor {q} has no prior proof"));
                }
                format!("pocklington:{q}")
            };
            Ok(factor(q, *exponent, &proof))
        })
        .collect::<Result<_, String>>()?;
    let product = factor_product(&factors)?;
    if product != &n - 1u8 {
        return Err(format!("{prime}: n-1 factorization product differs"));
    }
    let nm1 = &n - 1u8;
    let mut witnesses = Vec::with_capacity(factors.len());
    for item in &factors {
        let q = number(&item.prime)?;
        let exponent = &nm1 / &q;
        let mut found = None;
        for candidate in 2u64..=65_537 {
            let a = BigUint::from(candidate);
            if a.modpow(&nm1, &n) != BigUint::one() {
                continue;
            }
            let power = a.modpow(&exponent, &n);
            let difference = if power.is_zero() {
                &n - 1u8
            } else {
                power - 1u8
            };
            if difference.gcd(&n).is_one() {
                found = Some(candidate);
                break;
            }
        }
        let witness = found.ok_or_else(|| {
            format!(
                "{prime}: no Pocklington witness found for factor {}",
                item.prime
            )
        })?;
        witnesses.push(PocklingtonWitness {
            prime_divisor: item.prime.clone(),
            witness: witness.to_string(),
        });
    }
    Ok(PrimeCertificate {
        prime: prime.into(),
        method: "Pocklington theorem with complete proved factorization of n-1".into(),
        n_minus_one_factors: factors,
        known_factor_product: product.to_string(),
        witnesses,
    })
}

fn build_prime_certificates() -> Result<Vec<PrimeCertificate>, String> {
    // The factors were discovered after protocol freeze. The native code does
    // not trust that source: it checks every product and supplies a recursive
    // Pocklington certificate down to deterministic-u64 primes.
    let specs: Vec<(&str, Vec<(&str, u32)>)> = vec![
        (
            "1869236796843064056413",
            vec![("2", 2), ("23", 1), ("3911", 1), ("5195037399650551", 1)],
        ),
        (
            "112600039977335252189018616651067202354555503111433",
            vec![
                ("2", 3),
                ("857", 1),
                ("131909", 1),
                ("1387801", 1),
                ("124778574733515931", 1),
                ("718995365596683143", 1),
            ],
        ),
        (
            "2624747550333869278416773953",
            vec![
                ("2", 6),
                ("3", 2),
                ("1297", 1),
                ("16879", 1),
                ("208150935158385979", 1),
            ],
        ),
        (
            "774023187263532362759620327192479577272145303",
            vec![
                ("2", 1),
                ("3", 2),
                ("2411", 1),
                ("34282281433", 1),
                ("11290956913871", 1),
                ("46076956964474543", 1),
            ],
        ),
        (
            "835945042244614951780389953367877943453916927241",
            vec![
                ("2", 3),
                ("3", 3),
                ("5", 1),
                ("774023187263532362759620327192479577272145303", 1),
            ],
        ),
        (
            "1428624589419343516204097",
            vec![("2", 6), ("1669", 1), ("13374631042347059581", 1)],
        ),
        (
            "46523541035814968339936406074986559003387",
            vec![
                ("2", 1),
                ("269", 1),
                ("643", 1),
                ("215531", 1),
                ("333814693", 1),
                ("1869236796843064056413", 1),
            ],
        ),
        (
            "3317349640749355357762425066592395746459685764401801118712075735758936647",
            vec![
                ("2", 1),
                ("3", 1),
                ("19", 1),
                ("1861", 1),
                ("2729", 1),
                ("250259", 1),
                ("203333173", 1),
                ("112600039977335252189018616651067202354555503111433", 1),
            ],
        ),
        (
            "115792089210356248762697446949407573529996955224135760342422259061068512044369",
            vec![
                ("2", 4),
                ("3", 1),
                ("71", 1),
                ("131", 1),
                ("373", 1),
                ("3407", 1),
                ("17449", 1),
                ("38189", 1),
                ("187019741", 1),
                ("622491383", 1),
                ("1002328039319", 1),
                ("2624747550333869278416773953", 1),
            ],
        ),
        (
            "115792089210356248762697446949407573530086143415290314195533631308867097853951",
            vec![
                ("2", 1),
                ("3", 1),
                ("5", 2),
                ("17", 1),
                ("257", 1),
                ("641", 1),
                ("1531", 1),
                ("65537", 1),
                ("490463", 1),
                ("6700417", 1),
                ("835945042244614951780389953367877943453916927241", 1),
            ],
        ),
    ];
    let mut proved = HashSet::new();
    let mut certificates = Vec::with_capacity(specs.len());
    for (prime, factors) in specs {
        let certificate = pocklington_certificate(prime, &factors, &proved)?;
        proved.insert(prime.to_string());
        certificates.push(certificate);
    }
    Ok(certificates)
}

fn proved_factorization(
    value: &BigUint,
    raw_factors: &[(&str, u32)],
    certificates: &[PrimeCertificate],
) -> Result<Factorization, String> {
    let proved_large: HashSet<&str> = certificates.iter().map(|c| c.prime.as_str()).collect();
    let factors: Vec<FactorPower> = raw_factors
        .iter()
        .map(|(prime, exponent)| {
            let n = number(prime)?;
            let proof = if let Some(small) = n.to_u64() {
                if !is_prime_u64(small) {
                    return Err(format!("factor {prime} is not prime"));
                }
                "deterministic_miller_rabin_u64".to_string()
            } else if proved_large.contains(prime) {
                format!("pocklington:{prime}")
            } else {
                return Err(format!("factor {prime} has no primality certificate"));
            };
            Ok(factor(prime, *exponent, &proof))
        })
        .collect::<Result<_, String>>()?;
    let product = factor_product(&factors)?;
    if &product != value {
        return Err(format!(
            "factorization product {product} differs from {value}"
        ));
    }
    let largest = factors
        .iter()
        .map(|entry| number(&entry.prime).map(|value| (value, entry.prime.clone())))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .max_by(|left, right| left.0.cmp(&right.0))
        .ok_or("empty factorization")?;
    Ok(Factorization {
        value: value.to_string(),
        factors,
        factor_product: product.to_string(),
        remaining_cofactor: "1".into(),
        status: "complete_proved".into(),
        verification_method: FACTOR_VERIFICATION_METHOD.into(),
        discovery_source: FACTOR_DISCOVERY_SOURCE.into(),
        largest_prime_factor: largest.1,
        largest_prime_factor_bits: largest.0.bits(),
    })
}

fn unit_factorization() -> Factorization {
    Factorization {
        value: "1".into(),
        factors: Vec::new(),
        factor_product: "1".into(),
        remaining_cofactor: "1".into(),
        status: "complete_proved".into(),
        verification_method: "empty product".into(),
        discovery_source: "derived from fundamental Frobenius discriminant".into(),
        largest_prime_factor: "1".into(),
        largest_prime_factor_bits: 1,
    }
}

fn class_traits() -> Result<ClassTraits, String> {
    let start = StartCurve::p256();
    let p = start.p.clone();
    let n = start.order.clone();
    let p_int = BigInt::from(p.clone());
    let trace = &p_int + 1u8 - BigInt::from(n.clone());
    let disc = &trace * &trace - &p_int * 4u8;
    if disc.sign() != num_bigint::Sign::Minus {
        return Err("P-256 Frobenius discriminant is not negative".into());
    }
    let disc_abs = disc.magnitude().clone();
    let twist = (&p << 1) + 2u8 - &n;
    let certificates = build_prime_certificates()?;
    let field_factorization = proved_factorization(
        &p,
        &[(
            "115792089210356248762697446949407573530086143415290314195533631308867097853951",
            1,
        )],
        &certificates,
    )?;
    let order_factorization = proved_factorization(
        &n,
        &[(
            "115792089210356248762697446949407573529996955224135760342422259061068512044369",
            1,
        )],
        &certificates,
    )?;
    let disc_factorization = proved_factorization(
        &disc_abs,
        &[
            ("3", 1),
            ("5", 1),
            ("456597257999", 1),
            ("1428624589419343516204097", 1),
            ("46523541035814968339936406074986559003387", 1),
        ],
        &certificates,
    )?;
    let twist_factorization = proved_factorization(
        &twist,
        &[
            ("3", 1),
            ("5", 1),
            ("13", 1),
            ("179", 1),
            (
                "3317349640749355357762425066592395746459685764401801118712075735758936647",
                1,
            ),
        ],
        &certificates,
    )?;
    if disc.mod_floor(&BigInt::from(4u8)) != BigInt::one()
        || disc_factorization
            .factors
            .iter()
            .any(|factor| factor.exponent != 1)
    {
        return Err("P-256 Frobenius discriminant is not fundamental".into());
    }
    let conductor_factorization = unit_factorization();
    let class = ClassInfo::compute(&start, &[ROW_DEGREE, SPINE_DEGREE]);
    let ell = class
        .ells
        .iter()
        .map(|(degree, kind, valuation)| {
            if *valuation != 0 || !matches!(kind, EllKind::Elkies) {
                return Err(format!(
                    "degree {degree}: expected split depth-zero P-256 graph"
                ));
            }
            let residue = disc
                .mod_floor(&BigInt::from(*degree))
                .to_u64()
                .ok_or("discriminant residue does not fit u64")?;
            Ok(EllClassTraits {
                ell: *degree,
                frobenius_discriminant_mod_ell: residue,
                splitting: "elkies_split".into(),
                v_ell_frobenius_discriminant: *valuation,
                volcano_depth: 0,
                ell_divides_frobenius_conductor: false,
                expected_fp_rational_neighbours: 2,
                direction_for_certified_edges: "horizontal".into(),
                direction_status: "proved".into(),
                reason: "v_ell(Delta_pi)=0, so every endomorphism order is ell-maximal; split Frobenius gives two horizontal F_p-rational neighbours".into(),
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let hasse_bound_satisfied = &trace * &trace <= &p_int * 4u8;
    let ordinary = trace.gcd(&p_int).is_one();
    let extension_orders = [2, 3, 4, 6]
        .into_iter()
        .map(|degree| ExtensionOrder {
            degree,
            order: p256_structural::extension_field_order_generic(&p, &trace, degree).to_string(),
        })
        .collect();
    Ok(ClassTraits {
        field_p: p.to_string(),
        field_prime_status: "proved_prime".into(),
        field_prime_certificate: format!("pocklington:{}", field_factorization.value),
        group_order: n.to_string(),
        group_order_factorization: order_factorization,
        group_structure: "cyclic_prime_order".into(),
        trace: trace.to_string(),
        hasse_bound_satisfied,
        ordinary,
        rational_torsion: "no nontrivial rational torsion of order below the prime group order"
            .into(),
        frobenius_polynomial: format!("X^2-({trace})*X+{p}"),
        frobenius_discriminant: disc.to_string(),
        frobenius_discriminant_abs_bits: disc_abs.bits(),
        frobenius_discriminant_factorization: disc_factorization,
        cm_field_discriminant: disc.to_string(),
        cm_field_discriminant_status: "proved_fundamental".into(),
        frobenius_order_conductor: "1".into(),
        frobenius_order_conductor_factorization: conductor_factorization,
        endomorphism_order_conductor_all_curves: "1".into(),
        endomorphism_order_status: "proved_maximal_for_entire_fp_isogeny_class".into(),
        conductor_gap: "none".into(),
        twist_order: twist.to_string(),
        twist_order_factorization: twist_factorization,
        extension_orders,
        embedding_degree: class.embedding_degree,
        embedding_degree_searched_to: super::super::walk::EMBEDDING_SEARCH,
        ell,
        prime_certificates: certificates,
    })
}

fn source_bindings(receipt: &JUnionReceipt) -> Vec<SourceBinding> {
    receipt
        .inputs
        .iter()
        .enumerate()
        .map(|(source_index, input)| SourceBinding {
            source_index: source_index as u32,
            schema: input.schema.clone(),
            source_commit: input.source_commit.clone(),
            target_curves: input.target_curves,
            curve_records: input.curve_records,
            bound_records: input.bound_records,
            record_chain_sha256: input.record_chain_sha256.clone(),
            stored_sha256: input.stored_sha256.clone(),
            stored_bytes: input.stored_bytes,
            content_sha256: input.content_sha256.clone(),
            content_bytes: input.content_bytes,
        })
        .collect()
}

fn open_lines(path: &Path) -> Result<Box<dyn BufRead>, String> {
    let file = File::open(path).map_err(|error| format!("open {}: {error}", path.display()))?;
    if path.extension().and_then(|value| value.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

#[derive(Clone, Copy)]
struct SourceShape {
    schema: &'static str,
    side_or_width: u32,
    y_start: u32,
    y_end: u32,
}

fn validate_p256_source_header(
    path: &Path,
    value: &Value,
    expected_target_curves: u64,
) -> Result<(), String> {
    let start = StartCurve::p256();
    let trace = BigInt::from(start.p.clone()) + 1u8 - BigInt::from(start.order.clone());
    let expected_p = format!("{:064x}", start.p);
    let expected_order = format!("{:064x}", start.order);
    let context = format!("{} header", path.display());
    let require_text = |field: &str, expected: &str| -> Result<(), String> {
        if value[field].as_str() != Some(expected) {
            return Err(format!(
                "{context}: {field} differs from frozen P-256 value"
            ));
        }
        Ok(())
    };
    let require_u64 = |field: &str, expected: u64| -> Result<(), String> {
        if value[field].as_u64() != Some(expected) {
            return Err(format!(
                "{context}: {field} differs from frozen P-256 value"
            ));
        }
        Ok(())
    };

    require_text("record", "header")?;
    require_text("root_name", &start.name)?;
    require_text("root_icv1_slug", P256_ICV1)?;
    require_text("field_p", &expected_p)?;
    require_text("group_order", &expected_order)?;
    require_text("trace", &trace.to_string())?;
    require_text("canonical_model_rule", curve::CANONICAL_RULE)?;
    require_text("generator_rule", curve::GENERATOR_RULE)?;
    require_text(
        "first_edge_rule",
        "choose the numerically smaller rational root",
    )?;
    require_text(
        "later_edge_rule",
        "choose the rational root unequal to the predecessor j",
    )?;
    require_u64("target_curves", expected_target_curves)?;
    require_u64("spine_degree", SPINE_DEGREE)?;
    require_u64("row_degree", ROW_DEGREE)?;
    require_u64("root_seed", ROOT_SEED)?;
    require_u64("construction_audit_points", 1)?;
    require_u64("construction_audit_seed_x", 0)?;

    let source_commit = value["source_commit"]
        .as_str()
        .ok_or_else(|| format!("{context}: missing source_commit"))?;
    if source_commit.len() != 40 || !source_commit.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(format!(
            "{context}: source_commit is not 40 hexadecimal digits"
        ));
    }
    for field in ["spine_phi_checks", "row_phi_checks"] {
        if value[field].as_u64().is_none_or(|checks| checks == 0) {
            return Err(format!("{context}: {field} is missing or zero"));
        }
    }
    Ok(())
}

fn source_shape(path: &Path) -> Result<SourceShape, String> {
    let mut reader = open_lines(path)?;
    let mut line = String::new();
    reader
        .read_line(&mut line)
        .map_err(|error| format!("read {}: {error}", path.display()))?;
    let value: Value = serde_json::from_str(line.trim_end())
        .map_err(|error| format!("parse {} header: {error}", path.display()))?;
    let schema = value["schema"]
        .as_str()
        .ok_or_else(|| format!("{} header has no schema", path.display()))?;
    match schema {
        SCHEMA => {
            let header: HeaderRecord = serde_json::from_value(value.clone())
                .map_err(|error| format!("parse {} grid header: {error}", path.display()))?;
            if header.side < 2 {
                return Err(format!("{} grid side is less than two", path.display()));
            }
            let target_curves = u64::from(header.side) * u64::from(header.side);
            validate_p256_source_header(path, &value, target_curves)?;
            Ok(SourceShape {
                schema: SCHEMA,
                side_or_width: header.side,
                y_start: 0,
                y_end: header.side,
            })
        }
        STRIP_SCHEMA => {
            let header: StripHeaderRecord = serde_json::from_value(value.clone())
                .map_err(|error| format!("parse {} strip header: {error}", path.display()))?;
            if header.width < 2 || header.height == 0 {
                return Err(format!(
                    "{} strip width is less than two or height is zero",
                    path.display()
                ));
            }
            let target_curves = u64::from(header.width) * u64::from(header.height);
            validate_p256_source_header(path, &value, target_curves)?;
            let expected_context_records = u64::from(header.y_start != 0);
            if header.context_records != expected_context_records {
                return Err(format!(
                    "{} strip context record count differs",
                    path.display()
                ));
            }
            Ok(SourceShape {
                schema: STRIP_SCHEMA,
                side_or_width: header.width,
                y_start: header.y_start,
                y_end: header
                    .y_start
                    .checked_add(header.height)
                    .ok_or("strip y range overflow")?,
            })
        }
        other => Err(format!("{}: unsupported schema {other}", path.display())),
    }
}

struct PathMetricsTable {
    width: usize,
    values: Vec<(u16, u16)>,
}

impl PathMetricsTable {
    fn build(max_x: u32, max_y: u32) -> Result<Self, String> {
        let width = max_x as usize + 1;
        let len = width
            .checked_mul(max_y as usize + 1)
            .ok_or("path metric table size overflow")?;
        let mut values = Vec::with_capacity(len);
        let mut power_13 = BigUint::one();
        let mut y_digits = 1u16;
        let mut y_next_power_10 = BigUint::from(10u8);
        for _y in 0..=max_y {
            let mut degree = power_13.clone();
            let mut digits = y_digits;
            let mut next_power_10 = y_next_power_10.clone();
            for _x in 0..=max_x {
                values.push((
                    degree
                        .bits()
                        .try_into()
                        .map_err(|_| "path degree bit length exceeds u16")?,
                    digits,
                ));
                degree *= ROW_DEGREE;
                while degree >= next_power_10 {
                    next_power_10 *= 10u8;
                    digits = digits.checked_add(1).ok_or("decimal digits exceed u16")?;
                }
            }
            power_13 *= SPINE_DEGREE;
            while power_13 >= y_next_power_10 {
                y_next_power_10 *= 10u8;
                y_digits = y_digits.checked_add(1).ok_or("decimal digits exceed u16")?;
            }
        }
        Ok(Self { width, values })
    }

    fn get(&self, coordinate: Coordinate) -> Result<PathTraits, String> {
        let index = coordinate.y as usize * self.width + coordinate.x as usize;
        let (bits, digits) = *self
            .values
            .get(index)
            .ok_or_else(|| format!("no path metrics for ({},{})", coordinate.x, coordinate.y))?;
        Ok(PathTraits {
            degree_notation: format!(
                "{}^{}*{}^{}",
                ROW_DEGREE, coordinate.x, SPINE_DEGREE, coordinate.y
            ),
            degree_11_exponent: coordinate.x,
            degree_13_exponent: coordinate.y,
            degree_bit_length: bits,
            decimal_digits: digits,
            log10_lower_inclusive: digits.saturating_sub(1),
            log10_upper_exclusive: digits,
        })
    }
}

fn model_traits(
    field: &Field,
    model: &Model,
    recorded_j: &str,
    branch: &str,
) -> Result<ModelTraits, String> {
    let denominator = model.discriminant(field);
    let equation_discriminant = field.neg(&field.mul(&field.from_u64(16), &denominator));
    let c4 = field.mul(&field.from_i64(-48), &model.a);
    let c6 = field.mul(&field.from_i64(-864), &model.b);
    let j = model
        .j(field)
        .ok_or("census encountered a singular model")?;
    let zero = field.zero();
    let j_1728_value = field.from_u64(1728);
    let j_zero = j == zero;
    let j_1728 = j == j_1728_value;
    Ok(ModelTraits {
        model_denominator: hex_fe(field, &denominator),
        equation_discriminant: hex_fe(field, &equation_discriminant),
        c4: hex_fe(field, &c4),
        c6: hex_fe(field, &c6),
        j_recomputed_equal: hex_fe(field, &j) == recorded_j,
        quadratic_characters: QuadraticCharacters {
            a: field.legendre(&model.a) as i8,
            b: field.legendre(&model.b) as i8,
            j: field.legendre(&j) as i8,
            model_denominator: field.legendre(&denominator) as i8,
            c4: field.legendre(&c4) as i8,
            c6: field.legendre(&c6) as i8,
        },
        j_zero,
        j_1728,
        automorphism_order_over_algebraic_closure: if j_zero {
            6
        } else if j_1728 {
            4
        } else {
            2
        },
        canonical_model_branch: branch.into(),
    })
}

#[derive(Default)]
struct CensusAccumulator {
    curves: u64,
    edges: u64,
    row_edges: u64,
    spine_edges: u64,
    horizontal: u64,
    ascending: u64,
    descending: u64,
    unknown: u64,
    conductor_one: u64,
    unresolved_conductor: u64,
    registered_root: u64,
    a_minus_3: u64,
    general: u64,
    j_zero: u64,
    j_1728: u64,
    characters: CharacterSummary,
    qr_histogram: Vec<u64>,
    min_a: u64,
    min_b: u64,
    maximum_path: Option<PathExtremum>,
}

fn validate_edge_traits(edge: &EdgeTraits) -> Result<(), String> {
    match (edge.degree, edge.construction_axis.as_str()) {
        (ROW_DEGREE, "row_degree_11") | (SPINE_DEGREE, "spine_degree_13") => {}
        _ => return Err("edge degree and construction axis differ".into()),
    }
    if edge.direction_reason.trim().is_empty() {
        return Err("edge direction requires a nonempty reason".into());
    }
    if edge.conductor_status.trim().is_empty() {
        return Err("edge conductor status must not be empty".into());
    }
    match edge.volcano_direction.as_str() {
        "horizontal" | "ascending" | "descending" => {
            if edge.direction_status != "proved"
                || edge.source_endomorphism_conductor.is_none()
                || edge.target_endomorphism_conductor.is_none()
                || edge.conductor_ratio.is_none()
            {
                return Err(
                    "proved edge direction requires proved status and both conductors".into(),
                );
            }
        }
        "unknown" => {
            if edge.direction_status == "proved" {
                return Err("unknown edge direction cannot have proved status".into());
            }
        }
        other => return Err(format!("unknown volcano direction {other}")),
    }
    Ok(())
}

impl CensusAccumulator {
    fn new() -> Self {
        Self {
            min_a: u64::MAX,
            min_b: u64::MAX,
            qr_histogram: vec![0; 65],
            ..Self::default()
        }
    }

    fn observe(&mut self, record: &TraitCurveRecord) -> Result<(), String> {
        self.curves += 1;
        if let Some(edge) = &record.incoming_edge {
            validate_edge_traits(edge)?;
            if edge.target_coordinate != record.coordinate
                || edge.target_curve_uid != record.curve_uid
            {
                return Err("incoming edge target differs from its curve record".into());
            }
            self.edges += 1;
            match edge.construction_axis.as_str() {
                "row_degree_11" => self.row_edges += 1,
                "spine_degree_13" => self.spine_edges += 1,
                other => return Err(format!("unknown construction axis {other}")),
            }
            match edge.volcano_direction.as_str() {
                "horizontal" => self.horizontal += 1,
                "ascending" => self.ascending += 1,
                "descending" => self.descending += 1,
                "unknown" => self.unknown += 1,
                _ => unreachable!("validated above"),
            }
        }
        if record.incoming_edge.as_ref().is_none_or(|edge| {
            edge.source_endomorphism_conductor.as_deref() == Some("1")
                && edge.target_endomorphism_conductor.as_deref() == Some("1")
        }) {
            self.conductor_one += 1;
        } else {
            self.unresolved_conductor += 1;
        }
        match record.model.canonical_model_branch.as_str() {
            "registered_root" => self.registered_root += 1,
            "a_minus_3" => self.a_minus_3 += 1,
            "general_j_model" => self.general += 1,
            other => return Err(format!("unknown canonical model branch {other}")),
        }
        self.j_zero += u64::from(record.model.j_zero);
        self.j_1728 += u64::from(record.model.j_1728);
        let q = &record.model.quadratic_characters;
        self.characters.a.observe(q.a)?;
        self.characters.b.observe(q.b)?;
        self.characters.j.observe(q.j)?;
        self.characters
            .model_denominator
            .observe(q.model_denominator)?;
        self.characters.c4.observe(q.c4)?;
        self.characters.c6.observe(q.c6)?;
        let qr = record.detectors.qr_prefix_64 as usize;
        if qr >= self.qr_histogram.len() {
            return Err(format!("qr_prefix_64 {qr} exceeds 64"));
        }
        self.qr_histogram[qr] += 1;
        self.min_a = self.min_a.min(record.detectors.a_signed_bits);
        self.min_b = self.min_b.min(record.detectors.b_signed_bits);
        if self.maximum_path.as_ref().is_none_or(|old| {
            (record.path.degree_bit_length, record.path.decimal_digits)
                > (old.bit_length, old.decimal_digits)
        }) {
            self.maximum_path = Some(PathExtremum {
                coordinate: record.coordinate,
                bit_length: record.path.degree_bit_length,
                decimal_digits: record.path.decimal_digits,
            });
        }
        Ok(())
    }

    fn summary(&self, chain: [u8; 32]) -> Result<TraitCensusSummary, String> {
        if self.edges
            != self
                .horizontal
                .checked_add(self.ascending)
                .and_then(|value| value.checked_add(self.descending))
                .and_then(|value| value.checked_add(self.unknown))
                .ok_or("direction count overflow")?
        {
            return Err("direction buckets do not account for every edge".into());
        }
        let qr_prefix_64_histogram = self
            .qr_histogram
            .iter()
            .enumerate()
            .filter(|(_, count)| **count != 0)
            .map(|(value, count)| HistogramBin {
                value: value as u32,
                count: *count,
            })
            .collect();
        Ok(TraitCensusSummary {
            record: "summary".into(),
            schema: TRAIT_CENSUS_SUMMARY_SCHEMA.into(),
            bound_records: self.curves + 1,
            record_chain_sha256: hex::encode(chain),
            curves: self.curves,
            edges: self.edges,
            row_degree_11_edges: self.row_edges,
            spine_degree_13_edges: self.spine_edges,
            horizontal_edges: self.horizontal,
            ascending_edges: self.ascending,
            descending_edges: self.descending,
            unknown_direction_edges: self.unknown,
            endomorphism_conductor_one_curves: self.conductor_one,
            unresolved_conductor_curves: self.unresolved_conductor,
            registered_root_models: self.registered_root,
            a_minus_3_models: self.a_minus_3,
            general_canonical_models: self.general,
            j_zero_curves: self.j_zero,
            j_1728_curves: self.j_1728,
            quadratic_characters: self.characters.clone(),
            qr_prefix_64_histogram,
            minimum_a_signed_bits: self.min_a,
            minimum_b_signed_bits: self.min_b,
            maximum_path_degree: self.maximum_path.clone().ok_or("census has no curves")?,
            discrepancies: Vec::new(),
            ecdlp_speedup: None,
        })
    }
}

fn value_u64(value: &Value, field: &str, context: &str) -> Result<u64, String> {
    value[field]
        .as_u64()
        .ok_or_else(|| format!("{context}: missing {field}"))
}

fn derive_records<W: Write>(
    output: &mut LineWriter<'_, W>,
    source_paths: &[PathBuf],
    bindings: &[SourceBinding],
    class: &ClassTraits,
    path_metrics: &PathMetricsTable,
) -> Result<CensusAccumulator, String> {
    let start = StartCurve::p256();
    let field = Field::new(&start.p).ok_or("P-256 field construction failed")?;
    let minus_three = field.from_i64(-3);
    let mut accumulator = CensusAccumulator::new();
    let mut ordinal = 0u64;

    for (source_index, (path, binding)) in source_paths.iter().zip(bindings).enumerate() {
        let shape = source_shape(path)?;
        if shape.schema != binding.schema {
            return Err(format!(
                "{} schema changed after union audit",
                path.display()
            ));
        }
        let mut reader = open_lines(path)?;
        let mut line_number = 0u64;
        let mut emitted = 0u64;
        let mut spine_uids = BTreeMap::<u32, String>::new();
        let mut last_row_coordinate = None;
        let mut last_row_uid = None::<String>;
        loop {
            let mut line = String::new();
            let bytes = reader
                .read_line(&mut line)
                .map_err(|error| format!("read {}: {error}", path.display()))?;
            if bytes == 0 {
                break;
            }
            line_number += 1;
            let context = format!("{} line {line_number}", path.display());
            let value: Value = serde_json::from_str(line.trim_end())
                .map_err(|error| format!("{context}: {error}"))?;
            let kind = value["record"]
                .as_str()
                .ok_or_else(|| format!("{context}: missing record type"))?;
            if matches!(kind, "header" | "summary") {
                continue;
            }
            let coordinate: Coordinate = serde_json::from_value(value["coordinate"].clone())
                .map_err(|error| format!("{context}: coordinate: {error}"))?;
            let curve: CurveRecord = serde_json::from_value(value["curve"].clone())
                .map_err(|error| format!("{context}: curve: {error}"))?;
            if kind == "boundary" {
                spine_uids.insert(coordinate.y, curve.curve_uid);
                continue;
            }
            if kind != "root" && kind != "node" {
                return Err(format!("{context}: unsupported record type {kind}"));
            }
            if coordinate.x >= shape.side_or_width
                || coordinate.y < shape.y_start
                || coordinate.y >= shape.y_end
            {
                return Err(format!(
                    "{context}: coordinate ({},{}) lies outside header bounds",
                    coordinate.x, coordinate.y
                ));
            }
            let source_record_index = value_u64(&value, "index", &context)?;
            let (model, recorded_j) = parse_model(&field, &curve, &context)?;
            let branch = if kind == "root" {
                "registered_root"
            } else if model.a == minus_three {
                "a_minus_3"
            } else {
                "general_j_model"
            };
            let model_traits = model_traits(&field, &model, &curve.j, branch)?;
            if !model_traits.j_recomputed_equal {
                return Err(format!("{context}: recorded j differs from model"));
            }
            let incoming_edge = if kind == "root" {
                if coordinate != (Coordinate { x: 0, y: 0 }) {
                    return Err(format!("{context}: root is not at (0,0)"));
                }
                None
            } else {
                let degree = value_u64(&value, "degree", &context)?;
                let source_coordinate = if coordinate.x == 0 {
                    Coordinate {
                        x: 0,
                        y: coordinate.y.checked_sub(1).ok_or("spine underflow")?,
                    }
                } else {
                    Coordinate {
                        x: coordinate.x - 1,
                        y: coordinate.y,
                    }
                };
                let source_uid = if coordinate.x == 0 {
                    spine_uids.get(&source_coordinate.y).cloned()
                } else if coordinate.x == 1 {
                    spine_uids.get(&coordinate.y).cloned()
                } else if last_row_coordinate == Some(source_coordinate) {
                    last_row_uid.clone()
                } else {
                    None
                }
                .ok_or_else(|| {
                    format!(
                        "{context}: source identity for ({},{}) is unavailable",
                        source_coordinate.x, source_coordinate.y
                    )
                })?;
                let expected_degree = if coordinate.x == 0 {
                    SPINE_DEGREE
                } else {
                    ROW_DEGREE
                };
                if degree != expected_degree {
                    return Err(format!(
                        "{context}: degree {degree} differs from expected {expected_degree}"
                    ));
                }
                let expected_parent = if shape.schema == SCHEMA {
                    node_parent(shape.side_or_width, coordinate.x, coordinate.y)
                } else {
                    (u64::from(source_coordinate.y) << 32) | u64::from(source_coordinate.x)
                };
                if value_u64(&value, "parent", &context)? != expected_parent {
                    return Err(format!("{context}: parent index differs"));
                }
                let ell = class
                    .ell
                    .iter()
                    .find(|traits| traits.ell == degree)
                    .ok_or_else(|| format!("{context}: no class traits for degree {degree}"))?;
                Some(EdgeTraits {
                    source_coordinate,
                    target_coordinate: coordinate,
                    source_curve_uid: source_uid,
                    target_curve_uid: curve.curve_uid.clone(),
                    degree,
                    construction_axis: if degree == ROW_DEGREE {
                        "row_degree_11".into()
                    } else {
                        "spine_degree_13".into()
                    },
                    predecessor_j: value["predecessor_j"].as_str().map(str::to_owned),
                    defined_over: "F_p".into(),
                    separable: true,
                    cyclic: true,
                    splitting: ell.splitting.clone(),
                    ell_divides_frobenius_conductor: ell.ell_divides_frobenius_conductor,
                    volcano_direction: ell.direction_for_certified_edges.clone(),
                    direction_status: ell.direction_status.clone(),
                    direction_reason: ell.reason.clone(),
                    source_endomorphism_conductor: Some("1".into()),
                    target_endomorphism_conductor: Some("1".into()),
                    conductor_ratio: Some("1".into()),
                    conductor_status: "proved_from_fundamental_frobenius_discriminant".into(),
                })
            };
            let record = TraitCurveRecord {
                record: "curve".into(),
                ordinal,
                source_index: source_index as u32,
                source_stored_sha256: binding.stored_sha256.clone(),
                source_record_index,
                coordinate,
                a: curve.a.clone(),
                b: curve.b.clone(),
                j: hex_fe(&field, &recorded_j),
                icv1_slug: curve.icv1_slug.clone(),
                icv1_identity: curve.icv1_identity.clone(),
                ec1_alias: curve.ec1_alias.clone(),
                curve_uid: curve.curve_uid.clone(),
                order_status: curve.order_status.clone(),
                detectors: curve.detectors.clone(),
                model: model_traits,
                path: path_metrics.get(coordinate)?,
                incoming_edge,
            };
            accumulator.observe(&record)?;
            output.write(&record, true)?;
            if coordinate.x == 0 {
                spine_uids.insert(coordinate.y, curve.curve_uid.clone());
            } else {
                last_row_coordinate = Some(coordinate);
                last_row_uid = Some(curve.curve_uid.clone());
            }
            emitted += 1;
            ordinal += 1;
        }
        if emitted != binding.curve_records {
            return Err(format!(
                "{} emitted {emitted} curves, union audit bound {}",
                path.display(),
                binding.curve_records
            ));
        }
    }
    Ok(accumulator)
}

pub fn generate_trait_census_jsonl<W: Write>(
    writer: &mut W,
    source_paths: &[PathBuf],
    source_commit: String,
) -> Result<TraitCensusReceipt, String> {
    if source_commit.len() != 40 || !source_commit.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err("source commit must be a 40-digit hexadecimal Git object id".into());
    }
    let started = Instant::now();
    let union = audit_j_union_paths(source_paths)?;
    let bindings = source_bindings(&union);
    let shapes = source_paths
        .iter()
        .map(|path| source_shape(path))
        .collect::<Result<Vec<_>, _>>()?;
    let max_x = shapes
        .iter()
        .map(|shape| shape.side_or_width - 1)
        .max()
        .ok_or("trait census requires at least one source")?;
    let max_y = shapes
        .iter()
        .map(|shape| shape.y_end - 1)
        .max()
        .ok_or("trait census requires at least one source")?;
    let path_metrics = PathMetricsTable::build(max_x, max_y)?;
    let class = class_traits()?;
    let header = TraitCensusHeader {
        record: "header".into(),
        schema: TRAIT_CENSUS_SCHEMA.into(),
        source_commit: source_commit.clone(),
        input_curves: union.input_curves,
        unique_input_j_invariants: union.unique_j_invariants,
        sources: bindings.clone(),
        class: class.clone(),
        completeness: "complete_under_p256.isogeny-grid-trait-census/v1; unknowns require typed status and reason".into(),
        integer_factor_scope: "class integers only; canonical field-element representatives are not integer-factorized".into(),
        path_degree_encoding: "exact lossless factorization 11^x*13^y plus exact bit length and decimal-digit interval".into(),
    };
    let mut output = LineWriter::new(writer);
    output.write(&header, true)?;
    let accumulator = derive_records(&mut output, source_paths, &bindings, &class, &path_metrics)?;
    if accumulator.curves != union.input_curves {
        return Err(format!(
            "derived {} curves, union audit reported {}",
            accumulator.curves, union.input_curves
        ));
    }
    let summary = accumulator.summary(output.chain)?;
    output.write(&summary, false)?;
    output
        .writer
        .flush()
        .map_err(|error| format!("flush trait census: {error}"))?;
    Ok(TraitCensusReceipt {
        schema: TRAIT_CENSUS_RECEIPT_SCHEMA.into(),
        status: "pass".into(),
        source_commit,
        elapsed_ms: started.elapsed().as_millis(),
        resources: resource_usage(),
        uncompressed_bytes: output.bytes,
        sources: bindings,
        summary,
    })
}

struct CompareWriter<R: Read> {
    expected: R,
    compared: u64,
}

impl<R: Read> Write for CompareWriter<R> {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        let mut expected = vec![0u8; buffer.len()];
        self.expected.read_exact(&mut expected)?;
        if expected != buffer {
            let offset = expected
                .iter()
                .zip(buffer)
                .position(|(left, right)| left != right)
                .unwrap_or(0);
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "trait census differs from replay at uncompressed byte {}",
                    self.compared + offset as u64
                ),
            ));
        }
        self.compared += buffer.len() as u64;
        Ok(buffer.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn census_reader(path: &Path) -> Result<Box<dyn Read>, String> {
    let file = File::open(path).map_err(|error| format!("open {}: {error}", path.display()))?;
    if path.extension().and_then(|value| value.to_str()) == Some("gz") {
        Ok(Box::new(GzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

fn census_source_commit(path: &Path) -> Result<String, String> {
    let mut reader = BufReader::new(census_reader(path)?);
    let mut line = String::new();
    reader
        .read_line(&mut line)
        .map_err(|error| format!("read {}: {error}", path.display()))?;
    let header: TraitCensusHeader = serde_json::from_str(line.trim_end())
        .map_err(|error| format!("parse {} header: {error}", path.display()))?;
    if header.record != "header" || header.schema != TRAIT_CENSUS_SCHEMA {
        return Err(format!("{} is not a trait census", path.display()));
    }
    Ok(header.source_commit)
}

pub fn verify_trait_census_paths(
    census_path: &Path,
    source_paths: &[PathBuf],
) -> Result<TraitCensusReceipt, String> {
    let started = Instant::now();
    let source_commit = census_source_commit(census_path)?;
    let mut compare = CompareWriter {
        expected: census_reader(census_path)?,
        compared: 0,
    };
    let mut receipt = generate_trait_census_jsonl(&mut compare, source_paths, source_commit)?;
    let mut tail = [0u8; 1];
    if compare
        .expected
        .read(&mut tail)
        .map_err(|error| format!("read trait census tail: {error}"))?
        != 0
    {
        return Err("trait census has trailing bytes after replay".into());
    }
    receipt.elapsed_ms = started.elapsed().as_millis();
    receipt.resources = resource_usage();
    receipt.status = "pass_independent_source_replay".into();
    Ok(receipt)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static SERIAL: AtomicU64 = AtomicU64::new(0);

    fn temporary(stem: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "p256-trait-census-{stem}-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn p256_factor_certificates_prove_maximal_cm_order() {
        let class = class_traits().unwrap();
        assert_eq!(class.cm_field_discriminant_status, "proved_fundamental");
        assert_eq!(class.frobenius_order_conductor, "1");
        assert_eq!(class.endomorphism_order_conductor_all_curves, "1");
        assert!(class
            .ell
            .iter()
            .all(|ell| ell.direction_for_certified_edges == "horizontal"));
        assert_eq!(class.frobenius_discriminant_factorization.factors.len(), 5);
        assert_eq!(class.twist_order_factorization.status, "complete_proved");
    }

    #[test]
    fn path_metrics_are_exact_on_small_grid() {
        let table = PathMetricsTable::build(7, 7).unwrap();
        for y in 0..=7 {
            for x in 0..=7 {
                let exact = BigUint::from(ROW_DEGREE).pow(x) * BigUint::from(SPINE_DEGREE).pow(y);
                let metric = table.get(Coordinate { x, y }).unwrap();
                assert_eq!(u64::from(metric.degree_bit_length), exact.bits());
                assert_eq!(
                    usize::from(metric.decimal_digits),
                    exact.to_str_radix(10).len()
                );
            }
        }
    }

    #[test]
    fn deterministic_u64_primality_rejects_pseudoprimes() {
        assert!(is_prime_u64(18_446_744_073_709_551_557));
        assert!(!is_prime_u64(341_550_071_728_321));
        assert!(!is_prime_u64(3_825_123_056_546_413_051));
    }

    #[test]
    fn unresolved_direction_requires_status_reason_and_complete_schema() {
        let edge = EdgeTraits {
            source_coordinate: Coordinate { x: 0, y: 0 },
            target_coordinate: Coordinate { x: 1, y: 0 },
            source_curve_uid: "source".into(),
            target_curve_uid: "target".into(),
            degree: ROW_DEGREE,
            construction_axis: "row_degree_11".into(),
            predecessor_j: None,
            defined_over: "F_p".into(),
            separable: true,
            cyclic: true,
            splitting: "unresolved".into(),
            ell_divides_frobenius_conductor: false,
            volcano_direction: "unknown".into(),
            direction_status: "unresolved".into(),
            direction_reason: "conductor evidence is unavailable".into(),
            source_endomorphism_conductor: None,
            target_endomorphism_conductor: None,
            conductor_ratio: None,
            conductor_status: "unresolved: no conductor certificate".into(),
        };
        validate_edge_traits(&edge).unwrap();

        let mut missing_reason = edge.clone();
        missing_reason.direction_reason.clear();
        assert!(validate_edge_traits(&missing_reason).is_err());

        let mut omitted = serde_json::to_value(&edge).unwrap();
        omitted.as_object_mut().unwrap().remove("direction_reason");
        assert!(serde_json::from_value::<EdgeTraits>(omitted).is_err());
    }

    #[test]
    fn adjacent_sources_replay_and_tampering_or_overlap_is_rejected() {
        let legacy_path = temporary("legacy.jsonl");
        let strip_path = temporary("strip.jsonl");
        let census_path = temporary("census.jsonl");
        let tampered_path = temporary("tampered.jsonl");
        let tampered_source_path = temporary("tampered-source.jsonl");
        let source_commit = "1".repeat(40);

        let mut legacy = Vec::new();
        generate_jsonl(
            &mut legacy,
            &GridConfig {
                side: 3,
                batch_rows: 2,
                audit_points: 1,
                audit_seed_x: 0,
                source_commit: source_commit.clone(),
            },
        )
        .unwrap();
        fs::write(&legacy_path, legacy).unwrap();

        let mut strip = Vec::new();
        generate_strip_jsonl(
            &mut strip,
            &StripConfig {
                width: 3,
                y_start: 3,
                height: 1,
                batch_rows: 1,
                audit_points: 1,
                audit_seed_x: 0,
                source_commit: source_commit.clone(),
            },
        )
        .unwrap();
        fs::write(&strip_path, strip).unwrap();

        let sources = vec![legacy_path.clone(), strip_path.clone()];
        let mut census = Vec::new();
        let generated = generate_trait_census_jsonl(&mut census, &sources, source_commit).unwrap();
        assert_eq!(generated.summary.curves, 12);
        assert_eq!(generated.summary.edges, 11);
        assert_eq!(generated.summary.horizontal_edges, 11);
        assert_eq!(generated.summary.ascending_edges, 0);
        assert_eq!(generated.summary.descending_edges, 0);
        assert_eq!(generated.summary.unknown_direction_edges, 0);
        assert_eq!(generated.summary.endomorphism_conductor_one_curves, 12);
        fs::write(&census_path, &census).unwrap();
        let replayed = verify_trait_census_paths(&census_path, &sources).unwrap();
        assert_eq!(replayed.summary, generated.summary);

        let needle = b"\"volcano_direction\":\"horizontal\"";
        let position = census
            .windows(needle.len())
            .position(|window| window == needle)
            .unwrap();
        let value_offset = position + b"\"volcano_direction\":\"".len();
        census[value_offset] = b'H';
        fs::write(&tampered_path, census).unwrap();
        assert!(verify_trait_census_paths(&tampered_path, &sources).is_err());

        let mut tampered_source = fs::read(&strip_path).unwrap();
        let source_needle = b"\"root_name\":\"P-256\"";
        let source_position = tampered_source
            .windows(source_needle.len())
            .position(|window| window == source_needle)
            .unwrap();
        tampered_source[source_position + b"\"root_name\":\"".len()] = b'Q';
        fs::write(&tampered_source_path, tampered_source).unwrap();
        assert!(verify_trait_census_paths(
            &census_path,
            &[legacy_path.clone(), tampered_source_path.clone()]
        )
        .is_err());

        let mut rejected = Vec::new();
        assert!(generate_trait_census_jsonl(
            &mut rejected,
            &[legacy_path.clone(), legacy_path.clone()],
            "2".repeat(40),
        )
        .is_err());

        for path in [
            legacy_path,
            strip_path,
            census_path,
            tampered_path,
            tampered_source_path,
        ] {
            fs::remove_file(path).unwrap();
        }
    }
}
