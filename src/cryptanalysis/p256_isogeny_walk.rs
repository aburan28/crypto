//! Deterministic, certificate-producing degree-11 walk in P-256's
//! `F_p`-isogeny class.
//!
//! The production-scale campaign remains out of scope. This module implements
//! the bounded Gate-1 experiment frozen in
//! `research/p256_isogeny_walk_20261004/PROTOCOL.md`: classify small degrees,
//! walk `2^12` degree-11 edges through the classical modular polynomial, and
//! disambiguate every target twist with an exact-order point witness.

use crate::cryptanalysis::p256_isogeny_campaign::P256_ICV1;
use crate::ecc::{curve::CurveParams, point::Point};
use crate::hash::sha256::sha256;
use crate::prime_hyperelliptic::FpPoly;
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

pub const CERTIFICATE_SCHEMA: &str = "p256.isogeny-walk-certificate/v1";
pub const PHI11_SOURCE: &str = "https://math.mit.edu/~drew/modpolys/jfiles/phi_j_11.txt";
pub const PHI11_SHA256: &str = "985c6a76ca7cc2802ba3f9fd48112fe0d82ffe11dbaa5f49c57526658338b352";
pub const WALK_DEGREE: u64 = 11;
pub const DEFAULT_STEPS: u64 = 1 << 12;
const PHI11_TEXT: &str = include_str!("data/phi_j_11.txt");

const SMALL_PRIMES: [u64; 11] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DegreeCensus {
    pub degree: u64,
    pub discriminant_mod_degree: u64,
    pub class: String,
    pub fp_edges: u8,
    pub non_backtracking_walk: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PointWitness {
    pub x: String,
    pub y: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CurveRecord {
    pub a: String,
    pub b: String,
    pub j: String,
    pub order_witness: PointWitness,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WalkHeader {
    pub record: String,
    pub schema: String,
    pub curve: String,
    pub field_p: String,
    pub group_order: String,
    pub trace: String,
    pub degree: u64,
    pub requested_steps: u64,
    pub phi11_source: String,
    pub phi11_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CensusRecord {
    pub record: String,
    pub entries: Vec<DegreeCensus>,
    pub selected_degree: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StartRecord {
    pub record: String,
    pub curve: CurveRecord,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IsogenyWitness {
    pub kind: String,
    pub degree: u64,
    pub phi11_sha256: String,
    pub phi_value: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EdgeRecord {
    pub record: String,
    pub index: u64,
    pub predecessor_j: Option<String>,
    pub source_j: String,
    pub target: CurveRecord,
    pub direction: String,
    pub isogeny: IsogenyWitness,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CycleEvent {
    pub at_vertex: u64,
    pub first_seen_vertex: u64,
    pub period: u64,
    pub j: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SummaryRecord {
    pub record: String,
    pub emitted_steps: u64,
    pub verified_steps: u64,
    pub distinct_j: u64,
    pub first_cycle: Option<CycleEvent>,
    pub generation_ms: u128,
    pub edges_per_second: f64,
    pub serialized_edge_bytes: u64,
    pub bytes_per_edge: f64,
    pub edge_stream_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WalkCertificate {
    pub header: WalkHeader,
    pub census: CensusRecord,
    pub start: StartRecord,
    pub edges: Vec<EdgeRecord>,
    pub summary: SummaryRecord,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DivisionPolynomialFixture {
    pub scalar: u64,
    pub kernel_coefficients: Vec<String>,
    pub codomain_a: String,
    pub codomain_b: String,
    pub codomain_j: String,
}

#[derive(Clone, Debug)]
struct Phi11 {
    terms: Vec<(usize, usize, BigInt)>,
}

impl Phi11 {
    fn load() -> Result<Self, String> {
        let actual = hex::encode(sha256(PHI11_TEXT.as_bytes()));
        if actual != PHI11_SHA256 {
            return Err(format!(
                "Phi_11 source digest changed: expected {PHI11_SHA256}, got {actual}"
            ));
        }
        let mut terms = Vec::new();
        for (line_number, line) in PHI11_TEXT.lines().enumerate() {
            let (indices, coefficient) = line
                .split_once(' ')
                .ok_or_else(|| format!("malformed Phi_11 line {}", line_number + 1))?;
            let indices = indices
                .strip_prefix('[')
                .and_then(|value| value.strip_suffix(']'))
                .ok_or_else(|| format!("malformed Phi_11 indices on line {}", line_number + 1))?;
            let (i, j) = indices
                .split_once(',')
                .ok_or_else(|| format!("malformed Phi_11 indices on line {}", line_number + 1))?;
            let i = i
                .parse::<usize>()
                .map_err(|error| format!("invalid Phi_11 exponent: {error}"))?;
            let j = j
                .parse::<usize>()
                .map_err(|error| format!("invalid Phi_11 exponent: {error}"))?;
            let coefficient = BigInt::parse_bytes(coefficient.as_bytes(), 10)
                .ok_or_else(|| format!("invalid Phi_11 coefficient on line {}", line_number + 1))?;
            terms.push((i, j, coefficient));
        }
        if terms.len() != 79
            || !terms
                .iter()
                .any(|(i, j, c)| *i == 12 && *j == 0 && c.is_one())
        {
            return Err("Phi_11 table does not have the expected 79-term monic shape".into());
        }
        Ok(Self { terms })
    }

    fn specialize(&self, j_invariant: &BigUint, p: &BigUint) -> FpPoly {
        let mut powers = Vec::with_capacity(13);
        powers.push(BigUint::one());
        for exponent in 1..=12 {
            powers.push((&powers[exponent - 1] * j_invariant) % p);
        }
        let mut coefficients = vec![BigUint::zero(); 13];
        for (i, j, coefficient) in &self.terms {
            let coefficient = signed_mod(coefficient, p);
            add_mod_assign(&mut coefficients[*j], &(&coefficient * &powers[*i] % p), p);
            if i != j {
                add_mod_assign(&mut coefficients[*i], &(&coefficient * &powers[*j] % p), p);
            }
        }
        FpPoly::from_coeffs(coefficients, p.clone())
    }

    fn evaluate(&self, source_j: &BigUint, target_j: &BigUint, p: &BigUint) -> BigUint {
        self.specialize(source_j, p).eval(target_j)
    }
}

fn signed_mod(value: &BigInt, p: &BigUint) -> BigUint {
    let residue = value.magnitude() % p;
    if value.sign() == Sign::Minus && !residue.is_zero() {
        p - residue
    } else {
        residue
    }
}

fn add_mod_assign(target: &mut BigUint, value: &BigUint, p: &BigUint) {
    *target = (&*target + value) % p;
}

fn parse_hex(name: &str, value: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.as_bytes(), 16).ok_or_else(|| format!("invalid hexadecimal {name}"))
}

fn hex_big(value: &BigUint) -> String {
    format!("{value:064x}")
}

fn p256_trace(curve: &CurveParams) -> BigUint {
    &curve.p + BigUint::one() - &curve.n
}

fn mod_u64(value: &BigUint, modulus: u64) -> u64 {
    (value % BigUint::from(modulus))
        .to_u64()
        .expect("remainder fits u64")
}

fn mul_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

fn pow_mod_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut accumulator = 1 % modulus;
    while exponent != 0 {
        if exponent & 1 == 1 {
            accumulator = mul_mod_u64(accumulator, base, modulus);
        }
        base = mul_mod_u64(base, base, modulus);
        exponent >>= 1;
    }
    accumulator
}

pub fn degree_census() -> Vec<DegreeCensus> {
    let curve = CurveParams::p256();
    let trace = p256_trace(&curve);
    SMALL_PRIMES
        .into_iter()
        .map(|degree| {
            if degree == 2 {
                return DegreeCensus {
                    degree,
                    discriminant_mod_degree: 1,
                    class: "no-rational-2-torsion".into(),
                    fp_edges: 0,
                    non_backtracking_walk: false,
                };
            }
            let t = mod_u64(&trace, degree);
            let p = mod_u64(&curve.p, degree);
            let discriminant =
                (mul_mod_u64(t, t, degree) + degree - mul_mod_u64(4 % degree, p, degree)) % degree;
            let symbol = if discriminant == 0 {
                0
            } else {
                pow_mod_u64(discriminant, (degree - 1) / 2, degree)
            };
            let (class, fp_edges, non_backtracking_walk) = match symbol {
                0 => ("ramified-single-edge", 1, false),
                1 => ("split-two-edges", 2, true),
                _ => ("inert-no-edge", 0, false),
            };
            DegreeCensus {
                degree,
                discriminant_mod_degree: discriminant,
                class: class.into(),
                fp_edges,
                non_backtracking_walk,
            }
        })
        .collect()
}

pub fn selected_walk_degree(census: &[DegreeCensus]) -> Option<u64> {
    census
        .iter()
        .find(|entry| entry.non_backtracking_walk)
        .map(|entry| entry.degree)
}

fn modular_roots(phi: &Phi11, j_invariant: &BigUint, p: &BigUint) -> Result<[BigUint; 2], String> {
    let polynomial = phi.specialize(j_invariant, p).monic();
    if polynomial.degree() != Some(12) {
        return Err("Phi_11 specialization is not degree 12".into());
    }
    let x = FpPoly::x(p.clone());
    let rational_part = polynomial.gcd(&x.pow_mod(p, &polynomial).sub(&x));
    if rational_part.degree() != Some(2) {
        return Err(format!(
            "Phi_11 specialization has {} F_p-root degree, expected 2",
            rational_part.degree().unwrap_or(0)
        ));
    }
    let rational_part = rational_part.monic();
    let linear = rational_part.coeff(1);
    let constant = rational_part.coeff(0);
    let four_constant = (BigUint::from(4u32) * constant) % p;
    let discriminant = ((&linear * &linear) % p + p - four_constant) % p;
    let square_root = sqrt_p256(&discriminant, p)
        .ok_or_else(|| "quadratic F_p-root factor has nonsquare discriminant".to_string())?;
    let inverse_two = (p + BigUint::one()) >> 1usize;
    let negative_linear = if linear.is_zero() {
        BigUint::zero()
    } else {
        p - linear
    };
    let mut roots = [
        ((&negative_linear + &square_root) * &inverse_two) % p,
        ((&negative_linear + p - &square_root) * &inverse_two) % p,
    ];
    roots.sort();
    if roots[0] == roots[1]
        || !polynomial.eval(&roots[0]).is_zero()
        || !polynomial.eval(&roots[1]).is_zero()
    {
        return Err("Phi_11 roots are not two distinct verified field elements".into());
    }
    Ok(roots)
}

fn sqrt_p256(value: &BigUint, p: &BigUint) -> Option<BigUint> {
    let value = value % p;
    if value.is_zero() {
        return Some(BigUint::zero());
    }
    if (p % BigUint::from(4u32)) != BigUint::from(3u32) {
        return None;
    }
    let root = value.modpow(&((p + BigUint::one()) >> 2usize), p);
    ((&root * &root) % p == value).then_some(root)
}

fn curve_j(a: &BigUint, b: &BigUint, p: &BigUint) -> Option<BigUint> {
    let a3 = a.modpow(&BigUint::from(3u32), p);
    let denominator = (BigUint::from(4u32) * &a3 + BigUint::from(27u32) * b * b) % p;
    if denominator.is_zero() {
        return None;
    }
    let inverse = denominator.modpow(&(p - BigUint::from(2u32)), p);
    Some((BigUint::from(6912u32) * a3 * inverse) % p)
}

fn as_curve(
    a: &BigUint,
    b: &BigUint,
    witness: &PointWitness,
) -> Result<(CurveParams, Point), String> {
    let p256 = CurveParams::p256();
    let x = parse_hex("witness x", &witness.x)?;
    let y = parse_hex("witness y", &witness.y)?;
    let curve = CurveParams {
        name: "P-256-isogenous",
        p: p256.p,
        a: a.clone(),
        b: b.clone(),
        gx: x.clone(),
        gy: y.clone(),
        n: p256.n,
        h: 1,
    };
    let point = Point::Affine {
        x: curve.fe(x),
        y: curve.fe(y),
    };
    Ok((curve, point))
}

fn point_witness(a: &BigUint, b: &BigUint, p: &BigUint, n: &BigUint) -> Option<PointWitness> {
    for x_small in 0u32..512 {
        let x = BigUint::from(x_small);
        let rhs = (x.modpow(&BigUint::from(3u32), p) + a * &x + b) % p;
        let Some(y) = sqrt_p256(&rhs, p) else {
            continue;
        };
        if y.is_zero() {
            continue;
        }
        let witness = PointWitness {
            x: hex_big(&x),
            y: hex_big(&y),
        };
        let Ok((curve, point)) = as_curve(a, b, &witness) else {
            continue;
        };
        if curve.is_on_curve(&point)
            && matches!(point.scalar_mul_vartime(n, &curve.a_fe()), Point::Infinity)
        {
            return Some(witness);
        }
    }
    None
}

fn model_for_j(j: &BigUint, p: &BigUint, n: &BigUint) -> Result<CurveRecord, String> {
    let exceptional_zero = j.is_zero();
    let exceptional_1728 = j == &BigUint::from(1728u32);
    if exceptional_zero || exceptional_1728 {
        return Err("degree-11 walk reached unsupported exceptional j-invariant".into());
    }
    let denominator = (BigUint::from(1728u32) + p - j) % p;
    let k = (j * denominator.modpow(&(p - BigUint::from(2u32)), p)) % p;
    let a = (BigUint::from(3u32) * &k) % p;
    let b = (BigUint::from(2u32) * &k) % p;
    for candidate_b in [
        b.clone(),
        if b.is_zero() { BigUint::zero() } else { p - &b },
    ] {
        if let Some(order_witness) = point_witness(&a, &candidate_b, p, n) {
            let actual_j = curve_j(&a, &candidate_b, p)
                .ok_or_else(|| "canonical target model is singular".to_string())?;
            if actual_j != *j {
                return Err("canonical target model changed j-invariant".into());
            }
            return Ok(CurveRecord {
                a: hex_big(&a),
                b: hex_big(&candidate_b),
                j: hex_big(j),
                order_witness,
            });
        }
    }
    Err("neither quadratic twist supplied a P-256-order witness".into())
}

fn p256_start() -> Result<CurveRecord, String> {
    let curve = CurveParams::p256();
    let j = curve_j(&curve.a, &curve.b, &curve.p)
        .ok_or_else(|| "registered P-256 is singular".to_string())?;
    Ok(CurveRecord {
        a: hex_big(&curve.a),
        b: hex_big(&curve.b),
        j: hex_big(&j),
        order_witness: PointWitness {
            x: hex_big(&curve.gx),
            y: hex_big(&curve.gy),
        },
    })
}

fn edge_lines(edges: &[EdgeRecord]) -> Result<Vec<u8>, String> {
    let mut bytes = Vec::new();
    for edge in edges {
        let mut line = serde_json::to_vec(edge).map_err(|error| error.to_string())?;
        bytes.append(&mut line);
        bytes.push(b'\n');
    }
    Ok(bytes)
}

pub fn generate(steps: u64) -> Result<WalkCertificate, String> {
    if steps == 0 {
        return Err("walk must emit at least one edge".into());
    }
    let started = Instant::now();
    let p256 = CurveParams::p256();
    let trace = p256_trace(&p256);
    let phi = Phi11::load()?;
    let census_entries = degree_census();
    let selected_degree = selected_walk_degree(&census_entries)
        .ok_or_else(|| "no split prime at most 31".to_string())?;
    if selected_degree != WALK_DEGREE {
        return Err(format!(
            "frozen census expected degree {WALK_DEGREE}, selected {selected_degree}"
        ));
    }

    let start_curve = p256_start()?;
    let mut current_j = parse_hex("start j", &start_curve.j)?;
    let mut predecessor_j: Option<BigUint> = None;
    let mut edges = Vec::with_capacity(steps as usize);
    let mut first_seen = HashMap::new();
    first_seen.insert(current_j.clone(), 0u64);
    let mut first_cycle = None;

    for index in 0..steps {
        let roots = modular_roots(&phi, &current_j, &p256.p)?;
        let (target_j, direction) = match &predecessor_j {
            None => (roots[0].clone(), "lower-j-root".to_string()),
            Some(predecessor) if roots[0] == *predecessor => {
                (roots[1].clone(), "exclude-predecessor".to_string())
            }
            Some(predecessor) if roots[1] == *predecessor => {
                (roots[0].clone(), "exclude-predecessor".to_string())
            }
            Some(_) => {
                return Err(format!(
                    "edge {index}: predecessor is not a degree-11 neighbour"
                ))
            }
        };
        let target = model_for_j(&target_j, &p256.p, &p256.n)?;
        let edge = EdgeRecord {
            record: "edge".into(),
            index,
            predecessor_j: predecessor_j.as_ref().map(hex_big),
            source_j: hex_big(&current_j),
            target,
            direction,
            isogeny: IsogenyWitness {
                kind: "classical-Phi_11-plus-exact-order-twist-witness".into(),
                degree: WALK_DEGREE,
                phi11_sha256: PHI11_SHA256.into(),
                phi_value: "0".into(),
            },
        };
        let target_vertex = index + 1;
        if first_cycle.is_none() {
            if let Some(first_vertex) = first_seen.get(&target_j) {
                first_cycle = Some(CycleEvent {
                    at_vertex: target_vertex,
                    first_seen_vertex: *first_vertex,
                    period: target_vertex - first_vertex,
                    j: hex_big(&target_j),
                });
            }
        }
        first_seen.entry(target_j.clone()).or_insert(target_vertex);
        predecessor_j = Some(current_j);
        current_j = target_j;
        edges.push(edge);
    }

    let serialized_edges = edge_lines(&edges)?;
    let generation_ms = started.elapsed().as_millis();
    let elapsed_seconds = generation_ms.max(1) as f64 / 1000.0;
    let edge_bytes = serialized_edges.len() as u64;
    let summary = SummaryRecord {
        record: "summary".into(),
        emitted_steps: steps,
        verified_steps: 0,
        distinct_j: first_seen.len() as u64,
        first_cycle,
        generation_ms,
        edges_per_second: steps as f64 / elapsed_seconds,
        serialized_edge_bytes: edge_bytes,
        bytes_per_edge: edge_bytes as f64 / steps as f64,
        edge_stream_sha256: hex::encode(sha256(&serialized_edges)),
    };
    let mut certificate = WalkCertificate {
        header: WalkHeader {
            record: "header".into(),
            schema: CERTIFICATE_SCHEMA.into(),
            curve: P256_ICV1.into(),
            field_p: hex_big(&p256.p),
            group_order: hex_big(&p256.n),
            trace: hex_big(&trace),
            degree: WALK_DEGREE,
            requested_steps: steps,
            phi11_source: PHI11_SOURCE.into(),
            phi11_sha256: PHI11_SHA256.into(),
        },
        census: CensusRecord {
            record: "census".into(),
            entries: census_entries,
            selected_degree,
        },
        start: StartRecord {
            record: "start".into(),
            curve: start_curve,
        },
        edges,
        summary,
    };
    verify_inner(&mut certificate, true)?;
    Ok(certificate)
}

fn verify_curve_record(record: &CurveRecord, p: &BigUint, n: &BigUint) -> Result<BigUint, String> {
    let a = parse_hex("curve a", &record.a)?;
    let b = parse_hex("curve b", &record.b)?;
    let claimed_j = parse_hex("curve j", &record.j)?;
    if a >= *p || b >= *p || claimed_j >= *p {
        return Err("curve field element is not canonical".into());
    }
    let actual_j = curve_j(&a, &b, p).ok_or_else(|| "curve is singular".to_string())?;
    if actual_j != claimed_j {
        return Err("curve j-invariant does not match its model".into());
    }
    let (curve, point) = as_curve(&a, &b, &record.order_witness)?;
    if matches!(point, Point::Infinity) || !curve.is_on_curve(&point) {
        return Err("order witness is not a nonidentity curve point".into());
    }
    if !matches!(point.scalar_mul_vartime(n, &curve.a_fe()), Point::Infinity) {
        return Err("order witness is not killed by the P-256 group order".into());
    }
    Ok(actual_j)
}

fn hasse_interval_forces_exact_order(p: &BigUint, n: &BigUint) -> bool {
    let (floor, exact) = crate::cryptanalysis::p256_structural::isqrt_with_exactness(p);
    let ceil = if exact { floor } else { floor + BigUint::one() };
    p + BigUint::one() + BigUint::from(2u32) * ceil < BigUint::from(2u32) * n
}

fn verify_inner(certificate: &mut WalkCertificate, allow_unsealed: bool) -> Result<(), String> {
    let p256 = CurveParams::p256();
    let trace = p256_trace(&p256);
    let phi = Phi11::load()?;
    if certificate.header.record != "header"
        || certificate.header.schema != CERTIFICATE_SCHEMA
        || certificate.header.curve != P256_ICV1
        || certificate.header.field_p != hex_big(&p256.p)
        || certificate.header.group_order != hex_big(&p256.n)
        || certificate.header.trace != hex_big(&trace)
        || certificate.header.degree != WALK_DEGREE
        || certificate.header.phi11_source != PHI11_SOURCE
        || certificate.header.phi11_sha256 != PHI11_SHA256
    {
        return Err("certificate header does not name the frozen P-256 walk".into());
    }
    if !hasse_interval_forces_exact_order(&p256.p, &p256.n) {
        return Err("Hasse interval does not make the order witness decisive".into());
    }
    let expected_census = degree_census();
    if certificate.census.record != "census"
        || certificate.census.entries != expected_census
        || certificate.census.selected_degree != WALK_DEGREE
    {
        return Err("degree census does not match the frozen native census".into());
    }
    let expected_start = p256_start()?;
    if certificate.start.record != "start" || certificate.start.curve != expected_start {
        return Err("certificate does not start at the exact registered P-256 model".into());
    }
    verify_curve_record(&certificate.start.curve, &p256.p, &p256.n)?;
    let unsealed = certificate.summary.verified_steps == 0;
    if certificate.header.requested_steps == 0
        || certificate.summary.record != "summary"
        || certificate.edges.len() as u64 != certificate.header.requested_steps
        || certificate.summary.emitted_steps != certificate.header.requested_steps
        || (certificate.summary.verified_steps != certificate.header.requested_steps
            && !(allow_unsealed && unsealed))
    {
        return Err("certificate did not emit and verify the requested edge count".into());
    }

    let mut current_j = parse_hex("start j", &certificate.start.curve.j)?;
    let first_roots = modular_roots(&phi, &current_j, &p256.p)?;
    let mut predecessor_j: Option<BigUint> = None;
    let mut first_seen = HashMap::new();
    first_seen.insert(current_j.clone(), 0u64);
    let mut first_cycle = None;
    for (position, edge) in certificate.edges.iter().enumerate() {
        let index = position as u64;
        if edge.record != "edge" || edge.index != index || edge.source_j != hex_big(&current_j) {
            return Err(format!("edge {index}: discontinuous source or index"));
        }
        if edge.predecessor_j != predecessor_j.as_ref().map(hex_big) {
            return Err(format!("edge {index}: predecessor mismatch"));
        }
        if edge.isogeny.degree != WALK_DEGREE
            || edge.isogeny.phi11_sha256 != PHI11_SHA256
            || edge.isogeny.phi_value != "0"
            || edge.isogeny.kind != "classical-Phi_11-plus-exact-order-twist-witness"
        {
            return Err(format!("edge {index}: invalid isogeny witness metadata"));
        }
        let target_j = verify_curve_record(&edge.target, &p256.p, &p256.n)?;
        if !phi.evaluate(&current_j, &target_j, &p256.p).is_zero() {
            return Err(format!("edge {index}: Phi_11 relation failed"));
        }
        if let Some(predecessor) = &predecessor_j {
            if target_j == *predecessor || edge.direction != "exclude-predecessor" {
                return Err(format!("edge {index}: immediate backtracking"));
            }
        } else if target_j != first_roots[0] || edge.direction != "lower-j-root" {
            return Err("edge 0 does not use the frozen lower-j direction".into());
        }
        let target_vertex = index + 1;
        if first_cycle.is_none() {
            if let Some(first_vertex) = first_seen.get(&target_j) {
                first_cycle = Some(CycleEvent {
                    at_vertex: target_vertex,
                    first_seen_vertex: *first_vertex,
                    period: target_vertex - first_vertex,
                    j: hex_big(&target_j),
                });
            }
        }
        first_seen.entry(target_j.clone()).or_insert(target_vertex);
        predecessor_j = Some(current_j);
        current_j = target_j;
    }

    let serialized_edges = edge_lines(&certificate.edges)?;
    let expected_bytes_per_edge =
        serialized_edges.len() as f64 / certificate.header.requested_steps as f64;
    if certificate.summary.serialized_edge_bytes != serialized_edges.len() as u64
        || certificate.summary.bytes_per_edge != expected_bytes_per_edge
        || certificate.summary.edge_stream_sha256 != hex::encode(sha256(&serialized_edges))
        || certificate.summary.distinct_j != first_seen.len() as u64
        || certificate.summary.first_cycle != first_cycle
    {
        return Err("certificate summary does not match the edge stream".into());
    }
    certificate.summary.verified_steps = certificate.edges.len() as u64;
    Ok(())
}

/// Replay a sealed certificate. An unverified (`verified_steps = 0`) record is
/// accepted only inside [`generate`], immediately before it is sealed.
pub fn verify(certificate: &mut WalkCertificate) -> Result<(), String> {
    verify_inner(certificate, false)
}

pub fn to_json_lines(certificate: &WalkCertificate) -> Result<Vec<u8>, String> {
    let mut output = Vec::new();
    for line in [
        serde_json::to_vec(&certificate.header).map_err(|error| error.to_string())?,
        serde_json::to_vec(&certificate.census).map_err(|error| error.to_string())?,
        serde_json::to_vec(&certificate.start).map_err(|error| error.to_string())?,
    ] {
        output.extend(line);
        output.push(b'\n');
    }
    output.extend(edge_lines(&certificate.edges)?);
    output.extend(serde_json::to_vec(&certificate.summary).map_err(|error| error.to_string())?);
    output.push(b'\n');
    Ok(output)
}

pub fn from_json_lines(bytes: &[u8]) -> Result<WalkCertificate, String> {
    let lines: Vec<&[u8]> = bytes
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
        .collect();
    if lines.len() < 4 {
        return Err("certificate has too few JSON Lines records".into());
    }
    let header = serde_json::from_slice(lines[0]).map_err(|error| error.to_string())?;
    let census = serde_json::from_slice(lines[1]).map_err(|error| error.to_string())?;
    let start = serde_json::from_slice(lines[2]).map_err(|error| error.to_string())?;
    let summary =
        serde_json::from_slice(lines[lines.len() - 1]).map_err(|error| error.to_string())?;
    let edges = lines[3..lines.len() - 1]
        .iter()
        .map(|line| serde_json::from_slice(line).map_err(|error| error.to_string()))
        .collect::<Result<Vec<EdgeRecord>, String>>()?;
    Ok(WalkCertificate {
        header,
        census,
        start,
        edges,
        summary,
    })
}

fn poly_pow_small(polynomial: &FpPoly, exponent: usize) -> FpPoly {
    let mut result = FpPoly::one(polynomial.p.clone());
    for _ in 0..exponent {
        result = result.mul(polynomial);
    }
    result
}

fn division_g(
    index: usize,
    a: &BigUint,
    b: &BigUint,
    p: &BigUint,
    cache: &mut HashMap<usize, FpPoly>,
) -> FpPoly {
    if let Some(value) = cache.get(&index) {
        return value.clone();
    }
    let x = FpPoly::x(p.clone());
    let f = poly_pow_small(&x, 3)
        .add(&x.scalar_mul(a))
        .add(&FpPoly::constant(b.clone(), p.clone()));
    let value = match index {
        0 => FpPoly::zero(p.clone()),
        1 | 2 => FpPoly::one(p.clone()),
        3 => {
            let mut coefficients = vec![BigUint::zero(); 5];
            coefficients[0] = if a.is_zero() {
                BigUint::zero()
            } else {
                p - (a * a % p)
            };
            coefficients[1] = (BigUint::from(12u32) * b) % p;
            coefficients[2] = (BigUint::from(6u32) * a) % p;
            coefficients[4] = BigUint::from(3u32);
            FpPoly::from_coeffs(coefficients, p.clone())
        }
        4 => {
            let a2 = a * a % p;
            let a3 = &a2 * a % p;
            let b2 = b * b % p;
            let mut coefficients = vec![BigUint::zero(); 7];
            coefficients[0] = if (&a3 + BigUint::from(8u32) * &b2) % p == BigUint::zero() {
                BigUint::zero()
            } else {
                p - ((BigUint::from(2u32) * (&a3 + BigUint::from(8u32) * b2)) % p)
            };
            let ab = a * b % p;
            coefficients[1] = if ab.is_zero() {
                BigUint::zero()
            } else {
                p - ((BigUint::from(8u32) * ab) % p)
            };
            coefficients[2] = if a2.is_zero() {
                BigUint::zero()
            } else {
                p - ((BigUint::from(10u32) * a2) % p)
            };
            coefficients[3] = (BigUint::from(40u32) * b) % p;
            coefficients[4] = (BigUint::from(10u32) * a) % p;
            coefficients[6] = BigUint::from(2u32);
            FpPoly::from_coeffs(coefficients, p.clone())
        }
        odd if !odd.is_multiple_of(2) => {
            let m = (odd - 1) / 2;
            let gm2 = division_g(m + 2, a, b, p, cache);
            let gm = division_g(m, a, b, p, cache);
            let gm1 = division_g(m - 1, a, b, p, cache);
            let gp1 = division_g(m + 1, a, b, p, cache);
            let first = gm2.mul(&poly_pow_small(&gm, 3));
            let second = gm1.mul(&poly_pow_small(&gp1, 3));
            let sixteen_f2 = poly_pow_small(&f, 2).scalar_mul(&BigUint::from(16u32));
            if !m.is_multiple_of(2) {
                first.sub(&sixteen_f2.mul(&second))
            } else {
                sixteen_f2.mul(&first).sub(&second)
            }
        }
        even => {
            let m = even / 2;
            let gm = division_g(m, a, b, p, cache);
            let gp2 = division_g(m + 2, a, b, p, cache);
            let gm1 = division_g(m - 1, a, b, p, cache);
            let gm2 = division_g(m - 2, a, b, p, cache);
            let gp1 = division_g(m + 1, a, b, p, cache);
            gm.mul(
                &gp2.mul(&poly_pow_small(&gm1, 2))
                    .sub(&gm2.mul(&poly_pow_small(&gp1, 2))),
            )
        }
    };
    cache.insert(index, value.clone());
    value
}

fn scalar_x_rational(
    scalar: usize,
    a: &BigUint,
    b: &BigUint,
    p: &BigUint,
    cache: &mut HashMap<usize, FpPoly>,
) -> (FpPoly, FpPoly) {
    let x = FpPoly::x(p.clone());
    let f = poly_pow_small(&x, 3)
        .add(&x.scalar_mul(a))
        .add(&FpPoly::constant(b.clone(), p.clone()));
    let g = division_g(scalar, a, b, p, cache);
    let denominator = if scalar.is_multiple_of(2) {
        poly_pow_small(&g, 2)
            .mul(&f)
            .scalar_mul(&BigUint::from(4u32))
    } else {
        poly_pow_small(&g, 2)
    };
    let adjacent =
        division_g(scalar - 1, a, b, p, cache).mul(&division_g(scalar + 1, a, b, p, cache));
    let adjacent = if !scalar.is_multiple_of(2) {
        adjacent.mul(&f).scalar_mul(&BigUint::from(4u32))
    } else {
        adjacent
    };
    (x.mul(&denominator).sub(&adjacent), denominator)
}

fn compose_mod(polynomial: &FpPoly, argument: &FpPoly, modulus: &FpPoly) -> FpPoly {
    let mut result = FpPoly::zero(polynomial.p.clone());
    for coefficient in polynomial.coeffs.iter().rev() {
        result = result
            .mul(argument)
            .add(&FpPoly::constant(coefficient.clone(), polynomial.p.clone()))
            .rem(modulus);
    }
    result
}

fn power_sums_1_2_3(monic: &FpPoly) -> Result<[BigUint; 3], String> {
    let degree = monic
        .degree()
        .ok_or_else(|| "zero kernel polynomial".to_string())?;
    if degree < 3 || monic.lead() != BigUint::one() {
        return Err("kernel polynomial must be monic of degree at least three".into());
    }
    let p = &monic.p;
    let a1 = monic.coeff(degree - 1);
    let a2 = monic.coeff(degree - 2);
    let a3 = monic.coeff(degree - 3);
    let s1 = if a1.is_zero() {
        BigUint::zero()
    } else {
        p - &a1
    };
    let s2_sum = (&a1 * &s1 + BigUint::from(2u32) * &a2) % p;
    let s2 = if s2_sum.is_zero() {
        BigUint::zero()
    } else {
        p - s2_sum
    };
    let s3_sum = (&a1 * &s2 + &a2 * &s1 + BigUint::from(3u32) * a3) % p;
    let s3 = if s3_sum.is_zero() {
        BigUint::zero()
    } else {
        p - s3_sum
    };
    Ok([s1, s2, s3])
}

fn velu_codomain_from_kernel(
    a: &BigUint,
    b: &BigUint,
    kernel: &FpPoly,
) -> Result<(BigUint, BigUint), String> {
    let p = &kernel.p;
    let degree = kernel.degree().ok_or_else(|| "empty kernel".to_string())?;
    let [s1, s2, s3] = power_sums_1_2_3(kernel)?;
    let sum_v = (BigUint::from(6u32) * s2 + BigUint::from(2u32) * a * BigUint::from(degree)) % p;
    let sum_uvx = (BigUint::from(10u32) * s3
        + BigUint::from(6u32) * a * s1
        + BigUint::from(4u32) * b * BigUint::from(degree))
        % p;
    let five_v = BigUint::from(5u32) * sum_v % p;
    let seven_uvx = BigUint::from(7u32) * sum_uvx % p;
    let codomain_a = if five_v.is_zero() {
        a.clone()
    } else {
        (a + p - five_v) % p
    };
    let codomain_b = if seven_uvx.is_zero() {
        b.clone()
    } else {
        (b + p - seven_uvx) % p
    };
    Ok((codomain_a, codomain_b))
}

/// Slow, algebraically independent check of the first modular-polynomial
/// step. It obtains the two kernels from `psi_11`, the Frobenius eigenvalue
/// equations and x-only scalar multiplication, then applies Vélu's sums.
pub fn division_polynomial_fixtures() -> Result<Vec<DivisionPolynomialFixture>, String> {
    let curve = CurveParams::p256();
    let mut cache = HashMap::new();
    let division = division_g(11, &curve.a, &curve.b, &curve.p, &mut cache).monic();
    if division.degree() != Some(60) {
        return Err("psi_11 does not have degree 60".into());
    }
    let x = FpPoly::x(curve.p.clone());
    let x_frobenius = x.pow_mod(&curve.p, &division);
    let mut fixtures = Vec::new();
    for scalar in [2usize, 3usize] {
        let (numerator, denominator) =
            scalar_x_rational(scalar, &curve.a, &curve.b, &curve.p, &mut cache);
        let equation = x_frobenius.mul(&denominator).sub(&numerator).rem(&division);
        let kernel = division.gcd(&equation).monic();
        if kernel.degree() != Some(5) || !division.rem(&kernel).is_zero() {
            return Err(format!(
                "scalar {scalar}: Frobenius kernel is not degree five"
            ));
        }
        let (x_num, x_den) = scalar_x_rational(2, &curve.a, &curve.b, &curve.p, &mut cache);
        let (_, inverse_den, _) = x_den.ext_gcd(&kernel);
        let doubled_x = x_num.mul(&inverse_den).rem(&kernel);
        if !compose_mod(&kernel, &doubled_x, &kernel).is_zero() {
            return Err(format!(
                "scalar {scalar}: kernel roots are not closed under doubling"
            ));
        }
        let (codomain_a, codomain_b) = velu_codomain_from_kernel(&curve.a, &curve.b, &kernel)?;
        let codomain_j = curve_j(&codomain_a, &codomain_b, &curve.p)
            .ok_or_else(|| "Vélu codomain is singular".to_string())?;
        fixtures.push(DivisionPolynomialFixture {
            scalar: scalar as u64,
            kernel_coefficients: kernel.coeffs.iter().map(hex_big).collect(),
            codomain_a: hex_big(&codomain_a),
            codomain_b: hex_big(&codomain_b),
            codomain_j: hex_big(&codomain_j),
        });
    }
    fixtures.sort_by(|left, right| left.codomain_j.cmp(&right.codomain_j));
    Ok(fixtures)
}

pub fn verify_independent_start_fixture() -> Result<Vec<DivisionPolynomialFixture>, String> {
    let phi = Phi11::load()?;
    let curve = CurveParams::p256();
    let start_j = curve_j(&curve.a, &curve.b, &curve.p)
        .ok_or_else(|| "registered P-256 is singular".to_string())?;
    let modular: BTreeSet<String> = modular_roots(&phi, &start_j, &curve.p)?
        .into_iter()
        .map(|value| hex_big(&value))
        .collect();
    let fixtures = division_polynomial_fixtures()?;
    let division: BTreeSet<String> = fixtures
        .iter()
        .map(|fixture| fixture.codomain_j.clone())
        .collect();
    if modular != division {
        return Err("division-polynomial/Vélu neighbours disagree with Phi_11".into());
    }
    Ok(fixtures)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn census_selects_eleven_after_two_ramified_dead_ends() {
        let census = degree_census();
        assert_eq!(selected_walk_degree(&census), Some(11));
        assert_eq!(
            census.iter().find(|entry| entry.degree == 3).unwrap().class,
            "ramified-single-edge"
        );
        assert_eq!(
            census.iter().find(|entry| entry.degree == 5).unwrap().class,
            "ramified-single-edge"
        );
        assert_eq!(
            census.iter().find(|entry| entry.degree == 7).unwrap().class,
            "inert-no-edge"
        );
        assert_eq!(
            census
                .iter()
                .find(|entry| entry.degree == 11)
                .unwrap()
                .fp_edges,
            2
        );
    }

    #[test]
    fn short_walk_roundtrips_and_rejects_a_wrong_twist() {
        let mut certificate = generate(4).unwrap();
        verify(&mut certificate).unwrap();
        assert_eq!(certificate.summary.verified_steps, 4);

        let mut header = serde_json::to_value(&certificate.header).unwrap();
        header["unexpected"] = serde_json::json!(true);
        assert!(serde_json::from_value::<WalkHeader>(header).is_err());

        let bytes = to_json_lines(&certificate).unwrap();
        let mut decoded = from_json_lines(&bytes).unwrap();
        verify(&mut decoded).unwrap();

        let p = CurveParams::p256().p;
        let first = &mut decoded.edges[0].target;
        let b = parse_hex("b", &first.b).unwrap();
        first.b = hex_big(&(p - b));
        assert!(verify(&mut decoded).unwrap_err().contains("order witness"));
    }

    #[test]
    fn division_polynomial_fixture_agrees_with_classical_phi_11() {
        let fixtures = verify_independent_start_fixture().unwrap();
        assert_eq!(fixtures.len(), 2);
        assert_eq!(
            fixtures
                .iter()
                .map(|fixture| fixture.codomain_j.as_str())
                .collect::<Vec<_>>(),
            vec![
                "1e83ce42e311c91e26e534220f6c8af76ef81a2c22246611462ef9b5b044ff88",
                "f17b32fb5796946d1d2b09fd642c7661428bf217216ae83fb6354a5c87fefbc0",
            ]
        );
        assert_eq!(
            fixtures[0].kernel_coefficients[0],
            "f136e7e36906a7ed31e3dc2ef549d5e0001c8c6237f7205b0d86d62ddfacb549"
        );
        assert_eq!(
            fixtures[1].kernel_coefficients[0],
            "c066675b2ac76605dbb2aa7c9a20755a5906ef4eb6babd4be3d19dcf720b014d"
        );
        assert!(fixtures
            .iter()
            .all(|fixture| fixture.kernel_coefficients.len() == 6));
    }
}
