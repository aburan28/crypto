//! Exact, certificate-oriented CM relation census for NIST P-192.
//!
//! The public surface exists for the dedicated binary and independent replay;
//! it is not a general-purpose map-cost API.  The exact lane is complete:
//! pinned curve/order checks, a recursive Pocklington proof for the large
//! discriminant cofactor, oriented prime-ideal fixtures, exhaustive bounded
//! vector enumeration, form replay, and independent two-dimensional ideal-HNF
//! replay.

use crate::cryptanalysis::isogeny_walk::field::is_probable_prime;
use crate::ecc::curve::CurveParams;
use crate::ecc::point::Point;
use crate::hash::sha256::sha256;
use crate::isogeny::class_group::{class_above_prime, BinaryQuadraticForm};
use blake3::Hasher;
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeSet, HashSet};
use std::sync::OnceLock;

pub const SCHEMA_ORDER: &str = "p192.cm_order/v1";
pub const SCHEMA_POCKLINGTON: &str = "p192.pocklington/v1";
pub const SCHEMA_GENERATORS: &str = "p192.cm_generators/v1";
pub const SCHEMA_WEIGHTS: &str = "p192.cm_weights/v1";
pub const SCHEMA_SEARCH: &str = "p192.cm_relation_search/v2";
pub const SCHEMA_VERIFY: &str = "p192.cm_relation_verification/v1";

pub const P_DEC: &str = "6277101735386680763835789423207666416083908700390324961279";
pub const N_DEC: &str = "6277101735386680763835789423176059013767194773182842284081";
pub const T_DEC: &str = "31607402316713927207482677199";
pub const D_DEC: &str = "-24109379060336110122544161233113975664949272517896865359515";
pub const C_DEC: &str = "14140398275856956083603613626459809774163796198179979683";
pub const Q_DEC: &str = "3654106800343397140285403541412623476381";
pub const NONSCALAR_BOUND_DEC: &str = "6027344765084027530636040308278493916237318129474216339879";
pub const SYMBOLIC_LAMBDA_DEC: &str = "6277101735386680763835789423144451611450480845975359606884";
const A_HEX: &str = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFC";
const B_HEX: &str = "64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1";
const GX_HEX: &str = "188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012";
const GY_HEX: &str = "07192B95FFC8DA78631011ED6B24CDD573F977A11E794811";

pub const MAX_DEGREE: u64 = 1u64 << 48;
pub const EXPECTED_ORIENTED: u64 = 9_948_061;
pub const EXPECTED_CANONICAL: u64 = 4_974_348;
pub const EXPECTED_CANONICAL_NONZERO: u64 = 4_974_347;
pub const HARD_STATE_CAP: u64 = 5_000_000;
pub const EXPECTED_SCALAR_INCLUDING_ZERO: u64 = 104;
pub const EXPECTED_SCALAR_NONZERO: u64 = 103;

const EXPERIMENT_INCOMPLETE_STATUS: &str = "incomplete-map-calibration-and-payoff-open";
const VERIFICATION_NO_HIT_VERDICT: &str =
    "incomplete-exact-algebra-pass-no-nonscalar-relation-in-frozen-box";
const VERIFICATION_CANDIDATE_VERDICT: &str = "incomplete-candidate-only-map-and-payoff-open";

// Boundary members use the same fixed 13*i16 little-endian vector encoding as
// the census.  Concatenating zero members is exactly the empty byte string.
const EMPTY_BOUNDARY_ENCODING: &[u8] = b"";
const EMPTY_BOUNDARY_ENCODING_RULE: &str =
    "concatenated 13xi16-le canonical vectors in forward enumeration order; empty set is zero bytes";

const RAMIFIED: [(u64, u64); 3] = [(5, 2), (11, 3), (31, 7)];
const SPLIT: [(u64, u64, u64); 10] = [
    (13, 2, 5),
    (23, 21, 22),
    (37, 12, 35),
    (43, 8, 26),
    (73, 60, 67),
    (89, 6, 83),
    (101, 17, 70),
    (103, 5, 36),
    (107, 56, 68),
    (113, 26, 42),
];

const C_MINUS_ONE: [(&str, u32); 8] = [
    ("2", 1),
    ("3", 2),
    ("7", 1),
    ("17", 1),
    ("139", 1),
    ("11471", 1),
    ("1133039", 1),
    (Q_DEC, 1),
];

const Q_MINUS_ONE: [(&str, u32); 5] = [
    ("2", 2),
    ("5", 1),
    ("28929853", 1),
    ("386058915559", 1),
    ("16358799425486714897", 1),
];

const P28929853_MINUS_ONE: [(&str, u32); 5] =
    [("2", 2), ("3", 3), ("7", 1), ("17", 1), ("2251", 1)];
const P386058915559_MINUS_ONE: [(&str, u32); 4] = [("2", 1), ("3", 3), ("17", 1), ("420543481", 1)];
const P16358799425486714897_MINUS_ONE: [(&str, u32); 4] =
    [("2", 4), ("7", 1), ("449", 1), ("325302247563767", 1)];
const P1133039_MINUS_ONE: [(&str, u32); 3] = [("2", 1), ("397", 1), ("1427", 1)];
const P420543481_MINUS_ONE: [(&str, u32); 6] = [
    ("2", 3),
    ("3", 1),
    ("5", 1),
    ("7", 2),
    ("37", 1),
    ("1933", 1),
];
const P325302247563767_MINUS_ONE: [(&str, u32); 5] =
    [("2", 1), ("11", 1), ("839", 1), ("39371", 1), ("447637", 1)];
const P447637_MINUS_ONE: [(&str, u32); 4] = [("2", 2), ("3", 1), ("7", 1), ("73", 2)];

fn bu(s: &str) -> BigUint {
    BigUint::parse_bytes(s.as_bytes(), 10).expect("frozen decimal integer")
}

fn bu_hex(s: &str) -> BigUint {
    BigUint::parse_bytes(s.as_bytes(), 16).expect("frozen hexadecimal integer")
}

fn bi(s: &str) -> BigInt {
    BigInt::parse_bytes(s.as_bytes(), 10).expect("frozen signed decimal integer")
}

fn mod_u64(n: &BigInt, modulus: u64) -> u64 {
    n.mod_floor(&BigInt::from(modulus)).to_u64().unwrap()
}

fn form_strings(form: &BinaryQuadraticForm) -> [String; 3] {
    [form.a.to_string(), form.b.to_string(), form.c.to_string()]
}

fn strings_form(values: &[String; 3], disc: &BigInt) -> Result<BinaryQuadraticForm, String> {
    BinaryQuadraticForm::new(bi(&values[0]), bi(&values[1]), bi(&values[2]), disc)
        .map(|form| form.reduce())
        .ok_or_else(|| {
            format!(
                "invalid quadratic form [{}, {}, {}]",
                values[0], values[1], values[2]
            )
        })
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PocklingtonWitness {
    pub q: String,
    pub a: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FactorRecord {
    pub prime: String,
    pub exponent: u32,
    pub proof: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PocklingtonNode {
    pub n: String,
    pub factorization_of_n_minus_one: Vec<FactorRecord>,
    pub witnesses: Vec<PocklingtonWitness>,
    pub recursive_children: Vec<PocklingtonNode>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PocklingtonCertificate {
    pub schema: String,
    pub theorem: String,
    pub root: PocklingtonNode,
}

fn factor_product(factors: &[(&str, u32)]) -> BigUint {
    factors.iter().fold(BigUint::one(), |acc, (q, exponent)| {
        acc * bu(q).pow(*exponent)
    })
}

fn gcd_biguint(mut a: BigUint, mut b: BigUint) -> BigUint {
    while !b.is_zero() {
        let r = &a % &b;
        a = b;
        b = r;
    }
    a
}

fn trial_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n % 2 == 0 {
        return n == 2;
    }
    let mut divisor = 3u64;
    while divisor <= n / divisor {
        if n % divisor == 0 {
            return false;
        }
        divisor += 2;
    }
    true
}

fn common_pocklington_witness(n: &str) -> Result<&'static str, String> {
    match n {
        C_DEC | Q_DEC | "447637" => Ok("2"),
        "28929853" => Ok("15"),
        "386058915559" => Ok("6"),
        "16358799425486714897" => Ok("3"),
        "1133039" => Ok("13"),
        "420543481" => Ok("11"),
        "325302247563767" => Ok("5"),
        _ => Err(format!("no frozen common Pocklington witness for {n}")),
    }
}

fn make_node(
    n: &str,
    factors: &[(&str, u32)],
    children: Vec<PocklingtonNode>,
) -> Result<PocklingtonNode, String> {
    let value = bu(n);
    if factor_product(factors) != &value - 1u8 {
        return Err(format!(
            "hardcoded factorization does not multiply to {n}-1"
        ));
    }
    let child_names: BTreeSet<&str> = children.iter().map(|entry| entry.n.as_str()).collect();
    let common_witness = common_pocklington_witness(n)?;
    let mut records = Vec::with_capacity(factors.len());
    let mut witnesses = Vec::with_capacity(factors.len());
    for (q_text, exponent) in factors {
        let q = bu(q_text);
        let proof = if child_names.contains(*q_text) {
            "recursive-pocklington"
        } else {
            let q64 = q
                .to_u64()
                .ok_or_else(|| format!("leaf prime {q_text} does not fit u64"))?;
            if q64 >= 65_536 || !trial_prime_u64(q64) {
                return Err(format!("hardcoded leaf {q_text} is not prime"));
            }
            "exact-trial-division-below-65536"
        };
        records.push(FactorRecord {
            prime: (*q_text).to_owned(),
            exponent: *exponent,
            proof: proof.to_owned(),
        });
        witnesses.push(PocklingtonWitness {
            q: (*q_text).to_owned(),
            a: common_witness.to_owned(),
        });
    }
    Ok(PocklingtonNode {
        n: n.to_owned(),
        factorization_of_n_minus_one: records,
        witnesses,
        recursive_children: children,
    })
}

pub fn build_pocklington_certificate() -> Result<PocklingtonCertificate, String> {
    let p447637 = make_node("447637", &P447637_MINUS_ONE, vec![])?;
    let p325302247563767 = make_node(
        "325302247563767",
        &P325302247563767_MINUS_ONE,
        vec![p447637],
    )?;
    let p420543481 = make_node("420543481", &P420543481_MINUS_ONE, vec![])?;
    let p1133039 = make_node("1133039", &P1133039_MINUS_ONE, vec![])?;
    let p28929853 = make_node("28929853", &P28929853_MINUS_ONE, vec![])?;
    let p386058915559 = make_node("386058915559", &P386058915559_MINUS_ONE, vec![p420543481])?;
    let p16358799425486714897 = make_node(
        "16358799425486714897",
        &P16358799425486714897_MINUS_ONE,
        vec![p325302247563767],
    )?;
    let q = make_node(
        Q_DEC,
        &Q_MINUS_ONE,
        vec![p28929853, p386058915559, p16358799425486714897],
    )?;
    let c = make_node(C_DEC, &C_MINUS_ONE, vec![p1133039, q])?;
    let certificate = PocklingtonCertificate {
        schema: SCHEMA_POCKLINGTON.to_owned(),
        theorem: "recursive Pocklington with complete n-1 factorizations at nine nodes; terminal primes below 65536 use exact trial division".to_owned(),
        root: c,
    };
    verify_pocklington_certificate(&certificate)?;
    Ok(certificate)
}

fn verify_node(node: &PocklingtonNode) -> Result<(), String> {
    let n = bu(&node.n);
    if n < BigUint::from(3u8) || n.is_even() {
        return Err(format!("invalid odd Pocklington candidate {n}"));
    }
    let child_names: BTreeSet<&str> = node
        .recursive_children
        .iter()
        .map(|child| child.n.as_str())
        .collect();
    for child in &node.recursive_children {
        verify_node(child)?;
    }
    let mut product = BigUint::one();
    let mut seen = BTreeSet::new();
    for factor in &node.factorization_of_n_minus_one {
        let q = bu(&factor.prime);
        if !seen.insert(factor.prime.as_str()) {
            return Err(format!("duplicate prime factor {}", factor.prime));
        }
        product *= q.pow(factor.exponent);
        if child_names.contains(factor.prime.as_str()) {
            if factor.proof != "recursive-pocklington" {
                return Err(format!("recursive factor {} mislabeled", factor.prime));
            }
        } else {
            let q64 = q
                .to_u64()
                .ok_or_else(|| format!("unproved non-u64 leaf {}", factor.prime))?;
            if factor.proof != "exact-trial-division-below-65536"
                || q64 >= 65_536
                || !trial_prime_u64(q64)
            {
                return Err(format!(
                    "invalid exact trial-division leaf proof for {}",
                    factor.prime
                ));
            }
        }
    }
    if product != &n - 1u8 {
        return Err(format!("factorization product mismatch for {n}"));
    }
    if node.witnesses.len() != node.factorization_of_n_minus_one.len() {
        return Err(format!("witness count mismatch for {n}"));
    }
    let nm1 = &n - 1u8;
    for factor in &node.factorization_of_n_minus_one {
        let witness = node
            .witnesses
            .iter()
            .find(|witness| witness.q == factor.prime)
            .ok_or_else(|| format!("missing witness for q={}", factor.prime))?;
        let q = bu(&factor.prime);
        let a = bu(&witness.a);
        if a < BigUint::from(2u8) || a >= n {
            return Err(format!("witness out of range for q={q}"));
        }
        if a.modpow(&nm1, &n) != BigUint::one() {
            return Err(format!("Fermat equality failed for n={n}, q={q}"));
        }
        let x = a.modpow(&(&nm1 / &q), &n);
        let delta = if x.is_zero() { &n - 1u8 } else { x - 1u8 };
        if !gcd_biguint(delta, n.clone()).is_one() {
            return Err(format!("Pocklington gcd failed for n={n}, q={q}"));
        }
    }
    Ok(())
}

pub fn verify_pocklington_certificate(certificate: &PocklingtonCertificate) -> Result<(), String> {
    if certificate.schema != SCHEMA_POCKLINGTON {
        return Err("wrong Pocklington schema".to_owned());
    }
    if certificate.root.n != C_DEC {
        return Err("Pocklington root is not the frozen C".to_owned());
    }
    verify_node(&certificate.root)
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OrderReport {
    pub schema: String,
    pub curve: String,
    pub provenance: String,
    pub p: String,
    pub n: String,
    pub trace_t: String,
    pub discriminant_d: String,
    pub c: String,
    pub pinned_standard_tuple_verified: bool,
    pub nonsingular: bool,
    pub p_probable_prime_20_bases: bool,
    pub n_probable_prime_20_bases: bool,
    pub c_proven_prime_pocklington: bool,
    pub generator_nonidentity: bool,
    pub generator_on_curve: bool,
    pub n_times_generator_is_identity: bool,
    pub hasse_n_in_interval: bool,
    pub hasse_unique_multiple: bool,
    pub d_factorization: Vec<FactorRecord>,
    pub d_squarefree: bool,
    pub d_fundamental: bool,
    pub frobenius_conductor: u32,
    pub endomorphism_order: String,
    pub nonscalar_norm_lower_bound: String,
    pub theorem_bound_verified: bool,
    pub status: String,
}

pub fn build_order_report(certificate: &PocklingtonCertificate) -> Result<OrderReport, String> {
    verify_pocklington_certificate(certificate)?;
    let curve = CurveParams::p192();
    let tuple_verified = curve.name == "P-192"
        && curve.p == bu(P_DEC)
        && curve.a == bu_hex(A_HEX)
        && curve.b == bu_hex(B_HEX)
        && curve.gx == bu_hex(GX_HEX)
        && curve.gy == bu_hex(GY_HEX)
        && curve.n == bu(N_DEC)
        && curve.h == 1;
    if !tuple_verified {
        return Err("repository P-192 tuple differs from frozen standard tuple".to_owned());
    }
    let p = bu(P_DEC);
    let n = bu(N_DEC);
    let t = bu(T_DEC);
    if &p + 1u8 - &n != t {
        return Err("trace does not equal p+1-n".to_owned());
    }
    let d = bi(D_DEC);
    let nonsingular = ((BigUint::from(4u8) * curve.a.modpow(&BigUint::from(3u8), &p))
        + (BigUint::from(27u8) * curve.b.modpow(&BigUint::from(2u8), &p)))
        % &p
        != BigUint::zero();
    let computed_d =
        BigInt::from_biguint(Sign::Plus, &t * &t) - BigInt::from_biguint(Sign::Plus, &p << 2usize);
    if computed_d != d {
        return Err("frozen discriminant is not t^2-4p".to_owned());
    }
    let generator = curve.generator();
    let generator_nonidentity = generator != Point::Infinity;
    let generator_on_curve = curve.is_on_curve(&generator);
    let n_times_generator_is_identity = generator.scalar_mul(&n, &curve.a_fe()) == Point::Infinity;
    let hasse_n_in_interval = &t * &t <= (&p << 2usize);
    // n is the only positive multiple of n in Hasse's interval: the upper
    // endpoint is below 2n.  Avoid square roots by squaring the positive gap.
    let two_n_gap = (&n << 1usize) - (&p + 1u8);
    let hasse_unique_multiple =
        two_n_gap > BigUint::zero() && &two_n_gap * &two_n_gap > (&p << 2usize);
    let abs_d = d.abs().to_biguint().unwrap();
    let expected_abs_d = BigUint::from(5u8) * 11u8 * 31u8 * bu(C_DEC);
    let d_squarefree = abs_d == expected_abs_d;
    let d_fundamental = d_squarefree && d.mod_floor(&BigInt::from(4u8)) == BigInt::one();
    let lower = bu(NONSCALAR_BOUND_DEC);
    let theorem_bound_verified = &lower * 4u8 >= abs_d && (&lower - 1u8) * 4u8 < abs_d;
    let p_probable = is_probable_prime(&p);
    let n_probable = is_probable_prime(&n);
    let all = p_probable
        && n_probable
        && tuple_verified
        && nonsingular
        && generator_nonidentity
        && generator_on_curve
        && n_times_generator_is_identity
        && hasse_n_in_interval
        && hasse_unique_multiple
        && d_squarefree
        && d_fundamental
        && theorem_bound_verified;
    Ok(OrderReport {
        schema: SCHEMA_ORDER.to_owned(),
        curve: "NIST P-192 / secp192r1".to_owned(),
        provenance: "SEC 2 v2 / FIPS pinned tuple; p and n receive a fixed 20-base probable-prime check, not a new proof".to_owned(),
        p: P_DEC.to_owned(),
        n: N_DEC.to_owned(),
        trace_t: T_DEC.to_owned(),
        discriminant_d: D_DEC.to_owned(),
        c: C_DEC.to_owned(),
        pinned_standard_tuple_verified: tuple_verified,
        nonsingular,
        p_probable_prime_20_bases: p_probable,
        n_probable_prime_20_bases: n_probable,
        c_proven_prime_pocklington: true,
        generator_nonidentity,
        generator_on_curve,
        n_times_generator_is_identity,
        hasse_n_in_interval,
        hasse_unique_multiple,
        d_factorization: vec![
            FactorRecord { prime: "5".into(), exponent: 1, proof: "exact-trial-division-below-65536".into() },
            FactorRecord { prime: "11".into(), exponent: 1, proof: "exact-trial-division-below-65536".into() },
            FactorRecord { prime: "31".into(), exponent: 1, proof: "exact-trial-division-below-65536".into() },
            FactorRecord { prime: C_DEC.into(), exponent: 1, proof: "recursive-pocklington".into() },
        ],
        d_squarefree,
        d_fundamental,
        frobenius_conductor: 1,
        endomorphism_order: "End(E)=Z[pi]=O_D (D fundamental)".to_owned(),
        nonscalar_norm_lower_bound: NONSCALAR_BOUND_DEC.to_owned(),
        theorem_bound_verified,
        status: if all { "pass" } else { "invalid" }.to_owned(),
    })
}

pub fn verify_order_report(
    report: &OrderReport,
    certificate: &PocklingtonCertificate,
) -> Result<(), String> {
    let expected = build_order_report(certificate)?;
    let left = serde_json::to_value(report).map_err(|error| error.to_string())?;
    let right = serde_json::to_value(expected).map_err(|error| error.to_string())?;
    if left != right {
        return Err("order report differs from native recomputation".to_owned());
    }
    if report.status != "pass" {
        return Err("order report did not pass".to_owned());
    }
    Ok(())
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GeneratorRecord {
    pub index: usize,
    pub ell: u64,
    pub kind: String,
    pub roots: Vec<u64>,
    pub positive_root: u64,
    pub negative_root: u64,
    pub positive_form: [String; 3],
    pub negative_form: [String; 3],
    pub class_above_prime_crosscheck: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ControlResult {
    pub name: String,
    pub passed: bool,
    pub detail: String,
}

fn expected_theorem_control() -> ControlResult {
    ControlResult {
        name: "p192-nonscalar-norm-bound".to_owned(),
        passed: bu(NONSCALAR_BOUND_DEC) > BigUint::from(MAX_DEGREE),
        detail: format!("every non-scalar alpha has norm >= {NONSCALAR_BOUND_DEC} > 2^48"),
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GeneratorReport {
    pub schema: String,
    pub curve: String,
    pub discriminant_d: String,
    pub trace_t: String,
    pub max_ell: u64,
    pub convention: String,
    pub generators: Vec<GeneratorRecord>,
    pub controls: Vec<ControlResult>,
    pub status: String,
}

fn characteristic_roots(ell: u64, t: &BigInt, p: &BigInt) -> Vec<u64> {
    let tm = mod_u64(t, ell) as u128;
    let pm = mod_u64(p, ell) as u128;
    let modulus = ell as u128;
    (0..ell)
        .filter(|root| {
            let r = *root as u128;
            (r * r + modulus - (tm * r) % modulus + pm) % modulus == 0
        })
        .collect()
}

fn form_for_root(
    disc: &BigInt,
    t: &BigInt,
    p: &BigInt,
    ell: u64,
    root: u64,
) -> Result<BinaryQuadraticForm, String> {
    let a = BigInt::from(ell);
    let r = BigInt::from(root);
    if (&r * &r - t * &r + p).mod_floor(&a) != BigInt::zero() {
        return Err(format!(
            "root {root} is not a characteristic root mod {ell}"
        ));
    }
    // For I=(ell,pi-r) and sqrt(D)=2pi-t, the associated form has
    // (-b+sqrt(D))/2 congruent to pi-r, hence b=2r-t (mod 2ell).
    // Reducing b to (-ell,ell] keeps the unreduced form small without
    // changing its oriented ideal class.
    let two_a = &a * 2;
    let mut b = (&r * BigInt::from(2u8) - t).mod_floor(&two_a);
    if b > a {
        b -= &two_a;
    }
    let numerator = &b * &b - disc;
    let denominator = &a * 4;
    if !numerator.mod_floor(&denominator).is_zero() {
        return Err(format!(
            "root/form parity mismatch at ell={ell}, root={root}"
        ));
    }
    let c = numerator / denominator;
    BinaryQuadraticForm::new(a, b, c, disc)
        .map(|form| form.reduce())
        .ok_or_else(|| format!("prime ideal form at ell={ell}, root={root} is invalid"))
}

fn generator_records() -> Result<Vec<GeneratorRecord>, String> {
    let d = bi(D_DEC);
    let t = bi(T_DEC);
    let p = bi(P_DEC);
    let mut records = Vec::with_capacity(13);
    for (index, (ell, root)) in RAMIFIED.iter().copied().enumerate() {
        let roots = characteristic_roots(ell, &t, &p);
        if roots != vec![root] {
            return Err(format!(
                "ramified root fixture mismatch at ell={ell}: {roots:?}"
            ));
        }
        let form = form_for_root(&d, &t, &p, ell, root)?;
        let generic = class_above_prime(&d, ell)
            .ok_or_else(|| format!("class_above_prime rejected ramified ell={ell}"))?;
        let cross = form == generic || form.inverse() == generic;
        records.push(GeneratorRecord {
            index,
            ell,
            kind: "ramified".to_owned(),
            roots,
            positive_root: root,
            negative_root: root,
            positive_form: form_strings(&form),
            negative_form: form_strings(&form.inverse()),
            class_above_prime_crosscheck: cross,
        });
    }
    for (offset, (ell, positive, negative)) in SPLIT.iter().copied().enumerate() {
        let roots = characteristic_roots(ell, &t, &p);
        if roots != vec![positive, negative] {
            return Err(format!(
                "split root fixture mismatch at ell={ell}: {roots:?}"
            ));
        }
        let pos = form_for_root(&d, &t, &p, ell, positive)?;
        let neg = form_for_root(&d, &t, &p, ell, negative)?;
        if pos.inverse() != neg {
            return Err(format!("oriented forms are not inverses at ell={ell}"));
        }
        let generic = class_above_prime(&d, ell)
            .ok_or_else(|| format!("class_above_prime rejected split ell={ell}"))?;
        let cross = pos == generic || neg == generic;
        records.push(GeneratorRecord {
            index: RAMIFIED.len() + offset,
            ell,
            kind: "split".to_owned(),
            roots,
            positive_root: positive,
            negative_root: negative,
            positive_form: form_strings(&pos),
            negative_form: form_strings(&neg),
            class_above_prime_crosscheck: cross,
        });
    }
    Ok(records)
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct IdealHnf {
    /// Columns `(a,0)` and `(b,d)` in the basis `{1,pi}`.
    pub a: String,
    pub b: String,
    pub d: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Hnf {
    a: BigInt,
    b: BigInt,
    d: BigInt,
}

impl Hnf {
    fn identity() -> Self {
        Self {
            a: BigInt::one(),
            b: BigInt::zero(),
            d: BigInt::one(),
        }
    }

    fn prime(ell: u64, root: u64) -> Self {
        let a = BigInt::from(ell);
        Self {
            a: a.clone(),
            b: (-BigInt::from(root)).mod_floor(&a),
            d: BigInt::one(),
        }
    }

    fn vectors(&self) -> [(BigInt, BigInt); 2] {
        [
            (self.a.clone(), BigInt::zero()),
            (self.b.clone(), self.d.clone()),
        ]
    }

    fn record(&self) -> IdealHnf {
        IdealHnf {
            a: self.a.to_string(),
            b: self.b.to_string(),
            d: self.d.to_string(),
        }
    }
}

fn extended_gcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    if b.is_zero() {
        let sign = if a.is_negative() {
            -BigInt::one()
        } else {
            BigInt::one()
        };
        return (a.abs(), sign, BigInt::zero());
    }
    let (g, x1, y1) = extended_gcd(b, &a.mod_floor(b));
    (g, y1.clone(), x1 - a.div_floor(b) * y1)
}

fn hnf_from_vectors(vectors: &[(BigInt, BigInt)]) -> Result<Hnf, String> {
    if vectors.len() < 2 {
        return Err("a rank-two HNF needs at least two generators".to_owned());
    }
    let mut index = BigInt::zero();
    for i in 0..vectors.len() {
        for j in (i + 1)..vectors.len() {
            let det = &vectors[i].0 * &vectors[j].1 - &vectors[j].0 * &vectors[i].1;
            index = index.gcd(&det.abs());
        }
    }
    if index.is_zero() {
        return Err("ideal generators have rank below two".to_owned());
    }
    let mut gy = BigInt::zero();
    let mut x_combo = BigInt::zero();
    for (x, y) in vectors {
        let (next, s, coefficient) = extended_gcd(&gy, y);
        x_combo = s * x_combo + coefficient * x;
        gy = next;
    }
    if gy.is_zero() || !index.mod_floor(&gy).is_zero() {
        return Err("inconsistent HNF y-gcd/index".to_owned());
    }
    let a = &index / &gy;
    let b = x_combo.mod_floor(&a);
    Ok(Hnf { a, b, d: gy })
}

fn multiply_hnf(left: &Hnf, right: &Hnf, p: &BigInt, t: &BigInt) -> Result<Hnf, String> {
    let mut products = Vec::with_capacity(4);
    for (x1, y1) in left.vectors() {
        for (x2, y2) in right.vectors() {
            let x = &x1 * &x2 - p * &y1 * &y2;
            let y = &x1 * &y2 + &x2 * &y1 + t * &y1 * &y2;
            products.push((x, y));
        }
    }
    hnf_from_vectors(&products)
}

fn principal_hnf(u: &BigInt, v: &BigInt, p: &BigInt, t: &BigInt) -> Result<Hnf, String> {
    hnf_from_vectors(&[(u.clone(), v.clone()), (-v * p, u + v * t)])
}

fn norm(u: &BigInt, v: &BigInt, p: &BigInt, t: &BigInt) -> BigInt {
    u * u + t * u * v + p * v * v
}

fn ideal_product(vector: &[i16; 13], records: &[GeneratorRecord]) -> Result<Hnf, String> {
    let p = bi(P_DEC);
    let t = bi(T_DEC);
    let mut accumulator = Hnf::identity();
    for (index, exponent) in vector.iter().copied().enumerate() {
        if index < RAMIFIED.len() && exponent < 0 {
            return Err("ramified exponents must be unsigned".to_owned());
        }
        let record = &records[index];
        let root = if exponent < 0 {
            record.negative_root
        } else {
            record.positive_root
        };
        let prime = Hnf::prime(record.ell, root);
        for _ in 0..exponent.unsigned_abs() {
            accumulator = multiply_hnf(&accumulator, &prime, &p, &t)?;
        }
    }
    Ok(accumulator)
}

fn controls(records: &[GeneratorRecord]) -> Result<Vec<ControlResult>, String> {
    let d = bi(D_DEC);
    let principal = BinaryQuadraticForm::principal(&d).unwrap();
    let p = bi(P_DEC);
    let t = bi(T_DEC);
    let mut out = Vec::new();
    for (index, (ell, _)) in RAMIFIED.iter().copied().enumerate() {
        let form = strings_form(&records[index].positive_form, &d)?;
        let form_ok = form.square() == principal;
        let mut vector = [0i16; 13];
        vector[index] = 2;
        let ideal = ideal_product(&vector, records)?;
        let alpha = BigInt::from(ell);
        let expected = principal_hnf(&alpha, &BigInt::zero(), &p, &t)?;
        out.push(ControlResult {
            name: format!("p192-ramified-{ell}-square"),
            passed: form_ok && ideal == expected,
            detail: format!("g{ell}^2 principal with alpha=({ell},0)"),
        });
    }

    // D=-23, pi^2-pi+6=0.  The oriented ideal (2,pi-1) cubed is
    // generated by alpha=1+pi, whose norm is 8.
    let toy_d = BigInt::from(-23);
    let toy_t = BigInt::one();
    let toy_p = BigInt::from(6u8);
    let toy_form = form_for_root(&toy_d, &toy_t, &toy_p, 2, 1)?;
    let toy_principal = BinaryQuadraticForm::principal(&toy_d).unwrap();
    let mut toy_ideal = Hnf::identity();
    for _ in 0..3 {
        toy_ideal = multiply_hnf(&toy_ideal, &Hnf::prime(2, 1), &toy_p, &toy_t)?;
    }
    let toy_alpha = principal_hnf(&BigInt::one(), &BigInt::one(), &toy_p, &toy_t)?;
    out.push(ControlResult {
        name: "d-minus-23-nonscalar".to_owned(),
        passed: toy_form.pow(3) == toy_principal
            && toy_ideal == toy_alpha
            && norm(&BigInt::one(), &BigInt::one(), &toy_p, &toy_t) == BigInt::from(8u8),
        detail: "g2^3 principal, alpha=(1,1) in {1,pi}, norm 8".to_owned(),
    });

    // Symbolic all-ramified relation, including the deliberately unbuildable
    // C-degree edge.  This is algebraic evidence only and cannot pass the map
    // or payoff gate.
    let c = bu(C_DEC).to_u64();
    let c_big = bu(C_DEC);
    let t_big = bu(T_DEC);
    let root_c_numerator = if t_big.is_odd() {
        &t_big + &c_big
    } else {
        t_big.clone()
    };
    let root_c = (root_c_numerator >> 1usize) % &c_big;
    let root_c_u = root_c.to_u64();
    // C does not fit u64, so construct its HNF directly with BigInts.
    let prime_c = Hnf {
        a: bi(C_DEC),
        b: (-BigInt::from_biguint(Sign::Plus, root_c)).mod_floor(&bi(C_DEC)),
        d: BigInt::one(),
    };
    let mut symbolic = Hnf::identity();
    for index in 0..3 {
        symbolic = multiply_hnf(
            &symbolic,
            &Hnf::prime(records[index].ell, records[index].positive_root),
            &p,
            &t,
        )?;
    }
    symbolic = multiply_hnf(&symbolic, &prime_c, &p, &t)?;
    let alpha_u = -t.clone();
    let alpha_v = BigInt::from(2u8);
    let symbolic_expected = principal_hnf(&alpha_u, &alpha_v, &p, &t)?;
    let symbolic_lambda = (&alpha_u + &alpha_v).mod_floor(&bi(N_DEC));
    out.push(ControlResult {
        name: "p192-symbolic-c-degree".to_owned(),
        passed: c.is_none()
            && root_c_u.is_none()
            && symbolic == symbolic_expected
            && norm(&alpha_u, &alpha_v, &p, &t) == d.abs()
            && symbolic_lambda == bi(SYMBOLIC_LAMBDA_DEC),
        detail: format!(
            "g5*g11*g31*gC principal with alpha=(-t,2), lambda={SYMBOLIC_LAMBDA_DEC}; C-degree map intentionally unavailable"
        ),
    });
    Ok(out)
}

pub fn build_generator_report(
    order: &OrderReport,
    max_ell: u64,
) -> Result<GeneratorReport, String> {
    if order.schema != SCHEMA_ORDER || order.status != "pass" || order.discriminant_d != D_DEC {
        return Err("generator construction requires a passing frozen order report".to_owned());
    }
    let certificate = build_pocklington_certificate()?;
    let expected_order = build_order_report(&certificate)?;
    if serde_json::to_value(order).map_err(|error| error.to_string())?
        != serde_json::to_value(expected_order).map_err(|error| error.to_string())?
    {
        return Err(
            "generator construction order input differs from native recomputation".to_owned(),
        );
    }
    if max_ell != 113 {
        return Err("the frozen experiment requires max-ell=113".to_owned());
    }
    let generators = generator_records()?;
    let controls = controls(&generators)?;
    let status = if generators.len() == 13
        && generators
            .iter()
            .all(|entry| entry.class_above_prime_crosscheck)
        && controls.iter().all(|control| control.passed)
    {
        "pass"
    } else {
        "invalid"
    };
    Ok(GeneratorReport {
        schema: SCHEMA_GENERATORS.to_owned(),
        curve: "NIST P-192 / secp192r1".to_owned(),
        discriminant_d: D_DEC.to_owned(),
        trace_t: T_DEC.to_owned(),
        max_ell,
        convention: "(ell,pi-r), smaller characteristic root positive; global conjugation fixes the first nonzero split exponent positive".to_owned(),
        generators,
        controls,
        status: status.to_owned(),
    })
}

pub fn verify_generator_report(report: &GeneratorReport) -> Result<(), String> {
    let certificate = build_pocklington_certificate()?;
    let order = build_order_report(&certificate)?;
    let expected = build_generator_report(&order, 113)?;
    let left = serde_json::to_value(report).map_err(|error| error.to_string())?;
    let right = serde_json::to_value(expected).map_err(|error| error.to_string())?;
    if left != right {
        return Err("generator report differs from native recomputation".to_owned());
    }
    if report.status != "pass" {
        return Err("generator report did not pass".to_owned());
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct WeightsReport {
    pub schema: String,
    pub status: String,
    pub records_required: u64,
    pub records_emitted: u64,
    pub seed: String,
    pub reason: String,
    pub hot: Option<serde_json::Value>,
    pub cold: Option<serde_json::Value>,
}

pub fn verify_unsupported_weights(report: &WeightsReport) -> Result<(), String> {
    let expected = unsupported_weights("0x50313932434d5231", 16, 11, 1);
    if report != &expected {
        return Err(
            "weights artifact is not the exact frozen unsupported-calibration record".to_owned(),
        );
    }
    Ok(())
}

pub fn unsupported_weights(
    seed: &str,
    endpoints: u64,
    repetitions: u64,
    warmup: u64,
) -> WeightsReport {
    WeightsReport {
        schema: SCHEMA_WEIGHTS.to_owned(),
        status: "unsupported_open".to_owned(),
        records_required: 23 * endpoints,
        records_emitted: 0,
        seed: seed.to_owned(),
        reason: format!(
            "explicit globally oriented P-192 map construction is not implemented; no HOT/COLD tuple is fabricated (requested warmup={warmup}, paired_repetitions={repetitions})"
        ),
        hot: None,
        cold: None,
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct SetDigest {
    xor: [u8; 32],
    sums: [u64; 4],
}

impl SetDigest {
    fn add(&mut self, vector: &[i16; 13]) {
        let mut hasher = Hasher::new();
        hasher.update(b"p192-cm-vector/v1\0");
        for exponent in vector {
            hasher.update(&exponent.to_le_bytes());
        }
        let digest = *hasher.finalize().as_bytes();
        for (left, right) in self.xor.iter_mut().zip(digest) {
            *left ^= right;
        }
        for (index, chunk) in digest.chunks_exact(8).enumerate() {
            self.sums[index] =
                self.sums[index].wrapping_add(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
    }

    fn render(&self) -> String {
        format!(
            "xor={};sum={:016x}{:016x}{:016x}{:016x}",
            hex::encode(self.xor),
            self.sums[0],
            self.sums[1],
            self.sums[2],
            self.sums[3]
        )
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EnumerationStats {
    pub oriented_vectors: u64,
    pub conjugation_fixed_vectors: u64,
    pub canonical_vectors_including_zero: u64,
    pub canonical_nonzero_vectors: u64,
    pub scalar_principal_vectors_including_zero: u64,
    pub exact_equal_2pow48_vectors: u64,
    pub exact_equal_2pow48_encoding_rule: String,
    pub exact_equal_2pow48_encoding_hex: String,
    pub exact_equal_2pow48_sha256: String,
    pub canonical_set_digest: String,
}

fn is_canonical(vector: &[i16; 13]) -> bool {
    for exponent in &vector[RAMIFIED.len()..] {
        match exponent.cmp(&0) {
            Ordering::Greater => return true,
            Ordering::Less => return false,
            Ordering::Equal => {}
        }
    }
    true
}

fn enumerate_all(order: &[usize; 13]) -> Result<EnumerationStats, String> {
    fn recurse(
        position: usize,
        order: &[usize; 13],
        degree: u64,
        vector: &mut [i16; 13],
        stats: &mut EnumerationStats,
        digest: &mut SetDigest,
    ) -> Result<(), String> {
        if position == order.len() {
            stats.oriented_vectors += 1;
            if vector[RAMIFIED.len()..]
                .iter()
                .all(|exponent| *exponent == 0)
            {
                stats.conjugation_fixed_vectors += 1;
                if vector[..RAMIFIED.len()]
                    .iter()
                    .all(|exponent| exponent % 2 == 0)
                {
                    stats.scalar_principal_vectors_including_zero += 1;
                }
            }
            if is_canonical(vector) {
                stats.canonical_vectors_including_zero += 1;
                if vector.iter().any(|exponent| *exponent != 0) {
                    stats.canonical_nonzero_vectors += 1;
                }
                if degree == MAX_DEGREE {
                    stats.exact_equal_2pow48_vectors += 1;
                }
                digest.add(vector);
                if stats.canonical_vectors_including_zero > HARD_STATE_CAP {
                    return Err("canonical state cap exceeded".to_owned());
                }
            }
            return Ok(());
        }
        let index = order[position];
        let ell = if index < RAMIFIED.len() {
            RAMIFIED[index].0
        } else {
            SPLIT[index - RAMIFIED.len()].0
        };
        let mut power = 1u64;
        let mut exponent = 0i16;
        loop {
            vector[index] = exponent;
            recurse(position + 1, order, degree * power, vector, stats, digest)?;
            if index >= RAMIFIED.len() && exponent > 0 {
                vector[index] = -exponent;
                recurse(position + 1, order, degree * power, vector, stats, digest)?;
            }
            if degree > MAX_DEGREE / power / ell {
                break;
            }
            power *= ell;
            exponent += 1;
        }
        vector[index] = 0;
        Ok(())
    }

    let mut stats = EnumerationStats {
        oriented_vectors: 0,
        conjugation_fixed_vectors: 0,
        canonical_vectors_including_zero: 0,
        canonical_nonzero_vectors: 0,
        scalar_principal_vectors_including_zero: 0,
        exact_equal_2pow48_vectors: 0,
        exact_equal_2pow48_encoding_rule: EMPTY_BOUNDARY_ENCODING_RULE.to_owned(),
        exact_equal_2pow48_encoding_hex: hex::encode(EMPTY_BOUNDARY_ENCODING),
        exact_equal_2pow48_sha256: hex::encode(sha256(EMPTY_BOUNDARY_ENCODING)),
        canonical_set_digest: String::new(),
    };
    let mut digest = SetDigest::default();
    recurse(0, order, 1, &mut [0i16; 13], &mut stats, &mut digest)?;
    stats.canonical_set_digest = digest.render();
    Ok(stats)
}

fn compute_exact_enumeration_stats() -> Result<(EnumerationStats, EnumerationStats), String> {
    let forward = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
    let reverse = [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0];
    let first = enumerate_all(&forward)?;
    let second = enumerate_all(&reverse)?;
    for stats in [&first, &second] {
        if stats.oriented_vectors != EXPECTED_ORIENTED
            || stats.conjugation_fixed_vectors != 635
            || stats.canonical_vectors_including_zero != EXPECTED_CANONICAL
            || stats.canonical_nonzero_vectors != EXPECTED_CANONICAL_NONZERO
            || stats.scalar_principal_vectors_including_zero != EXPECTED_SCALAR_INCLUDING_ZERO
            || stats.exact_equal_2pow48_vectors != 0
            || stats.exact_equal_2pow48_encoding_rule != EMPTY_BOUNDARY_ENCODING_RULE
            || stats.exact_equal_2pow48_encoding_hex != hex::encode(EMPTY_BOUNDARY_ENCODING)
            || stats.exact_equal_2pow48_sha256 != hex::encode(sha256(EMPTY_BOUNDARY_ENCODING))
        {
            return Err(format!(
                "frozen enumeration count mismatch: oriented={}, canonical={}, nonzero={}",
                stats.oriented_vectors,
                stats.canonical_vectors_including_zero,
                stats.canonical_nonzero_vectors
            ));
        }
    }
    if first != second {
        return Err("enumeration-order tail check changed the exact state set".to_owned());
    }
    Ok((first, second))
}

pub fn exact_enumeration_stats() -> Result<(EnumerationStats, EnumerationStats), String> {
    static CACHE: OnceLock<Result<(EnumerationStats, EnumerationStats), String>> = OnceLock::new();
    CACHE.get_or_init(compute_exact_enumeration_stats).clone()
}

#[derive(Clone)]
struct PowerChoice {
    exponent: i16,
    degree: u64,
    form: BinaryQuadraticForm,
}

fn power_choices(records: &[GeneratorRecord]) -> Result<Vec<Vec<PowerChoice>>, String> {
    let d = bi(D_DEC);
    let identity = BinaryQuadraticForm::principal(&d).unwrap();
    let mut all = Vec::with_capacity(records.len());
    for (index, record) in records.iter().enumerate() {
        let positive = strings_form(&record.positive_form, &d)?;
        let negative = strings_form(&record.negative_form, &d)?;
        let mut choices = vec![PowerChoice {
            exponent: 0,
            degree: 1,
            form: identity.clone(),
        }];
        let mut degree = 1u64;
        let mut exponent = 0i16;
        while degree <= MAX_DEGREE / record.ell {
            degree *= record.ell;
            exponent += 1;
            choices.push(PowerChoice {
                exponent,
                degree,
                form: positive.pow(exponent as u64),
            });
            if index >= RAMIFIED.len() {
                choices.push(PowerChoice {
                    exponent: -exponent,
                    degree,
                    form: negative.pow(exponent as u64),
                });
            }
        }
        all.push(choices);
    }
    Ok(all)
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RelationRecord {
    pub exponent_vector: Vec<i16>,
    pub oriented_ideals: Vec<String>,
    pub degree: String,
    pub log2_degree: f64,
    pub degree_linear_proxy: u64,
    pub reduced_class_product: [String; 3],
    pub classification: String,
    pub alpha_u: Option<String>,
    pub alpha_v: Option<String>,
    pub lambda_mod_n: Option<String>,
    pub alpha_norm: Option<String>,
    pub ideal_product_hnf: Option<IdealHnf>,
    pub principal_hnf: Option<IdealHnf>,
    pub hnf_verified: bool,
    pub abstract_candidate_only: bool,
    pub map_status: String,
    pub payoff_status: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct WordRadiusControl {
    pub requested_radius: u32,
    pub meet_in_middle_half_radius: u32,
    pub half_ball_states: u64,
    pub duplicate_class_count: u64,
    pub scalar_ramified_squares_quotiented: bool,
    pub passed_no_extra_relation: bool,
    pub method: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct FormKey {
    a: Vec<u8>,
    b: Vec<u8>,
    c: Vec<u8>,
}

fn form_key(form: &BinaryQuadraticForm) -> FormKey {
    let reduced = form.reduce();
    FormKey {
        a: reduced.a.to_signed_bytes_be(),
        b: reduced.b.to_signed_bytes_be(),
        c: reduced.c.to_signed_bytes_be(),
    }
}

fn compute_word_radius_control(
    records: &[GeneratorRecord],
    radius: u32,
) -> Result<WordRadiusControl, String> {
    if radius != 16 {
        return Err("the frozen word-radius control is exactly 16".to_owned());
    }
    let half = radius / 2;
    let d = bi(D_DEC);
    let principal = BinaryQuadraticForm::principal(&d).unwrap();
    let choices = power_choices(records)?;
    let mut seen = HashSet::<FormKey>::new();
    let mut states = 0u64;
    let mut duplicates = 0u64;

    fn recurse(
        index: usize,
        remaining: u32,
        form: &BinaryQuadraticForm,
        choices: &[Vec<PowerChoice>],
        seen: &mut HashSet<FormKey>,
        states: &mut u64,
        duplicates: &mut u64,
    ) -> Result<(), String> {
        if index == 13 {
            *states += 1;
            if *states > 10_000_000 {
                return Err("word-radius half-ball unexpectedly exceeded 10M states".to_owned());
            }
            if !seen.insert(form_key(form)) {
                *duplicates += 1;
            }
            return Ok(());
        }
        for choice in &choices[index] {
            let cost = choice.exponent.unsigned_abs() as u32;
            if cost > remaining {
                continue;
            }
            // A ramified class is self-inverse and its square is the known
            // scalar relation, so its exact Cayley coordinate is one bit.
            if index < RAMIFIED.len() && !(choice.exponent == 0 || choice.exponent == 1) {
                continue;
            }
            let next = if choice.exponent == 0 {
                form.clone()
            } else {
                form.compose(&choice.form)
            };
            recurse(
                index + 1,
                remaining - cost,
                &next,
                choices,
                seen,
                states,
                duplicates,
            )?;
        }
        Ok(())
    }

    recurse(
        0,
        half,
        &principal,
        &choices,
        &mut seen,
        &mut states,
        &mut duplicates,
    )?;
    Ok(WordRadiusControl {
        requested_radius: radius,
        meet_in_middle_half_radius: half,
        half_ball_states: states,
        duplicate_class_count: duplicates,
        scalar_ramified_squares_quotiented: true,
        passed_no_extra_relation: duplicates == 0,
        method: "exact collision test on the complete radius-8 Cayley half-ball; a relation of length <=16 splits into two half-ball words".to_owned(),
    })
}

fn word_radius_control(
    records: &[GeneratorRecord],
    radius: u32,
) -> Result<WordRadiusControl, String> {
    if radius != 16 {
        return Err("the frozen word-radius control is exactly 16".to_owned());
    }
    static CACHE: OnceLock<Result<WordRadiusControl, String>> = OnceLock::new();
    CACHE
        .get_or_init(|| compute_word_radius_control(records, radius))
        .clone()
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SearchReport {
    pub schema: String,
    pub curve: String,
    pub discriminant_d: String,
    pub boundary: String,
    pub max_ell: u64,
    pub max_states: u64,
    pub word_radius_control: u32,
    pub enumeration_forward: EnumerationStats,
    pub enumeration_reverse: EnumerationStats,
    pub zero_vectors_excluded: u64,
    pub inverse_or_dual_backtracks_excluded: String,
    pub principal_relations: Vec<RelationRecord>,
    pub scalar_ramified_states_including_zero: u64,
    pub scalar_ramified_relation_count: u64,
    pub non_scalar_relation_count: u64,
    pub word_radius_control_result: WordRadiusControl,
    pub theorem_control: ControlResult,
    pub controls: Vec<ControlResult>,
    pub empirical_weights: WeightsReport,
    pub empirical_weights_status: String,
    pub exact_search_status: String,
    pub experiment_status: String,
    pub claim_ceiling: String,
}

fn vector_degree(vector: &[i16; 13], records: &[GeneratorRecord]) -> Result<u64, String> {
    let mut degree = 1u64;
    for (exponent, record) in vector.iter().zip(records) {
        for _ in 0..exponent.unsigned_abs() {
            degree = degree.checked_mul(record.ell).ok_or("degree overflow")?;
        }
    }
    Ok(degree)
}

fn expected_scalar_vectors(records: &[GeneratorRecord]) -> Result<Vec<[i16; 13]>, String> {
    let mut vectors = Vec::new();
    let mut e5 = 0i16;
    let mut d5 = 1u64;
    loop {
        let mut e11 = 0i16;
        let mut d11 = d5;
        loop {
            let mut e31 = 0i16;
            let mut degree = d11;
            loop {
                if e5 != 0 || e11 != 0 || e31 != 0 {
                    let mut vector = [0i16; 13];
                    vector[0] = e5;
                    vector[1] = e11;
                    vector[2] = e31;
                    if vector_degree(&vector, records)? != degree {
                        return Err("scalar-vector degree construction mismatch".to_owned());
                    }
                    vectors.push(vector);
                }
                let ell_squared = records[2]
                    .ell
                    .checked_mul(records[2].ell)
                    .ok_or("ell squared overflow")?;
                if degree > MAX_DEGREE / ell_squared {
                    break;
                }
                degree *= ell_squared;
                e31 += 2;
            }
            let ell_squared = records[1]
                .ell
                .checked_mul(records[1].ell)
                .ok_or("ell squared overflow")?;
            if d11 > MAX_DEGREE / ell_squared {
                break;
            }
            d11 *= ell_squared;
            e11 += 2;
        }
        let ell_squared = records[0]
            .ell
            .checked_mul(records[0].ell)
            .ok_or("ell squared overflow")?;
        if d5 > MAX_DEGREE / ell_squared {
            break;
        }
        d5 *= ell_squared;
        e5 += 2;
    }
    if vectors.len() as u64 != EXPECTED_SCALAR_NONZERO {
        return Err(format!(
            "derived {} scalar nonzero vectors, expected {EXPECTED_SCALAR_NONZERO}",
            vectors.len()
        ));
    }
    Ok(vectors)
}

#[cfg(test)]
fn expected_scalar_relations(records: &[GeneratorRecord]) -> Result<Vec<RelationRecord>, String> {
    let principal = BinaryQuadraticForm::principal(&bi(D_DEC)).unwrap();
    expected_scalar_vectors(records)?
        .iter()
        .map(|vector| relation_record(vector, &principal, records))
        .collect()
}

fn verify_relation_inventory(
    relations: &[RelationRecord],
    records: &[GeneratorRecord],
) -> Result<(), String> {
    if relations.len() as u64 != EXPECTED_SCALAR_NONZERO {
        return Err(format!(
            "relation inventory has {} entries, expected {EXPECTED_SCALAR_NONZERO}",
            relations.len()
        ));
    }
    let expected: HashSet<[i16; 13]> = expected_scalar_vectors(records)?.into_iter().collect();
    let mut seen = HashSet::with_capacity(relations.len());
    for relation in relations {
        verify_relation(relation, records)?;
        if relation.classification != "scalar-ramified" || !relation.hnf_verified {
            return Err("relation inventory contains a non-scalar or unverified entry".to_owned());
        }
        let vector = array13(&relation.exponent_vector)?;
        if !seen.insert(vector) {
            return Err("relation inventory contains a duplicate exponent vector".to_owned());
        }
    }
    if seen != expected {
        return Err("relation inventory is not the complete frozen scalar set".to_owned());
    }
    Ok(())
}

fn relation_record(
    vector: &[i16; 13],
    form: &BinaryQuadraticForm,
    records: &[GeneratorRecord],
) -> Result<RelationRecord, String> {
    let degree = vector_degree(vector, records)?;
    let p = bi(P_DEC);
    let t = bi(T_DEC);
    let n = bi(N_DEC);
    let ideal = ideal_product(vector, records)?;
    let sqrt = (degree as f64).sqrt() as u64;
    let exact_sqrt = if sqrt.checked_mul(sqrt) == Some(degree) {
        Some(sqrt)
    } else if (sqrt + 1).checked_mul(sqrt + 1) == Some(degree) {
        Some(sqrt + 1)
    } else {
        None
    };
    let split_nonzero = vector[RAMIFIED.len()..]
        .iter()
        .any(|exponent| *exponent != 0);
    let scalar = !split_nonzero
        && vector[..RAMIFIED.len()]
            .iter()
            .all(|exponent| exponent % 2 == 0)
        && exact_sqrt.is_some();
    let (classification, alpha_u, alpha_v, lambda, alpha_norm, principal, hnf_verified) = if scalar
    {
        let u = BigInt::from(exact_sqrt.unwrap());
        let v = BigInt::zero();
        let principal = principal_hnf(&u, &v, &p, &t)?;
        let norm_value = norm(&u, &v, &p, &t);
        (
            "scalar-ramified".to_owned(),
            Some(u.to_string()),
            Some(v.to_string()),
            Some(u.mod_floor(&n).to_string()),
            Some(norm_value.to_string()),
            Some(principal.clone()),
            ideal == principal,
        )
    } else {
        (
            "non-scalar-abstract".to_owned(),
            None,
            None,
            None,
            None,
            None,
            false,
        )
    };
    let oriented_ideals = vector
        .iter()
        .zip(records)
        .filter(|(exponent, _)| **exponent != 0)
        .map(|(exponent, record)| {
            let root = if *exponent < 0 {
                record.negative_root
            } else {
                record.positive_root
            };
            format!("({},{},pi-{})", exponent, record.ell, root)
        })
        .collect();
    let degree_linear_proxy = vector
        .iter()
        .zip(records)
        .map(|(exponent, record)| exponent.unsigned_abs() as u64 * record.ell)
        .sum();
    Ok(RelationRecord {
        exponent_vector: vector.to_vec(),
        oriented_ideals,
        degree: degree.to_string(),
        log2_degree: (degree as f64).log2(),
        degree_linear_proxy,
        reduced_class_product: form_strings(form),
        classification,
        alpha_u,
        alpha_v,
        lambda_mod_n: lambda,
        alpha_norm,
        ideal_product_hnf: Some(ideal.record()),
        principal_hnf: principal.map(|entry| entry.record()),
        hnf_verified,
        abstract_candidate_only: true,
        map_status: "not-built".to_owned(),
        payoff_status: "not-evaluable-without-HOT-COLD-map-costs".to_owned(),
    })
}

fn search_principal_relations(
    records: &[GeneratorRecord],
    max_states: u64,
) -> Result<Vec<RelationRecord>, String> {
    struct SearchState<'a> {
        records: &'a [GeneratorRecord],
        choices: Vec<Vec<PowerChoice>>,
        principal: BinaryQuadraticForm,
        relations: Vec<RelationRecord>,
        states: u64,
        max_states: u64,
    }

    fn recurse(
        index: usize,
        degree: u64,
        first_split_seen: bool,
        form: &BinaryQuadraticForm,
        vector: &mut [i16; 13],
        state: &mut SearchState<'_>,
    ) -> Result<(), String> {
        if index == 13 {
            state.states += 1;
            if state.states > state.max_states {
                return Err(format!("hard state cap {} exceeded", state.max_states));
            }
            if vector.iter().all(|exponent| *exponent == 0) {
                return Ok(());
            }
            if form == &state.principal {
                let relation = relation_record(vector, form, state.records)?;
                state.relations.push(relation);
            }
            return Ok(());
        }
        for choice in state.choices[index].clone() {
            if degree > MAX_DEGREE / choice.degree {
                continue;
            }
            if index < RAMIFIED.len() && choice.exponent < 0 {
                continue;
            }
            if index >= RAMIFIED.len() && !first_split_seen && choice.exponent < 0 {
                continue;
            }
            vector[index] = choice.exponent;
            let next_seen = first_split_seen || (index >= RAMIFIED.len() && choice.exponent != 0);
            let next_form = if choice.exponent == 0 {
                form.clone()
            } else {
                form.compose(&choice.form)
            };
            recurse(
                index + 1,
                degree * choice.degree,
                next_seen,
                &next_form,
                vector,
                state,
            )?;
        }
        vector[index] = 0;
        Ok(())
    }

    let d = bi(D_DEC);
    let principal = BinaryQuadraticForm::principal(&d).unwrap();
    let mut state = SearchState {
        records,
        choices: power_choices(records)?,
        principal: principal.clone(),
        relations: Vec::new(),
        states: 0,
        max_states,
    };
    recurse(0, 1, false, &principal, &mut [0i16; 13], &mut state)?;
    if state.states != EXPECTED_CANONICAL {
        return Err(format!(
            "form search visited {} canonical states, expected {}",
            state.states, EXPECTED_CANONICAL
        ));
    }
    Ok(state.relations)
}

pub fn run_exact_search(
    generators: &GeneratorReport,
    weights: Option<&WeightsReport>,
    max_states: u64,
    word_radius: u32,
) -> Result<SearchReport, String> {
    verify_generator_report(generators)?;
    if max_states != HARD_STATE_CAP {
        return Err(format!(
            "frozen exact search requires max-states={HARD_STATE_CAP}"
        ));
    }
    if word_radius != 16 {
        return Err("frozen exact search requires word-radius-control=16".to_owned());
    }
    let (forward, reverse) = exact_enumeration_stats()?;
    let relations = search_principal_relations(&generators.generators, max_states)?;
    let scalar_count = relations
        .iter()
        .filter(|entry| entry.classification == "scalar-ramified")
        .count() as u64;
    let non_scalar_count = relations.len() as u64 - scalar_count;
    if relations
        .iter()
        .any(|entry| entry.classification == "scalar-ramified" && !entry.hnf_verified)
    {
        return Err("a scalar relation failed independent ideal-HNF replay".to_owned());
    }
    if scalar_count != EXPECTED_SCALAR_NONZERO {
        return Err(format!(
            "found {scalar_count} nonzero scalar relations, expected {EXPECTED_SCALAR_NONZERO}"
        ));
    }
    if non_scalar_count != 0 {
        return Err(format!(
            "found {non_scalar_count} unexpected principal non-scalar vectors below the certified norm bound; premises are inconsistent and no incomplete alpha/HNF record will be emitted"
        ));
    }
    verify_relation_inventory(&relations, &generators.generators)?;
    let radius_control = word_radius_control(&generators.generators, word_radius)?;
    if !radius_control.passed_no_extra_relation {
        return Err(format!(
            "word-radius control found {} unexpected half-ball class collisions",
            radius_control.duplicate_class_count
        ));
    }
    let weights = weights.ok_or_else(|| {
        "frozen exact search requires the explicit unsupported-calibration artifact".to_owned()
    })?;
    verify_unsupported_weights(weights)?;
    let empirical_status = weights.status.clone();
    Ok(SearchReport {
        schema: SCHEMA_SEARCH.to_owned(),
        curve: "NIST P-192 / secp192r1".to_owned(),
        discriminant_d: D_DEC.to_owned(),
        boundary: "product ell^|e_ell| <= 2^48".to_owned(),
        max_ell: 113,
        max_states,
        word_radius_control: word_radius,
        enumeration_forward: forward,
        enumeration_reverse: reverse,
        zero_vectors_excluded: 1,
        inverse_or_dual_backtracks_excluded:
            "net exponent-vector encoding; zero vector separately excluded".to_owned(),
        principal_relations: relations,
        scalar_ramified_states_including_zero: EXPECTED_SCALAR_INCLUDING_ZERO,
        scalar_ramified_relation_count: scalar_count,
        non_scalar_relation_count: non_scalar_count,
        word_radius_control_result: radius_control,
        theorem_control: expected_theorem_control(),
        controls: generators.controls.clone(),
        empirical_weights: weights.clone(),
        empirical_weights_status: empirical_status,
        exact_search_status: "complete-exact-algebra-box".to_owned(),
        experiment_status: EXPERIMENT_INCOMPLETE_STATUS.to_owned(),
        claim_ceiling: if non_scalar_count == 0 {
            "NO_WEAKNESS_FOUND_WITHIN_SCOPE"
        } else {
            "CANDIDATE_ONLY_REQUIRES_MAP_AND_PAYOFF_REPLAY"
        }
        .to_owned(),
    })
}

fn array13(values: &[i16]) -> Result<[i16; 13], String> {
    values
        .try_into()
        .map_err(|_| "exponent vector must contain exactly 13 entries".to_owned())
}

fn verify_relation(entry: &RelationRecord, records: &[GeneratorRecord]) -> Result<(), String> {
    let vector = array13(&entry.exponent_vector)?;
    if !is_canonical(&vector) || vector.iter().all(|exponent| *exponent == 0) {
        return Err("relation vector is zero or not conjugacy-canonical".to_owned());
    }
    let degree = vector_degree(&vector, records)?;
    if entry.degree != degree.to_string() || degree > MAX_DEGREE {
        return Err("relation degree mismatch".to_owned());
    }
    let d = bi(D_DEC);
    let principal = BinaryQuadraticForm::principal(&d).unwrap();
    let mut product = principal.clone();
    for (index, exponent) in vector.iter().copied().enumerate() {
        let record = &records[index];
        let form = if exponent < 0 {
            strings_form(&record.negative_form, &d)?
        } else {
            strings_form(&record.positive_form, &d)?
        };
        product = product.compose(&form.pow(exponent.unsigned_abs() as u64));
    }
    if product != principal || entry.reduced_class_product != form_strings(&product) {
        return Err("relation form replay is non-principal or mismatched".to_owned());
    }
    let recomputed = relation_record(&vector, &product, records)?;
    let left = serde_json::to_value(entry).map_err(|error| error.to_string())?;
    let right = serde_json::to_value(recomputed).map_err(|error| error.to_string())?;
    if left != right {
        return Err("relation certificate differs from independent replay".to_owned());
    }
    Ok(())
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VerificationReport {
    pub schema: String,
    pub certificate_schema: String,
    pub exact_enumeration_replayed: bool,
    pub independent_ideal_hnf_replayed: bool,
    pub relation_count: u64,
    pub rejected_relation_count: u64,
    pub non_scalar_relation_count: u64,
    pub map_replay_status: String,
    pub calibration_status: String,
    pub payoff_gate_passed: bool,
    pub experiment_status: String,
    pub verdict: String,
    pub claim_ceiling: String,
}

fn verify_search_envelope(
    report: &SearchReport,
    records: &[GeneratorRecord],
) -> Result<(), String> {
    if report.schema != SCHEMA_SEARCH
        || report.curve != "NIST P-192 / secp192r1"
        || report.discriminant_d != D_DEC
        || report.boundary != "product ell^|e_ell| <= 2^48"
        || report.max_ell != 113
        || report.max_states != HARD_STATE_CAP
        || report.word_radius_control != 16
        || report.zero_vectors_excluded != 1
        || report.inverse_or_dual_backtracks_excluded
            != "net exponent-vector encoding; zero vector separately excluded"
        || report.scalar_ramified_states_including_zero != EXPECTED_SCALAR_INCLUDING_ZERO
        || report.scalar_ramified_relation_count != EXPECTED_SCALAR_NONZERO
        || report.non_scalar_relation_count != 0
        || report.exact_search_status != "complete-exact-algebra-box"
        || report.experiment_status != EXPERIMENT_INCOMPLETE_STATUS
        || report.claim_ceiling != "NO_WEAKNESS_FOUND_WITHIN_SCOPE"
    {
        return Err("one or more frozen search envelope fields differ".to_owned());
    }
    if report.controls != controls(records)? || report.theorem_control != expected_theorem_control()
    {
        return Err("top-level controls differ from native recomputation".to_owned());
    }
    verify_unsupported_weights(&report.empirical_weights)?;
    if report.empirical_weights_status != report.empirical_weights.status
        || report.empirical_weights_status != "unsupported_open"
    {
        return Err("empirical weight status is unbound or fabricated".to_owned());
    }
    Ok(())
}

pub fn verify_search_report(
    report: &SearchReport,
    independent_hnf: bool,
    replay_maps: bool,
) -> Result<VerificationReport, String> {
    if !independent_hnf {
        return Err("frozen verification requires --independent-ideal-hnf".to_owned());
    }
    let pocklington = build_pocklington_certificate()?;
    let order = build_order_report(&pocklington)?;
    verify_order_report(&order, &pocklington)?;
    let generators = build_generator_report(&order, 113)?;
    verify_generator_report(&generators)?;
    let records = generators.generators;
    verify_search_envelope(report, &records)?;
    let (forward, reverse) = exact_enumeration_stats()?;
    if report.enumeration_forward != forward || report.enumeration_reverse != reverse {
        return Err("search enumeration manifest failed independent replay".to_owned());
    }
    verify_relation_inventory(&report.principal_relations, &records)?;
    let radius_control = word_radius_control(&records, report.word_radius_control)?;
    if report.word_radius_control_result != radius_control {
        return Err("word-radius relation control failed independent replay".to_owned());
    }
    let non_scalar = 0;
    let map_status = if !replay_maps {
        "not-requested"
    } else if non_scalar == 0 {
        "not-applicable-no-nonscalar-hit"
    } else {
        "unsupported-open"
    };
    Ok(VerificationReport {
        schema: SCHEMA_VERIFY.to_owned(),
        certificate_schema: report.schema.clone(),
        exact_enumeration_replayed: true,
        independent_ideal_hnf_replayed: true,
        relation_count: report.principal_relations.len() as u64,
        rejected_relation_count: 0,
        non_scalar_relation_count: non_scalar,
        map_replay_status: map_status.to_owned(),
        calibration_status: report.empirical_weights_status.clone(),
        payoff_gate_passed: false,
        experiment_status: EXPERIMENT_INCOMPLETE_STATUS.to_owned(),
        verdict: if non_scalar == 0 {
            VERIFICATION_NO_HIT_VERDICT
        } else {
            VERIFICATION_CANDIDATE_VERDICT
        }
        .to_owned(),
        claim_ceiling: "NO_WEAKNESS_FOUND_WITHIN_SCOPE".to_owned(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixtures() -> (PocklingtonCertificate, OrderReport, GeneratorReport) {
        let certificate = build_pocklington_certificate().unwrap();
        let order = build_order_report(&certificate).unwrap();
        let generators = build_generator_report(&order, 113).unwrap();
        (certificate, order, generators)
    }

    fn search_report_fixture() -> SearchReport {
        let (_, _, generators) = fixtures();
        let (forward, reverse) = exact_enumeration_stats().unwrap();
        let weights = unsupported_weights("0x50313932434d5231", 16, 11, 1);
        SearchReport {
            schema: SCHEMA_SEARCH.to_owned(),
            curve: "NIST P-192 / secp192r1".to_owned(),
            discriminant_d: D_DEC.to_owned(),
            boundary: "product ell^|e_ell| <= 2^48".to_owned(),
            max_ell: 113,
            max_states: HARD_STATE_CAP,
            word_radius_control: 16,
            enumeration_forward: forward,
            enumeration_reverse: reverse,
            zero_vectors_excluded: 1,
            inverse_or_dual_backtracks_excluded:
                "net exponent-vector encoding; zero vector separately excluded".to_owned(),
            principal_relations: expected_scalar_relations(&generators.generators).unwrap(),
            scalar_ramified_states_including_zero: EXPECTED_SCALAR_INCLUDING_ZERO,
            scalar_ramified_relation_count: EXPECTED_SCALAR_NONZERO,
            non_scalar_relation_count: 0,
            word_radius_control_result: WordRadiusControl {
                requested_radius: 16,
                meet_in_middle_half_radius: 8,
                half_ball_states: 0,
                duplicate_class_count: 0,
                scalar_ramified_squares_quotiented: true,
                passed_no_extra_relation: true,
                method: "fixture-not-used-before-inventory-rejection".to_owned(),
            },
            theorem_control: expected_theorem_control(),
            controls: generators.controls,
            empirical_weights_status: weights.status.clone(),
            empirical_weights: weights,
            exact_search_status: "complete-exact-algebra-box".to_owned(),
            experiment_status: EXPERIMENT_INCOMPLETE_STATUS.to_owned(),
            claim_ceiling: "NO_WEAKNESS_FOUND_WITHIN_SCOPE".to_owned(),
        }
    }

    #[test]
    fn pocklington_tree_and_order_certificate_replay() {
        let (certificate, order, _) = fixtures();
        verify_pocklington_certificate(&certificate).unwrap();
        verify_order_report(&order, &certificate).unwrap();
        assert!(certificate
            .root
            .recursive_children
            .iter()
            .any(|child| child.n == Q_DEC));
        assert!(order.d_fundamental);
        assert_eq!(order.frobenius_conductor, 1);
    }

    #[test]
    fn pocklington_mutation_is_rejected() {
        let mut certificate = build_pocklington_certificate().unwrap();
        certificate.root.witnesses[0].a = "1".to_owned();
        assert!(verify_pocklington_certificate(&certificate).is_err());
    }

    #[test]
    fn generator_roots_forms_and_controls_replay() {
        let (_, _, generators) = fixtures();
        verify_generator_report(&generators).unwrap();
        assert_eq!(generators.generators.len(), 13);
        assert!(generators.controls.iter().all(|control| control.passed));
    }

    #[test]
    fn root_and_discriminant_mutations_are_rejected() {
        let (certificate, mut order, _) = fixtures();
        order.p = "17".to_owned();
        assert!(verify_order_report(&order, &certificate).is_err());
        assert!(build_generator_report(&order, 113).is_err());

        let (_, _, mut generators) = fixtures();
        generators.generators[3].positive_root += 1;
        assert!(verify_generator_report(&generators).is_err());
        let (_, _, mut generators) = fixtures();
        generators.discriminant_d = "-24".to_owned();
        assert!(verify_generator_report(&generators).is_err());
    }

    #[test]
    fn exact_boundary_counts_match_frozen_values_in_both_orders() {
        let (forward, reverse) = exact_enumeration_stats().unwrap();
        assert_eq!(forward, reverse);
        assert_eq!(forward.oriented_vectors, EXPECTED_ORIENTED);
        assert_eq!(forward.conjugation_fixed_vectors, 635);
        assert_eq!(forward.canonical_vectors_including_zero, EXPECTED_CANONICAL);
        assert_eq!(
            forward.canonical_nonzero_vectors,
            EXPECTED_CANONICAL_NONZERO
        );
        assert_eq!(forward.exact_equal_2pow48_vectors, 0);
        assert_eq!(
            forward.scalar_principal_vectors_including_zero,
            EXPECTED_SCALAR_INCLUDING_ZERO
        );
        assert_eq!(
            forward.exact_equal_2pow48_sha256,
            hex::encode(sha256(EMPTY_BOUNDARY_ENCODING))
        );
        assert_eq!(
            forward.exact_equal_2pow48_sha256,
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
    }

    #[test]
    fn full_exact_search_reaches_frozen_relation_inventory() {
        let (_, _, generators) = fixtures();
        let weights = unsupported_weights("0x50313932434d5231", 16, 11, 1);
        let report = run_exact_search(&generators, Some(&weights), HARD_STATE_CAP, 16).unwrap();
        assert_eq!(
            report.enumeration_forward.canonical_vectors_including_zero,
            EXPECTED_CANONICAL
        );
        assert_eq!(
            report.scalar_ramified_relation_count,
            EXPECTED_SCALAR_NONZERO
        );
        assert_eq!(report.non_scalar_relation_count, 0);
        assert_eq!(report.exact_search_status, "complete-exact-algebra-box");
        verify_search_report(&report, true, true).unwrap();
    }

    #[test]
    fn hnf_controls_reject_alpha_and_root_mutations() {
        let (_, _, generators) = fixtures();
        let mut vector = [0i16; 13];
        vector[0] = 2;
        let got = ideal_product(&vector, &generators.generators).unwrap();
        let expected =
            principal_hnf(&BigInt::from(5u8), &BigInt::zero(), &bi(P_DEC), &bi(T_DEC)).unwrap();
        assert_eq!(got, expected);
        let wrong_alpha =
            principal_hnf(&BigInt::from(6u8), &BigInt::zero(), &bi(P_DEC), &bi(T_DEC)).unwrap();
        assert_ne!(got, wrong_alpha);
        let wrong_root = Hnf::prime(5, 3);
        let squared = multiply_hnf(&wrong_root, &wrong_root, &bi(P_DEC), &bi(T_DEC)).unwrap();
        assert_ne!(squared, expected);
    }

    #[test]
    fn unsupported_calibration_is_explicit() {
        let weights = unsupported_weights("0x50313932434d5231", 16, 11, 1);
        assert_eq!(weights.status, "unsupported_open");
        assert_eq!(weights.records_required, 368);
        assert_eq!(weights.records_emitted, 0);
        assert!(weights.hot.is_none() && weights.cold.is_none());
    }

    #[test]
    fn fabricated_calibration_is_rejected() {
        let mut weights = unsupported_weights("0x50313932434d5231", 16, 11, 1);
        weights.status = "pass".to_owned();
        weights.records_emitted = 368;
        assert!(verify_unsupported_weights(&weights).is_err());
    }

    #[test]
    fn relation_inventory_rejects_subsets_and_duplicates() {
        let mut subset = search_report_fixture();
        subset.principal_relations.pop();
        assert!(verify_search_report(&subset, true, true).is_err());

        let mut duplicate = search_report_fixture();
        let first = duplicate.principal_relations[0].clone();
        *duplicate.principal_relations.last_mut().unwrap() = first;
        assert!(verify_search_report(&duplicate, true, true).is_err());
    }

    #[test]
    fn verifier_binds_envelope_controls_and_boundary_digest() {
        let records = generator_records().unwrap();
        let mut report = search_report_fixture();
        report.curve = "P-192-ish".to_owned();
        assert!(verify_search_envelope(&report, &records).is_err());

        let mut report = search_report_fixture();
        report.controls[0].detail.push_str(" mutated");
        assert!(verify_search_envelope(&report, &records).is_err());

        let mut report = search_report_fixture();
        report.enumeration_reverse.exact_equal_2pow48_sha256 = "00".repeat(32);
        assert!(verify_search_report(&report, true, true).is_err());

        let mut report = search_report_fixture();
        report.scalar_ramified_states_including_zero = 103;
        assert!(verify_search_envelope(&report, &records).is_err());
    }

    #[test]
    fn verification_labels_preserve_the_open_experiment_gate() {
        assert_eq!(
            EXPERIMENT_INCOMPLETE_STATUS,
            "incomplete-map-calibration-and-payoff-open"
        );
        assert!(VERIFICATION_NO_HIT_VERDICT.starts_with("incomplete-"));
        assert!(VERIFICATION_CANDIDATE_VERDICT.starts_with("incomplete-"));
    }
}
