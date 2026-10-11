use std::cmp::Ordering;

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::{params, Result};

const CAND_DOMAIN: &[u8] = b"P192-WCM-CAND-v1\0";
const FACT_DOMAIN: &[u8] = b"P192-WCM-FACT-v1\0";
const CANDIDATE_SHARD_DOMAIN: &[u8] = b"P192-WCM-CANDIDATE-SHARD-v1\0";
const SHELL_ROOT_DOMAIN: &[u8] = b"P192-WCM-SHELL-ROOT-v1\0";
const DISPOSITION_SHARD_DOMAIN: &[u8] = b"P192-WCM-DISPOSITION-SHARD-v1\0";
const DISPOSITION_SHELL_ROOT_DOMAIN: &[u8] = b"P192-WCM-DISPOSITION-SHELL-ROOT-v1\0";

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RationalFactor {
    pub prime: BigUint,
    pub exponent: BigUint,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmallEntry {
    pub ell: u32,
    /// Zero is split and one is ramified.
    pub kind: u8,
    pub exponent: BigInt,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LargePrimeEntry {
    pub q: u32,
    pub smaller_root: u32,
    pub larger_root: u32,
    /// One is the smaller/positive root and two is the larger/negative root.
    pub sign: u8,
    pub multiplicity: BigUint,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CandidateCertificate {
    /// Zero complete, one one-LP, two two-LP.
    pub candidate_type: u8,
    pub curve_uid: String,
    pub factor_base_sha256: [u8; 32],
    pub shell_id: String,
    pub candidate_index: u64,
    pub v_box: u64,
    pub x_box: u64,
    pub u: BigInt,
    pub v: BigInt,
    pub norm: BigUint,
    pub scalar_content: BigUint,
    pub rational_factors: Vec<RationalFactor>,
    pub small_entries: Vec<SmallEntry>,
    pub lp_entries: Vec<LargePrimeEntry>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FactorizationCertificate {
    pub candidate_sha256: [u8; 32],
    pub norm: BigUint,
    pub rational_factors: Vec<RationalFactor>,
    pub residual: BigUint,
    pub residual_type: u8,
    pub lp_entries: Vec<LargePrimeEntry>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShardDigest {
    pub shard_index: u64,
    pub first_candidate_index: u64,
    pub record_count: u64,
    pub sha256: [u8; 32],
}

pub fn sha256(bytes: &[u8]) -> [u8; 32] {
    Sha256::digest(bytes).into()
}

pub fn sha256_hex(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

fn push_len_prefixed(out: &mut Vec<u8>, bytes: &[u8]) -> Result<()> {
    let len = u32::try_from(bytes.len()).map_err(|_| "length does not fit u32".to_owned())?;
    out.extend_from_slice(&len.to_be_bytes());
    out.extend_from_slice(bytes);
    Ok(())
}

pub fn encode_shell_id(shell_id: &str) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(4 + shell_id.len());
    push_len_prefixed(&mut out, shell_id.as_bytes())?;
    Ok(out)
}

pub fn encode_nat(value: &BigUint, out: &mut Vec<u8>) -> Result<()> {
    let bytes = value.to_bytes_be();
    let len = u32::try_from(bytes.len()).map_err(|_| "nat length does not fit u32".to_owned())?;
    out.extend_from_slice(&len.to_be_bytes());
    out.extend_from_slice(&bytes);
    Ok(())
}

pub fn encode_sint(value: &BigInt, out: &mut Vec<u8>) -> Result<()> {
    let (sign, magnitude) = value.to_bytes_be();
    out.push(match sign {
        Sign::NoSign => 0,
        Sign::Plus => 1,
        Sign::Minus => 2,
    });
    encode_nat(&BigUint::from_bytes_be(&magnitude), out)
}

fn encode_rational_factors(factors: &[RationalFactor], out: &mut Vec<u8>) -> Result<()> {
    validate_rational_factors(factors)?;
    for factor in factors {
        encode_nat(&factor.prime, out)?;
        encode_nat(&factor.exponent, out)?;
    }
    Ok(())
}

fn encode_lp_entries(entries: &[LargePrimeEntry], out: &mut Vec<u8>) -> Result<()> {
    validate_lp_entries(entries)?;
    for entry in entries {
        out.extend_from_slice(&entry.q.to_be_bytes());
        out.extend_from_slice(&entry.smaller_root.to_be_bytes());
        out.extend_from_slice(&entry.larger_root.to_be_bytes());
        out.push(entry.sign);
        encode_nat(&entry.multiplicity, out)?;
    }
    Ok(())
}

fn validate_rational_factors(factors: &[RationalFactor]) -> Result<()> {
    let mut previous: Option<&BigUint> = None;
    for factor in factors {
        if factor.prime < BigUint::from(2u8) || factor.exponent.is_zero() {
            return Err("invalid rational factor".to_owned());
        }
        if previous.is_some_and(|p| p >= &factor.prime) {
            return Err("rational factors are not strictly increasing".to_owned());
        }
        previous = Some(&factor.prime);
    }
    Ok(())
}

fn validate_small_entries(entries: &[SmallEntry]) -> Result<()> {
    let mut previous = 0u32;
    for entry in entries {
        if entry.ell <= previous || entry.kind > 1 || entry.exponent.is_zero() {
            return Err("invalid or unordered small entry".to_owned());
        }
        if entry.kind == 1 && entry.exponent != BigInt::from(1u8) {
            return Err("ramified small exponent is not the canonical bit one".to_owned());
        }
        previous = entry.ell;
    }
    Ok(())
}

fn validate_lp_entries(entries: &[LargePrimeEntry]) -> Result<()> {
    let mut previous = 0u32;
    for entry in entries {
        if entry.q <= previous
            || entry.smaller_root >= entry.larger_root
            || entry.larger_root >= entry.q
            || u64::from(entry.q) <= params::ALGEBRAIC_BOUND
            || u64::from(entry.q) > params::LARGE_PRIME_BOUND
            || !matches!(entry.sign, 1 | 2)
            || entry.multiplicity.is_zero()
        {
            return Err("invalid or unordered LP entry".to_owned());
        }
        let modulus = u128::from(entry.q);
        let trace = params::t()
            .mod_floor(&BigInt::from(entry.q))
            .to_u128()
            .ok_or_else(|| "trace residue does not fit u128".to_owned())?;
        let field_norm = params::p()
            .mod_floor(&BigInt::from(entry.q))
            .to_u128()
            .ok_or_else(|| "field-norm residue does not fit u128".to_owned())?;
        for root in [entry.smaller_root, entry.larger_root] {
            let root = u128::from(root);
            if (root * root + modulus - (trace * root) % modulus + field_norm) % modulus != 0 {
                return Err("LP endpoint is not a Frobenius-polynomial root".to_owned());
            }
        }
        previous = entry.q;
    }
    Ok(())
}

fn bounded_power(base: &BigUint, exponent: &BigUint, name: &str) -> Result<BigUint> {
    let exponent = exponent
        .to_u32()
        .filter(|value| *value <= 256)
        .ok_or_else(|| format!("{name} exponent exceeds the P-192 certificate bound"))?;
    Ok(base.pow(exponent))
}

fn product_rational_factors(factors: &[RationalFactor]) -> Result<BigUint> {
    factors.iter().try_fold(BigUint::one(), |product, factor| {
        let power = bounded_power(&factor.prime, &factor.exponent, "rational factor")?;
        Ok(product * power)
    })
}

fn lp_residual(entries: &[LargePrimeEntry]) -> Result<(BigUint, u64)> {
    let mut multiplicity_total = 0u64;
    let product = entries
        .iter()
        .try_fold(BigUint::one(), |product, entry| -> Result<BigUint> {
            let multiplicity = entry
                .multiplicity
                .to_u64()
                .filter(|value| *value <= 2)
                .ok_or_else(|| "LP multiplicity exceeds two".to_owned())?;
            multiplicity_total = multiplicity_total
                .checked_add(multiplicity)
                .ok_or_else(|| "LP multiplicity total overflow".to_owned())?;
            let power = bounded_power(&BigUint::from(entry.q), &entry.multiplicity, "large-prime")?;
            Ok(product * power)
        })?;
    Ok((product, multiplicity_total))
}

fn validate_type_and_residual(candidate_type: u8, entries: &[LargePrimeEntry]) -> Result<BigUint> {
    let (residual, multiplicity) = lp_residual(entries)?;
    if u64::from(candidate_type) != multiplicity
        || (candidate_type == 0 && (!entries.is_empty() || residual != BigUint::one()))
    {
        return Err("candidate type does not equal LP multiplicity".to_owned());
    }
    Ok(residual)
}

fn validate_candidate_semantics(certificate: &CandidateCertificate) -> Result<()> {
    if certificate.scalar_content != BigUint::one()
        || certificate.v <= BigInt::zero()
        || certificate.v.to_u64() != Some(certificate.v_box)
        || certificate.u.abs().gcd(&certificate.v) != BigInt::one()
    {
        return Err("candidate is not a primitive scalar-content-one row".to_owned());
    }
    let centered_x = &certificate.u * 2u8 + params::t() * &certificate.v;
    if centered_x.to_u64() != Some(certificate.x_box) {
        return Err("candidate centered coordinate does not match u,v".to_owned());
    }
    let computed_norm = &certificate.u * &certificate.u
        + params::t() * &certificate.u * &certificate.v
        + params::p() * &certificate.v * &certificate.v;
    if computed_norm != BigInt::from_biguint(Sign::Plus, certificate.norm.clone()) {
        return Err("candidate norm does not match alpha".to_owned());
    }
    if certificate.rational_factors.len() != certificate.small_entries.len() {
        return Err("candidate rational/small factor counts differ".to_owned());
    }
    for (factor, entry) in certificate
        .rational_factors
        .iter()
        .zip(&certificate.small_entries)
    {
        if factor.prime.to_u32() != Some(entry.ell)
            || u64::from(entry.ell) > params::ALGEBRAIC_BOUND
        {
            return Err("candidate rational/small factor primes differ".to_owned());
        }
        let magnitude = entry
            .exponent
            .abs()
            .to_biguint()
            .ok_or_else(|| "small exponent magnitude conversion failed".to_owned())?;
        if magnitude != factor.exponent
            || (entry.kind == 1
                && (entry.exponent != BigInt::one() || factor.exponent != BigUint::one()))
        {
            return Err("candidate rational/small exponents differ".to_owned());
        }
    }
    let residual = validate_type_and_residual(certificate.candidate_type, &certificate.lp_entries)?;
    if product_rational_factors(&certificate.rational_factors)? * residual != certificate.norm {
        return Err("candidate factor product does not equal its norm".to_owned());
    }
    Ok(())
}

fn validate_factorization_semantics(certificate: &FactorizationCertificate) -> Result<()> {
    let residual = validate_type_and_residual(certificate.residual_type, &certificate.lp_entries)?;
    if residual != certificate.residual
        || certificate.rational_factors.iter().any(|factor| {
            factor
                .prime
                .to_u64()
                .is_none_or(|prime| prime > params::ALGEBRAIC_BOUND)
        })
        || product_rational_factors(&certificate.rational_factors)? * &certificate.residual
            != certificate.norm
    {
        return Err("factorization certificate norm/residual/type mismatch".to_owned());
    }
    Ok(())
}

pub fn encode_candidate_certificate(certificate: &CandidateCertificate) -> Result<Vec<u8>> {
    if certificate.candidate_type > 2 {
        return Err("candidate type is outside 0..=2".to_owned());
    }
    validate_rational_factors(&certificate.rational_factors)?;
    validate_small_entries(&certificate.small_entries)?;
    validate_lp_entries(&certificate.lp_entries)?;
    validate_candidate_semantics(certificate)?;
    let mut out = Vec::new();
    out.extend_from_slice(CAND_DOMAIN);
    out.push(certificate.candidate_type);
    push_len_prefixed(&mut out, certificate.curve_uid.as_bytes())?;
    out.extend_from_slice(&certificate.factor_base_sha256);
    push_len_prefixed(&mut out, certificate.shell_id.as_bytes())?;
    out.extend_from_slice(&certificate.candidate_index.to_be_bytes());
    out.extend_from_slice(&certificate.v_box.to_be_bytes());
    out.extend_from_slice(&certificate.x_box.to_be_bytes());
    encode_sint(&certificate.u, &mut out)?;
    encode_sint(&certificate.v, &mut out)?;
    encode_nat(&certificate.norm, &mut out)?;
    encode_nat(&certificate.scalar_content, &mut out)?;
    out.extend_from_slice(
        &u32::try_from(certificate.rational_factors.len())
            .map_err(|_| "too many rational factors".to_owned())?
            .to_be_bytes(),
    );
    encode_rational_factors(&certificate.rational_factors, &mut out)?;
    out.extend_from_slice(
        &u32::try_from(certificate.small_entries.len())
            .map_err(|_| "too many small entries".to_owned())?
            .to_be_bytes(),
    );
    for entry in &certificate.small_entries {
        out.extend_from_slice(&entry.ell.to_be_bytes());
        out.push(entry.kind);
        encode_sint(&entry.exponent, &mut out)?;
    }
    out.extend_from_slice(
        &u32::try_from(certificate.lp_entries.len())
            .map_err(|_| "too many LP entries".to_owned())?
            .to_be_bytes(),
    );
    encode_lp_entries(&certificate.lp_entries, &mut out)?;
    Ok(out)
}

pub fn encode_factorization_certificate(certificate: &FactorizationCertificate) -> Result<Vec<u8>> {
    if certificate.residual_type > 2 {
        return Err("residual type is outside 0..=2".to_owned());
    }
    validate_rational_factors(&certificate.rational_factors)?;
    validate_lp_entries(&certificate.lp_entries)?;
    validate_factorization_semantics(certificate)?;
    let mut out = Vec::new();
    out.extend_from_slice(FACT_DOMAIN);
    out.extend_from_slice(&certificate.candidate_sha256);
    encode_nat(&certificate.norm, &mut out)?;
    out.extend_from_slice(
        &u32::try_from(certificate.rational_factors.len())
            .map_err(|_| "too many rational factors".to_owned())?
            .to_be_bytes(),
    );
    encode_rational_factors(&certificate.rational_factors, &mut out)?;
    encode_nat(&certificate.residual, &mut out)?;
    out.push(certificate.residual_type);
    out.extend_from_slice(
        &u32::try_from(certificate.lp_entries.len())
            .map_err(|_| "too many LP entries".to_owned())?
            .to_be_bytes(),
    );
    encode_lp_entries(&certificate.lp_entries, &mut out)?;
    Ok(out)
}

pub fn encode_candidate_record(v: u64, x: u64) -> [u8; 16] {
    let mut out = [0u8; 16];
    out[..8].copy_from_slice(&v.to_be_bytes());
    out[8..].copy_from_slice(&x.to_be_bytes());
    out
}

pub fn encode_disposition_record(
    status: u8,
    candidate_index: u64,
    certificate_hashes: Option<([u8; 32], [u8; 32])>,
) -> Result<Vec<u8>> {
    if status > 6 {
        return Err("disposition status is outside 0..=6".to_owned());
    }
    let retained = matches!(status, 1..=3);
    if retained != certificate_hashes.is_some() {
        return Err("noncanonical disposition certificate presence".to_owned());
    }
    let mut out = Vec::with_capacity(if retained { 74 } else { 10 });
    out.push(status);
    out.extend_from_slice(&candidate_index.to_be_bytes());
    out.push(u8::from(retained));
    if let Some((candidate, factorization)) = certificate_hashes {
        out.extend_from_slice(&candidate);
        out.extend_from_slice(&factorization);
    }
    Ok(out)
}

pub fn candidate_shard_digest(
    shell_id: &str,
    shard_index: u64,
    record_count: u64,
    records: &[u8],
) -> Result<[u8; 32]> {
    let expected = usize::try_from(record_count)
        .ok()
        .and_then(|n| n.checked_mul(16))
        .ok_or_else(|| "candidate record byte count overflow".to_owned())?;
    if records.len() != expected {
        return Err("candidate shard byte length mismatch".to_owned());
    }
    let mut hasher = Sha256::new();
    hasher.update(CANDIDATE_SHARD_DOMAIN);
    hasher.update(encode_shell_id(shell_id)?);
    hasher.update(shard_index.to_be_bytes());
    hasher.update(record_count.to_be_bytes());
    hasher.update(records);
    Ok(hasher.finalize().into())
}

pub fn candidate_shell_digest(shell_id: &str, shards: &[ShardDigest]) -> Result<[u8; 32]> {
    validate_shards(shards, false)?;
    let mut hasher = Sha256::new();
    hasher.update(SHELL_ROOT_DOMAIN);
    hasher.update(encode_shell_id(shell_id)?);
    hasher.update(
        u64::try_from(shards.len())
            .map_err(|_| "too many shards".to_owned())?
            .to_be_bytes(),
    );
    for shard in shards {
        hasher.update(shard.shard_index.to_be_bytes());
        hasher.update(shard.record_count.to_be_bytes());
        hasher.update(shard.sha256);
    }
    Ok(hasher.finalize().into())
}

pub fn disposition_shard_digest(
    shell_id: &str,
    shard_index: u64,
    first_candidate_index: u64,
    record_count: u64,
    records: &[u8],
) -> Result<[u8; 32]> {
    let mut hasher = Sha256::new();
    hasher.update(DISPOSITION_SHARD_DOMAIN);
    hasher.update(encode_shell_id(shell_id)?);
    hasher.update(shard_index.to_be_bytes());
    hasher.update(first_candidate_index.to_be_bytes());
    hasher.update(record_count.to_be_bytes());
    hasher.update(records);
    Ok(hasher.finalize().into())
}

pub fn disposition_shell_digest(shell_id: &str, shards: &[ShardDigest]) -> Result<[u8; 32]> {
    validate_shards(shards, true)?;
    let mut hasher = Sha256::new();
    hasher.update(DISPOSITION_SHELL_ROOT_DOMAIN);
    hasher.update(encode_shell_id(shell_id)?);
    hasher.update(
        u64::try_from(shards.len())
            .map_err(|_| "too many shards".to_owned())?
            .to_be_bytes(),
    );
    for shard in shards {
        hasher.update(shard.shard_index.to_be_bytes());
        hasher.update(shard.first_candidate_index.to_be_bytes());
        hasher.update(shard.record_count.to_be_bytes());
        hasher.update(shard.sha256);
    }
    Ok(hasher.finalize().into())
}

fn validate_shards(shards: &[ShardDigest], include_first: bool) -> Result<()> {
    let mut next_first = 0u64;
    for (expected_index, shard) in shards.iter().enumerate() {
        if shard.shard_index != expected_index as u64 {
            return Err("noncontiguous shard index".to_owned());
        }
        if include_first && shard.first_candidate_index != next_first {
            return Err("noncontiguous disposition candidate interval".to_owned());
        }
        next_first = next_first
            .checked_add(shard.record_count)
            .ok_or_else(|| "shard record count overflow".to_owned())?;
    }
    Ok(())
}

/// RFC 8785 bytes for the integer/string-only schemas used by this protocol.
/// Floating-point values are rejected rather than normalized.
pub fn canonical_json(value: &Value) -> Result<Vec<u8>> {
    fn utf16_cmp(left: &str, right: &str) -> Ordering {
        left.encode_utf16()
            .collect::<Vec<_>>()
            .cmp(&right.encode_utf16().collect::<Vec<_>>())
    }
    fn write_string(value: &str, out: &mut String) {
        out.push('"');
        for ch in value.chars() {
            match ch {
                '"' => out.push_str("\\\""),
                '\\' => out.push_str("\\\\"),
                '\u{0008}' => out.push_str("\\b"),
                '\u{0009}' => out.push_str("\\t"),
                '\u{000a}' => out.push_str("\\n"),
                '\u{000c}' => out.push_str("\\f"),
                '\u{000d}' => out.push_str("\\r"),
                ch if ch <= '\u{001f}' => {
                    out.push_str(&format!("\\u{:04x}", u32::from(ch)));
                }
                ch => out.push(ch),
            }
        }
        out.push('"');
    }
    fn write(value: &Value, out: &mut String) -> Result<()> {
        match value {
            Value::Null => out.push_str("null"),
            Value::Bool(value) => out.push_str(if *value { "true" } else { "false" }),
            Value::Number(value) if value.is_i64() || value.is_u64() => {
                out.push_str(&value.to_string());
            }
            Value::Number(_) => return Err("floating-point JSON is not admitted".to_owned()),
            Value::String(value) => write_string(value, out),
            Value::Array(values) => {
                out.push('[');
                for (index, value) in values.iter().enumerate() {
                    if index != 0 {
                        out.push(',');
                    }
                    write(value, out)?;
                }
                out.push(']');
            }
            Value::Object(values) => {
                let mut entries: Vec<_> = values.iter().collect();
                entries.sort_by(|left, right| utf16_cmp(left.0, right.0));
                out.push('{');
                for (index, (key, value)) in entries.into_iter().enumerate() {
                    if index != 0 {
                        out.push(',');
                    }
                    write_string(key, out);
                    out.push(':');
                    write(value, out)?;
                }
                out.push('}');
            }
        }
        Ok(())
    }

    let mut out = String::new();
    write(value, &mut out)?;
    Ok(out.into_bytes())
}

pub fn parse_canonical_json(bytes: &[u8]) -> Result<Value> {
    if bytes.starts_with(&[0xef, 0xbb, 0xbf]) {
        return Err("JSON has a byte-order mark".to_owned());
    }
    let value: Value = serde_json::from_slice(bytes).map_err(|error| error.to_string())?;
    if canonical_json(&value)? != bytes {
        return Err("JSON is not the exact RFC-8785 canonical form".to_owned());
    }
    Ok(value)
}
