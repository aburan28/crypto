use std::{fs, path::Path};

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{Signed, ToPrimitive};
use serde::{Deserialize, Serialize};

use crate::{
    encoding::{parse_canonical_json, sha256},
    params::{self, ALGEBRAIC_BOUND, CURVE_UID, D_DEC, MAPPABLE_BOUND, P_DEC, T_DEC},
    Result,
};

pub const FACTOR_BASE_SCHEMA: &str = "p192-wcm-factor-base-v1";

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FactorBaseEntry {
    pub ell: u64,
    /// Zero split, one ramified.
    pub kind: u8,
    pub kronecker: i8,
    pub roots: Vec<u64>,
    pub positive_root_index: u8,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FactorBaseManifest {
    pub schema: String,
    pub curve_uid: String,
    pub p: String,
    pub t: String,
    #[serde(rename = "D")]
    pub discriminant: String,
    pub algebraic_bound: u64,
    pub mappable_bound: u64,
    pub source_commit: String,
    pub entry_count: u64,
    pub entries: Vec<FactorBaseEntry>,
}

#[derive(Clone, Debug)]
pub struct VerifiedFactorBase {
    pub manifest: FactorBaseManifest,
    pub canonical_bytes: Vec<u8>,
    pub sha256: [u8; 32],
}

fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((u128::from(left) * u128::from(right)) % u128::from(modulus)) as u64
}

fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut accumulator = 1 % modulus;
    base %= modulus;
    while exponent != 0 {
        if exponent & 1 == 1 {
            accumulator = mul_mod(accumulator, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exponent >>= 1;
    }
    accumulator
}

fn mod_u64(value: &BigInt, modulus: u64) -> u64 {
    value
        .mod_floor(&BigInt::from(modulus))
        .to_u64()
        .expect("residue is below u64 modulus")
}

fn sqrt_mod_prime(value: u64, prime: u64) -> Option<u64> {
    debug_assert!(prime > 2 && prime % 2 == 1);
    let value = value % prime;
    if value == 0 {
        return Some(0);
    }
    if pow_mod(value, (prime - 1) / 2, prime) != 1 {
        return None;
    }
    if prime % 4 == 3 {
        let root = pow_mod(value, (prime + 1) / 4, prime);
        return Some(root.min(prime - root));
    }

    let mut odd = prime - 1;
    let mut power = 0u32;
    while odd.is_multiple_of(2) {
        odd /= 2;
        power += 1;
    }
    let non_residue =
        (2..prime).find(|candidate| pow_mod(*candidate, (prime - 1) / 2, prime) == prime - 1)?;
    let mut c = pow_mod(non_residue, odd, prime);
    let mut x = pow_mod(value, odd.div_ceil(2), prime);
    let mut t = pow_mod(value, odd, prime);
    let mut m = power;
    while t != 1 {
        let mut i = 1u32;
        let mut probe = mul_mod(t, t, prime);
        while probe != 1 {
            probe = mul_mod(probe, probe, prime);
            i += 1;
            if i == m {
                return None;
            }
        }
        let b = pow_mod(c, 1u64 << (m - i - 1), prime);
        x = mul_mod(x, b, prime);
        let b_squared = mul_mod(b, b, prime);
        t = mul_mod(t, b_squared, prime);
        c = b_squared;
        m = i;
    }
    Some(x.min(prime - x))
}

pub fn primes_through(limit: u64) -> Result<Vec<u64>> {
    let len = usize::try_from(limit)
        .map_err(|_| "prime bound does not fit usize".to_owned())?
        .checked_add(1)
        .ok_or_else(|| "prime bound length overflow".to_owned())?;
    let mut flags = vec![true; len];
    if len > 0 {
        flags[0] = false;
    }
    if len > 1 {
        flags[1] = false;
    }
    let mut prime = 2usize;
    while prime <= (len.saturating_sub(1)) / prime {
        if flags[prime] {
            let mut multiple = prime * prime;
            while multiple < len {
                flags[multiple] = false;
                multiple += prime;
            }
        }
        prime += 1;
    }
    Ok(flags
        .into_iter()
        .enumerate()
        .filter_map(|(value, is_prime)| is_prime.then_some(value as u64))
        .collect())
}

pub fn kronecker_at_prime(discriminant: &BigInt, prime: u64) -> i8 {
    if prime == 2 {
        if discriminant.is_even() {
            return 0;
        }
        return match mod_u64(discriminant, 8) {
            1 | 7 => 1,
            3 | 5 => -1,
            _ => unreachable!("odd residue modulo eight"),
        };
    }
    let residue = mod_u64(discriminant, prime);
    if residue == 0 {
        0
    } else if pow_mod(residue, (prime - 1) / 2, prime) == 1 {
        1
    } else {
        -1
    }
}

pub fn roots_for_prime(
    prime: u64,
    trace: &BigInt,
    field_norm: &BigInt,
    discriminant: &BigInt,
) -> Result<Vec<u64>> {
    if prime == 2 {
        let mut roots = Vec::new();
        for root in 0..2 {
            let value =
                (root * root + 2 - (mod_u64(trace, 2) * root) % 2 + mod_u64(field_norm, 2)) % 2;
            if value == 0 {
                roots.push(root);
            }
        }
        return Ok(roots);
    }
    let symbol = kronecker_at_prime(discriminant, prime);
    if symbol < 0 {
        return Ok(Vec::new());
    }
    let sqrt_d = sqrt_mod_prime(mod_u64(discriminant, prime), prime)
        .ok_or_else(|| format!("Kronecker/root disagreement at prime {prime}"))?;
    let trace = mod_u64(trace, prime);
    let inverse_two = prime.div_ceil(2);
    let mut roots = vec![
        mul_mod((trace + sqrt_d) % prime, inverse_two, prime),
        mul_mod((trace + prime - sqrt_d) % prime, inverse_two, prime),
    ];
    roots.sort_unstable();
    roots.dedup();
    let field_norm = mod_u64(field_norm, prime);
    if roots.iter().any(|root| {
        let modulus = u128::from(prime);
        let root = u128::from(*root);
        (root * root + modulus - (u128::from(trace) * root) % modulus + u128::from(field_norm))
            % modulus
            != 0
    }) {
        return Err(format!("invalid characteristic root at prime {prime}"));
    }
    Ok(roots)
}

pub fn regenerate_entries() -> Result<Vec<FactorBaseEntry>> {
    let discriminant = params::d();
    let trace = params::t();
    let field_norm = params::p();
    let mut entries = Vec::new();
    for ell in primes_through(ALGEBRAIC_BOUND)? {
        let kronecker = kronecker_at_prime(&discriminant, ell);
        if kronecker < 0 {
            continue;
        }
        let roots = roots_for_prime(ell, &trace, &field_norm, &discriminant)?;
        let expected_len = if kronecker == 0 { 1 } else { 2 };
        if roots.len() != expected_len {
            return Err(format!(
                "expected {expected_len} root(s) at ell={ell}, found {}",
                roots.len()
            ));
        }
        entries.push(FactorBaseEntry {
            ell,
            kind: u8::from(kronecker == 0),
            kronecker,
            roots,
            positive_root_index: 0,
        });
    }
    Ok(entries)
}

fn minimal_decimal(text: &str, signed: bool) -> Result<()> {
    if text.is_empty() || text.starts_with('+') || text == "-0" {
        return Err("noncanonical decimal string".to_owned());
    }
    if !signed && text.starts_with('-') {
        return Err("negative unsigned decimal string".to_owned());
    }
    let parsed = if signed {
        BigInt::parse_bytes(text.as_bytes(), 10)
            .ok_or_else(|| "invalid signed decimal string".to_owned())?
            .to_string()
    } else {
        BigUint::parse_bytes(text.as_bytes(), 10)
            .ok_or_else(|| "invalid unsigned decimal string".to_owned())?
            .to_string()
    };
    if parsed != text {
        return Err("nonminimal decimal string".to_owned());
    }
    Ok(())
}

fn lowercase_hex(text: &str, len: usize) -> bool {
    text.len() == len
        && text
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

pub fn verify_manifest(manifest: &FactorBaseManifest) -> Result<()> {
    if manifest.schema != FACTOR_BASE_SCHEMA
        || manifest.curve_uid != CURVE_UID
        || manifest.p != P_DEC
        || manifest.t != T_DEC
        || manifest.discriminant != D_DEC
        || manifest.algebraic_bound != ALGEBRAIC_BOUND
        || manifest.mappable_bound != MAPPABLE_BOUND
    {
        return Err("factor-base envelope differs from the frozen P-192 contract".to_owned());
    }
    minimal_decimal(&manifest.p, false)?;
    minimal_decimal(&manifest.t, false)?;
    minimal_decimal(&manifest.discriminant, true)?;
    if !lowercase_hex(&manifest.source_commit, 40) {
        return Err("factor-base source_commit is not lowercase 40-hex".to_owned());
    }
    if manifest.entry_count != manifest.entries.len() as u64 {
        return Err("factor-base entry_count mismatch".to_owned());
    }
    let p = params::p();
    let t = params::t();
    let d = params::d();
    if &t * &t - &p * 4 != d || !d.is_negative() {
        return Err("frozen D=t^2-4p identity failed".to_owned());
    }
    let expected = regenerate_entries()?;
    if manifest.entries != expected {
        return Err("factor-base entries are not all-and-only the regenerated base".to_owned());
    }
    Ok(())
}

pub fn verify_factor_base(source_dir: &Path) -> Result<VerifiedFactorBase> {
    let json_path = source_dir.join("factor-base.json");
    let sidecar_path = source_dir.join("factor-base.sha256");
    let canonical_bytes =
        fs::read(&json_path).map_err(|error| format!("read {}: {error}", json_path.display()))?;
    let sidecar = fs::read(&sidecar_path)
        .map_err(|error| format!("read {}: {error}", sidecar_path.display()))?;
    verify_factor_base_bytes(&canonical_bytes, &sidecar)
}

pub fn verify_factor_base_bytes(
    canonical_bytes: &[u8],
    sidecar: &[u8],
) -> Result<VerifiedFactorBase> {
    let value = parse_canonical_json(canonical_bytes)?;
    let manifest: FactorBaseManifest =
        serde_json::from_value(value).map_err(|error| format!("factor-base schema: {error}"))?;
    verify_manifest(&manifest)?;
    let digest = sha256(canonical_bytes);
    let expected_sidecar = format!("{}  factor-base.json\n", hex::encode(digest));
    if sidecar != expected_sidecar.as_bytes() {
        return Err("factor-base.sha256 does not bind exact factor-base.json bytes".to_owned());
    }
    Ok(VerifiedFactorBase {
        manifest,
        canonical_bytes: canonical_bytes.to_vec(),
        sha256: digest,
    })
}
