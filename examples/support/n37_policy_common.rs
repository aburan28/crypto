//! Frozen inputs shared by the n37 four-policy producer and independent replay.
//! Neither caller is given a target scalar.

use crypto_lib::binary_ecc::{F2mElement, F2mPoly, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_field_basis::BinaryFieldBasis;
use crypto_lib::cryptanalysis::binary_velu::{
    elt, isogeny_from_kernel, to_u64, BinaryIsogeny, Curve,
};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::hash::sha256::sha256;
use num_traits::ToPrimitive;
use serde_json::Value;
use std::fs;
use std::io::Read;

pub const SOURCE_PATH: &str =
    "research/notes/ecc2k130/n37_four_policy_support_20261003/SOURCE42.jsonl";
pub const SOURCE_SHA: &str = "0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75";
pub const SOURCE_BASE_HASH: &str =
    "8423b135df3515b0e284126d9eeb95f71e3901ae12908b3cad8eb03d1f31d4bc";
pub const NATIVE_PATH: &str =
    "research/notes/ecc2k130/n37_native_basis_bridge_20261002/NATIVE42.json";
pub const NATIVE_SHA: &str = "bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c";
pub const ARCHIVE_PATH: &str = "research/koblitz_isogeny_descent_37_results_20260925/raw.json.gz";
pub const ARCHIVE_SHA: &str = "eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90";
pub const ORDER: u64 = 137_439_487_532;
pub const R: u64 = 230_603_167;
pub const K: usize = 42;
pub const N: usize = 37;

pub struct Context {
    pub kc: KoblitzCurve,
    pub control: Curve,
    pub archived_source: Curve,
    pub leaf: Curve,
    pub bridge: BinaryFieldBasis,
    pub isogeny: BinaryIsogeny,
    pub source: Value,
    pub native: Value,
    pub lambda: u64,
}

pub fn verified_bytes(path: &str, expected: &str) -> Result<Vec<u8>, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {path}: {error}"))?;
    let actual = hex::encode(sha256(&bytes));
    if actual != expected {
        return Err(format!("SHA-256 mismatch for {path}: {actual}"));
    }
    Ok(bytes)
}

pub fn point_from_value(value: &Value) -> Result<FastPoint, String> {
    let words: [u64; 2] = serde_json::from_value(value.clone())
        .map_err(|error| format!("point coordinates: {error}"))?;
    Ok(FastPoint::affine(words[0], words[1]))
}

pub fn point_words(point: FastPoint) -> Result<[u64; 2], String> {
    if point.infinity {
        return Err("support factor mapped to infinity".into());
    }
    Ok([point.x, point.y])
}

fn archive_field_element(expression: &str) -> Result<u64, String> {
    let mut value = 0u64;
    for monomial in expression.trim_matches(['(', ')']).split(" + ") {
        let bit = match monomial {
            "1" => 0,
            "a" => 1,
            power => power
                .strip_prefix("a^")
                .ok_or_else(|| format!("invalid archived field monomial {power}"))?
                .parse::<u32>()
                .map_err(|error| format!("archived field exponent: {error}"))?,
        };
        if bit >= 37 {
            return Err("archived field coefficient exceeds degree 37".into());
        }
        value ^= 1u64 << bit;
    }
    Ok(value)
}

fn archive_x_degree(suffix: &str) -> Result<usize, String> {
    if suffix.is_empty() {
        Ok(1)
    } else {
        suffix
            .strip_prefix('^')
            .ok_or("invalid archived x power".to_string())?
            .parse()
            .map_err(|error| format!("archived x exponent: {error}"))
    }
}

fn archive_kernel(expression: &str) -> Result<F2mPoly, String> {
    let mut terms = Vec::new();
    let mut start = 0;
    let mut depth = 0i32;
    for (index, symbol) in expression.char_indices() {
        match symbol {
            '(' => depth += 1,
            ')' => depth -= 1,
            '+' if depth == 0 => {
                terms.push(expression[start..index].trim());
                start = index + 1;
            }
            _ => {}
        }
        if depth < 0 {
            return Err("unbalanced archived kernel".into());
        }
    }
    if depth != 0 {
        return Err("unbalanced archived kernel".into());
    }
    terms.push(expression[start..].trim());
    let mut coeffs = vec![F2mElement::zero(37); 37];
    for term in terms {
        let (coefficient, degree) = if let Some((coefficient, suffix)) = term.rsplit_once("*x") {
            (
                archive_field_element(coefficient)?,
                archive_x_degree(suffix)?,
            )
        } else if let Some(suffix) = term.strip_prefix('x') {
            (1, archive_x_degree(suffix)?)
        } else {
            (archive_field_element(term)?, 0)
        };
        if degree > 36 {
            return Err("archived kernel exceeds degree 36".into());
        }
        coeffs[degree] = elt(to_u64(&coeffs[degree]) ^ coefficient, 37);
    }
    let kernel = F2mPoly::from_coeffs(coeffs, 37);
    if kernel.degree() != Some(36) || to_u64(&kernel.lead()) != 1 {
        return Err("archived kernel is not monic degree 36".into());
    }
    Ok(kernel)
}

impl Context {
    pub fn load() -> Result<Self, String> {
        let source_bytes = verified_bytes(SOURCE_PATH, SOURCE_SHA)?;
        let source_line = source_bytes
            .split(|&byte| byte == b'\n')
            .next()
            .ok_or("empty source base header")?;
        let source: Value = serde_json::from_slice(source_line)
            .map_err(|error| format!("source base header: {error}"))?;
        if source["kind"] != "point_defined_factor_base"
            || source["n"].as_u64() != Some(37)
            || source["a"].as_u64() != Some(0)
            || source["subgroup_order"].as_u64() != Some(R)
            || source["orbit_columns"].as_u64() != Some(K as u64)
            || source["factor_base_points"].as_u64() != Some((2 * N * K) as u64)
            || source["base_hash"].as_str() != Some(SOURCE_BASE_HASH)
        {
            return Err("frozen source base metadata changed".into());
        }
        let native: Value = serde_json::from_slice(&verified_bytes(NATIVE_PATH, NATIVE_SHA)?)
            .map_err(|error| format!("native manifest: {error}"))?;
        if native["rows"].as_array().map(Vec::len) != Some(K)
            || native["group_order"].as_u64() != Some(ORDER)
            || native["subgroup_order"].as_u64() != Some(R)
            || native["source_generator_image"].as_u64() != Some(10_156_182_909)
        {
            return Err("frozen native manifest metadata changed".into());
        }
        let kc = KoblitzCurve::new(0, 37).ok_or("source Koblitz curve")?;
        let lambda = kc.lambda.to_u64().ok_or("lambda does not fit u64")?;
        if kc.subgroup_order.to_u64() != Some(R)
            || kc.curve.irreducible.low_terms != vec![0, 1, 4, 6]
        {
            return Err("registered source curve changed".into());
        }
        let control =
            Curve::with_order(37, &kc.curve.irreducible, 0, 1, ORDER).ok_or("control curve")?;
        let archive_modulus = IrreduciblePoly {
            degree: 37,
            low_terms: vec![0, 1, 2, 3, 4, 5],
        };
        let archived_source =
            Curve::with_order(37, &archive_modulus, 0, 1, ORDER).ok_or("archive source curve")?;
        let bridge = BinaryFieldBasis::from_generator_image(
            &kc.curve.irreducible,
            &archive_modulus,
            10_156_182_909,
        )?;
        let archive_bytes = verified_bytes(ARCHIVE_PATH, ARCHIVE_SHA)?;
        let mut decompressed = String::new();
        flate2::read::GzDecoder::new(&archive_bytes[..])
            .read_to_string(&mut decompressed)
            .map_err(|error| format!("archive gzip: {error}"))?;
        let archive: Value = serde_json::from_str(&decompressed)
            .map_err(|error| format!("archive JSON: {error}"))?;
        let certificate = &archive["isogeny_certificate"];
        if certificate["degree"].as_u64() != Some(73)
            || certificate["field_modulus"] != "x^37 + x^5 + x^4 + x^3 + x^2 + x + 1"
        {
            return Err("archived degree-73 certificate changed".into());
        }
        let kernel = archive_kernel(
            certificate["kernel_polynomial"]
                .as_str()
                .ok_or("archived kernel polynomial")?,
        )?;
        let isogeny = isogeny_from_kernel(&archived_source, kernel, 73);
        let leaf_a6 = archive_field_element(
            certificate["target_normalized_coefficients"]["a6"]
                .as_str()
                .ok_or("archived leaf a6")?,
        )?;
        if isogeny.a6_codomain != leaf_a6 {
            return Err("reconstructed leaf coefficient mismatch".into());
        }
        let leaf =
            Curve::with_order(37, &archive_modulus, 0, leaf_a6, ORDER).ok_or("leaf curve")?;
        Ok(Self {
            kc,
            control,
            archived_source,
            leaf,
            bridge,
            isogeny,
            source,
            native,
            lambda,
        })
    }

    pub fn source_seeds(&self) -> Result<Vec<FastPoint>, String> {
        let rows = self.source["factor_base_representatives"]
            .as_array()
            .ok_or("source base representatives")?;
        if rows.len() != K {
            return Err("source base does not have 42 representatives".into());
        }
        rows.iter().map(point_from_value).collect()
    }

    pub fn native_seeds(&self) -> Result<Vec<FastPoint>, String> {
        self.native["rows"]
            .as_array()
            .ok_or("native rows")?
            .iter()
            .map(|row| point_from_value(&row["leaf"]))
            .collect()
    }

    // Replay deliberately recomputes this independently instead of calling it.
    #[allow(dead_code)]
    pub fn lambda_powers(&self) -> Vec<u64> {
        let mut power = 1u64;
        let mut powers = Vec::with_capacity(N);
        for _ in 0..N {
            powers.push(power);
            power = ((power as u128 * self.lambda as u128) % R as u128) as u64;
        }
        powers
    }
}
