//! Canonical JSON and the identities hashed from it.
//!
//! The repository's identity convention (`docs/curve-identities.md`,
//! `research/ic_candidate_tournament_20260915/identity.py`) hashes
//! sorted-key, compact JSON with no floats.  This is the same rule in
//! Rust, written out rather than left to `serde_json`'s map order, so an
//! identity cannot change because a dependency enabled `preserve_order`.

use serde_json::Value;

use crate::hash::sha256::sha256;

/// Lower-case hex SHA-256 of `bytes`.
pub fn sha256_hex(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

/// Canonical bytes of `v`: keys sorted, no whitespace, non-ASCII escaped.
///
/// Fails on a float: a float's decimal form is not stable across
/// producers, so it never enters an identity.  Integers and strings do.
pub fn canonical(v: &Value) -> Result<String, String> {
    let mut out = String::new();
    write(v, &mut out)?;
    Ok(out)
}

fn write(v: &Value, out: &mut String) -> Result<(), String> {
    match v {
        Value::Null => out.push_str("null"),
        Value::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Value::Number(n) => {
            if n.is_f64() {
                return Err(format!("float {n} in an identity; write it as a string"));
            }
            out.push_str(&n.to_string());
        }
        Value::String(s) => write_str(s, out),
        Value::Array(items) => {
            out.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write(item, out)?;
            }
            out.push(']');
        }
        Value::Object(map) => {
            let mut keys: Vec<&String> = map.keys().collect();
            keys.sort();
            out.push('{');
            for (i, k) in keys.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write_str(k, out);
                out.push(':');
                write(&map[*k], out)?;
            }
            out.push('}');
        }
    }
    Ok(())
}

/// JSON string escaping as Python's `json.dumps(..., ensure_ascii=True)`
/// writes it, so a Python reader hashing the same object agrees.
fn write_str(s: &str, out: &mut String) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c if (c as u32) < 0x80 => out.push(c),
            c => {
                let mut buf = [0u16; 2];
                for unit in c.encode_utf16(&mut buf) {
                    out.push_str(&format!("\\u{:04x}", unit));
                }
            }
        }
    }
    out.push('"');
}

/// Full SHA-256 of the canonical form.
pub fn digest(v: &Value) -> Result<String, String> {
    Ok(sha256_hex(canonical(v)?.as_bytes()))
}

/// `prefix` + `h` + the first 12 hex digits of the canonical digest, and
/// the full digest beside it.  The short form names; the full form is
/// what a join or an audit checks.
pub fn short_id(prefix: &str, v: &Value) -> Result<(String, String), String> {
    let full = digest(v)?;
    Ok((format!("{prefix}h{}", &full[..12]), full))
}

/// `prefix` + the first 12 hex digits, with no `h`: the workload form,
/// so a run id composes as `<method>W<12 hex>R<n>`, the repository's
/// `{candidate}W{workload}R{number}` convention.
pub fn bare_id(prefix: &str, v: &Value) -> Result<(String, String), String> {
    let full = digest(v)?;
    Ok((format!("{prefix}{}", &full[..12]), full))
}

/// SplitMix64: the one mixing function every derived seed goes through,
/// so a seed is reproducible from its inputs in any language.
pub fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// A 64-bit value derived from a domain label and integers, by SHA-256 of
/// their canonical JSON.  Used for planted logarithms and algorithm
/// seeds, so neither depends on an RNG crate's stream.
pub fn derive_u64(domain: &str, parts: &[u64]) -> u64 {
    let d = derive_bytes(domain, parts);
    u64::from_be_bytes(d[..8].try_into().expect("eight bytes"))
}

/// The whole digest [`derive_u64`] takes its first eight bytes from, for
/// a derivation that needs more than a word (a scalar on a two-word
/// subgroup takes sixteen).
pub fn derive_bytes(domain: &str, parts: &[u64]) -> [u8; 32] {
    let v = serde_json::json!({"domain": domain, "parts": parts});
    sha256(canonical(&v).expect("integers only").as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn canonical_sorts_and_compacts() {
        let v = json!({"b": 1, "a": [true, null, "x"], "c": {"z": 0, "y": "é"}});
        assert_eq!(
            canonical(&v).unwrap(),
            "{\"a\":[true,null,\"x\"],\"b\":1,\"c\":{\"y\":\"\\u00e9\",\"z\":0}}"
        );
    }

    #[test]
    fn floats_are_refused() {
        assert!(canonical(&json!({"x": 1.5})).is_err());
    }

    #[test]
    fn short_id_is_a_prefix_of_the_digest() {
        let (id, full) = short_id("W", &json!({"k": 1})).unwrap();
        assert_eq!(&id[2..], &full[..12]);
        assert_eq!(full.len(), 64);
    }
}
