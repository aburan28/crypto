//! Records a walk writes: a small value tree with exact integers, its YAML
//! and JSON forms, the canonical JSON the identities hash, and the EC1
//! identity of a prime-field representation.
//!
//! The YAML follows `docs/curves/ic/curves.yaml` (schema
//! `docs/curves/ic/curves.schema.json`); integers of any size are written
//! bare, as that file writes them.  JSON output writes field elements and
//! other wide integers as decimal strings so a reader without arbitrary
//! precision parses it exactly.

use num_bigint::{BigInt, BigUint};

use crate::hash::sha256::sha256;

/// A record value.  `Int` holds the decimal digits of an integer of any
/// size; `Map` keeps insertion order for YAML and is sorted for canonical
/// JSON.
#[derive(Clone, Debug, PartialEq)]
pub enum V {
    Null,
    Bool(bool),
    Int(String),
    Str(String),
    Seq(Vec<V>),
    Map(Vec<(String, V)>),
}

impl V {
    pub fn int<T: std::fmt::Display>(v: T) -> V {
        V::Int(v.to_string())
    }
    pub fn s(v: impl Into<String>) -> V {
        V::Str(v.into())
    }
    pub fn big(v: &BigUint) -> V {
        V::Int(v.to_string())
    }
    pub fn bigi(v: &BigInt) -> V {
        V::Int(v.to_string())
    }
    pub fn map(pairs: Vec<(&str, V)>) -> V {
        V::Map(pairs.into_iter().map(|(k, v)| (k.to_string(), v)).collect())
    }
    pub fn opt(v: Option<V>) -> V {
        v.unwrap_or(V::Null)
    }

    /// Canonical JSON: sorted keys, no whitespace, integers bare.  The
    /// identities in this repository hash exactly this form
    /// (`json.dumps(sort_keys=True, separators=(",", ":"))`).
    pub fn canonical_json(&self) -> String {
        let mut out = String::new();
        self.json_into(&mut out, true, false);
        out
    }

    /// Pretty JSON in insertion order; integers wider than fifteen digits
    /// are written as decimal strings.
    pub fn json(&self) -> String {
        let mut out = String::new();
        self.json_pretty(&mut out, 0);
        out.push('\n');
        out
    }

    fn json_into(&self, out: &mut String, sort: bool, ints_as_strings: bool) {
        match self {
            V::Null => out.push_str("null"),
            V::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            V::Int(d) => {
                if ints_as_strings {
                    out.push_str(&json_string(d))
                } else {
                    out.push_str(d)
                }
            }
            V::Str(s) => out.push_str(&json_string(s)),
            V::Seq(items) => {
                out.push('[');
                for (i, it) in items.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    it.json_into(out, sort, ints_as_strings);
                }
                out.push(']');
            }
            V::Map(pairs) => {
                let mut refs: Vec<&(String, V)> = pairs.iter().collect();
                if sort {
                    refs.sort_by(|a, b| a.0.cmp(&b.0));
                }
                out.push('{');
                for (i, (k, v)) in refs.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    out.push_str(&json_string(k));
                    out.push(':');
                    v.json_into(out, sort, ints_as_strings);
                }
                out.push('}');
            }
        }
    }

    fn json_pretty(&self, out: &mut String, indent: usize) {
        let pad = |n: usize| " ".repeat(n);
        match self {
            V::Int(d) if d.trim_start_matches('-').len() > 15 => out.push_str(&json_string(d)),
            V::Seq(items) if !items.is_empty() => {
                if items.iter().all(|i| !matches!(i, V::Seq(_) | V::Map(_))) {
                    // Wide integers become strings even inline.
                    let mut s2 = String::from("[");
                    for (i, it) in items.iter().enumerate() {
                        if i > 0 {
                            s2.push_str(", ");
                        }
                        let mut one = String::new();
                        it.json_pretty(&mut one, 0);
                        s2.push_str(&one);
                    }
                    s2.push(']');
                    out.push_str(&s2);
                    return;
                }
                out.push_str("[\n");
                for (i, it) in items.iter().enumerate() {
                    out.push_str(&pad(indent + 2));
                    it.json_pretty(out, indent + 2);
                    if i + 1 < items.len() {
                        out.push(',');
                    }
                    out.push('\n');
                }
                out.push_str(&pad(indent));
                out.push(']');
            }
            V::Map(pairs) if !pairs.is_empty() => {
                out.push_str("{\n");
                for (i, (k, v)) in pairs.iter().enumerate() {
                    out.push_str(&pad(indent + 2));
                    out.push_str(&json_string(k));
                    out.push_str(": ");
                    v.json_pretty(out, indent + 2);
                    if i + 1 < pairs.len() {
                        out.push(',');
                    }
                    out.push('\n');
                }
                out.push_str(&pad(indent));
                out.push('}');
            }
            other => other.json_into(out, false, false),
        }
    }

    /// Block-style YAML, two-space indentation, as `curves.yaml` is written.
    pub fn yaml(&self) -> String {
        let mut out = String::new();
        match self {
            V::Map(pairs) if !pairs.is_empty() => yaml_map(&mut out, pairs, 0),
            other => {
                out.push_str(&yaml_scalar(other));
                out.push('\n');
            }
        }
        out
    }
}

fn json_string(s: &str) -> String {
    let mut out = String::from("\"");
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn plain_ok(s: &str) -> bool {
    const RESERVED: [&str; 11] = [
        "true", "false", "null", "yes", "no", "on", "off", "y", "n", "~", "",
    ];
    let first = s.chars().next();
    first.is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || "_./+-=^*()".contains(c))
        && !RESERVED.contains(&s.to_ascii_lowercase().as_str())
}

fn yaml_key(k: &str) -> String {
    if plain_ok(k) {
        k.to_string()
    } else {
        format!("'{}'", k.replace('\'', "''"))
    }
}

fn yaml_scalar(v: &V) -> String {
    match v {
        V::Null => "null".into(),
        V::Bool(b) => b.to_string(),
        V::Int(d) => d.clone(),
        V::Str(s) => {
            if plain_ok(s) {
                s.clone()
            } else {
                json_string(s)
            }
        }
        V::Seq(_) => "[]".into(),
        V::Map(_) => "{}".into(),
    }
}

fn yaml_map(out: &mut String, pairs: &[(String, V)], indent: usize) {
    for (k, v) in pairs {
        out.push_str(&" ".repeat(indent));
        out.push_str(&yaml_key(k));
        out.push(':');
        match v {
            V::Map(m) if !m.is_empty() => {
                out.push('\n');
                yaml_map(out, m, indent + 2);
            }
            V::Seq(s) if !s.is_empty() => {
                out.push('\n');
                yaml_seq(out, s, indent + 2);
            }
            other => {
                out.push(' ');
                out.push_str(&yaml_scalar(other));
                out.push('\n');
            }
        }
    }
}

fn yaml_seq(out: &mut String, items: &[V], indent: usize) {
    for it in items {
        out.push_str(&" ".repeat(indent));
        out.push('-');
        match it {
            V::Map(m) if !m.is_empty() => {
                // First pair on the dash line, the rest aligned under it.
                let mut inner = String::new();
                yaml_map(&mut inner, m, indent + 2);
                out.push(' ');
                out.push_str(inner.trim_start());
            }
            V::Seq(s) if !s.is_empty() => {
                out.push('\n');
                yaml_seq(out, s, indent + 2);
            }
            other => {
                out.push(' ');
                out.push_str(&yaml_scalar(other));
                out.push('\n');
            }
        }
    }
}

pub fn sha256_hex(s: &str) -> String {
    sha256(s.as_bytes())
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// The EC1 field and curve records of a prime-field representation, in
/// the encoding `scripts/build_curve_registry.py::prime_record` writes.
pub fn ec1_records(
    p: &BigUint,
    a: &BigUint,
    b: &BigUint,
    subgroup_order: &BigUint,
    cofactor: &BigUint,
    g: (&BigUint, &BigUint),
) -> (V, V) {
    let field = V::map(vec![
        ("characteristic", V::big(p)),
        ("degree", V::int(1)),
        ("representation", V::s("prime")),
        (
            "element_encoding",
            V::s("hex integer modulo the characteristic"),
        ),
    ]);
    let curve = V::map(vec![
        ("model", V::s("short Weierstrass")),
        (
            "coefficients",
            V::Seq(vec![V::int(0), V::int(0), V::int(0), V::big(a), V::big(b)]),
        ),
        ("subgroup_order", V::s(subgroup_order.to_string())),
        ("cofactor", V::big(cofactor)),
        (
            "generator",
            V::Seq(vec![
                V::s(format!("0x{:x}", g.0)),
                V::s(format!("0x{:x}", g.1)),
            ]),
        ),
        ("target_group", V::s("prime-order subgroup")),
    ]);
    (field, curve)
}

/// `(EC1 alias, curve UID)` for a field and curve record, as
/// `tools/curve_identity.py::curve_identity` computes them.
pub fn ec1_identity(p: &BigUint, field: &V, curve: &V, tag: &str) -> (String, String) {
    let digest = sha256_hex(
        &V::map(vec![("field", field.clone()), ("curve", curve.clone())]).canonical_json(),
    );
    (
        format!("EC1P{}C{tag}h{}", p.bits(), &digest[..12]),
        format!("urn:ec-record:1:sha256:{digest}"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecc::curve::CurveParams;

    #[test]
    fn p256_ec1_matches_the_registry() {
        let c = CurveParams::p256();
        let (field, curve) =
            ec1_records(&c.p, &c.a, &c.b, &c.n, &BigUint::from(1u8), (&c.gx, &c.gy));
        let (id, uid) = ec1_identity(&c.p, &field, &curve, "p256");
        assert_eq!(id, "EC1P256Cp256h0523b774e066");
        assert_eq!(
            uid,
            "urn:ec-record:1:sha256:0523b774e0666cf37ff4dc9cf452eae024192a7f89bdad10ca226d6e64d4dd70"
        );
    }

    #[test]
    fn yaml_quotes_what_would_misparse() {
        let v = V::map(vec![
            ("j", V::s("0x1")),
            ("263", V::map(vec![("level", V::int(0))])),
            ("flag", V::s("true")),
            ("plain", V::s("short_weierstrass")),
            ("empty", V::Seq(vec![])),
        ]);
        let y = v.yaml();
        assert!(y.contains("j: \"0x1\"\n"));
        assert!(y.contains("'263':\n  level: 0\n"));
        assert!(y.contains("flag: \"true\"\n"));
        assert!(y.contains("plain: short_weierstrass\n"));
        assert!(y.contains("empty: []\n"));
    }
}
