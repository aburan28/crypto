//! Curves from outside the registry: isogeny-walk route files and the
//! IC crosswalk (`docs/curves/ic/curves.yaml`).
//!
//! The manifest `docs/curves/trait-sources.json` names each file and how
//! to read it.  A walked curve's ICV1 slug is recomputed from its model
//! and order (`curve_id::prime`) and must equal the one the walker
//! recorded; such slugs are computed identities, not registrations.  A
//! curve whose slug is already in the registry is skipped: the registry
//! record wins.  The IC crosswalk entries that point at registry models
//! are skipped for the same reason; the others are read in full.

use std::collections::BTreeMap;

use num_bigint::{BigInt, BigUint};
use num_traits::Zero;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use super::{irreducible, GeneratorCoords, Model, RegistryCurve, Representation, Source};
use crate::binary_ecc::F2mElement;
use crate::cryptanalysis::curve_id;

/// One entry of `docs/curves/trait-sources.json`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SourceSpec {
    /// `isogeny-walk` or `ic-crosswalk`.
    pub kind: String,
    /// Repository-relative path.
    pub path: String,
    /// The family label given to the curves it adds.
    pub family: String,
}

#[derive(Deserialize)]
struct Manifest {
    sources: Vec<SourceSpec>,
}

pub fn read_manifest(text: &str) -> Result<Vec<SourceSpec>, String> {
    let m: Manifest = serde_json::from_str(text).map_err(|e| format!("trait sources: {e}"))?;
    for s in &m.sources {
        if s.kind != "isogeny-walk" && s.kind != "ic-crosswalk" {
            return Err(format!("{}: unknown source kind {}", s.path, s.kind));
        }
    }
    Ok(m.sources)
}

fn uint(v: &Value, what: &str) -> Result<BigUint, String> {
    match v {
        Value::String(s) => s.parse().ok(),
        Value::Number(n) => n.as_u64().map(BigUint::from),
        _ => None,
    }
    .ok_or_else(|| format!("{what}: not an unsigned integer"))
}

fn int(v: &Value, what: &str) -> Result<BigInt, String> {
    match v {
        Value::String(s) => s.parse().ok(),
        Value::Number(n) => n.as_i64().map(BigInt::from),
        _ => None,
    }
    .ok_or_else(|| format!("{what}: not an integer"))
}

fn source(spec: &SourceSpec, uid: Option<&str>) -> Source {
    Source {
        kind: spec.kind.clone(),
        path: Some(spec.path.clone()),
        curve_uid: uid.map(str::to_string),
        registered: false,
    }
}

/// Levels written as `{"ℓ": level}` or `{"ℓ": {"level": l, "status": "proved"}}`;
/// only proved ones are kept.
fn levels(v: &Value, what: &str) -> Result<BTreeMap<u64, u32>, String> {
    let mut out = BTreeMap::new();
    let Some(map) = v.as_object() else {
        return Ok(out);
    };
    for (ell, l) in map {
        let ell: u64 = ell
            .parse()
            .map_err(|_| format!("{what}: volcano prime {ell}"))?;
        let level = match l {
            Value::Object(o) if o.get("status").and_then(Value::as_str) == Some("proved") => {
                o.get("level")
            }
            Value::Object(_) => continue,
            other => Some(other),
        };
        if let Some(level) = level {
            let level = uint(level, what)?;
            out.insert(ell, level.try_into().map_err(|_| format!("{what}: level"))?);
        }
    }
    Ok(out)
}

/// The `curve_nodes` of an isogeny-walk `isogeny_routes.json`.
pub fn read_walk(spec: &SourceSpec, text: &str) -> Result<Vec<RegistryCurve>, String> {
    let doc: Value = serde_json::from_str(text).map_err(|e| format!("{}: {e}", spec.path))?;
    let nodes = doc["curve_nodes"]
        .as_array()
        .ok_or_else(|| format!("{}: no curve_nodes", spec.path))?;
    let mut out = Vec::with_capacity(nodes.len());
    for node in nodes {
        let slug = node["icv1_slug"]
            .as_str()
            .ok_or_else(|| format!("{}: a node has no icv1_slug", spec.path))?;
        let what = format!("{}: {slug}", spec.path);
        let coeffs = node["coefficients_a1_a2_a3_a4_a6"]
            .as_array()
            .filter(|c| c.len() == 5)
            .ok_or_else(|| format!("{what}: coefficients"))?;
        for c in &coeffs[..3] {
            if !uint(c, &what)?.is_zero() {
                return Err(format!("{what}: not a short Weierstrass model"));
            }
        }
        let p = uint(&node["field_characteristic"], &what)?;
        let (a, b) = (uint(&coeffs[3], &what)?, uint(&coeffs[4], &what)?);
        let order = uint(&node["group_order"], &what)?;
        let id = curve_id::prime(&p, &a, &b, &order)
            .ok_or_else(|| format!("{what}: not a curve with that order"))?;
        if id.slug != slug {
            return Err(format!("{what}: the model and order give {}", id.slug));
        }
        let g = node["generator_G"]
            .as_array()
            .filter(|g| g.len() == 2)
            .ok_or_else(|| format!("{what}: generator"))?;
        out.push(RegistryCurve {
            slug: slug.to_string(),
            family: spec.family.clone(),
            trace: BigInt::from(p.clone()) + 1u8 - BigInt::from(order.clone()),
            model: Model::Prime {
                p: p.clone(),
                a: a % &p,
                b: b % &p,
            },
            order,
            end: None,
            representations: vec![Representation {
                generator: GeneratorCoords::Affine(uint(&g[0], &what)?, uint(&g[1], &what)?),
                subgroup_order: uint(&node["subgroup_order"], &what)?,
                cofactor: uint(&node["cofactor"], &what)?,
            }],
            levels: levels(&node["proved_volcano_levels"], &what)?,
            source: source(spec, node["curve_uid"].as_str()),
        });
    }
    Ok(out)
}

// ── A YAML subset ─────────────────────────────────────────────────────

/// Parse the block-style YAML subset the crosswalk is written in: nested
/// mappings, `- scalar` sequences, `[]`/`{}`, `#` comments, plain or
/// single/double-quoted scalars.  Every scalar but `null`, `true` and
/// `false` is kept as a string, so large integers stay exact.  Anything
/// outside the subset is an error, never a silent misreading.  A key
/// repeated in one mapping keeps its first value, and its dotted path is
/// returned so a reader can refuse an entry that has one.
pub fn parse_yaml(text: &str) -> Result<(Value, Vec<String>), String> {
    let lines: Vec<(usize, usize, &str)> = text
        .lines()
        .enumerate()
        .filter_map(|(i, l)| {
            let body = l.trim_end();
            let content = body.trim_start();
            if content.is_empty() || content.starts_with('#') {
                return None;
            }
            Some((i + 1, body.len() - content.len(), content))
        })
        .collect();
    let mut pos = 0;
    let mut dups = Vec::new();
    let v = block(&lines, &mut pos, 0, "", &mut dups)?;
    if pos != lines.len() {
        return Err(format!("line {}: unexpected indentation", lines[pos].0));
    }
    Ok((v, dups))
}

fn scalar(s: &str, line: usize) -> Result<Value, String> {
    let s = s.trim();
    Ok(match s {
        "null" | "~" | "" => Value::Null,
        "true" => Value::Bool(true),
        "false" => Value::Bool(false),
        "[]" => Value::Array(Vec::new()),
        "{}" => Value::Object(Default::default()),
        _ if s.len() >= 2
            && ((s.starts_with('\'') && s.ends_with('\''))
                || (s.starts_with('"') && s.ends_with('"'))) =>
        {
            Value::String(s[1..s.len() - 1].to_string())
        }
        _ if s.starts_with(['[', '{', '|', '>', '&', '*', '!', '\'', '"']) => {
            return Err(format!("line {line}: unsupported YAML scalar {s}"))
        }
        _ => Value::String(s.to_string()),
    })
}

fn split_key(content: &str, line: usize) -> Result<(String, &str), String> {
    let (key, rest) = if let Some(q) = content.strip_prefix('\'') {
        let end = q
            .find('\'')
            .ok_or(format!("line {line}: unterminated key"))?;
        (&q[..end], q[end + 1..].strip_prefix(':'))
    } else {
        match content.find(':') {
            Some(i) => (&content[..i], Some(&content[i + 1..])),
            None => (content, None),
        }
    };
    let rest = rest.ok_or(format!("line {line}: expected `key:`"))?;
    if !rest.is_empty() && !rest.starts_with(' ') {
        return Err(format!("line {line}: expected a space after `:`"));
    }
    Ok((key.to_string(), rest.trim()))
}

fn block(
    lines: &[(usize, usize, &str)],
    pos: &mut usize,
    indent: usize,
    path: &str,
    dups: &mut Vec<String>,
) -> Result<Value, String> {
    let Some(&(_, first, content)) = lines.get(*pos) else {
        return Ok(Value::Null);
    };
    if first < indent {
        return Ok(Value::Null);
    }
    let indent = first;
    if content.starts_with("- ") || content == "-" {
        let mut seq = Vec::new();
        while let Some(&(line, ind, c)) = lines.get(*pos) {
            if ind != indent {
                break;
            }
            let item = c
                .strip_prefix('-')
                .ok_or(format!("line {line}: expected `- ` in a sequence"))?
                .trim();
            let quoted = item.starts_with(['\'', '"']);
            if item.is_empty() || (!quoted && (item.contains(": ") || item.ends_with(':'))) {
                return Err(format!(
                    "line {line}: nested sequence items are not supported"
                ));
            }
            seq.push(scalar(item, line)?);
            *pos += 1;
        }
        return Ok(Value::Array(seq));
    }
    let mut map = serde_json::Map::new();
    while let Some(&(line, ind, c)) = lines.get(*pos) {
        if ind < indent {
            break;
        }
        if ind > indent {
            return Err(format!("line {line}: unexpected indentation"));
        }
        let (key, rest) = split_key(c, line)?;
        *pos += 1;
        let child = format!("{path}{key}.");
        let value = if rest.is_empty() {
            match lines.get(*pos) {
                Some(&(_, next, nc))
                    if next > indent || (next == indent && nc.starts_with("- ")) =>
                {
                    block(lines, pos, next, &child, dups)?
                }
                _ => Value::Null,
            }
        } else {
            scalar(rest, line)?
        };
        if map.contains_key(&key) {
            dups.push(format!("{path}{key}"));
        } else {
            map.insert(key, value);
        }
    }
    Ok(Value::Object(map))
}

// ── The IC crosswalk ──────────────────────────────────────────────────

/// The crosswalk entries that are not registry models.  A binary model
/// `y² + xy = x³ + a₂x² + a₄x + a₆` is taken to the registry form by
/// `y ↦ y + a₄`, which gives `b = a₆ + a₄²` and moves `G` to
/// `(x, y + a₄)`.
pub fn read_crosswalk(spec: &SourceSpec, text: &str) -> Result<Vec<RegistryCurve>, String> {
    let (doc, dups) = parse_yaml(text).map_err(|e| format!("{}: {e}", spec.path))?;
    let curves = doc["curves"]
        .as_object()
        .ok_or_else(|| format!("{}: no curves", spec.path))?;
    let mut out = Vec::new();
    for (name, c) in curves {
        let what = format!("{}: {name}", spec.path);
        if c["icv1_identity"]["slug"].is_string() {
            continue;
        }
        let prefix = format!("curves.{name}.");
        if let Some(d) = dups.iter().find(|d| d.starts_with(&prefix)) {
            return Err(format!("{what}: duplicate key {d}"));
        }
        let field = &c["field"];
        let curve = &c["curve"];
        if field["p"].as_str() != Some("2") {
            return Err(format!("{what}: only binary crosswalk entries are read"));
        }
        let n: u32 = uint(&field["n"], &what)?
            .try_into()
            .map_err(|_| format!("{what}: n"))?;
        let terms = field["defining_modulus_low_terms"]
            .as_array()
            .ok_or_else(|| format!("{what}: modulus"))?;
        let mut modulus = BigUint::zero();
        for t in terms {
            let e: u64 = uint(t, &what)?
                .try_into()
                .map_err(|_| format!("{what}: modulus term"))?;
            modulus.set_bit(e, true);
        }
        let irr = irreducible(&modulus);
        let coeffs: Vec<BigUint> = curve["coefficients_a1_a2_a3_a4_a6"]
            .as_array()
            .filter(|v| v.len() == 5)
            .ok_or_else(|| format!("{what}: coefficients"))?
            .iter()
            .map(|v| uint(v, &what))
            .collect::<Result<_, _>>()?;
        if coeffs[0] != BigUint::from(1u8) || !coeffs[2].is_zero() {
            return Err(format!(
                "{what}: only a₁ = 1, a₃ = 0 binary models are read"
            ));
        }
        let el = |v: &BigUint| F2mElement::from_biguint(v, n);
        let a4 = el(&coeffs[3]);
        let b = el(&coeffs[4]).add(&a4.square(&irr));
        let order = uint(&curve["curve_order"], &what)?;
        let trace = int(&curve["trace"], &what)?;
        if BigInt::from(BigUint::from(1u8) << n) + 1u8 - &trace != BigInt::from(order.clone()) {
            return Err(format!("{what}: order ≠ q + 1 − t"));
        }
        let a = coeffs[1].clone();
        let id = curve_id::binary(n, &modulus, &a, &b.to_biguint(), &order, None)
            .ok_or_else(|| format!("{what}: not a curve with that order"))?;
        let g = curve["generator_G"]
            .as_array()
            .filter(|g| g.len() == 2)
            .ok_or_else(|| format!("{what}: generator"))?;
        let gy = el(&uint(&g[1], &what)?).add(&a4).to_biguint();
        let endo = &c["endomorphism"];
        let end = match &c["trait_status"]["endomorphism_discriminant"] {
            d if d["status"].as_str() == Some("proved") => Some(int(&d["value"], &what)?),
            _ => None,
        };
        out.push(RegistryCurve {
            slug: id.slug,
            family: spec.family.clone(),
            model: Model::Binary {
                n,
                modulus,
                irr,
                a: el(&a),
                b,
            },
            trace,
            order,
            end,
            representations: vec![Representation {
                generator: GeneratorCoords::Affine(uint(&g[0], &what)?, gy),
                subgroup_order: uint(&curve["subgroup_order"], &what)?,
                cofactor: uint(&curve["cofactor"], &what)?,
            }],
            levels: levels(&endo["volcano_levels"], &what)?,
            source: source(spec, c["curve_uid"].as_str()),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_yaml_subset_round_trips_what_the_crosswalk_uses() {
        let (v, dups) = parse_yaml(
            "# c\na: 1\nb:\n  c: 'x y'\n  d:\n    - 3\n    - 4\n  e: []\n  '263':\n    level: 1\nf: null\nb:\n  c: 2\n",
        )
        .unwrap();
        assert_eq!(dups, ["b"]);
        assert_eq!(v["a"], "1");
        assert_eq!(v["b"]["c"], "x y");
        assert_eq!(v["b"]["d"][1], "4");
        assert_eq!(v["b"]["e"], Value::Array(vec![]));
        assert_eq!(v["b"]["263"]["level"], "1");
        assert!(v["f"].is_null());
        assert!(parse_yaml("a: |\n  x\n").is_err());
        assert!(parse_yaml("a:\n  - b: 1\n").is_err());
    }
}
