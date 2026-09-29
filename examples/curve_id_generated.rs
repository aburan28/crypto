//! Rebuild the curves the repository generated from a seed, and record
//! each one's model, its ICV1 identity and the name it was reported under.
//!
//! Frozen reports name these curves by handles the generators used before
//! ICV1 (`generated-<bits>bit-<p>`, `random-binary-n<n>-b<b>`,
//! `E_{a,b}/GF(2^k) over GF(2^n)`), and most record only the seed, not the
//! coefficients.  This example regenerates them so that the curve registry
//! (`scripts/build_curve_registry.py`) can resolve those handles.  Its
//! output is committed as `docs/curves/sources/generated.json`; see
//! `docs/curves/README.md` for the command.
//!
//! Specifications, one per argument:
//!
//! - `prime:<bits>:<seed>` — `find_prime_order_curve(bits, seed)`;
//! - `char2:<n>:<seed>:<max cofactor>` — `random_binary_instance`;
//! - `subfield:<k>:<n>:<a index>:<b index>` — `KoblitzCurve::subfield`.

use crypto_lib::cryptanalysis::curve_id::modulus_integer;
use crypto_lib::cryptanalysis::ic_boundary::{find_prime_order_curve, random_binary_instance};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use serde_json::{json, Value};

fn field(n: u32, low_terms: &[u32]) -> Value {
    json!({"kind": "binary", "degree": n, "polynomial_low_terms": low_terms})
}

fn one(spec: &str) -> Result<Value, String> {
    let parts: Vec<&str> = spec.split(':').collect();
    let num = |i: usize| -> Result<u64, String> {
        parts
            .get(i)
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| format!("{spec}: argument {i} is not a number"))
    };
    match parts.first().copied() {
        Some("prime") => {
            let (bits, seed) = (num(1)? as u32, num(2)?);
            let inst = find_prime_order_curve(bits, seed);
            let id = inst.curve_id();
            Ok(json!({
                "name": format!("generated-{bits}bit-{}", inst.curve.p),
                "generator_call": format!("find_prime_order_curve({bits}, {seed})"),
                "p": inst.curve.p, "a": inst.curve.a, "b": inst.curve.b,
                "group_order": inst.group_order, "subgroup_order": inst.r,
                "icv1": id.icv1, "slug": id.slug,
            }))
        }
        Some("char2") => {
            let (n, seed, cofactor) = (num(1)? as u32, num(2)?, num(3)?);
            let inst = random_binary_instance(n, seed, cofactor)
                .ok_or_else(|| format!("{spec}: no curve"))?;
            let id = inst.curve_id();
            Ok(json!({
                "name": format!("random-binary-n{n}-b{:x}", inst.b),
                "generator_call": format!("random_binary_instance({n}, {seed}, {cofactor})"),
                "field": field(n, &inst.irreducible.low_terms),
                "a": format!("0x{:x}", inst.a), "b": format!("0x{:x}", inst.b),
                "group_order": inst.group_order, "subgroup_order": inst.r,
                "koblitz": false, "icv1": id.icv1, "slug": id.slug,
            }))
        }
        Some("subfield") => {
            let (k, n, a, b) = (num(1)? as u32, num(2)? as u32, num(3)?, num(4)?);
            let c = KoblitzCurve::subfield(k, n, a, b).ok_or_else(|| format!("{spec}: no curve"))?;
            let id = c.curve_id();
            Ok(json!({
                "name": format!("E_{{{a},{b}}}/GF(2^{k}) over GF(2^{n})"),
                "generator_call": format!("KoblitzCurve::subfield({k}, {n}, {a}, {b})"),
                "field": field(n, &c.curve.irreducible.low_terms),
                "modulus": format!("0x{:x}", modulus_integer(&c.curve.irreducible)),
                "a": format!("0x{:x}", c.curve.a.to_biguint()),
                "b": format!("0x{:x}", c.curve.b.to_biguint()),
                "group_order": c.group_order.to_string(),
                "subgroup_order": c.subgroup_order.to_string(),
                "koblitz": false, "icv1": id.icv1, "slug": id.slug,
            }))
        }
        _ => Err(format!("{spec}: expected prime:, char2: or subfield:")),
    }
}

fn main() {
    let specs: Vec<String> = std::env::args().skip(1).collect();
    let mut curves = Vec::new();
    for spec in &specs {
        match one(spec) {
            Ok(v) => curves.push(v),
            Err(e) => {
                eprintln!("{e}");
                std::process::exit(1);
            }
        }
    }
    let doc = json!({
        "what_this_is": "Curves the repository generated from a seed, rebuilt by \
            examples/curve_id_generated.rs, with the handle each was reported under \
            before ICV1 and its ICV1 identity computed in Rust.",
        "specs": specs,
        "curves": curves,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
}
