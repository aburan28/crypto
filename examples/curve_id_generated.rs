//! Rebuild the curves the repository's constructors build, and record each
//! one's exact representation (modulus, coefficients, subgroup, cofactor,
//! generator), its ICV1 identity and the name it was reported under.
//!
//! Frozen reports name these curves by handles the generators used before
//! ICV1 (`generated-<bits>bit-<p>`, `random-binary-n<n>-b<b>`,
//! `E_{a,b}/GF(2^k) over GF(2^n)`, `K_a / GF(2^n)`, `bench-<bits>bit`), and
//! most record only the seed or the parameters, not the coefficients or the
//! generator.  This example rebuilds them so that the curve registry
//! (`scripts/build_curve_registry.py`) can resolve those handles and attach
//! each representation's EC1 identity (`docs/curve-identities.md`), which
//! needs the generator.  Its output is committed as
//! `docs/curves/sources/generated.json`; see `docs/curves/README.md` for the
//! command.
//!
//! Specifications, one per argument:
//!
//! - `prime:<bits>:<seed>` — `find_prime_order_curve(bits, seed)`;
//! - `char2:<n>:<seed>:<max cofactor>` — `random_binary_instance`;
//! - `subfield:<k>:<n>:<a index>:<b index>` — `KoblitzCurve::subfield`;
//! - `koblitz:<a>:<n>` — `KoblitzCurve::new(a, n)`;
//! - `roster:<bits>` — `roster_prime_instance(bits)`, the bench roster.

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::curve_id::modulus_integer;
use crypto_lib::cryptanalysis::ic_boundary::{
    find_prime_order_curve, random_binary_instance, roster_prime_instance, PrimeInstance,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use serde_json::{json, Value};

fn field(n: u32, low_terms: &[u32]) -> Value {
    json!({"kind": "binary", "degree": n, "polynomial_low_terms": low_terms})
}

fn prime_record(inst: &PrimeInstance, name: String, call: String) -> Value {
    let id = inst.curve_id();
    json!({
        "name": name, "generator_call": call,
        "p": inst.curve.p, "a": inst.curve.a, "b": inst.curve.b,
        "group_order": inst.group_order, "subgroup_order": inst.r, "cofactor": inst.cofactor,
        "generator": {"x": inst.generator.0, "y": inst.generator.1},
        "icv1": id.icv1, "slug": id.slug,
    })
}

fn koblitz_record(c: &KoblitzCurve, name: String, call: String) -> Result<Value, String> {
    let id = c.curve_id();
    let BinaryPoint::Affine { x, y } = c.generator() else {
        return Err(format!("{call}: generator at infinity"));
    };
    Ok(json!({
        "name": name, "generator_call": call,
        "field": field(c.n, &c.curve.irreducible.low_terms),
        "a": format!("0x{:x}", c.curve.a.to_biguint()),
        "b": format!("0x{:x}", c.curve.b.to_biguint()),
        "group_order": c.group_order.to_string(),
        "subgroup_order": c.subgroup_order.to_string(),
        "cofactor": c.cofactor.to_string(),
        "generator": {"x": format!("0x{:x}", x.to_biguint()), "y": format!("0x{:x}", y.to_biguint())},
        "koblitz": c.k == 1, "icv1": id.icv1, "slug": id.slug,
    }))
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
            Ok(prime_record(
                &inst,
                format!("generated-{bits}bit-{}", inst.curve.p),
                format!("find_prime_order_curve({bits}, {seed})"),
            ))
        }
        Some("roster") => {
            let bits = num(1)? as u32;
            let inst = roster_prime_instance(bits).ok_or_else(|| format!("{spec}: no curve"))?;
            Ok(prime_record(
                &inst,
                format!("bench-{bits}bit"),
                format!("roster_prime_instance({bits})"),
            ))
        }
        Some("koblitz") => {
            let (a, n) = (num(1)? as u8, num(2)? as u32);
            let c = KoblitzCurve::new(a, n).ok_or_else(|| format!("{spec}: no curve"))?;
            koblitz_record(
                &c,
                format!("K_{a} / GF(2^{n})"),
                format!("KoblitzCurve::new({a}, {n})"),
            )
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
                "cofactor": inst.cofactor,
                "generator": {"x": format!("0x{:x}", inst.generator.x),
                              "y": format!("0x{:x}", inst.generator.y)},
                "koblitz": false, "icv1": id.icv1, "slug": id.slug,
            }))
        }
        Some("subfield") => {
            let (k, n, a, b) = (num(1)? as u32, num(2)? as u32, num(3)?, num(4)?);
            let c =
                KoblitzCurve::subfield(k, n, a, b).ok_or_else(|| format!("{spec}: no curve"))?;
            let mut v = koblitz_record(
                &c,
                format!("E_{{{a},{b}}}/GF(2^{k}) over GF(2^{n})"),
                format!("KoblitzCurve::subfield({k}, {n}, {a}, {b})"),
            )?;
            v["modulus"] = json!(format!("0x{:x}", modulus_integer(&c.curve.irreducible)));
            Ok(v)
        }
        _ => Err(format!(
            "{spec}: expected prime:, roster:, char2:, subfield: or koblitz:"
        )),
    }
}

fn main() {
    let specs: Vec<String> = std::env::args().skip(1).collect();
    let (mut curves, mut unbuilt) = (Vec::new(), Vec::new());
    for spec in &specs {
        match one(spec) {
            Ok(v) => curves.push(v),
            Err(e) => unbuilt.push(json!({"spec": spec, "why": e})),
        }
    }
    let doc = json!({
        "what_this_is": "Curves the repository's constructors build, rebuilt by \
            examples/curve_id_generated.rs: each one's exact representation, the handle \
            it was reported under before ICV1, and its ICV1 identity computed in Rust.",
        "specs": specs,
        "curves": curves,
        "unbuilt": unbuilt,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
}
