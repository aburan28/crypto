//! **ICV1 curve identities** — the one way this repository names a curve.
//!
//! The specification is `docs/curves/ICV1.md` and the reference
//! implementation is `scripts/curve_id.py`; this module reproduces it byte
//! for byte, and `tests/curve_id.rs` pins vectors that file produced.  An
//! identity is
//!
//! ```text
//! ICV1:<field>:<trace>:<order>:<j>:<end>:<level>:<path>:<model12>
//! ```
//!
//! with a slug `icv1-<field tag>-t<trace>-<model8>` (a negative trace is
//! written `tm<|t|>`), which is the name a curve carries in text, tables and
//! reports.  `<model12>` is the SHA-256 of the curve's canonical model JSON,
//! so two models of one curve — another modulus, another Weierstrass form —
//! are two identities.
//!
//! Reports written before ICV1 name curves by the spellings they used then
//! (`K_0 / GF(2^41)`, `bench-24bit`, `random-binary-n27-b845462`, …).
//! Those files are frozen.  The registry (`docs/curves/registry.json`) maps
//! every such spelling to its identity; [`resolve`] and [`same_curve`] read
//! the alias map generated from it (`curve_aliases.json`, beside this file),
//! so code that replays a frozen report, or looks a curve up in a table keyed
//! the old way, accepts either name.

use std::collections::HashMap;
use std::sync::OnceLock;

use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};

use crate::binary_ecc::IrreduciblePoly;
use crate::hash::sha256::sha256;

/// Version of the canonical model JSON.
pub const VERSION: &str = "1";

/// Every name the registry knows, normalised, to its slug: generated from
/// `docs/curves/registry.json` by `scripts/build_curve_registry.py` and kept
/// under `src/` so that a checkout of `src/` alone still builds.
const ALIASES: &str = include_str!("curve_aliases.json");

/// A curve's ICV1 identity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CurveId {
    /// `ICV1:<field>:<trace>:<order>:<j>:<end>:<level>:<path>:<model12>`.
    pub icv1: String,
    /// `icv1-<field tag>-t<trace>-<model8>`, the name used in text.
    pub slug: String,
    /// The canonical model JSON the model hash is taken over.
    pub model_json: String,
}

fn hex_digest(s: &str) -> String {
    sha256(s.as_bytes())
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn hex(v: &BigUint) -> String {
    format!("0x{v:x}")
}

/// `x^m + Σ x^i` over the low terms, as the integer whose bit `i` is the
/// coefficient of `x^i`.
pub fn modulus_integer(p: &IrreduciblePoly) -> BigUint {
    let mut f = BigUint::zero();
    f.set_bit(u64::from(p.degree), true);
    for &i in &p.low_terms {
        f.set_bit(u64::from(i), true);
    }
    f
}

fn deg(f: &BigUint) -> u64 {
    f.bits().saturating_sub(1)
}

fn poly_mod(mut a: BigUint, f: &BigUint) -> BigUint {
    let df = deg(f);
    while !a.is_zero() && deg(&a) >= df {
        a ^= f << (deg(&a) - df);
    }
    a
}

fn clmul(a: &BigUint, b: &BigUint) -> BigUint {
    let mut r = BigUint::zero();
    for i in 0..b.bits() {
        if b.bit(i) {
            r ^= a << i;
        }
    }
    r
}

/// Inverse of a non-zero `a` modulo the irreducible `f` over `GF(2)`.
fn poly_inverse(a: &BigUint, f: &BigUint) -> Option<BigUint> {
    let (mut r0, mut r1) = (f.clone(), poly_mod(a.clone(), f));
    let (mut s0, mut s1) = (BigUint::zero(), BigUint::one());
    if r1.is_zero() {
        return None;
    }
    while !r1.is_one() {
        let mut q = BigUint::zero();
        let mut r = r0.clone();
        while !r.is_zero() && deg(&r) >= deg(&r1) {
            let shift = deg(&r) - deg(&r1);
            q.set_bit(shift, true);
            r ^= &r1 << shift;
        }
        if r.is_zero() {
            return None; // f was not irreducible
        }
        let s = &s0 ^ clmul(&q, &s1);
        r0 = std::mem::replace(&mut r1, r);
        s0 = std::mem::replace(&mut s1, s);
    }
    Some(poly_mod(s1, f))
}

/// `|q + 1 − #E| ≤ 2√q`: false for a subgroup order passed as the group's.
fn hasse_ok(q: &BigUint, order: &BigUint) -> bool {
    let t = BigInt::from(q.clone()) + 1 - BigInt::from(order.clone());
    &t * &t <= BigInt::from(q.clone()) * 4
}

fn identity(
    field: &str,
    field_tag: &str,
    trace: &BigInt,
    order: &BigUint,
    j: &str,
    end: &str,
    model_json: String,
) -> CurveId {
    let model = hex_digest(&model_json);
    let icv1 = format!(
        "ICV1:{field}:{trace}:{order}:{j}:{end}:unk:r:{}",
        &model[..12]
    );
    let t = if trace < &BigInt::zero() {
        format!("tm{}", -trace)
    } else {
        format!("t{trace}")
    };
    let slug = format!("icv1-{field_tag}-{t}-{}", &model[..8]);
    CurveId {
        icv1,
        slug,
        model_json,
    }
}

/// `y² + xy = x³ + a·x² + b` over `GF(2^m) = GF(2)[x]/(modulus)`.
///
/// `a`, `b` are polynomial-basis elements (bit `i` is `x^i`) and `order`
/// is `#E(GF(2^m))`, the whole group.  `end` is the certified
/// discriminant of `End(E)`, or `None` for `unk`.  Returns `None` for a
/// singular curve, a reducible modulus or an order outside the Hasse
/// interval (a subgroup order passed by mistake).
pub fn binary(
    m: u32,
    modulus: &BigUint,
    a: &BigUint,
    b: &BigUint,
    order: &BigUint,
    end: Option<i64>,
) -> Option<CurveId> {
    if deg(modulus) != u64::from(m) {
        return None;
    }
    let (a, b) = (poly_mod(a.clone(), modulus), poly_mod(b.clone(), modulus));
    let q = BigUint::one() << m;
    if b.is_zero() || !hasse_ok(&q, order) {
        return None;
    }
    let modhash = hex_digest(&format!("f2m-modulus:{}", hex(modulus)));
    let field = format!("f2m-{m}-{}", &modhash[..8]);
    let model_json = format!(
        r#"{{"a":"{}","b":"{}","field":"{field}","form":"y^2+xy=x^3+a*x^2+b","modulus":"{}","v":"{VERSION}"}}"#,
        hex(&a),
        hex(&b),
        hex(modulus)
    );
    let trace = BigInt::from(q) + 1 - BigInt::from(order.clone());
    let j = hex(&poly_inverse(&b, modulus)?);
    let end = end.map_or_else(|| "unk".to_string(), |d| d.to_string());
    Some(identity(
        &field,
        &format!("f2m{m}"),
        &trace,
        order,
        &j,
        &end,
        model_json,
    ))
}

/// The Koblitz curve `K_a : y² + xy = x³ + a·x² + 1` over `GF(2^n)`.
///
/// `End(E) ⊇ Z[τ]` with `τ² − t₁τ + 2 = 0`, `t₁ = ±1`, whose
/// discriminant `−7` is fundamental: `Z[τ]` is already the maximal order
/// of `Q(√−7)`, so `End(E) = Z[τ]` and the identity certifies `end = −7`.
pub fn koblitz(a: u8, n: u32, modulus: &BigUint, order: &BigUint) -> Option<CurveId> {
    if a > 1 {
        return None;
    }
    binary(
        n,
        modulus,
        &BigUint::from(a),
        &BigUint::one(),
        order,
        Some(-7),
    )
}

/// `y² = x³ + a·x + b` over `GF(p)`, `p` an odd prime above 3.
pub fn prime(p: &BigUint, a: &BigUint, b: &BigUint, order: &BigUint) -> Option<CurveId> {
    let (a, b) = (a % p, b % p);
    let a3 = a.modpow(&BigUint::from(3u8), p);
    let disc = (BigUint::from(4u8) * &a3 + BigUint::from(27u8) * &b * &b) % p;
    if disc.is_zero() || !hasse_ok(p, order) {
        return None;
    }
    let field = format!("fp-{p}");
    let model_json = format!(
        r#"{{"a":"{a}","b":"{b}","field":"{field}","form":"y^2=x^3+a*x+b","p":"{p}","v":"{VERSION}"}}"#
    );
    let trace = BigInt::from(p.clone()) + 1 - BigInt::from(order.clone());
    let inv = disc.modpow(&(p - BigUint::from(2u8)), p);
    let j = (BigUint::from(1728u32 * 4) * a3 * inv) % p;
    Some(identity(
        &field,
        &format!("fp{}", p.bits()),
        &trace,
        order,
        &j.to_string(),
        "unk",
        model_json,
    ))
}

/// The comparison key for a spelling: case, whitespace, braces,
/// underscores, surrounding backticks and subscript digits folded away.
/// Identical to `normalise_alias` in `scripts/curve_id.py`.
pub fn normalise(name: &str) -> String {
    name.trim()
        .trim_matches('`')
        .chars()
        .filter(|c| !c.is_whitespace() && !matches!(c, '{' | '}' | '_'))
        .map(|c| match c {
            '₀'..='₉' => char::from(b'0' + (c as u32 - '₀' as u32) as u8),
            _ => c.to_ascii_lowercase(),
        })
        .collect()
}

/// Every name the registry knows, normalised, mapped to its curve's slug.
fn aliases() -> &'static HashMap<String, String> {
    static MAP: OnceLock<HashMap<String, String>> = OnceLock::new();
    MAP.get_or_init(|| {
        let doc: serde_json::Value =
            serde_json::from_str(ALIASES).expect("curve_aliases.json parses");
        doc["aliases"]
            .as_object()
            .into_iter()
            .flatten()
            .filter_map(|(k, v)| Some((k.clone(), v.as_str()?.to_string())))
            .collect()
    })
}

/// The slug of the curve `name` denotes: a slug, an ICV1 string, a
/// standard name or a registered legacy spelling.  `None` when the
/// registry does not know the name.
pub fn resolve(name: &str) -> Option<&'static str> {
    aliases().get(&normalise(name)).map(String::as_str)
}

/// Whether two names denote one curve: equal as written, or resolved by
/// the registry to one slug.  This is how a rebuilt instance is matched
/// against the name a frozen report recorded.
pub fn same_curve(a: &str, b: &str) -> bool {
    if a == b {
        return true;
    }
    let ra = resolve(a).or_else(|| a.starts_with("icv1-").then_some(a));
    let rb = resolve(b).or_else(|| b.starts_with("icv1-").then_some(b));
    matches!((ra, rb), (Some(x), Some(y)) if x == y)
}

/// The degree of the binary field a curve name carries: `m` of an ICV1
/// slug's `f2m<m>`, `n` of a legacy `…-n<n>-…` or `…GF(2^<n>)`.
pub fn binary_degree(name: &str) -> Option<u32> {
    let digits = |s: &str| -> Option<u32> {
        let d: String = s.chars().take_while(char::is_ascii_digit).collect();
        d.parse().ok()
    };
    if let Some(rest) = name.strip_prefix("icv1-f2m") {
        return digits(rest);
    }
    if let Some(at) = name.find("2^") {
        return digits(&name[at + 2..]);
    }
    let at = name.find("-n")? + 2;
    digits(&name[at..])
}
