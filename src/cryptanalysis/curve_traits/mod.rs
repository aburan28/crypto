//! # Curve traits: size-independent signatures for registered curves
//!
//! The registry (`docs/curves/registry.json`, `AGENTS.md` §11) names every
//! curve the repository uses and records its model and group order.  It
//! does not say which curves are *alike*: that a family of Koblitz curves
//! shares one CM field and one descended Frobenius at every degree, that a
//! curve over `GF(2^18)` is the quadratic twist of one, or how far `Z[π]`
//! sits below the maximal order.  This module derives those traits from
//! each model and order, records how well each is established, and turns
//! them into grouping keys, so curves can be grouped by structure across
//! sizes (`curve_traits group`) or ranked by likeness to one curve
//! (`curve_traits similar`).
//!
//! ## What is computed
//!
//! - **Order check** ([`certify`]): exhaustive count on small fields, a
//!   generator certificate otherwise.
//! - **Frobenius discriminant** `Δ = t² − 4q = v²·d_K`: its factorisation,
//!   the CM discriminant `d_K` of `Q(π)`, the conductor `v = [O_K : Z[π]]`
//!   (the "conductor gap" between `Z[π]` and the maximal order) and
//!   `h(d_K)` for small `|d_K|`.  `End(E)` lies between `Z[π]` and `O_K`;
//!   which order it is needs a volcano walk and is not decided here.
//! - **Small primes** `ℓ ≤ 31`: the depth `v_ℓ(v)` of the `ℓ`-isogeny
//!   volcano and how `ℓ` splits in `Q(π)`.  Exact without factoring `Δ`.
//! - **Subfields** ([`subfield`], binary only): the field of `j`, the
//!   smallest field of definition and the trace of the descended model.
//! - **Subgroup, embedding degree, twist**: the prime subgroup and
//!   cofactor, the embedding degree `ord_r(q)`, the twist's cofactor.
//!
//! ## Status of a value
//!
//! Each derived value carries a [`Status`].  Factoring is budgeted, so a
//! large discriminant may stay partly unfactored; the record then holds a
//! bound or nothing, never a guess.  The trace itself is the registry's;
//! the order check says how well it is established.
//!
//! Records are written to `docs/curves/traits.json` by
//! `cargo run --release --bin curve_traits -- build`; see
//! `docs/curves/TRAITS.md`.

pub mod arith;
pub mod certify;
pub mod dissect;
pub mod query;
pub mod sources;
pub mod subfield;
pub mod torsion;

use std::collections::{BTreeMap, HashMap};
use std::sync::Mutex;

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use arith::Factored;

/// How well a value is established.  Ordered from strongest to weakest,
/// so the weaker of two is their maximum.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Status {
    /// Exact, and every primality claim under it is deterministic.
    Proved,
    /// Exact if every probable prime under it (Miller–Rabin, bases
    /// 2…53, above 3.3·10²⁴) is prime.
    Probable,
    /// Only a bound is known; the record says which.
    Bounded,
    /// Attempted and not determined within the budget.
    Unknown,
    /// Not attempted: beyond a declared limit.
    NotEvaluated,
    /// Does not apply to this curve.
    NotApplicable,
}

/// The rho iterations each composite gets by default.
pub const DEFAULT_BUDGET: u64 = 1 << 18;

/// `h(d_K)` is computed for `|d_K|` up to this bound.
pub const CLASS_NUMBER_MAX_DISC: u64 = 10_000_000;

/// A CM discriminant at most this large in absolute value is its own
/// grouping key; above it, keys say `large`.
pub const SMALL_CM_DISC: u64 = 1_000_000;

/// An embedding degree at most this is its own grouping key.
pub const SMALL_EMBEDDING_DEGREE: u64 = 20;

/// Powers of `q` tried when `r − 1` does not factor.
const EMBEDDING_SEARCH: u64 = 1000;

/// The small primes `ℓ` whose local data every record carries.
pub const SMALL_PRIMES: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];

/// A registered curve model.
#[derive(Clone, Debug)]
pub enum Model {
    /// `y² + xy = x³ + ax² + b` over `GF(2^n) = GF(2)[z]/(f)`.
    Binary {
        n: u32,
        modulus: BigUint,
        irr: IrreduciblePoly,
        a: F2mElement,
        b: F2mElement,
    },
    /// `y² = x³ + ax + b` over `GF(p)`.
    Prime { p: BigUint, a: BigUint, b: BigUint },
    /// A recognized catalogue model whose coordinate arithmetic is unavailable.
    /// Retain its exact field so order/trace identities can still be checked.
    Other { p: BigUint, degree: u32 },
    /// Short Weierstrass model over `GF(p^k)`, in the registered polynomial basis.
    Extension {
        p: BigUint,
        k: u32,
        q: BigUint,
        field_id: String,
        modulus: Vec<BigUint>,
        a: Vec<BigUint>,
        b: Vec<BigUint>,
    },
}

impl Model {
    /// The field size `q`.
    pub fn q(&self) -> BigUint {
        match self {
            Model::Binary { n, .. } => BigUint::one() << *n,
            Model::Prime { p, .. } => p.clone(),
            Model::Other { p, degree } => p.pow(*degree),
            Model::Extension { q, .. } => q.clone(),
        }
    }

    /// `n` for `GF(2^n)`, the bit length of `p` for `GF(p)`.
    pub fn size_bits(&self) -> u64 {
        match self {
            Model::Binary { n, .. } => u64::from(*n),
            Model::Prime { p, .. } => p.bits(),
            Model::Other { .. } => self.q().bits(),
            Model::Extension { q, .. } => q.bits(),
        }
    }
}

/// Generator coordinates in the model's declared field representation.
#[derive(Clone, Debug)]
pub enum GeneratorCoords {
    Affine(BigUint, BigUint),
    Extension(Vec<BigUint>, Vec<BigUint>),
}

/// One recorded representation, mapped into `Model` coordinates: a generator
/// of the prime subgroup.
#[derive(Clone, Debug)]
pub struct Representation {
    pub generator: GeneratorCoords,
    pub subgroup_order: BigUint,
    pub cofactor: BigUint,
}

/// A registry entry, parsed losslessly, or a curve from another source
/// (see [`sources`]) in the same shape.
#[derive(Clone, Debug)]
pub struct RegistryCurve {
    pub slug: String,
    pub family: String,
    pub model: Model,
    pub trace: BigInt,
    pub order: BigUint,
    /// The certified `End(E)` discriminant, when the source has one.
    pub end: Option<BigInt>,
    pub representations: Vec<Representation>,
    /// Volcano levels `v_ℓ([O_K : End(E)])` the source proves, by `ℓ`.
    pub levels: BTreeMap<u64, u32>,
    pub source: Source,
}

/// Where a record's curve came from.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Source {
    /// `registry`, `isogeny-walk` or `ic-crosswalk`.
    pub kind: String,
    /// The file it was read from, when not the registry.
    pub path: Option<String>,
    /// The source's representation identifier, when it records one.
    pub curve_uid: Option<String>,
    /// Whether the slug is a registered ICV1 identity (`AGENTS.md` §11);
    /// a walked or crosswalk curve's slug is computed, not registered.
    pub registered: bool,
}

impl Source {
    fn registry() -> Self {
        Source {
            kind: "registry".into(),
            path: None,
            curve_uid: None,
            registered: true,
        }
    }
}

/// Factoring shared by every record of one build: the rho budget, the
/// committed factor hints, and a memo, since isogenous curves share
/// `t² − 4q`, twist orders and `r − 1`.
pub struct Factoring {
    pub budget: u64,
    pub hints: Vec<BigUint>,
    memo: Mutex<HashMap<BigUint, Factored>>,
}

impl Factoring {
    pub fn new(budget: u64, hints: Vec<BigUint>) -> Self {
        Factoring {
            budget,
            hints,
            memo: Mutex::new(HashMap::new()),
        }
    }

    pub fn factor(&self, n: &BigUint) -> Factored {
        if let Some(f) = self.memo.lock().expect("memo").get(n) {
            return f.clone();
        }
        let f = arith::factor_hinted(n, self.budget, &self.hints);
        self.memo.lock().expect("memo").insert(n.clone(), f.clone());
        f
    }

    /// Every composite left unsplit so far, ascending.
    pub fn leftovers(&self) -> Vec<BigUint> {
        let memo = self.memo.lock().expect("memo");
        let mut out: Vec<BigUint> = memo
            .values()
            .flat_map(|f| f.composites.iter().map(|(c, _)| c.clone()))
            .collect();
        out.sort();
        out.dedup();
        out
    }
}

/// Parse `docs/curves/factor-hints.txt`: one prime per line, each above
/// the trial-division bound and passing [`arith::is_prime`].
pub fn read_hints(text: &str) -> Result<Vec<BigUint>, String> {
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let p: BigUint = line
            .parse()
            .map_err(|_| format!("factor hints line {}: not an integer", i + 1))?;
        if p.bits() <= 16 || !arith::is_prime(&p) {
            return Err(format!(
                "factor hints line {}: {p} is not a prime above 2^16",
                i + 1
            ));
        }
        out.push(p);
    }
    out.sort();
    out.dedup();
    Ok(out)
}

fn parse_uint(s: &str) -> Option<BigUint> {
    match s.strip_prefix("0x") {
        Some(h) => BigUint::parse_bytes(h.as_bytes(), 16),
        None => s.parse().ok(),
    }
}

fn value_uint(value: &Value) -> Option<BigUint> {
    match value {
        Value::String(s) => parse_uint(s),
        Value::Number(n) => parse_uint(&n.to_string()),
        _ => None,
    }
}

fn extension_coordinates(value: &Value, k: u32, p: &BigUint) -> Option<Vec<BigUint>> {
    let coefficients = value.as_array()?;
    if coefficients.len() != k as usize {
        return None;
    }
    let coordinates: Vec<BigUint> = coefficients.iter().map(value_uint).collect::<Option<_>>()?;
    coordinates
        .iter()
        .all(|coordinate| coordinate < p)
        .then_some(coordinates)
}

fn str_field<'a>(v: &'a Value, key: &str) -> Result<&'a str, String> {
    v[key]
        .as_str()
        .ok_or_else(|| format!("missing string field {key}"))
}

fn irreducible(modulus: &BigUint) -> IrreduciblePoly {
    let degree = (modulus.bits() - 1) as u32;
    IrreduciblePoly {
        degree,
        low_terms: (0..degree).filter(|&i| modulus.bit(u64::from(i))).collect(),
    }
}

/// Coordinate conversion from a registered prime-field model to the
/// short Weierstrass model used for order certification.
#[derive(Clone, Debug)]
enum CoordinateMap {
    Identity,
    Montgomery {
        a_over_three: BigUint,
        inv_b: BigUint,
    },
    Edwards {
        a_over_three: BigUint,
        inv_b: BigUint,
        inv_scale: BigUint,
    },
}

fn inverse_mod(value: &BigUint, p: &BigUint) -> Option<BigUint> {
    let modulus = BigInt::from(p.clone());
    let egcd = BigInt::from(value % p).extended_gcd(&modulus);
    if egcd.gcd != BigInt::one() {
        return None;
    }
    egcd.x.mod_floor(&modulus).to_biguint()
}

fn subtract_mod(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
    (a + p - (b % p)) % p
}

/// B*v² = u³ + A*u² + u, with X=(u+A/3)/B and Y=v/B.
fn montgomery_short(
    p: &BigUint,
    a: &BigUint,
    b: &BigUint,
) -> Result<(BigUint, BigUint, BigUint, BigUint), String> {
    if p <= &BigUint::from(3u8) {
        return Err("Montgomery conversion needs characteristic above three".into());
    }
    let a = a % p;
    let b = b % p;
    let inv_three =
        inverse_mod(&BigUint::from(3u8), p).ok_or("three is not invertible in the prime field")?;
    let inv_b = inverse_mod(&b, p).ok_or("Montgomery B is not invertible")?;
    let a_over_three = &a * &inv_three % p;
    let a2 = &a * &a % p;
    let a3 = &a2 * &a % p;
    let inv_b2 = &inv_b * &inv_b % p;
    let inv_b3 = &inv_b2 * &inv_b % p;
    let inv_27 = &inv_three * &inv_three % p * &inv_three % p;
    let short_a = subtract_mod(&BigUint::one(), &(&a2 * &inv_three % p), p) * &inv_b2 % p;
    let short_b =
        subtract_mod(&(BigUint::from(2u8) * &a3 * inv_27 % p), &a_over_three, p) * &inv_b3 % p;
    Ok((short_a, short_b, a_over_three, inv_b))
}

fn edwards_short(
    p: BigUint,
    a: BigUint,
    d: BigUint,
    scale: BigUint,
) -> Result<(Model, CoordinateMap), String> {
    if p <= BigUint::from(3u8) {
        return Err("Edwards conversion needs characteristic above three".into());
    }
    let inv_scale = inverse_mod(&scale, &p).ok_or("Edwards scale is not invertible")?;
    // x'=x/c, y'=y/c turns x²+y²=c²(1+d*x²*y²)
    // into a*x'²+y'²=1+(d*c⁴)*x'²*y'².
    let d = d * scale.modpow(&BigUint::from(4u8), &p) % &p;
    let a = a % &p;
    let inv_a_minus_d =
        inverse_mod(&subtract_mod(&a, &d, &p), &p).ok_or("Edwards a-d is not invertible")?;
    let mont_a = BigUint::from(2u8) * (&a + &d) * &inv_a_minus_d % &p;
    let mont_b = BigUint::from(4u8) * inv_a_minus_d % &p;
    let (short_a, short_b, a_over_three, inv_b) = montgomery_short(&p, &mont_a, &mont_b)?;
    Ok((
        Model::Prime {
            p,
            a: short_a,
            b: short_b,
        },
        CoordinateMap::Edwards {
            a_over_three,
            inv_b,
            inv_scale,
        },
    ))
}

fn parse_model(model_json: &str) -> Result<(Model, CoordinateMap), String> {
    let m: Value = serde_json::from_str(model_json).map_err(|e| format!("model JSON: {e}"))?;
    let num = |key: &str| {
        str_field(&m, key).and_then(|s| parse_uint(s).ok_or_else(|| format!("bad {key}: {s}")))
    };
    match str_field(&m, "form")? {
        "y^2+xy=x^3+a*x^2+b" => {
            let modulus = num("modulus")?;
            let n = (modulus.bits() - 1) as u32;
            let field = str_field(&m, "field")?;
            if !field.starts_with(&format!("f2m-{n}-")) {
                return Err(format!("field {field} is not of degree {n}"));
            }
            Ok((
                Model::Binary {
                    n,
                    irr: irreducible(&modulus),
                    a: F2mElement::from_biguint(&num("a")?, n),
                    b: F2mElement::from_biguint(&num("b")?, n),
                    modulus,
                },
                CoordinateMap::Identity,
            ))
        }
        "y^2=x^3+a*x+b" if m["a"].is_array() => {
            let p = num("p")?;
            let k: u32 = str_field(&m, "k")?
                .parse()
                .map_err(|_| "bad extension degree k")?;
            if p <= BigUint::one() || k < 2 {
                return Err("invalid extension field size".into());
            }
            let field_id = str_field(&m, "field")?.to_string();
            if !field_id.starts_with(&format!("fpk-{p}-{k}-")) {
                return Err(format!("field {field_id} is not GF({p}^{k})"));
            }
            let coeffs = |key: &str| -> Result<Vec<BigUint>, String> {
                m[key]
                    .as_array()
                    .ok_or_else(|| format!("missing coefficient array {key}"))?
                    .iter()
                    .map(|value| {
                        value
                            .as_str()
                            .and_then(parse_uint)
                            .ok_or_else(|| format!("bad coefficient in {key}"))
                    })
                    .collect()
            };
            let modulus = coeffs("modulus")?;
            let a = coeffs("a")?;
            let b = coeffs("b")?;
            if [modulus.len(), a.len(), b.len()]
                .iter()
                .any(|&len| len != k as usize)
                || modulus
                    .iter()
                    .chain(&a)
                    .chain(&b)
                    .any(|coefficient| coefficient >= &p)
            {
                return Err("extension coefficients are not canonical k-vectors".into());
            }
            let q = p.pow(k);
            Ok((
                Model::Extension {
                    p,
                    k,
                    q,
                    field_id,
                    modulus,
                    a,
                    b,
                },
                CoordinateMap::Identity,
            ))
        }
        "y^2=x^3+a*x+b" => Ok((
            Model::Prime {
                p: num("p")?,
                a: num("a")?,
                b: num("b")?,
            },
            CoordinateMap::Identity,
        )),
        "B*y^2=x^3+A*x^2+x" => {
            let p = num("p")?;
            let (short_a, short_b, a_over_three, inv_b) =
                montgomery_short(&p, &num("A")?, &num("B")?)?;
            Ok((
                Model::Prime {
                    p,
                    a: short_a,
                    b: short_b,
                },
                CoordinateMap::Montgomery {
                    a_over_three,
                    inv_b,
                },
            ))
        }
        "a*x^2+y^2=1+d*x^2*y^2" => edwards_short(num("p")?, num("a")?, num("d")?, BigUint::one()),
        "x^2+y^2=c^2*(1+d*x^2*y^2)" => {
            edwards_short(num("p")?, BigUint::one(), num("d")?, num("c")?)
        }
        other => Err(format!("unknown form {other}")),
    }
}

fn montgomery_coordinates(
    u: &BigUint,
    v: &BigUint,
    p: &BigUint,
    a_over_three: &BigUint,
    inv_b: &BigUint,
) -> (BigUint, BigUint) {
    ((u + a_over_three) * inv_b % p, v * inv_b % p)
}

fn short_coordinates(
    x: BigUint,
    y: BigUint,
    p: &BigUint,
    map: &CoordinateMap,
) -> Option<(BigUint, BigUint)> {
    match map {
        CoordinateMap::Identity => Some((x, y)),
        CoordinateMap::Montgomery {
            a_over_three,
            inv_b,
        } => Some(montgomery_coordinates(&x, &y, p, a_over_three, inv_b)),
        CoordinateMap::Edwards {
            a_over_three,
            inv_b,
            inv_scale,
        } => {
            let x = x * inv_scale % p;
            let y = y * inv_scale % p;
            let inv_one_minus_y = inverse_mod(&subtract_mod(&BigUint::one(), &y, p), p)?;
            let u = (&y + BigUint::one()) * inv_one_minus_y % p;
            let v = &u * inverse_mod(&x, p)? % p;
            Some(montgomery_coordinates(&u, &v, p, a_over_three, inv_b))
        }
    }
}

/// A representation's generator, mapped from its registered model to the
/// short Weierstrass model when its prime-field form requires conversion.
/// Big integers in representations are bare JSON numbers, which this
/// parser would round, so only string-typed and small fields are read.
fn parse_representation(rep: &Value, model: &Model, map: &CoordinateMap) -> Option<Representation> {
    let field = &rep["field"];
    let curve = &rep["curve"];
    match model {
        Model::Binary { modulus, .. } => {
            let exps = field["modulus_exponents"].as_array()?;
            let mut m = BigUint::zero();
            for e in exps {
                m.set_bit(e.as_u64()?, true);
            }
            if m != *modulus || !field["representation"].as_str()?.starts_with("polynomial") {
                return None;
            }
        }
        Model::Prime { .. } => {
            if field["degree"].as_u64()? != 1 {
                return None;
            }
        }
        Model::Other { .. } => return None,
        Model::Extension { p, k, modulus, .. } => {
            if field["degree"].as_u64()? != u64::from(*k)
                || value_uint(&field["characteristic"])? != *p
                || !field["representation"].as_str()?.starts_with("polynomial")
            {
                return None;
            }
            let field_modulus = field["modulus"].as_array()?;
            if field_modulus.len() != *k as usize + 1
                || value_uint(field_modulus.last()?)? != BigUint::one()
                || !field_modulus[..*k as usize]
                    .iter()
                    .zip(modulus)
                    .all(|(value, coefficient)| value_uint(value).as_ref() == Some(coefficient))
            {
                return None;
            }
        }
    }
    let gen = curve["generator"].as_array()?;
    let cofactor = value_uint(&curve["cofactor"])?;
    let generator = match model {
        Model::Binary { .. } => GeneratorCoords::Affine(
            parse_uint(gen.first()?.as_str()?)?,
            parse_uint(gen.get(1)?.as_str()?)?,
        ),
        Model::Prime { p, .. } => {
            let x = parse_uint(gen.first()?.as_str()?)?;
            let y = parse_uint(gen.get(1)?.as_str()?)?;
            let (x, y) = short_coordinates(x, y, p, map)?;
            GeneratorCoords::Affine(x, y)
        }
        Model::Extension { p, k, .. } => GeneratorCoords::Extension(
            extension_coordinates(gen.first()?, *k, p)?,
            extension_coordinates(gen.get(1)?, *k, p)?,
        ),
        Model::Other { .. } => return None,
    };
    Some(Representation {
        generator,
        subgroup_order: parse_uint(curve["subgroup_order"].as_str()?)?,
        cofactor,
    })
}

/// Parse `docs/curves/registry.json`.  The trace is read from the ICV1
/// string and checked against the recorded order.
pub fn read_registry(text: &str) -> Result<Vec<RegistryCurve>, String> {
    let doc: Value = serde_json::from_str(text).map_err(|e| format!("registry JSON: {e}"))?;
    let curves = doc["curves"].as_array().ok_or("registry has no curves")?;
    let mut out = Vec::with_capacity(curves.len());
    for c in curves {
        let slug = str_field(c, "slug")?.to_string();
        let ctx = |e: String| format!("{slug}: {e}");
        let (model, coord_map) = parse_model(str_field(c, "model_json")?).map_err(ctx)?;
        let icv1: Vec<&str> = str_field(c, "icv1")?.split(':').collect();
        let trace: BigInt = icv1
            .get(2)
            .and_then(|t| t.parse().ok())
            .ok_or_else(|| ctx("ICV1 string has no trace".into()))?;
        let order = parse_uint(str_field(c, "order")?).ok_or_else(|| ctx("bad order".into()))?;
        if BigInt::from(model.q()) + 1 - &trace != BigInt::from(order.clone()) {
            return Err(ctx(format!("order {order} ≠ q + 1 − t for t = {trace}")));
        }
        if &trace * &trace > BigInt::from(model.q()) * 4u8 {
            return Err(ctx(format!("|t| = |{trace}| exceeds the Hasse bound 2√q")));
        }
        let end = match str_field(c, "end")? {
            "unk" => None,
            d => Some(d.parse().map_err(|_| ctx(format!("bad end {d}")))?),
        };
        let representations = c["representations"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|r| parse_representation(r, &model, &coord_map))
            .collect();
        out.push(RegistryCurve {
            slug,
            family: str_field(c, "family")?.to_string(),
            model,
            trace,
            order,
            end,
            representations,
            levels: BTreeMap::new(),
            source: Source::registry(),
        });
    }
    Ok(out)
}

// ── The record ────────────────────────────────────────────────────────

/// Integers too large for JSON numbers are written as decimal strings.
pub type Factors = Vec<(String, u32)>;

fn factors(f: &[(BigUint, u32)]) -> Factors {
    f.iter().map(|(p, e)| (p.to_string(), *e)).collect()
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FieldInfo {
    /// `binary` or `prime`.
    pub kind: String,
    /// `n` of `GF(2^n)`; `1` for `GF(p)`.
    pub degree: u32,
    /// `n` for `GF(2^n)`, the bit length of `p` for `GF(p)`.
    pub size_bits: u64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Checked {
    pub status: Status,
    pub method: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Frobenius {
    /// `t² − 4q`.
    pub disc: String,
    pub disc_primes: Factors,
    /// Composites of `|t² − 4q|` the budget did not split.
    pub disc_unfactored: Factors,
    pub disc_status: Status,
    /// `d_K`, the fundamental discriminant of `Q(π)`.
    pub cm_disc: Option<String>,
    /// `v = [O_K : Z[π]]`; a lower bound when `status` is `bounded`.
    pub conductor: String,
    pub conductor_bits: u64,
    pub conductor_primes: Factors,
    pub conductor_unfactored: Factors,
    /// `log v / log √|Δ|`: 0 when `Z[π]` is maximal, near 1 when `d_K` is
    /// tiny beside `Δ`.  Absent when `v` is only bounded.
    pub conductor_fraction: Option<f64>,
    /// Of `cm_disc`, `conductor` and `conductor_fraction`.
    pub status: Status,
    /// `h(d_K)`.
    pub class_number: Option<u64>,
    pub class_number_status: Status,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SmallPrime {
    pub ell: u64,
    /// `v_ℓ(v)`: the depth of the `ℓ`-isogeny volcano.
    pub depth: u32,
    /// How `ℓ` splits in `Q(π)`: `split`, `inert` or `ramified`.
    pub splitting: String,
    /// Rational `ℓ`-isogenies: `1 + (d_K/ℓ)` when the depth is 0;
    /// absent when it depends on the volcano level or `ℓ` is the
    /// characteristic.
    pub rational_isogenies: Option<u32>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SubfieldRecord {
    pub j_field_degree: Option<u32>,
    pub definition_degree: Option<u32>,
    /// Trace of the descended model over `GF(2^k)`; `|t_k|` when
    /// `base_trace_sign_free`.
    pub base_trace: Option<String>,
    pub base_trace_sign_free: bool,
    pub status: Status,
    pub method: String,
}

/// The volcano level of `End(E)` at one prime `ℓ | v`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LevelRecord {
    pub ell: String,
    /// `v_ℓ(v)`: the depth of the `ℓ`-volcano.
    pub depth: u32,
    /// `v_ℓ([O_K : End(E)])`: 0 on the crater, `depth` on the floor.
    pub level: Option<u32>,
    /// The level is at most this.
    pub level_at_most: u32,
    pub method: String,
}

/// What is known of `End(E)`, an order of `Q(π)` between `Z[π]` and
/// `O_K`.  For an ordinary curve every geometric endomorphism is defined
/// over the working field, so this is also the ring over that field.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EndomorphismRecord {
    /// `f₀` with the conductor `[O_K : End(E)]` dividing it: the conductor
    /// `v` of `Z[π]`, sharpened by every other order known to lie in
    /// `End(E)` and by the volcano levels below.
    pub conductor_divides: Option<String>,
    /// `[O_K : End(E)]`, when every level is resolved.
    pub conductor: Option<String>,
    /// One entry per prime of `v`.
    pub levels: Vec<LevelRecord>,
    /// `true` when `End(E) = O_K` is established; `false` when a
    /// certificate fixes a smaller order; absent when unresolved.
    pub maximal: Option<bool>,
    /// [`Status::Bounded`] when only `conductor_divides` is known.
    pub status: Status,
    /// The orders known to lie in `End(E)`.
    pub contains: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Subgroup {
    pub order: Option<String>,
    pub bits: Option<u64>,
    pub cofactor: Option<String>,
    pub status: Status,
    pub source: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Embedding {
    /// `ord_r(q)`.
    pub degree: Option<String>,
    /// When only `degree > lower_bound` is known.
    pub lower_bound: Option<u64>,
    /// Bit length of `(r − 1)/k`.
    pub complement_bits: Option<u64>,
    pub status: Status,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Twist {
    /// `q + 1 + t`.
    pub order: String,
    pub largest_prime_bits: Option<u64>,
    pub cofactor: Option<String>,
    pub status: Status,
}

/// Everything derived about one curve.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurveTraits {
    pub slug: String,
    pub family: String,
    pub source: Source,
    pub field: FieldInfo,
    pub trace: String,
    pub order: String,
    pub order_check: Checked,
    pub ordinary: bool,
    /// `#E = q`.
    pub anomalous: bool,
    /// `t / 2√q ∈ [−1, 1]`.
    pub trace_ratio: f64,
    pub frobenius: Frobenius,
    pub endomorphism: EndomorphismRecord,
    pub small_primes: Vec<SmallPrime>,
    pub subfield: SubfieldRecord,
    pub subgroup: Subgroup,
    pub embedding: Embedding,
    pub twist: Twist,
    /// The DiSSECT families ([`dissect`]); absent for models whose
    /// arithmetic this module does not have.
    pub dissect: Option<dissect::DissectRecord>,
    /// The grouping keys of [`KEYS`], precomputed.
    pub keys: BTreeMap<String, String>,
}

// ── Computing it ──────────────────────────────────────────────────────

fn frobenius(c: &RegistryCurve, fx: &Factoring) -> Result<(Frobenius, BigInt), String> {
    let q = BigInt::from(c.model.q());
    let delta: BigInt = &c.trace * &c.trace - BigInt::from(4u8) * &q;
    let fac = fx.factor(delta.magnitude());
    let split = arith::split_discriminant(&delta, &fac);
    if let (Some(end), Some(dk)) = (&c.end, &split.cm_disc) {
        // A certified End(E) discriminant is d_K·f² for some f | v.
        let (quot, rem) = end.div_rem(dk);
        let f = quot.sqrt();
        if !rem.is_zero() || &f * &f != quot || !(&split.conductor % f.magnitude()).is_zero() {
            return Err(format!(
                "registry certifies End(E) of discriminant {end}, not an order of Q(√{dk})"
            ));
        }
    }
    let fraction = (split.cm_disc.is_some()).then(|| {
        let half_log = arith::log2_abs(&delta) / 2.0;
        let v = BigInt::from(split.conductor.clone());
        let log_v = if v.is_one() { 0.0 } else { arith::log2_abs(&v) };
        arith::round_to(log_v / half_log, 6)
    });
    let (class_number, class_number_status) = match &split.cm_disc {
        Some(dk) if dk.magnitude() <= &BigUint::from(CLASS_NUMBER_MAX_DISC) => (
            Some(crate::isogeny::class_group::class_number(dk)),
            split.status,
        ),
        Some(_) => (None, Status::NotEvaluated),
        None => (None, Status::Unknown),
    };
    Ok((
        Frobenius {
            disc: delta.to_string(),
            disc_primes: factors(&fac.primes),
            disc_unfactored: factors(&fac.composites),
            disc_status: fac.status(),
            cm_disc: split.cm_disc.as_ref().map(BigInt::to_string),
            conductor: split.conductor.to_string(),
            conductor_bits: split.conductor.bits(),
            conductor_primes: factors(&split.conductor_primes),
            conductor_unfactored: factors(&split.conductor_composites),
            conductor_fraction: fraction,
            status: split.status,
            class_number,
            class_number_status,
        },
        delta,
    ))
}

fn small_primes(delta: &BigInt, characteristic: &BigUint) -> Vec<SmallPrime> {
    SMALL_PRIMES
        .iter()
        .map(|&ell| {
            let loc = arith::ell_local(delta, ell);
            let is_char = *characteristic == BigUint::from(ell);
            SmallPrime {
                ell,
                depth: loc.depth,
                splitting: match loc.splitting {
                    1 => "split",
                    -1 => "inert",
                    _ => "ramified",
                }
                .to_string(),
                rational_isogenies: (loc.depth == 0 && !is_char)
                    .then(|| (1 + loc.splitting) as u32),
            }
        })
        .collect()
}

fn subfield_record(c: &RegistryCurve) -> Result<SubfieldRecord, String> {
    if matches!(&c.model, Model::Extension { .. }) {
        return Ok(SubfieldRecord {
            j_field_degree: None,
            definition_degree: None,
            base_trace: None,
            base_trace_sign_free: false,
            status: Status::NotEvaluated,
            method: "extension-field subfield analysis is not implemented".into(),
        });
    }
    let Model::Binary { n, irr, a, b, .. } = &c.model else {
        return Ok(SubfieldRecord {
            j_field_degree: None,
            definition_degree: None,
            base_trace: None,
            base_trace_sign_free: false,
            status: Status::NotApplicable,
            method: "a prime field has no proper subfield".into(),
        });
    };
    let s = subfield::analyse(*n, irr, a, b, &c.trace)?;
    Ok(SubfieldRecord {
        j_field_degree: Some(s.j_field_degree),
        definition_degree: Some(s.definition_degree),
        base_trace: s.base_trace.as_ref().map(BigInt::to_string),
        base_trace_sign_free: s.sign_free,
        status: s.status,
        method: s.method.into(),
    })
}

fn subgroup(c: &RegistryCurve, fx: &Factoring) -> Result<(Subgroup, Option<BigUint>), String> {
    if let Some(rep) = c.representations.first() {
        let r = rep.subgroup_order.clone();
        if &rep.cofactor * &r != c.order || !arith::is_prime(&r) {
            return Err(format!(
                "representation records subgroup order {r} and cofactor {}, not a prime-order \
                 subgroup of #E = {}",
                rep.cofactor, c.order
            ));
        }
        return Ok((
            Subgroup {
                order: Some(r.to_string()),
                bits: Some(r.bits()),
                cofactor: Some(rep.cofactor.to_string()),
                status: arith::prime_status(&r),
                source: "registry representation".into(),
            },
            Some(r),
        ));
    }
    let fac = fx.factor(&c.order);
    Ok(match fac.largest_prime() {
        Some(r) => {
            let r = r.clone();
            (
                Subgroup {
                    order: Some(r.to_string()),
                    bits: Some(r.bits()),
                    cofactor: Some((&c.order / &r).to_string()),
                    status: fac.status(),
                    source: "largest prime factor of the order".into(),
                },
                Some(r),
            )
        }
        None => (
            Subgroup {
                order: None,
                bits: None,
                cofactor: None,
                status: Status::Unknown,
                source: "the order did not factor within the budget".into(),
            },
            None,
        ),
    })
}

fn embedding(q: &BigUint, r: Option<&BigUint>, fx: &Factoring) -> Embedding {
    let unknown = |status| Embedding {
        degree: None,
        lower_bound: None,
        complement_bits: None,
        status,
    };
    let Some(r) = r else {
        return unknown(Status::Unknown);
    };
    if r.is_one() || (q % r).is_zero() {
        return unknown(Status::NotApplicable);
    }
    let rm1 = r - 1u8;
    let fac: Factored = fx.factor(&rm1);
    let r_status = arith::prime_status(r);
    let exact = |k: BigUint, status| Embedding {
        complement_bits: Some((&rm1 / &k).bits()),
        degree: Some(k.to_string()),
        lower_bound: None,
        status,
    };
    if fac.complete() {
        let k = arith::multiplicative_order(q, r, &fac);
        return exact(k, arith::weaker(r_status, fac.status()));
    }
    match arith::order_at_most(q, r, EMBEDDING_SEARCH) {
        Some(k) => exact(BigUint::from(k), r_status),
        None => Embedding {
            degree: None,
            lower_bound: Some(EMBEDDING_SEARCH),
            complement_bits: None,
            status: Status::Bounded,
        },
    }
}

fn twist(q: &BigUint, trace: &BigInt, fx: &Factoring) -> Twist {
    let order = (BigInt::from(q.clone()) + BigInt::one() + trace)
        .to_biguint()
        .expect("Hasse: q + 1 + t > 0");
    let fac = fx.factor(&order);
    let (bits, cofactor) = match fac.largest_prime() {
        Some(p) => (Some(p.bits()), Some((&order / p).to_string())),
        None => (None, None),
    };
    Twist {
        order: order.to_string(),
        largest_prime_bits: bits,
        cofactor,
        status: if fac.complete() {
            fac.status()
        } else {
            Status::Unknown
        },
    }
}

/// Bound the conductor of `End(E)` by every order known to lie in it.
/// `Z[π] ⊂ End(E)` gives `f | v`.  A field of definition `GF(2^k)` puts
/// its Frobenius `π_k` in `End(E)` (`π = π_k^{n/k}`, and `Q(π_k) = Q(π)`),
/// so `f | v_k` with `t_k² − 4·2^k = v_k²·d_K`.  `j = 0` and `j = 1728`
/// on a prime-field model put `Z[ζ₃]` and `Z[i]`, both maximal, in
/// `End(E)`.  A registry certificate fixes `End(E)` outright.  `Err` when
/// one of these contradicts `d_K`.
fn endomorphism(
    c: &RegistryCurve,
    frob: &Frobenius,
    sub: &SubfieldRecord,
    ordinary: bool,
) -> Result<EndomorphismRecord, String> {
    let record = |f0: Option<String>, maximal, status, contains: Vec<String>| EndomorphismRecord {
        conductor_divides: f0,
        conductor: None,
        levels: Vec::new(),
        maximal,
        status,
        contains,
    };
    if !ordinary {
        return Ok(record(None, None, Status::NotApplicable, Vec::new()));
    }
    let (Some(dk), Status::Proved | Status::Probable) = (&frob.cm_disc, frob.status) else {
        return Ok(record(None, None, Status::Unknown, Vec::new()));
    };
    let dk: BigInt = dk.parse().expect("written by frobenius()");
    let mut f0: BigUint = frob.conductor.parse().expect("written by frobenius()");
    let mut contains = vec!["Z[π]".to_string()];
    if let (Model::Binary { n, .. }, Some(k), Some(tk)) =
        (&c.model, sub.definition_degree, &sub.base_trace)
    {
        if k < *n {
            let tk: BigInt = tk.parse().expect("written by subfield_record()");
            let disc_k = (BigInt::one() << (k + 2)) - &tk * &tk;
            let (sq, rem) = disc_k.div_rem(&-&dk);
            let vk = sq.sqrt();
            if !rem.is_zero() || &vk * &vk != sq {
                return Err(format!(
                    "t_k² − 4·2^{k} = {} is not v²·d_K for d_K = {dk}",
                    -disc_k
                ));
            }
            f0 = f0.gcd(vk.magnitude());
            contains.push(format!("Z[π_{k}] (Frobenius of GF(2^{k}))"));
        }
    }
    if let Model::Prime { a, b, .. } = &c.model {
        for (special, d, ring) in [(a, -3, "Z[ζ₃] (j = 0)"), (b, -4, "Z[i] (j = 1728)")] {
            if special.is_zero() {
                if dk != BigInt::from(d) {
                    return Err(format!("{ring} ⊂ End(E), but d_K = {dk}"));
                }
                f0 = BigUint::one();
                contains.push(ring.to_string());
            }
        }
    }
    let status = frob.status;
    let v_primes: Vec<(BigUint, u32)> = frob
        .conductor_primes
        .iter()
        .map(|(p, e)| (p.parse().expect("written by frobenius()"), *e))
        .collect();
    let val = |n: &BigUint, ell: &BigUint| {
        let (mut n, mut e) = (n.clone(), 0u32);
        while !n.is_zero() && (&n % ell).is_zero() {
            n /= ell;
            e += 1;
        }
        e
    };
    for (&ell, &level) in &c.levels {
        if level > 0 && !v_primes.iter().any(|(p, _)| *p == BigUint::from(ell)) {
            return Err(format!(
                "the source certifies level {level} at ℓ = {ell}, which does not divide v"
            ));
        }
    }
    if let Some(end) = &c.end {
        let f = (end / &dk).sqrt();
        if !(&f0 % f.magnitude()).is_zero() {
            return Err(format!(
                "certified End(E) of conductor {f}, which does not divide {f0}"
            ));
        }
        contains.push(format!("End(E) of discriminant {end} (source certificate)"));
        let f = f.magnitude().clone();
        let mut levels = Vec::new();
        for (ell, depth) in &v_primes {
            let level = val(&f, ell);
            if let Some(&l) = ell.to_u64().and_then(|e| c.levels.get(&e)) {
                if l != level {
                    return Err(format!(
                        "level {l} at ℓ = {ell} contradicts End(E) of conductor {f}"
                    ));
                }
            }
            levels.push(LevelRecord {
                ell: ell.to_string(),
                depth: *depth,
                level: Some(level),
                level_at_most: level,
                method: "certified End(E)".into(),
            });
        }
        return Ok(EndomorphismRecord {
            conductor_divides: Some(f.to_string()),
            conductor: Some(f.to_string()),
            levels,
            maximal: Some(f.is_one()),
            status,
            contains,
        });
    }
    // Resolve each prime of v: a certified level, else the scalar test
    // (torsion.rs): Frobenius scalar on E[ℓ] iff ℓ | [End(E) : Z[π]].
    let mut levels = Vec::new();
    let mut bound = BigUint::one();
    let mut exact = BigUint::one();
    let mut resolved = true;
    for (ell, depth) in &v_primes {
        let e = val(&f0, ell);
        let small = ell.to_u64();
        let (level, at_most, method): (Option<u32>, u32, String) = if e == 0 {
            (Some(0), 0, "f₀ is prime to ℓ".into())
        } else if let Some(&l) = small.and_then(|s| c.levels.get(&s)) {
            if l > e {
                return Err(format!(
                    "the source certifies level {l} at ℓ = {ell}, above v_ℓ(f₀) = {e}"
                ));
            }
            // Where the scalar test runs, it must agree.
            let tested = small.and_then(|s| torsion::frobenius_is_scalar(&c.model, &c.trace, s));
            if tested.is_some_and(|scalar| scalar != (l < *depth)) {
                return Err(format!(
                    "the source certifies level {l} of {depth} at ℓ = {ell}; the scalar test on \
                     E[ℓ] says otherwise"
                ));
            }
            let how = if tested.is_some() {
                ", agreeing with the scalar test"
            } else {
                ""
            };
            (Some(l), l, format!("certified by {}{how}", c.source.kind))
        } else {
            match small.and_then(|l| torsion::frobenius_is_scalar(&c.model, &c.trace, l)) {
                Some(true) => {
                    let m = e.min(depth - 1);
                    (
                        (m == 0).then_some(0),
                        m,
                        "Frobenius is a scalar on E[ℓ]: above the floor".into(),
                    )
                }
                Some(false) => {
                    if *depth > e {
                        return Err(format!(
                            "Frobenius is not a scalar on E[{ell}], so End(E) is on the floor, \
                             but v_ℓ(f₀) = {e} < depth {depth}"
                        ));
                    }
                    (
                        Some(*depth),
                        *depth,
                        "Frobenius is not a scalar on E[ℓ]: the floor".into(),
                    )
                }
                None => (
                    None,
                    e,
                    "not tested (ℓ too large for the division-polynomial test)".into(),
                ),
            }
        };
        bound *= ell.pow(at_most);
        match level {
            Some(l) => exact *= ell.pow(l),
            None => resolved = false,
        }
        levels.push(LevelRecord {
            ell: ell.to_string(),
            depth: *depth,
            level,
            level_at_most: at_most,
            method,
        });
    }
    if levels
        .iter()
        .any(|l| l.level.is_some() && l.method.starts_with("Frobenius"))
    {
        contains.push("the levels the scalar test on E[ℓ] resolves".into());
    }
    Ok(if resolved {
        EndomorphismRecord {
            conductor_divides: Some(exact.to_string()),
            conductor: Some(exact.to_string()),
            maximal: Some(exact.is_one()),
            levels,
            status,
            contains,
        }
    } else {
        EndomorphismRecord {
            conductor_divides: Some(bound.to_string()),
            conductor: None,
            maximal: None,
            levels,
            status: Status::Bounded,
            contains,
        }
    })
}

/// Derive every trait of one curve.  `Err` when the registry entry
/// contradicts itself (an order the model refutes, a certified `End(E)`
/// outside `Q(π)`).
pub fn compute(c: &RegistryCurve, fx: &Factoring) -> Result<CurveTraits, String> {
    let ctx = |e: String| format!("{}: {e}", c.slug);
    let q = c.model.q();
    let characteristic = match &c.model {
        Model::Binary { .. } => BigUint::from(2u8),
        Model::Prime { p, .. } | Model::Other { p, .. } | Model::Extension { p, .. } => p.clone(),
    };
    let (status, method) =
        certify::check_order(&c.model, &c.order, &c.representations).map_err(ctx)?;
    let (frob, delta) = frobenius(c, fx).map_err(ctx)?;
    let ordinary = !(&c.trace % BigInt::from(characteristic.clone())).is_zero();
    let trace_ratio = {
        let r = arith::log2_abs(&c.trace) - 1.0 - arith::log2_abs(&BigInt::from(q.clone())) / 2.0;
        let mag = if c.trace.is_zero() { 0.0 } else { r.exp2() };
        arith::round_to(if c.trace.is_negative() { -mag } else { mag }, 6)
    };
    let (subgroup, r) = subgroup(c, fx).map_err(ctx)?;
    let sub = subfield_record(c).map_err(ctx)?;
    let endo = endomorphism(c, &frob, &sub, ordinary).map_err(ctx)?;
    let small = small_primes(&delta, &characteristic);
    let dissect = match &c.model {
        Model::Prime { .. } | Model::Binary { .. } if ordinary => {
            let scalar = |ell: u64| -> Option<bool> {
                let depth = small.iter().find(|s| s.ell == ell).map(|s| s.depth)?;
                if depth == 0 {
                    return Some(false);
                }
                let lv = endo.levels.iter().find(|l| l.ell == ell.to_string())?;
                match lv.level {
                    Some(l) => Some(l < depth),
                    None if lv.level_at_most < depth => Some(true),
                    None => None,
                }
            };
            Some(dissect::compute(
                fx,
                &c.model,
                &characteristic,
                &c.trace,
                &c.order,
                r.as_ref(),
                status,
                &scalar,
            ))
        }
        _ => None,
    };
    let mut out = CurveTraits {
        slug: c.slug.clone(),
        family: c.family.clone(),
        source: c.source.clone(),
        field: FieldInfo {
            kind: match c.model {
                Model::Binary { .. } => "binary",
                Model::Prime { .. } | Model::Other { .. } => "prime",
                Model::Extension { .. } => "extension",
            }
            .into(),
            degree: match &c.model {
                Model::Binary { n, .. } => *n,
                Model::Prime { .. } => 1,
                Model::Other { degree, .. } => *degree,
                Model::Extension { k, .. } => *k,
            },
            size_bits: c.model.size_bits(),
        },
        trace: c.trace.to_string(),
        order: c.order.to_string(),
        order_check: Checked {
            status,
            method: method.into(),
        },
        ordinary,
        anomalous: c.trace.is_one(),
        trace_ratio,
        small_primes: small,
        frobenius: frob,
        endomorphism: endo,
        subfield: sub,
        embedding: embedding(&q, r.as_ref(), fx),
        twist: twist(&q, &c.trace, fx),
        subgroup,
        dissect,
        keys: BTreeMap::new(),
    };
    out.keys = keys(&out);
    Ok(out)
}

// ── Grouping keys ─────────────────────────────────────────────────────

/// Every grouping key and what it means.  Each is discrete and does not
/// grow with the field, so equal keys can join curves of any size.
pub const KEYS: &[(&str, &str)] = &[
    ("char", "field characteristic: 2 or p"),
    ("family", "the registry's family label"),
    ("ordinary", "yes when p ∤ t"),
    (
        "cm",
        "d_K of Q(π) when |d_K| ≤ 10^6; 'large' above; 'unknown' when Δ did not factor",
    ),
    (
        "descent",
        "binary: k<deg>|t|=<|t_k|> for the smallest subfield GF(2^k) the curve is defined over (Koblitz: k1|t|=1); 'none'",
    ),
    (
        "descent_signed",
        "descent with the sign of t_k, or ± when both twists over GF(2^k) descend",
    ),
    ("jfield", "k<deg> when j lies in a proper subfield GF(2^k); 'full'"),
    ("cofactor", "#E / r"),
    ("twist_cofactor", "#E' / (largest prime of #E') for the quadratic twist"),
    (
        "embedding",
        "k=<deg> for an embedding degree ≤ 20; 'large' above; 'unknown'",
    ),
    (
        "conductor",
        "'1' when Z[π] is maximal; 'smooth' when every prime of v is below 2^16; 'rough'; 'unknown'",
    ),
    ("split", "how 3, 5, 7, 11, 13 split in Q(π): s split, i inert, r ramified"),
    ("depth", "v_ℓ(v) for ℓ = 2, 3, 5, 7: the ℓ-volcano depths"),
    (
        "end_order",
        "'maximal' when End(E) = O_K is established; 'non-maximal' when a smaller order is established; 'unresolved'; 'unknown'",
    ),
    (
        "isogenies",
        "rational ℓ-isogenies over GF(q) for ℓ = 2, 3, 5, 7 ('?' unknown, '-' the characteristic)",
    ),
    ("source", "registry, isogeny-walk or ic-crosswalk"),
];

/// The keys a curve's default signature is made of.
pub const DEFAULT_SIGNATURE: &[&str] = &["char", "ordinary", "cm", "descent", "cofactor"];

fn keys(t: &CurveTraits) -> BTreeMap<String, String> {
    let mut k = BTreeMap::new();
    let mut put = |key: &str, value: String| {
        k.insert(key.to_string(), value);
    };
    put(
        "char",
        if t.field.kind == "binary" { "2" } else { "p" }.into(),
    );
    put("family", t.family.clone());
    put("source", t.source.kind.clone());
    put(
        "isogenies",
        [2u64, 3, 5, 7]
            .iter()
            .map(|&l| {
                let n = match &t.dissect {
                    Some(d) => match d.torsion.iter().find(|r| r.ell == l) {
                        Some(r) => r
                            .rational_isogenies
                            .map_or_else(|| "?".to_string(), |n| n.to_string()),
                        None => "-".into(),
                    },
                    None => "?".into(),
                };
                format!("{l}:{n}")
            })
            .collect::<Vec<_>>()
            .join(","),
    );
    put("ordinary", if t.ordinary { "yes" } else { "no" }.into());
    put(
        "end_order",
        match (t.endomorphism.maximal, t.endomorphism.status) {
            (Some(true), _) => "maximal",
            (Some(false), _) => "non-maximal",
            (None, Status::Bounded) => "unresolved",
            (None, Status::NotApplicable) => "n/a",
            _ => "unknown",
        }
        .into(),
    );
    put(
        "cm",
        match &t.frobenius.cm_disc {
            Some(d)
                if d.parse::<i64>()
                    .is_ok_and(|v| v.unsigned_abs() <= SMALL_CM_DISC) =>
            {
                d.clone()
            }
            Some(_) => "large".into(),
            None => "unknown".into(),
        },
    );
    let s = &t.subfield;
    let (descent, signed) = match (s.definition_degree, &s.base_trace) {
        (Some(k), Some(bt)) if Some(k) != Some(t.field.degree) => {
            let abs = bt.trim_start_matches('-');
            let signed = if s.base_trace_sign_free {
                format!("k{k}|t|=±{abs}")
            } else {
                format!("k{k}|t|={bt}")
            };
            (format!("k{k}|t|={abs}"), signed)
        }
        (Some(k), None) if k != t.field.degree => (format!("k{k}|t|=?"), format!("k{k}|t|=?")),
        _ => ("none".into(), "none".into()),
    };
    put("descent", descent);
    put("descent_signed", signed);
    put(
        "jfield",
        match s.j_field_degree {
            Some(k) if k != t.field.degree => format!("k{k}"),
            _ => "full".into(),
        },
    );
    put(
        "cofactor",
        t.subgroup
            .cofactor
            .clone()
            .unwrap_or_else(|| "unknown".into()),
    );
    put(
        "twist_cofactor",
        t.twist.cofactor.clone().unwrap_or_else(|| "unknown".into()),
    );
    put(
        "embedding",
        match (&t.embedding.degree, t.embedding.lower_bound) {
            (Some(d), _) => match d.parse::<u64>() {
                Ok(v) if v <= SMALL_EMBEDDING_DEGREE => format!("k={v}"),
                _ => "large".into(),
            },
            (None, Some(lb)) if lb >= SMALL_EMBEDDING_DEGREE => "large".into(),
            _ => "unknown".into(),
        },
    );
    let f = &t.frobenius;
    put(
        "conductor",
        if f.status == Status::Bounded {
            "unknown".into()
        } else if f.conductor == "1" {
            "1".into()
        } else if f.conductor_unfactored.is_empty()
            && f.conductor_primes
                .iter()
                .all(|(p, _)| p.parse::<u64>().is_ok_and(|p| p < 1 << 16))
        {
            "smooth".into()
        } else {
            "rough".into()
        },
    );
    let local = |ell: u64| t.small_primes.iter().find(|s| s.ell == ell);
    put(
        "split",
        [3u64, 5, 7, 11, 13]
            .iter()
            .filter_map(|&l| local(l))
            .map(|s| format!("{}{}", s.ell, &s.splitting[..1]))
            .collect::<Vec<_>>()
            .join(""),
    );
    put(
        "depth",
        [2u64, 3, 5, 7]
            .iter()
            .filter_map(|&l| local(l))
            .map(|s| format!("{}:{}", s.ell, s.depth))
            .collect::<Vec<_>>()
            .join(","),
    );
    k
}

/// The curve's signature over `keys`: `key=value` joined by `;`.
pub fn signature(t: &CurveTraits, keys: &[&str]) -> String {
    keys.iter()
        .map(|k| format!("{k}={}", t.keys.get(*k).map_or("?", String::as_str)))
        .collect::<Vec<_>>()
        .join(";")
}

// ── The trait file ────────────────────────────────────────────────────

pub const SCHEMA: &str = "curve-traits/2";

/// One file a build read, by path and SHA-256.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct InputDigest {
    pub path: String,
    pub sha256: String,
}

/// A file handed to [`TraitFile::build_from`].
pub struct Input {
    pub path: String,
    pub text: String,
}

impl Input {
    fn digest(&self) -> InputDigest {
        let d = crate::hash::sha256::sha256(self.text.as_bytes());
        InputDigest {
            path: self.path.clone(),
            sha256: d.iter().map(|b| format!("{b:02x}")).collect(),
        }
    }
}

/// The DiSSECT families no record carries, and why.
pub const NOT_EVALUATED_FAMILIES: &[(&str, &str)] = &[
    (
        "brainpool_overlap",
        "defined on Brainpool's seed-derived bit strings; the registry records no seeds",
    ),
    (
        "multiples_x",
        "DiSSECT's definition depends on a choice of base point this module does not reproduce",
    ),
    (
        "q_torsion",
        "torsion of a lift to Q: the lift is a choice DiSSECT makes and this module does not",
    ),
    (
        "square_4p1",
        "the square part of 4q − 1 needs factoring numbers of the size of q; not attempted",
    ),
    (
        "class_number",
        "only exact h(d_K) for |d_K| ≤ 10^7 (frobenius.class_number); no analytic bounds",
    ),
    (
        "weierstrass",
        "the coefficients are in the registry's model_json; not copied",
    ),
];

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TraitFile {
    pub schema: String,
    pub generated_by: String,
    pub what_this_is: String,
    /// SHA-256 of the registry bytes the records were derived from.
    pub registry_sha256: String,
    /// Every file read: registry, factor hints, source manifest, sources.
    pub inputs: Vec<InputDigest>,
    /// Rho iterations per composite.
    pub factor_budget: u64,
    pub default_signature: Vec<String>,
    pub not_evaluated_families: BTreeMap<String, String>,
    pub curves: Vec<CurveTraits>,
}

impl TraitFile {
    /// The registry alone: no hints, no other sources.
    pub fn build(registry_text: &str, budget: u64) -> Result<Self, String> {
        let registry = Input {
            path: "docs/curves/registry.json".into(),
            text: registry_text.to_string(),
        };
        Ok(Self::build_from(&registry, None, None, &[], budget)?.0)
    }

    /// The registry, the factor hints, and the curves of every source the
    /// manifest names (`sources[i]` is the file at the manifest's `i`-th
    /// path).  Returns the [`Factoring`] too, for its leftovers.
    pub fn build_from(
        registry: &Input,
        hints: Option<&Input>,
        manifest: Option<&Input>,
        sources: &[Input],
        budget: u64,
    ) -> Result<(Self, Factoring), String> {
        use rayon::prelude::*;
        let mut curves = read_registry(&registry.text)?;
        let mut inputs = vec![registry.digest()];
        let hint_primes = match hints {
            Some(h) => {
                inputs.push(h.digest());
                read_hints(&h.text)?
            }
            None => Vec::new(),
        };
        if let Some(m) = manifest {
            inputs.push(m.digest());
            let specs = sources::read_manifest(&m.text)?;
            if specs.len() != sources.len()
                || specs.iter().zip(sources).any(|(s, i)| s.path != i.path)
            {
                return Err("the source files do not match the manifest".into());
            }
            let mut seen: std::collections::BTreeSet<String> =
                curves.iter().map(|c| c.slug.clone()).collect();
            for (spec, input) in specs.iter().zip(sources) {
                inputs.push(input.digest());
                let read = match spec.kind.as_str() {
                    "isogeny-walk" => sources::read_walk(spec, &input.text)?,
                    _ => sources::read_crosswalk(spec, &input.text)?,
                };
                curves.extend(read.into_iter().filter(|c| seen.insert(c.slug.clone())));
            }
        }
        let fx = Factoring::new(budget, hint_primes);
        let records: Result<Vec<_>, String> = curves.par_iter().map(|c| compute(c, &fx)).collect();
        let file = TraitFile {
            schema: SCHEMA.into(),
            generated_by: "cargo run --release --bin curve_traits -- build".into(),
            what_this_is: "Traits derived from each curve's model and order (the registry, then \
                           the sources in docs/curves/trait-sources.json), with the status of \
                           each value; see docs/curves/TRAITS.md."
                .into(),
            registry_sha256: inputs[0].sha256.clone(),
            inputs,
            factor_budget: budget,
            default_signature: DEFAULT_SIGNATURE.iter().map(|s| s.to_string()).collect(),
            not_evaluated_families: NOT_EVALUATED_FAMILIES
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect(),
            curves: records?,
        };
        Ok((file, fx))
    }

    /// Pretty JSON with a trailing newline, as committed.
    pub fn to_json(&self) -> String {
        let mut s = serde_json::to_string_pretty(self).expect("serialisable");
        s.push('\n');
        s
    }

    pub fn find(&self, name: &str) -> Option<&CurveTraits> {
        let slug = crate::cryptanalysis::curve_id::resolve(name).unwrap_or(name);
        self.curves.iter().find(|c| c.slug == slug)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn coordinate_models_keep_exact_fields_and_available_order_certificates() {
        use super::*;
        let (montgomery, _) = parse_model(
            r#"{"field":"fp-101","form":"B*y^2=x^3+A*x^2+x","p":"101","A":"3","B":"1"}"#,
        )
        .unwrap();
        assert_eq!(montgomery.q(), BigUint::from(101u64));
        // The prime-field Montgomery adapter now supports exact counting;
        // field size and curve order are distinct quantities.
        assert_eq!(
            certify::check_order(&montgomery, &BigUint::from(96u64), &[]).unwrap(),
            (Status::Proved, "exhaustive point count")
        );
        assert!(certify::check_order(&montgomery, &BigUint::from(101u64), &[]).is_err());

        let (extension, _) = parse_model(
            r#"{"field":"fpk-37-2-example","form":"y^2=x^3+a*x+b","p":"37","k":"2","modulus":["2","0"],"a":["1","0"],"b":["1","0"]}"#,
        )
        .unwrap();
        assert_eq!(extension.q(), BigUint::from(1369u64));
        assert_eq!(
            certify::check_order(&extension, &BigUint::from(1369u64), &[])
                .unwrap()
                .0,
            Status::NotEvaluated
        );
    }
    use super::*;

    const REGISTRY: &str = include_str!("../../../docs/curves/registry.json");

    fn fx() -> Factoring {
        Factoring::new(DEFAULT_BUDGET, Vec::new())
    }

    fn curve(slug: &str) -> RegistryCurve {
        read_registry(REGISTRY)
            .unwrap()
            .into_iter()
            .find(|c| c.slug == slug)
            .unwrap()
    }

    #[test]
    fn every_registry_entry_parses_with_its_order_identity() {
        let curves = read_registry(REGISTRY).unwrap();
        assert!(curves.len() >= 100);
        let with_gen = curves
            .iter()
            .filter(|c| !c.representations.is_empty())
            .count();
        assert!(with_gen >= 80, "only {with_gen} representations parsed");
    }

    #[test]
    fn extension_registry_preserves_field_and_subgroup_records() {
        let curves = read_registry(REGISTRY).unwrap();
        let extensions: Vec<_> = curves
            .iter()
            .filter(|curve| matches!(&curve.model, Model::Extension { .. }))
            .collect();
        assert_eq!(extensions.len(), 1324);
        for curve in extensions {
            let Model::Extension {
                p,
                k,
                q,
                field_id,
                modulus,
                a,
                b,
            } = &curve.model
            else {
                unreachable!()
            };
            assert_eq!(q, &p.pow(*k));
            assert!(field_id.starts_with(&format!("fpk-{p}-{k}-")));
            assert_eq!(
                (modulus.len(), a.len(), b.len()),
                (*k as usize, *k as usize, *k as usize)
            );
            let Some(rep) = curve.representations.first() else {
                // Catalogue entries may record field data without a
                // generator representation.
                continue;
            };
            assert_eq!(&rep.subgroup_order * &rep.cofactor, curve.order);
            assert!(matches!(
                &rep.generator,
                GeneratorCoords::Extension(x, y) if x.len() == *k as usize && y.len() == *k as usize
            ));
            assert_eq!(
                certify::check_order(&curve.model, &curve.order, &curve.representations)
                    .unwrap()
                    .0,
                Status::NotEvaluated
            );
        }
    }

    #[test]
    fn registered_montgomery_and_edwards_generators_certify() {
        let doc: Value = serde_json::from_str(REGISTRY).unwrap();
        let curves = read_registry(REGISTRY).unwrap();
        let mut checked = 0;
        for raw in doc["curves"].as_array().unwrap() {
            let model: Value = serde_json::from_str(raw["model_json"].as_str().unwrap()).unwrap();
            let form = model["form"].as_str().unwrap();
            if !matches!(
                form,
                "B*y^2=x^3+A*x^2+x" | "a*x^2+y^2=1+d*x^2*y^2" | "x^2+y^2=c^2*(1+d*x^2*y^2)"
            ) || raw["representations"].as_array().unwrap().is_empty()
            {
                continue;
            }
            let slug = raw["slug"].as_str().unwrap();
            let curve = curves.iter().find(|c| c.slug == slug).unwrap();
            assert!(!curve.representations.is_empty(), "{slug}");
            let (status, _) =
                certify::check_order(&curve.model, &curve.order, &curve.representations).unwrap();
            assert!(status <= Status::Bounded, "{slug}: {status:?}");
            checked += 1;
        }
        assert_eq!(checked, 20);
    }

    #[test]
    fn a_small_koblitz_curve() {
        let t = compute(&curve("icv1-f2m7-t13-616700dd"), &fx()).unwrap();
        assert_eq!(t.order_check.status, Status::Proved);
        assert_eq!(t.frobenius.cm_disc.as_deref(), Some("-7"));
        // 13² − 4·128 = −343 = −7·7²: v = 7.
        assert_eq!(t.frobenius.conductor, "7");
        assert_eq!(t.frobenius.class_number, Some(1));
        assert_eq!(t.keys["descent"], "k1|t|=1");
        assert_eq!(t.keys["descent_signed"], "k1|t|=-1");
        assert_eq!(t.keys["cm"], "-7");
        assert_eq!(t.keys["cofactor"], "4");
        let seven = t.small_primes.iter().find(|s| s.ell == 7).unwrap();
        assert_eq!((seven.depth, seven.splitting.as_str()), (1, "ramified"));
    }

    #[test]
    fn ecc2k_130_conductor_is_exact() {
        let t = compute(&curve("icv1-f2m131-tm22283658519494248867-115e0dc5"), &fx()).unwrap();
        assert_eq!(t.frobenius.cm_disc.as_deref(), Some("-7"));
        assert_eq!(
            t.frobenius.conductor_primes,
            vec![
                ("263".to_string(), 1),
                ("146505763881528721".to_string(), 1)
            ]
        );
        assert_eq!(t.order_check.status, Status::Proved);
        assert!(t.order_check.method.starts_with("descent count"));
        assert_eq!(t.keys["descent"], "k1|t|=1");
    }

    #[test]
    fn endomorphism_rings_from_contained_orders() {
        let get = |slug: &str| compute(&curve(slug), &fx()).unwrap();
        // K_0 over GF(2^7): Z[π] has conductor 7, but Z[τ] is maximal.
        let t = get("icv1-f2m7-t13-616700dd");
        assert_eq!(t.frobenius.conductor, "7");
        assert_eq!(t.endomorphism.maximal, Some(true));
        assert_eq!(t.keys["end_order"], "maximal");
        // The GF(2^18) twist: π_2 with t₂ = 3, 9 − 16 = −7, is maximal.
        let t = get("icv1-f2m18-t999-40751283");
        assert_eq!(t.frobenius.conductor, "85");
        assert_eq!(t.endomorphism.maximal, Some(true));
        assert!(t
            .endomorphism
            .contains
            .iter()
            .any(|o| o.starts_with("Z[π_2]")));
        // secp256k1: j = 0 puts Z[ζ₃] in End(E) though v has 128 bits.
        let t = get("icv1-fp256-t432420386565659656852420866390673177327-76dadd18");
        assert_eq!(t.endomorphism.maximal, Some(true));
        // A generic prime curve with v = 7 at depth 1: Frobenius is a
        // scalar on E[7], so End(E) is above the floor: O_K.
        let t = get("icv1-fp10-t5-192cb216");
        assert_eq!(t.endomorphism.conductor.as_deref(), Some("1"));
        assert_eq!(t.endomorphism.maximal, Some(true));
        assert_eq!(t.keys["end_order"], "maximal");
        // v = 271 is beyond the scalar test: End(E) is Z[π] or O_K, unresolved.
        let t = get("icv1-f2m17-tm101-8da609f8");
        assert_eq!(t.endomorphism.conductor_divides.as_deref(), Some("271"));
        assert_eq!(t.endomorphism.maximal, None);
        assert_eq!(t.endomorphism.status, Status::Bounded);
        assert_eq!(t.keys["end_order"], "unresolved");
    }

    #[test]
    fn the_gf2_18_twist_groups_with_koblitz_by_cm_not_by_descent() {
        let t = compute(&curve("icv1-f2m18-t999-40751283"), &fx()).unwrap();
        assert_eq!(t.keys["cm"], "-7");
        assert_eq!(t.keys["jfield"], "k1");
        assert_eq!(t.keys["descent"], "k2|t|=3");
    }
}
