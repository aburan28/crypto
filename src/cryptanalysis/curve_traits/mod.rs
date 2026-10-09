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
pub mod query;
pub mod subfield;

use std::collections::BTreeMap;

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use arith::{factor, Factored};

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
}

impl Model {
    /// The field size `q`.
    pub fn q(&self) -> BigUint {
        match self {
            Model::Binary { n, .. } => BigUint::one() << *n,
            Model::Prime { p, .. } => p.clone(),
        }
    }

    /// `n` for `GF(2^n)`, the bit length of `p` for `GF(p)`.
    pub fn size_bits(&self) -> u64 {
        match self {
            Model::Binary { n, .. } => u64::from(*n),
            Model::Prime { p, .. } => p.bits(),
        }
    }
}

/// One recorded representation, mapped into `Model` coordinates: a generator
/// of the prime subgroup.
#[derive(Clone, Debug)]
pub struct Representation {
    pub generator: (BigUint, BigUint),
    pub subgroup_order: BigUint,
    pub cofactor: BigUint,
}

/// A registry entry, parsed losslessly.
#[derive(Clone, Debug)]
pub struct RegistryCurve {
    pub slug: String,
    pub family: String,
    pub model: Model,
    pub trace: BigInt,
    pub order: BigUint,
    /// The registry's certified `End(E)` discriminant, when it has one.
    pub end: Option<BigInt>,
    pub representations: Vec<Representation>,
}

fn parse_uint(s: &str) -> Option<BigUint> {
    match s.strip_prefix("0x") {
        Some(h) => BigUint::parse_bytes(h.as_bytes(), 16),
        None => s.parse().ok(),
    }
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
    }
    let gen = curve["generator"].as_array()?;
    let coord = |v: &Value| parse_uint(v.as_str()?);
    let cofactor = match &curve["cofactor"] {
        Value::String(s) => parse_uint(s)?,
        v => BigUint::from(v.as_u64()?),
    };
    let x = coord(gen.first()?)?;
    let y = coord(gen.get(1)?)?;
    let generator = match model {
        Model::Binary { .. } => (x, y),
        Model::Prime { p, .. } => short_coordinates(x, y, p, map)?,
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

/// What is known of `End(E)`, an order of `Q(π)` between `Z[π]` and
/// `O_K`.  For an ordinary curve every geometric endomorphism is defined
/// over the working field, so this is also the ring over that field.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EndomorphismRecord {
    /// `f₀` with the conductor `[O_K : End(E)]` dividing it: the conductor
    /// `v` of `Z[π]`, sharpened by every other order known to lie in
    /// `End(E)`.
    pub conductor_divides: Option<String>,
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

/// Everything derived about one registered curve.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurveTraits {
    pub slug: String,
    pub family: String,
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
    /// The grouping keys of [`KEYS`], precomputed.
    pub keys: BTreeMap<String, String>,
}

// ── Computing it ──────────────────────────────────────────────────────

fn frobenius(c: &RegistryCurve, budget: u64) -> Result<(Frobenius, BigInt), String> {
    let q = BigInt::from(c.model.q());
    let delta: BigInt = &c.trace * &c.trace - BigInt::from(4u8) * &q;
    let fac = factor(delta.magnitude(), budget);
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

fn subgroup(c: &RegistryCurve, budget: u64) -> Result<(Subgroup, Option<BigUint>), String> {
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
    let fac = factor(&c.order, budget);
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

fn embedding(q: &BigUint, r: Option<&BigUint>, budget: u64) -> Embedding {
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
    let fac: Factored = factor(&rm1, budget);
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

fn twist(q: &BigUint, trace: &BigInt, budget: u64) -> Twist {
    let order = (BigInt::from(q.clone()) + BigInt::one() + trace)
        .to_biguint()
        .expect("Hasse: q + 1 + t > 0");
    let fac = factor(&order, budget);
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
    if let Some(end) = &c.end {
        let f = (end / &dk).sqrt();
        if !(&f0 % f.magnitude()).is_zero() {
            return Err(format!(
                "registry certifies End(E) of conductor {f}, which does not divide {f0}"
            ));
        }
        contains.push(format!(
            "End(E) of discriminant {end} (registry certificate)"
        ));
        return Ok(record(
            Some(f.to_string()),
            Some(f.is_one()),
            status,
            contains,
        ));
    }
    Ok(if f0.is_one() {
        record(Some("1".into()), Some(true), status, contains)
    } else {
        record(Some(f0.to_string()), None, Status::Bounded, contains)
    })
}

/// Derive every trait of one curve.  `Err` when the registry entry
/// contradicts itself (an order the model refutes, a certified `End(E)`
/// outside `Q(π)`).
pub fn compute(c: &RegistryCurve, budget: u64) -> Result<CurveTraits, String> {
    let ctx = |e: String| format!("{}: {e}", c.slug);
    let q = c.model.q();
    let characteristic = match &c.model {
        Model::Binary { .. } => BigUint::from(2u8),
        Model::Prime { p, .. } => p.clone(),
    };
    let (status, method) =
        certify::check_order(&c.model, &c.order, &c.representations).map_err(ctx)?;
    let (frob, delta) = frobenius(c, budget).map_err(ctx)?;
    let ordinary = !(&c.trace % BigInt::from(characteristic.clone())).is_zero();
    let trace_ratio = {
        let r = arith::log2_abs(&c.trace) - 1.0 - arith::log2_abs(&BigInt::from(q.clone())) / 2.0;
        let mag = if c.trace.is_zero() { 0.0 } else { r.exp2() };
        arith::round_to(if c.trace.is_negative() { -mag } else { mag }, 6)
    };
    let (subgroup, r) = subgroup(c, budget).map_err(ctx)?;
    let sub = subfield_record(c).map_err(ctx)?;
    let endo = endomorphism(c, &frob, &sub, ordinary).map_err(ctx)?;
    let mut out = CurveTraits {
        slug: c.slug.clone(),
        family: c.family.clone(),
        field: FieldInfo {
            kind: match c.model {
                Model::Binary { .. } => "binary",
                Model::Prime { .. } => "prime",
            }
            .into(),
            degree: match &c.model {
                Model::Binary { n, .. } => *n,
                Model::Prime { .. } => 1,
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
        small_primes: small_primes(&delta, &characteristic),
        frobenius: frob,
        endomorphism: endo,
        subfield: sub,
        embedding: embedding(&q, r.as_ref(), budget),
        twist: twist(&q, &c.trace, budget),
        subgroup,
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
        "'maximal' when End(E) = O_K is established; 'non-maximal' when certified smaller; 'unresolved'; 'unknown'",
    ),
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

pub const SCHEMA: &str = "curve-traits/1";

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TraitFile {
    pub schema: String,
    pub generated_by: String,
    pub what_this_is: String,
    /// SHA-256 of the registry bytes the records were derived from.
    pub registry_sha256: String,
    /// Rho iterations per composite.
    pub factor_budget: u64,
    pub default_signature: Vec<String>,
    pub curves: Vec<CurveTraits>,
}

impl TraitFile {
    pub fn build(registry_text: &str, budget: u64) -> Result<Self, String> {
        use rayon::prelude::*;
        let curves = read_registry(registry_text)?;
        let records: Result<Vec<_>, String> =
            curves.par_iter().map(|c| compute(c, budget)).collect();
        let digest = crate::hash::sha256::sha256(registry_text.as_bytes());
        Ok(TraitFile {
            schema: SCHEMA.into(),
            generated_by: "cargo run --release --bin curve_traits -- build".into(),
            what_this_is: "Traits derived from each registered curve's model and order, with the \
                           status of each value; see docs/curves/TRAITS.md."
                .into(),
            registry_sha256: digest.iter().map(|b| format!("{b:02x}")).collect(),
            factor_budget: budget,
            default_signature: DEFAULT_SIGNATURE.iter().map(|s| s.to_string()).collect(),
            curves: records?,
        })
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
    use super::*;

    const REGISTRY: &str = include_str!("../../../docs/curves/registry.json");

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
        let t = compute(&curve("icv1-f2m7-t13-616700dd"), DEFAULT_BUDGET).unwrap();
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
        let t = compute(
            &curve("icv1-f2m131-tm22283658519494248867-115e0dc5"),
            DEFAULT_BUDGET,
        )
        .unwrap();
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
        let get = |slug: &str| compute(&curve(slug), DEFAULT_BUDGET).unwrap();
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
        // A generic prime curve with v = 7: End(E) is Z[π] or O_K, unresolved.
        let t = get("icv1-fp10-t5-192cb216");
        assert_eq!(t.endomorphism.conductor_divides.as_deref(), Some("7"));
        assert_eq!(t.endomorphism.maximal, None);
        assert_eq!(t.endomorphism.status, Status::Bounded);
        assert_eq!(t.keys["end_order"], "unresolved");
    }

    #[test]
    fn the_gf2_18_twist_groups_with_koblitz_by_cm_not_by_descent() {
        let t = compute(&curve("icv1-f2m18-t999-40751283"), DEFAULT_BUDGET).unwrap();
        assert_eq!(t.keys["cm"], "-7");
        assert_eq!(t.keys["jfield"], "k1");
        assert_eq!(t.keys["descent"], "k2|t|=3");
    }
}
