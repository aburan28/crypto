//! # Jacobian transfer discovery across the curve registry.
//!
//! For every curve in `docs/curves/registry.json` this module enumerates the
//! known ways of moving its discrete logarithm into the Jacobian of a curve
//! of higher genus, prices each route with one frozen formula-level cost
//! model ([`cost::CostModel`]), and reports the ratio of that price to the
//! Pollard-rho reference on the original subgroup.  The result is
//! `docs/curves/jacobian-transfers.json`, rebuilt and checked by
//! `cargo run --release --bin jacobian_transfers`; see
//! `docs/curves/JACOBIAN_TRANSFERS.md`.
//!
//! Routes, by field:
//!
//! | field | route | what is exact | what is estimated |
//! |:--|:--|:--|:--|
//! | `GF(2^N)` | GHS Weil descent, every `N = n * l` | magic number, type, genus | solver cost |
//! | `GF(2^N)` | GHS isogeny extension (Galbraith–Hess–Smart) | annihilator-degree census | class membership, walk cost |
//! | `GF(2^N)` | same-field genus-3 cover (`COVERS.md`) | genus, degree | solver cost |
//! | `GF(p)` | same-field genus-2 cover (`COVERS.md`) | genus, degree | solver cost |
//! | `GF(p)` | base change to `GF(p^k)` then a genus-`k` descent | nothing (existence open) | solver cost at the most favourable genus |
//! | `GF(p^k)` | Gaudry–Diem decomposition on `E(F_{p^k})` | — | relation and decomposition cost |
//! | `GF(p^k)` | descent to a genus-`k` cover over `GF(p)` | nothing (existence open) | solver cost at the most favourable genus |
//!
//! Every row is a formula.  A `formula_advantage` row says the literature's
//! complexity formulas, with the stated constants, put the route below rho on
//! that subgroup by at least the margin; it is a lead with named open
//! obligations (`AGENTS.md`, curve-structure rules), not a speedup, a
//! measured result or a security claim.  A `boundary` row says the formulas
//! themselves place the route at or above rho, so no constant-factor work on
//! it can cross rho.

pub mod binary;
pub mod cost;
pub mod report;

use std::collections::BTreeMap;

use num_bigint::{BigInt, BigUint};
use num_traits::{ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::curve_traits::arith::{factor, log2_abs};
use crate::cryptanalysis::ec_trapdoor::DescentRow;
use cost::{CostModel, DestinationCost, DestinationSolver};

pub const SCHEMA: &str = "jacobian-transfers/1";

/// Pollard-rho iterations spent factoring an order whose subgroup the
/// registry does not record.
pub const DEFAULT_FACTOR_BUDGET: u64 = 200_000;

// ── Registry input ────────────────────────────────────────────────────

/// The field a registered model lives over.
#[derive(Clone, Debug)]
pub enum Field {
    /// `GF(2^n)` with its reduction polynomial and the curve's `a`, `b`.
    Binary {
        n: u32,
        irr: IrreduciblePoly,
        a: F2mElement,
        b: F2mElement,
    },
    /// `GF(p)`.
    Prime { p: BigUint },
    /// `GF(p^k)`, `k >= 2`.
    Extension { p: BigUint, k: u32 },
}

/// A registry entry, read for transfer screening.
#[derive(Clone, Debug)]
pub struct RegistryCurve {
    pub slug: String,
    pub family: String,
    pub icv1: String,
    pub standard_names: Vec<String>,
    pub field: Field,
    pub trace: BigInt,
    pub order: BigUint,
    /// The prime subgroup order the registry records, if any.
    pub subgroup_order: Option<BigUint>,
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

fn parse_field(model: &Value, label: &str) -> Result<Field, String> {
    let num = |key: &str| {
        str_field(model, key)
            .and_then(|s| parse_uint(s).ok_or_else(|| format!("bad {key}: {s}")))
    };
    if label.starts_with("f2m-") {
        let modulus = num("modulus")?;
        let n = (modulus.bits() - 1) as u32;
        if !label.starts_with(&format!("f2m-{n}-")) {
            return Err(format!("field {label} is not of degree {n}"));
        }
        let irr = IrreduciblePoly {
            degree: n,
            low_terms: (0..n).filter(|&i| modulus.bit(u64::from(i))).collect(),
        };
        let a = F2mElement::from_biguint(&num("a")?, n);
        let b = F2mElement::from_biguint(&num("b")?, n);
        Ok(Field::Binary { n, irr, a, b })
    } else if label.starts_with("fpk-") {
        let p = num("p")?;
        let k: u32 = str_field(model, "k")?
            .parse()
            .map_err(|_| "bad extension degree".to_string())?;
        if k < 2 {
            return Err("extension label with degree below 2".into());
        }
        Ok(Field::Extension { p, k })
    } else if label.starts_with("fp-") {
        Ok(Field::Prime { p: num("p")? })
    } else {
        Err(format!("unknown field label {label}"))
    }
}

/// Parse `docs/curves/registry.json`, every form of model, keeping what the
/// transfer screen needs.  The recorded order is checked against the ICV1
/// trace and the Hasse bound.
pub fn read_registry(text: &str) -> Result<Vec<RegistryCurve>, String> {
    let doc: Value = serde_json::from_str(text).map_err(|e| format!("registry JSON: {e}"))?;
    let curves = doc["curves"].as_array().ok_or("registry has no curves")?;
    let mut out = Vec::with_capacity(curves.len());
    for c in curves {
        let slug = str_field(c, "slug")?.to_string();
        let ctx = |e: String| format!("{slug}: {e}");
        let icv1 = str_field(c, "icv1")?.to_string();
        let parts: Vec<&str> = icv1.split(':').collect();
        let label = parts.get(1).ok_or_else(|| ctx("ICV1 has no field".into()))?;
        let model: Value = serde_json::from_str(str_field(c, "model_json")?)
            .map_err(|e| ctx(format!("model JSON: {e}")))?;
        let field = parse_field(&model, label).map_err(ctx)?;
        let trace: BigInt = parts
            .get(2)
            .and_then(|t| t.parse().ok())
            .ok_or_else(|| ctx("ICV1 string has no trace".into()))?;
        let order = parse_uint(str_field(c, "order")?).ok_or_else(|| ctx("bad order".into()))?;
        let q = field_size(&field);
        if BigInt::from(q.clone()) + 1 - &trace != BigInt::from(order.clone()) {
            return Err(ctx(format!("order {order} ≠ q + 1 − t for t = {trace}")));
        }
        if &trace * &trace > BigInt::from(q) * 4u8 {
            return Err(ctx(format!("|t| = |{trace}| exceeds the Hasse bound 2√q")));
        }
        let subgroup_order = c["representations"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|r| r["curve"]["subgroup_order"].as_str().and_then(parse_uint))
            .max();
        let standard_names = c["standard_names"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|v| v.as_str().map(str::to_string))
            .collect();
        out.push(RegistryCurve {
            slug,
            family: str_field(c, "family")?.to_string(),
            icv1,
            standard_names,
            field,
            trace,
            order,
            subgroup_order,
        });
    }
    Ok(out)
}

fn field_size(f: &Field) -> BigUint {
    match f {
        Field::Binary { n, .. } => BigUint::from(1u32) << *n,
        Field::Prime { p } => p.clone(),
        Field::Extension { p, k } => p.pow(*k),
    }
}

/// Round a stored value to three decimals and at most 12 significant digits,
/// so the committed file does not depend on the last bits of the platform's
/// `log2`/`exp2`/`ln` (CI rebuilds on Linux what was written on macOS).
pub fn round_to(x: f64, digits: i32) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let scale = 10f64.powi(digits);
    let r = (x * scale).round() / scale;
    if r.abs() >= 1e9 {
        format!("{r:.11e}").parse().expect("formatted float")
    } else {
        r
    }
}

/// `log2 x` for a positive integer, to double precision.
pub fn log2_biguint(x: &BigUint) -> f64 {
    if x.is_zero() {
        return f64::NEG_INFINITY;
    }
    let bits = x.bits();
    if bits <= 53 {
        return x.to_f64().expect("fits").log2();
    }
    let top = (x >> (bits - 53)).to_f64().expect("53 bits");
    top.log2() + (bits - 53) as f64
}

// ── The record ────────────────────────────────────────────────────────

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FieldInfo {
    /// `binary`, `prime` or `extension`.
    pub kind: String,
    /// `GF(2^n)`: `n`; `GF(p)`: 1; `GF(p^k)`: `k`.
    pub degree: u32,
    /// `log2` of the field size.
    pub log2_size: f64,
    /// `log2 p` for an extension field, else equal to `log2_size`.
    pub log2_base: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SubgroupInfo {
    /// Decimal prime order `r`, or the bound used when it is not established.
    pub order: String,
    pub bits: f64,
    /// `registry_representation`, `largest_prime_of_order`, or
    /// `largest_unfactored_composite_bound` (an upper bound on `r`).
    pub status: String,
}

/// How a route's destination relates to the subgroup.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Transport {
    /// The subgroup is carried whenever `r` is coprime to the map degree
    /// (a same-field cover pulls back injectively on the `r`-torsion).
    CarriedIfCoprimeToDegree,
    /// The Jacobian is large enough to contain the subgroup; carrying it is a
    /// necessary-condition check only.
    SizeAdmits,
    /// The Jacobian cannot contain a subgroup of order `r`.
    TooSmall,
    /// The descent curve is `E` itself over the subfield: nothing is gained.
    Degenerate,
    /// Not applicable: the route stays on `E`.
    SameGroup,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum RouteClass {
    /// The formulas, with the stated constants, place the route at least
    /// the margin below rho and the subgroup can be carried.
    FormulaAdvantage,
    /// The formulas place the route below rho but by less than the margin:
    /// the exponents favour it and the constants decide.
    FormulaMarginal,
    /// One bound of an unresolved invariant (the GHS magic number when the
    /// trace of `a` leaves it at `m_b` or `m_b + 1`) puts the route below rho
    /// and the other does not: the verdict waits on that invariant.
    Unresolved,
    /// The formulas place the route at or above rho.
    Boundary,
    /// The subgroup is not carried, so the route solves a different problem.
    NotCarried,
}

/// How much of a row is established.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Evidence {
    /// Geometry computed exactly on this curve; only the solver cost is a formula.
    ExactInvariant,
    /// The destination's existence is established by a construction family
    /// but not exercised; the cost is a formula.
    Construction,
    /// The destination's existence is not established; the row prices the
    /// most favourable case the literature allows.
    Proposed,
    /// Counting and class-membership heuristics enter the price.
    Heuristic,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Destination {
    /// `Jac(C)` or `E(F_{p^k})`.
    pub object: String,
    /// Decimal genus.
    pub genus: String,
    pub log2_genus: f64,
    /// `log2` of the field the destination is defined over.
    pub field_bits: f64,
    /// `log2 #Jac`, taken as `genus * field_bits`.
    pub log2_order: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RouteCost {
    /// Finding the route (an isogeny walk); absent when none is needed.
    pub log2_search: Option<f64>,
    /// Evaluating the transfer map once; absent when the route stays on `E`.
    pub log2_transfer: Option<f64>,
    pub destination: DestinationCost,
    /// `log2(2^search + 2^transfer + 2^destination)`.
    pub log2_total: f64,
    /// `reference - total`: positive means the formulas favour the route.
    pub ratio_log2: f64,
    /// `0.5 log2 r - (leading term of the destination solver)`: the ratio with
    /// every constant stripped, so exponent and constant effects can be told apart.
    pub exponent_ratio_log2: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TransferRoute {
    pub kind: String,
    pub label: String,
    pub evidence: Evidence,
    /// Exact or estimated structural data, by name.
    pub invariants: BTreeMap<String, Value>,
    pub destination: Destination,
    pub transport: Transport,
    pub cost: RouteCost,
    pub class: RouteClass,
    pub notes: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurveTransfers {
    pub slug: String,
    pub icv1: String,
    pub standard_names: Vec<String>,
    pub family: String,
    pub field: FieldInfo,
    pub trace: String,
    pub subgroup: SubgroupInfo,
    /// Pollard rho on the subgroup of `E`, in the model's unit.
    pub reference_log2: f64,
    /// Best class over the routes.
    pub verdict: RouteClass,
    /// Index into `routes` of the cheapest carried route.
    pub best_route: Option<usize>,
    pub routes: Vec<TransferRoute>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Summary {
    pub curves: usize,
    pub routes: usize,
    pub by_class: BTreeMap<String, usize>,
    pub by_route_kind: BTreeMap<String, usize>,
    /// Slugs whose verdict is `formula_advantage` or `formula_marginal`, with
    /// the route that gives it.
    pub formula_advantage: Vec<(String, String)>,
    /// Slugs whose verdict is `unresolved`, with the favourable-bound route.
    pub unresolved: Vec<(String, String)>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TransferFile {
    pub schema: String,
    pub generated_by: String,
    pub what_this_is: String,
    pub registry_sha256: String,
    pub factor_budget: u64,
    pub cost_model: CostModel,
    pub summary: Summary,
    pub curves: Vec<CurveTransfers>,
}

// ── The scan ──────────────────────────────────────────────────────────

fn subgroup_info(c: &RegistryCurve, budget: u64) -> SubgroupInfo {
    if let Some(r) = &c.subgroup_order {
        return SubgroupInfo {
            order: r.to_string(),
            bits: round_to(log2_biguint(r), 3),
            status: "registry_representation".into(),
        };
    }
    let f = factor(&c.order, budget);
    let largest_prime = f.primes.iter().map(|(p, _)| p.clone()).max();
    let largest_composite = f.composites.iter().map(|(c, _)| c.clone()).max();
    match (largest_prime, largest_composite) {
        (Some(p), Some(comp)) if comp > p => SubgroupInfo {
            order: comp.to_string(),
            bits: round_to(log2_biguint(&comp), 3),
            status: "largest_unfactored_composite_bound".into(),
        },
        (Some(p), _) => SubgroupInfo {
            order: p.to_string(),
            bits: round_to(log2_biguint(&p), 3),
            status: "largest_prime_of_order".into(),
        },
        (None, Some(comp)) => SubgroupInfo {
            order: comp.to_string(),
            bits: round_to(log2_biguint(&comp), 3),
            status: "largest_unfactored_composite_bound".into(),
        },
        (None, None) => SubgroupInfo {
            order: c.order.to_string(),
            bits: round_to(log2_biguint(&c.order), 3),
            status: "largest_prime_of_order".into(),
        },
    }
}

fn field_info(f: &Field) -> FieldInfo {
    match f {
        Field::Binary { n, .. } => FieldInfo {
            kind: "binary".into(),
            degree: *n,
            log2_size: f64::from(*n),
            log2_base: f64::from(*n),
        },
        Field::Prime { p } => {
            let b = round_to(log2_biguint(p), 3);
            FieldInfo {
                kind: "prime".into(),
                degree: 1,
                log2_size: b,
                log2_base: b,
            }
        }
        Field::Extension { p, k } => {
            let b = log2_biguint(p);
            FieldInfo {
                kind: "extension".into(),
                degree: *k,
                log2_size: round_to(b * f64::from(*k), 3),
                log2_base: round_to(b, 3),
            }
        }
    }
}

struct RouteBuilder<'a> {
    model: &'a CostModel,
    subgroup_bits: f64,
    reference_log2: f64,
    source_field_bits: f64,
}

impl RouteBuilder<'_> {
    fn jacobian_route(
        &self,
        kind: &str,
        label: String,
        evidence: Evidence,
        invariants: BTreeMap<String, Value>,
        genus: &BigUint,
        field_bits: f64,
        transport: Transport,
        log2_search: f64,
        notes: Vec<String>,
    ) -> TransferRoute {
        let log2_genus = log2_biguint(genus).max(0.0);
        let log2_order = log2_genus.exp2() * field_bits;
        let destination_cost = self
            .model
            .jacobian_destination(self.subgroup_bits, log2_genus, field_bits);
        // One evaluation of the transfer map: at least `g` operations over the
        // source field (the GHS conorm–norm map writes a genus-`g` divisor).
        let log2_transfer = log2_genus + self.model.log2_elliptic_op(self.source_field_bits);
        self.finish(
            kind,
            label,
            evidence,
            invariants,
            Destination {
                object: "Jac(C)".into(),
                genus: genus.to_string(),
                log2_genus: round_to(log2_genus, 3),
                field_bits: round_to(field_bits, 3),
                log2_order: round_to(log2_order, 3),
            },
            transport,
            log2_search,
            log2_transfer,
            destination_cost,
            notes,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        kind: &str,
        label: String,
        evidence: Evidence,
        invariants: BTreeMap<String, Value>,
        destination: Destination,
        transport: Transport,
        log2_search: f64,
        log2_transfer: f64,
        destination_cost: DestinationCost,
        notes: Vec<String>,
    ) -> TransferRoute {
        let log2_total = binary::log2_total(&[log2_search, log2_transfer, destination_cost.log2_cost]);
        let ratio_log2 = self.reference_log2 - log2_total;
        let exponent_ratio_log2 = 0.5 * self.subgroup_bits - destination_cost.log2_leading;
        let carried = matches!(
            transport,
            Transport::CarriedIfCoprimeToDegree | Transport::SizeAdmits | Transport::SameGroup
        );
        let class = if !carried {
            RouteClass::NotCarried
        } else if ratio_log2 >= self.model.advantage_margin_bits {
            RouteClass::FormulaAdvantage
        } else if ratio_log2 > 0.0 {
            RouteClass::FormulaMarginal
        } else {
            RouteClass::Boundary
        };
        let destination_cost = DestinationCost {
            solver: destination_cost.solver,
            log2_cost: round_to(destination_cost.log2_cost, 3),
            log2_leading: round_to(destination_cost.log2_leading, 3),
            considered: destination_cost
                .considered
                .iter()
                .map(|(s, c)| (*s, round_to(*c, 3)))
                .collect(),
            enge_gaudry_hypothesis_holds: destination_cost.enge_gaudry_hypothesis_holds,
        };
        TransferRoute {
            kind: kind.into(),
            label,
            evidence,
            invariants,
            destination,
            transport,
            cost: RouteCost {
                log2_search: log2_search.is_finite().then(|| round_to(log2_search, 3)),
                log2_transfer: log2_transfer.is_finite().then(|| round_to(log2_transfer, 3)),
                destination: destination_cost,
                log2_total: round_to(log2_total, 3),
                ratio_log2: round_to(ratio_log2, 3),
                exponent_ratio_log2: round_to(exponent_ratio_log2, 3),
            },
            class,
            notes,
        }
    }

    fn transport_by_size(&self, log2_genus: f64, field_bits: f64) -> Transport {
        if log2_genus.exp2() * field_bits + 1.0 >= self.subgroup_bits {
            Transport::SizeAdmits
        } else {
            Transport::TooSmall
        }
    }
}

fn inv<K: Into<String>, V: Into<Value>>(pairs: Vec<(K, V)>) -> BTreeMap<String, Value> {
    pairs.into_iter().map(|(k, v)| (k.into(), v.into())).collect()
}

fn binary_routes(rb: &RouteBuilder, c: &RegistryCurve, n_total: u32, irr: &IrreduciblePoly, a: &F2mElement, b: &F2mElement) -> Vec<TransferRoute> {
    let mut routes = Vec::new();
    let disc = &c.trace * &c.trace - BigInt::from(BigUint::from(1u32) << n_total) * 4;
    let log2_abs_disc = log2_abs(&disc);
    for audit in binary::ghs_rows(n_total, irr, a, b) {
        let l = audit.l;
        let magic = binary::ghs_magic(n_total, l, irr, a, b);
        let field_bits = f64::from(l);
        // When Tr(a) = 1 and n is even the reduced a may add one generator:
        // the GHS magic number is magic_b or magic_b + 1, and both are priced.
        let variants: Vec<(u32, &str)> = if magic.a_may_add_one {
            vec![(magic.magic_b, "lower"), (magic.magic_b + 1, "upper")]
        } else {
            vec![(magic.magic_b, "exact")]
        };
        let mut row = DescentRow::default();
        for (m, bound) in variants {
            let type_i = magic.orbit_length == m;
            let genus = crate::cryptanalysis::ec_trapdoor::ghs_genus_for_type(m, type_i);
            row = DescentRow { n: audit.n, l, magic_m: m, genus, type_i };
            let log2_genus = log2_biguint(&row.genus).max(0.0);
            let transport = if row.magic_m <= 1 {
                Transport::Degenerate
            } else {
                rb.transport_by_size(log2_genus, field_bits)
            };
            let mut notes = vec![
                "GHS Weil descent: magic number dim Span{(1, σ^i(√b))} and genus as in research/cover_route_prime_binary_20261008/ghs_magic.rs (Hess); the cover C/GF(2^l) is not constructed here".into(),
            ];
            if row.magic_m <= 1 {
                notes.push("sqrt(b) lies in the subfield: E is defined over GF(2^l) and the descent is the trace map onto E(GF(2^l)), which cannot carry the subgroup".into());
            }
            if bound != "exact" {
                notes.push(format!(
                    "Tr(a) = 1 with n even: a cannot be reduced into F_2 and may add one generator, so the magic number is {} or {}; this row prices the {bound} value and is heuristic until the reduced a is computed",
                    magic.magic_b,
                    magic.magic_b + 1
                ));
            }
            if audit.magic_m != magic.magic_b && audit.magic_m != magic.magic_b + u32::from(magic.a_may_add_one) {
                notes.push(format!(
                    "ec_trapdoor::audit_curve reports magic number {} here: it adds the Frobenius differences of an unreduced a as generators, while a matters only modulo {{α²+α}}",
                    audit.magic_m
                ));
            }
            let (label_suffix, evidence) = match bound {
                "exact" => (String::new(), Evidence::ExactInvariant),
                b => (format!(" [m = {m}, {b} bound]"), Evidence::Heuristic),
            };
            routes.push(rb.jacobian_route(
                "ghs_weil_descent",
                format!("GHS descent to Jac(C)/GF(2^{l}), N = {} * {l}{label_suffix}", row.n),
                evidence,
                inv(vec![
                    ("relative_degree", Value::from(row.n)),
                    ("base_degree", Value::from(l)),
                    ("magic_number", Value::from(row.magic_m)),
                    ("magic_number_b", Value::from(magic.magic_b)),
                    ("magic_number_bound", Value::from(bound)),
                    ("trace_a", Value::from(magic.trace_a)),
                    ("audit_curve_magic_number", Value::from(audit.magic_m)),
                    ("type_i", Value::from(row.type_i)),
                    ("cover_degree_log2", Value::from(row.magic_m)),
                ]),
                &row.genus,
                field_bits,
                transport,
                f64::NEG_INFINITY,
                notes,
            ));
        }
        let row = DescentRow { magic_m: magic.magic_b, ..row };
        // Galbraith–Hess–Smart isogeny extension, priced heuristically.
        let estimates = binary::isogeny_extension_estimates(rb.model, n_total, l, row.n, row.magic_m, log2_abs_disc);
        let mut chosen: Vec<&binary::IsogenyExtensionEstimate> = Vec::new();
        if let Some(first) = estimates.first() {
            chosen.push(first);
        }
        let best_plausible = estimates
            .iter()
            .filter(|e| e.log2_expected_isogenous >= 0.0)
            .min_by(|x, y| {
                let cx = route_cost_preview(rb, x, field_bits);
                let cy = route_cost_preview(rb, y, field_bits);
                cx.partial_cmp(&cy).unwrap_or(std::cmp::Ordering::Equal)
            });
        if let Some(bp) = best_plausible {
            if !chosen.iter().any(|e| e.annihilator_degree == bp.annihilator_degree) {
                chosen.push(bp);
            }
        }
        for e in chosen {
            let genus = binary::generic_genus(e.magic_bound);
            let log2_genus = log2_biguint(&genus).max(0.0);
            let transport = rb.transport_by_size(log2_genus, field_bits);
            let plausible = e.log2_expected_isogenous >= 0.0;
            routes.push(rb.jacobian_route(
                "ghs_isogeny_extension",
                format!(
                    "isogenous curve with magic number <= {} over the same tower, descent to Jac(C)/GF(2^{l})",
                    e.magic_bound
                ),
                Evidence::Heuristic,
                inv(vec![
                    ("relative_degree", Value::from(row.n)),
                    ("base_degree", Value::from(l)),
                    ("direct_magic_number", Value::from(row.magic_m)),
                    ("target_annihilator_degree", Value::from(e.annihilator_degree)),
                    ("target_magic_number_bound", Value::from(e.magic_bound)),
                    ("log2_candidate_curves", Value::from(round_to(e.log2_candidate_curves, 3))),
                    ("log2_isogeny_class_size", Value::from(round_to(0.5 * log2_abs_disc - std::f64::consts::PI.log2(), 3))),
                    ("log2_expected_isogenous_candidates", Value::from(round_to(e.log2_expected_isogenous, 3))),
                    ("isogenous_candidate_expected", Value::from(plausible)),
                ]),
                &genus,
                field_bits,
                transport,
                e.log2_search_cost,
                vec![
                    "heuristic: counts curves whose sqrt(b') has a small Frobenius annihilator and assumes isogeny classes are random of size sqrt(|t^2-4q|)/pi; no isogenous curve is exhibited".into(),
                    if plausible {
                        "at least one such isogenous curve is expected; the search cost prices a random walk in the class until one is hit".into()
                    } else {
                        "fewer than one such isogenous curve is expected in the class; the route is priced for the record".into()
                    },
                ],
            ));
        }
    }
    // Same-field genus-3 cover from COVERS.md (ordinary binary, b != 0).
    if !b.is_zero() {
        routes.push(rb.jacobian_route(
            "same_field_cover",
            format!("same-field genus-3 cover H -> E of degree 3 over GF(2^{n_total}) (ordinary_binary_cubic_pullback)"),
            Evidence::Construction,
            inv(vec![("map_degree", Value::from(3u32)), ("construction", Value::from("docs/curves/COVERS.md ordinary binary"))]),
            &BigUint::from(3u32),
            f64::from(n_total),
            Transport::CarriedIfCoprimeToDegree,
            f64::NEG_INFINITY,
            vec!["the Jacobian over the same field is larger than E: index calculus there is q^(4/3) against rho's q^(1/2)".into()],
        ));
    }
    routes
}

fn route_cost_preview(rb: &RouteBuilder, e: &binary::IsogenyExtensionEstimate, field_bits: f64) -> f64 {
    let genus = binary::generic_genus(e.magic_bound);
    let log2_genus = log2_biguint(&genus).max(0.0);
    let dest = rb.model.jacobian_destination(rb.subgroup_bits, log2_genus, field_bits);
    binary::log2_total(&[e.log2_search_cost, dest.log2_cost])
}

fn prime_routes(rb: &RouteBuilder, p: &BigUint) -> Vec<TransferRoute> {
    let p_bits = log2_biguint(p);
    let mut routes = Vec::new();
    if p > &BigUint::from(3u32) {
        routes.push(rb.jacobian_route(
            "same_field_cover",
            "same-field genus-2 cover H -> E of degree 2 over GF(p) (prime_quadratic_pullback_v1)".into(),
            Evidence::Construction,
            inv(vec![("map_degree", Value::from(2u32)), ("construction", Value::from("docs/curves/COVERS.md prime characteristic > 3"))]),
            &BigUint::from(2u32),
            p_bits,
            Transport::CarriedIfCoprimeToDegree,
            f64::NEG_INFINITY,
            vec!["the Jacobian over the same field is larger than E: its best solver is rho or q^(2-2/g) index calculus, both at or above rho on E".into()],
        ));
    }
    for k in [2u32, 3] {
        routes.push(rb.jacobian_route(
            "base_change_descent",
            format!("base change to GF(p^{k}) then descent to a genus-{k} cover over GF(p)"),
            Evidence::Proposed,
            inv(vec![("extension_degree", Value::from(k)), ("genus_lower_bound", Value::from(k))]),
            &BigUint::from(k),
            p_bits,
            Transport::SameGroup,
            f64::NEG_INFINITY,
            vec![
                "GF(p) has no proper subfield; any descent must first base-change, and a cover over GF(p) of the Weil restriction has genus >= k (existence for this curve not established)".into(),
                "priced at the most favourable genus the literature allows, g = k, so the row is a lower bound on the route".into(),
            ],
        ));
    }
    routes
}

fn extension_routes(rb: &RouteBuilder, c: &RegistryCurve, p: &BigUint, k: u32) -> Vec<TransferRoute> {
    let p_bits = log2_biguint(p);
    let two_torsion = (&c.order % 2u32).is_zero();
    let mut routes = Vec::new();
    // Gaudry–Diem decomposition on E(F_{p^k}) with a factor base from F_p.
    let dest_log2 = rb.model.log2_extension_ic(p_bits, k);
    let considered = vec![(DestinationSolver::ExtensionIndexCalculus, dest_log2)];
    routes.push(rb.finish(
        "extension_index_calculus",
        format!("Gaudry–Diem decomposition index calculus on E(GF(p^{k})) with factor base from GF(p)"),
        Evidence::Construction,
        inv(vec![
            ("extension_degree", Value::from(k)),
            ("relations_exponent", Value::from(round_to(2.0 - 2.0 / f64::from(k), 3))),
            ("decomposition_bezout_log2", Value::from(k * (k - 1))),
        ]),
        Destination {
            object: format!("E(GF(p^{k}))"),
            genus: "1".into(),
            log2_genus: 0.0,
            field_bits: round_to(p_bits * f64::from(k), 3),
            log2_order: round_to(p_bits * f64::from(k), 3),
        },
        Transport::SameGroup,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        DestinationCost {
            solver: DestinationSolver::ExtensionIndexCalculus,
            log2_cost: dest_log2,
            log2_leading: (2.0 - 2.0 / f64::from(k)) * p_bits + f64::from(k * (k - 1)),
            considered,
            enge_gaudry_hypothesis_holds: false,
        },
        vec![
            "no transfer: the method works on E itself with a subfield factor base; the decomposition price is the Bézout bound of the Semaev system, a stand-in for the Gröbner solve".into(),
            "the repository's gaudry_cubic measurements at k = 3 found the solve constant far above this formula; a formula advantage here is a lead, not a measurement".into(),
        ],
    ));
    // Descent to a genus-k cover over GF(p).
    let genus = BigUint::from(k);
    let log2_genus = log2_biguint(&genus);
    let transport = rb.transport_by_size(log2_genus, p_bits);
    routes.push(rb.jacobian_route(
        "cover_descent",
        format!("descent to a genus-{k} cover C/GF(p) of E/GF(p^{k})"),
        Evidence::Proposed,
        inv(vec![
            ("extension_degree", Value::from(k)),
            ("genus_lower_bound", Value::from(k)),
            ("rational_two_torsion_point", Value::from(two_torsion)),
        ]),
        &genus,
        p_bits,
        transport,
        f64::NEG_INFINITY,
        vec![
            "existence of a genus-k cover over GF(p) is not established for this curve; odd-characteristic GHS and Diem's genus-3 covers need extra structure (for k = 3, rational 2-torsion is the usual precondition, recorded as an exact invariant)".into(),
            "priced at the most favourable genus, g = k, so the row is a lower bound on the route".into(),
        ],
    ));
    routes
}

/// For each GHS tower priced at both magic-number bounds, the pair's
/// classes are replaced by `unresolved` when they disagree on whether the
/// route beats rho; when they agree, each keeps its own class.
fn resolve_bound_pairs(routes: &mut [TransferRoute]) {
    let favourable = |c: RouteClass| matches!(c, RouteClass::FormulaAdvantage | RouteClass::FormulaMarginal);
    let pairs: Vec<(usize, usize)> = routes
        .iter()
        .enumerate()
        .filter(|(_, r)| r.invariants.get("magic_number_bound").and_then(Value::as_str) == Some("lower"))
        .filter_map(|(i, r)| {
            let l = r.invariants.get("base_degree")?.clone();
            let j = routes.iter().position(|x| {
                x.kind == r.kind
                    && x.invariants.get("base_degree") == Some(&l)
                    && x.invariants.get("magic_number_bound").and_then(Value::as_str) == Some("upper")
            })?;
            Some((i, j))
        })
        .collect();
    for (i, j) in pairs {
        if favourable(routes[i].class) != favourable(routes[j].class) {
            let (win, lose) = if favourable(routes[i].class) { (i, j) } else { (j, i) };
            let other = class_name(routes[lose].class);
            routes[win].class = RouteClass::Unresolved;
            routes[win].notes.push(format!(
                "unresolved: this magic-number bound beats rho and the other bound is {other}; compute the reduced a (Hess) to decide"
            ));
            routes[lose].notes.push(
                "the other magic-number bound of this tower beats rho; the pair is unresolved".into(),
            );
        }
    }
}

/// Screen one registry curve.
pub fn screen(c: &RegistryCurve, model: &CostModel, budget: u64) -> CurveTransfers {
    let field = field_info(&c.field);
    let subgroup = subgroup_info(c, budget);
    let reference_log2 = model.log2_rho_elliptic(subgroup.bits, field.log2_size);
    let rb = RouteBuilder {
        model,
        subgroup_bits: subgroup.bits,
        reference_log2,
        source_field_bits: field.log2_size,
    };
    let routes = match &c.field {
        Field::Binary { n, irr, a, b } => binary_routes(&rb, c, *n, irr, a, b),
        Field::Prime { p } => prime_routes(&rb, p),
        Field::Extension { p, k } => extension_routes(&rb, c, p, *k),
    };
    let mut routes = routes;
    resolve_bound_pairs(&mut routes);
    let best_route = routes
        .iter()
        .enumerate()
        .filter(|(_, r)| r.class != RouteClass::NotCarried)
        .min_by(|(_, x), (_, y)| {
            x.cost
                .log2_total
                .partial_cmp(&y.cost.log2_total)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i);
    let verdict = routes
        .iter()
        .map(|r| r.class)
        .min()
        .unwrap_or(RouteClass::NotCarried);
    CurveTransfers {
        slug: c.slug.clone(),
        icv1: c.icv1.clone(),
        standard_names: c.standard_names.clone(),
        family: c.family.clone(),
        field,
        trace: c.trace.to_string(),
        subgroup,
        reference_log2: round_to(reference_log2, 3),
        verdict,
        best_route,
        routes,
    }
}

impl TransferFile {
    pub fn build(registry_text: &str, model: CostModel, budget: u64) -> Result<Self, String> {
        use rayon::prelude::*;
        let curves = read_registry(registry_text)?;
        let records: Vec<CurveTransfers> = curves.par_iter().map(|c| screen(c, &model, budget)).collect();
        let digest = crate::hash::sha256::sha256(registry_text.as_bytes());
        let mut by_class = BTreeMap::new();
        let mut by_route_kind = BTreeMap::new();
        let mut formula_advantage = Vec::new();
        let mut unresolved = Vec::new();
        for r in &records {
            *by_class.entry(class_name(r.verdict)).or_insert(0) += 1;
            for route in &r.routes {
                *by_route_kind.entry(route.kind.clone()).or_insert(0) += 1;
            }
            if matches!(r.verdict, RouteClass::FormulaAdvantage | RouteClass::FormulaMarginal) {
                let best = r.best_route.map(|i| r.routes[i].label.clone()).unwrap_or_default();
                formula_advantage.push((r.slug.clone(), best));
            } else if r.verdict == RouteClass::Unresolved {
                let best = r.best_route.map(|i| r.routes[i].label.clone()).unwrap_or_default();
                unresolved.push((r.slug.clone(), best));
            }
        }
        Ok(TransferFile {
            schema: SCHEMA.into(),
            generated_by: "cargo run --release --bin jacobian_transfers -- build".into(),
            what_this_is: "Every known route from each registered curve's discrete logarithm into a \
                           Jacobian, priced by one frozen formula-level cost model against Pollard rho \
                           on the same subgroup; see docs/curves/JACOBIAN_TRANSFERS.md. Rows are \
                           formulas and leads, never measurements or security claims."
                .into(),
            registry_sha256: digest.iter().map(|b| format!("{b:02x}")).collect(),
            factor_budget: budget,
            cost_model: model,
            summary: Summary {
                curves: records.len(),
                routes: records.iter().map(|r| r.routes.len()).sum(),
                by_class,
                by_route_kind,
                formula_advantage,
                unresolved,
            },
            curves: records,
        })
    }

    /// Pretty JSON with a trailing newline, as committed.
    pub fn to_json(&self) -> String {
        let mut s = serde_json::to_string_pretty(self).expect("serialisable");
        s.push('\n');
        s
    }

    pub fn find(&self, name: &str) -> Option<&CurveTransfers> {
        let slug = crate::cryptanalysis::curve_id::resolve(name).unwrap_or(name);
        self.curves
            .iter()
            .find(|c| c.slug == slug || c.standard_names.iter().any(|s| s == name))
    }
}

pub fn class_name(c: RouteClass) -> String {
    serde_json::to_value(c)
        .ok()
        .and_then(|v| v.as_str().map(str::to_string))
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    const REGISTRY: &str = include_str!("../../../docs/curves/registry.json");

    fn registry() -> Vec<RegistryCurve> {
        read_registry(REGISTRY).unwrap()
    }

    fn by_name(name: &str) -> RegistryCurve {
        let slug = crate::cryptanalysis::curve_id::resolve(name).unwrap_or(name);
        registry()
            .into_iter()
            .find(|c| c.slug == slug || c.standard_names.iter().any(|s| s == name))
            .unwrap_or_else(|| panic!("{name} not registered"))
    }

    #[test]
    fn every_registry_entry_parses() {
        let curves = registry();
        assert!(curves.len() >= 300);
        assert!(curves.iter().any(|c| matches!(c.field, Field::Extension { .. })));
        assert!(curves.iter().any(|c| matches!(c.field, Field::Binary { .. })));
        assert!(curves.iter().any(|c| matches!(c.field, Field::Prime { .. })));
    }

    #[test]
    fn oakley_group_three_reproduces_the_cover_route_census() {
        // research/cover_route_prime_binary_20261008: m = 5 over GF(2^31),
        // m = 26 over GF(2^5), m = 126 over GF(2).
        let c = by_name("Oakley Group 3");
        let r = screen(&c, &CostModel::default(), DEFAULT_FACTOR_BUDGET);
        let magic = |l: u32| {
            r.routes
                .iter()
                .find(|x| x.kind == "ghs_weil_descent" && x.invariants["base_degree"] == l)
                .unwrap()
                .invariants["magic_number"]
                .as_u64()
                .unwrap()
        };
        assert_eq!(magic(31), 5);
        assert_eq!(magic(5), 26);
        assert_eq!(magic(1), 126);
        let g16 = r
            .routes
            .iter()
            .find(|x| x.kind == "ghs_weil_descent" && x.invariants["base_degree"] == 31)
            .unwrap();
        // Type I: genus 2^4 - 1.
        assert_eq!(g16.destination.genus, "15");
        assert_eq!(g16.transport, Transport::SizeAdmits);
        // Genus 15 over GF(2^31): relation collection 15! * 2^31 = 2^71.25
        // against rho's 2^76 leaves an exponent ratio near 4.75 bits, which the
        // genus constants (10 g^2 per Jacobian operation) consume.
        assert!(g16.cost.ratio_log2 < cost::CostModel::default().advantage_margin_bits);
        assert_ne!(g16.class, RouteClass::FormulaAdvantage);
        assert_ne!(r.verdict, RouteClass::FormulaAdvantage);
    }

    #[test]
    fn unresolved_trace_of_a_prices_both_magic_bounds() {
        // c2pnb304w1: Tr(a) = 1 with n = 152 even over GF(4); magic_b = 8.
        let r = screen(&by_name("c2pnb304w1"), &CostModel::default(), DEFAULT_FACTOR_BUDGET);
        let rows: Vec<_> = r
            .routes
            .iter()
            .filter(|x| x.kind == "ghs_weil_descent" && x.invariants["base_degree"] == 2)
            .collect();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].invariants["magic_number"], 8);
        assert_eq!(rows[1].invariants["magic_number"], 9);
        assert!(rows.iter().all(|x| x.evidence == Evidence::Heuristic));
        // The lower bound cannot carry the 288-bit subgroup; the upper bound
        // beats rho in formula terms, so the pair is unresolved.
        assert_eq!(rows[0].transport, Transport::TooSmall);
        assert_eq!(rows[0].class, RouteClass::NotCarried);
        assert_eq!(rows[1].class, RouteClass::Unresolved);
        assert_eq!(r.verdict, RouteClass::Unresolved);
        assert_eq!(&r.routes[r.best_route.unwrap()], rows[1]);
        // c2pnb368w1 has Tr(a) = 0: one exact row per tower.
        let r = screen(&by_name("c2pnb368w1"), &CostModel::default(), DEFAULT_FACTOR_BUDGET);
        let g7 = r
            .routes
            .iter()
            .filter(|x| x.kind == "ghs_weil_descent" && x.invariants["base_degree"] == 92)
            .collect::<Vec<_>>();
        assert_eq!(g7.len(), 1);
        assert_eq!(g7[0].evidence, Evidence::ExactInvariant);
        assert_eq!(g7[0].destination.genus, "7");
    }

    #[test]
    fn koblitz_curves_degenerate_and_random_b_sits_at_the_maximum() {
        let k = screen(&by_name("sect163k1"), &CostModel::default(), DEFAULT_FACTOR_BUDGET);
        let ghs: Vec<_> = k.routes.iter().filter(|x| x.kind == "ghs_weil_descent").collect();
        assert_eq!(ghs.len(), 1);
        assert_eq!(ghs[0].invariants["magic_number"], 1);
        assert_eq!(ghs[0].transport, Transport::Degenerate);
        assert_eq!(k.verdict, RouteClass::Boundary);
        let b = screen(&by_name("sect113r1"), &CostModel::default(), DEFAULT_FACTOR_BUDGET);
        let ghs = b.routes.iter().find(|x| x.kind == "ghs_weil_descent").unwrap();
        assert_eq!(ghs.invariants["magic_number"], 113);
        let two_112 = BigUint::from(1u32) << 112u32;
        assert!(ghs.destination.genus == two_112.to_string()
            || ghs.destination.genus == (&two_112 - BigUint::from(1u32)).to_string());
        assert_eq!(b.verdict, RouteClass::Boundary);
    }

    #[test]
    fn prime_field_curves_are_boundary_rows_only() {
        for name in ["P-256", "secp256k1", "Curve25519"] {
            let r = screen(&by_name(name), &CostModel::default(), DEFAULT_FACTOR_BUDGET);
            assert_eq!(r.verdict, RouteClass::Boundary, "{name}");
            assert!(r.routes.iter().all(|x| x.cost.ratio_log2 < 0.0), "{name}");
            assert!(r.routes.iter().any(|x| x.kind == "same_field_cover"));
            assert_eq!(r.routes.iter().filter(|x| x.kind == "base_change_descent").count(), 2);
        }
    }

    #[test]
    fn extension_curves_get_both_extension_routes() {
        let c = registry()
            .into_iter()
            .find(|c| matches!(c.field, Field::Extension { .. }))
            .unwrap();
        let r = screen(&c, &CostModel::default(), DEFAULT_FACTOR_BUDGET);
        assert_eq!(r.field.kind, "extension");
        assert!(r.routes.iter().any(|x| x.kind == "extension_index_calculus"));
        assert!(r.routes.iter().any(|x| x.kind == "cover_descent"));
    }

    #[test]
    fn unrecorded_subgroups_fall_back_to_the_largest_prime() {
        let c = registry()
            .into_iter()
            .find(|c| c.subgroup_order.is_none())
            .unwrap();
        let r = screen(&c, &CostModel::default(), DEFAULT_FACTOR_BUDGET);
        assert_ne!(r.subgroup.status, "registry_representation");
        assert!(r.subgroup.bits > 0.0);
    }

    #[test]
    fn the_whole_registry_builds_and_round_trips() {
        let file = TransferFile::build(REGISTRY, CostModel::default(), DEFAULT_FACTOR_BUDGET).unwrap();
        assert_eq!(file.curves.len(), registry().len());
        let json = file.to_json();
        assert!(!json.contains("null,") || !json.contains("\"log2_cost\": null"), "non-finite cost");
        let back: TransferFile = serde_json::from_str(&json).unwrap();
        // The file parses back; values agree to 12 significant digits (the
        // JSON parser is not correctly rounded above 2^53, and `check` rebuilds
        // from the registry rather than round-tripping these numbers).
        assert_eq!(back.summary, file.summary);
        assert_eq!(back.curves.len(), file.curves.len());
        for (x, y) in back.curves.iter().zip(&file.curves) {
            assert_eq!(x.slug, y.slug);
            assert_eq!(x.verdict, y.verdict);
            assert_eq!(x.best_route, y.best_route);
            for (rx, ry) in x.routes.iter().zip(&y.routes) {
                let (a, b) = (rx.cost.log2_total, ry.cost.log2_total);
                assert!((a - b).abs() <= 1e-11 * a.abs().max(1.0), "{} {a} {b}", x.slug);
            }
        }
        // Rebuilding is deterministic.
        let again = TransferFile::build(REGISTRY, CostModel::default(), DEFAULT_FACTOR_BUDGET).unwrap();
        assert_eq!(again.to_json(), json);
        assert_eq!(file.summary.curves, file.curves.len());
        assert!(file.find("Oakley Group 3").is_some());
    }
}
