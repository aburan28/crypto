//! # ECDLP across every NIST curve: prime, binary random, Koblitz, and the
//! anomalous branch.
//!
//! One solver front end for the fifteen FIPS 186-4 curves — `P-192 … P-521`
//! over `F_p`, `B-163 … B-571` (random) and `K-163 … K-571` (Koblitz, the
//! "anomalous binary curves") over `F_{2^m}` — plus the structural attacks
//! that would apply to a curve of each kind if its parameters allowed them:
//!
//! | branch | when | implementation |
//! |---|---|---|
//! | Smart–Semaev–Satoh–Araki | `#E(F_p) = p` (trace 1) | [`anomalous`] — `p`-adic lift, formal log, polynomial time |
//! | Pohlig–Hellman | generator order composite | [`pohlig`] |
//! | MOV / Frey–Rück | embedding degree small | detected in the audit; the transfer itself is in `mov_attack` (supersingular `k = 2`) |
//! | BSGS / kangaroo | scalar known to lie in an interval | [`interval`] |
//! | Pollard rho | otherwise | [`rho`] — parallel, distinguished points, negation map, Frobenius folding on Koblitz curves |
//!
//! The audit ([`audit`]) runs each test against every curve.  None of the
//! fifteen is anomalous, has a usable embedding degree, or a composite
//! generator order, and every binary extension degree is prime (no
//! Weil-descent route), so on the real parameters the only branch that
//! runs is the generic one, whose expected cost — `√(πn/4)` additions with
//! the negation map, a further `√m` lower with Frobenius folding on
//! `K-m` — is what the audit reports.  The solvers nonetheless run *on the
//! real curves*: a planted scalar in an interval of `2^20` is recovered on
//! each of the fifteen in the tests, which exercises the full-width
//! arithmetic and the coefficient bookkeeping; the whole-group rho and the
//! Smart attack are demonstrated on small or constructed curves
//! ([`toy`], [`anomalous::generate`]).
//!
//! **Claim hygiene.**  Nothing here recovers a full-width key on a NIST
//! curve, and the reports say so: a solve without an interval on a real
//! curve is refused as infeasible unless forced, and then runs only to its
//! iteration budget.  Costs are stated as expected group additions, never
//! as security claims about deployed parameters.

use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use serde::Serialize;
use std::fmt;
use std::str::FromStr;

pub mod anomalous;
pub mod group;
pub mod interval;
pub mod pohlig;
pub mod rho;
pub mod toy;

pub use anomalous::{smart_attack, SmartReport};
pub use group::{BinaryGroup, EcdlpGroup, PrimeGroup};
pub use interval::{bsgs, kangaroo, IntervalReport, KangarooOptions};
pub use pohlig::{pohlig_hellman, PohligHellmanReport};
pub use rho::{pollard_rho, RhoOptions, RhoReport};

use crate::binary_ecc::curve::BinaryPoint;
use crate::cryptanalysis::curve_catalog::{self, CurveObject};
use crate::cryptanalysis::curve_traits::arith::{factor, is_prime};
use crate::cryptanalysis::mov_attack::embedding_degree;
use crate::ecc::point::Point;

/// `#E_a(F_{2^m})` for the Koblitz curve `y² + xy = x³ + ax² + 1`, from
/// the Lucas recurrence `V_0 = 2, V_1 = t, V_{k+1} = t V_k − 2 V_{k−1}`
/// with `t = (−1)^{1−a}` the trace over `F_2`: `#E = 2^m + 1 − V_m`.
pub fn koblitz_group_order(a: u8, m: u32) -> BigUint {
    let t: i128 = if a == 0 { -1 } else { 1 };
    let (mut v0, mut v1) = (BigInt::from(2), BigInt::from(t));
    for _ in 1..m {
        let next = BigInt::from(t) * &v1 - BigInt::from(2) * &v0;
        v0 = v1;
        v1 = next;
    }
    let order = (BigInt::one() << m) + BigInt::one() - v1;
    order.to_biguint().expect("positive group order")
}

/// The three NIST families.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Family {
    /// `P-xxx`: short Weierstrass over `F_p`.
    Prime,
    /// `K-xxx`: Koblitz (anomalous binary) curve over `F_{2^m}`.
    Koblitz,
    /// `B-xxx`: pseudo-random binary curve over `F_{2^m}`.
    BinaryRandom,
}

impl Family {
    pub fn tag(self) -> &'static str {
        match self {
            Family::Prime => "prime",
            Family::Koblitz => "koblitz",
            Family::BinaryRandom => "binary",
        }
    }
}

/// The NIST name and the `curve_catalog` entry that carries its constants.
/// `B-163` is `sect163r2` (not `r1`), per FIPS 186-4 D.1.3.1.
pub const NIST_CURVES: &[(&str, &str, Family)] = &[
    ("P-192", "p192", Family::Prime),
    ("P-224", "p224", Family::Prime),
    ("P-256", "p256", Family::Prime),
    ("P-384", "p384", Family::Prime),
    ("P-521", "p521", Family::Prime),
    ("K-163", "sect163k1", Family::Koblitz),
    ("K-233", "sect233k1", Family::Koblitz),
    ("K-283", "sect283k1", Family::Koblitz),
    ("K-409", "sect409k1", Family::Koblitz),
    ("K-571", "sect571k1", Family::Koblitz),
    ("B-163", "sect163r2", Family::BinaryRandom),
    ("B-233", "sect233r1", Family::BinaryRandom),
    ("B-283", "sect283r1", Family::BinaryRandom),
    ("B-409", "sect409r1", Family::BinaryRandom),
    ("B-571", "sect571r1", Family::BinaryRandom),
];

/// The concrete group behind a [`Curve`].
#[derive(Clone)]
pub enum Group {
    /// Both variants are boxed: a curve carries several multi-word
    /// constants, and the enum is held by value in [`Curve`].
    Prime(Box<PrimeGroup>),
    Binary(Box<BinaryGroup>),
}

/// A point on either kind of curve.
#[derive(Clone, Debug, PartialEq)]
pub enum AnyPoint {
    Prime(Point),
    Binary(BinaryPoint),
}

/// A curve the solver runs on: one of the fifteen NIST curves, a toy
/// curve, or a constructed anomalous curve.
#[derive(Clone)]
pub struct Curve {
    /// Display name (`P-256`, `K-163`, `toy-k23a1`, `anomalous-128-3`).
    pub name: String,
    /// `curve_catalog` name when the curve comes from the catalog.
    pub catalog_name: Option<&'static str>,
    pub family: Family,
    /// Where the parameters come from.
    pub standard: String,
    pub group: Group,
    /// `⌈log₂ q⌉` of the base field.
    pub field_bits: u64,
    /// `q`: `p`, or `2^m`.
    pub field_size: BigUint,
    /// `#E(F_q) = h · n`.
    pub group_order: BigUint,
    pub cofactor: BigUint,
}

impl Curve {
    fn from_catalog(
        nist_name: &str,
        catalog_name: &'static str,
        family: Family,
    ) -> Result<Self, String> {
        let entry = curve_catalog::by_name(catalog_name)
            .ok_or_else(|| format!("{catalog_name} is not in the curve catalog"))?;
        let field = entry.field();
        let cofactor = entry.cofactor();
        let group_order = entry.group_order();
        let (group, field_size) = match entry.object {
            CurveObject::Prime(cp) => {
                let fs = cp.p.clone();
                (Group::Prime(Box::new(PrimeGroup::new(cp))), fs)
            }
            CurveObject::Binary(bc) => {
                let fs = BigUint::one() << bc.m;
                (
                    Group::Binary(Box::new(BinaryGroup::new(nist_name, *bc)?)),
                    fs,
                )
            }
            CurveObject::Char3(_) => {
                return Err("characteristic-three curves are out of scope".into())
            }
        };
        Ok(Curve {
            name: nist_name.to_string(),
            catalog_name: Some(catalog_name),
            family,
            standard: entry.standard.to_string(),
            group,
            field_bits: field.bits,
            field_size,
            group_order,
            cofactor,
        })
    }

    fn from_prime_group(name: &str, standard: &str, g: PrimeGroup) -> Self {
        let c = &g.curve;
        Curve {
            name: name.to_string(),
            catalog_name: None,
            family: Family::Prime,
            standard: standard.to_string(),
            field_bits: c.p.bits(),
            field_size: c.p.clone(),
            group_order: &c.n * BigUint::from(c.h),
            cofactor: BigUint::from(c.h),
            group: Group::Prime(Box::new(g)),
        }
    }

    fn from_binary_group(name: &str, standard: &str, g: BinaryGroup) -> Self {
        let c = &g.curve;
        let family =
            if g.frobenius.is_some() || (c.b == crate::binary_ecc::f2m::F2mElement::one(c.m)) {
                Family::Koblitz
            } else {
                Family::BinaryRandom
            };
        Curve {
            name: name.to_string(),
            catalog_name: None,
            family,
            standard: standard.to_string(),
            field_bits: c.m as u64,
            field_size: BigUint::one() << c.m,
            group_order: &c.order * &c.cofactor,
            cofactor: c.cofactor.clone(),
            group: Group::Binary(Box::new(g)),
        }
    }

    /// The generator's (prime, on NIST curves) order `n`.
    pub fn order(&self) -> &BigUint {
        match &self.group {
            Group::Prime(g) => g.order(),
            Group::Binary(g) => g.order(),
        }
    }

    pub fn generator(&self) -> AnyPoint {
        match &self.group {
            Group::Prime(g) => AnyPoint::Prime(g.generator()),
            Group::Binary(g) => AnyPoint::Binary(g.generator()),
        }
    }

    pub fn is_identity(&self, p: &AnyPoint) -> bool {
        match p {
            AnyPoint::Prime(p) => matches!(p, Point::Infinity),
            AnyPoint::Binary(p) => matches!(p, BinaryPoint::Infinity),
        }
    }

    /// `[k]P`.
    pub fn mul(&self, p: &AnyPoint, k: &BigUint) -> Result<AnyPoint, String> {
        match (&self.group, p) {
            (Group::Prime(g), AnyPoint::Prime(p)) => Ok(AnyPoint::Prime(g.mul(p, k))),
            (Group::Binary(g), AnyPoint::Binary(p)) => Ok(AnyPoint::Binary(g.mul(p, k))),
            _ => Err("point does not belong to this curve's field".into()),
        }
    }

    /// `Q = [k]G`: the known-answer target for a planted instance.
    pub fn plant(&self, k: &BigUint) -> AnyPoint {
        self.mul(&self.generator(), k).expect("generator matches")
    }

    /// Parse hex coordinates (with or without `0x`) into a point of this
    /// curve, checking the curve equation.
    pub fn parse_point(&self, x_hex: &str, y_hex: &str) -> Result<AnyPoint, String> {
        let parse = |s: &str| {
            let s = s.trim().trim_start_matches("0x").trim_start_matches("0X");
            BigUint::parse_bytes(s.as_bytes(), 16)
                .ok_or_else(|| format!("bad hex coordinate {s:?}"))
        };
        let (x, y) = (parse(x_hex)?, parse(y_hex)?);
        let pt = match &self.group {
            Group::Prime(g) => AnyPoint::Prime(g.point(&x, &y)),
            Group::Binary(g) => AnyPoint::Binary(g.point(&x, &y)),
        };
        self.validate_target(&pt)?;
        Ok(pt)
    }

    /// Hex coordinates of a point; `None` for the identity.
    pub fn point_hex(&self, p: &AnyPoint) -> Option<(String, String)> {
        match (&self.group, p) {
            (Group::Prime(g), AnyPoint::Prime(p)) => g.coords_hex(p),
            (Group::Binary(g), AnyPoint::Binary(p)) => g.coords_hex(p),
            _ => None,
        }
    }

    /// On the curve and in the generator's subgroup (`[n]Q = O`).
    pub fn validate_target(&self, p: &AnyPoint) -> Result<(), String> {
        let on_curve = match (&self.group, p) {
            (Group::Prime(g), AnyPoint::Prime(p)) => g.is_on_curve(p),
            (Group::Binary(g), AnyPoint::Binary(p)) => g.is_on_curve(p),
            _ => return Err("point does not belong to this curve's field".into()),
        };
        if !on_curve {
            return Err("point is not on the curve".into());
        }
        if !self.cofactor.is_one() {
            let n_q = self.mul(p, self.order())?;
            if !self.is_identity(&n_q) {
                return Err(format!(
                    "point is not in the order-{} subgroup (cofactor {})",
                    self.order(),
                    self.cofactor
                ));
            }
        }
        Ok(())
    }

    /// Machine-readable parameters.
    pub fn describe(&self) -> serde_json::Value {
        let (gx, gy) = self.point_hex(&self.generator()).unwrap_or_default();
        serde_json::json!({
            "name": self.name,
            "catalog_name": self.catalog_name,
            "family": self.family.tag(),
            "standard": self.standard,
            "field_bits": self.field_bits,
            "order": self.order().to_string(),
            "order_bits": self.order().bits(),
            "cofactor": self.cofactor.to_string(),
            "group_order": self.group_order.to_string(),
            "generator": {"x": gx, "y": gy},
        })
    }
}

/// The fifteen FIPS 186-4 curves, in the order of [`NIST_CURVES`].
pub fn nist_curves() -> Vec<Curve> {
    NIST_CURVES
        .iter()
        .map(|(nist, cat, fam)| {
            Curve::from_catalog(nist, cat, *fam).unwrap_or_else(|e| panic!("{nist}: {e}"))
        })
        .collect()
}

/// Every name [`curve_by_name`] accepts, for listings.
pub fn curve_names() -> Vec<String> {
    let mut v: Vec<String> = NIST_CURVES.iter().map(|(n, _, _)| n.to_string()).collect();
    v.extend(toy::toy_names().into_iter().map(String::from));
    v.push("anomalous-<bits>[-<seed>]".to_string());
    v
}

/// Resolve a curve by NIST name (`P-256`, `k163`, `B-571`), catalog
/// alias (`secp256r1`, `sect233k1`), toy name (`toy-k23a1`,
/// `toy-p10039`, `toy-k17a0-full`) or `anomalous-<bits>[-<seed>]`.
pub fn curve_by_name(name: &str) -> Result<Curve, String> {
    let q = name.trim().to_ascii_lowercase().replace('_', "-");
    // NIST spellings: "p-256", "p256", "nistp256", "k-163", "b163".
    let compact = q.replace('-', "");
    for (nist, cat, fam) in NIST_CURVES {
        let nist_l = nist.to_ascii_lowercase();
        if q == nist_l
            || compact == nist_l.replace('-', "")
            || compact == format!("nist{}", nist_l.replace('-', ""))
            || q == *cat
        {
            return Curve::from_catalog(nist, cat, *fam);
        }
    }
    // Catalog aliases (secp256r1, prime256v1, …) that map onto a NIST curve.
    if let Some(entry) = curve_catalog::by_name(&q) {
        if let Some((nist, cat, fam)) = NIST_CURVES.iter().find(|(_, c, _)| *c == entry.name) {
            return Curve::from_catalog(nist, cat, *fam);
        }
        return Err(format!(
            "{name} is in the catalog but is not a NIST curve; this solver covers {}",
            NIST_CURVES
                .iter()
                .map(|(n, _, _)| *n)
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    match q.as_str() {
        "toy-p10039" => return Ok(Curve::from_prime_group(&q, "toy", toy::prime_small())),
        "toy-p98893" => return Ok(Curve::from_prime_group(&q, "toy", toy::prime_mid())),
        _ => {}
    }
    if let Some(rest) = q.strip_prefix("toy-k") {
        let (body, full) = match rest.strip_suffix("-full") {
            Some(b) => (b, true),
            None => (rest, false),
        };
        let (m, a) = body
            .split_once('a')
            .ok_or_else(|| format!("toy Koblitz name must look like toy-k23a1, got {name}"))?;
        let m: u32 = m.parse().map_err(|_| format!("bad m in {name}"))?;
        let a: u8 = a.parse().map_err(|_| format!("bad a in {name}"))?;
        let g = if full {
            toy::koblitz_full_group(m, a)?
        } else {
            toy::koblitz(m, a)?
        };
        return Ok(Curve::from_binary_group(&q, "toy", g));
    }
    if let Some(rest) = q.strip_prefix("anomalous-") {
        let (bits, seed) = match rest.split_once('-') {
            Some((b, s)) => (
                b,
                s.parse::<u64>()
                    .map_err(|_| format!("bad seed in {name}"))?,
            ),
            None => (rest, 0u64),
        };
        let bits: u32 = bits
            .parse()
            .map_err(|_| format!("bad bit size in {name}"))?;
        let cp = anomalous::generate(bits, seed)?;
        let g = PrimeGroup::new(cp);
        return Ok(Curve::from_prime_group(
            &format!("anomalous-{bits}-{seed}"),
            "constructed (CM, trace 1)",
            g,
        ));
    }
    Err(format!(
        "unknown curve {name:?}; known: {}",
        curve_names().join(", ")
    ))
}

// ── Audit ───────────────────────────────────────────────────────────────────

/// The solver branch a curve calls for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Method {
    /// Let the audit decide.
    Auto,
    /// Smart–Semaev–Satoh–Araki (anomalous prime curve).
    Smart,
    /// Pohlig–Hellman (composite generator order).
    PohligHellman,
    /// Baby-step giant-step on an interval.
    Bsgs,
    /// Parallel kangaroos on an interval.
    Kangaroo,
    /// Parallel rho on the whole group.
    Rho,
}

impl Method {
    pub fn tag(self) -> &'static str {
        match self {
            Method::Auto => "auto",
            Method::Smart => "smart",
            Method::PohligHellman => "pohlig-hellman",
            Method::Bsgs => "bsgs",
            Method::Kangaroo => "kangaroo",
            Method::Rho => "rho",
        }
    }
}

impl FromStr for Method {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        Ok(match s.to_ascii_lowercase().as_str() {
            "auto" => Method::Auto,
            "smart" | "anomalous" => Method::Smart,
            "pohlig-hellman" | "pohlig_hellman" | "ph" => Method::PohligHellman,
            "bsgs" => Method::Bsgs,
            "kangaroo" | "lambda" => Method::Kangaroo,
            "rho" => Method::Rho,
            _ => {
                return Err(format!(
                    "unknown method {s:?}; one of auto, smart, pohlig-hellman, bsgs, kangaroo, rho"
                ))
            }
        })
    }
}

impl fmt::Display for Method {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.tag())
    }
}

/// Embedding degrees up to this bound are searched; every NIST curve has
/// none below it (and in fact none below ~`n/…`, the audit just certifies
/// the bound).
pub const EMBEDDING_DEGREE_BOUND: u32 = 100;

/// Structural facts about a curve that decide which attack applies.
#[derive(Clone, Debug, Serialize)]
pub struct Audit {
    pub name: String,
    pub family: Family,
    pub standard: String,
    pub field_bits: u64,
    pub order_bits: u64,
    pub cofactor: String,
    /// `t = q + 1 − #E(F_q)`.
    pub trace: String,
    /// `#E(F_q) = q`, i.e. `t = 1`: Smart's attack applies.
    pub anomalous: bool,
    /// Miller–Rabin on the generator's order.
    pub order_is_prime: bool,
    /// `(q, e)` of the generator's order when composite.
    pub order_factors: Vec<(String, u32)>,
    /// Smallest `k ≤ bound` with `n | q^k − 1`, if any: MOV / Frey–Rück.
    pub embedding_degree: Option<u32>,
    pub embedding_degree_bound: u32,
    /// Binary curves: whether `m` is prime (composite `m` opens Weil
    /// descent / GHS).  `None` on prime fields.
    pub extension_degree_prime: Option<bool>,
    /// Koblitz: whether `#E_a(F_{2^m})` from the Lucas recurrence equals
    /// the catalog's `h · n`.
    pub koblitz_order_consistent: Option<bool>,
    /// Frobenius eigenvalue on the subgroup, Koblitz only.
    pub frobenius_eigenvalue: Option<String>,
    /// Automorphisms folded by the rho: `2` (negation) or `2m` (Koblitz).
    pub rho_class_size: u32,
    /// `log₂ √(πn / (2·class))`: expected rho additions with folding.
    pub rho_log2_expected: f64,
    /// `log₂ √(πn / 2)`: the same with no folding.
    pub rho_log2_expected_unfolded: f64,
    pub recommended: Method,
    pub notes: Vec<String>,
}

/// Run every structural check on a curve.
pub fn audit(curve: &Curve) -> Audit {
    let n = curve.order();
    let q = &curve.field_size;
    let trace = BigInt::from(q.clone()) + BigInt::one() - BigInt::from(curve.group_order.clone());
    let anomalous = curve.group_order == *q;
    let order_is_prime = is_prime(n);
    let order_factors = if order_is_prime {
        Vec::new()
    } else {
        let f = factor(n, 1 << 20);
        let mut v: Vec<(String, u32)> = f.primes.iter().map(|(p, e)| (p.to_string(), *e)).collect();
        v.extend(
            f.composites
                .iter()
                .map(|(c, e)| (format!("{c} (composite)"), *e)),
        );
        v
    };
    let embedding = embedding_degree(q, n, EMBEDDING_DEGREE_BOUND);
    let mut notes = Vec::new();
    let (extension_degree_prime, koblitz_order_consistent, frobenius_eigenvalue, class) =
        match &curve.group {
            Group::Prime(_) => (None, None, None, 2u32),
            Group::Binary(g) => {
                let m = g.curve.m;
                let m_prime = is_prime(&BigUint::from(m));
                if !m_prime {
                    notes.push(format!("extension degree m = {m} is composite: Weil descent (GHS) must be assessed"));
                }
                let (consistent, eig, class) = match &g.frobenius {
                    Some(fr) => {
                        let a: u8 = if g.curve.a.is_zero() { 0 } else { 1 };
                        let lucas = koblitz_group_order(a, m);
                        let ok = lucas == curve.group_order;
                        if !ok {
                            notes.push(format!(
                            "Koblitz order from the Lucas recurrence ({lucas}) differs from the catalog's h·n ({})",
                            curve.group_order
                        ));
                        }
                        (Some(ok), Some(fr.lambda.to_string()), 2 * m)
                    }
                    None => (None, None, 2),
                };
                (Some(m_prime), consistent, eig, class)
            }
        };
    let log2_expected = |class: u32| rho::expected_rho_iterations(n, class).log2();
    let recommended = if anomalous && matches!(curve.group, Group::Prime(_)) {
        notes.push(
            "trace 1: #E(F_p) = p — Smart's attack solves the ECDLP in polynomial time".into(),
        );
        Method::Smart
    } else if !order_is_prime {
        notes.push(
            "generator order is composite: Pohlig–Hellman reduces to its largest prime factor"
                .into(),
        );
        Method::PohligHellman
    } else {
        Method::Rho
    };
    if let Some(k) = embedding {
        notes.push(format!(
            "embedding degree {k} ≤ {EMBEDDING_DEGREE_BOUND}: MOV/Frey–Rück transfers the ECDLP to F_q^{k}"
        ));
    }
    if anomalous && matches!(curve.group, Group::Binary(_)) {
        notes.push("#E = 2^m on a binary curve (supersingular-type); the p-adic Smart attack is for prime fields".into());
    }
    Audit {
        name: curve.name.clone(),
        family: curve.family,
        standard: curve.standard.clone(),
        field_bits: curve.field_bits,
        order_bits: n.bits(),
        cofactor: curve.cofactor.to_string(),
        trace: trace.to_string(),
        anomalous,
        order_is_prime,
        order_factors,
        embedding_degree: embedding,
        embedding_degree_bound: EMBEDDING_DEGREE_BOUND,
        extension_degree_prime,
        koblitz_order_consistent,
        frobenius_eigenvalue,
        rho_class_size: class,
        rho_log2_expected: log2_expected(class),
        rho_log2_expected_unfolded: log2_expected(1),
        recommended,
        notes,
    }
}

// ── Dispatcher ──────────────────────────────────────────────────────────────

/// How to run [`solve`].
#[derive(Clone, Debug)]
pub struct SolveOptions {
    pub method: Method,
    /// `k ∈ [lo, lo + width)` when known.
    pub interval: Option<(BigUint, BigUint)>,
    pub threads: usize,
    pub seed: u64,
    /// Rho / kangaroo iteration budget; `None` runs to completion.
    pub max_iterations: Option<u64>,
    pub dp_bits: Option<u32>,
    /// Fold negation and Frobenius in the rho.
    pub fold_automorphisms: bool,
    /// Largest `log₂(expected additions)` a whole-group rho may start with
    /// unless [`force`](Self::force) is set.
    pub feasibility_log2: f64,
    /// Start an infeasible rho anyway (it stops at `max_iterations`).
    pub force: bool,
    /// Intervals up to `2^bsgs_max_bits` wide use BSGS, wider ones kangaroo.
    pub bsgs_max_bits: u64,
}

impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            method: Method::Auto,
            interval: None,
            threads: 1,
            seed: 0x5EED,
            max_iterations: None,
            dp_bits: None,
            fold_automorphisms: true,
            feasibility_log2: 40.0,
            force: false,
            bsgs_max_bits: 30,
        }
    }
}

/// The outcome of a solve.
#[derive(Clone, Debug, Serialize)]
pub struct SolveReport {
    pub curve: String,
    pub family: Family,
    /// The branch that ran (never `Auto`).
    pub method: Method,
    /// `k`, decimal, verified by `[k]G = Q`.
    pub scalar: Option<String>,
    pub scalar_hex: Option<String>,
    #[serde(skip)]
    pub scalar_value: Option<BigUint>,
    pub verified: bool,
    /// Group additions performed (rho/kangaroo iterations, BSGS steps,
    /// Pohlig–Hellman inner work); `0` for the Smart attack, which does no
    /// group work on `E(F_p)` beyond verification.
    pub group_ops: u64,
    pub elapsed_ms: u128,
    /// The textbook expectation for the branch that ran.
    pub expected_ops: f64,
    pub table_size: usize,
    pub threads: usize,
    pub notes: Vec<String>,
    pub failure: Option<String>,
    pub audit: Audit,
}

fn pick_method(a: &Audit, curve: &Curve, opts: &SolveOptions) -> Method {
    if opts.method != Method::Auto {
        return opts.method;
    }
    if a.anomalous && matches!(curve.group, Group::Prime(_)) {
        return Method::Smart;
    }
    if !a.order_is_prime {
        return Method::PohligHellman;
    }
    if let Some((_, width)) = &opts.interval {
        return if width.bits() <= opts.bsgs_max_bits {
            Method::Bsgs
        } else {
            Method::Kangaroo
        };
    }
    Method::Rho
}

fn run_generic<G: EcdlpGroup>(
    g: &G,
    target: &G::Elt,
    method: Method,
    opts: &SolveOptions,
    audit: &Audit,
    report: &mut SolveReport,
) {
    match method {
        Method::Bsgs | Method::Kangaroo => {
            let (lo, width) = match &opts.interval {
                Some(iv) => iv.clone(),
                None if g.order().bits() <= 40 => (BigUint::zero(), g.order().clone()),
                None => {
                    report.failure = Some(format!(
                        "{method} needs an interval (lo, width); none given and the group is too large to use [0, n)"
                    ));
                    return;
                }
            };
            if method == Method::Bsgs {
                let r = bsgs(g, target, &lo, &width);
                report.scalar_value = r.scalar;
                report.group_ops = r.group_ops;
                report.expected_ops = r.expected_ops;
                report.table_size = r.table_size;
                report.threads = 1;
            } else {
                let r = kangaroo(
                    g,
                    target,
                    &lo,
                    &width,
                    &KangarooOptions {
                        threads: opts.threads,
                        seed: opts.seed,
                        dp_bits: opts.dp_bits,
                        max_ops: opts.max_iterations,
                        ..KangarooOptions::default()
                    },
                );
                report.scalar_value = r.scalar;
                report.group_ops = r.group_ops;
                report.expected_ops = r.expected_ops;
                report.table_size = r.table_size;
                report.threads = r.threads;
            }
            report.notes.push(format!(
                "interval [{lo}, {lo} + {width}) — a bounded, known-answer instance"
            ));
            if report.scalar_value.is_none() {
                report.failure = Some("no scalar in the interval (or budget exhausted)".into());
            }
        }
        Method::Rho => {
            let class = if opts.fold_automorphisms {
                audit.rho_class_size
            } else {
                1
            };
            let expected = rho::expected_rho_iterations(g.order(), class);
            report.expected_ops = expected;
            if expected.log2() > opts.feasibility_log2 && !opts.force {
                report.failure = Some(format!(
                    "whole-group rho on {} expects 2^{:.1} additions (class size {class}); above the 2^{:.0} feasibility bound — pass an interval or force a budgeted run",
                    g.name(),
                    expected.log2(),
                    opts.feasibility_log2
                ));
                return;
            }
            let r = pollard_rho(
                g,
                target,
                &RhoOptions {
                    threads: opts.threads,
                    seed: opts.seed,
                    dp_bits: opts.dp_bits,
                    max_iterations: opts.max_iterations,
                    fold_automorphisms: opts.fold_automorphisms,
                    ..RhoOptions::default()
                },
            );
            report.scalar_value = r.scalar;
            report.group_ops = r.iterations;
            report.table_size = r.distinguished_points as usize;
            report.threads = r.threads;
            report.notes.push(format!(
                "rho: class size {}, dp_bits {}, {} walks, {} abandoned, {} cycle escapes",
                r.class_size, r.dp_bits, r.walks, r.abandoned_walks, r.cycle_escapes
            ));
            if report.scalar_value.is_none() {
                report.failure = Some(format!(
                    "rho budget of {} additions exhausted",
                    r.iterations
                ));
            }
        }
        Method::PohligHellman => {
            let r = pohlig_hellman(
                g,
                target,
                &RhoOptions {
                    threads: opts.threads,
                    seed: opts.seed,
                    ..RhoOptions::default()
                },
            );
            report.scalar_value = r.scalar;
            report.group_ops = r.group_ops;
            report.expected_ops = r
                .factors
                .iter()
                .map(|(q, e)| 2.0 * rho::biguint_to_f64(q).sqrt() * *e as f64)
                .sum();
            report.notes.push(format!(
                "order factors: {}",
                r.factors
                    .iter()
                    .map(|(q, e)| format!("{q}^{e}"))
                    .collect::<Vec<_>>()
                    .join(" · ")
            ));
            report.failure = r.failure;
        }
        Method::Smart | Method::Auto => unreachable!("dispatched before run_generic"),
    }
}

/// Solve `Q = [k]G` on `curve` by the branch the audit (or `opts.method`)
/// selects.  Every returned scalar is verified by `[k]G = Q`.
pub fn solve(curve: &Curve, target: &AnyPoint, opts: &SolveOptions) -> SolveReport {
    let t0 = std::time::Instant::now();
    let a = audit(curve);
    let method = pick_method(&a, curve, opts);
    let mut report = SolveReport {
        curve: curve.name.clone(),
        family: curve.family,
        method,
        scalar: None,
        scalar_hex: None,
        scalar_value: None,
        verified: false,
        group_ops: 0,
        elapsed_ms: 0,
        expected_ops: 0.0,
        table_size: 0,
        threads: opts.threads.max(1),
        notes: Vec::new(),
        failure: None,
        audit: a.clone(),
    };
    if let Err(e) = curve.validate_target(target) {
        report.failure = Some(e);
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }
    if curve.is_identity(target) {
        report.scalar_value = Some(BigUint::zero());
    } else {
        match (&curve.group, target, method) {
            (Group::Prime(g), AnyPoint::Prime(t), Method::Smart) => {
                let r = anomalous::smart_attack_group(g, t);
                report.scalar_value = r.scalar;
                report.failure = r.failure;
                report.notes.push(format!(
                    "Smart attack: {} lift(s) of (a, b) tried",
                    r.lifts_tried
                ));
                report.expected_ops = 0.0;
            }
            (Group::Binary(_), _, Method::Smart) => {
                report.failure = Some("Smart's attack needs an anomalous prime-field curve".into());
            }
            (Group::Prime(g), AnyPoint::Prime(t), m) => {
                run_generic(g.as_ref(), t, m, opts, &a, &mut report)
            }
            (Group::Binary(g), AnyPoint::Binary(t), m) => {
                run_generic(g.as_ref(), t, m, opts, &a, &mut report)
            }
            _ => report.failure = Some("point does not belong to this curve's field".into()),
        }
    }
    if let Some(k) = &report.scalar_value {
        let check = curve
            .mul(&curve.generator(), k)
            .map(|p| p == *target)
            .unwrap_or(false);
        report.verified = check;
        if check {
            report.scalar = Some(k.to_string());
            report.scalar_hex = Some(format!("{k:x}"));
        } else {
            report.failure = Some("solver returned a scalar that fails [k]G = Q".into());
            report.scalar_value = None;
        }
    }
    report.elapsed_ms = t0.elapsed().as_millis();
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_fifteen_nist_curves_load_and_verify() {
        let curves = nist_curves();
        assert_eq!(curves.len(), 15);
        let mut fams = std::collections::HashMap::new();
        for c in &curves {
            *fams.entry(c.family).or_insert(0) += 1;
            let g = c.generator();
            assert!(c.validate_target(&g).is_ok(), "{}: generator", c.name);
            assert!(
                c.is_identity(&c.mul(&g, c.order()).unwrap()),
                "{}: [n]G = O",
                c.name
            );
        }
        assert_eq!(fams[&Family::Prime], 5);
        assert_eq!(fams[&Family::Koblitz], 5);
        assert_eq!(fams[&Family::BinaryRandom], 5);
    }

    #[test]
    fn audit_finds_no_structural_weakness_on_nist_curves() {
        for c in nist_curves() {
            let a = audit(&c);
            assert!(!a.anomalous, "{}: anomalous", c.name);
            assert!(a.order_is_prime, "{}: composite order", c.name);
            assert!(
                a.embedding_degree.is_none(),
                "{}: embedding degree {:?}",
                c.name,
                a.embedding_degree
            );
            assert_eq!(a.recommended, Method::Rho);
            match c.family {
                Family::Prime => {
                    assert_eq!(a.cofactor, "1");
                    assert_eq!(a.rho_class_size, 2);
                }
                Family::Koblitz => {
                    assert_eq!(a.extension_degree_prime, Some(true));
                    assert_eq!(a.koblitz_order_consistent, Some(true), "{}", c.name);
                    assert_eq!(a.rho_class_size, 2 * c.field_bits as u32);
                    // FIPS 186-4: K-163 has a = 1 (h = 2); the other four have a = 0 (h = 4).
                    let expected_h = if c.name == "K-163" { "2" } else { "4" };
                    assert_eq!(a.cofactor, expected_h, "{}", c.name);
                    assert!(a.rho_log2_expected < a.rho_log2_expected_unfolded - 3.0);
                }
                Family::BinaryRandom => {
                    assert_eq!(a.extension_degree_prime, Some(true));
                    assert_eq!(a.cofactor, "2");
                    assert_eq!(a.rho_class_size, 2);
                }
            }
            // Trace is small relative to q (Hasse) and never 1.
            assert_ne!(a.trace, "1");
        }
    }

    #[test]
    fn koblitz_lucas_order_matches_known_small_values() {
        // From the Lucas recurrence by hand: #E_1(F_2^17) = 131174, #E_0(F_2^19) = 523492.
        assert_eq!(koblitz_group_order(1, 17), BigUint::from(131_174u32));
        assert_eq!(koblitz_group_order(0, 19), BigUint::from(523_492u32));
        // And at full size, against FIPS 186-4: K-163 has h = 2, n given in the catalog.
        let k163 = curve_by_name("K-163").unwrap();
        assert_eq!(koblitz_group_order(1, 163), k163.group_order);
    }

    #[test]
    fn bounded_instances_solve_on_every_nist_curve() {
        let lo = BigUint::from(1u32) << 40u32;
        let width = BigUint::from(1u32) << 16u32;
        for (i, c) in nist_curves().iter().enumerate() {
            let k = &lo + BigUint::from(1000u32 + 37 * i as u32);
            let q = c.plant(&k);
            for method in [Method::Bsgs, Method::Kangaroo] {
                let rep = solve(
                    c,
                    &q,
                    &SolveOptions {
                        method,
                        interval: Some((lo.clone(), width.clone())),
                        threads: 2,
                        ..SolveOptions::default()
                    },
                );
                assert_eq!(
                    rep.scalar_value.as_ref(),
                    Some(&k),
                    "{} {method}: {:?}",
                    c.name,
                    rep.failure
                );
                assert!(rep.verified);
            }
            // Auto picks BSGS for a 16-bit interval.
            let rep = solve(
                c,
                &q,
                &SolveOptions {
                    interval: Some((lo.clone(), width.clone())),
                    ..SolveOptions::default()
                },
            );
            assert_eq!(rep.method, Method::Bsgs);
            assert_eq!(rep.scalar_value, Some(k));
        }
    }

    #[test]
    fn whole_group_rho_on_a_nist_curve_is_refused_without_force() {
        let c = curve_by_name("P-256").unwrap();
        let q = c.plant(&BigUint::from(12345u32));
        let rep = solve(&c, &q, &SolveOptions::default());
        assert_eq!(rep.method, Method::Rho);
        assert!(rep.scalar.is_none());
        assert!(rep.failure.unwrap().contains("feasibility"));
        // Forced with a tiny budget: runs, stops, reports honestly.
        let rep = solve(
            &c,
            &q,
            &SolveOptions {
                force: true,
                max_iterations: Some(2000),
                dp_bits: Some(4),
                ..SolveOptions::default()
            },
        );
        assert!(rep.scalar.is_none());
        assert!(
            rep.group_ops >= 2000 && rep.group_ops < 6000,
            "{}",
            rep.group_ops
        );
    }

    #[test]
    fn dispatcher_routes_anomalous_toy_and_full_group_curves() {
        let an = curve_by_name("anomalous-96-2").unwrap();
        let a = audit(&an);
        assert!(a.anomalous);
        assert_eq!(a.recommended, Method::Smart);
        let k = an.order() / 3u32;
        let rep = solve(&an, &an.plant(&k), &SolveOptions::default());
        assert_eq!(rep.method, Method::Smart);
        assert_eq!(rep.scalar_value, Some(k));

        let toy_k = curve_by_name("toy-k23a1").unwrap();
        let k = BigUint::from(4_000_000u32);
        let rep = solve(&toy_k, &toy_k.plant(&k), &SolveOptions::default());
        assert_eq!(rep.method, Method::Rho);
        assert_eq!(rep.scalar_value, Some(k));

        let full = curve_by_name("toy-k19a0-full").unwrap();
        let k = BigUint::from(123_456u32);
        let rep = solve(&full, &full.plant(&k), &SolveOptions::default());
        assert_eq!(rep.method, Method::PohligHellman);
        assert_eq!(rep.scalar_value, Some(k));
    }

    #[test]
    fn names_resolve_in_every_spelling() {
        for s in [
            "P-256",
            "p256",
            "secp256r1",
            "prime256v1",
            "nistp256",
            "P_256",
        ] {
            assert_eq!(curve_by_name(s).unwrap().name, "P-256", "{s}");
        }
        for s in ["K-163", "k163", "sect163k1"] {
            assert_eq!(curve_by_name(s).unwrap().name, "K-163", "{s}");
        }
        assert_eq!(
            curve_by_name("B-163").unwrap().catalog_name,
            Some("sect163r2")
        );
        assert!(curve_by_name("secp256k1").is_err(), "not a NIST curve");
        assert!(curve_by_name("nope").is_err());
    }

    #[test]
    fn external_points_are_validated() {
        let c = curve_by_name("K-233").unwrap();
        let (gx, gy) = c.point_hex(&c.generator()).unwrap();
        assert!(c.parse_point(&gx, &gy).is_ok());
        assert!(c.parse_point(&gx, "1").is_err());
    }
}
