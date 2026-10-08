//! Streaming certificate for a million-curve P-256 isogeny grid.
//!
//! The ordinary breadth-first walker deliberately records the whole explored
//! multigraph and every root route.  Those are useful records at modest size,
//! but route duplication is quadratic on long paths.  This module records the
//! smaller object needed for a large enumeration: one kernel-certified parent
//! edge per new canonical curve in a fixed two-generator product grid.

mod census;
mod window;

pub use census::{
    generate_trait_census_jsonl, verify_trait_census_paths, TraitCensusReceipt, TraitCensusSummary,
    TRAIT_CENSUS_RECEIPT_SCHEMA, TRAIT_CENSUS_SCHEMA, TRAIT_CENSUS_SUMMARY_SCHEMA,
};
pub use window::{
    audit_j_union_paths, generate_strip_jsonl, verify_strip_jsonl, JUnionReceipt, StripConfig,
    StripHeaderRecord, UnionInputReceipt, J_UNION_SCHEMA, STRIP_RECEIPT_SCHEMA, STRIP_SCHEMA,
    STRIP_SUMMARY_SCHEMA,
};

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};
use std::time::Instant;

use num_bigint::{BigInt, BigUint};
use rayon::prelude::*;
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use super::curve::{self, Model, OrderAudit};
use super::field::{is_probable_prime, Fe, Field};
use super::kernel;
use super::modpoly::{self, ModPoly};
use super::poly::{self, Poly};
use super::record;
use super::walk::{ClassInfo, EllKind, StartCurve};
use crate::cryptanalysis::curve_id;
use crate::hash::sha256::sha256;

pub const SCHEMA: &str = "p256.isogeny-grid-certificate/v1";
pub const SUMMARY_SCHEMA: &str = "p256.isogeny-grid-summary/v1";
pub const RECEIPT_SCHEMA: &str = "p256.isogeny-grid-receipt/v1";
pub const SPINE_DEGREE: u64 = 13;
pub const ROW_DEGREE: u64 = 11;
pub const ROOT_SEED: u64 = 1;
pub const PRODUCTION_SIDE: u32 = 1_000;
pub const PRIMARY_QR_THRESHOLD: u32 = 54;
pub const COMPARABILITY_QR_THRESHOLD: u32 = 52;
pub const COEFFICIENT_BITS_THRESHOLD: u64 = 224;

#[derive(Clone, Debug)]
pub struct GridConfig {
    pub side: u32,
    pub batch_rows: usize,
    pub audit_points: usize,
    pub audit_seed_x: u64,
    pub source_commit: String,
}

impl GridConfig {
    pub fn production(source_commit: String) -> Self {
        Self {
            side: PRODUCTION_SIDE,
            batch_rows: (rayon::current_num_threads() * 2).max(1),
            audit_points: 1,
            audit_seed_x: 0,
            source_commit,
        }
    }

    fn validate(&self) -> Result<(), String> {
        if self.side < 2 {
            return Err("grid side must be at least two".into());
        }
        self.side
            .checked_mul(self.side)
            .ok_or("grid side squared does not fit u32")?;
        if self.batch_rows == 0 {
            return Err("row batch must be positive".into());
        }
        if self.audit_points == 0 {
            return Err("audit point count must be positive".into());
        }
        if self.audit_points != 1 || self.audit_seed_x != 0 {
            return Err("frozen construction uses one order-audit point from x = 0".into());
        }
        if self.source_commit.len() != 40
            || !self.source_commit.bytes().all(|b| b.is_ascii_hexdigit())
        {
            return Err("source commit must be a 40-digit hexadecimal Git object id".into());
        }
        Ok(())
    }

    fn target_curves(&self) -> u64 {
        u64::from(self.side) * u64::from(self.side)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Coordinate {
    pub x: u32,
    pub y: u32,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DetectorRecord {
    pub non_singular: bool,
    pub generator_valid: bool,
    pub a_minus_3_model: bool,
    pub qr_prefix_64: u32,
    pub a_signed_bits: u64,
    pub b_signed_bits: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CurveRecord {
    pub a: String,
    pub b: String,
    pub j: String,
    pub generator_x: String,
    pub generator_y: String,
    pub icv1_slug: String,
    pub icv1_identity: String,
    pub ec1_alias: String,
    pub curve_uid: String,
    pub order_status: String,
    pub detectors: DetectorRecord,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HeaderRecord {
    pub record: String,
    pub schema: String,
    pub source_commit: String,
    pub root_name: String,
    pub root_icv1_slug: String,
    pub field_p: String,
    pub group_order: String,
    pub trace: String,
    pub side: u32,
    pub target_curves: u64,
    pub spine_degree: u64,
    pub row_degree: u64,
    pub root_seed: u64,
    pub construction_audit_points: usize,
    pub construction_audit_seed_x: u64,
    pub canonical_model_rule: String,
    pub generator_rule: String,
    pub first_edge_rule: String,
    pub later_edge_rule: String,
    pub spine_phi_checks: usize,
    pub row_phi_checks: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RootRecord {
    pub record: String,
    pub index: u64,
    pub coordinate: Coordinate,
    pub curve: CurveRecord,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NodeRecord {
    pub record: String,
    pub index: u64,
    pub coordinate: Coordinate,
    pub parent: u64,
    pub degree: u64,
    pub predecessor_j: Option<String>,
    pub kernel_polynomial: Vec<String>,
    pub isomorphism_u2: String,
    pub curve: CurveRecord,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CurveLocator {
    pub index: u64,
    pub coordinate: Coordinate,
    pub icv1_slug: String,
    pub icv1_identity: String,
    pub ec1_alias: String,
    pub curve_uid: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Extremum {
    pub value: u64,
    pub curve: CurveLocator,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScreenHit {
    pub curve: CurveLocator,
    pub qr_prefix_64: u32,
    pub a_signed_bits: u64,
    pub b_signed_bits: u64,
    pub a_minus_3_model: bool,
    pub reasons: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SlugCollision {
    pub slug: String,
    pub first_index: u64,
    pub repeated_index: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HistogramBin {
    pub value: u32,
    pub count: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GridSummary {
    pub record: String,
    pub schema: String,
    pub bound_records: u64,
    pub record_chain_sha256: String,
    pub unique_curves: u64,
    pub parent_edges: u64,
    pub kernel_certificates: u64,
    pub order_proved_prime: u64,
    pub non_singular: u64,
    pub generator_valid: u64,
    pub unique_j_invariants: u64,
    pub unique_canonical_models: u64,
    pub unique_icv1_identities: u64,
    pub unique_ec1_aliases: u64,
    pub unique_curve_uids: u64,
    pub unique_icv1_slugs: u64,
    pub icv1_slug_collisions: Vec<SlugCollision>,
    pub a_minus_3_models: u64,
    pub qr_prefix_64_min: u32,
    pub qr_prefix_64_max: u32,
    pub qr_prefix_64_histogram: Vec<HistogramBin>,
    pub primary_screen_hits: Vec<ScreenHit>,
    pub qr_52_comparability_hits: Vec<ScreenHit>,
    pub minimum_b_signed_bits: Extremum,
    pub minimum_a_signed_bits_without_a_minus_3: Option<Extremum>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenerationReceipt {
    pub schema: String,
    pub status: String,
    pub source_commit: String,
    pub elapsed_ms: u128,
    pub resources: ResourceUsage,
    pub uncompressed_bytes: u64,
    pub summary: GridSummary,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VerificationReceipt {
    pub schema: String,
    pub status: String,
    pub source_commit: String,
    pub audit_points: usize,
    pub audit_seed_x: u64,
    pub elapsed_ms: u128,
    pub resources: ResourceUsage,
    pub uncompressed_bytes: u64,
    pub summary: GridSummary,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResourceUsage {
    pub status: String,
    pub user_cpu_micros: Option<u64>,
    pub system_cpu_micros: Option<u64>,
    pub peak_rss_bytes: Option<u64>,
}

#[cfg(unix)]
fn resource_usage() -> ResourceUsage {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::zeroed();
    // SAFETY: `getrusage` initializes the supplied `rusage` on success, and
    // the pointer is valid for the duration of the call.
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
        return ResourceUsage {
            status: "unavailable_getrusage_failed".into(),
            user_cpu_micros: None,
            system_cpu_micros: None,
            peak_rss_bytes: None,
        };
    }
    // SAFETY: the successful call above initialized every field.
    let usage = unsafe { usage.assume_init() };
    let micros = |time: libc::timeval| {
        (time.tv_sec.max(0) as u64)
            .saturating_mul(1_000_000)
            .saturating_add(time.tv_usec.max(0) as u64)
    };
    #[cfg(target_os = "macos")]
    let rss_multiplier = 1u64;
    #[cfg(not(target_os = "macos"))]
    let rss_multiplier = 1_024u64;
    ResourceUsage {
        status: "measured_getrusage_self".into(),
        user_cpu_micros: Some(micros(usage.ru_utime)),
        system_cpu_micros: Some(micros(usage.ru_stime)),
        peak_rss_bytes: Some((usage.ru_maxrss.max(0) as u64).saturating_mul(rss_multiplier)),
    }
}

#[cfg(not(unix))]
fn resource_usage() -> ResourceUsage {
    ResourceUsage {
        status: "unavailable_non_unix".into(),
        user_cpu_micros: None,
        system_cpu_micros: None,
        peak_rss_bytes: None,
    }
}

#[derive(Clone, Copy)]
struct GridPoint {
    model: Model,
    j: Fe,
}

struct GeneratedNode {
    record: NodeRecord,
    model: Model,
    j: Fe,
}

fn hex_fe(f: &Field, value: &Fe) -> String {
    format!("{:064x}", f.to_big(value))
}

fn parse_fe(f: &Field, value: &str, what: &str) -> Result<Fe, String> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(format!("{what}: expected 64 lowercase hexadecimal digits"));
    }
    let integer = BigUint::parse_bytes(value.as_bytes(), 16)
        .ok_or_else(|| format!("{what}: invalid hexadecimal field element"))?;
    if &integer >= f.modulus() {
        return Err(format!("{what}: field element is not canonical"));
    }
    Ok(f.from_big(&integer))
}

fn coefficient_bits(f: &Field, value: &Fe) -> u64 {
    let positive = f.to_big(value);
    let negative = f.modulus() - &positive;
    positive.min(negative).bits()
}

fn describe_curve(
    f: &Field,
    start: &StartCurve,
    model: &Model,
    is_root: bool,
    order_is_prime: bool,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<CurveRecord, String> {
    let j = model.j(f).ok_or("singular curve")?;
    let (a, b, j_big) = (f.to_big(&model.a), f.to_big(&model.b), f.to_big(&j));
    let generator = if is_root {
        (f.from_big(&start.gx), f.from_big(&start.gy))
    } else {
        curve::deterministic_generator(f, model, &start.cofactor)
    };
    let audit = curve::audit_order(
        f,
        model,
        &start.order,
        order_is_prime,
        audit_points,
        audit_seed_x,
    );
    if audit != OrderAudit::ProvedPrime {
        return Err(format!("order audit returned {audit:?}"));
    }
    let non_singular = !f.is_zero(&model.discriminant(f));
    let generator_valid = model.on_curve(f, &generator.0, &generator.1)
        && curve::scalar_mul(f, model, &generator.0, &generator.1, &start.subgroup_order).is_none();
    if !non_singular || !generator_valid {
        return Err("curve detector rejected the model or generator".into());
    }
    let identity = curve_id::prime_with_verified_j(&start.p, &a, &b, &start.order, &j_big)
        .ok_or("failed to construct the ICV1 identity")?;
    let (gx, gy) = (f.to_big(&generator.0), f.to_big(&generator.1));
    let (field_record, curve_record) = record::ec1_records(
        &start.p,
        &a,
        &b,
        &start.subgroup_order,
        &start.cofactor,
        (&gx, &gy),
    );
    let tag = if is_root { start.tag.as_str() } else { "fp" };
    let (ec1_alias, curve_uid) = record::ec1_identity(&start.p, &field_record, &curve_record, tag);
    Ok(CurveRecord {
        a: hex_fe(f, &model.a),
        b: hex_fe(f, &model.b),
        j: hex_fe(f, &j),
        generator_x: hex_fe(f, &generator.0),
        generator_y: hex_fe(f, &generator.1),
        icv1_slug: identity.slug,
        icv1_identity: identity.icv1,
        ec1_alias,
        curve_uid,
        order_status: "proved_prime".into(),
        detectors: DetectorRecord {
            non_singular,
            generator_valid,
            a_minus_3_model: !curve::a_minus_3_models(f, model).is_empty(),
            qr_prefix_64: curve::qr_prefix_64(f, model),
            a_signed_bits: coefficient_bits(f, &model.a),
            b_signed_bits: coefficient_bits(f, &model.b),
        },
    })
}

fn rational_roots(f: &Field, phi: &ModPoly, model: &Model) -> Result<[Fe; 2], String> {
    let j = model.j(f).ok_or("singular source curve")?;
    let roots = poly::roots(f, &phi.at_x(f, &j), ROOT_SEED ^ phi.ell);
    match roots.as_slice() {
        [left, right] if left != right => Ok([*left, *right]),
        _ => Err(format!(
            "Phi_{} has {} distinct rational roots, expected two",
            phi.ell,
            roots.len()
        )),
    }
}

fn selected_target(
    f: &Field,
    phi: &ModPoly,
    source: &Model,
    predecessor: Option<Fe>,
) -> Result<Fe, String> {
    let roots = rational_roots(f, phi, source)?;
    match predecessor {
        None => Ok(roots[0]),
        Some(previous) if roots[0] == previous => Ok(roots[1]),
        Some(previous) if roots[1] == previous => Ok(roots[0]),
        Some(_) => Err(format!(
            "the predecessor is not a degree-{} neighbour",
            phi.ell
        )),
    }
}

fn make_step(
    f: &Field,
    phi: &ModPoly,
    source: GridPoint,
    predecessor: Option<Fe>,
) -> Result<(Model, Fe, Poly, Fe), String> {
    let target_j = selected_target(f, phi, &source.model, predecessor)?;
    let isogeny = kernel::elkies_isogeny(f, &source.model, phi, &target_j)
        .map_err(|error| format!("degree-{} edge construction: {error:?}", phi.ell))?;
    let (target, canonical_u2) = curve::canonical_model(f, &isogeny.codomain)
        .ok_or("isogeny reached an exceptional target j-invariant")?;
    let linked_u2 = kernel::link_to_target(f, &isogeny.codomain, &target)
        .map_err(|error| format!("degree-{} target link: {error:?}", phi.ell))?;
    if canonical_u2 != linked_u2 || target.j(f) != Some(target_j) {
        return Err("canonical target or isomorphism changed the selected j-invariant".into());
    }
    Ok((target, target_j, isogeny.kernel, linked_u2))
}

fn node_index(side: u32, x: u32, y: u32) -> u64 {
    if x == 0 {
        u64::from(y)
    } else {
        u64::from(side) + u64::from(y) * u64::from(side - 1) + u64::from(x - 1)
    }
}

fn node_parent(side: u32, x: u32, y: u32) -> u64 {
    if x == 0 {
        u64::from(y - 1)
    } else if x == 1 {
        u64::from(y)
    } else {
        node_index(side, x - 1, y)
    }
}

fn make_node_record(
    f: &Field,
    start: &StartCurve,
    side: u32,
    x: u32,
    y: u32,
    phi: &ModPoly,
    source: GridPoint,
    predecessor: Option<Fe>,
    order_is_prime: bool,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<GeneratedNode, String> {
    let (model, j, kernel, u2) = make_step(f, phi, source, predecessor)?;
    let curve = describe_curve(
        f,
        start,
        &model,
        false,
        order_is_prime,
        audit_points,
        audit_seed_x,
    )?;
    Ok(GeneratedNode {
        record: NodeRecord {
            record: "node".into(),
            index: node_index(side, x, y),
            coordinate: Coordinate { x, y },
            parent: node_parent(side, x, y),
            degree: phi.ell,
            predecessor_j: predecessor.map(|value| hex_fe(f, &value)),
            kernel_polynomial: kernel.iter().map(|value| hex_fe(f, value)).collect(),
            isomorphism_u2: hex_fe(f, &u2),
            curve,
        },
        model,
        j,
    })
}

fn locator(index: u64, coordinate: Coordinate, curve: &CurveRecord) -> CurveLocator {
    CurveLocator {
        index,
        coordinate,
        icv1_slug: curve.icv1_slug.clone(),
        icv1_identity: curve.icv1_identity.clone(),
        ec1_alias: curve.ec1_alias.clone(),
        curve_uid: curve.curve_uid.clone(),
    }
}

fn screen_hit(index: u64, coordinate: Coordinate, curve: &CurveRecord) -> ScreenHit {
    let detector = &curve.detectors;
    let mut reasons = Vec::new();
    if detector.qr_prefix_64 >= PRIMARY_QR_THRESHOLD {
        reasons.push(format!("qr_prefix_64>={PRIMARY_QR_THRESHOLD}"));
    }
    if detector.b_signed_bits <= COEFFICIENT_BITS_THRESHOLD {
        reasons.push(format!("b_signed_bits<={COEFFICIENT_BITS_THRESHOLD}"));
    }
    if !detector.a_minus_3_model && detector.a_signed_bits <= COEFFICIENT_BITS_THRESHOLD {
        reasons.push(format!(
            "non_a_minus_3_a_signed_bits<={COEFFICIENT_BITS_THRESHOLD}"
        ));
    }
    ScreenHit {
        curve: locator(index, coordinate, curve),
        qr_prefix_64: detector.qr_prefix_64,
        a_signed_bits: detector.a_signed_bits,
        b_signed_bits: detector.b_signed_bits,
        a_minus_3_model: detector.a_minus_3_model,
        reasons,
    }
}

#[derive(Default)]
struct Accumulator {
    seen: u64,
    j_invariants: HashSet<Fe>,
    models: HashSet<(Fe, Fe)>,
    icv1_identities: HashSet<String>,
    ec1_aliases: HashSet<String>,
    curve_uids: HashSet<String>,
    slug_first: HashMap<String, u64>,
    slug_collisions: Vec<SlugCollision>,
    qr_histogram: Vec<u64>,
    qr_min: u32,
    qr_max: u32,
    a_minus_3: u64,
    primary_hits: Vec<ScreenHit>,
    comparability_hits: Vec<ScreenHit>,
    minimum_b: Option<Extremum>,
    minimum_non_a3_a: Option<Extremum>,
}

impl Accumulator {
    fn new() -> Self {
        Self {
            qr_min: 64,
            qr_histogram: vec![0; 65],
            ..Self::default()
        }
    }

    fn observe(
        &mut self,
        index: u64,
        coordinate: Coordinate,
        model: &Model,
        j: Fe,
        curve: &CurveRecord,
    ) -> Result<(), String> {
        if !self.j_invariants.insert(j) {
            return Err(format!(
                "grid collision at ({},{}): duplicate j {}",
                coordinate.x, coordinate.y, curve.j
            ));
        }
        if !self.models.insert((model.a, model.b)) {
            return Err(format!(
                "grid collision at ({},{}): duplicate canonical model",
                coordinate.x, coordinate.y
            ));
        }
        if !self.icv1_identities.insert(curve.icv1_identity.clone()) {
            return Err(format!("node {index}: duplicate ICV1 identity"));
        }
        if !self.ec1_aliases.insert(curve.ec1_alias.clone()) {
            return Err(format!(
                "node {index}: duplicate 48-bit EC1 alias {}",
                curve.ec1_alias
            ));
        }
        if !self.curve_uids.insert(curve.curve_uid.clone()) {
            return Err(format!("node {index}: duplicate full curve UID"));
        }
        if let Some(first_index) = self.slug_first.get(&curve.icv1_slug) {
            self.slug_collisions.push(SlugCollision {
                slug: curve.icv1_slug.clone(),
                first_index: *first_index,
                repeated_index: index,
            });
        } else {
            self.slug_first.insert(curve.icv1_slug.clone(), index);
        }
        if curve.order_status != "proved_prime"
            || !curve.detectors.non_singular
            || !curve.detectors.generator_valid
        {
            return Err(format!("node {index}: failed order or detector gate"));
        }
        let detector = &curve.detectors;
        self.qr_histogram[detector.qr_prefix_64 as usize] += 1;
        self.qr_min = self.qr_min.min(detector.qr_prefix_64);
        self.qr_max = self.qr_max.max(detector.qr_prefix_64);
        self.a_minus_3 += u64::from(detector.a_minus_3_model);
        let here = locator(index, coordinate, curve);
        if self
            .minimum_b
            .as_ref()
            .is_none_or(|old| detector.b_signed_bits < old.value)
        {
            self.minimum_b = Some(Extremum {
                value: detector.b_signed_bits,
                curve: here.clone(),
            });
        }
        if !detector.a_minus_3_model
            && self
                .minimum_non_a3_a
                .as_ref()
                .is_none_or(|old| detector.a_signed_bits < old.value)
        {
            self.minimum_non_a3_a = Some(Extremum {
                value: detector.a_signed_bits,
                curve: here,
            });
        }
        let hit = screen_hit(index, coordinate, curve);
        if !hit.reasons.is_empty() {
            self.primary_hits.push(hit.clone());
        }
        if detector.qr_prefix_64 >= COMPARABILITY_QR_THRESHOLD {
            self.comparability_hits.push(hit);
        }
        self.seen += 1;
        Ok(())
    }

    fn summary(&self, chain: [u8; 32]) -> Result<GridSummary, String> {
        let minimum_b = self.minimum_b.clone().ok_or("no b-bit minimum")?;
        let histogram = self
            .qr_histogram
            .iter()
            .enumerate()
            .filter(|(_, count)| **count != 0)
            .map(|(value, count)| HistogramBin {
                value: value as u32,
                count: *count,
            })
            .collect();
        Ok(GridSummary {
            record: "summary".into(),
            schema: SUMMARY_SCHEMA.into(),
            bound_records: self.seen + 1,
            record_chain_sha256: hex::encode(chain),
            unique_curves: self.seen,
            parent_edges: self.seen.saturating_sub(1),
            kernel_certificates: self.seen.saturating_sub(1),
            order_proved_prime: self.seen,
            non_singular: self.seen,
            generator_valid: self.seen,
            unique_j_invariants: self.j_invariants.len() as u64,
            unique_canonical_models: self.models.len() as u64,
            unique_icv1_identities: self.icv1_identities.len() as u64,
            unique_ec1_aliases: self.ec1_aliases.len() as u64,
            unique_curve_uids: self.curve_uids.len() as u64,
            unique_icv1_slugs: self.slug_first.len() as u64,
            icv1_slug_collisions: self.slug_collisions.clone(),
            a_minus_3_models: self.a_minus_3,
            qr_prefix_64_min: self.qr_min,
            qr_prefix_64_max: self.qr_max,
            qr_prefix_64_histogram: histogram,
            primary_screen_hits: self.primary_hits.clone(),
            qr_52_comparability_hits: self.comparability_hits.clone(),
            minimum_b_signed_bits: minimum_b,
            minimum_a_signed_bits_without_a_minus_3: self.minimum_non_a3_a.clone(),
        })
    }
}

fn update_chain(chain: &mut [u8; 32], line: &[u8]) {
    let mut input = Vec::with_capacity(chain.len() + line.len());
    input.extend_from_slice(chain);
    input.extend_from_slice(line);
    *chain = sha256(&input);
}

struct LineWriter<'a, W: Write> {
    writer: &'a mut W,
    chain: [u8; 32],
    bytes: u64,
}

impl<'a, W: Write> LineWriter<'a, W> {
    fn new(writer: &'a mut W) -> Self {
        Self {
            writer,
            chain: [0; 32],
            bytes: 0,
        }
    }

    fn write<T: Serialize>(&mut self, value: &T, bind: bool) -> Result<(), String> {
        let mut line = serde_json::to_vec(value).map_err(|error| error.to_string())?;
        line.push(b'\n');
        if bind {
            update_chain(&mut self.chain, &line);
        }
        self.writer
            .write_all(&line)
            .map_err(|error| format!("write certificate: {error}"))?;
        self.bytes += line.len() as u64;
        Ok(())
    }
}

fn build_context() -> Result<(StartCurve, Field, ModPoly, ModPoly), String> {
    let start = StartCurve::p256();
    let field = Field::new(&start.p).ok_or("P-256 field failed primality validation")?;
    let class = ClassInfo::compute(&start, &[ROW_DEGREE, SPINE_DEGREE]);
    for degree in [ROW_DEGREE, SPINE_DEGREE] {
        let info = class
            .ells
            .iter()
            .find(|entry| entry.0 == degree)
            .ok_or_else(|| format!("missing class information for degree {degree}"))?;
        if info.1 != EllKind::Elkies
            || class.depth(degree) != 0
            || class.expected_neighbours(degree) != Some(2)
        {
            return Err(format!(
                "degree {degree} is not a split, two-neighbour, depth-zero action"
            ));
        }
    }
    let built: Vec<Result<ModPoly, String>> = [SPINE_DEGREE, ROW_DEGREE]
        .par_iter()
        .map(|degree| modpoly::modular_polynomial(&field, *degree))
        .collect();
    let mut built = built.into_iter();
    let spine = built.next().expect("two modular polynomials")?;
    let row = built.next().expect("two modular polynomials")?;
    Ok((start, field, spine, row))
}

fn make_header(
    config: &GridConfig,
    start: &StartCurve,
    root: &CurveRecord,
    spine: &ModPoly,
    row: &ModPoly,
) -> HeaderRecord {
    let trace: BigInt =
        BigInt::from(start.p.clone()) + BigInt::from(1u8) - BigInt::from(start.order.clone());
    HeaderRecord {
        record: "header".into(),
        schema: SCHEMA.into(),
        source_commit: config.source_commit.clone(),
        root_name: start.name.clone(),
        root_icv1_slug: root.icv1_slug.clone(),
        field_p: format!("{:064x}", start.p),
        group_order: format!("{:064x}", start.order),
        trace: trace.to_string(),
        side: config.side,
        target_curves: config.target_curves(),
        spine_degree: SPINE_DEGREE,
        row_degree: ROW_DEGREE,
        root_seed: ROOT_SEED,
        construction_audit_points: config.audit_points,
        construction_audit_seed_x: config.audit_seed_x,
        canonical_model_rule: curve::CANONICAL_RULE.into(),
        generator_rule: curve::GENERATOR_RULE.into(),
        first_edge_rule: "choose the numerically smaller rational root".into(),
        later_edge_rule: "choose the rational root unequal to the predecessor j".into(),
        spine_phi_checks: spine.checks_passed,
        row_phi_checks: row.checks_passed,
    }
}

pub fn generate_jsonl<W: Write>(
    writer: &mut W,
    config: &GridConfig,
) -> Result<GenerationReceipt, String> {
    config.validate()?;
    let started = Instant::now();
    let (start, field, spine_phi, row_phi) = build_context()?;
    let order_is_prime = is_probable_prime(&start.order);
    if !order_is_prime {
        return Err("the registered P-256 group order is not prime".into());
    }
    let root_model = Model {
        a: field.from_big(&start.a),
        b: field.from_big(&start.b),
    };
    let root_j = root_model.j(&field).ok_or("registered P-256 is singular")?;
    let root_curve = describe_curve(
        &field,
        &start,
        &root_model,
        true,
        order_is_prime,
        config.audit_points,
        config.audit_seed_x,
    )?;
    let header = make_header(config, &start, &root_curve, &spine_phi, &row_phi);
    let mut output = LineWriter::new(writer);
    output.write(&header, true)?;
    let root_record = RootRecord {
        record: "root".into(),
        index: 0,
        coordinate: Coordinate { x: 0, y: 0 },
        curve: root_curve,
    };
    let mut accumulator = Accumulator::new();
    accumulator.observe(
        0,
        root_record.coordinate,
        &root_model,
        root_j,
        &root_record.curve,
    )?;
    output.write(&root_record, true)?;

    let mut spine = Vec::with_capacity(config.side as usize);
    spine.push(GridPoint {
        model: root_model,
        j: root_j,
    });
    for y in 1..config.side {
        let predecessor = if y == 1 {
            None
        } else {
            Some(spine[(y - 2) as usize].j)
        };
        let generated = make_node_record(
            &field,
            &start,
            config.side,
            0,
            y,
            &spine_phi,
            spine[(y - 1) as usize],
            predecessor,
            order_is_prime,
            config.audit_points,
            config.audit_seed_x,
        )?;
        accumulator.observe(
            generated.record.index,
            generated.record.coordinate,
            &generated.model,
            generated.j,
            &generated.record.curve,
        )?;
        output.write(&generated.record, true)?;
        spine.push(GridPoint {
            model: generated.model,
            j: generated.j,
        });
    }
    eprintln!(
        "p256-isogeny-million: generated {} spine curves",
        spine.len()
    );

    let batch_rows = config.batch_rows.min(config.side as usize).max(1);
    for batch_start in (0..config.side as usize).step_by(batch_rows) {
        let batch_end = (batch_start + batch_rows).min(config.side as usize);
        let rows: Vec<Result<Vec<GeneratedNode>, String>> = (batch_start..batch_end)
            .into_par_iter()
            .map(|y| {
                let mut source = spine[y];
                let mut predecessor = None;
                let mut row = Vec::with_capacity((config.side - 1) as usize);
                for x in 1..config.side {
                    let generated = make_node_record(
                        &field,
                        &start,
                        config.side,
                        x,
                        y as u32,
                        &row_phi,
                        source,
                        predecessor,
                        order_is_prime,
                        config.audit_points,
                        config.audit_seed_x,
                    )?;
                    predecessor = Some(source.j);
                    source = GridPoint {
                        model: generated.model,
                        j: generated.j,
                    };
                    row.push(generated);
                }
                Ok(row)
            })
            .collect();
        for row in rows {
            for generated in row? {
                accumulator.observe(
                    generated.record.index,
                    generated.record.coordinate,
                    &generated.model,
                    generated.j,
                    &generated.record.curve,
                )?;
                output.write(&generated.record, true)?;
            }
        }
        eprintln!(
            "p256-isogeny-million: rows {batch_end}/{}; {} unique curves",
            config.side, accumulator.seen
        );
    }
    if accumulator.seen != config.target_curves() {
        return Err(format!(
            "generated {} curves, expected {}",
            accumulator.seen,
            config.target_curves()
        ));
    }
    let summary = accumulator.summary(output.chain)?;
    output.write(&summary, false)?;
    output
        .writer
        .flush()
        .map_err(|error| format!("flush certificate: {error}"))?;
    Ok(GenerationReceipt {
        schema: RECEIPT_SCHEMA.into(),
        status: "pass".into(),
        source_commit: config.source_commit.clone(),
        elapsed_ms: started.elapsed().as_millis(),
        resources: resource_usage(),
        uncompressed_bytes: output.bytes,
        summary,
    })
}

struct LineReader<R: BufRead> {
    reader: R,
    chain: [u8; 32],
    bytes: u64,
    line_number: u64,
}

impl<R: BufRead> LineReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            chain: [0; 32],
            bytes: 0,
            line_number: 0,
        }
    }

    fn read<T: DeserializeOwned>(&mut self, bind: bool) -> Result<T, String> {
        let mut line = Vec::new();
        let bytes = self
            .reader
            .read_until(b'\n', &mut line)
            .map_err(|error| format!("read certificate: {error}"))?;
        self.line_number += 1;
        if bytes == 0 {
            return Err(format!(
                "line {}: unexpected end of certificate",
                self.line_number
            ));
        }
        if line.len() > 16 * 1024 * 1024 {
            return Err(format!("line {}: record exceeds 16 MiB", self.line_number));
        }
        if line.last() != Some(&b'\n') || line.ends_with(b"\r\n") {
            return Err(format!(
                "line {}: records require one LF terminator",
                self.line_number
            ));
        }
        if bind {
            update_chain(&mut self.chain, &line);
        }
        self.bytes += line.len() as u64;
        serde_json::from_slice(&line[..line.len() - 1])
            .map_err(|error| format!("line {}: {error}", self.line_number))
    }

    fn finish(mut self) -> Result<u64, String> {
        let mut trailing = Vec::new();
        self.reader
            .read_to_end(&mut trailing)
            .map_err(|error| format!("read trailing certificate bytes: {error}"))?;
        if !trailing.is_empty() {
            return Err("certificate has trailing records or bytes".into());
        }
        Ok(self.bytes)
    }
}

fn validate_header(
    header: &HeaderRecord,
    start: &StartCurve,
    root: &CurveRecord,
    spine: &ModPoly,
    row: &ModPoly,
) -> Result<GridConfig, String> {
    let config = GridConfig {
        side: header.side,
        batch_rows: 1,
        audit_points: header.construction_audit_points,
        audit_seed_x: header.construction_audit_seed_x,
        source_commit: header.source_commit.clone(),
    };
    config.validate()?;
    let expected = make_header(&config, start, root, spine, row);
    if header != &expected {
        return Err("certificate header differs from the frozen P-256 grid header".into());
    }
    Ok(config)
}

fn parse_model(f: &Field, curve_record: &CurveRecord, what: &str) -> Result<(Model, Fe), String> {
    let model = Model {
        a: parse_fe(f, &curve_record.a, &format!("{what} a"))?,
        b: parse_fe(f, &curve_record.b, &format!("{what} b"))?,
    };
    let j = model
        .j(f)
        .ok_or_else(|| format!("{what}: singular curve"))?;
    let recorded_j = parse_fe(f, &curve_record.j, &format!("{what} j"))?;
    if j != recorded_j {
        return Err(format!("{what}: recorded j differs from the model"));
    }
    Ok((model, j))
}

fn verify_curve_record(
    f: &Field,
    start: &StartCurve,
    record: &CurveRecord,
    model: &Model,
    is_root: bool,
    order_is_prime: bool,
    audit_points: usize,
    audit_seed_x: u64,
    what: &str,
) -> Result<(), String> {
    let expected = describe_curve(
        f,
        start,
        model,
        is_root,
        order_is_prime,
        audit_points,
        audit_seed_x,
    )?;
    if record != &expected {
        return Err(format!(
            "{what}: curve identity, generator, order, or detector record differs from replay"
        ));
    }
    Ok(())
}

fn verify_node(
    f: &Field,
    start: &StartCurve,
    side: u32,
    record: NodeRecord,
    phi: &ModPoly,
    source: GridPoint,
    predecessor: Option<Fe>,
    order_is_prime: bool,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<GeneratedNode, String> {
    let coordinate = record.coordinate;
    let expected_index = node_index(side, coordinate.x, coordinate.y);
    if record.record != "node"
        || record.index != expected_index
        || record.parent != node_parent(side, coordinate.x, coordinate.y)
        || record.degree != phi.ell
        || record.predecessor_j != predecessor.map(|value| hex_fe(f, &value))
    {
        return Err(format!(
            "node {}: index, coordinate, parent, degree, or predecessor differs from the grid",
            record.index
        ));
    }
    let selected_j = selected_target(f, phi, &source.model, predecessor)?;
    let (model, j) = parse_model(f, &record.curve, &format!("node {}", record.index))?;
    if j != selected_j {
        return Err(format!(
            "node {}: target j is not the frozen non-backtracking choice",
            record.index
        ));
    }
    let kernel: Vec<Fe> = record
        .kernel_polynomial
        .iter()
        .enumerate()
        .map(|(position, value)| {
            parse_fe(
                f,
                value,
                &format!("node {} kernel coefficient {position}", record.index),
            )
        })
        .collect::<Result<_, _>>()?;
    let velu = kernel::verify_kernel(f, &source.model, phi.ell, &kernel)
        .map_err(|error| format!("node {} kernel: {error:?}", record.index))?;
    let linked_u2 = kernel::link_to_target(f, &velu, &model)
        .map_err(|error| format!("node {} target link: {error:?}", record.index))?;
    let recorded_u2 = parse_fe(
        f,
        &record.isomorphism_u2,
        &format!("node {} isomorphism", record.index),
    )?;
    let (canonical, canonical_u2) = curve::canonical_model(f, &velu)
        .ok_or_else(|| format!("node {} exceptional target", record.index))?;
    if linked_u2 != recorded_u2 || canonical_u2 != recorded_u2 || canonical != model {
        return Err(format!(
            "node {}: target is not the canonical Velu codomain model",
            record.index
        ));
    }
    verify_curve_record(
        f,
        start,
        &record.curve,
        &model,
        false,
        order_is_prime,
        audit_points,
        audit_seed_x,
        &format!("node {}", record.index),
    )?;
    Ok(GeneratedNode { record, model, j })
}

pub fn verify_jsonl<R: BufRead>(
    reader: R,
    audit_points: usize,
    audit_seed_x: u64,
    batch_rows: usize,
) -> Result<VerificationReceipt, String> {
    if audit_points == 0 || batch_rows == 0 {
        return Err("audit points and row batch must be positive".into());
    }
    let started = Instant::now();
    let (start, field, spine_phi, row_phi) = build_context()?;
    let order_is_prime = is_probable_prime(&start.order);
    if !order_is_prime {
        return Err("the registered P-256 group order is not prime".into());
    }
    let mut input = LineReader::new(reader);
    let header: HeaderRecord = input.read(true)?;
    let root_record: RootRecord = input.read(true)?;
    if root_record.record != "root"
        || root_record.index != 0
        || root_record.coordinate != (Coordinate { x: 0, y: 0 })
    {
        return Err("the root record is not grid coordinate (0,0) at index zero".into());
    }
    let (root_model, root_j) = parse_model(&field, &root_record.curve, "root")?;
    if field.to_big(&root_model.a) != start.a || field.to_big(&root_model.b) != start.b {
        return Err("the root is not the registered P-256 model".into());
    }
    verify_curve_record(
        &field,
        &start,
        &root_record.curve,
        &root_model,
        true,
        order_is_prime,
        audit_points,
        audit_seed_x,
        "root",
    )?;
    let config = validate_header(&header, &start, &root_record.curve, &spine_phi, &row_phi)?;
    let mut accumulator = Accumulator::new();
    accumulator.observe(
        0,
        root_record.coordinate,
        &root_model,
        root_j,
        &root_record.curve,
    )?;
    let mut spine = Vec::with_capacity(config.side as usize);
    spine.push(GridPoint {
        model: root_model,
        j: root_j,
    });
    for y in 1..config.side {
        let record: NodeRecord = input.read(true)?;
        if record.coordinate != (Coordinate { x: 0, y }) {
            return Err(format!("spine record {y} has the wrong coordinate"));
        }
        let predecessor = if y == 1 {
            None
        } else {
            Some(spine[(y - 2) as usize].j)
        };
        let checked = verify_node(
            &field,
            &start,
            config.side,
            record,
            &spine_phi,
            spine[(y - 1) as usize],
            predecessor,
            order_is_prime,
            audit_points,
            audit_seed_x,
        )?;
        accumulator.observe(
            checked.record.index,
            checked.record.coordinate,
            &checked.model,
            checked.j,
            &checked.record.curve,
        )?;
        spine.push(GridPoint {
            model: checked.model,
            j: checked.j,
        });
    }
    eprintln!(
        "p256-isogeny-million: replayed {} spine curves",
        spine.len()
    );

    let batch_rows = batch_rows.min(config.side as usize).max(1);
    for batch_start in (0..config.side as usize).step_by(batch_rows) {
        let batch_end = (batch_start + batch_rows).min(config.side as usize);
        let mut encoded_rows = Vec::with_capacity(batch_end - batch_start);
        for y in batch_start..batch_end {
            let mut row = Vec::with_capacity((config.side - 1) as usize);
            for x in 1..config.side {
                let record: NodeRecord = input.read(true)?;
                if record.coordinate != (Coordinate { x, y: y as u32 }) {
                    return Err(format!(
                        "node {}: expected coordinate ({x},{y})",
                        record.index
                    ));
                }
                row.push(record);
            }
            encoded_rows.push(row);
        }
        let checked_rows: Vec<Result<Vec<GeneratedNode>, String>> = encoded_rows
            .into_par_iter()
            .enumerate()
            .map(|(offset, records)| {
                let y = batch_start + offset;
                let mut source = spine[y];
                let mut predecessor = None;
                let mut checked = Vec::with_capacity(records.len());
                for record in records {
                    let node = verify_node(
                        &field,
                        &start,
                        config.side,
                        record,
                        &row_phi,
                        source,
                        predecessor,
                        order_is_prime,
                        audit_points,
                        audit_seed_x,
                    )?;
                    predecessor = Some(source.j);
                    source = GridPoint {
                        model: node.model,
                        j: node.j,
                    };
                    checked.push(node);
                }
                Ok(checked)
            })
            .collect();
        for row in checked_rows {
            for checked in row? {
                accumulator.observe(
                    checked.record.index,
                    checked.record.coordinate,
                    &checked.model,
                    checked.j,
                    &checked.record.curve,
                )?;
            }
        }
        eprintln!(
            "p256-isogeny-million: replayed rows {batch_end}/{}; {} unique curves",
            config.side, accumulator.seen
        );
    }
    let recorded_summary: GridSummary = input.read(false)?;
    let computed_summary = accumulator.summary(input.chain)?;
    if recorded_summary != computed_summary {
        return Err("summary or record-chain digest differs from complete replay".into());
    }
    if computed_summary.unique_curves != config.target_curves() {
        return Err(format!(
            "certificate has {} curves, expected {}",
            computed_summary.unique_curves,
            config.target_curves()
        ));
    }
    let bytes = input.finish()?;
    Ok(VerificationReceipt {
        schema: RECEIPT_SCHEMA.into(),
        status: "pass".into(),
        source_commit: config.source_commit,
        audit_points,
        audit_seed_x,
        elapsed_ms: started.elapsed().as_millis(),
        resources: resource_usage(),
        uncompressed_bytes: bytes,
        summary: computed_summary,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn small_certificate() -> Vec<u8> {
        let mut bytes = Vec::new();
        generate_jsonl(
            &mut bytes,
            &GridConfig {
                side: 3,
                batch_rows: 2,
                audit_points: 1,
                audit_seed_x: 0,
                source_commit: "0".repeat(40),
            },
        )
        .expect("generate small grid");
        bytes
    }

    fn mutate(bytes: &[u8], line: usize, path: &[&str], value: serde_json::Value) -> Vec<u8> {
        let mut lines: Vec<Vec<u8>> = bytes
            .split(|byte| *byte == b'\n')
            .filter(|line| !line.is_empty())
            .map(<[u8]>::to_vec)
            .collect();
        let mut record: serde_json::Value =
            serde_json::from_slice(&lines[line]).expect("record parses");
        let mut selected = &mut record;
        for key in &path[..path.len() - 1] {
            selected = &mut selected[*key];
        }
        selected[path[path.len() - 1]] = value;
        lines[line] = serde_json::to_vec(&record).expect("record serializes");
        let mut out = Vec::new();
        for mut line in lines {
            out.append(&mut line);
            out.push(b'\n');
        }
        out
    }

    #[test]
    fn small_grid_replays_and_tampering_is_rejected() {
        let bytes = small_certificate();
        let receipt = verify_jsonl(Cursor::new(&bytes), 2, 7, 2).expect("replay small grid");
        assert_eq!(receipt.summary.unique_curves, 9);
        assert_eq!(receipt.summary.parent_edges, 8);

        let cases = [
            mutate(
                &bytes,
                2,
                &["curve", "j"],
                serde_json::Value::String("0".repeat(64)),
            ),
            mutate(&bytes, 2, &["parent"], serde_json::Value::from(7)),
            mutate(&bytes, 2, &["coordinate", "y"], serde_json::Value::from(2)),
            mutate(
                &bytes,
                2,
                &["curve", "curve_uid"],
                serde_json::Value::String("urn:ec-record:1:sha256:00".into()),
            ),
            mutate(
                &bytes,
                2,
                &["curve", "detectors", "qr_prefix_64"],
                serde_json::Value::from(64),
            ),
        ];
        for tampered in cases {
            assert!(verify_jsonl(Cursor::new(tampered), 2, 7, 2).is_err());
        }

        let mut kernel_tamper: serde_json::Value = serde_json::from_slice(
            bytes
                .split(|byte| *byte == b'\n')
                .filter(|line| !line.is_empty())
                .nth(2)
                .unwrap(),
        )
        .unwrap();
        kernel_tamper["kernel_polynomial"][0] = serde_json::Value::String("0".repeat(64));
        let changed = mutate(
            &bytes,
            2,
            &["kernel_polynomial"],
            kernel_tamper["kernel_polynomial"].clone(),
        );
        assert!(verify_jsonl(Cursor::new(changed), 2, 7, 2).is_err());

        let mut duplicate: serde_json::Value = serde_json::from_slice(
            bytes
                .split(|byte| *byte == b'\n')
                .filter(|line| !line.is_empty())
                .nth(2)
                .unwrap(),
        )
        .unwrap();
        let later: serde_json::Value = serde_json::from_slice(
            bytes
                .split(|byte| *byte == b'\n')
                .filter(|line| !line.is_empty())
                .nth(3)
                .unwrap(),
        )
        .unwrap();
        duplicate["index"] = later["index"].clone();
        duplicate["coordinate"] = later["coordinate"].clone();
        duplicate["parent"] = later["parent"].clone();
        let changed = mutate(&bytes, 3, &["curve"], duplicate["curve"].clone());
        assert!(verify_jsonl(Cursor::new(changed), 2, 7, 2).is_err());
    }
}
