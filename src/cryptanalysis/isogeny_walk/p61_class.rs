//! Bounded-memory certificate stream for the selected complete p61 class.
//!
//! The certified class-group target is
//! `Z/179047152 x Z/2`.  Degree 13 generates a subgroup of order
//! 19,894,128; degree 37 adds a relative factor 9 and degree 59 adds a
//! relative factor 2.  The producer therefore walks eighteen independent
//! degree-13 lanes.  A segment retains one current model and predecessor j
//! per lane, plus a bounded batch of records.  Checkpoints make every segment
//! restartable without retaining the full graph or materializing root paths.

use std::io::{BufRead, Write};
use std::time::Instant;

use num_bigint::BigUint;
use rayon::prelude::*;
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use super::curve::{self, Model, OrderAudit};
use super::field::{Fe, Field};
use super::kernel;
use super::modpoly::{self, ModPoly};
use super::poly::{self, Poly};
use super::walk::{ClassInfo, EllKind, StartCurve};
use crate::hash::sha256::sha256;

pub const SEGMENT_SCHEMA: &str = "ec-hardness/p61-class-segment/v1";
pub const SUMMARY_SCHEMA: &str = "ec-hardness/p61-class-segment-summary/v1";
pub const CHECKPOINT_SCHEMA: &str = "ec-hardness/p61-class-checkpoint/v1";
pub const RECEIPT_SCHEMA: &str = "ec-hardness/p61-class-segment-receipt/v1";

pub const ROW_DEGREE: u64 = 13;
pub const COSET_DEGREE: u64 = 37;
pub const TOP_DEGREE: u64 = 59;
pub const ROOT_SEED: u64 = 20261014;
pub const ROW_ORDER: u64 = 19_894_128;
pub const COSET_ROWS: u8 = 9;
pub const TOP_ROWS: u8 = 2;
pub const LANES: usize = COSET_ROWS as usize * TOP_ROWS as usize;
pub const EXPECTED_CLASS_VERTICES: u64 = 358_094_304;

const P: &str = "2305843009213693951";
const A: &str = "1055";
const B: &str = "16361";
const GROUP_ORDER: &str = "2305843009918553405";
const SUBGROUP_ORDER: &str = "6883113462443443";
const COFACTOR: &str = "335";
const GX: &str = "843688558021621073";
const GY: &str = "105039246913253664";
const TRACE: i64 = -704_859_453;
const D_K: &str = "-8726545188371316595";
const CLASS_GROUP_FACTORS: [u64; 2] = [179_047_152, 2];
const ROW_CLASS_VECTOR: [u64; 2] = [81_997_155, 0];
const COSET_CLASS_VECTOR: [u64; 2] = [36_835_963, 0];
const TOP_CLASS_VECTOR: [u64; 2] = [315_398, 1];

#[derive(Clone, Debug)]
pub struct SegmentConfig {
    pub x_start: u64,
    pub length: u64,
    pub batch_steps: usize,
    pub audit_points: usize,
    pub audit_seed_x: u64,
    pub source_commit: String,
    pub class_number_receipt_sha256: String,
    pub previous_checkpoint_sha256: Option<String>,
}

impl SegmentConfig {
    pub fn validate(&self) -> Result<(), String> {
        if self.length == 0 {
            return Err("segment length must be positive".into());
        }
        let end = self
            .x_start
            .checked_add(self.length)
            .ok_or("segment range overflows u64")?;
        if end > ROW_ORDER {
            return Err(format!(
                "segment end {end} exceeds frozen row order {ROW_ORDER}"
            ));
        }
        if self.batch_steps == 0 {
            return Err("batch steps must be positive".into());
        }
        if self.audit_points == 0 {
            return Err("audit point count must be positive".into());
        }
        if !is_lower_hex(&self.source_commit, 40) {
            return Err("source commit must be 40 lowercase hexadecimal digits".into());
        }
        if !is_lower_hex(&self.class_number_receipt_sha256, 64) {
            return Err(
                "class-number receipt digest must be 64 lowercase hexadecimal digits".into(),
            );
        }
        match (self.x_start, &self.previous_checkpoint_sha256) {
            (0, None) => {}
            (0, Some(_)) => {
                return Err("the first segment cannot have a predecessor checkpoint".into())
            }
            (_, None) => {
                return Err("continuation segments require a predecessor checkpoint".into())
            }
            (_, Some(digest)) if !is_lower_hex(digest, 64) => {
                return Err("predecessor checkpoint digest must be lowercase SHA-256".into())
            }
            _ => {}
        }
        Ok(())
    }

    pub fn x_end(&self) -> u64 {
        self.x_start + self.length
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Coordinate {
    pub x: u64,
    pub y: u8,
    pub z: u8,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CurveRecord {
    pub a: String,
    pub b: String,
    pub j: String,
    pub order_status: String,
    pub qr_prefix_64: u32,
    pub a_signed_bits: u64,
    pub b_signed_bits: u64,
    pub a_minus_3_model: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EdgeCertificate {
    pub parent: Coordinate,
    pub degree: u64,
    pub predecessor_j: Option<String>,
    pub kernel_polynomial: Vec<String>,
    pub isomorphism_u2: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NodeRecord {
    pub record: String,
    pub index: u64,
    pub coordinate: Coordinate,
    pub edge: Option<EdgeCertificate>,
    pub curve: CurveRecord,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClosureRecord {
    pub record: String,
    pub source: Coordinate,
    pub target: Coordinate,
    pub degree: u64,
    pub predecessor_j: String,
    pub kernel_polynomial: Vec<String>,
    pub isomorphism_u2: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SegmentHeader {
    pub record: String,
    pub schema: String,
    pub source_commit: String,
    pub class_number_receipt_sha256: String,
    pub previous_checkpoint_sha256: Option<String>,
    pub field_p: String,
    pub curve_a: String,
    pub curve_b: String,
    pub group_order: String,
    pub subgroup_order: String,
    pub cofactor: String,
    pub trace: i64,
    pub fundamental_discriminant: String,
    pub class_group_factors: [u64; 2],
    pub row_class_vector: [u64; 2],
    pub coset_class_vector: [u64; 2],
    pub top_class_vector: [u64; 2],
    pub row_order: u64,
    pub coset_rows: u8,
    pub top_rows: u8,
    pub expected_class_vertices: u64,
    pub x_start: u64,
    pub x_end: u64,
    pub batch_steps: usize,
    pub construction_audit_points: usize,
    pub construction_audit_seed_x: u64,
    pub row_degree: u64,
    pub coset_degree: u64,
    pub top_degree: u64,
    pub row_phi_checks: usize,
    pub coset_phi_checks: usize,
    pub top_phi_checks: usize,
    pub canonical_model_rule: String,
    pub selection_rule: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SegmentSummary {
    pub record: String,
    pub schema: String,
    pub x_start: u64,
    pub x_end: u64,
    pub node_records: u64,
    pub parent_edges: u64,
    pub closure_edges: u64,
    pub order_consistent_unproved: u64,
    pub record_chain_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LaneCheckpoint {
    pub y: u8,
    pub z: u8,
    pub predecessor_j: Option<String>,
    pub current: CurveRecord,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Checkpoint {
    pub schema: String,
    pub source_commit: String,
    pub class_number_receipt_sha256: String,
    pub previous_checkpoint_sha256: Option<String>,
    pub segment_sha256: String,
    pub x_start: u64,
    pub next_x: u64,
    pub record_chain_sha256: String,
    pub lanes: Vec<LaneCheckpoint>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenerationReceipt {
    pub schema: String,
    pub status: String,
    pub source_commit: String,
    pub elapsed_ms: u128,
    pub uncompressed_bytes: u64,
    pub summary: SegmentSummary,
    pub checkpoint: Checkpoint,
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
    pub uncompressed_bytes: u64,
    pub summary: SegmentSummary,
}

#[derive(Clone, Copy)]
pub(crate) struct Point {
    pub(crate) model: Model,
    pub(crate) j: Fe,
}

#[derive(Clone)]
struct LaneState {
    y: u8,
    z: u8,
    predecessor_j: Option<Fe>,
    current: Point,
}

pub(crate) struct Context {
    pub(crate) start: StartCurve,
    pub(crate) field: Field,
    pub(crate) row_phi: ModPoly,
    pub(crate) coset_phi: ModPoly,
    pub(crate) top_phi: ModPoly,
}

struct GeneratedNode {
    record: NodeRecord,
    point: Point,
}

pub(crate) fn is_lower_hex(value: &str, length: usize) -> bool {
    value.len() == length
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn parse_big(value: &str, what: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.as_bytes(), 10).ok_or_else(|| format!("invalid {what}"))
}

fn start_curve() -> Result<StartCurve, String> {
    Ok(StartCurve {
        name: "p61-i0025".into(),
        tag: "fp".into(),
        p: parse_big(P, "p")?,
        a: parse_big(A, "a")?,
        b: parse_big(B, "b")?,
        order: parse_big(GROUP_ORDER, "group order")?,
        subgroup_order: parse_big(SUBGROUP_ORDER, "subgroup order")?,
        cofactor: parse_big(COFACTOR, "cofactor")?,
        gx: parse_big(GX, "generator x")?,
        gy: parse_big(GY, "generator y")?,
    })
}

pub(crate) fn build_context() -> Result<Context, String> {
    if ROW_ORDER
        .checked_mul(u64::from(COSET_ROWS))
        .and_then(|value| value.checked_mul(u64::from(TOP_ROWS)))
        != Some(EXPECTED_CLASS_VERTICES)
    {
        return Err("frozen lattice dimensions do not multiply to the class number".into());
    }
    let start = start_curve()?;
    let field = Field::new(&start.p).ok_or("p61 field failed primality validation")?;
    let root = Model {
        a: field.from_big(&start.a),
        b: field.from_big(&start.b),
    };
    let gx = field.from_big(&start.gx);
    let gy = field.from_big(&start.gy);
    if !root.on_curve(&field, &gx, &gy)
        || curve::scalar_mul(&field, &root, &gx, &gy, &start.subgroup_order).is_some()
    {
        return Err("frozen root generator failed validation".into());
    }
    let class = ClassInfo::compute(&start, &[ROW_DEGREE, COSET_DEGREE, TOP_DEGREE]);
    for degree in [ROW_DEGREE, COSET_DEGREE, TOP_DEGREE] {
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
    let built: Vec<Result<ModPoly, String>> = [ROW_DEGREE, COSET_DEGREE, TOP_DEGREE]
        .par_iter()
        .map(|degree| modpoly::modular_polynomial(&field, *degree))
        .collect();
    let mut built = built.into_iter();
    Ok(Context {
        start,
        field,
        row_phi: built.next().expect("three modular polynomials")?,
        coset_phi: built.next().expect("three modular polynomials")?,
        top_phi: built.next().expect("three modular polynomials")?,
    })
}

pub(crate) fn hex_fe(field: &Field, value: &Fe) -> String {
    format!("{:016x}", field.to_big(value))
}

fn parse_fe(field: &Field, value: &str, what: &str) -> Result<Fe, String> {
    if !is_lower_hex(value, 16) {
        return Err(format!("{what}: expected 16 lowercase hexadecimal digits"));
    }
    let integer = BigUint::parse_bytes(value.as_bytes(), 16)
        .ok_or_else(|| format!("{what}: invalid hexadecimal field element"))?;
    if &integer >= field.modulus() {
        return Err(format!("{what}: field element is not canonical"));
    }
    Ok(field.from_big(&integer))
}

pub(crate) fn coefficient_bits(field: &Field, value: &Fe) -> u64 {
    let positive = field.to_big(value);
    let negative = field.modulus() - &positive;
    positive.min(negative).bits()
}

fn describe_curve(
    context: &Context,
    model: &Model,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<CurveRecord, String> {
    let field = &context.field;
    let j = model.j(field).ok_or("singular curve")?;
    let audit = curve::audit_order(
        field,
        model,
        &context.start.order,
        false,
        audit_points,
        audit_seed_x,
    );
    if audit == OrderAudit::Refuted {
        return Err("order audit refuted the frozen group order".into());
    }
    Ok(CurveRecord {
        a: hex_fe(field, &model.a),
        b: hex_fe(field, &model.b),
        j: hex_fe(field, &j),
        order_status: "consistent_unproved".into(),
        qr_prefix_64: curve::qr_prefix_64(field, model),
        a_signed_bits: coefficient_bits(field, &model.a),
        b_signed_bits: coefficient_bits(field, &model.b),
        a_minus_3_model: !curve::a_minus_3_models(field, model).is_empty(),
    })
}

pub(crate) fn parse_curve(
    context: &Context,
    record: &CurveRecord,
    what: &str,
) -> Result<Point, String> {
    let model = Model {
        a: parse_fe(&context.field, &record.a, &format!("{what} a"))?,
        b: parse_fe(&context.field, &record.b, &format!("{what} b"))?,
    };
    let j = model
        .j(&context.field)
        .ok_or_else(|| format!("{what}: singular curve"))?;
    if j != parse_fe(&context.field, &record.j, &format!("{what} j"))? {
        return Err(format!("{what}: recorded j differs from the model"));
    }
    Ok(Point { model, j })
}

fn rational_roots(context: &Context, phi: &ModPoly, model: &Model) -> Result<[Fe; 2], String> {
    let j = model.j(&context.field).ok_or("singular source curve")?;
    let roots = poly::roots(
        &context.field,
        &phi.at_x(&context.field, &j),
        ROOT_SEED ^ phi.ell,
    );
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
    context: &Context,
    phi: &ModPoly,
    source: &Model,
    predecessor: Option<Fe>,
) -> Result<Fe, String> {
    let roots = rational_roots(context, phi, source)?;
    match predecessor {
        None => Ok(roots[0]),
        Some(previous) if roots[0] == previous => Ok(roots[1]),
        Some(previous) if roots[1] == previous => Ok(roots[0]),
        Some(_) => Err(format!("predecessor is not a degree-{} neighbour", phi.ell)),
    }
}

fn make_step(
    context: &Context,
    phi: &ModPoly,
    source: Point,
    predecessor: Option<Fe>,
) -> Result<(Point, Poly, Fe), String> {
    let target_j = selected_target(context, phi, &source.model, predecessor)?;
    let isogeny = kernel::elkies_isogeny(&context.field, &source.model, phi, &target_j)
        .map_err(|error| format!("degree-{} edge construction: {error:?}", phi.ell))?;
    let (target, canonical_u2) = curve::canonical_model(&context.field, &isogeny.codomain)
        .ok_or("isogeny reached an exceptional target j-invariant")?;
    let linked_u2 = kernel::link_to_target(&context.field, &isogeny.codomain, &target)
        .map_err(|error| format!("degree-{} target link: {error:?}", phi.ell))?;
    if canonical_u2 != linked_u2 || target.j(&context.field) != Some(target_j) {
        return Err("canonical target or isomorphism changed the selected j-invariant".into());
    }
    Ok((
        Point {
            model: target,
            j: target_j,
        },
        isogeny.kernel,
        linked_u2,
    ))
}

pub fn coordinate_index(coordinate: Coordinate) -> Result<u64, String> {
    if coordinate.x >= ROW_ORDER || coordinate.y >= COSET_ROWS || coordinate.z >= TOP_ROWS {
        return Err("coordinate lies outside the frozen class lattice".into());
    }
    Ok(
        (u64::from(coordinate.z) * u64::from(COSET_ROWS) + u64::from(coordinate.y)) * ROW_ORDER
            + coordinate.x,
    )
}

fn node_record(
    context: &Context,
    coordinate: Coordinate,
    parent: Coordinate,
    phi: &ModPoly,
    source: Point,
    predecessor: Option<Fe>,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<GeneratedNode, String> {
    let (point, kernel, u2) = make_step(context, phi, source, predecessor)?;
    Ok(GeneratedNode {
        record: NodeRecord {
            record: "node".into(),
            index: coordinate_index(coordinate)?,
            coordinate,
            edge: Some(EdgeCertificate {
                parent,
                degree: phi.ell,
                predecessor_j: predecessor.map(|value| hex_fe(&context.field, &value)),
                kernel_polynomial: kernel
                    .iter()
                    .map(|value| hex_fe(&context.field, value))
                    .collect(),
                isomorphism_u2: hex_fe(&context.field, &u2),
            }),
            curve: describe_curve(context, &point.model, audit_points, audit_seed_x)?,
        },
        point,
    })
}

fn root_record(
    context: &Context,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<GeneratedNode, String> {
    let point = Point {
        model: Model {
            a: context.field.from_big(&context.start.a),
            b: context.field.from_big(&context.start.b),
        },
        j: Model {
            a: context.field.from_big(&context.start.a),
            b: context.field.from_big(&context.start.b),
        }
        .j(&context.field)
        .ok_or("frozen root is singular")?,
    };
    let coordinate = Coordinate { x: 0, y: 0, z: 0 };
    Ok(GeneratedNode {
        record: NodeRecord {
            record: "node".into(),
            index: 0,
            coordinate,
            edge: None,
            curve: describe_curve(context, &point.model, audit_points, audit_seed_x)?,
        },
        point,
    })
}

fn lane_seeds(
    context: &Context,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<Vec<GeneratedNode>, String> {
    let root = root_record(context, audit_points, audit_seed_x)?;
    let mut lanes = Vec::with_capacity(LANES);
    lanes.push(root);
    for z in 0..TOP_ROWS {
        let mut source = if z == 0 {
            lanes[0].point
        } else {
            let generated = node_record(
                context,
                Coordinate { x: 0, y: 0, z },
                Coordinate { x: 0, y: 0, z: 0 },
                &context.top_phi,
                lanes[0].point,
                None,
                audit_points,
                audit_seed_x,
            )?;
            let point = generated.point;
            lanes.push(generated);
            point
        };
        let mut predecessor = None;
        for y in 1..COSET_ROWS {
            let generated = node_record(
                context,
                Coordinate { x: 0, y, z },
                Coordinate { x: 0, y: y - 1, z },
                &context.coset_phi,
                source,
                predecessor,
                audit_points,
                audit_seed_x,
            )?;
            predecessor = Some(source.j);
            source = generated.point;
            lanes.push(generated);
        }
    }
    lanes.sort_by_key(|lane| (lane.record.coordinate.z, lane.record.coordinate.y));
    if lanes.len() != LANES {
        return Err("lane seed construction produced the wrong count".into());
    }
    Ok(lanes)
}

fn header(config: &SegmentConfig, context: &Context) -> SegmentHeader {
    SegmentHeader {
        record: "header".into(),
        schema: SEGMENT_SCHEMA.into(),
        source_commit: config.source_commit.clone(),
        class_number_receipt_sha256: config.class_number_receipt_sha256.clone(),
        previous_checkpoint_sha256: config.previous_checkpoint_sha256.clone(),
        field_p: P.into(),
        curve_a: A.into(),
        curve_b: B.into(),
        group_order: GROUP_ORDER.into(),
        subgroup_order: SUBGROUP_ORDER.into(),
        cofactor: COFACTOR.into(),
        trace: TRACE,
        fundamental_discriminant: D_K.into(),
        class_group_factors: CLASS_GROUP_FACTORS,
        row_class_vector: ROW_CLASS_VECTOR,
        coset_class_vector: COSET_CLASS_VECTOR,
        top_class_vector: TOP_CLASS_VECTOR,
        row_order: ROW_ORDER,
        coset_rows: COSET_ROWS,
        top_rows: TOP_ROWS,
        expected_class_vertices: EXPECTED_CLASS_VERTICES,
        x_start: config.x_start,
        x_end: config.x_end(),
        batch_steps: config.batch_steps,
        construction_audit_points: config.audit_points,
        construction_audit_seed_x: config.audit_seed_x,
        row_degree: ROW_DEGREE,
        coset_degree: COSET_DEGREE,
        top_degree: TOP_DEGREE,
        row_phi_checks: context.row_phi.checks_passed,
        coset_phi_checks: context.coset_phi.checks_passed,
        top_phi_checks: context.top_phi.checks_passed,
        canonical_model_rule: curve::CANONICAL_RULE.into(),
        selection_rule: "first edge chooses the numerically smaller rational root; later edges choose the root unequal to the predecessor j".into(),
    }
}

fn validate_checkpoint(
    config: &SegmentConfig,
    checkpoint: &Checkpoint,
    checkpoint_sha256: &str,
) -> Result<(), String> {
    if checkpoint.schema != CHECKPOINT_SCHEMA
        || checkpoint.source_commit != config.source_commit
        || checkpoint.class_number_receipt_sha256 != config.class_number_receipt_sha256
        || checkpoint.next_x != config.x_start
        || checkpoint.lanes.len() != LANES
        || config.previous_checkpoint_sha256.as_deref() != Some(checkpoint_sha256)
        || !is_lower_hex(&checkpoint.segment_sha256, 64)
        || !is_lower_hex(&checkpoint.record_chain_sha256, 64)
    {
        return Err("predecessor checkpoint does not match the requested segment".into());
    }
    for (index, lane) in checkpoint.lanes.iter().enumerate() {
        let expected_z = (index / COSET_ROWS as usize) as u8;
        let expected_y = (index % COSET_ROWS as usize) as u8;
        if lane.y != expected_y || lane.z != expected_z {
            return Err("predecessor checkpoint lane order differs from the frozen order".into());
        }
    }
    Ok(())
}

fn states_from_checkpoint(
    context: &Context,
    checkpoint: &Checkpoint,
) -> Result<Vec<LaneState>, String> {
    checkpoint
        .lanes
        .iter()
        .map(|lane| {
            Ok(LaneState {
                y: lane.y,
                z: lane.z,
                predecessor_j: lane
                    .predecessor_j
                    .as_deref()
                    .map(|value| parse_fe(&context.field, value, "checkpoint predecessor j"))
                    .transpose()?,
                current: parse_curve(context, &lane.current, "checkpoint current curve")?,
            })
        })
        .collect()
}

fn checkpoint_from_states(
    config: &SegmentConfig,
    states: &[LaneState],
    context: &Context,
    chain: [u8; 32],
) -> Result<Checkpoint, String> {
    let lanes = states
        .iter()
        .map(|lane| {
            Ok(LaneCheckpoint {
                y: lane.y,
                z: lane.z,
                predecessor_j: lane
                    .predecessor_j
                    .map(|value| hex_fe(&context.field, &value)),
                current: describe_curve(
                    context,
                    &lane.current.model,
                    config.audit_points,
                    config.audit_seed_x,
                )?,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    Ok(Checkpoint {
        schema: CHECKPOINT_SCHEMA.into(),
        source_commit: config.source_commit.clone(),
        class_number_receipt_sha256: config.class_number_receipt_sha256.clone(),
        previous_checkpoint_sha256: config.previous_checkpoint_sha256.clone(),
        segment_sha256: "0".repeat(64),
        x_start: config.x_start,
        next_x: config.x_end(),
        record_chain_sha256: hex::encode(chain),
        lanes,
    })
}

fn update_chain(chain: &mut [u8; 32], line: &[u8]) {
    let mut input = Vec::with_capacity(32 + line.len());
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
            .map_err(|error| format!("write segment: {error}"))?;
        self.bytes = self
            .bytes
            .checked_add(line.len() as u64)
            .ok_or("segment byte count overflow")?;
        Ok(())
    }
}

fn advance_lane(
    context: &Context,
    state: &mut LaneState,
    x_start: u64,
    x_end: u64,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<Vec<NodeRecord>, String> {
    let mut records = Vec::with_capacity((x_end - x_start) as usize);
    for x in x_start..x_end {
        let parent = Coordinate {
            x: x - 1,
            y: state.y,
            z: state.z,
        };
        let generated = node_record(
            context,
            Coordinate {
                x,
                y: state.y,
                z: state.z,
            },
            parent,
            &context.row_phi,
            state.current,
            state.predecessor_j,
            audit_points,
            audit_seed_x,
        )?;
        state.predecessor_j = Some(state.current.j);
        state.current = generated.point;
        records.push(generated.record);
    }
    Ok(records)
}

fn closure_record(
    context: &Context,
    state: &LaneState,
    seed: &GeneratedNode,
) -> Result<ClosureRecord, String> {
    let predecessor = state
        .predecessor_j
        .ok_or("final row state lacks a predecessor")?;
    let (target, kernel, u2) =
        make_step(context, &context.row_phi, state.current, Some(predecessor))?;
    if target.j != seed.point.j || target.model != seed.point.model {
        return Err(format!(
            "degree-{ROW_DEGREE} row ({},{}) does not close at x={ROW_ORDER}",
            state.y, state.z
        ));
    }
    Ok(ClosureRecord {
        record: "closure".into(),
        source: Coordinate {
            x: ROW_ORDER - 1,
            y: state.y,
            z: state.z,
        },
        target: Coordinate {
            x: 0,
            y: state.y,
            z: state.z,
        },
        degree: ROW_DEGREE,
        predecessor_j: hex_fe(&context.field, &predecessor),
        kernel_polynomial: kernel
            .iter()
            .map(|value| hex_fe(&context.field, value))
            .collect(),
        isomorphism_u2: hex_fe(&context.field, &u2),
    })
}

pub fn generate_segment_jsonl<W: Write>(
    writer: &mut W,
    config: &SegmentConfig,
    previous: Option<(&Checkpoint, &str)>,
) -> Result<GenerationReceipt, String> {
    config.validate()?;
    let started = Instant::now();
    let context = build_context()?;
    let mut output = LineWriter::new(writer);
    output.write(&header(config, &context), true)?;

    let mut parent_edges = 0u64;
    let mut states;
    if config.x_start == 0 {
        if previous.is_some() {
            return Err("the first segment cannot consume a checkpoint".into());
        }
        let seeds = lane_seeds(&context, config.audit_points, config.audit_seed_x)?;
        for seed in &seeds {
            parent_edges += u64::from(seed.record.edge.is_some());
            output.write(&seed.record, true)?;
        }
        states = seeds
            .iter()
            .map(|seed| LaneState {
                y: seed.record.coordinate.y,
                z: seed.record.coordinate.z,
                predecessor_j: None,
                current: seed.point,
            })
            .collect();
    } else {
        let (checkpoint, digest) = previous.ok_or("continuation segment lacks its checkpoint")?;
        validate_checkpoint(config, checkpoint, digest)?;
        states = states_from_checkpoint(&context, checkpoint)?;
    }

    let first_new_x = if config.x_start == 0 {
        1
    } else {
        config.x_start
    };
    let mut block_start = first_new_x;
    while block_start < config.x_end() {
        let block_end = config.x_end().min(block_start + config.batch_steps as u64);
        let rows: Vec<Result<Vec<NodeRecord>, String>> = states
            .par_iter_mut()
            .map(|state| {
                advance_lane(
                    &context,
                    state,
                    block_start,
                    block_end,
                    config.audit_points,
                    config.audit_seed_x,
                )
            })
            .collect();
        for row in rows {
            for record in row? {
                parent_edges += 1;
                output.write(&record, true)?;
            }
        }
        eprintln!(
            "p61-isogeny-class: segment {}..{} advanced lanes through x={block_end}",
            config.x_start,
            config.x_end()
        );
        block_start = block_end;
    }

    let mut closure_edges = 0u64;
    if config.x_end() == ROW_ORDER {
        let seeds = lane_seeds(&context, config.audit_points, config.audit_seed_x)?;
        for (state, seed) in states.iter().zip(seeds.iter()) {
            output.write(&closure_record(&context, state, seed)?, true)?;
            closure_edges += 1;
        }
    }

    let node_records = config
        .length
        .checked_mul(LANES as u64)
        .ok_or("node record count overflow")?;
    let summary = SegmentSummary {
        record: "summary".into(),
        schema: SUMMARY_SCHEMA.into(),
        x_start: config.x_start,
        x_end: config.x_end(),
        node_records,
        parent_edges,
        closure_edges,
        order_consistent_unproved: node_records,
        record_chain_sha256: hex::encode(output.chain),
    };
    output.write(&summary, false)?;
    let checkpoint = checkpoint_from_states(config, &states, &context, output.chain)?;
    Ok(GenerationReceipt {
        schema: RECEIPT_SCHEMA.into(),
        status: "pass".into(),
        source_commit: config.source_commit.clone(),
        elapsed_ms: started.elapsed().as_millis(),
        uncompressed_bytes: output.bytes,
        summary,
        checkpoint,
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
        let read = self
            .reader
            .read_until(b'\n', &mut line)
            .map_err(|error| format!("read segment: {error}"))?;
        self.line_number += 1;
        if read == 0 {
            return Err(format!(
                "line {}: unexpected end of segment",
                self.line_number
            ));
        }
        if line.last() != Some(&b'\n') || line.ends_with(b"\r\n") {
            return Err(format!("line {}: records require one LF", self.line_number));
        }
        if line.len() > 4 * 1024 * 1024 {
            return Err(format!("line {}: record exceeds 4 MiB", self.line_number));
        }
        if bind {
            update_chain(&mut self.chain, &line);
        }
        self.bytes = self
            .bytes
            .checked_add(line.len() as u64)
            .ok_or("segment byte count overflow")?;
        serde_json::from_slice(&line).map_err(|error| format!("line {}: {error}", self.line_number))
    }

    fn finish(mut self) -> Result<u64, String> {
        let mut trailing = Vec::new();
        self.reader
            .read_to_end(&mut trailing)
            .map_err(|error| format!("read segment tail: {error}"))?;
        if !trailing.is_empty() {
            return Err("segment has trailing bytes".into());
        }
        Ok(self.bytes)
    }
}

fn verify_curve(
    context: &Context,
    record: &CurveRecord,
    audit_points: usize,
    audit_seed_x: u64,
    what: &str,
) -> Result<Point, String> {
    let point = parse_curve(context, record, what)?;
    let expected = describe_curve(context, &point.model, audit_points, audit_seed_x)?;
    if record != &expected {
        return Err(format!(
            "{what}: detector or order record differs from replay"
        ));
    }
    Ok(point)
}

fn verify_edge(
    context: &Context,
    phi: &ModPoly,
    source: Point,
    predecessor: Option<Fe>,
    record: &NodeRecord,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<Point, String> {
    let edge = record
        .edge
        .as_ref()
        .ok_or("non-root node lacks an edge certificate")?;
    if edge.degree != phi.ell
        || edge.predecessor_j != predecessor.map(|value| hex_fe(&context.field, &value))
    {
        return Err(format!("node {} edge metadata differs", record.index));
    }
    let expected_j = selected_target(context, phi, &source.model, predecessor)?;
    let point = verify_curve(
        context,
        &record.curve,
        audit_points,
        audit_seed_x,
        &format!("node {}", record.index),
    )?;
    if point.j != expected_j {
        return Err(format!(
            "node {} target is not the selected neighbour",
            record.index
        ));
    }
    let kernel: Vec<Fe> = edge
        .kernel_polynomial
        .iter()
        .enumerate()
        .map(|(index, value)| {
            parse_fe(
                &context.field,
                value,
                &format!("node {} kernel coefficient {index}", record.index),
            )
        })
        .collect::<Result<_, _>>()?;
    let velu = kernel::verify_kernel(&context.field, &source.model, phi.ell, &kernel)
        .map_err(|error| format!("node {} kernel: {error:?}", record.index))?;
    let linked = kernel::link_to_target(&context.field, &velu, &point.model)
        .map_err(|error| format!("node {} target link: {error:?}", record.index))?;
    let recorded_u2 = parse_fe(
        &context.field,
        &edge.isomorphism_u2,
        &format!("node {} isomorphism", record.index),
    )?;
    let (canonical, canonical_u2) = curve::canonical_model(&context.field, &velu)
        .ok_or_else(|| format!("node {} exceptional target", record.index))?;
    if linked != recorded_u2 || canonical_u2 != recorded_u2 || canonical != point.model {
        return Err(format!(
            "node {} canonical target link differs",
            record.index
        ));
    }
    Ok(point)
}

fn verify_seed_records<R: BufRead>(
    input: &mut LineReader<R>,
    context: &Context,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<Vec<LaneState>, String> {
    let expected = lane_seeds(context, 1, 0)?;
    let mut states = Vec::with_capacity(LANES);
    let mut replayed = Vec::with_capacity(LANES);
    for expected_seed in expected {
        let record: NodeRecord = input.read(true)?;
        if record.record != "node"
            || record.index != expected_seed.record.index
            || record.coordinate != expected_seed.record.coordinate
            || record
                .edge
                .as_ref()
                .map(|edge| (&edge.parent, edge.degree, &edge.predecessor_j))
                != expected_seed
                    .record
                    .edge
                    .as_ref()
                    .map(|edge| (&edge.parent, edge.degree, &edge.predecessor_j))
        {
            return Err("lane seed coordinate or parent metadata differs".into());
        }
        let point = if let Some(edge) = &record.edge {
            let parent = replayed
                .iter()
                .find(|(coordinate, _): &&(Coordinate, Point)| *coordinate == edge.parent)
                .map(|(_, point)| *point)
                .ok_or("lane seed parent was not replayed first")?;
            let predecessor = if edge.degree == COSET_DEGREE {
                if record.coordinate.y < 2 {
                    None
                } else {
                    let predecessor_coordinate = Coordinate {
                        x: 0,
                        y: record.coordinate.y - 2,
                        z: record.coordinate.z,
                    };
                    Some(
                        replayed
                            .iter()
                            .find(|(coordinate, _)| *coordinate == predecessor_coordinate)
                            .map(|(_, point)| point.j)
                            .ok_or("lane seed predecessor was not replayed first")?,
                    )
                }
            } else {
                None
            };
            let phi = if edge.degree == COSET_DEGREE {
                &context.coset_phi
            } else if edge.degree == TOP_DEGREE {
                &context.top_phi
            } else {
                return Err("lane seed uses an unexpected degree".into());
            };
            verify_edge(
                context,
                phi,
                parent,
                predecessor,
                &record,
                audit_points,
                audit_seed_x,
            )?
        } else {
            if record.coordinate != (Coordinate { x: 0, y: 0, z: 0 }) {
                return Err("only the frozen root may omit a parent edge".into());
            }
            let point = verify_curve(context, &record.curve, audit_points, audit_seed_x, "root")?;
            if context.field.to_big(&point.model.a) != context.start.a
                || context.field.to_big(&point.model.b) != context.start.b
            {
                return Err("segment root differs from the frozen root".into());
            }
            point
        };
        replayed.push((record.coordinate, point));
        states.push(LaneState {
            y: record.coordinate.y,
            z: record.coordinate.z,
            predecessor_j: None,
            current: point,
        });
    }
    states.sort_by_key(|state| (state.z, state.y));
    Ok(states)
}

fn verify_row_record(
    context: &Context,
    state: &mut LaneState,
    x: u64,
    record: NodeRecord,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<(), String> {
    let coordinate = Coordinate {
        x,
        y: state.y,
        z: state.z,
    };
    let parent = Coordinate {
        x: x - 1,
        y: state.y,
        z: state.z,
    };
    if record.record != "node"
        || record.coordinate != coordinate
        || record.index != coordinate_index(coordinate)?
        || record.edge.as_ref().map(|edge| edge.parent) != Some(parent)
    {
        return Err(format!(
            "row node ({x},{},{}) metadata differs",
            state.y, state.z
        ));
    }
    let point = verify_edge(
        context,
        &context.row_phi,
        state.current,
        state.predecessor_j,
        &record,
        audit_points,
        audit_seed_x,
    )?;
    state.predecessor_j = Some(state.current.j);
    state.current = point;
    Ok(())
}

pub fn verify_segment_jsonl<R: BufRead>(
    reader: R,
    checkpoint: &Checkpoint,
    previous: Option<(&Checkpoint, &str)>,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<VerificationReceipt, String> {
    if audit_points == 0 {
        return Err("audit points must be positive".into());
    }
    let started = Instant::now();
    let context = build_context()?;
    let mut input = LineReader::new(reader);
    let header_record: SegmentHeader = input.read(true)?;
    if header_record.schema != SEGMENT_SCHEMA || header_record.record != "header" {
        return Err("segment header schema differs".into());
    }
    let length = header_record
        .x_end
        .checked_sub(header_record.x_start)
        .ok_or("segment header has x_end before x_start")?;
    let config = SegmentConfig {
        x_start: header_record.x_start,
        length,
        batch_steps: header_record.batch_steps,
        audit_points: header_record.construction_audit_points,
        audit_seed_x: header_record.construction_audit_seed_x,
        source_commit: header_record.source_commit.clone(),
        class_number_receipt_sha256: header_record.class_number_receipt_sha256.clone(),
        previous_checkpoint_sha256: header_record.previous_checkpoint_sha256.clone(),
    };
    config.validate()?;
    if header_record != header(&config, &context) {
        return Err("segment header differs from the frozen p61 header".into());
    }

    let mut states = if config.x_start == 0 {
        if previous.is_some() {
            return Err("first segment cannot have a predecessor checkpoint".into());
        }
        verify_seed_records(&mut input, &context, audit_points, audit_seed_x)?
    } else {
        let (prior, digest) = previous.ok_or("continuation verification lacks a predecessor")?;
        validate_checkpoint(&config, prior, digest)?;
        states_from_checkpoint(&context, prior)?
    };
    let mut parent_edges = if config.x_start == 0 {
        (LANES - 1) as u64
    } else {
        0
    };
    let first_new_x = if config.x_start == 0 {
        1
    } else {
        config.x_start
    };
    let mut block_start = first_new_x;
    while block_start < config.x_end() {
        let block_end = config.x_end().min(block_start + config.batch_steps as u64);
        let mut rows = Vec::with_capacity(states.len());
        for _ in 0..states.len() {
            let mut row = Vec::with_capacity((block_end - block_start) as usize);
            for _ in block_start..block_end {
                let record: NodeRecord = input.read(true)?;
                row.push(record);
            }
            rows.push(row);
        }
        let verified: Vec<Result<u64, String>> = states
            .par_iter_mut()
            .zip(rows.into_par_iter())
            .map(|(state, row)| {
                for (x, record) in (block_start..block_end).zip(row) {
                    verify_row_record(&context, state, x, record, audit_points, audit_seed_x)?;
                }
                Ok(block_end - block_start)
            })
            .collect();
        for count in verified {
            parent_edges = parent_edges
                .checked_add(count?)
                .ok_or("verification parent-edge count overflow")?;
        }
        block_start = block_end;
    }

    let mut closure_edges = 0u64;
    if config.x_end() == ROW_ORDER {
        let seeds = lane_seeds(&context, 1, 0)?;
        for (state, seed) in states.iter().zip(seeds.iter()) {
            let closure: ClosureRecord = input.read(true)?;
            let expected = closure_record(&context, state, seed)?;
            if closure != expected {
                return Err("row closure certificate differs from independent replay".into());
            }
            closure_edges += 1;
        }
    }
    let summary: SegmentSummary = input.read(false)?;
    let expected_summary = SegmentSummary {
        record: "summary".into(),
        schema: SUMMARY_SCHEMA.into(),
        x_start: config.x_start,
        x_end: config.x_end(),
        node_records: config.length * LANES as u64,
        parent_edges,
        closure_edges,
        order_consistent_unproved: config.length * LANES as u64,
        record_chain_sha256: hex::encode(input.chain),
    };
    if summary != expected_summary {
        return Err("segment summary or record-chain digest differs".into());
    }
    let chain = input.chain;
    let bytes = input.finish()?;
    let expected_checkpoint = checkpoint_from_states(&config, &states, &context, chain)?;
    if checkpoint.schema != expected_checkpoint.schema
        || checkpoint.source_commit != expected_checkpoint.source_commit
        || checkpoint.class_number_receipt_sha256 != expected_checkpoint.class_number_receipt_sha256
        || checkpoint.previous_checkpoint_sha256 != expected_checkpoint.previous_checkpoint_sha256
        || checkpoint.x_start != expected_checkpoint.x_start
        || checkpoint.next_x != expected_checkpoint.next_x
        || checkpoint.record_chain_sha256 != expected_checkpoint.record_chain_sha256
        || checkpoint.lanes != expected_checkpoint.lanes
        || !is_lower_hex(&checkpoint.segment_sha256, 64)
    {
        return Err("checkpoint differs from replayed terminal lane states".into());
    }
    Ok(VerificationReceipt {
        schema: RECEIPT_SCHEMA.into(),
        status: "pass".into(),
        source_commit: config.source_commit,
        audit_points,
        audit_seed_x,
        elapsed_ms: started.elapsed().as_millis(),
        uncompressed_bytes: bytes,
        summary,
    })
}

pub fn plan_record() -> serde_json::Value {
    serde_json::json!({
        "schema": "ec-hardness/p61-class-plan/v1",
        "root": "p61-i0025",
        "p": P,
        "trace": TRACE,
        "fundamental_discriminant": D_K,
        "class_group_factors": CLASS_GROUP_FACTORS,
        "expected_class_vertices": EXPECTED_CLASS_VERTICES,
        "lattice": {
            "row_degree": ROW_DEGREE,
            "row_class_vector": ROW_CLASS_VECTOR,
            "row_order": ROW_ORDER,
            "coset_degree": COSET_DEGREE,
            "coset_class_vector": COSET_CLASS_VECTOR,
            "coset_rows": COSET_ROWS,
            "top_degree": TOP_DEGREE,
            "top_class_vector": TOP_CLASS_VECTOR,
            "top_rows": TOP_ROWS,
            "lanes": LANES,
        },
        "full_adjacency_degrees_after_vertex_census": [5, 13, 37, 59],
        "memory_contract": "18 lane states plus 18*batch_steps node records",
        "completion_gate": "certified class number, 358094304 externally unique j-invariants, every segment replayed, and 18 row closures replayed",
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lattice_dimensions_match_the_selected_class() {
        assert_eq!(ROW_ORDER, CLASS_GROUP_FACTORS[0] / 9);
        assert_eq!(
            ROW_ORDER * u64::from(COSET_ROWS) * u64::from(TOP_ROWS),
            EXPECTED_CLASS_VERTICES
        );
        assert_eq!(LANES, 18);
    }

    #[test]
    fn coordinate_index_is_dense_and_rejects_out_of_range_values() {
        assert_eq!(coordinate_index(Coordinate { x: 0, y: 0, z: 0 }), Ok(0));
        assert_eq!(
            coordinate_index(Coordinate {
                x: ROW_ORDER - 1,
                y: COSET_ROWS - 1,
                z: TOP_ROWS - 1
            }),
            Ok(EXPECTED_CLASS_VERTICES - 1)
        );
        assert!(coordinate_index(Coordinate {
            x: ROW_ORDER,
            y: 0,
            z: 0
        })
        .is_err());
    }

    #[test]
    fn segment_config_requires_a_hash_bound_predecessor() {
        let base = SegmentConfig {
            x_start: 0,
            length: 10,
            batch_steps: 4,
            audit_points: 1,
            audit_seed_x: 0,
            source_commit: "a".repeat(40),
            class_number_receipt_sha256: "b".repeat(64),
            previous_checkpoint_sha256: None,
        };
        assert!(base.validate().is_ok());
        let mut continuation = base.clone();
        continuation.x_start = 10;
        assert!(continuation.validate().is_err());
        continuation.previous_checkpoint_sha256 = Some("c".repeat(64));
        assert!(continuation.validate().is_ok());
    }

    #[test]
    fn record_chain_binds_order_and_bytes() {
        let mut first = [0u8; 32];
        update_chain(&mut first, b"one\n");
        update_chain(&mut first, b"two\n");
        let mut second = [0u8; 32];
        update_chain(&mut second, b"two\n");
        update_chain(&mut second, b"one\n");
        assert_ne!(first, second);
    }
}
