//! Disjoint global-row continuations of the frozen P-256 product grid.
//!
//! A strip fixes `x = 0..width` and a half-open global `y` interval. The
//! degree-13 spine prefix before the interval is recomputed, not counted. For
//! a nonzero row start, one chain-bound boundary curve supplies the source of
//! the first emitted degree-13 edge. Exact cross-certificate uniqueness is a
//! separate streaming audit over canonical 256-bit j-invariants.

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::time::Instant;

use flate2::read::GzDecoder;
use serde_json::Value;

use super::*;
use crate::hash::sha256::Sha256;

use super::super::million_store::{digest_file, FileDigest};

pub const STRIP_SCHEMA: &str = "p256.isogeny-grid-strip-certificate/v1";
pub const STRIP_SUMMARY_SCHEMA: &str = "p256.isogeny-grid-strip-summary/v1";
pub const STRIP_RECEIPT_SCHEMA: &str = "p256.isogeny-grid-strip-receipt/v1";
pub const J_UNION_SCHEMA: &str = "p256.isogeny-grid-j-union/v1";

#[derive(Clone, Debug)]
pub struct StripConfig {
    pub width: u32,
    pub y_start: u32,
    pub height: u32,
    pub batch_rows: usize,
    pub audit_points: usize,
    pub audit_seed_x: u64,
    pub source_commit: String,
}

impl StripConfig {
    fn validate(&self) -> Result<(), String> {
        if self.width < 2 {
            return Err("strip width must be at least two".into());
        }
        if self.height == 0 {
            return Err("strip height must be positive".into());
        }
        self.y_start
            .checked_add(self.height)
            .ok_or("strip row interval exceeds u32")?;
        self.target_curves()
            .checked_add(2)
            .ok_or("strip record count overflow")?;
        if self.batch_rows == 0 {
            return Err("row batch must be positive".into());
        }
        if self.audit_points == 0 {
            return Err("audit point count must be positive".into());
        }
        if self.audit_points != 1 || self.audit_seed_x != 0 {
            return Err("frozen strip construction uses one order-audit point from x = 0".into());
        }
        if self.source_commit.len() != 40
            || !self
                .source_commit
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit())
        {
            return Err("source commit must be a 40-digit hexadecimal Git object id".into());
        }
        Ok(())
    }

    fn y_end(&self) -> u32 {
        self.y_start + self.height
    }

    fn target_curves(&self) -> u64 {
        u64::from(self.width) * u64::from(self.height)
    }

    fn context_records(&self) -> u64 {
        u64::from(self.y_start != 0)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StripHeaderRecord {
    pub record: String,
    pub schema: String,
    pub source_commit: String,
    pub root_name: String,
    pub root_icv1_slug: String,
    pub field_p: String,
    pub group_order: String,
    pub trace: String,
    pub width: u32,
    pub y_start: u32,
    pub height: u32,
    pub target_curves: u64,
    pub context_records: u64,
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
struct BoundaryRecord {
    record: String,
    coordinate: Coordinate,
    predecessor_j: Option<String>,
    curve: CurveRecord,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UnionInputReceipt {
    pub input: String,
    pub schema: String,
    pub source_commit: String,
    pub target_curves: u64,
    pub curve_records: u64,
    pub bound_records: u64,
    pub record_chain_sha256: String,
    pub stored_sha256: String,
    pub stored_bytes: u64,
    pub content_sha256: String,
    pub content_bytes: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JUnionReceipt {
    pub schema: String,
    pub status: String,
    pub elapsed_ms: u128,
    pub resources: ResourceUsage,
    pub inputs: Vec<UnionInputReceipt>,
    pub input_curves: u64,
    pub unique_j_invariants: u64,
}

#[derive(Clone, Copy)]
struct SeenJ {
    input: usize,
    index: u64,
    coordinate: Coordinate,
}

fn coordinate_index(coordinate: Coordinate) -> u64 {
    (u64::from(coordinate.y) << 32) | u64::from(coordinate.x)
}

fn parent_coordinate(coordinate: Coordinate) -> Result<Coordinate, String> {
    if coordinate.x != 0 {
        Ok(Coordinate {
            x: coordinate.x - 1,
            y: coordinate.y,
        })
    } else if coordinate.y != 0 {
        Ok(Coordinate {
            x: 0,
            y: coordinate.y - 1,
        })
    } else {
        Err("the registered root has no parent".into())
    }
}

fn make_strip_header(
    config: &StripConfig,
    start: &StartCurve,
    root: &CurveRecord,
    spine: &ModPoly,
    row: &ModPoly,
) -> StripHeaderRecord {
    let trace: BigInt =
        BigInt::from(start.p.clone()) + BigInt::from(1u8) - BigInt::from(start.order.clone());
    StripHeaderRecord {
        record: "header".into(),
        schema: STRIP_SCHEMA.into(),
        source_commit: config.source_commit.clone(),
        root_name: start.name.clone(),
        root_icv1_slug: root.icv1_slug.clone(),
        field_p: format!("{:064x}", start.p),
        group_order: format!("{:064x}", start.order),
        trace: trace.to_string(),
        width: config.width,
        y_start: config.y_start,
        height: config.height,
        target_curves: config.target_curves(),
        context_records: config.context_records(),
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

fn validate_strip_header(
    header: &StripHeaderRecord,
    start: &StartCurve,
    root: &CurveRecord,
    spine: &ModPoly,
    row: &ModPoly,
) -> Result<StripConfig, String> {
    let config = StripConfig {
        width: header.width,
        y_start: header.y_start,
        height: header.height,
        batch_rows: 1,
        audit_points: header.construction_audit_points,
        audit_seed_x: header.construction_audit_seed_x,
        source_commit: header.source_commit.clone(),
    };
    config.validate()?;
    if header != &make_strip_header(&config, start, root, spine, row) {
        return Err("certificate header differs from the frozen P-256 strip header".into());
    }
    Ok(config)
}

fn root_point(field: &Field, start: &StartCurve) -> Result<GridPoint, String> {
    let model = Model {
        a: field.from_big(&start.a),
        b: field.from_big(&start.b),
    };
    let j = model.j(field).ok_or("registered P-256 is singular")?;
    Ok(GridPoint { model, j })
}

fn spine_boundary(
    field: &Field,
    phi: &ModPoly,
    root: GridPoint,
    y_start: u32,
) -> Result<(GridPoint, Option<GridPoint>), String> {
    if y_start == 0 {
        return Ok((root, None));
    }
    let mut source = root;
    let mut predecessor = None;
    for _y in 1..y_start {
        let (model, j, _, _) = make_step(
            field,
            phi,
            source,
            predecessor.map(|point: GridPoint| point.j),
        )?;
        predecessor = Some(source);
        source = GridPoint { model, j };
    }
    Ok((source, predecessor))
}

#[allow(clippy::too_many_arguments)]
fn make_strip_node(
    field: &Field,
    start: &StartCurve,
    coordinate: Coordinate,
    phi: &ModPoly,
    source: GridPoint,
    predecessor: Option<Fe>,
    order_is_prime: bool,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<GeneratedNode, String> {
    let (model, j, kernel, u2) = make_step(field, phi, source, predecessor)?;
    let curve = describe_curve(
        field,
        start,
        &model,
        false,
        order_is_prime,
        audit_points,
        audit_seed_x,
    )?;
    let parent = parent_coordinate(coordinate)?;
    Ok(GeneratedNode {
        record: NodeRecord {
            record: "node".into(),
            index: coordinate_index(coordinate),
            coordinate,
            parent: coordinate_index(parent),
            degree: phi.ell,
            predecessor_j: predecessor.map(|value| hex_fe(field, &value)),
            kernel_polynomial: kernel.iter().map(|value| hex_fe(field, value)).collect(),
            isomorphism_u2: hex_fe(field, &u2),
            curve,
        },
        model,
        j,
    })
}

fn strip_summary(
    accumulator: &Accumulator,
    chain: [u8; 32],
    config: &StripConfig,
) -> Result<GridSummary, String> {
    let mut summary = accumulator.summary(chain)?;
    summary.schema = STRIP_SUMMARY_SCHEMA.into();
    summary.bound_records = config
        .target_curves()
        .checked_add(1 + config.context_records())
        .ok_or("strip bound-record count overflow")?;
    let certified = config
        .target_curves()
        .saturating_sub(u64::from(config.y_start == 0));
    summary.parent_edges = certified;
    summary.kernel_certificates = certified;
    Ok(summary)
}

pub fn generate_strip_jsonl<W: Write>(
    writer: &mut W,
    config: &StripConfig,
) -> Result<GenerationReceipt, String> {
    config.validate()?;
    let started = Instant::now();
    let (start, field, spine_phi, row_phi) = build_context()?;
    let order_is_prime = is_probable_prime(&start.order);
    if !order_is_prime {
        return Err("the registered P-256 group order is not prime".into());
    }
    let root = root_point(&field, &start)?;
    let root_curve = describe_curve(
        &field,
        &start,
        &root.model,
        true,
        order_is_prime,
        config.audit_points,
        config.audit_seed_x,
    )?;
    let header = make_strip_header(config, &start, &root_curve, &spine_phi, &row_phi);
    let mut output = LineWriter::new(writer);
    output.write(&header, true)?;
    let mut accumulator = Accumulator::new();
    let mut spine = Vec::with_capacity(config.height as usize);

    if config.y_start == 0 {
        let root_record = RootRecord {
            record: "root".into(),
            index: 0,
            coordinate: Coordinate { x: 0, y: 0 },
            curve: root_curve,
        };
        accumulator.observe(
            root_record.index,
            root_record.coordinate,
            &root.model,
            root.j,
            &root_record.curve,
        )?;
        output.write(&root_record, true)?;
        spine.push(root);
        let mut source = root;
        let mut predecessor = None;
        for y in 1..config.y_end() {
            let generated = make_strip_node(
                &field,
                &start,
                Coordinate { x: 0, y },
                &spine_phi,
                source,
                predecessor.map(|point: GridPoint| point.j),
                order_is_prime,
                config.audit_points,
                config.audit_seed_x,
            )?;
            predecessor = Some(source);
            source = GridPoint {
                model: generated.model,
                j: generated.j,
            };
            accumulator.observe(
                generated.record.index,
                generated.record.coordinate,
                &generated.model,
                generated.j,
                &generated.record.curve,
            )?;
            output.write(&generated.record, true)?;
            spine.push(source);
        }
    } else {
        let (mut source, mut predecessor) =
            spine_boundary(&field, &spine_phi, root, config.y_start)?;
        let boundary_y = config.y_start - 1;
        let boundary_curve = describe_curve(
            &field,
            &start,
            &source.model,
            boundary_y == 0,
            order_is_prime,
            config.audit_points,
            config.audit_seed_x,
        )?;
        output.write(
            &BoundaryRecord {
                record: "boundary".into(),
                coordinate: Coordinate {
                    x: 0,
                    y: boundary_y,
                },
                predecessor_j: predecessor.map(|point| hex_fe(&field, &point.j)),
                curve: boundary_curve,
            },
            true,
        )?;
        for y in config.y_start..config.y_end() {
            let generated = make_strip_node(
                &field,
                &start,
                Coordinate { x: 0, y },
                &spine_phi,
                source,
                predecessor.map(|point| point.j),
                order_is_prime,
                config.audit_points,
                config.audit_seed_x,
            )?;
            predecessor = Some(source);
            source = GridPoint {
                model: generated.model,
                j: generated.j,
            };
            accumulator.observe(
                generated.record.index,
                generated.record.coordinate,
                &generated.model,
                generated.j,
                &generated.record.curve,
            )?;
            output.write(&generated.record, true)?;
            spine.push(source);
        }
    }
    eprintln!(
        "p256-isogeny-million: generated strip spine rows {}..{}",
        config.y_start,
        config.y_end()
    );

    let batch_rows = config.batch_rows.min(config.height as usize).max(1);
    for batch_start in (0..config.height as usize).step_by(batch_rows) {
        let batch_end = (batch_start + batch_rows).min(config.height as usize);
        let rows: Vec<Result<Vec<GeneratedNode>, String>> = (batch_start..batch_end)
            .into_par_iter()
            .map(|offset| {
                let y = config.y_start + offset as u32;
                let mut source = spine[offset];
                let mut predecessor = None;
                let mut row = Vec::with_capacity((config.width - 1) as usize);
                for x in 1..config.width {
                    let generated = make_strip_node(
                        &field,
                        &start,
                        Coordinate { x, y },
                        &row_phi,
                        source,
                        predecessor.map(|point: GridPoint| point.j),
                        order_is_prime,
                        config.audit_points,
                        config.audit_seed_x,
                    )?;
                    predecessor = Some(source);
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
            "p256-isogeny-million: strip rows {}/{}; {} unique curves",
            batch_end, config.height, accumulator.seen
        );
    }
    if accumulator.seen != config.target_curves() {
        return Err(format!(
            "generated {} strip curves, expected {}",
            accumulator.seen,
            config.target_curves()
        ));
    }
    let summary = strip_summary(&accumulator, output.chain, config)?;
    output.write(&summary, false)?;
    output
        .writer
        .flush()
        .map_err(|error| format!("flush strip certificate: {error}"))?;
    Ok(GenerationReceipt {
        schema: STRIP_RECEIPT_SCHEMA.into(),
        status: "pass".into(),
        source_commit: config.source_commit.clone(),
        elapsed_ms: started.elapsed().as_millis(),
        resources: resource_usage(),
        uncompressed_bytes: output.bytes,
        summary,
    })
}

#[allow(clippy::too_many_arguments)]
fn verify_strip_node(
    field: &Field,
    start: &StartCurve,
    coordinate: Coordinate,
    record: NodeRecord,
    phi: &ModPoly,
    source: GridPoint,
    predecessor: Option<Fe>,
    order_is_prime: bool,
    audit_points: usize,
    audit_seed_x: u64,
) -> Result<GeneratedNode, String> {
    let parent = parent_coordinate(coordinate)?;
    if record.record != "node"
        || record.coordinate != coordinate
        || record.index != coordinate_index(coordinate)
        || record.parent != coordinate_index(parent)
        || record.degree != phi.ell
        || record.predecessor_j != predecessor.map(|value| hex_fe(field, &value))
    {
        return Err(format!(
            "strip node ({},{}): index, coordinate, parent, degree, or predecessor differs",
            coordinate.x, coordinate.y
        ));
    }
    let selected_j = selected_target(field, phi, &source.model, predecessor)?;
    let (model, j) = parse_model(
        field,
        &record.curve,
        &format!("strip node ({},{})", coordinate.x, coordinate.y),
    )?;
    if j != selected_j {
        return Err(format!(
            "strip node ({},{}): target j is not the frozen choice",
            coordinate.x, coordinate.y
        ));
    }
    let kernel: Vec<Fe> = record
        .kernel_polynomial
        .iter()
        .enumerate()
        .map(|(position, value)| {
            parse_fe(
                field,
                value,
                &format!(
                    "strip node ({},{}) kernel coefficient {position}",
                    coordinate.x, coordinate.y
                ),
            )
        })
        .collect::<Result<_, _>>()?;
    let velu = kernel::verify_kernel(field, &source.model, phi.ell, &kernel).map_err(|error| {
        format!(
            "strip node ({},{}) kernel: {error:?}",
            coordinate.x, coordinate.y
        )
    })?;
    let linked_u2 = kernel::link_to_target(field, &velu, &model).map_err(|error| {
        format!(
            "strip node ({},{}) target link: {error:?}",
            coordinate.x, coordinate.y
        )
    })?;
    let recorded_u2 = parse_fe(
        field,
        &record.isomorphism_u2,
        &format!("strip node ({},{}) isomorphism", coordinate.x, coordinate.y),
    )?;
    let (canonical, canonical_u2) = curve::canonical_model(field, &velu).ok_or_else(|| {
        format!(
            "strip node ({},{}) exceptional target",
            coordinate.x, coordinate.y
        )
    })?;
    if linked_u2 != recorded_u2 || canonical_u2 != recorded_u2 || canonical != model {
        return Err(format!(
            "strip node ({},{}): target is not the canonical Velu codomain",
            coordinate.x, coordinate.y
        ));
    }
    verify_curve_record(
        field,
        start,
        &record.curve,
        &model,
        false,
        order_is_prime,
        audit_points,
        audit_seed_x,
        &format!("strip node ({},{})", coordinate.x, coordinate.y),
    )?;
    Ok(GeneratedNode { record, model, j })
}

pub fn verify_strip_jsonl<R: BufRead>(
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
    let root = root_point(&field, &start)?;
    let mut input = LineReader::new(reader);
    let header: StripHeaderRecord = input.read(true)?;
    let construction_root = describe_curve(
        &field,
        &start,
        &root.model,
        true,
        order_is_prime,
        header.construction_audit_points,
        header.construction_audit_seed_x,
    )?;
    let config = validate_strip_header(&header, &start, &construction_root, &spine_phi, &row_phi)?;
    let mut accumulator = Accumulator::new();
    let mut spine = Vec::with_capacity(config.height as usize);

    if config.y_start == 0 {
        let root_record: RootRecord = input.read(true)?;
        if root_record.record != "root"
            || root_record.index != 0
            || root_record.coordinate != (Coordinate { x: 0, y: 0 })
        {
            return Err("the strip root is not P-256 at global coordinate (0,0)".into());
        }
        let (model, j) = parse_model(&field, &root_record.curve, "strip root")?;
        if model != root.model || j != root.j {
            return Err("the strip root differs from registered P-256".into());
        }
        verify_curve_record(
            &field,
            &start,
            &root_record.curve,
            &model,
            true,
            order_is_prime,
            audit_points,
            audit_seed_x,
            "strip root",
        )?;
        accumulator.observe(0, root_record.coordinate, &model, j, &root_record.curve)?;
        spine.push(root);
        let mut source = root;
        let mut predecessor = None;
        for y in 1..config.y_end() {
            let record: NodeRecord = input.read(true)?;
            let checked = verify_strip_node(
                &field,
                &start,
                Coordinate { x: 0, y },
                record,
                &spine_phi,
                source,
                predecessor.map(|point: GridPoint| point.j),
                order_is_prime,
                audit_points,
                audit_seed_x,
            )?;
            predecessor = Some(source);
            source = GridPoint {
                model: checked.model,
                j: checked.j,
            };
            accumulator.observe(
                checked.record.index,
                checked.record.coordinate,
                &checked.model,
                checked.j,
                &checked.record.curve,
            )?;
            spine.push(source);
        }
    } else {
        let (mut source, mut predecessor) =
            spine_boundary(&field, &spine_phi, root, config.y_start)?;
        let boundary: BoundaryRecord = input.read(true)?;
        let boundary_coordinate = Coordinate {
            x: 0,
            y: config.y_start - 1,
        };
        if boundary.record != "boundary"
            || boundary.coordinate != boundary_coordinate
            || boundary.predecessor_j != predecessor.map(|point| hex_fe(&field, &point.j))
        {
            return Err("strip boundary coordinate or predecessor differs".into());
        }
        let (boundary_model, boundary_j) = parse_model(&field, &boundary.curve, "strip boundary")?;
        if boundary_model != source.model || boundary_j != source.j {
            return Err("strip boundary curve differs from the replayed spine prefix".into());
        }
        verify_curve_record(
            &field,
            &start,
            &boundary.curve,
            &boundary_model,
            boundary_coordinate.y == 0,
            order_is_prime,
            audit_points,
            audit_seed_x,
            "strip boundary",
        )?;
        for y in config.y_start..config.y_end() {
            let record: NodeRecord = input.read(true)?;
            let checked = verify_strip_node(
                &field,
                &start,
                Coordinate { x: 0, y },
                record,
                &spine_phi,
                source,
                predecessor.map(|point| point.j),
                order_is_prime,
                audit_points,
                audit_seed_x,
            )?;
            predecessor = Some(source);
            source = GridPoint {
                model: checked.model,
                j: checked.j,
            };
            accumulator.observe(
                checked.record.index,
                checked.record.coordinate,
                &checked.model,
                checked.j,
                &checked.record.curve,
            )?;
            spine.push(source);
        }
    }
    eprintln!(
        "p256-isogeny-million: replayed strip spine rows {}..{}",
        config.y_start,
        config.y_end()
    );

    let batch_rows = batch_rows.min(config.height as usize).max(1);
    for batch_start in (0..config.height as usize).step_by(batch_rows) {
        let batch_end = (batch_start + batch_rows).min(config.height as usize);
        let mut encoded_rows = Vec::with_capacity(batch_end - batch_start);
        for offset in batch_start..batch_end {
            let y = config.y_start + offset as u32;
            let mut row = Vec::with_capacity((config.width - 1) as usize);
            for x in 1..config.width {
                let record: NodeRecord = input.read(true)?;
                if record.coordinate != (Coordinate { x, y }) {
                    return Err(format!(
                        "strip record {}: expected coordinate ({x},{y})",
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
            .map(|(batch_offset, records)| {
                let offset = batch_start + batch_offset;
                let y = config.y_start + offset as u32;
                let mut source = spine[offset];
                let mut predecessor = None;
                let mut checked = Vec::with_capacity(records.len());
                for (record_offset, record) in records.into_iter().enumerate() {
                    let coordinate = Coordinate {
                        x: record_offset as u32 + 1,
                        y,
                    };
                    let node = verify_strip_node(
                        &field,
                        &start,
                        coordinate,
                        record,
                        &row_phi,
                        source,
                        predecessor.map(|point: GridPoint| point.j),
                        order_is_prime,
                        audit_points,
                        audit_seed_x,
                    )?;
                    predecessor = Some(source);
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
            "p256-isogeny-million: replayed strip rows {}/{}; {} unique curves",
            batch_end, config.height, accumulator.seen
        );
    }
    let recorded_summary: GridSummary = input.read(false)?;
    let computed_summary = strip_summary(&accumulator, input.chain, &config)?;
    if recorded_summary != computed_summary {
        return Err("strip summary or record-chain digest differs from replay".into());
    }
    if computed_summary.unique_curves != config.target_curves() {
        return Err(format!(
            "strip has {} curves, expected {}",
            computed_summary.unique_curves,
            config.target_curves()
        ));
    }
    let bytes = input.finish()?;
    Ok(VerificationReceipt {
        schema: STRIP_RECEIPT_SCHEMA.into(),
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

fn parse_fixed_hex(value: &str, what: &str) -> Result<[u8; 32], String> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(format!(
            "{what} must be exactly 64 lowercase hexadecimal digits"
        ));
    }
    let mut out = [0u8; 32];
    hex::decode_to_slice(value, &mut out).map_err(|error| format!("{what}: {error}"))?;
    Ok(out)
}

fn json_u64(value: &Value, field: &str, input: &str, line: u64) -> Result<u64, String> {
    value[field]
        .as_u64()
        .ok_or_else(|| format!("{input} line {line}: missing integer {field}"))
}

fn audit_j_readers(
    mut inputs: Vec<(String, Box<dyn BufRead>, FileDigest)>,
) -> Result<JUnionReceipt, String> {
    if inputs.is_empty() {
        return Err("j-union audit requires at least one certificate".into());
    }
    let started = Instant::now();
    let names: Vec<String> = inputs.iter().map(|(name, _, _)| name.clone()).collect();
    let mut seen = HashMap::<[u8; 32], SeenJ>::new();
    let mut receipts = Vec::with_capacity(inputs.len());
    let mut input_curves = 0u64;

    for (input_index, (name, reader, stored)) in inputs.iter_mut().enumerate() {
        let mut chain = [0u8; 32];
        let mut content_hash = Sha256::new();
        let mut content_bytes = 0u64;
        let mut line_number = 0u64;
        let mut bound_records = 0u64;
        let mut curve_records = 0u64;
        let mut header_schema = None;
        let mut source_commit = None;
        let mut target_curves = None;
        let mut field_p = None;
        let mut recorded_summary = None;

        loop {
            let mut line = Vec::new();
            let bytes = reader
                .read_until(b'\n', &mut line)
                .map_err(|error| format!("read {name}: {error}"))?;
            if bytes == 0 {
                break;
            }
            line_number += 1;
            if line.len() > 16 * 1024 * 1024 {
                return Err(format!("{name} line {line_number}: record exceeds 16 MiB"));
            }
            if line.last() != Some(&b'\n') || line.ends_with(b"\r\n") {
                return Err(format!("{name} line {line_number}: records require one LF"));
            }
            content_hash.update(&line);
            content_bytes = content_bytes
                .checked_add(line.len() as u64)
                .ok_or("j-union content byte count overflow")?;
            let value: Value = serde_json::from_slice(&line[..line.len() - 1])
                .map_err(|error| format!("{name} line {line_number}: {error}"))?;
            let record = value["record"]
                .as_str()
                .ok_or_else(|| format!("{name} line {line_number}: missing record type"))?;
            if recorded_summary.is_some() {
                return Err(format!("{name}: records follow the summary"));
            }
            if record == "summary" {
                let summary: GridSummary = serde_json::from_value(value)
                    .map_err(|error| format!("{name} line {line_number}: {error}"))?;
                recorded_summary = Some(summary);
                continue;
            }
            update_chain(&mut chain, &line);
            bound_records += 1;

            match record {
                "header" => {
                    if line_number != 1 || header_schema.is_some() {
                        return Err(format!("{name}: header is not the first unique record"));
                    }
                    let schema = value["schema"]
                        .as_str()
                        .ok_or_else(|| format!("{name}: header has no schema"))?;
                    if schema != SCHEMA && schema != STRIP_SCHEMA {
                        return Err(format!("{name}: unsupported certificate schema {schema}"));
                    }
                    header_schema = Some(schema.to_owned());
                    source_commit = Some(
                        value["source_commit"]
                            .as_str()
                            .ok_or_else(|| format!("{name}: header has no source commit"))?
                            .to_owned(),
                    );
                    target_curves = Some(json_u64(&value, "target_curves", name, line_number)?);
                    field_p = Some(parse_fixed_hex(
                        value["field_p"]
                            .as_str()
                            .ok_or_else(|| format!("{name}: header has no field modulus"))?,
                        &format!("{name} field modulus"),
                    )?);
                }
                "boundary" => {
                    if header_schema.as_deref() != Some(STRIP_SCHEMA) {
                        return Err(format!(
                            "{name}: boundary record in a non-strip certificate"
                        ));
                    }
                }
                "root" | "node" => {
                    let modulus = field_p.ok_or_else(|| format!("{name}: curve before header"))?;
                    let j_text = value["curve"]["j"]
                        .as_str()
                        .ok_or_else(|| format!("{name} line {line_number}: curve has no j"))?;
                    let j = parse_fixed_hex(j_text, &format!("{name} line {line_number} j"))?;
                    if j >= modulus {
                        return Err(format!("{name} line {line_number}: j is outside F_p"));
                    }
                    let index = json_u64(&value, "index", name, line_number)?;
                    let coordinate = Coordinate {
                        x: json_u64(&value["coordinate"], "x", name, line_number)?
                            .try_into()
                            .map_err(|_| format!("{name} line {line_number}: x exceeds u32"))?,
                        y: json_u64(&value["coordinate"], "y", name, line_number)?
                            .try_into()
                            .map_err(|_| format!("{name} line {line_number}: y exceeds u32"))?,
                    };
                    if let Some(first) = seen.insert(
                        j,
                        SeenJ {
                            input: input_index,
                            index,
                            coordinate,
                        },
                    ) {
                        return Err(format!(
                            "duplicate j {j_text}: first {} index {} coordinate ({},{}), repeated {name} index {index} coordinate ({},{})",
                            names[first.input],
                            first.index,
                            first.coordinate.x,
                            first.coordinate.y,
                            coordinate.x,
                            coordinate.y
                        ));
                    }
                    curve_records += 1;
                }
                other => return Err(format!("{name}: unsupported record type {other}")),
            }
        }
        let schema = header_schema.ok_or_else(|| format!("{name}: missing header"))?;
        let summary = recorded_summary.ok_or_else(|| format!("{name}: missing summary"))?;
        let expected_summary_schema = if schema == SCHEMA {
            SUMMARY_SCHEMA
        } else {
            STRIP_SUMMARY_SCHEMA
        };
        if summary.schema != expected_summary_schema
            || summary.record_chain_sha256 != hex::encode(chain)
            || summary.bound_records != bound_records
            || summary.unique_curves != curve_records
            || summary.unique_j_invariants != curve_records
            || target_curves != Some(curve_records)
        {
            return Err(format!(
                "{name}: summary, chain, bound count, or curve population differs"
            ));
        }
        input_curves = input_curves
            .checked_add(curve_records)
            .ok_or("j-union curve count overflow")?;
        receipts.push(UnionInputReceipt {
            input: name.clone(),
            schema,
            source_commit: source_commit.ok_or_else(|| format!("{name}: missing source commit"))?,
            target_curves: target_curves.expect("checked above"),
            curve_records,
            bound_records,
            record_chain_sha256: summary.record_chain_sha256,
            stored_sha256: stored.sha256.clone(),
            stored_bytes: stored.bytes,
            content_sha256: hex::encode(content_hash.finalize()),
            content_bytes,
        });
    }
    if seen.len() as u64 != input_curves {
        return Err("j-union internal count differs from the input population".into());
    }
    Ok(JUnionReceipt {
        schema: J_UNION_SCHEMA.into(),
        status: "pass".into(),
        elapsed_ms: started.elapsed().as_millis(),
        resources: resource_usage(),
        inputs: receipts,
        input_curves,
        unique_j_invariants: seen.len() as u64,
    })
}

pub fn audit_j_union_paths(paths: &[PathBuf]) -> Result<JUnionReceipt, String> {
    let mut inputs = Vec::with_capacity(paths.len());
    for path in paths {
        let digest = digest_file(path)?;
        let file = File::open(path).map_err(|error| format!("open {}: {error}", path.display()))?;
        let reader: Box<dyn BufRead> =
            if path.extension().and_then(|value| value.to_str()) == Some("gz") {
                Box::new(BufReader::new(GzDecoder::new(file)))
            } else {
                Box::new(BufReader::new(file))
            };
        inputs.push((path.display().to_string(), reader, digest));
    }
    audit_j_readers(inputs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    use super::super::super::million_store::digest_reader;

    fn certificate(y_start: u32, height: u32) -> Vec<u8> {
        let mut bytes = Vec::new();
        generate_strip_jsonl(
            &mut bytes,
            &StripConfig {
                width: 4,
                y_start,
                height,
                batch_rows: 2,
                audit_points: 1,
                audit_seed_x: 0,
                source_commit: "0".repeat(40),
            },
        )
        .expect("generate strip");
        bytes
    }

    fn audit_inputs(named: Vec<(&str, Vec<u8>)>) -> Result<JUnionReceipt, String> {
        let mut inputs = Vec::new();
        for (name, bytes) in named {
            let digest = digest_reader(Cursor::new(&bytes))?;
            inputs.push((
                name.to_owned(),
                Box::new(Cursor::new(bytes)) as Box<dyn BufRead>,
                digest,
            ));
        }
        audit_j_readers(inputs)
    }

    #[test]
    fn adjacent_strips_replay_and_union_while_overlap_is_rejected() {
        let first = certificate(0, 3);
        let second = certificate(3, 2);
        let first_receipt =
            verify_strip_jsonl(Cursor::new(&first), 2, 7, 2).expect("replay first strip");
        let second_receipt =
            verify_strip_jsonl(Cursor::new(&second), 2, 7, 2).expect("replay second strip");
        assert_eq!(first_receipt.summary.unique_curves, 12);
        assert_eq!(first_receipt.summary.parent_edges, 11);
        assert_eq!(second_receipt.summary.unique_curves, 8);
        assert_eq!(second_receipt.summary.parent_edges, 8);

        let union = audit_inputs(vec![("first", first.clone()), ("second", second)])
            .expect("audit adjacent strips");
        assert_eq!(union.input_curves, 20);
        assert_eq!(union.unique_j_invariants, 20);

        let error = audit_inputs(vec![("first", first.clone()), ("overlap", first)])
            .expect_err("overlap must fail");
        assert!(error.contains("duplicate j"));
    }

    #[test]
    fn strip_replay_rejects_parent_tampering() {
        let bytes = certificate(3, 2);
        let mut lines: Vec<Vec<u8>> = bytes
            .split(|byte| *byte == b'\n')
            .filter(|line| !line.is_empty())
            .map(<[u8]>::to_vec)
            .collect();
        let mut first_node: Value = serde_json::from_slice(&lines[2]).expect("node parses");
        first_node["parent"] = Value::from(0);
        lines[2] = serde_json::to_vec(&first_node).expect("node serializes");
        let mut tampered = Vec::new();
        for mut line in lines {
            tampered.append(&mut line);
            tampered.push(b'\n');
        }
        assert!(verify_strip_jsonl(Cursor::new(tampered), 2, 7, 2).is_err());
    }
}
