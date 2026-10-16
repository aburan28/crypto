//! Complete selected-degree adjacency certificates for the frozen p61 class.
//!
//! Vertex census segments are immutable inputs.  Each adjacency shard binds one
//! independently verified census segment by SHA-256, selects a contiguous
//! record window, and emits every rational degree-5, 13, 37, and 59 isogeny
//! from every selected source.  Degree 5 is ramified and contributes one edge;
//! the split degrees contribute two each.  Targets are recorded by canonical
//! model and j-invariant now and are joined to dense vertex indices after the
//! complete vertex census has produced the global j-index.

use std::io::{BufRead, Write};
use std::time::Instant;

use rayon::prelude::*;
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use super::curve::{self, Model};
use super::field::{Fe, Field};
use super::kernel;
use super::modpoly::{self, ModPoly};
use super::p61_class::{
    build_context, coefficient_bits, coordinate_index, hex_fe, is_lower_hex, parse_curve, parse_fe,
    Context, Coordinate, CurveRecord, NodeRecord, SegmentHeader, EXPECTED_CLASS_VERTICES, LANES,
    ROOT_SEED, SEGMENT_SCHEMA,
};
use super::poly;
use super::walk::{ClassInfo, EllKind};
use crate::hash::sha256::sha256;

pub const ADJACENCY_SCHEMA: &str = "ec-hardness/p61-selected-adjacency/v1";
pub const ADJACENCY_SUMMARY_SCHEMA: &str = "ec-hardness/p61-selected-adjacency-summary/v1";
pub const ADJACENCY_RECEIPT_SCHEMA: &str = "ec-hardness/p61-selected-adjacency-receipt/v1";
pub const RAMIFIED_DEGREE: u64 = 5;
pub const SELECTED_DEGREES: [u64; 4] = [5, 13, 37, 59];
pub const EXPECTED_ROOT_COUNTS: [u8; 4] = [1, 2, 2, 2];
pub const DIRECTED_EDGES_PER_VERTEX: u64 = 7;
pub const KERNEL_COEFFICIENTS_PER_VERTEX: u64 = 115;

#[derive(Clone, Debug)]
pub struct AdjacencyConfig {
    pub source_record_offset: u64,
    pub source_record_count: u64,
    pub batch_vertices: usize,
    pub producer_source_commit: String,
    pub class_number_receipt_sha256: String,
    pub source_segment_sha256: String,
    pub source_verification_receipt_sha256: String,
}

impl AdjacencyConfig {
    pub fn validate(&self) -> Result<(), String> {
        if self.source_record_count == 0 {
            return Err("source record count must be positive".into());
        }
        self.source_record_offset
            .checked_add(self.source_record_count)
            .ok_or("source record window overflows u64")?;
        if self.batch_vertices == 0 {
            return Err("batch vertices must be positive".into());
        }
        if !is_lower_hex(&self.producer_source_commit, 40) {
            return Err("producer source commit must be 40 lowercase hexadecimal digits".into());
        }
        for (value, name) in [
            (&self.class_number_receipt_sha256, "class-number receipt"),
            (&self.source_segment_sha256, "source segment"),
            (
                &self.source_verification_receipt_sha256,
                "source verification receipt",
            ),
        ] {
            if !is_lower_hex(value, 64) {
                return Err(format!(
                    "{name} digest must be 64 lowercase hexadecimal digits"
                ));
            }
        }
        Ok(())
    }

    pub fn source_record_end(&self) -> u64 {
        self.source_record_offset + self.source_record_count
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DegreeChecks {
    pub degree: u64,
    pub rational_roots: u8,
    pub modular_polynomial_checks: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AdjacencyHeader {
    pub record: String,
    pub schema: String,
    pub root: String,
    pub field_p: String,
    pub expected_class_vertices: u64,
    pub selected_degrees: [u64; 4],
    pub expected_rational_roots: [u8; 4],
    pub directed_edges_per_vertex: u64,
    pub producer_source_commit: String,
    pub vertex_source_commit: String,
    pub class_number_receipt_sha256: String,
    pub source_segment_schema: String,
    pub source_segment_sha256: String,
    pub source_verification_receipt_sha256: String,
    pub source_x_start: u64,
    pub source_x_end: u64,
    pub source_segment_node_records: u64,
    pub source_record_offset: u64,
    pub source_record_count: u64,
    pub batch_vertices: usize,
    pub degree_checks: Vec<DegreeChecks>,
    pub canonical_model_rule: String,
    pub target_index_status: String,
    pub certificate_rule: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EdgeGeometry {
    pub source_target_j_xor_popcount: u32,
    pub source_target_j_absolute_delta_bits: u64,
    pub source_target_a_absolute_delta_bits: u64,
    pub source_target_b_absolute_delta_bits: u64,
    pub kernel_coefficient_count: u64,
    pub kernel_nonzero_coefficients: u64,
    pub kernel_signed_bits_min: u64,
    pub kernel_signed_bits_max: u64,
    pub kernel_signed_bits_sum: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AdjacencyEdgeRecord {
    pub record: String,
    pub edge_identity_sha256: String,
    pub source_index: u64,
    pub source_coordinate: Coordinate,
    pub degree: u64,
    pub root_ordinal: u8,
    pub rational_root_count: u8,
    pub direction: String,
    pub volcano_level: u32,
    pub source_curve: CurveRecord,
    pub target_curve: CurveRecord,
    pub velu_codomain_a: String,
    pub velu_codomain_b: String,
    pub kernel_polynomial: Vec<String>,
    pub isomorphism_u2: String,
    pub geometry: EdgeGeometry,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DegreeCount {
    pub degree: u64,
    pub directed_edges: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AdjacencySummary {
    pub record: String,
    pub schema: String,
    pub source_record_offset: u64,
    pub source_record_count: u64,
    pub directed_edges: u64,
    pub degree_counts: Vec<DegreeCount>,
    pub kernel_coefficients: u64,
    pub generation_failures: u64,
    pub record_chain_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AdjacencyGenerationReceipt {
    pub schema: String,
    pub status: String,
    pub producer_source_commit: String,
    pub elapsed_ms: u128,
    pub uncompressed_bytes: u64,
    pub summary: AdjacencySummary,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AdjacencyVerificationReceipt {
    pub schema: String,
    pub status: String,
    pub producer_source_commit: String,
    pub verifier_source_commit: String,
    pub elapsed_ms: u128,
    pub uncompressed_bytes: u64,
    pub summary: AdjacencySummary,
}

struct AdjacencyContext {
    base: Context,
    ramified_phi: ModPoly,
}

impl AdjacencyContext {
    fn build() -> Result<Self, String> {
        let base = build_context()?;
        let class = ClassInfo::compute(&base.start, &[RAMIFIED_DEGREE]);
        let degree = class
            .ells
            .iter()
            .find(|entry| entry.0 == RAMIFIED_DEGREE)
            .ok_or("missing degree-5 class information")?;
        if degree.1 != EllKind::Divides
            || class.depth(RAMIFIED_DEGREE) != 0
            || class.expected_neighbours(RAMIFIED_DEGREE) != Some(1)
        {
            return Err("degree 5 is not the frozen ramified one-neighbour action".into());
        }
        let ramified_phi = modpoly::modular_polynomial(&base.field, RAMIFIED_DEGREE)?;
        Ok(Self { base, ramified_phi })
    }

    fn phis(&self) -> [&ModPoly; 4] {
        [
            &self.ramified_phi,
            &self.base.row_phi,
            &self.base.coset_phi,
            &self.base.top_phi,
        ]
    }
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
            .map_err(|error| format!("write adjacency shard: {error}"))?;
        self.bytes = self
            .bytes
            .checked_add(line.len() as u64)
            .ok_or("adjacency shard byte count overflow")?;
        Ok(())
    }
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
            .map_err(|error| format!("read JSONL: {error}"))?;
        self.line_number += 1;
        if read == 0 {
            return Err(format!(
                "line {}: unexpected end of stream",
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
            .ok_or("JSONL byte count overflow")?;
        serde_json::from_slice(&line).map_err(|error| format!("line {}: {error}", self.line_number))
    }

    fn finish(mut self) -> Result<u64, String> {
        let mut trailing = Vec::new();
        self.reader
            .read_to_end(&mut trailing)
            .map_err(|error| format!("read adjacency tail: {error}"))?;
        if !trailing.is_empty() {
            return Err("adjacency shard has trailing bytes".into());
        }
        Ok(self.bytes)
    }
}

fn source_node_count(header: &SegmentHeader) -> Result<u64, String> {
    let columns = header
        .x_end
        .checked_sub(header.x_start)
        .ok_or("source header has x_end before x_start")?;
    columns
        .checked_mul(LANES as u64)
        .ok_or("source node count overflow".into())
}

fn validate_source_header(header: &SegmentHeader, config: &AdjacencyConfig) -> Result<u64, String> {
    if header.record != "header" || header.schema != SEGMENT_SCHEMA {
        return Err("source segment header schema differs".into());
    }
    if header.class_number_receipt_sha256 != config.class_number_receipt_sha256 {
        return Err("source segment uses a different class-number receipt".into());
    }
    let nodes = source_node_count(header)?;
    if config.source_record_end() > nodes {
        return Err(format!(
            "source record window {}..{} exceeds {nodes} node records",
            config.source_record_offset,
            config.source_record_end()
        ));
    }
    Ok(nodes)
}

fn header(
    context: &AdjacencyContext,
    source: &SegmentHeader,
    source_nodes: u64,
    config: &AdjacencyConfig,
) -> AdjacencyHeader {
    let degree_checks = context
        .phis()
        .into_iter()
        .zip(EXPECTED_ROOT_COUNTS)
        .map(|(phi, rational_roots)| DegreeChecks {
            degree: phi.ell,
            rational_roots,
            modular_polynomial_checks: phi.checks_passed,
        })
        .collect();
    AdjacencyHeader {
        record: "header".into(),
        schema: ADJACENCY_SCHEMA.into(),
        root: "p61-i0025".into(),
        field_p: context.base.start.p.to_string(),
        expected_class_vertices: EXPECTED_CLASS_VERTICES,
        selected_degrees: SELECTED_DEGREES,
        expected_rational_roots: EXPECTED_ROOT_COUNTS,
        directed_edges_per_vertex: DIRECTED_EDGES_PER_VERTEX,
        producer_source_commit: config.producer_source_commit.clone(),
        vertex_source_commit: source.source_commit.clone(),
        class_number_receipt_sha256: config.class_number_receipt_sha256.clone(),
        source_segment_schema: source.schema.clone(),
        source_segment_sha256: config.source_segment_sha256.clone(),
        source_verification_receipt_sha256: config
            .source_verification_receipt_sha256
            .clone(),
        source_x_start: source.x_start,
        source_x_end: source.x_end,
        source_segment_node_records: source_nodes,
        source_record_offset: config.source_record_offset,
        source_record_count: config.source_record_count,
        batch_vertices: config.batch_vertices,
        degree_checks,
        canonical_model_rule: curve::CANONICAL_RULE.into(),
        target_index_status: "pending_complete_census_j_to_index_join".into(),
        certificate_rule: "kernel polynomial is squarefree, divides psi_ell, is closed under the half-group generator, and its Velu codomain is F_p-isomorphic to the canonical target".into(),
    }
}

fn curve_record(field: &Field, model: &Model, order_status: &str) -> Result<CurveRecord, String> {
    let j = model.j(field).ok_or("singular adjacency endpoint")?;
    Ok(CurveRecord {
        a: hex_fe(field, &model.a),
        b: hex_fe(field, &model.b),
        j: hex_fe(field, &j),
        order_status: order_status.into(),
        qr_prefix_64: curve::qr_prefix_64(field, model),
        a_signed_bits: coefficient_bits(field, &model.a),
        b_signed_bits: coefficient_bits(field, &model.b),
        a_minus_3_model: !curve::a_minus_3_models(field, model).is_empty(),
    })
}

fn fe_u64(field: &Field, value: &Fe) -> u64 {
    field
        .to_big(value)
        .to_u64_digits()
        .first()
        .copied()
        .unwrap_or(0)
}

fn absolute_delta_bits(field: &Field, left: &Fe, right: &Fe) -> u64 {
    let left = field.to_big(left);
    let right = field.to_big(right);
    if left >= right {
        (left - right).bits()
    } else {
        (right - left).bits()
    }
}

fn edge_geometry(
    field: &Field,
    source: &Model,
    source_j: &Fe,
    target: &Model,
    target_j: &Fe,
    kernel_polynomial: &[Fe],
) -> EdgeGeometry {
    let signed_bits: Vec<u64> = kernel_polynomial
        .iter()
        .map(|value| coefficient_bits(field, value))
        .collect();
    EdgeGeometry {
        source_target_j_xor_popcount: (fe_u64(field, source_j) ^ fe_u64(field, target_j))
            .count_ones(),
        source_target_j_absolute_delta_bits: absolute_delta_bits(field, source_j, target_j),
        source_target_a_absolute_delta_bits: absolute_delta_bits(field, &source.a, &target.a),
        source_target_b_absolute_delta_bits: absolute_delta_bits(field, &source.b, &target.b),
        kernel_coefficient_count: kernel_polynomial.len() as u64,
        kernel_nonzero_coefficients: kernel_polynomial
            .iter()
            .filter(|value| !field.is_zero(value))
            .count() as u64,
        kernel_signed_bits_min: signed_bits.iter().copied().min().unwrap_or(0),
        kernel_signed_bits_max: signed_bits.iter().copied().max().unwrap_or(0),
        kernel_signed_bits_sum: signed_bits.iter().sum(),
    }
}

fn rational_roots(
    context: &AdjacencyContext,
    phi: &ModPoly,
    model: &Model,
    expected: usize,
) -> Result<Vec<Fe>, String> {
    let j = model
        .j(&context.base.field)
        .ok_or("singular source curve")?;
    let mut roots = poly::roots(
        &context.base.field,
        &phi.at_x(&context.base.field, &j),
        ROOT_SEED ^ phi.ell,
    );
    roots.sort_by_key(|value| context.base.field.to_big(value));
    roots.dedup();
    if roots.len() != expected {
        return Err(format!(
            "Phi_{} has {} distinct rational roots, expected {expected}",
            phi.ell,
            roots.len()
        ));
    }
    Ok(roots)
}

fn edge_identity(
    source_index: u64,
    degree: u64,
    root_ordinal: u8,
    source_j: &str,
    target_j: &str,
    kernel: &[String],
) -> String {
    let mut input = b"ec-hardness/p61-selected-edge/v1\0".to_vec();
    input.extend_from_slice(&source_index.to_le_bytes());
    input.extend_from_slice(&degree.to_le_bytes());
    input.push(root_ordinal);
    input.extend_from_slice(source_j.as_bytes());
    input.push(0);
    input.extend_from_slice(target_j.as_bytes());
    for coefficient in kernel {
        input.push(0);
        input.extend_from_slice(coefficient.as_bytes());
    }
    hex::encode(sha256(&input))
}

fn expand_node(
    context: &AdjacencyContext,
    node: &NodeRecord,
) -> Result<Vec<AdjacencyEdgeRecord>, String> {
    if node.record != "node" || node.index != coordinate_index(node.coordinate)? {
        return Err(format!("source node {} metadata differs", node.index));
    }
    let source = parse_curve(
        &context.base,
        &node.curve,
        &format!("source node {}", node.index),
    )?;
    let recomputed_source =
        curve_record(&context.base.field, &source.model, &node.curve.order_status)?;
    if recomputed_source != node.curve {
        return Err(format!("source node {} detector fields differ", node.index));
    }

    let mut out = Vec::with_capacity(DIRECTED_EDGES_PER_VERTEX as usize);
    for (phi, expected_roots) in context.phis().into_iter().zip(EXPECTED_ROOT_COUNTS) {
        let roots = rational_roots(context, phi, &source.model, expected_roots as usize)?;
        for (root_ordinal, target_j) in roots.into_iter().enumerate() {
            let isogeny =
                kernel::elkies_isogeny(&context.base.field, &source.model, phi, &target_j)
                    .map_err(|error| {
                        format!(
                            "source node {} degree-{} edge construction: {error:?}",
                            node.index, phi.ell
                        )
                    })?;
            let (target, canonical_u2) =
                curve::canonical_model(&context.base.field, &isogeny.codomain)
                    .ok_or("isogeny reached an exceptional target j-invariant")?;
            let linked_u2 = kernel::link_to_target(&context.base.field, &isogeny.codomain, &target)
                .map_err(|error| {
                    format!(
                        "source node {} degree-{} target link: {error:?}",
                        node.index, phi.ell
                    )
                })?;
            if canonical_u2 != linked_u2 || target.j(&context.base.field) != Some(target_j) {
                return Err(format!(
                    "source node {} degree-{} canonical target changed the selected j-invariant",
                    node.index, phi.ell
                ));
            }
            let source_curve = node.curve.clone();
            let target_curve = curve_record(
                &context.base.field,
                &target,
                "transported_by_verified_isogeny",
            )?;
            let kernel_polynomial: Vec<String> = isogeny
                .kernel
                .iter()
                .map(|value| hex_fe(&context.base.field, value))
                .collect();
            let edge_identity_sha256 = edge_identity(
                node.index,
                phi.ell,
                root_ordinal as u8,
                &source_curve.j,
                &target_curve.j,
                &kernel_polynomial,
            );
            out.push(AdjacencyEdgeRecord {
                record: "edge".into(),
                edge_identity_sha256,
                source_index: node.index,
                source_coordinate: node.coordinate,
                degree: phi.ell,
                root_ordinal: root_ordinal as u8,
                rational_root_count: expected_roots,
                direction: "horizontal".into(),
                volcano_level: 0,
                source_curve,
                target_curve,
                velu_codomain_a: hex_fe(&context.base.field, &isogeny.codomain.a),
                velu_codomain_b: hex_fe(&context.base.field, &isogeny.codomain.b),
                geometry: edge_geometry(
                    &context.base.field,
                    &source.model,
                    &source.j,
                    &target,
                    &target_j,
                    &isogeny.kernel,
                ),
                kernel_polynomial,
                isomorphism_u2: hex_fe(&context.base.field, &linked_u2),
            });
        }
    }
    if out.len() != DIRECTED_EDGES_PER_VERTEX as usize {
        return Err(format!(
            "source node {} produced {} edges, expected {}",
            node.index,
            out.len(),
            DIRECTED_EDGES_PER_VERTEX
        ));
    }
    Ok(out)
}

fn verify_node_edges(
    context: &AdjacencyContext,
    node: &NodeRecord,
    records: &[AdjacencyEdgeRecord],
) -> Result<(), String> {
    if node.record != "node" || node.index != coordinate_index(node.coordinate)? {
        return Err(format!("source node {} metadata differs", node.index));
    }
    if records.len() != DIRECTED_EDGES_PER_VERTEX as usize {
        return Err(format!(
            "source node {} has {} recorded edges, expected {}",
            node.index,
            records.len(),
            DIRECTED_EDGES_PER_VERTEX
        ));
    }
    let source = parse_curve(
        &context.base,
        &node.curve,
        &format!("source node {}", node.index),
    )?;
    let recomputed_source =
        curve_record(&context.base.field, &source.model, &node.curve.order_status)?;
    if recomputed_source != node.curve {
        return Err(format!("source node {} detector fields differ", node.index));
    }

    let mut position = 0usize;
    for (phi, expected_roots) in context.phis().into_iter().zip(EXPECTED_ROOT_COUNTS) {
        let roots = rational_roots(context, phi, &source.model, expected_roots as usize)?;
        for (root_ordinal, target_j) in roots.into_iter().enumerate() {
            let recorded = &records[position];
            position += 1;
            if recorded.record != "edge"
                || recorded.source_index != node.index
                || recorded.source_coordinate != node.coordinate
                || recorded.degree != phi.ell
                || recorded.root_ordinal != root_ordinal as u8
                || recorded.rational_root_count != expected_roots
                || recorded.direction != "horizontal"
                || recorded.volcano_level != 0
                || recorded.source_curve != node.curve
            {
                return Err(format!(
                    "source node {} degree-{} edge metadata differs",
                    node.index, phi.ell
                ));
            }
            let target = parse_curve(
                &context.base,
                &recorded.target_curve,
                &format!("source node {} degree-{} target", node.index, phi.ell),
            )?;
            if target.j != target_j {
                return Err(format!(
                    "source node {} degree-{} target is not rational root {}",
                    node.index, phi.ell, root_ordinal
                ));
            }
            let recomputed_target = curve_record(
                &context.base.field,
                &target.model,
                "transported_by_verified_isogeny",
            )?;
            if recorded.target_curve != recomputed_target {
                return Err(format!(
                    "source node {} degree-{} target detector fields differ",
                    node.index, phi.ell
                ));
            }
            let kernel_polynomial: Vec<Fe> = recorded
                .kernel_polynomial
                .iter()
                .enumerate()
                .map(|(index, value)| {
                    parse_fe(
                        &context.base.field,
                        value,
                        &format!(
                            "source node {} degree-{} kernel coefficient {index}",
                            node.index, phi.ell
                        ),
                    )
                })
                .collect::<Result<_, _>>()?;
            let velu = kernel::verify_kernel(
                &context.base.field,
                &source.model,
                phi.ell,
                &kernel_polynomial,
            )
            .map_err(|error| {
                format!(
                    "source node {} degree-{} kernel replay: {error:?}",
                    node.index, phi.ell
                )
            })?;
            if recorded.velu_codomain_a != hex_fe(&context.base.field, &velu.a)
                || recorded.velu_codomain_b != hex_fe(&context.base.field, &velu.b)
            {
                return Err(format!(
                    "source node {} degree-{} Velu codomain differs",
                    node.index, phi.ell
                ));
            }
            let (canonical, canonical_u2) = curve::canonical_model(&context.base.field, &velu)
                .ok_or("kernel replay reached an exceptional target")?;
            let linked_u2 = kernel::link_to_target(&context.base.field, &velu, &target.model)
                .map_err(|error| {
                    format!(
                        "source node {} degree-{} replay target link: {error:?}",
                        node.index, phi.ell
                    )
                })?;
            if canonical != target.model
                || canonical_u2 != linked_u2
                || recorded.isomorphism_u2 != hex_fe(&context.base.field, &linked_u2)
            {
                return Err(format!(
                    "source node {} degree-{} canonical target link differs",
                    node.index, phi.ell
                ));
            }
            let geometry = edge_geometry(
                &context.base.field,
                &source.model,
                &source.j,
                &target.model,
                &target.j,
                &kernel_polynomial,
            );
            let identity = edge_identity(
                node.index,
                phi.ell,
                root_ordinal as u8,
                &node.curve.j,
                &recorded.target_curve.j,
                &recorded.kernel_polynomial,
            );
            if recorded.geometry != geometry || recorded.edge_identity_sha256 != identity {
                return Err(format!(
                    "source node {} degree-{} geometry or identity differs",
                    node.index, phi.ell
                ));
            }
        }
    }
    Ok(())
}

fn expected_degree_counts(source_vertices: u64) -> Result<Vec<DegreeCount>, String> {
    SELECTED_DEGREES
        .into_iter()
        .zip(EXPECTED_ROOT_COUNTS)
        .map(|(degree, roots)| {
            Ok(DegreeCount {
                degree,
                directed_edges: source_vertices
                    .checked_mul(u64::from(roots))
                    .ok_or("degree edge count overflow")?,
            })
        })
        .collect()
}

fn summary(config: &AdjacencyConfig, chain: [u8; 32]) -> Result<AdjacencySummary, String> {
    Ok(AdjacencySummary {
        record: "summary".into(),
        schema: ADJACENCY_SUMMARY_SCHEMA.into(),
        source_record_offset: config.source_record_offset,
        source_record_count: config.source_record_count,
        directed_edges: config
            .source_record_count
            .checked_mul(DIRECTED_EDGES_PER_VERTEX)
            .ok_or("directed edge count overflow")?,
        degree_counts: expected_degree_counts(config.source_record_count)?,
        kernel_coefficients: config
            .source_record_count
            .checked_mul(KERNEL_COEFFICIENTS_PER_VERTEX)
            .ok_or("kernel coefficient count overflow")?,
        generation_failures: 0,
        record_chain_sha256: hex::encode(chain),
    })
}

fn skip_source_records<R: BufRead>(source: &mut LineReader<R>, count: u64) -> Result<(), String> {
    for _ in 0..count {
        let _: NodeRecord = source.read(false)?;
    }
    Ok(())
}

pub fn generate_adjacency_jsonl<R: BufRead, W: Write>(
    source_reader: R,
    writer: &mut W,
    config: &AdjacencyConfig,
) -> Result<AdjacencyGenerationReceipt, String> {
    config.validate()?;
    let started = Instant::now();
    let context = AdjacencyContext::build()?;
    let mut source = LineReader::new(source_reader);
    let source_header: SegmentHeader = source.read(false)?;
    let source_nodes = validate_source_header(&source_header, config)?;
    skip_source_records(&mut source, config.source_record_offset)?;

    let mut output = LineWriter::new(writer);
    output.write(
        &header(&context, &source_header, source_nodes, config),
        true,
    )?;
    let mut consumed = 0u64;
    while consumed < config.source_record_count {
        let count =
            (config.source_record_count - consumed).min(config.batch_vertices as u64) as usize;
        let mut nodes = Vec::with_capacity(count);
        for _ in 0..count {
            nodes.push(source.read::<NodeRecord>(false)?);
        }
        let expanded: Vec<Result<Vec<AdjacencyEdgeRecord>, String>> = nodes
            .par_iter()
            .map(|node| expand_node(&context, node))
            .collect();
        for edges in expanded {
            for edge in edges? {
                output.write(&edge, true)?;
            }
        }
        consumed += count as u64;
        eprintln!(
            "p61-isogeny-adjacency: source records {}..{} of {} completed",
            config.source_record_offset,
            config.source_record_offset + consumed,
            config.source_record_end()
        );
    }
    let summary = summary(config, output.chain)?;
    output.write(&summary, false)?;
    Ok(AdjacencyGenerationReceipt {
        schema: ADJACENCY_RECEIPT_SCHEMA.into(),
        status: "pass".into(),
        producer_source_commit: config.producer_source_commit.clone(),
        elapsed_ms: started.elapsed().as_millis(),
        uncompressed_bytes: output.bytes,
        summary,
    })
}

pub fn verify_adjacency_jsonl<SR: BufRead, AR: BufRead>(
    source_reader: SR,
    adjacency_reader: AR,
    verifier_source_commit: &str,
    expected_source_segment_sha256: &str,
    expected_source_verification_receipt_sha256: &str,
) -> Result<AdjacencyVerificationReceipt, String> {
    if !is_lower_hex(verifier_source_commit, 40) {
        return Err("verifier source commit must be 40 lowercase hexadecimal digits".into());
    }
    let started = Instant::now();
    let context = AdjacencyContext::build()?;
    let mut source = LineReader::new(source_reader);
    let source_header: SegmentHeader = source.read(false)?;
    let mut adjacency = LineReader::new(adjacency_reader);
    let recorded_header: AdjacencyHeader = adjacency.read(true)?;
    if recorded_header.source_segment_sha256 != expected_source_segment_sha256
        || recorded_header.source_verification_receipt_sha256
            != expected_source_verification_receipt_sha256
    {
        return Err("adjacency header is bound to different source evidence".into());
    }
    let config = AdjacencyConfig {
        source_record_offset: recorded_header.source_record_offset,
        source_record_count: recorded_header.source_record_count,
        batch_vertices: recorded_header.batch_vertices,
        producer_source_commit: recorded_header.producer_source_commit.clone(),
        class_number_receipt_sha256: recorded_header.class_number_receipt_sha256.clone(),
        source_segment_sha256: recorded_header.source_segment_sha256.clone(),
        source_verification_receipt_sha256: recorded_header
            .source_verification_receipt_sha256
            .clone(),
    };
    config.validate()?;
    let source_nodes = validate_source_header(&source_header, &config)?;
    if recorded_header != header(&context, &source_header, source_nodes, &config) {
        return Err("adjacency header differs from the frozen p61 contract".into());
    }
    skip_source_records(&mut source, config.source_record_offset)?;

    let mut consumed = 0u64;
    while consumed < config.source_record_count {
        let count =
            (config.source_record_count - consumed).min(config.batch_vertices as u64) as usize;
        let mut nodes = Vec::with_capacity(count);
        for _ in 0..count {
            nodes.push(source.read::<NodeRecord>(false)?);
        }
        let mut recorded = Vec::with_capacity(count);
        for _ in 0..count {
            let mut edges = Vec::with_capacity(DIRECTED_EDGES_PER_VERTEX as usize);
            for _ in 0..DIRECTED_EDGES_PER_VERTEX {
                edges.push(adjacency.read::<AdjacencyEdgeRecord>(true)?);
            }
            recorded.push(edges);
        }
        let verified: Vec<Result<(), String>> = nodes
            .par_iter()
            .zip(recorded.par_iter())
            .map(|(node, edges)| verify_node_edges(&context, node, edges))
            .collect();
        for result in verified {
            result?;
        }
        consumed += count as u64;
    }
    let recorded_summary: AdjacencySummary = adjacency.read(false)?;
    let expected_summary = summary(&config, adjacency.chain)?;
    if recorded_summary != expected_summary {
        return Err("adjacency summary or record-chain digest differs".into());
    }
    let bytes = adjacency.finish()?;
    Ok(AdjacencyVerificationReceipt {
        schema: ADJACENCY_RECEIPT_SCHEMA.into(),
        status: "pass".into(),
        producer_source_commit: config.producer_source_commit,
        verifier_source_commit: verifier_source_commit.into(),
        elapsed_ms: started.elapsed().as_millis(),
        uncompressed_bytes: bytes,
        summary: recorded_summary,
    })
}

pub fn plan_record() -> serde_json::Value {
    serde_json::json!({
        "schema": "ec-hardness/p61-selected-adjacency-plan/v1",
        "root": "p61-i0025",
        "vertices": EXPECTED_CLASS_VERTICES,
        "selected_degrees": SELECTED_DEGREES,
        "rational_roots_per_degree": EXPECTED_ROOT_COUNTS,
        "directed_edges_per_vertex": DIRECTED_EDGES_PER_VERTEX,
        "directed_edges": EXPECTED_CLASS_VERTICES * DIRECTED_EDGES_PER_VERTEX,
        "kernel_coefficients_per_vertex": KERNEL_COEFFICIENTS_PER_VERTEX,
        "source_contract": "only independently replayed and admitted vertex segments",
        "target_join": "target j-invariant is joined to the dense vertex index after the complete census j-index is sealed",
        "completion_gate": "every admitted vertex occurs in exactly one adjacency source window; every shard independently replays; the final j-index resolves every target; reciprocity and degree counts pass globally",
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    fn config() -> AdjacencyConfig {
        AdjacencyConfig {
            source_record_offset: 0,
            source_record_count: 10,
            batch_vertices: 4,
            producer_source_commit: "a".repeat(40),
            class_number_receipt_sha256: "b".repeat(64),
            source_segment_sha256: "c".repeat(64),
            source_verification_receipt_sha256: "d".repeat(64),
        }
    }

    #[test]
    fn selected_degrees_have_seven_directed_edges() {
        assert_eq!(
            EXPECTED_ROOT_COUNTS
                .iter()
                .copied()
                .map(u64::from)
                .sum::<u64>(),
            7
        );
        assert_eq!(EXPECTED_CLASS_VERTICES * 7, 2_506_660_128);
        assert_eq!(KERNEL_COEFFICIENTS_PER_VERTEX, 3 + 2 * 7 + 2 * 19 + 2 * 30);
    }

    #[test]
    fn config_rejects_empty_or_overflowing_windows() {
        let base = config();
        assert!(base.validate().is_ok());
        let mut empty = base.clone();
        empty.source_record_count = 0;
        assert!(empty.validate().is_err());
        let mut overflow = base;
        overflow.source_record_offset = u64::MAX;
        overflow.source_record_count = 2;
        assert!(overflow.validate().is_err());
    }

    #[test]
    fn summary_has_exact_degree_and_kernel_counts() {
        let summary = summary(&config(), [7; 32]).unwrap();
        assert_eq!(summary.directed_edges, 70);
        assert_eq!(summary.kernel_coefficients, 1_150);
        assert_eq!(summary.degree_counts[0].directed_edges, 10);
        assert_eq!(summary.degree_counts[3].directed_edges, 20);
    }

    #[test]
    fn absolute_delta_bit_length_handles_zero() {
        let field = Field::new(&BigUint::from(101u32)).unwrap();
        assert_eq!(
            absolute_delta_bits(&field, &field.from_u64(9), &field.from_u64(9)),
            0
        );
        assert_eq!(
            absolute_delta_bits(&field, &field.from_u64(9), &field.from_u64(1)),
            4
        );
    }
}
