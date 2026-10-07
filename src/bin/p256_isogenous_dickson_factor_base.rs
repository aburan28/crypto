//! Staged Dickson factor-base screen across a verified P-256 isogeny walk.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::curve_id;
use crypto_lib::cryptanalysis::ecbench::canonical::short_id;
use crypto_lib::cryptanalysis::isogeny_walk::curve::{self, Model};
use crypto_lib::cryptanalysis::isogeny_walk::field::{Fe, Field};
use crypto_lib::cryptanalysis::isogeny_walk::record::{ec1_identity, ec1_records};
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    self, WideCurveFacts, WideFactorBaseDump, WideFactorBaseFacts, WideFbPoint, DUMP_SCHEMA,
    FACTOR_BASE_SCHEMA, POINT_KEY_BYTES, POINT_KEY_ENCODING,
};
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use flate2::read::GzDecoder;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use rayon::prelude::*;
use serde::Serialize;
use serde_json::{json, Value};

const CERTIFICATE_SHA256: &str = "4c4798fe791c8f603a248782a861f3488216af7476c814bff59f007f87c0d162";
const ROUND27_SHA256: &str = "256271ebb6f1c3c736a01da4358f01a6822c4315bb2081adcba1dd05b1d847a5";
const ROOT_CURVE: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const ROOT_FACTOR_BASE: &str = "FB1h2f8621cda105";
const ROOT_EXPONENT_HEX: &str = "2b6fdc73dc04e7667129";
const DEPTH: u32 = 18;
const FIBRE_SIZE: usize = 1 << DEPTH;
const EXPECTED_EDGES: usize = 4_096;
const RHO_S: f64 = 1.3;

#[derive(Parser)]
#[command(about = "Screen the selected Dickson fibre on verified P-256 isogenous models")]
struct Cli {
    #[arg(long)]
    certificate: PathBuf,
    #[arg(long)]
    round27: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
    #[arg(long)]
    dump: Option<PathBuf>,
}

#[derive(Clone)]
struct CurveModel {
    node: usize,
    a: BigUint,
    b: BigUint,
    j: BigUint,
    gx: BigUint,
    gy: BigUint,
}

#[derive(Clone, PartialEq, Eq)]
struct Fp2 {
    re: BigUint,
    im: BigUint,
}

impl Fp2 {
    fn one() -> Self {
        Self {
            re: BigUint::one(),
            im: BigUint::zero(),
        }
    }

    fn mul(&self, other: &Self, p: &BigUint) -> Self {
        let ac = (&self.re * &other.re) % p;
        let bd = (&self.im * &other.im) % p;
        let ad = (&self.re * &other.im) % p;
        let bc = (&self.im * &other.re) % p;
        Self {
            re: if ac >= bd { ac - bd } else { ac + p - bd },
            im: (ad + bc) % p,
        }
    }

    fn pow(&self, exponent: &BigUint, p: &BigUint) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();
        let mut value = exponent.clone();
        while !value.is_zero() {
            if value.bit(0) {
                result = result.mul(&base, p);
            }
            base = base.mul(&base, p);
            value >>= 1usize;
        }
        result
    }
}

#[derive(Clone, Serialize)]
struct StageRow {
    node: usize,
    sample_depth: u32,
    samples: u64,
    liftable: u64,
    zero_rhs: u64,
}

#[derive(Serialize)]
struct Stage {
    depth: u32,
    models_evaluated: u64,
    samples_per_model: u64,
    selected_models: Vec<usize>,
    rows_sha256: String,
    rows: Vec<StageRow>,
}

#[derive(Clone, Serialize)]
struct FullCandidate {
    node: usize,
    icv1_slug: String,
    icv1: String,
    j: String,
    columns: u64,
    signed_points: u64,
    poisson_mean_m17: f64,
    poisson_success_m17: f64,
    cross_colour_s: f64,
    ratio_to_rho: f64,
    projected_138031_rows_log2_operations: f64,
    projected_cost_per_usable_relation_log2_operations: f64,
    sparse_linear_algebra_lower_bound_log2_operations: f64,
    complete_projected_lower_bound_log2_operations: f64,
    projected_peak_two_list_log2_bytes: f64,
    zero_rhs: u64,
    selected_by_ladder: bool,
}

#[derive(Default, Serialize)]
struct OperationCounts {
    rhs_field_multiplications: u64,
    legendre_tests: u64,
    legendre_squarings: u64,
    legendre_multiplications: u64,
    sqrt_exponentiations: u64,
    sqrt_squarings: u64,
    sqrt_multiplications: u64,
    sqrt_verification_squarings: u64,
    order_witness_scalar_multiplications: u64,
}

#[derive(Serialize)]
struct NativeFactorBase {
    fb_id: String,
    fb_sha256: String,
    family: String,
    curve_slug: String,
    curve_icv1: String,
    curve_ec1: String,
    curve_uid: String,
    walk_node: usize,
    columns: u64,
    signed_points: u64,
    points_sha256: String,
    point_key_encoding: String,
    all_points_on_curve: bool,
    all_points_nonidentity: bool,
    unique_up_to_sign: bool,
    point_lift_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    rebuild_identity_verified: bool,
    dump_schema: String,
    dump_bytes: u64,
    dump_sha256: String,
    logical_peak_bytes: u64,
}

#[derive(Serialize)]
struct CertificateAudit {
    compressed_sha256: String,
    header_schema: String,
    root_curve: String,
    degree: u64,
    requested_steps: u64,
    edges: u64,
    models: u64,
    continuity_failures: u64,
    missing_order_witnesses: u64,
    invalid_order_witnesses: u64,
    distinct_j_invariants: u64,
    generic_j_models: u64,
    j_zero_or_1728_models: u64,
    all_checks_passed: bool,
}

#[derive(Serialize)]
struct DegreeBoundary {
    short_weierstrass_s3_total_degree: u32,
    universal_degree_four_terms: Vec<String>,
    all_models_retain_universal_leading_support: bool,
    isogeny_degree: u32,
    isogeny_pullback_is_degree_reduction: bool,
    structured_residual_degree_of_regularity: Option<u32>,
    degree_gate_passed: bool,
    classification: String,
}

#[derive(Serialize)]
struct PromotionGates {
    zero_false_positives: bool,
    zero_false_negatives: bool,
    exact_point_lifts_and_identity_rebuild: bool,
    structured_residual_degree_at_most_5: bool,
    relation_collection_below_2_120: bool,
    per_usable_row_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    complete_cost_below_rho: bool,
    exact_non_generic_log_transport: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    root_curve: String,
    root_factor_base: String,
    relation_model: String,
    dependencies: BTreeMap<String, String>,
    certificate_audit: CertificateAudit,
    fibre: BTreeMap<String, String>,
    stage8: Stage,
    stage12: Stage,
    full_depth_candidates: Vec<FullCandidate>,
    best_candidate: FullCandidate,
    native_factor_base: NativeFactorBase,
    degree_boundary: DegreeBoundary,
    operations: OperationCounts,
    promotion_gates: PromotionGates,
    full_depth_screen_exhaustive_over_all_models: bool,
    discarded_prefix_branches_counted_as_exhaustive: bool,
    full_depth_unplanted_attempted: bool,
    dominant_obstruction: String,
    decision: String,
}

fn parse_hex(value: &str) -> Result<BigUint, String> {
    BigUint::parse_bytes(value.trim_start_matches("0x").as_bytes(), 16)
        .ok_or_else(|| format!("invalid hex integer `{value}`"))
}

fn file_sha256(path: &Path) -> Result<String, String> {
    fs::read(path)
        .map(|bytes| hex::encode(sha256(&bytes)))
        .map_err(|error| format!("{}: {error}", path.display()))
}

fn require_hash(path: &Path, expected: &str, label: &str) -> Result<String, String> {
    let digest = file_sha256(path)?;
    if digest != expected {
        return Err(format!(
            "{label} hash mismatch: expected {expected}, got {digest}"
        ));
    }
    Ok(digest)
}

fn target_curve(value: &Value, node: usize) -> Result<CurveModel, String> {
    let target = value
        .get("target")
        .or_else(|| value.get("curve"))
        .ok_or_else(|| format!("node {node} has no curve object"))?;
    let witness = target
        .get("order_witness")
        .ok_or_else(|| format!("node {node} has no order witness"))?;
    Ok(CurveModel {
        node,
        a: parse_hex(
            target
                .get("a")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("node {node} has no a"))?,
        )?,
        b: parse_hex(
            target
                .get("b")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("node {node} has no b"))?,
        )?,
        j: parse_hex(
            target
                .get("j")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("node {node} has no j"))?,
        )?,
        gx: parse_hex(
            witness
                .get("x")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("node {node} witness has no x"))?,
        )?,
        gy: parse_hex(
            witness
                .get("y")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("node {node} witness has no y"))?,
        )?,
    })
}

fn audit_certificate(
    path: &Path,
    field: &Field,
    order: &BigUint,
) -> Result<(CertificateAudit, Vec<CurveModel>), String> {
    let decoder = GzDecoder::new(
        fs::File::open(path).map_err(|error| format!("{}: {error}", path.display()))?,
    );
    let mut header_schema = String::new();
    let mut root_curve = String::new();
    let mut degree = 0u64;
    let mut requested_steps = 0u64;
    let mut models = Vec::with_capacity(EXPECTED_EDGES + 1);
    let mut edge_count = 0usize;
    let mut continuity_failures = 0u64;
    let mut previous_j: Option<String> = None;
    let mut census_seen = false;
    for line in BufReader::new(decoder).lines() {
        let line = line.map_err(|error| error.to_string())?;
        let value: Value = serde_json::from_str(&line).map_err(|error| error.to_string())?;
        match value.get("record").and_then(Value::as_str) {
            Some("header") => {
                header_schema = value
                    .get("schema")
                    .and_then(Value::as_str)
                    .unwrap_or("")
                    .into();
                root_curve = value
                    .get("curve")
                    .and_then(Value::as_str)
                    .unwrap_or("")
                    .into();
                degree = value.get("degree").and_then(Value::as_u64).unwrap_or(0);
                requested_steps = value
                    .get("requested_steps")
                    .and_then(Value::as_u64)
                    .unwrap_or(0);
            }
            Some("census") => census_seen = true,
            Some("start") => {
                let model = target_curve(&value, 0)?;
                previous_j = Some(format!("{:064x}", model.j));
                models.push(model);
            }
            Some("edge") => {
                let index = value
                    .get("index")
                    .and_then(Value::as_u64)
                    .ok_or("edge has no index")? as usize;
                if index != edge_count {
                    continuity_failures += 1;
                }
                let source_j = value.get("source_j").and_then(Value::as_str).unwrap_or("");
                if previous_j.as_deref() != Some(source_j) {
                    continuity_failures += 1;
                }
                let model = target_curve(&value, index + 1)?;
                previous_j = Some(format!("{:064x}", model.j));
                models.push(model);
                edge_count += 1;
            }
            Some("footer") | Some("summary") => {}
            Some(other) => return Err(format!("unexpected certificate record `{other}`")),
            None => return Err("certificate row has no record tag".into()),
        }
    }
    if header_schema != "p256.isogeny-walk-certificate/v1"
        || root_curve != ROOT_CURVE
        || degree != 11
        || requested_steps != EXPECTED_EDGES as u64
        || edge_count != EXPECTED_EDGES
        || models.len() != EXPECTED_EDGES + 1
        || !census_seen
    {
        return Err("certificate shape differs from the frozen walk".into());
    }
    let invalid_order_witnesses = models
        .par_iter()
        .filter(|entry| {
            let model = Model {
                a: field.from_big(&entry.a),
                b: field.from_big(&entry.b),
            };
            let x = field.from_big(&entry.gx);
            let y = field.from_big(&entry.gy);
            !model.on_curve(field, &x, &y)
                || curve::scalar_mul(field, &model, &x, &y, order).is_some()
        })
        .count() as u64;
    let distinct_j = models
        .iter()
        .map(|entry| entry.j.clone())
        .collect::<BTreeSet<_>>()
        .len() as u64;
    let special_j = models
        .iter()
        .filter(|entry| entry.j.is_zero() || entry.j == BigUint::from(1728u32))
        .count() as u64;
    let all_checks_passed = continuity_failures == 0
        && invalid_order_witnesses == 0
        && distinct_j == models.len() as u64
        && special_j == 0;
    if !all_checks_passed {
        return Err("certificate audit failed".into());
    }
    Ok((
        CertificateAudit {
            compressed_sha256: CERTIFICATE_SHA256.into(),
            header_schema,
            root_curve,
            degree,
            requested_steps,
            edges: edge_count as u64,
            models: models.len() as u64,
            continuity_failures,
            missing_order_witnesses: 0,
            invalid_order_witnesses,
            distinct_j_invariants: distinct_j,
            generic_j_models: models.len() as u64 - special_j,
            j_zero_or_1728_models: special_j,
            all_checks_passed,
        },
        models,
    ))
}

fn fibre_abscissae(p: &BigUint) -> Result<Vec<BigUint>, String> {
    let seed = Fp2 {
        re: BigUint::one(),
        im: BigUint::from(5u8),
    };
    let torus_generator = seed.pow(&(p - 1u8), p);
    let odd_part = (p + 1u8) >> 96usize;
    let g96 = torus_generator.pow(&odd_part, p);
    if g96.pow(&(BigUint::one() << 96usize), p) != Fp2::one()
        || g96.pow(&(BigUint::one() << 95usize), p) == Fp2::one()
    {
        return Err("norm-one generator did not have exact order 2^96".into());
    }
    let root_exponent = parse_hex(ROOT_EXPONENT_HEX)?;
    let root = g96.pow(&root_exponent, p);
    let step = g96.pow(&(BigUint::one() << (96usize - DEPTH as usize)), p);
    let mut z = root;
    let mut out = Vec::with_capacity(FIBRE_SIZE);
    for _ in 0..FIBRE_SIZE {
        out.push((BigUint::from(2u8) * &z.re) % p);
        z = z.mul(&step, p);
    }
    if out.iter().cloned().collect::<BTreeSet<_>>().len() != FIBRE_SIZE {
        return Err("Dickson fibre abscissae are not distinct".into());
    }
    Ok(out)
}

fn sample_indices(depth: u32) -> Vec<usize> {
    let stride = 1usize << (DEPTH - depth);
    (0..FIBRE_SIZE).step_by(stride).collect()
}

fn lift_count(field: &Field, model: &CurveModel, xs: &[Fe], indices: &[usize]) -> (u64, u64) {
    let curve = Model {
        a: field.from_big(&model.a),
        b: field.from_big(&model.b),
    };
    let mut liftable = 0u64;
    let mut zero = 0u64;
    for index in indices {
        let rhs = curve.rhs(field, &xs[*index]);
        match field.legendre(&rhs) {
            0 => {
                liftable += 1;
                zero += 1;
            }
            1 => liftable += 1,
            _ => {}
        }
    }
    (liftable, zero)
}

fn digest_rows(rows: &[StageRow]) -> String {
    let mut bytes = Vec::with_capacity(rows.len() * 32);
    for row in rows {
        bytes.extend((row.node as u64).to_be_bytes());
        bytes.extend(u64::from(row.sample_depth).to_be_bytes());
        bytes.extend(row.samples.to_be_bytes());
        bytes.extend(row.liftable.to_be_bytes());
        bytes.extend(row.zero_rhs.to_be_bytes());
    }
    hex::encode(sha256(&bytes))
}

fn stage(
    depth: u32,
    retain: usize,
    models: &[CurveModel],
    model_indices: &[usize],
    field: &Field,
    xs: &[Fe],
) -> Stage {
    let indices = sample_indices(depth);
    let mut rows = model_indices
        .par_iter()
        .map(|node| {
            let (liftable, zero_rhs) = lift_count(field, &models[*node], xs, &indices);
            StageRow {
                node: *node,
                sample_depth: depth,
                samples: indices.len() as u64,
                liftable,
                zero_rhs,
            }
        })
        .collect::<Vec<_>>();
    rows.sort_by_key(|row| row.node);
    let mut ranked = rows.clone();
    ranked.sort_by(|left, right| {
        right
            .liftable
            .cmp(&left.liftable)
            .then_with(|| left.node.cmp(&right.node))
    });
    let selected_models = ranked.iter().take(retain).map(|row| row.node).collect();
    Stage {
        depth,
        models_evaluated: rows.len() as u64,
        samples_per_model: indices.len() as u64,
        selected_models,
        rows_sha256: digest_rows(&rows),
        rows,
    }
}

fn disjoint_probability(columns: u64) -> f64 {
    (0..9).fold(1.0, |probability, offset| {
        probability * (columns - 8 - offset) as f64 / (columns - offset) as f64
    })
}

fn full_candidate(
    model: &CurveModel,
    columns: u64,
    zero_rhs: u64,
    selected: bool,
    p: &BigUint,
    n: &BigUint,
) -> Result<FullCandidate, String> {
    let identity = curve_id::prime(p, &model.a, &model.b, n)
        .ok_or_else(|| format!("node {} has no ICV1 identity", model.node))?;
    let metrics = p256_dickson_factor_base::relation_metrics(columns, 17)?;
    let cross_colour_s =
        (std::f64::consts::PI * columns as f64 / disjoint_probability(columns)).sqrt();
    let log2_sqrt_n = n
        .to_f64()
        .ok_or("P-256 order did not convert to f64")?
        .log2()
        / 2.0;
    let relation_rows = 138_031f64;
    let rows_s = (std::f64::consts::PI * relation_rows / disjoint_probability(columns)).sqrt();
    let rows_log2 = log2_sqrt_n + rows_s.log2();
    let per_relation_log2 = rows_log2 - relation_rows.log2();
    let classes = columns as f64;
    let nonzeros = relation_rows * 17.0;
    let sparse_la_log2 = (2.0 * classes * nonzeros + classes * classes).log2();
    let maximum = rows_log2.max(sparse_la_log2);
    let complete_log2 =
        maximum + (2f64.powf(rows_log2 - maximum) + 2f64.powf(sparse_la_log2 - maximum)).log2();
    let peak_bytes_log2 =
        log2_sqrt_n + (classes / disjoint_probability(columns)).sqrt().log2() + 6.0;
    Ok(FullCandidate {
        node: model.node,
        icv1_slug: identity.slug,
        icv1: identity.icv1,
        j: model.j.to_string(),
        columns,
        signed_points: 2 * columns,
        poisson_mean_m17: metrics.poisson_mean,
        poisson_success_m17: metrics.poisson_success,
        cross_colour_s,
        ratio_to_rho: cross_colour_s / RHO_S,
        projected_138031_rows_log2_operations: rows_log2,
        projected_cost_per_usable_relation_log2_operations: per_relation_log2,
        sparse_linear_algebra_lower_bound_log2_operations: sparse_la_log2,
        complete_projected_lower_bound_log2_operations: complete_log2,
        projected_peak_two_list_log2_bytes: peak_bytes_log2,
        zero_rhs,
        selected_by_ladder: selected,
    })
}

fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn point_key_bytes(x: &BigUint, high_sign: bool) -> Result<[u8; POINT_KEY_BYTES], String> {
    let key = ((x + BigUint::one()) << 1usize) + BigUint::from(high_sign);
    let raw = key.to_bytes_be();
    if raw.len() > POINT_KEY_BYTES {
        return Err("P-256 point key exceeded 33 bytes".into());
    }
    let mut out = [0u8; POINT_KEY_BYTES];
    out[POINT_KEY_BYTES - raw.len()..].copy_from_slice(&raw);
    Ok(out)
}

fn materialize(
    model: &CurveModel,
    candidate: &FullCandidate,
    p: &BigUint,
    n: &BigUint,
    field: &Field,
    xs_big: &[BigUint],
    xs: &[Fe],
    dump_path: Option<&Path>,
    operations: &mut OperationCounts,
) -> Result<NativeFactorBase, String> {
    let curve = Model {
        a: field.from_big(&model.a),
        b: field.from_big(&model.b),
    };
    let sqrt_exponent = (p + 1u8) >> 2usize;
    let mut columns = Vec::with_capacity(candidate.columns as usize);
    let mut failures = 0u64;
    for (index, x) in xs.iter().enumerate() {
        let rhs = curve.rhs(field, x);
        match field.legendre(&rhs) {
            -1 => continue,
            0 => {
                failures += 1;
                continue;
            }
            _ => {}
        }
        let y = field.pow(&rhs, &sqrt_exponent);
        operations.sqrt_exponentiations += 1;
        if field.sqr(&y) != rhs {
            failures += 1;
            continue;
        }
        operations.sqrt_verification_squarings += 1;
        let y_big = field.to_big(&y);
        let low_y = y_big.clone().min(p - &y_big);
        columns.push((xs_big[index].clone(), low_y));
    }
    if failures != 0 || columns.len() != candidate.columns as usize {
        return Err(format!(
            "native materialisation had {failures} failures and {} columns, expected {}",
            columns.len(),
            candidate.columns
        ));
    }
    columns.sort_by(|left, right| left.0.cmp(&right.0));
    let unique_up_to_sign = columns.windows(2).all(|pair| pair[0].0 < pair[1].0);
    if !unique_up_to_sign {
        return Err("materialised factor base has duplicate abscissae".into());
    }
    let mut point_blob = Vec::with_capacity(columns.len() * 2 * POINT_KEY_BYTES);
    let mut points = Vec::with_capacity(columns.len() * 2);
    for (column, (x, low_y)) in columns.iter().enumerate() {
        point_blob.extend(point_key_bytes(x, false)?);
        point_blob.extend(point_key_bytes(x, true)?);
        points.push(WideFbPoint {
            x: lower_hex(x),
            y: lower_hex(low_y),
            col: column as u64,
            coef: "1".into(),
        });
        points.push(WideFbPoint {
            x: lower_hex(x),
            y: lower_hex(&(p - low_y)),
            col: column as u64,
            coef: (n - 1u8).to_string(),
        });
    }
    let points_sha256 = hex::encode(sha256(&point_blob));
    let identity =
        curve_id::prime(p, &model.a, &model.b, n).ok_or("best model ICV1 reconstruction failed")?;
    let (field_record, curve_record) = ec1_records(
        p,
        &model.a,
        &model.b,
        n,
        &BigUint::one(),
        (&model.gx, &model.gy),
    );
    let tag = if model.node == 0 { "p256" } else { "fp" };
    let (ec1, curve_uid) = ec1_identity(p, &field_record, &curve_record, tag);
    let mut params = BTreeMap::new();
    params.insert("depth".into(), DEPTH.to_string());
    params.insert("isogeny_degree".into(), "11".into());
    params.insert("isogeny_walk_node".into(), model.node.to_string());
    params.insert("negation_folded".into(), "true".into());
    params.insert("point_key_encoding".into(), POINT_KEY_ENCODING.into());
    params.insert("root_exponent".into(), format!("0x{ROOT_EXPONENT_HEX}"));
    let fb_identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": identity.slug,
        "family": "isogenous-dickson-torus",
        "params": params,
        "columns": candidate.columns,
        "signed_points": candidate.signed_points,
        "points_sha256": points_sha256,
    });
    let (fb_id, fb_sha256) = short_id("FB1", &fb_identity)?;
    let facts = WideFactorBaseFacts {
        fb_id: fb_id.clone(),
        fb_sha256: fb_sha256.clone(),
        family: "isogenous-dickson-torus".into(),
        params,
        description: format!(
            "depth-18 selected Dickson fibre on degree-11 P-256 isogeny-walk node {}",
            model.node
        ),
        signed_points: candidate.signed_points,
        abscissae: candidate.columns,
        columns: candidate.columns,
        dimension: None,
        points_sha256: points_sha256.clone(),
    };
    let curve_facts = WideCurveFacts {
        slug: identity.slug.clone(),
        icv1: identity.icv1.clone(),
        family: "prime".into(),
        construction: format!("verified degree-11 P-256 isogeny-walk node {}", model.node),
        field_degree: None,
        field_bits: p.bits() as u32,
        group_order: n.to_string(),
        r: n.to_string(),
        cofactor: "1".into(),
        generator: [lower_hex(&model.gx), lower_hex(&model.gy)],
        automorphisms_available: 2,
        registered: model.node == 0,
        ec1: ec1.clone(),
        curve_uid: curve_uid.clone(),
    };
    let dump = WideFactorBaseDump {
        schema: DUMP_SCHEMA.into(),
        curve: curve_facts,
        factor_base: facts,
        build_adds: 0,
        build_doubles: 0,
        points,
    };
    let mut dump_bytes = serde_json::to_vec_pretty(&dump).map_err(|error| error.to_string())?;
    dump_bytes.push(b'\n');
    let dump_sha256 = hex::encode(sha256(&dump_bytes));
    let rebuilt: WideFactorBaseDump =
        serde_json::from_slice(&dump_bytes).map_err(|error| error.to_string())?;
    if rebuilt != dump {
        return Err("wide dump serialization round trip failed".into());
    }
    if let Some(path) = dump_path {
        fs::write(path, &dump_bytes).map_err(|error| format!("{}: {error}", path.display()))?;
    }
    let rebuilt_identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": dump.curve.slug,
        "family": dump.factor_base.family,
        "params": dump.factor_base.params,
        "columns": dump.factor_base.columns,
        "signed_points": dump.factor_base.signed_points,
        "points_sha256": dump.factor_base.points_sha256,
    });
    let (rebuilt_id, rebuilt_sha) = short_id("FB1", &rebuilt_identity)?;
    let rebuild_identity_verified = rebuilt_id == fb_id && rebuilt_sha == fb_sha256;
    if !rebuild_identity_verified {
        return Err("FB1 identity rebuild failed".into());
    }
    Ok(NativeFactorBase {
        fb_id,
        fb_sha256,
        family: "isogenous-dickson-torus".into(),
        curve_slug: identity.slug,
        curve_icv1: identity.icv1,
        curve_ec1: ec1,
        curve_uid,
        walk_node: model.node,
        columns: candidate.columns,
        signed_points: candidate.signed_points,
        points_sha256,
        point_key_encoding: POINT_KEY_ENCODING.into(),
        all_points_on_curve: true,
        all_points_nonidentity: true,
        unique_up_to_sign,
        point_lift_failures: failures,
        false_positives: 0,
        false_negatives: 0,
        rebuild_identity_verified,
        dump_schema: DUMP_SCHEMA.into(),
        dump_bytes: dump_bytes.len() as u64,
        dump_sha256,
        logical_peak_bytes: (xs_big.len() * 64 + dump_bytes.len()) as u64,
    })
}

fn pow_counts(exponent: &BigUint) -> (u64, u64) {
    (
        exponent.bits(),
        exponent
            .to_bytes_le()
            .iter()
            .map(|byte| u64::from(byte.count_ones()))
            .sum(),
    )
}

fn run(cli: &Cli) -> Result<ResultFile, String> {
    let certificate_hash = require_hash(&cli.certificate, CERTIFICATE_SHA256, "certificate")?;
    let round27_hash = require_hash(&cli.round27, ROUND27_SHA256, "round 27")?;
    let p256 = CurveParams::p256();
    let field = Field::new(&p256.p).ok_or("P-256 field construction failed")?;
    let (certificate_audit, models) = audit_certificate(&cli.certificate, &field, &p256.n)?;
    let xs_big = fibre_abscissae(&p256.p)?;
    let xs = xs_big.iter().map(|x| field.from_big(x)).collect::<Vec<_>>();

    let all_models = (0..models.len()).collect::<Vec<_>>();
    let stage8 = stage(8, 64, &models, &all_models, &field, &xs);
    let stage12 = stage(12, 8, &models, &stage8.selected_models, &field, &xs);
    let mut full_nodes = stage12.selected_models.clone();
    if !full_nodes.contains(&0) {
        full_nodes.push(0);
    }
    full_nodes.sort_unstable();
    let all_indices = (0..FIBRE_SIZE).collect::<Vec<_>>();
    let full_counts = full_nodes
        .par_iter()
        .map(|node| {
            let (columns, zero_rhs) = lift_count(&field, &models[*node], &xs, &all_indices);
            (*node, columns, zero_rhs)
        })
        .collect::<Vec<_>>();
    if full_counts.iter().any(|(_, _, zero)| *zero != 0) {
        return Err("a certified odd-order curve produced rational 2-torsion".into());
    }
    let selected_set = stage12
        .selected_models
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let mut full_depth_candidates = full_counts
        .iter()
        .map(|(node, columns, zero_rhs)| {
            full_candidate(
                &models[*node],
                *columns,
                *zero_rhs,
                selected_set.contains(node),
                &p256.p,
                &p256.n,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    full_depth_candidates.sort_by_key(|row| row.node);
    let best_candidate = full_depth_candidates
        .iter()
        .max_by(|left, right| {
            left.columns
                .cmp(&right.columns)
                .then_with(|| right.node.cmp(&left.node))
        })
        .ok_or("full-depth candidate list is empty")?
        .clone();
    let best_model = &models[best_candidate.node];

    let mut operations = OperationCounts {
        order_witness_scalar_multiplications: models.len() as u64,
        ..OperationCounts::default()
    };
    let total_legendre = stage8.models_evaluated * stage8.samples_per_model
        + stage12.models_evaluated * stage12.samples_per_model
        + full_depth_candidates.len() as u64 * FIBRE_SIZE as u64
        + FIBRE_SIZE as u64;
    let legendre_exponent = (&p256.p - 1u8) >> 1usize;
    let (legendre_squares, legendre_multiplies) = pow_counts(&legendre_exponent);
    operations.rhs_field_multiplications = 2 * total_legendre;
    operations.legendre_tests = total_legendre;
    operations.legendre_squarings = total_legendre * legendre_squares;
    operations.legendre_multiplications = total_legendre * legendre_multiplies;
    let sqrt_exponent = (&p256.p + 1u8) >> 2usize;
    let (sqrt_squares, sqrt_multiplies) = pow_counts(&sqrt_exponent);
    operations.sqrt_squarings = best_candidate.columns * sqrt_squares;
    operations.sqrt_multiplications = best_candidate.columns * sqrt_multiplies;
    let native_factor_base = materialize(
        best_model,
        &best_candidate,
        &p256.p,
        &p256.n,
        &field,
        &xs_big,
        &xs,
        cli.dump.as_deref(),
        &mut operations,
    )?;
    if operations.sqrt_exponentiations != best_candidate.columns {
        return Err("sqrt operation accounting mismatch".into());
    }

    let mut dependencies = BTreeMap::new();
    dependencies.insert("isogeny_certificate_sha256".into(), certificate_hash);
    dependencies.insert("round27_sha256".into(), round27_hash);
    let mut fibre = BTreeMap::new();
    fibre.insert("depth".into(), DEPTH.to_string());
    fibre.insert("abscissae".into(), FIBRE_SIZE.to_string());
    fibre.insert("root_exponent".into(), format!("0x{ROOT_EXPONENT_HEX}"));
    fibre.insert(
        "abscissae_sha256".into(),
        hex::encode(sha256(
            &xs_big
                .iter()
                .flat_map(|x| {
                    let mut out = [0u8; 32];
                    let raw = x.to_bytes_be();
                    out[32 - raw.len()..].copy_from_slice(&raw);
                    out
                })
                .collect::<Vec<_>>(),
        )),
    );
    let degree_boundary = DegreeBoundary {
        short_weierstrass_s3_total_degree: 4,
        universal_degree_four_terms: vec![
            "x1^2*x3^2".into(),
            "-2*x1*x2*x3^2".into(),
            "x2^2*x3^2".into(),
            "x1^2*x2^2".into(),
        ],
        all_models_retain_universal_leading_support: true,
        isogeny_degree: 11,
        isogeny_pullback_is_degree_reduction: false,
        structured_residual_degree_of_regularity: None,
        degree_gate_passed: false,
        classification: "Changing a,b changes lower S3 coefficients but not its universal degree-four leading support. Pulling membership through a rational degree-11 isogeny adds map equations and denominator clearing; no completed dreg<=5 measurement exists on the selected model.".into(),
    };
    let correctness = certificate_audit.all_checks_passed
        && native_factor_base.point_lift_failures == 0
        && native_factor_base.rebuild_identity_verified;
    let promotion_gates = PromotionGates {
        zero_false_positives: true,
        zero_false_negatives: true,
        exact_point_lifts_and_identity_rebuild: correctness,
        structured_residual_degree_at_most_5: false,
        relation_collection_below_2_120: false,
        per_usable_row_below_2_103: false,
        materialized_storage_below_2_50: false,
        complete_cost_below_rho: false,
        exact_non_generic_log_transport: false,
        promoted: false,
    };
    Ok(ResultFile {
        schema: "p256.isogenous_dickson_factor_base_screen/v1".into(),
        root_curve: ROOT_CURVE.into(),
        root_factor_base: ROOT_FACTOR_BASE.into(),
        relation_model: "17 distinct variable columns; signed 8+9 cross-colour collision modulo global negation".into(),
        dependencies,
        certificate_audit,
        fibre,
        stage8,
        stage12,
        full_depth_candidates,
        best_candidate,
        native_factor_base,
        degree_boundary,
        operations,
        promotion_gates,
        full_depth_screen_exhaustive_over_all_models: false,
        discarded_prefix_branches_counted_as_exhaustive: false,
        full_depth_unplanted_attempted: false,
        dominant_obstruction: "Isogenous model selection can change quadratic-residue lift density, but it does not remove the universal degree-four S3 support or create factor-base logarithm transport. The full-column cross-colour boundary remains hundreds of times rho before a measured solver.".into(),
        decision: "No screened isogenous Dickson base is promoted. Any lift-density improvement is a factor-base inventory result only; degree, quotient, collection, storage, and complete-cost gates remain unmet.".into(),
    })
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(result) => {
            let mut bytes = match serde_json::to_vec_pretty(&result) {
                Ok(bytes) => bytes,
                Err(error) => {
                    eprintln!("serialization failed: {error}");
                    return ExitCode::FAILURE;
                }
            };
            bytes.push(b'\n');
            if let Some(path) = cli.out {
                if let Err(error) = fs::write(&path, &bytes) {
                    eprintln!("{}: {error}", path.display());
                    return ExitCode::FAILURE;
                }
            } else {
                print!("{}", String::from_utf8_lossy(&bytes));
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("P-256 isogenous Dickson screen failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nested_samples_have_frozen_sizes() {
        assert_eq!(sample_indices(8).len(), 256);
        assert_eq!(sample_indices(12).len(), 4_096);
        assert!(sample_indices(8)
            .iter()
            .all(|index| sample_indices(12).contains(index)));
    }

    #[test]
    fn cross_colour_boundary_matches_round22_scale() {
        let columns = 131_458;
        let s = (std::f64::consts::PI * columns as f64 / disjoint_probability(columns)).sqrt();
        assert!(s / RHO_S > 494.0);
        assert!(s / RHO_S < 495.0);
    }

    #[test]
    fn selected_root_exponent_fits_depth_domain() {
        assert!(parse_hex(ROOT_EXPONENT_HEX).unwrap() < (BigUint::one() << (96 - DEPTH)));
    }

    #[test]
    fn universal_s3_terms_have_nonzero_p256_coefficients() {
        let p = CurveParams::p256().p;
        assert_ne!(BigUint::one() % &p, BigUint::zero());
        assert_ne!((&p - 2u8) % &p, BigUint::zero());
    }
}
