//! Pair-table representation and timing audit for the P-256 round-33 selector.

use std::collections::BTreeSet;
use std::fs;
use std::hint::black_box;
use std::mem::size_of;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::time::Instant;

use clap::Parser;
use crypto_lib::ct_bignum::U256;
use crypto_lib::ecc::p256_field::P256FieldElement;
use crypto_lib::ecc::p256_point::{P256ProjectivePoint, A_FE, B_FE};
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::Value;

const ROUND33_SHA256: &str = "9931bccbd9e2821f4498f465ce65f83efb400e7fb3740a239bdb40d3275ffbd8";
const CURVE_SLUG: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const PAIR_ENTRIES: u64 = 34_562_148_612;
const PAIRS_PER_SEGMENT: usize = 8;
const FULL_REPETITIONS: usize = 9;
const POOL_SIZE: usize = 4_096;
const CHECK_SEGMENTS: usize = 4_096;
const TARGET_NS: u128 = 50_000_000;

#[derive(Parser)]
#[command(about = "Audit pair-table layouts for the P-256 biased restart selector")]
struct Cli {
    #[arg(long)]
    round33: PathBuf,
    #[arg(long)]
    out: Option<PathBuf>,
    #[arg(long)]
    quick: bool,
}

#[derive(Clone, Serialize)]
struct Dependency {
    round: u64,
    path: String,
    sha256: String,
    schema: String,
}

#[derive(Clone, Serialize)]
struct ParityBudget {
    cutoff_capacity: u64,
    mean_capacity: f64,
    local_oracle_ratio_to_rho: f64,
    round33_ratio_to_rho_with_target_correction: f64,
    allowed_extra_additions_per_sample: f64,
    allowed_extra_additions_per_segment: f64,
}

#[derive(Clone, Serialize)]
struct StorageRow {
    layout: String,
    entries: u64,
    bytes_per_entry: u64,
    projected_bytes: u64,
    projected_log2_bytes: f64,
    below_2_50_bytes: bool,
    directly_addable: bool,
}

#[derive(Clone, Serialize)]
struct Correctness {
    pool_points: u64,
    roundtrips_checked: u64,
    segments_checked: u64,
    compressed_failures: u64,
    affine_failures: u64,
    projective_failures: u64,
    coefficient_failures: u64,
    invalid_sqrt_or_curve: u64,
    projective_size_bytes: u64,
    projective_size_is_96: bool,
    pool_sha256: String,
    index_sha256: String,
    output_sha256: String,
}

#[derive(Clone, Serialize)]
struct OpTiming {
    wall_ns: u64,
    process_cpu_ns: u64,
    segments: u64,
    checksum: u64,
}

#[derive(Clone, Serialize)]
struct TimingRepetition {
    repetition: u64,
    execution_order: Vec<String>,
    control: OpTiming,
    direct_eight_additions: OpTiming,
    layout_total: OpTiming,
}

#[derive(Clone, Serialize)]
struct Interval {
    lower: f64,
    median: f64,
    upper: f64,
}

#[derive(Clone, Serialize)]
struct TimingRow {
    layout: String,
    table_entries: u64,
    entry_bytes: u64,
    materialized_bytes: u64,
    index_stream_sha256: String,
    repetitions: Vec<TimingRepetition>,
    direct_wall_ns_per_segment: Interval,
    layout_wall_ns_per_segment: Interval,
    incremental_addition_equivalents: Interval,
    projected_ratio_to_rho: Interval,
    upper_endpoint_within_parity_budget: bool,
}

#[derive(Clone, Serialize)]
struct BestLayoutRow {
    layout: String,
    table_entries: u64,
    median_incremental_addition_equivalents: f64,
    upper_incremental_addition_equivalents: f64,
    median_projected_ratio_to_rho: f64,
    upper_projected_ratio_to_rho: f64,
    storage_gate: bool,
    correctness_gate: bool,
    local_timing_gate: bool,
    survives_optimistic_local_lower_bound: bool,
}

#[derive(Clone, Serialize)]
struct LargestMeasuredRow {
    layout: String,
    table_entries: u64,
    materialized_bytes: u64,
    median_incremental_addition_equivalents: f64,
    upper_incremental_addition_equivalents: f64,
    median_projected_ratio_to_rho: f64,
    upper_projected_ratio_to_rho: f64,
    upper_endpoint_within_parity_budget: bool,
}

#[derive(Clone, Serialize)]
struct Host {
    rustc: String,
    target_arch: String,
    target_os: String,
    cpu_model: String,
    release_profile: bool,
    timing_repetitions: u64,
    confidence_interval_rule: String,
}

#[derive(Clone, Serialize)]
struct SemanticEvidence {
    dependency_sha256: String,
    curve: String,
    pair_entries: u64,
    parity_budget: ParityBudget,
    storage: Vec<StorageRow>,
    correctness: Correctness,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    family: String,
    dependency: Dependency,
    parity_budget: ParityBudget,
    storage: Vec<StorageRow>,
    correctness: Correctness,
    host: Host,
    timing_rows: Vec<TimingRow>,
    best_by_layout: Vec<BestLayoutRow>,
    largest_measured_by_layout: Vec<LargestMeasuredRow>,
    semantic_evidence_sha256: String,
    raw_timing_sha256: String,
    any_optimistic_local_lower_bound_survives: bool,
    full_table_materialized_or_measured: bool,
    complete_time_gate: bool,
    imported_end_to_end_gates: ImportedGates,
    full_depth_unplanted_attempted: bool,
    interpretation: String,
    decision: String,
}

#[derive(Serialize)]
struct ImportedGates {
    proved_p256_usable_relation_probability: bool,
    structured_residual_degree_at_most_5: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    non_generic_end_to_end: bool,
    promoted: bool,
}

#[derive(Copy, Clone)]
struct PointRecord {
    compressed: [u8; 33],
    affine: [u8; 64],
    projective: P256ProjectivePoint,
    coefficient: [u8; 32],
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
enum Layout {
    Compressed,
    Affine,
    Projective,
    Coefficient,
}

impl Layout {
    fn name(self) -> &'static str {
        match self {
            Layout::Compressed => "compressed33",
            Layout::Affine => "affine64",
            Layout::Projective => "projective96",
            Layout::Coefficient => "coefficient32",
        }
    }

    fn entry_bytes(self) -> usize {
        match self {
            Layout::Compressed => 33,
            Layout::Affine => 64,
            Layout::Projective => size_of::<P256ProjectivePoint>(),
            Layout::Coefficient => 32,
        }
    }

    fn maximum_power(self) -> u32 {
        match self {
            Layout::Compressed => 23,
            Layout::Affine | Layout::Projective => 22,
            Layout::Coefficient => 23,
        }
    }
}

enum Table {
    Compressed(Vec<[u8; 33]>),
    Affine(Vec<[u8; 64]>),
    Projective(Vec<P256ProjectivePoint>),
    Coefficient(Vec<[u8; 32]>),
}

impl Table {
    fn build(layout: Layout, entries: usize, pool: &[PointRecord]) -> Self {
        match layout {
            Layout::Compressed => Table::Compressed(
                (0..entries)
                    .map(|index| pool[index & (pool.len() - 1)].compressed)
                    .collect(),
            ),
            Layout::Affine => Table::Affine(
                (0..entries)
                    .map(|index| pool[index & (pool.len() - 1)].affine)
                    .collect(),
            ),
            Layout::Projective => Table::Projective(
                (0..entries)
                    .map(|index| pool[index & (pool.len() - 1)].projective)
                    .collect(),
            ),
            Layout::Coefficient => Table::Coefficient(
                (0..entries)
                    .map(|index| pool[index & (pool.len() - 1)].coefficient)
                    .collect(),
            ),
        }
    }
}

fn json_at<'a>(value: &'a Value, path: &[&str]) -> Result<&'a Value, String> {
    let mut cursor = value;
    for key in path {
        cursor = cursor
            .get(*key)
            .ok_or_else(|| format!("missing JSON path {}", path.join(".")))?;
    }
    Ok(cursor)
}

fn json_u64(value: &Value, path: &[&str]) -> Result<u64, String> {
    json_at(value, path)?
        .as_u64()
        .ok_or_else(|| format!("JSON path {} is not a u64", path.join(".")))
}

fn json_f64(value: &Value, path: &[&str]) -> Result<f64, String> {
    json_at(value, path)?
        .as_f64()
        .ok_or_else(|| format!("JSON path {} is not numeric", path.join(".")))
}

fn json_str<'a>(value: &'a Value, path: &[&str]) -> Result<&'a str, String> {
    json_at(value, path)?
        .as_str()
        .ok_or_else(|| format!("JSON path {} is not a string", path.join(".")))
}

fn load_dependency(path: &Path) -> Result<(Dependency, ParityBudget), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND33_SHA256 {
        return Err(format!(
            "round-33 dependency hash mismatch: expected {ROUND33_SHA256}, got {digest}"
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("invalid JSON: {error}"))?;
    if json_str(&value, &["curve"])? != CURVE_SLUG {
        return Err("round-33 curve mismatch".into());
    }
    let cutoff = json_u64(&value, &["cutoff_sweep", "selected", "cutoff_capacity"])?;
    let mean = json_f64(&value, &["cutoff_sweep", "selected", "mean_capacity"])?;
    let local = json_f64(&value, &["cutoff_sweep", "local_oracle_ratio_to_rho"])?;
    let ratio = json_f64(
        &value,
        &[
            "cutoff_sweep",
            "selected",
            "optimistic_ratio_to_rho_with_target_correction",
        ],
    )?;
    if cutoff != 219
        || (mean - 224.361_249_307_510_55).abs() > 1e-12
        || (local - 0.964_336_477_130_181).abs() > 1e-15
        || (ratio - 0.998_569_034_150_285_9).abs() > 1e-15
        || json_u64(&value, &["pair_table_accounting", "lookups_per_segment"])? != 8
        || json_u64(
            &value,
            &["pair_table_accounting", "setup_group_additions_per_segment"],
        )? != 8
        || json_u64(
            &value,
            &[
                "pair_table_accounting",
                "omitted_target_correction_additions_per_segment",
            ],
        )? != 1
    {
        return Err("round-33 frozen values do not match the registered protocol".into());
    }
    let allowed_per_sample = (1.0 - ratio) / local;
    let allowed_per_segment = (mean + 1.0) * allowed_per_sample;
    Ok((
        Dependency {
            round: 33,
            path: path.display().to_string(),
            sha256: digest,
            schema: json_str(&value, &["schema"])?.into(),
        },
        ParityBudget {
            cutoff_capacity: cutoff,
            mean_capacity: mean,
            local_oracle_ratio_to_rho: local,
            round33_ratio_to_rho_with_target_correction: ratio,
            allowed_extra_additions_per_sample: allowed_per_sample,
            allowed_extra_additions_per_segment: allowed_per_segment,
        },
    ))
}

fn fixed_32(value: &BigUint) -> [u8; 32] {
    let bytes = value.to_bytes_be();
    assert!(bytes.len() <= 32);
    let mut out = [0u8; 32];
    out[32 - bytes.len()..].copy_from_slice(&bytes);
    out
}

fn sqrt_exponent(curve: &CurveParams) -> [u8; 32] {
    fixed_32(&((&curve.p + BigUint::one()) >> 2usize))
}

fn pow_public(base: &P256FieldElement, exponent: &[u8; 32]) -> P256FieldElement {
    let mut result = P256FieldElement::ONE;
    for byte in exponent {
        for bit in (0..8).rev() {
            result = result.sqr();
            if (byte >> bit) & 1 == 1 {
                result = result.mul(base);
            }
        }
    }
    result
}

fn decompress(encoded: &[u8; 33], exponent: &[u8; 32]) -> Option<P256ProjectivePoint> {
    if encoded[0] != 2 && encoded[0] != 3 {
        return None;
    }
    let mut x_bytes = [0u8; 32];
    x_bytes.copy_from_slice(&encoded[1..]);
    let x_u = U256::from_bytes_be(&x_bytes);
    let x = P256FieldElement::from_canonical(&x_u);
    let rhs = x.sqr().mul(&x).add(&A_FE.mul(&x)).add(&B_FE);
    let mut y = pow_public(&rhs, exponent);
    if y.sqr() != rhs {
        return None;
    }
    let parity = y.to_bytes_be()[31] & 1;
    if parity != encoded[0] - 2 {
        y = y.neg();
    }
    Some(P256ProjectivePoint {
        x,
        y,
        z: P256FieldElement::ONE,
    })
}

fn affine_decode(encoded: &[u8; 64]) -> P256ProjectivePoint {
    let mut x = [0u8; 32];
    let mut y = [0u8; 32];
    x.copy_from_slice(&encoded[..32]);
    y.copy_from_slice(&encoded[32..]);
    P256ProjectivePoint::from_affine(&U256::from_bytes_be(&x), &U256::from_bytes_be(&y))
}

fn projective_eq(left: &P256ProjectivePoint, right: &P256ProjectivePoint) -> bool {
    let left_identity = bool::from(left.is_identity());
    let right_identity = bool::from(right.is_identity());
    if left_identity || right_identity {
        return left_identity == right_identity;
    }
    left.x.mul(&right.z) == right.x.mul(&left.z) && left.y.mul(&right.z) == right.y.mul(&left.z)
}

fn point_checksum(point: &P256ProjectivePoint) -> u64 {
    point.x.as_montgomery().0[0] ^ point.y.as_montgomery().0[1] ^ point.z.as_montgomery().0[2]
}

fn make_pool(count: usize) -> Vec<PointRecord> {
    assert!(count.is_power_of_two());
    let curve = CurveParams::p256();
    let generator = P256ProjectivePoint::from_textbook(&curve.generator());
    let mut pool = Vec::with_capacity(count);
    for index in 0..count {
        let mut input = b"p256-pair-table-round34-pool".to_vec();
        input.extend_from_slice(&(index as u64).to_be_bytes());
        let mut scalar = BigUint::from_bytes_be(&sha256(&input)) % &curve.n;
        if scalar.is_zero() {
            scalar = BigUint::one();
        }
        let point = generator.scalar_mul_ct(&scalar, curve.order_bits());
        let (x, y) = point.to_affine().expect("nonzero scalar below group order");
        let x_bytes = x.to_bytes_be();
        let y_bytes = y.to_bytes_be();
        let mut compressed = [0u8; 33];
        compressed[0] = 2 + (y_bytes[31] & 1);
        compressed[1..].copy_from_slice(&x_bytes);
        let mut affine = [0u8; 64];
        affine[..32].copy_from_slice(&x_bytes);
        affine[32..].copy_from_slice(&y_bytes);
        pool.push(PointRecord {
            compressed,
            affine,
            projective: point,
            coefficient: fixed_32(&scalar),
        });
    }
    pool
}

fn make_indices(entries: usize, segments: usize, domain: &[u8]) -> Vec<usize> {
    assert!(entries.is_power_of_two());
    let mut indices = Vec::with_capacity(segments * PAIRS_PER_SEGMENT);
    for counter in 0..segments * PAIRS_PER_SEGMENT {
        let mut input = Vec::with_capacity(domain.len() + 8);
        input.extend_from_slice(domain);
        input.extend_from_slice(&(counter as u64).to_be_bytes());
        let digest = sha256(&input);
        let mut word = [0u8; 8];
        word.copy_from_slice(&digest[..8]);
        indices.push((u64::from_be_bytes(word) as usize) & (entries - 1));
    }
    indices
}

fn digest_indices(indices: &[usize]) -> String {
    let mut bytes = Vec::with_capacity(indices.len() * 8);
    for index in indices {
        bytes.extend_from_slice(&(*index as u64).to_be_bytes());
    }
    hex::encode(sha256(&bytes))
}

fn run_control(indices: &[usize]) -> u64 {
    let mut state = 0u64;
    for chunk in indices.chunks_exact(PAIRS_PER_SEGMENT) {
        for index in chunk {
            state = state.rotate_left(7) ^ *index as u64;
        }
    }
    black_box(state)
}

fn run_direct(
    segments: usize,
    generator: &P256ProjectivePoint,
    points: &[P256ProjectivePoint; PAIRS_PER_SEGMENT],
) -> u64 {
    let mut checksum = 0u64;
    for segment in 0..segments {
        let mut acc = *generator;
        for offset in 0..PAIRS_PER_SEGMENT {
            acc = acc.add(black_box(
                &points[(segment + offset) & (PAIRS_PER_SEGMENT - 1)],
            ));
        }
        checksum ^= point_checksum(black_box(&acc));
    }
    black_box(checksum)
}

fn run_table(
    table: &Table,
    indices: &[usize],
    exponent: &[u8; 32],
    curve: &CurveParams,
    generator: &P256ProjectivePoint,
) -> u64 {
    let mut checksum = 0u64;
    match table {
        Table::Compressed(entries) => {
            for chunk in indices.chunks_exact(PAIRS_PER_SEGMENT) {
                let mut acc = *generator;
                for index in chunk {
                    let point = decompress(black_box(&entries[*index]), exponent)
                        .expect("materialized table contains valid compressed points");
                    acc = acc.add(&point);
                }
                checksum ^= point_checksum(black_box(&acc));
            }
        }
        Table::Affine(entries) => {
            for chunk in indices.chunks_exact(PAIRS_PER_SEGMENT) {
                let mut acc = *generator;
                for index in chunk {
                    let point = affine_decode(black_box(&entries[*index]));
                    acc = acc.add(&point);
                }
                checksum ^= point_checksum(black_box(&acc));
            }
        }
        Table::Projective(entries) => {
            for chunk in indices.chunks_exact(PAIRS_PER_SEGMENT) {
                let mut acc = *generator;
                for index in chunk {
                    acc = acc.add(black_box(&entries[*index]));
                }
                checksum ^= point_checksum(black_box(&acc));
            }
        }
        Table::Coefficient(entries) => {
            for chunk in indices.chunks_exact(PAIRS_PER_SEGMENT) {
                let mut scalar = BigUint::one();
                for index in chunk {
                    scalar += BigUint::from_bytes_be(black_box(&entries[*index]));
                }
                scalar %= &curve.n;
                let point = generator.scalar_mul_ct(&scalar, curve.order_bits());
                checksum ^= point_checksum(black_box(&point));
            }
        }
    }
    black_box(checksum)
}

fn process_cpu_ns() -> u64 {
    let mut timestamp = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    // SAFETY: `timestamp` is a valid writable timespec and the clock id is
    // supported on the Linux benchmark host.  A failure is represented by 0.
    let status = unsafe { libc::clock_gettime(libc::CLOCK_PROCESS_CPUTIME_ID, &mut timestamp) };
    if status != 0 {
        0
    } else {
        (timestamp.tv_sec as u64) * 1_000_000_000 + timestamp.tv_nsec as u64
    }
}

fn time_operation<F>(segments: usize, mut operation: F) -> OpTiming
where
    F: FnMut() -> u64,
{
    let cpu_start = process_cpu_ns();
    let wall_start = Instant::now();
    let checksum = operation();
    let wall_ns = wall_start.elapsed().as_nanos() as u64;
    let cpu_ns = process_cpu_ns().saturating_sub(cpu_start);
    OpTiming {
        wall_ns,
        process_cpu_ns: cpu_ns,
        segments: segments as u64,
        checksum,
    }
}

fn calibrated_segments<F>(mut operation: F, maximum: usize) -> usize
where
    F: FnMut(usize) -> u64,
{
    let trial = 64usize;
    let start = Instant::now();
    black_box(operation(trial));
    let elapsed = start.elapsed().as_nanos().max(1);
    let scaled = ((TARGET_NS * trial as u128).div_ceil(elapsed)) as usize;
    scaled.clamp(trial, maximum)
}

fn sorted_interval(values: &[f64]) -> Interval {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let (lower_index, median_index, upper_index) = if sorted.len() == 9 {
        (1, 4, 7)
    } else {
        (0, sorted.len() / 2, sorted.len() - 1)
    };
    Interval {
        lower: sorted[lower_index],
        median: sorted[median_index],
        upper: sorted[upper_index],
    }
}

fn timing_row(
    layout: Layout,
    entries: usize,
    pool: &[PointRecord],
    exponent: &[u8; 32],
    curve: &CurveParams,
    budget: &ParityBudget,
    repetitions: usize,
) -> TimingRow {
    let table = Table::build(layout, entries, pool);
    let generator = P256ProjectivePoint::from_textbook(&curve.generator());
    let direct_points: [P256ProjectivePoint; PAIRS_PER_SEGMENT] =
        std::array::from_fn(|index| pool[index].projective);
    let seed_indices = make_indices(entries, 131_072, layout.name().as_bytes());
    let direct_segments = calibrated_segments(
        |segments| run_direct(segments, &generator, &direct_points),
        131_072,
    );
    let layout_maximum = match layout {
        Layout::Compressed => 16_384,
        Layout::Coefficient => 32_768,
        Layout::Affine | Layout::Projective => 131_072,
    };
    let layout_segments = calibrated_segments(
        |segments| {
            run_table(
                &table,
                &seed_indices[..segments * PAIRS_PER_SEGMENT],
                exponent,
                curve,
                &generator,
            )
        },
        layout_maximum,
    );
    let needed = direct_segments.max(layout_segments);
    let indices = &seed_indices[..needed * PAIRS_PER_SEGMENT];

    black_box(run_control(&indices[..layout_segments * PAIRS_PER_SEGMENT]));
    black_box(run_direct(direct_segments, &generator, &direct_points));
    black_box(run_table(
        &table,
        &indices[..layout_segments * PAIRS_PER_SEGMENT],
        exponent,
        curve,
        &generator,
    ));

    let mut raw = Vec::with_capacity(repetitions);
    for repetition in 0..repetitions {
        let mut control = None;
        let mut direct = None;
        let mut total = None;
        let mut order = Vec::with_capacity(3);
        for offset in 0..3 {
            match (repetition + offset + layout as usize) % 3 {
                0 => {
                    order.push("control".into());
                    control = Some(time_operation(layout_segments, || {
                        run_control(&indices[..layout_segments * PAIRS_PER_SEGMENT])
                    }));
                }
                1 => {
                    order.push("direct_eight_additions".into());
                    direct = Some(time_operation(direct_segments, || {
                        run_direct(direct_segments, &generator, &direct_points)
                    }));
                }
                _ => {
                    order.push("layout_total".into());
                    total = Some(time_operation(layout_segments, || {
                        run_table(
                            &table,
                            &indices[..layout_segments * PAIRS_PER_SEGMENT],
                            exponent,
                            curve,
                            &generator,
                        )
                    }));
                }
            }
        }
        raw.push(TimingRepetition {
            repetition: repetition as u64,
            execution_order: order,
            control: control.expect("control timing"),
            direct_eight_additions: direct.expect("direct timing"),
            layout_total: total.expect("layout timing"),
        });
    }

    let direct_values: Vec<f64> = raw
        .iter()
        .map(|row| row.direct_eight_additions.wall_ns as f64 / direct_segments as f64)
        .collect();
    let layout_values: Vec<f64> = raw
        .iter()
        .map(|row| row.layout_total.wall_ns as f64 / layout_segments as f64)
        .collect();
    let direct_interval = sorted_interval(&direct_values);
    let layout_interval = sorted_interval(&layout_values);
    let incremental = Interval {
        lower: (PAIRS_PER_SEGMENT as f64 * (layout_interval.lower / direct_interval.upper - 1.0))
            .max(0.0),
        median: (PAIRS_PER_SEGMENT as f64
            * (layout_interval.median / direct_interval.median - 1.0))
            .max(0.0),
        upper: (PAIRS_PER_SEGMENT as f64 * (layout_interval.upper / direct_interval.lower - 1.0))
            .max(0.0),
    };
    let project_ratio = |extra: f64| {
        budget.round33_ratio_to_rho_with_target_correction
            + budget.local_oracle_ratio_to_rho * extra / (budget.mean_capacity + 1.0)
    };
    let projected = Interval {
        lower: project_ratio(incremental.lower),
        median: project_ratio(incremental.median),
        upper: project_ratio(incremental.upper),
    };
    TimingRow {
        layout: layout.name().into(),
        table_entries: entries as u64,
        entry_bytes: layout.entry_bytes() as u64,
        materialized_bytes: entries as u64 * layout.entry_bytes() as u64,
        index_stream_sha256: digest_indices(indices),
        repetitions: raw,
        direct_wall_ns_per_segment: direct_interval,
        layout_wall_ns_per_segment: layout_interval,
        incremental_addition_equivalents: incremental,
        projected_ratio_to_rho: projected.clone(),
        upper_endpoint_within_parity_budget: projected.upper <= 1.0,
    }
}

fn storage_rows() -> Vec<StorageRow> {
    [
        Layout::Compressed,
        Layout::Affine,
        Layout::Projective,
        Layout::Coefficient,
    ]
    .into_iter()
    .map(|layout| {
        let bytes = PAIR_ENTRIES * layout.entry_bytes() as u64;
        StorageRow {
            layout: layout.name().into(),
            entries: PAIR_ENTRIES,
            bytes_per_entry: layout.entry_bytes() as u64,
            projected_bytes: bytes,
            projected_log2_bytes: (bytes as f64).log2(),
            below_2_50_bytes: bytes < (1u64 << 50),
            directly_addable: layout == Layout::Projective,
        }
    })
    .collect()
}

fn append_point(bytes: &mut Vec<u8>, point: &P256ProjectivePoint) {
    if let Some((x, y)) = point.to_affine() {
        bytes.extend_from_slice(&x.to_bytes_be());
        bytes.extend_from_slice(&y.to_bytes_be());
    } else {
        bytes.extend_from_slice(&[0u8; 64]);
    }
}

fn correctness(pool: &[PointRecord], exponent: &[u8; 32], curve: &CurveParams) -> Correctness {
    let mut compressed_failures = 0u64;
    let mut affine_failures = 0u64;
    let mut projective_failures = 0u64;
    let mut coefficient_failures = 0u64;
    let mut invalid = 0u64;
    let mut pool_bytes = Vec::with_capacity(pool.len() * (33 + 64 + 32));
    for record in pool {
        pool_bytes.extend_from_slice(&record.compressed);
        pool_bytes.extend_from_slice(&record.affine);
        pool_bytes.extend_from_slice(&record.coefficient);
        match decompress(&record.compressed, exponent) {
            Some(point) if projective_eq(&point, &record.projective) => {}
            Some(_) => compressed_failures += 1,
            None => {
                compressed_failures += 1;
                invalid += 1;
            }
        }
        let affine = affine_decode(&record.affine);
        if !projective_eq(&affine, &record.projective) {
            affine_failures += 1;
        }
        let x = affine.x;
        let rhs = x.sqr().mul(&x).add(&A_FE.mul(&x)).add(&B_FE);
        if affine.y.sqr() != rhs {
            invalid += 1;
        }
    }

    let indices = make_indices(pool.len(), CHECK_SEGMENTS, b"round34-correctness");
    let generator = P256ProjectivePoint::from_textbook(&curve.generator());
    let mut output_bytes = Vec::with_capacity(CHECK_SEGMENTS * 64 * 4);
    for chunk in indices.chunks_exact(PAIRS_PER_SEGMENT) {
        let mut reference = generator;
        let mut compressed = generator;
        let mut affine = generator;
        let mut projective = generator;
        let mut scalar = BigUint::one();
        for index in chunk {
            let record = &pool[*index];
            reference = reference.add(&record.projective);
            match decompress(&record.compressed, exponent) {
                Some(point) => compressed = compressed.add(&point),
                None => invalid += 1,
            }
            affine = affine.add(&affine_decode(&record.affine));
            projective = projective.add(&record.projective);
            scalar += BigUint::from_bytes_be(&record.coefficient);
        }
        scalar %= &curve.n;
        let coefficient = generator.scalar_mul_ct(&scalar, curve.order_bits());
        if !projective_eq(&compressed, &reference) {
            compressed_failures += 1;
        }
        if !projective_eq(&affine, &reference) {
            affine_failures += 1;
        }
        if !projective_eq(&projective, &reference) {
            projective_failures += 1;
        }
        if !projective_eq(&coefficient, &reference) {
            coefficient_failures += 1;
        }
        append_point(&mut output_bytes, &compressed);
        append_point(&mut output_bytes, &affine);
        append_point(&mut output_bytes, &projective);
        append_point(&mut output_bytes, &coefficient);
    }
    Correctness {
        pool_points: pool.len() as u64,
        roundtrips_checked: pool.len() as u64,
        segments_checked: CHECK_SEGMENTS as u64,
        compressed_failures,
        affine_failures,
        projective_failures,
        coefficient_failures,
        invalid_sqrt_or_curve: invalid,
        projective_size_bytes: size_of::<P256ProjectivePoint>() as u64,
        projective_size_is_96: size_of::<P256ProjectivePoint>() == 96,
        pool_sha256: hex::encode(sha256(&pool_bytes)),
        index_sha256: digest_indices(&indices),
        output_sha256: hex::encode(sha256(&output_bytes)),
    }
}

fn host(repetitions: usize) -> Host {
    let rustc = Command::new("rustc")
        .arg("--version")
        .output()
        .ok()
        .filter(|output| output.status.success())
        .map(|output| String::from_utf8_lossy(&output.stdout).trim().to_string())
        .unwrap_or_else(|| "unavailable".into());
    let cpu_model = fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|text| {
            text.lines()
                .find(|line| line.starts_with("model name"))
                .and_then(|line| line.split_once(':'))
                .map(|(_, value)| value.trim().to_string())
        })
        .unwrap_or_else(|| "unavailable".into());
    Host {
        rustc,
        target_arch: std::env::consts::ARCH.into(),
        target_os: std::env::consts::OS.into(),
        cpu_model,
        release_profile: !cfg!(debug_assertions),
        timing_repetitions: repetitions as u64,
        confidence_interval_rule: if repetitions == 9 {
            "sorted ranks 1,4,7 (zero-indexed) from nine repetitions".into()
        } else {
            "quick-mode minimum, median, maximum".into()
        },
    }
}

fn execute(cli: Cli) -> Result<ResultFile, String> {
    let (dependency, budget) = load_dependency(&cli.round33)?;
    let curve = CurveParams::p256();
    let exponent = sqrt_exponent(&curve);
    let pool_count = if cli.quick { 64 } else { POOL_SIZE };
    let pool = make_pool(pool_count);
    let correctness = correctness(&pool, &exponent, &curve);
    let storage = storage_rows();
    let repetitions = if cli.quick { 3 } else { FULL_REPETITIONS };
    let mut timing_rows = Vec::new();
    for layout in [
        Layout::Compressed,
        Layout::Affine,
        Layout::Projective,
        Layout::Coefficient,
    ] {
        let powers: BTreeSet<u32> = if cli.quick {
            [6, 10].into_iter().collect()
        } else {
            [10, 15, 20, layout.maximum_power()].into_iter().collect()
        };
        for power in powers {
            timing_rows.push(timing_row(
                layout,
                1usize << power,
                &pool,
                &exponent,
                &curve,
                &budget,
                repetitions,
            ));
        }
    }
    let correctness_gate = correctness.compressed_failures == 0
        && correctness.affine_failures == 0
        && correctness.projective_failures == 0
        && correctness.coefficient_failures == 0
        && correctness.invalid_sqrt_or_curve == 0
        && correctness.projective_size_is_96;
    let mut best_by_layout = Vec::new();
    for layout in [
        Layout::Compressed,
        Layout::Affine,
        Layout::Projective,
        Layout::Coefficient,
    ] {
        let best = timing_rows
            .iter()
            .filter(|row| row.layout == layout.name())
            .min_by(|left, right| {
                left.incremental_addition_equivalents
                    .median
                    .total_cmp(&right.incremental_addition_equivalents.median)
            })
            .expect("each layout has timing rows");
        let storage_gate = storage
            .iter()
            .find(|row| row.layout == layout.name())
            .expect("each layout has storage")
            .below_2_50_bytes;
        let local_timing_gate = best.upper_endpoint_within_parity_budget;
        best_by_layout.push(BestLayoutRow {
            layout: layout.name().into(),
            table_entries: best.table_entries,
            median_incremental_addition_equivalents: best.incremental_addition_equivalents.median,
            upper_incremental_addition_equivalents: best.incremental_addition_equivalents.upper,
            median_projected_ratio_to_rho: best.projected_ratio_to_rho.median,
            upper_projected_ratio_to_rho: best.projected_ratio_to_rho.upper,
            storage_gate,
            correctness_gate,
            local_timing_gate,
            survives_optimistic_local_lower_bound: storage_gate
                && correctness_gate
                && local_timing_gate,
        });
    }
    let mut largest_measured_by_layout = Vec::new();
    for layout in [
        Layout::Compressed,
        Layout::Affine,
        Layout::Projective,
        Layout::Coefficient,
    ] {
        let row = timing_rows
            .iter()
            .filter(|row| row.layout == layout.name())
            .max_by_key(|row| row.table_entries)
            .expect("each layout has timing rows");
        largest_measured_by_layout.push(LargestMeasuredRow {
            layout: layout.name().into(),
            table_entries: row.table_entries,
            materialized_bytes: row.materialized_bytes,
            median_incremental_addition_equivalents: row.incremental_addition_equivalents.median,
            upper_incremental_addition_equivalents: row.incremental_addition_equivalents.upper,
            median_projected_ratio_to_rho: row.projected_ratio_to_rho.median,
            upper_projected_ratio_to_rho: row.projected_ratio_to_rho.upper,
            upper_endpoint_within_parity_budget: row.upper_endpoint_within_parity_budget,
        });
    }
    let semantic = SemanticEvidence {
        dependency_sha256: dependency.sha256.clone(),
        curve: CURVE_SLUG.into(),
        pair_entries: PAIR_ENTRIES,
        parity_budget: budget.clone(),
        storage: storage.clone(),
        correctness: correctness.clone(),
    };
    let semantic_sha = hex::encode(sha256(
        &serde_json::to_vec(&semantic).map_err(|error| error.to_string())?,
    ));
    let timing_sha = hex::encode(sha256(
        &serde_json::to_vec(&timing_rows).map_err(|error| error.to_string())?,
    ));
    let any_survives = best_by_layout
        .iter()
        .any(|row| row.survives_optimistic_local_lower_bound);
    Ok(ResultFile {
        schema: "p256-pair-table-audit-v1".into(),
        curve: CURVE_SLUG.into(),
        family: "round-33 biased-restart pair-table representation audit".into(),
        dependency,
        parity_budget: budget,
        storage,
        correctness,
        host: host(repetitions),
        timing_rows,
        best_by_layout,
        largest_measured_by_layout,
        semantic_evidence_sha256: semantic_sha,
        raw_timing_sha256: timing_sha,
        any_optimistic_local_lower_bound_survives: any_survives,
        full_table_materialized_or_measured: false,
        complete_time_gate: false,
        imported_end_to_end_gates: ImportedGates {
            proved_p256_usable_relation_probability: false,
            structured_residual_degree_at_most_5: false,
            relation_collection_below_2_120: false,
            per_usable_relation_below_2_103: false,
            non_generic_end_to_end: false,
            promoted: false,
        },
        full_depth_unplanted_attempted: false,
        interpretation: if any_survives {
            "At least one faithful cache-resident representation fits the storage gate and its candidate-favouring local lower-bound interval preserves the narrow round-33 stage margin. This is not a complete-time or end-to-end parity result: the full table was not materialized, the largest measured working sets fail the timing budget, and achieved P-256 coverage, degree, collection and non-generic status remain open."
                .into()
        } else {
            "No faithful representation preserves the round-33 cutoff-219 margin under the registered local timing screen. This rejects that concrete pair-table implementation, not all factor bases or selectors."
                .into()
        },
        decision: if any_survives {
            "preserve only the projective layout's cache-resident lower bound as a stage candidate; require a locality construction before promotion or any full-depth relation"
                .into()
        } else {
            "reject the round-33 pair-table implementation; do not promote or run a full-depth relation"
                .into()
        },
    })
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let out = cli.out.clone();
    match execute(cli).and_then(|result| {
        let bytes = serde_json::to_vec_pretty(&result).map_err(|error| error.to_string())?;
        if let Some(path) = out {
            fs::write(&path, &bytes).map_err(|error| format!("{}: {error}", path.display()))?;
        } else {
            println!("{}", String::from_utf8_lossy(&bytes));
        }
        Ok(())
    }) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p256_pair_table_accounting: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compressed_and_affine_roundtrip() {
        let curve = CurveParams::p256();
        let exponent = sqrt_exponent(&curve);
        let pool = make_pool(16);
        for record in &pool {
            assert!(projective_eq(
                &decompress(&record.compressed, &exponent).unwrap(),
                &record.projective
            ));
            assert!(projective_eq(
                &affine_decode(&record.affine),
                &record.projective
            ));
        }
    }

    #[test]
    fn storage_gate_and_projective_size() {
        assert_eq!(size_of::<P256ProjectivePoint>(), 96);
        assert!(storage_rows().iter().all(|row| row.below_2_50_bytes));
    }

    #[test]
    fn parity_budget_is_narrow() {
        let mean = 224.361_249_307_510_55;
        let local = 0.964_336_477_130_181;
        let ratio = 0.998_569_034_150_285_9;
        let budget = (mean + 1.0) * (1.0 - ratio) / local;
        assert!(budget > 0.33 && budget < 0.34);
    }
}
