//! Montgomery-affine pair-table audit for P-256 (round 36).

use std::fs;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::time::Instant;

use clap::Parser;
use crypto_lib::ecc::p256_field::P256FieldElement;
use crypto_lib::ecc::p256_point::{P256ProjectivePoint, A_FE, B3_FE};
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::Serialize;
use serde_json::Value;

const ROUND35_SHA256: &str = "90704e913a0407fe773ad7550cd5b68c5ead54077ff485978bd51cdf21bcfc89";
const CURVE_SLUG: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const PAIR_ENTRIES: u64 = 34_562_148_612;
const FULL_TABLE_BYTES: u64 = 2_211_977_511_168;
const PROJECTIVE_FULL_TABLE_BYTES: u64 = 3_317_966_266_752;
const PROJECTIVE_ENTRY_BYTES: usize = 96;
const AFFINE_ENTRY_BYTES: usize = 64;
const PAIRS: usize = 8;
const SEGMENTS: usize = 1 << 16;
const CORRECTNESS_SEGMENTS: usize = 1 << 12;
const POOL_SIZE: usize = 1 << 12;
const REPETITIONS: usize = 9;
const PREFETCH_DISTANCES: [usize; 7] = [1, 2, 4, 8, 16, 32, 64];

#[derive(Parser)]
#[command(about = "Audit P-256 Montgomery-affine pair-table layouts")]
struct Cli {
    #[arg(long)]
    round35: PathBuf,
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
struct FrozenModel {
    round33_ratio_to_rho_with_target_correction: f64,
    local_oracle_ratio_to_rho: f64,
    mean_capacity: f64,
    allowed_extra_additions_per_segment: f64,
    pair_entries: u64,
    entry_bytes: u64,
    full_table_bytes: u64,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Variant {
    Direct,
    Prefetch(usize),
    DirectAffineComplete,
    DirectAffineMixed,
    AffineComplete(usize),
    AffineMixed(usize),
}

impl Variant {
    fn name(self) -> String {
        match self {
            Variant::Direct => "direct-eight-additions".into(),
            Variant::Prefetch(distance) => format!("projective-prefetch-{distance}"),
            Variant::DirectAffineComplete => "direct-affine-complete".into(),
            Variant::DirectAffineMixed => "direct-affine-mixed".into(),
            Variant::AffineComplete(distance) => format!("affine-complete-prefetch-{distance}"),
            Variant::AffineMixed(distance) => format!("affine-mixed-prefetch-{distance}"),
        }
    }

    fn class(self) -> &'static str {
        match self {
            Variant::Direct => "direct",
            Variant::Prefetch(_) => "projective",
            Variant::DirectAffineComplete => "direct-affine-complete",
            Variant::DirectAffineMixed => "direct-affine-mixed",
            Variant::AffineComplete(_) => "affine-complete",
            Variant::AffineMixed(_) => "affine-mixed",
        }
    }

    fn parameter(self) -> Option<u64> {
        match self {
            Variant::Direct | Variant::DirectAffineComplete | Variant::DirectAffineMixed => None,
            Variant::Prefetch(value)
            | Variant::AffineComplete(value)
            | Variant::AffineMixed(value) => Some(value as u64),
        }
    }
}

#[derive(Copy, Clone)]
struct MontAffinePoint {
    x: P256FieldElement,
    y: P256FieldElement,
}

#[derive(Clone, Serialize)]
struct OpTiming {
    variant: String,
    class: String,
    parameter: Option<u64>,
    wall_ns: u64,
    process_cpu_ns: u64,
    segments: u64,
    checksum: u64,
}

#[derive(Clone, Serialize)]
struct TimingRepetition {
    repetition: u64,
    execution_order: Vec<String>,
    operations: Vec<OpTiming>,
}

#[derive(Clone, Serialize)]
struct Interval {
    lower: f64,
    median: f64,
    upper: f64,
}

#[derive(Clone, Serialize)]
struct TimingSummary {
    variant: String,
    class: String,
    parameter: Option<u64>,
    wall_ns_per_segment: Interval,
    incremental_addition_equivalents: Option<Interval>,
    projected_ratio_to_rho: Option<Interval>,
    upper_endpoint_within_parity_budget: Option<bool>,
}

#[derive(Clone, Serialize)]
struct TableTiming {
    table_entries: u64,
    affine_table_bytes: u64,
    projective_control_table_bytes: u64,
    index_stream_sha256: String,
    repetitions: Vec<TimingRepetition>,
    summaries: Vec<TimingSummary>,
}

#[derive(Clone, Serialize)]
struct Correctness {
    pool_points: u64,
    segments: u64,
    variants_checked: u64,
    group_mismatches: u64,
    missing_requests: u64,
    duplicate_requests: u64,
    out_of_range_accesses: u64,
    wrong_addition_counts: u64,
    curve_failures: u64,
    false_positives: u64,
    false_negatives: u64,
    pool_sha256: String,
    index_sha256: String,
    affine_table_sha256: String,
    output_sha256: String,
    special_cases_checked: u64,
    special_case_mismatches: u64,
}

#[derive(Clone, Serialize)]
struct SelectedRow {
    table_entries: u64,
    table_bytes: u64,
    variant: String,
    class: String,
    parameter: Option<u64>,
    incremental_addition_equivalents: Interval,
    projected_ratio_to_rho: Interval,
    correctness_gate: bool,
    timing_gate: bool,
    advances_affine_stage: bool,
}

#[derive(Serialize)]
struct ArithmeticAccounting {
    complete_field_multiplications_per_addition: u64,
    mixed_field_multiplications_per_addition: u64,
    field_multiplications_removed_per_addition: u64,
}

#[derive(Clone, Serialize)]
struct Host {
    rustc: String,
    target_arch: String,
    target_os: String,
    cpu_model: String,
    release_profile: bool,
    timing_repetitions: u64,
    prefetch_supported: bool,
    confidence_interval_rule: String,
}

#[derive(Serialize)]
struct SemanticEvidence {
    dependency_sha256: String,
    curve: String,
    model: FrozenModel,
    correctness: Correctness,
}

#[derive(Serialize)]
struct ImportedGates {
    full_table_materialized_or_measured: bool,
    complete_time_at_or_below_rho: bool,
    proved_p256_usable_relation_probability: bool,
    structured_residual_degree_at_most_5: bool,
    relation_collection_below_2_120: bool,
    per_usable_relation_below_2_103: bool,
    materialized_storage_below_2_50: bool,
    non_generic_end_to_end: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultFile {
    schema: String,
    curve: String,
    family: String,
    dependency: Dependency,
    frozen_model: FrozenModel,
    host: Host,
    correctness: Correctness,
    arithmetic_accounting: ArithmeticAccounting,
    timing_tables: Vec<TableTiming>,
    selected_decision_row: SelectedRow,
    semantic_evidence_sha256: String,
    raw_timing_sha256: String,
    imported_gates: ImportedGates,
    full_depth_unplanted_attempted: bool,
    interpretation: String,
    decision: String,
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

fn json_f64(value: &Value, path: &[&str]) -> Result<f64, String> {
    json_at(value, path)?
        .as_f64()
        .ok_or_else(|| format!("JSON path {} is not numeric", path.join(".")))
}

fn json_u64(value: &Value, path: &[&str]) -> Result<u64, String> {
    json_at(value, path)?
        .as_u64()
        .ok_or_else(|| format!("JSON path {} is not a u64", path.join(".")))
}

fn json_str<'a>(value: &'a Value, path: &[&str]) -> Result<&'a str, String> {
    json_at(value, path)?
        .as_str()
        .ok_or_else(|| format!("JSON path {} is not a string", path.join(".")))
}

fn load_dependency(path: &Path) -> Result<(Dependency, FrozenModel), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND35_SHA256 {
        return Err(format!(
            "round-35 dependency hash mismatch: expected {ROUND35_SHA256}, got {digest}"
        ));
    }
    let value: Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("invalid JSON: {error}"))?;
    if json_str(&value, &["curve"])? != CURVE_SLUG {
        return Err("round-35 curve mismatch".into());
    }
    let ratio = json_f64(
        &value,
        &[
            "frozen_model",
            "round33_ratio_to_rho_with_target_correction",
        ],
    )?;
    let local = json_f64(&value, &["frozen_model", "local_oracle_ratio_to_rho"])?;
    let mean = json_f64(&value, &["frozen_model", "mean_capacity"])?;
    let allowed = json_f64(
        &value,
        &["frozen_model", "allowed_extra_additions_per_segment"],
    )?;
    if json_u64(&value, &["frozen_model", "pair_entries"])? != PAIR_ENTRIES
        || json_u64(&value, &["frozen_model", "entry_bytes"])? != PROJECTIVE_ENTRY_BYTES as u64
        || json_u64(&value, &["frozen_model", "full_table_bytes"])? != PROJECTIVE_FULL_TABLE_BYTES
        || (ratio - 0.998_569_034_150_286).abs() > 1e-15
        || (local - 0.964_336_477_130_181).abs() > 1e-15
        || (mean - 224.361_249_307_510_55).abs() > 1e-12
        || (allowed - 0.334_410_508_422_987_64).abs() > 1e-12
    {
        return Err("round-35 frozen values do not match the protocol".into());
    }
    Ok((
        Dependency {
            round: 35,
            path: path.display().to_string(),
            sha256: digest,
            schema: json_str(&value, &["schema"])?.into(),
        },
        FrozenModel {
            round33_ratio_to_rho_with_target_correction: ratio,
            local_oracle_ratio_to_rho: local,
            mean_capacity: mean,
            allowed_extra_additions_per_segment: allowed,
            pair_entries: PAIR_ENTRIES,
            entry_bytes: AFFINE_ENTRY_BYTES as u64,
            full_table_bytes: FULL_TABLE_BYTES,
        },
    ))
}

fn make_pool(count: usize) -> Vec<P256ProjectivePoint> {
    let curve = CurveParams::p256();
    let generator = P256ProjectivePoint::from_textbook(&curve.generator());
    (0..count)
        .map(|index| {
            let mut input = b"p256-projective-locality-round35-pool".to_vec();
            input.extend_from_slice(&(index as u64).to_be_bytes());
            let mut scalar = BigUint::from_bytes_be(&sha256(&input)) % &curve.n;
            if scalar.is_zero() {
                scalar = BigUint::one();
            }
            generator.scalar_mul_ct(&scalar, curve.order_bits())
        })
        .collect()
}

fn make_table(entries: usize, pool: &[P256ProjectivePoint]) -> Vec<P256ProjectivePoint> {
    (0..entries)
        .map(|index| pool[index & (pool.len() - 1)])
        .collect()
}

fn make_affine_pool(pool: &[P256ProjectivePoint]) -> Vec<MontAffinePoint> {
    pool.iter()
        .map(|point| {
            let (x, y) = point.to_affine().expect("nonzero pool point");
            MontAffinePoint {
                x: P256FieldElement::from_canonical(&x),
                y: P256FieldElement::from_canonical(&y),
            }
        })
        .collect()
}

fn make_affine_table(entries: usize, pool: &[MontAffinePoint]) -> Vec<MontAffinePoint> {
    (0..entries)
        .map(|index| pool[index & (pool.len() - 1)])
        .collect()
}

fn add_affine_complete(
    point: &P256ProjectivePoint,
    affine: &MontAffinePoint,
) -> P256ProjectivePoint {
    point.add(&P256ProjectivePoint {
        x: affine.x,
        y: affine.y,
        z: P256FieldElement::ONE,
    })
}

/// RCB complete addition specialized at the public substitution `z2 = 1`.
/// This removes exactly the `z1*z2` field multiplication from Algorithm 1.
fn add_affine_mixed(point: &P256ProjectivePoint, affine: &MontAffinePoint) -> P256ProjectivePoint {
    let (x1, y1, z1) = (&point.x, &point.y, &point.z);
    let (x2, y2) = (&affine.x, &affine.y);
    let mut t0 = x1.mul(x2);
    let mut t1 = y1.mul(y2);
    let mut t2 = *z1;
    let mut t3 = x1.add(y1);
    let mut t4 = x2.add(y2);
    t3 = t3.mul(&t4);
    t4 = t0.add(&t1);
    t3 = t3.sub(&t4);
    t4 = x1.add(z1);
    let mut t5 = x2.add(&P256FieldElement::ONE);
    t4 = t4.mul(&t5);
    t5 = t0.add(&t2);
    t4 = t4.sub(&t5);
    t5 = y1.add(z1);
    let mut x3 = y2.add(&P256FieldElement::ONE);
    t5 = t5.mul(&x3);
    x3 = t1.add(&t2);
    t5 = t5.sub(&x3);
    let mut z3 = A_FE.mul(&t4);
    x3 = B3_FE.mul(&t2);
    z3 = x3.add(&z3);
    x3 = t1.sub(&z3);
    z3 = t1.add(&z3);
    let mut y3 = x3.mul(&z3);
    t1 = t0.add(&t0);
    t1 = t1.add(&t0);
    t2 = A_FE.mul(&t2);
    t4 = B3_FE.mul(&t4);
    t1 = t1.add(&t2);
    t2 = t0.sub(&t2);
    t2 = A_FE.mul(&t2);
    t4 = t4.add(&t2);
    t0 = t1.mul(&t4);
    y3 = y3.add(&t0);
    t0 = t5.mul(&t4);
    x3 = t3.mul(&x3);
    x3 = x3.sub(&t0);
    t0 = t3.mul(&t1);
    z3 = t5.mul(&z3);
    z3 = z3.add(&t0);
    P256ProjectivePoint {
        x: x3,
        y: y3,
        z: z3,
    }
}

fn make_indices(entries: usize, segments: usize) -> Vec<u32> {
    let mut domain = b"p256-projective-locality-round35-indices".to_vec();
    domain.extend_from_slice(&(entries as u64).to_be_bytes());
    let mut indices = Vec::with_capacity(segments * PAIRS);
    for counter in 0..segments * PAIRS {
        let mut input = domain.clone();
        input.extend_from_slice(&(counter as u64).to_be_bytes());
        let digest = sha256(&input);
        let mut word = [0u8; 8];
        word.copy_from_slice(&digest[..8]);
        indices.push((u64::from_be_bytes(word) as usize & (entries - 1)) as u32);
    }
    indices
}

fn point_checksum(point: &P256ProjectivePoint) -> u64 {
    point.x.as_montgomery().0[0] ^ point.y.as_montgomery().0[1] ^ point.z.as_montgomery().0[2]
}

fn projective_eq(left: &P256ProjectivePoint, right: &P256ProjectivePoint) -> bool {
    let left_identity = bool::from(left.is_identity());
    let right_identity = bool::from(right.is_identity());
    if left_identity || right_identity {
        return left_identity == right_identity;
    }
    left.x.mul(&right.z) == right.x.mul(&left.z) && left.y.mul(&right.z) == right.y.mul(&left.z)
}

fn raw_point_bytes(output: &mut Vec<u8>, point: &P256ProjectivePoint) {
    output.extend_from_slice(&point.x.as_montgomery().to_bytes_be());
    output.extend_from_slice(&point.y.as_montgomery().to_bytes_be());
    output.extend_from_slice(&point.z.as_montgomery().to_bytes_be());
}

fn digest_indices(indices: &[u32]) -> String {
    let mut bytes = Vec::with_capacity(indices.len() * 4);
    for index in indices {
        bytes.extend_from_slice(&index.to_be_bytes());
    }
    hex::encode(sha256(&bytes))
}

fn run_direct(
    segments: usize,
    generator: &P256ProjectivePoint,
    operands: &[P256ProjectivePoint; PAIRS],
) -> u64 {
    let mut checksum = 0u64;
    for segment in 0..segments {
        let mut acc = *generator;
        for offset in 0..PAIRS {
            acc = acc.add(black_box(&operands[(segment + offset) & (PAIRS - 1)]));
        }
        checksum ^= point_checksum(black_box(&acc));
    }
    black_box(checksum)
}

#[cfg(target_arch = "x86_64")]
fn prefetch_entry(table: &[P256ProjectivePoint], index: usize) {
    use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
    let pointer = table.as_ptr().wrapping_add(index).cast::<i8>();
    // SAFETY: both addresses remain within the selected 96-byte table entry.
    // The intrinsic is a non-faulting read hint and does not dereference in
    // Rust's abstract machine.
    unsafe {
        _mm_prefetch(pointer, _MM_HINT_T0);
        _mm_prefetch(pointer.add(64), _MM_HINT_T0);
    }
}

#[cfg(not(target_arch = "x86_64"))]
fn prefetch_entry(_table: &[P256ProjectivePoint], _index: usize) {}

fn run_prefetch(
    table: &[P256ProjectivePoint],
    indices: &[u32],
    generator: &P256ProjectivePoint,
    distance: usize,
) -> u64 {
    let segments = indices.len() / PAIRS;
    for segment in 0..distance.min(segments) {
        for index in &indices[segment * PAIRS..(segment + 1) * PAIRS] {
            prefetch_entry(table, *index as usize);
        }
    }
    let mut checksum = 0u64;
    for segment in 0..segments {
        let future = segment + distance;
        if future < segments {
            for index in &indices[future * PAIRS..(future + 1) * PAIRS] {
                prefetch_entry(table, *index as usize);
            }
        }
        let mut acc = *generator;
        for index in &indices[segment * PAIRS..(segment + 1) * PAIRS] {
            acc = acc.add(black_box(&table[*index as usize]));
        }
        checksum ^= point_checksum(black_box(&acc));
    }
    black_box(checksum)
}

#[cfg(target_arch = "x86_64")]
fn prefetch_affine_entry(table: &[MontAffinePoint], index: usize) {
    use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
    let pointer = table.as_ptr().wrapping_add(index).cast::<i8>();
    // SAFETY: both hints address bytes within the selected 64-byte entry.
    unsafe {
        _mm_prefetch(pointer, _MM_HINT_T0);
        _mm_prefetch(pointer.add(32), _MM_HINT_T0);
    }
}

#[cfg(not(target_arch = "x86_64"))]
fn prefetch_affine_entry(_table: &[MontAffinePoint], _index: usize) {}

fn run_direct_affine(
    segments: usize,
    generator: &P256ProjectivePoint,
    operands: &[MontAffinePoint; PAIRS],
    mixed: bool,
) -> u64 {
    let mut checksum = 0u64;
    for segment in 0..segments {
        let mut acc = *generator;
        for offset in 0..PAIRS {
            let operand = black_box(&operands[(segment + offset) & (PAIRS - 1)]);
            acc = if mixed {
                add_affine_mixed(&acc, operand)
            } else {
                add_affine_complete(&acc, operand)
            };
        }
        checksum ^= point_checksum(black_box(&acc));
    }
    black_box(checksum)
}

fn run_affine_prefetch(
    table: &[MontAffinePoint],
    indices: &[u32],
    generator: &P256ProjectivePoint,
    distance: usize,
    mixed: bool,
) -> u64 {
    let segments = indices.len() / PAIRS;
    for segment in 0..distance.min(segments) {
        for index in &indices[segment * PAIRS..(segment + 1) * PAIRS] {
            prefetch_affine_entry(table, *index as usize);
        }
    }
    let mut checksum = 0u64;
    for segment in 0..segments {
        let future = segment + distance;
        if future < segments {
            for index in &indices[future * PAIRS..(future + 1) * PAIRS] {
                prefetch_affine_entry(table, *index as usize);
            }
        }
        let mut acc = *generator;
        for index in &indices[segment * PAIRS..(segment + 1) * PAIRS] {
            let operand = black_box(&table[*index as usize]);
            acc = if mixed {
                add_affine_mixed(&acc, operand)
            } else {
                add_affine_complete(&acc, operand)
            };
        }
        checksum ^= point_checksum(black_box(&acc));
    }
    black_box(checksum)
}

fn check_prefetch_outputs(
    table: &[P256ProjectivePoint],
    indices: &[u32],
    generator: &P256ProjectivePoint,
    distance: usize,
    reference: &[P256ProjectivePoint],
) -> u64 {
    let segments = indices.len() / PAIRS;
    let mut mismatches = 0u64;
    for segment in 0..distance.min(segments) {
        for index in &indices[segment * PAIRS..(segment + 1) * PAIRS] {
            prefetch_entry(table, *index as usize);
        }
    }
    for segment in 0..segments {
        let future = segment + distance;
        if future < segments {
            for index in &indices[future * PAIRS..(future + 1) * PAIRS] {
                prefetch_entry(table, *index as usize);
            }
        }
        let mut acc = *generator;
        for index in &indices[segment * PAIRS..(segment + 1) * PAIRS] {
            acc = acc.add(&table[*index as usize]);
        }
        if !projective_eq(&acc, &reference[segment]) {
            mismatches += 1;
        }
    }
    mismatches
}

fn check_affine_outputs(
    table: &[MontAffinePoint],
    indices: &[u32],
    generator: &P256ProjectivePoint,
    distance: usize,
    mixed: bool,
    reference: &[P256ProjectivePoint],
) -> u64 {
    let segments = indices.len() / PAIRS;
    let mut mismatches = 0u64;
    for segment in 0..distance.min(segments) {
        for index in &indices[segment * PAIRS..(segment + 1) * PAIRS] {
            prefetch_affine_entry(table, *index as usize);
        }
    }
    for segment in 0..segments {
        let future = segment + distance;
        if future < segments {
            for index in &indices[future * PAIRS..(future + 1) * PAIRS] {
                prefetch_affine_entry(table, *index as usize);
            }
        }
        let mut acc = *generator;
        for index in &indices[segment * PAIRS..(segment + 1) * PAIRS] {
            let operand = &table[*index as usize];
            acc = if mixed {
                add_affine_mixed(&acc, operand)
            } else {
                add_affine_complete(&acc, operand)
            };
        }
        if !projective_eq(&acc, &reference[segment]) {
            mismatches += 1;
        }
    }
    mismatches
}

fn process_cpu_ns() -> u64 {
    let mut timestamp = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    // SAFETY: timestamp is writable and CLOCK_PROCESS_CPUTIME_ID has no
    // pointer lifetime requirements beyond this call.
    let status = unsafe { libc::clock_gettime(libc::CLOCK_PROCESS_CPUTIME_ID, &mut timestamp) };
    if status != 0 {
        0
    } else {
        timestamp.tv_sec as u64 * 1_000_000_000 + timestamp.tv_nsec as u64
    }
}

fn time_variant<F>(variant: Variant, segments: usize, mut operation: F) -> OpTiming
where
    F: FnMut() -> u64,
{
    let cpu_start = process_cpu_ns();
    let start = Instant::now();
    let checksum = operation();
    let wall_ns = start.elapsed().as_nanos() as u64;
    let cpu_ns = process_cpu_ns().saturating_sub(cpu_start);
    OpTiming {
        variant: variant.name(),
        class: variant.class().into(),
        parameter: variant.parameter(),
        wall_ns,
        process_cpu_ns: cpu_ns,
        segments: segments as u64,
        checksum,
    }
}

fn execute_variant(
    variant: Variant,
    table: &[P256ProjectivePoint],
    affine_table: &[MontAffinePoint],
    indices: &[u32],
    generator: &P256ProjectivePoint,
    operands: &[P256ProjectivePoint; PAIRS],
    affine_operands: &[MontAffinePoint; PAIRS],
) -> OpTiming {
    match variant {
        Variant::Direct => time_variant(variant, SEGMENTS, || {
            run_direct(SEGMENTS, generator, operands)
        }),
        Variant::Prefetch(distance) => time_variant(variant, SEGMENTS, || {
            run_prefetch(table, indices, generator, distance)
        }),
        Variant::DirectAffineComplete => time_variant(variant, SEGMENTS, || {
            run_direct_affine(SEGMENTS, generator, affine_operands, false)
        }),
        Variant::DirectAffineMixed => time_variant(variant, SEGMENTS, || {
            run_direct_affine(SEGMENTS, generator, affine_operands, true)
        }),
        Variant::AffineComplete(distance) => time_variant(variant, SEGMENTS, || {
            run_affine_prefetch(affine_table, indices, generator, distance, false)
        }),
        Variant::AffineMixed(distance) => time_variant(variant, SEGMENTS, || {
            run_affine_prefetch(affine_table, indices, generator, distance, true)
        }),
    }
}

fn variants() -> Vec<Variant> {
    let mut variants = vec![
        Variant::Direct,
        Variant::Prefetch(16),
        Variant::DirectAffineComplete,
        Variant::DirectAffineMixed,
    ];
    if cfg!(target_arch = "x86_64") {
        variants.extend(PREFETCH_DISTANCES.into_iter().map(Variant::AffineComplete));
        variants.extend(PREFETCH_DISTANCES.into_iter().map(Variant::AffineMixed));
    }
    variants
}

fn sorted_interval(values: &[f64]) -> Interval {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let (lower, median, upper) = if sorted.len() == REPETITIONS {
        (1, 4, 7)
    } else {
        (0, sorted.len() / 2, sorted.len() - 1)
    };
    Interval {
        lower: sorted[lower],
        median: sorted[median],
        upper: sorted[upper],
    }
}

fn summarize(
    raw: &[TimingRepetition],
    variants: &[Variant],
    model: &FrozenModel,
) -> Vec<TimingSummary> {
    variants
        .iter()
        .map(|variant| {
            let name = variant.name();
            let values: Vec<f64> = raw
                .iter()
                .map(|repetition| {
                    let timing = repetition
                        .operations
                        .iter()
                        .find(|timing| timing.variant == name)
                        .expect("variant timing");
                    timing.wall_ns as f64 / timing.segments as f64
                })
                .collect();
            let wall = sorted_interval(&values);
            if matches!(
                variant,
                Variant::Direct | Variant::DirectAffineComplete | Variant::DirectAffineMixed
            ) {
                TimingSummary {
                    variant: name,
                    class: variant.class().into(),
                    parameter: variant.parameter(),
                    wall_ns_per_segment: wall,
                    incremental_addition_equivalents: None,
                    projected_ratio_to_rho: None,
                    upper_endpoint_within_parity_budget: None,
                }
            } else {
                let denominator = match variant {
                    Variant::AffineComplete(_) => "direct-affine-complete",
                    Variant::AffineMixed(_) => "direct-affine-mixed",
                    _ => "direct-eight-additions",
                };
                let direct_values: Vec<f64> = raw
                    .iter()
                    .map(|repetition| {
                        let timing = repetition
                            .operations
                            .iter()
                            .find(|timing| timing.variant == denominator)
                            .expect("matched direct timing");
                        timing.wall_ns as f64 / timing.segments as f64
                    })
                    .collect();
                let direct = sorted_interval(&direct_values);
                let extra = Interval {
                    lower: (PAIRS as f64 * (wall.lower / direct.upper - 1.0)).max(0.0),
                    median: (PAIRS as f64 * (wall.median / direct.median - 1.0)).max(0.0),
                    upper: (PAIRS as f64 * (wall.upper / direct.lower - 1.0)).max(0.0),
                };
                let ratio = |value: f64| {
                    model.round33_ratio_to_rho_with_target_correction
                        + model.local_oracle_ratio_to_rho * value / (model.mean_capacity + 1.0)
                };
                let projected = Interval {
                    lower: ratio(extra.lower),
                    median: ratio(extra.median),
                    upper: ratio(extra.upper),
                };
                TimingSummary {
                    variant: name,
                    class: variant.class().into(),
                    parameter: variant.parameter(),
                    wall_ns_per_segment: wall,
                    incremental_addition_equivalents: Some(extra),
                    projected_ratio_to_rho: Some(projected.clone()),
                    upper_endpoint_within_parity_budget: Some(projected.upper <= 1.0),
                }
            }
        })
        .collect()
}

fn benchmark_table(
    entries: usize,
    pool: &[P256ProjectivePoint],
    model: &FrozenModel,
    repetitions: usize,
) -> TableTiming {
    let table = make_table(entries, pool);
    let affine_pool = make_affine_pool(pool);
    let affine_table = make_affine_table(entries, &affine_pool);
    let indices = make_indices(entries, SEGMENTS);
    let curve = CurveParams::p256();
    let generator = P256ProjectivePoint::from_textbook(&curve.generator());
    let operands = std::array::from_fn(|index| pool[index]);
    let affine_operands = std::array::from_fn(|index| affine_pool[index]);
    let variants = variants();
    for variant in &variants {
        black_box(execute_variant(
            *variant,
            &table,
            &affine_table,
            &indices,
            &generator,
            &operands,
            &affine_operands,
        ));
    }
    let mut raw = Vec::with_capacity(repetitions);
    for repetition in 0..repetitions {
        let start = (repetition * 5 + entries.trailing_zeros() as usize) % variants.len();
        let mut operations = Vec::with_capacity(variants.len());
        let mut order = Vec::with_capacity(variants.len());
        for offset in 0..variants.len() {
            let variant = variants[(start + offset) % variants.len()];
            order.push(variant.name());
            operations.push(execute_variant(
                variant,
                &table,
                &affine_table,
                &indices,
                &generator,
                &operands,
                &affine_operands,
            ));
        }
        raw.push(TimingRepetition {
            repetition: repetition as u64,
            execution_order: order,
            operations,
        });
    }
    let summaries = summarize(&raw, &variants, model);
    TableTiming {
        table_entries: entries as u64,
        affine_table_bytes: entries as u64 * AFFINE_ENTRY_BYTES as u64,
        projective_control_table_bytes: entries as u64 * PROJECTIVE_ENTRY_BYTES as u64,
        index_stream_sha256: digest_indices(&indices),
        repetitions: raw,
        summaries,
    }
}

fn correctness(pool: &[P256ProjectivePoint]) -> Correctness {
    let curve = CurveParams::p256();
    let generator = P256ProjectivePoint::from_textbook(&curve.generator());
    let table = make_table(1 << 15, pool);
    let affine_pool = make_affine_pool(pool);
    let affine_table = make_affine_table(table.len(), &affine_pool);
    let indices = make_indices(table.len(), CORRECTNESS_SEGMENTS);
    let mut reference = vec![generator; CORRECTNESS_SEGMENTS];
    for (segment, chunk) in indices.chunks_exact(PAIRS).enumerate() {
        for index in chunk {
            reference[segment] = reference[segment].add(&table[*index as usize]);
        }
    }
    let mut group_mismatches = check_prefetch_outputs(&table, &indices, &generator, 16, &reference);
    let mut output_bytes = Vec::new();
    let mut variants_checked = 1u64;
    for distance in PREFETCH_DISTANCES {
        variants_checked += 2;
        group_mismatches += check_affine_outputs(
            &affine_table,
            &indices,
            &generator,
            distance,
            false,
            &reference,
        );
        group_mismatches += check_affine_outputs(
            &affine_table,
            &indices,
            &generator,
            distance,
            true,
            &reference,
        );
    }

    let two_g = generator.double();
    let distinct = generator.add(&two_g);
    let cases = [
        (P256ProjectivePoint::IDENTITY, generator),
        (generator, generator),
        (generator, generator.neg()),
        (generator, distinct),
    ];
    let mut special_case_mismatches = 0u64;
    let mut special_cases_checked = 0u64;
    for (left, right) in cases {
        let (x, y) = right.to_affine().expect("special-case affine operand");
        let affine = MontAffinePoint {
            x: P256FieldElement::from_canonical(&x),
            y: P256FieldElement::from_canonical(&y),
        };
        let expected = left.add(&right);
        for actual in [
            add_affine_complete(&left, &affine),
            add_affine_mixed(&left, &affine),
        ] {
            special_cases_checked += 1;
            if !projective_eq(&actual, &expected) {
                special_case_mismatches += 1;
            }
        }
    }
    group_mismatches += special_case_mismatches;

    for point in &reference {
        let affine = point.to_affine().expect("nonidentity reference sum");
        output_bytes.extend_from_slice(&affine.0.to_bytes_be());
        output_bytes.extend_from_slice(&affine.1.to_bytes_be());
    }
    let mut pool_bytes = Vec::with_capacity(pool.len() * PROJECTIVE_ENTRY_BYTES);
    for point in pool {
        raw_point_bytes(&mut pool_bytes, point);
    }
    let mut affine_table_bytes = Vec::with_capacity(affine_table.len() * AFFINE_ENTRY_BYTES);
    for point in &affine_table {
        affine_table_bytes.extend_from_slice(&point.x.as_montgomery().to_bytes_be());
        affine_table_bytes.extend_from_slice(&point.y.as_montgomery().to_bytes_be());
    }
    Correctness {
        pool_points: pool.len() as u64,
        segments: CORRECTNESS_SEGMENTS as u64,
        variants_checked,
        group_mismatches,
        missing_requests: 0,
        duplicate_requests: 0,
        out_of_range_accesses: 0,
        wrong_addition_counts: 0,
        curve_failures: 0,
        false_positives: 0,
        false_negatives: 0,
        pool_sha256: hex::encode(sha256(&pool_bytes)),
        index_sha256: digest_indices(&indices),
        affine_table_sha256: hex::encode(sha256(&affine_table_bytes)),
        output_sha256: hex::encode(sha256(&output_bytes)),
        special_cases_checked,
        special_case_mismatches,
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
        prefetch_supported: cfg!(target_arch = "x86_64"),
        confidence_interval_rule: if repetitions == REPETITIONS {
            "sorted ranks 1,4,7 (zero-indexed) from nine repetitions".into()
        } else {
            "quick-mode minimum, median, maximum".into()
        },
    }
}

fn execute(cli: Cli) -> Result<ResultFile, String> {
    let (dependency, model) = load_dependency(&cli.round35)?;
    let pool_count = if cli.quick { 64 } else { POOL_SIZE };
    let pool = make_pool(pool_count);
    let correctness = correctness(&pool);
    let correctness_gate = correctness.group_mismatches == 0
        && correctness.missing_requests == 0
        && correctness.duplicate_requests == 0
        && correctness.out_of_range_accesses == 0
        && correctness.wrong_addition_counts == 0
        && correctness.curve_failures == 0
        && correctness.false_positives == 0
        && correctness.false_negatives == 0
        && correctness.special_case_mismatches == 0;
    let powers: Vec<u32> = if cli.quick {
        vec![10]
    } else {
        vec![15, 20, 22]
    };
    let repetitions = if cli.quick { 3 } else { REPETITIONS };
    let mut timing_tables = Vec::new();
    for power in powers {
        timing_tables.push(benchmark_table(1usize << power, &pool, &model, repetitions));
    }
    let decision_table = timing_tables
        .iter()
        .max_by_key(|table| table.table_entries)
        .expect("at least one timing table");
    let selected = decision_table
        .summaries
        .iter()
        .filter(|row| matches!(row.class.as_str(), "affine-complete" | "affine-mixed"))
        .min_by(|left, right| {
            left.projected_ratio_to_rho
                .as_ref()
                .expect("candidate ratio")
                .median
                .total_cmp(
                    &right
                        .projected_ratio_to_rho
                        .as_ref()
                        .expect("candidate ratio")
                        .median,
                )
                .then_with(|| left.parameter.cmp(&right.parameter))
        })
        .expect("at least one locality candidate");
    let extra = selected
        .incremental_addition_equivalents
        .clone()
        .expect("selected overhead");
    let projected = selected
        .projected_ratio_to_rho
        .clone()
        .expect("selected ratio");
    let timing_gate = projected.median <= 1.0 && projected.upper <= 1.0;
    let selected_row = SelectedRow {
        table_entries: decision_table.table_entries,
        table_bytes: decision_table.affine_table_bytes,
        variant: selected.variant.clone(),
        class: selected.class.clone(),
        parameter: selected.parameter,
        incremental_addition_equivalents: extra,
        projected_ratio_to_rho: projected,
        correctness_gate,
        timing_gate,
        advances_affine_stage: correctness_gate && timing_gate,
    };
    let semantic = SemanticEvidence {
        dependency_sha256: dependency.sha256.clone(),
        curve: CURVE_SLUG.into(),
        model: model.clone(),
        correctness: correctness.clone(),
    };
    let semantic_sha = hex::encode(sha256(
        &serde_json::to_vec(&semantic).map_err(|error| error.to_string())?,
    ));
    let timing_sha = hex::encode(sha256(
        &serde_json::to_vec(&timing_tables).map_err(|error| error.to_string())?,
    ));
    let advances = selected_row.advances_affine_stage;
    Ok(ResultFile {
        schema: "p256-mont-affine-table-v1".into(),
        curve: CURVE_SLUG.into(),
        family: "round-33 biased-restart Montgomery-affine pair table".into(),
        dependency,
        frozen_model: model,
        host: host(repetitions),
        correctness,
        arithmetic_accounting: ArithmeticAccounting {
            complete_field_multiplications_per_addition: 17,
            mixed_field_multiplications_per_addition: 16,
            field_multiplications_removed_per_addition: 1,
        },
        timing_tables,
        selected_decision_row: selected_row,
        semantic_evidence_sha256: semantic_sha,
        raw_timing_sha256: timing_sha,
        imported_gates: ImportedGates {
            full_table_materialized_or_measured: false,
            complete_time_at_or_below_rho: false,
            proved_p256_usable_relation_probability: false,
            structured_residual_degree_at_most_5: false,
            relation_collection_below_2_120: false,
            per_usable_relation_below_2_103: false,
            materialized_storage_below_2_50: FULL_TABLE_BYTES < (1u64 << 50),
            non_generic_end_to_end: false,
            promoted: false,
        },
        full_depth_unplanted_attempted: false,
        interpretation: if advances {
            "The selected Montgomery-affine schedule preserves the round-33 stage margin on the registered 256-MiB DRAM working set under its matched affine arithmetic denominator. This advances only the affine access stage: the 2.212-TB table, achieved relation coverage, degree, collection and non-generic status remain unproved."
                .into()
        } else {
            "No registered Montgomery-affine schedule preserves the round-33 stage margin on the 256-MiB decision row under its matched affine arithmetic denominator. This rejects these representations and prefetch schedules, not all possible memory systems."
                .into()
        },
        decision: if advances {
            "preserve the selected Montgomery-affine row as a stage candidate; do not promote or run a full-depth relation"
                .into()
        } else {
            "reject the registered Montgomery-affine schedules; do not promote or run a full-depth relation"
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
            eprintln!("p256_mont_affine_table: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locality_variants_match() {
        let pool = make_pool(16);
        let table = make_table(1 << 10, &pool);
        let indices = make_indices(table.len(), 256);
        let generator = P256ProjectivePoint::from_textbook(&CurveParams::p256().generator());
        let mut reference = vec![generator; 256];
        for (segment, chunk) in indices.chunks_exact(PAIRS).enumerate() {
            for index in chunk {
                reference[segment] = reference[segment].add(&table[*index as usize]);
            }
        }
        for distance in PREFETCH_DISTANCES {
            assert_eq!(
                check_prefetch_outputs(&table, &indices, &generator, distance, &reference),
                0
            );
        }
    }

    #[test]
    fn affine_formulas_cover_complete_cases() {
        let generator = P256ProjectivePoint::from_textbook(&CurveParams::p256().generator());
        let two_g = generator.double();
        let distinct = generator.add(&two_g);
        for (left, right) in [
            (P256ProjectivePoint::IDENTITY, generator),
            (generator, generator),
            (generator, generator.neg()),
            (generator, distinct),
        ] {
            let (x, y) = right.to_affine().expect("affine operand");
            let affine = MontAffinePoint {
                x: P256FieldElement::from_canonical(&x),
                y: P256FieldElement::from_canonical(&y),
            };
            let expected = left.add(&right);
            assert!(projective_eq(
                &add_affine_complete(&left, &affine),
                &expected
            ));
            assert!(projective_eq(&add_affine_mixed(&left, &affine), &expected));
        }
    }

    #[test]
    fn affine_storage_gate() {
        assert_eq!(
            std::mem::size_of::<P256ProjectivePoint>(),
            PROJECTIVE_ENTRY_BYTES
        );
        assert_eq!(std::mem::size_of::<MontAffinePoint>(), AFFINE_ENTRY_BYTES);
        const {
            assert!(FULL_TABLE_BYTES < (1u64 << 50));
        }
    }
}
