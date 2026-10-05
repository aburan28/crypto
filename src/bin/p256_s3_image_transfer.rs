//! Transfer the indexed final-S3 image gate to the committed P-256 factor base.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{self, CURVE_SLUG};
use crypto_lib::ecc::p256_field::P256FieldElement as Fe;
use crypto_lib::ecc::{CurveParams, Point};
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

const SPEC: &str = "dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129";
const FB_ID: &str = "FB1h2f8621cda105";
const FB_SHA256: &str = "2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42";
const POINTS_SHA256: &str = "70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1";
const TERMINAL: &str = "0x5b17195299a3158b93389ad04c776fff2a8bb23ca8659b5b0c2b75f9009b65e5";
const COLUMNS: u64 = 131_458;
const SIGNED_POINTS: u64 = 262_916;
const ROUND13_SHA256: &str = "2008dcb659d3120157480b6096a4873d1f9c23ee30123f4dd3a32d450f53ba1d";
const ROUND14_SHA256: &str = "6eefb6b27768023af5850cd785d75ef3729484ac85e5defef9d1d596834ebccc";
const ROUND15_SHA256: &str = "4220dafa066613630338886a916218602cd61ca914bd66ac7994cde163a53ad5";
const ROUND16_SHA256: &str = "af2874ef5f0bdabec75b9bf707067232d0d6e3985f0cce4eba6834c2b5790e44";
const ROUND17_SHA256: &str = "57cf49ffbab65356b4ef69848e7824f0377bc058ddd57fe0bf6a8f7b77089ce7";
const ROUND6_SHA256: &str = "71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4";
const ATOM_COLUMNS: usize = 16;
const COLUMN_INDEX_BITS: usize = 18;
const PACKED_ATOM_BYTES: usize = ATOM_COLUMNS * COLUMN_INDEX_BITS / 8;
const OUTER_OFFSET: u64 = 90_322;
const OUTER_STRIDE: u64 = 23_509;
const OUTER_MAX_DEPTH: u32 = 10;
const OUTER_DEPTHS: [u32; 7] = [4, 5, 6, 7, 8, 9, 10];
const PLANTED_OUTER_POSITION: u64 = 7;
const PLANTED_OUTER_COLUMN: u64 = 123_427;
const OUTER_PUBLIC_TARGET_PREIMAGE: &str = concat!(
    "icv1-fp256-t89188191154553853111372247798585809583-f188c491",
    "/s17-outer-scan-round18/public-target-0"
);
const TARGET_PREIMAGE: &str = concat!(
    "icv1-fp256-t89188191154553853111372247798585809583-f188c491",
    "/s3-image-transfer-round8/target/0"
);

#[derive(Parser)]
#[command(about = "Validate the indexed final-S3 gate on the selected P-256 factor base")]
struct Cli {
    /// Compact deterministic JSON output; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
    /// Run the round-14 P-256 sixteen-leaf width transfer after hash-checking round 13.
    #[arg(long)]
    width_round13: Option<PathBuf>,
    /// Run the round-15 target-indexed compression after hash-checking round 14.
    #[arg(long)]
    compress_round14: Option<PathBuf>,
    /// Run the round-16 packed-atom compression after hash-checking round 15.
    #[arg(long)]
    pack_round15: Option<PathBuf>,
    /// Run the round-17 batched-inversion candidate after hash-checking round 16.
    #[arg(long)]
    batch_round16: Option<PathBuf>,
    /// Run the round-18 charged outer scan after hash-checking round 17.
    #[arg(long)]
    outer_scan_round17: Option<PathBuf>,
    /// Round-6 residual-degree receipt required by --outer-scan-round17.
    #[arg(long)]
    round6: Option<PathBuf>,
}

#[derive(Clone)]
struct Column {
    x: Fe,
    low: Point,
}

#[derive(Clone)]
struct SignedRow {
    point: Point,
    col: u64,
}

#[derive(Clone, Copy, Default, Serialize)]
struct MultiplicationCounts {
    coefficients: u64,
    discriminants: u64,
    square_roots: u64,
    inversions: u64,
    root_construction: u64,
    total: u64,
}

impl MultiplicationCounts {
    fn finish(&mut self) {
        self.total = self.coefficients
            + self.discriminants
            + self.square_roots
            + self.inversions
            + self.root_construction;
    }

    fn combined(self, other: Self) -> Self {
        let mut combined = Self {
            coefficients: self.coefficients + other.coefficients,
            discriminants: self.discriminants + other.discriminants,
            square_roots: self.square_roots + other.square_roots,
            inversions: self.inversions + other.inversions,
            root_construction: self.root_construction + other.root_construction,
            total: 0,
        };
        combined.finish();
        combined
    }

    fn add_assign(&mut self, other: Self) {
        self.coefficients += other.coefficients;
        self.discriminants += other.discriminants;
        self.square_roots += other.square_roots;
        self.inversions += other.inversions;
        self.root_construction += other.root_construction;
        self.finish();
    }
}

#[derive(Serialize)]
struct FactorBaseReceipt {
    spec: String,
    fb_id: String,
    fb_sha256: String,
    points_sha256: String,
    terminal: String,
    columns: u64,
    signed_points: u64,
    complete_rebuild_verified: bool,
}

#[derive(Serialize)]
struct TargetResult {
    kind: String,
    scalar: Option<String>,
    target: [String; 2],
    rejected_ordered_column_pairs: u64,
    quadratic_solves: u64,
    quadratic_roots_returned: u64,
    algebra_index_lookups: u64,
    linear_degeneracies: u64,
    universal_degeneracies: u64,
    algebra_pairs: usize,
    reference_pairs: usize,
    false_negatives: usize,
    false_positives: usize,
    multiplication_counts: MultiplicationCounts,
    reference_group_subtractions: u64,
    reference_index_lookups: u64,
    verification_group_additions: u64,
    all_emitted_pairs_verified: bool,
    exact: bool,
    pair_sha256: String,
    pairs: Vec<[u64; 2]>,
}

#[derive(Serialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    factor_base: FactorBaseReceipt,
    target_preimage: String,
    targets: Vec<TargetResult>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct P256ImageState {
    affine: BTreeMap<[u8; 32], Fe>,
    identity: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReferenceImageState {
    affine: BTreeSet<[u8; 32]>,
    identity: bool,
}

#[derive(Serialize)]
struct WidthSampleResult {
    kind: String,
    columns: Vec<u64>,
    x_coordinates: Vec<String>,
    two_leaf_widths: Vec<usize>,
    four_leaf_widths: Vec<usize>,
    eight_leaf_widths: Vec<usize>,
    sixteen_leaf_width: usize,
    identity_images_by_level: [u64; 4],
    quadratic_solves: u64,
    quadratic_roots_returned: u64,
    linear_degeneracies: u64,
    universal_degeneracies: u64,
    multiplication_counts: MultiplicationCounts,
    reference_group_additions: u64,
    exact: bool,
    final_image_sha256: String,
}

#[derive(Serialize)]
struct WidthStorageModel {
    generic_sixteen_leaf_x_entries: u64,
    bytes_per_x: u64,
    generic_sixteen_leaf_raw_bytes: u64,
    factor_base_columns: u64,
    unordered_factor_base_pairs: u64,
    two_root_pair_image_entry_upper_bound: u64,
    two_root_pair_image_raw_byte_upper_bound: u64,
    identity_aware_pair_image_affine_entries: u64,
    identity_aware_pair_image_raw_bytes: u64,
}

#[derive(Serialize)]
struct WidthExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    round13_sha256: String,
    factor_base: FactorBaseReceipt,
    local_maximum_degree: u32,
    generic_widths: [u64; 4],
    samples: Vec<WidthSampleResult>,
    storage_model: WidthStorageModel,
}

#[derive(Serialize)]
struct CompressedTargetResult {
    kind: String,
    scalar: Option<String>,
    target: [String; 2],
    reference_positive: bool,
    candidate_positive: bool,
    algebra_hits: u64,
    quadratic_solves_warm: u64,
    quadratic_solves_cold: u64,
    quadratic_roots_returned: u64,
    algebra_index_lookups: u64,
    linear_degeneracies: u64,
    universal_degeneracies: u64,
    multiplication_counts_warm: MultiplicationCounts,
    multiplication_counts_cold: MultiplicationCounts,
    false_negative: bool,
    false_positive: bool,
    exact: bool,
    hit_sha256: String,
}

#[derive(Serialize)]
struct CompressionSampleResult {
    kind: String,
    columns: Vec<u64>,
    x_coordinates: Vec<String>,
    left_eight_leaf_entries: usize,
    right_eight_leaf_entries: usize,
    retained_entries: usize,
    materialized_sixteen_leaf_entries: usize,
    retained_raw_bytes: usize,
    materialized_raw_bytes: usize,
    retained_to_materialized_ratio: f64,
    build_quadratic_solves: u64,
    full_materialization_quadratic_solves: u64,
    intermediate_reference_group_additions: u64,
    full_reference_group_additions: u64,
    build_linear_degeneracies: u64,
    build_universal_degeneracies: u64,
    build_multiplication_counts: MultiplicationCounts,
    left_image_sha256: String,
    right_image_sha256: String,
    targets: Vec<CompressedTargetResult>,
}

#[derive(Serialize)]
struct CompressionExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    round14_sha256: String,
    factor_base: FactorBaseReceipt,
    local_maximum_degree: u32,
    samples: Vec<CompressionSampleResult>,
}

#[derive(Serialize)]
struct PackedAtomSampleResult {
    kind: String,
    packed_atom_hex: String,
    packed_atom_sha256: String,
    bits_per_column: usize,
    packet_bytes: usize,
    round15_retained_raw_bytes: usize,
    materialized_sixteen_leaf_raw_bytes: usize,
    additional_persistent_reduction: f64,
    full_materialized_reduction: f64,
    columns: Vec<u64>,
    x_coordinates: Vec<String>,
    decoded_exact: bool,
    reencoded_exact: bool,
    dependency_sample_exact: bool,
    left_image_sha256: String,
    right_image_sha256: String,
    transient_reconstruction_raw_bytes: usize,
    build_quadratic_solves: u64,
    one_cold_query_quadratic_solves: u64,
    two_query_one_rebuild_quadratic_solves: u64,
    intermediate_reference_group_additions: u64,
    full_reference_group_additions: u64,
    build_linear_degeneracies: u64,
    build_universal_degeneracies: u64,
    build_multiplication_counts: MultiplicationCounts,
    targets: Vec<CompressedTargetResult>,
}

#[derive(Serialize)]
struct PackedAtomExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    round15_sha256: String,
    factor_base: FactorBaseReceipt,
    local_maximum_degree: u32,
    factor_base_column_bits: usize,
    packed_atom_bytes: usize,
    samples: Vec<PackedAtomSampleResult>,
}

#[derive(Clone, Default, Serialize)]
struct BatchInversionProfile {
    batch_sizes: Vec<usize>,
    batch_inversions: u64,
    scalar_fallbacks: u64,
}

#[derive(Serialize)]
struct BatchedTargetResult {
    kind: String,
    scalar: Option<String>,
    target: [String; 2],
    reference_positive: bool,
    scalar_positive: bool,
    batched_positive: bool,
    algebra_hits: u64,
    quadratic_solves: u64,
    quadratic_roots_returned: u64,
    algebra_index_lookups: u64,
    linear_degeneracies: u64,
    universal_degeneracies: u64,
    scalar_multiplication_counts_warm: MultiplicationCounts,
    batched_multiplication_counts_warm: MultiplicationCounts,
    scalar_multiplication_counts_cold: MultiplicationCounts,
    batched_multiplication_counts_cold: MultiplicationCounts,
    warm_multiplication_reduction: f64,
    cold_multiplication_reduction: f64,
    inversion_profile: BatchInversionProfile,
    false_negative: bool,
    false_positive: bool,
    exact: bool,
    hit_sha256: String,
}

#[derive(Serialize)]
struct BatchedInversionSampleResult {
    kind: String,
    packed_atom_hex: String,
    columns: Vec<u64>,
    decoded_exact: bool,
    left_image_sha256: String,
    right_image_sha256: String,
    scalar_build_multiplication_counts: MultiplicationCounts,
    batched_build_multiplication_counts: MultiplicationCounts,
    build_multiplication_reduction: f64,
    scalar_build_quadratic_solves: u64,
    batched_build_quadratic_solves: u64,
    build_inversion_profile: BatchInversionProfile,
    two_target_scalar_multiplications: u64,
    two_target_batched_multiplications: u64,
    two_target_multiplication_reduction: f64,
    transient_reconstruction_raw_bytes: usize,
    intermediate_reference_group_additions: u64,
    full_reference_group_additions: u64,
    targets: Vec<BatchedTargetResult>,
    exact: bool,
}

#[derive(Serialize)]
struct BatchedInversionExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    round16_sha256: String,
    factor_base: FactorBaseReceipt,
    local_maximum_degree: u32,
    modeled_relation_success_unchanged: f64,
    end_to_end_fraction_in_oracle: Option<f64>,
    end_to_end_speedup: Option<f64>,
    samples: Vec<BatchedInversionSampleResult>,
}

#[derive(Clone, Default, Serialize)]
struct BatchSizeTotals {
    sizes: BTreeMap<usize, u64>,
    batch_inversions: u64,
    scalar_fallbacks: u64,
}

impl BatchSizeTotals {
    fn record(&mut self, profile: &BatchInversionProfile) {
        for &size in &profile.batch_sizes {
            *self.sizes.entry(size).or_default() += 1;
        }
        self.batch_inversions += profile.batch_inversions;
        self.scalar_fallbacks += profile.scalar_fallbacks;
    }
}

struct BatchedQueryObservation {
    positive: bool,
    algebra_hits: u64,
    quadratic_solves: u64,
    quadratic_roots_returned: u64,
    algebra_index_lookups: u64,
    linear_degeneracies: u64,
    universal_degeneracies: u64,
    multiplication_counts: MultiplicationCounts,
    inversion_profile: BatchInversionProfile,
    hit_sha256: String,
}

#[derive(Clone, Serialize)]
struct OuterTargetIdentity {
    kind: String,
    scalar: Option<String>,
    target: [String; 2],
    planted: bool,
}

#[derive(Clone, Serialize)]
struct OuterRelation {
    target_kind: String,
    stream_position: u64,
    outer_column: u64,
    outer_negative: bool,
    atom_negative_mask: u32,
    direct_verified: bool,
}

#[derive(Clone, Serialize)]
struct OuterTargetCheckpoint {
    kind: String,
    signed_branch_cells: u64,
    oracle_calls: u64,
    quadratic_solves: u64,
    quadratic_roots_returned: u64,
    algebra_index_lookups: u64,
    algebra_hits: u64,
    candidate_positive_trials: u64,
    reference_positive_trials: u64,
    false_negatives: u64,
    false_positives: u64,
    multiplication_counts: MultiplicationCounts,
    target_adjustment_group_additions: u64,
    reference_index_lookups: u64,
    direct_relation_verification_additions: u64,
    batch_profile: BatchSizeTotals,
    relation_count: usize,
    relation_sha256: String,
    trial_stream_sha256: String,
    exact: bool,
}

#[derive(Serialize)]
struct OuterCheckpoint {
    depth: u32,
    stream_positions: u64,
    skipped_atom_columns: u64,
    eligible_outer_columns: u64,
    targets: Vec<OuterTargetCheckpoint>,
}

#[derive(Serialize)]
struct OuterGrowthFit {
    depths: Vec<u32>,
    public_field_multiplications: Vec<u64>,
    log2_multiplication_slope_per_depth: f64,
    outer_calls_growth_per_depth: f64,
}

#[derive(Serialize)]
struct ResidualDegreeEvidence {
    round6_sha256: String,
    maximum_by_residual_depth: BTreeMap<u32, u32>,
    every_component_complete_and_correct: bool,
    promotion_degree_gate_at_most_five: bool,
}

#[derive(Serialize)]
struct OuterMemoryModel {
    persistent_atom_bytes: u64,
    reconstructed_image_raw_bytes: u64,
    candidate_peak_logical_bytes: u64,
    factor_base_point_raw_bytes: u64,
    reference_entries: u64,
    reference_final_logical_bytes: u64,
    reference_build_peak_logical_bytes: u64,
    reference_memory_excluded_from_candidate: bool,
}

#[derive(Serialize)]
struct OuterCollectionProjection {
    classification: String,
    distinct_sixteen_column_atoms: String,
    distinct_sixteen_column_atoms_log2: f64,
    distinct_seventeen_column_sets: String,
    distinct_seventeen_column_sets_log2: f64,
    eligible_outer_columns_per_atom: u64,
    signed_oracle_calls_per_complete_atom_scan: u64,
    signed_candidate_domain_per_atom: String,
    poisson_mean_per_atom: f64,
    poisson_success_per_atom: f64,
    expected_atom_scans_per_relation: f64,
    expected_atom_scans_log2: f64,
    measured_public_multiplications_per_oracle_call: f64,
    projected_complete_atom_scan_multiplications: f64,
    projected_complete_atom_scan_multiplications_log2: f64,
    projected_multiplications_per_relation: f64,
    projected_multiplications_per_relation_log2: f64,
    whole_factor_base_signed_domain: String,
    whole_factor_base_poisson_mean: f64,
    whole_factor_base_poisson_success: f64,
    relation_rows: u64,
    projected_targets_for_relation_rows: f64,
    projected_collection_multiplications: f64,
    projected_collection_multiplications_log2: f64,
    below_2_pow_120: bool,
    below_2_pow_128: bool,
    bits_above_2_pow_128: f64,
    optimistic_lower_projection: bool,
}

#[derive(Serialize)]
struct SparseLinearAlgebraProjection {
    rows: u64,
    columns: u64,
    row_weight: u64,
    nonzeros: u64,
    csr_entry_bytes: u64,
    csr_row_offset_bytes: u64,
    csr_total_bytes: u64,
    wiedemann_sparse_matvecs: u64,
    wiedemann_nonzero_additions: u64,
    berlekamp_massey_scalar_operations: u64,
    field_vectors: u64,
    field_vector_bytes: u64,
    modeled_working_bytes: u64,
    operation_unit_separate_from_collection: bool,
}

#[derive(Serialize)]
struct OuterScanExperimentResult {
    schema: String,
    curve: String,
    field_prime: String,
    curve_a: String,
    curve_b: String,
    round17_sha256: String,
    factor_base: FactorBaseReceipt,
    atom_kind: String,
    packed_atom_hex: String,
    atom_columns: Vec<u64>,
    outer_offset: u64,
    outer_stride: u64,
    checkpoint_depths: Vec<u32>,
    planted_outer_position: u64,
    planted_outer_column: u64,
    targets: Vec<OuterTargetIdentity>,
    scalar_build_multiplication_counts: MultiplicationCounts,
    batched_build_multiplication_counts: MultiplicationCounts,
    build_inversion_profile: BatchInversionProfile,
    reference_signed_sums: u64,
    reference_build_group_additions: u64,
    intermediate_reference_group_additions: u64,
    checkpoints: Vec<OuterCheckpoint>,
    growth_fit: OuterGrowthFit,
    residual_degree_evidence: ResidualDegreeEvidence,
    local_maximum_degree: u32,
    memory: OuterMemoryModel,
    collection_projection: OuterCollectionProjection,
    sparse_linear_algebra_projection: SparseLinearAlgebraProjection,
    relations: Vec<OuterRelation>,
    exact: bool,
    attack_promotion_gate: bool,
}

#[derive(Clone, Default)]
struct OuterTargetAccumulator {
    signed_branch_cells: u64,
    oracle_calls: u64,
    quadratic_solves: u64,
    quadratic_roots_returned: u64,
    algebra_index_lookups: u64,
    algebra_hits: u64,
    candidate_positive_trials: u64,
    reference_positive_trials: u64,
    false_negatives: u64,
    false_positives: u64,
    multiplication_counts: MultiplicationCounts,
    target_adjustment_group_additions: u64,
    reference_index_lookups: u64,
    direct_relation_verification_additions: u64,
    batch_profile: BatchSizeTotals,
    relations: Vec<OuterRelation>,
    trial_stream: String,
}

struct OuterTargetSpec {
    kind: String,
    scalar: Option<BigUint>,
    point: Point,
    planted: bool,
}

#[derive(Clone, Copy)]
struct QuadraticRoots {
    roots: [Fe; 2],
    len: usize,
    linear: bool,
    universal: bool,
}

impl QuadraticRoots {
    fn empty() -> Self {
        Self {
            roots: [Fe::ZERO; 2],
            len: 0,
            linear: false,
            universal: false,
        }
    }
}

#[derive(Default)]
struct CountedField {
    multiplications: u64,
}

impl CountedField {
    fn mul(&mut self, left: Fe, right: Fe) -> Fe {
        self.multiplications += 1;
        left.mul(&right)
    }

    fn square(&mut self, value: Fe) -> Fe {
        self.multiplications += 1;
        value.sqr()
    }

    /// Right-to-left binary exponentiation, with every executed field
    /// multiplication charged.  The final unused base square is deliberately
    /// charged, matching the native loop that was run.
    fn pow(&mut self, mut base: Fe, exponent: &BigUint) -> Fe {
        let mut result = Fe::ONE;
        for bit in 0..exponent.bits() {
            if exponent.bit(bit) {
                result = self.mul(result, base);
            }
            base = self.square(base);
        }
        result
    }
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("not a non-negative integer: `{value}`"))
}

fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn point(x: &BigUint, y: &BigUint, curve: &CurveParams) -> Point {
    Point::Affine {
        x: curve.fe(x.clone()),
        y: curve.fe(y.clone()),
    }
}

fn affine_coordinates(value: &Point) -> Option<(&BigUint, &BigUint)> {
    match value {
        Point::Infinity => None,
        Point::Affine { x, y } => Some((&x.value, &y.value)),
    }
}

fn affine_key(value: &Point) -> Option<([u8; 32], [u8; 32])> {
    let (x, y) = affine_coordinates(value)?;
    Some((big_key(x), big_key(y)))
}

fn big_key(value: &BigUint) -> [u8; 32] {
    let bytes = value.to_bytes_be();
    assert!(bytes.len() <= 32, "P-256 coordinate exceeded 32 bytes");
    let mut raw = [0; 32];
    raw[32 - bytes.len()..].copy_from_slice(&bytes);
    raw
}

fn fe_key(value: Fe) -> [u8; 32] {
    value.to_bytes_be()
}

fn negate(value: &Point, curve: &CurveParams) -> Point {
    match value {
        Point::Infinity => Point::Infinity,
        Point::Affine { x, y } => Point::Affine {
            x: x.clone(),
            y: curve.fe(if y.value.is_zero() {
                BigUint::zero()
            } else {
                &curve.p - &y.value
            }),
        },
    }
}

fn coefficients(u: Fe, t: Fe, a: Fe, b: Fe, count: &mut CountedField) -> (Fe, Fe, Fe) {
    let two_b = b.add(&b);
    let t2 = count.square(t);
    let u2 = count.square(u);
    let a2 = count.square(a);
    let au = count.mul(a, u);
    let qa = count.square(t.sub(&u));

    let mut qb_inner = count.mul(t2, u);
    qb_inner = qb_inner.add(&count.mul(t, u2));
    qb_inner = qb_inner.add(&count.mul(t, a));
    qb_inner = qb_inner.add(&au);
    qb_inner = qb_inner.add(&two_b);
    let qb = qb_inner.add(&qb_inner).neg();

    let mut qc = count.mul(t2, u2);
    let t_inner = count.mul(t, au.add(&two_b));
    qc = qc.sub(&t_inner.add(&t_inner));
    qc = qc.add(&a2);
    let bu = count.mul(b, u);
    qc = qc.sub(&bu.add(&bu).add(&bu.add(&bu)));
    (qa, qb, qc)
}

fn solve_quadratic(
    qa: Fe,
    qb: Fe,
    qc: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    counts: &mut MultiplicationCounts,
) -> QuadraticRoots {
    if qa == Fe::ZERO {
        if qb == Fe::ZERO {
            return QuadraticRoots {
                universal: qc == Fe::ZERO,
                ..QuadraticRoots::empty()
            };
        }
        let mut inverse_field = CountedField::default();
        let inverse = inverse_field.pow(qb, inverse_exponent);
        counts.inversions += inverse_field.multiplications;
        let mut root_field = CountedField::default();
        let root = root_field.mul(qc.neg(), inverse);
        counts.root_construction += root_field.multiplications;
        return QuadraticRoots {
            roots: [root, Fe::ZERO],
            len: 1,
            linear: true,
            universal: false,
        };
    }

    let mut discriminant_field = CountedField::default();
    let qb2 = discriminant_field.square(qb);
    let ac = discriminant_field.mul(qa, qc);
    let four_ac = ac.add(&ac).add(&ac.add(&ac));
    let discriminant = qb2.sub(&four_ac);
    counts.discriminants += discriminant_field.multiplications;

    let mut sqrt_field = CountedField::default();
    let sqrt = sqrt_field.pow(discriminant, sqrt_exponent);
    let sqrt_check = sqrt_field.square(sqrt);
    counts.square_roots += sqrt_field.multiplications;
    if sqrt_check != discriminant {
        return QuadraticRoots::empty();
    }

    let mut inverse_field = CountedField::default();
    let inverse = inverse_field.pow(qa.add(&qa), inverse_exponent);
    counts.inversions += inverse_field.multiplications;

    let mut root_field = CountedField::default();
    let neg_qb = qb.neg();
    let first = root_field.mul(neg_qb.add(&sqrt), inverse);
    let second = root_field.mul(neg_qb.sub(&sqrt), inverse);
    counts.root_construction += root_field.multiplications;
    if first == second {
        QuadraticRoots {
            roots: [first, Fe::ZERO],
            len: 1,
            linear: false,
            universal: false,
        }
    } else {
        QuadraticRoots {
            roots: [first, second],
            len: 2,
            linear: false,
            universal: false,
        }
    }
}

fn roots_from_inverse(
    qb: Fe,
    sqrt: Fe,
    inverse: Fe,
    counts: &mut MultiplicationCounts,
) -> QuadraticRoots {
    let mut root_field = CountedField::default();
    let neg_qb = qb.neg();
    let first = root_field.mul(neg_qb.add(&sqrt), inverse);
    let second = root_field.mul(neg_qb.sub(&sqrt), inverse);
    counts.root_construction += root_field.multiplications;
    if first == second {
        QuadraticRoots {
            roots: [first, Fe::ZERO],
            len: 1,
            linear: false,
            universal: false,
        }
    } else {
        QuadraticRoots {
            roots: [first, second],
            len: 2,
            linear: false,
            universal: false,
        }
    }
}

fn batch_invert(
    denominators: &[Fe],
    inverse_exponent: &BigUint,
    counts: &mut MultiplicationCounts,
) -> Result<Vec<Fe>, String> {
    if denominators.len() < 2 {
        return Err("batch inversion requires at least two denominators".into());
    }
    if denominators.contains(&Fe::ZERO) {
        return Err("batch inversion received a zero denominator".into());
    }
    let mut field = CountedField::default();
    let mut prefixes = Vec::with_capacity(denominators.len());
    let mut product = Fe::ONE;
    for &denominator in denominators {
        prefixes.push(product);
        product = field.mul(product, denominator);
    }
    let mut inverse_product = field.pow(product, inverse_exponent);
    let mut inverses = vec![Fe::ZERO; denominators.len()];
    for index in (0..denominators.len()).rev() {
        inverses[index] = field.mul(inverse_product, prefixes[index]);
        if index != 0 {
            inverse_product = field.mul(inverse_product, denominators[index]);
        }
    }
    counts.inversions += field.multiplications;
    Ok(inverses)
}

fn solve_quadratic_batch(
    coefficients: &[(Fe, Fe, Fe)],
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    counts: &mut MultiplicationCounts,
    profile: &mut BatchInversionProfile,
) -> Result<Vec<QuadraticRoots>, String> {
    if coefficients.len() <= 1 {
        profile.scalar_fallbacks += coefficients.len() as u64;
        return Ok(coefficients
            .iter()
            .map(|&(qa, qb, qc)| {
                solve_quadratic(qa, qb, qc, sqrt_exponent, inverse_exponent, counts)
            })
            .collect());
    }

    let mut roots = vec![QuadraticRoots::empty(); coefficients.len()];
    let mut valid = Vec::new();
    for (index, &(qa, qb, qc)) in coefficients.iter().enumerate() {
        if qa == Fe::ZERO {
            profile.scalar_fallbacks += 1;
            roots[index] = solve_quadratic(qa, qb, qc, sqrt_exponent, inverse_exponent, counts);
            continue;
        }

        let mut discriminant_field = CountedField::default();
        let qb2 = discriminant_field.square(qb);
        let ac = discriminant_field.mul(qa, qc);
        let four_ac = ac.add(&ac).add(&ac.add(&ac));
        let discriminant = qb2.sub(&four_ac);
        counts.discriminants += discriminant_field.multiplications;

        let mut sqrt_field = CountedField::default();
        let sqrt = sqrt_field.pow(discriminant, sqrt_exponent);
        let sqrt_check = sqrt_field.square(sqrt);
        counts.square_roots += sqrt_field.multiplications;
        if sqrt_check == discriminant {
            valid.push((index, qb, sqrt, qa.add(&qa)));
        }
    }

    if valid.is_empty() {
        return Ok(roots);
    }
    if valid.len() == 1 {
        profile.scalar_fallbacks += 1;
        let (index, qb, sqrt, denominator) = valid[0];
        let mut inverse_field = CountedField::default();
        let inverse = inverse_field.pow(denominator, inverse_exponent);
        counts.inversions += inverse_field.multiplications;
        roots[index] = roots_from_inverse(qb, sqrt, inverse, counts);
        return Ok(roots);
    }

    let denominators: Vec<Fe> = valid.iter().map(|entry| entry.3).collect();
    let inverses = batch_invert(&denominators, inverse_exponent, counts)?;
    profile.batch_sizes.push(valid.len());
    profile.batch_inversions += 1;
    for ((index, qb, sqrt, _), inverse) in valid.into_iter().zip(inverses) {
        roots[index] = roots_from_inverse(qb, sqrt, inverse, counts);
    }
    Ok(roots)
}

fn algebra_pairs(
    target: &Point,
    columns: &[Column],
    x_index: &HashMap<[u8; 32], u64>,
    curve: &CurveParams,
) -> Result<(BTreeSet<(u64, u64)>, u64, u64, u64, MultiplicationCounts), String> {
    let (target_x, _) = affine_coordinates(target).ok_or("target is infinity")?;
    let t = Fe::from_biguint(target_x);
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut pairs = BTreeSet::new();
    let mut roots_returned = 0u64;
    let mut linear = 0u64;
    let mut universal = 0u64;
    let mut counts = MultiplicationCounts::default();
    for (left_col, column) in columns.iter().enumerate() {
        let mut coefficient_field = CountedField::default();
        let (qa, qb, qc) = coefficients(column.x, t, a, b, &mut coefficient_field);
        counts.coefficients += coefficient_field.multiplications;
        let roots = solve_quadratic(qa, qb, qc, &sqrt_exponent, &inverse_exponent, &mut counts);
        linear += u64::from(roots.linear);
        universal += u64::from(roots.universal);
        if roots.universal {
            return Err(format!(
                "universal quadratic at left column {left_col} has no bounded image"
            ));
        }
        roots_returned += roots.len as u64;
        for root in roots.roots[..roots.len].iter().copied() {
            if let Some(&right_col) = x_index.get(&fe_key(root)) {
                pairs.insert((left_col as u64, right_col));
            }
        }
    }
    counts.finish();
    Ok((pairs, roots_returned, linear, universal, counts))
}

fn reference_pairs(
    target: &Point,
    signed_rows: &[SignedRow],
    signed_index: &HashMap<([u8; 32], [u8; 32]), u64>,
    curve: &CurveParams,
) -> BTreeSet<(u64, u64)> {
    let a = curve.a_fe();
    let mut pairs = BTreeSet::new();
    for row in signed_rows {
        let right = target.add_vartime(&negate(&row.point, curve), &a);
        if let Some(key) = affine_key(&right) {
            if let Some(&right_col) = signed_index.get(&key) {
                pairs.insert((row.col, right_col));
            }
        }
    }
    pairs
}

fn verify_pairs(
    pairs: &BTreeSet<(u64, u64)>,
    target: &Point,
    columns: &[Column],
    curve: &CurveParams,
) -> Result<u64, String> {
    let a = curve.a_fe();
    let neg_target = negate(target, curve);
    let mut additions = 0u64;
    for &(left, right) in pairs {
        let left_point = &columns[left as usize].low;
        let right_point = &columns[right as usize].low;
        let mut verified = false;
        for left_sign in [false, true] {
            for right_sign in [false, true] {
                let lp = if left_sign {
                    negate(left_point, curve)
                } else {
                    left_point.clone()
                };
                let rp = if right_sign {
                    negate(right_point, curve)
                } else {
                    right_point.clone()
                };
                let sum = lp.add_vartime(&rp, &a);
                additions += 1;
                verified |= sum == *target || sum == neg_target;
            }
        }
        if !verified {
            return Err(format!(
                "emitted column pair ({left}, {right}) has no signed sum +/-Q"
            ));
        }
    }
    Ok(additions)
}

fn pair_digest(pairs: &BTreeSet<(u64, u64)>) -> String {
    let mut text = String::new();
    for (left, right) in pairs {
        writeln!(text, "{left},{right}").expect("writing to String cannot fail");
    }
    hex::encode(sha256(text.as_bytes()))
}

fn run_target(
    kind: &str,
    scalar: Option<&BigUint>,
    target: &Point,
    columns: &[Column],
    signed_rows: &[SignedRow],
    x_index: &HashMap<[u8; 32], u64>,
    signed_index: &HashMap<([u8; 32], [u8; 32]), u64>,
    curve: &CurveParams,
) -> Result<TargetResult, String> {
    if matches!(target, Point::Infinity) || !curve.is_on_curve(target) {
        return Err(format!("{kind} target is infinity or off curve"));
    }
    let (algebra, roots_returned, linear, universal, multiplication_counts) =
        algebra_pairs(target, columns, x_index, curve)?;
    let reference = reference_pairs(target, signed_rows, signed_index, curve);
    let false_negatives = reference.difference(&algebra).count();
    let false_positives = algebra.difference(&reference).count();
    let verification_group_additions = verify_pairs(&algebra, target, columns, curve)?;
    let (x, y) = affine_coordinates(target).expect("target was checked affine");
    Ok(TargetResult {
        kind: kind.into(),
        scalar: scalar.map(lower_hex),
        target: [lower_hex(x), lower_hex(y)],
        rejected_ordered_column_pairs: COLUMNS * COLUMNS,
        quadratic_solves: COLUMNS,
        quadratic_roots_returned: roots_returned,
        algebra_index_lookups: roots_returned,
        linear_degeneracies: linear,
        universal_degeneracies: universal,
        algebra_pairs: algebra.len(),
        reference_pairs: reference.len(),
        false_negatives,
        false_positives,
        multiplication_counts,
        reference_group_subtractions: SIGNED_POINTS,
        reference_index_lookups: SIGNED_POINTS,
        verification_group_additions,
        all_emitted_pairs_verified: true,
        exact: false_negatives == 0 && false_positives == 0,
        pair_sha256: pair_digest(&algebra),
        pairs: algebra.into_iter().map(|(l, r)| [l, r]).collect(),
    })
}

fn build_indexes(
    curve: &CurveParams,
) -> Result<
    (
        FactorBaseReceipt,
        Vec<Column>,
        Vec<SignedRow>,
        HashMap<[u8; 32], u64>,
        HashMap<([u8; 32], [u8; 32]), u64>,
    ),
    String,
> {
    let built = p256_dickson_factor_base::build(SPEC)?;
    let dump = &built.dump;
    if dump.factor_base.fb_id != FB_ID
        || dump.factor_base.fb_sha256 != FB_SHA256
        || dump.factor_base.points_sha256 != POINTS_SHA256
        || dump.factor_base.columns != COLUMNS
        || dump.factor_base.signed_points != SIGNED_POINTS
        || dump.factor_base.params.get("terminal").map(String::as_str) != Some(TERMINAL)
    {
        return Err("rebuilt factor base does not match the frozen FB1 identity".into());
    }
    p256_dickson_factor_base::verify(dump)?;

    let mut columns: Vec<Option<Column>> = vec![None; COLUMNS as usize];
    let mut signed_rows = Vec::with_capacity(SIGNED_POINTS as usize);
    let mut x_index = HashMap::with_capacity(COLUMNS as usize);
    let mut signed_index = HashMap::with_capacity(SIGNED_POINTS as usize);
    for row in &dump.points {
        let x = parse_big(&row.x)?;
        let y = parse_big(&row.y)?;
        let affine = point(&x, &y, curve);
        if !curve.is_on_curve(&affine) {
            return Err(format!(
                "factor-base row in column {} is off curve",
                row.col
            ));
        }
        let key = affine_key(&affine).expect("factor-base rows are affine");
        if signed_index.insert(key, row.col).is_some() {
            return Err("duplicate signed factor-base point".into());
        }
        signed_rows.push(SignedRow {
            point: affine.clone(),
            col: row.col,
        });
        if row.coef == "1" {
            let x_fe = Fe::from_biguint(&x);
            if x_index.insert(fe_key(x_fe), row.col).is_some() {
                return Err("duplicate factor-base abscissa".into());
            }
            if columns[row.col as usize]
                .replace(Column {
                    x: x_fe,
                    low: affine,
                })
                .is_some()
            {
                return Err(format!("duplicate low-y row for column {}", row.col));
            }
        }
    }
    let columns: Vec<Column> = columns
        .into_iter()
        .enumerate()
        .map(|(index, column)| column.ok_or_else(|| format!("missing column {index}")))
        .collect::<Result<_, _>>()?;
    if columns.len() as u64 != COLUMNS
        || signed_rows.len() as u64 != SIGNED_POINTS
        || x_index.len() as u64 != COLUMNS
        || signed_index.len() as u64 != SIGNED_POINTS
    {
        return Err("factor-base index cardinality mismatch".into());
    }
    Ok((
        FactorBaseReceipt {
            spec: SPEC.into(),
            fb_id: FB_ID.into(),
            fb_sha256: FB_SHA256.into(),
            points_sha256: POINTS_SHA256.into(),
            terminal: TERMINAL.into(),
            columns: COLUMNS,
            signed_points: SIGNED_POINTS,
            complete_rebuild_verified: true,
        },
        columns,
        signed_rows,
        x_index,
        signed_index,
    ))
}

fn full_point_key(value: &Point) -> [u8; 65] {
    let mut key = [0u8; 65];
    if let Some((x, y)) = affine_coordinates(value) {
        key[0] = 1;
        key[1..33].copy_from_slice(&big_key(x));
        key[33..].copy_from_slice(&big_key(y));
    }
    key
}

fn reference_image(
    points: &[Point],
    curve: &CurveParams,
) -> Result<(ReferenceImageState, u64), String> {
    let a = curve.a_fe();
    let mut sums = BTreeMap::from([(full_point_key(&Point::Infinity), Point::Infinity)]);
    let mut additions = 0u64;
    for point in points {
        let negated = negate(point, curve);
        if negated == *point {
            return Err("factor-base leaf has only one signed lift".into());
        }
        let mut next = BTreeMap::new();
        for sum in sums.values() {
            for signed in [point, &negated] {
                let value = sum.add_vartime(signed, &a);
                additions += 1;
                next.insert(full_point_key(&value), value);
            }
        }
        sums = next;
    }
    let mut affine = BTreeSet::new();
    let mut identity = false;
    for point in sums.values() {
        match affine_coordinates(point) {
            Some((x, _)) => {
                affine.insert(big_key(x));
            }
            None => identity = true,
        }
    }
    Ok((ReferenceImageState { affine, identity }, additions))
}

fn image_matches_reference(candidate: &P256ImageState, reference: &ReferenceImageState) -> bool {
    candidate.identity == reference.identity
        && candidate.affine.len() == reference.affine.len()
        && candidate
            .affine
            .keys()
            .zip(reference.affine.iter())
            .all(|(left, right)| left == right)
}

#[allow(clippy::too_many_arguments)]
fn compose_p256_images(
    left: &P256ImageState,
    right: &P256ImageState,
    a: Fe,
    b: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    solves: &mut u64,
    roots_returned: &mut u64,
    linear: &mut u64,
    universal: &mut u64,
    counts: &mut MultiplicationCounts,
) -> Result<P256ImageState, String> {
    let mut affine = BTreeMap::new();
    if left.identity {
        affine.extend(right.affine.iter().map(|(key, value)| (*key, *value)));
    }
    if right.identity {
        affine.extend(left.affine.iter().map(|(key, value)| (*key, *value)));
    }
    for &u in left.affine.values() {
        for &v in right.affine.values() {
            *solves += 1;
            let mut coefficient_field = CountedField::default();
            let (qa, qb, qc) = coefficients(u, v, a, b, &mut coefficient_field);
            counts.coefficients += coefficient_field.multiplications;
            let roots = solve_quadratic(qa, qb, qc, sqrt_exponent, inverse_exponent, counts);
            *linear += u64::from(roots.linear);
            *universal += u64::from(roots.universal);
            if roots.universal {
                return Err("universal quadratic has no bounded image".into());
            }
            *roots_returned += roots.len as u64;
            for root in roots.roots[..roots.len].iter().copied() {
                affine.insert(fe_key(root), root);
            }
        }
    }
    let shared_affine = left.affine.keys().any(|key| right.affine.contains_key(key));
    Ok(P256ImageState {
        affine,
        identity: (left.identity && right.identity) || shared_affine,
    })
}

#[allow(clippy::too_many_arguments)]
fn compose_p256_images_batched(
    left: &P256ImageState,
    right: &P256ImageState,
    a: Fe,
    b: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    solves: &mut u64,
    roots_returned: &mut u64,
    linear: &mut u64,
    universal: &mut u64,
    counts: &mut MultiplicationCounts,
    profile: &mut BatchInversionProfile,
) -> Result<P256ImageState, String> {
    let mut affine = BTreeMap::new();
    if left.identity {
        affine.extend(right.affine.iter().map(|(key, value)| (*key, *value)));
    }
    if right.identity {
        affine.extend(left.affine.iter().map(|(key, value)| (*key, *value)));
    }
    let mut coefficient_rows = Vec::with_capacity(left.affine.len() * right.affine.len());
    for &u in left.affine.values() {
        for &v in right.affine.values() {
            let mut coefficient_field = CountedField::default();
            coefficient_rows.push(coefficients(u, v, a, b, &mut coefficient_field));
            counts.coefficients += coefficient_field.multiplications;
        }
    }
    *solves += coefficient_rows.len() as u64;
    let roots = solve_quadratic_batch(
        &coefficient_rows,
        sqrt_exponent,
        inverse_exponent,
        counts,
        profile,
    )?;
    for result in roots {
        *linear += u64::from(result.linear);
        *universal += u64::from(result.universal);
        if result.universal {
            return Err("universal quadratic has no bounded image".into());
        }
        *roots_returned += result.len as u64;
        for root in result.roots[..result.len].iter().copied() {
            affine.insert(fe_key(root), root);
        }
    }
    let shared_affine = left.affine.keys().any(|key| right.affine.contains_key(key));
    Ok(P256ImageState {
        affine,
        identity: (left.identity && right.identity) || shared_affine,
    })
}

fn image_digest(image: &P256ImageState) -> String {
    let mut bytes = Vec::with_capacity(1 + 32 * image.affine.len());
    bytes.push(u8::from(image.identity));
    for key in image.affine.keys() {
        bytes.extend_from_slice(key);
    }
    hex::encode(sha256(&bytes))
}

fn verify_width_node(
    kind: &str,
    level: usize,
    node: usize,
    candidate: &P256ImageState,
    points: &[Point],
    curve: &CurveParams,
    reference_group_additions: &mut u64,
) -> Result<(), String> {
    let (reference, additions) = reference_image(points, curve)?;
    *reference_group_additions += additions;
    if !image_matches_reference(candidate, &reference) {
        return Err(format!(
            "{kind} {level}-leaf node {node} differs from exact signed group image"
        ));
    }
    Ok(())
}

fn run_width_sample(
    kind: &str,
    column_indices: Vec<u64>,
    columns: &[Column],
    curve: &CurveParams,
) -> Result<WidthSampleResult, String> {
    if column_indices.len() != 16
        || column_indices
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .len()
            != 16
    {
        return Err(format!("{kind} does not contain 16 distinct columns"));
    }
    let selected: Vec<&Column> = column_indices
        .iter()
        .map(|&index| {
            columns
                .get(index as usize)
                .ok_or_else(|| format!("{kind} column {index} is out of range"))
        })
        .collect::<Result<_, _>>()?;
    let points: Vec<Point> = selected.iter().map(|column| column.low.clone()).collect();
    let x_coordinates: Vec<String> = selected
        .iter()
        .map(|column| lower_hex(&BigUint::from_bytes_be(&column.x.to_bytes_be())))
        .collect();
    let mut images: Vec<P256ImageState> = selected
        .iter()
        .map(|column| P256ImageState {
            affine: BTreeMap::from([(fe_key(column.x), column.x)]),
            identity: false,
        })
        .collect();
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut solves = 0u64;
    let mut roots_returned = 0u64;
    let mut linear = 0u64;
    let mut universal = 0u64;
    let mut counts = MultiplicationCounts::default();
    let mut reference_group_additions = 0u64;
    let mut widths_by_level: Vec<Vec<usize>> = Vec::new();
    let mut identity_images_by_level = [0u64; 4];

    for (level_index, block_size) in [2usize, 4, 8, 16].into_iter().enumerate() {
        let mut next = Vec::with_capacity(images.len() / 2);
        for (node, pair) in images.chunks_exact(2).enumerate() {
            let image = compose_p256_images(
                &pair[0],
                &pair[1],
                a,
                b,
                &sqrt_exponent,
                &inverse_exponent,
                &mut solves,
                &mut roots_returned,
                &mut linear,
                &mut universal,
                &mut counts,
            )?;
            let start = node * block_size;
            verify_width_node(
                kind,
                block_size,
                node,
                &image,
                &points[start..start + block_size],
                curve,
                &mut reference_group_additions,
            )?;
            identity_images_by_level[level_index] += u64::from(image.identity);
            next.push(image);
        }
        widths_by_level.push(next.iter().map(|image| image.affine.len()).collect());
        images = next;
    }
    if images.len() != 1 {
        return Err(format!("{kind} did not reduce to one sixteen-leaf image"));
    }
    counts.finish();
    let final_image = &images[0];
    Ok(WidthSampleResult {
        kind: kind.into(),
        columns: column_indices,
        x_coordinates,
        two_leaf_widths: widths_by_level[0].clone(),
        four_leaf_widths: widths_by_level[1].clone(),
        eight_leaf_widths: widths_by_level[2].clone(),
        sixteen_leaf_width: widths_by_level[3][0],
        identity_images_by_level,
        quadratic_solves: solves,
        quadratic_roots_returned: roots_returned,
        linear_degeneracies: linear,
        universal_degeneracies: universal,
        multiplication_counts: counts,
        reference_group_additions,
        exact: true,
        final_image_sha256: image_digest(final_image),
    })
}

fn hash_sample_indices(sample: usize) -> Vec<u64> {
    let mut indices = Vec::with_capacity(16);
    for leaf in 0..16 {
        let mut counter = 0u64;
        loop {
            let preimage = format!(
                "{CURVE_SLUG}/s17-image-width-round14/sample/{sample}/leaf/{leaf}/counter/{counter}"
            );
            let digest = sha256(preimage.as_bytes());
            let index = (BigUint::from_bytes_be(&digest) % BigUint::from(COLUMNS))
                .to_u64()
                .expect("reduced column index fits u64");
            if !indices.contains(&index) {
                indices.push(index);
                break;
            }
            counter += 1;
        }
    }
    indices
}

fn round14_sample_specs() -> Vec<(String, Vec<u64>)> {
    vec![
        ("prefix".to_string(), (0..16).collect::<Vec<u64>>()),
        ("hash-0".to_string(), hash_sample_indices(0)),
        ("hash-1".to_string(), hash_sample_indices(1)),
        ("hash-2".to_string(), hash_sample_indices(2)),
    ]
}

fn run_width(round13_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round13_path).map_err(|error| error.to_string())?;
    let round13_sha256 = hex::encode(sha256(&bytes));
    if round13_sha256 != ROUND13_SHA256 {
        return Err(format!(
            "round-13 SHA-256 mismatch: expected {ROUND13_SHA256}, got {round13_sha256}"
        ));
    }
    let curve = CurveParams::p256();
    let (factor_base, columns, _, _, _) = build_indexes(&curve)?;
    let sample_specs = round14_sample_specs();
    let mut samples = Vec::new();
    for (kind, indices) in sample_specs {
        samples.push(run_width_sample(&kind, indices, &columns, &curve)?);
    }
    let unordered_factor_base_pairs = COLUMNS * (COLUMNS + 1) / 2;
    let two_root_pair_image_entry_upper_bound = 2 * unordered_factor_base_pairs;
    let identity_aware_pair_image_affine_entries = COLUMNS * COLUMNS;
    let result = WidthExperimentResult {
        schema: "p256.s17_image_width_transfer/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        round13_sha256,
        factor_base,
        local_maximum_degree: 2,
        generic_widths: [2, 8, 128, 32_768],
        samples,
        storage_model: WidthStorageModel {
            generic_sixteen_leaf_x_entries: 32_768,
            bytes_per_x: 32,
            generic_sixteen_leaf_raw_bytes: 32_768 * 32,
            factor_base_columns: COLUMNS,
            unordered_factor_base_pairs,
            two_root_pair_image_entry_upper_bound,
            two_root_pair_image_raw_byte_upper_bound: two_root_pair_image_entry_upper_bound * 32,
            identity_aware_pair_image_affine_entries,
            identity_aware_pair_image_raw_bytes: identity_aware_pair_image_affine_entries * 32,
        },
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn build_verified_eight_image(
    label: &str,
    xs: &[Fe],
    points: &[Point],
    curve: &CurveParams,
    a: Fe,
    b: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    solves: &mut u64,
    roots_returned: &mut u64,
    linear: &mut u64,
    universal: &mut u64,
    counts: &mut MultiplicationCounts,
    reference_group_additions: &mut u64,
) -> Result<P256ImageState, String> {
    if xs.len() != 8 || points.len() != 8 {
        return Err(format!("{label} requires exactly eight leaves"));
    }
    let mut images: Vec<P256ImageState> = xs
        .iter()
        .copied()
        .map(|x| P256ImageState {
            affine: BTreeMap::from([(fe_key(x), x)]),
            identity: false,
        })
        .collect();
    for block_size in [2usize, 4, 8] {
        let mut next = Vec::with_capacity(images.len() / 2);
        for (node, pair) in images.chunks_exact(2).enumerate() {
            let image = compose_p256_images(
                &pair[0],
                &pair[1],
                a,
                b,
                sqrt_exponent,
                inverse_exponent,
                solves,
                roots_returned,
                linear,
                universal,
                counts,
            )?;
            let start = node * block_size;
            verify_width_node(
                label,
                block_size,
                node,
                &image,
                &points[start..start + block_size],
                curve,
                reference_group_additions,
            )?;
            next.push(image);
        }
        images = next;
    }
    if images.len() != 1 {
        return Err(format!("{label} did not reduce to one eight-leaf image"));
    }
    Ok(images.pop().expect("one image was checked"))
}

#[allow(clippy::too_many_arguments)]
fn build_verified_eight_image_batched(
    label: &str,
    xs: &[Fe],
    points: &[Point],
    curve: &CurveParams,
    a: Fe,
    b: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    solves: &mut u64,
    roots_returned: &mut u64,
    linear: &mut u64,
    universal: &mut u64,
    counts: &mut MultiplicationCounts,
    profile: &mut BatchInversionProfile,
    reference_group_additions: &mut u64,
) -> Result<P256ImageState, String> {
    if xs.len() != 8 || points.len() != 8 {
        return Err(format!("{label} requires exactly eight leaves"));
    }
    let mut images: Vec<P256ImageState> = xs
        .iter()
        .copied()
        .map(|x| P256ImageState {
            affine: BTreeMap::from([(fe_key(x), x)]),
            identity: false,
        })
        .collect();
    for block_size in [2usize, 4, 8] {
        let mut next = Vec::with_capacity(images.len() / 2);
        for (node, pair) in images.chunks_exact(2).enumerate() {
            let image = compose_p256_images_batched(
                &pair[0],
                &pair[1],
                a,
                b,
                sqrt_exponent,
                inverse_exponent,
                solves,
                roots_returned,
                linear,
                universal,
                counts,
                profile,
            )?;
            let start = node * block_size;
            verify_width_node(
                label,
                block_size,
                node,
                &image,
                &points[start..start + block_size],
                curve,
                reference_group_additions,
            )?;
            next.push(image);
        }
        images = next;
    }
    if images.len() != 1 {
        return Err(format!("{label} did not reduce to one eight-leaf image"));
    }
    Ok(images.pop().expect("one image was checked"))
}

#[allow(clippy::too_many_arguments)]
fn query_compressed_target(
    kind: &str,
    scalar: Option<&BigUint>,
    target: &Point,
    reference_positive: bool,
    left: &P256ImageState,
    right: &P256ImageState,
    curve: &CurveParams,
    build_solves: u64,
    build_counts: MultiplicationCounts,
) -> Result<CompressedTargetResult, String> {
    if matches!(target, Point::Infinity) || !curve.is_on_curve(target) {
        return Err(format!("{kind} target is infinity or off curve"));
    }
    let (target_x, target_y) = affine_coordinates(target).expect("target was checked affine");
    let t = Fe::from_biguint(target_x);
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut solves = 0u64;
    let mut roots_returned = 0u64;
    let mut linear = 0u64;
    let mut universal = 0u64;
    let mut counts = MultiplicationCounts::default();
    let mut hits = BTreeSet::new();
    for (&left_key, &left_x) in &left.affine {
        solves += 1;
        let mut coefficient_field = CountedField::default();
        let (qa, qb, qc) = coefficients(left_x, t, a, b, &mut coefficient_field);
        counts.coefficients += coefficient_field.multiplications;
        let roots = solve_quadratic(qa, qb, qc, &sqrt_exponent, &inverse_exponent, &mut counts);
        linear += u64::from(roots.linear);
        universal += u64::from(roots.universal);
        if roots.universal {
            return Err(format!("{kind} query produced a universal quadratic"));
        }
        roots_returned += roots.len as u64;
        for root in roots.roots[..roots.len].iter().copied() {
            let right_key = fe_key(root);
            if right.affine.contains_key(&right_key) {
                let mut record = Vec::with_capacity(65);
                record.push(0);
                record.extend_from_slice(&left_key);
                record.extend_from_slice(&right_key);
                hits.insert(record);
            }
        }
    }
    let target_key = big_key(target_x);
    if left.identity && right.affine.contains_key(&target_key) {
        let mut record = Vec::with_capacity(33);
        record.push(1);
        record.extend_from_slice(&target_key);
        hits.insert(record);
    }
    if right.identity && left.affine.contains_key(&target_key) {
        let mut record = Vec::with_capacity(33);
        record.push(2);
        record.extend_from_slice(&target_key);
        hits.insert(record);
    }
    counts.finish();
    let candidate_positive = !hits.is_empty();
    let false_negative = reference_positive && !candidate_positive;
    let false_positive = !reference_positive && candidate_positive;
    let mut hit_bytes = Vec::new();
    for hit in &hits {
        hit_bytes.extend_from_slice(hit);
    }
    let cold_counts = build_counts.combined(counts);
    Ok(CompressedTargetResult {
        kind: kind.into(),
        scalar: scalar.map(lower_hex),
        target: [lower_hex(target_x), lower_hex(target_y)],
        reference_positive,
        candidate_positive,
        algebra_hits: hits.len() as u64,
        quadratic_solves_warm: solves,
        quadratic_solves_cold: build_solves + solves,
        quadratic_roots_returned: roots_returned,
        algebra_index_lookups: roots_returned,
        linear_degeneracies: linear,
        universal_degeneracies: universal,
        multiplication_counts_warm: counts,
        multiplication_counts_cold: cold_counts,
        false_negative,
        false_positive,
        exact: !false_negative && !false_positive,
        hit_sha256: hex::encode(sha256(&hit_bytes)),
    })
}

#[allow(clippy::too_many_arguments)]
fn observe_batched_target(
    target: &Point,
    left: &P256ImageState,
    right: &P256ImageState,
    curve: &CurveParams,
    a: Fe,
    b: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
) -> Result<BatchedQueryObservation, String> {
    if matches!(target, Point::Infinity) || !curve.is_on_curve(target) {
        return Err("batched observation target is infinity or off curve".into());
    }
    let (target_x, _) = affine_coordinates(target).expect("target was checked affine");
    let t = Fe::from_biguint(target_x);
    let mut counts = MultiplicationCounts::default();
    let mut profile = BatchInversionProfile::default();
    let mut left_keys = Vec::with_capacity(left.affine.len());
    let mut coefficient_rows = Vec::with_capacity(left.affine.len());
    for (&left_key, &left_x) in &left.affine {
        let mut coefficient_field = CountedField::default();
        coefficient_rows.push(coefficients(left_x, t, a, b, &mut coefficient_field));
        counts.coefficients += coefficient_field.multiplications;
        left_keys.push(left_key);
    }
    let results = solve_quadratic_batch(
        &coefficient_rows,
        sqrt_exponent,
        inverse_exponent,
        &mut counts,
        &mut profile,
    )?;
    let mut roots_returned = 0u64;
    let mut linear = 0u64;
    let mut universal = 0u64;
    let mut hits = BTreeSet::new();
    for (left_key, result) in left_keys.into_iter().zip(results) {
        linear += u64::from(result.linear);
        universal += u64::from(result.universal);
        if result.universal {
            return Err("batched observation produced a universal quadratic".into());
        }
        roots_returned += result.len as u64;
        for root in result.roots[..result.len].iter().copied() {
            let right_key = fe_key(root);
            if right.affine.contains_key(&right_key) {
                let mut record = Vec::with_capacity(65);
                record.push(0);
                record.extend_from_slice(&left_key);
                record.extend_from_slice(&right_key);
                hits.insert(record);
            }
        }
    }
    let target_key = big_key(target_x);
    if left.identity && right.affine.contains_key(&target_key) {
        let mut record = Vec::with_capacity(33);
        record.push(1);
        record.extend_from_slice(&target_key);
        hits.insert(record);
    }
    if right.identity && left.affine.contains_key(&target_key) {
        let mut record = Vec::with_capacity(33);
        record.push(2);
        record.extend_from_slice(&target_key);
        hits.insert(record);
    }
    counts.finish();
    let mut hit_bytes = Vec::new();
    for hit in &hits {
        hit_bytes.extend_from_slice(hit);
    }
    Ok(BatchedQueryObservation {
        positive: !hits.is_empty(),
        algebra_hits: hits.len() as u64,
        quadratic_solves: coefficient_rows.len() as u64,
        quadratic_roots_returned: roots_returned,
        algebra_index_lookups: roots_returned,
        linear_degeneracies: linear,
        universal_degeneracies: universal,
        multiplication_counts: counts,
        inversion_profile: profile,
        hit_sha256: hex::encode(sha256(&hit_bytes)),
    })
}

#[allow(clippy::too_many_arguments)]
fn query_batched_target(
    kind: &str,
    scalar: Option<&BigUint>,
    target: &Point,
    reference_positive: bool,
    left: &P256ImageState,
    right: &P256ImageState,
    curve: &CurveParams,
    scalar_result: &CompressedTargetResult,
    scalar_build_counts: MultiplicationCounts,
    batched_build_counts: MultiplicationCounts,
) -> Result<BatchedTargetResult, String> {
    if matches!(target, Point::Infinity) || !curve.is_on_curve(target) {
        return Err(format!("{kind} target is infinity or off curve"));
    }
    let (target_x, target_y) = affine_coordinates(target).expect("target was checked affine");
    let t = Fe::from_biguint(target_x);
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut counts = MultiplicationCounts::default();
    let mut profile = BatchInversionProfile::default();
    let mut left_keys = Vec::with_capacity(left.affine.len());
    let mut coefficient_rows = Vec::with_capacity(left.affine.len());
    for (&left_key, &left_x) in &left.affine {
        let mut coefficient_field = CountedField::default();
        coefficient_rows.push(coefficients(left_x, t, a, b, &mut coefficient_field));
        counts.coefficients += coefficient_field.multiplications;
        left_keys.push(left_key);
    }
    let results = solve_quadratic_batch(
        &coefficient_rows,
        &sqrt_exponent,
        &inverse_exponent,
        &mut counts,
        &mut profile,
    )?;
    let mut roots_returned = 0u64;
    let mut linear = 0u64;
    let mut universal = 0u64;
    let mut hits = BTreeSet::new();
    for (left_key, result) in left_keys.into_iter().zip(results) {
        linear += u64::from(result.linear);
        universal += u64::from(result.universal);
        if result.universal {
            return Err(format!(
                "{kind} batched query produced a universal quadratic"
            ));
        }
        roots_returned += result.len as u64;
        for root in result.roots[..result.len].iter().copied() {
            let right_key = fe_key(root);
            if right.affine.contains_key(&right_key) {
                let mut record = Vec::with_capacity(65);
                record.push(0);
                record.extend_from_slice(&left_key);
                record.extend_from_slice(&right_key);
                hits.insert(record);
            }
        }
    }
    let target_key = big_key(target_x);
    if left.identity && right.affine.contains_key(&target_key) {
        let mut record = Vec::with_capacity(33);
        record.push(1);
        record.extend_from_slice(&target_key);
        hits.insert(record);
    }
    if right.identity && left.affine.contains_key(&target_key) {
        let mut record = Vec::with_capacity(33);
        record.push(2);
        record.extend_from_slice(&target_key);
        hits.insert(record);
    }
    counts.finish();
    let batched_positive = !hits.is_empty();
    let false_negative = reference_positive && !batched_positive;
    let false_positive = !reference_positive && batched_positive;
    let mut hit_bytes = Vec::new();
    for hit in &hits {
        hit_bytes.extend_from_slice(hit);
    }
    let hit_sha256 = hex::encode(sha256(&hit_bytes));
    let batched_cold = batched_build_counts.combined(counts);
    let scalar_cold = scalar_build_counts.combined(scalar_result.multiplication_counts_warm);
    let target_coordinates = [lower_hex(target_x), lower_hex(target_y)];
    let scalar_text = scalar.map(lower_hex);
    let exact = !false_negative
        && !false_positive
        && scalar_result.kind == kind
        && scalar_result.scalar == scalar_text
        && scalar_result.target == target_coordinates
        && scalar_result.reference_positive == reference_positive
        && scalar_result.candidate_positive == batched_positive
        && scalar_result.algebra_hits == hits.len() as u64
        && scalar_result.quadratic_solves_warm == coefficient_rows.len() as u64
        && scalar_result.quadratic_roots_returned == roots_returned
        && scalar_result.algebra_index_lookups == roots_returned
        && scalar_result.linear_degeneracies == linear
        && scalar_result.universal_degeneracies == universal
        && scalar_result.hit_sha256 == hit_sha256;
    let warm_multiplication_reduction =
        scalar_result.multiplication_counts_warm.total as f64 / counts.total as f64;
    let cold_multiplication_reduction = scalar_cold.total as f64 / batched_cold.total as f64;
    if !exact || warm_multiplication_reduction < 2.0 || cold_multiplication_reduction < 2.0 {
        return Err(format!(
            "{kind} batched target gate failed: exact={exact}, warm={warm_multiplication_reduction:.6}, cold={cold_multiplication_reduction:.6}"
        ));
    }
    Ok(BatchedTargetResult {
        kind: kind.into(),
        scalar: scalar_text,
        target: target_coordinates,
        reference_positive,
        scalar_positive: scalar_result.candidate_positive,
        batched_positive,
        algebra_hits: hits.len() as u64,
        quadratic_solves: coefficient_rows.len() as u64,
        quadratic_roots_returned: roots_returned,
        algebra_index_lookups: roots_returned,
        linear_degeneracies: linear,
        universal_degeneracies: universal,
        scalar_multiplication_counts_warm: scalar_result.multiplication_counts_warm,
        batched_multiplication_counts_warm: counts,
        scalar_multiplication_counts_cold: scalar_cold,
        batched_multiplication_counts_cold: batched_cold,
        warm_multiplication_reduction,
        cold_multiplication_reduction,
        inversion_profile: profile,
        false_negative,
        false_positive,
        exact,
        hit_sha256,
    })
}

fn verify_round14_samples(
    dependency: &serde_json::Value,
    sample_specs: &[(String, Vec<u64>)],
    columns: &[Column],
) -> Result<(), String> {
    let samples = dependency
        .get("samples")
        .and_then(serde_json::Value::as_array)
        .ok_or("round-14 dependency has no samples array")?;
    if samples.len() != sample_specs.len() {
        return Err("round-14 dependency sample count changed".into());
    }
    for (sample, (expected_kind, expected_columns)) in samples.iter().zip(sample_specs) {
        if sample.get("kind").and_then(serde_json::Value::as_str) != Some(expected_kind) {
            return Err(format!("round-14 sample name changed for {expected_kind}"));
        }
        let stored_columns: Vec<u64> = sample
            .get("columns")
            .and_then(serde_json::Value::as_array)
            .ok_or_else(|| format!("round-14 {expected_kind} columns are missing"))?
            .iter()
            .map(|value| {
                value
                    .as_u64()
                    .ok_or_else(|| format!("round-14 {expected_kind} has a non-u64 column"))
            })
            .collect::<Result<_, _>>()?;
        if &stored_columns != expected_columns {
            return Err(format!("round-14 {expected_kind} column selection changed"));
        }
        let stored_xs: Vec<&str> = sample
            .get("x_coordinates")
            .and_then(serde_json::Value::as_array)
            .ok_or_else(|| format!("round-14 {expected_kind} x-coordinates are missing"))?
            .iter()
            .map(|value| {
                value
                    .as_str()
                    .ok_or_else(|| format!("round-14 {expected_kind} has a non-string x"))
            })
            .collect::<Result<_, _>>()?;
        let rebuilt_xs: Vec<String> = expected_columns
            .iter()
            .map(|&index| {
                lower_hex(&BigUint::from_bytes_be(
                    &columns[index as usize].x.to_bytes_be(),
                ))
            })
            .collect();
        if stored_xs
            .iter()
            .zip(&rebuilt_xs)
            .any(|(stored, rebuilt)| *stored != rebuilt)
        {
            return Err(format!("round-14 {expected_kind} x-coordinates changed"));
        }
    }
    Ok(())
}

fn run_compression_sample(
    kind: &str,
    column_indices: Vec<u64>,
    columns: &[Column],
    curve: &CurveParams,
) -> Result<CompressionSampleResult, String> {
    let selected: Vec<&Column> = column_indices
        .iter()
        .map(|&index| {
            columns
                .get(index as usize)
                .ok_or_else(|| format!("{kind} column {index} is out of range"))
        })
        .collect::<Result<_, _>>()?;
    if selected.len() != 16 {
        return Err(format!("{kind} requires exactly sixteen columns"));
    }
    let xs: Vec<Fe> = selected.iter().map(|column| column.x).collect();
    let points: Vec<Point> = selected.iter().map(|column| column.low.clone()).collect();
    let x_coordinates: Vec<String> = xs
        .iter()
        .map(|x| lower_hex(&BigUint::from_bytes_be(&x.to_bytes_be())))
        .collect();
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut build_solves = 0u64;
    let mut build_roots = 0u64;
    let mut build_linear = 0u64;
    let mut build_universal = 0u64;
    let mut build_counts = MultiplicationCounts::default();
    let mut intermediate_reference_group_additions = 0u64;
    let left = build_verified_eight_image(
        &format!("{kind}/left"),
        &xs[..8],
        &points[..8],
        curve,
        a,
        b,
        &sqrt_exponent,
        &inverse_exponent,
        &mut build_solves,
        &mut build_roots,
        &mut build_linear,
        &mut build_universal,
        &mut build_counts,
        &mut intermediate_reference_group_additions,
    )?;
    let right = build_verified_eight_image(
        &format!("{kind}/right"),
        &xs[8..],
        &points[8..],
        curve,
        a,
        b,
        &sqrt_exponent,
        &inverse_exponent,
        &mut build_solves,
        &mut build_roots,
        &mut build_linear,
        &mut build_universal,
        &mut build_counts,
        &mut intermediate_reference_group_additions,
    )?;
    build_counts.finish();
    if left.affine.len() != 128
        || right.affine.len() != 128
        || left.identity
        || right.identity
        || build_solves != 152
        || build_roots != 304
    {
        return Err(format!("{kind} half-image shape or count changed"));
    }

    let (full_reference, full_reference_group_additions) = reference_image(&points, curve)?;
    let curve_a = curve.a_fe();
    let planted = points.iter().fold(Point::Infinity, |sum, point| {
        sum.add_vartime(point, &curve_a)
    });
    if matches!(planted, Point::Infinity) || !curve.is_on_curve(&planted) {
        return Err(format!("{kind} planted target is infinity or off curve"));
    }
    let public_preimage = format!("{CURVE_SLUG}/s17-target-join-round15/sample/{kind}/public");
    let mut public_scalar = BigUint::from_bytes_be(&sha256(public_preimage.as_bytes())) % &curve.n;
    if public_scalar.is_zero() {
        public_scalar = BigUint::one();
    }
    let public = curve
        .generator()
        .scalar_mul_vartime(&public_scalar, &curve_a);

    let reference_membership = |target: &Point| -> Result<bool, String> {
        let (x, _) = affine_coordinates(target).ok_or("compressed target is infinity")?;
        Ok(full_reference.affine.contains(&big_key(x)))
    };
    let targets = vec![
        query_compressed_target(
            "planted-positive",
            None,
            &planted,
            reference_membership(&planted)?,
            &left,
            &right,
            curve,
            build_solves,
            build_counts,
        )?,
        query_compressed_target(
            "hash-public",
            Some(&public_scalar),
            &public,
            reference_membership(&public)?,
            &left,
            &right,
            curve,
            build_solves,
            build_counts,
        )?,
    ];
    if targets.iter().any(|target| !target.exact) {
        return Err(format!("{kind} compressed target classification mismatch"));
    }
    let retained_entries = left.affine.len() + right.affine.len();
    let materialized_sixteen_leaf_entries = 32_768usize;
    Ok(CompressionSampleResult {
        kind: kind.into(),
        columns: column_indices,
        x_coordinates,
        left_eight_leaf_entries: left.affine.len(),
        right_eight_leaf_entries: right.affine.len(),
        retained_entries,
        materialized_sixteen_leaf_entries,
        retained_raw_bytes: retained_entries * 32,
        materialized_raw_bytes: materialized_sixteen_leaf_entries * 32,
        retained_to_materialized_ratio: retained_entries as f64
            / materialized_sixteen_leaf_entries as f64,
        build_quadratic_solves: build_solves,
        full_materialization_quadratic_solves: 16_536,
        intermediate_reference_group_additions,
        full_reference_group_additions,
        build_linear_degeneracies: build_linear,
        build_universal_degeneracies: build_universal,
        build_multiplication_counts: build_counts,
        left_image_sha256: image_digest(&left),
        right_image_sha256: image_digest(&right),
        targets,
    })
}

fn run_compression(round14_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round14_path).map_err(|error| error.to_string())?;
    let round14_sha256 = hex::encode(sha256(&bytes));
    if round14_sha256 != ROUND14_SHA256 {
        return Err(format!(
            "round-14 SHA-256 mismatch: expected {ROUND14_SHA256}, got {round14_sha256}"
        ));
    }
    let dependency: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    let curve = CurveParams::p256();
    let (factor_base, columns, _, _, _) = build_indexes(&curve)?;
    let sample_specs = round14_sample_specs();
    verify_round14_samples(&dependency, &sample_specs, &columns)?;
    let mut samples = Vec::new();
    for (kind, indices) in sample_specs {
        samples.push(run_compression_sample(&kind, indices, &columns, &curve)?);
    }
    let result = CompressionExperimentResult {
        schema: "p256.s17_target_indexed_compression/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        round14_sha256,
        factor_base,
        local_maximum_degree: 2,
        samples,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn validate_atom_columns(column_indices: &[u64]) -> Result<(), String> {
    if column_indices.len() != ATOM_COLUMNS {
        return Err(format!(
            "packed atom requires {ATOM_COLUMNS} columns, got {}",
            column_indices.len()
        ));
    }
    let mut seen = BTreeSet::new();
    for &column in column_indices {
        if column >= COLUMNS {
            return Err(format!(
                "packed atom column {column} is outside the {COLUMNS}-column dictionary"
            ));
        }
        if !seen.insert(column) {
            return Err(format!("packed atom repeats column {column}"));
        }
    }
    Ok(())
}

fn pack_atom_columns(column_indices: &[u64]) -> Result<[u8; PACKED_ATOM_BYTES], String> {
    validate_atom_columns(column_indices)?;
    let mut packet = [0u8; PACKED_ATOM_BYTES];
    let mut cursor = 0usize;
    for &column in column_indices {
        for shift in (0..COLUMN_INDEX_BITS).rev() {
            let bit = ((column >> shift) & 1) as u8;
            packet[cursor / 8] |= bit << (7 - cursor % 8);
            cursor += 1;
        }
    }
    if cursor != PACKED_ATOM_BYTES * 8 {
        return Err("packed atom bit count changed".into());
    }
    Ok(packet)
}

fn unpack_atom_columns(packet: &[u8]) -> Result<Vec<u64>, String> {
    if packet.len() != PACKED_ATOM_BYTES {
        return Err(format!(
            "packed atom requires {PACKED_ATOM_BYTES} bytes, got {}",
            packet.len()
        ));
    }
    let mut columns = Vec::with_capacity(ATOM_COLUMNS);
    let mut cursor = 0usize;
    for _ in 0..ATOM_COLUMNS {
        let mut column = 0u64;
        for _ in 0..COLUMN_INDEX_BITS {
            let bit = (packet[cursor / 8] >> (7 - cursor % 8)) & 1;
            column = (column << 1) | u64::from(bit);
            cursor += 1;
        }
        columns.push(column);
    }
    if cursor != packet.len() * 8 {
        return Err("packed atom decoder did not consume the packet".into());
    }
    validate_atom_columns(&columns)?;
    Ok(columns)
}

fn round15_sample<'a>(
    dependency: &'a serde_json::Value,
    kind: &str,
) -> Result<&'a serde_json::Value, String> {
    dependency["samples"]
        .as_array()
        .ok_or("round-15 samples are missing")?
        .iter()
        .find(|sample| sample["kind"].as_str() == Some(kind))
        .ok_or_else(|| format!("round-15 sample {kind} is missing"))
}

fn run_packed_atom_sample(
    kind: &str,
    expected_columns: Vec<u64>,
    dependency: &serde_json::Value,
    columns: &[Column],
    curve: &CurveParams,
) -> Result<PackedAtomSampleResult, String> {
    let stored_sample = round15_sample(dependency, kind)?;
    let stored_columns: Vec<u64> = serde_json::from_value(stored_sample["columns"].clone())
        .map_err(|error| format!("round-15 {kind} columns: {error}"))?;
    if stored_columns != expected_columns {
        return Err(format!("round-15 {kind} column selection changed"));
    }

    let packet = pack_atom_columns(&stored_columns)?;
    let decoded = unpack_atom_columns(&packet)?;
    let decoded_exact = decoded == stored_columns;
    let reencoded_exact = pack_atom_columns(&decoded)? == packet;
    if !decoded_exact || !reencoded_exact {
        return Err(format!("{kind} packed atom failed canonical round trip"));
    }

    let recomputed = run_compression_sample(kind, decoded, columns, curve)?;
    let recomputed_value = serde_json::to_value(&recomputed).map_err(|error| error.to_string())?;
    if &recomputed_value != stored_sample {
        return Err(format!(
            "{kind} reconstructed round-15 sample differs from its frozen receipt"
        ));
    }
    if recomputed.targets.len() != 2 || recomputed.targets.iter().any(|target| !target.exact) {
        return Err(format!("{kind} reconstructed targets are not both exact"));
    }
    let one_cold_query_quadratic_solves = recomputed.targets[0].quadratic_solves_cold;
    if recomputed
        .targets
        .iter()
        .any(|target| target.quadratic_solves_cold != one_cold_query_quadratic_solves)
    {
        return Err(format!("{kind} cold-query solve counts differ"));
    }
    let two_query_one_rebuild_quadratic_solves = recomputed.build_quadratic_solves
        + recomputed
            .targets
            .iter()
            .map(|target| target.quadratic_solves_warm)
            .sum::<u64>();
    if packet.len() != 36
        || recomputed.retained_raw_bytes != 8_192
        || recomputed.materialized_raw_bytes != 1_048_576
        || one_cold_query_quadratic_solves != 280
        || two_query_one_rebuild_quadratic_solves != 408
        || recomputed.build_linear_degeneracies != 0
        || recomputed.build_universal_degeneracies != 0
    {
        return Err(format!("{kind} packed-atom boundary changed"));
    }
    let additional_persistent_reduction =
        recomputed.retained_raw_bytes as f64 / packet.len() as f64;
    if additional_persistent_reduction < 30.0 {
        return Err(format!(
            "{kind} compression gate failed: {additional_persistent_reduction:.6}x"
        ));
    }

    Ok(PackedAtomSampleResult {
        kind: recomputed.kind,
        packed_atom_hex: hex::encode(packet),
        packed_atom_sha256: hex::encode(sha256(&packet)),
        bits_per_column: COLUMN_INDEX_BITS,
        packet_bytes: packet.len(),
        round15_retained_raw_bytes: recomputed.retained_raw_bytes,
        materialized_sixteen_leaf_raw_bytes: recomputed.materialized_raw_bytes,
        additional_persistent_reduction,
        full_materialized_reduction: recomputed.materialized_raw_bytes as f64 / packet.len() as f64,
        columns: recomputed.columns,
        x_coordinates: recomputed.x_coordinates,
        decoded_exact,
        reencoded_exact,
        dependency_sample_exact: true,
        left_image_sha256: recomputed.left_image_sha256,
        right_image_sha256: recomputed.right_image_sha256,
        transient_reconstruction_raw_bytes: recomputed.retained_raw_bytes,
        build_quadratic_solves: recomputed.build_quadratic_solves,
        one_cold_query_quadratic_solves,
        two_query_one_rebuild_quadratic_solves,
        intermediate_reference_group_additions: recomputed.intermediate_reference_group_additions,
        full_reference_group_additions: recomputed.full_reference_group_additions,
        build_linear_degeneracies: recomputed.build_linear_degeneracies,
        build_universal_degeneracies: recomputed.build_universal_degeneracies,
        build_multiplication_counts: recomputed.build_multiplication_counts,
        targets: recomputed.targets,
    })
}

fn run_packed_compression(round15_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round15_path).map_err(|error| error.to_string())?;
    let round15_sha256 = hex::encode(sha256(&bytes));
    if round15_sha256 != ROUND15_SHA256 {
        return Err(format!(
            "round-15 SHA-256 mismatch: expected {ROUND15_SHA256}, got {round15_sha256}"
        ));
    }
    let dependency: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if dependency["schema"].as_str() != Some("p256.s17_target_indexed_compression/v1")
        || dependency["curve"].as_str() != Some(CURVE_SLUG)
        || dependency["factor_base"]["fb_id"].as_str() != Some(FB_ID)
        || dependency["factor_base"]["fb_sha256"].as_str() != Some(FB_SHA256)
        || dependency["factor_base"]["points_sha256"].as_str() != Some(POINTS_SHA256)
    {
        return Err("round-15 identity receipt changed".into());
    }

    let curve = CurveParams::p256();
    let (factor_base, columns, _, _, _) = build_indexes(&curve)?;
    let mut samples = Vec::new();
    for (kind, expected_columns) in round14_sample_specs() {
        samples.push(run_packed_atom_sample(
            &kind,
            expected_columns,
            &dependency,
            &columns,
            &curve,
        )?);
    }
    let result = PackedAtomExperimentResult {
        schema: "p256.s17_packed_atom/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        round15_sha256,
        factor_base,
        local_maximum_degree: 2,
        factor_base_column_bits: COLUMN_INDEX_BITS,
        packed_atom_bytes: PACKED_ATOM_BYTES,
        samples,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn run_batched_inversion_sample(
    kind: &str,
    dependency: &serde_json::Value,
    columns: &[Column],
    curve: &CurveParams,
) -> Result<BatchedInversionSampleResult, String> {
    let stored_sample = dependency["samples"]
        .as_array()
        .ok_or("round-16 samples are missing")?
        .iter()
        .find(|sample| sample["kind"].as_str() == Some(kind))
        .ok_or_else(|| format!("round-16 sample {kind} is missing"))?;
    let packed_atom_hex = stored_sample["packed_atom_hex"]
        .as_str()
        .ok_or_else(|| format!("round-16 {kind} packet is missing"))?;
    let packet =
        hex::decode(packed_atom_hex).map_err(|error| format!("round-16 {kind} packet: {error}"))?;
    let decoded = unpack_atom_columns(&packet)?;
    let stored_columns: Vec<u64> = serde_json::from_value(stored_sample["columns"].clone())
        .map_err(|error| format!("round-16 {kind} columns: {error}"))?;
    let decoded_exact =
        decoded == stored_columns && hex::encode(pack_atom_columns(&decoded)?) == packed_atom_hex;
    if !decoded_exact {
        return Err(format!("round-16 {kind} packet failed exact decode"));
    }

    let scalar = run_compression_sample(kind, decoded.clone(), columns, curve)?;
    let scalar_targets =
        serde_json::to_value(&scalar.targets).map_err(|error| error.to_string())?;
    if stored_sample["columns"]
        != serde_json::to_value(&scalar.columns).map_err(|error| error.to_string())?
        || stored_sample["x_coordinates"]
            != serde_json::to_value(&scalar.x_coordinates).map_err(|error| error.to_string())?
        || stored_sample["left_image_sha256"].as_str() != Some(&scalar.left_image_sha256)
        || stored_sample["right_image_sha256"].as_str() != Some(&scalar.right_image_sha256)
        || stored_sample["targets"] != scalar_targets
    {
        return Err(format!(
            "round-16 {kind} receipt differs from scalar replay"
        ));
    }

    let selected: Vec<&Column> = decoded
        .iter()
        .map(|&index| {
            columns
                .get(index as usize)
                .ok_or_else(|| format!("{kind} column {index} is out of range"))
        })
        .collect::<Result<_, _>>()?;
    let xs: Vec<Fe> = selected.iter().map(|column| column.x).collect();
    let points: Vec<Point> = selected.iter().map(|column| column.low.clone()).collect();
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut build_solves = 0u64;
    let mut build_roots = 0u64;
    let mut build_linear = 0u64;
    let mut build_universal = 0u64;
    let mut build_counts = MultiplicationCounts::default();
    let mut build_profile = BatchInversionProfile::default();
    let mut intermediate_reference_group_additions = 0u64;
    let left = build_verified_eight_image_batched(
        &format!("{kind}/left/batched"),
        &xs[..8],
        &points[..8],
        curve,
        a,
        b,
        &sqrt_exponent,
        &inverse_exponent,
        &mut build_solves,
        &mut build_roots,
        &mut build_linear,
        &mut build_universal,
        &mut build_counts,
        &mut build_profile,
        &mut intermediate_reference_group_additions,
    )?;
    let right = build_verified_eight_image_batched(
        &format!("{kind}/right/batched"),
        &xs[8..],
        &points[8..],
        curve,
        a,
        b,
        &sqrt_exponent,
        &inverse_exponent,
        &mut build_solves,
        &mut build_roots,
        &mut build_linear,
        &mut build_universal,
        &mut build_counts,
        &mut build_profile,
        &mut intermediate_reference_group_additions,
    )?;
    build_counts.finish();
    if image_digest(&left) != scalar.left_image_sha256
        || image_digest(&right) != scalar.right_image_sha256
        || build_solves != scalar.build_quadratic_solves
        || build_roots != 304
        || build_linear != 0
        || build_universal != 0
        || build_counts.inversions != 5_802
        || build_counts.total != 52_010
        || build_profile.batch_sizes != [4, 4, 64, 4, 4, 64]
        || build_profile.batch_inversions != 6
        || build_profile.scalar_fallbacks != 8
    {
        return Err(format!("{kind} batched build boundary changed"));
    }
    let build_multiplication_reduction =
        scalar.build_multiplication_counts.total as f64 / build_counts.total as f64;
    if build_multiplication_reduction < 2.0 {
        return Err(format!("{kind} build multiplication gate failed"));
    }

    let curve_a = curve.a_fe();
    let planted = points.iter().fold(Point::Infinity, |sum, point| {
        sum.add_vartime(point, &curve_a)
    });
    let public_preimage = format!("{CURVE_SLUG}/s17-target-join-round15/sample/{kind}/public");
    let mut public_scalar = BigUint::from_bytes_be(&sha256(public_preimage.as_bytes())) % &curve.n;
    if public_scalar.is_zero() {
        public_scalar = BigUint::one();
    }
    let public = curve
        .generator()
        .scalar_mul_vartime(&public_scalar, &curve_a);
    let targets = vec![
        query_batched_target(
            "planted-positive",
            None,
            &planted,
            scalar.targets[0].reference_positive,
            &left,
            &right,
            curve,
            &scalar.targets[0],
            scalar.build_multiplication_counts,
            build_counts,
        )?,
        query_batched_target(
            "hash-public",
            Some(&public_scalar),
            &public,
            scalar.targets[1].reference_positive,
            &left,
            &right,
            curve,
            &scalar.targets[1],
            scalar.build_multiplication_counts,
            build_counts,
        )?,
    ];
    let two_target_scalar_multiplications = scalar.build_multiplication_counts.total
        + scalar
            .targets
            .iter()
            .map(|target| target.multiplication_counts_warm.total)
            .sum::<u64>();
    let two_target_batched_multiplications = build_counts.total
        + targets
            .iter()
            .map(|target| target.batched_multiplication_counts_warm.total)
            .sum::<u64>();
    let two_target_multiplication_reduction =
        two_target_scalar_multiplications as f64 / two_target_batched_multiplications as f64;
    let exact = targets.iter().all(|target| target.exact)
        && two_target_scalar_multiplications == 280_704
        && two_target_batched_multiplications == 131_368
        && two_target_multiplication_reduction >= 2.0;
    if !exact {
        return Err(format!("{kind} batched two-target gate failed"));
    }

    Ok(BatchedInversionSampleResult {
        kind: kind.into(),
        packed_atom_hex: packed_atom_hex.into(),
        columns: decoded,
        decoded_exact,
        left_image_sha256: image_digest(&left),
        right_image_sha256: image_digest(&right),
        scalar_build_multiplication_counts: scalar.build_multiplication_counts,
        batched_build_multiplication_counts: build_counts,
        build_multiplication_reduction,
        scalar_build_quadratic_solves: scalar.build_quadratic_solves,
        batched_build_quadratic_solves: build_solves,
        build_inversion_profile: build_profile,
        two_target_scalar_multiplications,
        two_target_batched_multiplications,
        two_target_multiplication_reduction,
        transient_reconstruction_raw_bytes: scalar.retained_raw_bytes,
        intermediate_reference_group_additions,
        full_reference_group_additions: scalar.full_reference_group_additions,
        targets,
        exact,
    })
}

fn run_batched_inversion(round16_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round16_path).map_err(|error| error.to_string())?;
    let round16_sha256 = hex::encode(sha256(&bytes));
    if round16_sha256 != ROUND16_SHA256 {
        return Err(format!(
            "round-16 SHA-256 mismatch: expected {ROUND16_SHA256}, got {round16_sha256}"
        ));
    }
    let dependency: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if dependency["schema"].as_str() != Some("p256.s17_packed_atom/v1")
        || dependency["curve"].as_str() != Some(CURVE_SLUG)
        || dependency["factor_base"]["fb_id"].as_str() != Some(FB_ID)
        || dependency["factor_base"]["fb_sha256"].as_str() != Some(FB_SHA256)
        || dependency["factor_base"]["points_sha256"].as_str() != Some(POINTS_SHA256)
    {
        return Err("round-16 identity receipt changed".into());
    }
    let curve = CurveParams::p256();
    let (factor_base, columns, _, _, _) = build_indexes(&curve)?;
    let mut samples = Vec::new();
    for (kind, _) in round14_sample_specs() {
        samples.push(run_batched_inversion_sample(
            &kind,
            &dependency,
            &columns,
            &curve,
        )?);
    }
    let result = BatchedInversionExperimentResult {
        schema: "p256.s17_batched_inversion/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        round16_sha256,
        factor_base,
        local_maximum_degree: 2,
        modeled_relation_success_unchanged: 0.964_000_18,
        end_to_end_fraction_in_oracle: None,
        end_to_end_speedup: None,
        samples,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn signed_atom_reference(
    points: &[Point],
    curve: &CurveParams,
) -> Result<(BTreeMap<[u8; 65], Vec<u32>>, u64), String> {
    if points.len() != ATOM_COLUMNS {
        return Err(format!("signed atom reference needs {ATOM_COLUMNS} points"));
    }
    let curve_a = curve.a_fe();
    let mut states = BTreeMap::from([(
        full_point_key(&Point::Infinity),
        (Point::Infinity, vec![0u32]),
    )]);
    let mut additions = 0u64;
    for (index, point) in points.iter().enumerate() {
        let negative = negate(point, curve);
        let mut next: BTreeMap<[u8; 65], (Point, Vec<u32>)> = BTreeMap::new();
        for (_, (sum, masks)) in states {
            for (negative_sign, signed) in [(false, point), (true, &negative)] {
                for &mask in &masks {
                    let value = sum.add_vartime(signed, &curve_a);
                    additions += 1;
                    let next_mask = if negative_sign {
                        mask | (1u32 << index)
                    } else {
                        mask
                    };
                    let key = full_point_key(&value);
                    let entry = next.entry(key).or_insert_with(|| (value, Vec::new()));
                    entry.1.push(next_mask);
                }
            }
        }
        states = next;
    }
    let reference: BTreeMap<[u8; 65], Vec<u32>> = states
        .into_iter()
        .map(|(key, (_, masks))| (key, masks))
        .collect();
    let signed_sums: usize = reference.values().map(Vec::len).sum();
    if signed_sums != 1usize << ATOM_COLUMNS {
        return Err(format!(
            "signed atom reference retained {signed_sums} masks, expected {}",
            1usize << ATOM_COLUMNS
        ));
    }
    Ok((reference, additions))
}

fn verify_outer_relation(
    mask: u32,
    atom_points: &[Point],
    outer_point: &Point,
    outer_negative: bool,
    target: &Point,
    curve: &CurveParams,
) -> (bool, u64) {
    let curve_a = curve.a_fe();
    let mut sum = Point::Infinity;
    let mut additions = 0u64;
    for (index, point) in atom_points.iter().enumerate() {
        let signed = if mask & (1u32 << index) == 0 {
            point.clone()
        } else {
            negate(point, curve)
        };
        sum = sum.add_vartime(&signed, &curve_a);
        additions += 1;
    }
    let signed_outer = if outer_negative {
        negate(outer_point, curve)
    } else {
        outer_point.clone()
    };
    sum = sum.add_vartime(&signed_outer, &curve_a);
    additions += 1;
    (sum == *target, additions)
}

fn relation_digest(relations: &[OuterRelation]) -> Result<String, String> {
    let bytes = serde_json::to_vec(relations).map_err(|error| error.to_string())?;
    Ok(hex::encode(sha256(&bytes)))
}

fn outer_target_checkpoint(
    kind: &str,
    accumulator: &OuterTargetAccumulator,
) -> Result<OuterTargetCheckpoint, String> {
    Ok(OuterTargetCheckpoint {
        kind: kind.into(),
        signed_branch_cells: accumulator.signed_branch_cells,
        oracle_calls: accumulator.oracle_calls,
        quadratic_solves: accumulator.quadratic_solves,
        quadratic_roots_returned: accumulator.quadratic_roots_returned,
        algebra_index_lookups: accumulator.algebra_index_lookups,
        algebra_hits: accumulator.algebra_hits,
        candidate_positive_trials: accumulator.candidate_positive_trials,
        reference_positive_trials: accumulator.reference_positive_trials,
        false_negatives: accumulator.false_negatives,
        false_positives: accumulator.false_positives,
        multiplication_counts: accumulator.multiplication_counts,
        target_adjustment_group_additions: accumulator.target_adjustment_group_additions,
        reference_index_lookups: accumulator.reference_index_lookups,
        direct_relation_verification_additions: accumulator.direct_relation_verification_additions,
        batch_profile: accumulator.batch_profile.clone(),
        relation_count: accumulator.relations.len(),
        relation_sha256: relation_digest(&accumulator.relations)?,
        trial_stream_sha256: hex::encode(sha256(accumulator.trial_stream.as_bytes())),
        exact: accumulator.false_negatives == 0 && accumulator.false_positives == 0,
    })
}

#[allow(clippy::too_many_arguments)]
fn run_outer_trial(
    target: &OuterTargetSpec,
    stream_position: u64,
    outer_column: u64,
    outer_negative: bool,
    outer_point: &Point,
    atom_points: &[Point],
    reference: &BTreeMap<[u8; 65], Vec<u32>>,
    left: &P256ImageState,
    right: &P256ImageState,
    curve: &CurveParams,
    a: Fe,
    b: Fe,
    sqrt_exponent: &BigUint,
    inverse_exponent: &BigUint,
    accumulator: &mut OuterTargetAccumulator,
) -> Result<(), String> {
    let curve_a = curve.a_fe();
    let signed_outer = if outer_negative {
        negate(outer_point, curve)
    } else {
        outer_point.clone()
    };
    let adjusted = target
        .point
        .add_vartime(&negate(&signed_outer, curve), &curve_a);
    accumulator.target_adjustment_group_additions += 1;
    if matches!(adjusted, Point::Infinity) {
        return Err(format!(
            "{} outer trial {stream_position}/{outer_column} adjusted to infinity",
            target.kind
        ));
    }
    let reference_masks = reference.get(&full_point_key(&adjusted));
    let reference_positive = reference_masks.is_some_and(|masks| !masks.is_empty());
    accumulator.reference_index_lookups += 1;
    let observation = observe_batched_target(
        &adjusted,
        left,
        right,
        curve,
        a,
        b,
        sqrt_exponent,
        inverse_exponent,
    )?;
    accumulator.signed_branch_cells += 1;
    accumulator.oracle_calls += 1;
    accumulator.quadratic_solves += observation.quadratic_solves;
    accumulator.quadratic_roots_returned += observation.quadratic_roots_returned;
    accumulator.algebra_index_lookups += observation.algebra_index_lookups;
    accumulator.algebra_hits += observation.algebra_hits;
    accumulator.candidate_positive_trials += u64::from(observation.positive);
    accumulator.reference_positive_trials += u64::from(reference_positive);
    accumulator.false_negatives += u64::from(reference_positive && !observation.positive);
    accumulator.false_positives += u64::from(!reference_positive && observation.positive);
    accumulator
        .multiplication_counts
        .add_assign(observation.multiplication_counts);
    accumulator
        .batch_profile
        .record(&observation.inversion_profile);
    writeln!(
        accumulator.trial_stream,
        "{},{stream_position},{outer_column},{},{},{},{},{},{},{}",
        target.kind,
        u8::from(outer_negative),
        u8::from(reference_positive),
        u8::from(observation.positive),
        observation.algebra_hits,
        observation.quadratic_roots_returned,
        observation.multiplication_counts.total,
        observation.hit_sha256
    )
    .expect("writing a String cannot fail");
    if observation.linear_degeneracies != 0 || observation.universal_degeneracies != 0 {
        return Err(format!(
            "{} outer trial produced {} linear and {} universal degeneracies",
            target.kind, observation.linear_degeneracies, observation.universal_degeneracies
        ));
    }
    if reference_positive != observation.positive {
        return Err(format!(
            "{} outer trial {stream_position}/{outer_column}/{} disagreed: reference={reference_positive}, candidate={}",
            target.kind,
            u8::from(outer_negative),
            observation.positive
        ));
    }
    if let Some(masks) = reference_masks {
        for &mask in masks {
            let (direct_verified, additions) = verify_outer_relation(
                mask,
                atom_points,
                outer_point,
                outer_negative,
                &target.point,
                curve,
            );
            accumulator.direct_relation_verification_additions += additions;
            if !direct_verified {
                return Err(format!(
                    "{} outer relation failed direct group verification",
                    target.kind
                ));
            }
            accumulator.relations.push(OuterRelation {
                target_kind: target.kind.clone(),
                stream_position,
                outer_column,
                outer_negative,
                atom_negative_mask: mask,
                direct_verified,
            });
        }
    }
    Ok(())
}

fn log2_slope(depths: &[u32], values: &[u64]) -> Result<f64, String> {
    if depths.len() != values.len() || depths.len() < 2 || values.contains(&0) {
        return Err("growth fit needs matching nonzero samples".into());
    }
    let count = depths.len() as f64;
    let mean_x = depths.iter().map(|&value| value as f64).sum::<f64>() / count;
    let logs: Vec<f64> = values.iter().map(|&value| (value as f64).log2()).collect();
    let mean_y = logs.iter().sum::<f64>() / count;
    let numerator: f64 = depths
        .iter()
        .map(|&value| value as f64)
        .zip(&logs)
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum();
    let denominator: f64 = depths
        .iter()
        .map(|&value| {
            let delta = value as f64 - mean_x;
            delta * delta
        })
        .sum();
    Ok(numerator / denominator)
}

fn binomial_big(n: u64, k: u32) -> BigUint {
    let k = u64::from(k).min(n - u64::from(k));
    let mut value = BigUint::one();
    for index in 0..k {
        value *= BigUint::from(n - index);
        value /= BigUint::from(index + 1);
    }
    value
}

fn load_residual_degree_evidence(path: &PathBuf) -> Result<ResidualDegreeEvidence, String> {
    let bytes = std::fs::read(path).map_err(|error| error.to_string())?;
    let digest = hex::encode(sha256(&bytes));
    if digest != ROUND6_SHA256 {
        return Err(format!(
            "round-6 SHA-256 mismatch: expected {ROUND6_SHA256}, got {digest}"
        ));
    }
    let dependency: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if dependency["schema"].as_str() != Some("p256.dickson_residual_scaling/v1")
        || dependency["curve"].as_str() != Some(CURVE_SLUG)
    {
        return Err("round-6 degree identity changed".into());
    }
    let cells = dependency["cells"]
        .as_array()
        .ok_or("round-6 cells are missing")?;
    let mut maximum_by_residual_depth = BTreeMap::new();
    let mut every_component_complete_and_correct = true;
    for cell in cells {
        let residual = cell["residual_depth"]
            .as_u64()
            .ok_or("round-6 cell residual depth is missing")? as u32;
        let degree = cell["max_solving_degree"]
            .as_u64()
            .ok_or("round-6 cell maximum degree is missing")? as u32;
        maximum_by_residual_depth
            .entry(residual)
            .and_modify(|value: &mut u32| *value = (*value).max(degree))
            .or_insert(degree);
        let components = cell["components"]
            .as_u64()
            .ok_or("round-6 component count is missing")?;
        every_component_complete_and_correct &= cell["complete_components"].as_u64()
            == Some(components)
            && cell["correct_components"].as_u64() == Some(components);
    }
    let expected = BTreeMap::from([(1u32, 3u32), (2, 3), (3, 4)]);
    if maximum_by_residual_depth != expected || !every_component_complete_and_correct {
        return Err(format!(
            "round-6 residual-degree boundary changed: {maximum_by_residual_depth:?}"
        ));
    }
    Ok(ResidualDegreeEvidence {
        round6_sha256: digest,
        maximum_by_residual_depth,
        every_component_complete_and_correct,
        promotion_degree_gate_at_most_five: true,
    })
}

fn run_outer_scan(
    round17_path: &PathBuf,
    round6_path: &PathBuf,
    out: Option<PathBuf>,
) -> Result<(), String> {
    let round17_bytes = std::fs::read(round17_path).map_err(|error| error.to_string())?;
    let round17_sha256 = hex::encode(sha256(&round17_bytes));
    if round17_sha256 != ROUND17_SHA256 {
        return Err(format!(
            "round-17 SHA-256 mismatch: expected {ROUND17_SHA256}, got {round17_sha256}"
        ));
    }
    let dependency: serde_json::Value =
        serde_json::from_slice(&round17_bytes).map_err(|error| error.to_string())?;
    if dependency["schema"].as_str() != Some("p256.s17_batched_inversion/v1")
        || dependency["curve"].as_str() != Some(CURVE_SLUG)
        || dependency["factor_base"]["fb_id"].as_str() != Some(FB_ID)
        || dependency["factor_base"]["fb_sha256"].as_str() != Some(FB_SHA256)
        || dependency["factor_base"]["points_sha256"].as_str() != Some(POINTS_SHA256)
    {
        return Err("round-17 identity receipt changed".into());
    }
    let stored_sample = dependency["samples"]
        .as_array()
        .ok_or("round-17 samples are missing")?
        .iter()
        .find(|sample| sample["kind"].as_str() == Some("hash-0"))
        .ok_or("round-17 hash-0 sample is missing")?;
    let packed_atom_hex = stored_sample["packed_atom_hex"]
        .as_str()
        .ok_or("round-17 hash-0 packet is missing")?;
    let packet = hex::decode(packed_atom_hex).map_err(|error| error.to_string())?;
    let atom_columns = unpack_atom_columns(&packet)?;
    let expected_atom_columns = vec![
        114_616, 101_224, 36_382, 129_678, 17_773, 72_570, 33_360, 118_255, 88_179, 58_271, 40_458,
        117_426, 13_544, 42_490, 88_617, 31_941,
    ];
    if atom_columns != expected_atom_columns
        || packed_atom_hex
            != "6fee18b6823879fa8e115b51b7a20941cdef561cce39f27829cab20d3a0a5fa568a47cc5"
    {
        return Err("round-17 hash-0 packet changed".into());
    }
    let residual_degree_evidence = load_residual_degree_evidence(round6_path)?;

    let curve = CurveParams::p256();
    let (factor_base, columns, _, _, _) = build_indexes(&curve)?;
    let scalar = run_compression_sample("hash-0", atom_columns.clone(), &columns, &curve)?;
    if stored_sample["left_image_sha256"].as_str() != Some(&scalar.left_image_sha256)
        || stored_sample["right_image_sha256"].as_str() != Some(&scalar.right_image_sha256)
    {
        return Err("round-17 hash-0 scalar image receipt changed".into());
    }
    let selected: Vec<&Column> = atom_columns
        .iter()
        .map(|&index| {
            columns
                .get(index as usize)
                .ok_or_else(|| format!("atom column {index} is out of range"))
        })
        .collect::<Result<_, _>>()?;
    let xs: Vec<Fe> = selected.iter().map(|column| column.x).collect();
    let atom_points: Vec<Point> = selected.iter().map(|column| column.low.clone()).collect();
    let a = Fe::from_biguint(&curve.a);
    let b = Fe::from_biguint(&curve.b);
    let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
    let inverse_exponent = &curve.p - BigUint::from(2u8);
    let mut build_solves = 0u64;
    let mut build_roots = 0u64;
    let mut build_linear = 0u64;
    let mut build_universal = 0u64;
    let mut build_counts = MultiplicationCounts::default();
    let mut build_profile = BatchInversionProfile::default();
    let mut intermediate_reference_group_additions = 0u64;
    let left = build_verified_eight_image_batched(
        "outer/hash-0/left",
        &xs[..8],
        &atom_points[..8],
        &curve,
        a,
        b,
        &sqrt_exponent,
        &inverse_exponent,
        &mut build_solves,
        &mut build_roots,
        &mut build_linear,
        &mut build_universal,
        &mut build_counts,
        &mut build_profile,
        &mut intermediate_reference_group_additions,
    )?;
    let right = build_verified_eight_image_batched(
        "outer/hash-0/right",
        &xs[8..],
        &atom_points[8..],
        &curve,
        a,
        b,
        &sqrt_exponent,
        &inverse_exponent,
        &mut build_solves,
        &mut build_roots,
        &mut build_linear,
        &mut build_universal,
        &mut build_counts,
        &mut build_profile,
        &mut intermediate_reference_group_additions,
    )?;
    build_counts.finish();
    if image_digest(&left) != scalar.left_image_sha256
        || image_digest(&right) != scalar.right_image_sha256
        || build_solves != 152
        || build_roots != 304
        || build_linear != 0
        || build_universal != 0
        || build_counts.total != 52_010
        || build_profile.batch_sizes != [4, 4, 64, 4, 4, 64]
        || build_profile.batch_inversions != 6
        || build_profile.scalar_fallbacks != 8
    {
        return Err("round-18 batched atom build boundary changed".into());
    }

    let (reference, reference_build_group_additions) = signed_atom_reference(&atom_points, &curve)?;
    let reference_signed_sums: u64 = reference.values().map(|masks| masks.len() as u64).sum();
    let curve_a = curve.a_fe();
    let atom_low_sum = atom_points.iter().fold(Point::Infinity, |sum, point| {
        sum.add_vartime(point, &curve_a)
    });
    let planted_outer = columns
        .get(PLANTED_OUTER_COLUMN as usize)
        .ok_or("planted outer column is out of range")?;
    let planted_target = atom_low_sum.add_vartime(&planted_outer.low, &curve_a);
    if matches!(planted_target, Point::Infinity) || !curve.is_on_curve(&planted_target) {
        return Err("planted outer target is infinity or off curve".into());
    }
    let mut public_scalar =
        BigUint::from_bytes_be(&sha256(OUTER_PUBLIC_TARGET_PREIMAGE.as_bytes())) % &curve.n;
    if public_scalar.is_zero() {
        public_scalar = BigUint::one();
    }
    let public_target = curve
        .generator()
        .scalar_mul_vartime(&public_scalar, &curve_a);
    if matches!(public_target, Point::Infinity) || !curve.is_on_curve(&public_target) {
        return Err("public outer target is infinity or off curve".into());
    }
    let targets = [
        OuterTargetSpec {
            kind: "planted-control".into(),
            scalar: None,
            point: planted_target,
            planted: true,
        },
        OuterTargetSpec {
            kind: "hash-public".into(),
            scalar: Some(public_scalar),
            point: public_target,
            planted: false,
        },
    ];
    let target_identities: Vec<OuterTargetIdentity> = targets
        .iter()
        .map(|target| {
            let (x, y) = affine_coordinates(&target.point).expect("targets were checked affine");
            OuterTargetIdentity {
                kind: target.kind.clone(),
                scalar: target.scalar.as_ref().map(lower_hex),
                target: [lower_hex(x), lower_hex(y)],
                planted: target.planted,
            }
        })
        .collect();

    let atom_set: BTreeSet<u64> = atom_columns.iter().copied().collect();
    let mut accumulators = vec![OuterTargetAccumulator::default(); targets.len()];
    let mut checkpoints = Vec::new();
    let mut seen_outer = BTreeSet::new();
    let mut skipped_atom_columns = 0u64;
    let mut eligible_outer_columns = 0u64;
    for stream_position in 0..(1u64 << OUTER_MAX_DEPTH) {
        let outer_column = (OUTER_OFFSET + OUTER_STRIDE * stream_position) % COLUMNS;
        if !seen_outer.insert(outer_column) {
            return Err(format!(
                "outer stream repeated column {outer_column} at position {stream_position}"
            ));
        }
        if stream_position == PLANTED_OUTER_POSITION && outer_column != PLANTED_OUTER_COLUMN {
            return Err("planted outer stream position changed".into());
        }
        if atom_set.contains(&outer_column) {
            skipped_atom_columns += 1;
        } else {
            eligible_outer_columns += 1;
            let outer_point = &columns[outer_column as usize].low;
            for (target, accumulator) in targets.iter().zip(&mut accumulators) {
                for outer_negative in [false, true] {
                    run_outer_trial(
                        target,
                        stream_position,
                        outer_column,
                        outer_negative,
                        outer_point,
                        &atom_points,
                        &reference,
                        &left,
                        &right,
                        &curve,
                        a,
                        b,
                        &sqrt_exponent,
                        &inverse_exponent,
                        accumulator,
                    )?;
                }
            }
        }
        let positions = stream_position + 1;
        if let Some(&depth) = OUTER_DEPTHS
            .iter()
            .find(|&&depth| positions == 1u64 << depth)
        {
            let target_rows = targets
                .iter()
                .zip(&accumulators)
                .map(|(target, accumulator)| outer_target_checkpoint(&target.kind, accumulator))
                .collect::<Result<_, _>>()?;
            checkpoints.push(OuterCheckpoint {
                depth,
                stream_positions: positions,
                skipped_atom_columns,
                eligible_outer_columns,
                targets: target_rows,
            });
        }
    }
    if checkpoints.len() != OUTER_DEPTHS.len() {
        return Err("outer scan missed a frozen checkpoint".into());
    }
    let public_field_multiplications: Vec<u64> = checkpoints
        .iter()
        .map(|checkpoint| build_counts.total + checkpoint.targets[1].multiplication_counts.total)
        .collect();
    let public_oracle_calls: Vec<u64> = checkpoints
        .iter()
        .map(|checkpoint| checkpoint.targets[1].oracle_calls)
        .collect();
    let growth_fit = OuterGrowthFit {
        depths: OUTER_DEPTHS.to_vec(),
        public_field_multiplications: public_field_multiplications.clone(),
        log2_multiplication_slope_per_depth: log2_slope(
            &OUTER_DEPTHS,
            &public_field_multiplications,
        )?,
        outer_calls_growth_per_depth: log2_slope(&OUTER_DEPTHS, &public_oracle_calls)?,
    };

    let final_public = &checkpoints
        .last()
        .expect("checkpoints are nonempty")
        .targets[1];
    let multiplications_per_oracle_call =
        final_public.multiplication_counts.total as f64 / final_public.oracle_calls as f64;
    let eligible_full = COLUMNS - ATOM_COLUMNS as u64;
    let full_calls = 2 * eligible_full;
    let signed_domain = BigUint::from(eligible_full) << 17usize;
    let distinct_sixteen_column_atoms = binomial_big(COLUMNS, 16);
    let distinct_seventeen_column_sets = binomial_big(COLUMNS, 17);
    let whole_factor_base_signed_domain = &distinct_seventeen_column_sets << 17usize;
    let group_order = curve
        .n
        .to_f64()
        .ok_or("P-256 subgroup order cannot be represented as f64")?;
    let poisson_mean_per_atom = signed_domain
        .to_f64()
        .ok_or("outer signed domain cannot be represented as f64")?
        / group_order;
    let poisson_success_per_atom = -(-poisson_mean_per_atom).exp_m1();
    let expected_atom_scans_per_relation = 1.0 / poisson_success_per_atom;
    let complete_atom_scan_multiplications =
        build_counts.total as f64 + full_calls as f64 * multiplications_per_oracle_call;
    let relation_rows = (COLUMNS * 105).div_ceil(100);
    let whole_factor_base_poisson_mean = whole_factor_base_signed_domain
        .to_f64()
        .ok_or("whole factor-base signed domain cannot be represented as f64")?
        / group_order;
    let whole_factor_base_poisson_success = -(-whole_factor_base_poisson_mean).exp_m1();
    let projected_targets = relation_rows as f64 / whole_factor_base_poisson_success;
    let projected_per_relation =
        complete_atom_scan_multiplications * expected_atom_scans_per_relation;
    let projected_collection = projected_targets * projected_per_relation;
    let collection_log2 = projected_collection.log2();
    let below_2_pow_120 = collection_log2 < 120.0;
    let below_2_pow_128 = collection_log2 < 128.0;
    let classification = if below_2_pow_120 {
        "promotion candidate"
    } else if below_2_pow_128 {
        "not promoted: less than eight bits of margin"
    } else {
        "rejected: enumeration exponent exceeds rho"
    };
    let collection_projection = OuterCollectionProjection {
        classification: classification.into(),
        distinct_sixteen_column_atoms: distinct_sixteen_column_atoms.to_string(),
        distinct_sixteen_column_atoms_log2: distinct_sixteen_column_atoms
            .to_f64()
            .ok_or("sixteen-column atom count cannot be represented as f64")?
            .log2(),
        distinct_seventeen_column_sets: distinct_seventeen_column_sets.to_string(),
        distinct_seventeen_column_sets_log2: distinct_seventeen_column_sets
            .to_f64()
            .ok_or("seventeen-column set count cannot be represented as f64")?
            .log2(),
        eligible_outer_columns_per_atom: eligible_full,
        signed_oracle_calls_per_complete_atom_scan: full_calls,
        signed_candidate_domain_per_atom: signed_domain.to_string(),
        poisson_mean_per_atom,
        poisson_success_per_atom,
        expected_atom_scans_per_relation,
        expected_atom_scans_log2: expected_atom_scans_per_relation.log2(),
        measured_public_multiplications_per_oracle_call: multiplications_per_oracle_call,
        projected_complete_atom_scan_multiplications: complete_atom_scan_multiplications,
        projected_complete_atom_scan_multiplications_log2: complete_atom_scan_multiplications
            .log2(),
        projected_multiplications_per_relation: projected_per_relation,
        projected_multiplications_per_relation_log2: projected_per_relation.log2(),
        whole_factor_base_signed_domain: whole_factor_base_signed_domain.to_string(),
        whole_factor_base_poisson_mean,
        whole_factor_base_poisson_success,
        relation_rows,
        projected_targets_for_relation_rows: projected_targets,
        projected_collection_multiplications: projected_collection,
        projected_collection_multiplications_log2: collection_log2,
        below_2_pow_120,
        below_2_pow_128,
        bits_above_2_pow_128: collection_log2 - 128.0,
        optimistic_lower_projection: true,
    };

    let row_weight = 17u64;
    let nonzeros = relation_rows * row_weight;
    let csr_entry_bytes = nonzeros * 5;
    let csr_row_offset_bytes = (relation_rows + 1) * 8;
    let csr_total_bytes = csr_entry_bytes + csr_row_offset_bytes;
    let wiedemann_sparse_matvecs = 2 * COLUMNS;
    let wiedemann_nonzero_additions = wiedemann_sparse_matvecs * nonzeros;
    let berlekamp_massey_scalar_operations = COLUMNS * COLUMNS;
    let field_vectors = 3u64;
    let field_vector_bytes = field_vectors * COLUMNS * 32;
    let sparse_linear_algebra_projection = SparseLinearAlgebraProjection {
        rows: relation_rows,
        columns: COLUMNS,
        row_weight,
        nonzeros,
        csr_entry_bytes,
        csr_row_offset_bytes,
        csr_total_bytes,
        wiedemann_sparse_matvecs,
        wiedemann_nonzero_additions,
        berlekamp_massey_scalar_operations,
        field_vectors,
        field_vector_bytes,
        modeled_working_bytes: csr_total_bytes + field_vector_bytes,
        operation_unit_separate_from_collection: true,
    };

    let relations: Vec<OuterRelation> = accumulators
        .iter()
        .flat_map(|accumulator| accumulator.relations.iter().cloned())
        .collect();
    let planted_witness = relations.iter().any(|relation| {
        relation.target_kind == "planted-control"
            && relation.stream_position == PLANTED_OUTER_POSITION
            && relation.outer_column == PLANTED_OUTER_COLUMN
            && !relation.outer_negative
            && relation.atom_negative_mask == 0
            && relation.direct_verified
    });
    let exact = planted_witness
        && checkpoints
            .iter()
            .flat_map(|checkpoint| &checkpoint.targets)
            .all(|target| target.exact)
        && residual_degree_evidence.every_component_complete_and_correct;
    if !exact {
        return Err("round-18 exactness or planted-witness gate failed".into());
    }
    let attack_promotion_gate = exact
        && residual_degree_evidence.promotion_degree_gate_at_most_five
        && collection_projection.below_2_pow_120;
    let memory = OuterMemoryModel {
        persistent_atom_bytes: PACKED_ATOM_BYTES as u64,
        reconstructed_image_raw_bytes: ((left.affine.len() + right.affine.len()) * 32) as u64,
        candidate_peak_logical_bytes: 75_008,
        factor_base_point_raw_bytes: COLUMNS * 96,
        reference_entries: reference_signed_sums,
        reference_final_logical_bytes: reference_signed_sums * (65 + 4),
        reference_build_peak_logical_bytes: ((1u64 << 15) + (1u64 << 16)) * (65 + 64 + 4),
        reference_memory_excluded_from_candidate: true,
    };
    let result = OuterScanExperimentResult {
        schema: "p256.s17_outer_scan/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        round17_sha256,
        factor_base,
        atom_kind: "hash-0".into(),
        packed_atom_hex: packed_atom_hex.into(),
        atom_columns,
        outer_offset: OUTER_OFFSET,
        outer_stride: OUTER_STRIDE,
        checkpoint_depths: OUTER_DEPTHS.to_vec(),
        planted_outer_position: PLANTED_OUTER_POSITION,
        planted_outer_column: PLANTED_OUTER_COLUMN,
        targets: target_identities,
        scalar_build_multiplication_counts: scalar.build_multiplication_counts,
        batched_build_multiplication_counts: build_counts,
        build_inversion_profile: build_profile,
        reference_signed_sums,
        reference_build_group_additions,
        intermediate_reference_group_additions,
        checkpoints,
        growth_fit,
        residual_degree_evidence,
        local_maximum_degree: 2,
        memory,
        collection_projection,
        sparse_linear_algebra_projection,
        relations,
        exact,
        attack_promotion_gate,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn run_transfer(out: Option<PathBuf>) -> Result<(), String> {
    let curve = CurveParams::p256();
    let (factor_base, columns, signed_rows, x_index, signed_index) = build_indexes(&curve)?;
    let a = curve.a_fe();
    let planted = columns[0].low.add_vartime(&columns[1].low, &a);
    if matches!(planted, Point::Infinity) || !curve.is_on_curve(&planted) {
        return Err("planted target is infinity or off curve".into());
    }

    let mut hash_scalar = BigUint::from_bytes_be(&sha256(TARGET_PREIMAGE.as_bytes())) % &curve.n;
    if hash_scalar.is_zero() {
        hash_scalar = BigUint::one();
    }
    let hash_target = curve
        .generator()
        .scalar_mul_vartime(&hash_scalar, &curve.a_fe());

    let targets = vec![
        run_target(
            "planted-positive",
            None,
            &planted,
            &columns,
            &signed_rows,
            &x_index,
            &signed_index,
            &curve,
        )?,
        run_target(
            "hash-public",
            Some(&hash_scalar),
            &hash_target,
            &columns,
            &signed_rows,
            &x_index,
            &signed_index,
            &curve,
        )?,
    ];
    if targets.iter().any(|target| !target.exact) {
        return Err("algebraic and signed-point pair sets differ".into());
    }
    let result = ExperimentResult {
        schema: "p256.s3_image_transfer/v1".into(),
        curve: CURVE_SLUG.into(),
        field_prime: lower_hex(&curve.p),
        curve_a: lower_hex(&curve.a),
        curve_b: lower_hex(&curve.b),
        factor_base,
        target_preimage: TARGET_PREIMAGE.into(),
        targets,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn run(cli: Cli) -> Result<(), String> {
    let Cli {
        out,
        width_round13,
        compress_round14,
        pack_round15,
        batch_round16,
        outer_scan_round17,
        round6,
    } = cli;
    let continuation_modes = usize::from(width_round13.is_some())
        + usize::from(compress_round14.is_some())
        + usize::from(pack_round15.is_some())
        + usize::from(batch_round16.is_some())
        + usize::from(outer_scan_round17.is_some());
    if continuation_modes > 1 {
        return Err("choose only one continuation mode".into());
    }
    if let Some(path) = width_round13 {
        return run_width(&path, out);
    }
    if let Some(path) = compress_round14 {
        return run_compression(&path, out);
    }
    if let Some(path) = pack_round15 {
        return run_packed_compression(&path, out);
    }
    if let Some(path) = batch_round16 {
        return run_batched_inversion(&path, out);
    }
    if let Some(path) = outer_scan_round17 {
        let round6 = round6.ok_or("--outer-scan-round17 requires --round6")?;
        return run_outer_scan(&path, &round6, out);
    }
    if round6.is_some() {
        return Err("--round6 is only valid with --outer-scan-round17".into());
    }
    run_transfer(out)
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s3(u: Fe, v: Fe, t: Fe, a: Fe, b: Fe) -> Fe {
        let q = u.mul(&v);
        let sum = u.add(&v);
        let difference = u.sub(&v);
        let first = t.sqr().mul(&difference.sqr());
        let inner = sum.mul(&q.add(&a)).add(&b.add(&b));
        let second = t.add(&t).mul(&inner).neg();
        let third = q.sub(&a).sqr().sub(&b.add(&b).add(&b.add(&b)).mul(&sum));
        first.add(&second).add(&third)
    }

    #[test]
    fn coefficients_equal_direct_s3() {
        let curve = CurveParams::p256();
        let a = Fe::from_biguint(&curve.a);
        let b = Fe::from_biguint(&curve.b);
        for (u, t) in [(1u64, 2u64), (17, 301), (500, 500), (1_000_003, 42)] {
            let u = Fe::from_biguint(&BigUint::from(u));
            let t = Fe::from_biguint(&BigUint::from(t));
            let mut count = CountedField::default();
            let (qa, qb, qc) = coefficients(u, t, a, b, &mut count);
            for v in [0u64, 1, 2, 99, 1024, 1_000_033] {
                let v = Fe::from_biguint(&BigUint::from(v));
                let polynomial = qa.mul(&v.sqr()).add(&qb.mul(&v)).add(&qc);
                assert_eq!(polynomial, s3(u, v, t, a, b));
            }
        }
    }

    #[test]
    fn quadratic_solver_recovers_constructed_roots_and_linear_case() {
        let curve = CurveParams::p256();
        let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
        let inverse_exponent = &curve.p - BigUint::from(2u8);
        let r1 = Fe::from_biguint(&BigUint::from(123u16));
        let r2 = Fe::from_biguint(&BigUint::from(456u16));
        let qa = Fe::ONE;
        let qb = r1.add(&r2).neg();
        let qc = r1.mul(&r2);
        let mut counts = MultiplicationCounts::default();
        let roots = solve_quadratic(qa, qb, qc, &sqrt_exponent, &inverse_exponent, &mut counts);
        let actual: BTreeSet<[u8; 32]> = roots.roots[..roots.len]
            .iter()
            .copied()
            .map(fe_key)
            .collect();
        let expected = BTreeSet::from([fe_key(r1), fe_key(r2)]);
        assert_eq!(actual, expected);

        let mut counts = MultiplicationCounts::default();
        let linear = solve_quadratic(
            Fe::ZERO,
            Fe::from_biguint(&BigUint::from(7u8)),
            Fe::from_biguint(&BigUint::from(35u8)),
            &sqrt_exponent,
            &inverse_exponent,
            &mut counts,
        );
        assert!(linear.linear);
        assert_eq!(linear.len, 1);
        let seven = Fe::from_biguint(&BigUint::from(7u8));
        let thirty_five = Fe::from_biguint(&BigUint::from(35u8));
        assert_eq!(seven.mul(&linear.roots[0]).add(&thirty_five), Fe::ZERO);
    }

    #[test]
    fn packed_atom_is_canonical_and_round_trips() {
        let columns: Vec<u64> = (0..ATOM_COLUMNS as u64)
            .map(|index| index * 8_001)
            .collect();
        let packet = pack_atom_columns(&columns).unwrap();
        assert_eq!(packet.len(), 36);
        assert_eq!(unpack_atom_columns(&packet).unwrap(), columns);

        let mut aggregate = BigUint::zero();
        for column in &columns {
            aggregate = (aggregate << COLUMN_INDEX_BITS) + BigUint::from(*column);
        }
        let encoded = aggregate.to_bytes_be();
        let mut reference = [0u8; PACKED_ATOM_BYTES];
        reference[PACKED_ATOM_BYTES - encoded.len()..].copy_from_slice(&encoded);
        assert_eq!(packet, reference);
    }

    #[test]
    fn packed_atom_rejects_noncanonical_columns() {
        let mut duplicate: Vec<u64> = (0..ATOM_COLUMNS as u64).collect();
        duplicate[15] = duplicate[0];
        assert!(pack_atom_columns(&duplicate).is_err());

        let mut out_of_range: Vec<u64> = (0..ATOM_COLUMNS as u64).collect();
        out_of_range[15] = COLUMNS;
        assert!(pack_atom_columns(&out_of_range).is_err());
        assert!(unpack_atom_columns(&[0u8; PACKED_ATOM_BYTES - 1]).is_err());
    }

    #[test]
    fn batched_quadratics_match_scalar_roots_and_charge_one_inversion() {
        let curve = CurveParams::p256();
        let sqrt_exponent = (&curve.p + BigUint::one()) >> 2usize;
        let inverse_exponent = &curve.p - BigUint::from(2u8);
        let rows: Vec<(Fe, Fe, Fe)> = [(3u64, 11u64), (17, 29), (101, 303)]
            .into_iter()
            .map(|(left, right)| {
                let left = Fe::from_biguint(&BigUint::from(left));
                let right = Fe::from_biguint(&BigUint::from(right));
                (Fe::ONE, left.add(&right).neg(), left.mul(&right))
            })
            .collect();
        let mut scalar_counts = MultiplicationCounts::default();
        let scalar: Vec<QuadraticRoots> = rows
            .iter()
            .map(|&(qa, qb, qc)| {
                solve_quadratic(
                    qa,
                    qb,
                    qc,
                    &sqrt_exponent,
                    &inverse_exponent,
                    &mut scalar_counts,
                )
            })
            .collect();
        let mut batched_counts = MultiplicationCounts::default();
        let mut profile = BatchInversionProfile::default();
        let batched = solve_quadratic_batch(
            &rows,
            &sqrt_exponent,
            &inverse_exponent,
            &mut batched_counts,
            &mut profile,
        )
        .unwrap();
        for (left, right) in scalar.iter().zip(&batched) {
            let left: BTreeSet<[u8; 32]> =
                left.roots[..left.len].iter().copied().map(fe_key).collect();
            let right: BTreeSet<[u8; 32]> = right.roots[..right.len]
                .iter()
                .copied()
                .map(fe_key)
                .collect();
            assert_eq!(left, right);
        }
        assert_eq!(profile.batch_sizes, [3]);
        assert_eq!(profile.batch_inversions, 1);
        assert_eq!(profile.scalar_fallbacks, 0);
        assert_eq!(scalar_counts.inversions, 3 * 384);
        assert_eq!(batched_counts.inversions, 384 + 3 * 3 - 1);
    }

    #[test]
    fn outer_stream_is_full_cycle_and_freezes_the_planted_column() {
        let mut seen = BTreeSet::new();
        for position in 0..COLUMNS {
            let column = (OUTER_OFFSET + OUTER_STRIDE * position) % COLUMNS;
            assert!(seen.insert(column));
            if position == PLANTED_OUTER_POSITION {
                assert_eq!(column, PLANTED_OUTER_COLUMN);
            }
        }
        assert_eq!(seen.len() as u64, COLUMNS);
        assert_eq!(binomial_big(5, 2), BigUint::from(10u8));
    }
}
