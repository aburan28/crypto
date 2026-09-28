#![recursion_limit = "256"]
#![allow(dead_code)]
//! Kuhn–Struik batched Pollard rho over the Koblitz public-fixture corpus —
//! **matched-arithmetic variant** of `koblitz_rho_batch_ks.rs`.
//!
//! This file is a byte-for-byte copy of `koblitz_rho_batch_ks.rs` except for
//! exactly two functional changes, made to remove an arithmetic asymmetry
//! against `koblitz_orbit_dlp_fast.rs` (see
//! `research/notes/index-calculus/RESEARCH_MATCHED_RHO_ORBIT_DLP_20260928.md`
//! for the frozen protocol this file exists to run):
//!
//! 1. `raw_canonicalize`'s `SignedFrobenius` branch no longer walks the orbit
//!    by calling `raw_square` on the point coordinates `curve.n` times. It
//!    builds a normal basis once (`NormalBasis`, copied from
//!    `koblitz_orbit_dlp_fast.rs`, over the same `Gf2` field the library uses
//!    for this exact curve) and reads every Frobenius image `x^(2^k)` as an
//!    O(1) bit rotation instead.
//! 2. `raw_inverse` no longer does Fermat's-little-theorem square-and-multiply
//!    (`~n` squarings, up to `~n` multiplications). It does Itoh–Tsujii
//!    inversion with the `a^(2^k)` steps done as the same normal-basis
//!    rotation, exactly as `koblitz_orbit_dlp_fast.rs`'s `invert()` does.
//!
//! CLI arguments, environment variables, JSON schema, walk logic, the
//! distinguished-point table, the jump table and negation handling are
//! otherwise identical to `koblitz_rho_batch_ks.rs`, so the two remain a
//! controlled comparison: same corpus, same seeds, same steps, same decisions,
//! different cost only where the two listed changes apply.
//!
//! Operation counting: every real field multiplication (`raw_mul_field`) and
//! squaring (`raw_square`) call increments a global counter, auditable at the
//! call site. Normal-basis rotations and the `Linear` byte-table applies used
//! to build/consult the basis are O(1) word operations and are not counted,
//! matching `koblitz_orbit_dlp_fast.rs`'s own accounting convention (see that
//! file and the protocol note). `NormalBasis::new`'s one-time O(n) library
//! `Gf2::sqr` calls (candidate search at startup) are excluded from the count
//! as well: they are a bounded, one-time setup cost, negligible next to the
//! millions of per-step calls in a full batch, and are called out here
//! explicitly rather than silently folded in or dropped.

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::collections::{BTreeSet, HashMap};
use std::fs;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const JUMPS: usize = 32;
const PRECOMPUTED: u32 = u32::MAX;

/// Real field multiplications (`raw_mul_field` calls), across the whole
/// process. The unit the matched-arithmetic comparison is scored in.
static MUL_CALLS: AtomicU64 = AtomicU64::new(0);
/// Real field squarings (`raw_square` calls), across the whole process.
static SQR_CALLS: AtomicU64 = AtomicU64::new(0);
/// Field-inversion invocations (`raw_inverse` calls). Diagnostic only: an
/// inversion's real cost is already counted through the `raw_mul_field` and
/// `raw_square` calls it makes internally, so this is not added to the op
/// total, only reported alongside it.
static INV_CALLS: AtomicU64 = AtomicU64::new(0);

fn field_ops_snapshot() -> (u64, u64, u64) {
    (
        MUL_CALLS.load(Ordering::Relaxed),
        SQR_CALLS.load(Ordering::Relaxed),
        INV_CALLS.load(Ordering::Relaxed),
    )
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Quotient {
    Ordinary,
    Negation,
    SignedFrobenius,
}

impl Quotient {
    fn parse(value: &str) -> Self {
        match value {
            "ordinary" => Self::Ordinary,
            "negation_only" => Self::Negation,
            "signed_frobenius" => Self::SignedFrobenius,
            _ => panic!("unknown quotient mode {value}"),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Ordinary => "ordinary",
            Self::Negation => "negation_only",
            Self::SignedFrobenius => "signed_frobenius",
        }
    }

    fn uses_negation(self) -> bool {
        self != Self::Ordinary
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RawPoint {
    Infinity,
    Affine { x: u64, y: u64 },
}

#[derive(Clone, Copy)]
struct RawState {
    point: RawPoint,
    a: u64,
    b: u64,
}

#[derive(Clone, Copy)]
struct RawJump {
    point: RawPoint,
    a: u64,
    b: u64,
}

#[derive(Default)]
struct Charges {
    group_additions: u64,
    scalar_multiplications: u64,
    canonicalizations: u64,
    /// Orbit maps applied during canonicalization. Under the matched
    /// arithmetic these are O(1) normal-basis rotations, not `raw_square`
    /// calls, so this remains a diagnostic orbit-size count, not a real-cost
    /// proxy; the real cost is `MUL_CALLS`/`SQR_CALLS` (see `field_ops` in
    /// the JSON output).
    frobenius_maps: u64,
    negations_examined: u64,
    partition_hashes: u64,
    table_queries: u64,
    table_inserts: u64,
    failed_collisions: u64,
    fruitless_cycle_restarts: u64,
}
fn raw_point(point: &BinaryPoint) -> RawPoint {
    match point {
        BinaryPoint::Infinity => RawPoint::Infinity,
        BinaryPoint::Affine { x, y } => RawPoint::Affine {
            x: x.raw_bits().first().copied().unwrap_or(0),
            y: y.raw_bits().first().copied().unwrap_or(0),
        },
    }
}

fn raw_key(point: RawPoint) -> (u8, u64, u64) {
    match point {
        RawPoint::Infinity => (0, 0, 0),
        RawPoint::Affine { x, y } => (1, x, y),
    }
}
fn raw_reduce(curve: &KoblitzCurve, mut wide: u128) -> u64 {
    if curve.n <= 31 && wide <= u64::MAX as u128 {
        let mut narrow = wide as u64;
        let mask = (1u64 << curve.n) - 1;
        while narrow >> curve.n != 0 {
            let high = narrow >> curve.n;
            narrow &= mask;
            for &term in &curve.curve.irreducible.low_terms {
                narrow ^= high << term;
            }
        }
        return narrow;
    }
    let mask = (1u128 << curve.n) - 1;
    while wide >> curve.n != 0 {
        let high = wide >> curve.n;
        wide &= mask;
        for &term in &curve.curve.irreducible.low_terms {
            wide ^= high << term;
        }
    }
    wide as u64
}

fn raw_square(curve: &KoblitzCurve, value: u64) -> u64 {
    SQR_CALLS.fetch_add(1, Ordering::Relaxed);
    if curve.n <= 31 {
        let mut wide = value;
        wide = (wide | (wide << 16)) & 0x0000_ffff_0000_ffff;
        wide = (wide | (wide << 8)) & 0x00ff_00ff_00ff_00ff;
        wide = (wide | (wide << 4)) & 0x0f0f_0f0f_0f0f_0f0f;
        wide = (wide | (wide << 2)) & 0x3333_3333_3333_3333;
        wide = (wide | (wide << 1)) & 0x5555_5555_5555_5555;
        return raw_reduce(curve, wide as u128);
    }
    raw_reduce(curve, carryless_product(value, value))
}

fn carryless_product_software(left: u64, right: u64) -> u128 {
    let mut product = 0u128;
    let mut value = right;
    while value != 0 {
        let bit = value.trailing_zeros();
        product ^= (left as u128) << bit;
        value &= value - 1;
    }
    product
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "aes")]
unsafe fn carryless_product_pmull(left: u64, right: u64) -> u128 {
    std::arch::aarch64::vmull_p64(left, right)
}

fn carryless_product(left: u64, right: u64) -> u128 {
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("aes") {
        // SAFETY: the runtime feature check above proves PMULL availability.
        return unsafe { carryless_product_pmull(left, right) };
    }
    carryless_product_software(left, right)
}

fn raw_mul_field(curve: &KoblitzCurve, left: u64, right: u64) -> u64 {
    MUL_CALLS.fetch_add(1, Ordering::Relaxed);
    raw_reduce(curve, carryless_product(left, right))
}

/// F2-linear map on n <= 64 bits applied by byte tables. Copied verbatim
/// from `examples/koblitz_orbit_dlp_fast.rs`.
struct Linear {
    tables: Vec<[u64; 256]>,
}

impl Linear {
    fn from_images(images: &[u64]) -> Self {
        let chunks = images.len().div_ceil(8);
        let mut tables = vec![[0u64; 256]; chunks];
        for (chunk, table) in tables.iter_mut().enumerate() {
            for byte in 1usize..256 {
                let low = byte & (byte - 1);
                let bit = chunk * 8 + byte.trailing_zeros() as usize;
                let image = images.get(bit).copied().unwrap_or(0);
                table[byte] = table[low] ^ image;
            }
        }
        Self { tables }
    }

    #[inline(always)]
    fn apply(&self, v: u64) -> u64 {
        let mut acc = 0u64;
        for (chunk, table) in self.tables.iter().enumerate() {
            acc ^= table[((v >> (8 * chunk)) & 0xff) as usize];
        }
        acc
    }
}

/// Columns `c_i` form matrix M (poly = M * normal). Returns the images of the
/// poly basis vectors under M^{-1}, or None when M is singular. Copied
/// verbatim from `examples/koblitz_orbit_dlp_fast.rs`.
fn invert_columns(columns: &[u64], n: usize) -> Option<Vec<u64>> {
    let mut rows: Vec<u128> = (0..n)
        .map(|r| {
            let mut row = 0u128;
            for (i, &column) in columns.iter().enumerate() {
                if (column >> r) & 1 == 1 {
                    row |= 1u128 << i;
                }
            }
            row | (1u128 << (64 + r))
        })
        .collect();
    for pivot in 0..n {
        let found = (pivot..n).find(|&r| (rows[r] >> pivot) & 1 == 1)?;
        rows.swap(pivot, found);
        for r in 0..n {
            if r != pivot && (rows[r] >> pivot) & 1 == 1 {
                rows[r] ^= rows[pivot];
            }
        }
    }
    let mut images = vec![0u64; n];
    for (i, row) in rows.iter().enumerate() {
        let high = (row >> 64) as u64;
        for (r, image) in images.iter_mut().enumerate() {
            if (high >> r) & 1 == 1 {
                *image |= 1u64 << i;
            }
        }
    }
    Some(images)
}

/// Frobenius as an O(1) bit rotation in a normal basis, built from `Gf2`
/// conjugates. Same construction and same field (`Gf2::new` on the curve's
/// own `IrreducblePoly`) as `koblitz_orbit_dlp_fast.rs`'s `NormalBasis`;
/// trimmed to the pieces `raw_canonicalize`/`raw_inverse` need
/// (`rotate`/`to_normal`/`to_poly`; the root-table `canonical()` helper is
/// unused here and omitted).
struct NormalBasis {
    n: u32,
    mask: u64,
    to_normal: Linear,
    to_poly: Linear,
}

impl NormalBasis {
    fn new(gf: &Gf2) -> Self {
        let n = gf.n as usize;
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        for candidate in 2u64.. {
            let mut conjugates = Vec::with_capacity(n);
            let mut value = candidate & mask;
            for _ in 0..n {
                conjugates.push(value);
                value = gf.sqr(value);
            }
            let Some(inverse_columns) = invert_columns(&conjugates, n) else {
                continue;
            };
            return Self {
                n: gf.n,
                mask,
                to_normal: Linear::from_images(&inverse_columns),
                to_poly: Linear::from_images(&conjugates),
            };
        }
        unreachable!()
    }

    #[inline(always)]
    fn rotate(&self, v: u64, k: u32) -> u64 {
        if k == 0 {
            return v;
        }
        ((v << k) | (v >> (self.n - k))) & self.mask
    }
}

/// Itoh–Tsujii inversion with the `a^(2^k)` steps done as normal-basis
/// rotations instead of Fermat's-little-theorem square-and-multiply. Same
/// algorithm as `koblitz_orbit_dlp_fast.rs`'s `invert()`, expressed with this
/// file's `raw_mul_field`/`raw_square` (so the same global counters price it).
fn raw_inverse(curve: &KoblitzCurve, nb: &NormalBasis, value: u64) -> u64 {
    assert_ne!(value, 0);
    INV_CALLS.fetch_add(1, Ordering::Relaxed);
    let e = curve.n - 1;
    let mut c = value;
    let mut k = 1u32;
    let bits = 32 - e.leading_zeros();
    for i in (0..bits - 1).rev() {
        let raised = nb.to_poly.apply(nb.rotate(nb.to_normal.apply(c), k));
        c = raw_mul_field(curve, c, raised);
        k *= 2;
        if (e >> i) & 1 == 1 {
            c = raw_mul_field(curve, raw_square(curve, c), value);
            k += 1;
        }
    }
    raw_square(curve, c)
}

fn raw_neg(point: RawPoint) -> RawPoint {
    match point {
        RawPoint::Infinity => RawPoint::Infinity,
        RawPoint::Affine { x, y } => RawPoint::Affine { x, y: y ^ x },
    }
}

fn raw_double(curve: &KoblitzCurve, nb: &NormalBasis, point: RawPoint) -> RawPoint {
    let RawPoint::Affine { x, y } = point else {
        return RawPoint::Infinity;
    };
    if x == 0 {
        return RawPoint::Infinity;
    }
    let lambda = x ^ raw_mul_field(curve, y, raw_inverse(curve, nb, x));
    let x3 = raw_square(curve, lambda) ^ lambda ^ curve.a as u64;
    let y3 = raw_square(curve, x) ^ raw_mul_field(curve, lambda ^ 1, x3);
    RawPoint::Affine { x: x3, y: y3 }
}

fn raw_add(curve: &KoblitzCurve, nb: &NormalBasis, left: RawPoint, right: RawPoint) -> RawPoint {
    match (left, right) {
        (RawPoint::Infinity, point) | (point, RawPoint::Infinity) => point,
        (RawPoint::Affine { x: x1, y: y1 }, RawPoint::Affine { x: x2, y: y2 }) => {
            if x1 == x2 {
                return if y1 ^ y2 == x1 {
                    RawPoint::Infinity
                } else {
                    raw_double(curve, nb, left)
                };
            }
            let lambda = raw_mul_field(curve, y1 ^ y2, raw_inverse(curve, nb, x1 ^ x2));
            let x3 = raw_square(curve, lambda) ^ lambda ^ x1 ^ x2 ^ curve.a as u64;
            let y3 = raw_mul_field(curve, lambda, x1 ^ x3) ^ x3 ^ y1;
            RawPoint::Affine { x: x3, y: y3 }
        }
    }
}

fn raw_scalar_mul(curve: &KoblitzCurve, nb: &NormalBasis, point: RawPoint, scalar: u64) -> RawPoint {
    let mut result = RawPoint::Infinity;
    for bit in (0..64 - scalar.leading_zeros()).rev() {
        result = raw_double(curve, nb, result);
        if (scalar >> bit) & 1 == 1 {
            result = raw_add(curve, nb, result, point);
        }
    }
    result
}

fn raw_on_curve(curve: &KoblitzCurve, point: RawPoint) -> bool {
    let RawPoint::Affine { x, y } = point else {
        return true;
    };
    let x_squared = raw_square(curve, x);
    let left = raw_square(curve, y) ^ raw_mul_field(curve, x, y);
    let right = raw_mul_field(curve, x_squared, x) ^ if curve.a == 1 { x_squared } else { 0 } ^ 1;
    left == right
}

fn public_point_targets(curve: &KoblitzCurve, nb: &NormalBasis, modulus: u64) -> Option<Vec<RawPoint>> {
    let path = std::env::var("KIC_RHO_TARGET_POINTS_JSONL").ok()?;
    let input = fs::read_to_string(path).expect("rho point target file must be readable");
    let limit = 1u64 << curve.n;
    let points = input
        .lines()
        .map(|line| {
            let [x, y]: [u64; 2] =
                serde_json::from_str(line).expect("rho point target must be JSON [x,y]");
            assert!(
                x < limit && y < limit,
                "rho point coordinates must be field elements"
            );
            let point = RawPoint::Affine { x, y };
            assert!(
                raw_on_curve(curve, point),
                "rho point target must be on the curve"
            );
            assert_eq!(
                raw_scalar_mul(curve, nb, point, modulus),
                RawPoint::Infinity,
                "rho point target must belong to the prime-order subgroup"
            );
            point
        })
        .collect::<Vec<_>>();
    Some(points)
}

fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}

fn signed_automorphism_size(lambda: u64, modulus: u64, n: u32) -> usize {
    let mut scalars = BTreeSet::new();
    let mut current = 1u64;
    for _ in 0..n {
        scalars.insert(current);
        scalars.insert((modulus - current) % modulus);
        current = mul_mod(current, lambda, modulus);
    }
    assert_eq!(current, 1);
    scalars.len()
}

fn sub_mod(left: u64, right: u64, modulus: u64) -> u64 {
    if left >= right {
        left - right
    } else {
        modulus - (right - left)
    }
}

fn inverse_mod(value: u64, modulus: u64) -> Option<u64> {
    if value == 0 {
        return None;
    }
    let (mut old_r, mut r) = (modulus as i128, value as i128);
    let (mut old_t, mut t) = (0i128, 1i128);
    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_t, t) = (t, old_t - quotient * t);
    }
    (old_r == 1).then_some(old_t.rem_euclid(modulus as i128) as u64)
}

/// Matched-arithmetic `raw_canonicalize`: the `SignedFrobenius` orbit is read
/// off a normal basis instead of walked by repeated `raw_square`. `point` for
/// orbit position `exponent` is computed directly from the ORIGINAL point's
/// normal coordinates via one O(1) rotation, rather than by chaining through
/// `exponent` real squarings from the previous position — mathematically the
/// same value (`x^(2^exponent)`), reached without the field arithmetic. Tie
/// -breaking, negation handling and the multiplier bookkeeping are otherwise
/// identical to `koblitz_rho_batch_ks.rs`.
fn raw_canonicalize(
    curve: &KoblitzCurve,
    nb: &NormalBasis,
    state: RawState,
    mode: Quotient,
    modulus: u64,
    lambda: u64,
    charges: &mut Charges,
) -> RawState {
    charges.canonicalizations += 1;
    if state.point == RawPoint::Infinity || mode == Quotient::Ordinary {
        return state;
    }
    let powers = if mode == Quotient::SignedFrobenius {
        curve.n
    } else {
        1
    };
    let RawPoint::Affine { x: x0, y: y0 } = state.point else {
        unreachable!("infinity handled above");
    };
    let x_normal = nb.to_normal.apply(x0);
    let y_normal = nb.to_normal.apply(y0);
    let mut multiplier = 1u64;
    let mut best_key = raw_key(state.point);
    let mut best_point = state.point;
    let mut best_multiplier = multiplier;
    for exponent in 0..powers {
        let point = if exponent == 0 {
            state.point
        } else {
            RawPoint::Affine {
                x: nb.to_poly.apply(nb.rotate(x_normal, exponent)),
                y: nb.to_poly.apply(nb.rotate(y_normal, exponent)),
            }
        };
        let key = raw_key(point);
        if key < best_key {
            best_key = key;
            best_point = point;
            best_multiplier = multiplier;
        }
        if mode.uses_negation() {
            charges.negations_examined += 1;
            let negative = raw_neg(point);
            if raw_key(negative) < best_key {
                best_key = raw_key(negative);
                best_point = negative;
                best_multiplier = modulus - multiplier;
            }
        }
        if exponent + 1 < powers {
            multiplier = mul_mod(multiplier, lambda, modulus);
            charges.frobenius_maps += 1;
        }
    }
    RawState {
        point: best_point,
        a: mul_mod(state.a, best_multiplier, modulus),
        b: mul_mod(state.b, best_multiplier, modulus),
    }
}
fn raw_partition(point: RawPoint) -> usize {
    let (_, x, y) = raw_key(point);
    let mut value = x ^ y.rotate_left(21) ^ 0x9e37_79b9_7f4a_7c15;
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^= value >> 31;
    value as usize % JUMPS
}

fn dp_hash(point: RawPoint) -> u64 {
    let (_, x, y) = raw_key(point);
    let mut v = x.rotate_left(7) ^ y ^ 0x2545_f491_4f6c_dd1d;
    v ^= v >> 33;
    v = v.wrapping_mul(0xff51_afd7_ed55_8ccd);
    v ^= v >> 33;
    v = v.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    v ^ (v >> 33)
}

#[derive(Clone, Copy)]
struct Trail {
    a: u64,
    b: u64,
    target: u32,
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        6,
        "usage: <n> <a> <mode> <fixtures> <batch_seed>"
    );
    let n: u32 = args[1].parse().unwrap();
    let a: u8 = args[2].parse().unwrap();
    let mode = Quotient::parse(&args[3]);
    let fixtures: u32 = args[4].parse().unwrap();
    let batch_seed: u64 = args[5].parse().unwrap();
    let dp_bits: u32 = std::env::var("KIC_RHO_DP_BITS")
        .map(|v| v.parse().expect("KIC_RHO_DP_BITS must be an integer"))
        .unwrap_or(8);
    assert!(dp_bits < 32);
    let shared_corpus = std::env::var("KIC_RHO_BATCH_CORPUS").ok();
    assert!(matches!(n, 7 | 11 | 13 | 17 | 19 | 23 | 37 | 41 | 53));
    assert!(fixtures > 0);

    let process_started = Instant::now();
    let curve = KoblitzCurve::new(a, n).expect("frozen exact rung must construct");
    // NormalBasis::new uses the library's own Gf2 (same irreducible
    // polynomial as curve.curve.irreducible, verified against raw_square by
    // this file's own `cargo test --example` property tests) to find a
    // normal-basis generator once, at startup. That one-time O(n) search is
    // not on the counted call path (see the module doc comment).
    let gf = Gf2::new(&curve.curve.irreducible);
    let nb = NormalBasis::new(&gf);
    let modulus = curve.subgroup_order.to_u64_digits()[0];
    let lambda = curve.lambda.to_u64_digits()[0];
    let automorphisms = match mode {
        Quotient::Ordinary => 1,
        Quotient::Negation => 2,
        Quotient::SignedFrobenius => signed_automorphism_size(lambda, modulus, n),
    };
    let generator = raw_point(curve.generator());
    let point_targets = public_point_targets(&curve, &nb, modulus);
    if let Some(points) = &point_targets {
        assert_eq!(
            points.len(),
            fixtures as usize,
            "rho point target count must equal fixtures argument"
        );
    }
    let mut charges = Charges::default();

    let setup_started = Instant::now();
    let field_ops_before_setup = field_ops_snapshot();
    let jump_digest =
        blake3::hash(format!("KIC-KS-BATCH-JUMPS-v1|{n}|{a}|{batch_seed}").as_bytes());
    let mut jump_rng = StdRng::seed_from_u64(u64::from_le_bytes(
        jump_digest.as_bytes()[..8].try_into().unwrap(),
    ));
    let jumps: Vec<RawJump> = (0..JUMPS)
        .map(|_| loop {
            let s = jump_rng.gen_range(1..modulus);
            charges.scalar_multiplications += 1;
            let point = raw_scalar_mul(&curve, &nb, generator, s);
            if point != RawPoint::Infinity {
                break RawJump { point, a: s, b: 0 };
            }
        })
        .collect();
    let setup_ms = setup_started.elapsed().as_secs_f64() * 1000.0;
    let field_ops_after_setup = field_ops_snapshot();

    let dp_mask = (1u64 << dp_bits) - 1;
    let walk_cap = 8u64 << dp_bits;
    let ideal_single =
        (std::f64::consts::PI * modulus as f64 / (2.0 * automorphisms as f64)).sqrt();
    let step_cap = (ideal_single.ceil() as u64)
        .saturating_mul(2_000)
        .max(1_000_000);
    let mut table: HashMap<(u8, u64, u64), Trail> = HashMap::new();
    let mut solved: Vec<u64> = Vec::with_capacity(fixtures as usize);

    // Bernstein–Lange precomputation: G-only walks whose distinguished points
    // have known logarithms. Charged separately from the target loop.
    let precompute_walks: u64 = std::env::var("KIC_RHO_PRECOMPUTE_WALKS")
        .map(|v| {
            v.parse()
                .expect("KIC_RHO_PRECOMPUTE_WALKS must be an integer")
        })
        .unwrap_or(0);
    let precompute_started = Instant::now();
    let mut precompute_steps = 0u64;
    if precompute_walks > 0 {
        let digest = blake3::hash(format!("KIC-KS-PRECOMPUTE-v1|{n}|{a}|{batch_seed}").as_bytes());
        let mut rng = StdRng::seed_from_u64(u64::from_le_bytes(
            digest.as_bytes()[..8].try_into().unwrap(),
        ));
        let stride_a = rng.gen_range(1..modulus);
        let stride = raw_scalar_mul(&curve, &nb, generator, stride_a);
        let mut cursor_a = rng.gen_range(1..modulus);
        let mut cursor = raw_scalar_mul(&curve, &nb, generator, cursor_a);
        'pre: for _ in 0..precompute_walks {
            let (start, start_a) = (cursor, cursor_a);
            cursor = raw_add(&curve, &nb, cursor, stride);
            cursor_a = (cursor_a + stride_a) % modulus;
            if start == RawPoint::Infinity {
                continue;
            }
            let mut state = raw_canonicalize(
                &curve,
                &nb,
                RawState {
                    point: start,
                    a: start_a,
                    b: 0,
                },
                mode,
                modulus,
                lambda,
                &mut charges,
            );
            let mut previous = [RawPoint::Infinity; 4];
            let mut length = 0u64;
            while dp_hash(state.point) & dp_mask != 0 {
                let jump = &jumps[raw_partition(state.point)];
                let next = raw_canonicalize(
                    &curve,
                    &nb,
                    RawState {
                        point: raw_add(&curve, &nb, state.point, jump.point),
                        a: (state.a + jump.a) % modulus,
                        b: 0,
                    },
                    mode,
                    modulus,
                    lambda,
                    &mut charges,
                );
                precompute_steps += 1;
                length += 1;
                if next.point == state.point || previous.contains(&next.point) || length > walk_cap
                {
                    continue 'pre;
                }
                previous = [state.point, previous[0], previous[1], previous[2]];
                state = next;
            }
            table.entry(raw_key(state.point)).or_insert(Trail {
                a: state.a,
                b: 0,
                target: PRECOMPUTED,
            });
        }
    }
    let precompute_ms = precompute_started.elapsed().as_secs_f64() * 1000.0;
    let field_ops_after_precompute = field_ops_snapshot();
    let freeze_table = std::env::var("KIC_RHO_FREEZE_TABLE").is_ok_and(|v| v == "1");
    let precompute_table_entries = table.len();
    let mut total_steps = 0u64;
    let mut cross_solves = 0u32;

    for index in 0..fixtures {
        let material = match &shared_corpus {
            Some(corpus) => {
                format!("KIC-SHARED-PUBLIC-FIXTURE-v1|{n}|{a}|{corpus}|{batch_seed}|{index}")
            }
            None => format!(
                "TASK-KIC-DIRECT-BATCH-20260910|rho|{n}|{a}|{}|{batch_seed}|{index}",
                mode.name()
            ),
        };
        let digest = blake3::hash(material.as_bytes());
        let seed = u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap());
        let mut rng = StdRng::seed_from_u64(seed);
        // Consume the same RNG word in both modes, so the walk schedule remains
        // frozen even though a point-only run has no target discrete-log label.
        let generated_d0 = rng.gen_range(1..modulus);
        let started = Instant::now();
        let field_ops_before_fixture = field_ops_snapshot();
        let q = if let Some(points) = &point_targets {
            points[index as usize]
        } else {
            charges.scalar_multiplications += 1;
            raw_scalar_mul(&curve, &nb, generator, generated_d0)
        };
        let table_before = table.len();
        // Successive walk starts step by a fixed stride: one addition per walk
        // instead of a fresh scalar multiplication.
        let stride_a = rng.gen_range(1..modulus);
        let stride = raw_scalar_mul(&curve, &nb, generator, stride_a);
        let mut cursor_a = rng.gen_range(0..modulus);
        let mut cursor = raw_add(
            &curve,
            &nb,
            raw_scalar_mul(&curve, &nb, generator, cursor_a),
            q,
        );
        charges.scalar_multiplications += 2;
        charges.group_additions += 1;

        let mut steps = 0u64;
        let mut walks = 0u64;
        let mut fruitless = 0u64;
        let mut capped = 0u64;
        let mut wasted_merges = 0u64;
        let mut recovered = None;
        let mut via_target = None;
        'walks: while steps < step_cap {
            walks += 1;
            let start = cursor;
            let start_a = cursor_a;
            cursor = raw_add(&curve, &nb, cursor, stride);
            cursor_a = (cursor_a + stride_a) % modulus;
            charges.group_additions += 1;
            if start == RawPoint::Infinity {
                continue;
            }
            let mut state = raw_canonicalize(
                &curve,
                &nb,
                RawState {
                    point: start,
                    a: start_a,
                    b: 1,
                },
                mode,
                modulus,
                lambda,
                &mut charges,
            );
            let mut previous = [RawPoint::Infinity; 4];
            let mut length = 0u64;
            while dp_hash(state.point) & dp_mask != 0 {
                let jump = &jumps[raw_partition(state.point)];
                charges.partition_hashes += 1;
                let next = RawState {
                    point: raw_add(&curve, &nb, state.point, jump.point),
                    a: (state.a + jump.a) % modulus,
                    b: state.b,
                };
                charges.group_additions += 1;
                let next = raw_canonicalize(&curve, &nb, next, mode, modulus, lambda, &mut charges);
                steps += 1;
                length += 1;
                if next.point == state.point || previous.contains(&next.point) {
                    fruitless += 1;
                    continue 'walks;
                }
                if length > walk_cap {
                    capped += 1;
                    continue 'walks;
                }
                previous = [state.point, previous[0], previous[1], previous[2]];
                state = next;
            }
            let key = raw_key(state.point);
            charges.table_queries += 1;
            let Some(&hit) = table.get(&key) else {
                if !freeze_table {
                    table.insert(
                        key,
                        Trail {
                            a: state.a,
                            b: state.b,
                            target: index,
                        },
                    );
                    charges.table_inserts += 1;
                }
                continue;
            };
            // Both trails reach this point: a + b·d = a' + b'·d' with d' known or = d.
            let candidate = if hit.target == index {
                let denominator = sub_mod(state.b, hit.b, modulus);
                inverse_mod(denominator, modulus)
                    .map(|inv| mul_mod(sub_mod(hit.a, state.a, modulus), inv, modulus))
            } else {
                let base = if hit.target == PRECOMPUTED {
                    0
                } else {
                    solved[hit.target as usize]
                };
                let known =
                    (hit.a as u128 + mul_mod(hit.b, base, modulus) as u128) % modulus as u128;
                inverse_mod(state.b, modulus)
                    .map(|inv| mul_mod(sub_mod(known as u64, state.a, modulus), inv, modulus))
            };
            let Some(candidate) = candidate else {
                wasted_merges += 1;
                continue;
            };
            charges.scalar_multiplications += 1;
            if raw_scalar_mul(&curve, &nb, generator, candidate) == q {
                recovered = Some(candidate);
                via_target = Some(hit.target);
                break;
            }
            charges.failed_collisions += 1;
        }
        let solve_ms = started.elapsed().as_secs_f64() * 1000.0;
        let recovered = recovered.expect("batched rho exceeded the per-target step cap");
        if point_targets.is_none() {
            assert_eq!(recovered, generated_d0);
        }
        let via = via_target.unwrap();
        if via != index {
            cross_solves += 1;
        }
        solved.push(recovered);
        total_steps += steps;
        let q_key = raw_key(q);
        let (fixture_mul_after, fixture_sqr_after, fixture_inv_after) = field_ops_snapshot();
        let (fixture_mul_before, fixture_sqr_before, fixture_inv_before) = field_ops_before_fixture;
        println!(
            "{}",
            json!({
                "kind":"rho_ks_batch_fixture",
                "evidence_class":"measured_rho_observation",
                "n":n,"a":a,"quotient_mode":mode.name(),"automorphism_size":automorphisms,
                "fixture_index":index,"fixture_seed":seed,"batch_seed":batch_seed,
                "target_source":if point_targets.is_some() {"explicit_public_points"} else {"derived_known_scalar"},
                "published_fixture_scalar":point_targets.is_none().then_some(generated_d0),
                "recovered_fixture_scalar":recovered,
                "published_q":[q_key.1,q_key.2],"verified":true,
                "solved_via_target":via,"solved_via_precomputed":via == PRECOMPUTED,"cross_target_solve":via != index,
                "walk_steps":steps,"walks":walks,"fruitless_two_cycles":fruitless,
                "capped_walks":capped,"wasted_merges":wasted_merges,
                "table_entries_before":table_before,"table_entries_after":table.len(),
                "total_ms":solve_ms,
                "ideal_independent_steps":ideal_single,
                "field_ops":{
                    "mul_calls":fixture_mul_after - fixture_mul_before,
                    "sqr_calls":fixture_sqr_after - fixture_sqr_before,
                    "inv_calls":fixture_inv_after - fixture_inv_before,
                },
            })
        );
    }
    let process_ms = process_started.elapsed().as_secs_f64() * 1000.0;
    let (total_mul_calls, total_sqr_calls, total_inv_calls) = field_ops_snapshot();
    let (setup_mul_before, setup_sqr_before, setup_inv_before) = field_ops_before_setup;
    let (setup_mul_after, setup_sqr_after, setup_inv_after) = field_ops_after_setup;
    let (precompute_mul_after, precompute_sqr_after, precompute_inv_after) =
        field_ops_after_precompute;
    println!(
        "{}",
        json!({
            "kind":"rho_ks_batch_summary","producer_version":"v4_point_input_matched_arith",
            "n":n,"a":a,"quotient_mode":mode.name(),"automorphism_size":automorphisms,
            "fixtures":fixtures,"batch_seed":batch_seed,"corpus":shared_corpus,"dp_bits":dp_bits,"jump_count":JUMPS,
            "target_source":if point_targets.is_some() {"explicit_public_points"} else {"derived_known_scalar"},
            "all_verified":true,"cross_target_solves":cross_solves,
            "total_walk_steps":total_steps,"table_entries":table.len(),
            "table_payload_lower_bound_bytes":table.len() * (1 + 4 * std::mem::size_of::<u64>()),
            "setup_ms":setup_ms,"in_process_ms":process_ms,
            "precompute_walks":precompute_walks,"freeze_table":freeze_table,"precompute_steps":precompute_steps,"precompute_ms":precompute_ms,"precompute_table_entries":precompute_table_entries,
            "charges":{
                "group_additions":charges.group_additions,
                "scalar_multiplications":charges.scalar_multiplications,
                "canonicalizations":charges.canonicalizations,
                "frobenius_maps":charges.frobenius_maps,
                "negations_examined":charges.negations_examined,
                "partition_hashes":charges.partition_hashes,
                "table_queries":charges.table_queries,
                "table_inserts":charges.table_inserts,
                "failed_collisions":charges.failed_collisions
            },
            "field_ops":{
                "unit":"gf2_mul_or_sqr_call",
                "note":"counts every raw_mul_field/raw_square call, i.e. the real carryless-multiply-and-reduce operations the matched-arithmetic run performs; normal-basis rotations and Linear::apply byte-table lookups are O(1) word ops and are not counted; NormalBasis::new's one-time O(n) library Gf2::sqr calls at startup are excluded as negligible (see module doc comment)",
                "setup_mul_calls":setup_mul_after - setup_mul_before,
                "setup_sqr_calls":setup_sqr_after - setup_sqr_before,
                "setup_inv_calls":setup_inv_after - setup_inv_before,
                "precompute_mul_calls":precompute_mul_after - setup_mul_after,
                "precompute_sqr_calls":precompute_sqr_after - setup_sqr_after,
                "precompute_inv_calls":precompute_inv_after - setup_inv_after,
                "targets_mul_calls":total_mul_calls - precompute_mul_after,
                "targets_sqr_calls":total_sqr_calls - precompute_sqr_after,
                "targets_inv_calls":total_inv_calls - precompute_inv_after,
                "total_mul_calls":total_mul_calls,
                "total_sqr_calls":total_sqr_calls,
                "total_inv_calls":total_inv_calls,
                "total_field_ops":total_mul_calls + total_sqr_calls,
            },
            "scope":if point_targets.is_some() {
                "public synthetic point-only targets; no scalar labels supplied to producer"
            } else {
                "published synthetic toy fixtures; no external point, unknown scalar, or production key"
            }
        })
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng_from(seed: u64) -> StdRng {
        StdRng::seed_from_u64(seed)
    }

    /// Correctness gate 1: `Gf2::mul`/`Gf2::sqr` (the library field this
    /// file's `NormalBasis` is built over) agree with this file's own
    /// `raw_mul_field`/`raw_square` on every input, for every admitted n.
    /// If this ever fails, the normal-basis substitution below is unsound:
    /// `Gf2` would not be the same field `KoblitzCurve`'s raw arithmetic uses.
    #[test]
    fn gf2_matches_raw_field_arithmetic() {
        for &n in &[7u32, 11, 13, 17, 19, 23, 37, 41, 53] {
            for a in [0u8, 1u8] {
                let Some(curve) = KoblitzCurve::new(a, n) else {
                    continue;
                };
                let gf = Gf2::new(&curve.curve.irreducible);
                let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
                let mut rng = rng_from(0x6b6f_626c_6974_7a00 ^ (n as u64) ^ ((a as u64) << 40));
                for _ in 0..2000 {
                    let x: u64 = rng.gen::<u64>() & mask;
                    let y: u64 = rng.gen::<u64>() & mask;
                    assert_eq!(
                        raw_square(&curve, x),
                        gf.sqr(x),
                        "raw_square disagrees with Gf2::sqr at n={n} a={a} x={x}"
                    );
                    assert_eq!(
                        raw_mul_field(&curve, x, y),
                        gf.mul(x, y),
                        "raw_mul_field disagrees with Gf2::mul at n={n} a={a} x={x} y={y}"
                    );
                }
            }
        }
    }

    /// Correctness gate 2: the normal-basis rotation reproduces `k`
    /// applications of `raw_square` for every k in 0..n, i.e. the exact
    /// value `raw_canonicalize`'s rewritten loop now reads off in O(1).
    #[test]
    fn normal_basis_rotation_matches_repeated_squaring() {
        for &n in &[7u32, 11, 13, 17, 19, 23, 37, 41, 53] {
            let Some(curve) = KoblitzCurve::new(0, n) else {
                continue;
            };
            let gf = Gf2::new(&curve.curve.irreducible);
            let nb = NormalBasis::new(&gf);
            let mask = (1u64 << n) - 1;
            let mut rng = rng_from(0x726f_7461_7465_0000 ^ (n as u64));
            for _ in 0..500 {
                let x: u64 = rng.gen::<u64>() & mask;
                let mut expected = x;
                let x_normal = nb.to_normal.apply(x);
                for k in 0..n {
                    assert_eq!(
                        nb.to_poly.apply(nb.rotate(x_normal, k)),
                        expected,
                        "rotation by {k} disagrees with {k} chained raw_square calls at n={n} x={x}"
                    );
                    expected = raw_square(&curve, expected);
                }
            }
        }
    }

    /// Correctness gate 3: the matched-arithmetic `raw_canonicalize`
    /// (`SignedFrobenius`) picks the exact same canonical representative
    /// point (and the same multiplier) as the original repeated-squaring
    /// implementation, for a batch of random points. The algorithm's output
    /// must not change; only its cost may.
    #[test]
    fn canonicalization_output_unchanged_by_the_matched_arithmetic() {
        fn original_canonicalize(
            curve: &KoblitzCurve,
            state: RawState,
            modulus: u64,
            lambda: u64,
        ) -> RawState {
            if state.point == RawPoint::Infinity {
                return state;
            }
            let powers = curve.n;
            let mut point = state.point;
            let mut multiplier = 1u64;
            let mut best_key = raw_key(point);
            let mut best_point = point;
            let mut best_multiplier = multiplier;
            for exponent in 0..powers {
                let key = raw_key(point);
                if key < best_key {
                    best_key = key;
                    best_point = point;
                    best_multiplier = multiplier;
                }
                let negative = raw_neg(point);
                if raw_key(negative) < best_key {
                    best_key = raw_key(negative);
                    best_point = negative;
                    best_multiplier = modulus - multiplier;
                }
                if exponent + 1 < powers {
                    point = match point {
                        RawPoint::Infinity => RawPoint::Infinity,
                        RawPoint::Affine { x, y } => RawPoint::Affine {
                            x: raw_square(curve, x),
                            y: raw_square(curve, y),
                        },
                    };
                    multiplier = mul_mod(multiplier, lambda, modulus);
                }
            }
            RawState {
                point: best_point,
                a: mul_mod(state.a, best_multiplier, modulus),
                b: mul_mod(state.b, best_multiplier, modulus),
            }
        }

        for &n in &[7u32, 11, 13, 17, 19, 23, 37, 41, 53] {
            let Some(curve) = KoblitzCurve::new(0, n) else {
                continue;
            };
            let gf = Gf2::new(&curve.curve.irreducible);
            let nb = NormalBasis::new(&gf);
            let modulus = curve.subgroup_order.to_u64_digits()[0];
            let lambda = curve.lambda.to_u64_digits()[0];
            let mask = (1u64 << n) - 1;
            let generator = raw_point(curve.generator());
            let mut rng = rng_from(0x6361_6e6f_6e5f_636b ^ (n as u64));
            let mut examined = 0;
            while examined < 200 {
                let scalar = rng.gen_range(1..modulus);
                let point = raw_scalar_mul(&curve, &nb, generator, scalar);
                if point == RawPoint::Infinity {
                    continue;
                }
                let RawPoint::Affine { x, y } = point else {
                    continue;
                };
                let _ = mask;
                let a = rng.gen_range(0..modulus);
                let b = rng.gen_range(0..modulus);
                let start = RawState {
                    point: RawPoint::Affine { x, y },
                    a,
                    b,
                };
                let expected = original_canonicalize(&curve, start, modulus, lambda);
                let mut charges = Charges::default();
                let got = raw_canonicalize(
                    &curve,
                    &nb,
                    start,
                    Quotient::SignedFrobenius,
                    modulus,
                    lambda,
                    &mut charges,
                );
                assert_eq!(
                    got.point, expected.point,
                    "canonical representative changed at n={n} scalar={scalar}"
                );
                assert_eq!(
                    got.a, expected.a,
                    "canonicalized coefficient a changed at n={n} scalar={scalar}"
                );
                assert_eq!(
                    got.b, expected.b,
                    "canonicalized coefficient b changed at n={n} scalar={scalar}"
                );
                examined += 1;
            }
        }
    }

    /// Correctness gate 4: the Itoh–Tsujii `raw_inverse` is a genuine
    /// multiplicative inverse under `raw_mul_field`, for many random nonzero
    /// field elements, at every admitted n.
    #[test]
    fn raw_inverse_is_a_correct_multiplicative_inverse() {
        for &n in &[7u32, 11, 13, 17, 19, 23, 37, 41, 53] {
            let Some(curve) = KoblitzCurve::new(0, n) else {
                continue;
            };
            let gf = Gf2::new(&curve.curve.irreducible);
            let nb = NormalBasis::new(&gf);
            let mask = (1u64 << n) - 1;
            let mut rng = rng_from(0x696e_7665_7274_0000 ^ (n as u64));
            let mut examined = 0;
            while examined < 500 {
                let x: u64 = rng.gen::<u64>() & mask;
                if x == 0 {
                    continue;
                }
                let inv = raw_inverse(&curve, &nb, x);
                assert_eq!(
                    raw_mul_field(&curve, x, inv),
                    1,
                    "raw_inverse({x}) is not a multiplicative inverse at n={n}"
                );
                examined += 1;
            }
        }
    }
}
