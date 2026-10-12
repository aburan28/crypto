#![recursion_limit = "256"]
#![allow(dead_code)]
//! Kuhn–Struik batched Pollard rho over the Koblitz public-fixture corpus.
//!
//! Arithmetic, canonicalization and partition are copied verbatim from
//! `koblitz_rho_fixture.rs` (packed backend) so per-step cost matches the
//! independent comparator. Jumps are multiples of G only and one
//! distinguished-point table persists across targets, so later targets can
//! finish on the trails of earlier, already-solved ones. Targets are derived
//! exactly as the frozen `packed <batch_seed>` independent runs derive them.
//! `KIC_RHO_POINT_INPUT` accepts public [x,y] JSONL without reading the
//! known-answer scalars. `KIC_RHO_GENERATE_ONLY=1` emits a separate fixture
//! stream for a frozen public-point corpus without running the rho search.

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::collections::{BTreeSet, HashMap};
use std::ops::Deref;
use std::time::Instant;

const JUMPS: usize = 32;
const PARALLEL_WALKS: usize = 32;

/// F2-linear map on n <= 64 bits applied by byte tables.
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
            let basis = Self {
                n: gf.n,
                mask,
                to_normal: Linear::from_images(&inverse_columns),
                to_poly: Linear::from_images(&conjugates),
            };
            return basis;
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

    /// Smallest rotation of `v` and the shift that produces it.
    #[inline(always)]
    fn canonical(&self, v: u64) -> (u64, u32) {
        let mut best = v;
        let mut best_shift = 0;
        let mut current = v;
        for shift in 1..self.n {
            current = ((current << 1) | (current >> (self.n - 1))) & self.mask;
            if current < best {
                best = current;
                best_shift = shift;
            }
        }
        (best, best_shift)
    }
}

/// Columns `c_i` form matrix M (poly = M * normal). Returns the images of the
/// poly basis vectors under M^{-1}, or None when M is singular.
fn invert_columns(columns: &[u64], n: usize) -> Option<Vec<u64>> {
    // Row r of [M | I]: low n bits are M[r][i] over normal index i.
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
    // Row i now reads x_i = sum_r Minv[i][r] y_r in its high half.
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

/// The rho walk uses the same table-reduced single-word field as the
/// compact-orbit producer. Construction is charged inside the process.
struct RhoCurve {
    base: KoblitzCurve,
    field: Gf2,
    normal_basis: Option<NormalBasis>,
    lambda_pows: Vec<u64>,
}

impl RhoCurve {
    fn new(a: u8, n: u32) -> Self {
        Self::with_normal_basis(a, n, false)
    }

    fn with_normal_basis(a: u8, n: u32, use_normal_basis: bool) -> Self {
        let base = KoblitzCurve::new(a, n).expect("frozen exact rung must construct");
        let field = Gf2::new(&base.curve.irreducible);
        let normal_basis = use_normal_basis.then(|| NormalBasis::new(&field));
        let modulus = base.subgroup_order.to_u64_digits()[0];
        let lambda = base.lambda.to_u64_digits()[0];
        let mut lambda_pows = Vec::with_capacity(n as usize);
        let mut value = 1u64;
        for _ in 0..n {
            lambda_pows.push(value);
            value = mul_mod(value, lambda, modulus);
        }
        Self {
            base,
            field,
            normal_basis,
            lambda_pows,
        }
    }
}

impl Deref for RhoCurve {
    type Target = KoblitzCurve;

    fn deref(&self) -> &Self::Target {
        &self.base
    }
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
    frobenius_maps: u64,
    negations_examined: u64,
    partition_hashes: u64,
    table_queries: u64,
    table_inserts: u64,
    failed_collisions: u64,
    fruitless_cycle_restarts: u64,
    batch_inversion_calls: u64,
    batch_inversion_inputs: u64,
    batch_fallback_additions: u64,
    normal_basis_transforms: u64,
    canonical_rotation_steps: u64,
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
fn raw_reduce(curve: &RhoCurve, mut wide: u128) -> u64 {
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

fn raw_square(curve: &RhoCurve, value: u64) -> u64 {
    curve.field.sqr(value)
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

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "pclmulqdq")]
unsafe fn carryless_product_pclmul(left: u64, right: u64) -> u128 {
    use std::arch::x86_64::*;
    let a = _mm_set_epi64x(0, left as i64);
    let b = _mm_set_epi64x(0, right as i64);
    let product = _mm_clmulepi64_si128::<0x00>(a, b);
    let low = _mm_cvtsi128_si64(product) as u64;
    let high = _mm_cvtsi128_si64(_mm_srli_si128::<8>(product)) as u64;
    ((high as u128) << 64) | low as u128
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "aes")]
unsafe fn carryless_product_pmull(left: u64, right: u64) -> u128 {
    std::arch::aarch64::vmull_p64(left, right)
}

fn field_product_backend() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("pclmulqdq") {
        return "pclmulqdq";
    }
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("aes") {
        return "pmull";
    }
    "portable"
}

fn carryless_product(left: u64, right: u64) -> u128 {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("pclmulqdq") {
        // SAFETY: guarded by the runtime CPU-feature check above.
        return unsafe { carryless_product_pclmul(left, right) };
    }
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("aes") {
        // SAFETY: the runtime feature check above proves PMULL availability.
        return unsafe { carryless_product_pmull(left, right) };
    }
    carryless_product_software(left, right)
}

fn raw_mul_field(curve: &RhoCurve, left: u64, right: u64) -> u64 {
    curve.field.mul(left, right)
}

fn raw_inverse(curve: &RhoCurve, value: u64) -> u64 {
    assert_ne!(value, 0);
    curve.field.inv(value)
}

fn raw_neg(point: RawPoint) -> RawPoint {
    match point {
        RawPoint::Infinity => RawPoint::Infinity,
        RawPoint::Affine { x, y } => RawPoint::Affine { x, y: y ^ x },
    }
}

fn raw_double(curve: &RhoCurve, point: RawPoint) -> RawPoint {
    let RawPoint::Affine { x, y } = point else {
        return RawPoint::Infinity;
    };
    if x == 0 {
        return RawPoint::Infinity;
    }
    let lambda = x ^ raw_mul_field(curve, y, raw_inverse(curve, x));
    let x3 = raw_square(curve, lambda) ^ lambda ^ curve.a as u64;
    let y3 = raw_square(curve, x) ^ raw_mul_field(curve, lambda ^ 1, x3);
    RawPoint::Affine { x: x3, y: y3 }
}

fn raw_add(curve: &RhoCurve, left: RawPoint, right: RawPoint) -> RawPoint {
    match (left, right) {
        (RawPoint::Infinity, point) | (point, RawPoint::Infinity) => point,
        (RawPoint::Affine { x: x1, y: y1 }, RawPoint::Affine { x: x2, y: y2 }) => {
            if x1 == x2 {
                return if y1 ^ y2 == x1 {
                    RawPoint::Infinity
                } else {
                    raw_double(curve, left)
                };
            }
            let lambda = raw_mul_field(curve, y1 ^ y2, raw_inverse(curve, x1 ^ x2));
            let x3 = raw_square(curve, lambda) ^ lambda ^ x1 ^ x2 ^ curve.a as u64;
            let y3 = raw_mul_field(curve, lambda, x1 ^ x3) ^ x3 ^ y1;
            RawPoint::Affine { x: x3, y: y3 }
        }
    }
}

/// Add one step for up to 32 independent walks with one field inversion.
/// The nondegenerate lanes use Montgomery's trick; exceptional pairs use
/// the ordinary complete group law, so a batch is bit-for-bit equivalent
/// to applying `raw_add` independently to each active lane.
fn raw_add_pairwise(
    curve: &RhoCurve,
    left: &[RawPoint; PARALLEL_WALKS],
    right: &[RawPoint; PARALLEL_WALKS],
    active: &[bool; PARALLEL_WALKS],
    charges: &mut Charges,
) -> [RawPoint; PARALLEL_WALKS] {
    let mut denoms = [0u64; PARALLEL_WALKS];
    let mut prefixes = [0u64; PARALLEL_WALKS];
    let mut inverses = [0u64; PARALLEL_WALKS];
    let mut ordinary = [false; PARALLEL_WALKS];
    let mut product = 1u64;
    let mut count = 0u64;
    for i in 0..PARALLEL_WALKS {
        if !active[i] {
            continue;
        }
        if let (RawPoint::Affine { x: x1, .. }, RawPoint::Affine { x: x2, .. }) =
            (left[i], right[i])
        {
            if x1 != x2 {
                ordinary[i] = true;
                denoms[i] = x1 ^ x2;
                prefixes[i] = product;
                product = raw_mul_field(curve, product, denoms[i]);
                count += 1;
            }
        }
    }
    if count != 0 {
        charges.batch_inversion_calls += 1;
        charges.batch_inversion_inputs += count;
        let mut inverse_product = raw_inverse(curve, product);
        for i in (0..PARALLEL_WALKS).rev() {
            if ordinary[i] {
                inverses[i] = raw_mul_field(curve, inverse_product, prefixes[i]);
                inverse_product = raw_mul_field(curve, inverse_product, denoms[i]);
            }
        }
        debug_assert_eq!(inverse_product, 1);
    }
    let mut out = [RawPoint::Infinity; PARALLEL_WALKS];
    for i in 0..PARALLEL_WALKS {
        if !active[i] {
            continue;
        }
        if ordinary[i] {
            let (RawPoint::Affine { x: x1, y: y1 }, RawPoint::Affine { x: x2, y: y2 }) =
                (left[i], right[i])
            else {
                unreachable!();
            };
            let lambda = raw_mul_field(curve, y1 ^ y2, inverses[i]);
            let x3 = raw_square(curve, lambda) ^ lambda ^ x1 ^ x2 ^ curve.a as u64;
            let y3 = raw_mul_field(curve, lambda, x1 ^ x3) ^ x3 ^ y1;
            out[i] = RawPoint::Affine { x: x3, y: y3 };
        } else {
            charges.batch_fallback_additions += 1;
            out[i] = raw_add(curve, left[i], right[i]);
        }
    }
    out
}

fn raw_scalar_mul(curve: &RhoCurve, point: RawPoint, scalar: u64) -> RawPoint {
    let mut result = RawPoint::Infinity;
    for bit in (0..64 - scalar.leading_zeros()).rev() {
        result = raw_double(curve, result);
        if (scalar >> bit) & 1 == 1 {
            result = raw_add(curve, result, point);
        }
    }
    result
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

fn raw_canonicalize_reference(
    curve: &RhoCurve,
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
            point = match point {
                RawPoint::Infinity => RawPoint::Infinity,
                RawPoint::Affine { x, y } => RawPoint::Affine {
                    x: raw_square(curve, x),
                    y: raw_square(curve, y),
                },
            };
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
/// Canonicalize in a normal basis: Frobenius is a word rotation.  The
/// complete `(x, min(y, x+y))` comparison covers short x-orbits too.
fn normal_canonicalize(
    curve: &RhoCurve,
    state: RawState,
    mode: Quotient,
    modulus: u64,
    charges: &mut Charges,
) -> RawState {
    let basis = curve.normal_basis.as_ref().unwrap();
    let RawPoint::Affine { x, y } = state.point else {
        unreachable!();
    };
    let powers = if mode == Quotient::SignedFrobenius {
        curve.n
    } else {
        1
    };
    let mut current_x = basis.to_normal.apply(x);
    let mut current_y = basis.to_normal.apply(y);
    let mut best = (u64::MAX, u64::MAX);
    let mut best_shift = 0u32;
    let mut best_negated = false;
    for shift in 0..powers {
        let negated_y = current_x ^ current_y;
        let (candidate_y, negated) = if mode.uses_negation() && negated_y < current_y {
            (negated_y, true)
        } else {
            (current_y, false)
        };
        let key = (current_x, candidate_y);
        if key < best {
            best = key;
            best_shift = shift;
            best_negated = negated;
        }
        if shift + 1 < powers {
            current_x = basis.rotate(current_x, 1);
            current_y = basis.rotate(current_y, 1);
        }
    }
    charges.canonicalizations += 1;
    charges.frobenius_maps += u64::from(powers - 1);
    charges.canonical_rotation_steps += 2 * u64::from(powers - 1);
    charges.normal_basis_transforms += 4;
    if mode.uses_negation() {
        charges.negations_examined += u64::from(powers);
    }
    let mut factor = curve.lambda_pows[best_shift as usize];
    if best_negated {
        factor = modulus - factor;
    }
    RawState {
        point: RawPoint::Affine {
            x: basis.to_poly.apply(best.0),
            y: basis.to_poly.apply(best.1),
        },
        a: mul_mod(state.a, factor, modulus),
        b: mul_mod(state.b, factor, modulus),
    }
}

/// Find the reference signed-Frobenius representative from the x orbit.
/// Ordinarily only the winning y power must be squared; an exceptional
/// short x orbit falls back to the full reference enumeration.
fn raw_canonicalize(
    curve: &RhoCurve,
    state: RawState,
    mode: Quotient,
    modulus: u64,
    lambda: u64,
    charges: &mut Charges,
) -> RawState {
    if state.point == RawPoint::Infinity || mode == Quotient::Ordinary {
        charges.canonicalizations += 1;
        return state;
    }
    if curve.normal_basis.is_some() {
        return normal_canonicalize(curve, state, mode, modulus, charges);
    }
    let powers = if mode == Quotient::SignedFrobenius {
        curve.n
    } else {
        1
    };
    let RawPoint::Affine { x: x0, y: y0 } = state.point else {
        unreachable!();
    };
    let mut x = x0;
    let mut best_x = x0;
    let mut best_k = 0u32;
    let mut tied = false;
    for exponent in 1..powers {
        x = raw_square(curve, x);
        if x < best_x {
            best_x = x;
            best_k = exponent;
            tied = false;
        } else if x == best_x {
            tied = true;
        }
    }
    if tied {
        return raw_canonicalize_reference(curve, state, mode, modulus, lambda, charges);
    }
    charges.canonicalizations += 1;
    charges.frobenius_maps += u64::from(powers - 1);
    if mode.uses_negation() {
        charges.negations_examined += u64::from(powers);
    }
    let mut y = y0;
    let mut factor = 1u64;
    for _ in 0..best_k {
        y = raw_square(curve, y);
        factor = mul_mod(factor, lambda, modulus);
    }
    if mode.uses_negation() && (best_x ^ y) < y {
        y ^= best_x;
        factor = modulus - factor;
    }
    RawState {
        point: RawPoint::Affine { x: best_x, y },
        a: mul_mod(state.a, factor, modulus),
        b: mul_mod(state.b, factor, modulus),
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

#[derive(Clone, Copy)]
struct WalkLane {
    state: RawState,
    previous: [RawPoint; 4],
    length: u64,
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
    let canonicalization_backend =
        std::env::var("KIC_RHO_CANON_BACKEND").unwrap_or_else(|_| "normal_basis".to_string());
    assert!(matches!(
        canonicalization_backend.as_str(),
        "normal_basis" | "poly_xfirst"
    ));
    let fixtures: u32 = args[4].parse().unwrap();
    let batch_seed: u64 = args[5].parse().unwrap();
    let dp_bits: u32 = std::env::var("KIC_RHO_DP_BITS")
        .map(|v| v.parse().expect("KIC_RHO_DP_BITS must be an integer"))
        .unwrap_or(8);
    assert!(dp_bits < 32);
    let shared_corpus = std::env::var("KIC_RHO_BATCH_CORPUS").ok();
    assert!(matches!(n, 7 | 11 | 13 | 17 | 19 | 23 | 37 | 41 | 53 | 61));
    assert!(fixtures > 0);

    let process_started = Instant::now();
    let curve = RhoCurve::with_normal_basis(a, n, canonicalization_backend == "normal_basis");
    let modulus = curve.subgroup_order.to_u64_digits()[0];
    let lambda = curve.lambda.to_u64_digits()[0];
    let automorphisms = match mode {
        Quotient::Ordinary => 1,
        Quotient::Negation => 2,
        Quotient::SignedFrobenius => signed_automorphism_size(lambda, modulus, n),
    };
    let generator = raw_point(curve.generator());
    let point_input = std::env::var("KIC_RHO_POINT_INPUT").ok().map(|path| {
        let input = std::fs::read_to_string(path).expect("read public point JSONL");
        let points: Vec<RawPoint> = input
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                let [x, y]: [u64; 2] =
                    serde_json::from_str(line).expect("public target must be [x,y]");
                assert!(x < (1u64 << n) && y < (1u64 << n));
                RawPoint::Affine { x, y }
            })
            .collect();
        assert_eq!(points.len(), fixtures as usize);
        points
    });
    let generate_only = std::env::var("KIC_RHO_GENERATE_ONLY")
        .map(|value| value == "1")
        .unwrap_or(false);
    if generate_only {
        assert!(
            point_input.is_none(),
            "generation and point input are exclusive"
        );
        assert!(
            shared_corpus.is_some(),
            "generation requires a named corpus"
        );
        for index in 0..fixtures {
            let material = format!(
                "KIC-SHARED-PUBLIC-FIXTURE-v1|{n}|{a}|{}|{batch_seed}|{index}",
                shared_corpus.as_ref().unwrap()
            );
            let digest = blake3::hash(material.as_bytes());
            let seed = u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap());
            let mut rng = StdRng::seed_from_u64(seed);
            let scalar = rng.gen_range(1..modulus);
            let q = raw_key(raw_scalar_mul(&curve, generator, scalar));
            println!(
                "{}",
                json!({
                    "kind":"rho_ks_public_fixture", "n":n, "a":a,
                    "fixture_index":index, "fixture_seed":seed,
                    "batch_seed":batch_seed, "corpus":shared_corpus,
                    "subgroup_order":modulus, "automorphism_size":automorphisms,
                    "field_modulus_low_terms":curve.curve.irreducible.low_terms,
                    "generator":[raw_key(generator).1,raw_key(generator).2],
                    "published_fixture_scalar":scalar, "published_q":[q.1,q.2],
                })
            );
        }
        return;
    }
    if point_input.is_some() {
        assert!(
            shared_corpus.is_some(),
            "public point input requires a named corpus"
        );
    }
    let mut charges = Charges::default();

    let setup_started = Instant::now();
    let jump_digest =
        blake3::hash(format!("KIC-KS-BATCH-JUMPS-v1|{n}|{a}|{batch_seed}").as_bytes());
    let mut jump_rng = StdRng::seed_from_u64(u64::from_le_bytes(
        jump_digest.as_bytes()[..8].try_into().unwrap(),
    ));
    let jumps: Vec<RawJump> = (0..JUMPS)
        .map(|_| loop {
            let s = jump_rng.gen_range(1..modulus);
            charges.scalar_multiplications += 1;
            let point = raw_scalar_mul(&curve, generator, s);
            if point != RawPoint::Infinity {
                break RawJump { point, a: s, b: 0 };
            }
        })
        .collect();
    let setup_ms = setup_started.elapsed().as_secs_f64() * 1000.0;

    let dp_mask = (1u64 << dp_bits) - 1;
    let walk_cap = 8u64 << dp_bits;
    let ideal_single =
        (std::f64::consts::PI * modulus as f64 / (2.0 * automorphisms as f64)).sqrt();
    let step_cap = (ideal_single.ceil() as u64)
        .saturating_mul(2_000)
        .max(1_000_000);
    let mut table: HashMap<(u8, u64, u64), Trail> = HashMap::new();
    let mut solved: Vec<u64> = Vec::with_capacity(fixtures as usize);
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
        let d0 = if point_input.is_some() {
            None
        } else {
            Some(rng.gen_range(1..modulus))
        };
        let started = Instant::now();
        let q = if let Some(points) = &point_input {
            points[index as usize]
        } else {
            charges.scalar_multiplications += 1;
            raw_scalar_mul(&curve, generator, d0.unwrap())
        };
        let table_before = table.len();
        // Successive walk starts step by a fixed stride: one addition per walk
        // instead of a fresh scalar multiplication.
        let stride_a = rng.gen_range(1..modulus);
        let stride = raw_scalar_mul(&curve, generator, stride_a);
        let mut cursor_a = rng.gen_range(0..modulus);
        let mut cursor = raw_add(&curve, raw_scalar_mul(&curve, generator, cursor_a), q);
        charges.scalar_multiplications += 2;
        charges.group_additions += 1;

        let mut steps = 0u64;
        let mut walks = 0u64;
        let mut fruitless = 0u64;
        let mut capped = 0u64;
        let mut wasted_merges = 0u64;
        let mut recovered = None;
        let mut via_target = None;
        let mut packed_batches = 0u64;
        let mut lanes: [Option<WalkLane>; PARALLEL_WALKS] = [None; PARALLEL_WALKS];
        'search: while steps < step_cap {
            // Fill vacant lanes with deterministic stride starts.  Each lane
            // is an independent walk; all lanes share the target's DP table.
            for lane in &mut lanes {
                if lane.is_some() {
                    continue;
                }
                walks += 1;
                let start = cursor;
                let start_a = cursor_a;
                cursor = raw_add(&curve, cursor, stride);
                cursor_a = (cursor_a + stride_a) % modulus;
                charges.group_additions += 1;
                if start == RawPoint::Infinity {
                    continue;
                }
                *lane = Some(WalkLane {
                    state: raw_canonicalize(
                        &curve,
                        RawState {
                            point: start,
                            a: start_a,
                            b: 1,
                        },
                        mode,
                        modulus,
                        lambda,
                        &mut charges,
                    ),
                    previous: [RawPoint::Infinity; 4],
                    length: 0,
                });
            }
            // Resolve distinguished points in lane order before stepping
            // again.  A collision against this or an earlier target can
            // yield the log; failed collisions remain charged and retained.
            for lane in &mut lanes {
                let Some(current) = *lane else {
                    continue;
                };
                if dp_hash(current.state.point) & dp_mask != 0 {
                    continue;
                }
                *lane = None;
                let state = current.state;
                let key = raw_key(state.point);
                charges.table_queries += 1;
                let Some(&hit) = table.get(&key) else {
                    table.insert(
                        key,
                        Trail {
                            a: state.a,
                            b: state.b,
                            target: index,
                        },
                    );
                    charges.table_inserts += 1;
                    continue;
                };
                let candidate = if hit.target == index {
                    let denominator = sub_mod(state.b, hit.b, modulus);
                    inverse_mod(denominator, modulus)
                        .map(|inv| mul_mod(sub_mod(hit.a, state.a, modulus), inv, modulus))
                } else {
                    let known = (hit.a as u128
                        + mul_mod(hit.b, solved[hit.target as usize], modulus) as u128)
                        % modulus as u128;
                    inverse_mod(state.b, modulus)
                        .map(|inv| mul_mod(sub_mod(known as u64, state.a, modulus), inv, modulus))
                };
                let Some(candidate) = candidate else {
                    wasted_merges += 1;
                    continue;
                };
                charges.scalar_multiplications += 1;
                if raw_scalar_mul(&curve, generator, candidate) == q {
                    recovered = Some(candidate);
                    via_target = Some(hit.target);
                    break 'search;
                }
                charges.failed_collisions += 1;
            }
            let mut left = [RawPoint::Infinity; PARALLEL_WALKS];
            let mut right = [RawPoint::Infinity; PARALLEL_WALKS];
            let mut active = [false; PARALLEL_WALKS];
            let mut jump_a = [0u64; PARALLEL_WALKS];
            let mut active_count = 0u64;
            for i in 0..PARALLEL_WALKS {
                if active_count >= step_cap - steps {
                    break;
                }
                let Some(lane) = lanes[i] else {
                    continue;
                };
                let jump = &jumps[raw_partition(lane.state.point)];
                charges.partition_hashes += 1;
                left[i] = lane.state.point;
                right[i] = jump.point;
                active[i] = true;
                jump_a[i] = jump.a;
                active_count += 1;
            }
            if active_count == 0 {
                continue;
            }
            let next_points = raw_add_pairwise(&curve, &left, &right, &active, &mut charges);
            packed_batches += 1;
            for i in 0..PARALLEL_WALKS {
                if !active[i] {
                    continue;
                }
                let lane = lanes[i].unwrap();
                let next = raw_canonicalize(
                    &curve,
                    RawState {
                        point: next_points[i],
                        a: (lane.state.a + jump_a[i]) % modulus,
                        b: lane.state.b,
                    },
                    mode,
                    modulus,
                    lambda,
                    &mut charges,
                );
                charges.group_additions += 1;
                steps += 1;
                let length = lane.length + 1;
                if next.point == lane.state.point || lane.previous.contains(&next.point) {
                    fruitless += 1;
                    lanes[i] = None;
                    continue;
                }
                if length > walk_cap {
                    capped += 1;
                    lanes[i] = None;
                    continue;
                }
                lanes[i] = Some(WalkLane {
                    state: next,
                    previous: [
                        lane.state.point,
                        lane.previous[0],
                        lane.previous[1],
                        lane.previous[2],
                    ],
                    length,
                });
            }
        }
        let solve_ms = started.elapsed().as_secs_f64() * 1000.0;
        let recovered = recovered.expect("batched rho exceeded the per-target step cap");
        if let Some(published) = d0 {
            assert_eq!(recovered, published);
        }
        let via = via_target.unwrap();
        if via != index {
            cross_solves += 1;
        }
        solved.push(recovered);
        total_steps += steps;
        let q_key = raw_key(q);
        println!(
            "{}",
            json!({
                "kind":"rho_ks_batch_fixture",
                "evidence_class":"measured_rho_observation",
                "n":n,"a":a,"quotient_mode":mode.name(),"automorphism_size":automorphisms,
                "fixture_index":index,"fixture_seed":seed,"batch_seed":batch_seed,
                "published_fixture_scalar":d0,"recovered_fixture_scalar":recovered,
                "published_q":[q_key.1,q_key.2],"verified":true,
                "target_source":if point_input.is_some() { "public_point_jsonl" } else { "generated_fixture" },
                "solved_via_target":via,"cross_target_solve":via != index,
                "walk_steps":steps,"walks":walks,"packed_batches":packed_batches,"fruitless_two_cycles":fruitless,
                "capped_walks":capped,"wasted_merges":wasted_merges,
                "table_entries_before":table_before,"table_entries_after":table.len(),
                "total_ms":solve_ms,
                "ideal_independent_steps":ideal_single,
            })
        );
    }
    let process_ms = process_started.elapsed().as_secs_f64() * 1000.0;
    println!(
        "{}",
        json!({
            "kind":"rho_ks_batch_summary","producer_version":"v3_parallel_walks_gf2",
            "n":n,"a":a,"quotient_mode":mode.name(),"automorphism_size":automorphisms,
            "fixtures":fixtures,"batch_seed":batch_seed,"corpus":shared_corpus,"dp_bits":dp_bits,"jump_count":JUMPS,
            "parallel_walks":PARALLEL_WALKS,
            "canonicalization_backend":canonicalization_backend,
            "target_source":if point_input.is_some() { "public_point_jsonl" } else { "generated_fixture" },
            "field_product_backend":curve.field.kernel_name(),
            "inversion_backend":"itoh_tsujii",
            "all_verified":true,"cross_target_solves":cross_solves,
            "total_walk_steps":total_steps,"table_entries":table.len(),
            "table_payload_lower_bound_bytes":table.len() * (1 + 4 * std::mem::size_of::<u64>()),
            "setup_ms":setup_ms,"in_process_ms":process_ms,
            "charges":{
                "group_additions":charges.group_additions,
                "scalar_multiplications":charges.scalar_multiplications,
                "canonicalizations":charges.canonicalizations,
                "partition_hashes":charges.partition_hashes,
                "table_queries":charges.table_queries,
                "table_inserts":charges.table_inserts,
                "failed_collisions":charges.failed_collisions,
                "batch_inversion_calls":charges.batch_inversion_calls,
                "batch_inversion_inputs":charges.batch_inversion_inputs,
                "batch_fallback_additions":charges.batch_fallback_additions,
                "normal_basis_transforms":charges.normal_basis_transforms,
                "canonical_rotation_steps":charges.canonical_rotation_steps
            },
            "scope":"published synthetic toy fixtures; no external point, unknown scalar, or production key"
        })
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runtime_carryless_product_matches_bit_serial_reference() {
        let mut state = 0x6a09_e667_f3bc_c909u64;
        for _ in 0..128 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let left = state;
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let right = state;
            assert_eq!(
                carryless_product(left, right),
                carryless_product_software(left, right)
            );
        }
    }

    #[test]
    fn x_first_canonicalization_matches_full_orbit_reference() {
        for n in [13, 37, 41, 53] {
            let curve = RhoCurve::new(0, n);
            let generator = raw_point(curve.generator());
            let modulus = curve.subgroup_order.to_u64_digits()[0];
            let lambda = curve.lambda.to_u64_digits()[0];
            let mut rng = StdRng::seed_from_u64(0x5be0_cd19_137e_2179 ^ u64::from(n));
            for mode in [Quotient::Negation, Quotient::SignedFrobenius] {
                for _ in 0..32 {
                    let state = RawState {
                        point: raw_scalar_mul(&curve, generator, rng.gen_range(1..modulus)),
                        a: rng.gen_range(0..modulus),
                        b: rng.gen_range(0..modulus),
                    };
                    let fast = raw_canonicalize(
                        &curve,
                        state,
                        mode,
                        modulus,
                        lambda,
                        &mut Charges::default(),
                    );
                    let reference = raw_canonicalize_reference(
                        &curve,
                        state,
                        mode,
                        modulus,
                        lambda,
                        &mut Charges::default(),
                    );
                    assert_eq!(raw_key(fast.point), raw_key(reference.point));
                    assert_eq!(fast.a, reference.a);
                    assert_eq!(fast.b, reference.b);
                }
                // The x=0 two-torsion point has a short x orbit and
                // exercises the complete-enumeration fallback.
                let exceptional = RawState {
                    point: RawPoint::Affine { x: 0, y: 1 },
                    a: 3,
                    b: 5,
                };
                let fast = raw_canonicalize(
                    &curve,
                    exceptional,
                    mode,
                    modulus,
                    lambda,
                    &mut Charges::default(),
                );
                let reference = raw_canonicalize_reference(
                    &curve,
                    exceptional,
                    mode,
                    modulus,
                    lambda,
                    &mut Charges::default(),
                );
                assert_eq!(raw_key(fast.point), raw_key(reference.point));
                assert_eq!(fast.a, reference.a);
                assert_eq!(fast.b, reference.b);
            }
        }
    }

    #[test]
    fn normal_basis_rotation_matches_frobenius_and_roundtrips() {
        for n in [13, 37, 41, 53] {
            let curve = RhoCurve::with_normal_basis(0, n, true);
            let basis = curve.normal_basis.as_ref().unwrap();
            let mask = (1u64 << n) - 1;
            let mut rng = StdRng::seed_from_u64(0x428a_2f98_d728_ae22 ^ u64::from(n));
            for _ in 0..128 {
                let value = rng.gen::<u64>() & mask;
                let normal = basis.to_normal.apply(value);
                assert_eq!(basis.to_poly.apply(normal), value);
                assert_eq!(
                    basis.to_poly.apply(basis.rotate(normal, 1)),
                    raw_square(&curve, value)
                );
            }
        }
    }

    #[test]
    fn normal_canonicalization_keeps_point_and_coefficients_in_signed_orbit() {
        for n in [13, 37, 41, 53] {
            let curve = RhoCurve::with_normal_basis(0, n, true);
            let generator = raw_point(curve.generator());
            let modulus = curve.subgroup_order.to_u64_digits()[0];
            let lambda = curve.lambda.to_u64_digits()[0];
            let mut rng = StdRng::seed_from_u64(0x7137_4491_23ef_65cd ^ u64::from(n));
            for _ in 0..32 {
                let state = RawState {
                    point: raw_scalar_mul(&curve, generator, rng.gen_range(1..modulus)),
                    a: rng.gen_range(0..modulus),
                    b: rng.gen_range(0..modulus),
                };
                let candidate = raw_canonicalize(
                    &curve,
                    state,
                    Quotient::SignedFrobenius,
                    modulus,
                    lambda,
                    &mut Charges::default(),
                );
                let again = raw_canonicalize(
                    &curve,
                    candidate,
                    Quotient::SignedFrobenius,
                    modulus,
                    lambda,
                    &mut Charges::default(),
                );
                assert_eq!(raw_key(candidate.point), raw_key(again.point));
                assert_eq!((candidate.a, candidate.b), (again.a, again.b));
                let mut point = state.point;
                let mut factor = 1u64;
                let mut found = false;
                for _ in 0..n {
                    for (orbit_point, orbit_factor) in
                        [(point, factor), (raw_neg(point), modulus - factor)]
                    {
                        if candidate.point == orbit_point
                            && candidate.a == mul_mod(state.a, orbit_factor, modulus)
                            && candidate.b == mul_mod(state.b, orbit_factor, modulus)
                        {
                            found = true;
                        }
                    }
                    point = match point {
                        RawPoint::Infinity => RawPoint::Infinity,
                        RawPoint::Affine { x, y } => RawPoint::Affine {
                            x: raw_square(&curve, x),
                            y: raw_square(&curve, y),
                        },
                    };
                    factor = mul_mod(factor, lambda, modulus);
                }
                assert!(found, "n={n}");
            }
        }
    }

    #[test]
    fn pairwise_batch_add_matches_complete_group_law() {
        for n in [13, 37, 41, 53] {
            let curve = RhoCurve::new(0, n);
            let generator = raw_point(curve.generator());
            let modulus = curve.subgroup_order.to_u64_digits()[0];
            let mut rng = StdRng::seed_from_u64(0x1f83_d9ab_fb41_bd6b ^ u64::from(n));
            for _ in 0..4 {
                let mut left = [RawPoint::Infinity; PARALLEL_WALKS];
                let mut right = [RawPoint::Infinity; PARALLEL_WALKS];
                let mut active = [true; PARALLEL_WALKS];
                for i in 0..PARALLEL_WALKS {
                    left[i] = raw_scalar_mul(&curve, generator, rng.gen_range(1..modulus));
                    right[i] = match i % 9 {
                        0 => RawPoint::Infinity,
                        1 => left[i],
                        2 => raw_neg(left[i]),
                        _ => raw_scalar_mul(&curve, generator, rng.gen_range(1..modulus)),
                    };
                    if i % 11 == 0 {
                        active[i] = false;
                    }
                }
                let mut charges = Charges::default();
                let batch = raw_add_pairwise(&curve, &left, &right, &active, &mut charges);
                for i in 0..PARALLEL_WALKS {
                    assert_eq!(
                        batch[i],
                        if active[i] {
                            raw_add(&curve, left[i], right[i])
                        } else {
                            RawPoint::Infinity
                        },
                        "n={n} lane={i}"
                    );
                }
                assert_eq!(charges.batch_inversion_calls, 1);
                assert!(charges.batch_inversion_inputs > 0);
                assert!(charges.batch_fallback_additions > 0);
            }
        }
    }

    #[test]
    fn itoh_inverse_matches_square_and_multiply_reference() {
        for n in [13, 37, 41, 53] {
            let curve = RhoCurve::new(0, n);
            let mask = (1u64 << n) - 1;
            let exponent = (1u64 << n) - 2;
            let mut state = 0x243f_6a88_85a3_08d3u64;
            for _ in 0..64 {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let value = (state & mask).max(1);
                let mut reference = 1u64;
                let mut base = value;
                for bit in 0..n {
                    if (exponent >> bit) & 1 == 1 {
                        reference = raw_mul_field(&curve, reference, base);
                    }
                    base = raw_square(&curve, base);
                }
                let candidate = raw_inverse(&curve, value);
                assert_eq!(candidate, reference, "n={n} value={value}");
                assert_eq!(raw_mul_field(&curve, value, candidate), 1);
            }
        }
    }
}
