#![recursion_limit = "256"]
//! Public-fixture Pollard-rho calibration for the Koblitz crossover task.
//!
//! This executable is deliberately limited to the exact toy rungs.  Every
//! discrete logarithm is a published synthetic fixture used to validate the
//! walk.  It never accepts an external point or a secret scalar.
//!
//! Backends (argument 5):
//!
//! * `strong` — the reference to measure index calculus against: the library
//!   [`crypto_lib::cryptanalysis::koblitz_strong_rho`] walk (distinguished
//!   points, normal-basis signed-Frobenius canonical form, library `Gf2`
//!   arithmetic, `KIC_RHO_LANES` lockstep walks with batched inversion,
//!   default 32; `KIC_RHO_DP_BITS`, default 4). `signed_frobenius` only.
//! * `reference`, `packed` — kept unchanged so archived stage runs reproduce.
//!   **Neither is a valid `vs_rho` reference**: both store every step in a
//!   table (no distinguished points) and canonicalize by an O(n)
//!   polynomial-basis scan; at n = 53 `packed` costs about 43× the
//!   instructions and 60× the wall time of `strong` on the same target
//!   (`docs/ic/BOUNDARY_TARGETS.md`, 2026-09-30 erratum). They print a note to
//!   stderr saying so.

use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{points_with_x, KoblitzCurve};
use crypto_lib::cryptanalysis::koblitz_strong_rho::{
    RawPoint as StrongPoint, StrongRho, StrongRhoCharges, StrongRhoParams,
};
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

const TASK_ID: &str = "TASK-KIC-SAT-RHO-CROSSOVER-20260909";
/// The exact toy rungs.  47, 57 and 61 are ledger §23's sizes the others miss,
/// so that the rule's comparison can check its rho against the strong walk at
/// every size it runs (`docs/ic/BOUNDARY_TARGETS.md`, 2026-10-01).
const RUNGS: [u32; 13] = [7, 11, 13, 17, 19, 23, 37, 41, 47, 53, 57, 59, 61];
const JUMPS: usize = 32;
const MAX_RESTARTS: u64 = 128;
/// Fruitless-collision restart budget for larger fields (n≥41).
const MAX_RESTARTS_LARGE: u64 = 100_000;

fn peak_rss_bytes() -> Option<u64> {
    #[cfg(unix)]
    {
        let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
        // SAFETY: getrusage initializes the complete value on success.
        if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
            return None;
        }
        // SAFETY: getrusage succeeded above.
        let rss = unsafe { usage.assume_init() }.ru_maxrss;
        if rss < 0 {
            return None;
        }
        #[cfg(target_os = "macos")]
        return Some(rss as u64);
        #[cfg(not(target_os = "macos"))]
        return (rss as u64).checked_mul(1024);
    }
    #[cfg(not(unix))]
    {
        None
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Quotient {
    Ordinary,
    Negation,
    SignedFrobenius,
}

#[derive(Clone, Copy)]
enum FixtureTarget {
    SeededScalar,
    ExplicitScalar(u64),
    PublicHash(u64),
}

impl FixtureTarget {
    fn parse(value: Option<&str>) -> Self {
        match value {
            None => Self::SeededScalar,
            Some(value) if value.starts_with("hash:") => Self::PublicHash(
                value[5..]
                    .parse()
                    .expect("public hash target seed must be a u64"),
            ),
            Some(value) => Self::ExplicitScalar(
                value
                    .parse()
                    .expect("explicit validation scalar must be a u64"),
            ),
        }
    }
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

#[derive(Clone)]
struct State {
    point: BinaryPoint,
    a: u64,
    b: u64,
}

#[derive(Clone)]
struct Jump {
    point: BinaryPoint,
    a: u64,
    b: u64,
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
}

fn point_key(point: &BinaryPoint) -> (u8, u64, u64) {
    match point {
        BinaryPoint::Infinity => (0, 0, 0),
        BinaryPoint::Affine { x, y } => (
            1,
            x.raw_bits().first().copied().unwrap_or(0),
            y.raw_bits().first().copied().unwrap_or(0),
        ),
    }
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
    raw_reduce(curve, carryless_product(left, right))
}

fn raw_inverse(curve: &KoblitzCurve, value: u64) -> u64 {
    assert_ne!(value, 0);
    let exponent = (1u64 << curve.n) - 2;
    let mut result = 1u64;
    let mut base = value;
    for bit in 0..curve.n {
        if (exponent >> bit) & 1 == 1 {
            result = raw_mul_field(curve, result, base);
        }
        base = raw_square(curve, base);
    }
    result
}

fn raw_neg(point: RawPoint) -> RawPoint {
    match point {
        RawPoint::Infinity => RawPoint::Infinity,
        RawPoint::Affine { x, y } => RawPoint::Affine { x, y: y ^ x },
    }
}

fn raw_double(curve: &KoblitzCurve, point: RawPoint) -> RawPoint {
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

fn raw_add(curve: &KoblitzCurve, left: RawPoint, right: RawPoint) -> RawPoint {
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

fn raw_scalar_mul(curve: &KoblitzCurve, point: RawPoint, scalar: u64) -> RawPoint {
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

fn canonicalize(
    curve: &KoblitzCurve,
    state: State,
    mode: Quotient,
    modulus: u64,
    lambda: u64,
    charges: &mut Charges,
) -> State {
    charges.canonicalizations += 1;
    if state.point == BinaryPoint::Infinity || mode == Quotient::Ordinary {
        return state;
    }

    let powers = if mode == Quotient::SignedFrobenius {
        curve.n
    } else {
        1
    };
    let mut point = state.point.clone();
    let mut multiplier = 1u64;
    let mut best_key = point_key(&point);
    let mut best_point = point.clone();
    let mut best_multiplier = multiplier;

    for exponent in 0..powers {
        let key = point_key(&point);
        if key < best_key {
            best_key = key;
            best_point = point.clone();
            best_multiplier = multiplier;
        }
        if mode.uses_negation() {
            charges.negations_examined += 1;
            let negative = point_neg(&point);
            let negative_key = point_key(&negative);
            if negative_key < best_key {
                best_key = negative_key;
                best_point = negative;
                best_multiplier = if multiplier == 0 {
                    0
                } else {
                    modulus - multiplier
                };
            }
        }
        if exponent + 1 < powers {
            point = curve.frobenius(&point);
            multiplier = mul_mod(multiplier, lambda, modulus);
            charges.frobenius_maps += 1;
        }
    }

    State {
        point: best_point,
        a: mul_mod(state.a, best_multiplier, modulus),
        b: mul_mod(state.b, best_multiplier, modulus),
    }
}

fn raw_canonicalize(
    curve: &KoblitzCurve,
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

fn partition(point: &BinaryPoint) -> usize {
    let (_, x, y) = point_key(point);
    let mut value = x ^ y.rotate_left(21) ^ 0x9e37_79b9_7f4a_7c15;
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^= value >> 31;
    value as usize % JUMPS
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

fn random_state(
    curve: &KoblitzCurve,
    q: &BinaryPoint,
    rng: &mut StdRng,
    modulus: u64,
    charges: &mut Charges,
) -> State {
    loop {
        let a = rng.gen_range(0..modulus);
        let b = rng.gen_range(0..modulus);
        charges.scalar_multiplications += 2;
        charges.group_additions += 1;
        let point = curve.add(
            &curve.mul(curve.generator(), &BigUint::from(a)),
            &curve.mul(q, &BigUint::from(b)),
        );
        if point != BinaryPoint::Infinity {
            return State { point, a, b };
        }
    }
}

fn make_jumps(
    curve: &KoblitzCurve,
    q: &BinaryPoint,
    rng: &mut StdRng,
    modulus: u64,
    charges: &mut Charges,
) -> Vec<Jump> {
    (0..JUMPS)
        .map(|_| loop {
            let state = random_state(curve, q, rng, modulus, charges);
            if state.point != BinaryPoint::Infinity {
                break Jump {
                    point: state.point,
                    a: state.a,
                    b: state.b,
                };
            }
        })
        .collect()
}

fn raw_random_state(
    curve: &KoblitzCurve,
    generator: RawPoint,
    q: RawPoint,
    rng: &mut StdRng,
    modulus: u64,
    charges: &mut Charges,
) -> RawState {
    loop {
        let a = rng.gen_range(0..modulus);
        let b = rng.gen_range(0..modulus);
        charges.scalar_multiplications += 2;
        charges.group_additions += 1;
        let point = raw_add(
            curve,
            raw_scalar_mul(curve, generator, a),
            raw_scalar_mul(curve, q, b),
        );
        if point != RawPoint::Infinity {
            return RawState { point, a, b };
        }
    }
}

fn raw_make_jumps(
    curve: &KoblitzCurve,
    generator: RawPoint,
    q: RawPoint,
    rng: &mut StdRng,
    modulus: u64,
    charges: &mut Charges,
) -> Vec<RawJump> {
    (0..JUMPS)
        .map(|_| {
            let state = raw_random_state(curve, generator, q, rng, modulus, charges);
            RawJump {
                point: state.point,
                a: state.a,
                b: state.b,
            }
        })
        .collect()
}

fn public_hash_target(curve: &KoblitzCurve, seed: u64) -> (BinaryPoint, u64) {
    const DOMAIN: &[u8] = b"ic-workflow-public-target-v1\0";
    let mask = (1u64 << curve.n) - 1;
    for counter in 0u64..1_000_000 {
        let mut hasher = blake3::Hasher::new();
        hasher.update(DOMAIN);
        hasher.update(&curve.n.to_le_bytes());
        hasher.update(&[curve.a]);
        hasher.update(&curve.k.to_le_bytes());
        hasher.update(&curve.b_index.to_le_bytes());
        hasher.update(&seed.to_le_bytes());
        hasher.update(&counter.to_le_bytes());
        let digest = hasher.finalize();
        let bytes = digest.as_bytes();
        let mut word = [0u8; 8];
        word.copy_from_slice(&bytes[..8]);
        let x = F2mElement::from_biguint(&BigUint::from(u64::from_le_bytes(word) & mask), curve.n);
        let mut lifts = points_with_x(&curve.curve, &x);
        lifts.sort_by_key(point_key);
        if lifts.is_empty() {
            continue;
        }
        let raw = lifts[usize::from(bytes[8] & 1) % lifts.len()].clone();
        let target = curve.mul(&raw, &curve.cofactor);
        if target == BinaryPoint::Infinity {
            continue;
        }
        assert_eq!(
            curve.mul(&target, &curve.subgroup_order),
            BinaryPoint::Infinity
        );
        return (target, counter);
    }
    panic!("public hash-to-curve target attempt cap exhausted");
}

fn solve_fixture(
    curve: &KoblitzCurve,
    mode: Quotient,
    fixture_index: u64,
    fixture_seed: u64,
    fixture_target: FixtureTarget,
) -> serde_json::Value {
    let modulus = curve.subgroup_order.to_u64_digits()[0];
    let lambda = curve.lambda.to_u64_digits()[0];
    let signed_size = signed_automorphism_size(lambda, modulus, curve.n);
    let fixture_seed = std::env::var("KIC_RHO_WALK_SEED")
        .ok()
        .map(|value| {
            value
                .parse::<u64>()
                .expect("KIC_RHO_WALK_SEED must be an integer")
        })
        .unwrap_or(fixture_seed);
    let fixed_target = std::env::var("KIC_RHO_FIXED_TARGET_SCALAR")
        .ok()
        .map(|value| {
            value
                .parse::<u64>()
                .expect("KIC_RHO_FIXED_TARGET_SCALAR must be an integer")
        });
    let mut rng = StdRng::seed_from_u64(fixture_seed);
    let d0 = fixed_target.unwrap_or_else(|| rng.gen_range(1..modulus));
    assert!(
        (1..modulus).contains(&d0),
        "fixed target scalar must be in [1,r)"
    );
    let generated_scalar = rng.gen_range(1..modulus);
    let target_generation_started = Instant::now();
    let (known_scalar, q, fixture_scalar_source, public_hash_seed, public_hash_counter) =
        match fixture_target {
            FixtureTarget::SeededScalar => (
                Some(generated_scalar),
                curve.mul(curve.generator(), &BigUint::from(generated_scalar)),
                "seeded_fixture_scalar",
                None,
                None,
            ),
            FixtureTarget::ExplicitScalar(scalar) => (
                Some(scalar),
                curve.mul(curve.generator(), &BigUint::from(scalar)),
                "explicit_public_validation_scalar",
                None,
                None,
            ),
            FixtureTarget::PublicHash(seed) => {
                let (target, counter) = public_hash_target(curve, seed);
                (
                    None,
                    target,
                    "public_hash_unknown_scalar",
                    Some(seed),
                    Some(counter),
                )
            }
        };
    let target_generation_ms = target_generation_started.elapsed().as_secs_f64() * 1000.0;
    let mut charges = Charges {
        scalar_multiplications: 1,
        ..Charges::default()
    };
    let started = Instant::now();
    let jumps = make_jumps(curve, &q, &mut rng, modulus, &mut charges);
    let setup_ms = started.elapsed().as_secs_f64() * 1000.0;
    let walk_started = Instant::now();
    let mut table: HashMap<(u8, u64, u64), (u64, u64)> = HashMap::new();
    let ideal_steps = (std::f64::consts::PI * modulus as f64
        / (2.0
            * match mode {
                Quotient::Ordinary => 1.0,
                Quotient::Negation => 2.0,
                Quotient::SignedFrobenius => signed_size as f64,
            }))
    .sqrt();
    // Toy rungs finish with 200×√(πr/2A); n≥41 needs more headroom —
    // distinguished-point density and restart waste grow with the field.
    let safety = if curve.n >= 41 { 2_000 } else { 200 };
    let max_steps = (ideal_steps.ceil() as u64)
        .saturating_mul(safety)
        .max(10_000);
    let mut steps = 0u64;
    let mut restarts = 0u64;
    let mut recovered = None;
    let restart_cap = if curve.n >= 41 {
        MAX_RESTARTS_LARGE
    } else {
        MAX_RESTARTS
    };

    'restart: while restarts <= restart_cap && steps < max_steps {
        let initial = random_state(curve, &q, &mut rng, modulus, &mut charges);
        let mut state = canonicalize(curve, initial, mode, modulus, lambda, &mut charges);
        loop {
            let key = point_key(&state.point);
            charges.table_queries += 1;
            if let Some(&(old_a, old_b)) = table.get(&key) {
                let numerator = sub_mod(old_a, state.a, modulus);
                let denominator = sub_mod(state.b, old_b, modulus);
                if let Some(inverse) = inverse_mod(denominator, modulus) {
                    let candidate = mul_mod(numerator, inverse, modulus);
                    charges.scalar_multiplications += 1;
                    if curve.mul(curve.generator(), &BigUint::from(candidate)) == q {
                        recovered = Some(candidate);
                        break 'restart;
                    }
                    charges.failed_collisions += 1;
                } else {
                    charges.failed_collisions += 1;
                }
                charges.fruitless_cycle_restarts += 1;
                restarts += 1;
                continue 'restart;
            }
            table.insert(key, (state.a, state.b));
            charges.table_inserts += 1;
            let jump_index = partition(&state.point);
            charges.partition_hashes += 1;
            let jump = &jumps[jump_index];
            state = State {
                point: curve.add(&state.point, &jump.point),
                a: (state.a + jump.a) % modulus,
                b: (state.b + jump.b) % modulus,
            };
            charges.group_additions += 1;
            state = canonicalize(curve, state, mode, modulus, lambda, &mut charges);
            steps += 1;
            if steps >= max_steps {
                break 'restart;
            }
        }
    }
    let walk_ms = walk_started.elapsed().as_secs_f64() * 1000.0;
    let recovered = recovered.expect("public rho fixture exceeded the frozen step/restart cap");
    if let Some(expected) = known_scalar {
        assert_eq!(recovered, expected);
    }
    let table_entries = table.len();
    let generator_point_key = point_key(curve.generator());
    let q_point_key = point_key(&q);

    json!({
        "schema_version":"1.0",
        "task_id":TASK_ID,
        "kind":"rho_public_fixture",
        "evidence_class":"measured_rho_observation",
        "n":curve.n,
        "a":curve.curve.a.raw_bits().first().copied().unwrap_or(0),
        "subgroup_order":modulus,
        "lambda":lambda,
        "quotient_mode":mode.name(),
        "arithmetic_backend":"reference_f2m",
        "automorphism_size":match mode {
            Quotient::Ordinary => 1,
            Quotient::Negation => 2,
            Quotient::SignedFrobenius => signed_size as u32,
        },
        "fixture_index":fixture_index,
        "fixture_seed":fixture_seed,
        "published_fixture_scalar":d0,
        "target_generation_ms_excluded":target_generation_ms,
        "published_fixture_scalar":known_scalar,
        "fixture_scalar_source":fixture_scalar_source,
        "target_scalar_constructed":known_scalar.is_some(),
        "target_kind":if known_scalar.is_some() {"known_scalar_multiple"} else {"public_hash_to_curve_cofactor"},
        "public_hash_seed":public_hash_seed,
        "public_hash_counter":public_hash_counter,
        "recovered_fixture_scalar":recovered,
        "generator":[generator_point_key.1,generator_point_key.2],
        "published_q":[q_point_key.1,q_point_key.2],
        "generator_point_key":[generator_point_key.1,generator_point_key.2],
        "published_q_point_key":[q_point_key.1,q_point_key.2],
        "field_modulus_low_terms":curve.curve.irreducible.low_terms,
        "verified":true,
        "ideal_steps":ideal_steps,
        "walk_steps":steps,
        "restarts":restarts,
        "jump_count":JUMPS,
        "distinguished_bits":0,
        "table_entries":table_entries,
        "table_payload_lower_bound_bytes":table_entries * (1 + 5 * std::mem::size_of::<u64>()),
        "setup_ms":setup_ms,
        "target_generation_ms":target_generation_ms,
        "walk_ms":walk_ms,
        "total_ms":target_generation_ms + setup_ms + walk_ms,
        "charges":{
            "group_additions":charges.group_additions,
            "scalar_multiplications":charges.scalar_multiplications,
            "canonicalizations":charges.canonicalizations,
            "frobenius_maps":charges.frobenius_maps,
            "negations_examined":charges.negations_examined,
            "partition_hashes":charges.partition_hashes,
            "table_queries":charges.table_queries,
            "table_inserts":charges.table_inserts,
            "failed_collisions":charges.failed_collisions,
            "fruitless_cycle_restarts":charges.fruitless_cycle_restarts
        },
        "scope":if known_scalar.is_some() {"published synthetic toy fixture; no external point or production key"} else {"public_hash_unknown_scalar"}
    })
}

fn solve_fixture_packed(
    curve: &KoblitzCurve,
    mode: Quotient,
    fixture_index: u64,
    fixture_seed: u64,
    fixture_target: FixtureTarget,
) -> serde_json::Value {
    let modulus = curve.subgroup_order.to_u64_digits()[0];
    let lambda = curve.lambda.to_u64_digits()[0];
    let signed_size = signed_automorphism_size(lambda, modulus, curve.n);
    let generator = raw_point(curve.generator());
    let fixture_seed = std::env::var("KIC_RHO_WALK_SEED")
        .ok()
        .map(|value| {
            value
                .parse::<u64>()
                .expect("KIC_RHO_WALK_SEED must be an integer")
        })
        .unwrap_or(fixture_seed);
    let fixed_target = std::env::var("KIC_RHO_FIXED_TARGET_SCALAR")
        .ok()
        .map(|value| {
            value
                .parse::<u64>()
                .expect("KIC_RHO_FIXED_TARGET_SCALAR must be an integer")
        });
    let mut rng = StdRng::seed_from_u64(fixture_seed);
    let d0 = fixed_target.unwrap_or_else(|| rng.gen_range(1..modulus));
    assert!(
        (1..modulus).contains(&d0),
        "fixed target scalar must be in [1,r)"
    );
    let generated_scalar = rng.gen_range(1..modulus);
    let target_generation_started = Instant::now();
    let (known_scalar, reference_q, fixture_scalar_source, public_hash_seed, public_hash_counter) =
        match fixture_target {
            FixtureTarget::SeededScalar => (
                Some(generated_scalar),
                curve.mul(curve.generator(), &BigUint::from(generated_scalar)),
                "seeded_fixture_scalar",
                None,
                None,
            ),
            FixtureTarget::ExplicitScalar(scalar) => (
                Some(scalar),
                curve.mul(curve.generator(), &BigUint::from(scalar)),
                "explicit_public_validation_scalar",
                None,
                None,
            ),
            FixtureTarget::PublicHash(seed) => {
                let (target, counter) = public_hash_target(curve, seed);
                (
                    None,
                    target,
                    "public_hash_unknown_scalar",
                    Some(seed),
                    Some(counter),
                )
            }
        };
    let target_generation_ms = target_generation_started.elapsed().as_secs_f64() * 1000.0;
    let mut charges = Charges::default();
    let q = raw_point(&reference_q);
    charges.scalar_multiplications += 1;
    let started = Instant::now();
    let jumps = raw_make_jumps(curve, generator, q, &mut rng, modulus, &mut charges);
    let setup_ms = started.elapsed().as_secs_f64() * 1000.0;
    let walk_started = Instant::now();
    let mut table: HashMap<(u8, u64, u64), (u64, u64)> = HashMap::new();
    let ideal_steps = (std::f64::consts::PI * modulus as f64
        / (2.0
            * match mode {
                Quotient::Ordinary => 1.0,
                Quotient::Negation => 2.0,
                Quotient::SignedFrobenius => signed_size as f64,
            }))
    .sqrt();
    // Toy rungs finish with 200×√(πr/2A); n≥41 needs more headroom —
    // distinguished-point density and restart waste grow with the field.
    let safety = if curve.n >= 41 { 2_000 } else { 200 };
    let max_steps = (ideal_steps.ceil() as u64)
        .saturating_mul(safety)
        .max(10_000);
    let mut steps = 0u64;
    let mut restarts = 0u64;
    let mut recovered = None;
    let restart_cap = if curve.n >= 41 {
        MAX_RESTARTS_LARGE
    } else {
        MAX_RESTARTS
    };

    'restart: while restarts <= restart_cap && steps < max_steps {
        let initial = raw_random_state(curve, generator, q, &mut rng, modulus, &mut charges);
        let mut state = raw_canonicalize(curve, initial, mode, modulus, lambda, &mut charges);
        loop {
            let key = raw_key(state.point);
            charges.table_queries += 1;
            if let Some(&(old_a, old_b)) = table.get(&key) {
                let numerator = sub_mod(old_a, state.a, modulus);
                let denominator = sub_mod(state.b, old_b, modulus);
                if let Some(inverse) = inverse_mod(denominator, modulus) {
                    let candidate = mul_mod(numerator, inverse, modulus);
                    charges.scalar_multiplications += 1;
                    if raw_scalar_mul(curve, generator, candidate) == q {
                        recovered = Some(candidate);
                        break 'restart;
                    }
                    charges.failed_collisions += 1;
                } else {
                    charges.failed_collisions += 1;
                }
                charges.fruitless_cycle_restarts += 1;
                restarts += 1;
                continue 'restart;
            }
            table.insert(key, (state.a, state.b));
            charges.table_inserts += 1;
            let jump = raw_make_jump_ref(&jumps, raw_partition(state.point));
            charges.partition_hashes += 1;
            state = RawState {
                point: raw_add(curve, state.point, jump.point),
                a: (state.a + jump.a) % modulus,
                b: (state.b + jump.b) % modulus,
            };
            charges.group_additions += 1;
            state = raw_canonicalize(curve, state, mode, modulus, lambda, &mut charges);
            steps += 1;
            if steps >= max_steps {
                break 'restart;
            }
        }
    }
    let walk_ms = walk_started.elapsed().as_secs_f64() * 1000.0;
    let recovered = recovered.expect("packed public rho fixture exceeded cap");
    let validation_started = Instant::now();
    if let Some(expected) = known_scalar {
        assert_eq!(recovered, expected);
    }
    assert_eq!(raw_point(&reference_q), q);
    assert_eq!(
        curve.mul(curve.generator(), &BigUint::from(recovered)),
        reference_q
    );
    let validation_ms = validation_started.elapsed().as_secs_f64() * 1000.0;
    let table_entries = table.len();
    let generator_point_key = raw_key(generator);
    let q_point_key = raw_key(q);

    json!({
        "schema_version":"1.0",
        "task_id":TASK_ID,
        "kind":"rho_public_fixture",
        "evidence_class":"measured_rho_observation",
        "n":curve.n,
        "a":curve.a,
        "subgroup_order":modulus,
        "lambda":lambda,
        "quotient_mode":mode.name(),
        "arithmetic_backend":"packed_u64_polynomial_basis",
        "automorphism_size":match mode {
            Quotient::Ordinary => 1,
            Quotient::Negation => 2,
            Quotient::SignedFrobenius => signed_size as u32,
        },
        "fixture_index":fixture_index,
        "fixture_seed":fixture_seed,
        "published_fixture_scalar":d0,
        "target_generation_ms_excluded":target_generation_ms,
        "published_fixture_scalar":known_scalar,
        "fixture_scalar_source":fixture_scalar_source,
        "target_scalar_constructed":known_scalar.is_some(),
        "target_kind":if known_scalar.is_some() {"known_scalar_multiple"} else {"public_hash_to_curve_cofactor"},
        "public_hash_seed":public_hash_seed,
        "public_hash_counter":public_hash_counter,
        "recovered_fixture_scalar":recovered,
        "generator":[generator_point_key.1,generator_point_key.2],
        "published_q":[q_point_key.1,q_point_key.2],
        "generator_point_key":[generator_point_key.1,generator_point_key.2],
        "published_q_point_key":[q_point_key.1,q_point_key.2],
        "field_modulus_low_terms":curve.curve.irreducible.low_terms,
        "verified":true,
        "reference_group_validation":true,
        "ideal_steps":ideal_steps,
        "walk_steps":steps,
        "restarts":restarts,
        "jump_count":JUMPS,
        "distinguished_bits":0,
        "table_entries":table_entries,
        "table_payload_lower_bound_bytes":table_entries * (1 + 5 * std::mem::size_of::<u64>()),
        "setup_ms":setup_ms,
        "target_generation_ms":target_generation_ms,
        "walk_ms":walk_ms,
        "validation_ms":validation_ms,
        "total_ms":target_generation_ms + setup_ms + walk_ms + validation_ms,
        "charges":{
            "group_additions":charges.group_additions,
            "scalar_multiplications":charges.scalar_multiplications,
            "canonicalizations":charges.canonicalizations,
            "frobenius_maps":charges.frobenius_maps,
            "negations_examined":charges.negations_examined,
            "partition_hashes":charges.partition_hashes,
            "table_queries":charges.table_queries,
            "table_inserts":charges.table_inserts,
            "failed_collisions":charges.failed_collisions,
            "fruitless_cycle_restarts":charges.fruitless_cycle_restarts
        },
        "scope":if known_scalar.is_some() {"published synthetic toy fixture; no external point or production key"} else {"public_hash_unknown_scalar"}
    })
}

/// `KIC_RHO_LANES` / `KIC_RHO_DP_BITS` overrides of the strong defaults.
fn strong_params() -> StrongRhoParams {
    let mut params = StrongRhoParams::default();
    if let Ok(value) = std::env::var("KIC_RHO_LANES") {
        params.lanes = value.parse().expect("KIC_RHO_LANES must be an integer");
        assert!(params.lanes >= 1, "KIC_RHO_LANES must be at least 1");
    }
    if let Ok(value) = std::env::var("KIC_RHO_DP_BITS") {
        params.dp_bits = value.parse().expect("KIC_RHO_DP_BITS must be an integer");
        assert!(params.dp_bits < 32, "KIC_RHO_DP_BITS must be below 32");
    }
    params
}

/// The `strong` backend: one target, same target derivation and JSON schema as
/// `solve_fixture_packed`. The jump table is drawn from
/// `blake3("KIC-KS-BATCH-JUMPS-v1|n|a|batch_seed")` and the walk starts from the
/// fixture RNG after the fixture scalar, exactly as
/// `examples/koblitz_rho_batch_ks_strong.rs` does at rung 3 with one fixture, so
/// the two walk the same trajectory for the same seeds.
fn solve_fixture_strong(
    curve: &KoblitzCurve,
    mode: Quotient,
    fixture_index: u64,
    fixture_seed: u64,
    fixture_target: FixtureTarget,
    batch_seed: Option<u64>,
) -> serde_json::Value {
    assert!(
        mode == Quotient::SignedFrobenius,
        "the strong backend implements signed_frobenius only"
    );
    let params = strong_params();
    let setup_started = Instant::now();
    let rho = StrongRho::new(curve);
    let modulus = rho.modulus();
    let lambda = rho.lambda();
    let precompute_ms = setup_started.elapsed().as_secs_f64() * 1000.0;
    let mut rng = StdRng::seed_from_u64(fixture_seed);
    let generated_scalar = rng.gen_range(1..modulus);
    let target_generation_started = Instant::now();
    let (known_scalar, reference_q, fixture_scalar_source, public_hash_seed, public_hash_counter) =
        match fixture_target {
            FixtureTarget::SeededScalar => (
                Some(generated_scalar),
                curve.mul(curve.generator(), &BigUint::from(generated_scalar)),
                "seeded_fixture_scalar",
                None,
                None,
            ),
            FixtureTarget::ExplicitScalar(scalar) => (
                Some(scalar),
                curve.mul(curve.generator(), &BigUint::from(scalar)),
                "explicit_public_validation_scalar",
                None,
                None,
            ),
            FixtureTarget::PublicHash(seed) => {
                let (target, counter) = public_hash_target(curve, seed);
                (
                    None,
                    target,
                    "public_hash_unknown_scalar",
                    Some(seed),
                    Some(counter),
                )
            }
        };
    let target_generation_ms = target_generation_started.elapsed().as_secs_f64() * 1000.0;
    let started = Instant::now();
    let q = StrongPoint::from_binary(&reference_q);
    let mut charges = StrongRhoCharges {
        scalar_multiplications: 1,
        ..StrongRhoCharges::default()
    };
    let jump_material = match batch_seed {
        Some(seed) => format!("KIC-KS-BATCH-JUMPS-v1|{}|{}|{seed}", curve.n, curve.a),
        None => format!("KIC-KS-BATCH-JUMPS-v1|{}|{}|{TASK_ID}", curve.n, curve.a),
    };
    let digest = blake3::hash(jump_material.as_bytes());
    let jump_seed = u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap());
    let jumps = rho.jumps(jump_seed, &mut charges);
    let setup_ms = precompute_ms + started.elapsed().as_secs_f64() * 1000.0;
    let walk_started = Instant::now();
    let outcome = rho
        .solve(q, &jumps, &mut rng, &params, charges)
        .expect("strong public rho fixture exceeded its step cap");
    let walk_ms = walk_started.elapsed().as_secs_f64() * 1000.0;
    let recovered = outcome.scalar;
    let validation_started = Instant::now();
    if let Some(expected) = known_scalar {
        assert_eq!(recovered, expected);
    }
    assert_eq!(
        curve.mul(curve.generator(), &BigUint::from(recovered)),
        reference_q
    );
    let validation_ms = validation_started.elapsed().as_secs_f64() * 1000.0;
    let generator_point_key = raw_key(raw_point(curve.generator()));
    let q_point_key = raw_key(raw_point(&reference_q));
    let c = outcome.charges;

    json!({
        "schema_version":"1.0",
        "task_id":TASK_ID,
        "kind":"rho_public_fixture",
        "evidence_class":"measured_rho_observation",
        "reference_grade":"strong",
        "n":curve.n,
        "a":curve.a,
        "subgroup_order":modulus,
        "lambda":lambda,
        "quotient_mode":mode.name(),
        "arithmetic_backend":"gf2_normal_basis_lockstep_lanes",
        "automorphism_size":outcome.automorphisms,
        "fixture_index":fixture_index,
        "fixture_seed":fixture_seed,
        "batch_seed":batch_seed,
        "published_fixture_scalar":known_scalar,
        "fixture_scalar_source":fixture_scalar_source,
        "target_scalar_constructed":known_scalar.is_some(),
        "target_kind":if known_scalar.is_some() {"known_scalar_multiple"} else {"public_hash_to_curve_cofactor"},
        "public_hash_seed":public_hash_seed,
        "public_hash_counter":public_hash_counter,
        "recovered_fixture_scalar":recovered,
        "generator":[generator_point_key.1,generator_point_key.2],
        "published_q":[q_point_key.1,q_point_key.2],
        "generator_point_key":[generator_point_key.1,generator_point_key.2],
        "published_q_point_key":[q_point_key.1,q_point_key.2],
        "field_modulus_low_terms":curve.curve.irreducible.low_terms,
        "verified":true,
        "reference_group_validation":true,
        "ideal_steps":outcome.ideal_steps,
        "walk_steps":outcome.walk_steps,
        "walks":outcome.walks,
        "restarts":outcome.fruitless + outcome.capped,
        "fruitless_cycles":outcome.fruitless,
        "capped_walks":outcome.capped,
        "wasted_merges":outcome.wasted_merges,
        "lanes":params.lanes,
        "jump_count":JUMPS,
        "distinguished_bits":params.dp_bits,
        "table_entries":outcome.table_entries,
        "table_payload_lower_bound_bytes":outcome.table_entries * (1 + 4 * std::mem::size_of::<u64>()),
        "peak_rss_bytes":peak_rss_bytes(),
        "setup_ms":setup_ms,
        "target_generation_ms":target_generation_ms,
        "walk_ms":walk_ms,
        "validation_ms":validation_ms,
        "total_ms":target_generation_ms + setup_ms + walk_ms + validation_ms,
        "charges":{
            "group_additions":c.group_additions,
            "scalar_multiplications":c.scalar_multiplications,
            "canonicalizations":c.canonicalizations,
            "frobenius_maps":0,
            "negations_examined":0,
            "partition_hashes":c.partition_hashes,
            "table_queries":c.table_queries,
            "table_inserts":c.table_inserts,
            "failed_collisions":c.failed_collisions,
            "fruitless_cycle_restarts":outcome.fruitless
        },
        "scope":if known_scalar.is_some() {"published synthetic toy fixture; no external point or production key"} else {"public_hash_unknown_scalar"}
    })
}

fn raw_make_jump_ref(jumps: &[RawJump], index: usize) -> &RawJump {
    &jumps[index]
}

// ── Wide (`u128`-word) packed rho backend for `64 < n ≤ 127` ─────
//
// Twin of the `u64` packed path with identical walk semantics; the
// `u64` code above is deliberately untouched so every `n ≤ 63`
// fixture stays byte-identical.  Raw point operations reuse the
// library's [`FastBinaryCurve128`] (itself group-law-tested against
// the general implementation); subgroup scalars stay `u64` (every
// admitted `r < 2^64` through `n = 127`); only field words widen.
mod wide {
    use super::{
        inverse_mod, mul_mod, signed_automorphism_size, sub_mod, Charges, Quotient, JUMPS,
        MAX_RESTARTS, MAX_RESTARTS_LARGE, TASK_ID,
    };
    use crypto_lib::binary_ecc::BinaryPoint;
    use crypto_lib::cryptanalysis::koblitz_fast_arith::{FastBinaryCurve128, FastPoint128};
    use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
    use crypto_lib::cryptanalysis::semaev_decomp::Gf2_128;
    use num_bigint::BigUint;
    use rand::{rngs::StdRng, Rng, SeedableRng};
    use std::collections::HashMap;
    use std::time::Instant;

    /// Raw point as two wide field words; `None` is infinity.
    pub(crate) type RawPoint128 = FastPoint128;

    #[derive(Clone, Copy)]
    pub(crate) struct RawState128 {
        pub(crate) point: RawPoint128,
        pub(crate) a: u64,
        pub(crate) b: u64,
    }

    #[derive(Clone, Copy)]
    pub(crate) struct RawJump128 {
        pub(crate) point: RawPoint128,
        pub(crate) a: u64,
        pub(crate) b: u64,
    }

    pub(crate) struct WideBackend {
        pub(crate) fast: FastBinaryCurve128,
        pub(crate) gf: Gf2_128,
        mask: u128,
    }

    impl WideBackend {
        pub(crate) fn new(curve: &KoblitzCurve) -> Option<Self> {
            let n = curve.n;
            if n <= 63 || n > 127 {
                return None;
            }
            let fast = FastBinaryCurve128::new(
                &curve.curve.irreducible,
                curve.curve.a.raw_bits().first().copied().unwrap_or(0) as u128,
            )?;
            let gf = Gf2_128::new(&curve.curve.irreducible);
            Some(Self {
                fast,
                gf,
                mask: (1u128 << n) - 1,
            })
        }

        #[inline(always)]
        pub(crate) fn norm(&self, w: u128) -> u128 {
            w & self.mask
        }
    }

    pub(crate) fn raw_point128(backend: &WideBackend, point: &BinaryPoint) -> RawPoint128 {
        match point {
            BinaryPoint::Infinity => None,
            BinaryPoint::Affine { x, y } => {
                let words = |e: &crypto_lib::binary_ecc::F2mElement| {
                    let raw = e.raw_bits();
                    backend.norm(
                        raw.first().copied().unwrap_or(0) as u128
                            | ((raw.get(1).copied().unwrap_or(0) as u128) << 64),
                    )
                };
                Some((words(x), words(y)))
            }
        }
    }

    pub(crate) fn raw_key128(point: RawPoint128) -> (u8, u128, u128) {
        match point {
            None => (0, 0, 0),
            Some((x, y)) => (1, x, y),
        }
    }

    fn point_coordinates_json128(key: (u8, u128, u128)) -> [String; 2] {
        [key.1.to_string(), key.2.to_string()]
    }

    fn parse_point_coordinates128(encoded: &str) -> (u128, u128) {
        let [x, y]: [String; 2] =
            serde_json::from_str(encoded).expect("public point must be a string [x,y] array");
        (
            x.parse().expect("public point x must be decimal"),
            y.parse().expect("public point y must be decimal"),
        )
    }

    pub(crate) fn raw_neg128(point: RawPoint128) -> RawPoint128 {
        FastBinaryCurve128::neg(point)
    }

    pub(crate) fn raw_double128(backend: &WideBackend, point: RawPoint128) -> RawPoint128 {
        match point {
            None => None,
            Some((x, y)) => backend.fast.double(x, y),
        }
    }

    pub(crate) fn raw_add128(
        backend: &WideBackend,
        left: RawPoint128,
        right: RawPoint128,
    ) -> RawPoint128 {
        backend.fast.add(left, right)
    }

    pub(crate) fn raw_scalar_mul128(
        backend: &WideBackend,
        point: RawPoint128,
        scalar: u64,
    ) -> RawPoint128 {
        let mut result = None;
        for bit in (0..64 - scalar.leading_zeros()).rev() {
            result = raw_double128(backend, result);
            if (scalar >> bit) & 1 == 1 {
                result = raw_add128(backend, result, point);
            }
        }
        result
    }

    pub(crate) fn raw_canonicalize128(
        backend: &WideBackend,
        curve: &KoblitzCurve,
        state: RawState128,
        mode: Quotient,
        modulus: u64,
        lambda: u64,
        charges: &mut Charges,
    ) -> RawState128 {
        charges.canonicalizations += 1;
        if state.point.is_none() || mode == Quotient::Ordinary {
            return state;
        }
        let powers = if mode == Quotient::SignedFrobenius {
            curve.n
        } else {
            1
        };
        let mut point = state.point;
        let mut multiplier = 1u64;
        let mut best_key = raw_key128(point);
        let mut best_point = point;
        let mut best_multiplier = multiplier;
        for exponent in 0..powers {
            let key = raw_key128(point);
            if key < best_key {
                best_key = key;
                best_point = point;
                best_multiplier = multiplier;
            }
            if mode.uses_negation() {
                charges.negations_examined += 1;
                let negative = raw_neg128(point);
                if raw_key128(negative) < best_key {
                    best_key = raw_key128(negative);
                    best_point = negative;
                    best_multiplier = modulus - multiplier;
                }
            }
            if exponent + 1 < powers {
                point = point.map(|(x, y)| (backend.gf.sqr(x), backend.gf.sqr(y)));
                multiplier = mul_mod(multiplier, lambda, modulus);
                charges.frobenius_maps += 1;
            }
        }
        RawState128 {
            point: best_point,
            a: mul_mod(state.a, best_multiplier, modulus),
            b: mul_mod(state.b, best_multiplier, modulus),
        }
    }

    pub(crate) fn raw_partition128(point: RawPoint128) -> usize {
        let (_, x, y) = raw_key128(point);
        let mut value = x ^ y.rotate_left(21) ^ 0x9E37_79B9_7F4A_7C15_9E37_79B9_7F4A_7C15;
        value ^= value >> 30;
        value = value.wrapping_mul(0xBF58_476D_1CE4_E5B9_BF58_476D_1CE4_E5B9);
        value ^= value >> 27;
        value = value.wrapping_mul(0x94D0_49BB_1331_11EB_94D0_49BB_1331_11EB);
        value ^= value >> 31;
        value as usize % JUMPS
    }

    pub(crate) fn raw_random_state128(
        backend: &WideBackend,
        _curve: &KoblitzCurve,
        generator: RawPoint128,
        q: RawPoint128,
        rng: &mut StdRng,
        modulus: u64,
        charges: &mut Charges,
    ) -> RawState128 {
        loop {
            let a = rng.gen_range(0..modulus);
            let b = rng.gen_range(0..modulus);
            charges.scalar_multiplications += 2;
            charges.group_additions += 1;
            let point = raw_add128(
                backend,
                raw_scalar_mul128(backend, generator, a),
                raw_scalar_mul128(backend, q, b),
            );
            if point.is_some() {
                return RawState128 { point, a, b };
            }
        }
    }

    pub(crate) fn raw_make_jumps128(
        backend: &WideBackend,
        curve: &KoblitzCurve,
        generator: RawPoint128,
        q: RawPoint128,
        rng: &mut StdRng,
        modulus: u64,
        charges: &mut Charges,
    ) -> Vec<RawJump128> {
        (0..JUMPS)
            .map(|_| {
                let state =
                    raw_random_state128(backend, curve, generator, q, rng, modulus, charges);
                RawJump128 {
                    point: state.point,
                    a: state.a,
                    b: state.b,
                }
            })
            .collect()
    }

    /// Deterministic per-fixture walk seed, same material shape as the
    /// `u64` flow's non-batch seed derivation.
    pub fn fixture_seed(n: u32, a: u8, mode_name: &str, fixture_index: u64) -> u64 {
        let material = format!("{TASK_ID}|rho|{n}|{a}|{mode_name}|{fixture_index}");
        let digest = blake3::hash(material.as_bytes());
        u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap())
    }

    /// Wide twin of `solve_fixture_packed`: same walk, caps, charges,
    /// and JSON shapes (field coordinates widen to `u128`); the
    /// recovered scalar is validated against the general path.
    pub fn solve_fixture_packed128(
        curve: &KoblitzCurve,
        mode: Quotient,
        fixture_index: u64,
        fixture_seed: u64,
    ) -> serde_json::Value {
        let backend = WideBackend::new(curve).expect("wide backend needs 64 < n ≤ 127");
        let modulus = curve.subgroup_order.to_u64_digits()[0];
        let lambda = curve.lambda.to_u64_digits()[0];
        let signed_size = signed_automorphism_size(lambda, modulus, curve.n);
        let generator = raw_point128(&backend, curve.generator());
        let fixture_seed = std::env::var("KIC_RHO_WALK_SEED")
            .ok()
            .map(|value| {
                value
                    .parse::<u64>()
                    .expect("KIC_RHO_WALK_SEED must be an integer")
            })
            .unwrap_or(fixture_seed);
        let fixed_target = std::env::var("KIC_RHO_FIXED_TARGET_SCALAR")
            .ok()
            .map(|value| {
                value
                    .parse::<u64>()
                    .expect("KIC_RHO_FIXED_TARGET_SCALAR must be an integer")
            });
        let public_target = std::env::var("KIC_RHO_PUBLIC_TARGET_POINT")
            .ok()
            .map(|encoded| parse_point_coordinates128(&encoded));
        assert!(
            public_target.is_none() || fixed_target.is_none(),
            "a public point run must not receive a known-answer scalar"
        );
        let mut rng = StdRng::seed_from_u64(fixture_seed);
        let d0 = if public_target.is_some() {
            None
        } else {
            Some(fixed_target.unwrap_or_else(|| rng.gen_range(1..modulus)))
        };
        if let Some(d0) = d0 {
            assert!(
                (1..modulus).contains(&d0),
                "fixed target scalar must be in [1,r)"
            );
        }
        let mut charges = Charges::default();
        let (q, target_generation_ms) = if let Some((x, y)) = public_target {
            (Some((x, y)), 0.0)
        } else {
            let target_generation_started = Instant::now();
            let q = raw_scalar_mul128(
                &backend,
                generator,
                d0.expect("generated target needs a scalar"),
            );
            let elapsed = target_generation_started.elapsed().as_secs_f64() * 1000.0;
            charges.scalar_multiplications += 1;
            (q, elapsed)
        };
        let started = Instant::now();
        let jumps = raw_make_jumps128(
            &backend,
            curve,
            generator,
            q,
            &mut rng,
            modulus,
            &mut charges,
        );
        let setup_ms = started.elapsed().as_secs_f64() * 1000.0;
        let walk_started = Instant::now();
        let mut table: HashMap<(u8, u128, u128), (u64, u64)> = HashMap::new();
        let ideal_steps = (std::f64::consts::PI * modulus as f64
            / (2.0
                * match mode {
                    Quotient::Ordinary => 1.0,
                    Quotient::Negation => 2.0,
                    Quotient::SignedFrobenius => signed_size as f64,
                }))
        .sqrt();
        let safety = if curve.n >= 41 { 2_000 } else { 200 };
        let max_steps = (ideal_steps.ceil() as u64)
            .saturating_mul(safety)
            .max(10_000);
        let mut steps = 0u64;
        let mut restarts = 0u64;
        let mut recovered = None;
        let restart_cap = if curve.n >= 41 {
            MAX_RESTARTS_LARGE
        } else {
            MAX_RESTARTS
        };

        'restart: while restarts <= restart_cap && steps < max_steps {
            let initial = raw_random_state128(
                &backend,
                curve,
                generator,
                q,
                &mut rng,
                modulus,
                &mut charges,
            );
            let mut state = raw_canonicalize128(
                &backend,
                curve,
                initial,
                mode,
                modulus,
                lambda,
                &mut charges,
            );
            loop {
                let key = raw_key128(state.point);
                charges.table_queries += 1;
                if let Some(&(old_a, old_b)) = table.get(&key) {
                    let numerator = sub_mod(old_a, state.a, modulus);
                    let denominator = sub_mod(state.b, old_b, modulus);
                    if let Some(inverse) = inverse_mod(denominator, modulus) {
                        let candidate = mul_mod(numerator, inverse, modulus);
                        charges.scalar_multiplications += 1;
                        if raw_scalar_mul128(&backend, generator, candidate) == q {
                            recovered = Some(candidate);
                            break 'restart;
                        }
                        charges.failed_collisions += 1;
                    } else {
                        charges.failed_collisions += 1;
                    }
                    charges.fruitless_cycle_restarts += 1;
                    restarts += 1;
                    continue 'restart;
                }
                table.insert(key, (state.a, state.b));
                charges.table_inserts += 1;
                let jump = &jumps[raw_partition128(state.point)];
                charges.partition_hashes += 1;
                state = RawState128 {
                    point: raw_add128(&backend, state.point, jump.point),
                    a: (state.a + jump.a) % modulus,
                    b: (state.b + jump.b) % modulus,
                };
                charges.group_additions += 1;
                state = raw_canonicalize128(
                    &backend,
                    curve,
                    state,
                    mode,
                    modulus,
                    lambda,
                    &mut charges,
                );
                steps += 1;
                if steps >= max_steps {
                    break 'restart;
                }
            }
        }
        let walk_ms = walk_started.elapsed().as_secs_f64() * 1000.0;
        let recovered = recovered.expect("wide packed public rho fixture exceeded cap");
        let validation_started = Instant::now();
        if let Some(d0) = d0 {
            assert_eq!(recovered, d0);
        }
        let reference_q = curve.mul(curve.generator(), &BigUint::from(recovered));
        assert_eq!(raw_point128(&backend, &reference_q), q);
        assert_eq!(
            curve.mul(curve.generator(), &BigUint::from(recovered)),
            reference_q
        );
        let validation_ms = validation_started.elapsed().as_secs_f64() * 1000.0;
        let table_entries = table.len();
        let generator_point_key = raw_key128(generator);
        let q_point_key = raw_key128(q);

        serde_json::json!({
            "schema_version":"1.0",
            "task_id":TASK_ID,
            "kind":"rho_public_fixture",
            "evidence_class":"measured_rho_observation",
            "n":curve.n,
            "a":curve.curve.a.raw_bits().first().copied().unwrap_or(0),
            "subgroup_order":modulus,
            "lambda":lambda,
            "quotient_mode":mode.name(),
            "arithmetic_backend":"packed_u128_polynomial_basis",
            "automorphism_size":match mode {
                Quotient::Ordinary => 1,
                Quotient::Negation => 2,
                Quotient::SignedFrobenius => signed_size as u32,
            },
            "fixture_index":fixture_index,
            "fixture_seed":fixture_seed,
            "published_fixture_scalar":d0,
            "target_input_kind":if public_target.is_some() { "public_point" } else { "fixture_scalar_generated_before_online" },
            "target_generation_ms_excluded":target_generation_ms,
            "recovered_fixture_scalar":recovered,
            "field_value_encoding":"base-10 strings for u128 field values",
            "generator":point_coordinates_json128(generator_point_key),
            "published_q":point_coordinates_json128(q_point_key),
            "generator_point_key":point_coordinates_json128(generator_point_key),
            "published_q_point_key":point_coordinates_json128(q_point_key),
            "field_modulus_low_terms":curve.curve.irreducible.low_terms,
            "verified":true,
            "reference_group_validation":true,
            "ideal_steps":ideal_steps,
            "walk_steps":steps,
            "restarts":restarts,
            "jump_count":JUMPS,
            "distinguished_bits":0,
            "table_entries":table_entries,
            "table_payload_lower_bound_bytes":table_entries * (1 + 2 * std::mem::size_of::<u128>() + 2 * std::mem::size_of::<u64>()),
            "setup_ms":setup_ms,
            "walk_ms":walk_ms,
            "validation_ms":validation_ms,
            "total_ms":target_generation_ms + setup_ms + walk_ms + validation_ms,
            "charges":{
                "group_additions":charges.group_additions,
                "scalar_multiplications":charges.scalar_multiplications,
                "canonicalizations":charges.canonicalizations,
                "frobenius_maps":charges.frobenius_maps,
                "negations_examined":charges.negations_examined,
                "partition_hashes":charges.partition_hashes,
                "table_queries":charges.table_queries,
                "table_inserts":charges.table_inserts,
                "failed_collisions":charges.failed_collisions,
                "fruitless_cycle_restarts":charges.fruitless_cycle_restarts
            },
            "scope":"published synthetic toy fixture; no external point, unknown scalar, or production key"
        })
    }

    #[cfg(test)]
    mod wide_json_tests {
        use super::{parse_point_coordinates128, point_coordinates_json128};

        #[test]
        fn wide_point_key_serializes_as_decimal_strings() {
            assert_eq!(
                point_coordinates_json128((1, (1u128 << 100) + 7, (1u128 << 80) + 9)),
                [
                    ((1u128 << 100) + 7).to_string(),
                    ((1u128 << 80) + 9).to_string()
                ]
            );
        }

        #[test]
        fn public_point_input_parses_decimal_string_coordinates() {
            assert_eq!(
                parse_point_coordinates128(
                    "[\"1267650600228229401496703205383\",\"1208925819614629174706185\"]"
                ),
                ((1u128 << 100) + 7, (1u128 << 80) + 9)
            );
        }
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert!(
        (5..=8).contains(&args.len()),
        "usage: <n> <a> <mode> <fixtures> [reference|packed|strong] [batch_seed] [explicit_fixture_scalar|hash:public_seed]"
    );
    let n: u32 = args[1].parse().unwrap();
    let a: u8 = args[2].parse().unwrap();
    let mode = Quotient::parse(&args[3]);
    let fixtures: u64 = args[4].parse().unwrap();
    let backend = args.get(5).map(String::as_str).unwrap_or("reference");
    let batch_seed = args.get(6).map(|value| value.parse::<u64>().unwrap());
    // Wide packed backend (`u128` words) for 64 < n ≤ 127.  The `u64`
    // flow below stays byte-identical; the wide path only implements
    // the packed walk (the beat control) and re-parses its own args.
    if n > 63 {
        assert_eq!(
            backend, "packed",
            "wide rungs only implement the packed walk"
        );
        let curve = KoblitzCurve::new(a, n).expect("frozen exact rung must construct");
        for fixture_index in 0..fixtures {
            // Same seed-derivation discipline as the `u64` batch path.
            let seed = match batch_seed {
                Some(bs) => {
                    let material = format!(
                        "TASK-KIC-DIRECT-BATCH-20260910|rho|{n}|{a}|{}|{bs}|{fixture_index}",
                        mode.name()
                    );
                    let digest = blake3::hash(material.as_bytes());
                    u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap())
                }
                None => wide::fixture_seed(n, a, mode.name(), fixture_index),
            };
            println!(
                "{}",
                wide::solve_fixture_packed128(&curve, mode, fixture_index, seed)
            );
        }
        return;
    }
    let shared_corpus = std::env::var("KIC_RHO_BATCH_CORPUS").ok();
    let shared_fixture_offset: u64 = std::env::var("KIC_RHO_FIXTURE_OFFSET")
        .ok()
        .map(|value| {
            value
                .parse()
                .expect("KIC_RHO_FIXTURE_OFFSET must be an integer")
        })
        .unwrap_or(0);
    assert!(
        shared_corpus.is_none() || batch_seed.is_some(),
        "KIC_RHO_BATCH_CORPUS requires the batch_seed argument"
    );
    assert!(matches!(n, 7 | 11 | 13 | 17 | 19 | 23 | 37 | 41 | 53 | 61));
    let fixture_target = FixtureTarget::parse(args.get(7).map(String::as_str));
    assert!(matches!(backend, "reference" | "packed" | "strong"));
    if backend != "strong" {
        eprintln!(
            "note: backend `{backend}` stores every step and canonicalizes by an O(n) scan; \
             it is not a valid vs_rho reference (docs/ic/BOUNDARY_TARGETS.md, 2026-09-30 \
             erratum). Use `strong`."
        );
    }
    assert!(
        RUNGS.contains(&n),
        "n = {n} is not one of the exact rungs {RUNGS:?}"
    );
    assert!(fixtures > 0);
    let curve = KoblitzCurve::new(a, n).expect("frozen exact rung must construct");
    let modulus = curve.subgroup_order.to_u64_digits()[0];
    if let FixtureTarget::ExplicitScalar(scalar) = fixture_target {
        assert!(fixtures == 1, "an explicit scalar requires one fixture");
        assert!(
            (1..modulus).contains(&scalar),
            "explicit scalar must be in 1..r"
        );
    }
    if matches!(fixture_target, FixtureTarget::PublicHash(_)) {
        assert!(fixtures == 1, "a public hash target requires one fixture");
    }
    for fixture_index in 0..fixtures {
        let material = if let (Some(corpus), Some(batch_seed)) = (&shared_corpus, batch_seed) {
            format!(
                "KIC-SHARED-PUBLIC-FIXTURE-v1|{n}|{a}|{corpus}|{batch_seed}|{}",
                shared_fixture_offset + fixture_index
            )
        } else if let Some(batch_seed) = batch_seed {
            format!(
                "TASK-KIC-DIRECT-BATCH-20260910|rho|{n}|{a}|{}|{batch_seed}|{fixture_index}",
                mode.name()
            )
        } else {
            format!("{TASK_ID}|rho|{n}|{a}|{}|{fixture_index}", mode.name())
        };
        let digest = blake3::hash(material.as_bytes());
        let seed = u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap());
        let result = if backend == "strong" {
            solve_fixture_strong(
                &curve,
                mode,
                fixture_index,
                seed,
                fixture_target,
                batch_seed,
            )
        } else if backend == "packed" {
            solve_fixture_packed(&curve, mode, fixture_index, seed, fixture_target)
        } else {
            solve_fixture(&curve, mode, fixture_index, seed, fixture_target)
        };
        println!("{}", result);
    }
}

#[cfg(test)]
mod packed_tests {
    use super::*;

    /// The strong backend recovers seeded and explicit scalars and reports the
    /// same JSON fields the stage verifiers read.
    #[test]
    fn strong_backend_recovers_fixture_scalars() {
        let mut ran = 0;
        for (n, a) in [(13, 0), (17, 1), (19, 1), (23, 1), (37, 0), (41, 0)] {
            let Some(curve) = KoblitzCurve::new(a, n) else {
                continue;
            };
            ran += 1;
            let modulus = curve.subgroup_order.to_u64_digits()[0];
            for (index, target) in [
                FixtureTarget::SeededScalar,
                FixtureTarget::ExplicitScalar(modulus / 5 + 3),
            ]
            .into_iter()
            .enumerate()
            {
                let row = solve_fixture_strong(
                    &curve,
                    Quotient::SignedFrobenius,
                    index as u64,
                    1_000 + n as u64,
                    target,
                    Some(531_310),
                );
                assert_eq!(row["verified"], true);
                assert_eq!(row["reference_grade"], "strong");
                assert_eq!(
                    row["recovered_fixture_scalar"], row["published_fixture_scalar"],
                    "n={n} a={a}"
                );
                for key in [
                    "published_q",
                    "walk_steps",
                    "restarts",
                    "setup_ms",
                    "table_entries",
                ] {
                    assert!(!row[key].is_null(), "missing {key}");
                }
            }
        }
        assert!(ran >= 4, "strong backend exercised only {ran} rungs");
    }

    /// Ledger §23's sizes are rungs, and a public target at one the earlier
    /// rungs missed is the point `ic workflow` hashed from the same seed:
    /// §23's `T01` at `icv1-f2m47-t22705043-f4e44623`.
    #[test]
    fn the_rule_comparisons_sizes_are_rungs_with_its_targets() {
        for n in [41, 47, 53, 57, 59, 61] {
            assert!(RUNGS.contains(&n), "n = {n}");
        }
        let curve = KoblitzCurve::new(1, 47).expect("§23's curve at n = 47 constructs");
        let (target, counter) = public_hash_target(&curve, 23_001);
        assert_eq!(
            point_key(&target),
            (1, 59_876_159_786_985, 136_992_923_299_992)
        );
        assert_eq!(counter, 0);
    }

    #[test]
    fn packed_rho_arithmetic_matches_reference() {
        for (n, a) in [
            (7, 1),
            (11, 1),
            (13, 0),
            (17, 1),
            (19, 1),
            (23, 1),
            (37, 0),
            (41, 0),
            (59, 1),
        ] {
            let curve = KoblitzCurve::new(a, n).unwrap();
            let generator = raw_point(curve.generator());
            for scalar in 0..128u64 {
                assert_eq!(
                    raw_scalar_mul(&curve, generator, scalar),
                    raw_point(&curve.mul(curve.generator(), &BigUint::from(scalar)))
                );
            }
            let left = raw_scalar_mul(&curve, generator, 37);
            let right = raw_scalar_mul(&curve, generator, 91);
            assert_eq!(
                raw_add(&curve, left, right),
                raw_point(&curve.add(
                    &curve.mul(curve.generator(), &BigUint::from(37u64)),
                    &curve.mul(curve.generator(), &BigUint::from(91u64)),
                ))
            );
            let RawPoint::Affine { x, .. } = left else {
                panic!("nonzero scalar must be affine");
            };
            let reference = crypto_lib::binary_ecc::F2mElement::from_biguint(&BigUint::from(x), n);
            assert_eq!(
                raw_square(&curve, x),
                reference
                    .square(&curve.curve.irreducible)
                    .raw_bits()
                    .first()
                    .copied()
                    .unwrap_or(0)
            );
        }
    }
}

#[cfg(test)]
mod wide_packed_tests {
    use super::wide::*;
    use super::*;

    fn backend71() -> (KoblitzCurve, WideBackend) {
        let curve = KoblitzCurve::new(0, 71).expect("K_0/F_2^71 must construct");
        let backend = WideBackend::new(&curve).expect("wide backend at n = 71");
        (curve, backend)
    }

    #[test]
    fn wide_raw_arithmetic_matches_reference() {
        let (curve, backend) = backend71();
        let generator = raw_point128(&backend, curve.generator());
        for scalar in 0..128u64 {
            assert_eq!(
                raw_scalar_mul128(&backend, generator, scalar),
                raw_point128(
                    &backend,
                    &curve.mul(curve.generator(), &BigUint::from(scalar))
                ),
                "scalar {scalar}"
            );
        }
        let left = raw_scalar_mul128(&backend, generator, 37);
        let right = raw_scalar_mul128(&backend, generator, 91);
        assert_eq!(
            raw_add128(&backend, left, right),
            raw_point128(
                &backend,
                &curve.add(
                    &curve.mul(curve.generator(), &BigUint::from(37u64)),
                    &curve.mul(curve.generator(), &BigUint::from(91u64)),
                )
            ),
            "add"
        );
        assert_eq!(
            raw_neg128(left),
            raw_point128(
                &backend,
                &point_neg(&curve.mul(curve.generator(), &BigUint::from(37u64)))
            ),
            "neg"
        );
        // Square agrees with the general implementation on the x-word.
        if let Some((x, _)) = left {
            let reference = crypto_lib::binary_ecc::F2mElement::from_biguint(&BigUint::from(x), 71);
            let back = backend
                .gf
                .from_element(&reference.square(&curve.curve.irreducible));
            assert_eq!(backend.fast.gf.sqr(x), back, "square");
        }
    }

    #[test]
    fn wide_partition_is_deterministic_in_range() {
        let (curve, backend) = backend71();
        let generator = raw_point128(&backend, curve.generator());
        let p = raw_scalar_mul128(&backend, generator, 12345);
        let first = raw_partition128(p);
        assert!(first < JUMPS);
        for scalar in [1u64, 2, 999983, 1 << 40] {
            let q = raw_scalar_mul128(&backend, generator, scalar);
            assert_eq!(raw_partition128(q), raw_partition128(q));
            assert!(raw_partition128(q) < JUMPS);
        }
        assert_eq!(raw_partition128(p), first);
    }

    #[test]
    fn wide_canonicalize_picks_minimum_signed_orbit_key() {
        let (curve, backend) = backend71();
        let modulus = curve.subgroup_order.to_u64_digits()[0];
        let lambda = curve.lambda.to_u64_digits()[0];
        let generator = raw_point128(&backend, curve.generator());
        let mut charges = Charges::default();
        for scalar in [7u64, 123456789, 1 << 50] {
            let point = raw_scalar_mul128(&backend, generator, scalar);
            let state = RawState128 { point, a: 1, b: 2 };
            let out = raw_canonicalize128(
                &backend,
                &curve,
                state,
                Quotient::SignedFrobenius,
                modulus,
                lambda,
                &mut charges,
            );
            // The canonical point is in the signed Frobenius orbit:
            // some (k, s) with out = s·λ^k·P reproduces the multiplier
            // relation on the (a, b) side.
            let mut found = false;
            let mut probe = point;
            let mut mult = 1u64;
            for _ in 0..curve.n {
                for (cand, m) in [(probe, mult), (raw_neg128(probe), modulus - mult)] {
                    if raw_key128(cand) == raw_key128(out.point)
                        && mul_mod(state.a, m, modulus) == out.a
                        && mul_mod(state.b, m, modulus) == out.b
                    {
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
                probe = match probe {
                    None => None,
                    Some((x, y)) => Some((backend.gf.sqr(x), backend.gf.sqr(y))),
                };
                mult = mul_mod(mult, lambda, modulus);
            }
            assert!(found, "canonical point must be a signed orbit image");
            // ... and its key is minimal over the whole signed orbit.
            let best = {
                let mut best = raw_key128(point);
                let mut probe = point;
                for _ in 0..curve.n {
                    best = best.min(raw_key128(probe));
                    best = best.min(raw_key128(raw_neg128(probe)));
                    probe = match probe {
                        None => None,
                        Some((x, y)) => Some((backend.gf.sqr(x), backend.gf.sqr(y))),
                    };
                }
                best
            };
            assert_eq!(raw_key128(out.point), best, "canonical key must be minimal");
        }
    }
}
