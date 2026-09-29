#![recursion_limit = "256"]
#![allow(dead_code)]
//! Kuhn–Struik batched Pollard rho over the Koblitz public-fixture corpus —
//! **reference-strengthening ladder** (`KIC_RHO_RUNG`), for the protocol in
//! `research/notes/index-calculus/RESEARCH_STRONG_RHO_LADDER_20260929.md`.
//!
//! Same jump table, distinguished-point table, target derivation, walk
//! schedule, cycle rules and JSON schema as `koblitz_rho_batch_ks.rs` and
//! `koblitz_rho_batch_ks_matched_arith.rs`, restricted to the mode the
//! comparison uses (`signed_frobenius`, derived known-scalar targets, no
//! Bernstein–Lange precomputation). Cumulative rungs:
//!
//! * `0` — `koblitz_rho_batch_ks_matched_arith.rs` reproduced: normal-basis
//!   rotation and Itoh–Tsujii inversion, but the rho file's own field
//!   multiply/square (software carry-less loop on x86-64).
//! * `1` — every field multiply/square is the library `Gf2::mul`/`Gf2::sqr`,
//!   i.e. the code the IC arm runs (hardware carry-less multiply, table
//!   reduction). Field values are identical, so the walk is bit-identical to
//!   rung 0.
//! * `2` — canonical orbit representative = least rotation of the x normal
//!   coordinates, sign chosen by the y coordinate, converted back to the
//!   polynomial basis once; orbit multiplier from a `λ^k` table; `mul_mod` by a
//!   float-quotient estimate. A different (still class-invariant)
//!   representative, so the walk changes.
//! * `3` — `KIC_RHO_LANES` walks advance in lockstep and their `x₁+x₂`
//!   denominators are inverted together by `Gf2::batch_inv`.
//!
//! Environment: `KIC_RHO_RUNG` (default 3), `KIC_RHO_LANES` (default 32, rung 3
//! only), `KIC_RHO_DP_BITS`, `KIC_RHO_BATCH_CORPUS`.
//! Usage: `<n> <a> signed_frobenius <fixtures> <batch_seed>`.

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

const JUMPS: usize = 32;

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
    partition_hashes: u64,
    table_queries: u64,
    table_inserts: u64,
    failed_collisions: u64,
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

// ---------------------------------------------------------------------------
// Rung-0 field arithmetic: the rho file's own (software carry-less on x86-64).
// ---------------------------------------------------------------------------

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

fn raw_square_soft(curve: &KoblitzCurve, value: u64) -> u64 {
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

/// F2-linear map on n <= 64 bits applied by byte tables. Copied verbatim from
/// `examples/koblitz_orbit_dlp_fast.rs`.
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

/// Copied verbatim from `examples/koblitz_orbit_dlp_fast.rs`.
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

/// Frobenius as an O(1) bit rotation in a normal basis built from `Gf2`
/// conjugates; same construction as `koblitz_orbit_dlp_fast.rs`.
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

/// Everything the group law needs. `rung >= 1` routes every field multiply and
/// square through the library `Gf2` (the code the IC arm runs).
struct Field<'a> {
    curve: &'a KoblitzCurve,
    gf: &'a Gf2,
    nb: &'a NormalBasis,
    rung: u8,
}

impl Field<'_> {
    #[inline(always)]
    fn mul(&self, left: u64, right: u64) -> u64 {
        if self.rung >= 1 {
            self.gf.mul(left, right)
        } else {
            raw_reduce(self.curve, carryless_product(left, right))
        }
    }

    #[inline(always)]
    fn sqr(&self, value: u64) -> u64 {
        if self.rung >= 1 {
            self.gf.sqr(value)
        } else {
            raw_square_soft(self.curve, value)
        }
    }

    /// Itoh–Tsujii with the `a^(2^k)` steps as normal-basis rotations, exactly
    /// as `koblitz_orbit_dlp_fast.rs`'s `invert()`.
    fn inv(&self, value: u64) -> u64 {
        assert_ne!(value, 0);
        let e = self.curve.n - 1;
        let mut c = value;
        let mut k = 1u32;
        let bits = 32 - e.leading_zeros();
        for i in (0..bits - 1).rev() {
            let raised = self
                .nb
                .to_poly
                .apply(self.nb.rotate(self.nb.to_normal.apply(c), k));
            c = self.mul(c, raised);
            k *= 2;
            if (e >> i) & 1 == 1 {
                c = self.mul(self.sqr(c), value);
                k += 1;
            }
        }
        self.sqr(c)
    }
}

fn point_neg(point: RawPoint) -> RawPoint {
    match point {
        RawPoint::Infinity => RawPoint::Infinity,
        RawPoint::Affine { x, y } => RawPoint::Affine { x, y: y ^ x },
    }
}

fn point_double(f: &Field, point: RawPoint) -> RawPoint {
    let RawPoint::Affine { x, y } = point else {
        return RawPoint::Infinity;
    };
    if x == 0 {
        return RawPoint::Infinity;
    }
    let lambda = x ^ f.mul(y, f.inv(x));
    let x3 = f.sqr(lambda) ^ lambda ^ f.curve.a as u64;
    let y3 = f.sqr(x) ^ f.mul(lambda ^ 1, x3);
    RawPoint::Affine { x: x3, y: y3 }
}

fn point_add(f: &Field, left: RawPoint, right: RawPoint) -> RawPoint {
    match (left, right) {
        (RawPoint::Infinity, point) | (point, RawPoint::Infinity) => point,
        (RawPoint::Affine { x: x1, y: y1 }, RawPoint::Affine { x: x2, y: y2 }) => {
            if x1 == x2 {
                return if y1 ^ y2 == x1 {
                    RawPoint::Infinity
                } else {
                    point_double(f, left)
                };
            }
            let lambda = f.mul(y1 ^ y2, f.inv(x1 ^ x2));
            add_with_slope(f, x1, y1, x2, lambda)
        }
    }
}

/// Finish an addition given the slope `lambda`.
#[inline(always)]
fn add_with_slope(f: &Field, x1: u64, y1: u64, x2: u64, lambda: u64) -> RawPoint {
    let x3 = f.sqr(lambda) ^ lambda ^ x1 ^ x2 ^ f.curve.a as u64;
    let y3 = f.mul(lambda, x1 ^ x3) ^ x3 ^ y1;
    RawPoint::Affine { x: x3, y: y3 }
}

fn scalar_mul(f: &Field, point: RawPoint, scalar: u64) -> RawPoint {
    let mut result = RawPoint::Infinity;
    for bit in (0..64 - scalar.leading_zeros()).rev() {
        result = point_double(f, result);
        if (scalar >> bit) & 1 == 1 {
            result = point_add(f, result, point);
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Scalar arithmetic modulo the subgroup order.
// ---------------------------------------------------------------------------

fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}

/// `a·b mod m` for `a, b < m < 2^50` by a floating-point quotient estimate:
/// the estimate is off by at most one, so the wrapped remainder lies in
/// `(−m, 2m)` and one conditional correction finishes it.
#[derive(Clone, Copy)]
struct FastMod {
    m: u64,
    inv: f64,
}

impl FastMod {
    fn new(m: u64) -> Self {
        assert!(m > 1 && m < (1u64 << 50), "FastMod needs 1 < m < 2^50");
        Self {
            m,
            inv: 1.0 / m as f64,
        }
    }

    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.m && b < self.m);
        let q = ((a as f64) * (b as f64) * self.inv) as u64;
        let r = a.wrapping_mul(b).wrapping_sub(q.wrapping_mul(self.m)) as i64;
        let m = self.m as i64;
        (if r < 0 {
            r + m
        } else if r >= m {
            r - m
        } else {
            r
        }) as u64
    }
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

// ---------------------------------------------------------------------------
// Canonicalization.
// ---------------------------------------------------------------------------

/// Rungs 0 and 1: the representative and multiplier of PR #955's
/// `raw_canonicalize` (least `(1, x, y)` in the polynomial basis over the whole
/// signed Frobenius orbit), the orbit read off the normal basis.
fn canonicalize_poly(
    f: &Field,
    state: RawState,
    modulus: u64,
    lambda: u64,
    charges: &mut Charges,
) -> RawState {
    charges.canonicalizations += 1;
    let RawPoint::Affine { x: x0, y: y0 } = state.point else {
        return state;
    };
    let nb = f.nb;
    let x_normal = nb.to_normal.apply(x0);
    let y_normal = nb.to_normal.apply(y0);
    let mut multiplier = 1u64;
    let mut best_key = raw_key(state.point);
    let mut best_point = state.point;
    let mut best_multiplier = multiplier;
    for exponent in 0..f.curve.n {
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
        let negative = point_neg(point);
        if raw_key(negative) < best_key {
            best_key = raw_key(negative);
            best_point = negative;
            best_multiplier = modulus - multiplier;
        }
        if exponent + 1 < f.curve.n {
            multiplier = mul_mod(multiplier, lambda, modulus);
        }
    }
    RawState {
        point: best_point,
        a: mul_mod(state.a, best_multiplier, modulus),
        b: mul_mod(state.b, best_multiplier, modulus),
    }
}

/// Reference definition of the rung-2 representative, by exhaustive search over
/// the `2n` orbit members ordered by `(x normal coords, y normal coords)`.
/// Returns `(rotation k, negated)`. Used to define ties and to test the fast
/// path against.
fn normal_coordinate_choice(nb: &NormalBasis, xc: u64, yc: u64) -> (u32, bool) {
    let mut best = (u64::MAX, u64::MAX);
    let mut choice = (0u32, false);
    for k in 0..nb.n {
        let xk = nb.rotate(xc, k);
        let yk = nb.rotate(yc, k);
        if (xk, yk) < best {
            best = (xk, yk);
            choice = (k, false);
        }
        if (xk, yk ^ xk) < best {
            best = (xk, yk ^ xk);
            choice = (k, true);
        }
    }
    choice
}

/// Rung 2: least rotation of the x normal coordinates by integer
/// rotate/compare, sign by the y coordinate, one conversion back.
struct FastCanon<'a> {
    nb: &'a NormalBasis,
    lam_pow: Vec<u64>,
    fm: FastMod,
    modulus: u64,
}

impl<'a> FastCanon<'a> {
    fn new(nb: &'a NormalBasis, lambda: u64, modulus: u64) -> Self {
        let mut lam_pow = Vec::with_capacity(nb.n as usize);
        let mut cur = 1u64;
        for _ in 0..nb.n {
            lam_pow.push(cur);
            cur = mul_mod(cur, lambda, modulus);
        }
        FastCanon {
            nb,
            lam_pow,
            fm: FastMod::new(modulus),
            modulus,
        }
    }

    #[inline]
    fn apply(&self, state: RawState, charges: &mut Charges) -> RawState {
        charges.canonicalizations += 1;
        let RawPoint::Affine { x, y } = state.point else {
            return state;
        };
        let nb = self.nb;
        let n = nb.n;
        let xc = nb.to_normal.apply(x);
        let yc = nb.to_normal.apply(y);
        let mut best = xc;
        let mut best_k = 0u32;
        let mut cur = xc;
        let mut tie = false;
        for k in 1..n {
            cur = ((cur << 1) | (cur >> (n - 1))) & nb.mask;
            if cur < best {
                best = cur;
                best_k = k;
                tie = false;
            } else if cur == best {
                tie = true;
            }
        }
        let (k, negated) = if tie {
            normal_coordinate_choice(nb, xc, yc)
        } else {
            let yk = nb.rotate(yc, best_k);
            (best_k, (yk ^ best) < yk)
        };
        let xk = nb.rotate(xc, k);
        let yk = nb.rotate(yc, k);
        let yk = if negated { yk ^ xk } else { yk };
        let lam = self.lam_pow[k as usize];
        let multiplier = if negated { self.modulus - lam } else { lam };
        RawState {
            point: RawPoint::Affine {
                x: nb.to_poly.apply(xk),
                y: nb.to_poly.apply(yk),
            },
            a: self.fm.mul(state.a, multiplier),
            b: self.fm.mul(state.b, multiplier),
        }
    }
}

enum Canon<'a> {
    Poly {
        f: &'a Field<'a>,
        modulus: u64,
        lambda: u64,
    },
    Fast(FastCanon<'a>),
}

impl Canon<'_> {
    #[inline]
    fn apply(&self, state: RawState, charges: &mut Charges) -> RawState {
        match self {
            Canon::Poly { f, modulus, lambda } => {
                canonicalize_poly(f, state, *modulus, *lambda, charges)
            }
            Canon::Fast(fast) => fast.apply(state, charges),
        }
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

struct RunCfg {
    n: u32,
    a: u8,
    fixtures: u32,
    batch_seed: u64,
    dp_bits: u32,
    corpus: Option<String>,
    rung: u8,
    lanes: usize,
}

struct RunResult {
    total_steps: u64,
    per_fixture_steps: Vec<u64>,
    recovered: Vec<u64>,
    planted: Vec<u64>,
    table_entries: usize,
    cross_solves: u32,
}

struct Lane {
    state: RawState,
    previous: [RawPoint; 4],
    length: u64,
    live: bool,
}

/// Resolve a distinguished point against the shared table. Returns
/// `Some((scalar, via_target))` when it yields a verified logarithm of `q`.
#[allow(clippy::too_many_arguments)]
fn resolve_distinguished_point(
    f: &Field,
    table: &mut HashMap<(u8, u64, u64), Trail>,
    solved: &[u64],
    state: RawState,
    index: u32,
    q: RawPoint,
    generator: RawPoint,
    modulus: u64,
    charges: &mut Charges,
    wasted_merges: &mut u64,
) -> Option<(u64, u32)> {
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
        return None;
    };
    // Both trails reach this point: a + b·d = a' + b'·d' with d' known or = d.
    let candidate = if hit.target == index {
        let denominator = sub_mod(state.b, hit.b, modulus);
        inverse_mod(denominator, modulus)
            .map(|inv| mul_mod(sub_mod(hit.a, state.a, modulus), inv, modulus))
    } else {
        let base = solved[hit.target as usize];
        let known = (hit.a as u128 + mul_mod(hit.b, base, modulus) as u128) % modulus as u128;
        inverse_mod(state.b, modulus)
            .map(|inv| mul_mod(sub_mod(known as u64, state.a, modulus), inv, modulus))
    };
    let Some(candidate) = candidate else {
        *wasted_merges += 1;
        return None;
    };
    charges.scalar_multiplications += 1;
    if scalar_mul(f, generator, candidate) == q {
        return Some((candidate, hit.target));
    }
    charges.failed_collisions += 1;
    None
}

fn run(cfg: &RunCfg, emit: &mut dyn FnMut(Value)) -> RunResult {
    let (n, a) = (cfg.n, cfg.a);
    let process_started = Instant::now();
    let curve = KoblitzCurve::new(a, n).expect("frozen exact rung must construct");
    let gf = Gf2::new(&curve.curve.irreducible);
    let nb = NormalBasis::new(&gf);
    let f = Field {
        curve: &curve,
        gf: &gf,
        nb: &nb,
        rung: cfg.rung,
    };
    let modulus = curve.subgroup_order.to_u64_digits()[0];
    let lambda = curve.lambda.to_u64_digits()[0];
    let automorphisms = signed_automorphism_size(lambda, modulus, n);
    let generator = raw_point(curve.generator());
    let canon = if cfg.rung >= 2 {
        Canon::Fast(FastCanon::new(&nb, lambda, modulus))
    } else {
        Canon::Poly {
            f: &f,
            modulus,
            lambda,
        }
    };
    let mut charges = Charges::default();

    let setup_started = Instant::now();
    let jump_digest =
        blake3::hash(format!("KIC-KS-BATCH-JUMPS-v1|{n}|{a}|{}", cfg.batch_seed).as_bytes());
    let mut jump_rng = StdRng::seed_from_u64(u64::from_le_bytes(
        jump_digest.as_bytes()[..8].try_into().unwrap(),
    ));
    let jumps: Vec<RawJump> = (0..JUMPS)
        .map(|_| loop {
            let s = jump_rng.gen_range(1..modulus);
            charges.scalar_multiplications += 1;
            let point = scalar_mul(&f, generator, s);
            if point != RawPoint::Infinity {
                break RawJump { point, a: s, b: 0 };
            }
        })
        .collect();
    let setup_ms = setup_started.elapsed().as_secs_f64() * 1000.0;

    let dp_mask = (1u64 << cfg.dp_bits) - 1;
    let walk_cap = 8u64 << cfg.dp_bits;
    let ideal_single =
        (std::f64::consts::PI * modulus as f64 / (2.0 * automorphisms as f64)).sqrt();
    let step_cap = (ideal_single.ceil() as u64)
        .saturating_mul(2_000)
        .max(1_000_000);
    let mut table: HashMap<(u8, u64, u64), Trail> = HashMap::new();
    let mut solved: Vec<u64> = Vec::with_capacity(cfg.fixtures as usize);
    let mut planted: Vec<u64> = Vec::with_capacity(cfg.fixtures as usize);
    let mut per_fixture_steps: Vec<u64> = Vec::with_capacity(cfg.fixtures as usize);
    let mut total_steps = 0u64;
    let mut cross_solves = 0u32;

    let mut denominators: Vec<u64> = vec![0; cfg.lanes.max(1)];
    let mut jump_index: Vec<usize> = vec![0; cfg.lanes.max(1)];
    let mut active: Vec<usize> = Vec::with_capacity(cfg.lanes.max(1));
    let mut scratch: Vec<u64> = Vec::new();

    for index in 0..cfg.fixtures {
        let material = match &cfg.corpus {
            Some(corpus) => format!(
                "KIC-SHARED-PUBLIC-FIXTURE-v1|{n}|{a}|{corpus}|{}|{index}",
                cfg.batch_seed
            ),
            None => format!(
                "TASK-KIC-DIRECT-BATCH-20260910|rho|{n}|{a}|signed_frobenius|{}|{index}",
                cfg.batch_seed
            ),
        };
        let digest = blake3::hash(material.as_bytes());
        let seed = u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap());
        let mut rng = StdRng::seed_from_u64(seed);
        let generated_d0 = rng.gen_range(1..modulus);
        let started = Instant::now();
        charges.scalar_multiplications += 1;
        let q = scalar_mul(&f, generator, generated_d0);
        let table_before = table.len();
        let stride_a = rng.gen_range(1..modulus);
        let stride = scalar_mul(&f, generator, stride_a);
        let mut cursor_a = rng.gen_range(0..modulus);
        let mut cursor = point_add(&f, scalar_mul(&f, generator, cursor_a), q);
        charges.scalar_multiplications += 2;
        charges.group_additions += 1;

        let mut steps = 0u64;
        let mut walks = 0u64;
        let mut fruitless = 0u64;
        let mut capped = 0u64;
        let mut wasted_merges = 0u64;
        let mut recovered = None;
        let mut via_target = None;

        if cfg.rung < 3 {
            'walks: while steps < step_cap {
                walks += 1;
                let start = cursor;
                let start_a = cursor_a;
                cursor = point_add(&f, cursor, stride);
                cursor_a = (cursor_a + stride_a) % modulus;
                charges.group_additions += 1;
                if start == RawPoint::Infinity {
                    continue;
                }
                let mut state = canon.apply(
                    RawState {
                        point: start,
                        a: start_a,
                        b: 1,
                    },
                    &mut charges,
                );
                let mut previous = [RawPoint::Infinity; 4];
                let mut length = 0u64;
                while dp_hash(state.point) & dp_mask != 0 {
                    let jump = &jumps[raw_partition(state.point)];
                    charges.partition_hashes += 1;
                    let next = RawState {
                        point: point_add(&f, state.point, jump.point),
                        a: (state.a + jump.a) % modulus,
                        b: state.b,
                    };
                    charges.group_additions += 1;
                    let next = canon.apply(next, &mut charges);
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
                if let Some((scalar, via)) = resolve_distinguished_point(
                    &f,
                    &mut table,
                    &solved,
                    state,
                    index,
                    q,
                    generator,
                    modulus,
                    &mut charges,
                    &mut wasted_merges,
                ) {
                    recovered = Some(scalar);
                    via_target = Some(via);
                    break;
                }
            }
        } else {
            let lanes_n = cfg.lanes.max(1);
            let mut lanes: Vec<Lane> = (0..lanes_n)
                .map(|_| Lane {
                    state: RawState {
                        point: RawPoint::Infinity,
                        a: 0,
                        b: 0,
                    },
                    previous: [RawPoint::Infinity; 4],
                    length: 0,
                    live: false,
                })
                .collect();
            'target: while steps < step_cap {
                // Phase A: start a walk in every idle lane.
                for lane in lanes.iter_mut().filter(|lane| !lane.live) {
                    loop {
                        walks += 1;
                        let start = cursor;
                        let start_a = cursor_a;
                        cursor = point_add(&f, cursor, stride);
                        cursor_a = (cursor_a + stride_a) % modulus;
                        charges.group_additions += 1;
                        if start == RawPoint::Infinity {
                            continue;
                        }
                        lane.state = canon.apply(
                            RawState {
                                point: start,
                                a: start_a,
                                b: 1,
                            },
                            &mut charges,
                        );
                        lane.previous = [RawPoint::Infinity; 4];
                        lane.length = 0;
                        lane.live = true;
                        break;
                    }
                }
                // Phase B: lanes standing on a distinguished point.
                for lane in lanes.iter_mut() {
                    if lane.live && dp_hash(lane.state.point) & dp_mask == 0 {
                        lane.live = false;
                        if let Some((scalar, via)) = resolve_distinguished_point(
                            &f,
                            &mut table,
                            &solved,
                            lane.state,
                            index,
                            q,
                            generator,
                            modulus,
                            &mut charges,
                            &mut wasted_merges,
                        ) {
                            recovered = Some(scalar);
                            via_target = Some(via);
                            break 'target;
                        }
                    }
                }
                // Phase C: one step for every live lane, denominators inverted
                // together.
                active.clear();
                for (li, lane) in lanes.iter().enumerate() {
                    if !lane.live {
                        continue;
                    }
                    let ji = raw_partition(lane.state.point);
                    charges.partition_hashes += 1;
                    let j = active.len();
                    jump_index[j] = ji;
                    denominators[j] = match (lane.state.point, jumps[ji].point) {
                        (RawPoint::Affine { x: x1, .. }, RawPoint::Affine { x: x2, .. }) => x1 ^ x2,
                        _ => 0,
                    };
                    active.push(li);
                }
                f.gf.batch_inv(&mut denominators[..active.len()], &mut scratch);
                for (j, &li) in active.iter().enumerate() {
                    let lane = &mut lanes[li];
                    let jump = &jumps[jump_index[j]];
                    let sum = match (lane.state.point, jump.point) {
                        (RawPoint::Affine { x: x1, y: y1 }, RawPoint::Affine { x: x2, y: y2 })
                            if denominators[j] != 0 =>
                        {
                            let slope = f.mul(y1 ^ y2, denominators[j]);
                            add_with_slope(&f, x1, y1, x2, slope)
                        }
                        _ => point_add(&f, lane.state.point, jump.point),
                    };
                    let next = RawState {
                        point: sum,
                        a: (lane.state.a + jump.a) % modulus,
                        b: lane.state.b,
                    };
                    charges.group_additions += 1;
                    let next = canon.apply(next, &mut charges);
                    steps += 1;
                    lane.length += 1;
                    if next.point == lane.state.point || lane.previous.contains(&next.point) {
                        fruitless += 1;
                        lane.live = false;
                        continue;
                    }
                    if lane.length > walk_cap {
                        capped += 1;
                        lane.live = false;
                        continue;
                    }
                    lane.previous = [
                        lane.state.point,
                        lane.previous[0],
                        lane.previous[1],
                        lane.previous[2],
                    ];
                    lane.state = next;
                }
            }
        }

        let solve_ms = started.elapsed().as_secs_f64() * 1000.0;
        let recovered = recovered.expect("batched rho exceeded the per-target step cap");
        assert_eq!(recovered, generated_d0, "recovered scalar is wrong");
        let via = via_target.unwrap();
        if via != index {
            cross_solves += 1;
        }
        solved.push(recovered);
        planted.push(generated_d0);
        per_fixture_steps.push(steps);
        total_steps += steps;
        let q_key = raw_key(q);
        emit(json!({
            "kind":"rho_ks_batch_fixture",
            "evidence_class":"measured_rho_observation",
            "n":n,"a":a,"quotient_mode":"signed_frobenius","automorphism_size":automorphisms,
            "fixture_index":index,"fixture_seed":seed,"batch_seed":cfg.batch_seed,
            "target_source":"derived_known_scalar",
            "published_fixture_scalar":generated_d0,
            "recovered_fixture_scalar":recovered,
            "published_q":[q_key.1,q_key.2],"verified":true,
            "solved_via_target":via,"cross_target_solve":via != index,
            "walk_steps":steps,"walks":walks,"fruitless_two_cycles":fruitless,
            "capped_walks":capped,"wasted_merges":wasted_merges,
            "table_entries_before":table_before,"table_entries_after":table.len(),
            "total_ms":solve_ms,
            "ideal_independent_steps":ideal_single,
        }));
    }
    let process_ms = process_started.elapsed().as_secs_f64() * 1000.0;
    emit(json!({
        "kind":"rho_ks_batch_summary","producer_version":"v5_strong_rho_ladder",
        "rung":cfg.rung,"lanes":if cfg.rung >= 3 { cfg.lanes } else { 1 },
        "n":n,"a":a,"quotient_mode":"signed_frobenius","automorphism_size":automorphisms,
        "fixtures":cfg.fixtures,"batch_seed":cfg.batch_seed,"corpus":cfg.corpus,
        "dp_bits":cfg.dp_bits,"jump_count":JUMPS,
        "target_source":"derived_known_scalar",
        "library_gf2_kernel":gf.kernel_name(),
        "rho_field_path":if cfg.rung >= 1 { "library Gf2 (same code as IC)" } else { "rho-local carry-less product and bit-serial reduction" },
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
            "failed_collisions":charges.failed_collisions
        },
        "scope":"published synthetic toy fixtures; no external point, unknown scalar, or production key",
    }));
    RunResult {
        total_steps,
        per_fixture_steps,
        recovered: solved,
        planted,
        table_entries: table.len(),
        cross_solves,
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        6,
        "usage: <n> <a> signed_frobenius <fixtures> <batch_seed>"
    );
    assert_eq!(
        args[3], "signed_frobenius",
        "only signed_frobenius is built"
    );
    let cfg = RunCfg {
        n: args[1].parse().unwrap(),
        a: args[2].parse().unwrap(),
        fixtures: args[4].parse().unwrap(),
        batch_seed: args[5].parse().unwrap(),
        dp_bits: std::env::var("KIC_RHO_DP_BITS")
            .map(|v| v.parse().expect("KIC_RHO_DP_BITS must be an integer"))
            .unwrap_or(8),
        corpus: std::env::var("KIC_RHO_BATCH_CORPUS").ok(),
        rung: std::env::var("KIC_RHO_RUNG")
            .map(|v| v.parse().expect("KIC_RHO_RUNG must be 0..=3"))
            .unwrap_or(3),
        lanes: std::env::var("KIC_RHO_LANES")
            .map(|v| v.parse().expect("KIC_RHO_LANES must be an integer"))
            .unwrap_or(32),
    };
    assert!(cfg.dp_bits < 32 && cfg.fixtures > 0 && cfg.rung <= 3);
    assert!(matches!(
        cfg.n,
        7 | 11 | 13 | 17 | 19 | 23 | 37 | 41 | 53 | 61
    ));
    run(&cfg, &mut |line| println!("{line}"));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng_from(seed: u64) -> StdRng {
        StdRng::seed_from_u64(seed)
    }

    const SIZES: [u32; 9] = [7, 11, 13, 17, 19, 23, 37, 41, 53];

    /// Every `(a, n)` in `SIZES × {0, 1}` that constructs a frozen exact rung.
    fn each_curve(mut body: impl FnMut(u8, u32, KoblitzCurve)) {
        let mut seen = 0;
        for &n in &SIZES {
            for a in [0u8, 1u8] {
                if let Some(curve) = KoblitzCurve::new(a, n) {
                    body(a, n, curve);
                    seen += 1;
                }
            }
        }
        assert!(seen >= 4, "too few constructible curves: {seen}");
    }

    /// R1 gate: the library `Gf2` and this file's own rung-0 arithmetic are the
    /// same field, so switching rungs cannot change a value.
    #[test]
    fn gf2_matches_rung0_field_arithmetic() {
        for &n in &SIZES {
            for a in [0u8, 1u8] {
                let Some(curve) = KoblitzCurve::new(a, n) else {
                    continue;
                };
                let gf = Gf2::new(&curve.curve.irreducible);
                let mask = (1u64 << n) - 1;
                let mut rng = rng_from(0x6b6f_626c_6974_7a00 ^ (n as u64) ^ ((a as u64) << 40));
                for _ in 0..2000 {
                    let x: u64 = rng.gen::<u64>() & mask;
                    let y: u64 = rng.gen::<u64>() & mask;
                    assert_eq!(raw_square_soft(&curve, x), gf.sqr(x), "sqr n={n} a={a}");
                    assert_eq!(
                        raw_reduce(&curve, carryless_product(x, y)),
                        gf.mul(x, y),
                        "mul n={n} a={a}"
                    );
                }
            }
        }
    }

    /// The Itoh–Tsujii inverse is a genuine inverse under every rung's field.
    #[test]
    fn inverse_is_correct_under_every_rung() {
        each_curve(|a, n, curve| {
            let gf = Gf2::new(&curve.curve.irreducible);
            let nb = NormalBasis::new(&gf);
            let mask = (1u64 << n) - 1;
            let mut rng = rng_from(0x696e_7665_7274_0000 ^ (n as u64) ^ ((a as u64) << 40));
            for rung in 0..=3u8 {
                let f = Field {
                    curve: &curve,
                    gf: &gf,
                    nb: &nb,
                    rung,
                };
                for _ in 0..200 {
                    let x: u64 = rng.gen::<u64>() & mask;
                    if x == 0 {
                        continue;
                    }
                    assert_eq!(f.mul(x, f.inv(x)), 1, "inverse n={n} a={a} rung={rung}");
                }
            }
        });
    }

    /// `FastMod` equals the `u128 %` reference on random and edge inputs.
    #[test]
    fn fast_mod_matches_u128_reference() {
        let mut rng = rng_from(0x6d6f_6472_6564_0001);
        let moduli = [
            3u64,
            97,
            (1 << 20) + 7,
            (1 << 33) - 9,
            (1 << 44) + 5,
            (1u64 << 49) - 81,
            (1u64 << 50) - 27,
        ];
        for &m in &moduli {
            let fm = FastMod::new(m);
            let edge = [0, 1, 2, m / 2, m - 2, m - 1];
            for &a in &edge {
                for &b in &edge {
                    assert_eq!(fm.mul(a, b), mul_mod(a, b, m), "edge m={m} a={a} b={b}");
                }
            }
            for _ in 0..200_000 {
                let (a, b) = (rng.gen_range(0..m), rng.gen_range(0..m));
                assert_eq!(fm.mul(a, b), mul_mod(a, b, m), "m={m} a={a} b={b}");
            }
        }
    }

    /// Rung-2 canonicalization: a class function (same result for every orbit
    /// member and its negative), the fast path equals the exhaustive
    /// definition, and the returned multiplier really carries the input point
    /// to the representative (`canon = [m]P`, replayed with the group law).
    #[test]
    fn fast_canonicalization_is_class_invariant_and_multiplier_correct() {
        each_curve(|a, n, curve| {
            let gf = Gf2::new(&curve.curve.irreducible);
            let nb = NormalBasis::new(&gf);
            let f = Field {
                curve: &curve,
                gf: &gf,
                nb: &nb,
                rung: 1,
            };
            let modulus = curve.subgroup_order.to_u64_digits()[0];
            let lambda = curve.lambda.to_u64_digits()[0];
            let fast = FastCanon::new(&nb, lambda, modulus);
            let generator = raw_point(curve.generator());
            let mut rng = rng_from(0x6661_7374_6361_6e00 ^ (n as u64) ^ ((a as u64) << 40));
            let mut charges = Charges::default();
            let mut examined = 0;
            while examined < 60 {
                let scalar = rng.gen_range(1..modulus);
                let p = scalar_mul(&f, generator, scalar);
                let RawPoint::Affine { x, y } = p else {
                    continue;
                };
                // Multiplier replay: a = 1, b = 0 exposes the multiplier.
                let probe = fast.apply(
                    RawState {
                        point: p,
                        a: 1,
                        b: 0,
                    },
                    &mut charges,
                );
                assert_eq!(
                    scalar_mul(&f, p, probe.a),
                    probe.point,
                    "n={n}: multiplier does not carry P to the representative"
                );
                // Fast path equals the exhaustive normal-coordinate definition.
                let (xc, yc) = (nb.to_normal.apply(x), nb.to_normal.apply(y));
                let (k, negated) = normal_coordinate_choice(&nb, xc, yc);
                let (xk, yk) = (nb.rotate(xc, k), nb.rotate(yc, k));
                let yk = if negated { yk ^ xk } else { yk };
                assert_eq!(
                    probe.point,
                    RawPoint::Affine {
                        x: nb.to_poly.apply(xk),
                        y: nb.to_poly.apply(yk)
                    },
                    "n={n}: fast path differs from the exhaustive definition"
                );
                // Class invariance over the whole signed orbit.
                let (mut ox, mut oy) = (x, y);
                for _ in 0..n {
                    for candidate in [
                        RawPoint::Affine { x: ox, y: oy },
                        point_neg(RawPoint::Affine { x: ox, y: oy }),
                    ] {
                        let got = fast.apply(
                            RawState {
                                point: candidate,
                                a: 1,
                                b: 0,
                            },
                            &mut charges,
                        );
                        assert_eq!(got.point, probe.point, "n={n}: not a class function");
                    }
                    ox = gf.sqr(ox);
                    oy = gf.sqr(oy);
                }
                examined += 1;
            }
        });
    }

    fn cfg(a: u8, n: u32, fixtures: u32, rung: u8, lanes: usize) -> RunCfg {
        RunCfg {
            n,
            a,
            fixtures,
            batch_seed: 531_310,
            dp_bits: 4,
            corpus: Some(format!("strong-test-a{a}-n{n}")),
            rung,
            lanes,
        }
    }

    /// G2 in miniature: rung 1 changes no field value, so the whole walk is
    /// bit-identical to rung 0 — per-fixture steps, table size, cross solves.
    #[test]
    fn rung1_walk_is_bit_identical_to_rung0() {
        let mut ran = 0;
        each_curve(|a, n, _| {
            if !(13..=23).contains(&n) {
                return;
            }
            ran += 1;
            let r0 = run(&cfg(a, n, 40, 0, 1), &mut |_| {});
            let r1 = run(&cfg(a, n, 40, 1, 1), &mut |_| {});
            assert_eq!(r0.per_fixture_steps, r1.per_fixture_steps, "n={n} a={a}");
            assert_eq!(r0.table_entries, r1.table_entries, "n={n} a={a}");
            assert_eq!(r0.cross_solves, r1.cross_solves, "n={n} a={a}");
            assert_eq!(r0.recovered, r1.recovered, "n={n} a={a}");
        });
        assert!(ran >= 2, "rung-1 identity exercised only {ran} curves");
    }

    /// End to end at small n: every rung recovers every planted scalar (the
    /// run asserts each recovered scalar equals the planted one and verified
    /// `[d]G = Q`).
    #[test]
    fn every_rung_recovers_every_target_end_to_end() {
        let mut ran = 0;
        each_curve(|a, n, _| {
            if !(13..=23).contains(&n) {
                return;
            }
            ran += 1;
            for (rung, lanes) in [(0u8, 1usize), (1, 1), (2, 1), (3, 8), (3, 32)] {
                let r = run(&cfg(a, n, 60, rung, lanes), &mut |_| {});
                assert_eq!(
                    r.recovered, r.planted,
                    "n={n} a={a} rung={rung} lanes={lanes}"
                );
                assert!(r.total_steps > 0);
            }
        });
        assert!(ran >= 2, "end-to-end exercised only {ran} curves");
    }
}
