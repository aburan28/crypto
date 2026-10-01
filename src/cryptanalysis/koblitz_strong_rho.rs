//! Strong single-target Pollard rho on Koblitz curves `E_a / F_{2^n}`: the
//! reference that `vs_rho` comparisons should be measured against.
//!
//! Each property below fixes a weakness of a reference that an earlier
//! comparison in this repository was measured against (see the 2026-09-29 and
//! 2026-09-30 errata in `docs/ic/BOUNDARY_TARGETS.md`):
//!
//! * **distinguished points**: only states whose hash has `dp_bits` low zero
//!   bits are stored, instead of a table holding every step;
//! * **signed-Frobenius quotient in a normal basis**: Frobenius is a rotation of
//!   normal coordinates, so the canonical representative of the orbit
//!   `{±φ^k(P)}` is the least rotation of the x normal coordinates (sign chosen
//!   by y), converted back once; the orbit multiplier `±λ^k` comes from a table;
//! * **library field arithmetic**: every multiply, square and inversion is
//!   [`Gf2`] (hardware carry-less multiply where the CPU has one), inversion by
//!   Itoh–Tsujii with its Frobenius powers taken as rotations;
//! * **lockstep lanes**: `lanes` independent walks advance together and their
//!   slope denominators are inverted with one field inversion
//!   ([`Gf2::batch_inv`], Montgomery's trick).
//!
//! The algorithm is rung 3 of `examples/koblitz_rho_batch_ks_strong.rs` at one
//! target. With the same seeds it walks the same trajectory step for step;
//! `koblitz_rho_fixture … strong` is checked against that example in
//! `research/notes/index-calculus/RESEARCH_FIXED_REFERENCE_RHO_20261001.md`.
//!
//! Walk construction (Kuhn–Struik style, single target): every walk starts at
//! `c·G + Q` (`b = 1`) for a fresh `c`, jumps are known multiples of `G`, and the
//! state `(P, a, b)` keeps `P = a·G + b·Q`; canonicalization multiplies `a` and
//! `b` by the orbit multiplier. Two walks meeting at a stored distinguished
//! point give `a + b·d = a' + b'·d`, solved when `b ≠ b'` and verified by
//! `[d]G = Q` before it is returned.

use crate::binary_ecc::BinaryPoint;
use crate::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crate::cryptanalysis::semaev_decomp::Gf2;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::{BTreeSet, HashMap};

/// Number of precomputed jumps (`r`-adding walk).
pub const JUMPS: usize = 32;

/// A point of `E(F_{2^n})` with coordinates packed in one `u64` each
/// (polynomial basis, as [`Gf2`] uses).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RawPoint {
    /// The point at infinity.
    Infinity,
    /// An affine point.
    Affine {
        /// x-coordinate.
        x: u64,
        /// y-coordinate.
        y: u64,
    },
}

impl RawPoint {
    /// Pack a library point (`n ≤ 64`).
    pub fn from_binary(point: &BinaryPoint) -> Self {
        match point {
            BinaryPoint::Infinity => RawPoint::Infinity,
            BinaryPoint::Affine { x, y } => RawPoint::Affine {
                x: x.raw_bits().first().copied().unwrap_or(0),
                y: y.raw_bits().first().copied().unwrap_or(0),
            },
        }
    }

    fn key(self) -> (u8, u64, u64) {
        match self {
            RawPoint::Infinity => (0, 0, 0),
            RawPoint::Affine { x, y } => (1, x, y),
        }
    }
}

/// Walk state: `point = a·G + b·Q`.
#[derive(Clone, Copy, Debug)]
pub struct WalkState {
    /// Current point.
    pub point: RawPoint,
    /// Coefficient of `G`.
    pub a: u64,
    /// Coefficient of `Q`.
    pub b: u64,
}

/// One precomputed jump `s·G`.
#[derive(Clone, Copy, Debug)]
pub struct Jump {
    /// The point `s·G`.
    pub point: RawPoint,
    /// Its discrete log `s`.
    pub a: u64,
}

/// Tunables. The defaults are the configuration measured as the strong
/// reference in PR #1090 (`lanes = 32`, `dp_bits = 4`).
#[derive(Clone, Copy, Debug)]
pub struct StrongRhoParams {
    /// Walks advanced in lockstep (batch-inversion width); at least 1.
    pub lanes: usize,
    /// A state is distinguished when `dp_hash & (2^dp_bits − 1) == 0`.
    pub dp_bits: u32,
    /// Give up after `step_cap_factor × ⌈√(πr / 2A)⌉` steps (at least 10^6).
    pub step_cap_factor: u64,
}

impl Default for StrongRhoParams {
    fn default() -> Self {
        Self {
            lanes: 32,
            dp_bits: 4,
            step_cap_factor: 2_000,
        }
    }
}

/// Operation counts for one solve (whole solve, setup included).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StrongRhoCharges {
    /// Group additions (walk steps and walk starts).
    pub group_additions: u64,
    /// Scalar multiplications (jumps, start stride, candidate checks).
    pub scalar_multiplications: u64,
    /// Orbit canonicalizations.
    pub canonicalizations: u64,
    /// Jump-index hashes.
    pub partition_hashes: u64,
    /// Distinguished-point table lookups.
    pub table_queries: u64,
    /// Distinguished-point table insertions.
    pub table_inserts: u64,
    /// Collisions whose candidate failed `[d]G = Q`.
    pub failed_collisions: u64,
}

/// Result of a successful solve.
#[derive(Clone, Debug)]
pub struct StrongRhoOutcome {
    /// The verified discrete logarithm `d` with `[d]G = Q`.
    pub scalar: u64,
    /// Group-operation steps taken by all walks.
    pub walk_steps: u64,
    /// Walks started.
    pub walks: u64,
    /// Walks abandoned in a fruitless cycle (length ≤ 4).
    pub fruitless: u64,
    /// Walks abandoned at the length cap `8·2^dp_bits`.
    pub capped: u64,
    /// Collisions with `b = b'` (no information).
    pub wasted_merges: u64,
    /// Distinguished points stored when the solve ended.
    pub table_entries: usize,
    /// `√(πr / 2A)`, the expected step count of one ideal walk.
    pub ideal_steps: f64,
    /// Size `A` of the signed-Frobenius automorphism group on the subgroup.
    pub automorphisms: usize,
    /// Operation counts.
    pub charges: StrongRhoCharges,
}

/// F2-linear map on `n ≤ 64` bits applied by byte tables.
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

/// Invert the `n × n` matrix over F2 whose columns are `columns`; `None` when
/// singular. Returns the columns of the inverse.
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

/// Normal basis built from the conjugates of the first normal element: in it,
/// squaring (Frobenius) is a one-bit rotation.
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
        unreachable!("every finite field has a normal basis")
    }

    #[inline(always)]
    fn rotate(&self, v: u64, k: u32) -> u64 {
        if k == 0 {
            return v;
        }
        ((v << k) | (v >> (self.n - k))) & self.mask
    }
}

/// `a·b mod m`: a floating-point quotient estimate when `m < 2^50` (off by at
/// most one, so one conditional correction finishes it), else `u128`.
#[derive(Clone, Copy)]
struct MulMod {
    m: u64,
    inv: f64,
    fast: bool,
}

impl MulMod {
    fn new(m: u64) -> Self {
        assert!(m > 1, "modulus must exceed 1");
        Self {
            m,
            inv: 1.0 / m as f64,
            fast: m < (1u64 << 50),
        }
    }

    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.m && b < self.m);
        if !self.fast {
            return ((a as u128 * b as u128) % self.m as u128) as u64;
        }
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

fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}

fn add_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 + right as u128) % modulus as u128) as u64
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

/// Jump index of a point (a fixed 64-bit mixer; part of the walk definition).
fn partition(point: RawPoint) -> usize {
    let (_, x, y) = point.key();
    let mut value = x ^ y.rotate_left(21) ^ 0x9e37_79b9_7f4a_7c15;
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^= value >> 31;
    value as usize % JUMPS
}

/// Distinguished-point hash (independent of [`partition`]).
fn dp_hash(point: RawPoint) -> u64 {
    let (_, x, y) = point.key();
    let mut v = x.rotate_left(7) ^ y ^ 0x2545_f491_4f6c_dd1d;
    v ^= v >> 33;
    v = v.wrapping_mul(0xff51_afd7_ed55_8ccd);
    v ^= v >> 33;
    v = v.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    v ^ (v >> 33)
}

#[derive(Clone, Copy)]
struct Lane {
    state: WalkState,
    previous: [RawPoint; 4],
    length: u64,
    live: bool,
}

/// Per-curve precomputation: field, normal basis, orbit-multiplier table.
pub struct StrongRho {
    n: u32,
    curve_a: u64,
    gf: Gf2,
    nb: NormalBasis,
    modulus: u64,
    lambda: u64,
    automorphisms: usize,
    generator: RawPoint,
    lam_pow: Vec<u64>,
    mm: MulMod,
}

impl StrongRho {
    /// Precompute for `curve`. Panics when `n > 63` or the subgroup order does
    /// not fit in a `u64` (the packed representation needs both).
    pub fn new(curve: &KoblitzCurve) -> Self {
        assert!(curve.n <= 63, "packed field elements need n <= 63");
        let digits = curve.subgroup_order.to_u64_digits();
        assert_eq!(digits.len(), 1, "subgroup order must fit in a u64");
        let modulus = digits[0];
        let lambda = curve.lambda.to_u64_digits().first().copied().unwrap_or(0);
        let gf = Gf2::new(&curve.curve.irreducible);
        let nb = NormalBasis::new(&gf);
        let mut lam_pow = Vec::with_capacity(curve.n as usize);
        let mut cur = 1u64;
        let mut orbit = BTreeSet::new();
        for _ in 0..curve.n {
            lam_pow.push(cur);
            orbit.insert(cur);
            orbit.insert((modulus - cur) % modulus);
            cur = mul_mod(cur, lambda, modulus);
        }
        assert_eq!(cur, 1, "lambda must have order dividing n");
        Self {
            n: curve.n,
            curve_a: curve.a as u64,
            gf,
            nb,
            modulus,
            lambda,
            automorphisms: orbit.len(),
            generator: RawPoint::from_binary(curve.generator()),
            lam_pow,
            mm: MulMod::new(modulus),
        }
    }

    /// Prime subgroup order `r`.
    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    /// Eigenvalue `λ` of Frobenius on the subgroup.
    pub fn lambda(&self) -> u64 {
        self.lambda
    }

    /// Size of the signed-Frobenius automorphism group on the subgroup.
    pub fn automorphisms(&self) -> usize {
        self.automorphisms
    }

    /// The generator `G`.
    pub fn generator(&self) -> RawPoint {
        self.generator
    }

    /// Expected steps of one ideal walk, `√(πr / 2A)`.
    pub fn ideal_steps(&self) -> f64 {
        (std::f64::consts::PI * self.modulus as f64 / (2.0 * self.automorphisms as f64)).sqrt()
    }

    #[inline(always)]
    fn mul(&self, left: u64, right: u64) -> u64 {
        self.gf.mul(left, right)
    }

    #[inline(always)]
    fn sqr(&self, value: u64) -> u64 {
        self.gf.sqr(value)
    }

    /// Itoh–Tsujii inversion with the `x^(2^k)` steps as normal-basis rotations.
    fn inv(&self, value: u64) -> u64 {
        assert_ne!(value, 0);
        let e = self.n - 1;
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

    fn double(&self, point: RawPoint) -> RawPoint {
        let RawPoint::Affine { x, y } = point else {
            return RawPoint::Infinity;
        };
        if x == 0 {
            return RawPoint::Infinity;
        }
        let lambda = x ^ self.mul(y, self.inv(x));
        let x3 = self.sqr(lambda) ^ lambda ^ self.curve_a;
        let y3 = self.sqr(x) ^ self.mul(lambda ^ 1, x3);
        RawPoint::Affine { x: x3, y: y3 }
    }

    /// Group law.
    pub fn add(&self, left: RawPoint, right: RawPoint) -> RawPoint {
        match (left, right) {
            (RawPoint::Infinity, point) | (point, RawPoint::Infinity) => point,
            (RawPoint::Affine { x: x1, y: y1 }, RawPoint::Affine { x: x2, y: y2 }) => {
                if x1 == x2 {
                    return if y1 ^ y2 == x1 {
                        RawPoint::Infinity
                    } else {
                        self.double(left)
                    };
                }
                let lambda = self.mul(y1 ^ y2, self.inv(x1 ^ x2));
                self.add_with_slope(x1, y1, x2, lambda)
            }
        }
    }

    #[inline(always)]
    fn add_with_slope(&self, x1: u64, y1: u64, x2: u64, lambda: u64) -> RawPoint {
        let x3 = self.sqr(lambda) ^ lambda ^ x1 ^ x2 ^ self.curve_a;
        let y3 = self.mul(lambda, x1 ^ x3) ^ x3 ^ y1;
        RawPoint::Affine { x: x3, y: y3 }
    }

    /// `k·P` by double-and-add.
    pub fn scalar_mul(&self, point: RawPoint, scalar: u64) -> RawPoint {
        let mut result = RawPoint::Infinity;
        for bit in (0..64 - scalar.leading_zeros()).rev() {
            result = self.double(result);
            if (scalar >> bit) & 1 == 1 {
                result = self.add(result, point);
            }
        }
        result
    }

    /// Exhaustive definition of the canonical representative: the least
    /// `(x, y)` in normal coordinates over the `2n` orbit members. Returns
    /// `(rotation k, negated)`. Used for ties and as the test reference.
    fn choice_exhaustive(&self, xc: u64, yc: u64) -> (u32, bool) {
        let nb = &self.nb;
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

    /// Canonical representative of the signed-Frobenius orbit of the state's
    /// point, with `a, b` multiplied by the orbit multiplier `±λ^k`.
    pub fn canonicalize(&self, state: WalkState) -> WalkState {
        let RawPoint::Affine { x, y } = state.point else {
            return state;
        };
        let nb = &self.nb;
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
            self.choice_exhaustive(xc, yc)
        } else {
            let yk = nb.rotate(yc, best_k);
            (best_k, (yk ^ best) < yk)
        };
        let xk = nb.rotate(xc, k);
        let yk = nb.rotate(yc, k);
        let yk = if negated { yk ^ xk } else { yk };
        let lam = self.lam_pow[k as usize];
        let multiplier = if negated { self.modulus - lam } else { lam };
        WalkState {
            point: RawPoint::Affine {
                x: nb.to_poly.apply(xk),
                y: nb.to_poly.apply(yk),
            },
            a: self.mm.mul(state.a, multiplier),
            b: self.mm.mul(state.b, multiplier),
        }
    }

    /// The jump table drawn from `jump_seed` (`s` uniform in `[1, r)`, `s·G ≠ O`).
    pub fn jumps(&self, jump_seed: u64, charges: &mut StrongRhoCharges) -> Vec<Jump> {
        let mut rng = StdRng::seed_from_u64(jump_seed);
        (0..JUMPS)
            .map(|_| loop {
                let s = rng.gen_range(1..self.modulus);
                charges.scalar_multiplications += 1;
                let point = self.scalar_mul(self.generator, s);
                if point != RawPoint::Infinity {
                    break Jump { point, a: s };
                }
            })
            .collect()
    }

    /// Solve `[d]G = Q` for one target `q` in the prime subgroup.
    ///
    /// `jumps` is the jump table (see [`StrongRho::jumps`]); walk starts draw
    /// `stride ∈ [1, r)` and `cursor ∈ [0, r)` from `start_rng` and then step
    /// `c·G + Q` by `stride·G`. Returns `None` when the step cap is reached.
    pub fn solve(
        &self,
        q: RawPoint,
        jumps: &[Jump],
        start_rng: &mut StdRng,
        params: &StrongRhoParams,
        mut charges: StrongRhoCharges,
    ) -> Option<StrongRhoOutcome> {
        assert_eq!(jumps.len(), JUMPS, "jump table must have JUMPS entries");
        assert!(params.dp_bits < 32, "dp_bits must be below 32");
        let modulus = self.modulus;
        let stride_a = start_rng.gen_range(1..modulus);
        let stride = self.scalar_mul(self.generator, stride_a);
        let mut cursor_a = start_rng.gen_range(0..modulus);
        let mut cursor = self.add(self.scalar_mul(self.generator, cursor_a), q);
        charges.scalar_multiplications += 2;
        charges.group_additions += 1;

        let dp_mask = (1u64 << params.dp_bits) - 1;
        let walk_cap = 8u64 << params.dp_bits;
        let ideal = self.ideal_steps();
        let step_cap = (ideal.ceil() as u64)
            .saturating_mul(params.step_cap_factor)
            .max(1_000_000);
        let lanes_n = params.lanes.max(1);
        let mut table: HashMap<(u8, u64, u64), (u64, u64)> = HashMap::new();
        let mut lanes = vec![
            Lane {
                state: WalkState {
                    point: RawPoint::Infinity,
                    a: 0,
                    b: 0,
                },
                previous: [RawPoint::Infinity; 4],
                length: 0,
                live: false,
            };
            lanes_n
        ];
        let mut denominators = vec![0u64; lanes_n];
        let mut jump_index = vec![0usize; lanes_n];
        let mut active: Vec<usize> = Vec::with_capacity(lanes_n);
        let mut scratch: Vec<u64> = Vec::new();
        let (mut steps, mut walks, mut fruitless, mut capped, mut wasted) = (0u64, 0u64, 0, 0, 0);

        while steps < step_cap {
            // Start a walk in every idle lane.
            for lane in lanes.iter_mut().filter(|lane| !lane.live) {
                loop {
                    walks += 1;
                    let start = cursor;
                    let start_a = cursor_a;
                    cursor = self.add(cursor, stride);
                    cursor_a = add_mod(cursor_a, stride_a, modulus);
                    charges.group_additions += 1;
                    if start == RawPoint::Infinity {
                        continue;
                    }
                    charges.canonicalizations += 1;
                    lane.state = self.canonicalize(WalkState {
                        point: start,
                        a: start_a,
                        b: 1,
                    });
                    lane.previous = [RawPoint::Infinity; 4];
                    lane.length = 0;
                    lane.live = true;
                    break;
                }
            }
            // Lanes standing on a distinguished point.
            for lane in lanes.iter_mut() {
                if !(lane.live && dp_hash(lane.state.point) & dp_mask == 0) {
                    continue;
                }
                lane.live = false;
                let key = lane.state.point.key();
                charges.table_queries += 1;
                let Some(&(hit_a, hit_b)) = table.get(&key) else {
                    table.insert(key, (lane.state.a, lane.state.b));
                    charges.table_inserts += 1;
                    continue;
                };
                let denominator = sub_mod(lane.state.b, hit_b, modulus);
                let Some(inverse) = inverse_mod(denominator, modulus) else {
                    wasted += 1;
                    continue;
                };
                let candidate = mul_mod(sub_mod(hit_a, lane.state.a, modulus), inverse, modulus);
                charges.scalar_multiplications += 1;
                if self.scalar_mul(self.generator, candidate) == q {
                    return Some(StrongRhoOutcome {
                        scalar: candidate,
                        walk_steps: steps,
                        walks,
                        fruitless,
                        capped,
                        wasted_merges: wasted,
                        table_entries: table.len(),
                        ideal_steps: ideal,
                        automorphisms: self.automorphisms,
                        charges,
                    });
                }
                charges.failed_collisions += 1;
            }
            // One step for every live lane, denominators inverted together.
            active.clear();
            for (li, lane) in lanes.iter().enumerate() {
                if !lane.live {
                    continue;
                }
                let ji = partition(lane.state.point);
                charges.partition_hashes += 1;
                let j = active.len();
                jump_index[j] = ji;
                denominators[j] = match (lane.state.point, jumps[ji].point) {
                    (RawPoint::Affine { x: x1, .. }, RawPoint::Affine { x: x2, .. }) => x1 ^ x2,
                    _ => 0,
                };
                active.push(li);
            }
            self.gf
                .batch_inv(&mut denominators[..active.len()], &mut scratch);
            for (j, &li) in active.iter().enumerate() {
                let lane = &mut lanes[li];
                let jump = &jumps[jump_index[j]];
                let sum = match (lane.state.point, jump.point) {
                    (RawPoint::Affine { x: x1, y: y1 }, RawPoint::Affine { x: x2, y: y2 })
                        if denominators[j] != 0 =>
                    {
                        let slope = self.mul(y1 ^ y2, denominators[j]);
                        self.add_with_slope(x1, y1, x2, slope)
                    }
                    _ => self.add(lane.state.point, jump.point),
                };
                charges.group_additions += 1;
                charges.canonicalizations += 1;
                let next = self.canonicalize(WalkState {
                    point: sum,
                    a: add_mod(lane.state.a, jump.a, modulus),
                    b: lane.state.b,
                });
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
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    /// Koblitz rungs that construct, small enough to solve in a unit test.
    fn curves() -> Vec<KoblitzCurve> {
        let mut out = Vec::new();
        for n in [7u32, 11, 13, 17, 19, 23, 29, 31, 37, 41] {
            for a in [0u8, 1] {
                if let Some(curve) = KoblitzCurve::new(a, n) {
                    out.push(curve);
                }
            }
        }
        assert!(out.len() >= 8, "only {} rungs constructed", out.len());
        out
    }

    fn lib_point(curve: &KoblitzCurve, k: u64) -> RawPoint {
        RawPoint::from_binary(&curve.mul(curve.generator(), &BigUint::from(k)))
    }

    #[test]
    fn group_law_matches_the_library_curve() {
        for curve in curves() {
            let rho = StrongRho::new(&curve);
            let g = rho.generator();
            for k in [0u64, 1, 2, 3, 37, 91, 1000, rho.modulus() - 1] {
                assert_eq!(
                    rho.scalar_mul(g, k),
                    lib_point(&curve, k),
                    "n={} k={k}",
                    curve.n
                );
            }
            let p = lib_point(&curve, 37);
            let r = lib_point(&curve, 91);
            assert_eq!(rho.add(p, r), lib_point(&curve, 128));
            assert_eq!(rho.add(p, p), lib_point(&curve, 74));
        }
    }

    #[test]
    fn inverse_is_an_inverse() {
        for curve in curves() {
            let rho = StrongRho::new(&curve);
            let mask = (1u64 << curve.n) - 1;
            let mut rng = StdRng::seed_from_u64(curve.n as u64);
            for _ in 0..200 {
                let v = rng.gen::<u64>() & mask;
                if v != 0 {
                    assert_eq!(rho.mul(v, rho.inv(v)), 1, "n={}", curve.n);
                }
            }
        }
    }

    #[test]
    fn mul_mod_paths_agree() {
        let mut rng = StdRng::seed_from_u64(7);
        for m in [3u64, 1 << 20, (1 << 49) + 1, (1 << 50) - 3, (1 << 57) + 13] {
            let mm = MulMod::new(m);
            for _ in 0..20_000 {
                let (a, b) = (rng.gen_range(0..m), rng.gen_range(0..m));
                assert_eq!(mm.mul(a, b), mul_mod(a, b, m), "m={m}");
            }
            assert_eq!(mm.mul(m - 1, m - 1), mul_mod(m - 1, m - 1, m));
        }
    }

    /// The representative is a function of the orbit, and the multiplier is
    /// right: canon(σP) is the same point for every signed Frobenius image σP,
    /// and if P = a·G + b·Q then canon(P) = a'·G + b'·Q with (a', b') returned.
    #[test]
    fn canonical_form_is_orbit_invariant_and_multiplier_is_correct() {
        for curve in curves() {
            let rho = StrongRho::new(&curve);
            let r = rho.modulus();
            let d = r / 3 + 1;
            let q = rho.scalar_mul(rho.generator(), d);
            let mut rng = StdRng::seed_from_u64(1000 + curve.n as u64);
            for _ in 0..40 {
                let (a, b) = (rng.gen_range(0..r), rng.gen_range(1..r));
                let p = rho.add(rho.scalar_mul(rho.generator(), a), rho.scalar_mul(q, b));
                if p == RawPoint::Infinity {
                    continue;
                }
                let canon = rho.canonicalize(WalkState { point: p, a, b });
                let expected = rho.add(
                    rho.scalar_mul(rho.generator(), canon.a),
                    rho.scalar_mul(q, canon.b),
                );
                assert_eq!(canon.point, expected, "multiplier wrong at n={}", curve.n);
                // every orbit member maps to the same representative
                let mut image = p;
                for _ in 0..curve.n {
                    let RawPoint::Affine { x, y } = image else {
                        unreachable!()
                    };
                    for member in [image, RawPoint::Affine { x, y: y ^ x }] {
                        let c = rho.canonicalize(WalkState {
                            point: member,
                            a: 0,
                            b: 1,
                        });
                        assert_eq!(c.point, canon.point, "not orbit-invariant at n={}", curve.n);
                    }
                    image = RawPoint::Affine {
                        x: rho.sqr(x),
                        y: rho.sqr(y),
                    };
                }
            }
        }
    }

    #[test]
    fn solves_planted_targets_at_every_lane_width() {
        for curve in curves().into_iter().filter(|c| (13..=31).contains(&c.n)) {
            let rho = StrongRho::new(&curve);
            for (lanes, dp_bits) in [(1usize, 2u32), (8, 3), (32, 4)] {
                let params = StrongRhoParams {
                    lanes,
                    dp_bits,
                    ..StrongRhoParams::default()
                };
                for seed in 0..4u64 {
                    let mut charges = StrongRhoCharges::default();
                    let jumps = rho.jumps(seed ^ 0x5eed, &mut charges);
                    let mut rng = StdRng::seed_from_u64(seed);
                    let d = rng.gen_range(1..rho.modulus());
                    let q = rho.scalar_mul(rho.generator(), d);
                    let out = rho
                        .solve(q, &jumps, &mut rng, &params, charges)
                        .unwrap_or_else(|| panic!("no solve n={} lanes={lanes}", curve.n));
                    assert_eq!(out.scalar, d, "n={} lanes={lanes} seed={seed}", curve.n);
                    assert!(out.walk_steps > 0 && out.table_entries > 0);
                }
            }
        }
    }

    /// A strong rho should need about √(πr/2A) steps; at n = 37/41 the median
    /// over a few seeds stays within a generous factor of that.
    #[test]
    fn step_counts_are_near_the_ideal() {
        for curve in curves().into_iter().filter(|c| c.n >= 37) {
            let rho = StrongRho::new(&curve);
            let mut ratios = Vec::new();
            for seed in 0..9u64 {
                let mut charges = StrongRhoCharges::default();
                let jumps = rho.jumps(seed + 77, &mut charges);
                let mut rng = StdRng::seed_from_u64(seed + 99);
                let d = rng.gen_range(1..rho.modulus());
                let q = rho.scalar_mul(rho.generator(), d);
                let out = rho
                    .solve(q, &jumps, &mut rng, &StrongRhoParams::default(), charges)
                    .expect("solve");
                assert_eq!(out.scalar, d);
                ratios.push(out.walk_steps as f64 / out.ideal_steps);
            }
            ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let median = ratios[ratios.len() / 2];
            assert!(
                (0.3..4.0).contains(&median),
                "n={} median steps/ideal {median}",
                curve.n
            );
        }
    }
}
