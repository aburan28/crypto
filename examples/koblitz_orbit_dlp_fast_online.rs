//! Compact-orbit relation extraction with a shared-factor-log DLP pass.
//!
//! Public synthetic Koblitz fixtures only. Reads a retained
//! `point_defined_factor_base` header, builds the Frobenius-quotiented S3
//! regular-root index once (no pair table, no edge selectors), then:
//!
//! 1. rank stage: decomposes random `[a]G` until the K orbit logs are fixed,
//! 2. target stage: decomposes each published target and recovers its log.
//!
//! Every stage is timed in this one process.
//!
//! Frobenius is a bit rotation in a normal basis, so a demanded partner
//! x-coordinate needs one canonical-rotation probe instead of n probes.
//!
//! Single-target online version of `koblitz_orbit_dlp_fast.rs` (which stays
//! frozen for the evidence that pins it). The relation search, rank stage and
//! descent are unchanged; only the target timers differ. Each target's online
//! interval is one chain of laps on a single clock, from the first operation on
//! the target point (`target_query_begin`) to the end of the `[d]G = Q` check
//! (`recovery_check_end`), so the five exclusive phases sum to the interval:
//! query hash, PDP (candidate generation, partner roots, index probes),
//! relation check (every `lift` call: base-point choice and group sum), descent
//! (log combination) and recovery check. Reusable setup (base, index, rank
//! stage, linear algebra) ends before the first target's interval starts.
//!
//! Usage: <base_header.jsonl> <target_points.jsonl|legacy_scalars.txt> <rank_seed> <out.jsonl>

#![recursion_limit = "256"]

use crypto_lib::cryptanalysis::ecbench::isolation::HwCounters;
use crypto_lib::cryptanalysis::koblitz_fast_arith::{s3_x_roots, FastBinaryCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::collections::HashMap;
use std::io::Write;
use std::process::ExitCode;
use std::time::{Duration, Instant};

/// Consecutive laps on one clock: the laps of an interval sum to its length.
struct PhaseClock {
    last: Instant,
}

impl PhaseClock {
    fn start() -> Self {
        Self {
            last: Instant::now(),
        }
    }

    fn lap(&mut self) -> Duration {
        let now = Instant::now();
        let elapsed = now - self.last;
        self.last = now;
        elapsed
    }
}

/// Exclusive PDP / relation-check time of one `extract` call.
#[derive(Default)]
struct ExtractSplit {
    pdp: Duration,
    relation_check: Duration,
    relation_checks: u64,
    root_calls: u64,
    skipped_symmetric_states: u64,
}

fn ms(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1000.0
}

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

/// Table-driven S3 root solver: S3(u, v, t) = 0 as a quadratic in v.
struct S3Solver {
    square: Linear,
    half_trace: Linear,
    trace_mask: u64,
    b: u64,
}

impl S3Solver {
    fn new(gf: &Gf2, b: u64) -> Self {
        let n = gf.n as usize;
        assert!(n % 2 == 1, "half-trace solver needs odd n");
        let unit: Vec<u64> = (0..n).map(|bit| 1u64 << bit).collect();
        let half_trace_images: Vec<u64> = unit
            .iter()
            .map(|&e| {
                let mut acc = e;
                let mut power = e;
                for _ in 0..(n - 1) / 2 {
                    power = gf.sqr(gf.sqr(power));
                    acc ^= power;
                }
                acc
            })
            .collect();
        let mut trace_mask = 0u64;
        for (bit, &e) in unit.iter().enumerate() {
            let mut acc = e;
            let mut power = e;
            for _ in 1..n {
                power = gf.sqr(power);
                acc ^= power;
            }
            if acc & 1 == 1 {
                trace_mask |= 1u64 << bit;
            }
        }
        Self {
            square: Linear::from_images(&unit.iter().map(|&e| gf.sqr(e)).collect::<Vec<_>>()),
            half_trace: Linear::from_images(&half_trace_images),
            trace_mask,
            b,
        }
    }

    /// Same roots, in the same order, as `s3_x_roots` on the generic branch.
    #[inline(always)]
    fn roots(&self, gf: &Gf2, basis: &NormalBasis, left: u64, right: u64) -> Option<[u64; 2]> {
        let s = left ^ right;
        let a = self.square.apply(s);
        let p = gf.mul(left, right);
        let ps = self.square.apply(p);
        if a == 0 || ps == 0 {
            return s3_x_roots(gf, self.b, left, right);
        }
        let inv_combined = invert(gf, basis, gf.mul(a, ps));
        let q = gf.mul(p, gf.mul(ps, inv_combined));
        let c = ps ^ self.b;
        let d = gf.mul(gf.mul(c, a), gf.mul(a, inv_combined));
        if (d & self.trace_mask).count_ones() & 1 == 1 {
            return None;
        }
        let first = gf.mul(q, self.half_trace.apply(d));
        Some([first, first ^ q])
    }
}

/// Itoh-Tsujii inversion with the `a^(2^k)` steps done as normal-basis rotations.
#[inline(always)]
fn invert(gf: &Gf2, basis: &NormalBasis, a: u64) -> u64 {
    let e = gf.n - 1;
    let mut c = a;
    let mut k = 1u32;
    let bits = 32 - e.leading_zeros();
    for i in (0..bits - 1).rev() {
        let raised = basis
            .to_poly
            .apply(basis.rotate(basis.to_normal.apply(c), k));
        c = gf.mul(c, raised);
        k *= 2;
        if (e >> i) & 1 == 1 {
            c = gf.mul(gf.sqr(c), a);
            k += 1;
        }
    }
    gf.sqr(c)
}

/// Open-addressing u64 -> u64 table; keys are < 2^63 so u64::MAX marks empty.
/// Key and value share a slot so a probe touches one cache line.
/// A negative answer skips an exact root-table probe. Every inserted key is
/// also inserted here, so the filter cannot discard a table witness.
struct RootBloom {
    words: Vec<u64>,
    bit_mask: usize,
}

impl RootBloom {
    fn with_capacity(keys: usize) -> Self {
        let bits = keys
            .checked_mul(8)
            .and_then(usize::checked_next_power_of_two)
            .expect("root Bloom size overflow")
            .max(64);
        Self {
            words: vec![0; bits / 64],
            bit_mask: bits - 1,
        }
    }

    #[inline(always)]
    fn positions(&self, key: u64) -> [usize; 3] {
        let first = mix_root_bloom(key ^ 0xA5A5_38F9_9328_6F3D);
        let step = mix_root_bloom(key.rotate_left(17) ^ 0xD4C9_3E72_B719_485B) | 1;
        [0, 1, 2].map(|index| first.wrapping_add(index * step) as usize & self.bit_mask)
    }

    #[inline(always)]
    fn insert(&mut self, key: u64) {
        for bit in self.positions(key) {
            self.words[bit >> 6] |= 1u64 << (bit & 63);
        }
    }

    #[inline(always)]
    fn contains(&self, key: u64) -> bool {
        self.positions(key)
            .into_iter()
            .all(|bit| self.words[bit >> 6] & (1u64 << (bit & 63)) != 0)
    }

    fn bytes(&self) -> usize {
        self.words.len() * std::mem::size_of::<u64>()
    }
}

#[inline(always)]
fn mix_root_bloom(mut value: u64) -> u64 {
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

struct RootTable {
    slots: Vec<(u64, u64)>,
    shift: u32,
    len: usize,
    bloom: Option<RootBloom>,
}

impl RootTable {
    fn with_capacity(entries: usize) -> Self {
        let slots = (entries * 2).next_power_of_two().max(16);
        Self {
            slots: vec![(u64::MAX, 0); slots],
            shift: 64 - slots.trailing_zeros(),
            len: 0,
            bloom: (std::env::var("KIC_ROOT_BLOOM").as_deref() == Ok("1"))
                .then(|| RootBloom::with_capacity(entries)),
        }
    }

    #[inline(always)]
    fn slot(&self, key: u64) -> usize {
        (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> self.shift) as usize
    }

    fn insert_if_absent(&mut self, key: u64, value: u64) {
        let mask = self.slots.len() - 1;
        let mut slot = self.slot(key);
        loop {
            if self.slots[slot].0 == u64::MAX {
                self.slots[slot] = (key, value);
                self.len += 1;
                if let Some(bloom) = &mut self.bloom {
                    bloom.insert(key);
                }
                return;
            }
            if self.slots[slot].0 == key {
                return;
            }
            slot = (slot + 1) & mask;
        }
    }

    #[inline(always)]
    fn get(&self, key: u64) -> Option<u64> {
        if self
            .bloom
            .as_ref()
            .is_some_and(|bloom| !bloom.contains(key))
        {
            return None;
        }
        let mask = self.slots.len() - 1;
        let mut slot = self.slot(key);
        loop {
            let (stored, value) = self.slots[slot];
            if stored == key {
                return Some(value);
            }
            if stored == u64::MAX {
                return None;
            }
            slot = (slot + 1) & mask;
        }
    }
}

#[cfg(test)]
mod root_bloom_table_tests {
    use super::{mix_root_bloom, RootBloom, RootTable};

    #[test]
    fn screened_lookup_matches_exact_table() {
        let mut reference = RootTable::with_capacity(4096);
        reference.bloom = None;
        let mut screened = RootTable::with_capacity(4096);
        screened.bloom = Some(RootBloom::with_capacity(4096));
        for i in 0..4096u64 {
            let key = mix_root_bloom(i) & (u64::MAX >> 1);
            reference.insert_if_absent(key, i);
            screened.insert_if_absent(key, i);
            screened.insert_if_absent(key, i + 1);
        }
        for i in 0..100_000u64 {
            let key = mix_root_bloom(i ^ 0xD1B5_4A32_D192_ED03) & (u64::MAX >> 1);
            assert_eq!(screened.get(key), reference.get(key));
        }
        for &(key, value) in &reference.slots {
            if key != u64::MAX {
                assert!(screened.bloom.as_ref().unwrap().contains(key));
                assert_eq!(screened.get(key), Some(value));
            }
        }
    }
}

#[derive(Clone, Copy)]
struct State {
    left: u16,
    right: u16,
    relative: u16,
    normal_roots: [u64; 2],
}

struct Index {
    states: Vec<State>,
    table: RootTable,
    shifted: Vec<Vec<u64>>,
    pair_candidates: usize,
    counterparts: Option<Vec<u32>>,
}

fn pack(left: u16, right: u16, relative: u16, shift: u32) -> u64 {
    (left as u64) | ((right as u64) << 16) | ((relative as u64) << 32) | ((shift as u64) << 48)
}

fn unpack(value: u64) -> (usize, usize, usize, u32) {
    (
        (value & 0xffff) as usize,
        ((value >> 16) & 0xffff) as usize,
        ((value >> 32) & 0xffff) as usize,
        ((value >> 48) & 0xffff) as u32,
    )
}

fn build_index(
    gf: &Gf2,
    basis: &NormalBasis,
    solver: &S3Solver,
    reps: &[u64],
    pair_symmetry: bool,
) -> Index {
    let n = gf.n as usize;
    let shifted: Vec<Vec<u64>> = reps
        .iter()
        .map(|&code| {
            let mut row = Vec::with_capacity(n);
            let mut value = code;
            for _ in 0..n {
                row.push(value);
                value = gf.sqr(value);
            }
            row
        })
        .collect();
    let width = reps.len();
    let slot = |left: usize, right: usize, relative: usize| (left * width + right) * n + relative;
    let mate_slot =
        |left: usize, right: usize, relative: usize| slot(right, left, (n - relative) % n);
    let mut pair_candidates = 0;
    let (states, counterparts) = if pair_symmetry {
        // Build each unordered pair orbit once, then restore its swapped
        // alias in the original scan position without another S3 solve.
        // This keeps the public-query scan origin identical to the control.
        let mut slots: Vec<Option<State>> = vec![None; width * width * n];
        for left in 0..width {
            for right in left..width {
                for relative in 0..n {
                    if left == right && relative > n - relative {
                        continue;
                    }
                    pair_candidates += 1;
                    if let Some(roots) =
                        solver.roots(gf, basis, shifted[left][0], shifted[right][relative])
                    {
                        let normal_roots = [
                            basis.to_normal.apply(roots[0]),
                            basis.to_normal.apply(roots[1]),
                        ];
                        let original = slot(left, right, relative);
                        let swapped = mate_slot(left, right, relative);
                        assert!(slots[original].is_none());
                        slots[original] = Some(State {
                            left: left as u16,
                            right: right as u16,
                            relative: relative as u16,
                            normal_roots,
                        });
                        if swapped != original {
                            assert!(slots[swapped].is_none());
                            slots[swapped] = Some(State {
                                left: right as u16,
                                right: left as u16,
                                relative: ((n - relative) % n) as u16,
                                normal_roots: [
                                    basis.rotate(normal_roots[0], ((n - relative) % n) as u32),
                                    basis.rotate(normal_roots[1], ((n - relative) % n) as u32),
                                ],
                            });
                        }
                    }
                }
            }
        }
        let mut positions = vec![u32::MAX; slots.len()];
        let mut states = Vec::new();
        for (candidate_slot, state) in slots.into_iter().enumerate() {
            if let Some(state) = state {
                positions[candidate_slot] = states.len() as u32;
                states.push(state);
            }
        }
        let counterparts = states
            .iter()
            .map(|state| {
                let mate = mate_slot(
                    state.left as usize,
                    state.right as usize,
                    state.relative as usize,
                );
                let position = positions[mate];
                assert_ne!(position, u32::MAX);
                position
            })
            .collect();
        (states, Some(counterparts))
    } else {
        let mut states = Vec::new();
        for left in 0..width {
            for right in 0..width {
                for relative in 0..n {
                    pair_candidates += 1;
                    if let Some(roots) =
                        solver.roots(gf, basis, shifted[left][0], shifted[right][relative])
                    {
                        states.push(State {
                            left: left as u16,
                            right: right as u16,
                            relative: relative as u16,
                            normal_roots: [
                                basis.to_normal.apply(roots[0]),
                                basis.to_normal.apply(roots[1]),
                            ],
                        });
                    }
                }
            }
        }
        (states, None)
    };
    let mut table = RootTable::with_capacity(pair_candidates * 2);
    for (index, state) in states.iter().enumerate() {
        if counterparts
            .as_ref()
            .is_some_and(|partners: &Vec<u32>| partners[index] < index as u32)
        {
            continue;
        }
        for &root in &state.normal_roots {
            let (canonical, shift) = basis.canonical(root);
            table.insert_if_absent(
                canonical,
                pack(state.left, state.right, state.relative, shift),
            );
        }
    }
    Index {
        states,
        table,
        shifted,
        pair_candidates,
        counterparts,
    }
}

struct Relation {
    point_indices: [usize; 4],
    x_codes: [u64; 4],
    intermediates: [u64; 2],
    probes: u64,
}

struct ExtractOutcome {
    relation: Option<Relation>,
    probes: u64,
    capped: bool,
}

#[derive(Clone, Copy)]
enum TargetInput {
    PublicPoint(FastPoint),
    KnownAnswerScalar(u64),
}

struct Base {
    points: Vec<FastPoint>,
    labels: Vec<(usize, u64)>,
    by_x: HashMap<u64, Vec<usize>>,
    columns: usize,
}

fn lift(
    fast: &FastBinaryCurve,
    base: &Base,
    codes: &[u64; 4],
    target: FastPoint,
) -> Option<[usize; 4]> {
    let choices: Vec<&Vec<usize>> = codes
        .iter()
        .map(|code| base.by_x.get(code))
        .collect::<Option<_>>()?;
    for &a in choices[0] {
        for &b in choices[1] {
            let ab = fast.add(base.points[a], base.points[b]);
            for &c in choices[2] {
                let abc = fast.add(ab, base.points[c]);
                for &d in choices[3] {
                    if fast.add(abc, base.points[d]) == target {
                        return Some([a, b, c, d]);
                    }
                }
            }
        }
    }
    None
}

fn extract(
    gf: &Gf2,
    fast: &FastBinaryCurve,
    basis: &NormalBasis,
    solver: &S3Solver,
    index: &Index,
    base: &Base,
    target: FastPoint,
    start: usize,
    probe_cap: Option<u64>,
    clock: &mut PhaseClock,
    split: &mut ExtractSplit,
) -> ExtractOutcome {
    let relation = extract_timed(
        gf, fast, basis, solver, index, base, target, start, probe_cap, clock, split,
    );
    split.pdp += clock.lap();
    relation
}

/// The search of `koblitz_orbit_dlp_fast.rs`; laps are taken only around
/// `lift`, which runs on index hits, so the probe loop carries no timer.
#[allow(clippy::too_many_arguments)]
fn extract_timed(
    gf: &Gf2,
    fast: &FastBinaryCurve,
    basis: &NormalBasis,
    solver: &S3Solver,
    index: &Index,
    base: &Base,
    target: FastPoint,
    start: usize,
    probe_cap: Option<u64>,
    clock: &mut PhaseClock,
    split: &mut ExtractSplit,
) -> ExtractOutcome {
    let Some((target_x, _)) = target else {
        return ExtractOutcome {
            relation: None,
            probes: 0,
            capped: false,
        };
    };
    let n = gf.n;
    let mut probes = 0u64;
    // Callers provide a rotating, target-dependent origin so rank-stage rows
    // do not concentrate on the same early columns.
    let len = index.states.len();
    let start = start % len.max(1);
    let first = index.states[start..]
        .iter()
        .enumerate()
        .map(|(offset, state)| (start + offset, state));
    let second = index.states[..start].iter().enumerate();
    for (state_index, state) in first.chain(second) {
        if let Some(counterparts) = &index.counterparts {
            let mate_index = counterparts[state_index] as usize;
            let distance = (state_index + len - start) % len;
            let mate_distance = (mate_index + len - start) % len;
            if mate_distance < distance {
                split.skipped_symmetric_states += 1;
                continue;
            }
        }
        for shift in 0..n {
            let left_x = index.shifted[state.left as usize][shift as usize];
            let right_x = index.shifted[state.right as usize]
                [(shift as usize + state.relative as usize) % n as usize];
            for &normal_root in &state.normal_roots {
                let absolute = basis.to_poly.apply(basis.rotate(normal_root, shift));
                split.root_calls += 1;
                let Some(partners) = solver.roots(gf, basis, absolute, target_x) else {
                    continue;
                };
                for partner in partners {
                    probes += 1;
                    let (canonical, partner_shift) =
                        basis.canonical(basis.to_normal.apply(partner));
                    let Some(value) = index.table.get(canonical) else {
                        if probe_cap.is_some_and(|cap| probes >= cap) {
                            return ExtractOutcome {
                                relation: None,
                                probes,
                                capped: true,
                            };
                        }
                        continue;
                    };
                    let (left2, right2, relative2, stored_shift) = unpack(value);
                    let left_shift2 = ((stored_shift + n - partner_shift) % n) as usize;
                    let right_shift2 = (left_shift2 + relative2) % n as usize;
                    let codes = [
                        left_x,
                        right_x,
                        index.shifted[left2][left_shift2],
                        index.shifted[right2][right_shift2],
                    ];
                    split.pdp += clock.lap();
                    let lifted = lift(fast, base, &codes, target);
                    split.relation_check += clock.lap();
                    split.relation_checks += 1;
                    if let Some(point_indices) = lifted {
                        return ExtractOutcome {
                            relation: Some(Relation {
                                point_indices,
                                x_codes: codes,
                                intermediates: [absolute, partner],
                                probes,
                            }),
                            probes,
                            capped: false,
                        };
                    }
                    if probe_cap.is_some_and(|cap| probes >= cap) {
                        return ExtractOutcome {
                            relation: None,
                            probes,
                            capped: true,
                        };
                    }
                }
            }
        }
    }
    ExtractOutcome {
        relation: None,
        probes,
        capped: false,
    }
}

/// Frobenius-closed signed-orbit base from a deterministic public x-scan.
///
/// Each accepted x lifts to a point, is projected into the order-r subgroup by
/// the cofactor, and contributes its whole signed Frobenius orbit. Point logs
/// are unknown; the label of member `phi^k(+-R_j)` is `(j, +-lambda^k)`.
fn construct_base(
    fast: &FastBinaryCurve,
    curve: &KoblitzCurve,
    b: u64,
    columns: usize,
    r: u64,
) -> (
    Vec<FastPoint>,
    Vec<(usize, u64)>,
    Vec<Option<[u64; 2]>>,
    u64,
) {
    let gf = &fast.gf;
    let n = fast.n as usize;
    let lambda = curve.lambda.to_u64().unwrap();
    let mut seen = std::collections::HashSet::new();
    let mut points = Vec::with_capacity(columns * 2 * n);
    let mut labels = Vec::with_capacity(columns * 2 * n);
    let mut reps = Vec::with_capacity(columns);
    let mut scanned = 0u64;
    let mut raw_x = 1u64;
    while reps.len() < columns {
        scanned += 1;
        for point in fast.points_with_x(b, raw_x) {
            let Some((px, py)) = fast.scalar_mul(point, &curve.cofactor) else {
                continue;
            };
            if px <= 1 {
                continue;
            }
            let mut key = px;
            let mut x = px;
            for _ in 1..n {
                x = gf.sqr(x);
                key = key.min(x);
            }
            if !seen.insert(key) {
                continue;
            }
            let column = reps.len();
            reps.push(Some([px, py]));
            let (mut cx, mut cy) = (px, py);
            let mut coefficient = 1u64;
            for _ in 0..n {
                points.push(Some((cx, cy)));
                labels.push((column, coefficient));
                points.push(Some((cx, cx ^ cy)));
                labels.push((column, (r - coefficient) % r));
                cx = gf.sqr(cx);
                cy = gf.sqr(cy);
                coefficient = mulmod(coefficient, lambda, r);
            }
            if reps.len() == columns {
                break;
            }
        }
        raw_x += 1;
    }
    let rep = reps[0].map(|[x, y]| (x, y));
    assert_eq!(
        fast.scalar_mul(rep, &BigUint::from(lambda)),
        rep.map(|(x, y)| (gf.sqr(x), gf.sqr(y))),
        "Frobenius must act as lambda on the subgroup"
    );
    (points, labels, reps, scanned)
}

fn mulmod(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 * b as u128) % r as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, r: u64) -> u64 {
    let mut acc = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            acc = mulmod(acc, a, r);
        }
        a = mulmod(a, a, r);
        e >>= 1;
    }
    acc
}

/// Incremental row echelon mod prime r over K unknowns plus a right-hand side.
struct Echelon {
    r: u64,
    columns: usize,
    pivots: Vec<Option<Vec<u64>>>,
    rank: usize,
}

impl Echelon {
    fn insert(&mut self, mut row: Vec<u64>) -> bool {
        for column in 0..self.columns {
            if row[column] == 0 {
                continue;
            }
            if let Some(pivot) = &self.pivots[column] {
                let factor = row[column];
                for (value, &p) in row.iter_mut().zip(pivot.iter()).skip(column) {
                    *value = (*value + self.r - mulmod(factor, p, self.r)) % self.r;
                }
            } else {
                let inverse = powmod(row[column], self.r - 2, self.r);
                for value in row.iter_mut().skip(column) {
                    *value = mulmod(*value, inverse, self.r);
                }
                self.pivots[column] = Some(row);
                self.rank += 1;
                return true;
            }
        }
        false
    }

    fn solve(&self) -> Vec<u64> {
        let mut solution = vec![0u64; self.columns];
        for column in (0..self.columns).rev() {
            let pivot = self.pivots[column].as_ref().expect("full rank");
            let mut value = pivot[self.columns];
            for other in column + 1..self.columns {
                value = (value + self.r - mulmod(pivot[other], solution[other], self.r)) % self.r;
            }
            solution[column] = value;
        }
        solution
    }
}

fn relation_row(base: &Base, relation: &Relation, rhs: u64, r: u64) -> Vec<u64> {
    let mut row = vec![0u64; base.columns + 1];
    for &index in &relation.point_indices {
        let (column, coefficient) = base.labels[index];
        row[column] = (row[column] + coefficient) % r;
    }
    row[base.columns] = rhs % r;
    row
}

fn peak_rss_bytes() -> Option<u64> {
    #[cfg(unix)]
    {
        let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
        // SAFETY: getrusage writes the complete rusage value on success.
        if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
            return None;
        }
        // SAFETY: the successful getrusage call initialized usage above.
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

// ── Wide (`u128`-word) compact-orbit pipeline for `64 < n ≤ 127` ──
//
// Twin of the `u64` path above with identical stage semantics and the same
// online-interval instrumentation (five phase timers on one clock, hardware
// counters around the whole target interval).  Subgroup scalars, orbit
// indices, packed witnesses, and the rank linear algebra stay `u64`; only
// field elements move to the wider word.
mod wide {
    use super::{mix_root_bloom, ms, pack, unpack, Echelon, ExtractSplit, HwCounters, PhaseClock};
    use crypto_lib::cryptanalysis::koblitz_fast_arith::{
        s3_x_roots_128, FastBinaryCurve128, FastPoint128,
    };
    use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
    use crypto_lib::cryptanalysis::semaev_decomp::Gf2_128;
    use num_bigint::BigUint;
    use num_traits::ToPrimitive;
    use serde_json::{json, Value};
    use std::collections::HashMap;
    use std::io::Write;
    use std::time::Instant;

    /// F2-linear map on n ≤ 128 bits applied by byte tables.
    struct Linear128 {
        tables: Vec<[u128; 256]>,
    }

    impl Linear128 {
        fn from_images(images: &[u128]) -> Self {
            let chunks = images.len().div_ceil(8);
            let mut tables = vec![[0u128; 256]; chunks];
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
        fn apply(&self, v: u128) -> u128 {
            let mut acc = 0u128;
            for (chunk, table) in self.tables.iter().enumerate() {
                acc ^= table[((v >> (8 * chunk)) & 0xff) as usize];
            }
            acc
        }
    }

    struct NormalBasis128 {
        n: u32,
        mask: u128,
        to_normal: Linear128,
        to_poly: Linear128,
        sparse_canonical: bool,
    }

    impl NormalBasis128 {
        fn new(gf: &Gf2_128) -> Self {
            let n = gf.n as usize;
            let mask = (1u128 << n) - 1;
            for candidate in 2u128.. {
                let mut conjugates = Vec::with_capacity(n);
                let mut value = candidate & mask;
                for _ in 0..n {
                    conjugates.push(value);
                    value = gf.sqr(value);
                }
                let Some(inverse_columns) = invert_columns128(&conjugates, n) else {
                    continue;
                };
                return Self {
                    n: gf.n,
                    mask,
                    to_normal: Linear128::from_images(&inverse_columns),
                    to_poly: Linear128::from_images(&conjugates),
                    sparse_canonical: std::env::var("KIC_SPARSE_CANONICAL").as_deref() == Ok("1"),
                };
            }
            unreachable!()
        }

        #[inline(always)]
        fn rotate(&self, v: u128, k: u32) -> u128 {
            if k == 0 {
                return v;
            }
            ((v << k) | (v >> (self.n - k))) & self.mask
        }

        #[inline(always)]
        fn canonical(&self, v: u128) -> (u128, u32) {
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

        /// Wide-word analogue of the exact cyclic zero-gap selector.
        #[inline(always)]
        fn canonical_sparse(&self, v: u128) -> (u128, u32) {
            let value = v & self.mask;
            if value == 0 || value == self.mask {
                return (value, 0);
            }

            let n = self.n as usize;
            let zeros = (!value) & self.mask;
            let mut zero_runs = [0u128; 7];
            zero_runs[0] = zeros;
            for level in 1..zero_runs.len() {
                let previous_length = 1usize << (level - 1);
                if previous_length >= n {
                    zero_runs[level] = 0;
                    continue;
                }
                zero_runs[level] = zero_runs[level - 1]
                    & self.rotate(zero_runs[level - 1], previous_length as u32);
            }

            let mut longest_gap = 0usize;
            let mut starts = self.mask;
            for level in (0..zero_runs.len()).rev() {
                let length = 1usize << level;
                if longest_gap + length >= n {
                    continue;
                }
                let eligible = starts & self.rotate(zero_runs[level], longest_gap as u32);
                if eligible != 0 {
                    starts = eligible;
                    longest_gap += length;
                }
            }

            let mut best_value = u128::MAX;
            let mut best_shift = u32::MAX;
            while starts != 0 {
                let start = starts.trailing_zeros() as usize;
                starts &= starts - 1;
                let shift = (n - 1 - start) as u32;
                let rotated = self.rotate(value, shift);
                if rotated < best_value || (rotated == best_value && shift < best_shift) {
                    best_value = rotated;
                    best_shift = shift;
                }
            }
            (best_value, best_shift)
        }

        #[inline(always)]
        fn canonical_selected(&self, v: u128) -> (u128, u32) {
            if self.sparse_canonical {
                self.canonical_sparse(v)
            } else {
                self.canonical(v)
            }
        }
    }

    /// Same Gauss-Jordan as the twin, with rows as `(m, i)` word pairs.
    fn invert_columns128(columns: &[u128], n: usize) -> Option<Vec<u128>> {
        let mut rows: Vec<(u128, u128)> = (0..n)
            .map(|r| {
                let mut m = 0u128;
                for (i, &column) in columns.iter().enumerate() {
                    if (column >> r) & 1 == 1 {
                        m |= 1u128 << i;
                    }
                }
                (m, 1u128 << r)
            })
            .collect();
        for pivot in 0..n {
            let found = (pivot..n).find(|&r| (rows[r].0 >> pivot) & 1 == 1)?;
            rows.swap(pivot, found);
            for r in 0..n {
                if r != pivot && (rows[r].0 >> pivot) & 1 == 1 {
                    rows[r].0 ^= rows[pivot].0;
                    rows[r].1 ^= rows[pivot].1;
                }
            }
        }
        let mut images = vec![0u128; n];
        for (i, &(_, ih)) in rows.iter().enumerate() {
            for r in 0..n {
                if (ih >> r) & 1 == 1 {
                    images[r] |= 1u128 << i;
                }
            }
        }
        Some(images)
    }

    /// Table-driven S3 root solver on wide words.
    struct S3Solver128 {
        square: Linear128,
        half_trace: Linear128,
        trace_mask: u128,
        b: u128,
    }

    impl S3Solver128 {
        fn new(gf: &Gf2_128, b: u128) -> Self {
            let n = gf.n as usize;
            assert!(n % 2 == 1, "half-trace solver needs odd n");
            let unit: Vec<u128> = (0..n).map(|bit| 1u128 << bit).collect();
            let half_trace_images: Vec<u128> = unit
                .iter()
                .map(|&e| {
                    let mut acc = e;
                    let mut power = e;
                    for _ in 0..(n - 1) / 2 {
                        power = gf.sqr(gf.sqr(power));
                        acc ^= power;
                    }
                    acc
                })
                .collect();
            let mut trace_mask = 0u128;
            for (bit, &e) in unit.iter().enumerate() {
                let mut acc = e;
                let mut power = e;
                for _ in 1..n {
                    power = gf.sqr(power);
                    acc ^= power;
                }
                if acc & 1 == 1 {
                    trace_mask |= 1u128 << bit;
                }
            }
            Self {
                square: Linear128::from_images(
                    &unit.iter().map(|&e| gf.sqr(e)).collect::<Vec<_>>(),
                ),
                half_trace: Linear128::from_images(&half_trace_images),
                trace_mask,
                b,
            }
        }

        /// Everything before the inversion: `a = (left + right)^2`,
        /// `p = left·right`, `ps = p^2`, and the denominator `a·ps`.
        #[inline(always)]
        fn prepare(&self, gf: &Gf2_128, left: u128, right: u128) -> PreparedRoot128 {
            let a = self.square.apply(left ^ right);
            let p = gf.mul(left, right);
            let ps = self.square.apply(p);
            let denominator = if a == 0 || ps == 0 { 0 } else { gf.mul(a, ps) };
            PreparedRoot128 {
                left,
                right,
                a,
                p,
                ps,
                denominator,
            }
        }

        /// Everything after the inversion, given `inverse = denominator^{-1}`.
        #[inline(always)]
        fn finish(
            &self,
            gf: &Gf2_128,
            prepared: PreparedRoot128,
            inverse: u128,
        ) -> Option<[u128; 2]> {
            if prepared.denominator == 0 {
                return s3_x_roots_128(gf, self.b, prepared.left, prepared.right);
            }
            let q = gf.mul(prepared.p, gf.mul(prepared.ps, inverse));
            let c = prepared.ps ^ self.b;
            let d = gf.mul(gf.mul(c, prepared.a), gf.mul(prepared.a, inverse));
            if (d & self.trace_mask).count_ones() & 1 == 1 {
                return None;
            }
            let first = gf.mul(q, self.half_trace.apply(d));
            Some([first, first ^ q])
        }

        /// Scalar form: one inversion per call, same root order as the
        /// batched form.
        #[inline(always)]
        fn roots(
            &self,
            gf: &Gf2_128,
            basis: &NormalBasis128,
            left: u128,
            right: u128,
        ) -> Option<[u128; 2]> {
            let prepared = self.prepare(gf, left, right);
            let inverse = if prepared.denominator == 0 {
                0
            } else {
                invert128(gf, basis, prepared.denominator)
            };
            self.finish(gf, prepared, inverse)
        }
    }

    #[derive(Clone, Copy)]
    struct PreparedRoot128 {
        left: u128,
        right: u128,
        a: u128,
        p: u128,
        ps: u128,
        denominator: u128,
    }

    /// Reusable scratch for solving many S3 root pairs with one shared
    /// inversion (Montgomery's trick over the non-exceptional
    /// denominators).
    struct RootBatch128 {
        prepared: Vec<PreparedRoot128>,
        inverses: Vec<u128>,
        prefixes: Vec<u128>,
        regular: Vec<usize>,
        results: Vec<Option<[u128; 2]>>,
    }

    impl RootBatch128 {
        fn with_capacity(entries: usize) -> Self {
            Self {
                prepared: Vec::with_capacity(entries),
                inverses: Vec::with_capacity(entries),
                prefixes: Vec::with_capacity(entries),
                regular: Vec::with_capacity(entries),
                results: Vec::with_capacity(entries),
            }
        }

        fn solve(
            &mut self,
            gf: &Gf2_128,
            basis: &NormalBasis128,
            solver: &S3Solver128,
            inputs: &[(u128, u128)],
        ) {
            self.prepared.clear();
            self.inverses.clear();
            self.prefixes.clear();
            self.regular.clear();
            self.results.clear();
            self.inverses.resize(inputs.len(), 0);
            self.prefixes.resize(inputs.len(), 0);
            for (index, &(left, right)) in inputs.iter().enumerate() {
                let prepared = solver.prepare(gf, left, right);
                if prepared.denominator != 0 {
                    self.regular.push(index);
                }
                self.prepared.push(prepared);
            }
            if !self.regular.is_empty() {
                let mut product = 1u128;
                for &index in &self.regular {
                    self.prefixes[index] = product;
                    product = gf.mul(product, self.prepared[index].denominator);
                }
                // A product of non-zero field elements is non-zero.
                let mut inverse_product = invert128(gf, basis, product);
                for &index in self.regular.iter().rev() {
                    self.inverses[index] = gf.mul(inverse_product, self.prefixes[index]);
                    inverse_product = gf.mul(inverse_product, self.prepared[index].denominator);
                }
            }
            for (prepared, &inverse) in self.prepared.iter().zip(&self.inverses) {
                self.results.push(solver.finish(gf, *prepared, inverse));
            }
        }
    }

    #[inline(always)]
    fn invert128(gf: &Gf2_128, basis: &NormalBasis128, a: u128) -> u128 {
        let e = gf.n - 1;
        let mut c = a;
        let mut k = 1u32;
        let bits = 32 - e.leading_zeros();
        for i in (0..bits - 1).rev() {
            let raised = basis
                .to_poly
                .apply(basis.rotate(basis.to_normal.apply(c), k));
            c = gf.mul(c, raised);
            k *= 2;
            if (e >> i) & 1 == 1 {
                c = gf.mul(gf.sqr(c), a);
                k += 1;
            }
        }
        gf.sqr(c)
    }

    /// Open-addressing `u128 -> u64` table; `u128::MAX` marks empty.
    /// Prefilter is the frozen wide example's two-word Bloom.
    struct RootBloom128 {
        words: Vec<u64>,
        bit_mask: usize,
    }

    impl RootBloom128 {
        const HASHES: usize = 3;

        fn with_capacity(keys: usize) -> Self {
            let bits = keys
                .checked_mul(8)
                .and_then(usize::checked_next_power_of_two)
                .expect("root Bloom size overflow")
                .max(64);
            Self {
                words: vec![0; bits / 64],
                bit_mask: bits - 1,
            }
        }

        #[inline(always)]
        fn hash_pair(key: u128) -> (u64, u64) {
            let low = key as u64;
            let high = (key >> 64) as u64;
            (
                mix_root_bloom(low ^ high.rotate_left(23) ^ 0xA5A5_38F9_9328_6F3D),
                mix_root_bloom(high ^ low.rotate_left(17) ^ 0xD4C9_3E72_B719_485B) | 1,
            )
        }

        #[inline(always)]
        fn insert(&mut self, key: u128) {
            let (first, step) = Self::hash_pair(key);
            for index in 0..Self::HASHES {
                let bit =
                    first.wrapping_add((index as u64).wrapping_mul(step)) as usize & self.bit_mask;
                self.words[bit >> 6] |= 1u64 << (bit & 63);
            }
        }

        #[inline(always)]
        fn contains(&self, key: u128) -> bool {
            let (first, step) = Self::hash_pair(key);
            for index in 0..Self::HASHES {
                let bit =
                    first.wrapping_add((index as u64).wrapping_mul(step)) as usize & self.bit_mask;
                if self.words[bit >> 6] & (1u64 << (bit & 63)) == 0 {
                    return false;
                }
            }
            true
        }

        fn bytes(&self) -> usize {
            self.words.len() * std::mem::size_of::<u64>()
        }
    }

    /// Open-addressing `u128 -> u64` table; `u128::MAX` marks empty.
    struct RootTable128 {
        slots: Vec<(u128, u64)>,
        shift: u32,
        len: usize,
        bloom: Option<RootBloom128>,
    }

    impl RootTable128 {
        fn with_capacity(entries: usize) -> Self {
            let slots = (entries * 2).next_power_of_two().max(16);
            Self {
                slots: vec![(u128::MAX, 0); slots],
                shift: 128 - slots.trailing_zeros(),
                len: 0,
                bloom: (std::env::var("KIC_ROOT_BLOOM").as_deref() == Ok("1"))
                    .then(|| RootBloom128::with_capacity(entries)),
            }
        }

        #[inline(always)]
        fn slot(&self, key: u128) -> usize {
            (key.wrapping_mul(0x9E37_79B9_7F4A_7C15_9E37_79B9_7F4A_7C15) >> self.shift) as usize
        }

        fn insert_if_absent(&mut self, key: u128, value: u64) {
            let mask = self.slots.len() - 1;
            let mut slot = self.slot(key);
            loop {
                if self.slots[slot].0 == u128::MAX {
                    self.slots[slot] = (key, value);
                    self.len += 1;
                    if let Some(bloom) = &mut self.bloom {
                        bloom.insert(key);
                    }
                    return;
                }
                if self.slots[slot].0 == key {
                    return;
                }
                slot = (slot + 1) & mask;
            }
        }

        #[inline(always)]
        fn get(&self, key: u128) -> Option<u64> {
            if self
                .bloom
                .as_ref()
                .is_some_and(|bloom| !bloom.contains(key))
            {
                return None;
            }
            let mask = self.slots.len() - 1;
            let mut slot = self.slot(key);
            loop {
                let (stored, value) = self.slots[slot];
                if stored == key {
                    return Some(value);
                }
                if stored == u128::MAX {
                    return None;
                }
                slot = (slot + 1) & mask;
            }
        }
    }

    struct State128 {
        left: u16,
        right: u16,
        relative: u16,
        normal_roots: [u128; 2],
    }

    struct Index128 {
        states: Vec<State128>,
        table: RootTable128,
        shifted: Vec<Vec<u128>>,
    }

    struct Relation128 {
        point_indices: [usize; 4],
        x_codes: [u128; 4],
        intermediates: [u128; 2],
        probes: u64,
    }

    #[derive(Clone, Copy)]
    enum TargetInput128 {
        PublicPoint(FastPoint128),
        KnownAnswerScalar(u64),
    }

    // serde_json's default Number representation cannot emit arbitrary
    // u128 field elements, so every wide field coordinate is a decimal
    // string on the wire.
    fn encode_point128(point: FastPoint128) -> Option<[String; 2]> {
        point.map(|(x, y)| [x.to_string(), y.to_string()])
    }

    fn encode_pair128(point: Option<[u128; 2]>) -> Option<[String; 2]> {
        point.map(|[x, y]| [x.to_string(), y.to_string()])
    }

    fn decode_word128(value: &Value) -> u128 {
        if let Some(decimal) = value.as_str() {
            decimal
                .parse()
                .expect("wide field coordinate must be decimal")
        } else if let Some(word) = value.as_u64() {
            word as u128
        } else {
            panic!("wide field coordinate must be a decimal string or u64")
        }
    }

    /// Parse a `[x, y]` target line directly from decimal text: JSON
    /// numbers beyond `u64` lose precision in `serde_json::Value`, so wide
    /// coordinates are read from the raw line (integers or quoted strings).
    fn parse_wide_pair_line(line: &str) -> [u128; 2] {
        let inner = line.trim().trim_start_matches('[').trim_end_matches(']');
        let mut words = inner.split(',');
        let parse = |raw: &str| -> u128 {
            raw.trim()
                .trim_matches('"')
                .parse()
                .expect("wide field coordinate must be a decimal integer or string")
        };
        let x = parse(words.next().expect("wide point needs two coordinates"));
        let y = parse(words.next().expect("wide point needs two coordinates"));
        assert!(
            words.next().is_none(),
            "wide point must have two coordinates"
        );
        [x, y]
    }

    fn decode_pair128(value: &Value) -> Option<[u128; 2]> {
        if value.is_null() {
            return None;
        }
        let pair = value.as_array().expect("wide point must be a [x,y] array");
        assert_eq!(pair.len(), 2, "wide point must have two coordinates");
        Some([decode_word128(&pair[0]), decode_word128(&pair[1])])
    }

    fn decode_point_list128(value: &Value) -> Vec<Option<[u128; 2]>> {
        value
            .as_array()
            .expect("wide point list must be an array")
            .iter()
            .map(decode_pair128)
            .collect()
    }

    struct Base128 {
        points: Vec<FastPoint128>,
        labels: Vec<(usize, u64)>,
        by_x: HashMap<u128, Vec<usize>>,
        columns: usize,
    }

    fn build_index128(
        gf: &Gf2_128,
        basis: &NormalBasis128,
        solver: &S3Solver128,
        reps: &[u128],
    ) -> Index128 {
        let n = gf.n as usize;
        let shifted: Vec<Vec<u128>> = reps
            .iter()
            .map(|&code| {
                let mut row = Vec::with_capacity(n);
                let mut value = code;
                for _ in 0..n {
                    row.push(value);
                    value = gf.sqr(value);
                }
                row
            })
            .collect();
        let mut states = Vec::new();
        for left in 0..reps.len() {
            for right in 0..reps.len() {
                for relative in 0..n {
                    if let Some(roots) =
                        solver.roots(gf, basis, shifted[left][0], shifted[right][relative])
                    {
                        states.push(State128 {
                            left: left as u16,
                            right: right as u16,
                            relative: relative as u16,
                            normal_roots: [
                                basis.to_normal.apply(roots[0]),
                                basis.to_normal.apply(roots[1]),
                            ],
                        });
                    }
                }
            }
        }
        let mut table = RootTable128::with_capacity(states.len() * 2);
        for state in &states {
            for &root in &state.normal_roots {
                let (canonical, shift) = basis.canonical_selected(root);
                table.insert_if_absent(
                    canonical,
                    pack(state.left, state.right, state.relative, shift),
                );
            }
        }
        Index128 {
            states,
            table,
            shifted,
        }
    }

    fn lift128(
        fast: &FastBinaryCurve128,
        base: &Base128,
        codes: &[u128; 4],
        target: FastPoint128,
    ) -> Option<[usize; 4]> {
        let choices: Vec<&Vec<usize>> = codes
            .iter()
            .map(|code| base.by_x.get(code))
            .collect::<Option<_>>()?;
        for &a in choices[0] {
            for &b in choices[1] {
                let ab = fast.add(base.points[a], base.points[b]);
                for &c in choices[2] {
                    let abc = fast.add(ab, base.points[c]);
                    for &d in choices[3] {
                        if fast.add(abc, base.points[d]) == target {
                            return Some([a, b, c, d]);
                        }
                    }
                }
            }
        }
        None
    }

    /// Per-state scratch for [`scan_state128`]: the `2n` candidate pairs
    /// and their batched root solutions.
    struct StateScan128 {
        inputs: Vec<(u128, u128)>,
        batch: RootBatch128,
    }

    impl StateScan128 {
        fn new(n: u32) -> Self {
            let entries = 2 * n as usize;
            Self {
                inputs: Vec::with_capacity(entries),
                batch: RootBatch128::with_capacity(entries),
            }
        }
    }

    /// Scan one index state against `target` in the sequential candidate
    /// order (shift outer, stored root inner) and return the first lifted
    /// relation.  Candidate order, probe count and published relation are
    /// exactly those of the frozen wide example; the laps around `lift128`
    /// feed the online PDP / relation-check split.
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn scan_state128(
        gf: &Gf2_128,
        fast: &FastBinaryCurve128,
        basis: &NormalBasis128,
        solver: &S3Solver128,
        index: &Index128,
        base: &Base128,
        target: FastPoint128,
        target_x: u128,
        state: &State128,
        scan: &mut StateScan128,
        probes: &mut u64,
        clock: &mut PhaseClock,
        split: &mut ExtractSplit,
    ) -> Option<Relation128> {
        let n = gf.n;
        scan.inputs.clear();
        for shift in 0..n {
            for &normal_root in &state.normal_roots {
                let absolute = basis.to_poly.apply(basis.rotate(normal_root, shift));
                scan.inputs.push((absolute, target_x));
            }
        }
        scan.batch.solve(gf, basis, solver, &scan.inputs);
        split.root_calls += scan.inputs.len() as u64;
        let mut candidate = 0usize;
        for shift in 0..n {
            let left_x = index.shifted[state.left as usize][shift as usize];
            let right_x = index.shifted[state.right as usize]
                [(shift as usize + state.relative as usize) % n as usize];
            for _ in 0..state.normal_roots.len() {
                let absolute = scan.inputs[candidate].0;
                let result = scan.batch.results[candidate];
                candidate += 1;
                let Some(partners) = result else {
                    continue;
                };
                for partner in partners {
                    *probes += 1;
                    let (canonical, partner_shift) =
                        basis.canonical(basis.to_normal.apply(partner));
                    let Some(value) = index.table.get(canonical) else {
                        continue;
                    };
                    let (left2, right2, relative2, stored_shift) = unpack(value);
                    let left_shift2 = ((stored_shift + n - partner_shift) % n) as usize;
                    let right_shift2 = (left_shift2 + relative2) % n as usize;
                    let codes = [
                        left_x,
                        right_x,
                        index.shifted[left2][left_shift2],
                        index.shifted[right2][right_shift2],
                    ];
                    split.pdp += clock.lap();
                    let lifted = lift128(fast, base, &codes, target);
                    split.relation_check += clock.lap();
                    split.relation_checks += 1;
                    if let Some(point_indices) = lifted {
                        return Some(Relation128 {
                            point_indices,
                            x_codes: codes,
                            intermediates: [absolute, partner],
                            probes: *probes,
                        });
                    }
                }
            }
        }
        None
    }

    /// The wide search with the same candidate order as the frozen wide
    /// example; laps are taken only around `lift128`, so the probe loop
    /// carries no timer.
    #[allow(clippy::too_many_arguments)]
    fn extract128_timed(
        gf: &Gf2_128,
        fast: &FastBinaryCurve128,
        basis: &NormalBasis128,
        solver: &S3Solver128,
        index: &Index128,
        base: &Base128,
        target: FastPoint128,
        start: usize,
        clock: &mut PhaseClock,
        split: &mut ExtractSplit,
    ) -> Option<Relation128> {
        let (target_x, _) = target?;
        let mut probes = 0u64;
        let mut scan = StateScan128::new(gf.n);
        let start = start % index.states.len().max(1);
        for state in index.states[start..].iter().chain(&index.states[..start]) {
            if let Some(relation) = scan_state128(
                gf,
                fast,
                basis,
                solver,
                index,
                base,
                target,
                target_x,
                state,
                &mut scan,
                &mut probes,
                clock,
                split,
            ) {
                return Some(relation);
            }
        }
        None
    }

    /// Exclusive PDP / relation-check time of one wide `extract` call.
    #[allow(clippy::too_many_arguments)]
    fn extract128(
        gf: &Gf2_128,
        fast: &FastBinaryCurve128,
        basis: &NormalBasis128,
        solver: &S3Solver128,
        index: &Index128,
        base: &Base128,
        target: FastPoint128,
        start: usize,
        clock: &mut PhaseClock,
        split: &mut ExtractSplit,
    ) -> Option<Relation128> {
        let relation = extract128_timed(
            gf, fast, basis, solver, index, base, target, start, clock, split,
        );
        split.pdp += clock.lap();
        relation
    }

    /// Frobenius-closed signed-orbit base from a deterministic public
    /// x-scan (wide twin: scanned abscissae are small integers, so the
    /// u64 scan counter converts losslessly).
    fn construct_base128(
        fast: &FastBinaryCurve128,
        curve: &KoblitzCurve,
        b: u128,
        columns: usize,
        r: u64,
    ) -> (
        Vec<FastPoint128>,
        Vec<(usize, u64)>,
        Vec<Option<[u128; 2]>>,
        u64,
    ) {
        let gf = &fast.gf;
        let n = fast.n as usize;
        let lambda = curve.lambda.to_u64().unwrap();
        let mut seen = std::collections::HashSet::new();
        let mut points = Vec::with_capacity(columns * 2 * n);
        let mut labels = Vec::with_capacity(columns * 2 * n);
        let mut reps = Vec::with_capacity(columns);
        let mut scanned = 0u64;
        let mut raw_x = 1u64;
        while reps.len() < columns {
            scanned += 1;
            for point in fast.points_with_x(b, raw_x as u128) {
                let Some((px, py)) = fast.scalar_mul(point, &curve.cofactor) else {
                    continue;
                };
                if px <= 1 {
                    continue;
                }
                let mut key = px;
                let mut x = px;
                for _ in 1..n {
                    x = gf.sqr(x);
                    key = key.min(x);
                }
                if !seen.insert(key) {
                    continue;
                }
                let column = reps.len();
                reps.push(Some([px, py]));
                let (mut cx, mut cy) = (px, py);
                let mut coefficient = 1u64;
                for _ in 0..n {
                    points.push(Some((cx, cy)));
                    labels.push((column, coefficient));
                    points.push(Some((cx, cx ^ cy)));
                    labels.push((column, (r - coefficient) % r));
                    cx = gf.sqr(cx);
                    cy = gf.sqr(cy);
                    coefficient = super::mulmod(coefficient, lambda, r);
                }
                if reps.len() == columns {
                    break;
                }
            }
            raw_x += 1;
        }
        let rep = reps[0].map(|[x, y]| (x, y));
        assert_eq!(
            fast.scalar_mul(rep, &BigUint::from(lambda)),
            rep.map(|(x, y)| (gf.sqr(x), gf.sqr(y))),
            "Frobenius must act as lambda on the subgroup"
        );
        (points, labels, reps, scanned)
    }

    fn relation_row128(base: &Base128, relation: &Relation128, rhs: u64, r: u64) -> Vec<u64> {
        let mut row = vec![0u64; base.columns + 1];
        for &index in &relation.point_indices {
            let (column, coefficient) = base.labels[index];
            row[column] = (row[column] + coefficient) % r;
        }
        row[base.columns] = rhs % r;
        row
    }

    /// Wide online pipeline for `64 < n ≤ 127`: same stages, record
    /// fields and JSON shapes as the frozen wide flow, plus the online
    /// interval instrumentation of the `u64` flow above (five phase
    /// timers on one clock, hardware counters, relation-check counting).
    pub fn run(arguments: &[String]) -> ! {
        assert_eq!(arguments.len(), 5);
        let process_started = Instant::now();
        let constructed: Option<(u32, u8, usize)> =
            arguments[1].strip_prefix("construct:").map(|spec| {
                let parts: Vec<&str> = spec.split(':').collect();
                assert_eq!(parts.len(), 3, "construct:<n>:<a>:<columns>");
                (
                    parts[0].parse().unwrap(),
                    parts[1].parse().unwrap(),
                    parts[2].parse().unwrap(),
                )
            });
        let header: Value = match constructed {
            Some((n, a, _)) => json!({"n":n, "a":a}),
            None => {
                let header_bytes = std::fs::read(&arguments[1]).expect("read base header");
                let header_line = header_bytes.split(|&b| b == b'\n').next().unwrap();
                let header: Value = serde_json::from_slice(header_line).expect("parse base header");
                assert_eq!(header["kind"], "point_defined_factor_base");
                header
            }
        };
        let n = header["n"].as_u64().unwrap() as u32;
        let a = header["a"].as_u64().unwrap() as u8;
        assert!((64..=127).contains(&n), "wide pipeline needs 64 < n ≤ 127");
        let target_inputs: Vec<TargetInput128> = std::fs::read_to_string(&arguments[2])
            .expect("read target points or legacy fixture scalars")
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                let line = line.trim();
                if line.starts_with('[') {
                    let [x, y] = parse_wide_pair_line(line);
                    TargetInput128::PublicPoint(Some((x, y)))
                } else {
                    TargetInput128::KnownAnswerScalar(line.parse().unwrap())
                }
            })
            .collect();
        let mut rank_seed: u64 = arguments[3].parse().unwrap();
        let mut out = std::fs::File::create(&arguments[4]).expect("create output");

        let setup_started = Instant::now();
        let curve = KoblitzCurve::new(a, n).expect("admitted Koblitz rung");
        let r = curve.subgroup_order.to_u64().unwrap();
        if constructed.is_none() {
            assert_eq!(header["subgroup_order"].as_u64(), Some(r));
            let low_terms: Vec<u64> =
                serde_json::from_value(header["field_modulus_low_terms"].clone()).unwrap();
            assert_eq!(
                low_terms,
                curve
                    .curve
                    .irreducible
                    .low_terms
                    .iter()
                    .map(|&t| t as u64)
                    .collect::<Vec<_>>(),
                "base header field must match the constructed curve"
            );
        }
        let fast = FastBinaryCurve128::new(&curve.curve.irreducible, a as u128).unwrap();
        let gf = Gf2_128::new(&curve.curve.irreducible);
        let b = gf.from_element(&curve.curve.b);
        let generator: FastPoint128 = match curve.generator() {
            crypto_lib::binary_ecc::BinaryPoint::Affine { x, y } => {
                Some((fast.word(x), fast.word(y)))
            }
            crypto_lib::binary_ecc::BinaryPoint::Infinity => panic!("generator must be affine"),
        };
        let mut scanned_x = None;
        let (points, labels, representatives): (
            Vec<FastPoint128>,
            Vec<(usize, u64)>,
            Vec<Option<[u128; 2]>>,
        ) = match constructed {
            Some((_, _, columns)) => {
                let (points, labels, reps, scanned) =
                    construct_base128(&fast, &curve, b, columns, r);
                scanned_x = Some(scanned);
                (points, labels, reps)
            }
            None => {
                let coordinates = decode_point_list128(&header["factor_base_point_coordinates"]);
                let representatives = decode_point_list128(&header["factor_base_representatives"]);
                (
                    coordinates.iter().map(|c| c.map(|[x, y]| (x, y))).collect(),
                    serde_json::from_value(header["factor_base_point_labels"].clone()).unwrap(),
                    representatives,
                )
            }
        };
        let base_hash = match constructed {
            Some(_) => {
                let mut hasher = blake3::Hasher::new();
                hasher.update(b"compact-orbit-constructed-base-v1");
                for (point, &(column, coefficient)) in points.iter().zip(&labels) {
                    let (x, y): (u128, u128) = point.unwrap();
                    for word in [x, y] {
                        hasher.update(&word.to_le_bytes());
                    }
                    hasher.update(&(column as u64).to_le_bytes());
                    hasher.update(&coefficient.to_le_bytes());
                }
                hasher.finalize().to_hex().to_string()
            }
            None => header["base_hash"].as_str().unwrap().to_owned(),
        };
        let mut by_x: HashMap<u128, Vec<usize>> = HashMap::new();
        for (index, point) in points.iter().enumerate() {
            if let Some((x, _)) = point {
                by_x.entry(*x).or_default().push(index);
            }
        }
        let base = Base128 {
            points,
            labels,
            by_x,
            columns: representatives.len(),
        };
        let reps: Vec<u128> = representatives
            .iter()
            .map(|p| p.expect("affine representative")[0])
            .collect();
        let base_elapsed = setup_started.elapsed();
        let base_load_ms = ms(base_elapsed);

        let basis_started = Instant::now();
        let basis = NormalBasis128::new(&gf);
        let solver = S3Solver128::new(&gf, b);
        for probe in 1..2048u128 {
            let x = probe.wrapping_mul(0x2545_F491_4F6C_DD1D) & basis.mask;
            assert_eq!(
                basis
                    .to_poly
                    .apply(basis.rotate(basis.to_normal.apply(x), 1)),
                gf.sqr(x)
            );
            let y = (probe * 0x9E37_79B9) & basis.mask;
            assert_eq!(
                solver.roots(&gf, &basis, x, y),
                s3_x_roots_128(&gf, b, x, y)
            );
        }
        let basis_elapsed = basis_started.elapsed();
        let basis_ms = ms(basis_elapsed);

        let index_started = Instant::now();
        let index = build_index128(&gf, &basis, &solver, &reps);
        let index_elapsed = index_started.elapsed();
        let index_ms = ms(index_elapsed);
        let setup_ms = ms(setup_started.elapsed());

        // Rank stage: rows sum(coefficient * L_column) = a for random [a]G,
        // with the same query / PDP / relation-check split as the `u64`
        // online flow so the cold-phase accounting stays complete.
        let rank_started = Instant::now();
        let mut echelon = Echelon {
            r,
            columns: base.columns,
            pivots: vec![None; base.columns],
            rank: 0,
        };
        let mut rank_attempts = 0u64;
        let mut rank_failures = 0u64;
        let mut rank_relations = 0u64;
        let mut rank_probes = 0u64;
        let mut rank_root_calls = 0u64;
        let mut rank_rows_without_gain = 0u64;
        let mut rank_query_elapsed = std::time::Duration::ZERO;
        let mut rank_pdp_elapsed = std::time::Duration::ZERO;
        let mut rank_relation_check_elapsed = std::time::Duration::ZERO;
        // Guided: decompose [a]G - R_j for a column j that has no pivot yet.
        while echelon.rank < base.columns {
            let query_started = Instant::now();
            rank_seed = rank_seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let scalar = (rank_seed >> 11) % (r - 1) + 1;
            let column = (0..base.columns)
                .find(|&c| echelon.pivots[c].is_none())
                .unwrap();
            let rep: FastPoint128 = representatives[column].map(|[x, y]| (x, y));
            let point = fast.add(
                fast.scalar_mul(generator, &BigUint::from(scalar)),
                FastBinaryCurve128::neg(rep),
            );
            rank_attempts += 1;
            let query_elapsed = query_started.elapsed();
            rank_query_elapsed += query_elapsed;
            let mut rank_clock = PhaseClock::start();
            let mut rank_split = ExtractSplit::default();
            let relation = extract128(
                &gf,
                &fast,
                &basis,
                &solver,
                &index,
                &base,
                point,
                (rank_seed >> 20) as usize,
                &mut rank_clock,
                &mut rank_split,
            );
            rank_pdp_elapsed += rank_split.pdp;
            rank_relation_check_elapsed += rank_split.relation_check;
            rank_root_calls += rank_split.root_calls;
            match relation {
                Some(relation) => {
                    rank_relations += 1;
                    rank_probes += relation.probes;
                    let mut row = relation_row128(&base, &relation, scalar, r);
                    row[column] = (row[column] + 1) % r;
                    if !echelon.insert(row) {
                        rank_rows_without_gain += 1;
                    }
                }
                None => rank_failures += 1,
            }
        }
        let rank_elapsed = rank_started.elapsed();
        let rank_ms = ms(rank_elapsed);
        let rank_matrix_elapsed = rank_elapsed
            .checked_sub(rank_query_elapsed + rank_pdp_elapsed + rank_relation_check_elapsed)
            .expect("rank phase accounting exceeds rank interval");
        let la_started = Instant::now();
        let logs = echelon.solve();
        let la_elapsed = la_started.elapsed();
        let la_ms = ms(la_elapsed);
        let setup_complete_elapsed = process_started.elapsed();
        let setup_complete_ns = setup_complete_elapsed.as_nanos() as u64;
        let input_control_elapsed = setup_complete_elapsed
            .checked_sub(base_elapsed + basis_elapsed + index_elapsed + rank_elapsed + la_elapsed)
            .expect("setup phase accounting exceeds setup interval");

        let targets_started = Instant::now();
        let mut solved = 0usize;
        let mut failed = 0usize;
        let mut target_query_ms = Vec::with_capacity(target_inputs.len());
        let mut single_target_phases = None;
        for (fixture_index, &target_input) in target_inputs.iter().enumerate() {
            // Fixture construction is outside the online interval: a real
            // query arrives as a point Q, not as its known-answer scalar.
            let (target, published, target_generation_ms) = match target_input {
                TargetInput128::PublicPoint(target) => (target, None, 0.0),
                TargetInput128::KnownAnswerScalar(scalar) => {
                    let target_generation_started = Instant::now();
                    let target = fast.scalar_mul(generator, &BigUint::from(scalar));
                    let elapsed = target_generation_started.elapsed().as_secs_f64() * 1000.0;
                    (target, Some(scalar), elapsed)
                }
            };
            let online_counters = HwCounters::open();
            online_counters.start();
            let mut clock = PhaseClock::start();
            let online_started = clock.last;
            // The scan origin must be a function of the public target point
            // only; the wide hash matches the frozen wide producer.
            let (target_x, target_y) = target.expect("published fixture target is not infinity");
            let mut query_hash =
                target_x ^ target_y.rotate_left(29) ^ 0x9E37_79B9_7F4A_7C15_9E37_79B9_7F4A_7C15;
            query_hash = (query_hash ^ (query_hash >> 30))
                .wrapping_mul(0xBF58_476D_1CE4_E5B9_BF58_476D_1CE4_E5B9);
            query_hash = (query_hash ^ (query_hash >> 27))
                .wrapping_mul(0x94D0_49BB_1331_11EB_94D0_49BB_1331_11EB);
            query_hash ^= query_hash >> 31;
            let start = query_hash as usize;
            let query_phase = clock.lap();
            let mut split = ExtractSplit::default();
            let relation = extract128(
                &gf, &fast, &basis, &solver, &index, &base, target, start, &mut clock, &mut split,
            );
            let recovered = relation.as_ref().map(|relation| {
                relation.point_indices.iter().fold(0u64, |acc, &index| {
                    let (column, coefficient) = base.labels[index];
                    (acc + super::mulmod(coefficient, logs[column], r)) % r
                })
            });
            let descent_phase = clock.lap();
            let verified =
                recovered.map(|d| fast.scalar_mul(generator, &BigUint::from(d)) == target);
            let recovery_check_phase = clock.lap();
            let online_stopped = clock.last;
            let online_hw_counts = online_counters.stop();
            let online = online_stopped - online_started;
            if target_inputs.len() == 1 {
                single_target_phases = Some((
                    query_phase,
                    split.pdp,
                    split.relation_check,
                    descent_phase,
                    recovery_check_phase,
                ));
            }
            let elapsed = ms(online);
            target_query_ms.push(elapsed);
            if verified == Some(true) {
                solved += 1;
            } else {
                failed += 1;
            }
            let target_exit_code = if verified == Some(true) { 0 } else { 1 };
            let record = json!({
                "kind":"compact_orbit_dlp_target",
                "n":n, "a":a, "fixture_index":fixture_index,
                "published_fixture_scalar":published,
                "field_value_encoding":"base-10 strings for u128 field values",
                "generator":encode_point128(generator),
                "target":encode_point128(target),
                "published_q":encode_point128(target),
                "exit_code":target_exit_code,
                "x_codes":relation.as_ref().map(|relation| relation.x_codes.map(|x| x.to_string())),
                "pinned_intermediates":relation.as_ref().map(|relation| relation.intermediates.map(|x| x.to_string())),
                "point_indices":relation.as_ref().map(|relation| relation.point_indices),
                "probes":relation.as_ref().map(|relation| relation.probes),
                "recovered_scalar":recovered,
                "recovered_matches_published":published.and_then(|scalar| recovered.map(|d| d == scalar % r)),
                "group_verified":verified,
                "relation_checks":split.relation_checks,
                "root_calls":split.root_calls,
                "skipped_symmetric_states":split.skipped_symmetric_states,
                "target_generation_ms_excluded":target_generation_ms,
                "online_start_event":"target_query_begin",
                "online_stop_event":"recovery_check_end",
                "online_start_ns":(online_started - process_started).as_nanos() as u64,
                "online_stop_ns":(online_stopped - process_started).as_nanos() as u64,
                "online_ms":elapsed,
                "online_hw_counts":online_hw_counts,
                "online_hw_scope":"calling thread, user space, target query through recovery check",
                "target_query_ms":ms(query_phase),
                "target_pdp_ms":ms(split.pdp),
                "target_relation_check_ms":ms(split.relation_check),
                "target_descent_ms":ms(descent_phase),
                "target_recovery_check_ms":ms(recovery_check_phase),
                "target_phase_sum_ms":ms(query_phase + split.pdp + split.relation_check + descent_phase + recovery_check_phase),
                "target_ms":elapsed,
            });
            writeln!(out, "{record}").unwrap();
        }
        let targets_ms = targets_started.elapsed().as_secs_f64() * 1000.0;
        let total_ms = process_started.elapsed().as_secs_f64() * 1000.0;
        let mut sorted = target_query_ms.clone();
        sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let median = sorted.get(sorted.len() / 2).copied();
        let (cold_in_process_ms, cold_phase_ms) = if let Some((
            target_query,
            target_pdp,
            target_relation_check,
            target_descent,
            target_recovery_check,
        )) = single_target_phases
        {
            let phases = [
                input_control_elapsed,
                base_elapsed,
                basis_elapsed,
                index_elapsed,
                rank_query_elapsed,
                rank_pdp_elapsed,
                rank_relation_check_elapsed,
                rank_matrix_elapsed,
                la_elapsed,
                target_query,
                target_pdp,
                target_relation_check,
                target_descent,
                target_recovery_check,
            ];
            let phase_sum: std::time::Duration = phases.iter().copied().sum();
            let cold_total = setup_complete_elapsed
                + target_query
                + target_pdp
                + target_relation_check
                + target_descent
                + target_recovery_check;
            assert_eq!(phase_sum, cold_total, "exclusive cold phases differ");
            (
                Some(ms(cold_total)),
                Some(json!({
                    "setup_input_control":ms(input_control_elapsed),
                    "isogeny":0.0,
                    "factor_base":ms(base_elapsed),
                    "precompute_basis_and_selftest":ms(basis_elapsed),
                    "precompute_index":ms(index_elapsed),
                    "rank_query":ms(rank_query_elapsed),
                    "rank_pdp":ms(rank_pdp_elapsed),
                    "rank_relation_check":ms(rank_relation_check_elapsed),
                    "rank_matrix_build":ms(rank_matrix_elapsed),
                    "relation_la":ms(la_elapsed),
                    "target_query":ms(target_query),
                    "target_pdp":ms(target_pdp),
                    "target_relation_check":ms(target_relation_check),
                    "target_descent":ms(target_descent),
                    "target_recovery_check":ms(target_recovery_check),
                })),
            )
        } else {
            (None, None)
        };
        let summary = json!({
            "kind":"compact_orbit_dlp_summary",
            "schema_version":"1.0",
            "producer_version":"compact_orbit_online_pair_symmetry_v1",
            "setup_complete_ns":setup_complete_ns,
            "cold_in_process_ms":cold_in_process_ms,
            "cold_phase_ms":cold_phase_ms,
            "n":n, "a":a,
            "base_hash":base_hash,
            "base_source":if constructed.is_some() { "constructed_in_process_x_scan" } else { "retained_header" },
            "base_scanned_x":scanned_x,
            "orbit_columns":base.columns,
            "factor_base_points":base.points.len(),
            "regular_states":index.states.len(),
            "root_table_entries":index.table.len,
            "root_bloom_enabled":index.table.bloom.is_some(),
            "root_bloom_bytes":index.table.bloom.as_ref().map_or(0, RootBloom128::bytes),
            "sparse_canonical_enabled":basis.sparse_canonical,
            "pair_table_entries":0,
            "edge_selectors":0,
            "timing_ms":{
                "base_load":base_load_ms,
                "normal_basis_and_selftest":basis_ms,
                "index_build":index_ms,
                "setup_total":setup_ms,
                "rank_stage":rank_ms,
                "linear_algebra":la_ms,
                "targets_total":targets_ms,
                "target_median":median,
                "process_total":total_ms,
            },
            "rank_attempts":rank_attempts,
            "rank_relations":rank_relations,
            "rank_failures":rank_failures,
            "rank_probes_total":rank_probes,
            "rank_root_calls_total":rank_root_calls,
            "rank_rows_without_gain":rank_rows_without_gain,
            "rank_policy":"guided: decompose [a]G - R_j for the first pivotless column j",
            "rank_probes_mean":if rank_relations > 0 { rank_probes as f64 / rank_relations as f64 } else { 0.0 },
            "rank":echelon.rank,
            "targets":target_inputs.len(),
            "targets_solved":solved,
            "targets_failed":failed,
            "peak_rss_bytes":super::peak_rss_bytes(),
            "threads":1,
            "target_threads":1,
            "scope":"public synthetic Koblitz fixtures; shared-factor-log DLP with every stage timed in one process; no external points or key recovery",
        });
        println!("{summary}");
        if let Ok(path) = std::env::var("KIC_DUMP_BASE") {
            let dump = json!({
                "kind":"point_defined_factor_base",
                "n":n, "a":a,
                "subgroup_order":r,
                "field_modulus_low_terms":curve.curve.irreducible.low_terms,
                "orbit_columns":base.columns,
                "base_hash":base_hash,
                "factor_base_points":base.points.len(),
                "field_value_encoding":"base-10 strings for u128 field values",
                "factor_base_point_coordinates":base.points.iter().map(|p| encode_point128(*p)).collect::<Vec<_>>(),
                "factor_base_point_labels":base.labels,
                "factor_base_representatives":representatives.iter().copied().map(encode_pair128).collect::<Vec<_>>(),
            });
            std::fs::write(path, format!("{dump}\n")).expect("write base dump");
        }
        std::process::exit(if failed == 0 { 0 } else { 1 });
    }
}

fn main() -> ExitCode {
    let arguments: Vec<String> = std::env::args().collect();
    assert_eq!(
        arguments.len(),
        5,
        "usage: <base_header.jsonl | construct:<n>:<a>:<columns>> <target_points_or_legacy_scalars.txt> <rank_seed> <out.jsonl>"
    );
    // Wide pipeline (`u128` words) for 64 < n ≤ 127.  The degree is
    // detected here so the `u64` flow below stays byte-identical; the wide
    // path re-parses its own arguments.
    {
        let wide_n: Option<u32> = match arguments[1].strip_prefix("construct:") {
            Some(spec) => spec.split(':').next().and_then(|n| n.parse().ok()),
            None => {
                let header_bytes = std::fs::read(&arguments[1]).expect("read base header");
                let header_line = header_bytes.split(|&b| b == b'\n').next().unwrap();
                let header: Value = serde_json::from_slice(header_line).expect("parse base header");
                header["n"].as_u64().map(|n| n as u32)
            }
        };
        if wide_n.is_some_and(|n| n > 63) {
            wide::run(&arguments);
        }
    }
    let process_started = Instant::now();
    let constructed: Option<(u32, u8, usize)> =
        arguments[1].strip_prefix("construct:").map(|spec| {
            let parts: Vec<&str> = spec.split(':').collect();
            assert_eq!(parts.len(), 3, "construct:<n>:<a>:<columns>");
            (
                parts[0].parse().unwrap(),
                parts[1].parse().unwrap(),
                parts[2].parse().unwrap(),
            )
        });
    let header: Value = match constructed {
        Some((n, a, _)) => json!({"n":n, "a":a}),
        None => {
            let header_bytes = std::fs::read(&arguments[1]).expect("read base header");
            let header_line = header_bytes.split(|&b| b == b'\n').next().unwrap();
            let header: Value = serde_json::from_slice(header_line).expect("parse base header");
            assert_eq!(header["kind"], "point_defined_factor_base");
            header
        }
    };
    let n = header["n"].as_u64().unwrap() as u32;
    let a = header["a"].as_u64().unwrap() as u8;
    let target_inputs: Vec<TargetInput> = std::fs::read_to_string(&arguments[2])
        .expect("read target points or legacy fixture scalars")
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let line = line.trim();
            if line.starts_with('[') {
                let [x, y]: [u64; 2] =
                    serde_json::from_str(line).expect("target point must be [x,y]");
                TargetInput::PublicPoint(Some((x, y)))
            } else {
                TargetInput::KnownAnswerScalar(line.parse().unwrap())
            }
        })
        .collect();
    let mut rank_seed: u64 = arguments[3].parse().unwrap();
    let rank_probe_cap = std::env::var("KIC_RANK_PROBE_CAP").ok().map(|raw| {
        let cap: u64 = raw
            .parse()
            .expect("KIC_RANK_PROBE_CAP must be a positive integer");
        assert!(cap > 0, "KIC_RANK_PROBE_CAP must be positive");
        cap
    });
    let pair_symmetry = match std::env::var("KIC_PAIR_SYMMETRY").as_deref() {
        Ok("1") => true,
        Ok("0") | Err(_) => false,
        Ok(_) => panic!("KIC_PAIR_SYMMETRY must be 0 or 1"),
    };
    let mut out = std::fs::File::create(&arguments[4]).expect("create output");

    let setup_started = Instant::now();
    let curve = KoblitzCurve::new(a, n).expect("admitted Koblitz rung");
    let r = curve.subgroup_order.to_u64().unwrap();
    if constructed.is_none() {
        assert_eq!(header["subgroup_order"].as_u64(), Some(r));
        let low_terms: Vec<u64> =
            serde_json::from_value(header["field_modulus_low_terms"].clone()).unwrap();
        assert_eq!(
            low_terms,
            curve
                .curve
                .irreducible
                .low_terms
                .iter()
                .map(|&t| t as u64)
                .collect::<Vec<_>>(),
            "base header field must match the constructed curve"
        );
    }
    let fast = FastBinaryCurve::new(&curve.curve.irreducible, a as u64).unwrap();
    let gf = Gf2::new(&curve.curve.irreducible);
    let b = gf.from_element(&curve.curve.b);
    let generator: FastPoint = match curve.generator() {
        crypto_lib::binary_ecc::BinaryPoint::Affine { x, y } => Some((fast.word(x), fast.word(y))),
        crypto_lib::binary_ecc::BinaryPoint::Infinity => panic!("generator must be affine"),
    };
    let mut scanned_x = None;
    let (points, labels, representatives): (
        Vec<FastPoint>,
        Vec<(usize, u64)>,
        Vec<Option<[u64; 2]>>,
    ) = match constructed {
        Some((_, _, columns)) => {
            let (points, labels, reps, scanned) = construct_base(&fast, &curve, b, columns, r);
            scanned_x = Some(scanned);
            (points, labels, reps)
        }
        None => {
            let coordinates: Vec<Option<[u64; 2]>> =
                serde_json::from_value(header["factor_base_point_coordinates"].clone()).unwrap();
            (
                coordinates.iter().map(|c| c.map(|[x, y]| (x, y))).collect(),
                serde_json::from_value(header["factor_base_point_labels"].clone()).unwrap(),
                serde_json::from_value(header["factor_base_representatives"].clone()).unwrap(),
            )
        }
    };
    let base_hash = match constructed {
        Some(_) => {
            let mut hasher = blake3::Hasher::new();
            hasher.update(b"compact-orbit-constructed-base-v1");
            for (point, &(column, coefficient)) in points.iter().zip(&labels) {
                let (x, y) = point.unwrap();
                for word in [x, y, column as u64, coefficient] {
                    hasher.update(&word.to_le_bytes());
                }
            }
            hasher.finalize().to_hex().to_string()
        }
        None => header["base_hash"].as_str().unwrap().to_owned(),
    };
    let mut by_x: HashMap<u64, Vec<usize>> = HashMap::new();
    for (index, point) in points.iter().enumerate() {
        if let Some((x, _)) = point {
            by_x.entry(*x).or_default().push(index);
        }
    }
    let base = Base {
        points,
        labels,
        by_x,
        columns: representatives.len(),
    };
    let reps: Vec<u64> = representatives
        .iter()
        .map(|p| p.expect("affine representative")[0])
        .collect();
    let base_elapsed = setup_started.elapsed();
    let base_load_ms = ms(base_elapsed);

    let basis_started = Instant::now();
    let basis = NormalBasis::new(&gf);
    let solver = S3Solver::new(&gf, b);
    for probe in 1..2048u64 {
        let x = probe.wrapping_mul(0x2545_F491_4F6C_DD1D) & basis.mask;
        assert_eq!(
            basis
                .to_poly
                .apply(basis.rotate(basis.to_normal.apply(x), 1)),
            gf.sqr(x)
        );
        let y = (probe * 0x9E37_79B9) & basis.mask;
        assert_eq!(solver.roots(&gf, &basis, x, y), s3_x_roots(&gf, b, x, y));
    }
    let basis_elapsed = basis_started.elapsed();
    let basis_ms = ms(basis_elapsed);

    let index_started = Instant::now();
    let index = build_index(&gf, &basis, &solver, &reps, pair_symmetry);
    let index_elapsed = index_started.elapsed();
    let index_ms = ms(index_elapsed);
    let setup_ms = ms(setup_started.elapsed());

    // Rank stage: rows sum(coefficient * L_column) = a for random [a]G.
    // The optional trace lets a separate implementation replay every attempted
    // group relation, rank transition, and final base-log solution.
    let rank_started = Instant::now();
    let mut rank_trace = std::env::var("KIC_DUMP_RANK")
        .ok()
        .map(|path| std::fs::File::create(path).expect("create rank trace"));
    if let Some(trace) = rank_trace.as_mut() {
        writeln!(
            trace,
            "{}",
            json!({
                "kind":"compact_orbit_rank_header", "schema_version":1,
                "n":n, "a":a, "subgroup_order":r, "base_hash":base_hash,
                "rank_probe_cap":rank_probe_cap,
                "pair_symmetry":pair_symmetry,
                "generator":generator.map(|(x, y)| [x, y]),
                "orbit_columns":base.columns, "factor_base_points":base.points.len(),
            })
        )
        .unwrap();
    }
    let mut echelon = Echelon {
        r,
        columns: base.columns,
        pivots: vec![None; base.columns],
        rank: 0,
    };
    let mut rank_attempts = 0u64;
    let mut rank_failures = 0u64;
    let mut rank_relations = 0u64;
    let mut rank_probes = 0u64;
    let mut rank_root_calls = 0u64;
    let mut rank_skipped_symmetric_states = 0u64;
    let mut rank_relation_probes = 0u64;
    let mut rank_capped_attempts = 0u64;
    let mut rank_capped_probes = 0u64;
    let mut rank_query_elapsed = Duration::ZERO;
    let mut rank_pdp_elapsed = Duration::ZERO;
    let mut rank_relation_check_elapsed = Duration::ZERO;
    // Guided: decompose [a]G - R_j for a column j that has no pivot yet, so the
    // row a = L_j + sum(...) always contains column j and raises the rank.
    let mut rank_rows_without_gain = 0u64;
    while echelon.rank < base.columns {
        let query_started = Instant::now();
        rank_seed = rank_seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let scalar = (rank_seed >> 11) % (r - 1) + 1;
        let column = (0..base.columns)
            .find(|&c| echelon.pivots[c].is_none())
            .unwrap();
        let rep: FastPoint = representatives[column].map(|[x, y]| (x, y));
        let point = fast.add(
            fast.scalar_mul(generator, &BigUint::from(scalar)),
            FastBinaryCurve::neg(rep),
        );
        rank_attempts += 1;
        let rank_before = echelon.rank;
        let query_elapsed = query_started.elapsed();
        rank_query_elapsed += query_elapsed;
        let mut rank_clock = PhaseClock::start();
        let mut rank_split = ExtractSplit::default();
        let outcome = extract(
            &gf,
            &fast,
            &basis,
            &solver,
            &index,
            &base,
            point,
            (rank_seed >> 20) as usize,
            rank_probe_cap,
            &mut rank_clock,
            &mut rank_split,
        );
        rank_pdp_elapsed += rank_split.pdp;
        rank_relation_check_elapsed += rank_split.relation_check;
        rank_root_calls += rank_split.root_calls;
        rank_skipped_symmetric_states += rank_split.skipped_symmetric_states;
        let ExtractOutcome {
            relation,
            probes,
            capped,
        } = outcome;
        rank_probes += probes;
        match relation {
            Some(relation) => {
                rank_relations += 1;
                rank_relation_probes += probes;
                let mut row = relation_row(&base, &relation, scalar, r);
                row[column] = (row[column] + 1) % r;
                let gained = echelon.insert(row.clone());
                if !gained {
                    rank_rows_without_gain += 1;
                }
                if let Some(trace) = rank_trace.as_mut() {
                    writeln!(
                        trace,
                        "{}",
                        json!({
                            "kind":"compact_orbit_rank_attempt", "attempt_index":rank_attempts - 1,
                            "found":true, "scalar":scalar, "pivotless_column":column,
                            "target":point.map(|(x, y)| [x, y]),
                            "point_indices":relation.point_indices, "x_codes":relation.x_codes,
                            "pinned_intermediates":relation.intermediates,
                            "probes":probes, "capped":false, "row":row,
                            "query_ms":ms(query_elapsed), "pdp_ms":ms(rank_split.pdp),
                            "relation_check_ms":ms(rank_split.relation_check),
                            "relation_checks":rank_split.relation_checks,
                            "root_calls":rank_split.root_calls,
                            "skipped_symmetric_states":rank_split.skipped_symmetric_states,
                            "rank_before":rank_before, "rank_after":echelon.rank, "gained":gained,
                        })
                    )
                    .unwrap();
                }
            }
            None => {
                rank_failures += 1;
                if capped {
                    rank_capped_attempts += 1;
                    rank_capped_probes += probes;
                }
                if let Some(trace) = rank_trace.as_mut() {
                    writeln!(
                        trace,
                        "{}",
                        json!({
                            "kind":"compact_orbit_rank_attempt", "attempt_index":rank_attempts - 1,
                            "found":false, "scalar":scalar, "pivotless_column":column,
                            "target":point.map(|(x, y)| [x, y]),
                            "probes":probes, "capped":capped,
                            "query_ms":ms(query_elapsed), "pdp_ms":ms(rank_split.pdp),
                            "relation_check_ms":ms(rank_split.relation_check),
                            "relation_checks":rank_split.relation_checks,
                            "root_calls":rank_split.root_calls,
                            "skipped_symmetric_states":rank_split.skipped_symmetric_states,
                            "rank_before":rank_before, "rank_after":echelon.rank,
                        })
                    )
                    .unwrap();
                }
            }
        }
    }
    let rank_elapsed = rank_started.elapsed();
    let rank_ms = ms(rank_elapsed);
    let rank_matrix_elapsed = rank_elapsed
        .checked_sub(rank_query_elapsed + rank_pdp_elapsed + rank_relation_check_elapsed)
        .expect("rank phase accounting exceeds rank interval");
    let la_started = Instant::now();
    let logs = echelon.solve();
    if let Some(trace) = rank_trace.as_mut() {
        writeln!(
            trace,
            "{}",
            json!({
                "kind":"compact_orbit_rank_solution", "rank":echelon.rank,
                "attempts":rank_attempts, "relations":rank_relations,
                "failures":rank_failures, "rows_without_gain":rank_rows_without_gain,
                "rank_probe_cap":rank_probe_cap, "probes_total":rank_probes,
                "capped_attempts":rank_capped_attempts, "capped_probes":rank_capped_probes,
                "logs":logs,
            })
        )
        .unwrap();
        trace.flush().unwrap();
    }
    let la_elapsed = la_started.elapsed();
    let la_ms = ms(la_elapsed);
    let setup_complete_elapsed = process_started.elapsed();
    let setup_complete_ns = setup_complete_elapsed.as_nanos() as u64;
    let input_control_elapsed = setup_complete_elapsed
        .checked_sub(base_elapsed + basis_elapsed + index_elapsed + rank_elapsed + la_elapsed)
        .expect("setup phase accounting exceeds setup interval");

    let targets_started = Instant::now();
    let mut solved = 0usize;
    let mut failed = 0usize;
    let mut target_query_ms = Vec::with_capacity(target_inputs.len());
    let mut single_target_phases = None;
    for (fixture_index, &target_input) in target_inputs.iter().enumerate() {
        // Fixture construction is outside the online interval: a real query
        // arrives as a point Q, not as its known-answer scalar.
        let (target, published, target_generation_ms) = match target_input {
            TargetInput::PublicPoint(target) => (target, None, 0.0),
            TargetInput::KnownAnswerScalar(scalar) => {
                let target_generation_started = Instant::now();
                let target = fast.scalar_mul(generator, &BigUint::from(scalar));
                let elapsed = target_generation_started.elapsed().as_secs_f64() * 1000.0;
                (target, Some(scalar), elapsed)
            }
        };
        let online_counters = HwCounters::open();
        online_counters.start();
        let mut clock = PhaseClock::start();
        let online_started = clock.last;
        // The scan origin must be a function of the public target point only;
        // using the fixture's known-answer scalar here would leak the DLP.
        let (target_x, target_y) = target.expect("published fixture target is not infinity");
        let mut query_hash = target_x ^ target_y.rotate_left(29) ^ 0x9E37_79B9_7F4A_7C15;
        query_hash = (query_hash ^ (query_hash >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        query_hash = (query_hash ^ (query_hash >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        query_hash ^= query_hash >> 31;
        let start = query_hash as usize;
        let query_phase = clock.lap();
        let mut split = ExtractSplit::default();
        let relation = extract(
            &gf, &fast, &basis, &solver, &index, &base, target, start, None, &mut clock, &mut split,
        )
        .relation;
        let recovered = relation.as_ref().map(|relation| {
            relation.point_indices.iter().fold(0u64, |acc, &index| {
                let (column, coefficient) = base.labels[index];
                (acc + mulmod(coefficient, logs[column], r)) % r
            })
        });
        let descent_phase = clock.lap();
        let verified = recovered.map(|d| fast.scalar_mul(generator, &BigUint::from(d)) == target);
        let recovery_check_phase = clock.lap();
        let online_stopped = clock.last;
        let online_hw_counts = online_counters.stop();
        let online = online_stopped - online_started;
        if target_inputs.len() == 1 {
            single_target_phases = Some((
                query_phase,
                split.pdp,
                split.relation_check,
                descent_phase,
                recovery_check_phase,
            ));
        }
        let elapsed = ms(online);
        target_query_ms.push(elapsed);
        if verified == Some(true) {
            solved += 1;
        } else {
            failed += 1;
        }
        let target_exit_code = if verified == Some(true) { 0 } else { 1 };
        let record = json!({
            "kind":"compact_orbit_dlp_target",
            "n":n, "a":a, "fixture_index":fixture_index,
            "published_fixture_scalar":published,
            "generator":generator.map(|(x, y)| [x, y]),
            "target":target.map(|(x, y)| [x, y]),
            "published_q":target.map(|(x, y)| [x, y]),
            "exit_code":target_exit_code,
            "x_codes":relation.as_ref().map(|relation| relation.x_codes),
            "pinned_intermediates":relation.as_ref().map(|relation| relation.intermediates),
            "point_indices":relation.as_ref().map(|relation| relation.point_indices),
            "probes":relation.as_ref().map(|relation| relation.probes),
            "recovered_scalar":recovered,
            "recovered_matches_published":published.and_then(|scalar| recovered.map(|d| d == scalar % r)),
            "group_verified":verified,
            "relation_checks":split.relation_checks,
            "root_calls":split.root_calls,
            "skipped_symmetric_states":split.skipped_symmetric_states,
            "target_generation_ms_excluded":target_generation_ms,
            "online_start_event":"target_query_begin",
            "online_stop_event":"recovery_check_end",
            "online_start_ns":(online_started - process_started).as_nanos() as u64,
            "online_stop_ns":(online_stopped - process_started).as_nanos() as u64,
            "online_ms":elapsed,
            "online_hw_counts":online_hw_counts,
            "online_hw_scope":"calling thread, user space, target query through recovery check",
            "target_query_ms":ms(query_phase),
            "target_pdp_ms":ms(split.pdp),
            "target_relation_check_ms":ms(split.relation_check),
            "target_descent_ms":ms(descent_phase),
            "target_recovery_check_ms":ms(recovery_check_phase),
            "target_phase_sum_ms":ms(query_phase + split.pdp + split.relation_check + descent_phase + recovery_check_phase),
            "target_ms":elapsed,
        });
        writeln!(out, "{record}").unwrap();
    }
    let targets_ms = targets_started.elapsed().as_secs_f64() * 1000.0;
    let total_ms = process_started.elapsed().as_secs_f64() * 1000.0;
    let mut sorted = target_query_ms.clone();
    sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let median = sorted.get(sorted.len() / 2).copied();
    let (cold_in_process_ms, cold_phase_ms) = if let Some((
        target_query,
        target_pdp,
        target_relation_check,
        target_descent,
        target_recovery_check,
    )) = single_target_phases
    {
        let phases = [
            input_control_elapsed,
            base_elapsed,
            basis_elapsed,
            index_elapsed,
            rank_query_elapsed,
            rank_pdp_elapsed,
            rank_relation_check_elapsed,
            rank_matrix_elapsed,
            la_elapsed,
            target_query,
            target_pdp,
            target_relation_check,
            target_descent,
            target_recovery_check,
        ];
        let phase_sum: Duration = phases.iter().copied().sum();
        let cold_total = setup_complete_elapsed
            + target_query
            + target_pdp
            + target_relation_check
            + target_descent
            + target_recovery_check;
        assert_eq!(phase_sum, cold_total, "exclusive cold phases differ");
        (
            Some(ms(cold_total)),
            Some(json!({
                "setup_input_control":ms(input_control_elapsed),
                "isogeny":0.0,
                "factor_base":ms(base_elapsed),
                "precompute_basis_and_selftest":ms(basis_elapsed),
                "precompute_index":ms(index_elapsed),
                "rank_query":ms(rank_query_elapsed),
                "rank_pdp":ms(rank_pdp_elapsed),
                "rank_relation_check":ms(rank_relation_check_elapsed),
                "rank_matrix_build":ms(rank_matrix_elapsed),
                "relation_la":ms(la_elapsed),
                "target_query":ms(target_query),
                "target_pdp":ms(target_pdp),
                "target_relation_check":ms(target_relation_check),
                "target_descent":ms(target_descent),
                "target_recovery_check":ms(target_recovery_check),
            })),
        )
    } else {
        (None, None)
    };
    let summary = json!({
        "kind":"compact_orbit_dlp_summary",
        "schema_version":"1.0",
        "producer_version":"compact_orbit_online_pair_symmetry_v1",
        "setup_complete_ns":setup_complete_ns,
        "cold_in_process_ms":cold_in_process_ms,
        "cold_phase_ms":cold_phase_ms,
        "n":n, "a":a,
        "base_hash":base_hash,
        "base_source":if constructed.is_some() { "constructed_in_process_x_scan" } else { "retained_header" },
        "base_scanned_x":scanned_x,
        "orbit_columns":base.columns,
        "factor_base_points":base.points.len(),
        "regular_states":index.states.len(),
        "index_pair_candidates":index.pair_candidates,
        "pair_symmetry":pair_symmetry,
        "root_table_entries":index.table.len,
        "root_bloom_enabled":index.table.bloom.is_some(),
        "root_bloom_bytes":index.table.bloom.as_ref().map_or(0, RootBloom::bytes),
        "pair_table_entries":0,
        "edge_selectors":0,
        "timing_ms":{
            "base_load":base_load_ms,
            "normal_basis_and_selftest":basis_ms,
            "index_build":index_ms,
            "setup_total":setup_ms,
            "rank_stage":rank_ms,
            "linear_algebra":la_ms,
            "targets_total":targets_ms,
            "target_median":median,
            "process_total":total_ms,
        },
        "rank_attempts":rank_attempts,
        "rank_relations":rank_relations,
        "rank_failures":rank_failures,
        "rank_probe_cap":rank_probe_cap,
        "rank_capped_attempts":rank_capped_attempts,
        "rank_capped_probes":rank_capped_probes,
        "rank_exhausted_attempts":rank_failures - rank_capped_attempts,
        "rank_probes_total":rank_probes,
        "rank_root_calls_total":rank_root_calls,
        "rank_skipped_symmetric_states":rank_skipped_symmetric_states,
        "rank_rows_without_gain":rank_rows_without_gain,
        "rank_policy":"guided: decompose [a]G - R_j for the first pivotless column j",
        "rank_restart_policy":if rank_probe_cap.is_some() { "abandon after cap support probes; draw next scalar" } else { "unbounded" },
        "rank_probes_mean":if rank_relations > 0 { rank_relation_probes as f64 / rank_relations as f64 } else { 0.0 },
        "rank_probes_mean_per_attempt":if rank_attempts > 0 { rank_probes as f64 / rank_attempts as f64 } else { 0.0 },
        "rank":echelon.rank,
        "targets":target_inputs.len(),
        "targets_solved":solved,
        "targets_failed":failed,
        "peak_rss_bytes":peak_rss_bytes(),
        "threads":1,
        "scope":"public synthetic Koblitz fixtures; shared-factor-log DLP with every stage timed in one process; no external points or key recovery",
    });
    println!("{summary}");
    if let Ok(path) = std::env::var("KIC_DUMP_BASE") {
        let dump = json!({
            "kind":"point_defined_factor_base",
            "n":n, "a":a,
            "subgroup_order":r,
            "field_modulus_low_terms":curve.curve.irreducible.low_terms,
            "orbit_columns":base.columns,
            "base_hash":base_hash,
            "factor_base_points":base.points.len(),
            "factor_base_point_coordinates":base.points.iter().map(|p| p.map(|(x, y)| [x, y])).collect::<Vec<_>>(),
            "factor_base_point_labels":base.labels,
            "factor_base_representatives":representatives,
        });
        std::fs::write(path, format!("{dump}\n")).expect("write base dump");
    }
    if failed == 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

#[cfg(test)]
mod pair_symmetry_tests {
    use super::*;

    #[test]
    fn unordered_pair_orbits_preserve_every_index_root() {
        let curve = KoblitzCurve::new(0, 13).unwrap();
        let gf = Gf2::new(&curve.curve.irreducible);
        let basis = NormalBasis::new(&gf);
        let solver = S3Solver::new(&gf, gf.from_element(&curve.curve.b));
        let reps = [2, 3, 5, 7, 11];
        let full = build_index(&gf, &basis, &solver, &reps, false);
        let folded = build_index(&gf, &basis, &solver, &reps, true);
        let keys = |index: &Index| {
            let mut roots: Vec<u64> = index
                .table
                .slots
                .iter()
                .filter_map(|&(key, _)| (key != u64::MAX).then_some(key))
                .collect();
            roots.sort_unstable();
            roots
        };
        assert_eq!(
            full.pair_candidates,
            reps.len() * reps.len() * gf.n as usize
        );
        assert_eq!(folded.pair_candidates, 165);
        let full_keys = keys(&full);
        assert_eq!(full_keys, keys(&folded));
        assert_eq!(folded.states.len(), full.states.len());
        for key in full_keys {
            assert_eq!(full.table.get(key), folded.table.get(key));
        }
        for (index, &mate) in folded.counterparts.as_ref().unwrap().iter().enumerate() {
            assert_eq!(
                folded.counterparts.as_ref().unwrap()[mate as usize],
                index as u32
            );
        }
    }

    #[test]
    fn n53_pair_quotient_preserves_the_complete_root_table() {
        let curve = KoblitzCurve::new(0, 53).unwrap();
        let fast = FastBinaryCurve::new(&curve.curve.irreducible, 0).unwrap();
        let gf = Gf2::new(&curve.curve.irreducible);
        let b = gf.from_element(&curve.curve.b);
        let (_, _, representatives, _) = construct_base(
            &fast,
            &curve,
            b,
            220,
            curve.subgroup_order.to_u64().unwrap(),
        );
        let reps: Vec<u64> = representatives
            .iter()
            .map(|point| point.unwrap()[0])
            .collect();
        let basis = NormalBasis::new(&gf);
        let solver = S3Solver::new(&gf, b);
        let full = build_index(&gf, &basis, &solver, &reps, false);
        let folded = build_index(&gf, &basis, &solver, &reps, true);
        assert_eq!(full.pair_candidates, 2_565_200);
        assert_eq!(folded.pair_candidates, 1_282_710);
        assert_eq!(full.states.len(), folded.states.len());
        assert_eq!(full.table.len, folded.table.len);
        for &(key, value) in &full.table.slots {
            if key != u64::MAX {
                assert_eq!(folded.table.get(key), Some(value));
            }
        }
    }
}
