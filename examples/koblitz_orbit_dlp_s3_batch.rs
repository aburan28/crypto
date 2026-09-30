//! Experimental swap/Frobenius-quotiented compact-orbit shared-log DLP.
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
//! Usage: <base_header.jsonl> <target_points.jsonl|legacy_scalars.txt> <rank_seed> <out.jsonl>

use crypto_lib::cryptanalysis::koblitz_fast_arith::{s3_x_roots, FastBinaryCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

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

    #[inline(always)]
    fn prepare(&self, gf: &Gf2, left: u64, right: u64) -> PreparedRoot {
        let a = self.square.apply(left ^ right);
        let p = gf.mul(left, right);
        let ps = self.square.apply(p);
        let denominator = if a == 0 || ps == 0 { 0 } else { gf.mul(a, ps) };
        PreparedRoot {
            left,
            right,
            a,
            p,
            ps,
            denominator,
        }
    }

    #[inline(always)]
    fn finish(&self, gf: &Gf2, prepared: PreparedRoot, inverse: u64) -> Option<[u64; 2]> {
        if prepared.denominator == 0 {
            return s3_x_roots(gf, self.b, prepared.left, prepared.right);
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

    /// Scalar control, with the same root order as the frozen swap producer.
    #[inline(always)]
    fn roots(&self, gf: &Gf2, basis: &NormalBasis, left: u64, right: u64) -> Option<[u64; 2]> {
        let prepared = self.prepare(gf, left, right);
        let inverse = if prepared.denominator == 0 {
            0
        } else {
            invert(gf, basis, prepared.denominator)
        };
        self.finish(gf, prepared, inverse)
    }

    #[inline(always)]
    fn roots_counted(
        &self,
        gf: &Gf2,
        basis: &NormalBasis,
        left: u64,
        right: u64,
        counts: &mut S3Counts,
    ) -> Option<[u64; 2]> {
        let prepared = self.prepare(gf, left, right);
        counts.calls += 1;
        let inverse = if prepared.denominator == 0 {
            counts.exceptional_calls += 1;
            0
        } else {
            counts.regular_calls += 1;
            counts.scalar_inversions += 1;
            invert(gf, basis, prepared.denominator)
        };
        let result = self.finish(gf, prepared, inverse);
        if result.is_some() {
            counts.returned_pairs += 1;
        } else {
            counts.no_root_returns += 1;
        }
        counts.consumed_candidates += 1;
        result
    }
}

#[derive(Clone, Copy)]
struct PreparedRoot {
    left: u64,
    right: u64,
    a: u64,
    p: u64,
    ps: u64,
    denominator: u64,
}

#[derive(Clone, Copy, Default)]
struct S3Counts {
    calls: u64,
    regular_calls: u64,
    exceptional_calls: u64,
    returned_pairs: u64,
    no_root_returns: u64,
    scalar_inversions: u64,
    batch_inversions: u64,
    batch_multiplications: u64,
    consumed_candidates: u64,
    discarded_candidates: u64,
    table_lookups: u64,
    table_lookup_probes: u64,
    table_insert_probes: u64,
    lift_attempts: u64,
    root_keys_considered: u64,
    prefilter_checks: u64,
    prefilter_definite_misses: u64,
    prefilter_false_positives: u64,
    prefilter_true_positives: u64,
}

impl S3Counts {
    fn validate(&self, window: usize) {
        assert_eq!(self.calls, self.regular_calls + self.exceptional_calls);
        assert_eq!(self.calls, self.returned_pairs + self.no_root_returns);
        assert_eq!(
            self.calls,
            self.consumed_candidates + self.discarded_candidates
        );
        assert!(self.table_lookup_probes >= self.table_lookups);
        assert_eq!(
            self.root_keys_considered,
            self.table_lookups + self.prefilter_definite_misses
        );
        assert_eq!(
            self.prefilter_checks,
            self.prefilter_definite_misses
                + self.prefilter_false_positives
                + self.prefilter_true_positives
        );
        assert!(self.prefilter_checks == 0 || self.prefilter_checks == self.root_keys_considered);
        if window == 1 {
            assert_eq!(self.scalar_inversions, self.regular_calls);
            assert_eq!(self.batch_inversions, 0);
            assert_eq!(self.discarded_candidates, 0);
        } else {
            assert_eq!(self.scalar_inversions, 0);
            assert_eq!(self.batch_multiplications, self.regular_calls * 3);
        }
    }

    fn add(&mut self, other: &Self) {
        self.calls += other.calls;
        self.regular_calls += other.regular_calls;
        self.exceptional_calls += other.exceptional_calls;
        self.returned_pairs += other.returned_pairs;
        self.no_root_returns += other.no_root_returns;
        self.scalar_inversions += other.scalar_inversions;
        self.batch_inversions += other.batch_inversions;
        self.batch_multiplications += other.batch_multiplications;
        self.consumed_candidates += other.consumed_candidates;
        self.discarded_candidates += other.discarded_candidates;
        self.table_lookups += other.table_lookups;
        self.table_lookup_probes += other.table_lookup_probes;
        self.table_insert_probes += other.table_insert_probes;
        self.lift_attempts += other.lift_attempts;
        self.root_keys_considered += other.root_keys_considered;
        self.prefilter_checks += other.prefilter_checks;
        self.prefilter_definite_misses += other.prefilter_definite_misses;
        self.prefilter_false_positives += other.prefilter_false_positives;
        self.prefilter_true_positives += other.prefilter_true_positives;
    }

    fn as_json(&self) -> Value {
        json!({
            "calls":self.calls, "regular_calls":self.regular_calls,
            "exceptional_calls":self.exceptional_calls,
            "returned_pairs":self.returned_pairs,
            "no_root_returns":self.no_root_returns,
            "scalar_inversions":self.scalar_inversions,
            "batch_inversions":self.batch_inversions,
            "batch_multiplications":self.batch_multiplications,
            "consumed_candidates":self.consumed_candidates,
            "discarded_candidates":self.discarded_candidates,
            "table_lookups":self.table_lookups,
            "table_lookup_probes":self.table_lookup_probes,
            "table_insert_probes":self.table_insert_probes,
            "lift_attempts":self.lift_attempts,
            "root_keys_considered":self.root_keys_considered,
            "prefilter_checks":self.prefilter_checks,
            "prefilter_definite_misses":self.prefilter_definite_misses,
            "prefilter_false_positives":self.prefilter_false_positives,
            "prefilter_true_positives":self.prefilter_true_positives,
        })
    }
}

#[derive(Default)]
struct RootBatchScratch {
    prepared: Vec<PreparedRoot>,
    inverses: Vec<u64>,
    prefixes: Vec<u64>,
    regular_indices: Vec<usize>,
    results: Vec<Option<[u64; 2]>>,
}

impl RootBatchScratch {
    fn solve(
        &mut self,
        gf: &Gf2,
        basis: &NormalBasis,
        solver: &S3Solver,
        inputs: &[(u64, u64)],
        counts: &mut S3Counts,
    ) {
        self.prepared.clear();
        self.inverses.clear();
        self.prefixes.clear();
        self.regular_indices.clear();
        self.results.clear();
        self.inverses.resize(inputs.len(), 0);
        self.prefixes.resize(inputs.len(), 0);
        for (index, &(left, right)) in inputs.iter().enumerate() {
            let prepared = solver.prepare(gf, left, right);
            counts.calls += 1;
            if prepared.denominator == 0 {
                counts.exceptional_calls += 1;
            } else {
                counts.regular_calls += 1;
                self.regular_indices.push(index);
            }
            self.prepared.push(prepared);
        }
        if !self.regular_indices.is_empty() {
            let mut product = 1;
            for &index in &self.regular_indices {
                self.prefixes[index] = product;
                product = gf.mul(product, self.prepared[index].denominator);
                counts.batch_multiplications += 1;
            }
            let mut inverse_product = invert(gf, basis, product);
            counts.batch_inversions += 1;
            for &index in self.regular_indices.iter().rev() {
                self.inverses[index] = gf.mul(inverse_product, self.prefixes[index]);
                inverse_product = gf.mul(inverse_product, self.prepared[index].denominator);
                counts.batch_multiplications += 2;
            }
        }
        for (prepared, &inverse) in self.prepared.iter().zip(&self.inverses) {
            let result = solver.finish(gf, *prepared, inverse);
            if result.is_some() {
                counts.returned_pairs += 1;
            } else {
                counts.no_root_returns += 1;
            }
            self.results.push(result);
        }
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
struct RootTable {
    slots: Vec<(u64, u64)>,
    shift: u32,
    len: usize,
}

impl RootTable {
    fn with_capacity(entries: usize) -> Self {
        let slots = (entries * 2).next_power_of_two().max(16);
        Self {
            slots: vec![(u64::MAX, 0); slots],
            shift: 64 - slots.trailing_zeros(),
            len: 0,
        }
    }

    #[inline(always)]
    fn slot(&self, key: u64) -> usize {
        (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> self.shift) as usize
    }

    fn insert_if_absent(&mut self, key: u64, value: u64) -> u64 {
        let mask = self.slots.len() - 1;
        let mut slot = self.slot(key);
        let mut probes = 0;
        loop {
            probes += 1;
            if self.slots[slot].0 == u64::MAX {
                self.slots[slot] = (key, value);
                self.len += 1;
                return probes;
            }
            if self.slots[slot].0 == key {
                return probes;
            }
            slot = (slot + 1) & mask;
        }
    }

    #[cfg(test)]
    #[inline(always)]
    fn get(&self, key: u64) -> Option<u64> {
        self.get_counted(key).0
    }

    #[inline(always)]
    fn get_counted(&self, key: u64) -> (Option<u64>, u64) {
        let mask = self.slots.len() - 1;
        let mut slot = self.slot(key);
        let mut probes = 0;
        loop {
            probes += 1;
            let (stored, value) = self.slots[slot];
            if stored == key {
                return (Some(value), probes);
            }
            if stored == u64::MAX {
                return (None, probes);
            }
            slot = (slot + 1) & mask;
        }
    }
}

/// Three hashes within one 64-byte block. Every inserted key sets all three
/// bits; a definite miss is therefore safe to skip, regardless of collisions.
struct BlockedBloom {
    blocks: Vec<[u64; 8]>,
}

impl BlockedBloom {
    fn with_capacity(expected_keys: usize) -> Self {
        let blocks = expected_keys
            .saturating_mul(12)
            .div_ceil(512)
            .next_power_of_two()
            .max(1);
        Self {
            blocks: vec![[0; 8]; blocks],
        }
    }

    #[inline(always)]
    fn positions(&self, key: u64) -> (usize, [usize; 3]) {
        let mut h = key.wrapping_add(0x9E37_79B9_7F4A_7C15);
        h = (h ^ (h >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        h = (h ^ (h >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        h ^= h >> 31;
        (
            (h as usize) & (self.blocks.len() - 1),
            [
                ((h >> 19) & 511) as usize,
                ((h >> 28) & 511) as usize,
                ((h >> 37) & 511) as usize,
            ],
        )
    }

    #[inline(always)]
    fn insert(&mut self, key: u64) {
        let (block, positions) = self.positions(key);
        for position in positions {
            self.blocks[block][position >> 6] |= 1u64 << (position & 63);
        }
    }

    #[inline(always)]
    fn may_contain(&self, key: u64) -> bool {
        let (block, positions) = self.positions(key);
        positions
            .iter()
            .all(|&position| self.blocks[block][position >> 6] & (1u64 << (position & 63)) != 0)
    }

    fn bytes(&self) -> usize {
        self.blocks.len() * 64
    }
}

struct State {
    left: u16,
    right: u16,
    relative: u16,
    normal_roots: [u64; 2],
}

struct Index {
    states: Vec<State>,
    table: RootTable,
    prefilter: Option<BlockedBloom>,
    shifted: Vec<Vec<u64>>,
    representative_candidates: usize,
}

/// Choose one state from the involution (l,r,t) <-> (r,l,-t mod n).
/// Its only fixed states at odd n are (l,l,0).
#[inline(always)]
fn retain_swap_representative(left: usize, right: usize, relative: usize, n: usize) -> bool {
    (left, right, relative) <= (right, left, (n - relative) % n)
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

fn append_index_batch(
    gf: &Gf2,
    basis: &NormalBasis,
    solver: &S3Solver,
    inputs: &[(u64, u64)],
    metadata: &[(usize, usize, usize)],
    scratch: &mut RootBatchScratch,
    states: &mut Vec<State>,
    counts: &mut S3Counts,
) {
    scratch.solve(gf, basis, solver, inputs, counts);
    counts.consumed_candidates += inputs.len() as u64;
    for (&(left, right, relative), roots) in metadata.iter().zip(&scratch.results) {
        if let Some(roots) = roots {
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

fn build_index(
    gf: &Gf2,
    basis: &NormalBasis,
    solver: &S3Solver,
    reps: &[u64],
    window: usize,
    use_prefilter: bool,
    counts: &mut S3Counts,
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
    let mut states = Vec::new();
    let mut representative_candidates = 0usize;
    let mut inputs = Vec::with_capacity(window);
    let mut metadata = Vec::with_capacity(window);
    let mut scratch = RootBatchScratch::default();
    for left in 0..reps.len() {
        for right in 0..reps.len() {
            for relative in 0..n {
                if !retain_swap_representative(left, right, relative, n) {
                    continue;
                }
                representative_candidates += 1;
                let pair = (shifted[left][0], shifted[right][relative]);
                if window == 1 {
                    if let Some(roots) = solver.roots_counted(gf, basis, pair.0, pair.1, counts) {
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
                } else {
                    inputs.push(pair);
                    metadata.push((left, right, relative));
                    if inputs.len() == window {
                        append_index_batch(
                            gf,
                            basis,
                            solver,
                            &inputs,
                            &metadata,
                            &mut scratch,
                            &mut states,
                            counts,
                        );
                        inputs.clear();
                        metadata.clear();
                    }
                }
            }
        }
    }
    if !inputs.is_empty() {
        append_index_batch(
            gf,
            basis,
            solver,
            &inputs,
            &metadata,
            &mut scratch,
            &mut states,
            counts,
        );
    }
    let mut table = RootTable::with_capacity(states.len() * 2);
    let mut prefilter = use_prefilter.then(|| BlockedBloom::with_capacity(states.len() * 2));
    for state in &states {
        for &root in &state.normal_roots {
            let (canonical, shift) = basis.canonical(root);
            if let Some(filter) = prefilter.as_mut() {
                filter.insert(canonical);
            }
            counts.table_insert_probes += table.insert_if_absent(
                canonical,
                pack(state.left, state.right, state.relative, shift),
            );
        }
    }
    Index {
        states,
        table,
        prefilter,
        shifted,
        representative_candidates,
    }
}

struct Relation {
    point_indices: [usize; 4],
    x_codes: [u64; 4],
    intermediates: [u64; 2],
    probes: u64,
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

#[derive(Clone, Copy)]
struct QueryCandidate {
    left_x: u64,
    right_x: u64,
    absolute: u64,
}

fn consume_query_batch(
    gf: &Gf2,
    fast: &FastBinaryCurve,
    basis: &NormalBasis,
    solver: &S3Solver,
    index: &Index,
    base: &Base,
    target: FastPoint,
    window: usize,
    inputs: &[(u64, u64)],
    candidates: &[QueryCandidate],
    scratch: &mut RootBatchScratch,
    counts: &mut S3Counts,
    probes: &mut u64,
) -> Option<Relation> {
    if window > 1 {
        scratch.solve(gf, basis, solver, inputs, counts);
    }
    let n = gf.n;
    for (candidate_index, &candidate) in candidates.iter().enumerate() {
        let roots = if window == 1 {
            solver.roots_counted(
                gf,
                basis,
                inputs[candidate_index].0,
                inputs[candidate_index].1,
                counts,
            )
        } else {
            counts.consumed_candidates += 1;
            scratch.results[candidate_index]
        };
        let Some(partners) = roots else {
            continue;
        };
        for partner in partners {
            *probes += 1;
            let (canonical, partner_shift) = basis.canonical(basis.to_normal.apply(partner));
            counts.root_keys_considered += 1;
            if let Some(filter) = &index.prefilter {
                counts.prefilter_checks += 1;
                if !filter.may_contain(canonical) {
                    counts.prefilter_definite_misses += 1;
                    continue;
                }
            }
            let (value, slot_probes) = index.table.get_counted(canonical);
            counts.table_lookups += 1;
            counts.table_lookup_probes += slot_probes;
            if index.prefilter.is_some() {
                if value.is_some() {
                    counts.prefilter_true_positives += 1;
                } else {
                    counts.prefilter_false_positives += 1;
                }
            }
            let Some(value) = value else {
                continue;
            };
            let (left2, right2, relative2, stored_shift) = unpack(value);
            let left_shift2 = ((stored_shift + n - partner_shift) % n) as usize;
            let right_shift2 = (left_shift2 + relative2) % n as usize;
            let codes = [
                candidate.left_x,
                candidate.right_x,
                index.shifted[left2][left_shift2],
                index.shifted[right2][right_shift2],
            ];
            counts.lift_attempts += 1;
            if let Some(point_indices) = lift(fast, base, &codes, target) {
                if window > 1 {
                    counts.discarded_candidates += (candidates.len() - candidate_index - 1) as u64;
                }
                return Some(Relation {
                    point_indices,
                    x_codes: codes,
                    intermediates: [candidate.absolute, partner],
                    probes: *probes,
                });
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
    window: usize,
    counts: &mut S3Counts,
) -> Option<Relation> {
    let (target_x, _) = target?;
    let n = gf.n;
    let mut probes = 0u64;
    let mut inputs = Vec::with_capacity(window);
    let mut candidates = Vec::with_capacity(window);
    let mut scratch = RootBatchScratch::default();
    // Keep the original state/shift/root order. A full batch is prepared
    // before consumption, and unused prefetched candidates remain charged.
    let start = start % index.states.len().max(1);
    for state in index.states[start..].iter().chain(&index.states[..start]) {
        for shift in 0..n {
            let left_x = index.shifted[state.left as usize][shift as usize];
            let right_x = index.shifted[state.right as usize]
                [(shift as usize + state.relative as usize) % n as usize];
            for &normal_root in &state.normal_roots {
                let absolute = basis.to_poly.apply(basis.rotate(normal_root, shift));
                inputs.push((absolute, target_x));
                candidates.push(QueryCandidate {
                    left_x,
                    right_x,
                    absolute,
                });
                if inputs.len() == window {
                    if let Some(relation) = consume_query_batch(
                        gf,
                        fast,
                        basis,
                        solver,
                        index,
                        base,
                        target,
                        window,
                        &inputs,
                        &candidates,
                        &mut scratch,
                        counts,
                        &mut probes,
                    ) {
                        return Some(relation);
                    }
                    inputs.clear();
                    candidates.clear();
                }
            }
        }
    }
    if !inputs.is_empty() {
        return consume_query_batch(
            gf,
            fast,
            basis,
            solver,
            index,
            base,
            target,
            window,
            &inputs,
            &candidates,
            &mut scratch,
            counts,
            &mut probes,
        );
    }
    None
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
    let output = std::process::Command::new("ps")
        .args(["-o", "rss=", "-p", &std::process::id().to_string()])
        .output()
        .ok()?;
    let kib: u64 = String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse()
        .ok()?;
    Some(kib * 1024)
}

fn main() {
    let arguments: Vec<String> = std::env::args().collect();
    assert_eq!(
        arguments.len(),
        5,
        "usage: <base_header.jsonl | construct:<n>:<a>:<columns>> <target_points_or_legacy_scalars.txt> <rank_seed> <out.jsonl>"
    );
    let process_started = Instant::now();
    let window: usize = std::env::var("KIC_S3_BATCH_WINDOW")
        .unwrap_or_else(|_| "1".to_string())
        .parse()
        .expect("batch window integer");
    assert!(
        matches!(window, 1 | 16 | 64),
        "batch window must be 1, 16, or 64"
    );
    let use_prefilter = match std::env::var("KIC_S3_PREFILTER").as_deref() {
        Ok("blocked") => true,
        Ok("off") | Err(_) => false,
        _ => panic!("KIC_S3_PREFILTER must be off or blocked"),
    };
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
    let base_load_ms = setup_started.elapsed().as_secs_f64() * 1000.0;

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
    let basis_ms = basis_started.elapsed().as_secs_f64() * 1000.0;

    let index_started = Instant::now();
    let mut index_counts = S3Counts::default();
    let index = build_index(
        &gf,
        &basis,
        &solver,
        &reps,
        window,
        use_prefilter,
        &mut index_counts,
    );
    index_counts.validate(window);
    let index_ms = index_started.elapsed().as_secs_f64() * 1000.0;
    let setup_ms = setup_started.elapsed().as_secs_f64() * 1000.0;

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
    let mut rank_counts = S3Counts::default();
    // Guided: decompose [a]G - R_j for a column j that has no pivot yet, so the
    // row a = L_j + sum(...) always contains column j and raises the rank.
    let mut rank_rows_without_gain = 0u64;
    while echelon.rank < base.columns {
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
        match extract(
            &gf,
            &fast,
            &basis,
            &solver,
            &index,
            &base,
            point,
            (rank_seed >> 20) as usize,
            window,
            &mut rank_counts,
        ) {
            Some(relation) => {
                rank_relations += 1;
                rank_probes += relation.probes;
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
                            "probes":relation.probes, "row":row,
                            "rank_before":rank_before, "rank_after":echelon.rank, "gained":gained,
                        })
                    )
                    .unwrap();
                }
            }
            None => {
                rank_failures += 1;
                if let Some(trace) = rank_trace.as_mut() {
                    writeln!(
                        trace,
                        "{}",
                        json!({
                            "kind":"compact_orbit_rank_attempt", "attempt_index":rank_attempts - 1,
                            "found":false, "scalar":scalar, "pivotless_column":column,
                            "target":point.map(|(x, y)| [x, y]),
                            "rank_before":rank_before, "rank_after":echelon.rank,
                        })
                    )
                    .unwrap();
                }
            }
        }
    }
    rank_counts.validate(window);
    let rank_ms = rank_started.elapsed().as_secs_f64() * 1000.0;
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
                "logs":logs,
            })
        )
        .unwrap();
        trace.flush().unwrap();
    }
    let la_ms = la_started.elapsed().as_secs_f64() * 1000.0;

    let targets_started = Instant::now();
    let mut solved = 0usize;
    let mut failed = 0usize;
    let mut target_query_ms = Vec::with_capacity(target_inputs.len());
    let mut target_counts = S3Counts::default();
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
        let query_started = Instant::now();
        // The scan origin must be a function of the public target point only;
        // using the fixture's known-answer scalar here would leak the DLP.
        let query_hash_started = Instant::now();
        let (target_x, target_y) = target.expect("published fixture target is not infinity");
        let mut query_hash = target_x ^ target_y.rotate_left(29) ^ 0x9E37_79B9_7F4A_7C15;
        query_hash = (query_hash ^ (query_hash >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        query_hash = (query_hash ^ (query_hash >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        query_hash ^= query_hash >> 31;
        let start = query_hash as usize;
        let target_query_stage_ms = query_hash_started.elapsed().as_secs_f64() * 1000.0;
        let target_pdp_started = Instant::now();
        let mut query_counts = S3Counts::default();
        let relation = extract(
            &gf,
            &fast,
            &basis,
            &solver,
            &index,
            &base,
            target,
            start,
            window,
            &mut query_counts,
        );
        query_counts.validate(window);
        target_counts.add(&query_counts);
        let target_pdp_and_relation_check_ms = target_pdp_started.elapsed().as_secs_f64() * 1000.0;
        let target_descent_started = Instant::now();
        let recovered = relation.as_ref().map(|relation| {
            relation.point_indices.iter().fold(0u64, |acc, &index| {
                let (column, coefficient) = base.labels[index];
                (acc + mulmod(coefficient, logs[column], r)) % r
            })
        });
        let target_descent_ms = target_descent_started.elapsed().as_secs_f64() * 1000.0;
        let target_recovery_check_started = Instant::now();
        let verified = recovered.map(|d| fast.scalar_mul(generator, &BigUint::from(d)) == target);
        let target_recovery_check_ms =
            target_recovery_check_started.elapsed().as_secs_f64() * 1000.0;
        let elapsed = query_started.elapsed().as_secs_f64() * 1000.0;
        target_query_ms.push(elapsed);
        if verified == Some(true) {
            solved += 1;
        } else {
            failed += 1;
        }
        let record = json!({
            "kind":"compact_orbit_dlp_target",
            "n":n, "a":a, "fixture_index":fixture_index,
            "published_fixture_scalar":published,
            "generator":generator.map(|(x, y)| [x, y]),
            "target":target.map(|(x, y)| [x, y]),
            "published_q":target.map(|(x, y)| [x, y]),
            "exit_code":0,
            "x_codes":relation.as_ref().map(|relation| relation.x_codes),
            "pinned_intermediates":relation.as_ref().map(|relation| relation.intermediates),
            "point_indices":relation.as_ref().map(|relation| relation.point_indices),
            "probes":relation.as_ref().map(|relation| relation.probes),
            "s3_counts":query_counts.as_json(),
            "recovered_scalar":recovered,
            "recovered_matches_published":published.and_then(|scalar| recovered.map(|d| d == scalar % r)),
            "group_verified":verified,
            "target_generation_ms_excluded":target_generation_ms,
            "target_query_ms":target_query_stage_ms,
            "target_pdp_and_relation_check_ms":target_pdp_and_relation_check_ms,
            "target_descent_ms":target_descent_ms,
            "target_recovery_check_ms":target_recovery_check_ms,
            "target_phase_sum_ms":target_query_stage_ms + target_pdp_and_relation_check_ms + target_descent_ms + target_recovery_check_ms,
            "target_ms":elapsed,
        });
        writeln!(out, "{record}").unwrap();
    }
    target_counts.validate(window);
    let targets_ms = targets_started.elapsed().as_secs_f64() * 1000.0;
    let total_ms = process_started.elapsed().as_secs_f64() * 1000.0;
    let mut sorted = target_query_ms.clone();
    sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let median = sorted.get(sorted.len() / 2).copied();
    let summary = json!({
        "kind":"compact_orbit_dlp_summary",
        "schema_version":"1.0",
        "n":n, "a":a,
        "base_hash":base_hash,
        "base_source":if constructed.is_some() { "constructed_in_process_x_scan" } else { "retained_header" },
        "base_scanned_x":scanned_x,
        "orbit_columns":base.columns,
        "factor_base_points":base.points.len(),
        "index_policy":"swap_frobenius_quotient",
        "s3_batch_window":window,
        "root_prefilter_policy":if use_prefilter { "blocked_bloom_512_3hash" } else { "off" },
        "root_prefilter_bytes":index.prefilter.as_ref().map_or(0, BlockedBloom::bytes),
        "root_prefilter_blocks":index.prefilter.as_ref().map_or(0, |filter| filter.blocks.len()),
        "index_s3_counts":index_counts.as_json(),
        "rank_s3_counts":rank_counts.as_json(),
        "target_s3_counts":target_counts.as_json(),
        "ordered_state_candidates":base.columns * base.columns * n as usize,
        "representative_state_candidates":index.representative_candidates,
        "regular_states":index.states.len(),
        "root_table_entries":index.table.len,
        "root_table_slots":index.table.slots.len(),
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
        "rank_rows_without_gain":rank_rows_without_gain,
        "rank_policy":"guided: decompose [a]G - R_j for the first pivotless column j",
        "rank_probes_mean":if rank_relations > 0 { rank_probes as f64 / rank_relations as f64 } else { 0.0 },
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
}

#[cfg(test)]
mod swap_quotient_tests {
    use super::*;
    use crypto_lib::binary_ecc::IrreduciblePoly;
    use std::collections::HashSet;

    #[test]
    fn exhaustive_gf32_scalar_and_batch_roots_match_in_original_order() {
        let gf = Gf2::new(&IrreduciblePoly {
            degree: 5,
            low_terms: vec![0, 2],
        });
        let basis = NormalBasis::new(&gf);
        let solver = S3Solver::new(&gf, 1);
        let inputs: Vec<(u64, u64)> = (0..32).flat_map(|x| (0..32).map(move |y| (x, y))).collect();
        let scalar: Vec<_> = inputs
            .iter()
            .map(|&(x, y)| solver.roots(&gf, &basis, x, y))
            .collect();
        let mut exceptional = 0;
        for window in [16, 64] {
            let mut scratch = RootBatchScratch::default();
            let mut counts = S3Counts::default();
            let mut batched = Vec::new();
            for chunk in inputs.chunks(window) {
                scratch.solve(&gf, &basis, &solver, chunk, &mut counts);
                counts.consumed_candidates += chunk.len() as u64;
                batched.extend_from_slice(&scratch.results);
            }
            counts.validate(window);
            exceptional = counts.exceptional_calls;
            assert_eq!(batched, scalar);
            assert_eq!(counts.calls, 1024);
            assert!(counts.batch_inversions > 0);
        }
        assert!(exceptional > 0);
        assert!(scalar.iter().any(Option::is_none));
        assert!(scalar.iter().any(Option::is_some));
    }

    #[test]
    fn query_batch_keeps_first_full_point_witness_and_charges_partial_tail() {
        let curve = KoblitzCurve::new(0, 41).unwrap();
        let fast = FastBinaryCurve::new(&curve.curve.irreducible, 0).unwrap();
        let gf = Gf2::new(&curve.curve.irreducible);
        let b = gf.from_element(&curve.curve.b);
        let basis = NormalBasis::new(&gf);
        let solver = S3Solver::new(&gf, b);
        let r = curve.subgroup_order.to_u64().unwrap();
        let (points, labels, representatives, _) = construct_base(&fast, &curve, b, 1, r);
        let reps = vec![representatives[0].unwrap()[0]];
        let mut by_x: HashMap<u64, Vec<usize>> = HashMap::new();
        for (index, point) in points.iter().enumerate() {
            by_x.entry(point.unwrap().0).or_default().push(index);
        }
        let base = Base {
            points,
            labels,
            by_x,
            columns: 1,
        };
        let mut index_counts = S3Counts::default();
        let index = build_index(&gf, &basis, &solver, &reps, 1, false, &mut index_counts);
        let mut filter_index_counts = S3Counts::default();
        let filtered = build_index(
            &gf,
            &basis,
            &solver,
            &reps,
            1,
            true,
            &mut filter_index_counts,
        );
        assert_eq!(filtered.table.slots, index.table.slots);
        assert_eq!(filter_index_counts.as_json(), index_counts.as_json());
        let mut found = None;
        for i in 0..8 {
            for j in 0..8 {
                let pair = fast.add(base.points[i], base.points[j]);
                let q = fast.add(pair, pair);
                let mut scalar_counts = S3Counts::default();
                if let Some(relation) = extract(
                    &gf,
                    &fast,
                    &basis,
                    &solver,
                    &index,
                    &base,
                    q,
                    0,
                    1,
                    &mut scalar_counts,
                ) {
                    found = Some((q, relation, scalar_counts));
                    break;
                }
            }
            if found.is_some() {
                break;
            }
        }
        let (q, expected, scalar_counts) =
            found.expect("a signed Frobenius orbit has a four-point witness");
        for window in [16, 64] {
            let mut counts = S3Counts::default();
            let actual = extract(
                &gf,
                &fast,
                &basis,
                &solver,
                &index,
                &base,
                q,
                0,
                window,
                &mut counts,
            )
            .unwrap();
            counts.validate(window);
            assert_eq!(actual.point_indices, expected.point_indices);
            assert_eq!(actual.x_codes, expected.x_codes);
            assert_eq!(actual.intermediates, expected.intermediates);
            assert_eq!(actual.probes, expected.probes);
            assert_eq!(
                counts.consumed_candidates,
                scalar_counts.consumed_candidates
            );
            assert!(counts.discarded_candidates < window as u64);
            let mut filter_counts = S3Counts::default();
            let filtered_relation = extract(
                &gf,
                &fast,
                &basis,
                &solver,
                &filtered,
                &base,
                q,
                0,
                window,
                &mut filter_counts,
            )
            .unwrap();
            filter_counts.validate(window);
            assert_eq!(filtered_relation.point_indices, actual.point_indices);
            assert_eq!(filtered_relation.x_codes, actual.x_codes);
            assert_eq!(filtered_relation.intermediates, actual.intermediates);
            assert_eq!(filtered_relation.probes, actual.probes);
            assert_eq!(
                filter_counts.root_keys_considered,
                counts.root_keys_considered
            );
            assert_eq!(filter_counts.calls, counts.calls);
        }
        // This deliberately invalid point cannot be a sum of curve points.
        // Its complete scan ends in a non-full batch, whose actual work is counted.
        let impossible = Some((0, 0));
        let mut scalar_counts = S3Counts::default();
        assert!(extract(
            &gf,
            &fast,
            &basis,
            &solver,
            &index,
            &base,
            impossible,
            0,
            1,
            &mut scalar_counts
        )
        .is_none());
        let mut batch_counts = S3Counts::default();
        assert!(extract(
            &gf,
            &fast,
            &basis,
            &solver,
            &index,
            &base,
            impossible,
            0,
            64,
            &mut batch_counts
        )
        .is_none());
        batch_counts.validate(64);
        assert_eq!(batch_counts.calls, scalar_counts.calls);
        assert_ne!(batch_counts.calls % 64, 0);
        assert_eq!(batch_counts.discarded_candidates, 0);
        let mut filtered_counts = S3Counts::default();
        assert!(extract(
            &gf,
            &fast,
            &basis,
            &solver,
            &filtered,
            &base,
            impossible,
            0,
            64,
            &mut filtered_counts
        )
        .is_none());
        filtered_counts.validate(64);
        assert_eq!(filtered_counts.calls, batch_counts.calls);
        assert_eq!(
            filtered_counts.root_keys_considered,
            batch_counts.root_keys_considered
        );
        assert!(filtered_counts.prefilter_definite_misses > 0);
    }

    #[test]
    fn exhaustive_gf32_root_filter_only_skips_absent_keys() {
        let gf = Gf2::new(&IrreduciblePoly {
            degree: 5,
            low_terms: vec![0, 2],
        });
        let basis = NormalBasis::new(&gf);
        let solver = S3Solver::new(&gf, 1);
        let mut table = RootTable::with_capacity(2048);
        let mut filter = BlockedBloom::with_capacity(2048);
        for left in 0..32 {
            for right in 0..32 {
                if let Some(roots) = solver.roots(&gf, &basis, left, right) {
                    for root in roots {
                        let (key, _) = basis.canonical(basis.to_normal.apply(root));
                        filter.insert(key);
                        table.insert_if_absent(key, left * 32 + right);
                    }
                }
            }
        }
        let mut skipped_absent = 0;
        for key in 0..4096u64 {
            let actual = table.get(key);
            if !filter.may_contain(key) {
                assert!(actual.is_none(), "false negative for canonical root {key}");
                skipped_absent += 1;
            }
        }
        assert!(table.len > 0);
        assert!(skipped_absent > 3000);
    }

    #[test]
    fn exhaustive_odd_degree_swap_preserves_roots_and_table_witnesses() {
        // All x in GF(2^5), including zero, equal inputs, and the solver's
        // a=0/ps=0 fallback branches. This checks the algebraic identity
        // independently of factor-base construction or target search order.
        let gf = Gf2::new(&IrreduciblePoly {
            degree: 5,
            low_terms: vec![0, 2],
        });
        let basis = NormalBasis::new(&gf);
        let solver = S3Solver::new(&gf, 1);
        let reps: Vec<u64> = (0..(1u64 << 5)).collect();
        let mut counts = S3Counts::default();
        let index = build_index(&gf, &basis, &solver, &reps, 1, false, &mut counts);
        counts.validate(1);
        let n = gf.n as usize;
        let ordered = reps.len() * reps.len() * n;
        assert_eq!(index.representative_candidates, (ordered + reps.len()) / 2);
        for left in 0..reps.len() {
            assert!(retain_swap_representative(left, left, 0, n));
        }

        let mut all_keys = HashSet::new();
        let mut regular_ordered = 0usize;
        let mut regular_fixed = 0usize;
        let mut exceptional_regular = 0usize;
        for left in 0..reps.len() {
            for right in 0..reps.len() {
                for relative in 0..n {
                    let a = index.shifted[left][0];
                    let b = index.shifted[right][relative];
                    let roots = solver.roots(&gf, &basis, a, b);
                    let mate_relative = (n - relative) % n;
                    let mate = solver.roots(
                        &gf,
                        &basis,
                        index.shifted[right][0],
                        index.shifted[left][mate_relative],
                    );
                    assert_eq!(roots.is_some(), mate.is_some());
                    if let (Some(mut original), Some(mut swapped)) = (roots, mate) {
                        regular_ordered += 1;
                        if left == right && relative == 0 {
                            regular_fixed += 1;
                        }
                        if a == b || a == 0 || b == 0 {
                            exceptional_regular += 1;
                        }
                        // Shifting the mate by relative aligns its unordered
                        // pair with the original state's unshifted pair.
                        for root in &mut swapped {
                            for _ in 0..relative {
                                *root = gf.sqr(*root);
                            }
                        }
                        original.sort_unstable();
                        swapped.sort_unstable();
                        assert_eq!(original, swapped, "({left},{right},{relative})");
                        for root in original {
                            let key = basis.canonical(basis.to_normal.apply(root)).0;
                            all_keys.insert(key);
                            assert!(index.table.get(key).is_some());
                        }
                    }
                }
            }
        }
        assert!(exceptional_regular > 0);
        assert_eq!(index.states.len(), (regular_ordered + regular_fixed) / 2);
        assert_eq!(index.table.len, all_keys.len());
        for key in all_keys {
            let encoded = index.table.get(key).unwrap();
            let (left, right, relative, stored_shift) = unpack(encoded);
            assert!(retain_swap_representative(left, right, relative, n));
            let roots = solver
                .roots(
                    &gf,
                    &basis,
                    index.shifted[left][0],
                    index.shifted[right][relative],
                )
                .unwrap();
            assert!(roots.into_iter().any(|root| {
                basis.canonical(basis.to_normal.apply(root)) == (key, stored_shift)
            }));
        }
    }
}
