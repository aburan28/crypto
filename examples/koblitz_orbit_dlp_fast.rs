//! Compact-orbit relation extraction with a shared-factor-log DLP pass.
//!
//! Public synthetic Koblitz fixtures only. Reads a retained
//! `point_defined_factor_base` header, builds the Frobenius-quotiented S3
//! regular-root index once (no pair table, no edge selectors), then
//!   1. rank stage: decomposes random `[a]G` until the K orbit logs are fixed,
//!   2. target stage: decomposes each published target and recovers its log.
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
        let chunks = (images.len() + 7) / 8;
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
        let raised = basis.to_poly.apply(basis.rotate(basis.to_normal.apply(c), k));
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

    fn insert_if_absent(&mut self, key: u64, value: u64) {
        let mask = self.slots.len() - 1;
        let mut slot = self.slot(key);
        loop {
            if self.slots[slot].0 == u64::MAX {
                self.slots[slot] = (key, value);
                self.len += 1;
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

fn build_index(gf: &Gf2, basis: &NormalBasis, solver: &S3Solver, reps: &[u64]) -> Index {
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
    for left in 0..reps.len() {
        for right in 0..reps.len() {
            for relative in 0..n {
                if let Some(roots) = solver.roots(gf, basis, shifted[left][0], shifted[right][relative]) {
                    states.push(State {
                        left: left as u16,
                        right: right as u16,
                        relative: relative as u16,
                        normal_roots: [basis.to_normal.apply(roots[0]), basis.to_normal.apply(roots[1])],
                    });
                }
            }
        }
    }
    let mut table = RootTable::with_capacity(states.len() * 2);
    for state in &states {
        for &root in &state.normal_roots {
            let (canonical, shift) = basis.canonical(root);
            table.insert_if_absent(canonical, pack(state.left, state.right, state.relative, shift));
        }
    }
    Index { states, table, shifted }
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

fn lift(fast: &FastBinaryCurve, base: &Base, codes: &[u64; 4], target: FastPoint) -> Option<[usize; 4]> {
    let choices: Vec<&Vec<usize>> = codes.iter().map(|code| base.by_x.get(code)).collect::<Option<_>>()?;
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
) -> Option<Relation> {
    let (target_x, _) = target?;
    let n = gf.n;
    let mut probes = 0u64;
    // Callers provide a rotating, target-dependent origin so rank-stage rows
    // do not concentrate on the same early columns.
    let start = start % index.states.len().max(1);
    for state in index.states[start..].iter().chain(&index.states[..start]) {
        for shift in 0..n {
            let left_x = index.shifted[state.left as usize][shift as usize];
            let right_x = index.shifted[state.right as usize][(shift as usize + state.relative as usize) % n as usize];
            for &normal_root in &state.normal_roots {
                let absolute = basis.to_poly.apply(basis.rotate(normal_root, shift));
                let Some(partners) = solver.roots(gf, basis, absolute, target_x) else {
                    continue;
                };
                for partner in partners {
                    probes += 1;
                    let (canonical, partner_shift) = basis.canonical(basis.to_normal.apply(partner));
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
                    if let Some(point_indices) = lift(fast, base, &codes, target) {
                        return Some(Relation { point_indices, x_codes: codes, intermediates: [absolute, partner], probes });
                    }
                }
            }
        }
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
) -> (Vec<FastPoint>, Vec<(usize, u64)>, Vec<Option<[u64; 2]>>, u64) {
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
    let kib: u64 = String::from_utf8_lossy(&output.stdout).trim().parse().ok()?;
    Some(kib * 1024)
}

// ── Wide (`u128`-word) compact-orbit pipeline for `64 < n ≤ 127` ──
//
// Twin of the `u64` path above with identical stage semantics; the
// `u64` code is deliberately untouched so every `n ≤ 63` fixture stays
// byte-identical.  Subgroup scalars, orbit indices, packed witnesses,
// and the rank linear algebra stay `u64` (every admitted `r < 2^64`
// through `n = 127`); only field elements move to the wider word.
mod wide {
    use super::{pack, unpack, Echelon, mulmod, peak_rss_bytes};
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
            let chunks = (images.len() + 7) / 8;
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
    }

    /// Same Gauss-Jordan as the twin, with rows as `(m, i)` word pairs
    /// (the M and I halves no longer share one word past 64 columns).
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
                square: Linear128::from_images(&unit.iter().map(|&e| gf.sqr(e)).collect::<Vec<_>>()),
                half_trace: Linear128::from_images(&half_trace_images),
                trace_mask,
                b,
            }
        }

        #[inline(always)]
        fn roots(
            &self,
            gf: &Gf2_128,
            basis: &NormalBasis128,
            left: u128,
            right: u128,
        ) -> Option<[u128; 2]> {
            let s = left ^ right;
            let a = self.square.apply(s);
            let p = gf.mul(left, right);
            let ps = self.square.apply(p);
            if a == 0 || ps == 0 {
                return s3_x_roots_128(gf, self.b, left, right);
            }
            let inv_combined = invert128(gf, basis, gf.mul(a, ps));
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

    #[inline(always)]
    fn invert128(gf: &Gf2_128, basis: &NormalBasis128, a: u128) -> u128 {
        let e = gf.n - 1;
        let mut c = a;
        let mut k = 1u32;
        let bits = 32 - e.leading_zeros();
        for i in (0..bits - 1).rev() {
            let raised = basis.to_poly.apply(basis.rotate(basis.to_normal.apply(c), k));
            c = gf.mul(c, raised);
            k *= 2;
            if (e >> i) & 1 == 1 {
                c = gf.mul(gf.sqr(c), a);
                k += 1;
            }
        }
        gf.sqr(c)
    }

    /// Open-addressing `u128 -> u64` table; `u128::MAX` marks empty
    /// (keys are `< 2^71`, so the marker is unreachable).
    struct RootTable128 {
        slots: Vec<(u128, u64)>,
        shift: u32,
        len: usize,
    }

    impl RootTable128 {
        fn with_capacity(entries: usize) -> Self {
            let slots = (entries * 2).next_power_of_two().max(16);
            Self {
                slots: vec![(u128::MAX, 0); slots],
                shift: 128 - slots.trailing_zeros(),
                len: 0,
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

    /// Deterministic rank-stage scalar / start derivation for the parallel
    /// guided rank: each column `j` uses `mix(rank_seed, j, attempt)`, so
    /// the set of extractions is independent of thread scheduling.  The
    /// resulting factor-log table is the unique full-rank solution either
    /// way, so the published target row is unchanged.
    fn rank_mix(seed: u64, column: u64, attempt: u64) -> u64 {
        let mut z = seed
            ^ column.wrapping_mul(0x9E37_79B9_7F4A_7C15)
            ^ attempt.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Order-preserving parallel target extraction: workers scan disjoint
    /// blocks of the same rotated state order the sequential scan uses, each
    /// stopping at its block's first hit; the merge picks the smallest
    /// global position.  The published relation (point indices, x-codes,
    /// intermediates) is therefore identical to the sequential first hit —
    /// only the wall time and the per-block probe diagnostic differ.
    /// Returns the relation and the summed relation-check nanoseconds.
    fn extract128_parallel(
        gf: &Gf2_128,
        fast: &FastBinaryCurve128,
        basis: &NormalBasis128,
        solver: &S3Solver128,
        index: &Index128,
        base: &Base128,
        target: FastPoint128,
        start: usize,
        threads: usize,
    ) -> (Option<Relation128>, u128) {
        let n_states = index.states.len();
        if threads <= 1 || n_states == 0 {
            let mut relation_check_ns = 0u128;
            let relation = extract128(
                gf,
                fast,
                basis,
                solver,
                index,
                base,
                target,
                start,
                &mut relation_check_ns,
            );
            return (relation, relation_check_ns);
        }
        let (target_x, _) = match target {
            Some((x, _)) => (x, ()),
            None => return (None, 0),
        };
        let n_u32 = gf.n;
        let n = n_u32 as usize;
        let chunk = (n_states + threads - 1) / threads;
        let (sender, receiver) = std::sync::mpsc::channel::<(usize, Relation128, u128)>();
        std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(threads);
            for thread in 0..threads {
                let sender = sender.clone();
                let gf = gf;
                let fast = fast;
                let basis = basis;
                let solver = solver;
                let index = index;
                let base = base;
                handles.push(scope.spawn(move || {
                    let pos_lo = thread * chunk;
                    let pos_hi = ((thread + 1) * chunk).min(n_states);
                    if pos_lo >= pos_hi {
                        return;
                    }
                    // Rotated position p maps to absolute state
                    // (start + p) mod n_states; scan in order.
                    let state_lo = (start + pos_lo) % n_states;
                    let length = pos_hi - pos_lo;
                    let mut relation_check_ns = 0u128;
                    let mut probes = 0u64;
                    let mut offset = 0usize;
                    'scan: while offset < length {
                        let state = &index.states[(state_lo + offset) % n_states];
                        for shift in 0..n_u32 {
                            let left_x = index.shifted[state.left as usize][shift as usize];
                            let right_x = index.shifted[state.right as usize]
                                [(shift as usize + state.relative as usize) % n];
                            for &normal_root in &state.normal_roots {
                                let absolute =
                                    basis.to_poly.apply(basis.rotate(normal_root, shift));
                                let Some(partners) =
                                    solver.roots(gf, basis, absolute, target_x)
                                else {
                                    continue;
                                };
                                for partner in partners {
                                    probes += 1;
                                    let (canonical, partner_shift) =
                                        basis.canonical(basis.to_normal.apply(partner));
                                    let Some(value) = index.table.get(canonical) else {
                                        continue;
                                    };
                                    let (left2, right2, relative2, stored_shift) = unpack(value);
                                    let left_shift2 =
                                        ((stored_shift + n_u32 - partner_shift) % n_u32) as usize;
                                    let right_shift2 = (left_shift2 + relative2 as usize) % n;
                                    let codes = [
                                        left_x,
                                        right_x,
                                        index.shifted[left2][left_shift2],
                                        index.shifted[right2][right_shift2],
                                    ];
                                    if let Some(point_indices) = lift128(
                                        fast,
                                        base,
                                        &codes,
                                        target,
                                        &mut relation_check_ns,
                                    ) {
                                        let position = pos_lo + offset;
                                        sender
                                            .send((
                                                position,
                                                Relation128 {
                                                    point_indices,
                                                    x_codes: codes,
                                                    intermediates: [absolute, partner],
                                                    probes,
                                                },
                                                relation_check_ns,
                                            ))
                                            .expect("target receiver is alive");
                                        break 'scan;
                                    }
                                }
                            }
                        }
                        offset += 1;
                    }
                }));
            }
            drop(sender);
            let mut best: Option<(usize, Relation128, u128)> = None;
            let mut check_ns_total = 0u128;
            for (position, relation, check_ns) in receiver {
                check_ns_total += check_ns;
                if best.as_ref().map_or(true, |(pos, _, _)| position < *pos) {
                    best = Some((position, relation, check_ns));
                }
            }
            for handle in handles {
                handle.join().expect("target worker thread");
            }
            (best.map(|(_, relation, _)| relation), check_ns_total)
        })
    }

    // serde_json's default Number representation cannot emit arbitrary
    // u128 field elements.  Keep the wide-only wire format lossless by
    // encoding every field coordinate as a decimal string; the legacy
    // u64 pipeline and its JSON fixtures remain unchanged.
    fn encode_point128(point: FastPoint128) -> Option<[String; 2]> {
        point.map(|(x, y)| [x.to_string(), y.to_string()])
    }

    fn encode_pair128(point: Option<[u128; 2]>) -> Option<[String; 2]> {
        point.map(|[x, y]| [x.to_string(), y.to_string()])
    }

    fn decode_word128(value: &Value) -> u128 {
        if let Some(decimal) = value.as_str() {
            decimal.parse().expect("wide field coordinate must be decimal")
        } else if let Some(word) = value.as_u64() {
            word as u128
        } else {
            panic!("wide field coordinate must be a decimal string or u64")
        }
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
                    if let Some(roots) = solver.roots(gf, basis, shifted[left][0], shifted[right][relative]) {
                        states.push(State128 {
                            left: left as u16,
                            right: right as u16,
                            relative: relative as u16,
                            normal_roots: [basis.to_normal.apply(roots[0]), basis.to_normal.apply(roots[1])],
                        });
                    }
                }
            }
        }
        let mut table = RootTable128::with_capacity(states.len() * 2);
        for state in &states {
            for &root in &state.normal_roots {
                let (canonical, shift) = basis.canonical(root);
                table.insert_if_absent(canonical, pack(state.left, state.right, state.relative, shift));
            }
        }
        Index128 { states, table, shifted }
    }

    fn lift128(
        fast: &FastBinaryCurve128,
        base: &Base128,
        codes: &[u128; 4],
        target: FastPoint128,
        relation_check_ns: &mut u128,
    ) -> Option<[usize; 4]> {
        let choices: Vec<&Vec<usize>> = codes.iter().map(|code| base.by_x.get(code)).collect::<Option<_>>()?;
        for &a in choices[0] {
            for &b in choices[1] {
                let ab = fast.add(base.points[a], base.points[b]);
                for &c in choices[2] {
                    let abc = fast.add(ab, base.points[c]);
                    for &d in choices[3] {
                        let check_started = Instant::now();
                        let matches = fast.add(abc, base.points[d]) == target;
                        *relation_check_ns += check_started.elapsed().as_nanos();
                        if matches {
                            return Some([a, b, c, d]);
                        }
                    }
                }
            }
        }
        None
    }

    fn extract128(
        gf: &Gf2_128,
        fast: &FastBinaryCurve128,
        basis: &NormalBasis128,
        solver: &S3Solver128,
        index: &Index128,
        base: &Base128,
        target: FastPoint128,
        start: usize,
        relation_check_ns: &mut u128,
    ) -> Option<Relation128> {
        let (target_x, _) = target?;
        let n = gf.n;
        let mut probes = 0u64;
        let start = start % index.states.len().max(1);
        for state in index.states[start..].iter().chain(&index.states[..start]) {
            for shift in 0..n {
                let left_x = index.shifted[state.left as usize][shift as usize];
                let right_x = index.shifted[state.right as usize][(shift as usize + state.relative as usize) % n as usize];
                for &normal_root in &state.normal_roots {
                    let absolute = basis.to_poly.apply(basis.rotate(normal_root, shift));
                    let Some(partners) = solver.roots(gf, basis, absolute, target_x) else {
                        continue;
                    };
                    for partner in partners {
                        probes += 1;
                        let (canonical, partner_shift) = basis.canonical(basis.to_normal.apply(partner));
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
                        if let Some(point_indices) = lift128(
                            fast,
                            base,
                            &codes,
                            target,
                            relation_check_ns,
                        ) {
                            return Some(Relation128 { point_indices, x_codes: codes, intermediates: [absolute, partner], probes });
                        }
                    }
                }
            }
        }
        None
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
    ) -> (Vec<FastPoint128>, Vec<(usize, u64)>, Vec<Option<[u128; 2]>>, u64) {
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

    /// Wide pipeline for `64 < n ≤ 127`: same stages, timers, and JSON
    /// shapes as the `u64` flow (field values widen to `u128`).
    pub fn run(arguments: &[String]) -> ! {
        assert_eq!(arguments.len(), 5);
        let process_started = Instant::now();
        let constructed: Option<(u32, u8, usize)> = arguments[1].strip_prefix("construct:").map(|spec| {
            let parts: Vec<&str> = spec.split(':').collect();
            assert_eq!(parts.len(), 3, "construct:<n>:<a>:<columns>");
            (parts[0].parse().unwrap(), parts[1].parse().unwrap(), parts[2].parse().unwrap())
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
                    let encoded: Value = serde_json::from_str(line)
                        .expect("target point must be a JSON [x,y] array");
                    let [x, y] = decode_pair128(&encoded)
                        .expect("target point must be affine");
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
                curve.curve.irreducible.low_terms.iter().map(|&t| t as u64).collect::<Vec<_>>(),
                "base header field must match the constructed curve"
            );
        }
        let fast = FastBinaryCurve128::new(&curve.curve.irreducible, a as u128).unwrap();
        let gf = Gf2_128::new(&curve.curve.irreducible);
        let b = gf.from_element(&curve.curve.b);
        let generator: FastPoint128 = match curve.generator() {
            crypto_lib::binary_ecc::BinaryPoint::Affine { x, y } => Some((fast.word(x), fast.word(y))),
            crypto_lib::binary_ecc::BinaryPoint::Infinity => panic!("generator must be affine"),
        };
        let mut scanned_x = None;
        let (points, labels, representatives): (Vec<FastPoint128>, Vec<(usize, u64)>, Vec<Option<[u128; 2]>>) =
            match constructed {
                Some((_, _, columns)) => {
                    let (points, labels, reps, scanned) = construct_base128(&fast, &curve, b, columns, r);
                    scanned_x = Some(scanned);
                    (points, labels, reps)
                }
                None => {
                    let coordinates =
                        decode_point_list128(&header["factor_base_point_coordinates"]);
                    let representatives =
                        decode_point_list128(&header["factor_base_representatives"]);
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
        let base = Base128 { points, labels, by_x, columns: representatives.len() };
        let reps: Vec<u128> = representatives.iter().map(|p| p.expect("affine representative")[0]).collect();
        let base_load_ms = setup_started.elapsed().as_secs_f64() * 1000.0;

        let basis_started = Instant::now();
        let basis = NormalBasis128::new(&gf);
        let solver = S3Solver128::new(&gf, b);
        for probe in 1..2048u128 {
            let x = probe.wrapping_mul(0x2545_F491_4F6C_DD1D) & basis.mask;
            assert_eq!(basis.to_poly.apply(basis.rotate(basis.to_normal.apply(x), 1)), gf.sqr(x));
            let y = (probe * 0x9E37_79B9) & basis.mask;
            assert_eq!(solver.roots(&gf, &basis, x, y), s3_x_roots_128(&gf, b, x, y));
        }
        let basis_ms = basis_started.elapsed().as_secs_f64() * 1000.0;

        let index_started = Instant::now();
        let index = build_index128(&gf, &basis, &solver, &reps);
        let index_ms = index_started.elapsed().as_secs_f64() * 1000.0;
        let setup_ms = setup_started.elapsed().as_secs_f64() * 1000.0;

        let rank_started = Instant::now();
        let mut echelon = super::Echelon { r, columns: base.columns, pivots: vec![None; base.columns], rank: 0 };
        let mut rank_attempts = 0u64;
        let mut rank_failures = 0u64;
        let mut rank_relations = 0u64;
        let mut rank_probes = 0u64;
        let mut rank_rows_without_gain = 0u64;
        let mut rank_unused_relation_check_ns = 0u128;
        let rank_threads: usize = std::env::var("KIC_RANK_THREADS")
            .ok()
            .and_then(|value| value.trim().parse::<usize>().ok())
            .filter(|&threads| threads > 0)
            .unwrap_or(1);
        let rank_policy_parallel = rank_threads > 1;
        if rank_policy_parallel {
            // Parallel guided rank: every column is anchored by its own
            // relation with a scheduling-independent scalar, so the row set
            // (and therefore the solved logs, the unique full-rank
            // solution) does not depend on thread interleaving.  Rows may
            // arrive covering already-pivoted columns; those count as
            // rows without gain, exactly as sequential duplicates would.
            let (sender, receiver) =
                std::sync::mpsc::channel::<(usize, u64, Relation128, u64, u128, u64)>();
            std::thread::scope(|scope| {
                let mut handles = Vec::with_capacity(rank_threads);
                for thread in 0..rank_threads {
                    let sender = sender.clone();
                    let gf = &gf;
                    let fast = &fast;
                    let basis = &basis;
                    let solver = &solver;
                    let index = &index;
                    let base = &base;
                    let generator = &generator;
                    let representatives = &representatives;
                    handles.push(scope.spawn(move || {
                        let mut thread_attempts = 0u64;
                        let mut thread_relation_check_ns = 0u128;
                        for column in (thread..base.columns).step_by(rank_threads) {
                            let rep: FastPoint128 = representatives[column].map(|[x, y]| (x, y));
                            let mut attempt = 0u64;
                            loop {
                                let mixed = rank_mix(rank_seed, column as u64, attempt);
                                let scalar = mixed % (r - 1) + 1;
                                let start = (mixed >> 20) as usize;
                                let point = fast.add(
                                    fast.scalar_mul(*generator, &BigUint::from(scalar)),
                                    FastBinaryCurve128::neg(rep),
                                );
                                thread_attempts += 1;
                                match extract128(
                                    gf,
                                    fast,
                                    basis,
                                    solver,
                                    index,
                                    base,
                                    point,
                                    start,
                                    &mut thread_relation_check_ns,
                                ) {
                                    Some(relation) => {
                                        let probes = relation.probes;
                                        sender
                                            .send((
                                                column,
                                                scalar,
                                                relation,
                                                probes,
                                                thread_relation_check_ns,
                                                thread_attempts,
                                            ))
                                            .expect("rank receiver is alive");
                                        thread_relation_check_ns = 0;
                                        thread_attempts = 0;
                                        break;
                                    }
                                    None => attempt += 1,
                                }
                            }
                        }
                    }));
                }
                drop(sender);
                for (column, scalar, relation, probes, relation_check_ns, attempts) in receiver {
                    rank_attempts += attempts;
                    rank_relations += 1;
                    rank_probes += probes;
                    rank_unused_relation_check_ns += relation_check_ns;
                    let mut row = relation_row128(&base, &relation, scalar, r);
                    row[column] = (row[column] + 1) % r;
                    if !echelon.insert(row) {
                        rank_rows_without_gain += 1;
                    }
                }
                for handle in handles {
                    handle.join().expect("rank worker thread");
                }
            });
        } else {
            while echelon.rank < base.columns {
                rank_seed = rank_seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                let scalar = (rank_seed >> 11) % (r - 1) + 1;
                let column = (0..base.columns).find(|&c| echelon.pivots[c].is_none()).unwrap();
                let rep: FastPoint128 = representatives[column].map(|[x, y]| (x, y));
                let point = fast.add(fast.scalar_mul(generator, &BigUint::from(scalar)), FastBinaryCurve128::neg(rep));
                rank_attempts += 1;
                let mut unused_relation_check_ns = 0u128;
                match extract128(
                    &gf,
                    &fast,
                    &basis,
                    &solver,
                    &index,
                    &base,
                    point,
                    (rank_seed >> 20) as usize,
                    &mut unused_relation_check_ns,
                ) {
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
        }
        let rank_ms = rank_started.elapsed().as_secs_f64() * 1000.0;
        let la_started = Instant::now();
        let logs = echelon.solve();
        let la_ms = la_started.elapsed().as_secs_f64() * 1000.0;

        let targets_started = Instant::now();
        let target_threads: usize = std::env::var("KIC_TARGET_THREADS")
            .ok()
            .and_then(|value| value.trim().parse::<usize>().ok())
            .filter(|&threads| threads > 0)
            .unwrap_or(1);
        let mut solved = 0usize;
        let mut failed = 0usize;
        let mut target_query_ms = Vec::with_capacity(target_inputs.len());
        for (fixture_index, &target_input) in target_inputs.iter().enumerate() {
            let (target, published, target_generation_ms) = match target_input {
                TargetInput128::PublicPoint(target) => (target, None, 0.0),
                TargetInput128::KnownAnswerScalar(scalar) => {
                    let target_generation_started = Instant::now();
                    let target = fast.scalar_mul(generator, &BigUint::from(scalar));
                    let elapsed = target_generation_started.elapsed().as_secs_f64() * 1000.0;
                    (target, Some(scalar), elapsed)
                }
            };
            let query_started = Instant::now();
            let query_hash_started = Instant::now();
            let (target_x, target_y) = target.expect("published fixture target is not infinity");
            let mut query_hash = target_x ^ target_y.rotate_left(29) ^ 0x9E37_79B9_7F4A_7C15_9E37_79B9_7F4A_7C15;
            query_hash = (query_hash ^ (query_hash >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9_BF58_476D_1CE4_E5B9);
            query_hash = (query_hash ^ (query_hash >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB_94D0_49BB_1331_11EB);
            query_hash ^= query_hash >> 31;
            let start = query_hash as usize;
            let target_query_stage_ms = query_hash_started.elapsed().as_secs_f64() * 1000.0;
            let target_pdp_started = Instant::now();
            let (relation, target_relation_check_ns) = extract128_parallel(
                &gf,
                &fast,
                &basis,
                &solver,
                &index,
                &base,
                target,
                start,
                target_threads,
            );
            let target_decomposition_total_ms = target_pdp_started.elapsed().as_secs_f64() * 1000.0;
            let target_relation_check_ms = target_relation_check_ns as f64 / 1_000_000.0;
            let target_pdp_ms = (target_decomposition_total_ms - target_relation_check_ms).max(0.0);
            let target_descent_started = Instant::now();
            let recovered = relation.as_ref().map(|relation| {
                relation.point_indices.iter().fold(0u64, |acc, &index| {
                    let (column, coefficient) = base.labels[index];
                    (acc + super::mulmod(coefficient, logs[column], r)) % r
                })
            });
            let target_descent_ms = target_descent_started.elapsed().as_secs_f64() * 1000.0;
            let target_recovery_check_started = Instant::now();
            let verified = recovered.map(|d| fast.scalar_mul(generator, &BigUint::from(d)) == target);
            let target_recovery_check_ms = target_recovery_check_started.elapsed().as_secs_f64() * 1000.0;
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
                "field_value_encoding":"base-10 strings for u128 field values",
                "generator":encode_point128(generator),
                "target":encode_point128(target),
                "published_q":encode_point128(target),
                "exit_code":0,
                "x_codes":relation.as_ref().map(|relation| relation.x_codes.map(|x| x.to_string())),
                "pinned_intermediates":relation.as_ref().map(|relation| relation.intermediates.map(|x| x.to_string())),
                "point_indices":relation.as_ref().map(|relation| relation.point_indices),
                "probes":relation.as_ref().map(|relation| relation.probes),
                "recovered_scalar":recovered,
                "recovered_matches_published":published.and_then(|scalar| recovered.map(|d| d == scalar % r)),
                "group_verified":verified,
                "target_generation_ms_excluded":target_generation_ms,
                "target_query_ms":target_query_stage_ms,
                "target_pdp_ms":target_pdp_ms,
                "target_relation_check_ms":target_relation_check_ms,
                "target_descent_ms":target_descent_ms,
                "target_recovery_check_ms":target_recovery_check_ms,
                "target_phase_sum_ms":target_query_stage_ms + target_pdp_ms + target_relation_check_ms + target_descent_ms + target_recovery_check_ms,
                "target_ms":elapsed,
            });
            writeln!(out, "{record}").unwrap();
        }
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
            "regular_states":index.states.len(),
            "root_table_entries":index.table.len,
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
            "rank_policy": if rank_policy_parallel {
                "guided parallel: every column anchors [a]G - R_j with a scheduling-independent scalar; logs are the unique full-rank solution"
            } else {
                "guided: decompose [a]G - R_j for the first pivotless column j"
            },
            "rank_probes_mean":if rank_relations > 0 { rank_probes as f64 / rank_relations as f64 } else { 0.0 },
            "rank":echelon.rank,
            "targets":target_inputs.len(),
            "targets_solved":solved,
            "targets_failed":failed,
            "peak_rss_bytes":super::peak_rss_bytes(),
            "threads":rank_threads,
            "target_threads":target_threads,
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
        std::process::exit(0);
    }

    #[cfg(test)]
    mod wide_json_tests {
        use super::{decode_pair128, decode_point_list128, encode_pair128};
        use serde_json::json;

        #[test]
        fn wide_point_coordinates_round_trip_as_decimal_strings() {
            let point = Some([(1u128 << 100) + 7, (1u128 << 80) + 9]);
            let encoded = serde_json::to_value(encode_pair128(point)).unwrap();
            assert_eq!(
                encoded,
                json!([((1u128 << 100) + 7).to_string(), ((1u128 << 80) + 9).to_string()])
            );
            assert_eq!(decode_pair128(&encoded), point);
        }

        #[test]
        fn wide_point_parser_accepts_small_legacy_numeric_words() {
            let encoded = json!([17u64, u64::MAX]);
            assert_eq!(decode_pair128(&encoded), Some([17, u64::MAX as u128]));
        }

        #[test]
        fn wide_point_list_preserves_infinity_entries() {
            let encoded = json!([["18446744073709551616", "9"], null]);
            assert_eq!(
                decode_point_list128(&encoded),
                vec![Some([1u128 << 64, 9]), None]
            );
        }
    }
}

fn main() {
    let arguments: Vec<String> = std::env::args().collect();
    assert_eq!(
        arguments.len(),
        5,
        "usage: <base_header.jsonl | construct:<n>:<a>:<columns>> <target_points_or_legacy_scalars.txt> <rank_seed> <out.jsonl>"
    );
    // Wide pipeline (`u128` words) for 64 < n ≤ 127.  The degree is
    // detected here so the `u64` flow below stays byte-identical;
    // the wide path re-parses its own arguments.
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
    let constructed: Option<(u32, u8, usize)> = arguments[1].strip_prefix("construct:").map(|spec| {
        let parts: Vec<&str> = spec.split(':').collect();
        assert_eq!(parts.len(), 3, "construct:<n>:<a>:<columns>");
        (parts[0].parse().unwrap(), parts[1].parse().unwrap(), parts[2].parse().unwrap())
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
                let [x, y]: [u64; 2] = serde_json::from_str(line).expect("target point must be [x,y]");
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
            curve.curve.irreducible.low_terms.iter().map(|&t| t as u64).collect::<Vec<_>>(),
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
    let (points, labels, representatives): (Vec<FastPoint>, Vec<(usize, u64)>, Vec<Option<[u64; 2]>>) =
        match constructed {
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
    let base = Base { points, labels, by_x, columns: representatives.len() };
    let reps: Vec<u64> = representatives.iter().map(|p| p.expect("affine representative")[0]).collect();
    let base_load_ms = setup_started.elapsed().as_secs_f64() * 1000.0;

    let basis_started = Instant::now();
    let basis = NormalBasis::new(&gf);
    let solver = S3Solver::new(&gf, b);
    for probe in 1..2048u64 {
        let x = probe.wrapping_mul(0x2545_F491_4F6C_DD1D) & basis.mask;
        assert_eq!(basis.to_poly.apply(basis.rotate(basis.to_normal.apply(x), 1)), gf.sqr(x));
        let y = (probe * 0x9E37_79B9) & basis.mask;
        assert_eq!(solver.roots(&gf, &basis, x, y), s3_x_roots(&gf, b, x, y));
    }
    let basis_ms = basis_started.elapsed().as_secs_f64() * 1000.0;

    let index_started = Instant::now();
    let index = build_index(&gf, &basis, &solver, &reps);
    let index_ms = index_started.elapsed().as_secs_f64() * 1000.0;
    let setup_ms = setup_started.elapsed().as_secs_f64() * 1000.0;

    // Rank stage: rows sum(coefficient * L_column) = a for random [a]G.
    let rank_started = Instant::now();
    let mut echelon = Echelon { r, columns: base.columns, pivots: vec![None; base.columns], rank: 0 };
    let mut rank_attempts = 0u64;
    let mut rank_failures = 0u64;
    let mut rank_relations = 0u64;
    let mut rank_probes = 0u64;
    // Guided: decompose [a]G - R_j for a column j that has no pivot yet, so the
    // row a = L_j + sum(...) always contains column j and raises the rank.
    let mut rank_rows_without_gain = 0u64;
    while echelon.rank < base.columns {
        rank_seed = rank_seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let scalar = (rank_seed >> 11) % (r - 1) + 1;
        let column = (0..base.columns).find(|&c| echelon.pivots[c].is_none()).unwrap();
        let rep: FastPoint = representatives[column].map(|[x, y]| (x, y));
        let point = fast.add(fast.scalar_mul(generator, &BigUint::from(scalar)), FastBinaryCurve::neg(rep));
        rank_attempts += 1;
        match extract(&gf, &fast, &basis, &solver, &index, &base, point, (rank_seed >> 20) as usize) {
            Some(relation) => {
                rank_relations += 1;
                rank_probes += relation.probes;
                let mut row = relation_row(&base, &relation, scalar, r);
                row[column] = (row[column] + 1) % r;
                if !echelon.insert(row) {
                    rank_rows_without_gain += 1;
                }
            }
            None => rank_failures += 1,
        }
    }
    let rank_ms = rank_started.elapsed().as_secs_f64() * 1000.0;
    let la_started = Instant::now();
    let logs = echelon.solve();
    let la_ms = la_started.elapsed().as_secs_f64() * 1000.0;

    let targets_started = Instant::now();
    let mut solved = 0usize;
    let mut failed = 0usize;
    let mut target_query_ms = Vec::with_capacity(target_inputs.len());
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
        let relation = extract(&gf, &fast, &basis, &solver, &index, &base, target, start);
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
        let target_recovery_check_ms = target_recovery_check_started.elapsed().as_secs_f64() * 1000.0;
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
        "regular_states":index.states.len(),
        "root_table_entries":index.table.len,
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
