//! Free-relation census: self-collisions in the S3 pair-root table.
//!
//! Public synthetic Koblitz fixtures only. Loads (or constructs) the same
//! `point_defined_factor_base` used by `koblitz_orbit_dlp_fast`, enumerates
//! the relative pair-states `(left, right, relative_shift)` and their two S3
//! roots exactly as `build_index` does there, but instead of keeping only the
//! first state per canonical root, keeps every `(canonical_root, state,
//! root_index, rotation)` entry, sorts it, and reads off every pair of
//! distinct pair-states that share a root up to Frobenius rotation.
//!
//! Each such collision is a candidate four-term relation among factor-base
//! points. Every candidate is verified as an exact group identity
//! (`sum of the four subgroup points == O`) before it is counted, and the
//! verified rows are inserted into the same incremental echelon (mod r) the
//! rank stage uses, so the reported rank is measured, not assumed.
//!
//! No claim is made here about ECC2K-130, no target is solved, and no
//! vs_rho comparison is implied. This is a stage census for the free-relation
//! prediction in `RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md` §3.
//!
//! Usage: <base_header.jsonl | construct:<n>:<a>:<columns>> <out.jsonl>

use crypto_lib::cryptanalysis::koblitz_fast_arith::FastBinaryCurve;
use crypto_lib::cryptanalysis::koblitz_fast_arith::FastPoint;
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
            return crypto_lib::cryptanalysis::koblitz_fast_arith::s3_x_roots(
                gf, self.b, left, right,
            );
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

#[derive(Clone, Copy)]
struct State {
    left: u16,
    right: u16,
    relative: u16,
}

/// One `(canonical root, state, root index, rotation to canonical)` entry.
/// `rotation` is the shift c with `rotate(root, c) == canonical`.
#[derive(Clone, Copy)]
struct Entry {
    key: u64,
    state: u32,
    #[allow(dead_code)]
    root_index: u8,
    rotation: u8,
}

struct Base {
    points: Vec<FastPoint>,
    labels: Vec<(usize, u64)>,
    by_x: HashMap<u64, Vec<usize>>,
    columns: usize,
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

/// Frobenius-closed signed-orbit base from a deterministic public x-scan,
/// identical to `koblitz_orbit_dlp_fast::construct_base`.
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

/// Incremental row echelon mod prime r over K unknowns (homogeneous rows).
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
}

fn main() {
    let arguments: Vec<String> = std::env::args().collect();
    assert_eq!(
        arguments.len(),
        3,
        "usage: <base_header.jsonl | construct:<n>:<a>:<columns>> <out.jsonl>"
    );
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
        Some((n, a, _)) => json!({"n": n, "a": a}),
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
    assert!(n <= 63, "this census covers the u64 path (n <= 63)");

    let curve = KoblitzCurve::new(a, n).expect("admitted Koblitz rung");
    let r = curve.subgroup_order.to_u64().unwrap();
    let fast = FastBinaryCurve::new(&curve.curve.irreducible, a as u64).unwrap();
    let gf = Gf2::new(&curve.curve.irreducible);
    let b = gf.from_element(&curve.curve.b);

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
        assert_eq!(
            solver.roots(&gf, &basis, x, y),
            crypto_lib::cryptanalysis::koblitz_fast_arith::s3_x_roots(&gf, b, x, y)
        );
    }

    // State scan: identical triple set to `build_index`, both roots kept.
    let scan_started = Instant::now();
    let n_usize = n as usize;
    let shifted: Vec<Vec<u64>> = reps
        .iter()
        .map(|&code| {
            let mut row = Vec::with_capacity(n_usize);
            let mut value = code;
            for _ in 0..n_usize {
                row.push(value);
                value = gf.sqr(value);
            }
            row
        })
        .collect();
    let mut states: Vec<State> = Vec::new();
    let mut entries: Vec<Entry> = Vec::new();
    for left in 0..reps.len() {
        for right in 0..reps.len() {
            for relative in 0..n_usize {
                if let Some(roots) =
                    solver.roots(&gf, &basis, shifted[left][0], shifted[right][relative])
                {
                    let state = State {
                        left: left as u16,
                        right: right as u16,
                        relative: relative as u16,
                    };
                    let state_index = states.len() as u32;
                    states.push(state);
                    for (root_index, &root) in roots.iter().enumerate() {
                        let normal_root = basis.to_normal.apply(root);
                        let (key, rotation) = basis.canonical(normal_root);
                        entries.push(Entry {
                            key,
                            state: state_index,
                            root_index: root_index as u8,
                            rotation: rotation as u8,
                        });
                    }
                }
            }
        }
    }
    let scan_ms = scan_started.elapsed().as_secs_f64() * 1000.0;

    // Sort by canonical key and walk collision runs.
    let sort_started = Instant::now();
    entries.sort_unstable_by_key(|entry| entry.key);
    let sort_ms = sort_started.elapsed().as_secs_f64() * 1000.0;

    let verify_started = Instant::now();
    let mut rows_dump = std::env::var("KIC_DUMP_ROWS")
        .ok()
        .map(|path| std::fs::File::create(path).expect("create rows dump"));
    let mut collision_pairs = 0u64;
    let mut candidate_pairs = 0u64;
    let mut trivial_only_pairs = 0u64;
    let mut no_sum_pairs = 0u64;
    let mut verified_relations = 0u64;
    let mut verified_rows: Vec<Vec<u64>> = Vec::new();
    let mut echelon = Echelon {
        r,
        columns: base.columns,
        pivots: vec![None; base.columns],
        rank: 0,
    };
    let mut rank_gain_rows = 0u64;
    let mut i = 0usize;
    while i < entries.len() {
        let mut j = i + 1;
        while j < entries.len() && entries[j].key == entries[i].key {
            j += 1;
        }
        let run = &entries[i..j];
        if run.len() > 1 {
            // Pair up every two entries in the run whose pair-states differ
            // beyond the trivial transpose. Every later key-run is walked
            // too: a fixed multiple of the column count is not a finished
            // census while new rows still raise the homogeneous rank.
            for x in 0..run.len() {
                for y in (x + 1)..run.len() {
                    let first = &run[x];
                    let second = &run[y];
                    let s1 = states[first.state as usize];
                    let s2 = states[second.state as usize];
                    // Trivial: same state, or the transposed pair with the
                    // mirrored relative shift (identical point pair).
                    let same_pair =
                        (s1.left == s2.left && s1.right == s2.right && s1.relative == s2.relative)
                            || (s1.left == s2.right
                                && s1.right == s2.left
                                && (s1.relative as usize + s2.relative as usize)
                                    .is_multiple_of(n_usize));
                    if same_pair {
                        continue;
                    }
                    collision_pairs += 1;
                    match verify_collision(
                        &fast,
                        &shifted,
                        &base,
                        r,
                        n,
                        s1,
                        first.rotation,
                        s2,
                        second.rotation,
                    ) {
                        Verified::NonTrivial(row) => {
                            candidate_pairs += 1;
                            if let Some(dump) = rows_dump.as_mut() {
                                writeln!(dump, "{}", json!({"kind":"free_relation_row","row":row}))
                                    .unwrap();
                            }
                            if echelon.insert(row.clone()) {
                                rank_gain_rows += 1;
                                verified_rows.push(row);
                            }
                            verified_relations += 1;
                        }
                        Verified::TrivialOnly => trivial_only_pairs += 1,
                        Verified::NoSum => no_sum_pairs += 1,
                    }
                }
            }
        }
        i = j;
    }
    let verify_ms = verify_started.elapsed().as_secs_f64() * 1000.0;

    let mut out = std::fs::File::create(&arguments[2]).expect("create output");
    writeln!(
        out,
        "{}",
        json!({
            "kind":"koblitz_free_relation_census",
            "schema_version":1,
            "n":n, "a":a, "subgroup_order":r,
            "orbit_columns":base.columns,
            "factor_base_points":base.points.len(),
            "scanned_x":scanned_x,
            "states":states.len(),
            "root_entries":entries.len(),
            "scan_ms":scan_ms,
            "sort_ms":sort_ms,
            "collision_pairs":collision_pairs,
            "nontrivial_relations":candidate_pairs,
            "trivial_only_pairs":trivial_only_pairs,
            "no_sum_pairs":no_sum_pairs,
            "verified_relations":verified_relations,
            "rank_gain_rows":rank_gain_rows,
            "homogeneous_rank":echelon.rank,
            "rank_deficit_vs_columns":base.columns - echelon.rank,
            "verify_ms":verify_ms,
            "total_ms":process_started.elapsed().as_secs_f64() * 1000.0,
            "claim_boundary":"Stage census only: free-relation yield and homogeneous rank on the frozen factor base. No target solved, no DLP recovered, no vs_rho comparison, no ledger promotion.",
        })
    )
    .unwrap();
}

/// Verify one candidate collision as an exact group identity and return the
/// homogeneous relation row over Z/rZ, or None when no sign combination of
/// the four x-codes sums to O.
#[allow(clippy::too_many_arguments)]
fn verify_collision(
    fast: &FastBinaryCurve,
    shifted: &[Vec<u64>],
    base: &Base,
    r: u64,
    n: u32,
    s1: State,
    rotation1: u8,
    s2: State,
    rotation2: u8,
) -> Verified {
    let n_usize = n as usize;
    // The shared root in normal basis for each state, at shift 0:
    // rotate(root, rotation) == key, so root == rotate(key, n - rotation).
    // The two states' shared x-code fixes the relative alignment: state 2 must
    // be evaluated at the shift that reproduces state 1's root value.
    let rho = ((rotation1 as i32 - rotation2 as i32).rem_euclid(n as i32)) as u32;
    let shift2 = (n - rho) % n;
    let codes = [
        shifted[s1.left as usize][0],
        shifted[s1.right as usize][s1.relative as usize],
        shifted[s2.left as usize][shift2 as usize],
        shifted[s2.right as usize][(shift2 as usize + s2.relative as usize) % n_usize],
    ];
    // Every x-code must belong to the factor base; then search the sign
    // combinations for an exact sum of O.
    let choices: Vec<&Vec<usize>> = codes
        .iter()
        .map(|code| base.by_x.get(code))
        .collect::<Option<_>>()
        .unwrap_or_default();
    if choices.len() < 4 {
        return Verified::NoSum;
    }
    // Search the sign combinations for exact sums of O, preferring a choice
    // whose relation row is nonzero (a trivial pairing P + (-P) + Q + (-Q)
    // also sums to O but carries no information).
    let mut best: Option<Vec<usize>> = None;
    let debug = std::env::var("KIC_DEBUG_COLLISIONS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    let mut debug_printed = 0usize;
    'combo: for &a in choices[0] {
        for &b in choices[1] {
            let ab = fast.add(base.points[a], base.points[b]);
            for &c in choices[2] {
                let abc = fast.add(ab, base.points[c]);
                for &d in choices[3] {
                    if fast.add(abc, base.points[d]).is_none() {
                        if debug_printed < debug {
                            debug_printed += 1;
                            eprintln!(
                                "collision sum-O combo: codes={codes:?} indices=({a},{b},{c},{d}) s1=({},{},{}) rot1={} s2=({},{},{}) rot2={}",
                                s1.left, s1.right, s1.relative, rotation1,
                                s2.left, s2.right, s2.relative, rotation2
                            );
                        }
                        let mut row = vec![0u64; base.columns];
                        for &index in &[a, b, c, d] {
                            let (column, coefficient) = base.labels[index];
                            row[column] = (row[column] + coefficient) % r;
                        }
                        if row.iter().any(|&value| value != 0) {
                            best = Some(vec![a, b, c, d]);
                            break 'combo;
                        }
                        best = Some(vec![a, b, c, d]);
                    }
                }
            }
        }
    }
    let indices = match best {
        Some(indices) => indices,
        None => return Verified::NoSum,
    };
    let mut row = vec![0u64; base.columns];
    for &index in &indices {
        let (column, coefficient) = base.labels[index];
        row[column] = (row[column] + coefficient) % r;
    }
    if row.iter().all(|&value| value == 0) {
        return Verified::TrivialOnly;
    }
    Verified::NonTrivial(row)
}

enum Verified {
    NonTrivial(Vec<u64>),
    TrivialOnly,
    NoSum,
}
