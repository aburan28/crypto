//! Exact structural screen for replaying Boolean Macaulay pivot traces.
//! Only deterministic generated quadratic systems are accepted.

use flate2::{read::GzDecoder, write::GzEncoder, Compression, GzBuilder};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

const PROTOCOL_BYTES: &[u8] = include_bytes!("protocol.json");
const SOURCE_BYTES: &[u8] = include_bytes!("worker.rs");
const MAGIC: &[u8] = b"BPIV1\0";
const FROZEN_PROTOCOL_SHA256: &str =
    "83b846a2ff0a364448494f9bacac31d0fe7429092654a7649ca3c097b21a0754";

#[derive(Clone, Deserialize)]
struct Protocol {
    schema_version: u32,
    variables: Vec<u8>,
    discovery_seeds: Vec<u64>,
    holdout_seeds: Vec<u64>,
    degree_bound: u8,
    row_cap: usize,
    column_cap: usize,
    screen: Screen,
}

#[derive(Clone, Deserialize)]
struct Screen {
    sizes_at_least: u8,
    minimum_fraction: f64,
    maximum_delta_rank_fraction: f64,
    require_identical_full_pivot_trace: bool,
    run_holdout_only_if_discovery_passes: bool,
}

#[derive(Clone)]
struct System {
    n: u8,
    polys: Vec<Vec<u32>>,
}

#[derive(Clone)]
struct Matrix {
    columns: Vec<u32>,
    labels: Vec<(usize, u32)>,
    rows: Vec<Vec<u64>>,
}

#[derive(Clone)]
struct Echelon {
    rows: Vec<Vec<u64>>,
    trace: Vec<Option<usize>>,
    rank: usize,
    xor_words: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
struct PairRecord {
    n: u8,
    seed: u64,
    branch_variable: u8,
    columns: usize,
    rows_zero: usize,
    rows_one: usize,
    matrix_payload_bytes_zero: usize,
    matrix_payload_bytes_one: usize,
    rank_zero: usize,
    rank_one: usize,
    rank_delta: Option<usize>,
    degree_stable: bool,
    same_pivot_columns: Option<bool>,
    same_full_trace: Option<bool>,
    replay_verified: Option<bool>,
    trace_zero: Vec<Option<usize>>,
    trace_one: Vec<Option<usize>>,
    row_digest_zero: String,
    row_digest_one: String,
    row_digest_delta: Option<String>,
    elimination_xor_words_zero: u64,
    elimination_xor_words_one: u64,
    elimination_xor_words_delta: Option<u64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
struct Summary {
    pairs: usize,
    eligible_pairs: usize,
    guarded_fallbacks: usize,
    nontrivial_eligible_large_pairs: usize,
    favorable_large_pairs: usize,
    favorable_fraction: Option<f64>,
    screen_pass: bool,
    performance_ratio: Option<f64>,
    full_solver_cost: Option<f64>,
    rho_ratio: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
struct ResultFile {
    schema_version: u32,
    phase: String,
    protocol_sha256: String,
    source_sha256: String,
    producer_binary_sha256: String,
    prior_discovery_sha256: Option<String>,
    host_os: String,
    host_arch: String,
    records: Vec<PairRecord>,
    summary: Summary,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct Manifest {
    schema_version: u32,
    files: BTreeMap<String, String>,
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn load_protocol() -> Protocol {
    assert_eq!(sha256(PROTOCOL_BYTES), FROZEN_PROTOCOL_SHA256);
    let p: Protocol = serde_json::from_slice(PROTOCOL_BYTES).expect("parse protocol");
    assert_eq!(p.schema_version, 1);
    assert_eq!(p.variables, [8, 12, 16, 20, 24]);
    assert_eq!(p.degree_bound, 3);
    assert_eq!(p.discovery_seeds.len(), 4);
    assert_eq!(p.holdout_seeds.len(), 4);
    assert_eq!(p.row_cap, 4096);
    assert_eq!(p.column_cap, 8192);
    assert_eq!(p.screen.sizes_at_least, 12);
    assert_eq!(p.screen.minimum_fraction, 0.8);
    assert_eq!(p.screen.maximum_delta_rank_fraction, 0.125);
    assert!(p.screen.require_identical_full_pivot_trace);
    assert!(p.screen.run_holdout_only_if_discovery_passes);
    p
}

fn next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn fixture(n: u8, seed: u64) -> System {
    let mut state = seed;
    let mut polys = Vec::with_capacity(n as usize);
    for _ in 0..n {
        let mut terms = BTreeSet::new();
        while terms.len() < usize::from(2 * n) {
            let a = (next(&mut state) % u64::from(n)) as u8;
            let b = (next(&mut state) % u64::from(n)) as u8;
            if a != b {
                terms.insert((1u32 << a) | (1u32 << b));
            }
        }
        polys.push(terms.into_iter().collect());
    }
    System { n, polys }
}

fn restrict(system: &System, variable: u8, value: bool) -> System {
    assert!(variable < system.n);
    let low = (1u32 << variable) - 1;
    let upper = !((1u32 << (variable + 1)) - 1);
    let bit = 1u32 << variable;
    let mut polys = Vec::with_capacity(system.polys.len());
    for poly in &system.polys {
        let mut terms = BTreeSet::new();
        for &term in poly {
            if term & bit != 0 && !value {
                continue;
            }
            let reduced = (term & low) | ((term & upper) >> 1);
            if !terms.insert(reduced) {
                terms.remove(&reduced);
            }
        }
        polys.push(terms.into_iter().collect());
    }
    System {
        n: system.n - 1,
        polys,
    }
}

fn degree(poly: &[u32]) -> Option<u8> {
    poly.iter().map(|m| m.count_ones() as u8).max()
}

fn monomials(n: u8, max_degree: u8) -> Vec<u32> {
    fn add(out: &mut Vec<u32>, n: u8, left: u8, start: u8, mask: u32) {
        if left == 0 {
            out.push(mask);
            return;
        }
        for bit in start..n {
            add(out, n, left - 1, bit + 1, mask | (1u32 << bit));
        }
    }
    let mut out = Vec::new();
    for k in 0..=max_degree.min(n) {
        add(&mut out, n, k, 0, 0);
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn labelled_matrix(system: &System, p: &Protocol) -> Matrix {
    let mut columns = monomials(system.n, p.degree_bound);
    columns.reverse();
    assert!(columns.len() <= p.column_cap, "column cap");
    let slots: BTreeMap<_, _> = columns.iter().enumerate().map(|(i, &m)| (m, i)).collect();
    let words = columns.len().div_ceil(64);
    let mut labels = Vec::new();
    let mut rows = Vec::new();
    for (generator, poly) in system.polys.iter().enumerate() {
        let Some(d) = degree(poly) else { continue };
        if d > p.degree_bound {
            continue;
        }
        for multiplier in monomials(system.n, p.degree_bound - d) {
            let mut row = vec![0u64; words];
            for &term in poly {
                let product = term | multiplier;
                let slot = slots[&product];
                row[slot / 64] ^= 1u64 << (slot % 64);
            }
            labels.push((generator, multiplier));
            rows.push(row);
            assert!(rows.len() <= p.row_cap, "row cap");
        }
    }
    Matrix {
        columns,
        labels,
        rows,
    }
}

// Independently reconstruct each product with set-parity cancellation and
// compare its monomial list with the bits of the packed constructor.
fn oracle_check(system: &System, matrix: &Matrix, p: &Protocol) {
    let mut expected_labels = Vec::new();
    let mut expected_rows = Vec::new();
    for (generator, poly) in system.polys.iter().enumerate() {
        let Some(d) = degree(poly) else { continue };
        if d > p.degree_bound {
            continue;
        }
        for multiplier in oracle_multipliers(system.n, p.degree_bound - d) {
            let mut parity = BTreeSet::new();
            for &term in poly {
                let product = term | multiplier;
                if !parity.insert(product) {
                    parity.remove(&product);
                }
            }
            expected_labels.push((generator, multiplier));
            expected_rows.push(parity);
        }
    }
    assert_eq!(matrix.labels, expected_labels, "row labels differ");
    assert_eq!(matrix.rows.len(), expected_rows.len());
    for (row, expected) in matrix.rows.iter().zip(&expected_rows) {
        let actual: BTreeSet<u32> = matrix
            .columns
            .iter()
            .enumerate()
            .filter(|(i, _)| row[i / 64] & (1u64 << (i % 64)) != 0)
            .map(|(_, &monomial)| monomial)
            .collect();
        assert_eq!(&actual, expected, "product parity differs");
    }
}

fn oracle_multipliers(n: u8, max_degree: u8) -> Vec<u32> {
    let mut levels = vec![vec![0u32]];
    for _ in 1..=max_degree {
        let mut next_level = BTreeSet::new();
        for &mask in levels.last().expect("level") {
            let start = if mask == 0 {
                0
            } else {
                32 - mask.leading_zeros()
            };
            for bit in start..u32::from(n) {
                next_level.insert(mask | (1u32 << bit));
            }
        }
        levels.push(next_level.into_iter().collect());
    }
    let mut out: Vec<u32> = levels.into_iter().flatten().collect();
    out.sort_unstable();
    out.dedup();
    out
}

fn bit(row: &[u64], column: usize) -> bool {
    row[column / 64] & (1u64 << (column % 64)) != 0
}

fn echelon(input: &Matrix) -> Echelon {
    let mut rows = input.rows.clone();
    let mut trace = Vec::with_capacity(input.columns.len());
    let mut pivot = 0;
    let mut xor_words = 0u64;
    for column in 0..input.columns.len() {
        let source = (pivot..rows.len()).find(|&r| bit(&rows[r], column));
        trace.push(source);
        if let Some(source) = source {
            rows.swap(pivot, source);
            for later in pivot + 1..rows.len() {
                if bit(&rows[later], column) {
                    for word in 0..rows[later].len() {
                        rows[later][word] ^= rows[pivot][word];
                        xor_words += 1;
                    }
                }
            }
            pivot += 1;
        }
    }
    Echelon {
        rows,
        trace,
        rank: pivot,
        xor_words,
    }
}

fn replay(input: &Matrix, trace: &[Option<usize>]) -> Option<Echelon> {
    if trace.len() != input.columns.len() {
        return None;
    }
    let mut rows = input.rows.clone();
    let mut pivot = 0;
    let mut xor_words = 0u64;
    for (column, &source) in trace.iter().enumerate() {
        match source {
            Some(source) => {
                if source < pivot || source >= rows.len() || !bit(&rows[source], column) {
                    return None;
                }
                if (pivot..source).any(|r| bit(&rows[r], column)) {
                    return None;
                }
                rows.swap(pivot, source);
                for later in pivot + 1..rows.len() {
                    if bit(&rows[later], column) {
                        for word in 0..rows[later].len() {
                            rows[later][word] ^= rows[pivot][word];
                            xor_words += 1;
                        }
                    }
                }
                pivot += 1;
            }
            None => {
                if (pivot..rows.len()).any(|r| bit(&rows[r], column)) {
                    return None;
                }
            }
        }
    }
    Some(Echelon {
        rows,
        trace: trace.to_vec(),
        rank: pivot,
        xor_words,
    })
}

fn matrix_digest(matrix: &Matrix) -> String {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(matrix.columns.len() as u64).to_le_bytes());
    bytes.extend_from_slice(&(matrix.rows.len() as u64).to_le_bytes());
    for &column in &matrix.columns {
        bytes.extend_from_slice(&column.to_le_bytes());
    }
    for row in &matrix.rows {
        for &word in row {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
    }
    sha256(&bytes)
}

fn matrix_payload_bytes(matrix: &Matrix) -> usize {
    matrix.columns.len() * std::mem::size_of::<u32>()
        + matrix.rows.len() * matrix.columns.len().div_ceil(64) * std::mem::size_of::<u64>()
}

fn encode_matrix(raw: &mut Vec<u8>, matrix: &Matrix) {
    raw.extend_from_slice(&(matrix.columns.len() as u32).to_le_bytes());
    raw.extend_from_slice(&(matrix.rows.len() as u32).to_le_bytes());
    raw.extend_from_slice(&(matrix.columns.len().div_ceil(64) as u32).to_le_bytes());
    for &column in &matrix.columns {
        raw.extend_from_slice(&column.to_le_bytes());
    }
    for row in &matrix.rows {
        for &word in row {
            raw.extend_from_slice(&word.to_le_bytes());
        }
    }
}

fn pair(n: u8, seed: u64, variable: u8, p: &Protocol, raw: &mut Vec<u8>) -> PairRecord {
    let original = fixture(n, seed);
    let zero = restrict(&original, variable, false);
    let one = restrict(&original, variable, true);
    let matrix_zero = labelled_matrix(&zero, p);
    let matrix_one = labelled_matrix(&one, p);
    oracle_check(&zero, &matrix_zero, p);
    oracle_check(&one, &matrix_one, p);
    assert_eq!(matrix_zero.columns, matrix_one.columns);
    let reduction_zero = echelon(&matrix_zero);
    let reduction_one = echelon(&matrix_one);
    raw.push(n);
    raw.extend_from_slice(&seed.to_le_bytes());
    raw.push(variable);
    encode_matrix(raw, &matrix_zero);
    encode_matrix(raw, &matrix_one);
    let stable = matrix_zero.labels == matrix_one.labels;
    let (rank_delta, same_columns, same_trace, replay_verified, delta_digest, delta_ops) = if stable
    {
        let mut delta = matrix_zero.clone();
        for (row, other) in delta.rows.iter_mut().zip(&matrix_one.rows) {
            for (left, right) in row.iter_mut().zip(other) {
                *left ^= right;
            }
        }
        let reduction_delta = echelon(&delta);
        let replayed = replay(&matrix_one, &reduction_zero.trace);
        let replay_verified = replayed
            .as_ref()
            .is_some_and(|r| r.rows == reduction_one.rows && r.rank == reduction_one.rank);
        let pivot_columns = |trace: &[Option<usize>]| {
            trace
                .iter()
                .enumerate()
                .filter_map(|(i, row)| row.map(|_| i))
                .collect::<Vec<_>>()
        };
        let same_columns =
            pivot_columns(&reduction_zero.trace) == pivot_columns(&reduction_one.trace);
        let same_trace = reduction_zero.trace == reduction_one.trace;
        assert_eq!(replay_verified, same_trace, "replay certification differs");
        (
            Some(reduction_delta.rank),
            Some(same_columns),
            Some(same_trace),
            Some(replay_verified),
            Some(matrix_digest(&delta)),
            Some(reduction_delta.xor_words),
        )
    } else {
        (None, None, None, None, None, None)
    };
    PairRecord {
        n,
        seed,
        branch_variable: variable,
        columns: matrix_zero.columns.len(),
        rows_zero: matrix_zero.rows.len(),
        rows_one: matrix_one.rows.len(),
        matrix_payload_bytes_zero: matrix_payload_bytes(&matrix_zero),
        matrix_payload_bytes_one: matrix_payload_bytes(&matrix_one),
        rank_zero: reduction_zero.rank,
        rank_one: reduction_one.rank,
        rank_delta,
        degree_stable: stable,
        same_pivot_columns: same_columns,
        same_full_trace: same_trace,
        replay_verified,
        trace_zero: reduction_zero.trace,
        trace_one: reduction_one.trace,
        row_digest_zero: matrix_digest(&matrix_zero),
        row_digest_one: matrix_digest(&matrix_one),
        row_digest_delta: delta_digest,
        elimination_xor_words_zero: reduction_zero.xor_words,
        elimination_xor_words_one: reduction_one.xor_words,
        elimination_xor_words_delta: delta_ops,
    }
}

fn summarize(records: &[PairRecord], p: &Protocol) -> Summary {
    let eligible_pairs = records.iter().filter(|r| r.degree_stable).count();
    let large: Vec<_> = records
        .iter()
        .filter(|r| {
            r.n >= p.screen.sizes_at_least
                && r.degree_stable
                && r.rank_zero > 0
                && r.rank_delta.is_some_and(|rank| rank > 0)
        })
        .collect();
    let favorable = large
        .iter()
        .filter(|r| {
            8 * r.rank_delta.unwrap() <= r.rank_zero
                && r.same_full_trace == Some(true)
                && r.replay_verified == Some(true)
        })
        .count();
    let fraction = (!large.is_empty()).then_some(favorable as f64 / large.len() as f64);
    Summary {
        pairs: records.len(),
        eligible_pairs,
        guarded_fallbacks: records.len() - eligible_pairs,
        nontrivial_eligible_large_pairs: large.len(),
        favorable_large_pairs: favorable,
        favorable_fraction: fraction,
        screen_pass: !large.is_empty() && 5 * favorable >= 4 * large.len(),
        performance_ratio: None,
        full_solver_cost: None,
        rho_ratio: None,
    }
}

fn compute(phase: &str, prior: Option<String>, binary_hash: String) -> (ResultFile, Vec<u8>) {
    let p = load_protocol();
    let seeds = match phase {
        "discovery" => &p.discovery_seeds,
        "holdout" => &p.holdout_seeds,
        _ => panic!("phase must be discovery or holdout"),
    };
    let mut records = Vec::new();
    let mut raw = MAGIC.to_vec();
    raw.extend_from_slice(&((p.variables.len() * seeds.len() * 4) as u32).to_le_bytes());
    for &n in &p.variables {
        for &seed in seeds {
            for variable in [0, n / 3, (2 * n) / 3, n - 1] {
                records.push(pair(n, seed, variable, &p, &mut raw));
            }
        }
    }
    assert_eq!(records.len(), 80);
    let summary = summarize(&records, &p);
    (
        ResultFile {
            schema_version: 1,
            phase: phase.to_string(),
            protocol_sha256: sha256(PROTOCOL_BYTES),
            source_sha256: sha256(SOURCE_BYTES),
            producer_binary_sha256: binary_hash,
            prior_discovery_sha256: prior,
            host_os: std::env::consts::OS.to_string(),
            host_arch: std::env::consts::ARCH.to_string(),
            records,
            summary,
        },
        raw,
    )
}

fn binary_hash() -> String {
    sha256(&fs::read(std::env::current_exe().expect("current executable")).expect("read binary"))
}

fn verify_result(path: &Path, raw_path: &Path) -> ResultFile {
    let result: ResultFile =
        serde_json::from_slice(&fs::read(path).expect("result")).expect("parse result");
    assert_eq!(result.schema_version, 1);
    assert_eq!(result.protocol_sha256, sha256(PROTOCOL_BYTES));
    assert_eq!(result.source_sha256, sha256(SOURCE_BYTES));
    let (mut replayed, expected_raw) = compute(
        &result.phase,
        result.prior_discovery_sha256.clone(),
        result.producer_binary_sha256.clone(),
    );
    // The machine metadata is provenance, not a field to regenerate on the
    // verifier's machine. All mathematical records are recomputed above.
    replayed.host_os.clone_from(&result.host_os);
    replayed.host_arch.clone_from(&result.host_arch);
    assert_eq!(result, replayed, "structural result replay differs");
    let mut decoder = GzDecoder::new(File::open(raw_path).expect("raw gzip"));
    let mut actual_raw = Vec::new();
    decoder.read_to_end(&mut actual_raw).expect("decode raw");
    assert_eq!(actual_raw, expected_raw, "raw matrix bytes differ");
    result
}

fn run(phase: &str, output: &Path, discovery: Option<&Path>) {
    assert!(!output.exists(), "refuse to overwrite result");
    let prior = if phase == "holdout" {
        let bundle = discovery.expect("holdout requires discovery bundle path");
        verify_bundle(bundle);
        let path = bundle.join("result.json");
        let prior_result = verify_result(&path, &bundle.join("raw.bin.gz"));
        assert_eq!(prior_result.phase, "discovery");
        assert!(prior_result.summary.screen_pass, "discovery did not pass");
        Some(sha256(&fs::read(path).expect("discovery bytes")))
    } else {
        assert!(discovery.is_none());
        None
    };
    let (result, raw) = compute(phase, prior, binary_hash());
    let raw_path = output.with_file_name("raw.bin.gz");
    assert!(!raw_path.exists(), "refuse to overwrite raw matrix corpus");
    let raw_file = File::create(&raw_path).expect("create raw");
    let mut encoder: GzEncoder<File> = GzBuilder::new()
        .mtime(0)
        .write(raw_file, Compression::default());
    encoder.write_all(&raw).expect("write raw");
    encoder.finish().expect("finish raw gzip");
    fs::write(
        output,
        serde_json::to_vec_pretty(&result).expect("serialize result"),
    )
    .expect("write result");
    verify_result(output, &raw_path);
}

fn seal(dir: &Path) {
    let manifest_path = dir.join("manifest.json");
    assert!(!manifest_path.exists(), "refuse to overwrite manifest");
    let mut files = BTreeMap::new();
    for name in [
        "Cargo.toml",
        "Cargo.lock",
        "PROTOCOL.md",
        "README.md",
        "protocol.json",
        "worker.rs",
        "run.sh",
        "boolean-pivot-screen",
        "result.json",
        "raw.bin.gz",
        "rustc.txt",
        "test.stdout",
        "test.stderr",
        "build.stdout",
        "build.stderr",
        "worker.stdout",
        "worker.stderr",
        "git-head.txt",
        "host.txt",
    ] {
        files.insert(
            name.to_string(),
            sha256(&fs::read(dir.join(name)).expect(name)),
        );
    }
    fs::write(
        manifest_path,
        serde_json::to_vec_pretty(&Manifest {
            schema_version: 1,
            files,
        })
        .expect("serialize manifest"),
    )
    .expect("write manifest");
}

fn verify_bundle(dir: &Path) {
    let manifest: Manifest =
        serde_json::from_slice(&fs::read(dir.join("manifest.json")).expect("manifest"))
            .expect("parse manifest");
    assert_eq!(manifest.schema_version, 1);
    assert_eq!(manifest.files.len(), 19);
    for (name, digest) in &manifest.files {
        assert_eq!(
            &sha256(&fs::read(dir.join(name)).expect("bundle member")),
            digest
        );
    }
    assert_eq!(fs::read(dir.join("protocol.json")).unwrap(), PROTOCOL_BYTES);
    assert_eq!(fs::read(dir.join("worker.rs")).unwrap(), SOURCE_BYTES);
    let result = verify_result(&dir.join("result.json"), &dir.join("raw.bin.gz"));
    assert_eq!(
        result.producer_binary_sha256,
        sha256(&fs::read(dir.join("boolean-pivot-screen")).unwrap())
    );
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    match args.as_slice() {
        [_, command, phase, output] if command == "--run" => {
            run(phase, &PathBuf::from(output), None)
        }
        [_, command, phase, output, discovery] if command == "--run" => {
            run(phase, &PathBuf::from(output), Some(Path::new(discovery)))
        }
        [_, command, dir] if command == "--seal" => seal(Path::new(dir)),
        [_, command, dir] if command == "--verify-bundle" => verify_bundle(Path::new(dir)),
        _ => panic!("usage: worker --run discovery OUT | --run holdout OUT DISCOVERY | --seal DIR | --verify-bundle DIR"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn independent_sparse_rank(matrix: &Matrix) -> usize {
        let mut pivots: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
        for row in &matrix.rows {
            let mut support: BTreeSet<u32> = matrix
                .columns
                .iter()
                .enumerate()
                .filter(|(i, _)| bit(row, *i))
                .map(|(_, &column)| column)
                .collect();
            while let Some(&lead) = support.iter().next_back() {
                if let Some(pivot) = pivots.get(&lead) {
                    for &term in pivot {
                        if !support.insert(term) {
                            support.remove(&term);
                        }
                    }
                } else {
                    pivots.insert(lead, support);
                    break;
                }
            }
        }
        pivots.len()
    }

    #[test]
    fn restriction_derivative_is_affine() {
        for n in [8, 12, 16, 20, 24] {
            let source = fixture(n, 20261003);
            for variable in [0, n / 3, (2 * n) / 3, n - 1] {
                let zero = restrict(&source, variable, false);
                let one = restrict(&source, variable, true);
                for (a, b) in zero.polys.iter().zip(&one.polys) {
                    let mut derivative = BTreeSet::new();
                    for &term in a.iter().chain(b) {
                        if !derivative.insert(term) {
                            derivative.remove(&term);
                        }
                    }
                    assert!(derivative.iter().all(|m| m.count_ones() <= 1));
                }
            }
        }
    }

    #[test]
    fn multipliers_agree_with_independent_enumerator() {
        for n in 1..=12 {
            for d in 0..=3 {
                assert_eq!(monomials(n, d), oracle_multipliers(n, d));
            }
        }
    }

    #[test]
    fn replay_requires_the_exact_schedule() {
        let p = load_protocol();
        let original = fixture(8, 17);
        let matrix = labelled_matrix(&restrict(&original, 0, false), &p);
        oracle_check(&restrict(&original, 0, false), &matrix, &p);
        let fresh = echelon(&matrix);
        let played = replay(&matrix, &fresh.trace).expect("self replay");
        assert_eq!(played.rows, fresh.rows);
        let mut broken = fresh.trace.clone();
        let pivot = broken.iter().position(Option::is_some).unwrap();
        broken[pivot] = None;
        assert!(replay(&matrix, &broken).is_none());
    }

    #[test]
    fn independent_sparse_rank_agrees_on_paired_matrices_and_deltas() {
        let p = load_protocol();
        for seed in [17, 20261003, 3141593, 2718281] {
            let source = fixture(8, seed);
            for variable in [0, 2, 5, 7] {
                let zero = labelled_matrix(&restrict(&source, variable, false), &p);
                let one = labelled_matrix(&restrict(&source, variable, true), &p);
                oracle_check(&restrict(&source, variable, false), &zero, &p);
                oracle_check(&restrict(&source, variable, true), &one, &p);
                assert_eq!(echelon(&zero).rank, independent_sparse_rank(&zero));
                assert_eq!(echelon(&one).rank, independent_sparse_rank(&one));
                if zero.labels == one.labels {
                    let mut delta = zero.clone();
                    for (left, right) in delta.rows.iter_mut().zip(&one.rows) {
                        for (a, b) in left.iter_mut().zip(right) {
                            *a ^= b;
                        }
                    }
                    assert_eq!(echelon(&delta).rank, independent_sparse_rank(&delta));
                }
            }
        }
    }

    #[test]
    fn degree_drop_changes_labelled_skeleton() {
        let p = load_protocol();
        let source = System {
            n: 4,
            polys: vec![vec![0b0011]],
        };
        let zero = restrict(&source, 0, false);
        let one = restrict(&source, 0, true);
        let m0 = labelled_matrix(&zero, &p);
        let m1 = labelled_matrix(&one, &p);
        oracle_check(&zero, &m0, &p);
        oracle_check(&one, &m1, &p);
        assert_ne!(m0.labels, m1.labels);
    }
}
