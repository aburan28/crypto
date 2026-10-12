//! Exact native screen of Boolean F5 selected-row signatures.
//! Inputs are deterministic generated public Boolean polynomials only.

use crypto_lib::cryptanalysis::matrix_f5_f2::F5Criterion;
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, BufWriter, Write};
use std::time::Instant;

const PROTOCOL_BYTES: &[u8] = include_bytes!("protocol.json");
const SOURCE_BYTES: &[u8] = include_bytes!("worker.rs");
const VERIFY_BYTES: &[u8] = include_bytes!("verify.rs");
const FROZEN_PROTOCOL_SHA256: &str =
    "73a954c813732a58bb3a2ddc8f9b2d8f2c041fc4e4656c8c493fb96d8312af64";

mod verify;

#[derive(Deserialize)]
struct Protocol {
    schema_version: u32,
    degree_bound: u32,
    variables: Vec<u8>,
    batches: Vec<usize>,
    families: Vec<String>,
    discovery_seeds: Vec<u64>,
    holdout_seeds: Vec<u64>,
    worker_seconds: u64,
    selected_row_budget_per_system: usize,
    evidence_cap_bytes: usize,
    primary_variables: Vec<u8>,
    primary_batch: usize,
    minimum_largest_signature_fraction: f64,
}

#[derive(Clone)]
struct System {
    n: u8,
    quadratic: Vec<Vec<u64>>,
    // Bit zero is the constant, bit v+1 is variable v.
    affine: Vec<u64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
struct Signature {
    selected_bits_hex: String,
    selected_bits_sha256: String,
    labels: usize,
    selected: usize,
    pruned_on_grid: usize,
    reported_pruned: u64,
    koszul_pruned: u64,
    frobenius_pruned: u64,
    lower_rows: u64,
    lower_zero_reductions: u64,
    criterion_word_xors: u64,
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn fixture_with_state(n: u8, seed: u64) -> (System, u64) {
    assert!([12, 16, 20, 24].contains(&n));
    let mut state = seed;
    let mut quadratic = Vec::with_capacity(n as usize);
    for _ in 0..n {
        let mut terms = BTreeSet::new();
        while terms.len() < usize::from(2 * n) {
            let a = next(&mut state) % u64::from(n);
            let b = next(&mut state) % u64::from(n);
            if a != b {
                terms.insert((1u64 << a) | (1u64 << b));
            }
        }
        quadratic.push(terms.into_iter().collect());
    }
    (
        System {
            n,
            quadratic,
            affine: vec![0; n as usize],
        },
        state,
    )
}

fn assignments(n: u8, seed: u64, batch: usize, family: &str) -> Vec<System> {
    assert!(["independent_affine", "walk_affine"].contains(&family));
    assert!([2, 8, 32].contains(&batch));
    let (base, mut state) = fixture_with_state(n, seed);
    let mut current = base.clone();
    let mut out = Vec::with_capacity(batch);
    out.push(base.clone());
    for _ in 1..batch {
        let mut input = base.clone();
        if family == "independent_affine" {
            for affine in &mut input.affine {
                for slot in 0..=n {
                    if next(&mut state) & 1 != 0 {
                        *affine |= 1u64 << slot;
                    }
                }
            }
        } else {
            for affine in &mut current.affine {
                let slot = next(&mut state) % u64::from(n + 1);
                *affine ^= 1u64 << slot;
            }
            input.affine.clone_from(&current.affine);
        }
        out.push(input);
    }
    out
}

fn monomials_up_to_two(n: u8) -> Vec<u64> {
    let mut out =
        Vec::with_capacity(1 + usize::from(n) + (usize::from(n) * (usize::from(n) - 1) / 2));
    out.push(0);
    for a in 0..n {
        out.push(1u64 << a);
        for b in a + 1..n {
            out.push((1u64 << a) | (1u64 << b));
        }
    }
    out.sort_unstable();
    out
}

fn polynomials(system: &System) -> Vec<F2BoolPoly> {
    system
        .quadratic
        .iter()
        .zip(&system.affine)
        .map(|(quadratic, &affine)| {
            let mut terms = quadratic
                .iter()
                .copied()
                .map(F2BoolMono::from_mask)
                .collect::<Vec<_>>();
            if affine & 1 != 0 {
                terms.push(F2BoolMono::from_mask(0));
            }
            for variable in 0..system.n {
                if affine & (1u64 << (variable + 1)) != 0 {
                    terms.push(F2BoolMono::from_mask(1u64 << variable));
                }
            }
            F2BoolPoly::from_monos(terms, system.n as usize)
        })
        .collect()
}

fn signature(system: &System, multipliers: &[u64], budget: usize) -> Signature {
    let polys = polynomials(system);
    let mask = (1u64 << system.n) - 1;
    let criterion = F5Criterion::new(&polys, system.n as usize, 4, mask);
    let labels = polys.len() * multipliers.len();
    assert!(labels <= budget, "selected-row budget; censored");
    let mut bits = vec![0u8; labels.div_ceil(8)];
    let mut selected = 0;
    for generator in 0..polys.len() {
        for (slot, &multiplier) in multipliers.iter().enumerate() {
            let index = generator * multipliers.len() + slot;
            if !criterion.prunes(generator, multiplier) {
                bits[index / 8] |= 1u8 << (index % 8);
                selected += 1;
            }
        }
    }
    let (koszul_pruned, frobenius_pruned) = criterion.pruned_by_part();
    let (lower_rows, lower_zero_reductions) = criterion.lower_level_rows();
    Signature {
        selected_bits_hex: hex::encode(&bits),
        selected_bits_sha256: sha256(&bits),
        labels,
        selected,
        pruned_on_grid: labels - selected,
        reported_pruned: criterion.pruned_count(),
        koszul_pruned,
        frobenius_pruned,
        lower_rows,
        lower_zero_reductions,
        criterion_word_xors: criterion.word_ops(),
    }
}

fn cell_records(
    n: u8,
    seed: u64,
    batch: usize,
    family: &str,
    budget: usize,
    split: &str,
) -> Vec<Value> {
    let systems = assignments(n, seed, batch, family);
    let multipliers = monomials_up_to_two(n);
    let cell = format!("n{n}-{split}-{seed}-{family}-b{batch}");
    let mut outputs = Vec::with_capacity(batch + 1);
    let mut frequencies = BTreeMap::<String, usize>::new();
    let mut signatures = Vec::with_capacity(batch);
    for (index, system) in systems.iter().enumerate() {
        let record = signature(system, &multipliers, budget);
        *frequencies
            .entry(record.selected_bits_hex.clone())
            .or_default() += 1;
        signatures.push(record.selected_bits_hex.clone());
        outputs.push(json!({
            "kind":"system", "cell":cell, "index":index,
            "n":n, "seed":seed, "family":family, "batch":batch,
            "quadratic":system.quadratic,
            "affine":system.affine,
            "signature":record
        }));
    }
    let largest_class = frequencies.values().copied().max().unwrap();
    let adjacent_equal = signatures
        .windows(2)
        .filter(|pair| pair[0] == pair[1])
        .count();
    let base_equal = signatures
        .iter()
        .filter(|bits| *bits == &signatures[0])
        .count();
    outputs.push(json!({
        "kind":"cell", "cell":cell, "n":n, "seed":seed,
        "family":family, "batch":batch,
        "distinct_signatures":frequencies.len(),
        "largest_class":largest_class,
        "adjacent_equal":adjacent_equal,
        "base_equal":base_equal,
        "gate":5*largest_class >= 4*batch
    }));
    outputs
}

fn emit<W: Write>(writer: &mut W, value: &Value, bytes: &mut usize, cap: usize) {
    let line = serde_json::to_vec(value).expect("serialize row");
    *bytes += line.len() + 1;
    assert!(*bytes <= cap, "evidence cap; censored");
    writer.write_all(&line).expect("write row");
    writer.write_all(b"\n").expect("write newline");
}

fn campaign(phase: &str, protocol_path: &str) {
    assert_eq!(sha256(PROTOCOL_BYTES), FROZEN_PROTOCOL_SHA256);
    let supplied = std::fs::read(protocol_path).expect("read protocol");
    assert_eq!(supplied, PROTOCOL_BYTES, "changed protocol bytes");
    let protocol: Protocol = serde_json::from_slice(PROTOCOL_BYTES).expect("protocol JSON");
    assert_eq!(protocol.schema_version, 1);
    assert_eq!(protocol.degree_bound, 4);
    assert_eq!(protocol.variables, [12, 16, 20, 24]);
    assert_eq!(protocol.batches, [2, 8, 32]);
    assert_eq!(protocol.families, ["independent_affine", "walk_affine"]);
    assert_eq!(protocol.primary_variables, [16, 20, 24]);
    assert_eq!(protocol.primary_batch, 32);
    assert_eq!(protocol.minimum_largest_signature_fraction, 0.8);
    let (split, seeds) = match phase {
        "discovery" => ("discovery", protocol.discovery_seeds),
        "holdout" => ("holdout", protocol.holdout_seeds),
        _ => panic!("phase must be discovery or holdout"),
    };
    let mut writer = BufWriter::new(io::stdout().lock());
    let mut bytes = 0;
    emit(
        &mut writer,
        &json!({
        "kind":"campaign", "phase":phase,
        "protocol_sha256":sha256(PROTOCOL_BYTES),
        "source_sha256":sha256(SOURCE_BYTES),
        "verifier_sha256":sha256(VERIFY_BYTES)
        }),
        &mut bytes,
        protocol.evidence_cap_bytes,
    );
    let started = Instant::now();
    for n in protocol.variables {
        for &seed in &seeds {
            for family in &protocol.families {
                for &batch in &protocol.batches {
                    for record in cell_records(
                        n,
                        seed,
                        batch,
                        family,
                        protocol.selected_row_budget_per_system,
                        split,
                    ) {
                        emit(
                            &mut writer,
                            &record,
                            &mut bytes,
                            protocol.evidence_cap_bytes,
                        );
                    }
                    assert!(
                        started.elapsed().as_secs() <= protocol.worker_seconds,
                        "worker cap; incomplete run censored"
                    );
                }
            }
        }
    }
    writer.flush().expect("flush raw");
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.len() == 3 && args[1] == "--seal" {
        verify::seal(&args[2]);
        return;
    }
    if args.len() == 3 && args[1] == "--verify-bundle" {
        verify::verify_bundle(&args[2]);
        return;
    }
    if args.len() == 4 && args[1] == "--failure" {
        verify::failure(&args[2], &args[3]);
        return;
    }
    if args.len() == 4 && args[1] == "--check-discovery" {
        verify::check_discovery(&args[2], &args[3]);
        return;
    }
    if args.len() == 5 && args[1] == "--verify" {
        verify::verify_run(&args[2], &args[3], &args[4]);
        return;
    }
    if args.len() == 4 && args[1] == "--campaign" {
        campaign(&args[2], &args[3]);
        return;
    }
    panic!("usage: --campaign PHASE PROTOCOL | --verify PHASE RAW RESULTS | --check-discovery BUNDLE BINDING | --failure BUNDLE REASON | --seal BUNDLE | --verify-bundle BUNDLE");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generated_affine_families_preserve_quadratic_support() {
        for n in [12, 16, 20, 24] {
            for family in ["independent_affine", "walk_affine"] {
                let systems = assignments(n, 17, 8, family);
                assert_eq!(systems.len(), 8);
                for system in &systems {
                    assert_eq!(system.quadratic, systems[0].quadratic);
                    assert!(system
                        .quadratic
                        .iter()
                        .all(|q| q.len() == usize::from(2 * n)));
                }
                if family == "walk_affine" {
                    for pair in systems.windows(2) {
                        for (&a, &b) in pair[0].affine.iter().zip(&pair[1].affine) {
                            assert_eq!((a ^ b).count_ones(), 1);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn signature_bits_and_counts_match_public_criterion() {
        let systems = assignments(12, 19, 2, "independent_affine");
        let multipliers = monomials_up_to_two(12);
        assert_eq!(multipliers.len(), 79);
        for system in systems {
            let got = signature(&system, &multipliers, 8192);
            assert_eq!(got.labels, 12 * 79);
            assert_eq!(got.selected + got.pruned_on_grid, got.labels);
            assert_eq!(
                sha256(&hex::decode(&got.selected_bits_hex).unwrap()),
                got.selected_bits_sha256
            );
            assert!(got.lower_rows > 0);
        }
    }
}
