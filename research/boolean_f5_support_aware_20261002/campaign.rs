//! Frozen same-binary producer and replay for support-aware Boolean F5 batches.
//! All fixtures are generated public Boolean systems; no curve or key input exists.

use super::*;
use crypto_lib::cryptanalysis::koblitz_groebner::matrix_f4_f2;
use crypto_lib::cryptanalysis::matrix_f5_f2::{F5Report, F5Timings};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::fs;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;
use std::thread;
use std::time::{Duration, Instant};

const PROTOCOL_BYTES: &[u8] = include_bytes!("protocol.json");
const WORKER_BYTES: &[u8] = include_bytes!("worker.rs");
const CAMPAIGN_BYTES: &[u8] = include_bytes!("campaign.rs");
const GF2_BYTES: &[u8] = include_bytes!("../../src/cryptanalysis/gf2_elim.rs");
const FROZEN_PROTOCOL_SHA256: &str =
    "317118981ba789248225609f92aa1bb6ae0834a97a9da9c9ce321a157676546a";

#[derive(Deserialize)]
struct Protocol {
    schema_version: u32,
    study: String,
    parent_discovery_run_id: u64,
    parent_source_head: String,
    parent_worker_sha256: String,
    parent_manifest_sha256: String,
    degree_bound: u32,
    output_form: String,
    high_prefix_rule: String,
    support_projection_guard: String,
    metric: String,
    arms: Vec<String>,
    primary_candidate: String,
    env: BTreeMap<String, String>,
    variables: Vec<u8>,
    batches: Vec<usize>,
    families: Vec<String>,
    discovery_seeds: Vec<u64>,
    holdout_seeds: Vec<u64>,
    repetitions: usize,
    bootstrap_resamples: usize,
    bootstrap_seed: u64,
    primary_variables: Vec<u8>,
    primary_batch: usize,
    primary_groups: usize,
    primary_median_greater_than: f64,
    primary_lower_bound_greater_than: f64,
    structural_primary_sparse_recovery_required: bool,
    support_vs_anchor_lower_bound_above_aa_floor: bool,
    nonregression_variables: Vec<u8>,
    nonregression_lower_bound_greater_than: f64,
    aa_noise_quantile: f64,
    retained_context_cap_bytes: usize,
    worker_seconds: u64,
    evidence_cap_bytes: usize,
    maximum_other_cpu_fraction: f64,
    maximum_psi_some_avg10: f64,
    full_ic_cost: Option<f64>,
    rho_ratio: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
struct Reference {
    report: F5Report,
    output_digest: String,
    returned_rows: usize,
    returned_terms: usize,
    direct_pack_used: bool,
    direct_unpack_used: bool,
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn route() -> Protocol {
    assert_eq!(sha256(PROTOCOL_BYTES), FROZEN_PROTOCOL_SHA256);
    let p: Protocol = serde_json::from_slice(PROTOCOL_BYTES).expect("protocol");
    assert_eq!(p.schema_version, 1);
    assert_eq!(p.study, "support_aware_exact_boolean_f5_batch");
    assert_eq!(p.parent_discovery_run_id, 37045285547);
    assert_eq!(
        p.parent_source_head,
        "53f26c669852108c8128c16341f722a1c49e4f5d"
    );
    assert_eq!(
        p.parent_worker_sha256,
        "a37fc763173eb45d0542c23a24eb9074287882584a6ba11c4566b053ec5ce3e6"
    );
    assert_eq!(
        p.parent_manifest_sha256,
        "651b2ae2932bcebd2bf14d012c41477671c3b69e9e320393762a90a2eececbe5"
    );
    assert_eq!(p.degree_bound, 4);
    assert_eq!(p.output_form, "Echelon");
    assert_eq!(p.high_prefix_rule, "floor(binomial(n,4)/64)*64");
    assert_eq!(
        p.support_projection_guard,
        "all_high_prefix_columns_present_and_exact_selected_row_labels"
    );
    assert_eq!(p.metric, "complete_cold_batch_wall_ns");
    assert_eq!(
        p.arms,
        [
            "aa_a",
            "aa_b",
            "fresh",
            "anchor",
            "support",
            "support_unpack",
            "support_unpack_layout"
        ]
    );
    assert_eq!(p.primary_candidate, "support_unpack_layout");
    assert_eq!(p.variables, [12, 16, 20, 24]);
    assert_eq!(p.batches, [2, 8, 32]);
    assert_eq!(p.families, ["independent_affine", "walk_affine"]);
    assert_eq!(p.discovery_seeds, [20261103, 3141753]);
    assert_eq!(p.holdout_seeds, [20261110, 4242601]);
    assert_eq!(p.repetitions, 7);
    assert_eq!(p.bootstrap_resamples, 4000);
    assert_eq!(p.primary_variables, [24]);
    assert_eq!(p.primary_batch, 32);
    assert_eq!(p.primary_groups, 4);
    assert_eq!(p.primary_median_greater_than, 2.0);
    assert_eq!(p.primary_lower_bound_greater_than, 2.0);
    assert!(p.structural_primary_sparse_recovery_required);
    assert!(p.support_vs_anchor_lower_bound_above_aa_floor);
    assert_eq!(p.nonregression_variables, [16, 20]);
    assert_eq!(p.nonregression_lower_bound_greater_than, 0.95);
    assert_eq!(p.aa_noise_quantile, 0.975);
    assert_eq!(p.retained_context_cap_bytes, 128 * 1024 * 1024);
    assert_eq!(p.worker_seconds, 2100);
    assert_eq!(p.evidence_cap_bytes, 64 * 1024 * 1024);
    assert_eq!(p.maximum_other_cpu_fraction, 0.1);
    assert_eq!(p.maximum_psi_some_avg10, 5.0);
    assert!(p.full_ic_cost.is_none() && p.rho_ratio.is_none());
    p
}

fn check_environment(expected: &BTreeMap<String, String>) {
    for (name, value) in expected {
        assert_eq!(std::env::var(name).as_deref(), Ok(value.as_str()), "{name}");
    }
    for (name, _) in std::env::vars() {
        if (name.starts_with("KIC_F5_") || name.starts_with("KIC_GF2_"))
            && !expected.contains_key(&name)
        {
            panic!("unregistered F5/GF2 option {name}");
        }
    }
}

fn check_host() {
    #[cfg(target_arch = "x86_64")]
    {
        assert_eq!(std::env::consts::OS, "linux");
        assert!(std::arch::is_x86_feature_detected!("avx2"));
    }
    #[cfg(not(target_arch = "x86_64"))]
    panic!("qualified campaign requires Linux x86-64 AVX2");
}

fn assignments(n: u8, seed: u64, batch: usize, family: &str) -> Vec<System> {
    assert!([2, 8, 32].contains(&batch));
    assert!(["independent_affine", "walk_affine"].contains(&family));
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
    let base = System {
        n,
        quadratic,
        affine: vec![0; n as usize],
    };
    let mut current = base.clone();
    let mut out = vec![base.clone()];
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

fn output_digest(rows: &[F2BoolPoly]) -> (String, usize) {
    let mut hasher = Sha256::new();
    hasher.update(b"boolean-f5-returned-rows-v1");
    hasher.update((rows.len() as u64).to_le_bytes());
    let mut terms = 0;
    for row in rows {
        hasher.update((row.n_vars as u64).to_le_bytes());
        hasher.update((row.terms.len() as u64).to_le_bytes());
        for mono in &row.terms {
            hasher.update(mono.mask.to_le_bytes());
            terms += 1;
        }
    }
    (format!("{:x}", hasher.finalize()), terms)
}

fn call_f5(polys: &[F2BoolPoly], n: u8) -> (Vec<F2BoolPoly>, F5Report, F5Timings) {
    matrix_f5_f2_with_form_timed(polys, n as usize, 4, F5OutputForm::Echelon)
        .expect("complete F5 step")
}

fn reference(polys: &[F2BoolPoly], n: u8) -> Reference {
    let (rows, report, timings) = call_f5(polys, n);
    let (output_digest, returned_terms) = output_digest(&rows);
    assert!(timings.direct_unpack_used);
    Reference {
        report,
        output_digest,
        returned_rows: rows.len(),
        returned_terms,
        direct_pack_used: timings.direct_pack_used,
        direct_unpack_used: timings.direct_unpack_used,
    }
}

fn f4_cross_check(polys: &[F2BoolPoly], n: u8) -> Option<u64> {
    if n > 16 {
        return None;
    }
    let f4 = matrix_f4_f2(polys, n as usize, 4).expect("small F4");
    let (f5, _, _) = call_f5(polys, n);
    let f4_fingerprint = canonical_row_space_fingerprint(&f4).expect("F4 fingerprint");
    assert_eq!(
        f4_fingerprint,
        canonical_row_space_fingerprint(&f5).unwrap()
    );
    Some(f4_fingerprint)
}

fn arm_order(repetition: usize) -> [&'static str; 7] {
    if repetition % 2 == 0 {
        [
            "aa_a",
            "aa_b",
            "fresh",
            "anchor",
            "support",
            "support_unpack",
            "support_unpack_layout",
        ]
    } else {
        [
            "support_unpack_layout",
            "support_unpack",
            "support",
            "anchor",
            "fresh",
            "aa_b",
            "aa_a",
        ]
    }
}

fn sparse_eligibility(systems: &[System]) -> Vec<bool> {
    let basis = Basis::new(systems[0].n);
    let first = selected_matrix(&systems[0], &basis);
    let prefix_word = basis.high_cols / 64;
    systems
        .iter()
        .map(|system| {
            let selected = selected_matrix(system, &basis);
            !selected.full_columns
                && selected.labels == first.labels
                && selected.occupied_columns[..prefix_word]
                    .iter()
                    .all(|&word| word == u64::MAX)
                && selected.rows.len() == first.rows.len()
                && selected
                    .rows
                    .iter()
                    .zip(&first.rows)
                    .all(|(row, base)| row[..prefix_word] == base[..prefix_word])
        })
        .collect()
}

enum CacheContext {
    High(HighCache),
    Layout(LayoutCache),
}

impl CacheContext {
    fn compile(arm: &str, base: &System) -> Result<Self, &'static str> {
        if arm == "support_unpack_layout" {
            LayoutCache::compile(base).map(Self::Layout)
        } else {
            HighCache::compile(base).map(Self::High)
        }
    }

    fn context_bytes(&self) -> usize {
        match self {
            Self::High(cache) => cache.context_bytes,
            Self::Layout(cache) => cache.context_bytes,
        }
    }

    fn compile_word_xors(&self) -> u64 {
        match self {
            Self::High(cache) => cache.compile_word_xors,
            Self::Layout(cache) => cache.high.compile_word_xors,
        }
    }

    fn apply(
        &mut self,
        arm: &str,
        input: &System,
    ) -> (Vec<F2BoolPoly>, bool, &'static str, u64, u64, bool) {
        match self {
            Self::High(cache) if arm == "anchor" => {
                let (rows, hit, reason, ops, tail) = cache.apply(input);
                (rows, hit, reason, ops, tail, false)
            }
            Self::High(cache) if arm == "support" => cache.apply_support(input, false),
            Self::High(cache) if arm == "support_unpack" => cache.apply_support(input, true),
            Self::Layout(cache) if arm == "support_unpack_layout" => cache.apply(input),
            _ => panic!("cache context/arm mismatch"),
        }
    }
}

fn timed_batch(
    systems: &[System],
    polys: &[Vec<F2BoolPoly>],
    references: &[Reference],
    repetition: usize,
    arm: &str,
) -> Value {
    let n = systems[0].n;
    let batch = systems.len();
    let candidate = !matches!(arm, "fresh" | "aa_a" | "aa_b");
    let order = (0..batch)
        .map(|position| (position + repetition) % batch)
        .collect::<Vec<_>>();
    let mut total_ns = 0u128;
    let mut validation_ns = 0u128;
    let mut compile_ns = 0u128;
    let mut context_bytes = 0usize;
    let mut compile_word_xors = 0u64;
    let mut compile_error = None;
    let mut cache = if candidate {
        let start = Instant::now();
        let result = CacheContext::compile(arm, &systems[order[0]]);
        compile_ns = start.elapsed().as_nanos();
        total_ns += compile_ns;
        match result {
            Ok(cache) => {
                context_bytes = cache.context_bytes();
                compile_word_xors = cache.compile_word_xors();
                assert!(context_bytes <= 128 * 1024 * 1024);
                Some(cache)
            }
            Err(reason) => {
                assert_ne!(reason, "context-cap", "resource cap is censored");
                compile_error = Some(reason);
                None
            }
        }
    } else {
        None
    };
    let mut calls = Vec::with_capacity(batch);
    for &index in &order {
        let start = Instant::now();
        let (rows, report, timing, hit, reason, candidate_word_xors, tail_word_xors, projected) =
            if candidate {
                let (rows, hit, reason, ops, tail, projected) = if let Some(ref mut cache) = cache {
                    cache.apply(arm, &systems[index])
                } else {
                    (
                        baseline(&systems[index]),
                        false,
                        compile_error.unwrap(),
                        0,
                        0,
                        false,
                    )
                };
                (rows, None, None, hit, reason, ops, tail, projected)
            } else {
                let (rows, report, timing) = call_f5(&polys[index], n);
                (
                    rows,
                    Some(report),
                    Some(timing),
                    false,
                    "baseline",
                    0,
                    0,
                    false,
                )
            };
        let call_ns = start.elapsed().as_nanos();
        total_ns += call_ns;
        let validation_start = Instant::now();
        let expected = &references[index];
        let (digest, terms) = output_digest(&rows);
        assert_eq!(digest, expected.output_digest);
        assert_eq!(rows.len(), expected.returned_rows);
        assert_eq!(terms, expected.returned_terms);
        // Keep validation work matched between arms. Every timed return is
        // pinned by its ordered-output digest; the full Vec comparison is
        // performed once per unique public fixture before the timed batches.
        let selected = selected_matrix(&systems[index], &Basis::new(n));
        assert_eq!(selected.rows_f4, expected.report.rows_f4);
        assert_eq!(selected.rows.len() as u64, expected.report.rows_built);
        assert_eq!(selected.criterion_rows, expected.report.criterion_rows);
        assert_eq!(
            selected.criterion_word_xors,
            expected.report.criterion_word_ops
        );
        let (selected_rows, rows_f4, criterion_rows, criterion_word_xors, transform_word_xors) =
            if candidate {
                let transform = if hit {
                    assert!(candidate_word_xors >= selected.criterion_word_xors + tail_word_xors);
                    candidate_word_xors - selected.criterion_word_xors - tail_word_xors
                } else {
                    0
                };
                (
                    selected.rows.len() as u64,
                    selected.rows_f4,
                    selected.criterion_rows,
                    selected.criterion_word_xors,
                    transform,
                )
            } else {
                assert_eq!(report.unwrap(), expected.report);
                let timing = timing.unwrap();
                assert_eq!(timing.direct_pack_used, expected.direct_pack_used);
                assert_eq!(timing.direct_unpack_used, expected.direct_unpack_used);
                (0, 0, 0, 0, 0)
            };
        let check_ns = validation_start.elapsed().as_nanos();
        validation_ns += check_ns;
        let drop_start = Instant::now();
        drop(rows);
        let destruction_ns = drop_start.elapsed().as_nanos();
        total_ns += destruction_ns;
        calls.push(json!({
            "input_index":index,"call_ns":call_ns,"destruction_ns":destruction_ns,
            "validation_ns":check_ns,"output_digest":digest,
            "returned_rows":expected.returned_rows,"returned_terms":terms,
            "report":report,"cache_hit":hit,"fallback_reason":reason,
            "support_projected":projected,
            "candidate_word_xors":candidate_word_xors,
            "candidate_transform_word_xors":transform_word_xors,
            "candidate_tail_word_xors":tail_word_xors,
            "selected_rows":selected_rows,"rows_f4":rows_f4,
            "criterion_rows":criterion_rows,"criterion_word_xors":criterion_word_xors
        }));
    }
    let drop_start = Instant::now();
    drop(cache);
    let cache_drop_ns = drop_start.elapsed().as_nanos();
    total_ns += cache_drop_ns;
    json!({
        "kind":"sample","repetition":repetition,"arm":arm,
        "total_ns":total_ns,"compile_ns":compile_ns,
        "cache_drop_ns":cache_drop_ns,"compile_word_xors":compile_word_xors,
        "context_bytes":context_bytes,"compile_error":compile_error,
        "validation_ns":validation_ns,"calls":calls
    })
}

fn cell(split: &str, n: u8, seed: u64, family: &str, batch: usize, repetitions: usize) {
    let systems = assignments(n, seed, batch, family);
    let polys = systems.iter().map(polynomials).collect::<Vec<_>>();
    let references = polys.iter().map(|p| reference(p, n)).collect::<Vec<_>>();
    let eligible_sparse = sparse_eligibility(&systems);
    let mut preflight_caches = [
        "anchor",
        "support",
        "support_unpack",
        "support_unpack_layout",
    ]
    .into_iter()
    .map(|arm| {
        let cache = match CacheContext::compile(arm, &systems[0]) {
            Ok(cache) => Some(cache),
            Err(reason) => {
                assert_ne!(reason, "context-cap", "resource cap is censored");
                None
            }
        };
        (arm, cache)
    })
    .collect::<Vec<_>>();
    for (index, system) in systems.iter().enumerate() {
        let (expected, report, _) = call_f5(&polys[index], n);
        assert_eq!(report, references[index].report);
        for (arm, cache) in &mut preflight_caches {
            let actual = if let Some(cache) = cache {
                cache.apply(arm, system).0
            } else {
                baseline(system)
            };
            assert_eq!(
                actual, expected,
                "direct preflight {arm} ordered output {index}"
            );
        }
    }
    drop(preflight_caches);
    let f4_fingerprint = f4_cross_check(&polys[0], n);
    let cell_id = format!("n{n}-{split}-{seed}-{family}-b{batch}");
    println!(
        "{}",
        json!({
            "kind":"fixture","cell":cell_id,"n":n,"seed":seed,
            "family":family,"batch":batch,"quadratic":systems[0].quadratic,
            "affine":systems.iter().map(|s| &s.affine).collect::<Vec<_>>(),
            "references":references,"small_f4_fingerprint":f4_fingerprint,
            "eligible_sparse":eligible_sparse,
            "direct_preflight":true
        })
    );
    for repetition in 0..repetitions {
        for arm in arm_order(repetition) {
            let mut sample = timed_batch(&systems, &polys, &references, repetition, arm);
            sample["cell"] = json!(cell_id);
            println!("{sample}");
        }
    }
}

fn campaign(phase: &str, protocol_path: &str) {
    let p = route();
    check_host();
    check_environment(&p.env);
    assert_eq!(
        fs::read(protocol_path).expect("protocol file"),
        PROTOCOL_BYTES
    );
    let seeds = match phase {
        "discovery" => &p.discovery_seeds,
        "holdout" => &p.holdout_seeds,
        _ => panic!("phase"),
    };
    let header = json!({
        "kind":"campaign","phase":phase,"protocol_sha256":sha256(PROTOCOL_BYTES),
        "worker_sha256":sha256(WORKER_BYTES),"campaign_sha256":sha256(CAMPAIGN_BYTES),
        "gf2_sha256":sha256(GF2_BYTES),"host_arch":"x86_64","route":p.env
    });
    let mut writer = BufWriter::new(io::stdout().lock());
    writeln!(writer, "{header}").unwrap();
    let mut bytes = header.to_string().len() + 1;
    let started = Instant::now();
    let executable = std::env::current_exe().expect("worker path");
    for &n in &p.variables {
        for &seed in seeds {
            for family in &p.families {
                for &batch in &p.batches {
                    let output = Command::new(&executable)
                        .args([
                            "--cell",
                            phase,
                            &n.to_string(),
                            &seed.to_string(),
                            family,
                            &batch.to_string(),
                            &p.repetitions.to_string(),
                        ])
                        .output()
                        .expect("cell process");
                    bytes += output.stdout.len();
                    assert!(bytes <= p.evidence_cap_bytes, "evidence cap; censored");
                    writer.write_all(&output.stdout).expect("raw records");
                    writer.flush().expect("flush cell");
                    assert!(
                        output.status.success(),
                        "cell failed: {}",
                        String::from_utf8_lossy(&output.stderr)
                    );
                    assert!(
                        started.elapsed().as_secs() <= p.worker_seconds,
                        "worker cap; censored"
                    );
                }
            }
        }
    }
}

pub(super) fn main_cli() {
    let args = std::env::args().collect::<Vec<_>>();
    match args
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>()
        .as_slice()
    {
        [_, "--cell", split, n, seed, family, batch, repetitions] => {
            let p = route();
            check_environment(&p.env);
            if *split != "development" {
                check_host();
            }
            cell(
                split,
                n.parse().unwrap(),
                seed.parse().unwrap(),
                family,
                batch.parse().unwrap(),
                repetitions.parse().unwrap(),
            );
        }
        [_, "--campaign", phase, path] => campaign(phase, path),
        [_, "--verify", phase, raw, results] => verify_run(phase, raw, results),
        [_, "--wait-quiet", path] => wait_quiet(path),
        [_, "--seal", path] => seal(path),
        [_, "--verify-bundle", path] => {
            verify_bundle(path);
        }
        [_, "--check-discovery", bundle, binding] => check_discovery(bundle, binding),
        [_, "--failure", dir, reason] => failure(dir, reason),
        _ => panic!("usage: --cell | --campaign | --verify | --seal | --verify-bundle"),
    }
}

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("read JSON")).expect("parse JSON")
}

fn raw_lines(path: &Path) -> Vec<Value> {
    BufReader::new(fs::File::open(path).expect("raw file"))
        .lines()
        .map(|line| serde_json::from_str(&line.expect("raw line")).expect("raw JSON"))
        .collect()
}

fn integer(value: &Value, field: &str) -> u64 {
    value[field]
        .as_u64()
        .unwrap_or_else(|| panic!("numeric {field}"))
}

fn median(values: &[f64]) -> f64 {
    assert!(!values.is_empty());
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let middle = sorted.len() / 2;
    if sorted.len().is_multiple_of(2) {
        (sorted[middle - 1] + sorted[middle]) / 2.0
    } else {
        sorted[middle]
    }
}

fn quantile(values: &[f64], fraction: f64) -> f64 {
    assert!(!values.is_empty() && (0.0..=1.0).contains(&fraction));
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let index = ((sorted.len() as f64 * fraction).ceil() as usize)
        .saturating_sub(1)
        .min(sorted.len() - 1);
    sorted[index]
}

fn bootstrap_lower(values: &[f64], seed: u64, resamples: usize) -> f64 {
    let mut state = seed;
    let mut medians = Vec::with_capacity(resamples);
    for _ in 0..resamples {
        let mut sample = Vec::with_capacity(values.len());
        for _ in values {
            sample.push(values[(next(&mut state) % values.len() as u64) as usize]);
        }
        medians.push(median(&sample));
    }
    quantile(&medians, 0.025)
}

fn parse_cpu_list(list: &str) -> Vec<u64> {
    let mut out = BTreeSet::new();
    for part in list.trim().split(',') {
        if let Some((start, end)) = part.split_once('-') {
            let start = start.parse::<u64>().expect("range start");
            let end = end.parse::<u64>().expect("range end");
            assert!(start <= end);
            out.extend(start..=end);
        } else {
            out.insert(part.parse().expect("CPU number"));
        }
    }
    assert!(!out.is_empty());
    out.into_iter().collect()
}

fn resource(path: &Path, phase: &str, p: &Protocol) -> Value {
    let lines = raw_lines(path);
    assert_eq!(lines.len(), 1);
    let receipt = &lines[0];
    assert_eq!(receipt["schema"], "isolated-bench/1");
    assert_eq!(receipt["mode"], "run");
    assert_eq!(receipt["label"], format!("f5-support-aware/{phase}"));
    let siblings =
        parse_cpu_list(&fs::read_to_string(path.with_file_name("cpu-siblings.txt")).unwrap());
    let reserved = receipt["reserved_cpus"]
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_u64().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(reserved, siblings);
    assert_eq!(receipt["run"]["contended"], false);
    assert_eq!(receipt["run"]["exit_status"], 0);
    assert_eq!(receipt["host"]["machine"], "x86_64");
    let settle = receipt["preflight"]["settle"]["seconds"].as_f64().unwrap();
    let settle_other = receipt["preflight"]["settle"]["other_cpu_seconds"]
        .as_f64()
        .unwrap();
    assert!(settle > 0.0 && settle_other / settle <= p.maximum_other_cpu_fraction);
    let psi = receipt["preflight"]["conditions"]["psi_cpu"]["some"]["avg10"]
        .as_f64()
        .unwrap();
    assert!(psi <= p.maximum_psi_some_avg10);
    let wall = receipt["run"]["wall_seconds"].as_f64().unwrap();
    let other = receipt["run"]["other_cpu_seconds"].as_f64().unwrap();
    assert!(wall > 0.0 && wall <= p.worker_seconds as f64);
    assert!(other / wall <= p.maximum_other_cpu_fraction);
    json!({"wall_seconds":wall,"other_cpu_seconds":other,"psi_some_avg10":psi,
        "max_rss_kib":integer(&receipt["run"],"max_rss_kib"),
        "reserved_cpus":receipt["reserved_cpus"],"contended":false})
}

fn readiness(path: &Path, p: &Protocol) -> Value {
    let record = read_json(path);
    assert_eq!(record["accepted"], true);
    let samples = record["samples"].as_array().expect("readiness samples");
    assert!(!samples.is_empty() && samples.len() <= 30);
    for (index, sample) in samples.iter().enumerate() {
        assert_eq!(integer(sample, "index"), index as u64);
        let psi = sample["psi_some_avg10"].as_f64().unwrap();
        assert_eq!(sample["accepted"], psi <= p.maximum_psi_some_avg10);
    }
    assert_eq!(samples.last().unwrap()["accepted"], true);
    record
}

fn verify_sample(
    line: &Value,
    cell_id: &str,
    repetition: usize,
    arm: &str,
    references: &[Reference],
    systems: &[System],
    eligible_sparse: &[bool],
    p: &Protocol,
) -> (f64, usize, usize, u128, usize, usize) {
    assert_eq!(line["kind"], "sample");
    assert_eq!(line["cell"], cell_id);
    assert_eq!(integer(line, "repetition"), repetition as u64);
    assert_eq!(line["arm"], arm);
    let candidate = !matches!(arm, "fresh" | "aa_a" | "aa_b");
    let support_arm = matches!(arm, "support" | "support_unpack" | "support_unpack_layout");
    let calls = line["calls"].as_array().expect("calls");
    assert_eq!(calls.len(), references.len());
    assert_eq!(eligible_sparse.len(), calls.len());
    let compile_ns = integer(line, "compile_ns") as u128;
    let cache_drop_ns = integer(line, "cache_drop_ns") as u128;
    let context_bytes = integer(line, "context_bytes") as usize;
    assert!(context_bytes <= p.retained_context_cap_bytes);
    if !candidate {
        assert_eq!(compile_ns, 0);
        assert_eq!(context_bytes, 0);
        assert_eq!(integer(line, "compile_word_xors"), 0);
        assert!(line["compile_error"].is_null());
    }
    let mut total = compile_ns + cache_drop_ns;
    let mut validation = 0u128;
    let mut hits = 0;
    let mut fallbacks = 0;
    let mut eligible = 0;
    let mut recovered = 0;
    let mut word_xors = if candidate {
        integer(line, "compile_word_xors") as u128
    } else {
        0
    };
    for (position, call) in calls.iter().enumerate() {
        let index = (position + repetition) % references.len();
        assert_eq!(integer(call, "input_index"), index as u64);
        let expected = &references[index];
        let projected = call["support_projected"]
            .as_bool()
            .expect("support_projected");
        if support_arm && eligible_sparse[index] {
            eligible += 1;
            recovered += usize::from(projected);
        }
        assert_eq!(call["output_digest"], expected.output_digest);
        assert_eq!(
            integer(call, "returned_rows"),
            expected.returned_rows as u64
        );
        assert_eq!(
            integer(call, "returned_terms"),
            expected.returned_terms as u64
        );
        let call_ns = integer(call, "call_ns") as u128;
        let destruction_ns = integer(call, "destruction_ns") as u128;
        let check_ns = integer(call, "validation_ns") as u128;
        assert!(call_ns > 0 && check_ns > 0);
        total += call_ns + destruction_ns;
        validation += check_ns;
        if candidate {
            assert!(call["report"].is_null());
            let selected = selected_matrix(&systems[index], &Basis::new(systems[index].n));
            assert_eq!(selected.rows_f4, expected.report.rows_f4);
            assert_eq!(selected.rows.len() as u64, expected.report.rows_built);
            assert_eq!(selected.criterion_rows, expected.report.criterion_rows);
            assert_eq!(
                selected.criterion_word_xors,
                expected.report.criterion_word_ops
            );
            assert_eq!(integer(call, "selected_rows"), selected.rows.len() as u64);
            assert_eq!(integer(call, "rows_f4"), selected.rows_f4);
            assert_eq!(integer(call, "criterion_rows"), selected.criterion_rows);
            assert_eq!(
                integer(call, "criterion_word_xors"),
                selected.criterion_word_xors
            );
            if call["cache_hit"] == true {
                let expected_reason = if arm == "anchor" {
                    assert!(!projected);
                    "hit"
                } else if projected {
                    "hit-projected"
                } else {
                    "hit-full"
                };
                assert_eq!(call["fallback_reason"], expected_reason);
                let total = integer(call, "candidate_word_xors");
                let transform = integer(call, "candidate_transform_word_xors");
                let tail = integer(call, "candidate_tail_word_xors");
                assert_eq!(total, selected.criterion_word_xors + transform + tail);
                word_xors += u128::from(total);
                hits += 1;
            } else {
                assert!(!projected);
                assert_eq!(integer(call, "candidate_word_xors"), 0);
                assert_eq!(integer(call, "candidate_transform_word_xors"), 0);
                assert_eq!(integer(call, "candidate_tail_word_xors"), 0);
                word_xors += u128::from(expected.report.word_ops());
                fallbacks += 1;
            }
        } else {
            assert_eq!(
                call["report"],
                serde_json::to_value(expected.report).unwrap()
            );
            assert_eq!(call["cache_hit"], false);
            assert!(!projected);
            assert_eq!(call["fallback_reason"], "baseline");
            assert_eq!(integer(call, "candidate_word_xors"), 0);
            for field in [
                "candidate_transform_word_xors",
                "candidate_tail_word_xors",
                "selected_rows",
                "rows_f4",
                "criterion_rows",
                "criterion_word_xors",
            ] {
                assert_eq!(integer(call, field), 0);
            }
            word_xors += u128::from(expected.report.word_ops());
        }
    }
    assert_eq!(integer(line, "total_ns") as u128, total);
    assert_eq!(integer(line, "validation_ns") as u128, validation);
    (
        total as f64,
        hits,
        fallbacks,
        word_xors,
        eligible,
        recovered,
    )
}

fn compute_report(phase: &str, raw: &Path, conditions: &Path) -> Value {
    let p = route();
    check_host();
    check_environment(&p.env);
    let seeds = match phase {
        "discovery" => &p.discovery_seeds,
        "holdout" => &p.holdout_seeds,
        _ => panic!("phase"),
    };
    let size = fs::metadata(raw).expect("raw size").len() as usize;
    assert!(size <= p.evidence_cap_bytes);
    let lines = raw_lines(raw);
    let mut cursor = 0;
    let mut take = || {
        let row = lines
            .get(cursor)
            .unwrap_or_else(|| panic!("truncated at {cursor}"));
        cursor += 1;
        row
    };
    let header = json!({
        "kind":"campaign","phase":phase,"protocol_sha256":sha256(PROTOCOL_BYTES),
        "worker_sha256":sha256(WORKER_BYTES),"campaign_sha256":sha256(CAMPAIGN_BYTES),
        "gf2_sha256":sha256(GF2_BYTES),"host_arch":"x86_64","route":p.env
    });
    assert_eq!(take(), &header);
    let mut groups = Vec::new();
    let mut primary_passed = 0;
    let mut nonregression_passed = 0;
    let mut cells = 0;
    let mut timed_calls = 0;
    let mut cache_hits = 0;
    let mut fallbacks = 0;
    for &n in &p.variables {
        for &seed in seeds {
            for family in &p.families {
                for &batch in &p.batches {
                    let systems = assignments(n, seed, batch, family);
                    let polys = systems.iter().map(polynomials).collect::<Vec<_>>();
                    let references = polys.iter().map(|p| reference(p, n)).collect::<Vec<_>>();
                    let eligible_sparse = sparse_eligibility(&systems);
                    let fingerprint = f4_cross_check(&polys[0], n);
                    let cell_id = format!("n{n}-{phase}-{seed}-{family}-b{batch}");
                    let fixture = json!({
                        "kind":"fixture","cell":cell_id,"n":n,"seed":seed,
                        "family":family,"batch":batch,"quadratic":systems[0].quadratic,
                        "affine":systems.iter().map(|s| &s.affine).collect::<Vec<_>>(),
                        "references":references,"small_f4_fingerprint":fingerprint,
                        "eligible_sparse":eligible_sparse,
                        "direct_preflight":true
                    });
                    assert_eq!(take(), &fixture, "fixture {cell_id}");
                    // Independent direct replay of each unique output. The
                    // producer preflight compares vectors, and every timed
                    // return carries an ordered-output digest.
                    let mut replay_caches = [
                        "anchor",
                        "support",
                        "support_unpack",
                        "support_unpack_layout",
                    ]
                    .into_iter()
                    .map(|arm| (arm, CacheContext::compile(arm, &systems[0]).ok()))
                    .collect::<Vec<_>>();
                    for (index, system) in systems.iter().enumerate() {
                        let (expected, report, _) = call_f5(&polys[index], n);
                        assert_eq!(report, references[index].report);
                        for (arm, cache) in &mut replay_caches {
                            let actual = if let Some(cache) = cache {
                                cache.apply(arm, system).0
                            } else {
                                baseline(system)
                            };
                            assert_eq!(
                                actual, expected,
                                "replayed {arm} ordered output {cell_id}/{index}"
                            );
                        }
                    }
                    let mut ratios = Vec::with_capacity(p.repetitions);
                    let mut support_anchor_ratios = Vec::with_capacity(p.repetitions);
                    let mut arm_ratios: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
                    let mut arm_counts: BTreeMap<&str, (usize, usize, usize, usize)> =
                        BTreeMap::new();
                    let mut counted_ratios = Vec::with_capacity(p.repetitions);
                    let mut aa_ratios = Vec::with_capacity(p.repetitions);
                    for repetition in 0..p.repetitions {
                        let mut times = BTreeMap::new();
                        let mut counted = BTreeMap::new();
                        for arm in arm_order(repetition) {
                            let sample = take();
                            let (time, hits, missed, words, eligible, recovered) = verify_sample(
                                sample,
                                &cell_id,
                                repetition,
                                arm,
                                &references,
                                &systems,
                                &eligible_sparse,
                                &p,
                            );
                            times.insert(arm, time);
                            counted.insert(arm, words);
                            let counts = arm_counts.entry(arm).or_default();
                            counts.0 += hits;
                            counts.1 += missed;
                            counts.2 += eligible;
                            counts.3 += recovered;
                            cache_hits += hits;
                            fallbacks += missed;
                            timed_calls += batch;
                        }
                        let aa = (times["aa_a"] / times["aa_b"]).max(times["aa_b"] / times["aa_a"]);
                        aa_ratios.push(aa);
                        ratios.push(times["fresh"] / times["support_unpack_layout"]);
                        support_anchor_ratios.push(times["anchor"] / times["support"]);
                        for arm in [
                            "anchor",
                            "support",
                            "support_unpack",
                            "support_unpack_layout",
                        ] {
                            arm_ratios
                                .entry(arm)
                                .or_default()
                                .push(times["fresh"] / times[arm]);
                        }
                        counted_ratios.push(
                            counted["fresh"] as f64 / counted["support_unpack_layout"] as f64,
                        );
                    }
                    let family_code = if family == "independent_affine" { 1 } else { 2 };
                    let lower = bootstrap_lower(
                        &ratios,
                        p.bootstrap_seed ^ (u64::from(n) << 32) ^ seed ^ family_code ^ batch as u64,
                        p.bootstrap_resamples,
                    );
                    let centre = median(&ratios);
                    let aa_floor = quantile(&aa_ratios, p.aa_noise_quantile);
                    let attribution_lower = bootstrap_lower(
                        &support_anchor_ratios,
                        p.bootstrap_seed
                            ^ (u64::from(n) << 32)
                            ^ seed
                            ^ family_code
                            ^ batch as u64
                            ^ 0x5355_5050_4f52_54,
                        p.bootstrap_resamples,
                    );
                    let sparse_expected =
                        eligible_sparse.iter().filter(|&&yes| yes).count() * p.repetitions;
                    let structural_pass = sparse_expected > 0
                        && ["support", "support_unpack", "support_unpack_layout"]
                            .iter()
                            .all(|arm| {
                                let counts = arm_counts[arm];
                                counts.2 == sparse_expected && counts.3 == sparse_expected
                            });
                    let primary = p.primary_variables.contains(&n) && batch == p.primary_batch;
                    let nonregression = p.nonregression_variables.contains(&n) && batch == 32;
                    let pass = if primary {
                        centre > p.primary_median_greater_than
                            && lower > p.primary_lower_bound_greater_than
                            && centre > aa_floor
                            && lower > aa_floor
                            && attribution_lower > aa_floor
                            && structural_pass
                    } else if nonregression {
                        lower > p.nonregression_lower_bound_greater_than
                    } else {
                        true
                    };
                    if primary {
                        primary_passed += usize::from(pass);
                    }
                    if nonregression {
                        nonregression_passed += usize::from(pass);
                    }
                    let mut ablations = serde_json::Map::new();
                    for (arm_index, arm) in [
                        "anchor",
                        "support",
                        "support_unpack",
                        "support_unpack_layout",
                    ]
                    .iter()
                    .enumerate()
                    {
                        let values = &arm_ratios[arm];
                        let bounds = bootstrap_lower(
                            values,
                            p.bootstrap_seed
                                ^ (u64::from(n) << 32)
                                ^ seed
                                ^ family_code
                                ^ batch as u64
                                ^ (arm_index as u64 + 1),
                            p.bootstrap_resamples,
                        );
                        let counts = arm_counts[arm];
                        ablations.insert(
                            (*arm).to_owned(),
                            json!({
                                "fresh_over_arm_median":median(values),
                                "bootstrap_lower_95":bounds,
                                "hits":counts.0,"fallbacks":counts.1,
                                "eligible_sparse":counts.2,"recovered_sparse":counts.3
                            }),
                        );
                    }
                    groups.push(json!({"n":n,"seed":seed,"family":family,"batch":batch,
                        "paired_median":centre,"bootstrap_lower_95":lower,
                        "counted_word_xor_ratio_median":median(&counted_ratios),
                        "aa_noise_floor_975":aa_floor,"primary":primary,
                        "support_over_anchor_median":median(&support_anchor_ratios),
                        "support_over_anchor_lower_95":attribution_lower,
                        "structural_pass":structural_pass,"sparse_expected_per_arm":sparse_expected,
                        "ablations":ablations,
                        "nonregression":nonregression,"pass":pass}));
                    cells += 1;
                }
            }
        }
    }
    assert_eq!(cursor, lines.len(), "trailing raw records");
    assert_eq!(cells, 48);
    assert_eq!(
        groups.iter().filter(|g| g["primary"] == true).count(),
        p.primary_groups
    );
    assert_eq!(
        groups.iter().filter(|g| g["nonregression"] == true).count(),
        8
    );
    let resource = resource(conditions, phase, &p);
    let readiness = readiness(&raw.with_file_name("readiness.json"), &p);
    let gate_pass = primary_passed == p.primary_groups && nonregression_passed == 8;
    json!({
        "schema_version":1,"phase":phase,"protocol_sha256":sha256(PROTOCOL_BYTES),
        "worker_sha256":sha256(WORKER_BYTES),"campaign_sha256":sha256(CAMPAIGN_BYTES),
        "gf2_sha256":sha256(GF2_BYTES),"raw_sha256":sha256(&fs::read(raw).unwrap()),
        "raw_bytes":size,"cells":cells,"timed_calls":timed_calls,
        "cache_hits":cache_hits,"fallbacks":fallbacks,"groups":groups,
        "primary_passed":primary_passed,"nonregression_passed":nonregression_passed,
        "gate_pass":gate_pass,"resource":resource,"readiness":readiness,
        "full_ic_cost":null,"rho_ratio":null
    })
}

fn verify_run(phase: &str, raw: &str, results: &str) {
    let output = Path::new(results);
    assert!(!output.exists(), "refuse to overwrite results");
    if phase == "holdout" {
        validate_binding(&output.with_file_name("binding.json"));
    }
    let report = compute_report(
        phase,
        Path::new(raw),
        &output.with_file_name("conditions.jsonl"),
    );
    fs::write(output, serde_json::to_vec_pretty(&report).unwrap()).expect("results");
}

fn validate_binding(path: &Path) {
    let binding = read_json(path);
    assert_eq!(binding["protocol_sha256"], sha256(PROTOCOL_BYTES));
    assert_eq!(binding["worker_sha256"], sha256(WORKER_BYTES));
    assert_eq!(binding["campaign_sha256"], sha256(CAMPAIGN_BYTES));
    assert_eq!(binding["gf2_sha256"], sha256(GF2_BYTES));
    for field in ["discovery_manifest_sha256", "discovery_results_sha256"] {
        assert_eq!(binding[field].as_str().unwrap().len(), 64);
    }
}

fn seal(dir: &str) {
    let root = Path::new(dir);
    let manifest = root.join("manifest.json");
    assert!(!manifest.exists(), "refuse to overwrite manifest");
    let mut files = BTreeMap::new();
    for entry in fs::read_dir(root).expect("bundle directory") {
        let entry = entry.expect("entry");
        assert!(!entry.file_type().unwrap().is_symlink());
        if !entry.file_type().unwrap().is_file() {
            continue;
        }
        let name = entry.file_name().into_string().expect("UTF-8 name");
        assert_ne!(name, "manifest.json");
        files.insert(
            name,
            sha256(&fs::read(entry.path()).expect("bundle member")),
        );
    }
    for required in ["worker.rs", "campaign.rs", "protocol.json", "raw.jsonl"] {
        assert!(files.contains_key(required), "missing {required}");
    }
    fs::write(
        manifest,
        serde_json::to_vec_pretty(&json!({"schema_version":1,"files":files})).unwrap(),
    )
    .expect("manifest");
}

fn verify_bundle(dir: &str) -> Value {
    let root = Path::new(dir);
    let manifest = read_json(&root.join("manifest.json"));
    assert_eq!(manifest["schema_version"], 1);
    let files = manifest["files"].as_object().expect("files");
    assert!(files.len() >= 8);
    for (name, digest) in files {
        assert!(!name.contains('/') && !name.contains(".."));
        assert_eq!(
            sha256(&fs::read(root.join(name)).expect("member")),
            digest.as_str().unwrap()
        );
    }
    assert_eq!(fs::read(root.join("worker.rs")).unwrap(), WORKER_BYTES);
    assert_eq!(fs::read(root.join("campaign.rs")).unwrap(), CAMPAIGN_BYTES);
    assert_eq!(
        fs::read(root.join("protocol.json")).unwrap(),
        PROTOCOL_BYTES
    );
    assert_eq!(fs::read(root.join("gf2_elim.rs")).unwrap(), GF2_BYTES);
    if root.join("failure.json").exists() {
        let failure = read_json(&root.join("failure.json"));
        assert_eq!(failure["performance_admitted"], false);
        return failure;
    }
    let old = read_json(&root.join("results.json"));
    let phase = old["phase"].as_str().expect("phase");
    if phase == "holdout" {
        validate_binding(&root.join("binding.json"));
    }
    let expected = compute_report(
        phase,
        &root.join("raw.jsonl"),
        &root.join("conditions.jsonl"),
    );
    assert_eq!(old, expected, "replay differs");
    old
}

fn check_discovery(bundle: &str, binding_path: &str) {
    let result = verify_bundle(bundle);
    assert_eq!(result["phase"], "discovery");
    assert_eq!(result["gate_pass"], true, "discovery gate failed");
    let target = Path::new(binding_path);
    assert!(!target.exists());
    let root = Path::new(bundle);
    fs::write(
        target,
        serde_json::to_vec_pretty(&json!({
            "protocol_sha256":sha256(PROTOCOL_BYTES),"worker_sha256":sha256(WORKER_BYTES),
            "campaign_sha256":sha256(CAMPAIGN_BYTES),"gf2_sha256":sha256(GF2_BYTES),
            "discovery_manifest_sha256":sha256(&fs::read(root.join("manifest.json")).unwrap()),
            "discovery_results_sha256":sha256(&fs::read(root.join("results.json")).unwrap())
        }))
        .unwrap(),
    )
    .expect("binding");
}

fn failure(dir: &str, reason: &str) {
    let path = Path::new(dir).join("failure.json");
    assert!(!path.exists());
    fs::write(
        path,
        serde_json::to_vec_pretty(&json!({
            "reason":reason,"performance_admitted":false,"full_ic_cost":null,"rho_ratio":null
        }))
        .unwrap(),
    )
    .expect("failure");
}

fn psi_some_avg10() -> f64 {
    let content = fs::read_to_string("/proc/pressure/cpu").expect("CPU PSI");
    content
        .lines()
        .find(|line| line.starts_with("some "))
        .expect("some PSI")
        .split_whitespace()
        .find_map(|field| field.strip_prefix("avg10="))
        .expect("avg10")
        .parse()
        .expect("numeric PSI")
}

fn wait_quiet(path: &str) {
    let target = Path::new(path);
    assert!(!target.exists());
    let p = route();
    let mut samples = Vec::new();
    for index in 0..30 {
        let psi = psi_some_avg10();
        samples.push(json!({"index":index,"psi_some_avg10":psi,
            "accepted":psi<=p.maximum_psi_some_avg10}));
        if psi <= p.maximum_psi_some_avg10 {
            fs::write(
                target,
                serde_json::to_vec_pretty(&json!({
                    "accepted":true,"samples":samples
                }))
                .unwrap(),
            )
            .expect("readiness");
            return;
        }
        thread::sleep(Duration::from_secs(2));
    }
    fs::write(
        target,
        serde_json::to_vec_pretty(&json!({
            "accepted":false,"samples":samples
        }))
        .unwrap(),
    )
    .expect("readiness refusal");
    panic!("CPU PSI remained above bound");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_fixture_recurrence_matches_prior_native_study() {
        let systems = assignments(12, 17, 8, "independent_affine");
        assert_eq!(systems[0], fixture(12, 17));
        assert_eq!(systems[0].affine, vec![0; 12]);
        assert_eq!(systems.len(), 8);
    }

    #[test]
    fn bootstrap_lower_is_deterministic() {
        let ratios = [1.2, 1.5, 1.8, 2.1, 2.4, 2.6, 3.0];
        assert_eq!(
            bootstrap_lower(&ratios, 17, 1000),
            bootstrap_lower(&ratios, 17, 1000)
        );
    }

    #[test]
    fn cpu_sibling_ranges_are_parsed_exactly() {
        assert_eq!(parse_cpu_list("2-3\n"), vec![2, 3]);
        assert_eq!(parse_cpu_list("0,2-3"), vec![0, 2, 3]);
    }

    #[test]
    fn development_receipts_replay_against_native_reference() {
        let systems = assignments(12, 17, 2, "independent_affine");
        let polys = systems.iter().map(polynomials).collect::<Vec<_>>();
        let references = polys
            .iter()
            .map(|input| reference(input, 12))
            .collect::<Vec<_>>();
        let eligible_sparse = sparse_eligibility(&systems);
        let p = route();
        for arm in [
            "fresh",
            "anchor",
            "support",
            "support_unpack",
            "support_unpack_layout",
        ] {
            let mut sample = timed_batch(&systems, &polys, &references, 0, arm);
            sample["cell"] = json!("development-n12");
            let (time, hits, misses, words, eligible, recovered) = verify_sample(
                &sample,
                "development-n12",
                0,
                arm,
                &references,
                &systems,
                &eligible_sparse,
                &p,
            );
            assert!(time > 0.0 && words > 0);
            if arm != "fresh" {
                assert_eq!(hits + misses, 2);
            } else {
                assert_eq!(hits + misses, 0);
            }
            assert!(recovered <= eligible);
        }
    }

    #[test]
    fn resource_gate_rejects_contended_receipt() {
        let root = std::env::temp_dir().join(format!("f5-support-resource-{}", std::process::id()));
        assert!(!root.exists());
        fs::create_dir(&root).unwrap();
        let receipt_path = root.join("conditions.jsonl");
        fs::write(root.join("cpu-siblings.txt"), b"2-3\n").unwrap();
        let mut receipt = json!({
            "schema":"isolated-bench/1","mode":"run","label":"f5-support-aware/discovery",
            "host":{"machine":"x86_64"},"reserved_cpus":[2,3],
            "preflight":{"settle":{"seconds":2.0,"other_cpu_seconds":0.02},
                "conditions":{"psi_cpu":{"some":{"avg10":1.0}}}},
            "run":{"contended":false,"exit_status":0,"wall_seconds":10.0,
                "other_cpu_seconds":0.2,"max_rss_kib":8192}
        });
        fs::write(&receipt_path, format!("{receipt}\n")).unwrap();
        resource(&receipt_path, "discovery", &route());
        receipt["run"]["contended"] = json!(true);
        fs::write(&receipt_path, format!("{receipt}\n")).unwrap();
        assert!(
            std::panic::catch_unwind(|| resource(&receipt_path, "discovery", &route())).is_err()
        );
        fs::remove_dir_all(root).unwrap();
    }
}
