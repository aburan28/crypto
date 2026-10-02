//! Native phase accounting, resource admission and sealed replay.

use super::*;
use serde_json::Value;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::thread;
use std::time::Duration;

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("read JSON")).expect("parse JSON")
}

fn raw_lines(path: &Path) -> Vec<Value> {
    let input = fs::File::open(path).expect("raw campaign");
    BufReader::new(input)
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
    assert!(!values.is_empty());
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let index = ((sorted.len() as f64 * fraction).ceil() as usize)
        .saturating_sub(1)
        .min(sorted.len() - 1);
    sorted[index]
}

fn bootstrap_upper(values: &[f64], seed: u64, resamples: usize) -> f64 {
    let mut state = seed;
    let mut medians = Vec::with_capacity(resamples);
    for _ in 0..resamples {
        let mut sample = Vec::with_capacity(values.len());
        for _ in values {
            sample.push(values[(next(&mut state) % values.len() as u64) as usize]);
        }
        medians.push(median(&sample));
    }
    quantile(&medians, 0.975)
}

fn resource(path: &Path, phase: &str, worker_seconds: u64) -> Value {
    let lines = raw_lines(path);
    assert_eq!(lines.len(), 1);
    let receipt = &lines[0];
    assert_eq!(receipt["schema"], "isolated-bench/1");
    assert_eq!(receipt["mode"], "run");
    assert_eq!(receipt["label"], format!("f5-ceiling/{phase}"));
    assert_eq!(receipt["reserved_cpus"].as_array().unwrap().len(), 1);
    assert_eq!(receipt["run"]["contended"], false);
    assert_eq!(receipt["run"]["exit_status"], 0);
    assert_eq!(receipt["host"]["machine"], "x86_64");
    let settle_seconds = receipt["preflight"]["settle"]["seconds"].as_f64().unwrap();
    let settle_other = receipt["preflight"]["settle"]["other_cpu_seconds"]
        .as_f64()
        .unwrap();
    assert!(settle_seconds > 0.0 && settle_other / settle_seconds <= 0.10);
    let psi = receipt["preflight"]["conditions"]["psi_cpu"]["some"]["avg10"]
        .as_f64()
        .unwrap();
    assert!(psi <= 5.0);
    let wall = receipt["run"]["wall_seconds"].as_f64().unwrap();
    let other = receipt["run"]["other_cpu_seconds"].as_f64().unwrap();
    assert!(wall > 0.0 && wall <= worker_seconds as f64 && other / wall <= 0.10);
    json!({
        "wall_seconds":wall,
        "other_cpu_seconds":other,
        "psi_some_avg10":psi,
        "max_rss_kib":integer(&receipt["run"],"max_rss_kib"),
        "contended":false,
        "reserved_cpus":receipt["reserved_cpus"]
    })
}

fn readiness(path: &Path) -> Value {
    let record = read_json(path);
    assert_eq!(record["accepted"], true);
    let samples = record["samples"].as_array().expect("samples");
    assert!(!samples.is_empty() && samples.len() <= 30);
    for (index, sample) in samples.iter().enumerate() {
        assert_eq!(integer(sample, "index"), index as u64);
        let psi = sample["psi_some_avg10"].as_f64().unwrap();
        assert_eq!(sample["accepted"], psi <= 5.0);
    }
    assert_eq!(samples.last().unwrap()["accepted"], true);
    record
}

fn verify_sample(
    line: &Value,
    cell_id: &str,
    repetition: usize,
    arm: usize,
    references: &[Reference],
    batch: usize,
) -> f64 {
    assert_eq!(line["kind"], "sample");
    assert_eq!(line["cell"], cell_id);
    assert_eq!(integer(line, "repetition"), repetition as u64);
    assert_eq!(line["arm"], if arm == 0 { "aa_a" } else { "aa_b" });
    let calls = line["calls"].as_array().expect("call records");
    assert_eq!(calls.len(), batch);
    let mut total = 0u128;
    let mut criterion = 0u128;
    let mut build = 0u128;
    let mut reduce = 0u128;
    let mut unpack = 0u128;
    let mut validation = 0u128;
    for (position, call) in calls.iter().enumerate() {
        let index = (position + repetition + arm) % batch;
        assert_eq!(integer(call, "input_index"), index as u64);
        assert_eq!(
            call["report"],
            serde_json::to_value(references[index].report).unwrap()
        );
        assert_eq!(call["output_digest"], references[index].output_digest);
        assert_eq!(
            integer(call, "returned_rows"),
            references[index].returned_rows as u64
        );
        assert_eq!(
            integer(call, "returned_terms"),
            references[index].returned_terms as u64
        );
        assert_eq!(call["direct_pack_used"], references[index].direct_pack_used);
        assert_eq!(call["direct_unpack_used"], true);
        let t = integer(call, "total_ns") as u128;
        let c = integer(call, "criterion_ns") as u128;
        let b = integer(call, "build_ns") as u128;
        let r = integer(call, "reduce_ns") as u128;
        let u = integer(call, "unpack_ns") as u128;
        let v = integer(call, "validation_ns") as u128;
        assert!(t > b + r && c + b + r + u + v <= t);
        total += t;
        criterion += c;
        build += b;
        reduce += r;
        unpack += u;
        validation += v;
    }
    assert_eq!(integer(line, "total_ns") as u128, total);
    assert_eq!(integer(line, "criterion_ns") as u128, criterion);
    assert_eq!(integer(line, "build_ns") as u128, build);
    assert_eq!(integer(line, "reduce_ns") as u128, reduce);
    assert_eq!(integer(line, "unpack_ns") as u128, unpack);
    assert_eq!(integer(line, "validation_ns") as u128, validation);
    let ceiling = total as f64 / (total - build - reduce) as f64;
    assert_eq!(
        line["optimistic_ceiling"].as_f64().unwrap().to_bits(),
        ceiling.to_bits()
    );
    ceiling
}

fn compute_report(phase: &str, raw: &Path, conditions: &Path) -> Value {
    let (route_env, protocol) = route();
    check_host();
    check_route_environment(&route_env);
    let (split, seeds) = match phase {
        "discovery" => ("discovery", protocol.discovery_seeds),
        "holdout" => ("holdout", protocol.holdout_seeds),
        _ => panic!("phase"),
    };
    let size = fs::metadata(raw).expect("raw size").len() as usize;
    assert!(size <= protocol.evidence_cap_bytes);
    let lines = raw_lines(raw);
    let mut cursor = 0usize;
    let mut take = || {
        let row = lines
            .get(cursor)
            .unwrap_or_else(|| panic!("truncated at {cursor}"));
        cursor += 1;
        row
    };
    let expected_header = json!({
        "kind":"campaign","phase":phase,
        "protocol_sha256":sha256(PROTOCOL_BYTES),
        "source_sha256":sha256(SOURCE_BYTES),
        "verifier_sha256":sha256(VERIFY_BYTES),
        "host_arch":"x86_64","route":route_env
    });
    assert_eq!(take(), &expected_header);
    let mut cells = 0;
    let mut calls = 0;
    let mut direct_pack_hits = 0;
    let mut primary = Vec::new();
    let mut gate = true;
    for &n in &protocol.variables {
        for &seed in &seeds {
            for family in &protocol.families {
                for &batch in &protocol.batches {
                    let systems = assignments(n, seed, batch, family);
                    let polys = systems.iter().map(polynomials).collect::<Vec<_>>();
                    let references = polys
                        .iter()
                        .map(|input| reference(input, n))
                        .collect::<Vec<_>>();
                    let f4_fingerprint = f4_cross_check(&polys[0], n);
                    let cell_id = format!("n{n}-{split}-{seed}-{family}-b{batch}");
                    let expected_fixture = json!({
                        "kind":"fixture","cell":cell_id,"n":n,"seed":seed,
                        "family":family,"batch":batch,
                        "quadratic":systems[0].quadratic,
                        "affine":systems.iter().map(|s|&s.affine).collect::<Vec<_>>(),
                        "references":references,
                        "small_f4_fingerprint":f4_fingerprint
                    });
                    assert_eq!(take(), &expected_fixture, "fixture {cell_id}");
                    let mut ceilings = Vec::new();
                    let mut aa_ratios = Vec::new();
                    for repetition in 0..protocol.repetitions {
                        let a = take();
                        let b = take();
                        let ua = verify_sample(a, &cell_id, repetition, 0, &references, batch);
                        let ub = verify_sample(b, &cell_id, repetition, 1, &references, batch);
                        direct_pack_hits += a["calls"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .chain(b["calls"].as_array().unwrap())
                            .filter(|call| call["direct_pack_used"] == true)
                            .count();
                        let ta = integer(a, "total_ns") as f64;
                        let tb = integer(b, "total_ns") as f64;
                        aa_ratios.push((ta / tb).max(tb / ta));
                        ceilings.push(ua);
                        ceilings.push(ub);
                        calls += batch * 2;
                    }
                    if protocol.primary_variables.contains(&n) && batch == protocol.primary_batch {
                        let family_code = if family == "independent_affine" { 1 } else { 2 };
                        let upper = bootstrap_upper(
                            &ceilings,
                            protocol.bootstrap_seed ^ (u64::from(n) << 32) ^ seed ^ family_code,
                            protocol.bootstrap_resamples,
                        );
                        let aa_floor = quantile(&aa_ratios, 0.975);
                        let pass = upper > protocol.upper_bound_greater_than;
                        gate &= pass;
                        primary.push(json!({
                            "n":n,"seed":seed,"family":family,"batch":batch,
                            "observations":ceilings.len(),
                            "median_optimistic_ceiling":median(&ceilings),
                            "bootstrap_upper_95":upper,
                            "aa_noise_floor_975":aa_floor,
                            "pass":pass
                        }));
                    }
                    cells += 1;
                }
            }
        }
    }
    assert_eq!(cursor, lines.len(), "trailing records");
    assert_eq!(cells, 48);
    assert_eq!(primary.len(), 4);
    let resource = resource(conditions, phase, protocol.worker_seconds);
    let readiness = readiness(&raw.with_file_name("readiness.json"));
    json!({
        "schema_version":1,
        "phase":phase,
        "protocol_sha256":sha256(PROTOCOL_BYTES),
        "source_sha256":sha256(SOURCE_BYTES),
        "verifier_sha256":sha256(VERIFY_BYTES),
        "raw_sha256":sha256(&fs::read(raw).unwrap()),
        "raw_bytes":size,
        "cells":cells,"timed_calls":calls,
        "direct_pack_hits":direct_pack_hits,
        "direct_pack_fallbacks":calls-direct_pack_hits,
        "resource":resource,
        "readiness":readiness,
        "primary_groups":primary,
        "primary_passed":primary.iter().filter(|row|row["pass"]==true).count(),
        "gate_pass":gate,
        "complete_f5_speedup":null,
        "full_ic_cost":null,
        "rho_ratio":null
    })
}

pub(super) fn verify_run(phase: &str, raw: &str, results: &str) {
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
    fs::write(output, serde_json::to_vec_pretty(&report).unwrap()).expect("write result");
}

fn validate_binding(path: &Path) {
    let binding = read_json(path);
    assert_eq!(binding["protocol_sha256"], sha256(PROTOCOL_BYTES));
    assert_eq!(binding["source_sha256"], sha256(SOURCE_BYTES));
    assert_eq!(binding["verifier_sha256"], sha256(VERIFY_BYTES));
    for field in ["discovery_manifest_sha256", "discovery_results_sha256"] {
        assert_eq!(binding[field].as_str().unwrap().len(), 64);
    }
}

pub(super) fn seal(dir: &str) {
    let root = Path::new(dir);
    let manifest = root.join("manifest.json");
    assert!(!manifest.exists(), "refuse to overwrite manifest");
    let mut files = BTreeMap::new();
    for entry in fs::read_dir(root).expect("bundle") {
        let entry = entry.expect("entry");
        let kind = entry.file_type().expect("file type");
        assert!(!kind.is_symlink());
        if !kind.is_file() {
            continue;
        }
        let name = entry.file_name().into_string().expect("UTF-8 member");
        assert_ne!(name, "manifest.json");
        files.insert(name, sha256(&fs::read(entry.path()).expect("member")));
    }
    for required in ["worker.rs", "verify.rs", "protocol.json", "raw.jsonl"] {
        assert!(files.contains_key(required), "missing {required}");
    }
    fs::write(
        manifest,
        serde_json::to_vec_pretty(&json!({
            "schema_version":1,"files":files
        }))
        .unwrap(),
    )
    .expect("manifest");
}

pub(super) fn verify_bundle(dir: &str) -> Value {
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
    assert_eq!(fs::read(root.join("worker.rs")).unwrap(), SOURCE_BYTES);
    assert_eq!(fs::read(root.join("verify.rs")).unwrap(), VERIFY_BYTES);
    assert_eq!(
        fs::read(root.join("protocol.json")).unwrap(),
        PROTOCOL_BYTES
    );
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
    assert_eq!(old, expected, "bundle replay differs");
    old
}

pub(super) fn check_discovery(bundle: &str, binding_path: &str) {
    let result = verify_bundle(bundle);
    assert_eq!(result["phase"], "discovery");
    assert_eq!(result["gate_pass"], true, "discovery ceiling gate failed");
    assert!(!Path::new(binding_path).exists());
    let root = Path::new(bundle);
    fs::write(
        binding_path,
        serde_json::to_vec_pretty(&json!({
            "protocol_sha256":sha256(PROTOCOL_BYTES),
            "source_sha256":sha256(SOURCE_BYTES),
            "verifier_sha256":sha256(VERIFY_BYTES),
            "discovery_manifest_sha256":sha256(&fs::read(root.join("manifest.json")).unwrap()),
            "discovery_results_sha256":sha256(&fs::read(root.join("results.json")).unwrap())
        }))
        .unwrap(),
    )
    .expect("binding");
}

pub(super) fn failure(dir: &str, reason: &str) {
    let path = Path::new(dir).join("failure.json");
    assert!(!path.exists());
    fs::write(
        path,
        serde_json::to_vec_pretty(&json!({
            "reason":reason,"performance_admitted":false,
            "complete_f5_speedup":null,"full_ic_cost":null,"rho_ratio":null
        }))
        .unwrap(),
    )
    .expect("failure");
}

fn psi_some_avg10() -> f64 {
    let contents = fs::read_to_string("/proc/pressure/cpu").expect("CPU PSI");
    let line = contents
        .lines()
        .find(|line| line.starts_with("some "))
        .expect("some PSI");
    line.split_whitespace()
        .find_map(|field| field.strip_prefix("avg10="))
        .expect("avg10")
        .parse()
        .expect("numeric PSI")
}

pub(super) fn wait_quiet(path: &str) {
    let target = Path::new(path);
    assert!(!target.exists(), "refuse to overwrite readiness");
    let mut samples = Vec::new();
    for index in 0..30 {
        let psi = psi_some_avg10();
        samples.push(json!({"index":index,"psi_some_avg10":psi,"accepted":psi<=5.0}));
        if psi <= 5.0 {
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
    panic!("CPU PSI remained above 5.0");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ceiling_bootstrap_is_deterministic() {
        let values = [1.2, 1.5, 1.7, 2.1, 2.4, 2.8, 3.0];
        assert_eq!(
            bootstrap_upper(&values, 17, 1000),
            bootstrap_upper(&values, 17, 1000)
        );
    }

    #[test]
    fn reported_ceiling_survives_json_roundtrip() {
        for value in [1.9999999999999998_f64, 2.0000000000000004_f64] {
            let text = serde_json::to_string(&value).unwrap();
            let parsed: f64 = serde_json::from_str(&text).unwrap();
            assert_eq!(value.to_bits(), parsed.to_bits());
        }
    }

    #[test]
    fn resource_gate_rejects_contended_receipt() {
        let path = std::env::temp_dir().join(format!("f5-ceiling-resource-{}", std::process::id()));
        assert!(!path.exists());
        let mut record = json!({
            "schema":"isolated-bench/1","mode":"run","label":"f5-ceiling/discovery",
            "host":{"machine":"x86_64"},"reserved_cpus":[1],
            "preflight":{"settle":{"seconds":2.0,"other_cpu_seconds":0.02},
                "conditions":{"psi_cpu":{"some":{"avg10":1.0}}}},
            "run":{"contended":false,"exit_status":0,"wall_seconds":10.0,
                "other_cpu_seconds":0.2,"max_rss_kib":8192}
        });
        fs::write(&path, format!("{record}\n")).unwrap();
        resource(&path, "discovery", 900);
        record["run"]["contended"] = json!(true);
        fs::write(&path, format!("{record}\n")).unwrap();
        assert!(std::panic::catch_unwind(|| resource(&path, "discovery", 900)).is_err());
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn sealed_failure_cannot_be_promoted() {
        let root = std::env::temp_dir().join(format!("f5-ceiling-failure-{}", std::process::id()));
        assert!(!root.exists());
        fs::create_dir(&root).unwrap();
        fs::write(root.join("worker.rs"), SOURCE_BYTES).unwrap();
        fs::write(root.join("verify.rs"), VERIFY_BYTES).unwrap();
        fs::write(root.join("protocol.json"), PROTOCOL_BYTES).unwrap();
        fs::write(root.join("raw.jsonl"), b"").unwrap();
        for name in ["Cargo.toml", "Cargo.lock", "PROTOCOL.md", "README.md"] {
            fs::write(root.join(name), b"retained").unwrap();
        }
        failure(root.to_str().unwrap(), "resource-refused");
        seal(root.to_str().unwrap());
        assert_eq!(
            verify_bundle(root.to_str().unwrap())["performance_admitted"],
            false
        );
        fs::remove_dir_all(root).unwrap();
    }
}
