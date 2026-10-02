//! Native replay, resource admission, and sealed evidence for graded reuse.

use super::*;
use serde_json::{json, Value};
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::thread;
use std::time::Duration;

fn hash(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("read JSON")).expect("parse JSON")
}

fn number(value: &Value, field: &str) -> u64 {
    value[field]
        .as_u64()
        .unwrap_or_else(|| panic!("missing numeric {field}"))
}

fn value_lines(path: &Path) -> Vec<Value> {
    let input = fs::File::open(path).expect("open raw campaign");
    BufReader::new(input)
        .lines()
        .map(|line| serde_json::from_str(&line.expect("raw line")).expect("raw JSON"))
        .collect()
}

fn expected_sample(
    arm: Arm,
    inputs: &[System],
    expected: &[(Vec<Row>, usize)],
    cap: usize,
) -> Option<Value> {
    let mut engine = Engine::compile(arm, &inputs[0])?;
    let retained_bytes = engine.retained_bytes();
    assert!(retained_bytes <= cap);
    let setup_work = engine.setup_work();
    let mut work = Work::default();
    let mut hasher = Sha256::new();
    hasher.update(b"graded-macaulay-output-v1");
    let mut hits = 0;
    let mut output_bytes = 0;
    for (input, (reference, reference_rank)) in inputs.iter().zip(expected) {
        let (rows, rank, part, hit) = engine.apply(input);
        assert_eq!(rank, *reference_rank);
        assert_eq!(&rows, reference);
        hits += usize::from(hit);
        work.high_word_xors += part.high_word_xors;
        work.low_word_xors += part.low_word_xors;
        work.tail_low_word_xors += part.tail_low_word_xors;
        work.row_swaps += part.row_swaps;
        output_bytes += hash_output(&mut hasher, &rows, rank);
    }
    Some(json!({
        "retained_bytes":retained_bytes,
        "setup_work":setup_work,
        "output_bytes":output_bytes,
        "output_digest":format!("{:x}",hasher.finalize()),
        "hits":hits,
        "fallbacks":inputs.len()-hits,
        "work":work
    }))
}

fn verify_sample(
    line: &Value,
    repetition: usize,
    order: usize,
    arm: Arm,
    expected_name: &str,
    context: &SampleContext<'_>,
) -> u64 {
    let cell = context.cell;
    let inputs = context.inputs;
    let expected = context.expected;
    let cap = context.cap;
    assert_eq!(line["kind"], "sample");
    assert_eq!(line["cell"], cell);
    assert_eq!(number(line, "repetition"), repetition as u64);
    assert_eq!(number(line, "order"), order as u64);
    assert_eq!(line["arm"], expected_name);
    let semantic = expected_sample(arm, inputs, expected, cap).expect("arm applicable");
    for key in [
        "retained_bytes",
        "setup_work",
        "output_bytes",
        "output_digest",
        "hits",
        "fallbacks",
        "work",
    ] {
        assert_eq!(
            line[key], semantic[key],
            "sample {cell} {expected_name} {key}"
        );
    }
    let total = number(line, "total_ns");
    let setup = number(line, "setup_ns");
    let apply = number(line, "apply_ns");
    let validate = number(line, "validate_ns");
    assert!(total > 0 && total >= setup + apply + validate);
    total
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
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    assert!(!sorted.is_empty());
    let index = ((sorted.len() as f64 * fraction).ceil() as usize)
        .saturating_sub(1)
        .min(sorted.len() - 1);
    sorted[index]
}

fn bootstrap_lower(values: &[f64], seed: u64) -> f64 {
    let mut state = seed;
    let mut medians = Vec::with_capacity(4000);
    for _ in 0..4000 {
        let mut sample = Vec::with_capacity(values.len());
        for _ in values {
            sample.push(values[(next(&mut state) % values.len() as u64) as usize]);
        }
        medians.push(median(&sample));
    }
    quantile(&medians, 0.025)
}

fn resource_receipt(phase: &str, path: &Path, campaign_seconds: u64) -> Value {
    let lines = value_lines(path);
    assert_eq!(lines.len(), 1, "one reserved worker required");
    let receipt = &lines[0];
    assert_eq!(receipt["schema"], "isolated-bench/1");
    assert_eq!(receipt["mode"], "run");
    assert_eq!(receipt["label"], format!("graded-tail/{phase}"));
    assert_eq!(receipt["reserved_cpus"].as_array().unwrap().len(), 1);
    assert_eq!(receipt["run"]["contended"], false);
    assert_eq!(receipt["run"]["exit_status"], 0);
    let settle = receipt["preflight"]["settle"]["seconds"].as_f64().unwrap();
    let other_settle = receipt["preflight"]["settle"]["other_cpu_seconds"]
        .as_f64()
        .unwrap();
    assert!(settle > 0.0 && other_settle / settle <= 0.10);
    let psi = receipt["preflight"]["conditions"]["psi_cpu"]["some"]["avg10"]
        .as_f64()
        .unwrap();
    assert!(psi <= 5.0);
    let wall = receipt["run"]["wall_seconds"].as_f64().unwrap();
    let other = receipt["run"]["other_cpu_seconds"].as_f64().unwrap();
    assert!(wall > 0.0 && wall <= campaign_seconds as f64 && other / wall <= 0.10);
    assert_eq!(receipt["host"]["machine"], "aarch64");
    json!({
        "wall_seconds":wall,
        "other_cpu_seconds":other,
        "psi_some_avg10":psi,
        "max_rss_kib":number(&receipt["run"],"max_rss_kib"),
        "contended":false,
        "reserved_cpus":receipt["reserved_cpus"]
    })
}

fn readiness_receipt(path: &Path) -> Value {
    let readiness = read_json(path);
    assert_eq!(readiness["accepted"], true);
    let samples = readiness["samples"].as_array().expect("readiness samples");
    assert!(!samples.is_empty() && samples.len() <= 30);
    for (index, sample) in samples.iter().enumerate() {
        assert_eq!(number(sample, "index"), index as u64);
        let psi = sample["psi_some_avg10"].as_f64().unwrap();
        assert_eq!(sample["accepted"], psi <= 5.0);
    }
    assert_eq!(samples.last().unwrap()["accepted"], true);
    readiness
}

fn compute_run(phase: &str, raw: &Path, conditions: &Path) -> Value {
    assert_eq!(hash(PROTOCOL_BYTES), FROZEN_PROTOCOL_SHA256);
    let protocol: Protocol = serde_json::from_slice(PROTOCOL_BYTES).expect("frozen protocol");
    let (split, seeds) = match phase {
        "discovery" => ("discovery", protocol.discovery_seeds),
        "full" => ("holdout", protocol.holdout_seeds),
        _ => panic!("phase"),
    };
    let lines = value_lines(raw);
    let mut cursor = 0;
    let mut take = || {
        let line = lines
            .get(cursor)
            .unwrap_or_else(|| panic!("truncated raw at {cursor}"));
        cursor += 1;
        line
    };
    let header = take();
    assert_eq!(header["kind"], "campaign");
    assert_eq!(header["phase"], phase);
    assert_eq!(header["protocol_sha256"], hash(PROTOCOL_BYTES));
    assert_eq!(header["source_sha256"], hash(SOURCE_BYTES));
    assert_eq!(header["verifier_sha256"], hash(VERIFY_BYTES));
    let mut groups = BTreeMap::<(u8, String), Vec<f64>>::new();
    let mut noise = BTreeMap::<(u8, String), Vec<f64>>::new();
    let mut fixtures = 0;
    let mut samples = 0;
    let mut verified_outputs = 0;
    for &n in &protocol.variables {
        for &seed in &seeds {
            for family in &protocol.families {
                for &batch in &protocol.batches {
                    let inputs = assignments(n, seed, batch, family);
                    let basis = Basis::new(n);
                    assert!(row_count(&inputs[0]) <= protocol.row_cap);
                    assert!(basis.high.len() + basis.low.len() <= protocol.column_cap);
                    let mut expected = Vec::with_capacity(batch);
                    for input in &inputs {
                        let matrix = construct_full(input, &basis);
                        oracle_check_input(input, &matrix, &basis);
                        let (rows, rank, _) = reduce_constructed(matrix, &basis);
                        expected.push((rows, rank));
                    }
                    let cell = format!("n{n}-{split}-{seed}-{family}-b{batch}");
                    let expected_fixture = json!({
                        "kind":"fixture", "cell":cell, "n":n, "seed":seed,
                        "batch":batch, "family":family,
                        "output_ranks":expected.iter().map(|(_,r)|*r).collect::<Vec<_>>(),
                        "quadratic_support":inputs[0].quadratic,
                        "affine_masks":inputs.iter().map(|input|&input.affine).collect::<Vec<_>>(),
                        "escape_indices":inputs.iter().enumerate().filter_map(|(i,input)|
                            (input.quadratic!=inputs[0].quadratic).then_some(i)).collect::<Vec<_>>()
                    });
                    assert_eq!(take(), &expected_fixture, "fixture {cell}");
                    fixtures += 1;
                    let context = SampleContext {
                        cell: &cell,
                        inputs: &inputs,
                        expected: &expected,
                        cap: protocol.retained_context_cap_bytes,
                    };
                    for repetition in 0..protocol.repetitions {
                        let mut aa = [0u64; 2];
                        for (order, value) in aa.iter_mut().enumerate() {
                            *value = verify_sample(
                                take(),
                                repetition,
                                order,
                                Arm::RankedDense,
                                if order == 0 { "aa_a" } else { "aa_b" },
                                &context,
                            );
                            samples += 1;
                            verified_outputs += batch;
                        }
                        let mut candidate = None;
                        let mut best_reference = u64::MAX;
                        for order in 0..Arm::ALL.len() {
                            let arm = Arm::ALL[(repetition + order) % Arm::ALL.len()];
                            let cost = verify_sample(
                                take(),
                                repetition,
                                order + 2,
                                arm,
                                arm.name(),
                                &context,
                            );
                            if arm == Arm::Graded {
                                candidate = Some(cost);
                            } else {
                                best_reference = best_reference.min(cost);
                            }
                            samples += 1;
                            verified_outputs += batch;
                        }
                        if batch == 64
                            && n >= 12
                            && (family == "independent_affine" || family == "walk_affine")
                        {
                            let key = (n, family.clone());
                            groups
                                .entry(key.clone())
                                .or_default()
                                .push(best_reference as f64 / candidate.unwrap() as f64);
                            let aa_ratio = aa[0] as f64 / aa[1] as f64;
                            noise
                                .entry(key)
                                .or_default()
                                .push(aa_ratio.max(1.0 / aa_ratio));
                        }
                    }
                }
            }
        }
    }
    assert_eq!(cursor, lines.len(), "trailing or missing raw records");
    assert_eq!(fixtures, 200);
    assert_eq!(groups.len(), 8);
    let mut group_rows = Vec::new();
    let threshold = if phase == "discovery" { 1.5 } else { 2.0 };
    let mut gate_pass = true;
    for ((n, family), ratios) in groups {
        assert_eq!(ratios.len(), 20);
        let aa = &noise[&(n, family.clone())];
        let family_code = if family == "independent_affine" { 1 } else { 2 };
        let lower = bootstrap_lower(&ratios, 20261003 ^ ((n as u64) << 32) ^ family_code);
        let floor = quantile(aa, 0.975);
        let pass = lower > threshold && lower > floor;
        gate_pass &= pass;
        group_rows.push(json!({
            "n":n,"family":family,"batch":64,
            "observations":ratios.len(),
            "paired_median_ratio":median(&ratios),
            "bootstrap_lower_95":lower,
            "aa_noise_floor_975":floor,
            "pass":pass
        }));
    }
    let resource = resource_receipt(phase, conditions, protocol.campaign_seconds);
    let readiness = readiness_receipt(&raw.with_file_name("readiness.json"));
    json!({
        "schema_version":1,
        "phase":phase,
        "protocol_sha256":hash(PROTOCOL_BYTES),
        "source_sha256":hash(SOURCE_BYTES),
        "verifier_sha256":hash(VERIFY_BYTES),
        "fixtures":fixtures,
        "samples":samples,
        "verified_outputs":verified_outputs,
        "resource":resource,
        "readiness":readiness,
        "groups":group_rows,
        "gate_pass":gate_pass,
        "full_solver_cost":null,
        "relation_yield":null,
        "rho_ratio":null
    })
}

pub(super) fn verify_run(phase: &str, raw: &str, conditions: &str, results: &str) {
    let output = Path::new(results);
    assert!(!output.exists(), "refuse to overwrite results");
    if phase == "full" {
        validate_binding(&output.with_file_name("binding.json"));
    }
    let report = compute_run(phase, Path::new(raw), Path::new(conditions));
    fs::write(output, serde_json::to_vec_pretty(&report).unwrap()).expect("write results");
}

fn validate_binding(path: &Path) {
    let binding = read_json(path);
    assert_eq!(binding["source_sha256"], hash(SOURCE_BYTES));
    assert_eq!(binding["verifier_sha256"], hash(VERIFY_BYTES));
    assert_eq!(binding["protocol_sha256"], hash(PROTOCOL_BYTES));
    for key in ["discovery_manifest_sha256", "discovery_results_sha256"] {
        assert_eq!(binding[key].as_str().unwrap().len(), 64);
    }
}

pub(super) fn seal(dir: &str) {
    let root = Path::new(dir);
    let manifest = root.join("manifest.json");
    assert!(!manifest.exists(), "refuse to overwrite manifest");
    let mut files = BTreeMap::new();
    for entry in fs::read_dir(root).expect("bundle directory") {
        let entry = entry.expect("directory entry");
        let kind = entry.file_type().expect("file type");
        assert!(!kind.is_symlink(), "no bundle symlinks");
        if !kind.is_file() {
            continue;
        }
        let name = entry.file_name().into_string().expect("UTF-8 file name");
        assert_ne!(name, "manifest.json");
        files.insert(name, hash(&fs::read(entry.path()).expect("member")));
    }
    assert!(
        files.contains_key("worker.rs")
            && files.contains_key("verify.rs")
            && files.contains_key("protocol.json")
    );
    fs::write(
        manifest,
        serde_json::to_vec_pretty(&json!({
            "schema_version":1,"files":files
        }))
        .unwrap(),
    )
    .expect("write manifest");
}

pub(super) fn verify_bundle(dir: &str) -> Value {
    let root = Path::new(dir);
    let manifest = read_json(&root.join("manifest.json"));
    assert_eq!(manifest["schema_version"], 1);
    let files = manifest["files"].as_object().expect("manifest files");
    assert!(files.len() >= 9);
    for (name, expected) in files {
        assert!(!name.contains('/') && !name.contains(".."));
        assert_eq!(
            hash(&fs::read(root.join(name)).expect("manifest member")),
            expected.as_str().unwrap(),
            "member {name}"
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
    if root.join("results.json").exists() {
        let old = read_json(&root.join("results.json"));
        let phase = old["phase"].as_str().expect("phase");
        if phase == "full" {
            validate_binding(&root.join("binding.json"));
        }
        let expected = compute_run(
            phase,
            &root.join("raw.jsonl"),
            &root.join("conditions.jsonl"),
        );
        assert_eq!(old, expected, "bundle replay differs");
        old
    } else {
        panic!("missing result or failure");
    }
}

pub(super) fn check_discovery(bundle: &str, binding: &str) {
    let result = verify_bundle(bundle);
    assert_eq!(result["phase"], "discovery");
    assert_eq!(result["gate_pass"], true, "discovery gate failed");
    assert_eq!(result["source_sha256"], hash(SOURCE_BYTES));
    assert_eq!(result["verifier_sha256"], hash(VERIFY_BYTES));
    assert_eq!(result["protocol_sha256"], hash(PROTOCOL_BYTES));
    assert!(!Path::new(binding).exists());
    let root = Path::new(bundle);
    fs::write(
        binding,
        serde_json::to_vec_pretty(&json!({
        "source_sha256":hash(SOURCE_BYTES),
        "verifier_sha256":hash(VERIFY_BYTES),
            "protocol_sha256":hash(PROTOCOL_BYTES),
            "discovery_manifest_sha256":hash(&fs::read(root.join("manifest.json")).unwrap()),
            "discovery_results_sha256":hash(&fs::read(root.join("results.json")).unwrap())
        }))
        .unwrap(),
    )
    .expect("binding");
}

pub(super) fn failure(dir: &str, reason: &str) {
    let path = Path::new(dir).join("failure.json");
    assert!(!path.exists(), "refuse to overwrite failure");
    fs::write(
        path,
        serde_json::to_vec_pretty(&json!({
            "reason":reason,
            "performance_admitted":false,
            "full_solver_cost":null,
            "rho_ratio":null
        }))
        .unwrap(),
    )
    .expect("write failure");
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
    let target = PathBuf::from(path);
    assert!(!target.exists(), "refuse to overwrite readiness");
    let mut samples = Vec::new();
    for index in 0..30 {
        let psi = psi_some_avg10();
        samples.push(json!({"index":index,"psi_some_avg10":psi,"accepted":psi<=5.0}));
        if psi <= 5.0 {
            fs::write(
                &target,
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
        &target,
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
    fn reported_medians_survive_json_roundtrip_exactly() {
        for value in [1.0247742717975405_f64, 1.4431176749304875_f64] {
            let bytes = serde_json::to_string(&value).unwrap();
            let restored: f64 = serde_json::from_str(&bytes).unwrap();
            assert_eq!(value.to_bits(), restored.to_bits());
        }
    }

    #[test]
    fn native_sample_replay_rejects_changed_digest() {
        let inputs = assignments(8, 41, 4, "support_escape");
        let basis = Basis::new(8);
        let expected = inputs
            .iter()
            .map(|input| {
                let raw = construct_full(input, &basis);
                oracle_check_input(input, &raw, &basis);
                let (rows, rank, _) = reduce_constructed(raw, &basis);
                (rows, rank)
            })
            .collect::<Vec<_>>();
        let context = SampleContext {
            cell: "development",
            inputs: &inputs,
            expected: &expected,
            cap: 64 * 1024 * 1024,
        };
        let observation = sample(Arm::Graded, "graded", 0, 2, &context).unwrap();
        let mut value = serde_json::to_value(observation).unwrap();
        verify_sample(&value, 0, 2, Arm::Graded, "graded", &context);
        value["output_digest"] = json!("00");
        assert!(std::panic::catch_unwind(|| verify_sample(
            &value,
            0,
            2,
            Arm::Graded,
            "graded",
            &context
        ))
        .is_err());
    }

    #[test]
    fn resource_gate_rejects_contention() {
        let path = std::env::temp_dir().join(format!(
            "graded-resource-test-{}-{}.jsonl",
            std::process::id(),
            17
        ));
        let mut receipt = json!({
            "schema":"isolated-bench/1", "mode":"run",
            "label":"graded-tail/discovery", "reserved_cpus":[3],
            "host":{"machine":"aarch64"},
            "preflight":{"settle":{"seconds":2.0,"other_cpu_seconds":0.02},
                "conditions":{"psi_cpu":{"some":{"avg10":1.0}}}},
            "run":{"contended":false,"exit_status":0,"wall_seconds":10.0,
                "other_cpu_seconds":0.2,"max_rss_kib":8192}
        });
        fs::write(&path, format!("{receipt}\n")).unwrap();
        resource_receipt("discovery", &path, 1200);
        receipt["run"]["contended"] = json!(true);
        fs::write(&path, format!("{receipt}\n")).unwrap();
        assert!(std::panic::catch_unwind(|| resource_receipt("discovery", &path, 1200)).is_err());
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn sealed_failure_with_partial_result_stays_unadmitted() {
        let root = std::env::temp_dir().join(format!(
            "graded-sealed-failure-{}-{}",
            std::process::id(),
            29
        ));
        assert!(!root.exists());
        fs::create_dir(&root).unwrap();
        fs::write(root.join("worker.rs"), SOURCE_BYTES).unwrap();
        fs::write(root.join("verify.rs"), VERIFY_BYTES).unwrap();
        fs::write(root.join("protocol.json"), PROTOCOL_BYTES).unwrap();
        for name in [
            "Cargo.toml",
            "Cargo.lock",
            "PROTOCOL.md",
            "README.md",
            "run.sh",
        ] {
            fs::write(root.join(name), b"retained").unwrap();
        }
        fs::write(root.join("results.json"), b"{\"gate_pass\":true}").unwrap();
        failure(root.to_str().unwrap(), "post-seal-replay-failed");
        seal(root.to_str().unwrap());
        let replay = verify_bundle(root.to_str().unwrap());
        assert_eq!(replay["performance_admitted"], false);
        assert_eq!(replay["reason"], "post-seal-replay-failed");
        fs::remove_dir_all(root).unwrap();
    }
}
