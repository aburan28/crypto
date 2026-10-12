//! Native replay of the frozen cell order, complete matrices and resource gate.
use super::*;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

fn field<'a>(row: &'a Value, key: &str) -> &'a Value {
    row.get(key).unwrap_or_else(|| panic!("missing {key}"))
}
fn integer(row: &Value, key: &str) -> u64 {
    field(row, key)
        .as_u64()
        .unwrap_or_else(|| panic!("invalid integer {key}"))
}
fn string<'a>(row: &'a Value, key: &str) -> &'a str {
    field(row, key)
        .as_str()
        .unwrap_or_else(|| panic!("invalid string {key}"))
}
fn number(row: &Value, key: &str) -> f64 {
    let value = field(row, key).as_f64().expect("finite resource number");
    assert!(value.is_finite(), "nonfinite {key}");
    value
}
fn sha256(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
fn next_json(lines: &mut impl Iterator<Item = std::io::Result<String>>) -> Value {
    let line = lines
        .next()
        .expect("missing JSONL row")
        .expect("read JSONL row");
    serde_json::from_str(&line).expect("valid JSONL row")
}
fn median(values: &[f64]) -> f64 {
    assert!(!values.is_empty());
    let mut values = values.to_vec();
    values.sort_by(f64::total_cmp);
    if values.len() % 2 == 1 {
        values[values.len() / 2]
    } else {
        (values[values.len() / 2 - 1] + values[values.len() / 2]) / 2.0
    }
}
fn interval(values: &[f64], seed: u64) -> [f64; 2] {
    assert_eq!(values.len(), 10);
    let mut state = seed;
    let mut medians = Vec::with_capacity(4000);
    for _ in 0..4000 {
        let draw: Vec<_> = (0..values.len())
            .map(|_| values[(next(&mut state) as usize) % values.len()])
            .collect();
        medians.push(median(&draw));
    }
    medians.sort_by(f64::total_cmp);
    [medians[100], medians[3899]]
}
fn oracle_digest(expected: &[Matrix]) -> String {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"boolean-matrix-v1");
    for matrix in expected {
        digest_matrix(&mut hasher, matrix);
    }
    hasher.finalize().to_hex().to_string()
}
fn independent_hit(
    variant: usize,
    input: &Input,
    base: &Input,
    envelope: &Input,
    actual: &Matrix,
    base_matrix: &Matrix,
) -> bool {
    match variant {
        0 => false,
        1 => actual.columns == base_matrix.columns,
        2 | 3 => input == base,
        4..=6 => {
            input.n == envelope.n
                && input.degree == envelope.degree
                && input.active == envelope.active
                && input.polys.len() == envelope.polys.len()
                && input
                    .polys
                    .iter()
                    .zip(&envelope.polys)
                    .all(|(poly, allowed)| {
                        poly.iter().all(|term| allowed.binary_search(term).is_ok())
                    })
        }
        _ => unreachable!(),
    }
}
fn expected_counts(
    variant: usize,
    inputs: &[Input],
    expected: &[Matrix],
    base: &Input,
    envelope: &Input,
    base_matrix: &Matrix,
) -> (u64, u64, u64) {
    let mut counts = (0, 0, 0);
    for (input, matrix) in inputs.iter().zip(expected) {
        let hit = independent_hit(variant, input, base, envelope, matrix, base_matrix);
        counts.0 += u64::from(hit);
        counts.1 += u64::from(variant != 0 && !hit);
        counts.2 += u64::from(hit && input != base);
    }
    counts
}
fn verify_sample(
    sample: &Value,
    cell: &str,
    rep: usize,
    order: usize,
    name: &str,
    variant: usize,
    base: &Input,
    envelope: &Input,
    inputs: &[Input],
    expected: &[Matrix],
    base_matrix: &Matrix,
    digest: &str,
) -> u64 {
    assert_eq!(string(sample, "type"), "sample");
    assert_eq!(string(sample, "cell"), cell);
    assert_eq!(string(sample, "variant"), name);
    assert_eq!(integer(sample, "rep"), rep as u64);
    assert_eq!(integer(sample, "order"), order as u64);
    assert_eq!(integer(sample, "verified_outputs"), inputs.len() as u64);
    assert_eq!(string(sample, "output_digest"), digest);
    let payload: usize = expected.iter().map(Matrix::payload).sum();
    assert_eq!(integer(sample, "output_bytes"), payload as u64);
    let retained = Cache::compile(variant, base, envelope).retained();
    assert_eq!(integer(sample, "retained_bytes"), retained as u64);
    let counts = expected_counts(variant, inputs, expected, base, envelope, base_matrix);
    assert_eq!(integer(sample, "hits"), counts.0);
    assert_eq!(integer(sample, "fallbacks"), counts.1);
    assert_eq!(integer(sample, "changed_hits"), counts.2);
    let total = integer(sample, "total_ns");
    assert!(total > 0);
    assert!(integer(sample, "apply_ns") > 0 && integer(sample, "validation_ns") > 0);
    assert!(
        integer(sample, "setup_ns")
            + integer(sample, "apply_ns")
            + integer(sample, "validation_ns")
            <= total
    );
    total
}

fn verify_resource(path: &str, phase: &str, protocol_path: &str, worker_cap: f64) -> Value {
    let raw = std::fs::read_to_string(path).expect("resource conditions");
    let mut lines = raw.lines();
    let conditions: Value = serde_json::from_str(lines.next().expect("one resource row")).unwrap();
    assert!(lines.next().is_none(), "multiple resource rows");
    assert_eq!(string(&conditions, "schema"), "isolated-bench/1");
    assert_eq!(string(&conditions, "mode"), "run");
    let pre = field(&conditions, "preflight");
    let settle = field(pre, "settle");
    assert_eq!(number(settle, "seconds"), 2.0);
    assert!((0.0..=0.2).contains(&number(settle, "other_cpu_seconds")));
    for key in ["psi_cpu", "psi_memory"] {
        assert!((0.0..=5.0).contains(&number(
            field(field(pre, "conditions"), key).get("some").unwrap(),
            "avg10"
        )));
    }
    let run = field(&conditions, "run");
    assert_eq!(integer(run, "exit_status"), 0);
    assert_eq!(field(run, "contended").as_bool(), Some(false));
    let wall = number(run, "wall_seconds");
    assert!(wall > 0.0 && wall <= worker_cap);
    assert!((0.0..=0.1 * wall).contains(&number(run, "other_cpu_seconds")));
    let cpus = field(&conditions, "reserved_cpus")
        .as_array()
        .expect("reserved cpus");
    assert!(!cpus.is_empty() && field(run, "cpus") == field(&conditions, "reserved_cpus"));
    let host = field(&conditions, "host");
    assert!(integer(host, "logical_cpus") > cpus.len() as u64);
    assert!(!string(host, "machine").is_empty() && !string(host, "kernel").is_empty());
    for when in ["before", "after"] {
        for key in ["psi_cpu", "psi_memory"] {
            assert!(field(field(run, when), key).is_object());
        }
    }
    for key in [
        "voluntary_switches",
        "involuntary_switches",
        "minor_faults",
        "major_faults",
        "max_rss_kib",
    ] {
        let _ = integer(run, key);
    }
    let command = field(run, "command").as_array().expect("worker command");
    assert_eq!(command.len(), 4);
    assert!(command[0]
        .as_str()
        .unwrap()
        .ends_with("boolean-envelope-indexed"));
    assert_eq!(command[1].as_str(), Some("--campaign"));
    assert_eq!(command[2].as_str(), Some(phase));
    assert_eq!(command[3].as_str(), Some(protocol_path));
    conditions
}

pub(super) fn run(
    phase: &str,
    protocol_path: &str,
    raw_path: &str,
    conditions_path: &str,
    result_path: &str,
) {
    assert!(["discovery", "full"].contains(&phase));
    let protocol_bytes = std::fs::read(protocol_path).expect("frozen protocol");
    assert_eq!(protocol_bytes, include_bytes!("protocol.json"));
    let protocol: StudyProtocol = serde_json::from_slice(&protocol_bytes).unwrap();
    assert_eq!(protocol.variables, [6, 8, 10, 12]);
    assert_eq!(
        protocol.families,
        ["repeat", "coefficients", "degree_cycle", "escape"]
    );
    assert_eq!(protocol.batches, [1, 4, 16, 64]);
    assert_eq!(protocol.repetitions, 10);
    assert_eq!(protocol.bootstrap.resamples, 4000);
    assert_eq!(protocol.dramatic_gate.groups, 32);
    assert_eq!(protocol.dramatic_gate.lower_bound_greater_than, 2.0);
    assert_eq!(protocol.noise.floor_quantile, 0.975);
    let names: Vec<_> = protocol
        .reference_arms
        .iter()
        .chain(&protocol.candidate_arms)
        .map(String::as_str)
        .collect();
    assert_eq!(names, ARM_NAMES);
    let (split, seeds) = if phase == "discovery" {
        ("discovery", &protocol.discovery_seeds)
    } else {
        ("holdout", &protocol.holdout_seeds)
    };
    let replay_started = Instant::now();
    let resource = verify_resource(
        conditions_path,
        phase,
        protocol_path,
        protocol.limits.worker_seconds as f64,
    );
    let file = File::open(raw_path).expect("raw campaign stream");
    let mut lines = BufReader::new(file).lines();
    let header = next_json(&mut lines);
    assert_eq!(string(&header, "type"), "campaign");
    assert_eq!(string(&header, "phase"), phase);
    assert_eq!(
        string(&header, "protocol_digest"),
        blake3::hash(&protocol_bytes).to_hex().as_str()
    );
    let mut cells = Vec::new();
    let mut gates = Vec::new();
    let mut aa_observations = 0usize;
    let mut ab_observations = 0usize;
    let mut verified_outputs = 0usize;
    let mut decisions = BTreeMap::new();
    for arm in &protocol.candidate_arms {
        decisions.insert(arm.clone(), (0usize, 0usize, 0usize));
    }
    for &n in &protocol.variables {
        for &seed in seeds {
            for family in &protocol.families {
                for &batch in &protocol.batches {
                    let cell = format!("n{n}-{split}-{seed}-{family}-b{batch}");
                    let fixture_row = next_json(&mut lines);
                    assert_eq!(string(&fixture_row, "type"), "fixture");
                    assert_eq!(string(&fixture_row, "cell"), cell);
                    assert_eq!(integer(&fixture_row, "n"), u64::from(n));
                    assert_eq!(integer(&fixture_row, "seed"), seed);
                    assert_eq!(integer(&fixture_row, "batch"), batch as u64);
                    assert_eq!(string(&fixture_row, "family"), family);
                    assert_eq!(integer(&fixture_row, "degree"), 3);
                    let base = fixture(n, seed);
                    let envelope = declared_envelope(&base);
                    let inputs = assignments(&base, &envelope, seed, batch, family);
                    let polys: Vec<_> = inputs.iter().map(|input| &input.polys).collect();
                    assert_eq!(integer(&fixture_row, "active"), u64::from(base.active));
                    assert_eq!(field(&fixture_row, "base"), &json!(base.polys));
                    assert_eq!(field(&fixture_row, "envelope"), &json!(envelope.polys));
                    assert_eq!(field(&fixture_row, "inputs"), &json!(polys));
                    assert_eq!(
                        integer(&fixture_row, "newly_required_multipliers"),
                        inputs
                            .iter()
                            .map(|input| newly_required(&base, input))
                            .sum::<usize>() as u64
                    );
                    let expected: Vec<_> = inputs
                        .iter()
                        .map(|input| oracle(input, CAPS).unwrap())
                        .collect();
                    let base_matrix = oracle(&base, CAPS).unwrap();
                    let digest = oracle_digest(&expected);
                    let mut noise = Vec::new();
                    let mut totals = vec![vec![0u64; ARM_NAMES.len()]; protocol.repetitions];
                    for rep in 0..protocol.repetitions {
                        let first = next_json(&mut lines);
                        let second = next_json(&mut lines);
                        let a = verify_sample(
                            &first,
                            &cell,
                            rep,
                            0,
                            "aa_a",
                            3,
                            &base,
                            &envelope,
                            &inputs,
                            &expected,
                            &base_matrix,
                            &digest,
                        );
                        let b = verify_sample(
                            &second,
                            &cell,
                            rep,
                            1,
                            "aa_b",
                            3,
                            &base,
                            &envelope,
                            &inputs,
                            &expected,
                            &base_matrix,
                            &digest,
                        );
                        assert!(a > 0 && b > 0);
                        noise.push((a as f64 / b as f64).max(b as f64 / a as f64));
                        aa_observations += 2;
                        verified_outputs += 2 * batch;
                        for order in 0..ARM_NAMES.len() {
                            let variant = (rep + order) % ARM_NAMES.len();
                            let sample = next_json(&mut lines);
                            totals[rep][variant] = verify_sample(
                                &sample,
                                &cell,
                                rep,
                                order,
                                ARM_NAMES[variant],
                                variant,
                                &base,
                                &envelope,
                                &inputs,
                                &expected,
                                &base_matrix,
                                &digest,
                            );
                            ab_observations += 1;
                            verified_outputs += batch;
                        }
                    }
                    noise.sort_by(f64::total_cmp);
                    let floor_index =
                        (protocol.noise.floor_quantile * (noise.len() - 1) as f64).ceil() as usize;
                    let noise_floor = noise[floor_index];
                    let reference: Vec<u64> = totals
                        .iter()
                        .map(|rep| *rep[..5].iter().min().unwrap())
                        .collect();
                    let mut medians = BTreeMap::new();
                    for (index, arm) in ARM_NAMES.iter().enumerate() {
                        medians.insert(
                            *arm,
                            median(
                                &totals
                                    .iter()
                                    .map(|row| row[index] as f64)
                                    .collect::<Vec<_>>(),
                            ),
                        );
                    }
                    if ["coefficients", "degree_cycle"].contains(&family.as_str())
                        && [16, 64].contains(&batch)
                    {
                        for index in 5..ARM_NAMES.len() {
                            let ratios: Vec<_> = (0..protocol.repetitions)
                                .map(|rep| reference[rep] as f64 / totals[rep][index] as f64)
                                .collect();
                            let mut tag = [0u8; 8];
                            tag.copy_from_slice(&blake3::hash(cell.as_bytes()).as_bytes()[..8]);
                            let interval = interval(
                                &ratios,
                                protocol.bootstrap.seed ^ u64::from_le_bytes(tag) ^ index as u64,
                            );
                            let dramatic = interval[0]
                                > protocol
                                    .dramatic_gate
                                    .lower_bound_greater_than
                                    .max(noise_floor);
                            let incremental = interval[0] > 1.0f64.max(noise_floor);
                            let counts = decisions.get_mut(ARM_NAMES[index]).unwrap();
                            counts.0 += usize::from(dramatic);
                            counts.1 += usize::from(incremental);
                            counts.2 += 1;
                            gates.push(json!({"cell":cell,"candidate":ARM_NAMES[index],"median":median(&ratios),"ci95":interval,"noise_floor":noise_floor,"dramatic":dramatic,"incremental":incremental,"observations":ratios.len()}));
                        }
                    }
                    cells.push(json!({"cell":cell,"noise_floor":noise_floor,"medians_ns":medians,"verified_outputs":batch * (ARM_NAMES.len() + 2) * protocol.repetitions}));
                }
            }
        }
    }
    assert!(lines.next().is_none(), "unexpected trailing campaign rows");
    assert_eq!(cells.len(), 128);
    assert_eq!(aa_observations, 2560);
    assert_eq!(ab_observations, 8960);
    assert_eq!(gates.len(), 2 * protocol.dramatic_gate.groups);
    let decisions: BTreeMap<_, _> = decisions.into_iter().map(|(arm, (dramatic, incremental, total))| {
        assert_eq!(total, protocol.dramatic_gate.groups);
        let pass = phase == "full" && dramatic == total;
        (arm, json!({"dramatic_groups":dramatic,"incremental_groups":incremental,"groups":total,"dramatic_pass":pass}))
    }).collect();
    let raw = std::fs::read(raw_path).unwrap();
    let conditions = std::fs::read(conditions_path).unwrap();
    let result_directory = Path::new(result_path).parent().unwrap();
    let binding = if phase == "full" {
        serde_json::from_slice::<Value>(
            &std::fs::read(result_directory.join("binding.json")).unwrap(),
        )
        .unwrap()
    } else {
        Value::Null
    };
    assert!(
        number(field(&resource, "run"), "wall_seconds") + replay_started.elapsed().as_secs_f64()
            <= protocol.limits.campaign_seconds as f64,
        "campaign/replay cap"
    );
    let result = json!({"schema_version":1,"phase":phase,"qualified":true,"all_complete":true,"cells":cells,"gates":gates,"decisions":decisions,"aa_observations":aa_observations,"ab_observations":ab_observations,"verified_outputs":verified_outputs,"protocol_sha256":sha256(&protocol_bytes),"raw_sha256":sha256(&raw),"conditions_sha256":sha256(&conditions),"discovery_binding":binding,"resource":resource,"full_ic_cost":Value::Null,"calibrated_operation_ratio":Value::Null,"rho_ratio":Value::Null});
    std::fs::write(result_path, serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    let report_path = Path::new(result_path).with_file_name("RESULT.md");
    std::fs::write(report_path, render(&result)).unwrap();
}

fn render(result: &Value) -> String {
    let phase = string(result, "phase");
    let mut note = format!("# Native indexed-envelope {phase} readback\n\n");
    note.push_str("This is a complete, resource-qualified **Boolean matrix-construction stage** comparison. It does not measure a complete index-calculus solve or a rho boundary. Full-IC and rho costs remain null.\n\n");
    note.push_str(&format!("All {} cells were complete: {} A/B observations, {} A/A observations, and {} oracle-verified output matrices. Per-call cold totals include compilation, every application, validation, output allocation and destruction.\n\n", field(result,"cells").as_array().unwrap().len(), integer(result,"ab_observations"), integer(result,"aa_observations"), integer(result,"verified_outputs")));
    note.push_str("| Candidate | Dramatic groups | Incremental groups | Full primary gate |\n|---|---:|---:|---|\n");
    for (arm, decision) in field(result, "decisions").as_object().unwrap() {
        let verdict = if phase == "discovery" {
            "discovery only"
        } else if field(decision, "dramatic_pass").as_bool() == Some(true) {
            "PASS; confirmation required"
        } else {
            "FAIL"
        };
        note.push_str(&format!(
            "| {arm} | {}/{} | {}/{} | {verdict} |\n",
            integer(decision, "dramatic_groups"),
            integer(decision, "groups"),
            integer(decision, "incremental_groups"),
            integer(decision, "groups")
        ));
    }
    note.push_str("\nThe primary ratio for each group uses the pointwise fastest old control, a paired 95% median-ratio interval, and the A/A noise floor. Every declared group is retained below.\n\n| Cell | Candidate | Median ratio | 95% lower | 95% upper | A/A floor | >2x gate |\n|---|---|---:|---:|---:|---:|---|\n");
    if phase == "full" {
        let binding = field(result, "discovery_binding");
        note.push_str(&format!(
            "Discovery binding: manifest SHA-256 `{}`, source commit `{}`.\n\n",
            string(binding, "manifest_sha256"),
            string(binding, "source_head")
        ));
    }
    for gate in field(result, "gates").as_array().unwrap() {
        let ci = field(gate, "ci95").as_array().unwrap();
        note.push_str(&format!(
            "| {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {} |\n",
            string(gate, "cell"),
            string(gate, "candidate"),
            number(gate, "median"),
            ci[0].as_f64().unwrap(),
            ci[1].as_f64().unwrap(),
            number(gate, "noise_floor"),
            if field(gate, "dramatic").as_bool() == Some(true) {
                "PASS"
            } else {
                "FAIL"
            }
        ));
    }
    let run = field(field(result, "resource"), "run");
    let host = field(field(result, "resource"), "host");
    note.push_str(&format!("\nHost: {} on {} logical CPUs; reserved CPUs {:?}. Whole-campaign wall time {:.3} s, other-process CPU {:.3} s, peak worker RSS {} KiB. The resource receipt covers the whole campaign; individual short calls have separate timers and A/A pairs, not individual resource certifications.\n\n",string(host,"machine"),integer(host,"logical_cpus"),field(field(result,"resource"),"reserved_cpus"),number(run,"wall_seconds"),number(run,"other_cpu_seconds"),integer(run,"max_rss_kib")));
    note.push_str("The frozen raw JSONL, input source, binaries, host records, resource receipt, native verifier result and this report are bound by the artifact manifest. A positive full result requires a further unchanged-source run on unused seeds.\n");
    note
}

pub(super) fn failure(path: &str, stage: &str) {
    assert!([
        "resource-run-failed",
        "verifier-failed",
        "discovery-binding-failed"
    ]
    .contains(&stage));
    let destination = Path::new(path).join("EXECUTION_STATUS.json");
    assert!(!destination.exists());
    std::fs::write(
        destination,
        serde_json::to_vec_pretty(
            &json!({"status":"FAILED_OR_UNQUALIFIED","stage":stage,"performance_admitted":false}),
        )
        .unwrap(),
    )
    .unwrap();
}

pub(super) fn seal(path: &str) {
    let root = Path::new(path);
    assert_eq!(
        std::fs::read(root.join("protocol.json")).unwrap(),
        include_bytes!("protocol.json")
    );
    let failed = root.join("EXECUTION_STATUS.json").exists();
    if !failed {
        let result: Value =
            serde_json::from_slice(&std::fs::read(root.join("results.json")).unwrap()).unwrap();
        assert_eq!(field(&result, "qualified").as_bool(), Some(true));
        assert_eq!(field(&result, "all_complete").as_bool(), Some(true));
        if string(&result, "phase") == "full" {
            assert!(root.join("binding.json").exists());
        }
    }
    let mut files = BTreeMap::<String, String>::new();
    for entry in std::fs::read_dir(root).unwrap() {
        let entry = entry.unwrap();
        if !entry.file_type().unwrap().is_file() {
            continue;
        }
        let name = entry.file_name().into_string().unwrap();
        assert_ne!(name, "manifest.json", "fresh output directory required");
        files.insert(name, sha256(&std::fs::read(entry.path()).unwrap()));
    }
    for name in [
        "protocol.json",
        "worker.rs",
        "baseline_worker.rs",
        "verify.rs",
        "Cargo.toml",
        "Cargo.lock",
        "run.sh",
        "boolean-envelope-indexed",
        "raw.jsonl",
        "worker.stderr",
        "git-head.txt",
        "rustc.txt",
        "cpuinfo.txt",
        "meminfo.txt",
        "baseline",
        "baseline-test.stdout",
        "test.stdout",
    ] {
        assert!(files.contains_key(name), "missing retained member {name}");
    }
    if !failed {
        for name in ["conditions.jsonl", "results.json", "RESULT.md"] {
            assert!(files.contains_key(name), "missing qualified member {name}");
        }
    }
    std::fs::write(root.join("manifest.json"), serde_json::to_vec_pretty(&json!({"schema_version":1,"status":if failed {"failed"} else {"qualified"},"files":files})).unwrap()).unwrap();
}

pub(super) fn check_discovery(path: &str, source: &str, binding_path: &str) {
    let root = Path::new(path);
    let manifest_bytes = std::fs::read(root.join("manifest.json")).unwrap();
    let manifest: Value = serde_json::from_slice(&manifest_bytes).unwrap();
    assert_eq!(string(&manifest, "status"), "qualified");
    for (name, digest) in field(&manifest, "files").as_object().unwrap() {
        assert_eq!(
            sha256(&std::fs::read(root.join(name)).unwrap()),
            digest.as_str().unwrap(),
            "changed discovery member {name}"
        );
    }
    assert!(!root.join("EXECUTION_STATUS.json").exists());
    let result: Value =
        serde_json::from_slice(&std::fs::read(root.join("results.json")).unwrap()).unwrap();
    assert_eq!(string(&result, "phase"), "discovery");
    assert_eq!(field(&result, "qualified").as_bool(), Some(true));
    assert_eq!(field(&result, "all_complete").as_bool(), Some(true));
    for name in [
        "worker.rs",
        "verify.rs",
        "baseline_worker.rs",
        "protocol.json",
        "Cargo.toml",
        "Cargo.lock",
        "run.sh",
    ] {
        assert_eq!(
            std::fs::read(root.join(name)).unwrap(),
            std::fs::read(Path::new(source).join(name)).unwrap(),
            "changed discovery source {name}"
        );
    }
    let source_head = std::fs::read_to_string(root.join("git-head.txt")).unwrap();
    std::fs::write(
        binding_path,
        serde_json::to_vec_pretty(&json!({"source_path":path,"manifest_sha256":sha256(&manifest_bytes),"source_head":source_head.trim()})).unwrap(),
    )
    .unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn intervals_are_deterministic_and_use_paired_values() {
        let ratios = [2.5f64; 10];
        assert_eq!(interval(&ratios, 17), [2.5, 2.5]);
        assert_eq!(median(&ratios), 2.5);
    }
    #[test]
    fn replay_rejects_changed_output_or_charged_work() {
        let base = fixture(6, 17);
        let envelope = declared_envelope(&base);
        let expected = oracle(&base, CAPS).unwrap();
        let digest = oracle_digest(std::slice::from_ref(&expected));
        let mut sample = json!({
            "type":"sample","cell":"unit","rep":0,"order":0,"variant":"aa_a",
            "setup_ns":1,"apply_ns":1,"validation_ns":1,"total_ns":3,
            "retained_bytes":Cache::compile(3,&base,&envelope).retained(),
            "output_bytes":expected.payload(),"output_digest":digest,
            "hits":1,"fallbacks":0,"changed_hits":0,"verified_outputs":1
        });
        let expected_rows = [expected.clone()];
        let inputs = [base.clone()];
        assert_eq!(
            verify_sample(
                &sample,
                "unit",
                0,
                0,
                "aa_a",
                3,
                &base,
                &envelope,
                &inputs,
                &expected_rows,
                &expected,
                &digest
            ),
            3
        );
        sample["output_digest"] = json!("wrong");
        assert!(std::panic::catch_unwind(|| verify_sample(
            &sample,
            "unit",
            0,
            0,
            "aa_a",
            3,
            &base,
            &envelope,
            &inputs,
            &expected_rows,
            &expected,
            &digest
        ))
        .is_err());
        sample["output_digest"] = json!(digest);
        sample["total_ns"] = json!(2);
        assert!(std::panic::catch_unwind(|| verify_sample(
            &sample,
            "unit",
            0,
            0,
            "aa_a",
            3,
            &base,
            &envelope,
            &inputs,
            &expected_rows,
            &expected,
            &digest
        ))
        .is_err());
    }
    #[test]
    fn resource_threshold_rejects_a_short_contended_window() {
        let path = std::env::temp_dir().join(format!(
            "boolean-envelope-resource-{}.jsonl",
            std::process::id()
        ));
        let mut sample = json!({
            "schema":"isolated-bench/1","mode":"run","reserved_cpus":[3],
            "host":{"logical_cpus":4,"machine":"aarch64","kernel":"test"},
            "preflight":{"settle":{"seconds":2.0,"other_cpu_seconds":0.0},
                         "conditions":{"psi_cpu":{"some":{"avg10":0.0}},
                                       "psi_memory":{"some":{"avg10":0.0}}}},
            "run":{"command":["/retained/boolean-envelope-indexed","--campaign","discovery","protocol.json"],
                   "cpus":[3],"exit_status":0,"contended":false,
                   "wall_seconds":0.075,"other_cpu_seconds":0.0,
                   "before":{"psi_cpu":{},"psi_memory":{}},
                   "after":{"psi_cpu":{},"psi_memory":{}},
                   "voluntary_switches":0,"involuntary_switches":0,
                   "minor_faults":0,"major_faults":0,"max_rss_kib":1024}
        });
        std::fs::write(&path, format!("{sample}\n")).unwrap();
        let path_text = path.to_str().unwrap();
        verify_resource(path_text, "discovery", "protocol.json", 900.0);
        sample["run"]["other_cpu_seconds"] = json!(0.01);
        sample["run"]["contended"] = json!(true);
        std::fs::write(&path, format!("{sample}\n")).unwrap();
        assert!(std::panic::catch_unwind(|| verify_resource(
            path_text,
            "discovery",
            "protocol.json",
            900.0
        ))
        .is_err());
        std::fs::remove_file(path).unwrap();
    }
}
