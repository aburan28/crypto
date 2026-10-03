//! Native paired complete-call benchmark for the fixed cut-19 F5 form.
//! Run under tools/isolated_bench.py reserve on Linux for qualified timings.

use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

const PRIMARY: &str = "f5_n24_m24_d4";
const CASES: [&str; 7] = [
    "f5_n12_m12_d4",
    "f5_n16_m16_d3",
    "f5_n16_m16_d4",
    "f5_n20_m20_d3",
    "f5_n20_m20_d4",
    "f5_n24_m24_d3",
    PRIMARY,
];
const TIMING: [&str; 5] = [
    "wall_ms",
    "criterion_ms",
    "f5_build_ms",
    "reduce_ms",
    "unpack_ms",
];

fn hash_file(path: impl AsRef<Path>) -> String {
    hex::encode(sha256(&fs::read(path).expect("hash input")))
}

fn write_receipt(path: &Path, value: &Value) {
    let temp = path.with_extension("json.tmp");
    fs::write(
        &temp,
        format!("{}\n", serde_json::to_string_pretty(value).unwrap()),
    )
    .unwrap();
    fs::rename(temp, path).unwrap();
}

fn command_output(program: &str, args: &[&str]) -> String {
    let output = Command::new(program)
        .args(args)
        .output()
        .expect("host command");
    assert!(output.status.success(), "{program} failed");
    String::from_utf8_lossy(&output.stdout).trim().to_owned()
}

fn seed_for(workload: &str) -> u64 {
    match workload {
        "frozen" => 0,
        "holdout_a" => 0xbadc_0de1,
        "holdout_b" => 0x5eed_2026,
        "holdout_c" => 0xf5c0_2a28,
        "confirm_a" => 0xc0d4_1ab3_8b87_2c0d,
        "confirm_b" => 0xae05_fd63_76bf_7a39,
        "confirm_c" => 0x467a_3257_c4d7_8bba,
        "confirm_d" => 0xff83_e4df_6657_b12f,
        _ => panic!("unknown workload: {workload}"),
    }
}

fn run_one(
    binary: &Path,
    workload: &str,
    seed: u64,
    threads: usize,
    mode: u8,
    candidate_rank_bits: u8,
    phase: &str,
    pair: usize,
    position: usize,
) -> Value {
    let mut child = Command::new(binary)
        .args(["1", "24", "f5", &format!("{seed:016x}")])
        .env("RAYON_NUM_THREADS", threads.to_string())
        .env("KIC_GF2_DEFER_ABOVE", "0")
        .env("KIC_GF2_WORD_BATCH", "0")
        .env("KIC_F5_ECHELON", if mode == 1 { "5" } else { "2" })
        .env("KIC_F5_COLUMN_SLACK", "512")
        .env("KIC_F5_FUSED_BUILD", "1")
        .env("KIC_GF2_FORCE_AVX2", "1")
        .env("KIC_F5_AVX512_UNPACK", "0")
        .env("KIC_F5_DIRECT_PACK", "1")
        .env("KIC_GF2_REUSE_TABLE", "1")
        .env("KIC_GF2_TABLES", "4")
        .env("KIC_GF2_SIMD", "1")
        .env("KIC_F5_UNPACK_DIRECT", "1")
        .env("KIC_GF2_AVX2_TABLE_BUILD", "0")
        .env("KIC_GF2_BRANCHLESS_STRIP", "1")
        .env("KIC_GF2_ROW_BASIS_RANK", if mode == 1 { "1" } else { "0" })
        .env(
            "KIC_GF2_ROW_BASIS_BITS",
            if mode == 1 { candidate_rank_bits } else { 8 }.to_string(),
        )
        .env(
            "KIC_GF2_ROW_BASIS_TABLES",
            if mode == 1 { "8" } else { "4" },
        )
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("benchmark spawn");
    let start = Instant::now();
    let timed_out = loop {
        if child.try_wait().unwrap().is_some() {
            break false;
        }
        if start.elapsed() >= Duration::from_secs(120) {
            child.kill().unwrap();
            break true;
        }
        thread::sleep(Duration::from_millis(20));
    };
    let output = child.wait_with_output().unwrap();
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    let mut record = json!({
        "phase": phase, "workload": workload, "seed_xor_hex": format!("{seed:016x}"),
        "pair": pair, "position": position, "mode": mode,
        "rank_table_bits_cap": if mode == 1 { candidate_rank_bits } else { 8 },
        "row_basis_rank_tables": if mode == 1 { 8 } else { 4 },
        "rayon_threads": threads, "process_wall_ms": start.elapsed().as_secs_f64() * 1000.0,
        "exit_code": output.status.code(), "stdout": stdout, "stderr": stderr,
        "status": if timed_out { "timeout" } else if output.status.success() { "ok" } else { "failure" },
    });
    if record["status"] != "ok" {
        return record;
    }
    let cases: Result<Vec<Value>, _> = stdout
        .lines()
        .filter(|s| !s.trim().is_empty())
        .map(serde_json::from_str)
        .collect();
    let Ok(cases) = cases else {
        record["status"] = json!("output_error");
        return record;
    };
    let parsed: BTreeMap<String, Value> = cases
        .iter()
        .filter_map(|case| {
            case["case"]
                .as_str()
                .map(|name| (name.to_owned(), case.clone()))
        })
        .collect();
    let expected: BTreeSet<_> = CASES.into_iter().map(str::to_owned).collect();
    if parsed.len() != cases.len() || parsed.keys().cloned().collect::<BTreeSet<_>>() != expected {
        record["status"] = json!("output_error");
    }
    record["cases"] = json!(parsed);
    record
}

fn non_timing(case: &Value) -> Value {
    let mut fields = case.as_object().unwrap().clone();
    for key in TIMING {
        fields.remove(key);
    }
    Value::Object(fields)
}

fn structural_errors(reference: &Value, candidate: &Value) -> Vec<String> {
    let mut errors = Vec::new();
    for name in CASES {
        let a = &reference["cases"][name];
        let b = &candidate["cases"][name];
        for key in [
            "row_space_fp",
            "rank",
            "criterion_word_ops",
            "rows_f4",
            "rows_built",
            "cols",
            "rows_pruned",
            "direct_pack_used",
            "direct_unpack_used",
        ] {
            if a[key] != b[key] {
                errors.push(format!("{name}: {key} mismatch"));
            }
        }
        if a["column_cert_attempted"] != false
            || a["selected_original_used"] != false
            || a["support_split_attempted"] != false
            || a["support_split_inner_vars"] != 0
        {
            errors.push(format!("{name}: reference certificate route active"));
        }
        if name == PRIMARY {
            if b["support_split_attempted"] != true
                || b["support_split_original_used"] != true
                || b["support_split_inner_vars"] != 19
            {
                errors.push(format!("{name}: support route failed"));
            }
            if b["column_cert_attempted"] != false
                || b["selected_original_used"] != false
                || b["selected_cols"] != 0
            {
                errors.push(format!("{name}: selected-column route active"));
            }
            let inner = b["support_split_inner_rows"].as_u64().unwrap_or(0);
            let outer = b["support_split_outer_rows"].as_u64().unwrap_or(0);
            if inner == 0 || outer == 0 || inner + outer > b["rows_built"].as_u64().unwrap_or(0) {
                errors.push(format!("{name}: invalid support row counts"));
            }
            let old_ops = a["reduce_word_ops"].as_u64().unwrap_or(0);
            let new_ops = b["reduce_word_ops"].as_u64().unwrap_or(u64::MAX);
            if new_ops.saturating_mul(100) > old_ops.saturating_mul(35) {
                errors.push(format!("{name}: counted work exceeded 35% gate"));
            }
        } else if non_timing(a) != non_timing(b) {
            errors.push(format!("{name}: inactive output or route mismatch"));
        }
    }
    errors
}

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    sorted[sorted.len() / 2]
}

fn bootstrap_bounds(values: &[f64]) -> [f64; 2] {
    assert_eq!(values.len(), 5);
    let mut sampled = Vec::with_capacity(3125);
    for n in 0..3125 {
        let mut k = n;
        let mut draw = [0.0; 5];
        for v in &mut draw {
            *v = values[k % 5];
            k /= 5;
        }
        sampled.push(median(&draw));
    }
    sampled.sort_by(f64::total_cmp);
    [
        sampled[((sampled.len() - 1) as f64 * 0.025) as usize],
        sampled[((sampled.len() - 1) as f64 * 0.975) as usize],
    ]
}

fn summary(runs: &[Value]) -> Value {
    let mut result = serde_json::Map::new();
    for name in CASES {
        let mut fields = serde_json::Map::new();
        for field in ["wall_ms", "f5_build_ms", "reduce_ms", "unpack_ms"] {
            let mut aa = Vec::new();
            let mut paired = Vec::new();
            for phase in ["aa", "paired"] {
                for pair in 0..5 {
                    let group: Vec<_> = runs
                        .iter()
                        .filter(|r| r["phase"] == phase && r["pair"] == pair)
                        .collect();
                    assert_eq!(group.len(), 2);
                    let (a, b) = if phase == "aa" {
                        (&group[0], &group[1])
                    } else if group[0]["mode"] == 0 {
                        (&group[0], &group[1])
                    } else {
                        (&group[1], &group[0])
                    };
                    let ratio = a["cases"][name][field].as_f64().unwrap()
                        / b["cases"][name][field].as_f64().unwrap();
                    if phase == "aa" {
                        aa.push(ratio);
                    } else {
                        paired.push(ratio);
                    }
                }
            }
            fields.insert(field.to_owned(), json!({"aa_ratios": aa, "aa_min": aa.iter().copied().fold(f64::INFINITY, f64::min),
                "aa_max": aa.iter().copied().fold(0.0, f64::max), "prior_new": {
                    "ratios": paired, "median": median(&paired), "bootstrap_95pct": bootstrap_bounds(&paired)
                }}));
        }
        result.insert(name.to_owned(), Value::Object(fields));
    }
    Value::Object(result)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        args.len() == 5 || args.len() == 6,
        "usage: f5_support_cut19_pair BINARY OUTPUT WORKLOAD THREADS [CANDIDATE_RANK_BITS]"
    );
    let binary = fs::canonicalize(&args[1]).unwrap();
    let output = PathBuf::from(&args[2]);
    let workload = &args[3];
    let threads: usize = args[4].parse().unwrap();
    assert!(threads == 1 || threads == 2);
    let candidate_rank_bits: u8 = args.get(5).map_or(8, |value| value.parse().unwrap());
    assert!(candidate_rank_bits == 6 || candidate_rank_bits == 8);
    let seed = seed_for(workload);
    fs::create_dir_all(output.parent().unwrap()).unwrap();
    let mut receipt = json!({
        "status": "running", "primary_case": PRIMARY, "workload": workload,
        "seed_xor_hex": format!("{seed:016x}"), "pairs_per_phase": 5,
        "rayon_threads": threads, "candidate_rank_bits_cap": candidate_rank_bits,
        "binary_sha256": hash_file(&binary),
        "benchmark_source_sha256": hash_file("examples/f4_f2_bench.rs"),
        "harness_source_sha256": hash_file("examples/f5_support_cut19_pair.rs"),
        "f5_source_sha256": hash_file("src/cryptanalysis/matrix_f5_f2.rs"),
        "row_builder_source_sha256": hash_file("src/cryptanalysis/koblitz_groebner.rs"),
        "gf2_source_sha256": hash_file("src/cryptanalysis/gf2_elim.rs"),
        "mono_source_sha256": hash_file("src/cryptanalysis/pq_groebner_f2.rs"),
        "git_sha": command_output("git", &["rev-parse", "HEAD"]),
        "head_ref": std::env::var("GITHUB_HEAD_REF").ok(),
        "host": {"os": std::env::consts::OS, "arch": std::env::consts::ARCH,
            "logical_cpus": thread::available_parallelism().unwrap().get(),
            "rustc": command_output("rustc", &["--version"]),
            "cpuinfo": fs::read_to_string("/proc/cpuinfo").ok().map(|s| s.lines().take(30).collect::<Vec<_>>().join("\n")),
            "affinity": command_output("sh", &["-c", "command -v taskset >/dev/null && taskset -pc $$ || true"]),
        },
        "qualified_isolation": false,
        "isolation_reason": "Use external isolated_bench reserve receipt; local Apple runs are nonpromoting",
        "runs": [],
    });
    write_receipt(&output, &receipt);
    let mut reference: Option<Value> = None;
    let mut signatures: BTreeMap<u8, BTreeMap<String, Value>> = BTreeMap::new();
    let mut sequence = vec![("warmup", 0, 0, 0), ("warmup", 0, 1, 1)];
    for pair in 0..5 {
        sequence.extend([("aa", pair, 0, 0), ("aa", pair, 1, 0)]);
    }
    for pair in 0..5 {
        sequence.extend([
            ("paired", pair, 0, (pair % 2) as u8),
            ("paired", pair, 1, (1 - pair % 2) as u8),
        ]);
    }
    for (phase, pair, position, mode) in sequence {
        let mut record = run_one(
            &binary,
            workload,
            seed,
            threads,
            mode,
            candidate_rank_bits,
            phase,
            pair,
            position,
        );
        if record["status"] == "ok" {
            let signature: BTreeMap<String, Value> = CASES
                .into_iter()
                .map(|name| (name.to_owned(), non_timing(&record["cases"][name])))
                .collect();
            match signatures.get(&mode) {
                Some(previous) if previous != &signature => {
                    record["status"] = json!("output_mismatch")
                }
                None => {
                    signatures.insert(mode, signature);
                }
                _ => {}
            }
            if record["status"] == "ok" {
                if record["cases"][PRIMARY]["direct_pack_used"] != true
                    || record["cases"][PRIMARY]["direct_unpack_used"] != true
                {
                    record["status"] = json!("route_mismatch");
                } else if mode == 0 {
                    if record["cases"][PRIMARY]["support_split_attempted"] != false
                        || record["cases"][PRIMARY]["support_split_inner_vars"] != 0
                    {
                        record["status"] = json!("reference_route_mismatch");
                    } else {
                        reference = Some(record.clone());
                    }
                } else {
                    let errors = structural_errors(reference.as_ref().unwrap(), &record);
                    if !errors.is_empty() {
                        record["status"] = json!("candidate_structural_mismatch");
                        record["errors"] = json!(errors);
                    }
                }
            }
        }
        let failed = record["status"] != "ok";
        receipt["runs"].as_array_mut().unwrap().push(record);
        receipt["output_signatures"] = json!(signatures);
        if failed {
            receipt["status"] =
                receipt["runs"].as_array().unwrap().last().unwrap()["status"].clone();
        }
        write_receipt(&output, &receipt);
        if failed {
            std::process::exit(1);
        }
    }
    receipt["summary"] = summary(receipt["runs"].as_array().unwrap());
    receipt["status"] = json!("complete");
    write_receipt(&output, &receipt);
    println!(
        "{}",
        serde_json::to_string_pretty(&receipt["summary"][PRIMARY]).unwrap()
    );
}
