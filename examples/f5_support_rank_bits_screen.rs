//! Native one-shot sensitivity screen for exact F5 row-basis table widths.

use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

const PRIMARY: &str = "f5_n24_m24_d4";
const SEEDS: [(&str, u64); 4] = [
    ("frozen", 0),
    ("holdout_a", 0xbadc_0de1),
    ("holdout_b", 0x5eed_2026),
    ("holdout_c", 0xf5c0_2a28),
];
const BITS: [usize; 4] = [8, 5, 6, 7];
const CASES: [&str; 7] = [
    "f5_n12_m12_d4",
    "f5_n16_m16_d3",
    "f5_n16_m16_d4",
    "f5_n20_m20_d3",
    "f5_n20_m20_d4",
    "f5_n24_m24_d3",
    PRIMARY,
];

fn hash_file(path: impl AsRef<Path>) -> String {
    hex::encode(sha256(&fs::read(path).unwrap()))
}

fn save(path: &Path, receipt: &Value) {
    let temp = path.with_extension("json.tmp");
    fs::write(
        &temp,
        format!("{}\n", serde_json::to_string_pretty(receipt).unwrap()),
    )
    .unwrap();
    fs::rename(temp, path).unwrap();
}

fn run_one(binary: &Path, name: &str, seed: u64, bits: usize) -> Value {
    let start = Instant::now();
    let mut child = Command::new(binary)
        .args(["1", "24", "f5", &format!("{seed:016x}")])
        .env("RAYON_NUM_THREADS", "1")
        .env("KIC_GF2_DEFER_ABOVE", "0")
        .env("KIC_GF2_WORD_BATCH", "0")
        .env("KIC_F5_ECHELON", "5")
        .env("KIC_GF2_ROW_BASIS_BITS", bits.to_string())
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
        .env("KIC_GF2_ROW_BASIS_RANK", "1")
        .env("KIC_GF2_ROW_BASIS_TABLES", "8")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
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
    let mut record = json!({
        "workload": name, "seed_xor_hex": format!("{seed:016x}"), "rank_table_bits": bits,
        "process_wall_ms": start.elapsed().as_secs_f64() * 1000.0,
        "exit_code": output.status.code(), "stdout": stdout,
        "stderr": String::from_utf8_lossy(&output.stderr),
        "status": if timed_out { "timeout" } else if output.status.success() { "ok" } else { "failure" },
    });
    if record["status"] != "ok" {
        return record;
    }
    let parsed: Result<Vec<Value>, _> = stdout
        .lines()
        .filter(|line| !line.is_empty())
        .map(serde_json::from_str)
        .collect();
    let Ok(parsed) = parsed else {
        record["status"] = json!("output_error");
        return record;
    };
    let cases: BTreeMap<String, Value> = parsed
        .iter()
        .filter_map(|case| {
            case["case"]
                .as_str()
                .map(|name| (name.to_owned(), case.clone()))
        })
        .collect();
    if parsed.len() != CASES.len()
        || cases.keys().cloned().collect::<BTreeSet<_>>()
            != CASES.into_iter().map(str::to_owned).collect()
    {
        record["status"] = json!("output_error");
    }
    record["cases"] = json!(cases);
    record
}

fn non_timing(case: &Value) -> Value {
    let mut fields = case.as_object().unwrap().clone();
    for key in [
        "wall_ms",
        "criterion_ms",
        "f5_build_ms",
        "reduce_ms",
        "unpack_ms",
    ] {
        fields.remove(key);
    }
    Value::Object(fields)
}

fn compare(reference: &Value, candidate: &Value) -> Vec<String> {
    let mut errors = Vec::new();
    let a = &reference["cases"][PRIMARY];
    let b = &candidate["cases"][PRIMARY];
    let mut old = non_timing(a).as_object().unwrap().clone();
    let mut new = non_timing(b).as_object().unwrap().clone();
    old.remove("reduce_word_ops");
    new.remove("reduce_word_ops");
    if old != new {
        errors.push("primary original-row output or invariant mismatch".to_owned());
    }
    if b["support_split_attempted"] != true
        || b["support_split_original_used"] != true
        || b["support_split_inner_vars"] != 19
    {
        errors.push("fixed cut-19 certificate failed".to_owned());
    }
    for name in CASES.into_iter().filter(|&name| name != PRIMARY) {
        if non_timing(&reference["cases"][name]) != non_timing(&candidate["cases"][name]) {
            errors.push(format!("{name}: inactive output mismatch"));
        }
    }
    errors
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        3,
        "usage: f5_support_rank_bits_screen BINARY RECEIPT"
    );
    let binary = fs::canonicalize(&args[1]).unwrap();
    let output = Path::new(&args[2]);
    fs::create_dir_all(output.parent().unwrap()).unwrap();
    let mut receipt = json!({
        "status":"running", "binary_sha256":hash_file(&binary),
        "driver_source_sha256":hash_file("examples/f5_support_rank_bits_screen.rs"),
        "benchmark_source_sha256":hash_file("examples/f4_f2_bench.rs"),
        "f5_source_sha256":hash_file("src/cryptanalysis/matrix_f5_f2.rs"),
        "gf2_source_sha256":hash_file("src/cryptanalysis/gf2_elim.rs"),
        "rank_table_bits":BITS, "seeds":SEEDS.iter().map(|(name, seed)| (name.to_string(),format!("{seed:016x}"))).collect::<BTreeMap<_,_>>(),
        "host":{"os":std::env::consts::OS,"arch":std::env::consts::ARCH,
            "rustc":String::from_utf8_lossy(&Command::new("rustc").arg("--version").output().unwrap().stdout).trim()},
        "runs":[],
    });
    save(output, &receipt);
    let mut all_ok = true;
    for &(name, seed) in &SEEDS {
        let mut reference = None;
        for bits in BITS {
            let mut run = run_one(&binary, name, seed, bits);
            if run["status"] == "ok" {
                if bits == 8 {
                    if run["cases"][PRIMARY]["support_split_original_used"] != true
                        || run["cases"][PRIMARY]["support_split_inner_vars"] != 19
                    {
                        run["status"] = json!("reference_certificate_failed");
                    } else {
                        reference = Some(run.clone());
                    }
                } else if let Some(reference) = &reference {
                    let errors = compare(reference, &run);
                    if !errors.is_empty() {
                        run["status"] = json!("structural_mismatch");
                        run["errors"] = json!(errors);
                    }
                } else {
                    run["status"] = json!("reference_missing");
                }
            }
            if run["status"] != "ok" {
                all_ok = false;
            }
            receipt["runs"].as_array_mut().unwrap().push(run);
            save(output, &receipt);
        }
    }
    if !all_ok {
        receipt["status"] = json!("failed");
        save(output, &receipt);
        std::process::exit(1);
    }
    let mut decisions = Vec::new();
    for bits in BITS.into_iter().filter(|&bits| bits != 8) {
        let mut cells = Vec::new();
        let mut advance = true;
        for &(name, _) in &SEEDS {
            let runs = receipt["runs"].as_array().unwrap();
            let reference = runs
                .iter()
                .find(|r| r["workload"] == name && r["rank_table_bits"] == 8)
                .unwrap();
            let candidate = runs
                .iter()
                .find(|r| r["workload"] == name && r["rank_table_bits"] == bits)
                .unwrap();
            let ref_ops = reference["cases"][PRIMARY]["reduce_word_ops"]
                .as_u64()
                .unwrap();
            let ops = candidate["cases"][PRIMARY]["reduce_word_ops"]
                .as_u64()
                .unwrap();
            let used = candidate["cases"][PRIMARY]["support_split_original_used"] == true;
            let ratio = ops as f64 / ref_ops as f64;
            if candidate["status"] != "ok" || !used || ratio > 0.90 {
                advance = false;
            }
            cells.push(json!({"workload":name,"used_original_rows":used,
                "inner_rows":candidate["cases"][PRIMARY]["support_split_inner_rows"],
                "outer_rows":candidate["cases"][PRIMARY]["support_split_outer_rows"],
                "reduce_word_ops":ops,"reference_word_ops":ref_ops,"work_ratio":ratio}));
        }
        decisions.push(json!({"rank_table_bits":bits,"advance":advance,"cells":cells}));
    }
    receipt["decisions"] = json!(decisions);
    receipt["status"] = json!(if all_ok { "complete" } else { "failed" });
    save(output, &receipt);
    println!(
        "{}",
        serde_json::to_string_pretty(&receipt["decisions"]).unwrap()
    );
    if !all_ok {
        std::process::exit(1);
    }
}
