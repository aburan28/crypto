//! Native custody, execution, and replay for the frozen N53 root-batch panel.
//! Usage: n53_s3_root_batch_run <generate|run|analyze> <ic> <rho> <fixture> <replay>

use crypto_lib::binary_ecc::curve::scalar_mul;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::fs::{self, File};
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

const STUDY: &str = "research/notes/ecc2k130/n53_s3_root_batch_20261012";
const HASH_SEEDS: [u64; 3] = [53281307, 53291307, 53301307];
const RANK_SEEDS: [u64; 3] = [532101, 532102, 532103];
const RHO_SEEDS: [u64; 3] = [531201, 531202, 531203];
const CAP_SECS: u64 = 90;
const RSS_CAP: u64 = 16 * 1024 * 1024 * 1024;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn study() -> PathBuf {
    root().join(STUDY)
}

fn sha(path: &Path) -> Result<String, String> {
    let result = Command::new("/usr/bin/shasum")
        .args(["-a", "256"])
        .arg(path)
        .output()
        .map_err(|e| format!("shasum {}: {e}", path.display()))?;
    if !result.status.success() {
        return Err(format!("shasum {} failed", path.display()));
    }
    String::from_utf8(result.stdout)
        .map_err(|e| e.to_string())?
        .split_whitespace()
        .next()
        .map(str::to_owned)
        .ok_or_else(|| "empty shasum result".into())
}

fn git(args: &[&str]) -> Result<String, String> {
    let result = Command::new("git")
        .args(args)
        .current_dir(root())
        .output()
        .map_err(|e| format!("git {args:?}: {e}"))?;
    if !result.status.success() {
        return Err(format!(
            "git {args:?}: {}",
            String::from_utf8_lossy(&result.stderr)
        ));
    }
    String::from_utf8(result.stdout).map_err(|e| e.to_string())
}

fn require_pushed_clean() -> Result<String, String> {
    let head = git(&["rev-parse", "HEAD"])?;
    let upstream = git(&["rev-parse", "@{u}"])?;
    if head != upstream || !git(&["status", "--porcelain"])?.is_empty() {
        return Err("source and frozen inputs must be committed and pushed, worktree clean".into());
    }
    Ok(head.trim().to_owned())
}

fn write_new(path: &Path, bytes: impl AsRef<[u8]>) -> Result<(), String> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    fs::write(path, bytes).map_err(|e| format!("{}: {e}", path.display()))
}

fn read_json(path: &Path) -> Result<Value, String> {
    serde_json::from_slice(&fs::read(path).map_err(|e| e.to_string())?)
        .map_err(|e| format!("{}: {e}", path.display()))
}

fn one_jsonl(path: &Path) -> Result<Value, String> {
    let data = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let rows: Vec<_> = data.lines().filter(|line| !line.is_empty()).collect();
    if rows.len() != 1 {
        return Err(format!("{}: expected one JSON row", path.display()));
    }
    serde_json::from_str(rows[0]).map_err(|e| e.to_string())
}

fn point(value: &Value) -> Result<BinaryPoint, String> {
    let [x, y]: [u64; 2] = serde_json::from_value(value.clone()).map_err(|e| e.to_string())?;
    Ok(BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(x), 53),
        y: F2mElement::from_biguint(&BigUint::from(y), 53),
    })
}

fn schedule() -> Vec<Value> {
    let arms = ["control", "candidate", "rho"];
    let mut cells = Vec::with_capacity(27);
    for target_index in 0..3 {
        for seed_index in 0..3 {
            for offset in 0..3 {
                cells.push(json!({
                    "target_index":target_index,
                    "seed_index":seed_index,
                    "arm":arms[(seed_index + offset) % 3],
                }));
            }
        }
    }
    cells
}

fn generate(
    fixture_bin: &Path,
    ic_bin: &Path,
    rho_bin: &Path,
    replay_bin: &Path,
) -> Result<(), String> {
    let source_head = require_pushed_clean()?;
    let curve = KoblitzCurve::new(0, 53).ok_or("N53 curve unavailable")?;
    let r = curve
        .subgroup_order
        .to_u64()
        .ok_or("subgroup outside u64")?;
    let mut points = Vec::new();
    let mut fixtures = Vec::new();
    for (target_index, seed) in HASH_SEEDS.iter().enumerate() {
        let args = ["53", "0", "signed_frobenius", "1", "strong", "700003"];
        let result = Command::new(fixture_bin)
            .args(args)
            .arg(format!("hash:{seed}"))
            .current_dir(root())
            .output()
            .map_err(|e| format!("fixture {seed}: {e}"))?;
        let input_dir = study().join("inputs");
        let fixture_path = input_dir.join(format!("fixture_{target_index}.jsonl"));
        let stderr_path = input_dir.join(format!("fixture_{target_index}.stderr"));
        write_new(&fixture_path, &result.stdout)?;
        write_new(&stderr_path, &result.stderr)?;
        if !result.status.success() {
            return Err(format!("fixture {seed} exit {:?}", result.status.code()));
        }
        let fixture = one_jsonl(&fixture_path)?;
        if fixture["n"] != 53
            || fixture["a"] != 0
            || fixture["target_kind"] != "public_hash_to_curve_cofactor"
            || fixture["public_hash_seed"] != *seed
            || fixture["published_fixture_scalar"] != Value::Null
            || fixture["verified"] != true
        {
            return Err(format!("fixture {seed} metadata failed"));
        }
        let q = point(&fixture["published_q"])?;
        let recovered = fixture["recovered_fixture_scalar"]
            .as_u64()
            .ok_or("fixture scalar missing")?;
        if !curve.curve.is_on_curve(&q)
            || scalar_mul(&curve.curve, &q, &BigUint::from(r)) != BinaryPoint::Infinity
            || scalar_mul(&curve.curve, curve.generator(), &BigUint::from(recovered)) != q
            || points.contains(&fixture["published_q"])
        {
            return Err(format!(
                "fixture {seed} general-law or uniqueness replay failed"
            ));
        }
        let compact = fixture["published_q"].to_string();
        let spaced = compact.replace(',', ", ");
        for pattern in [&compact, &spaced] {
            let found = Command::new("git")
                .args([
                    "grep",
                    "-l",
                    "-F",
                    "-e",
                    pattern,
                    "HEAD",
                    "--",
                    "research",
                    "experiments",
                    "docs",
                ])
                .current_dir(root())
                .output()
                .map_err(|e| e.to_string())?;
            if found.status.success() {
                return Err(format!(
                    "fixture {seed} duplicates committed point: {}",
                    String::from_utf8_lossy(&found.stdout)
                ));
            }
            if found.status.code() != Some(1) {
                return Err(format!(
                    "git grep point check failed: {}",
                    String::from_utf8_lossy(&found.stderr)
                ));
            }
        }
        write_new(
            &input_dir.join(format!("point_{target_index}.jsonl")),
            format!("{compact}\n"),
        )?;
        points.push(fixture["published_q"].clone());
        fixtures.push(json!({
            "target_index":target_index, "public_hash_seed":seed,
            "public_q":fixture["published_q"],
            "fixture_sha256":sha(&fixture_path)?,
            "point_sha256":sha(&input_dir.join(format!("point_{target_index}.jsonl")))?,
            "independent_general_law_replay":true
        }));
    }
    let frozen = json!({
        "schema":"n53-s3-root-batch-frozen-v1",
        "source_head":source_head,
        "curve_id":"EC1N53Ce0hb097de99be9a",
        "subgroup_order":r,
        "base_hash":"7af2460c8b5a2c29f9d1aa7fefecbcde3a6ce761293d0dab0cecfc3bc980b973",
        "rank_seeds":RANK_SEEDS, "rho_seeds":RHO_SEEDS,
        "schedule":schedule(),
        "wall_cap_seconds":CAP_SECS, "rss_cap_bytes":RSS_CAP,
        "source_sha256":sha(&root().join("examples/koblitz_s3_root_batch_experiment.rs"))?,
        "runner_sha256":sha(&root().join("examples/n53_s3_root_batch_run.rs"))?,
        "rho_source_sha256":sha(&root().join("examples/koblitz_rho_batch_ks_strong_online.rs"))?,
        "fixture_source_sha256":sha(&root().join("examples/koblitz_rho_fixture.rs"))?,
        "replay_source_sha256":sha(&root().join("examples/koblitz_one_target_replay.rs"))?,
        "ic_binary_sha256":sha(ic_bin)?, "rho_binary_sha256":sha(rho_bin)?,
        "fixture_binary_sha256":sha(fixture_bin)?, "replay_binary_sha256":sha(replay_bin)?,
        "targets":fixtures
    });
    write_new(&study().join("FROZEN.json"), format!("{frozen:#}\n"))?;
    println!("generated 3 disjoint public Q; commit and push inputs/FROZEN.json before run");
    Ok(())
}

fn ensure_frozen(ic_bin: &Path, rho_bin: &Path, replay_bin: &Path) -> Result<Value, String> {
    require_pushed_clean()?;
    let frozen = read_json(&study().join("FROZEN.json"))?;
    for (key, path) in [
        (
            "source_sha256",
            root().join("examples/koblitz_s3_root_batch_experiment.rs"),
        ),
        (
            "runner_sha256",
            root().join("examples/n53_s3_root_batch_run.rs"),
        ),
        (
            "rho_source_sha256",
            root().join("examples/koblitz_rho_batch_ks_strong_online.rs"),
        ),
        (
            "replay_source_sha256",
            root().join("examples/koblitz_one_target_replay.rs"),
        ),
        ("ic_binary_sha256", ic_bin.to_path_buf()),
        ("rho_binary_sha256", rho_bin.to_path_buf()),
        ("replay_binary_sha256", replay_bin.to_path_buf()),
    ] {
        if frozen[key].as_str() != Some(sha(&path)?.as_str()) {
            return Err(format!("frozen {key} differs from {}", path.display()));
        }
    }
    for (target_index, target) in frozen["targets"]
        .as_array()
        .ok_or("targets absent")?
        .iter()
        .enumerate()
    {
        let path = study().join(format!("inputs/point_{target_index}.jsonl"));
        if target["point_sha256"].as_str() != Some(sha(&path)?.as_str()) {
            return Err(format!("frozen target {target_index} differs"));
        }
    }
    Ok(frozen)
}

fn run_cell(
    target_index: usize,
    seed_index: usize,
    arm: &str,
    ic_bin: &Path,
    rho_bin: &Path,
    frozen: &Value,
) -> Result<(), String> {
    let directory = study().join(format!("runs/q{target_index}/s{seed_index}/{arm}"));
    if directory.exists() {
        return Err(format!("refusing to overwrite {}", directory.display()));
    }
    fs::create_dir_all(&directory).map_err(|e| e.to_string())?;
    let q = frozen["targets"][target_index]["public_q"]
        .as_array()
        .ok_or("public point absent")?;
    let point_file = study().join(format!("inputs/point_{target_index}.jsonl"));
    let (binary, args, output_name) = if arm == "rho" {
        (
            rho_bin,
            vec![
                "53".to_owned(),
                "0".to_owned(),
                "signed_frobenius".to_owned(),
                "1".to_owned(),
                RHO_SEEDS[seed_index].to_string(),
            ],
            "rho.jsonl",
        )
    } else {
        (
            ic_bin,
            vec![
                "construct:53:0:220".to_owned(),
                point_file.display().to_string(),
                RANK_SEEDS[seed_index].to_string(),
                directory.join("targets.jsonl").display().to_string(),
            ],
            "summary.jsonl",
        )
    };
    let intent = json!({
        "schema":"n53-s3-root-batch-intent-v1",
        "target_index":target_index, "seed_index":seed_index, "arm":arm,
        "rank_seed":RANK_SEEDS[seed_index], "rho_seed":RHO_SEEDS[seed_index],
        "command":{"binary":binary,"args":args},
        "environment":if arm == "rho" {
            json!({"KIC_RHO_TARGET_POINT":format!("{},{}",q[0],q[1]),
                "KIC_RHO_RUNG":"3","KIC_RHO_LANES":"32","KIC_RHO_DP_BITS":"4"})
        } else {
            json!({"KIC_PAIR_SYMMETRY":"0","KIC_ROOT_BLOOM":"0",
                "KIC_S3_ROOT_BATCH":if arm == "candidate" {"1"} else {"0"},
                "KIC_RANK_PROBE_CAP":Value::Null})
        },
        "point_sha256":frozen["targets"][target_index]["point_sha256"],
        "binary_sha256":sha(binary)?, "frozen_sha256":sha(&study().join("FROZEN.json"))?,
        "wall_cap_seconds":CAP_SECS, "rss_cap_bytes":RSS_CAP
    });
    write_new(&directory.join("intent.json"), format!("{intent:#}\n"))?;
    let stdout = File::create(directory.join(output_name)).map_err(|e| e.to_string())?;
    let stderr = File::create(directory.join("producer.stderr")).map_err(|e| e.to_string())?;
    let mut command = Command::new("/usr/bin/time");
    command
        .arg("-l")
        .arg(binary)
        .args(&args)
        .current_dir(root());
    command
        .stdout(Stdio::from(stdout))
        .stderr(Stdio::from(stderr));
    for key in [
        "KIC_RANK_PROBE_CAP",
        "KIC_RHO_EXPLICIT_SCALAR",
        "KIC_RHO_BATCH_CORPUS",
    ] {
        command.env_remove(key);
    }
    if arm == "rho" {
        command.env("KIC_RHO_TARGET_POINT", format!("{},{}", q[0], q[1]));
        command
            .env("KIC_RHO_RUNG", "3")
            .env("KIC_RHO_LANES", "32")
            .env("KIC_RHO_DP_BITS", "4");
    } else {
        command.env_remove("KIC_RHO_TARGET_POINT");
        command.env("KIC_PAIR_SYMMETRY", "0");
        command.env("KIC_ROOT_BLOOM", "0");
        command.env(
            "KIC_S3_ROOT_BATCH",
            if arm == "candidate" { "1" } else { "0" },
        );
        command.env("KIC_DUMP_BASE", directory.join("base.jsonl"));
        command.env("KIC_DUMP_RANK", directory.join("rank.jsonl"));
    }
    // SAFETY: setpgid runs in the child before exec; the parent may then kill
    // the full time-wrapper/process group when the frozen wall cap fires.
    unsafe {
        command.pre_exec(|| {
            if libc::setpgid(0, 0) == 0 {
                Ok(())
            } else {
                Err(std::io::Error::last_os_error())
            }
        });
    }
    let began = Instant::now();
    let mut child = command.spawn().map_err(|e| e.to_string())?;
    let (exit_code, timed_out) = loop {
        if let Some(status) = child.try_wait().map_err(|e| e.to_string())? {
            break (status.code(), false);
        }
        if began.elapsed() >= Duration::from_secs(CAP_SECS) {
            // SAFETY: child is the leader of its own process group.
            unsafe {
                libc::kill(-(child.id() as i32), libc::SIGKILL);
            }
            let status = child.wait().map_err(|e| e.to_string())?;
            break (status.code(), true);
        }
        thread::sleep(Duration::from_millis(100));
    };
    let wall_ms = began.elapsed().as_secs_f64() * 1000.0;
    let stderr_text =
        fs::read_to_string(directory.join("producer.stderr")).map_err(|e| e.to_string())?;
    let peak_rss_bytes = stderr_text.lines().find_map(|line| {
        line.contains("maximum resident set size")
            .then(|| line.split_whitespace().next()?.parse::<u64>().ok())
            .flatten()
    });
    let status = json!({
        "schema":"n53-s3-root-batch-status-v1", "arm":arm,
        "exit_code":exit_code, "timed_out":timed_out, "wall_ms":wall_ms,
        "peak_rss_bytes":peak_rss_bytes,
        "status":if !timed_out && exit_code == Some(0) && peak_rss_bytes.is_some_and(|rss| rss <= RSS_CAP) { "COMPLETED" } else { "INCOMPLETE" },
        "stdout_sha256":sha(&directory.join(output_name))?,
        "stderr_sha256":sha(&directory.join("producer.stderr"))?
    });
    write_new(&directory.join("status.json"), format!("{status:#}\n"))?;
    println!(
        "q{target_index} s{seed_index} {arm}: {} ({wall_ms:.1} ms)",
        status["status"]
    );
    Ok(())
}

fn run(ic_bin: &Path, rho_bin: &Path, replay_bin: &Path) -> Result<(), String> {
    let frozen = ensure_frozen(ic_bin, rho_bin, replay_bin)?;
    if frozen["schedule"] != json!(schedule()) {
        return Err("frozen schedule differs from native runner".into());
    }
    for cell in frozen["schedule"].as_array().ok_or("schedule absent")? {
        let target_index = cell["target_index"].as_u64().ok_or("target index absent")? as usize;
        let seed_index = cell["seed_index"].as_u64().ok_or("seed index absent")? as usize;
        let arm = cell["arm"].as_str().ok_or("arm absent")?;
        if target_index >= 3 || seed_index >= 3 {
            return Err("schedule index outside frozen grid".into());
        }
        run_cell(target_index, seed_index, arm, ic_bin, rho_bin, &frozen)?;
    }
    Ok(())
}

fn selected_trace(trace: &Path) -> Result<Vec<Value>, String> {
    let text = fs::read_to_string(trace).map_err(|e| e.to_string())?;
    let mut rows = Vec::new();
    for line in text.lines().filter(|line| !line.is_empty()) {
        let mut row: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
        for key in [
            "query_ms",
            "pdp_ms",
            "relation_check_ms",
            "root_calls",
            "root_results_consumed",
            "root_ops",
            "root_scratch_peak_bytes",
            "s3_root_batch",
        ] {
            row.as_object_mut()
                .ok_or("trace row not object")?
                .remove(key);
        }
        rows.push(row);
    }
    Ok(rows)
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(|a, b| a.total_cmp(b));
    values[values.len() / 2]
}

fn target_block_interval(values: &[f64]) -> Option<[f64; 2]> {
    if values.len() != 9 {
        return None;
    }
    let mut state = 0x53_BA_7C_12_2026u64;
    let mut draws = Vec::with_capacity(10_000);
    for _ in 0..10_000 {
        let mut sample = Vec::with_capacity(9);
        for _ in 0..3 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let target = (state % 3) as usize;
            sample.extend_from_slice(&values[target * 3..target * 3 + 3]);
        }
        draws.push(median(sample));
    }
    draws.sort_by(|a, b| a.total_cmp(b));
    Some([draws[250], draws[9750]])
}

fn u64_at(value: &Value, key: &str) -> Result<u64, String> {
    value[key].as_u64().ok_or_else(|| format!("missing {key}"))
}

fn normalized_rho(
    rho_dir: &Path,
    target_index: usize,
    seed_index: usize,
    status: &Value,
    frozen: &Value,
) -> Result<PathBuf, String> {
    let raw = fs::read_to_string(rho_dir.join("rho.jsonl")).map_err(|e| e.to_string())?;
    let rows: Vec<Value> = raw
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| serde_json::from_str(line).map_err(|e| e.to_string()))
        .collect::<Result<_, _>>()?;
    if rows.len() != 2 {
        return Err("rho must have exactly one fixture and one summary".into());
    }
    let item = &rows[0];
    let summary = &rows[1];
    let fixture = one_jsonl(&study().join(format!("inputs/fixture_{target_index}.jsonl")))?;
    if item["kind"] != "rho_ks_batch_fixture"
        || summary["kind"] != "rho_ks_batch_summary"
        || item["n"] != 53
        || item["a"] != 0
        || item["published_q"] != frozen["targets"][target_index]["public_q"]
        || item["published_q"] != fixture["published_q"]
        || item["batch_seed"] != RHO_SEEDS[seed_index]
        || summary["batch_seed"] != RHO_SEEDS[seed_index]
        || item["target_source"] != "public_point"
        || summary["target_source"] != "public_point"
        || item["published_fixture_scalar"] != Value::Null
        || item["verified"] != true
        || summary["all_verified"] != true
        || summary["rung"] != 3
        || summary["lanes"] != 32
        || summary["dp_bits"] != 4
        || summary["jump_count"] != 32
        || summary["fixtures"] != 1
        || item["recovered_fixture_scalar"] != fixture["recovered_fixture_scalar"]
    {
        return Err("raw rho settings, point, or scalar differ from frozen fixture".into());
    }
    let online = item["online_ms"].as_f64().ok_or("rho online missing")?;
    let walk = item["walk_ms"].as_f64().ok_or("rho walk missing")?;
    let collision = item["collision_ms"]
        .as_f64()
        .ok_or("rho collision missing")?;
    let recovery = item["recovery_check_ms"]
        .as_f64()
        .ok_or("rho recovery missing")?;
    let setup = summary["setup_ms"].as_f64().ok_or("rho setup missing")?;
    let excluded = item["target_generation_ms_excluded"]
        .as_f64()
        .ok_or("rho target preparation missing")?;
    if (online - (walk + collision + recovery)).abs() > 1e-6 {
        return Err("raw rho online phases differ".into());
    }
    let normal = json!({
        "n":53,"a":0,"subgroup_order":frozen["subgroup_order"],
        "field_modulus_low_terms":fixture["field_modulus_low_terms"],
        "published_q":item["published_q"],
        "published_fixture_scalar":Value::Null,
        "recovered_fixture_scalar":item["recovered_fixture_scalar"],
        "verified":true,"reference_grade":"strong",
        "quotient_mode":"signed_frobenius",
        "target_kind":"public_hash_to_curve_cofactor",
        "public_hash_seed":HASH_SEEDS[target_index],
        "batch_seed":RHO_SEEDS[seed_index],
        "peak_rss_bytes":status["peak_rss_bytes"],
        "walk_ms":walk,"validation_ms":collision+recovery,
        "setup_ms":setup,"target_generation_ms":excluded,
        "total_ms":setup+online+excluded,
        "raw_rho_sha256":sha(&rho_dir.join("rho.jsonl"))?
    });
    let path = rho_dir.join("rho_normalized.jsonl");
    write_new(&path, format!("{normal}\n"))?;
    Ok(path)
}

fn analyze(ic_bin: &Path, rho_bin: &Path, replay_bin: &Path) -> Result<(), String> {
    let frozen = read_json(&study().join("FROZEN.json"))?;
    for (key, path) in [
        (
            "source_sha256",
            root().join("examples/koblitz_s3_root_batch_experiment.rs"),
        ),
        (
            "runner_sha256",
            root().join("examples/n53_s3_root_batch_run.rs"),
        ),
        ("ic_binary_sha256", ic_bin.to_path_buf()),
        ("rho_binary_sha256", rho_bin.to_path_buf()),
        ("replay_binary_sha256", replay_bin.to_path_buf()),
    ] {
        if frozen[key].as_str() != Some(sha(&path)?.as_str()) {
            return Err(format!("analysis source/binary mismatch at {key}"));
        }
    }
    let mut pairs = Vec::new();
    let mut missing = Vec::new();
    let mut cold_ratios = Vec::new();
    let mut online_ratios = Vec::new();
    for target_index in 0..3 {
        for seed_index in 0..3 {
            let dir = study().join(format!("runs/q{target_index}/s{seed_index}"));
            let control = dir.join("control");
            let candidate = dir.join("candidate");
            let rho = dir.join("rho");
            let statuses = [
                control.join("status.json"),
                candidate.join("status.json"),
                rho.join("status.json"),
            ];
            if statuses.iter().any(|path| !path.exists()) {
                missing.push(json!({"target_index":target_index,"seed_index":seed_index,"reason":"missing cell"}));
                continue;
            }
            let status_values: Vec<Value> = statuses
                .iter()
                .map(|p| read_json(p))
                .collect::<Result<_, _>>()?;
            if status_values
                .iter()
                .any(|status| status["status"] != "COMPLETED")
            {
                missing.push(json!({"target_index":target_index,"seed_index":seed_index,"reason":"incomplete cell","statuses":status_values}));
                continue;
            }
            if selected_trace(&control.join("rank.jsonl"))?
                != selected_trace(&candidate.join("rank.jsonl"))?
            {
                missing.push(json!({"target_index":target_index,"seed_index":seed_index,"reason":"rank trace differs"}));
                continue;
            }
            let control_summary = one_jsonl(&control.join("summary.jsonl"))?;
            let candidate_summary = one_jsonl(&candidate.join("summary.jsonl"))?;
            let control_target = one_jsonl(&control.join("targets.jsonl"))?;
            let candidate_target = one_jsonl(&candidate.join("targets.jsonl"))?;
            if control_summary["base_hash"] != frozen["base_hash"]
                || candidate_summary["base_hash"] != frozen["base_hash"]
                || control_summary["root_table_entries"] != candidate_summary["root_table_entries"]
                || sha(&control.join("base.jsonl"))? != sha(&candidate.join("base.jsonl"))?
                || control_summary["rank"] != 220
                || candidate_summary["rank"] != 220
                || control_summary["rank_relations"] != 220
                || candidate_summary["rank_relations"] != 220
                || control_summary["rank_probes_total"] != candidate_summary["rank_probes_total"]
                || control_target["recovered_scalar"] != candidate_target["recovered_scalar"]
                || control_target["point_indices"] != candidate_target["point_indices"]
            {
                missing.push(json!({"target_index":target_index,"seed_index":seed_index,"reason":"base/rank/target mismatch"}));
                continue;
            }
            let mut replays = Vec::new();
            let normalized_rho = match normalized_rho(
                &rho,
                target_index,
                seed_index,
                &status_values[2],
                &frozen,
            ) {
                Ok(path) => path,
                Err(error) => {
                    missing.push(json!({"target_index":target_index,"seed_index":seed_index,"reason":"rho normalization failed","error":error}));
                    continue;
                }
            };
            for (arm, arm_dir) in [("control", &control), ("candidate", &candidate)] {
                let result = Command::new(replay_bin)
                    .args([
                        arm_dir.join("base.jsonl").display().to_string(),
                        arm_dir.join("rank.jsonl").display().to_string(),
                        arm_dir.join("targets.jsonl").display().to_string(),
                        arm_dir.join("summary.jsonl").display().to_string(),
                        normalized_rho.display().to_string(),
                        study()
                            .join(format!("inputs/point_{target_index}.jsonl"))
                            .display()
                            .to_string(),
                        HASH_SEEDS[target_index].to_string(),
                        RHO_SEEDS[seed_index].to_string(),
                    ])
                    .current_dir(root())
                    .output()
                    .map_err(|e| e.to_string())?;
                write_new(&arm_dir.join("replay.stdout"), &result.stdout)?;
                write_new(&arm_dir.join("replay.stderr"), &result.stderr)?;
                if !result.status.success() {
                    missing.push(json!({"target_index":target_index,"seed_index":seed_index,"arm":arm,"reason":"independent replay failed"}));
                    break;
                }
                replays.push(
                    serde_json::from_slice::<Value>(&result.stdout).map_err(|e| e.to_string())?,
                );
            }
            if replays.len() != 2 {
                continue;
            }
            let control_cold = control_summary["cold_in_process_ms"]
                .as_f64()
                .ok_or("control cold missing")?;
            let candidate_cold = candidate_summary["cold_in_process_ms"]
                .as_f64()
                .ok_or("candidate cold missing")?;
            let control_online = control_target["online_ms"]
                .as_f64()
                .ok_or("control online missing")?;
            let candidate_online = candidate_target["online_ms"]
                .as_f64()
                .ok_or("candidate online missing")?;
            let cold_ratio = candidate_cold / control_cold;
            let online_ratio = candidate_online / control_online;
            let control_ops = &control_summary["rank_root_ops"];
            let candidate_ops = &candidate_summary["rank_root_ops"];
            let control_inversions = u64_at(control_ops, "inversion_primitives")?;
            let candidate_inversions = u64_at(candidate_ops, "inversion_primitives")?;
            let count_consistent = [
                (&control_summary, control_ops),
                (&candidate_summary, candidate_ops),
            ]
            .iter()
            .all(|(summary, ops)| {
                let calls = summary["rank_root_calls_total"].as_u64();
                let regular = ops["regular_calls"].as_u64();
                let degenerate = ops["degenerate_calls"].as_u64();
                matches!((calls, regular, degenerate), (Some(c), Some(r), Some(d)) if c == r + d)
            });
            let consumed_same = control_summary["rank_root_results_consumed"]
                == candidate_summary["rank_root_results_consumed"];
            let inversion_gate = candidate_inversions * 2 <= control_inversions;
            let rho_cold = replays[0]["rho_cold_in_process_ms"]
                .as_f64()
                .ok_or("rho cold missing")?;
            let rho_online = replays[0]["rho_online_ms"]
                .as_f64()
                .ok_or("rho online missing")?;
            cold_ratios.push(cold_ratio);
            online_ratios.push(online_ratio);
            pairs.push(json!({
                "target_index":target_index,"seed_index":seed_index,
                "control_cold_ms":control_cold,"candidate_cold_ms":candidate_cold,
                "control_online_ms":control_online,"candidate_online_ms":candidate_online,
                "candidate_over_control_cold":cold_ratio,
                "candidate_over_control_online":online_ratio,
                "control_over_rho_cold":control_cold/rho_cold,
                "candidate_over_rho_cold":candidate_cold/rho_cold,
                "control_over_rho_online":control_online/rho_online,
                "candidate_over_rho_online":candidate_online/rho_online,
                "control_root_ops":control_ops,
                "candidate_root_ops":candidate_ops,
                "control_target_root_ops":control_target["root_ops"],
                "candidate_target_root_ops":candidate_target["root_ops"],
                "operation_counts_consistent":count_consistent,
                "rank_root_results_consumed_equal":consumed_same,
                "inversion_gate":inversion_gate,
                "control_rank_root_calls":control_summary["rank_root_calls_total"],
                "candidate_rank_root_calls":candidate_summary["rank_root_calls_total"],
                "control_rank_probes":control_summary["rank_probes_total"],
                "candidate_rank_probes":candidate_summary["rank_probes_total"],
                "control_replay":replays[0], "candidate_replay":replays[1],
                "rho_online_ms":rho_online,
                "rho_cold_ms":rho_cold,
                "status":"VERIFIED"
            }));
        }
    }
    let analysis = json!({
        "schema":"n53-s3-root-batch-analysis-v1",
        "frozen_sha256":sha(&study().join("FROZEN.json"))?,
        "ic_binary_sha256":sha(ic_bin)?,"rho_binary_sha256":sha(rho_bin)?,
        "replay_binary_sha256":sha(replay_bin)?,
        "requested_pairs":9,"verified_pairs":pairs.len(),"missing_or_invalid":missing,
        "pairs":pairs,
        "candidate_over_control_cold_median":if cold_ratios.is_empty() { Value::Null } else { json!(median(cold_ratios.clone())) },
        "candidate_over_control_online_median":if online_ratios.is_empty() { Value::Null } else { json!(median(online_ratios.clone())) },
        "candidate_over_control_cold_target_block_interval":target_block_interval(&cold_ratios),
        "candidate_over_control_online_target_block_interval":target_block_interval(&online_ratios),
        "exact_gate":pairs.len()==9 && pairs.iter().all(|pair| pair["operation_counts_consistent"]==true && pair["rank_root_results_consumed_equal"]==true),
        "inversion_gate":pairs.len()==9 && pairs.iter().all(|pair| pair["inversion_gate"]==true),
        "descriptive_progression_gate":pairs.len()==9 && median(cold_ratios.clone()) <= 0.90 && median(online_ratios.clone()) <= 1.10,
        "host_timing_class":"macOS_shared_L0_descriptive"
    });
    write_new(&study().join("ANALYSIS.json"), format!("{analysis:#}\n"))?;
    println!(
        "{}",
        json!({"verified_pairs":analysis["verified_pairs"],"missing_or_invalid":analysis["missing_or_invalid"],"cold_median":analysis["candidate_over_control_cold_median"],"online_median":analysis["candidate_over_control_online_median"]})
    );
    Ok(())
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 6 {
        eprintln!(
            "usage: n53_s3_root_batch_run <generate|run|analyze> <ic> <rho> <fixture> <replay>"
        );
        std::process::exit(2);
    }
    let ic = Path::new(&args[2]);
    let rho = Path::new(&args[3]);
    let fixture = Path::new(&args[4]);
    let replay = Path::new(&args[5]);
    let result = match args[1].as_str() {
        "generate" => generate(fixture, ic, rho, replay),
        "run" => run(ic, rho, replay),
        "analyze" => analyze(ic, rho, replay),
        _ => Err("mode must be generate, run, or analyze".into()),
    };
    if let Err(error) = result {
        eprintln!("n53 S3 root-batch {0}: {error}", args[1]);
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_schedule_has_all_cells_and_balanced_positions() {
        let cells = schedule();
        assert_eq!(cells.len(), 27);
        for target in 0..3 {
            for seed in 0..3 {
                let group = &cells[(target * 3 + seed) * 3..(target * 3 + seed + 1) * 3];
                let mut arms: Vec<_> = group
                    .iter()
                    .map(|cell| cell["arm"].as_str().unwrap())
                    .collect();
                arms.sort_unstable();
                assert_eq!(arms, ["candidate", "control", "rho"]);
            }
            for position in 0..3 {
                let mut arms: Vec<_> = (0..3)
                    .map(|seed| {
                        cells[(target * 3 + seed) * 3 + position]["arm"]
                            .as_str()
                            .unwrap()
                    })
                    .collect();
                arms.sort_unstable();
                assert_eq!(arms, ["candidate", "control", "rho"]);
            }
        }
    }

    #[test]
    fn target_block_interval_keeps_constant_ratio() {
        assert_eq!(target_block_interval(&[0.75; 9]), Some([0.75, 0.75]));
    }
}
