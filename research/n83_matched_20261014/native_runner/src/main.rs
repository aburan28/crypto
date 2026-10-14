//! Write-once, predeclared matched IC/rho campaigns for public synthetic points.
mod monitor;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const N: usize = 32;
const MEMORY: u64 = 12 * 1024 * 1024 * 1024;
const BUDGET: u64 = 2100;
fn err(s: impl Into<String>) -> Box<dyn std::error::Error> {
    s.into().into()
}
fn epoch_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_millis()
}
fn read_json(p: &Path) -> Result<Value> {
    Ok(serde_json::from_reader(File::open(p)?)?)
}
fn read_jsonl(p: &Path) -> Result<Vec<Value>> {
    BufReader::new(File::open(p)?)
        .lines()
        .filter_map(|l| match l {
            Ok(l) if l.trim().is_empty() => None,
            other => Some(other),
        })
        .enumerate()
        .map(|(i, l)| {
            Ok(serde_json::from_str(&l?)
                .map_err(|e| err(format!("{}:{}:{e}", p.display(), i + 1)))?)
        })
        .collect()
}
fn hash(p: &Path) -> Result<String> {
    let o = Command::new("/usr/bin/sha256sum")
        .arg("--")
        .arg(p)
        .output()?;
    if !o.status.success() {
        return Err(err(format!("hash failed {}", p.display())));
    }
    Ok(String::from_utf8(o.stdout)?
        .split_whitespace()
        .next()
        .ok_or_else(|| err("empty SHA256"))?
        .into())
}
fn hash_bytes(b: &[u8]) -> Result<String> {
    let mut c = Command::new("/usr/bin/sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    c.stdin
        .take()
        .ok_or_else(|| err("hash stdin missing"))?
        .write_all(b)?;
    let o = c.wait_with_output()?;
    if !o.status.success() {
        return Err(err("source hash failed"));
    }
    Ok(String::from_utf8(o.stdout)?
        .split_whitespace()
        .next()
        .ok_or_else(|| err("empty SHA256"))?
        .into())
}
fn new_json(p: &Path, v: &Value) -> Result<()> {
    let mut f = OpenOptions::new().write(true).create_new(true).open(p)?;
    serde_json::to_writer_pretty(&mut f, v)?;
    f.write_all(b"\n")?;
    f.sync_all()?;
    fs::set_permissions(p, fs::Permissions::from_mode(0o444))?;
    Ok(())
}
fn new_jsonl(p: &Path, rows: &[Value]) -> Result<()> {
    let mut f = OpenOptions::new().write(true).create_new(true).open(p)?;
    for r in rows {
        serde_json::to_writer(&mut f, r)?;
        f.write_all(b"\n")?;
    }
    f.sync_all()?;
    fs::set_permissions(p, fs::Permissions::from_mode(0o444))?;
    Ok(())
}
fn path(v: &Value, k: &str) -> Result<PathBuf> {
    v[k].as_str()
        .map(PathBuf::from)
        .ok_or_else(|| err(format!("protocol {k} missing")))
}
fn scalar(v: &Value) -> Option<String> {
    match v {
        Value::String(s) => Some(s.clone()),
        Value::Number(n) => Some(n.to_string()),
        _ => None,
    }
}
fn point(v: &Value) -> Option<Vec<String>> {
    let a = v.as_array()?;
    if a.len() != 2 {
        return None;
    }
    a.iter().map(scalar).collect()
}
fn policy() -> Value {
    json!({"backend":"fold","scalar_kernel":"reference","lanes":256,"dp_bits":6,"escape_fruitless":false,"step_cap_factor":2000})
}
fn check_hash(v: &Value, key: &str, file: &Path) -> Result<()> {
    let expected = v[key]
        .as_str()
        .ok_or_else(|| err(format!("required digest {key} missing")))?;
    if hash(file)? != expected {
        return Err(err(format!("{key} mismatch before child execution")));
    }
    Ok(())
}
#[derive(Clone)]
struct Fixture {
    id: String,
    point: Value,
    scalar: String,
}
struct Protocol {
    raw: Value,
    digest: String,
    cwd: PathBuf,
    rho: PathBuf,
    ic: PathBuf,
    base: PathBuf,
    base_hash: String,
    fixtures: Vec<Fixture>,
    bindings: BTreeMap<String, String>,
}
impl Protocol {
    fn load(file: &Path) -> Result<Self> {
        let v = read_json(file)?;
        if v["schema"] != "n83.matched-campaign/v1" || v["protocol_id"].as_str().is_none() {
            return Err(err("matched protocol identity missing"));
        }
        if v["resource"]["cpu"] != 3
            || v["resource"]["memory_bytes"] != MEMORY
            || v["resource"]["timeout_s"] != BUDGET
        {
            return Err(err("common envelope must be CPU3/one thread/12GiB/2100s"));
        }
        let cwd = path(&v, "cwd")?;
        let rho = path(&v, "rho_binary")?;
        let ic = path(&v, "ic_binary")?;
        let base = path(&v, "base_file")?;
        let sidecar = path(&v, "validation_sidecar")?;
        let generator = path(&v, "generator_receipt")?;
        let human = path(&v, "human_protocol")?;
        let mut bindings = BTreeMap::new();
        for (key, p) in [
            ("expected_rho_binary_sha256", &rho),
            ("expected_ic_binary_sha256", &ic),
            ("expected_base_file_sha256", &base),
            ("expected_validation_sidecar_sha256", &sidecar),
            ("expected_generator_receipt_sha256", &generator),
            ("human_protocol_sha256", &human),
        ] {
            check_hash(&v, key, p)?;
            bindings.insert(p.to_string_lossy().into(), hash(p)?);
        }
        if let Some(sources) = v["relevant_source_sha256"].as_object() {
            for (name, digest) in sources {
                let candidate = Path::new(name);
                let source = if candidate.is_absolute() {
                    candidate.to_path_buf()
                } else {
                    cwd.join(candidate)
                };
                if hash(&source)?.as_str()
                    != digest
                        .as_str()
                        .ok_or_else(|| err("invalid relevant-source hash"))?
                {
                    return Err(err(format!("relevant source mismatch: {name}")));
                }
                bindings.insert(source.to_string_lossy().into(), hash(&source)?);
            }
        }
        if v["build_receipt"].is_string() {
            let receipt = path(&v, "build_receipt")?;
            check_hash(&v, "build_receipt_sha256", &receipt)?;
            bindings.insert(receipt.to_string_lossy().into(), hash(&receipt)?);
        }
        let g = read_json(&generator)?;
        if g["status"] != "PASS" || g["target_count"] != N {
            return Err(err("fixture generation not PASS for32 targets"));
        }
        let s = read_json(&sidecar)?;
        if s["must_not_be_passed_to_measured_children"] != true {
            return Err(err("sidecar must be validation only"));
        }
        let entries = s["entries"]
            .as_array()
            .ok_or_else(|| err("sidecar entries missing"))?;
        if entries.len() != N {
            return Err(err("exactly32 sidecar entries required"));
        }
        let dir = path(&v, "fixtures_dir")?;
        let mut fixtures = Vec::new();
        let mut unique = BTreeSet::new();
        for i in 1..=N {
            let id = format!("target_{i:02}");
            let f = dir.join(format!("{id}.jsonl"));
            let rows = read_jsonl(&f)?;
            if rows.len() != 1
                || point(&rows[0]).is_none()
                || !unique.insert(point(&rows[0]).unwrap())
            {
                return Err(err(format!(
                    "{id}: exactly one distinct public point required"
                )));
            }
            let es: Vec<_> = entries.iter().filter(|e| e["fixture_id"] == id).collect();
            if es.len() != 1 || es[0]["point"] != rows[0] {
                return Err(err(format!("{id}: sidecar identity mismatch")));
            }
            let expected = v["expected_fixture_file_sha256"][&id]
                .as_str()
                .ok_or_else(|| err(format!("{id}: digest missing")))?;
            if hash(&f)? != expected {
                return Err(err(format!("{id}: fixture digest mismatch")));
            }
            let k = scalar(&es[0]["scalar_validation_only"])
                .ok_or_else(|| err("sidecar scalar missing"))?;
            let parsed: u64 = k.parse()?;
            if parsed == 0 || parsed >= 8_569_786_107_849_059 {
                return Err(err("sidecar scalar outside subgroup"));
            }
            bindings.insert(f.to_string_lossy().into(), expected.into());
            fixtures.push(Fixture {
                id,
                point: rows[0].clone(),
                scalar: k,
            });
        }
        let base_hash = v["expected_base_hash"]
            .as_str()
            .ok_or_else(|| err("expected base identity missing"))?
            .into();
        Ok(Self {
            raw: v,
            digest: hash(file)?,
            cwd,
            rho,
            ic,
            base,
            base_hash,
            fixtures,
            bindings,
        })
    }
    fn fixture(&self, id: &str) -> Result<&Fixture> {
        self.fixtures
            .iter()
            .find(|f| f.id == id)
            .ok_or_else(|| err(format!("unknown target {id}")))
    }
    fn provenance(&self) -> Result<Value> {
        let exe = std::env::current_exe()?;
        let cpu = fs::read_to_string("/proc/cpuinfo")?;
        let u = Command::new("/usr/bin/uname").arg("-a").output()?;
        Ok(
            json!({"controller_executable":exe,"controller_executable_sha256":hash(&exe)?,"compiled_source_sha256":{"main":hash_bytes(include_bytes!("main.rs"))?,"monitor":hash_bytes(include_bytes!("monitor.rs"))?,"manifest":hash_bytes(include_bytes!("../Cargo.toml"))?,"lock":hash_bytes(include_bytes!("../Cargo.lock"))?},"compiled_source_revision":option_env!("N83_SOURCE_COMMIT"),"input_binary_bindings":self.bindings,"fresh_host":{"epoch_ms":epoch_ms(),"uname":String::from_utf8_lossy(&u.stdout).trim(),"cpu_model":cpu.lines().find(|l|l.starts_with("model name")),"cpuinfo_sha256":hash_bytes(cpu.as_bytes())?},"scalar_sidecar_opened_by_parent_only":true}),
        )
    }
}
#[derive(Clone, Debug)]
struct Run {
    id: String,
    kind: String,
    arm: String,
    cohort: usize,
    pair: usize,
    leg: String,
    targets: Vec<String>,
    seeds: Vec<u64>,
}
fn matrix() -> Vec<Run> {
    let mut runs = Vec::new();
    for pair in 1..=5 {
        for leg in ["A", "B"] {
            runs.push(Run {
                id: format!("AA{pair:02}-{leg}"),
                kind: "aa".into(),
                arm: "rho".into(),
                cohort: 0,
                pair,
                leg: leg.into(),
                targets: vec!["target_01".into()],
                seeds: vec![202_610_146_410 + pair as u64],
            });
        }
    }
    for cohort in 1..=3 {
        let offset = (cohort - 1) * 11;
        let mut indices: Vec<_> = (0..N).map(|j| (j + offset) % N).collect();
        if cohort == 2 {
            indices.reverse();
        }
        let ids: Vec<_> = indices
            .iter()
            .map(|i| format!("target_{:02}", i + 1))
            .collect();
        let seeds: Vec<_> = indices
            .iter()
            .map(|i| 202_610_146_501 + ((cohort - 1) * N + i) as u64)
            .collect();
        let arms = if cohort % 2 == 1 {
            vec!["rho", "ic"]
        } else {
            vec!["ic", "rho"]
        };
        for arm in arms {
            runs.push(Run {
                id: format!("C{cohort:02}-{arm}"),
                kind: "cohort".into(),
                arm: arm.into(),
                cohort,
                pair: 0,
                leg: String::new(),
                targets: ids.clone(),
                seeds: seeds.clone(),
            });
        }
        let subset = match cohort {
            1 => 1..=3,
            2 => 4..=6,
            _ => 7..=8,
        };
        for pair in subset {
            let arms = if pair % 2 == 1 {
                vec!["ic", "rho"]
            } else {
                vec!["rho", "ic"]
            };
            for arm in arms {
                runs.push(Run {
                    id: format!("K{pair:02}-{arm}"),
                    kind: "cold".into(),
                    arm: arm.into(),
                    cohort: 0,
                    pair,
                    leg: String::new(),
                    targets: vec![format!("target_{pair:02}")],
                    seeds: vec![202_610_146_600 + pair as u64],
                });
            }
        }
    }
    runs
}
fn run_json(r: &Run, root: &Path) -> Value {
    json!({"run_id":r.id,"kind":r.kind,"arm":r.arm,"cohort_id":r.cohort,"pair_id":r.pair,"leg":r.leg,"target_ids":r.targets,"walk_seeds":if r.arm=="rho"{Some(&r.seeds)}else{None},"rank_seed":if r.arm=="ic"{Some(29)}else{None},"target_file":root.join("inputs").join(format!("{}-points.jsonl",r.id)),"seed_file":if r.arm=="rho"{Some(root.join("inputs").join(format!("{}-seeds.json",r.id)))}else{None},"out_dir":root.join("runs").join(&r.id)})
}
fn prepare(p: &Protocol, root: &Path) -> Result<i32> {
    fs::create_dir(root)?;
    fs::create_dir(root.join("inputs"))?;
    fs::create_dir(root.join("runs"))?;
    let config = root.join("inputs/rho_config.json");
    new_json(&config, &policy())?;
    let runs = matrix();
    let mut generated = BTreeMap::new();
    generated.insert(config.to_string_lossy().to_string(), hash(&config)?);
    for r in &runs {
        let points = r
            .targets
            .iter()
            .map(|id| p.fixture(id).map(|f| f.point.clone()))
            .collect::<Result<Vec<_>>>()?;
        let file = root.join("inputs").join(format!("{}-points.jsonl", r.id));
        new_jsonl(&file, &points)?;
        generated.insert(file.to_string_lossy().into(), hash(&file)?);
        if r.arm == "rho" {
            let f = root.join("inputs").join(format!("{}-seeds.json", r.id));
            new_json(&f, &json!(r.seeds))?;
            generated.insert(f.to_string_lossy().into(), hash(&f)?);
        }
    }
    let manifest = json!({"schema":"n83.matched-manifest/v1","protocol":p.raw,"protocol_sha256":p.digest,"frozen_before_any_child_epoch_ms":epoch_ms(),"provenance":p.provenance()?,"generated_input_sha256":generated,"ordered_runs":runs.iter().map(|r|run_json(r,root)).collect::<Vec<_>>(),"requested_child_processes":32,"requested_target_answers":218,"target_count":N,"shared_setup_repetitions":3,"true_cold_target_subset_count":8,"aa_pairs":5,"rho_policy":policy(),"common_resource_envelope":{"cpu":3,"os_threads":1,"memory_bytes":MEMORY,"timeout_s":BUDGET},"measurement_status":"FROZEN_NOT_RUN","failure_policy":"retain every row; no adaptive target/config/seed replacement; complete-panel ratios require all verified uncensored rows"});
    new_json(&root.join("campaign_manifest.json"), &manifest)?;
    println!(
        "{}",
        json!({"status":"PREPARED_NOT_RUN","manifest":root.join("campaign_manifest.json"),"manifest_sha256":hash(&root.join("campaign_manifest.json"))?,"requested_child_processes":runs.len()})
    );
    Ok(0)
}
fn validate_frozen(p: &Protocol, root: &Path) -> Result<Value> {
    let m = read_json(&root.join("campaign_manifest.json"))?;
    if m["schema"] != "n83.matched-manifest/v1"
        || m["protocol_sha256"] != p.digest
        || m["protocol"] != p.raw
        || m["rho_policy"] != policy()
    {
        return Err(err("frozen campaign identity mismatch"));
    }
    let expected: Vec<_> = matrix().iter().map(|r| run_json(r, root)).collect();
    if m["ordered_runs"] != json!(expected) {
        return Err(err("frozen ordered run matrix mismatch"));
    }
    if m["provenance"]["compiled_source_sha256"] != p.provenance()?["compiled_source_sha256"]
        || m["provenance"]["controller_executable_sha256"] != hash(&std::env::current_exe()?)?
    {
        return Err(err(
            "controller executable/source differs from frozen campaign",
        ));
    }
    if m["provenance"]["input_binary_bindings"] != json!(p.bindings) {
        return Err(err("bound input/executable identity changed"));
    }
    for (file, digest) in m["generated_input_sha256"]
        .as_object()
        .ok_or_else(|| err("generated bindings missing"))?
    {
        if hash(Path::new(file))?.as_str()
            != digest
                .as_str()
                .ok_or_else(|| err("invalid generated digest"))?
        {
            return Err(err(format!("frozen generated input changed: {file}")));
        }
    }
    Ok(m)
}
fn target_record_validate(row: &Value, f: &Fixture, arm: &str, seed: Option<u64>) -> bool {
    let observed = row
        .get("target")
        .or_else(|| row.get("published_q"))
        .and_then(point);
    let k = row.get("recovered_scalar").and_then(scalar);
    let timing = row
        .get("target_ms")
        .and_then(Value::as_f64)
        .is_some_and(|n| n.is_finite() && n >= 0.);
    let basic = observed == point(&f.point)
        && k.as_deref() == Some(&f.scalar)
        && row["group_verified"] == true
        && timing;
    if arm == "rho" {
        basic
            && monitor::rho_timing_valid(row)
            && row["schema"] == "n83-rho-matched-target/v1"
            && row["kind"] == "n83_rho_matched_target"
            && row["status"] == "PASS_VERIFIED_PUBLIC_TARGET"
            && row["verified"] == true
            && row["exit_code"] == 0
            && row["known_answer_input"] == false
            && row["validation_sidecar_opened"] == false
            && row["optimization_config"] == policy()
            && row["walk_seed"].as_u64() == seed
            && row["fresh_target_dependent_DP_walk_and_rng_state"] == true
            && row["cross_target_collision_sharing"] == false
    } else {
        basic
            && row["published_fixture_scalar"].is_null()
            && [
                "target_query_ms",
                "target_pdp_ms",
                "target_relation_check_ms",
                "target_descent_ms",
                "target_recovery_check_ms",
            ]
            .iter()
            .all(|k| row[k].as_f64().is_some_and(|n| n.is_finite() && n >= 0.))
    }
}
fn observe(p: &Protocol, r: &Run, dir: &Path) -> Result<Value> {
    let receipt = read_json(&dir.join("execution.json"))?;
    let semantic = dir.join("targets.jsonl");
    let (target_rows, semantic_parse_error) = if semantic.is_file() {
        match read_jsonl(&semantic) {
            Ok(rows) => (rows, None),
            Err(e) => (Vec::new(), Some(e.to_string())),
        }
    } else {
        (Vec::new(), Some("semantic file was not produced".into()))
    };
    let mut rows = Vec::new();
    let mut all_correct = target_rows.len() == r.targets.len();
    for (index, id) in r.targets.iter().enumerate() {
        let f = p.fixture(id)?;
        let matching: Vec<_> = target_rows
            .iter()
            .filter(|v| {
                v.get("target")
                    .or_else(|| v.get("published_q"))
                    .and_then(point)
                    == point(&f.point)
            })
            .collect();
        let raw = if matching.len() == 1 {
            matching[0].clone()
        } else {
            Value::Null
        };
        let valid = matching.len() == 1
            && raw["fixture_index"].as_u64() == Some(index as u64)
            && target_record_validate(&raw, f, &r.arm, Some(r.seeds[index]));
        let wrong = raw
            .get("recovered_scalar")
            .and_then(scalar)
            .is_some_and(|s| s != f.scalar);
        all_correct &= valid;
        rows.push(json!({"run_id":r.id,"kind":r.kind,"target_id":id,"cohort_id":r.cohort,"pair_id":r.pair,"leg":r.leg,"arm":r.arm,"seed":if r.arm=="rho"{Some(r.seeds[index])}else{None},"rank_seed":if r.arm=="ic"{Some(29)}else{None},"status":if wrong{"WRONG_ANSWER"}else if valid{"COMPLETED_VERIFIED"}else{"TARGET_NOT_COMPLETED_OR_UNVERIFIED"},"verified":valid,"independently_verified":false,"online_wall_s":raw["target_ms"].as_f64().map(|n|n/1000.),"rawtargetrecord":raw}));
    }
    let (stdout, stdout_parse_error) = match read_jsonl(&dir.join("stdout.jsonl")) {
        Ok(rows) => (rows, None),
        Err(e) => (Vec::new(), Some(e.to_string())),
    };
    let summaries: Vec<_> = stdout
        .iter()
        .filter(|v| {
            if r.arm == "ic" {
                v["kind"] == "compact_orbit_dlp_summary"
            } else {
                v["kind"] == "n83_rho_matched_summary"
            }
        })
        .collect();
    let summary = if summaries.len() == 1 {
        summaries[0].clone()
    } else {
        Value::Null
    };
    let summary_ok = if r.arm == "ic" {
        summary["rank"] == 600
            && summary["orbit_columns"] == 600
            && summary["base_hash"] == p.base_hash
    } else {
        summary["optimization_config"] == policy()
            && summary["target_dependent_state_shared"] == false
            && summary["target_independent_curve_normal_basis_constructions"] == 1
    };
    let counts_ok = summary["targets"].as_u64() == Some(r.targets.len() as u64)
        && summary["targets_solved"].as_u64() == Some(r.targets.len() as u64)
        && summary["targets_failed"] == 0;
    let complete =
        receipt["status"] == "COMPLETED_VERIFIED" && all_correct && summary_ok && counts_ok;
    let t = &summary["timing_ms"];
    let prep = if r.arm == "ic" {
        t["setup_total"]
            .as_f64()
            .zip(t["rank_stage"].as_f64())
            .zip(t["linear_algebra"].as_f64())
            .map(|((a, b), c)| a + b + c)
    } else {
        t["setup_total"].as_f64()
    };
    let shared_prep = if r.arm == "ic" {
        prep
    } else {
        t["setup_total"]
            .as_f64()
            .zip(t["all_G_only_jumps"].as_f64())
            .map(|(a, b)| a - b)
    };
    let per_target_prep = if r.arm == "ic" {
        Some(0.)
    } else {
        t["all_G_only_jumps"]
            .as_f64()
            .map(|n| n / r.targets.len() as f64)
    };
    let status = if complete {
        "COMPLETED_VERIFIED"
    } else if rows.iter().any(|r| r["status"] == "WRONG_ANSWER") {
        "WRONG_ANSWER"
    } else if receipt["status"] == "COMPLETED_VERIFIED" {
        "SEMANTIC_VALIDATION_FAILED"
    } else {
        receipt["status"]
            .as_str()
            .unwrap_or("INCOMPLETE_OR_UNVERIFIED")
    };
    Ok(
        json!({"schema":"n83.matched-observation/v1","run_id":r.id,"kind":r.kind,"arm":r.arm,"cohort_id":r.cohort,"pair_id":r.pair,"leg":r.leg,"target_id":if r.targets.len()==1{Some(&r.targets[0])}else{None},"target_ids":r.targets,"target_count":r.targets.len(),"status":status,"all_verified":complete,"independently_verified":false,"parent_known_answer_match_all":all_correct,"summary_ok":summary_ok,"summary_counts_ok":counts_ok,"semantic_parse_error":semantic_parse_error,"stdout_parse_error":stdout_parse_error,"elapsed_process_wall_s":receipt["elapsed_process_wall_s"],"cpu_user_s":receipt["cpu_user_s"],"cpu_system_s":receipt["cpu_system_s"],"cpu_total_s":receipt["cpu_total_s"],"peak_rss_bytes":receipt["peak_rss_bytes"],"producer_preparation_ms":prep,"producer_shared_preparation_ms":shared_prep,"producer_per_target_preparation_ms":per_target_prep,"summary":summary,"target_rows":rows,"raw_process_receipt":receipt,"execution_receipt_sha256":hash(&dir.join("execution.json"))?,"campaign_manifest_sha256":hash(&dir.parent().unwrap().parent().unwrap().join("campaign_manifest.json"))?}),
    )
}
fn validate_existing(p: &Protocol, r: &Run, dir: &Path) -> Result<()> {
    let observation = read_json(&dir.join("observation.json"))?;
    let receipt = read_json(&dir.join("execution.json"))?;
    if observation["execution_receipt_sha256"] != hash(&dir.join("execution.json"))? {
        return Err(err("resume observation/execution hash mismatch"));
    }
    for (key, name) in [
        ("stdout_sha256", "stdout.jsonl"),
        ("stderr_sha256", "stderr.txt"),
    ] {
        if receipt[key].as_str() != Some(hash(&dir.join(name))?.as_str()) {
            return Err(err(format!("resume {name} digest mismatch")));
        }
    }
    let semantic = dir.join("targets.jsonl");
    if semantic.is_file() {
        if receipt["semantic_sha256"].as_str() != Some(hash(&semantic)?.as_str()) {
            return Err(err("resume semantic digest mismatch"));
        }
    } else if !receipt["semantic_sha256"].is_null() {
        return Err(err("resume missing semantic file"));
    }
    let points = dir
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("inputs")
        .join(format!("{}-points.jsonl", r.id));
    if receipt["target_sha256"] != hash(&points)?
        || receipt["binary_sha256"] != hash(if r.arm == "ic" { &p.ic } else { &p.rho })?
        || receipt["resource_envelope"]["timeout_s"].as_f64() != Some(BUDGET as f64)
        || receipt["resource_envelope"]["cpu_affinity"] != json!([3])
        || receipt["resource_envelope"]["address_space_limit_bytes"] != MEMORY
    {
        return Err(err("resume target/binary/common-resource binding mismatch"));
    }
    if observe(p, r, dir)? != observation {
        return Err(err(
            "resume raw records do not reproduce the frozen observation",
        ));
    }
    Ok(())
}
fn run_campaign(p: &Protocol, root: &Path) -> Result<i32> {
    validate_frozen(p, root)?;
    for r in matrix() {
        let dir = root.join("runs").join(&r.id);
        if dir.exists() {
            validate_existing(p, &r, &dir)?;
            println!(
                "{}",
                json!({"run_id":r.id,"status":"EXISTING_FROZEN_RECEIPT_RETAINED"})
            );
            continue;
        }
        let points = root.join("inputs").join(format!("{}-points.jsonl", r.id));
        let semantic = dir.join("targets.jsonl");
        let mut args = vec![
            "--arm".into(),
            if r.arm == "ic" && r.targets.len() > 1 {
                "ic_batch".into()
            } else if r.arm == "rho" && r.targets.len() > 1 {
                "rho_batch".into()
            } else {
                r.arm.clone()
            },
            "--target-count".into(),
            r.targets.len().to_string(),
            "--out-dir".into(),
            dir.to_string_lossy().into(),
            "--target-file".into(),
            points.to_string_lossy().into(),
            "--semantic-file".into(),
            semantic.to_string_lossy().into(),
            "--cwd".into(),
            p.cwd.to_string_lossy().into(),
            "--timeout-s".into(),
            BUDGET.to_string(),
            "--memory-bytes".into(),
            MEMORY.to_string(),
            "--cpu".into(),
            "3".into(),
        ];
        if r.arm == "ic" {
            args.extend(["--expected-base-hash".into(), p.base_hash.clone()]);
        }
        args.push("--".into());
        if r.arm == "ic" {
            args.extend([
                p.ic.to_string_lossy().into(),
                p.base.to_string_lossy().into(),
                points.to_string_lossy().into(),
                "29".into(),
                semantic.to_string_lossy().into(),
            ]);
        } else {
            args.extend([
                p.rho.to_string_lossy().into(),
                points.to_string_lossy().into(),
                root.join("inputs")
                    .join(format!("{}-seeds.json", r.id))
                    .to_string_lossy()
                    .into(),
                root.join("inputs/rho_config.json").to_string_lossy().into(),
                semantic.to_string_lossy().into(),
            ]);
        }
        println!(
            "{}",
            json!({"starting_run":r.id,"kind":r.kind,"arm":r.arm,"target_count":r.targets.len(),"epoch_ms":epoch_ms()})
        );
        let _ = monitor::launch(&args)?;
        let observation = observe(p, &r, &dir)?;
        new_json(&dir.join("observation.json"), &observation)?;
        println!(
            "{}",
            json!({"finished_run":r.id,"status":observation["status"],"process_wall_s":observation["elapsed_process_wall_s"]})
        );
        if monitor::was_interrupted() {
            return Ok(130);
        }
    }
    println!(
        "{}",
        json!({"status":"ALL_DECLARED_CHILDREN_ATTEMPTED","manifest_sha256":hash(&root.join("campaign_manifest.json"))?,"independent_replay_required":true})
    );
    Ok(0)
}
fn collect(p: &Protocol, root: &Path, out: &Path, replay: Option<&Path>) -> Result<i32> {
    validate_frozen(p, root)?;
    let manifest_hash = hash(&root.join("campaign_manifest.json"))?;
    let independently_verified = if let Some(f) = replay {
        let r = read_json(f)?;
        if r["status"] != "PASS"
            || r["checked_target_rows"] != 218
            || r["invalid_rows"] != 0
            || r["checked_aa_semantic_pairs"] != 5
            || r["aa_semantic_gate"] != true
            || r["campaign_manifest_sha256"] != manifest_hash
            || r["protocol_sha256"] != p.digest
        {
            return Err(err(
                "independent replay receipt does not cover all218 bound target rows and five semantic A/A pairs",
            ));
        }
        true
    } else {
        false
    };
    let mut processes = Vec::new();
    let mut targets = Vec::new();
    let mut auxiliary = Vec::new();
    for r in matrix() {
        let dir = root.join("runs").join(&r.id);
        validate_existing(p, &r, &dir)?;
        let mut o = read_json(&dir.join("observation.json"))?;
        if o["campaign_manifest_sha256"] != manifest_hash
            || o["run_id"] != r.id
            || o["arm"] != r.arm
            || o["target_ids"] != json!(r.targets)
        {
            return Err(err("observation/matrix binding mismatch"));
        }
        let mut rows = o["target_rows"]
            .as_array()
            .ok_or_else(|| err("observation target rows missing"))?
            .clone();
        for row in &mut rows {
            row["independently_verified"] = json!(independently_verified);
            if !independently_verified {
                row["verified"] = json!(false);
            }
        }
        if r.kind == "cohort" {
            targets.extend(rows);
        } else {
            auxiliary.extend(rows);
        }
        o.as_object_mut().unwrap().remove("target_rows");
        o["seed"] = if r.arm == "rho" && r.seeds.len() == 1 {
            json!(r.seeds[0])
        } else {
            Value::Null
        };
        o["rank_seed"] = if r.arm == "ic" {
            json!(29)
        } else {
            Value::Null
        };
        o["independently_verified"] = json!(independently_verified);
        if !independently_verified {
            o["all_verified"] = json!(false);
        }
        processes.push(o);
    }
    new_json(
        out,
        &json!({"schema":"n83_matched_statistics_input_v1","process_rows":processes,"target_rows":targets,"auxiliary_target_rows":auxiliary,"independent_verification":replay.map(read_json).transpose()?,"metadata":{"target_ids":p.fixtures.iter().map(|f|&f.id).collect::<Vec<_>>(),"cold_target_ids":p.fixtures.iter().take(8).map(|f|&f.id).collect::<Vec<_>>(),"campaign_manifest_sha256":manifest_hash,"protocol_sha256":p.digest,"independent_replay_required":true,"independent_replay_receipt":replay,"independent_replay_receipt_sha256":replay.map(hash).transpose()?,"independent_verification_receipt_sha256":replay.map(hash).transpose()?,"independently_verified":independently_verified,"requested_processes":32,"requested_target_answers":218,"target_count":32,"cohort_repeats":3,"cold_subset_count":8,"aa_pairs":5,"ICrank_seed":29,"common_resource_envelope":{"cpu":3,"os_threads":1,"memory_bytes":MEMORY,"timeout_s":BUDGET},"IConline_excludes_input_validation_that_rho_includes":true,"ICfive_phase_timer_residuals_retained":true,"no_common_calibrated_operation_unit":true}}),
    )?;
    println!(
        "{}",
        json!({"status":"COLLECTED","out":out,"independently_verified":independently_verified})
    );
    Ok(0)
}
fn cli() -> Result<i32> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let command = args.first().ok_or_else(|| err("missing command"))?;
    let mut opts = BTreeMap::new();
    let mut i = 1;
    while i < args.len() {
        let k = args[i].clone();
        let v = args
            .get(i + 1)
            .ok_or_else(|| err("missing option value"))?
            .clone();
        if opts.insert(k, v).is_some() {
            return Err(err("duplicate option"));
        }
        i += 2;
    }
    let req = |k: &str| {
        opts.get(k)
            .map(PathBuf::from)
            .ok_or_else(|| err(format!("missing {k}")))
    };
    let p = Protocol::load(&req("--protocol")?)?;
    let root = req("--out-dir")?;
    match command.as_str(){"prepare"=>prepare(&p,&root),"run"=>run_campaign(&p,&root),"collect"=>collect(&p,&root,&req("--out")?,opts.get("--independent-replay-receipt").map(Path::new)),_=>Err(err("commands prepare|run|collect --protocol JSON --out-dir DIR [--out NEW_JSON] [--independent-replay-receipt JSON]"))}
}
fn main() {
    match cli() {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            eprintln!("matched runner error: {e}");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn complete_matrix_and_seed_domains_are_predeclared() {
        let m = matrix();
        assert_eq!(m.len(), 32);
        assert_eq!(m.iter().map(|r| r.targets.len()).sum::<usize>(), 218);
        let cohort_rho: Vec<_> = m
            .iter()
            .filter(|r| r.kind == "cohort" && r.arm == "rho")
            .collect();
        let seeds: BTreeSet<_> = cohort_rho.iter().flat_map(|r| r.seeds.iter()).collect();
        assert_eq!(seeds.len(), 96);
        assert_eq!(seeds.first(), Some(&&202_610_146_501));
        assert_eq!(seeds.last(), Some(&&202_610_146_596));
    }
    #[test]
    fn wrong_scalar_point_policy_or_seed_rejects_a_target() {
        let f = Fixture {
            id: "target_01".into(),
            point: json!(["1", "2"]),
            scalar: "3".into(),
        };
        let good = json!({"schema":"n83-rho-matched-target/v1","kind":"n83_rho_matched_target","status":"PASS_VERIFIED_PUBLIC_TARGET","target":["1","2"],"recovered_scalar":"3","group_verified":true,"verified":true,"exit_code":0,"known_answer_input":false,"validation_sidecar_opened":false,"target_ms":1.,"online_total_ms":1.,"online_interval_ns":1000000,"online_phase_sum_ns":999990,"online_timer_boundary_residual_ns":10,"online_phase_ns":{"target_query_and_validation":20,"walk_and_collision":999950,"generic_scalar_replay":20},"optimization_config":policy(),"walk_seed":4,"fresh_target_dependent_DP_walk_and_rng_state":true,"cross_target_collision_sharing":false});
        assert!(target_record_validate(&good, &f, "rho", Some(4)));
        for (k, v) in [
            ("recovered_scalar", json!("8")),
            ("target", json!(["5", "6"])),
            ("walk_seed", json!(5)),
            ("cross_target_collision_sharing", json!(true)),
            ("online_phase_ns", Value::Null),
            ("online_phase_sum_ns", json!(1000000)),
            ("online_timer_boundary_residual_ns", json!(11)),
            ("verified", json!(false)),
        ] {
            let mut bad = good.clone();
            bad[k] = v;
            assert!(!target_record_validate(&bad, &f, "rho", Some(4)));
        }
    }
    #[test]
    fn repeated_aa_seed_is_within_pair_only_and_orders_alternate() {
        let m = matrix();
        for pair in 1..=5 {
            let rows: Vec<_> = m
                .iter()
                .filter(|r| r.kind == "aa" && r.pair == pair)
                .collect();
            assert_eq!(rows.len(), 2);
            assert_eq!(rows[0].seeds, rows[1].seeds);
        }
        assert_eq!(m[10].arm, "rho");
        assert_eq!(m[18].arm, "ic");
        assert_eq!(m[26].arm, "rho");
    }
    #[test]
    fn a_launch_failure_is_preserved_and_can_be_resumed_without_rerunning() {
        let root = std::env::temp_dir().join(format!(
            "n83-matched-launch-failure-{}-{}",
            std::process::id(),
            epoch_ms()
        ));
        fs::create_dir(&root).unwrap();
        fs::create_dir(root.join("inputs")).unwrap();
        fs::create_dir(root.join("runs")).unwrap();
        new_json(
            &root.join("campaign_manifest.json"),
            &json!({"test_only":true}),
        )
        .unwrap();
        let binary = root.join("missing-interpreter");
        fs::write(&binary, "#!/missing/n83-matched-interpreter\nexit 0\n").unwrap();
        fs::set_permissions(&binary, fs::Permissions::from_mode(0o755)).unwrap();
        let r = matrix()[0].clone();
        let points = root.join("inputs").join(format!("{}-points.jsonl", r.id));
        new_jsonl(&points, &[json!(["1", "2"])]).unwrap();
        let dir = root.join("runs").join(&r.id);
        let semantic = dir.join("targets.jsonl");
        let args = vec![
            "--arm".into(),
            "rho".into(),
            "--out-dir".into(),
            dir.to_string_lossy().into(),
            "--target-file".into(),
            points.to_string_lossy().into(),
            "--semantic-file".into(),
            semantic.to_string_lossy().into(),
            "--cwd".into(),
            root.to_string_lossy().into(),
            "--timeout-s".into(),
            BUDGET.to_string(),
            "--memory-bytes".into(),
            MEMORY.to_string(),
            "--cpu".into(),
            "3".into(),
            "--".into(),
            binary.to_string_lossy().into(),
        ];
        assert_eq!(monitor::launch(&args).unwrap(), 1);
        let p = Protocol {
            raw: json!({}),
            digest: String::new(),
            cwd: root.clone(),
            rho: binary.clone(),
            ic: binary,
            base: points,
            base_hash: String::new(),
            fixtures: vec![Fixture {
                id: "target_01".into(),
                point: json!(["1", "2"]),
                scalar: "3".into(),
            }],
            bindings: BTreeMap::new(),
        };
        let observation = observe(&p, &r, &dir).unwrap();
        assert_eq!(observation["status"], "LAUNCH_FAILED");
        assert_eq!(observation["all_verified"], false);
        new_json(&dir.join("observation.json"), &observation).unwrap();
        validate_existing(&p, &r, &dir).unwrap();
        assert_eq!(observation["cpu_total_s"], Value::Null);
        assert_eq!(observation["peak_rss_bytes"], Value::Null);
    }
}
