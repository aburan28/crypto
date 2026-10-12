//! Frozen finite optimization stages for public synthetic legacy n83 fixtures.
mod monitor;
mod statistics;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const RHO_CAP: f64 = 60.;
const MEMORY: u64 = 12 * 1024 * 1024 * 1024;
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
    let mut rows = Vec::new();
    for (i, l) in BufReader::new(File::open(p)?).lines().enumerate() {
        let l = l?;
        if l.trim().is_empty() {
            continue;
        }
        rows.push(
            serde_json::from_str(&l).map_err(|e| err(format!("{}:{}:{e}", p.display(), i + 1)))?,
        );
    }
    Ok(rows)
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
fn path(v: &Value, k: &str) -> Result<PathBuf> {
    v[k].as_str()
        .map(PathBuf::from)
        .ok_or_else(|| err(format!("protocol {k} missing")))
}
#[derive(Clone, Debug)]
struct Fixture {
    id: String,
    path: PathBuf,
    point: Value,
    expected_scalar: String,
}
#[derive(Clone, Debug)]
struct Policy {
    backend: String,
    lanes: u64,
    dp: u64,
    escape: bool,
}
impl Policy {
    fn new(backend: &str, lanes: u64, dp: u64, escape: bool) -> Self {
        Self {
            backend: backend.into(),
            lanes,
            dp,
            escape,
        }
    }
    fn id(&self) -> String {
        format!(
            "{}-L{}-DP{}-E{}",
            self.backend,
            self.lanes,
            self.dp,
            u8::from(self.escape)
        )
    }
    fn json(&self) -> Value {
        json!({"backend":self.backend,"lanes":self.lanes,"dp_bits":self.dp,"escape_fruitless":self.escape,"step_cap_factor":2000})
    }
    fn from_json(v: &Value) -> Result<Self> {
        let b = v["backend"]
            .as_str()
            .ok_or_else(|| err("selected backend missing"))?;
        let l = v["lanes"]
            .as_u64()
            .ok_or_else(|| err("selected lanes missing"))?;
        let d = v["dp_bits"]
            .as_u64()
            .ok_or_else(|| err("selected DP missing"))?;
        let e = v["escape_fruitless"]
            .as_bool()
            .ok_or_else(|| err("selected escape missing"))?;
        if !["wide", "fold"].contains(&b)
            || ![16, 32, 64].contains(&l)
            || ![4, 6].contains(&d)
            || v["step_cap_factor"] != 2000
        {
            return Err(err("invalid selected policy"));
        }
        Ok(Self::new(b, l, d, e))
    }
}
struct Protocol {
    raw: Value,
    path: PathBuf,
    digest: String,
    cwd: PathBuf,
    rho: PathBuf,
    ic: PathBuf,
    base: PathBuf,
    base_hash: String,
    fixtures: PathBuf,
    sidecar: PathBuf,
    generator: PathBuf,
    training: Vec<Fixture>,
    holdout: Vec<Fixture>,
}
impl Protocol {
    fn load(p: &Path) -> Result<Self> {
        let v = read_json(p)?;
        if v["schema"] != "n83.optimization-campaign/v1" || v["protocol_id"].as_str().is_none() {
            return Err(err("protocol schema/identity missing"));
        }
        if let Some(expected)=v["human_protocol_sha256"].as_str(){
            let human=path(&v,"human_protocol")?;
            if hash(&human)?!=expected{return Err(err("human protocol SHA256 mismatch before any child"));}
        }
        if v["resource"]["cpu"] != 3
            || v["resource"]["memory_bytes"].as_u64() != Some(MEMORY)
            || v["resource"]["rho_timeout_s"] != 60
            || v["resource"]["ic_timeout_s"] != 1800
        {
            return Err(err("protocol must freeze CPU3,12GiB,rho60s,IC1800s"));
        }
        let cwd = path(&v, "cwd")?;
        let rho = path(&v, "rho_binary")?;
        let ic = path(&v, "ic_binary")?;
        let base = path(&v, "base_file")?;
        let fixtures = path(&v, "fixtures_dir")?;
        let sidecar = path(&v, "validation_sidecar")?;
        let generator = path(&v, "generator_receipt")?;
        let base_hash = v["expected_base_hash"]
            .as_str()
            .ok_or_else(|| err("base hash missing"))?
            .to_string();
        let gen = read_json(&generator)?;
        if gen["status"] != "PASS"
            || gen["training_target_count"] != 3
            || gen["holdout_target_count"] != 5
        {
            return Err(err("fixture prerequisite not PASS for3training/5holdout"));
        }
        if gen["retained_base_label_hash"].as_str() != Some(base_hash.as_str())
            || gen["retained_base_file_sha256"].as_str() != Some(hash(&base)?.as_str())
        {
            return Err(err(
                "retained base label/file hash differsfrom independentfixture receipt",
            ));
        }
        let s = read_json(&sidecar)?;
        if s["must_not_be_passed_to_measured_children"] != true {
            return Err(err("sidecar validation-only declaration missing"));
        }
        let entries = s["entries"]
            .as_array()
            .ok_or_else(|| err("sidecar entries missing"))?;
        let mut points = BTreeSet::new();
        let mut training = Vec::new();
        let mut holdout = Vec::new();
        for (split, count, out) in [("training", 3, &mut training), ("holdout", 5, &mut holdout)] {
            for i in 1..=count {
                let id = format!("{split}_{i}");
                let file = fixtures.join(format!("{id}.jsonl"));
                let rows = read_jsonl(&file)?;
                if rows.len() != 1 || point(&rows[0]).is_none() {
                    return Err(err(format!("{id} is not one point")));
                }
                let matches: Vec<_> = entries.iter().filter(|e| e["fixture_id"] == id).collect();
                if matches.len() != 1
                    || matches[0]["split"] != split
                    || matches[0]["point"] != rows[0]
                {
                    return Err(err(format!("{id} sidecar point/identity mismatch")));
                }
                let expected = matches[0]["scalar_validation_only"]
                    .as_str()
                    .ok_or_else(|| err("sidecar scalar missing"))?
                    .to_string();
                if expected.parse::<u64>().is_err() || !points.insert(rows[0].to_string()) {
                    return Err(err("invalid scalar or duplicate fixture point"));
                }
                out.push(Fixture {
                    id,
                    path: file,
                    point: rows[0].clone(),
                    expected_scalar: expected,
                });
            }
        }
        let batch = read_jsonl(&fixtures.join("holdout_5_targets.jsonl"))?;
        if batch != holdout.iter().map(|f| f.point.clone()).collect::<Vec<_>>() {
            return Err(err("ICbatch does not match ordered5holdout points"));
        }
        Ok(Self {
            raw: v,
            path: p.into(),
            digest: hash(p)?,
            cwd,
            rho,
            ic,
            base,
            base_hash,
            fixtures,
            sidecar,
            generator,
            training,
            holdout,
        })
    }
    fn fixture(&self, id: &str) -> Result<&Fixture> {
        self.training
            .iter()
            .chain(self.holdout.iter())
            .find(|f| f.id == id)
            .ok_or_else(|| err(format!("unknown fixture {id}")))
    }
    fn provenance(&self) -> Result<Value> {
        let exe = std::env::current_exe()?;
        let mut fixture_hashes = BTreeMap::new();
        for f in self.training.iter().chain(self.holdout.iter()) {
            fixture_hashes.insert(f.id.clone(), hash(&f.path)?);
        }
        let uname = Command::new("/usr/bin/uname").arg("-a").output()?;
        let cpu_info = fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
        let parent_host = self.raw["host_context_file"].as_str().map(PathBuf::from);
        Ok(
            json!({"protocol_path":self.path,"protocol_sha256":self.digest,"executable":exe,"executable_sha256":hash(&exe)?,
            "compiled_source_sha256":{"main":hash_bytes(include_bytes!("main.rs"))?,"monitor":hash_bytes(include_bytes!("monitor.rs"))?,"statistics":hash_bytes(include_bytes!("statistics.rs"))?,"manifest":hash_bytes(include_bytes!("../Cargo.toml"))?},
            "rho_binary":self.rho,"rho_binary_sha256":hash(&self.rho)?,"ic_binary":self.ic,"ic_binary_sha256":hash(&self.ic)?,"base_file_sha256":hash(&self.base)?,
            "fixtures_sha256":fixture_hashes,"holdout_batch_sha256":hash(&self.fixtures.join("holdout_5_targets.jsonl"))?,"validation_sidecar_sha256":hash(&self.sidecar)?,"generator_receipt_sha256":hash(&self.generator)?,
            "fresh_host":{"unix_epoch_ms":epoch_ms(),"uname":String::from_utf8_lossy(&uname.stdout).trim(),"cpu_model":cpu_info.lines().find(|l|l.starts_with("model name")),"cpuinfo_sha256":hash_bytes(cpu_info.as_bytes())?},
            "previous_host_context":if let Some(p)=parent_host.as_ref(){Some(read_json(p)?)}else{None},"previous_host_context_sha256":if let Some(p)=parent_host.as_ref(){Some(hash(p)?)}else{None},
            "monitor_poll_ms":5,"monitor_delay":"Nominal0-5ms plus scheduling/poll-body; no strict bound under arbitrary scheduling","scalar_sidecar_opened_by_parent_only":true}),
        )
    }
}
#[derive(Clone)]
struct Run {
    id: String,
    target: String,
    seed: u64,
    arm: String,
    policy: Policy,
}
impl Run {
    fn json(&self, root: &Path) -> Value {
        json!({"run_id":self.id,"target_id":self.target,"seed":self.seed,"arm":self.arm,"candidate_id":self.policy.id(),"config":self.policy.json(),"config_file":root.join("configs").join(format!("{}.json",self.policy.id())),"out_dir":root.join("runs").join(&self.id)})
    }
}
fn training_cases(p: &Protocol) -> Vec<(&Fixture, u64)> {
    p.training
        .iter()
        .enumerate()
        .flat_map(|(q, f)| (0..3).map(move |s| (f, 202610128311 + q as u64 * 3 + s)))
        .collect()
}
fn holdout_cases(p: &Protocol) -> Vec<(&Fixture, u64)> {
    p.holdout
        .iter()
        .enumerate()
        .flat_map(|(q, f)| (0..3).map(move |s| (f, 202610128341 + q as u64 * 3 + s)))
        .collect()
}
fn stage_plan(
    stage: &str,
    p: &Protocol,
    backend: Option<&str>,
    selected: Option<&Policy>,
) -> Result<Vec<Run>> {
    let mut runs = Vec::new();
    match stage {
        "kernel" => {
            for (case, (f, seed)) in training_cases(p).into_iter().enumerate() {
                let policies = [
                    Policy::new("wide", 32, 4, false),
                    Policy::new("fold", 32, 4, false),
                ];
                let order = if case % 2 == 0 { [0, 1] } else { [1, 0] };
                for i in order {
                    let policy = policies[i].clone();
                    runs.push(Run {
                        id: format!("K{:02}-{}", case + 1, policy.id()),
                        target: f.id.clone(),
                        seed,
                        arm: policy.id(),
                        policy,
                    });
                }
            }
        }
        "policy" => {
            let b = backend.ok_or_else(|| err("policy stage backend missing"))?;
            let mut policies = Vec::new();
            for lanes in [16, 32, 64] {
                for dp in [4, 6] {
                    for escape in [false, true] {
                        policies.push(Policy::new(b, lanes, dp, escape));
                    }
                }
            }
            for (case, (f, seed)) in training_cases(p).into_iter().enumerate() {
                for offset in 0..policies.len() {
                    let policy = policies[(case + offset) % policies.len()].clone();
                    runs.push(Run {
                        id: format!("P{:02}-{}", case + 1, policy.id()),
                        target: f.id.clone(),
                        seed,
                        arm: policy.id(),
                        policy,
                    });
                }
            }
        }
        "holdout" => {
            let chosen = selected.ok_or_else(|| err("holdout selected policy missing"))?;
            let default = Policy::new("wide", 32, 4, false);
            for (case, (f, seed)) in holdout_cases(p).into_iter().enumerate() {
                let arms = [("default", default.clone()), ("selected", chosen.clone())];
                let order = if case % 2 == 0 { [0, 1] } else { [1, 0] };
                for i in order {
                    let (arm, policy) = arms[i].clone();
                    runs.push(Run {
                        id: format!("H{:02}-{arm}", case + 1),
                        target: f.id.clone(),
                        seed,
                        arm: arm.into(),
                        policy,
                    });
                }
            }
        }
        _ => return Err(err("unknown stage")),
    }
    Ok(runs)
}
fn freeze_stage(
    root: &Path,
    stage: &str,
    p: &Protocol,
    runs: &[Run],
    selection_inputs: &[(PathBuf, String)],
) -> Result<()> {
    fs::create_dir(root)?;
    fs::create_dir(root.join("configs"))?;
    fs::create_dir(root.join("runs"))?;
    let mut configs = BTreeMap::new();
    for r in runs {
        configs.insert(r.policy.id(), r.policy.json());
    }
    for (id, c) in &configs {
        new_json(&root.join("configs").join(format!("{id}.json")), c)?;
    }
    new_json(&root.join("protocol.json"), &p.raw)?;
    new_json(
        &root.join("stage_manifest.json"),
        &json!({"schema":"n83.optimization-stage/v1","stage":stage,"frozen_before_launch_epoch_ms":epoch_ms(),"planned_launches":runs.len(),"full_block_only_no_adaptive_stopping":true,
        "resource_envelope":p.raw["resource"],"selection_input_sha256":selection_inputs.iter().map(|(p,h)|json!({"path":p,"sha256":h})).collect::<Vec<_>>(),"source_and_input_provenance":p.provenance()?,
        "configs":configs,"ordered_runs":runs.iter().map(|r|r.json(root)).collect::<Vec<_>>(),"selection_rule":"all9verified uncapped completions; wronganswer rejects; arithmetic mean min(wall,60) primary; default wins if within2percent best; noeligible blocks nextstage"}),
    )?;
    Ok(())
}
fn run_rho(p: &Protocol, root: &Path, r: &Run) -> Result<Value> {
    let f = p.fixture(&r.target)?;
    let dir = root.join("runs").join(&r.id);
    let cfg = root.join("configs").join(format!("{}.json", r.policy.id()));
    let args = vec![
        "--arm".into(),
        "rho".into(),
        "--out-dir".into(),
        dir.to_string_lossy().into(),
        "--target-file".into(),
        f.path.to_string_lossy().into(),
        "--cwd".into(),
        p.cwd.to_string_lossy().into(),
        "--timeout-s".into(),
        "60".into(),
        "--memory-bytes".into(),
        MEMORY.to_string(),
        "--cpu".into(),
        "3".into(),
        "--".into(),
        p.rho.to_string_lossy().into(),
        f.path.to_string_lossy().into(),
        r.seed.to_string(),
        cfg.to_string_lossy().into(),
        dir.join("producer.json").to_string_lossy().into(),
    ];
    eprintln!(
        "{}",
        json!({"starting_run":r.id,"target_id":r.target,"seed":r.seed,"candidate":r.policy.id(),"epoch_ms":epoch_ms()})
    );
    let _ = monitor::launch(&args)?;
    let receipt = read_json(&dir.join("execution.json"))?;
    let recovered = receipt["semantic"]["recovered_scalar"].as_str();
    let wrong = recovered.is_some_and(|s| s != f.expected_scalar);
    let valid = receipt["status"] == "COMPLETED_VERIFIED"
        && receipt["semantic"]["valid"] == true
        && recovered == Some(f.expected_scalar.as_str())
        && receipt["semantic"]["selected_row"]["optimization_config"] == r.policy.json();
    let status = if wrong {
        "WRONG_ANSWER"
    } else if valid {
        "COMPLETED_VERIFIED"
    } else if receipt["status"] == "COMPLETED_VERIFIED" {
        "SEMANTIC_VALIDATION_FAILED"
    } else {
        receipt["status"].as_str().unwrap_or("PROCESS_FAILED")
    };
    let wall = receipt["elapsed_process_wall_s"].as_f64();
    let charge = if valid {
        wall.unwrap_or(RHO_CAP).min(RHO_CAP)
    } else {
        RHO_CAP
    };
    let record = json!({"schema":"n83.optimization-observation/v1","run_id":r.id,"candidate_id":r.policy.id(),"config":r.policy.json(),"target_id":r.target,"seed":r.seed,"arm":r.arm,"status":status,
        "parent_known_answer_scalar_match":recovered==Some(f.expected_scalar.as_str()),"wrong_answer":wrong,"verified":valid,"capped_charge_wall_s":charge,"capped_charge_rule":"verified min(actualwall,60), everyunverified/cap60; nozeroimputation",
        "elapsed_process_wall_s":wall,"producer_online_ms":receipt["semantic"]["producer_online_ms"],"producer_preparation_ms":receipt["semantic"]["producer_preparation_ms"],"cpu_total_s":receipt["cpu_total_s"],"peak_rss_bytes":receipt["peak_rss_bytes"],
        "producer_phase_gap_ms":receipt["semantic"]["producer_phase_gap_ms"],"recovered_scalar":recovered,"trajectory":{"scalar":recovered,"group_verified":receipt["semantic"]["group_verified"],"producer_status":receipt["semantic"]["selected_row"]["status"],"rho_policy":receipt["semantic"]["selected_row"]["rho_policy"],"outcome":receipt["semantic"]["selected_row"]["outcome"],"preparation_charges":receipt["semantic"]["selected_row"]["preparation_charges"]},
        "process_receipt_path":dir.join("execution.json"),"process_receipt_sha256":hash(&dir.join("execution.json"))?,"target_file_sha256":hash(&f.path)?,"config_file_sha256":hash(&cfg)?,"resource_envelope":receipt["resource_envelope"],"raw_process_receipt":receipt});
    new_json(&dir.join("observation.json"), &record)?;
    Ok(record)
}
fn candidate_scores(
    rows: &[Value],
    default_id: &str,
    stage: &str,
) -> (Vec<Value>, Option<Value>, Vec<Value>) {
    let mut grouped: BTreeMap<String, Vec<&Value>> = BTreeMap::new();
    for row in rows {
        grouped
            .entry(row["candidate_id"].as_str().unwrap_or("").into())
            .or_default()
            .push(row);
    }
    let mut trajectory_checks = Vec::new();
    let mut rejected = BTreeSet::new();
    if stage == "kernel" {
        let mut cases: BTreeMap<(String, u64), BTreeMap<String, &Value>> = BTreeMap::new();
        for row in rows {
            cases
                .entry((
                    row["target_id"].as_str().unwrap_or("").into(),
                    row["seed"].as_u64().unwrap_or(0),
                ))
                .or_default()
                .insert(row["config"]["backend"].as_str().unwrap_or("").into(), row);
        }
        for ((target, seed), case) in cases {
            let wide = case.get("wide");
            let fold = case.get("fold");
            let equal = wide.zip(fold).is_some_and(|(w, f)| {
                w["status"] == "COMPLETED_VERIFIED"
                    && f["status"] == "COMPLETED_VERIFIED"
                    && w["trajectory"] == f["trajectory"]
            });
            if !equal {
                rejected.insert("fold-L32-DP4-E0".to_string());
            }
            trajectory_checks.push(json!({"target_id":target,"seed":seed,"same_verified_outcome_and_all_charges":equal,"wide_status":wide.map(|r|r["status"].clone()),"fold_status":fold.map(|r|r["status"].clone())}));
        }
        if trajectory_checks.len() != 9 {
            rejected.insert("fold-L32-DP4-E0".into());
        }
    }
    let mut scores = Vec::new();
    for (id, values) in grouped {
        let n = values.len();
        let verified = values
            .iter()
            .filter(|r| r["status"] == "COMPLETED_VERIFIED" && r["verified"] == true)
            .count();
        let wrong = values.iter().any(|r| r["wrong_answer"] == true);
        let mean = if n > 0 {
            Some(
                values
                    .iter()
                    .map(|r| r["capped_charge_wall_s"].as_f64().unwrap_or(RHO_CAP))
                    .sum::<f64>()
                    / n as f64,
            )
        } else {
            None
        };
        let online: Vec<_> = values
            .iter()
            .filter(|r| r["verified"] == true)
            .filter_map(|r| r["producer_online_ms"].as_f64())
            .collect();
        let mut identities = BTreeSet::new();
        for r in &values {
            identities.insert((r["target_id"].to_string(), r["seed"].to_string()));
        }
        let eligible =
            n == 9 && verified == 9 && !wrong && !rejected.contains(&id) && identities.len() == 9;
        scores.push(json!({"candidate_id":id,"config":values.first().map(|r|r["config"].clone()),"observations":n,"verified_completions":verified,"wrong_answer":wrong,"trajectory_rejected":rejected.contains(&id),"eligible":eligible,
            "primary_capped_arithmetic_mean_process_wall_s":mean,"secondary_arithmetic_mean_online_ms":if online.len()==9{Some(online.iter().sum::<f64>()/9.)}else{None},"all_case_rows":values}));
    }
    let eligible: Vec<_> = scores.iter().filter(|s| s["eligible"] == true).collect();
    let best = eligible.iter().min_by(|a, b| {
        a["primary_capped_arithmetic_mean_process_wall_s"]
            .as_f64()
            .unwrap()
            .total_cmp(
                &b["primary_capped_arithmetic_mean_process_wall_s"]
                    .as_f64()
                    .unwrap(),
            )
            .then_with(|| a["candidate_id"].as_str().cmp(&b["candidate_id"].as_str()))
    });
    let chosen = best.map(|best| {
        let min = best["primary_capped_arithmetic_mean_process_wall_s"]
            .as_f64()
            .unwrap();
        eligible
            .iter()
            .find(|s| {
                s["candidate_id"] == default_id
                    && s["primary_capped_arithmetic_mean_process_wall_s"]
                        .as_f64()
                        .unwrap()
                        <= min * 1.02
            })
            .copied()
            .unwrap_or(best)
            .clone()
    });
    (scores, chosen, trajectory_checks)
}
fn selection(
    p: &Protocol,
    root: &Path,
    stage: &str,
    rows: &[Value],
    default: &Policy,
) -> Result<Value> {
    let (scores, chosen, checks) = candidate_scores(rows, &default.id(), stage);
    let available = chosen.is_some();
    let receipt = json!({"schema":"n83.optimization-selection/v1","stage":stage,"status":if available{"SELECTED_FROZEN"}else{"FAILED_NO_COMPLETE_CANDIDATE"},"selection_epoch_ms":epoch_ms(),"protocol_sha256":p.digest,
        "stage_manifest_sha256":hash(&root.join("stage_manifest.json"))?,"stage_records_sha256":hash(&root.join("records.jsonl"))?,"selection_rule":"all9verified uncapped first; thenmin capped arithmetic mean fullprocess; default ifwithin2percentbest; wronganswersreject",
        "default_candidate_id":default.id(),"selected_candidate_id":chosen.as_ref().map(|s|s["candidate_id"].clone()),"selected_config":chosen.as_ref().map(|s|s["config"].clone()),"kernel_trajectory_checks":checks,"candidate_scores":scores,"no_best_seed_or_per_target_selection":true});
    new_json(&root.join("selection.json"), &receipt)?;
    Ok(receipt)
}
fn load_selection(p: &Protocol, path: &Path, stage: &str) -> Result<Policy> {
    let s = read_json(path)?;
    if s["schema"] != "n83.optimization-selection/v1"
        || s["stage"] != stage
        || s["status"] != "SELECTED_FROZEN"
        || s["protocol_sha256"] != p.digest
    {
        return Err(err(format!(
            "{stage} selection noteligible or differentprotocol"
        )));
    }
    Policy::from_json(&s["selected_config"])
}
fn stage(
    p: &Protocol,
    name: &str,
    root: &Path,
    kernel: Option<&Path>,
    policy: Option<&Path>,
) -> Result<i32> {
    let mut inputs = Vec::new();
    let mut selected = None;
    let mut backend = None;
    if name != "kernel" {
        let kp = kernel.ok_or_else(|| err("--kernel-selection required"))?;
        let k = load_selection(p, kp, "kernel")?;
        backend = Some(k.backend);
        inputs.push((kp.to_path_buf(), hash(kp)?));
    }
    if name == "holdout" {
        let pp = policy.ok_or_else(|| err("--policy-selection required"))?;
        let s = load_selection(p, pp, "policy")?;
        if Some(s.backend.as_str()) != backend.as_deref() {
            return Err(err("policy selectedbackend differsfrom kernelselection"));
        }
        selected = Some(s);
        inputs.push((pp.to_path_buf(), hash(pp)?));
    }
    let plan = stage_plan(name, p, backend.as_deref(), selected.as_ref())?;
    freeze_stage(root, name, p, &plan, &inputs)?;
    let mut records = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(root.join("records.jsonl"))?;
    let mut rows = Vec::new();
    for run in &plan {
        let row = run_rho(p, root, run)?;
        writeln!(records, "{}", row)?;
        records.flush()?;
        rows.push(row);
        if monitor::was_interrupted() {
            break;
        }
    }
    records.sync_all()?;
    drop(records);
    fs::set_permissions(
        root.join("records.jsonl"),
        fs::Permissions::from_mode(0o444),
    )?;
    let all = rows.len() == plan.len();
    let disposition = if all {
        "FULL_BLOCK_COMPLETED"
    } else {
        "INTERRUPTED_PARTIAL_BLOCK"
    };
    new_json(
        &root.join("completion.json"),
        &json!({"schema":"n83.optimization-stage-completion/v1","status":disposition,"planned_runs":plan.len(),"retained_runs":rows.len(),"verified_runs":rows.iter().filter(|r|r["status"]=="COMPLETED_VERIFIED").count(),"finished_epoch_ms":epoch_ms()}),
    )?;
    if !all {
        return Ok(1);
    }
    if name == "kernel" || name == "policy" {
        let d = Policy::new(
            if name == "kernel" {
                "wide"
            } else {
                backend.as_deref().unwrap()
            },
            32,
            4,
            false,
        );
        let s = selection(p, root, name, &rows, &d)?;
        println!(
            "{}",
            json!({"stage":name,"status":s["status"],"selected_candidate_id":s["selected_candidate_id"],"selection_path":root.join("selection.json"),"selection_sha256":hash(&root.join("selection.json"))?})
        );
        return Ok(if s["status"] == "SELECTED_FROZEN" {
            0
        } else {
            1
        });
    }
    println!(
        "{}",
        json!({"stage":name,"out_dir":root,"full_block_completed":true,"verified_rows":rows.iter().filter(|r|r["status"]=="COMPLETED_VERIFIED").count()})
    );
    Ok(0)
}
fn ic_batch(p: &Protocol, root: &Path) -> Result<i32> {
    fs::create_dir(root)?;
    let batch = p.fixtures.join("holdout_5_targets.jsonl");
    new_json(
        &root.join("stage_manifest.json"),
        &json!({"schema":"n83.optimization-ic-batch/v1","frozen_before_launch_epoch_ms":epoch_ms(),"protocol":p.raw,"provenance":p.provenance()?,"target_count":5,"rank_seed":29,"shared_preparation_not_cold_single_target":true,"stock_extract_finite_index_scan_no20millioncap":true}),
    )?;
    let dir = root.join("process");
    let target_file = dir.join("targets.jsonl");
    let args = vec![
        "--arm".into(),
        "ic_batch".into(),
        "--target-count".into(),
        "5".into(),
        "--out-dir".into(),
        dir.to_string_lossy().into(),
        "--target-file".into(),
        batch.to_string_lossy().into(),
        "--semantic-file".into(),
        target_file.to_string_lossy().into(),
        "--expected-base-hash".into(),
        p.base_hash.clone(),
        "--cwd".into(),
        p.cwd.to_string_lossy().into(),
        "--timeout-s".into(),
        "1800".into(),
        "--memory-bytes".into(),
        MEMORY.to_string(),
        "--cpu".into(),
        "3".into(),
        "--".into(),
        p.ic.to_string_lossy().into(),
        p.base.to_string_lossy().into(),
        batch.to_string_lossy().into(),
        "29".into(),
        target_file.to_string_lossy().into(),
    ];
    let _ = monitor::launch(&args)?;
    let receipt = read_json(&dir.join("execution.json"))?;
    let target_records = if target_file.is_file() {
        read_jsonl(&target_file)?
    } else {
        Vec::new()
    };
    let mut rows = Vec::new();
    for f in &p.holdout {
        let matches: Vec<_> = target_records
            .iter()
            .filter(|r| {
                r.get("target")
                    .or_else(|| r.get("published_q"))
                    .and_then(point)
                    == point(&f.point)
            })
            .collect();
        let r = if matches.len() == 1 {
            matches[0].clone()
        } else {
            Value::Null
        };
        let recovered = r.get("recovered_scalar").and_then(scalar);
        let correct = recovered.as_deref() == Some(f.expected_scalar.as_str());
        rows.push(json!({"target_id":f.id,"target":f.point,"target_ms":r["target_ms"],"group_verified":r["group_verified"],"recovered_scalar":recovered,"parent_known_answer_scalar_match":correct,
            "status":if r.is_null(){"TARGET_NOT_COMPLETED"}else if !correct{"WRONG_ANSWER"}else if r["group_verified"]==true{"COMPLETED_VERIFIED"}else{"SEMANTIC_VALIDATION_FAILED"},"raw_target_record":r}));
    }
    let complete = receipt["status"] == "COMPLETED_VERIFIED"
        && receipt["semantic"]["valid"] == true
        && rows.iter().all(|r| r["status"] == "COMPLETED_VERIFIED");
    let output = json!({"schema":"n83.optimization-ic-batch-observations/v1","status":if complete{"COMPLETED_VERIFIED"}else{"INCOMPLETE_OR_UNVERIFIED"},"shared_preparation_not_cold_single_target":true,
        "producer_preparation_ms":receipt["semantic"]["producer_preparation_ms"],"elapsed_process_wall_s":receipt["elapsed_process_wall_s"],"cpu_total_s":receipt["cpu_total_s"],"peak_rss_bytes":receipt["peak_rss_bytes"],"rows":rows,"raw_process_receipt":receipt,
        "process_receipt_sha256":hash(&dir.join("execution.json"))?,"ICcap1800s_rhocap60s_unequal":true,"no_common_calibrated_operation_unit":true});
    new_json(&root.join("observations.json"), &output)?;
    println!("{}", output);
    Ok(if complete { 0 } else { 1 })
}
fn analyze(holdout: &Path, ic_dir: &Path, out: &Path) -> Result<i32> {
    let rows = read_jsonl(&holdout.join("records.jsonl"))?;
    let manifest = read_json(&holdout.join("stage_manifest.json"))?;
    let completion = read_json(&holdout.join("completion.json"))?;
    let ic = read_json(&ic_dir.join("observations.json"))?;
    let ic_manifest = read_json(&ic_dir.join("stage_manifest.json"))?;
    if manifest["stage"] != "holdout" || completion["planned_runs"] != 30 {
        return Err(err("not frozen30-run holdout stage"));
    }
    let ic_raw = ic["rows"].as_array().ok_or_else(|| err("ICrows missing"))?;
    let mut ic_admitted = ic_raw.clone();
    if ic["status"] != "COMPLETED_VERIFIED" {
        for row in &mut ic_admitted {
            row["group_verified"] = json!(false);
            row["batch_complete_verified"] = json!(false);
        }
    }
    let mut result = statistics::summarize(&rows, &ic_admitted);
    result["input_provenance"] = json!({"holdout_manifest_sha256":hash(&holdout.join("stage_manifest.json"))?,"holdout_records_sha256":hash(&holdout.join("records.jsonl"))?,"ic_observations_sha256":hash(&ic_dir.join("observations.json"))?,"ic_manifest_sha256":hash(&ic_dir.join("stage_manifest.json"))?});
    result["holdout_completion_raw"] = completion;
    result["holdout_manifest_raw"] = manifest;
    result["ic_manifest_raw"] = ic_manifest;
    result["ic_batch_raw"] = ic;
    result["analysis_source_sha256"] = json!(hash_bytes(include_bytes!("statistics.rs"))?);
    new_json(out, &result)?;
    println!("{}", result);
    Ok(0)
}
fn main_run() -> Result<i32> {
    let a: Vec<_> = std::env::args().skip(1).collect();
    let cmd = a.first().ok_or_else(|| err("commandmissing"))?;
    if (a.len() - 1) % 2 != 0 {
        return Err(err("options need values"));
    }
    let mut opts = BTreeMap::new();
    for x in a[1..].chunks(2) {
        if opts.insert(x[0].as_str(), x[1].as_str()).is_some() {
            return Err(err("duplicateoption"));
        }
    }
    let req = |k: &str| {
        opts.get(k)
            .map(PathBuf::from)
            .ok_or_else(|| err(format!("missing{k}")))
    };
    if cmd == "analyze" {
        return analyze(&req("--holdout-dir")?, &req("--ic-dir")?, &req("--out")?);
    }
    let p = Protocol::load(&req("--protocol")?)?;
    let root = req("--out-dir")?;
    match cmd.as_str() {
        "kernel" | "policy" | "holdout" => stage(
            &p,
            cmd,
            &root,
            opts.get("--kernel-selection").map(Path::new),
            opts.get("--policy-selection").map(Path::new),
        ),
        "ic-batch" => ic_batch(&p, &root),
        _ => Err(err("unknown command")),
    }
}
fn main() {
    match main_run() {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            eprintln!("optimizationrunnererror:{e}\nCommands kernel/policy/holdout/ic-batch --protocol JSON --out-dir NEWDIR [--kernel-selection JSON] [--policy-selection JSON]; analyze --holdout-dir DIR --ic-dir DIR --out NEWJSON");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn row(id: &str, case: u64, cost: f64, verified: bool) -> Value {
        json!({"candidate_id":id,"target_id":format!("training_{}",case/3+1),"seed":case+1,"status":if verified{"COMPLETED_VERIFIED"}else{"TIMEOUT"},"verified":verified,"wrong_answer":false,"capped_charge_wall_s":cost,"producer_online_ms":cost*1000.,"config":Policy::new("wide",32,4,false).json()})
    }
    #[test]
    fn incomplete_candidate_never_beats_complete_candidate() {
        let mut rows = Vec::new();
        for i in 0..9 {
            rows.push(row("complete", i, 20., true));
            rows.push(row("incomplete", i, if i == 8 { 60. } else { 1. }, i != 8));
        }
        let (_, chosen, _) = candidate_scores(&rows, "none", "policy");
        assert_eq!(chosen.unwrap()["candidate_id"], "complete");
    }
    #[test]
    fn default_wins_declared_two_percent_tie() {
        let mut rows = Vec::new();
        for i in 0..9 {
            rows.push(row("default", i, 10.1, true));
            rows.push(row("fast", i, 10., true));
        }
        let (_, chosen, _) = candidate_scores(&rows, "default", "policy");
        assert_eq!(chosen.unwrap()["candidate_id"], "default");
    }
    #[test]
    fn wrong_answer_rejects_otherwise_complete_candidate() {
        let mut rows = Vec::new();
        for i in 0..9 {
            let mut r = row("wrong", i, 1., true);
            if i == 0 {
                r["wrong_answer"] = json!(true);
            }
            rows.push(r);
        }
        let (_, chosen, _) = candidate_scores(&rows, "wrong", "policy");
        assert!(chosen.is_none());
    }
    #[test]
    fn incomplete_kernel_trajectory_blocks_fold() {
        let mut rows = Vec::new();
        for i in 0..9 {
            let mut r = row("fold-L32-DP4-E0", i, 1., true);
            r["config"] = Policy::new("fold", 32, 4, false).json();
            rows.push(r);
        }
        let (_, chosen, checks) = candidate_scores(&rows, "wide-L32-DP4-E0", "kernel");
        assert!(chosen.is_none());
        assert_eq!(checks.len(), 9);
    }
}
