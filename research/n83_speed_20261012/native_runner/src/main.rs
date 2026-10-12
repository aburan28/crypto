//! Frozen finite scalar-kernel and policy stages for public synthetic degree-83 fixtures.
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
    scalar_kernel: String,
    lanes: u64,
    dp: u64,
    escape: bool,
}
impl Policy {
    fn new(scalar_kernel: &str, lanes: u64, dp: u64, escape: bool) -> Self {
        Self {
            scalar_kernel: scalar_kernel.into(),
            lanes,
            dp,
            escape,
        }
    }
    fn id(&self) -> String {
        format!(
            "{}-L{}-DP{}-E{}",
            self.scalar_kernel,
            self.lanes,
            self.dp,
            u8::from(self.escape)
        )
    }
    fn json(&self) -> Value {
        json!({"backend":"fold","scalar_kernel":self.scalar_kernel,"lanes":self.lanes,"dp_bits":self.dp,"escape_fruitless":self.escape,"step_cap_factor":2000})
    }
    fn from_json(v: &Value) -> Result<Self> {
        let b = v["scalar_kernel"]
            .as_str()
            .ok_or_else(|| err("selected scalar kernel missing"))?;
        let l = v["lanes"]
            .as_u64()
            .ok_or_else(|| err("selected lanes missing"))?;
        let d = v["dp_bits"]
            .as_u64()
            .ok_or_else(|| err("selected DP missing"))?;
        let e = v["escape_fruitless"]
            .as_bool()
            .ok_or_else(|| err("selected escape missing"))?;
        if v["backend"] != "fold"
            || !["reference", "barrett"].contains(&b)
            || ![64, 128, 256].contains(&l)
            || ![6, 8, 10].contains(&d)
            || e
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
    fixtures: PathBuf,
    sidecar: PathBuf,
    generator: PathBuf,
    training: Vec<Fixture>,
    holdout: Vec<Fixture>,
}
impl Protocol {
    fn load(p: &Path) -> Result<Self> {
        let v = read_json(p)?;
        if v["schema"] != "n83.speed-campaign/v1" || v["protocol_id"].as_str().is_none() {
            return Err(err("protocol schema/identity missing"));
        }
        if let Some(expected) = v["human_protocol_sha256"].as_str() {
            let human = path(&v, "human_protocol")?;
            if hash(&human)? != expected {
                return Err(err("human protocol SHA256 mismatch before any child"));
            }
        }
        if v["resource"]["cpu"] != 3
            || v["resource"]["memory_bytes"].as_u64() != Some(MEMORY)
            || v["resource"]["rho_timeout_s"] != 60
        {
            return Err(err("protocol must freeze CPU3,12GiB,rho60s"));
        }
        let cwd = path(&v, "cwd")?;
        let rho = path(&v, "rho_binary")?;
        if let Some(expected) = v["expected_rho_binary_sha256"].as_str() {
            if hash(&rho)? != expected {
                return Err(err("rho binary differs from frozen protocol"));
            }
        }
        let fixtures = path(&v, "fixtures_dir")?;
        let sidecar = path(&v, "validation_sidecar")?;
        let generator = path(&v, "generator_receipt")?;
        for (key, file) in [
            ("expected_validation_sidecar_sha256", &sidecar),
            ("expected_generator_receipt_sha256", &generator),
        ] {
            if let Some(expected) = v[key].as_str() {
                if hash(file)? != expected {
                    return Err(err(format!("{key} mismatch before any child")));
                }
            }
        }
        let gen = read_json(&generator)?;
        if gen["status"] != "PASS"
            || gen["training_target_count"] != 3
            || gen["holdout_target_count"] != 5
        {
            return Err(err("fixture prerequisite not PASS for3training/5holdout"));
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
                if !expected
                    .parse::<u64>()
                    .is_ok_and(|d| d > 0 && d < 8569786107849059)
                    || !points.insert(rows[0].to_string())
                {
                    return Err(err("invalid scalar or duplicate fixture point"));
                }
                if let Some(expected) = v["expected_fixture_file_sha256"][&id].as_str() {
                    if hash(&file)? != expected {
                        return Err(err(format!(
                            "fixture hash differs from frozen protocol: {id}"
                        )));
                    }
                }
                out.push(Fixture {
                    id,
                    path: file,
                    point: rows[0].clone(),
                    expected_scalar: expected,
                });
            }
        }
        Ok(Self {
            raw: v,
            path: p.into(),
            digest: hash(p)?,
            cwd,
            rho,
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
            "rho_binary":self.rho,"rho_binary_sha256":hash(&self.rho)?,
            "fixtures_dir":self.fixtures,"fixtures_sha256":fixture_hashes,"fixture_points":self.training.iter().chain(self.holdout.iter()).map(|f|json!({"fixture_id":f.id,"point":f.point})).collect::<Vec<_>>(),"validation_sidecar_sha256":hash(&self.sidecar)?,"generator_receipt_sha256":hash(&self.generator)?,
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
        .flat_map(|(q, f)| (0..3).map(move |s| (f, 202610126211 + q as u64 * 3 + s)))
        .collect()
}
fn holdout_cases(p: &Protocol) -> Vec<(&Fixture, u64)> {
    p.holdout
        .iter()
        .enumerate()
        .flat_map(|(q, f)| (0..3).map(move |s| (f, 202610126241 + q as u64 * 3 + s)))
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
        "aa" => {
            let target = &p.training[0];
            for pair in 1..=5 {
                for arm in if pair % 2 == 1 {
                    ["a", "b"]
                } else {
                    ["b", "a"]
                } {
                    runs.push(Run {
                        id: format!("AA{pair:02}-{arm}"),
                        target: target.id.clone(),
                        seed: 202610126211,
                        arm: arm.into(),
                        policy: Policy::new("reference", 64, 6, false),
                    });
                }
            }
        }
        "arithmetic" => {
            for (case, (f, seed)) in training_cases(p).into_iter().enumerate() {
                let policies = [
                    Policy::new("reference", 64, 6, false),
                    Policy::new("barrett", 64, 6, false),
                ];
                let order = if case % 2 == 0 { [0, 1] } else { [1, 0] };
                for i in order {
                    let policy = policies[i].clone();
                    runs.push(Run {
                        id: format!("A{:02}-{}", case + 1, policy.id()),
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
            for lanes in [64, 128, 256] {
                for dp in [6, 8, 10] {
                    policies.push(Policy::new(b, lanes, dp, false));
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
            let default = Policy::new("reference", 64, 6, false);
            for (case, (f, seed)) in holdout_cases(p).into_iter().enumerate() {
                let arms = [("current", default.clone()), ("selected", chosen.clone())];
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
        &json!({"schema":"n83.speed-stage/v1","stage":stage,"frozen_before_launch_epoch_ms":epoch_ms(),"planned_launches":runs.len(),"full_block_only_no_adaptive_stopping":true,
        "resource_envelope":p.raw["resource"],"selection_input_sha256":selection_inputs.iter().map(|(p,h)|json!({"path":p,"sha256":h})).collect::<Vec<_>>(),"source_and_input_provenance":p.provenance()?,
"configs":configs,"ordered_runs":runs.iter().map(|r|r.json(root)).collect::<Vec<_>>(),"selection_rule":match stage{"aa"=>"five same-target same-seed A/A pairs; verify exported counters and retain all noise ratios; no selection","holdout"=>"chosen config frozen using separate training cases; all15paired cases retained; no retuning",_=>"all9verified uncapped completions; wronganswer rejects; arithmetic mean min(wall,60) primary; default wins if within2percent best; noeligible blocks nextstage"}}),
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
    let wall = receipt["elapsed_process_wall_s"].as_f64();
    let valid = receipt["status"] == "COMPLETED_VERIFIED"
        && receipt["semantic"]["valid"] == true
        && recovered == Some(f.expected_scalar.as_str())
        && wall.is_some_and(|w| w.is_finite() && w > 0. && w <= RHO_CAP)
        && receipt["semantic"]["producer_online_ms"]
            .as_f64()
            .is_some_and(|t| t.is_finite() && t > 0.)
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
    let charge = if valid {
        wall.unwrap_or(RHO_CAP).min(RHO_CAP)
    } else {
        RHO_CAP
    };
    let record = json!({"schema":"n83.speed-observation/v1","run_id":r.id,"candidate_id":r.policy.id(),"config":r.policy.json(),"target_id":r.target,"seed":r.seed,"arm":r.arm,"status":status,
        "parent_known_answer_scalar_match":recovered==Some(f.expected_scalar.as_str()),"wrong_answer":wrong,"verified":valid,"capped_charge_wall_s":charge,"capped_charge_rule":"verified min(actualwall,60), everyunverified/cap60; nozeroimputation",
        "elapsed_process_wall_s":wall,"producer_online_ms":receipt["semantic"]["producer_online_ms"],"producer_preparation_ms":receipt["semantic"]["producer_preparation_ms"],"cpu_total_s":receipt["cpu_total_s"],"peak_rss_bytes":receipt["peak_rss_bytes"],
"producer_phase_gap_ms":receipt["semantic"]["producer_phase_gap_ms"],"recovered_scalar":recovered,"trajectory":{"walk_seed":receipt["semantic"]["selected_row"]["walk_seed"],"jump_seed":receipt["semantic"]["selected_row"]["jump_seed"],"scalar":recovered,"group_verified":receipt["semantic"]["group_verified"],"producer_status":receipt["semantic"]["selected_row"]["status"],"rho_policy":receipt["semantic"]["selected_row"]["rho_policy"],"outcome":receipt["semantic"]["selected_row"]["outcome"],"preparation_charges":receipt["semantic"]["selected_row"]["preparation_charges"]},
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
    if stage == "arithmetic" {
        let mut cases: BTreeMap<(String, u64), BTreeMap<String, &Value>> = BTreeMap::new();
        for row in rows {
            cases
                .entry((
                    row["target_id"].as_str().unwrap_or("").into(),
                    row["seed"].as_u64().unwrap_or(0),
                ))
                .or_default()
                .insert(
                    row["config"]["scalar_kernel"].as_str().unwrap_or("").into(),
                    row,
                );
        }
        for ((target, seed), case) in cases {
            let wide = case.get("reference");
            let fold = case.get("barrett");
            let equal = wide.zip(fold).is_some_and(|(w, f)| {
                w["status"] == "COMPLETED_VERIFIED"
                    && f["status"] == "COMPLETED_VERIFIED"
                    && w["verified"] == true
                    && f["verified"] == true
                    && w["trajectory"].is_object()
                    && w["trajectory"]["outcome"].is_object()
                    && w["trajectory"] == f["trajectory"]
            });
            if !equal {
                rejected.insert("barrett-L64-DP6-E0".to_string());
            }
            trajectory_checks.push(json!({"target_id":target,"seed":seed,"same_verified_outcome_and_all_charges":equal,"reference_status":wide.map(|r|r["status"].clone()),"barrett_status":fold.map(|r|r["status"].clone()),"comparison_scope":"exported outcome, counters, preparation charges and policies; no full step trace is exported"}));
        }
        if trajectory_checks.len() != 9 {
            rejected.insert("barrett-L64-DP6-E0".into());
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
        let eligible = n == 9
            && verified == 9
            && !wrong
            && !rejected.contains(&id)
            && identities.len() == 9
            && values.iter().all(|r| {
                r["capped_charge_wall_s"]
                    .as_f64()
                    .is_some_and(|c| c.is_finite() && c > 0. && c <= RHO_CAP)
            });
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
    let receipt = json!({"schema":"n83.speed-selection/v1","stage":stage,"status":if available{"SELECTED_FROZEN"}else{"FAILED_NO_COMPLETE_CANDIDATE"},"selection_epoch_ms":epoch_ms(),"protocol_sha256":p.digest,
        "stage_manifest_sha256":hash(&root.join("stage_manifest.json"))?,"stage_records_sha256":hash(&root.join("records.jsonl"))?,"selection_rule":"all9verified uncapped first; thenmin capped arithmetic mean fullprocess; default ifwithin2percentbest; wronganswersreject",
        "default_candidate_id":default.id(),"selected_candidate_id":chosen.as_ref().map(|s|s["candidate_id"].clone()),"selected_config":chosen.as_ref().map(|s|s["config"].clone()),"arithmetic_counter_equivalence_checks":checks,"candidate_scores":scores,"no_best_seed_or_per_target_selection":true});
    new_json(&root.join("selection.json"), &receipt)?;
    Ok(receipt)
}
fn load_selection(p: &Protocol, path: &Path, stage: &str) -> Result<Policy> {
    let s = read_json(path)?;
    if s["schema"] != "n83.speed-selection/v1"
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
    if name != "arithmetic" && name != "aa" {
        let kp = kernel.ok_or_else(|| err("--arithmetic-selection required"))?;
        let k = load_selection(p, kp, "arithmetic")?;
        backend = Some(k.scalar_kernel);
        inputs.push((kp.to_path_buf(), hash(kp)?));
    }
    if name == "holdout" {
        let pp = policy.ok_or_else(|| err("--policy-selection required"))?;
        let s = load_selection(p, pp, "policy")?;
        if Some(s.scalar_kernel.as_str()) != backend.as_deref() {
            return Err(err(
                "policy selected scalar kernel differs from arithmetic selection",
            ));
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
        &json!({"schema":"n83.speed-stage-completion/v1","status":disposition,"planned_runs":plan.len(),"retained_runs":rows.len(),"verified_runs":rows.iter().filter(|r|r["status"]=="COMPLETED_VERIFIED").count(),"finished_epoch_ms":epoch_ms(),"charged_actual_process_wall_s":rows.iter().filter_map(|r|r["elapsed_process_wall_s"].as_f64()).sum::<f64>(),"capped_training_selection_charge_wall_s":rows.iter().filter_map(|r|r["capped_charge_wall_s"].as_f64()).sum::<f64>()}),
    )?;
    if !all {
        return Ok(1);
    }
    if name == "aa" {
        let summary = statistics::summarize_aa(&rows);
        let valid = summary["all_complete_verified_and_counter_equal"] == true;
        new_json(&root.join("noise_summary.json"), &summary)?;
        println!(
            "{}",
            json!({"stage":"aa","valid":valid,"noise_summary_path":root.join("noise_summary.json"),"noise_summary_sha256":hash(&root.join("noise_summary.json"))?})
        );
        return Ok(if valid { 0 } else { 1 });
    }
    if name == "arithmetic" || name == "policy" {
        let d = Policy::new(
            if name == "arithmetic" {
                "reference"
            } else {
                backend.as_deref().unwrap()
            },
            64,
            6,
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
fn analyze(holdout: &Path, out: &Path) -> Result<i32> {
    let rows = read_jsonl(&holdout.join("records.jsonl"))?;
    let manifest = read_json(&holdout.join("stage_manifest.json"))?;
    let completion = read_json(&holdout.join("completion.json"))?;
    if manifest["stage"] != "holdout"
        || completion["planned_runs"] != 30
        || completion["retained_runs"] != 30
        || completion["status"] != "FULL_BLOCK_COMPLETED"
    {
        return Err(err("not a complete frozen30-run holdout stage"));
    }
    let planned = manifest["ordered_runs"]
        .as_array()
        .ok_or_else(|| err("frozen ordered runs missing"))?;
    if planned.len() != 30 || rows.len() != 30 {
        return Err(err("holdout requires all30rows"));
    }
    for (declared, row) in planned.iter().zip(&rows) {
        for k in [
            "run_id",
            "target_id",
            "seed",
            "arm",
            "candidate_id",
            "config",
        ] {
            if declared[k] != row[k] {
                return Err(err(format!(
                    "holdout row identity differs from frozen plan: {k}"
                )));
            }
        }
    }
    let mut result = statistics::summarize(&rows);
    result["input_provenance"] = json!({"holdout_manifest_sha256":hash(&holdout.join("stage_manifest.json"))?,"holdout_records_sha256":hash(&holdout.join("records.jsonl"))?});
    result["holdout_completion_raw"] = completion;
    result["holdout_manifest_raw"] = manifest;
    result["analysis_source_sha256"] = json!(hash_bytes(include_bytes!("statistics.rs"))?);
    if let Some(parent) = holdout.parent() {
        let mut costs = Vec::new();
        let mut complete = true;
        for (name, expected) in [
            ("aa", 10),
            ("arithmetic", 18),
            ("policy", 81),
            ("holdout", 30),
        ] {
            let stage = parent.join(name);
            let record_path = stage.join("records.jsonl");
            let done_path = stage.join("completion.json");
            if !record_path.is_file() || !done_path.is_file() {
                complete = false;
                continue;
            }
            let records = read_jsonl(&record_path)?;
            let done = read_json(&done_path)?;
            complete &= records.len() == expected && done["status"] == "FULL_BLOCK_COMPLETED";
            costs.push(json!({"stage":name,"planned_launches":expected,"observed_launches":records.len(),"verified_launches":records.iter().filter(|r|r["verified"]==true).count(),
                "actual_child_process_wall_s":records.iter().filter_map(|r|r["elapsed_process_wall_s"].as_f64()).sum::<f64>(),"actual_child_cpu_total_s":records.iter().filter_map(|r|r["cpu_total_s"].as_f64()).sum::<f64>(),
                "capped_selection_charge_wall_s":records.iter().filter_map(|r|r["capped_charge_wall_s"].as_f64()).sum::<f64>(),"records_sha256":hash(&record_path)?,"completion_sha256":hash(&done_path)?}));
        }
        result["campaign_cost"] = json!({"complete_139_launch_campaign":complete &&costs.len()==4,"stages":costs,"actual_child_process_wall_s":if complete &&costs.len()==4 {Some(costs.iter().filter_map(|c|c["actual_child_process_wall_s"].as_f64()).sum::<f64>())}else{None},
            "actual_child_cpu_total_s":if complete &&costs.len()==4 {Some(costs.iter().filter_map(|c|c["actual_child_cpu_total_s"].as_f64()).sum::<f64>())}else{None},"scope":"Sum of actual monitored child process costs across A/A,arithmetic,policy,holdout; excludes builds, external verification and parent inter-run orchestration. Failed actual costs remain; capped selection cost is separate."});
    }
    new_json(out, &result)?;
    println!(
        "{}",
        json!({"statistics_path":out,"primary":result["primary_current_over_selected_process_wall"],"all_verified":result["all_verified"]})
    );
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
        return analyze(&req("--holdout-dir")?, &req("--out")?);
    }
    let p = Protocol::load(&req("--protocol")?)?;
    let root = req("--out-dir")?;
    match cmd.as_str() {
        "aa" | "arithmetic" | "policy" | "holdout" => stage(
            &p,
            cmd,
            &root,
            opts.get("--arithmetic-selection").map(Path::new),
            opts.get("--policy-selection").map(Path::new),
        ),
        _ => Err(err("unknown command")),
    }
}
fn main() {
    match main_run() {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            eprintln!("speedrunnererror:{e}\nCommands aa/arithmetic/policy/holdout --protocol JSON --out-dir NEWDIR [--arithmetic-selection JSON] [--policy-selection JSON]; analyze --holdout-dir DIR --out NEWJSON");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn row(id: &str, case: u64, cost: f64, verified: bool) -> Value {
        json!({"candidate_id":id,"target_id":format!("training_{}",case/3+1),"seed":case+1,"status":if verified{"COMPLETED_VERIFIED"}else{"TIMEOUT"},"verified":verified,"wrong_answer":false,"capped_charge_wall_s":cost,"producer_online_ms":cost*1000.,"config":Policy::new("reference",64,6,false).json()})
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
    fn incomplete_arithmetic_equivalence_blocks_barrett() {
        let mut rows = Vec::new();
        for i in 0..9 {
            let mut r = row("barrett-L64-DP6-E0", i, 1., true);
            r["config"] = Policy::new("barrett", 64, 6, false).json();
            rows.push(r);
        }
        let (_, chosen, checks) = candidate_scores(&rows, "reference-L64-DP6-E0", "arithmetic");
        assert!(chosen.is_none());
        assert_eq!(checks.len(), 9);
    }
    #[test]
    fn changed_arithmetic_counter_rejects_barrett_even_if_faster() {
        let mut rows = Vec::new();
        for i in 0..9 {
            for (kernel, cost) in [("reference", 2.), ("barrett", 1.)] {
                let policy = Policy::new(kernel, 64, 6, false);
                let mut r = row(&policy.id(), i, cost, true);
                r["config"] = policy.json();
                r["trajectory"] = json!({"outcome":{"walk_steps":100+i}});
                if kernel == "barrett" && i == 0 {
                    r["trajectory"]["outcome"]["walk_steps"] = json!(999);
                }
                rows.push(r);
            }
        }
        let (_, chosen, checks) = candidate_scores(&rows, "reference-L64-DP6-E0", "arithmetic");
        assert_eq!(chosen.unwrap()["candidate_id"], "reference-L64-DP6-E0");
        assert_eq!(
            checks
                .iter()
                .filter(|c| c["same_verified_outcome_and_all_charges"] == false)
                .count(),
            1
        );
    }
    #[test]
    fn policy_grid_cannot_accept_prior_default_or_escaped_policy() {
        assert!(Policy::from_json(&Policy::new("reference", 32, 4, false).json()).is_err());
        assert!(Policy::from_json(&Policy::new("barrett", 64, 6, true).json()).is_err());
        assert!(Policy::from_json(&Policy::new("barrett", 256, 10, false).json()).is_ok());
    }
    #[test]
    fn zero_cost_is_not_an_eligible_completion() {
        let rows = (0..9).map(|i| row("zero", i, 0., true)).collect::<Vec<_>>();
        assert!(candidate_scores(&rows, "zero", "policy").1.is_none());
    }
}
