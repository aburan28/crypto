//! Frozen, one-use execution of the disclosed synthetic prepared F5 control.
//! Fresh sampling, ordinary collection and cross-method claims are excluded.
use super::{f5_target, prepared_f5, sat_control::native};
use crypto_lib::cryptanalysis::prepared_sat_control::{canonical_sha, sha256};
use native::{load, read, require, save};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{fs, path::Path, process::Command, time::Instant};

const PREPARATION: &str = "research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json";
const SCOPE: &str = "disclosed-n17-native-prepared-f5-control-v1";
const STATE: &str = "edbff76da6442b9f2e5e8235682c9bf1052465f310c765ba8a37d018a60bf107";
const WORKER: &str = "ic_tournament_worker";
const FILES: [&str; 4] = [
    "config.json",
    "preparation.json",
    "job.json",
    "host-context.json",
];

#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Control {
    schema_version: u64,
    question: String,
    target: [u64; 2],
    algorithm_seed: u64,
    max_queries: u64,
    controller_timeout_ms: u64,
}
impl Control {
    fn validate(&self) -> Result<(), String> {
        require(
            self.schema_version == 1
                && self.question == "native-f5-known-input-transport-and-recovery-control-v1"
                && self.target == [52411, 72106]
                && (1..=8).contains(&self.max_queries)
                && (1_000..=900_000).contains(&self.controller_timeout_ms),
            "outside the disclosed n17 F5 control envelope",
        )
    }
}
fn config(capsule: &Path) -> Result<Control, String> {
    let cfg: Control =
        serde_json::from_value(load(&capsule.join("config.json"))?).map_err(|e| e.to_string())?;
    cfg.validate()?;
    Ok(cfg)
}
fn require_executable_registration(registration: &Value) -> Result<(), String> {
    require(
        registration["validation_only"] == false,
        "F5 build-validation capsule cannot dispatch or admit a scientific job",
    )
}
fn environment() -> Vec<(String, String)> {
    [
        ("LC_ALL", "C"),
        ("RAYON_NUM_THREADS", "1"),
        ("IC_ARTIFACT_CACHE", "off"),
        ("IC_F2_BACKEND", "cpu"),
    ]
    .into_iter()
    .map(|(k, v)| (k.into(), v.into()))
    .collect()
}
fn reject_ancestor_cargo_configuration(source: &Path) -> Result<(), String> {
    // Cargo merges ancestor configuration even when CARGO_HOME is empty.
    // Only the inventoried snapshot's own .cargo/config.toml is admissible.
    for parent in source.ancestors().skip(1) {
        for name in [".cargo/config", ".cargo/config.toml"] {
            require(
                !parent.join(name).try_exists().map_err(|e| e.to_string())?,
                "unregistered ancestor Cargo configuration; choose an isolated capsule parent",
            )?;
        }
    }
    Ok(())
}
fn mathematical_job(prep: &Value, cfg: &Control) -> Result<Value, String> {
    cfg.validate()?;
    prepared_f5::verify(prep)?;
    let coordinate = |v: &Value| -> Result<String, String> {
        v.as_u64()
            .or_else(|| v.as_str().and_then(|s| s.parse::<u64>().ok()))
            .map(|v| v.to_string())
            .ok_or_else(|| "invalid prepared coordinate".into())
    };
    let point = |p: &Value| -> Result<[String; 2], String> {
        let p = p
            .as_array()
            .filter(|p| p.len() == 2)
            .ok_or("invalid prepared point")?;
        Ok([coordinate(&p[0])?, coordinate(&p[1])?])
    };
    let state = &prep["record"];
    let geometry = prep["certificate"]["inputs"]["base"]
        .as_array()
        .ok_or("missing geometry")?
        .iter()
        .map(point)
        .collect::<Result<Vec<_>, _>>()?;
    let columns = state["column_logs"]
        .as_array()
        .ok_or("missing column logs")?
        .iter()
        .map(|column| {
            Ok(json!({"point":point(&column["point"])? ,
            "log":column["log"].as_u64().ok_or("invalid column log")?.to_string()}))
        })
        .collect::<Result<Vec<Value>, String>>()?;
    Ok(
        json!({"mode":"ic","degree":17,"curve_a":1,"target_seeds":[],
        "public_targets":[["52411","72106"]],"algorithm_seed":cfg.algorithm_seed,
        "factor_base":{"kind":"standard_subspace","dimension":6},"exclusive_phases":true,
        "config":{"solver":"f5","linear_algebra":"dense","summands":3,
            "groebner_degree":3,"node_budget":8192,"conflict_budget":100000,"batch_trials":1,
            "max_trials":cfg.max_queries,"collection_window":null,"factor_base":null,
            "factor_base_cube_root":false,"factor_base_orbits":null,"rho_parallel_walks":32,
            "sparse":crypto_lib::cryptanalysis::koblitz_sparse_la::SparseSolveOptions::default()},
        "prepared":{"mathematical_state_sha256":STATE,"factor_base":geometry,"columns":columns}}),
    )
}
fn check_capsule(capsule: &Path) -> Result<Value, String> {
    let registration = load(&capsule.join("registration.json"))?;
    let seal = load(&capsule.join("seal.json"))?;
    require(
        seal["registration_sha256"] == canonical_sha(&registration)?,
        "F5 registration seal differs",
    )?;
    require(
        registration["schema_version"] == 1
            && registration["scope"] == SCOPE
            && registration["hardware"] == json!({"os":"macos","architecture":"aarch64"})
            && registration["fresh_paired_qualification"] == false
            && registration["headline_eligible"] == false
            && registration["online_speedup"].is_null()
            && registration["validation_only"].is_boolean()
            && registration["runtime_environment"] == json!(environment()),
        "F5 registration scope differs",
    )?;
    native::check_tree(&capsule.join("immutable"), &registration["immutable_files"])?;
    for name in FILES {
        require(
            registration[name] == sha256(&read(&capsule.join(name), 16 * 1024 * 1024)?),
            "F5 capsule input differs",
        )?;
    }
    let cfg = config(capsule)?;
    let prep = load(&capsule.join("preparation.json"))?;
    require(
        load(&capsule.join("job.json"))? == mathematical_job(&prep, &cfg)?,
        "F5 worker input is not mathematical-only",
    )?;
    require(
        registration["preparation_admission"] == prepared_f5::verify(&prep)?,
        "F5 preparation admission differs",
    )?;
    let bin = capsule.join("immutable/bin");
    require(
        registration["worker_sha256"] == sha256(&read(&bin.join(WORKER), 128 * 1024 * 1024)?)
            && registration["auditor_sha256"]
                == sha256(&read(&bin.join("icprog"), 128 * 1024 * 1024)?),
        "F5 binary binding differs",
    )?;
    Ok(registration)
}

/// Freeze builds from retained sources and offline dependencies, without a job.
pub fn freeze(
    root: &Path,
    out: &Path,
    config_file: &Path,
    cargo: &Path,
    rustc: &Path,
    host: &Path,
    validation_only: bool,
) -> Result<String, String> {
    native::enforce_hardware()?;
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let cargo = cargo.canonicalize().map_err(|e| e.to_string())?;
    let rustc = rustc.canonicalize().map_err(|e| e.to_string())?;
    let cfg: Control = serde_json::from_value(load(config_file)?).map_err(|e| e.to_string())?;
    cfg.validate()?;
    let prep = load(&root.join(PREPARATION))?;
    let admission = prepared_f5::verify(&prep)?;
    let host = load(host)?;
    require(host["os"] == "macos" && host["architecture"] == "aarch64"
        && host["physical_hardware"] == true && host["calibrated_performance_environment"] == false,
        "F5 host context must identify physical disclosed-control hardware without calibration claims")?;
    let head = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(&root)
        .output()
        .map_err(|e| e.to_string())?;
    require(head.status.success(), "source checkout lacks HEAD")?;
    let dirty = Command::new("git")
        .args([
            "status",
            "--porcelain",
            "--untracked-files=all",
            "--",
            "src",
            "examples/ic_tournament_worker.rs",
            "Cargo.toml",
            "docs/ic/calibration.json",
            "docs/ecbench/schema.sql",
            "docs/curves/registry.json",
        ])
        .current_dir(&root)
        .output()
        .map_err(|e| e.to_string())?;
    require(
        dirty.status.success() && dirty.stdout.is_empty(),
        "commit every F5 worker/controller/auditor source before freezing",
    )?;
    fs::create_dir(out).map_err(|e| e.to_string())?;
    let out = out.canonicalize().map_err(|e| e.to_string())?;
    let immutable = out.join("immutable");
    fs::create_dir(&immutable).map_err(|e| e.to_string())?;
    let source = immutable.join("source");
    fs::create_dir(&source).map_err(|e| e.to_string())?;
    native::copy_tree(&root.join("src"), &source.join("src"))?;
    for name in [
        "Cargo.toml",
        "Cargo.lock",
        "examples/ic_tournament_worker.rs",
        "docs/ic/calibration.json",
        "docs/ecbench/schema.sql",
        "docs/curves/registry.json",
    ] {
        let path = source.join(name);
        fs::create_dir_all(path.parent().ok_or("missing source parent")?)
            .map_err(|e| e.to_string())?;
        native::create(&path, &read(&root.join(name), 16 * 1024 * 1024)?)?;
    }
    let receipts = immutable.join("build-receipts");
    fs::create_dir(&receipts).map_err(|e| e.to_string())?;
    let path = std::env::var("PATH").map_err(|e| e.to_string())?;
    let vendor_env = vec![("PATH".into(), path.clone()), ("LC_ALL".into(), "C".into())];
    super::sat_control::native_step(
        &cargo,
        &[
            "vendor".into(),
            "--locked".into(),
            "--offline".into(),
            source.join("vendor").to_string_lossy().into_owned(),
        ],
        &root,
        &receipts.join("vendor.log"),
        &vendor_env,
    )?;
    fs::create_dir(source.join(".cargo")).map_err(|e| e.to_string())?;
    native::create(&source.join(".cargo/config.toml"), b"[source.crates-io]\nreplace-with = \"vendored-sources\"\n[source.vendored-sources]\ndirectory = \"vendor\"\n")?;
    let home = out.join("build-home");
    fs::create_dir(&home).map_err(|e| e.to_string())?;
    let target = out.join("build-target");
    let snapshot = native::inventory(&source)?;
    let source_sha = canonical_sha(&snapshot)?;
    let rustc_sha = sha256(&read(&rustc, 128 * 1024 * 1024)?);
    let cargo_sha = sha256(&read(&cargo, 128 * 1024 * 1024)?);
    let build = json!({"schema_version":1,"source_manifest_sha256":source_sha,
        "rustc_sha256":rustc_sha,"cargo_sha256":cargo_sha,"profile":"release","features":"default",
        "hardware":{"os":"macos","architecture":"aarch64"},"algorithm_flags":"none"});
    let build_sha = canonical_sha(&build)?;
    save(&receipts.join("build-identity.json"), &build)?;
    reject_ancestor_cargo_configuration(&source)?;
    let env = vec![
        ("PATH".into(), path),
        ("LC_ALL".into(), "C".into()),
        ("RUSTC".into(), rustc.to_string_lossy().into_owned()),
        ("CARGO_HOME".into(), home.to_string_lossy().into_owned()),
        (
            "CARGO_TARGET_DIR".into(),
            target.to_string_lossy().into_owned(),
        ),
        (
            "IC_GENERIC_SOURCE_MANIFEST_SHA256".into(),
            source_sha.clone(),
        ),
        ("IC_GENERIC_BUILD_SHA256".into(), build_sha.clone()),
    ];
    super::sat_control::native_step(
        &rustc,
        &["--version".into(), "--verbose".into()],
        &source,
        &receipts.join("rustc.log"),
        &env,
    )?;
    super::sat_control::native_step(
        &cargo,
        &["--version".into(), "--verbose".into()],
        &source,
        &receipts.join("cargo.log"),
        &env,
    )?;
    super::sat_control::native_step(
        &cargo,
        &[
            "build".into(),
            "--locked".into(),
            "--offline".into(),
            "--release".into(),
            "--example".into(),
            WORKER.into(),
            "--bin".into(),
            "icprog".into(),
        ],
        &source,
        &receipts.join("build.log"),
        &env,
    )?;
    native::check_tree(&source, &snapshot)?;
    let bin = immutable.join("bin");
    fs::create_dir(&bin).map_err(|e| e.to_string())?;
    for (name, path) in [
        (WORKER, target.join("release/examples").join(WORKER)),
        ("icprog", target.join("release/icprog")),
    ] {
        native::create(&bin.join(name), &read(&path, 128 * 1024 * 1024)?)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(bin.join(name), fs::Permissions::from_mode(0o755))
                .map_err(|e| e.to_string())?;
        }
    }
    let identity = json!({"schema_version":1,"source_manifest_sha256":source_sha,"build_sha256":build_sha,
        "target_arch":"aarch64","target_os":"macos"});
    super::sat_control::native_step(
        &bin.join(WORKER),
        &["--build-identity".into()],
        &out,
        &receipts.join("worker-identity.log"),
        &environment(),
    )?;
    require(
        serde_json::from_slice::<Value>(&read(&receipts.join("worker-identity.log"), 64 * 1024)?)
            .map_err(|e| e.to_string())?
            == identity,
        "frozen worker compiled identity differs",
    )?;
    require(
        sha256(&read(&rustc, 128 * 1024 * 1024)?) == rustc_sha
            && sha256(&read(&cargo, 128 * 1024 * 1024)?) == cargo_sha,
        "F5 compiler or Cargo changed during the frozen build",
    )?;
    save(&out.join("config.json"), &json!(cfg))?;
    save(&out.join("preparation.json"), &prep)?;
    save(&out.join("job.json"), &mathematical_job(&prep, &cfg)?)?;
    save(&out.join("host-context.json"), &host)?;
    let mut registration = json!({"schema_version":1,"scope":SCOPE,
        "source_commit":String::from_utf8(head.stdout).map_err(|e|e.to_string())?.trim(),
        "hardware":{"os":"macos","architecture":"aarch64"},"immutable_files":native::inventory(&immutable)?,
        "worker_sha256":sha256(&read(&bin.join(WORKER),128*1024*1024)?),
        "auditor_sha256":sha256(&read(&bin.join("icprog"),128*1024*1024)?),"worker_build_identity":identity,
        "preparation_admission":admission,"runtime_environment":environment(),
        "primary_question":"disclosed synthetic known-input transport and recovery; no fresh-target speed claim",
        "ordinary_queries_executed":0,"fresh_paired_qualification":false,"headline_eligible":false,
        "promotion_eligible":false,"candidate_id":null,"workload_id":null,"run_id":null,"validation_only":validation_only,
        "reference":null,"operation_calibration":null,"online_speedup":null,"status":"registered-not-dispatched"});
    for name in FILES {
        registration[name] = json!(sha256(&read(&out.join(name), 16 * 1024 * 1024)?));
    }
    save(&out.join("registration.json"), &registration)?;
    let seal = json!({"registration_sha256":canonical_sha(&registration)?});
    save(&out.join("seal.json"), &seal)?;
    check_capsule(&out)?;
    serde_json::to_string_pretty(&seal).map_err(|e| e.to_string())
}

/// The claim is durable before the sole native scientific launch. No resume path.
pub fn execute(capsule: &Path, execution: &Path, expected: &str) -> Result<String, String> {
    native::enforce_hardware()?;
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    let registration = check_capsule(&capsule)?;
    require_executable_registration(&registration)?;
    let registration_sha = canonical_sha(&registration)?;
    require(
        registration_sha == expected,
        "F5 dispatch lacks the exact external registration seal",
    )?;
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    require(
        registration["auditor_sha256"] == sha256(&read(&own, 128 * 1024 * 1024)?),
        "F5 dispatch must use the frozen controller",
    )?;
    let cfg = config(&capsule)?;
    let job = read(&capsule.join("job.json"), 1_048_576)?;
    fs::create_dir(execution).map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    save(
        &capsule.join("consumed.json"),
        &json!({"registration_sha256":registration_sha,
        "execution":execution,"status":"consumed-before-launch"}),
    )?;
    let child = native::measured_child_request(native::ChildRequest {
        program: &capsule.join("immutable/bin").join(WORKER),
        args: &[],
        cwd: &execution,
        stem: &execution.join("worker"),
        deadline_ms: cfg.controller_timeout_ms,
        ledger: &execution.join("worker-pids"),
        helper: false,
        input: Some(&job),
        environment: &environment(),
    });
    let drain = native::drain_ledger(&execution.join("worker-pids"));
    let gate = check_capsule(&capsule);
    let mut terminal = json!({"schema_version":1,"registration_sha256":registration_sha,
        "source_gate_passed":gate.is_ok(),"source_gate_error":gate.err(),
        "child_drain_passed":drain.is_ok(),"child_drain_error":drain.err()});
    match child {
        Ok(child) => {
            terminal["worker"] = json!({"exit_code":child.exit_code,"timed_out":child.timed_out,"receipt":child.receipt})
        }
        Err(error) => terminal["worker_transport_error"] = json!(error),
    }
    save(&execution.join("terminal.json"), &terminal)?;
    serde_json::to_string_pretty(&terminal).map_err(|e| e.to_string())
}

fn verify_transport(
    registration: &Value,
    cfg: &Control,
    capsule: &Path,
    execution: &Path,
    terminal: &Value,
) -> Result<Value, String> {
    require(
        terminal["schema_version"] == 1
            && terminal["registration_sha256"] == canonical_sha(registration)?
            && terminal["source_gate_passed"] == true
            && terminal["child_drain_passed"] == true,
        "F5 terminal source/drain gate failed",
    )?;
    let claim = load(&capsule.join("consumed.json"))?;
    require(
        claim["registration_sha256"] == canonical_sha(registration)?
            && claim["status"] == "consumed-before-launch"
            && claim["execution"] == json!(execution.canonicalize().map_err(|e| e.to_string())?),
        "F5 audit differs from the one-use claim",
    )?;
    let receipt = load(&execution.join("worker.receipt.json"))?;
    let output = read(&execution.join("worker.stdout"), 8 * 1024 * 1024)?;
    let error = read(&execution.join("worker.stderr"), 8 * 1024 * 1024)?;
    let input = read(&capsule.join("job.json"), 1_048_576)?;
    let env = environment()
        .into_iter()
        .collect::<std::collections::BTreeMap<_, _>>();
    require(
        receipt["argv"] == json!([capsule.join("immutable/bin").join(WORKER)])
            && receipt["cwd"] == json!(execution.canonicalize().map_err(|e| e.to_string())?)
            && receipt["deadline_ms"] == cfg.controller_timeout_ms
            && receipt["environment"] == json!(env)
            && receipt["executable_sha256_before"] == registration["worker_sha256"]
            && receipt["executable_sha256_after"] == registration["worker_sha256"]
            && receipt["stdin_sha256"] == sha256(&input)
            && receipt["stdin_bytes"] == input.len()
            && receipt["stdin_written_bytes"] == input.len()
            && receipt["stdin_write_error"].is_null()
            && receipt["stdout_sha256"] == sha256(&output)
            && receipt["stdout_bytes"] == output.len()
            && receipt["stderr_sha256"] == sha256(&error)
            && receipt["stderr_bytes"] == error.len()
            && receipt["timed_out"] == false
            && receipt["output_limit"] == false
            && receipt["process_group_drain_requested"] == true
            && receipt["process_group_drain_confirmed"] == true
            && terminal["worker"]["receipt"] == receipt
            && terminal["worker"]["timed_out"] == false
            && terminal["worker"]["exit_code"] == receipt["exit_code"],
        "F5 exact worker invocation/output receipt differs",
    )?;
    let pid = receipt["pid"]
        .as_u64()
        .filter(|&pid| pid > 1 && pid <= u32::MAX as u64)
        .ok_or("invalid F5 child PID")?;
    let ledger = String::from_utf8(read(&execution.join("worker-pids"), 64 * 1024)?)
        .map_err(|e| e.to_string())?;
    require(
        ledger == format!("start {pid}\ndone {pid}\n"),
        "F5 child ledger is incomplete or contains extra launches",
    )?;
    let report: Value = serde_json::from_slice(&output).map_err(|e| e.to_string())?;
    require(
        report["generic_build"] == registration["worker_build_identity"]
            && !report
                .as_object()
                .ok_or("F5 output is not an object")?
                .contains_key("test_fixture"),
        "F5 producer build or fixture status differs",
    )?;
    let complete = report["status"] == "complete";
    require(
        receipt["exit_code"] == if complete { 0 } else { 2 },
        "F5 producer exit disagrees with completion",
    )?;
    Ok(report)
}

/// Frozen independent checker: sources, transport, every target attempt and scalar.
/// No solver or archived producer is ever called by this entry point.
pub fn audit(capsule: &Path, execution: &Path, out: &Path) -> Result<String, String> {
    let started = Instant::now();
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    let registration = check_capsule(&capsule)?;
    require_executable_registration(&registration)?;
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    let before = sha256(&read(&own, 128 * 1024 * 1024)?);
    require(
        registration["auditor_sha256"] == before,
        "F5 audit must use its preexecution-frozen independent checker",
    )?;
    let cfg = config(&capsule)?;
    let terminal = load(&execution.join("terminal.json"))?;
    let result =
        verify_transport(&registration, &cfg, &capsule, &execution, &terminal).and_then(|report| {
            f5_target::verify(
                &load(&capsule.join("preparation.json"))?,
                cfg.algorithm_seed,
                cfg.max_queries,
                &report,
            )
        });
    let mut receipt = match result {
        Ok(mathematics) => json!({"status":"PASS_NATIVE_SOURCE_BOUND_F5_CONTROL_AUDIT",
            "source_bound_execution_admitted":true,"target_complete":mathematics["target_complete"],"mathematics":mathematics}),
        Err(error) => {
            json!({"status":"TERMINAL_F5_CONTROL_NOT_ADMITTED","source_bound_execution_admitted":false,
            "target_complete":false,"error":error,"terminal":terminal,"complete_online_diagnostic_ns":null})
        }
    };
    let gate = check_capsule(&capsule);
    require(
        gate.is_ok() && before == sha256(&read(&own, 128 * 1024 * 1024)?),
        "F5 audit sources changed during verification",
    )?;
    receipt["schema_version"] = json!(1);
    receipt["registration_sha256"] = json!(canonical_sha(&registration)?);
    receipt["checker_binary_sha256"] = json!(before);
    receipt["independent_external_audit_wall_ns"] = json!(started.elapsed().as_nanos().to_string());
    receipt["audit_timing_policy"] = json!(
        "independent external replay outside original producer interval; no headline online claim"
    );
    receipt["native_children_executed_by_auditor"] = json!(0);
    receipt["fresh_paired_qualification"] = json!(false);
    receipt["headline_eligible"] = json!(false);
    receipt["promotion_eligible"] = json!(false);
    receipt["online_speedup"] = Value::Null;
    receipt["full_goal_complete"] = json!(false);
    receipt["ordinary_queries_executed"] = json!(0);
    save(out, &receipt)?;
    serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn cfg() -> Control {
        Control {
            schema_version: 1,
            question: "native-f5-known-input-transport-and-recovery-control-v1".into(),
            target: [52411, 72106],
            algorithm_seed: 2026100310,
            max_queries: 8,
            controller_timeout_ms: 900000,
        }
    }
    fn prep() -> Value {
        serde_json::from_str(include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/f5-preparation.json")).unwrap()
    }
    fn temp(label: &str) -> std::path::PathBuf {
        static SERIAL: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "native-f5-controller-{label}-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        path
    }
    #[test]
    fn validation_capsules_cannot_dispatch_or_admit_scientific_jobs() {
        require_executable_registration(&json!({"validation_only":false})).unwrap();
        for document in [
            json!({"validation_only":true}),
            json!({"validation_only":null}),
            json!({}),
        ] {
            assert!(require_executable_registration(&document).is_err());
        }
    }
    #[test]
    fn strict_control_rejects_other_targets_budgets_questions_and_unknown_fields() {
        cfg().validate().unwrap();
        for (key, value) in [
            ("target", json!([1, 2])),
            ("max_queries", json!(9)),
            ("max_queries", json!(0)),
            ("controller_timeout_ms", json!(900001)),
            ("question", json!("fresh-paired")),
            ("schema_version", json!(2)),
        ] {
            let mut document = json!(cfg());
            document[key] = value;
            let candidate: Control = serde_json::from_value(document).unwrap();
            assert!(candidate.validate().is_err());
        }
        let mut document = json!(cfg());
        document["known_scalar"] = json!(24886);
        assert!(serde_json::from_value::<Control>(document).is_err());
    }
    #[test]
    fn mathematical_stdin_contains_only_accepted_geometry_logs_point_and_limits() {
        let document = prep();
        let job = mathematical_job(&document, &cfg()).unwrap();
        let expected: Value = serde_json::from_str(include_str!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json")).unwrap();
        assert_eq!(job["prepared"], expected);
        assert_eq!(job["public_targets"], json!([["52411", "72106"]]));
        assert_eq!(job["target_seeds"], json!([]));
        assert_eq!(job["algorithm_seed"], 2026100310u64);
        assert_eq!(job["prepared"]["factor_base"].as_array().unwrap().len(), 63);
        assert_eq!(job["prepared"]["columns"].as_array().unwrap().len(), 29);
        for key in [
            "certificate",
            "attempts",
            "recovered",
            "scalar",
            "solutions",
            "provenance",
        ] {
            assert!(job.get(key).is_none());
            assert!(job["prepared"].get(key).is_none());
        }
        let mut changed = document;
        changed["record"]["column_logs"][0]["log"] = json!(0);
        assert!(mathematical_job(&changed, &cfg()).is_err());
    }
    #[test]
    fn empty_cargo_home_does_not_authorize_ancestor_build_configuration() {
        let p = temp("ancestor");
        let source = p.join("source");
        fs::create_dir(&source).unwrap();
        fs::create_dir(p.join(".cargo")).unwrap();
        fs::write(
            p.join(".cargo/config.toml"),
            b"[build]\nrustflags = [\"--cfg\", \"changed\"]\n",
        )
        .unwrap();
        assert!(reject_ancestor_cargo_configuration(&source).is_err());
        fs::remove_dir_all(p).unwrap();
    }
    #[test]
    fn simulated_transport_rejects_modified_environment_input_output_claim_and_drain() {
        // Data fixture only: neither this placeholder worker nor any solver runs.
        // Public audit separately requires the actual frozen checker/source tree.
        let p = temp("transport-fixture");
        let execution = p.join("execution");
        fs::create_dir(&execution).unwrap();
        let execution = execution.canonicalize().unwrap();
        let cfg = cfg();
        let worker = p.join("immutable/bin").join(WORKER);
        let job = serde_json::to_vec(&mathematical_job(&prep(), &cfg).unwrap()).unwrap();
        native::create(&p.join("job.json"), &job).unwrap();
        let raw = include_bytes!("../../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-v3-control-v1/controls/pipeline.stdout");
        let old: Value = serde_json::from_slice(raw).unwrap();
        native::create(&execution.join("worker.stdout"), raw).unwrap();
        native::create(&execution.join("worker.stderr"), b"").unwrap();
        native::create(&execution.join("worker-pids"), b"start 42\ndone 42\n").unwrap();
        let registration = json!({"test_fixture":true,"worker_sha256":"placeholder-never-executed", "worker_build_identity":old["generic_build"]});
        save(
            &p.join("consumed.json"),
            &json!({"registration_sha256":canonical_sha(&registration).unwrap(),
            "status":"consumed-before-launch","execution":execution}),
        )
        .unwrap();
        let receipt = json!({"argv":[worker],"cwd":execution,"pid":42,"deadline_ms":900000,
            "environment":environment().into_iter().collect::<std::collections::BTreeMap<_,_>>(),
            "executable_sha256_before":"placeholder-never-executed","executable_sha256_after":"placeholder-never-executed",
            "stdin_sha256":sha256(&job),"stdin_bytes":job.len(),"stdin_written_bytes":job.len(),"stdin_write_error":null,
            "stdout_sha256":sha256(raw),"stdout_bytes":raw.len(),"stderr_sha256":sha256(b""),"stderr_bytes":0,
            "exit_code":0,"timed_out":false,"output_limit":false,"process_group_drain_requested":true,"process_group_drain_confirmed":true});
        save(&execution.join("worker.receipt.json"), &receipt).unwrap();
        let terminal = json!({"schema_version":1,"registration_sha256":canonical_sha(&registration).unwrap(),
            "source_gate_passed":true,"child_drain_passed":true,"worker":{"exit_code":0,"timed_out":false,"receipt":receipt}});
        assert_eq!(
            verify_transport(&registration, &cfg, &p, &execution, &terminal).unwrap(),
            old
        );
        // The subsequent mathematics gate rejects this old seed under the new config.
        assert!(f5_target::verify(&prep(), cfg.algorithm_seed, cfg.max_queries, &old).is_err());
        for (pointer, value) in [
            ("/environment/RAYON_NUM_THREADS", json!("2")),
            ("/stdin_written_bytes", json!(0)),
            ("/stdin_write_error", json!("broken pipe")),
            ("/stdout_sha256", json!("wrong")),
            ("/process_group_drain_confirmed", json!(false)),
            ("/argv/0", json!("/unregistered/worker")),
            ("/exit_code", json!(2)),
        ] {
            let mut changed = receipt.clone();
            *changed.pointer_mut(pointer).unwrap() = value;
            fs::write(
                execution.join("worker.receipt.json"),
                serde_json::to_vec(&changed).unwrap(),
            )
            .unwrap();
            let mut forged = terminal.clone();
            forged["worker"]["receipt"] = changed;
            assert!(
                verify_transport(&registration, &cfg, &p, &execution, &forged).is_err(),
                "accepted {pointer}"
            );
        }
        fs::write(
            execution.join("worker.receipt.json"),
            serde_json::to_vec(&receipt).unwrap(),
        )
        .unwrap();
        fs::write(p.join("consumed.json"), b"{}").unwrap();
        assert!(verify_transport(&registration, &cfg, &p, &execution, &terminal).is_err());
        fs::remove_dir_all(p).unwrap();
    }
}
