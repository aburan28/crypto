//! Source-bound native development registration and independent execution audit.
use super::{json as records, oracle::Curve, prepared_sat, sat_query_law::QueryLaw, sat_source};
use crypto_lib::cryptanalysis::prepared_sat_control::{
    canonical_sha, sha256, ControlConfig, PHASES, PREPARATION_SHA256, TARGET,
};
use serde_json::{json, Value};
use std::{
    fs,
    path::Path,
    process::{Command, Stdio},
    time::Instant,
};
#[path = "../prepared_sat_worker/native.rs"]
#[allow(dead_code)] // The common transport is also compiled into the producer.
pub(super) mod native;
use native::{load, read, require, save};
const PREP:&str="research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/sat-preparation.json";
const ASSETS:&str="research/ic_candidate_tournament_20260915/goal_20260924/static-sat-runtime-v3/native-inputs-macos-arm64";
fn native_step(
    program: &Path,
    args: &[String],
    cwd: &Path,
    log: &Path,
    environment: &[(String, String)],
) -> Result<Value, String> {
    let out = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(log)
        .map_err(|e| e.to_string())?;
    let err = out.try_clone().map_err(|e| e.to_string())?;
    let before = sha256(&read(program, 128 * 1024 * 1024)?);
    let started = Instant::now();
    let status = Command::new(program)
        .args(args)
        .current_dir(cwd)
        .env_clear()
        .envs(environment.iter().cloned())
        .stdin(Stdio::null())
        .stdout(out)
        .stderr(err)
        .status()
        .map_err(|e| e.to_string())?;
    let after = sha256(&read(program, 128 * 1024 * 1024)?);
    let receipt = json!({"program":program,"argv":args,"environment":environment,"cwd":cwd,"exit_code":status.code(),
        "wall_ns":started.elapsed().as_nanos().to_string(),"program_sha256_before":before,"program_sha256_after":after,
        "log_sha256":sha256(&read(log,16*1024*1024)?)});
    save(&log.with_extension("receipt.json"), &receipt)?;
    require(
        status.success() && before == after,
        "native build step failed; original log retained",
    )?;
    Ok(receipt)
}
/// Freeze and build a controller and auditor from their retained sources and offline vendor tree.
/// This launches no exporter, SAT solver or scientific target computation.
pub fn freeze(
    root: &Path,
    out: &Path,
    config_file: &Path,
    cargo: &Path,
    rustc: &Path,
) -> Result<String, String> {
    native::enforce_hardware()?;
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let cargo = cargo.canonicalize().map_err(|e| e.to_string())?;
    let rustc = rustc.canonicalize().map_err(|e| e.to_string())?;
    let config: ControlConfig =
        serde_json::from_value(load(config_file)?).map_err(|e| e.to_string())?;
    config.validate()?;
    let prep = load(&root.join(PREP))?;
    require(
        canonical_sha(&prep)? == PREPARATION_SHA256,
        "preparation seal differs",
    )?;
    let preparation_admission = prepared_sat::verify(&records::parse(&prep.to_string())?)?;
    let head = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(&root)
        .output()
        .map_err(|e| e.to_string())?;
    require(head.status.success(), "source checkout has no Git HEAD")?;
    let dirty = Command::new("git")
        .args([
            "status",
            "--porcelain",
            "--untracked-files=all",
            "--",
            "src",
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
        "commit the complete controller sources before freezing",
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
        "docs/ic/calibration.json",
        "docs/ecbench/schema.sql",
        "docs/curves/registry.json",
    ] {
        let path = source.join(name);
        fs::create_dir_all(path.parent().ok_or("missing source parent")?)
            .map_err(|e| e.to_string())?;
        native::create(&path, &read(&root.join(name), 16 * 1024 * 1024)?)?;
    }
    let native_assets = native::extract_assets(&root.join(ASSETS), &immutable.join("assets"))?;
    let receipts = immutable.join("build-receipts");
    fs::create_dir(&receipts).map_err(|e| e.to_string())?;
    let path = std::env::var("PATH").map_err(|e| e.to_string())?;
    let vendor_env = vec![("PATH".into(), path.clone()), ("LC_ALL".into(), "C".into())];
    native_step(
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
    // Build only from the frozen tree. Ambient Cargo configuration/home and flags are excluded.
    fs::create_dir(source.join(".cargo")).map_err(|e| e.to_string())?;
    native::create(&source.join(".cargo/config.toml"),b"[source.crates-io]\nreplace-with = \"vendored-sources\"\n[source.vendored-sources]\ndirectory = \"vendor\"\n")?;
    let build_home = out.join("build-home");
    fs::create_dir(&build_home).map_err(|e| e.to_string())?;
    let target = out.join("build-target");
    let build_env = vec![
        ("PATH".into(), path),
        ("LC_ALL".into(), "C".into()),
        ("RUSTC".into(), rustc.to_string_lossy().into_owned()),
        (
            "CARGO_HOME".into(),
            build_home.to_string_lossy().into_owned(),
        ),
        (
            "CARGO_TARGET_DIR".into(),
            target.to_string_lossy().into_owned(),
        ),
    ];
    native_step(
        &rustc,
        &["--version".into(), "--verbose".into()],
        &source,
        &receipts.join("rustc.log"),
        &build_env,
    )?;
    native_step(
        &cargo,
        &["--version".into(), "--verbose".into()],
        &source,
        &receipts.join("cargo.log"),
        &build_env,
    )?;
    let snapshot = native::inventory(&source)?;
    native_step(
        &cargo,
        &[
            "build".into(),
            "--locked".into(),
            "--offline".into(),
            "--release".into(),
            "--bin".into(),
            "prepared_sat_worker".into(),
            "--bin".into(),
            "icprog".into(),
        ],
        &source,
        &receipts.join("build.log"),
        &build_env,
    )?;
    native::check_tree(&source, &snapshot)?;
    let bin = immutable.join("bin");
    fs::create_dir(&bin).map_err(|e| e.to_string())?;
    for name in ["prepared_sat_worker", "icprog"] {
        native::create(
            &bin.join(name),
            &read(&target.join("release").join(name), 128 * 1024 * 1024)?,
        )?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(bin.join(name), fs::Permissions::from_mode(0o755))
                .map_err(|e| e.to_string())?;
        }
    }
    save(&out.join("config.json"), &json!(config))?;
    save(&out.join("preparation.json"), &prep)?;
    let registration = json!({"schema_version":1,"scope":"disclosed-n17-native-prepared-control","source_commit":String::from_utf8(head.stdout).map_err(|e|e.to_string())?.trim(),
        "hardware":{"os":"macos","architecture":"aarch64"},"immutable_files":native::inventory(&immutable)?,
        "worker_sha256":sha256(&read(&bin.join("prepared_sat_worker"),128*1024*1024)?),
        "auditor_sha256":sha256(&read(&bin.join("icprog"),128*1024*1024)?),
        "config.json":sha256(&read(&out.join("config.json"),16*1024*1024)?),"preparation.json":sha256(&read(&out.join("preparation.json"),16*1024*1024)?),
        "preparation_admission":preparation_admission,"native_assets":native_assets,
        "primary_question":"known-input source-bound transport and recovery control; no fresh-target speed claim",
        "ordinary_queries_executed":0,"fresh_paired_qualification":false,"headline_eligible":false,
        "accounting":"five exclusive target-dependent phases; independent external replay timed separately",
        "reference":null,"operation_calibration":null,"online_speedup":null,"status":"registered-not-dispatched"});
    save(&out.join("registration.json"), &registration)?;
    let seal = json!({"registration_sha256":canonical_sha(&registration)?});
    save(&out.join("seal.json"), &seal)?;
    native::check_capsule(&out)?;
    serde_json::to_string_pretty(&seal).map_err(|e| e.to_string())
}
pub fn execute(
    capsule: &Path,
    execution: &Path,
    expected_registration: &str,
) -> Result<String, String> {
    native::enforce_hardware()?;
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    let registration = native::check_capsule(&capsule)?;
    let registration_sha = canonical_sha(&registration)?;
    require(
        registration_sha == expected_registration,
        "dispatch registration is not the explicitly frozen seal",
    )?;
    let dispatcher = std::env::current_exe().map_err(|e| e.to_string())?;
    require(
        registration["auditor_sha256"] == sha256(&read(&dispatcher, 128 * 1024 * 1024)?),
        "dispatch must use the frozen controller/checker binary",
    )?;
    let prep = load(&capsule.join("preparation.json"))?;
    prepared_sat::verify(&records::parse(&prep.to_string())?)?;
    fs::create_dir(execution).map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    // Consumed before any scientific native launch, including crashes and interrupted execution.
    save(
        &capsule.join("consumed.json"),
        &json!({"registration_sha256":registration_sha,"execution":execution,"status":"consumed-before-launch"}),
    )?;
    let cfg = native::config(&capsule)?;
    let result = native::measured_child(
        &capsule.join("immutable/bin/prepared_sat_worker"),
        &[
            "run".into(),
            "--capsule".into(),
            capsule.to_string_lossy().into_owned(),
            "--execution".into(),
            execution.to_string_lossy().into_owned(),
        ],
        &execution,
        &execution.join("controller"),
        cfg.controller_timeout_ms,
        &execution.join("controller-pids"),
        false,
    );
    // Outer kill reaches inherited helpers; the durable ledger reaches children that escaped.
    let child_drain = native::drain_ledger(&execution.join("child-pids"));
    let source_gate = native::check_capsule(&capsule);
    let terminal = match result {
        Ok(child) => {
            json!({"schema_version":1,"registration_sha256":registration_sha,"controller":child,"source_gate_passed":source_gate.is_ok(),"source_gate_error":source_gate.err(),"child_drain_passed":child_drain.is_ok(),"child_drain_error":child_drain.err()})
        }
        Err(e) => {
            json!({"schema_version":1,"registration_sha256":registration_sha,"controller_transport_error":e,"source_gate_passed":source_gate.is_ok(),"source_gate_error":source_gate.err(),"child_drain_passed":child_drain.is_ok(),"child_drain_error":child_drain.err()})
        }
    };
    save(&execution.join("terminal.json"), &terminal)?;
    serde_json::to_string_pretty(&terminal).map_err(|e| e.to_string())
}
fn number(v: &Value) -> Result<u64, String> {
    v.as_u64()
        .or_else(|| v.as_str().and_then(|s| s.parse().ok()))
        .ok_or("invalid independent scalar".into())
}
fn power(mut v: u64, mut e: u64) -> u64 {
    let mut out = 1;
    while e != 0 {
        if e & 1 != 0 {
            out = out * v % 65587;
        }
        v = v * v % 65587;
        e >>= 1;
    }
    out
}

/// Arithmetic audit consumes outputs only; it can never dispatch native search.
pub fn audit(capsule: &Path, execution: &Path, out: &Path) -> Result<String, String> {
    let before = Instant::now();
    let registration = native::check_capsule(capsule)?;
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    let own_hash = sha256(&read(&own, 128 * 1024 * 1024)?);
    require(
        registration["auditor_sha256"] == own_hash,
        "audit must use the frozen independent checker",
    )?;
    let terminal = load(&execution.join("terminal.json"))?;
    require(
        terminal["registration_sha256"] == canonical_sha(&registration)?
            && terminal["source_gate_passed"] == true
            && terminal["child_drain_passed"] == true,
        "terminal source gate failed",
    )?;
    let claim = load(&capsule.join("consumed.json"))?;
    require(
        claim["registration_sha256"] == canonical_sha(&registration)?
            && claim["execution"] == json!(execution.canonicalize().map_err(|e| e.to_string())?),
        "audit execution differs from one-use claim",
    )?;
    let prep = load(&capsule.join("preparation.json"))?;
    let prepared_record = records::parse(&prep.to_string())?;
    let prep_admission = prepared_sat::verify(&prepared_record)?;
    let producer_path = execution.join("producer.json");
    if !producer_path.exists() {
        let report = json!({"status":"TERMINAL_INCOMPLETE_OR_INTERRUPTED","registration_sha256":canonical_sha(&registration)?,"terminal":terminal,
        "source_bound_execution_admitted":false,"online_wall_ns":null,"online_speedup":null,"promotion_eligible":false,"native_children_executed_by_auditor":0});
        save(out, &report)?;
        return serde_json::to_string_pretty(&report).map_err(|e| e.to_string());
    }
    require(
        terminal["controller"]["exit_code"] == 0 && terminal["controller"]["timed_out"] == false,
        "producer lacks a clean terminal controller",
    )?;
    require(
        terminal["controller"]["receipt"]["process_group_drain_confirmed"] == true,
        "controller descendants did not drain",
    )?;
    let report = load(&producer_path)?;
    require(
        !contains_mock_receipt(&report),
        "mock process receipts cannot admit scientific execution",
    )?;
    let cfg = native::config(capsule)?;
    let (recovered, audited, total) =
        verify_report(&prep, &cfg, &report, execution, &registration)?;
    let phases = &report["online_phases_ns"];
    native::check_capsule(capsule)?;
    require(
        sha256(&read(&own, 128 * 1024 * 1024)?) == own_hash,
        "checker changed during independent audit",
    )?;
    require(
        total
            <= terminal["controller"]["receipt"]["child_wall_ns"]
                .as_u64()
                .ok_or("controller timing missing")?,
        "online interval exceeds controller interval",
    )?;
    let admitted = json!({"schema_version":1,"status":"PASS_NATIVE_SOURCE_BOUND_CONTROL_AUDIT","registration_sha256":canonical_sha(&registration)?,
        "producer_sha256":sha256(&read(&producer_path,16*1024*1024)?),"auditor_binary_sha256":own_hash,"preparation_admission":prep_admission,
        "target_input":TARGET,"audited_attempts":audited,"recovered_scalar":recovered,"scalar_verified":recovered.is_some(),
        "source_bound_execution_admitted":true,"target_complete":recovered.is_some(),"online_attempt_wall_ns":total,"online_wall_ns":recovered.map(|_|total),
        "online_phases_ns":phases,"independent_replay_wall_ns":before.elapsed().as_nanos().to_string(),
        "independent_replay_inside_online_interval":false,"native_children_executed_by_auditor":0,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,"online_speedup":null});
    save(out, &admitted)?;
    serde_json::to_string_pretty(&admitted).map_err(|e| e.to_string())
}

fn contains_mock_receipt(value: &Value) -> bool {
    match value {
        Value::Object(map) => {
            map.get("test_fixture") == Some(&Value::Bool(true))
                || map.values().any(contains_mock_receipt)
        }
        Value::Array(values) => values.iter().any(contains_mock_receipt),
        _ => false,
    }
}

pub(super) fn verify_report(
    prep: &Value,
    cfg: &ControlConfig,
    report: &Value,
    execution: &Path,
    registration: &Value,
) -> Result<(Option<u64>, Vec<Value>, u64), String> {
    cfg.validate()?;
    let prepared_record = records::parse(&prep.to_string())?;
    require(
        report["schema_version"] == 1
            && report["method"] == "native-prepared-sat-control-v1"
            && report["ordinary_queries_executed"] == 0
            && report["mathematical_state_sha256"]
                == crypto_lib::cryptanalysis::prepared_sat_control::STATE_SHA256
            && report["config"] == json!(cfg)
            && report["target_input"] == json!(TARGET)
            && report["preparation_sha256"] == PREPARATION_SHA256,
        "producer input differs from frozen workload",
    )?;
    let input = prepared_record.at("certificate")?.at("inputs")?;
    let curve = Curve::new(input.at("fixture")?)?;
    let base = input
        .at("base")?
        .as_arr()
        .ok_or("missing geometric base")?
        .iter()
        .map(|p| curve.decode(p))
        .collect::<Result<Vec<_>, _>>()?;
    let (_, projection) = prepared_sat::projected_columns(&curve, &base)?;
    let pairs = prepared_sat::pair_sums(&curve, &base);
    let target = curve.decode(&records::parse(&json!(TARGET).to_string())?)?;
    let logs = prep["record"]["column_logs"]
        .as_array()
        .ok_or("missing logs")?
        .iter()
        .map(|v| number(&v["log"]))
        .collect::<Result<Vec<_>, _>>()?;
    let attempts = report["target_attempts"]
        .as_array()
        .filter(|a| !a.is_empty() && a.len() <= cfg.max_queries)
        .ok_or("invalid target chronology")?;
    let mut law = QueryLaw::new(cfg.descent_query_seed);
    let mut recovered = None;
    let mut audited = Vec::new();
    for (trial, row) in attempts.iter().enumerate() {
        require(
            recovered.is_none() && row["trial"] == trial,
            "attempt follows success or trial order differs",
        )?;
        let (a, b) = (law.scalar(), law.scalar());
        require(row["a"] == a && row["b"] == b, "target query law differs")?;
        let point = curve.add(
            curve.mul(Some(curve.g), a as u128),
            curve.mul(target, b as u128),
        );
        require(
            row["public_point"] == json!(point),
            "target query point differs",
        )?;
        let status = row["status"].as_str().ok_or("missing target outcome")?;
        let dir = execution.join(format!("query-{trial:02}"));
        if point.is_none() {
            require(status == "IDENTITY_QUERY", "identity query dispatched")?;
            continue;
        }
        require(
            load(&dir.join("query.json"))? == json!({"trial":trial,"point":point}),
            "native query input differs",
        )?;
        if status == "QUERY_TRANSPORT_FAILURE" {
            require(
                row["candidate_scalar"].is_null()
                    && row["witness_indices"].is_null()
                    && row["source_model_valid"] == false,
                "transport failure contributes scalar",
            )?;
            audited.push(json!({"trial":trial,"status":status,"verified_relation":false}));
            continue;
        }
        let p = point.ok_or("identity dispatch")?;
        let (manifest, anf, cnf, files) =
            native::validate_exports(&dir.join("instance"), [p.0, p.1])?;
        require(
            row["source_receipt"]["files"] == files,
            "source receipt differs from retained inputs",
        )?;
        let receipt = load(&dir.join("cms.receipt.json"))?;
        let exporter_receipt = load(&dir.join("exporter.receipt.json"))?;
        require(
            exporter_receipt["exit_code"] == 0
                && exporter_receipt["timed_out"] == false
                && exporter_receipt["deadline_ms"] == cfg.exporter_timeout_ms
                && receipt["deadline_ms"] == cfg.solver_timeout_ms
                && exporter_receipt["process_group_drain_confirmed"] == true
                && receipt["process_group_drain_confirmed"] == true
                && row["source_receipt"]["exporter"] == exporter_receipt,
            "native deadlines/export receipt differ",
        )?;
        for (role, pin, child) in [
            ("exporter", native::EXPORTER_SHA, &exporter_receipt),
            ("cms", native::CMS_SHA, &receipt),
        ] {
            let launch = load(&dir.join(format!("{role}.launch.json")))?;
            require(
                launch["role"] == role
                    && launch["binary_sha256_before_exec"] == pin
                    && launch["config_sha256"] == canonical_sha(&json!(cfg))?
                    && launch["query_sha256"] == canonical_sha(&load(&dir.join("query.json"))?)?
                    && launch["pid"] == child["pid"]
                    && launch["environment"] == json!({"LC_ALL":"C"})
                    && child["executable_sha256_before"] == registration["worker_sha256"]
                    && child["executable_sha256_after"] == registration["worker_sha256"],
                "native launch source/query identity differs",
            )?;
        }
        let stdout = read(&dir.join("cms.stdout"), 8 * 1024 * 1024)?;
        require(
            row["native"]["receipt"] == receipt
                && row["native"]["exit_code"] == receipt["exit_code"]
                && row["native"]["timed_out"] == receipt["timed_out"]
                && receipt["stdout_sha256"] == sha256(&stdout)
                && receipt["stderr_sha256"]
                    == sha256(&read(&dir.join("cms.stderr"), 8 * 1024 * 1024)?)
                && row["native"]["stdout"]
                    == std::str::from_utf8(&stdout).map_err(|e| e.to_string())?,
            "native output/receipt differs",
        )?;
        let stdout = std::str::from_utf8(&stdout).map_err(|e| e.to_string())?;
        let statuses = stdout
            .lines()
            .filter(|l| l.starts_with("s "))
            .collect::<Vec<_>>();
        let native_status = if receipt["timed_out"] == true {
            "TIMEOUT"
        } else {
            match (receipt["exit_code"].as_i64(), statuses.as_slice()) {
                (Some(10), ["s SATISFIABLE"]) => "SAT_MODEL",
                (Some(20), ["s UNSATISFIABLE"]) => "SOURCE_UNSAT",
                (Some(15), ["s INDETERMINATE"]) => "CONFLICT_BUDGET_INCONCLUSIVE",
                (Some(0), ["s UNKNOWN"]) => "UNKNOWN_INCONCLUSIVE",
                _ => "NATIVE_ERROR",
            }
        };
        if native_status == "SAT_MODEL" {
            let validated = sat_source::verify_native_model(&anf, &cnf, stdout);
            if status == "INVALID_OR_NONLIFTING_SOURCE_MODEL" {
                require(
                    row["source_model_valid"] == false
                        && row["witness_indices"].is_null()
                        && row["candidate_scalar"].is_null(),
                    "invalid model contributes a relation",
                )?;
                if let Ok(model) = validated {
                    let xs: [u64; 3] = std::array::from_fn(|i| {
                        (0..6).fold(0, |v, k| v | ((model[6 * i + k] as u64) << k))
                    });
                    let exists =
                        base.iter()
                            .filter(|p| p.is_some_and(|p| p.0 == xs[0]))
                            .any(|p0| {
                                base.iter()
                                    .filter(|p| p.is_some_and(|p| p.0 == xs[1]))
                                    .any(|p1| {
                                        base.iter()
                                            .filter(|p| p.is_some_and(|p| p.0 == xs[2]))
                                            .any(|p2| curve.add(curve.add(*p0, *p1), *p2) == point)
                                    })
                            });
                    require(!exists, "claimed nonlifting model has a witness")?;
                }
            } else {
                let model = validated?;
                require(
                    row["source_model_valid"] == true
                        && row["source_model_sha256"]
                            == sha256(
                                &model[..51].iter().map(|&b| u8::from(b)).collect::<Vec<_>>(),
                            )
                        && row["expanded_model_sha256"]
                            == sha256(&model.iter().map(|&b| u8::from(b)).collect::<Vec<_>>()),
                    "source model digest/validity differs",
                )?;
                require(
                    status == "VALID_POINT_WITNESS",
                    "validated model has an invalid producer status",
                )?;
                require(
                    model.len()
                        == manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
                            .as_u64()
                            .ok_or("missing model width")? as usize,
                    "model width differs",
                )?;
                let indices = row["witness_indices"]
                    .as_array()
                    .filter(|a| a.len() == 3)
                    .ok_or("missing witness indices")?
                    .iter()
                    .map(|i| {
                        i.as_u64()
                            .filter(|&i| i < 63)
                            .map(|i| i as usize)
                            .ok_or("bad witness index")
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let mut sum = None;
                let mut projected = 0;
                for (j, &i) in indices.iter().enumerate() {
                    let expected = (0..6).fold(0, |v, k| v | ((model[6 * j + k] as u64) << k));
                    require(
                        base[i].is_some_and(|p| p.0 == expected),
                        "witness is unrelated to model",
                    )?;
                    sum = curve.add(sum, base[i]);
                    if let Some((column, coefficient)) = projection[i] {
                        projected = (projected + coefficient * logs[column] % 65587) % 65587;
                    }
                }
                require(sum == point, "target relation does not readd")?;
                let scalar = (projected + 65587 - 2 * a % 65587) % 65587
                    * power(2 * b % 65587, 65585)
                    % 65587;
                require(
                    row["candidate_scalar"] == scalar
                        && curve.mul(Some(curve.g), scalar as u128) == target,
                    "recovered scalar fails independent replay",
                )?;
                recovered = Some(scalar);
            }
        } else {
            require(
                status == native_status
                    && row["source_model_valid"] == false
                    && row["witness_indices"].is_null()
                    && row["candidate_scalar"].is_null(),
                "failed outcome contributes relation or status differs",
            )?;
            if status == "SOURCE_UNSAT" {
                require(
                    !prepared_sat::has_three_sum(&curve, &base, &pairs, point),
                    "claimed negative has a group decomposition",
                )?;
            }
        }
        audited.push(json!({"trial":trial,"status":status,"source_sha256":files}));
    }
    require(
        recovered.is_some() || attempts.len() == cfg.max_queries,
        "producer stopped early without success",
    )?;
    require(
        report["recovered_scalar"] == json!(recovered)
            && report["scalar_verified"] == recovered.is_some()
            && report["status"]
                == if recovered.is_some() {
                    "complete"
                } else {
                    "incomplete"
                },
        "completion fields differ",
    )?;
    let phases = report["online_phases_ns"]
        .as_object()
        .filter(|p| p.len() == 5)
        .ok_or("online phase set differs")?;
    let total = PHASES
        .iter()
        .map(|name| {
            phases
                .get(*name)
                .and_then(Value::as_u64)
                .ok_or("online phase missing")
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .try_fold(0u64, |sum, v| sum.checked_add(v).ok_or("timing overflow"))?;
    require(
        total > 0
            && report["online_attempt_wall_ns"] == total
            && report["online_wall_ns"] == json!(recovered.map(|_| total)),
        "online timing interval differs",
    )?;
    require(
        report["headline_eligible"] == false
            && report["fresh_paired_qualification"] == false
            && report["promotion_eligible"] == false
            && report["online_speedup"].is_null(),
        "known-input control promoted",
    )?;
    Ok((recovered, audited, total))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_lib::cryptanalysis::prepared_sat_control::{
        self as producer, NativeOutput, OnlineClock, QueryBackend, QueryOutput,
    };
    struct Fixture {
        root: std::path::PathBuf,
        execution: std::path::PathBuf,
    }
    impl QueryBackend for Fixture {
        fn query(
            &mut self,
            trial: usize,
            point: [u64; 2],
            cfg: &ControlConfig,
            _: &mut OnlineClock,
        ) -> Result<QueryOutput, String> {
            let dir = self.execution.join(format!("query-{trial:02}"));
            fs::create_dir(&dir).unwrap();
            let query = json!({"trial":trial,"point":point});
            save(&dir.join("query.json"), &query).unwrap();
            if trial == 0 {
                return Err("deliberate transport control; no native search".into());
            }
            assert_eq!(point, [62577, 27783]);
            let retained=self.root.join("research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check");
            fs::create_dir(dir.join("instance")).unwrap();
            for name in [
                "manifest.json",
                "instance.anf",
                "instance.xor.cnf",
                "instance.magma",
            ] {
                native::create(
                    &dir.join("instance").join(name),
                    &read(&retained.join(name), 2 * 1024 * 1024).unwrap(),
                )
                .unwrap();
            }
            let source = load(&retained.join("result.json")).unwrap();
            let model = source["cnf_assignment"].as_array().unwrap();
            let stdout = format!(
                "s SATISFIABLE\nv {} 0\n",
                model
                    .iter()
                    .enumerate()
                    .map(|(i, v)| if v.as_bool().unwrap() {
                        (i as i32 + 1).to_string()
                    } else {
                        (-(i as i32 + 1)).to_string()
                    })
                    .collect::<Vec<_>>()
                    .join(" ")
            );
            let mut cms_receipt = Value::Null;
            let mut exporter_receipt = Value::Null;
            for (role, pin) in [("exporter", native::EXPORTER_SHA), ("cms", native::CMS_SHA)] {
                let content = if role == "cms" {
                    stdout.as_bytes()
                } else {
                    b"retained control fixture\n"
                };
                native::create(&dir.join(format!("{role}.stdout")), content).unwrap();
                native::create(&dir.join(format!("{role}.stderr")), b"").unwrap();
                let receipt = json!({"test_fixture":true,"exit_code":if role=="cms"{10}else{0},"timed_out":false,"pid":42,
                    "deadline_ms":if role=="cms"{cfg.solver_timeout_ms}else{cfg.exporter_timeout_ms},"process_group_drain_confirmed":true,
                    "executable_sha256_before":"mock-worker","executable_sha256_after":"mock-worker",
                    "stdout_sha256":sha256(content),"stderr_sha256":sha256(b"")});
                save(&dir.join(format!("{role}.receipt.json")), &receipt).unwrap();
                save(&dir.join(format!("{role}.launch.json")),&json!({"test_fixture":true,"role":role,"binary_sha256_before_exec":pin,
                    "pid":42,"environment":{"LC_ALL":"C"},"config_sha256":canonical_sha(&json!(cfg)).unwrap(),"query_sha256":canonical_sha(&query).unwrap()})).unwrap();
                if role == "cms" {
                    cms_receipt = receipt;
                } else {
                    exporter_receipt = receipt;
                }
            }
            let (manifest, anf, cnf, files) =
                native::validate_exports(&dir.join("instance"), point).unwrap();
            Ok(QueryOutput {
                native: NativeOutput {
                    exit_code: Some(10),
                    timed_out: false,
                    stdout,
                    receipt: cms_receipt,
                },
                manifest,
                anf,
                cnf,
                source_receipt: json!({"files":files,"files_unchanged_before_after":true,"exporter":exporter_receipt}),
            })
        }
        fn progress(&mut self, _: &Value) -> Result<(), String> {
            Ok(())
        }
    }
    #[test]
    fn independent_complete_flow_and_mutation_controls_without_search() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf();
        let execution =
            std::env::temp_dir().join(format!("sat-complete-audit-{}", std::process::id()));
        fs::create_dir(&execution).unwrap();
        let prep = load(&root.join(PREP)).unwrap();
        let cfg = ControlConfig {
            schema_version: 1,
            question: "native-prepared-known-input-control-v1".into(),
            target: TARGET,
            descent_query_seed: 2026100102,
            export_nonce: 2026100103,
            conflict_budget: 100_000,
            max_queries: 2,
            exporter_timeout_ms: 30_000,
            solver_timeout_ms: 60_000,
            controller_timeout_ms: 180_000,
        };
        let registration = json!({"worker_sha256":"mock-worker"});
        let mut fixture = Fixture {
            root,
            execution: execution.clone(),
        };
        let report = producer::solve(&prep, &cfg, &mut fixture).unwrap();
        assert!(contains_mock_receipt(&report));
        let (scalar, rows, _) =
            verify_report(&prep, &cfg, &report, &execution, &registration).unwrap();
        assert_eq!(scalar, Some(24886));
        assert_eq!(rows.len(), 2);
        // These controls call the arithmetic/receipt verifier, never the source-bound admission wrapper.
        for (path, value) in [
            ("scalar", json!(24887)),
            ("a", json!(0)),
            ("model-valid", json!(false)),
            ("model-hash", json!("wrong")),
            ("witness", json!([29, 51, 3])),
            ("time", json!(0)),
            ("promotion", json!(true)),
            ("status", json!("SOURCE_UNSAT")),
        ] {
            let mut bad = report.clone();
            match path {
                "scalar" => bad["recovered_scalar"] = value,
                "a" => bad["target_attempts"][1]["a"] = value,
                "model-valid" => bad["target_attempts"][1]["source_model_valid"] = value,
                "model-hash" => bad["target_attempts"][1]["expanded_model_sha256"] = value,
                "witness" => bad["target_attempts"][1]["witness_indices"] = value,
                "time" => bad["online_attempt_wall_ns"] = value,
                "promotion" => bad["headline_eligible"] = value,
                "status" => bad["target_attempts"][1]["status"] = value,
                _ => unreachable!(),
            }
            assert!(
                verify_report(&prep, &cfg, &bad, &execution, &registration).is_err(),
                "accepted {path}"
            );
        }
        let source = execution.join("query-01/instance/instance.anf");
        let original = read(&source, 2 * 1024 * 1024).unwrap();
        fs::write(&source, b"changed").unwrap();
        assert!(verify_report(&prep, &cfg, &report, &execution, &registration).is_err());
        fs::write(&source, original).unwrap();
        fs::remove_dir_all(execution).unwrap();
    }
}
