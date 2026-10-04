//! New target-free source freeze and one-use native controller; never old controls.
use super::{ordinary_preparation, sat_control::native};
#[path = "../prepared_ordinary_worker/contract.rs"]
#[allow(dead_code)] // Shared typed capsule/claim/query contract with the worker.
mod capsule;
use crypto_lib::cryptanalysis::{
    prepared_ordinary::Family,
    prepared_sat_control::{canonical_sha, sha256},
};
use native::{load, read, require, save};
use serde_json::{json, Value};
use std::{fs, path::Path, process::Command, time::Instant};
const ASSETS: &str = "research/ic_candidate_tournament_20260915/goal_20260924/static-sat-runtime-v3/native-inputs-macos-arm64";

fn reject_ancestor_config(source: &Path) -> Result<(), String> {
    for parent in source.ancestors().skip(1) {
        for name in [".cargo/config", ".cargo/config.toml"] {
            require(
                !parent.join(name).try_exists().map_err(|e| e.to_string())?,
                "unregistered ancestor Cargo configuration; use an isolated capsule parent",
            )?;
        }
    }
    Ok(())
}
fn copy_file(root: &Path, source: &Path, name: &str) -> Result<(), String> {
    let target = source.join(name);
    fs::create_dir_all(target.parent().ok_or("source file lacks parent")?)
        .map_err(|e| e.to_string())?;
    native::create(&target, &read(&root.join(name), 256 * 1024 * 1024)?)
}
fn executable_copy(from: &Path, to: &Path) -> Result<(), String> {
    native::create(to, &read(from, 128 * 1024 * 1024)?)?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(to, fs::Permissions::from_mode(0o755)).map_err(|e| e.to_string())?;
    }
    Ok(())
}
/// Builds only. The worker identity command executes no arithmetic or solver.
pub fn freeze(
    root: &Path,
    out: &Path,
    config_file: &Path,
    cargo: &Path,
    rustc: &Path,
    host_file: &Path,
    validation_only: bool,
) -> Result<String, String> {
    native::enforce_hardware()?;
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let cargo = cargo.canonicalize().map_err(|e| e.to_string())?;
    let rustc = rustc.canonicalize().map_err(|e| e.to_string())?;
    let cfg: capsule::Config =
        serde_json::from_slice(&read(config_file, 65536)?).map_err(|e| e.to_string())?;
    cfg.validate()?;
    let host = load(host_file)?;
    require(
        host["os"] == "macos"
            && host["architecture"] == "aarch64"
            && host["physical_hardware"] == true
            && host["calibrated_performance_environment"] == false,
        "ordinary host context must not claim calibrated CPU performance",
    )?;
    let head = Command::new("git")
        .args(["-c", "gc.auto=0", "rev-parse", "HEAD"])
        .current_dir(&root)
        .output()
        .map_err(|e| e.to_string())?;
    require(head.status.success(), "ordinary source checkout lacks HEAD")?;
    let dirty = Command::new("git")
        .args([
            "-c",
            "gc.auto=0",
            "status",
            "--porcelain",
            "--untracked-files=all",
            "--",
            "src",
            "examples",
            "benches",
            "tests",
            "docs",
            "Cargo.toml",
            "build.rs",
        ])
        .current_dir(&root)
        .output()
        .map_err(|e| e.to_string())?;
    require(
        dirty.status.success() && dirty.stdout.is_empty(),
        "commit every ordinary-worker/controller source and embedded document before freezing",
    )?;
    fs::create_dir(out).map_err(|e| e.to_string())?;
    let out = out.canonicalize().map_err(|e| e.to_string())?;
    let immutable = out.join("immutable");
    fs::create_dir(&immutable).map_err(|e| e.to_string())?;
    let source = immutable.join("source");
    fs::create_dir(&source).map_err(|e| e.to_string())?;
    for name in ["src", "examples", "benches", "tests", "docs"] {
        native::copy_tree(&root.join(name), &source.join(name))?;
    }
    for name in ["Cargo.toml", "Cargo.lock"] {
        copy_file(&root, &source, name)?;
    }
    if root
        .join("build.rs")
        .try_exists()
        .map_err(|e| e.to_string())?
    {
        copy_file(&root, &source, "build.rs")?;
    }
    let receipts = immutable.join("build-receipts");
    fs::create_dir(&receipts).map_err(|e| e.to_string())?;
    if cfg.plan.family == Family::Cryptominisat {
        let archive = immutable.join("native-archive");
        fs::create_dir(&archive).map_err(|e| e.to_string())?;
        for name in ["assets.tar.gz", "manifest.json"] {
            native::create(
                &archive.join(name),
                &read(&root.join(ASSETS).join(name), 8 * 1024 * 1024)?,
            )?;
        }
        let assets = native::extract_assets(&archive, &immutable.join("assets"))?;
        save(&receipts.join("native-assets.json"), &assets)?;
    }
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
    native::create(&source.join(".cargo/config.toml"),b"[source.crates-io]\nreplace-with = \"vendored-sources\"\n[source.vendored-sources]\ndirectory = \"vendor\"\n")?;
    reject_ancestor_config(&source)?;
    let snapshot = native::inventory(&source)?;
    let source_sha = canonical_sha(&snapshot)?;
    let cargo_sha = sha256(&read(&cargo, 128 * 1024 * 1024)?);
    let rustc_sha = sha256(&read(&rustc, 128 * 1024 * 1024)?);
    let build = json!({"schema_version":1,"scope":capsule::SCOPE,"source_manifest_sha256":source_sha,
        "cargo_sha256":cargo_sha,"rustc_sha256":rustc_sha,"profile":"release","features":"default",
        "hardware":{"os":"macos","architecture":"aarch64"},"algorithm_flags":"none",
        "native_archive_sha256":if cfg.plan.family==Family::Cryptominisat{Some(native::ARCHIVE_SHA)}else{None}});
    let build_sha = canonical_sha(&build)?;
    save(&receipts.join("build-identity.json"), &build)?;
    let home = out.join("build-home");
    fs::create_dir(&home).map_err(|e| e.to_string())?;
    let target = out.join("build-target");
    let env = vec![
        ("PATH".into(), path),
        ("LC_ALL".into(), "C".into()),
        ("RUSTC".into(), rustc.to_string_lossy().into_owned()),
        ("CARGO_HOME".into(), home.to_string_lossy().into_owned()),
        (
            "CARGO_TARGET_DIR".into(),
            target.to_string_lossy().into_owned(),
        ),
        ("CARGO_INCREMENTAL".into(), "0".into()),
        (
            "IC_ORDINARY_SOURCE_MANIFEST_SHA256".into(),
            source_sha.clone(),
        ),
        ("IC_ORDINARY_BUILD_SHA256".into(), build_sha.clone()),
    ];
    for (program, name) in [(&rustc, "rustc"), (&cargo, "cargo")] {
        super::sat_control::native_step(
            program,
            &["--version".into(), "--verbose".into()],
            &source,
            &receipts.join(format!("{name}.log")),
            &env,
        )?;
    }
    super::sat_control::native_step(
        &cargo,
        &[
            "build".into(),
            "--locked".into(),
            "--offline".into(),
            "--release".into(),
            "--bin".into(),
            capsule::WORKER.into(),
            "--bin".into(),
            "icprog".into(),
        ],
        &source,
        &receipts.join("build.log"),
        &env,
    )?;
    native::check_tree(&source, &snapshot)?;
    require(
        sha256(&read(&cargo, 128 * 1024 * 1024)?) == cargo_sha
            && sha256(&read(&rustc, 128 * 1024 * 1024)?) == rustc_sha,
        "ordinary compiler/Cargo changed during freeze",
    )?;
    let bin = immutable.join("bin");
    fs::create_dir(&bin).map_err(|e| e.to_string())?;
    for name in [capsule::WORKER, "icprog"] {
        executable_copy(&target.join("release").join(name), &bin.join(name))?;
    }
    let identity = json!({"schema_version":1,"source_manifest_sha256":source_sha,"build_sha256":build_sha,
        "target_os":"macos","target_arch":"aarch64"});
    let run_env = capsule::environment().into_iter().collect::<Vec<_>>();
    super::sat_control::native_step(
        &bin.join(capsule::WORKER),
        &["build-identity".into()],
        &out,
        &receipts.join("worker-identity.log"),
        &run_env,
    )?;
    capsule::require_identity(
        &serde_json::from_slice(&read(&receipts.join("worker-identity.log"), 65536)?)
            .map_err(|e| e.to_string())?,
        &identity,
    )?;
    save(&out.join("config.json"), &json!(cfg))?;
    save(&out.join("host-context.json"), &host)?;
    let registration = capsule::Registration {
        schema_version: 1,
        scope: capsule::SCOPE.into(),
        source_commit: String::from_utf8(head.stdout)
            .map_err(|e| e.to_string())?
            .trim()
            .into(),
        hardware: json!({"os":"macos","architecture":"aarch64"}),
        immutable_files: native::inventory(&immutable)?,
        source_manifest_sha256: source_sha,
        config_sha256: sha256(&read(&out.join("config.json"), 65536)?),
        host_context_sha256: sha256(&read(&out.join("host-context.json"), 65536)?),
        worker_sha256: sha256(&read(&bin.join(capsule::WORKER), 128 * 1024 * 1024)?),
        auditor_sha256: sha256(&read(&bin.join("icprog"), 128 * 1024 * 1024)?),
        worker_build_identity: identity,
        runtime_environment: capsule::environment(),
        validation_only,
    };
    save(&out.join("registration.json"), &json!(registration))?;
    let seal = json!({"registration_sha256":canonical_sha(&json!(registration))?});
    save(&out.join("seal.json"), &seal)?;
    let (record, _) = capsule::check_capsule(&out)?;
    verify_build(&out, &record)?;
    if cfg.plan.family == Family::Cryptominisat {
        verify_native_assets(&out)?;
    }
    serde_json::to_string_pretty(&seal).map_err(|e| e.to_string())
}
fn expected_seal(record: &capsule::Registration, hash: &str, expected: &str) -> Result<(), String> {
    require(
        !record.validation_only
            && hash == expected
            && expected.len() == 64
            && expected
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b)),
        "ordinary registration is validation-only or differs from external seal",
    )
}
fn frozen_self(capsule: &Path, record: &capsule::Registration) -> Result<String, String> {
    capsule::require_identity(&capsule::identity(), &record.worker_build_identity)?;
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    let hash = sha256(&read(&own, 128 * 1024 * 1024)?);
    require(
        own.canonicalize().map_err(|e| e.to_string())?
            == capsule
                .join("immutable/bin/icprog")
                .canonicalize()
                .map_err(|e| e.to_string())?
            && hash == record.auditor_sha256,
        "ordinary execution/audit requires the frozen controller",
    )?;
    Ok(hash)
}
fn worker_args(capsule: &Path, execution: &Path) -> Vec<String> {
    vec![
        "run".into(),
        "--capsule".into(),
        capsule.to_string_lossy().into_owned(),
        "--execution".into(),
        execution.to_string_lossy().into_owned(),
    ]
}
/// Durably consume before the sole launch; always retain transport/source/drain outcomes.
pub fn execute(capsule: &Path, execution: &Path, expected: &str) -> Result<String, String> {
    native::enforce_hardware()?;
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    let (record, hash) = capsule::check_capsule(&capsule)?;
    expected_seal(&record, &hash, expected)?;
    let own = frozen_self(&capsule, &record)?;
    let cfg = capsule::config(&capsule)?;
    fs::create_dir(execution).map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    save(
        &capsule.join("consumed.json"),
        &json!(capsule::Claim {
            schema_version: 1,
            scope: capsule::SCOPE.into(),
            registration_sha256: hash.clone(),
            execution: execution.clone(),
            worker_sha256: record.worker_sha256.clone()
        }),
    )?;
    let env = capsule::environment().into_iter().collect::<Vec<_>>();
    let args = worker_args(&capsule, &execution);
    let child = native::measured_child_request(native::ChildRequest {
        program: &capsule.join("immutable/bin").join(capsule::WORKER),
        args: &args,
        cwd: &execution,
        stem: &execution.join("worker"),
        deadline_ms: cfg.worker_timeout_ms,
        ledger: &execution.join("worker-pids"),
        helper: false,
        input: None,
        environment: &env,
    });
    let worker_drain = native::drain_ledger(&execution.join("worker-pids"));
    let role_drain = native::drain_ledger(&execution.join("child-pids"));
    let gate = capsule::check_capsule(&capsule);
    let evidence = native::inventory(&execution);
    let own_after = std::env::current_exe()
        .map_err(|e| e.to_string())
        .and_then(|p| read(&p, 128 * 1024 * 1024))
        .map(|b| sha256(&b));
    let unchanged = own_after.as_ref().is_ok_and(|hash| hash == &own);
    let mut terminal = json!({"schema_version":1,"scope":capsule::SCOPE,"registration_sha256":hash,
        "controller_sha256":own,"controller_unchanged":unchanged,"controller_read_error":own_after.err(),"source_gate_passed":gate.is_ok(),"source_gate_error":gate.err(),
        "worker_drain_passed":worker_drain.is_ok(),"worker_drain_error":worker_drain.err(),
        "role_drain_passed":role_drain.is_ok(),"role_drain_error":role_drain.err(),
        "execution_files":evidence.as_ref().ok(),"execution_inventory_error":evidence.as_ref().err(),
        "source_bound_execution_admitted":false,"full_goal_complete":false});
    match child {
        Ok(child) => {
            terminal["worker"] = json!({"exit_code":child.exit_code,"timed_out":child.timed_out,"receipt":child.receipt})
        }
        Err(error) => terminal["worker_transport_error"] = json!(error),
    }
    save(&execution.join("terminal.json"), &terminal)?;
    serde_json::to_string_pretty(&terminal).map_err(|e| e.to_string())
}

/// Data-only chronology of a killed panel. No solver, no runtime admission.
pub fn prefix(execution: &Path, cfg: &capsule::Config) -> Result<Value, String> {
    cfg.validate()?;
    let mut records = Vec::new();
    let mut pending = None;
    let mut names = std::collections::BTreeSet::new();
    for entry in fs::read_dir(execution).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let name = entry.file_name().to_string_lossy().into_owned();
        if name.starts_with("start-") || name.starts_with("progress-") {
            names.insert(name);
        }
    }
    let mut expected_names = std::collections::BTreeSet::new();
    for trial in 0..cfg.plan.planned_queries {
        let start = format!("start-{trial:03}.json");
        let done = format!("progress-{trial:03}.json");
        if !names.contains(&start) {
            break;
        }
        expected_names.insert(start.clone());
        let start_record = load(&execution.join(&start))?;
        let scalar = crypto_lib::cryptanalysis::koblitz_index_calculus::probe_scalar(
            cfg.plan.algorithm_seed,
            trial as u64,
            65587,
        );
        require(
            start_record == json!({"trial":trial,"scalar":scalar}),
            "durable ordinary start differs from frozen law",
        )?;
        if !names.contains(&done) {
            pending = Some(start_record);
            break;
        }
        expected_names.insert(done.clone());
        let record = load(&execution.join(done))?;
        require(
            record["attempt"]["trial"] == trial && record["attempt"]["scalar"] == scalar,
            "durable ordinary completion differs from start",
        )?;
        records.push(record);
    }
    require(
        names == expected_names,
        "ordinary prefix has extra, gapped or out-of-panel records",
    )?;
    Ok(
        json!({"completed_queries":records.len(),"records":records,"pending_start":pending,
        "panel_count_complete":expected_names.len()==cfg.plan.planned_queries*2,
        "source_bound_execution_admitted":false,"full_goal_complete":false}),
    )
}

/// Inspect a retained prefix and independently audit any completed mathematical report.
/// Actual source/model runtime admission still requires the complete receipt gate.
pub fn inspect(
    capsule: &Path,
    execution: &Path,
    expected: &str,
    out: &Path,
) -> Result<String, String> {
    let started = Instant::now();
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    let (record, hash) = capsule::check_capsule(&capsule)?;
    expected_seal(&record, &hash, expected)?;
    let own = frozen_self(&capsule, &record)?;
    capsule::claim(&capsule, &execution, &record, &hash)?;
    let cfg = capsule::config(&capsule)?;
    let retained = prefix(&execution, &cfg)?;
    let report = if execution
        .join("producer.json")
        .try_exists()
        .map_err(|e| e.to_string())?
    {
        Some(load(&execution.join("producer.json"))?)
    } else {
        None
    };
    let mathematics = report
        .as_ref()
        .map(|r| {
            require(
                r["plan"] == json!(cfg.plan)
                    && r["query_records"] == retained["records"]
                    && r["mathematical_input"]["attempts"]
                        == json!(retained["records"]
                            .as_array()
                            .ok_or("missing prefix records")?
                            .iter()
                            .map(|r| r["attempt"].clone())
                            .collect::<Vec<_>>()),
                "ordinary prefix and optional producer transcript differ",
            )?;
            ordinary_preparation::audit(&r["mathematical_input"])
        })
        .transpose()?;
    capsule::check_capsule(&capsule)?;
    require(
        sha256(&read(
            &std::env::current_exe().map_err(|e| e.to_string())?,
            128 * 1024 * 1024,
        )?) == own,
        "ordinary prefix checker changed during inspection",
    )?;
    let result = json!({"schema_version":1,"status":"ORDINARY_PREFIX_INSPECTED_RUNTIME_NOT_ADMITTED",
        "registration_sha256":hash,"prefix":retained,"mathematics":mathematics,
        "terminal":load(&execution.join("terminal.json"))?,"checker_sha256":own,
        "inspection_wall_ns":started.elapsed().as_nanos().to_string(),"native_children_executed_by_inspector":0,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,"headline_eligible":false,
        "promotion_eligible":false,"full_goal_complete":false,"online_wall_ns":null,"online_speedup":null});
    save(out, &result)?;
    serde_json::to_string_pretty(&result).map_err(|e| e.to_string())
}

fn integer(value: &Value) -> Result<u64, String> {
    value
        .as_u64()
        .ok_or("missing unsigned timing/counter".into())
}
fn pid(receipt: &Value) -> Result<u32, String> {
    receipt["pid"]
        .as_u64()
        .filter(|&v| v > 1 && v <= u32::MAX as u64)
        .map(|v| v as u32)
        .ok_or("invalid ordinary child PID".into())
}
fn verify_receipt(
    execution: &Path,
    stem: &Path,
    receipt: &Value,
    argv: &[String],
    program_hash: &str,
    deadline: u64,
    env: &Value,
) -> Result<u32, String> {
    let stdout = read(&stem.with_extension("stdout"), 8 * 1024 * 1024)?;
    let stderr = read(&stem.with_extension("stderr"), 8 * 1024 * 1024)?;
    require(
        receipt["argv"] == json!(argv)
            && receipt["cwd"] == json!(execution)
            && receipt["executable_sha256_before"] == program_hash
            && receipt["executable_sha256_after"] == program_hash
            && receipt["deadline_ms"] == deadline
            && receipt["environment"] == *env
            && receipt["stdout_sha256"] == sha256(&stdout)
            && receipt["stdout_bytes"] == stdout.len()
            && receipt["stderr_sha256"] == sha256(&stderr)
            && receipt["stderr_bytes"] == stderr.len()
            && receipt["output_limit"] == false
            && receipt["process_group_drain_requested"] == true
            && receipt["process_group_drain_confirmed"] == true,
        "ordinary exact invocation/output/drain receipt differs",
    )?;
    require(
        integer(&receipt["child_wall_ns"])? > 0,
        "ordinary child lacks observed interval",
    )?;
    pid(receipt)
}
fn verify_costs(report: &Value, count: usize) -> Result<u64, String> {
    let mut phases = std::collections::BTreeMap::<String, u64>::new();
    let mut windows = 0u64;
    let mut add = |clock: &Value| -> Result<(), String> {
        const PHASES: [&str; 14] = [
            "setup",
            "factor_base",
            "precompute",
            "queries",
            "pdp",
            "relation_check",
            "matrix_build",
            "relation_la",
            "target_query",
            "target_pdp",
            "target_relation_check",
            "target_descent",
            "recovery_check",
            "rho_solve",
        ];
        const ONLINE: [&str; 6] = [
            "target_query",
            "target_pdp",
            "target_relation_check",
            "target_descent",
            "recovery_check",
            "rho_solve",
        ];
        require(
            clock["schema_version"] == 1
                && clock.as_object().is_some_and(|o| o.len() == 5)
                && clock["phases_ns"]
                    .as_object()
                    .is_some_and(|p| p.len() == 14 && PHASES.iter().all(|k| p.contains_key(*k)))
                && clock["online_phases_ns"]
                    .as_object()
                    .is_some_and(|p| p.len() == 6 && ONLINE.iter().all(|k| p.contains_key(*k))),
            "preparation snapshot is incomplete or has undeclared phases",
        )?;
        require(
            clock["online_wall_ns"].is_null()
                && clock["online_phases_ns"]
                    .as_object()
                    .is_some_and(|p| p.values().all(Value::is_null)),
            "preparation unexpectedly entered online timing",
        )?;
        let mut sum = 0u64;
        for (name, value) in clock["phases_ns"]
            .as_object()
            .ok_or("missing preparation phases")?
        {
            if !value.is_null() {
                let value = integer(value)?;
                sum = sum.checked_add(value).ok_or("phase overflow")?;
                let total = phases.entry(name.clone()).or_default();
                *total = total.checked_add(value).ok_or("phase overflow")?;
            }
        }
        let observed = integer(&clock["observed_wall_ns"])?;
        require(sum == observed, "preparation phase window does not close")?;
        windows = windows.checked_add(observed).ok_or("window overflow")?;
        Ok(())
    };
    let initial = report["phase_windows"]
        .as_array()
        .filter(|p| p.len() == 2)
        .ok_or("preparation window list differs")?;
    require(
        initial[0]["stage"] == "reusable_initialization"
            && initial[1]["stage"] == "final_relation_la_and_group_replay",
        "preparation stage order differs",
    )?;
    add(&initial[0]["clock"])?;
    let records = report["query_records"]
        .as_array()
        .filter(|r| r.len() == count)
        .ok_or("preparation costs omit attempt")?;
    for record in records {
        require(
            record["costs"]["phases_ns"]["queries"].is_u64()
                && record["costs"]["phases_ns"]["pdp"].is_u64()
                && record["costs"]["phases_ns"]["relation_check"].is_u64(),
            "ordinary attempt lacks its executed phase costs",
        )?;
        add(&record["costs"])?;
    }
    add(&initial[1]["clock"])?;
    let controller = integer(&report["controller_setup_ns"])?;
    let total = phases.entry("setup".into()).or_default();
    *total = total.checked_add(controller).ok_or("controller overflow")?;
    let wall = integer(&report["preparation_wall_ns"])?;
    require(
        windows.checked_add(controller) == Some(wall)
            && report["exclusive_preparation_phases_ns"] == json!(phases)
            && report["online_wall_ns"].is_null()
            && report["online_speedup"].is_null()
            && report["source_bound_execution_admitted"] == false
            && report["full_goal_complete"] == false,
        "preparation aggregate, unknown online costs or producer admission flags differ",
    )?;
    Ok(wall)
}
fn verify_matrix(report: &Value, mathematics: &Value) -> Result<(), String> {
    let matrix = &report["relation_matrix"];
    require(
        matrix["modulus"] == "65587"
            && matrix["solver"] == "dense-gauss"
            && matrix["sparse_options"].is_null(),
        "ordinary final relation-LA contract differs",
    )?;
    let parse = |v: &Value| {
        v.as_str()
            .and_then(|s| s.parse::<u64>().ok())
            .filter(|&v| v < 65587)
            .ok_or("invalid relation residue".to_string())
    };
    let columns = matrix["column_points"]
        .as_array()
        .filter(|c| c.len() == 29)
        .ok_or("ordinary matrix columns differ")?
        .iter()
        .map(|p| {
            let p = p
                .as_array()
                .filter(|p| p.len() == 2)
                .ok_or("invalid matrix point")?;
            p.iter()
                .map(|v| {
                    v.as_str()
                        .and_then(|s| s.parse::<u64>().ok())
                        .filter(|&v| v < 1 << 17)
                        .ok_or("invalid matrix coordinate".to_string())
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    require(
        json!(columns) == mathematics["column_points"],
        "producer matrix column order differs from independent quotient",
    )?;
    let rows = matrix["rows"]
        .as_array()
        .filter(|r| r.len() <= 512)
        .ok_or("ordinary matrix rows exceed panel")?
        .iter()
        .map(|row| {
            let mut coefficients = vec![0u64; 29];
            let mut previous = None;
            for entry in row["entries"]
                .as_array()
                .ok_or("missing matrix row entries")?
            {
                let entry = entry
                    .as_array()
                    .filter(|e| e.len() == 2)
                    .ok_or("invalid matrix entry")?;
                let j = entry[0]
                    .as_u64()
                    .filter(|&v| v < 29)
                    .ok_or("matrix index outside columns")? as usize;
                require(
                    previous.is_none_or(|old| old < j),
                    "matrix entries repeated or unsorted",
                )?;
                previous = Some(j);
                coefficients[j] = parse(&entry[1])?;
                require(coefficients[j] != 0, "explicit zero sparse entry")?;
            }
            Ok(json!({"coefficients":coefficients,"rhs":parse(&row["rhs"])?}))
        })
        .collect::<Result<Vec<Value>, String>>()?;
    require(
        rows.len() as u64 == integer(&mathematics["accepted_rows"])?
            && canonical_sha(&json!(rows))? == mathematics["rows_sha256"]
            && report["log_report"]["relations"] == mathematics["accepted_rows"]
            && report["log_report"]["duplicate_relations"] == mathematics["duplicate_relations"],
        "producer relation matrix/counts differ from independently reconstructed rows",
    )
}
fn verify_cms_rows(
    root: &Path,
    execution: &Path,
    cfg: &capsule::Config,
    record: &capsule::Registration,
    report: &Value,
) -> Result<usize, String> {
    let curve = super::oracle::Curve::new(&super::json::parse(
        &report["mathematical_input"]["fixture"].to_string(),
    )?)?;
    let base = report["mathematical_input"]["geometric_base"]
        .as_array()
        .ok_or("missing audited base")?
        .iter()
        .map(|p| curve.decode(&super::json::parse(&p.to_string())?))
        .collect::<Result<Vec<_>, String>>()?;
    let records = report["query_records"]
        .as_array()
        .ok_or("missing ordinary source records")?;
    let mut ledger = String::new();
    for (trial, row) in records.iter().enumerate() {
        let dir = execution.join(format!("query-{trial:03}"));
        let q = capsule::query(&dir, cfg)?;
        let p = curve
            .mul(Some(curve.g), q.scalar as u128)
            .ok_or("ordinary query is identity")?;
        require(
            q.point == [p.0, p.1],
            "native query differs from independent group arithmetic",
        )?;
        let front = &row["evidence"]["frontend"];
        let outcome = row["attempt"]["outcome"]
            .as_str()
            .ok_or("missing ordinary outcome")?;
        require(
            outcome != "transport_failure",
            "source transport failure retained but not runtime-admitted by this complete gate",
        )?;
        let (manifest, anf, cnf, files) = native::validate_exports(&dir.join("instance"), q.point)?;
        require(
            front["public_point"] == json!(q.point)
                && front["manifest"] == manifest
                && front["anf_sha256"] == sha256(anf.as_bytes())
                && front["anf_bytes"] == anf.len()
                && front["cnf_sha256"] == sha256(cnf.as_bytes())
                && front["cnf_bytes"] == cnf.len()
                && front["source_receipt"]["files"] == files
                && front["source_receipt"]["files_unchanged_before_after"] == true,
            "ordinary original source/hash receipt differs",
        )?;
        let mut cms = None;
        for (role, deadline) in [
            ("exporter", cfg.exporter_timeout_ms),
            ("cms", cfg.solver_timeout_ms),
        ] {
            let child = load(&dir.join(format!("{role}.receipt.json")))?;
            let argv = vec![
                root.join("immutable/bin")
                    .join(capsule::WORKER)
                    .to_string_lossy()
                    .into_owned(),
                "child".into(),
                "--capsule".into(),
                root.to_string_lossy().into_owned(),
                "--attempt".into(),
                dir.to_string_lossy().into_owned(),
                "--role".into(),
                role.into(),
            ];
            let child_pid = verify_receipt(
                &dir,
                &dir.join(role),
                &child,
                &argv,
                &record.worker_sha256,
                deadline,
                &json!({"LC_ALL":"C"}),
            )?;
            ledger.push_str(&format!("start {child_pid}\ndone {child_pid}\n"));
            let (program, args, pin) = capsule::role(root, cfg, &q, role)?;
            let launch = load(&dir.join(format!("{role}.launch.json")))?;
            require(
                launch
                    == json!({"role":role,"program":program,"argv":args,"binary_sha256_before_exec":pin,
                "config_sha256":canonical_sha(&json!(cfg))?,"query_sha256":canonical_sha(&json!(q))?,
                "environment":{"LC_ALL":"C"},"cwd":dir,"pid":child_pid}),
                "ordinary native exec source/argv receipt differs",
            )?;
            if role == "exporter" {
                require(
                    child["exit_code"] == 0
                        && child["timed_out"] == false
                        && front["source_receipt"]["exporter"] == child,
                    "ordinary exporter failed or receipt differs",
                )?;
            } else {
                cms = Some(child);
            }
        }
        let cms = cms.ok_or("missing native CMS receipt")?;
        let stdout = read(&dir.join("cms.stdout"), 8 * 1024 * 1024)?;
        require(
            front["native"]["receipt"] == cms
                && front["native"]["exit_code"] == cms["exit_code"]
                && front["native"]["timed_out"] == cms["timed_out"]
                && front["native"]["stdout_sha256"] == sha256(&stdout)
                && front["native"]["stdout_bytes"] == stdout.len(),
            "ordinary original native output differs",
        )?;
        let stdout = std::str::from_utf8(&stdout).map_err(|e| e.to_string())?;
        let statuses = stdout
            .lines()
            .filter(|l| l.starts_with("s "))
            .collect::<Vec<_>>();
        let status = if cms["timed_out"] == true {
            "TIMEOUT"
        } else {
            match (cms["exit_code"].as_i64(), statuses.as_slice()) {
                (Some(10), ["s SATISFIABLE"]) => "SAT_MODEL",
                (Some(20), ["s UNSATISFIABLE"]) => "SOURCE_UNSAT",
                (Some(15), ["s INDETERMINATE"]) => "CONFLICT_BUDGET_INCONCLUSIVE",
                (Some(0), ["s UNKNOWN"]) => "UNKNOWN_INCONCLUSIVE",
                _ => "NATIVE_ERROR",
            }
        };
        require(
            front["native_status"] == status,
            "ordinary producer native classification differs",
        )?;
        match status {
            "SAT_MODEL" => {
                let checked = super::sat_source::verify_native_model(&anf, &cnf, stdout);
                if outcome == "witness" {
                    let model = checked?;
                    let indices: [usize; 3] =
                        serde_json::from_value(row["attempt"]["indices"].clone())
                            .map_err(|e| e.to_string())?;
                    for (i, &index) in indices.iter().enumerate() {
                        let x = (0..6).fold(0, |v, b| v | ((model[i * 6 + b] as u64) << b));
                        require(
                            base.get(index).is_some_and(|p| p.is_some_and(|p| p.0 == x)),
                            "ordinary witness differs from native source assignment",
                        )?;
                    }
                } else {
                    require(
                        outcome == "invalid_model" && row["attempt"]["indices"].is_null(),
                        "SAT model has incompatible outcome",
                    )?;
                    if let Ok(model) = checked {
                        let xs: [u64; 3] = std::array::from_fn(|i| {
                            (0..6).fold(0, |v, b| v | ((model[i * 6 + b] as u64) << b))
                        });
                        let lifts =
                            base.iter()
                                .filter(|p| p.is_some_and(|p| p.0 == xs[0]))
                                .any(|a| {
                                    base.iter().filter(|p| p.is_some_and(|p| p.0 == xs[1])).any(
                                        |b| {
                                            base.iter()
                                                .filter(|p| p.is_some_and(|p| p.0 == xs[2]))
                                                .any(|c| {
                                                    curve.add(curve.add(*a, *b), *c) == Some(p)
                                                })
                                        },
                                    )
                                });
                        require(!lifts, "claimed invalid source model actually lifts")?;
                    }
                }
            }
            "SOURCE_UNSAT" => require(
                outcome == "proved_unsat",
                "native negative classified differently",
            )?,
            "TIMEOUT" => require(
                outcome == "timeout",
                "native timeout classified differently",
            )?,
            "CONFLICT_BUDGET_INCONCLUSIVE" | "UNKNOWN_INCONCLUSIVE" => require(
                outcome == "incomplete",
                "native budget stop classified differently",
            )?,
            _ => {
                return Err("native error retained but not admitted by complete source gate".into())
            }
        }
    }
    require(
        std::str::from_utf8(&read(&execution.join("child-pids"), 65536)?)
            .map_err(|e| e.to_string())?
            == ledger,
        "ordinary native role ledger has missing or extra launches",
    )?;
    Ok(records.len() * 2)
}
fn verify_native_assets(root: &Path) -> Result<(), String> {
    let temp = std::env::temp_dir().join(format!("ordinary-assets-audit-{}", std::process::id()));
    // The accepted fixed archive is decoded as data; no native file is executed.
    let expected = native::extract_assets(&root.join("immutable/native-archive"), &temp)?;
    let result = (|| {
        require(
            expected == load(&root.join("immutable/build-receipts/native-assets.json"))?
                && expected["files"] == native::inventory(&root.join("immutable/assets"))?,
            "ordinary native files differ from independently decoded accepted archive",
        )
    })();
    let cleanup = fs::remove_dir_all(&temp).map_err(|e| e.to_string());
    result?;
    cleanup
}
fn verify_execution_inventory(execution: &Path, terminal: &Value) -> Result<(), String> {
    let mut files = native::inventory(execution)?;
    files
        .as_object_mut()
        .ok_or("invalid execution inventory")?
        .remove("terminal.json");
    require(
        terminal["execution_inventory_error"].is_null()
            && terminal["execution_files"].is_object()
            && files == terminal["execution_files"],
        "ordinary retained execution files changed after drain",
    )
}
fn verify_build(root: &Path, record: &capsule::Registration) -> Result<(), String> {
    let dir = root.join("immutable/build-receipts");
    let build = load(&dir.join("build-identity.json"))?;
    require(
        canonical_sha(&build)? == record.worker_build_identity["build_sha256"]
            && build["source_manifest_sha256"] == record.source_manifest_sha256
            && build["scope"] == capsule::SCOPE
            && build["profile"] == "release"
            && build["features"] == "default"
            && build["algorithm_flags"] == "none"
            && build["hardware"] == record.hardware,
        "ordinary compiled build identity/flags differ",
    )?;
    let cfg = capsule::config(root)?;
    require(
        build["native_archive_sha256"]
            == if cfg.plan.family == Family::Cryptominisat {
                json!(native::ARCHIVE_SHA)
            } else {
                Value::Null
            },
        "ordinary compiled native asset binding differs",
    )?;
    let build_env = load(&dir.join("build.receipt.json"))?["environment"].clone();
    for name in ["vendor", "rustc", "cargo", "build", "worker-identity"] {
        let receipt = load(&dir.join(format!("{name}.receipt.json")))?;
        let bytes = read(&dir.join(format!("{name}.log")), 16 * 1024 * 1024)?;
        let pin = if name == "worker-identity" {
            json!(record.worker_sha256)
        } else if name == "rustc" {
            build["rustc_sha256"].clone()
        } else {
            build["cargo_sha256"].clone()
        };
        require(
            receipt["exit_code"] == 0
                && receipt["program_sha256_before"] == pin
                && receipt["program_sha256_after"] == pin
                && receipt["log_sha256"] == sha256(&bytes),
            "ordinary build tool/log receipt differs",
        )?;
        if name == "worker-identity" {
            require(
                serde_json::from_slice::<Value>(&bytes).map_err(|e| e.to_string())?
                    == record.worker_build_identity
                    && receipt["argv"] == json!(["build-identity"])
                    && receipt["environment"]
                        == json!(capsule::environment().into_iter().collect::<Vec<_>>()),
                "ordinary observed compiled worker identity differs",
            )?;
        }
        if name == "rustc" || name == "cargo" {
            require(
                receipt["argv"] == json!(["--version", "--verbose"])
                    && receipt["environment"] == build_env,
                "ordinary build tool version arguments/environment differ",
            )?;
        }
        if name == "vendor" {
            require(
                receipt["argv"]
                    == json!([
                        "vendor",
                        "--locked",
                        "--offline",
                        root.join("immutable/source/vendor")
                    ]),
                "ordinary vendor arguments differ",
            )?;
            let env: Vec<(String, String)> = serde_json::from_value(receipt["environment"].clone())
                .map_err(|e| e.to_string())?;
            require(
                env.len() == 2 && env[0].0 == "PATH" && env[1] == ("LC_ALL".into(), "C".into()),
                "ordinary vendor environment differs",
            )?;
        }
        if name == "build" {
            require(
                receipt["argv"]
                    == json!([
                        "build",
                        "--locked",
                        "--offline",
                        "--release",
                        "--bin",
                        capsule::WORKER,
                        "--bin",
                        "icprog"
                    ]),
                "ordinary build arguments differ",
            )?;
            let pairs: Vec<(String, String)> =
                serde_json::from_value(receipt["environment"].clone())
                    .map_err(|e| e.to_string())?;
            let env: std::collections::BTreeMap<String, String> = pairs.iter().cloned().collect();
            require(
                pairs.len() == 8
                    && env.len() == 8
                    && env.get("IC_ORDINARY_SOURCE_MANIFEST_SHA256")
                        == Some(&record.source_manifest_sha256)
                    && env.get("IC_ORDINARY_BUILD_SHA256").map(String::as_str)
                        == record.worker_build_identity["build_sha256"].as_str()
                    && env.get("CARGO_INCREMENTAL").map(String::as_str) == Some("0")
                    && env.get("LC_ALL").map(String::as_str) == Some("C")
                    && ["PATH", "RUSTC", "CARGO_HOME", "CARGO_TARGET_DIR"]
                        .iter()
                        .all(|k| env.contains_key(*k)),
                "ordinary build environment/flags differ",
            )?;
        }
    }
    Ok(())
}

/// Complete frozen data/source audit; never launches a producer or solver.
pub fn audit(root: &Path, execution: &Path, expected: &str, out: &Path) -> Result<String, String> {
    let started = Instant::now();
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    let own_sha = sha256(&read(&own, 128 * 1024 * 1024)?);
    let result = (|| -> Result<Value, String> {
        let root = root.canonicalize().map_err(|e| e.to_string())?;
        let execution = execution.canonicalize().map_err(|e| e.to_string())?;
        let (record, hash) = capsule::check_capsule(&root)?;
        expected_seal(&record, &hash, expected)?;
        require(
            frozen_self(&root, &record)? == own_sha,
            "ordinary audit checker differs",
        )?;
        capsule::claim(&root, &execution, &record, &hash)?;
        verify_build(&root, &record)?;
        let cfg = capsule::config(&root)?;
        let terminal = load(&execution.join("terminal.json"))?;
        require(
            terminal["registration_sha256"] == hash
                && terminal["scope"] == capsule::SCOPE
                && terminal["controller_sha256"] == own_sha
                && terminal["controller_unchanged"] == true
                && terminal["source_gate_passed"] == true
                && terminal["worker_drain_passed"] == true
                && terminal["role_drain_passed"] == true,
            "ordinary terminal source/drain gate failed",
        )?;
        verify_execution_inventory(&execution, &terminal)?;
        let receipt = load(&execution.join("worker.receipt.json"))?;
        let argv = std::iter::once(
            root.join("immutable/bin")
                .join(capsule::WORKER)
                .to_string_lossy()
                .into_owned(),
        )
        .chain(worker_args(&root, &execution))
        .collect::<Vec<_>>();
        let worker_pid = verify_receipt(
            &execution,
            &execution.join("worker"),
            &receipt,
            &argv,
            &record.worker_sha256,
            cfg.worker_timeout_ms,
            &json!(capsule::environment()),
        )?;
        require(
            receipt["exit_code"] == 0
                && receipt["timed_out"] == false
                && terminal["worker"]["receipt"] == receipt
                && terminal["worker"]["exit_code"] == 0
                && terminal["worker"]["timed_out"] == false,
            "ordinary worker is incomplete or interrupted; prefix remains diagnostic",
        )?;
        require(
            std::str::from_utf8(&read(&execution.join("worker-pids"), 65536)?)
                .map_err(|e| e.to_string())?
                == format!("start {worker_pid}\ndone {worker_pid}\n"),
            "ordinary worker ledger differs",
        )?;
        require(
            load(&execution.join("worker-started.json"))?
                == json!({"registration_sha256":hash,
            "worker_sha256":record.worker_sha256,"worker_build_identity":record.worker_build_identity,"pid":worker_pid}),
            "ordinary worker start identity differs",
        )?;
        let worker_terminal = load(&execution.join("worker-terminal.json"))?;
        require(
            worker_terminal["child_drain_confirmed"] == true
                && worker_terminal["immutable_binding_unchanged"] == true
                && worker_terminal["source_bound_execution_admitted"] == false,
            "ordinary producer terminal binding differs",
        )?;
        let retained = prefix(&execution, &cfg)?;
        require(
            retained["panel_count_complete"] == true && retained["pending_start"].is_null(),
            "ordinary planned panel not complete",
        )?;
        let bytes = read(&execution.join("producer.json"), 16 * 1024 * 1024)?;
        let report: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
        require(
            report["schema_version"] == 1
                && report["question"] == "native-target-free-preparation-producer-n17-v1"
                && report["plan"] == json!(cfg.plan)
                && report["query_records"] == retained["records"]
                && report["mathematical_input"]
                    == load(&execution.join("mathematical-input.json"))?,
            "ordinary producer, durable prefix and mathematical input differ",
        )?;
        let rows = report["query_records"]
            .as_array()
            .ok_or("missing ordinary records")?;
        require(
            report["mathematical_input"]["attempts"]
                == json!(rows
                    .iter()
                    .map(|r| r["attempt"].clone())
                    .collect::<Vec<_>>())
                && report["mathematical_input"]["plan"]["family"] == json!(cfg.plan.family)
                && report["mathematical_input"]["plan"]["algorithm_seed"]
                    == cfg.plan.algorithm_seed
                && report["mathematical_input"]["plan"]["planned_queries"]
                    == cfg.plan.planned_queries,
            "ordinary mathematical query plan differs",
        )?;
        let mathematics = ordinary_preparation::audit(&report["mathematical_input"])?;
        verify_matrix(&report, &mathematics)?;
        let wall = verify_costs(&report, cfg.plan.planned_queries)?;
        require(
            wall <= integer(&receipt["child_wall_ns"])? && report["panel_complete"] == true,
            "preparation interval exceeds original worker interval",
        )?;
        let roles = if cfg.plan.family == Family::Cryptominisat {
            verify_native_assets(&root)?;
            verify_cms_rows(&root, &execution, &cfg, &record, &report)?
        } else {
            require(
                !execution.join("child-pids").exists(),
                "native F5 unexpectedly launched external roles",
            )?;
            0
        };
        capsule::check_capsule(&root)?;
        Ok(
            json!({"schema_version":1,"status":"PASS_NATIVE_SOURCE_BOUND_ORDINARY_PREPARATION_AUDIT",
            "registration_sha256":hash,"source_bound_execution_admitted":true,
            "mathematical_preparation_complete":mathematics["mathematical_preparation_complete"],"mathematics":mathematics,
            "producer_sha256":sha256(&bytes),"worker_pid":worker_pid,"audited_worker_calls":1,
            "audited_role_calls":roles,"audited_ordinary_queries":rows.len(),"preparation_wall_ns":wall,
            "exclusive_preparation_phases_ns":report["exclusive_preparation_phases_ns"],
            "outer_worker_wall_ns":receipt["child_wall_ns"],"worker_setup_and_tail_ns":integer(&receipt["child_wall_ns"])?-wall,
            "wall_timing_class":"ordinary-host-stage-diagnostic","peak_memory":null,"operation_counts":null}),
        )
    })();
    let mut receipt = match result {
        Ok(result) => result,
        Err(error) => json!({"schema_version":1,
        "status":"TERMINAL_ORDINARY_PREPARATION_NOT_ADMITTED","registration_sha256":expected,
        "source_bound_execution_admitted":false,"mathematical_preparation_complete":false,"error":error}),
    };
    require(
        sha256(&read(&own, 128 * 1024 * 1024)?) == own_sha,
        "ordinary auditor changed during verification",
    )?;
    receipt["checker_sha256"] = json!(own_sha);
    receipt["independent_audit_wall_ns"] = json!(started.elapsed().as_nanos().to_string());
    receipt["native_children_executed_by_auditor"] = json!(0);
    receipt["fresh_paired_qualification"] = json!(false);
    receipt["headline_eligible"] = json!(false);
    receipt["promotion_eligible"] = json!(false);
    receipt["full_goal_complete"] = json!(false);
    receipt["online_wall_ns"] = Value::Null;
    receipt["online_speedup"] = Value::Null;
    save(out, &receipt)?;
    serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    fn config() -> capsule::Config {
        capsule::Config {
            schema_version: 1,
            question: capsule::SCOPE.into(),
            plan: crypto_lib::cryptanalysis::prepared_ordinary::Plan {
                schema_version: 1,
                question: "native-target-free-preparation-n17-v1".into(),
                family: Family::MatrixF5,
                algorithm_seed: 2026100311,
                planned_queries: 3,
            },
            export_nonce: 1,
            conflict_budget: 100000,
            exporter_timeout_ms: 30000,
            solver_timeout_ms: 30000,
            worker_timeout_ms: 120000,
        }
    }
    fn temp(label: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!(
            "ordinary-controller-{label}-{}",
            std::process::id()
        ));
        fs::create_dir(&p).unwrap();
        p
    }
    #[test]
    fn partial_prefix_retains_final_start_and_rejects_gaps_and_changed_law() {
        let p = temp("prefix");
        let cfg = config();
        let scalar = |t| {
            crypto_lib::cryptanalysis::koblitz_index_calculus::probe_scalar(
                cfg.plan.algorithm_seed,
                t,
                65587,
            )
        };
        save(
            &p.join("start-000.json"),
            &json!({"trial":0,"scalar":scalar(0)}),
        )
        .unwrap();
        save(
            &p.join("progress-000.json"),
            &json!({"attempt":{"trial":0,"scalar":scalar(0),"outcome":"timeout","indices":null}}),
        )
        .unwrap();
        save(
            &p.join("start-001.json"),
            &json!({"trial":1,"scalar":scalar(1)}),
        )
        .unwrap();
        let retained = prefix(&p, &cfg).unwrap();
        assert_eq!(retained["completed_queries"], 1);
        assert_eq!(retained["pending_start"]["trial"], 1);
        assert_eq!(retained["panel_count_complete"], false);
        assert_eq!(retained["source_bound_execution_admitted"], false);
        save(&p.join("progress-002.json"), &json!({})).unwrap();
        assert!(prefix(&p, &cfg).is_err());
        fs::remove_file(p.join("progress-002.json")).unwrap();
        fs::write(
            p.join("start-001.json"),
            serde_json::to_vec(&json!({"trial":1,"scalar":scalar(1)+1})).unwrap(),
        )
        .unwrap();
        assert!(prefix(&p, &cfg).is_err());
        fs::remove_dir_all(p).unwrap();
    }
    #[test]
    fn ancestor_configuration_and_worker_argv_are_explicit() {
        let p = temp("ancestor");
        let source = p.join("source");
        fs::create_dir(&source).unwrap();
        reject_ancestor_config(&source).unwrap();
        fs::create_dir(p.join(".cargo")).unwrap();
        fs::write(p.join(".cargo/config.toml"), b"[build]\n").unwrap();
        assert!(reject_ancestor_config(&source).is_err());
        assert_eq!(
            worker_args(Path::new("/capsule"), Path::new("/execution")),
            vec!["run", "--capsule", "/capsule", "--execution", "/execution"]
        );
        fs::remove_dir_all(p).unwrap();
    }
    // These numbers are accounting controls, never measured performance.
    fn snapshot(charges: &[(&str, u64)]) -> Value {
        let names = [
            "setup",
            "factor_base",
            "precompute",
            "queries",
            "pdp",
            "relation_check",
            "matrix_build",
            "relation_la",
            "target_query",
            "target_pdp",
            "target_relation_check",
            "target_descent",
            "recovery_check",
            "rho_solve",
        ];
        let mut phases = serde_json::Map::new();
        for name in names {
            phases.insert(name.into(), Value::Null);
        }
        for &(name, cost) in charges {
            phases.insert(name.into(), json!(cost));
        }
        let mut online = serde_json::Map::new();
        for name in &names[8..] {
            online.insert((*name).into(), Value::Null);
        }
        json!({"schema_version":1,"phases_ns":phases,"observed_wall_ns":charges.iter().map(|v|v.1).sum::<u64>(),
            "online_phases_ns":online,"online_wall_ns":null})
    }
    #[test]
    fn exclusive_costs_reject_missing_attempt_phases_and_inconsistent_closure() {
        let report = json!({"phase_windows":[{"stage":"reusable_initialization","clock":snapshot(&[("setup",2)])},
            {"stage":"final_relation_la_and_group_replay","clock":snapshot(&[("relation_la",5)])}],
            "query_records":[{"costs":snapshot(&[("queries",1),("pdp",3),("relation_check",1)])}],
            "controller_setup_ns":1,"preparation_wall_ns":13,
            "exclusive_preparation_phases_ns":{"setup":3,"queries":1,"pdp":3,"relation_check":1,"relation_la":5},
            "online_wall_ns":null,"online_speedup":null,"source_bound_execution_admitted":false,"full_goal_complete":false});
        assert_eq!(verify_costs(&report, 1).unwrap(), 13);
        let mut changed = report.clone();
        changed["preparation_wall_ns"] = json!(12);
        assert!(verify_costs(&changed, 1).is_err());
        let mut changed = report.clone();
        changed["query_records"][0]["costs"]["phases_ns"]
            .as_object_mut()
            .unwrap()
            .remove("pdp");
        assert!(verify_costs(&changed, 1).is_err());
        let mut changed = report.clone();
        changed["query_records"][0]["costs"]["phases_ns"]["pdp"] = Value::Null;
        assert!(verify_costs(&changed, 1).is_err());
        let mut changed = report.clone();
        changed["query_records"][0]["costs"]["online_phases_ns"]["target_pdp"] = json!(0);
        assert!(verify_costs(&changed, 1).is_err());
        let mut changed = report.clone();
        changed["exclusive_preparation_phases_ns"]["pdp"] = json!(2);
        assert!(verify_costs(&changed, 1).is_err());
        assert!(verify_costs(&report, 2).is_err());
    }
    #[test]
    fn matrix_representation_rejects_wrong_coefficients_columns_and_duplicate_entries() {
        // Synthetic representation check. Independent group reconstruction is a separate gate.
        let columns = (0..29u64).map(|i| [i, i + 1]).collect::<Vec<_>>();
        let mut coefficients = vec![0u64; 29];
        coefficients[5] = 1;
        let math = json!({"column_points":columns,"accepted_rows":1,"duplicate_relations":0,
            "rows_sha256":canonical_sha(&json!([{"coefficients":coefficients,"rhs":2}])).unwrap()});
        let report = json!({"relation_matrix":{"modulus":"65587","solver":"dense-gauss","sparse_options":null,
            "column_points":columns.iter().map(|p|p.map(|v|v.to_string())).collect::<Vec<_>>(),
            "rows":[{"entries":[[5,"1"]],"rhs":"2"}]},"log_report":{"relations":1,"duplicate_relations":0}});
        verify_matrix(&report, &math).unwrap();
        let mut changed = report.clone();
        changed["relation_matrix"]["rows"][0]["entries"][0][1] = json!("2");
        assert!(verify_matrix(&changed, &math).is_err());
        let mut changed = report.clone();
        changed["relation_matrix"]["column_points"]
            .as_array_mut()
            .unwrap()
            .swap(0, 1);
        assert!(verify_matrix(&changed, &math).is_err());
        let mut changed = report.clone();
        changed["relation_matrix"]["rows"][0]["entries"] = json!([[5, "1"], [5, "1"]]);
        assert!(verify_matrix(&changed, &math).is_err());
    }
    #[test]
    fn exact_receipt_rejects_environment_drain_and_output_drift() {
        let p = temp("receipt");
        let stem = p.join("worker");
        fs::write(stem.with_extension("stdout"), b"data").unwrap();
        fs::write(stem.with_extension("stderr"), b"").unwrap();
        let argv = vec!["/frozen/worker".into(), "run".into()];
        let env = json!(capsule::environment());
        let receipt = json!({"argv":argv,"cwd":p,"executable_sha256_before":"pin","executable_sha256_after":"pin",
            "deadline_ms":10,"environment":env,"stdout_sha256":sha256(b"data"),"stdout_bytes":4,
            "stderr_sha256":sha256(b""),"stderr_bytes":0,"output_limit":false,
            "process_group_drain_requested":true,"process_group_drain_confirmed":true,"child_wall_ns":7,"pid":123});
        assert_eq!(
            verify_receipt(&p, &stem, &receipt, &argv, "pin", 10, &env).unwrap(),
            123
        );
        let mut changed = receipt.clone();
        changed["environment"]["RAYON_NUM_THREADS"] = json!("2");
        assert!(verify_receipt(&p, &stem, &changed, &argv, "pin", 10, &env).is_err());
        let mut changed = receipt.clone();
        changed["process_group_drain_confirmed"] = json!(false);
        assert!(verify_receipt(&p, &stem, &changed, &argv, "pin", 10, &env).is_err());
        fs::write(stem.with_extension("stdout"), b"drift").unwrap();
        assert!(verify_receipt(&p, &stem, &receipt, &argv, "pin", 10, &env).is_err());
        fs::remove_dir_all(p).unwrap();
    }
    #[test]
    fn drained_execution_inventory_rejects_modified_and_added_evidence() {
        let p = temp("inventory");
        fs::write(p.join("producer.json"), b"{}").unwrap();
        let terminal = json!({"execution_files":native::inventory(&p).unwrap(),"execution_inventory_error":null});
        save(&p.join("terminal.json"), &terminal).unwrap();
        verify_execution_inventory(&p, &terminal).unwrap();
        fs::write(p.join("extra"), b"changed").unwrap();
        assert!(verify_execution_inventory(&p, &terminal).is_err());
        fs::remove_file(p.join("extra")).unwrap();
        fs::write(p.join("producer.json"), b"[]").unwrap();
        assert!(verify_execution_inventory(&p, &terminal).is_err());
        fs::remove_dir_all(p).unwrap();
    }
    #[test]
    fn build_receipts_accept_actual_pair_encoding_and_reject_flag_or_log_drift() {
        let p = temp("build-receipts");
        let dir = p.join("immutable/build-receipts");
        fs::create_dir_all(&dir).unwrap();
        save(&p.join("config.json"), &json!(config())).unwrap();
        let build = json!({"schema_version":1,"scope":capsule::SCOPE,"source_manifest_sha256":"sourcepin",
            "cargo_sha256":"cargopin","rustc_sha256":"rustcpin","profile":"release","features":"default",
            "algorithm_flags":"none","hardware":{"os":"macos","architecture":"aarch64"},"native_archive_sha256":null});
        save(&dir.join("build-identity.json"), &build).unwrap();
        let identity = json!({"source_manifest_sha256":"sourcepin","build_sha256":canonical_sha(&build).unwrap()});
        let record = capsule::Registration {
            schema_version: 1,
            scope: capsule::SCOPE.into(),
            source_commit: "commit".into(),
            hardware: build["hardware"].clone(),
            immutable_files: json!({}),
            source_manifest_sha256: "sourcepin".into(),
            config_sha256: "config".into(),
            host_context_sha256: "host".into(),
            worker_sha256: "workerpin".into(),
            auditor_sha256: "auditor".into(),
            worker_build_identity: identity.clone(),
            runtime_environment: capsule::environment(),
            validation_only: true,
        };
        let env = json!([
            ["PATH", "/tools"],
            ["LC_ALL", "C"],
            ["RUSTC", "/tools/rustc"],
            ["CARGO_HOME", "/home"],
            ["CARGO_TARGET_DIR", "/target"],
            ["CARGO_INCREMENTAL", "0"],
            ["IC_ORDINARY_SOURCE_MANIFEST_SHA256", "sourcepin"],
            ["IC_ORDINARY_BUILD_SHA256", canonical_sha(&build).unwrap()]
        ]);
        for name in ["vendor", "rustc", "cargo", "build", "worker-identity"] {
            let bytes = if name == "worker-identity" {
                serde_json::to_vec(&identity).unwrap()
            } else {
                format!("{name} control log").into_bytes()
            };
            fs::write(dir.join(format!("{name}.log")), &bytes).unwrap();
            let pin = match name {
                "rustc" => "rustcpin",
                "worker-identity" => "workerpin",
                _ => "cargopin",
            };
            let args = match name {
                "vendor" => json!([
                    "vendor",
                    "--locked",
                    "--offline",
                    p.join("immutable/source/vendor")
                ]),
                "build" => json!([
                    "build",
                    "--locked",
                    "--offline",
                    "--release",
                    "--bin",
                    capsule::WORKER,
                    "--bin",
                    "icprog"
                ]),
                "worker-identity" => json!(["build-identity"]),
                _ => json!(["--version", "--verbose"]),
            };
            let run_env = match name {
                "vendor" => json!([["PATH", "/tools"], ["LC_ALL", "C"]]),
                "worker-identity" => json!(capsule::environment().into_iter().collect::<Vec<_>>()),
                _ => env.clone(),
            };
            save(&dir.join(format!("{name}.receipt.json")),&json!({"exit_code":0,"program_sha256_before":pin,
                "program_sha256_after":pin,"log_sha256":sha256(&bytes),"argv":args,"environment":run_env})).unwrap();
        }
        verify_build(&p, &record).unwrap();
        let path = dir.join("build.receipt.json");
        let original = load(&path).unwrap();
        let mut changed = original.clone();
        changed["environment"][5][1] = json!("1");
        fs::write(&path, serde_json::to_vec(&changed).unwrap()).unwrap();
        assert!(verify_build(&p, &record).is_err());
        fs::write(&path, serde_json::to_vec(&original).unwrap()).unwrap();
        fs::write(dir.join("rustc.log"), b"different version").unwrap();
        assert!(verify_build(&p, &record).is_err());
        fs::remove_dir_all(p).unwrap();
    }
}
