//! Native composer/verifier for Stage 193 WDSat and CryptoMiniSat controls.

use crypto_lib::hash::sha256::sha256;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

const RESULT_SCHEMA: &str = "koblitz_stage193_sat_comparators.v1";
const VERIFICATION_SCHEMA: &str = "koblitz_stage193_sat_comparators_verification.v1";
const METER_SCHEMA: &str = "koblitz_native_process_meter.v1";
const RECEIPT_SCHEMA: &str = "koblitz_native_meter_receipt.v1";
const SOURCE_INSTANCE_ID: &str = "954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7";
const ANF_SHA256: &str = "34c516d2568e327a3c5b403d6871845c80d523f0b7f012bd183be0bbc1469508";
const XOR_SHA256: &str = "a529c4d8bb1e1a4b5904134b1e5d099adfb80806215f1ff128088d7d518d4bb8";
const MANIFEST_SHA256: &str = "188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c";
const WDSAT_CONFIG_SHA256: &str =
    "4a73c3b5a14ded98f749d282b597594ae3ac05355a39cadcb9345a90bf735352";
const WDSAT_SOURCE_COMMIT: &str = "55d55b2620d768d9f7c78dcd8990a0689533c1d0";
const CMS_VERSION: &str = "5.14.7";
const CMS_BINARY_SHA256: &str = "a3f85c3709b5e2a040bf82a4a604d1c7b9f10219bbf180a9e0f72319a2e892ac";
const METER_BINARY_PATH: &str =
    "/Volumes/SSD990/kic-stage192-target-e575d53f1/release/examples/koblitz_f4_stage188";
const METER_BINARY_SHA256: &str =
    "971ddd3b0a77b86aa1388632a201747cec42d570e4d9ba9536c7c61e8361c524";
const STAGE192_SHA256: &str = "e1849a1b2ae97bdc3633025eaadcb18438ff2a49448820a936a890db24904a48";
const METER_DESCRIPTION: &str = "fresh child process wait4/getrusage; all child threads included";
const THREAD_CAPS: [&str; 7] = [
    "RAYON_NUM_THREADS=1",
    "VECLIB_MAXIMUM_THREADS=1",
    "OPENBLAS_NUM_THREADS=1",
    "OMP_NUM_THREADS=1",
    "MKL_NUM_THREADS=1",
    "BLIS_NUM_THREADS=1",
    "NUMEXPR_NUM_THREADS=1",
];
const WDSAT_SOURCE_INVENTORY: [(&str, &str); 16] = [
    (
        "cnf.c",
        "ec86f8f08d9272579b20d6dfb9513ca69bdfca3c490e4b11e75878f13745a9ee",
    ),
    (
        "cnf.h",
        "f1c5479b60c4a49192c84569a74a9838b855c44411e35c066817076c90879a62",
    ),
    (
        "dimacs.c",
        "d165dea3b5c2919c8d29a05a037e0b98ed7aac4ff9f7c8c073a7a9c2151d9ffb",
    ),
    (
        "dimacs.h",
        "f40ede80353b747851463b0fd376c3b9f7e394367b21df115ab8f5f797b30cfd",
    ),
    (
        "fp.c",
        "c51eca3b78beb801810ab82a9ae5d97af688338d567ec4c645873ec4adf1fd06",
    ),
    (
        "fp.h",
        "cff83916044c1370c0c07f271fa8eea6446ff1e29032618821cdda3c75b68e97",
    ),
    (
        "main.c",
        "0ab5b48557aa3a5a6229124bda3dea36a1a771d2a8f92721e6d44a4291b0973d",
    ),
    (
        "makefile",
        "3cbeb6537587d54d28c37937e17d841b551d774f7a9946aeee91264e9413b48d",
    ),
    (
        "wdsat.c",
        "1f17ba03459bbdda89539e542bbda75d8860037082a25d7c073dc0d90095a732",
    ),
    (
        "wdsat.h",
        "476b9ae9fcf66b140625e6b8f44179a121bcc08cb425f63dd7b6caf02c1d6e8a",
    ),
    (
        "wdsat_utils.c",
        "430b2d23de167854782cfd46dd0e08987abcd1d69cfde69d82482dd81899740e",
    ),
    (
        "wdsat_utils.h",
        "15cd7a3921902d157e7c83bc333c516ec8ba685e91170c897e06c5fe3654479b",
    ),
    (
        "xorgauss.c",
        "e9c62dda272a3801999b03da2ada8282539796caf98befd6768ce79027e526ec",
    ),
    (
        "xorgauss.h",
        "02014f27602a0e02553467acde33d14c69c7aa0f8a5fc01f16deece57b438445",
    ),
    (
        "xorset.c",
        "584ef01e43561df607a7e65d340dcd30736121834c4b66f81f9f208b6d924c0d",
    ),
    (
        "xorset.h",
        "0cb90df166194472d21bcb7740d5ab2c1887db89eea9e2522a4d06e977099db9",
    ),
];

type AnyResult<T> = Result<T, String>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Artifact {
    path: String,
    bytes: u64,
    sha256: String,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct ProcessMetrics {
    schema: String,
    meter: String,
    command: Vec<String>,
    environment: BTreeMap<String, String>,
    pid: u32,
    returncode: i32,
    timed_out: bool,
    watchdog_seconds: u64,
    wall_seconds: f64,
    user_seconds: Option<f64>,
    system_seconds: Option<f64>,
    total_core_seconds: Option<f64>,
    single_core_seconds: Option<f64>,
    peak_rss_bytes: Option<u64>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct MeterReceipt {
    schema: String,
    process: ProcessMetrics,
    stdout: Artifact,
    stderr: Artifact,
    metrics: Artifact,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct SolverOutcome {
    solver: String,
    status: String,
    terminal: Option<String>,
    conflicts: Option<u64>,
    wall_seconds: f64,
    total_core_seconds: f64,
    single_core_seconds: f64,
    peak_rss_bytes: u64,
    timed_out: bool,
    returncode: i32,
    receipt: MeterReceipt,
    stdout: Artifact,
    stderr: Artifact,
    binary: Artifact,
    input: Artifact,
    executed_input: Artifact,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct WdsatProvenance {
    source_commit: String,
    worktree: MeterReceipt,
    source_check: MeterReceipt,
    config_install: MeterReceipt,
    clean: MeterReceipt,
    build: MeterReceipt,
    installed_config: Artifact,
    source_inventory: Vec<Artifact>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct CmsProvenance {
    version: String,
    version_check: MeterReceipt,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Charge {
    components: usize,
    wall_seconds_sum: f64,
    total_core_seconds_sum: f64,
    peak_rss_bytes_max: u64,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct ResultRecord {
    schema: String,
    source_commit: String,
    finalizer_commit: String,
    source_instance_id: String,
    factor_base_contract: Value,
    meter_binary: Artifact,
    manifest: Artifact,
    wdsat_anf: Artifact,
    cryptominisat_xor_dimacs: Artifact,
    wdsat_config: Artifact,
    wdsat_provenance: WdsatProvenance,
    cryptominisat_provenance: CmsProvenance,
    wdsat: SolverOutcome,
    cryptominisat: SolverOutcome,
    inherited_stage192_result: Artifact,
    inherited_stage192_result_sha256: String,
    measured_stage_charge: Charge,
    cumulative_measured_lower_bound: Charge,
    complete_campaign_cost: Option<Charge>,
    licensed_magma_available: bool,
    all_seven_gates_passed: bool,
    koblitz_index_calculus_sota: bool,
    claim_boundary: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct Verification {
    schema: String,
    result_sha256: String,
    checks: usize,
    passed: usize,
    failed: Vec<String>,
}

fn main() {
    if let Err(error) = real_main() {
        eprintln!("{error}");
        std::process::exit(2);
    }
}

fn real_main() -> AnyResult<()> {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("compose") if args.len() == 5 => compose(
            Path::new(&args[1]),
            Path::new(&args[2]),
            &args[3],
            &args[4],
        ),
        Some("verify") if args.len() == 3 => {
            let verification = verify(Path::new(&args[1]))?;
            write_json(Path::new(&args[2]), &verification)?;
            if verification.failed.is_empty() {
                Ok(())
            } else {
                Err(format!(
                    "verification failed: {}",
                    verification.failed.join("; ")
                ))
            }
        }
        _ => Err("usage: koblitz_stage193_sat_comparators compose STAGE RESULT SOURCE_COMMIT FINALIZER_COMMIT\n       koblitz_stage193_sat_comparators verify RESULT VERIFICATION".into()),
    }
}

fn compose(stage: &Path, out: &Path, source_commit: &str, finalizer_commit: &str) -> AnyResult<()> {
    validate_commit("source", source_commit)?;
    validate_commit("finalizer", finalizer_commit)?;
    let stage = stage
        .canonicalize()
        .map_err(|e| format!("{}: {e}", stage.display()))?;
    let manifest = artifact(&stage, &stage.join("input/manifest.json"))?;
    let meter_binary = artifact(&stage, Path::new(METER_BINARY_PATH))?;
    let anf = artifact(&stage, &stage.join("input/instance.anf"))?;
    let xor = artifact(&stage, &stage.join("input/instance.xor.cnf"))?;
    let config = artifact(&stage, &stage.join("input/wdsat-config.h"))?;
    require_hash("manifest", &manifest, MANIFEST_SHA256)?;
    require_hash("native process meter", &meter_binary, METER_BINARY_SHA256)?;
    require_hash("ANF", &anf, ANF_SHA256)?;
    require_hash("XOR-DIMACS", &xor, XOR_SHA256)?;
    require_hash("WDSat config", &config, WDSAT_CONFIG_SHA256)?;
    validate_manifest(&stage.join("input/manifest.json"))?;
    let wdsat_provenance = load_wdsat_provenance(&stage)?;
    let wdsat_binary = wdsat_root(&wdsat_provenance)?.join("wdsat_solver");
    let cryptominisat_provenance = load_cms_provenance(&stage)?;
    let wdsat = solver_outcome(
        &stage,
        "wdsat",
        "development/wdsat-run",
        &anf,
        &wdsat_binary,
    )?;
    let cms = solver_outcome(
        &stage,
        "cryptominisat",
        "development/cryptominisat-run",
        &xor,
        Path::new("/opt/homebrew/bin/cryptominisat5"),
    )?;
    let inherited_path = stage
        .parent()
        .ok_or_else(|| "stage has no continuation parent".to_string())?
        .join("stage-192-selected-parallel-boundary-20261002/result.json");
    let inherited = artifact(&stage, &inherited_path)?;
    require_hash("Stage 192 result", &inherited, STAGE192_SHA256)?;
    let measured_stage_charge = measured_charge(&stage.join("development"))?;
    let cumulative_measured_lower_bound = Charge {
        components: 582 + measured_stage_charge.components,
        wall_seconds_sum: 23_910.746_620_539 + measured_stage_charge.wall_seconds_sum,
        total_core_seconds_sum: 61_559.696_526 + measured_stage_charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(measured_stage_charge.peak_rss_bytes_max),
    };
    write_json(
        out,
        &ResultRecord {
            schema: RESULT_SCHEMA.into(),
            source_commit: source_commit.into(),
            finalizer_commit: finalizer_commit.into(),
            source_instance_id: SOURCE_INSTANCE_ID.into(),
            factor_base_contract: serde_json::json!({
                "definition":"span_F2(1,z,...,z^8)",
                "target_subgroup_enumerated":false,
                "discrete_log_labels_used":false,
            }),
            meter_binary,
            manifest,
            wdsat_anf: anf,
            cryptominisat_xor_dimacs: xor,
            wdsat_config: config,
            wdsat_provenance,
            cryptominisat_provenance,
            wdsat,
            cryptominisat: cms,
            inherited_stage192_result: inherited,
            inherited_stage192_result_sha256: STAGE192_SHA256.into(),
            measured_stage_charge,
            cumulative_measured_lower_bound,
            complete_campaign_cost: None,
            licensed_magma_available: false,
            all_seven_gates_passed: false,
            koblitz_index_calculus_sota: false,
            claim_boundary: "Same opened n=59 decomposition instance. WDSat and CryptoMiniSat terminal/resource controls only; UNKNOWN and watchdog timeout are censored. Licensed Magma, full index calculus, rho crossover, external reproduction, novelty and SOTA remain open.".into(),
        },
    )
}

fn verify(path: &Path) -> AnyResult<Verification> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let result: ResultRecord =
        serde_json::from_slice(&bytes).map_err(|e| format!("result JSON: {e}"))?;
    let stage_root = find_stage_root(path)?;
    let stage = stage_root
        .canonicalize()
        .map_err(|e| format!("{}: {e}", stage_root.display()))?;
    let mut checks = 0usize;
    let mut failures = Vec::new();
    check(
        &mut checks,
        &mut failures,
        result.schema == RESULT_SCHEMA && result.source_instance_id == SOURCE_INSTANCE_ID,
        "result identity",
    );
    check(
        &mut checks,
        &mut failures,
        validate_manifest(&stage.join("input/manifest.json")).is_ok(),
        "manifest contract",
    );
    for (name, artifact, expected) in [
        (
            "native process meter",
            &result.meter_binary,
            METER_BINARY_SHA256,
        ),
        ("manifest", &result.manifest, MANIFEST_SHA256),
        ("ANF", &result.wdsat_anf, ANF_SHA256),
        ("XOR", &result.cryptominisat_xor_dimacs, XOR_SHA256),
        ("WDSat config", &result.wdsat_config, WDSAT_CONFIG_SHA256),
        (
            "Stage 192",
            &result.inherited_stage192_result,
            STAGE192_SHA256,
        ),
    ] {
        checks += 2;
        match rehash_artifact(&stage, artifact) {
            Ok(actual) => {
                if actual != *artifact {
                    failures.push(format!("{name} artifact replay"));
                }
                if artifact.sha256 != expected {
                    failures.push(format!("{name} expected hash"));
                }
            }
            Err(error) => failures.push(format!("{name}: {error}")),
        }
    }
    let wdsat_provenance = load_wdsat_provenance(&stage)?;
    check(
        &mut checks,
        &mut failures,
        wdsat_provenance == result.wdsat_provenance,
        "WDSat provenance replay",
    );
    let cms_provenance = load_cms_provenance(&stage)?;
    check(
        &mut checks,
        &mut failures,
        cms_provenance == result.cryptominisat_provenance,
        "CryptoMiniSat provenance replay",
    );
    let wdsat_binary = wdsat_root(&wdsat_provenance)?.join("wdsat_solver");
    let wdsat = solver_outcome(
        &stage,
        "wdsat",
        "development/wdsat-run",
        &result.wdsat_anf,
        &wdsat_binary,
    )?;
    let cms = solver_outcome(
        &stage,
        "cryptominisat",
        "development/cryptominisat-run",
        &result.cryptominisat_xor_dimacs,
        Path::new("/opt/homebrew/bin/cryptominisat5"),
    )?;
    check(
        &mut checks,
        &mut failures,
        wdsat == result.wdsat,
        "WDSat replay",
    );
    check(
        &mut checks,
        &mut failures,
        cms == result.cryptominisat,
        "CryptoMiniSat replay",
    );
    let charge = measured_charge(&stage.join("development"))?;
    check(
        &mut checks,
        &mut failures,
        charges_match(&charge, &result.measured_stage_charge),
        "stage charge",
    );
    let cumulative = Charge {
        components: 582 + charge.components,
        wall_seconds_sum: 23_910.746_620_539 + charge.wall_seconds_sum,
        total_core_seconds_sum: 61_559.696_526 + charge.total_core_seconds_sum,
        peak_rss_bytes_max: 6_310_576_128u64.max(charge.peak_rss_bytes_max),
    };
    check(
        &mut checks,
        &mut failures,
        charges_match(&cumulative, &result.cumulative_measured_lower_bound),
        "cumulative charge",
    );
    check(
        &mut checks,
        &mut failures,
        result.factor_base_contract["definition"] == "span_F2(1,z,...,z^8)"
            && result.factor_base_contract["target_subgroup_enumerated"] == Value::Bool(false)
            && result.factor_base_contract["discrete_log_labels_used"] == Value::Bool(false)
            && result.wdsat_provenance.source_commit == WDSAT_SOURCE_COMMIT
            && result.cryptominisat_provenance.version == CMS_VERSION
            && result.complete_campaign_cost.is_none()
            && !result.licensed_magma_available
            && !result.all_seven_gates_passed
            && !result.koblitz_index_calculus_sota,
        "claim boundary",
    );
    Ok(Verification {
        schema: VERIFICATION_SCHEMA.into(),
        result_sha256: hex::encode(sha256(&bytes)),
        checks,
        passed: checks.saturating_sub(failures.len()),
        failed: failures,
    })
}

fn solver_outcome(
    stage: &Path,
    solver: &str,
    relative: &str,
    input: &Artifact,
    expected_binary: &Path,
) -> AnyResult<SolverOutcome> {
    let directory = stage.join(relative);
    let receipt: MeterReceipt = read_json(&directory.join("receipt.json"))?;
    validate_receipt(&directory, &receipt)?;
    let (binary_path, executed_input_path) =
        validate_solver_command(stage, solver, &receipt, expected_binary)?;
    let stdout = artifact(&stage, &directory.join(&receipt.stdout.path))?;
    let stderr = artifact(&stage, &directory.join(&receipt.stderr.path))?;
    let output = fs::read_to_string(directory.join(&receipt.stdout.path))
        .map_err(|e| format!("{solver} stdout: {e}"))?;
    let (mut status, terminal, conflicts) = if solver == "cryptominisat" {
        parse_cms(&output, receipt.process.timed_out)
    } else {
        parse_wdsat(&output, receipt.process.timed_out)
    };
    if !receipt.process.timed_out && receipt.process.returncode != 0 {
        status = "process_failure".into();
    }
    validate_solver_terminal(solver, &receipt.process, terminal.as_deref(), conflicts)?;
    let binary = artifact(stage, &binary_path)?;
    if solver == "cryptominisat" {
        require_hash("CryptoMiniSat executable", &binary, CMS_BINARY_SHA256)?;
    }
    let executed_input = artifact(stage, &executed_input_path)?;
    if executed_input.bytes != input.bytes || executed_input.sha256 != input.sha256 {
        return Err(format!("{solver} executed input differs from frozen input"));
    }
    let core = receipt
        .process
        .total_core_seconds
        .ok_or_else(|| format!("{solver} core time unavailable"))?;
    let rss = receipt
        .process
        .peak_rss_bytes
        .ok_or_else(|| format!("{solver} RSS unavailable"))?;
    Ok(SolverOutcome {
        solver: solver.into(),
        status,
        terminal,
        conflicts,
        wall_seconds: receipt.process.wall_seconds,
        total_core_seconds: core,
        single_core_seconds: core,
        peak_rss_bytes: rss,
        timed_out: receipt.process.timed_out,
        returncode: receipt.process.returncode,
        receipt,
        stdout,
        stderr,
        binary,
        input: input.clone(),
        executed_input,
    })
}

fn parse_cms(output: &str, outer_timeout: bool) -> (String, Option<String>, Option<u64>) {
    let terminal = output
        .lines()
        .rev()
        .find_map(|line| line.trim().strip_prefix("s ").map(str::to_owned));
    let conflicts = output.lines().rev().find_map(|line| {
        let line = line.trim();
        if !line.starts_with("c conflicts") {
            return None;
        }
        line.split(':')
            .nth(1)
            .and_then(|tail| tail.split_whitespace().next())
            .and_then(|value| value.parse().ok())
    });
    let status = if outer_timeout {
        "timeout_inconclusive"
    } else {
        match terminal.as_deref() {
            Some("SATISFIABLE") => "sat_unverified",
            Some("UNSATISFIABLE") => "unsat",
            Some("INDETERMINATE" | "UNKNOWN") => "timeout_inconclusive",
            _ => "unknown_inconclusive",
        }
    };
    (status.into(), terminal, conflicts)
}

fn parse_wdsat(output: &str, timed_out: bool) -> (String, Option<String>, Option<u64>) {
    let lines: Vec<&str> = output.lines().map(str::trim).collect();
    let terminal_index = lines
        .iter()
        .rposition(|line| *line == "UNSAT" || line.starts_with("SAT:"));
    let terminal = terminal_index.map(|index| lines[index].to_owned());
    let conflicts = terminal_index.and_then(|index| {
        lines[index + 1..].iter().find_map(|line| {
            line.strip_prefix("conf:")
                .unwrap_or(line)
                .parse::<u64>()
                .ok()
        })
    });
    let status = if timed_out && terminal.is_none() {
        "timeout_inconclusive"
    } else {
        match terminal.as_deref() {
            Some("UNSAT") => "unsat",
            Some(value) if value.starts_with("SAT:") => "sat_unverified",
            _ => "unknown_inconclusive",
        }
    };
    (status.into(), terminal, conflicts)
}

fn validate_solver_command(
    stage: &Path,
    solver: &str,
    receipt: &MeterReceipt,
    expected_binary: &Path,
) -> AnyResult<(PathBuf, PathBuf)> {
    let command = &receipt.process.command;
    let prefix = env_prefix();
    if !command.starts_with(&prefix) {
        return Err(format!("{solver} thread-cap prefix differs from protocol"));
    }
    let tail = &command[prefix.len()..];
    let expected_binary = expected_binary
        .canonicalize()
        .map_err(|e| format!("{}: {e}", expected_binary.display()))?;
    let actual_binary = tail
        .first()
        .ok_or_else(|| format!("{solver} command lacks binary"))?;
    let actual_binary = Path::new(actual_binary)
        .canonicalize()
        .map_err(|e| format!("{actual_binary}: {e}"))?;
    if actual_binary != expected_binary {
        return Err(format!("{solver} executable differs from protocol"));
    }
    let (expected_watchdog, input) = if solver == "wdsat" {
        if tail.len() != 5 || tail[1] != "-i" || tail[3] != "-g" || tail[4] != wdsat_generators() {
            return Err("WDSat arguments differ from protocol".into());
        }
        let input = PathBuf::from(&tail[2]);
        if input.file_name().and_then(|name| name.to_str()) != Some("instance.anf")
            || !input.starts_with("/tmp/kic193")
        {
            return Err("WDSat input is not the frozen short-path copy".into());
        }
        (120, input)
    } else if solver == "cryptominisat" {
        let expected_options = [
            "--threads",
            "1",
            "--verb",
            "1",
            "--verbstat",
            "1",
            "--printsol",
            "0",
            "--zero-exit-status",
            "--maxtime",
            "120",
        ];
        if tail.len() != expected_options.len() + 2
            || !tail[1..tail.len() - 1]
                .iter()
                .map(String::as_str)
                .eq(expected_options)
        {
            return Err("CryptoMiniSat arguments differ from protocol".into());
        }
        let input = PathBuf::from(tail.last().expect("length checked"));
        let expected_input = stage
            .join("input/instance.xor.cnf")
            .canonicalize()
            .map_err(|e| format!("{}: {e}", stage.join("input/instance.xor.cnf").display()))?;
        let actual_input = input
            .canonicalize()
            .map_err(|e| format!("{}: {e}", input.display()))?;
        if actual_input != expected_input {
            return Err("CryptoMiniSat input path differs from protocol".into());
        }
        (135, input)
    } else {
        return Err(format!("unsupported solver {solver}"));
    };
    if receipt.process.watchdog_seconds != expected_watchdog {
        return Err(format!("{solver} watchdog differs from protocol"));
    }
    Ok((actual_binary, input))
}

fn validate_solver_terminal(
    solver: &str,
    process: &ProcessMetrics,
    terminal: Option<&str>,
    conflicts: Option<u64>,
) -> AnyResult<()> {
    if process.timed_out {
        if terminal.is_some() || conflicts.is_some() || process.returncode == 0 {
            return Err(format!("{solver} outer-timeout semantics are inconsistent"));
        }
        return Ok(());
    }
    if process.returncode != 0 {
        if terminal.is_some() {
            return Err(format!("{solver} failed after emitting a terminal result"));
        }
        return Ok(());
    }
    if terminal.is_none() {
        return Err(format!("{solver} clean exit lacks terminal status"));
    }
    if conflicts.is_none() {
        return Err(format!("{solver} terminal status lacks conflict count"));
    }
    Ok(())
}

fn load_wdsat_provenance(stage: &Path) -> AnyResult<WdsatProvenance> {
    let worktree = load_receipt(stage, "development/wdsat-worktree")?;
    require_success("WDSat worktree", &worktree.process, 60)?;
    let command = &worktree.process.command;
    if command.len() != 8
        || command[0] != "/usr/bin/git"
        || command[1] != "-C"
        || command[2] != "/Volumes/SSD990/WDSat"
        || command[3] != "worktree"
        || command[4] != "add"
        || command[5] != "--detach"
        || command[7] != WDSAT_SOURCE_COMMIT
    {
        return Err("WDSat worktree command differs from protocol".into());
    }
    let root = PathBuf::from(&command[6]);
    if !root.is_absolute()
        || root.file_name().and_then(|name| name.to_str()) != Some("wdsat-stage193-55d55b2")
    {
        return Err("WDSat worktree destination differs from protocol".into());
    }

    let source_check = load_receipt(stage, "development/wdsat-source-check")?;
    require_success("WDSat source check", &source_check.process, 30)?;
    let expected_source_command = vec![
        "/usr/bin/git".into(),
        "-C".into(),
        root.display().to_string(),
        "rev-parse".into(),
        "HEAD".into(),
    ];
    if source_check.process.command != expected_source_command {
        return Err("WDSat source-check command differs from protocol".into());
    }
    let source_stdout = fs::read_to_string(
        stage
            .join("development/wdsat-source-check")
            .join(&source_check.stdout.path),
    )
    .map_err(|e| format!("WDSat source-check stdout: {e}"))?;
    if source_stdout.trim() != WDSAT_SOURCE_COMMIT {
        return Err("WDSat source commit differs from protocol".into());
    }

    let config_install = load_receipt(stage, "development/wdsat-config-install")?;
    require_success("WDSat config install", &config_install.process, 30)?;
    let expected_config_command = vec![
        "/bin/cp".into(),
        stage.join("input/wdsat-config.h").display().to_string(),
        root.join("src/config.h").display().to_string(),
    ];
    if config_install.process.command != expected_config_command {
        return Err("WDSat config-install command differs from protocol".into());
    }

    let clean = load_receipt(stage, "development/wdsat-clean")?;
    require_success("WDSat clean", &clean.process, 60)?;
    let mut expected_clean = env_prefix();
    expected_clean.extend([
        "/usr/bin/make".into(),
        "-C".into(),
        root.join("src").display().to_string(),
        "clean".into(),
    ]);
    if clean.process.command != expected_clean {
        return Err("WDSat clean command differs from protocol".into());
    }

    let build = load_receipt(stage, "development/wdsat-build")?;
    require_success("WDSat build", &build.process, 120)?;
    let mut expected_build = env_prefix();
    expected_build.extend([
        "/usr/bin/make".into(),
        "-C".into(),
        root.join("src").display().to_string(),
        "all".into(),
    ]);
    if build.process.command != expected_build {
        return Err("WDSat build command differs from protocol".into());
    }

    let installed_config = artifact(stage, &root.join("src/config.h"))?;
    require_hash(
        "installed WDSat config",
        &installed_config,
        WDSAT_CONFIG_SHA256,
    )?;
    let mut source_inventory = Vec::with_capacity(WDSAT_SOURCE_INVENTORY.len());
    for (name, expected_hash) in WDSAT_SOURCE_INVENTORY {
        let item = artifact(stage, &root.join("src").join(name))?;
        require_hash(name, &item, expected_hash)?;
        source_inventory.push(item);
    }
    Ok(WdsatProvenance {
        source_commit: WDSAT_SOURCE_COMMIT.into(),
        worktree,
        source_check,
        config_install,
        clean,
        build,
        installed_config,
        source_inventory,
    })
}

fn load_cms_provenance(stage: &Path) -> AnyResult<CmsProvenance> {
    let version_check = load_receipt(stage, "development/cryptominisat-version")?;
    require_success("CryptoMiniSat version", &version_check.process, 30)?;
    let mut expected = env_prefix();
    expected.extend([
        "/opt/homebrew/bin/cryptominisat5".into(),
        "--version".into(),
    ]);
    if version_check.process.command != expected {
        return Err("CryptoMiniSat version command differs from protocol".into());
    }
    let binary = artifact(stage, Path::new("/opt/homebrew/bin/cryptominisat5"))?;
    require_hash("CryptoMiniSat executable", &binary, CMS_BINARY_SHA256)?;
    let output = fs::read_to_string(
        stage
            .join("development/cryptominisat-version")
            .join(&version_check.stdout.path),
    )
    .map_err(|e| format!("CryptoMiniSat version stdout: {e}"))?;
    if !output
        .lines()
        .any(|line| line.trim() == format!("c CryptoMiniSat version {CMS_VERSION}"))
    {
        return Err("CryptoMiniSat version output differs from protocol".into());
    }
    Ok(CmsProvenance {
        version: CMS_VERSION.into(),
        version_check,
    })
}

fn load_receipt(stage: &Path, relative: &str) -> AnyResult<MeterReceipt> {
    let directory = stage.join(relative);
    let receipt = read_json(&directory.join("receipt.json"))?;
    validate_receipt(&directory, &receipt)?;
    Ok(receipt)
}

fn wdsat_root(provenance: &WdsatProvenance) -> AnyResult<PathBuf> {
    provenance
        .worktree
        .process
        .command
        .get(6)
        .map(PathBuf::from)
        .ok_or_else(|| "WDSat provenance lacks worktree root".into())
}

fn require_success(name: &str, process: &ProcessMetrics, watchdog: u64) -> AnyResult<()> {
    if process.returncode == 0 && !process.timed_out && process.watchdog_seconds == watchdog {
        Ok(())
    } else {
        Err(format!("{name} did not complete under the frozen watchdog"))
    }
}

fn env_prefix() -> Vec<String> {
    std::iter::once("/usr/bin/env".to_string())
        .chain(THREAD_CAPS.into_iter().map(str::to_string))
        .collect()
}

fn wdsat_generators() -> String {
    (1..=27)
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn validate_manifest(path: &Path) -> AnyResult<()> {
    let manifest: Value = read_json(path)?;
    let expected_basis = (0..9)
        .map(|bit| Value::String((1u16 << bit).to_string()))
        .collect::<Vec<_>>();
    let requirements = [
        ("/n", serde_json::json!(59)),
        ("/ell", serde_json::json!(9)),
        ("/m", serde_json::json!(3)),
        ("/source_variables", serde_json::json!(78)),
        ("/source_equations", serde_json::json!(110)),
        ("/target_mode", serde_json::json!("explicit_affine")),
        (
            "/source_instance/id_blake3",
            serde_json::json!(SOURCE_INSTANCE_ID),
        ),
        (
            "/factor_base_predicate/kind",
            serde_json::json!("polynomial_subspace"),
        ),
        (
            "/factor_base_predicate/definition",
            serde_json::json!("x = sum(c_i z^i), c_i in F_2, 0 <= i < ell"),
        ),
        (
            "/factor_base_predicate/enumerates_target_subgroup",
            serde_json::json!(false),
        ),
        (
            "/factor_base_predicate/uses_discrete_log_labels",
            serde_json::json!(false),
        ),
    ];
    for (pointer, expected) in requirements {
        if manifest.pointer(pointer) != Some(&expected) {
            return Err(format!("manifest {pointer} differs from protocol"));
        }
    }
    if manifest.get("factor_base_basis_bitmasks") != Some(&Value::Array(expected_basis)) {
        return Err("manifest factor-base basis differs from protocol".into());
    }
    Ok(())
}

fn validate_receipt(directory: &Path, receipt: &MeterReceipt) -> AnyResult<()> {
    if receipt.schema != RECEIPT_SCHEMA || receipt.process.schema != METER_SCHEMA {
        return Err(format!("{} meter schema", directory.display()));
    }
    if receipt.process.meter != METER_DESCRIPTION
        || !receipt.process.environment.is_empty()
        || receipt.process.command.is_empty()
        || receipt.process.pid == 0
        || receipt.process.watchdog_seconds == 0
        || !receipt.process.wall_seconds.is_finite()
        || receipt.process.wall_seconds < 0.0
        || receipt.process.single_core_seconds.is_some()
    {
        return Err(format!("{} process metadata", directory.display()));
    }
    let user = receipt
        .process
        .user_seconds
        .ok_or_else(|| format!("{} lacks user time", directory.display()))?;
    let system = receipt
        .process
        .system_seconds
        .ok_or_else(|| format!("{} lacks system time", directory.display()))?;
    let core = receipt
        .process
        .total_core_seconds
        .ok_or_else(|| format!("{} lacks core time", directory.display()))?;
    if !user.is_finite()
        || !system.is_finite()
        || !core.is_finite()
        || user < 0.0
        || system < 0.0
        || core < 0.0
        || !close(user + system, core)
        || receipt.process.peak_rss_bytes.is_none()
    {
        return Err(format!("{} resource accounting", directory.display()));
    }
    if receipt.stdout.path != "stdout.txt"
        || receipt.stderr.path != "stderr.txt"
        || receipt.metrics.path != "metrics.json"
    {
        return Err(format!("{} receipt paths", directory.display()));
    }
    for artifact in [&receipt.stdout, &receipt.stderr, &receipt.metrics] {
        let actual = artifact_from_file(&directory.join(&artifact.path), artifact.path.clone())?;
        if actual != *artifact {
            return Err(format!("{} artifact mismatch", directory.display()));
        }
    }
    let metrics: ProcessMetrics = read_json(&directory.join("metrics.json"))?;
    if metrics != receipt.process {
        return Err(format!("{} metrics replay", directory.display()));
    }
    Ok(())
}

fn measured_charge(development: &Path) -> AnyResult<Charge> {
    let mut files = Vec::new();
    collect_named_files(development, "metrics.json", &mut files)?;
    files.sort();
    let mut receipt_files = Vec::new();
    collect_named_files(development, "receipt.json", &mut receipt_files)?;
    receipt_files.sort();
    let mut covered_metrics = BTreeSet::new();
    for receipt_file in receipt_files {
        let directory = receipt_file
            .parent()
            .ok_or_else(|| format!("{} has no parent", receipt_file.display()))?;
        let receipt: MeterReceipt = read_json(&receipt_file)?;
        validate_receipt(directory, &receipt)?;
        covered_metrics.insert(directory.join(&receipt.metrics.path));
    }
    let measured_metrics = files.iter().cloned().collect::<BTreeSet<_>>();
    if covered_metrics != measured_metrics {
        return Err("development metrics and authenticated receipts differ".into());
    }
    let mut charge = Charge {
        components: files.len(),
        wall_seconds_sum: 0.0,
        total_core_seconds_sum: 0.0,
        peak_rss_bytes_max: 0,
    };
    for file in files {
        let metrics: ProcessMetrics = read_json(&file)?;
        if metrics.schema != METER_SCHEMA {
            return Err(format!("{} wrong meter schema", file.display()));
        }
        charge.wall_seconds_sum += metrics.wall_seconds;
        charge.total_core_seconds_sum += metrics
            .total_core_seconds
            .ok_or_else(|| format!("{} lacks core time", file.display()))?;
        charge.peak_rss_bytes_max = charge.peak_rss_bytes_max.max(
            metrics
                .peak_rss_bytes
                .ok_or_else(|| format!("{} lacks RSS", file.display()))?,
        );
    }
    Ok(charge)
}

fn collect_named_files(directory: &Path, name: &str, out: &mut Vec<PathBuf>) -> AnyResult<()> {
    for entry in fs::read_dir(directory).map_err(|e| format!("{}: {e}", directory.display()))? {
        let entry = entry.map_err(|e| e.to_string())?;
        let kind = entry.file_type().map_err(|e| e.to_string())?;
        if kind.is_symlink() {
            return Err(format!("refusing symlink {}", entry.path().display()));
        }
        if kind.is_dir() {
            collect_named_files(&entry.path(), name, out)?;
        } else if kind.is_file() && entry.file_name() == name {
            out.push(entry.path());
        }
    }
    Ok(())
}

fn artifact(stage: &Path, path: &Path) -> AnyResult<Artifact> {
    let shown = path
        .strip_prefix(stage)
        .map(Path::to_path_buf)
        .unwrap_or_else(|_| path.to_path_buf());
    artifact_from_file(path, shown.display().to_string())
}

fn artifact_from_file(path: &Path, shown: String) -> AnyResult<Artifact> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Artifact {
        path: shown,
        bytes: bytes.len() as u64,
        sha256: hex::encode(sha256(&bytes)),
    })
}

fn rehash_artifact(stage: &Path, artifact: &Artifact) -> AnyResult<Artifact> {
    let path = Path::new(&artifact.path);
    let path = if path.is_absolute() {
        path.to_path_buf()
    } else {
        stage.join(path)
    };
    artifact_from_file(&path, artifact.path.clone())
}

fn require_hash(name: &str, artifact: &Artifact, expected: &str) -> AnyResult<()> {
    if artifact.sha256 == expected {
        Ok(())
    } else {
        Err(format!(
            "{name} hash mismatch: {} != {expected}",
            artifact.sha256
        ))
    }
}

fn charges_match(left: &Charge, right: &Charge) -> bool {
    left.components == right.components
        && close(left.wall_seconds_sum, right.wall_seconds_sum)
        && close(left.total_core_seconds_sum, right.total_core_seconds_sum)
        && left.peak_rss_bytes_max == right.peak_rss_bytes_max
}

fn close(left: f64, right: f64) -> bool {
    left.is_finite()
        && right.is_finite()
        && (left - right).abs() <= 1e-12 * left.abs().max(right.abs()).max(1.0)
}

fn find_stage_root(path: &Path) -> AnyResult<&Path> {
    path.ancestors()
        .skip(1)
        .find(|candidate| {
            candidate.join("PROTOCOL.md").is_file() && candidate.join("development").is_dir()
        })
        .ok_or_else(|| format!("{} has no Stage 193 root", path.display()))
}

fn validate_commit(name: &str, commit: &str) -> AnyResult<()> {
    if commit.len() == 40 && commit.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        Ok(())
    } else {
        Err(format!("{name} commit must be a full SHA"))
    }
}

fn check(checks: &mut usize, failures: &mut Vec<String>, condition: bool, name: &str) {
    *checks += 1;
    if !condition {
        failures.push(name.into());
    }
}

fn read_json<T: for<'de> Deserialize<'de>>(path: &Path) -> AnyResult<T> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    serde_json::from_slice(&bytes).map_err(|e| format!("{}: {e}", path.display()))
}

fn write_json(path: &Path, value: &impl Serialize) -> AnyResult<()> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    let mut file = File::create(path).map_err(|e| format!("{}: {e}", path.display()))?;
    serde_json::to_writer_pretty(&mut file, value).map_err(|e| e.to_string())?;
    file.write_all(b"\n").map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_cryptominisat_terminal_and_conflicts() {
        let text = "c conflicts                : 12345\ns INDETERMINATE\n";
        let got = parse_cms(text, false);
        assert_eq!(got.0, "timeout_inconclusive");
        assert_eq!(got.1.as_deref(), Some("INDETERMINATE"));
        assert_eq!(got.2, Some(12345));
    }

    #[test]
    fn wdsat_timeout_without_terminal_is_censored() {
        let got = parse_wdsat("", true);
        assert_eq!(got.0, "timeout_inconclusive");
        assert_eq!(got.1, None);
        assert_eq!(got.2, None);
    }

    #[test]
    fn wdsat_unsat_retains_conflicts() {
        let got = parse_wdsat("UNSAT\n231\n", false);
        assert_eq!(got.0, "unsat");
        assert_eq!(got.1.as_deref(), Some("UNSAT"));
        assert_eq!(got.2, Some(231));
    }
}
