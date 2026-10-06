//! Disclosed three-way CMS transport control on three closed, fixed instances.
//! The accepted file/stdin roles and new post-buffer-marker role each run once.
//! No target is fresh, no yield is inferred, and no online speed is measured.
#[path = "../src/bin/prepared_sat_worker/native.rs"]
#[allow(dead_code)]
mod native;

use crypto_lib::cryptanalysis::{
    koblitz_fast::FastPoint,
    prepared_sat_control::{
        native_status, parse_model, sha256, validate_manifest, validate_source_syntax, verify_anf,
        verify_cnf, NativeOutput, PreparedState, ORDER,
    },
};
use serde_json::{json, Value};
use std::{
    collections::HashMap,
    env, fs,
    path::{Path, PathBuf},
    process::ExitCode,
};

const INPUT_ROOT: &str = "research/ic_candidate_tournament_20260915/goal_20260924/native-sat-control-registration-v1/result-v1/data";
const PREPARATION_SHA: &str = "67319e5bb0c65e9ae9a10bfbd929487afd5a7e4e8c885b7d123da79b385a9d49";
const ACCEPTED_CMS_SHA: &str = "6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af";
const ACCEPTED_BUNDLE_SHA: &str =
    "a1fd5bd49c80076f3b64fd5cb51d891b278afde765b4d3e39853692ac318bd96";
const SOURCE_SHA: &str = "467b1c3d00a7d6e893332b4d8b42c6326301974d22885aa745f1f926da050323";
const CADICAL_SHA: &str = "8264713f3dc1c4455162d2912238712bd8030fceab0f4b430d106b5a58058d";
const CADIBACK_SHA: &str = "e0aa8f5d67c04527135dde5fe5f943e672d6af11a6f7f924e9f0c07ad3bffba0";
const PATCH_SHA: &str = "b335c52f6673fd6a17d32d0c5506699058b3ee04f0089200305f2c89fe9fde82";
const BUILD_SCRIPT: &str = "research/ic_candidate_tournament_20260915/goal_20260924/native-sat-target-transport-v1/build-marked-cms.sh";
const ACCEPTED_READY: &str = "c Reading from standard input... Use '-h' or '--help' for help.";
const MARKED_READY: &str = "c PREPARED_STDIN_READY_v1";
const DEADLINE_MS: u64 = 120_000;
const CONFLICTS: u64 = 1_000_000;

struct Control {
    point: [u64; 2],
    expected: &'static str,
    manifest: &'static str,
    anf: &'static str,
    cnf: &'static str,
    magma: &'static str,
}
const CONTROLS: [Control; 3] = [
    Control {
        point: [40991, 73355],
        expected: "SOURCE_UNSAT",
        manifest: "d95a7101a3e69759cb1be39363fb305049159c23c5595e0ef82b6511447d7840",
        anf: "a7867583d28fe9593fd17dabdf6545d5bf79b0fda08888e5172f051a3d918eda",
        cnf: "05b0913d27231ae82601b44f042765cf9ba2727aeb66543467e92127b3a2d92f",
        magma: "1e7609332d721678d763b1eb48658b505eae39506b0f6e4d0edff8261595eaab",
    },
    Control {
        point: [73003, 104622],
        expected: "SOURCE_UNSAT",
        manifest: "58ffc396a551c424e19aa38af3addce803c6aaa7df3de899732d094fd663ca62",
        anf: "eab39a019e5298286e92169efcbe0eb381ce77af2d46749748d0ef202e1478ff",
        cnf: "03a017d1763f5be128f8069f78158723a456b0803757eaa03ffbb647d9c6d960",
        magma: "c5f7baa8ef80be4d749fb3be5ac6c3d957978166d534826090455378ec6e5d7e",
    },
    Control {
        point: [59775, 2910],
        expected: "SAT_MODEL",
        manifest: "3ecd7c38a81d15b6849e538c89c1673f3e0a979d38e7365f2abb7cd2bfaa9e1e",
        anf: "1783698ff0fda61fc056b24b19d12e421f895c03ebbbc63fd399390b82def031",
        cnf: "2077d013463e6a58270a38f29752ec0395099d1db0abfcfaadb7f6b8970d784c",
        magma: "958701bb9e154d2bea03b83e8243aa2b5dfc90278f421ddf3e21ac1d88159a0f",
    },
];

fn require(ok: bool, why: &str) -> Result<(), String> {
    native::require(ok, why)
}
fn pinned(path: &Path, limit: usize, expected: &str) -> Result<Vec<u8>, String> {
    let bytes = native::read(path, limit)?;
    require(
        sha256(&bytes) == expected,
        "disclosed control input digest differs",
    )?;
    Ok(bytes)
}
fn cms_args(file: Option<&Path>) -> Vec<String> {
    let mut args = vec![
        "--verb".into(),
        "1".into(),
        "--threads".into(),
        "1".into(),
        "--random".into(),
        "1".into(),
        "--maxsol".into(),
        "1".into(),
        "--maxconfl".into(),
        CONFLICTS.to_string(),
    ];
    if let Some(file) = file {
        args.push(file.to_string_lossy().into_owned());
    }
    args
}
fn geometric_witness(state: &PreparedState, point: FastPoint) -> Option<[usize; 3]> {
    let mut pairs = HashMap::new();
    for (i, &a) in state.geometry.iter().enumerate() {
        for (j, &b) in state.geometry.iter().enumerate() {
            pairs.entry(state.curve.add(a, b)).or_insert((i, j));
        }
    }
    for (k, &c) in state.geometry.iter().enumerate() {
        let wanted = state.curve.add(point, state.curve.neg(c));
        if let Some(&(i, j)) = pairs.get(&wanted) {
            if state
                .curve
                .add(state.curve.add(state.geometry[i], state.geometry[j]), c)
                == point
            {
                return Some([i, j, k]);
            }
        }
    }
    None
}
fn model_check(
    output: &NativeOutput,
    anf: &str,
    cnf: &str,
    state: &PreparedState,
    point: FastPoint,
) -> Result<Value, String> {
    if native_status(output) != "SAT_MODEL" {
        return Ok(Value::Null);
    }
    let width = cnf
        .lines()
        .find(|line| line.starts_with("p cnf "))
        .and_then(|line| line.split_whitespace().nth(2))
        .and_then(|value| value.parse::<usize>().ok())
        .ok_or("disclosed CNF variable count is missing")?;
    let model = parse_model(&output.stdout, width)?;
    verify_anf(anf, model.get(..51).ok_or("short disclosed SAT model")?)?;
    verify_cnf(cnf, &model)?;
    let indices = state
        .lift(&model, point)
        .ok_or("source-valid disclosed SAT model has no full-point lift")?;
    require(
        indices.iter().fold(FastPoint::INFINITY, |sum, &index| {
            state.curve.add(sum, state.geometry[index])
        }) == point,
        "disclosed SAT witness does not readd to the public point",
    )?;
    Ok(json!({"source_model_valid":true,"model_sha256":sha256(
        &model.iter().map(|&bit| u8::from(bit)).collect::<Vec<_>>()),
        "full_point_witness_indices":indices}))
}
fn run_file(
    cms: &Path,
    dir: &Path,
    cnf: &[u8],
    env: &[(String, String)],
) -> Result<NativeOutput, String> {
    let copied = dir.join("file-input.xor.cnf");
    native::create(&copied, cnf)?;
    let args = cms_args(Some(&copied));
    native::measured_child_request(native::ChildRequest {
        program: cms,
        args: &args,
        cwd: dir,
        stem: &dir.join("accepted-file"),
        deadline_ms: DEADLINE_MS,
        ledger: &dir.join("accepted-file-pids"),
        helper: false,
        input: None,
        environment: env,
    })
}
fn run_stdin(
    cms: &Path,
    pin: &str,
    dir: &Path,
    role: &str,
    marker: &str,
    cnf: &[u8],
    env: &[(String, String)],
) -> Result<NativeOutput, String> {
    let args = cms_args(None);
    native::PreparedChild::launch(native::PreparedLaunch {
        program: cms,
        expected_sha256: pin,
        args: &args,
        cwd: dir,
        stem: &dir.join(role),
        ledger: &dir.join(format!("{role}-pids")),
        environment: env,
        marker,
        ready_deadline_ms: 30_000,
    })?
    .deliver(cnf, DEADLINE_MS)
}
fn role_result(
    dir: &Path,
    role: &str,
    result: Result<NativeOutput, String>,
    anf: &str,
    cnf: &str,
    state: &PreparedState,
    point: FastPoint,
) -> Value {
    let ledger = native::drain_ledger(&dir.join(format!("{role}-pids")));
    match result {
        Ok(output) => {
            let status = native_status(&output);
            let model = model_check(&output, anf, cnf, state, point);
            json!({"role":role,"status":status,"exit_code":output.exit_code,
                "timed_out":output.timed_out,"receipt":output.receipt,
                "model_check":model.as_ref().ok(),"model_error":model.as_ref().err(),
                "process_group_drain_passed":ledger.is_ok(),"drain_error":ledger.err(),
                "source_bound_target_admitted":false})
        }
        Err(error) => json!({"role":role,"error":error,"process_group_drain_passed":ledger.is_ok(),
            "drain_error":ledger.err(),"source_bound_target_admitted":false}),
    }
}
fn input(root: &Path, trial: usize, control: &Control) -> Result<(String, String), String> {
    let dir = root
        .join(INPUT_ROOT)
        .join(format!("execution/query-{trial:02}/instance"));
    let manifest_bytes = pinned(&dir.join("manifest.json"), 1_048_576, control.manifest)?;
    let manifest: Value = serde_json::from_slice(&manifest_bytes).map_err(|e| e.to_string())?;
    validate_manifest(&manifest, control.point)?;
    let anf = pinned(&dir.join("instance.anf"), 16 * 1024 * 1024, control.anf)?;
    let cnf = pinned(&dir.join("instance.xor.cnf"), 16 * 1024 * 1024, control.cnf)?;
    pinned(&dir.join("instance.magma"), 2 * 1024 * 1024, control.magma)?;
    let anf = String::from_utf8(anf).map_err(|e| e.to_string())?;
    let cnf = String::from_utf8(cnf).map_err(|e| e.to_string())?;
    let width = manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
        .as_u64()
        .ok_or("disclosed source variable count missing")? as usize;
    validate_source_syntax(&anf, &cnf, width)?;
    Ok((anf, cnf))
}
fn run(root: &Path, accepted: &Path, marked_build: &Path, out: &Path) -> Result<Value, String> {
    native::enforce_hardware()?;
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let accepted = accepted.canonicalize().map_err(|e| e.to_string())?;
    let marked_build = marked_build.canonicalize().map_err(|e| e.to_string())?;
    let marked = marked_build.join("bin/prepared-cms");
    let build = native::load(&marked_build.join("build-receipt.json"))?;
    let terminal = native::load(&marked_build.join("terminal.json"))?;
    let inputs = native::load(&marked_build.join("build-inputs.json"))?;
    require(
        build["status"] == "BUILT_UNVALIDATED"
            && terminal["status"] == "BUILT_UNVALIDATED"
            && terminal["exit_code"] == 0
            && build["transport_parity_passed"] == false
            && inputs["accepted_asset_bundle_sha256"] == ACCEPTED_BUNDLE_SHA
            && inputs["source_archive_sha256"] == SOURCE_SHA
            && inputs["cadical_archive_sha256"] == CADICAL_SHA
            && inputs["cadiback_archive_sha256"] == CADIBACK_SHA
            && inputs["postbuffer_patch_sha256"] == PATCH_SHA
            && inputs["build_script_sha256"]
                == sha256(&native::read(&root.join(BUILD_SCRIPT), 64 * 1024)?)
            && sha256(&native::read(
                &marked_build.join("source.tar"),
                8 * 1024 * 1024,
            )?) == SOURCE_SHA
            && sha256(&native::read(
                &marked_build.join("cadical.tar"),
                8 * 1024 * 1024,
            )?) == CADICAL_SHA
            && sha256(&native::read(
                &marked_build.join("cadiback.tar"),
                8 * 1024 * 1024,
            )?) == CADIBACK_SHA
            && sha256(&native::read(
                &marked_build.join("postbuffer.patch"),
                64 * 1024,
            )?) == PATCH_SHA
            && build["prepared_cms_sha256"] == sha256(&native::read(&marked, 16 * 1024 * 1024)?)
            && sha256(&native::read(&accepted, 16 * 1024 * 1024)?) == ACCEPTED_CMS_SHA,
        "disclosed CMS binary or build receipt differs",
    )?;
    let marked_pin = build["prepared_cms_sha256"]
        .as_str()
        .ok_or("marked CMS binary pin missing")?;
    let preparation = pinned(
        &root.join(INPUT_ROOT).join("preparation.json"),
        16 * 1024 * 1024,
        PREPARATION_SHA,
    )?;
    let preparation: Value = serde_json::from_slice(&preparation).map_err(|e| e.to_string())?;
    let state = PreparedState::load(&preparation)?;
    let mut sources = Vec::new();
    for (trial, control) in CONTROLS.iter().enumerate() {
        let (anf, cnf) = input(&root, trial, control)?;
        let point = FastPoint::affine(control.point[0], control.point[1]);
        require(
            state.curve.is_on_curve(point) && state.curve.mul_u64(point, ORDER).infinity,
            "disclosed public point is outside the declared subgroup",
        )?;
        let geometric = geometric_witness(&state, point);
        require(
            geometric.is_some() == (control.expected == "SAT_MODEL"),
            "disclosed geometric class differs before transport roles",
        )?;
        sources.push((anf, cnf, geometric));
    }
    let output_parent = out
        .parent()
        .ok_or("disclosed control output lacks parent")?
        .canonicalize()
        .map_err(|e| e.to_string())?;
    require(
        !output_parent.starts_with(&root) && !output_parent.starts_with(&marked_build),
        "disclosed control output would alter original source or build evidence",
    )?;
    fs::create_dir(out).map_err(|e| e.to_string())?;
    let out = out.canonicalize().map_err(|e| e.to_string())?;
    let env = vec![("LC_ALL".to_string(), "C".to_string())];
    let mut rows = Vec::new();
    for (trial, control) in CONTROLS.iter().enumerate() {
        let dir = out.join(format!("query-{trial:02}"));
        fs::create_dir(&dir).map_err(|e| e.to_string())?;
        let (anf, cnf, geometric) = &sources[trial];
        let point = FastPoint::affine(control.point[0], control.point[1]);
        let file = role_result(
            &dir,
            "accepted-file",
            run_file(&accepted, &dir, cnf.as_bytes(), &env),
            anf,
            cnf,
            &state,
            point,
        );
        let accepted_stdin = role_result(
            &dir,
            "accepted-stdin",
            run_stdin(
                &accepted,
                ACCEPTED_CMS_SHA,
                &dir,
                "accepted-stdin",
                ACCEPTED_READY,
                cnf.as_bytes(),
                &env,
            ),
            anf,
            cnf,
            &state,
            point,
        );
        let marked_stdin = role_result(
            &dir,
            "marked-stdin",
            run_stdin(
                &marked,
                marked_pin,
                &dir,
                "marked-stdin",
                MARKED_READY,
                cnf.as_bytes(),
                &env,
            ),
            anf,
            cnf,
            &state,
            point,
        );
        let source_unchanged = input(&root, trial, control).is_ok_and(|(after_anf, after_cnf)| {
            after_anf == anf.as_str() && after_cnf == cnf.as_str()
        });
        let binaries_unchanged = native::read(&accepted, 16 * 1024 * 1024)
            .is_ok_and(|bytes| sha256(&bytes) == ACCEPTED_CMS_SHA)
            && native::read(&marked, 16 * 1024 * 1024)
                .is_ok_and(|bytes| sha256(&bytes) == marked_pin);
        let roles = [&file, &accepted_stdin, &marked_stdin];
        let pass = source_unchanged
            && binaries_unchanged
            && roles.iter().all(|row| {
                row["status"] == control.expected
                    && row["timed_out"] == false
                    && row["process_group_drain_passed"] == true
                    && row["model_error"].is_null()
                    && (control.expected != "SAT_MODEL"
                        || row["model_check"]["source_model_valid"] == true)
            });
        let row = json!({"trial":trial,"point":control.point,"expected_status":control.expected,
            "manifest_sha256":control.manifest,"anf_sha256":control.anf,"cnf_sha256":control.cnf,
            "magma_sha256":control.magma,"geometric_witness":geometric,
            "source_unchanged":source_unchanged,"binaries_unchanged":binaries_unchanged,
            "accepted_file":file,"accepted_stdin":accepted_stdin,"marked_stdin":marked_stdin,
            "three_way_transport_pass":pass,"source_bound_target_admitted":false});
        native::save(&dir.join("result.json"), &row)?;
        rows.push(row);
    }
    let pass = rows
        .iter()
        .all(|row| row["three_way_transport_pass"] == true);
    require(
        sha256(&native::read(&accepted, 16 * 1024 * 1024)?) == ACCEPTED_CMS_SHA
            && sha256(&native::read(&marked, 16 * 1024 * 1024)?) == marked_pin,
        "CMS binary changed during disclosed transport control",
    )?;
    Ok(
        json!({"schema_version":1,"question":"marked-cms-three-way-disclosed-control-v1",
        "status":if pass {"THREE_WAY_PARITY_REPORTED_REPLAY_PENDING"} else {"TRANSPORT_PARITY_FAILED"},
        "rows":rows,"accepted_cms_sha256":ACCEPTED_CMS_SHA,"marked_cms_sha256":marked_pin,
        "disclosed_points_only":true,"natural_relation_yield":null,"target_count":0,
        "independent_data_only_replay_passed":false,"source_bound_target_admitted":false,
        "fresh_paired_qualification":false,"online_speedup":null}),
    )
}
fn main() -> ExitCode {
    let args = env::args_os().collect::<Vec<_>>();
    if args.len() != 5 {
        eprintln!(
            "usage: ic_marked_cms_control SOURCE_REPO ACCEPTED_CMS MARKED_BUILD NEW_OUTPUT_DIR"
        );
        return ExitCode::FAILURE;
    }
    let paths = args[1..].iter().map(PathBuf::from).collect::<Vec<_>>();
    if !paths.iter().all(|path| path.is_absolute()) {
        eprintln!("all disclosed control paths must be absolute");
        return ExitCode::FAILURE;
    }
    let out = &paths[3];
    let result = run(&paths[0], &paths[1], &paths[2], out);
    if out.is_dir() {
        let report = result.as_ref().cloned().unwrap_or_else(|error| {
            json!({"schema_version":1,"question":"marked-cms-three-way-disclosed-control-v1",
                "status":"CONTROL_INTERRUPTED_OR_FAILED","error":error,"natural_relation_yield":null,
                "independent_data_only_replay_passed":false,"source_bound_target_admitted":false,
                "online_speedup":null})
        });
        if let Err(error) = native::save(&out.join("terminal.json"), &report) {
            eprintln!("cannot retain disclosed control terminal: {error}");
            return ExitCode::FAILURE;
        }
    }
    match result {
        Ok(report) if report["status"] == "THREE_WAY_PARITY_REPORTED_REPLAY_PENDING" => {
            ExitCode::SUCCESS
        }
        Ok(_) => ExitCode::FAILURE,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}
