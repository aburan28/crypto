//! Data-only replay of the nine original disclosed marked-CMS control roles.
//! A transport pass is narrow: it cannot admit a fresh one-target IC solver.
use super::{json as strict_json, oracle::Curve, sat_control::native, sat_source, target_math};
use crypto_lib::cryptanalysis::prepared_sat_control::sha256;
use serde_json::{json, Value};
use std::{collections::BTreeSet, path::Path, time::Instant};

const INPUT_ROOT: &str = "research/ic_candidate_tournament_20260915/goal_20260924/native-sat-control-registration-v1/result-v1/data";
const BUILD_SCRIPT: &str = "research/ic_candidate_tournament_20260915/goal_20260924/native-sat-target-transport-v1/build-marked-cms.sh";
const PREPARATION_SHA: &str = "67319e5bb0c65e9ae9a10bfbd929487afd5a7e4e8c885b7d123da79b385a9d49";
const ACCEPTED_CMS_SHA: &str = "6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af";
const ACCEPTED_BUNDLE_SHA: &str =
    "a1fd5bd49c80076f3b64fd5cb51d891b278afde765b4d3e39853692ac318bd96";
const SOURCE_SHA: &str = "467b1c3d00a7d6e893332b4d8b42c6326301974d22885aa745f1f926da050323";
const CADICAL_SHA: &str = "8264713f3dc1c4455162d2912238712bd8030fceab0f4b430d106b5a58058d";
const CADIBACK_SHA: &str = "e0aa8f5d67c04527135dde5fe5f943e672d6af11a6f7f924e9f0c07ad3bffba0";
const PATCH_SHA: &str = "b335c52f6673fd6a17d32d0c5506699058b3ee04f0089200305f2c89fe9fde82";
const ACCEPTED_READY: &str = "c Reading from standard input... Use '-h' or '--help' for help.";
const MARKED_READY: &str = "c PREPARED_STDIN_READY_v1";

struct Control {
    point: [u64; 2],
    expected: &'static str,
    files: [&'static str; 4], // manifest, ANF, CNF, Magma
}
const CONTROLS: [Control; 3] = [
    Control {
        point: [40991, 73355],
        expected: "SOURCE_UNSAT",
        files: [
            "d95a7101a3e69759cb1be39363fb305049159c23c5595e0ef82b6511447d7840",
            "a7867583d28fe9593fd17dabdf6545d5bf79b0fda08888e5172f051a3d918eda",
            "05b0913d27231ae82601b44f042765cf9ba2727aeb66543467e92127b3a2d92f",
            "1e7609332d721678d763b1eb48658b505eae39506b0f6e4d0edff8261595eaab",
        ],
    },
    Control {
        point: [73003, 104622],
        expected: "SOURCE_UNSAT",
        files: [
            "58ffc396a551c424e19aa38af3addce803c6aaa7df3de899732d094fd663ca62",
            "eab39a019e5298286e92169efcbe0eb381ce77af2d46749748d0ef202e1478ff",
            "03a017d1763f5be128f8069f78158723a456b0803757eaa03ffbb647d9c6d960",
            "c5f7baa8ef80be4d749fb3be5ac6c3d957978166d534826090455378ec6e5d7e",
        ],
    },
    Control {
        point: [59775, 2910],
        expected: "SAT_MODEL",
        files: [
            "3ecd7c38a81d15b6849e538c89c1673f3e0a979d38e7365f2abb7cd2bfaa9e1e",
            "1783698ff0fda61fc056b24b19d12e421f895c03ebbbc63fd399390b82def031",
            "2077d013463e6a58270a38f29752ec0395099d1db0abfcfaadb7f6b8970d784c",
            "958701bb9e154d2bea03b83e8243aa2b5dfc90278f421ddf3e21ac1d88159a0f",
        ],
    },
];

fn require(ok: bool, why: &str) -> Result<(), String> {
    native::require(ok, why)
}
fn load(path: &Path) -> Result<Value, String> {
    target_math::parse(&native::read(path, 16 * 1024 * 1024)?)
}
fn pinned(path: &Path, max: usize, hash: &str) -> Result<Vec<u8>, String> {
    let bytes = native::read(path, max)?;
    require(
        sha256(&bytes) == hash,
        "disclosed source/build byte hash differs",
    )?;
    Ok(bytes)
}
fn args(file: Option<&Path>) -> Vec<String> {
    let mut out = vec![
        "--verb".into(),
        "1".into(),
        "--threads".into(),
        "1".into(),
        "--random".into(),
        "1".into(),
        "--maxsol".into(),
        "1".into(),
        "--maxconfl".into(),
        "1000000".into(),
    ];
    if let Some(file) = file {
        out.push(file.to_string_lossy().into_owned());
    }
    out
}
fn status(receipt: &Value, stdout: &str) -> Result<&'static str, String> {
    let lines = stdout
        .lines()
        .filter(|line| line.starts_with("s "))
        .collect::<Vec<_>>();
    if receipt["timed_out"] == true {
        return Ok("TIMEOUT");
    }
    match (receipt["exit_code"].as_i64(), lines.as_slice()) {
        (Some(10), ["s SATISFIABLE"]) => Ok("SAT_MODEL"),
        (Some(20), ["s UNSATISFIABLE"]) => Ok("SOURCE_UNSAT"),
        (Some(15), ["s INDETERMINATE"]) => Ok("CONFLICT_BUDGET_INCONCLUSIVE"),
        (Some(0), ["s UNKNOWN"]) => Ok("UNKNOWN_INCONCLUSIVE"),
        _ => Err("disclosed CMS native exit/status is contradictory".into()),
    }
}
fn marker_prefix(bytes: &[u8], marker: &str) -> Result<String, String> {
    let marker = format!("{marker}\n").into_bytes();
    let mut end = 0;
    let mut found = None;
    for line in bytes.split_inclusive(|&byte| byte == b'\n') {
        end += line.len();
        if line == marker {
            require(
                found.replace(end).is_none(),
                "disclosed READY marker repeated",
            )?;
        }
    }
    Ok(sha256(
        &bytes[..found.ok_or("disclosed READY marker absent")?],
    ))
}
fn curve_and_base(preparation: &Value) -> Result<(Curve, Vec<super::oracle::Point>), String> {
    let fixture = json!({"degree":17,"curve_a":1,"subgroup_order":65587,"group_order":131174,
        "cofactor":2,"generator":[43693,23339],"lambda":17184,
        "irreducible":{"degree":17,"low_terms":[0,3]},"targets":[],
        "target_seeds":[],"target_scalar_constructed":false});
    let curve = Curve::new(&strict_json::parse(&fixture.to_string())?)?;
    let base = preparation["record"]["factor_base"]["points"]
        .as_array()
        .ok_or("disclosed preparation lacks geometric points")?
        .iter()
        .map(|value| curve.decode(&strict_json::parse(&value.to_string())?))
        .collect::<Result<Vec<_>, _>>()?;
    require(base.len() == 63, "disclosed geometric base count differs")?;
    Ok((curve, base))
}
fn has_three_sum(
    curve: &Curve,
    base: &[super::oracle::Point],
    target: super::oracle::Point,
) -> bool {
    base.iter().any(|&a| {
        base.iter().any(|&b| {
            base.iter()
                .any(|&c| curve.add(curve.add(a, b), c) == target)
        })
    })
}
fn model_witness(
    model: &[bool],
    row: &Value,
    curve: &Curve,
    base: &[super::oracle::Point],
    target: super::oracle::Point,
) -> Result<(), String> {
    let indices: [usize; 3] =
        serde_json::from_value(row["model_check"]["full_point_witness_indices"].clone())
            .map_err(|e| e.to_string())?;
    require(model.len() >= 51, "disclosed CMS model is short")?;
    let xs: [u64; 3] = std::array::from_fn(|i| {
        (0..6).fold(0, |value, bit| value | ((model[i * 6 + bit] as u64) << bit))
    });
    require(
        indices.iter().enumerate().all(|(i, &index)| {
            base.get(index)
                .is_some_and(|p| p.is_some_and(|p| p.0 == xs[i]))
        }) && indices
            .iter()
            .fold(None, |sum, &index| curve.add(sum, base[index]))
            == target,
        "disclosed SAT model witness fails independent x/full-point replay",
    )
}
fn build_pin(root: &Path, build: &Path, accepted: &Path) -> Result<String, String> {
    let inputs = load(&build.join("build-inputs.json"))?;
    let receipt = load(&build.join("build-receipt.json"))?;
    let terminal = load(&build.join("terminal.json"))?;
    require(
        inputs["accepted_asset_bundle_sha256"] == ACCEPTED_BUNDLE_SHA
            && inputs["source_archive_sha256"] == SOURCE_SHA
            && inputs["cadical_archive_sha256"] == CADICAL_SHA
            && inputs["cadiback_archive_sha256"] == CADIBACK_SHA
            && inputs["postbuffer_patch_sha256"] == PATCH_SHA
            && inputs["build_script_sha256"]
                == sha256(&native::read(&root.join(BUILD_SCRIPT), 65536)?)
            && terminal["status"] == "BUILT_UNVALIDATED"
            && terminal["exit_code"] == 0
            && receipt["status"] == "BUILT_UNVALIDATED"
            && receipt["transport_parity_passed"] == false,
        "marked CMS source/build receipt differs",
    )?;
    for (name, hash, max) in [
        ("source.tar", SOURCE_SHA, 8 * 1024 * 1024),
        ("cadical.tar", CADICAL_SHA, 8 * 1024 * 1024),
        ("cadiback.tar", CADIBACK_SHA, 8 * 1024 * 1024),
        ("postbuffer.patch", PATCH_SHA, 65536),
    ] {
        pinned(&build.join(name), max, hash)?;
    }
    for (name, hash) in [
        (
            "src/main.cpp",
            "b6c65e963442cc63df10b0d4888737693aeebccb94b6d2379b3fb2bef37160b2",
        ),
        (
            "src/dimacsparser.h",
            "a3919c1761e77243ed526c3be38b620f111da071a7ee098595e88ab3980a2355",
        ),
        (
            "src/streambuffer.h",
            "e58dce81a60271641ab6f2c8923e3a1c02d945a75463c9d666d791ec0a077c0a",
        ),
    ] {
        pinned(&build.join("source").join(name), 2 * 1024 * 1024, hash)?;
    }
    for (name, field) in [
        ("configure.stdout", "configure_stdout"),
        ("configure.stderr", "configure_stderr"),
        ("build.stdout", "build_stdout"),
        ("build.stderr", "build_stderr"),
    ] {
        require(
            receipt["logs_sha256"][field]
                == sha256(&native::read(&build.join(name), 16 * 1024 * 1024)?),
            "marked CMS retained build log differs",
        )?;
    }
    let pin = receipt["prepared_cms_sha256"]
        .as_str()
        .ok_or("marked CMS pin absent")?
        .to_string();
    require(
        pin.len() == 64
            && pin.bytes().all(|b| b.is_ascii_hexdigit())
            && sha256(&native::read(
                &build.join("bin/prepared-cms"),
                16 * 1024 * 1024,
            )?) == pin
            && sha256(&native::read(accepted, 16 * 1024 * 1024)?) == ACCEPTED_CMS_SHA,
        "disclosed accepted/marked executable differs",
    )?;
    Ok(pin)
}
fn role(
    dir: &Path,
    name: &str,
    binary: &Path,
    pin: &str,
    cnf: &[u8],
    row: &Value,
    curve: &Curve,
    base: &[super::oracle::Point],
    target: super::oracle::Point,
    anf: &str,
) -> Result<&'static str, String> {
    let receipt = load(&dir.join(format!("{name}.receipt.json")))?;
    let stdout = native::read(&dir.join(format!("{name}.stdout")), 8 * 1024 * 1024)?;
    let stderr = native::read(&dir.join(format!("{name}.stderr")), 8 * 1024 * 1024)?;
    let text = std::str::from_utf8(&stdout).map_err(|e| e.to_string())?;
    let file = (name == "accepted-file").then(|| dir.join("file-input.xor.cnf"));
    let argv = std::iter::once(binary.to_string_lossy().into_owned())
        .chain(args(file.as_deref()))
        .collect::<Vec<_>>();
    let pid = receipt["pid"]
        .as_u64()
        .filter(|pid| *pid > 1)
        .ok_or("disclosed role PID absent")?;
    require(
        row["receipt"] == receipt
            && row["process_group_drain_passed"] == true
            && receipt["argv"] == json!(argv)
            && receipt["cwd"] == json!(dir)
            && receipt["environment"] == json!({"LC_ALL":"C"})
            && receipt["deadline_ms"] == 120000
            && receipt["executable_sha256_before"] == pin
            && receipt["executable_sha256_after"] == pin
            && receipt["stdout_bytes"] == stdout.len()
            && receipt["stdout_sha256"] == sha256(&stdout)
            && receipt["stderr_bytes"] == stderr.len()
            && receipt["stderr_sha256"] == sha256(&stderr)
            && receipt["process_group_drain_confirmed"] == true
            && native::read(&dir.join(format!("{name}-pids")), 65536)?
                == format!("start {pid}\ndone {pid}\n").as_bytes(),
        "original disclosed role identity/output/PID drain differs",
    )?;
    if let Some(file) = file {
        require(
            native::read(&file, 16 * 1024 * 1024)? == cnf && receipt["stdin_bytes"].is_null(),
            "disclosed file-mode CNF or stdin differs",
        )?;
    } else {
        let ready = load(&dir.join(format!("{name}.ready.json")))?;
        let marker = if name == "marked-stdin" {
            MARKED_READY
        } else {
            ACCEPTED_READY
        };
        require(
            ready["state"] == "ready-without-input"
                && ready["pid"] == pid
                && ready["argv"] == json!(argv)
                && ready["cwd"] == json!(dir)
                && ready["environment"] == json!({"LC_ALL":"C"})
                && ready["executable_sha256"] == pin
                && ready["marker"] == marker
                && ready["stdin_written_bytes"] == 0
                && ready["stdout_at_ready_sha256"] == marker_prefix(&stdout, marker)?
                && receipt["state"] == "used"
                && receipt["marker"] == marker
                && receipt["readiness_marker_count"] == 1
                && receipt["launch_to_ready_ns"] == ready["ready_after_launch_ns"]
                && receipt["stdin_bytes"] == cnf.len()
                && receipt["stdin_sha256"] == sha256(cnf)
                && receipt["stdin_written_bytes"] == cnf.len()
                && receipt["stdin_write_error"].is_null(),
            "original disclosed zero-input READY/CNF delivery differs",
        )?;
    }
    let observed = status(&receipt, text)?;
    require(
        row["status"] == observed
            && row["exit_code"] == receipt["exit_code"]
            && row["timed_out"] == receipt["timed_out"],
        "disclosed role status differs from original native bytes",
    )?;
    if observed == "SAT_MODEL" {
        let model = sat_source::verify_native_model(
            anf,
            std::str::from_utf8(cnf).map_err(|e| e.to_string())?,
            text,
        )?;
        require(
            row["model_error"].is_null()
                && row["model_check"]["source_model_valid"] == true
                && row["model_check"]["model_sha256"]
                    == sha256(&model.iter().map(|&b| u8::from(b)).collect::<Vec<_>>()),
            "disclosed source-valid model digest differs",
        )?;
        model_witness(&model, row, curve, base, target)?;
    } else {
        require(
            row["model_check"].is_null() && row["model_error"].is_null(),
            "non-model disclosed role claims a source model",
        )?;
    }
    Ok(observed)
}
fn verify(root: &Path, accepted: &Path, build: &Path, execution: &Path) -> Result<Value, String> {
    let marked_pin = build_pin(root, build, accepted)?;
    let inventory = native::inventory(execution)?;
    let observed = inventory
        .as_object()
        .ok_or("disclosed execution inventory is malformed")?
        .keys()
        .cloned()
        .collect::<BTreeSet<_>>();
    let mut expected = BTreeSet::from(["terminal.json".to_string()]);
    for trial in 0..CONTROLS.len() {
        let dir = format!("query-{trial:02}");
        for name in ["result.json", "file-input.xor.cnf"] {
            expected.insert(format!("{dir}/{name}"));
        }
        for role in ["accepted-file", "accepted-stdin", "marked-stdin"] {
            for suffix in ["stdout", "stderr", "receipt.json"] {
                expected.insert(format!("{dir}/{role}.{suffix}"));
            }
            expected.insert(format!("{dir}/{role}-pids"));
            if role != "accepted-file" {
                expected.insert(format!("{dir}/{role}.ready.json"));
            }
        }
    }
    require(
        observed == expected,
        "disclosed original execution has missing or extra role files",
    )?;
    let prep = target_math::parse(&pinned(
        &root.join(INPUT_ROOT).join("preparation.json"),
        16 * 1024 * 1024,
        PREPARATION_SHA,
    )?)?;
    let (curve, base) = curve_and_base(&prep)?;
    let terminal = load(&execution.join("terminal.json"))?;
    let rows = terminal["rows"]
        .as_array()
        .ok_or("disclosed control terminal lacks rows")?;
    require(
        rows.len() == 3
            && terminal["question"] == "marked-cms-three-way-disclosed-control-v1"
            && terminal["source_bound_target_admitted"] == false
            && terminal["independent_data_only_replay_passed"] == false
            && terminal["natural_relation_yield"].is_null()
            && terminal["online_speedup"].is_null()
            && terminal["accepted_cms_sha256"] == ACCEPTED_CMS_SHA
            && terminal["marked_cms_sha256"] == marked_pin,
        "disclosed control terminal makes an unadmitted claim",
    )?;
    let mut outcomes = Vec::new();
    let mut pids = BTreeSet::new();
    for (trial, control) in CONTROLS.iter().enumerate() {
        let dir = execution.join(format!("query-{trial:02}"));
        let row = load(&dir.join("result.json"))?;
        require(
            row == rows[trial]
                && row["trial"] == trial
                && row["point"] == json!(control.point)
                && row["expected_status"] == control.expected
                && row["manifest_sha256"] == control.files[0]
                && row["anf_sha256"] == control.files[1]
                && row["cnf_sha256"] == control.files[2]
                && row["magma_sha256"] == control.files[3]
                && row["source_bound_target_admitted"] == false,
            "disclosed row differs from fixed input or terminal",
        )?;
        let src = root
            .join(INPUT_ROOT)
            .join(format!("execution/query-{trial:02}/instance"));
        let manifest = target_math::parse(&pinned(
            &src.join("manifest.json"),
            1_048_576,
            control.files[0],
        )?)?;
        require(
            manifest["target"]["x"] == control.point[0].to_string()
                && manifest["target"]["y"] == control.point[1].to_string(),
            "disclosed source manifest point differs",
        )?;
        let anf = String::from_utf8(pinned(
            &src.join("instance.anf"),
            16 * 1024 * 1024,
            control.files[1],
        )?)
        .map_err(|e| e.to_string())?;
        let cnf = pinned(
            &src.join("instance.xor.cnf"),
            16 * 1024 * 1024,
            control.files[2],
        )?;
        pinned(
            &src.join("instance.magma"),
            2 * 1024 * 1024,
            control.files[3],
        )?;
        let target = curve.decode(&strict_json::parse(&json!(control.point).to_string())?)?;
        require(
            target.is_some()
                && has_three_sum(&curve, &base, target) == (control.expected == "SAT_MODEL"),
            "disclosed full-point geometric class differs",
        )?;
        if control.expected == "SAT_MODEL" {
            let indices: [usize; 3] = serde_json::from_value(row["geometric_witness"].clone())
                .map_err(|e| e.to_string())?;
            require(
                indices.iter().all(|&index| index < base.len())
                    && indices
                        .iter()
                        .fold(None, |sum, &index| curve.add(sum, base[index]))
                        == target,
                "disclosed geometric control witness differs",
            )?;
        } else {
            require(
                row["geometric_witness"].is_null(),
                "geometric negative row claims a witness",
            )?;
        }
        let file = role(
            &dir,
            "accepted-file",
            accepted,
            ACCEPTED_CMS_SHA,
            &cnf,
            &row["accepted_file"],
            &curve,
            &base,
            target,
            &anf,
        )?;
        let stdin = role(
            &dir,
            "accepted-stdin",
            accepted,
            ACCEPTED_CMS_SHA,
            &cnf,
            &row["accepted_stdin"],
            &curve,
            &base,
            target,
            &anf,
        )?;
        let marked = role(
            &dir,
            "marked-stdin",
            &build.join("bin/prepared-cms"),
            &marked_pin,
            &cnf,
            &row["marked_stdin"],
            &curve,
            &base,
            target,
            &anf,
        )?;
        for name in ["accepted_file", "accepted_stdin", "marked_stdin"] {
            let pid = row[name]["receipt"]["pid"]
                .as_u64()
                .ok_or("disclosed PID missing")?;
            require(pids.insert(pid), "disclosed role PID reused")?;
        }
        let pass = [file, stdin, marked].iter().all(|&v| v == control.expected)
            && row["source_unchanged"] == true
            && row["binaries_unchanged"] == true;
        require(
            row["three_way_transport_pass"] == pass,
            "disclosed producer parity verdict differs from original roles",
        )?;
        outcomes.push(
            json!({"trial":trial,"point":control.point,"statuses":[file,stdin,marked],
            "three_way_transport_pass":pass}),
        );
    }
    let pass = outcomes
        .iter()
        .all(|row| row["three_way_transport_pass"] == true);
    require(
        terminal["status"]
            == if pass {
                "THREE_WAY_PARITY_REPORTED_REPLAY_PENDING"
            } else {
                "TRANSPORT_PARITY_FAILED"
            },
        "disclosed terminal parity status differs",
    )?;
    Ok(
        json!({"schema_version":1,"status":if pass {"PASS_DISCLOSED_CMS_TRANSPORT_PARITY"} else {"AUDITED_DISCLOSED_CMS_TRANSPORT_FAILURE"},
        "audited_roles":9,"rows":outcomes,"source_binary_pins_checked":true,
        "disclosed_transport_parity_admitted":pass,"source_bound_target_admitted":false,
        "fresh_paired_qualification":false,"natural_relation_yield":null,"online_speedup":null,
        "native_solver_calls_by_auditor":0}),
    )
}
pub(super) fn run(
    root: &Path,
    accepted: &Path,
    build: &Path,
    execution: &Path,
    out: &Path,
) -> Result<String, String> {
    let started = Instant::now();
    let root = root.canonicalize().map_err(|e| e.to_string())?;
    let accepted = accepted.canonicalize().map_err(|e| e.to_string())?;
    let build = build.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    let parent = out
        .parent()
        .ok_or("disclosed audit output lacks parent")?
        .canonicalize()
        .map_err(|e| e.to_string())?;
    require(
        !parent.starts_with(&root)
            && !parent.starts_with(&build)
            && !parent.starts_with(&execution),
        "disclosed audit output would alter original source/build/execution",
    )?;
    let before = native::inventory(&execution)?;
    let result = verify(&root, &accepted, &build, &execution);
    require(
        native::inventory(&execution)? == before,
        "disclosed execution changed during data-only replay",
    )?;
    let error = result.as_ref().err().cloned();
    let mut receipt=result.unwrap_or_else(|reason|json!({"schema_version":1,"status":"DISCLOSED_CMS_TRANSPORT_NOT_VERIFIED",
        "error":reason,"disclosed_transport_parity_admitted":false,"source_bound_target_admitted":false,
        "natural_relation_yield":null,"online_speedup":null,"native_solver_calls_by_auditor":0}));
    receipt["audit_wall_ns"] = json!(started.elapsed().as_nanos().to_string());
    native::save(out, &receipt)?;
    if let Some(reason) = error {
        return Err(format!(
            "disclosed CMS data-only audit failed; receipt retained: {reason}"
        ));
    }
    serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strict_native_status_rejects_contradictions() {
        assert_eq!(
            status(
                &json!({"timed_out":false,"exit_code":20}),
                "s UNSATISFIABLE\n"
            )
            .unwrap(),
            "SOURCE_UNSAT"
        );
        assert!(status(
            &json!({"timed_out":false,"exit_code":20}),
            "s SATISFIABLE\n"
        )
        .is_err());
        assert!(status(
            &json!({"timed_out":false,"exit_code":10}),
            "s SATISFIABLE\ns UNSATISFIABLE\n"
        )
        .is_err());
    }
}
