//! Supplement ecbench's standard audit, which selects verified records only.
//! Re-execute failures using the same native ecbench child and exact inputs.
use crypto_lib::cryptanalysis::ecbench::{
    canonical::sha256_hex,
    methods::{MethodSpec, PhaseRecord, SolveReport},
    record::{json_roundtrip, ChildInput, ChildOutput, Record, CHILD_INPUT_SCHEMA},
    runner::read_records,
    spec::{plan, Spec},
    workload::CurveSpec,
};
use serde_json::{json, Value};
use std::{
    fs,
    io::Write,
    path::Path,
    process::{Command, Stdio},
};

fn phase_projection(phases: &[PhaseRecord]) -> Value {
    json!(phases.iter().map(|p|json!({"name":p.name,"adds":p.adds,"doubles":p.doubles,
        "scalar_mults":p.scalar_mults,"native":p.native,"gae_bits":json_roundtrip(p.gae).to_bits()})).collect::<Vec<_>>())
}

fn differences(r: &Record, rep: &SolveReport) -> Vec<&'static str> {
    let mut diffs = Vec::new();
    if rep.recovered.map(|n| n.to_string()) != r.outcome.recovered {
        diffs.push("recovered");
    }
    if rep.exhausted != (r.outcome.status == "exhausted") {
        diffs.push("exhausted");
    }
    if Some(json_roundtrip(rep.total_gae).to_bits())
        != r.cost.total_gae.map(|n| json_roundtrip(n).to_bits())
    {
        diffs.push("total_gae");
    }
    if phase_projection(&rep.phases) != phase_projection(&r.phases) {
        diffs.push("phases");
    }
    if rep.counters != r.counters {
        diffs.push("counters");
    }
    if rep.factor_base != r.factor_base {
        diffs.push("factor_base");
    }
    if rep.field_ops != r.field_ops {
        diffs.push("field_ops");
    }
    if rep.unpriced != r.cost.unpriced
        || rep.deterministic != r.cost.deterministic
        || Some(rep.automorphisms_used) != r.boundaries.automorphisms_used
    {
        diffs.push("pricing_or_determinism");
    }
    diffs
}

fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let [binary, study] = args.as_slice() else {
        return Err("usage: ecbench_open_questions_replay_failures ECBENCH STUDY_DIR".into());
    };
    let root = Path::new(study);
    let output = root.join("failed-replays-v2.json");
    if output.exists() {
        return Err("refusing to overwrite failed replay receipt".into());
    }
    let binary = fs::canonicalize(binary).map_err(|e| e.to_string())?;
    let digest = sha256_hex(&fs::read(&binary).map_err(|e| e.to_string())?);
    let evidence = root.join("failed-replay-exec-v2");
    if evidence.exists() {
        return Err("refusing to reuse replay evidence directory".into());
    }
    fs::create_dir(&evidence).map_err(|e| e.to_string())?;
    let mut dirs: Vec<_> = fs::read_dir(root.join("sessions"))
        .map_err(|e| e.to_string())?
        .map(|e| e.map(|e| e.path()))
        .collect::<Result<_, _>>()
        .map_err(|e| e.to_string())?;
    dirs.sort();
    let mut results = Vec::new();
    for dir in dirs {
        if !dir.is_dir() {
            continue;
        }
        let audit: Value = serde_json::from_slice(
            &fs::read(root.join(format!(
                "{}.audit.json",
                dir.file_name().unwrap().to_str().unwrap()
            )))
            .map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        if audit["ok"] != true {
            return Err("standard session audit must pass before supplement".into());
        }
        let session: Value =
            serde_json::from_slice(&fs::read(dir.join("session.json")).map_err(|e| e.to_string())?)
                .map_err(|e| e.to_string())?;
        if session["binary_sha256"] != digest {
            return Err("replay binary differs from measured binary".into());
        }
        let s: Spec =
            serde_json::from_slice(&fs::read(dir.join("spec.json")).map_err(|e| e.to_string())?)
                .map_err(|e| e.to_string())?;
        let p = plan(s)?;
        for r in read_records(&dir)?
            .into_iter()
            .filter(|r| !r.warmup && r.outcome.status != "verified")
        {
            let ex = &p.executions[r.seq as usize];
            let w = &p.workloads[ex.workload];
            let input = ChildInput {
                schema: CHILD_INPUT_SCHEMA.into(),
                curve: CurveSpec::explicit(p.instance(ex.workload))
                    .unwrap_or_else(|| w.curve_spec.clone()),
                target_seed: w.target_seed,
                target_index: w.target_index,
                target_kind: w.kind(),
                expected_workload_id: w.workload_id.clone(),
                method: MethodSpec {
                    id: r.method.id.clone(),
                    params: r.method.params.clone(),
                },
                expected_method_id: r.method.method_id.clone(),
                algorithm_seed: ex.algorithm_seed,
            };
            let input_bytes = serde_json::to_vec(&input).map_err(|e| e.to_string())?;
            let stem = format!("{}-{}", r.session_id, r.seq);
            fs::write(evidence.join(format!("{stem}.input.json")), &input_bytes)
                .map_err(|e| e.to_string())?;
            let mut child = Command::new(&binary)
                .arg("exec")
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| e.to_string())?;
            child
                .stdin
                .take()
                .ok_or("missing stdin")?
                .write_all(&input_bytes)
                .map_err(|e| e.to_string())?;
            let out = child.wait_with_output().map_err(|e| e.to_string())?;
            fs::write(evidence.join(format!("{stem}.output.json")), &out.stdout)
                .map_err(|e| e.to_string())?;
            if !out.stderr.is_empty() {
                fs::write(evidence.join(format!("{stem}.stderr")), &out.stderr)
                    .map_err(|e| e.to_string())?;
            }
            let parsed = serde_json::from_slice::<ChildOutput>(&out.stdout);
            let mut diffs = Vec::<String>::new();
            if !out.status.success() {
                diffs.push(format!("child exit {:?}", out.status.code()));
            }
            match parsed {
                Err(e) => diffs.push(format!("child output: {e}")),
                Ok(out) => {
                    if (out.report.is_some()
                        && (out.workload_id.as_ref() != Some(&r.workload.workload_id)
                            || out.method_id.as_ref() != Some(&r.method.method_id)))
                        || out
                            .workload_id
                            .as_ref()
                            .is_some_and(|id| id != &r.workload.workload_id)
                        || out
                            .method_id
                            .as_ref()
                            .is_some_and(|id| id != &r.method.method_id)
                    {
                        diffs.push("identities".into());
                    }
                    match out.report {
                        Some(rep) => {
                            diffs.extend(differences(&r, &rep).into_iter().map(str::to_string))
                        }
                        None => {
                            if r.outcome.status != "error" || out.error != r.outcome.error {
                                diffs.push("error_outcome".into());
                            }
                        }
                    }
                }
            }
            results.push(json!({"record_id":r.record_id,"session_id":r.session_id,"seq":r.seq,"status":r.outcome.status,
                "reproduced":diffs.is_empty(),"differences":diffs,"input_sha256":sha256_hex(&input_bytes),
                "output_sha256":sha256_hex(&out.stdout),"evidence_stem":format!("failed-replay-exec-v2/{stem}")}));
        }
    }
    let ok = results.iter().all(|r| r["reproduced"] == true);
    fs::write(
        output,
        format!(
            "{}\n",
            serde_json::to_string_pretty(&json!({"schema":"ecbench-failed-replay-supplement/v1",
        "binary_sha256":digest,"ok":ok,"same_host":true,"replays":results,
        "scope":"failure status and stable counts only; no independent-host or wall-time claim"}))
            .unwrap()
        ),
    )
    .map_err(|e| e.to_string())?;
    if ok {
        Ok(())
    } else {
        Err("failed replay differences; receipt and outputs preserved".into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn wall_time_does_not_mask_a_changed_integer_count() {
        let a = PhaseRecord {
            name: "setup".into(),
            adds: 3,
            doubles: 2,
            scalar_mults: 1,
            gae: 5.,
            wall_ns: Some(100),
            native: Default::default(),
        };
        let mut b = a.clone();
        b.wall_ns = Some(900);
        assert_eq!(
            phase_projection(&[a.clone()]),
            phase_projection(&[b.clone()])
        );
        b.adds += 1;
        assert_ne!(phase_projection(&[a]), phase_projection(&[b]));
    }
}
