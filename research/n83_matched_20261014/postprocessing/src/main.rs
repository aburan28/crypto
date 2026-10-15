use n83_matched_supplement::diagnostic_source_equal;
use serde_json::{json, Value};
use std::{
    env,
    fs::{self, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
    process::Command,
};

fn read(path: &Path) -> Result<Value, Box<dyn std::error::Error>> {
    Ok(serde_json::from_slice(&fs::read(path)?)?)
}
fn read_lines(path: &Path) -> Result<Vec<Value>, Box<dyn std::error::Error>> {
    fs::read_to_string(path)?
        .lines()
        .filter(|s| !s.trim().is_empty())
        .map(|line| Ok(serde_json::from_str(line)?))
        .collect()
}
fn hash(path: &Path) -> Result<String, Box<dyn std::error::Error>> {
    let result = Command::new("sha256sum").arg("--").arg(path).output()?;
    if !result.status.success() {
        return Err(format!("sha256sum failed for {}", path.display()).into());
    }
    let text = String::from_utf8(result.stdout)?;
    let digest = text.split_whitespace().next().ok_or("missing digest")?;
    if digest.len() != 64 || !digest.chars().all(|c| c.is_ascii_hexdigit()) {
        return Err("invalid digest output".into());
    }
    Ok(digest.to_owned())
}
fn main() {
    let args = env::args().collect::<Vec<_>>();
    if !(5..=6).contains(&args.len()) {
        eprintln!("usage: n83-matched-supplement INPUT_JSON FROZEN_STATISTICS_JSON CAMPAIGN_ROOT NEW_OUTPUT_JSON [COMMITTED_CAMPAIGN_PREFIX]");
        std::process::exit(2);
    }
    let result = (|| -> Result<(), Box<dyn std::error::Error>> {
        let input_path = Path::new(&args[1]);
        let stats_path = Path::new(&args[2]);
        let root = Path::new(&args[3]);
        let prefix = args
            .get(5)
            .cloned()
            .unwrap_or_else(|| root.to_string_lossy().into_owned());
        let mut input = read(input_path)?;
        let mut stats = read(stats_path)?;
        n83_matched_supplement::validate_gate(&input, &stats)?;
        let manifest_hash = hash(&root.join("campaign_manifest.json"))?;
        if input["metadata"]["campaign_manifest_sha256"] != manifest_hash {
            return Err("actual campaign manifest hash differs".into());
        }
        let mut source_bindings = Vec::new();
        for row in input["process_rows"]
            .as_array_mut()
            .ok_or("missing processes")?
        {
            let id = row["run_id"].as_str().ok_or("missing run id")?.to_owned();
            if id.is_empty() || !id.chars().all(|c| c.is_ascii_alphanumeric() || c == '-') {
                return Err("unsafe run id".into());
            }
            let file = root.join("runs").join(&id).join("observation.json");
            let observed = read(&file)?;
            if !diagnostic_source_equal(&observed["summary"], &row["summary"])
                || !diagnostic_source_equal(
                    &observed["raw_process_receipt"],
                    &row["raw_process_receipt"],
                )
                || !diagnostic_source_equal(
                    &json!({"producer_preparation_ms":observed["producer_preparation_ms"]}),
                    &json!({"producer_preparation_ms":row["producer_preparation_ms"]}),
                )
            {
                return Err(format!("process source binding differs {id}").into());
            }
            for key in [
                "run_id",
                "kind",
                "arm",
                "target_count",
                "cohort_id",
                "target_id",
                "elapsed_process_wall_s",
                "cpu_user_s",
                "cpu_system_s",
                "cpu_total_s",
                "peak_rss_bytes",
            ] {
                if !diagnostic_source_equal(&json!({key:observed[key]}), &json!({key:row[key]})) {
                    return Err(format!("process source field {key} differs for {id}").into());
                }
            }
            let execution = read(&root.join("runs").join(&id).join("execution.json"))?;
            if !diagnostic_source_equal(&execution, &observed["raw_process_receipt"]) {
                return Err(format!("execution source differs {id}").into());
            }
            let stdout = read_lines(&root.join("runs").join(&id).join("stdout.jsonl"))?;
            let summaries = stdout
                .iter()
                .filter(|v| {
                    matches!(
                        v["kind"].as_str(),
                        Some("compact_orbit_dlp_summary" | "n83_rho_matched_summary")
                    )
                })
                .collect::<Vec<_>>();
            if summaries.len() != 1 || !diagnostic_source_equal(summaries[0], &observed["summary"])
            {
                return Err(format!("producer summary differs {id}").into());
            }
            row["summary"] = summaries[0].clone();
            row["raw_process_receipt"] = execution.clone();
            for key in [
                "elapsed_process_wall_s",
                "cpu_user_s",
                "cpu_system_s",
                "cpu_total_s",
                "peak_rss_bytes",
            ] {
                row[key] = execution[key].clone();
            }
            let timing = &row["summary"]["timing_ms"];
            let prep = if row["arm"] == "ic" {
                timing["setup_total"]
                    .as_f64()
                    .zip(timing["rank_stage"].as_f64())
                    .zip(timing["linear_algebra"].as_f64())
                    .map(|((a, b), c)| a + b + c)
            } else {
                timing["setup_total"].as_f64()
            };
            row["producer_preparation_ms"] = json!(prep);
            let source = json!({"file":format!("{prefix}/runs/{id}/observation.json"),"sha256":hash(&file)?,"json_pointer":"/","summary_json_pointer":"/summary/timing_ms","execution_file":format!("{prefix}/runs/{id}/execution.json"),"execution_sha256":hash(&root.join("runs").join(&id).join("execution.json"))?});
            source_bindings.push(source.clone());
            row["source"] = source;
        }
        for field in ["target_rows", "auxiliary_target_rows"] {
            for row in input[field].as_array_mut().ok_or("missing target rows")? {
                let id = row["run_id"].as_str().ok_or("missing run id")?.to_owned();
                let file = root.join("runs").join(&id).join("observation.json");
                let observed = read(&file)?;
                let rows = observed["target_rows"]
                    .as_array()
                    .ok_or("source target rows missing")?;
                let matching = rows
                    .iter()
                    .enumerate()
                    .filter(|(_, r)| r["target_id"] == row["target_id"])
                    .collect::<Vec<_>>();
                if matching.len() != 1
                    || !diagnostic_source_equal(
                        &matching[0].1["rawtargetrecord"],
                        &row["rawtargetrecord"],
                    )
                {
                    return Err(
                        format!("target source binding differs {id}/{}", row["target_id"]).into(),
                    );
                }
                let index = matching[0].0;
                if row["rawtargetrecord"]["fixture_index"].as_u64() != Some(index as u64) {
                    return Err(format!("producer line index binding differs {id}").into());
                }
                let producer = read_lines(&root.join("runs").join(&id).join("targets.jsonl"))?;
                let source = producer.get(index).ok_or("producer target line missing")?;
                if !diagnostic_source_equal(source, &matching[0].1["rawtargetrecord"]) {
                    return Err(format!("producer target differs {id}/{index}").into());
                }
                row["rawtargetrecord"] = source.clone();
                row["online_wall_s"] = json!(source["target_ms"].as_f64().map(|n| n / 1000.));
                row["source"] = json!({"file":format!("{prefix}/runs/{id}/observation.json"),"sha256":hash(&file)?,"json_pointer":format!("/target_rows/{index}/rawtargetrecord"),"producer_file":format!("{prefix}/runs/{id}/targets.jsonl"),"producer_sha256":hash(&root.join("runs").join(&id).join("targets.jsonl"))?,"producer_line_1_based":index+1});
            }
        }
        // Validate original frozen rows first; only source metadata is appended.
        stats["process_rows"] = input["process_rows"].clone();
        stats["target_rows"] = input["target_rows"].clone();
        stats["auxiliary_target_rows"] = input["auxiliary_target_rows"].clone();
        let mut output = n83_matched_supplement::supplemental(&input, &stats)
            .map_err(|s| format!("supplement rejected: {s}"))?;
        output["source_bindings"] = json!(source_bindings);
        output["provenance"] = json!({"normalized_input":{"path":input_path,"sha256":hash(input_path)?},"frozen_statistics":{"path":stats_path,"sha256":hash(stats_path)?,"modified":false,"recomputed":false},"campaign_manifest_sha256":manifest_hash,"source":{"src/lib.rs":hash(&PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/lib.rs"))?,"src/main.rs":hash(&PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/main.rs"))?,"Cargo.lock":hash(&PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("Cargo.lock"))?},"executable_sha256":hash(&env::current_exe()?)?,"native_calculations_only":true,"calculation_authority":"Original producer targets.jsonl/stdout.jsonl and execution.json; normalization drift at named diagnostic f64 clocks permits at most1ULP. No integer/type/provenance relaxation.","source_pointer_prefix":prefix});
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&args[4])?;
        serde_json::to_writer_pretty(&mut file, &output)?;
        file.write_all(b"\n")?;
        println!(
            "{}",
            json!({"status":"PASS_SUPPLEMENTAL_PRESENTATION_DATA","output":args[4],"counts":output["counts"],"frozen_primary_recomputed":false})
        );
        Ok(())
    })();
    if let Err(e) = result {
        eprintln!("{e}");
        std::process::exit(1);
    }
}
