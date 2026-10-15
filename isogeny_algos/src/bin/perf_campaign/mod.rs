//! Five A/A pairs and five interleaved A/B pairs, with optional native Linux isolation.
use super::{load_cases, text_field, Json};
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

pub fn run(args: &[String]) -> Result<(), String> {
    if args.len() != 5 && args.len() != 7 {
        return Err("campaign BASE_BINARY CANDIDATE_BINARY BASE_SHA CANDIDATE_SHA FRESH_DIR [ISOLATED_BENCH CPU_LIST]".into());
    }
    let base = std::fs::canonicalize(&args[0]).map_err(|e| e.to_string())?;
    let candidate = std::fs::canonicalize(&args[1]).map_err(|e| e.to_string())?;
    let isolator = if args.len() == 7 {
        Some(std::fs::canonicalize(&args[5]).map_err(|e| e.to_string())?)
    } else {
        None
    };
    let dir = Path::new(&args[4]);
    std::fs::create_dir(dir).map_err(|e| e.to_string())?;
    let mut attempts = vec![];
    // The portable matrix intentionally starts with two degrees. Degree 101/257
    // counts are a separate series; they cannot silently replace these fixtures.
    for phase in ["aa", "ab"] {
        for round in 0..5 {
            let order = if phase == "ab" && round % 2 == 1 {
                ["b", "a"]
            } else {
                ["a", "b"]
            };
            for arm in order {
                let is_base = phase == "aa" || arm == "a";
                let binary = if is_base { &base } else { &candidate };
                let revision = &args[if is_base { 2 } else { 3 }];
                let name = format!("{phase}-{round}-{arm}");
                let output = dir.join(format!("{name}.jsonl"));
                let isolation = dir.join(format!("{name}-isolation.jsonl"));
                let mut cmd = Command::new(isolator.as_ref().unwrap_or(binary));
                if isolator.is_some() {
                    cmd.args(["run", "--cpus", &args[6], "--out"])
                        .arg(&isolation)
                        .args(["--label", &name, "--"])
                        .arg(binary);
                }
                cmd.arg("--out").arg(&output).args([
                    "--revision",
                    revision,
                    "--bits",
                    "192,224,256,511",
                    "--ells",
                    "7,31",
                ]);
                if isolator.is_some() {
                    cmd.args(["--diagnostic-timing", "--samples", "5"]);
                }
                cmd.env("RAYON_NUM_THREADS", "1");
                cmd.stdout(Stdio::from(
                    std::fs::File::create(dir.join(format!("{name}.stdout")))
                        .map_err(|e| e.to_string())?,
                ));
                cmd.stderr(Stdio::from(
                    std::fs::File::create(dir.join(format!("{name}.stderr")))
                        .map_err(|e| e.to_string())?,
                ));
                let started = std::time::Instant::now();
                let status = cmd.status();
                let exit = status.as_ref().ok().and_then(|s| s.code());
                let qualified = isolator.is_some() && exit == Some(0) && isolation_pass(&isolation);
                let valid = exit == Some(0)
                    && load_cases(output.to_str().ok_or("non-UTF8 output path")?).is_ok();
                attempts.push(Json::obj(vec![
                    ("name", Json::str(&name)),
                    ("phase", Json::str(phase)),
                    ("round", Json::Num(round)),
                    ("arm", Json::str(arm)),
                    ("baseline_binary", Json::Bool(is_base)),
                    (
                        "output",
                        Json::str(output.file_name().unwrap().to_string_lossy()),
                    ),
                    (
                        "isolation_receipt",
                        if isolator.is_some() {
                            Json::str(isolation.file_name().unwrap().to_string_lossy())
                        } else {
                            Json::Null
                        },
                    ),
                    (
                        "child_exit",
                        exit.map(|e| Json::Num(e as i64)).unwrap_or(Json::Null),
                    ),
                    (
                        "failure",
                        status
                            .err()
                            .map(|e| Json::str(e.to_string()))
                            .unwrap_or(Json::Null),
                    ),
                    (
                        "elapsed_all_phases_ms",
                        Json::str(started.elapsed().as_millis()),
                    ),
                    ("valid", Json::Bool(valid)),
                    ("qualified_timing", Json::Bool(qualified)),
                ]));
                write_json(
                    &dir.join(format!("{name}-attempt.json")),
                    attempts.last().unwrap(),
                )?;
            }
        }
    }
    let mut pairs = vec![];
    let mut reference = None;
    let mut exact = true;
    let mut count_references: [Option<Vec<Json>>; 2] = [None, None];
    let mut repeated_counts = true;
    for attempt in &attempts {
        if attempt.get("valid") != Some(&Json::Bool(true)) {
            exact = false;
            continue;
        }
        let (_, rows) = load_cases(dir.join(text_field(attempt, "output")?).to_str().unwrap())?;
        let arm = if attempt.get("baseline_binary") == Some(&Json::Bool(true)) {
            0
        } else {
            1
        };
        if let Some(previous) = &count_references[arm] {
            if rows.len() != previous.len()
                || rows
                    .iter()
                    .zip(previous)
                    .any(|(x, y)| x.get("operations") != y.get("operations"))
            {
                repeated_counts = false;
            }
        } else {
            count_references[arm] = Some(rows.clone());
        }
        if let Some(ref original) = reference {
            let original: &Vec<Json> = original;
            if rows.len() != original.len()
                || rows.iter().zip(original).any(|(x, y)| {
                    ["bits", "ell", "algorithm", "input_sha256", "output_sha256"]
                        .iter()
                        .any(|k| x.get(k) != y.get(k))
                })
            {
                exact = false;
            }
        } else {
            reference = Some(rows);
        }
    }
    for phase in ["aa", "ab"] {
        for round in 0..5 {
            let a = attempts
                .iter()
                .find(|a| text_field(a, "name").ok() == Some(format!("{phase}-{round}-a").as_str()))
                .unwrap();
            let b = attempts
                .iter()
                .find(|a| text_field(a, "name").ok() == Some(format!("{phase}-{round}-b").as_str()))
                .unwrap();
            if [a, b]
                .iter()
                .any(|a| a.get("qualified_timing") != Some(&Json::Bool(true)))
            {
                continue;
            }
            let (_, ar) = load_cases(dir.join(text_field(a, "output")?).to_str().unwrap())?;
            let (_, br) = load_cases(dir.join(text_field(b, "output")?).to_str().unwrap())?;
            for (x, y) in ar.iter().zip(&br) {
                let am = median_samples(x)?;
                let bm = median_samples(y)?;
                pairs.push(Json::obj(vec![
                    ("phase", Json::str(phase)),
                    ("round", Json::Num(round)),
                    ("bits", x.get("bits").unwrap().clone()),
                    ("ell", x.get("ell").unwrap().clone()),
                    ("algorithm", x.get("algorithm").unwrap().clone()),
                    ("a_median_ns", Json::Float(am)),
                    ("b_median_ns", Json::Float(bm)),
                    ("a_over_b", Json::Float(am / bm)),
                ]));
            }
        }
    }
    let valid = attempts
        .iter()
        .filter(|a| a.get("valid") == Some(&Json::Bool(true)))
        .count();
    let qualified = attempts
        .iter()
        .filter(|a| a.get("qualified_timing") == Some(&Json::Bool(true)))
        .count();
    let result=Json::obj(vec![("schema",Json::str("isogeny-perf.campaign/v1")),
        ("status",Json::str(if exact && repeated_counts && valid==20 {"PASS"}else{"CENSORED_OR_MISMATCH"})),
        ("exact_outputs_all_runs",Json::Bool(exact && valid==20)),("valid_runs",Json::Num(valid as i64)),
        ("repeatable_counts_each_binary",Json::Bool(repeated_counts && valid==20)),
        ("qualified_timing_runs",Json::Num(qualified as i64)),("requested_runs",Json::Num(20)),
        ("wall_speedup_eligible",Json::Bool(exact && repeated_counts && qualified==20)),
        ("timing_note",Json::str("Report five round ratios and A/A spread per case. No pooled processor or field-width claim; VM host interference remains outside the reservation.")),
        ("attempts",Json::Arr(attempts)),("qualified_pairs",Json::Arr(pairs))]);
    write_json(&dir.join("SUMMARY.json"), &result)?;
    if !exact || !repeated_counts || valid != 20 {
        return Err("campaign has censored or mismatched runs; see SUMMARY.json".into());
    }
    Ok(())
}
fn write_json(path: &Path, j: &Json) -> Result<(), String> {
    let mut f = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| e.to_string())?;
    writeln!(f, "{}", j.dump()).map_err(|e| e.to_string())
}
fn isolation_pass(path: &Path) -> bool {
    let Ok(raw) = std::fs::read_to_string(path) else {
        return false;
    };
    let rows: Vec<_> = raw.lines().filter(|l| !l.trim().is_empty()).collect();
    if rows.len() != 1 {
        return false;
    }
    let Ok(row) = Json::parse(rows[0]) else {
        return false;
    };
    if row.get("schema").and_then(Json::as_str) != Some("isolated-bench/1")
        || row.get("mode").and_then(Json::as_str) != Some("run")
        || row
            .get("tool")
            .and_then(|t| t.get("implementation"))
            .and_then(Json::as_str)
            != Some("native")
    {
        return false;
    }
    let Some(result) = row.get("run") else {
        return false;
    };
    result.get("contended") == Some(&Json::Bool(false))
        && result.get("exit_status").and_then(Json::as_f64) == Some(0.0)
}
fn median_samples(j: &Json) -> Result<f64, String> {
    let Some(Json::Arr(values)) = j.get("diagnostic_ns") else {
        return Err("missing samples".into());
    };
    if values.len() != 5 {
        return Err("expected five samples".into());
    }
    let mut samples: Vec<f64> = values
        .iter()
        .map(|v| {
            v.as_str()
                .ok_or("sample not integer")?
                .parse::<u64>()
                .map(|n| n as f64)
                .map_err(|_| "invalid sample")
        })
        .collect::<Result<_, _>>()?;
    if samples.iter().any(|s| *s <= 0.0) {
        return Err("zero duration".into());
    }
    samples.sort_by(f64::total_cmp);
    Ok(samples[2])
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn odd_median_and_reject_unrecorded_time() {
        let row = Json::obj(vec![(
            "diagnostic_ns",
            Json::Arr([5, 1, 3, 2, 4].into_iter().map(Json::str).collect()),
        )]);
        assert_eq!(median_samples(&row).unwrap(), 3.0);
        assert!(median_samples(&Json::Null).is_err());
    }
}
