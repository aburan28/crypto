//! Recheck downloaded campaigns with the strict native reader and retained executables.
use isogeny_algos::{json::Json, sha256::sha256_hex};
use std::{collections::BTreeMap, io::Write, path::Path, process::Command};

fn s<'a>(j: &'a Json, k: &str) -> Result<&'a str, String> {
    j.get(k)
        .and_then(Json::as_str)
        .ok_or_else(|| format!("missing {k}"))
}
fn read(path: &Path) -> Result<Json, String> {
    Json::parse(&std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?)
}
fn run_file(
    path: &Path,
    binary: &str,
    revision: &str,
    ells: &str,
) -> Result<(Json, Vec<Json>), String> {
    let bytes = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut rows = bytes
        .lines()
        .map(Json::parse)
        .collect::<Result<Vec<_>, _>>()?;
    if rows.len() < 3 {
        return Err("incomplete CPU run".into());
    }
    let header = rows.remove(0);
    if s(&header, "binary_sha256")? != binary
        || s(&header, "source_revision_supplied_by_launcher")? != revision
    {
        return Err(format!("wrong executable or source in {}", path.display()));
    }
    let Some(Json::Arr(args)) = header.get("args") else {
        return Err("missing arguments".into());
    };
    let value = |flag| {
        args.windows(2)
            .find(|p| p[0].as_str() == Some(flag))
            .and_then(|p| p[1].as_str())
    };
    if value("--bits") != Some("192,224,256,511") || value("--ells") != Some(ells) {
        return Err("unexpected CPU workload".into());
    }
    let exit = rows.pop().ok_or("missing exit")?;
    if s(&exit, "schema")? != "isogeny-perf.exit/v1" || s(&exit, "status")? != "PASS" {
        return Err("incomplete run".into());
    }
    Ok((header, rows))
}
fn median(row: &Json) -> Result<f64, String> {
    let Some(Json::Arr(values)) = row.get("diagnostic_ns") else {
        return Err("missing native samples".into());
    };
    let mut ns = values
        .iter()
        .map(|v| {
            v.as_str()
                .ok_or("non-string sample")?
                .parse::<u64>()
                .map_err(|_| "bad sample")
        })
        .collect::<Result<Vec<_>, _>>()?;
    if ns.len() != 5 || ns.contains(&0) {
        return Err("expected five positive native samples".into());
    }
    ns.sort();
    Ok(ns[2] as f64)
}
fn audit(args: &[String], out: &Path) -> Result<Json, String> {
    let reader = &args[1];
    let root = Path::new(&args[2]);
    let hashes = ["baseline", "candidate"]
        .map(|name| std::fs::read(root.join("binaries").join(name)).map(|b| sha256_hex(&b)));
    let hashes = [
        hashes[0].as_ref().map_err(|e| e.to_string())?,
        hashes[1].as_ref().map_err(|e| e.to_string())?,
    ];
    let compare = |a: &Path, b: &Path, label: &str, equal_counts: bool| -> Result<(), String> {
        let result = out.join(format!("{label}.json"));
        if !Command::new(reader)
            .arg("compare")
            .arg(a)
            .arg(b)
            .arg(&result)
            .status()
            .map_err(|e| e.to_string())?
            .success()
        {
            return Err(format!("strict native comparison rejected {label}"));
        }
        let j = read(&result)?;
        if s(&j, "status")? != "PASS" {
            return Err("comparison did not pass".into());
        }
        if equal_counts {
            let Some(Json::Arr(cases)) = j.get("cases") else {
                return Err("missing cases".into());
            };
            for c in cases {
                for op in ["add_sub", "mul", "inv", "neg"] {
                    let v = c
                        .get("operations")
                        .and_then(|v| v.get(op))
                        .ok_or("missing counts")?;
                    if s(v, "baseline")? != s(v, "candidate")? {
                        return Err(format!("nonrepeatable counts in {label}"));
                    }
                }
            }
        }
        Ok(())
    };
    compare(
        &root.join("baseline.jsonl"),
        &root.join("candidate.jsonl"),
        "full-matrix",
        false,
    )?;
    run_file(
        &root.join("baseline.jsonl"),
        hashes[0],
        &args[3],
        "7,31,101",
    )?;
    run_file(
        &root.join("candidate.jsonl"),
        hashes[1],
        &args[4],
        "7,31,101",
    )?;
    let campaign = root.join("campaign");
    let summary_path = campaign.join("SUMMARY.json");
    let summary = read(&summary_path)?;
    if s(&summary, "status")? != "PASS" {
        return Err("producer campaign is censored or mismatched; full-matrix comparison retained separately".into());
    }
    let Some(Json::Arr(attempts)) = summary.get("attempts") else {
        return Err("missing attempts".into());
    };
    if attempts.len() != 20 {
        return Err("expected twenty attempts".into());
    }
    let mut data = BTreeMap::new();
    let mut host = None;
    let mut qualified = 0;
    for phase in ["aa", "ab"] {
        for round in 0..5 {
            for arm in ["a", "b"] {
                let name = format!("{phase}-{round}-{arm}");
                let matching: Vec<_> = attempts
                    .iter()
                    .filter(|a| a.get("name").and_then(Json::as_str) == Some(&name))
                    .collect();
                if matching.len() != 1 {
                    return Err("missing or repeated attempt".into());
                }
                let a = matching[0];
                let base = phase == "aa" || arm == "a";
                if a.get("valid") != Some(&Json::Bool(true))
                    || a.get("child_exit").and_then(Json::as_f64) != Some(0.0)
                    || a.get("baseline_binary") != Some(&Json::Bool(base))
                {
                    return Err(format!("invalid arm {name}"));
                }
                let path = campaign.join(format!("{name}.jsonl"));
                if s(a, "output")? != format!("{name}.jsonl") {
                    return Err("foreign output path".into());
                }
                let (header, rows) = run_file(
                    &path,
                    hashes[usize::from(!base)],
                    &args[if base { 3 } else { 4 }],
                    "7,31",
                )?;
                let current_host = header.get("host").ok_or("missing host")?;
                if host.as_ref().is_some_and(|h| h != current_host) {
                    return Err("mixed campaign hosts".into());
                }
                host = Some(current_host.clone());
                let reference = campaign.join(if base { "aa-0-a.jsonl" } else { "ab-0-b.jsonl" });
                compare(&reference, &path, &name, true)?;
                let mut q = false;
                if let Some(receipt) = a.get("isolation_receipt").and_then(Json::as_str) {
                    if receipt != format!("{name}-isolation.jsonl") {
                        return Err("foreign isolation path".into());
                    }
                    let raw = std::fs::read_to_string(campaign.join(receipt))
                        .map_err(|e| e.to_string())?;
                    let lines: Vec<_> = raw.lines().filter(|l| !l.trim().is_empty()).collect();
                    if lines.len() == 1 {
                        let iso = Json::parse(lines[0])?;
                        if let Some(run) = iso.get("run") {
                            let command_matches = match (run.get("command"), header.get("args")) {
                                (Some(Json::Arr(cmd)), Some(Json::Arr(argv))) => {
                                    cmd.get(1..) == Some(argv.as_slice())
                                }
                                _ => false,
                            };
                            q = s(&iso, "schema")? == "isolated-bench/1"
                                && s(&iso, "mode")? == "run"
                                && s(&iso, "label")? == name
                                && iso
                                    .get("tool")
                                    .and_then(|t| t.get("implementation"))
                                    .and_then(Json::as_str)
                                    == Some("native")
                                && run.get("contended") == Some(&Json::Bool(false))
                                && run.get("exit_status").and_then(Json::as_f64) == Some(0.0)
                                && command_matches;
                        }
                    }
                }
                if a.get("qualified_timing") != Some(&Json::Bool(q)) {
                    return Err(format!("unbound timing receipt {name}"));
                }
                qualified += usize::from(q);
                data.insert(name, rows);
            }
        }
    }
    compare(
        &campaign.join("aa-0-a.jsonl"),
        &campaign.join("ab-0-b.jsonl"),
        "cross-arm-output",
        false,
    )?;
    let eligible = qualified == 20;
    if summary.get("qualified_timing_runs").and_then(Json::as_f64) != Some(qualified as f64)
        || summary.get("wall_speedup_eligible") != Some(&Json::Bool(eligible))
    {
        return Err("incorrect timing admission summary".into());
    }
    if eligible {
        let Some(Json::Arr(pairs)) = summary.get("qualified_pairs") else {
            return Err("missing pairs".into());
        };
        if pairs.len() != 640 {
            return Err("expected 640 paired case medians".into());
        }
        for phase in ["aa", "ab"] {
            for round in 0..5 {
                let ar = &data[&format!("{phase}-{round}-a")];
                let br = &data[&format!("{phase}-{round}-b")];
                for a in ar {
                    let same_cell = |b: &&Json| {
                        ["bits", "ell", "algorithm"]
                            .iter()
                            .all(|k| a.get(k) == b.get(k))
                    };
                    let b = br.iter().find(same_cell).ok_or("missing paired cell")?;
                    let matching: Vec<_> = pairs
                        .iter()
                        .filter(|p| {
                            p.get("phase").and_then(Json::as_str) == Some(phase)
                                && p.get("round").and_then(Json::as_f64) == Some(round as f64)
                                && ["bits", "ell", "algorithm"]
                                    .iter()
                                    .all(|k| a.get(k) == p.get(k))
                        })
                        .collect();
                    if matching.len() != 1 {
                        return Err("missing or duplicate median pair".into());
                    }
                    let p = matching[0];
                    let (am, bm) = (median(a)?, median(b)?);
                    if p.get("a_median_ns").and_then(Json::as_f64) != Some(am)
                        || p.get("b_median_ns").and_then(Json::as_f64) != Some(bm)
                        || p.get("a_over_b").and_then(Json::as_f64).is_none_or(|r| {
                            !r.is_finite() || (r - am / bm).abs() > 1e-12 * (am / bm)
                        })
                    {
                        return Err("summary medians differ from raw native samples".into());
                    }
                }
            }
        }
    }
    Ok(Json::obj(vec![
        ("schema", Json::str("isogeny-perf.artifact-audit/v1")),
        ("status", Json::str("PASS")),
        ("valid_runs", Json::Num(20)),
        ("qualified_timing_runs", Json::Num(qualified as i64)),
        ("wall_speedup_eligible", Json::Bool(eligible)),
        ("actual_host", host.unwrap()),
        ("baseline_binary_sha256", Json::str(hashes[0])),
        ("candidate_binary_sha256", Json::str(hashes[1])),
        (
            "summary_sha256",
            Json::str(sha256_hex(
                &std::fs::read(summary_path).map_err(|e| e.to_string())?,
            )),
        ),
    ]))
}
fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 6 {
        return Err("audit-cpu NATIVE_READER ARTIFACT_DIR BASE_SHA CANDIDATE_SHA FRESH_DIR".into());
    }
    let out = Path::new(&args[5]);
    std::fs::create_dir(out).map_err(|e| e.to_string())?;
    let result = audit(&args, out);
    let value = match &result {
        Ok(j) => j.clone(),
        Err(e) => Json::obj(vec![
            ("schema", Json::str("isogeny-perf.artifact-audit/v1")),
            ("status", Json::str("CENSORED")),
            ("reason", Json::str(e)),
        ]),
    };
    writeln!(
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(out.join("AUDIT.json"))
            .map_err(|e| e.to_string())?,
        "{}",
        value.dump()
    )
    .map_err(|e| e.to_string())?;
    result.map(|_| ())
}
