//! Native validation and statistical comparison of the frozen hash study.
//! `cargo run --release --example f4_f5_hash_compare -- RECORDS_DIRECTORY`
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

type Result<T> = std::result::Result<T, String>;

fn number(v: &Value, field: &str) -> Result<u64> {
    v[field]
        .as_u64()
        .ok_or_else(|| format!("missing/invalid {field}"))
}

fn key(v: &Value, repeat: bool) -> String {
    if repeat {
        json!([v["case"], v["seed"], v["repeat"], v["arm"]]).to_string()
    } else {
        json!([v["case"], v["seed"], v["arm"]]).to_string()
    }
}

fn comparable(v: &Value) -> Value {
    let mut v = v.clone();
    for field in [
        "repeat",
        "allocations",
        "allocated_bytes",
        "live_at_entry_bytes",
        "peak_live_bytes",
    ] {
        v.as_object_mut().unwrap().remove(field);
    }
    v
}

fn expected_keys(repeats: u64) -> BTreeSet<String> {
    let mut cases = Vec::new();
    for case in ["n8_d3", "n8_d4", "n10_d4", "n10_d5"] {
        for seed in [0x2545_f491_4f6c_dd1du64, 0x9e37_79b9_7f4a_7c15] {
            cases.push((case, seed));
        }
    }
    for case in ["holdout_n12_d3", "holdout_n12_d4", "holdout_shifted_d4"] {
        cases.push((case, 0xd1b5_4a32_d192_ed03));
    }
    cases
        .into_iter()
        .flat_map(|(case, seed)| {
            (0..repeats).flat_map(move |repeat| {
                ["f4", "f5", "criterion"]
                    .into_iter()
                    .map(move |arm| json!([case, seed, repeat, arm]).to_string())
            })
        })
        .collect()
}

fn parse_dump(text: &str) -> Result<Option<(String, u64)>> {
    let Some(label) = text
        .lines()
        .find_map(|line| line.strip_prefix("desc: Trigger: Client Request: "))
    else {
        return Ok(None);
    };
    let label: Value = serde_json::from_str(label).map_err(|e| e.to_string())?;
    if !label.is_array() || label.as_array().unwrap().len() != 4 {
        return Err("invalid request label".into());
    }
    // The command uses --separate-threads=no. Refuse separately attributed
    // dumps rather than presenting a main-thread slice as complete work.
    if text.lines().any(|line| line.starts_with("thread:")) {
        return Err("separate-thread dump: worker-inclusive merged profile required".into());
    }
    let events: Vec<_> = text
        .lines()
        .find_map(|line| line.strip_prefix("events: "))
        .ok_or("events missing")?
        .split_whitespace()
        .collect();
    let ir = events
        .iter()
        .position(|&event| event == "Ir")
        .ok_or("Ir missing")?;
    // After ZERO_STATS, Callgrind 3.24's incremental `summary` can wrap.
    // The file's `totals` is the sum of its actual self-cost records.
    let counts: Vec<u64> = text
        .lines()
        .find_map(|line| line.strip_prefix("totals: "))
        .ok_or("totals missing")?
        .split_whitespace()
        .map(|n| n.parse::<u64>().map_err(|e| e.to_string()))
        .collect::<Result<_>>()?;
    if counts.len() != events.len() || counts[ir] == 0 {
        return Err("invalid/zero instruction summary".into());
    }
    Ok(Some((label.to_string(), counts[ir])))
}

struct Run {
    rows: BTreeMap<String, Value>,
    ir: BTreeMap<String, u64>,
    sources: Vec<Value>,
}

fn read_run(root: &Path, name: &str, repeats: u64) -> Result<Run> {
    let directory = root.join(name);
    if std::fs::read_to_string(directory.join("exit-status.txt"))
        .map_err(|e| e.to_string())?
        .trim()
        != "0"
    {
        return Err(format!("{name}: run did not complete successfully"));
    }
    let records = std::fs::read(directory.join("records.jsonl")).map_err(|e| e.to_string())?;
    let mut rows = BTreeMap::new();
    for line in std::str::from_utf8(&records)
        .map_err(|e| e.to_string())?
        .lines()
    {
        let row: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
        let peak = number(&row, "peak_live_bytes")?;
        if peak < number(&row, "live_at_entry_bytes")? {
            return Err("peak below entry live bytes".into());
        }
        number(&row, "allocations")?;
        number(&row, "allocated_bytes")?;
        if row["input_digest"].as_str().is_none() {
            return Err("input digest missing".into());
        }
        if row["arm"] == "criterion" {
            if row["criterion_digest"].as_str().is_none() {
                return Err("criterion membership digest missing".into());
            }
            number(&row, "word_ops")?;
            number(&row, "pruned")?;
            for field in ["lower_rows", "pruned_by_part"] {
                let a = row[field]
                    .as_array()
                    .ok_or("criterion accounting missing")?;
                if a.len() != 2 || a.iter().any(|v| v.as_u64().is_none()) {
                    return Err("invalid criterion accounting".into());
                }
            }
        } else {
            if row["output_digest"].as_str().is_none() {
                return Err("output digest missing".into());
            }
            if row["arm"] == "f4" {
                number(&row, "rank")?;
                number(&row, "word_ops")?;
            } else {
                for field in [
                    "degree",
                    "rows_f4",
                    "rows_built",
                    "rows_pruned",
                    "cols",
                    "rank",
                    "zero_reductions",
                    "criterion_rows",
                    "criterion_word_ops",
                    "reduce_word_ops",
                ] {
                    number(&row["report"], field)?;
                }
            }
        }
        if rows.insert(key(&row, true), row).is_some() {
            return Err("duplicate native record".into());
        }
    }
    if rows.keys().cloned().collect::<BTreeSet<_>>() != expected_keys(repeats) {
        return Err(format!(
            "{name}: incomplete or changed frozen workload coverage"
        ));
    }
    let mut ir = BTreeMap::new();
    let mut sources = vec![
        json!({"path": format!("{name}/records.jsonl"), "blake3":blake3::hash(&records).to_hex().to_string()}),
    ];
    let mut paths: Vec<_> = std::fs::read_dir(&directory)
        .map_err(|e| e.to_string())?
        .map(|entry| entry.map(|e| e.path()).map_err(|e| e.to_string()))
        .collect::<Result<_>>()?;
    paths.sort();
    for path in paths {
        if !path
            .file_name()
            .unwrap()
            .to_string_lossy()
            .starts_with("callgrind.out")
        {
            continue;
        }
        let bytes = std::fs::read(&path).map_err(|e| e.to_string())?;
        if let Some((key, count)) =
            parse_dump(std::str::from_utf8(&bytes).map_err(|e| e.to_string())?)?
        {
            let text = std::str::from_utf8(&bytes).map_err(|e| e.to_string())?;
            if self_costs(text, &BTreeMap::new())?.values().sum::<u64>() != count {
                return Err(format!(
                    "self costs do not sum to recorded totals at {}",
                    path.display()
                ));
            }
            if ir.insert(key, count).is_some() {
                return Err("duplicate Callgrind request dump".into());
            }
            sources.push(json!({"path":format!("{name}/{}",path.file_name().unwrap().to_string_lossy()),"blake3":blake3::hash(&bytes).to_hex().to_string()}));
        }
    }
    if ir.keys().ne(rows.keys()) {
        return Err(format!("{name}: missing/extra request dumps"));
    }
    let mut seen = BTreeMap::new();
    let mut outputs = BTreeMap::<String, String>::new();
    let mut inputs = BTreeMap::<String, String>::new();
    for row in rows.values() {
        let input_key = json!([row["case"], row["seed"], row["repeat"]]).to_string();
        let input = row["input_digest"].as_str().unwrap().to_string();
        if inputs
            .insert(input_key, input.clone())
            .is_some_and(|old| old != input)
        {
            return Err("F4/F5/criterion used different inputs".into());
        }
        let case_key = key(row, false);
        if seen
            .insert(case_key.clone(), comparable(row))
            .is_some_and(|old| old != comparable(row))
        {
            return Err(format!("unstable output/accounting at {case_key}"));
        }
        if row["arm"] != "criterion" {
            let input_key = json!([row["case"], row["seed"], row["repeat"]]).to_string();
            let digest = row["output_digest"].as_str().unwrap().to_string();
            if outputs
                .insert(input_key, digest.clone())
                .is_some_and(|old| old != digest)
            {
                return Err("F4/F5 reduced outputs differ".into());
            }
        }
    }
    Ok(Run { rows, ir, sources })
}

fn check_against(reference: &Run, run: &Run) -> Result<()> {
    let expected: BTreeMap<_, _> = reference
        .rows
        .values()
        .map(|row| (key(row, false), comparable(row)))
        .collect();
    for row in run.rows.values() {
        if expected.get(&key(row, false)) != Some(&comparable(row)) {
            return Err(format!(
                "baseline output/accounting mismatch at {}",
                key(row, true)
            ));
        }
    }
    Ok(())
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    let middle = v.len() / 2;
    if v.len().is_multiple_of(2) {
        (v[middle - 1] + v[middle]) / 2.0
    } else {
        v[middle]
    }
}

fn stats(a: &[u64], b: &[u64], ap: &[u64]) -> Value {
    assert_eq!(a.len(), 5);
    assert_eq!(b.len(), 5);
    assert_eq!(ap.len(), 5);
    assert!(a.iter().chain(b).chain(ap).all(|&n| n > 0));
    let d: Vec<_> = a
        .iter()
        .zip(b)
        .map(|(&a, &b)| (a as f64 / b as f64).ln())
        .collect();
    let noise: Vec<_> = a
        .iter()
        .zip(ap)
        .map(|(&a, &ap)| (a as f64 / ap as f64).ln().abs())
        .collect();
    let sigma = 1.4826 * median(noise);
    let observed = median(d.clone());
    let mut state = 0x243f_6a88_85a3_08d3u64;
    let mut boot = Vec::with_capacity(2000);
    for _ in 0..2000 {
        let mut sample = Vec::with_capacity(5);
        for _ in 0..5 {
            // Xorshift64's nonzero cycle has 2^64-1 states, divisible by 5.
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            sample.push(d[(state % 5) as usize]);
        }
        boot.push(median(sample));
    }
    boot.sort_by(f64::total_cmp);
    let interval = [boot[49], boot[1949]];
    let threshold = (2.0 * sigma).max(1.02f64.ln());
    let verdict = if interval[0] > 0.0 && observed > threshold {
        "fewer_instructions"
    } else if interval[1] < 0.0 && -observed > threshold {
        "more_instructions"
    } else {
        "neutral"
    };
    let base = median(a.iter().map(|&n| n as f64).collect());
    let candidate = median(b.iter().map(|&n| n as f64).collect());
    json!({"baseline_ir":a,"candidate_ir":b,"repeat_baseline_ir":ap,
        "baseline_median_ir":base,"candidate_median_ir":candidate,
        "median_log_ratio":observed,"instruction_ratio":observed.exp(),
        "instruction_reduction_percent":100.0*(1.0-(-observed).exp()),
        "log_ratio_ci95":interval,"instruction_reduction_ci95_percent":interval.map(|x|100.0*(1.0-(-x).exp())),
        "aa_sigma":sigma,"log_threshold":threshold,"verdict":verdict})
}

fn total(run: &Run, arm: &str) -> u64 {
    run.rows
        .iter()
        .filter(|(_, v)| v["arm"] == arm)
        .map(|(key, _)| run.ir[key])
        .sum()
}

fn function_names(text: &str) -> Result<BTreeMap<String, String>> {
    let mut names = BTreeMap::new();
    for line in text.lines() {
        if let Some(function) = line
            .strip_prefix("fn=")
            .or_else(|| line.strip_prefix("cfn="))
        {
            if let Some(rest) = function.strip_prefix('(') {
                let (id, name) = rest.split_once(')').ok_or("malformed function id")?;
                if !name.trim().is_empty() {
                    names.insert(id.to_string(), name.trim().to_string());
                }
            }
        }
    }
    Ok(names)
}

fn self_costs(text: &str, known_names: &BTreeMap<String, String>) -> Result<BTreeMap<String, u64>> {
    let mut names = known_names.clone();
    let mut costs = BTreeMap::<String, u64>::new();
    let mut current = String::new();
    let mut called_cost = false;
    for line in text.lines() {
        if let Some(function) = line
            .strip_prefix("fn=")
            .or_else(|| line.strip_prefix("cfn="))
        {
            let (id, name) = if let Some(rest) = function.strip_prefix('(') {
                let (id, rest) = rest.split_once(')').ok_or("malformed function id")?;
                (id.to_string(), rest.trim())
            } else {
                (function.to_string(), function)
            };
            if !name.is_empty() {
                names.insert(id.clone(), name.to_string());
            }
            if line.starts_with("fn=") {
                current = id;
                called_cost = false;
            }
        } else if line.starts_with("calls=") {
            called_cost = true;
        } else if !current.is_empty()
            && line
                .as_bytes()
                .first()
                .is_some_and(|b| b.is_ascii_digit() || b"+-*".contains(b))
        {
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() != 2 {
                return Err("self-cost inspection requires positions: line and events: Ir".into());
            }
            let cost: u64 = fields[1].parse::<u64>().map_err(|e| e.to_string())?;
            if !called_cost {
                *costs.entry(current.clone()).or_default() += cost;
            }
            called_cost = false;
        }
    }
    let mut result = BTreeMap::new();
    for (id, cost) in costs {
        *result
            .entry(
                names
                    .get(&id)
                    .cloned()
                    .unwrap_or_else(|| format!("<function {id}>")),
            )
            .or_default() += cost;
    }
    Ok(result)
}

fn inspect(root: &Path, name: &str, repeats: u64) -> Result<Value> {
    let run = read_run(root, name, repeats)?;
    let mut profile = BTreeMap::<String, u64>::new();
    let mut criterion_summary = 0u64;
    let mut names = BTreeMap::new();
    for source in &run.sources {
        let path = source["path"].as_str().unwrap();
        if !path.ends_with("records.jsonl") {
            let text = std::fs::read_to_string(root.join(path)).map_err(|e| e.to_string())?;
            names.extend(function_names(&text)?);
        }
    }
    for source in &run.sources {
        let path = source["path"].as_str().unwrap();
        if path.ends_with("records.jsonl") {
            continue;
        }
        let text = std::fs::read_to_string(root.join(path)).map_err(|e| e.to_string())?;
        let (label, ir) = parse_dump(&text)?.ok_or("request missing")?;
        let label: Value = serde_json::from_str(&label).map_err(|e| e.to_string())?;
        if label[3] != "criterion" {
            continue;
        }
        let costs = self_costs(&text, &names)?;
        if costs.values().sum::<u64>() != ir {
            return Err(format!(
                "self costs {} do not sum to Ir {ir} at {path}",
                costs.values().sum::<u64>()
            ));
        }
        criterion_summary += ir;
        for (function, cost) in costs {
            *profile.entry(function).or_default() += cost;
        }
    }
    let mut functions: Vec<_> = profile.into_iter().collect();
    functions.sort_by_key(|(_, cost)| std::cmp::Reverse(*cost));
    let hashing_ir: u64 = functions
        .iter()
        .filter(|(name, _)| {
            name.contains("::hash::sip::")
                || (name.contains("RandomState") && name.contains("hash_one"))
        })
        .map(|(_, ir)| ir)
        .sum();
    let arms: Vec<_> = ["f4", "f5", "criterion"]
        .into_iter()
        .map(|arm| json!({"arm":arm,"calls":11*repeats,"ir":total(&run,arm)}))
        .collect();
    Ok(
        json!({"run":name,"repeats":repeats,"systems":11,"merged_worker_profile":true,
        "arms":arms,
        "criterion_self_cost_sum":criterion_summary,"standard_hash_self_ir":hashing_ir,
        "standard_hash_self_percent":100.0*hashing_ir as f64/criterion_summary as f64,
        "criterion_hotspots":functions.into_iter().take(20).map(|(function,ir)|json!({"function":function,"ir":ir})).collect::<Vec<_>>(),"sources":run.sources}),
    )
}

fn analyze(root: &Path) -> Result<Value> {
    let a = read_run(root, "primary-a", 3)?;
    let b = read_run(root, "primary-b", 3)?;
    let aa = read_run(root, "aa-a", 3)?;
    check_against(&a, &b)?;
    check_against(&a, &aa)?;
    let mut rounds = Vec::new();
    for r in 0..5 {
        let a1 = read_run(root, &format!("round-{r}-a"), 1)?;
        let b1 = read_run(root, &format!("round-{r}-b"), 1)?;
        let ap = read_run(root, &format!("round-{r}-aprime"), 1)?;
        for run in [&a1, &b1, &ap] {
            check_against(&a, run)?;
        }
        rounds.push((a1, b1, ap));
    }
    let mut memory_regressions = Vec::new();
    let mut per_case = Vec::new();
    for (label, ar) in &a.rows {
        let br = &b.rows[label];
        let ma = number(ar, "peak_live_bytes")? - number(ar, "live_at_entry_bytes")?;
        let mb = number(br, "peak_live_bytes")? - number(br, "live_at_entry_bytes")?;
        if number(br, "allocations")? > number(ar, "allocations")?
            || number(br, "allocated_bytes")? > number(ar, "allocated_bytes")?
            || number(br, "peak_live_bytes")? > number(ar, "peak_live_bytes")?
            || mb > ma
        {
            memory_regressions.push(label.clone());
        }
        per_case.push(json!({"key":label,"case":ar["case"],"seed":ar["seed"],"repeat":ar["repeat"],"arm":ar["arm"],
            "baseline_ir":a.ir[label],"candidate_ir":b.ir[label],
            "baseline_allocations":ar["allocations"],"candidate_allocations":br["allocations"],
            "baseline_allocated_bytes":ar["allocated_bytes"],"candidate_allocated_bytes":br["allocated_bytes"],
            "baseline_live_at_entry_bytes":ar["live_at_entry_bytes"],"candidate_live_at_entry_bytes":br["live_at_entry_bytes"],
            "baseline_peak_live_bytes":ar["peak_live_bytes"],"candidate_peak_live_bytes":br["peak_live_bytes"],
            "baseline_peak_additional_live_bytes":ma,"candidate_peak_additional_live_bytes":mb}));
    }
    let mut memory_cases = Vec::new();
    let mut keys_seen = BTreeSet::new();
    for row in a.rows.values() {
        let case_key = key(row, false);
        if !keys_seen.insert(case_key.clone()) {
            continue;
        }
        let mut result = json!({"key":case_key,"case":row["case"],"seed":row["seed"],"arm":row["arm"],"samples":3});
        for (role, run) in [("baseline", &a), ("candidate", &b)] {
            for field in [
                "allocations",
                "allocated_bytes",
                "live_at_entry_bytes",
                "peak_live_bytes",
            ] {
                let values: Vec<_> = run
                    .rows
                    .values()
                    .filter(|r| key(r, false) == case_key)
                    .map(|r| number(r, field).unwrap() as f64)
                    .collect();
                result[format!("{role}_{field}")] = json!(median(values));
            }
            let values: Vec<_> = run
                .rows
                .values()
                .filter(|r| key(r, false) == case_key)
                .map(|r| {
                    (number(r, "peak_live_bytes").unwrap()
                        - number(r, "live_at_entry_bytes").unwrap()) as f64
                })
                .collect();
            result[format!("{role}_peak_additional_live_bytes")] = json!(median(values));
        }
        result["additional_live_byte_reduction"] = json!(
            result["baseline_peak_additional_live_bytes"]
                .as_f64()
                .unwrap()
                - result["candidate_peak_additional_live_bytes"]
                    .as_f64()
                    .unwrap()
        );
        memory_cases.push(result);
    }
    for (round, (ar, br, _)) in rounds.iter().enumerate() {
        for (label, a) in &ar.rows {
            let b = &br.rows[label];
            let ma = number(a, "peak_live_bytes")? - number(a, "live_at_entry_bytes")?;
            let mb = number(b, "peak_live_bytes")? - number(b, "live_at_entry_bytes")?;
            if number(b, "allocations")? > number(a, "allocations")?
                || number(b, "allocated_bytes")? > number(a, "allocated_bytes")?
                || number(b, "peak_live_bytes")? > number(a, "peak_live_bytes")?
                || mb > ma
            {
                memory_regressions.push(format!("round-{round}:{label}"));
            }
        }
    }
    let mut arms = Vec::new();
    let mut case_verdicts = Vec::new();
    for arm in ["f4", "f5", "criterion"] {
        let mut av = Vec::new();
        let mut bv = Vec::new();
        let mut apv = Vec::new();
        for (ar, br, ap) in &rounds {
            av.push(total(ar, arm));
            bv.push(total(br, arm));
            apv.push(total(ap, arm));
        }
        let mut result = stats(&av, &bv, &apv);
        result["arm"] = json!(arm);
        result["primary_baseline_ir"] = json!(total(&a, arm));
        result["primary_candidate_ir"] = json!(total(&b, arm));
        result["aa_primary_ir"] = json!(total(&aa, arm));
        result["calls_per_round"] = json!(11);
        result["primary_calls"] = json!(33);
        arms.push(result);
    }
    for (label, row) in rounds[0].0.rows.iter() {
        let mut av = Vec::new();
        let mut bv = Vec::new();
        let mut apv = Vec::new();
        for (ar, br, ap) in &rounds {
            av.push(ar.ir[label]);
            bv.push(br.ir[label]);
            apv.push(ap.ir[label]);
        }
        let mut result = stats(&av, &bv, &apv);
        result["key"] = json!(label);
        result["case"] = row["case"].clone();
        result["seed"] = row["seed"].clone();
        result["arm"] = row["arm"].clone();
        case_verdicts.push(result);
    }
    let f4 = arms.iter().find(|a| a["arm"] == "f4").unwrap();
    let f5 = arms.iter().find(|a| a["arm"] == "f5").unwrap();
    let accepted = f5["verdict"] == "fewer_instructions"
        && f4["verdict"] != "more_instructions"
        && memory_regressions.is_empty();
    let mut sources = Vec::new();
    for run in [&a, &b, &aa] {
        sources.extend_from_slice(&run.sources);
    }
    for (ar, br, ap) in &rounds {
        for run in [ar, br, ap] {
            sources.extend_from_slice(&run.sources);
        }
    }
    Ok(
        json!({"schema":"f4-f5.hash-profile/v1","systems":11,"primary_repeats":3,"paired_rounds":5,
        "bootstrap_resamples":2000,"bootstrap_seed":"0x243f6a8885a308d3","exact_outputs_and_counts":true,
        "measurement":"Callgrind Ir, instrumented complete calls, merged worker coverage",
        "runtime_speedup":null,"rss":null,"acceptance_met":accepted,"memory_regressions":memory_regressions,
        "arms":arms,"case_verdicts":case_verdicts,"primary_cases":per_case,"memory_cases":memory_cases,"sources":sources}),
    )
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let result = if args.get(1).is_some_and(|s| s == "inspect") {
        assert_eq!(args.len(), 5, "inspect RECORDS_DIRECTORY RUN_NAME REPEATS");
        inspect(Path::new(&args[2]), &args[3], args[4].parse().unwrap())
    } else {
        assert_eq!(args.len(), 2, "RECORDS_DIRECTORY");
        analyze(Path::new(&args[1]))
    };
    match result {
        Ok(value) => println!("{}", serde_json::to_string_pretty(&value).unwrap()),
        Err(error) => {
            eprintln!("INVALID comparison: {error}");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn instruction_gate_requires_noise_control_and_reports_regressions() {
        assert_eq!(
            stats(&[1000; 5], &[800; 5], &[1000; 5])["verdict"],
            "fewer_instructions"
        );
        assert_eq!(
            stats(&[1000; 5], &[1200; 5], &[1000; 5])["verdict"],
            "more_instructions"
        );
        assert_eq!(
            stats(&[1000; 5], &[990; 5], &[1000; 5])["verdict"],
            "neutral"
        );
        assert_eq!(
            stats(&[1000; 5], &[900; 5], &[700; 5])["verdict"],
            "neutral"
        );
        assert_eq!(
            stats(
                &[800, 1200, 1000, 800, 1200],
                &[1000; 5],
                &[800, 1200, 1000, 800, 1200]
            )["verdict"],
            "neutral"
        );
    }
    #[test]
    fn profile_parser_rejects_missing_zero_and_partial_thread_counts() {
        let header = "desc: Trigger: Client Request: [\"x\",1,0,\"f5\"]\n";
        assert!(parse_dump(&format!("{header}events: Ir\ntotals: 0\n")).is_err());
        assert!(parse_dump(&format!("{header}events: Ir\n")).is_err());
        assert!(parse_dump(&format!("{header}thread: 1\nevents: Ir\ntotals: 100\n")).is_err());
        assert_eq!(
            parse_dump(&format!("{header}events: Dr Ir\ntotals: 12 100\n"))
                .unwrap()
                .unwrap()
                .1,
            100
        );
    }
    #[test]
    fn frozen_keys_include_every_seed_holdout_arm_and_repetition() {
        assert_eq!(expected_keys(3).len(), 99);
        assert!(expected_keys(1).contains(
            &json!([
                "holdout_shifted_d4",
                0xd1b5_4a32_d192_ed03u64,
                0,
                "criterion"
            ])
            .to_string()
        ));
    }
    #[test]
    fn self_costs_exclude_inclusive_call_edges_and_resolve_late_names() {
        let text = "fn=(1) caller\n0 10\ncfn=(2)\ncalls=1 0\n0 30\n0 5\nfn=(2) callee\n0 30\n";
        let costs = self_costs(text, &BTreeMap::new()).unwrap();
        assert_eq!(costs["caller"], 15);
        assert_eq!(costs["callee"], 30);
    }
    #[test]
    fn wrapped_incremental_summary_does_not_replace_the_file_totals() {
        let text = "desc: Trigger: Client Request: [\"x\",1,0,\"f5\"]\nevents: Ir\nsummary: 18446744073709499096\ntotals: 13432\n";
        assert_eq!(parse_dump(text).unwrap().unwrap().1, 13432);
    }
}
