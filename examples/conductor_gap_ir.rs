//! Freeze, execute and independently validate a whole-solve instruction census.
//! No curve solver is implemented here: every job uses the archived ecbench.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crypto_lib::cryptanalysis::ecbench::{
    callgrind, canonical, record::ChildInput, workload::Workload,
};
use serde_json::{json, Value};

const CELLS: [&str; 2] = ["ell8209-q29", "ell131101-q37"];
const CRATERS: [&str; 2] = ["icv1-fp29-t47-14392892", "icv1-fp37-t277-91fd6e29"];
const FLOORS: [&str; 2] = ["icv1-fp29-t47-1f934956", "icv1-fp37-t277-090dc414"];
const BINARY_HASH: &str = "429de171c022b71427759b47d411d8f8d6dd52e9f2daa819de23f4a67536194a";

fn read(p: &Path) -> Result<Vec<u8>, String> {
    fs::read(p).map_err(|e| format!("{}: {e}", p.display()))
}
fn value(p: &Path) -> Result<Value, String> {
    serde_json::from_slice(&read(p)?).map_err(|e| e.to_string())
}
fn write(p: &Path, v: &Value) -> Result<(), String> {
    let mut bytes = serde_json::to_vec_pretty(v).map_err(|e| e.to_string())?;
    bytes.push(b'\n');
    fs::write(p, bytes).map_err(|e| e.to_string())
}
fn string<'a>(v: &'a Value, key: &str) -> Result<&'a str, String> {
    v[key]
        .as_str()
        .ok_or_else(|| format!("missing string {key}"))
}
fn hash(p: &Path) -> Result<String, String> {
    Ok(canonical::sha256_hex(&read(p)?))
}
fn binary_check(p: &Path) -> Result<(), String> {
    if hash(p)? != BINARY_HASH {
        return Err("ecbench differs from the frozen original executable".into());
    }
    Ok(())
}

fn freeze(study: &Path, archive: &Path, binary: &Path) -> Result<(), String> {
    binary_check(binary)?;
    let inputs = study.join("inputs");
    fs::create_dir(&inputs).map_err(|e| format!("refusing existing input directory: {e}"))?;
    let mut jobs = Vec::new();
    let mut sources = Vec::new();
    for cell in CELLS {
        let session = archive.join("sessions").join(cell);
        let spec = session.join("spec.json");
        let records = session.join("records.jsonl");
        sources.push(json!({"cell":cell,"records_sha256":hash(&records)?,"spec_sha256":hash(&spec)?,"session_sha256":hash(&session.join("session.json"))?}));
        let text = String::from_utf8(read(&records)?).map_err(|e| e.to_string())?;
        let mut controls = Vec::new();
        let mut n = 0;
        for line in text.lines() {
            let row: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
            if row["warmup"] != false {
                continue;
            }
            if row["outcome"]["status"] != "verified"
                || row["outcome"]["matches_target"] != true
                || row["cost"]["deterministic"] != true
            {
                return Err("archive has an unverified/nondeterministic measured row".into());
            }
            let seq = row["seq"].as_u64().ok_or("missing sequence")?;
            let id = format!("{cell}-seq-{seq}");
            let out = Command::new(binary)
                .args(["profile-input", "--spec"])
                .arg(&spec)
                .args(["--seq", &seq.to_string()])
                .output()
                .map_err(|e| e.to_string())?;
            if !out.status.success() {
                return Err(format!(
                    "export {id}: {}",
                    String::from_utf8_lossy(&out.stderr)
                ));
            }
            let input: Value = serde_json::from_slice(&out.stdout).map_err(|e| e.to_string())?;
            if input["expected_workload_id"] != row["workload"]["workload_id"]
                || input["expected_method_id"] != row["method"]["method_id"]
                || input["algorithm_seed"]
                    .as_u64()
                    .map(|x| x.to_string())
                    .as_deref()
                    != row["algorithm_seed"].as_str()
            {
                return Err(format!("export changed archived identities: {id}"));
            }
            let path = inputs.join(format!("{id}.json"));
            write(&path, &input)?;
            let job = json!({"id":id,"cell":cell,"seq":seq,"arm":row["arm"],"control_repeat":false,"input_sha256":hash(&path)?,"record":row});
            if job["arm"] == "ic-mitm" && job["record"]["workload"]["target_index"] == 0 {
                controls.push(job.clone());
            }
            jobs.push(job);
            n += 1;
        }
        if n != 48 || controls.len() != 2 {
            return Err(format!("unexpected frozen workload count at {cell}"));
        }
        for mut job in controls {
            job["id"] = Value::String(format!("{}-ic-aa", string(&job, "id")?));
            job["control_repeat"] = json!(true);
            let original = inputs.join(format!("{cell}-seq-{}.json", job["seq"]));
            fs::copy(
                original,
                inputs.join(format!("{}.json", string(&job, "id")?)),
            )
            .map_err(|e| e.to_string())?;
            jobs.push(job);
        }
    }
    if jobs.len() != 100 {
        return Err("expected 100 frozen profiles".into());
    }
    write(
        &study.join("jobs.json"),
        &json!({"schema":"conductor_gap.callgrind_jobs/v1","binary_sha256":BINARY_HASH,"sources":sources,"jobs":jobs}),
    )
}

fn checked_jobs(study: &Path, archive: &Path) -> Result<Vec<Value>, String> {
    let frozen = value(&study.join("jobs.json"))?;
    if frozen["schema"] != "conductor_gap.callgrind_jobs/v1"
        || frozen["binary_sha256"] != BINARY_HASH
    {
        return Err("wrong manifest contract".into());
    }
    let mut archived = BTreeMap::new();
    for src in frozen["sources"].as_array().ok_or("missing sources")? {
        let cell = string(src, "cell")?;
        let dir = archive.join("sessions").join(cell);
        for (file, key) in [
            ("records.jsonl", "records_sha256"),
            ("spec.json", "spec_sha256"),
            ("session.json", "session_sha256"),
        ] {
            if hash(&dir.join(file))? != string(src, key)? {
                return Err(format!("changed frozen source {cell}/{file}"));
            }
        }
        for line in String::from_utf8(read(&dir.join("records.jsonl"))?)
            .map_err(|e| e.to_string())?
            .lines()
        {
            let row: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
            archived.insert(
                (
                    cell.to_owned(),
                    row["seq"].as_u64().ok_or("missing archive seq")?,
                ),
                row,
            );
        }
    }
    let jobs = frozen["jobs"].as_array().ok_or("missing jobs")?.clone();
    if jobs.len() != 100 {
        return Err("missing frozen jobs".into());
    }
    let mut ids = BTreeSet::new();
    for job in &jobs {
        let id = string(job, "id")?;
        if !ids.insert(id) {
            return Err("duplicate frozen job".into());
        }
        let key = (
            string(job, "cell")?.to_owned(),
            job["seq"].as_u64().ok_or("missing job seq")?,
        );
        if archived.get(&key) != Some(&job["record"])
            || job["record"]["warmup"] != false
            || job["arm"] != job["record"]["arm"]
        {
            return Err(format!("job {id} differs from archive"));
        }
        if hash(&study.join("inputs").join(format!("{id}.json")))? != string(job, "input_sha256")? {
            return Err(format!("input hash changed: {id}"));
        }
    }
    Ok(jobs)
}

fn census(study: &Path, archive: &Path, binary: &Path, valgrind: &Path) -> Result<(), String> {
    binary_check(binary)?;
    let jobs = checked_jobs(study, archive)?;
    let raw = study.join("raw");
    fs::create_dir(&raw).map_err(|e| format!("refusing existing census directory: {e}"))?;
    for (i, job) in jobs.iter().enumerate() {
        let id = string(job, "id")?;
        let dir = raw.join(id);
        fs::create_dir(&dir).map_err(|e| e.to_string())?;
        fs::copy(
            study.join("inputs").join(format!("{id}.json")),
            dir.join("input.json"),
        )
        .map_err(|e| e.to_string())?;
        let start = std::time::Instant::now();
        let result = Command::new("timeout")
            .args(["--signal=TERM", "--kill-after=15", "1800"])
            .arg(valgrind)
            .args(["--tool=callgrind", "--cache-sim=no", "--branch-sim=no"])
            .arg(format!(
                "--callgrind-out-file={}",
                dir.join("profile").display()
            ))
            .arg(binary)
            .arg("exec")
            .env("ECBENCH_CALLGRIND_SOLVE", "1")
            .env("RAYON_NUM_THREADS", "1")
            .env_remove("ECBENCH_CALLGRIND_TARGET")
            .stdin(Stdio::from(
                File::open(dir.join("input.json")).map_err(|e| e.to_string())?,
            ))
            .stdout(Stdio::from(
                File::create(dir.join("child.json")).map_err(|e| e.to_string())?,
            ))
            .stderr(Stdio::from(
                File::create(dir.join("valgrind.log")).map_err(|e| e.to_string())?,
            ))
            .status();
        let (exit, error) = match result {
            Ok(s) => (s.code(), None),
            Err(e) => (None, Some(e.to_string())),
        };
        write(
            &dir.join("execution.json"),
            &json!({"exit_code":exit,"launch_error":error,"watchdog_seconds":1800,"elapsed_seconds_descriptive":start.elapsed().as_secs_f64()}),
        )?;
        let receipt = callgrind::solve_ir(&dir.join("profile"));
        match receipt {
            Ok(r) => write(
                &dir.join("ir.json"),
                &serde_json::to_value(r).map_err(|e| e.to_string())?,
            )?,
            Err(e) => write(&dir.join("ir-error.json"), &json!({"error":e}))?,
        }
        let mut paths = fs::read_dir(&dir)
            .map_err(|e| e.to_string())?
            .map(|e| e.map(|v| v.path()))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| e.to_string())?;
        paths.sort();
        let mut sums = String::new();
        for p in paths {
            sums.push_str(&format!(
                "{}  {}\n",
                hash(&p)?,
                p.file_name().unwrap().to_string_lossy()
            ));
        }
        fs::write(dir.join("SHA256SUMS"), sums).map_err(|e| e.to_string())?;
        println!(
            "{}/{} {} exit={:?} elapsed={:.1}s",
            i + 1,
            jobs.len(),
            id,
            exit,
            start.elapsed().as_secs_f64()
        );
        use std::io::Write;
        std::io::stdout().flush().map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn check_files(dir: &Path) -> Result<(), String> {
    let text = String::from_utf8(read(&dir.join("SHA256SUMS"))?).map_err(|e| e.to_string())?;
    let mut names = BTreeSet::new();
    for line in text.lines() {
        let (sha, name) = line.split_once("  ").ok_or("malformed hash receipt")?;
        if Path::new(name).components().count() != 1
            || !names.insert(name.to_owned())
            || hash(&dir.join(name))? != sha
        {
            return Err(format!("raw hash/filename mismatch: {name}"));
        }
    }
    for name in [
        "input.json",
        "child.json",
        "ir.json",
        "valgrind.log",
        "execution.json",
    ] {
        if !names.contains(name) {
            return Err(format!("missing raw receipt {name}"));
        }
    }
    for entry in fs::read_dir(dir).map_err(|e| e.to_string())? {
        let p = entry.map_err(|e| e.to_string())?.path();
        let n = p.file_name().unwrap().to_string_lossy();
        if n != "SHA256SUMS" && !names.contains(n.as_ref()) {
            return Err(format!("unhashed raw file {n}"));
        }
    }
    Ok(())
}

fn validate(job: &Value, dir: &Path) -> Result<u64, String> {
    check_files(dir)?;
    let input_bytes = read(&dir.join("input.json"))?;
    if canonical::sha256_hex(&input_bytes) != string(job, "input_sha256")? {
        return Err("raw input differs from frozen input".into());
    }
    let input: ChildInput = serde_json::from_slice(&input_bytes).map_err(|e| e.to_string())?;
    let child = value(&dir.join("child.json"))?;
    let old = &job["record"];
    let execution = value(&dir.join("execution.json"))?;
    if execution["exit_code"] != 0
        || !execution["launch_error"].is_null()
        || !child["error"].is_null()
        || child["report"]["exhausted"] != false
    {
        return Err("profile failed or exhausted".into());
    }
    if child["workload_id"] != old["workload"]["workload_id"]
        || child["method_id"] != old["method"]["method_id"]
        || input.expected_workload_id != string(&old["workload"], "workload_id")?
        || input.expected_method_id != string(&old["method"], "method_id")?
        || input.algorithm_seed.to_string() != string(old, "algorithm_seed")?
    {
        return Err("profile changed workload/method/seed".into());
    }
    let recovered = child["report"]["recovered"]
        .as_u64()
        .ok_or("no recovered scalar")?;
    if recovered.to_string() != string(&old["outcome"], "recovered")? {
        return Err("wrong archived scalar".into());
    }
    let (workload, instance) = Workload::build(
        &input.curve,
        input.target_seed,
        input.target_index,
        input.target_kind,
    )?;
    if workload.workload_id != input.expected_workload_id
        || instance.mul_generator_hex(u128::from(recovered)).as_ref() != Some(&workload.target)
    {
        return Err("independent point verification failed".into());
    }
    let report = &child["report"];
    if report["total_gae"] != old["cost"]["total_gae"]
        || report["counters"] != old["counters"]
        || report["unpriced"] != old["cost"]["unpriced"]
        || report["deterministic"] != true
    {
        return Err("deterministic counters changed".into());
    }
    let phases = report["phases"].as_array().ok_or("missing child phases")?;
    let old_phases = old["phases"].as_array().ok_or("missing archived phases")?;
    if phases.len() != old_phases.len() {
        return Err("phase count changed".into());
    }
    for (a, b) in phases.iter().zip(old_phases) {
        for k in ["name", "adds", "doubles", "scalar_mults", "gae", "native"] {
            if a[k] != b[k] {
                return Err(format!("phase counter changed: {k}"));
            }
        }
    }
    let raw = callgrind::solve_ir(&dir.join("profile"))?;
    let receipt = value(&dir.join("ir.json"))?;
    if serde_json::to_value(&raw).map_err(|e| e.to_string())? != receipt {
        return Err("Ir receipt differs from raw profiles".into());
    }
    Ok(raw.solve_ir)
}

fn next(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}
fn ratio(a: &[f64], b: &[f64], paired: bool) -> Value {
    let point = a.iter().sum::<f64>() / a.len() as f64 / (b.iter().sum::<f64>() / b.len() as f64);
    let mut state = 20261014;
    let mut samples = Vec::with_capacity(20000);
    for _ in 0..20000 {
        let mut sa = 0.;
        let mut sb = 0.;
        for _ in 0..a.len() {
            let i = (next(&mut state) % a.len() as u64) as usize;
            sa += a[i];
            if paired {
                sb += b[i];
            }
        }
        if !paired {
            for _ in 0..b.len() {
                sb += b[(next(&mut state) % b.len() as u64) as usize];
            }
        }
        samples.push(sa / a.len() as f64 / (sb / b.len() as f64));
    }
    samples.sort_by(f64::total_cmp);
    json!({"ratio_of_means":point,"bootstrap_95":[samples[499],samples[19499]],"bootstrap_99":[samples[99],samples[19899]],"paired":paired,"resamples":20000,"seed":20261014})
}

fn analyze(study: &Path, archive: &Path) -> Result<(), String> {
    let jobs = checked_jobs(study, archive)?;
    let mut rows = Vec::new();
    let mut all_ok = true;
    let mut groups: BTreeMap<(String, String), BTreeMap<String, u64>> = BTreeMap::new();
    for job in &jobs {
        let id = string(job, "id")?;
        let result = validate(job, &study.join("raw").join(id));
        let (ir, error) = match result {
            Ok(n) => (Some(n), None),
            Err(e) => {
                all_ok = false;
                (None, Some(e))
            }
        };
        let arm = if job["control_repeat"] == true {
            "ic-aa"
        } else {
            string(job, "arm")?
        };
        let old = &job["record"];
        rows.push(json!({"id":id,"cell":job["cell"],"arm":arm,"slug":old["workload"]["curve"]["slug"],"r":old["workload"]["curve"]["r"],"workload_id":old["workload"]["workload_id"],"target_index":old["workload"]["target_index"],"source_run_id":old["run_id"],"source_record_id":old["record_id"],"source_total_gae":old["cost"]["total_gae"],"solve_ir":ir,"error":error}));
        if let Some(n) = ir {
            let key = (
                string(&old["workload"]["curve"], "slug")?.to_owned(),
                string(&old["workload"], "workload_id")?.to_owned(),
            );
            if groups
                .entry(key)
                .or_default()
                .insert(arm.to_owned(), n)
                .is_some()
            {
                return Err("duplicate target/arm".into());
            }
        }
    }
    let mut aa_max = 0f64;
    let mut aa = Vec::new();
    for ((slug, w), arms) in &groups {
        for (base, copy) in [("rho-neg", "rho-neg-aa"), ("ic-mitm", "ic-aa")] {
            if let (Some(a), Some(b)) = (arms.get(base), arms.get(copy)) {
                let deviation = (*b as f64 / *a as f64 - 1.).abs();
                aa_max = aa_max.max(deviation);
                aa.push(json!({"slug":slug,"workload_id":w,"base":base,"copy":copy,"ratio":*b as f64 / *a as f64,"absolute_relative_deviation":deviation}));
            }
        }
    }
    let mut summaries = Vec::new();
    let mut comparisons = Vec::new();
    if all_ok {
        if aa.len() != 36 || groups.len() != 32 {
            return Err("missing A/A or primary workloads".into());
        }
        for i in 0..2 {
            for slug in [CRATERS[i], FLOORS[i]] {
                let g: Vec<_> = groups.iter().filter(|((s, _), _)| s == slug).collect();
                if g.len() != 8 {
                    return Err("missing fixed-curve targets".into());
                }
                let r = rows
                    .iter()
                    .find(|v| v["slug"] == slug)
                    .and_then(|v| v["r"].as_u64())
                    .ok_or("missing subgroup order")?;
                for arm in ["rho-neg", "rho-neg-aa", "ic-mitm"] {
                    let counts = g
                        .iter()
                        .map(|(_, a)| a.get(arm).copied().ok_or("missing primary arm"))
                        .collect::<Result<Vec<_>, _>>()?;
                    let mean = counts.iter().map(|&n| n as f64).sum::<f64>() / 8.;
                    summaries.push(json!({"cell":CELLS[i],"slug":slug,"arm":arm,"r":r,"n":8,"mean_ir":mean,"mean_ir_over_sqrt_r":mean/(r as f64).sqrt(),"min_ir":counts.iter().min(),"max_ir":counts.iter().max()}));
                }
                let ic = g
                    .iter()
                    .map(|(_, a)| a["ic-mitm"] as f64)
                    .collect::<Vec<_>>();
                let rho = g
                    .iter()
                    .map(|(_, a)| a["rho-neg"] as f64)
                    .collect::<Vec<_>>();
                comparisons.push(json!({"cell":CELLS[i],"slug":slug,"kind":"ic_over_rho","estimate":ratio(&ic,&rho,true)}));
            }
            let sample = |slug: &str| {
                groups
                    .iter()
                    .filter(|((s, _), _)| s == slug)
                    .map(|(_, a)| a["ic-mitm"] as f64)
                    .collect::<Vec<_>>()
            };
            comparisons.push(json!({"cell":CELLS[i],"kind":"floor_over_crater_ic","estimate":ratio(&sample(FLOORS[i]),&sample(CRATERS[i]),false)}));
        }
    }
    let leads = comparisons
        .iter()
        .filter(|v| v["kind"] == "floor_over_crater_ic")
        .map(|v| &v["estimate"])
        .collect::<Vec<_>>();
    let direction = |e: &Value| -> i32 {
        let p = e["ratio_of_means"].as_f64().unwrap();
        let lo = e["bootstrap_99"][0].as_f64().unwrap();
        let hi = e["bootstrap_99"][1].as_f64().unwrap();
        if lo > 1. && p > 1.1 {
            1
        } else if hi < 1. && p < 0.9 {
            -1
        } else {
            0
        }
    };
    let repeated_lead = all_ok
        && aa_max <= 0.02
        && leads.len() == 2
        && direction(leads[0]) != 0
        && direction(leads[0]) == direction(leads[1]);
    let result = json!({"schema":"conductor_gap.callgrind_analysis/v1","unit":"callgrind.Ir","expected_profiles":100,"verified_profiles":rows.iter().filter(|v|v["error"].is_null()).count(),"all_profiles_verified":all_ok,"jobs_sha256":hash(&study.join("jobs.json"))?,"aa_pairs":aa,"aa_max_deviation":aa_max,"aa_within_two_percent":all_ok && aa.len()==36 && aa_max<=0.02,"comparison_status":if !all_ok {"inadmissible_missing_or_failed_profiles"} else if aa_max>0.02 {"exploratory_aa_gate_failed"} else {"admissible_implementation_diagnostic"},"summaries":summaries,"comparisons":comparisons,"repeated_instruction_lead":repeated_lead,"rows":rows,"claim_limit":"complete simulated user-space implementation instruction census; external independent verification excluded; no native wall-time, intrinsic-hardness, asymptotic or cryptographic-scale conclusion"});
    write(&study.join("analysis.json"), &result)?;
    let mut md=String::from("# Complete implementation instruction census\n\nUnit: Callgrind simulated user-space instructions, inside `methods::solve`. External independent verification is outside this interval. No native wall-time claim.\n\n");
    md.push_str(&format!("Verified profiles: {}/100. A/A maximum relative deviation: {:.6}%. Frozen repeated instruction-lead criterion: {}.\n\n",result["verified_profiles"],aa_max*100.,repeated_lead));
    md.push_str("| Curve | Arm | n | Mean Ir | Mean Ir/sqrt(r) | Min Ir | Max Ir |\n|---|---|---:|---:|---:|---:|---:|\n");
    for s in result["summaries"].as_array().unwrap() {
        md.push_str(&format!(
            "| {} | {} | 8 | {:.0} | {:.6} | {} | {} |\n",
            string(s, "slug")?,
            string(s, "arm")?,
            s["mean_ir"].as_f64().unwrap(),
            s["mean_ir_over_sqrt_r"].as_f64().unwrap(),
            s["min_ir"],
            s["max_ir"]
        ));
    }
    md.push_str("\n| Comparison | Curve/gap | Ratio | 95% interval | 99% interval |\n|---|---|---:|---|---|\n");
    for c in result["comparisons"].as_array().unwrap() {
        let e = &c["estimate"];
        md.push_str(&format!(
            "| {} | {} | {:.6} | {:.6} to {:.6} | {:.6} to {:.6} |\n",
            string(c, "kind")?,
            c["slug"].as_str().unwrap_or(c["cell"].as_str().unwrap()),
            e["ratio_of_means"].as_f64().unwrap(),
            e["bootstrap_95"][0].as_f64().unwrap(),
            e["bootstrap_95"][1].as_f64().unwrap(),
            e["bootstrap_99"][0].as_f64().unwrap(),
            e["bootstrap_99"][1].as_f64().unwrap()
        ));
    }
    md.push_str("\nAll individual profiles, source run/record IDs and validation failures are in analysis.json. Missing profiles invalidate inference; no successful-subset analysis is used. Intervals containing one do not prove equivalence.\n");
    fs::write(study.join("analysis.md"), md).map_err(|e| e.to_string())?;
    if !all_ok {
        return Err("analysis retained failures; inference is inadmissible".into());
    }
    Ok(())
}

fn run() -> Result<(), String> {
    let a = std::env::args().collect::<Vec<_>>();
    if a.len() < 4 {
        return Err(
            "usage: conductor_gap_ir freeze|run|analyze STUDY ARCHIVE [ECBENCH [VALGRIND]]".into(),
        );
    }
    let study = PathBuf::from(&a[2]);
    let archive = PathBuf::from(&a[3]);
    match a[1].as_str() {
        "freeze" if a.len() == 5 => freeze(&study, &archive, Path::new(&a[4])),
        "run" if a.len() == 6 => census(&study, &archive, Path::new(&a[4]), Path::new(&a[5])),
        "analyze" if a.len() == 4 => analyze(&study, &archive),
        _ => Err("invalid command/argument count".into()),
    }
}
fn main() {
    if let Err(e) = run() {
        eprintln!("conductor_gap_ir: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn independent_null_and_separated_controls() {
        let null = ratio(&[12.; 8], &[12.; 8], false);
        assert_eq!(null["ratio_of_means"], 1.);
        assert_eq!(null["bootstrap_99"], json!([1., 1.]));
        let separation = ratio(&[24.; 8], &[12.; 8], false);
        assert_eq!(separation["bootstrap_99"], json!([2., 2.]));
    }
    #[test]
    fn paired_resampling_preserves_target_dependence() {
        let estimate = ratio(&[2., 8., 20., 60.], &[1., 4., 10., 30.], true);
        assert_eq!(estimate["bootstrap_99"], json!([2., 2.]));
    }
}
