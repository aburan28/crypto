//! Read-only preparation of F4/F5 plot data. No benchmark execution.
use serde::Deserialize;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

#[allow(dead_code)]
#[path = "../../../src/hash/sha256.rs"]
mod sha256;

type Result<T> = std::result::Result<T, String>;

fn digest(bytes: &[u8]) -> String {
    sha256::sha256(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn read_json(path: &Path) -> Result<(Value, String)> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let value = serde_json::from_slice(&bytes).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok((value, digest(&bytes)))
}

fn text<'a>(v: &'a Value, key: &str) -> Result<&'a str> {
    v.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| format!("missing nonempty string {key}"))
}

fn number(v: &Value, key: &str) -> Result<u64> {
    v.get(key)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("missing unsigned integer {key}"))
}

fn revision(s: &str) -> Result<()> {
    if s.len() == 40 && s.bytes().all(|b| b.is_ascii_hexdigit()) {
        Ok(())
    } else {
        Err(format!("invalid full revision {s:?}"))
    }
}

fn grouped(n: u64) -> String {
    let s = n.to_string();
    s.chars()
        .enumerate()
        .fold(String::new(), |mut out, (i, c)| {
            if i > 0 && (s.len() - i).is_multiple_of(3) {
                out.push(',');
            }
            out.push(c);
            out
        })
}

#[derive(Deserialize)]
struct Spec {
    schema: String,
    groups: Vec<Group>,
}

#[derive(Deserialize)]
struct Group {
    id: String,
    title: String,
    workload: String,
    series: Vec<Series>,
}

#[derive(Deserialize)]
struct Series {
    id: String,
    label: String,
    color: String,
}

struct Snapshot {
    value: Value,
    source: String,
    sha256: String,
    utc_label: String,
    conditions: String,
    kernels: BTreeMap<String, Value>,
}

fn load_snapshots(dir: &Path) -> Result<Vec<Snapshot>> {
    let mut paths: Vec<PathBuf> = fs::read_dir(dir)
        .map_err(|e| e.to_string())?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|e| e == "json"))
        .collect();
    // The historical producer names snapshots by UTC committed timestamp, then revision.
    paths.sort();
    if paths.is_empty() {
        return Err("no history snapshots found".into());
    }
    let mut out = Vec::new();
    let mut revisions = BTreeSet::new();
    for path in paths {
        let (value, hash) = read_json(&path)?;
        if text(&value, "schema")? != "perf.history/v1" {
            return Err(format!("{}: unsupported history schema", path.display()));
        }
        let rev = text(&value, "rev")?;
        revision(rev)?;
        if !revisions.insert(rev.to_string()) {
            return Err(format!("duplicate revision {rev}"));
        }
        let file = path
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or("non-UTF8 snapshot name")?;
        let expected_suffix = format!("-{}.json", &rev[..12]);
        if !file.is_ascii()
            || file.len() != 34
            || !file.ends_with(&expected_suffix)
            || file.as_bytes()[8] != b'T'
            || file.as_bytes()[15] != b'Z'
            || !file[..8].bytes().all(|b| b.is_ascii_digit())
            || !file[9..15].bytes().all(|b| b.is_ascii_digit())
        {
            return Err(format!(
                "snapshot name must contain UTC commit timestamp and revision: {file}"
            ));
        }
        let utc_label = format!("{}-{}-{}", &file[..4], &file[4..6], &file[6..8]);
        text(&value, "committed_at")?;
        text(&value, "recorded_at")?;
        let harness = text(&value, "harness_sha256")?;
        let host = value
            .get("host")
            .filter(|h| h.is_object())
            .ok_or("snapshot host missing")?;
        let conditions = digest(
            json!({"host":host,"harness":harness,"valgrind":value.get("valgrind")})
                .to_string()
                .as_bytes(),
        );
        let rows = value
            .get("kernels")
            .and_then(Value::as_array)
            .ok_or("kernels array missing")?;
        let mut kernels = BTreeMap::new();
        for row in rows {
            let id = text(row, "id")?.to_string();
            text(row, "fingerprint")?;
            if number(row, "ir")? == 0 {
                return Err(format!("{file}: zero Ir for {id}"));
            }
            if kernels.insert(id.clone(), row.clone()).is_some() {
                return Err(format!("duplicate kernel {id}"));
            }
        }
        out.push(Snapshot {
            value,
            source: file.to_string(),
            sha256: hash,
            utc_label,
            conditions,
            kernels,
        });
    }
    Ok(out)
}

fn history_points(snaps: &[Snapshot], id: &str, metric: &str) -> Result<Vec<Value>> {
    let mut points = Vec::new();
    let mut previous: Option<(String, String)> = None;
    for (x, snap) in snaps.iter().enumerate() {
        let row = snap.kernels.get(id);
        let measured = row.and_then(|r| r.get(metric)).and_then(Value::as_u64);
        let fingerprint = row
            .and_then(|r| r.get("fingerprint"))
            .and_then(Value::as_str);
        let mut reasons = Vec::new();
        let mut display = Value::Null;
        let mut label = "missing".to_string();
        if let (Some(n), Some(fp)) = (measured, fingerprint) {
            if n == 0 {
                return Err(format!("{}: zero {metric} for {id}", snap.source));
            }
            if metric == "wall_ns" && number(&snap.value, "wall_samples")? == 0 {
                return Err(format!("{}: wall time without samples", snap.source));
            }
            if let Some((old_fp, old_conditions)) = &previous {
                if old_fp != fp {
                    reasons.push("output fingerprint changed");
                }
                if old_conditions != &snap.conditions {
                    reasons.push("host, compiler, harness or profiler changed");
                }
            } else if x > 0 {
                reasons.push("previous measurement missing");
            }
            if metric == "wall_ns" {
                let ms = n as f64 / 1_000_000.0;
                display = json!(ms);
                label = format!("{ms:.3}");
            } else {
                display = json!(n);
                label = grouped(n);
            }
            previous = Some((fp.to_string(), snap.conditions.clone()));
        } else {
            reasons.push("measurement missing");
            previous = None;
        }
        let rev = text(&snap.value, "rev")?;
        points.push(json!({
            "x":x,"date":snap.utc_label,"revision":rev,"label":format!("{}\n{}",snap.utc_label,&rev[..8]),
            "value":display,"value_label":label,"raw_value":measured,"fingerprint":fingerprint,
            "break_before":!reasons.is_empty(),"break_reasons":reasons,"source":snap.source,
            "source_sha256":snap.sha256,"conditions_sha256":snap.conditions,
            "host":snap.value["host"],"harness_sha256":snap.value["harness_sha256"],
            "wall_samples":snap.value.get("wall_samples")
        }));
    }
    Ok(points)
}

fn history_groups(spec: Spec, snaps: &[Snapshot]) -> Result<Vec<Value>> {
    if spec.schema != "groebner.groups/v1" {
        return Err("unsupported groups schema".into());
    }
    let mut groups = Vec::new();
    let mut ids = BTreeSet::new();
    for group in spec.groups {
        if !ids.insert(group.id.clone()) {
            return Err(format!("duplicate group {}", group.id));
        }
        if group.series.len() < 2 {
            return Err(format!(
                "{}: at least two algorithm series required",
                group.id
            ));
        }
        let mut charts = Vec::new();
        for (metric, suffix, title, unit) in [
            (
                "ir",
                "instructions",
                "Recorded profiler instructions",
                "Recorded Ir (instructions)",
            ),
            (
                "wall_ns",
                "runtime",
                "Historical CI runtime diagnostics",
                "Recorded median runtime (ms)",
            ),
        ] {
            let mut series = Vec::new();
            for s in &group.series {
                series.push(json!({"id":s.id,"label":s.label,"color":s.color,"points":history_points(snaps,&s.id,metric)?}));
            }
            let caption = if metric == "ir" {
                format!("{} archived revisions; recorded historical profiler scope, worker coverage incomplete.\nGaps mark changed outputs or measurement conditions. Sources: plot-data.json.", snaps.len())
            } else {
                format!("{} archived revisions; recorded medians, no stored spread or isolated timing eligibility.\nCI runtime diagnostics; sample counts and condition changes retained in plot-data.json.", snaps.len())
            };
            charts.push(json!({"id":format!("{}-{suffix}",group.id),"title":title,"unit":unit,"series":series,"caption":caption}));
        }
        groups.push(json!({
            "id":group.id,"title":group.title,"workload":group.workload,"charts":charts,
            "notes":[
                "Historical Python-orchestrated perf.history/v1 snapshots; no new performance run or speedup claim.",
                "Runtime: one recorded median per kernel and revision; sample counts retained per point, no stored spread, no isolated timing eligibility. Runtime is diagnostic.",
                "Profiler Ir is the recorded historical scope. Worker-thread coverage is incomplete for these kernels; it is not a complete solver operation count.",
                "Lines stop at missing values, changed output fingerprints, or changed host/compiler/harness/profiler. Different algorithm fingerprints are expected and are not equated.",
                "X positions are measured revision order, labelled with UTC commit date and revision; spacing is not elapsed time."
            ]
        }));
    }
    Ok(groups)
}

fn row_key(row: &Value) -> Result<String> {
    Ok(json!([
        text(row, "case")?,
        number(row, "seed")?,
        number(row, "repeat")?,
        text(row, "arm")?
    ])
    .to_string())
}

fn read_rows(path: &Path) -> Result<(BTreeMap<String, Value>, String)> {
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let contents = std::str::from_utf8(&bytes).map_err(|e| e.to_string())?;
    let mut rows = BTreeMap::new();
    for line in contents.lines().filter(|s| !s.trim().is_empty()) {
        let row: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
        let key = row_key(&row)?;
        if rows.insert(key.clone(), row).is_some() {
            return Err(format!("duplicate allocation row {key}"));
        }
    }
    if rows.is_empty() {
        return Err("empty allocation records".into());
    }
    Ok((rows, digest(&bytes)))
}

fn comparable_row(row: &Value) -> Value {
    let mut value = row.clone();
    if let Some(object) = value.as_object_mut() {
        object.remove("allocated_bytes");
        object.remove("allocations");
    }
    value
}

// Validate the original records as well as the frozen native comparison.
// A hand-edited total or a changed output must not become a plotted point.
fn validate_allocation(
    comparison: &Value,
    baseline: &BTreeMap<String, Value>,
    candidate: &BTreeMap<String, Value>,
) -> Result<(BTreeMap<String, [u64; 5]>, String)> {
    if comparison["acceptance_met"] != true || comparison["exact_outputs_and_counts"] != true {
        return Err("allocation comparison was not accepted with exact outputs".into());
    }
    let repeats = number(comparison, "repeats")?;
    let systems = number(comparison, "systems")?;
    let expected_rows = systems
        .checked_mul(repeats)
        .and_then(|n| n.checked_mul(3))
        .ok_or("allocation workload/repetition contract overflow")?;
    if repeats == 0 || systems == 0 || baseline.len() as u64 != expected_rows {
        return Err("allocation workload/repetition contract differs from raw coverage".into());
    }
    if baseline.keys().ne(candidate.keys()) {
        return Err("allocation record keys differ".into());
    }
    let mut totals: BTreeMap<String, [u64; 5]> = BTreeMap::new();
    let mut inputs = BTreeMap::new();
    let expected: BTreeMap<String, Value> = comparison["comparisons"]
        .as_array()
        .ok_or("comparisons missing")?
        .iter()
        .map(|v| Ok((text(v, "key")?.to_string(), v["counts"].clone())))
        .collect::<Result<_>>()?;
    if expected.len() != baseline.len() {
        return Err("comparison count differs from raw rows".into());
    }
    for (key, b) in baseline {
        let c = &candidate[key];
        if number(b, "repeat")? >= repeats {
            return Err("allocation repetition outside the frozen contract".into());
        }
        let arm = text(b, "arm")?;
        if !["f4", "f5", "criterion"].contains(&arm) {
            return Err(format!("unknown allocation arm {arm}"));
        }
        if arm == "criterion" {
            number(b, "word_ops")?;
            number(b, "pruned")?;
            for field in ["lower_rows", "pruned_by_part"] {
                let values = b[field]
                    .as_array()
                    .ok_or_else(|| format!("criterion {field} missing"))?;
                if values.is_empty() || values.iter().any(|v| v.as_u64().is_none()) {
                    return Err(format!("invalid criterion {field}"));
                }
            }
        } else {
            text(b, "output_digest")?;
            if arm == "f4" {
                number(b, "rank")?;
                number(b, "word_ops")?;
            } else {
                let report = b
                    .get("report")
                    .filter(|r| r.is_object())
                    .ok_or("F5 report missing")?;
                for field in [
                    "cols",
                    "criterion_rows",
                    "criterion_word_ops",
                    "degree",
                    "rank",
                    "reduce_word_ops",
                    "rows_built",
                    "rows_f4",
                    "rows_pruned",
                    "zero_reductions",
                ] {
                    number(report, field)?;
                }
            }
        }
        if comparable_row(b) != comparable_row(c) {
            return Err(format!("output or accounting differs at {key}"));
        }
        let counts = [
            number(b, "allocations")?,
            number(c, "allocations")?,
            number(b, "allocated_bytes")?,
            number(c, "allocated_bytes")?,
        ];
        if expected.get(key) != Some(&json!(counts)) {
            return Err(format!("frozen counts differ at {key}"));
        }
        let entry = totals.entry(arm.to_string()).or_insert([0; 5]);
        for (i, n) in counts.into_iter().chain([1]).enumerate() {
            entry[i] = entry[i].checked_add(n).ok_or("allocation total overflow")?;
        }
        let input_key = json!([text(b, "case")?, number(b, "seed")?]).to_string();
        let input = text(b, "input_digest")?.to_string();
        if inputs
            .insert(input_key.clone(), input.clone())
            .is_some_and(|old| old != input)
        {
            return Err(format!("algorithms used different inputs for {input_key}"));
        }
    }
    if inputs.len() as u64 != systems {
        return Err("allocation input count differs from the frozen systems".into());
    }
    let arms = comparison["arms"]
        .as_array()
        .ok_or("allocation arms missing")?;
    if arms.len() != 3 || totals.len() != 3 {
        return Err("three allocation arms required".into());
    }
    let mut arm_ids = BTreeSet::new();
    for arm in arms {
        let id = text(arm, "arm")?;
        if !arm_ids.insert(id) {
            return Err("duplicate allocation arm".into());
        }
        let want = [
            number(arm, "baseline_allocations")?,
            number(arm, "candidate_allocations")?,
            number(arm, "baseline_allocated_bytes")?,
            number(arm, "candidate_allocated_bytes")?,
            number(arm, "calls")?,
        ];
        if totals.get(id) != Some(&want) {
            return Err(format!("allocation totals differ for {id}"));
        }
    }
    let suite = digest(
        json!({"inputs":inputs,"repeats":comparison["repeats"],"systems":comparison["systems"]})
            .to_string()
            .as_bytes(),
    );
    Ok((totals, suite))
}

fn allocation_groups(root: &Path, index: &Value) -> Result<(Vec<Value>, Vec<Value>)> {
    if text(index, "schema")? != "f4-f5.allocation-history/v1" {
        return Err("unsupported allocation index".into());
    }
    let studies = index["studies"]
        .as_array()
        .ok_or("allocation studies missing")?;
    if studies.is_empty() {
        return Err("no allocation studies".into());
    }
    let mut points: BTreeMap<(String, String), Vec<Value>> = BTreeMap::new();
    let mut sources = Vec::new();
    let mut previous_suite: Option<String> = None;
    let mut workload = String::new();
    let mut allocation_scope = String::new();
    for (study_pos, study) in studies.iter().enumerate() {
        let comparison_path = text(study, "comparison")?;
        let (comparison, comparison_hash) = read_json(&root.join(comparison_path))?;
        let (baseline, baseline_hash) = read_rows(&root.join(text(study, "baseline_records")?))?;
        let (candidate, candidate_hash) = read_rows(&root.join(text(study, "candidate_records")?))?;
        let (totals, suite) = validate_allocation(&comparison, &baseline, &candidate)?;
        if study_pos == 0 {
            workload = format!("Matching frozen Boolean allocation suite: {} seeded systems, {} repetitions, {} calls per arm and revision. Polynomial dimensions, degrees and seeds are retained in the source audit protocol.",
                number(&comparison,"systems")?,number(&comparison,"repeats")?,totals["f4"][4]);
            allocation_scope = format!("Frozen native audit: {} seeded systems, {} repetitions, {} calls per arm and revision.",
                number(&comparison,"systems")?,number(&comparison,"repeats")?,totals["f4"][4]);
        }
        if previous_suite.as_ref().is_some_and(|old| old != &suite) {
            return Err(
                "allocation workload/repetition contract changed: use a separate history group"
                    .into(),
            );
        }
        previous_suite = Some(suite.clone());
        for (arm, counts) in totals {
            for (side, role) in ["baseline", "candidate"].into_iter().enumerate() {
                let rev = text(study, &format!("{role}_revision"))?;
                revision(rev)?;
                for (metric, n) in [("requests", counts[side]), ("bytes", counts[2 + side])] {
                    points.entry((arm.clone(),metric.to_string())).or_default().push(json!({
                        "x":2*study_pos+side,"date":text(study,"date")?,"revision":rev,
                        "label":format!("{}\n{role} {}",text(study,"date")?,&rev[..8]),
                        "value":n,"value_label":grouped(n),"raw_value":n,"calls":counts[4],
                        "break_before":false,"break_reasons":[],"source":comparison_path,
                        "source_sha256":comparison_hash,"run_id":text(study,&format!("{role}_run"))?,
                        "suite_sha256":suite
                    }));
                }
            }
        }
        sources.push(json!({"study":study["id"],"comparison":comparison_path,"comparison_sha256":comparison_hash,
            "baseline_records":study["baseline_records"],"baseline_sha256":baseline_hash,
            "candidate_records":study["candidate_records"],"candidate_sha256":candidate_hash,"suite_sha256":suite}));
    }
    let mut groups = Vec::new();
    for (id, title, arms) in [
        (
            "complete-allocations",
            "Complete Boolean calls: allocation history",
            vec![("f4", "F4", "#326eab"), ("f5", "F5", "#d17b35")],
        ),
        (
            "criterion-allocations",
            "F5 criterion: allocation history",
            vec![("criterion", "F5 criterion", "#7d54a3")],
        ),
    ] {
        let mut charts = Vec::new();
        for (metric, title, unit) in [
            (
                "bytes",
                "Cumulative requested allocation bytes",
                "Requested bytes (cumulative)",
            ),
            (
                "requests",
                "Allocation and reallocation requests",
                "Allocation / reallocation requests",
            ),
        ] {
            let series: Vec<Value> = arms.iter().map(|(arm,label,color)| json!({
                "id":arm,"label":label,"color":color,"points":points[&(arm.to_string(),metric.to_string())]
            })).collect();
            let caption = format!("{allocation_scope}\n{} accepted audit(s); baseline/candidate points, not a daily sequence. Sources and run IDs: plot-data.json.",studies.len());
            charts.push(
                json!({"id":format!("{id}-{metric}"),"title":title,"unit":unit,"series":series,"caption":caption}),
            );
        }
        groups.push(json!({
            "id":id,"title":title,"charts":charts,
            "workload":workload,
            "notes":[
                "Native allocation audit. Original baseline/candidate rows, accepted comparison, exact outputs and accounting, shared inputs and all totals revalidated by the Rust plot-data tool.",
                "One baseline and one candidate point per accepted audit, not a daily sequence. A/A and replication pairs are controls, not additional algorithm revisions.",
                "Bytes are cumulative requested traffic, not peak memory. No timing is measured by this allocation series.",
                "Criterion construction is a separate part of F5 and is never compared as a complete solver.",
                "Audit revisions come from their frozen records. They are not inserted as measurements of historical CI workloads; the workloads and metrics remain separate."
            ]
        }));
    }
    Ok((groups, sources))
}

fn run() -> Result<()> {
    let mut args = std::env::args().skip(1);
    if args.next().as_deref() != Some("prepare") {
        return Err("usage: groebner-history prepare --history DIR --groups FILE --allocation-index FILE --root DIR --out FILE".into());
    }
    let mut options = BTreeMap::new();
    while let Some(key) = args.next() {
        if ![
            "--history",
            "--groups",
            "--allocation-index",
            "--root",
            "--out",
        ]
        .contains(&key.as_str())
        {
            return Err(format!("unknown option {key}"));
        }
        let value = args
            .next()
            .ok_or_else(|| format!("{key}: value required"))?;
        if options.insert(key.clone(), value).is_some() {
            return Err(format!("duplicate option {key}"));
        }
    }
    let option = |key: &str| {
        options
            .get(key)
            .map(PathBuf::from)
            .ok_or_else(|| format!("{key} required"))
    };
    let snaps = load_snapshots(&option("--history")?)?;
    let (spec_value, spec_hash) = read_json(&option("--groups")?)?;
    let spec: Spec = serde_json::from_value(spec_value).map_err(|e| e.to_string())?;
    let mut groups = history_groups(spec, &snaps)?;
    let (index, index_hash) = read_json(&option("--allocation-index")?)?;
    let (allocations, allocation_sources) = allocation_groups(&option("--root")?, &index)?;
    groups.extend(allocations);
    let sources: Vec<Value> = snaps
        .iter()
        .map(|s| {
            json!({
                "file":s.source,"sha256":s.sha256,"revision":s.value["rev"],
                "committed_at":s.value["committed_at"],"recorded_at":s.value["recorded_at"],
                "host":s.value["host"],"harness_sha256":s.value["harness_sha256"],
                "valgrind":s.value["valgrind"],"wall_samples":s.value["wall_samples"]
            })
        })
        .collect();
    let output = json!({
        "schema":"groebner.plot-data/v1","producer":"tools/groebner-history (native Rust, read-only evidence preparation)",
        "snapshot_count":snaps.len(),"groups_sha256":spec_hash,"allocation_index_sha256":index_hash,
        "historical_sources":sources,"allocation_sources":allocation_sources,"groups":groups,
        "runtime_conversion":"Recorded wall_ns / 1,000,000 = milliseconds; no runtime ratio or confidence interval computed.",
        "claims":"Historical diagnostic values and native allocation counts only; no new performance measurement or speedup claim."
    });
    let out = option("--out")?;
    if let Some(parent) = out.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let mut bytes = serde_json::to_vec_pretty(&output).map_err(|e| e.to_string())?;
    bytes.push(b'\n');
    fs::write(&out, bytes).map_err(|e| e.to_string())?;
    println!(
        "Prepared {} historical snapshots and {} plot groups: {}",
        snaps.len(),
        output["groups"].as_array().unwrap().len(),
        out.display()
    );
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("groebner-history: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn snap(fp: &str, conditions: &str, wall: Option<u64>) -> Snapshot {
        Snapshot {
            value: json!({"rev":"0123456789012345678901234567890123456789","wall_samples":5,"host":{},"harness_sha256":"h"}),
            source: "fixture".into(),
            sha256: "hash".into(),
            utc_label: "2026-10-11".into(),
            conditions: conditions.into(),
            kernels: BTreeMap::from([(
                "k".into(),
                json!({"id":"k","fingerprint":fp,"ir":100,"wall_ns":wall}),
            )]),
        }
    }

    #[test]
    fn fingerprint_and_machine_changes_break_history_connections() {
        let s = [
            snap("a", "host1", Some(1_500_000)),
            snap("a", "host1", Some(1_500_000)),
            snap("b", "host1", Some(1_500_000)),
            snap("b", "host2", Some(1_500_000)),
        ];
        let points = history_points(&s, "k", "ir").unwrap();
        assert_eq!(
            points
                .iter()
                .map(|p| p["break_before"].as_bool().unwrap())
                .collect::<Vec<_>>(),
            [false, false, true, true]
        );
        assert!(points[2]["break_reasons"]
            .to_string()
            .contains("fingerprint"));
        assert!(points[3]["break_reasons"].to_string().contains("host"));
    }

    #[test]
    fn missing_runtime_stays_missing_and_breaks_the_next_connection() {
        let points = history_points(
            &[
                snap("a", "h", Some(1_500_000)),
                snap("a", "h", None),
                snap("a", "h", Some(2_500_000)),
            ],
            "k",
            "wall_ns",
        )
        .unwrap();
        assert_eq!(points[0]["value"], 1.5);
        assert!(points[1]["value"].is_null());
        assert_eq!(points[2]["value"], 2.5);
        assert_eq!(points[2]["break_before"], true);
    }

    fn allocation_fixture() -> (Value, BTreeMap<String, Value>, BTreeMap<String, Value>) {
        let mut baseline = BTreeMap::new();
        let mut candidate = BTreeMap::new();
        let mut comparisons = Vec::new();
        let mut arms = Vec::new();
        for arm in ["f4", "f5", "criterion"] {
            let mut b = json!({"arm":arm,"case":"same","seed":7,"repeat":0,"input_digest":"i","output_digest":"o","rank":2,"word_ops":7,"allocations":10,"allocated_bytes":100});
            if arm == "criterion" {
                b.as_object_mut().unwrap().remove("output_digest");
                b.as_object_mut().unwrap().remove("rank");
                b["lower_rows"] = json!([2, 3]);
                b["pruned_by_part"] = json!([0, 1]);
                b["pruned"] = json!(1);
            } else if arm == "f5" {
                b.as_object_mut().unwrap().remove("rank");
                b.as_object_mut().unwrap().remove("word_ops");
                b["report"] = json!({"cols":3,"criterion_rows":1,"criterion_word_ops":0,"degree":3,"rank":2,"reduce_word_ops":7,"rows_built":2,"rows_f4":3,"rows_pruned":1,"zero_reductions":0});
            }
            let mut c = b.clone();
            c["allocations"] = json!(5);
            c["allocated_bytes"] = json!(50);
            let key = row_key(&b).unwrap();
            comparisons.push(json!({"key":key,"counts":[10,5,100,50]}));
            arms.push(json!({"arm":arm,"baseline_allocations":10,"candidate_allocations":5,"baseline_allocated_bytes":100,"candidate_allocated_bytes":50,"calls":1}));
            baseline.insert(key.clone(), b);
            candidate.insert(key, c);
        }
        (
            json!({"acceptance_met":true,"exact_outputs_and_counts":true,"systems":1,"repeats":1,"comparisons":comparisons,"arms":arms}),
            baseline,
            candidate,
        )
    }

    #[test]
    fn allocation_outputs_and_frozen_totals_are_checked_against_raw_rows() {
        let (mut comparison, b, c) = allocation_fixture();
        assert!(validate_allocation(&comparison, &b, &c).is_ok());
        comparison["arms"][0]["candidate_allocated_bytes"] = json!(49);
        assert!(validate_allocation(&comparison, &b, &c)
            .unwrap_err()
            .contains("totals"));
    }

    #[test]
    fn changed_output_cannot_be_plotted_as_a_valid_allocation_improvement() {
        let (comparison, b, mut c) = allocation_fixture();
        c.values_mut().find(|r| r["arm"] == "f4").unwrap()["output_digest"] = json!("different");
        assert!(validate_allocation(&comparison, &b, &c)
            .unwrap_err()
            .contains("output"));
    }

    #[test]
    fn unaccepted_comparisons_and_changed_inputs_are_rejected() {
        let (mut comparison, b, c) = allocation_fixture();
        comparison["acceptance_met"] = json!(false);
        assert!(validate_allocation(&comparison, &b, &c).is_err());
        let (comparison, mut b, mut c) = allocation_fixture();
        let key = b.keys().next().unwrap().clone();
        b.get_mut(&key).unwrap()["input_digest"] = json!("other");
        c.get_mut(&key).unwrap()["input_digest"] = json!("other");
        assert!(validate_allocation(&comparison, &b, &c)
            .unwrap_err()
            .contains("different inputs"));
    }

    #[test]
    fn declared_repetitions_must_match_the_raw_record_coverage() {
        let (mut comparison, b, c) = allocation_fixture();
        comparison["repeats"] = json!(2);
        assert!(validate_allocation(&comparison, &b, &c)
            .unwrap_err()
            .contains("contract"));
    }
}
