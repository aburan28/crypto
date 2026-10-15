//! Exploratory coverage census. Preserve false flags; admit no speedup ratios.
use isogeny_algos::{json::Json, sha256::sha256_hex};
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Write as _,
    path::Path,
};
fn s<'a>(j: &'a Json, k: &str) -> Result<&'a str, String> {
    j.get(k)
        .and_then(Json::as_str)
        .ok_or_else(|| format!("missing {k}"))
}
fn census(root: &Path, label: &str) -> Result<(Json, String), String> {
    let raw = std::fs::read(root.join("family-baselines.jsonl")).map_err(|e| e.to_string())?;
    let text = std::str::from_utf8(&raw).map_err(|e| e.to_string())?;
    let mut groups: BTreeMap<String, (usize, usize, BTreeSet<String>)> = BTreeMap::new();
    let mut false_rows = vec![];
    let mut rows = 0usize;
    for (index, line) in text.lines().enumerate() {
        let j = Json::parse(line)?;
        let group = s(&j, "group")?.to_string();
        let algo = s(&j, "algo")?.to_string();
        let flag = match j.get("verified") {
            Some(Json::Bool(b)) => *b,
            _ => return Err(format!("missing boolean at row {}", index + 1)),
        };
        let entry = groups.entry(group).or_default();
        entry.0 += 1;
        entry.1 += usize::from(flag);
        entry.2.insert(algo);
        if !flag {
            false_rows.push(Json::obj(vec![
                ("line", Json::str(index + 1)),
                ("record", j),
            ]));
        }
        rows += 1;
    }
    if rows == 0 {
        return Err("empty legacy output".into());
    }
    let mut csv =
        String::from("profile,group,rows,verified_true,verified_false,distinct_algorithm_labels\n");
    let groups_json = groups
        .into_iter()
        .map(|(group, (count, verified, algos))| {
            writeln!(
                csv,
                "{label},{group},{count},{verified},{},{}",
                count - verified,
                algos.len()
            )
            .unwrap();
            Json::obj(vec![
                ("group", Json::str(group)),
                ("rows", Json::str(count)),
                ("verified_true", Json::str(verified)),
                ("verified_false", Json::str(count - verified)),
                (
                    "algorithm_labels",
                    Json::Arr(algos.into_iter().map(Json::str).collect()),
                ),
            ])
        })
        .collect();
    let isolation = root.join("family-baselines-isolation.jsonl");
    let admitted = if isolation.exists() {
        let raw = std::fs::read_to_string(&isolation).map_err(|e| e.to_string())?;
        let lines: Vec<_> = raw.lines().collect();
        if lines.len() != 1 {
            return Err("legacy isolation record count mismatch".into());
        }
        let j = Json::parse(lines[0])?;
        j.get("schema").and_then(Json::as_str) == Some("isolated-bench/1")
            && j.get("tool")
                .and_then(|t| t.get("implementation"))
                .and_then(Json::as_str)
                == Some("native")
            && j.get("run")
                .and_then(|r| r.get("exit_status"))
                .and_then(Json::as_f64)
                == Some(0.0)
            && j.get("run").and_then(|r| r.get("contended")) == Some(&Json::Bool(false))
    } else {
        false
    };
    let binary = root.join("bin/bench");
    let mut missing = vec![
        Json::str("full exact-output fingerprints"),
        Json::str("matched A/A and A/B controls"),
    ];
    let binary_hash = if binary.is_file() {
        Json::str(sha256_hex(
            &std::fs::read(binary).map_err(|e| e.to_string())?,
        ))
    } else {
        missing.push(Json::str(
            "retained legacy executable in original hosted artifacts",
        ));
        Json::Null
    };
    Ok((
        Json::obj(vec![
            ("schema", Json::str("isogeny-perf.legacy-census/v1")),
            ("status", Json::str("EXPLORATORY_ONLY")),
            ("profile", Json::str(label)),
            ("source_jsonl_sha256", Json::str(sha256_hex(&raw))),
            ("rows", Json::str(rows)),
            (
                "native_outer_isolation_record_admitted",
                Json::Bool(admitted),
            ),
            ("qualified_optimization_ratio", Json::Bool(false)),
            ("retained_legacy_binary_sha256", binary_hash),
            ("missing_evidence", Json::Arr(missing)),
            ("groups", Json::Arr(groups_json)),
            ("verified_false_rows", Json::Arr(false_rows)),
        ]),
        csv,
    ))
}
fn main() -> Result<(), String> {
    let a: Vec<_> = std::env::args().collect();
    if a.len() != 4 {
        return Err("audit-families RAW_DIR PROFILE FRESH_OUTPUT_STEM".into());
    }
    let (j, csv) = census(Path::new(&a[1]), &a[2])?;
    for (ext, data) in [("json", j.dump()), ("csv", csv)] {
        use std::io::Write;
        let mut out = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(format!("{}.{ext}", a[3]))
            .map_err(|e| e.to_string())?;
        writeln!(out, "{data}").map_err(|e| e.to_string())?;
    }
    Ok(())
}
