//! Native reports and dashboard additions from frozen ecbench records.
//! Historical metadata is read from an exact commit, never a Python rerun.
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Write as _,
    fs,
    io::Write as _,
    path::Path,
    process::{Command, Stdio},
};

#[path = "ecbench_leads_findings.rs"]
mod lead_findings;

const BASE: &str = "a80b15bbc5c1ea4f2db2b294a9bcddc31f06fdd7";
const INCOMING_BASE: &str = BASE;
const STUDY: &str = "research/ecbench_open_questions_20261014";
const LEADS: &str = "research/ecbench_leads_20261015";

fn studies() -> Result<Vec<(String, Value)>, String> {
    let mut studies = vec![(STUDY.into(), read(Path::new(STUDY).join("summary.json"))?)];
    if Path::new(LEADS).join("summary.json").exists() {
        studies.push((LEADS.into(), read(Path::new(LEADS).join("summary.json"))?));
    }
    Ok(studies)
}

fn read(path: impl AsRef<Path>) -> Result<Value, String> {
    serde_json::from_slice(&fs::read(path).map_err(|e| e.to_string())?).map_err(|e| e.to_string())
}
fn baseline(path: &str) -> Result<String, String> {
    baseline_at(BASE, path)
}
fn baseline_at(commit: &str, path: &str) -> Result<String, String> {
    let out = Command::new("git")
        .args(["show", &format!("{commit}:{path}")])
        .output()
        .map_err(|e| e.to_string())?;
    if !out.status.success() {
        return Err(format!("cannot read frozen metadata {path}"));
    }
    String::from_utf8(out.stdout).map_err(|e| e.to_string())
}

fn catalog_derivative_append_only(path: &str) -> Result<bool, String> {
    if !matches!(
        path,
        "docs/curves/covers.json" | "docs/curves/cover-links.yaml"
    ) {
        return Ok(false);
    }
    let old: Value = serde_json::from_str(&baseline(path)?).map_err(|e| e.to_string())?;
    let actual = read(path)?;
    let additions = read("docs/curves/sources/ecbench-ten-questions.json")?;
    let allowed: BTreeSet<_> = additions["curves"]
        .as_array()
        .ok_or("missing model additions")?
        .iter()
        .map(|r| r["slug"].as_str().unwrap())
        .collect();
    let mut expected = old.clone();
    expected["registry_sha256"] = json!(sha256_hex(
        &fs::read("docs/curves/registry.json").map_err(|e| e.to_string())?
    ));
    let old_rows = old["curves"]
        .as_array()
        .ok_or("missing old catalog derivatives")?;
    let new_rows = actual["curves"]
        .as_array()
        .ok_or("missing current catalog derivatives")?;
    if !old_rows.iter().all(|row| new_rows.contains(row)) {
        return Ok(false);
    }
    let extra: Vec<_> = new_rows
        .iter()
        .filter(|row| !old_rows.contains(row))
        .collect();
    if extra.len() != allowed.len()
        || !extra.iter().all(|r| {
            allowed.contains(if path.ends_with(".json") {
                r["slug"].as_str().unwrap_or("")
            } else {
                r["icv1_identity"]["slug"].as_str().unwrap_or("")
            })
        })
    {
        return Ok(false);
    }
    expected["curves"] = actual["curves"].clone();
    if path.ends_with(".json") {
        let total = |status: &str| new_rows.iter().filter(|r| r["status"] == status).count();
        expected["summary"] = json!({"verified":total("verified_over_declared_field"),"unsupported":total("unsupported"),"invalid_input":total("invalid_input")});
    } else {
        for (kind, uid) in [("covers", "cover_uid"), ("maps", "map_uid")] {
            let a = old[kind].as_object().ok_or("missing old graph objects")?;
            let b = actual[kind].as_object().ok_or("missing graph objects")?;
            if !a.iter().all(|(k, v)| b.get(k) == Some(v)) {
                return Ok(false);
            }
            let ids: BTreeSet<_> = extra
                .iter()
                .flat_map(|r| {
                    r["cover_links"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|l| l[uid].as_str().unwrap())
                })
                .collect();
            let new_ids: BTreeSet<_> = b
                .keys()
                .filter(|k| !a.contains_key(*k))
                .map(String::as_str)
                .collect();
            if ids != new_ids {
                return Ok(false);
            }
            expected[kind] = actual[kind].clone();
        }
    }
    // Construction validity and all graph identities are separately
    // rederived by the unchanged native curve_cover_check --check in CI.
    Ok(expected == actual)
}
fn guard_historical_sources() -> Result<(), String> {
    let browser: Value =
        serde_json::from_str(&baseline("docs/browser/data.json")?).map_err(|e| e.to_string())?;
    let board: Value =
        serde_json::from_str(&baseline("docs/ic/leaderboard.json")?).map_err(|e| e.to_string())?;
    let mut files = BTreeMap::new();
    for (path, hash) in browser["sources"]
        .as_object()
        .ok_or("missing historical browser sources")?
    {
        if path != "docs/ic/leaderboard.json" {
            files.insert(
                path.clone(),
                hash.as_str().ok_or("invalid source hash")?.to_string(),
            );
        }
    }
    for source in board["sources"]
        .as_object()
        .ok_or("missing historical board sources")?
        .values()
    {
        if let (Some(path), Some(hash)) = (source["path"].as_str(), source["sha256"].as_str()) {
            files.insert(path.to_string(), hash.to_string());
        }
    }
    // Never silently lose another campaign's newly committed sessions.
    let mut indexed: BTreeSet<String> = browser["sources"]
        .as_object()
        .unwrap()
        .keys()
        .cloned()
        .collect();
    for (study, summary) in studies()? {
        for session in summary["sessions"]
            .as_array()
            .ok_or("missing campaign sessions")?
        {
            indexed.insert(format!(
                "{study}/sessions/{}/records.jsonl",
                session["name"].as_str().ok_or("missing session name")?
            ));
        }
    }
    let tracked = Command::new("git")
        .args([
            "ls-files",
            "--",
            "research/ecbench_*/sessions/*/records.jsonl",
            "research/isogeny_conductor_gap_*/sessions/*/records.jsonl",
            "research/fiber_ic_scaling_*/sessions/*/records.jsonl",
        ])
        .output()
        .map_err(|e| e.to_string())?;
    if !tracked.status.success() {
        return Err("cannot enumerate indexed session inputs".into());
    }
    for path in std::str::from_utf8(&tracked.stdout)
        .map_err(|e| e.to_string())?
        .lines()
    {
        // Git's wildcard also matches slashes; the frozen browser globs use
        // exactly research/<study>/sessions/<session>/records.jsonl.
        let components: Vec<_> = path.split('/').collect();
        if components.len() != 5 || components[2] != "sessions" {
            continue;
        }
        if !indexed.contains(path) {
            return Err(format!(
                "new session needs native metadata ingestion: {path}"
            ));
        }
    }
    for (path, expected) in files {
        let actual =
            sha256_hex(&fs::read(&path).map_err(|e| format!("historical input {path}: {e}"))?);
        if actual != expected {
            if catalog_derivative_append_only(&path)? {
                continue;
            }
            if path == "docs/curves/registry.json" {
                let mut expected: Value =
                    serde_json::from_str(&baseline(&path)?).map_err(|e| e.to_string())?;
                let additions = read("docs/curves/sources/ecbench-ten-questions.json")?;
                for row in additions["curves"]
                    .as_array()
                    .ok_or("missing registry additions")?
                {
                    if !expected["curves"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .any(|r| r["slug"] == row["slug"])
                    {
                        expected["curves"].as_array_mut().unwrap().push(row.clone());
                    }
                }
                if read(&path)? == expected {
                    continue;
                }
            }
            return Err(format!(
                "historical source changed: {path}; extend native regeneration before publication"
            ));
        }
    }
    Ok(())
}
fn pretty(v: &Value) -> String {
    let mut out = Vec::new();
    let mut s = serde_json::Serializer::with_formatter(
        &mut out,
        serde_json::ser::PrettyFormatter::with_indent(b" "),
    );
    serde::Serialize::serialize(v, &mut s).unwrap();
    String::from_utf8(out).unwrap()
}
fn whitespace(bytes: &[u8], mut i: usize) -> usize {
    while i < bytes.len() && bytes[i].is_ascii_whitespace() {
        i += 1;
    }
    i
}
fn end_string(bytes: &[u8], mut i: usize) -> usize {
    i += 1;
    while i < bytes.len() {
        match bytes[i] {
            b'\\' => i += 2,
            b'"' => return i + 1,
            _ => i += 1,
        }
    }
    i
}
fn end_value(bytes: &[u8], i: usize) -> usize {
    if bytes[i] == b'"' {
        return end_string(bytes, i);
    }
    if bytes[i] == b'[' || bytes[i] == b'{' {
        let mut depth = 1usize;
        let mut j = i + 1;
        while j < bytes.len() {
            match bytes[j] {
                b'"' => {
                    j = end_string(bytes, j);
                    continue;
                }
                b'[' | b'{' => depth += 1,
                b']' | b'}' => {
                    depth -= 1;
                    if depth == 0 {
                        return j + 1;
                    }
                }
                _ => {}
            }
            j += 1;
        }
        return j;
    }
    let mut j = i;
    while j < bytes.len() && !b",]} \r\n\t".contains(&bytes[j]) {
        j += 1;
    }
    j
}
/// Validate through serde_json; retain byte formatting of unchanged rows.
fn preserve(old: &str, new: &Value) -> String {
    fn render(old: &str, new: &Value, depth: usize) -> String {
        let previous = serde_json::from_str::<Value>(old).ok();
        if previous.as_ref() == Some(new) {
            return old.to_string();
        }
        let mut fields = Vec::<(String, &str)>::new();
        let mut elements = Vec::<&str>::new();
        if let Some(previous) = &previous {
            if previous.is_object() || previous.is_array() {
                let bytes = old.as_bytes();
                let mut i = whitespace(bytes, 0) + 1;
                let last = old.trim_end().len() - 1;
                while {
                    i = whitespace(bytes, i);
                    i < last
                } {
                    if bytes[i] == b',' {
                        i += 1;
                        continue;
                    }
                    let key = if previous.is_object() {
                        let stop = end_string(bytes, i);
                        let key = serde_json::from_str::<String>(&old[i..stop]).unwrap();
                        i = whitespace(bytes, stop) + 1;
                        i = whitespace(bytes, i);
                        Some(key)
                    } else {
                        None
                    };
                    let stop = end_value(bytes, i);
                    if let Some(key) = key {
                        fields.push((key, &old[i..stop]));
                    } else {
                        elements.push(&old[i..stop]);
                    }
                    i = stop;
                }
            }
        }
        let indent = " ".repeat(depth + 1);
        match new {
            Value::Object(map) if !map.is_empty() => {
                let mut keys: Vec<_> = fields
                    .iter()
                    .map(|(k, _)| k.clone())
                    .filter(|k| map.contains_key(k))
                    .collect();
                for key in map.keys() {
                    if !keys.contains(key) {
                        keys.push(key.clone());
                    }
                }
                let lines: Vec<_> = keys
                    .iter()
                    .map(|key| {
                        let old = fields
                            .iter()
                            .find(|(k, _)| k == key)
                            .map_or("", |(_, v)| *v);
                        format!(
                            "{indent}{}: {}",
                            serde_json::to_string(key).unwrap(),
                            render(old, &map[key], depth + 1)
                        )
                    })
                    .collect();
                format!("{{\n{}\n{}}}", lines.join(",\n"), " ".repeat(depth))
            }
            Value::Array(values) if !values.is_empty() => {
                let lines: Vec<_> = values
                    .iter()
                    .enumerate()
                    .map(|(i, v)| {
                        // Appended rows with the same object schema inherit the
                        // established field order used by catalog renderers.
                        let template = elements.get(i).copied().or_else(|| {
                            elements.last().copied().filter(|old| {
                                let Ok(Value::Object(previous)) = serde_json::from_str(old) else {
                                    return false;
                                };
                                v.as_object()
                                    .is_some_and(|next| previous.keys().eq(next.keys()))
                            })
                        });
                        format!("{indent}{}", render(template.unwrap_or(""), v, depth + 1))
                    })
                    .collect();
                format!("[\n{}\n{}]", lines.join(",\n"), " ".repeat(depth))
            }
            _ => new.to_string(),
        }
    }
    let mut output = render(old.trim_end(), new, 0);
    if old.ends_with('\n') {
        output.push('\n');
    }
    output
}
fn write_or_check(path: &Path, bytes: &[u8], check: bool) -> Result<(), String> {
    if check {
        if fs::read(path).map_err(|e| e.to_string())? != bytes {
            return Err(format!("stale native generated file {}", path.display()));
        }
    } else {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        fs::write(path, bytes).map_err(|e| e.to_string())?;
    }
    Ok(())
}
fn update_json(path: &str, value: &Value, check: bool) -> Result<(), String> {
    let out = preserve(&baseline(path)?, value);
    if serde_json::from_str::<Value>(&out).map_err(|e| e.to_string())? != *value {
        return Err(format!("native render changed JSON data: {path}"));
    }
    write_or_check(Path::new(path), out.as_bytes(), check)
}
fn catalog_alias_layout(check: bool) -> Result<(), String> {
    let path = Path::new("src/cryptanalysis/curve_aliases.json");
    let doc = read(path)?;
    let mut aliases = Vec::new();
    let mut serializer = serde_json::Serializer::with_formatter(
        &mut aliases,
        serde_json::ser::PrettyFormatter::with_indent(b""),
    );
    serde::Serialize::serialize(&doc["aliases"], &mut serializer).map_err(|e| e.to_string())?;
    let text = format!("{{\n\"schema_version\": {},\n\"generated_from\": {},\n\"generated_by\": {},\n\"aliases\": {}\n}}\n",
        doc["schema_version"],doc["generated_from"],doc["generated_by"],
        String::from_utf8(aliases).map_err(|e|e.to_string())?);
    if serde_json::from_str::<Value>(&text).map_err(|e| e.to_string())? != doc {
        return Err("alias layout changed metadata".into());
    }
    write_or_check(path, text.as_bytes(), check)
}
fn dashboard_presentation(input: &str, check: bool) -> Result<(), String> {
    // The existing source-pinned renderer owns the overview and its canonical
    // timeline embed. Render in memory so --check never mutates a working tree.
    let mut child = Command::new("node")
        .args(["tools/render_ic_dashboard.mjs", "--stdin", "--json"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| e.to_string())?;
    child
        .stdin
        .take()
        .ok_or("missing renderer input")?
        .write_all(input.as_bytes())
        .map_err(|e| e.to_string())?;
    let output = child.wait_with_output().map_err(|e| e.to_string())?;
    if !output.status.success() {
        return Err(format!(
            "dashboard presentation failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    let rendered: Value = serde_json::from_slice(&output.stdout).map_err(|e| e.to_string())?;
    write_or_check(
        Path::new("docs/index-calculus-scoreboard.html"),
        rendered["page"]
            .as_str()
            .ok_or("missing rendered page")?
            .as_bytes(),
        check,
    )?;
    write_or_check(
        Path::new("docs/ic/dashboard-overview-data.json"),
        rendered["data_text"]
            .as_str()
            .ok_or("missing rendered data")?
            .as_bytes(),
        check,
    )
}
fn html(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
fn number(v: &Value) -> String {
    v.as_f64()
        .map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| "—".into())
}
fn union(value: &mut Value, item: Value) {
    let array = value.as_array_mut().unwrap();
    if !array.contains(&item) {
        array.push(item);
    }
}
fn upsert<'a>(rows: &'a mut Value, key: &str, id: &Value, create: Value) -> &'a mut Value {
    let rows = rows.as_array_mut().unwrap();
    let index = rows.iter().position(|r| &r[key] == id).unwrap_or_else(|| {
        rows.push(create);
        rows.len() - 1
    });
    &mut rows[index]
}
fn yield_row(r: &Value) -> Value {
    let phase = |name| {
        r["phases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|p| p["name"] == name)
            .cloned()
            .unwrap_or(Value::Null)
    };
    let rel = phase("relations");
    let la = phase("linear_algebra");
    let fb = &r["factor_base"];
    let native = &rel["native"];
    let trials = native["trials"].as_f64();
    let count = native["relations"]
        .as_f64()
        .or_else(|| native["hits"].as_f64());
    let lookups = native["lookups"].as_f64();
    let ratio =
        |a: Option<f64>, b: Option<f64>| a.zip(b).and_then(|(a, b)| (b > 0.).then_some(a / b));
    let mut gae = serde_json::Map::new();
    for p in [
        "factor_base",
        "oracle_setup",
        "relations",
        "linear_algebra",
        "verify",
    ] {
        gae.insert(p.into(), phase(p)["gae"].clone());
    }
    let params = &r["method"]["params"];
    let dec = &r["detail"]["decomposition"];
    let mut result = json!({"run_id":r["run_id"],"record_id":r["record_id"],"session_id":r["session_id"],"arm":r["arm"],"round":r["round"],
        "status":r["outcome"]["status"],"curve":r["workload"]["curve"]["slug"],
        "log2_r":r["workload"]["curve"]["r"].as_f64().map(f64::log2),"target_index":r["workload"]["target_index"],"target":r["workload"]["target"],
        "method_id":r["method"]["method_id"],"fb_id":fb["fb_id"],"fb_family":fb["family"],"fb_params":fb["params"],"fb_columns":fb["columns"],"fb_signed_points":fb["signed_points"],
        "oracle":params["oracle"].as_str().or_else(||r["method"]["id"].as_str()),"solver_name":params["solver"].as_str().filter(|s|!s.is_empty()).or_else(||dec["solver"].as_str()),
        "summands":dec["summands"],"hit_rate":dec["hit_rate"],"trials":native["trials"],"relations":count,"yield":ratio(count,trials),"lookups":native["lookups"],"lookups_per_relation":ratio(lookups,count),
        "walk_steps":native["walk_steps"],"lift_failures":native["lift_failures"],"unliftable_systems":native["unliftable_systems"],"frobfold_mismatches":native["frobfold_mismatches"],
        "matrix_rows":la["native"]["rows"],"matrix_columns":la["native"]["columns"],"matrix_rank":la["native"]["rank"],"row_ops":la["native"]["row_ops"],
        "gae":gae,"total_gae":r["cost"]["total_gae"],"s":r["cost"]["s"],"deterministic":r["cost"]["deterministic"],"solver":r["solver"]});
    if r["method"]["id"] == "ic.orbit_rank_v1" {
        let c = &r["counters"];
        result["summands"] = json!(3);
        result["trials"] = c["trials"].clone();
        result["relations"] = c["point_checked_relations"].clone();
        result["yield"] = json!(ratio(
            c["point_checked_relations"].as_f64().or(Some(0.)),
            c["trials"].as_f64()
        ));
        result["lookups"] = c["table_probes_uncharged"].clone();
        result["matrix_rank"] = c["rank"].clone();
        result["matrix_columns"] = c["factor_base_columns"]
            .as_u64()
            .map(|x| json!(x + 1))
            .unwrap_or(Value::Null);
        result["matrix_rows"] = c["point_checked_relations"]
            .as_u64()
            .map_or(json!(0), |x| json!(x));
        result["row_ops"] = c["matrix_multiply_subtracts_uncharged"].clone();
        result["native_counters"] = c.clone();
        for p in r["phases"].as_array().unwrap() {
            result["gae"][p["name"].as_str().unwrap()] = p["gae"].clone();
        }
    }
    result
}
fn dashboard_browser(_root: &Path, _summary: &Value, check: bool) -> Result<(), String> {
    let mut data: Value =
        serde_json::from_str(&baseline("docs/browser/data.json")?).map_err(|e| e.to_string())?;
    let additions = read("docs/curves/sources/ecbench-ten-questions.json")?;
    for row in additions["curves"]
        .as_array()
        .ok_or("missing registered models")?
    {
        if data["curves"]
            .as_array()
            .unwrap()
            .iter()
            .any(|c| c["slug"] == row["slug"])
        {
            continue;
        }
        let rep = &row["representations"][0];
        let order: f64 = row["order"]
            .as_str()
            .unwrap()
            .parse()
            .map_err(|_| "invalid model order")?;
        let r: f64 = rep["curve"]["subgroup_order"]
            .as_str()
            .unwrap()
            .parse()
            .map_err(|_| "invalid subgroup")?;
        let p: u64 = row["params"]["p"]
            .as_str()
            .unwrap()
            .parse()
            .map_err(|_| "invalid prime")?;
        data["curves"].as_array_mut().unwrap().push(json!({"slug":row["slug"],"icv1":row["icv1"],
            "family":"prime","field":{"characteristic":"p","p":row["params"]["p"],"bits":64-p.leading_zeros(),"degree":1,"modulus":null},
            "a":row["params"]["a"],"b":row["params"]["b"],"trace":row["trace"],"order":row["order"],
            "j":row["j"],"endomorphism_discriminant":null,"construction":row["params"]["generator_call"],
            "order_bits":(order.log2()*1000.).round()/1000.,"r":rep["curve"]["subgroup_order"],"r_bits":(r.log2()*1000.).round()/1000.,
            "cofactor":rep["curve"]["cofactor"],"generator":rep["curve"]["generator"],"ec1":[rep["ec1"]],"curve_uid":[rep["curve_uid"]],
            "ec1_unresolved":null,"model_json":row["model_json"],"retired_names":row["aliases"],"standard_names":[],
            "standards_provenance":[],"sources":row["sources"],"target_group":"prime-order subgroup", "leaderboard":[],
            "on_leaderboard":false,"ecbench":[],"factor_bases":[],"tournament_cells":[],"hyperelliptic_cover":null,"cover_links":[]}));
    }
    for (study, summary) in studies()? {
        let root = Path::new(&study);
        for s in summary["sessions"].as_array().unwrap() {
            let name = s["name"].as_str().unwrap();
            let dir = root.join("sessions").join(name);
            let session = read(dir.join("session.json"))?;
            let bytes = fs::read(dir.join("records.jsonl")).map_err(|e| e.to_string())?;
            let mut status = BTreeMap::<String, u64>::new();
            let mut curves = BTreeSet::new();
            let mut arms = BTreeMap::new();
            for line in std::str::from_utf8(&bytes)
                .map_err(|e| e.to_string())?
                .lines()
            {
                let r: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
                let slug = &r["workload"]["curve"]["slug"];
                let m = &r["method"];
                let id = &m["method_id"];
                *status
                    .entry(r["outcome"]["status"].as_str().unwrap().into())
                    .or_default() += 1;
                curves.insert(slug.as_str().unwrap().to_string());
                arms.insert(
                    r["arm"].as_str().unwrap().to_string(),
                    json!({"arm":r["arm"],"role":r["role"],"method":m["id"],"method_id":id}),
                );
                let entry = upsert(
                    &mut data["methods"],
                    "method_id",
                    id,
                    json!({"method_id":id,"method_sha256":m["method_sha256"],"method":m["id"],"family":m["family"],"params":m["params"],"sessions":[],"curves":[],"runs":0,"verified":0}),
                );
                union(&mut entry["sessions"], session["session_id"].clone());
                union(&mut entry["curves"], slug.clone());
                entry["runs"] = json!(entry["runs"].as_u64().unwrap() + 1);
                if r["outcome"]["status"] == "verified" {
                    entry["verified"] = json!(entry["verified"].as_u64().unwrap() + 1);
                }
                let fb = &r["factor_base"];
                if !fb.is_null() {
                    let mut new = fb.clone();
                    new["curve"] = slug.clone();
                    new["sessions"] = json!([]);
                    new["methods"] = json!([]);
                    let f = upsert(&mut data["factor_bases"], "fb_id", &fb["fb_id"], new);
                    union(&mut f["sessions"], session["session_id"].clone());
                    union(&mut f["methods"], id.clone());
                    if let Some(c) = data["curves"]
                        .as_array_mut()
                        .unwrap()
                        .iter_mut()
                        .find(|c| &c["slug"] == slug)
                    {
                        union(&mut c["factor_bases"], fb["fb_id"].clone());
                    }
                }
                if r["warmup"] == false && m["family"] == "ic" {
                    data["yields"].as_array_mut().unwrap().push(yield_row(&r));
                }
            }
            for row in summary["rows"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|r| r["session"] == name)
            {
                let slug = &row["curve"];
                let m = arms.get(row["arm"].as_str().unwrap()).unwrap();
                if let Some(c) = data["curves"]
                    .as_array_mut()
                    .unwrap()
                    .iter_mut()
                    .find(|c| &c["slug"] == slug)
                {
                    c["ecbench"].as_array_mut().unwrap().push(json!({"session_id":session["session_id"],"arm":row["arm"],"method":m["method"],"method_id":m["method_id"],"runs":row["records"],"verified":row["verified"],"mean_s":row["mean_s"],"mean_ratio_to_floor":null,"scope":"counted-work diagnostic; unpriced native work; see campaign"}));
                }
            }
            data["sessions"].as_array_mut().unwrap().push(json!({"session_id":session["session_id"],"dir":format!("{study}/sessions/{name}"),"label":session["label"],"status":session["status"],"spec_id":session["spec_id"],"env_class_id":session["env_class_id"],"binary_sha256":session["binary_sha256"],"git_commit":session["git_commit"],"records":s["records"],"status_counts":status,"started_unix_ms":session["started_unix_ms"],"curves":curves,"arms":arms.values().collect::<Vec<_>>()}));
            data["sources"][format!("{study}/sessions/{name}/records.jsonl")] =
                json!(sha256_hex(&bytes));
        }
    }
    let covers = read("docs/curves/covers.json")?;
    let graph = read("docs/curves/cover-links.yaml")?;
    for curve in data["curves"].as_array_mut().unwrap() {
        let slug = &curve["slug"];
        curve["hyperelliptic_cover"] = covers["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| &c["slug"] == slug)
            .ok_or("missing native cover finding")?
            .clone();
        let node = graph["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|n| n["icv1_identity"]["slug"] == curve["slug"])
            .ok_or("missing native cover graph node")?;
        let links: Vec<_> = node["cover_links"]
            .as_array()
            .unwrap()
            .iter()
            .map(|link| {
                let mut value = link.clone();
                value["cover_id"] =
                    graph["covers"][link["cover_uid"].as_str().unwrap()]["cover_id"].clone();
                value["map_id"] =
                    graph["maps"][link["map_uid"].as_str().unwrap()]["map_id"].clone();
                value
            })
            .collect();
        curve["cover_links"] = json!(links);
    }
    for key in [
        "curves",
        "methods",
        "factor_bases",
        "candidates",
        "rounds",
        "sessions",
        "yields",
    ] {
        data["counts"][key] = json!(data[key].as_array().unwrap().len());
    }
    for (kind,name,what) in [
        ("oracles","pdp3-koblitz-v2 / pdp3-binary-v2","Versioned complete chained-system shape; solver calls and old adapter counts preserved."),
        ("oracles","ic.orbit_rank_v1","Width-generic exact direct or signed-Frobenius quotient pair tables, full-point witness checks and modular target pinning; native work unpriced."),
        ("factor_base_families","projected-raw-x-orbits-v1","Raw-x lifts projected by the cofactor, retaining distinct signed Frobenius orbits with unknown logarithms.")
    ] {
        if !data["vocabulary"][kind].as_array().unwrap().iter().any(|v|v["name"]==name) {
            data["vocabulary"][kind].as_array_mut().unwrap().push(json!({"name":name,"where":"native ecbench","what":what}));
        }
    }
    data["sources"]["docs/ic/leaderboard.json"] = json!(sha256_hex(
        &fs::read("docs/ic/leaderboard.json").map_err(|e| e.to_string())?
    ));
    data["sources"]["docs/curves/registry.json"] = json!(sha256_hex(
        &fs::read("docs/curves/registry.json").map_err(|e| e.to_string())?
    ));
    for path in ["docs/curves/covers.json", "docs/curves/cover-links.yaml"] {
        data["sources"][path] = json!(sha256_hex(&fs::read(path).map_err(|e| e.to_string())?));
    }
    data["generated_by"]=json!("examples/ecbench_open_questions_publish.rs (native augmentation of frozen historical metadata)");
    data["native_campaign"] = json!({"baseline_commit":BASE,"incoming_main_commit":INCOMING_BASE,"scope":"all campaign/follow-up sessions and incoming-main sessions indexed; historical scientific analyses retained"});
    update_json("docs/browser/data.json", &data, check)
}
fn table(rows: &[Value]) -> String {
    let mut out=String::from("| Curve | Arm | Verified / attempts | Mean cold S | Counted IC / rho | Status |\n|---|---|---:|---:|---:|---|\n");
    for r in rows {
        let _ = writeln!(
            out,
            "| {} | {} | {} / {} | {} | {} | {} |",
            r["curve"].as_str().unwrap(),
            r["arm"].as_str().unwrap(),
            r["verified"],
            r["records"],
            number(&r["mean_s"]),
            number(&r["paired_complete_counted_ratio_candidate_over_rho"]),
            r["statuses"]
        );
    }
    out
}
fn section(rows: &[Value], id: &str, title: &str, intro: &str, study: &str) -> String {
    let mut out=format!("<section id=\"{}\"><h2>{}</h2><p>{} Count quotients do not bound full-cost ratios. L0 wall measurements are inadmissible. <a href=\"../{study}/REPORT.md\">Report, sources and remaining obligations</a>.</p><div class=\"scroll\"><table><thead><tr><th>Curve</th><th>Arm</th><th>Verified / attempts</th><th>Mean cold S</th><th>Counted IC / rho</th><th>Status</th></tr></thead><tbody>",html(id),html(title),html(intro));
    for r in rows {
        let _ = write!(
            out,
            "<tr><td>{}</td><td>{}</td><td>{}/{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
            html(r["curve"].as_str().unwrap()),
            html(r["arm"].as_str().unwrap()),
            r["verified"],
            r["records"],
            number(&r["mean_s"]),
            number(&r["paired_complete_counted_ratio_candidate_over_rho"]),
            html(&r["statuses"].to_string())
        );
    }
    out.push_str("</tbody></table></div></section>\n");
    out.replace(
        "class=\"scroll\"",
        "class=\"table-scroll\" style=\"overflow-x:auto\"",
    )
}

fn generate(check: bool) -> Result<(), String> {
    guard_historical_sources()?;
    catalog_alias_layout(check)?;
    let root = Path::new(STUDY);
    let summary = read(root.join("summary.json"))?;
    let literature = read(root.join("literature/questions.json"))?;
    let controls = read(root.join("controls-v3.json"))?;
    let failed = read(root.join("failed-replays-v2.json"))?;
    if failed["ok"] != true
        || summary["sessions"]
            .as_array()
            .unwrap()
            .iter()
            .any(|s| s["audit_ok"] != true)
    {
        return Err("all replay gates must pass before publication".into());
    }
    let mut attempts = 0u64;
    let mut verified = 0u64;
    let mut statuses = BTreeMap::<String, u64>::new();
    for r in summary["rows"].as_array().unwrap() {
        attempts += r["records"].as_u64().unwrap();
        verified += r["verified"].as_u64().unwrap();
        for (k, v) in r["statuses"].as_object().unwrap() {
            *statuses.entry(k.clone()).or_default() += v.as_u64().unwrap();
        }
    }
    let remaining=[
        "A scalable prime-field base with executable cheap membership/decomposition and a fully priced held-out crossover.",
        "Growing-dimension solving-degree bounds; correct chained-system shape metadata and held-out satisfiable/infeasible queries.",
        "Compact algebraic witnesses for relative Frobenius phases, including membership, splitting and lift costs.",
        "Scalable nonlinear algebraic membership/solving; gamma weight is not a normal-basis Hamming ideal.",
        "A useful prime-degree parameter family, subgroup yield law and complete measured IC integration.",
        "Growing arity with full splitting/lifting/matrix costs; literal ten-summand toy results do not transfer to N83.",
        "Rank and descent reliability laws, optimized large-prime policies and fully priced reusable setup/online windows.",
        "Explicit paid isogeny transport and a Hecke-rank theorem or certified counterexample in its stated model.",
        "Independent exact source-order proofs, a general Jacobian norm and a cheaper higher-genus relation/descent solver.",
        "Wide IC at N83, N131 transfer justification, isolated A/A, independent host, complete pricing, held-out scaling and matched GPU/annealing hardware.",
    ];
    let experiments = [
        vec!["q01-prime"],
        vec!["q02-binary", "q02-solvers"],
        vec!["q03-frobenius"],
        vec!["q04-nonlinear"],
        vec![],
        vec!["q06-many-summand"],
        vec!["q07-rank", "q06-many-summand"],
        vec!["q08-isogeny"],
        vec![],
        vec!["q10-n83-gate"],
    ];
    let mut requirements = Vec::new();
    for (i, q) in literature["questions"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
    {
        requirements.push(json!({"question":i+1,"title":q["title"],"original_question":q["question"],"original_resolution":q["resolution"],"sources":q["sources"],"status":"partial","diagnostic_sessions":experiments[i],"remaining":remaining[i],"mathematical_question_resolved":false}));
    }
    let results = json!({"schema":"ecbench-ten-questions-results/v1","attempts":attempts,"verified_recoveries":verified,"statuses":statuses,
        "standard_audits":summary["sessions"].as_array().unwrap().len(),"standard_verified_replays":verified,
        "supplemental_failed_replays":failed["replays"].as_array().unwrap().len(),"supplemental_failed_replays_ok":failed["ok"],
        "requirements":requirements,"controls":controls,"rows":summary["rows"],
        "claims":{"wall_time":false,"fully_priced_speedup":false,"asymptotic_advance":false,"independent_host":false,"gpu":false,"quantum_hardware":false}});
    write_or_check(
        &root.join("results.json"),
        format!("{}\n", pretty(&results)).as_bytes(),
        check,
    )?;
    let mut report=format!("# All ten questions: ecbench execution and remaining obligations\n\nExecuted nine ecbench sessions with **{attempts} attempts**, **{verified} verified scalar recoveries**, **{} exhausted searches** and **{} applicability error**. Nine session audits passed; all {verified} recoveries replayed exactly, and the {} failure outcomes/counts passed a separate native replay. Every original mathematical question remains open.\n\nThe native harness is aburan28/crypto; aburan28/ecbench is empty. [Frozen protocol](PROTOCOL-v1.md), [structural amendment](CONTROLS-AMENDMENT-v2.md), [extended-control/replay amendment](CONTROLS-AMENDMENT-v3.md), [source manifest](SOURCE-MANIFEST.json), [commands](execution-manifest.json), [summary](summary.json), [complete outcomes](results.json). Measured binary cc5bba39d33465ddc1b1d31af8ac25c6680d41193fc76012d9a1019b6046f02d, compiled from clean commit c716e462934128c75e99dd458748e08623a2708c. Host ECBENV2h687984d4e73d, virtualized AMD EPYC 9V45, earned L0. Public target law hash_to_subgroup_v1; seed20261014.\n\nAll comparisons use one public target per workload with paired rho and A/A controls. Setup is charged in every cold execution. Primary online windows and their exclusive phases are retained in the raw records, but their L0 clocks are inadmissible for a performance claim. Native work remains unpriced. S is counted group-addition equivalents / sqrt(r); its quotient is a counted-work diagnostic, **not a bound on the ratio of complete costs**. The harness floor is a birthday-style reference constant. Successful-only mean S excludes failures; paired ratios are null whenever any execution in that curve/arm is censored. All failures retain their recorded work and resource vector. No scaling exponent was fitted.\n\n",statuses.get("exhausted").unwrap_or(&0),statuses.get("error").unwrap_or(&0),failed["replays"].as_array().unwrap().len());
    report.push_str("## Results by original question\n\n");
    for (i, q) in literature["questions"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
    {
        let evidence = if experiments[i].is_empty() {
            if i == 4 {
                "controls-v3.json / quasi_subfield".into()
            } else {
                "controls-v3.json / covers; prime-cover-replay.json".into()
            }
        } else {
            experiments[i].join(", ")
        };
        let citations = q["sources"]
            .as_array()
            .unwrap()
            .iter()
            .map(|id| {
                let s = literature["sources"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|s| s["id"] == *id)
                    .unwrap();
                format!("[{}]({})", id.as_str().unwrap(), s["url"].as_str().unwrap())
            })
            .collect::<Vec<_>>()
            .join("; ");
        let _=writeln!(report,"{}. **{}** — partial.\n\n   Evidence: {}. Remaining: {}\n\n   Literature/repo basis: {}\n",i+1,q["title"].as_str().unwrap(),evidence,remaining[i],citations);
    }
    report.push_str("## Structural observations\n\nThe complete quasi-subfield coefficient searches examined32,16384,2097152 and2048 vectors at (n,n0)=(5,2),(7,3),(7,4),(11,3), finding0,127,127 and0 splitting polynomials. Retained witnesses passed independent root enumeration; their actual point/subgroup support is recorded. This finite census is not a general family or complete IC attack.\n\nAll six gamma-weight base constructions passed on-curve, subgroup-annihilation, Frobenius-closure and coefficient checks. The cap is eight projected orbits; raw weight domains and projection counts are retained. Hamming ideals for syndrome decoding do not automatically apply to this base.\n\nThe degree controls measure m2 Macaulay matrices through degree5, at dimensions3/4, and exhaust their Boolean assignments. They also record actual m3 chain sizes. Ecbench PDP3 shape currently reports only summand coordinates; intermediate field variables and the second equation block must be included before interpreting solver shape/growth. First fall degree is not solving degree. Pipeline trial-cap exhaustion does not by itself imply an individual algebraic query hit its node budget.\n\nThe degree7 and exact degree83 E0/E1 genus3 pullbacks passed generator/image annihilation, sampled scalar homomorphism, typed push-pull, infinity, ramified and off-curve controls. An additional degree7 public-point Jacobian BSGS control recovered the scalar and replayed it on the elliptic curve. This is functional generic recovery, with no measured speed claim. The F101 genus2 cover certificate replay passed; two-infinity Jacobian arithmetic/pullback/norm remain unsupported. E1 degree83 uses catalogued r=8569786107849059 and cofactor1128547018; it does not replace the E0 N83 subgroup. The first E1 cofactor2 constructor failure is preserved.\n\nThe QUBO encodes the repository's descended Semaev system. All8192 original assignments passed faithful monomial/slack completion checks; a small independent fixture exhausts every auxiliary assignment. This does not enumerate every auxiliary assignment of the full encoding, prove an x-only zero-energy solution lifts to a point relation, or run hardware annealing.\n\n## Session ledger and counted costs\n\n| Session | Attempts | Verified | Exact verified replays | Audit |\n|---|---:|---:|---:|---|\n");
    let mut render_sessions = Vec::new();
    for s in summary["sessions"].as_array().unwrap() {
        let name = s["name"].as_str().unwrap();
        let _ = writeln!(
            report,
            "| [{}](sessions/{name}/session.json) | {} | {} | {} | {} |",
            name, s["records"], s["verified_records"], s["replays"], s["audit_ok"]
        );
        let session = read(root.join("sessions").join(name).join("session.json"))?;
        let rows: Vec<_> = summary["rows"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|r| r["session"] == name)
            .cloned()
            .collect();
        let epoch = session["started_unix_ms"].as_u64().unwrap_or(0) / 1000;
        let minute = epoch / 60 % 60;
        let hour = epoch / 3600 % 24;
        let path = format!("reports/run-20261014T{hour:02}{minute:02}Z-{name}");
        let run_report=format!("# {name}: native ecbench run\n\nSession {}, spec{}, source{}, binary{}. Host {}, L0; wall clocks inadmissible. Public single-target workloads, seeds20261014, no planted log, zero warmup and cold setup per attempt. [Frozen spec](../../specs/{name}.json), [raw session](../../sessions/{name}/session.json), [audit](../../{name}.audit.json), [resource vector](../../{name}.resources.json).\n\n{} attempts, {} verified; {} standard exact replays. Failed deterministic outcomes have a separate [supplemental replay](../../failed-replays.json). All failures are retained. The N83 applicability error is an implementation gap, not a mathematical obstruction.\n\nMean/min/max S use verified recoveries only. A partial count quotient does not bound full cost. No speed, wall, exponent or hardware claim.\n\n{}\n\n![Recorded S and success ranges](cost.svg)\n",s["session_id"],session["spec_id"],session["git_commit"],session["binary_sha256"],session["env_class_id"],s["records"],s["verified_records"],s["replays"],table(&rows));
        let run_report =
            run_report.replace("../../failed-replays.json", "../../failed-replays-v2.json");
        write_or_check(
            &root.join(&path).join("REPORT.md"),
            run_report.as_bytes(),
            check,
        )?;
        render_sessions
            .push(json!({"session":name,"path":path,"metadata":session,"ledger":s,"rows":rows}));
    }
    report.push('\n');
    report.push_str(&table(summary["rows"].as_array().unwrap()));
    report.push_str("\n## Evidence limits and delivery\n\nAll ten diagnostic tracks were attempted; all ten original resolution criteria remain partial. No IC speedup, new asymptotic result, exact universal lower-bound crossing or new mathematical theorem is established. N83 IC was actually attempted and returned the missing-wide-implementation error. Degree31 has ord_31(2)=5 and cannot replace degree83/131 structure. Independent-host/L2/GPU/annealer gates remain blocked by missing configured runner access; [access problem](../../problems/ecbench-ten-question-runner-access.md), [open resolution gates](../../problems/ecbench-ten-question-resolution-gates.md). Local unreviewed evidence is not a merged or independently accepted result.\n\nBuild the native ecbench_open_questions, ecbench_open_questions_run, ecbench_open_questions_replay_failures and ecbench_open_questions_publish examples plus the ecbench binary. Freeze a **new directory**, run the native dispatcher, audit/supplement failures, analyze, and render the recorded output. Never reuse or overwrite a sealed session. Exact commands are in execution-manifest.json and COMMANDS.md.\n");
    write_or_check(&root.join("REPORT.md"), report.as_bytes(), check)?;
    write_or_check(
        &root.join("render-data.json"),
        format!(
            "{}\n",
            pretty(&json!({"sessions":render_sessions,"results":results,"literature":literature}))
        )
        .as_bytes(),
        check,
    )?;
    let rows = summary["rows"].as_array().unwrap();
    let mut block=section(rows,"ecbench-ten-questions-20261014","Ten open questions: native ecbench diagnostics, 2026-10-14","937 attempts; historical counted pipeline refused wide IC. All original resolutions remain partial.",STUDY);
    for (study, summary) in studies()?.into_iter().filter(|(p, _)| p != STUDY) {
        let attempts: u64 = summary["sessions"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s["records"].as_u64().unwrap())
            .sum();
        block.push_str(&section(summary["rows"].as_array().unwrap(),"ecbench-leads-20261015","Pursuing the leads: versioned shape and width-generic IC, 2026-10-15",
            &format!("{attempts} fresh public-target attempts; exact E0 degree83 bounded IC is executable. Successful N83 recovery and the mathematical gates remain open."),&study));
    }
    let output = baseline_at(INCOMING_BASE, "docs/index-calculus-scoreboard.html")?.replacen(
        "</body>",
        &format!("{block}</body>"),
        1,
    );
    write_or_check(
        Path::new("docs/ic/ecbench-diagnostics.fragment.html"),
        block.as_bytes(),
        check,
    )?;
    let mut leaderboard: Value =
        serde_json::from_str(&baseline("docs/ic/leaderboard.json")?).map_err(|e| e.to_string())?;
    let additions = read("docs/curves/sources/ecbench-ten-questions.json")?;
    for row in additions["curves"].as_array().unwrap() {
        if leaderboard["roster"]
            .as_array()
            .unwrap()
            .iter()
            .any(|r| r["slug"] == row["slug"])
        {
            continue;
        }
        let p = &row["params"];
        let order: u64 = row["order"]
            .as_str()
            .unwrap()
            .parse()
            .map_err(|_| "bad registry order")?;
        leaderboard["roster"].as_array_mut().unwrap().push(json!({"slug":row["slug"],"family":"prime",
            "field":format!("GF(p), p = {}",p["p"].as_str().unwrap()),
            "coefficients":format!("a = {}, b = {}",p["a"].as_str().unwrap(),p["b"].as_str().unwrap()),
            "order_bits":64-order.leading_zeros(),"standard":row["standard_names"],
            "ec1":row["representations"].as_array().unwrap().iter().map(|r|r["ec1"].clone()).collect::<Vec<_>>(),
            "ec1_unresolved":null,"legacy":row["aliases"],"on_board":false,
            "curves_yaml_key":null,"factor_base_link_status":"no_curves_yaml_record"}));
    }
    leaderboard["diagnostic_campaigns"] = json!(studies()?.iter().map(|(p,s)|Ok(json!({"id":p,"scope":"partial counted-work diagnostics, no fully priced ratios","source":format!("{p}/summary.json"),"sha256":sha256_hex(&fs::read(Path::new(p).join("summary.json")).map_err(|e|e.to_string())?),"rows":s["rows"]}))).collect::<Result<Vec<Value>,String>>()?);
    leaderboard["sources"]["registry"]["sha256"] = json!(sha256_hex(
        &fs::read("docs/curves/registry.json").map_err(|e| e.to_string())?
    ));
    leaderboard["native_augmentation"] = json!({"generator":"examples/ecbench_open_questions_publish.rs","baseline_commit":BASE,"historical_analyses_recomputed":false});
    update_json("docs/ic/leaderboard.json", &leaderboard, check)?;
    let mut timeline: Value = serde_json::from_str(&baseline("docs/ic/progress-timeline.json")?)
        .map_err(|e| e.to_string())?;
    let points:Vec<_>=rows.iter().filter(|r|r["arm"]!="rho" && r["arm"]!="rho-aa" && !r["paired_complete_counted_ratio_candidate_over_rho"].is_null()).map(|r|json!({"date":"2026-10-14","value":r["paired_complete_counted_ratio_candidate_over_rho"],"reference":"matched","label":format!("{} / {}",r["curve"].as_str().unwrap(),r["arm"].as_str().unwrap()),"detail":"paired counted-work quotient only; no full-cost bound or wall result; failed arms remain in the linked ledger","class":"accounting","anchor":"#ecbench-ten-questions-20261014","source":format!("https://github.com/aburan28/crypto/blob/research/ecbench-ten-questions-20261014/{STUDY}/REPORT.md")})).collect();
    timeline["series"].as_array_mut().unwrap().push(json!({"id":"ec-ten-question-counts","label":"Ten-question campaign: counted diagnostics","role":"diagnostic","barrier_test":false,"points":points}));
    for (study, summary) in studies()?.into_iter().filter(|(p, _)| p != STUDY) {
        let points:Vec<_>=summary["rows"].as_array().unwrap().iter().filter(|r|!r["arm"].as_str().unwrap().starts_with("rho") && !r["paired_complete_counted_ratio_candidate_over_rho"].is_null())
            .map(|r|json!({"date":"2026-10-15","value":r["paired_complete_counted_ratio_candidate_over_rho"],"reference":"matched","label":format!("{} / {}",r["curve"].as_str().unwrap(),r["arm"].as_str().unwrap()),"detail":"bounded counted quotient only; full costs and wall result unknown","class":"accounting","anchor":"#ecbench-leads-20261015","source":format!("https://github.com/aburan28/crypto/blob/research/ecbench-ten-questions-20261014/{study}/REPORT.md")})).collect();
        timeline["series"].as_array_mut().unwrap().push(json!({"id":"ec-ten-question-leads","label":"Follow-up: counted diagnostics","role":"diagnostic","barrier_test":false,"points":points}));
    }
    timeline["updated"] = json!(if Path::new(LEADS).join("summary.json").exists() {
        "2026-10-15"
    } else {
        "2026-10-14"
    });
    update_json("docs/ic/progress-timeline.json", &timeline, check)?;
    dashboard_presentation(&output, check)?;
    dashboard_browser(root, &summary, check)?;
    Ok(())
}
pub fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args == ["--followup-findings"] {
        return lead_findings::main();
    }
    if args == ["--historical-inputs"] {
        return guard_historical_sources();
    }
    if !args.is_empty() && args != ["--check"] {
        return Err("usage: ecbench_open_questions_publish [--check|--historical-inputs]".into());
    }
    generate(!args.is_empty())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn preserve_unrelated_arrays_and_string_escapes() {
        let old="{\n \"historic\": [{\"x\": \"[}\\\"\"}],\n \"rows\": [{\"n\": 1}],\n \"counts\": {\"rows\":1}\n}\n";
        let mut value: Value = serde_json::from_str(old).unwrap();
        value["rows"].as_array_mut().unwrap().push(json!({"n":2}));
        value["counts"]["rows"] = json!(2);
        value["native"] = json!(true);
        let output = preserve(old, &value);
        assert!(output.contains("\"historic\": [{\"x\": \"[}\\\"\"}]"));
        assert_eq!(serde_json::from_str::<Value>(&output).unwrap(), value);
    }
}
