//! Prime Degree Frontier: fresh exact replay, immutable rules and separate curve scores.
use super::{catalog, is_prime, provenance, replay, sha256_hex, supervisor, Options};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    fs,
    io::Write,
    path::{Path, PathBuf},
};

const RULES: &[u8] = include_bytes!("../challenge.json");
fn read_bytes(path: &Path, maximum: u64) -> Result<Vec<u8>, String> {
    let meta = fs::symlink_metadata(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if !meta.is_file() || meta.len() > maximum {
        return Err("expected a regular file within the published byte limit".into());
    }
    // Bound the read itself as well as its initial metadata.
    use std::io::Read;
    let mut bytes = Vec::new();
    fs::File::open(path)
        .map_err(|e| e.to_string())?
        .take(maximum + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    if bytes.len() as u64 > maximum {
        return Err("file grew beyond its byte limit".into());
    }
    Ok(bytes)
}
fn read_json(path: &Path, maximum: u64) -> Result<(Vec<u8>, Value), String> {
    let bytes = read_bytes(path, maximum)?;
    let value = serde_json::from_slice(&bytes).map_err(|e| format!("invalid JSON: {e}"))?;
    Ok((bytes, value))
}
fn save(path: &Path, value: &Value) -> Result<(), String> {
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| e.to_string())?;
    writeln!(
        file,
        "{}",
        serde_json::to_string_pretty(value).map_err(|e| e.to_string())?
    )
    .map_err(|e| e.to_string())
}
fn rules() -> Value {
    serde_json::from_slice(RULES).unwrap()
}
fn season() -> String {
    sha256_hex(format!("{}:{}", sha256_hex(RULES), env!("ISOGENY_REGISTRY_SHA256")).as_bytes())
}
fn standard_prime(row: &Value) -> bool {
    row["family"] == "prime"
        && row["standard_names"]
            .as_array()
            .is_some_and(|names| !names.is_empty())
}
fn curves() -> Vec<Value> {
    catalog::registry()["curves"].as_array().unwrap().iter().filter(|r| standard_prime(r)).map(|r| {
        json!({"source_icv1":r["slug"],"standard_names":r["standard_names"],"aliases":r["aliases"],
            "capabilities":catalog::capabilities(r),"best_prime_degree":Value::Null,"status":"UNSCORED"})
    }).collect()
}
fn policy() -> Value {
    let c = curves();
    json!({"schema":"isogeny-challenge-policy/v1","status":"PASS","season":season(),"rules":rules(),
        "rules_sha256":sha256_hex(RULES),"catalogue_sha256":env!("ISOGENY_REGISTRY_SHA256"),
        "tool":provenance(),"curves":c,"curve_count":c.len(),
        "degree_limit_note":"the bundled search adapter stops at 1000000; larger full-map submissions may be replayed within verifier and artifact limits; degree encoding is u64"})
}
fn resolve(name: &str) -> Result<&'static Value, String> {
    let row = catalog::resolve(name)?;
    if !standard_prime(row) {
        return Err("the main track requires a registered standard prime-field curve".into());
    }
    Ok(row)
}
fn plan(o: &Options) -> Result<Value, String> {
    let lower = o.number("from", Some(10454))?;
    let upper = o.number("to", Some(65537))?;
    let max_order = o.number("max-order", Some(16))?;
    if lower < 3 || upper < lower || upper > 1_000_000 || max_order == 0 || max_order > 16 {
        return Err("plan needs 3 <= from <= to <= 1000000 and 1 <= max-order <= 16".into());
    }
    let selected: Vec<_> = if o.get("curve").unwrap_or("all") == "all" {
        catalog::registry()["curves"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|r| standard_prime(r))
            .collect()
    } else {
        vec![resolve(o.need("curve")?)?]
    };
    let primes: Vec<_> = (lower..=upper).rev().filter(|&n| is_prime(n)).collect();
    let mut rows = Vec::new();
    for row in selected {
        let caps = catalog::capabilities(row);
        let mut candidates = Vec::new();
        if !row["order"].is_null() {
            let q = catalog::field_order(row)?;
            let t = &(&q + &super::Int::one()) - &catalog::integer(&row["order"]);
            for &ell in &primes {
                if let Some(roots) = super::eigenvalues_from_trace(&q, &t, ell) {
                    let orders = roots.map(|a| super::multiplicative_order(a, ell));
                    if *orders.iter().min().unwrap() <= max_order {
                        let name = catalog::label(row);
                        let kernel = caps["construction"] == true
                            && ell > 1009
                            && super::constructor::selection(name, ell).is_some();
                        candidates.push(json!({"ell":ell.to_string(),"eigenvalues":roots,"eigenvalue_orders":orders,
                            "constructor_route":if caps["construction"] != true {"UNAVAILABLE"} else if kernel {"kernel"} else {"modpoly"},"proof_status":"CANDIDATE_ONLY"}));
                    }
                }
            }
        }
        rows.push(json!({"source_icv1":row["slug"],"standard_names":row["standard_names"],"capabilities":caps,
            "screen_status":if row["order"].is_null() {"MISSING_ORDER"} else {"SCREENED"},"candidates":candidates}));
    }
    Ok(
        json!({"schema":"isogeny-challenge-plan/v1","status":"PASS","season":season(),"lower":lower,"upper":upper,
        "max_order":max_order,"ordering":"descending degree within each curve","curves":rows,
        "scope":"small eigenvalue-order shortlist, not exhaustive existence testing or constructed maps"}),
    )
}
struct Submission {
    label: String,
    slug: String,
    ell: u64,
    summary: Vec<u8>,
    map: Vec<u8>,
    stderr: Vec<u8>,
}
fn preflight(input: &Path) -> Result<Submission, String> {
    let (summary_bytes, summary) = read_json(&input.join("search.json"), 1 << 20)?;
    let attempts = summary["attempts"].as_array().ok_or("missing attempts")?;
    if attempts.len() != 1 {
        return Err("submit exactly one source and prime degree per directory".into());
    }
    let a = &attempts[0];
    let label = a["curve"].as_str().ok_or("missing curve")?;
    let source = resolve(label)?;
    if label != catalog::label(source) {
        return Err("use the constructor's canonical curve label".into());
    }
    if catalog::capabilities(source)["independent_replay"] != true {
        return Err("this curve is listed but its replay adapter is unavailable".into());
    }
    let ell = a["ell"]
        .as_u64()
        .ok_or("degree must be an exact unsigned integer")?;
    if ell < 3 || !is_prime(ell) {
        return Err(
            "only odd prime degrees score; chain products and composites are excluded".into(),
        );
    }
    let name = format!("ell-{ell}.json");
    if a["receipt"]["stdout"] != name
        || a["receipt"]["status"] != "PASS"
        || a["receipt"]["exit_code"] != 0
    {
        return Err("missing completed construction receipt".into());
    }
    let directory = input.join(label);
    let m = fs::symlink_metadata(&directory).map_err(|e| e.to_string())?;
    if !m.is_dir() {
        return Err("curve directory must be a real directory".into());
    }
    // The parent copies bounded raw bytes; potentially expensive JSON parsing
    // and algebra happen only inside the memory/time-limited replay child.
    let map = read_bytes(&directory.join(&name), 256 << 20)?;
    let stderr_path = directory.join(format!("ell-{ell}.stderr.txt"));
    let stderr = read_bytes(&stderr_path, 8 << 20)?;
    if sha256_hex(&map) != a["receipt"]["stdout_sha256"].as_str().unwrap_or("")
        || sha256_hex(&stderr) != a["receipt"]["stderr_sha256"].as_str().unwrap_or("")
    {
        return Err("construction bytes differ from receipt hashes".into());
    }
    Ok(Submission {
        label: label.into(),
        slug: source["slug"].as_str().unwrap().into(),
        ell,
        summary: summary_bytes,
        map,
        stderr,
    })
}
fn validate_document(s: &Submission) -> Result<(), String> {
    let source = resolve(&s.slug)?;
    let label = s.label.as_str();
    let ell = s.ell;
    let doc: Value =
        serde_json::from_slice(&s.map).map_err(|e| format!("invalid map JSON: {e}"))?;
    if doc["status"] != "PASS"
        || doc["curve"]["preset"] != label
        || doc["result"]["ell"] != ell
        || doc["result"]["source_icv1"]["slug"] != source["slug"]
    {
        return Err("source or degree substitution".into());
    }
    for k in ["p", "a", "b"] {
        if doc["curve"][k].is_null()
            || catalog::integer(&doc["curve"][k]) != catalog::integer(&source["params"][k])
        {
            return Err("source coefficients differ from registered model".into());
        }
    }
    let n = usize::try_from(ell)
        .ok()
        .and_then(|x| x.checked_add(1))
        .ok_or("degree exceeds platform representation")?;
    let maps = doc["result"]["isogenies"]
        .as_array()
        .ok_or("missing explicit maps")?;
    if maps.is_empty() || maps.len() > 2 {
        return Err("expected one or two explicit eigenline maps".into());
    }
    for m in maps {
        for (value, expected) in [
            (&m["kernel"], (n - 2) / 2 + 1),
            (&m["x_map"]["num"], n),
            (&m["x_map"]["den"], n - 1),
        ] {
            if value.as_array().map(Vec::len) != Some(expected) {
                return Err(
                    "explicit coefficient count does not match the claimed prime degree".into(),
                );
            }
        }
    }
    if maps.len() == 2 {
        let p = num_bigint::BigInt::parse_bytes(
            catalog::integer(&source["params"]["p"])
                .to_string()
                .as_bytes(),
            10,
        )
        .unwrap();
        let residue = |v: &Value| {
            let n = num_bigint::BigInt::parse_bytes(catalog::integer(v).to_string().as_bytes(), 10)
                .unwrap();
            ((n % &p) + &p) % &p
        };
        let left = maps[0]["kernel"].as_array().unwrap();
        let right = maps[1]["kernel"].as_array().unwrap();
        if left
            .iter()
            .zip(right)
            .all(|(a, b)| residue(a) == residue(b))
        {
            return Err("duplicate kernels cannot substantiate complete eigenline coverage".into());
        }
        if !two_distinct_eigenlines(source, ell)? {
            return Err("complete two-line coverage requires distinct Frobenius roots".into());
        }
    }
    Ok(())
}
fn two_distinct_eigenlines(source: &Value, ell: u64) -> Result<bool, String> {
    let q = catalog::field_order(source)?;
    let trace = &(&q + &super::Int::one()) - &catalog::integer(&source["order"]);
    // Called after coefficient-count validation: the artifact byte cap already
    // bounds ell far below the u64 modular-product limit of this helper.
    Ok(super::eigenvalues_from_trace(&q, &trace, ell).is_some())
}
fn snapshot(s: &Submission, out: &Path) -> Result<(), String> {
    fs::create_dir(out).map_err(|e| e.to_string())?;
    fs::create_dir(out.join(&s.label)).map_err(|e| e.to_string())?;
    fs::write(out.join("search.json"), &s.summary).map_err(|e| e.to_string())?;
    fs::write(
        out.join(&s.label).join(format!("ell-{}.json", s.ell)),
        &s.map,
    )
    .map_err(|e| e.to_string())?;
    fs::write(
        out.join(&s.label).join(format!("ell-{}.stderr.txt", s.ell)),
        &s.stderr,
    )
    .map_err(|e| e.to_string())
}
fn worker(input: &Path) -> Result<Value, String> {
    let s = preflight(input)?;
    validate_document(&s)?;
    let mut certificate = replay::verify(input, false);
    if certificate["degree_coverage"] == "PARTIAL"
        && !two_distinct_eigenlines(resolve(&s.slug)?, s.ell)?
    {
        certificate["expected_total_eigenlines"] = Value::Null;
        certificate["enumeration_note"] = json!("partial map verified; total eigenline count not established for repeated Frobenius roots");
    }
    let rows = certificate["records"]
        .as_array()
        .ok_or("missing fresh verification records")?;
    if certificate["status"] != "PASS" || rows.is_empty() {
        return Err("fresh independent replay did not pass".into());
    }
    for r in rows {
        if r["status"] != "PASS"
            || r["ell"] != s.ell
            || r["source_icv1"] != s.slug
            || r["public_scalar_transport_checks"] != 20
        {
            return Err(
                "fresh verification does not bind this source, degree and twenty transports".into(),
            );
        }
    }
    Ok(
        json!({"schema":"isogeny-challenge-verified-map/v1","status":"PASS","season":season(),"source_icv1":s.slug,
        "ell":s.ell.to_string(),"map_sha256":sha256_hex(&s.map),"fresh_certificate":certificate}),
    )
}
fn judge(input: &Path, out: &Path, timeout: u64) -> Result<(Value, i32), String> {
    fs::create_dir(out).map_err(|e| format!("judgement needs a fresh directory: {e}"))?;
    let result: Result<(Value, i32), String> = (|| {
        let s = preflight(input)?;
        let artifact = out.join("artifact");
        snapshot(&s, &artifact)?;
        let artifact = artifact.canonicalize().map_err(|e| e.to_string())?;
        let args = vec![
            "challenge".into(),
            "__judge".into(),
            "--input".into(),
            artifact.to_string_lossy().into_owned(),
            "--rss-ready".into(),
            artifact
                .parent()
                .unwrap()
                .join("resources-ready")
                .to_string_lossy()
                .into_owned(),
        ];
        let exe = std::env::current_exe().map_err(|e| e.to_string())?;
        let receipt = supervisor::run(&exe, out, &args, "verification", timeout);
        let receipt: Value = serde_json::from_str(&receipt.dump()).map_err(|e| e.to_string())?;
        if receipt["status"] != "PASS" {
            let invalid = receipt["status"] == "FAIL" && receipt["exit_code"] != 3;
            return Ok((
                json!({"schema":"isogeny-challenge-entry/v1","status":if invalid {"REJECTED"} else {"CENSORED"},"reason":receipt["status"],"score":Value::Null,"verification":receipt}),
                if invalid { 1 } else { 3 },
            ));
        }
        if receipt["sampled_peak_resident_bytes"]
            .as_str()
            .and_then(|n| n.parse::<u64>().ok())
            .unwrap_or(0)
            == 0
        {
            return Ok((
                json!({"schema":"isogeny-challenge-entry/v1","status":"CENSORED",
                "reason":"RESOURCE_MONITOR_UNAVAILABLE","score":Value::Null,"verification":receipt}),
                3,
            ));
        }
        let (_, verified) = read_json(&out.join("verification.json"), 16 << 20)?;
        let expected_ell = s.ell.to_string();
        if verified["status"] != "PASS"
            || verified["season"] != season()
            || verified["ell"] != expected_ell.as_str()
            || verified["source_icv1"] != s.slug
            || verified["map_sha256"] != sha256_hex(&s.map)
        {
            return Err("verification output binding failed".into());
        }
        Ok((
            json!({"schema":"isogeny-challenge-entry/v1","status":"ACCEPTED","season":season(),"source_icv1":s.slug,
            "prime_degree":s.ell.to_string(),"score":s.ell.to_string(),"representation":"expanded rational map", "degree_coverage":verified["fresh_certificate"]["degree_coverage"],
            "map_sha256":sha256_hex(&s.map),"snapshot_search_sha256":sha256_hex(&s.summary),"source_submission":input.to_string_lossy(),
            "judge_binary_sha256":sha256_hex(&fs::read(exe).map_err(|e|e.to_string())?),"tool":provenance(),"verification":receipt,
            "construction_accounting":"external constructor receipts retained; costs are not independently attested by this judge"}),
            0,
        ))
    })();
    let (mut result, code) = match result {
        Ok(r) => r,
        Err(e) => (
            json!({"schema":"isogeny-challenge-entry/v1","status":"REJECTED","reason":e,"score":Value::Null}),
            1,
        ),
    };
    result["season"] = json!(season());
    save(&out.join("entry.json"), &result)?;
    Ok((result, code))
}
pub fn execute(args: &[String]) -> Result<(Value, i32), String> {
    match args.first().map(String::as_str).unwrap_or("policy") {
        "policy" => {
            if args.len() > 1 {
                return Err("policy takes no arguments".into());
            }
            Ok((policy(), 0))
        }
        "plan" => {
            let o = Options::parse(&args[1..], &["curve", "from", "to", "max-order"])?;
            Ok((plan(&o)?, 0))
        }
        "__judge" => {
            let o = Options::parse(&args[1..], &["input", "rss-ready"])?;
            let ready = Path::new(o.need("rss-ready")?);
            let started = std::time::Instant::now();
            while !ready.is_file() {
                if started.elapsed() >= std::time::Duration::from_secs(10) {
                    return Ok((
                        json!({"status":"CENSORED","reason":"RESOURCE_HANDSHAKE_UNAVAILABLE"}),
                        3,
                    ));
                }
                std::thread::sleep(std::time::Duration::from_millis(5));
            }
            let result = std::panic::catch_unwind(|| worker(Path::new(o.need("input").unwrap())));
            // All parsing and algebra have completed, including failure paths.
            // The parent can stop RSS polling during process-exit cleanup.
            fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(ready.with_extension("done"))
                .map_err(|e| e.to_string())?;
            match result {
                Ok(result) => Ok((result?, 0)),
                Err(_) => Err("independent replay rejected the submitted map".into()),
            }
        }
        "judge" => {
            let o = Options::parse(&args[1..], &["input", "out", "timeout"])?;
            let timeout = o.number("timeout", Some(3600))?;
            if timeout == 0 || timeout > 3600 {
                return Err("judge timeout must be 1..3600 seconds".into());
            }
            judge(
                Path::new(o.need("input")?),
                Path::new(o.need("out")?),
                timeout,
            )
        }
        "run" => {
            let o = Options::parse(&args[1..], &["curve", "ell", "out", "timeout"])?;
            let row = resolve(o.need("curve")?)?;
            let ell = o.number("ell", None)?;
            let timeout = o.number("timeout", Some(1800))?;
            if !(3..=1_000_000).contains(&ell) || !is_prime(ell) {
                return Err("bundled construction needs an odd prime degree <= 1000000".into());
            }
            if timeout == 0 || timeout > 1800 {
                return Err("construction timeout must be 1..1800 seconds".into());
            }
            if catalog::capabilities(row)["independent_replay"] != true {
                return Err("curve is listed but the full-map adapter is unavailable".into());
            }
            let name = catalog::label(row);
            let rep = catalog::representation(row).ok_or("missing registered subgroup")?;
            if catalog::integer(&rep["curve"]["subgroup_order"]).mod_u64(ell) == 0 {
                return Err("degree must be coprime to the registered subgroup order".into());
            }
            if super::eigenvalues(name, ell).is_none() {
                return Err(
                    "bundled construction needs a split prime with distinct Frobenius eigenvalues"
                        .into(),
                );
            }
            let out = Path::new(o.need("out")?);
            fs::create_dir(out).map_err(|e| e.to_string())?;
            let artifact = out.join("construction");
            let command = vec![
                "search".into(),
                "--curve".into(),
                catalog::label(row).into(),
                "--ell".into(),
                ell.to_string(),
                "--construct-only".into(),
                "--timeout".into(),
                timeout.to_string(),
                "--out".into(),
                artifact.to_string_lossy().into_owned(),
            ];
            let (construction, code) = super::execute(&command)?;
            save(&out.join("construction-result.json"), &construction)?;
            let (result, code) = if code == 0 {
                judge(&artifact, &out.join("judgement"), 3600)?
            } else {
                (
                    json!({"status":"CENSORED","score":Value::Null,"construction":construction}),
                    code,
                )
            };
            save(&out.join("RESULT.json"), &result)?;
            Ok((result, code))
        }
        "leaderboard" => {
            let o = Options::parse(&args[1..], &["input", "out"])?;
            let list = Path::new(o.need("input")?);
            let (_, manifest) = read_json(list, 1 << 20)?;
            let inputs = manifest
                .as_array()
                .ok_or("input must be a JSON array of raw submission directory paths")?;
            if inputs.len() > 256 {
                return Err("at most 256 submissions per judgement batch".into());
            }
            if inputs.iter().any(|p| !p.is_string()) {
                return Err("each submission must be a path string, never a saved score".into());
            }
            let out = Path::new(o.need("out")?);
            fs::create_dir(out).map_err(|e| e.to_string())?;
            let mut attempts = vec![];
            let mut best: BTreeMap<String, (u64, Vec<Value>)> = BTreeMap::new();
            for (i, item) in inputs.iter().enumerate() {
                let name = item
                    .as_str()
                    .ok_or("each submission must be a path string, never a saved score")?;
                let input = PathBuf::from(name);
                let input = if input.is_absolute() {
                    input
                } else {
                    list.parent().unwrap_or(Path::new(".")).join(input)
                };
                let (entry, code) = judge(&input, &out.join(format!("entry-{i}")), 3600)?;
                if code == 0 {
                    let slug = entry["source_icv1"].as_str().unwrap().to_owned();
                    let ell = entry["prime_degree"]
                        .as_str()
                        .unwrap()
                        .parse::<u64>()
                        .unwrap();
                    let row = best.entry(slug).or_insert((ell, vec![]));
                    if ell > row.0 {
                        *row = (ell, vec![]);
                    }
                    if ell == row.0 {
                        row.1.push(entry.clone());
                    }
                }
                attempts.push(entry);
            }
            let mut rows = curves();
            for row in &mut rows {
                if let Some((ell, winners)) = best.get(row["source_icv1"].as_str().unwrap()) {
                    row["best_prime_degree"] = json!(ell.to_string());
                    row["status"] = json!("VERIFIED");
                    row["co_leaders"] = json!(winners);
                }
            }
            let result = json!({"schema":"isogeny-challenge-leaderboard/v1","status":"PASS","season":season(),"curves":rows,"attempts":attempts,"ranking":"separate per exact curve; equal degrees share rank; every entry replayed in this batch"});
            save(&out.join("leaderboard.json"), &result)?;
            Ok((result, 0))
        }
        _ => Err("challenge commands: policy, plan, run, judge, leaderboard".into()),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn catalogue_has_every_standard_prime_model_and_no_research_descendants() {
        let rows = curves();
        let canonical: Value = serde_json::from_str(crate::REGISTRY).unwrap();
        let expected: std::collections::BTreeSet<_> = canonical["curves"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|r| {
                r["family"] == "prime"
                    && r["standard_names"]
                        .as_array()
                        .is_some_and(|names| !names.is_empty())
            })
            .map(|r| r["slug"].as_str().unwrap())
            .collect();
        let actual: std::collections::BTreeSet<_> = rows
            .iter()
            .map(|r| r["source_icv1"].as_str().unwrap())
            .collect();
        assert_eq!(actual, expected);
        assert!(rows.len() > 20);
        for r in rows {
            assert_eq!(
                resolve(r["source_icv1"].as_str().unwrap()).unwrap()["family"],
                "prime"
            );
            assert!(!r["standard_names"].as_array().unwrap().is_empty());
        }
        let research = catalog::registry()["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["family"] == "prime" && !standard_prime(r))
            .unwrap();
        assert!(resolve(research["slug"].as_str().unwrap()).is_err());
    }
    #[test]
    fn policy_identifies_rules_and_catalogue() {
        let p = policy();
        assert_eq!(p["rules_sha256"], sha256_hex(RULES));
        assert_eq!(p["season"].as_str().unwrap().len(), 64);
        assert_eq!(p["rules"]["replay_timeout_seconds"], 3600);
    }
    #[test]
    fn plan_rejects_unbounded_ranges_and_nonprime_models() {
        assert!(execute(&["plan".into(), "--to".into(), "1000001".into()]).is_err());
        assert!(resolve("sect163k1").is_err());
    }
}
