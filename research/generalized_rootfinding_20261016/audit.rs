//! Native integrity and recorded-output audit; never launches historical Python.
use flate2::read::GzDecoder;
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, io::Read, path::{Component, Path}};

#[allow(dead_code)]
#[path = "../../src/hash/sha256.rs"]
mod sha256;

fn digest(bytes: &[u8]) -> String {
    sha256::sha256(bytes).iter().map(|b| format!("{b:02x}")).collect()
}

fn safe_relative(name: &str) -> Result<&Path, String> {
    let path = Path::new(name);
    if name.is_empty() || !path.components().all(|c| matches!(c, Component::Normal(_))) {
        return Err(format!("unsafe relative path: {name}"));
    }
    Ok(path)
}

fn text<'a>(v: &'a Value, key: &str) -> Result<&'a str, String> {
    v[key].as_str().ok_or_else(|| format!("missing string {key}"))
}

fn decode(root: &Path, entry: &Value) -> Result<Vec<u8>, String> {
    let stored = safe_relative(text(entry, "stored")?)?;
    let encoded = fs::read(root.join(stored)).map_err(|e| e.to_string())?;
    if digest(&encoded) != text(entry, "transport_sha256")? {
        return Err(format!("transport digest mismatch: {}", stored.display()));
    }
    let bytes = match text(entry, "encoding")? {
        "identity" => encoded,
        "gzip" => {
            let mut out = Vec::new();
            GzDecoder::new(encoded.as_slice()).take(512 * 1024 * 1024 + 1)
                .read_to_end(&mut out).map_err(|e| e.to_string())?;
            out
        }
        other => return Err(format!("unsupported encoding: {other}")),
    };
    if entry["bytes"].as_u64() != Some(bytes.len() as u64)
        || digest(&bytes) != text(entry, "sha256")? {
        return Err(format!("original length/digest mismatch: {}", stored.display()));
    }
    Ok(bytes)
}

fn audit(root: &Path) -> Result<Value, String> {
    let manifest: Value = serde_json::from_slice(&fs::read(root.join("TRANSPORT.json"))
        .map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    let entries = manifest["entries"].as_array().ok_or("missing entries")?;
    let mut expected = BTreeMap::new();
    let mut results = BTreeMap::new();
    let mut source_manifests = Vec::new();
    for entry in entries {
        let name = text(entry, "original")?;
        safe_relative(name)?;
        if expected.insert(name.to_owned(), text(entry, "sha256")?.to_owned()).is_some() {
            return Err(format!("duplicate original path: {name}"));
        }
        let bytes = decode(root, entry)?;
        if name.ends_with("/SOURCE_MANIFEST.json") {
            source_manifests.push(serde_json::from_slice::<Value>(&bytes).map_err(|e| e.to_string())?);
        }
        if name.ends_with("/scientific_result.json") {
            let result: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
            if result["schema"] != "scientific_result_v1" || result["status"] != "passed"
                || result["evidence"]["release_gate_claim"] != false {
                return Err(format!("unexpected result contract: {name}"));
            }
            results.insert(text(&result, "experiment_id")?.to_owned(), result);
        }
    }
    if expected.len() != 112 || results.len() != 3 {
        return Err("expected frozen 112-file snapshot and three results".into());
    }
    for source in source_manifests {
        let pins = source["source_and_input_sha256"].as_object().ok_or("missing source pins")?;
        for (name, hash) in pins {
            if expected.get(name).map(String::as_str) != hash.as_str() {
                return Err(format!("source closure mismatch: {name}"));
            }
        }
    }
    let batch = results.get("batch-separator-manyfield-exact-coldpairs-20261016-v1")
        .ok_or("missing batch result")?;
    let evidence = &batch["evidence"];
    if evidence["distinct_fields"] != 20 || evidence["curve_contexts"] != 38 {
        return Err("field/context coverage mismatch".into());
    }
    let records = evidence["records"].as_array().ok_or("missing cold records")?;
    if records.len() != 6 { return Err("expected three cold pairs".into()); }
    let mut reference: BTreeMap<String, Value> = BTreeMap::new();
    let mut output_comparisons = 0;
    for record in records {
        let rows = record["result"]["rows"].as_array().ok_or("missing context rows")?;
        if rows.len() != 38 { return Err("incomplete cold context panel".into()); }
        for row in rows {
            let case = text(row, "case")?.to_owned();
            let answers = &row["answers"];
            if !answers.is_array() { return Err("missing recorded all-root answers".into()); }
            if let Some(previous) = reference.get(&case) {
                if previous != answers { return Err(format!("recorded output mismatch: {case}")); }
                output_comparisons += 1;
            } else { reference.insert(case, answers.clone()); }
        }
    }
    if reference.len() != 38 { return Err("duplicate/missing contexts".into()); }
    Ok(json!({"schema":"historical-rootfinding-audit/v1", "status":"passed",
        "payloads_verified":expected.len(), "scientific_results":results.len(),
        "recorded_output_comparisons":output_comparisons, "independent_resolve":false,
        "qualified_performance_claim":false, "python_executed":false}))
}

fn main() {
    let root = std::env::args().nth(1).unwrap_or_else(|| env!("CARGO_MANIFEST_DIR").into());
    match audit(Path::new(&root)) {
        Ok(value) => println!("{value}"),
        Err(error) => { eprintln!("audit failed: {error}"); std::process::exit(1); }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_escape_paths() {
        for name in ["", "/tmp/x", "../x", "a/../x", "./x"] {
            assert!(safe_relative(name).is_err(), "{name}");
        }
        assert!(safe_relative("historical/valid.json").is_ok());
    }
    #[test]
    fn sha256_known_answer() {
        assert_eq!(digest(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    }
    #[test]
    fn rejects_modified_transport() {
        let root = std::env::temp_dir().join(format!("rootfinding-audit-{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        fs::write(root.join("changed"), b"changed").unwrap();
        let entry = json!({"stored":"changed", "transport_sha256":digest(b"original"),
            "encoding":"identity", "bytes":8, "sha256":digest(b"original")});
        assert!(decode(&root, &entry).unwrap_err().contains("transport digest mismatch"));
        fs::remove_file(root.join("changed")).unwrap();
        fs::remove_dir(root).unwrap();
    }
}
