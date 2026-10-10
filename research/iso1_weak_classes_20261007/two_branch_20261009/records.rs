//! Canonical identities for independently replayed source, route and literal models.
use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, path::Path};
fn hash(s: &[u8]) -> String {
    hex::encode(sha256(s))
}
fn kv(s: &str) -> BTreeMap<&str, &str> {
    s.split('|')
        .skip(1)
        .filter_map(|v| v.split_once('='))
        .collect()
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 2, "STUDY REPO");
    let base = Path::new(&a[0]);
    let repo = Path::new(&a[1]);
    let mut raw = fs::read_to_string(base.join("replay_run1/stdout.txt")).unwrap();
    assert!(fs::read(base.join("replay_run1/stderr.txt"))
        .unwrap()
        .is_empty());
    assert!(raw.lines().last().unwrap().starts_with("ALGEBRA_COMPLETE|"));
    let general = base.join("general_replay_run2/stdout.txt");
    if general.exists() {
        assert!(fs::read(base.join("general_replay_run2/stderr.txt"))
            .unwrap()
            .is_empty());
        let text = fs::read_to_string(general).unwrap();
        assert!(text
            .lines()
            .last()
            .unwrap()
            .starts_with("ALGEBRA_COMPLETE|"));
        raw.push_str(&text);
    }
    let mut records = Vec::new();
    let mut entries = Vec::new();
    for line in raw.lines().filter(|l| l.starts_with("MODEL|")) {
        let m = kv(line);
        let p = m["p"];
        let degree = m.get("degree").copied().unwrap_or("6");
        let degree_number: usize = degree.parse().unwrap();
        let modulus: Value = serde_json::from_str(m["modulus"]).unwrap();
        let av: Value = serde_json::from_str(m["short_a"]).unwrap();
        let bv: Value = serde_json::from_str(m["short_b"]).unwrap();
        let j: Value = serde_json::from_str(m["j"]).unwrap();
        let mods = modulus
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect::<Vec<_>>()
            .join(",");
        let field = format!(
            "fpk-{p}-{degree}-{}",
            &hash(format!("fpk-modulus:{p}:{mods}").as_bytes())[..8]
        );
        let model = json!({"a":av,"b":bv,"field":field,"form":"y^2=x^3+a*x+b","k":degree,"modulus":modulus,"p":p,"v":"1"});
        let canonical = serde_json::to_string(&model).unwrap();
        let h = hash(canonical.as_bytes());
        let trace = m["trace"];
        let order = m["order"];
        let signed = if let Some(t) = trace.strip_prefix('-') {
            format!("tm{t}")
        } else {
            format!("t{trace}")
        };
        let jj = j
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect::<Vec<_>>()
            .join(",");
        let prime: u64 = p.parse().unwrap();
        let bits = 64 - prime.leading_zeros();
        let slug = format!("icv1-fp{bits}k{degree}-{signed}-{}", &h[..8]);
        let icv1 = format!("ICV1:{field}:{trace}:{order}:{jj}:unk:unk:r:{}", &h[..12]);
        records.push(json!({"slug":slug,"icv1":icv1,"model_json":canonical,"model_sha256":h,"trace":trace,"order":order,"j_coefficients":j,"endomorphism_order_conductor":null,"volcano_level":null,"subgroup_order":null,"generator":null,"verification":"independent recorded-model point count, saved-route kernel replay and point transport","raw_source":if degree_number==6{"replay_run1/stdout.txt"}else{"general_replay_run2/stdout.txt"}}));
        entries.push(json!({"slug":slug,"icv1":icv1,"family":"extension","params":{"p":p,"k":degree_number,"modulus":modulus,"a":av,"b":bv,"generator_call":{"law":"frozen independent source, kernel-certified route, or literal weak-model conversion; higher-degree rows are separate constructed controls","source":"recorded coordinate replay"}},"model_json":canonical,"trace":trace,"order":order,"j":jj,"end":"unk","aliases":[],"standard_names":[],"sources":["research/iso1_weak_classes_20261007/two_branch_20261009/curve_records.json"],"representations":[]}));
    }
    fs::write(base.join("curve_records.json"),serde_json::to_vec_pretty(&json!({"schema":"iso1.two_branch_curve_records/v1","replay_sha256":hash(raw.as_bytes()),"replay_hash_scope":"concatenation of replay_run1 and, when present, general_replay_run2 stdout bytes","records":records})).unwrap()).unwrap();
    let path = repo.join("docs/curves/registry.json");
    let old = fs::read_to_string(&path).unwrap();
    let marker = "\n ],\n \"standards_source\":";
    let at = old.rfind(marker).unwrap();
    let mut appendix = String::new();
    let mut added = 0;
    for e in &entries {
        if !old.contains(&format!("\"slug\": \"{}\"", e["slug"].as_str().unwrap())) {
            appendix.push_str(",\n");
            appendix.push_str(&serde_json::to_string_pretty(e).unwrap());
            added += 1;
        }
    }
    fs::write(path, format!("{}{}{}", &old[..at], appendix, &old[at..])).unwrap();
    fs::write(repo.join("docs/curves/sources/iso1_two_branch_20261009.json"),serde_json::to_vec_pretty(&json!({"schema":"iso1-curve-sources/v1","generated_by":"iso1_two_branch_records","curves":entries})).unwrap()).unwrap();
    let specs = repo.join("docs/curves/sources/specs.txt");
    let mut s = fs::read_to_string(&specs).unwrap();
    let entry = "# Native seed/source manifest: docs/curves/sources/iso1_two_branch_20261009.json";
    if !s.lines().any(|v| v == entry) {
        s.push_str(&format!("\n{entry}\n"));
        fs::write(specs, s).unwrap();
    }
    println!(
        "{} independently replayed models; {added} registered",
        records.len()
    );
}
