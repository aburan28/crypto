//! Research companion for a future frozen campaign. No target is embedded.
//! Parsing and arithmetic are confined to an RSS-acknowledged worker.
#![allow(dead_code)]
#[path = "../verification_field.rs"]
mod verification_field;
mod field {
    pub type Field = crate::verification_field::PrimeField<3>;
    pub type Fe = crate::verification_field::Element<3>;
    pub const FIELD_BITS: usize = 192;
}
#[path = "../../../../src/cryptanalysis/isogeny_walk/curve.rs"]
mod curve;
#[path = "../../../../research/p192_p224_large_degree_isogenies_20261009/verification_kernel.rs"]
mod kernel;
#[path = "../verification_modpoly.rs"]
mod modpoly;
#[path = "../../../../research/p192_p224_large_degree_isogenies_20261009/verification_poly.rs"]
mod poly;
#[path = "../supervisor.rs"]
mod supervisor;

use field::{Fe, Field};
use num_bigint::BigUint;
use poly::{FixedModulus, Poly};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

const P192: &str = "6277101735386680763835789423207666416083908700390324961279";
const SEED: &str = "isogeny-fixed-modulus-v1";
const POLICY: &str = if cfg!(feature = "cached-fixed-modulus-reduction") {
    "cached-fixed-modulus-v1"
} else {
    "reference-rem-v1"
};
fn sha(bytes: &[u8]) -> String {
    hex::encode(Sha256::digest(bytes))
}
fn fresh(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| e.to_string())?;
    file.write_all(bytes)
        .and_then(|_| file.sync_all())
        .map_err(|e| e.to_string())
}
fn write_json(path: &Path, value: &impl Serialize) -> Result<(), String> {
    let mut bytes = serde_json::to_vec_pretty(value).map_err(|e| e.to_string())?;
    bytes.push(b'\n');
    fresh(path, &bytes)
}
fn encode(f: &Field, p: &Poly) -> Vec<String> {
    p.iter().map(|c| f.to_big(c).to_string()).collect()
}
fn decode(f: &Field, p: &[String]) -> Result<Poly, String> {
    p.iter()
        .map(|s| {
            if s.is_empty()
                || (s.len() > 1 && s.starts_with('0'))
                || !s.bytes().all(|b| b.is_ascii_digit())
            {
                return Err("noncanonical decimal coefficient".into());
            }
            let x = BigUint::parse_bytes(s.as_bytes(), 10).ok_or("invalid coefficient")?;
            if x >= *f.modulus() {
                return Err("coefficient outside field".into());
            }
            Ok(f.from_big(&x))
        })
        .collect()
}
fn output_digest(f: &Field, results: &[Poly]) -> String {
    sha(&serde_json::to_vec(&results.iter().map(|p| encode(f, p)).collect::<Vec<_>>()).unwrap())
}
#[derive(Serialize, Deserialize, Clone, Debug)]
#[serde(deny_unknown_fields)]
struct Case {
    id: String,
    kind: String,
    modulus: Vec<String>,
    pairs: Vec<[Vec<String>; 2]>,
    outer: Vec<String>,
    inner: Vec<String>,
    expected_sha256: String,
}
#[derive(Serialize, Deserialize, Debug)]
#[serde(deny_unknown_fields)]
struct Dataset {
    schema: String,
    field: String,
    seed: String,
    input_map_sha256: String,
    extracted_kernel_sha256: String,
    extracted_kernel: Vec<String>,
    cases: Vec<Case>,
}
fn stream(f: &Field, case: &str, role: &str, ordinal: usize, length: usize) -> Poly {
    (0..length)
        .map(|index| {
            // Length-prefixed components make the domains unambiguous.
            let mut h = Sha256::new();
            for part in [SEED.as_bytes(), case.as_bytes(), role.as_bytes()] {
                h.update((part.len() as u64).to_le_bytes());
                h.update(part);
            }
            h.update((ordinal as u64).to_le_bytes());
            h.update((index as u64).to_le_bytes());
            f.from_big(&(BigUint::from_bytes_be(&h.finalize()) % f.modulus()))
        })
        .collect()
}
fn make_case(f: &Field, id: &str, kind: &str, m: Poly, count: usize) -> Case {
    let d = m.len() - 1;
    let pairs = if kind == "product" {
        (0..count)
            .map(|i| {
                [
                    encode(f, &stream(f, id, "left", i, d)),
                    encode(f, &stream(f, id, "right", i, d)),
                ]
            })
            .collect()
    } else {
        vec![]
    };
    let (outer, inner) = if kind == "composition" {
        (
            encode(f, &stream(f, id, "outer", 0, d + 1)),
            encode(f, &stream(f, id, "inner", 0, d)),
        )
    } else {
        (vec![], vec![])
    };
    let mut c = Case {
        id: id.into(),
        kind: kind.into(),
        modulus: encode(f, &m),
        pairs,
        outer,
        inner,
        expected_sha256: String::new(),
    };
    c.expected_sha256 = output_digest(
        f,
        &evaluate_reference(f, &c).expect("generated reference case"),
    );
    c
}
fn evaluate_reference(f: &Field, c: &Case) -> Result<Vec<Poly>, String> {
    let m = decode(f, &c.modulus)?;
    FixedModulus::new(f, &m)?; // Validate; this setup is not an operation timer.
    if c.kind == "product" {
        c.pairs
            .iter()
            .map(|[a, b]| Ok(poly::mulmod(f, &decode(f, a)?, &decode(f, b)?, &m)))
            .collect()
    } else if c.kind == "composition" {
        Ok(vec![kernel::compose_mod(
            f,
            &decode(f, &c.outer)?,
            &decode(f, &c.inner)?,
            &m,
        )])
    } else {
        Err("unknown case kind".into())
    }
}
fn generate(input: &Path, out: &Path) -> Result<Value, String> {
    let start = Instant::now();
    let raw = fs::read(input).map_err(|e| e.to_string())?;
    let document: Value = serde_json::from_slice(&raw).map_err(|e| e.to_string())?;
    if document["curve"]["p"].as_str() != Some(P192) {
        return Err("P-192 map required".into());
    }
    let f = Field::new(&BigUint::parse_bytes(P192.as_bytes(), 10).unwrap()).unwrap();
    let kernel: Vec<String> =
        serde_json::from_value(document["result"]["isogenies"][0]["kernel"].clone())
            .map_err(|e| e.to_string())?;
    let h = decode(&f, &kernel)?;
    if h.len() < 2 || h.last() != Some(&f.one()) {
        return Err("nonconstant monic kernel required".into());
    }
    let parsing_ms = start.elapsed().as_millis();
    let mut cases = vec![];
    for (degree, count) in [(127, 256), (1023, 32)] {
        let id = format!("synthetic-{degree}");
        let mut m = stream(&f, &id, "modulus", 0, degree + 1);
        m[degree] = f.one();
        cases.push(make_case(
            &f,
            &format!("{id}-products"),
            "product",
            m.clone(),
            count,
        ));
        cases.push(make_case(
            &f,
            &format!("{id}-cold"),
            "product",
            m.clone(),
            1,
        ));
        cases.push(make_case(
            &f,
            &format!("{id}-composition"),
            "composition",
            m,
            1,
        ));
    }
    cases.push(make_case(
        &f,
        "actual-kernel-products",
        "product",
        h.clone(),
        4,
    ));
    cases.push(make_case(&f, "actual-kernel-cold", "product", h, 1));
    let dataset = Dataset {
        schema: "isogeny-fixed-modulus-dataset/v1".into(),
        field: P192.into(),
        seed: SEED.into(),
        input_map_sha256: sha(&raw),
        extracted_kernel_sha256: sha(&serde_json::to_vec(&kernel).unwrap()),
        extracted_kernel: kernel,
        cases,
    };
    let generation_ms = start.elapsed().as_millis() - parsing_ms;
    write_json(&out.join("dataset.json"), &dataset)?;
    Ok(
        json!({"schema":"fixed-modulus-generation/v1","status":"PASS","compiled_policy":POLICY,
        "input_map_sha256":dataset.input_map_sha256,"extracted_kernel_sha256":dataset.extracted_kernel_sha256,
        "dataset_sha256":sha(&fs::read(out.join("dataset.json")).map_err(|e|e.to_string())?),
        "parsing_ms_operational_only":parsing_ms,"generation_and_reference_verification_ms_operational_only":generation_ms,
        "total_ms_operational_only":start.elapsed().as_millis(),"case_count":dataset.cases.len(),
        "cost_scope":"generation includes reference expected-output calculation; no measured comparison"}),
    )
}
fn measure(input: &Path, id: &str, out: &Path) -> Result<Value, String> {
    let total = Instant::now();
    let bytes = fs::read(input).map_err(|e| e.to_string())?;
    let dataset: Dataset = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    if dataset.schema != "isogeny-fixed-modulus-dataset/v1"
        || dataset.field != P192
        || dataset.seed != SEED
    {
        return Err("unsupported dataset contract".into());
    }
    let f = Field::new(&BigUint::parse_bytes(P192.as_bytes(), 10).unwrap()).unwrap();
    let matches: Vec<_> = dataset.cases.iter().filter(|c| c.id == id).collect();
    if matches.len() != 1 {
        return Err("case must exist exactly once".into());
    }
    let c = matches[0];
    let m = decode(&f, &c.modulus)?;
    let pairs: Vec<_> = c
        .pairs
        .iter()
        .map(|[a, b]| Ok((decode(&f, a)?, decode(&f, b)?)))
        .collect::<Result<_, String>>()?;
    let outer = decode(&f, &c.outer)?;
    let inner = decode(&f, &c.inner)?;
    if !matches!(c.kind.as_str(), "product" | "composition")
        || (c.kind == "product" && pairs.is_empty())
    {
        return Err("invalid operation case".into());
    }
    if m.len() < 2 || m.last() != Some(&f.one()) {
        return Err("nonconstant monic modulus required".into());
    }
    let parsing = total.elapsed();
    let setup_start = Instant::now();
    #[cfg(feature = "cached-fixed-modulus-reduction")]
    let context = FixedModulus::new(&f, &m)?;
    let setup = setup_start.elapsed();
    let op_start = Instant::now();
    let results: Vec<Poly> = if c.kind == "product" {
        pairs
            .iter()
            .map(|(a, b)| {
                #[cfg(feature = "cached-fixed-modulus-reduction")]
                {
                    context.mulmod(std::hint::black_box(a), std::hint::black_box(b))
                }
                #[cfg(not(feature = "cached-fixed-modulus-reduction"))]
                {
                    poly::mulmod(&f, std::hint::black_box(a), std::hint::black_box(b), &m)
                }
            })
            .collect()
    } else {
        #[cfg(feature = "cached-fixed-modulus-reduction")]
        {
            vec![kernel::compose_mod_cached(
                std::hint::black_box(&outer),
                std::hint::black_box(&inner),
                &context,
            )]
        }
        #[cfg(not(feature = "cached-fixed-modulus-reduction"))]
        {
            vec![kernel::compose_mod(
                &f,
                std::hint::black_box(&outer),
                std::hint::black_box(&inner),
                &m,
            )]
        }
    };
    std::hint::black_box(&results);
    let operation = op_start.elapsed();
    let encoding_start = Instant::now();
    let digest = output_digest(&f, &results);
    let encoding = encoding_start.elapsed();
    let verify_start = Instant::now();
    if digest != c.expected_sha256 {
        return Err("reference output digest mismatch".into());
    }
    let verification = verify_start.elapsed();
    let record = json!({"schema":"fixed-modulus-measurement/v1","status":"PASS","case":id,
        "compiled_policy":POLICY,"dataset_sha256":sha(&bytes),"output_sha256":digest,
        "operation_count":results.len(),"modulus_degree":m.len()-1,
        "parsing_ns":parsing.as_nanos(),"context_setup_ns":setup.as_nanos(),"operation_ns":operation.as_nanos(),
        "encoding_checksum_ns":encoding.as_nanos(),"verification_ns":verification.as_nanos(),
        "total_ns":total.elapsed().as_nanos(),"cold_operation_ns":(setup+operation).as_nanos(),
        "cpu_time":"requires external OS user/system receipt; stage CPU is not estimated",
        "git_commit":env!("ISOGENY_GIT_COMMIT"),"git_dirty":env!("ISOGENY_GIT_DIRTY")});
    write_json(&out.join("measurement.json"), &record)?;
    Ok(record)
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct AnalysisManifest {
    schema: String,
    cases: Vec<AnalysisCase>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct AnalysisCase {
    id: String,
    metric: String,
    // Ordered paths, not precomputed times. Future frozen manifests must bind
    // their hashes. Every receipt is independently checked below.
    aa: Vec<[BoundMeasurement; 2]>,
    abba: Vec<[BoundMeasurement; 4]>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BoundMeasurement {
    path: PathBuf,
    sha256: String,
}
fn observed(
    r: &BoundMeasurement,
    id: &str,
    metric: &str,
    policy: &str,
) -> Result<(f64, String, String), String> {
    if !r.path.is_absolute() {
        return Err("measurement path must be absolute".into());
    }
    let bytes = fs::read(&r.path).map_err(|e| e.to_string())?;
    if sha(&bytes) != r.sha256 {
        return Err("measurement hash mismatch".into());
    }
    let v: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    if v["schema"] != "fixed-modulus-measurement/v1"
        || v["status"] != "PASS"
        || v["case"] != id
        || v["compiled_policy"] != policy
    {
        return Err("measurement identity/status/policy mismatch".into());
    }
    if !matches!(metric, "total_ns" | "cold_operation_ns" | "operation_ns") {
        return Err("unsupported wall metric; CPU requires a separately bound OS receipt".into());
    }
    let n = v[metric].as_u64().ok_or("missing nonnegative duration")?;
    if n == 0 {
        return Err("zero duration cannot be adjudicated".into());
    }
    Ok((
        n as f64,
        v["dataset_sha256"]
            .as_str()
            .ok_or("missing dataset hash")?
            .into(),
        v["output_sha256"]
            .as_str()
            .ok_or("missing output hash")?
            .into(),
    ))
}
fn adjudicate(aa: &[[f64; 2]], blocks: &[[f64; 4]]) -> Result<Value, String> {
    if aa.len() != 8 || blocks.len() != 8 {
        return Err("requires all eight A/A pairs and eight ABBA blocks".into());
    }
    if aa
        .iter()
        .flatten()
        .chain(blocks.iter().flatten())
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err("all samples must be finite and positive".into());
    }
    let mut noise: Vec<f64> = aa
        .iter()
        .map(|[a, b]| 2.0 * (a - b).abs() / (a + b))
        .collect();
    let original_noise = noise.clone();
    noise.sort_by(f64::total_cmp);
    let rank = (9 * noise.len()).div_ceil(10); // one based nearest rank: 8 of 8
    let p90 = noise[rank - 1];
    let ratios: Vec<f64> = blocks
        .iter()
        .map(|[a, b, c, d]| ((b / a) * (c / d)).sqrt())
        .collect();
    let mut sorted = ratios.clone();
    sorted.sort_by(f64::total_cmp);
    let median = (sorted[3] + sorted[4]) / 2.0;
    Ok(
        json!({"aa_relative_differences":original_noise,"aa_nearest_rank_p90":p90,"aa_rank":rank,
        "abba_block_ratios":ratios,"median_paired_b_over_a":median,
        "scoped_microcase_gate_passed":p90<=0.10 && ratios.iter().all(|r|*r<1.0) && median<=0.90,
        "wall_metric_scope":"descriptive unless host isolation/thermal/background gates are separately satisfied",
        "end_to_end_replay_claim":false}),
    )
}
fn analyze(input: &Path, out: &Path) -> Result<Value, String> {
    let start = Instant::now();
    let bytes = fs::read(input).map_err(|e| e.to_string())?;
    let manifest: AnalysisManifest = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    if manifest.schema != "fixed-modulus-analysis-manifest/v1" || manifest.cases.is_empty() {
        return Err("invalid analysis manifest".into());
    }
    let mut results = vec![];
    let mut used = std::collections::BTreeSet::new();
    for case in manifest.cases {
        let mut identity = None;
        let mut load = |b: &BoundMeasurement, policy: &str| -> Result<f64, String> {
            if !used.insert(b.path.clone()) {
                return Err("measurement reused in manifest".into());
            }
            let (t, digest, output) = observed(b, &case.id, &case.metric, policy)?;
            match &identity {
                None => identity = Some((digest, output)),
                Some(old) if *old == (digest, output) => {}
                _ => return Err("dataset/output differs between arms".into()),
            }
            Ok(t)
        };
        let mut aa = vec![];
        let mut abba = vec![];
        for [a, b] in &case.aa {
            aa.push([load(a, "reference-rem-v1")?, load(b, "reference-rem-v1")?]);
        }
        for [a, b, c, d] in &case.abba {
            abba.push([
                load(a, "reference-rem-v1")?,
                load(b, "cached-fixed-modulus-v1")?,
                load(c, "cached-fixed-modulus-v1")?,
                load(d, "reference-rem-v1")?,
            ]);
        }
        results.push(
            json!({"case":case.id,"metric":case.metric,"adjudication":adjudicate(&aa,&abba)?}),
        );
    }
    let record = json!({"schema":"fixed-modulus-analysis/v1","status":"PASS","manifest_sha256":sha(&bytes),
        "cases":results,"total_ms_operational_only":start.elapsed().as_millis(),
        "limits":"CPU adjudication and external host/resource gates remain independently required; this is not public admission"});
    write_json(&out.join("analysis.json"), &record)?;
    Ok(record)
}
fn options(
    args: &[String],
    allowed: &[&str],
) -> Result<std::collections::BTreeMap<String, String>, String> {
    if !args.len().is_multiple_of(2) {
        return Err("options need values".into());
    }
    let mut o = std::collections::BTreeMap::new();
    for pair in args.chunks_exact(2) {
        let key = pair[0].strip_prefix("--").ok_or("option required")?;
        if !allowed.contains(&key) || o.insert(key.into(), pair[1].clone()).is_some() {
            return Err("unknown or duplicate option".into());
        }
    }
    Ok(o)
}
fn unprotected(args: &[String]) -> Result<Value, String> {
    let operation = args.first().ok_or("operation required")?;
    let allowed = if operation == "measure" {
        vec!["input", "out", "case"]
    } else {
        vec!["input", "out"]
    };
    let o = options(&args[1..], &allowed)?;
    let input = Path::new(o.get("input").ok_or("--input required")?);
    let out = Path::new(o.get("out").ok_or("--out required")?);
    match operation.as_str() {
        "generate" => generate(input, out),
        "measure" => measure(input, o.get("case").ok_or("--case required")?, out),
        "analyze" => analyze(input, out),
        _ => Err("unknown operation".into()),
    }
}
fn execute(args: &[String]) -> Result<(Value, i32), String> {
    if args.first().is_some_and(|s| s == "fixed-modulus-bench")
        && args.get(1).is_some_and(|s| s == "__worker")
    {
        if args.len() < 5 || args[args.len() - 2] != "--rss-ready" {
            return Err("positive RSS acknowledgement required".into());
        }
        let ready = Path::new(&args[args.len() - 1]);
        if !ready.is_absolute() {
            return Err("absolute RSS acknowledgement path required".into());
        }
        let start = Instant::now();
        while !ready.is_file() {
            if start.elapsed() >= Duration::from_secs(10) {
                return Ok((
                    json!({"status":"CENSORED","reason":"RESOURCE_HANDSHAKE_UNAVAILABLE"}),
                    3,
                ));
            }
            std::thread::sleep(Duration::from_millis(5));
        }
        let result = std::panic::catch_unwind(|| unprotected(&args[2..args.len() - 2]));
        fresh(&ready.with_extension("done"), b"")?;
        return match result {
            Ok(Ok(v)) => Ok((json!({"status":"PASS","record":v}), 0)),
            Ok(Err(e)) => Ok((json!({"status":"FAIL","error":e}), 1)),
            Err(_) => Ok((
                json!({"status":"FAIL","error":"native companion worker panicked"}),
                1,
            )),
        };
    }
    let op = args.first().map(String::as_str).unwrap_or("help");
    if op == "version" {
        return Ok((
            json!({"schema":"fixed-modulus-bench/v1","compiled_policy":POLICY,
        "git_commit":env!("ISOGENY_GIT_COMMIT"),"git_dirty":env!("ISOGENY_GIT_DIRTY")}),
            0,
        ));
    }
    if !matches!(op, "generate" | "measure" | "analyze") {
        return Err("fixed-modulus-bench generate --input MAP --out FRESH; measure --input DATASET --case CASE --out FRESH; analyze --input MANIFEST --out FRESH".into());
    }
    let allowed = if op == "measure" {
        vec!["input", "out", "case"]
    } else {
        vec!["input", "out"]
    };
    let o = options(&args[1..], &allowed)?;
    for key in ["input", "out"] {
        if !o.contains_key(key) {
            return Err(format!("--{key} required"));
        }
    }
    if op == "measure" && !o.contains_key("case") {
        return Err("--case required".into());
    }
    let out = Path::new(o.get("out").unwrap());
    fs::create_dir(out).map_err(|e| e.to_string())?;
    let out = fs::canonicalize(out).map_err(|e| e.to_string())?;
    let input = fs::canonicalize(o.get("input").unwrap()).map_err(|e| e.to_string())?;
    let ready = out.join("resources-ready");
    let mut worker = vec![
        "fixed-modulus-bench".into(),
        "__worker".into(),
        op.into(),
        "--input".into(),
        input.to_string_lossy().into_owned(),
        "--out".into(),
        out.to_string_lossy().into_owned(),
    ];
    if let Some(id) = o.get("case") {
        worker.extend(["--case".into(), id.clone()]);
    }
    worker.extend(["--rss-ready".into(), ready.to_string_lossy().into_owned()]);
    let receipt = supervisor::run(
        &std::env::current_exe().map_err(|e| e.to_string())?,
        &out,
        &worker,
        "native",
        if op == "measure" { 300 } else { 120 },
    );
    let value = serde_json::from_str::<Value>(&receipt.dump()).map_err(|e| e.to_string())?;
    let code = match value["status"].as_str() {
        Some("PASS") => 0,
        Some("TIMEOUT" | "MEMORY_LIMIT" | "RESOURCE_MONITOR_UNAVAILABLE") => 3,
        _ => 1,
    };
    Ok((
        json!({"schema":"fixed-modulus-supervised/v1","status":if code==0{"COMPLETED"}else{"FAILED"},
        "compiled_policy":POLICY,"resource_receipt":value,"public_admission":false}),
        code,
    ))
}
fn main() {
    match execute(&std::env::args().skip(1).collect::<Vec<_>>()) {
        Ok((value, code)) => {
            println!("{}", serde_json::to_string_pretty(&value).unwrap());
            std::process::exit(code)
        }
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn p192_and_small_prime_reduction_controls() {
        let p = BigUint::parse_bytes(P192.as_bytes(), 10).unwrap();
        poly::fixed_modulus_control_suite(&Field::new(&p).unwrap());
        kernel::fixed_modulus_kernel_controls(&Field::new(&p).unwrap());
        for p in [101u64, 107] {
            poly::fixed_modulus_control_suite(&Field::new(&BigUint::from(p)).unwrap());
        }
    }
    #[test]
    fn native_percentile_uses_maximum_of_eight_and_refuses_missing_samples() {
        let aa = vec![[100.0, 100.0]; 8];
        let abba = vec![[100.0, 80.0, 80.0, 100.0]; 8];
        assert_eq!(
            adjudicate(&aa, &abba).unwrap()["scoped_microcase_gate_passed"],
            true
        );
        let mut noisy = aa.clone();
        noisy[7] = [100.0, 120.0];
        assert_eq!(
            adjudicate(&noisy, &abba).unwrap()["scoped_microcase_gate_passed"],
            false
        );
        assert_eq!(adjudicate(&noisy, &abba).unwrap()["aa_rank"], 8);
        assert!(adjudicate(&aa[..7], &abba).is_err());
        let mut invalid = aa;
        invalid[0][0] = f64::NAN;
        assert!(adjudicate(&invalid, &abba).is_err());
    }
    #[test]
    fn deterministic_dataset_controls_are_bounded_and_match_both_reducers() {
        let f = Field::new(&BigUint::parse_bytes(P192.as_bytes(), 10).unwrap()).unwrap();
        for degree in [7, 127] {
            let mut m = stream(&f, "control", "modulus", 0, degree + 1);
            m[degree] = f.one();
            let c = make_case(&f, "control", "product", m.clone(), 2);
            let context = FixedModulus::new(&f, &m).unwrap();
            let result: Vec<Poly> = c
                .pairs
                .iter()
                .map(|[a, b]| context.mulmod(&decode(&f, a).unwrap(), &decode(&f, b).unwrap()))
                .collect();
            assert_eq!(output_digest(&f, &result), c.expected_sha256);
            assert_eq!(
                stream(&f, "control", "left", 0, 3),
                stream(&f, "control", "left", 0, 3)
            );
            assert_ne!(
                stream(&f, "control", "left", 0, 3),
                stream(&f, "control", "right", 0, 3)
            );
        }
    }
}
