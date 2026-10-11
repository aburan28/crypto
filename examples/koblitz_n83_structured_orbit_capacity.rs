//! Bounded construction-only gate for the stored N83 structured-orbit base.
//! The driver imposes a 240-second wall deadline and a 4 GiB, zero-swap
//! Linux cgroup; this process reads both limits back before model construction.
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::binary_semaev_chain_sat::{
    frobenius_orbit_domain_size, ChainedS3Encoding, FrobeniusOrbitDomainSize,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use flate2::read::GzDecoder;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::error::Error;
use std::fs;
use std::io::{BufRead, BufReader};
use std::time::Instant;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const OBJECT_BLAKE3: &str = "aea278aff2d79a592d5164d7437ddde8f8e9377fa481058dfc16afa176a31389";
const SET_BLAKE3: &str = "db3e7f25877c75dfd6f53db46447972ce6a44434cb7ea1e53e450d21af0fe79c";
const CURVE: &str = "icv1-f2m83-tm6151469093347-debefd74";
const LIMIT: u64 = 4 * 1024 * 1024 * 1024;
const SOURCE: &[u8] = include_bytes!("koblitz_n83_structured_orbit_capacity.rs");

fn decimal(v: &Value) -> Result<BigUint> {
    Ok(v.as_str().ok_or("decimal coordinate string")?.parse()?)
}

fn affine(v: &Value) -> Result<BinaryPoint> {
    let pair = v.as_array().ok_or("point pair")?;
    if pair.len() != 2 {
        return Err("point pair arity".into());
    }
    let x = decimal(&pair[0])?;
    let y = decimal(&pair[1])?;
    if x.bits() > 83 || y.bits() > 83 {
        return Err("coordinate outside F2^83".into());
    }
    Ok(BinaryPoint::Affine {
        x: F2mElement::from_biguint(&x, 83),
        y: F2mElement::from_biguint(&y, 83),
    })
}

fn cgroup_memory() -> Result<(u64, u64)> {
    #[cfg(not(target_os = "linux"))]
    return Err("Linux cgroup memory gate required".into());
    #[cfg(target_os = "linux")]
    {
        let maximum: u64 = fs::read_to_string("/sys/fs/cgroup/memory.max")?
            .trim()
            .parse()?;
        let swap: u64 = fs::read_to_string("/sys/fs/cgroup/memory.swap.max")?
            .trim()
            .parse()?;
        if maximum != LIMIT || swap != 0 {
            return Err("cgroup memory/swap readback mismatch".into());
        }
        Ok((maximum, swap))
    }
}

fn cgroup_peak_bytes() -> Result<u64> {
    #[cfg(not(target_os = "linux"))]
    return Err("Linux cgroup memory gate required".into());
    #[cfg(target_os = "linux")]
    return Ok(fs::read_to_string("/sys/fs/cgroup/memory.peak")?
        .trim()
        .parse()?);
}

fn run() -> Result<Value> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: capacity OBJECT.jsonl.gz MANIFEST.json PROBE-CORPUS.json".into());
    }
    let start = Instant::now();
    let manifest_bytes = fs::read(&args[2])?;
    let manifest: Value = serde_json::from_slice(&manifest_bytes)?;
    let compressed = fs::read(&args[1])?;
    if manifest["schema"] != "n83.structured-orbit-panel/v1"
        || manifest["compressed_blake3"] != OBJECT_BLAKE3
        || manifest["compressed_bytes"].as_u64() != Some(compressed.len() as u64)
        || manifest["historical_set_blake3"] != SET_BLAKE3
        || blake3::hash(&compressed).to_hex().as_str() != OBJECT_BLAKE3
    {
        return Err("stored object/manifest binding".into());
    }
    let mut lines = BufReader::new(GzDecoder::new(&compressed[..]));
    let mut header_line = String::new();
    lines.read_line(&mut header_line)?;
    let header: Value = serde_json::from_str(&header_line)?;
    if header["schema"] != "n83.structured-orbit-base/v1"
        || header["curve_slug"] != CURVE
        || header["orbit_columns"] != 2001
        || header["point_count"] != 332166
        || header["historical_set_blake3"] != SET_BLAKE3
    {
        return Err("structured-orbit header".into());
    }
    let rep_rows = header["representatives"]
        .as_array()
        .ok_or("representatives")?;
    if rep_rows.len() != 2001 {
        return Err("representative count".into());
    }
    let mut reps = Vec::with_capacity(rep_rows.len());
    let mut distinct = HashSet::new();
    for row in rep_rows {
        let BinaryPoint::Affine { x, .. } = affine(row)? else {
            return Err("affine representative".into());
        };
        let code = x.to_biguint();
        if !distinct.insert(code.clone()) {
            return Err("duplicate representative x".into());
        }
        reps.push(code);
    }
    let corpus: Value = serde_json::from_slice(&fs::read(&args[3])?)?;
    if corpus["schema"] != "n83.public-probe-corpus/v1" {
        return Err("public corpus schema".into());
    }
    let mut selected = corpus["targets"]
        .as_array()
        .ok_or("public targets")?
        .iter()
        .filter(|row| row["a"] == 0 && row["fixture"] == 0);
    let target = affine(&selected.next().ok_or("primary target")?["point"])?;
    if selected.next().is_some() {
        return Err("ambiguous primary target".into());
    }
    let kc = KoblitzCurve::known_n83_k0().ok_or("pinned K0")?;
    if kc.label() != CURVE
        || !kc.curve.is_on_curve(&target)
        || kc.mul(&target, &kc.subgroup_order) != BinaryPoint::Infinity
    {
        return Err("public target group check".into());
    }
    let BinaryPoint::Affine { x: target_x, .. } = &target else {
        return Err("affine target".into());
    };
    let import_ms = start.elapsed().as_secs_f64() * 1000.0;
    let expected = FrobeniusOrbitDomainSize {
        variables: 26_265,
        clauses: 872_260,
        xor_rows: 3_320,
    };
    if frobenius_orbit_domain_size(83, 5, reps.len()) != Some(expected) {
        return Err("domain-size preflight".into());
    }
    let (memory_limit, swap_limit) = cgroup_memory()?;
    let model_start = Instant::now();
    let mut model = ChainedS3Encoding::build(
        5,
        &kc.curve.irreducible,
        &kc.curve.b,
        Some(target_x),
        0,
        150_000,
    )
    .map_err(|_| "chained-S3 build failed")?;
    let core = (
        model.solver.n_vars(),
        model.solver.n_clauses(),
        model.solver.n_xors(),
    );
    if model
        .constrain_frobenius_orbit_summands(&reps, 150_000, 12_000_000)
        .map_err(|_| "Frobenius-domain install failed")?
        != expected
    {
        return Err("domain-size install".into());
    }
    let model_ms = model_start.elapsed().as_secs_f64() * 1000.0;
    if model.solver.n_vars() as u64 - core.0 as u64 != expected.variables
        || model.solver.n_clauses() - core.1 != expected.clauses
        || model.solver.n_xors() - core.2 != expected.xor_rows
    {
        return Err("installed model counts".into());
    }
    let source_commit = std::env::var("N83_SOURCE_COMMIT")?;
    if source_commit.len() != 40 || !source_commit.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err("source commit environment".into());
    }
    Ok(json!({
        "schema":"n83.structured-orbit-capacity-worker/v1",
        "status":"PASS_model_construction_only",
        "source_commit":source_commit,
        "source_blake3":blake3::hash(SOURCE).to_hex().to_string(),
        "chain_source_blake3":blake3::hash(include_bytes!("../src/cryptanalysis/binary_semaev_chain_sat.rs")).to_hex().to_string(),
        "object_compressed_blake3":OBJECT_BLAKE3,
        "manifest_blake3":blake3::hash(&manifest_bytes).to_hex().to_string(),
        "public_corpus_canonical_json_blake3":blake3::hash(&serde_json::to_vec(&corpus)?).to_hex().to_string(),
        "curve_slug":CURVE, "summands":5, "identity_mask":0,
        "orbit_columns":reps.len(), "domain_x_coordinates":166083,
        "max_variables":150000, "max_domain_clauses":12000000,
        "memory_cgroup_limit_bytes":memory_limit,
        "memory_cgroup_swap_limit_bytes":swap_limit,
        "memory_cgroup_peak_bytes":cgroup_peak_bytes()?,
        "core_variables":core.0, "core_clauses":core.1, "core_xors":core.2,
        "domain_variables":expected.variables, "domain_clauses":expected.clauses,
        "domain_xors":expected.xor_rows,
        "sat_variables":model.solver.n_vars(),
        "sat_clauses":model.solver.n_clauses(),
        "sat_xors":model.solver.n_xors(),
        "import_elapsed_ms_L0":import_ms,
        "model_construction_elapsed_ms_L0":model_ms,
        "total_elapsed_ms_L0":start.elapsed().as_secs_f64()*1000.0,
        "solver_search_executed":false,
        "relation_stage_executed":false,
        "rank_stage_executed":false,
        "total_index_calculus_runtime_ms":Value::Null,
        "selected_best_total_runtime":Value::Null
    }))
}

fn main() -> Result<()> {
    println!("{}", run()?);
    Ok(())
}
