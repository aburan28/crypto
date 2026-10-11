//! Import the frozen S3 structured-orbit object into the generic IC driver.
//!
//! `verify` checks the complete archive and factor-base maps without solver
//! search. `run` is a separate, hard-wall-capped experiment entry point; it
//! must be launched by a budgeted supervisor in a 4 GiB, zero-swap cgroup.
use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_explicit_frobenius_point_factor_base,
    koblitz_index_calculus_dlp_with_factor_base_and_progress, ChainS3Limits, DecompositionStrategy,
    FrobeniusFactorBase, KoblitzCurve, KoblitzIcOptions,
};
use flate2::read::GzDecoder;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::error::Error;
use std::fs::{self, OpenOptions};
use std::io::{Read, Write};
use std::path::Path;
use std::time::{Duration, Instant};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const STUDY: &str = "koblitz_n83_factor_base_sweep_20261008";
const CURVE: &str = "icv1-f2m83-tm6151469093347-debefd74";
const OBJECT_HASH: &str = "aea278aff2d79a592d5164d7437ddde8f8e9377fa481058dfc16afa176a31389";
const PLAIN_HASH: &str = "c137a82dada17a78f2f5cd2e8b48a1731a76fa4ca9ecfc8971b8f92b9821bc61";
const SET_HASH: &str = "db3e7f25877c75dfd6f53db46447972ce6a44434cb7ea1e53e450d21af0fe79c";
const MANIFEST_HASH: &str = "d6adad9bec74350f37bf917e737006869bc0f343692f668fa4dc893afe332745";
const OBJECT_NAME: &str =
    "objects/aea278aff2d79a592d5164d7437ddde8f8e9377fa481058dfc16afa176a31389.jsonl.gz";
const S3_URI: &str = "s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v3-structured-orbit/a0/objects/aea278aff2d79a592d5164d7437ddde8f8e9377fa481058dfc16afa176a31389.jsonl.gz";
const MANIFEST: &[u8] = include_bytes!("../research/koblitz_n83_factor_base_sweep_20261008/verification/structured-d12-orbit-20261010/manifest.json");
const EXPORTER: &[u8] = include_bytes!("koblitz_n83_structured_orbit_export.rs");
const DESIGN: &[u8] =
    include_bytes!("../research/koblitz_n83_factor_base_sweep_20261008/structured-orbit-v3.json");

fn decimal(v: &Value) -> Result<BigUint> {
    let s = v
        .as_str()
        .ok_or("decimal coordinate or coefficient required")?;
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return Err("noncanonical decimal field".into());
    }
    let n: BigUint = s.parse()?;
    if n.to_string() != s {
        return Err("noncanonical decimal field".into());
    }
    Ok(n)
}

fn point(v: &Value) -> Result<BinaryPoint> {
    let pair = v.as_array().ok_or("affine point pair required")?;
    if pair.len() != 2 {
        return Err("affine point arity".into());
    }
    let (x, y) = (decimal(&pair[0])?, decimal(&pair[1])?);
    if x.bits() > 83 || y.bits() > 83 {
        return Err("coordinate outside F2^83".into());
    }
    Ok(BinaryPoint::Affine {
        x: F2mElement::from_biguint(&x, 83),
        y: F2mElement::from_biguint(&y, 83),
    })
}

fn words(p: &BinaryPoint) -> Result<(u128, u128)> {
    let BinaryPoint::Affine { x, y } = p else {
        return Err("point at infinity in factor base".into());
    };
    let u128_of = |value: &F2mElement| -> Result<u128> {
        let limbs = value.to_biguint().to_u64_digits();
        if limbs.len() > 2 {
            return Err("wide factor-base coordinate".into());
        }
        Ok(u128::from(*limbs.first().unwrap_or(&0))
            | (u128::from(*limbs.get(1).unwrap_or(&0)) << 64))
    };
    Ok((u128_of(x)?, u128_of(y)?))
}

fn set_hash(points: &HashSet<(u128, u128)>) -> String {
    let mut sorted: Vec<_> = points.iter().copied().collect();
    sorted.sort_unstable_by_key(|&(x, y)| {
        let mut tag = [0u8; 33];
        tag[0] = 1;
        tag[1..17].copy_from_slice(&x.to_le_bytes());
        tag[17..33].copy_from_slice(&y.to_le_bytes());
        tag
    });
    let mut h = blake3::Hasher::new();
    for (x, y) in sorted {
        h.update(&[1]);
        h.update(&x.to_le_bytes());
        h.update(&y.to_le_bytes());
    }
    h.finalize().to_hex().to_string()
}

fn bound_archive(manifest_bytes: &[u8], compressed: &[u8]) -> Result<Vec<u8>> {
    if blake3::hash(MANIFEST).to_hex().as_str() != MANIFEST_HASH
        || manifest_bytes != MANIFEST
        || blake3::hash(compressed).to_hex().as_str() != OBJECT_HASH
    {
        return Err("pinned manifest or compressed S3 object mismatch".into());
    }
    let manifest: Value = serde_json::from_slice(manifest_bytes)?;
    if manifest["schema"] != "n83.structured-orbit-panel/v1"
        || manifest["status"] != "completed_factor_base_object"
        || manifest["study"] != STUDY
        || manifest["object"] != OBJECT_NAME
        || manifest["s3_uri"] != S3_URI
        || manifest["compressed_blake3"] != OBJECT_HASH
        || manifest["plain_blake3"] != PLAIN_HASH
        || manifest["historical_set_blake3"] != SET_HASH
        || manifest["compressed_bytes"].as_u64() != Some(compressed.len() as u64)
        || manifest["source_blake3"] != blake3::hash(EXPORTER).to_hex().as_str()
        || manifest["design_blake3"] != blake3::hash(DESIGN).to_hex().as_str()
        || manifest["orbit_columns"] != 2001
        || manifest["point_records"] != 332166
    {
        return Err("manifest binding".into());
    }
    let mut plain = Vec::new();
    GzDecoder::new(compressed).read_to_end(&mut plain)?;
    if blake3::hash(&plain).to_hex().as_str() != PLAIN_HASH {
        return Err("uncompressed object mismatch".into());
    }
    Ok(plain)
}

fn decode_base(plain: &[u8], kc: &KoblitzCurve) -> Result<FrobeniusFactorBase> {
    let text = std::str::from_utf8(plain)?;
    let mut lines = text.lines();
    let header: Value = serde_json::from_str(lines.next().ok_or("missing header")?)?;
    if header["schema"] != "n83.structured-orbit-base/v1"
        || header["study"] != STUDY
        || header["curve_slug"] != CURVE
        || header["n"] != 83
        || header["curve_a"] != 0
        || header["curve_b"] != "1"
        || header["modulus_low_terms"] != json!([0, 1, 2, 45])
        || header["subgroup_order"] != kc.subgroup_order.to_string()
        || header["cofactor"] != kc.cofactor.to_string()
        || header["frobenius_lambda"] != kc.lambda.to_string()
        || header["source_dimension"] != 12
        || header["closure"] != "signed_frobenius"
        || header["orbit_columns"] != 2001
        || header["point_count"] != 332166
        || header["historical_set_blake3"] != SET_HASH
    {
        return Err("header or curve binding".into());
    }
    let reps = header["representatives"]
        .as_array()
        .ok_or("representatives")?;
    let sources = header["source_points_before_cofactor_clearing"]
        .as_array()
        .ok_or("sources")?;
    let codes = header["source_x_codes"].as_array().ok_or("source codes")?;
    if reps.len() != 2001 || sources.len() != reps.len() || codes.len() != reps.len() {
        return Err("representative and source lengths".into());
    }
    let mut points = Vec::with_capacity(332166);
    let mut seen = HashSet::with_capacity(332166);
    let mut previous_code = None;
    for (column, ((rep_row, source_row), code_row)) in
        reps.iter().zip(sources).zip(codes).enumerate()
    {
        let source = point(source_row)?;
        let rep = point(rep_row)?;
        let code = code_row.as_u64().ok_or("source code")?;
        if code >= 4096
            || previous_code.is_some_and(|old| old >= code)
            || words(&source)?.0 != u128::from(code)
            || !kc.curve.is_on_curve(&source)
            || kc.mul(&source, &kc.cofactor) != rep
            || !kc.curve.is_on_curve(&rep)
            || kc.mul(&rep, &kc.subgroup_order) != BinaryPoint::Infinity
            || kc.mul(&rep, &kc.lambda) != kc.frobenius(&rep)
        {
            return Err("source projection or subgroup certificate".into());
        }
        previous_code = Some(code);
        let mut current = rep.clone();
        let mut coefficient = BigUint::from(1u8);
        for phase in 0..83u64 {
            for negative in [false, true] {
                let row: Value = serde_json::from_str(lines.next().ok_or("missing point row")?)?;
                let expected = if negative {
                    point_neg(&current)
                } else {
                    current.clone()
                };
                let label = if negative {
                    &kc.subgroup_order - &coefficient
                } else {
                    coefficient.clone()
                };
                if point(&row["point"])? != expected
                    || row["column"].as_u64() != Some(column as u64)
                    || row["phase"].as_u64() != Some(phase)
                    || row["negative"].as_bool() != Some(negative)
                    || decimal(&row["coefficient"])? != label
                    || !seen.insert(words(&expected)?)
                {
                    return Err("point, orbit, coefficient or uniqueness mismatch".into());
                }
                points.push(expected);
            }
            current = kc.frobenius(&current);
            coefficient = (&coefficient * &kc.lambda) % &kc.subgroup_order;
        }
        if current != rep {
            return Err("orbit did not close".into());
        }
    }
    if lines.next().is_some() || points.len() != 332166 || set_hash(&seen) != SET_HASH {
        return Err("trailing rows, point count or point-set hash".into());
    }
    let fb = build_explicit_frobenius_point_factor_base(kc, points)
        .ok_or("generic orbit reconstruction failed")?;
    if fb.signed_orbits.len() != 2001
        || fb.orbits.len() != 4002
        || fb.subspace.len() != 166083
        || fb
            .signed_orbit_of
            .iter()
            .enumerate()
            .any(|(i, &actual)| actual != (i / 166, ((i % 166) / 2) as u32, i % 2 == 1))
    {
        return Err("generic orbit map differs from archive labels".into());
    }
    Ok(fb)
}

fn import(manifest_path: &Path, object_path: &Path) -> Result<(KoblitzCurve, FrobeniusFactorBase)> {
    let manifest_bytes = fs::read(manifest_path)?;
    let compressed = fs::read(object_path)?;
    let plain = bound_archive(&manifest_bytes, &compressed)?;
    let kc = KoblitzCurve::known_n83_k0().ok_or("pinned K0 curve")?;
    if kc.label() != CURVE || kc.cofactor != BigUint::from(4u8) {
        return Err("pinned curve mismatch".into());
    }
    let fb = decode_base(&plain, &kc)?;
    Ok((kc, fb))
}

fn public_target(corpus_path: &Path, kc: &KoblitzCurve) -> Result<(BinaryPoint, String)> {
    let corpus: Value = serde_json::from_slice(&fs::read(corpus_path)?)?;
    if corpus["schema"] != "n83.public-probe-corpus/v1" || corpus["study"] != STUDY {
        return Err("public corpus schema or study".into());
    }
    let selected: Vec<_> = corpus["targets"]
        .as_array()
        .ok_or("public targets")?
        .iter()
        .filter(|row| row["a"] == 0 && row["fixture"] == 0)
        .collect();
    if selected.len() != 1 {
        return Err("ambiguous primary target".into());
    }
    let target = point(&selected[0]["point"])?;
    if !kc.curve.is_on_curve(&target)
        || kc.mul(&target, &kc.subgroup_order) != BinaryPoint::Infinity
    {
        return Err("target subgroup check".into());
    }
    Ok((
        target,
        blake3::hash(&serde_json::to_vec(&corpus)?)
            .to_hex()
            .to_string(),
    ))
}

fn validate_answer(validation_path: &Path, corpus_hash: &str, answer: &BigUint) -> Result<bool> {
    let record: Value = serde_json::from_slice(&fs::read(validation_path)?)?;
    if record["schema"] != "n83.public-probe-validation/v1"
        || record["corpus_canonical_json_blake3"] != corpus_hash
        || record["oracle_receives_known_answer_scalars"] != false
    {
        return Err("validation/corpus binding".into());
    }
    let selected: Vec<_> = record["validation"]
        .as_array()
        .ok_or("validation rows")?
        .iter()
        .filter(|row| row["a"] == 0 && row["fixture"] == 0)
        .collect();
    if selected.len() != 1 {
        return Err("ambiguous validation target".into());
    }
    Ok(decimal(&selected[0]["known_answer_scalar"])? == *answer)
}

fn write_new(path: &Path, value: &Value) -> Result<()> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    serde_json::to_writer_pretty(&mut file, value)?;
    file.write_all(b"\n")?;
    file.sync_all()?;
    Ok(())
}

fn cgroup_gate() -> Result<()> {
    #[cfg(not(target_os = "linux"))]
    return Err("run requires a Linux cgroup".into());
    #[cfg(target_os = "linux")]
    {
        if fs::read_to_string("/sys/fs/cgroup/memory.max")?.trim() != "4294967296"
            || fs::read_to_string("/sys/fs/cgroup/memory.swap.max")?.trim() != "0"
        {
            return Err("run requires exactly 4 GiB memory and zero swap".into());
        }
        Ok(())
    }
}

fn cgroup_peak_bytes() -> Result<u64> {
    #[cfg(not(target_os = "linux"))]
    return Err("run requires a Linux cgroup".into());
    #[cfg(target_os = "linux")]
    {
        let peak: u64 = fs::read_to_string("/sys/fs/cgroup/memory.peak")?
            .trim()
            .parse()?;
        if peak == 0 || peak > 4 * 1024 * 1024 * 1024 {
            return Err("cgroup peak outside declared memory limit".into());
        }
        Ok(peak)
    }
}

fn run(args: &[String]) -> Result<()> {
    if args.len() != 10 {
        return Err("usage: run MANIFEST OBJECT CORPUS VALIDATION BUDGET_SECONDS MAX_TRIALS CONFLICT_BUDGET MAX_MODELS OUTPUT_DIR".into());
    }
    let budget: u64 = args[5].parse()?;
    let max_trials: usize = args[6].parse()?;
    let conflicts: u64 = args[7].parse()?;
    let models: usize = args[8].parse()?;
    if !(1..=7200).contains(&budget)
        || max_trials == 0
        || !(1..=10_000_000).contains(&conflicts)
        || !(1..=1024).contains(&models)
    {
        return Err("invalid hard work cap".into());
    }
    cgroup_gate()?;
    let start = Instant::now();
    let output = Path::new(&args[9]);
    fs::create_dir(output)?;
    write_new(
        &output.join("config.json"),
        &json!({
            "schema":"n83.structured-orbit-cold-config/v1", "study":STUDY,
            "object_blake3":OBJECT_HASH, "manifest_blake3":MANIFEST_HASH,
            "public_corpus_path":args[3], "validation_path":args[4],
            "budget_seconds":budget, "max_trials":max_trials,
            "conflict_budget":conflicts, "max_models":models,
            "summands":5, "strategy":"ChainedS3FrobeniusOrbit",
            "memory_cgroup_limit_bytes":4294967296u64,
            "memory_cgroup_swap_limit_bytes":0,
            "source_blake3":blake3::hash(include_bytes!("koblitz_n83_structured_orbit_worker.rs")).to_hex().to_string()
        }),
    )?;
    let cap_path = output.join("cap.json");
    std::thread::spawn(move || {
        std::thread::sleep(Duration::from_secs(budget));
        let _ = write_new(
            &cap_path,
            &json!({
                "schema":"n83.structured-orbit-cold-cap/v1",
                "status":"UNKNOWN_budget", "budget_seconds":budget,
                "elapsed_process_wall_ms":start.elapsed().as_secs_f64()*1000.0
            }),
        );
        std::process::exit(124);
    });
    let (kc, fb) = import(Path::new(&args[1]), Path::new(&args[2]))?;
    let import_ms = start.elapsed().as_secs_f64() * 1000.0;
    let (target, corpus_hash) = public_target(Path::new(&args[3]), &kc)?;
    let mut opts = KoblitzIcOptions::default();
    opts.m = 5;
    opts.strategy = DecompositionStrategy::ChainedS3FrobeniusOrbit;
    opts.max_trials = max_trials;
    opts.seed = 2026100901;
    opts.allow_direct_relation = false;
    opts.chain_s3_limits = ChainS3Limits {
        max_variables: 150_000,
        max_domain_clauses: 12_000_000,
        max_models: models,
        conflict_budget: conflicts,
    };
    let mut events = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output.join("events.jsonl"))?;
    let mut write_failed = false;
    let report = koblitz_index_calculus_dlp_with_factor_base_and_progress(
        &kc,
        &target,
        &fb,
        &opts,
        &mut |event| {
            if writeln!(
                events,
                "{}",
                json!({
                    "elapsed_process_wall_ms":start.elapsed().as_secs_f64()*1000.0,
                    "event":format!("{event:?}")
                })
            )
            .and_then(|_| events.flush())
            .is_err()
            {
                write_failed = true;
            }
        },
    )
    .ok_or("generic IC driver refused imported base")?;
    if write_failed {
        return Err("progress-event write failed".into());
    }
    let valid = match &report.log {
        Some(answer) => {
            kc.mul(kc.generator(), answer) == target
                && validate_answer(Path::new(&args[4]), &corpus_hash, answer)?
        }
        None => false,
    };
    let bad = report.inconsistent_relations != 0
        || report.verification_failures != 0
        || report.sat_invalid_models != 0
        || report.direct_relation
        || (report.log.is_some() && !valid);
    let memory_peak = cgroup_peak_bytes()?;
    let status = if bad {
        "PRODUCER_FAILURE"
    } else if valid {
        "PASS_verified_target"
    } else if report.sat_unknowns != 0 {
        "UNKNOWN_solver_cap"
    } else {
        "UNKNOWN_trial_cap"
    };
    let attempts: Vec<_> = report
        .attempt_records
        .iter()
        .map(|r| {
            json!({
                "trial":r.trial, "a":r.coefficient_a.to_string(),
                "b":r.coefficient_b.to_string(), "target":words(&r.target).ok(),
                "disposition":format!("{:?}",r.disposition),
                "decomposition_indices":r.decomposition_indices,
                "solver_calls":r.solver_calls, "models":r.models,
                "conflicts":r.conflicts, "implied_rows":r.implied_rows
            })
        })
        .collect();
    let relations: Vec<_> = report
        .relation_matrix
        .iter()
        .map(|r| {
            json!({
                "a":r.coef_a.to_string(), "b":r.coef_b.to_string(),
                "summands":r.summands, "summand_negated":r.summand_negated,
                "row":r.row.iter().map(ToString::to_string).collect::<Vec<_>>()
            })
        })
        .collect();
    write_new(
        &output.join("result.json"),
        &json!({
            "schema":"n83.structured-orbit-cold-result/v1", "study":STUDY,
            "status":status, "object_blake3":OBJECT_HASH,
            "manifest_blake3":MANIFEST_HASH,
            "public_corpus_canonical_json_blake3":corpus_hash,
            "factor_base_size":report.factor_base_size, "orbit_count":report.orbit_count,
            "trials":report.trials, "relations":report.relations,
            "independent_relations":report.independent_relations,
            "dependent_relations":report.dependent_relations,
            "inconsistent_relations":report.inconsistent_relations,
            "verification_failures":report.verification_failures,
            "sat_calls":report.sat_calls, "sat_unknowns":report.sat_unknowns,
            "sat_invalid_models":report.sat_invalid_models,
            "sat_models":report.sat_models, "sat_conflicts":report.sat_conflicts,
            "matrix_rows":report.matrix_rows, "matrix_columns":report.matrix_columns,
            "terminal_matrix_rank":report.terminal_matrix_rank,
            "linear_solve_attempts":report.linear_solve_attempts,
            "rank_history":report.rank_history.iter().map(|r| json!({
                "rows":r.rows,"columns":r.columns,"rank":r.rank,
                "candidate_produced":r.candidate_produced,
                "candidate_verified":r.candidate_verified
            })).collect::<Vec<_>>(),
            "attempt_records":attempts, "relation_matrix":relations,
            "candidate_log":report.log.as_ref().map(ToString::to_string),
            "candidate_verified":valid,
            "memory_cgroup_peak_bytes":memory_peak,
            "base_import_elapsed_ms_L0":import_ms,
            "relation_collection_ns_L0":report.relation_collection_ns.to_string(),
            "linear_algebra_ns_L0":report.linear_algebra_ns.to_string(),
            "worker_wall_ms_L0":start.elapsed().as_secs_f64()*1000.0,
            "total_index_calculus_runtime_ms":Value::Null,
            "selected_best_total_runtime":Value::Null
        }),
    )?;
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("verify") if args.len() == 3 => {
            let (kc, fb) = import(Path::new(&args[1]), Path::new(&args[2]))?;
            println!("{}", json!({
                "schema":"n83.structured-orbit-import/v1", "status":"PASS",
                "curve_slug":kc.label(), "object_blake3":OBJECT_HASH,
                "points":fb.points.len(), "signed_orbits":fb.signed_orbits.len(),
                "x_coordinates":fb.subspace.len(),
                "solver_search_executed":false
            }));
            Ok(())
        }
        Some("run") => run(&args),
        _ => Err("usage: verify MANIFEST OBJECT | run MANIFEST OBJECT CORPUS VALIDATION BUDGET_SECONDS MAX_TRIALS CONFLICT_BUDGET MAX_MODELS OUTPUT_DIR".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_manifest_and_archive_require_exact_bytes() {
        assert_eq!(blake3::hash(MANIFEST).to_hex().as_str(), MANIFEST_HASH);
        assert!(bound_archive(MANIFEST, b"changed gzip object").is_err());
        let mut changed = MANIFEST.to_vec();
        changed.push(b' ');
        assert!(bound_archive(&changed, &[]).is_err());
    }

    #[test]
    fn decimal_and_affine_fields_reject_noncanonical_values() {
        assert!(decimal(&json!("01")).is_err());
        assert!(decimal(&json!("-1")).is_err());
        assert!(point(&json!(["1", "1"])).is_ok());
        assert!(point(&json!(["9671406556917033397649408", "0"])).is_err());
    }
}
