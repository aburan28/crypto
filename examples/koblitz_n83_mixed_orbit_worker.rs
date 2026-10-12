//! Replay two pinned N83 bases and collect bounded one/two-residual relations.
//!
//! `verify` and `verify-alt` import pinned S3 objects and check their exact
//! orbit-labelled prefix. `run` and `run-alt` additionally require a 4 GiB, zero-swap Linux cgroup and
//! an external budget supervisor. Neither mode downloads an object.
//! `verify-control` and `control-run` use a pinned known-answer witness only;
//! their outcomes are not natural-target relation-yield measurements.

use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_explicit_frobenius_point_factor_base, ChainS3Limits, FrobeniusFactorBase, KoblitzCurve,
};
use crypto_lib::cryptanalysis::koblitz_large_prime::{
    build_exact_residual_suffix_factor_base, FeedOutcome, LargePrimeLimits, MixedOrbitDomains,
    WideLargePrimeEliminator,
};
use flate2::read::GzDecoder;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
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
const PREFIX: &str = "s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/";
const PILOT_PREFIX: &str =
    "s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/";
const SMALL_MANIFEST: &[u8] = include_bytes!("../research/koblitz_n83_factor_base_sweep_20261008/verification/matched-k1182-k2048-900s-20261011/evidence-01/panel-k1182/manifest.json");
const SMALL_600_MANIFEST: &[u8] =
    include_bytes!("../research/koblitz_n83_factor_base_sweep_20261008/pilot-01/manifest.json");
const ENVELOPE_MANIFEST: &[u8] = include_bytes!("../research/koblitz_n83_factor_base_sweep_20261008/verification/remaining-budget-k2048-20261010/manifest.json");
const SMALL_HASH: &str = "76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92";
const SMALL_600_HASH: &str = "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e";
const ENVELOPE_HASH: &str = "82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851";
const PUBLIC_CORPUS_HASH: &str = "831eaa7be2387634fd610f6706103b2af9b4e274f0db59923c0280af3cbed74b";
const CONTROL_FIXTURE: &[u8] = include_bytes!("../research/koblitz_n83_factor_base_sweep_20261008/verification/mixed-orbit-known-witness-source-20261012/fixture.json");

struct PinnedBase {
    manifest: &'static [u8],
    object_hash: &'static str,
    columns: usize,
    panel_schema: &'static str,
    panel_status: &'static str,
    panel_count: usize,
    s3_prefix: &'static str,
}

const SMALL: PinnedBase = PinnedBase {
    manifest: SMALL_MANIFEST,
    object_hash: SMALL_HASH,
    columns: 1182,
    panel_schema: "n83.factor-base-panel/v2-size-frontier",
    panel_status: "completed_factor_base_object",
    panel_count: 1,
    s3_prefix: PREFIX,
};
const SMALL_600: PinnedBase = PinnedBase {
    manifest: SMALL_600_MANIFEST,
    object_hash: SMALL_600_HASH,
    columns: 600,
    panel_schema: "n83.factor-base-panel/v1",
    panel_status: "completed_factor_base_panel",
    panel_count: 54,
    s3_prefix: PILOT_PREFIX,
};
const ENVELOPE: PinnedBase = PinnedBase {
    manifest: ENVELOPE_MANIFEST,
    object_hash: ENVELOPE_HASH,
    columns: 2048,
    panel_schema: "n83.factor-base-panel/v2-size-frontier",
    panel_status: "completed_factor_base_object",
    panel_count: 1,
    s3_prefix: PREFIX,
};

fn decimal(value: &Value) -> Result<BigUint> {
    let raw = value.as_str().ok_or("decimal string required")?;
    if raw.is_empty() || !raw.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err("invalid decimal encoding".into());
    }
    let number: BigUint = raw.parse()?;
    if number.to_string() != raw {
        return Err("noncanonical decimal encoding".into());
    }
    Ok(number)
}

fn point(value: &Value) -> Result<BinaryPoint> {
    let pair = value.as_array().ok_or("affine pair required")?;
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

fn point_json(point: &BinaryPoint) -> Result<Value> {
    let BinaryPoint::Affine { x, y } = point else {
        return Err("identity is not an affine point".into());
    };
    Ok(json!([
        x.to_biguint().to_string(),
        y.to_biguint().to_string()
    ]))
}

fn word_pair(point: &BinaryPoint) -> Result<(u128, u128)> {
    let BinaryPoint::Affine { x, y } = point else {
        return Err("identity in factor base".into());
    };
    Ok((
        x.to_biguint().to_u128().ok_or("wide x coordinate")?,
        y.to_biguint().to_u128().ok_or("wide y coordinate")?,
    ))
}

fn point_set_hash(mut points: Vec<(u128, u128)>) -> String {
    points.sort_unstable();
    let mut hash = blake3::Hasher::new();
    hash.update(b"n83.sorted-polynomial-points/v1");
    for (x, y) in points {
        hash.update(&x.to_le_bytes());
        hash.update(&y.to_le_bytes());
    }
    hash.finalize().to_hex().to_string()
}

fn bound_archive(
    spec: &PinnedBase,
    manifest: &[u8],
    compressed: &[u8],
) -> Result<(Value, Vec<u8>)> {
    if manifest != spec.manifest || blake3::hash(compressed).to_hex().as_str() != spec.object_hash {
        return Err("pinned manifest or compressed object mismatch".into());
    }
    let header: Value = serde_json::from_slice(manifest)?;
    let rows = header["bases"].as_array().ok_or("manifest bases")?;
    let object = format!("objects/{}.jsonl.gz", spec.object_hash);
    if header["schema"] != spec.panel_schema
        || header["status"] != spec.panel_status
        || header["study"] != STUDY
        || header["completed_base_count"] != spec.panel_count
        || rows.len() != spec.panel_count
        || (spec.panel_count > 1 && header["expected_base_count"] != spec.panel_count)
    {
        return Err("manifest panel binding".into());
    }
    let selected: Vec<_> = rows
        .iter()
        .filter(|row| row["compressed_blake3"] == spec.object_hash)
        .collect();
    if selected.len() != 1 {
        return Err("manifest selected row binding".into());
    }
    let row = selected[0];
    if row["a"] != 0
        || row["policy"] != "public_x_hash"
        || row["seed"] != 2026100801u64
        || row["columns"] != spec.columns
        || row["points"] != spec.columns * 166
        || row["object"] != object
        || row["compressed_blake3"] != spec.object_hash
        || row["s3_uri"] != format!("{}{object}", spec.s3_prefix)
        || row["bytes"].as_u64() != Some(compressed.len() as u64)
    {
        return Err("manifest row binding".into());
    }
    let mut plain = Vec::new();
    GzDecoder::new(compressed).read_to_end(&mut plain)?;
    if blake3::hash(&plain).to_hex().as_str() != row["plain_blake3"].as_str().ok_or("plain hash")? {
        return Err("plain object hash mismatch".into());
    }
    Ok((row.clone(), plain))
}

fn decode_base(row: &Value, plain: &[u8], curve: &KoblitzCurve) -> Result<FrobeniusFactorBase> {
    let text = std::str::from_utf8(plain)?;
    let mut lines = text.lines();
    let header: Value = serde_json::from_str(lines.next().ok_or("missing base header")?)?;
    let columns = row["columns"].as_u64().ok_or("columns")? as usize;
    if header["schema"] != "n83.factor-base/v1"
        || header["study"] != STUDY
        || header["curve_slug"] != CURVE
        || header["n"] != 83
        || header["curve_a"] != 0
        || header["curve_b"] != "1"
        || header["modulus_low_terms"] != json!([0, 1, 2, 45])
        || header["policy"] != row["policy"]
        || header["seed"] != row["seed"]
        || header["orbit_columns"] != row["columns"]
        || header["point_count"] != row["points"]
        || header["closure"] != "signed_frobenius"
        || header["subgroup_order"] != curve.subgroup_order.to_string()
        || header["cofactor"] != curve.cofactor.to_string()
        || header["frobenius_lambda"] != curve.lambda.to_string()
        || point(&header["generator"])? != *curve.generator()
    {
        return Err("base header and curve mismatch".into());
    }
    let reps = header["representatives"]
        .as_array()
        .ok_or("representatives")?;
    let sources = header["source_points_before_cofactor_clearing"]
        .as_array()
        .ok_or("source points")?;
    let indices = header["accepted_candidate_indices"]
        .as_array()
        .ok_or("candidate indices")?;
    if reps.len() != columns || sources.len() != columns || indices.len() != columns {
        return Err("representative length mismatch".into());
    }
    let mut points = Vec::with_capacity(columns * 166);
    let mut seen = HashSet::with_capacity(columns * 166);
    let mut previous_index = 0;
    for (column, ((rep, source), index)) in reps.iter().zip(sources).zip(indices).enumerate() {
        let rep = point(rep)?;
        let source = point(source)?;
        let candidate = index.as_u64().ok_or("candidate index")?;
        if candidate <= previous_index
            || !curve.curve.is_on_curve(&source)
            || curve.mul(&source, &curve.cofactor) != rep
            || !curve.curve.is_on_curve(&rep)
            || curve.mul(&rep, &curve.subgroup_order) != BinaryPoint::Infinity
            || curve.mul(&rep, &curve.lambda) != curve.frobenius(&rep)
        {
            return Err("source projection or representative mismatch".into());
        }
        previous_index = candidate;
        let mut current = rep.clone();
        let mut coefficient = BigUint::from(1u8);
        for phase in 0..83u64 {
            for negative in [false, true] {
                let archive_row: Value =
                    serde_json::from_str(lines.next().ok_or("missing point row")?)?;
                let expected = if negative {
                    point_neg(&current)
                } else {
                    current.clone()
                };
                let label = if negative {
                    &curve.subgroup_order - &coefficient
                } else {
                    coefficient.clone()
                };
                if point(&archive_row["point"])? != expected
                    || archive_row["column"].as_u64() != Some(column as u64)
                    || archive_row["phase"].as_u64() != Some(phase)
                    || archive_row["negative"].as_bool() != Some(negative)
                    || decimal(&archive_row["coefficient"])? != label
                    || !seen.insert(word_pair(&expected)?)
                {
                    return Err("point row or orbit label mismatch".into());
                }
                points.push(expected);
            }
            current = curve.frobenius(&current);
            coefficient = (&coefficient * &curve.lambda) % &curve.subgroup_order;
        }
        if current != rep {
            return Err("Frobenius orbit did not close".into());
        }
    }
    if lines.next().is_some()
        || previous_index != row["raw_x_candidates"].as_u64().ok_or("candidate count")?
        || points.len() != columns * 166
        || point_set_hash(seen.into_iter().collect())
            != row["point_set_blake3"].as_str().ok_or("point set hash")?
    {
        return Err("trailing rows, candidate count or point set mismatch".into());
    }
    let base = build_explicit_frobenius_point_factor_base(curve, points)
        .ok_or("generic orbit reconstruction failed")?;
    if base.signed_orbits.len() != columns
        || base.orbits.len() != columns * 2
        || base
            .signed_orbit_of
            .iter()
            .enumerate()
            .any(|(index, &label)| {
                label != (index / 166, ((index % 166) / 2) as u32, index % 2 == 1)
            })
    {
        return Err("generic orbit map differs from archive".into());
    }
    Ok(base)
}

fn import(
    spec: &PinnedBase,
    manifest: &Path,
    object: &Path,
    curve: &KoblitzCurve,
) -> Result<FrobeniusFactorBase> {
    let (row, plain) = bound_archive(spec, &fs::read(manifest)?, &fs::read(object)?)?;
    decode_base(&row, &plain, curve)
}

fn import_pair(
    args: &[String],
    small_spec: &PinnedBase,
) -> Result<(KoblitzCurve, FrobeniusFactorBase, FrobeniusFactorBase)> {
    let curve = KoblitzCurve::known_n83_k0().ok_or("pinned Koblitz curve")?;
    if curve.label() != CURVE || curve.cofactor != BigUint::from(4u8) {
        return Err("curve identity mismatch".into());
    }
    let small = import(small_spec, Path::new(&args[1]), Path::new(&args[2]), &curve)?;
    let envelope = import(&ENVELOPE, Path::new(&args[3]), Path::new(&args[4]), &curve)?;
    let residual = build_exact_residual_suffix_factor_base(&curve, &small, &envelope)?;
    let expected_residual_orbits = ENVELOPE.columns - small_spec.columns;
    if residual.signed_orbits.len() != expected_residual_orbits
        || residual.points.len() != expected_residual_orbits * 166
    {
        return Err("residual suffix size mismatch".into());
    }
    Ok((curve, small, residual))
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
            return Err("run requires 4 GiB memory and zero swap".into());
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
            return Err("invalid cgroup memory peak".into());
        }
        Ok(peak)
    }
}

fn public_target(path: &Path, curve: &KoblitzCurve) -> Result<(BinaryPoint, String)> {
    let corpus: Value = serde_json::from_slice(&fs::read(path)?)?;
    let digest = blake3::hash(&serde_json::to_vec(&corpus)?)
        .to_hex()
        .to_string();
    if digest != PUBLIC_CORPUS_HASH {
        return Err("public corpus differs from matched N83 fixture".into());
    }
    if corpus["schema"] != "n83.public-probe-corpus/v1" || corpus["study"] != STUDY {
        return Err("public corpus mismatch".into());
    }
    let rows = corpus["targets"].as_array().ok_or("public targets")?;
    let selected: Vec<_> = rows
        .iter()
        .filter(|row| row["a"] == 0 && row["fixture"] == 0)
        .collect();
    if selected.len() != 1 {
        return Err("ambiguous public target".into());
    }
    let target = point(&selected[0]["point"])?;
    if !curve.curve.is_on_curve(&target)
        || curve.mul(&target, &curve.subgroup_order) != BinaryPoint::Infinity
    {
        return Err("target subgroup check".into());
    }
    Ok((target, digest))
}

fn trial_coefficient(label: &[u8], trial: u64, modulus: &BigUint) -> BigUint {
    let mut hash = blake3::Hasher::new();
    hash.update(b"n83.mixed-orbit-partial-trial/v1");
    hash.update(label);
    hash.update(&trial.to_le_bytes());
    BigUint::from_bytes_le(hash.finalize().as_bytes()) % modulus
}

fn parse_control_bytes(bytes: &[u8]) -> Result<Value> {
    if bytes != CONTROL_FIXTURE {
        return Err("known-answer fixture differs from pinned bytes".into());
    }
    let fixture: Value = serde_json::from_slice(bytes)?;
    if fixture["schema"] != "n83.mixed-orbit-known-witness/v1"
        || fixture["scope"] != "known_answer_control_only"
        || fixture["curve_slug"] != CURVE
        || fixture["equation"]["a"] != "0"
        || fixture["equation"]["b"] != "1"
        || fixture["equation"]["small_summands"] != 3
        || fixture["equation"]["residual_summands"] != 2
        || fixture["archive"]["small"]["s3_uri"] != format!("{PREFIX}objects/{SMALL_HASH}.jsonl.gz")
        || fixture["archive"]["envelope"]["s3_uri"]
            != format!("{PREFIX}objects/{ENVELOPE_HASH}.jsonl.gz")
    {
        return Err("known-answer fixture metadata".into());
    }
    Ok(fixture)
}

fn control_target(
    path: &Path,
    curve: &KoblitzCurve,
    small: &FrobeniusFactorBase,
    residual: &FrobeniusFactorBase,
) -> Result<BinaryPoint> {
    let fixture = parse_control_bytes(&fs::read(path)?)?;
    if fixture["subgroup_order"] != curve.subgroup_order.to_string() {
        return Err("known-answer subgroup order".into());
    }
    let terms = fixture["terms"].as_array().ok_or("known-answer terms")?;
    let indices = [0usize, 166, 332, 196_212, 196_378];
    if terms.len() != indices.len() {
        return Err("known-answer term count".into());
    }
    let mut sum = BinaryPoint::Infinity;
    for (&index, term) in indices.iter().zip(terms) {
        let selected = if index < 196_212 {
            &small.points[index]
        } else {
            &residual.points[index - 196_212]
        };
        if term["global_point_index"] != index
            || term["domain"] != if index < 196_212 { "small" } else { "residual" }
            || term["column"] != index / 166
            || term["phase"] != 0
            || term["negative"] != false
            || term["coefficient"] != "1"
            || point(&term["point"])? != *selected
        {
            return Err("known-answer selected point differs from imported base".into());
        }
        sum = curve.add(&sum, selected);
    }
    let target = point(&fixture["equation"]["target"])?;
    if target == BinaryPoint::Infinity
        || !curve.curve.is_on_curve(&target)
        || curve.mul(&target, &curve.subgroup_order) != BinaryPoint::Infinity
        || sum != target
    {
        return Err("known-answer group sum mismatch".into());
    }
    Ok(target)
}

fn control_run(args: &[String]) -> Result<()> {
    if args.len() != 10 {
        return Err("usage: control-run SMALL_MANIFEST SMALL_OBJECT ENVELOPE_MANIFEST ENVELOPE_OBJECT FIXTURE BUDGET_SECONDS CONFLICT_BUDGET MAX_MODELS OUTPUT_DIR".into());
    }
    let budget: u64 = args[6].parse()?;
    let conflicts: u64 = args[7].parse()?;
    let models: usize = args[8].parse()?;
    if !(1..=7200).contains(&budget)
        || !(1..=10_000_000).contains(&conflicts)
        || !(1..=1024).contains(&models)
    {
        return Err("invalid bounded control configuration".into());
    }
    cgroup_gate()?;
    let start = Instant::now();
    let output = Path::new(&args[9]);
    fs::create_dir(output)?;
    write_new(
        &output.join("config.json"),
        &json!({
            "schema":"n83.mixed-orbit-known-witness-config/v1",
            "scope":"known_answer_control_only", "study":STUDY,
            "fixture_blake3":blake3::hash(CONTROL_FIXTURE).to_hex().to_string(),
            "small_object_blake3":SMALL_HASH, "envelope_object_blake3":ENVELOPE_HASH,
            "budget_seconds":budget, "conflict_budget":conflicts, "max_models":models,
            "small_summands":3, "residual_summands":2,
            "memory_cgroup_limit_bytes":4294967296u64,
            "memory_cgroup_swap_limit_bytes":0
        }),
    )?;
    let cap_path = output.join("cap.json");
    std::thread::spawn(move || {
        std::thread::sleep(Duration::from_secs(budget));
        let _ = write_new(
            &cap_path,
            &json!({
                "schema":"n83.mixed-orbit-known-witness-cap/v1",
                "status":"UNKNOWN_budget", "budget_seconds":budget,
                "elapsed_process_wall_ms":start.elapsed().as_secs_f64()*1000.0
            }),
        );
        std::process::exit(124);
    });
    let (curve, small, residual) = import_pair(args, &SMALL)?;
    let target = control_target(Path::new(&args[5]), &curve, &small, &residual)?;
    let import_ms = start.elapsed().as_secs_f64() * 1000.0;
    let domains = MixedOrbitDomains::new(&curve, &small, &residual)?;
    let limits = ChainS3Limits {
        max_variables: 150_000,
        max_domain_clauses: 12_000_000,
        max_models: models,
        conflict_budget: conflicts,
    };
    let (candidate, stats) = domains.search(
        1,
        &BigUint::from(0u8),
        &BigUint::from(1u8),
        &target,
        3,
        2,
        limits,
        true,
    );
    let candidate_record = candidate.as_ref().map(|row| -> Result<Value> {
        let mut sum = BinaryPoint::Infinity;
        for &index in &row.small_point_indices {
            sum = curve.add(&sum, small.points.get(index).ok_or("candidate small index")?);
        }
        for point in &row.residual_points {
            sum = curve.add(&sum, point);
        }
        if sum != target || row.coef_a != BigUint::from(0u8)
            || row.coef_b != BigUint::from(1u8)
        {
            return Err("known-answer candidate failed group replay".into());
        }
        Ok(json!({
            "a":row.coef_a.to_string(), "b":row.coef_b.to_string(),
            "small_point_indices":row.small_point_indices,
            "residual_points":row.residual_points.iter().map(point_json).collect::<Result<Vec<_>>>()?
        }))
    }).transpose()?;
    let status = if stats.invalid_models != 0 {
        "PRODUCER_FAILURE_invalid_model"
    } else if candidate_record.is_some() {
        "PASS_known_witness_solver_access"
    } else if stats.exhausted && stats.patterns_attempted == 0 {
        "UNKNOWN_preflight_cap"
    } else if stats.exhausted {
        "UNKNOWN_search_cap"
    } else if stats.refuted {
        "PRODUCER_FAILURE_known_witness_refuted"
    } else {
        "PRODUCER_FAILURE_unclassified_solver_result"
    };
    write_new(
        &output.join("result.json"),
        &json!({
            "schema":"n83.mixed-orbit-known-witness-result/v1",
            "scope":"known_answer_control_only", "status":status,
            "candidate":candidate_record, "search":stats,
            "target":point_json(&target)?,
            "base_import_elapsed_ms_L0":import_ms,
            "worker_wall_ms_L0":start.elapsed().as_secs_f64()*1000.0,
            "memory_cgroup_peak_bytes":cgroup_peak_bytes()?,
            "natural_relations":Value::Null,
            "total_index_calculus_runtime_ms":Value::Null,
            "selected_best_total_runtime":Value::Null
        }),
    )?;
    Ok(())
}

fn run(args: &[String], small_spec: &PinnedBase) -> Result<()> {
    if args.len() != 13 {
        return Err("usage: run SMALL_MANIFEST SMALL_OBJECT ENVELOPE_MANIFEST ENVELOPE_OBJECT CORPUS BUDGET_SECONDS MAX_TRIALS CONFLICT_BUDGET MAX_MODELS SMALL_SUMMANDS RESIDUAL_SUMMANDS OUTPUT_DIR".into());
    }
    let budget: u64 = args[6].parse()?;
    let trials: u64 = args[7].parse()?;
    let conflicts: u64 = args[8].parse()?;
    let models: usize = args[9].parse()?;
    let small_summands: usize = args[10].parse()?;
    let residual_summands: usize = args[11].parse()?;
    if !(1..=7200).contains(&budget)
        || !(1..=1_000_000).contains(&trials)
        || !(1..=10_000_000).contains(&conflicts)
        || !(1..=1024).contains(&models)
        || small_summands < 3
        || !(1..=2).contains(&residual_summands)
        || !(5..=6).contains(&(small_summands + residual_summands))
    {
        return Err("invalid bounded search configuration".into());
    }
    cgroup_gate()?;
    let start = Instant::now();
    let output = Path::new(&args[12]);
    fs::create_dir(output)?;
    write_new(
        &output.join("config.json"),
        &json!({
            "schema":"n83.mixed-orbit-partial-config/v1", "study":STUDY,
            "small_object_blake3":small_spec.object_hash, "envelope_object_blake3":ENVELOPE_HASH,
            "small_manifest_blake3":blake3::hash(small_spec.manifest).to_hex().to_string(),
            "envelope_manifest_blake3":blake3::hash(ENVELOPE_MANIFEST).to_hex().to_string(),
            "source_blake3":blake3::hash(include_bytes!("koblitz_n83_mixed_orbit_worker.rs")).to_hex().to_string(),
            "public_corpus_path":args[5], "public_corpus_canonical_json_blake3":PUBLIC_CORPUS_HASH,
            "budget_seconds":budget,
            "max_trials":trials, "conflict_budget_per_trial":conflicts,
            "max_models_per_trial":models, "small_summands":small_summands,
            "residual_summands":residual_summands, "order_first_pair":true,
            "max_variables":150000, "max_domain_clauses":12000000,
            "graph_limits":{"max_merge_depth":32,"max_row_terms":4096,"max_pivots":1000000},
            "memory_cgroup_limit_bytes":4294967296u64,
            "memory_cgroup_swap_limit_bytes":0
        }),
    )?;
    let cap_path = output.join("cap.json");
    std::thread::spawn(move || {
        std::thread::sleep(Duration::from_secs(budget));
        let _ = write_new(
            &cap_path,
            &json!({
                "schema":"n83.mixed-orbit-partial-cap/v1", "status":"UNKNOWN_budget",
                "budget_seconds":budget,
                "elapsed_process_wall_ms":start.elapsed().as_secs_f64()*1000.0
            }),
        );
        std::process::exit(124);
    });
    let (curve, small, residual) = import_pair(args, small_spec)?;
    let import_ms = start.elapsed().as_secs_f64() * 1000.0;
    let (target, corpus_hash) = public_target(Path::new(&args[5]), &curve)?;
    let domains = MixedOrbitDomains::new(&curve, &small, &residual)?;
    let mut graph =
        WideLargePrimeEliminator::new(&curve, &small, target.clone(), LargePrimeLimits::default())?;
    let limits = ChainS3Limits {
        max_variables: 150_000,
        max_domain_clauses: 12_000_000,
        max_models: models,
        conflict_budget: conflicts,
    };
    let mut events = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(output.join("events.jsonl"))?;
    let mut completed = 0u64;
    let mut attempted = 0u64;
    let mut exhausted = 0u64;
    let mut refuted = 0u64;
    let mut solver_calls = 0u64;
    let mut solver_conflicts = 0u64;
    let mut solver_models = 0u64;
    let mut invalid_models = 0u64;
    let mut group_rejected_models = 0u64;
    let mut preflight_cap = false;
    for trial in 0..trials {
        let a = trial_coefficient(b"a", trial, &curve.subgroup_order);
        let b = trial_coefficient(b"b", trial, &curve.subgroup_order);
        let (candidate, stats) = domains.search(
            trial + 1,
            &a,
            &b,
            &target,
            small_summands,
            residual_summands,
            limits,
            true,
        );
        attempted += 1;
        exhausted += u64::from(stats.exhausted);
        refuted += u64::from(stats.refuted);
        solver_calls += stats.solver_calls as u64;
        solver_conflicts += stats.conflicts;
        solver_models += stats.models as u64;
        invalid_models += stats.invalid_models as u64;
        group_rejected_models += stats.group_rejected_models as u64;
        preflight_cap = stats.exhausted && stats.patterns_attempted == 0;
        let candidate_record = candidate.as_ref().map(|row| -> Result<Value> {
            Ok(json!({
                "source_id":row.source_id, "a":row.coef_a.to_string(),
                "b":row.coef_b.to_string(),
                "small_point_indices":row.small_point_indices,
                "residual_points":row.residual_points.iter().map(point_json).collect::<Result<Vec<_>>>()?
            }))
        }).transpose()?;
        let graph_outcome = if let Some(candidate) = candidate {
            let fed = graph.feed(candidate)?;
            match fed {
                FeedOutcome::PivotStored { merge_depth } => {
                    json!({"kind":"pivot","merge_depth":merge_depth})
                }
                FeedOutcome::UnknownCap {
                    reason,
                    merge_depth,
                } => {
                    json!({"kind":"UNKNOWN_graph_cap","reason":format!("{reason:?}"),"merge_depth":merge_depth})
                }
                FeedOutcome::Full(row) => {
                    completed += 1;
                    json!({"kind":"full","a":row.coef_a.to_string(),"b":row.coef_b.to_string(),
                        "small":row.small.iter().map(|(k,v)| (k.to_string(),v.to_string())).collect::<Vec<_>>(),
                        "source_coeffs":row.source_coeffs.iter().map(|(k,v)| (k.to_string(),v.to_string())).collect::<Vec<_>>(),
                        "merge_depth":row.merge_depth})
                }
            }
        } else {
            Value::Null
        };
        serde_json::to_writer(
            &mut events,
            &json!({
                "trial":trial, "a":a.to_string(), "b":b.to_string(),
                "elapsed_process_wall_ms":start.elapsed().as_secs_f64()*1000.0,
                "search":stats, "candidate":candidate_record, "graph":graph_outcome
            }),
        )?;
        events.write_all(b"\n")?;
        events.flush()?;
        if preflight_cap || invalid_models != 0 {
            break;
        }
    }
    events.sync_all()?;
    let graph_stats = graph.stats();
    let status = if invalid_models != 0 {
        "PRODUCER_FAILURE_invalid_model"
    } else if preflight_cap {
        "UNKNOWN_preflight_cap"
    } else if graph_stats.capped_rows != 0 {
        "UNKNOWN_graph_cap"
    } else if exhausted != 0 {
        "UNKNOWN_search_cap"
    } else {
        "PASS_bounded_partial_search"
    };
    write_new(
        &output.join("result.json"),
        &json!({
            "schema":"n83.mixed-orbit-partial-result/v1", "study":STUDY,
            "status":status, "small_object_blake3":small_spec.object_hash,
            "envelope_object_blake3":ENVELOPE_HASH,
            "public_corpus_canonical_json_blake3":corpus_hash,
            "small_signed_orbits":small.signed_orbits.len(),
            "residual_signed_orbits":residual.signed_orbits.len(),
            "trials_attempted":attempted, "solver_exhausted_trials":exhausted,
            "solver_refuted_trials":refuted, "completed_rows":completed,
            "solver_calls":solver_calls, "solver_conflicts":solver_conflicts,
            "solver_models":solver_models, "invalid_models":invalid_models,
            "group_rejected_models":group_rejected_models,
            "graph_accepted_candidates":graph_stats.accepted_candidates,
            "graph_pivot_rows":graph_stats.pivot_rows,
            "graph_full_from_partials":graph_stats.full_from_partials,
            "graph_capped_rows":graph_stats.capped_rows,
            "graph_pivot_reductions":graph_stats.pivot_reductions,
            "base_import_elapsed_ms_L0":import_ms,
            "worker_wall_ms_L0":start.elapsed().as_secs_f64()*1000.0,
            "memory_cgroup_peak_bytes":cgroup_peak_bytes()?,
            "matrix_rank":Value::Null, "total_index_calculus_runtime_ms":Value::Null,
            "selected_best_total_runtime":Value::Null
        }),
    )?;
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("verify" | "verify-alt") if args.len() == 5 => {
            let small_spec = if args[0] == "verify-alt" { &SMALL_600 } else { &SMALL };
            let (curve, small, residual) = import_pair(&args, small_spec)?;
            println!("{}", json!({
                "schema":"n83.mixed-orbit-import/v1", "status":"PASS",
                "curve_slug":curve.label(), "small_object_blake3":small_spec.object_hash,
                "envelope_object_blake3":ENVELOPE_HASH,
                "small_points":small.points.len(), "small_signed_orbits":small.signed_orbits.len(),
                "residual_points":residual.points.len(), "residual_signed_orbits":residual.signed_orbits.len(),
                "solver_search_executed":false
            }));
            Ok(())
        }
        Some("run") => run(&args, &SMALL),
        Some("run-alt") => run(&args, &SMALL_600),
        Some("verify-control") if args.len() == 6 => {
            let (curve, small, residual) = import_pair(&args, &SMALL)?;
            let target = control_target(Path::new(&args[5]), &curve, &small, &residual)?;
            println!("{}", json!({"schema":"n83.mixed-orbit-known-witness-verify/v1",
                "status":"PASS", "scope":"known_answer_control_only",
                "target":point_json(&target)?, "solver_search_executed":false}));
            Ok(())
        }
        Some("control-run") => control_run(&args),
        _ => Err("usage: verify|verify-alt SMALL_MANIFEST SMALL_OBJECT ENVELOPE_MANIFEST ENVELOPE_OBJECT | verify-control SMALL_MANIFEST SMALL_OBJECT ENVELOPE_MANIFEST ENVELOPE_OBJECT FIXTURE | run|run-alt SMALL_MANIFEST SMALL_OBJECT ENVELOPE_MANIFEST ENVELOPE_OBJECT CORPUS BUDGET_SECONDS MAX_TRIALS CONFLICT_BUDGET MAX_MODELS SMALL_SUMMANDS RESIDUAL_SUMMANDS OUTPUT_DIR | control-run SMALL_MANIFEST SMALL_OBJECT ENVELOPE_MANIFEST ENVELOPE_OBJECT FIXTURE BUDGET_SECONDS CONFLICT_BUDGET MAX_MODELS OUTPUT_DIR".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_witness_fixture_is_pinned_and_distinct_from_natural_corpus() {
        let fixture = parse_control_bytes(CONTROL_FIXTURE).unwrap();
        assert_eq!(fixture["scope"], "known_answer_control_only");
        let mut altered = CONTROL_FIXTURE.to_vec();
        altered[0] ^= 1;
        assert!(parse_control_bytes(&altered).is_err());
    }

    #[test]
    fn pinned_manifests_reject_wrong_bytes_or_archive() {
        assert!(bound_archive(&SMALL, ENVELOPE_MANIFEST, b"").is_err());
        assert!(bound_archive(&SMALL, SMALL_MANIFEST, b"changed archive").is_err());
        assert!(bound_archive(&SMALL_600, SMALL_MANIFEST, b"").is_err());
        assert!(bound_archive(&SMALL_600, SMALL_600_MANIFEST, b"changed archive").is_err());
        assert!(bound_archive(&ENVELOPE, ENVELOPE_MANIFEST, b"changed archive").is_err());
    }

    #[test]
    fn trial_coefficients_are_distinct_and_replayable() {
        let modulus = BigUint::from(97u8);
        let a = trial_coefficient(b"a", 7, &modulus);
        assert_eq!(a, trial_coefficient(b"a", 7, &modulus));
        assert_ne!(a, trial_coefficient(b"a", 8, &modulus));
        assert_ne!(a, trial_coefficient(b"b", 7, &modulus));
    }

    #[test]
    fn matched_public_corpus_digest_is_pinned() {
        let corpus: Value = serde_json::from_slice(include_bytes!(
            "../research/koblitz_n83_factor_base_sweep_20261008/pilot-01/probe-corpus.json"
        ))
        .unwrap();
        assert_eq!(
            blake3::hash(&serde_json::to_vec(&corpus).unwrap())
                .to_hex()
                .as_str(),
            PUBLIC_CORPUS_HASH
        );
    }
}
