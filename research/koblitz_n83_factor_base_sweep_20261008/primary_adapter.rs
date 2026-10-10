//! Replay-bound bridge from the stored primary N83 points to the generic IC base.
//!
//! The library's explicit-orbit constructor currently requires a u64 field.
//! This bridge uses the verified JSONL orbit order directly and checks every
//! point and signed-Frobenius label before building the public structure.

use super::{
    compact_cold, curve, frozen_v2_source_commit, general, load_object, parse_point, point_json,
    point_set_hash, words, write_new_json, Result, SEEDS, STUDY,
};
use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::binary_semaev_chain_sat::{
    chain_s3_max_variables, finite_domain_clause_count, ChainedS3Encoding,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_union_s4_encoding, koblitz_index_calculus_dlp_with_factor_base_and_progress,
    koblitz_point_count, point_key, DecompositionStrategy, FactorBaseDomain, FrobeniusFactorBase,
    KoblitzCurve, KoblitzIcOptions, SatDecompositionOptions,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::time::{Duration, Instant};

fn checked_cgroup_memory_limit(memory_mib: u64) -> Result<u64> {
    let bytes = memory_mib
        .checked_mul(1024 * 1024)
        .filter(|&value| value > 0)
        .ok_or("invalid capacity memory limit")?;
    let actual: u64 = fs::read_to_string("/sys/fs/cgroup/memory.max")?
        .trim()
        .parse()
        .map_err(|_| "capacity worker requires a finite cgroup memory.max")?;
    let swap: u64 = fs::read_to_string("/sys/fs/cgroup/memory.swap.max")?
        .trim()
        .parse()
        .map_err(|_| "capacity worker requires a finite cgroup memory.swap.max")?;
    if actual != bytes || swap != 0 {
        return Err("capacity worker requires the declared hard memory cap and zero swap".into());
    }
    Ok(bytes)
}

fn process_usage() -> Value {
    #[cfg(unix)]
    {
        let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
        if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } == 0 {
            let usage = unsafe { usage.assume_init() };
            let seconds = |value: libc::timeval| {
                value.tv_sec as f64 + value.tv_usec as f64 / 1_000_000.0
            };
            let multiplier = if cfg!(target_os = "macos") { 1 } else { 1024 };
            return json!({
                "cpu_seconds":seconds(usage.ru_utime)+seconds(usage.ru_stime),
                "peak_rss_bytes":(usage.ru_maxrss.max(0) as u64).saturating_mul(multiplier)
            });
        }
    }
    json!({"cpu_seconds":null,"peak_rss_bytes":null})
}

pub(super) struct PrimaryBase {
    pub(super) curve: KoblitzCurve,
    pub(super) factor_base: FrobeniusFactorBase,
}

fn from_record(row: &Value, bytes: &[u8]) -> Result<PrimaryBase> {
    let mut lines = std::str::from_utf8(bytes)?.lines();
    let header: Value = serde_json::from_str(lines.next().ok_or("missing base header")?)?;
    let columns = row["columns"].as_u64().ok_or("column count")? as usize;
    if columns == 0 || header["orbit_columns"].as_u64() != Some(columns as u64) {
        return Err("primary orbit-column count mismatch".into());
    }
    let pinned = curve(0)?;
    if row["a"].as_u64() != Some(0)
        || header["curve_a"].as_u64() != Some(0)
        || header["schema"] != "n83.factor-base/v1"
        || header["study"] != STUDY
        || header["curve_slug"] != "icv1-f2m83-tm6151469093347-debefd74"
        || header["n"] != 83
        || header["policy"] != row["policy"]
        || header["seed"] != row["seed"]
        || header["closure"] != "signed_frobenius"
        || header["modulus_low_terms"] != json!([0, 1, 2, 45])
        || header["subgroup_order"] != pinned.subgroup_order.to_string()
        || header["cofactor"] != pinned.cofactor.to_string()
        || header["frobenius_lambda"] != pinned.lambda.to_string()
        || header["generator"] != point_json(words(pinned.generator()))
        || header["point_count"].as_u64() != Some((166 * columns) as u64)
        || row["points"].as_u64() != Some((166 * columns) as u64)
    {
        return Err("primary curve or manifest binding mismatch".into());
    }
    let reps = header["representatives"]
        .as_array()
        .ok_or("representatives")?;
    if reps.len() != columns {
        return Err("representative count mismatch".into());
    }

    let mut points = Vec::with_capacity(166 * columns);
    let mut point_words = Vec::with_capacity(166 * columns);
    let mut orbits = Vec::with_capacity(2 * columns);
    let mut orbit_of = Vec::with_capacity(166 * columns);
    let mut signed_orbits = Vec::with_capacity(columns);
    let mut signed_orbit_of = Vec::with_capacity(166 * columns);
    let mut x_values = BTreeMap::new();
    let mut seen = HashSet::with_capacity(166 * columns);

    for (column, value) in reps.iter().enumerate() {
        let rep_word = parse_point(value)?.ok_or("identity representative")?;
        let rep = general(Some(rep_word));
        if !pinned.curve.is_on_curve(&rep)
            || pinned.mul(&rep, &pinned.subgroup_order) != BinaryPoint::Infinity
            || pinned.mul(&rep, &pinned.lambda) != pinned.frobenius(&rep)
        {
            return Err("primary representative subgroup or eigenvalue mismatch".into());
        }
        let mut current = rep.clone();
        let mut coefficient = BigUint::one();
        let mut positive_orbit = Vec::with_capacity(83);
        let mut negative_orbit = Vec::with_capacity(83);
        let mut signed_members = Vec::with_capacity(166);
        for phase in 0..83u32 {
            if let BinaryPoint::Affine { x, .. } = &current {
                x_values.entry(x.to_biguint()).or_insert_with(|| x.clone());
            } else {
                return Err("identity inside primary orbit".into());
            }
            for negative in [false, true] {
                let point = if negative {
                    point_neg(&current)
                } else {
                    current.clone()
                };
                let entry: Value = serde_json::from_str(lines.next().ok_or("missing point row")?)?;
                let expected_coefficient = if negative {
                    (&pinned.subgroup_order - &coefficient).to_string()
                } else {
                    coefficient.to_string()
                };
                let index = points.len();
                if entry["point"] != point_json(words(&point))
                    || entry["column"].as_u64() != Some(column as u64)
                    || entry["phase"].as_u64() != Some(phase as u64)
                    || entry["negative"].as_bool() != Some(negative)
                    || entry["coefficient"] != expected_coefficient
                    || !seen.insert(point_key(&point))
                {
                    return Err("primary point, label or distinctness mismatch".into());
                }
                point_words.push(words(&point));
                points.push(point);
                orbit_of.push((2 * column + usize::from(negative), phase));
                signed_orbit_of.push((column, phase, negative));
                signed_members.push(index);
                if negative {
                    negative_orbit.push(index);
                } else {
                    positive_orbit.push(index);
                }
            }
            current = pinned.frobenius(&current);
            coefficient = (&coefficient * &pinned.lambda) % &pinned.subgroup_order;
        }
        if current != rep {
            return Err("primary orbit did not close".into());
        }
        orbits.push(positive_orbit);
        orbits.push(negative_orbit);
        signed_orbits.push(signed_members);
    }
    if lines.next().is_some()
        || points.len() != 166 * columns
        || x_values.len() != 83 * columns
        || row["point_set_blake3"] != point_set_hash(point_words)?
    {
        return Err("primary point count or set hash mismatch".into());
    }

    let curve = KoblitzCurve {
        a: 0,
        n: 83,
        trace: -1,
        group_order: koblitz_point_count(0, 83),
        subgroup_order: pinned.subgroup_order,
        cofactor: pinned.cofactor,
        lambda: pinned.lambda,
        frobenius_is_endomorphism: true,
        curve: pinned.curve,
    };
    let factor_base = FrobeniusFactorBase {
        domain: FactorBaseDomain::ExplicitFrobeniusOrbits {
            representatives: columns,
        },
        ell: 83,
        f_j: 0,
        linearised_exponents: Vec::new(),
        subspace: x_values.into_values().collect(),
        subspace_basis: (0..83)
            .map(|bit| F2mElement::from_bit_positions(&[bit], 83))
            .collect(),
        points,
        orbits,
        orbit_of,
        signed_orbits,
        signed_orbit_of,
    };
    Ok(PrimaryBase { curve, factor_base })
}

fn load_panel_selected(
    root: &Path,
    columns: usize,
    policy: &str,
    seed: u64,
) -> Result<(PrimaryBase, Value, String)> {
    if !["public_x_sequential", "public_x_hash", "public_x_gray_prefix"].contains(&policy)
        || !SEEDS.contains(&seed)
    {
        return Err("primary base policy or seed outside frozen panel".into());
    }
    let manifest_bytes = fs::read(root.join("manifest.json"))?;
    let manifest_hash = blake3::hash(&manifest_bytes).to_hex().to_string();
    let manifest: Value = serde_json::from_slice(&manifest_bytes)?;
    let replay: Value = serde_json::from_slice(&fs::read(root.join("replay.json"))?)?;
    if manifest["schema"] != "n83.factor-base-panel/v1"
        || manifest["study"] != STUDY
        || manifest["status"] != "completed_factor_base_panel"
        || manifest["completed_base_count"] != 54
        || replay["schema"] != "n83.factor-base-replay/v1"
        || replay["status"] != "PASS"
        || replay["panel_manifest_blake3"] != manifest_hash
    {
        return Err("panel or replay receipt is not complete and bound".into());
    }
    let rows = manifest["bases"].as_array().ok_or("panel base rows")?;
    if rows.len() != 54 {
        return Err("panel base count mismatch".into());
    }
    let mut rows = rows.iter().filter(|row| {
        row["a"].as_u64() == Some(0)
            && row["policy"] == policy
            && row["seed"].as_u64() == Some(seed)
            && row["columns"].as_u64() == Some(columns as u64)
    });
    let row = rows.next().ok_or("requested primary base absent")?;
    if rows.next().is_some() {
        return Err("ambiguous primary base selection".into());
    }
    let checks = replay["checks"].as_array().ok_or("replay checks")?;
    if !checks
        .iter()
        .any(|check| check["object"] == row["object"] && check["status"] == "PASS")
    {
        return Err("selected object lacks a successful replay".into());
    }
    let base = from_record(row, &load_object(root, row)?)?;
    Ok((base, row.clone(), manifest_hash))
}

fn load_panel(root: &Path, columns: usize) -> Result<(PrimaryBase, Value, String)> {
    load_panel_selected(root, columns, "public_x_hash", SEEDS[0])
}

pub(super) fn check_panel(root: &Path, columns: usize) -> Result<Value> {
    let (base, row, _) = load_panel(root, columns)?;
    Ok(json!({
        "schema": "n83.primary-factor-base-adapter/v1",
        "study": STUDY,
        "status": "PASS",
        "curve_a": 0,
        "subgroup_order": base.curve.subgroup_order.to_string(),
        "object": row["object"],
        "point_set_blake3": row["point_set_blake3"],
        "orbit_columns": columns,
        "points": base.factor_base.points.len(),
        "frobenius_orbits": base.factor_base.orbits.len(),
        "signed_frobenius_orbits": base.factor_base.signed_orbits.len(),
        "solver_stage_executed": false,
        "total_index_calculus_runtime_ms": Value::Null,
        "selected_best_total_runtime": Value::Null
    }))
}

/// Convert a checked four-summand witness into exact signed-orbit
/// coefficients. The retained primary base consists entirely of order-r
/// points, so replaying the row against its representatives must reproduce
/// the original target without a cofactor projection.
fn checked_primary_relation_row(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    indices: &[usize; 4],
    target: &BinaryPoint,
) -> Result<Vec<BigUint>> {
    let r = &kc.subgroup_order;
    if r <= &BigUint::one() || fb.signed_orbit_of.len() != fb.points.len() {
        return Err("invalid primary relation domain".into());
    }
    let mut row = vec![BigUint::zero(); fb.unknowns()];
    let mut point_sum = BinaryPoint::Infinity;
    for &index in indices {
        let point = fb
            .points
            .get(index)
            .ok_or("relation point index out of range")?;
        let &(column, phase, negative) = fb
            .signed_orbit_of
            .get(index)
            .ok_or("relation point label absent")?;
        let rep_index = *fb
            .signed_orbits
            .get(column)
            .and_then(|orbit| orbit.first())
            .ok_or("relation orbit representative absent")?;
        let representative = fb
            .points
            .get(rep_index)
            .ok_or("relation representative index out of range")?;
        if phase >= kc.n || kc.mul(representative, r) != BinaryPoint::Infinity {
            return Err("relation representative outside the primary subgroup".into());
        }
        let positive = kc.lambda.modpow(&BigUint::from(phase), r);
        let coefficient = if negative && !positive.is_zero() {
            r - &positive
        } else {
            positive
        };
        if kc.mul(representative, &coefficient) != *point {
            return Err("relation point and signed-Frobenius label disagree".into());
        }
        row[column] = (&row[column] + coefficient) % r;
        point_sum = kc.add(&point_sum, point);
    }
    if point_sum != *target {
        return Err("relation point sum differs from target".into());
    }
    let mut row_sum = BinaryPoint::Infinity;
    for (column, coefficient) in row.iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let rep_index = fb.signed_orbits[column][0];
        row_sum = kc.add(&row_sum, &kc.mul(&fb.points[rep_index], coefficient));
    }
    if row_sum != *target {
        return Err("relation row failed independent group replay".into());
    }
    Ok(row)
}

/// The container sees only a supervisor-made source snapshot, never the
/// checkout's Git metadata or unrelated files. The outer supervisor must
/// independently attest a clean checkout at this commit and retain the
/// snapshot hashes; this worker checks that the three files it executes
/// match its compiled bytes.
fn checked_primary_probe_source(root: &Path) -> Result<()> {
    for (path, compiled) in [
        (
            "examples/koblitz_n83_factor_base_export.rs",
            include_bytes!("../../examples/koblitz_n83_factor_base_export.rs").as_slice(),
        ),
        (
            "research/koblitz_n83_factor_base_sweep_20261008/primary_adapter.rs",
            include_bytes!("primary_adapter.rs").as_slice(),
        ),
        (
            "research/koblitz_n83_factor_base_sweep_20261008/compact_cold.rs",
            include_bytes!("compact_cold.rs").as_slice(),
        ),
    ] {
        if blake3::hash(&fs::read(root.join(path))?) != blake3::hash(compiled) {
            return Err("compiled S3 probe source differs from frozen source".into());
        }
    }
    Ok(())
}

fn frozen_primary_probe_source() -> Result<(String, &'static str)> {
    if let Some(path) = std::env::var_os("ICV1_FROZEN_SOURCE_DIR") {
        let root = std::path::PathBuf::from(path);
        let attestation: Value =
            serde_json::from_slice(&fs::read(root.join("attestation.json"))?)?;
        let commit = attestation["source_commit"]
            .as_str()
            .ok_or("source snapshot lacks a commit")?;
        if attestation["schema"] != "n83.primary-probe-source-attestation/v1"
            || attestation["status_clean"] != true
            || commit.len() != 40
            || !commit.bytes().all(|byte| byte.is_ascii_hexdigit())
        {
            return Err("invalid supervised source attestation".into());
        }
        checked_primary_probe_source(&root)?;
        Ok((commit.to_owned(), "supervised_snapshot"))
    } else {
        let commit = frozen_v2_source_commit()?;
        checked_primary_probe_source(Path::new(env!("CARGO_MANIFEST_DIR")))?;
        Ok((commit, "clean_checkout"))
    }
}

/// Pin the retained-base importer and both layers of the chained-S3 model.
/// The outer guard freezes these sources at a clean commit before launching
/// the Linux worker. The worker compares the snapshot to its compiled bytes.
fn frozen_chain_source() -> Result<(String, &'static str)> {
    let files: [(&str, &[u8]); 4] = [
        ("examples/koblitz_n83_factor_base_export.rs", include_bytes!("../../examples/koblitz_n83_factor_base_export.rs")),
        ("research/koblitz_n83_factor_base_sweep_20261008/primary_adapter.rs", include_bytes!("primary_adapter.rs")),
        ("src/cryptanalysis/binary_semaev_chain_sat.rs", include_bytes!("../../src/cryptanalysis/binary_semaev_chain_sat.rs")),
        ("src/cryptanalysis/koblitz_index_calculus.rs", include_bytes!("../../src/cryptanalysis/koblitz_index_calculus.rs")),
    ];
    let (root, commit, mode) = if let Some(path) = std::env::var_os("ICV1_FROZEN_SOURCE_DIR") {
        let root = std::path::PathBuf::from(path);
        let attestation: Value = serde_json::from_slice(&fs::read(root.join("attestation.json"))?)?;
        let commit = attestation["source_commit"].as_str().ok_or("chain source snapshot lacks a commit")?;
        if attestation["schema"] != "n83.chain-s3-source-attestation/v1"
            || attestation["status_clean"] != true
            || commit.len() != 40
            || !commit.bytes().all(|byte| byte.is_ascii_hexdigit()) {
            return Err("invalid chained-S3 source attestation".into());
        }
        (root, commit.to_owned(), "supervised_snapshot")
    } else {
        (std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")), frozen_v2_source_commit()?, "clean_checkout")
    };
    for (path, compiled) in files {
        if blake3::hash(&fs::read(root.join(path))?) != blake3::hash(compiled) {
            return Err("compiled chained-S3 source differs from frozen source".into());
        }
    }
    Ok((commit, mode))
}

/// Construct one identity-pattern model from an exact replayed primary base.
/// This is a capacity diagnostic only: no SAT search or relation claim occurs.
/// A hard cgroup memory cap is mandatory here; the caller supplies the wall cap.
#[allow(clippy::too_many_arguments)]
pub(super) fn chain_capacity_cli(
    root: &Path, columns: usize, policy: &str, seed: u64, m: usize,
    identity_mask: u32, max_variables: u32, max_domain_clauses: usize,
    memory_mib: u64, output: &Path,
) -> Result<Value> {
    if ![64, 256, 600].contains(&columns) || !(5..=6).contains(&m)
        || identity_mask >= (1u32 << (m - 2))
        || max_variables == 0 || max_domain_clauses == 0 {
        return Err("chained-S3 capacity parameters outside supported domain".into());
    }
    let memory_limit_bytes = checked_cgroup_memory_limit(memory_mib)?;
    let (source_commit, source_attestation_mode) = frozen_chain_source()?;
    let start = Instant::now();
    let (base, row, manifest_hash) = load_panel_selected(root, columns, policy, seed)?;
    let base_import_ms = start.elapsed().as_secs_f64() * 1000.0;
    let (target, corpus_hash) = public_target(root, &base.curve)?;
    let target_validation_ms = start.elapsed().as_secs_f64() * 1000.0 - base_import_ms;
    let mut codes: Vec<BigUint> = base.factor_base.subspace.iter().map(F2mElement::to_biguint).collect();
    codes.sort_unstable();
    codes.dedup();
    if codes.len() != 83 * columns {
        return Err("primary coordinate domain does not contain 83 distinct x values per column".into());
    }
    let required_variables = chain_s3_max_variables(83, m).ok_or("unsupported chained-S3 arity")?;
    let required_domain_clauses = finite_domain_clause_count(&codes, 83, m).ok_or("invalid chained-S3 domain")?;
    let model_start = Instant::now();
    let mut enc = None;
    let status = if required_variables > u64::from(max_variables) {
        "UNKNOWN_variable_cap"
    } else if required_domain_clauses > max_domain_clauses {
        "UNKNOWN_domain_clause_cap"
    } else {
        let target_x = match &target { BinaryPoint::Affine { x, .. } => Some(x), BinaryPoint::Infinity => None };
        let mut built = ChainedS3Encoding::build(
            m, &base.curve.curve.irreducible, &base.curve.curve.b,
            target_x, identity_mask, max_variables,
        ).map_err(|_| "chained-S3 model build failed after preflight")?;
        built.constrain_summands(&codes, max_domain_clauses)
            .map_err(|_| "chained-S3 finite-domain install failed after preflight")?;
        enc = Some(built);
        "PASS_model_construction_only"
    };
    let model_construction_ms = model_start.elapsed().as_secs_f64() * 1000.0;
    let cgroup_peak_bytes: u64 = fs::read_to_string("/sys/fs/cgroup/memory.peak")?.trim().parse()?;
    let mut receipt = json!({
        "schema":"n83.chain-s3-capacity-worker/v1", "study":STUDY, "status":status,
        "curve_a":0, "fixture":0, "orbit_columns":columns, "policy":policy, "seed":seed,
        "summands":m, "identity_mask":identity_mask, "object":row["object"],
        "point_set_blake3":row["point_set_blake3"], "panel_manifest_blake3":manifest_hash,
        "public_corpus_canonical_json_blake3":corpus_hash,
        "source_commit":source_commit, "source_attestation_mode":source_attestation_mode,
        "source_exporter_blake3":blake3::hash(include_bytes!("../../examples/koblitz_n83_factor_base_export.rs")).to_hex().to_string(),
        "source_adapter_blake3":blake3::hash(include_bytes!("primary_adapter.rs")).to_hex().to_string(),
        "source_chain_blake3":blake3::hash(include_bytes!("../../src/cryptanalysis/binary_semaev_chain_sat.rs")).to_hex().to_string(),
        "source_index_calculus_blake3":blake3::hash(include_bytes!("../../src/cryptanalysis/koblitz_index_calculus.rs")).to_hex().to_string()
    });
    let measurements = json!({
        "memory_cgroup_limit_bytes":memory_limit_bytes, "memory_cgroup_swap_limit_bytes":0,
        "memory_cgroup_peak_bytes":cgroup_peak_bytes,
        "max_variables":max_variables, "max_domain_clauses":max_domain_clauses,
        "required_max_variables":required_variables, "required_domain_clauses":required_domain_clauses,
        "legal_x_coordinates":codes.len(),
        "sat_variables":enc.as_ref().map(|e| e.solver.n_vars()),
        "sat_clauses":enc.as_ref().map(|e| e.solver.n_clauses()),
        "sat_xor_rows":enc.as_ref().map(|e| e.solver.n_xors()),
        "sat_and_gates":enc.as_ref().map(|e| e.n_and_gates),
        "sat_s3_nodes":enc.as_ref().map(|e| e.n_s3_nodes),
        "base_import_ms":base_import_ms, "target_validation_ms":target_validation_ms,
        "model_construction_ms":model_construction_ms,
        "process_wall_ms":start.elapsed().as_secs_f64()*1000.0, "process_usage":process_usage(),
        "solver_search_executed":false, "relation_stage_executed":false,
        "rank_stage_executed":false, "total_index_calculus_runtime_ms":Value::Null,
        "selected_best_total_runtime":Value::Null
    });
    receipt.as_object_mut().ok_or("invalid chained-S3 receipt")?
        .extend(measurements.as_object().ok_or("invalid chained-S3 measurements")?.clone());
    write_new_json(output, &receipt)?;
    Ok(receipt)
}

/// Probe the published primary point with a replayed retained base through
/// the point-only wide compact index. A hit is also converted to a full-width
/// relation row and independently replayed. This is one relation diagnostic,
/// not a rank or total-runtime result. A Linux cgroup with a hard memory limit and
/// zero swap is mandatory; the caller must also impose a process-wall cap.
pub(super) fn s3_probe_cli(
    root: &Path,
    columns: usize,
    policy: &str,
    seed: u64,
    unordered_pairs: bool,
    max_candidate_states: usize,
    memory_mib: u64,
    output: &Path,
) -> Result<Value> {
    if ![64, 256, 600].contains(&columns) {
        return Err("S3 probe supports retained primary K=64/256/600 bases".into());
    }
    let memory_limit_bytes = checked_cgroup_memory_limit(memory_mib)?;
    let (source_commit, source_attestation_mode) = frozen_primary_probe_source()?;
    let started = Instant::now();
    let (base, row, manifest_hash) = load_panel_selected(root, columns, policy, seed)?;
    let base_import_ms = started.elapsed().as_secs_f64() * 1000.0;
    let (target, corpus_hash) = public_target(root, &base.curve)?;
    let target_validation_ms = started.elapsed().as_secs_f64() * 1000.0 - base_import_ms;
    let representatives: Vec<_> = base
        .factor_base
        .signed_orbits
        .iter()
        .map(|orbit| base.factor_base.points[orbit[0]].clone())
        .collect();
    let probe_started = Instant::now();
    let probe = compact_cold::probe_checked_four_sum(
        &base.curve,
        &base.factor_base.points,
        &representatives,
        &target,
        unordered_pairs,
        max_candidate_states,
    )?;
    let probe_ms = probe_started.elapsed().as_secs_f64() * 1000.0;
    let row_started = Instant::now();
    let relation_row = if probe["status"] == "HIT" {
        let indices: [usize; 4] = probe["point_indices"]
            .as_array()
            .ok_or("compact hit lacks point indices")?
            .iter()
            .map(|value| {
                value
                    .as_u64()
                    .and_then(|index| usize::try_from(index).ok())
                    .ok_or("compact hit has invalid point index")
            })
            .collect::<std::result::Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| "compact hit is not a four-summand witness")?;
        let row = checked_primary_relation_row(&base.curve, &base.factor_base, &indices, &target)?;
        let dense_decimal: Vec<_> = row.iter().map(ToString::to_string).collect();
        let nonzero: Vec<_> = row
            .iter()
            .enumerate()
            .filter(|(_, coefficient)| !coefficient.is_zero())
            .map(|(column, coefficient)| {
                json!({"column":column,"coefficient_decimal":coefficient.to_string()})
            })
            .collect();
        Some(json!({
            "schema":"n83.primary-wide-relation-row/v1",
            "modulus_decimal":base.curve.subgroup_order.to_string(),
            "orbit_columns":columns,
            "nonzero":nonzero,
            "dense_decimal_blake3":blake3::hash(&serde_json::to_vec(&dense_decimal)?).to_hex().to_string(),
            "group_verified":true
        }))
    } else if probe["status"] == "MISS" {
        None
    } else {
        return Err("compact probe returned an unknown status".into());
    };
    let relation_row_stage_executed = relation_row.is_some();
    let relation_row_ms = row_started.elapsed().as_secs_f64() * 1000.0;
    let cgroup_peak_bytes: u64 = fs::read_to_string("/sys/fs/cgroup/memory.peak")?
        .trim()
        .parse()
        .map_err(|_| "S3 probe requires a numeric cgroup memory.peak")?;
    let receipt = json!({
        "schema":"n83.primary-compact-s3-probe/v1",
        "study":STUDY,
        "status":probe["status"],
        "curve_a":0,
        "fixture":0,
        "orbit_columns":columns,
        "policy":policy,
        "seed":seed,
        "object":row["object"],
        "point_set_blake3":row["point_set_blake3"],
        "panel_manifest_blake3":manifest_hash,
        "public_corpus_canonical_json_blake3":corpus_hash,
        "source_commit":source_commit,
        "source_attestation_mode":source_attestation_mode,
        "source_exporter_blake3":blake3::hash(include_bytes!("../../examples/koblitz_n83_factor_base_export.rs")).to_hex().to_string(),
        "source_adapter_blake3":blake3::hash(include_bytes!("primary_adapter.rs")).to_hex().to_string(),
        "source_compact_blake3":blake3::hash(include_bytes!("compact_cold.rs")).to_hex().to_string(),
        "memory_cgroup_limit_bytes":memory_limit_bytes,
        "memory_cgroup_swap_limit_bytes":0,
        "memory_cgroup_peak_bytes":cgroup_peak_bytes,
        "max_candidate_states":max_candidate_states,
        "base_import_ms":base_import_ms,
        "target_validation_ms":target_validation_ms,
        "point_probe_ms":probe_ms,
        "relation_row_ms":relation_row_ms,
        "process_wall_ms":started.elapsed().as_secs_f64()*1000.0,
        "process_usage":process_usage(),
        "probe":probe,
        "relation_row":relation_row,
        "relation_row_stage_executed":relation_row_stage_executed,
        "rank_stage_executed":false,
        "column_log_verification":false,
        "total_index_calculus_runtime_ms":Value::Null,
        "selected_best_total_runtime":Value::Null
    });
    write_new_json(output, &receipt)?;
    Ok(receipt)
}

/// Construct the exact K=64/256/600 finite-domain factored S4 model and
/// stop before search. The caller must impose and independently verify a
/// hard process-wall limit; this worker refuses an absent cgroup memory cap.
pub(super) fn build_capacity_cli(
    root: &Path,
    columns: usize,
    memory_mib: u64,
    output: &Path,
) -> Result<Value> {
    if ![64, 256, 600].contains(&columns) {
        return Err("capacity gate supports only retained primary K=64/256/600 bases".into());
    }
    let memory_limit_bytes = checked_cgroup_memory_limit(memory_mib)?;
    let start = Instant::now();
    let (base, row, manifest_hash) = load_panel(root, columns)?;
    let base_import_ms = start.elapsed().as_secs_f64() * 1000.0;
    let (target, corpus_hash) = public_target(root, &base.curve)?;
    let target_validation_ms = start.elapsed().as_secs_f64() * 1000.0 - base_import_ms;
    let model_start = Instant::now();
    let enc = build_union_s4_encoding(
        &base.curve,
        &base.factor_base,
        &target,
        SatDecompositionOptions {
            factored_s4: true,
            ..Default::default()
        },
    )?;
    let model_construction_ms = model_start.elapsed().as_secs_f64() * 1000.0;
    let cgroup_peak_bytes: u64 = fs::read_to_string("/sys/fs/cgroup/memory.peak")?
        .trim()
        .parse()
        .map_err(|_| "capacity worker requires a numeric cgroup memory.peak")?;
    let receipt = json!({
        "schema":"n83.factored-s4-capacity-worker/v1",
        "study":STUDY,
        "status":"PASS_model_construction_only",
        "curve_a":0,
        "fixture":0,
        "orbit_columns":columns,
        "object":row["object"],
        "point_set_blake3":row["point_set_blake3"],
        "panel_manifest_blake3":manifest_hash,
        "public_corpus_canonical_json_blake3":corpus_hash,
        "source_adapter_blake3":blake3::hash(include_bytes!("primary_adapter.rs")).to_hex().to_string(),
        "source_union_builder_blake3":blake3::hash(include_bytes!("../../src/cryptanalysis/koblitz_index_calculus.rs")).to_hex().to_string(),
        "source_factored_encoder_blake3":blake3::hash(include_bytes!("../../src/cryptanalysis/semaev_sat.rs")).to_hex().to_string(),
        "memory_cgroup_limit_bytes":memory_limit_bytes,
        "memory_cgroup_swap_limit_bytes":0,
        "memory_cgroup_peak_bytes":cgroup_peak_bytes,
        "base_points":base.factor_base.points.len(),
        "legal_x_coordinates":base.factor_base.subspace.len(),
        "sat_variables":enc.solver.n_vars(),
        "sat_clauses":enc.solver.n_clauses(),
        "sat_xor_rows":enc.solver.n_xors(),
        "sat_x_variables":enc.n_x_vars,
        "sat_e_variables":enc.n_e_vars,
        "sat_aux_variables":enc.n_aux_vars,
        "trivially_unsat":enc.trivially_unsat,
        "base_import_ms":base_import_ms,
        "target_validation_ms":target_validation_ms,
        "model_construction_ms":model_construction_ms,
        "process_wall_ms":start.elapsed().as_secs_f64()*1000.0,
        "process_usage":process_usage(),
        "solver_search_executed":false,
        "model_lifting_executed":false,
        "total_index_calculus_runtime_ms":Value::Null,
        "selected_best_total_runtime":Value::Null
    });
    write_new_json(output, &receipt)?;
    Ok(receipt)
}

fn public_target(root: &Path, kc: &KoblitzCurve) -> Result<(BinaryPoint, String)> {
    let corpus: Value = serde_json::from_slice(&fs::read(root.join("probe-corpus.json"))?)?;
    if corpus["schema"] != "n83.public-probe-corpus/v1" || corpus["study"] != STUDY {
        return Err("public target corpus mismatch".into());
    }
    let targets = corpus["targets"].as_array().ok_or("public targets")?;
    let mut selected = targets
        .iter()
        .filter(|row| row["a"] == 0 && row["fixture"] == 0);
    let row = selected
        .next()
        .ok_or("primary public fixture zero absent")?;
    if selected.next().is_some() {
        return Err("ambiguous primary public fixture".into());
    }
    let encoded = parse_point(&row["point"])?.ok_or("identity public target")?;
    let target = general(Some(encoded));
    if !kc.curve.is_on_curve(&target)
        || kc.mul(&target, &kc.subgroup_order) != BinaryPoint::Infinity
    {
        return Err("primary public target fails group validation".into());
    }
    let corpus_hash = blake3::hash(&serde_json::to_vec(&corpus)?)
        .to_hex()
        .to_string();
    Ok((target, corpus_hash))
}

fn candidate_matches_validation(
    root: &Path,
    corpus_hash: &str,
    candidate: &BigUint,
) -> Result<bool> {
    let validation: Value = serde_json::from_slice(&fs::read(root.join("probe-validation.json"))?)?;
    if validation["schema"] != "n83.public-probe-validation/v1"
        || validation["corpus_canonical_json_blake3"] != corpus_hash
        || validation["oracle_receives_known_answer_scalars"] != false
    {
        return Err("public validation sidecar mismatch".into());
    }
    let entries = validation["validation"]
        .as_array()
        .ok_or("validation entries")?;
    let mut selected = entries
        .iter()
        .filter(|row| row["a"] == 0 && row["fixture"] == 0);
    let row = selected.next().ok_or("primary validation fixture absent")?;
    if selected.next().is_some() {
        return Err("ambiguous primary validation fixture".into());
    }
    let expected: BigUint = row["known_answer_scalar"]
        .as_str()
        .ok_or("known-answer encoding")?
        .parse()?;
    Ok(*candidate == expected)
}

/// Run one bounded primary workload. The public target is the solver's only
/// target input; its known-answer sidecar is opened only after a candidate log.
#[allow(clippy::too_many_arguments)]
pub(super) fn run_cli(
    root: &Path,
    columns: usize,
    m: usize,
    strategy_name: &str,
    max_trials: usize,
    budget_seconds: u64,
    run_dir: &Path,
) -> Result<Value> {
    if ![64, 256, 600].contains(&columns) || !(2..=6).contains(&m) {
        return Err("unsupported primary K or summand count".into());
    }
    let strategy = match strategy_name {
        "enumerate" => DecompositionStrategy::Enumerate,
        "sat-m3" => return Err("primary SAT S4 model replay and resource-capacity gate is incomplete for F2^83".into()),
        _ => return Err("primary strategy must be enumerate".into()),
    };
    if budget_seconds == 0 || budget_seconds > 86_400 {
        return Err("primary wall cap must be 1..86400 seconds".into());
    }
    let start = Instant::now();
    fs::create_dir(run_dir)?;
    write_new_json(
        &run_dir.join("config.json"),
        &json!({
            "schema":"n83.primary-cold-config/v1", "study":STUDY,
            "panel_dir":root.to_string_lossy(),
            "curve_a":0, "fixture":0, "orbit_columns":columns,
            "summands":m, "strategy":strategy_name, "max_trials":max_trials,
            "budget_seconds":budget_seconds, "allow_direct_relation":false,
            "source_adapter_blake3":blake3::hash(include_bytes!("primary_adapter.rs")).to_hex().to_string()
        }),
    )?;
    let cap_path = run_dir.join("cap.json");
    let cap_strategy = strategy_name.to_owned();
    let deadline = start + Duration::from_secs(budget_seconds);
    std::thread::spawn(move || {
        std::thread::sleep(deadline.saturating_duration_since(Instant::now()));
        let cap = json!({
            "schema":"n83.primary-cold-cap/v1", "status":"UNKNOWN_budget",
            "orbit_columns":columns, "summands":m, "strategy":cap_strategy,
            "max_trials":max_trials, "budget_seconds":budget_seconds,
            "elapsed_process_wall_ms":start.elapsed().as_secs_f64()*1000.0,
            "selected_best_total_runtime":Value::Null
        });
        if write_new_json(&cap_path, &cap).is_err() {
            std::process::exit(125);
        }
        std::process::exit(124);
    });

    let import_start = Instant::now();
    let (base, row, manifest_hash) = load_panel(root, columns)?;
    let base_import_ms = import_start.elapsed().as_secs_f64() * 1000.0;
    let target_start = Instant::now();
    let (target, corpus_hash) = public_target(root, &base.curve)?;
    let target_validation_ms = target_start.elapsed().as_secs_f64() * 1000.0;

    let mut options = KoblitzIcOptions::default();
    options.m = m;
    options.strategy = strategy;
    options.max_trials = max_trials;
    options.seed = 2026100901;
    options.allow_direct_relation = false;
    let mut events = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(run_dir.join("events.jsonl"))?;
    let solver_start = Instant::now();
    let mut event_write_failed = false;
    let report = koblitz_index_calculus_dlp_with_factor_base_and_progress(
        &base.curve,
        &target,
        &base.factor_base,
        &options,
        &mut |event| {
            let line = json!({
                "elapsed_process_wall_ms":start.elapsed().as_secs_f64()*1000.0,
                "event":format!("{event:?}")
            });
            if writeln!(events, "{line}")
                .and_then(|_| events.flush())
                .is_err()
            {
                event_write_failed = true;
            }
        },
    );
    let solver_ms = solver_start.elapsed().as_secs_f64() * 1000.0;
    if event_write_failed {
        return Err("primary progress-event write failed".into());
    }
    let Some(report) = report else {
        let summary = json!({
            "schema":"n83.primary-cold-result/v1", "status":"PRODUCER_FAILURE",
            "reason":"generic solver rejected the primary base or strategy",
            "solver_stage_executed":false, "selected_best_total_runtime":Value::Null
        });
        write_new_json(&run_dir.join("summary.json"), &summary)?;
        return Err("generic primary solver rejected the imported base or strategy".into());
    };
    let verification_start = Instant::now();
    let candidate_valid = match &report.log {
        Some(candidate) => {
            base.curve.mul(base.curve.generator(), candidate) == target
                && candidate_matches_validation(root, &corpus_hash, candidate).unwrap_or(false)
        }
        None => false,
    };
    let post_solver_validation_ms = verification_start.elapsed().as_secs_f64() * 1000.0;
    let status = if report.inconsistent_relations != 0
        || report.verification_failures != 0
        || report.sat_invalid_models != 0
        || report.direct_relation
        || (report.log.is_some() && !candidate_valid)
    {
        "PRODUCER_FAILURE"
    } else if candidate_valid {
        "PASS_verified_target_only"
    } else if !report.m_cofactor_admissible {
        "INADMISSIBLE_cofactor_class"
    } else if max_trials == 0 {
        "PREFLIGHT_ONLY"
    } else {
        "UNKNOWN_trial_cap"
    };
    let summary = json!({
        "schema":"n83.primary-cold-result/v1", "study":STUDY,
        "status":status, "curve_a":0, "fixture":0,
        "object":row["object"], "point_set_blake3":row["point_set_blake3"],
        "panel_manifest_blake3":manifest_hash,
        "public_corpus_canonical_json_blake3":corpus_hash,
        "orbit_columns":columns, "summands":m, "strategy":strategy_name,
        "max_trials":max_trials, "sampler_seed":options.seed,
        "budget_seconds":budget_seconds,
        "base_import_ms":base_import_ms, "target_validation_ms":target_validation_ms,
        "solver_ms":solver_ms, "post_solver_validation_ms":post_solver_validation_ms,
        "pre_summary_process_wall_ms":start.elapsed().as_secs_f64()*1000.0,
        "solver_stage_executed":report.trials>0,
        "report":{
            "factor_base_size":report.factor_base_size,
            "orbit_count":report.orbit_count,
            "trials":report.trials, "relations":report.relations,
            "independent_relations":report.independent_relations,
            "dependent_relations":report.dependent_relations,
            "inconsistent_relations":report.inconsistent_relations,
            "verification_failures":report.verification_failures,
            "sat_unknowns":report.sat_unknowns,
            "sat_invalid_models":report.sat_invalid_models,
            "sat_group_rejected_models":report.sat_group_rejected_models,
            "linear_solve_attempts":report.linear_solve_attempts,
            "relation_collection_ns":report.relation_collection_ns.to_string(),
            "linear_algebra_ns":report.linear_algebra_ns.to_string(),
            "m_cofactor_admissible":report.m_cofactor_admissible,
            "direct_relations_skipped":report.direct_relations_skipped
        },
        "verified_log":if candidate_valid { report.log.as_ref().map(ToString::to_string) } else { None },
        "column_log_verification":false,
        "total_index_calculus_runtime_ms":Value::Null,
        "selected_best_total_runtime":Value::Null
    });
    write_new_json(&run_dir.join("summary.json"), &summary)?;
    if status == "PRODUCER_FAILURE" {
        Err("primary solver result failed verification; inspect summary.json".into())
    } else {
        Ok(summary)
    }
}

#[cfg(test)]
mod tests {
    use super::super::construct;
    use super::*;
    use crypto_lib::cryptanalysis::koblitz_index_calculus::{
        enumerate_decompose, koblitz_index_calculus_dlp_with_factor_base, KoblitzIcOptions,
    };
    use crypto_lib::cryptanalysis::koblitz_relation_solver::WideRankTracker;
    use std::collections::HashMap;
    use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

    #[test]
    fn supervised_source_snapshot_rejects_a_changed_compiled_file() {
        let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
        let snapshot = std::env::temp_dir().join(format!("n83-probe-source-{nonce}"));
        let checkout = Path::new(env!("CARGO_MANIFEST_DIR"));
        for relative in [
            "examples/koblitz_n83_factor_base_export.rs",
            "research/koblitz_n83_factor_base_sweep_20261008/primary_adapter.rs",
            "research/koblitz_n83_factor_base_sweep_20261008/compact_cold.rs",
        ] {
            let destination = snapshot.join(relative);
            fs::create_dir_all(destination.parent().unwrap()).unwrap();
            fs::copy(checkout.join(relative), destination).unwrap();
        }
        assert!(checked_primary_probe_source(&snapshot).is_ok());
        fs::write(
            snapshot.join("research/koblitz_n83_factor_base_sweep_20261008/primary_adapter.rs"),
            b"changed source",
        )
        .unwrap();
        assert!(checked_primary_probe_source(&snapshot).is_err());
        fs::remove_dir_all(snapshot).unwrap();
    }

    #[test]
    fn small_primary_export_maps_to_generic_orbits_and_zero_trial_preflight() {
        let (bytes, row) = construct(
            0,
            "public_x_hash",
            2,
            17,
            Instant::now() + Duration::from_secs(60),
        )
        .unwrap();
        let base = from_record(&row, &bytes).unwrap();
        assert_eq!(base.factor_base.points.len(), 332);
        assert_eq!(base.factor_base.orbits.len(), 4);
        assert_eq!(base.factor_base.signed_orbits.len(), 2);
        assert_eq!(base.factor_base.subspace.len(), 166);
        for (index, point) in base.factor_base.points.iter().enumerate() {
            let (column, phase, negative) = base.factor_base.signed_orbit_of[index];
            assert_eq!(column, index / 166);
            assert_eq!(phase as usize, (index % 166) / 2);
            assert_eq!(negative, index % 2 == 1);
            assert_eq!(
                base.factor_base.orbit_of[index].0,
                2 * column + usize::from(negative)
            );
            assert_eq!(
                point,
                &base.factor_base.points[base.factor_base.signed_orbits[column][index % 166]]
            );
        }
        let kc = &base.curve;
        let fb = &base.factor_base;
        let indices = [0, 2, 165, 168];
        let target = indices.iter().fold(BinaryPoint::Infinity, |sum, &index| {
            kc.add(&sum, &fb.points[index])
        });
        let relation = checked_primary_relation_row(kc, fb, &indices, &target).unwrap();
        let lambda_82 = kc.lambda.modpow(&BigUint::from(82u32), &kc.subgroup_order);
        assert_eq!(
            relation[0],
            (BigUint::one() + &kc.lambda + &kc.subgroup_order - lambda_82) % &kc.subgroup_order
        );
        assert_eq!(relation[1], kc.lambda);
        assert!(relation.iter().any(|coefficient| coefficient.bits() > 64));
        assert!(
            checked_primary_relation_row(kc, fb, &[usize::MAX, 2, 165, 168], &target).is_err()
        );
        assert!(checked_primary_relation_row(kc, fb, &indices, &point_neg(&target)).is_err());
        let mut bad_labels = fb.clone();
        bad_labels.signed_orbit_of[2] = (1, 1, false);
        assert!(checked_primary_relation_row(kc, &bad_labels, &indices, &target).is_err());
        let mut rank = WideRankTracker::new(&kc.subgroup_order, fb.unknowns());
        for (column, indices) in [(0, [0; 4]), (1, [166; 4])] {
            let target = indices.iter().fold(BinaryPoint::Infinity, |sum, &index| {
                kc.add(&sum, &fb.points[index])
            });
            let row = checked_primary_relation_row(kc, fb, &indices, &target).unwrap();
            assert_eq!(rank.insert(row.clone()), column + 1);
            assert_eq!(rank.insert(row), column + 1);
        }
        let mut opts = KoblitzIcOptions::default();
        opts.max_trials = 0;
        opts.allow_direct_relation = false;
        let report = koblitz_index_calculus_dlp_with_factor_base(
            &base.curve,
            base.curve.generator(),
            &base.factor_base,
            &opts,
        )
        .unwrap();
        assert_eq!(report.factor_base_size, 332);
        assert_eq!(report.orbit_count, 2);
        assert_eq!(report.trials, 0);
        assert!(report.log.is_none());

        let mut entries: Vec<Value> = std::str::from_utf8(&bytes)
            .unwrap()
            .lines()
            .map(|line| serde_json::from_str(line).unwrap())
            .collect();
        let mut wrong_header = entries.clone();
        wrong_header[0]["seed"] = json!(18);
        let bad_header = wrong_header
            .iter()
            .map(Value::to_string)
            .collect::<Vec<_>>()
            .join("\n");
        assert!(from_record(&row, bad_header.as_bytes()).is_err());
        entries[1]["coefficient"] = json!("0");
        let bad = entries
            .iter()
            .map(Value::to_string)
            .collect::<Vec<_>>()
            .join("\n");
        assert!(from_record(&row, bad.as_bytes()).is_err());
    }

    #[test]
    fn wide_enumerator_matches_generic_witness_order_on_primary_curve() {
        let (bytes, row) = construct(
            0,
            "public_x_hash",
            2,
            17,
            Instant::now() + Duration::from_secs(60),
        )
        .unwrap();
        let base = from_record(&row, &bytes).unwrap();
        let kc = &base.curve;
        let fb = &base.factor_base;
        let index: HashMap<_, _> = fb
            .points
            .iter()
            .enumerate()
            .map(|(i, p)| (point_key(p), i))
            .collect();
        fn reference(
            kc: &KoblitzCurve,
            fb: &FrobeniusFactorBase,
            index: &HashMap<(BigUint, BigUint), usize>,
            target: &BinaryPoint,
            m: usize,
            start: usize,
        ) -> Option<Vec<usize>> {
            if m == 0 {
                return (target == &BinaryPoint::Infinity).then(Vec::new);
            }
            if m == 1 {
                let i = *index.get(&point_key(target))?;
                return (i >= start).then(|| vec![i]);
            }
            for i in start..fb.points.len() {
                let rest = kc.add(target, &point_neg(&fb.points[i]));
                if let Some(mut tail) = reference(kc, fb, index, &rest, m - 1, i) {
                    let mut result = vec![i];
                    result.append(&mut tail);
                    return Some(result);
                }
            }
            None
        }
        let p0 = &fb.points[0];
        let p2 = &fb.points[2];
        let two = kc.add(p0, p2);
        let three = kc.add(&kc.add(p0, p0), p2);
        let cases = [
            (BinaryPoint::Infinity, 0),
            (p0.clone(), 1),
            (BinaryPoint::Infinity, 2),
            (two, 2),
            (three, 3),
        ];
        for (target, m) in cases {
            assert_eq!(
                enumerate_decompose(kc, fb, &index, &target, m),
                reference(kc, fb, &index, &target, m, 0),
                "wide witness at m={m}"
            );
        }
    }

    #[test]
    fn public_target_and_post_solver_validation_bind_the_same_fixture() {
        let (bytes, row) = construct(
            0,
            "public_x_hash",
            2,
            17,
            Instant::now() + Duration::from_secs(60),
        )
        .unwrap();
        let base = from_record(&row, &bytes).unwrap();
        let suffix = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "n83-primary-fixture-{}-{suffix}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        let target = base.curve.generator().clone();
        let corpus = json!({
            "schema":"n83.public-probe-corpus/v1", "study":STUDY,
            "targets":[{"a":0,"fixture":0,"point":point_json(words(&target))}]
        });
        write_new_json(&root.join("probe-corpus.json"), &corpus).unwrap();
        let (selected, hash) = public_target(&root, &base.curve).unwrap();
        assert_eq!(selected, target);
        let validation = json!({
            "schema":"n83.public-probe-validation/v1",
            "corpus_canonical_json_blake3":hash,
            "oracle_receives_known_answer_scalars":false,
            "validation":[{"a":0,"fixture":0,"known_answer_scalar":"1"}]
        });
        write_new_json(&root.join("probe-validation.json"), &validation).unwrap();
        assert!(candidate_matches_validation(&root, &hash, &BigUint::one()).unwrap());
        assert!(!candidate_matches_validation(&root, &hash, &BigUint::from(2u8)).unwrap());
        let duplicate = json!({
            "schema":"n83.public-probe-corpus/v1", "study":STUDY,
            "targets":[
                {"a":0,"fixture":0,"point":point_json(words(&target))},
                {"a":0,"fixture":0,"point":point_json(words(&target))}
            ]
        });
        fs::write(
            root.join("probe-corpus.json"),
            serde_json::to_vec(&duplicate).unwrap(),
        )
        .unwrap();
        assert!(public_target(&root, &base.curve).is_err());
        let rejected_run = root.join("rejected-sat-m3");
        let error = run_cli(&root, 64, 3, "sat-m3", 1, 1, &rejected_run).unwrap_err();
        assert!(error.to_string().contains("gate is incomplete for F2^83"));
        assert!(!rejected_run.exists());
        fs::remove_dir_all(root).unwrap();
    }
}
