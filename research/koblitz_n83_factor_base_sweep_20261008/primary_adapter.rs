//! Replay-bound bridge from the stored primary N83 points to the generic IC base.
//!
//! The library's explicit-orbit constructor currently requires a u64 field.
//! This bridge uses the verified JSONL orbit order directly and checks every
//! point and signed-Frobenius label before building the public structure.

use super::{
    curve, general, load_object, parse_point, point_json, point_set_hash, words, write_new_json,
    Result, SEEDS, STUDY,
};
use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    koblitz_index_calculus_dlp_with_factor_base_and_progress, koblitz_point_count, point_key,
    DecompositionStrategy, FactorBaseDomain, FrobeniusFactorBase, KoblitzCurve, KoblitzIcOptions,
};
use num_bigint::BigUint;
use num_traits::One;
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::time::{Duration, Instant};

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

fn load_panel(root: &Path, columns: usize) -> Result<(PrimaryBase, Value, String)> {
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
            && row["policy"] == "public_x_hash"
            && row["seed"].as_u64() == Some(SEEDS[0])
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
        koblitz_index_calculus_dlp_with_factor_base, KoblitzIcOptions,
    };
    use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

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
