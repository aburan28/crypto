//! Analyze the sealed Round 314 exact-log-transport ecbench sessions.
//!
//! This is deliberately an operation-count analysis.  Every session is L0,
//! so wall time is retained in the records but is never promoted to a speed
//! comparison.  The P-256 transfer decision is read from the existing exact
//! FB1 transport census rather than inferred from the Koblitz experiments.

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ecbench::record::Record;
use crypto_lib::cryptanalysis::ecbench::runner::{read_records, read_session};
use crypto_lib::cryptanalysis::ecbench::spec::Level;
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

const ROOT: &str = "research/ecbench_known_log_transport_round314_20261014";
const P256_TRANSPORT: &str = "research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_fb_log_transport_round21_20261006/transport-result.json";
const CURVE: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
const FIXED_FB: &str = "FB1h2f8621cda105";
const ARMS: [&str; 5] = [
    "rho-strong",
    "transport-hash",
    "transport-hash-aa",
    "transport-sequential",
    "claw-c1",
];

struct Input {
    tag: &'static str,
    session: &'static str,
    audit: &'static str,
    expected_session: &'static str,
}

const INPUTS: [Input; 4] = [
    Input {
        tag: "n37",
        session: "sessions/iteration-n37",
        audit: "ITERATION_N37_AUDIT.json",
        expected_session: "ECBS1h9245b38c90f5",
    },
    Input {
        tag: "n43",
        session: "sessions/iteration-n43",
        audit: "ITERATION_N43_AUDIT.json",
        expected_session: "ECBS1hadc6338f9904",
    },
    Input {
        tag: "n47",
        session: "sessions/iteration-n47",
        audit: "ITERATION_N47_AUDIT.json",
        expected_session: "ECBS1h02a93eb07dc0",
    },
    // This curve has the largest subgroup even though its field degree is 41.
    Input {
        tag: "n41",
        session: "sessions/iteration-n41",
        audit: "ITERATION_N41_AUDIT.json",
        expected_session: "ECBS1ha4734a5bf3e0",
    },
];

fn bytes(path: &Path) -> Result<Vec<u8>, String> {
    fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}

fn value(path: &Path) -> Result<Value, String> {
    serde_json::from_slice(&bytes(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

fn mean(values: impl Iterator<Item = f64>) -> f64 {
    let values: Vec<_> = values.collect();
    values.iter().sum::<f64>() / values.len() as f64
}

fn phase_signature(record: &Record) -> Vec<Value> {
    record
        .phases
        .iter()
        .map(|p| {
            json!({
                "name": p.name,
                "adds": p.adds,
                "doubles": p.doubles,
                "scalar_mults": p.scalar_mults,
                "gae": p.gae,
                "native": p.native,
            })
        })
        .collect()
}

fn checked_audit(path: &Path, dir: &Path, session_id: &str) -> Result<Value, String> {
    let audit = value(path)?;
    if audit["ok"] != true
        || audit["session_status"] != "complete"
        || audit["session_id"] != session_id
        || audit["records"] != 240
        || audit["verified_records"] != 240
        || audit["replays"].as_array().map(Vec::len) != Some(200)
        || audit["replays"]
            .as_array()
            .is_none_or(|r| r.iter().any(|x| x["reproduced"] != true))
    {
        return Err(format!(
            "{} is not a complete 200-replay audit",
            path.display()
        ));
    }
    for name in [
        "host.json",
        "plan.json",
        "records.jsonl",
        "session.json",
        "spec.json",
    ] {
        let actual = sha256_hex(&bytes(&dir.join(name))?);
        if audit["files"][name] != actual {
            return Err(format!("{} has a stale {name} hash", path.display()));
        }
    }
    Ok(audit)
}

fn checked_comparison(dir: &Path, a: &str, b: &str) -> Result<Value, String> {
    let path = dir.join(format!("comparisons/{a}__{b}.json"));
    let comparison = value(&path)?;
    if comparison["ops"]["status"] != "ok"
        || comparison["ops"]["pairs"] != 40
        || comparison["ops"]["same_seeds"] != true
        || comparison["ops"]["same_session"] != true
        || comparison["ops"]["bounded"] != true
        || comparison["wall"]["status"] != "descriptive"
    {
        return Err(format!(
            "{} is not the expected paired L0 comparison",
            path.display()
        ));
    }
    Ok(comparison)
}

fn checked_record(record: &Record, session_id: &str) -> Result<(), String> {
    if record.session_id != session_id
        || record.warmup != (record.round == 0)
        || record.outcome.status != "verified"
        || record.outcome.matches_target != Some(true)
        || !record.cost.lower_bound
        || !record.cost.deterministic
        || record.cost.unpriced.is_empty()
        || record.isolation.level != Level::L0
    {
        return Err(format!(
            "{} fails the status/accounting gate",
            record.record_id
        ));
    }
    let charged = record.cost.total_gae.ok_or("record lacks total GAE")?;
    let phase_sum: f64 = record.phases.iter().map(|p| p.gae).sum();
    if (charged - phase_sum).abs() > 1e-7_f64.max(charged * 1e-11) {
        return Err(format!(
            "{} has nonexclusive charged phases",
            record.record_id
        ));
    }
    Ok(())
}

fn arm_summary(records: &[&Record]) -> Result<Value, String> {
    if records.len() != 40 {
        return Err(format!(
            "arm has {} measured records, expected 40",
            records.len()
        ));
    }
    let mut phase_totals = BTreeMap::<String, f64>::new();
    let mut native_totals = BTreeMap::<String, u64>::new();
    let mut adds = 0u64;
    let mut doubles = 0u64;
    let mut scalar_mults = 0u64;
    let mut max_rss_kib = 0u64;
    for record in records {
        max_rss_kib = max_rss_kib.max(record.time.max_rss_kib);
        for phase in &record.phases {
            *phase_totals.entry(phase.name.clone()).or_default() += phase.gae;
            adds += phase.adds;
            doubles += phase.doubles;
            scalar_mults += phase.scalar_mults;
            for (name, count) in &phase.native {
                *native_totals.entry(name.clone()).or_default() += count;
            }
        }
    }
    let divisor = records.len() as f64;
    let phase_means: BTreeMap<_, _> = phase_totals
        .into_iter()
        .map(|(name, total)| (name, total / divisor))
        .collect();
    let mean_counter = |name: &str| {
        mean(
            records
                .iter()
                .map(|r| r.counters.get(name).copied().unwrap_or(0) as f64),
        )
    };
    Ok(json!({
        "method": records[0].method.id,
        "measured_runs": records.len(),
        "verified_runs": records.len(),
        "mean_cold_gae_lower_bound": mean(records.iter().map(|r| r.cost.total_gae.unwrap())),
        "mean_cold_s_lower_bound": mean(records.iter().map(|r| r.cost.s.unwrap())),
        "mean_ratio_to_generic_floor_lower_bound": mean(records.iter().map(|r| r.boundaries.ratio_to_floor.expect("checked ecbench boundary"))),
        "mean_group_ops": {
            "adds": adds as f64 / divisor,
            "doubles": doubles as f64 / divisor,
            "scalar_mults": scalar_mults as f64 / divisor,
        },
        "mean_target_attempts": mean_counter("target_attempts"),
        "mean_pair_table_entries": mean_counter("pair_table_entries"),
        "mean_factor_base_bytes_floor": mean_counter("factor_base_vec_bytes_floor_uncharged"),
        "peak_process_rss_kib": max_rss_kib,
        "phase_mean_gae_lower_bound": phase_means,
        "phase_native_totals": native_totals,
    }))
}

fn fit(points: &[(String, f64, f64)]) -> Value {
    let n = points.len() as f64;
    let mean_x = points.iter().map(|p| p.1).sum::<f64>() / n;
    let mean_y = points.iter().map(|p| p.2).sum::<f64>() / n;
    let sxx = points.iter().map(|p| (p.1 - mean_x).powi(2)).sum::<f64>();
    let slope = points
        .iter()
        .map(|p| (p.1 - mean_x) * (p.2 - mean_y))
        .sum::<f64>()
        / sxx;
    let intercept = mean_y - slope * mean_x;
    let residuals: Vec<Value> = points
        .iter()
        .map(|(tag, x, y)| {
            let predicted = intercept + slope * x;
            json!({
                "panel": tag,
                "log2_r": x,
                "observed_log2": y,
                "predicted_log2": predicted,
                "residual_log2": y - predicted,
            })
        })
        .collect();
    let sse = points
        .iter()
        .map(|(_, x, y)| (y - intercept - slope * x).powi(2))
        .sum::<f64>();
    let syy = points.iter().map(|p| (p.2 - mean_y).powi(2)).sum::<f64>();
    let standard_error = (sse / (n - 2.0) / sxx).sqrt();
    // Student-t(0.975, 2), because this frozen scaling panel has four points.
    let t95 = 4.302_652_729_696_142;
    json!({
        "model": "log2(metric) = intercept + exponent*log2(r)",
        "points": points.len(),
        "exponent": slope,
        "intercept": intercept,
        "r_squared": if syy == 0.0 { 1.0 } else { 1.0 - sse / syy },
        "slope_standard_error": standard_error,
        "slope_ci95_t_df2": [slope - t95 * standard_error, slope + t95 * standard_error],
        "residuals": residuals,
    })
}

fn extract(row: &Value, pointer: &str) -> Result<f64, String> {
    row.pointer(pointer)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("missing numeric result field {pointer}"))
}

fn run(root: &Path, p256_path: &Path, out_path: &Path) -> Result<(), String> {
    let mut rows = Vec::new();
    let mut input_hashes = Map::new();
    let expected_arms: BTreeSet<_> = ARMS.into_iter().collect();
    let mut total_records = 0usize;
    let mut total_measured = 0usize;
    let mut total_replays = 0usize;

    for input in &INPUTS {
        let dir = root.join(input.session);
        let audit_path = root.join(input.audit);
        let session = read_session(&dir)?;
        if session.session_id != input.expected_session || session.status != "complete" {
            return Err(format!(
                "{} has wrong session identity/status",
                dir.display()
            ));
        }
        let audit = checked_audit(&audit_path, &dir, input.expected_session)?;
        let records = read_records(&dir)?;
        if records.len() != 240 || records.iter().filter(|r| r.warmup).count() != 40 {
            return Err(format!(
                "{} lacks 40 warm-ups plus 200 measured runs",
                dir.display()
            ));
        }
        total_records += records.len();
        total_measured += records.iter().filter(|r| !r.warmup).count();
        total_replays += audit["replays"].as_array().unwrap().len();
        let mut by_arm = BTreeMap::<String, Vec<&Record>>::new();
        let mut by_pair = BTreeMap::<(String, u32), BTreeMap<String, &Record>>::new();
        for record in &records {
            checked_record(record, input.expected_session)?;
            if !expected_arms.contains(record.arm.as_str()) {
                return Err(format!(
                    "{} has unexpected arm {}",
                    dir.display(),
                    record.arm
                ));
            }
            if !record.warmup {
                by_arm.entry(record.arm.clone()).or_default().push(record);
                by_pair
                    .entry((record.workload.workload_id.clone(), record.round))
                    .or_default()
                    .insert(record.arm.clone(), record);
            }
        }
        if by_pair.len() != 40 || by_arm.len() != ARMS.len() {
            return Err(format!(
                "{} does not contain the complete paired panel",
                dir.display()
            ));
        }
        for ((workload, round), arms) in &by_pair {
            if arms.keys().map(String::as_str).collect::<BTreeSet<_>>() != expected_arms {
                return Err(format!("{workload} round {round} lacks an arm"));
            }
            let reference = arms["rho-strong"];
            for record in arms.values() {
                if record.workload.target != reference.workload.target
                    || record.algorithm_seed != reference.algorithm_seed
                    || record.outcome.recovered != reference.outcome.recovered
                {
                    return Err(format!(
                        "{workload} round {round} is not a same-point solve"
                    ));
                }
            }
            let hash = arms["transport-hash"];
            let aa = arms["transport-hash-aa"];
            if hash.cost.total_gae != aa.cost.total_gae
                || hash.cost.s != aa.cost.s
                || hash.counters != aa.counters
                || hash.factor_base != aa.factor_base
                || phase_signature(hash) != phase_signature(aa)
            {
                return Err(format!("{workload} round {round} failed exact counted A/A"));
            }
        }
        let rho_hash = checked_comparison(&dir, "rho-strong", "transport-hash")?;
        let seq_hash = checked_comparison(&dir, "transport-sequential", "transport-hash")?;
        let aa = checked_comparison(&dir, "transport-hash", "transport-hash-aa")?;
        let rho_seq = checked_comparison(&dir, "rho-strong", "transport-sequential")?;
        let rho_claw = checked_comparison(&dir, "rho-strong", "claw-c1")?;
        if aa["ops"]["ratio_b_over_a"] != 1.0 || aa["ops"]["ci95"] != json!([1.0, 1.0]) {
            return Err(format!(
                "{} has nonidentity A/A counted work",
                dir.display()
            ));
        }
        let first = records.first().ok_or("empty records")?;
        let r = first.workload.curve.r as f64;
        let mut arms = Map::new();
        for arm in ARMS {
            arms.insert(arm.into(), arm_summary(&by_arm[arm])?);
        }
        rows.push(json!({
            "panel": input.tag,
            "session_id": session.session_id,
            "spec_id": session.spec_id,
            "binary_sha256": session.binary_sha256,
            "curve": first.workload.curve.slug,
            "field_degree": first.workload.curve.field_degree,
            "subgroup_order": first.workload.curve.r,
            "log2_subgroup_order": r.log2(),
            "target_count": 8,
            "measured_rounds_per_target": 5,
            "all_records_verified": true,
            "all_measured_replayed_identically": true,
            "isolation_level": "L0",
            "arms": arms,
            "paired_count_ratios": {
                "hash_over_rho": rho_hash["ops"]["ratio_b_over_a"],
                "hash_over_rho_ci95": rho_hash["ops"]["ci95"],
                "hash_over_sequential": seq_hash["ops"]["ratio_b_over_a"],
                "hash_over_sequential_ci95": seq_hash["ops"]["ci95"],
                "sequential_over_rho": rho_seq["ops"]["ratio_b_over_a"],
                "claw_over_rho": rho_claw["ops"]["ratio_b_over_a"],
                "hash_aa": aa["ops"]["ratio_b_over_a"],
            },
            "wall_time_admissible": false,
            "audit_sha256": sha256_hex(&bytes(&audit_path)?),
            "records_sha256": sha256_hex(&bytes(&dir.join("records.jsonl"))?),
        }));
        input_hashes.insert(
            input.audit.into(),
            Value::String(sha256_hex(&bytes(&audit_path)?)),
        );
        input_hashes.insert(
            format!("{}/records.jsonl", input.session),
            Value::String(sha256_hex(&bytes(&dir.join("records.jsonl"))?)),
        );
    }

    rows.sort_by(|a, b| {
        a["subgroup_order"]
            .as_u64()
            .cmp(&b["subgroup_order"].as_u64())
    });
    let points = |pointer: &str| -> Result<Vec<(String, f64, f64)>, String> {
        rows.iter()
            .map(|row| {
                Ok((
                    row["panel"].as_str().unwrap().to_string(),
                    extract(row, "/log2_subgroup_order")?,
                    extract(row, pointer)?.log2(),
                ))
            })
            .collect()
    };
    let work_fit = fit(&points("/arms/transport-hash/mean_cold_gae_lower_bound")?);
    let rho_fit = fit(&points("/arms/rho-strong/mean_cold_gae_lower_bound")?);
    let table_fit = fit(&points("/arms/transport-hash/mean_pair_table_entries")?);
    let rss_fit = fit(&points("/arms/transport-hash/peak_process_rss_kib")?);
    let attempts_fit = fit(&points("/arms/transport-hash/mean_target_attempts")?);

    let p256_bytes = bytes(p256_path)?;
    let p256: Value = serde_json::from_slice(&p256_bytes).map_err(|e| e.to_string())?;
    if p256["schema"] != "p256.fb_log_transport/v1"
        || p256["curve"] != CURVE
        || p256["factor_base"]["fb_id"] != FIXED_FB
        || p256["quotient"]["total_independent_relation_rank"] != 0
        || p256["quotient"]["quotient_dimension"] != 131_458
        || p256["false_positives"] != 0
        || p256["false_negatives"] != 0
        || p256["exact_group_replay"] != true
    {
        return Err("P-256 transfer certificate does not match the frozen negative result".into());
    }
    let p256_block_images: u64 = p256["block_censuses"]
        .as_array()
        .ok_or("P-256 block censuses missing")?
        .iter()
        .map(|row| {
            row["registered_images"]
                .as_u64()
                .ok_or("P-256 block image count missing")
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .sum();
    input_hashes.insert(
        P256_TRANSPORT.into(),
        Value::String(sha256_hex(&p256_bytes)),
    );

    let degree37 = rows
        .iter()
        .find(|r| r["field_degree"] == 37)
        .ok_or("degree-37 row missing")?;
    let degree37_ratio = extract(degree37, "/paired_count_ratios/hash_over_rho")?;
    let exponent = work_fit["exponent"].as_f64().unwrap();
    let correctness_gate = total_records == 960 && total_measured == 800 && total_replays == 800;
    let aa_gate = rows
        .iter()
        .all(|r| r["paired_count_ratios"]["hash_aa"] == 1.0);
    let output = json!({
        "schema": "ecbench.known_log_transport.round314.analysis/v1",
        "date": "2026-10-14",
        "status": "verified_l0_operation_count_negative_result",
        "question": "Can a fully charged exact-log signed-Frobenius base reach cold one-target parity with strong Pollard rho, and can that mechanism transfer to the fixed P-256 Dickson base?",
        "comparison_unit": "ecbench.gae; S = cold total GAE / sqrt(subgroup order)",
        "rows": rows,
        "scaling_fits": {
            "hash_candidate_cold_work": work_fit,
            "strong_rho_cold_work": rho_fit,
            "hash_candidate_pair_entries": table_fit,
            "hash_candidate_peak_rss_kib": rss_fit,
            "hash_candidate_target_attempts": attempts_fit,
        },
        "correctness": {
            "records": total_records,
            "measured_records": total_measured,
            "measured_records_replayed_identically": total_replays,
            "wrong_answers": 0,
            "exhausted_runs": 0,
            "crashes": 0,
            "timeouts": 0,
            "observed_false_positives": 0,
            "observed_false_negatives_on_checked_instances": 0,
            "exact_counted_aa_pairs": 160,
        },
        "p256_transfer": {
            "curve": CURVE,
            "factor_base": FIXED_FB,
            "certificate_schema": p256["schema"],
            "certificate_sha256": sha256_hex(&p256_bytes),
            "small_multiple_images": p256["multiplier_census"]["images"],
            "dickson_block_images": p256_block_images,
            "transport_relations": p256["quotient"]["total_independent_relation_rank"],
            "remaining_log_dimensions": p256["quotient"]["quotient_dimension"],
            "structured_residual_degree": p256["residual_degree_evidence"]["structured_residual_maximum"],
            "unsplit_s17_degree_of_regularity": p256["residual_degree_evidence"]["unsplit_s17_degree_of_regularity"],
            "transfer_passes": false,
            "reason": "Known multiples have exact logs only by changing the point set. The complete registered census finds no admitted log transport inside FB1h2f8621cda105; obtaining a fixed FB1 point's scalar is the ECDLP being solved.",
        },
        "promotion_gates": {
            "all_records_verified_and_replayed": correctness_gate,
            "exact_counted_aa": aa_gate,
            "degree37_hash_over_rho_at_most_1_30": degree37_ratio <= 1.30,
            "candidate_work_exponent_at_most_0_55": exponent <= 0.55,
            "cold_setup_included": true,
            "p256_exact_membership_and_log_transport": false,
            "all_conjunctive_gates_pass": false,
        },
        "decision": {
            "classification": "engineering improvement; negative parity and transfer result",
            "parity_achieved": false,
            "full_p256_attempt_authorized": false,
            "dominant_obstruction": "The exact-log construction is a generic scalar-labelled base, not the fixed Dickson base. Even on applicable Koblitz analogues its fully charged cold count remains above strong rho, while the exact P-256 census leaves all 131,458 FB1 log dimensions independent.",
            "round_stop": "The preregistered single post-holdout iteration has been consumed; another algorithmic change requires a new preregistration and fresh hidden targets.",
        },
        "limits": [
            "Every session earned L0, so wall-time ratios are descriptive and excluded.",
            "Field arithmetic, hashing, allocation, and modular combination remain explicitly unpriced; all reported operation ratios are lower-bound diagnostics.",
            "Four finite scaling points do not prove an asymptotic exponent.",
            "The transported-log method bypasses a polynomial system; it does not establish a new Groebner degree. The fixed P-256 base retains its independently certified structured residual degree 4 and unknown unsplit S17 degree.",
            "No Koblitz measurement is extrapolated into a P-256 attack cost after the exact transfer gate fails.",
        ],
        "input_sha256": input_hashes,
    });
    let mut out = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(out_path)
        .map_err(|e| format!("{}: {e}", out_path.display()))?;
    serde_json::to_writer_pretty(&mut out, &output).map_err(|e| e.to_string())?;
    out.write_all(b"\n").map_err(|e| e.to_string())?;
    Ok(())
}

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let out = args.next().ok_or(
        "usage: known_log_transport_round314_analyze OUT [RESEARCH_ROOT [P256_CERTIFICATE]]",
    )?;
    let root = args.next().unwrap_or_else(|| ROOT.into());
    let p256 = args.next().unwrap_or_else(|| P256_TRANSPORT.into());
    if args.next().is_some() {
        return Err("too many arguments".into());
    }
    run(Path::new(&root), Path::new(&p256), &PathBuf::from(out))
}
