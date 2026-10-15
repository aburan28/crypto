//! Additional presentation fields, calculated only after full independent replay.
//! This never selects policies or recalculates any primary estimator/bootstrap.
use serde_json::{json, Value};
use std::collections::BTreeSet;

pub type Result<T> = std::result::Result<T, String>;
const IC_PHASES: [&str; 5] = [
    "target_query_ms",
    "target_pdp_ms",
    "target_relation_check_ms",
    "target_descent_ms",
    "target_recovery_check_ms",
];
const RHO_COUNTERS: [&str; 7] = [
    "walk_steps",
    "walks",
    "fruitless_cycles",
    "fruitless_escapes",
    "capped_walks",
    "wasted_merges",
    "table_entries",
];
fn rho_counters(raw: &Value) -> Result<Value> {
    RHO_COUNTERS
        .iter()
        .map(|key| {
            Ok((
                key.to_string(),
                json!(raw["outcome"][*key]
                    .as_u64()
                    .ok_or_else(|| format!("missing exact rho work counter {key}"))?),
            ))
        })
        .collect::<Result<serde_json::Map<_, _>>>()
        .map(Value::Object)
}
fn counter_totals(rows: &[Value]) -> Result<Value> {
    let mut totals = serde_json::Map::new();
    for key in RHO_COUNTERS {
        let mut sum = 0u64;
        for row in rows {
            sum = sum
                .checked_add(
                    row["work_counters"][key]
                        .as_u64()
                        .ok_or("missing exact exported counter")?,
                )
                .ok_or("rho counter sum overflow")?;
        }
        totals.insert(key.into(), json!(sum));
    }
    Ok(
        json!({"target_rows":rows.len(),"totals":totals,"scope":"Exact integer sums of exported work counters, not distinct simultaneous table memory. Internal capped walks/restarts are charged algorithmic work, separate from process/global timeouts."}),
    )
}
fn number(v: &Value, name: &str) -> Result<f64> {
    v.as_f64()
        .filter(|n| n.is_finite())
        .ok_or_else(|| format!("missing/nonfinite {name}"))
}
fn nonnegative(v: &Value, name: &str) -> Result<f64> {
    let n = number(v, name)?;
    if n < 0. {
        Err(format!("negative {name}"))
    } else {
        Ok(n)
    }
}
fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 16. * f64::EPSILON * a.abs().max(b.abs()).max(1.)
}

fn diagnostic_clock_key(key: &str) -> bool {
    matches!(
        key,
        "target_ms"
            | "target_query_ms"
            | "target_pdp_ms"
            | "target_relation_check_ms"
            | "target_descent_ms"
            | "target_recovery_check_ms"
            | "target_phase_sum_ms"
            | "target_generation_ms_excluded"
            | "online_total_ms"
            | "online_ms"
            | "phase_gap_ms"
            | "phase_sum_ms"
            | "online_wall_s"
            | "producer_preparation_ms"
            | "producer_shared_preparation_ms"
            | "producer_per_target_preparation_ms"
            | "elapsed_process_wall_s"
            | "cpu_user_s"
            | "cpu_system_s"
            | "cpu_total_s"
    )
}

/// JSON structure, integer types and all non-clock values remain exact.
/// At named diagnostic binary64 clock fields, tolerate one representable value
/// of parser/serializer round-trip drift, never a larger timing discrepancy.
pub fn diagnostic_source_equal(a: &Value, b: &Value) -> bool {
    fn one_ulp(a: &Value, b: &Value) -> bool {
        let (Some(an), Some(bn)) = (a.as_number(), b.as_number()) else {
            return false;
        };
        if !an.is_f64() || !bn.is_f64() {
            return false;
        }
        let (Some(a), Some(b)) = (an.as_f64(), bn.as_f64()) else {
            return false;
        };
        if !a.is_finite() || !b.is_finite() {
            return false;
        }
        let order = |n: f64| {
            let bits = n.to_bits();
            if bits >> 63 == 0 {
                bits | (1u64 << 63)
            } else {
                !bits
            }
        };
        order(a).abs_diff(order(b)) <= 1
    }
    fn walk(a: &Value, b: &Value, clock: bool, timing: bool) -> bool {
        if a == b {
            return true;
        }
        match (a, b) {
            (Value::Object(a), Value::Object(b)) => {
                a.len() == b.len()
                    && a.iter().all(|(k, v)| {
                        b.get(k).is_some_and(|w| {
                            walk(v, w, diagnostic_clock_key(k) || timing, k == "timing_ms")
                        })
                    })
            }
            (Value::Array(a), Value::Array(b)) => {
                a.len() == b.len() && a.iter().zip(b).all(|(v, w)| walk(v, w, false, false))
            }
            (Value::Number(_), Value::Number(_)) => clock && one_ulp(a, b),
            _ => false,
        }
    }
    walk(a, b, false, false)
}
pub fn descriptive(values: &[f64]) -> Value {
    let mut v = values.to_vec();
    v.sort_by(f64::total_cmp);
    let q = |p: f64| -> Option<f64> {
        if v.is_empty() {
            return None;
        }
        let i = p * (v.len() - 1) as f64;
        let lo = i.floor() as usize;
        let hi = i.ceil() as usize;
        Some(v[lo] + (v[hi] - v[lo]) * (i - lo as f64))
    };
    json!({"n":v.len(),"sum":if v.is_empty(){None}else{Some(v.iter().sum::<f64>())},"mean":if v.is_empty(){None}else{Some(v.iter().sum::<f64>()/v.len()as f64)},"median":q(0.5),"p90":q(0.9),"p95":q(0.95),"min":v.first(),"max":v.last(),"negative_rows":v.iter().filter(|x|**x<0.).count(),"quantile_rule":"linear interpolation at p*(n-1), midpoint median for even n","scope":"descriptive recorded rows only; no bootstrap or population-tail assertion"})
}

pub fn validate_gate(input: &Value, stats: &Value) -> Result<()> {
    let a = &input["independent_verification"];
    if input["schema"] != "n83_matched_statistics_input_v1"
        || input["metadata"]["independently_verified"] != true
        || stats["schema"] != "n83_matched_statistics_v1"
        || stats["all_required_rows_valid"] != true
        || a["status"] != "PASS"
        || a["checked_target_rows"] != 218
        || a["invalid_rows"] != 0
        || a["checked_aa_semantic_pairs"] != 5
        || a["aa_semantic_gate"] != true
        || a["campaign_manifest_sha256"] != input["metadata"]["campaign_manifest_sha256"]
        || a["protocol_sha256"] != input["metadata"]["protocol_sha256"]
        || !diagnostic_source_equal(&stats["independent_verification"], a)
        || stats["input_metadata"] != input["metadata"]
    {
        return Err("supplement requires matching full218-row independent PASS, five AA semantic pairs, and unchanged successful frozen analysis".into());
    }
    if !diagnostic_source_equal(&stats["process_rows"], &input["process_rows"])
        || !diagnostic_source_equal(&stats["target_rows"], &input["target_rows"])
        || !diagnostic_source_equal(
            &stats["auxiliary_target_rows"],
            &input["auxiliary_target_rows"],
        )
    {
        return Err("frozen analysis and normalized input rows differ".into());
    }
    Ok(())
}

pub fn ic_phase_row(row: &Value) -> Result<Value> {
    let raw = &row["rawtargetrecord"];
    let total = nonnegative(&raw["target_ms"], "IC target_ms")?;
    if total == 0. {
        return Err("zero IC target_ms".into());
    }
    let declared = nonnegative(&raw["target_phase_sum_ms"], "IC declared phase sum")?;
    let mut phases = serde_json::Map::new();
    let mut phase_sum = 0.;
    for key in IC_PHASES {
        let n = nonnegative(&raw[key], key)?;
        phase_sum += n;
        phases.insert(key.into(), json!(n));
    }
    if !close(phase_sum, declared) {
        return Err("IC five-phase sum is inconsistent with exported phase sum".into());
    }
    let residual = total - declared;
    Ok(
        json!({"run_id":row["run_id"],"kind":row["kind"],"target_id":row["target_id"],"arm":"ic","cohort_id":row["cohort_id"],"target_ms":total,"phases_ms":phases,"five_phase_sum_ms":declared,"recomputed_five_phase_sum_ms":phase_sum,"target_minus_phase_sum_ms":residual,"target_minus_phase_sum_fraction":residual/total,"group_verified":raw["group_verified"],"source":row["source"],"residual_scope":"Signed target_ms minus exported target_phase_sum_ms; retained without absolute value, clipping, reassignment or adding nested timings twice."}),
    )
}

pub fn rho_phase_row(row: &Value) -> Result<Value> {
    let raw = &row["rawtargetrecord"];
    let integer = |key: &str| {
        raw[key]
            .as_u64()
            .ok_or_else(|| format!("missing/noninteger rho {key}"))
    };
    let total = integer("online_interval_ns")?;
    let sum = integer("online_phase_sum_ns")?;
    let residual = integer("online_timer_boundary_residual_ns")?;
    if sum.checked_add(residual) != Some(total) {
        return Err("rho integer online/phase/residual identity failed".into());
    }
    let fields = [
        "target_query_and_validation",
        "walk_and_collision",
        "generic_scalar_replay",
    ];
    let mut phases = serde_json::Map::new();
    let mut recomputed = 0u64;
    for key in fields {
        let n = raw["online_phase_ns"][key]
            .as_u64()
            .ok_or_else(|| format!("missing rho phase {key}"))?;
        recomputed = recomputed.checked_add(n).ok_or("rho phase overflow")?;
        phases.insert(key.into(), json!(n as f64 / 1e6));
    }
    if recomputed != sum {
        return Err("rho integer phase components fail exact sum".into());
    }
    let counters = rho_counters(raw)?;
    Ok(
        json!({"run_id":row["run_id"],"kind":row["kind"],"target_id":row["target_id"],"arm":"rho","cohort_id":row["cohort_id"],"online_ms":total as f64/1e6,"phases_ms":phases,"phase_sum_ms":sum as f64/1e6,"boundary_residual_ms":residual as f64/1e6,"boundary_residual_ns":residual,"work_counters":counters,"source":row["source"],"scope":"Includes all114rho outputs, with ten A/A observations kept separately by kind; exact integer clock identity checked."}),
    )
}

pub fn process_row(row: &Value) -> Result<Value> {
    let wall = nonnegative(&row["elapsed_process_wall_s"], "process wall")?;
    let user = nonnegative(&row["cpu_user_s"], "user CPU")?;
    let system = nonnegative(&row["cpu_system_s"], "system CPU")?;
    let cpu = nonnegative(&row["cpu_total_s"], "total CPU")?;
    if wall == 0. || !close(user + system, cpu) {
        return Err("process CPU identity/positive wall failed".into());
    }
    let rss = row["peak_rss_bytes"]
        .as_u64()
        .ok_or("missing exact peak RSS bytes")?;
    let native = row["raw_process_receipt"]["peak_rss_native_kib"]
        .as_u64()
        .ok_or("missing native peak RSS KiB")?;
    if native.checked_mul(1024) != Some(rss) {
        return Err("process RSS KiB/bytes identity failed".into());
    }
    let prep = nonnegative(&row["producer_preparation_ms"], "preparation")? / 1000.;
    let timings = &row["summary"]["timing_ms"];
    let producer = nonnegative(&timings["process_total"], "producer total")? / 1000.;
    let targets = nonnegative(&timings["targets_total"], "producer targets total")? / 1000.;
    let start = row["raw_process_receipt"]["started_epoch_ms"].as_u64();
    let finish = row["raw_process_receipt"]["finished_epoch_ms"].as_u64();
    Ok(
        json!({"run_id":row["run_id"],"kind":row["kind"],"arm":row["arm"],"target_id":row["target_id"],"cohort_id":row["cohort_id"],"target_count":row["target_count"],"wall_s":wall,"cpu_user_s":user,"cpu_system_s":system,"cpu_total_s":cpu,"cpu_over_wall":cpu/wall,"peak_rss_bytes":rss,"peak_rss_mib":rss as f64/(1024.*1024.),"producer_preparation_s":prep,"producer_targets_total_s":targets,"producer_process_total_s":producer,"process_wall_minus_producer_total_s":wall-producer,"producer_total_minus_preparation_and_targets_s":producer-prep-targets,"started_epoch_ms":start,"finished_epoch_ms":finish,"producer_phase_ms":timings,"source":row["source"],"scope":"Actual whole-process wall/CPU/RSS. Signed residuals retained; setup_total already contains nested load/basis/index stages and is never added to them twice. Epoch millisecond stamps are scheduling records, not a replacement for monotonic wall."}),
    )
}

fn descriptions(rows: &[Value], field: &str) -> Value {
    descriptive(
        &rows
            .iter()
            .filter_map(|r| r[field].as_f64())
            .collect::<Vec<_>>(),
    )
}
fn phase_summaries(rows: &[Value], fields: &[&str]) -> Value {
    fields
        .iter()
        .map(|key| {
            (
                key.to_string(),
                descriptive(
                    &rows
                        .iter()
                        .filter_map(|r| r["phases_ms"][*key].as_f64())
                        .collect::<Vec<_>>(),
                ),
            )
        })
        .collect::<serde_json::Map<_, _>>()
        .into()
}

pub fn supplemental(input: &Value, stats: &Value) -> Result<Value> {
    validate_gate(input, stats)?;
    let mut ic = Vec::new();
    let mut rho = Vec::new();
    let mut seen = BTreeSet::new();
    let targets = input["target_rows"]
        .as_array()
        .ok_or("missing cohort target rows")?;
    let auxiliary = input["auxiliary_target_rows"]
        .as_array()
        .ok_or("missing auxiliary target rows")?;
    for row in targets.iter().chain(auxiliary) {
        let key = (
            row["run_id"].as_str().ok_or("missing run id")?,
            row["target_id"].as_str().ok_or("missing target id")?,
        );
        if !seen.insert(key) {
            return Err("duplicate supplemental target identity".into());
        }
        if row["verified"] != true || row["independently_verified"] != true {
            return Err("supplement requires each target independently verified".into());
        }
        match row["arm"].as_str() {
            Some("ic") => ic.push(ic_phase_row(row)?),
            Some("rho") => rho.push(rho_phase_row(row)?),
            _ => return Err("unknown target arm".into()),
        }
    }
    if ic.len() != 104 || rho.len() != 114 {
        return Err(format!(
            "expected104IC/114rho targets, got{}/{}",
            ic.len(),
            rho.len()
        ));
    }
    let mut processes = Vec::new();
    for row in input["process_rows"]
        .as_array()
        .ok_or("missing process rows")?
    {
        processes.push(process_row(row)?);
    }
    if processes.len() != 32 {
        return Err("expected32process rows".into());
    }
    let mut groups = Vec::new();
    for arm in ["ic", "rho"] {
        for kind in ["cohort", "cold", "aa"] {
            let ps = processes
                .iter()
                .filter(|r| r["arm"] == arm && r["kind"] == kind)
                .cloned()
                .collect::<Vec<_>>();
            let ts = if arm == "ic" { &ic } else { &rho }
                .iter()
                .filter(|r| r["kind"] == kind)
                .cloned()
                .collect::<Vec<_>>();
            let work = if arm == "rho" {
                counter_totals(&ts)?
            } else {
                Value::Null
            };
            groups.push(json!({"arm":arm,"kind":kind,"processes":ps.len(),"targets":ts.len(),"wall_s":descriptions(&ps,"wall_s"),"cpu_total_s":descriptions(&ps,"cpu_total_s"),"cpu_user_s":descriptions(&ps,"cpu_user_s"),"cpu_system_s":descriptions(&ps,"cpu_system_s"),"cpu_over_wall":descriptions(&ps,"cpu_over_wall"),"peak_rss_bytes":descriptions(&ps,"peak_rss_bytes"),"peak_rss_mib":descriptions(&ps,"peak_rss_mib"),"process_wall_minus_producer_total_s":descriptions(&ps,"process_wall_minus_producer_total_s"),"producer_total_minus_preparation_and_targets_s":descriptions(&ps,"producer_total_minus_preparation_and_targets_s"),"online_ms":descriptions(&ts,if arm=="ic"{"target_ms"}else{"online_ms"}),"phase_ms":phase_summaries(&ts,if arm=="ic"{&IC_PHASES}else{&["target_query_and_validation","walk_and_collision","generic_scalar_replay"]}),"online_boundary_residual_ms":descriptions(&ts,if arm=="ic"{"target_minus_phase_sum_ms"}else{"boundary_residual_ms"}),"rho_work_counter_totals":work}));
        }
    }
    let main_rho = rho
        .iter()
        .filter(|r| r["kind"] != "aa")
        .cloned()
        .collect::<Vec<_>>();
    let mut rho_by_process = Vec::new();
    for process in processes.iter().filter(|p| p["arm"] == "rho") {
        let rows = rho
            .iter()
            .filter(|r| r["run_id"] == process["run_id"])
            .cloned()
            .collect::<Vec<_>>();
        rho_by_process.push(json!({"run_id":process["run_id"],"kind":process["kind"],"cohort_id":process["cohort_id"],"work":counter_totals(&rows)?}));
    }
    Ok(
        json!({"schema":"n83_matched_supplement_v1","purpose":"Native presentation-only phase, residual, CPU and memory fields omitted from frozen main analysis; primary statistics/bootstrap/configuration unchanged.","all_ic_targets":ic,"all_rho_targets":rho,"all_processes":processes,"group_summaries":groups,"ic104_phase_ms":phase_summaries(&ic,&IC_PHASES),"ic104_target_minus_phase_sum_ms":descriptions(&ic,"target_minus_phase_sum_ms"),"ic104_target_minus_phase_sum_fraction":descriptions(&ic,"target_minus_phase_sum_fraction"),"all114_rho_work_counter_totals":counter_totals(&rho)?,"main104_rho_work_counter_totals":counter_totals(&main_rho)?,"rho_work_counters_by_process":rho_by_process,"counts":{"ic_targets":104,"rho_targets":114,"all_targets":218,"processes":32},"independent_verification":input["independent_verification"],"limitations":["Supplemental descriptive fields only; no new primary estimator, bootstrap, optimization or ledger admission.","Phase boundaries differ across IC and rho. IC phase residual cannot be reassigned to a particular algorithm stage.","Process peak RSS is not per-target peak memory; it is not summed as simultaneously required memory.","All signed timing residuals are retained. Producer process wall and supervisor whole-process wall have different boundaries.","Zero process timeouts does not imply zero internal capped walks/restarts. Ten A/A rows remain separate from104 main rho answers."]}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ic() -> Value {
        json!({"run_id":"mock","target_id":"mock","rawtargetrecord":{"target_ms":11.,"target_query_ms":1.,"target_pdp_ms":2.,"target_relation_check_ms":3.,"target_descent_ms":1.,"target_recovery_check_ms":3.,"target_phase_sum_ms":10.,"group_verified":true}})
    }
    #[test]
    fn ic_residual_is_signed_and_not_a_sixth_phase() {
        let r = ic_phase_row(&ic()).unwrap();
        assert_eq!(r["target_minus_phase_sum_ms"], 1.);
        let mut row = ic();
        row["rawtargetrecord"]["target_ms"] = json!(9.);
        assert_eq!(
            ic_phase_row(&row).unwrap()["target_minus_phase_sum_ms"],
            -1.
        );
    }
    #[test]
    fn missing_and_inconsistent_phase_never_become_zero() {
        let mut r = ic();
        r["rawtargetrecord"]
            .as_object_mut()
            .unwrap()
            .remove("target_pdp_ms");
        assert!(ic_phase_row(&r).is_err());
        let mut r = ic();
        r["rawtargetrecord"]["target_phase_sum_ms"] = json!(11.);
        assert!(ic_phase_row(&r).is_err());
    }
    #[test]
    fn even_median_and_negative_counts() {
        let d = descriptive(&[-2., 1., 3., 8.]);
        assert_eq!(d["median"], 2.);
        assert_eq!(d["negative_rows"], 1);
        assert!(descriptive(&[])["mean"].is_null());
    }
    #[test]
    fn rho_integer_identity_is_exact() {
        let mut r = json!({"rawtargetrecord":{"online_interval_ns":20,"online_phase_sum_ns":19,"online_timer_boundary_residual_ns":1,"online_phase_ns":{"target_query_and_validation":2,"walk_and_collision":15,"generic_scalar_replay":2}}});
        r["rawtargetrecord"]["outcome"] = RHO_COUNTERS
            .iter()
            .map(|k| (k.to_string(), json!(0)))
            .collect::<serde_json::Map<_, _>>()
            .into();
        assert!(rho_phase_row(&r).is_ok());
        r["rawtargetrecord"]["online_timer_boundary_residual_ns"] = json!(2);
        assert!(rho_phase_row(&r).is_err());
    }
    #[test]
    fn cpu_and_rss_identity_not_stage_double_counting() {
        let mut r = json!({"elapsed_process_wall_s":4.,"cpu_user_s":2.,"cpu_system_s":1.,"cpu_total_s":3.,"peak_rss_bytes":2048,"raw_process_receipt":{"peak_rss_native_kib":2},"producer_preparation_ms":1000.,"summary":{"timing_ms":{"process_total":3000.,"targets_total":1000.,"setup_total":1000.,"base_load":400.,"normal_basis_and_selftest":100.,"index_build":500.}}});
        let s = process_row(&r).unwrap();
        assert_eq!(s["producer_total_minus_preparation_and_targets_s"], 1.);
        assert_eq!(s["cpu_over_wall"], 0.75);
        r["peak_rss_bytes"] = json!(2049);
        assert!(process_row(&r).is_err());
    }
    #[test]
    fn failed_or_missing_independent_gate_blocks_execution() {
        assert!(validate_gate(&json!({}), &json!({})).is_err());
    }

    #[test]
    fn counter_sums_remain_exact_and_overflow_fails() {
        let counters = |n: u64| {
            RHO_COUNTERS
                .iter()
                .map(|k| (k.to_string(), json!(n)))
                .collect::<serde_json::Map<_, _>>()
        };
        let rows = vec![
            json!({"work_counters":counters(7)}),
            json!({"work_counters":counters(11)}),
        ];
        assert_eq!(counter_totals(&rows).unwrap()["totals"]["capped_walks"], 18);
        let rows = vec![
            json!({"work_counters":counters(u64::MAX)}),
            json!({"work_counters":counters(1)}),
        ];
        assert!(counter_totals(&rows).is_err());
    }

    #[test]
    fn only_named_float_clock_fields_allow_one_ulp() {
        let a = json!({"target_ms":1.0,"counter":1,"hash":"original","online_interval_ns":1});
        let mut b = a.clone();
        b["target_ms"] = json!(f64::from_bits(1.0f64.to_bits() + 1));
        assert!(diagnostic_source_equal(&a, &b));
        b["target_ms"] = json!(f64::from_bits(1.0f64.to_bits() + 2));
        assert!(!diagnostic_source_equal(&a, &b));
        b = a.clone();
        b["counter"] = json!(1.0);
        assert!(!diagnostic_source_equal(&a, &b));
        b = a.clone();
        b["hash"] = json!("changed");
        assert!(!diagnostic_source_equal(&a, &b));
        b = a.clone();
        b["online_interval_ns"] = json!(2);
        assert!(!diagnostic_source_equal(&a, &b));
        let x = json!({"timing_ms":{"rank_stage":1.0}});
        let y = json!({"timing_ms":{"rank_stage":f64::from_bits(1.0f64.to_bits()+1)}});
        assert!(diagnostic_source_equal(&x, &y));
    }
}
