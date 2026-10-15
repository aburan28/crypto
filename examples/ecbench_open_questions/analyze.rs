//! Analyze frozen records only; never recompute or alter a session.
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
};

fn read(path: &Path) -> Result<Value, String> {
    serde_json::from_slice(&fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?)
        .map_err(|e| e.to_string())
}

pub fn run(root: &Path) -> Result<(), String> {
    let mut dirs: Vec<_> = fs::read_dir(root.join("sessions"))
        .map_err(|e| e.to_string())?
        .map(|e| e.map(|e| e.path()))
        .collect::<Result<_, _>>()
        .map_err(|e| e.to_string())?;
    dirs.sort();
    let mut all = Vec::new();
    let mut rows = Vec::new();
    for dir in dirs {
        if !dir.is_dir() {
            continue;
        }
        let name = dir.file_name().unwrap().to_str().unwrap();
        let audit = read(&root.join(format!("{name}.audit.json")))?;
        let session = read(&dir.join("session.json"))?;
        let text = fs::read_to_string(dir.join("records.jsonl")).map_err(|e| e.to_string())?;
        let records: Vec<Value> = text
            .lines()
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()
            .map_err(|e| e.to_string())?;
        let measured: Vec<_> = records.iter().filter(|r| r["warmup"] == false).collect();
        let mut groups: BTreeMap<(String, String), Vec<&Value>> = BTreeMap::new();
        for r in &measured {
            let slug = r["workload"]["curve"]["slug"]
                .as_str()
                .ok_or("missing curve slug")?;
            let arm = r["arm"].as_str().ok_or("missing arm")?;
            groups
                .entry((slug.into(), arm.into()))
                .or_default()
                .push(*r);
        }
        let references: BTreeMap<_, _> = measured
            .iter()
            .filter(|r| r["role"] == "reference")
            .map(|r| {
                (
                    (
                        r["workload"]["workload_id"].to_string(),
                        r["round"].to_string(),
                    ),
                    *r,
                )
            })
            .collect();
        for ((curve, arm), r) in groups {
            let mut statuses = BTreeMap::<String, usize>::new();
            let mut unpriced = BTreeSet::new();
            let mut s = Vec::new();
            let mut gae = Vec::new();
            let mut pairs = Vec::new();
            let mut errors = BTreeSet::new();
            let mut phases: BTreeMap<String, Vec<f64>> = BTreeMap::new();
            let mut counters: BTreeMap<String, u128> = BTreeMap::new();
            let mut fbids = BTreeSet::new();
            for row in &r {
                *statuses
                    .entry(
                        row["outcome"]["status"]
                            .as_str()
                            .ok_or("missing status")?
                            .into(),
                    )
                    .or_default() += 1;
                if let Some(e) = row["outcome"]["error"].as_str() {
                    errors.insert(e.to_string());
                }
                for (k, v) in row["counters"].as_object().ok_or("missing counters")? {
                    *counters.entry(k.clone()).or_default() +=
                        v.as_u64().ok_or("non-integer counter")? as u128;
                }
                if let Some(u) = row["cost"]["unpriced"].as_array() {
                    for v in u {
                        unpriced.insert(v.clone().to_string());
                    }
                }
                if let Some(id) = row["factor_base"]["fb_id"].as_str() {
                    fbids.insert(id.to_string());
                }
                if row["outcome"]["matches_target"] == true
                    && row["outcome"]["status"] == "verified"
                {
                    if let Some(v) = row["cost"]["s"].as_f64() {
                        s.push(v);
                    }
                    if let Some(v) = row["cost"]["total_gae"].as_f64() {
                        gae.push(v);
                    }
                    if let Some(ps) = row["phases"].as_array() {
                        for p in ps {
                            if let Some(v) = p["gae"].as_f64() {
                                phases
                                    .entry(p["name"].as_str().unwrap().into())
                                    .or_default()
                                    .push(v);
                            }
                        }
                    }
                    let key = (
                        row["workload"]["workload_id"].to_string(),
                        row["round"].to_string(),
                    );
                    if let Some(reference) = references.get(&key) {
                        if reference["outcome"]["matches_target"] == true
                            && reference["outcome"]["status"] == "verified"
                        {
                            if let (Some(a), Some(b)) = (
                                reference["cost"]["total_gae"].as_f64(),
                                row["cost"]["total_gae"].as_f64(),
                            ) {
                                pairs.push((a, b));
                            }
                        }
                    }
                }
            }
            let ref_sum: f64 = pairs.iter().map(|p| p.0).sum();
            let cand_sum: f64 = pairs.iter().map(|p| p.1).sum();
            let ratio = (pairs.len() == r.len() && !pairs.is_empty() && ref_sum > 0.)
                .then_some(cand_sum / ref_sum);
            let phase_means: BTreeMap<_, _> =
                phases.iter().map(|(k, v)| (k.clone(), mean(v))).collect();
            rows.push(json!({"session":name,"session_id":session["session_id"],"curve":curve,"arm":arm,
                "records":r.len(),"statuses":statuses,"verified":s.len(),"mean_s":mean(&s),
                "min_s":s.iter().copied().reduce(f64::min),"max_s":s.iter().copied().reduce(f64::max),
                "mean_cold_gae":mean(&gae),"paired_complete_counted_ratio_candidate_over_rho":ratio,
                "paired_successes":pairs.len(),"ratio_null_if_censored":true,"unpriced":unpriced,
                "fully_priced_cost_ratio":null,"wall_time_result":null,"errors":errors,
                "phase_mean_gae":phase_means,"counter_sums":counters,"factor_base_ids":fbids,
                "audit_ok":audit["ok"],"scope":"counted-work diagnostic; partial counts are not a bound on full-cost ratios"}));
        }
        all.push(json!({"name":name,"session_id":session["session_id"],"records":records.len(),"measured":measured.len(),
            "audit_ok":audit["ok"],"replays":audit["replays"].as_array().map(Vec::len),
            "verified_records":audit["verified_records"],"problems":audit["problems"]}));
    }
    super::save(
        &root.join("summary.json"),
        &json!({"schema":"ecbench-ten-questions-summary/v1","sessions":all,"rows":rows,
        "wall_claim":false,"end_to_end_speedup_claim":false,"asymptotic_claim":false,
        "new_mathematical_theorem_claim":false,"source":"sealed records and replay receipts only"}),
    )
}

fn mean(xs: &[f64]) -> Option<f64> {
    (!xs.is_empty()).then(|| xs.iter().sum::<f64>() / xs.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn absent_cost_is_missing_not_zero() {
        assert_eq!(mean(&[]), None);
        assert_eq!(mean(&[1., 3.]), Some(2.));
    }
}
