//! Read-only native aggregation and public-target replay input generation.
use super::{write, ROOT};
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, path::Path};

pub fn run() {
    let root = Path::new(ROOT);
    let mut counts = BTreeMap::<String, u64>::new();
    let mut records = Vec::new();
    let mut vectors = BTreeMap::<String, BTreeMap<String, Value>>::new();
    let mut replay = String::new();
    let mut metrics = BTreeMap::<u64, BTreeMap<String, Value>>::new();
    let mut sealed = Vec::new();
    let mut summaries = Vec::new();
    let mut replication_metrics = Vec::new();
    let mut dirs: Vec<_> = fs::read_dir(root.join("sessions"))
        .unwrap()
        .map(|v| v.unwrap().path())
        .filter(|v| v.is_dir())
        .collect();
    dirs.sort();
    for dir in dirs {
        let name = dir.file_name().unwrap().to_str().unwrap();
        if !(name.starts_with("primary-")
            || name.starts_with("rep-")
            || name.starts_with("recurrence-pilot-")
            || name.starts_with("control"))
        {
            continue;
        }
        let file = dir.join("records.jsonl");
        if !file.exists() {
            continue;
        }
        let seal_path = dir.join("session.json");
        if !seal_path.exists() || super::read(&seal_path)["status"] != "complete" {
            continue;
        }
        let seal = super::read(&seal_path);
        let raw = fs::read(&file).unwrap();
        assert_eq!(super::sha256_hex(&raw), seal["records_sha256"]);
        sealed.push(name.to_string());
        let mut arms = BTreeMap::<String, Vec<Value>>::new();
        for line in fs::read_to_string(file).unwrap().lines() {
            let v: Value = serde_json::from_str(line).unwrap();
            arms.entry(v["arm"].as_str().unwrap().to_string())
                .or_default()
                .push(v.clone());
            let cohort = if name.starts_with("primary-") {
                "primary"
            } else if name.starts_with("rep-") {
                "replication"
            } else {
                "control"
            };
            let arm = v["arm"].as_str().unwrap();
            let status = v["outcome"]["status"].as_str().unwrap();
            if cohort == "primary" && arm != "rho" {
                let p = v["workload"]["curve_spec"]["p"].as_u64().unwrap();
                let a = v["workload"]["curve_spec"]["a"].as_u64().unwrap();
                metrics.entry(p).or_default().insert(format!("{a}/{arm}"), json!({"a":a,"arm":arm,"queries":v["detail"]["native"]["queries"],"cpu":v["detail"]["native"]["cpu_seconds"],"status":status}));
            }
            *counts
                .entry(format!("{cohort}/{arm}/{status}"))
                .or_default() += 1;
            if (cohort == "primary" || cohort == "replication")
                && (arm == "velu" || arm == "scalar")
            {
                let key = format!(
                    "{cohort}/{}/{}",
                    v["workload"]["workload_id"].as_str().unwrap(),
                    v["round"]
                );
                vectors.entry(key).or_default().insert(
                    arm.to_string(),
                    json!({
                        "base":v["detail"]["replay"]["base_sha256"],
                        "relations":v["detail"]["relations"],
                        "queries":v["detail"]["native"]["queries"],
                        "rank":v["detail"]["native"]["rank"],
                        "recovered":v["outcome"]["recovered"],
                        "seed":v["algorithm_seed"], "status":status
                    }),
                );
            }
            if cohort == "replication" && arm != "rho" {
                let label = name.strip_prefix("rep-").unwrap();
                let family = label
                    .rsplit_once('-')
                    .map(|(prefix, _)| prefix)
                    .unwrap_or(label);
                replication_metrics.push(json!({"group":family,"session":name,"workload":v["workload"]["workload_id"],"p":v["workload"]["curve_spec"]["p"],"a":v["workload"]["curve_spec"]["a"],"round":v["round"],"target_index":v["workload"]["target_index"],"arm":arm,"status":status,"queries":v["detail"]["native"]["queries"],"native_backend_cpu_seconds":v["detail"]["native"]["cpu_seconds"],"process_cpu_seconds":(v["time"]["user_ns"].as_u64().unwrap() as f64+v["time"]["sys_ns"].as_u64().unwrap() as f64)/1e9}));
            }
            if status == "verified" {
                let c = &v["workload"]["curve_spec"];
                let d = v["outcome"]["recovered"].as_str().unwrap();
                let target = &v["workload"]["target"];
                let qx =
                    u64::from_str_radix(target[0].as_str().unwrap().trim_start_matches("0x"), 16)
                        .unwrap();
                let qy =
                    u64::from_str_radix(target[1].as_str().unwrap().trim_start_matches("0x"), 16)
                        .unwrap();
                replay += &format!(
                    "[{},{},{},{},[{gox},{goy}],[{qx},{qy}],{}]\n",
                    c["p"],
                    c["a"],
                    c["b"],
                    c["r"],
                    d,
                    gox = c["gx"],
                    goy = c["gy"],
                    qx = qx,
                    qy = qy
                );
            }
            records.push(json!({"session":name,"record_id":v["record_id"],"arm":arm,"status":status,
                "workload":v["workload"]["workload_id"],"native_cpu_seconds":v["detail"]["native"]["cpu_seconds"],
                "native_queries":v["detail"]["native"]["queries"],"online":v["online"],"process_cpu_user_ns":v["time"]["user_ns"],
                "complete_common_price":false,"isolation":v["isolation"]["level"]}));
        }
        let mut cells = Vec::new();
        for (arm, rows) in arms {
            let mut statuses = BTreeMap::<String, u64>::new();
            let mut cpu = 0u128;
            let mut phases = BTreeMap::<String, (u128, u64)>::new();
            for row in &rows {
                *statuses
                    .entry(row["outcome"]["status"].as_str().unwrap().to_string())
                    .or_default() += 1;
                cpu += row["time"]["user_ns"].as_u64().unwrap() as u128
                    + row["time"]["sys_ns"].as_u64().unwrap() as u128;
                if let Some(values) = row["online"]["phases_ns"].as_object() {
                    for (k, v) in values {
                        let pair = phases.entry(k.clone()).or_default();
                        pair.0 += v.as_u64().unwrap() as u128;
                        pair.1 += 1;
                    }
                }
            }
            let phases:Vec<_>=phases.into_iter().map(|(name,(total,n))|json!({"phase":name,"observations":n,"mean_wall_seconds":total as f64/n as f64/1e9})).collect();
            cells.push(json!({"arm":arm,"attempts":rows.len(),"statuses":statuses,"mean_recorded_process_cpu_seconds":cpu as f64/rows.len() as f64/1e9,"online_phases":phases}));
        }
        summaries.push(json!({"session":name,"records_sha256":seal["records_sha256"],"binary_sha256":seal["binary_sha256"],"cells":cells,"isolation":"L0","common_price":"incomplete; no admitted IC/rho ratio","scope":"All attempts retained; process CPU means include every status. Phase means include only available online receipts; their observation counts are explicit."}));
    }
    let mut matched = 0;
    let mut replication_matched = 0;
    for (key, pair) in &vectors {
        if let (Some(a), Some(b)) = (pair.get("velu"), pair.get("scalar")) {
            assert_eq!(a, b, "same-crater scientific output mismatch: {key}");
            if key.starts_with("primary/") {
                matched += 1;
            } else {
                replication_matched += 1;
            }
        }
    }
    write(
        root.join("results/campaign-progress.json"),
        &json!({"schema":"conductor-gap.progress/v1",
        "scope":"snapshot of complete sessions only; each record-file hash checked against its seal",
        "counts":counts,"same_crater_exact_solve_matches":matched,"replication_same_crater_exact_solve_matches":replication_matched,"session_json_present":sealed,
        "records":records,"interpretation":"L0; unpriced work; no admitted IC/rho speedup or hardness claim"}),
    );
    let mut paired = Vec::new();
    for (p, entries) in metrics {
        let Some(crater) = entries.values().find(|v| v["arm"] == "velu") else {
            continue;
        };
        let Some(floor) = entries
            .values()
            .find(|v| v["a"] != crater["a"] && v["arm"] == "scalar")
        else {
            continue;
        };
        paired.push(json!({"p":p,"crater":crater,"floor":floor}));
    }
    let complete = paired.len() == 100
        && paired
            .iter()
            .all(|v| v["crater"]["status"] == "verified" && v["floor"]["status"] == "verified");
    let mut statistics = json!({"paired_endpoints_observed":paired.len(),"all_100_pairs_verified":complete,"paired":paired,"bootstrap_unit":"whole conductor-gap pair; 10000 resamples; ratio of arithmetic means","scope":"L0 native CPU/query diagnostics, no common-price IC/rho comparison"});
    if complete {
        for metric in ["queries", "cpu"] {
            let xs: Vec<_> = paired
                .iter()
                .map(|v| {
                    (
                        v["crater"][metric].as_f64().unwrap(),
                        v["floor"][metric].as_f64().unwrap(),
                    )
                })
                .collect();
            let ratio = |indices: &[usize]| {
                indices.iter().map(|&i| xs[i].1).sum::<f64>()
                    / indices.iter().map(|&i| xs[i].0).sum::<f64>()
            };
            let mut seed = 20261018u64;
            let mut draws = Vec::new();
            for _ in 0..10000 {
                let mut indices = Vec::new();
                for _ in 0..xs.len() {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    indices.push((seed % xs.len() as u64) as usize);
                }
                draws.push(ratio(&indices));
            }
            draws.sort_by(f64::total_cmp);
            statistics[format!("floor_over_crater_{metric}")] = json!({"ratio_of_means":ratio(&(0..xs.len()).collect::<Vec<_>>()),"paired_percentile_95":[draws[250],draws[9749]]});
        }
    }
    write(root.join("results/paired-primary.json"), &statistics);
    write(
        root.join("results/session-summaries.json"),
        &json!({"schema":"conductor-gap.session-summaries/v1","sessions":summaries}),
    );
    let replication_comparisons = replication_comparisons(&replication_metrics);
    write(
        root.join("results/replication-metrics.json"),
        &json!({"schema":"conductor-gap.replication-metrics/v1","scope":"All 155 IC replication attempts, retained individually; same-crater exact matches checked separately. L0 CPU and query diagnostics; no common-price IC/rho claim.","comparisons":replication_comparisons,"attempts":replication_metrics}),
    );
    fs::write(root.join("prepared/target-replay.gpvec"), replay).unwrap();
}

fn replication_comparisons(rows: &[Value]) -> Vec<Value> {
    let mut out = Vec::new();
    for family in ["fresh", "isomorphic", "d43", "pilot-five-rounds"] {
        let crater_group = if family == "pilot-five-rounds" {
            family.to_string()
        } else {
            format!("{family}-crater")
        };
        let floor_group = if family == "pilot-five-rounds" {
            family.to_string()
        } else {
            format!("{family}-floor")
        };
        let mut pairs = Vec::new();
        for root in rows
            .iter()
            .filter(|v| v["group"] == crater_group && v["arm"] == "velu")
        {
            let floor: Vec<_> = rows
                .iter()
                .filter(|v| {
                    v["group"] == floor_group
                        && v["arm"] == "scalar"
                        && v["p"] == root["p"]
                        && v["a"] != root["a"]
                        && v["round"] == root["round"]
                        && v["target_index"] == root["target_index"]
                })
                .collect();
            assert_eq!(floor.len(), 1, "replication pairing missing/ambiguous");
            assert_eq!(root["status"], "verified");
            assert_eq!(floor[0]["status"], "verified");
            pairs.push((root, floor[0]));
        }
        let mut clusters = BTreeMap::<u64, Vec<usize>>::new();
        for (i, (root, _)) in pairs.iter().enumerate() {
            clusters
                .entry(root["p"].as_u64().unwrap())
                .or_default()
                .push(i);
        }
        let mut entry = json!({"family":family,"endpoint_pairs":pairs.len(),"curve_gap_clusters":clusters.len(),"bootstrap":"10000 whole-gap cluster resamples; all targets/rounds of each sampled gap retained","scope":"L0 diagnostics; only one curve pair in the query-seed control, so no curve-population interval there"});
        let clusters: Vec<_> = clusters.into_values().collect();
        for metric in [
            "queries",
            "native_backend_cpu_seconds",
            "process_cpu_seconds",
        ] {
            let ratio = |indices: &[usize]| {
                indices
                    .iter()
                    .map(|&i| pairs[i].1[metric].as_f64().unwrap())
                    .sum::<f64>()
                    / indices
                        .iter()
                        .map(|&i| pairs[i].0[metric].as_f64().unwrap())
                        .sum::<f64>()
            };
            let estimate = ratio(&(0..pairs.len()).collect::<Vec<_>>());
            let mut result =
                json!({"floor_over_crater_ratio_of_means":estimate,"cluster_percentile_95":null});
            if clusters.len() > 1 {
                let mut seed = 20261019u64;
                let mut draws = Vec::new();
                for _ in 0..10000 {
                    let mut indices = Vec::new();
                    for _ in 0..clusters.len() {
                        seed ^= seed << 13;
                        seed ^= seed >> 7;
                        seed ^= seed << 17;
                        indices
                            .extend_from_slice(&clusters[(seed % clusters.len() as u64) as usize]);
                    }
                    draws.push(ratio(&indices));
                }
                draws.sort_by(f64::total_cmp);
                result["cluster_percentile_95"] = json!([draws[250], draws[9749]]);
            }
            entry[metric] = result;
        }
        out.push(entry);
    }
    out
}
