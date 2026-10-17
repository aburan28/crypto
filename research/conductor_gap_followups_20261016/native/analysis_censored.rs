use super::*;
use rand::Rng;
use std::collections::BTreeMap;
fn rows() -> Vec<Value> {
    let root = Path::new(ROOT);
    let mut rows = Vec::new();
    for cohort in [
        "crossover",
        "yield",
        "rho",
        "rho-coordinate",
        "orbit",
        "dose",
        "overlap",
    ] {
        let p = root.join(format!("prepared/{cohort}.json"));
        if !p.exists() {
            continue;
        }
        let jobs = read(p);
        for (index, j) in jobs.as_array().unwrap().iter().enumerate() {
            let label = format!("{cohort}-{index:04}");
            let receipt = root.join(format!("results/receipts/{label}.json"));
            if !receipt.exists() {
                continue;
            }
            let seal = read(receipt);
            assert_eq!(
                seal["input_sha256"],
                sha256_hex(serde_json::to_string(j).unwrap().as_bytes())
            );
            let dir = root.join("runs").join(&label);
            let mut row = json!({"cohort":cohort,"index":index,"pair":j["pair"],"rep":j["rep"],"arm":j["arm"],"curve":j["curve"],"partition_seed":j["partition_seed"],"completed":seal["completed"],"exit_code":seal["exit_code"],"selection":j["selection"]});
            if seal["completed"] != true {
                rows.push(row);
                continue;
            }
            let bytes = fs::read(dir.join("result.json")).unwrap();
            assert_eq!(sha256_hex(&bytes), seal["result_sha256"]);
            let v: Value = serde_json::from_slice(&bytes).unwrap();
            let input = fs::read(dir.join("input.json")).unwrap();
            assert_eq!(sha256_hex(&input), v["input_sha256"]);
            assert_eq!(serde_json::from_slice::<Value>(&input).unwrap(), *j);
            let time = seal["process_time"].as_str().unwrap();
            let time: Value =
                serde_json::from_str(time.lines().find(|s| s.starts_with('{')).unwrap()).unwrap();
            row["process_cpu_seconds"] = json!(
                time["user_seconds"].as_f64().unwrap() + time["sys_seconds"].as_f64().unwrap()
            );
            row["process_wall_seconds"] = time["wall_seconds"].clone();
            row["max_rss_kib"] = time["max_rss_kib"].clone();
            if let Some(n) = v["result"].get("native") {
                assert_eq!(read(dir.join("stdout.json")), *n);
                let ev = fs::read(dir.join("events.jsonl")).unwrap();
                assert_eq!(sha256_hex(&ev), v["result"]["replay"]["events_sha256"]);
                row["native"] = n.clone();
                row["base_sha256"] = v["result"]["replay"]["base_sha256"].clone();
                row["events_sha256"] = v["result"]["replay"]["events_sha256"].clone();
                row["verified"] = n["verified"].clone();
                row["queries"] = n["queries"].clone();
                row["query_cap"] = j["query_limit"].clone();
                row["observed_query_lower"] = n["queries"].clone();
                row["restricted_query_upper"] =
                    if n["verified"] == true || n["status"] == "QUERY_LIMIT" {
                        n["queries"].clone()
                    } else {
                        j["query_limit"].clone()
                    };
                row["native_cpu_seconds"] = n["cpu_seconds"].clone();
            } else {
                let n = &v["result"]["result"];
                row["verified"] = n["verified"].clone();
                row["rho"] = v["result"].clone();
                let ops = &n["group_ops"];
                let adds = u(ops, "adds") + u(ops, "doubles");
                row["group_addition_equivalents"] = json!(adds);
                row["group_ops_per_sqrt_r"] =
                    json!(adds as f64 / (u(&j["curve"], "N") as f64).sqrt());
                row["process_cpu_per_group_op"] =
                    json!(row["process_cpu_seconds"].as_f64().unwrap() / adds as f64);
            }
            rows.push(row);
        }
    }
    rows
}
fn grouped(rows: &[Value]) -> Value {
    let mut groups = BTreeMap::<String, Vec<&Value>>::new();
    for row in rows {
        groups
            .entry(format!(
                "{}/{}",
                row["cohort"].as_str().unwrap(),
                row["arm"].as_str().unwrap()
            ))
            .or_default()
            .push(row);
    }
    let mut out = BTreeMap::new();
    for (name, rs) in groups {
        let complete: Vec<_> = rs
            .iter()
            .filter(|r| r["completed"] == true)
            .copied()
            .collect();
        let mean = |key: &str| {
            let v: Vec<f64> = complete.iter().filter_map(|r| r[key].as_f64()).collect();
            if v.is_empty() {
                Value::Null
            } else {
                json!(v.iter().sum::<f64>() / v.len() as f64)
            }
        };
        let mut statuses = BTreeMap::<String, usize>::new();
        for r in &rs {
            let status = if r["completed"] != true {
                "execution-failure"
            } else if r["native"].is_object() {
                r["native"]["status"].as_str().unwrap()
            } else if r["verified"] == true {
                "PASS"
            } else {
                "OPERATION_LIMIT"
            };
            *statuses.entry(status.into()).or_default() += 1;
        }
        out.insert(name,json!({"attempts":rs.len(),"completed":complete.len(),"verified":rs.iter().filter(|r|r["verified"]==true).count(),"statuses":statuses,"mean_process_cpu_seconds":mean("process_cpu_seconds"),"mean_native_cpu_seconds":mean("native_cpu_seconds"),"mean_observed_queries":mean("queries"),"restricted_mean_query_lower":mean("observed_query_lower"),"restricted_mean_query_upper":mean("restricted_query_upper"),"mean_group_ops_per_sqrt_r":mean("group_ops_per_sqrt_r"),"mean_process_cpu_per_group_op":mean("process_cpu_per_group_op")}));
    }
    json!(out)
}
fn contrast(
    rows: &[Value],
    cohort: &str,
    numerator: &str,
    denominator: &str,
    metric: &str,
) -> Value {
    if rows.iter().any(|r| {
        r["cohort"] == cohort
            && (r["arm"] == numerator || r["arm"] == denominator)
            && r["completed"] != true
    }) {
        return json!({"blocked":"execution failures in the comparison panel; do not compute a successful-only ratio","metric":metric});
    }
    let mut pairs = BTreeMap::<u64, (Vec<f64>, Vec<f64>)>::new();
    for r in rows {
        if r["cohort"] != cohort || r["completed"] != true {
            continue;
        }
        let Some(value) = r[metric].as_f64() else {
            continue;
        };
        let entry = pairs.entry(u(r, "pair")).or_default();
        if r["arm"] == numerator {
            entry.0.push(value)
        } else if r["arm"] == denominator {
            entry.1.push(value)
        }
    }
    let pairs: Vec<_> = pairs
        .into_values()
        .filter(|(a, b)| !a.is_empty() && !b.is_empty())
        .map(|(a, b)| {
            (
                a.iter().sum::<f64>() / a.len() as f64,
                b.iter().sum::<f64>() / b.len() as f64,
            )
        })
        .collect();
    if pairs.is_empty() {
        return Value::Null;
    }
    let ratio = pairs.iter().map(|x| x.0).sum::<f64>() / pairs.iter().map(|x| x.1).sum::<f64>();
    let mut rng = StdRng::seed_from_u64(2026101602);
    let mut draws = Vec::new();
    for _ in 0..10000 {
        let (mut a, mut b) = (0.0, 0.0);
        for _ in &pairs {
            let p = pairs[rng.gen_range(0..pairs.len())];
            a += p.0;
            b += p.1;
        }
        draws.push(a / b)
    }
    draws.sort_by(f64::total_cmp);
    json!({"clusters":pairs.len(),"numerator":numerator,"denominator":denominator,"metric":metric,"ratio_of_cluster_means":ratio,"bootstrap_95":[draws[249],draws[9749]],"draws":10000,"seed":2026101602u64,"interpretation":"Observed budget-limited attempt work; inspect censors and execution failures before interpreting as solve work. No successful-only speed ratio."})
}
pub fn progress() {
    let r = rows();
    write(
        Path::new(ROOT).join("results/progress.json"),
        &json!({"schema":"conductor-followup.progress/v1","sealed_attempts":r.len(),"groups":grouped(&r),"complete_primary_expected":2160,"additional_overlap_expected":40,"additional_coordinates_expected":16}),
    );
    println!("{} sealed attempts", r.len());
}
pub fn analyze() {
    let root = Path::new(ROOT);
    assert!(
        !root.join("results/campaign.lock").exists(),
        "do not bootstrap while the timed campaign is active"
    );
    let r = rows();
    assert_eq!(
        r.len(),
        2216,
        "analyze requires every scheduled receipt; use progress for an incomplete campaign"
    );
    let mut contrasts = BTreeMap::new();
    for (cohort, num, den, metrics) in [
        (
            "crossover",
            "floor-recurrence",
            "crater-recurrence",
            vec![
                "process_cpu_seconds",
                "observed_query_lower",
                "restricted_query_upper",
                "native_cpu_seconds",
            ],
        ),
        (
            "crossover",
            "crater-cm",
            "crater-recurrence",
            vec![
                "process_cpu_seconds",
                "observed_query_lower",
                "native_cpu_seconds",
            ],
        ),
        (
            "yield",
            "floor-recurrence",
            "crater-recurrence",
            vec![
                "observed_query_lower",
                "restricted_query_upper",
                "process_cpu_seconds",
            ],
        ),
        (
            "rho",
            "floor-rho",
            "crater-rho",
            vec![
                "group_ops_per_sqrt_r",
                "process_cpu_seconds",
                "process_cpu_per_group_op",
            ],
        ),
    ] {
        for metric in metrics {
            contrasts.insert(
                format!("{cohort}/{num}/{metric}"),
                contrast(&r, cohort, num, den, metric),
            );
        }
    }
    let mut exact = Vec::new();
    for pair in 0..100 {
        for rep in 0..4 {
            let cm = r.iter().find(|x| {
                x["cohort"] == "crossover"
                    && x["pair"] == pair
                    && x["rep"] == rep
                    && x["arm"] == "crater-cm"
                    && x["completed"] == true
            });
            let re = r.iter().find(|x| {
                x["cohort"] == "crossover"
                    && x["pair"] == pair
                    && x["rep"] == rep
                    && x["arm"] == "crater-recurrence"
                    && x["completed"] == true
            });
            if let (Some(a), Some(b)) = (cm, re) {
                assert_eq!(a["base_sha256"], b["base_sha256"]);
                let load = |row: &Value| -> Vec<Value> {
                    fs::read_to_string(root.join(format!(
                        "runs/crossover-{:04}/events.jsonl",
                        u(row, "index")
                    )))
                    .unwrap()
                    .lines()
                    .map(|s| serde_json::from_str(s).unwrap())
                    .collect()
                };
                let ae = load(a);
                let be = load(b);
                let receipt = compare_prefix(a, b, &ae, &be);
                exact.push(json!({"pair":pair,"rep":rep,"ok":true,"comparison":receipt}));
            }
        }
    }
    let mut audit = Vec::new();
    for row in &r {
        if row["cohort"] == "orbit" && row["completed"] == true {
            let n = &row["native"];
            let b = u(n, "B");
            audit.push(json!({"pair":row["pair"],"conductor":row["curve"]["conductor"],"points":b,"distinct_pair_abscissae":n["pair_table_keys"],"pair_abscissa_collisions":b*b-u(n,"pair_table_keys"),"signed_point_collisions":0,"replay":"all points independently checked; repeated abscissae cause a preserved failed attempt"}));
        }
    }
    let mut ranks = Vec::new();
    for row in &r {
        if row["cohort"] != "yield" || row["completed"] != true {
            continue;
        }
        let dir = root.join(format!("runs/yield-{:04}", u(row, "index")));
        let events: Vec<Value> = fs::read_to_string(dir.join("events.jsonl"))
            .unwrap()
            .lines()
            .map(|s| serde_json::from_str(s).unwrap())
            .collect();
        for cutoff in [1000000u64, 10000000, 30000000, 100000000, 300000000] {
            let exposed = row["verified"] == true || u(row, "queries") >= cutoff;
            let rank = events
                .iter()
                .filter(|e| u(e, "query") < cutoff)
                .map(|e| u(e, "rank"))
                .max()
                .unwrap_or(0);
            ranks.push(json!({"pair":row["pair"],"rep":row["rep"],"arm":row["arm"],"query_cutoff":cutoff,"rank":if exposed{json!(rank)}else{Value::Null},"exposed_or_already_solved":exposed}));
        }
    }
    let dose: Vec<_> = r
        .iter()
        .filter(|x| x["cohort"] == "dose")
        .cloned()
        .collect();
    for d in [-19, -43] {
        let family: Vec<_> = dose
            .iter()
            .filter(|x| x["curve"]["D"] == d)
            .map(|x| {
                let mut v = x.clone();
                v["pair"] = v["curve"]["class"].clone();
                let method = if v["arm"] == "dose-rho" { "rho" } else { "ic" };
                v["arm"] = json!(format!("{method}-stratum-{}", v["curve"]["stratum"]));
                v
            })
            .collect();
        for method in ["ic", "rho"] {
            for stratum in 2..=4 {
                for metric in if method == "ic" {
                    vec![
                        "process_cpu_seconds",
                        "observed_query_lower",
                        "restricted_query_upper",
                    ]
                } else {
                    vec!["process_cpu_seconds", "group_ops_per_sqrt_r"]
                } {
                    let num = format!("{method}-stratum-{stratum}");
                    let den = format!("{method}-stratum-1");
                    contrasts.insert(
                        format!("dose/D{d}/{num}/{metric}"),
                        contrast(&family, "dose", &num, &den, metric),
                    );
                }
            }
        }
    }
    write(
        root.join("results/analysis.json"),
        &json!({"schema":"conductor-followup.analysis/v1","sealed_attempts":r.len(),"expected_total":2216,"groups":grouped(&r),"contrasts":contrasts,"same_crater_checks":exact,"orbit_audit":audit,"yield_rank_snapshots":ranks,"dose_records":dose,"records":r,"scope":"53-bit synthetic public targets; no general hardness or fully priced IC/rho speedup claim. Bootstrap intervals are descriptive and unadjusted for multiple comparisons."}),
    );
}
