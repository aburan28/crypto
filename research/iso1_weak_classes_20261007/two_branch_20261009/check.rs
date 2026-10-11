//! Compare exact two-branch class support and construction receipts without censoring failures.
use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, path::Path};
fn fields(s: &str) -> BTreeMap<&str, &str> {
    s.split('|')
        .skip(1)
        .filter_map(|f| f.split_once('='))
        .collect()
}
fn num(m: &BTreeMap<&str, &str>, k: &str) -> u64 {
    m[k].parse().unwrap()
}
fn verify_receipt(dir: &Path) -> String {
    let r: Value = serde_json::from_slice(&fs::read(dir.join("receipt.json")).unwrap()).unwrap();
    assert_eq!(r["status"], "PASS");
    for file in ["stdout", "stderr"] {
        assert_eq!(
            r[format!("{file}_sha256")],
            hex::encode(sha256(&fs::read(dir.join(format!("{file}.txt"))).unwrap()))
        );
    }
    assert!(fs::read(dir.join("stderr.txt")).unwrap().is_empty());
    let freeze: Value =
        serde_json::from_slice(&fs::read(dir.join("source_freeze.json")).unwrap()).unwrap();
    assert_eq!(
        freeze["script_sha256"],
        hex::encode(sha256(&fs::read(dir.join("executed_script.gp")).unwrap()))
    );
    if dir.join("executed_driver.txt").exists() {
        assert_eq!(
            freeze["driver_source_sha256"],
            hex::encode(sha256(&fs::read(dir.join("executed_driver.txt")).unwrap()))
        );
    } else {
        assert_eq!(
            hex::encode(sha256(
                &fs::read(dir.join("executed_driver_v1.txt")).unwrap()
            )),
            "4e97c67f96d8c2adc12280a539ec086853e57a777784bc65eb39bedbe50feb87"
        );
    }
    fs::read_to_string(dir.join("stdout.txt")).unwrap()
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 1);
    let base = Path::new(&args[0]);
    let old = fs::read_to_string(base.join("../p7_exact_absolute.csv")).unwrap();
    let expected: BTreeMap<i64, bool> = old
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<_> = l.split(',').collect();
            (c[2].parse().unwrap(), c[3].parse::<u64>().unwrap() > 0)
        })
        .collect();
    let mut class_tables = BTreeMap::new();
    let mut summaries = Vec::new();
    for p in [7u64, 11, 13, 17, 19, 29, 37] {
        let dir = base.join(format!("p{p}_run1"));
        if !dir.join("receipt.json").exists() {
            continue;
        }
        let raw = verify_receipt(&dir);
        if [13, 37].contains(&p) {
            let old = fs::read_to_string(base.join(format!("../p{p}_twist_derived.csv"))).unwrap();
            let mut weighted = BTreeMap::<i64, u64>::new();
            for line in raw.lines().filter(|l| l.starts_with("PLUS|")) {
                let m = fields(line);
                let t = m["trace"].parse::<i64>().unwrap().abs();
                *weighted.entry(t).or_default() += num(&m, "orbit_size");
            }
            for line in old.lines().skip(1) {
                let c: Vec<_> = line.split(',').collect();
                let t: i64 = c[2].parse().unwrap();
                if t % p as i64 != 0 {
                    assert_eq!(
                        c[3].parse::<u64>().unwrap(),
                        weighted.get(&t.abs()).copied().unwrap_or(0)
                    );
                }
            }
        }
        let mut table = BTreeMap::new();
        let mut summary = None;
        for line in raw.lines() {
            let m = fields(line);
            if line.starts_with("CLASS|") {
                let t: i64 = m["trace"].parse().unwrap();
                let a = m["plus"] == "1";
                let b = m["minus"] == "1";
                assert_eq!(m["combined"] == "1", a || b);
                if p == 7 {
                    assert_eq!(a, expected[&t]);
                }
                table.insert(t, (a, b));
            }
            if line.starts_with("TWO_BRANCH_COMPLETE|") {
                summary = Some(
                    m.iter()
                        .map(|(k, v)| ((*k).to_string(), (*v).to_string()))
                        .collect::<BTreeMap<_, _>>(),
                );
            }
        }
        let s = summary.unwrap();
        assert_eq!(table.len(), s["ordinary_classes"].parse::<usize>().unwrap());
        assert_eq!(
            table.values().filter(|v| v.0).count(),
            s["plus_support"].parse::<usize>().unwrap()
        );
        assert_eq!(
            table.values().filter(|v| v.1).count(),
            s["minus_support"].parse::<usize>().unwrap()
        );
        assert_eq!(
            table.values().filter(|v| v.0 && v.1).count(),
            s["overlap"].parse::<usize>().unwrap()
        );
        assert_eq!(
            table.values().filter(|v| !v.0 && !v.1).count(),
            s["combined_zero"].parse::<usize>().unwrap()
        );
        let mut csv = String::from("p,trace,cubic,quadratic,combined\n");
        for (t, (a, b)) in &table {
            csv.push_str(&format!(
                "{p},{t},{},{},{}\n",
                u8::from(*a),
                u8::from(*b),
                u8::from(*a || *b)
            ));
        }
        fs::write(base.join(format!("p{p}_classes.csv")), csv).unwrap();
        summaries.push(json!(s));
        class_tables.insert(p, table);
    }
    let raw = verify_receipt(&base.join("construction_p7_run1"));
    let mut seed_traces = BTreeMap::new();
    let mut costs = Vec::new();
    let mut plus = 0;
    let mut minus = 0;
    let mut zeros = 0;
    let mut missed = 0;
    let mut routes = 0;
    for line in raw.lines() {
        let m = fields(line);
        if line.starts_with("SOURCE|") {
            seed_traces.insert(num(&m, "seed"), m["trace"].parse::<i64>().unwrap());
        }
        if line.starts_with("CONSTRUCTION|") {
            let t = seed_traces[&num(&m, "seed")];
            let label = class_tables[&7][&t];
            if m["status"] == "VERIFIED_WITNESS" {
                assert!(label.0 || label.1);
                if m["branch"] == "PLUS" {
                    assert!(label.0);
                    plus += 1;
                } else {
                    assert_eq!(m["branch"], "MINUS");
                    assert!(label.1);
                    minus += 1;
                }
                routes += num(&m, "route_length");
            } else if label.0 || label.1 {
                missed += 1;
            } else {
                zeros += 1;
            }
            costs.push(json!(m
                .iter()
                .map(|(k, v)| ((*k).to_string(), (*v).to_string()))
                .collect::<BTreeMap<_, _>>()));
        }
    }
    assert_eq!(costs.len(), 64);
    assert_eq!(seed_traces.len(), 64);
    let mut followed = BTreeMap::new();
    let mut attempts = Vec::new();
    for run in [
        "construction_p7_run1",
        "construction_p7_odd_run3",
        "construction_p7_followup_run1",
        "construction_p7_targeted_run1",
        "construction_p7_degree19_run1",
    ] {
        let dir = base.join(run);
        if !dir.join("receipt.json").exists() {
            continue;
        }
        let text = verify_receipt(&dir);
        for line in text.lines().filter(|l| l.starts_with("CONSTRUCTION|")) {
            let m = fields(line);
            let seed = num(&m, "seed");
            let t = seed_traces[&seed];
            let label = class_tables[&7][&t];
            if m["status"] == "VERIFIED_WITNESS" {
                assert!(label.0 || label.1);
                followed.insert(seed, m["branch"].to_string());
            }
            attempts.push(json!({"run":run,"trace":t,"cost":m.iter().map(|(k,v)|((*k).to_string(),(*v).to_string())).collect::<BTreeMap<_,_>>()}));
        }
    }
    let mut large_rows = Vec::new();
    let mut large_counts = BTreeMap::new();
    let batch = base.join("construction_large_run1");
    if batch.join("batch_receipt.json").exists() {
        let receipt: Value =
            serde_json::from_slice(&fs::read(batch.join("batch_receipt.json")).unwrap()).unwrap();
        assert_eq!(receipt["completed"], 384);
        assert_eq!(receipt["failed"], 0);
        for i in 1..=384 {
            let raw = verify_receipt(&batch.join(format!("source{i:04}")));
            let mut input = None;
            for line in raw.lines() {
                let m = fields(line);
                if line.starts_with("SOURCE|") {
                    input = Some(
                        m.iter()
                            .map(|(k, v)| ((*k).to_string(), (*v).to_string()))
                            .collect::<BTreeMap<_, _>>(),
                    );
                }
                if line.starts_with("CONSTRUCTION|") {
                    let key = (m["p"].to_string(), m["status"].to_string());
                    *large_counts.entry(key).or_insert(0usize) += 1;
                    large_rows.push(json!({"source":input.take().unwrap(),"cost":m.iter().map(|(k,v)|((*k).to_string(),(*v).to_string())).collect::<BTreeMap<_,_>>()}));
                }
            }
        }
        assert_eq!(large_rows.len(), 384);
    }
    let general = base.join("construction_general_run1");
    if general.join("batch_receipt.json").exists() {
        let receipt: Value =
            serde_json::from_slice(&fs::read(general.join("batch_receipt.json")).unwrap()).unwrap();
        assert_eq!(receipt["completed"], 128);
        assert_eq!(receipt["failed"], 0);
        for i in 1..=128 {
            let raw = verify_receipt(&general.join(format!("source{i:04}")));
            let mut input = None;
            for line in raw.lines() {
                let m = fields(line);
                if line.starts_with("SOURCE|") {
                    input = Some(
                        m.iter()
                            .map(|(k, v)| ((*k).to_string(), (*v).to_string()))
                            .collect::<BTreeMap<_, _>>(),
                    );
                }
                if line.starts_with("CONSTRUCTION|") {
                    let key = (m["p"].to_string(), m["status"].to_string());
                    *large_counts.entry(key).or_insert(0usize) += 1;
                    large_rows.push(json!({"source":input.take().unwrap(),"cost":m.iter().map(|(k,v)|((*k).to_string(),(*v).to_string())).collect::<BTreeMap<_,_>>()}));
                }
            }
        }
        assert_eq!(large_rows.len(), 512);
    }
    let cm = verify_receipt(&base.join("cm_controls_run1"));
    let micro = verify_receipt(&base.join("cm_microfactor_run1"));
    let degree6 = verify_receipt(&base.join("cm_degree6_run1"));
    let mut expected_cm = BTreeMap::new();
    for line in cm.lines().filter(|l| l.starts_with("CM_CLASS|")) {
        let m = fields(line);
        expected_cm.insert(num(&m, "trace"), num(&m, "minus_parameters"));
        assert_eq!(
            num(&m, "minus_parameters") > 0,
            class_tables[&7][&(num(&m, "trace") as i64)].1
        );
    }
    for output in [&micro, &degree6] {
        let rows: Vec<_> = output
            .lines()
            .filter(|l| l.starts_with("MICROFACTOR_CLASS|"))
            .collect();
        assert_eq!(rows.len(), expected_cm.len());
        for line in rows {
            let m = fields(line);
            assert_eq!(num(&m, "minus_parameters"), expected_cm[&num(&m, "trace")]);
        }
    }
    let large_counts: Vec<_> = large_counts
        .into_iter()
        .map(|((p, status), count)| json!({"p":p,"status":status,"count":count}))
        .collect();
    let summary = json!({"schema":"iso1.two_branch_summary/v1","censuses":summaries,"p7_construction":{"sources":64,"cubic_witnesses":plus,"quadratic_witnesses":minus,"exact_zero_classes":zeros,"exact_positive_missed":missed,"verified_route_edges":routes,"cost_rows":costs,"all_attempts":attempts,"constructed_after_followups":followed.len(),"final_sources":followed},"large_counts":large_counts,"large_rows":large_rows,"timings":"shared-host resource diagnostics"});
    fs::write(
        base.join("summary.json"),
        serde_json::to_vec_pretty(&summary).unwrap(),
    )
    .unwrap();
    println!("{}",serde_json::to_string(&json!({"censuses":summary["censuses"],"p7":{"plus":plus,"minus":minus,"zero":zeros,"missed":missed,"route_edges":routes}})).unwrap());
}
