//! Conservative floor-distance bounds; never changes a frozen graph database.
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use rusqlite::{Connection, OpenFlags};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, VecDeque},
    env, fs,
    path::Path,
};

fn distances(adj: &[Vec<usize>], seeds: impl Iterator<Item = usize>) -> Vec<Option<u32>> {
    let mut out = vec![None; adj.len()];
    let mut queue = VecDeque::new();
    for i in seeds {
        out[i] = Some(0);
        queue.push_back(i);
    }
    while let Some(i) = queue.pop_front() {
        let d = out[i].unwrap() + 1;
        for &j in &adj[i] {
            if out[j].is_none() {
                out[j] = Some(d);
                queue.push_back(j);
            }
        }
    }
    out
}

fn certify(path: &Path, cutoff: u64) -> Value {
    let sha = sha256_hex(&fs::read(path).unwrap());
    let db = Connection::open_with_flags(
        format!(
            "file:{}?immutable=1",
            path.canonicalize().unwrap().display()
        ),
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_URI,
    )
    .unwrap();
    let meta = |key: &str| -> Value {
        let v: String = db
            .query_row("SELECT value FROM metadata WHERE key=?", [key], |r| {
                r.get(0)
            })
            .unwrap();
        serde_json::from_str(&v).unwrap()
    };
    let plan = meta("plan");
    let replay = meta("verification");
    assert_eq!(replay["status"], "REPLAYED");
    assert!(plan["fundamental_discriminant"].as_i64().unwrap() < -4);
    let mut st = db
        .prepare("SELECT id,expanded,conductor FROM nodes ORDER BY id")
        .unwrap();
    let nodes: Vec<(u64, bool, Option<u64>)> = st
        .query_map([], |r| Ok((r.get(0)?, r.get::<_, u64>(1)? == 1, r.get(2)?)))
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    for (i, (id, _, _)) in nodes.iter().enumerate() {
        assert_eq!(*id, i as u64);
    }
    let known: Vec<_> = nodes
        .iter()
        .map(|(id, expanded, _)| *expanded && *id < cutoff)
        .collect();
    let mut exact = vec![Some(1u64); nodes.len()];
    let mut heights: Vec<Vec<Value>> = vec![vec![]; nodes.len()];
    let mut bounds_checked = 0;
    let mut labels_checked = 0;
    for pair in plan["conductor_factors"].as_array().unwrap() {
        let ell = pair[0].as_u64().unwrap();
        let depth = pair[1].as_u64().unwrap() as u32;
        assert!(depth > 0);
        assert!(plan["selected_degrees"]
            .as_array()
            .or_else(|| plan["degrees"].as_array())
            .unwrap()
            .iter()
            .any(|x| x.as_u64() == Some(ell)));
        let mut adj = vec![vec![]; nodes.len()];
        let mut degree = vec![0u64; nodes.len()];
        let mut edges = db
            .prepare("SELECT source,target FROM edges WHERE ell=? ORDER BY source,target,kernel")
            .unwrap();
        for edge in edges
            .query_map([ell], |r| {
                Ok((r.get::<_, usize>(0)?, r.get::<_, usize>(1)?))
            })
            .unwrap()
        {
            let (i, j) = edge.unwrap();
            if !known[i] {
                continue;
            }
            degree[i] += 1;
            adj[i].push(j);
            adj[j].push(i);
        }
        for i in 0..nodes.len() {
            if known[i] {
                assert!(
                    degree[i] == 1 || degree[i] == ell + 1,
                    "unexpected complete kernel count"
                );
            }
        }
        let floor: Vec<_> = (0..nodes.len())
            .map(|i| known[i] && degree[i] == 1)
            .collect();
        let upper = distances(&adj, (0..nodes.len()).filter(|&i| floor[i]));
        let lower = distances(&adj, (0..nodes.len()).filter(|&i| floor[i] || !known[i]));
        for i in 0..nodes.len() {
            let height = match (lower[i], upper[i]) {
                (Some(l), Some(u)) if l == u && u <= depth => Some(u),
                _ => None,
            };
            if let Some(h) = height {
                exact[i] = exact[i].map(|f| f.checked_mul(ell.pow(depth - h)).unwrap());
            } else {
                exact[i] = None;
            }
            heights[i].push(json!({"ell":ell,"depth":depth,"known_floor":floor[i],"complete_kernel_count":if known[i]{Some(degree[i])}else{None},"lower_height":lower[i],"upper_height":upper[i],"exact_height":height}));
            if let Some(mut truth) = nodes[i].2 {
                let mut exponent = 0;
                while truth % ell == 0 {
                    exponent += 1;
                    truth /= ell;
                }
                assert!(exponent <= depth);
                let true_height = depth - exponent;
                if let Some(l) = lower[i] {
                    assert!(l <= true_height, "lower bound exceeds truth");
                    bounds_checked += 1;
                }
                if let Some(u) = upper[i] {
                    assert!(u >= true_height, "upper bound below truth");
                    bounds_checked += 1;
                }
                if let Some(h) = height {
                    assert_eq!(h, true_height);
                    labels_checked += 1;
                }
            }
        }
    }
    let mut population = BTreeMap::<u64, u64>::new();
    let mut rows = vec![];
    let mut labeled = 0;
    for (i, value) in exact.iter().enumerate() {
        if let Some(f) = value {
            if let Some(truth) = nodes[i].2 {
                assert_eq!(*f, truth);
            }
            labeled += 1;
            *population.entry(*f).or_default() += 1;
        }
        rows.push(json!({"vertex":i,"expanded_visible":known[i],"conductor":value,
            "status":if value.is_some(){"EXACT_LOCAL_DISTANCE_BOUNDS"}else{"UNDETERMINED"},"prime_bounds":heights[i]}));
    }
    json!({"schema":"volcano-local-strata/v1","db_sha256":sha,"cutoff":cutoff,
        "scope":"Only local exact labels; does not establish complete class or global stratum populations",
        "map_replay":replay,"plan":plan,"vertices":nodes.len(),"vertices_exact":labeled,
        "vertices_unknown":nodes.len()-labeled,"strata_exact_observed":population.len(),
        "exact_populations":population.into_iter().map(|(f,n)|json!({"conductor":f,"vertices":n})).collect::<Vec<_>>(),
        "validation":{"all_bounds_and_labels_match_available_complete_truth":true,"bounds_checked":bounds_checked,"prime_labels_checked":labels_checked},"rows":rows})
}

fn main() {
    let a: Vec<_> = env::args().collect();
    if a.get(1).map(String::as_str) == Some("certify") {
        assert_eq!(
            a.len(),
            4,
            "local-strata certify PREFIX_DB FRESH_CERTIFICATE"
        );
        assert!(!Path::new(&a[3]).exists());
        let result = certify(Path::new(&a[2]), u64::MAX);
        fs::write(&a[3], serde_json::to_vec_pretty(&result).unwrap()).unwrap();
        println!(
            "{}",
            json!({"vertices":result["vertices"],"exact":result["vertices_exact"],
            "unknown":result["vertices_unknown"],"strata":result["strata_exact_observed"],
            "source_conductor":result["rows"][0]["conductor"]})
        );
        return;
    }
    if a.get(1).map(String::as_str) == Some("summarize") {
        assert_eq!(
            a.len(),
            5,
            "local-strata summarize FF_RESULTS C_CERTIFICATE FRESH_STATS"
        );
        assert!(!Path::new(&a[4]).exists());
        let raw = fs::read(&a[2]).unwrap();
        let ff: Value = serde_json::from_slice(&raw).unwrap();
        let raw_cert = fs::read(&a[3]).unwrap();
        let cert: Value = serde_json::from_slice(&raw_cert).unwrap();
        assert_eq!(ff["prefix"]["sha256"], cert["db_sha256"]);
        let mut groups = BTreeMap::<(Option<u64>, String), Vec<&Value>>::new();
        let mut identities = std::collections::BTreeSet::new();
        let mut exact_records = 0;
        let mut unknown_records = 0;
        for row in ff["rows"].as_array().unwrap() {
            let id = row["graph_vertex"].as_u64().unwrap();
            assert!(identities.insert((id, row["target_index"].as_u64().unwrap())));
            let c = &cert["rows"][id as usize];
            assert_eq!(c["vertex"], id);
            let conductor = c["conductor"].as_u64();
            if conductor.is_some() {
                exact_records += 1;
            } else {
                unknown_records += 1;
            }
            groups
                .entry((conductor, row["status"].as_str().unwrap().into()))
                .or_default()
                .push(row);
        }
        assert_eq!(
            identities.len() as u64,
            cert["vertices"].as_u64().unwrap() * 2
        );
        let mut summaries = vec![];
        for ((conductor, status), rows) in groups {
            let mut metrics = serde_json::Map::new();
            let mut keys: Vec<(String, String)> =
                vec![("end_to_end_ns".into(), "/end_to_end_ns".into())];
            for key in ff["summary"]["metrics"].as_object().unwrap().keys() {
                keys.push((key.clone(), format!("/relation_report/{key}")));
            }
            for key in ff["summary"]["phase_operations"]
                .as_object()
                .unwrap()
                .keys()
            {
                let (phase, primitive) = key.split_once('/').unwrap();
                keys.push((
                    key.clone(),
                    format!("/phases/{phase}/operations/{primitive}"),
                ));
            }
            for (name, pointer) in keys {
                let values: Vec<f64> = rows
                    .iter()
                    .filter_map(|r| {
                        r.pointer(&pointer).and_then(|v| {
                            v.as_f64()
                                .or_else(|| v.as_str().and_then(|s| s.parse().ok()))
                        })
                    })
                    .collect();
                let value = if values.is_empty() {
                    json!({"n":0,"min":null,"max":null,"mean":null})
                } else {
                    json!({"n":values.len(),"min":values.iter().copied().reduce(f64::min),"max":values.iter().copied().reduce(f64::max),"sum":values.iter().sum::<f64>(),"mean":values.iter().sum::<f64>()/values.len() as f64})
                };
                metrics.insert(name, value);
            }
            summaries.push(json!({"conductor":conductor,"status":status,"records":rows.len(),"metrics":metrics}));
        }
        let result = json!({"schema":"volcano-local-stratum-ff-summary/v1","ff_results_sha256":sha256_hex(&raw),"certificate_sha256":sha256_hex(&raw_cert),"db_sha256":cert["db_sha256"],
            "scope":"Locally exact strata only; unknown conductor kept as a separate null bucket; every status separated; present-value n explicit; not full C class or qualified timing",
            "records":identities.len(),"records_exact_conductor":exact_records,"records_unknown_conductor":unknown_records,"groups":summaries});
        fs::write(&a[4], serde_json::to_vec_pretty(&result).unwrap()).unwrap();
        return;
    }
    assert_eq!(
        a.len(),
        4,
        "local-strata COMPLETE_GRAPH_DIR C_PREFIX FRESH_OUT_DIR"
    );
    let out = Path::new(&a[3]);
    fs::create_dir(out).unwrap();
    let mut calibrations = vec![];
    for label in ["two", "three", "mixed"] {
        for cutoff in [0, 1, 17, 1024, u64::MAX] {
            let result = certify(&Path::new(&a[1]).join(format!("{label}.sqlite")), cutoff);
            calibrations.push(json!({"graph":label,"cutoff":cutoff,"db_sha256":result["db_sha256"],"vertices":result["vertices"],"exact":result["vertices_exact"],"validation":result["validation"]}));
        }
    }
    fs::write(
        out.join("calibration.json"),
        serde_json::to_vec_pretty(&json!({"ok":true,"calibrations":calibrations})).unwrap(),
    )
    .unwrap();
    let result = certify(Path::new(&a[2]), u64::MAX);
    assert_eq!(
        result["db_sha256"],
        "468ce8320b6b7deae37198ca3afc07534fa4fe0e950f621e45ad322f23a2c881"
    );
    fs::write(
        out.join("C-local-strata.json"),
        serde_json::to_vec_pretty(&result).unwrap(),
    )
    .unwrap();
    println!(
        "{}",
        json!({"vertices":result["vertices"],"exact":result["vertices_exact"],"unknown":result["vertices_unknown"],"strata":result["strata_exact_observed"],"source_conductor":result["rows"][0]["conductor"]})
    );
}
