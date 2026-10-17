//! Aggregate recorded values only; absent counters/timeouts stay absent.
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use rusqlite::Connection;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[derive(Default)]
struct Stats {
    n: u64,
    sum: f64,
    min: Option<f64>,
    max: Option<f64>,
}
impl Stats {
    fn add(&mut self, x: f64) {
        assert!(x.is_finite());
        self.n += 1;
        self.sum += x;
        self.min = Some(self.min.map_or(x, |v| v.min(x)));
        self.max = Some(self.max.map_or(x, |v| v.max(x)));
    }
    fn value(&self) -> Value {
        json!({"n":self.n,"sum":self.sum,"min":self.min,"max":self.max,"mean":if self.n>0 {Some(self.sum/self.n as f64)}else{None}})
    }
}
fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        return Err("usage: summarize FRESH_OUT CAMPAIGN...".into());
    }
    let out = Path::new(&a[1]);
    fs::create_dir(out)?;
    let mut obs = File::create(out.join("observations.jsonl"))?;
    let mut groups: BTreeMap<(String, u64, String, String), BTreeMap<String, Stats>> =
        BTreeMap::new();
    let mut statuses: BTreeMap<(String, u64, String, String), u64> = BTreeMap::new();
    let mut outcome_totals: BTreeMap<String, u64> = BTreeMap::new();
    let mut campaigns = Vec::new();
    let mut total = 0;
    let mut verified = 0;
    let mut trials = 0u64;
    let mut relations = 0u64;
    let mut measured_yield = 0u64;
    for path in &a[2..] {
        let root = Path::new(path);
        let manifest: Value = serde_json::from_slice(&fs::read(root.join("campaign.json"))?)?;
        let c = Connection::open_with_flags(
            root.join("ledger.sqlite"),
            rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
        )?;
        let mut jobs =
            c.prepare("SELECT graph,vertex,conductor,state,session FROM jobs ORDER BY ordinal")?;
        let mut coverage: BTreeMap<(String, String), u64> = BTreeMap::new();
        let mut strata: BTreeMap<(String, u64, String), u64> = BTreeMap::new();
        let mut bytes = 0u64;
        let mut sessions = 0;
        for row in jobs.query_map([], |r| {
            Ok((
                r.get::<_, String>(0)?,
                r.get::<_, u64>(1)?,
                r.get::<_, u64>(2)?,
                r.get::<_, String>(3)?,
                r.get::<_, Option<String>>(4)?,
            ))
        })? {
            let (graph, vertex, conductor, state, session) = row?;
            *coverage.entry((graph.clone(), state.clone())).or_default() += 1;
            *strata
                .entry((graph.clone(), conductor, state.clone()))
                .or_default() += 1;
            if state != "audited" {
                continue;
            }
            let session = session.unwrap();
            sessions += 1;
            let records = root.join(&session).join("records.jsonl");
            bytes += fs::metadata(&records)?.len();
            for line in BufReader::new(File::open(&records)?).lines() {
                let v: Value = serde_json::from_str(&line?)?;
                if v["warmup"] == true {
                    continue;
                }
                let arm = v["arm"].as_str().unwrap().to_string();
                let status = v["outcome"]["status"].as_str().unwrap().to_string();
                total += 1;
                *outcome_totals.entry(status.clone()).or_default() += 1;
                if status == "verified" {
                    verified += 1;
                }
                let mut m: BTreeMap<String, f64> = BTreeMap::new();
                for (name, pointer) in [
                    ("cold_gae_lower_bound", "/cost/total_gae"),
                    ("S_lower_bound", "/cost/s"),
                    ("raw_unqualified_solve_seconds", "/time/solve_wall_ns"),
                    ("process_seconds", "/time/process_wall_ns"),
                    ("peak_rss_kib", "/time/max_rss_kib"),
                    ("relation_trials", "/counters/targets_tried"),
                    ("relations_found", "/counters/relations_found"),
                    ("matrix_rank", "/counters/matrix_rank"),
                    ("dependent_relations", "/counters/matrix_dependent"),
                ] {
                    if let Some(x) = v.pointer(pointer).and_then(Value::as_f64) {
                        m.insert(
                            name.into(),
                            if name.ends_with("seconds") {
                                x / 1e9
                            } else {
                                x
                            },
                        );
                    }
                }
                if let (Some(t), Some(r)) = (
                    v["counters"]["targets_tried"].as_u64(),
                    v["counters"]["relations_found"].as_u64(),
                ) {
                    trials += t;
                    relations += r;
                    measured_yield += 1;
                    if t > 0 {
                        m.insert("relation_yield".into(), r as f64 / t as f64);
                    }
                }
                let mut native = BTreeMap::new();
                if let Some(phases) = v["phases"].as_array() {
                    for p in phases {
                        let name = p["name"].as_str().unwrap();
                        for field in ["adds", "doubles", "scalar_mults", "gae", "wall_ns"] {
                            if let Some(x) = p[field].as_f64() {
                                m.insert(format!("phase.{name}.{field}"), x);
                            }
                        }
                        if let Some(obj) = p["native"].as_object() {
                            for (k, x) in obj {
                                if let Some(n) = x.as_f64() {
                                    native.insert(format!("{name}.{k}"), n);
                                    m.insert(format!("native.{name}.{k}"), n);
                                }
                            }
                        }
                    }
                }
                if let Some(counters) = v["counters"].as_object() {
                    for (k, x) in counters {
                        if let Some(n) = x.as_f64() {
                            m.insert(format!("counter.{k}"), n);
                        }
                    }
                }
                // Keep graph totals as well as each actual conductor stratum.
                for stratum in [0, conductor] {
                    let key = (graph.clone(), stratum, arm.clone(), status.clone());
                    *statuses.entry(key.clone()).or_default() += 1;
                    let group = groups.entry(key).or_default();
                    for (k, x) in &m {
                        group.entry(k.clone()).or_default().add(*x);
                    }
                }
                writeln!(
                    obs,
                    "{}",
                    json!({"campaign":root,"session":session,"graph":graph,"vertex":vertex,"conductor":conductor,"record_id":v["record_id"],"run_id":v["run_id"],"workload_id":v["workload"]["workload_id"],"arm":arm,"status":status,"isolation":v["isolation"]["level"],"metrics":m,"native":native,"unpriced":v["cost"]["unpriced"]})
                )?;
            }
        }
        let coverage: Vec<Value> = coverage
            .into_iter()
            .map(|((g, s), n)| json!({"graph":g,"state":s,"vertices":n}))
            .collect();
        let strata_audited = strata.keys().filter(|(_, _, s)| s == "audited").count();
        let stratum_coverage: Vec<Value> = strata
            .into_iter()
            .map(|((g, f, s), n)| json!({"graph":g,"conductor":f,"state":s,"vertices":n}))
            .collect();
        campaigns.push(json!({"path":root,"manifest_sha256":sha256_hex(&fs::read(root.join("campaign.json"))?),"mode":manifest["mode"],"planned_vertices":manifest["vertices_planned"],"coverage":coverage,"stratum_coverage":stratum_coverage,"strata_audited":strata_audited,"audited_sessions":sessions,"audited_records_bytes":bytes}));
    }
    obs.sync_all()?;
    let mut csv = File::create(out.join("statistics.csv"))?;
    writeln!(
        csv,
        "graph,conductor_0_means_all,arm,status,metric,n,min,max,mean"
    )?;
    let mut summaries = Vec::new();
    for ((graph, conductor, arm, status), metrics) in groups {
        let mut output = BTreeMap::new();
        for (k, s) in metrics {
            writeln!(
                csv,
                "{graph},{conductor},{arm},{status},{k},{},{},{},{}",
                s.n,
                s.min.unwrap(),
                s.max.unwrap(),
                s.sum / s.n as f64
            )?;
            output.insert(k, s.value());
        }
        let n = statuses[&(graph.clone(), conductor, arm.clone(), status.clone())];
        summaries.push(json!({"graph":graph,"conductor_0_means_all":conductor,"arm":arm,"status":status,"records":n,"metrics":output}));
    }
    let summary = json!({"schema":"volcano-ic-summary/v1","wall_claim":"UNQUALIFIED; L0/debug diagnostic; raw observations only","aggregation":"Each status separate. Means over present values only; n explicit. Conductor 0 denotes graph aggregate, not a stratum. Timeout counts are unknown, not zero.","campaigns":campaigns,"measured_records":total,"verified_records":verified,"outcome_totals":outcome_totals,"relation_observations":measured_yield,"relation_trials":trials,"relations_found":relations,"weighted_relation_yield":if trials>0{Some(relations as f64/trials as f64)}else{None},"groups":summaries});
    serde_json::to_writer_pretty(File::create(out.join("summary.json"))?, &summary)?;
    println!(
        "{}",
        json!({"measured_records":total,"verified_records":verified,"relation_trials":trials,"relations_found":relations,"out":out})
    );
    Ok(())
}
