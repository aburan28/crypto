//! Native descriptive aggregation and accounting checks; never pool replay rows.
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    env, fs,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};
#[derive(Default)]
struct Stat {
    n: u64,
    min: f64,
    max: f64,
    sum: f64,
}
impl Stat {
    fn add(&mut self, x: f64) {
        if self.n == 0 {
            self.min = x;
            self.max = x;
        } else {
            self.min = self.min.min(x);
            self.max = self.max.max(x);
        }
        self.n += 1;
        self.sum += x;
    }
    fn value(&self) -> Value {
        json!({"n":self.n,"min":self.min,"max":self.max,"sum":self.sum,"mean":self.sum/self.n as f64})
    }
}
fn main() {
    let a: Vec<_> = env::args().collect();
    assert_eq!(a.len(), 3, "summarize RUN_DIR FRESH_OUT");
    let run = Path::new(&a[1]);
    let out = Path::new(&a[2]);
    assert!(!out.exists());
    fs::create_dir(out).unwrap();
    let meta: Value = serde_json::from_slice(&fs::read(run.join("meta.json")).unwrap()).unwrap();
    let mut groups = BTreeMap::<(String, Option<u64>, String), BTreeMap<String, Stat>>::new();
    let mut global = BTreeMap::<(String, String), BTreeMap<String, Stat>>::new();
    let mut pairs = BTreeMap::<(u64, u64), Vec<Value>>::new();
    let mut outcomes = BTreeMap::<String, u64>::new();
    let mut totals = BTreeMap::<String, u64>::new();
    let mut count = 0u64;
    let (mut smallest, mut largest) = (u64::MAX, 0u64);
    for line in BufReader::new(fs::File::open(run.join("records.jsonl")).unwrap()).lines() {
        let row: Value = serde_json::from_str(&line.unwrap()).unwrap();
        count += 1;
        let method = row["method"].as_str().unwrap().to_owned();
        assert!([
            "ec.ph_bsgs.full_point_v1",
            "ff.k1_ph_relations.full_point_v1"
        ]
        .contains(&method.as_str()));
        let status = row["status"].as_str().unwrap().to_owned();
        let conductor = row["conductor"].as_u64();
        let vertex = row["vertex"].as_u64().unwrap();
        let target = row["target_index"].as_u64().unwrap();
        assert!(target < 2);
        *outcomes.entry(format!("{method}/{status}")).or_default() += 1;
        let mut metrics = BTreeMap::<String, f64>::new();
        if let Some(clock) = row["end_to_end_ns"].as_str() {
            metrics.insert("raw_end_to_end_ns".into(), clock.parse().unwrap());
        }
        if let Some(phases) = row["trace"]["phases"].as_object().filter(|m| !m.is_empty()) {
            for p in phases.values() {
                for (name, pointer) in [
                    ("ec_adds", "/ec_group/adds"),
                    ("ec_doubles", "/ec_group/doubles"),
                    ("ec_scalar_mult_calls", "/ec_group/scalar_mults"),
                    ("ec_field_muls", "/ec_field/muls"),
                    ("ec_field_squares", "/ec_field/sqrs"),
                    ("ec_field_inversions_unexpanded", "/ec_field/invs"),
                    ("explicit_fp_muls", "/explicit_fp/multiplications"),
                    ("explicit_fp_squares", "/explicit_fp/squarings"),
                    (
                        "explicit_fp_inversions_already_expanded",
                        "/explicit_fp/inversions",
                    ),
                    ("sqrt_calls_unexpanded", "/sqrt_calls_unexpanded"),
                    ("miller_steps", "/miller_steps"),
                ] {
                    let n = p
                        .pointer(pointer)
                        .and_then(Value::as_u64)
                        .expect("phase counter missing");
                    *metrics.entry(name.into()).or_default() += n as f64;
                }
            }
            metrics.insert(
                "ec_adds_and_doubles".into(),
                metrics["ec_adds"] + metrics["ec_doubles"],
            );
            metrics.insert(
                "explicit_fp_muls_and_squares".into(),
                metrics["explicit_fp_muls"] + metrics["explicit_fp_squares"],
            );
        }
        if let Some(counters) = row["trace"]["counters"].as_object() {
            for (k, v) in counters {
                let n = v.as_u64().unwrap();
                metrics.insert(k.clone(), n as f64);
                *totals.entry(format!("{method}/{k}")).or_default() += n;
            }
            for tag in ["base_relations", "target_relations"] {
                if let Some(q) = counters
                    .get(&format!("{tag}_queries"))
                    .and_then(Value::as_u64)
                    .filter(|q| *q > 0)
                {
                    let r = counters[&format!("{tag}_relations")].as_u64().unwrap_or(0);
                    assert!(r <= q);
                    metrics.insert(format!("{tag}_equations_per_query"), r as f64 / q as f64);
                }
            }
        }
        if status == "verified" {
            assert_eq!(row["known_answer_checked"], true);
            let mut m: u64 = row["point_order"].as_str().unwrap().parse().unwrap();
            assert!(m >= 42 && m % 42 == 0);
            smallest = smallest.min(m);
            largest = largest.max(m);
            let mut digits = 0;
            for q in [2, 3, 7] {
                while m % q == 0 {
                    m /= q;
                    digits += 1;
                }
            }
            assert_eq!(m, 1);
            let details = row["trace"]["details"].as_array().unwrap();
            assert_eq!(
                details.iter().filter(|d| d["kind"] == "digit").count(),
                digits
            );
            if method.starts_with("ff.") {
                let anchors = details
                    .iter()
                    .filter(|d| d["kind"] == "anchor_relation")
                    .count() as u64;
                let targets = details
                    .iter()
                    .filter(|d| d["kind"] == "target_relation")
                    .count() as u64;
                assert_eq!(targets, digits as u64);
                assert_eq!(
                    row["trace"]["counters"]["target_relations_relations"],
                    targets
                );
                assert_eq!(
                    row["trace"]["counters"]["base_relations_relations"],
                    anchors
                );
                let bases: Vec<_> = details.iter().filter(|d| d["kind"] == "base").collect();
                assert_eq!(bases.len(), 3);
                let ranks: u64 = bases.iter().map(|d| d["rank"].as_u64().unwrap()).sum();
                assert_eq!(ranks, 4);
                assert_eq!(
                    anchors,
                    ranks
                        + row["trace"]["counters"]["matrix_dependent_rows"]
                            .as_u64()
                            .unwrap()
                );
                assert_eq!(
                    details
                        .iter()
                        .filter(
                            |d| d["kind"] == "transfer" && d["exact_image_order_checked"] == true
                        )
                        .count(),
                    3
                );
            }
        }
        pairs.entry((vertex,target)).or_default().push(json!({"method":method,"status":status,"workload":row["workload_sha256"],"order":row["point_order"],"recovered":row["recovered"]}));
        for (name, value) in metrics {
            groups
                .entry((method.clone(), conductor, status.clone()))
                .or_default()
                .entry(name.clone())
                .or_default()
                .add(value);
            global
                .entry((method.clone(), status.clone()))
                .or_default()
                .entry(name)
                .or_default()
                .add(value);
        }
    }
    assert_eq!(count, meta["method_records_requested"].as_u64().unwrap());
    let vertices = meta["vertices_requested"].as_u64().unwrap();
    assert_eq!(pairs.len() as u64, vertices * 2);
    let mut paired_verified = 0u64;
    for ((vertex, _), p) in &pairs {
        if meta["mode"] == "prefix" {
            assert!(*vertex < vertices);
        } else {
            assert_eq!(*vertex, 0);
        }
        assert_eq!(p.len(), 2);
        assert_ne!(p[0]["method"], p[1]["method"]);
        assert_eq!(p[0]["workload"], p[1]["workload"]);
        assert_eq!(p[0]["order"], p[1]["order"]);
        if p[0]["status"] == "verified" && p[1]["status"] == "verified" {
            assert_eq!(p[0]["recovered"], p[1]["recovered"]);
            paired_verified += 1;
        }
    }
    let mut csv = BufWriter::new(fs::File::create(out.join("aggregate-statistics.csv")).unwrap());
    writeln!(csv, "method,stratum,status,metric,n,min,max,mean").unwrap();
    let mut rows = vec![];
    for ((method, conductor, status), metrics) in groups {
        let mut data = serde_json::Map::new();
        for (name, s) in metrics {
            writeln!(
                csv,
                "{method},{},{status},{name},{},{},{},{}",
                conductor.map(|f| f.to_string()).unwrap_or("unknown".into()),
                s.n,
                s.min,
                s.max,
                s.sum / s.n as f64
            )
            .unwrap();
            data.insert(name, s.value());
        }
        rows.push(json!({"method":method,"conductor":conductor,"status":status,"metrics":data}));
    }
    let overall:Vec<_>=global.into_iter().map(|((m,s),v)|json!({"method":m,"status":s,"metrics":v.into_iter().map(|(k,x)|(k,x.value())).collect::<BTreeMap<_,_>>()})).collect();
    let doc = json!({"schema":"C-full-point-native-summary/v1","records":count,"paired_verified":paired_verified,"outcomes":outcomes,"point_order_min":smallest.to_string(),"point_order_max":largest.to_string(),"totals":totals,"overall":overall,"strata":rows,"records_sha256":sha256_hex(&fs::read(run.join("records.jsonl")).unwrap()),"aggregation_source_sha256":sha256_hex(include_bytes!("summarize.rs")),"scope":"Native sums/min/max/means. Primitive categories remain separate, not a priced total. EC inversions and square-root internals unexpanded. Raw clock observations are unqualified. Replay rows are not pooled. Accounting conservation and metadata pairing checks are not independent arithmetic replay."});
    fs::write(
        out.join("summary.json"),
        serde_json::to_vec_pretty(&doc).unwrap(),
    )
    .unwrap();
    println!(
        "{}",
        json!({"records":count,"paired_verified":paired_verified,"outcomes":outcomes,"point_order_min":smallest.to_string(),"point_order_max":largest.to_string()})
    );
}
