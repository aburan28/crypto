//! Native input preparation and analysis. Python is used only for rendering.
mod certificates;
mod index;
use crypto_lib::cryptanalysis::{
    ecbench::{
        canonical::{canonical, sha256_hex},
        workload::CurveSpec,
    },
    ic_boundary::{GroupOps, PrimeCurve, PrimeInstance},
    ic_framework::cm_orbit::derive_map,
};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, HashMap},
    fs,
    path::Path,
};
const ROOT: &str = "research/conductor_gap_cm_orbits_20261015";
fn read(p: impl AsRef<Path>) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn write(p: impl AsRef<Path>, v: &Value) {
    fs::write(p, serde_json::to_string_pretty(v).unwrap() + "\n").unwrap();
}
fn n(v: &Value, k: &str) -> u64 {
    v[k].as_u64().unwrap()
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    match args[1].as_str() {
        "prepare" => prepare(&args[2]),
        "index" => index::run(),
        "browser" => index::browser(),
        "leaderboard" => index::leaderboard(),
        "certificates" => certificates::run(),
        "summarize" => summarize(&args[2], &args[3]),
        _ => panic!("prepare BACKEND | summarize JSONL OUT"),
    }
}
fn prepare(backend: &str) {
    let root = Path::new(ROOT);
    fs::create_dir_all(root.join("prepared")).unwrap();
    let pairs = read(root.join("inputs/pairs.json"));
    let hash = sha256_hex(&fs::read(backend).unwrap());
    let mut maps = Vec::new();
    let mut params = String::new();
    let mut kernels = String::new();
    let mut catalog = Vec::new();
    let mut orbit = String::from("job,pair,role,p,r,a,b,Gx,Gy,Qx,Qy,seed,repeat\n");
    let csv = fs::read_to_string(root.join("inputs/workloads.csv")).unwrap();
    let targets: HashMap<_, _> = csv
        .lines()
        .skip(1)
        .map(|l| {
            let v: Vec<u64> = l.split(',').map(|x| x.parse().unwrap()).collect();
            ((v[1], v[2], v[12]), v)
        })
        .collect();
    let mut shards: Vec<Vec<Value>> = vec![Vec::new(); 8];
    for v in pairs.as_array().unwrap() {
        let id = n(v, "index");
        let p = n(v, "p");
        let r = n(v, "N");
        let g0 = (v["G0"][0].as_u64().unwrap(), v["G0"][1].as_u64().unwrap());
        let inst = PrimeInstance {
            name: format!("CM gap {id}"),
            curve: PrimeCurve {
                p,
                a: n(v, "a0"),
                b: n(v, "b0"),
            },
            group_order: r,
            r,
            cofactor: 1,
            generator: g0,
        };
        let map = derive_map(&inst, 19, &mut GroupOps::default())
            .unwrap_or_else(|e| panic!("map {id}: {e}"));
        params += &format!(
            "{p} {r} {id} {} {} {} {} {} {} {} {} {}\n",
            map.lambda,
            map.inv_c,
            map.inv_cu,
            map.kernel_x[0],
            map.kernel_x[1],
            map.kernel_t[0],
            map.kernel_t[1],
            map.kernel_rhs[0],
            map.kernel_rhs[1]
        );
        maps.push(json!({"pair":id,"lambda":map.lambda,"degree":map.degree,"kernel_x":map.kernel_x,"inverse_c":map.inv_c,"inverse_cu":map.inv_cu,"generator_action_checked":true}));
        kernels += &format!(
            "{id} {} {p} {} {} {} {} {} {} {} {} {}\n",
            n(v, "ell"),
            n(v, "twist_a"),
            n(v, "twist_b"),
            v["kernel_G"][0],
            v["kernel_G"][1],
            n(v, "twist_d"),
            n(v, "a0"),
            n(v, "b0"),
            n(v, "a1"),
            n(v, "b1")
        );
        for role in 0..=1 {
            let w = &targets[&(id, role, 0)];
            let spec = json!({"kind":"prime_public_fixture","p":p,"a":w[5],"b":w[6],"group_order":r,"r":r,"gx":w[7],"gy":w[8],"qx":w[9],"qy":w[10]});
            let cs: CurveSpec = serde_json::from_value(spec.clone()).unwrap();
            let i = cs.build().unwrap();
            let facts = i.facts(&cs);
            catalog.push(json!({"pair":id,"role":role,"facts":facts,"spec":spec,"conductor":if role==0{1}else{n(v,"ell")}}));
            shards[(id % 4) as usize + 4 * role as usize].push(spec);
            for arm in if role == 0 { vec![0, 2, 3] } else { vec![1] } {
                orbit += &format!(
                    "{},{id},{arm},{p},{r},{},{},{},{},{},{},{},0\n",
                    id * 3 + arm,
                    w[5],
                    w[6],
                    w[7],
                    w[8],
                    w[9],
                    w[10],
                    w[11]
                );
            }
        }
    }
    fs::write(root.join("prepared/cm-params.txt"), params).unwrap();
    fs::write(root.join("prepared/kernels.txt"), kernels).unwrap();
    fs::write(root.join("prepared/orbits.csv"), orbit).unwrap();
    write(root.join("prepared/cm-maps.json"), &json!(maps));
    write(root.join("prepared/catalog.json"), &json!(catalog));
    for (i, curves) in shards.into_iter().enumerate() {
        let method = |mode: &str| json!({"id":"ic.cm_folded_native_v1","params":{"backend":backend,"backend_sha256":hash,"mode":mode}});
        let mut arms = vec![
            json!({"name":"scalar","role":"baseline","method":method("scalar")}),
            json!({"name":"rho","role":"reference","method":{"id":"rho.negation"}}),
        ];
        if i < 4 {
            arms.push(json!({"name":"velu","role":"candidate","method":method("velu")}));
        }
        write(
            root.join(format!("prepared/shard-{i}.json")),
            &json!({"schema":"ecbench.spec/v1","label":format!("CM orbit public regression shard {i}"),"description":"One target per workload; L0 descriptive native costs, no eligible speedup claim.","workloads":{"curves":curves,"targets_per_curve":1,"target_seed":20261015,"target_kind":"public"},"arms":arms,"measurement":{"rounds":1,"warmup":0,"order":"shuffle","seed":20261015,"isolation_required":"L0","timeout_seconds":1900}}),
        );
    }
    write(
        root.join("prepared/input-manifest.json"),
        &json!({"schema":"cm-orbits.inputs/v1","backend_sha256":hash,"pairs":100,"endpoints":200,"CM_maps":maps.len(),"orbit_arms":400,"canonical_CM_trace":1,"primary_targets":"original 200 public regression fixtures; validation scalars excluded","prep_source_sha256":sha256_hex(include_bytes!("study.rs"))}),
    );
    // A separate source catalog is append-only provenance; main registry is updated
    // by the native registration command after generated artifacts are validated.
    println!("prepared {} maps, {} endpoints", maps.len(), catalog.len());
}
fn summarize(input: &str, out: &str) {
    let text = fs::read_to_string(input).unwrap();
    let rows: Vec<Value> = text
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let mut groups: BTreeMap<(u64, u64), Vec<f64>> = BTreeMap::new();
    let mut fingerprints: BTreeMap<(u64, u64, u64), String> = BTreeMap::new();
    for v in &rows {
        groups
            .entry((n(v, "pair"), n(v, "role")))
            .or_default()
            .push(v["orbit_cpu_seconds"].as_f64().unwrap());
        fingerprints.insert(
            (n(v, "pair"), n(v, "role"), n(v, "repeat")),
            v["point_record_fnv64"].as_str().unwrap().into(),
        );
    }
    let mut identical = 0;
    for ((pair, role, rep), h) in &fingerprints {
        if *role == 0 {
            assert_eq!(Some(h), fingerprints.get(&(*pair, 2, *rep)));
            identical += 1;
        }
    }
    let mean = |a: &[f64]| a.iter().sum::<f64>() / a.len() as f64;
    let data:Vec<_>=groups.iter().map(|((pair,role),v)|json!({"pair":pair,"role":role,"repetitions":v.len(),"mean_cpu_seconds":mean(v),"min_cpu_seconds":v.iter().cloned().fold(f64::INFINITY,f64::min),"max_cpu_seconds":v.iter().cloned().fold(0.,f64::max)})).collect();
    let mut aggregate = Vec::new();
    for role in 0..4 {
        let vals: Vec<f64> = rows
            .iter()
            .filter(|v| n(v, "role") == role)
            .map(|v| v["orbit_cpu_seconds"].as_f64().unwrap())
            .collect();
        aggregate.push(json!({"role":role,"mean_cpu_seconds":mean(&vals),"runs":vals.len()}));
    }
    write(
        out,
        &json!({"schema":"cm-orbits.summary/v1","source_sha256":sha256_hex(text.as_bytes()),"runs":rows.len(),"identical_same_crater_pairs":identical,"fingerprint":"FNV64 over all serialized point/column/coefficient records; noncryptographic equality diagnostic","isolation":"L0; descriptive CPU observations only","per_pair":data,"aggregate":aggregate,"not_measured":"ECDLP hardness or a complete online IC/rho comparison"}),
    );
    println!(
        "{} construction runs, {} matching crater comparisons",
        rows.len(),
        identical
    );
    let _ = canonical(&json!({"v":1}));
}
