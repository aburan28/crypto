//! Native input preparation and analysis. Python is used only for rendering.
mod campaign_analysis;
mod certificates;
mod index;
mod replay_campaign;
mod replication;
mod rho_analysis;
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
        "campaign-analysis" => campaign_analysis::run(),
        "rho-analysis" => rho_analysis::run(),
        "replay-campaign" => replay_campaign::run(),
        "cold-batches" => cold_batches(),
        "index" => index::run(),
        "dashboard" => index::dashboard(),
        "certificates" => certificates::run(),
        "replication" => replication::run(&args[2], &args[3]),
        "upper-ids" => replication::upper_ids(),
        "constructed-floor-input" => {
            let raw = read(&args[2]);
            let searched =
                fs::read_to_string(Path::new(ROOT).join("prepared/search-certificates.jsonl"))
                    .unwrap();
            let crater: Value = searched
                .lines()
                .map(|s| serde_json::from_str::<Value>(s).unwrap())
                .find(|c| c["id"] == raw["id"])
                .expect("unknown construction id");
            assert_eq!(raw["p"], crater["p"]);
            assert_eq!(
                raw["ell"].as_u64().unwrap().to_string(),
                crater["ell"].as_str().unwrap()
            );
            assert_eq!(raw["reference_supplied"], false);
            let decimal = |v: &Value| v.as_str().unwrap().to_string();
            fs::write(
                &args[3],
                format!(
                    "[{},{},{},{},{},{},{},{},{}]\n",
                    raw["id"],
                    decimal(&raw["p"]),
                    decimal(&raw["a"]),
                    decimal(&raw["b"]),
                    decimal(&crater["N"]),
                    decimal(&crater["ell"]),
                    crater["D"],
                    decimal(&crater["a"]),
                    decimal(&crater["b"])
                ),
            )
            .unwrap();
            if args.len() > 4 {
                assert_eq!(
                    args.len(),
                    6,
                    "optional torsion output requires kernel input"
                );
                let kernel = fs::read_to_string(&args[5]).unwrap();
                let f: Vec<_> = kernel.split_whitespace().collect();
                assert_eq!(f.len(), 12);
                assert_eq!(f[0], raw["id"].to_string());
                assert_eq!(f[1], decimal(&crater["ell"]));
                assert_eq!(f[2], decimal(&raw["p"]));
                assert_eq!(f[8], decimal(&crater["a"]));
                assert_eq!(f[9], decimal(&crater["b"]));
                let q = (1 - crater["D"].as_i64().unwrap()) / 4;
                assert!(q == 5 || q == 11);
                fs::write(
                    &args[4],
                    format!(
                        "[{},{},{},{},{},{},{},{},{},{},{},{},{},{}]\n",
                        f[0],
                        f[1],
                        f[2],
                        decimal(&crater["N"]),
                        q,
                        f[3],
                        f[4],
                        f[5],
                        f[6],
                        f[7],
                        f[8],
                        f[9],
                        decimal(&raw["a"]),
                        decimal(&raw["b"])
                    ),
                )
                .unwrap();
            }
        }
        "controls" => replication::controls(&args[2], &args[3]),
        "recurrence-inputs" => replication::recurrence_inputs(),
        "torsion-inputs" => replication::torsion_inputs(),
        "certified-floor-inputs" => {
            let input =
                fs::read_to_string(Path::new(ROOT).join("prepared/torsion-structure.gpvec"))
                    .unwrap();
            let mut out = String::new();
            for line in input.lines() {
                let v: Value = serde_json::from_str(line).unwrap();
                out += &format!(
                    "[{},{},{},{},{},{},{},{},{}]\n",
                    v[0],
                    v[2],
                    v[12],
                    v[13],
                    v[3],
                    v[1],
                    1 - 4 * v[4].as_i64().unwrap(),
                    v[10],
                    v[11]
                );
            }
            fs::write(
                Path::new(ROOT).join("prepared/all-certified-floors.gpvec"),
                out,
            )
            .unwrap();
        }
        "recurrence-controls" => replication::recurrence_controls(&args[2]),
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
    let mut recurrence_matches = 0;
    for ((pair, role, rep), h) in &fingerprints {
        if *role == 0 {
            assert_eq!(Some(h), fingerprints.get(&(*pair, 2, *rep)));
            identical += 1;
        }
        if *role == 4 || *role == 5 {
            let direct = if *role == 4 { 0 } else { 1 };
            assert_eq!(Some(h), fingerprints.get(&(*pair, direct, *rep)));
            recurrence_matches += 1;
        }
    }
    let mean = |a: &[f64]| a.iter().sum::<f64>() / a.len() as f64;
    let data:Vec<_>=groups.iter().map(|((pair,role),v)|json!({"pair":pair,"role":role,"repetitions":v.len(),"mean_cpu_seconds":mean(v),"min_cpu_seconds":v.iter().cloned().fold(f64::INFINITY,f64::min),"max_cpu_seconds":v.iter().cloned().fold(0.,f64::max)})).collect();
    let mut aggregate = Vec::new();
    for role in 0..6 {
        let vals: Vec<f64> = rows
            .iter()
            .filter(|v| n(v, "role") == role)
            .map(|v| v["orbit_cpu_seconds"].as_f64().unwrap())
            .collect();
        if vals.is_empty() {
            continue;
        }
        aggregate.push(json!({"role":role,"mean_cpu_seconds":mean(&vals),"runs":vals.len()}));
        let last = aggregate.last_mut().unwrap();
        last["points_per_construction"] = json!(8192);
        last["orbit_steps_per_construction"] = json!(8188);
        last["amortized_cpu_seconds_per_orbit_step"] = json!(mean(&vals) / 8188.0);
        let first: Vec<_> = rows
            .iter()
            .filter(|v| n(v, "role") == role && n(v, "repeat") == 0)
            .map(|v| v["orbit_cpu_seconds"].as_f64().unwrap())
            .collect();
        let later: Vec<_> = rows
            .iter()
            .filter(|v| n(v, "role") == role && n(v, "repeat") > 0)
            .map(|v| v["orbit_cpu_seconds"].as_f64().unwrap())
            .collect();
        let last = aggregate.last_mut().unwrap();
        last["first_pass_mean_cpu_seconds"] = json!(mean(&first));
        last["later_pass_mean_cpu_seconds"] = json!(mean(&later));
        last["pass_scope"]=json!("First observed model/arm pass versus later shuffled passes in one fresh process; not a fresh process or cache flush per arm.");
    }
    write(
        out,
        &json!({"schema":"cm-orbits.summary/v1","source_sha256":sha256_hex(text.as_bytes()),"runs":rows.len(),"identical_same_crater_pairs":identical,"recurrence_same_endpoint_matches":recurrence_matches,"fingerprint":"FNV64 supplementary diagnostic; backend separately compares complete vectors in memory","isolation":"L0; descriptive CPU observations only","per_pair":data,"aggregate":aggregate,"not_measured":"ECDLP hardness or a complete online IC/rho comparison"}),
    );
    println!(
        "{} construction runs, {} matching crater comparisons",
        rows.len(),
        identical
    );
    let _ = canonical(&json!({"v":1}));
}

fn cold_batches() {
    let root = Path::new(ROOT);
    let mut rows = Vec::new();
    let mut orbit_cpu = BTreeMap::<u64, f64>::new();
    for line in fs::read_to_string(root.join("results/cold-orbits.jsonl"))
        .unwrap()
        .lines()
    {
        let v: Value = serde_json::from_str(line).unwrap();
        *orbit_cpu.entry(n(&v, "repeat")).or_default() += v["orbit_cpu_seconds"].as_f64().unwrap();
    }
    for pass in 0..20 {
        let path = root.join(format!("results/cold-passes/{pass}-time.txt"));
        let text = fs::read_to_string(&path).unwrap();
        let field = |key: &str| {
            text.lines()
                .find_map(|line| line.trim().strip_prefix(key))
                .unwrap()
                .trim()
                .to_string()
        };
        assert_eq!(field("Exit status:"), "0");
        let cpu = field("User time (seconds):").parse::<f64>().unwrap()
            + field("System time (seconds):").parse::<f64>().unwrap();
        let elapsed_line = text
            .lines()
            .find(|line| line.trim().starts_with("Elapsed (wall clock)"))
            .unwrap();
        let parts: Vec<f64> = elapsed_line
            .rsplit_once(": ")
            .unwrap()
            .1
            .split(':')
            .map(|v| v.parse().unwrap())
            .collect();
        let wall = parts.into_iter().fold(0.0, |sum, v| 60.0 * sum + v);
        rows.push(json!({"pass":pass,"process_cpu_seconds":cpu,"wall_seconds":wall,"native_orbit_cpu_seconds":orbit_cpu[&pass],"outside_orbit_cpu_seconds":cpu-orbit_cpu[&pass],"time_source":path.to_str().unwrap(),"time_source_sha256":sha256_hex(text.as_bytes())}));
    }
    let mean =
        |key: &str| rows.iter().map(|r| r[key].as_f64().unwrap()).sum::<f64>() / rows.len() as f64;
    write(
        root.join("results/cold-batches.json"),
        &json!({"schema":"conductor-gap.cold-batches/v1","processes":20,"jobs_per_process":600,"mean_process_cpu_seconds":mean("process_cpu_seconds"),"mean_native_orbit_cpu_seconds":mean("native_orbit_cpu_seconds"),"mean_outside_orbit_cpu_seconds":mean("outside_orbit_cpu_seconds"),"mean_wall_seconds":mean("wall_seconds"),"scope":"Fresh process per 600-job batch, not per arm. Precomputed CM-map parameter files; excludes prior map derivation. Outside-orbit CPU includes parsing, hashing, equality checks, serialization and thread/process housekeeping. L0 co-scheduled; no causal cold/warm timing or ECDLP speed claim.","batches":rows}),
    );
}
