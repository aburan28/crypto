//! Freeze replication inputs without exposing target logarithms.
use super::*;
use std::collections::HashSet;
fn lines(path: &str) -> Vec<Value> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|x| !x.is_empty())
        .map(|x| serde_json::from_str(x).unwrap())
        .collect()
}
fn textn(v: &Value, k: &str) -> u64 {
    v[k].as_str().unwrap().parse().unwrap()
}
fn mul(a: u64, b: u64, p: u64) -> u64 {
    (a as u128 * b as u128 % p as u128) as u64
}
fn explicit(v: &Value, role: u64) -> Value {
    let g = &v[format!("G{role}")];
    json!({"kind":"prime_explicit","p":v["p"],"a":v[format!("a{role}")],"b":v[format!("b{role}")],"group_order":v["N"],"r":v["N"],"gx":g[0],"gy":g[1]})
}
fn fromcert(v: &Value) -> Value {
    json!({"kind":"prime_explicit","p":textn(v,"p"),"a":textn(v,"a"),"b":textn(v,"b"),"group_order":textn(v,"N"),"r":textn(v,"N"),"gx":textn(v,"Gx"),"gy":textn(v,"Gy")})
}
fn freeze(
    name: &str,
    curves: Vec<Value>,
    root: bool,
    d43: bool,
    targets: u64,
    rounds: u64,
    public_fixture: bool,
    backend: &str,
    hash: &str,
) {
    let method = |mode: &str| json!({"id":if d43{"ic.cm_folded_native_d43_v1"}else{"ic.cm_folded_native_v1"},"params":{"backend":backend,"backend_sha256":hash,"mode":mode}});
    let mut arms = vec![
        json!({"name":"scalar","role":"baseline","method":method("scalar")}),
        json!({"name":"rho","role":"reference","method":{"id":"rho.negation"}}),
    ];
    if root {
        arms.push(json!({"name":"velu","role":"candidate","method":method("velu")}));
    }
    let path = Path::new(ROOT).join("prepared/replication");
    fs::create_dir_all(&path).unwrap();
    write(
        path.join(format!("{name}.json")),
        &json!({"schema":"ecbench.spec/v1","label":format!("Conductor-gap replication {name}"),"description":if public_fixture{"Matched explicit public regression point; validation scalar excluded"}else{"Fresh public hash-to-subgroup targets; no logarithm supplied or computed during preparation"},"workloads":{"curves":curves,"targets_per_curve":targets,"target_seed":20261016,"target_kind":"public"},"arms":arms,"measurement":{"rounds":rounds,"warmup":0,"order":"alternate","seed":20261017,"isolation_required":"L0","timeout_seconds":1900}}),
    );
}
fn shards(
    label: &str,
    curves: Vec<Value>,
    root: bool,
    d43: bool,
    targets: u64,
    backend: &str,
    hash: &str,
) {
    for shard in 0..4 {
        let c: Vec<_> = curves
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 4 == shard)
            .map(|(_, v)| v.clone())
            .collect();
        if !c.is_empty() {
            let fixtures = c.iter().all(|v| v["kind"] == "prime_public_fixture");
            freeze(
                &format!("{label}-{shard}"),
                c,
                root,
                d43,
                targets,
                1,
                fixtures,
                backend,
                hash,
            );
        }
    }
}
pub fn run(q5: &str, q11: &str) {
    let root = Path::new(ROOT);
    let pairs = read(root.join("inputs/pairs.json"));
    let data = pairs.as_array().unwrap();
    let csv = fs::read_to_string(root.join("inputs/workloads.csv")).unwrap();
    let rows: Vec<Vec<u64>> = csv
        .lines()
        .skip(1)
        .map(|l| l.split(',').map(|x| x.parse().unwrap()).collect())
        .collect();
    let mut selected = Vec::new();
    for w in &rows {
        if w[12] == 0 && !selected.contains(&w[1]) {
            selected.push(w[1]);
            if selected.len() == 10 {
                break;
            }
        }
    }
    let q5hash = sha256_hex(&fs::read(q5).unwrap());
    let q11hash = sha256_hex(&fs::read(q11).unwrap());
    let mut fresh = [Vec::new(), Vec::new()];
    let mut iso = [Vec::new(), Vec::new()];
    let mut model_records = Vec::new();
    let mut public_inputs = Vec::new();
    for id in &selected {
        let v = &data[*id as usize];
        for role in 0..=1 {
            let cs = explicit(v, role);
            let spec: CurveSpec = serde_json::from_value(cs.clone()).unwrap();
            spec.build().unwrap();
            fresh[role as usize].push(cs.clone());
            let w = rows
                .iter()
                .find(|w| w[1] == *id && w[2] == role && w[12] == 0)
                .unwrap();
            let mut isomorphic = cs.clone();
            let p = n(&cs, "p");
            isomorphic["kind"] = json!("prime_public_fixture");
            isomorphic["a"] = json!(mul(16, n(&cs, "a"), p));
            isomorphic["b"] = json!(mul(64, n(&cs, "b"), p));
            isomorphic["gx"] = json!(mul(4, n(&cs, "gx"), p));
            isomorphic["gy"] = json!(mul(8, n(&cs, "gy"), p));
            isomorphic["qx"] = json!(mul(4, w[9], p));
            isomorphic["qy"] = json!(mul(8, w[10], p));
            let ispec: CurveSpec = serde_json::from_value(isomorphic.clone()).unwrap();
            ispec.build().unwrap();
            model_records.push(json!({"p":p,"a":isomorphic["a"],"b":isomorphic["b"],"order":isomorphic["group_order"],"r":isomorphic["r"],"generator":[format!("0x{:x}",n(&isomorphic,"gx")),format!("0x{:x}",n(&isomorphic,"gy"))]}));
            public_inputs.push(json!({"parent_pair":id,"role":role,"coordinate_isomorphism_u":2,"source":cs,"target":isomorphic}));
            iso[role as usize].push(isomorphic);
        }
    }
    for role in 0..=1 {
        shards(
            &format!("fresh-{}", if role == 0 { "crater" } else { "floor" }),
            fresh[role].clone(),
            role == 0,
            false,
            2,
            q5,
            &q5hash,
        );
        shards(
            &format!("isomorphic-{}", if role == 0 { "crater" } else { "floor" }),
            iso[role].clone(),
            role == 0,
            false,
            1,
            q5,
            &q5hash,
        );
        let w = rows
            .iter()
            .find(|w| w[1] == 46 && w[2] == role as u64 && w[12] == 0)
            .unwrap();
        freeze(
            &format!("pilot-five-rounds-{role}"),
            vec![
                json!({"kind":"prime_public_fixture","p":w[3],"a":w[5],"b":w[6],"group_order":w[4],"r":w[4],"gx":w[7],"gy":w[8],"qx":w[9],"qy":w[10]}),
            ],
            role == 0,
            false,
            1,
            5,
            true,
            q5,
            &q5hash,
        );
    }
    let extra = lines(&format!("{ROOT}/results/extra-floor-certificates.jsonl"));
    let meta = lines(&format!("{ROOT}/prepared/extra-kernel-catalog.jsonl"));
    let mut seen: HashMap<u64, HashSet<String>> = data
        .iter()
        .map(|v| (n(v, "index"), HashSet::from([v["j1"].to_string()])))
        .collect();
    let mut extra_models = Vec::new();
    for f in extra {
        let id = n(&f, "id");
        let m = meta.iter().find(|m| n(m, "id") == id).unwrap();
        let parent = n(m, "parent_pair");
        assert!(selected.contains(&parent));
        assert!(
            seen.get_mut(&parent)
                .unwrap()
                .insert(f["j"].as_str().unwrap().into()),
            "repeated floor j vertex"
        );
        let cs = fromcert(&f);
        serde_json::from_value::<CurveSpec>(cs.clone())
            .unwrap()
            .build()
            .unwrap();
        model_records.push(json!({"p":cs["p"],"a":cs["a"],"b":cs["b"],"order":cs["group_order"],"r":cs["r"],"generator":[format!("0x{:x}",n(&cs,"gx")),format!("0x{:x}",n(&cs,"gy"))]}));
        extra_models.push(cs);
    }
    shards(
        "additional-floor",
        extra_models,
        false,
        false,
        1,
        q5,
        &q5hash,
    );
    let searched = lines(&format!("{ROOT}/prepared/search-certificates.jsonl"));
    let floors = lines(&format!("{ROOT}/results/floor-certificates-v2.jsonl"));
    let mut craters = Vec::new();
    let mut fmodels = Vec::new();
    for f in floors.iter().filter(|f| f["D"] == -43) {
        let cert = searched.iter().find(|c| c["id"] == f["id"]).unwrap();
        let c = fromcert(cert);
        let cs: CurveSpec = serde_json::from_value(c.clone()).unwrap();
        let inst = cs.build().unwrap();
        let crypto_lib::cryptanalysis::ecbench::workload::Instance::Prime(ref inst) = inst else {
            unreachable!()
        };
        derive_map(inst, 43, &mut GroupOps::default()).unwrap();
        let floor = fromcert(f);
        for v in [&c, &floor] {
            model_records.push(json!({"p":v["p"],"a":v["a"],"b":v["b"],"order":v["group_order"],"r":v["r"],"generator":[format!("0x{:x}",n(v,"gx")),format!("0x{:x}",n(v,"gy"))]}));
        }
        craters.push(c);
        fmodels.push(floor);
    }
    shards("d43-crater", craters, true, true, 1, q11, &q11hash);
    shards("d43-floor", fmodels, false, true, 1, q11, &q11hash);
    write(
        root.join("prepared/replication_curve_ids.json"),
        &json!(model_records),
    );
    write(
        root.join("prepared/replication-cohort.json"),
        &json!({"selected_pairs_in_original_shuffled_order":selected,"fresh_public_targets":40,"pilot_five_round_endpoints":2,"isomorphic_models":20,"distinct_additional_floor_vertices":20,"D43_pairs":10,"query_cap":1000000000u64,"cold_wall_watchdog_seconds":1800,"reference":"independent same-public-point rho for every IC attempt","q5_backend_sha256":q5hash,"q11_backend_sha256":q11hash,"isomorphism_checks":public_inputs}),
    );
    println!("froze 265 replication attempts across 30 specs; 60 additional models");
}

/// Register upper search models separately from the solve-ready 53-bit cohort.
pub fn upper_ids() {
    use num_bigint::BigUint;
    let mut models = Vec::new();
    let mut certificates = vec![
        (
            "prepared/search-certificates.jsonl",
            "crater; certified prime order",
        ),
        (
            "results/floor-certificates-v2.jsonl",
            "floor; independently certified",
        ),
    ];
    let upper_floor = "results/precomputed-96-floor-certificate.jsonl";
    if Path::new(ROOT).join(upper_floor).exists() {
        certificates.push((upper_floor, "floor; independently certified"));
    }
    for (path, status) in certificates {
        for c in lines(&format!("{ROOT}/{path}")) {
            let gx: BigUint = c["Gx"].as_str().unwrap().parse().unwrap();
            let gy: BigUint = c["Gy"].as_str().unwrap().parse().unwrap();
            models.push(json!({"p":c["p"],"a":c["a"],"b":c["b"],"order":c["N"],"r":c["N"],"generator":[format!("0x{gx:x}"),format!("0x{gy:x}")],"status":status,"source_certificate":path,"frobenius_conductor":c["ell"],"fundamental_discriminant":c["D"]}));
        }
    }
    write(
        Path::new(ROOT).join("prepared/upper_curve_ids.json"),
        &json!(models),
    );
}

pub fn controls(q5: &str, q11: &str) {
    let root = Path::new(ROOT);
    let mut records = Vec::new();
    let q11hash = sha256_hex(&fs::read(q11).unwrap());
    let cr = lines(&format!("{ROOT}/prepared/d43-control-crater.jsonl"))[0].clone();
    let fl = lines(&format!(
        "{ROOT}/results/d43-control-floor-certificate.jsonl"
    ))[0]
        .clone();
    for (role, c) in [(0, &cr), (1, &fl)] {
        let cs = fromcert(c);
        freeze(
            &format!("control-d43-{role}"),
            vec![cs.clone()],
            role == 0,
            true,
            1,
            2,
            false,
            q11,
            &q11hash,
        );
        let path = root.join(format!("prepared/replication/control-d43-{role}.json"));
        let mut spec = read(&path);
        for a in spec["arms"].as_array_mut().unwrap() {
            if a["name"] != "rho" {
                a["method"]["params"]["size"] = json!("64");
                a["method"]["params"]["query_limit"] = json!("500000");
            }
        }
        write(path, &spec);
        records.push(json!({"p":cs["p"],"a":cs["a"],"b":cs["b"],"order":cs["group_order"],"r":cs["r"],"generator":[format!("0x{:x}",n(&cs,"gx")),format!("0x{:x}",n(&cs,"gy"))]}));
    }
    let mut floor = read(root.join("prepared/control.json"));
    floor["workloads"]["curves"] = json!([{"kind":"prime_public_fixture","p":343073,"a":43182,"b":57354,"group_order":342833,"r":342833,"gx":2,"gy":20047,"qx":134265,"qy":91951}]);
    floor["arms"]
        .as_array_mut()
        .unwrap()
        .retain(|a| a["name"] != "velu");
    floor["label"] = json!("D=-19 native floor correctness control");
    write(root.join("prepared/control-floor.json"), &floor);
    for cs in [
        read(root.join("prepared/control.json"))["workloads"]["curves"][0].clone(),
        floor["workloads"]["curves"][0].clone(),
    ] {
        records.push(json!({"p":cs["p"],"a":cs["a"],"b":cs["b"],"order":cs["group_order"],"r":cs["r"],"generator":[format!("0x{:x}",n(&cs,"gx")),format!("0x{:x}",n(&cs,"gy"))]}));
    }
    write(
        root.join("prepared/control_curve_ids.json"),
        &json!(records),
    );
    assert_eq!(
        sha256_hex(&fs::read(q5).unwrap()),
        floor["arms"][0]["method"]["params"]["backend_sha256"]
            .as_str()
            .unwrap_or("")
    );
}

pub fn recurrence_inputs() {
    let root = Path::new(ROOT);
    let csv = fs::read_to_string(root.join("prepared/orbits.csv")).unwrap();
    let mut out = String::from("job,pair,role,p,r,a,b,Gx,Gy,Qx,Qy,seed,repeat\n");
    for line in csv.lines().skip(1) {
        let mut v: Vec<u64> = line.split(',').map(|x| x.parse().unwrap()).collect();
        v[0] = 6 * v[1] + v[2];
        out += &(v.iter().map(u64::to_string).collect::<Vec<_>>().join(",") + "\n");
        if v[2] == 0 || v[2] == 1 {
            v[2] = if v[2] == 0 { 4 } else { 5 };
            v[0] = 6 * v[1] + v[2];
            out += &(v.iter().map(u64::to_string).collect::<Vec<_>>().join(",") + "\n");
        }
    }
    fs::write(root.join("prepared/recurrence-orbits.csv"), out).unwrap();
}

pub fn torsion_inputs() {
    let root = Path::new(ROOT);
    let pairs = read(root.join("inputs/pairs.json"));
    let data = pairs.as_array().unwrap();
    let mut out = String::new();
    let append = |out: &mut String, v: Vec<String>| *out += &format!("[{}]\n", v.join(","));
    let val = |v: &Value| {
        v.as_str()
            .map(str::to_string)
            .unwrap_or_else(|| v.to_string())
    };
    for v in data {
        append(
            &mut out,
            vec![
                val(&v["index"]),
                val(&v["ell"]),
                val(&v["p"]),
                val(&v["N"]),
                "5".into(),
                val(&v["twist_a"]),
                val(&v["twist_b"]),
                val(&v["kernel_G"][0]),
                val(&v["kernel_G"][1]),
                val(&v["twist_d"]),
                val(&v["a0"]),
                val(&v["b0"]),
                val(&v["a1"]),
                val(&v["b1"]),
            ],
        );
    }
    let kernel_text = fs::read_to_string(root.join("prepared/new-kernels.txt")).unwrap()
        + &fs::read_to_string(root.join("prepared/d43-control-kernel.txt")).unwrap();
    let kernels: HashMap<u64, Vec<String>> = kernel_text
        .lines()
        .map(|l| {
            let a: Vec<String> = l.split_whitespace().map(str::to_string).collect();
            (a[0].parse().unwrap(), a)
        })
        .collect();
    for file in [
        "results/floor-certificates-v2.jsonl",
        "results/d43-control-floor-certificate.jsonl",
    ] {
        for f in lines(&format!("{ROOT}/{file}")) {
            let k = &kernels[&n(&f, "id")];
            let q = if f["D"] == -19 { "5" } else { "11" };
            append(
                &mut out,
                vec![
                    k[0].clone(),
                    k[1].clone(),
                    k[2].clone(),
                    val(&f["N"]),
                    q.into(),
                    k[3].clone(),
                    k[4].clone(),
                    k[5].clone(),
                    k[6].clone(),
                    k[7].clone(),
                    k[8].clone(),
                    k[9].clone(),
                    val(&f["a"]),
                    val(&f["b"]),
                ],
            );
        }
    }
    let meta = lines(&format!("{ROOT}/prepared/extra-kernel-catalog.jsonl"));
    for f in lines(&format!("{ROOT}/results/extra-floor-certificates.jsonl")) {
        let m = meta.iter().find(|m| m["id"] == f["id"]).unwrap();
        let p = &data[n(m, "parent_pair") as usize];
        assert_eq!(p["p"].to_string(), f["p"].as_str().unwrap());
        append(
            &mut out,
            vec![
                val(&f["id"]),
                val(&f["ell"]),
                val(&f["p"]),
                val(&f["N"]),
                "5".into(),
                val(&p["twist_a"]),
                val(&p["twist_b"]),
                val(&p["kernel_G"][0]),
                val(&p["kernel_G"][1]),
                val(&p["twist_d"]),
                val(&p["a0"]),
                val(&p["b0"]),
                val(&f["a"]),
                val(&f["b"]),
            ],
        );
    }
    fs::write(root.join("prepared/torsion-structure.gpvec"), out).unwrap();
}

pub fn recurrence_controls(backend: &str) {
    let root = Path::new(ROOT);
    let hash = sha256_hex(&fs::read(backend).unwrap());
    for (role, file) in [(0, "control.json"), (1, "control-floor.json")] {
        let mut s = read(root.join(format!("prepared/{file}")));
        let mut arm = s["arms"]
            .as_array()
            .unwrap()
            .iter()
            .find(|a| a["name"] == "scalar")
            .unwrap()
            .clone();
        arm["name"] = json!("recurrence");
        arm["role"] = json!("candidate");
        arm["method"]["id"] = json!("ic.cm_folded_recurrence_v1");
        arm["method"]["params"]["backend"] = json!(backend);
        arm["method"]["params"]["backend_sha256"] = json!(hash);
        let rho = s["arms"]
            .as_array()
            .unwrap()
            .iter()
            .find(|a| a["name"] == "rho")
            .unwrap()
            .clone();
        s["arms"] = json!([arm, rho]);
        s["label"] = json!(format!("Generic recurrence small control endpoint {role}"));
        write(
            root.join(format!("prepared/control-recurrence-{role}.json")),
            &s,
        );
    }
    for role in 0..=1 {
        let mut s = read(root.join(format!("prepared/shard-{}.json", 2 + 4 * role)));
        s["workloads"]["curves"]
            .as_array_mut()
            .unwrap()
            .retain(|v| n(v, "p") == 5947698030651343);
        assert_eq!(s["workloads"]["curves"].as_array().unwrap().len(), 1);
        let mut arm = s["arms"]
            .as_array()
            .unwrap()
            .iter()
            .find(|a| a["name"] == "scalar")
            .unwrap()
            .clone();
        arm["name"] = json!("recurrence");
        arm["role"] = json!("candidate");
        arm["method"]["id"] = json!("ic.cm_folded_recurrence_v1");
        arm["method"]["params"]["backend"] = json!(backend);
        arm["method"]["params"]["backend_sha256"] = json!(hash);
        let rho = s["arms"]
            .as_array()
            .unwrap()
            .iter()
            .find(|a| a["name"] == "rho")
            .unwrap()
            .clone();
        s["arms"] = json!([arm, rho]);
        s["label"] = json!(format!(
            "Generic recurrence pair-46 endpoint {role}, same primary public point and query seed"
        ));
        write(
            root.join(format!("prepared/recurrence-pilot-{role}.json")),
            &s,
        );
    }
}
