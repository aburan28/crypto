//! Native per-vertex input certification. No DLP solver or timing runs here.
use crypto_lib::cryptanalysis::{
    curve_traits::arith,
    ecbench::canonical::{canonical, sha256_hex},
    ic_boundary::{CountedGroup, GroupOps, PrimeCurve, PrimeInstance, PrimePoint},
};
use rusqlite::{Connection, OpenFlags};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    env,
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
};

fn main() {
    let args: Vec<_> = env::args().collect();
    assert_eq!(
        args.len(),
        4,
        "inventory REPOSITORY GRAPH_DIRECTORY FRESH_OUTPUT"
    );
    let repo = Path::new(&args[1]);
    let graphs = Path::new(&args[2]);
    let out = Path::new(&args[3]);
    fs::create_dir(out).expect("output must be fresh");
    let traits: Value =
        serde_json::from_slice(&fs::read(repo.join("docs/curves/traits.json")).unwrap()).unwrap();
    let mut summaries = vec![];
    for (label, expected_n, expected_sha) in [
        (
            "two",
            8192,
            "8ec54bf16456fd0cf6a62f8758a68425e0151157a44c11e6c2e00ddff67a8f64",
        ),
        (
            "three",
            6561,
            "a72ec162a221914f12e9c8806652e865b99565b7ddbcba860af0bae923de76b2",
        ),
        (
            "mixed",
            10070,
            "621d236a1a604ef391978aabb8b419f10db6daa54977831e37fb515d9bdc9d96",
        ),
    ] {
        let graph_path = graphs.join(format!("{label}.sqlite"));
        assert_eq!(sha256_hex(&fs::read(&graph_path).unwrap()), expected_sha);
        let study = repo.join("research/full_volcano_classes_20261017");
        let plan: Value = serde_json::from_slice(
            &fs::read(study.join(format!("evidence/{label}-plan.json"))).unwrap(),
        )
        .unwrap();
        let graph_audit: Value = serde_json::from_slice(
            &fs::read(study.join(format!("evidence/{label}-verify.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(graph_audit["vertex_class_complete"], true);
        assert_eq!(graph_audit["replay_receipt"]["status"], "REPLAYED");
        let p = plan["input"]["p"].as_u64().unwrap();
        let order = plan["order"].as_u64().unwrap();
        let source = plan["source_identity"]["slug"].as_str().unwrap();
        let trait_row = traits["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|v| v["slug"] == source)
            .unwrap();
        assert_eq!(trait_row["subgroup"]["status"], "proved");
        let r: u64 = trait_row["subgroup"]["order"]
            .as_str()
            .unwrap()
            .parse()
            .unwrap();
        assert!(arith::is_prime(&r.into()));
        assert_eq!(order % r, 0);
        let h = order / r;
        assert_ne!(
            h % r,
            0,
            "this inventory requires r to divide N exactly once"
        );
        let r1_factors = arith::factor(&(r - 1).into(), 1_000_000);
        assert!(r1_factors.complete());
        let k = arith::multiplicative_order(&p.into(), &r.into(), &r1_factors);
        let k: u64 = k.to_string().parse().unwrap();
        let expected_k: u64 = trait_row["embedding"]["degree"]
            .as_str()
            .unwrap()
            .parse()
            .unwrap();
        assert_eq!(k, expected_k);
        let ff = json!({"embedding_degree":k,"embedding_degree_status":"EXACT_NATIVE_MULTIPLICATIVE_ORDER",
            "r_minus_1_factorization":r1_factors.primes.iter().map(|(q,e)|json!([q.to_string(),e])).collect::<Vec<_>>(),
            "destination_field":"F_(p^k)","mathematical_reduction":"NOT_EXCLUDED_BY_LARGE_K",
            "existing_backend":"mov_attack_supersingular_k2; BSGS/PH destination solver, not finite-field IC",
            "finite_field_index_calculus":"UNIMPLEMENTED_FOR_THIS_FIELD_AND_DEGREE","pairing_execution":"NOT_RUN"});
        let db =
            Connection::open_with_flags(&graph_path, OpenFlags::SQLITE_OPEN_READ_ONLY).unwrap();
        let mut stmt = db
            .prepare("SELECT id,a,b,conductor,expanded FROM nodes ORDER BY id")
            .unwrap();
        let mut rows = stmt.query([]).unwrap();
        let mut file =
            BufWriter::new(File::create(out.join(format!("{label}-workloads.jsonl"))).unwrap());
        let mut populations = BTreeMap::<u64, u64>::new();
        let mut count = 0;
        while let Some(row) = rows.next().unwrap() {
            let id: u64 = row.get(0).unwrap();
            let a: u64 = row.get(1).unwrap();
            let b: u64 = row.get(2).unwrap();
            let conductor: u64 = row.get(3).unwrap();
            assert_eq!(row.get::<_, u64>(4).unwrap(), 1);
            let curve = PrimeCurve { p, a, b };
            let mut construction = GroupOps::default();
            let mut found = None;
            for x in 1..=10000 {
                if let Some(y) = curve.sqrt(curve.rhs(x)) {
                    let g = curve.mul(&mut construction, PrimePoint::affine(x, y), h);
                    if !g.infinity {
                        found = Some(g);
                        break;
                    }
                }
            }
            let g = found.expect("deterministic generator budget exhausted");
            assert!(curve.is_on_curve(g));
            let mut verification = GroupOps::default();
            assert!(curve.mul(&mut verification, g, r).infinity);
            let instance = PrimeInstance {
                name: label.into(),
                curve,
                group_order: order,
                r,
                cofactor: h,
                generator: (g.x, g.y),
            };
            let identity = instance.curve_id();
            let explicit = json!({"kind":"prime_explicit","p":p,"a":a,"b":b,"group_order":order,"r":r,"gx":g.x,"gy":g.y});
            let key = json!({"schema":"volcano-workload-input/v1","graph_sha256":expected_sha,"vertex":id,"conductor":conductor,"curve":explicit});
            let key_hash = sha256_hex(canonical(&key).unwrap().as_bytes());
            writeln!(file,"{}",json!({"input_id":format!("VWI1h{key_hash}"),"input":key,"model_slug":identity.slug,"icv1":identity.icv1,
                "ec1":null,"curve_uid":null,"representation_registration":"PENDING_EXACT_SUBGROUP_GENERATOR_REGISTRATION",
                "generator_certificate":{"on_curve":true,"r_prime":true,"r_times_G":"INFINITY","G":"NONZERO","order":r},
                "finite_field_applicability":ff,"direct_ic":"PENDING_BACKEND_BASE_VALIDATION_AND_EXECUTION",
                "solver_outcome":"NOT_RUN","relation_trials":null,"relations_found":null,"cold_time_ns":null,"total_operations":null})).unwrap();
            *populations.entry(conductor).or_default() += 1;
            count += 1;
        }
        file.flush().unwrap();
        assert_eq!(count, expected_n);
        assert_eq!(populations.len(), plan["strata"].as_array().unwrap().len());
        for s in plan["strata"].as_array().unwrap() {
            assert_eq!(
                populations[&s["conductor"].as_u64().unwrap()],
                s["expected_vertices"].as_u64().unwrap()
            );
        }
        let path = out.join(format!("{label}-workloads.jsonl"));
        summaries.push(json!({"case":label,"source_model":source,"vertices_certified":count,"strata_certified":populations.len(),"largest_prime_subgroup":r,
            "cofactor":h,"embedding_degree":k,"finite_field_applicability":ff,"graph_sha256":expected_sha,
            "workload_file":path.file_name().unwrap().to_str().unwrap(),"workload_bytes":fs::metadata(&path).unwrap().len(),
            "workload_sha256":sha256_hex(&fs::read(&path).unwrap()),"per_conductor_vertices":populations,"solver_runs":0}));
        eprintln!("{label}: {count} certified workload inputs; zero solver runs");
    }
    fs::write(out.join("inventory.json"),format!("{}\n",serde_json::to_string_pretty(&json!({"schema":"volcano-graphwide-inventory/v1","cases":summaries,"solver_execution":"NOT_RUN","graph_vertex_coverage":"ALL_THREE_COMPLETE_CLASSES"})).unwrap())).unwrap();
}
