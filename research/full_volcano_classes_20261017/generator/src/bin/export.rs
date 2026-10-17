//! Export verified graph records and representative models without Python.
use rusqlite::{Connection, OpenFlags};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::{
    env,
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
    process::Command,
};

fn save(path: &Path, value: &Value) {
    fs::write(
        path,
        format!("{}\n", serde_json::to_string_pretty(value).unwrap()),
    )
    .unwrap();
}
fn main() {
    let args: Vec<_> = env::args().collect();
    assert_eq!(args.len(), 4, "export STUDY GRAPHDIR CLI");
    let study = Path::new(&args[1]);
    let graphs = Path::new(&args[2]);
    let mut counts = vec![];
    for label in ["two", "three", "mixed"] {
        let plan: Value = serde_json::from_slice(
            &fs::read(study.join(format!("evidence/{label}-plan.json"))).unwrap(),
        )
        .unwrap();
        let audit: Value = serde_json::from_slice(
            &fs::read(study.join(format!("evidence/{label}-verify.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(audit["vertex_class_complete"], true);
        assert_eq!(audit["selected_degree_edges_complete"], true);
        let db = Connection::open_with_flags(
            graphs.join(format!("{label}.sqlite")),
            OpenFlags::SQLITE_OPEN_READ_ONLY,
        )
        .unwrap();
        let dir = graphs.join(label);
        fs::create_dir_all(&dir).unwrap();
        let mut nodes = BufWriter::new(File::create(dir.join("nodes.csv")).unwrap());
        writeln!(nodes, "id,a,b,j,endomorphism_conductor").unwrap();
        let mut stmt = db
            .prepare("SELECT id,a,b,j,conductor FROM nodes ORDER BY id")
            .unwrap();
        let mut rows = stmt.query([]).unwrap();
        while let Some(row) = rows.next().unwrap() {
            let data = (0..5)
                .map(|i| row.get::<_, u64>(i).unwrap().to_string())
                .collect::<Vec<_>>();
            writeln!(nodes, "{}", data.join(",")).unwrap();
        }
        let mut edges = BufWriter::new(File::create(dir.join("edges.jsonl")).unwrap());
        let mut stratum_edges: BTreeMap<(u64, u64, u64), u64> = BTreeMap::new();
        let mut stmt = db.prepare("SELECT e.source,e.ell,e.kernel,e.target,e.u,s.conductor,t.conductor FROM edges e JOIN nodes s ON s.id=e.source JOIN nodes t ON t.id=e.target ORDER BY e.source,e.ell,e.kernel").unwrap();
        let mut rows = stmt.query([]).unwrap();
        while let Some(row) = rows.next().unwrap() {
            let source: u64 = row.get(0).unwrap();
            let ell: u64 = row.get(1).unwrap();
            let kernel: Value = serde_json::from_str(&row.get::<_, String>(2).unwrap()).unwrap();
            let target: u64 = row.get(3).unwrap();
            let u: u64 = row.get(4).unwrap();
            let sf: u64 = row.get(5).unwrap();
            let tf: u64 = row.get(6).unwrap();
            *stratum_edges.entry((sf, tf, ell)).or_default() += 1;
            let direction = if sf == tf {
                "horizontal"
            } else if tf == sf * ell {
                "descending"
            } else if sf == tf * ell {
                "ascending"
            } else {
                panic!("invalid conductor step")
            };
            writeln!(edges, "{}", json!({"source":source,"target":target,"degree":ell,"kernel_polynomial_ascending":kernel,"target_isomorphism_u":u,"direction":direction})).unwrap();
        }
        let stratum_edges: Vec<_> = stratum_edges.into_iter().map(|((source,target,degree),count)| json!({"source_conductor":source,"target_conductor":target,"degree":degree,"directed_maps":count})).collect();
        save(
            &study.join(format!("evidence/{label}-stratum-edges.json")),
            &json!(stratum_edges),
        );
        save(&dir.join("stratum-edges.json"), &json!(stratum_edges));
        let mut strata = vec![];
        let mut stmt = db.prepare("SELECT conductor,count(*),min(id) FROM nodes GROUP BY conductor ORDER BY conductor").unwrap();
        let data = stmt
            .query_map([], |r| {
                Ok((
                    r.get::<_, u64>(0)?,
                    r.get::<_, u64>(1)?,
                    r.get::<_, u64>(2)?,
                ))
            })
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(data.len(), plan["strata"].as_array().unwrap().len());
        for ((f, count, id), predicted) in data.iter().zip(plan["strata"].as_array().unwrap()) {
            assert_eq!(*f, predicted["conductor"].as_u64().unwrap());
            assert_eq!(*count, predicted["expected_vertices"].as_u64().unwrap());
            let result = Command::new(&args[3])
                .args([
                    "vertex",
                    "--input",
                    graphs.join(format!("{label}.sqlite")).to_str().unwrap(),
                    "--id",
                    &id.to_string(),
                ])
                .output()
                .unwrap();
            assert!(result.status.success());
            let vertex: Value = serde_json::from_slice(&result.stdout).unwrap();
            assert_eq!(vertex["endomorphism_conductor"], *f);
            save(
                &dir.join(format!("stratum-{f}-representative.json")),
                &vertex,
            );
            if *f == plan["frobenius_conductor"].as_u64().unwrap() {
                save(
                    &study.join(format!("inputs/{label}-floor.json")),
                    &vertex["model"],
                );
                let floor_input = serde_json::from_value(vertex["model"].clone()).unwrap();
                let floor_plan = volcano_cli::plan(floor_input, "auto").unwrap();
                save(
                    &study.join(format!("evidence/{label}-floor-plan.json")),
                    &serde_json::to_value(floor_plan).unwrap(),
                );
            }
            strata.push(json!({"conductor":f,"vertices":count,"representative_vertex":id,"representative":format!("{label}/stratum-{f}-representative.json")}));
        }
        save(
            &study.join(format!("evidence/{label}-strata.json")),
            &json!(strata),
        );
        save(&dir.join("strata.json"), &json!(strata));
        counts.push(json!({"case":label,"source_identity":plan["source_identity"],"field_p":plan["input"]["p"],"fundamental_discriminant":plan["fundamental_discriminant"],"maximum_conductor":plan["frobenius_conductor"],"vertices":audit["vertices"],"directed_maps":audit["directed_kernel_edges"],"degrees":audit["selected_degrees"],"strata":data.len(),"status":audit["status"],"independent_implementation_replay":"NOT_RUN","dlp_measurements":"NOT_RUN"}));
    }
    save(&study.join("evidence/coverage.json"), &json!(counts));
    println!(
        "{}",
        json!({"status":"EXPORTED_VERIFIED_CLASSES","cases":counts})
    );
}
