//! Reconstruct saved maps without generation seeds or graph search.
use std::{collections::BTreeMap, fmt::Write, fs, path::Path};
fn kv(s: &str) -> BTreeMap<&str, &str> {
    s.split('|')
        .skip(1)
        .filter_map(|v| v.split_once('='))
        .collect()
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 1);
    let base = Path::new(&a[0]);
    let mut gp = String::new();
    let mut routes = 0;
    let mut literals = 0;
    for run in [
        "construction_p7_run1",
        "construction_p7_odd_run3",
        "construction_p7_targeted_run1",
        "construction_p7_degree19_run1",
    ] {
        let path = base.join(run).join("receipt.json");
        if !path.exists() {
            continue;
        }
        let receipt: serde_json::Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        assert_eq!(receipt["status"], "PASS");
        let raw = fs::read_to_string(base.join(run).join("stdout.txt")).unwrap();
        let mut p = "";
        let mut modulus = "";
        let mut order = "";
        let mut ea = "";
        let mut eb = "";
        for line in raw.lines() {
            let m = kv(line);
            if line.starts_with("SOURCE|") {
                p = m["p"];
                modulus = m["modulus"];
                order = m["order"];
                writeln!(
                    gp,
                    "check_model({p},{modulus},{},{},{order});",
                    m["short_a"], m["short_b"]
                )
                .unwrap();
            }
            if line.starts_with("ENDPOINT|") {
                ea = m["short_a"];
                eb = m["short_b"];
                writeln!(gp, "check_model({p},{modulus},{ea},{eb},{order});").unwrap();
            }
            if line.starts_with("ROUTE|") {
                let (mode, kernel) = if let Some(poly) = m.get("kernel_polynomial") {
                    (1, *poly)
                } else {
                    (0, m["kernel"])
                };
                writeln!(
                    gp,
                    "check_route({p},{modulus},{},{},{},{},{kernel},{mode},{},{order});",
                    m["source_a"], m["source_b"], m["target_a"], m["target_b"], m["degree"]
                )
                .unwrap();
                routes += 1;
            }
            if line.starts_with("LITERAL_PLUS|") {
                writeln!(
                    gp,
                    "check_plus({p},{modulus},{ea},{eb},{},{},{},{order});",
                    m["alpha"], m["h"], m["origin"]
                )
                .unwrap();
                literals += 1;
            }
            if line.starts_with("LITERAL_MINUS|") {
                writeln!(
                    gp,
                    "check_minus({p},{modulus},{ea},{eb},{},{},{},{},{},{},{order});",
                    m["alpha"], m["d"], m["twist"], m["h"], m["c3"], m["a2"]
                )
                .unwrap();
                literals += 1;
            }
        }
    }
    writeln!(
        gp,
        "print(\"ALGEBRA_COMPLETE|saved_routes={routes}|literal_models={literals}\");"
    )
    .unwrap();
    fs::write(base.join("replay_inputs.gp"), gp).unwrap();
    println!("emitted {routes} routes and {literals} literal models");
}
