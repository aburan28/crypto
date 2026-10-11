//! Emit independent GP reconstruction from the saved documentary identities.
use serde_json::Value;
use std::{collections::BTreeMap, fmt::Write, fs, path::Path};
fn json(v: &Value) -> String {
    serde_json::to_string(v).unwrap()
}
fn main() {
    let a = std::env::args().skip(1).collect::<Vec<_>>();
    assert_eq!(a.len(), 1);
    let base = Path::new(&a[0]);
    let mut input = String::new();
    let mut models = 0;
    let mut routes = 0;
    for run in [
        "validation_p7_run1",
        "pilot_run1",
        "positive_controls_run1",
        "population_run1",
    ] {
        let doc: Value =
            serde_json::from_slice(&fs::read(base.join(run).join("curve_records.json")).unwrap())
                .unwrap();
        for record in doc["records"].as_array().unwrap() {
            let p = record["characteristic"].as_str().unwrap();
            let k = record["field_degree"].as_u64().unwrap();
            let modulus = &record["modulus_low_coefficients"];
            let mut roles = BTreeMap::new();
            for curve in record["curves"].as_array().unwrap() {
                let role = curve["role"].as_str().unwrap();
                roles.insert(role, curve);
                let model: Value =
                    serde_json::from_str(curve["model_json"].as_str().unwrap()).unwrap();
                let weak = if role == "source" {
                    i32::from(record["direct_weak"].as_bool().unwrap())
                } else if role == "witness" || role == "control_known_weak" {
                    1
                } else {
                    -1
                };
                let four = if role == "full4" { 1 } else { -1 };
                writeln!(
                    input,
                    "verify_record({p},{k},{},{},{},{},{},{},{weak},{four});",
                    json(modulus),
                    json(&model["a"]),
                    json(&model["b"]),
                    json(&curve["j_coefficients"]),
                    json(&curve["original_root_coordinates"]),
                    curve["group_order"].as_str().unwrap()
                )
                .unwrap();
                models += 1;
            }
            for route in record["route"].as_array().unwrap() {
                let parent = route["parent"].as_str().unwrap();
                let child = route["child"].as_str().unwrap();
                let source_role = if parent == "1" {
                    if record["selected_trace_sign"].as_str() == Some("-1") {
                        "selected_twist".to_string()
                    } else {
                        "source".to_string()
                    }
                } else {
                    format!("route{parent}")
                };
                let target_role = format!("route{child}");
                let source = roles[source_role.as_str()];
                let target = roles[target_role.as_str()];
                writeln!(
                    input,
                    "verify_route({p},{k},{},{},{},{},{},{});",
                    json(modulus),
                    json(&source["original_root_coordinates"]),
                    json(&target["original_root_coordinates"]),
                    route["degree"].as_str().unwrap(),
                    json(&route["kernel_x"]),
                    source["group_order"].as_str().unwrap()
                )
                .unwrap();
                routes += 1;
            }
        }
    }
    writeln!(
        input,
        "print(\"REPLAY_COMPLETE|models={models}|routes={routes}\");"
    )
    .unwrap();
    fs::write(base.join("model_inputs.gp"), input).unwrap();
    println!("emitted {models} models and {routes} routes");
}
