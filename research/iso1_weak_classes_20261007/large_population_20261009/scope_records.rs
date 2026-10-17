//! Freeze exact quadratic-branch controls and canonical curve identities.
use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, path::Path};
fn sha(s: &[u8]) -> String {
    hex::encode(sha256(s))
}
fn kv(line: &str) -> BTreeMap<&str, &str> {
    line.split('|')
        .skip(1)
        .filter_map(|s| s.split_once('='))
        .collect()
}
fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    assert_eq!(args.len(), 1);
    let s = Path::new(&args[0]);
    let raw = fs::read_to_string(s.join("quadratic_scope_check.txt")).unwrap();
    assert!(fs::read(s.join("quadratic_scope_check.stderr"))
        .unwrap()
        .is_empty());
    assert!(raw.lines().last().unwrap().ends_with("|all_x_controls=3"));
    let meta = kv(raw.lines().next().unwrap());
    let modulus: Value = serde_json::from_str(meta["modulus"]).unwrap();
    let d: Value = serde_json::from_str(meta["d"]).unwrap();
    let mods = modulus
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect::<Vec<_>>()
        .join(",");
    let field = format!(
        "fpk-7-6-{}",
        &sha(format!("fpk-modulus:7:{mods}").as_bytes())[..8]
    );
    let mut records = Vec::new();
    let mut inputs = String::new();
    for line in raw.lines().filter(|l| l.starts_with("ALL_X|")) {
        let v = kv(line);
        assert_eq!(v["ellcard"], v["quartic_count"]);
        let a: Value = serde_json::from_str(v["short_a"]).unwrap();
        let b: Value = serde_json::from_str(v["short_b"]).unwrap();
        let j: Value = serde_json::from_str(v["j"]).unwrap();
        let alpha: Value = serde_json::from_str(v["alpha"]).unwrap();
        let model = json!({"a":a,"b":b,"field":field,"form":"y^2=x^3+a*x+b","k":"6","modulus":modulus,"p":"7","v":"1"});
        let canonical = serde_json::to_string(&model).unwrap();
        let hash = sha(canonical.as_bytes());
        let trace = v["trace"];
        let order = v["ellcard"];
        let signed = if let Some(t) = trace.strip_prefix('-') {
            format!("tm{t}")
        } else {
            format!("t{trace}")
        };
        let jj = j
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect::<Vec<_>>()
            .join(",");
        records.push(json!({"slug":format!("icv1-fp3k6-{signed}-{}",&hash[..8]),"icv1":format!("ICV1:{field}:{trace}:{order}:{jj}:unk:unk:r:{}",&hash[..12]),"model_json":canonical,"model_sha256":hash,"trace":trace,"group_order":order,"j_coefficients":j,"alpha":alpha,"base_field_nonsquare_d":d,"modulus_low_coefficients":modulus,"field_degree":6,"characteristic":"7","sample_number":v["sample"],"seed":meta["seed"],"fundamental_discriminant":v["D_K"],"frobenius_order_conductor":v["f_pi"],"endomorphism_order_conductor":null,"volcano_level":null,"subgroup_order":null,"quartic_model":"y^2=(x^2-d)*(x-alpha)*(x-alpha^49)","quartic_origin":"(alpha,0)","quartic_to_cubic_map":"X=c3/(x-alpha);Y=c3*y/(x-alpha)^2;c3=f'(alpha);c2=f''(alpha)/2;c1=f'''(alpha)/6","class_family":"Joux-Vitse quadratic h branch","verification":"ellcard plus complete affine-x count and two infinity points; explicit coordinate-map point replay"}));
        inputs.push_str(&format!(
            "scope_replay({}, {}, {trace}, {order}, {}, {}, {});\n",
            serde_json::to_string(&modulus).unwrap(),
            serde_json::to_string(&alpha).unwrap(),
            serde_json::to_string(&d).unwrap(),
            serde_json::to_string(&a).unwrap(),
            serde_json::to_string(&b).unwrap()
        ));
    }
    assert_eq!(records.len(), 3);
    fs::write(s.join("quadratic_curve_records.json"),serde_json::to_string_pretty(&json!({"schema":"iso1-quadratic-scope/v1","raw_sha256":sha(raw.as_bytes()),"script_sha256":sha(&fs::read(s.join("quadratic_scope_check.gp")).unwrap()),"input_law":"100 constructed quartics with a fixed nonsplit h; three exact all-x controls","records":records})).unwrap()+"\n").unwrap();
    fs::write(s.join("quadratic_inputs.gp"), inputs).unwrap();
    println!("three exact quadratic controls frozen with canonical identities");
}
