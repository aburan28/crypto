//! Public walk replay controls. Fixtures are copies of sealed campaign steps.
use serde_json::{json, Value};
use std::{fs, path::Path, process::Command};
#[path = "../../src/hash/sha256.rs"]
#[allow(dead_code)]
mod sha256;
fn hash(b: &[u8]) -> String {
    sha256::sha256(b)
        .iter()
        .map(|x| format!("{x:02x}"))
        .collect()
}
fn load(p: &Path) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn save(p: &Path, v: &Value) {
    fs::write(p, serde_json::to_vec_pretty(v).unwrap()).unwrap();
}
fn copy(a: &Path, b: &Path) {
    fs::create_dir_all(b).unwrap();
    for item in fs::read_dir(a).unwrap() {
        let p = item.unwrap().path();
        let q = b.join(p.file_name().unwrap());
        if p.is_dir() {
            copy(&p, &q);
        } else {
            fs::copy(p, q).unwrap();
        }
    }
}
fn invoke(binary: &str, dir: &Path, expected: bool) -> Value {
    let o = Command::new(binary)
        .args(["verify-walk", "--input", dir.to_str().unwrap()])
        .output()
        .unwrap();
    assert_eq!(
        o.status.success(),
        expected,
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let stdout: Value = serde_json::from_slice(&o.stdout).unwrap();
    fs::write(dir.join("control.stdout.json"), &o.stdout).unwrap();
    fs::write(dir.join("control.stderr.txt"), &o.stderr).unwrap();
    json!({"exit_code":o.status.code(),"stdout":stdout,"stderr_sha256":hash(&o.stderr)})
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        4,
        "walk-controls BINARY SEALED_WALK FRESH_OUTPUT"
    );
    let original = Path::new(&args[2]);
    let out = Path::new(&args[3]);
    assert!(!out.exists());
    fs::create_dir_all(out).unwrap();
    let mut m = load(&original.join("walk.json"));
    assert!(m["verified_maps"].as_u64().unwrap() >= 2);
    m["steps"].as_array_mut().unwrap().truncate(2);
    m["verified_maps"] = json!(2);
    m["requested_maps"] = json!(2);
    m["status"] = json!("PASS");
    let valid = out.join("valid");
    for index in 0..2 {
        copy(
            &original.join(format!("step-{index:03}")),
            &valid.join(format!("step-{index:03}")),
        );
    }
    save(&valid.join("walk.json"), &m);
    let pass = invoke(&args[1], &valid, true);
    let source = out.join("changed-source");
    copy(&valid, &source);
    let mut input = load(&source.join("step-001/model.json"));
    input["a"] = json!("1");
    save(&source.join("step-001/model.json"), &input);
    let mut changed = m.clone();
    changed["steps"][1]["source_model_sha256"] =
        json!(hash(&fs::read(source.join("step-001/model.json")).unwrap()));
    save(&source.join("walk.json"), &changed);
    let source_fail = invoke(&args[1], &source, false);
    let gen = out.join("changed-generator");
    copy(&valid, &gen);
    let mut target = load(&gen.join("step-000/curves.json"));
    target["curves"][0]["generator"][0] = json!("0x1");
    save(&gen.join("step-000/curves.json"), &target);
    let mut changed = m.clone();
    changed["steps"][0]["curves_sha256"] =
        json!(hash(&fs::read(gen.join("step-000/curves.json")).unwrap()));
    save(&gen.join("walk.json"), &changed);
    let generator_fail = invoke(&args[1], &gen, false);
    let empty = out.join("empty-chain");
    copy(&valid, &empty);
    let mut empty_manifest = m.clone();
    empty_manifest["steps"] = json!([]);
    empty_manifest["verified_maps"] = json!(0);
    empty_manifest["requested_maps"] = json!(0);
    save(&empty.join("walk.json"), &empty_manifest);
    let empty_fail = invoke(&args[1], &empty, false);
    save(
        &out.join("CONTROLS.json"),
        &json!({"schema":"isogeny-walk-controls/v1","status":"PASS",
        "positive_two_step_chain":pass,"changed_source_with_updated_hash_rejected":source_fail,
        "changed_mapped_generator_with_updated_hash_rejected":generator_fail,
        "empty_chain_rejected":empty_fail}),
    );
}
