//! Prime-degree walks anchored in embedded public curves, replayed before advancing.
use crate::{catalog, provenance, replay, sha256_hex, standard_order, supervisor};
use serde_json::{json, Value};
use std::{collections::HashSet, fs, path::Path, time::Instant};

fn write(path: &Path, value: &Value) {
    let temp = path.with_extension("pending");
    fs::write(
        &temp,
        format!("{}\n", serde_json::to_string_pretty(value).unwrap()),
    )
    .unwrap();
    fs::rename(temp, path).unwrap();
}
fn anchor(name: &str) -> Value {
    json!({"curves":[catalog::resolve(name).unwrap().clone()]})
}
fn advance(registry: &mut Value, target: &Value) {
    let row = json!({"slug":target["name"],"icv1":target["icv1"],"family":"prime",
        "params":{"p":target["p"],"a":target["a"],"b":target["b"]},
        "order":target["group_order"], "representations":[{
        "ec1":target["target_ec1"],"curve_uid":target["target_curve_uid"],
        "curve":{"model":"short Weierstrass","generator":target["generator"],
        "subgroup_order":target["subgroup_order"],"cofactor":target["cofactor"]}}]});
    assert!(
        !registry["curves"]
            .as_array()
            .unwrap()
            .iter()
            .any(|r| r["slug"] == row["slug"]),
        "walk revisited a model; this is not a new source"
    );
    registry["curves"].as_array_mut().unwrap().push(row);
}
fn model(row: &Value) -> Value {
    json!({"p":row["params"]["p"],"a":row["params"]["a"],"b":row["params"]["b"]})
}
pub fn run(
    name: &str,
    ell: u64,
    count: u64,
    out: &Path,
    timeout: u64,
) -> Result<(Value, i32), String> {
    if !(1010..=1_000_000).contains(&ell)
        || crate::constructor::selection(name, ell).is_none()
        || !(1..=1000).contains(&count)
        || timeout == 0
    {
        return Err("walk needs P192/P224, a supported prime degree above 1009, count 1..1000, and positive timeout".into());
    }
    if out.exists() {
        return Err("walk --out must be a fresh directory".into());
    }
    fs::create_dir_all(out).map_err(|e| e.to_string())?;
    let exe = std::env::current_exe().map_err(|e| e.to_string())?;
    let digest = sha256_hex(&fs::read(&exe).map_err(|e| e.to_string())?);
    let mut registry = anchor(name);
    let mut seen = HashSet::new();
    let mut manifest = json!({"schema":"isogeny-prime-walk/v1","status":"IN_PROGRESS",
        "tool":provenance(),"executable_sha256":digest,"curve":name,"ell":ell,
        "requested_maps":count,"verified_maps":0,"degree_coverage":"PARTIAL",
        "protocol_sha256":sha256_hex(crate::KERNEL_PROTOCOL),"timeout_seconds":timeout,
        "resource_envelope":{"construction_rss_limit_bytes":8589934592u64,
            "construction_timeout_seconds":timeout,"walk_concurrency":1},
        "steps":[]});
    write(&out.join("walk.json"), &manifest);
    for index in 0..count {
        let start = Instant::now();
        let step_name = format!("step-{index:03}");
        let step = out.join(&step_name);
        let dir = step.join(name);
        fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
        let source = registry["curves"].as_array().unwrap().last().unwrap();
        let input = model(source);
        let source_slug = source["slug"].clone();
        write(&step.join("model.json"), &input);
        let args = vec![
            "kernel-first".into(),
            "--curve".into(),
            name.into(),
            "--ell".into(),
            ell.to_string(),
            "--order".into(),
            standard_order(name).to_string(),
            "--seed".into(),
            "1".into(),
            "--model".into(),
            step.join("model.json")
                .canonicalize()
                .unwrap()
                .to_string_lossy()
                .into_owned(),
        ];
        let receipt: Value = serde_json::from_str(
            &supervisor::run(&exe, &dir, &args, &format!("ell-{ell}"), timeout).dump(),
        )
        .unwrap();
        write(
            &step.join("search.json"),
            &json!({"schema":"isogeny-search/v1",
            "source_commit":env!("ISOGENY_GIT_COMMIT"),"tool":provenance(),
            "executable_sha256":digest,"reference_protocol_sha256":sha256_hex(crate::KERNEL_PROTOCOL),
            "timeout_seconds":timeout,"attempts":[{"curve":name,"ell":ell,"receipt":receipt}]}),
        );
        if receipt["status"] != "PASS" {
            manifest["status"] = receipt["status"].clone();
            manifest["failed_step"] = json!(step_name);
            write(&out.join("walk.json"), &manifest);
            return Ok((manifest, 3));
        }
        let verify_start = Instant::now();
        let certificate = replay::verify_with_registry(&step, true, &registry);
        let verification_ms = verify_start.elapsed().as_millis();
        assert_eq!(certificate["records"].as_array().unwrap().len(), 1);
        let record = &certificate["records"][0];
        assert_eq!(record["source_icv1"], source_slug);
        let key = format!("{}:{}", source_slug, record["certificate_sha256"]);
        assert!(seen.insert(key), "duplicate map");
        let targets: Value =
            serde_json::from_slice(&fs::read(step.join("curves.json")).unwrap()).unwrap();
        let target = &targets["curves"][0];
        advance(&mut registry, target);
        let row = json!({"index":index,"directory":step_name,"source_icv1":source_slug,
            "target_icv1":target["name"],"ell":ell,"status":"PASS",
            "source_model_sha256":sha256_hex(&fs::read(step.join("model.json")).unwrap()),
            "search_sha256":sha256_hex(&fs::read(step.join("search.json")).unwrap()),
            "replay_sha256":sha256_hex(&fs::read(step.join("replay.json")).unwrap()),
            "curves_sha256":sha256_hex(&fs::read(step.join("curves.json")).unwrap()),
            "map_sha256":record["certificate_sha256"],"construction_elapsed_ms_operational_only":receipt["elapsed_ms_operational_only"],
            "verification_elapsed_ms":verification_ms,"invocation_elapsed_ms":start.elapsed().as_millis()});
        manifest["steps"].as_array_mut().unwrap().push(row);
        manifest["verified_maps"] = json!(index + 1);
        write(&out.join("walk.json"), &manifest);
        eprintln!(
            "walk {name} ell={ell}: {}/{} distinct maps independently verified",
            index + 1,
            count
        );
    }
    manifest["status"] = json!("PASS");
    write(&out.join("walk.json"), &manifest);
    write(&out.join("verified-registry.json"), &registry);
    Ok((manifest, 0))
}
pub fn verify(out: &Path) -> Result<Value, String> {
    let manifest: Value =
        serde_json::from_slice(&fs::read(out.join("walk.json")).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;
    if manifest["schema"] != "isogeny-prime-walk/v1" {
        return Err("unsupported walk schema".into());
    }
    let name = manifest["curve"].as_str().ok_or("missing anchor")?;
    let ell = manifest["ell"].as_u64().ok_or("missing degree")?;
    if !(1010..=1_000_000).contains(&ell) || crate::constructor::selection(name, ell).is_none() {
        return Err("unsupported walk degree or anchor".into());
    }
    let steps = manifest["steps"].as_array().ok_or("missing steps")?;
    let requested = manifest["requested_maps"]
        .as_u64()
        .ok_or("missing requested count")?;
    if !(1..=1000).contains(&requested) || steps.is_empty() || steps.len() as u64 > requested {
        return Err(
            "walk replay requires a nonempty verified prefix within requested count 1..1000".into(),
        );
    }
    assert_eq!(
        steps.len() as u64,
        manifest["verified_maps"].as_u64().unwrap()
    );
    let mut registry = anchor(name);
    let mut maps = HashSet::new();
    for (index, row) in steps.iter().enumerate() {
        assert_eq!(row["index"], index);
        assert_eq!(row["directory"], format!("step-{index:03}"));
        let step = out.join(format!("step-{index:03}"));
        for (file, key) in [
            ("model.json", "source_model_sha256"),
            ("search.json", "search_sha256"),
            ("replay.json", "replay_sha256"),
            ("curves.json", "curves_sha256"),
        ] {
            assert_eq!(row[key], sha256_hex(&fs::read(step.join(file)).unwrap()));
        }
        let source = registry["curves"].as_array().unwrap().last().unwrap();
        assert_eq!(row["source_icv1"], source["slug"]);
        let input: Value =
            serde_json::from_slice(&fs::read(step.join("model.json")).unwrap()).unwrap();
        assert_eq!(input, model(source));
        let summary: Value =
            serde_json::from_slice(&fs::read(step.join("search.json")).unwrap()).unwrap();
        assert_eq!(summary["attempts"][0]["ell"], ell);
        assert_eq!(summary["attempts"][0]["curve"], name);
        let (certificate, expected_targets) = replay::verify_detail(&step, false, &registry);
        let retained: Value =
            serde_json::from_slice(&fs::read(step.join("replay.json")).unwrap()).unwrap();
        assert_eq!(certificate, retained);
        let record = &certificate["records"][0];
        assert_eq!(record["source_icv1"], row["source_icv1"]);
        assert_eq!(record["target_icv1"], row["target_icv1"]);
        assert_eq!(record["certificate_sha256"], row["map_sha256"]);
        assert!(maps.insert(format!(
            "{}:{}",
            record["source_icv1"], record["certificate_sha256"]
        )));
        let targets: Value =
            serde_json::from_slice(&fs::read(step.join("curves.json")).unwrap()).unwrap();
        // Bind the transported generator and model to exact replay.
        assert_eq!(targets, expected_targets);
        advance(&mut registry, &targets["curves"][0]);
    }
    if manifest["status"] == "PASS" {
        assert_eq!(manifest["requested_maps"], manifest["verified_maps"]);
    }
    Ok(
        json!({"schema":"isogeny-prime-walk-replay/v1","status":"PASS","curve":name,
        "ell":ell,"verified_maps":maps.len(),"degree_coverage":"PARTIAL",
        "requested_maps":requested,"construction_status":manifest["status"],
        "complete":manifest["status"]=="PASS" && maps.len() as u64==requested,
        "scope":"same-host independent implementation","walk_sha256":sha256_hex(&fs::read(out.join("walk.json")).unwrap())}),
    )
}
