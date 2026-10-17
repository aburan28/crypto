//! Native write-once stage controller. No target generation or scientific analysis.
mod monitor;
use serde_json::{json, Value};
use std::{
    fs::{self, File, OpenOptions},
    io::Write,
    path::Path,
    process::Command,
};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
fn hash(p: &Path) -> Result<String> {
    let o = Command::new("/usr/bin/sha256sum")
        .arg("--")
        .arg(p)
        .output()?;
    if !o.status.success() {
        return Err(format!("hash failed {}", p.display()).into());
    }
    Ok(String::from_utf8(o.stdout)?
        .split_whitespace()
        .next()
        .ok_or("empty hash")?
        .into())
}
fn write_new(p: &Path, v: &Value) -> Result<()> {
    let mut f = OpenOptions::new().write(true).create_new(true).open(p)?;
    serde_json::to_writer_pretty(&mut f, v)?;
    writeln!(f)?;
    f.sync_all()?;
    Ok(())
}
fn read(p: &str) -> Result<Value> {
    Ok(serde_json::from_reader(File::open(p)?)?)
}
fn string(v: &Value, k: &str) -> Result<String> {
    Ok(v[k].as_str().ok_or_else(|| format!("missing {k}"))?.into())
}
fn plan(config: &str, out: &str) -> Result<()> {
    let c = read(config)?;
    let root = Path::new(out);
    fs::create_dir(root)?;
    let fixtures = string(&c, "training_fixtures")?;
    let ic = string(&c, "ic_binary")?;
    let rho = string(&c, "rho_binary")?;
    let mut bindings = std::collections::BTreeMap::<String, String>::new();
    for key in ["rho_config", "protocol", "candidate_manifest"] {
        let p = string(&c, key)?;
        bindings.insert(p.clone(), hash(Path::new(&p))?);
    }
    for candidate in c["candidates"].as_array().ok_or("candidates missing")? {
        let p = string(candidate, "base_spec")?;
        if !p.starts_with("construct:") {
            bindings.insert(p.clone(), hash(Path::new(&p))?);
        }
    }
    for name in ["generator_receipt.json", "targets.jsonl"] {
        let p = Path::new(&fixtures).join(name);
        bindings.insert(p.to_string_lossy().into(), hash(&p)?);
    }
    let mut rows = Vec::new();
    let mut add = |id: String,
                   method: &str,
                   candidate: Option<&Value>,
                   point_index: usize,
                   seed: u64|
     -> Result<()> {
        let points = Path::new(&fixtures).join(format!("target_{point_index:02}.jsonl"));
        let seeds = root.join(format!("{id}.seeds.json"));
        if method == "rho" {
            write_new(&seeds, &json!([seed]))?;
            bindings.insert(seeds.to_string_lossy().into(), hash(&seeds)?);
        }
        rows.push(json!({"id":id,"method":method,"candidate":candidate,"target_file":points,"target_sha256":hash(&points)?,"target_count":1,"rank_seed":if method=="ic"{Some(if point_index==0{29u64}else{43u64})}else{None},"seed_file":if method=="rho"{Some(seeds)}else{None}}));
        Ok(())
    };
    for pair in 0..2 {
        for repeat in 0..2 {
            add(
                format!("aa_{pair}_{repeat}"),
                "rho",
                None,
                0,
                202610156411 + pair,
            )?;
        }
    }
    for target in 0..2 {
        let mut cs = c["candidates"]
            .as_array()
            .ok_or("missing candidates")?
            .clone();
        if cs.len() != 5 {
            return Err("five candidates required".into());
        }
        if target == 1 {
            cs.reverse()
        };
        for candidate in &cs {
            add(
                format!(
                    "train_{target}_{}",
                    candidate["id"].as_str().ok_or("candidate id")?
                ),
                "ic",
                Some(candidate),
                target,
                0,
            )?;
        }
        add(
            format!("train_{target}_rho"),
            "rho",
            None,
            target,
            202610156801 + target as u64,
        )?;
    }
    for repeat in 0..2 {
        add(format!("aa_2_{repeat}"), "rho", None, 0, 202610156413)?;
    }
    let matrix = json!({"kind":"n83_basis_training_matrix","input_bindings":bindings,"config":c,"config_sha256":hash(Path::new(config))?,"ic_binary_sha256":hash(Path::new(&ic))?,"rho_binary_sha256":hash(Path::new(&rho))?,"controller_sha256":hash(&std::env::current_exe()?)?,"rows":rows,"processes":18,"answers":18,"holdouts_generated":false});
    write_new(&root.join("matrix.json"), &matrix)
}
fn run(matrix_path: &str, out: &str) -> Result<()> {
    let matrix = read(matrix_path)?;
    let c = &matrix["config"];
    let ic = string(c, "ic_binary")?;
    let rho = string(c, "rho_binary")?;
    if hash(Path::new(&ic))? != string(&matrix, "ic_binary_sha256")?
        || hash(Path::new(&rho))? != string(&matrix, "rho_binary_sha256")?
        || hash(&std::env::current_exe()?)? != string(&matrix, "controller_sha256")?
    {
        return Err("frozen binary mismatch".into());
    }
    for (p, h) in matrix["input_bindings"]
        .as_object()
        .ok_or("input bindings missing")?
    {
        if hash(Path::new(p))? != h.as_str().ok_or("binding hash missing")? {
            return Err(format!("frozen input mismatch {p}").into());
        }
    }
    let root = Path::new(out);
    fs::create_dir(root)?;
    write_new(&root.join("matrix.json"), &matrix)?;
    let mut failures = 0;
    for row in matrix["rows"].as_array().ok_or("rows missing")? {
        if monitor::was_interrupted() {
            return Err("interrupted; campaign not resumed".into());
        }
        let id = string(row, "id")?;
        let method = string(row, "method")?;
        let points = string(row, "target_file")?;
        if hash(Path::new(&points))? != string(row, "target_sha256")? {
            return Err("target hash mismatch".into());
        }
        let dir = root.join(&id);
        let answers = dir.join("answers.jsonl");
        let mut args: Vec<String> = vec![
            "--arm".into(),
            method.clone(),
            "--out-dir".into(),
            dir.to_string_lossy().into(),
            "--target-file".into(),
            points.clone(),
            "--target-count".into(),
            "1".into(),
            "--cwd".into(),
            string(c, "cwd")?,
            "--timeout-s".into(),
            "2100".into(),
            "--memory-bytes".into(),
            "12884901888".into(),
            "--cpu".into(),
            "3".into(),
            "--semantic-file".into(),
            answers.to_string_lossy().into(),
        ];
        let mut command = Vec::new();
        if method == "ic" {
            let candidate = &row["candidate"];
            let base = string(candidate, "base_spec")?;
            args.extend([
                "--expected-base-hash".into(),
                string(candidate, "base_hash")?,
            ]);
            for (k, v) in [
                (
                    "KIC_DUMP_BASE",
                    dir.join("base.jsonl").to_string_lossy().into_owned(),
                ),
                (
                    "KIC_DUMP_RANK",
                    dir.join("rank.jsonl").to_string_lossy().into_owned(),
                ),
                (
                    "KIC_PACKED_ROOT_TABLE",
                    if candidate["dense"] == true { "1" } else { "0" }.into(),
                ),
                (
                    "KIC_PACKED_ROOT_TABLE_DENSE",
                    if candidate["dense"] == true { "1" } else { "0" }.into(),
                ),
            ] {
                args.extend(["--env".into(), format!("{k}={v}")]);
            }
            command.extend([
                ic.clone(),
                base,
                points,
                row["rank_seed"]
                    .as_u64()
                    .ok_or("rank seed missing")?
                    .to_string(),
                answers.to_string_lossy().into(),
            ]);
        } else {
            command.extend([
                rho.clone(),
                points,
                string(row, "seed_file")?,
                string(c, "rho_config")?,
                answers.to_string_lossy().into(),
            ]);
        }
        args.push("--".into());
        args.extend(command);
        let status = monitor::launch(&args)?;
        if status != 0 {
            failures += 1;
        }
        println!(
            "{}",
            json!({"id":id,"status":status,"failures_so_far":failures})
        );
    }
    write_new(
        &root.join("stage_complete.json"),
        &json!({"status":"ALL_DECLARED_PROCESSES_ATTEMPTED","failures":failures,"winner":null,"independent_verification_pending":true}),
    )
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    let result = match a.first().map(String::as_str) {
        Some("plan-training") if a.len() == 3 => plan(&a[1], &a[2]),
        Some("run") if a.len() == 3 => run(&a[1], &a[2]),
        Some("run-arm") => monitor::launch(&a[1..]).and_then(|code| {
            if code == 0 {
                Ok(())
            } else {
                Err(format!("arm exit {code}").into())
            }
        }),
        _ => Err(
            "usage: runner plan-training CONFIG NEW_DIR | run MATRIX NEW_DIR | run-arm ARGS".into(),
        ),
    };
    if let Err(e) = result {
        eprintln!("ERROR: {e}");
        std::process::exit(1)
    }
}

#[cfg(test)]
mod plan_tests {
    use super::*;
    fn scratch() -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!(
            "n83-basis-plan-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&p).unwrap();
        p
    }
    #[test]
    fn training_matrix_is_complete_ordered_seeded_and_hash_bound() {
        let root = scratch();
        let points = root.join("fixtures");
        fs::create_dir(&points).unwrap();
        for name in ["target_00.jsonl", "target_01.jsonl", "targets.jsonl"] {
            fs::write(points.join(name), "[\"1\",\"2\"]\n").unwrap();
        }
        fs::write(points.join("generator_receipt.json"), "{}").unwrap();
        let artifact = root.join("artifact.json");
        fs::write(&artifact, "{}").unwrap();
        let exe = std::env::current_exe().unwrap();
        let candidates:Vec<_>=(0..5).map(|i|json!({"id":format!("C{i}"),"base_spec":artifact,"base_hash":"test-only","columns":600,"dense":i>0})).collect();
        let c = json!({"training_fixtures":points,"ic_binary":exe,"rho_binary":exe,"rho_config":artifact,"protocol":artifact,"candidate_manifest":artifact,"candidates":candidates,"cwd":root});
        let config = root.join("config.json");
        write_new(&config, &c).unwrap();
        let dir = root.join("plan");
        plan(config.to_str().unwrap(), dir.to_str().unwrap()).unwrap();
        let matrix_path = dir.join("matrix.json");
        let matrix = read(matrix_path.to_str().unwrap()).unwrap();
        let rows = matrix["rows"].as_array().unwrap();
        assert_eq!(rows.len(), 18);
        assert_eq!(rows.iter().filter(|r| r["method"] == "ic").count(), 10);
        assert_eq!(
            rows.iter()
                .filter(|r| r["id"].as_str().unwrap().starts_with("aa_"))
                .count(),
            6
        );
        assert_eq!(rows[4]["id"], "train_0_C0");
        assert_eq!(rows[10]["id"], "train_1_C4");
        assert_eq!(rows[4]["rank_seed"], 29);
        assert_eq!(rows[10]["rank_seed"], 43);
        assert!(
            plan(config.to_str().unwrap(), dir.to_str().unwrap()).is_err(),
            "write-once plan"
        );
        fs::write(&artifact, "changed").unwrap();
        let campaign = root.join("should_not_launch");
        assert!(
            run(matrix_path.to_str().unwrap(), campaign.to_str().unwrap())
                .unwrap_err()
                .to_string()
                .contains("frozen input mismatch")
        );
        assert!(!campaign.exists(), "reject before any measured launch");
        fs::remove_dir_all(&root).unwrap();
    }
}
