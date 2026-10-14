//! Thin native dispatcher: all measurement and auditing belongs to ecbench.
use serde_json::{json, Value};
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::{SystemTime, UNIX_EPOCH},
};

fn execute(binary: &Path, args: &[String], log: &Path) -> Result<Value, String> {
    let out = fs::File::create(log).map_err(|e| e.to_string())?;
    let started = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| e.to_string())?
        .as_secs();
    let status = Command::new(binary)
        .args(args)
        .stdout(Stdio::from(out.try_clone().map_err(|e| e.to_string())?))
        .stderr(Stdio::from(out))
        .status()
        .map_err(|e| e.to_string())?;
    Ok(
        json!({"binary":binary,"arguments":args,"log":log,"started_unix_utc":started,"exit_code":status.code(),"success":status.success()}),
    )
}

fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let [binary, study] = args.as_slice() else {
        return Err("usage: ecbench_open_questions_run ECBENCH STUDY_DIR".into());
    };
    let binary = fs::canonicalize(binary).map_err(|e| e.to_string())?;
    let root = PathBuf::from(study);
    let manifest = root.join("execution-manifest.json");
    if manifest.exists() {
        return Err("refusing to repeat/overwrite a campaign; freeze a new study directory".into());
    }
    fs::create_dir_all(root.join("sessions")).map_err(|e| e.to_string())?;
    fs::create_dir_all(root.join("comparisons")).map_err(|e| e.to_string())?;
    let mut specs: Vec<_> = fs::read_dir(root.join("specs"))
        .map_err(|e| e.to_string())?
        .map(|p| p.map(|p| p.path()))
        .collect::<Result<_, _>>()
        .map_err(|e| e.to_string())?;
    specs.retain(|p| p.extension().is_some_and(|e| e == "json"));
    specs.sort();
    let mut records = Vec::new();
    for spec in specs {
        let name = spec.file_stem().unwrap().to_str().unwrap();
        let dir = root.join("sessions").join(name);
        if dir.exists() {
            return Err(format!(
                "refusing to touch existing session {}",
                dir.display()
            ));
        }
        let doc: Value = serde_json::from_slice(&fs::read(&spec).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;
        let cmd = vec![
            "run".into(),
            "--spec".into(),
            spec.to_string_lossy().into(),
            "--out".into(),
            dir.to_string_lossy().into(),
            "--cpus".into(),
            "none".into(),
            "--allow-busy".into(),
            "--wait".into(),
        ];
        println!("starting {name}");
        records.push(execute(
            &binary,
            &cmd,
            &root.join(format!("{name}.run.log")),
        )?);
        if dir.join("session.json").exists() {
            let cmd = vec![
                "verify".into(),
                "--dir".into(),
                dir.to_string_lossy().into(),
                "--replay-all".into(),
                "--exit-code".into(),
                "--out".into(),
                root.join(format!("{name}.audit.json"))
                    .to_string_lossy()
                    .into(),
            ];
            records.push(execute(
                &binary,
                &cmd,
                &root.join(format!("{name}.verify.log")),
            )?);
            records.push(execute(
                &binary,
                &["table".into(), "--dir".into(), dir.to_string_lossy().into()],
                &root.join(format!("{name}.table.md")),
            )?);
            records.push(execute(
                &binary,
                &[
                    "resources".into(),
                    "--dir".into(),
                    dir.to_string_lossy().into(),
                    "--out".into(),
                    root.join(format!("{name}.resources.json"))
                        .to_string_lossy()
                        .into(),
                ],
                &root.join(format!("{name}.resources.log")),
            )?);
            let has_rho = doc["arms"]
                .as_array()
                .unwrap()
                .iter()
                .any(|a| a["name"] == "rho");
            if has_rho {
                for arm in doc["arms"].as_array().unwrap() {
                    let arm = arm["name"].as_str().unwrap();
                    if arm == "rho" {
                        continue;
                    }
                    records.push(execute(
                        &binary,
                        &[
                            "compare".into(),
                            "--dir".into(),
                            dir.to_string_lossy().into(),
                            "--a".into(),
                            "rho".into(),
                            "--b".into(),
                            arm.into(),
                            "--resamples".into(),
                            "4000".into(),
                            "--seed".into(),
                            "20261014".into(),
                            "--json".into(),
                        ],
                        &root
                            .join("comparisons")
                            .join(format!("{name}__rho__{arm}.json")),
                    )?);
                }
            }
        }
        fs::write(&manifest,format!("{}\n",serde_json::to_string_pretty(&json!({"schema":"ecbench-ten-questions-execution/v1","completed_through":name,"commands":records})).unwrap())).map_err(|e|e.to_string())?;
        println!("finished {name}; failures retained in manifest");
    }
    Ok(())
}
