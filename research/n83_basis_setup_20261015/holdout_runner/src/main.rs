//! Fresh-holdout controller, versioned separately from measured training code.
mod monitor;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
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
fn read(p: &Path) -> Result<Value> {
    Ok(serde_json::from_reader(File::open(p)?)?)
}
fn string(v: &Value, k: &str) -> Result<String> {
    Ok(v[k].as_str().ok_or_else(|| format!("missing {k}"))?.into())
}
fn new_json(p: &Path, v: &Value) -> Result<()> {
    let mut f = OpenOptions::new().write(true).create_new(true).open(p)?;
    serde_json::to_writer_pretty(&mut f, v)?;
    writeln!(f)?;
    f.sync_all()?;
    Ok(())
}
fn new_points(p: &Path, rows: &[Value]) -> Result<()> {
    let mut f = OpenOptions::new().write(true).create_new(true).open(p)?;
    for v in rows {
        serde_json::to_writer(&mut f, v)?;
        writeln!(f)?;
    }
    f.sync_all()?;
    Ok(())
}
fn public_points(p: &Path, n: usize) -> Result<Vec<Value>> {
    let text = fs::read_to_string(p)?;
    let rows: Vec<Value> = text
        .lines()
        .map(serde_json::from_str)
        .collect::<std::result::Result<_, _>>()?;
    if rows.len() != n {
        return Err("public target count mismatch".into());
    }
    let mut seen = BTreeSet::new();
    for row in &rows {
        let a = row
            .as_array()
            .filter(|a| a.len() == 2)
            .ok_or("public coordinates required")?;
        for v in a {
            let s = v.as_str().ok_or("decimal string coordinates required")?;
            let x: u128 = s.parse()?;
            if x >= (1u128 << 83) {
                return Err("coordinate exceeds field".into());
            }
        }
        if !seen.insert(row.to_string()) {
            return Err("duplicate public target".into());
        }
    }
    Ok(rows)
}
const COLD_SEEDS: [u64; 8] = [101, 103, 107, 109, 113, 127, 131, 137];
const ORDERS: [[&str; 3]; 6] = [
    ["control", "selected", "rho"],
    ["selected", "rho", "control"],
    ["rho", "control", "selected"],
    ["control", "rho", "selected"],
    ["rho", "selected", "control"],
    ["selected", "control", "rho"],
];
fn append(
    rows: &mut Vec<Value>,
    bindings: &mut BTreeMap<String, String>,
    root: &Path,
    id: String,
    stage: &str,
    role: &str,
    points: PathBuf,
    indices: Vec<usize>,
    candidate: Option<&Value>,
    rank_seed: Option<u64>,
    seeds: Vec<u64>,
) -> Result<()> {
    let seed_file = if role == "rho" {
        let p = root.join(format!("{id}.seeds.json"));
        new_json(&p, &json!(seeds))?;
        bindings.insert(p.to_string_lossy().into(), hash(&p)?);
        Some(p)
    } else {
        None
    };
    bindings.insert(points.to_string_lossy().into(), hash(&points)?);
    rows.push(json!({"id":id,"stage":stage,"role":role,"method":if role=="rho"{"rho"}else{"ic"},"candidate":candidate,"rank_seed":rank_seed,"target_count":indices.len(),"target_indices":indices,"target_file":points,"target_sha256":hash(&points)?,"seed_file":seed_file}));
    Ok(())
}
fn plan(config: &Path, root: &Path) -> Result<()> {
    let c = read(config)?;
    let selection = read(Path::new(&string(&c, "selection_receipt")?))?;
    if selection["status"] != "TRAINING_REVIEWED" || selection["holdouts_generated"] != false {
        return Err("selection is not frozen training-only evidence".into());
    }
    let selected = selection["selected_candidate"].as_str().unwrap_or("C0");
    let candidates = c["candidates"].as_array().ok_or("candidates missing")?;
    let control = candidates
        .iter()
        .find(|v| v["id"] == "C0")
        .ok_or("C0 missing")?;
    let chosen = candidates
        .iter()
        .find(|v| v["id"] == selected)
        .ok_or("selected candidate absent")?;
    let fixtures = PathBuf::from(string(&c, "holdout_fixtures")?);
    let points = public_points(&fixtures.join("targets.jsonl"), 8)?;
    let training = PathBuf::from(string(&c, "training_fixtures")?);
    let training_q = public_points(&training.join("target_00.jsonl"), 1)?;
    fs::create_dir(root)?;
    let mut bindings = BTreeMap::new();
    for key in [
        "protocol",
        "candidate_manifest",
        "selection_receipt",
        "selection_git_binding",
        "rho_config",
    ] {
        let p = string(&c, key)?;
        bindings.insert(p.clone(), hash(Path::new(&p))?);
    }
    for dir in [&fixtures, &training] {
        let p = dir.join("generator_receipt.json");
        bindings.insert(p.to_string_lossy().into(), hash(&p)?);
    }
    for candidate in [control, chosen] {
        let p = string(candidate, "base_spec")?;
        if !p.starts_with("construct:") {
            bindings.insert(p.clone(), hash(Path::new(&p))?);
        }
    }
    let mut rows = Vec::new();
    for target in 0..8 {
        for role in ORDERS[target % 6] {
            let candidate = if role == "control" {
                Some(control)
            } else if role == "selected" {
                Some(chosen)
            } else {
                None
            };
            append(
                &mut rows,
                &mut bindings,
                root,
                format!("cold_{target}_{role}"),
                "cold",
                role,
                fixtures.join(format!("target_{target:02}.jsonl")),
                vec![target],
                candidate,
                if role == "rho" {
                    None
                } else {
                    Some(COLD_SEEDS[target])
                },
                vec![202610156901 + target as u64],
            )?;
        }
    }
    // Pair3 ends the cold stage; its point is training target0, never a holdout.
    let aa_point = root.join("AA_training_target_00.jsonl");
    new_points(&aa_point, &training_q)?;
    for repeat in 0..2 {
        append(
            &mut rows,
            &mut bindings,
            root,
            format!("aa_3_{repeat}"),
            "cold",
            "rho",
            aa_point.clone(),
            vec![0],
            None,
            None,
            vec![202610156414],
        )?;
    }
    for cohort in 0..3 {
        let rotation = [0, 3, 5][cohort];
        let mut order: Vec<usize> = (0..8).map(|i| (i + rotation) % 8).collect();
        if cohort == 1 {
            order.reverse();
        }
        let file = root.join(format!("cohort_{cohort}.points.jsonl"));
        new_points(
            &file,
            &order.iter().map(|&i| points[i].clone()).collect::<Vec<_>>(),
        )?;
        for role in ORDERS[cohort] {
            let candidate = if role == "control" {
                Some(control)
            } else if role == "selected" {
                Some(chosen)
            } else {
                None
            };
            append(
                &mut rows,
                &mut bindings,
                root,
                format!("cohort_{cohort}_{role}"),
                "cohort",
                role,
                file.clone(),
                order.clone(),
                candidate,
                if role == "rho" {
                    None
                } else {
                    Some([149, 151, 157][cohort])
                },
                (0..8)
                    .map(|i| 202610157001 + 8 * cohort as u64 + i)
                    .collect(),
            )?;
        }
    }
    for repeat in 0..2 {
        append(
            &mut rows,
            &mut bindings,
            root,
            format!("aa_4_{repeat}"),
            "cohort",
            "rho",
            aa_point.clone(),
            vec![0],
            None,
            None,
            vec![202610156415],
        )?;
    }
    if rows.len() != 37
        || rows
            .iter()
            .map(|r| r["target_count"].as_u64().unwrap())
            .sum::<u64>()
            != 100
    {
        return Err("matrix dimensions violated".into());
    }
    let ic = string(&c, "ic_binary")?;
    let rho = string(&c, "rho_binary")?;
    new_json(
        &root.join("matrix.json"),
        &json!({"kind":"n83_basis_holdout_matrix/v1","config":c,"config_sha256":hash(config)?,"input_bindings":bindings,"ic_binary_sha256":hash(Path::new(&ic))?,"rho_binary_sha256":hash(Path::new(&rho))?,"controller_sha256":hash(&std::env::current_exe()?)?,"selected_candidate":selected,"processes":37,"answers":100,"rows":rows}),
    )
}
fn run(matrix_path: &Path, root: &Path, stage: &str) -> Result<()> {
    if !["cold", "cohort"].contains(&stage) {
        return Err("stage must be cold or cohort".into());
    }
    let m = read(matrix_path)?;
    let c = &m["config"];
    let ic = string(c, "ic_binary")?;
    let rho = string(c, "rho_binary")?;
    for (p, h) in m["input_bindings"].as_object().ok_or("bindings missing")? {
        if hash(Path::new(p))? != h.as_str().ok_or("binding hash missing")? {
            return Err(format!("frozen input mismatch {p}").into());
        }
    }
    if hash(Path::new(&ic))? != string(&m, "ic_binary_sha256")?
        || hash(Path::new(&rho))? != string(&m, "rho_binary_sha256")?
        || hash(&std::env::current_exe()?)? != string(&m, "controller_sha256")?
    {
        return Err("binary freeze mismatch".into());
    }
    fs::create_dir(root)?;
    new_json(&root.join("matrix.json"), &m)?;
    let mut failures = 0;
    let mut attempted = 0;
    for row in m["rows"]
        .as_array()
        .ok_or("rows missing")?
        .iter()
        .filter(|r| r["stage"] == stage)
    {
        if monitor::was_interrupted() {
            return Err("interrupted; no directory retry".into());
        }
        let id = string(row, "id")?;
        let method = string(row, "method")?;
        let n = row["target_count"].as_u64().ok_or("target count missing")?;
        let point = string(row, "target_file")?;
        let dir = root.join(&id);
        let answers = dir.join("answers.jsonl");
        let mut args: Vec<String> = vec![
            "--arm".into(),
            format!("{}{}", method, if n == 8 { "_batch" } else { "" }),
            "--out-dir".into(),
            dir.to_string_lossy().into(),
            "--target-file".into(),
            point.clone(),
            "--target-count".into(),
            n.to_string(),
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
        let command: Vec<String> = if method == "ic" {
            let candidate = &row["candidate"];
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
            vec![
                ic.clone(),
                string(candidate, "base_spec")?,
                point,
                row["rank_seed"]
                    .as_u64()
                    .ok_or("rank seed missing")?
                    .to_string(),
                answers.to_string_lossy().into(),
            ]
        } else {
            vec![
                rho.clone(),
                point,
                string(row, "seed_file")?,
                string(c, "rho_config")?,
                answers.to_string_lossy().into(),
            ]
        };
        args.push("--".into());
        args.extend(command);
        let status = monitor::launch(&args)?;
        attempted += 1;
        if status != 0 {
            failures += 1;
        }
        println!(
            "{}",
            json!({"id":id,"status":status,"failures_so_far":failures})
        );
    }
    new_json(
        &root.join("stage_complete.json"),
        &json!({"status":"ALL_DECLARED_STAGE_PROCESSES_ATTEMPTED","stage":stage,"attempted":attempted,"failures":failures,"winner":null,"independent_verification_pending":true}),
    )
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    let r = match a.first().map(String::as_str) {
        Some("plan") if a.len() == 3 => plan(Path::new(&a[1]), Path::new(&a[2])),
        Some("run") if a.len() == 4 => run(Path::new(&a[1]), Path::new(&a[2]), &a[3]),
        _ => {
            Err("usage: holdout_runner plan CONFIG NEW_DIR | run MATRIX NEW_DIR cold|cohort".into())
        }
    };
    if let Err(e) = r {
        eprintln!("ERROR: {e}");
        std::process::exit(1)
    }
}

#[cfg(test)]
mod protocol_tests {
    use super::*;
    #[test]
    fn complete_plan_preserves_counts_orders_seeds_and_bindings() {
        let root = std::env::temp_dir().join(format!(
            "n83-holdout-plan-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&root).unwrap();
        let fixtures = root.join("holdouts");
        let training = root.join("training");
        fs::create_dir(&fixtures).unwrap();
        fs::create_dir(&training).unwrap();
        let q: Vec<_> = (0..8)
            .map(|i| json!([(i + 1).to_string(), (i + 2).to_string()]))
            .collect();
        new_points(&fixtures.join("targets.jsonl"), &q).unwrap();
        for i in 0..8 {
            new_points(&fixtures.join(format!("target_{i:02}.jsonl")), &q[i..i + 1]).unwrap();
        }
        new_points(&training.join("target_00.jsonl"), &[json!(["100", "101"])]).unwrap();
        for dir in [&fixtures, &training] {
            new_json(&dir.join("generator_receipt.json"), &json!({})).unwrap();
        }
        let selection = root.join("selection.json");
        new_json(&selection,&json!({"status":"TRAINING_REVIEWED","selected_candidate":"C3","holdouts_generated":false})).unwrap();
        let artifact = root.join("artifact.json");
        new_json(&artifact, &json!({})).unwrap();
        let exe = std::env::current_exe().unwrap();
        let candidates = json!([{"id":"C0","base_spec":artifact,"base_hash":"control","columns":600,"dense":false},{"id":"C3","base_spec":"construct:83:1:1271","base_hash":"selected","columns":1271,"dense":true}]);
        let cfg = root.join("config.json");
        new_json(&cfg,&json!({"selection_receipt":selection,"selection_git_binding":artifact,"protocol":artifact,"candidate_manifest":artifact,"rho_config":artifact,"ic_binary":exe,"rho_binary":exe,"holdout_fixtures":fixtures,"training_fixtures":training,"candidates":candidates,"cwd":root})).unwrap();
        let output = root.join("plan");
        plan(&cfg, &output).unwrap();
        let matrix = read(&output.join("matrix.json")).unwrap();
        let rows = matrix["rows"].as_array().unwrap();
        assert_eq!(rows.len(), 37);
        assert_eq!(
            rows.iter()
                .map(|r| r["target_count"].as_u64().unwrap())
                .sum::<u64>(),
            100
        );
        assert_eq!(rows.iter().filter(|r| r["stage"] == "cold").count(), 26);
        assert_eq!(rows.iter().filter(|r| r["method"] == "ic").count(), 22);
        assert_eq!(rows[0]["id"], "cold_0_control");
        assert_eq!(rows[3]["id"], "cold_1_selected");
        assert_eq!(rows[0]["rank_seed"], 101);
        assert_eq!(rows[21]["rank_seed"], 137);
        assert_eq!(rows[26]["id"], "cohort_0_control");
        assert_eq!(rows[29]["target_indices"], json!([2, 1, 0, 7, 6, 5, 4, 3]));
        assert_eq!(rows[32]["target_indices"], json!([5, 6, 7, 0, 1, 2, 3, 4]));
        let seeds = read(Path::new(rows[32]["seed_file"].as_str().unwrap())).unwrap();
        assert_eq!(
            seeds,
            json!([
                202610157017u64,
                202610157018u64,
                202610157019u64,
                202610157020u64,
                202610157021u64,
                202610157022u64,
                202610157023u64,
                202610157024u64
            ])
        );
        assert!(plan(&cfg, &output).is_err());
        fs::write(&artifact, "changed").unwrap();
        let campaign = root.join("should-not-launch");
        assert!(run(&output.join("matrix.json"), &campaign, "cold")
            .unwrap_err()
            .to_string()
            .contains("frozen input mismatch"));
        assert!(!campaign.exists());
        fs::remove_dir_all(root).unwrap();
    }
}
