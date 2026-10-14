//! Native supervision and integrity snapshots. Algebraic replay stays in isogeny.
use serde_json::{json, Value};
use std::{
    collections::HashSet,
    fmt::Write,
    fs,
    path::Path,
    process::{Command, Stdio},
    time::Instant,
};
#[path = "../../src/hash/sha256.rs"]
#[allow(dead_code)]
mod sha256;
#[allow(dead_code)]
mod visuals {
    use super::*;
    include!("../isogeny_catalogue_setup_20261009/visuals.rs");
    pub fn report(out: &Path, counts: &[Value], state: &str) {
        let mut ops = vec![
            Op::Text(42., 55., 28., "Additional public prime-degree maps".into()),
            Op::Text(
                42.,
                90.,
                18.,
                format!("Campaign {state}; target 2000 additional distinct maps"),
            ),
            Op::Text(
                42.,
                120.,
                15.,
                "Integrity snapshot of exact independent replay receipts; same host".into(),
            ),
        ];
        for (i, row) in counts.iter().enumerate() {
            let y = 205. + i as f64 * 200.;
            ops.push(Op::Text(
                42.,
                y,
                22.,
                format!(
                    "{}: prime degree {}",
                    row["curve"].as_str().unwrap(),
                    row["ell"]
                ),
            ));
            ops.push(Op::Text(
                42.,
                y + 35.,
                17.,
                format!(
                    "{} sealed steps; {} old maps excluded; {} additional",
                    row["verified_steps"], row["prior_duplicates"], row["additional_maps"]
                ),
            ));
            ops.push(Op::Rect(42., y + 55., 950., 35., "#e5eaf0".into()));
            ops.push(Op::Rect(
                42.,
                y + 55.,
                950. * row["additional_maps"].as_u64().unwrap() as f64 / 1000.,
                35.,
                "#247f67".into(),
            ));
            ops.push(Op::Text(
                42.,
                y + 115.,
                15.,
                "Bar: additional independently verified maps, target 1000 per class".into(),
            ));
        }
        ops.push(Op::Text(
            42.,
            660.,
            16.,
            "One Frobenius eigenline per step; coverage PARTIAL; two prime degrees".into(),
        ));
        ops.push(Op::Text(
            42.,
            695.,
            15.,
            "Proposed remaining steps are excluded from the verified count.".into(),
        ));
        fs::write(out.join("SEARCH.svg"), svg(&ops)).unwrap();
        let mut body = String::new();
        text(
            &mut body,
            42.,
            755.,
            22.,
            "At least 2000 additional prime-degree maps",
        );
        let mut y = 712.;
        for p in [format!("Status: {state}. The request stays open until at least 2000 additional distinct maps are sealed and delivered."),
            "Two fresh 1050-step walks start on the published P192/P224 models. The earlier 50-step prefix of each class is excluded by exact source, degree and monic-kernel identity.".into(),
            "Each sealed step passed independent kernel, codomain, exact rational curve-equation, subgroup order, and 20 public scalar-transport checks inside the dedicated executable.".into(),
            "This snapshot checks linked file hashes and replay receipt contents. It does not repeat the algebraic verification. All arithmetic replay is same-host and uses separate limb and polynomial code.".into(),
            "Evidence: PROTOCOL.md, SUMMARY.json, per-class walk.json, raw kernel/map JSON, replay.json, curves.json and per-walk exit receipts. Counts and operational phase sums are in SUMMARY.json.".into(),
            "The IC performance ledger remains unchanged. Canonical model graphs are updated only with archived verified findings; unfinished steps are not registered.".into()] {
            wrapped(&mut body, &mut y, &p, 11.);
        }
        emit_pdf(
            &out.join("REPORT.pdf"),
            vec![(612, 792, body), (842, 595, vector_pdf(&ops))],
        );
    }
}
fn hash(bytes: &[u8]) -> String {
    sha256::sha256(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
fn load(p: &Path) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn save(p: &Path, v: &Value) {
    let tmp = p.with_extension("pending");
    fs::write(
        &tmp,
        format!("{}\n", serde_json::to_string_pretty(v).unwrap()),
    )
    .unwrap();
    fs::rename(tmp, p).unwrap();
}
fn inventory(root: &Path, dir: &Path, rows: &mut Vec<Value>) {
    let mut paths: Vec<_> = fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .collect();
    paths.sort();
    for path in paths {
        // The live restart is durable, but mutable. Frozen inventories include
        // its archived checkpoints only after the campaign seals them.
        if path == root.join("runs") {
            continue;
        }
        if path.is_dir() {
            inventory(root, &path, rows);
        } else if path.file_name().unwrap() != "FILES.json" {
            let bytes = fs::read(&path).unwrap();
            rows.push(json!({"path":path.strip_prefix(root).unwrap().to_str().unwrap(),"bytes":bytes.len(),"sha256":hash(&bytes)}));
        }
    }
}
fn copy_tree(source: &Path, target: &Path) {
    fs::create_dir_all(target).unwrap();
    for entry in fs::read_dir(source).unwrap() {
        let source = entry.unwrap().path();
        let target = target.join(source.file_name().unwrap());
        if source.is_dir() {
            copy_tree(&source, &target);
        } else {
            fs::copy(source, target).unwrap();
        }
    }
}
fn check_long_prefix(binary: &Path, old: &Path, out: &Path) {
    assert!(!out.exists());
    fs::create_dir_all(out).unwrap();
    let walk = old.join("evidence/p224");
    let mut m = load(&walk.join("walk.json"));
    for index in 0..2 {
        let name = format!("step-{index:03}");
        copy_tree(&walk.join(&name), &out.join(name));
    }
    m["steps"].as_array_mut().unwrap().truncate(2);
    m["verified_maps"] = json!(2);
    m["requested_maps"] = json!(1050);
    m["status"] = json!("IN_PROGRESS");
    save(&out.join("walk.json"), &m);
    let result = Command::new(binary)
        .args(["verify-walk", "--input", out.to_str().unwrap()])
        .output()
        .unwrap();
    fs::write(out.join("stdout.json"), &result.stdout).unwrap();
    fs::write(out.join("stderr.txt"), &result.stderr).unwrap();
    assert!(result.status.success());
    let replay: Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(replay["verified_maps"], 2);
    assert_eq!(replay["requested_maps"], 1050);
    assert_eq!(replay["complete"], false);
    save(
        &out.join("CONTROL.json"),
        &json!({"status":"PASS","exit_code":result.status.code(),"assertion":"large request admits an exact verified prefix without claiming completion"}),
    );
}
fn identity(doc: &Value) -> String {
    let ell = doc["result"]["ell"].as_u64().unwrap();
    let maps = doc["result"]["isogenies"].as_array().unwrap();
    assert_eq!(maps.len(), 1);
    let kernel = maps[0]["kernel"].as_array().unwrap();
    assert_eq!(kernel.len() as u64, (ell - 1) / 2 + 1);
    let p = doc["curve"]["p"]
        .as_str()
        .expect("decimal field characteristic");
    assert!(
        !p.is_empty() && p.bytes().all(|c| c.is_ascii_digit()) && !p.starts_with('0'),
        "canonical positive characteristic required"
    );
    for coefficient in kernel {
        let s = coefficient.as_str().unwrap();
        assert!(
            !s.is_empty()
                && s.bytes().all(|c| c.is_ascii_digit())
                && (s == "0" || !s.starts_with('0')),
            "canonical decimal coefficient required"
        );
        assert!(
            s.len() < p.len() || (s.len() == p.len() && s < p),
            "coefficient must be reduced modulo the field characteristic"
        );
    }
    assert_eq!(kernel.last().unwrap(), "1", "kernel must be monic");
    hash(&serde_json::to_vec(&json!([doc["result"]["source_icv1"]["icv1"], ell, kernel])).unwrap())
}
fn raw(dir: &Path, curve: &str, ell: u64) -> std::path::PathBuf {
    dir.join(curve).join(format!("ell-{ell}.json"))
}
fn prior_keys(old: &Path) -> HashSet<String> {
    let mut keys = HashSet::new();
    for curve in ["p192", "p224"] {
        let dir = old.join("evidence").join(curve);
        let walk = load(&dir.join("walk.json"));
        assert_eq!(walk["status"], "PASS");
        let ell = walk["ell"].as_u64().unwrap();
        for row in walk["steps"].as_array().unwrap() {
            assert!(keys.insert(identity(&load(&raw(
                &dir.join(row["directory"].as_str().unwrap()),
                curve,
                ell
            )))));
        }
    }
    assert_eq!(keys.len(), 100);
    keys
}
fn snapshot(old: &Path, out: &Path) -> Value {
    let prior = prior_keys(old);
    let mut all = HashSet::new();
    let mut counts = vec![];
    let mut records = vec![];
    let mut total = 0;
    let mut complete = true;
    for (curve, ell) in [("p192", 2341), ("p224", 1117)] {
        let dir = out.join(curve);
        let mut additional = 0;
        let mut duplicates = 0;
        let mut construction = 0u64;
        let mut replay = 0u64;
        let mut invocation = 0u64;
        let mut steps = 0;
        let mut state = "NOT_STARTED".to_owned();
        if dir.join("walk.json").exists() {
            let m = load(&dir.join("walk.json"));
            assert_eq!(m["curve"], curve);
            assert_eq!(m["ell"], ell);
            assert_eq!(m["requested_maps"], 1050);
            assert_eq!(m["degree_coverage"], "PARTIAL");
            state = m["status"].as_str().unwrap().to_owned();
            let old_m = load(&old.join("evidence").join(curve).join("walk.json"));
            let mut source = old_m["steps"][0]["source_icv1"].clone();
            for (index, row) in m["steps"].as_array().unwrap().iter().enumerate() {
                assert_eq!(row["index"], index);
                assert_eq!(row["directory"], format!("step-{index:03}"));
                assert_eq!(row["source_icv1"], source);
                assert_eq!(row["status"], "PASS");
                let step = dir.join(row["directory"].as_str().unwrap());
                for (file, key) in [
                    ("model.json", "source_model_sha256"),
                    ("search.json", "search_sha256"),
                    ("replay.json", "replay_sha256"),
                    ("curves.json", "curves_sha256"),
                ] {
                    assert_eq!(row[key], hash(&fs::read(step.join(file)).unwrap()));
                }
                let bytes = fs::read(raw(&step, curve, ell)).unwrap();
                assert_eq!(row["map_sha256"], hash(&bytes));
                let doc: Value = serde_json::from_slice(&bytes).unwrap();
                let certificate = load(&step.join("replay.json"));
                assert_eq!(certificate["status"], "PASS");
                let r = &certificate["records"][0];
                for key in [
                    "status",
                    "kernel_check",
                    "codomain_check",
                    "subgroup_check",
                    "exact_rational_map_check",
                ] {
                    assert_eq!(r[key], "PASS");
                }
                assert_eq!(r["public_scalar_transport_checks"], 20);
                assert_eq!(r["certificate_sha256"], row["map_sha256"]);
                assert_eq!(r["source_icv1"], source);
                assert_eq!(r["target_icv1"], row["target_icv1"]);
                assert_eq!(r["ell"], ell);
                let search = load(&step.join("search.json"));
                let receipt = &search["attempts"][0]["receipt"];
                assert_eq!(receipt["status"], "PASS");
                assert_eq!(receipt["exit_code"], 0);
                assert_eq!(receipt["stdout_sha256"], row["map_sha256"]);
                assert_eq!(
                    receipt["stderr_sha256"],
                    hash(
                        &fs::read(step.join(curve).join(format!("ell-{ell}.stderr.txt"))).unwrap()
                    )
                );
                let id = identity(&doc);
                assert!(all.insert(id.clone()), "repeated map within campaign");
                let is_prior = prior.contains(&id);
                if is_prior {
                    duplicates += 1;
                } else {
                    additional += 1;
                    total += 1;
                }
                source = row["target_icv1"].clone();
                construction += row["construction_elapsed_ms_operational_only"]
                    .as_u64()
                    .unwrap();
                replay += row["verification_elapsed_ms"].as_u64().unwrap();
                invocation += row["invocation_elapsed_ms"].as_u64().unwrap();
                records.push(json!({"curve":curve,"step":index,"directory":format!("{curve}/step-{index:03}"),"normalized_map_identity":id,"previously_counted":is_prior,"independent_replay":r}));
                steps += 1;
            }
            assert_eq!(m["verified_maps"], steps);
        }
        let exit = out.join(format!("{curve}-exit.json"));
        let child_pass = exit.exists() && load(&exit)["exit_code"] == 0;
        complete &= state == "PASS" && steps == 1050 && child_pass;
        counts.push(json!({"curve":curve,"ell":ell,"requested_steps":1050,"verified_steps":steps,"additional_maps":additional,"prior_duplicates":duplicates,"construction_status":state,"child_exit_pass":child_pass,
            "construction_elapsed_ms_operational_only":construction,"independent_replay_elapsed_ms":replay,"summed_invocation_elapsed_ms":invocation}));
    }
    json!({"schema":"additional-prime-maps/v1","status":if complete && total>=2000 {"PASS"} else {"IN_PROGRESS"},"requested_additional_maps":2000,"verified_additional_maps":total,"remaining_additional_maps":2000u64.saturating_sub(total),"distinct_prime_degrees":2,"degree_coverage":"PARTIAL","integrity_audit":"PASS","audit_scope":"linked hashes and exact-replay receipts; algebraic replay occurs in the walk executable","counts":counts,"records":records})
}
fn launch(binary: &Path, out: &Path, protocol: &Path, old: &Path) {
    assert_eq!(prior_keys(old).len(), 100);
    assert!(!out.exists(), "launch needs a fresh output root");
    fs::create_dir_all(out).unwrap();
    fs::copy(binary, out.join("isogeny")).unwrap();
    fs::copy(protocol, out.join("PROTOCOL.md")).unwrap();
    let exe = out.join("isogeny").canonicalize().unwrap();
    let version = Command::new(&exe).arg("--version").output().unwrap();
    assert!(version.status.success());
    let v: Value = serde_json::from_slice(&version.stdout).unwrap();
    assert_eq!(v["tool"]["git_dirty"], "false");
    save(
        &out.join("LAUNCH.json"),
        &json!({"schema":"additional-prime-launch/v1","supervisor_pid":std::process::id(),"tool":v["tool"],"binary_sha256":hash(&fs::read(&exe).unwrap()),"supervisor_binary_sha256":hash(&fs::read(std::env::current_exe().unwrap()).unwrap()),"supervisor_source_commit":option_env!("CAMPAIGN_SOURCE_COMMIT").unwrap_or("unknown"),"protocol_sha256":hash(&fs::read(protocol).unwrap()),"prior_manifest_sha256":hash(&fs::read(old.join("MANIFEST.json")).unwrap()),"per_curve_count":1050,"concurrent_walks":2}),
    );
    let prevent_sleep = Command::new("caffeinate")
        .args(["-i", "-w", &std::process::id().to_string()])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .ok();
    let handles: Vec<_> = [("p192", "2341"), ("p224", "1117")].into_iter().map(|(curve, ell)| {
        let exe = exe.clone(); let out = out.to_owned();
        std::thread::spawn(move || {
            let args = vec!["walk".to_owned(),"--curve".into(),curve.into(),"--ell".into(),ell.into(),"--count".into(),"1050".into(),"--timeout".into(),"1800".into(),"--out".into(),out.join(curve).to_str().unwrap().into()];
            let start = Instant::now();
            let status = Command::new(exe).args(&args).stdout(fs::File::create(out.join(format!("{curve}.stdout.json"))).unwrap()).stderr(fs::File::create(out.join(format!("{curve}.stderr.log"))).unwrap()).status().unwrap();
            save(&out.join(format!("{curve}-exit.json")), &json!({"schema":"additional-prime-child-exit/v1","curve":curve,"args":args,"exit_code":status.code(),"elapsed_ms_operational_only":start.elapsed().as_millis(),"stdout_sha256":hash(&fs::read(out.join(format!("{curve}.stdout.json"))).unwrap()),"stderr_sha256":hash(&fs::read(out.join(format!("{curve}.stderr.log"))).unwrap())}));
            status.success()
        })
    }).collect();
    let ok = handles
        .into_iter()
        .map(|h| h.join().unwrap())
        .fold(true, |a, b| a & b);
    drop(prevent_sleep);
    save(
        &out.join("SUPERVISOR_EXIT.json"),
        &json!({"status":if ok {"PASS"} else {"FAIL"},"pid":std::process::id()}),
    );
    // Preserve the actual outcome even if no chat turn is active at completion.
    let final_summary = snapshot(old, out);
    save(&out.join("FINAL_SUMMARY.json"), &final_summary);
    if ok && final_summary["status"] != "PASS" {
        std::process::exit(4);
    }
    if !ok {
        std::process::exit(3);
    }
}
fn main() {
    let a: Vec<_> = std::env::args().collect();
    match a.get(1).map(String::as_str) {
        Some("manifest") => {
            assert_eq!(a.len(), 3, "manifest STUDY");
            let root = Path::new(&a[2]);
            let mut files = vec![];
            inventory(root, root, &mut files);
            save(
                &root.join("FILES.json"),
                &json!({"schema":"additional-prime-artifacts/v1","scope":"frozen study artifacts; live output is separate until archived","files":files}),
            );
        }
        Some("check-long-prefix") => {
            assert_eq!(
                a.len(),
                5,
                "check-long-prefix BINARY OLD_STUDY FRESH_OUTPUT"
            );
            check_long_prefix(Path::new(&a[2]), Path::new(&a[3]), Path::new(&a[4]));
        }
        Some("launch") => {
            assert_eq!(a.len(), 6, "launch BINARY FRESH_OUTPUT PROTOCOL OLD_STUDY");
            launch(
                Path::new(&a[2]),
                Path::new(&a[3]),
                Path::new(&a[4]),
                Path::new(&a[5]),
            );
        }
        Some("status") => {
            assert_eq!(a.len(), 4, "status OLD_STUDY OUTPUT");
            let mut v = snapshot(Path::new(&a[2]), Path::new(&a[3]));
            v.as_object_mut().unwrap().remove("records");
            println!("{}", serde_json::to_string_pretty(&v).unwrap());
        }
        Some("report") => {
            assert_eq!(a.len(), 5, "report OLD_STUDY OUTPUT FRESH_REPORT");
            let v = snapshot(Path::new(&a[2]), Path::new(&a[3]));
            let out = Path::new(&a[4]);
            assert!(!out.exists());
            fs::create_dir_all(out).unwrap();
            save(&out.join("SUMMARY.json"), &v);
            let mut md = format!("# Additional large prime-degree maps\n\nCampaign status: **{}**. Independently verified additional maps: **{} / 2000**. Remaining: **{}**. Two distinct prime degrees, with coverage PARTIAL at each step.\n\n| Class | Prime degree | Sealed steps | Prior maps excluded | Additional maps |\n|---|---:|---:|---:|---:|\n", v["status"].as_str().unwrap(),v["verified_additional_maps"],v["remaining_additional_maps"]);
            for r in v["counts"].as_array().unwrap() {
                writeln!(
                    md,
                    "| {} | {} | {} | {} | {} |",
                    r["curve"].as_str().unwrap(),
                    r["ell"],
                    r["verified_steps"],
                    r["prior_duplicates"],
                    r["additional_maps"]
                )
                .unwrap();
            }
            md.push_str("\nEach sealed step passed exact independent kernel, codomain, rational curve-equation and subgroup checks and 20 public scalar-transport checks. SUMMARY.json binds all exact ICV1/EC1 identities and normalized map keys. This integrity snapshot checks linked file hashes and replay receipts; algebraic replay occurs inside the dedicated walk executable before advancing. Repeated prefix work is included in phase sums but excluded from new-map yield.\n\n![Verified additional map counts](SEARCH.svg)\n\nThe original hundred-map evidence remains frozen. Raw evidence is in the launch directory until archived and committed. The requested 2000-map outcome remains open while the count is below 2000 or either walk is unfinished. Canonical graphs and models remain at the earlier archived campaign until new results are durably saved. No IC measurement or performance ledger row is promoted.\n");
            fs::write(out.join("RESULTS.md"), md).unwrap();
            visuals::report(
                out,
                v["counts"].as_array().unwrap(),
                v["status"].as_str().unwrap(),
            );
        }
        _ => panic!("use launch, status or report"),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn novelty_excludes_receipt_changes_but_detects_kernel_changes() {
        let mut doc = json!({"curve":{"p":"101"},"result":{"ell":5,"source_icv1":{"icv1":"ICV1:test-source"},"isogenies":[{"kernel":["2","3","1"]}]}});
        let first = identity(&doc);
        doc["elapsed_ms"] = json!(9000);
        assert_eq!(identity(&doc), first);
        doc["result"]["isogenies"][0]["kernel"][0] = json!("4");
        assert_ne!(identity(&doc), first);
        doc["result"]["source_icv1"]["icv1"] = json!("ICV1:other-source");
        assert_ne!(identity(&doc), first);
    }
    #[test]
    #[should_panic(expected = "kernel must be monic")]
    fn malformed_kernel_does_not_count() {
        identity(
            &json!({"curve":{"p":"101"},"result":{"ell":5,"source_icv1":{"icv1":"ICV1:test-source"},"isogenies":[{"kernel":["2","3","2"]}]}}),
        );
    }
    #[test]
    #[should_panic(expected = "canonical decimal coefficient required")]
    fn alternate_decimal_spelling_cannot_inflate_novelty() {
        identity(
            &json!({"curve":{"p":"101"},"result":{"ell":5,"source_icv1":{"icv1":"ICV1:test-source"},"isogenies":[{"kernel":["02","3","1"]}]}}),
        );
    }
    #[test]
    #[should_panic(expected = "coefficient must be reduced")]
    fn unreduced_coefficient_cannot_inflate_novelty() {
        identity(
            &json!({"curve":{"p":"101"},"result":{"ell":5,"source_icv1":{"icv1":"ICV1:test-source"},"isogenies":[{"kernel":["103","3","1"]}]}}),
        );
    }
}
