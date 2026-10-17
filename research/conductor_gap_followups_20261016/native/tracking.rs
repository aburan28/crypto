use super::*;
use std::{collections::BTreeMap, io::Write};
fn digest(path: impl AsRef<Path>) -> String {
    sha256_hex(&fs::read(path).unwrap())
}
pub fn freeze() {
    let root = Path::new(ROOT);
    assert!(!root.join("prepared/freeze.json").exists());
    for (file, count) in [
        ("dose-independent-certificate", 32),
        ("orbit-order-certificate", 100),
    ] {
        let text = fs::read_to_string(root.join(format!("results/{file}.jsonl"))).unwrap();
        let rows: Vec<Value> = text
            .lines()
            .map(|s| serde_json::from_str(s).unwrap())
            .collect();
        assert_eq!(rows.len(), count);
        assert!(rows.iter().all(|r| r["ok"] == true));
        let errors = fs::read_to_string(root.join(format!("results/{file}.stderr"))).unwrap();
        for forbidden in ["syntax error", "at top-level", "incorrect type", "Error"] {
            assert!(!errors.contains(forbidden));
        }
    }
    let mut files = BTreeMap::new();
    for file in fs::read_dir(root.join("prepared")).unwrap() {
        let p = file.unwrap().path();
        if p.is_file() {
            files.insert(p.to_string_lossy().into_owned(), digest(&p));
        }
    }
    for p in [
        "/workspace/native-deps/cm-followups-v1",
        "/workspace/native-deps/cm-folded-events-q5",
        "/workspace/native-deps/cm-folded-events-q11",
    ] {
        files.insert(p.to_string(), digest(p));
    }
    for dir in ["native", "protocols"] {
        for file in fs::read_dir(root.join(dir)).unwrap() {
            let p = file.unwrap().path();
            if p.is_file() {
                files.insert(p.to_string_lossy().into_owned(), digest(&p));
            }
        }
    }
    for p in ["group.hpp", "campaign_params.hpp"] {
        let p = Path::new(OLD).join("native").join(p);
        files.insert(p.to_string_lossy().into_owned(), digest(p));
    }
    for file in [
        "dose-independent-certificate.jsonl",
        "orbit-order-certificate.jsonl",
        "dose-construction-v3.jsonl",
    ] {
        let p = root.join("results").join(file);
        files.insert(p.to_string_lossy().into_owned(), digest(&p));
    }
    write(
        root.join("prepared/freeze.json"),
        &json!({"schema":"conductor-followup.freeze/v1","files":files,"sampling_seed":2026101601u64,"bootstrap_seed":2026101602u64,"overlap_holdout_seed":2026101603u64,"description":"Pre-run inputs, source and executable digests; one measured process at a time, pinned CPU 0, L0 local isolation only."}),
    );
}
pub fn instrumentation() {
    let root = Path::new(ROOT);
    let mut records = Vec::new();
    for control in [0, 2, 4] {
        let source = root.join(format!("runs/control-{control:04}"));
        let j = read(source.join("input.json"));
        let degree = u(&j, "degree");
        for repeat in 0..20 {
            let mut arms = Vec::new();
            let order = if repeat % 2 == 0 {
                ["noevents", "events"]
            } else {
                ["events", "noevents"]
            };
            for kind in order {
                let dir = root.join(format!(
                    "calibration/instrumentation-{control}-{repeat}-{kind}"
                ));
                fs::create_dir_all(dir.join("factor-bases")).unwrap();
                fs::create_dir(dir.join("relations")).unwrap();
                let start = Instant::now();
                let out =
                    Command::new(format!("/workspace/native-deps/cm-folded-{kind}-q{degree}"))
                        .arg(source.join("work.csv"))
                        .args(["64", "300000", "180", "1"])
                        .arg(&dir)
                        .arg(root.join("prepared/control-params.txt"))
                        .output()
                        .unwrap();
                let wall = start.elapsed().as_secs_f64();
                fs::write(dir.join("stdout.json"), &out.stdout).unwrap();
                fs::write(dir.join("stderr.log"), &out.stderr).unwrap();
                assert!(out.status.success());
                let v: Value = serde_json::from_slice(&out.stdout).unwrap();
                let mut rows: Vec<Value> = fs::read_to_string(dir.join("relations/job-0.jsonl"))
                    .unwrap()
                    .lines()
                    .map(|s| serde_json::from_str(s).unwrap())
                    .collect();
                for r in &mut rows {
                    r.as_object_mut().unwrap().remove("query");
                }
                arms.push((
                    kind,
                    v,
                    digest(dir.join("factor-bases/job-0.csv")),
                    rows,
                    wall,
                ));
            }
            for key in [
                "queries",
                "hits",
                "rank",
                "recovered",
                "verified",
                "pair_table_keys",
            ] {
                assert_eq!(arms[0].1[key], arms[1].1[key]);
            }
            assert_eq!(arms[0].2, arms[1].2);
            assert_eq!(arms[0].3, arms[1].3);
            records.push(json!({"control":control,"repeat":repeat,"equality_verified":true,"arms":arms.iter().map(|a|json!({"kind":a.0,"native":a.1,"base_sha256":a.2,"cold_wall_seconds":a.4})).collect::<Vec<_>>()}));
        }
    }
    write(
        root.join("results/instrumentation-controls.json"),
        &json!({"schema":"conductor-followup.instrumentation/v1","trials":records.len(),"records":records,"scope":"60 paired small-fixture process comparisons; exact arithmetic/row/recovery equality; wall diagnostics include process launch and parsing, not a production overhead estimate"}),
    );
}
fn check_freeze() {
    let f = read(Path::new(ROOT).join("prepared/freeze.json"));
    for (p, h) in f["files"].as_object().unwrap() {
        assert_eq!(digest(p), h.as_str().unwrap(), "frozen file changed: {p}");
    }
}
fn journal(row: &Value) {
    let p = Path::new(ROOT).join("results/execution.jsonl");
    let mut f = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(p)
        .unwrap();
    writeln!(f, "{}", serde_json::to_string(row).unwrap()).unwrap();
    f.sync_all().unwrap();
}
fn timestamp() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs()
}
fn local_guard() {
    let status = fs::read_to_string("/proc/self/status").unwrap();
    assert!(
        status.lines().any(|line| line == "Cpus_allowed_list:\t0"),
        "launch the campaign through taskset -c 0"
    );
    let own = std::process::id();
    for entry in fs::read_dir("/proc").unwrap().filter_map(Result::ok) {
        let Ok(pid) = entry.file_name().to_string_lossy().parse::<u32>() else {
            continue;
        };
        if pid == own {
            continue;
        }
        let Ok(comm) = fs::read_to_string(entry.path().join("comm")) else {
            continue;
        };
        let comm = comm.trim();
        assert!(
            ![
                "rustc",
                "cargo",
                "cc1plus",
                "cc1",
                "cm-folded",
                "cm-folded-event",
                "cm-followups-v1"
            ]
            .contains(&comm)
                && !comm.starts_with("crypto_lib-"),
            "local compiler/solver overlaps launch: pid {pid} {comm}"
        );
    }
}
fn selection() {
    let root = Path::new(ROOT);
    let mut scores = BTreeMap::<u64, u64>::new();
    for index in 0..800 {
        let dir = root.join(format!("runs/orbit-{index:04}"));
        let j = read(dir.join("input.json"));
        if u(&j, "B") != 8192 {
            continue;
        }
        let v = read(dir.join("result.json"));
        let keys = u(&v["result"]["native"], "pair_table_keys");
        let score = 8192u64.pow(2) - keys;
        *scores.entry(u(&j, "pair")).or_default() += score;
    }
    assert_eq!(scores.len(), 100);
    let mut ranked: Vec<_> = scores.iter().map(|(&p, &s)| (p, s)).collect();
    ranked.sort_by_key(|&(p, s)| (std::cmp::Reverse(s), seed(&format!("overlap-tie/{p}"))));
    let selected: Vec<_> = ranked[..10].iter().map(|x| x.0).collect();
    let mut heldout: Vec<_> = (0..100u64).filter(|p| !selected.contains(p)).collect();
    heldout.shuffle(&mut StdRng::seed_from_u64(2026101603));
    heldout.truncate(10);
    let pairs = read(Path::new(OLD).join("inputs/pairs.json"));
    let mut jobs = Vec::new();
    for (&pair, kind) in selected
        .iter()
        .map(|p| (p, "selected"))
        .chain(heldout.iter().map(|p| (p, "heldout")))
    {
        for role in 0..2 {
            let v = endpoint(&pairs[pair as usize], role);
            let i = inst(&v);
            let q = target(&i.curve, i.r, &format!("overlap/{pair}/{role}"));
            let mut j = job(
                "overlap",
                pair,
                0,
                if role == 0 {
                    "crater-recurrence"
                } else {
                    "floor-recurrence"
                },
                &v,
                q,
                seed(&format!("overlap-query/{pair}")),
                0,
            );
            j["selection"] = json!(kind);
            jobs.push(j);
        }
    }
    jobs.shuffle(&mut StdRng::seed_from_u64(2026101604));
    write(root.join("prepared/overlap.json"), &json!(jobs));
    write(
        root.join("prepared/overlap-selection.json"),
        &json!({"ranked":ranked,"selected":selected,"heldout":heldout,"jobs_sha256":digest(root.join("prepared/overlap.json")),"rule":"Ten highest summed pair-table abscissa collisions at length 2048; deterministic hash tie-break. Ten randomly drawn remaining gaps are a held-out comparison, not an independent population replication."}),
    );
}
pub fn campaign() {
    let root = Path::new(ROOT);
    check_freeze();
    local_guard();
    let mut host = BTreeMap::new();
    for p in [
        "/proc/self/status",
        "/proc/pressure/cpu",
        "/proc/pressure/memory",
        "/sys/fs/cgroup/cpu.max",
        "/sys/fs/cgroup/cpu.stat",
        "/sys/fs/cgroup/memory.max",
        "/sys/fs/cgroup/memory.events",
        "/sys/fs/cgroup/cpuset.cpus.effective",
    ] {
        if let Ok(text) = fs::read_to_string(p) {
            host.insert(p, text);
        }
    }
    write(
        root.join("results/host-start.json"),
        &json!({"pid":std::process::id(),"unix_seconds":timestamp(),"local_observations":host,"scope":"Container-visible local launch guard; no physical-host or frequency isolation attestation."}),
    );
    let lock = root.join("results/campaign.lock");
    let mut guard = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&lock)
        .expect("campaign lock exists: inspect process before removing");
    writeln!(guard, "{}", std::process::id()).unwrap();
    journal(&json!({"event":"campaign-start","pid":std::process::id(),"unix_seconds":timestamp()}));
    for cohort in [
        "crossover",
        "yield",
        "rho",
        "rho-coordinate",
        "orbit",
        "dose",
        "overlap",
    ] {
        if cohort == "overlap" && !root.join("prepared/overlap.json").exists() {
            selection()
        }
        let jobs = read(root.join(format!("prepared/{cohort}.json")));
        for index in 0..jobs.as_array().unwrap().len() {
            let label = format!("{cohort}-{index:04}");
            let receipt = root.join(format!("results/receipts/{label}.json"));
            if receipt.exists() {
                continue;
            }
            check_freeze();
            local_guard();
            let logs = root.join("results/processes");
            fs::create_dir_all(&logs).unwrap();
            fs::create_dir_all(root.join("results/receipts")).unwrap();
            let timefile = logs.join(format!("{label}-time.json"));
            let out = fs::File::create(logs.join(format!("{label}.stdout"))).unwrap();
            let err = fs::File::create(logs.join(format!("{label}.stderr"))).unwrap();
            journal(
                &json!({"event":"run-start","cohort":cohort,"index":index,"unix_seconds":timestamp()}),
            );
            let status=Command::new("/workspace/native-deps/root/usr/bin/time").arg("-f").arg("{\"user_seconds\":%U,\"sys_seconds\":%S,\"wall_seconds\":%e,\"max_rss_kib\":%M,\"exit_code\":%x}").arg("-o").arg(&timefile).arg("timeout").args(["--signal=TERM","--kill-after=10s","900s"]).arg("/workspace/native-deps/cm-followups-v1").arg("run-one").arg(cohort).arg(index.to_string()).stdout(out).stderr(err).status().unwrap();
            let result = root.join(format!("runs/{label}/result.json"));
            let record = json!({"schema":"conductor-followup.execution/v1","cohort":cohort,"index":index,"exit_code":status.code(),"process_time":if timefile.exists(){json!(fs::read_to_string(&timefile).unwrap())}else{Value::Null},"input_sha256":sha256_hex(serde_json::to_string(&jobs[index]).unwrap().as_bytes()),"result_sha256":if result.exists(){json!(digest(&result))}else{Value::Null},"unix_seconds":timestamp(),"completed":status.success()&&result.exists()});
            write(&receipt, &record);
            journal(&record);
        }
    }
    journal(&json!({"event":"campaign-finished","unix_seconds":timestamp()}));
    fs::remove_file(lock).unwrap();
}
