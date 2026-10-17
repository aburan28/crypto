//! Native checkpointed orchestration; ecbench owns solving, counts and audits.
use crypto_lib::cryptanalysis::ecbench::canonical::{canonical, sha256_hex};
use crypto_lib::cryptanalysis::ecbench::{
    methods::MethodSpec,
    record::{self, ChildInput, CHILD_INPUT_SCHEMA},
    runner,
};
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps, PrimeCurve, PrimePoint};
use rusqlite::{params, Connection};
use serde_json::{json, Value};
use std::collections::{BTreeSet, VecDeque};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::os::fd::AsRawFd;
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
static STOP: AtomicBool = AtomicBool::new(false);
extern "C" {
    fn flock(fd: i32, operation: i32) -> i32;
    fn signal(sig: i32, handler: usize) -> usize;
    fn kill(pid: i32, sig: i32) -> i32;
}
extern "C" fn stop(_: i32) {
    STOP.store(true, Ordering::Relaxed);
}
fn now() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_millis()
}
fn hash(p: &Path) -> Result<String> {
    Ok(sha256_hex(&fs::read(p)?))
}
fn write_json(p: &Path, v: &Value) -> Result<()> {
    let mut f = File::create(p)?;
    serde_json::to_writer_pretty(&mut f, v)?;
    f.write_all(b"\n")?;
    f.sync_all()?;
    Ok(())
}
fn read_json(p: &Path) -> Result<Value> {
    Ok(serde_json::from_slice(&fs::read(p)?)?)
}
fn bytes(p: &Path) -> Result<u64> {
    let mut n = 0;
    for e in fs::read_dir(p)? {
        let e = e?;
        let t = e.file_type()?;
        if t.is_dir() {
            n += bytes(&e.path())?;
        } else if t.is_file() {
            n += e.metadata()?.len();
        }
    }
    Ok(n)
}
fn db(root: &Path) -> Result<Connection> {
    let c = Connection::open(root.join("ledger.sqlite"))?;
    c.busy_timeout(Duration::from_secs(10))?;
    c.execute_batch("PRAGMA journal_mode=WAL; PRAGMA synchronous=FULL; CREATE TABLE IF NOT EXISTS jobs(id TEXT PRIMARY KEY, ordinal INTEGER UNIQUE NOT NULL, graph TEXT NOT NULL, vertex INTEGER NOT NULL, conductor INTEGER NOT NULL, state TEXT NOT NULL DEFAULT 'queued', attempts INTEGER NOT NULL DEFAULT 0, session TEXT, audit TEXT, summary TEXT); CREATE TABLE IF NOT EXISTS attempts(job TEXT NOT NULL, attempt INTEGER NOT NULL, state TEXT NOT NULL, session TEXT NOT NULL, started_ms TEXT NOT NULL, finished_ms TEXT, PRIMARY KEY(job,attempt));")?;
    Ok(c)
}
fn arm(name: &str, role: &str, id: &str, params: Value) -> Value {
    json!({"name":name,"role":role,"method":{"id":id,"params":params}})
}
fn spec(row: &Value, mode: &str) -> Value {
    let ctl = mode == "control";
    let full = mode == "full";
    let mut arms = vec![
        arm("rho", "reference", "rho.negation", json!({})),
        arm("bsgs", "reference", "bsgs.textbook", json!({})),
        arm("rho-aa", "control", "rho.negation", json!({})),
    ];
    let sizes = if ctl {
        vec![2]
    } else if full {
        vec![16, 32, 64, 128]
    } else {
        vec![16, 64]
    };
    let trials = if ctl {
        1000
    } else if full {
        1_000_000
    } else {
        4096
    };
    let cap = if ctl { 64 } else { 20000 };
    for size in sizes {
        for oracle in ["mitm", "subtract"] {
            arms.push(arm(&format!("projected-{size}-{oracle}"),"candidate","ic.pipeline_counted",json!({"factor_base":format!("prime-projected-abscissa:size={size},scan_cap={cap}"),"oracle":oracle,"targets":"walk","linalg":"incremental-gauss-full-rank-checked","max_trials":trials.to_string()})));
        }
    }
    arms.push(arm("raw-control","baseline","ic.pipeline_counted",json!({"factor_base":format!("prime-abscissa:size={}",if ctl {2}else{64}),"oracle":"mitm","targets":"walk","linalg":"incremental-gauss-full-rank-checked","max_trials":trials.to_string()})));
    json!({"schema":"ecbench.spec/v1","label":format!("volcano-{mode}-{}",row["id"].as_str().unwrap()),"description":"Execution addendum: bounded diagnostic unless full. Public one-target workloads; no runtime promotion.","workloads":{"curves":[row["input"]["curve"].clone()],"targets_per_curve":if full {8}else{2},"target_seed":20261017,"target_kind":"public"},"arms":arms,"measurement":{"rounds":if full {5}else{1},"warmup":if full {1}else{0},"order":"alternate","seed":1,"isolation_required":"L2","timeout_seconds":if ctl {10}else if full {600}else{3}}})
}
fn control() -> Value {
    let c = PrimeCurve {
        p: 113,
        a: 19,
        b: 2,
    };
    let mut ops = GroupOps::default();
    for x in 1..113 {
        if let Some(y) = c.sqrt(c.rhs(x)) {
            let g = c.mul(&mut ops, PrimePoint::affine(x, y), 16);
            if !g.infinity && c.mul(&mut ops, g, 7).infinity {
                return json!({"id":"control-v0","graph":"control","input":{"vertex":0,"conductor":1,"curve":{"kind":"prime_explicit","p":113,"a":19,"b":2,"group_order":112,"r":7,"gx":g.x,"gy":g.y}}});
            }
        }
    }
    panic!("missing fixture");
}
fn plan(root: &Path, input: &Path, binary: &Path, mode: &str) -> Result<()> {
    if !["control", "pilot", "survey", "full"].contains(&mode) {
        return Err("mode must be control, pilot, survey or full".into());
    }
    fs::create_dir(root)?;
    fs::create_dir(root.join("sessions"))?;
    let binary = fs::canonicalize(binary)?;
    fs::copy(&binary, root.join("ecbench"))?;
    let mut q = File::create(root.join("queue.jsonl"))?;
    let mut sources = Vec::new();
    let mut c = db(root)?;
    let tx = c.transaction()?;
    let mut count = 0u64;
    let mut emit = |row: Value| -> Result<()> {
        writeln!(q, "{}", canonical(&row)?)?;
        tx.execute(
            "INSERT INTO jobs(id,ordinal,graph,vertex,conductor) VALUES(?1,?2,?3,?4,?5)",
            params![
                row["id"].as_str().unwrap(),
                count,
                row["graph"].as_str().unwrap(),
                row["input"]["vertex"].as_u64().unwrap(),
                row["input"]["conductor"].as_u64().unwrap()
            ],
        )?;
        count += 1;
        Ok(())
    };
    if mode == "control" {
        emit(control())?;
    } else {
        let mut priority = Vec::new();
        let mut remainder = Vec::new();
        for graph in ["two", "three", "mixed"] {
            let path = fs::canonicalize(input.join(format!("{graph}-workloads.jsonl")))?;
            sources.push(json!({"path":path,"sha256":hash(&path)?}));
            let mut seen = BTreeSet::new();
            let mut first = VecDeque::new();
            let mut rest = VecDeque::new();
            for line in BufReader::new(File::open(&path)?).lines() {
                let mut row: Value = serde_json::from_str(&line?)?;
                row["graph"] = json!(graph);
                row["id"] = json!(format!("{graph}-v{}", row["input"]["vertex"]));
                let conductor = row["input"]["conductor"]
                    .as_u64()
                    .ok_or("missing conductor")?;
                if seen.insert(conductor) {
                    first.push_back(row);
                } else {
                    rest.push_back(row);
                }
                if mode == "pilot" {
                    break;
                }
            }
            priority.push(first);
            remainder.push(rest);
        }
        // Representatives of all 51 strata first, then every other vertex.
        // Each phase preserves round-robin graph order; no vertex is omitted.
        for mut phase in [priority, remainder] {
            loop {
                let mut any = false;
                for queue in &mut phase {
                    if let Some(row) = queue.pop_front() {
                        any = true;
                        emit(row)?;
                    }
                }
                if !any {
                    break;
                }
            }
        }
    }
    drop(emit);
    q.sync_all()?;
    tx.commit()?;
    let source = Command::new("git").args(["rev-parse", "HEAD"]).output()?;
    let diff = Command::new("git").args(["diff", "--binary"]).output()?;
    fs::write(root.join("source.patch"), &diff.stdout)?;
    fs::write(root.join("driver.rs"), include_str!("driver.rs"))?;
    fs::write(
        root.join("execution-protocol.md"),
        include_str!("PROTOCOL.md"),
    )?;
    write_json(
        &root.join("campaign.json"),
        &json!({"schema":"volcano-ic-campaign/v1","created_ms":now().to_string(),"mode":mode,"vertices_planned":count,"queue_sha256":hash(&root.join("queue.jsonl"))?,"ecbench_sha256":hash(&root.join("ecbench"))?,"driver_binary_sha256":hash(&std::env::current_exe()?)?,"driver_source_sha256":hash(&root.join("driver.rs"))?,"protocol_sha256":hash(&root.join("execution-protocol.md"))?,"source_head":String::from_utf8_lossy(&source.stdout).trim(),"source_patch_sha256":hash(&root.join("source.patch"))?,"inputs":sources,"wall_claim":"UNQUALIFIED_DIAGNOSTIC","representation_registration":"SEE_INPUT; NOT_PROMOTED","finite_field_ic":"UNIMPLEMENTED_FOR_ACTUAL_EMBEDDING_DEGREES"}),
    )?;
    status(root)
}
fn command(bin: &Path, args: &[&str], log: &Path) -> Result<i32> {
    let f = File::create(log)?;
    let mut child = Command::new(bin)
        .args(args)
        .env("RAYON_NUM_THREADS", "1")
        .stdout(Stdio::from(f.try_clone()?))
        .stderr(Stdio::from(f))
        .spawn()?;
    loop {
        if let Some(s) = child.try_wait()? {
            return Ok(s.code().unwrap_or(-1));
        }
        if STOP.load(Ordering::Relaxed) {
            unsafe {
                kill(child.id() as i32, 15);
            }
            return Ok(child.wait()?.code().unwrap_or(-1));
        }
        std::thread::sleep(Duration::from_millis(100));
    }
}
fn replay_exhausted(session: &Path, out: &Path) -> Result<bool> {
    if out.exists() {
        return Err("exhaustion replay receipt already exists".into());
    }
    let mut checks = Vec::new();
    for r in runner::read_records(session)? {
        if r.warmup || !r.cost.deterministic || r.outcome.status != "exhausted" {
            continue;
        }
        let input = ChildInput {
            schema: CHILD_INPUT_SCHEMA.into(),
            curve: r.workload.curve_spec.clone(),
            target_seed: r.workload.target_seed,
            target_index: r.workload.target_index,
            target_kind: r.workload.kind(),
            expected_workload_id: r.workload.workload_id.clone(),
            method: MethodSpec {
                id: r.method.id.clone(),
                params: r.method.params.clone(),
            },
            expected_method_id: r.method.method_id.clone(),
            algorithm_seed: r.algorithm_seed.parse()?,
        };
        let output = record::child_main(&serde_json::to_string(&input)?);
        let mut diffs = Vec::new();
        if let Some(rep) = output.report {
            if !rep.exhausted || rep.recovered.map(|k| k.to_string()) != r.outcome.recovered {
                diffs.push("outcome");
            }
            let shape = |phases: &[crypto_lib::cryptanalysis::ecbench::methods::PhaseRecord]| {
                phases
                    .iter()
                    .map(|p| {
                        (
                            p.name.clone(),
                            p.adds,
                            p.doubles,
                            p.scalar_mults,
                            record::json_roundtrip(p.gae).to_bits(),
                            p.native.clone(),
                        )
                    })
                    .collect::<Vec<_>>()
            };
            if shape(&rep.phases) != shape(&r.phases) {
                diffs.push("phases");
            }
            if Some(record::json_roundtrip(rep.total_gae).to_bits())
                != r.cost.total_gae.map(f64::to_bits)
            {
                diffs.push("total_gae");
            }
            if rep.counters != r.counters {
                diffs.push("counters");
            }
            if rep.factor_base != r.factor_base {
                diffs.push("factor_base");
            }
            if rep.field_ops != r.field_ops {
                diffs.push("field_ops");
            }
            if rep.unpriced != r.cost.unpriced
                || rep.deterministic != r.cost.deterministic
                || Some(rep.automorphisms_used) != r.boundaries.automorphisms_used
            {
                diffs.push("accounting_identity");
            }
        } else {
            diffs.push("missing_report");
        }
        checks.push(json!({"record_id":r.record_id,"seq":r.seq,"reproduced":diffs.is_empty(),"differences":diffs,"error":output.error}));
    }
    let ok = checks.iter().all(|v| v["reproduced"] == true);
    write_json(
        out,
        &json!({"schema":"volcano-exhaustion-replay/v1","session":session,"records_sha256":hash(&session.join("records.jsonl"))?,"replayer_binary_sha256":hash(&std::env::current_exe()?)?,"ok":ok,"replays":checks,"scope":"Measured deterministic exhausted records; successful records replayed by ecbench verify; timeouts have no deterministic final report."}),
    )?;
    Ok(ok)
}
fn run(root: &Path, max_jobs: u64, max_seconds: u64, max_bytes: u64) -> Result<()> {
    let lock = OpenOptions::new()
        .create(true)
        .truncate(false)
        .write(true)
        .open(root.join("writer.lock"))?;
    if unsafe { flock(lock.as_raw_fd(), 2 | 4) } != 0 {
        return Err("another runner holds this campaign".into());
    }
    unsafe {
        signal(2, stop as *const () as usize);
        signal(15, stop as *const () as usize);
    }
    let manifest = read_json(&root.join("campaign.json"))?;
    let mode = manifest["mode"].as_str().unwrap();
    if mode == "full" {
        return Err("full launch gated: exact registration and qualified storage/time pilot are still required; use separately labeled survey".into());
    }
    for (field, path) in [
        ("queue_sha256", root.join("queue.jsonl")),
        ("ecbench_sha256", root.join("ecbench")),
        ("driver_binary_sha256", std::env::current_exe()?),
    ] {
        if manifest[field].as_str() != Some(&hash(&path)?) {
            return Err(format!("fingerprint changed: {field}").into());
        }
    }
    let c = db(root)?;
    c.execute(
        "UPDATE jobs SET state='interrupted' WHERE state='running'",
        [],
    )?;
    c.execute(
        "UPDATE attempts SET state='interrupted' WHERE state='running'",
        [],
    )?;
    let start = Instant::now();
    let mut completed = 0;
    let mut used = bytes(root)?;
    write_json(
        &root.join("worker.json"),
        &json!({"pid":std::process::id(),"started_ms":now().to_string(),"max_jobs":max_jobs,"max_seconds":max_seconds,"max_bytes":max_bytes}),
    )?;
    let mut reason = "queue_finished";
    for line in BufReader::new(File::open(root.join("queue.jsonl"))?).lines() {
        let row: Value = serde_json::from_str(&line?)?;
        let id = row["id"].as_str().unwrap();
        let (state, attempts): (String, u64) =
            c.query_row("SELECT state,attempts FROM jobs WHERE id=?1", [id], |r| {
                Ok((r.get(0)?, r.get(1)?))
            })?;
        if state == "audited" {
            continue;
        }
        if STOP.load(Ordering::Relaxed) || root.join("STOP").exists() {
            reason = "stop_requested";
            break;
        }
        if completed >= max_jobs {
            reason = "job_limit";
            break;
        }
        if start.elapsed().as_secs() >= max_seconds {
            reason = "elapsed_limit";
            break;
        }
        if used.saturating_add(64 * 1024 * 1024) >= max_bytes {
            reason = "storage_guard";
            break;
        }
        let attempt = attempts + 1;
        let rel = format!("sessions/{id}-attempt-{attempt:03}");
        let session = root.join(&rel);
        let prefix = root.join(format!("{id}-attempt-{attempt:03}"));
        let specpath = prefix.with_extension("spec.json");
        write_json(&specpath, &spec(&row, mode))?;
        c.execute("INSERT INTO attempts(job,attempt,state,session,started_ms) VALUES(?1,?2,'running',?3,?4)",params![id,attempt,rel,now().to_string()])?;
        c.execute(
            "UPDATE jobs SET state='running',attempts=?2,session=?3 WHERE id=?1",
            params![id, attempt, rel],
        )?;
        println!("{} RUN {id} attempt {attempt}", now());
        std::io::stdout().flush()?;
        let bin = root.join("ecbench");
        let log = prefix.with_extension("run.log");
        let rc = command(
            &bin,
            &[
                "run",
                "--spec",
                specpath.to_str().unwrap(),
                "--out",
                session.to_str().unwrap(),
                "--cpus",
                "auto",
                "--allow-busy",
                "--wait",
                "--settle",
                "0.1",
                "--quiet",
            ],
            &log,
        )?;
        let audit = prefix.with_extension("audit.json");
        let mut ok = false;
        if rc == 0 && !STOP.load(Ordering::Relaxed) {
            let code = command(
                &bin,
                &[
                    "verify",
                    "--dir",
                    session.to_str().unwrap(),
                    "--replay-all",
                    "--out",
                    audit.to_str().unwrap(),
                    "--exit-code",
                ],
                &prefix.with_extension("audit.log"),
            )?;
            ok = code == 0 && read_json(&audit).map(|v| v["ok"] == true).unwrap_or(false);
            if ok {
                ok = replay_exhausted(&session, &prefix.with_extension("exhaustion-replay.json"))?;
            }
        }
        let summary = if session.join("session.json").exists() {
            read_json(&session.join("session.json"))?
        } else {
            json!({"status":"not_started","exit":rc})
        };
        if summary["status_counts"]["wrong_answer"]
            .as_u64()
            .unwrap_or(0)
            > 0
        {
            ok = false;
        }
        let newstate = if ok {
            "audited"
        } else if STOP.load(Ordering::Relaxed) {
            "interrupted"
        } else {
            "failed"
        };
        c.execute(
            "UPDATE jobs SET state=?2,audit=?3,summary=?4 WHERE id=?1",
            params![
                id,
                newstate,
                audit.to_str(),
                serde_json::to_string(&summary)?
            ],
        )?;
        c.execute(
            "UPDATE attempts SET state=?3,finished_ms=?4 WHERE job=?1 AND attempt=?2",
            params![id, attempt, newstate, now().to_string()],
        )?;
        println!("{} {newstate} {id}: {}", now(), summary["status_counts"]);
        std::io::stdout().flush()?;
        if !ok {
            reason = "attempt_failed_or_interrupted";
            break;
        }
        completed += 1;
        used = bytes(root)?;
    }
    write_json(
        &root.join("worker-exit.json"),
        &json!({"pid":std::process::id(),"finished_ms":now().to_string(),"reason":reason,"vertices_this_invocation":completed,"elapsed_seconds":start.elapsed().as_secs(),"spool_bytes":used}),
    )?;
    c.execute_batch("PRAGMA wal_checkpoint(TRUNCATE)")?;
    status(root)?;
    if reason == "attempt_failed_or_interrupted" && !STOP.load(Ordering::Relaxed) {
        return Err("attempt failed; inspect preserved logs before resume".into());
    }
    Ok(())
}
fn status(root: &Path) -> Result<()> {
    let c = db(root)?;
    let mut st = c.prepare(
        "SELECT graph,state,count(*) FROM jobs GROUP BY graph,state ORDER BY graph,state",
    )?;
    let rows=st.query_map([],|r|Ok(json!({"graph":r.get::<_,String>(0)?,"state":r.get::<_,String>(1)?,"vertices":r.get::<_,u64>(2)?})))?.collect::<std::result::Result<Vec<_>,_>>()?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({"campaign":root,"coverage":rows}))?
    );
    Ok(())
}
fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(String::as_str) {
        Some("plan") if a.len()==7 && a[6]=="freeze"=>plan(Path::new(&a[2]),Path::new(&a[3]),Path::new(&a[4]),&a[5]),
        Some("run") if a.len()==6=>run(Path::new(&a[2]),a[3].parse()?,a[4].parse()?,a[5].parse()?),
        Some("status") if a.len()==3=>status(Path::new(&a[2])),
        Some("replay-exhausted") if a.len()==4=>if replay_exhausted(Path::new(&a[2]),Path::new(&a[3]))? {Ok(())}else{Err("exhaustion replay mismatch".into())},
        _=>Err("usage: driver plan OUT INPUTS ECBENCH control|pilot|survey|full freeze | driver run CAMPAIGN MAX_JOBS MAX_SECONDS MAX_BYTES | driver status CAMPAIGN | driver replay-exhausted SESSION FRESH_RECEIPT".into())
    }
}
