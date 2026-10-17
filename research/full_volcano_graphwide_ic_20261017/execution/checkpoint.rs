//! Read-only receipt/seal inventory, then a byte-preserving stopped checkpoint.
//! This checks receipt bindings; it is not an independent solver replay.
use crypto_lib::cryptanalysis::ecbench::{canonical::sha256_hex, record::check_line_seal};
use rusqlite::Connection;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader};
use std::os::fd::AsRawFd;
use std::path::{Component, Path};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
extern "C" {
    fn flock(fd: i32, operation: i32) -> i32;
}
fn read(p: &Path) -> Result<Value> {
    Ok(serde_json::from_slice(&fs::read(p)?)?)
}
fn hash(p: &Path) -> Result<String> {
    Ok(sha256_hex(&fs::read(p)?))
}
fn write(p: &Path, v: &Value) -> Result<()> {
    let f = File::create(p)?;
    serde_json::to_writer_pretty(&f, v)?;
    f.sync_all()?;
    Ok(())
}
fn safe(s: &str) -> Result<&Path> {
    let p = Path::new(s);
    if p.components().all(|c| matches!(c, Component::Normal(_))) {
        Ok(p)
    } else {
        Err("non-relative capsule path".into())
    }
}
fn require(ok: bool, what: &str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(what.into())
    }
}
fn copy(src: &Path, dst: &Path) -> Result<()> {
    fs::create_dir(dst)?;
    for e in fs::read_dir(src)? {
        let e = e?;
        let n = e.file_name();
        let name = n.to_string_lossy();
        if name == "STOP"
            || name == "writer.lock"
            || name.ends_with("-wal")
            || name.ends_with("-shm")
        {
            continue;
        }
        if e.file_type()?.is_dir() {
            copy(&e.path(), &dst.join(n))?;
        } else {
            fs::copy(e.path(), dst.join(n))?;
        }
    }
    Ok(())
}
fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: checkpoint STOPPED_CAMPAIGN FROZEN_DRIVER FRESH_OUT".into());
    }
    let root = Path::new(&args[1]);
    let driver = Path::new(&args[2]);
    let out = Path::new(&args[3]);
    let lock = OpenOptions::new()
        .write(true)
        .open(root.join("writer.lock"))?;
    require(
        unsafe { flock(lock.as_raw_fd(), 2 | 4) } == 0,
        "campaign still has a writer; pause before snapshot",
    )?;
    if let Ok(m) = fs::metadata(root.join("ledger.sqlite-wal")) {
        require(
            m.len() == 0,
            "nonempty WAL: checkpoint must finish before snapshot",
        )?;
    }
    let manifest = read(&root.join("campaign.json"))?;
    for (key, path) in [
        ("queue_sha256", root.join("queue.jsonl")),
        ("ecbench_sha256", root.join("ecbench")),
        ("driver_binary_sha256", driver.to_path_buf()),
    ] {
        require(
            manifest[key].as_str() == Some(&hash(&path)?),
            &format!("fingerprint mismatch: {key}"),
        )?;
    }
    let mut queue = BTreeMap::new();
    let mut planned: BTreeMap<String, BTreeMap<u64, u64>> = BTreeMap::new();
    for line in BufReader::new(File::open(root.join("queue.jsonl"))?).lines() {
        let row: Value = serde_json::from_str(&line?)?;
        let id = row["id"].as_str().ok_or("queue id")?.to_string();
        let graph = row["graph"].as_str().ok_or("graph")?.to_string();
        let f = row["input"]["conductor"].as_u64().ok_or("conductor")?;
        *planned
            .entry(graph.clone())
            .or_default()
            .entry(f)
            .or_default() += 1;
        let old = queue.insert(
            id,
            (
                graph,
                row["input"]["vertex"].as_u64().ok_or("vertex")?,
                f,
                row["input"]["curve"]["r"].as_u64().ok_or("r")?,
            ),
        );
        require(old.is_none(), "duplicate queue id")?;
    }
    require(
        queue.len() as u64
            == manifest["vertices_planned"]
                .as_u64()
                .ok_or("planned size")?,
        "queue population mismatch",
    )?;
    let conn = Connection::open_with_flags(
        root.join("ledger.sqlite"),
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
    )?;
    let mut jobs =
        conn.prepare("SELECT id,graph,vertex,conductor,state,session FROM jobs ORDER BY ordinal")?;
    let mut seen = BTreeSet::new();
    let mut done: BTreeMap<String, BTreeMap<u64, u64>> = BTreeMap::new();
    let mut states: BTreeMap<String, u64> = BTreeMap::new();
    let mut outcomes: BTreeMap<String, u64> = BTreeMap::new();
    let mut capsules = Vec::new();
    let mut measured = 0;
    let mut replayed_ok = 0;
    let mut replayed_exhausted = 0;
    for job in jobs.query_map([], |r| {
        Ok((
            r.get::<_, String>(0)?,
            r.get::<_, String>(1)?,
            r.get::<_, u64>(2)?,
            r.get::<_, u64>(3)?,
            r.get::<_, String>(4)?,
            r.get::<_, Option<String>>(5)?,
        ))
    })? {
        let (id, graph, vertex, f, state, session) = job?;
        let expected = queue.get(&id).ok_or("ledger job missing from queue")?;
        require(
            expected.0 == graph && expected.1 == vertex && expected.2 == f,
            "ledger input identity mismatch",
        )?;
        require(seen.insert(id.clone()), "duplicate ledger job")?;
        *states.entry(state.clone()).or_default() += 1;
        require(state != "running", "unfinalized job")?;
        if state != "audited" {
            continue;
        }
        let rel = session.ok_or("audited job without session")?;
        let path = root.join(safe(&rel)?);
        let name = path
            .file_name()
            .ok_or("session name")?
            .to_str()
            .ok_or("session UTF8")?;
        let receipt_path = root.join(format!("{name}.audit.json"));
        let receipt = read(&receipt_path)?;
        require(receipt["ok"] == true, "failed standard audit")?;
        require(
            receipt["auditor_binary_sha256"] == manifest["ecbench_sha256"],
            "auditor binary mismatch",
        )?;
        for (file, expected_hash) in receipt["files"].as_object().ok_or("audit file inventory")? {
            require(
                expected_hash.as_str() == Some(&hash(&path.join(safe(file)?))?),
                &format!("audit file hash mismatch: {id}/{file}"),
            )?;
        }
        let s = read(&path.join("session.json"))?;
        require(
            s["status"] == "complete" && s["session_id"] == receipt["session_id"],
            "session status/identity",
        )?;
        let exhaustion_path = root.join(format!("{name}.exhaustion-replay.json"));
        let exhaustion = read(&exhaustion_path)?;
        require(exhaustion["ok"] == true, "failed exhaustion replay")?;
        require(
            exhaustion["records_sha256"].as_str() == Some(&hash(&path.join("records.jsonl"))?),
            "exhaustion receipt binding",
        )?;
        require(
            exhaustion["replayer_binary_sha256"] == manifest["driver_binary_sha256"],
            "exhaustion replayer identity",
        )?;
        let mut expected_success = BTreeSet::new();
        let mut expected_exhausted = BTreeSet::new();
        let mut lines = 0;
        for line in BufReader::new(File::open(path.join("records.jsonl"))?).lines() {
            let line = line?;
            let r: Value = serde_json::from_str(&line)?;
            let rid = r["record_id"].as_str().ok_or("record id")?;
            require(check_line_seal(&line, rid), "record seal mismatch")?;
            require(
                r["session_id"] == s["session_id"],
                "record session mismatch",
            )?;
            require(
                r["workload"]["curve"]["r"].as_u64() == Some(expected.3),
                "subgroup mismatch",
            )?;
            lines += 1;
            if r["warmup"] == true {
                continue;
            }
            measured += 1;
            let outcome = r["outcome"]["status"].as_str().ok_or("outcome")?;
            *outcomes.entry(outcome.to_string()).or_default() += 1;
            require(outcome != "wrong_answer", "wrong answer in audited session")?;
            if r["cost"]["deterministic"] == true {
                if outcome == "verified" {
                    require(
                        r["outcome"]["matches_target"] == true,
                        "missing verified multiplication",
                    )?;
                    expected_success.insert(rid.to_string());
                } else if outcome == "exhausted" {
                    expected_exhausted.insert(rid.to_string());
                }
            }
        }
        require(
            s["records_written"].as_u64() == Some(lines)
                && receipt["records"].as_u64() == Some(lines),
            "record population mismatch",
        )?;
        let collect = |v: &Value| -> Result<BTreeSet<String>> {
            let mut set = BTreeSet::new();
            for r in v.as_array().ok_or("replay inventory")? {
                require(r["reproduced"] == true, "replay mismatch")?;
                require(
                    set.insert(
                        r["record_id"]
                            .as_str()
                            .ok_or("replay record id")?
                            .to_string(),
                    ),
                    "duplicate replay id",
                )?;
            }
            Ok(set)
        };
        let success = collect(&receipt["replays"])?;
        let exhausted = collect(&exhaustion["replays"])?;
        require(
            success == expected_success && exhausted == expected_exhausted,
            "missing or extra deterministic replay",
        )?;
        replayed_ok += success.len();
        replayed_exhausted += exhausted.len();
        *done.entry(graph.clone()).or_default().entry(f).or_default() += 1;
        capsules.push(json!({"id":id,"graph":graph,"vertex":vertex,"conductor":f,"session":rel,"standard_audit_sha256":hash(&receipt_path)?,"exhaustion_receipt_sha256":hash(&exhaustion_path)?,"records_sha256":hash(&path.join("records.jsonl"))?,"successful_replays":success.len(),"exhausted_replays":exhausted.len()}));
    }
    require(
        seen.len() == queue.len(),
        "ledger population differs from frozen queue",
    )?;
    let mut graphs = Vec::new();
    let mut planned_strata = 0;
    let mut audited_strata = 0;
    for (graph, strata) in planned {
        let audited = done.get(&graph).cloned().unwrap_or_default();
        let expected_set: BTreeSet<_> = strata.keys().copied().collect();
        let actual_set: BTreeSet<_> = audited.keys().copied().collect();
        require(
            actual_set.is_subset(&expected_set),
            "unexpected audited conductor",
        )?;
        planned_strata += expected_set.len();
        audited_strata += actual_set.len();
        let strata_rows:Vec<_>=strata.iter().map(|(f,n)|json!({"conductor":f,"vertices_planned":n,"vertices_audited":audited.get(f).copied().unwrap_or(0)})).collect();
        graphs.push(json!({"graph":graph,"vertices_planned":strata.values().sum::<u64>(),"vertices_audited":audited.values().sum::<u64>(),"strata_planned":expected_set.len(),"strata_audited":actual_set.len(),"max_planned_conductor":expected_set.last(),"max_audited_conductor":actual_set.last(),"strata":strata_rows}));
    }
    let summary = json!({"schema":"volcano-checkpoint-inventory/v1","ok":true,"scope":"Receipt/file hashes, record seals, ledger/queue identities and exact replay coverage; not a fresh or independent solver replay.","campaign_manifest_sha256":hash(&root.join("campaign.json"))?,"vertices_planned":queue.len(),"vertices_audited":capsules.len(),"states":states,"strata_planned":planned_strata,"strata_audited":audited_strata,"measured_records":measured,"outcome_totals":outcomes,"successful_replays":replayed_ok,"exhausted_replays":replayed_exhausted,"graphs":graphs,"capsules":capsules});
    fs::create_dir(out)?;
    write(&out.join("inventory.json"), &summary)?;
    copy(root, &out.join("campaign"))?;
    fs::create_dir(out.join("bin"))?;
    fs::copy(driver, out.join("bin/volcano-ic-driver"))?;
    fs::write(
        out.join("checkpoint-source.rs"),
        include_str!("checkpoint.rs"),
    )?;
    println!(
        "{}",
        json!({"ok":true,"vertices":summary["vertices_audited"],"strata":summary["strata_audited"],"measured_records":measured,"outcomes":summary["outcome_totals"],"out":out})
    );
    Ok(())
}
