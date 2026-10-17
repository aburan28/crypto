//! Compare completed deterministic accounting across frozen compiler profiles.
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use rusqlite::Connection;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
type Key = (String, String, u64, String, bool);
fn records(root: &Path) -> Result<BTreeMap<Key, Value>> {
    let conn = Connection::open_with_flags(
        root.join("ledger.sqlite"),
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
    )?;
    let mut st = conn.prepare("SELECT state,session FROM jobs ORDER BY ordinal")?;
    let mut out = BTreeMap::new();
    for job in st.query_map([], |r| {
        Ok((r.get::<_, String>(0)?, r.get::<_, Option<String>>(1)?))
    })? {
        let (state, session) = job?;
        if state != "audited" {
            return Err("all planned control/pilot jobs must be audited before matching".into());
        }
        for line in BufReader::new(File::open(
            root.join(session.ok_or("missing session")?)
                .join("records.jsonl"),
        )?)
        .lines()
        {
            let r: Value = serde_json::from_str(&line?)?;
            let key = (
                r["workload"]["workload_id"]
                    .as_str()
                    .ok_or("workload")?
                    .to_string(),
                r["arm"].as_str().ok_or("arm")?.to_string(),
                r["round"].as_u64().ok_or("round")?,
                r["algorithm_seed"].as_str().ok_or("seed")?.to_string(),
                r["warmup"].as_bool().ok_or("warmup")?,
            );
            if out.insert(key, r).is_some() {
                return Err("duplicate match key".into());
            }
        }
    }
    Ok(out)
}
fn counts(r: &Value) -> Value {
    let phases: Vec<_> = r["phases"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| {
            let mut v = p.clone();
            v.as_object_mut().unwrap().remove("wall_ns");
            v
        })
        .collect();
    json!({"status":r["outcome"]["status"],"recovered":r["outcome"]["recovered"],"matches_target":r["outcome"]["matches_target"],"cost":r["cost"],"phases":phases,"counters":r["counters"],"factor_base":r["factor_base"],"field_ops":r["field_ops"],"automorphisms_used":r["boundaries"]["automorphisms_used"]})
}
fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 4 {
        return Err("usage: match_profiles BASE_CAMPAIGN CANDIDATE_CAMPAIGN FRESH_RECEIPT".into());
    }
    let out = Path::new(&a[3]);
    if out.exists() {
        return Err("receipt already exists".into());
    }
    let base = records(Path::new(&a[1]))?;
    let candidate = records(Path::new(&a[2]))?;
    let bk: BTreeSet<_> = base.keys().collect();
    let ck: BTreeSet<_> = candidate.keys().collect();
    if bk != ck {
        return Err("different workload/arm/round/seed populations".into());
    }
    let mut checks = Vec::new();
    let mut compared = 0;
    let mut censored = 0;
    let mut ok = true;
    for (key, b) in &base {
        let c = &candidate[key];
        let mut diff = Vec::new();
        if b["method"]["method_id"] != c["method"]["method_id"]
            || b["method"]["params"] != c["method"]["params"]
        {
            diff.push("method_identity");
        }
        let final_report = |r: &Value| {
            r["cost"]["deterministic"] == true
                && r["cost"]["total_gae"].is_number()
                && ["verified", "exhausted"]
                    .iter()
                    .any(|s| r["outcome"]["status"] == *s)
        };
        let classification = if final_report(b) && final_report(c) {
            compared += 1;
            if counts(b) != counts(c) {
                diff.push("deterministic_accounting");
            }
            "completed_count_comparison"
        } else {
            censored += 1;
            if b["outcome"]["status"] != "timeout" && c["outcome"]["status"] != "timeout" {
                diff.push("unexpected_nonreport_outcome");
            }
            "censored_not_count_comparable"
        };
        if b["outcome"]["status"] == "wrong_answer" || c["outcome"]["status"] == "wrong_answer" {
            diff.push("wrong_answer");
        }
        ok &= diff.is_empty();
        checks.push(json!({"workload":key.0,"arm":key.1,"round":key.2,"seed":key.3,"warmup":key.4,"baseline_record":b["record_id"],"candidate_record":c["record_id"],"baseline_status":b["outcome"]["status"],"candidate_status":c["outcome"]["status"],"classification":classification,"differences":diff}));
    }
    ok &= compared > 0;
    let hash = |p: &Path| -> Result<String> { Ok(sha256_hex(&fs::read(p)?)) };
    let receipt = json!({"schema":"volcano-profile-accounting-match/v1","ok":ok,"baseline_campaign":a[1],"candidate_campaign":a[2],"baseline_manifest_sha256":hash(&Path::new(&a[1]).join("campaign.json"))?,"candidate_manifest_sha256":hash(&Path::new(&a[2]).join("campaign.json"))?,"completed_pairs_compared":compared,"censored_pairs_not_compared":censored,"scope":"Exact deterministic outcomes, GAE, phase group/native counts, field operations, counters, base identity and accounting flags. Elapsed values are excluded. No equality claim for censored pairs.","checks":checks});
    serde_json::to_writer_pretty(File::create(out)?, &receipt)?;
    println!(
        "{}",
        json!({"ok":ok,"compared":compared,"censored":censored})
    );
    if !ok {
        return Err("profile-accounting mismatch; campaign launch must halt".into());
    }
    Ok(())
}
