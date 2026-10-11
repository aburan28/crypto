//! Read-side operations shared by the CLI and the MCP server. Each returns JSON.

use anyhow::{bail, Context, Result};
use serde_json::{json, Value};
use std::path::Path;

use crate::controller::Controller;
use crate::record::*;
use crate::spec::{Placement, Resources};
use crate::util::{now, sha256_file};

pub fn overview(c: &Controller) -> Result<Value> {
    let t = now();
    let mut counts = std::collections::BTreeMap::<String, usize>::new();
    let mut running = Vec::new();
    let mut hourly = 0.0;
    let mut spent = 0.0;
    for id in c.job_ids()? {
        let Ok(r) = c.record(&id) else { continue };
        *counts.entry(r.state.as_str().into()).or_default() += 1;
        spent += r.cost_usd(t);
        if r.state == JobState::Running {
            if let Some(a) = r.current() {
                hourly += a.hourly_usd;
                running.push(json!({
                    "job_id": r.job_id, "name": r.name, "attempt": r.attempt,
                    "provider": a.instance.provider, "instance_type": a.instance.instance_type,
                    "zone": a.instance.zone, "hourly_usd": a.hourly_usd,
                    "running_for_h": (t - a.launched_at) / 3600.0,
                }));
            }
        }
    }
    Ok(json!({
        "store": c.store.url(),
        "providers": c.providers.iter().map(|p| p.name()).collect::<Vec<_>>(),
        "jobs_by_state": counts,
        "running": running,
        "burn_usd_per_hour": hourly,
        "estimated_spend_usd_all_jobs": spent,
        "budget": {"max_running": c.cfg.budget.max_running, "max_hourly_usd": c.cfg.budget.max_hourly_usd},
    }))
}

pub fn offers(c: &Controller, r: &Resources, p: &Placement, limit: usize) -> Value {
    let (offers, errors) = c.offers(r, p);
    json!({"offers": offers.into_iter().take(limit).collect::<Vec<_>>(), "errors": errors})
}

pub fn status(c: &Controller, id: &str) -> Result<Value> {
    let rec = c.record(id)?;
    let hb: Option<Heartbeat> = if rec.attempt > 0 {
        c.store
            .get_json(&attempt_key(id, rec.attempt, "heartbeat.json"))?
    } else {
        None
    };
    let seqs = checkpoint_seqs(&c.store, id)?;
    Ok(json!({
        "record": rec,
        "heartbeat": hb,
        "checkpoints": seqs,
        "cost_usd_estimate": rec.cost_usd(now()),
        "has_result": c.store.exists(&job_key(id, "result.json"))?,
    }))
}

pub fn jobs(c: &Controller, state: Option<&str>) -> Result<Value> {
    let t = now();
    let mut out = Vec::new();
    for id in c.job_ids()? {
        let Ok(r) = c.record(&id) else { continue };
        if state.is_some_and(|s| s != r.state.as_str()) {
            continue;
        }
        out.push(json!({
            "job_id": r.job_id, "name": r.name, "state": r.state.as_str(), "attempt": r.attempt,
            "reason": r.reason, "cost_usd_estimate": r.cost_usd(t), "created_at": r.created_at,
        }));
    }
    Ok(Value::Array(out))
}

pub fn result(c: &Controller, id: &str) -> Result<Value> {
    c.store
        .get_json::<Value>(&job_key(id, "result.json"))?
        .with_context(|| format!("job {id} has no result yet; see status"))
}

pub fn logs(
    c: &Controller,
    id: &str,
    attempt: Option<u32>,
    stream: &str,
    max_bytes: usize,
) -> Result<Value> {
    let rec = c.record(id)?;
    let n = attempt.unwrap_or(rec.attempt);
    if n == 0 {
        bail!("job {id} has not started");
    }
    let name = match stream {
        "stdout" | "stderr" | "agent" => format!("{stream}.log"),
        other => bail!("stream must be stdout, stderr or agent, not {other}"),
    };
    let text = match c.store.get(&attempt_key(id, n, &name))? {
        Some(b) => {
            let start = b.len().saturating_sub(max_bytes);
            String::from_utf8_lossy(&b[start..]).into_owned()
        }
        None => String::new(),
    };
    Ok(
        json!({"job_id": id, "attempt": n, "stream": stream, "text": text,
        "note": "logs are uploaded when an attempt ends"}),
    )
}

/// Download the job's outputs into `dest`, verifying hashes against the result.
pub fn fetch(c: &Controller, id: &str, dest: &Path) -> Result<Value> {
    let res: ResultDoc = c
        .store
        .get_json(&job_key(id, "result.json"))?
        .with_context(|| format!("job {id} has no result yet"))?;
    let mut files = Vec::new();
    for f in &res.outputs {
        crate::util::safe_rel(&f.path)?;
        let p = dest.join(&f.path);
        if !c
            .store
            .get_file(&job_key(id, &format!("output/{}", f.path)), &p)?
        {
            bail!("output {} missing from the store", f.path);
        }
        let got = sha256_file(&p)?;
        if got != f.sha256 {
            bail!(
                "output {}: sha256 {got} does not match the result's {}",
                f.path,
                f.sha256
            );
        }
        files.push(p.display().to_string());
    }
    Ok(json!({"job_id": id, "dest": dest.display().to_string(), "files": files, "verified": true}))
}
