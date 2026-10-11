//! `spotlab mcp`: a Model Context Protocol server over stdio
//! (newline-delimited JSON-RPC 2.0). It also runs the controller's reconcile
//! loop in the background while it is up, so a Claude Code session that
//! submits work keeps resuming it across preemptions without anything else
//! running.

use anyhow::{bail, Result};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Duration;

use crate::controller::Controller;
use crate::ops;
use crate::spec::{Input, JobSpec, Placement, Resources};

pub fn serve(c: Controller, reconcile_every_s: u64) -> Result<()> {
    let c = Arc::new(c);
    let lock = Arc::new(Mutex::new(()));
    if reconcile_every_s > 0 {
        let (c2, l2) = (c.clone(), lock.clone());
        std::thread::spawn(move || loop {
            {
                let _g = l2.lock().unwrap();
                match c2.reconcile(reconcile_every_s as f64 * 3.0) {
                    Ok(ev) => ev.iter().for_each(|e| eprintln!("spotlab reconcile: {e}")),
                    Err(e) => eprintln!("spotlab reconcile: {e:#}"),
                }
            }
            std::thread::sleep(Duration::from_secs(reconcile_every_s));
        });
    }
    let stdin = std::io::stdin();
    let mut stdout = std::io::stdout();
    for line in stdin.lock().lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let msg: Value = match serde_json::from_str(&line) {
            Ok(v) => v,
            Err(e) => {
                write(
                    &mut stdout,
                    &json!({"jsonrpc":"2.0","id":null,"error":{"code":-32700,"message":e.to_string()}}),
                )?;
                continue;
            }
        };
        let Some(id) = msg.get("id").cloned() else {
            continue;
        }; // notification
        let method = msg["method"].as_str().unwrap_or("");
        let reply = match method {
            "initialize" => json!({"jsonrpc":"2.0","id":id,"result":{
                "protocolVersion": msg["params"]["protocolVersion"].as_str().unwrap_or("2025-06-18"),
                "capabilities": {"tools": {}},
                "serverInfo": {"name": "spotlab", "version": env!("CARGO_PKG_VERSION")},
                "instructions": INSTRUCTIONS,
            }}),
            "ping" => json!({"jsonrpc":"2.0","id":id,"result":{}}),
            "tools/list" => json!({"jsonrpc":"2.0","id":id,"result":{"tools": tools()}}),
            "tools/call" => {
                let name = msg["params"]["name"].as_str().unwrap_or("");
                let args = msg["params"]["arguments"].clone();
                let out = call(&c, &lock, name, &args);
                let (text, is_error) = match out {
                    Ok(v) => (serde_json::to_string_pretty(&v)?, false),
                    Err(e) => (format!("{e:#}"), true),
                };
                json!({"jsonrpc":"2.0","id":id,"result":{"content":[{"type":"text","text":text}],"isError":is_error}})
            }
            _ => {
                json!({"jsonrpc":"2.0","id":id,"error":{"code":-32601,"message":format!("unknown method {method}")}})
            }
        };
        write(&mut stdout, &reply)?;
    }
    Ok(())
}

fn write(out: &mut impl Write, v: &Value) -> Result<()> {
    writeln!(out, "{}", serde_json::to_string(v)?)?;
    out.flush()?;
    Ok(())
}

const INSTRUCTIONS: &str = "spotlab runs resumable jobs on the cheapest spot capacity across AWS and GCP. \
Call spotlab_overview first, spotlab_offers to see prices, spotlab_submit to start a job. A job must keep \
its resumable state in $SPOTLAB_CHECKPOINT_DIR (atomic writes, save on SIGTERM, restore at start) and its \
results in $SPOTLAB_OUTPUT_DIR. Results are throughput runs, never performance measurements.";

fn tools() -> Value {
    let obj =
        |props: Value, req: &[&str]| json!({"type":"object","properties":props,"required":req});
    let resources = json!({
        "vcpus": {"type":"integer","description":"minimum vCPUs"},
        "memory_gb": {"type":"number"},
        "arch": {"type":"string","enum":["x86_64","aarch64","any"]},
        "cpu_flags": {"type":"array","items":{"type":"string"},"description":"e.g. avx512, vpclmulqdq"},
        "providers": {"type":"array","items":{"type":"string","enum":["aws","gcp","local"]}},
        "max_hourly_usd": {"type":"number"},
    });
    let mut submit_props = resources.as_object().unwrap().clone();
    for (k, v) in json!({
        "spec": {"type":"object","description":"a full spotlab.job/v1 document; when given, the flat fields are ignored"},
        "name": {"type":"string"},
        "argv": {"type":"array","items":{"type":"string"}},
        "env": {"type":"object","additionalProperties":{"type":"string"}},
        "inputs": {"type":"array","items":{"type":"object","properties":{
            "path":{"type":"string"},"content":{"type":"string"},"local_file":{"type":"string"},"executable":{"type":"boolean"}}},
            "description":"files placed in the working directory; local_file is read on this machine and uploaded"},
        "checkpoint_interval_s": {"type":"integer"},
        "max_total_hours": {"type":"number"},
        "max_cost_usd": {"type":"number"},
        "labels": {"type":"object","additionalProperties":{"type":"string"}},
        "launch_now": {"type":"boolean","description":"run a reconcile pass immediately (default true)"},
    }).as_object().unwrap() {
        submit_props.insert(k.clone(), v.clone());
    }
    json!([
        {"name":"spotlab_overview","description":"Store, providers, jobs by state, running instances and the current burn rate. Call first.","inputSchema":obj(json!({}),&[])},
        {"name":"spotlab_offers","description":"Cheapest spot offers that fit the given resources, across providers (AWS prices are live).","inputSchema":obj(resources.clone(),&[])},
        {"name":"spotlab_submit","description":"Submit a resumable job and (by default) launch it on the cheapest offer.","inputSchema":obj(Value::Object(submit_props),&[])},
        {"name":"spotlab_status","description":"A job's record, latest heartbeat and checkpoints.","inputSchema":obj(json!({"job_id":{"type":"string"}}),&["job_id"])},
        {"name":"spotlab_jobs","description":"List jobs, optionally by state (queued, running, succeeded, failed, cancelled, dead).","inputSchema":obj(json!({"state":{"type":"string"}}),&[])},
        {"name":"spotlab_result","description":"The result document of a finished job: attempts, preemptions, cost, outputs.","inputSchema":obj(json!({"job_id":{"type":"string"}}),&["job_id"])},
        {"name":"spotlab_logs","description":"Tail of an attempt's stdout, stderr or agent log (uploaded when the attempt ends).","inputSchema":obj(json!({"job_id":{"type":"string"},"attempt":{"type":"integer"},"stream":{"type":"string","enum":["stdout","stderr","agent"]},"max_bytes":{"type":"integer"}}),&["job_id"])},
        {"name":"spotlab_fetch","description":"Download a finished job's outputs to a local directory, verifying hashes.","inputSchema":obj(json!({"job_id":{"type":"string"},"dest":{"type":"string"}}),&["job_id","dest"])},
        {"name":"spotlab_cancel","description":"Cancel a job; its instance is terminated on the next reconcile.","inputSchema":obj(json!({"job_id":{"type":"string"},"reconcile_now":{"type":"boolean"}}),&["job_id"])},
        {"name":"spotlab_reconcile","description":"Run one controller pass now: settle ended attempts, relaunch preempted jobs, launch queued ones.","inputSchema":obj(json!({}),&[])},
    ])
}

fn resources_from(a: &Value) -> (Resources, Placement) {
    let mut r = Resources::default();
    if let Some(v) = a["vcpus"].as_u64() {
        r.vcpus = v as u32;
    }
    if let Some(v) = a["memory_gb"].as_f64() {
        r.memory_gb = v;
    }
    if let Some(v) = a["arch"].as_str() {
        r.arch = v.into();
    }
    if let Some(v) = a["cpu_flags"].as_array() {
        r.cpu_flags = v
            .iter()
            .filter_map(|x| x.as_str().map(String::from))
            .collect();
    }
    let mut p = Placement::default();
    if let Some(v) = a["providers"].as_array() {
        p.providers = v
            .iter()
            .filter_map(|x| x.as_str().map(String::from))
            .collect();
    }
    p.max_hourly_usd = a["max_hourly_usd"].as_f64();
    (r, p)
}

fn str_map(v: &Value) -> BTreeMap<String, String> {
    v.as_object()
        .map(|m| {
            m.iter()
                .filter_map(|(k, v)| v.as_str().map(|s| (k.clone(), s.to_string())))
                .collect()
        })
        .unwrap_or_default()
}

fn call(c: &Controller, lock: &Mutex<()>, name: &str, a: &Value) -> Result<Value> {
    let job = || {
        a["job_id"]
            .as_str()
            .map(String::from)
            .ok_or_else(|| anyhow::anyhow!("job_id is required"))
    };
    match name {
        "spotlab_overview" => ops::overview(c),
        "spotlab_offers" => {
            let (r, p) = resources_from(a);
            Ok(ops::offers(c, &r, &p, 20))
        }
        "spotlab_submit" => {
            let spec = if a["spec"].is_object() {
                JobSpec::parse(&a["spec"].to_string())?
            } else {
                let argv: Vec<String> = a["argv"]
                    .as_array()
                    .map(|v| {
                        v.iter()
                            .filter_map(|x| x.as_str().map(String::from))
                            .collect()
                    })
                    .unwrap_or_default();
                let mut spec = JobSpec::parse(&json!({"name": a["name"].as_str().unwrap_or("job"), "command": {"argv": argv}}).to_string())?;
                spec.command.env = str_map(&a["env"]);
                let (r, p) = resources_from(a);
                spec.resources = r;
                spec.placement = p;
                spec.labels = str_map(&a["labels"]);
                if let Some(v) = a["checkpoint_interval_s"].as_u64() {
                    spec.checkpoint.interval_s = v;
                }
                if let Some(v) = a["max_total_hours"].as_f64() {
                    spec.limits.max_total_hours = v;
                }
                spec.limits.max_cost_usd = a["max_cost_usd"].as_f64();
                for i in a["inputs"].as_array().cloned().unwrap_or_default() {
                    spec.inputs.push(Input {
                        path: i["path"].as_str().unwrap_or("").into(),
                        content: i["content"].as_str().map(String::from),
                        store_key: None,
                        local_file: i["local_file"].as_str().map(String::from),
                        sha256: None,
                        executable: i["executable"].as_bool().unwrap_or(false),
                    });
                }
                spec
            };
            let job_id = c.submit(spec, &std::env::current_dir()?)?;
            let events = if a["launch_now"].as_bool().unwrap_or(true) {
                let _g = lock.lock().unwrap();
                c.reconcile(180.0)?
            } else {
                vec![]
            };
            Ok(json!({"job_id": job_id, "events": events, "status": ops::status(c, &job_id)?}))
        }
        "spotlab_status" => ops::status(c, &job()?),
        "spotlab_jobs" => ops::jobs(c, a["state"].as_str()),
        "spotlab_result" => ops::result(c, &job()?),
        "spotlab_logs" => ops::logs(
            c,
            &job()?,
            a["attempt"].as_u64().map(|x| x as u32),
            a["stream"].as_str().unwrap_or("stdout"),
            a["max_bytes"].as_u64().unwrap_or(8192) as usize,
        ),
        "spotlab_fetch" => {
            let Some(d) = a["dest"].as_str() else {
                bail!("dest is required")
            };
            ops::fetch(c, &job()?, &PathBuf::from(d))
        }
        "spotlab_cancel" => {
            let id = job()?;
            c.cancel(&id)?;
            let events = if a["reconcile_now"].as_bool().unwrap_or(true) {
                let _g = lock.lock().unwrap();
                c.reconcile(180.0)?
            } else {
                vec![]
            };
            Ok(json!({"job_id": id, "cancel_requested": true, "events": events}))
        }
        "spotlab_reconcile" => {
            let _g = lock.lock().unwrap();
            Ok(json!({"events": c.reconcile(180.0)?}))
        }
        other => bail!("unknown tool {other}"),
    }
}
