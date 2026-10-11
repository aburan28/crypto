//! spotlab: resumable experiments on the cheapest AWS or GCP spot capacity.
//! See README.md.

mod agent;
mod controller;
mod demo;
mod mcp;
mod ops;
mod provider;
mod record;
mod spec;
mod store;
mod util;

use anyhow::{bail, Context, Result};
use serde_json::Value;
use std::path::{Path, PathBuf};

use controller::Controller;
use spec::{Config, Input, JobSpec, Placement, Resources};

const USAGE: &str = "\
spotlab — resumable experiments on the cheapest AWS/GCP spot capacity

Controller (reads spotlab.json, $SPOTLAB_CONFIG or --config PATH):
  spotlab offers [--vcpus N] [--memory-gb G] [--arch A] [--flags a,b] [--provider P]... [--max-hourly X]
  spotlab submit SPEC.json [--no-launch]
  spotlab run [--name N] [--vcpus N] [--memory-gb G] [--arch A] [--flags a,b] [--provider P]...
              [--max-hourly X] [--interval S] [--max-hours H] [--max-cost X]
              [--input LOCAL[:DEST]]... [--env K=V]... [--no-launch] -- CMD ARGS...
  spotlab reconcile [--watch SECONDS]
  spotlab overview | jobs [--state S] | status JOB | result JOB
  spotlab logs JOB [--attempt N] [--stream stdout|stderr|agent]
  spotlab fetch JOB --dest DIR | cancel JOB
  spotlab publish-agent [--binary PATH] [--arch x86_64|aarch64]
  spotlab mcp [--reconcile-every SECONDS]      MCP server over stdio
  spotlab local-preempt JOB                    simulate an interruption (local provider)

On the instance:
  spotlab agent --store URL --job ID --attempt N --provider aws|gcp|local [--workdir DIR] [--preempt-file F]

Reference workload:
  spotlab demo-job --to N [--every K] [--step-us U]
";

struct Args(Vec<String>);

impl Args {
    fn take(&mut self, flag: &str) -> Option<String> {
        let i = self.0.iter().position(|a| a == flag)?;
        if i + 1 >= self.0.len() {
            return None;
        }
        let v = self.0.remove(i + 1);
        self.0.remove(i);
        Some(v)
    }
    fn take_all(&mut self, flag: &str) -> Vec<String> {
        let mut out = Vec::new();
        while let Some(v) = self.take(flag) {
            out.push(v);
        }
        out
    }
    fn flag(&mut self, flag: &str) -> bool {
        match self.0.iter().position(|a| a == flag) {
            Some(i) => {
                self.0.remove(i);
                true
            }
            None => false,
        }
    }
    fn pos(&mut self, what: &str) -> Result<String> {
        if self.0.is_empty() || self.0[0].starts_with("--") {
            bail!("missing {what}\n\n{USAGE}");
        }
        Ok(self.0.remove(0))
    }
    fn done(&self) -> Result<()> {
        if !self.0.is_empty() {
            bail!("unexpected arguments: {}", self.0.join(" "));
        }
        Ok(())
    }
}

fn print(v: &Value) {
    println!("{}", serde_json::to_string_pretty(v).unwrap_or_default());
}

fn resources(a: &mut Args) -> Result<(Resources, Placement)> {
    let mut r = Resources::default();
    if let Some(v) = a.take("--vcpus") {
        r.vcpus = v.parse().context("--vcpus")?;
    }
    if let Some(v) = a.take("--memory-gb") {
        r.memory_gb = v.parse().context("--memory-gb")?;
    }
    if let Some(v) = a.take("--arch") {
        r.arch = v;
    }
    if let Some(v) = a.take("--flags") {
        r.cpu_flags = v
            .split(',')
            .filter(|s| !s.is_empty())
            .map(String::from)
            .collect();
    }
    let mut p = Placement {
        providers: a.take_all("--provider"),
        ..Default::default()
    };
    if let Some(v) = a.take("--max-hourly") {
        p.max_hourly_usd = Some(v.parse().context("--max-hourly")?);
    }
    Ok((r, p))
}

fn main() {
    if let Err(e) = real_main() {
        eprintln!("spotlab: {e:#}");
        std::process::exit(1);
    }
}

fn real_main() -> Result<()> {
    let mut a = Args(std::env::args().skip(1).collect());
    let config = a.take("--config");
    let cmd = if a.0.is_empty() {
        "help".to_string()
    } else {
        a.0.remove(0)
    };
    let ctl = |config: &Option<String>| -> Result<Controller> {
        Controller::new(Config::locate(config.as_deref())?)
    };
    match cmd.as_str() {
        "help" | "-h" | "--help" => print!("{USAGE}"),
        "agent" => {
            let args = agent::AgentArgs {
                store: a.take("--store").context("--store")?,
                job: a.take("--job").context("--job")?,
                attempt: a.take("--attempt").context("--attempt")?.parse()?,
                provider: a.take("--provider").unwrap_or_else(|| "local".into()),
                workdir: a.take("--workdir").map(PathBuf::from),
                preempt_file: a.take("--preempt-file").map(PathBuf::from),
            };
            a.done()?;
            agent::run_agent(args)?;
        }
        "demo-job" => {
            let to: u64 = a.take("--to").context("--to")?.parse()?;
            let every: u64 = a
                .take("--every")
                .map(|v| v.parse())
                .transpose()?
                .unwrap_or(1000);
            let step_us: u64 = a
                .take("--step-us")
                .map(|v| v.parse())
                .transpose()?
                .unwrap_or(0);
            a.done()?;
            demo::run(to, every, step_us)?;
        }
        "offers" => {
            let (r, p) = resources(&mut a)?;
            a.done()?;
            print(&ops::offers(&ctl(&config)?, &r, &p, 30));
        }
        "submit" => {
            let path = a.pos("SPEC.json")?;
            let no_launch = a.flag("--no-launch");
            a.done()?;
            let text = std::fs::read_to_string(&path).with_context(|| format!("read {path}"))?;
            let spec = JobSpec::parse(&text)?;
            let base = Path::new(&path)
                .parent()
                .map(Path::to_path_buf)
                .unwrap_or_default();
            let c = ctl(&config)?;
            submit_and_report(&c, spec, &base, no_launch)?;
        }
        "run" => {
            let sep =
                a.0.iter()
                    .position(|x| x == "--")
                    .context("missing -- before the command")?;
            let argv: Vec<String> = a.0.split_off(sep)[1..].to_vec();
            let name = a
                .take("--name")
                .unwrap_or_else(|| argv.first().cloned().unwrap_or_default());
            let (r, p) = resources(&mut a)?;
            let interval = a.take("--interval");
            let max_hours = a.take("--max-hours");
            let max_cost = a.take("--max-cost");
            let inputs = a.take_all("--input");
            let envs = a.take_all("--env");
            let no_launch = a.flag("--no-launch");
            a.done()?;
            let mut spec = JobSpec::parse(
                &serde_json::json!({"name": name, "command": {"argv": argv}}).to_string(),
            )?;
            spec.resources = r;
            spec.placement = p;
            if let Some(v) = interval {
                spec.checkpoint.interval_s = v.parse()?;
            }
            if let Some(v) = max_hours {
                spec.limits.max_total_hours = v.parse()?;
            }
            if let Some(v) = max_cost {
                spec.limits.max_cost_usd = Some(v.parse()?);
            }
            for kv in envs {
                let (k, v) = kv.split_once('=').context("--env K=V")?;
                spec.command.env.insert(k.into(), v.into());
            }
            for i in inputs {
                let (local, dest) = match i.split_once(':') {
                    Some((l, d)) => (l.to_string(), d.to_string()),
                    None => (
                        i.clone(),
                        Path::new(&i)
                            .file_name()
                            .context("--input")?
                            .to_string_lossy()
                            .into_owned(),
                    ),
                };
                spec.inputs.push(Input {
                    path: dest,
                    content: None,
                    store_key: None,
                    local_file: Some(local),
                    sha256: None,
                    executable: false,
                });
            }
            spec.validate()?;
            let c = ctl(&config)?;
            submit_and_report(&c, spec, &std::env::current_dir()?, no_launch)?;
        }
        "reconcile" => {
            let watch: Option<u64> = a.take("--watch").map(|v| v.parse()).transpose()?;
            a.done()?;
            let c = ctl(&config)?;
            loop {
                let ttl = watch.map(|w| w as f64 * 3.0).unwrap_or(120.0);
                for e in c.reconcile(ttl)? {
                    println!("{e}");
                }
                match watch {
                    Some(w) => std::thread::sleep(std::time::Duration::from_secs(w.max(5))),
                    None => {
                        c.release_lease()?;
                        break;
                    }
                }
            }
        }
        "overview" => print(&ops::overview(&ctl(&config)?)?),
        "jobs" => {
            let st = a.take("--state");
            a.done()?;
            print(&ops::jobs(&ctl(&config)?, st.as_deref())?);
        }
        "status" => {
            let id = a.pos("JOB")?;
            print(&ops::status(&ctl(&config)?, &id)?);
        }
        "result" => {
            let id = a.pos("JOB")?;
            print(&ops::result(&ctl(&config)?, &id)?);
        }
        "logs" => {
            let id = a.pos("JOB")?;
            let attempt = a.take("--attempt").map(|v| v.parse()).transpose()?;
            let stream = a.take("--stream").unwrap_or_else(|| "stdout".into());
            a.done()?;
            let v = ops::logs(&ctl(&config)?, &id, attempt, &stream, 1 << 20)?;
            print!("{}", v["text"].as_str().unwrap_or(""));
        }
        "fetch" => {
            let id = a.pos("JOB")?;
            let dest = a.take("--dest").context("--dest")?;
            a.done()?;
            print(&ops::fetch(&ctl(&config)?, &id, Path::new(&dest))?);
        }
        "cancel" => {
            let id = a.pos("JOB")?;
            a.done()?;
            ctl(&config)?.cancel(&id)?;
            println!("cancel requested for {id}; the next reconcile terminates its instance");
        }
        "publish-agent" => {
            let bin = a
                .take("--binary")
                .map(PathBuf::from)
                .unwrap_or(std::env::current_exe()?);
            let arch = a
                .take("--arch")
                .unwrap_or_else(|| std::env::consts::ARCH.to_string());
            a.done()?;
            if !matches!(arch.as_str(), "x86_64" | "aarch64") {
                bail!("--arch must be x86_64 or aarch64");
            }
            let c = ctl(&config)?;
            let key = format!("bin/spotlab-{arch}");
            c.store.put_file(&key, &bin)?;
            println!(
                "uploaded {} as {}/{key} (sha256 {})",
                bin.display(),
                c.store.url(),
                util::sha256_file(&bin)?
            );
        }
        "mcp" => {
            let every: u64 = a
                .take("--reconcile-every")
                .map(|v| v.parse())
                .transpose()?
                .unwrap_or(60);
            a.done()?;
            mcp::serve(ctl(&config)?, every)?;
        }
        "local-preempt" => {
            let id = a.pos("JOB")?;
            a.done()?;
            let c = ctl(&config)?;
            let lc = c
                .cfg
                .local
                .clone()
                .context("no local provider in the config")?;
            let rec = c
                .store
                .get_json::<record::JobRecord>(&record::job_key(&id, "record.json"))?
                .context("no such job")?;
            let f = provider::local::Local { cfg: lc }.preempt_file(&id, rec.attempt);
            std::fs::write(&f, b"preempt\n")?;
            println!("wrote {}", f.display());
        }
        other => bail!("unknown command {other}\n\n{USAGE}"),
    }
    Ok(())
}

fn submit_and_report(c: &Controller, spec: JobSpec, base: &Path, no_launch: bool) -> Result<()> {
    let id = c.submit(spec, base)?;
    println!("submitted {id}");
    if !no_launch {
        for e in c.reconcile(120.0)? {
            println!("{e}");
        }
        c.release_lease()?;
    }
    Ok(())
}
