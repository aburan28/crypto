//! ephemeral (Rust): run one command on a fresh EC2 instance that terminates
//! itself. A single-binary port of `../ephemeral` (bash) with the same
//! subcommands, flags, S3 layout, SSM config, tags and exit statuses, driving
//! the same `../runner.sh` (embedded at build time). Jobs launched by either
//! CLI can be followed, fetched and killed by the other.

mod cloud;
mod git;
mod job;

use anyhow::{Context, Result, bail};
use aws_sdk_ec2::types::{
    BlockDeviceMapping, EbsBlockDevice, Filter, InstanceInterruptionBehavior,
    InstanceMarketOptionsRequest, InstanceType, LaunchTemplateSpecification, MarketType,
    ResourceType, SpotInstanceType, SpotMarketOptions, Tag, TagSpecification, VolumeType,
};
use base64::Engine;
use clap::{Args, Parser, Subcommand};
use cloud::{Cloud, err};
use job::JobSpec;
use rand::seq::SliceRandom;
use serde_json::Value;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

/// The runner the instance executes. EPHEMERAL_RUNNER overrides it with a file.
const RUNNER: &str = include_str!("../../runner.sh");

#[derive(Parser)]
#[command(
    name = "ephemeral",
    version,
    about = "Run one command on a fresh EC2 instance that terminates itself.",
    after_help = "ID may be `last` for the most recent job launched from this machine.\n\n\
Exit status of `run` and `logs -f`: the command's own status when it ran, 124 on timeout, \
75 when spot reclaimed the instance, 70 when the runner failed before or around the command, \
130 when you detach with Ctrl-C.\n\n\
Environment: EPHEMERAL_CONFIG (SSM parameter, default /ephemeral/config), EPHEMERAL_REGION \
(default AWS_REGION, then us-east-1), AWS_PROFILE, HTTPS_PROXY, AWS_CA_BUNDLE."
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Launch, follow the log, fetch outputs
    Run(Box<RunArgs>),
    /// Live runners and recent jobs
    Ls {
        #[arg(default_value_t = 10)]
        recent: usize,
    },
    /// Phase, instance, exit code
    Status { id: String },
    /// Print (or follow with -f) the job log
    Logs {
        id: String,
        #[arg(short)]
        f: bool,
    },
    /// Download everything the job left in S3
    Fetch { id: String, dir: Option<PathBuf> },
    /// Terminate the runner now
    Kill { id: String },
    /// Session Manager shell on the runner (needs the aws CLI and session-manager-plugin)
    Shell { id: String },
    /// The settings Terraform published
    Config,
    /// Bake the rust toolchain into an AMI (administrator)
    Bake {
        #[arg(long, default_value = "x86_64")]
        arch: String,
        #[arg(long = "type")]
        instance_type: Option<String>,
    },
    /// List the git and cargo caches in S3
    Cache,
}

#[derive(Args, Clone)]
struct RunArgs {
    /// Instance type. GPU types (g*, p*) get the NVIDIA Deep Learning AMI, Graviton types arm64
    #[arg(long = "type", default_value = "c7i.4xlarge")]
    instance_type: String,
    /// Spot capacity; an interruption ends the job as `interrupted`
    #[arg(long)]
    spot: bool,
    /// Hard lifetime, e.g. 90m, 6h, 2d (default from config). The command stops 10 min before it
    #[arg(long)]
    ttl: Option<String>,
    /// Root volume size in GB
    #[arg(long, default_value_t = 100)]
    disk: i32,
    /// none | build | rust
    #[arg(long, default_value = "rust", value_parser = ["none", "build", "rust"])]
    toolchain: String,
    /// Do not restore or save the git object cache or (rust) the cargo cache
    #[arg(long)]
    no_cache: bool,
    /// rust: the directory holding Cargo.lock (default: the repo root)
    #[arg(long, default_value = ".")]
    cargo_dir: String,
    /// Shell command run in the checkout before the job
    #[arg(long, default_value = "")]
    setup: String,
    /// KEY=VALUE for the job (repeatable). Stored in S3: no long-lived secrets
    #[arg(long = "env", value_name = "K=V")]
    envs: Vec<String>,
    /// File or glob, relative to the checkout, copied to outputs (repeatable)
    #[arg(long = "artifact", value_name = "PATH")]
    artifacts: Vec<String>,
    /// Sparse checkout: fetch only these directories (comma-separated)
    #[arg(long, value_delimiter = ',')]
    paths: Vec<String>,
    /// Ship tracked changes only, not new files
    #[arg(long)]
    no_untracked: bool,
    /// Run a pushed commit instead of the working tree
    #[arg(long = "ref", value_name = "REF")]
    git_ref: Option<String>,
    /// Do not check anything out
    #[arg(long)]
    no_source: bool,
    /// Label for the job id (default: the command's first word)
    #[arg(long)]
    name: Option<String>,
    /// Print the job id and return; follow later with `logs -f`
    #[arg(long)]
    detach: bool,
    /// Leave the instance up after the job, until the TTL
    #[arg(long)]
    keep: bool,
    /// Where outputs land (default ~/.cache/ephemeral/runs/ID)
    #[arg(long)]
    out: Option<PathBuf>,
    /// Explicit AMI (default: resolved for the instance type and toolchain)
    #[arg(long)]
    ami: Option<String>,
    /// The command. One argument is a shell snippet; several are an argv
    #[arg(last = true, required = true, value_name = "CMD")]
    command: Vec<String>,
}

#[tokio::main]
async fn main() -> ExitCode {
    // Die quietly on a closed pipe (`ephemeral ls | head`), like a shell tool.
    // SAFETY: restoring the default disposition of one signal at startup.
    unsafe { libc::signal(libc::SIGPIPE, libc::SIG_DFL) };
    let cli = Cli::parse();
    match dispatch(cli.cmd).await {
        Ok(code) => ExitCode::from(code.clamp(0, 255) as u8),
        Err(e) => {
            eprintln!("ephemeral: {e:#}");
            ExitCode::from(1)
        }
    }
}

async fn dispatch(cmd: Cmd) -> Result<i32> {
    let cloud = Cloud::new().await?;
    let ctx = Ctx::load(cloud).await?;
    match cmd {
        Cmd::Run(a) => {
            let (id, out) = ctx.launch(&a).await?;
            if a.detach {
                println!("{id}");
                return Ok(0);
            }
            Ok(ctx
                .follow(
                    &id,
                    Some(&out.unwrap_or_else(|| cache_dir().join("runs").join(&id))),
                )
                .await?
                .code())
        }
        Cmd::Ls { recent } => ctx.ls(recent).await.map(|_| 0),
        Cmd::Status { id } => ctx.status(&resolve_id(&id)?).await.map(|_| 0),
        Cmd::Logs { id, f } => {
            let id = resolve_id(&id)?;
            if f {
                return Ok(ctx
                    .follow(&id, Some(&cache_dir().join("runs").join(&id)))
                    .await?
                    .code());
            }
            let log = ctx
                .cloud
                .get(&ctx.bucket, &format!("jobs/{id}/job.log"), None)
                .await;
            let log = log.with_context(|| format!("no log yet for {id}"))?;
            use std::io::Write;
            std::io::stdout().write_all(&log)?;
            Ok(0)
        }
        Cmd::Fetch { id, dir } => {
            let id = resolve_id(&id)?;
            let dir = dir.unwrap_or_else(|| cache_dir().join("runs").join(&id));
            ctx.cloud
                .sync_down(&ctx.bucket, &format!("jobs/{id}/"), &dir)
                .await?;
            println!("{}", dir.display());
            Ok(0)
        }
        Cmd::Kill { id } => ctx.kill(&resolve_id(&id)?).await.map(|_| 0),
        Cmd::Shell { id } => {
            let id = resolve_id(&id)?;
            let (iid, _) = ctx
                .cloud
                .instance_of(&id)
                .await
                .with_context(|| format!("no instance for {id}"))?;
            let status = std::process::Command::new("aws")
                .args([
                    "ssm",
                    "start-session",
                    "--region",
                    &ctx.cloud.region,
                    "--target",
                    &iid,
                ])
                .status()
                .context("shell needs the aws CLI and session-manager-plugin")?;
            Ok(status.code().unwrap_or(1))
        }
        Cmd::Config => {
            println!("{}", serde_json::to_string_pretty(&ctx.config)?);
            Ok(0)
        }
        Cmd::Bake {
            arch,
            instance_type,
        } => ctx.bake(&arch, instance_type).await.map(|_| 0),
        Cmd::Cache => ctx.cache().await.map(|_| 0),
    }
}

fn cache_dir() -> PathBuf {
    let base = std::env::var_os("XDG_CACHE_HOME")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(std::env::var_os("HOME").unwrap_or_default()).join(".cache")
        });
    base.join("ephemeral")
}

fn resolve_id(id: &str) -> Result<String> {
    let id = if id == "last" {
        std::fs::read_to_string(cache_dir().join("last"))
            .context("no job launched from this machine yet")?
            .trim()
            .to_string()
    } else {
        id.to_string()
    };
    if !job::valid_id(&id) {
        bail!("bad job id: {id}");
    }
    Ok(id)
}

fn now_secs() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs()
}

fn utc_now() -> String {
    chrono::Utc::now().format("%Y-%m-%dT%H:%M:%SZ").to_string()
}

fn utc_of(secs: u64) -> String {
    chrono::DateTime::from_timestamp(secs as i64, 0)
        .map(|t| t.format("%Y-%m-%dT%H:%M:%SZ").to_string())
        .unwrap_or_default()
}

/// `id -un`@`hostname -s`
fn owner() -> String {
    let run = |c: &str, a: &[&str]| {
        std::process::Command::new(c)
            .args(a)
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
    };
    let user = run("id", &["-un"])
        .or_else(|| std::env::var("USER").ok())
        .unwrap_or_else(|| "unknown".into());
    let host = run("hostname", &["-s"]).unwrap_or_else(|| "unknown".into());
    format!("{user}@{host}")
}

fn tag(k: &str, v: &str) -> Tag {
    Tag::builder().key(k).value(v).build()
}

/// How following a job ended.
enum Followed {
    Final { phase: String, code: Option<i64> },
    Detached,
}

impl Followed {
    fn code(&self) -> i32 {
        match self {
            Followed::Final { phase, code } => job::exit_status(phase, *code),
            Followed::Detached => 130,
        }
    }
}

struct Ctx {
    cloud: Cloud,
    config: Value,
    bucket: String,
}

impl Ctx {
    async fn load(cloud: Cloud) -> Result<Ctx> {
        let param =
            std::env::var("EPHEMERAL_CONFIG").unwrap_or_else(|_| "/ephemeral/config".into());
        let raw = cloud.parameter(&param).await.with_context(|| {
            format!(
                "cannot read SSM parameter {param} in {} (deployed? right region/profile?)",
                cloud.region
            )
        })?;
        let config: Value =
            serde_json::from_str(&raw).context("the config parameter is not JSON")?;
        let bucket = config["bucket"]
            .as_str()
            .context("config has no bucket")?
            .to_string();
        Ok(Ctx {
            cloud,
            config,
            bucket,
        })
    }

    fn cfg(&self, key: &str) -> Result<String> {
        match &self.config[key] {
            Value::String(s) => Ok(s.clone()),
            Value::Number(n) => Ok(n.to_string()),
            _ => bail!("config has no {key}"),
        }
    }

    /// SSM parameter naming the baked rust image for an architecture.
    fn baked_param(&self, arch: &str) -> Option<String> {
        let prefix = self.config["baked_ami_prefix"]
            .as_str()
            .filter(|p| !p.is_empty())?;
        Some(format!("{prefix}/rust-{arch}"))
    }

    async fn put(&self, id: &str, name: &str, body: impl Into<Vec<u8>>) -> Result<()> {
        self.cloud
            .put(&self.bucket, &format!("jobs/{id}/{name}"), body.into())
            .await
    }

    /// Stage the job in S3 and launch its instance. Returns the id and --out.
    async fn launch(&self, a: &RunArgs) -> Result<(String, Option<PathBuf>)> {
        if a.no_source && a.git_ref.is_some() {
            bail!("--ref and --no-source contradict each other");
        }
        if a.disk < 1 {
            bail!("--disk wants a number of GB");
        }
        for kv in &a.envs {
            if !kv.contains('=') {
                bail!("--env wants KEY=VALUE");
            }
        }
        let cargo_dir = job::check_cargo_dir(&a.cargo_dir)?;

        let ttl = match &a.ttl {
            Some(t) => t.clone(),
            None => self.cfg("default_ttl")?,
        };
        let ttl_s = job::span_seconds(&ttl)?;
        let max_h: u64 = self.cfg("max_ttl_hours")?.parse()?;
        if ttl_s > max_h * 3600 {
            bail!("--ttl {ttl} is over the {max_h}h ceiling");
        }
        if ttl_s < 900 {
            bail!("--ttl must be at least 15m (boot, setup and upload need time)");
        }

        let name = a
            .name
            .clone()
            .unwrap_or_else(|| job::default_name(&a.command[0]));
        let stamp = chrono::Utc::now().format("%Y%m%d-%H%M%S").to_string();
        let id = job::job_id(&stamp, &name, rand::random());
        let out = match &a.out {
            Some(o) => {
                std::fs::create_dir_all(o)?;
                Some(std::fs::canonicalize(o)?)
            }
            None => None,
        };

        // Source: the commit a runner can fetch, and the patch on top of it.
        let (mut repo_url, mut commit, mut head, mut patch) =
            (String::new(), String::new(), String::new(), None);
        if !a.no_source {
            let repo = git::Repo::discover()?;
            repo_url = match repo.try_git(&["remote", "get-url", "origin"]) {
                Some(u) => job::origin_https(&u),
                None => self.cfg("repo_url")?,
            };
            head = repo.git(&["rev-parse", "HEAD"])?;
            if let Some(r) = &a.git_ref {
                commit = repo.pushed_ref(r)?;
            } else {
                commit = repo.pushed_base()?;
                let tmp = tempdir()?;
                let diff = repo.worktree_patch(&commit, !a.no_untracked, &tmp);
                let _ = std::fs::remove_dir_all(&tmp);
                let diff = diff?;
                if !diff.is_empty() {
                    let n = repo.git(&["rev-list", "--count", &format!("{commit}..HEAD")])?;
                    eprintln!(
                        "ephemeral: shipping {n} unpushed commit(s) + local changes as a {} KiB patch on {}",
                        diff.len() / 1024,
                        &commit[..12]
                    );
                    if diff.len() >= 200 << 20 {
                        bail!(
                            "patch is {} MiB: commit and push large files, or use --ref",
                            diff.len() >> 20
                        );
                    }
                    patch = Some(diff);
                }
            }
        }

        // AMI by architecture and GPU.
        let it = self
            .cloud
            .ec2
            .describe_instance_types()
            .instance_types(InstanceType::from(a.instance_type.as_str()))
            .send()
            .await
            .map_err(|_| {
                anyhow::anyhow!(
                    "unknown instance type {} in {}",
                    a.instance_type,
                    self.cloud.region
                )
            })?;
        let info = it
            .instance_types()
            .first()
            .context("unknown instance type")?;
        let arch = info
            .processor_info()
            .and_then(|p| p.supported_architectures().first())
            .map(|a| a.as_str().to_string())
            .unwrap_or_default();
        let gpu = info
            .gpu_info()
            .and_then(|g| g.gpus().first())
            .and_then(|g| g.manufacturer())
            .unwrap_or("");
        let ami_ref = if let Some(ami) = &a.ami {
            ami.clone()
        } else if gpu == "NVIDIA" {
            if arch != "x86_64" {
                bail!("{}: arm64 GPU instances are not supported", a.instance_type);
            }
            format!("resolve:ssm:{}", self.cfg("gpu_ami_param")?)
        } else if let Some(baked) = self.baked_param(&arch).filter(|_| a.toolchain != "none")
            && self.cloud.parameter(&baked).await.is_ok()
        {
            format!("resolve:ssm:{baked}")
        } else if arch == "arm64" {
            format!("resolve:ssm:{}", self.cfg("arm_ami_param")?)
        } else {
            format!("resolve:ssm:{}", self.cfg("cpu_ami_param")?)
        };

        let deadline = now_secs() + ttl_s;
        let owner = owner();
        let display = a.command.join(" ");
        let spec = JobSpec {
            id: id.clone(),
            created: utc_now(),
            owner: owner.clone(),
            instance_type: a.instance_type.clone(),
            spot: a.spot,
            ttl: ttl.clone(),
            deadline,
            disk_gb: a.disk,
            ami: ami_ref.clone(),
            toolchain: a.toolchain.clone(),
            setup: a.setup.clone(),
            repo_url,
            commit: commit.clone(),
            local_head: head,
            patch: patch.is_some(),
            keep: a.keep,
            command_display: display.clone(),
            paths: a.paths.iter().filter(|p| !p.is_empty()).cloned().collect(),
            artifacts: a
                .artifacts
                .iter()
                .filter(|p| !p.is_empty())
                .cloned()
                .collect(),
            cache: !a.no_cache,
            cargo_dir,
        };
        let runner = match std::env::var("EPHEMERAL_RUNNER") {
            Ok(p) => std::fs::read_to_string(&p)
                .with_context(|| format!("reading EPHEMERAL_RUNNER {p}"))?,
            Err(_) => RUNNER.to_string(),
        };
        let status = format!(
            "{{\"job\":\"{id}\",\"phase\":\"launching\",\"updated\":\"{}\",\"deadline\":{deadline},\"exit_code\":null}}\n",
            utc_now()
        );
        self.put(&id, "job.json", serde_json::to_vec(&spec)?)
            .await?;
        self.put(&id, "cmd.sh", job::cmd_script(&a.command)).await?;
        self.put(&id, "env", job::env_file(&a.envs)).await?;
        self.put(&id, "runner.sh", runner).await?;
        self.put(&id, "status.json", status).await?;
        if let Some(p) = patch {
            self.put(&id, "patch.diff", p).await?;
        }

        let user_data = job::user_data(&id, &self.bucket, &self.cloud.region, deadline, a.keep);
        let managed = self.cfg("managed_tag")?;
        let instance_tags = TagSpecification::builder()
            .resource_type(ResourceType::Instance)
            .tags(tag("Name", &format!("ephemeral-{id}")))
            .tags(tag("ephemeral:job", &id))
            .tags(tag("ephemeral:deadline", &deadline.to_string()))
            .tags(tag(&managed, "true"))
            .tags(tag("ephemeral:owner", &owner))
            .build();
        let volume_tags = TagSpecification::builder()
            .resource_type(ResourceType::Volume)
            .tags(tag(&managed, "true"))
            .tags(tag("ephemeral:job", &id))
            .build();
        let market = a.spot.then(|| {
            InstanceMarketOptionsRequest::builder()
                .market_type(MarketType::Spot)
                .spot_options(
                    SpotMarketOptions::builder()
                        .spot_instance_type(SpotInstanceType::OneTime)
                        .instance_interruption_behavior(InstanceInterruptionBehavior::Terminate)
                        .build(),
                )
                .build()
        });

        std::fs::create_dir_all(cache_dir())?;
        std::fs::write(cache_dir().join("last"), format!("{id}\n"))?;

        // Try each configured subnet (one per AZ) until one has capacity.
        let mut subnets: Vec<String> = self.config["subnet_ids"]
            .as_array()
            .map(|v| {
                v.iter()
                    .filter_map(|s| s.as_str().map(String::from))
                    .collect()
            })
            .unwrap_or_default();
        subnets.shuffle(&mut rand::rng());
        let mut iid = None;
        for subnet in &subnets {
            let res = self
                .cloud
                .ec2
                .run_instances()
                .launch_template(
                    LaunchTemplateSpecification::builder()
                        .launch_template_id(self.cfg("launch_template")?)
                        .version("$Default")
                        .build(),
                )
                .image_id(&ami_ref)
                .instance_type(InstanceType::from(a.instance_type.as_str()))
                .subnet_id(subnet)
                .security_group_ids(self.cfg("security_group")?)
                .block_device_mappings(
                    BlockDeviceMapping::builder()
                        .device_name("/dev/xvda")
                        .ebs(
                            EbsBlockDevice::builder()
                                .volume_size(a.disk)
                                .volume_type(VolumeType::Gp3)
                                .encrypted(true)
                                .delete_on_termination(true)
                                .build(),
                        )
                        .build(),
                )
                .user_data(base64::engine::general_purpose::STANDARD.encode(&user_data))
                .set_instance_market_options(market.clone())
                .tag_specifications(instance_tags.clone())
                .tag_specifications(volume_tags.clone())
                .min_count(1)
                .max_count(1)
                .send()
                .await;
            match res {
                Ok(out) => {
                    iid = out
                        .instances()
                        .first()
                        .and_then(|i| i.instance_id())
                        .map(String::from);
                    break;
                }
                Err(e) => {
                    let msg = err(e).to_string();
                    if [
                        "InsufficientInstanceCapacity",
                        "Unsupported",
                        "capacity-not-available",
                    ]
                    .iter()
                    .any(|s| msg.contains(s))
                    {
                        let short: String = msg.chars().take(200).collect();
                        eprintln!("ephemeral: {subnet}: {short}; trying the next zone");
                        continue;
                    }
                    return Err(self
                        .launch_failed(&id, deadline, &format!("run-instances failed: {msg}"))
                        .await);
                }
            }
        }
        let Some(iid) = iid else {
            let why = format!(
                "no zone could launch {} (reasons above); try another type, --spot or not, or later",
                a.instance_type
            );
            return Err(self.launch_failed(&id, deadline, &why).await);
        };

        eprintln!("ephemeral: job {id}");
        eprintln!(
            "ephemeral:   {}{} {iid} ({ami_ref}), ttl {ttl} (terminates by {})",
            a.instance_type,
            if a.spot { " spot" } else { "" },
            utc_of(deadline)
        );
        let source = if commit.is_empty() {
            String::new()
        } else {
            format!(
                "source {}{}, ",
                &commit[..12],
                if spec.patch { " + local patch" } else { "" }
            )
        };
        eprintln!("ephemeral:   {source}command: {display}");
        eprintln!("ephemeral:   s3://{}/jobs/{id}/", self.bucket);
        Ok((id, out))
    }

    /// A job that never got an instance must not sit at "launching" forever.
    async fn launch_failed(&self, id: &str, deadline: u64, reason: &str) -> anyhow::Error {
        let clean: String = reason
            .chars()
            .filter(|c| !matches!(c, '"' | '\\' | '\n'))
            .take(300)
            .collect();
        let status = serde_json::json!({
            "job": id, "phase": "error", "reason": format!("launch failed: {clean}"),
            "updated": utc_now(), "deadline": deadline, "exit_code": null,
        });
        let _ = self.put(id, "status.json", format!("{status}\n")).await;
        anyhow::anyhow!("{reason}")
    }

    /// Print the log as it grows, return when the job reaches a final phase,
    /// then download its outputs. Ctrl-C detaches; the job keeps running.
    async fn follow(&self, id: &str, out: Option<&Path>) -> Result<Followed> {
        let poll = Duration::from_secs(
            std::env::var("EPHEMERAL_POLL")
                .ok()
                .and_then(|p| p.parse().ok())
                .unwrap_or(10),
        );
        let log_key = format!("jobs/{id}/job.log");
        let mut offset = 0usize;
        let mut tick = 0u64;
        let mut st = None;
        let ctrl_c = tokio::signal::ctrl_c();
        tokio::pin!(ctrl_c);
        let phase = loop {
            offset += self.print_log(&log_key, offset).await;
            st = self.cloud.status(&self.bucket, id).await.or(st);
            let phase = st
                .as_ref()
                .and_then(|s| s["phase"].as_str())
                .unwrap_or("")
                .to_string();
            if job::is_final(&phase) {
                break phase;
            }
            // A runner that died without reporting: check EC2 about once a minute.
            tick += 1;
            if tick.is_multiple_of(6)
                && let Some((_, state)) = self.cloud.instance_of(id).await
                && (state == "terminated" || state == "shutting-down")
            {
                tokio::time::sleep(poll).await;
                st = self.cloud.status(&self.bucket, id).await.or(st);
                let phase = st
                    .as_ref()
                    .and_then(|s| s["phase"].as_str())
                    .unwrap_or("")
                    .to_string();
                break if job::is_final(&phase) && phase != "lost" {
                    phase
                } else {
                    "lost".into()
                };
            }
            tokio::select! {
                _ = tokio::time::sleep(poll) => {}
                _ = &mut ctrl_c => {
                    eprintln!();
                    eprintln!("ephemeral: detached; {id} keeps running. ephemeral logs -f {id} | ephemeral kill {id}");
                    return Ok(Followed::Detached);
                }
            }
        };
        self.print_log(&log_key, offset).await;

        let code = st.as_ref().and_then(|s| s["exit_code"].as_i64());
        if let Some(out) = out {
            if let Err(e) = self
                .cloud
                .sync_down(&self.bucket, &format!("jobs/{id}/"), out)
                .await
            {
                eprintln!("ephemeral: fetching outputs: {e:#}");
            }
            eprintln!("ephemeral: outputs: {}", out.display());
        }
        let shown = code.map(|c| format!(" (exit {c})")).unwrap_or_default();
        eprintln!("ephemeral: job {id}: {phase}{shown}");
        Ok(Followed::Final { phase, code })
    }

    /// Print the log from `offset`; returns the number of bytes printed.
    async fn print_log(&self, key: &str, offset: usize) -> usize {
        use std::io::Write;
        match self
            .cloud
            .get(&self.bucket, key, Some(format!("bytes={offset}-")))
            .await
        {
            Some(chunk) => {
                let mut o = std::io::stdout();
                let _ = o.write_all(&chunk);
                let _ = o.flush();
                chunk.len()
            }
            None => 0,
        }
    }

    async fn ls(&self, recent: usize) -> Result<()> {
        println!("live runners:");
        let states = ["pending", "running", "stopping", "shutting-down"];
        let out = self
            .cloud
            .ec2
            .describe_instances()
            .filters(
                Filter::builder()
                    .name(format!("tag:{}", self.cfg("managed_tag")?))
                    .values("true")
                    .build(),
            )
            .filters(
                Filter::builder()
                    .name("instance-state-name")
                    .set_values(Some(states.map(String::from).to_vec()))
                    .build(),
            )
            .send()
            .await
            .map_err(err)?;
        for i in out.reservations().iter().flat_map(|r| r.instances()) {
            let job = i
                .tags()
                .iter()
                .find(|t| t.key() == Some("ephemeral:job"))
                .and_then(|t| t.value())
                .unwrap_or("None");
            println!(
                "  {job}\t{}\t{}\t{}\t{}\t{}",
                i.instance_id().unwrap_or(""),
                i.instance_type().map(|t| t.as_str()).unwrap_or(""),
                i.state()
                    .and_then(|s| s.name())
                    .map(|n| n.as_str())
                    .unwrap_or(""),
                i.instance_lifecycle()
                    .map(|l| l.as_str())
                    .unwrap_or("on-demand"),
                i.launch_time()
                    .map(|t| utc_of(t.secs() as u64))
                    .unwrap_or_default(),
            );
        }
        println!("recent jobs:");
        let mut ids = Vec::new();
        let mut pages = self
            .cloud
            .s3
            .list_objects_v2()
            .bucket(&self.bucket)
            .prefix("jobs/")
            .delimiter("/")
            .into_paginator()
            .send();
        while let Some(page) = pages.next().await {
            for p in page.map_err(err)?.common_prefixes() {
                if let Some(p) = p.prefix() {
                    ids.push(
                        p.trim_start_matches("jobs/")
                            .trim_end_matches('/')
                            .to_string(),
                    );
                }
            }
        }
        ids.sort();
        for id in ids.iter().skip(ids.len().saturating_sub(recent)) {
            let st = self.cloud.status(&self.bucket, id).await;
            let phase = st.as_ref().and_then(|s| s["phase"].as_str()).unwrap_or("");
            let code = st
                .as_ref()
                .map(|s| &s["exit_code"])
                .filter(|c| !c.is_null())
                .map(|c| c.to_string());
            println!("  {id:<48} {phase} {}", code.unwrap_or_default());
        }
        Ok(())
    }

    async fn status(&self, id: &str) -> Result<()> {
        let st = self
            .cloud
            .status(&self.bucket, id)
            .await
            .with_context(|| format!("no status for {id}"))?;
        println!("{}", serde_json::to_string_pretty(&st)?);
        match self.cloud.instance_of(id).await {
            Some((iid, state)) => println!("instance: {iid}\t{state}"),
            None => println!("instance: "),
        }
        Ok(())
    }

    async fn kill(&self, id: &str) -> Result<()> {
        let (iid, _) = self
            .cloud
            .instance_of(id)
            .await
            .with_context(|| format!("no instance for {id}"))?;
        let out = self
            .cloud
            .ec2
            .terminate_instances()
            .instance_ids(&iid)
            .send()
            .await
            .map_err(err)?;
        let state = out
            .terminating_instances()
            .first()
            .and_then(|t| t.current_state())
            .and_then(|s| s.name())
            .map(|n| n.as_str().to_string())
            .unwrap_or_default();
        println!("{state}");
        let status = serde_json::json!({
            "job": id, "phase": "killed", "reason": format!("ephemeral kill from {}", owner()),
            "instance_id": iid, "updated": utc_now(), "exit_code": null,
        });
        self.put(id, "status.json", format!("{status}\n")).await?;
        eprintln!("ephemeral: terminated {iid} ({id})");
        Ok(())
    }

    async fn cache(&self) -> Result<()> {
        let mut any = false;
        let mut pages = self
            .cloud
            .s3
            .list_objects_v2()
            .bucket(&self.bucket)
            .prefix("cache/")
            .into_paginator()
            .send();
        while let Some(page) = pages.next().await {
            for o in page.map_err(err)?.contents() {
                any = true;
                let when = o
                    .last_modified()
                    .map(|t| utc_of(t.secs() as u64))
                    .unwrap_or_default();
                println!(
                    "{when} {:>12} {}",
                    job::human_size(o.size().unwrap_or(0)),
                    o.key().unwrap_or("")
                );
            }
        }
        if !any {
            println!("no caches yet");
        }
        Ok(())
    }

    /// Bake the rust preset into an AMI: run a --keep job on the stock image,
    /// snapshot it once the toolchain is installed, and point the SSM parameter
    /// at it. Needs ec2:CreateImage, DeregisterImage, DeleteSnapshot and
    /// ssm:PutParameter on top of the launcher policy.
    async fn bake(&self, arch: &str, instance_type: Option<String>) -> Result<()> {
        let (default_type, stock) = match arch {
            "x86_64" => ("c7i.xlarge", self.cfg("cpu_ami_param")?),
            "arm64" => ("c7g.xlarge", self.cfg("arm_ami_param")?),
            _ => bail!("--arch is x86_64 or arm64"),
        };
        let param = self
            .baked_param(arch)
            .context("this deployment publishes no baked_ami_prefix: terraform apply first")?;
        let script = r#"set -e; . "$HOME/.cargo/env"; rustc -V; cargo -V
         printf "toolchain=rust\nbaked=%s\nrustc=%s\n" "$(date -u +%FT%TZ)" "$(rustc -V)" >/etc/ephemeral-image
         dnf clean all -q; sync"#;
        let args = RunArgs::try_parse_from_args(&[
            "--type",
            &instance_type.unwrap_or_else(|| default_type.into()),
            "--toolchain",
            "rust",
            "--no-source",
            "--keep",
            "--detach",
            "--ttl",
            "1h",
            "--disk",
            "20",
            "--ami",
            &format!("resolve:ssm:{stock}"),
            "--name",
            &format!("bake-{arch}"),
            "--",
            script,
        ])?;
        let (id, _) = self
            .launch(&args)
            .await
            .context("bake job did not launch")?;
        // Whatever happens next, the bake instance does not outlive this command.
        let res = self.bake_image(&id, arch, &param).await;
        if let Some((iid, state)) = self.cloud.instance_of(&id).await
            && state != "terminated"
        {
            let _ = self
                .cloud
                .ec2
                .terminate_instances()
                .instance_ids(&iid)
                .send()
                .await;
        }
        let ami = res?;
        println!("{ami}");
        Ok(())
    }

    async fn bake_image(&self, id: &str, arch: &str, param: &str) -> Result<String> {
        match self.follow(id, None).await? {
            Followed::Final { phase, .. } if phase == "done" => {}
            _ => bail!("bake job {id} did not finish cleanly"),
        }
        let (iid, _) = self
            .cloud
            .instance_of(id)
            .await
            .with_context(|| format!("no instance for {id}"))?;
        let name = format!(
            "ephemeral-rust-{arch}-{}",
            chrono::Utc::now().format("%Y%m%d-%H%M%S")
        );
        eprintln!("ephemeral: creating image {name} from {iid}");
        let baked = format!("rust-{arch}");
        let ami = self
            .cloud
            .ec2
            .create_image()
            .instance_id(&iid)
            .name(&name)
            .no_reboot(true)
            .description(format!(
                "ephemeral runner: AL2023 + build tools + rustup stable (job {id})"
            ))
            .tag_specifications(
                TagSpecification::builder()
                    .resource_type(ResourceType::Image)
                    .tags(tag("ephemeral:baked", &baked))
                    .tags(tag("ephemeral:job", id))
                    .build(),
            )
            .tag_specifications(
                TagSpecification::builder()
                    .resource_type(ResourceType::Snapshot)
                    .tags(tag("ephemeral:baked", &baked))
                    .build(),
            )
            .send()
            .await
            .map_err(err)
            .context("create-image failed")?
            .image_id()
            .context("create-image returned no id")?
            .to_string();
        let mut state = String::from("pending");
        for _ in 0..120 {
            tokio::time::sleep(Duration::from_secs(15)).await;
            let out = self
                .cloud
                .ec2
                .describe_images()
                .image_ids(&ami)
                .send()
                .await
                .map_err(err)?;
            state = out
                .images()
                .first()
                .and_then(|i| i.state())
                .map(|s| s.as_str().to_string())
                .unwrap_or_default();
            if state != "pending" {
                break;
            }
        }
        if state != "available" {
            bail!("image {ami} is {state}");
        }
        self.cloud
            .ec2
            .terminate_instances()
            .instance_ids(&iid)
            .send()
            .await
            .map_err(err)?;
        self.cloud
            .ssm
            .put_parameter()
            .name(param)
            .r#type(aws_sdk_ssm::types::ParameterType::String)
            .data_type("aws:ec2:image")
            .value(&ami)
            .overwrite(true)
            .send()
            .await
            .map_err(err)?;
        eprintln!("ephemeral: {param} -> {ami}");

        // Keep the newest two bakes; jobs already launching hold their own copy.
        let out = self
            .cloud
            .ec2
            .describe_images()
            .owners("self")
            .filters(
                Filter::builder()
                    .name("tag:ephemeral:baked")
                    .values(&baked)
                    .build(),
            )
            .send()
            .await
            .map_err(err)?;
        let mut images: Vec<_> = out.images().to_vec();
        images.sort_by(|a, b| a.creation_date().cmp(&b.creation_date()));
        let n = images.len().saturating_sub(2);
        for old in &images[..n] {
            let Some(old_id) = old.image_id() else {
                continue;
            };
            self.cloud
                .ec2
                .deregister_image()
                .image_id(old_id)
                .send()
                .await
                .map_err(err)?;
            for snap in old
                .block_device_mappings()
                .iter()
                .filter_map(|m| m.ebs().and_then(|e| e.snapshot_id()))
            {
                self.cloud
                    .ec2
                    .delete_snapshot()
                    .snapshot_id(snap)
                    .send()
                    .await
                    .map_err(err)?;
            }
            eprintln!("ephemeral: deregistered {old_id}");
        }
        Ok(ami)
    }
}

impl RunArgs {
    fn try_parse_from_args(args: &[&str]) -> Result<RunArgs> {
        #[derive(Parser)]
        struct Wrap {
            #[command(flatten)]
            run: RunArgs,
        }
        let argv = std::iter::once("ephemeral-run").chain(args.iter().copied());
        Ok(Wrap::try_parse_from(argv)?.run)
    }
}

fn tempdir() -> Result<PathBuf> {
    let dir = std::env::temp_dir().join(format!(
        "ephemeral-{}-{:x}",
        std::process::id(),
        rand::random::<u32>()
    ));
    std::fs::create_dir(&dir)?;
    Ok(dir)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_args_need_the_separator() {
        let a = RunArgs::try_parse_from_args(&[
            "--spot",
            "--paths",
            "taskq,src",
            "--env",
            "A=1",
            "--",
            "echo",
            "hi",
        ])
        .unwrap();
        assert!(a.spot);
        assert_eq!(a.paths, ["taskq", "src"]);
        assert_eq!(a.envs, ["A=1"]);
        assert_eq!(a.command, ["echo", "hi"]);
        assert_eq!(a.instance_type, "c7i.4xlarge");
        assert!(RunArgs::try_parse_from_args(&["echo", "hi"]).is_err());
        assert!(RunArgs::try_parse_from_args(&["--toolchain", "go", "--", "x"]).is_err());
    }

    #[test]
    fn runner_is_embedded() {
        assert!(RUNNER.starts_with("#!"));
        assert!(RUNNER.contains("status.json"));
    }
}
