//! The job specification (`spotlab.job/v1`) and the controller config.

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;

use crate::util::{safe_rel, sha256_hex};

pub const JOB_SCHEMA: &str = "spotlab.job/v1";
pub const RESULT_SCHEMA: &str = "spotlab.result/v1";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct JobSpec {
    #[serde(default = "job_schema")]
    pub schema: String,
    pub name: String,
    pub command: Command,
    #[serde(default)]
    pub inputs: Vec<Input>,
    #[serde(default)]
    pub resources: Resources,
    #[serde(default)]
    pub placement: Placement,
    #[serde(default)]
    pub checkpoint: CheckpointPolicy,
    #[serde(default)]
    pub limits: Limits,
    #[serde(default)]
    pub labels: BTreeMap<String, String>,
}

fn job_schema() -> String {
    JOB_SCHEMA.to_string()
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Command {
    pub argv: Vec<String>,
    #[serde(default)]
    pub env: BTreeMap<String, String>,
    /// Relative to the job's working directory.
    #[serde(default)]
    pub cwd: Option<String>,
}

/// A file placed in the working directory before the command starts.
/// Exactly one of `content`, `store_key` or `local_file` is set; the CLI and
/// the MCP server upload a `local_file` to the store and rewrite it as a
/// `store_key` with its `sha256` before the spec is saved.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub content: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub store_key: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub local_file: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
    #[serde(default)]
    pub executable: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Resources {
    #[serde(default = "one")]
    pub vcpus: u32,
    #[serde(default = "one_f")]
    pub memory_gb: f64,
    /// `x86_64`, `aarch64`, or `any`.
    #[serde(default = "any")]
    pub arch: String,
    /// CPU features every candidate instance type must list (`avx512`, `vpclmulqdq`, ...).
    #[serde(default)]
    pub cpu_flags: Vec<String>,
}

impl Default for Resources {
    fn default() -> Self {
        Resources {
            vcpus: 1,
            memory_gb: 1.0,
            arch: "any".into(),
            cpu_flags: vec![],
        }
    }
}

fn one() -> u32 {
    1
}
fn one_f() -> f64 {
    1.0
}
fn any() -> String {
    "any".into()
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Default)]
#[serde(deny_unknown_fields)]
pub struct Placement {
    /// Providers this job may run on; empty means every configured provider.
    #[serde(default)]
    pub providers: Vec<String>,
    /// Never launch an offer above this hourly price (USD).
    #[serde(default)]
    pub max_hourly_usd: Option<f64>,
    /// Only these instance types, if set.
    #[serde(default)]
    pub instance_types: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct CheckpointPolicy {
    /// How often the agent uploads the checkpoint directory when it changed.
    #[serde(default = "default_interval")]
    pub interval_s: u64,
    /// Seconds the command gets after SIGTERM on a preemption notice before
    /// SIGKILL. GCP gives 30 s of notice in total and the final upload must
    /// fit too, so the default is conservative.
    #[serde(default = "default_grace")]
    pub grace_s: u64,
    /// Checkpoints kept in the store; older ones are deleted.
    #[serde(default = "default_keep")]
    pub keep: usize,
}

impl Default for CheckpointPolicy {
    fn default() -> Self {
        CheckpointPolicy {
            interval_s: default_interval(),
            grace_s: default_grace(),
            keep: default_keep(),
        }
    }
}

fn default_interval() -> u64 {
    300
}
fn default_grace() -> u64 {
    15
}
fn default_keep() -> usize {
    2
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Limits {
    /// Launches (first run plus resumes) before the job is declared dead.
    #[serde(default = "default_attempts")]
    pub max_attempts: u32,
    /// Wall-clock budget summed over all attempts.
    #[serde(default = "default_hours")]
    pub max_total_hours: f64,
    /// Spend ceiling for this job, estimated from launch prices.
    #[serde(default)]
    pub max_cost_usd: Option<f64>,
    /// An attempt whose agent has not written a heartbeat for this long is
    /// treated as lost and its instance terminated.
    #[serde(default = "default_stale")]
    pub heartbeat_stale_s: u64,
}

impl Default for Limits {
    fn default() -> Self {
        Limits {
            max_attempts: default_attempts(),
            max_total_hours: default_hours(),
            max_cost_usd: None,
            heartbeat_stale_s: default_stale(),
        }
    }
}

fn default_attempts() -> u32 {
    20
}
fn default_hours() -> f64 {
    24.0
}
fn default_stale() -> u64 {
    900
}

impl JobSpec {
    pub fn parse(text: &str) -> Result<JobSpec> {
        let spec: JobSpec = serde_json::from_str(text).context("parse job spec")?;
        spec.validate()?;
        Ok(spec)
    }

    pub fn validate(&self) -> Result<()> {
        if self.schema != JOB_SCHEMA {
            bail!("schema must be {JOB_SCHEMA}, got {}", self.schema);
        }
        if self.name.trim().is_empty() {
            bail!("name is required");
        }
        if self.command.argv.is_empty() {
            bail!("command.argv is empty");
        }
        if let Some(cwd) = &self.command.cwd {
            safe_rel(cwd)?;
        }
        for k in self.command.env.keys() {
            if k.starts_with("SPOTLAB_") {
                bail!("env {k}: the SPOTLAB_ prefix is reserved");
            }
        }
        for i in &self.inputs {
            safe_rel(&i.path)?;
            let n = i.content.is_some() as u8
                + i.store_key.is_some() as u8
                + i.local_file.is_some() as u8;
            if n != 1 {
                bail!(
                    "input {}: set exactly one of content, store_key, local_file",
                    i.path
                );
            }
        }
        if !matches!(self.resources.arch.as_str(), "x86_64" | "aarch64" | "any") {
            bail!("resources.arch must be x86_64, aarch64 or any");
        }
        if self.resources.vcpus == 0 {
            bail!("resources.vcpus must be at least 1");
        }
        for p in &self.placement.providers {
            if !matches!(p.as_str(), "aws" | "gcp" | "local") {
                bail!("unknown provider {p}");
            }
        }
        if self.checkpoint.interval_s < 10 {
            bail!("checkpoint.interval_s must be at least 10");
        }
        if self.limits.max_attempts == 0 {
            bail!("limits.max_attempts must be at least 1");
        }
        Ok(())
    }

    /// The canonical form (every default written out) and its hash.
    pub fn canonical(&self) -> Result<(String, String)> {
        let v = serde_json::to_value(self)?;
        let text = serde_json::to_string(&v)?; // serde_json maps are sorted (BTreeMap)
        let hash = sha256_hex(text.as_bytes());
        Ok((text, hash))
    }
}

/// Controller configuration: where state lives and what may be launched.
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Config {
    /// `s3://bucket/prefix`, `gs://bucket/prefix`, or a local directory.
    pub store: String,
    #[serde(default)]
    pub aws: Option<crate::provider::aws::AwsConfig>,
    #[serde(default)]
    pub gcp: Option<crate::provider::gcp::GcpConfig>,
    #[serde(default)]
    pub local: Option<crate::provider::local::LocalConfig>,
    #[serde(default)]
    pub budget: Budget,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Budget {
    /// Instances running at once across every job.
    #[serde(default = "default_running")]
    pub max_running: usize,
    /// Sum of hourly prices of everything running.
    #[serde(default = "default_hourly")]
    pub max_hourly_usd: f64,
}

impl Default for Budget {
    fn default() -> Self {
        Budget {
            max_running: default_running(),
            max_hourly_usd: default_hourly(),
        }
    }
}

fn default_running() -> usize {
    4
}
fn default_hourly() -> f64 {
    5.0
}

impl Config {
    pub fn load(path: &Path) -> Result<Config> {
        let text =
            std::fs::read_to_string(path).with_context(|| format!("read {}", path.display()))?;
        let cfg: Config =
            serde_json::from_str(&text).with_context(|| format!("parse {}", path.display()))?;
        if cfg.store.trim().is_empty() {
            bail!("config: store is required");
        }
        Ok(cfg)
    }

    /// Resolve the config: `--config`, then `SPOTLAB_CONFIG`, then `./spotlab.json`.
    pub fn locate(explicit: Option<&str>) -> Result<Config> {
        let path = explicit
            .map(String::from)
            .or_else(|| std::env::var("SPOTLAB_CONFIG").ok())
            .unwrap_or_else(|| "spotlab.json".into());
        Config::load(Path::new(&path))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn minimal() -> &'static str {
        r#"{"name":"t","command":{"argv":["true"]}}"#
    }

    #[test]
    fn defaults_fill_and_hash_is_stable() {
        let a = JobSpec::parse(minimal()).unwrap();
        let (text, h1) = a.canonical().unwrap();
        let b = JobSpec::parse(&text).unwrap();
        let (_, h2) = b.canonical().unwrap();
        assert_eq!(h1, h2);
        assert_eq!(a.checkpoint.grace_s, 15);
        assert_eq!(a.limits.max_attempts, 20);
    }

    #[test]
    fn rejects_bad_specs() {
        assert!(JobSpec::parse(r#"{"name":"t","command":{"argv":[]}}"#).is_err());
        assert!(JobSpec::parse(
            r#"{"name":"t","command":{"argv":["x"]},"inputs":[{"path":"../a","content":"x"}]}"#
        )
        .is_err());
        assert!(
            JobSpec::parse(r#"{"name":"t","command":{"argv":["x"]},"inputs":[{"path":"a"}]}"#)
                .is_err()
        );
        assert!(
            JobSpec::parse(r#"{"name":"t","command":{"argv":["x"],"env":{"SPOTLAB_X":"1"}}}"#)
                .is_err()
        );
        assert!(JobSpec::parse(r#"{"name":"t","command":{"argv":["x"]},"bogus":1}"#).is_err());
    }
}
