//! Submission to an independent runner: an `isolab.job/v1` for a spec.
//!
//! isolab (`isolab/README.md`) places a job on a lab worker: a cgroup v2
//! isolated partition of whole cores on one NUMA node, other threads and
//! IRQs moved off, frequency policy set, quiet gated by PSI, and a result
//! committed once with the worker's own observations.  `ecbench` inside
//! the job then runs with `--cpus inherit`: it takes the placement isolab
//! made, pins its measured children inside it, and grades every run from
//! its own observations as it would anywhere else.  The session directory
//! comes back as hashed artifacts, and `verify` re-derives it on the
//! worker after the clock stops.
//!
//! Two ways to get the binary there:
//!
//! - **binary** (preferred): upload a binary built for the worker's
//!   architecture.  Its hash is the code's identity, the job needs no
//!   network, and nothing is resolved on the worker.
//! - **source**: export the repository at a full commit and build it on
//!   the worker.  The repository does not track `Cargo.lock`, so the
//!   build resolves dependencies afresh and needs network; the job says
//!   so in its labels, and the binary's hash in every record remains the
//!   identity a comparison checks.

use serde_json::{json, Value};

use crate::cryptanalysis::ecbench::spec::Plan;

/// How the binary reaches the worker.
#[derive(Clone, Debug)]
pub enum Delivery {
    /// A local binary, uploaded by the isolab client.
    Binary { path: String },
    /// A git URL and full 40-hex commit, built on the worker.
    Source { url: String, commit: String },
}

#[derive(Clone, Debug)]
pub struct JobOptions {
    pub pool: String,
    pub delivery: Delivery,
    /// Logical CPUs to place: 2 is one core with its sibling under SMT;
    /// the measured child takes one, the runner sleeps on the other.
    pub cpus: u32,
    pub memory_mb: Option<u64>,
    pub image: Option<String>,
    pub timeout_seconds: u64,
    /// `strict`, `standard` or `best_effort` (isolab's fidelity policy).
    pub policy: String,
}

/// The job document.  `spec_text` travels inline as `spec.json`.
pub fn job(plan: &Plan, spec_text: &str, o: &JobOptions) -> Result<Value, String> {
    if !matches!(o.policy.as_str(), "strict" | "standard" | "best_effort") {
        return Err(format!(
            "policy `{}` is not strict, standard or best_effort",
            o.policy
        ));
    }
    let mut inputs = vec![json!({"path": "spec.json", "content": spec_text, "encoding": "utf8"})];
    let mut build = Vec::new();
    let mut labels = json!({
        "ecbench_spec_id": plan.spec_id,
        "ecbench_spec_sha256": plan.spec_sha256,
        "ecbench_label": plan.spec.label.chars().take(100).collect::<String>(),
    });
    let (bin, network) = match &o.delivery {
        Delivery::Binary { path } => {
            inputs.push(json!({"path": "ecbench", "local_file": path, "mode": "755"}));
            labels["ecbench_delivery"] = json!("binary");
            ("./ecbench".to_string(), "none")
        }
        Delivery::Source { url, commit } => {
            if commit.len() != 40 || !commit.chars().all(|c| c.is_ascii_hexdigit()) {
                return Err("source delivery needs a full 40-hex commit".into());
            }
            inputs.push(json!({"path": "src", "git": {"url": url, "commit": commit}}));
            build.push(json!({
                "argv": ["cargo", "build", "--release", "--bin", "ecbench"],
                "cwd": "src",
                "timeout_s": 3600,
            }));
            labels["ecbench_delivery"] = json!("source");
            labels["git_commit"] = json!(commit);
            labels["ecbench_dependencies"] =
                json!("resolved on the worker: Cargo.lock is not tracked");
            ("src/target/release/ecbench".to_string(), "bridge")
        }
    };
    let timeout = plan.executions.len() as u64 * plan.spec.measurement.timeout_seconds.max(1) + 600;
    let mut doc = json!({
        "schema": "isolab.job/v1",
        "name": format!("ecbench {}", plan.spec.label).chars().take(120).collect::<String>(),
        "pool": o.pool,
        "runtime": {"backend": "auto", "network": network},
        "inputs": inputs,
        "build": build,
        "command": {
            "argv": [bin, "run", "--spec", "spec.json", "--out", "${ISOLAB_OUTPUT_DIR}/session",
                     "--cpus", "inherit", "--allow-busy"],
            "env": {"RAYON_NUM_THREADS": "1"},
        },
        "measure": {"warmups": 0, "repeats": 1},
        "resources": {"cpus": o.cpus, "numa_node": "single", "smt": "isolate"},
        "fidelity": {
            "policy": o.policy,
            "require": ["cpu_partition", "numa_bind", "evicted", "no_steal"],
        },
        "limits": {"timeout_s": o.timeout_seconds.max(timeout)},
        "verify": {
            "argv": [bin_for_verify(&o.delivery), "verify", "--dir", "${ISOLAB_OUTPUT_DIR}/session",
                     "--replay", "3", "--exit-code"],
            "timeout_s": 3600,
        },
        "labels": labels,
    });
    if let Some(m) = o.memory_mb {
        doc["resources"]["memory_mb"] = json!(m);
    }
    if let Some(img) = &o.image {
        doc["runtime"]["image"] = json!(img);
    }
    Ok(doc)
}

fn bin_for_verify(d: &Delivery) -> &'static str {
    match d {
        Delivery::Binary { .. } => "./ecbench",
        Delivery::Source { .. } => "src/target/release/ecbench",
    }
}

/// `${VAR}` in a path, expanded from the environment: isolab hands the
/// output directory as `$ISOLAB_OUTPUT_DIR`, and argv is not a shell.
pub fn expand_env(s: &str) -> Result<String, String> {
    let mut out = String::new();
    let mut rest = s;
    while let Some(start) = rest.find("${") {
        out.push_str(&rest[..start]);
        let after = &rest[start + 2..];
        let end = after
            .find('}')
            .ok_or_else(|| format!("unclosed ${{ in `{s}`"))?;
        let name = &after[..end];
        let val = std::env::var(name)
            .map_err(|_| format!("`{s}` names ${{{name}}}, which is not set"))?;
        out.push_str(&val);
        rest = &after[end + 1..];
    }
    out.push_str(rest);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn env_expansion() {
        std::env::set_var("ECBENCH_TEST_DIR", "/x/y");
        assert_eq!(
            expand_env("${ECBENCH_TEST_DIR}/session").unwrap(),
            "/x/y/session"
        );
        assert_eq!(expand_env("plain").unwrap(), "plain");
        assert!(expand_env("${ECBENCH_SURELY_UNSET_VAR}").is_err());
    }
}
