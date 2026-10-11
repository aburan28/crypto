//! A provider whose "instances" are agent processes on this machine. It
//! exists for tests and rehearsals: preemption is simulated by creating a
//! file (`spotlab local-preempt JOB`), which the agent treats exactly like a
//! cloud's interruption notice.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::os::unix::process::CommandExt;
use std::path::PathBuf;
use std::process::{Command, Stdio};

use super::{static_offers, CatalogEntry, InstState, LaunchRequest, Offer, Provider};
use crate::record::InstanceRef;
use crate::spec::{Placement, Resources};

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct LocalConfig {
    pub catalog: Vec<CatalogEntry>,
    /// The `spotlab` binary to run as the agent; defaults to the running one.
    #[serde(default)]
    pub agent_binary: Option<String>,
    /// Where agents keep their working directories, logs and preempt files.
    pub work_root: String,
}

pub struct Local {
    pub cfg: LocalConfig,
}

impl Local {
    pub fn preempt_file(&self, job_id: &str, attempt: u32) -> PathBuf {
        PathBuf::from(&self.cfg.work_root).join(format!("{job_id}-a{attempt:04}.preempt"))
    }
}

fn pid_of(inst: &InstanceRef) -> Option<i32> {
    inst.instance_id.strip_prefix("pid-")?.parse().ok()
}

impl Provider for Local {
    fn name(&self) -> &str {
        "local"
    }

    fn offers(&self, r: &Resources, p: &Placement) -> Result<Vec<Offer>> {
        Ok(static_offers(
            "local",
            &["local".to_string()],
            &self.cfg.catalog,
            r,
            p,
        ))
    }

    fn launch(&self, offer: &Offer, req: &LaunchRequest) -> Result<InstanceRef> {
        let root = PathBuf::from(&self.cfg.work_root);
        std::fs::create_dir_all(&root)?;
        let bin = match &self.cfg.agent_binary {
            Some(b) => PathBuf::from(b),
            None => std::env::current_exe()?,
        };
        let log = std::fs::File::create(
            root.join(format!("{}-a{:04}.agent.log", req.job_id, req.attempt)),
        )?;
        let child = Command::new(&bin)
            .args(["agent", "--store", &req.store_url, "--job", &req.job_id])
            .args(["--attempt", &req.attempt.to_string(), "--provider", "local"])
            .arg("--workdir")
            .arg(root.join(format!("{}-a{:04}", req.job_id, req.attempt)))
            .arg("--preempt-file")
            .arg(self.preempt_file(&req.job_id, req.attempt))
            .stdin(Stdio::null())
            .stdout(log.try_clone()?)
            .stderr(log)
            .process_group(0)
            .spawn()
            .with_context(|| format!("spawn agent {}", bin.display()))?;
        Ok(InstanceRef {
            provider: "local".into(),
            region: "local".into(),
            zone: offer.zone.clone(),
            instance_type: offer.instance_type.clone(),
            instance_id: format!("pid-{}", child.id()),
        })
    }

    fn state(&self, inst: &InstanceRef) -> Result<InstState> {
        let Some(pid) = pid_of(inst) else {
            return Ok(InstState::Gone("bad id".into()));
        };
        // Reap it if it is our child; otherwise ask whether it still exists.
        let mut status = 0;
        let r = unsafe { libc::waitpid(pid, &mut status, libc::WNOHANG) };
        if r == pid {
            return Ok(InstState::Gone("exited".into()));
        }
        let alive = unsafe { libc::kill(pid, 0) } == 0;
        Ok(if alive {
            InstState::Running
        } else {
            InstState::Gone("exited".into())
        })
    }

    fn terminate(&self, inst: &InstanceRef) -> Result<()> {
        if let Some(pid) = pid_of(inst) {
            unsafe {
                libc::kill(-pid, libc::SIGKILL);
                libc::kill(pid, libc::SIGKILL);
                libc::waitpid(pid, std::ptr::null_mut(), libc::WNOHANG);
            }
        }
        Ok(())
    }
}
