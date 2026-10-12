//! Documents written to the store: job record, heartbeats, exits,
//! checkpoint manifests and the result.

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::store::Store;
use crate::util::{now, sha256_file, walk_files};

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum JobState {
    /// Waiting for an instance (first launch or after a preemption).
    Queued,
    /// An instance was launched for the current attempt.
    Running,
    Succeeded,
    Failed,
    Cancelled,
    /// Attempts, hours or cost exhausted, or the spec cannot be placed.
    Dead,
}

impl JobState {
    pub fn terminal(self) -> bool {
        matches!(
            self,
            JobState::Succeeded | JobState::Failed | JobState::Cancelled | JobState::Dead
        )
    }
    pub fn as_str(self) -> &'static str {
        match self {
            JobState::Queued => "queued",
            JobState::Running => "running",
            JobState::Succeeded => "succeeded",
            JobState::Failed => "failed",
            JobState::Cancelled => "cancelled",
            JobState::Dead => "dead",
        }
    }
}

/// Where an attempt runs.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct InstanceRef {
    pub provider: String,
    pub region: String,
    pub zone: String,
    pub instance_type: String,
    pub instance_id: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct AttemptRecord {
    pub attempt: u32,
    pub instance: InstanceRef,
    pub hourly_usd: f64,
    pub launched_at: f64,
    #[serde(default)]
    pub ended_at: Option<f64>,
    /// `completed`, `preempted`, `lost`, `stale`, `cancelled`, `deadline`, `fenced`.
    #[serde(default)]
    pub end_reason: Option<String>,
}

impl AttemptRecord {
    pub fn hours(&self, at: f64) -> f64 {
        (self.ended_at.unwrap_or(at) - self.launched_at).max(0.0) / 3600.0
    }
    pub fn cost(&self, at: f64) -> f64 {
        self.hours(at) * self.hourly_usd
    }
}

/// The job's state. Created by `submit`; afterwards written only by the
/// controller, so it never needs compare-and-swap.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct JobRecord {
    pub schema: String,
    pub job_id: String,
    pub name: String,
    pub spec_sha256: String,
    pub state: JobState,
    pub created_at: f64,
    pub updated_at: f64,
    /// The attempt currently allowed to run; agents of older attempts are fenced.
    pub attempt: u32,
    #[serde(default)]
    pub attempts: Vec<AttemptRecord>,
    /// Why the job is waiting, or how it ended.
    #[serde(default)]
    pub reason: Option<String>,
    #[serde(default)]
    pub labels: std::collections::BTreeMap<String, String>,
}

impl JobRecord {
    pub fn current(&self) -> Option<&AttemptRecord> {
        self.attempts
            .iter()
            .rev()
            .find(|a| a.attempt == self.attempt)
    }
    pub fn current_mut(&mut self) -> Option<&mut AttemptRecord> {
        let n = self.attempt;
        self.attempts.iter_mut().rev().find(|a| a.attempt == n)
    }
    pub fn total_hours(&self, at: f64) -> f64 {
        self.attempts.iter().map(|a| a.hours(at)).sum()
    }
    pub fn cost_usd(&self, at: f64) -> f64 {
        self.attempts.iter().map(|a| a.cost(at)).sum()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Heartbeat {
    pub attempt: u32,
    pub at: f64,
    pub phase: String,
    pub elapsed_s: f64,
    #[serde(default)]
    pub last_checkpoint_seq: Option<u64>,
}

/// How an attempt ended, from its agent's point of view.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ExitDoc {
    pub attempt: u32,
    /// `completed` (the command exited), `preempted`, `fenced` or `error`.
    pub reason: String,
    #[serde(default)]
    pub exit_code: Option<i32>,
    #[serde(default)]
    pub signal: Option<i32>,
    pub started_at: f64,
    pub ended_at: f64,
    pub wall_s: f64,
    pub user_cpu_s: f64,
    pub sys_cpu_s: f64,
    pub max_rss_kb: i64,
    /// The checkpoint this attempt resumed from, if any.
    #[serde(default)]
    pub restored_seq: Option<u64>,
    /// The last checkpoint this attempt uploaded.
    #[serde(default)]
    pub last_checkpoint_seq: Option<u64>,
    #[serde(default)]
    pub error: Option<String>,
    #[serde(default)]
    pub outputs: Vec<FileEntry>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct FileEntry {
    pub path: String,
    pub sha256: String,
    pub bytes: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct CheckpointManifest {
    pub seq: u64,
    pub attempt: u32,
    pub at: f64,
    pub files: Vec<FileEntry>,
    /// Hash over the sorted (path, sha256) list: identical contents, identical digest.
    pub digest: String,
}

pub fn job_key(id: &str, rest: &str) -> String {
    format!("jobs/{id}/{rest}")
}
pub fn attempt_key(id: &str, attempt: u32, rest: &str) -> String {
    format!("jobs/{id}/attempts/{attempt:04}/{rest}")
}
pub fn ckpt_prefix(id: &str) -> String {
    format!("jobs/{id}/ckpt/")
}
fn ckpt_key(id: &str, seq: u64, rest: &str) -> String {
    format!("jobs/{id}/ckpt/{seq:08}/{rest}")
}

/// Hash files under `dir`, in sorted order.
pub fn hash_dir(dir: &Path) -> Result<(Vec<FileEntry>, String)> {
    let mut files = Vec::new();
    let mut agg = String::new();
    for rel in walk_files(dir)? {
        let p = dir.join(&rel);
        let sha = sha256_file(&p)?;
        let bytes = std::fs::metadata(&p)?.len();
        agg.push_str(&rel);
        agg.push('\0');
        agg.push_str(&sha);
        agg.push('\n');
        files.push(FileEntry {
            path: rel,
            sha256: sha,
            bytes,
        });
    }
    Ok((files, crate::util::sha256_hex(agg.as_bytes())))
}

/// Sequence numbers of complete checkpoints (those with a manifest), ascending.
pub fn checkpoint_seqs(store: &Store, id: &str) -> Result<Vec<u64>> {
    let prefix = ckpt_prefix(id);
    let mut seqs: Vec<u64> = store
        .list(&prefix)?
        .iter()
        .filter_map(|k| k.strip_prefix(&prefix))
        .filter_map(|rest| rest.strip_suffix("/manifest.json"))
        .filter_map(|n| n.parse().ok())
        .collect();
    seqs.sort_unstable();
    seqs.dedup();
    Ok(seqs)
}

/// Upload `dir` as checkpoint `seq`: files first, manifest last.
pub fn upload_checkpoint(
    store: &Store,
    id: &str,
    seq: u64,
    attempt: u32,
    dir: &Path,
) -> Result<CheckpointManifest> {
    let (files, digest) = hash_dir(dir)?;
    for f in &files {
        store.put_file(
            &ckpt_key(id, seq, &format!("files/{}", f.path)),
            &dir.join(&f.path),
        )?;
    }
    // Re-hash: a file that changed while it was uploading would make the
    // manifest disagree with the bytes in the store.
    let (_, after) = hash_dir(dir)?;
    if after != digest {
        bail!("checkpoint directory changed during upload; write checkpoint files atomically (rename)");
    }
    let m = CheckpointManifest {
        seq,
        attempt,
        at: now(),
        files,
        digest,
    };
    store.put_json(&ckpt_key(id, seq, "manifest.json"), &m)?;
    Ok(m)
}

/// Download checkpoint `seq` into `dir`, verifying every file's hash.
pub fn restore_checkpoint(
    store: &Store,
    id: &str,
    seq: u64,
    dir: &Path,
) -> Result<CheckpointManifest> {
    let m: CheckpointManifest = store
        .get_json(&ckpt_key(id, seq, "manifest.json"))?
        .with_context(|| format!("checkpoint {seq} has no manifest"))?;
    for f in &m.files {
        crate::util::safe_rel(&f.path)?;
        let dst = dir.join(&f.path);
        if !store.get_file(&ckpt_key(id, seq, &format!("files/{}", f.path)), &dst)? {
            bail!("checkpoint {seq}: {} missing from the store", f.path);
        }
        let got = sha256_file(&dst)?;
        if got != f.sha256 {
            bail!(
                "checkpoint {seq}: {} hash {got} does not match manifest {}",
                f.path,
                f.sha256
            );
        }
    }
    Ok(m)
}

/// Delete all but the newest `keep` checkpoints.
pub fn prune_checkpoints(store: &Store, id: &str, keep: usize) -> Result<()> {
    let seqs = checkpoint_seqs(store, id)?;
    if seqs.len() <= keep {
        return Ok(());
    }
    for seq in &seqs[..seqs.len() - keep] {
        // Manifest first, so a half-deleted checkpoint is never mistaken for a complete one.
        store.delete(&ckpt_key(id, *seq, "manifest.json"))?;
        for k in store.list(&ckpt_key(id, *seq, ""))? {
            store.delete(&k)?;
        }
    }
    Ok(())
}

/// The result document, composed by the controller when a job ends.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ResultDoc {
    pub schema: String,
    pub job_id: String,
    pub name: String,
    pub spec_sha256: String,
    pub status: String,
    pub exit_code: Option<i32>,
    pub reason: Option<String>,
    pub attempts: Vec<AttemptSummary>,
    /// More than one attempt ran and a later one restored a checkpoint.
    pub resumed: bool,
    pub preemptions: u32,
    pub total_hours: f64,
    pub cost_usd_estimate: f64,
    pub outputs: Vec<FileEntry>,
    pub evidence: String,
    pub labels: std::collections::BTreeMap<String, String>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct AttemptSummary {
    pub attempt: u32,
    pub instance: InstanceRef,
    pub hourly_usd: f64,
    pub launched_at: f64,
    pub ended_at: Option<f64>,
    pub end_reason: Option<String>,
    pub exit: Option<ExitDoc>,
}

/// What a spot run is worth as evidence, stated in every result.
pub const EVIDENCE: &str = "throughput run on spot virtual machines: wall-clock and CPU times are \
not performance measurements (shared hosts, no isolation, and a resumed run's segments ran on \
different machines); counted units the job reports itself are admissible only if the job \
accumulates them across checkpoints";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checkpoint_round_trip_and_prune() {
        let base = std::env::temp_dir().join(format!("spotlab-ck-{}", std::process::id()));
        let store = Store::open(base.join("store").to_str().unwrap()).unwrap();
        let src = base.join("src");
        std::fs::create_dir_all(src.join("sub")).unwrap();
        std::fs::write(src.join("state.bin"), b"hello").unwrap();
        std::fs::write(src.join("sub/x"), b"world").unwrap();
        for seq in 1..=3 {
            upload_checkpoint(&store, "j1", seq, 1, &src).unwrap();
        }
        assert_eq!(checkpoint_seqs(&store, "j1").unwrap(), vec![1, 2, 3]);
        prune_checkpoints(&store, "j1", 2).unwrap();
        assert_eq!(checkpoint_seqs(&store, "j1").unwrap(), vec![2, 3]);
        let dst = base.join("dst");
        let m = restore_checkpoint(&store, "j1", 3, &dst).unwrap();
        assert_eq!(m.files.len(), 2);
        assert_eq!(std::fs::read(dst.join("sub/x")).unwrap(), b"world");
        // Tampering is caught.
        store
            .put("jobs/j1/ckpt/00000003/files/state.bin", b"evil")
            .unwrap();
        assert!(restore_checkpoint(&store, "j1", 3, &base.join("dst2")).is_err());
        std::fs::remove_dir_all(base).ok();
    }
}
