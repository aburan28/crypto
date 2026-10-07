//! Native control-plane contracts for a proposed `2^40`-record search over
//! curves isogenous to P-256.
//!
//! This module deliberately does **not** implement the isogeny walk.  A bounded
//! native prefix now exists in [`super::p256_isogeny_walk`], with its portable
//! one-unit wrapper in [`super::p256_isogeny_task`], but neither supplies the
//! independently rooted shards this production layout requires.  An S3 receipt
//! is also not evidence that its referenced curves are isogenous.  The types
//! below cover deterministic shard ownership, immutable S3 attempt keys,
//! authoritative completion markers, and the compact audit-reference artifact
//! sent through Cairn's existing commit-reveal transport.
//!
//! The scientific and deployment gates are frozen in
//! `research/p256_isogeny_s3_cairn_20261004/PROTOCOL.md`.

use serde::{Deserialize, Serialize};
use serde_json::{json, Value as Json};

use super::pollard_collab::cairn::commitment_hash;

pub const RUN_SCHEMA: &str = "p256.isogeny.s3.run/v1";
pub const COMPLETE_SCHEMA: &str = "p256.isogeny.s3.shard-complete/v1";
pub const RECEIPT_SCHEMA: &str = "p256.isogeny.cairn.shard-receipt/v1";
pub const P256_ICV1: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
pub const S3_PREFIX: &str = "p256-isogeny-search";

/// The immutable layout of one campaign.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RunManifest {
    pub schema: String,
    pub run_id: String,
    pub curve: String,
    pub shard_bits: u8,
    pub records_per_shard_bits: u8,
    pub config_sha256: String,
    /// Prefix screening is discovery only; it is never a speedup claim.
    pub result_class: String,
}

impl RunManifest {
    pub fn validate(&self) -> Result<(), String> {
        if self.schema != RUN_SCHEMA {
            return Err(format!("schema must be {RUN_SCHEMA}"));
        }
        validate_segment("run_id", &self.run_id)?;
        if self.curve != P256_ICV1 {
            return Err(format!("curve must be registered P-256 ({P256_ICV1})"));
        }
        if self.result_class != "stage-diagnostic" {
            return Err("result_class must be stage-diagnostic".into());
        }
        validate_sha256("config_sha256", &self.config_sha256)?;
        self.shard_count()?;
        self.records_per_shard()?;
        self.total_records()?;
        Ok(())
    }

    /// The production proposal is exactly `2^20` shards of `2^20` emitted
    /// records.  Records are not silently relabelled as globally unique
    /// isomorphism classes.
    pub fn validate_2p40(&self) -> Result<(), String> {
        self.validate()?;
        if self.shard_bits != 20 || self.records_per_shard_bits != 20 {
            return Err("the production layout must be 2^20 shards x 2^20 records".into());
        }
        Ok(())
    }

    pub fn shard_count(&self) -> Result<u64, String> {
        1u64.checked_shl(u32::from(self.shard_bits))
            .ok_or_else(|| "shard_bits is too large".into())
    }

    pub fn records_per_shard(&self) -> Result<u64, String> {
        1u64.checked_shl(u32::from(self.records_per_shard_bits))
            .ok_or_else(|| "records_per_shard_bits is too large".into())
    }

    pub fn total_records(&self) -> Result<u128, String> {
        let bits = u32::from(self.shard_bits) + u32::from(self.records_per_shard_bits);
        1u128
            .checked_shl(bits)
            .ok_or_else(|| "total record exponent is too large".into())
    }

    pub fn run_root(&self) -> Result<String, String> {
        validate_segment("run_id", &self.run_id)?;
        Ok(format!("{S3_PREFIX}/runs/{}", self.run_id))
    }

    pub fn shard_prefix(&self, shard_id: u64) -> Result<String, String> {
        self.require_shard(shard_id)?;
        Ok(format!("{}/shards/{shard_id:08x}", self.run_root()?))
    }

    pub fn complete_key(&self, shard_id: u64) -> Result<String, String> {
        Ok(format!("{}/complete.json", self.shard_prefix(shard_id)?))
    }

    pub fn attempt_prefix(&self, shard_id: u64, attempt_id: &str) -> Result<String, String> {
        validate_segment("attempt_id", attempt_id)?;
        Ok(format!(
            "{}/attempts/{attempt_id}",
            self.shard_prefix(shard_id)?
        ))
    }

    pub fn assignment(&self, task_index: u64, task_count: u64) -> Result<ShardAssignment, String> {
        if task_count == 0 {
            return Err("task_count must be positive".into());
        }
        if task_index >= task_count {
            return Err("task_index must be smaller than task_count".into());
        }
        Ok(ShardAssignment {
            first: task_index,
            step: task_count,
            shard_count: self.shard_count()?,
        })
    }

    fn require_shard(&self, shard_id: u64) -> Result<(), String> {
        if shard_id >= self.shard_count()? {
            return Err(format!("shard_id {shard_id} is outside the run"));
        }
        Ok(())
    }
}

/// A scheduler-neutral strided assignment.  AWS Batch array children can use
/// `(array_index, array_size)` without racing for a mutable queue object.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShardAssignment {
    first: u64,
    step: u64,
    shard_count: u64,
}

impl ShardAssignment {
    pub fn first(&self) -> u64 {
        self.first
    }

    pub fn step(&self) -> u64 {
        self.step
    }

    pub fn shard_count(&self) -> u64 {
        self.shard_count
    }

    pub fn len(&self) -> u64 {
        if self.first >= self.shard_count {
            0
        } else {
            1 + (self.shard_count - 1 - self.first) / self.step
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn shard_at(&self, position: u64) -> Option<u64> {
        let shard = self.first.checked_add(position.checked_mul(self.step)?)?;
        (shard < self.shard_count).then_some(shard)
    }

    pub fn contains(&self, shard_id: u64) -> bool {
        shard_id < self.shard_count
            && shard_id >= self.first
            && (shard_id - self.first).is_multiple_of(self.step)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ObjectRef {
    pub key: String,
    pub sha256: String,
    pub bytes: u64,
}

impl ObjectRef {
    fn validate(&self, label: &str, exact_key: &str) -> Result<(), String> {
        if self.key != exact_key {
            return Err(format!("{label}.key must be {exact_key}"));
        }
        validate_sha256(&format!("{label}.sha256"), &self.sha256)?;
        if self.bytes == 0 {
            return Err(format!("{label}.bytes must be positive"));
        }
        Ok(())
    }
}

/// The only authoritative record that a shard finished.  Writers first upload
/// immutable attempt objects, verify their hashes, then create this marker with
/// S3 `If-None-Match: *`.  Losing attempts remain collectable but never become
/// run state.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ShardComplete {
    pub schema: String,
    pub run_id: String,
    pub config_sha256: String,
    pub shard_id: u64,
    pub attempt_id: String,
    pub records: u64,
    pub distinct_j_within_shard: u64,
    pub raw: ObjectRef,
    pub stage1: ObjectRef,
    pub merkle_root_sha256: String,
}

impl ShardComplete {
    pub fn validate(&self, run: &RunManifest) -> Result<(), String> {
        run.validate()?;
        if self.schema != COMPLETE_SCHEMA {
            return Err(format!("schema must be {COMPLETE_SCHEMA}"));
        }
        if self.run_id != run.run_id {
            return Err("marker run_id does not match RUN.json".into());
        }
        if self.config_sha256 != run.config_sha256 {
            return Err("marker config_sha256 does not match RUN.json".into());
        }
        run.require_shard(self.shard_id)?;
        validate_segment("attempt_id", &self.attempt_id)?;
        let expected = run.records_per_shard()?;
        if self.records != expected {
            return Err(format!(
                "marker has {} records, expected {expected}",
                self.records
            ));
        }
        if self.distinct_j_within_shard != self.records {
            return Err("committed shards must contain no repeated j-invariant".into());
        }
        validate_sha256("merkle_root_sha256", &self.merkle_root_sha256)?;
        let attempt = run.attempt_prefix(self.shard_id, &self.attempt_id)?;
        self.raw
            .validate("raw", &format!("{attempt}/raw.jsonl.zst"))?;
        self.stage1
            .validate("stage1", &format!("{attempt}/stage1.jsonl.zst"))?;
        Ok(())
    }

    /// Produce the deliberately compact Cairn artifact.  This is an
    /// audit/reference receipt, not proof of the S3 object's contents or of
    /// membership in the P-256 isogeny class.
    pub fn cairn_receipt(
        &self,
        run: &RunManifest,
        complete_sha256: &str,
    ) -> Result<CairnShardReceipt, String> {
        self.validate(run)?;
        validate_sha256("complete_sha256", complete_sha256)?;
        Ok(CairnShardReceipt {
            schema: RECEIPT_SCHEMA.into(),
            claim_class: "coordination-receipt".into(),
            verification: "s3-content-address-only".into(),
            run_id: self.run_id.clone(),
            curve: run.curve.clone(),
            shard_id: self.shard_id,
            records: self.records,
            config_sha256: self.config_sha256.clone(),
            complete_key: run.complete_key(self.shard_id)?,
            complete_sha256: complete_sha256.into(),
            raw_sha256: self.raw.sha256.clone(),
            stage1_sha256: self.stage1.sha256.clone(),
            merkle_root_sha256: self.merkle_root_sha256.clone(),
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CairnShardReceipt {
    pub schema: String,
    pub claim_class: String,
    pub verification: String,
    pub run_id: String,
    pub curve: String,
    pub shard_id: u64,
    pub records: u64,
    pub config_sha256: String,
    pub complete_key: String,
    pub complete_sha256: String,
    pub raw_sha256: String,
    pub stage1_sha256: String,
    pub merkle_root_sha256: String,
}

impl CairnShardReceipt {
    pub fn validate(&self, run: &RunManifest) -> Result<(), String> {
        if self.schema != RECEIPT_SCHEMA
            || self.claim_class != "coordination-receipt"
            || self.verification != "s3-content-address-only"
        {
            return Err("receipt overstates or changes its verification class".into());
        }
        if self.run_id != run.run_id || self.curve != run.curve {
            return Err("receipt identity does not match the run".into());
        }
        run.require_shard(self.shard_id)?;
        if self.records != run.records_per_shard()? {
            return Err("receipt record count does not match the run".into());
        }
        if self.config_sha256 != run.config_sha256 {
            return Err("receipt config hash does not match the run".into());
        }
        if self.complete_key != run.complete_key(self.shard_id)? {
            return Err("receipt complete key does not match the shard".into());
        }
        for (name, value) in [
            ("complete_sha256", &self.complete_sha256),
            ("raw_sha256", &self.raw_sha256),
            ("stage1_sha256", &self.stage1_sha256),
            ("merkle_root_sha256", &self.merkle_root_sha256),
        ] {
            validate_sha256(name, value)?;
        }
        Ok(())
    }

    pub fn artifact(&self) -> Json {
        json!({
            "schema": self.schema,
            "claim_class": self.claim_class,
            "verification": self.verification,
            "run_id": self.run_id,
            "curve": self.curve,
            "shard_id": self.shard_id,
            "records": self.records,
            "config_sha256": self.config_sha256,
            "complete_key": self.complete_key,
            "complete_sha256": self.complete_sha256,
            "raw_sha256": self.raw_sha256,
            "stage1_sha256": self.stage1_sha256,
            "merkle_root_sha256": self.merkle_root_sha256,
        })
    }

    /// Stable local idempotency key.  A retry with a different attempt cannot
    /// mint a second local submission for the same run/shard pair.
    pub fn claim_key(&self) -> String {
        format!("{}:{:08x}", self.run_id, self.shard_id)
    }

    /// Use Cairn's consensus-critical commitment encoding.  Signing,
    /// persistence, HTTP submission and epoch-delayed reveal remain the job of
    /// `pollard_collab::cairn::CairnTransport`.
    pub fn commitment_hash(
        &self,
        objective_id: &str,
        submitter: &str,
        nonce: &str,
    ) -> Result<String, String> {
        commitment_hash(objective_id, submitter, &self.artifact(), nonce)
    }
}

fn validate_sha256(name: &str, value: &str) -> Result<(), String> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
    {
        return Err(format!(
            "{name} must be 64 lowercase hexadecimal characters"
        ));
    }
    Ok(())
}

fn validate_segment(name: &str, value: &str) -> Result<(), String> {
    let safe = !value.is_empty()
        && value.len() <= 128
        && value != "."
        && value != ".."
        && value
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, b'.' | b'_' | b'-'));
    if !safe {
        return Err(format!("{name} is not a safe S3 path segment"));
    }
    Ok(())
}
