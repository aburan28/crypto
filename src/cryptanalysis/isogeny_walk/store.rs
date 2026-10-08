//! S3 storage for walk and trait outputs, write-once and hash-checked.
//!
//! The layout follows PR #1330's contract
//! (`src/cryptanalysis/p256_isogeny_campaign.rs`), generalised from its
//! fixed `2^40` campaign to any walk:
//!
//! ```text
//! s3://<bucket>/<prefix>/runs/<run_id>/
//!   attempts/<attempt_id>/<file>.gz      immutable, If-None-Match: *
//!   complete.json                        the only authoritative record
//!   traits/<of>/shard-<i>/attempts/…, complete.json
//! ```
//!
//! - **Write-once objects.** Every object is created with
//!   `If-None-Match: *`, so nothing is ever overwritten, whether or not the
//!   bucket keeps versions.
//! - **Checked uploads.** Each upload's SHA-256, as S3 reports it, must equal
//!   the local hash before the marker is written.
//! - **One winner.** The marker lists every object by key, stored SHA-256,
//!   byte count and uncompressed SHA-256.  If two attempts race, one marker
//!   wins and the loser's objects stay collectable, never authoritative.
//! - **Checked downloads.** [`fetch`] recomputes both hashes of every object
//!   before writing it.
//!
//! The transport is the AWS CLI (`aws s3api`), run as a subprocess with the
//! caller's credentials; no SDK is linked.

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use base64::Engine;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use super::record::{sha256_hex, V};
use crate::hash::sha256::sha256;

pub const COMPLETE_SCHEMA: &str = "isogeny-walk.s3.complete/v1";
pub const COMPLETE_FILE: &str = "complete.json";

/// `s3://bucket/prefix`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct S3Loc {
    pub bucket: String,
    pub prefix: String,
}

impl S3Loc {
    pub fn parse(uri: &str) -> Result<Self, String> {
        let rest = uri
            .strip_prefix("s3://")
            .ok_or_else(|| format!("{uri}: expected s3://bucket/prefix"))?;
        let (bucket, prefix) = rest.split_once('/').unwrap_or((rest, ""));
        if bucket.is_empty()
            || !bucket
                .bytes()
                .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == b'-' || c == b'.')
        {
            return Err(format!("{uri}: bad bucket name"));
        }
        Ok(S3Loc {
            bucket: bucket.to_string(),
            prefix: prefix.trim_matches('/').to_string(),
        })
    }

    /// The key under this location's prefix.
    pub fn key(&self, rel: &str) -> String {
        if self.prefix.is_empty() {
            rel.to_string()
        } else {
            format!("{}/{rel}", self.prefix)
        }
    }

    pub fn uri(&self, rel: &str) -> String {
        format!("s3://{}/{}", self.bucket, self.key(rel))
    }
}

/// A path segment: lower-case letters, digits, `-`, `_`, `.`, no `..`.
pub fn segment(s: &str) -> Result<&str, String> {
    if s.is_empty()
        || s.contains("..")
        || !s
            .bytes()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || b"-_.".contains(&c))
    {
        return Err(format!("{s:?} is not a safe key segment"));
    }
    Ok(s)
}

/// A fresh attempt id: the taskq task id when running under taskq, else
/// the host and the time.
pub fn attempt_id() -> String {
    let base = std::env::var("TASKQ_TASK_ID")
        .ok()
        .or_else(|| std::env::var("HOSTNAME").ok())
        .unwrap_or_else(|| "local".into());
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let clean: String = base
        .to_ascii_lowercase()
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' {
                c
            } else {
                '-'
            }
        })
        .collect();
    format!("{}-{nanos}", clean.trim_matches('-'))
}

fn aws(args: &[&str]) -> Result<(bool, String, String), String> {
    let out = Command::new("aws")
        .args(args)
        .output()
        .map_err(|e| format!("cannot run the AWS CLI: {e}"))?;
    Ok((
        out.status.success(),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    ))
}

fn b64_sha256(bytes: &[u8]) -> String {
    base64::engine::general_purpose::STANDARD.encode(sha256(bytes))
}

/// What a write-once put found.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Put {
    Created,
    /// The key already existed; nothing was written.
    Exists,
}

/// Metadata S3 reports for an object uploaded with a SHA-256 checksum.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Head {
    pub bytes: u64,
    pub sha256: String,
}

fn missing(err: &str) -> bool {
    err.contains("NoSuchKey") || err.contains("Not Found") || err.contains("404")
}

fn checksum_hex(value: &str) -> Result<String, String> {
    let bytes = base64::engine::general_purpose::STANDARD
        .decode(value)
        .map_err(|e| format!("invalid S3 SHA-256 checksum: {e}"))?;
    if bytes.len() != 32 {
        return Err(format!(
            "S3 SHA-256 checksum decoded to {} bytes, expected 32",
            bytes.len()
        ));
    }
    Ok(hex::encode(bytes))
}

/// Read an object's exact byte count and stored SHA-256 checksum.
pub fn head(loc: &S3Loc, key: &str) -> Result<Option<Head>, String> {
    let (ok, out, err) = aws(&[
        "s3api",
        "head-object",
        "--bucket",
        &loc.bucket,
        "--key",
        key,
        "--checksum-mode",
        "ENABLED",
        "--output",
        "json",
    ])?;
    if !ok {
        if missing(&err) {
            return Ok(None);
        }
        return Err(format!("head s3://{}/{key}: {}", loc.bucket, err.trim()));
    }
    let value: serde_json::Value = serde_json::from_str(&out).map_err(|e| e.to_string())?;
    let bytes = value["ContentLength"]
        .as_u64()
        .ok_or("S3 head response has no ContentLength")?;
    let encoded = value["ChecksumSHA256"]
        .as_str()
        .ok_or("S3 head response has no ChecksumSHA256")?;
    Ok(Some(Head {
        bytes,
        sha256: checksum_hex(encoded)?,
    }))
}

/// Upload an existing file without buffering it, conditional on the key not
/// existing. An existing content-addressed object is accepted only when its
/// byte count and S3-reported SHA-256 match the caller's descriptor.
pub fn put_file_once(
    loc: &S3Loc,
    key: &str,
    path: &Path,
    expected_sha256: &str,
    expected_bytes: u64,
) -> Result<Put, String> {
    let metadata = std::fs::metadata(path).map_err(|e| format!("stat {}: {e}", path.display()))?;
    if metadata.len() != expected_bytes {
        return Err(format!(
            "{} has {} bytes, expected {expected_bytes}",
            path.display(),
            metadata.len()
        ));
    }
    let digest = hex::decode(expected_sha256)
        .map_err(|e| format!("invalid expected SHA-256 {expected_sha256}: {e}"))?;
    if digest.len() != 32 {
        return Err("expected SHA-256 must contain 32 bytes".into());
    }
    let expected_b64 = base64::engine::general_purpose::STANDARD.encode(digest);
    let body = path.to_string_lossy().into_owned();
    let (ok, out, err) = aws(&[
        "s3api",
        "put-object",
        "--bucket",
        &loc.bucket,
        "--key",
        key,
        "--body",
        &body,
        "--if-none-match",
        "*",
        "--checksum-algorithm",
        "SHA256",
        "--checksum-sha256",
        &expected_b64,
        "--output",
        "json",
    ])?;
    if !ok {
        if !err.contains("PreconditionFailed") && !err.contains("ConditionalRequestConflict") {
            return Err(format!("put s3://{}/{key}: {}", loc.bucket, err.trim()));
        }
        let existing = head(loc, key)?.ok_or_else(|| format!("{key}: object vanished"))?;
        if existing.bytes != expected_bytes || existing.sha256 != expected_sha256 {
            return Err(format!(
                "s3://{}/{key}: existing object differs from its content address",
                loc.bucket
            ));
        }
        return Ok(Put::Exists);
    }
    let value: serde_json::Value = serde_json::from_str(&out).map_err(|e| e.to_string())?;
    if value["ChecksumSHA256"].as_str() != Some(expected_b64.as_str()) {
        return Err(format!(
            "put s3://{}/{key}: S3 reports SHA-256 {:?}, expected {expected_b64}",
            loc.bucket, value["ChecksumSHA256"]
        ));
    }
    Ok(Put::Created)
}

/// Download an object directly to a new local file. Returns `false` when the
/// key does not exist and never exposes a successful partial download.
pub fn get_to_file(loc: &S3Loc, key: &str, output: &Path) -> Result<bool, String> {
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()));
    }
    if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent).map_err(|e| format!("create {}: {e}", parent.display()))?;
    }
    let outfile = output.to_string_lossy().into_owned();
    let (ok, _, err) = aws(&[
        "s3api",
        "get-object",
        "--bucket",
        &loc.bucket,
        "--key",
        key,
        &outfile,
    ])?;
    if !ok {
        let _ = std::fs::remove_file(output);
        if missing(&err) {
            return Ok(false);
        }
        return Err(format!("get s3://{}/{key}: {}", loc.bucket, err.trim()));
    }
    File::open(output)
        .and_then(|file| file.sync_all())
        .map_err(|e| format!("sync {}: {e}", output.display()))?;
    Ok(true)
}

/// Create `key` with `bytes` unless it exists, and check S3's SHA-256 of
/// what it stored against ours.
pub fn put_once(loc: &S3Loc, key: &str, bytes: &[u8], scratch: &Path) -> Result<Put, String> {
    std::fs::create_dir_all(scratch).map_err(|e| e.to_string())?;
    let tmp = scratch.join(format!("put-{}", sha256_hex(key)));
    std::fs::write(&tmp, bytes).map_err(|e| e.to_string())?;
    let body = tmp.to_string_lossy().into_owned();
    let (ok, out, err) = aws(&[
        "s3api",
        "put-object",
        "--bucket",
        &loc.bucket,
        "--key",
        key,
        "--body",
        &body,
        "--if-none-match",
        "*",
        "--checksum-algorithm",
        "SHA256",
        "--output",
        "json",
    ])?;
    let _ = std::fs::remove_file(&tmp);
    if !ok {
        if err.contains("PreconditionFailed") {
            return Ok(Put::Exists);
        }
        return Err(format!("put s3://{}/{key}: {}", loc.bucket, err.trim()));
    }
    let v: serde_json::Value = serde_json::from_str(&out).map_err(|e| e.to_string())?;
    let want = b64_sha256(bytes);
    if v["ChecksumSHA256"].as_str() != Some(want.as_str()) {
        return Err(format!(
            "put s3://{}/{key}: S3 reports SHA-256 {:?}, expected {want}",
            loc.bucket, v["ChecksumSHA256"]
        ));
    }
    Ok(Put::Created)
}

/// The bytes at `key`, or `None` if it does not exist.
pub fn get(loc: &S3Loc, key: &str, scratch: &Path) -> Result<Option<Vec<u8>>, String> {
    std::fs::create_dir_all(scratch).map_err(|e| e.to_string())?;
    let tmp = scratch.join(format!("get-{}", sha256_hex(key)));
    let outfile = tmp.to_string_lossy().into_owned();
    let (ok, _, err) = aws(&[
        "s3api",
        "get-object",
        "--bucket",
        &loc.bucket,
        "--key",
        key,
        &outfile,
    ])?;
    if !ok {
        if err.contains("NoSuchKey") || err.contains("Not Found") || err.contains("404") {
            return Ok(None);
        }
        return Err(format!("get s3://{}/{key}: {}", loc.bucket, err.trim()));
    }
    let bytes = std::fs::read(&tmp).map_err(|e| e.to_string())?;
    let _ = std::fs::remove_file(&tmp);
    Ok(Some(bytes))
}

fn gzip(bytes: &[u8]) -> Result<Vec<u8>, String> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(bytes).map_err(|e| e.to_string())?;
    enc.finish().map_err(|e| e.to_string())
}

fn gunzip(bytes: &[u8]) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    GzDecoder::new(bytes)
        .read_to_end(&mut out)
        .map_err(|e| e.to_string())?;
    Ok(out)
}

/// What [`publish`] did.
pub struct Published {
    /// The marker now at `complete.json` (ours, or the winner's).
    pub marker: V,
    pub marker_uri: String,
    /// `true` if our attempt wrote the marker.
    pub won: bool,
}

/// Upload `files` (local paths, stored gzip-compressed under their file
/// names) as one attempt of `run_rel`, then create `run_rel/complete.json`.
/// `meta` is recorded in the marker.
pub fn publish(
    loc: &S3Loc,
    run_rel: &str,
    files: &[PathBuf],
    meta: V,
    scratch: &Path,
) -> Result<Published, String> {
    // A completed run is final: skip the upload when its marker exists.
    // A racing attempt can still lose at the conditional marker write below.
    let marker_key = loc.key(&format!("{run_rel}/{COMPLETE_FILE}"));
    if let Some(theirs) = get(loc, &marker_key, scratch)? {
        let theirs: serde_json::Value =
            serde_json::from_slice(&theirs).map_err(|e| e.to_string())?;
        return Ok(Published {
            marker: V::s(theirs.to_string()),
            marker_uri: format!("s3://{}/{marker_key}", loc.bucket),
            won: false,
        });
    }
    let attempt = attempt_id();
    let mut objects = Vec::new();
    for path in files {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .ok_or_else(|| format!("{}: no file name", path.display()))?;
        segment(name)?;
        let content = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let stored = gzip(&content)?;
        let key = loc.key(&format!("{run_rel}/attempts/{attempt}/{name}.gz"));
        if put_once(loc, &key, &stored, scratch)? == Put::Exists {
            return Err(format!("{key} already exists; attempt ids must be fresh"));
        }
        objects.push(V::map(vec![
            ("name", V::s(name)),
            ("key", V::s(key)),
            ("encoding", V::s("gzip")),
            ("sha256", V::s(sha256_hex_bytes(&stored))),
            ("bytes", V::int(stored.len())),
            ("content_sha256", V::s(sha256_hex_bytes(&content))),
            ("content_bytes", V::int(content.len())),
        ]));
    }
    let marker = V::map(vec![
        ("schema", V::s(COMPLETE_SCHEMA)),
        ("run", V::s(loc.uri(run_rel))),
        ("attempt_id", V::s(attempt)),
        ("objects", V::Seq(objects)),
        ("meta", meta),
    ]);
    let text = marker.json();
    let marker_uri = format!("s3://{}/{marker_key}", loc.bucket);
    match put_once(loc, &marker_key, text.as_bytes(), scratch)? {
        Put::Created => Ok(Published {
            marker,
            marker_uri,
            won: true,
        }),
        Put::Exists => {
            let theirs = get(loc, &marker_key, scratch)?.ok_or("marker vanished")?;
            let theirs: serde_json::Value =
                serde_json::from_slice(&theirs).map_err(|e| e.to_string())?;
            Ok(Published {
                marker: V::s(theirs.to_string()),
                marker_uri,
                won: false,
            })
        }
    }
}

fn sha256_hex_bytes(b: &[u8]) -> String {
    sha256(b).iter().map(|x| format!("{x:02x}")).collect()
}

/// Download the completed run at `run_rel` into `out`, checking every
/// object's stored and uncompressed SHA-256.  Returns the marker.
pub fn fetch(
    loc: &S3Loc,
    run_rel: &str,
    out: &Path,
    scratch: &Path,
) -> Result<serde_json::Value, String> {
    let marker_key = loc.key(&format!("{run_rel}/{COMPLETE_FILE}"));
    let text = get(loc, &marker_key, scratch)?
        .ok_or_else(|| format!("s3://{}/{marker_key}: no completed run", loc.bucket))?;
    let marker: serde_json::Value = serde_json::from_slice(&text).map_err(|e| e.to_string())?;
    if marker["schema"].as_str() != Some(COMPLETE_SCHEMA) {
        return Err(format!("{marker_key}: not a {COMPLETE_SCHEMA} marker"));
    }
    std::fs::create_dir_all(out).map_err(|e| e.to_string())?;
    for o in marker["objects"].as_array().ok_or("objects")? {
        let name = segment(o["name"].as_str().ok_or("name")?)?.to_string();
        let key = o["key"].as_str().ok_or("key")?;
        if !key.starts_with(&loc.key(&format!("{run_rel}/attempts/"))) {
            return Err(format!("{key}: outside this run"));
        }
        let stored = get(loc, key, scratch)?.ok_or_else(|| format!("{key}: missing"))?;
        if Some(sha256_hex_bytes(&stored).as_str()) != o["sha256"].as_str() {
            return Err(format!("{key}: stored SHA-256 differs from the marker"));
        }
        let content = match o["encoding"].as_str() {
            Some("gzip") => gunzip(&stored)?,
            Some("identity") => stored,
            other => return Err(format!("{key}: unknown encoding {other:?}")),
        };
        if Some(sha256_hex_bytes(&content).as_str()) != o["content_sha256"].as_str() {
            return Err(format!("{key}: content SHA-256 differs from the marker"));
        }
        std::fs::write(out.join(&name), &content).map_err(|e| e.to_string())?;
    }
    Ok(marker)
}

/// `runs/<run_id>` for a walk.
pub fn walk_rel(run_id: &str) -> Result<String, String> {
    Ok(format!("runs/{}", segment(run_id)?))
}

/// `runs/<run_id>/traits/<of>/shard-<i>` for a trait shard.
pub fn shard_rel(run_id: &str, shard: usize, of: usize) -> Result<String, String> {
    Ok(format!(
        "{}/traits/{of}/shard-{shard:04}",
        walk_rel(run_id)?
    ))
}

/// `runs/<run_id>/traits/<of>/collected` for merged traits.
pub fn collected_rel(run_id: &str, of: usize) -> Result<String, String> {
    Ok(format!("{}/traits/{of}/collected", walk_rel(run_id)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locations_and_segments() {
        let l = S3Loc::parse("s3://crypto-autoresearcher/isogeny-walk/").unwrap();
        assert_eq!(l.bucket, "crypto-autoresearcher");
        assert_eq!(
            l.key("runs/x/complete.json"),
            "isogeny-walk/runs/x/complete.json"
        );
        assert!(S3Loc::parse("https://x").is_err());
        assert!(segment("../x").is_err() && segment("A").is_err() && segment("p256-ab12").is_ok());
        assert_eq!(
            shard_rel("p256-ab", 3, 8).unwrap(),
            "runs/p256-ab/traits/8/shard-0003"
        );
        let data = b"isogeny".repeat(100);
        assert_eq!(gunzip(&gzip(&data).unwrap()).unwrap(), data);
        // S3's ChecksumSHA256 is base64 of the raw digest.
        assert_eq!(
            b64_sha256(b"probe\n"),
            "Jb4yNVba03ertX/n7IxLmaZSf0iN2ijQybaGUoZZyQk="
        );
    }
}
