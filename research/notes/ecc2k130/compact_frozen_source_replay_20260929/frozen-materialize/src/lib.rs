//! Materialize the frozen source tree of a compact-orbit evaluation.
//!
//! A compact-orbit `FROZEN.json` pins the SHA-256 of the files its evaluation
//! builds and runs, but only the files the evaluation's own PR touched. The
//! rest of the build (the other ~500 files of `src/`, the data files `src/`
//! embeds, and the Cargo lock the run used) was whatever `main` held at the
//! time. Replaying from the live checkout therefore mixes frozen files with
//! newer ones, and the mix stops compiling as soon as `main` moves: a pinned
//! `Cargo.toml` without `redis` beside a live module that imports it, a pinned
//! index-calculus crate beside a live module that calls a method added later.
//!
//! The materialized tree is built in three steps, each fail-closed:
//!
//! 1. copy the checkout, without `.git` and `target`;
//! 2. replace the unpinned inputs (`src/`, the embedded data files, the lock)
//!    by the committed archive of the base commit the freeze names, whose
//!    SHA-256 is recorded in `BASES.json`;
//! 3. replace every pinned file whose live bytes drifted by its committed
//!    gzip snapshot, then re-hash every pinned file against the freeze.
//!
//! This is the native replacement of `materialize.py` (step 3 is the same
//! contract; step 2 is new). `materialize.py` stays for the older consumers
//! that import it.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

/// Largest snapshot or base archive member read into memory. The cap guards
/// against a corrupted gzip stream, not a size any source approaches.
pub const MAX_SNAPSHOT_BYTES: u64 = 8 * 1024 * 1024;
/// Cap on a whole base archive after decompression (src is about 12 MB).
pub const MAX_BASE_BYTES: u64 = 256 * 1024 * 1024;

pub const SNAPSHOTS_SCHEMA: &str = "compact-frozen-source-snapshots-v1";
pub const BASES_SCHEMA: &str = "compact-frozen-source-bases-v1";
pub const RECEIPT_SCHEMA: &str = "compact-frozen-source-materialization-v2";

/// The inputs of a build that no freeze pins, restored from the base commit.
/// A path that does not exist at the base commit is removed from the tree, so
/// a live file that did not exist then cannot leak in.
pub const BASE_PATHS: &[&str] = &[
    "src",
    "docs/ecbench/schema.sql",
    "docs/curves/registry.json",
    "docs/ic/boundary_targets.json",
    "research/notes/ecc2k130/compact_shared_log_20260925/Cargo.lock",
];

type Res<T> = Result<T, String>;

pub fn sha_hex(data: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(data);
    h.finalize().iter().map(|b| format!("{b:02x}")).collect()
}

fn read(path: &Path) -> Res<Vec<u8>> {
    fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}

fn read_json(path: &Path) -> Res<Value> {
    serde_json::from_slice(&read(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

fn str_field<'a>(v: &'a Value, key: &str, ctx: &str) -> Res<&'a str> {
    v.get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("{ctx}: missing string field {key:?}"))
}

fn u64_field(v: &Value, key: &str, ctx: &str) -> Res<u64> {
    v.get(key)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("{ctx}: missing integer field {key:?}"))
}

/// Gunzip with a hard cap on the decompressed size.
fn gunzip_capped(packed: &[u8], cap: u64, what: &str) -> Res<Vec<u8>> {
    let mut out = Vec::new();
    GzDecoder::new(packed)
        .take(cap + 1)
        .read_to_end(&mut out)
        .map_err(|e| format!("{what}: {e}"))?;
    if out.len() as u64 > cap {
        return Err(format!("{what}: larger than {cap} bytes"));
    }
    Ok(out)
}

// ───────────────────────────── pins and snapshots ─────────────────────────────

/// Every pinned path and digest of the given freezes. Two freezes that pin one
/// path to different digests are refused.
pub fn load_pins(freezes: &[PathBuf]) -> Res<BTreeMap<String, String>> {
    let mut pins: BTreeMap<String, String> = BTreeMap::new();
    for freeze in freezes {
        let doc = read_json(freeze)?;
        let table = doc
            .get("source_sha256")
            .and_then(Value::as_object)
            .ok_or_else(|| format!("{}: no source_sha256", freeze.display()))?;
        for (path, digest) in table {
            let digest = digest
                .as_str()
                .ok_or_else(|| format!("{path}: pin is not a string"))?;
            if let Some(prev) = pins.insert(path.clone(), digest.to_string()) {
                if prev != digest {
                    return Err(format!("{path}: conflicting pins"));
                }
            }
        }
    }
    Ok(pins)
}

pub struct Snapshots {
    dir: PathBuf,
    by_digest: BTreeMap<String, Value>,
}

impl Snapshots {
    pub fn load(replay_dir: &Path) -> Res<Self> {
        let manifest = read_json(&replay_dir.join("SNAPSHOTS.json"))?;
        if manifest.get("schema").and_then(Value::as_str) != Some(SNAPSHOTS_SCHEMA) {
            return Err("SNAPSHOTS.json: unexpected schema".into());
        }
        let mut by_digest = BTreeMap::new();
        for entry in manifest
            .get("snapshots")
            .and_then(Value::as_array)
            .ok_or("SNAPSHOTS.json: no snapshots")?
        {
            by_digest.insert(
                str_field(entry, "sha256", "snapshot")?.to_string(),
                entry.clone(),
            );
        }
        Ok(Self {
            dir: replay_dir.join("historical_snapshots"),
            by_digest,
        })
    }

    /// The frozen bytes of `path`, accepted only if the manifest lists the
    /// snapshot under that path and its gzip hash, length and content hash all
    /// match.
    pub fn bytes(&self, path: &str, digest: &str) -> Res<Vec<u8>> {
        let entry = self
            .by_digest
            .get(digest)
            .ok_or_else(|| format!("{path}: drifted and no snapshot of {digest}"))?;
        let listed = str_field(entry, "path", path)?;
        if listed != path {
            return Err(format!("{path}: snapshot {digest} is {listed}"));
        }
        let file = self.dir.join(format!("{digest}.gz"));
        let packed = read(&file)?;
        if sha_hex(&packed) != str_field(entry, "gzip_sha256", path)? {
            return Err(format!("{path}: snapshot gzip hash"));
        }
        let raw = gunzip_capped(&packed, MAX_SNAPSHOT_BYTES, path)?;
        if raw.len() as u64 != u64_field(entry, "bytes", path)? {
            return Err(format!("{path}: snapshot length"));
        }
        if sha_hex(&raw) != digest {
            return Err(format!("{path}: snapshot content hash"));
        }
        Ok(raw)
    }
}

// ───────────────────────────────── base trees ─────────────────────────────────

pub struct Base {
    pub commit: String,
    pub basis: String,
    pub archive: PathBuf,
    pub archive_sha256: String,
    pub archive_bytes: u64,
    pub freezes: Vec<String>,
}

fn load_bases(replay_dir: &Path) -> Res<Vec<Base>> {
    let doc = read_json(&replay_dir.join("BASES.json"))?;
    if doc.get("schema").and_then(Value::as_str) != Some(BASES_SCHEMA) {
        return Err("BASES.json: unexpected schema".into());
    }
    let mut out = Vec::new();
    for b in doc
        .get("bases")
        .and_then(Value::as_array)
        .ok_or("BASES.json: no bases")?
    {
        let ctx = "BASES.json base";
        out.push(Base {
            commit: str_field(b, "commit", ctx)?.to_string(),
            basis: str_field(b, "basis", ctx)?.to_string(),
            archive: replay_dir.join(str_field(b, "archive", ctx)?),
            archive_sha256: str_field(b, "archive_sha256", ctx)?.to_string(),
            archive_bytes: u64_field(b, "archive_bytes", ctx)?,
            freezes: b
                .get("freezes")
                .and_then(Value::as_array)
                .ok_or("BASES.json: base without freezes")?
                .iter()
                .filter_map(|f| f.as_str().map(str::to_string))
                .collect(),
        });
    }
    Ok(out)
}

/// The one base every given freeze agrees on. A freeze without a base, or two
/// freezes with different bases, fail closed.
fn base_for<'a>(bases: &'a [Base], freeze_rels: &[String]) -> Res<&'a Base> {
    let mut chosen: Option<&Base> = None;
    for rel in freeze_rels {
        let hit = bases
            .iter()
            .find(|b| b.freezes.iter().any(|f| f == rel))
            .ok_or_else(|| {
                format!(
                    "{rel}: no base tree in BASES.json (run `frozen-materialize snapshot-base`)"
                )
            })?;
        match chosen {
            Some(prev) if prev.commit != hit.commit => {
                return Err(format!(
                    "freezes disagree on the base tree: {} vs {}",
                    prev.commit, hit.commit
                ))
            }
            _ => chosen = Some(hit),
        }
    }
    chosen.ok_or_else(|| "no freeze given".to_string())
}

fn is_base_path(path: &str) -> bool {
    BASE_PATHS
        .iter()
        .any(|p| path == *p || path.starts_with(&format!("{p}/")))
}

/// Replace the unpinned inputs of `dest` by the base archive. Returns the
/// number of files restored.
fn restore_base(base: &Base, dest: &Path) -> Res<usize> {
    let packed = read(&base.archive)?;
    if packed.len() as u64 != base.archive_bytes {
        return Err(format!("{}: archive length", base.archive.display()));
    }
    if sha_hex(&packed) != base.archive_sha256 {
        return Err(format!(
            "{}: archive hash differs from BASES.json",
            base.archive.display()
        ));
    }
    let raw = gunzip_capped(&packed, MAX_BASE_BYTES, "base archive")?;
    // Wipe first, so a live file the base commit did not have cannot remain.
    for path in BASE_PATHS {
        let target = dest.join(path);
        match fs::symlink_metadata(&target) {
            Ok(m) if m.is_dir() => fs::remove_dir_all(&target),
            Ok(_) => fs::remove_file(&target),
            Err(_) => Ok(()),
        }
        .map_err(|e| format!("{}: {e}", target.display()))?;
    }
    let mut restored = 0usize;
    let mut archive = tar::Archive::new(raw.as_slice());
    for entry in archive.entries().map_err(|e| e.to_string())? {
        let mut entry = entry.map_err(|e| e.to_string())?;
        let kind = entry.header().entry_type();
        if kind.is_pax_global_extensions() || kind.is_pax_local_extensions() {
            continue;
        }
        let path = entry
            .path()
            .map_err(|e| e.to_string())?
            .to_string_lossy()
            .trim_end_matches('/')
            .to_string();
        if path.is_empty() {
            continue;
        }
        // `git archive` lists the parent directories of a file path too
        // (`docs/`, `docs/curves/`); those hold nothing and are allowed.
        let ancestor_dir = kind.is_dir()
            && BASE_PATHS
                .iter()
                .any(|p| p.starts_with(&format!("{path}/")));
        if !is_base_path(&path) && !ancestor_dir {
            return Err(format!(
                "base archive member outside the base paths: {path}"
            ));
        }
        if kind.is_file() {
            // `unpack_in` refuses `..` and absolute members.
            if !entry.unpack_in(dest).map_err(|e| format!("{path}: {e}"))? {
                return Err(format!("base archive member refused: {path}"));
            }
            restored += 1;
        } else if kind.is_dir() {
            fs::create_dir_all(dest.join(&path)).map_err(|e| format!("{path}: {e}"))?;
        } else {
            return Err(format!(
                "base archive member {path} is not a file or directory"
            ));
        }
    }
    Ok(restored)
}

// ─────────────────────────────── materialization ───────────────────────────────

fn copy_tree(src: &Path, dst: &Path) -> Res<()> {
    fs::create_dir_all(dst).map_err(|e| format!("{}: {e}", dst.display()))?;
    for item in fs::read_dir(src).map_err(|e| format!("{}: {e}", src.display()))? {
        let item = item.map_err(|e| e.to_string())?;
        let name = item.file_name();
        if name == ".git" || name == "target" {
            continue;
        }
        let from = item.path();
        let to = dst.join(&name);
        let kind = item.file_type().map_err(|e| e.to_string())?;
        if kind.is_symlink() {
            let link = fs::read_link(&from).map_err(|e| e.to_string())?;
            std::os::unix::fs::symlink(link, &to).map_err(|e| format!("{}: {e}", to.display()))?;
        } else if kind.is_dir() {
            copy_tree(&from, &to)?;
        } else {
            fs::copy(&from, &to).map_err(|e| format!("{}: {e}", from.display()))?;
        }
    }
    Ok(())
}

/// Build the frozen tree for `freezes` (absolute paths under `root`) in
/// `dest`, which must not exist. `replay_dir` holds `SNAPSHOTS.json`,
/// `BASES.json`, `historical_snapshots/` and the base archives.
pub fn materialize(root: &Path, freezes: &[PathBuf], dest: &Path, replay_dir: &Path) -> Res<Value> {
    if dest.exists() {
        return Err("never overwrite a materialized tree".into());
    }
    let pins = load_pins(freezes)?;
    let snapshots = Snapshots::load(replay_dir)?;
    let bases = load_bases(replay_dir)?;
    let rels: Vec<String> = freezes
        .iter()
        .map(|f| {
            f.strip_prefix(root)
                .map(|p| p.to_string_lossy().into_owned())
                .map_err(|_| format!("{}: not under {}", f.display(), root.display()))
        })
        .collect::<Res<_>>()?;
    let base = base_for(&bases, &rels)?;

    copy_tree(root, dest)?;
    let restored = restore_base(base, dest)?;

    let mut replayed = serde_json::Map::new();
    let mut live: Vec<&String> = Vec::new();
    let mut from_base: Vec<&String> = Vec::new();
    for (path, digest) in &pins {
        let target = dest.join(path);
        let current = fs::read(&target).ok().map(|b| sha_hex(&b));
        if current.as_deref() == Some(digest.as_str()) {
            if is_base_path(path) {
                from_base.push(path);
            } else {
                live.push(path);
            }
            continue;
        }
        if let Some(dir) = target.parent() {
            fs::create_dir_all(dir).map_err(|e| e.to_string())?;
        }
        fs::write(&target, snapshots.bytes(path, digest)?)
            .map_err(|e| format!("{}: {e}", target.display()))?;
        replayed.insert(
            path.clone(),
            json!({"frozen_sha256": digest, "live_sha256": current}),
        );
    }
    for (path, digest) in &pins {
        let got = sha_hex(&read(&dest.join(path))?);
        if got != *digest {
            return Err(format!(
                "{path}: materialized bytes hash to {got}, frozen {digest}"
            ));
        }
    }
    Ok(json!({
        "schema": RECEIPT_SCHEMA,
        "freezes": rels,
        "pinned_files": pins.len(),
        "live": live,
        "pinned_matching_base": from_base,
        "replayed_from_snapshot": replayed,
        "base": {
            "commit": base.commit,
            "basis": base.basis,
            "archive_sha256": base.archive_sha256,
            "files_restored": restored,
        },
    }))
}

// ───────────────────────────── base archive generation ─────────────────────────────

fn git(git_dir: Option<&Path>, args: &[&str]) -> Res<Vec<u8>> {
    let mut cmd = Command::new("git");
    if let Some(dir) = git_dir {
        cmd.arg("-C").arg(dir);
    }
    let out = cmd.args(args).output().map_err(|e| format!("git: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "git {}: {}",
            args.join(" "),
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }
    Ok(out.stdout)
}

/// Write the base archive of `commit` into `replay_dir/bases/` and record it
/// in `BASES.json` for `freeze_rel`. Developer tool: run from a full checkout,
/// commit the outputs. Deterministic for a given commit and `flate2` version.
pub fn snapshot_base(
    repo: &Path,
    replay_dir: &Path,
    freeze_rel: &str,
    commit: &str,
    basis: &str,
) -> Res<Value> {
    let full = String::from_utf8(git(
        Some(repo),
        &["rev-parse", "--verify", &format!("{commit}^{{commit}}")],
    )?)
    .map_err(|e| e.to_string())?
    .trim()
    .to_string();
    let existing: Vec<&str> = BASE_PATHS
        .iter()
        .copied()
        .filter(|p| git(Some(repo), &["cat-file", "-e", &format!("{full}:{p}")]).is_ok())
        .collect();
    if !existing.contains(&"src") {
        return Err(format!("{full}: no src/"));
    }
    let mut args = vec!["archive", "--format=tar", &full, "--"];
    args.extend(existing.iter());
    let tarball = git(Some(repo), &args)?;
    let mut enc = GzEncoder::new(Vec::new(), Compression::best());
    std::io::Write::write_all(&mut enc, &tarball).map_err(|e| e.to_string())?;
    let packed = enc.finish().map_err(|e| e.to_string())?;
    let rel = format!("bases/{}.tar.gz", &full[..12]);
    let file = replay_dir.join(&rel);
    fs::create_dir_all(replay_dir.join("bases")).map_err(|e| e.to_string())?;
    fs::write(&file, &packed).map_err(|e| format!("{}: {e}", file.display()))?;

    // Provenance a reviewer can re-derive with `git rev-parse <commit>:<path>`.
    let mut tree_oids = serde_json::Map::new();
    for p in &existing {
        let oid = String::from_utf8(git(Some(repo), &["rev-parse", &format!("{full}:{p}")])?)
            .map_err(|e| e.to_string())?;
        tree_oids.insert((*p).to_string(), json!(oid.trim()));
    }
    let entry = |freezes: Vec<Value>| {
        json!({
            "commit": full,
            "basis": basis,
            "archive": rel,
            "archive_sha256": sha_hex(&packed),
            "archive_bytes": packed.len(),
            "archived_paths": existing,
            "absent_at_commit": BASE_PATHS.iter().filter(|p| !existing.contains(p)).collect::<Vec<_>>(),
            "git_oids": tree_oids,
            "freezes": freezes,
        })
    };
    let bases_path = replay_dir.join("BASES.json");
    let mut doc = if bases_path.exists() {
        read_json(&bases_path)?
    } else {
        json!({"schema": BASES_SCHEMA, "bases": []})
    };
    let list = doc["bases"].as_array_mut().ok_or("BASES.json: no bases")?;
    // A freeze belongs to one base: drop it from any other entry first.
    for b in list.iter_mut() {
        if b["commit"] != json!(full) {
            if let Some(f) = b["freezes"].as_array_mut() {
                f.retain(|x| x.as_str() != Some(freeze_rel));
            }
        }
    }
    list.retain(|b| {
        b["freezes"].as_array().map_or(true, |f| !f.is_empty()) || b["commit"] == json!(full)
    });
    let mut freezes: BTreeSet<String> = BTreeSet::new();
    freezes.insert(freeze_rel.to_string());
    if let Some(pos) = list.iter().position(|b| b["commit"] == json!(full)) {
        for f in list[pos]["freezes"].as_array().into_iter().flatten() {
            if let Some(s) = f.as_str() {
                freezes.insert(s.to_string());
            }
        }
        list[pos] = entry(freezes.iter().map(|s| json!(s)).collect());
    } else {
        list.push(entry(freezes.iter().map(|s| json!(s)).collect()));
    }
    list.sort_by(|a, b| a["commit"].as_str().cmp(&b["commit"].as_str()));
    let text = serde_json::to_string_pretty(&doc).map_err(|e| e.to_string())? + "\n";
    fs::write(&bases_path, text).map_err(|e| e.to_string())?;
    Ok(entry(freezes.iter().map(|s| json!(s)).collect()))
}

#[cfg(test)]
mod tests;
