//! # The artifact store
//!
//! Factor bases, relation tables, DP archives and other files too large for
//! Git live in object storage and are fetched when something needs them.
//! Code names a file, not a location:
//!
//! ```text
//! ic run --factor-base store:factor-bases/icv1-f2m41-.../fb1-....json
//! ```
//!
//! and [`resolve_path`] turns that name into a verified local path,
//! downloading it into a content-addressed cache on first use.
//!
//! The pieces:
//!
//! * **Objects** are immutable and keyed by their SHA-256
//!   (`sha256/ab/abcdef…`). A cached copy therefore never goes stale, a
//!   remote object is never overwritten with different bytes, and every
//!   download is checked against the digest the index promised.
//! * **The index** (`store/index.json`, tracked in Git) maps a logical name
//!   to digest, size, description and provenance: the manifest AGENTS.md
//!   asks for when output is too large for Git. A name that already holds
//!   other content is not repointed without `--replace`, and the old digest
//!   is kept as superseded.
//! * **The cache** (`$CRYPTO_STORE_CACHE`, default `~/.cache/crypto-store`)
//!   holds verified objects, read-only, written atomically.
//! * **The remote** is `s3://bucket/prefix` (through the `aws` CLI, so
//!   credentials, SSO and multipart transfers are the CLI's) or a
//!   `file:///directory`, set by `CRYPTO_STORE_REMOTE` or
//!   `store/config.json`.
//!
//! The store prices nothing and makes no performance claim. It moves bytes
//! and checks them; the evidence rules in AGENTS.md apply to what is built
//! on the bytes. The remote is never deleted from: `rm` edits the index
//! only, and `gc` prunes the local cache only.

pub mod backend;
pub mod cli;
pub mod index;

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use serde::Deserialize;

use crate::hash::sha256::Sha256;
use backend::{object_key, Backend, DirBackend, Meta, RemoteSpec, S3Backend};
use index::{validate_name, validate_sha256, Entry, Index};

/// Where the store keeps things, resolved from the environment, the
/// optional `store/config.json`, and the repository layout.
#[derive(Clone, Debug)]
pub struct Config {
    pub index_path: PathBuf,
    pub cache_dir: PathBuf,
    pub remote: Option<RemoteSpec>,
    pub endpoint: Option<String>,
    pub aws: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FileConfig {
    remote: Option<String>,
    endpoint: Option<String>,
}

impl Config {
    /// Configuration for the repository containing the working directory.
    pub fn discover() -> Result<Self, String> {
        let cwd = std::env::current_dir().map_err(|e| e.to_string())?;
        Self::resolve(&cwd, &|k| std::env::var(k).ok().filter(|v| !v.is_empty()))
    }

    /// `env` is injected so the precedence rules are tested without
    /// touching the process environment.
    pub fn resolve(start: &Path, env: &dyn Fn(&str) -> Option<String>) -> Result<Self, String> {
        let root = find_root(start);
        let index_path = match env("CRYPTO_STORE_INDEX") {
            Some(p) => PathBuf::from(p),
            None => root.join("store/index.json"),
        };
        let file = root.join("store/config.json");
        let from_file: Option<FileConfig> = match fs::read(&file) {
            Ok(bytes) => Some(
                serde_json::from_slice(&bytes)
                    .map_err(|e| format!("{}: invalid store config: {e}", file.display()))?,
            ),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => None,
            Err(e) => return Err(format!("{}: {e}", file.display())),
        };
        let remote_text = env("CRYPTO_STORE_REMOTE")
            .or_else(|| from_file.as_ref().and_then(|f| f.remote.clone()));
        let remote = remote_text
            .map(|t| RemoteSpec::parse(&t))
            .transpose()
            .map_err(|e| format!("store remote: {e}"))?;
        let endpoint = env("CRYPTO_STORE_ENDPOINT").or_else(|| from_file.and_then(|f| f.endpoint));
        let cache_dir = match env("CRYPTO_STORE_CACHE") {
            Some(p) => PathBuf::from(p),
            None => match (env("XDG_CACHE_HOME"), env("HOME")) {
                (Some(x), _) => PathBuf::from(x).join("crypto-store"),
                (None, Some(h)) => PathBuf::from(h).join(".cache/crypto-store"),
                (None, None) => std::env::temp_dir().join("crypto-store-cache"),
            },
        };
        let aws = env("CRYPTO_STORE_AWS").unwrap_or_else(|| "aws".to_string());
        Ok(Self {
            index_path,
            cache_dir,
            remote,
            endpoint,
            aws,
        })
    }
}

/// The repository root: the nearest ancestor with a `.git` or a
/// `store/index.json`, else the starting directory.
fn find_root(start: &Path) -> PathBuf {
    for dir in start.ancestors() {
        if dir.join(".git").exists() || dir.join("store/index.json").is_file() {
            return dir.to_path_buf();
        }
    }
    start.to_path_buf()
}

/// SHA-256 and length of a file, streamed (the store never holds a whole
/// object in memory).
pub fn hash_file(path: &Path) -> Result<(String, u64), String> {
    let mut file = fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut total = 0u64;
    loop {
        let n = file
            .read(&mut buf)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
        total += n as u64;
    }
    Ok((hex::encode(hasher.finalize()), total))
}

/// What a caller gets back from [`Store::ensure`].
#[derive(Clone, Debug)]
pub struct Resolved {
    pub path: PathBuf,
    pub sha256: String,
    pub bytes: u64,
    /// True when this call downloaded the object.
    pub fetched: bool,
    pub name: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct PushMeta {
    pub description: String,
    pub tags: Vec<String>,
    pub provenance: String,
    pub replace: bool,
    /// Do not copy the pushed file into the local cache.
    pub no_cache: bool,
}

#[derive(Clone, Debug)]
pub struct PushReport {
    pub name: String,
    pub sha256: String,
    pub bytes: u64,
    /// False when the remote already held this exact object.
    pub uploaded: bool,
    pub cached: bool,
    pub replaced: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CacheState {
    Ok,
    Missing,
    Corrupt(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RemoteState {
    Ok,
    Absent,
    Mismatch(String),
}

#[derive(Clone, Debug)]
pub struct VerifyRow {
    pub name: String,
    pub sha256: String,
    pub bytes: u64,
    pub cache: CacheState,
    pub remote: Option<RemoteState>,
}

pub struct Store {
    pub cfg: Config,
    pub index: Index,
    backend: Option<Box<dyn Backend>>,
}

fn make_backend(cfg: &Config) -> Option<Box<dyn Backend>> {
    cfg.remote.as_ref().map(|r| match r {
        RemoteSpec::Dir(p) => Box::new(DirBackend::new(p.clone())) as Box<dyn Backend>,
        RemoteSpec::S3 { bucket, prefix } => Box::new(S3Backend {
            bucket: bucket.clone(),
            prefix: prefix.clone(),
            endpoint: cfg.endpoint.clone(),
            aws: cfg.aws.clone(),
        }),
    })
}

fn unique_suffix() -> String {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.subsec_nanos())
        .unwrap_or(0);
    format!("{}.{}", std::process::id(), nanos)
}

#[cfg(unix)]
fn make_read_only(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    let _ = fs::set_permissions(path, fs::Permissions::from_mode(0o444));
}
#[cfg(not(unix))]
fn make_read_only(_path: &Path) {}

impl Store {
    pub fn open(cfg: Config) -> Result<Self, String> {
        let index = Index::load(&cfg.index_path)?;
        let backend = make_backend(&cfg);
        Ok(Self {
            cfg,
            index,
            backend,
        })
    }

    pub fn open_default() -> Result<Self, String> {
        Self::open(Config::discover()?)
    }

    pub fn cache_path(&self, sha256: &str) -> PathBuf {
        self.cfg
            .cache_dir
            .join("sha256")
            .join(&sha256[..2])
            .join(sha256)
    }

    pub fn remote(&self) -> Result<&dyn Backend, String> {
        self.backend.as_deref().ok_or_else(|| {
            "no remote configured: set CRYPTO_STORE_REMOTE (s3://bucket/prefix or \
             file:///directory) or put {\"remote\": \"...\"} in store/config.json"
                .to_string()
        })
    }

    pub fn save_index(&self) -> Result<(), String> {
        self.index.save(&self.cfg.index_path)
    }

    /// Turn `store:NAME`, `NAME` or `sha256:DIGEST` into
    /// `(name, digest, size when the index knows it)`.
    pub fn lookup(&self, spec: &str) -> Result<(Option<String>, String, Option<u64>), String> {
        let spec = spec.strip_prefix("store:").unwrap_or(spec);
        if let Some(digest) = spec.strip_prefix("sha256:") {
            if digest.len() == 64 {
                validate_sha256(digest)?;
                let known = self.index.entries.iter().find(|(_, e)| e.sha256 == digest);
                return Ok((
                    known.map(|(n, _)| n.clone()),
                    digest.to_string(),
                    known.map(|(_, e)| e.bytes),
                ));
            }
            let hits = self.index.by_sha256(digest);
            return match hits.as_slice() {
                [(name, e)] => Ok((Some((*name).clone()), e.sha256.clone(), Some(e.bytes))),
                [] => Err(format!(
                    "no indexed entry has a digest starting {digest:?} (a prefix needs 8+ \
                     digits and an index entry; use the full 64 for an unindexed object)"
                )),
                many => Err(format!(
                    "digest prefix {digest:?} is ambiguous: {}",
                    many.iter()
                        .map(|(n, _)| n.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                )),
            };
        }
        validate_name(spec)?;
        match self.index.get(spec) {
            Some(e) => Ok((Some(spec.to_string()), e.sha256.clone(), Some(e.bytes))),
            None => {
                let near: Vec<&str> = self
                    .index
                    .entries
                    .keys()
                    .filter(|n| n.contains(spec) || spec.contains(n.as_str()))
                    .take(5)
                    .map(String::as_str)
                    .collect();
                Err(format!(
                    "{spec:?} is not in {}{}",
                    self.cfg.index_path.display(),
                    if near.is_empty() {
                        String::new()
                    } else {
                        format!(" (similar: {})", near.join(", "))
                    }
                ))
            }
        }
    }

    /// A verified local path for a named object, downloading it when it is
    /// not cached. A cached object is trusted by size only (it was hashed
    /// when it entered the cache and is read-only); `verify` re-hashes.
    pub fn ensure(&self, spec: &str) -> Result<Resolved, String> {
        let (name, sha256, expected) = self.lookup(spec)?;
        let path = self.cache_path(&sha256);
        if let Ok(meta) = fs::metadata(&path) {
            if expected.is_none_or(|b| b == meta.len()) {
                return Ok(Resolved {
                    path,
                    sha256,
                    bytes: meta.len(),
                    fetched: false,
                    name,
                });
            }
            // Wrong size: the cache entry is damaged. Replace it.
            let _ = fs::remove_file(&path);
        }
        let remote = self.remote().map_err(|e| {
            format!(
                "{} is not in the local cache ({}) and cannot be fetched: {e}",
                name.as_deref().unwrap_or(&sha256),
                self.cfg.cache_dir.display()
            )
        })?;
        let key = object_key(&sha256)?;
        let meta = remote
            .head(&key)?
            .ok_or_else(|| format!("{key} is absent from {}", remote.describe()))?;
        if let Some(b) = expected {
            if meta.bytes != b {
                return Err(format!(
                    "{key}: the remote object is {} bytes, the index says {b}",
                    meta.bytes
                ));
            }
        }
        let tmp = self
            .cfg
            .cache_dir
            .join("tmp")
            .join(format!("{sha256}.{}", unique_suffix()));
        remote.get(&key, &tmp)?;
        let (got, bytes) = match hash_file(&tmp) {
            Ok(v) => v,
            Err(e) => {
                let _ = fs::remove_file(&tmp);
                return Err(e);
            }
        };
        if got != sha256 {
            let _ = fs::remove_file(&tmp);
            return Err(format!(
                "{key}: downloaded bytes hash to {got}, expected {sha256}; discarded"
            ));
        }
        make_read_only(&tmp);
        if let Some(dir) = path.parent() {
            fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        }
        fs::rename(&tmp, &path).map_err(|e| {
            let _ = fs::remove_file(&tmp);
            format!("{}: {e}", path.display())
        })?;
        Ok(Resolved {
            path,
            sha256,
            bytes,
            fetched: true,
            name,
        })
    }

    /// Upload a file (unless the remote already holds it), cache it, and
    /// record it in the index under `name`.
    pub fn push(&mut self, name: &str, file: &Path, meta: PushMeta) -> Result<PushReport, String> {
        validate_name(name)?;
        let m = fs::metadata(file).map_err(|e| format!("{}: {e}", file.display()))?;
        if !m.is_file() {
            return Err(format!("{} is not a regular file", file.display()));
        }
        let (sha256, bytes) = hash_file(file)?;
        let entry = Entry {
            sha256: sha256.clone(),
            bytes,
            description: meta.description.clone(),
            tags: meta.tags.clone(),
            provenance: meta.provenance.clone(),
            superseded: vec![],
        };
        // Settle the index first: a name conflict must fail before any upload.
        let replaced = self.index.get(name).is_some_and(|old| old.sha256 != sha256);
        let mut next = self.index.clone();
        next.put(name, entry, meta.replace)?;

        let remote = self.remote()?;
        let key = object_key(&sha256)?;
        let present: Option<Meta> = remote.head(&key)?;
        let uploaded = match present {
            Some(m) if m.bytes == bytes => false,
            _ => {
                remote.put(&key, file, &sha256)?;
                let after = remote
                    .head(&key)?
                    .ok_or_else(|| format!("{key}: absent from the remote right after upload"))?;
                if after.bytes != bytes {
                    return Err(format!(
                        "{key}: the remote holds {} bytes after upload, expected {bytes}",
                        after.bytes
                    ));
                }
                if let Some(s) = &after.sha256 {
                    if *s != sha256 {
                        return Err(format!(
                            "{key}: the remote recorded digest {s}, expected {sha256}"
                        ));
                    }
                }
                true
            }
        };
        let mut cached = false;
        if !meta.no_cache {
            let path = self.cache_path(&sha256);
            if !path.is_file() {
                backend::copy_into(file, &path)?;
                make_read_only(&path);
            }
            cached = true;
        }
        self.index = next;
        self.save_index()?;
        Ok(PushReport {
            name: name.to_string(),
            sha256,
            bytes,
            uploaded,
            cached,
            replaced,
        })
    }

    /// Check indexed objects: re-hash the cached copy and, with `remote`,
    /// confirm the remote object exists with the indexed size.
    pub fn verify(&self, names: &[String], remote: bool) -> Result<Vec<VerifyRow>, String> {
        let chosen: Vec<(&String, &Entry)> = if names.is_empty() {
            self.index.entries.iter().collect()
        } else {
            let mut v = Vec::new();
            for n in names {
                let (name, _, _) = self.lookup(n)?;
                let name = name.ok_or_else(|| format!("{n}: not an indexed name"))?;
                let (k, e) = self.index.entries.get_key_value(&name).expect("looked up");
                v.push((k, e));
            }
            v
        };
        let mut rows = Vec::new();
        for (name, e) in chosen {
            let path = self.cache_path(&e.sha256);
            let cache = if !path.exists() {
                CacheState::Missing
            } else {
                match hash_file(&path) {
                    Ok((h, b)) if h == e.sha256 && b == e.bytes => CacheState::Ok,
                    Ok((h, b)) => CacheState::Corrupt(format!("hashes to {h} ({b} bytes)")),
                    Err(err) => CacheState::Corrupt(err),
                }
            };
            let remote_state = if remote {
                let be = self.remote()?;
                Some(match be.head(&object_key(&e.sha256)?)? {
                    None => RemoteState::Absent,
                    Some(m) if m.bytes != e.bytes => {
                        RemoteState::Mismatch(format!("{} bytes, index says {}", m.bytes, e.bytes))
                    }
                    Some(m) => match m.sha256 {
                        Some(s) if s != e.sha256 => {
                            RemoteState::Mismatch(format!("recorded digest {s}"))
                        }
                        _ => RemoteState::Ok,
                    },
                })
            } else {
                None
            };
            rows.push(VerifyRow {
                name: name.clone(),
                sha256: e.sha256.clone(),
                bytes: e.bytes,
                cache,
                remote: remote_state,
            });
        }
        Ok(rows)
    }

    /// Cached objects the index no longer names, and leftover partial
    /// downloads. Removes them when `apply`. The remote is never touched.
    pub fn gc(&self, apply: bool) -> Result<Vec<(PathBuf, u64)>, String> {
        let keep: std::collections::BTreeSet<&str> = self
            .index
            .entries
            .values()
            .map(|e| e.sha256.as_str())
            .collect();
        let mut doomed = Vec::new();
        let objects = self.cfg.cache_dir.join("sha256");
        if let Ok(shards) = fs::read_dir(&objects) {
            for shard in shards.flatten() {
                if let Ok(items) = fs::read_dir(shard.path()) {
                    for item in items.flatten() {
                        let name = item.file_name().to_string_lossy().into_owned();
                        if !keep.contains(name.as_str()) {
                            let len = item.metadata().map(|m| m.len()).unwrap_or(0);
                            doomed.push((item.path(), len));
                        }
                    }
                }
            }
        }
        if let Ok(items) = fs::read_dir(self.cfg.cache_dir.join("tmp")) {
            for item in items.flatten() {
                let len = item.metadata().map(|m| m.len()).unwrap_or(0);
                doomed.push((item.path(), len));
            }
        }
        doomed.sort();
        if apply {
            for (p, _) in &doomed {
                fs::remove_file(p).map_err(|e| format!("{}: {e}", p.display()))?;
            }
        }
        Ok(doomed)
    }
}

/// True for the spellings [`resolve_path`] treats as store references.
pub fn is_store_ref(s: &str) -> bool {
    s.starts_with("store:") || s.starts_with("sha256:")
}

/// The one call code makes to accept a store reference where it accepted a
/// path: a plain path is returned unchanged; `store:NAME` or `sha256:DIGEST`
/// becomes a verified local path, fetched on first use.
pub fn resolve_path(arg: &Path) -> Result<PathBuf, String> {
    let text = match arg.to_str() {
        Some(t) if is_store_ref(t) => t,
        _ => return Ok(arg.to_path_buf()),
    };
    let store = Store::open_default()?;
    let got = store.ensure(text)?;
    if got.fetched {
        eprintln!(
            "store: fetched {} ({} bytes, sha256 {})",
            got.name.as_deref().unwrap_or(text),
            got.bytes,
            &got.sha256[..16]
        );
    }
    Ok(got.path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU32, Ordering};

    static N: AtomicU32 = AtomicU32::new(0);

    /// A scratch directory and a store whose remote and cache are inside it.
    fn scratch() -> (PathBuf, Store) {
        let root = std::env::temp_dir().join(format!(
            "store-test-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::SeqCst)
        ));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let cfg = Config {
            index_path: root.join("repo/store/index.json"),
            cache_dir: root.join("cache"),
            remote: Some(RemoteSpec::Dir(root.join("remote"))),
            endpoint: None,
            aws: "aws".into(),
        };
        (root, Store::open(cfg).unwrap())
    }

    fn file(root: &Path, name: &str, bytes: &[u8]) -> PathBuf {
        let p = root.join(name);
        fs::write(&p, bytes).unwrap();
        p
    }

    #[test]
    fn hashing_matches_the_known_digest_across_buffer_sizes() {
        let root = std::env::temp_dir().join(format!("store-hash-{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        let abc = file(&root, "abc", b"abc");
        assert_eq!(
            hash_file(&abc).unwrap(),
            (
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad".to_string(),
                3
            )
        );
        // Longer than the 1 MiB read buffer, to cross a refill.
        let big: Vec<u8> = (0..(3 << 19) + 7).map(|i| (i % 251) as u8).collect();
        let p = file(&root, "big", &big);
        assert_eq!(
            hash_file(&p).unwrap(),
            (
                hex::encode(crate::hash::sha256::sha256(&big)),
                big.len() as u64
            )
        );
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn push_then_pull_on_a_clean_cache_round_trips() {
        let (root, mut s) = scratch();
        let src = file(&root, "fb.bin", b"factor base bytes");
        let r = s
            .push(
                "factor-bases/demo",
                &src,
                PushMeta {
                    description: "demo".into(),
                    tags: vec!["fb".into()],
                    provenance: "test".into(),
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(r.uploaded && r.cached && !r.replaced);

        // A second machine: same index and remote, empty cache.
        let mut cfg = s.cfg.clone();
        cfg.cache_dir = root.join("cache2");
        let t = Store::open(cfg).unwrap();
        let got = t.ensure("store:factor-bases/demo").unwrap();
        assert!(got.fetched);
        assert_eq!(fs::read(&got.path).unwrap(), b"factor base bytes");
        let again = t.ensure("factor-bases/demo").unwrap();
        assert!(!again.fetched);
        assert_eq!(again.path, got.path);
        assert!(t
            .ensure(&format!("sha256:{}", &r.sha256[..12]))
            .unwrap()
            .path
            .ends_with(&r.sha256));
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn pushing_twice_uploads_once_and_a_conflict_uploads_nothing() {
        let (root, mut s) = scratch();
        let a = file(&root, "a", b"one");
        let b = file(&root, "b", b"two");
        assert!(s.push("n", &a, PushMeta::default()).unwrap().uploaded);
        assert!(!s.push("n", &a, PushMeta::default()).unwrap().uploaded);
        let objects_before = s.remote().unwrap().list("").unwrap().len();
        let err = s.push("n", &b, PushMeta::default()).unwrap_err();
        assert!(err.contains("--replace"), "{err}");
        assert_eq!(s.remote().unwrap().list("").unwrap().len(), objects_before);
        let r = s
            .push(
                "n",
                &b,
                PushMeta {
                    replace: true,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(r.replaced && r.uploaded);
        assert_eq!(s.index.get("n").unwrap().superseded.len(), 1);
        // The old object is still there: replacing never deletes.
        assert_eq!(
            s.remote().unwrap().list("").unwrap().len(),
            objects_before + 1
        );
        // The index on disk is what a reader sees.
        assert_eq!(Index::load(&s.cfg.index_path).unwrap(), s.index);
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn a_corrupted_remote_object_is_refused_and_not_cached() {
        let (root, mut s) = scratch();
        let src = file(&root, "x", b"original");
        let r = s.push("x", &src, PushMeta::default()).unwrap();
        // Corrupt the remote object keeping its length.
        let key = object_key(&r.sha256).unwrap();
        fs::remove_file(root.join("remote").join(&key)).unwrap();
        fs::write(root.join("remote").join(&key), b"tampered").unwrap();
        let mut cfg = s.cfg.clone();
        cfg.cache_dir = root.join("fresh");
        let t = Store::open(cfg).unwrap();
        let err = t.ensure("x").unwrap_err();
        assert!(err.contains("expected"), "{err}");
        assert!(!t.cache_path(&r.sha256).exists());
        assert!(t.gc(false).unwrap().is_empty(), "no partial file left over");
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn a_damaged_cache_entry_is_replaced_and_verify_reports_it() {
        let (root, mut s) = scratch();
        let src = file(&root, "y", b"payload-1234");
        let r = s.push("y", &src, PushMeta::default()).unwrap();
        let cached = s.cache_path(&r.sha256);
        // Same size, wrong bytes: only `verify` can see it.
        fs::remove_file(&cached).unwrap();
        fs::write(&cached, b"payload-9999").unwrap();
        let rows = s.verify(&[], true).unwrap();
        assert!(matches!(rows[0].cache, CacheState::Corrupt(_)));
        assert_eq!(rows[0].remote, Some(RemoteState::Ok));
        // Wrong size: `ensure` notices and refetches.
        fs::remove_file(&cached).unwrap();
        fs::write(&cached, b"short").unwrap();
        let got = s.ensure("y").unwrap();
        assert!(got.fetched);
        assert_eq!(fs::read(&got.path).unwrap(), b"payload-1234");
        assert_eq!(s.verify(&[], false).unwrap()[0].cache, CacheState::Ok);
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn no_remote_and_no_cache_says_exactly_what_to_set() {
        let (root, mut s) = scratch();
        let src = file(&root, "z", b"zzz");
        s.push("z", &src, PushMeta::default()).unwrap();
        let mut cfg = s.cfg.clone();
        cfg.cache_dir = root.join("empty-cache");
        cfg.remote = None;
        let mut t = Store::open(cfg).unwrap();
        let err = t.ensure("z").unwrap_err();
        assert!(err.contains("CRYPTO_STORE_REMOTE"), "{err}");
        assert!(t.push("z2", &src, PushMeta::default()).is_err());
        let err = t.ensure("nope").unwrap_err();
        assert!(err.contains("not in"), "{err}");
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn gc_prunes_unindexed_cache_objects_only() {
        let (root, mut s) = scratch();
        let a = file(&root, "a", b"keep me");
        let b = file(&root, "b", b"drop me");
        s.push("keep", &a, PushMeta::default()).unwrap();
        let rb = s.push("drop", &b, PushMeta::default()).unwrap();
        s.index.entries.remove("drop");
        s.save_index().unwrap();
        let listed = s.gc(false).unwrap();
        assert_eq!(listed.len(), 1);
        assert!(listed[0].0.ends_with(&rb.sha256));
        assert!(s.cache_path(&rb.sha256).exists(), "dry run removes nothing");
        s.gc(true).unwrap();
        assert!(!s.cache_path(&rb.sha256).exists());
        assert!(s.ensure("keep").is_ok());
        // The remote still has both objects.
        assert_eq!(s.remote().unwrap().list("").unwrap().len(), 2);
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn configuration_precedence_is_env_then_file_then_default() {
        let root = std::env::temp_dir().join(format!("store-cfg-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("repo/.git")).unwrap();
        fs::create_dir_all(root.join("repo/store")).unwrap();
        fs::write(
            root.join("repo/store/config.json"),
            r#"{"remote":"s3://from-file/p","endpoint":"http://file:9000"}"#,
        )
        .unwrap();
        let sub = root.join("repo/deep/er");
        fs::create_dir_all(&sub).unwrap();
        let none = |_: &str| None;
        let c = Config::resolve(&sub, &none).unwrap();
        assert_eq!(c.index_path, root.join("repo/store/index.json"));
        assert_eq!(
            c.remote,
            Some(RemoteSpec::S3 {
                bucket: "from-file".into(),
                prefix: "p".into()
            })
        );
        assert_eq!(c.endpoint.as_deref(), Some("http://file:9000"));
        assert_eq!(c.aws, "aws");
        let env = |k: &str| match k {
            "CRYPTO_STORE_REMOTE" => Some("file:///srv/x".to_string()),
            "CRYPTO_STORE_CACHE" => Some("/c".to_string()),
            "CRYPTO_STORE_AWS" => Some("/bin/myaws".to_string()),
            _ => None,
        };
        let c = Config::resolve(&sub, &env).unwrap();
        assert_eq!(c.remote, Some(RemoteSpec::Dir("/srv/x".into())));
        assert_eq!(c.cache_dir, PathBuf::from("/c"));
        assert_eq!(c.aws, "/bin/myaws");
        let bad = |k: &str| (k == "CRYPTO_STORE_REMOTE").then(|| "nonsense".to_string());
        assert!(Config::resolve(&sub, &bad).is_err());
        fs::write(root.join("repo/store/config.json"), r#"{"remote":1}"#).unwrap();
        assert!(Config::resolve(&sub, &none).is_err());
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn plain_paths_pass_through_resolve_unchanged() {
        assert_eq!(
            resolve_path(Path::new("some/file.json")).unwrap(),
            PathBuf::from("some/file.json")
        );
        assert!(is_store_ref("store:a") && is_store_ref("sha256:ab"));
        assert!(!is_store_ref("./store:a") && !is_store_ref("a"));
    }
}
