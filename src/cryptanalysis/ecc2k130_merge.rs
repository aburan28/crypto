//! The ECC2K-130 distinguished-point merge: a port of
//! `ecc2k130/aws/merge.py` that must agree with it byte for byte.
//!
//! Workers' corpora are bucketed by a mixed hash of the orbit key; a pass
//! sorts every bucket that received records, drops exact duplicates, and
//! reports adjacent equal keys with distinct seeds as collisions, which the
//! host client then re-walks and solves.  Everything the Python writes — the
//! bucket files, `state.json`, `campaign.lock.json`, `solution.json` and the
//! JSON summary on stdout — is reproduced exactly, which is what lets the
//! certification suite and `test_merge_parity.py` hold the two side by side.
//!
//! The parts of Python that decide those bytes are reproduced rather than
//! approximated: `json.dumps` (ASCII escaping, key order, indentation),
//! dict insertion order and mutation, `np.lexsort`'s order, and the
//! campaign-id hash of `protocol.campaignContract`.

use std::collections::BTreeSet;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::os::unix::fs::OpenOptionsExt;
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use num_bigint::BigUint;
use rayon::prelude::*;

use super::ecc2k130_pyjson::{dumps, parse_json, py_splitlines, py_strip, Json, Obj, Style};
use crate::hash::sha256::{sha256, Sha256};

// ── Storage protocol (`aws/protocol.py`) ───────────────────────────

pub const PROTOCOL: &str = "ecc2k-seed-orbit-v1";
pub const RECORD_BYTES: u64 = 32;
const RECORD_V2_BYTES: u64 = 72;
const DP_MAGIC_V2: &[u8; 8] = b"ECC2KDP2";
const DP_MAGIC_TABLE3: &[u8; 8] = b"ECC2KDT3";
const DP_HEADER_BYTES: u64 = 16;

/// `protocol.WALKS`; part of the campaign id, so it must stay byte-equal.
fn walk_description(walk: &str) -> Option<&'static str> {
    match walk {
        "sigma" => Some("sigma^(3+((normal-weight(x)>>1)&7))(R)+R"),
        "table" => Some(concat!(
            "R+(-1)^eps(R)*sigma^k(R)(T[(normal-weight(x)>>1)&7]);",
            "k=frobenius-phase(x);eps=pivot-coordinate(y);",
            "cycle-rule=v3-raw-cycle8-canonical-eligible-anchor"
        )),
        _ => None,
    }
}

fn hex_sha256(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

/// `protocol.sha256File`.
pub fn sha256_file(path: &Path) -> io::Result<String> {
    let mut f = File::open(path)?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(hex::encode(h.finalize()))
}

/// `protocol.syncDirectory`.
fn sync_dir(path: &Path) -> io::Result<()> {
    File::open(if path.as_os_str().is_empty() {
        Path::new(".")
    } else {
        path
    })?
    .sync_all()
}

fn absolute(path: &Path) -> io::Result<PathBuf> {
    if path.is_absolute() {
        Ok(path.to_path_buf())
    } else {
        Ok(std::env::current_dir()?.join(path))
    }
}

/// `protocol.atomicJson`: indent 2, sorted keys, newline, fsync, rename,
/// fsync the directory; the file is mode 0600, as `mkstemp` makes it.
pub fn atomic_json(path: &Path, v: &Json) -> Result<(), String> {
    let body = dumps(v, Style::Indent(2), true)? + "\n";
    let parent = absolute(path)
        .map_err(|e| e.to_string())?
        .parent()
        .map(Path::to_path_buf)
        .ok_or("no parent directory")?;
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let tmp = parent.join(format!(".commit-{}-{stamp}", std::process::id()));
    let written = (|| -> io::Result<()> {
        let mut f = OpenOptions::new()
            .write(true)
            .create_new(true)
            .mode(0o600)
            .open(&tmp)?;
        f.write_all(body.as_bytes())?;
        f.sync_all()?;
        drop(f);
        fs::rename(&tmp, path)?;
        sync_dir(&parent)
    })();
    if tmp.exists() {
        let _ = fs::remove_file(&tmp);
    }
    written.map_err(|e| format!("{}: {e}", path.display()))
}

/// What `protocol.campaignContract` returns: `{"id": ..., "contract": ...}`.
pub fn campaign_contract(config: &Obj) -> Result<Json, String> {
    if config.get("storageProtocol") != Some(&Json::str(PROTOCOL)) {
        return Err(format!(
            "strict storage requires storageProtocol={PROTOCOL}"
        ));
    }
    if config.get("extraArgs").is_some_and(Json::truthy) {
        return Err("strict campaigns forbid extraArgs overriding protocol parameters".into());
    }
    let mut c = Obj::new();
    for key in [
        "curve",
        "dpWeight",
        "maxIters",
        "packed",
        "workers",
        "batch",
        "binarySha256",
        "hostBinarySha256",
        "sourceSha256",
    ] {
        let v = config
            .get(key)
            .ok_or_else(|| format!("campaign config has no {key}"))?;
        c.set(key, v.clone());
    }
    let live = c.get("maxIters").cloned().unwrap_or(Json::Null);
    let guard = config
        .get("contractMaxIters")
        .cloned()
        .unwrap_or(live.clone());
    c.set("maxIters", guard.clone());
    let (Some(live), Some(guard)) = (live.as_int(), guard.as_int()) else {
        return Err("invalid maxIters".into());
    };
    if guard < 0 {
        return Err("invalid maxIters".into());
    }
    if guard == 0 && live != 0 {
        return Err(format!(
            "the campaign was created with no guard (maxIters 0); maxIters {live} would cut trails already under way"
        ));
    }
    if live < guard {
        return Err(format!(
            "maxIters {live} is below the campaign's guard {guard}; only a raise keeps trails already under way"
        ));
    }
    let curve = c.get("curve").and_then(Json::as_int);
    if !matches!(curve, Some(23 | 41 | 83 | 131)) {
        return Err("strict storage currently supports normal-basis curves 23/41/83/131".into());
    }
    let curve = curve.expect("checked above");
    for key in ["dpWeight", "maxIters", "workers", "batch"] {
        if c.get(key).and_then(Json::as_int).is_none_or(|v| v < 0) {
            return Err(format!("invalid {key}"));
        }
    }
    let int = |key: &str| c.get(key).and_then(Json::as_int).expect("checked above");
    if !(0..=curve).contains(&int("dpWeight")) || int("workers") == 0 || int("batch") == 0 {
        return Err("invalid cutoff or worker geometry".into());
    }
    match c.get("packed") {
        Some(Json::Bool(packed)) if !*packed || curve == 131 => {}
        _ => return Err("invalid packed backend".into()),
    }
    for key in ["binarySha256", "hostBinarySha256", "sourceSha256"] {
        let pinned = c.get(key).and_then(Json::as_str).is_some_and(|s| {
            s.len() == 64
                && s.bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        });
        if !pinned {
            return Err(format!("strict campaigns require a pinned {key}"));
        }
    }
    let walk = match config.get("walk") {
        None => "sigma".to_string(),
        Some(Json::Str(s)) => s.clone(),
        Some(other) => {
            return Err(format!(
                "unknown walk {}; one of ['sigma', 'table']",
                dumps(other, Style::Compact, false).unwrap_or_default()
            ))
        }
    };
    let described = walk_description(&walk)
        .ok_or_else(|| format!("unknown walk '{walk}'; one of ['sigma', 'table']"))?;
    c.set("protocol", Json::str(PROTOCOL));
    c.set("recordBytes", Json::Int(RECORD_BYTES as i128));
    c.set(
        "key",
        Json::str("min-normal-basis-x-over-frobenius;negation-quotient"),
    );
    c.set("walk", Json::str(described));
    c.set(
        "seed",
        Json::str("run16-walk32-counter16;splitmix64;128-frobenius-terms"),
    );
    c.set(
        "coefficients",
        Json::str("absent;recover-by-seed-replay;verify-kP-equals-Q"),
    );
    let raw = dumps(&Json::Obj(c.clone()), Style::Compact, true)?;
    Ok(Json::Obj(
        Obj::new()
            .with("id", Json::Str(hex_sha256(raw.as_bytes())))
            .with("contract", Json::Obj(c)),
    ))
}

fn campaign_id(campaign: &Json) -> &str {
    campaign
        .as_obj()
        .and_then(|o| o.get("id"))
        .and_then(Json::as_str)
        .unwrap_or("")
}

/// `protocol.verifyEnvelope(path, manifest, campaign, "dp")`.
fn verify_envelope(path: &Path, manifest: &Json, campaign: &Json) -> Result<(), String> {
    let size = fs::metadata(path)
        .map_err(|e| format!("{}: {e}", path.display()))?
        .len();
    if size % RECORD_BYTES != 0 {
        return Err("partial DP record".into());
    }
    let digest = sha256_file(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let expected = [
        ("protocol", Json::str(PROTOCOL)),
        ("campaignId", Json::str(campaign_id(campaign))),
        ("kind", Json::str("dp")),
        ("bytes", Json::Int(size as i128)),
        ("sha256", Json::Str(digest)),
        ("records", Json::Int((size / RECORD_BYTES) as i128)),
    ];
    let meta = manifest.as_obj();
    let mismatch = expected
        .iter()
        .any(|(k, v)| meta.and_then(|m| m.get(k)).unwrap_or(&Json::Null) != v);
    if mismatch {
        return Err(format!(
            "manifest mismatch, corruption or incompatible campaign: {}",
            path.display()
        ));
    }
    Ok(())
}

/// `protocol.bindDirectory`.
fn bind_directory(root: &Path, campaign: &Json) -> Result<(), String> {
    let path = root.join("campaign.lock.json");
    if path.exists() {
        if read_json(&path)? != *campaign {
            return Err("directory belongs to a different campaign".into());
        }
        return Ok(());
    }
    if ["walk.ck", "dp.bin", "state.json"]
        .iter()
        .any(|name| root.join(name).exists())
    {
        return Err("legacy data requires an explicit audited migration or new directory".into());
    }
    atomic_json(&path, campaign)
}

fn read_json(path: &Path) -> Result<Json, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse_json(&text).map_err(|e| format!("{}: {e}", path.display()))
}

// ── Corpus formats ─────────────────────────────────────────────────

/// `(header bytes, record bytes, witness layout?)` for a corpus file, by
/// its magic.  The witness (v2) layout is `seed, iters, k0..k2, counts[8]`.
fn corpus_format(path: &Path) -> io::Result<(u64, u64, bool)> {
    let mut magic = Vec::with_capacity(8);
    File::open(path)?.take(8).read_to_end(&mut magic)?;
    Ok(if magic == DP_MAGIC_TABLE3 {
        (DP_HEADER_BYTES, RECORD_BYTES, false)
    } else if magic == DP_MAGIC_V2 {
        (DP_HEADER_BYTES, RECORD_V2_BYTES, true)
    } else {
        (0, RECORD_BYTES, false)
    })
}

fn word(b: &[u8], at: usize) -> u64 {
    u64::from_le_bytes(b[at..at + 8].try_into().expect("8 bytes"))
}

/// `(seed, k0, k1, k2)` of one record in either layout.
fn key_record(rec: &[u8], v2: bool) -> [u64; 4] {
    if v2 {
        [word(rec, 0), word(rec, 16), word(rec, 24), word(rec, 32)]
    } else {
        [word(rec, 0), word(rec, 8), word(rec, 16), word(rec, 24)]
    }
}

fn record_bytes(r: &[u64; 4]) -> [u8; 32] {
    let mut out = [0u8; 32];
    for (i, w) in r.iter().enumerate() {
        out[i * 8..i * 8 + 8].copy_from_slice(&w.to_le_bytes());
    }
    out
}

/// The bucket of an orbit key: raw low bits of sparse canonical `x` keys
/// are not uniform, so all three words are mixed first.
pub fn bucket_of(k0: u64, k1: u64, k2: u64, buckets: u64) -> u64 {
    let mut h = k0.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    h ^= k1
        .wrapping_add(0x632B_E59B_D9B4_E019)
        .wrapping_mul(0xBF58_476D_1CE4_E5B9);
    h ^= k2
        .wrapping_add(0x94D0_49BB_1331_11EB)
        .wrapping_mul(0x2545_F491_4F6C_DD1D);
    (h ^ (h >> 29)) & buckets.wrapping_sub(1)
}

/// `*.bin` under `root`, relative paths sorted as Python sorts `os.walk`'s
/// joined paths.  Symlinked directories are listed but not entered.
fn source_files(root: &Path) -> Vec<String> {
    fn walk(dir: &Path, rel: &str, out: &mut Vec<String>) {
        let Ok(entries) = fs::read_dir(dir) else {
            return;
        };
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            let path = entry.path();
            let sub = if rel.is_empty() {
                name.clone()
            } else {
                format!("{rel}/{name}")
            };
            let Ok(ft) = entry.file_type() else { continue };
            let is_dir = if ft.is_symlink() {
                fs::metadata(&path).is_ok_and(|m| m.is_dir())
            } else {
                ft.is_dir()
            };
            if is_dir {
                if !ft.is_symlink() {
                    walk(&path, &sub, out);
                }
            } else if name.ends_with(".bin") {
                out.push(sub);
            }
        }
    }
    let mut out = Vec::new();
    walk(root, "", &mut out);
    out.sort();
    out
}

fn bucket_name(b: u64) -> String {
    format!("{b:04x}.bin")
}

/// 8 MiB of records at a time, as the Python reads them.
const READ_BYTES: u64 = 8 * 1024 * 1024;

/// `fh.read(n)` on a regular file: up to `n` bytes, short only at EOF.
fn read_up_to(f: &mut File, n: usize) -> io::Result<Vec<u8>> {
    let mut buf = Vec::with_capacity(n);
    (&mut *f).take(n as u64).read_to_end(&mut buf)?;
    Ok(buf)
}

fn obj_mut<'a>(state: &'a mut Obj, key: &str) -> Result<&'a mut Obj, String> {
    match state.get_mut(key) {
        Some(Json::Obj(o)) => Ok(o),
        _ => Err(format!("merge state has no {key} table")),
    }
}

fn int_list(state: &Obj, key: &str) -> Result<Vec<i128>, String> {
    match state.get(key) {
        Some(Json::List(items)) => items
            .iter()
            .map(|v| {
                v.as_int()
                    .ok_or_else(|| format!("merge state {key} holds a non-integer"))
            })
            .collect(),
        _ => Err(format!("merge state has no {key} list")),
    }
}

/// `merge.ingest`: append every unread whole record of every source file to
/// its bucket.  Returns the records added.
pub fn ingest(
    state: &mut Obj,
    root: &Path,
    work: &Path,
    campaign: Option<&Json>,
) -> Result<u64, String> {
    let nb = state
        .get("buckets")
        .and_then(Json::as_int)
        .ok_or("merge state has no bucket count")? as u64;
    let bucket_dir = work.join("buckets");
    fs::create_dir_all(&bucket_dir).map_err(|e| e.to_string())?;
    let mut dirty: BTreeSet<i128> = int_list(state, "dirty")?.into_iter().collect();
    let mut added = 0u64;
    for rel in source_files(root) {
        let path = root.join(&rel);
        if let Some(campaign) = campaign {
            let meta_path = PathBuf::from(format!("{}.json", path.display()));
            if !meta_path.exists() {
                log(&format!("pending manifest: {rel}"));
                continue;
            }
            let meta = read_json(&meta_path)?;
            verify_envelope(&path, &meta, campaign)?;
            let digest = meta
                .as_obj()
                .and_then(|m| m.get("sha256"))
                .cloned()
                .unwrap_or(Json::Null);
            let digests = obj_mut(state, "digests")?;
            if digests.get(&rel).is_some_and(|d| *d != digest) {
                return Err(format!("immutable source changed: {rel}"));
            }
            digests.set(&rel, digest);
        }
        let (head, stride, v2) = corpus_format(&path).map_err(|e| format!("{rel}: {e}"))?;
        let offsets = obj_mut(state, "offsets")?;
        let committed = match offsets.get(&rel) {
            Some(v) => v
                .as_int()
                .ok_or_else(|| format!("offset of {rel} is not an integer"))?,
            None => head as i128,
        };
        let done = committed.max(head as i128) as u64;
        let size = fs::metadata(&path).map_err(|e| e.to_string())?.len();
        let whole = if size > head {
            size - (size - head) % stride
        } else {
            head
        };
        if whole < done {
            return Err(format!("source shrank below committed offset: {rel}"));
        }
        if whole == done {
            continue;
        }
        let mut f = File::open(&path).map_err(|e| e.to_string())?;
        f.seek(SeekFrom::Start(done)).map_err(|e| e.to_string())?;
        let chunk = (READ_BYTES / stride) * stride;
        let mut remaining = whole - done;
        while remaining > 0 {
            let data =
                read_up_to(&mut f, remaining.min(chunk) as usize).map_err(|e| e.to_string())?;
            if data.is_empty() || !(data.len() as u64).is_multiple_of(stride) {
                return Err(format!("source truncated during ingest: {rel}"));
            }
            remaining -= data.len() as u64;
            let mut grouped: std::collections::BTreeMap<u64, Vec<u8>> = Default::default();
            let mut n = 0u64;
            for rec in data.chunks_exact(stride as usize) {
                let r = key_record(rec, v2);
                grouped
                    .entry(bucket_of(r[1], r[2], r[3], nb))
                    .or_default()
                    .extend_from_slice(&record_bytes(&r));
                n += 1;
            }
            for (b, bytes) in grouped {
                let bucket_path = bucket_dir.join(bucket_name(b));
                if bucket_path
                    .metadata()
                    .is_ok_and(|m| m.len() % RECORD_BYTES != 0)
                {
                    return Err(
                        "torn bucket; rebuild derived merge directory from immutable corpora"
                            .into(),
                    );
                }
                let mut out = OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&bucket_path)
                    .map_err(|e| e.to_string())?;
                out.write_all(&bytes).map_err(|e| e.to_string())?;
                out.sync_all().map_err(|e| e.to_string())?;
                dirty.insert(b as i128);
            }
            added += n;
        }
        obj_mut(state, "offsets")?.set(&rel, Json::Int(whole as i128));
    }
    sync_dir(&bucket_dir).map_err(|e| e.to_string())?;
    state.set(
        "dirty",
        Json::List(dirty.into_iter().map(Json::Int).collect()),
    );
    Ok(added)
}

/// One collision, as `merge.detect` builds it.
fn collision(a: &[u64; 4], b: &[u64; 4], bucket: i128) -> Obj {
    Obj::new()
        .with("seedA", Json::Str(format!("{:016x}", a[0])))
        .with("seedB", Json::Str(format!("{:016x}", b[0])))
        .with(
            "key",
            Json::Str(format!("{:016x}{:016x}{:016x}", a[3], a[2], a[1])),
        )
        .with("bucket", Json::Int(bucket))
}

/// `merge.detect`: sort each dirty bucket by `(k0, k1, k2, seed)`, drop exact
/// duplicates, rewrite it, and return the adjacent equal-key pairs and the
/// records kept.
pub fn detect(state: &mut Obj, work: &Path) -> Result<(Vec<Obj>, u64), String> {
    let bucket_dir = work.join("buckets");
    let mut found = Vec::new();
    let mut total = 0u64;
    for b in int_list(state, "dirty")? {
        let path = bucket_dir.join(bucket_name(b as u64));
        let meta = path
            .metadata()
            .map_err(|_| format!("committed bucket missing: {}", path.display()))?;
        if meta.len() % RECORD_BYTES != 0 {
            return Err(format!("torn bucket: {}", path.display()));
        }
        let bytes = fs::read(&path).map_err(|e| e.to_string())?;
        if bytes.is_empty() {
            continue;
        }
        let mut recs: Vec<[u64; 4]> = bytes
            .chunks_exact(RECORD_BYTES as usize)
            .map(|r| key_record(r, false))
            .collect();
        recs.sort_unstable_by_key(|r| (r[1], r[2], r[3], r[0]));
        recs.dedup();
        for pair in recs.windows(2) {
            if pair[0][1..] == pair[1][1..] {
                found.push(collision(&pair[0], &pair[1], b));
            }
        }
        let tmp = PathBuf::from(format!("{}.tmp", path.display()));
        let mut out = File::create(&tmp).map_err(|e| e.to_string())?;
        let mut body = Vec::with_capacity(recs.len() * RECORD_BYTES as usize);
        for r in &recs {
            body.extend_from_slice(&record_bytes(r));
        }
        out.write_all(&body).map_err(|e| e.to_string())?;
        out.sync_all().map_err(|e| e.to_string())?;
        drop(out);
        fs::rename(&tmp, &path).map_err(|e| e.to_string())?;
        sync_dir(&bucket_dir).map_err(|e| e.to_string())?;
        total += recs.len() as u64;
    }
    state.set("dirty", Json::List(Vec::new()));
    Ok((found, total))
}

/// Bucket file names, sorted.
fn bucket_files(work: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(work.join("buckets"))
        .map(|it| {
            it.flatten()
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .filter(|n| n.ends_with(".bin"))
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    names
}

/// `merge.bucketTotal`: records across every bucket file.
pub fn bucket_total(work: &Path) -> u64 {
    let dir = work.join("buckets");
    bucket_files(work)
        .iter()
        .map(|n| fs::metadata(dir.join(n)).map_or(0, |m| m.len()))
        .sum::<u64>()
        / RECORD_BYTES
}

/// `merge.verifyBuckets`.  Buckets are hashed in parallel; the failure
/// reported is still the first in `state.json` order, as in the Python.
fn verify_buckets(state: &Obj, work: &Path) -> Result<(), String> {
    let Some(Json::Obj(hashes)) = state.get("bucketHashes") else {
        return Ok(());
    };
    let entries: Vec<&(String, Json)> = hashes.iter().collect();
    let intact: Vec<bool> = entries
        .par_iter()
        .map(|(name, digest)| {
            let path = work.join("buckets").join(name);
            path.is_file() && sha256_file(&path).is_ok_and(|d| Json::Str(d) == *digest)
        })
        .collect();
    match intact.iter().position(|ok| !ok) {
        Some(i) => Err(format!(
            "bucket integrity failure; rebuild merge state from immutable source corpora: {}",
            entries[i].0
        )),
        None => Ok(()),
    }
}

/// `merge.hashBuckets`, in parallel.
fn hash_buckets(state: &mut Obj, work: &Path) -> Result<(), String> {
    let dir = work.join("buckets");
    if !dir.is_dir() {
        return Err(format!("{}: no such directory", dir.display()));
    }
    let names = bucket_files(work);
    let digests: Vec<io::Result<String>> = names
        .par_iter()
        .map(|name| sha256_file(&dir.join(name)))
        .collect();
    let mut hashes = Obj::new();
    for (name, digest) in names.iter().zip(digests) {
        hashes.set(name, Json::Str(digest.map_err(|e| e.to_string())?));
    }
    state.set("bucketHashes", Json::Obj(hashes));
    Ok(())
}

fn default_state() -> Obj {
    Obj::new()
        .with("version", Json::Int(2))
        .with("buckets", Json::Int(4096))
        .with("offsets", Json::Obj(Obj::new()))
        .with("digests", Json::Obj(Obj::new()))
        .with("bucketHashes", Json::Obj(Obj::new()))
        .with("dirty", Json::List(Vec::new()))
        .with("collisions", Json::List(Vec::new()))
        .with("solved", Json::Null)
}

/// `merge.loadState`.
fn load_state(path: &Path) -> Result<Obj, String> {
    if !path.exists() {
        return Ok(default_state());
    }
    match read_json(path)? {
        Json::Obj(o) => Ok(o),
        _ => Err(format!("{}: not a JSON object", path.display())),
    }
}

// ── Solving ────────────────────────────────────────────────────────

/// What `merge.solve` learns from the client.
struct Solved {
    returncode: i128,
    k: Option<String>,
    verified: bool,
    matches_published: Option<String>,
    tail: Vec<String>,
}

/// An exception `runMerge` catches: its Python class name is recorded.
struct SolveFailure(&'static str);

fn os_error_name(e: &io::Error) -> &'static str {
    match e.kind() {
        io::ErrorKind::NotFound => "FileNotFoundError",
        io::ErrorKind::PermissionDenied => "PermissionError",
        io::ErrorKind::AlreadyExists => "FileExistsError",
        _ => "OSError",
    }
}

/// `re.search(r"matches the (?:published|planted) (?:solution|discrete log): (\w+)", t)`.
fn matches_published(t: &str) -> Option<String> {
    const PREFIXES: [&str; 4] = [
        "matches the published solution: ",
        "matches the published discrete log: ",
        "matches the planted solution: ",
        "matches the planted discrete log: ",
    ];
    for (at, _) in t.char_indices() {
        for p in PREFIXES {
            if t[at..].starts_with(p) {
                let w: String = t[at + p.len()..]
                    .chars()
                    .take_while(|c| c.is_alphanumeric() || *c == '_')
                    .collect();
                if !w.is_empty() {
                    return Some(w);
                }
            }
        }
    }
    None
}

/// Run `cmd` to completion or `timeout`, capturing both streams as text.
fn run_captured(cmd: &[String], timeout: Duration) -> Result<(i128, String), SolveFailure> {
    let mut child = Command::new(&cmd[0])
        .args(&cmd[1..])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| SolveFailure(os_error_name(&e)))?;
    let mut out = child.stdout.take().expect("piped stdout");
    let mut err = child.stderr.take().expect("piped stderr");
    let out_reader = std::thread::spawn(move || {
        let mut buf = Vec::new();
        let _ = out.read_to_end(&mut buf);
        buf
    });
    let err_reader = std::thread::spawn(move || {
        let mut buf = Vec::new();
        let _ = err.read_to_end(&mut buf);
        buf
    });
    let start = Instant::now();
    let status = loop {
        match child.try_wait() {
            Ok(Some(status)) => break status,
            Ok(None) => {}
            Err(e) => return Err(SolveFailure(os_error_name(&e))),
        }
        let elapsed = start.elapsed();
        if elapsed >= timeout {
            let _ = child.kill();
            let _ = child.wait();
            let _ = out_reader.join();
            let _ = err_reader.join();
            return Err(SolveFailure("TimeoutExpired"));
        }
        std::thread::sleep((timeout - elapsed).min(Duration::from_millis(5)));
    };
    let stdout = out_reader.join().unwrap_or_default();
    let stderr = err_reader.join().unwrap_or_default();
    let code = status
        .code()
        .map(i128::from)
        .or_else(|| status.signal().map(|s| -i128::from(s)))
        .unwrap_or(-1);
    let mut text = String::from_utf8_lossy(&stdout).into_owned();
    text.push_str(&String::from_utf8_lossy(&stderr));
    Ok((code, text))
}

fn pair_field<'a>(pair: &'a Obj, key: &str) -> Result<&'a str, String> {
    pair.get(key)
        .and_then(Json::as_str)
        .ok_or_else(|| format!("collision record has no {key}"))
}

/// `merge.solve`: hand the two colliding records to the host client and read
/// its answer.  `Err(Ok(..))` is a failure `runMerge` records and retries;
/// `Err(Err(..))` ends the run, as an uncaught exception does.
fn solve(
    pair: &Obj,
    client: &Path,
    curve: i128,
    work: &Path,
    extra: &[String],
    timeout: Duration,
    walk: &str,
) -> Result<Solved, Result<SolveFailure, String>> {
    let (a, b) = (
        pair_field(pair, "seedA").map_err(Err)?,
        pair_field(pair, "seedB").map_err(Err)?,
    );
    let seed = |s: &str| {
        u64::from_str_radix(s, 16)
            .map_err(|_| Err(format!("invalid literal for int() with base 16: '{s}'")))
    };
    let (seed_a, seed_b) = (seed(a)?, seed(b)?);
    let key = pair_field(pair, "key").map_err(Err)?;
    let key = BigUint::parse_bytes(key.as_bytes(), 16)
        .ok_or_else(|| Err(format!("invalid literal for int() with base 16: '{key}'")))?;
    let digits = key.to_u64_digits();
    let w = |i: usize| digits.get(i).copied().unwrap_or(0);
    let words = [w(0), w(1), w(2)];
    let pair_file = work.join(format!("pair-{a}-{b}.bin"));
    let mut body = Vec::new();
    if walk == "table" {
        body.extend_from_slice(DP_MAGIC_TABLE3);
        body.extend_from_slice(&3u32.to_le_bytes());
        body.extend_from_slice(&(RECORD_BYTES as u32).to_le_bytes());
    }
    for s in [seed_a, seed_b] {
        body.extend_from_slice(&record_bytes(&[s, words[0], words[1], words[2]]));
    }
    fs::write(&pair_file, &body).map_err(|e| Ok(SolveFailure(os_error_name(&e))))?;
    let mut cmd: Vec<String> = [
        client.display().to_string(),
        "--curve".into(),
        curve.to_string(),
        "--threads".into(),
        "1".into(),
        "--steps".into(),
        "1".into(),
        "--launches".into(),
        "1".into(),
        "--verify".into(),
        "0".into(),
        "--run-id".into(),
        "65535".into(),
        "--load".into(),
        pair_file.display().to_string(),
    ]
    .into();
    cmd.extend(extra.iter().cloned());
    log(&format!("solving: {}", cmd.join(" ")));
    let (returncode, out) = run_captured(&cmd, timeout).map_err(Ok)?;
    let lines = py_splitlines(py_strip(&out));
    let tail = lines[lines.len().saturating_sub(12)..].to_vec();
    let mut solved = Solved {
        returncode,
        k: None,
        verified: false,
        matches_published: None,
        tail,
    };
    for line in py_splitlines(&out) {
        let t = py_strip(&line);
        if let Some(k) = t.strip_prefix("k = ") {
            solved.k = Some(py_strip(k).to_string());
        }
        if t.contains("verified [k]P == Q") {
            solved.verified = true;
        }
        if let Some(m) = matches_published(t) {
            solved.matches_published = Some(m);
        }
    }
    solved.verified = solved.verified
        && solved.k.as_deref().is_some_and(|k| !k.is_empty())
        && solved.returncode == 0
        && solved.matches_published.as_deref() != Some("NO");
    Ok(solved)
}

// ── The merge ──────────────────────────────────────────────────────

/// `merge.py`'s command line, after argument parsing.
#[derive(Clone, Debug)]
pub struct MergeArgs {
    pub work: PathBuf,
    pub local: Option<PathBuf>,
    pub s3: Option<String>,
    pub client: PathBuf,
    pub curve: i128,
    pub walk: String,
    pub buckets: i128,
    pub detect_only: bool,
    pub campaign: Option<PathBuf>,
    pub legacy: bool,
    pub solve_timeout: f64,
    pub dp_weight: Option<i128>,
    pub client_arg: Vec<String>,
}

/// How a run ends: argparse's exit 2, an uncaught exception's exit 1.
#[derive(Debug)]
pub enum MergeError {
    Usage(String),
    Fatal(String),
}

impl MergeError {
    pub fn exit_code(&self) -> i32 {
        match self {
            MergeError::Usage(_) => 2,
            MergeError::Fatal(_) => 1,
        }
    }
}

impl From<String> for MergeError {
    fn from(s: String) -> Self {
        MergeError::Fatal(s)
    }
}

/// Progress goes to stderr, so stdout is exactly one JSON summary.
fn log(msg: &str) {
    let tm = broken_down_now(false);
    eprintln!("{:02}:{:02}:{:02} {msg}", tm.tm_hour, tm.tm_min, tm.tm_sec);
}

/// `time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())`.
fn utc_stamp() -> String {
    let tm = broken_down_now(true);
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        tm.tm_year + 1900,
        tm.tm_mon + 1,
        tm.tm_mday,
        tm.tm_hour,
        tm.tm_min,
        tm.tm_sec
    )
}

fn broken_down_now(utc: bool) -> libc::tm {
    let t = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs()) as libc::time_t;
    // SAFETY: gmtime_r and localtime_r only write the `tm` they are given.
    let mut tm: libc::tm = unsafe { std::mem::zeroed() };
    unsafe {
        if utc {
            libc::gmtime_r(&t, &mut tm);
        } else {
            libc::localtime_r(&t, &mut tm);
        }
    }
    tm
}

/// `merge.splitMaxIters`.
fn split_max_iters(args: &Json) -> (Vec<Json>, Option<i128>) {
    let list = match args {
        Json::List(l) => l.clone(),
        _ => Vec::new(),
    };
    match list.iter().position(|a| a.as_str() == Some("--max-iters")) {
        None => (list, None),
        Some(i) => {
            let value = list
                .get(i + 1)
                .and_then(Json::as_str)
                .and_then(|s| s.trim().parse::<i128>().ok());
            let mut rest = list[..i].to_vec();
            rest.extend_from_slice(list.get(i + 2..).unwrap_or(&[]));
            (rest, value)
        }
    }
}

/// `merge.bindingCompatible`: everything must match except that the solver's
/// `--max-iters` may rise.
pub fn binding_compatible(old: &Obj, new: &Obj) -> bool {
    if old == new {
        return true;
    }
    let none = Json::Null;
    let (old_args, before) = split_max_iters(old.get("clientArgs").unwrap_or(&none));
    let (new_args, after) = split_max_iters(new.get("clientArgs").unwrap_or(&none));
    let mut a = old.clone();
    a.set("clientArgs", Json::List(old_args));
    let mut b = new.clone();
    b.set("clientArgs", Json::List(new_args));
    a == b && matches!((before, after), (Some(x), Some(y)) if y > x)
}

/// `json.dumps(v, indent=1)` and a newline, as `print` writes it.
fn print_json(out: &mut dyn Write, v: &Json) -> Result<(), MergeError> {
    let text = dumps(v, Style::Indent(1), false)?;
    writeln!(out, "{text}").map_err(|e| MergeError::Fatal(e.to_string()))
}

fn collisions(state: &Obj) -> Result<&Vec<Json>, String> {
    match state.get("collisions") {
        Some(Json::List(l)) => Ok(l),
        _ => Err("merge state has no collisions list".into()),
    }
}

fn collisions_mut(state: &mut Obj) -> Result<&mut Vec<Json>, String> {
    match state.get_mut("collisions") {
        Some(Json::List(l)) => Ok(l),
        _ => Err("merge state has no collisions list".into()),
    }
}

fn verified(pair: &Json) -> bool {
    pair.as_obj()
        .and_then(|p| p.get("result"))
        .and_then(Json::as_obj)
        .and_then(|r| r.get("verified"))
        .is_some_and(Json::truthy)
}

fn opt_str(v: &Option<String>) -> Json {
    v.as_ref().map_or(Json::Null, |s| Json::Str(s.clone()))
}

fn py_repr(v: &Option<String>) -> String {
    v.clone().unwrap_or_else(|| "None".into())
}

/// `merge.main` after parsing: validate, bind, lock and run one pass.
pub fn run(mut args: MergeArgs, stdout: &mut dyn Write) -> Result<i32, MergeError> {
    if args.campaign.is_some() == args.legacy {
        return Err(MergeError::Usage(
            "choose exactly one of --campaign or --legacy".into(),
        ));
    }
    if args.buckets < 1 || args.buckets & (args.buckets - 1) != 0 {
        return Err(MergeError::Usage(
            "--buckets must be a positive power of two".into(),
        ));
    }
    if args.solve_timeout <= 0.0 {
        return Err(MergeError::Usage("--solve-timeout must be positive".into()));
    }
    let mut campaign = None;
    if let Some(path) = args.campaign.clone() {
        if !args.client_arg.is_empty() {
            return Err(MergeError::Usage(
                "strict merge forbids --client-arg overrides".into(),
            ));
        }
        let config = match read_json(&path)? {
            Json::Obj(o) => o,
            _ => {
                return Err(MergeError::Fatal(format!(
                    "{}: not a JSON object",
                    path.display()
                )))
            }
        };
        campaign = Some(campaign_contract(&config)?);
        let host = sha256_file(&args.client)
            .map_err(|e| MergeError::Fatal(format!("{}: {e}", args.client.display())))?;
        if config.get("hostBinarySha256") != Some(&Json::Str(host)) {
            return Err(MergeError::Usage(
                "solver binary hash differs from pinned campaign hostBinarySha256".into(),
            ));
        }
        let int = |key: &str| {
            config
                .get(key)
                .and_then(Json::as_int)
                .ok_or_else(|| MergeError::Fatal(format!("campaign config has no integer {key}")))
        };
        args.curve = int("curve")?;
        args.dp_weight = Some(int("dpWeight")?);
        args.walk = match config.get("walk") {
            Some(Json::Str(w)) => w.clone(),
            _ => "sigma".into(),
        };
        if let Some(cap) = config.get("maxIters").filter(|v| v.truthy()) {
            args.client_arg = vec!["--max-iters".into(), dumps(cap, Style::Compact, false)?];
        }
    }
    if let Some(dp) = args.dp_weight {
        let mut with = vec!["--dp-weight".to_string(), dp.to_string()];
        with.append(&mut args.client_arg);
        args.client_arg = with;
    }
    fs::create_dir_all(&args.work).map_err(|e| MergeError::Fatal(e.to_string()))?;
    let lock = OpenOptions::new()
        .create(true)
        .append(true)
        .read(true)
        .open(args.work.join("merge.lock"))
        .map_err(|e| MergeError::Fatal(e.to_string()))?;
    use std::os::unix::io::AsRawFd;
    // SAFETY: flock on a descriptor this function owns for the whole run.
    if unsafe { libc::flock(lock.as_raw_fd(), libc::LOCK_EX | libc::LOCK_NB) } != 0 {
        return Err(MergeError::Fatal(format!(
            "{}: {}",
            args.work.join("merge.lock").display(),
            io::Error::last_os_error()
        )));
    }
    let code = run_merge(&args, campaign.as_ref(), stdout);
    drop(lock);
    code
}

fn run_merge(
    args: &MergeArgs,
    campaign: Option<&Json>,
    stdout: &mut dyn Write,
) -> Result<i32, MergeError> {
    if let Some(campaign) = campaign {
        bind_directory(&args.work, campaign)?;
    }
    let state_path = args.work.join("state.json");
    let mut state = load_state(&state_path)?;
    if state.get("version") != Some(&Json::Int(2)) {
        return Err(MergeError::Fatal(
            "old bucket layout; use a fresh merge directory (preserve source corpora)".into(),
        ));
    }
    verify_buckets(&state, &args.work)?;
    let mut binding = Obj::new()
        .with(
            "campaignId",
            campaign.map_or(Json::Null, |c| Json::str(campaign_id(c))),
        )
        .with("curve", Json::Int(args.curve))
        .with("dpWeight", args.dp_weight.map_or(Json::Null, Json::Int))
        .with(
            "clientArgs",
            Json::List(args.client_arg.iter().map(|a| Json::str(a)).collect()),
        );
    if args.walk == "table" {
        binding.set("walk", Json::str("table-v3"));
    }
    if let Some(Json::Obj(old)) = state.get("binding") {
        if !binding_compatible(old, &binding) {
            return Err(MergeError::Fatal(
                "merge state belongs to different campaign/solve parameters".into(),
            ));
        }
    } else if state.contains("binding") {
        return Err(MergeError::Fatal(
            "merge state belongs to different campaign/solve parameters".into(),
        ));
    }
    state.set("binding", Json::Obj(binding));
    if !state_path.exists() {
        state.set("buckets", Json::Int(args.buckets));
    }
    if let Some(solved) = state.get("solved").filter(|s| s.truthy()) {
        let k = solved
            .as_obj()
            .and_then(|s| s.get("k"))
            .map(|k| match k {
                Json::Str(s) => s.clone(),
                other => dumps(other, Style::Compact, false).unwrap_or_default(),
            })
            .unwrap_or_default();
        log(&format!("already solved: k = {k}"));
        print_json(stdout, solved)?;
        return Ok(0);
    }

    let root = match &args.s3 {
        Some(uri) => {
            let cache = args.work.join("s3cache");
            fs::create_dir_all(&cache).map_err(|e| MergeError::Fatal(e.to_string()))?;
            let cmd = [
                "aws",
                "s3",
                "sync",
                uri,
                &cache.display().to_string(),
                "--only-show-errors",
            ];
            log(&cmd.join(" "));
            let status = Command::new(cmd[0])
                .args(&cmd[1..])
                .status()
                .map_err(|e| MergeError::Fatal(format!("aws: {e}")))?;
            if !status.success() {
                return Err(MergeError::Fatal(format!(
                    "Command {cmd:?} returned non-zero exit status {}",
                    status.code().unwrap_or(-1)
                )));
            }
            cache
        }
        None => args.local.clone().expect("--local or --s3"),
    };
    let t0 = Instant::now();
    let added = ingest(&mut state, &root, &args.work, campaign)?;
    hash_buckets(&mut state, &args.work)?;
    atomic_json(&state_path, &Json::Obj(state.clone()))?;
    log(&format!(
        "ingested {added} new records in {:.1} s ({} buckets to re-sort)",
        t0.elapsed().as_secs_f64(),
        int_list(&state, "dirty")?.len()
    ));
    let t1 = Instant::now();
    let (found, sorted) = detect(&mut state, &args.work)?;
    let known: BTreeSet<(String, String)> = collisions(&state)?
        .iter()
        .filter_map(|c| {
            let c = c.as_obj()?;
            Some((
                c.get("seedA")?.as_str()?.to_string(),
                c.get("seedB")?.as_str()?.to_string(),
            ))
        })
        .collect();
    let new: Vec<Obj> = found
        .into_iter()
        .filter(|c| {
            let k = (
                c.get("seedA")
                    .and_then(Json::as_str)
                    .unwrap_or("")
                    .to_string(),
                c.get("seedB")
                    .and_then(Json::as_str)
                    .unwrap_or("")
                    .to_string(),
            );
            !known.contains(&k)
        })
        .collect();
    let first_new = collisions(&state)?.len();
    collisions_mut(&mut state)?.extend(new.iter().cloned().map(Json::Obj));
    hash_buckets(&mut state, &args.work)?;
    atomic_json(&state_path, &Json::Obj(state.clone()))?;
    let pending: Vec<usize> = collisions(&state)?
        .iter()
        .enumerate()
        .filter(|(_, c)| !verified(c))
        .map(|(i, _)| i)
        .collect();
    let corpus = bucket_total(&args.work);
    log(&format!(
        "sorted {sorted} records in {:.1} s; corpus {corpus} records; {} collision(s), {} new, {} unsolved",
        t1.elapsed().as_secs_f64(),
        collisions(&state)?.len(),
        new.len(),
        pending.len()
    ));
    let n_collisions = collisions(&state)?.len();
    let (n_new, n_pending) = (new.len(), pending.len());
    // `summary["new"]` aliases the pairs it lists, so a solve that marks
    // one shows in the summary; it is rendered from the state each time.
    let summary = |state: &Obj, solution: Json| -> Result<Json, String> {
        let all = collisions(state)?;
        Ok(Json::Obj(
            Obj::new()
                .with("added", Json::Int(added as i128))
                .with("corpus", Json::Int(corpus as i128))
                .with("collisions", Json::Int(n_collisions as i128))
                .with(
                    "new",
                    Json::List(all[first_new..first_new + n_new].to_vec()),
                )
                .with("pending", Json::Int(n_pending as i128))
                .with("solution", solution),
        ))
    };
    if args.detect_only || pending.is_empty() {
        print_json(stdout, &summary(&state, Json::Null)?)?;
        return Ok(0);
    }
    let timeout = Duration::from_secs_f64(args.solve_timeout.min(1e9));
    for i in pending {
        let pair = match &collisions(&state)?[i] {
            Json::Obj(o) => o.clone(),
            _ => {
                return Err(MergeError::Fatal(
                    "collision record is not an object".into(),
                ))
            }
        };
        let res = solve(
            &pair,
            &args.client,
            args.curve,
            &args.work,
            &args.client_arg,
            timeout,
            &args.walk,
        );
        let solved = match res {
            Ok(s) => s,
            Err(Ok(SolveFailure(name))) => {
                if let Json::Obj(p) = &mut collisions_mut(&mut state)?[i] {
                    p.set(
                        "result",
                        Json::Obj(
                            Obj::new()
                                .with("verified", Json::Bool(false))
                                .with("error", Json::str(name)),
                        ),
                    );
                }
                atomic_json(&state_path, &Json::Obj(state.clone()))?;
                continue;
            }
            Err(Err(fatal)) => return Err(MergeError::Fatal(fatal)),
        };
        let tail = Json::List(solved.tail.iter().map(|l| Json::str(l)).collect());
        let pair_now = {
            let Json::Obj(p) = &mut collisions_mut(&mut state)?[i] else {
                unreachable!("checked above")
            };
            p.set(
                "result",
                Json::Obj(
                    Obj::new()
                        .with("k", opt_str(&solved.k))
                        .with("verified", Json::Bool(solved.verified))
                        .with("returncode", Json::Int(solved.returncode))
                        .with("tail", tail.clone()),
                ),
            );
            p.clone()
        };
        atomic_json(&state_path, &Json::Obj(state.clone()))?;
        log(&format!(
            "seeds {} / {}: k = {}, verified {}",
            pair_field(&pair_now, "seedA").unwrap_or(""),
            pair_field(&pair_now, "seedB").unwrap_or(""),
            py_repr(&solved.k),
            if solved.verified { "True" } else { "False" }
        ));
        if solved.k.as_deref().is_some_and(|k| !k.is_empty()) && solved.verified {
            let res = Json::Obj(
                Obj::new()
                    .with("pair", Json::Obj(pair_now))
                    .with("returncode", Json::Int(solved.returncode))
                    .with("k", opt_str(&solved.k))
                    .with("verified", Json::Bool(solved.verified))
                    .with("matchesPublished", opt_str(&solved.matches_published))
                    .with("tail", tail)
                    .with("when", Json::Str(utc_stamp())),
            );
            state.set("solved", res.clone());
            atomic_json(&state_path, &Json::Obj(state.clone()))?;
            atomic_json(&args.work.join("solution.json"), &res)?;
            if let Some(uri) = &args.s3 {
                let dest = format!(
                    "{}/solution.json",
                    uri.trim_end_matches('/')
                        .rsplit_once('/')
                        .map_or(uri.trim_end_matches('/'), |(head, _)| head)
                );
                let _ = Command::new("aws")
                    .args([
                        "s3",
                        "cp",
                        &args.work.join("solution.json").display().to_string(),
                        &dest,
                        "--only-show-errors",
                    ])
                    .status();
            }
            print_json(stdout, &summary(&state, res)?)?;
            return Ok(0);
        }
    }
    print_json(stdout, &summary(&state, Json::Null)?)?;
    Ok(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `merge.py`'s numpy hash on 2026-10-06.
    #[test]
    fn buckets_match_numpy() {
        let keys = [
            (0u64, 0u64, 0u64),
            (1, 2, 3),
            (u64::MAX, 0x1234_5678_90AB_CDEF, 7),
            (0x83E0_790D_2530_939F, 0xCDD9_D226_E697_7EDB, 5),
        ];
        for (nb, want) in [
            (16u64, [9u64, 10, 10, 11]),
            (4096, [1465, 2826, 858, 1883]),
            (1 << 20, [628153, 854794, 156506, 329563]),
        ] {
            let got: Vec<u64> = keys
                .iter()
                .map(|&(a, b, c)| bucket_of(a, b, c, nb))
                .collect();
            assert_eq!(got, want, "{nb} buckets");
        }
    }

    fn config(extra: &[(&str, Json)]) -> Obj {
        let mut c = parse_json(&format!(
            r#"{{"storageProtocol": "{PROTOCOL}", "curve": 41, "dpWeight": 13, "maxIters": 100000,
                "packed": false, "workers": 1, "batch": 32, "extraArgs": [],
                "binarySha256": "{}", "hostBinarySha256": "{}", "sourceSha256": "{}"}}"#,
            "a".repeat(64),
            "b".repeat(64),
            "c".repeat(64)
        ))
        .unwrap();
        let Json::Obj(ref mut o) = c else { panic!() };
        for (k, v) in extra {
            o.set(k, v.clone());
        }
        o.clone()
    }

    /// `protocol.campaignContract(...)["id"]` on 2026-10-06.
    #[test]
    fn campaign_ids_match_protocol_py() {
        let id = |c: &Obj| campaign_id(&campaign_contract(c).unwrap()).to_string();
        assert_eq!(
            id(&config(&[])),
            "321ccd25adc2af730f41fe0b94b4455e41d1194ad45f5516c66f5dd74e63804e"
        );
        assert_eq!(
            id(&config(&[("walk", Json::str("table"))])),
            "51ebe68681944bb90b761906ce82b4a78000ce200d70fe85931dbfe994c2322f"
        );
        assert_eq!(
            id(&config(&[
                ("maxIters", Json::Int(1 << 32)),
                ("contractMaxIters", Json::Int(100000))
            ])),
            "321ccd25adc2af730f41fe0b94b4455e41d1194ad45f5516c66f5dd74e63804e"
        );
        assert_eq!(
            id(&config(&[
                ("curve", Json::Int(131)),
                ("dpWeight", Json::Int(32)),
                ("packed", Json::Bool(true)),
                ("workers", Json::Int(385024)),
                ("batch", Json::Int(16)),
            ])),
            "dee1ad8798583a72cff02a03f169ee3fe7ee4e2000f807a86b3930ddaa87d945"
        );
    }

    #[test]
    fn the_contract_refuses_what_protocol_py_refuses() {
        for (field, value) in [
            ("extraArgs", parse_json(r#"["--dp-weight", "3"]"#).unwrap()),
            ("binarySha256", Json::str("")),
            ("curve", Json::Bool(true)),
            ("curve", Json::Int(97)),
            ("packed", Json::Bool(true)),
            ("dpWeight", Json::Int(42)),
            ("walk", Json::str("random")),
            ("maxIters", Json::Int(99)),
            ("storageProtocol", Json::str("other")),
        ] {
            let mut c = config(&[]);
            if field == "maxIters" {
                c.set("contractMaxIters", Json::Int(100000));
            }
            c.set(field, value);
            assert!(campaign_contract(&c).is_err(), "{field}");
        }
    }

    #[test]
    fn file_hash_agrees_with_the_reference_sha256() {
        let dir = scratch("filehash");
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        for len in [0usize, 1, 55, 56, 64, 1000, 65_537, (1 << 20) + 7] {
            let data: Vec<u8> = (0..len)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    x as u8
                })
                .collect();
            let path = dir.join(format!("{len}.bin"));
            fs::write(&path, &data).unwrap();
            assert_eq!(
                sha256_file(&path).unwrap(),
                hex::encode(sha256(&data)),
                "{len} bytes"
            );
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn the_solver_verdict_is_read_like_the_regex() {
        assert_eq!(
            matches_published("  matches the planted discrete log: yes"),
            Some("yes".into())
        );
        assert_eq!(
            matches_published(
                "matches the published solution: !! matches the planted solution: NO"
            ),
            Some("NO".into())
        );
        assert_eq!(matches_published("matches the published solution:"), None);
    }

    #[test]
    fn max_iters_may_rise_and_nothing_else_may_move() {
        let b = |args: &[&str]| -> Obj {
            Obj::new()
                .with("campaignId", Json::str("c"))
                .with("curve", Json::Int(41))
                .with("dpWeight", Json::Int(13))
                .with(
                    "clientArgs",
                    Json::List(args.iter().map(|a| Json::str(a)).collect()),
                )
        };
        let old = b(&["--dp-weight", "13", "--max-iters", "1073741824"]);
        assert!(binding_compatible(&old, &old));
        assert!(binding_compatible(
            &old,
            &b(&["--dp-weight", "13", "--max-iters", "4294967296"])
        ));
        assert!(!binding_compatible(
            &old,
            &b(&["--dp-weight", "13", "--max-iters", "1024"])
        ));
        assert!(!binding_compatible(
            &old,
            &b(&["--dp-weight", "14", "--max-iters", "4294967296"])
        ));
        assert!(!binding_compatible(&old, &b(&["--dp-weight", "13"])));
        assert!(!binding_compatible(&b(&["--dp-weight", "13"]), &old));
    }

    fn rec(seed: u64, key: [u64; 3]) -> [u8; 32] {
        record_bytes(&[seed, key[0], key[1], key[2]])
    }

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "ecc2k-merge-{name}-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(dir.join("corpus")).unwrap();
        dir
    }

    #[test]
    fn all_three_layouts_reach_the_same_buckets() {
        let dir = scratch("layouts");
        let key = [3u64, 4, 5];
        fs::write(dir.join("corpus/v1.bin"), rec(7, key)).unwrap();
        let mut v2 = DP_MAGIC_V2.to_vec();
        v2.extend_from_slice(&[2, 0, 0, 0, 72, 0, 0, 0]);
        v2.extend_from_slice(&8u64.to_le_bytes());
        v2.extend_from_slice(&999u64.to_le_bytes());
        for w in key {
            v2.extend_from_slice(&w.to_le_bytes());
        }
        v2.extend_from_slice(&[1u8; 32]);
        fs::create_dir_all(dir.join("corpus/sub")).unwrap();
        fs::write(dir.join("corpus/sub/v2.bin"), v2).unwrap();
        let mut t3 = DP_MAGIC_TABLE3.to_vec();
        t3.extend_from_slice(&[3, 0, 0, 0, 32, 0, 0, 0]);
        t3.extend_from_slice(&rec(9, key));
        fs::write(dir.join("corpus/t3.bin"), t3).unwrap();
        let mut state = default_state().with("buckets", Json::Int(16));
        assert_eq!(
            ingest(&mut state, &dir.join("corpus"), &dir, None).unwrap(),
            3
        );
        let (found, total) = detect(&mut state, &dir).unwrap();
        assert_eq!(total, 3);
        let seeds: Vec<(&str, &str)> = found
            .iter()
            .map(|c| {
                (
                    c.get("seedA").unwrap().as_str().unwrap(),
                    c.get("seedB").unwrap().as_str().unwrap(),
                )
            })
            .collect();
        assert_eq!(
            seeds,
            [
                ("0000000000000007", "0000000000000008"),
                ("0000000000000008", "0000000000000009")
            ]
        );
        assert_eq!(
            found[0].get("key").unwrap().as_str().unwrap(),
            "000000000000000500000000000000040000000000000003"
        );
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn a_crash_between_append_and_state_only_duplicates() {
        let dir = scratch("crash");
        fs::write(
            dir.join("corpus/a.bin"),
            [rec(1, [5, 6, 7]), rec(2, [5, 6, 7])].concat(),
        )
        .unwrap();
        let fresh = default_state().with("buckets", Json::Int(16));
        let mut state = fresh.clone();
        ingest(&mut state, &dir.join("corpus"), &dir, None).unwrap();
        let mut before = fresh;
        ingest(&mut before, &dir.join("corpus"), &dir, None).unwrap();
        let (found, total) = detect(&mut before, &dir).unwrap();
        assert_eq!((found.len(), total), (1, 2));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn torn_and_shrunk_inputs_fail_closed() {
        let dir = scratch("torn");
        let src = dir.join("corpus/a.bin");
        fs::write(&src, rec(1, [5, 6, 7])).unwrap();
        let mut state = default_state().with("buckets", Json::Int(16));
        ingest(&mut state, &dir.join("corpus"), &dir, None).unwrap();
        let bucket = dir
            .join("buckets")
            .join(bucket_name(bucket_of(5, 6, 7, 16)));
        OpenOptions::new()
            .append(true)
            .open(&bucket)
            .unwrap()
            .write_all(b"torn")
            .unwrap();
        assert!(detect(&mut state.clone(), &dir)
            .unwrap_err()
            .contains("torn bucket"));
        fs::write(&src, b"").unwrap();
        assert!(ingest(&mut state, &dir.join("corpus"), &dir, None)
            .unwrap_err()
            .contains("shrank"));
        fs::remove_dir_all(dir).unwrap();
    }
}
