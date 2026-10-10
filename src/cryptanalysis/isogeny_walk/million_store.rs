//! Deterministic local cache and immutable S3 descriptors for the million-grid certificate.

use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use flate2::{read::GzDecoder, Compression, GzBuilder};
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use super::million::{
    Coordinate, GenerationReceipt, GridSummary, HeaderRecord, NodeRecord, RootRecord,
    VerificationReceipt, RECEIPT_SCHEMA, ROW_DEGREE, SPINE_DEGREE,
};
use super::store;
use crate::cryptanalysis::p256_isogeny_campaign::P256_ICV1;
use crate::hash::sha256::{sha256, Sha256};

pub const CACHE_SCHEMA: &str = "p256.isogeny-grid-cache/v1";
pub const COMPLETE_SCHEMA: &str = "p256.isogeny-grid-s3-complete/v1";
pub const CAIRN_RECEIPT_SCHEMA: &str = "p256.isogeny-grid-cairn-receipt/v1";
pub const CACHE_MANIFEST_FILE: &str = "manifest.json";
pub const DEFAULT_CHUNK_ROWS: u32 = 8;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FileDigest {
    pub sha256: String,
    pub bytes: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EncodedDigest {
    pub encoding: String,
    pub stored_sha256: String,
    pub stored_bytes: u64,
    pub content_sha256: String,
    pub content_bytes: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RowRange {
    pub start: u32,
    pub end_exclusive: u32,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CachePart {
    pub name: String,
    pub file: String,
    pub role: String,
    pub rows: Option<RowRange>,
    pub records: u64,
    pub digest: EncodedDigest,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CacheManifest {
    pub schema: String,
    pub run_id: String,
    pub result_class: String,
    pub source_commit: String,
    pub root_icv1_slug: String,
    pub side: u32,
    pub unique_curves: u64,
    pub chunk_rows: u32,
    pub canonical_file: String,
    pub canonical: EncodedDigest,
    pub generation_receipt: FileDigest,
    pub verification_receipt: FileDigest,
    pub summary: GridSummary,
    pub parts: Vec<CachePart>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StoredObject {
    pub key: String,
    pub encoding: String,
    pub stored_sha256: String,
    pub stored_bytes: u64,
    pub content_sha256: String,
    pub content_bytes: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StoredPart {
    pub name: String,
    pub role: String,
    pub rows: Option<RowRange>,
    pub records: u64,
    pub object: StoredObject,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CompleteMarker {
    pub schema: String,
    pub run_id: String,
    pub run: String,
    pub result_class: String,
    pub source_commit: String,
    pub root_icv1_slug: String,
    pub side: u32,
    pub unique_curves: u64,
    pub chunk_rows: u32,
    pub canonical: StoredObject,
    pub generation_receipt: StoredObject,
    pub verification_receipt: StoredObject,
    pub summary: GridSummary,
    pub parts: Vec<StoredPart>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CairnReceipt {
    pub schema: String,
    pub claim_class: String,
    pub verification: String,
    pub run_id: String,
    pub result_class: String,
    pub source_commit: String,
    pub root_icv1_slug: String,
    pub unique_curves: u64,
    pub complete_key: String,
    pub complete_sha256: String,
    pub canonical_key: String,
    pub canonical_sha256: String,
    pub canonical_content_sha256: String,
    pub record_chain_sha256: String,
}

impl CairnReceipt {
    pub fn artifact(&self) -> Result<serde_json::Value, String> {
        serde_json::to_value(self).map_err(|error| error.to_string())
    }

    pub fn claim_key(&self) -> String {
        self.run_id.clone()
    }
}

pub struct Published {
    pub marker: CompleteMarker,
    pub marker_uri: String,
    pub marker_sha256: String,
    pub created: bool,
    pub cairn_receipt: CairnReceipt,
}

fn hex_digest(digest: [u8; 32]) -> String {
    hex::encode(digest)
}

pub fn digest_reader<R: Read>(mut reader: R) -> Result<FileDigest, String> {
    let mut state = Sha256::new();
    let mut bytes = 0u64;
    let mut buffer = [0u8; 1024 * 1024];
    loop {
        let read = reader
            .read(&mut buffer)
            .map_err(|error| format!("read while hashing: {error}"))?;
        if read == 0 {
            break;
        }
        state.update(&buffer[..read]);
        bytes = bytes
            .checked_add(read as u64)
            .ok_or("byte count overflow while hashing")?;
    }
    Ok(FileDigest {
        sha256: hex_digest(state.finalize()),
        bytes,
    })
}

pub fn digest_file(path: &Path) -> Result<FileDigest, String> {
    let file = File::open(path).map_err(|error| format!("open {}: {error}", path.display()))?;
    digest_reader(BufReader::new(file))
}

pub fn digest_gzip(path: &Path) -> Result<EncodedDigest, String> {
    let stored = digest_file(path)?;
    let file = File::open(path).map_err(|error| format!("open {}: {error}", path.display()))?;
    let content = digest_reader(GzDecoder::new(BufReader::new(file)))?;
    Ok(EncodedDigest {
        encoding: "gzip".into(),
        stored_sha256: stored.sha256,
        stored_bytes: stored.bytes,
        content_sha256: content.sha256,
        content_bytes: content.bytes,
    })
}

fn read_json<T: DeserializeOwned>(path: &Path) -> Result<T, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    serde_json::from_slice(&bytes).map_err(|error| format!("parse {}: {error}", path.display()))
}

fn validate_receipts(
    generation: &GenerationReceipt,
    verification: &VerificationReceipt,
) -> Result<(), String> {
    if generation.schema != RECEIPT_SCHEMA
        || verification.schema != RECEIPT_SCHEMA
        || generation.status != "pass"
        || verification.status != "pass"
    {
        return Err(
            "generation and verification receipts must both be passing grid receipts".into(),
        );
    }
    if generation.source_commit != verification.source_commit
        || generation.uncompressed_bytes != verification.uncompressed_bytes
        || generation.summary != verification.summary
    {
        return Err("generation and verification receipts disagree".into());
    }
    if verification.audit_points != 2 || verification.audit_seed_x != 7 {
        return Err("the independent receipt must use two audit points from x = 7".into());
    }
    Ok(())
}

fn partial_path(path: &Path) -> PathBuf {
    let mut name: OsString = path.as_os_str().to_owned();
    name.push(".partial");
    PathBuf::from(name)
}

fn read_line<R: BufRead>(reader: &mut R, line_number: &mut u64) -> Result<Vec<u8>, String> {
    let mut line = Vec::new();
    let read = reader
        .read_until(b'\n', &mut line)
        .map_err(|error| format!("read certificate: {error}"))?;
    *line_number += 1;
    if read == 0 {
        return Err(format!(
            "line {}: unexpected end of certificate",
            line_number
        ));
    }
    if line.len() > 16 * 1024 * 1024 {
        return Err(format!("line {}: record exceeds 16 MiB", line_number));
    }
    if line.last() != Some(&b'\n') || line.ends_with(b"\r\n") {
        return Err(format!(
            "line {}: records require one LF terminator",
            line_number
        ));
    }
    Ok(line)
}

fn parse_line<T: DeserializeOwned>(line: &[u8], line_number: u64) -> Result<T, String> {
    serde_json::from_slice(line).map_err(|error| format!("line {line_number}: {error}"))
}

fn bind_line(chain: &mut [u8; 32], line: &[u8]) {
    let mut input = Vec::with_capacity(32 + line.len());
    input.extend_from_slice(chain);
    input.extend_from_slice(line);
    *chain = sha256(&input);
}

fn observe_line(raw: &mut Sha256, raw_bytes: &mut u64, line: &[u8]) -> Result<(), String> {
    raw.update(line);
    *raw_bytes = raw_bytes
        .checked_add(line.len() as u64)
        .ok_or("uncompressed byte count overflow")?;
    Ok(())
}

fn write_part(
    dir: &Path,
    name: String,
    role: &str,
    rows: Option<RowRange>,
    records: u64,
    content: &[u8],
) -> Result<CachePart, String> {
    store::segment(&name)?;
    let file_name = format!("{name}.jsonl.gz");
    let path = dir.join(&file_name);
    let file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&path)
        .map_err(|error| format!("create {}: {error}", path.display()))?;
    let buffered = BufWriter::new(file);
    let mut encoder = GzBuilder::new()
        .mtime(0)
        .operating_system(255)
        .write(buffered, Compression::new(6));
    encoder
        .write_all(content)
        .map_err(|error| format!("write {}: {error}", path.display()))?;
    let mut buffered = encoder
        .finish()
        .map_err(|error| format!("finish {}: {error}", path.display()))?;
    buffered
        .flush()
        .map_err(|error| format!("flush {}: {error}", path.display()))?;
    buffered
        .get_ref()
        .sync_all()
        .map_err(|error| format!("sync {}: {error}", path.display()))?;
    let stored = digest_file(&path)?;
    Ok(CachePart {
        name,
        file: file_name,
        role: role.into(),
        rows,
        records,
        digest: EncodedDigest {
            encoding: "gzip".into(),
            stored_sha256: stored.sha256,
            stored_bytes: stored.bytes,
            content_sha256: hex_digest(sha256(content)),
            content_bytes: content.len() as u64,
        },
    })
}

fn append_bound(
    content: &mut Vec<u8>,
    raw: &mut Sha256,
    raw_bytes: &mut u64,
    chain: &mut [u8; 32],
    line: &[u8],
) -> Result<(), String> {
    bind_line(chain, line);
    observe_line(raw, raw_bytes, line)?;
    content.extend_from_slice(line);
    Ok(())
}

/// Split a verified canonical certificate into independently fetchable row
/// slices. Concatenating the decompressed parts in manifest order reproduces
/// the exact canonical uncompressed bytes.
pub fn stage_cache(
    input: &Path,
    generation_path: &Path,
    verification_path: &Path,
    run_id: &str,
    output: &Path,
    chunk_rows: u32,
) -> Result<CacheManifest, String> {
    store::segment(run_id)?;
    if chunk_rows == 0 {
        return Err("chunk_rows must be positive".into());
    }
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()));
    }
    let partial = partial_path(output);
    if partial.exists() {
        return Err(format!(
            "incomplete cache already exists at {}",
            partial.display()
        ));
    }
    fs::create_dir_all(&partial)
        .map_err(|error| format!("create {}: {error}", partial.display()))?;

    let generation: GenerationReceipt = read_json(generation_path)?;
    let verification: VerificationReceipt = read_json(verification_path)?;
    validate_receipts(&generation, &verification)?;
    let generation_digest = digest_file(generation_path)?;
    let verification_digest = digest_file(verification_path)?;
    let canonical_stored = digest_file(input)?;

    let file = File::open(input).map_err(|error| format!("open {}: {error}", input.display()))?;
    let mut reader = BufReader::new(GzDecoder::new(file));
    let mut line_number = 0u64;
    let mut raw = Sha256::new();
    let mut raw_bytes = 0u64;
    let mut chain = [0u8; 32];
    let mut parts = Vec::new();

    let header_line = read_line(&mut reader, &mut line_number)?;
    let header: HeaderRecord = parse_line(&header_line, line_number)?;
    if header.record != "header"
        || header.source_commit != generation.source_commit
        || header.spine_degree != SPINE_DEGREE
        || header.row_degree != ROW_DEGREE
        || header.target_curves != u64::from(header.side) * u64::from(header.side)
        || header.target_curves != generation.summary.unique_curves
    {
        return Err("certificate header disagrees with the verified receipts".into());
    }
    let mut preamble = Vec::new();
    append_bound(
        &mut preamble,
        &mut raw,
        &mut raw_bytes,
        &mut chain,
        &header_line,
    )?;
    let root_line = read_line(&mut reader, &mut line_number)?;
    let root: RootRecord = parse_line(&root_line, line_number)?;
    if root.record != "root" || root.index != 0 || root.coordinate != (Coordinate { x: 0, y: 0 }) {
        return Err("certificate root is not coordinate (0,0) at index zero".into());
    }
    append_bound(
        &mut preamble,
        &mut raw,
        &mut raw_bytes,
        &mut chain,
        &root_line,
    )?;
    for y in 1..header.side {
        let line = read_line(&mut reader, &mut line_number)?;
        let node: NodeRecord = parse_line(&line, line_number)?;
        if node.record != "node"
            || node.index != u64::from(y)
            || node.coordinate != (Coordinate { x: 0, y })
            || node.degree != SPINE_DEGREE
        {
            return Err(format!("line {line_number}: malformed spine position {y}"));
        }
        append_bound(&mut preamble, &mut raw, &mut raw_bytes, &mut chain, &line)?;
    }
    parts.push(write_part(
        &partial,
        "part-0000-preamble".into(),
        "preamble",
        None,
        u64::from(header.side) + 1,
        &preamble,
    )?);

    let mut part_index = 1u32;
    for row_start in (0..header.side).step_by(chunk_rows as usize) {
        let row_end = row_start.saturating_add(chunk_rows).min(header.side);
        let mut content = Vec::new();
        let mut records = 0u64;
        for y in row_start..row_end {
            for x in 1..header.side {
                let line = read_line(&mut reader, &mut line_number)?;
                let node: NodeRecord = parse_line(&line, line_number)?;
                let expected_index = u64::from(header.side)
                    + u64::from(y) * u64::from(header.side - 1)
                    + u64::from(x - 1);
                if node.record != "node"
                    || node.index != expected_index
                    || node.coordinate != (Coordinate { x, y })
                    || node.degree != ROW_DEGREE
                {
                    return Err(format!(
                        "line {line_number}: malformed grid position ({x},{y})"
                    ));
                }
                append_bound(&mut content, &mut raw, &mut raw_bytes, &mut chain, &line)?;
                records += 1;
            }
        }
        parts.push(write_part(
            &partial,
            format!(
                "part-{part_index:04}-rows-{row_start:04}-{last:04}",
                last = row_end - 1
            ),
            "rows",
            Some(RowRange {
                start: row_start,
                end_exclusive: row_end,
            }),
            records,
            &content,
        )?);
        part_index += 1;
    }

    let summary_line = read_line(&mut reader, &mut line_number)?;
    let summary: GridSummary = parse_line(&summary_line, line_number)?;
    if summary != generation.summary || summary.record_chain_sha256 != hex_digest(chain) {
        return Err("certificate summary or record chain differs from the receipts".into());
    }
    observe_line(&mut raw, &mut raw_bytes, &summary_line)?;
    parts.push(write_part(
        &partial,
        format!("part-{part_index:04}-summary"),
        "summary",
        None,
        1,
        &summary_line,
    )?);
    let mut trailing = Vec::new();
    reader
        .read_to_end(&mut trailing)
        .map_err(|error| format!("read certificate tail: {error}"))?;
    if !trailing.is_empty() {
        return Err("certificate has bytes after its summary".into());
    }
    let raw_sha256 = hex_digest(raw.finalize());
    if raw_bytes != generation.uncompressed_bytes {
        return Err(format!(
            "certificate has {raw_bytes} uncompressed bytes, receipt says {}",
            generation.uncompressed_bytes
        ));
    }

    let canonical_file = input
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| format!("{} has no UTF-8 file name", input.display()))?
        .to_string();
    let manifest = CacheManifest {
        schema: CACHE_SCHEMA.into(),
        run_id: run_id.into(),
        result_class: "bounded-structural-screen-no-ecdlp-speedup".into(),
        source_commit: generation.source_commit,
        root_icv1_slug: header.root_icv1_slug,
        side: header.side,
        unique_curves: summary.unique_curves,
        chunk_rows,
        canonical_file,
        canonical: EncodedDigest {
            encoding: "gzip".into(),
            stored_sha256: canonical_stored.sha256,
            stored_bytes: canonical_stored.bytes,
            content_sha256: raw_sha256,
            content_bytes: raw_bytes,
        },
        generation_receipt: generation_digest,
        verification_receipt: verification_digest,
        summary,
        parts,
    };
    let mut manifest_bytes =
        serde_json::to_vec_pretty(&manifest).map_err(|error| error.to_string())?;
    manifest_bytes.push(b'\n');
    let manifest_path = partial.join(CACHE_MANIFEST_FILE);
    fs::write(&manifest_path, manifest_bytes)
        .map_err(|error| format!("write {}: {error}", manifest_path.display()))?;
    fs::rename(&partial, output).map_err(|error| {
        format!(
            "commit cache {} -> {}: {error}",
            partial.display(),
            output.display()
        )
    })?;
    Ok(manifest)
}

pub fn read_cache(path: &Path) -> Result<CacheManifest, String> {
    let manifest: CacheManifest = read_json(&path.join(CACHE_MANIFEST_FILE))?;
    if manifest.schema != CACHE_SCHEMA {
        return Err(format!("cache schema must be {CACHE_SCHEMA}"));
    }
    store::segment(&manifest.run_id)?;
    if manifest.chunk_rows == 0 || manifest.parts.is_empty() {
        return Err("cache manifest has no usable parts".into());
    }
    Ok(manifest)
}

pub fn object_key(sha256: &str) -> Result<String, String> {
    if sha256.len() != 64
        || !sha256
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err("object SHA-256 must be 64 lowercase hexadecimal characters".into());
    }
    Ok(format!("objects/sha256/{}/{}", &sha256[..2], sha256))
}

fn identity_digest(digest: &FileDigest) -> EncodedDigest {
    EncodedDigest {
        encoding: "identity".into(),
        stored_sha256: digest.sha256.clone(),
        stored_bytes: digest.bytes,
        content_sha256: digest.sha256.clone(),
        content_bytes: digest.bytes,
    }
}

fn stored_object(loc: &store::S3Loc, digest: &EncodedDigest) -> Result<StoredObject, String> {
    Ok(StoredObject {
        key: loc.key(&object_key(&digest.stored_sha256)?),
        encoding: digest.encoding.clone(),
        stored_sha256: digest.stored_sha256.clone(),
        stored_bytes: digest.stored_bytes,
        content_sha256: digest.content_sha256.clone(),
        content_bytes: digest.content_bytes,
    })
}

fn validate_encoded(path: &Path, expected: &EncodedDigest) -> Result<(), String> {
    let actual = match expected.encoding.as_str() {
        "gzip" => digest_gzip(path)?,
        "identity" => {
            let digest = digest_file(path)?;
            identity_digest(&digest)
        }
        other => return Err(format!("unsupported local encoding {other}")),
    };
    if &actual != expected {
        return Err(format!("{} differs from its cache digest", path.display()));
    }
    Ok(())
}

fn upload(loc: &store::S3Loc, path: &Path, digest: &EncodedDigest) -> Result<StoredObject, String> {
    validate_encoded(path, digest)?;
    let object = stored_object(loc, digest)?;
    store::put_file_once(
        loc,
        &object.key,
        path,
        &object.stored_sha256,
        object.stored_bytes,
    )?;
    Ok(object)
}

impl CompleteMarker {
    pub fn validate(&self, loc: &store::S3Loc) -> Result<(), String> {
        if self.schema != COMPLETE_SCHEMA {
            return Err(format!("completion schema must be {COMPLETE_SCHEMA}"));
        }
        store::segment(&self.run_id)?;
        if self.run != loc.uri(&format!("runs/{}", self.run_id)) {
            return Err("completion marker run URI disagrees with its location".into());
        }
        if self.result_class != "bounded-structural-screen-no-ecdlp-speedup" {
            return Err("completion marker overstates its result class".into());
        }
        if self.root_icv1_slug != P256_ICV1
            || self.source_commit.len() != 40
            || !self
                .source_commit
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit())
        {
            return Err("completion marker does not name registered P-256 source code".into());
        }
        if self.side == 0
            || self.unique_curves != u64::from(self.side) * u64::from(self.side)
            || self.summary.unique_curves != self.unique_curves
            || self.summary.parent_edges + 1 != self.unique_curves
            || self.summary.kernel_certificates != self.summary.parent_edges
            || self.summary.order_proved_prime != self.unique_curves
            || self.summary.non_singular != self.unique_curves
            || self.summary.generator_valid != self.unique_curves
            || self.summary.unique_j_invariants != self.unique_curves
            || self.summary.unique_canonical_models != self.unique_curves
            || self.summary.unique_icv1_identities != self.unique_curves
            || self.summary.unique_ec1_aliases != self.unique_curves
            || self.summary.unique_curve_uids != self.unique_curves
        {
            return Err("completion marker has inconsistent grid counts".into());
        }
        if self.chunk_rows == 0 || self.parts.len() < 3 {
            return Err("completion marker has no usable row slices".into());
        }
        let mut names = HashSet::new();
        let mut expected_row = 0u32;
        let mut records = 0u64;
        let mut content_bytes = 0u64;
        for (index, part) in self.parts.iter().enumerate() {
            if !names.insert(&part.name) {
                return Err(format!("duplicate cache part {}", part.name));
            }
            match (index, part.role.as_str(), &part.rows) {
                (0, "preamble", None) if part.records == u64::from(self.side) + 1 => {}
                (last, "summary", None) if last + 1 == self.parts.len() && part.records == 1 => {}
                (_, "rows", Some(rows))
                    if rows.start == expected_row
                        && rows.end_exclusive > rows.start
                        && rows.end_exclusive <= self.side
                        && rows.end_exclusive - rows.start <= self.chunk_rows
                        && part.records
                            == u64::from(rows.end_exclusive - rows.start)
                                * u64::from(self.side - 1) =>
                {
                    expected_row = rows.end_exclusive;
                }
                _ => return Err(format!("cache part {} breaks the row partition", part.name)),
            }
            records = records
                .checked_add(part.records)
                .ok_or("cache record count overflow")?;
            content_bytes = content_bytes
                .checked_add(part.object.content_bytes)
                .ok_or("cache content byte count overflow")?;
            validate_stored(loc, &part.object)?;
        }
        if expected_row != self.side
            || records != self.unique_curves + 2
            || content_bytes != self.canonical.content_bytes
        {
            return Err("cache parts do not exactly cover the canonical certificate".into());
        }
        validate_stored(loc, &self.canonical)?;
        validate_stored(loc, &self.generation_receipt)?;
        validate_stored(loc, &self.verification_receipt)?;
        Ok(())
    }

    pub fn cairn_receipt(
        &self,
        loc: &store::S3Loc,
        complete_key: &str,
        complete_sha256: &str,
    ) -> Result<CairnReceipt, String> {
        self.validate(loc)?;
        object_key(complete_sha256)?;
        let expected_key = loc.key(&format!("runs/{}/complete.json", self.run_id));
        if complete_key != expected_key {
            return Err("Cairn receipt complete key disagrees with the run".into());
        }
        Ok(CairnReceipt {
            schema: CAIRN_RECEIPT_SCHEMA.into(),
            claim_class: "coordination-receipt".into(),
            verification: "s3-content-address-only".into(),
            run_id: self.run_id.clone(),
            result_class: self.result_class.clone(),
            source_commit: self.source_commit.clone(),
            root_icv1_slug: self.root_icv1_slug.clone(),
            unique_curves: self.unique_curves,
            complete_key: complete_key.into(),
            complete_sha256: complete_sha256.into(),
            canonical_key: self.canonical.key.clone(),
            canonical_sha256: self.canonical.stored_sha256.clone(),
            canonical_content_sha256: self.canonical.content_sha256.clone(),
            record_chain_sha256: self.summary.record_chain_sha256.clone(),
        })
    }
}

fn validate_stored(loc: &store::S3Loc, object: &StoredObject) -> Result<(), String> {
    let expected = loc.key(&object_key(&object.stored_sha256)?);
    if object.key != expected {
        return Err(format!(
            "{} is not the object's content-addressed key",
            object.key
        ));
    }
    object_key(&object.content_sha256)?;
    if object.stored_bytes == 0 || object.content_bytes == 0 {
        return Err(format!("{} has an empty byte count", object.key));
    }
    match object.encoding.as_str() {
        "gzip" => {}
        "identity" => {
            if object.stored_sha256 != object.content_sha256
                || object.stored_bytes != object.content_bytes
            {
                return Err("identity object has different stored and content digests".into());
            }
        }
        other => return Err(format!("{} has unknown encoding {other}", object.key)),
    }
    Ok(())
}

/// Validate every local byte, upload content-addressed objects, and create the
/// immutable completion marker last. A retry reuses matching objects and the
/// same deterministic marker.
pub fn publish_cache(
    loc: &store::S3Loc,
    cache_dir: &Path,
    canonical_path: &Path,
    generation_path: &Path,
    verification_path: &Path,
    scratch: &Path,
) -> Result<Published, String> {
    let cache = read_cache(cache_dir)?;
    if canonical_path.file_name().and_then(|name| name.to_str())
        != Some(cache.canonical_file.as_str())
    {
        return Err("canonical artifact file name differs from the cache manifest".into());
    }
    let canonical = upload(loc, canonical_path, &cache.canonical)?;
    let generation = upload(
        loc,
        generation_path,
        &identity_digest(&cache.generation_receipt),
    )?;
    let verification = upload(
        loc,
        verification_path,
        &identity_digest(&cache.verification_receipt),
    )?;
    let mut parts = Vec::with_capacity(cache.parts.len());
    for part in &cache.parts {
        parts.push(StoredPart {
            name: part.name.clone(),
            role: part.role.clone(),
            rows: part.rows.clone(),
            records: part.records,
            object: upload(loc, &cache_dir.join(&part.file), &part.digest)?,
        });
    }
    let marker = CompleteMarker {
        schema: COMPLETE_SCHEMA.into(),
        run_id: cache.run_id.clone(),
        run: loc.uri(&format!("runs/{}", cache.run_id)),
        result_class: cache.result_class,
        source_commit: cache.source_commit,
        root_icv1_slug: cache.root_icv1_slug,
        side: cache.side,
        unique_curves: cache.unique_curves,
        chunk_rows: cache.chunk_rows,
        canonical,
        generation_receipt: generation,
        verification_receipt: verification,
        summary: cache.summary,
        parts,
    };
    marker.validate(loc)?;
    let mut bytes = serde_json::to_vec_pretty(&marker).map_err(|error| error.to_string())?;
    bytes.push(b'\n');
    let marker_sha256 = hex_digest(sha256(&bytes));
    let complete_key = loc.key(&format!("runs/{}/complete.json", marker.run_id));
    let created = match store::put_once(loc, &complete_key, &bytes, scratch)? {
        store::Put::Created => true,
        store::Put::Exists => {
            let existing = store::get(loc, &complete_key, scratch)?
                .ok_or_else(|| format!("{complete_key}: completion marker vanished"))?;
            if existing != bytes {
                return Err(format!(
                    "s3://{}/{complete_key}: a different completion marker already won",
                    loc.bucket
                ));
            }
            false
        }
    };
    let cairn_receipt = marker.cairn_receipt(loc, &complete_key, &marker_sha256)?;
    Ok(Published {
        marker,
        marker_uri: format!("s3://{}/{complete_key}", loc.bucket),
        marker_sha256,
        created,
        cairn_receipt,
    })
}

pub fn fetch_marker(
    loc: &store::S3Loc,
    run_id: &str,
    scratch: &Path,
) -> Result<(CompleteMarker, String, String), String> {
    store::segment(run_id)?;
    let key = loc.key(&format!("runs/{run_id}/complete.json"));
    let bytes = store::get(loc, &key, scratch)?
        .ok_or_else(|| format!("s3://{}/{key}: no completed run", loc.bucket))?;
    let marker: CompleteMarker =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    marker.validate(loc)?;
    if marker.run_id != run_id {
        return Err("completion marker run id differs from the requested run".into());
    }
    Ok((marker, key, hex_digest(sha256(&bytes))))
}

/// Fetch the canonical certificate (`part = None`) or one named cache part,
/// validate both stored and decompressed hashes, then atomically expose it.
pub fn fetch_object(
    loc: &store::S3Loc,
    run_id: &str,
    part: Option<&str>,
    output: &Path,
    scratch: &Path,
) -> Result<StoredObject, String> {
    let (marker, _, _) = fetch_marker(loc, run_id, scratch)?;
    let object = match part {
        None => marker.canonical.clone(),
        Some(name) => marker
            .parts
            .iter()
            .find(|candidate| candidate.name == name)
            .map(|candidate| candidate.object.clone())
            .ok_or_else(|| format!("run {run_id} has no cache part {name}"))?,
    };
    let expected = EncodedDigest {
        encoding: object.encoding.clone(),
        stored_sha256: object.stored_sha256.clone(),
        stored_bytes: object.stored_bytes,
        content_sha256: object.content_sha256.clone(),
        content_bytes: object.content_bytes,
    };
    if output.exists() {
        validate_encoded(output, &expected)?;
        return Ok(object);
    }
    let partial = partial_path(output);
    if partial.exists() {
        return Err(format!(
            "incomplete download already exists at {}",
            partial.display()
        ));
    }
    if !store::get_to_file(loc, &object.key, &partial)? {
        return Err(format!(
            "s3://{}/{}: object is missing",
            loc.bucket, object.key
        ));
    }
    validate_encoded(&partial, &expected).map_err(|error| {
        format!(
            "{error}; preserved the failed download at {}",
            partial.display()
        )
    })?;
    fs::rename(&partial, output).map_err(|error| {
        format!(
            "commit download {} -> {}: {error}",
            partial.display(),
            output.display()
        )
    })?;
    Ok(object)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::isogeny_walk::million::{generate_jsonl, verify_jsonl, GridConfig};
    use std::sync::atomic::{AtomicU64, Ordering};

    static SERIAL: AtomicU64 = AtomicU64::new(0);

    fn temporary() -> PathBuf {
        std::env::temp_dir().join(format!(
            "p256-million-store-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn streaming_digest_matches_one_shot() {
        let bytes: Vec<u8> = (0..2_000_003).map(|i| (i * 19) as u8).collect();
        let digest = digest_reader(bytes.as_slice()).unwrap();
        assert_eq!(digest.sha256, hex::encode(sha256(&bytes)));
        assert_eq!(digest.bytes, bytes.len() as u64);
    }

    #[test]
    fn content_keys_are_sharded_and_strict() {
        let digest = "01abcdef0123456789abcdef0123456789abcdef0123456789abcdef01234567";
        assert_eq!(
            object_key(digest).unwrap(),
            format!("objects/sha256/01/{digest}")
        );
        assert!(object_key("ABC").is_err());
        assert!(object_key(&"g".repeat(64)).is_err());
    }

    #[test]
    fn staged_parts_reconstruct_the_exact_verified_certificate() {
        std::thread::Builder::new()
            .stack_size(8 * 1024 * 1024)
            .spawn(staged_parts_reconstruct_the_exact_verified_certificate_inner)
            .unwrap()
            .join()
            .unwrap();
    }

    fn staged_parts_reconstruct_the_exact_verified_certificate_inner() {
        let root = temporary();
        fs::create_dir_all(&root).unwrap();
        let mut raw = Vec::new();
        let generated = generate_jsonl(
            &mut raw,
            &GridConfig {
                side: 3,
                batch_rows: 2,
                audit_points: 1,
                audit_seed_x: 0,
                source_commit: "0".repeat(40),
            },
        )
        .unwrap();
        let verified = verify_jsonl(BufReader::new(raw.as_slice()), 2, 7, 2).unwrap();
        let artifact = root.join("grid.jsonl.gz");
        let mut encoder = GzBuilder::new()
            .mtime(0)
            .operating_system(255)
            .write(Vec::new(), Compression::new(6));
        encoder.write_all(&raw).unwrap();
        fs::write(&artifact, encoder.finish().unwrap()).unwrap();
        let generation = root.join("GENERATE.json");
        let verification = root.join("VERIFY.json");
        fs::write(&generation, serde_json::to_vec_pretty(&generated).unwrap()).unwrap();
        fs::write(&verification, serde_json::to_vec_pretty(&verified).unwrap()).unwrap();
        let cache = root.join("cache");
        let manifest = stage_cache(
            &artifact,
            &generation,
            &verification,
            "p256-grid-test",
            &cache,
            2,
        )
        .unwrap();
        assert_eq!(manifest.parts.len(), 4);
        let mut rebuilt = Vec::new();
        for part in &manifest.parts {
            GzDecoder::new(File::open(cache.join(&part.file)).unwrap())
                .read_to_end(&mut rebuilt)
                .unwrap();
        }
        assert_eq!(rebuilt, raw);
        assert_eq!(manifest.canonical.content_sha256, hex::encode(sha256(&raw)));
        let loc = store::S3Loc::parse("s3://example/p256").unwrap();
        let marker = CompleteMarker {
            schema: COMPLETE_SCHEMA.into(),
            run_id: manifest.run_id.clone(),
            run: loc.uri(&format!("runs/{}", manifest.run_id)),
            result_class: manifest.result_class.clone(),
            source_commit: manifest.source_commit.clone(),
            root_icv1_slug: manifest.root_icv1_slug.clone(),
            side: manifest.side,
            unique_curves: manifest.unique_curves,
            chunk_rows: manifest.chunk_rows,
            canonical: stored_object(&loc, &manifest.canonical).unwrap(),
            generation_receipt: stored_object(&loc, &identity_digest(&manifest.generation_receipt))
                .unwrap(),
            verification_receipt: stored_object(
                &loc,
                &identity_digest(&manifest.verification_receipt),
            )
            .unwrap(),
            summary: manifest.summary.clone(),
            parts: manifest
                .parts
                .iter()
                .map(|part| StoredPart {
                    name: part.name.clone(),
                    role: part.role.clone(),
                    rows: part.rows.clone(),
                    records: part.records,
                    object: stored_object(&loc, &part.digest).unwrap(),
                })
                .collect(),
        };
        marker.validate(&loc).unwrap();
        let complete_key = loc.key("runs/p256-grid-test/complete.json");
        let receipt = marker
            .cairn_receipt(&loc, &complete_key, &"a".repeat(64))
            .unwrap();
        assert_eq!(receipt.claim_class, "coordination-receipt");
        assert_eq!(receipt.verification, "s3-content-address-only");
        assert_eq!(receipt.claim_key(), "p256-grid-test");
        let mut bad_marker = marker;
        bad_marker.parts[1].rows.as_mut().unwrap().start = 1;
        assert!(bad_marker.validate(&loc).is_err());
        let first = cache.join(&manifest.parts[0].file);
        let mut damaged = fs::read(&first).unwrap();
        damaged.push(0);
        fs::write(&first, damaged).unwrap();
        assert!(validate_encoded(&first, &manifest.parts[0].digest).is_err());
        fs::remove_dir_all(root).unwrap();
    }
}
