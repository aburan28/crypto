//! External-memory ledger and uniqueness audit for p61 class-census segments.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, ErrorKind, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::isogeny_walk::million_store::{digest_file, FileDigest};
use crypto_lib::cryptanalysis::isogeny_walk::p61_class::{
    coordinate_index, Checkpoint, ClosureRecord, GenerationReceipt, NodeRecord, SegmentHeader,
    SegmentSummary, VerificationReceipt, CHECKPOINT_SCHEMA, COSET_ROWS, EXPECTED_CLASS_VERTICES,
    LANES, RECEIPT_SCHEMA, ROW_ORDER, SEGMENT_SCHEMA, SUMMARY_SCHEMA, TOP_ROWS,
};
use crypto_lib::hash::sha256::{sha256, Sha256};
use flate2::read::GzDecoder;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use serde_json::json;

const LEDGER_SCHEMA: &str = "ec-hardness/p61-class-ledger/v1";
const AUDIT_SCHEMA: &str = "ec-hardness/p61-class-uniqueness-audit/v1";
const INDEX_MAGIC: &[u8; 8] = b"P61JIDX1";
const FIELD_P: u64 = 2_305_843_009_213_693_951;
const MAX_DUPLICATE_EXAMPLES: usize = 32;

#[derive(Debug, Parser)]
#[command(about = "Build and externally audit a p61 class-census ledger")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Build an immutable ordered ledger from generation and verification receipts.
    BuildLedger {
        #[arg(long)]
        generation_receipts_dir: PathBuf,
        #[arg(long)]
        verification_receipts_dir: PathBuf,
        #[arg(long)]
        class_number_receipt: PathBuf,
        #[arg(long)]
        class_certification_state: String,
        #[arg(long)]
        auditor_source_commit: String,
        #[arg(long)]
        output: PathBuf,
    },
    /// Extract every j-invariant, externally sort it, and prove uniqueness or retain duplicates.
    AuditLedger {
        #[arg(long)]
        ledger: PathBuf,
        #[arg(long)]
        work_dir: PathBuf,
        #[arg(long)]
        sorted_index_output: PathBuf,
        #[arg(long)]
        receipt_output: PathBuf,
        #[arg(long, default_value_t = 5_000_000)]
        chunk_records: usize,
        #[arg(long)]
        auditor_source_commit: String,
    },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct LedgerEntry {
    x_start: u64,
    x_end: u64,
    node_records: u64,
    segment: PathBuf,
    segment_digest: FileDigest,
    checkpoint: PathBuf,
    checkpoint_digest: FileDigest,
    generation_receipt: PathBuf,
    generation_receipt_digest: FileDigest,
    verification_receipt: PathBuf,
    verification_receipt_digest: FileDigest,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Ledger {
    schema: String,
    status: String,
    builder_source_commit: String,
    source_commit: String,
    class_number_receipt: PathBuf,
    class_number_receipt_digest: FileDigest,
    class_certification_state: String,
    x_start: u64,
    x_end: u64,
    expected_node_records: u64,
    expected_complete_class_vertices: u64,
    entries: Vec<LedgerEntry>,
}

#[derive(Debug, Deserialize)]
struct GenerationEnvelope {
    status: String,
    segment: PathBuf,
    segment_digest: FileDigest,
    checkpoint: PathBuf,
    checkpoint_digest: FileDigest,
    receipt: GenerationReceipt,
}

#[derive(Debug, Deserialize)]
struct VerificationEnvelope {
    status: String,
    segment: PathBuf,
    segment_digest: FileDigest,
    checkpoint: PathBuf,
    checkpoint_digest: FileDigest,
    receipt: VerificationReceipt,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct IndexEntry {
    j: u64,
    index: u64,
}

#[derive(Clone, Debug, Serialize)]
struct DurableWrite {
    digest: FileDigest,
    sync_status: String,
}

#[derive(Clone, Debug, Serialize)]
struct ChunkMeta {
    path: PathBuf,
    records: u64,
    min_j: String,
    max_j: String,
    digest: FileDigest,
    sync_status: String,
}

#[derive(Clone, Debug, Serialize)]
struct SegmentAudit {
    x_start: u64,
    x_end: u64,
    node_records: u64,
    closure_records: u64,
    record_chain_sha256: String,
    segment_digest: FileDigest,
    checkpoint_digest: FileDigest,
    generation_receipt_digest: FileDigest,
    verification_receipt_digest: FileDigest,
}

#[derive(Clone, Debug, Serialize)]
struct DuplicateExample {
    j: String,
    first_index: u64,
    duplicate_index: u64,
}

#[derive(Clone, Debug, Serialize)]
struct AuditReceipt {
    schema: String,
    status: String,
    auditor_source_commit: String,
    ledger_builder_source_commit: String,
    source_commit: String,
    class_number_receipt_sha256: String,
    class_certification_state: String,
    ledger: PathBuf,
    ledger_digest: FileDigest,
    x_start: u64,
    x_end: u64,
    records: u64,
    unique_j_invariants: u64,
    duplicate_records: u64,
    duplicate_groups: u64,
    duplicate_examples: Vec<DuplicateExample>,
    complete_class_coverage: bool,
    chunks: Vec<ChunkMeta>,
    sorted_index: PathBuf,
    sorted_index_digest: FileDigest,
    sorted_index_sync_status: String,
    segments: Vec<SegmentAudit>,
    elapsed_ms: u128,
}

fn partial_path(path: &Path) -> PathBuf {
    let mut name: OsString = path.as_os_str().to_owned();
    name.push(".partial");
    PathBuf::from(name)
}

fn is_lower_hex(value: &str, length: usize) -> bool {
    value.len() == length
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn sync_file(file: &File, path: &Path) -> Result<String, String> {
    match file.sync_all() {
        Ok(()) => Ok("sync_all_completed".into()),
        Err(error)
            if error.kind() == ErrorKind::Unsupported || error.raw_os_error() == Some(45) =>
        {
            Ok("sync_all_unsupported_by_filesystem".into())
        }
        Err(error) => Err(format!("sync {}: {error}", path.display())),
    }
}

fn read_json<T: DeserializeOwned>(path: &Path) -> Result<T, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    serde_json::from_slice(&bytes).map_err(|error| format!("parse {}: {error}", path.display()))
}

fn write_json_new<T: Serialize>(path: &Path, value: &T) -> Result<DurableWrite, String> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    if let Some(parent) = path.parent().filter(|value| !value.as_os_str().is_empty()) {
        fs::create_dir_all(parent)
            .map_err(|error| format!("create {}: {error}", parent.display()))?;
    }
    let partial = partial_path(path);
    let file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&partial)
        .map_err(|error| format!("create {}: {error}", partial.display()))?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, value).map_err(|error| error.to_string())?;
    writer.write_all(b"\n").map_err(|error| error.to_string())?;
    writer.flush().map_err(|error| error.to_string())?;
    let sync_status = sync_file(writer.get_ref(), &partial)?;
    fs::rename(&partial, path).map_err(|error| {
        format!(
            "commit {} -> {}: {error}",
            partial.display(),
            path.display()
        )
    })?;
    Ok(DurableWrite {
        digest: digest_file(path)?,
        sync_status,
    })
}

fn json_files(directory: &Path) -> Result<Vec<PathBuf>, String> {
    let mut paths = fs::read_dir(directory)
        .map_err(|error| format!("read {}: {error}", directory.display()))?
        .map(|entry| {
            entry
                .map(|value| value.path())
                .map_err(|error| format!("read {} entry: {error}", directory.display()))
        })
        .collect::<Result<Vec<_>, _>>()?;
    paths.retain(|path| {
        path.extension()
            .is_some_and(|extension| extension == "json")
    });
    paths.sort();
    Ok(paths)
}

fn ensure_digest(path: &Path, expected: &FileDigest, what: &str) -> Result<FileDigest, String> {
    let actual = digest_file(path)?;
    if &actual != expected {
        return Err(format!(
            "{what} digest differs for {}: expected {:?}, observed {:?}",
            path.display(),
            expected,
            actual
        ));
    }
    Ok(actual)
}

fn validate_generation_envelope(path: &Path, envelope: &GenerationEnvelope) -> Result<(), String> {
    if envelope.status != "pass"
        || envelope.receipt.status != "pass"
        || envelope.receipt.schema != RECEIPT_SCHEMA
        || envelope.receipt.checkpoint.schema != CHECKPOINT_SCHEMA
        || envelope.receipt.summary.schema != SUMMARY_SCHEMA
        || envelope.receipt.summary.record != "summary"
        || envelope.receipt.checkpoint.segment_sha256 != envelope.segment_digest.sha256
    {
        return Err(format!(
            "generation receipt {} is not an admitted p61 pass",
            path.display()
        ));
    }
    ensure_digest(&envelope.segment, &envelope.segment_digest, "segment")?;
    ensure_digest(
        &envelope.checkpoint,
        &envelope.checkpoint_digest,
        "checkpoint",
    )?;
    Ok(())
}

fn validate_verification_envelope(
    path: &Path,
    envelope: &VerificationEnvelope,
) -> Result<(), String> {
    if envelope.status != "pass"
        || envelope.receipt.status != "pass"
        || envelope.receipt.schema != RECEIPT_SCHEMA
        || envelope.receipt.summary.schema != SUMMARY_SCHEMA
        || envelope.receipt.summary.record != "summary"
    {
        return Err(format!(
            "verification receipt {} is not an admitted p61 pass",
            path.display()
        ));
    }
    ensure_digest(
        &envelope.segment,
        &envelope.segment_digest,
        "verified segment",
    )?;
    ensure_digest(
        &envelope.checkpoint,
        &envelope.checkpoint_digest,
        "verified checkpoint",
    )?;
    Ok(())
}

fn build_ledger(
    generation_dir: &Path,
    verification_dir: &Path,
    class_number_receipt: &Path,
    class_certification_state: String,
    auditor_source_commit: String,
    output: &Path,
) -> Result<(), String> {
    if class_certification_state.trim().is_empty() {
        return Err("class certification state must be explicit".into());
    }
    if !is_lower_hex(&auditor_source_commit, 40) {
        return Err("auditor source commit must be 40 lowercase hexadecimal digits".into());
    }
    let class_number_receipt_digest = digest_file(class_number_receipt)?;
    let class_number_value: serde_json::Value = read_json(class_number_receipt)?;
    let recorded_certification_state = class_number_value
        .pointer("/selection/class_number_status")
        .and_then(|value| value.as_str())
        .ok_or("class-number receipt lacks selection.class_number_status")?;
    let recorded_vertices = class_number_value
        .pointer("/selection/expected_vertices")
        .and_then(|value| value.as_u64())
        .ok_or("class-number receipt lacks selection.expected_vertices")?;
    if recorded_certification_state != class_certification_state
        || recorded_vertices != EXPECTED_CLASS_VERTICES
    {
        return Err(format!(
            "declared class certification state or vertex count differs from {}",
            class_number_receipt.display()
        ));
    }
    let verification_paths = json_files(verification_dir)?;
    let mut verifications = HashMap::new();
    for path in verification_paths {
        let envelope: VerificationEnvelope = read_json(&path)?;
        validate_verification_envelope(&path, &envelope)?;
        if verifications
            .insert(envelope.segment_digest.sha256.clone(), (path, envelope))
            .is_some()
        {
            return Err("multiple verification receipts bind the same segment".into());
        }
    }

    let mut entries = Vec::new();
    let mut source_commit: Option<String> = None;
    for generation_path in json_files(generation_dir)? {
        let envelope: GenerationEnvelope = read_json(&generation_path)?;
        validate_generation_envelope(&generation_path, &envelope)?;
        let verification = verifications
            .remove(&envelope.segment_digest.sha256)
            .ok_or_else(|| {
                format!(
                    "segment {} has no matching verification receipt",
                    envelope.segment.display()
                )
            })?;
        let (verification_path, verification_envelope) = verification;
        if verification_envelope.segment != envelope.segment
            || verification_envelope.segment_digest != envelope.segment_digest
            || verification_envelope.checkpoint != envelope.checkpoint
            || verification_envelope.checkpoint_digest != envelope.checkpoint_digest
            || verification_envelope.receipt.source_commit != envelope.receipt.source_commit
            || verification_envelope.receipt.summary != envelope.receipt.summary
        {
            return Err(format!(
                "generation and verification receipts disagree for {}",
                envelope.segment.display()
            ));
        }
        if envelope.receipt.checkpoint.class_number_receipt_sha256
            != class_number_receipt_digest.sha256
        {
            return Err(format!(
                "{} is bound to a different class-number receipt",
                generation_path.display()
            ));
        }
        match &source_commit {
            None => source_commit = Some(envelope.receipt.source_commit.clone()),
            Some(expected) if expected != &envelope.receipt.source_commit => {
                return Err("generation receipts mix source commits".into())
            }
            _ => {}
        }
        entries.push(LedgerEntry {
            x_start: envelope.receipt.summary.x_start,
            x_end: envelope.receipt.summary.x_end,
            node_records: envelope.receipt.summary.node_records,
            segment: envelope.segment,
            segment_digest: envelope.segment_digest,
            checkpoint: envelope.checkpoint,
            checkpoint_digest: envelope.checkpoint_digest,
            generation_receipt_digest: digest_file(&generation_path)?,
            generation_receipt: generation_path,
            verification_receipt_digest: digest_file(&verification_path)?,
            verification_receipt: verification_path,
        });
    }
    if entries.is_empty() {
        return Err("generation receipt directory contains no admitted entries".into());
    }
    if !verifications.is_empty() {
        return Err(
            "verification directory contains receipts absent from generation ledger".into(),
        );
    }
    entries.sort_by_key(|entry| (entry.x_start, entry.x_end));
    let x_start = entries[0].x_start;
    let mut next_x = x_start;
    let mut expected_node_records = 0u64;
    let mut prior_checkpoint_digest: Option<String> = None;
    for entry in &entries {
        if entry.x_start != next_x || entry.x_end <= entry.x_start || entry.x_end > ROW_ORDER {
            return Err(format!(
                "ledger is not a positive contiguous partition at {}..{} (expected start {})",
                entry.x_start, entry.x_end, next_x
            ));
        }
        let expected = (entry.x_end - entry.x_start)
            .checked_mul(LANES as u64)
            .ok_or("node count overflow")?;
        if entry.node_records != expected {
            return Err(format!(
                "segment {}..{} has {} nodes, expected {}",
                entry.x_start, entry.x_end, entry.node_records, expected
            ));
        }
        let checkpoint: Checkpoint = read_json(&entry.checkpoint)?;
        if checkpoint.source_commit != source_commit.as_deref().unwrap_or_default()
            || checkpoint.class_number_receipt_sha256 != class_number_receipt_digest.sha256
            || checkpoint.x_start != entry.x_start
            || checkpoint.next_x != entry.x_end
            || checkpoint.segment_sha256 != entry.segment_digest.sha256
            || checkpoint.previous_checkpoint_sha256 != prior_checkpoint_digest
        {
            return Err(format!(
                "checkpoint continuity differs at {}..{}",
                entry.x_start, entry.x_end
            ));
        }
        prior_checkpoint_digest = Some(entry.checkpoint_digest.sha256.clone());
        expected_node_records = expected_node_records
            .checked_add(expected)
            .ok_or("ledger node count overflow")?;
        next_x = entry.x_end;
    }
    let ledger = Ledger {
        schema: LEDGER_SCHEMA.into(),
        status: "sealed".into(),
        builder_source_commit: auditor_source_commit,
        source_commit: source_commit.ok_or("missing source commit")?,
        class_number_receipt: class_number_receipt.to_path_buf(),
        class_number_receipt_digest,
        class_certification_state,
        x_start,
        x_end: next_x,
        expected_node_records,
        expected_complete_class_vertices: EXPECTED_CLASS_VERTICES,
        entries,
    };
    let written = write_json_new(output, &ledger)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "status": "pass",
            "ledger": output,
            "ledger_digest": written.digest,
            "ledger_sync_status": written.sync_status,
            "x_start": ledger.x_start,
            "x_end": ledger.x_end,
            "entries": ledger.entries.len(),
            "expected_node_records": ledger.expected_node_records,
        }))
        .map_err(|error| error.to_string())?
    );
    Ok(())
}

fn update_record_chain(chain: &mut [u8; 32], line: &[u8]) {
    let mut input = Vec::with_capacity(32 + line.len());
    input.extend_from_slice(chain);
    input.extend_from_slice(line);
    *chain = sha256(&input);
}

fn parse_j(value: &str) -> Result<u64, String> {
    if value.len() != 16
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err("j-invariant is not sixteen lowercase hexadecimal digits".into());
    }
    let parsed = u64::from_str_radix(value, 16).map_err(|error| error.to_string())?;
    if parsed >= FIELD_P {
        return Err("j-invariant is not a canonical field element".into());
    }
    Ok(parsed)
}

fn local_ordinal(node: &NodeRecord, x_start: u64, x_end: u64) -> Result<usize, String> {
    if node.coordinate.x < x_start
        || node.coordinate.x >= x_end
        || node.coordinate.y >= COSET_ROWS
        || node.coordinate.z >= TOP_ROWS
    {
        return Err(format!(
            "node {} coordinate lies outside segment {}..{}",
            node.index, x_start, x_end
        ));
    }
    if coordinate_index(node.coordinate)? != node.index {
        return Err(format!(
            "node {} coordinate/index mapping differs",
            node.index
        ));
    }
    let width = usize::try_from(x_end - x_start).map_err(|error| error.to_string())?;
    let lane =
        usize::from(node.coordinate.z) * usize::from(COSET_ROWS) + usize::from(node.coordinate.y);
    let x = usize::try_from(node.coordinate.x - x_start).map_err(|error| error.to_string())?;
    lane.checked_mul(width)
        .and_then(|value| value.checked_add(x))
        .ok_or("local node ordinal overflow".into())
}

fn set_once(bits: &mut [u64], ordinal: usize) -> Result<(), String> {
    let word = bits
        .get_mut(ordinal / 64)
        .ok_or("local node ordinal lies outside bitmap")?;
    let mask = 1u64 << (ordinal % 64);
    if *word & mask != 0 {
        return Err(format!("duplicate local node ordinal {ordinal}"));
    }
    *word |= mask;
    Ok(())
}

fn read_line_strict<R: BufRead>(
    reader: &mut R,
    line_number: u64,
) -> Result<Option<Vec<u8>>, String> {
    let mut line = Vec::new();
    let read = reader
        .read_until(b'\n', &mut line)
        .map_err(|error| format!("line {line_number}: {error}"))?;
    if read == 0 {
        return Ok(None);
    }
    if line.last() != Some(&b'\n') || line.ends_with(b"\r\n") {
        return Err(format!("line {line_number}: records require one LF"));
    }
    if line.len() > 4 * 1024 * 1024 {
        return Err(format!("line {line_number}: record exceeds 4 MiB"));
    }
    Ok(Some(line))
}

fn flush_chunk(
    entries: &mut Vec<IndexEntry>,
    work_dir: &Path,
    chunks: &mut Vec<ChunkMeta>,
) -> Result<(), String> {
    if entries.is_empty() {
        return Ok(());
    }
    entries.sort_unstable();
    let path = work_dir.join(format!("j-index-chunk-{:06}.bin", chunks.len()));
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    let partial = partial_path(&path);
    let file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&partial)
        .map_err(|error| format!("create {}: {error}", partial.display()))?;
    let mut writer = BufWriter::new(file);
    let mut hasher = Sha256::new();
    let mut bytes = 0u64;
    for entry in entries.iter() {
        let mut encoded = [0u8; 16];
        encoded[..8].copy_from_slice(&entry.j.to_le_bytes());
        encoded[8..].copy_from_slice(&entry.index.to_le_bytes());
        writer
            .write_all(&encoded)
            .map_err(|error| error.to_string())?;
        hasher.update(&encoded);
        bytes += encoded.len() as u64;
    }
    writer.flush().map_err(|error| error.to_string())?;
    let sync_status = sync_file(writer.get_ref(), &partial)?;
    fs::rename(&partial, &path).map_err(|error| {
        format!(
            "commit {} -> {}: {error}",
            partial.display(),
            path.display()
        )
    })?;
    let digest = FileDigest {
        sha256: hex::encode(hasher.finalize()),
        bytes,
    };
    ensure_digest(&path, &digest, "sorted chunk")?;
    chunks.push(ChunkMeta {
        path,
        records: entries.len() as u64,
        min_j: format!("{:016x}", entries.first().expect("nonempty chunk").j),
        max_j: format!("{:016x}", entries.last().expect("nonempty chunk").j),
        digest,
        sync_status,
    });
    entries.clear();
    Ok(())
}

fn validate_entry_artifacts(
    ledger: &Ledger,
    entry: &LedgerEntry,
    prior_checkpoint_digest: Option<&str>,
) -> Result<Checkpoint, String> {
    ensure_digest(&entry.segment, &entry.segment_digest, "ledger segment")?;
    ensure_digest(
        &entry.checkpoint,
        &entry.checkpoint_digest,
        "ledger checkpoint",
    )?;
    ensure_digest(
        &entry.generation_receipt,
        &entry.generation_receipt_digest,
        "generation receipt",
    )?;
    ensure_digest(
        &entry.verification_receipt,
        &entry.verification_receipt_digest,
        "verification receipt",
    )?;
    let generation: GenerationEnvelope = read_json(&entry.generation_receipt)?;
    let verification: VerificationEnvelope = read_json(&entry.verification_receipt)?;
    validate_generation_envelope(&entry.generation_receipt, &generation)?;
    validate_verification_envelope(&entry.verification_receipt, &verification)?;
    if generation.segment != entry.segment
        || generation.segment_digest != entry.segment_digest
        || generation.checkpoint != entry.checkpoint
        || generation.checkpoint_digest != entry.checkpoint_digest
        || verification.segment != entry.segment
        || verification.segment_digest != entry.segment_digest
        || verification.checkpoint != entry.checkpoint
        || verification.checkpoint_digest != entry.checkpoint_digest
        || generation.receipt.source_commit != ledger.source_commit
        || verification.receipt.source_commit != ledger.source_commit
        || generation.receipt.summary != verification.receipt.summary
        || generation.receipt.summary.x_start != entry.x_start
        || generation.receipt.summary.x_end != entry.x_end
        || generation.receipt.summary.node_records != entry.node_records
    {
        return Err(format!(
            "ledger artifacts disagree at {}..{}",
            entry.x_start, entry.x_end
        ));
    }
    let checkpoint: Checkpoint = read_json(&entry.checkpoint)?;
    if checkpoint.schema != CHECKPOINT_SCHEMA
        || checkpoint.source_commit != ledger.source_commit
        || checkpoint.class_number_receipt_sha256 != ledger.class_number_receipt_digest.sha256
        || checkpoint.x_start != entry.x_start
        || checkpoint.next_x != entry.x_end
        || checkpoint.segment_sha256 != entry.segment_digest.sha256
        || checkpoint.previous_checkpoint_sha256.as_deref() != prior_checkpoint_digest
    {
        return Err(format!(
            "checkpoint chain differs at {}..{}",
            entry.x_start, entry.x_end
        ));
    }
    Ok(checkpoint)
}

fn audit_segment(
    ledger: &Ledger,
    entry: &LedgerEntry,
    buffer: &mut Vec<IndexEntry>,
    chunk_records: usize,
    work_dir: &Path,
    chunks: &mut Vec<ChunkMeta>,
) -> Result<SegmentAudit, String> {
    let file = File::open(&entry.segment)
        .map_err(|error| format!("open {}: {error}", entry.segment.display()))?;
    let mut reader = BufReader::new(GzDecoder::new(file));
    let width = usize::try_from(entry.x_end - entry.x_start).map_err(|error| error.to_string())?;
    let expected_nodes = width
        .checked_mul(LANES)
        .ok_or("segment bitmap size overflow")?;
    let mut bits = vec![0u64; expected_nodes.div_ceil(64)];
    let mut line_number = 0u64;
    let mut header: Option<SegmentHeader> = None;
    let mut summary: Option<SegmentSummary> = None;
    let mut nodes = 0u64;
    let mut closures = 0u64;
    let mut record_chain = [0u8; 32];

    loop {
        line_number += 1;
        let Some(line) = read_line_strict(&mut reader, line_number)? else {
            break;
        };
        if summary.is_some() {
            return Err(format!(
                "{} has records after its summary",
                entry.segment.display()
            ));
        }
        let value: serde_json::Value = serde_json::from_slice(&line)
            .map_err(|error| format!("{} line {line_number}: {error}", entry.segment.display()))?;
        let record = value
            .get("record")
            .and_then(|value| value.as_str())
            .ok_or_else(|| {
                format!(
                    "{} line {line_number}: missing record tag",
                    entry.segment.display()
                )
            })?;
        match record {
            "header" => {
                if header.is_some() || line_number != 1 {
                    return Err("segment header is not the unique first record".into());
                }
                let parsed: SegmentHeader =
                    serde_json::from_value(value).map_err(|error| error.to_string())?;
                if parsed.schema != SEGMENT_SCHEMA
                    || parsed.source_commit != ledger.source_commit
                    || parsed.class_number_receipt_sha256
                        != ledger.class_number_receipt_digest.sha256
                    || parsed.x_start != entry.x_start
                    || parsed.x_end != entry.x_end
                    || parsed.expected_class_vertices != EXPECTED_CLASS_VERTICES
                    || parsed.row_order != ROW_ORDER
                    || usize::from(parsed.coset_rows) * usize::from(parsed.top_rows) != LANES
                {
                    return Err(format!(
                        "{} header differs from the ledger and frozen p61 lattice",
                        entry.segment.display()
                    ));
                }
                update_record_chain(&mut record_chain, &line);
                header = Some(parsed);
            }
            "node" => {
                if header.is_none() {
                    return Err("node precedes segment header".into());
                }
                let node: NodeRecord =
                    serde_json::from_value(value).map_err(|error| error.to_string())?;
                let ordinal = local_ordinal(&node, entry.x_start, entry.x_end)?;
                set_once(&mut bits, ordinal)?;
                buffer.push(IndexEntry {
                    j: parse_j(&node.curve.j)?,
                    index: node.index,
                });
                if buffer.len() >= chunk_records {
                    flush_chunk(buffer, work_dir, chunks)?;
                }
                nodes += 1;
                update_record_chain(&mut record_chain, &line);
            }
            "closure" => {
                if header.is_none() || entry.x_end != ROW_ORDER {
                    return Err("closure appears outside the final row segment".into());
                }
                let _: ClosureRecord =
                    serde_json::from_value(value).map_err(|error| error.to_string())?;
                closures += 1;
                update_record_chain(&mut record_chain, &line);
            }
            "summary" => {
                let parsed: SegmentSummary =
                    serde_json::from_value(value).map_err(|error| error.to_string())?;
                summary = Some(parsed);
            }
            other => return Err(format!("unexpected segment record {other}")),
        }
    }

    let summary = summary.ok_or("segment lacks a summary")?;
    if nodes != expected_nodes as u64
        || nodes != entry.node_records
        || bits
            .iter()
            .map(|word| word.count_ones() as u64)
            .sum::<u64>()
            != nodes
        || summary.schema != SUMMARY_SCHEMA
        || summary.x_start != entry.x_start
        || summary.x_end != entry.x_end
        || summary.node_records != nodes
        || summary.closure_edges != closures
        || summary.record_chain_sha256 != hex::encode(record_chain)
    {
        return Err(format!(
            "{} coverage, summary, or record chain differs",
            entry.segment.display()
        ));
    }
    if (entry.x_end == ROW_ORDER && closures != LANES as u64)
        || (entry.x_end != ROW_ORDER && closures != 0)
    {
        return Err("row closure count differs from the frozen lattice".into());
    }
    Ok(SegmentAudit {
        x_start: entry.x_start,
        x_end: entry.x_end,
        node_records: nodes,
        closure_records: closures,
        record_chain_sha256: summary.record_chain_sha256,
        segment_digest: entry.segment_digest.clone(),
        checkpoint_digest: entry.checkpoint_digest.clone(),
        generation_receipt_digest: entry.generation_receipt_digest.clone(),
        verification_receipt_digest: entry.verification_receipt_digest.clone(),
    })
}

struct ChunkReader {
    reader: BufReader<File>,
    expected: u64,
    read: u64,
}

impl ChunkReader {
    fn open(meta: &ChunkMeta) -> Result<Self, String> {
        ensure_digest(&meta.path, &meta.digest, "merge chunk")?;
        let file = File::open(&meta.path)
            .map_err(|error| format!("open {}: {error}", meta.path.display()))?;
        Ok(Self {
            reader: BufReader::new(file),
            expected: meta.records,
            read: 0,
        })
    }

    fn next(&mut self) -> Result<Option<IndexEntry>, String> {
        let mut encoded = [0u8; 16];
        let first = self
            .reader
            .read(&mut encoded)
            .map_err(|error| format!("read sorted chunk: {error}"))?;
        if first == 0 {
            if self.read != self.expected {
                return Err(format!(
                    "sorted chunk ended after {} records, expected {}",
                    self.read, self.expected
                ));
            }
            return Ok(None);
        }
        self.reader
            .read_exact(&mut encoded[first..])
            .map_err(|error| format!("truncated sorted chunk record: {error}"))?;
        self.read += 1;
        if self.read > self.expected {
            return Err("sorted chunk contains more records than its manifest".into());
        }
        Ok(Some(IndexEntry {
            j: u64::from_le_bytes(encoded[..8].try_into().expect("eight bytes")),
            index: u64::from_le_bytes(encoded[8..].try_into().expect("eight bytes")),
        }))
    }
}

struct MergeResult {
    unique_j: u64,
    duplicate_records: u64,
    duplicate_groups: u64,
    duplicate_examples: Vec<DuplicateExample>,
    sorted_index: DurableWrite,
}

fn merge_chunks(
    chunks: &[ChunkMeta],
    output: &Path,
    total_records: u64,
) -> Result<MergeResult, String> {
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()));
    }
    if let Some(parent) = output
        .parent()
        .filter(|value| !value.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .map_err(|error| format!("create {}: {error}", parent.display()))?;
    }
    let partial = partial_path(output);
    let file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&partial)
        .map_err(|error| format!("create {}: {error}", partial.display()))?;
    let mut writer = BufWriter::new(file);
    let mut hasher = Sha256::new();
    let mut header = [0u8; 16];
    header[..8].copy_from_slice(INDEX_MAGIC);
    header[8..].copy_from_slice(&total_records.to_le_bytes());
    writer
        .write_all(&header)
        .map_err(|error| error.to_string())?;
    hasher.update(&header);

    let mut readers = chunks
        .iter()
        .map(ChunkReader::open)
        .collect::<Result<Vec<_>, _>>()?;
    let mut heap = BinaryHeap::new();
    for (index, reader) in readers.iter_mut().enumerate() {
        if let Some(entry) = reader.next()? {
            heap.push(Reverse((entry, index)));
        }
    }
    let mut observed = 0u64;
    let mut unique_j = 0u64;
    let mut duplicate_records = 0u64;
    let mut duplicate_groups = 0u64;
    let mut duplicate_examples = Vec::new();
    let mut previous: Option<IndexEntry> = None;
    let mut in_duplicate_group = false;
    let mut duplicate_group_first: Option<IndexEntry> = None;
    while let Some(Reverse((entry, source))) = heap.pop() {
        let mut encoded = [0u8; 16];
        encoded[..8].copy_from_slice(&entry.j.to_le_bytes());
        encoded[8..].copy_from_slice(&entry.index.to_le_bytes());
        writer
            .write_all(&encoded)
            .map_err(|error| error.to_string())?;
        hasher.update(&encoded);
        observed += 1;
        match previous {
            Some(prior) if prior.j == entry.j => {
                duplicate_records += 1;
                if !in_duplicate_group {
                    duplicate_groups += 1;
                    in_duplicate_group = true;
                    duplicate_group_first = Some(prior);
                }
                if duplicate_examples.len() < MAX_DUPLICATE_EXAMPLES {
                    duplicate_examples.push(DuplicateExample {
                        j: format!("{:016x}", entry.j),
                        first_index: duplicate_group_first
                            .expect("duplicate group has a first record")
                            .index,
                        duplicate_index: entry.index,
                    });
                }
            }
            _ => {
                unique_j += 1;
                in_duplicate_group = false;
                duplicate_group_first = None;
            }
        }
        previous = Some(entry);
        if let Some(next) = readers[source].next()? {
            heap.push(Reverse((next, source)));
        }
    }
    if observed != total_records {
        return Err(format!(
            "merge produced {observed} records, expected {total_records}"
        ));
    }
    writer.flush().map_err(|error| error.to_string())?;
    let sync_status = sync_file(writer.get_ref(), &partial)?;
    fs::rename(&partial, output).map_err(|error| {
        format!(
            "commit {} -> {}: {error}",
            partial.display(),
            output.display()
        )
    })?;
    let digest = FileDigest {
        sha256: hex::encode(hasher.finalize()),
        bytes: 16u64
            .checked_add(total_records.checked_mul(16).ok_or("index size overflow")?)
            .ok_or("index size overflow")?,
    };
    ensure_digest(output, &digest, "merged j-index")?;
    Ok(MergeResult {
        unique_j,
        duplicate_records,
        duplicate_groups,
        duplicate_examples,
        sorted_index: DurableWrite {
            digest,
            sync_status,
        },
    })
}

fn audit_ledger(
    ledger_path: &Path,
    work_dir: &Path,
    sorted_index_output: &Path,
    receipt_output: &Path,
    chunk_records: usize,
    auditor_source_commit: String,
) -> Result<(), String> {
    if chunk_records == 0 {
        return Err("chunk records must be positive".into());
    }
    if !is_lower_hex(&auditor_source_commit, 40) {
        return Err("auditor source commit must be 40 lowercase hexadecimal digits".into());
    }
    if sorted_index_output.exists() || receipt_output.exists() {
        return Err("refusing to overwrite a sorted index or audit receipt".into());
    }
    fs::create_dir_all(work_dir)
        .map_err(|error| format!("create {}: {error}", work_dir.display()))?;
    let started = Instant::now();
    let ledger_digest = digest_file(ledger_path)?;
    let ledger: Ledger = read_json(ledger_path)?;
    if ledger.schema != LEDGER_SCHEMA
        || ledger.status != "sealed"
        || ledger.entries.is_empty()
        || !is_lower_hex(&ledger.builder_source_commit, 40)
        || ledger.x_start != ledger.entries[0].x_start
        || ledger.x_end != ledger.entries.last().map(|entry| entry.x_end).unwrap_or(0)
        || ledger.expected_complete_class_vertices != EXPECTED_CLASS_VERTICES
    {
        return Err("ledger header is malformed or does not target the frozen p61 class".into());
    }
    ensure_digest(
        &ledger.class_number_receipt,
        &ledger.class_number_receipt_digest,
        "class-number receipt",
    )?;
    let mut next_x = ledger.x_start;
    let mut prior_checkpoint_digest: Option<String> = None;
    let mut entries_buffer = Vec::with_capacity(chunk_records.min(5_000_000));
    let mut chunks = Vec::new();
    let mut segments = Vec::new();
    let mut total_records = 0u64;
    for entry in &ledger.entries {
        if entry.x_start != next_x || entry.x_end <= entry.x_start {
            return Err("ledger entries are not a positive contiguous partition".into());
        }
        validate_entry_artifacts(&ledger, entry, prior_checkpoint_digest.as_deref())?;
        segments.push(audit_segment(
            &ledger,
            entry,
            &mut entries_buffer,
            chunk_records,
            work_dir,
            &mut chunks,
        )?);
        total_records = total_records
            .checked_add(entry.node_records)
            .ok_or("audit record count overflow")?;
        prior_checkpoint_digest = Some(entry.checkpoint_digest.sha256.clone());
        next_x = entry.x_end;
    }
    if next_x != ledger.x_end || total_records != ledger.expected_node_records {
        return Err("ledger terminal range or node count differs".into());
    }
    flush_chunk(&mut entries_buffer, work_dir, &mut chunks)?;
    let merge = merge_chunks(&chunks, sorted_index_output, total_records)?;
    let complete_class_coverage = ledger.x_start == 0
        && ledger.x_end == ROW_ORDER
        && total_records == EXPECTED_CLASS_VERTICES;
    let status = if merge.duplicate_records > 0 {
        "duplicate_j_invariants"
    } else if complete_class_coverage && ledger.class_certification_state == "unconditional" {
        "complete_class_pass"
    } else if complete_class_coverage {
        "full_data_pass_certificate_blocked"
    } else {
        "partial_census_pass"
    };
    let receipt = AuditReceipt {
        schema: AUDIT_SCHEMA.into(),
        status: status.into(),
        auditor_source_commit,
        ledger_builder_source_commit: ledger.builder_source_commit,
        source_commit: ledger.source_commit,
        class_number_receipt_sha256: ledger.class_number_receipt_digest.sha256,
        class_certification_state: ledger.class_certification_state,
        ledger: ledger_path.to_path_buf(),
        ledger_digest,
        x_start: ledger.x_start,
        x_end: ledger.x_end,
        records: total_records,
        unique_j_invariants: merge.unique_j,
        duplicate_records: merge.duplicate_records,
        duplicate_groups: merge.duplicate_groups,
        duplicate_examples: merge.duplicate_examples,
        complete_class_coverage,
        chunks,
        sorted_index: sorted_index_output.to_path_buf(),
        sorted_index_digest: merge.sorted_index.digest,
        sorted_index_sync_status: merge.sorted_index.sync_status,
        segments,
        elapsed_ms: started.elapsed().as_millis(),
    };
    let written = write_json_new(receipt_output, &receipt)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "status": receipt.status,
            "receipt": receipt_output,
            "receipt_digest": written.digest,
            "receipt_sync_status": written.sync_status,
            "records": receipt.records,
            "unique_j_invariants": receipt.unique_j_invariants,
            "duplicate_records": receipt.duplicate_records,
            "sorted_index": receipt.sorted_index,
            "sorted_index_digest": receipt.sorted_index_digest,
        }))
        .map_err(|error| error.to_string())?
    );
    if receipt.duplicate_records > 0 {
        return Err(format!(
            "global uniqueness failed; evidence retained at {}",
            receipt_output.display()
        ));
    }
    Ok(())
}

fn run() -> Result<(), String> {
    match Cli::parse().command {
        Command::BuildLedger {
            generation_receipts_dir,
            verification_receipts_dir,
            class_number_receipt,
            class_certification_state,
            auditor_source_commit,
            output,
        } => build_ledger(
            &generation_receipts_dir,
            &verification_receipts_dir,
            &class_number_receipt,
            class_certification_state,
            auditor_source_commit,
            &output,
        ),
        Command::AuditLedger {
            ledger,
            work_dir,
            sorted_index_output,
            receipt_output,
            chunk_records,
            auditor_source_commit,
        } => audit_ledger(
            &ledger,
            &work_dir,
            &sorted_index_output,
            &receipt_output,
            chunk_records,
            auditor_source_commit,
        ),
    }
}

fn main() {
    if let Err(error) = run() {
        eprintln!("p61-isogeny-audit: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn j_parser_is_canonical_and_strict() {
        assert_eq!(parse_j("0000000000000000"), Ok(0));
        assert_eq!(parse_j("1ffffffffffffffe"), Ok(FIELD_P - 1));
        assert!(parse_j("1fffffffffffffff").is_err());
        assert!(parse_j("ABCDEF0000000000").is_err());
        assert!(parse_j("1").is_err());
    }

    #[test]
    fn local_bitmap_rejects_duplicates() {
        let mut bits = vec![0u64; 2];
        assert!(set_once(&mut bits, 65).is_ok());
        assert!(set_once(&mut bits, 65).is_err());
        assert_eq!(bits.iter().map(|word| word.count_ones()).sum::<u32>(), 1);
    }

    #[test]
    fn index_entry_binary_layout_round_trips() {
        let entry = IndexEntry {
            j: 0x0123_4567_89ab_cdef,
            index: 0xfedc_ba98_7654_3210,
        };
        let mut encoded = [0u8; 16];
        encoded[..8].copy_from_slice(&entry.j.to_le_bytes());
        encoded[8..].copy_from_slice(&entry.index.to_le_bytes());
        assert_eq!(
            IndexEntry {
                j: u64::from_le_bytes(encoded[..8].try_into().unwrap()),
                index: u64::from_le_bytes(encoded[8..].try_into().unwrap()),
            },
            entry
        );
    }

    #[test]
    fn external_merge_counts_duplicate_j_invariants() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let temporary = std::env::temp_dir().join(format!(
            "p61-isogeny-audit-test-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&temporary).unwrap();
        let mut entries = vec![
            IndexEntry { j: 7, index: 9 },
            IndexEntry { j: 3, index: 1 },
            IndexEntry { j: 7, index: 4 },
        ];
        let mut chunks = Vec::new();
        flush_chunk(&mut entries, &temporary, &mut chunks).unwrap();
        let output = temporary.join("index.bin");
        let merged = merge_chunks(&chunks, &output, 3).unwrap();
        assert_eq!(merged.unique_j, 2);
        assert_eq!(merged.duplicate_records, 1);
        assert_eq!(merged.duplicate_groups, 1);
        assert_eq!(merged.duplicate_examples.len(), 1);
        assert_eq!(merged.sorted_index.digest.bytes, 64);
        fs::remove_dir_all(&temporary).unwrap();
    }
}
