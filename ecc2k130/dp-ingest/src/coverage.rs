//! How many records each slot has actually produced, from what it uploaded.

use crate::keys::{classify_key, slot_of_key, KeyKind};
use crate::RECORD_BYTES;
use std::collections::BTreeMap;

/// One stream's contribution from a bucket listing.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoverageRow {
    pub slot: u32,
    /// Orbit stream id, or `None` for the legacy key shape.
    pub stream: Option<String>,
    pub records: u64,
    /// Furthest `offset / 32 + records` any object in the stream reaches.
    pub covered: u64,
}

/// `(slot, stream, records, covered)` per stream from a bucket listing.
///
/// `objects` are `(key, records, _)` as `s3Objects` returns them.
pub fn coverage_rows(objects: &[(String, u64)]) -> Vec<CoverageRow> {
    let mut groups: BTreeMap<(u32, Option<String>), (u64, u64)> = BTreeMap::new();
    for (key, records) in objects {
        let Some(slot) = slot_of_key(key) else {
            continue;
        };
        let c = classify_key(key);
        let stream = match c.kind {
            KeyKind::Orbit => c.stream.map(str::to_string),
            _ => None,
        };
        let end = match (c.kind, c.offset) {
            (KeyKind::Orbit, Some(offset)) => offset / RECORD_BYTES as u64 + *records,
            _ => 0,
        };
        let acc = groups.entry((slot, stream)).or_insert((0, 0));
        acc.0 += *records;
        acc.1 = acc.1.max(end);
    }
    groups
        .into_iter()
        .map(|((slot, stream), (records, covered))| CoverageRow {
            slot,
            stream,
            records,
            covered,
        })
        .collect()
}

/// Records each slot has actually produced.
///
/// Within one stream offsets never overlap, so the furthest end is the count.
/// The legacy key shape carries no stream and is summed.
pub fn covered_records(rows: &[CoverageRow]) -> BTreeMap<u32, u64> {
    let mut out = BTreeMap::new();
    for row in rows {
        let n = if row.stream.is_some() {
            row.covered
        } else {
            row.records
        };
        *out.entry(row.slot).or_insert(0) += n;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coverage_from_a_listing_dedupes_re_uploaded_stretches() {
        let sid = "a8b4133d5c5f414ebd1337f15603588d";
        let key = |off: u64, sha: &str| {
            format!(
                "dp/slot-90001/{sid}-{off:016}-{sha}.bin",
                sha = sha.repeat(32)
            )
        };
        let legacy = "dp/slot-00002/1789311001-0000000000000000.bin".to_string();
        let objects = vec![
            (key(0, "aa"), 100),
            (key(3200, "bb"), 50),
            (key(0, "cc"), 150),
            (legacy, 7),
            (
                format!("dp/slot-00002/{sid}-{:016}-{:}.bin", 0, "ee".repeat(32)),
                5,
            ),
        ];
        let rows = coverage_rows(&objects);
        let covered = covered_records(&rows);
        assert_eq!(covered.get(&90001), Some(&150));
        assert_eq!(covered.get(&2), Some(&12));
    }
}
