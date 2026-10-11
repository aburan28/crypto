//! The two `dp/` key shapes and how each names its upload time.

use regex::Regex;
use std::sync::OnceLock;

/// `dp/slot-00002/1789311001-0000000000000000.bin` — uploader clock in the name.
pub static LEGACY_KEY_RE: OnceLock<Regex> = OnceLock::new();
/// `dp/slot-00140/<32-hex>-<offset>-<64-hex>.bin` — content hash, no clock.
pub static ORBIT_KEY_RE: OnceLock<Regex> = OnceLock::new();
static SLOT_OF_KEY_RE: OnceLock<Regex> = OnceLock::new();

fn legacy() -> &'static Regex {
    LEGACY_KEY_RE.get_or_init(|| Regex::new(r"^dp/(slot-\d+)/(\d+)-(\d+)\.bin$").unwrap())
}

fn orbit() -> &'static Regex {
    ORBIT_KEY_RE.get_or_init(|| {
        Regex::new(r"^dp/(slot-\d+)/([0-9a-f]{32})-(\d+)-([0-9a-f]{64})\.bin$").unwrap()
    })
}

fn slot_re() -> &'static Regex {
    SLOT_OF_KEY_RE.get_or_init(|| Regex::new(r"^dp/slot-(\d+)/").unwrap())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KeyKind {
    Legacy,
    Orbit,
    Envelope,
    Unrecognised,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Classified<'a> {
    pub kind: KeyKind,
    pub key: &'a str,
    /// Orbit stream id (32 hex), when the key is orbit-shaped.
    pub stream: Option<&'a str>,
    /// Byte offset from the name, when present.
    pub offset: Option<u64>,
}

/// Classify one object key under `dp/`.
pub fn classify_key(key: &str) -> Classified<'_> {
    if key.ends_with(".bin.json") {
        return Classified {
            kind: KeyKind::Envelope,
            key,
            stream: None,
            offset: None,
        };
    }
    if let Some(c) = legacy().captures(key) {
        return Classified {
            kind: KeyKind::Legacy,
            key,
            stream: None,
            offset: c.get(3).and_then(|m| m.as_str().parse().ok()),
        };
    }
    if let Some(c) = orbit().captures(key) {
        return Classified {
            kind: KeyKind::Orbit,
            key,
            stream: c.get(2).map(|m| m.as_str()),
            offset: c.get(3).and_then(|m| m.as_str().parse().ok()),
        };
    }
    Classified {
        kind: KeyKind::Unrecognised,
        key,
        stream: None,
        offset: None,
    }
}

/// The identifier the existing rows use: the object key, slashes to dashes.
pub fn worker_id(key: &str) -> String {
    let rest = key.strip_prefix("dp/").unwrap_or(key);
    format!("dp-{}", rest.replace('/', "-"))
}

/// Upload time of one dp object, or `None` if the key is not a point object.
///
/// Legacy names carry the uploader's clock; orbit names need `last_modified`
/// (Unix seconds). Reading only the first shape is how five hours of a
/// 133-worker fleet went into `dp/` without a row reaching the store.
pub fn found_at(key: &str, last_modified: Option<i64>) -> Option<i64> {
    if let Some(c) = legacy().captures(key) {
        return c.get(2)?.as_str().parse().ok();
    }
    if orbit().is_match(key) {
        return last_modified;
    }
    None
}

/// Slot number from a `dp/slot-N/…` key.
pub fn slot_of_key(key: &str) -> Option<u32> {
    slot_re()
        .captures(key)
        .and_then(|c| c.get(1)?.as_str().parse().ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    const LEGACY: &str = "dp/slot-00002/1789311001-0000000000000000.bin";
    const ORBIT: &str =
        "dp/slot-00140/a8b4133d5c5f414ebd1337f15603588d-0000000000000000-bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb.bin";

    #[test]
    fn legacy_key_takes_found_at_from_the_name() {
        assert_eq!(found_at(LEGACY, None), Some(1_789_311_001));
        assert_eq!(classify_key(LEGACY).kind, KeyKind::Legacy);
    }

    #[test]
    fn orbit_key_takes_found_at_from_last_modified() {
        assert_eq!(found_at(ORBIT, Some(1_700_000_000)), Some(1_700_000_000));
        assert_eq!(found_at(ORBIT, None), None);
    }

    #[test]
    fn orbit_key_is_a_point_object() {
        let c = classify_key(ORBIT);
        assert_eq!(c.kind, KeyKind::Orbit);
        assert_eq!(c.stream, Some("a8b4133d5c5f414ebd1337f15603588d"));
        assert_eq!(c.offset, Some(0));
    }

    #[test]
    fn envelopes_and_foreign_keys_are_not_points() {
        assert_eq!(
            classify_key(&(LEGACY.to_string() + ".json")).kind,
            KeyKind::Envelope
        );
        assert_eq!(classify_key("dp/README.txt").kind, KeyKind::Unrecognised);
        assert_eq!(found_at("dp/README.txt", Some(1)), None);
    }

    #[test]
    fn worker_id_turns_slashes_to_dashes() {
        assert_eq!(worker_id("dp/slot-00002/1-0.bin"), "dp-slot-00002-1-0.bin");
    }

    #[test]
    fn slot_of_key_reads_both_shapes() {
        assert_eq!(slot_of_key(ORBIT), Some(140));
        assert_eq!(slot_of_key(LEGACY), Some(2));
        assert_eq!(slot_of_key("ckpt/slot-00002.ck"), None);
    }
}
