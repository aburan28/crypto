//! The index: the small, git-tracked map from a logical name to the
//! content hash, size and provenance of a large file held in the store.
//!
//! AGENTS.md asks that output too large for Git be represented by "a
//! manifest with its content hash, byte count, durable accessible location,
//! extraction command". This is that manifest, one for the whole repository.
//! The bytes live in the remote under their SHA-256, so an entry can never
//! point at different content than it was written for.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use serde::{Deserialize, Serialize};

pub const SCHEMA: &str = "crypto.store-index/v1";

/// What a name pointed at before `push --replace` moved it. Kept, not
/// deleted: a superseded figure stays visible as the "was".
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Superseded {
    pub sha256: String,
    pub bytes: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Entry {
    pub sha256: String,
    pub bytes: u64,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub description: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tags: Vec<String>,
    /// How the file was produced: the command, or what it was extracted from.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub provenance: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub superseded: Vec<Superseded>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Index {
    pub schema: String,
    pub entries: BTreeMap<String, Entry>,
}

impl Default for Index {
    fn default() -> Self {
        Self {
            schema: SCHEMA.to_string(),
            entries: BTreeMap::new(),
        }
    }
}

/// A logical name: path-like, printable, no way out of its namespace.
pub fn validate_name(name: &str) -> Result<(), String> {
    if name.is_empty() || name.len() > 200 {
        return Err("name must be 1 to 200 characters".into());
    }
    if name.starts_with('/') || name.ends_with('/') || name.contains("//") {
        return Err(format!(
            "name {name:?} must not start or end with '/' or contain '//'"
        ));
    }
    if name.split('/').any(|part| part == "." || part == "..") {
        return Err(format!("name {name:?} must not contain '.' or '..' parts"));
    }
    if !name
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'-' | b'/' | b'+' | b'='))
    {
        return Err(format!(
            "name {name:?} may only use letters, digits and . _ - / + ="
        ));
    }
    Ok(())
}

/// A SHA-256 digest as the store spells it: 64 lowercase hex digits.
pub fn validate_sha256(s: &str) -> Result<(), String> {
    if s.len() == 64 && s.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f')) {
        Ok(())
    } else {
        Err(format!("{s:?} is not a lowercase hex SHA-256"))
    }
}

impl Index {
    pub fn load(path: &Path) -> Result<Self, String> {
        let bytes = match fs::read(path) {
            Ok(b) => b,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(Self::default()),
            Err(e) => return Err(format!("{}: {e}", path.display())),
        };
        let index: Index = serde_json::from_slice(&bytes)
            .map_err(|e| format!("{}: invalid store index: {e}", path.display()))?;
        if index.schema != SCHEMA {
            return Err(format!(
                "{}: schema {:?}, expected {SCHEMA:?}",
                path.display(),
                index.schema
            ));
        }
        for (name, entry) in &index.entries {
            validate_name(name).map_err(|e| format!("{}: {e}", path.display()))?;
            validate_sha256(&entry.sha256)
                .map_err(|e| format!("{}: {name}: {e}", path.display()))?;
        }
        Ok(index)
    }

    /// Write sorted, pretty JSON with a trailing newline, atomically, so a
    /// diff of the index shows exactly the entries that changed.
    pub fn save(&self, path: &Path) -> Result<(), String> {
        if let Some(dir) = path.parent() {
            fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        }
        let mut text = serde_json::to_string_pretty(self).map_err(|e| e.to_string())?;
        text.push('\n');
        let tmp = path.with_extension(format!("json.tmp{}", std::process::id()));
        fs::write(&tmp, text).map_err(|e| format!("{}: {e}", tmp.display()))?;
        fs::rename(&tmp, path).map_err(|e| format!("{}: {e}", path.display()))
    }

    /// Find an entry by exact name.
    pub fn get(&self, name: &str) -> Option<&Entry> {
        self.entries.get(name)
    }

    /// Find the entries holding a digest (full, or a unique prefix of at
    /// least eight digits).
    pub fn by_sha256(&self, digest: &str) -> Vec<(&String, &Entry)> {
        if digest.len() < 8 {
            return Vec::new();
        }
        self.entries
            .iter()
            .filter(|(_, e)| e.sha256.starts_with(digest))
            .collect()
    }

    /// Insert or move a name. A name that already holds different content is
    /// refused unless `replace`, which keeps the old digest as superseded.
    pub fn put(&mut self, name: &str, mut entry: Entry, replace: bool) -> Result<(), String> {
        validate_name(name)?;
        validate_sha256(&entry.sha256)?;
        if let Some(old) = self.entries.get(name) {
            if old.sha256 == entry.sha256 {
                // Same content: keep what was recorded, fill in what was not.
                if entry.description.is_empty() {
                    entry.description = old.description.clone();
                }
                if entry.provenance.is_empty() {
                    entry.provenance = old.provenance.clone();
                }
                for tag in &old.tags {
                    if !entry.tags.contains(tag) {
                        entry.tags.push(tag.clone());
                    }
                }
                entry.superseded = old.superseded.clone();
            } else if replace {
                let mut was = old.superseded.clone();
                was.push(Superseded {
                    sha256: old.sha256.clone(),
                    bytes: old.bytes,
                });
                entry.superseded = was;
            } else {
                return Err(format!(
                    "{name} already holds sha256 {} ({} bytes); pass --replace to move it \
                     (the old digest is kept as superseded)",
                    old.sha256, old.bytes
                ));
            }
        }
        entry.tags.sort();
        entry.tags.dedup();
        self.entries.insert(name.to_string(), entry);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const A: &str = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    const B: &str = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";

    fn entry(sha: &str, bytes: u64) -> Entry {
        Entry {
            sha256: sha.into(),
            bytes,
            description: String::new(),
            tags: vec![],
            provenance: String::new(),
            superseded: vec![],
        }
    }

    #[test]
    fn names_cannot_escape_or_surprise() {
        for ok in ["fb/icv1-x/fb1-1.bin", "a", "a-b_c.d/e+f=g"] {
            assert!(validate_name(ok).is_ok(), "{ok}");
        }
        for bad in [
            "", "/a", "a/", "a//b", "../a", "a/../b", "a/./b", "a b", "a\\b", "é", "a\nb",
        ] {
            assert!(validate_name(bad).is_err(), "{bad:?}");
        }
        assert!(validate_name(&"x".repeat(201)).is_err());
    }

    #[test]
    fn digests_must_be_lowercase_hex_of_full_length() {
        assert!(validate_sha256(A).is_ok());
        assert!(validate_sha256(&A.to_uppercase()).is_err());
        assert!(validate_sha256(&A[..63]).is_err());
        assert!(validate_sha256(&format!("{}g", &A[..63])).is_err());
    }

    #[test]
    fn a_name_is_not_silently_repointed() {
        let mut ix = Index::default();
        ix.put("fb/one", entry(A, 10), false).unwrap();
        let err = ix.put("fb/one", entry(B, 20), false).unwrap_err();
        assert!(err.contains("--replace"), "{err}");
        assert_eq!(ix.get("fb/one").unwrap().sha256, A);
        ix.put("fb/one", entry(B, 20), true).unwrap();
        let now = ix.get("fb/one").unwrap();
        assert_eq!(now.sha256, B);
        assert_eq!(
            now.superseded,
            vec![Superseded {
                sha256: A.into(),
                bytes: 10
            }]
        );
    }

    #[test]
    fn republishing_the_same_content_keeps_the_record() {
        let mut ix = Index::default();
        let mut first = entry(A, 10);
        first.description = "kept".into();
        first.tags = vec!["fb".into()];
        ix.put("n", first, false).unwrap();
        let mut again = entry(A, 10);
        again.tags = vec!["extra".into()];
        ix.put("n", again, false).unwrap();
        let e = ix.get("n").unwrap();
        assert_eq!(e.description, "kept");
        assert_eq!(e.tags, vec!["extra".to_string(), "fb".to_string()]);
    }

    #[test]
    fn digest_prefix_lookup_needs_eight_digits() {
        let mut ix = Index::default();
        ix.put("n", entry(A, 1), false).unwrap();
        assert_eq!(ix.by_sha256("aaaaaaaa").len(), 1);
        assert!(ix.by_sha256("aaaaaaa").is_empty());
        assert!(ix.by_sha256("bbbbbbbb").is_empty());
    }

    #[test]
    fn save_then_load_round_trips_and_is_sorted() {
        let dir = std::env::temp_dir().join(format!("store-index-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        let path = dir.join("nested/index.json");
        let mut ix = Index::default();
        ix.put("z", entry(B, 2), false).unwrap();
        ix.put("a", entry(A, 1), false).unwrap();
        ix.save(&path).unwrap();
        let text = fs::read_to_string(&path).unwrap();
        assert!(text.ends_with("}\n"));
        assert!(text.find("\"a\"").unwrap() < text.find("\"z\"").unwrap());
        assert_eq!(Index::load(&path).unwrap(), ix);
        assert_eq!(
            Index::load(&dir.join("absent.json")).unwrap(),
            Index::default()
        );
        fs::write(&path, "{\"schema\":\"other\",\"entries\":{}}").unwrap();
        assert!(Index::load(&path).is_err());
        let _ = fs::remove_dir_all(&dir);
    }
}
