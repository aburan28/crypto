//! The object store every party shares: a local directory, S3 or GCS.
//!
//! Layout under the store root:
//!
//! ```text
//! bin/spotlab-<arch>                       agent binaries the bootstrap downloads
//! blobs/<sha256>                           uploaded inputs, content-addressed
//! jobs/<id>/spec.json                      the normalised spec (written once)
//! jobs/<id>/record.json                    job state (written by the controller only)
//! jobs/<id>/cancel                         cancel request marker
//! jobs/<id>/ckpt/<seq>/files/...           checkpoint contents
//! jobs/<id>/ckpt/<seq>/manifest.json       written last; a checkpoint exists iff this does
//! jobs/<id>/attempts/<n>/heartbeat.json    agent liveness and progress
//! jobs/<id>/attempts/<n>/exit.json         how the attempt ended (written once by its agent)
//! jobs/<id>/attempts/<n>/{stdout,stderr,agent}.log
//! jobs/<id>/output/...                     $SPOTLAB_OUTPUT_DIR of the finishing attempt
//! jobs/<id>/result.json                    the result (written once by the controller)
//! ```
//!
//! S3 and GCS are driven through the `aws` and `gcloud` command-line tools,
//! which on a cloud instance pick up its attached identity, so the agent
//! never handles a credential.

use anyhow::{bail, Context, Result};
use std::path::{Path, PathBuf};

use crate::util::{run, run_ok, s, write_atomic};

#[derive(Clone, Debug)]
pub enum Store {
    Local(PathBuf),
    S3(String),
    Gcs(String),
}

impl Store {
    pub fn open(url: &str) -> Result<Store> {
        let trimmed = url.trim_end_matches('/');
        if let Some(rest) = trimmed.strip_prefix("s3://") {
            if rest.is_empty() {
                bail!("store {url}: bucket missing");
            }
            Ok(Store::S3(trimmed.to_string()))
        } else if let Some(rest) = trimmed.strip_prefix("gs://") {
            if rest.is_empty() {
                bail!("store {url}: bucket missing");
            }
            Ok(Store::Gcs(trimmed.to_string()))
        } else {
            let p = PathBuf::from(trimmed);
            std::fs::create_dir_all(&p)
                .with_context(|| format!("create store dir {}", p.display()))?;
            Ok(Store::Local(p.canonicalize()?))
        }
    }

    pub fn url(&self) -> String {
        match self {
            Store::Local(p) => p.display().to_string(),
            Store::S3(u) | Store::Gcs(u) => u.clone(),
        }
    }

    fn full(&self, key: &str) -> String {
        match self {
            Store::Local(p) => p.join(key).display().to_string(),
            Store::S3(u) | Store::Gcs(u) => format!("{u}/{key}"),
        }
    }

    pub fn put(&self, key: &str, bytes: &[u8]) -> Result<()> {
        match self {
            Store::Local(p) => write_atomic(&p.join(key), bytes),
            Store::S3(_) => cli_ok(
                &s(&[
                    "aws",
                    "s3",
                    "cp",
                    "--only-show-errors",
                    "-",
                    &self.full(key),
                ]),
                Some(bytes),
            ),
            Store::Gcs(_) => cli_ok(
                &s(&["gcloud", "storage", "cp", "-", &self.full(key)]),
                Some(bytes),
            ),
        }
    }

    pub fn get(&self, key: &str) -> Result<Option<Vec<u8>>> {
        match self {
            Store::Local(p) => match std::fs::read(p.join(key)) {
                Ok(b) => Ok(Some(b)),
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
                Err(e) => Err(e.into()),
            },
            Store::S3(_) | Store::Gcs(_) => {
                let argv = match self {
                    Store::S3(_) => s(&[
                        "aws",
                        "s3",
                        "cp",
                        "--only-show-errors",
                        &self.full(key),
                        "-",
                    ]),
                    _ => s(&["gcloud", "storage", "cat", &self.full(key)]),
                };
                let out = run(&argv, None)?;
                if out.ok() {
                    Ok(Some(out.stdout))
                } else if is_missing(&out.stderr) {
                    Ok(None)
                } else {
                    bail!("{} failed: {}", argv.join(" "), out.stderr.trim())
                }
            }
        }
    }

    pub fn get_json<T: serde::de::DeserializeOwned>(&self, key: &str) -> Result<Option<T>> {
        match self.get(key)? {
            Some(b) => Ok(Some(
                serde_json::from_slice(&b).with_context(|| format!("parse {key}"))?,
            )),
            None => Ok(None),
        }
    }

    pub fn put_json<T: serde::Serialize>(&self, key: &str, v: &T) -> Result<()> {
        self.put(key, &serde_json::to_vec_pretty(v)?)
    }

    pub fn exists(&self, key: &str) -> Result<bool> {
        match self {
            Store::Local(p) => Ok(p.join(key).is_file()),
            _ => Ok(self.list(key)?.iter().any(|k| k == key)),
        }
    }

    pub fn put_file(&self, key: &str, path: &Path) -> Result<()> {
        match self {
            Store::Local(p) => {
                let dst = p.join(key);
                if let Some(d) = dst.parent() {
                    std::fs::create_dir_all(d)?;
                }
                let tmp = dst.with_extension(format!("tmp.{}", std::process::id()));
                std::fs::copy(path, &tmp)?;
                std::fs::rename(tmp, dst)?;
                Ok(())
            }
            Store::S3(_) => cli_ok(
                &s(&[
                    "aws",
                    "s3",
                    "cp",
                    "--only-show-errors",
                    &path.display().to_string(),
                    &self.full(key),
                ]),
                None,
            ),
            Store::Gcs(_) => cli_ok(
                &s(&[
                    "gcloud",
                    "storage",
                    "cp",
                    &path.display().to_string(),
                    &self.full(key),
                ]),
                None,
            ),
        }
    }

    /// Download `key` to `path`; false if it does not exist.
    pub fn get_file(&self, key: &str, path: &Path) -> Result<bool> {
        if let Some(d) = path.parent() {
            std::fs::create_dir_all(d)?;
        }
        match self {
            Store::Local(p) => {
                let src = p.join(key);
                if !src.is_file() {
                    return Ok(false);
                }
                std::fs::copy(src, path)?;
                Ok(true)
            }
            _ => {
                let argv = match self {
                    Store::S3(_) => s(&[
                        "aws",
                        "s3",
                        "cp",
                        "--only-show-errors",
                        &self.full(key),
                        &path.display().to_string(),
                    ]),
                    _ => s(&[
                        "gcloud",
                        "storage",
                        "cp",
                        &self.full(key),
                        &path.display().to_string(),
                    ]),
                };
                let out = run(&argv, None)?;
                if out.ok() {
                    Ok(true)
                } else if is_missing(&out.stderr) {
                    Ok(false)
                } else {
                    bail!("{} failed: {}", argv.join(" "), out.stderr.trim())
                }
            }
        }
    }

    /// Every key under `prefix`, sorted.
    pub fn list(&self, prefix: &str) -> Result<Vec<String>> {
        let mut keys: Vec<String> = match self {
            Store::Local(p) => {
                // List the deepest existing directory, then filter by the prefix.
                let base = match prefix.rfind('/') {
                    Some(i) => &prefix[..i],
                    None => "",
                };
                let dir = p.join(base);
                crate::util::walk_files(&dir)?
                    .into_iter()
                    .map(|rel| {
                        if base.is_empty() {
                            rel
                        } else {
                            format!("{base}/{rel}")
                        }
                    })
                    .filter(|k| k.starts_with(prefix) && !k.contains(".tmp."))
                    .collect()
            }
            Store::S3(u) => {
                let out = run(
                    &s(&["aws", "s3", "ls", "--recursive", &format!("{u}/{prefix}")]),
                    None,
                )?;
                if !out.ok() {
                    if out.stdout.is_empty() && (out.status == 1 || is_missing(&out.stderr)) {
                        return Ok(vec![]);
                    }
                    bail!("aws s3 ls failed: {}", out.stderr.trim());
                }
                // Lines are "date time size key"; keys are bucket-relative.
                let root = u.trim_start_matches("s3://");
                let root_prefix = match root.find('/') {
                    Some(i) => format!("{}/", &root[i + 1..]),
                    None => String::new(),
                };
                out.stdout_str()
                    .lines()
                    .filter_map(s3_ls_key)
                    .filter_map(|k| k.strip_prefix(&root_prefix).map(str::to_string))
                    .filter(|k| k.starts_with(prefix))
                    .collect()
            }
            Store::Gcs(u) => {
                let out = run(
                    &s(&["gcloud", "storage", "ls", &format!("{u}/{prefix}**")]),
                    None,
                )?;
                if !out.ok() {
                    if is_missing(&out.stderr) {
                        return Ok(vec![]);
                    }
                    bail!("gcloud storage ls failed: {}", out.stderr.trim());
                }
                let root = format!("{u}/");
                out.stdout_str()
                    .lines()
                    .filter_map(|l| l.trim().strip_prefix(&root).map(str::to_string))
                    .filter(|k| !k.ends_with('/'))
                    .collect()
            }
        };
        keys.sort();
        Ok(keys)
    }

    pub fn delete(&self, key: &str) -> Result<()> {
        match self {
            Store::Local(p) => match std::fs::remove_file(p.join(key)) {
                Ok(()) => Ok(()),
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
                Err(e) => Err(e.into()),
            },
            Store::S3(_) => cli_ok(
                &s(&["aws", "s3", "rm", "--only-show-errors", &self.full(key)]),
                None,
            ),
            Store::Gcs(_) => {
                let out = run(&s(&["gcloud", "storage", "rm", &self.full(key)]), None)?;
                if out.ok() || is_missing(&out.stderr) {
                    Ok(())
                } else {
                    bail!("gcloud storage rm failed: {}", out.stderr.trim())
                }
            }
        }
    }
}

fn cli_ok(argv: &[String], stdin: Option<&[u8]>) -> Result<()> {
    match stdin {
        None => run_ok(argv).map(|_| ()),
        Some(data) => {
            let out = run(argv, Some(data))?;
            if !out.ok() {
                bail!(
                    "{} exited {}: {}",
                    argv.join(" "),
                    out.status,
                    out.stderr.trim()
                );
            }
            Ok(())
        }
    }
}

/// The key from an `aws s3 ls --recursive` line: "date time size key".
fn s3_ls_key(line: &str) -> Option<String> {
    let mut rest = line;
    for _ in 0..3 {
        rest = rest.trim_start();
        let i = rest.find(char::is_whitespace)?;
        rest = &rest[i..];
    }
    let key = rest.trim_start();
    (!key.is_empty()).then(|| key.to_string())
}

fn is_missing(stderr: &str) -> bool {
    let e = stderr.to_ascii_lowercase();
    e.contains("404")
        || e.contains("not found")
        || e.contains("nosuchkey")
        || e.contains("no urls matched")
        || e.contains("does not exist")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_s3_listing_lines() {
        assert_eq!(
            s3_ls_key("2026-10-09 12:00:00       1234 p/jobs/a b.json").unwrap(),
            "p/jobs/a b.json"
        );
        assert!(s3_ls_key("garbage").is_none());
    }

    #[test]
    fn local_store_round_trip_and_list() {
        let dir = std::env::temp_dir().join(format!("spotlab-store-{}", std::process::id()));
        let st = Store::open(dir.to_str().unwrap()).unwrap();
        st.put("jobs/a/x.json", b"1").unwrap();
        st.put("jobs/a/ckpt/00000001/manifest.json", b"2").unwrap();
        st.put("jobs/b/x.json", b"3").unwrap();
        assert_eq!(st.get("jobs/a/x.json").unwrap().unwrap(), b"1");
        assert!(st.get("jobs/a/none").unwrap().is_none());
        assert_eq!(
            st.list("jobs/a/").unwrap(),
            vec!["jobs/a/ckpt/00000001/manifest.json", "jobs/a/x.json"]
        );
        assert_eq!(st.list("jobs/").unwrap().len(), 3);
        st.delete("jobs/a/x.json").unwrap();
        assert!(!st.exists("jobs/a/x.json").unwrap());
        std::fs::remove_dir_all(dir).ok();
    }
}
