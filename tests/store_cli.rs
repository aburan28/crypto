//! End to end through the `crypto store` command and the S3 backend, against
//! a stand-in `aws` executable that keeps a bucket in a directory and models
//! the behaviours the real transport depends on: absent keys answer 404,
//! `If-None-Match: *` refuses to overwrite, and a PUT whose bytes do not
//! match the SHA-256 it declared is rejected.
//!
//! The stand-in is a test fixture, not repository tooling: it exists so the
//! argument vectors in `src/store/backend.rs` are exercised through a real
//! subprocess boundary without a network or credentials. Nothing here is a
//! claim about real S3 behaviour beyond what those flags are documented to do.

#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const FAKE_AWS: &str = r#"#!/bin/sh
root="$FAKE_S3_ROOT"
cmd="$1"; op="$2"; shift 2
bucket=""; key=""; body=""; sum=""; meta=""; prefix=""; n=0; pos1=""; pos2=""
while [ $# -gt 0 ]; do
  case "$1" in
    --bucket) bucket="$2"; shift 2;;
    --key) key="$2"; shift 2;;
    --body) body="$2"; shift 2;;
    --checksum-sha256) sum="$2"; shift 2;;
    --metadata) meta="$2"; shift 2;;
    --prefix) prefix="$2"; shift 2;;
    --checksum-mode|--checksum-algorithm|--output|--endpoint-url|--if-none-match) shift 2;;
    --no-progress|--only-show-errors) shift;;
    *) n=$((n+1)); if [ "$n" = 1 ]; then pos1="$1"; else pos2="$1"; fi; shift;;
  esac
done
echo "$cmd $op" >> "$root/calls.log"
b64() { openssl dgst -sha256 -binary "$1" | base64; }
missing() { echo "An error occurred (404) when calling the HeadObject operation: Not Found" >&2; exit 254; }
case "$cmd $op" in
  "s3api head-object")
    f="$root/$bucket/$key"; [ -f "$f" ] || missing
    size=$(wc -c < "$f" | tr -d ' ')
    if [ -f "$f.sum" ]; then
      printf '{"ContentLength": %s, "ChecksumSHA256": "%s"}\n' "$size" "$(cat "$f.sum")"
    else
      printf '{"ContentLength": %s, "Metadata": {"sha256": "%s"}}\n' "$size" "$(cat "$f.meta" 2>/dev/null)"
    fi;;
  "s3api put-object")
    f="$root/$bucket/$key"
    if [ -f "$f" ]; then echo "An error occurred (PreconditionFailed) when calling the PutObject operation" >&2; exit 254; fi
    if [ "$(b64 "$body")" != "$sum" ]; then echo "An error occurred (BadDigest) when calling the PutObject operation" >&2; exit 254; fi
    mkdir -p "$(dirname "$f")"; cp "$body" "$f"; printf '%s' "$sum" > "$f.sum"
    printf '{"ChecksumSHA256": "%s"}\n' "$sum";;
  "s3 cp")
    case "$pos1" in
      s3://*) rest="${pos1#s3://}"; f="$root/$rest"; [ -f "$f" ] || { echo "An error occurred (404) when calling the HeadObject operation: Not Found" >&2; exit 1; }
              cp "$f" "$pos2";;
      *) rest="${pos2#s3://}"; f="$root/$rest"; mkdir -p "$(dirname "$f")"; cp "$pos1" "$f"; printf '%s' "${meta#sha256=}" > "$f.meta";;
    esac;;
  "s3api list-objects-v2")
    cd "$root/$bucket" 2>/dev/null || exit 0
    first=1; printf '{"Contents":['
    for f in $(find . -type f ! -name '*.sum' ! -name '*.meta' | sed 's|^\./||' | sort); do
      case "$f" in "$prefix"*) [ "$first" = 1 ] || printf ','; first=0
        printf '{"Key":"%s","Size":%s}' "$f" "$(wc -c < "$f" | tr -d ' ')";; esac
    done; printf ']}\n';;
  *) echo "fake aws: unsupported: $cmd $op" >&2; exit 2;;
esac
"#;

struct Env {
    root: PathBuf,
}

impl Env {
    fn new(tag: &str) -> Self {
        let root = std::env::temp_dir().join(format!("store-cli-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("s3")).unwrap();
        fs::write(root.join("aws"), FAKE_AWS).unwrap();
        fs::set_permissions(root.join("aws"), fs::Permissions::from_mode(0o755)).unwrap();
        Self { root }
    }

    /// Run `crypto store …` as one machine: its own cache, the shared index
    /// file and the shared fake bucket.
    fn store(&self, cache: &str, args: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_crypto"))
            .arg("store")
            .args(args)
            .current_dir(&self.root)
            .env("CRYPTO_STORE_INDEX", self.root.join("index.json"))
            .env("CRYPTO_STORE_CACHE", self.root.join(cache))
            .env("CRYPTO_STORE_REMOTE", "s3://test-bucket/store")
            .env("CRYPTO_STORE_AWS", self.root.join("aws"))
            .env("FAKE_S3_ROOT", self.root.join("s3"))
            .output()
            .unwrap()
    }

    fn calls(&self) -> String {
        fs::read_to_string(self.root.join("s3/calls.log")).unwrap_or_default()
    }

    fn object(&self, sha: &str) -> PathBuf {
        self.root
            .join(format!("s3/test-bucket/store/sha256/{}/{sha}", &sha[..2]))
    }
}

impl Drop for Env {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn text(o: &Output) -> String {
    format!(
        "{}{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    )
}

fn sha_of(path: &Path) -> String {
    let out = Command::new("sha256sum").arg(path).output().unwrap();
    String::from_utf8(out.stdout)
        .unwrap()
        .split_whitespace()
        .next()
        .unwrap()
        .to_string()
}

#[test]
fn push_on_one_machine_pull_on_another_through_s3() {
    let env = Env::new("roundtrip");
    let payload: Vec<u8> = (0..300_000u32).map(|i| (i % 253) as u8).collect();
    let file = env.root.join("fb.bin");
    fs::write(&file, &payload).unwrap();
    let sha = sha_of(&file);

    let push = env.store(
        "cache-a",
        &[
            "push",
            file.to_str().unwrap(),
            "--name",
            "factor-bases/demo/fb.bin",
            "--description",
            "demo factor base",
            "--tag",
            "fb",
        ],
    );
    assert!(push.status.success(), "{}", text(&push));
    assert!(text(&push).contains("uploaded"), "{}", text(&push));
    assert!(env.object(&sha).is_file());
    let index = fs::read_to_string(env.root.join("index.json")).unwrap();
    assert!(index.contains(&sha) && index.contains("factor-bases/demo/fb.bin"));
    assert!(env.calls().contains("s3api put-object"));

    // Pushing the same bytes again does not upload.
    let puts = env.calls().matches("put-object").count();
    let again = env.store(
        "cache-a",
        &[
            "push",
            file.to_str().unwrap(),
            "--name",
            "factor-bases/demo/fb.bin",
        ],
    );
    assert!(again.status.success(), "{}", text(&again));
    assert!(text(&again).contains("already in the remote"));
    assert_eq!(env.calls().matches("put-object").count(), puts);

    // Another machine, empty cache: pull downloads and verifies.
    let pull = env.store("cache-b", &["pull", "store:factor-bases/demo/fb.bin"]);
    assert!(pull.status.success(), "{}", text(&pull));
    let path = String::from_utf8(pull.stdout).unwrap().trim().to_string();
    assert_eq!(fs::read(&path).unwrap(), payload);
    assert!(path.ends_with(&sha));
    assert!(env.calls().contains("s3 cp"));
    let mode = fs::metadata(&path).unwrap().permissions().mode() & 0o777;
    assert_eq!(mode, 0o444, "cached objects are read-only");

    // A second pull is served from the cache without touching S3.
    let calls = env.calls().lines().count();
    let cached = env.store("cache-b", &["pull", &format!("sha256:{}", &sha[..16])]);
    assert!(cached.status.success(), "{}", text(&cached));
    assert_eq!(env.calls().lines().count(), calls);

    // --out copies, and refuses to overwrite.
    let out = env.root.join("copy.bin");
    let o = env.store(
        "cache-b",
        &[
            "pull",
            "factor-bases/demo/fb.bin",
            "--out",
            out.to_str().unwrap(),
        ],
    );
    assert!(o.status.success(), "{}", text(&o));
    assert_eq!(fs::read(&out).unwrap(), payload);
    let o = env.store(
        "cache-b",
        &[
            "pull",
            "factor-bases/demo/fb.bin",
            "--out",
            out.to_str().unwrap(),
        ],
    );
    assert!(!o.status.success() && text(&o).contains("not overwriting"));

    // verify, locally and against the remote.
    let v = env.store("cache-b", &["verify", "--remote"]);
    assert!(v.status.success(), "{}", text(&v));
    assert!(text(&v).contains("cache ok, remote ok"), "{}", text(&v));

    // ls, info and status read the same index.
    let ls = env.store("cache-b", &["ls", "factor-bases/"]);
    assert!(text(&ls).contains("factor-bases/demo/fb.bin"));
    let status = env.store("cache-b", &["status"]);
    assert!(
        text(&status).contains("ok (1 objects)"),
        "{}",
        text(&status)
    );
}

#[test]
fn a_tampered_remote_object_is_rejected_and_never_cached() {
    let env = Env::new("tamper");
    let file = env.root.join("t.bin");
    fs::write(&file, b"genuine bytes").unwrap();
    let sha = sha_of(&file);
    assert!(env
        .store("cache-a", &["push", file.to_str().unwrap(), "--name", "t"])
        .status
        .success());
    // Same length, different content.
    fs::write(env.object(&sha), b"forged bytes!").unwrap();
    let pull = env.store("cache-b", &["pull", "t"]);
    assert!(!pull.status.success());
    assert!(text(&pull).contains("expected"), "{}", text(&pull));
    assert!(!env
        .root
        .join("cache-b/sha256")
        .join(&sha[..2])
        .join(&sha)
        .exists());
    // And `verify --remote` reports the digest S3 recorded still matches the
    // index (S3 cannot see the forgery), which is why every download is hashed.
    let v = env.store("cache-b", &["verify", "--remote"]);
    assert!(v.status.success(), "{}", text(&v));
}

#[test]
fn a_name_is_not_repointed_without_replace_and_nothing_is_uploaded() {
    let env = Env::new("replace");
    let a = env.root.join("a.bin");
    let b = env.root.join("b.bin");
    fs::write(&a, b"first").unwrap();
    fs::write(&b, b"second").unwrap();
    assert!(env
        .store("c", &["push", a.to_str().unwrap(), "--name", "n"])
        .status
        .success());
    let puts = env.calls().matches("put-object").count();
    let clash = env.store("c", &["push", b.to_str().unwrap(), "--name", "n"]);
    assert!(!clash.status.success());
    assert!(text(&clash).contains("--replace"), "{}", text(&clash));
    assert_eq!(env.calls().matches("put-object").count(), puts);
    let moved = env.store(
        "c",
        &["push", b.to_str().unwrap(), "--name", "n", "--replace"],
    );
    assert!(moved.status.success(), "{}", text(&moved));
    let index = fs::read_to_string(env.root.join("index.json")).unwrap();
    assert!(index.contains("superseded"), "{index}");
    // rm edits the index only.
    let rm = env.store("c", &["rm", "n"]);
    assert!(rm.status.success(), "{}", text(&rm));
    assert!(env.object(&sha_of(&a)).is_file() && env.object(&sha_of(&b)).is_file());
}

#[test]
fn a_missing_aws_cli_is_reported_not_hidden() {
    let env = Env::new("noaws");
    let f = env.root.join("f");
    fs::write(&f, b"x").unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_crypto"))
        .args(["store", "push", f.to_str().unwrap(), "--name", "f"])
        .current_dir(&env.root)
        .env("CRYPTO_STORE_INDEX", env.root.join("index.json"))
        .env("CRYPTO_STORE_CACHE", env.root.join("c"))
        .env("CRYPTO_STORE_REMOTE", "s3://test-bucket/store")
        .env("CRYPTO_STORE_AWS", env.root.join("no-such-aws"))
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(text(&out).contains("CRYPTO_STORE_AWS"), "{}", text(&out));
    assert!(
        !env.root.join("index.json").exists(),
        "index untouched on failure"
    );
}

#[test]
fn a_store_reference_works_wherever_a_path_is_accepted() {
    // `ic --factor-base store:NAME` resolves through the same call. A name
    // that is not indexed fails with the store's message, before any other work.
    let env = Env::new("ic");
    let out = Command::new(env!("CARGO_BIN_EXE_ic"))
        .args(["run", "--factor-base", "store:absent/name"])
        .current_dir(&env.root)
        .env("CRYPTO_STORE_INDEX", env.root.join("index.json"))
        .env("CRYPTO_STORE_CACHE", env.root.join("c"))
        .output()
        .unwrap();
    assert!(
        text(&out).contains("is not in"),
        "unexpected output: {}",
        text(&out)
    );
}
