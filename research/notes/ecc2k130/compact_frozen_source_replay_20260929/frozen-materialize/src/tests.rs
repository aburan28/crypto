//! Fail-closed checks for the materializer. The first four are the contract of
//! `test_materialize.py`; the rest cover the base-tree restore.

use super::*;
use std::sync::atomic::{AtomicU32, Ordering};

static N: AtomicU32 = AtomicU32::new(0);

fn gz(raw: &[u8]) -> Vec<u8> {
    let mut e = GzEncoder::new(Vec::new(), Compression::default());
    std::io::Write::write_all(&mut e, raw).unwrap();
    e.finish().unwrap()
}

fn tar_of(files: &[(&str, &[u8])]) -> Vec<u8> {
    let mut b = tar::Builder::new(Vec::new());
    for (path, data) in files {
        let mut h = tar::Header::new_gnu();
        h.set_size(data.len() as u64);
        h.set_mode(0o644);
        h.set_cksum();
        b.append_data(&mut h, path, *data).unwrap();
    }
    b.into_inner().unwrap()
}

struct Fx {
    tmp: PathBuf,
    root: PathBuf,
    replay: PathBuf,
    freeze: PathBuf,
}

const FROZEN: &[u8] = b"frozen bytes\n";
const FREEZE_REL: &str = "FROZEN.json";

impl Fx {
    /// A checkout with one drifted pinned file (`src/a.rs`), one unchanged
    /// pinned file (`res/b.txt`), a live-only file in `src/`, and a replay
    /// directory holding the snapshot of the drifted file and a base archive
    /// whose `src/` has an unpinned file the live tree lacks.
    fn new() -> Self {
        let tmp = std::env::temp_dir().join(format!(
            "frozen-materialize-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::SeqCst)
        ));
        let _ = fs::remove_dir_all(&tmp);
        let root = tmp.join("root");
        let replay = tmp.join("replay");
        fs::create_dir_all(root.join("src")).unwrap();
        fs::create_dir_all(root.join("res")).unwrap();
        fs::create_dir_all(replay.join("historical_snapshots")).unwrap();
        fs::create_dir_all(replay.join("bases")).unwrap();
        fs::write(root.join("src/a.rs"), b"drifted bytes\n").unwrap();
        fs::write(root.join("src/live_only.rs"), b"newer than the freeze\n").unwrap();
        fs::write(root.join("res/b.txt"), b"unchanged\n").unwrap();
        let freeze = root.join(FREEZE_REL);
        fs::write(
            &freeze,
            json!({"source_sha256": {"src/a.rs": sha_hex(FROZEN), "res/b.txt": sha_hex(b"unchanged\n")}})
                .to_string(),
        )
        .unwrap();
        let fx = Self {
            tmp,
            root,
            replay,
            freeze,
        };
        fx.write_snapshot(FROZEN, "src/a.rs", None);
        fx.write_base(&[("src/a.rs", FROZEN), ("src/base_only.rs", b"base file\n")]);
        fx
    }

    fn write_snapshot(&self, raw: &[u8], path: &str, digest: Option<String>) {
        let digest = digest.unwrap_or_else(|| sha_hex(raw));
        let packed = gz(raw);
        fs::write(
            self.replay
                .join(format!("historical_snapshots/{digest}.gz")),
            &packed,
        )
        .unwrap();
        fs::write(
            self.replay.join("SNAPSHOTS.json"),
            json!({"schema": SNAPSHOTS_SCHEMA, "snapshots": [{
                "path": path, "sha256": digest, "bytes": raw.len(),
                "gzip_sha256": sha_hex(&packed)}]})
            .to_string(),
        )
        .unwrap();
    }

    fn write_base(&self, files: &[(&str, &[u8])]) {
        let packed = gz(&tar_of(files));
        fs::write(self.replay.join("bases/base.tar.gz"), &packed).unwrap();
        fs::write(
            self.replay.join("BASES.json"),
            json!({"schema": BASES_SCHEMA, "bases": [{
                "commit": "0123456789abcdef", "basis": "source_base_main",
                "archive": "bases/base.tar.gz", "archive_sha256": sha_hex(&packed),
                "archive_bytes": packed.len(), "freezes": [FREEZE_REL]}]})
            .to_string(),
        )
        .unwrap();
    }

    fn run(&self) -> Res<Value> {
        materialize(
            &self.root,
            &[self.freeze.clone()],
            &self.tmp.join("out"),
            &self.replay,
        )
    }
}

impl Drop for Fx {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.tmp);
    }
}

#[test]
fn a_parent_directory_member_of_a_base_path_is_accepted() {
    let fx = Fx::new();
    let mut b = tar::Builder::new(Vec::new());
    let mut dir = tar::Header::new_gnu();
    dir.set_entry_type(tar::EntryType::Directory);
    dir.set_size(0);
    dir.set_mode(0o755);
    dir.set_cksum();
    b.append_data(&mut dir, "docs/", &[][..]).unwrap();
    let mut dir2 = dir.clone();
    b.append_data(&mut dir2, "docs/ecbench/", &[][..]).unwrap();
    let data = b"schema";
    let mut f = tar::Header::new_gnu();
    f.set_size(data.len() as u64);
    f.set_mode(0o644);
    f.set_cksum();
    b.append_data(&mut f, "docs/ecbench/schema.sql", &data[..])
        .unwrap();
    let mut f2 = tar::Header::new_gnu();
    f2.set_size(FROZEN.len() as u64);
    f2.set_mode(0o644);
    f2.set_cksum();
    b.append_data(&mut f2, "src/a.rs", FROZEN).unwrap();
    let packed = gz(&b.into_inner().unwrap());
    fs::write(fx.replay.join("bases/base.tar.gz"), &packed).unwrap();
    fs::write(
        fx.replay.join("BASES.json"),
        json!({"schema": BASES_SCHEMA, "bases": [{
            "commit": "0123456789abcdef", "basis": "source_base_main",
            "archive": "bases/base.tar.gz", "archive_sha256": sha_hex(&packed),
            "archive_bytes": packed.len(), "freezes": [FREEZE_REL]}]})
        .to_string(),
    )
    .unwrap();
    fx.run().unwrap();
    assert_eq!(
        fs::read(fx.tmp.join("out/docs/ecbench/schema.sql")).unwrap(),
        b"schema"
    );
}

#[test]
fn a_pinned_file_is_restored_to_its_frozen_bytes() {
    let fx = Fx::new();
    let r = fx.run().unwrap();
    let out = fx.tmp.join("out");
    assert_eq!(fs::read(out.join("src/a.rs")).unwrap(), FROZEN);
    assert_eq!(
        r["replayed_from_snapshot"].as_object().unwrap().len(),
        0,
        "the base already holds the frozen bytes: {r}"
    );
    assert_eq!(r["schema"], "compact-frozen-source-materialization-v1");
    assert_eq!(r["pinned_files"], 2);
    assert_eq!(r["live"], json!(["res/b.txt", "src/a.rs"]));
    assert_eq!(r["base"]["commit"], "0123456789abcdef");
    // The checkout itself is never modified.
    assert_eq!(
        fs::read(fx.root.join("src/a.rs")).unwrap(),
        b"drifted bytes\n"
    );
}

#[test]
fn a_pinned_file_the_base_does_not_hold_comes_from_the_snapshot() {
    let fx = Fx::new();
    fx.write_base(&[("src/base_only.rs", b"base file\n")]);
    let r = fx.run().unwrap();
    assert_eq!(fs::read(fx.tmp.join("out/src/a.rs")).unwrap(), FROZEN);
    assert_eq!(
        r["replayed_from_snapshot"]
            .as_object()
            .unwrap()
            .keys()
            .collect::<Vec<_>>(),
        ["src/a.rs"]
    );
}

#[test]
fn missing_snapshot_fails_closed() {
    let fx = Fx::new();
    fx.write_base(&[("src/base_only.rs", b"base file\n")]);
    fs::write(fx.root.join("res/b.txt"), b"also drifted\n").unwrap();
    let e = fx.run().unwrap_err();
    assert!(e.contains("no snapshot"), "{e}");
}

#[test]
fn tampered_snapshot_fails_closed() {
    let fx = Fx::new();
    fx.write_base(&[("src/base_only.rs", b"base file\n")]);
    fx.write_snapshot(b"not the frozen bytes\n", "src/a.rs", Some(sha_hex(FROZEN)));
    assert!(fx.run().is_err());
}

#[test]
fn snapshot_for_another_path_is_rejected() {
    let fx = Fx::new();
    fx.write_base(&[("src/base_only.rs", b"base file\n")]);
    fx.write_snapshot(FROZEN, "src/other.rs", None);
    let e = fx.run().unwrap_err();
    assert!(e.contains("src/other.rs"), "{e}");
}

#[test]
fn the_base_tree_replaces_src_and_drops_files_the_base_never_had() {
    let fx = Fx::new();
    fx.run().unwrap();
    let out = fx.tmp.join("out");
    assert_eq!(
        fs::read(out.join("src/base_only.rs")).unwrap(),
        b"base file\n"
    );
    assert!(
        !out.join("src/live_only.rs").exists(),
        "a live-only file must not leak in"
    );
    // The checkout keeps its own file.
    assert!(fx.root.join("src/live_only.rs").exists());
}

#[test]
fn a_tampered_base_archive_fails_closed() {
    let fx = Fx::new();
    let mut packed = fs::read(fx.replay.join("bases/base.tar.gz")).unwrap();
    let last = packed.len() / 2;
    packed[last] ^= 0x55;
    fs::write(fx.replay.join("bases/base.tar.gz"), packed).unwrap();
    let e = fx.run().unwrap_err();
    assert!(e.contains("archive"), "{e}");
}

#[test]
fn a_freeze_without_a_base_fails_closed() {
    let fx = Fx::new();
    let doc = json!({"schema": BASES_SCHEMA, "bases": []});
    fs::write(fx.replay.join("BASES.json"), doc.to_string()).unwrap();
    let e = fx.run().unwrap_err();
    assert!(e.contains("no base tree"), "{e}");
}

#[test]
fn an_archive_member_outside_the_base_paths_is_refused() {
    let fx = Fx::new();
    fx.write_base(&[
        ("src/a.rs", FROZEN),
        ("docs/other.md", b"not a base path\n"),
    ]);
    let e = fx.run().unwrap_err();
    assert!(e.contains("outside the base paths"), "{e}");
}

#[test]
fn an_existing_destination_is_never_overwritten() {
    let fx = Fx::new();
    fs::create_dir_all(fx.tmp.join("out")).unwrap();
    assert!(fx.run().unwrap_err().contains("never overwrite"));
}

#[test]
fn conflicting_pins_are_refused() {
    let fx = Fx::new();
    let other = fx.root.join("OTHER.json");
    fs::write(
        &other,
        json!({"source_sha256": {"src/a.rs": "00"}}).to_string(),
    )
    .unwrap();
    let e = materialize(
        &fx.root,
        &[fx.freeze.clone(), other],
        &fx.tmp.join("out"),
        &fx.replay,
    )
    .unwrap_err();
    assert!(e.contains("conflicting pins"), "{e}");
}
