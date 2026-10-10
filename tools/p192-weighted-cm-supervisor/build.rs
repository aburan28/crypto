use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn full_lower_hex(value: &str) -> bool {
    value.len() == 40
        && value.bytes().any(|byte| byte != b'0')
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn git_output(root: &Path, arguments: &[&str]) -> Option<String> {
    let output = Command::new("/usr/bin/git")
        .env_clear()
        .env("LC_ALL", "C")
        .env("GIT_CONFIG_NOSYSTEM", "1")
        .env("GIT_CONFIG_GLOBAL", "/dev/null")
        .env("GIT_NO_REPLACE_OBJECTS", "1")
        .env("GIT_GRAFT_FILE", "/dev/null")
        .env("GIT_NO_LAZY_FETCH", "1")
        .env("GIT_OPTIONAL_LOCKS", "0")
        .env("GIT_TERMINAL_PROMPT", "0")
        .env("GIT_PAGER", "cat")
        .env("GIT_LITERAL_PATHSPECS", "1")
        .arg("--no-pager")
        .arg("--no-replace-objects")
        .args(["-c", "core.fsmonitor=false"])
        .args(["-c", "core.commitGraph=false"])
        .args(["-c", "protocol.ext.allow=never"])
        .args(["-c", "core.hooksPath=/dev/null"])
        .arg("-C")
        .arg(root)
        .args(arguments)
        .output()
        .ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_owned())
}

fn emit_rerun_path(path: &Path, label: &str) {
    let text = path
        .to_str()
        .unwrap_or_else(|| panic!("{label} must be UTF-8 for Cargo provenance tracking"));
    assert!(
        !text.chars().any(char::is_control),
        "{label} contains a control character unsafe for Cargo build-script directives"
    );
    println!("cargo:rerun-if-changed={text}");
}

fn assert_release_package_snapshot(repository: &Path, package: &Path, commit: &str) {
    let relative = package
        .strip_prefix(repository)
        .expect("supervisor package is outside its repository")
        .to_str()
        .expect("supervisor package path must be UTF-8");
    let flags = git_output(repository, &["ls-files", "-v", "--", relative])
        .expect("inspect supervisor package index flags");
    assert!(
        flags.lines().all(|line| line.starts_with("H ")),
        "supervisor package has assume-unchanged, skip-worktree, or noncanonical index state"
    );
    let tree = git_output(
        repository,
        &[
            "ls-tree",
            "-r",
            "--format=%(objectmode) %(objectname) %(path)",
            commit,
            "--",
            relative,
        ],
    )
    .expect("read supervisor package tree");
    let mut expected = BTreeSet::new();
    for row in tree.lines() {
        let mut fields = row.splitn(3, ' ');
        let mode = fields.next().expect("tree mode");
        let object = fields.next().expect("tree object");
        let path = fields.next().expect("tree path");
        assert!(matches!(mode, "100644" | "100755"));
        assert!(expected.insert(path.to_owned()), "duplicate package path");
        let digest = git_output(repository, &["hash-object", "--no-filters", "--", path])
            .expect("hash supervisor working source");
        assert_eq!(
            digest, object,
            "supervisor source differs from {commit}: {path}"
        );
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let executable = fs::metadata(repository.join(path))
                .expect("inspect supervisor source mode")
                .permissions()
                .mode()
                & 0o111
                != 0;
            assert_eq!(executable, mode == "100755", "source mode differs: {path}");
        }
    }
    assert!(!expected.is_empty(), "supervisor package tree is empty");

    fn collect(repository: &Path, package: &Path, current: &Path, files: &mut BTreeSet<String>) {
        for entry in fs::read_dir(current).expect("read supervisor package") {
            let entry = entry.expect("read supervisor package entry");
            let path = entry.path();
            let metadata = fs::symlink_metadata(&path).expect("inspect supervisor package path");
            assert!(
                !metadata.file_type().is_symlink(),
                "package symlink: {}",
                path.display()
            );
            if current == package && entry.file_name() == "target" {
                assert!(
                    metadata.is_dir(),
                    "supervisor target exclusion is not a directory"
                );
                continue;
            }
            if metadata.is_dir() {
                collect(repository, package, &path, files);
            } else {
                assert!(
                    metadata.is_file(),
                    "non-file package object: {}",
                    path.display()
                );
                files.insert(
                    path.strip_prefix(repository)
                        .expect("supervisor package prefix")
                        .to_str()
                        .expect("supervisor source path must be UTF-8")
                        .replace('\\', "/"),
                );
            }
        }
    }
    let mut actual = BTreeSet::new();
    collect(repository, package, package, &mut actual);
    assert_eq!(
        actual, expected,
        "supervisor package inventory differs from {commit}"
    );
}

fn main() {
    for variable in [
        "P192_WCM_PROTOCOL_COMMIT",
        "P192_WCM_SOURCE_COMMIT",
        "P192_WCM_VERIFIER_COMMIT",
    ] {
        println!("cargo:rerun-if-env-changed={variable}");
    }
    for input in [".gitignore", "Cargo.lock", "Cargo.toml", "build.rs", "src"] {
        println!("cargo:rerun-if-changed={input}");
    }

    let protocol_commit = env::var("P192_WCM_PROTOCOL_COMMIT")
        .expect("P192_WCM_PROTOCOL_COMMIT is required for every supervisor build");
    let source_commit = env::var("P192_WCM_SOURCE_COMMIT")
        .expect("P192_WCM_SOURCE_COMMIT is required for every supervisor build");
    let verifier_commit = env::var("P192_WCM_VERIFIER_COMMIT")
        .expect("P192_WCM_VERIFIER_COMMIT is required for every supervisor build");
    for (label, value) in [
        ("P192_WCM_PROTOCOL_COMMIT", &protocol_commit),
        ("P192_WCM_SOURCE_COMMIT", &source_commit),
        ("P192_WCM_VERIFIER_COMMIT", &verifier_commit),
    ] {
        assert!(
            full_lower_hex(value),
            "{label} must be nonzero and exactly 40 lowercase hexadecimal characters"
        );
    }

    let manifest_dir = PathBuf::from(
        env::var_os("CARGO_MANIFEST_DIR").expect("Cargo supplies CARGO_MANIFEST_DIR"),
    );
    let profile = env::var("PROFILE").expect("Cargo supplies PROFILE");
    let (git_metadata_present, dirty) =
        if let Some(repository) = git_output(&manifest_dir, &["rev-parse", "--show-toplevel"]) {
            let repository = PathBuf::from(repository);
            let head = git_output(&repository, &["rev-parse", "HEAD"])
                .expect("Git metadata exists but HEAD cannot be resolved");
            assert_eq!(
                source_commit, head,
                "P192_WCM_SOURCE_COMMIT differs from the containing repository HEAD"
            );
            if profile == "release" {
                assert_release_package_snapshot(&repository, &manifest_dir, &source_commit);
            }
            let status = git_output(
                &repository,
                &["status", "--porcelain=v1", "--untracked-files=normal"],
            )
            .expect("Git metadata exists but cleanliness cannot be checked");
            let git_dir = PathBuf::from(
                git_output(&repository, &["rev-parse", "--absolute-git-dir"])
                    .expect("Git metadata exists but the Git directory cannot be resolved"),
            );
            emit_rerun_path(&git_dir.join("HEAD"), "Git HEAD path");
            emit_rerun_path(&git_dir.join("index"), "Git index path");
            (true, !status.is_empty())
        } else {
            (false, true)
        };
    if profile == "release" && !git_metadata_present {
        panic!("release supervisor builds require Git metadata");
    }
    if profile == "release" && dirty {
        panic!("release supervisor builds require a clean Git worktree");
    }
    if !git_metadata_present {
        println!("cargo:warning=Git metadata unavailable; supervisor is not admission eligible");
    } else if dirty {
        println!(
            "cargo:warning=dirty non-release supervisor build is test-only and not admission eligible"
        );
    }

    println!("cargo:rustc-env=P192_WCM_EMBEDDED_PROTOCOL_COMMIT={protocol_commit}");
    println!("cargo:rustc-env=P192_WCM_EMBEDDED_SUPERVISOR_COMMIT={source_commit}");
    println!("cargo:rustc-env=P192_WCM_EMBEDDED_VERIFIER_COMMIT={verifier_commit}");
    println!("cargo:rustc-env=P192_WCM_EMBEDDED_BUILD_PROFILE={profile}");
    println!(
        "cargo:rustc-env=P192_WCM_EMBEDDED_SUPERVISOR_DIRTY={}",
        if dirty { "true" } else { "false" }
    );
    println!(
        "cargo:rustc-env=P192_WCM_EMBEDDED_GIT_METADATA_PRESENT={}",
        if git_metadata_present {
            "true"
        } else {
            "false"
        }
    );
}
