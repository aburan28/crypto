use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn full_lower_hex(value: &str) -> bool {
    value.len() == 40
        && value.bytes().any(|byte| byte != b'0')
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn run_git(root: &Path, arguments: &[&str]) -> std::io::Result<Output> {
    Command::new("/usr/bin/git")
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
}

fn command_output(root: &Path, arguments: &[&str]) -> String {
    let output = run_git(root, arguments)
        .unwrap_or_else(|error| panic!("cannot execute git {:?}: {error}", arguments));
    if !output.status.success() {
        panic!(
            "git {:?} failed: {}",
            arguments,
            String::from_utf8_lossy(&output.stderr).trim()
        );
    }
    String::from_utf8(output.stdout)
        .expect("git output must be UTF-8")
        .trim()
        .to_owned()
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

fn enforce_release_tree(repository: &Path, package: &Path, revision: &str) {
    let prefix = package
        .strip_prefix(repository)
        .expect("producer package must be inside its repository")
        .to_str()
        .expect("producer package path must be UTF-8");
    let index_rows = command_output(repository, &["ls-files", "-v", "--", prefix]);
    assert!(
        index_rows.lines().all(|row| row.starts_with("H ")),
        "producer package index contains assume-unchanged or skip-worktree entries"
    );
    let committed = command_output(
        repository,
        &[
            "ls-tree",
            "-r",
            "--format=%(objectmode) %(objectname) %(path)",
            revision,
            "--",
            prefix,
        ],
    );
    let mut committed_paths = BTreeSet::new();
    for description in committed.lines() {
        let pieces = description.splitn(3, ' ').collect::<Vec<_>>();
        assert_eq!(pieces.len(), 3, "malformed producer tree row");
        assert!(matches!(pieces[0], "100644" | "100755"));
        assert!(committed_paths.insert(pieces[2].to_owned()));
        let working_object = command_output(
            repository,
            &["hash-object", "--no-filters", "--", pieces[2]],
        );
        assert_eq!(
            working_object, pieces[1],
            "producer bytes differ from {revision}: {}",
            pieces[2]
        );
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mode = fs::metadata(repository.join(pieces[2]))
                .expect("inspect producer source mode")
                .permissions()
                .mode();
            assert_eq!(mode & 0o111 != 0, pieces[0] == "100755");
        }
    }
    assert!(
        !committed_paths.is_empty(),
        "producer package tree is empty"
    );

    fn walk(repository: &Path, package: &Path, directory: &Path, seen: &mut BTreeSet<String>) {
        let mut entries = fs::read_dir(directory)
            .expect("read producer source directory")
            .collect::<std::result::Result<Vec<_>, _>>()
            .expect("read producer source entry");
        entries.sort_by_key(|entry| entry.file_name());
        for entry in entries {
            let candidate = entry.path();
            let metadata = fs::symlink_metadata(&candidate).expect("inspect producer source path");
            assert!(!metadata.file_type().is_symlink());
            if directory == package && entry.file_name() == "target" {
                assert!(
                    metadata.is_dir(),
                    "producer target exclusion is not a directory"
                );
                continue;
            }
            if metadata.is_dir() {
                walk(repository, package, &candidate, seen);
            } else {
                assert!(metadata.is_file());
                seen.insert(
                    candidate
                        .strip_prefix(repository)
                        .expect("producer source prefix")
                        .to_str()
                        .expect("producer source path must be UTF-8")
                        .replace('\\', "/"),
                );
            }
        }
    }
    let mut working_paths = BTreeSet::new();
    walk(repository, package, package, &mut working_paths);
    assert_eq!(
        working_paths, committed_paths,
        "producer package inventory differs"
    );
}

fn main() {
    println!("cargo:rerun-if-env-changed=P192_WCM_PROTOCOL_COMMIT");
    println!("cargo:rerun-if-env-changed=P192_WCM_SOURCE_COMMIT");
    for input in [".gitignore", "Cargo.lock", "Cargo.toml", "build.rs", "src"] {
        println!("cargo:rerun-if-changed={input}");
    }

    let protocol_commit = env::var("P192_WCM_PROTOCOL_COMMIT")
        .expect("P192_WCM_PROTOCOL_COMMIT is required at build time");
    let source_commit = env::var("P192_WCM_SOURCE_COMMIT")
        .expect("P192_WCM_SOURCE_COMMIT is required at build time");
    if !full_lower_hex(&protocol_commit) {
        panic!("P192_WCM_PROTOCOL_COMMIT must be full lower-case 40-hex");
    }
    if !full_lower_hex(&source_commit) {
        panic!("P192_WCM_SOURCE_COMMIT must be full lower-case 40-hex");
    }

    let manifest_dir = PathBuf::from(
        env::var_os("CARGO_MANIFEST_DIR").expect("Cargo supplies CARGO_MANIFEST_DIR"),
    );
    let profile = env::var("PROFILE").expect("Cargo did not supply PROFILE to build.rs");

    // The environment values are the provenance inputs.  Git is only an
    // independent consistency gate when metadata accompanies the source.
    let git_metadata_present = run_git(&manifest_dir, &["rev-parse", "--is-inside-work-tree"])
        .map(|output| output.status.success())
        .unwrap_or(false);
    let dirty = if git_metadata_present {
        let head = command_output(&manifest_dir, &["rev-parse", "HEAD"]);
        if head != source_commit {
            panic!(
                "P192_WCM_SOURCE_COMMIT differs from Git HEAD: env={source_commit}, HEAD={head}"
            );
        }
        if profile == "release" {
            let repository = PathBuf::from(command_output(
                &manifest_dir,
                &["rev-parse", "--show-toplevel"],
            ));
            enforce_release_tree(&repository, &manifest_dir, &source_commit);
        }
        let status = command_output(
            &manifest_dir,
            &["status", "--porcelain=v1", "--untracked-files=normal"],
        );
        let dirty = !status.is_empty();
        let git_dir = PathBuf::from(command_output(
            &manifest_dir,
            &["rev-parse", "--absolute-git-dir"],
        ));
        emit_rerun_path(&git_dir.join("HEAD"), "Git HEAD path");
        emit_rerun_path(&git_dir.join("index"), "Git index path");
        dirty
    } else {
        false
    };
    if !git_metadata_present && profile == "release" {
        panic!("release producer builds require Git metadata for HEAD and clean-tree checks");
    }
    if dirty && profile == "release" {
        panic!("release producer builds require a clean Git worktree");
    }
    if !git_metadata_present {
        println!("cargo:warning=producer build has no Git metadata and is not admission-eligible");
    }
    if dirty {
        println!(
            "cargo:warning=dirty non-release producer build is not admission-eligible and preflight will refuse it"
        );
    }

    println!("cargo:rustc-env=P192_WCM_PROTOCOL_COMMIT={protocol_commit}");
    println!("cargo:rustc-env=P192_WCM_SOURCE_COMMIT={source_commit}");
    println!("cargo:rustc-env=P192_WCM_BUILD_PROFILE={profile}");
    println!(
        "cargo:rustc-env=P192_WCM_GIT_METADATA_PRESENT={}",
        if git_metadata_present {
            "true"
        } else {
            "false"
        }
    );
    println!(
        "cargo:rustc-env=P192_WCM_SOURCE_DIRTY={}",
        if dirty { "true" } else { "false" }
    );
}
