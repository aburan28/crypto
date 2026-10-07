use std::{
    collections::BTreeSet,
    env, fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};

fn is_lower_hex(value: &str, length: usize) -> bool {
    value.len() == length
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        && value.bytes().any(|byte| byte != b'0')
}

fn git_raw(root: &Path, args: &[&str]) -> Option<Output> {
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
        .args(args)
        .output()
        .ok()
}

fn git_output(root: &Path, args: &[&str]) -> Option<String> {
    let output = git_raw(root, args)?;
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

fn verify_release_package(repository: &Path, package: &Path, commit: &str) {
    let subtree = package
        .strip_prefix(repository)
        .expect("verifier package must be inside its repository")
        .to_str()
        .expect("verifier package path must be UTF-8");
    let flag_rows = git_output(repository, &["ls-files", "-v", "--", subtree])
        .expect("inspect verifier index flags");
    assert!(
        flag_rows
            .lines()
            .all(|row| row.as_bytes().starts_with(b"H ")),
        "verifier package index has special flags"
    );
    let snapshot = git_output(
        repository,
        &[
            "ls-tree",
            "-r",
            "--format=%(objectmode) %(objectname) %(path)",
            commit,
            "--",
            subtree,
        ],
    )
    .expect("load verifier commit tree");
    let mut expected = BTreeSet::new();
    for row in snapshot.lines() {
        let (mode, tail) = row.split_once(' ').expect("verifier tree mode");
        let (object, relative) = tail.split_once(' ').expect("verifier tree object");
        assert!(mode == "100644" || mode == "100755");
        assert!(expected.insert(relative.to_owned()));
        let observed = git_output(repository, &["hash-object", "--no-filters", "--", relative])
            .expect("hash verifier working file");
        assert_eq!(observed, object, "verifier bytes changed: {relative}");
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let permissions = fs::metadata(repository.join(relative))
                .expect("verifier source metadata")
                .permissions()
                .mode();
            assert_eq!(permissions & 0o111 != 0, mode == "100755");
        }
    }
    assert!(!expected.is_empty());

    fn enumerate(
        repository: &Path,
        package: &Path,
        directory: &Path,
        output: &mut BTreeSet<String>,
    ) {
        let entries = fs::read_dir(directory).expect("read verifier package");
        for item in entries {
            let item = item.expect("read verifier package entry");
            let path = item.path();
            let kind = fs::symlink_metadata(&path).expect("verifier path metadata");
            assert!(!kind.file_type().is_symlink());
            if directory == package && item.file_name() == "target" {
                assert!(
                    kind.is_dir(),
                    "verifier target exclusion is not a directory"
                );
                continue;
            }
            if kind.is_dir() {
                enumerate(repository, package, &path, output);
            } else {
                assert!(kind.is_file());
                output.insert(
                    path.strip_prefix(repository)
                        .expect("verifier source prefix")
                        .to_str()
                        .expect("verifier source path must be UTF-8")
                        .replace('\\', "/"),
                );
            }
        }
    }
    let mut observed = BTreeSet::new();
    enumerate(repository, package, package, &mut observed);
    assert_eq!(observed, expected, "verifier package inventory changed");
}

fn main() {
    println!("cargo:rerun-if-env-changed=P192_WCM_VERIFIER_COMMIT");
    for input in [
        ".gitignore",
        "Cargo.lock",
        "Cargo.toml",
        "README.md",
        "build.rs",
        "src",
    ] {
        println!("cargo:rerun-if-changed={input}");
    }

    let commit = env::var("P192_WCM_VERIFIER_COMMIT").unwrap_or_else(|_| {
        panic!("P192_WCM_VERIFIER_COMMIT is required for every verifier build")
    });
    assert!(
        is_lower_hex(&commit, 40),
        "P192_WCM_VERIFIER_COMMIT must be nonzero and exactly 40 lowercase hexadecimal characters"
    );

    let manifest_dir = PathBuf::from(
        env::var_os("CARGO_MANIFEST_DIR").expect("Cargo supplies CARGO_MANIFEST_DIR"),
    );
    let profile = env::var("PROFILE").expect("Cargo supplies PROFILE");
    if let Some(repository) = git_output(&manifest_dir, &["rev-parse", "--show-toplevel"]) {
        let repository = PathBuf::from(repository);
        let head = git_output(&repository, &["rev-parse", "HEAD"])
            .expect("Git metadata exists but HEAD cannot be resolved");
        assert_eq!(
            commit, head,
            "P192_WCM_VERIFIER_COMMIT does not equal the containing repository HEAD"
        );
        if profile == "release" {
            verify_release_package(&repository, &manifest_dir, &commit);
        }
        let git_dir = PathBuf::from(
            git_output(&repository, &["rev-parse", "--absolute-git-dir"])
                .expect("Git metadata exists but the Git directory cannot be resolved"),
        );
        emit_rerun_path(&git_dir.join("HEAD"), "Git HEAD path");
        emit_rerun_path(&git_dir.join("index"), "Git index path");
        let untracked = git_raw(
            &repository,
            &["ls-files", "--others", "--exclude-standard", "-z"],
        )
        .expect("Git metadata exists but untracked files cannot be enumerated");
        if !untracked.status.success() {
            panic!(
                "Git untracked-file enumeration failed: {}",
                String::from_utf8_lossy(&untracked.stderr).trim()
            );
        }
        let dirty = !git_raw(&repository, &["diff-index", "--quiet", "HEAD", "--"])
            .expect("Git metadata exists but cleanliness cannot be checked")
            .status
            .success()
            || !untracked.stdout.is_empty();
        if dirty && profile == "release" {
            panic!("release verifier builds require a clean worktree");
        }
        if dirty {
            println!(
                "cargo:warning=dirty non-release verifier build is test-only and not admission eligible"
            );
        }
        println!(
            "cargo:rustc-env=P192_WCM_EMBEDDED_VERIFIER_DIRTY={}",
            if dirty { "true" } else { "false" }
        );
        println!("cargo:rustc-env=P192_WCM_EMBEDDED_BUILD_PROFILE={profile}");
    } else {
        // Admission builds are required to run from Git metadata.  Mark an
        // uninspectable development build dirty so runtime fails closed.
        println!("cargo:warning=Git metadata unavailable; verifier is not admission eligible");
        println!("cargo:rustc-env=P192_WCM_EMBEDDED_VERIFIER_DIRTY=true");
        println!("cargo:rustc-env=P192_WCM_EMBEDDED_BUILD_PROFILE=unknown");
    }
    println!("cargo:rustc-env=P192_WCM_EMBEDDED_VERIFIER_COMMIT={commit}");
}
