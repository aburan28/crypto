use std::{
    env,
    path::{Path, PathBuf},
    process::Command,
};

fn is_lower_hex(value: &str, length: usize) -> bool {
    value.len() == length
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        && value.bytes().any(|byte| byte != b'0')
}

fn git_output(root: &Path, args: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(args)
        .output()
        .ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_owned())
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
    if let Some(repository) = git_output(&manifest_dir, &["rev-parse", "--show-toplevel"]) {
        let repository = PathBuf::from(repository);
        let head = git_output(&repository, &["rev-parse", "HEAD"])
            .expect("Git metadata exists but HEAD cannot be resolved");
        assert_eq!(
            commit, head,
            "P192_WCM_VERIFIER_COMMIT does not equal the containing repository HEAD"
        );
        let dirty = !Command::new("git")
            .arg("-C")
            .arg(&repository)
            .args(["diff-index", "--quiet", "HEAD", "--"])
            .status()
            .expect("Git metadata exists but cleanliness cannot be checked")
            .success()
            || !git_output(&repository, &["ls-files", "--others", "--exclude-standard"])
                .unwrap_or_default()
                .is_empty();
        let profile = env::var("PROFILE").expect("Cargo supplies PROFILE");
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
