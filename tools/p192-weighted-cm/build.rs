use std::env;
use std::path::PathBuf;
use std::process::Command;

fn full_lower_hex(value: &str) -> bool {
    value.len() == 40
        && value.bytes().any(|byte| byte != b'0')
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn command_output(arguments: &[&str]) -> String {
    let output = Command::new("git")
        .args(arguments)
        .output()
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

fn main() {
    println!("cargo:rerun-if-env-changed=P192_WCM_PROTOCOL_COMMIT");
    println!("cargo:rerun-if-env-changed=P192_WCM_SOURCE_COMMIT");

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

    // The environment values are the provenance inputs.  Git is only an
    // independent consistency gate when metadata accompanies the source.
    let git_metadata_present = Command::new("git")
        .args(["rev-parse", "--is-inside-work-tree"])
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false);
    let dirty = if git_metadata_present {
        let head = command_output(&["rev-parse", "HEAD"]);
        if head != source_commit {
            panic!(
                "P192_WCM_SOURCE_COMMIT differs from Git HEAD: env={source_commit}, HEAD={head}"
            );
        }
        let status = command_output(&["status", "--porcelain=v1", "--untracked-files=normal"]);
        let dirty = !status.is_empty();
        let git_dir = PathBuf::from(command_output(&["rev-parse", "--absolute-git-dir"]));
        println!("cargo:rerun-if-changed={}", git_dir.join("HEAD").display());
        println!("cargo:rerun-if-changed={}", git_dir.join("index").display());
        dirty
    } else {
        false
    };
    let profile = env::var("PROFILE").expect("Cargo did not supply PROFILE to build.rs");
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
