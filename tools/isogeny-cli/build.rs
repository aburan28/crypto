use std::{env, path::Path, process::Command};

fn git(repo: &Path, args: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(repo)
        .args(args)
        .output()
        .ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())
}

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    println!("cargo:rerun-if-env-changed=CRYPTO_BUILD_GIT_COMMIT");
    println!("cargo:rerun-if-env-changed=CRYPTO_BUILD_GIT_DIRTY");
    for reference in [
        "HEAD".to_string(),
        git(&root, &["symbolic-ref", "HEAD"]).unwrap_or_else(|| "HEAD".into()),
    ] {
        if let Some(path) = git(&root, &["rev-parse", "--git-path", &reference]) {
            println!("cargo:rerun-if-changed={path}");
        }
    }
    let revision = env::var("CRYPTO_BUILD_GIT_COMMIT")
        .ok()
        .or_else(|| git(&root, &["rev-parse", "HEAD"]))
        .unwrap_or_else(|| "unknown".into());
    let dirty = env::var("CRYPTO_BUILD_GIT_DIRTY").ok().unwrap_or_else(|| {
        git(
            &root,
            &["status", "--porcelain", "--untracked-files=normal"],
        )
        .map(|s| (!s.is_empty()).to_string())
        .unwrap_or_else(|| "unknown".into())
    });
    println!("cargo:rustc-env=ISOGENY_GIT_COMMIT={revision}");
    println!("cargo:rustc-env=ISOGENY_GIT_DIRTY={dirty}");
}
