use std::process::Command;

fn main() {
    println!("cargo:rerun-if-env-changed=IC_BUILD_GIT_COMMIT");
    let root = std::env::var("CARGO_MANIFEST_DIR").ok();
    if let Some(root) = &root {
        let head = std::path::Path::new(root).join(".git/HEAD");
        println!("cargo:rerun-if-changed={}", head.display());
        if let Ok(text) = std::fs::read_to_string(&head) {
            if let Some(reference) = text.trim().strip_prefix("ref: ") {
                // `.git/HEAD` itself does not change when a branch advances.
                // Track the resolved branch ref so a post-commit rebuild
                // cannot retain the previous source hash.
                println!(
                    "cargo:rerun-if-changed={}",
                    std::path::Path::new(root)
                        .join(".git")
                        .join(reference)
                        .display()
                );
            }
        }
    }
    let commit = std::env::var("IC_BUILD_GIT_COMMIT").ok().or_else(|| {
        let root = root.as_ref()?;
        let output = Command::new("git")
            .args(["-C", root.as_str(), "rev-parse", "HEAD"])
            .output()
            .ok()?;
        output
            .status
            .success()
            .then(|| String::from_utf8_lossy(&output.stdout).trim().to_owned())
    });
    if let Some(commit) = commit.filter(|s| !s.is_empty()) {
        println!("cargo:rustc-env=IC_BUILD_GIT_COMMIT={commit}");
    }
}
