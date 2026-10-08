//! Build-time source provenance shared by the research binaries.
//!
//! Release and installation jobs set `CRYPTO_BUILD_GIT_COMMIT` explicitly.
//! GitHub Actions supplies `GITHUB_SHA` as a repository-neutral fallback.
//! Keeping these as compile-time values avoids consulting an unrelated
//! worktree after a binary has been copied or installed.

/// The source commit embedded in this build, when the builder supplied one.
pub fn git_commit() -> Option<&'static str> {
    option_env!("CRYPTO_BUILD_GIT_COMMIT").or(option_env!("GITHUB_SHA"))
}

/// Whether the embedded source tree was dirty, when the builder supplied it.
pub fn git_dirty() -> Option<bool> {
    option_env!("CRYPTO_BUILD_GIT_DIRTY").and_then(parse_bool)
}

/// Which compile-time input supplied [`git_commit`].
pub fn source() -> Option<&'static str> {
    if option_env!("CRYPTO_BUILD_GIT_COMMIT").is_some() {
        Some("CRYPTO_BUILD_GIT_COMMIT")
    } else if option_env!("GITHUB_SHA").is_some() {
        Some("GITHUB_SHA")
    } else {
        None
    }
}

fn parse_bool(value: &str) -> Option<bool> {
    match value {
        "1" | "true" | "TRUE" => Some(true),
        "0" | "false" | "FALSE" => Some(false),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::parse_bool;

    #[test]
    fn dirty_flag_is_strict() {
        assert_eq!(parse_bool("1"), Some(true));
        assert_eq!(parse_bool("true"), Some(true));
        assert_eq!(parse_bool("0"), Some(false));
        assert_eq!(parse_bool("false"), Some(false));
        assert_eq!(parse_bool("yes"), None);
    }
}
