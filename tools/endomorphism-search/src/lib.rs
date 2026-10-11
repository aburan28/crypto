//! Exact bounded searches. No secret-key inputs or discrete-log solver.
pub mod curves;
pub mod orders;
pub mod scalar;

// Reuse the study library's implementation for metadata and evidence hashes.
#[allow(dead_code)]
#[path = "../../../src/hash/sha256.rs"]
pub mod sha256;

pub const SCHEMA: &str = "endomorphism-search/native-v1";
pub type Result<T> = std::result::Result<T, String>;

pub fn digest(bytes: &[u8]) -> String {
    sha256::sha256(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
