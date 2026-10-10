//! Build the repository's identity checker without unrelated library modules.

#[path = "../../../src/hash/sha256.rs"]
pub mod sha256_source;
pub mod hash {
    pub use crate::sha256_source as sha256;
}

#[path = "../../../src/cryptanalysis/identity_certificate.rs"]
pub mod identity_certificate_source;
pub mod cryptanalysis {
    pub use crate::identity_certificate_source as identity_certificate;
}
