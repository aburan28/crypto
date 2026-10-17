//! Public synthetic-input research only. Statistical candidates are not breaks.
pub mod cipher;
pub mod diagnostics;
pub mod search;
pub mod stats;
pub mod trail;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use sha2::{Digest, Sha256};

pub fn digest(bytes: &[u8]) -> String {
    hex(&Sha256::digest(bytes))
}

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|x| format!("{x:02x}")).collect()
}

pub fn unhex(s: &str) -> Result<Vec<u8>, String> {
    if !s.is_ascii() || !s.len().is_multiple_of(2) {
        return Err("expected an even number of ASCII hexadecimal digits".into());
    }
    s.as_bytes()
        .chunks_exact(2)
        .map(|c| {
            let digit = |x: u8| match x {
                b'0'..=b'9' => Ok(x - b'0'),
                b'a'..=b'f' => Ok(x - b'a' + 10),
                b'A'..=b'F' => Ok(x - b'A' + 10),
                _ => Err("invalid hexadecimal digit".to_string()),
            };
            Ok((digit(c[0])? << 4) | digit(c[1])?)
        })
        .collect()
}

pub fn block(s: &str) -> Result<cipher::Block, String> {
    unhex(s)?
        .try_into()
        .map_err(|_| "expected a 16-byte block".into())
}

/// Deterministic, domain-separated public fixtures; not secret key generation.
pub fn stream(seed: u64, domain: &str, index: u64) -> ChaCha20Rng {
    let mut h = Sha256::new();
    h.update(b"aes-discover/public-fixtures/v1\0");
    h.update(seed.to_le_bytes());
    h.update((domain.len() as u64).to_le_bytes());
    h.update(domain.as_bytes());
    h.update(index.to_le_bytes());
    ChaCha20Rng::from_seed(h.finalize().into())
}
