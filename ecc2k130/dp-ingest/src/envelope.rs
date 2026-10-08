//! Commit markers beside a corpus object, and the table-v3 header.

use crate::decode::RECORD_BYTES;
use sha2::{Digest, Sha256};

pub const DP_MAGIC_V2: &[u8; 8] = b"ECC2KDP2";
pub const DP_MAGIC_TABLE3: &[u8; 8] = b"ECC2KDT3";
pub const ENVELOPE_SUFFIX: &str = ".bin.json";

/// Strip the v3 table header. Signatures are checked over the original bytes.
pub fn table_record_body(body: &[u8]) -> Result<&[u8], String> {
    if body.len() < 8 || &body[..8] != DP_MAGIC_TABLE3.as_slice() {
        return Ok(body);
    }
    if body.len() < 16 {
        return Err("invalid table-v3 corpus header".into());
    }
    let version = u32::from_le_bytes(body[8..12].try_into().unwrap());
    let stride = u32::from_le_bytes(body[12..16].try_into().unwrap());
    if version != 3 || stride != RECORD_BYTES as u32 {
        return Err("invalid table-v3 corpus header".into());
    }
    Ok(&body[16..])
}

/// Refuse a body that does not match its own commit marker.
///
/// `body` is the raw GET bytes **before** [`table_record_body`]: the table-v3
/// header is inside the sha256 the worker wrote.
pub fn check_envelope(
    key: &str,
    body: &[u8],
    sha256_hex: Option<&str>,
    records: Option<u64>,
) -> Result<(), String> {
    if let Some(want) = sha256_hex {
        let got = hex::encode(Sha256::digest(body));
        if got != want {
            return Err(format!(
                "{key}: body sha256 {got} does not match envelope {want}"
            ));
        }
    }
    if let Some(want) = records {
        let records_body = table_record_body(body)?;
        if !(records_body.len() as u64).is_multiple_of(RECORD_BYTES as u64) {
            return Err(format!(
                "{key}: body length {} is not a multiple of {RECORD_BYTES}",
                records_body.len()
            ));
        }
        let got = records_body.len() as u64 / RECORD_BYTES as u64;
        if got != want {
            return Err(format!(
                "{key}: envelope records={want} but body holds {got}"
            ));
        }
    }
    Ok(())
}

/// Whether this body is a WITNESS=1 corpus the campaign must refuse.
pub fn is_witness_v2(body: &[u8]) -> bool {
    body.len() >= 8 && &body[..8] == DP_MAGIC_V2.as_slice()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_v3_header_is_stripped() {
        let mut body = DP_MAGIC_TABLE3.to_vec();
        body.extend_from_slice(&3u32.to_le_bytes());
        body.extend_from_slice(&(RECORD_BYTES as u32).to_le_bytes());
        body.extend_from_slice(&[7u8; 64]);
        assert_eq!(table_record_body(&body).unwrap(), &[7u8; 64]);
    }

    #[test]
    fn plain_body_passes_through() {
        let body = [1u8; 64];
        assert_eq!(table_record_body(&body).unwrap(), &body);
    }

    #[test]
    fn a_body_matching_its_manifest_is_accepted() {
        let body = [9u8; 64];
        let digest = hex::encode(Sha256::digest(body));
        check_envelope("dp/x.bin", &body, Some(&digest), Some(2)).unwrap();
    }

    #[test]
    fn a_corrupted_body_is_refused() {
        let body = [9u8; 64];
        let err = check_envelope("dp/x.bin", &body, Some(&"00".repeat(32)), Some(2)).unwrap_err();
        assert!(err.contains("sha256"), "{err}");
    }

    #[test]
    fn a_record_count_that_disagrees_is_refused() {
        let body = [9u8; 64];
        let digest = hex::encode(Sha256::digest(body));
        let err = check_envelope("dp/x.bin", &body, Some(&digest), Some(3)).unwrap_err();
        assert!(err.contains("records=3"), "{err}");
    }

    #[test]
    fn witness_v2_is_recognised() {
        let mut body = DP_MAGIC_V2.to_vec();
        body.extend_from_slice(&[0u8; 72]);
        assert!(is_witness_v2(&body));
        assert!(!is_witness_v2(&[0u8; 32]));
    }

    #[test]
    fn table_v3_hash_covers_the_header() {
        let mut body = DP_MAGIC_TABLE3.to_vec();
        body.extend_from_slice(&3u32.to_le_bytes());
        body.extend_from_slice(&(RECORD_BYTES as u32).to_le_bytes());
        body.extend_from_slice(&[1u8; 32]);
        let digest = hex::encode(Sha256::digest(&body));
        // Hashing the stripped body would disagree with the envelope.
        let stripped = table_record_body(&body).unwrap();
        assert_ne!(hex::encode(Sha256::digest(stripped)), digest);
        check_envelope("dp/t.bin", &body, Some(&digest), Some(1)).unwrap();
    }
}
