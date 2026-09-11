//! **TLS 1.3 record layer** (RFC 8446 §5).
//!
//! Each record on the wire looks like:
//!
//! ```text
//!     struct {
//!         ContentType opaque_type = application_data;  /* always 23 for encrypted */
//!         ProtocolVersion legacy_record_version = 0x0303;
//!         uint16 length;       /* of (encrypted_record + auth_tag) */
//!         opaque encrypted_record[length];
//!     } TLSCiphertext;
//! ```
//!
//! Plaintext records (only ClientHello / ServerHello before keys are
//! installed) use:
//!
//! ```text
//!     struct {
//!         ContentType type;          /* handshake = 22, app_data = 23, alert = 21 */
//!         ProtocolVersion version = 0x0303;
//!         uint16 length;
//!         opaque fragment[length];
//!     } TLSPlaintext;
//! ```
//!
//! ## AEAD per-record nonce
//!
//! RFC 8446 §5.3: the per-record nonce is the **xor of the IV** (the
//! `*_traffic_iv` from the key schedule) with the 64-bit big-endian
//! sequence number, padded on the left with zeros to 12 bytes:
//!
//! ```text
//!     nonce_i = iv ⊕ (00 00 00 00 || BE64(seq_i))
//! ```
//!
//! ## AAD
//!
//! TLS 1.3 uses the cleartext record header as AAD for the AEAD:
//!
//! ```text
//!     aad = opaque_type || legacy_record_version || length
//! ```
//!
//! `length` is the length of the encrypted record (= payload +
//! content_type byte + auth tag).

use crate::symmetric::aes::{aes_gcm_decrypt, aes_gcm_encrypt, AesKey};

/// Content types (RFC 8446 §B.1).
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContentType {
    ChangeCipherSpec = 20,
    Alert = 21,
    Handshake = 22,
    ApplicationData = 23,
}

impl ContentType {
    pub fn from_u8(b: u8) -> Option<Self> {
        match b {
            20 => Some(ContentType::ChangeCipherSpec),
            21 => Some(ContentType::Alert),
            22 => Some(ContentType::Handshake),
            23 => Some(ContentType::ApplicationData),
            _ => None,
        }
    }
}

/// **Encode a plaintext TLS record** for the wire.
pub fn encode_tls_plaintext(ct: ContentType, payload: &[u8]) -> Vec<u8> {
    assert!(payload.len() <= 0xFFFF, "record fragment too large");
    let mut out = Vec::with_capacity(5 + payload.len());
    out.push(ct as u8);
    out.extend_from_slice(&[0x03, 0x03]); // legacy_record_version
    out.extend_from_slice(&(payload.len() as u16).to_be_bytes());
    out.extend_from_slice(payload);
    out
}

/// Parsed TLS record.
#[derive(Clone, Debug)]
pub struct ParsedRecord {
    pub content_type: ContentType,
    pub legacy_version: [u8; 2],
    pub payload: Vec<u8>,
    /// Total wire bytes consumed.
    pub consumed: usize,
}

/// **Decode** a TLS record from a byte slice.  Returns `None` if the
/// slice is too short or the content type is unknown.
pub fn parse_tls_record(buf: &[u8]) -> Option<ParsedRecord> {
    if buf.len() < 5 {
        return None;
    }
    let ct = ContentType::from_u8(buf[0])?;
    let legacy_version = [buf[1], buf[2]];
    let length = u16::from_be_bytes([buf[3], buf[4]]) as usize;
    if buf.len() < 5 + length {
        return None;
    }
    Some(ParsedRecord {
        content_type: ct,
        legacy_version,
        payload: buf[5..5 + length].to_vec(),
        consumed: 5 + length,
    })
}

/// **Compute the per-record AEAD nonce** per RFC 8446 §5.3.
pub fn derive_nonce(iv: &[u8; 12], seq: u64) -> [u8; 12] {
    let seq_be = seq.to_be_bytes();
    let mut nonce = *iv;
    // XOR seq into the last 8 bytes, big-endian.
    for i in 0..8 {
        nonce[4 + i] ^= seq_be[i];
    }
    nonce
}

/// **Encrypt a TLS 1.3 record** with AES-128-GCM.  Wraps the
/// `inner_plaintext = content_data || content_type || zero_padding`
/// pattern from RFC 8446 §5.2.
///
/// Returns the full wire bytes (including the 5-byte record header).
pub fn encrypt_record(
    key: &AesKey,
    iv: &[u8; 12],
    seq: u64,
    inner_content_type: ContentType,
    plaintext: &[u8],
) -> Vec<u8> {
    // Build TLSInnerPlaintext: content || content_type || 0..0 (no padding).
    let mut inner = plaintext.to_vec();
    inner.push(inner_content_type as u8);
    // Compute encrypted record length = inner.len() + tag(16).
    let ct_total_len = inner.len() + 16;
    // AAD = record header bytes (opaque_type=23, version=0x0303, length).
    let mut aad = [0u8; 5];
    aad[0] = ContentType::ApplicationData as u8;
    aad[1] = 0x03;
    aad[2] = 0x03;
    aad[3..5].copy_from_slice(&(ct_total_len as u16).to_be_bytes());
    let nonce = derive_nonce(iv, seq);
    let ct_with_tag = aes_gcm_encrypt(&inner, key, &nonce, &aad);
    debug_assert_eq!(ct_with_tag.len(), ct_total_len);
    let mut wire = aad.to_vec();
    wire.extend_from_slice(&ct_with_tag);
    wire
}

/// **Decrypt a TLS 1.3 record** given its full wire bytes.
/// Returns `(inner_content_type, plaintext)`.
pub fn decrypt_record(
    key: &AesKey,
    iv: &[u8; 12],
    seq: u64,
    wire: &[u8],
) -> Option<(ContentType, Vec<u8>)> {
    if wire.len() < 5 + 16 {
        return None;
    }
    if wire[0] != ContentType::ApplicationData as u8 {
        return None; // TLS 1.3 always uses opaque_type = 23 for encrypted
    }
    let length = u16::from_be_bytes([wire[3], wire[4]]) as usize;
    if wire.len() < 5 + length {
        return None;
    }
    let ct_with_tag = &wire[5..5 + length];
    let aad = &wire[0..5];
    let nonce = derive_nonce(iv, seq);
    let mut inner = aes_gcm_decrypt(ct_with_tag, key, &nonce, aad).ok()?;
    // Strip zero padding from the END, then read the trailing
    // content-type byte.
    while inner.last() == Some(&0) {
        inner.pop();
    }
    let content_type_byte = inner.pop()?;
    let content_type = ContentType::from_u8(content_type_byte)?;
    Some((content_type, inner))
}

// ── Tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Plaintext record round-trips.
    #[test]
    fn plaintext_record_round_trips() {
        let payload = b"hello tls".to_vec();
        let wire = encode_tls_plaintext(ContentType::Handshake, &payload);
        let parsed = parse_tls_record(&wire).expect("parse");
        assert_eq!(parsed.content_type, ContentType::Handshake);
        assert_eq!(parsed.payload, payload);
        assert_eq!(parsed.consumed, wire.len());
    }

    /// Encrypted record round-trips with matching key/iv/seq.
    #[test]
    fn encrypted_record_round_trips() {
        let key = AesKey::Aes128([0x42u8; 16]);
        let iv = [0xAAu8; 12];
        let plaintext = b"GET / HTTP/1.1\r\n\r\n".to_vec();
        let wire = encrypt_record(&key, &iv, 0, ContentType::ApplicationData, &plaintext);
        // Wire header: opaque_type=23, version=0x0303, length=...
        assert_eq!(wire[0], 23);
        assert_eq!(&wire[1..3], &[0x03, 0x03]);
        let (ct, pt) = decrypt_record(&key, &iv, 0, &wire).expect("decrypt");
        assert_eq!(ct, ContentType::ApplicationData);
        assert_eq!(pt, plaintext);
    }

    /// Wrong sequence number → decryption fails (the per-record
    /// nonce changed).
    #[test]
    fn wrong_seq_fails_decrypt() {
        let key = AesKey::Aes128([0u8; 16]);
        let iv = [0u8; 12];
        let wire = encrypt_record(&key, &iv, 5, ContentType::ApplicationData, b"data");
        assert!(decrypt_record(&key, &iv, 5, &wire).is_some());
        assert!(decrypt_record(&key, &iv, 6, &wire).is_none());
    }

    /// Tampered ciphertext → AEAD rejects.
    #[test]
    fn tampered_ct_fails() {
        let key = AesKey::Aes128([0u8; 16]);
        let iv = [0u8; 12];
        let mut wire = encrypt_record(&key, &iv, 0, ContentType::ApplicationData, b"secret");
        wire[8] ^= 1;
        assert!(decrypt_record(&key, &iv, 0, &wire).is_none());
    }

    /// **RFC 8446 §5.3 example**: `iv ⊕ seq_be_padded` gives the right
    /// per-record nonce.  We check at `seq = 0xDEADBEEF`.
    #[test]
    fn nonce_derivation_matches_rfc_8446_section_5_3() {
        let mut iv = [0u8; 12];
        for i in 0..12 {
            iv[i] = i as u8;
        }
        let n = derive_nonce(&iv, 0xDEAD_BEEFu64);
        // Expected: bytes 0..4 unchanged; bytes 4..12 = iv[4..12] xor BE64(0xDEADBEEF)
        // BE64(0xDEADBEEF) = 00 00 00 00 DE AD BE EF.
        let expect = [
            iv[0], iv[1], iv[2], iv[3],
            iv[4] ^ 0x00, iv[5] ^ 0x00, iv[6] ^ 0x00, iv[7] ^ 0x00,
            iv[8] ^ 0xDE, iv[9] ^ 0xAD, iv[10] ^ 0xBE, iv[11] ^ 0xEF,
        ];
        assert_eq!(n, expect);
    }
}
