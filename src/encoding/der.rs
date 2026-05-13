//! Minimal DER encoder/decoder for ECDSA signatures.
//!
//! ECDSA signatures over the wire are DER-encoded `SEQUENCE { r INTEGER, s INTEGER }`.
//! The full DER specification (X.690) is large; we ship just what's
//! needed to round-trip ECDSA signatures from Bitcoin / TLS / X.509.
//!
//! # Wire format
//!
//! ```text
//!   0x30 [total-length]
//!     0x02 [r-length] [r-bytes]
//!     0x02 [s-length] [s-bytes]
//! ```
//!
//! Each INTEGER is **big-endian, two's-complement**.  Critically
//! for our use:
//!
//! - If the high bit of the magnitude is set, a `0x00` prefix
//!   byte must be added (so the value is unambiguously
//!   positive).  Bitcoin's "low-S" rule depends on this.
//! - Leading `0x00` bytes are forbidden unless required by the
//!   above rule (Bitcoin BIP-66 strict-DER).
//!
//! We implement the strict form: encoding produces canonical
//! output, decoding rejects non-canonical inputs with
//! [`DerError::NonCanonical`].
//!
//! # What this is NOT
//!
//! - A general DER parser.  We don't handle indefinite-length,
//!   tagged context-specific types, OBJECT IDENTIFIERs, or any
//!   of the X.509 universe.  See `der` / `pkcs1` / `x509-parser`
//!   crates for that.
//! - A BER parser.  BER permits non-minimal forms; DER does not.

use num_bigint::BigUint;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DerError {
    /// Truncated input.
    Truncated,
    /// Wrong tag where SEQUENCE / INTEGER expected.
    UnexpectedTag(u8),
    /// Length encoding malformed (e.g., long-form length we don't support).
    BadLength,
    /// Trailing bytes after the SEQUENCE end.
    TrailingBytes,
    /// Encoded value is not in canonical DER form (leading zeros,
    /// missing zero-prefix when high bit is set, total length
    /// mismatch, …).
    NonCanonical(&'static str),
    /// `r` or `s` is zero — invalid ECDSA signature.
    ZeroComponent,
}

/// Encode `(r, s)` as a strict-DER ECDSA signature.
///
/// Returns the variable-length DER encoding (typically 70-72
/// bytes for 256-bit `n`).  Both components must be non-zero.
pub fn encode_ecdsa_signature(r: &BigUint, s: &BigUint) -> Result<Vec<u8>, DerError> {
    if r.bits() == 0 || s.bits() == 0 {
        return Err(DerError::ZeroComponent);
    }
    let r_bytes = encode_integer(r);
    let s_bytes = encode_integer(s);
    let inner_len = 2 + r_bytes.len() + 2 + s_bytes.len();
    if inner_len > 0x7F {
        // Real ECDSA signatures fit comfortably in short-form
        // length encoding.  Anything > 127 bytes is almost
        // certainly a bug.
        return Err(DerError::BadLength);
    }
    let mut out = Vec::with_capacity(2 + inner_len);
    out.push(0x30); // SEQUENCE
    out.push(inner_len as u8);
    out.push(0x02); // INTEGER
    out.push(r_bytes.len() as u8);
    out.extend_from_slice(&r_bytes);
    out.push(0x02); // INTEGER
    out.push(s_bytes.len() as u8);
    out.extend_from_slice(&s_bytes);
    Ok(out)
}

/// Encode a positive `BigUint` as the body of a DER INTEGER.
/// Adds a `0x00` prefix iff the high bit of the magnitude is
/// set (so the integer reads as positive in two's complement).
fn encode_integer(n: &BigUint) -> Vec<u8> {
    let mut bytes = n.to_bytes_be();
    if bytes.is_empty() {
        bytes.push(0x00);
    }
    if bytes[0] & 0x80 != 0 {
        let mut out = Vec::with_capacity(bytes.len() + 1);
        out.push(0x00);
        out.extend_from_slice(&bytes);
        out
    } else {
        bytes
    }
}

/// Decode a strict-DER ECDSA signature into `(r, s)`.  Rejects
/// non-canonical encodings (per BIP-66).
pub fn decode_ecdsa_signature(input: &[u8]) -> Result<(BigUint, BigUint), DerError> {
    if input.len() < 2 {
        return Err(DerError::Truncated);
    }
    if input[0] != 0x30 {
        return Err(DerError::UnexpectedTag(input[0]));
    }
    // Short-form length only (BIP-66 forbids long-form for
    // signatures, since they can't be > 127 bytes).
    let total_len = input[1];
    if total_len & 0x80 != 0 {
        return Err(DerError::BadLength);
    }
    let body = &input[2..];
    if body.len() != total_len as usize {
        return Err(DerError::NonCanonical("SEQUENCE length mismatch"));
    }

    let (r, rest) = decode_integer(body)?;
    let (s, rest) = decode_integer(rest)?;
    if !rest.is_empty() {
        return Err(DerError::TrailingBytes);
    }
    if r.bits() == 0 || s.bits() == 0 {
        return Err(DerError::ZeroComponent);
    }
    Ok((r, s))
}

fn decode_integer(input: &[u8]) -> Result<(BigUint, &[u8]), DerError> {
    if input.len() < 2 {
        return Err(DerError::Truncated);
    }
    if input[0] != 0x02 {
        return Err(DerError::UnexpectedTag(input[0]));
    }
    let len = input[1];
    if len & 0x80 != 0 {
        return Err(DerError::BadLength);
    }
    let len = len as usize;
    if input.len() < 2 + len {
        return Err(DerError::Truncated);
    }
    let bytes = &input[2..2 + len];
    let rest = &input[2 + len..];

    if bytes.is_empty() {
        return Err(DerError::NonCanonical("zero-length INTEGER"));
    }
    // Reject high-bit set (would be negative in DER).
    if bytes[0] & 0x80 != 0 {
        return Err(DerError::NonCanonical(
            "INTEGER has high bit set (negative)",
        ));
    }
    // Reject leading zero unless required by next-byte high bit.
    if bytes.len() > 1 && bytes[0] == 0x00 && bytes[1] & 0x80 == 0 {
        return Err(DerError::NonCanonical("non-minimal INTEGER encoding"));
    }
    Ok((BigUint::from_bytes_be(bytes), rest))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Round-trip a small (r, s) pair.
    #[test]
    fn roundtrip_small() {
        let r = BigUint::from(42u32);
        let s = BigUint::from(7u32);
        let der = encode_ecdsa_signature(&r, &s).unwrap();
        let (r2, s2) = decode_ecdsa_signature(&der).unwrap();
        assert_eq!(r, r2);
        assert_eq!(s, s2);
    }

    /// Round-trip 256-bit values.
    #[test]
    fn roundtrip_256bit() {
        let r = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364140",
            16,
        )
        .unwrap();
        let s = BigUint::from(1u32);
        let der = encode_ecdsa_signature(&r, &s).unwrap();
        // r has high bit set ⇒ should have 0x00 prefix in encoding ⇒
        // r-length byte = 33.
        assert_eq!(der[2], 0x02); // INTEGER tag
        assert_eq!(der[3], 33); // length with prefix byte
        assert_eq!(der[4], 0x00); // sign-prefix
        let (r2, s2) = decode_ecdsa_signature(&der).unwrap();
        assert_eq!(r, r2);
        assert_eq!(s, s2);
    }

    /// **Real Bitcoin signature decode.**  This is a standard
    /// BIP-66 test vector (a well-formed mainnet signature in
    /// strict-DER form).
    #[test]
    fn decodes_real_bitcoin_signature() {
        // Signature from a real Bitcoin transaction (well-formed
        // strict-DER, BIP-66 compliant).
        let der = hex_decode(
            "30440220181522ec8eca07de4860a4acdd12909d831cc56cbbac46\
             22082221a8868d260902200ed99fcfaffe9c25b8eecf95efa6dab\
             dde6acff05bf60ce0eb45e10ada04b41401",
        );
        // The trailing 0x01 is the sighash-type byte appended by
        // Bitcoin — not part of the DER signature itself.  Strip
        // it before decoding.
        let sig_only = &der[..der.len() - 1];
        let (r, s) = decode_ecdsa_signature(sig_only).unwrap();
        assert!(r.bits() > 0 && r.bits() <= 256);
        assert!(s.bits() > 0 && s.bits() <= 256);
    }

    /// Strict-DER rejection: non-minimal encoding.
    #[test]
    fn rejects_non_minimal_integer() {
        // SEQUENCE { INTEGER 00 00 01, INTEGER 01 } — non-minimal
        // (leading 00 not required since 0x00 has high bit clear).
        let bad = vec![0x30, 0x08, 0x02, 0x03, 0x00, 0x00, 0x01, 0x02, 0x01, 0x01];
        let result = decode_ecdsa_signature(&bad);
        assert!(matches!(result, Err(DerError::NonCanonical(_))));
    }

    /// Strict-DER rejection: high-bit-set INTEGER without 0x00 prefix.
    #[test]
    fn rejects_negative_integer() {
        // SEQUENCE { INTEGER 80, INTEGER 01 } — 0x80 has high bit
        // set, would be negative without the 0x00 prefix.
        let bad = vec![0x30, 0x06, 0x02, 0x01, 0x80, 0x02, 0x01, 0x01];
        let result = decode_ecdsa_signature(&bad);
        assert!(matches!(result, Err(DerError::NonCanonical(_))));
    }

    /// Strict-DER rejection: zero r component.
    #[test]
    fn rejects_zero_component() {
        // SEQUENCE { INTEGER 00, INTEGER 01 }
        let bad = vec![0x30, 0x06, 0x02, 0x01, 0x00, 0x02, 0x01, 0x01];
        let result = decode_ecdsa_signature(&bad);
        // Either ZeroComponent or NonCanonical (depending on which
        // gate trips first; both are correct refusals).
        assert!(result.is_err());
    }

    /// Trailing bytes are rejected.
    #[test]
    fn rejects_trailing_bytes() {
        let r = BigUint::from(1u32);
        let s = BigUint::from(1u32);
        let mut der = encode_ecdsa_signature(&r, &s).unwrap();
        der.push(0xFF);
        // After a clean encoding, an extra byte makes the SEQUENCE
        // length mismatch the inner contents.  The decoder should
        // reject.
        assert!(decode_ecdsa_signature(&der).is_err());
    }

    fn hex_decode(s: &str) -> Vec<u8> {
        let s: String = s.chars().filter(|c| !c.is_whitespace()).collect();
        (0..s.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
            .collect()
    }
}
