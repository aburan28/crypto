//! Byte/hex/BigUint encoding utilities.

use num_bigint::BigUint;

/// Encode bytes as lowercase hex string.
pub fn to_hex(bytes: &[u8]) -> String {
    hex::encode(bytes)
}

/// Decode a hex string into bytes, returning an error on bad input.
pub fn from_hex(s: &str) -> Result<Vec<u8>, hex::FromHexError> {
    hex::decode(s)
}

/// Encode a `BigUint` as a fixed-width big-endian byte array.
/// Pads with leading zeros if the value is shorter than `len`.
pub fn bigint_to_bytes_be(n: &BigUint, len: usize) -> Vec<u8> {
    let raw = n.to_bytes_be();
    if raw.len() >= len {
        raw[raw.len() - len..].to_vec()
    } else {
        let mut out = vec![0u8; len];
        out[len - raw.len()..].copy_from_slice(&raw);
        out
    }
}

/// Decode a big-endian byte slice into a `BigUint`.
pub fn bytes_be_to_bigint(bytes: &[u8]) -> BigUint {
    BigUint::from_bytes_be(bytes)
}
