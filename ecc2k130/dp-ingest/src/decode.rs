//! One 32-byte distinguished-point record → the columns of one store row.

pub const RECORD_BYTES: usize = 32;
/// 17-byte big-endian mod 2^131, per `rho_campaigns.meta`.
pub const COEFF_BYTES: usize = 17;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Decoded {
    /// The client's three little-endian uint64 keys, as packed.
    pub point_key: [u8; 24],
    /// Seed as 17 big-endian bytes.
    pub a: [u8; COEFF_BYTES],
    /// Always seventeen zero bytes here.
    pub b: [u8; COEFF_BYTES],
    /// Seed as big-endian *integer* bytes: no leading zeros (zero → `[0]`).
    pub walk_seed: Vec<u8>,
}

/// One 32-byte record → the column values for one row.
///
/// # Panics
///
/// When `record` is not exactly [`RECORD_BYTES`] long.
pub fn decode(record: &[u8]) -> Decoded {
    assert_eq!(
        record.len(),
        RECORD_BYTES,
        "record must be {RECORD_BYTES} bytes"
    );
    let seed = u64::from_le_bytes(record[0..8].try_into().unwrap());
    let mut point_key = [0u8; 24];
    point_key.copy_from_slice(&record[8..RECORD_BYTES]);
    let mut a = [0u8; COEFF_BYTES];
    a[COEFF_BYTES - 8..].copy_from_slice(&seed.to_be_bytes());
    let walk_seed = {
        let full = seed.to_be_bytes();
        let trimmed: Vec<u8> = full.iter().copied().skip_while(|&b| b == 0).collect();
        if trimmed.is_empty() {
            vec![0]
        } else {
            trimmed
        }
    };
    Decoded {
        point_key,
        a,
        b: [0u8; COEFF_BYTES],
        walk_seed,
    }
}

/// The walk seed as an integer, however wide the stored bytes happen to be.
pub fn seed_of(coeff: &[u8]) -> u128 {
    coeff
        .iter()
        .fold(0u128, |acc, &b| (acc << 8) | u128::from(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_key_is_the_last_24_bytes_of_the_record() {
        let record: Vec<u8> = (0..32).collect();
        assert_eq!(decode(&record).point_key.as_slice(), &record[8..]);
    }

    #[test]
    fn walk_seed_is_the_seed_big_endian_without_leading_zeros() {
        let mut record = [0u8; 32];
        record[..8].copy_from_slice(&7u64.to_le_bytes());
        assert_eq!(decode(&record).walk_seed, [0x07]);
    }

    #[test]
    fn coefficients_are_17_bytes() {
        let record: Vec<u8> = (0..32).collect();
        let r = decode(&record);
        assert_eq!((r.a.len(), r.b.len()), (17, 17));
        assert_eq!(r.b, [0u8; 17]);
    }

    #[test]
    fn seed_of_ignores_leading_width() {
        assert_eq!(seed_of(&[0x07]), 7);
        assert_eq!(seed_of(&[0x00, 0x00, 0x07]), 7);
        assert_eq!(seed_of(&[]), 0);
        assert_eq!(seed_of(&[0xff; 8]), u128::from(u64::MAX));
    }

    #[test]
    fn zero_seed_is_a_single_zero_byte() {
        assert_eq!(decode(&[0u8; 32]).walk_seed, [0]);
    }
}
