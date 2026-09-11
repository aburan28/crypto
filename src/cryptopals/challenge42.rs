//! # Challenge 42 — Bleichenbacher's e=3 forgery
//!
//! With a sloppy PKCS#1-v1.5 signature verifier that doesn't check
//! padding bytes to the end of the modulus, you can forge an RSA
//! signature without knowing the private key.  Construct a value
//! whose cube starts with `00 01 FF .. 00 ASN.1 H 0…` — the integer
//! cube root of `2^(modbits - …)` works.

use crate::cryptopals::Report;
use crate::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::One;

fn int_cube_root_ceil(n: &BigUint) -> BigUint {
    // Returns smallest x with x^3 >= n.
    let three = BigUint::from(3u32);
    let two = BigUint::from(2u32);
    let mut x = BigUint::one() << ((n.bits() / 3) as usize + 2);
    loop {
        let next = (&x * &two + n / (&x * &x)) / &three;
        if next >= x {
            // Adjust to ceil.
            if &x * &x * &x < *n {
                return &x + BigUint::one();
            }
            return x;
        }
        x = next;
    }
}

/// Sloppy verifier: parses 0x00 0x01 FF... 00 <ASN.1+hash> then
/// returns success without scanning to the end.
pub fn sloppy_verify(s: &BigUint, n: &BigUint, e: &BigUint, msg: &[u8]) -> bool {
    let cube = s.modpow(e, n);
    let modulus_bytes = (n.bits() as usize + 7) / 8;
    let mut bytes = cube.to_bytes_be();
    while bytes.len() < modulus_bytes {
        bytes.insert(0, 0);
    }
    if bytes[0] != 0x00 || bytes[1] != 0x01 {
        return false;
    }
    // Skip FFs.
    let mut i = 2;
    while i < bytes.len() && bytes[i] == 0xff {
        i += 1;
    }
    if i >= bytes.len() || bytes[i] != 0x00 {
        return false;
    }
    i += 1;
    let asn1: [u8; 19] = [
        0x30, 0x31, 0x30, 0x0d, 0x06, 0x09, 0x60, 0x86, 0x48, 0x01, 0x65, 0x03, 0x04,
        0x02, 0x01, 0x05, 0x00, 0x04, 0x20,
    ];
    if bytes.len() < i + asn1.len() + 32 {
        return false;
    }
    if &bytes[i..i + asn1.len()] != asn1 {
        return false;
    }
    i += asn1.len();
    let h = sha256(msg);
    bytes[i..i + 32] == h
}

pub fn run() -> Report {
    let mut r = Report::new(42, "Bleichenbacher e=3 signature forgery");
    // Pick a 3072-bit modulus so the cube root forgery has plenty of room.
    let kp = crate::asymmetric::rsa::RsaKeyPair::generate(3072);
    let n = kp.public.n.clone();
    let msg = b"hi mom";
    // Construct a "block" that the lax parser accepts but with the
    // hash near the beginning, padded with zeros on the right.
    let modulus_bytes = (n.bits() as usize + 7) / 8;
    let asn1: [u8; 19] = [
        0x30, 0x31, 0x30, 0x0d, 0x06, 0x09, 0x60, 0x86, 0x48, 0x01, 0x65, 0x03, 0x04,
        0x02, 0x01, 0x05, 0x00, 0x04, 0x20,
    ];
    let h = sha256(msg);
    let mut block = vec![0u8; modulus_bytes];
    block[0] = 0x00;
    block[1] = 0x01;
    // Minimal pad: a few 0xff bytes then 0x00.
    block[2] = 0xff;
    block[3] = 0xff;
    block[4] = 0xff;
    block[5] = 0xff;
    block[6] = 0xff;
    block[7] = 0xff;
    block[8] = 0xff;
    block[9] = 0x00;
    block[10..10 + asn1.len()].copy_from_slice(&asn1);
    block[10 + asn1.len()..10 + asn1.len() + 32].copy_from_slice(&h);
    // The rest stays zero — sloppy parser doesn't check.
    let target = BigUint::from_bytes_be(&block);
    let forged = int_cube_root_ceil(&target);
    let ok = sloppy_verify(&forged, &n, &BigUint::from(3u32), msg);
    r.line(format!("Sloppy verify accepts forged sig (e=3): {}", ok));
    if !ok {
        // Fall back: try +1 and -1 if rounding tripped us.
        let alt1 = &forged + BigUint::one();
        let alt2 = &forged - BigUint::one();
        let _ = (alt1, alt2);
    }
    assert!(ok);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn forges_e3_signature() {
        assert!(super::run().success);
    }
}
