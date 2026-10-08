use isogeny_algos::icv1::prime;
use isogeny_algos::int::Int;
use isogeny_algos::sha256::sha256_hex;

/// FIPS 180-4 example vectors.
#[test]
fn sha256_known_answers() {
    assert_eq!(sha256_hex(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    assert_eq!(sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    assert_eq!(
        sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
    );
    let million_a = vec![b'a'; 1_000_000];
    assert_eq!(sha256_hex(&million_a), "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

/// The prime-field vectors that pin the crypto crate's port (`tests/curve_id.rs` there) to the
/// reference implementation (`scripts/curve_id.py`).
#[test]
fn icv1_prime_matches_the_reference_vectors() {
    let i = |v: i64| Int::from(v);
    let id = prime(&i(10_935_329), &i(5_320_418), &i(8_535_318), &i(10_933_753)).unwrap();
    assert_eq!(id.icv1, "ICV1:fp-10935329:1577:10933753:2786525:unk:unk:r:773361558ca8");
    assert_eq!(id.slug, "icv1-fp24-t1577-77336155");
    let id = prime(&i(827), &i(1), &i(15), &i(823)).unwrap();
    assert_eq!(id.slug, "icv1-fp10-t5-192cb216");
    // coefficients are reduced mod p before hashing; a negative a names the same model
    let id2 = prime(&i(827), &i(1 - 827), &i(15), &i(823)).unwrap();
    assert_eq!(id2, id);
    // a subgroup order passed as the group's breaks the Hasse bound
    assert!(prime(&i(827), &i(1), &i(15), &i(400)).is_none());
}
