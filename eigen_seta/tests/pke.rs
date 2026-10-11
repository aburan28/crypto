use eigen_seta::kem::{EigenSeta, PublicKey};
use eigen_seta::params::Toy16;
use eigen_seta::pke::{EigenSetaPke, PkeCiphertext, NONCE_BYTES, TAG_BYTES};
use rand_core::OsRng;

#[test]
fn pke_round_trip_and_rejection() {
    let mut rng = OsRng;
    let (pk, sk) = EigenSeta::<Toy16, 4>::keygen(&mut rng);
    let aad = b"message context";
    assert!(
        EigenSetaPke::<Toy16, 4>::decrypt(&sk, &PkeCiphertext { bytes: Vec::new() }, aad).is_none()
    );
    let invalid_pk = PublicKey {
        bytes: vec![0; EigenSeta::<Toy16, 4>::pk_bytes()],
    };
    assert!(EigenSetaPke::<Toy16, 4>::encrypt(&invalid_pk, b"test", aad, &mut rng).is_none());
    for message in [&b""[..], &b"Eigen-Seta end-to-end PKE"[..]] {
        let ct = EigenSetaPke::<Toy16, 4>::encrypt(&pk, message, aad, &mut rng).unwrap();
        assert_eq!(
            ct.bytes.len(),
            EigenSeta::<Toy16, 4>::ct_bytes() + NONCE_BYTES + message.len() + TAG_BYTES
        );
        assert_eq!(
            EigenSetaPke::<Toy16, 4>::decrypt(&sk, &ct, aad).as_deref(),
            Some(message)
        );
        assert!(EigenSetaPke::<Toy16, 4>::decrypt(&sk, &ct, b"wrong context").is_none());
        for index in [0, EigenSeta::<Toy16, 4>::ct_bytes(), ct.bytes.len() - 1] {
            let mut tampered = ct.clone();
            tampered.bytes[index] ^= 1;
            assert!(EigenSetaPke::<Toy16, 4>::decrypt(&sk, &tampered, aad).is_none());
        }
        let truncated = PkeCiphertext {
            bytes: ct.bytes[..ct.bytes.len() - 1].to_vec(),
        };
        assert!(EigenSetaPke::<Toy16, 4>::decrypt(&sk, &truncated, aad).is_none());
    }
}
