use eigen_seta::kem::{Ciphertext, EigenSeta};
use eigen_seta::params::{Toy16, Toy32};
use rand_core::OsRng;
use std::time::Instant;

fn roundtrip<P: eigen_seta::params::ParamSet<N>, const N: usize>(trials: usize) {
    let mut rng = OsRng;
    let t0 = Instant::now();
    let (pk, sk) = EigenSeta::<P, N>::keygen(&mut rng);
    eprintln!(
        "[{}] keygen {:?}; pk {} B, ct {} B; secret degree has {} bits",
        P::NAME,
        t0.elapsed(),
        pk.bytes.len(),
        EigenSeta::<P, N>::ct_bytes(),
        sk.secret_degree_bits()
    );
    for i in 0..trials {
        let t1 = Instant::now();
        let (ct, k1) = EigenSeta::<P, N>::encaps(&pk, &mut rng).unwrap();
        let te = t1.elapsed();
        let t2 = Instant::now();
        let k2 = EigenSeta::<P, N>::decaps(&sk, &ct);
        let td = t2.elapsed();
        assert_eq!(k1, k2, "decapsulation mismatch in trial {i}");
        eprintln!("[{}] trial {i}: encaps {:?}, decaps {:?}", P::NAME, te, td);
        // tampering: flip one byte of the ciphertext -> implicit rejection key
        let mut bad = ct.bytes.clone();
        let j = bad.len() - 3;
        bad[j] ^= 1;
        let k3 = EigenSeta::<P, N>::decaps(&sk, &Ciphertext { bytes: bad });
        assert_ne!(k3, k1, "tampered ciphertext accepted");
    }
}

#[test]
fn kem_toy16() {
    roundtrip::<Toy16, 4>(4);
}

#[test]
#[ignore]
fn kem_toy32() {
    roundtrip::<Toy32, 6>(2);
}
