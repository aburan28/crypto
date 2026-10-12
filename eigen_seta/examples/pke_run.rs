//! cargo run --release --example pke_run -- <Toy16|Toy32|L1|L3|L5>
use eigen_seta::kem::EigenSeta;
use eigen_seta::params::*;
use eigen_seta::pke::EigenSetaPke;
use rand_core::OsRng;
use std::time::Instant;

fn go<P: ParamSet<N>, const N: usize>() {
    let mut rng = OsRng;
    let message = b"Eigen-Seta authenticated PKE round trip";
    let aad = b"pke-run-v1";
    let t0 = Instant::now();
    let (pk, sk) = EigenSeta::<P, N>::keygen(&mut rng);
    let keygen_s = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let ct = EigenSetaPke::<P, N>::encrypt(&pk, message, aad, &mut rng).unwrap();
    let encrypt_s = t1.elapsed().as_secs_f64();
    let t2 = Instant::now();
    let recovered = EigenSetaPke::<P, N>::decrypt(&sk, &ct, aad).unwrap();
    let decrypt_s = t2.elapsed().as_secs_f64();
    let correct = recovered == message;
    println!(
        "{} keygen_s {:.3} encrypt_s {:.4} decrypt_s {:.3} pk_bytes {} pke_ct_bytes {} message_bytes {} correct {}",
        P::NAME,
        keygen_s,
        encrypt_s,
        decrypt_s,
        pk.bytes.len(),
        ct.bytes.len(),
        message.len(),
        correct
    );
    assert!(correct);
}

fn main() {
    let set = std::env::args().nth(1).unwrap_or_else(|| "Toy16".into());
    match set.as_str() {
        "Toy16" => go::<Toy16, 4>(),
        "Toy32" => go::<Toy32, 6>(),
        "L1" => go::<L1, 21>(),
        "L3" => go::<L3, 32>(),
        "L5" => go::<L5, 42>(),
        _ => panic!("unknown set"),
    }
}
