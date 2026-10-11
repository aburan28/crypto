//! Reproducible small-parameter Eigen-Séta experiment.
//! Run: cargo run --release --example eigen_seta_toy
use isogeny_algos::eigen_seta_toy::{
    toy_decaps, toy_encaps, toy_kani_first_half, toy_keygen, toy_public_recover, ToyParameters,
};
use isogeny_algos::field::Rng;

fn main() {
    let params = ToyParameters {
        prime: 16127,
        two_exponent: 6,
        three_exponent: 2,
    };
    let mut rng = Rng::new(0x44c0);
    let sk = toy_keygen(params, &mut rng).expect("CM toy key generation");
    let (ct, sent) = toy_encaps(&sk.public, &mut rng).expect("toy encapsulation");
    let received = toy_decaps(&sk, &ct).expect("toy decapsulation");
    let public = toy_public_recover(&sk.public, &ct).expect("public exhaustive recovery");
    assert_eq!(sent, received);
    assert_eq!(sent, public);
    let kani = toy_kani_first_half(&sk, &ct).expect("first Kani half");
    println!("toy KEM round trip: PASS");
    println!("public exhaustive recovery: PASS (insecure toy)");
    println!(
        "first Kani half: {} product steps + {} theta steps, split={:?}",
        kani.product_steps,
        kani.first_half.nulls.len(),
        kani.first_half.split
    );
}
