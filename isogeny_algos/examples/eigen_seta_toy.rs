//! Reproducible small-parameter Eigen-Séta experiment.
//! Run: cargo run --release --example eigen_seta_toy
use isogeny_algos::eigen_seta_toy::{
    toy_decaps, toy_encaps, toy_kani_dual_first_half, toy_kani_first_half, toy_keygen,
    toy_public_recover, toy_trapdoor_decaps, ToyParameters,
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
    let oracle = toy_trapdoor_decaps(&sk, &ct).expect("two-halves toy oracle decapsulation");
    assert_eq!(sent, received);
    assert_eq!(sent, public);
    assert_eq!(sent, oracle);
    let kani = toy_kani_first_half(&sk, &ct).expect("first Kani half");
    let dual = toy_kani_dual_first_half(&sk, &ct).expect("dual first half");
    println!("toy KEM round trip: PASS");
    println!("public exhaustive recovery: PASS (insecure toy)");
    println!("two-halves odd-torsion oracle recovery: PASS (exponential toy)");
    println!(
        "first Kani half: {} product steps + {} theta steps, split={:?}",
        kani.product_steps,
        kani.first_half.nulls.len(),
        kani.first_half.split
    );
    println!(
        "dual first half: {} product steps + {} theta steps, split={:?}",
        dual.product_steps,
        dual.first_half.nulls.len(),
        dual.first_half.split
    );

    let special = ToyParameters {
        prime: 21503,
        two_exponent: 7,
        three_exponent: 1,
    };
    let special_sk = toy_keygen(special, &mut rng).expect("special-key toy keygen");
    let (special_ct, special_sent) =
        toy_encaps(&special_sk.public, &mut rng).expect("special-key toy encapsulation");
    let special_received = toy_trapdoor_decaps(&special_sk, &special_ct)
        .expect("secret-dependent single-chain toy decapsulation");
    assert_eq!(special_sent, special_received);
    println!("secret-dependent single-chain toy round trip: PASS (insecure CM key)");
}
