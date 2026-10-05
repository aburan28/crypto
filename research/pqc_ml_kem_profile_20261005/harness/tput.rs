//! Throughput loop matching pqc_tput.c: same key, N calls, total rdtsc / N.
use crypto_lib::pqc::ml_kem::ML_KEM_768;
use mlkem_dev::ml_kem::*;
use std::hint::black_box;
fn rd() -> u64 { unsafe { core::arch::x86_64::_rdtsc() } }
const N: u64 = 20000;
fn main() {
    let p = ML_KEM_768;
    let mut coins = [7u8; 32];
    let (ek, dk) = ml_kem_keygen_internal(&p, &coins, &coins);
    let (ct, _) = ml_kem_encaps_internal(&p, &ek, &coins).unwrap();
    for _ in 0..1000 { black_box(ml_kem_encaps_internal(&p, &ek, &coins)); }
    let t = rd();
    for i in 0..N { coins[0] = i as u8; black_box(ml_kem_keygen_internal(&p, &coins, &coins)); }
    let kg = (rd() - t) / N;
    let t = rd();
    for i in 0..N { coins[0] = i as u8; black_box(ml_kem_encaps_internal(&p, &ek, &coins)); }
    let en = (rd() - t) / N;
    let t = rd();
    for _ in 0..N { black_box(ml_kem_decaps(&p, &dk, &ct)); }
    let de = (rd() - t) / N;
    println!("keygen\t{kg}\nencaps\t{en}\ndecaps\t{de}");
}
