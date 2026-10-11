//! Constant-time probe.  Generates keys and ciphertexts from fixed seeds, then
//! runs the secret phase of one decapsulation in `probe_core` (never inlined),
//! so that `valgrind --tool=callgrind --toggle-collect='*probe_core*'` reports
//! the instruction count of that phase alone.  Equal counts across secret keys
//! and across valid vs. FO-rejected ciphertexts show that control flow does
//! not depend on secret data (memory addresses are covered separately).
//!
//!   ct_probe <key_seed> <ct_seed> <valid|foreign> [Toy16|Toy32|L1|L3|L5]
use eigen_seta::kem::{Ciphertext, DecapsPrep, EigenSeta, SecretKey};
use eigen_seta::params::*;
use rand_core::{CryptoRng, RngCore};

/// SplitMix64: deterministic test RNG (NOT for production keys).
struct TestRng(u64);
impl RngCore for TestRng {
    fn next_u32(&mut self) -> u32 { self.next_u64() as u32 }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn fill_bytes(&mut self, d: &mut [u8]) { for c in d.chunks_mut(8) { let v = self.next_u64().to_le_bytes(); c.copy_from_slice(&v[..c.len()]); } }
    fn try_fill_bytes(&mut self, d: &mut [u8]) -> Result<(), rand_core::Error> { self.fill_bytes(d); Ok(()) }
}
impl CryptoRng for TestRng {}

type Fl<P, const N: usize> = eigen_seta::fp2::Fp2<<P as ParamSet<N>>::Cfg, N>;

/// Only the secret phase is measured: the public preparation (parsing,
/// point lifting, canonical bases; variable time by design) runs outside.
#[inline(never)]
fn probe_core<P: ParamSet<N>, const N: usize>(sk: &SecretKey, prep: &DecapsPrep<Fl<P, N>>, ct: &Ciphertext) -> [u8; 32] {
    EigenSeta::<P, N>::decaps_core(std::hint::black_box(sk), std::hint::black_box(prep), std::hint::black_box(ct))
}

fn go<P: ParamSet<N>, const N: usize>(ks: u64, cs: u64, mode: &str) {
    type K<P, const N: usize> = EigenSeta<P, N>;
    let (pk, sk) = K::<P, N>::keygen(&mut TestRng(1000 + ks));
    let mut crng = TestRng(2000 + cs);
    let (ct, k) = match mode {
        "valid" => K::<P, N>::encaps(&pk, &mut crng).unwrap(),
        // a well-formed ciphertext for a different key: passes parsing and runs
        // the full trapdoor, then fails the FO re-encryption check
        "foreign" => {
            let (pk2, _) = K::<P, N>::keygen(&mut TestRng(9000 + ks));
            K::<P, N>::encaps(&pk2, &mut crng).unwrap()
        }
        _ => panic!("mode"),
    };
    let prep = K::<P, N>::decaps_prepare(&sk, &ct).expect("well-formed ciphertext");
    let k2 = probe_core::<P, N>(&sk, &prep, &ct);
    println!("{} key {ks} ct {cs} {mode} accepted {}", P::NAME, k == k2);
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let ks: u64 = a[1].parse().unwrap();
    let cs: u64 = a[2].parse().unwrap();
    let mode = a[3].as_str();
    match a.get(4).map(|s| s.as_str()).unwrap_or("Toy16") {
        "Toy16" => go::<Toy16, 4>(ks, cs, mode),
        "Toy32" => go::<Toy32, 6>(ks, cs, mode),
        "L1" => go::<L1, 21>(ks, cs, mode),
        "L3" => go::<L3, 32>(ks, cs, mode),
        "L5" => go::<L5, 42>(ks, cs, mode),
        _ => panic!("unknown set"),
    }
}
