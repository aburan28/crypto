//! Known-answer tests: key generation, encapsulation and decapsulation from
//! a fixed deterministic RNG seed must reproduce pinned SHA3-256 digests of
//! the public key, ciphertext and shared secret.  These are regression KATs
//! for this implementation (no second implementation exists yet); a change
//! to any of them is a change of the scheme's wire behaviour and must be
//! deliberate.  Set KAT_PRINT=1 to print fresh values.
mod common;
use common::testrng::TestRng;
use eigen_seta::kem::EigenSeta;
use eigen_seta::params::{ParamSet, Toy16, Toy32, L1};
use sha3::{Digest, Sha3_256};

fn hexd(b: &[u8]) -> String {
    Sha3_256::digest(b).iter().map(|x| format!("{x:02x}")).collect()
}

fn kat<P: ParamSet<N>, const N: usize>(seed: u64, want: [&str; 3]) {
    let mut rng = TestRng(seed);
    let (pk, sk) = EigenSeta::<P, N>::keygen(&mut rng);
    let (ct, ss) = EigenSeta::<P, N>::encaps(&pk, &mut rng).unwrap();
    assert_eq!(EigenSeta::<P, N>::decaps(&sk, &ct), ss, "{}: KAT decapsulation", P::NAME);
    let got = [hexd(&pk.bytes), hexd(&ct.bytes), hexd(&ss)];
    if std::env::var("KAT_PRINT").is_ok() {
        eprintln!("KAT {} seed {seed:#x}: pk {} ct {} ss {}", P::NAME, got[0], got[1], got[2]);
        return;
    }
    for (i, what) in ["pk", "ct", "ss"].iter().enumerate() {
        assert_eq!(got[i], want[i], "{} KAT mismatch in {what}", P::NAME);
    }
}

#[test]
fn kat_toy16() {
    kat::<Toy16, 4>(0x4b41_5431, ["0572c4e2914a0731d28e79492498940c1f68d475b28bf1b869611a662b834c50", "23274e0fd49da2d853cbcaa8585f1b25fe0b8bd0855ba93ab16c9e75010c1768", "7f4d278b93b08600cd729b74a8620e3c679d50954667d1d292c5f0b5b38b70c4"]);
}

#[test]
fn kat_toy32() {
    kat::<Toy32, 6>(0x4b41_5432, ["b2c77814e60221e216aa07ee7f9e982af7aa6e12b6de4077d1d0e4b504c05504", "1c4f32b175779233a7b2abf8c70134a934d9af85a49d58e958525df0ae01604d", "1c7d20a8f878ffbe8a08a3defb3fa5f3d82d6107be87f20689a2b8609490771e"]);
}

#[test]
#[ignore] // about a minute in release mode
fn kat_l1() {
    kat::<L1, 21>(0x4b41_5433, ["f225592743478273085120c6bdbd8946b4351526215220a24ac213ae7864a7e1", "568aa2c8fc225150610cca3523dfabe547c75c921e2422e0727a2af14b3b624c", "db17a4189fa2cad1c9597433e954849b568ebf2c0496d627d6e9ad74ad8a95dc"]);
}
