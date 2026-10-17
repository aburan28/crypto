//! Long polynomial products over multi-limb prime fields: Karatsuba (`poly::mul_raw` with the NTT
//! disabled) against the multi-prime NTT (`FpM::conv_long`), checked equal, timed per length.
//!
//! `cargo run --release --example profile_ntt`
use isogeny_algos::field::*;
use isogeny_algos::fpm::{FpM, NTT_N};
use isogeny_algos::poly;
use std::time::Instant;

fn time<T>(reps: usize, mut f: impl FnMut() -> T) -> (f64, T) {
    let mut best = f64::MAX;
    let mut out = None;
    for _ in 0..reps {
        let t = Instant::now();
        let r = f();
        best = best.min(t.elapsed().as_secs_f64() * 1e6);
        out = Some(r);
    }
    (best, out.unwrap())
}

fn run<const N: usize>(name: &str, f: &FpM<N>, lens: &[usize]) {
    let mut rng = Rng::new(7);
    println!("# {name}: len  karatsuba_us  ntt_us  ratio");
    for &n in lens {
        let a: Vec<[u64; N]> = (0..n).map(|_| f.random(&mut rng)).collect();
        let b: Vec<[u64; N]> = (0..n).map(|_| f.random(&mut rng)).collect();
        let reps = if n < 512 {
            30
        } else if n < 4096 {
            5
        } else {
            2
        };
        NTT_N.with(|k| k.set(Some(usize::MAX)));
        let (tk, rk) = time(reps, || poly::mul_raw(f, &a, &b));
        NTT_N.with(|k| k.set(Some(0)));
        let (tn, rn) = time(reps, || poly::mul_raw(f, &a, &b));
        NTT_N.with(|k| k.set(None));
        assert_eq!(rk, rn, "{name} len {n}: NTT product differs");
        println!("{n:6} {tk:12.1} {tn:9.1} {:6.2}", tk / tn);
    }
}

fn main() {
    let lens = [
        16, 24, 32, 40, 48, 64, 80, 96, 128, 160, 192, 256, 320, 512, 1024,
    ];
    run(
        "FpM1 (2^61-1)",
        &FpM::<1>::from_dec("2305843009213693951"),
        &lens,
    );
    run(
        "FpM2 (2^127-1)",
        &FpM::<2>::from_dec("170141183460469231731687303715884105727"),
        &lens,
    );
    run(
        "FpM4 (2^255-19)",
        &FpM::<4>::from_dec(
            "57896044618658097711785492504343953926634992332820282019728792003956564819949",
        ),
        &lens,
    );
    let csidh = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";
    run("FpM8 (CSIDH-512)", &FpM::<8>::from_dec(csidh), &lens);
}
