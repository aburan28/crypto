//! Squaring against the general product of two equal vectors (FpM2, FpM4): ns per call.
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::poly;
use std::hint::black_box;
use std::time::Instant;

fn run<const N: usize>(name: &str, f: &FpM<N>) {
    let mut rng = Rng::new(1);
    println!("# {name}: n  mul_us  sqr_us  conv_us  conv_sq_us");
    for n in [8usize, 16, 24, 32, 48, 64, 90, 128, 160, 192] {
        let a: Vec<[u64; N]> = (0..n).map(|_| f.random(&mut rng)).collect();
        let b = a.clone();
        let reps = 300;
        let t = |g: &dyn Fn() -> Vec<[u64; N]>| {
            let mut best = f64::MAX;
            for _ in 0..5 {
                let s = Instant::now();
                for _ in 0..reps {
                    black_box(g());
                }
                best = best.min(s.elapsed().as_secs_f64() * 1e6 / reps as f64);
            }
            best
        };
        let m = t(&|| poly::mul_raw(f, &a, &b));
        let s = t(&|| poly::sqr_raw(f, &a));
        let c = t(&|| f.conv(&a, &b));
        let cs = t(&|| f.conv_sq(&a));
        println!("{n:4} {m:8.2} {s:8.2} {c:8.2} {cs:8.2}");
    }
}

fn main() {
    run(
        "FpM2",
        &FpM::<2>::from_dec("170141183460469231731687303715884105727"),
    );
    run(
        "FpM4",
        &FpM::<4>::from_dec(
            "57896044618658097711785492504343953926634992332820282019728792003956564819949",
        ),
    );
    let csidh = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";
    run("FpM8", &FpM::<8>::from_dec(csidh));
}
