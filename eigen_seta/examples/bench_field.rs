//! Micro-benchmark: field mul, Fermat inversion vs safegcd inversion.
use eigen_seta::fp::{Fp, FpConfig};
use eigen_seta::params::*;
use std::hint::black_box;
use std::time::Instant;

fn bench<C: FpConfig<N>, const N: usize>() {
    type F<C, const N: usize> = Fp<C, N>;
    let mut x = F::<C, N>::from_u64(123456789);
    let it = 20000;
    let mut mul = f64::MAX;
    for _ in 0..15 {
        let t = Instant::now();
        for _ in 0..it { x = black_box(x * x + F::<C, N>::one()); }
        mul = mul.min(t.elapsed().as_nanos() as f64 / it as f64);
    }
    let it = 200;
    let t = Instant::now();
    for _ in 0..it { x = black_box(x.invert_fermat() + F::<C, N>::one()); }
    let fer = t.elapsed().as_nanos() as f64 / it as f64;
    let t = Instant::now();
    for _ in 0..it { x = black_box(x.invert() + F::<C, N>::one()); }
    let sg = t.elapsed().as_nanos() as f64 / it as f64;
    println!("{:6} limbs {:2}  mul {:8.0} ns  fermat {:10.0} ns  safegcd {:9.0} ns  ratio {:5.1}  safegcd/mul {:5.0}",
        C::NAME, N, mul, fer, sg, fer / sg, sg / mul);
}

fn main() {
    bench::<Toy16Fp, 4>();
    bench::<Toy32Fp, 6>();
    bench::<L1Fp, 21>();
    bench::<L3Fp, 32>();
    bench::<L5Fp, 42>();
}
