use eigen_seta::fp::{Fp, FpConfig};
use std::hint::black_box;
use std::time::Instant;
macro_rules! cfg {
    ($name:ident, $n:expr) => {
        #[derive(Clone, Copy, Debug)]
        struct $name;
        impl FpConfig<$n> for $name {
            const MODULUS: [u64; $n] = {
                let mut m = [u64::MAX; $n];
                m[$n - 1] = u64::MAX >> 1;
                m[0] = u64::MAX - 58;
                m
            };
            const NAME: &'static str = stringify!($name);
        }
    };
}
cfg!(C16, 16);
cfg!(C20, 20);
cfg!(C21, 21);
cfg!(C22, 22);
cfg!(C24, 24);
cfg!(C26, 26);
cfg!(C28, 28);
cfg!(C32, 32);
cfg!(C42, 42);
fn b<C: FpConfig<N>, const N: usize>() {
    let mut x = Fp::<C, N>::from_u64(1234567);
    let mut best = f64::MAX;
    for _ in 0..10 {
        let t = Instant::now();
        for _ in 0..10000 {
            x = black_box(x * x + Fp::<C, N>::one());
        }
        best = best.min(t.elapsed().as_nanos() as f64 / 10000.0);
    }
    println!(
        "N={:2} mul {:6.0} ns  ns/N^2 {:.3}",
        N,
        best,
        best / (N * N) as f64
    );
}
fn main() {
    b::<C16, 16>();
    b::<C20, 20>();
    b::<C21, 21>();
    b::<C22, 22>();
    b::<C24, 24>();
    b::<C26, 26>();
    b::<C28, 28>();
    b::<C32, 32>();
    b::<C42, 42>();
}
