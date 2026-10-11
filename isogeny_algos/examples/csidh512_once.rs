//! One CSIDH-512 action with exponents in [-5, 5]^74 (for profiling).
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::path::csidh::Csidh;

fn main() {
    let cs = Csidh::<FpM<8>>::csidh512();
    let mut rng = Rng::new(10_600);
    let e: Vec<i32> = (0..74).map(|_| rng.below(11) as i32 - 5).collect();
    let a = cs.action_fast(cs.fp.zero(), &e, &mut rng);
    println!("{:?}", cs.fp.to_big(&a).to_dec());
}
