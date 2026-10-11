//! sqrt-Velu (Montgomery) codomain at l = 10007 over a 511-bit prime, repeated (for profiling).
use isogeny_algos::field::*;
use isogeny_algos::kernel::montgomery as mg;
use isogeny_algos::kernel::sqrt_velu_mont::SqrtVeluMont;

fn main() {
    let mut rng = Rng::new(5);
    let w = isogeny_algos::testdata::supersingular_workload::<8>(&[3, 10007], 511, &mut rng);
    let f = &w.f;
    let a = f.zero();
    let k24 = mg::proj24(f, a);
    let cof = f.modulus().add_small(1).divrem_small(10007).0;
    let kp = loop {
        let k = mg::ladder_p(f, k24, (f.random(&mut rng), f.one()), &cof);
        if !f.is_zero(k.1) {
            break k;
        }
    };
    let mut acc = f.zero();
    for _ in 0..5 {
        let sv = SqrtVeluMont::new(f, a, kp, 10007).unwrap();
        acc = f.add(acc, sv.codomain(f));
    }
    println!("{:?}", f.to_big(&acc).to_dec().len());
}
