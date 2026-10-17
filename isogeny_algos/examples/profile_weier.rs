//! Weierstrass sqrt-Velu (pairs) vs x-only Velu at large l (scratch profiling).
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::kernel::sqrt_velu::sqrt_velu_fast;
use isogeny_algos::kernel::xonly::velu_xonly_fast;
use std::time::Instant;

fn med(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn run<const N: usize>(bits: usize, ells: &[u64]) {
    let mut rng = Rng::new(14_000);
    let w = isogeny_algos::testdata::supersingular_workload::<N>(ells, bits, &mut rng);
    let f: &FpM<N> = &w.f;
    for &(ell, p) in &w.pts {
        if ell < 100 {
            continue;
        }
        let Pt::Aff(x0, _) = p else { unreachable!() };
        let t = random_point_f(f, &w.e, &mut rng);
        let reps = if ell > 3000 { 7 } else { 21 };
        let (mut a, mut b, mut c, mut d) = (vec![], vec![], vec![], vec![]);
        for _ in 0..reps {
            let s = Instant::now();
            let v = velu_xonly_fast(f, &w.e, x0, ell).unwrap();
            a.push(s.elapsed().as_nanos() as f64);
            let s = Instant::now();
            std::hint::black_box(v.eval(f, &t));
            b.push(s.elapsed().as_nanos() as f64);
            let s = Instant::now();
            let q = sqrt_velu_fast(f, &w.e, x0, ell).unwrap();
            c.push(s.elapsed().as_nanos() as f64);
            let s = Instant::now();
            std::hint::black_box(q.eval(f, &t));
            d.push(s.elapsed().as_nanos() as f64);
        }
        let (a, b, c, d) = (med(a) / 1e3, med(b) / 1e3, med(c) / 1e3, med(d) / 1e3);
        println!(
            "{bits} l={ell:6}: codomain xonly {a:8.0} sqrt {c:8.0} us {:5.2}x | image xonly {b:7.0} sqrt {d:7.0} us {:5.2}x",
            a / c,
            b / d
        );
    }
}

fn main() {
    let ells = [3u64, 101, 401, 1009, 4001, 10007];
    run::<4>(256, &ells);
    run::<8>(511, &ells);
}
