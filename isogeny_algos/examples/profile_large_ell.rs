//! Montgomery sqrt-Velu vs projective Velu at large l: codomain (setup + codomain) and one point
//! image (scratch profiling).
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::kernel::montgomery as mg;
use isogeny_algos::kernel::sqrt_velu_mont::SqrtVeluMont;
use std::time::Instant;

fn med(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn run<const N: usize>(bits: usize, ells: &[u64]) {
    let mut rng = Rng::new(5);
    let w = isogeny_algos::testdata::supersingular_workload::<N>(ells, bits, &mut rng);
    let f: &FpM<N> = &w.f;
    let a = f.zero();
    let k24 = mg::proj24(f, a);
    for &ell in ells {
        if ell < 100 {
            continue;
        }
        let cof = f.modulus().add_small(1).divrem_small(ell).0;
        let kp = loop {
            let k = mg::ladder_p(f, k24, (f.random(&mut rng), f.one()), &cof);
            if !f.is_zero(k.1) {
                break k;
            }
        };
        let u = f.random(&mut rng);
        let reps = if ell > 3000 { 9 } else { 31 };
        let (mut tv, mut tve, mut tnew, mut tcod, mut tev) =
            (vec![], vec![], vec![], vec![], vec![]);
        for _ in 0..reps {
            let d = ((ell - 1) / 2) as usize;
            let t = Instant::now();
            let ms = mg::multiples_p(f, k24, kp, d);
            let pre = mg::kernel_pre(f, &ms);
            std::hint::black_box(mg::velu_codomain_p(f, k24, &pre, ell));
            tv.push(t.elapsed().as_nanos() as f64);
            let t = Instant::now();
            std::hint::black_box(mg::isog_xz_pre(f, &pre, (u, f.one())));
            tve.push(t.elapsed().as_nanos() as f64);
            let t = Instant::now();
            let sv = SqrtVeluMont::new(f, a, kp, ell).unwrap();
            tnew.push(t.elapsed().as_nanos() as f64);
            let t = Instant::now();
            std::hint::black_box(sv.codomain_proj(f));
            tcod.push(t.elapsed().as_nanos() as f64);
            std::hint::black_box(sv.eval_nd(f, u));
            let t = Instant::now();
            std::hint::black_box(sv.eval_nd(f, u));
            tev.push(t.elapsed().as_nanos() as f64);
        }
        let (v, n, c) = (med(tv) / 1e3, med(tnew) / 1e3, med(tcod) / 1e3);
        let (ve, ev) = (med(tve) / 1e3, med(tev) / 1e3);
        println!(
            "{bits} l={ell:6}: codomain velu {v:8.0} us, sqrt {:7.0} (new {n:6.0} + cod {c:6.0}) {:5.2}x | image velu {ve:6.0} sqrt {ev:6.0} {:5.2}x",
            n + c,
            v / (n + c),
            ve / ev
        );
    }
}

fn main() {
    let ells = [3u64, 101, 401, 1009, 4001, 10007];
    run::<4>(256, &ells);
    run::<8>(511, &ells);
}
