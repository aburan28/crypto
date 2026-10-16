//! Montgomery sqrt-Velu: pair-by-pair products vs remainder trees at large l (scratch profiling).
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
    let mut rng = Rng::new(9);
    let w = isogeny_algos::testdata::supersingular_workload::<N>(ells, bits, &mut rng);
    let f: &FpM<N> = &w.f;
    let a = f.zero();
    let k24 = mg::proj24(f, a);
    for &ell in ells {
        if ell < 1000 {
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
        let reps = if ell > 30000 { 5 } else { 9 };
        let mut res = vec![];
        for tree in [false, true] {
            let (mut tc, mut te) = (vec![], vec![]);
            for _ in 0..reps {
                let mut sv = SqrtVeluMont::new(f, a, kp, ell).unwrap();
                sv.force_tree = Some(tree);
                let t = Instant::now();
                std::hint::black_box(sv.codomain_proj(f));
                tc.push(t.elapsed().as_nanos() as f64);
                std::hint::black_box(sv.eval_nd(f, u));
                let t = Instant::now();
                std::hint::black_box(sv.eval_nd(f, u));
                te.push(t.elapsed().as_nanos() as f64);
            }
            res.push((med(tc) / 1e3, med(te) / 1e3));
        }
        let (b, bp) = isogeny_algos::kernel::sqrt_velu_mont::sizes(ell);
        println!(
            "{bits} l={ell:7} b={b} b'={bp}: codomain direct {:8.0} tree {:8.0} us ({:.2}) | image direct {:8.0} tree {:8.0} us ({:.2})",
            res[0].0, res[1].0, res[0].0 / res[1].0, res[0].1, res[1].1, res[0].1 / res[1].1
        );
    }
}

fn main() {
    let args: Vec<u64> = std::env::args()
        .skip(1)
        .map(|s| s.parse().unwrap())
        .collect();
    let ells: Vec<u64> = if args.is_empty() {
        vec![3, 1009, 10007, 40009, 100003]
    } else {
        std::iter::once(3).chain(args).collect()
    };
    run::<4>(256, &ells);
    run::<8>(511, &ells);
}
