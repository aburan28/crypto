//! Phase profile of SEA point counting over the 127-bit prime 2^127 - 1: time to compute Phi_l
//! per l, then per l the time of each phase SEA runs (specialise Phi_l(j, Y), Y^q mod it, roots,
//! Elkies codomain, BMSS kernel, eigenvalue, Atkin degree), and the whole `sea_opts` call.
//!
//! `cargo run --release --example profile_sea [max_ell] [curves]`
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::find::{bmss, elkies, sea};
use isogeny_algos::fpm::FpM;
use isogeny_algos::poly;
use std::collections::HashMap;
use std::time::Instant;

fn ms(t: Instant) -> f64 {
    t.elapsed().as_secs_f64() * 1e3
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let max_ell: usize = args.first().and_then(|s| s.parse().ok()).unwrap_or(89);
    let curves: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(2);
    let f = FpM::<2>::from_dec("170141183460469231731687303715884105727");
    let ells: Vec<usize> = (3..=max_ell).filter(|&l| is_prime(l as u64)).collect();
    let t = Instant::now();
    let mut phis: HashMap<usize, Phi<_>> = ells
        .iter()
        .copied()
        .zip(Phi::compute_many(&f, &ells))
        .collect();
    println!("phi total (compute_many, l <= {max_ell}) {:.1} ms", ms(t));
    let mut rng = Rng::new(19_500);
    for c in 0..curves {
        let e = loop {
            let e = Curve::new(f.random(&mut rng), f.random(&mut rng));
            let j = jinv(&f, &e);
            if is_smooth(&f, &e) && j != f.zero() && j != f.from_u64(1728) {
                break e;
            }
        };
        let j = jinv(&f, &e);
        println!("# curve {c}: per-l phases, ms (spec powq roots codom bmss eig atkin)");
        let mut tot = [0.0f64; 7];
        for l in 3..=max_ell {
            if !is_prime(l as u64) {
                continue;
            }
            let phi = &phis[&l];
            let t = Instant::now();
            let g = poly::monic(&f, &phi.y_poly(&f, j));
            let t_spec = ms(t);
            let t = Instant::now();
            let xi = poly::powmod_big(&f, &poly::x_poly(&f), &f.q(), &g);
            let t_pow = ms(t);
            let t = Instant::now();
            let mut roots = vec![];
            poly::split_roots(
                &f,
                &poly::gcd(&f, &g, &poly::sub(&f, &xi, &poly::x_poly(&f))),
                &mut rng,
                &mut roots,
            );
            let t_roots = ms(t);
            let (mut t_cod, mut t_bmss, mut t_eig, mut t_atk) = (0.0, 0.0, 0.0, 0.0);
            let kind;
            if let Some(&jt) = roots.first() {
                kind = "E";
                let t = Instant::now();
                let et = elkies::elkies_codomain(&f, phi, &e, jt);
                t_cod = ms(t);
                if let Some(et) = et {
                    let t = Instant::now();
                    let iso = bmss::isogeny(&f, bmss::Method::FastElkiesPrime, &e, &et, l, None);
                    t_bmss = ms(t);
                    if let Some(iso) = iso {
                        let t = Instant::now();
                        let _ = sea::elkies_eigenvalue(&f, &e, &iso.ker, l as u64);
                        t_eig = ms(t);
                    }
                }
            } else {
                kind = "A";
                let t = Instant::now();
                let _ = sea::atkin_degree_with(&f, &g, &xi);
                t_atk = ms(t);
            }
            let row = [t_spec, t_pow, t_roots, t_cod, t_bmss, t_eig, t_atk];
            for (a, b) in tot.iter_mut().zip(row) {
                *a += b;
            }
            println!(
                "l={l:3} {kind} {:8.2} {:8.2} {:8.2} {:8.2} {:8.2} {:8.2} {:8.2}",
                row[0], row[1], row[2], row[3], row[4], row[5], row[6]
            );
        }
        println!(
            "sum   {:8.1} {:8.1} {:8.1} {:8.1} {:8.1} {:8.1} {:8.1}",
            tot[0], tot[1], tot[2], tot[3], tot[4], tot[5], tot[6]
        );
        let t = Instant::now();
        let r = sea::sea_opts(&f, &e, max_ell, &mut phis, &mut rng, 0);
        let st = r.as_ref().map(|x| x.1.clone()).unwrap_or_default();
        println!(
            "sea_opts {:.1} ms ok={} elkies={} atkin={} candidates={}",
            ms(t),
            r.is_some(),
            st.elkies.len(),
            st.atkin.len(),
            st.candidates
        );
    }
}
