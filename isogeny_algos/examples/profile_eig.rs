//! Split of the Elkies eigenvalue phase at 2^127 - 1: X^q mod h, f^((q-1)/2) mod h, and the
//! scan over lambda (whole `elkies_eigenvalue` minus the two powers).
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::find::{bmss, elkies, sea};
use isogeny_algos::fpm::FpM;
use isogeny_algos::poly;
use std::time::Instant;

fn main() {
    let f = FpM::<2>::from_dec("170141183460469231731687303715884105727");
    let ells = [29usize, 47, 53, 59, 67, 73, 79, 83];
    let phis = Phi::compute_many(&f, &ells);
    let mut rng = Rng::new(5);
    println!("l  d  xq_ms  yq_ms  scan_ms  lambda");
    let mut done = 0;
    while done < 12 {
        let e = Curve::new(f.random(&mut rng), f.random(&mut rng));
        let j = jinv(&f, &e);
        if !is_smooth(&f, &e) || j == f.zero() {
            continue;
        }
        for (phi, &l) in phis.iter().zip(&ells) {
            let g = poly::monic(&f, &phi.y_poly(&f, j));
            let xi = poly::powmod_big(&f, &poly::x_poly(&f), &f.q(), &g);
            let mut roots = vec![];
            poly::split_roots(
                &f,
                &poly::gcd(&f, &g, &poly::sub(&f, &xi, &poly::x_poly(&f))),
                &mut rng,
                &mut roots,
            );
            let Some(&jt) = roots.first() else { continue };
            let Some(et) = elkies::elkies_codomain(&f, phi, &e, jt) else {
                continue;
            };
            let Some(iso) = bmss::isogeny(&f, bmss::Method::FastElkiesPrime, &e, &et, l, None)
            else {
                continue;
            };
            let h = &iso.ker;
            let q = f.q();
            let t = Instant::now();
            let _ = poly::powmod_big(&f, &poly::x_poly(&f), &q, h);
            let txq = t.elapsed().as_secs_f64() * 1e3;
            let fxh = poly::rem(&f, &vec![e.b, e.a, f.zero(), f.one()], h);
            let t = Instant::now();
            let _ = poly::powmod_big(&f, &fxh, &q.sub_small(1).shr(1), h);
            let tyq = t.elapsed().as_secs_f64() * 1e3;
            let t = Instant::now();
            let lam = sea::elkies_eigenvalue(&f, &e, h, l as u64);
            let tall = t.elapsed().as_secs_f64() * 1e3;
            println!(
                "{l:3} {:3} {txq:6.2} {tyq:6.2} {:7.2}  {lam:?}",
                (l - 1) / 2,
                tall - txq - tyq
            );
            done += 1;
        }
    }
}
