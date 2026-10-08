//! One KLPT run with timings of the stages (for profiling).
use isogeny_algos::field::Rng;
use isogeny_algos::int::Int;
use isogeny_algos::quat::klpt::*;
use isogeny_algos::quat::*;
use std::time::Instant;

fn main() {
    let mut rng = Rng::new(902);
    let p = Int::from_big(&isogeny_algos::bigint::Big::from_dec(
        &std::env::args().nth(1).unwrap_or("1000003".into()),
    ));
    // ideal norm: argument 2 ("small" = 3^5 * 5, else a prime near p)
    let generic = std::env::args().nth(2).as_deref() == Some("generic");
    eprintln!("p = {p}");
    let alg = Alg::new(&p);
    eprintln!("alg");
    let o0 = alg.o0();
    eprintln!("o0 {:?}", o0);
    let n = if generic {
        let mut m = &p + &Int::from(2i64);
        while !m.is_probable_prime() {
            m = &m + &Int::from(2i64);
        }
        m
    } else {
        Int::from(3i64.pow(5) * 5)
    };
    let b = o0.basis();
    let t0 = Instant::now();
    let a = if generic {
        // alpha = a + b i + c j + d k with n | Nrd(alpha): a = sqrt(-(b^2 + p(c^2 + d^2))) mod n
        loop {
            let (bb, cc, dd) = (
                Int::from(rng.next() >> 1),
                Int::from(rng.next() >> 1),
                Int::from(rng.next() >> 1),
            );
            let t = -&(&(&bb * &bb) + &(&p * &(&(&cc * &cc) + &(&dd * &dd))));
            if let Some(aa) = Int::sqrt_mod_prime(&t, &n) {
                let x = Quat::new([aa, bb, cc, dd], Int::one());
                let i = ideal_from(&alg, &o0, &n, &x);
                if i.norm(&o0) == (n.clone(), Int::one()) {
                    break x;
                }
            }
        }
    } else {
        loop {
            let mut x = Quat::zero();
            for e in &b {
                let c = rng.below(2001) as i64 - 1000;
                x = x.add(&e.scale(&Int::from(c), &Int::one()));
            }
            let (na, _) = alg.nrd(&x);
            if na.modulo(&n).is_zero() && !x.is_zero() {
                let i = ideal_from(&alg, &o0, &n, &x);
                if i.norm(&o0) == (n.clone(), Int::one()) {
                    break x;
                }
            }
        }
    };
    eprintln!("ideal found in {:?}", t0.elapsed());
    let i = ideal_from(&alg, &o0, &n, &a);
    let t0 = Instant::now();
    let rb = i.reduced_basis(&alg);
    eprintln!(
        "LLL {:?}: {:?}",
        t0.elapsed(),
        rb.iter().map(|v| alg.qf(v)).collect::<Vec<_>>()
    );
    let t0 = Instant::now();
    let sv = i.short_vectors(&alg, &(&alg.qf(&rb[0]) * &Int::from(16i64)), 2000);
    eprintln!("short vectors {:?}: {}", t0.elapsed(), sv.len());
    let t0 = Instant::now();
    let res = klpt(&alg, &i, 2, &mut rng);
    eprintln!(
        "klpt {:?}: {:?}",
        t0.elapsed(),
        res.map(|r| (r.e, r.e0, r.e1, r.n_prime))
    );
}
