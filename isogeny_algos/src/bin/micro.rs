//! Micro-benchmarks of the arithmetic primitives (ns per operation). Usage: micro [--out FILE]
use isogeny_algos::field::*;
use isogeny_algos::poly;
use isogeny_algos::series;
use std::hint::black_box;
use std::io::Write;
use std::time::Instant;

fn per_op<T>(reps: usize, mut f: impl FnMut() -> T) -> f64 {
    // warm-up
    for _ in 0..reps / 10 + 1 {
        black_box(f());
    }
    let mut best = f64::INFINITY;
    for _ in 0..5 {
        let t = Instant::now();
        for _ in 0..reps {
            black_box(f());
        }
        best = best.min(t.elapsed().as_nanos() as f64 / reps as f64);
    }
    best
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out = args
        .iter()
        .position(|a| a == "--out")
        .map(|i| args[i + 1].clone());
    let mut lines = vec![];
    let mut rec = |name: &str, size: usize, ns: f64| {
        let s = format!("{{\"op\":\"{name}\",\"size\":{size},\"ns\":{ns:.2}}}");
        println!("{s}");
        lines.push(s);
    };
    let mut rng = Rng::new(1);
    for bits in [30u32, 61] {
        let fp = Zp::new(next_prime(1u64 << (bits - 1)));
        let xs: Vec<u64> = (0..1024).map(|_| fp.random(&mut rng)).collect();
        let mut i = 0usize;
        rec(
            &format!("zp{bits}_mul"),
            1,
            per_op(200_000, || {
                i = (i + 1) & 1023;
                fp.mul(xs[i], xs[(i + 7) & 1023])
            }),
        );
        rec(
            &format!("zp{bits}_inv"),
            1,
            per_op(20_000, || {
                i = (i + 1) & 1023;
                fp.inv(xs[i] | 1)
            }),
        );
        rec(
            &format!("zp{bits}_sqrt"),
            1,
            per_op(2_000, || {
                i = (i + 1) & 1023;
                Field::sqrt(&fp, xs[i])
            }),
        );
        let fm = isogeny_algos::fpm::FpM::<1>::from_u64_modulus(fp.p);
        let ms: Vec<[u64; 1]> = (0..1024).map(|_| fm.random(&mut rng)).collect();
        rec(
            &format!("fpm1_{bits}_mul"),
            1,
            per_op(200_000, || {
                i = (i + 1) & 1023;
                fm.mul(ms[i], ms[(i + 7) & 1023])
            }),
        );
        let mut chain = ms[0];
        rec(
            &format!("fpm1_{bits}_mul_latency"),
            1,
            per_op(200_000, || {
                chain = fm.mul(chain, ms[5]);
                chain
            }),
        );
        let mut chainz = xs[0];
        rec(
            &format!("zp{bits}_mul_latency"),
            1,
            per_op(200_000, || {
                chainz = fp.mul(chainz, xs[5]);
                chainz
            }),
        );
        rec(
            &format!("fpm1_{bits}_inv"),
            1,
            per_op(20_000, || {
                i = (i + 1) & 1023;
                fm.inv(ms[i])
            }),
        );
        let f2 = Zp2::new(fp.p);
        let ys: Vec<(u64, u64)> = (0..1024).map(|_| f2.random(&mut rng)).collect();
        rec(
            &format!("zp2_{bits}_mul"),
            1,
            per_op(200_000, || {
                i = (i + 1) & 1023;
                f2.mul(ys[i], ys[(i + 7) & 1023])
            }),
        );
        rec(
            &format!("zp2_{bits}_inv"),
            1,
            per_op(20_000, || {
                i = (i + 1) & 1023;
                f2.inv(ys[i])
            }),
        );
        for n in [16usize, 64, 256, 1024] {
            let a: Vec<u64> = (0..n).map(|_| fp.random(&mut rng)).collect();
            let b: Vec<u64> = (0..n).map(|_| fp.random(&mut rng)).collect();
            let reps = (4_000_000 / (n * n)).max(3);
            rec(
                &format!("poly{bits}_mul"),
                n,
                per_op(reps, || poly::mul(&fp, &a, &b)),
            );
            let c: Vec<u64> = (0..2 * n).map(|_| fp.random(&mut rng)).collect();
            let mut m: Vec<u64> = (0..n).map(|_| fp.random(&mut rng)).collect();
            m.push(1);
            rec(
                &format!("poly{bits}_rem_2n_by_n"),
                n,
                per_op(reps, || poly::rem(&fp, &c, &m)),
            );
            if n <= 256 {
                let x = poly::x_poly(&fp);
                rec(
                    &format!("poly{bits}_powmod_x_p"),
                    n,
                    per_op((reps / 64).max(2), || {
                        poly::powmod(&fp, &x, fp.p as u128, &m)
                    }),
                );
                rec(
                    &format!("poly{bits}_gcd"),
                    n,
                    per_op(reps.max(2), || poly::gcd(&fp, &a, &m)),
                );
            }
            let mut s = a.clone();
            s[0] = 1;
            rec(
                &format!("series{bits}_inv"),
                n,
                per_op(reps, || series::inv(&fp, &s, n)),
            );
            rec(
                &format!("series{bits}_mul"),
                n,
                per_op(reps, || series::mul(&fp, &a, &b, n)),
            );
        }
    }
    if let Some(o) = out {
        let mut f = std::fs::File::create(o).unwrap();
        for l in lines {
            writeln!(f, "{l}").unwrap();
        }
    }
}
