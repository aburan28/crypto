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
    // multi-limb Montgomery fields: P-256 (N = 4) and the CSIDH-512 prime (N = 8)
    fn big_field<const N: usize>(
        name: &str,
        f: &isogeny_algos::fpm::FpM<N>,
        rng: &mut Rng,
        rec: &mut impl FnMut(&str, usize, f64),
    ) {
        let ms: Vec<[u64; N]> = (0..1024).map(|_| f.random(rng)).collect();
        let mut i = 0usize;
        rec(
            &format!("{name}_mul"),
            1,
            per_op(200_000, || {
                i = (i + 1) & 1023;
                f.mul(ms[i], ms[(i + 7) & 1023])
            }),
        );
        let mut chain = ms[0];
        rec(
            &format!("{name}_mul_latency"),
            1,
            per_op(200_000, || {
                chain = f.mul(chain, ms[5]);
                chain
            }),
        );
        rec(
            &format!("{name}_sq"),
            1,
            per_op(200_000, || {
                i = (i + 1) & 1023;
                f.sq(ms[i])
            }),
        );
        rec(
            &format!("{name}_add"),
            1,
            per_op(1_000_000, || {
                i = (i + 1) & 1023;
                f.add(ms[i], ms[(i + 7) & 1023])
            }),
        );
        rec(
            &format!("{name}_inv"),
            1,
            per_op(2_000, || {
                i = (i + 1) & 1023;
                f.inv(ms[i])
            }),
        );
        let e = f.modulus().sub_small(1).shr(1);
        rec(
            &format!("{name}_legendre_pow"),
            1,
            per_op(2_000, || {
                i = (i + 1) & 1023;
                f.pow_big(ms[i], &e)
            }),
        );
        rec(
            &format!("{name}_sqrt"),
            1,
            per_op(500, || {
                i = (i + 1) & 1023;
                f.sqrt(f.sq(ms[i]))
            }),
        );
    }
    // polynomial products over multi-limb fields, by Karatsuba threshold
    fn big_poly<const N: usize>(
        name: &str,
        f: &isogeny_algos::fpm::FpM<N>,
        rng: &mut Rng,
        rec: &mut impl FnMut(&str, usize, f64),
    ) {
        for n in [8usize, 16, 32, 64, 128, 256] {
            let a: Vec<[u64; N]> = (0..n).map(|_| f.random(rng)).collect();
            let b: Vec<[u64; N]> = (0..n).map(|_| f.random(rng)).collect();
            let reps = (400_000 / (n * n)).max(3);
            for th in [4usize, 8, 16, 32, 1 << 20] {
                isogeny_algos::fpm::KARATSUBA_N.with(|k| k.set(Some(th)));
                let t = per_op(reps, || poly::mul(f, &a, &b));
                rec(&format!("{name}_polymul_kth{th}"), n, t);
            }
            isogeny_algos::fpm::KARATSUBA_N.with(|k| k.set(None));
        }
    }
    let p256 = isogeny_algos::fpm::FpM::<4>::from_dec(
        "115792089210356248762697446949407573530086143415290314195533631308867097853951",
    );
    big_field("fpm4_p256", &p256, &mut rng, &mut rec);
    big_poly("fpm4_p256", &p256, &mut rng, &mut rec);
    let c512 = isogeny_algos::path::csidh::Csidh::<isogeny_algos::fpm::FpM<8>>::csidh512().fp;
    big_field("fpm8_csidh512", &c512, &mut rng, &mut rec);
    big_poly("fpm8_csidh512", &c512, &mut rng, &mut rec);
    if let Some(o) = out {
        let mut f = std::fs::File::create(o).unwrap();
        for l in lines {
            writeln!(f, "{l}").unwrap();
        }
    }
}
