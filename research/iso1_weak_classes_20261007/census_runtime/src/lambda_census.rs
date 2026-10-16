//! Exhaustive Legendre census over F_{p^6} for small p (added 2026-10-16).
//!
//! For every lambda outside {0, 1}, records the certified trace of
//! L_lambda: y^2 = x(x-1)(x-lambda), the 2-Sylow type (a, b) of L_lambda(F_Q)
//! with a <= b, whether lambda lies in Omega = (F_Q^*)^(2^e) with
//! e = v2(Q-1), and whether it lies in the norm-one torus T of order q^2+q+1.
//! Output rows aggregate counts by (trace, a, b, in_omega, in_t).
use crypto_lib::certified::certified_order;
use crypto_lib::field_kernel::{EllE, Fld, Fq3, PtE6, E2, E6};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::Write;

fn sylow2(ec: &EllE, n: u128, rng: &mut StdRng) -> (u32, u32) {
    let v = n.trailing_zeros();
    let m = n >> v;
    // b = exponent of the 2-Sylow subgroup: the maximum order over 24 random
    // points (a smaller exponent survives with probability <= 2^-24 per curve).
    let mut b = 0;
    for _ in 0..24 {
        let mut r: PtE6 = ec.mul_u128(&ec.random_point(rng), m);
        let mut j = 0;
        while !r.inf {
            r = ec.add(&r, &r);
            j += 1;
        }
        b = b.max(j);
    }
    (v - b, b)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let p: u64 = args.next().expect("p").parse().unwrap();
    let out = args.next().expect("output csv");
    let f = Fq3::new(p);
    let q = p * p;
    let q3 = (q as u128).pow(3);
    // ISO1_OMEGA_EXP=k replaces e by min(k, v2(Q-1)): Omega becomes the 2^k-th powers.
    let e = match std::env::var("ISO1_OMEGA_EXP") {
        Ok(k) => k.parse::<u32>().unwrap().min((q3 - 1).trailing_zeros()),
        Err(_) => (q3 - 1).trailing_zeros(),
    };
    let e2 = |k: u64| E2([k % p, k / p]);
    let e6 = |k: u64| E6([e2(k % q), e2((k / q) % q), e2(k / (q * q))]);
    // ISO1_SKIP_SYLOW=1 records a = b = 0 and skips the 2-Sylow sampling.
    let skip_sylow = std::env::var("ISO1_SKIP_SYLOW").is_ok_and(|v| v == "1");
    let nonsquare = (2..).map(e6).find(|x| !f.is_square(x)).unwrap();
    // ISO1_OMEGA_ONLY=1 enumerates Omega = <g^(2^e)> multiplicatively instead
    // of all of F_Q; rows then all have in_omega = 1.
    if std::env::var("ISO1_OMEGA_ONLY").is_ok_and(|v| v == "1") {
        let n = q3 - 1;
        let mut primes = Vec::new();
        let (mut m, mut d) = (n, 2u128);
        while d * d <= m {
            if m % d == 0 {
                primes.push(d);
                while m % d == 0 {
                    m /= d;
                }
            }
            d += 1;
        }
        if m > 1 {
            primes.push(m);
        }
        let g = (2..)
            .map(e6)
            .find(|x| primes.iter().all(|&r| f.pow(x, n / r) != E6::ONE))
            .unwrap();
        let h = f.pow(&g, 1u128 << e);
        let size = (n >> e) as u64;
        let block = 4096u64;
        // ISO1_CHECKPOINT=FILE: each finished block appends `b|t:in_t:count,...`;
        // a restarted run reloads them and skips those blocks.
        let checkpoint = std::env::var("ISO1_CHECKPOINT").ok();
        let mut resumed = BTreeMap::<(i128, u32, u32, bool, bool), u64>::new();
        let mut done = std::collections::HashSet::<u64>::new();
        if let Some(path) = &checkpoint {
            if let Ok(text) = std::fs::read_to_string(path) {
                for line in text.lines() {
                    let Some((b, tail)) = line.split_once('|') else {
                        continue;
                    };
                    if !tail.ends_with(';') || !done.insert(b.parse().unwrap()) {
                        continue; // torn final line, or a duplicate block
                    }
                    for item in tail
                        .trim_end_matches(';')
                        .split(',')
                        .filter(|x| !x.is_empty())
                    {
                        let v: Vec<i128> = item.split(':').map(|x| x.parse().unwrap()).collect();
                        *resumed.entry((v[0], 0, 0, true, v[1] == 1)).or_default() += v[2] as u64;
                    }
                }
            }
            eprintln!("checkpoint {path}: resumed {} blocks", done.len());
        }
        let log = checkpoint.as_ref().map(|path| {
            std::sync::Mutex::new(
                std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(path)
                    .unwrap(),
            )
        });
        let starts: Vec<u64> = (0..size.div_ceil(block))
            .filter(|b| !done.contains(b))
            .collect();
        let agg = starts
            .into_par_iter()
            .map(|b| {
                let f = Fq3::new(p);
                let mut rng = StdRng::seed_from_u64(0x0e6a ^ b ^ (p << 40));
                let mut m = BTreeMap::<(i128, u32, u32, bool, bool), u64>::new();
                let mut lam = f.pow(&h, (b * block) as u128);
                for _ in b * block..((b + 1) * block).min(size) {
                    if lam != E6::ONE {
                        let a2 = f.neg(&f.add(&E6::ONE, &lam));
                        let n = certified_order(&f, a2, lam, &nonsquare, &mut rng, 32)
                            .expect("order not isolated")
                            .order;
                        let in_t = f.pow(&lam, (q * q + q + 1) as u128) == E6::ONE;
                        *m.entry((q3 as i128 + 1 - n as i128, 0, 0, true, in_t))
                            .or_default() += 1;
                    }
                    lam = f.mul(&lam, &h);
                }
                if let Some(log) = &log {
                    let items: Vec<String> = m
                        .iter()
                        .map(|((t, _, _, _, tt), c)| format!("{t}:{}:{c}", u8::from(*tt)))
                        .collect();
                    let line = format!("{b}|{};\n", items.join(","));
                    log.lock().unwrap().write_all(line.as_bytes()).unwrap();
                }
                m
            })
            .reduce(BTreeMap::new, |mut x, y| {
                for (k, v) in y {
                    *x.entry(k).or_default() += v;
                }
                x
            });
        let mut agg = agg;
        for (k, v) in resumed {
            *agg.entry(k).or_default() += v;
        }
        let mut w = std::io::BufWriter::new(std::fs::File::create(&out).unwrap());
        writeln!(w, "p,e,trace,a,b,in_omega,in_t,count").unwrap();
        for ((t, a, b, o, tt), c) in agg {
            writeln!(w, "{p},{e},{t},{a},{b},{},{},{c}", o as u8, tt as u8).unwrap();
        }
        return;
    }
    let chunks: Vec<u64> = (0..q * q).collect();
    let agg = chunks
        .into_par_iter()
        .map(|hi| {
            let f = Fq3::new(p);
            let mut rng = StdRng::seed_from_u64(0x1a4b_da ^ hi ^ (p << 40));
            let mut m = BTreeMap::<(i128, u32, u32, bool, bool), u64>::new();
            for lo in 0..q {
                let lam = e6(hi * q + lo);
                if lam == E6::ZERO || lam == E6::ONE {
                    continue;
                }
                // y^2 = x^3 - (1+lambda) x^2 + lambda x
                let a2 = f.neg(&f.add(&E6::ONE, &lam));
                let a4 = lam;
                let n = certified_order(&f, a2, a4, &nonsquare, &mut rng, 32)
                    .expect("order not isolated")
                    .order;
                let (a, b) = if skip_sylow {
                    (0, 0)
                } else {
                    sylow2(&EllE::from_a2_a4(&f, a2, a4), n, &mut rng)
                };
                let in_omega = f.pow(&lam, (q3 - 1) >> e) == E6::ONE;
                let in_t = f.pow(&lam, (q * q + q + 1) as u128) == E6::ONE;
                *m.entry((q3 as i128 + 1 - n as i128, a, b, in_omega, in_t))
                    .or_default() += 1;
            }
            m
        })
        .reduce(BTreeMap::new, |mut x, y| {
            for (k, v) in y {
                *x.entry(k).or_default() += v;
            }
            x
        });
    let mut w = std::io::BufWriter::new(std::fs::File::create(&out).unwrap());
    writeln!(w, "p,e,trace,a,b,in_omega,in_t,count").unwrap();
    for ((t, a, b, o, tt), c) in agg {
        writeln!(w, "{p},{e},{t},{a},{b},{},{},{c}", o as u8, tt as u8).unwrap();
    }
}
