//! Interleaved harness. Sides:
//!   a = crypto_lib::pqc::fast::ml_kem (baseline, main)
//!   b = mlkem_dev::ml_kem             (candidate, unprepared API)
//!   c = mlkem_dev::ml_kem prepared    (candidate, public-data reuse)
//! Deterministic inputs (rotating seeds), so every side does the same work and
//! callgrind counts are stable.
//!
//!   mbench cycles [rounds]
//!   mbench count <a|b|c> <512|768|1024> <keygen|encaps|decaps> <iters>
use crypto_lib::pqc::fast::ml_kem as base;
use crypto_lib::pqc::ml_kem::{
    MlKemDecapsKey, MlKemEncapsKey, MlKemParams, ML_KEM_1024, ML_KEM_512, ML_KEM_768,
};
use mlkem_dev::ml_kem as cand;
use std::hint::black_box;

fn seeds(n: usize, tag: u8) -> Vec<[u8; 32]> {
    let mut x = 0x9e37_79b9_7f4a_7c15u64 ^ ((tag as u64) << 40);
    (0..n)
        .map(|_| {
            let mut s = [0u8; 32];
            for b in s.iter_mut() {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                *b = (x >> 24) as u8;
            }
            s
        })
        .collect()
}

#[inline]
fn rdtsc() -> u64 {
    unsafe { core::arch::x86_64::_rdtsc() }
}

fn params(s: &str) -> MlKemParams {
    match s {
        "512" => ML_KEM_512,
        "768" => ML_KEM_768,
        "1024" => ML_KEM_1024,
        _ => panic!("set"),
    }
}

const NS: usize = 16;

struct Fix {
    p: MlKemParams,
    d: Vec<[u8; 32]>,
    z: Vec<[u8; 32]>,
    m: Vec<[u8; 32]>,
    ek: Vec<MlKemEncapsKey>,
    dk: Vec<MlKemDecapsKey>,
    ct: Vec<Vec<u8>>,
    pe: Vec<cand::MlKemPreparedEncapsKey>,
    pd: Vec<cand::MlKemPreparedDecapsKey>,
}

fn fixture(p: MlKemParams) -> Fix {
    let d = seeds(NS, 1);
    let z = seeds(NS, 2);
    let m = seeds(NS, 3);
    let (mut ek, mut dk, mut ct, mut pe, mut pd) = (vec![], vec![], vec![], vec![], vec![]);
    for i in 0..NS {
        let (e, k) = base::ml_kem_keygen_internal(&p, &d[i], &z[i]);
        let (c, ss) = base::ml_kem_encaps_internal(&p, &e, &m[i]).unwrap();
        let a = cand::MlKemPreparedEncapsKey::new(&p, &e).unwrap();
        let b = cand::MlKemPreparedDecapsKey::new(&p, &k).unwrap();
        // Every side must produce the same bytes before any of it is timed.
        assert_eq!(cand::ml_kem_encaps_internal(&p, &e, &m[i]).unwrap(), (c.clone(), ss));
        assert_eq!(a.encaps_internal(&m[i]), (c.clone(), ss));
        assert_eq!(b.decaps(&c).unwrap(), ss);
        assert_eq!(cand::ml_kem_decaps(&p, &k, &c).unwrap(), ss);
        let (e2, k2) = cand::ml_kem_keygen_internal(&p, &d[i], &z[i]);
        assert!(e2 == e && k2 == k);
        ek.push(e);
        dk.push(k);
        ct.push(c);
        pe.push(a);
        pd.push(b);
    }
    Fix { p, d, z, m, ek, dk, ct, pe, pd }
}

fn run(side: char, op: &str, f: &Fix, i: usize) {
    let p = &f.p;
    let i = i % NS;
    match (side, op) {
        ('a', "keygen") => {
            black_box(base::ml_kem_keygen_internal(p, &f.d[i], &f.z[i]));
        }
        ('b', "keygen") => {
            black_box(cand::ml_kem_keygen_internal(p, &f.d[i], &f.z[i]));
        }
        ('c', "keygen") => {
            black_box(cand::ml_kem_keygen_prepared_internal(p, &f.d[i], &f.z[i]));
        }
        ('a', "encaps") => {
            black_box(base::ml_kem_encaps_internal(p, &f.ek[i], &f.m[i]));
        }
        ('b', "encaps") => {
            black_box(cand::ml_kem_encaps_internal(p, &f.ek[i], &f.m[i]));
        }
        ('c', "encaps") => {
            black_box(f.pe[i].encaps_internal(&f.m[i]));
        }
        ('a', "decaps") => {
            black_box(base::ml_kem_decaps(p, &f.dk[i], &f.ct[i]));
        }
        ('b', "decaps") => {
            black_box(cand::ml_kem_decaps(p, &f.dk[i], &f.ct[i]));
        }
        ('c', "decaps") => {
            black_box(f.pd[i].decaps(&f.ct[i]));
        }
        ('b', "prep-ek") => {
            black_box(cand::MlKemPreparedEncapsKey::new(p, &f.ek[i]));
        }
        ('b', "prep-dk") => {
            black_box(cand::MlKemPreparedDecapsKey::new(p, &f.dk[i]));
        }
        _ => panic!("op"),
    }
}

fn median(v: &mut [u64]) -> f64 {
    v.sort_unstable();
    v[v.len() / 2] as f64 / 1e3
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(|s| s.as_str()) == Some("count") {
        let f = fixture(params(&args[3]));
        let side = args[2].chars().next().unwrap();
        let iters: usize = args[5].parse().unwrap();
        for i in 0..iters {
            run(side, &args[4], &f, i);
        }
        return;
    }
    let rounds: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(5);
    let reps = 200;
    println!("kilocycles (rdtsc), median of {} samples per cell, ABC interleaved", rounds * reps);
    println!("| set | op | base | cand | prepared | base/cand | base/prepared |");
    println!("|---|---|---:|---:|---:|---:|---:|");
    for set in ["512", "768", "1024"] {
        let f = fixture(params(set));
        for op in ["keygen", "encaps", "decaps"] {
            let mut s: [Vec<u64>; 3] = [vec![], vec![], vec![]];
            for _ in 0..20 {
                for side in ['a', 'b', 'c'] {
                    run(side, op, &f, 0);
                }
            }
            for _ in 0..rounds {
                for (k, side) in ['a', 'b', 'c'].into_iter().enumerate() {
                    for r in 0..reps {
                        let t0 = rdtsc();
                        run(side, op, &f, r);
                        let t1 = rdtsc();
                        s[k].push(t1 - t0);
                    }
                }
            }
            let a = median(&mut s[0]);
            let b = median(&mut s[1]);
            let c = median(&mut s[2]);
            println!(
                "| {set} | {op} | {a:.1} | {b:.1} | {c:.1} | {:.2}x | {:.2}x |",
                a / b,
                a / c
            );
        }
        for op in ["prep-ek", "prep-dk"] {
            let mut v = vec![];
            for r in 0..reps * rounds {
                let t0 = rdtsc();
                run('b', op, &f, r);
                v.push(rdtsc() - t0);
            }
            println!("| {set} | {op} | | {:.1} | | | |", median(&mut v));
        }
    }
}
