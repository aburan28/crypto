//! dudect-style timing test of the secret phase of decapsulation (Toy16).
//!
//! One key.  Class 0 = one fixed ciphertext, class 1 = a pool of fresh
//! ciphertexts (half of them FO-rejected: made for another key).  Classes are
//! interleaved in a random order; the public preparation runs outside the
//! timed region.  Welch's t statistic is reported on the raw samples and on
//! samples cropped at several percentiles (Reparaz–Balasch–Verbauwhede).
//! |t| > 4.5 is the usual "leak detected" threshold.
//!
//!   cargo run --release --example dudect_decaps -- [samples_per_class]
use eigen_seta::kem::{EigenSeta, SecretKey, DecapsPrep, Ciphertext};
use eigen_seta::params::{Toy16, Toy16Fp};
use rand_core::{OsRng, RngCore};
use std::time::Instant;

type K = EigenSeta<Toy16, 4>;
type Prep = DecapsPrep<eigen_seta::fp2::Fp2<Toy16Fp, 4>>;

#[inline(never)]
fn timed(sk: &SecretKey, prep: &Prep, ct: &Ciphertext) -> u64 {
    let t = Instant::now();
    let k = K::decaps_core(std::hint::black_box(sk), std::hint::black_box(prep), std::hint::black_box(ct));
    std::hint::black_box(k);
    t.elapsed().as_nanos() as u64
}

fn welch(a: &[f64], b: &[f64]) -> f64 {
    let m = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let v = |x: &[f64], mu: f64| x.iter().map(|y| (y - mu) * (y - mu)).sum::<f64>() / (x.len() as f64 - 1.0);
    let (ma, mb) = (m(a), m(b));
    (ma - mb) / (v(a, ma) / a.len() as f64 + v(b, mb) / b.len() as f64).sqrt()
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(300);
    let mut rng = OsRng;
    let (pk, sk) = K::keygen(&mut rng);
    let (pk2, _) = K::keygen(&mut rng);
    let pool_n = 16;
    let fixed = K::encaps(&pk, &mut rng).unwrap().0;
    let fixed_prep = K::decaps_prepare(&sk, &fixed).unwrap();
    let pool: Vec<(Ciphertext, Prep)> = (0..pool_n)
        .map(|i| {
            let ct = if i % 2 == 0 { K::encaps(&pk, &mut rng) } else { K::encaps(&pk2, &mut rng) }.unwrap().0;
            let prep = K::decaps_prepare(&sk, &ct).unwrap();
            (ct, prep)
        })
        .collect();
    // warm-up
    for _ in 0..5 {
        timed(&sk, &fixed_prep, &fixed);
    }
    let mut s0 = Vec::with_capacity(n);
    let mut s1 = Vec::with_capacity(n);
    let mut j = 0;
    while s0.len() < n || s1.len() < n {
        let class = (rng.next_u32() & 1) as usize;
        if class == 0 && s0.len() < n {
            s0.push(timed(&sk, &fixed_prep, &fixed) as f64);
        } else if class == 1 && s1.len() < n {
            let (ct, prep) = &pool[j % pool_n];
            j += 1;
            s1.push(timed(&sk, prep, ct) as f64);
        }
    }
    let med = |x: &[f64]| { let mut y = x.to_vec(); y.sort_by(|a, b| a.partial_cmp(b).unwrap()); y[y.len() / 2] };
    println!("samples/class {n}; median ns: fixed {:.0}, random {:.0}", med(&s0), med(&s1));
    println!("t(raw) = {:+.2}", welch(&s0, &s1));
    let mut all: Vec<f64> = s0.iter().chain(s1.iter()).copied().collect();
    all.sort_by(|a, b| a.partial_cmp(b).unwrap());
    for pct in [0.5, 0.75, 0.9] {
        let thr = all[((all.len() as f64) * pct) as usize];
        let a: Vec<f64> = s0.iter().copied().filter(|x| *x <= thr).collect();
        let b: Vec<f64> = s1.iter().copied().filter(|x| *x <= thr).collect();
        println!("t(crop <= p{:.0}) = {:+.2}  (n = {}, {})", pct * 100.0, welch(&a, &b), a.len(), b.len());
    }
}
