//! Our side of the external comparison (`external/`): the same workloads, inputs and output
//! digests as `external/flint_bench.c` and `external/pari_bench.gp`, one JSON line per case.
//!
//! Inputs (shared by all three programs), full-size residues: with G = 11400714819323198485
//! (the 64-bit golden-ratio constant), a_i = G^(i + 16), b_i = G^(i + 1000003) mod p, and the
//! monic modulus g of degree d has g_i = G^(i + 2000003) (i < d). Digest of a
//! polynomial: its value at r = 123456789 mod p, in decimal. Curves for SEA:
//! y^2 = x^3 + (k + 3) x + (7 k + 11), k = 0, 1, ..
//!
//! ```text
//! cargo run --release --example ext_compare -- poly  <p61|p127|p255|p511>
//! cargo run --release --example ext_compare -- powx  <p61|p127|p255|p511>
//! cargo run --release --example ext_compare -- phi   p127 <max_ell>
//! cargo run --release --example ext_compare -- sea   <p61|p127> <curves>
//! ```
use isogeny_algos::bigint::Big;
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::find::sea;
use isogeny_algos::fpm::FpM;
use isogeny_algos::poly::{self, Poly};
use std::collections::HashMap;
use std::time::Instant;

const P61: &str = "2305843009213693951";
const P127: &str = "170141183460469231731687303715884105727";
const P255: &str = "57896044618658097711785492504343953926634992332820282019728792003956564819949";
const P511: &str = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";
const LENS: [usize; 11] = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384];
const DEGS: [usize; 7] = [32, 64, 90, 128, 256, 512, 1024];

const G: u64 = 11_400_714_819_323_198_485;

/// G^(i + s) mod p for i < n.
fn seq<const N: usize>(f: &FpM<N>, n: usize, s: u64) -> Vec<[u64; N]> {
    let g = f.from_u64(G);
    let mut x = f.pow(g, s as u128);
    (0..n)
        .map(|_| {
            let v = x;
            x = f.mul(x, g);
            v
        })
        .collect()
}

fn digest<const N: usize>(f: &FpM<N>, a: &Poly<FpM<N>>) -> String {
    f.to_big(&poly::eval(f, a, f.from_u64(123_456_789)))
        .to_dec()
}

/// (best, median) of `reps` timed calls, in microseconds, and the last result.
fn time<T>(reps: usize, mut g: impl FnMut() -> T) -> (f64, f64, T) {
    let mut ts = vec![];
    let mut out = None;
    for _ in 0..reps {
        let t = Instant::now();
        let r = g();
        ts.push(t.elapsed().as_secs_f64() * 1e6);
        out = Some(r);
    }
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (ts[0], ts[ts.len() / 2], out.unwrap())
}

fn reps_for(n: usize) -> usize {
    if n <= 256 {
        50
    } else if n <= 2048 {
        9
    } else {
        5
    }
}

fn poly_case<const N: usize>(tag: &str, f: &FpM<N>) {
    for &n in &LENS {
        if N >= 8 && n > 8192 {
            continue;
        }
        let a = seq(f, n, 16);
        let b = seq(f, n, 1_000_003);
        let reps = reps_for(n);
        let (best, med, r) = time(reps, || poly::mul(f, &a, &b));
        println!(
            "{{\"tool\":\"isogeny_algos\",\"op\":\"mul\",\"p\":\"{tag}\",\"n\":{n},\"reps\":{reps},\"best_us\":{best:.2},\"median_us\":{med:.2},\"digest\":\"{}\"}}",
            digest(f, &r)
        );
        let (best, med, r) = time(reps, || poly::mul(f, &a, &a));
        println!(
            "{{\"tool\":\"isogeny_algos\",\"op\":\"sqr\",\"p\":\"{tag}\",\"n\":{n},\"reps\":{reps},\"best_us\":{best:.2},\"median_us\":{med:.2},\"digest\":\"{}\"}}",
            digest(f, &r)
        );
    }
}

fn powx_case<const N: usize>(tag: &str, f: &FpM<N>) {
    let x = poly::x_poly(f);
    for &d in &DEGS {
        let mut g = seq(f, d, 2_000_003);
        g.push(f.one());
        let reps = if d <= 128 { 9 } else { 3 };
        let (best, med, r) = time(reps, || poly::powmod_big(f, &x, &f.q(), &g));
        println!(
            "{{\"tool\":\"isogeny_algos\",\"op\":\"powx\",\"p\":\"{tag}\",\"n\":{d},\"reps\":{reps},\"best_us\":{best:.2},\"median_us\":{med:.2},\"digest\":\"{}\"}}",
            digest(f, &r)
        );
    }
}

fn primes_upto(m: usize) -> Vec<usize> {
    (3..=m).filter(|&l| is_prime(l as u64)).collect()
}

fn phi_case<const N: usize>(tag: &str, f: &FpM<N>, max_ell: usize) {
    let ells = primes_upto(max_ell);
    let t = Instant::now();
    let phis = Phi::compute_many(f, &ells);
    let total = t.elapsed().as_secs_f64() * 1e6;
    let (j0, r) = (f.from_u64(987_654_321), f.from_u64(123_456_789));
    for (l, phi) in ells.iter().zip(&phis) {
        println!(
            "{{\"tool\":\"isogeny_algos\",\"op\":\"phi_digest\",\"p\":\"{tag}\",\"ell\":{l},\"digest\":\"{}\"}}",
            f.to_big(&phi.eval(f, j0, r)).to_dec()
        );
    }
    println!(
        "{{\"tool\":\"isogeny_algos\",\"op\":\"phi_all\",\"p\":\"{tag}\",\"max_ell\":{max_ell},\"total_us\":{total:.0}}}"
    );
}

fn sea_case<const N: usize>(tag: &str, f: &FpM<N>, max_ell: usize, curves: usize) {
    let ells = primes_upto(max_ell);
    let t = Instant::now();
    let mut phis: HashMap<usize, Phi<FpM<N>>> = ells
        .iter()
        .copied()
        .zip(Phi::compute_many(f, &ells))
        .collect();
    let pre = t.elapsed().as_secs_f64() * 1e6;
    println!(
        "{{\"tool\":\"isogeny_algos\",\"op\":\"sea_phi_precompute\",\"p\":\"{tag}\",\"max_ell\":{max_ell},\"total_us\":{pre:.0}}}"
    );
    let mut rng = Rng::new(42);
    for k in 0..curves as u64 {
        let e = Curve::new(f.from_u64(k + 3), f.from_u64(7 * k + 11));
        if !is_smooth(f, &e) {
            continue;
        }
        let t = Instant::now();
        let r = sea::sea_opts(f, &e, max_ell, &mut phis, &mut rng, 0);
        let us = t.elapsed().as_secs_f64() * 1e6;
        let order = r.map(|x| x.0.to_dec()).unwrap_or_else(|| "none".into());
        println!(
            "{{\"tool\":\"isogeny_algos\",\"op\":\"sea\",\"p\":\"{tag}\",\"k\":{k},\"us\":{us:.0},\"order\":\"{order}\"}}"
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (op, tag) = (args[0].as_str(), args[1].as_str());
    let arg = |i: usize, d: usize| args.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
    macro_rules! with_field {
        ($run:ident $(, $extra:expr)*) => {
            match tag {
                "p61" => $run(tag, &FpM::<1>::from_dec(P61) $(, $extra)*),
                "p127" => $run(tag, &FpM::<2>::from_dec(P127) $(, $extra)*),
                "p255" => $run(tag, &FpM::<4>::from_dec(P255) $(, $extra)*),
                "p511" => $run(tag, &FpM::<8>::from_dec(P511) $(, $extra)*),
                _ => panic!("unknown prime {tag}"),
            }
        };
    }
    match op {
        "poly" => with_field!(poly_case),
        "powx" => with_field!(powx_case),
        "phi" => with_field!(phi_case, arg(2, 89)),
        "sea" => {
            let max_ell = if tag == "p61" { 61 } else { 89 };
            with_field!(sea_case, max_ell, arg(2, 10))
        }
        _ => panic!("unknown op {op}"),
    }
    let _ = Big::from_u64(0);
}
