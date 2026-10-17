//! Our side of the per-family external comparison (`external/families.sh`, against
//! `external/pari_families.gp`): kernel -> codomain, kernel -> rational maps, and l^n chains,
//! on the same inputs with the same outputs, one JSON line per case.
//!
//! Inputs: the supersingular curve E: y^2 = x^3 + x over p = 4 l^n k - 1 (p = 3 mod 4, k the
//! least integer making p prime with the stated number of bits), so E(F_p) is cyclic of order
//! p + 1 and has a rational point of order l^n: P = [(p + 1) / l^n] P0, P0 = (x0, y0) with the
//! least x0 >= 2 that works. `gen` prints the cases as GP code (with the point, so both tools
//! use the same y); the other commands rebuild them.
//!
//! Outputs: the codomain (a4, a6) of the Vélu-normalised isogeny (PARI's `ellisogeny`
//! normalisation, equal to ours), and for `maps` the x-map at r = 123456789.
//!
//! ```text
//! cargo run --release --example ext_families -- gen > cases.gp
//! cargo run --release --example ext_families -- codomain   (x-only Vélu and sqrt-Vélu)
//! cargo run --release --example ext_families -- maps       (kernel polynomial + Kohel)
//! cargo run --release --example ext_families -- chain      (l^n, optimal strategy)
//! ```
use isogeny_algos::api::{with_field, FieldTask, PrimeFieldInt};
use isogeny_algos::curve::*;
use isogeny_algos::int::Int;
use isogeny_algos::kernel::chain::{ell_power_isogeny, Strategy};
use isogeny_algos::kernel::kohel::kohel;
use isogeny_algos::kernel::sqrt_velu::sqrt_velu_fast;
use isogeny_algos::kernel::xonly::{multiples_x_fast, velu_xonly_fast};
use isogeny_algos::poly;
use std::time::Instant;

/// (bits, l, n): codomain and maps cases have n = 1
const KERNEL_CASES: [(usize, u64); 12] = [
    (61, 101),
    (61, 1009),
    (61, 10007),
    (127, 101),
    (127, 1009),
    (127, 10007),
    (255, 101),
    (255, 1009),
    (255, 10007),
    (511, 101),
    (511, 1009),
    (511, 10007),
];
const CHAIN_CASES: [(usize, u64, u32); 6] = [
    (127, 2, 60),
    (127, 3, 38),
    (255, 2, 120),
    (255, 3, 76),
    (511, 2, 240),
    (511, 3, 152),
];
const MAPS_MAX_ELL: u64 = 1009;

struct Case {
    bits: usize,
    ell: u64,
    n: u32,
    p: Int,
    x: Int,
    y: Int,
}

fn case(bits: usize, ell: u64, n: u32) -> Case {
    let ln = Int::from(ell).pow(n);
    let four_ln = &Int::from(4u64) * &ln;
    let mut k = Int::from(2u64).pow(bits as u32 - 1).div_floor(&four_ln);
    let p = loop {
        k = &k + &Int::one();
        let p = &(&four_ln * &k) - &Int::one();
        if p.bits() == bits && p.is_probable_prime() {
            break p;
        }
        assert!(p.bits() <= bits, "no prime of {bits} bits");
    };
    struct T<'a>(&'a Int, &'a Int, u64);
    impl FieldTask for T<'_> {
        type Out = (Int, Int);
        fn run<F: PrimeFieldInt>(self, f: &F) -> (Int, Int) {
            let e = Curve::new(f.one(), f.zero());
            let cof = (&(self.0 + &Int::one())).div_floor(self.1);
            let mut x0 = 2u64;
            loop {
                let x = f.from_u64(x0);
                if let Some(y) = f.sqrt(rhs(f, &e, x)) {
                    let q = pmul_big(f, &e, &Pt::Aff(x, y), cof.mag());
                    let low = pmul_big(f, &e, &q, self.1.div_floor(&Int::from(self.2)).mag());
                    if let (Pt::Aff(qx, qy), false) = (q, low == Pt::Inf) {
                        return (f.int(qx), f.int(qy));
                    }
                }
                x0 += 1;
            }
        }
    }
    let (x, y) = with_field(&p, T(&p, &ln, ell));
    Case {
        bits,
        ell,
        n,
        p,
        x,
        y,
    }
}

fn kernel_cases() -> Vec<Case> {
    KERNEL_CASES.iter().map(|&(b, l)| case(b, l, 1)).collect()
}
fn chain_cases() -> Vec<Case> {
    CHAIN_CASES.iter().map(|&(b, l, n)| case(b, l, n)).collect()
}

/// (best, median) microseconds of `reps` calls, and the last result.
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

fn line(op: &str, c: &Case, reps: usize, t: (f64, f64), out: &str) {
    println!(
        "{{\"tool\":\"ours\",\"op\":\"{op}\",\"bits\":{},\"ell\":{},\"n\":{},\"reps\":{reps},\"best_us\":{:.2},\"median_us\":{:.2},{out}}}",
        c.bits, c.ell, c.n, t.0, t.1
    );
}

fn reps(ell: u64, bits: usize) -> usize {
    if ell <= 101 {
        51
    } else if ell <= 1009 && bits <= 255 {
        11
    } else {
        5
    }
}

struct Run<'a>(&'a str, &'a Case);
impl FieldTask for Run<'_> {
    type Out = ();
    fn run<F: PrimeFieldInt>(self, f: &F) {
        let Run(op, c) = self;
        let e = Curve::new(f.one(), f.zero());
        let (x, y) = (f.elem(&c.x), f.elem(&c.y));
        let cod =
            |cv: &Curve<F::E>| format!("\"a4\":\"{}\",\"a6\":\"{}\"", f.int(cv.a), f.int(cv.b));
        match op {
            "codomain" => {
                let r = reps(c.ell, c.bits);
                let t = time(r, || velu_xonly_fast(f, &e, x, c.ell).unwrap().cod);
                line("codomain_xonly", c, r, (t.0, t.1), &cod(&t.2));
                let t = time(r, || sqrt_velu_fast(f, &e, x, c.ell).unwrap().cod);
                line("codomain_sqrtvelu", c, r, (t.0, t.1), &cod(&t.2));
            }
            "maps" => {
                if c.ell > MAPS_MAX_ELL {
                    return;
                }
                let r = reps(c.ell, c.bits);
                let t = time(r, || {
                    let n = ((c.ell - 1) / 2) as usize;
                    let xs = multiples_x_fast(f, &e, x, n).unwrap();
                    kohel(f, &e, &poly::from_roots(f, &xs), c.ell)
                });
                let at = f.from_u64(123_456_789);
                let v = f.div(poly::eval(f, &t.2.num, at), poly::eval(f, &t.2.den, at));
                let out = format!("{},\"x_map_at_r\":\"{}\"", cod(&t.2.cod), f.int(v));
                line("maps_kohel", c, r, (t.0, t.1), &out);
            }
            "chain" => {
                let l = c.ell as f64;
                let strat = Strategy::Optimal {
                    mul: 12.0 * l.log2().ceil(),
                    eval: 4.0 * l,
                    iso: 6.0 * l,
                };
                let r = 5;
                let t = time(r, || {
                    ell_power_isogeny(f, &e, &Pt::Aff(x, y), c.ell, c.n as usize, &[], strat)
                        .0
                        .cod
                });
                line("chain_optimal", c, r, (t.0, t.1), &cod(&t.2));
            }
            _ => unreachable!(),
        }
    }
}

fn main() {
    let op = std::env::args().nth(1).unwrap_or_default();
    match op.as_str() {
        "gen" => {
            for (name, cs) in [("KCASES", kernel_cases()), ("CCASES", chain_cases())] {
                let rows: Vec<String> = cs
                    .iter()
                    .map(|c| {
                        format!(
                            "[{}, {}, {}, {}, {}, {}]",
                            c.bits, c.ell, c.n, c.p, c.x, c.y
                        )
                    })
                    .collect();
                println!("{name} = [{}];", rows.join(", "));
            }
        }
        "codomain" | "maps" => {
            for c in kernel_cases() {
                with_field(&c.p, Run(&op, &c));
            }
        }
        "chain" => {
            for c in chain_cases() {
                with_field(&c.p, Run("chain", &c));
            }
        }
        _ => eprintln!("usage: ext_families gen|codomain|maps|chain"),
    }
}
