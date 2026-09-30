//! Brute-force validation of the F4 verdict on random small Boolean systems.

use crate::boolpoly::{Mono, Poly, W};
use crate::f4::{linear_solutions, truncated_groebner, Outcome};

pub fn run(seed: u64, trials: u32) {
    let only: Option<u32> = std::env::var("TRIAL").ok().map(|v| v.parse().unwrap());
    let mut s = seed | 1;
    let mut next = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    let mut counts = [0u32; 4];
    for t in 0..trials {
        let nv = 4 + (next() % 5) as usize; // 4..8
        let ne = nv - 1 + (next() % 4) as usize;
        let mut eqs = Vec::new();
        for _ in 0..ne {
            let mut terms = Vec::new();
            for i in 0..nv {
                if next() % 3 == 0 {
                    terms.push(Mono::var(i));
                }
                for j in i + 1..nv {
                    if next() % 4 == 0 {
                        terms.push(Mono::var(i).mul(&Mono::var(j)));
                    }
                }
            }
            if next() % 2 == 0 {
                terms.push(Mono::ONE);
            }
            eqs.push(Poly::from_terms(terms));
        }
        // brute force
        let mut sols: Vec<u64> = Vec::new();
        for a in 0..(1u64 << nv) {
            let pt = [a, 0, 0, 0, 0, 0, 0, 0];
            if eqs.iter().all(|e| !e.eval(&pt)) {
                sols.push(a);
            }
        }
        if let Some(o) = only { if o != t { continue; } }
        let affine = is_affine(&sols);
        let (out, st) = truncated_groebner(&eqs, nv as u32, u64::MAX, u64::MAX);
        let verdict = match out {
            Outcome::Inconsistent => {
                counts[0] += 1;
                assert!(sols.is_empty(), "trial {t}: F4 says inconsistent but {} solutions", sols.len());
                "inconsistent"
            }
            Outcome::Linear(rref) => {
                counts[1] += 1;
                let ls = linear_solutions(&rref, nv, 20).unwrap();
                let mut got: Vec<u64> = ls.iter().map(|p| p[0]).collect();
                got.sort();
                if got != sols {
                    for e in &eqs { eprintln!("eq: {e}"); }
                    for r in &rref { eprintln!("rref: {r}"); }
                    panic!("trial {t}: linear solution set differs: got {got:?} want {sols:?} (nv={nv})");
                }
                "linear"
            }
            Outcome::Wild => {
                counts[2] += 1;
                assert!(!affine, "trial {t}: complete GB wild but solutions {:?} affine (nv={nv}, ne={ne})", sols);
                "wild"
            }
            Outcome::Budget => {
                counts[3] += 1;
                "budget"
            }
        };
        let _ = (verdict, st);
        let _ = W;
    }
    eprintln!("selftest ok: inconsistent={} linear={} wild={} budget={}", counts[0], counts[1], counts[2], counts[3]);
}

fn is_affine(sols: &[u64]) -> bool {
    if sols.is_empty() {
        return true;
    }
    let base = sols[0];
    let set: std::collections::HashSet<u64> = sols.iter().copied().collect();
    for &a in sols {
        for &b in sols {
            if !set.contains(&(a ^ b ^ base)) {
                return false;
            }
        }
    }
    true
}
