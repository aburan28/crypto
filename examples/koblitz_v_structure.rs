//! **Why is the monomial subspace cheap?**  Equation structure and solver
//! work of the m=2 decomposition system for several factor-base subspaces V.
//!
//! Families: `mono` = <1..z^{l-1}>; `shiftK` = z^K·mono; `frobJ` = mono^(2^J)
//! (Frobenius image, same curve-independent structure); `rebasis` = mono
//! with a random GL(l) change of basis (same V, different variables);
//! `randS` = random V (seed S).  Same random x(R) per probe across families.
//!
//! Run: `cargo run --release --example koblitz_v_structure <n> <l> <probes> fam,...`
//! Prints one JSON line per family.

use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    solve_boolean_system_filtered, FieldStructure, SolveOptions, SolverEngine, SplitRule,
};
use crypto_lib::cryptanalysis::koblitz_isogeny_cost::{
    factor_base_basis, field_for, geometric_basis,
};
use crypto_lib::cryptanalysis::polynomial_reuse::build_decomposition_system_reusing;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::{BTreeMap, HashSet};
use std::time::Instant;

fn el(v: u64, n: u32) -> F2mElement {
    let b: Vec<u32> = (0..n).filter(|i| (v >> i) & 1 == 1).collect();
    F2mElement::from_bit_positions(&b, n)
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let n: u32 = a[0].parse().unwrap();
    let l: u32 = a[1].parse().unwrap();
    let probes: usize = a[2].parse().unwrap();
    let solve = std::env::var("VS_NOSOLVE").is_err();
    let irr = field_for(n).unwrap();
    let st = FieldStructure::new(n, &irr);
    let mono = factor_base_basis(n, l, None);
    let mut rng = StdRng::seed_from_u64(0x5eed ^ n as u64);
    let xs: Vec<u64> = (0..probes)
        .map(|_| rng.gen::<u64>() & ((1 << n) - 1))
        .collect();
    let b = el(1, n); // a6 = 1 representative; structure does not depend on it
    for fam in a[3].split(',') {
        let basis: Vec<F2mElement> = if fam == "mono" {
            mono.clone()
        } else if let Some(k) = fam.strip_prefix("shift") {
            let zk = el(1, n).mul(&el(2, n), &irr); // z
            let mut s = el(1, n);
            for _ in 0..k.parse::<u32>().unwrap() {
                s = s.mul(&zk, &irr);
            }
            mono.iter().map(|e| e.mul(&s, &irr)).collect()
        } else if let Some(j) = fam.strip_prefix("frob") {
            mono.iter()
                .map(|e| e.square_k_times(j.parse().unwrap(), &irr))
                .collect()
        } else if fam == "rebasis" {
            // random invertible combination of the monomial basis
            loop {
                let rows: Vec<u64> = (0..l).map(|_| rng.gen::<u64>() & ((1 << l) - 1)).collect();
                let mut ech: Vec<u64> = vec![];
                let mut ok = true;
                for &r in &rows {
                    let mut x = r;
                    for &e in &ech {
                        x = x.min(x ^ e);
                    }
                    if x == 0 {
                        ok = false;
                        break;
                    }
                    ech.push(x);
                    ech.sort_unstable_by(|a, b| b.cmp(a));
                }
                if ok {
                    break rows.iter().map(|&r| el(r, n)).collect();
                }
            }
        } else if let Some(s) = fam.strip_prefix("geom") {
            geometric_basis(n, l, s.parse().unwrap(), &irr)
        } else if let Some(s) = fam.strip_prefix("rand") {
            factor_base_basis(n, l, Some(s.parse().unwrap()))
        } else {
            panic!("family {fam}")
        };
        let opts = SolveOptions {
            engine: SolverEngine::MatrixF4 { max_degree: 4 },
            max_solutions: usize::MAX,
            node_budget: 20_000,
            split_rule: SplitRule::default(),
        };
        let (mut terms, mut eqs, mut mono_union, mut red, mut ns, mut sols, mut exh) =
            (0usize, 0usize, 0usize, 0usize, 0u128, 0usize, 0usize);
        let mut deg_hist: BTreeMap<u32, usize> = BTreeMap::new();
        let mut maxdeg_built = 0u32;
        let mut qrank = 0usize;
        for &x in &xs {
            let sys = build_decomposition_system_reusing(&basis, &el(x, n), &b, 2, &st).unwrap();
            let mut u = HashSet::new();
            for e in &sys.equations {
                eqs += 1;
                terms += e.terms.len();
                for t in &e.terms {
                    u.insert(t.mask);
                    *deg_hist.entry(t.degree()).or_default() += 1;
                }
            }
            mono_union += u.len();
            // rank of the quadratic parts: n_eqs - rank = linear equations
            // obtainable by linear algebra alone (degree falls at degree 2)
            let quads: Vec<u64> = {
                let mut q: Vec<u64> = u.iter().copied().filter(|m| m.count_ones() >= 2).collect();
                q.sort_unstable();
                q
            };
            let rows: Vec<Vec<bool>> = sys
                .equations
                .iter()
                .map(|e| {
                    let mut r = vec![false; quads.len()];
                    for t in &e.terms {
                        if let Ok(i) = quads.binary_search(&t.mask) {
                            r[i] ^= true;
                        }
                    }
                    r
                })
                .collect();
            qrank += rank(rows);
            if !solve {
                continue;
            }
            let t0 = Instant::now();
            let (r, s) =
                solve_boolean_system_filtered(&sys.equations, sys.n_vars, &opts, |_| false);
            ns += t0.elapsed().as_nanos();
            sols += r.len();
            red += s.reductions;
            exh += s.exhausted as usize;
            maxdeg_built = maxdeg_built.max(s.max_degree_built);
        }
        let p = probes as f64;
        println!(
            "{{\"n\":{n},\"l\":{l},\"family\":\"{fam}\",\"probes\":{probes},\"eqs_per_sys\":{:.1},\"terms_per_eq\":{:.2},\"distinct_monos_per_sys\":{:.1},\"deg_hist\":{:?},\"reductions_per_call\":{:.1},\"ms_per_call\":{:.2},\"exhausted\":{exh},\"max_degree_built\":{maxdeg_built},\"roots\":{sols},\"quad_rank\":{:.2}}}",
            eqs as f64 / p, terms as f64 / eqs as f64, mono_union as f64 / p, deg_hist,
            red as f64 / p, ns as f64 / p / 1e6, qrank as f64 / p
        );
    }
}

fn rank(mut rows: Vec<Vec<bool>>) -> usize {
    let cols = rows.first().map_or(0, |r| r.len());
    let mut r = 0;
    for c in 0..cols {
        let Some(p) = (r..rows.len()).find(|&i| rows[i][c]) else {
            continue;
        };
        rows.swap(r, p);
        for i in 0..rows.len() {
            if i != r && rows[i][c] {
                let (a, b) = if i < r {
                    let (x, y) = rows.split_at_mut(r);
                    (&mut x[i], &y[0])
                } else {
                    let (x, y) = rows.split_at_mut(i);
                    (&mut y[0], &x[r])
                };
                for k in 0..cols {
                    a[k] ^= b[k];
                }
            }
        }
        r += 1;
    }
    r
}
