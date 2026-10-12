//! E6 of `research/fall_degree_bounds_20261008/PREREGISTRATION.md`
//! (Amendment 4): `d_F` beyond the ladder's reach.
//!
//! - `newcells`: Part B. Draw cells exactly as `dreg_ladder` does
//!   (`ladder_sample` with `cell_seed(seed, n, ℓ, m)`), keep the first `k`
//!   draws with no solutions, and measure `d_F` (the first `c` with
//!   `1 ∈ V_{F,c}`) up to `cmax`.
//! - `sat`: Part C. On committed satisfiable `m = 3` draws, enumerate the
//!   solution set of the chained system `S₃(x₁, x₂, e₁)`, `S₃(e₁, x₃, x_R)`
//!   exactly, check it against every equation, and compute `d_last` as
//!   `lfd_fast.sage` defines it (closure dimensions against
//!   `dim (I(S) ∩ R_{≤c})`, threshold `DG`).
//!
//! The closure engine is `examples/hky_last_fall_degree.rs`, included as a
//! module, so this measures with the validated engine.
#[path = "hky_last_fall_degree.rs"]
#[allow(dead_code)]
mod engine;

use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_bench::{chained_s3_solution_count, chained_solution_count, ladder_sample};
use crypto_lib::cryptanalysis::koblitz_groebner::{build_decomposition_system, FieldStructure};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use engine::{closure_dims, eval, ideal_side, monomials_deglex, quad_roots, ArtinSchreier, Field};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::{json, Value};
use std::time::Instant;

/// `dreg_ladder`'s cell seed (examples/dreg_ladder.rs).
fn cell_seed(seed: u64, n: u32, ell: usize, m: usize) -> u64 {
    let base = seed ^ (u64::from(n) << 40) ^ ((ell as u64) << 32);
    if m == 3 {
        base
    } else {
        base ^ ((m as u64) << 48)
    }
}

fn word(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

fn system(m: usize, n: u32, basis: &[F2mElement], x_r: &F2mElement) -> (Vec<Vec<u64>>, u32) {
    let irr = find_irreducible_sparse(n).unwrap();
    let st = FieldStructure::new(n, &irr);
    let sys = build_decomposition_system(basis, x_r, &F2mElement::one(n), m, &st).expect("system builds");
    let eqs = sys
        .equations
        .iter()
        .map(|p| p.terms.iter().map(|t| t.mask).collect())
        .collect();
    (eqs, sys.n_vars as u32)
}

fn refutation_degree(eqs: &[Vec<u64>], nv: u32, cmax: u32) -> (Option<u32>, serde_json::Map<String, Value>) {
    let mut dims = serde_json::Map::new();
    for c in 2..=cmax {
        let full = monomials_deglex(nv, c).len();
        let (_, closed) = closure_dims(eqs, nv, c, full);
        dims.insert(c.to_string(), json!(closed));
        if closed == full {
            return (Some(c), dims);
        }
    }
    (None, dims)
}

fn field_of(n: u32) -> Field {
    let irr = find_irreducible_sparse(n).unwrap();
    let mut poly = 1u64 << n;
    for &t in &irr.low_terms {
        poly |= 1u64 << t;
    }
    Field { n, poly }
}

/// All solutions of the `m = 3` chain, as variable assignments in
/// `build_decomposition_system`'s layout: summand `i` owns bits
/// `[iℓ, (i+1)ℓ)`, then `e₁`'s `n` bits.
fn chain3_solutions(n: u32, basis: &[u64], x_r: u64) -> Vec<u64> {
    let f = field_of(n);
    let asl = ArtinSchreier::new(&f);
    let ell = basis.len();
    let elem = |c: u64| -> u64 {
        (0..ell)
            .filter(|j| c >> j & 1 == 1)
            .fold(0, |acc, j| acc ^ basis[j])
    };
    // S₃(a, c, d) = (a + c)² d² + a c d + (a c)² + 1, as a quadratic in d
    let s3_roots = |a: u64, c: u64| -> Vec<u64> {
        let ac = f.mul(a, c);
        if a ^ c == 0 && ac == 0 {
            return vec![]; // constant 1: no root
        }
        quad_roots(&f, &asl, f.sqr(a ^ c), ac, f.sqr(ac) ^ 1)
    };
    let mut out = vec![];
    for c1 in 0..(1u64 << ell) {
        for c2 in 0..(1u64 << ell) {
            for e1 in s3_roots(elem(c1), elem(c2)) {
                for c3 in 0..(1u64 << ell) {
                    let (a, c) = (e1, elem(c3));
                    let ac = f.mul(a, c);
                    let v = f.mul(f.sqr(a ^ c), f.sqr(x_r)) ^ f.mul(ac, x_r) ^ f.sqr(ac) ^ 1;
                    if v == 0 {
                        out.push(c1 | c2 << ell | c3 << (2 * ell) | e1 << (3 * ell));
                    }
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| args.iter().position(|a| a == k).map(|i| args[i + 1].clone());
    match args.get(1).map(String::as_str) {
        Some("newcells") => {
            // --cells m:n:ell,... --seed S --keep K --cmax C
            let seed: u64 = get("--seed").unwrap().parse().unwrap();
            let keep: usize = get("--keep").map_or(4, |v| v.parse().unwrap());
            let cmax: u32 = get("--cmax").map_or(4, |v| v.parse().unwrap());
            for cell in get("--cells").unwrap().split(',') {
                let p: Vec<usize> = cell.split(':').map(|x| x.parse().unwrap()).collect();
                let (m, n, ell) = (p[0], p[1] as u32, p[2]);
                let mut rng = StdRng::seed_from_u64(cell_seed(seed, n, ell, m));
                let irr = find_irreducible_sparse(n).unwrap();
                let (mut kept, mut drawn) = (0, 0);
                while kept < keep && drawn < 1000 {
                    let (basis, x_r) = ladder_sample(n, ell, &mut rng);
                    let draw = drawn;
                    drawn += 1;
                    let b1 = F2mElement::one(n);
                    let sols = if m == 3 {
                        chained_s3_solution_count(&basis, &x_r, &b1, &irr)
                    } else {
                        chained_solution_count(m, &basis, &x_r, &b1, &irr)
                    };
                    if sols > 0 {
                        continue;
                    }
                    kept += 1;
                    let t0 = Instant::now();
                    let (eqs, nv) = system(m, n, &basis, &x_r);
                    let (d_f, dims) = refutation_degree(&eqs, nv, cmax);
                    println!("{}", json!({
                        "schema": "hky_grid_newcell/v1", "m": m, "n": n, "ell": ell, "draw": draw, "seed": seed,
                        "v_basis": basis.iter().map(word).collect::<Vec<_>>(), "x_r": word(&x_r),
                        "n_vars": nv, "n_eqs": eqs.len(), "solutions": 0,
                        "d_f": d_f, "d_f_at_least": if d_f.is_none() { Some(cmax + 1) } else { None },
                        "closure_dims": dims, "secs": (t0.elapsed().as_secs_f64() * 10.0).round() / 10.0,
                    }));
                }
            }
        }
        Some("sat") => {
            let cmax: u32 = get("--cmax").map_or(5, |v| v.parse().unwrap());
            let keep: usize = get("--keep").map_or(4, |v| v.parse().unwrap());
            let text = std::fs::read_to_string(&args[2]).unwrap();
            let mut kept = 0;
            for line in text.lines().filter(|l| !l.trim().is_empty()) {
                let rec: Value = serde_json::from_str(line).unwrap();
                if rec.get("draw").is_none() || rec["solutions"].as_u64().unwrap_or(0) == 0 || kept >= keep {
                    continue;
                }
                kept += 1;
                let t0 = Instant::now();
                let n = rec["n"].as_u64().unwrap() as u32;
                let words: Vec<u64> = rec["v_basis"].as_array().unwrap().iter().map(|w| w.as_u64().unwrap()).collect();
                let x_r = rec["x_r"].as_u64().unwrap();
                let el = |w: u64| F2mElement::from_biguint(&num_bigint::BigUint::from(w), n);
                let basis: Vec<F2mElement> = words.iter().map(|&w| el(w)).collect();
                let (eqs, nv) = system(3, n, &basis, &el(x_r));
                let sols = chain3_solutions(n, &words, x_r);
                // every enumerated solution satisfies every equation
                for &s in &sols {
                    for e in &eqs {
                        assert!(!e.iter().fold(false, |acc, &m| acc ^ eval(m, s)), "solution fails an equation");
                    }
                }
                let (idim, dg) = ideal_side(nv, &sols, cmax);
                let mut clo = std::collections::HashMap::new();
                let mut e0 = None;
                let mut dims = serde_json::Map::new();
                for c in 2..=cmax {
                    let (_, closed) = closure_dims(&eqs, nv, c, idim[c as usize]);
                    clo.insert(c, closed);
                    dims.insert(c.to_string(), json!([closed, idim[c as usize]]));
                    if closed == idim[c as usize] && c >= dg {
                        e0 = Some(c);
                        break;
                    }
                }
                let d_last = e0.map(|e| {
                    let mut dl = e;
                    for c in (2..e).rev() {
                        if clo.get(&c) == Some(&idim[c as usize]) {
                            dl = c;
                        } else {
                            break;
                        }
                    }
                    dl
                });
                println!("{}", json!({
                    "schema": "hky_grid_sat/v1", "cell": rec["cell"], "n": n, "ell": rec["ell"], "draw": rec["draw"],
                    "solutions_recorded": rec["solutions"], "solutions_enumerated": sols.len(),
                    "n_vars": nv, "DG": dg, "e0": e0, "d_last": d_last,
                    "d_last_at_least": if d_last.is_none() { Some(cmax + 1) } else { None },
                    "closure_vs_ideal_dims": dims, "secs": (t0.elapsed().as_secs_f64() * 10.0).round() / 10.0,
                }));
            }
        }
        _ => eprintln!("usage: hky_grid newcells --cells m:n:ell,... --seed S [--keep K] [--cmax C] | sat LADDER.jsonl [--keep K] [--cmax C]"),
    }
}
