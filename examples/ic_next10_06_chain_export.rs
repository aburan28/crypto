//! IC next-10 item 6: export one unsatisfiable draw of the chained `m`-summand
//! Weil-descent system (the system `dreg_ladder --m M` measures) as an msolve
//! input over GF(2), with the field equations `x_i^2 + x_i` appended.
//!
//! The draw is selected exactly as `dreg_ladder --unsat-index K` selects it
//! (same cell seed, same sampling, same exact solution count), so draw `K` of a
//! cell is the same system the committed Macaulay runs measured.
//!
//! ```sh
//! cargo run --release --example ic_next10_06_chain_export -- \
//!     --m 4 --n 6 --ell 2 --seed 20261006 --unsat-index 0 --out sys.ms
//! ```
//!
//! Prints one JSON line (draw index, basis, x_R, unknowns, equations,
//! solution count) on stdout.  `--sat-ok` exports the K-th draw regardless of
//! satisfiability (for satisfiable controls), still reporting its exact count;
//! `--no-count` skips the count.

use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_bench::{
    chained_s3_solution_count, chained_solution_count, ladder_sample,
};
use crypto_lib::cryptanalysis::koblitz_groebner::{build_decomposition_system, FieldStructure};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::io::Write;

fn cell_seed(seed: u64, n: u32, ell: usize, m: usize) -> u64 {
    // Identical to dreg_ladder's cell_seed.
    let base = seed ^ (u64::from(n) << 40) ^ ((ell as u64) << 32);
    if m == 3 {
        base
    } else {
        base ^ ((m as u64) << 48)
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let flag = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let num = |name: &str, d: u64| flag(name).and_then(|v| v.parse().ok()).unwrap_or(d);
    let m = num("--m", 4) as usize;
    let n = num("--n", 6) as u32;
    let ell = num("--ell", 2) as usize;
    let seed = num("--seed", 20261006);
    let k = num("--unsat-index", 0) as usize;
    let max_draws = num("--max-draws", 4096) as usize;
    let sat_ok = args.iter().any(|a| a == "--sat-ok");
    // `--no-count`: skip the exact count (cells where it is too dear); the
    // JSON then reports `"solutions":-1` and the GB itself decides.
    let no_count = args.iter().any(|a| a == "--no-count");
    let out = flag("--out").expect("--out FILE");

    let irr = find_irreducible_sparse(n).expect("field");
    let st = FieldStructure::new(n, &irr);
    let b = F2mElement::one(n);
    let mut rng = StdRng::seed_from_u64(cell_seed(seed, n, ell, m));
    let mut seen = 0usize;
    for drawn in 0..max_draws {
        let (basis, x_r) = ladder_sample(n, ell, &mut rng);
        let solutions = if no_count {
            0
        } else if m == 3 {
            chained_s3_solution_count(&basis, &x_r, &b, &irr)
        } else {
            chained_solution_count(m, &basis, &x_r, &b, &irr)
        };
        if solutions > 0 && !sat_ok {
            continue;
        }
        if seen < k {
            seen += 1;
            continue;
        }
        let sys = build_decomposition_system(&basis, &x_r, &b, m, &st).expect("system builds");
        let nv = sys.n_vars;
        let mut f = std::fs::File::create(&out).expect("create out");
        let vars: Vec<String> = (0..nv).map(|i| format!("x{i}")).collect();
        writeln!(f, "{}", vars.join(",")).unwrap();
        writeln!(f, "2").unwrap();
        let mut polys: Vec<String> = Vec::new();
        for eq in &sys.equations {
            if eq.terms.is_empty() {
                continue;
            }
            let terms: Vec<String> = eq
                .terms
                .iter()
                .map(|t| {
                    if t.mask == 0 {
                        "1".to_string()
                    } else {
                        (0..nv)
                            .filter(|i| (t.mask >> i) & 1 == 1)
                            .map(|i| format!("x{i}"))
                            .collect::<Vec<_>>()
                            .join("*")
                    }
                })
                .collect();
            polys.push(terms.join("+"));
        }
        let n_eqs = polys.len();
        for i in 0..nv {
            polys.push(format!("x{i}^2+x{i}"));
        }
        writeln!(f, "{}", polys.join(",\n")).unwrap();
        let word = |e: &F2mElement| e.raw_bits().first().copied().unwrap_or(0);
        let max_deg = sys
            .equations
            .iter()
            .flat_map(|e| e.terms.iter().map(|t| t.mask.count_ones()))
            .max()
            .unwrap_or(0);
        println!(
            r#"{{"m":{m},"n":{n},"ell":{ell},"seed":{seed},"unsat_index":{k},"draw":{drawn},"n_vars":{nv},"n_eqs":{n_eqs},"max_eq_degree":{max_deg},"v_basis":{:?},"x_r":{},"solutions":{},"sat_ok":{sat_ok}}}"#,
            basis.iter().map(word).collect::<Vec<_>>(),
            word(&x_r),
            if no_count { -1i64 } else { solutions as i64 }
        );
        return;
    }
    eprintln!("no unsatisfiable draw {k} within {max_draws} draws");
    std::process::exit(2);
}
