//! Measure what each modelling choice costs the SAT pipeline for
//! binary-Semaev index calculus.
//!
//! Three axes, measured independently:
//!
//! 1. **Factor base** — confining the unknowns to an `l`-dimensional
//!    subspace instead of all of `F_{2ⁿ}`.
//! 2. **Parity encoding** — native XOR rows (Gauss-Jordan inside the
//!    solver) versus Tseitin expansion into CNF.
//! 3. **Decomposition size** — `S₃` (two points) versus symmetrised
//!    `S₄` (three points), which is where index calculus starts to pay.
//!
//! ```bash
//! cargo run --release --example semaev_sat_bench
//! ```

use crypto_lib::binary_ecc::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_semaev::binary_semaev_s3;
use crypto_lib::cryptanalysis::binary_semaev_s4::{elementary_symmetric_3, symmetrised_s4_eval};
use crypto_lib::cryptanalysis::sat::SolveResult;
use crypto_lib::cryptanalysis::semaev_sat::{
    decode_x1_x2, encode_semaev_s3_subspace, encode_semaev_s4, XorEncoding,
};
use std::time::Instant;

fn irr_for(n: u32) -> IrreduciblePoly {
    let low = match n {
        4 => vec![0, 1],
        5 => vec![0, 2],
        6 => vec![0, 1],
        7 => vec![0, 1],
        8 => vec![0, 2, 3, 4],
        9 => vec![0, 1],
        10 => vec![0, 3],
        11 => vec![0, 2],
        12 => vec![0, 1, 2, 3],
        19 => vec![0, 1, 2, 5],
        _ => panic!("no irreducible tabulated for n = {n}"),
    };
    IrreduciblePoly {
        degree: n,
        low_terms: low,
    }
}

/// Build a `b` that puts the symmetric point `X₁ = X₂ = x₃` on the
/// `S₃` variety, so the instance is known-SAT.
fn planted_s3(n: u32, irr: &IrreduciblePoly, x3: &F2mElement) -> F2mElement {
    let x3_2 = x3.square(irr);
    let x3_3 = x3_2.mul(x3, irr);
    let x3_4 = x3_3.mul(x3, irr);
    x3_3.add(&x3_4)
}

fn run_s3(n: u32, l: u32, encoding: XorEncoding, budget: u64) {
    let irr = irr_for(n);
    let x3 = F2mElement::from_bit_positions(&[1, 3], n);
    let b = planted_s3(n, &irr, &x3);

    let t0 = Instant::now();
    let mut enc = encode_semaev_s3_subspace(n, l, &irr, &b, &x3, encoding);
    let build = t0.elapsed();
    let (vars, clauses, xors) = (
        enc.solver.n_vars(),
        enc.solver.n_clauses(),
        enc.solver.n_xors(),
    );

    enc.solver.conflict_budget = budget;
    let t1 = Instant::now();
    let res = enc.solver.solve();
    let solve = t1.elapsed();

    let verdict = match res {
        SolveResult::Sat => {
            let (x1, x2) = decode_x1_x2(&enc);
            if binary_semaev_s3(&x1, &x2, &x3, &b, &irr).is_zero() {
                "SAT (verified)"
            } else {
                "SAT (BAD DECODE)"
            }
        }
        SolveResult::Unsat => "UNSAT",
        SolveResult::Unknown => "budget hit",
    };

    println!(
        "  n={n:<3} l={l:<3} {:<7} {vars:>6} {clauses:>8} {xors:>6}  {:>8.1?} {:>10.1?}  {verdict}",
        match encoding {
            XorEncoding::Native => "native",
            XorEncoding::Cnf => "cnf",
        },
        build,
        solve
    );
}

fn run_s4(n: u32, l: u32, x_r: &F2mElement, encoding: XorEncoding, budget: u64, solve_it: bool) {
    let irr = irr_for(n);
    let b = F2mElement::one(n);

    let t0 = Instant::now();
    let mut enc = encode_semaev_s4(n, l, &irr, &b, x_r, encoding);
    let build = t0.elapsed();
    let (vars, clauses, xors) = (
        enc.solver.n_vars(),
        enc.solver.n_clauses(),
        enc.solver.n_xors(),
    );

    if !solve_it {
        println!(
            "  n={n:<3} l={l:<3} {:<7} {vars:>6} {clauses:>8} {xors:>6}  {:>8.1?} {:>10}  (not solved)",
            match encoding {
                XorEncoding::Native => "native",
                XorEncoding::Cnf => "cnf",
            },
            build,
            "-"
        );
        return;
    }

    enc.solver.conflict_budget = budget;
    let t1 = Instant::now();
    let res = enc.solver.solve();
    let solve = t1.elapsed();

    let verdict = match res {
        SolveResult::Sat => {
            let xs = enc.decode();
            let (e1, e2, e3) = elementary_symmetric_3(&xs[0], &xs[1], &xs[2], &irr);
            if symmetrised_s4_eval(&e1, &e2, &e3, x_r, &irr).is_zero() {
                "SAT (verified)".to_string()
            } else {
                "SAT (BAD DECODE)".to_string()
            }
        }
        SolveResult::Unsat => "UNSAT".to_string(),
        SolveResult::Unknown => "budget hit".to_string(),
    };

    println!(
        "  n={n:<3} l={l:<3} {:<7} {vars:>6} {clauses:>8} {xors:>6}  {:>8.1?} {:>10.1?}  {verdict}",
        match encoding {
            XorEncoding::Native => "native",
            XorEncoding::Cnf => "cnf",
        },
        build,
        solve
    );
}

fn header() {
    println!(
        "  {:<11} {:<7} {:>6} {:>8} {:>6}  {:>8} {:>10}  {}",
        "params", "parity", "vars", "clauses", "xors", "build", "solve", "result"
    );
}

fn main() {
    println!("\n=== S₃, unrestricted (l = n): n equations in 2n unknowns ===");
    header();
    for n in [4u32, 5, 6, 7] {
        run_s3(n, n, XorEncoding::Cnf, 2_000_000);
        run_s3(n, n, XorEncoding::Native, 2_000_000);
    }

    println!("\n=== S₃, factor base l = ⌈n/2⌉: n equations in ~n unknowns ===");
    header();
    for n in [6u32, 8, 10, 12] {
        let l = n.div_ceil(2);
        run_s3(n, l, XorEncoding::Cnf, 2_000_000);
        run_s3(n, l, XorEncoding::Native, 2_000_000);
    }

    println!("\n=== symmetrised S₄, factor base l = 6, n = 19 (corpus instance) ===");
    header();
    let n = 19;
    let x_r = F2mElement::from_bit_positions(&[1, 2, 5, 6, 7, 8, 10, 11, 12, 16], n);
    run_s4(n, 6, &x_r, XorEncoding::Cnf, 0, false);
    run_s4(n, 6, &x_r, XorEncoding::Native, 8_000_000, true);
    println!();
}
