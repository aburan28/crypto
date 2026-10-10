//! How much of the Macaulay matrix work in the binary Semaev decomposition
//! systems is reduction to zero?
//!
//! The Joux–Vitse F4-trace idea (replay the useful row multiples recorded
//! on one target for every later target of the same shape) can only save
//! the rows that reduce to zero, so its ceiling is `rows / rank` at each
//! degree. This probe measures that ratio, per degree, for the ledger's
//! landed decomposition shape: Koblitz `K_0` over `F_{2^31}`, a
//! 16-dimensional subspace factor base, `m = 2` and `m = 3` summands,
//! on several random targets, and reports whether the ratio is stable
//! across targets.
//!
//! ```bash
//! cargo run --release --example macaulay_syzygy_probe -- [--n 31] [--dim 16] [--targets 4]
//! ```

use std::time::Instant;

use crypto_lib::binary_ecc::f2m::F2mElement;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, macaulay_profile, FieldStructure,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{find_irreducible, KoblitzCurve};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn elem_from_bits(bits: u64, n: u32) -> F2mElement {
    let positions: Vec<u32> = (0..n).filter(|k| (bits >> k) & 1 == 1).collect();
    F2mElement::from_bit_positions(&positions, n)
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut n = 31u32;
    let mut dim = 16u32;
    let mut targets = 4usize;
    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "--n" => {
                i += 1;
                n = argv[i].parse().unwrap_or(31);
            }
            "--dim" => {
                i += 1;
                dim = argv[i].parse().unwrap_or(16);
            }
            "--targets" => {
                i += 1;
                targets = argv[i].parse().unwrap_or(4);
            }
            _ => {}
        }
        i += 1;
    }
    let irr = find_irreducible(n).expect("irreducible");
    let k0 = KoblitzCurve::new(0, n);
    let b = k0.as_ref().map(|c| c.curve.b.clone()).unwrap_or_else(|| F2mElement::one(n));
    let st = FieldStructure::new(n, &irr);
    let basis: Vec<F2mElement> = (0..dim).map(|k| F2mElement::from_bit_positions(&[k], n)).collect();
    let mut rng = StdRng::seed_from_u64(31);
    println!("# Macaulay rows vs rank, K_0 over F_2^{n}, subspace dim {dim}\n");
    println!("| m | vars | eqs | target | degree | rows | cols | rank | rows/rank | syzygies | seconds |");
    println!("|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|");
    for m in [2usize, 3] {
        for t in 0..targets {
            let xr = elem_from_bits(rng.gen_range(1..(1u64 << n)), n);
            let Some(sys) = build_decomposition_system(&basis, &xr, &b, m, &st) else {
                println!("| {m} | — | — | {t} | — | system too large | | | | | |");
                continue;
            };
            let max_deg = if m == 2 { 4 } else { 3 };
            for d in 2..=max_deg {
                let t0 = Instant::now();
                match macaulay_profile(&sys.equations, sys.n_vars, d) {
                    Some(p) => println!(
                        "| {m} | {} | {} | {t} | {d} | {} | {} | {} | {:.3} | {} | {:.2} |",
                        sys.n_vars,
                        sys.equations.len(),
                        p.rows,
                        p.cols,
                        p.rank,
                        p.rows as f64 / p.rank.max(1) as f64,
                        p.syzygies(),
                        t0.elapsed().as_secs_f64()
                    ),
                    None => {
                        println!("| {m} | {} | {} | {t} | {d} | exceeds size limit | | | | | |", sys.n_vars, sys.equations.len());
                        break;
                    }
                }
            }
        }
    }
}
