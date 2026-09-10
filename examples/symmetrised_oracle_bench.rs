//! Paired decomposition-oracle benchmark: the production `x`-system against
//! the symmetrised Artin–Schreier-frame system, on the same curve and the
//! same invariant subspace, with every verdict gated by enumeration.
//!
//! ```bash
//! cargo run --release --example symmetrised_oracle_bench                 # ladder up to n = 15, plus n = 21 at m = 2
//! cargo run --release --example symmetrised_oracle_bench -- --full       # also n = 21 at m = 3 (minutes)
//! cargo run --release --example symmetrised_oracle_bench -- --targets 12 --no-sat
//! ```
//!
//! Columns: Boolean unknowns / equations / degree of the arm's system;
//! found / refuted / inconclusive targets (mixes differ between the two
//! bases, so medians are split by verdict and never compared across
//! regimes); `via T` = relations that lifted through `R + T`; median
//! milliseconds per verdict; median effort (F4 splits or SAT conflicts);
//! first fall degree; and the enumeration gate.

use crypto_lib::cryptanalysis::koblitz_symmetrised::{format_paired, paired_bench, PairedOptions};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let flag = |f: &str| args.iter().any(|a| a == f);
    let value = |f: &str| {
        args.iter()
            .position(|a| a == f)
            .and_then(|i| args.get(i + 1))
            .and_then(|v| v.parse::<usize>().ok())
    };
    let full = flag("--full");
    let mut opts = PairedOptions::default();
    if let Some(t) = value("--targets") {
        opts.targets = t;
    }
    if flag("--no-sat") {
        opts.sat = false;
    }
    if flag("--no-direct") {
        opts.direct_x = false;
    }
    if flag("--whole-group") {
        opts.targets_in_subgroup = false;
    }

    println!("=== Paired oracle benchmark: x-system vs symmetrised system ===");
    println!(
        "targets = {}, node budget = {}, SAT conflict budget = {}, targets in <G> = {}",
        opts.targets, opts.node_budget, opts.conflict_budget, opts.targets_in_subgroup
    );
    println!();

    let mut ladder: Vec<(u8, u32, usize)> = vec![
        (0, 9, 2),
        (1, 9, 2),
        (0, 9, 3),
        (1, 9, 3),
        (0, 15, 2),
        (1, 15, 2),
        (0, 15, 3),
        (1, 15, 3),
        (0, 21, 2),
        (1, 21, 2),
    ];
    if full {
        ladder.push((0, 21, 3));
        ladder.push((1, 21, 3));
    }
    for (a, n, m) in ladder {
        let t0 = Instant::now();
        // The direct x-S4 arm is Boolean degree 7; past n = 9 its Macaulay
        // matrices exceed the engine's limits and it only splits.  Keep it
        // where it can finish.
        let mut o = opts.clone();
        if n > 9 {
            o.direct_x = false;
        }
        match paired_bench(a, n, m, &o) {
            Some(b) => {
                print!("{}", format_paired(&b));
                println!(
                    "   divisor indices {:?}  ({:.1} s)",
                    b.divisor,
                    t0.elapsed().as_secs_f64()
                );
            }
            None => println!("== K_{a}/F_2^{n} m = {m}: no instance (dimension or size limit)"),
        }
        println!();
    }
}
