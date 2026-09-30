//! Where one degree of one ladder draw's refutation spends its work.
//!
//! ```sh
//! cargo run --release --example refute_profile -- \
//!     --n 8 --basis 232,164,245,177 --xr 200 --degree 6 [--f5-count]
//! ```
//!
//! `--basis` and `--xr` are the `v_basis` and `x_r` fields of a ladder
//! JSON line.  Prints the Macaulay shape, the build time, and per column
//! degree band the pivots, row additions, words written and the nonzeros
//! left, so a speedup can be aimed at the band that costs.
//!
//! A profiling tool: it reports no degree, and its timings are
//! practicality notes on whatever machine it runs.

use crypto_lib::cryptanalysis::koblitz_bench::{ladder_f5_row_counts, ladder_refutation_profile};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let flag = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let n: u32 = flag("--n").expect("--n").parse().expect("n");
    let basis: Vec<u64> = flag("--basis")
        .expect("--basis a,b,..")
        .split(',')
        .map(|w| w.trim().parse().expect("basis word"))
        .collect();
    let x_r: u64 = flag("--xr").expect("--xr").parse().expect("x_r");
    let degree: u32 = flag("--degree").expect("--degree").parse().expect("degree");

    // `--f5-count`: only count the rows, and the rows the F5 criterion keeps.
    if args.iter().any(|a| a == "--f5-count") {
        let (full, kept, ms) = ladder_f5_row_counts(n, &basis, x_r, degree).expect("draw builds");
        println!(
            "n={n} ell={} degree={degree} rows {full} f5_kept {kept} pruned {} ({:.1}%) in {ms:.0} ms",
            basis.len(),
            full - kept,
            100.0 * (full - kept) as f64 / full.max(1) as f64
        );
        return;
    }

    let p = ladder_refutation_profile(n, &basis, x_r, degree).expect("draw builds");
    println!(
        "n={n} ell={} degree={degree} vars={} matrix {} x {} nnz {} build {:.0} ms",
        basis.len(),
        p.n_vars,
        p.rows,
        p.cols,
        p.start_nnz,
        p.build_ms
    );
    println!("| band degree | columns | pivots | row adds | words written | nnz left | ms |");
    println!("|--:|--:|--:|--:|--:|--:|--:|");
    for b in &p.bands {
        println!(
            "| {} | {} | {} | {} | {} | {} | {:.0} |",
            b.degree, b.columns, b.pivots, b.xors, b.xor_output, b.active_nnz_after, b.ms
        );
    }
    println!(
        "high rank {} vanished {} linear rows {} max weight {}",
        p.high_rank, p.vanished, p.linear_rows, p.max_weight
    );
}
