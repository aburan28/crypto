//! Bounded Macaulay diagnostics (legacy executable name).
//!
//! Reports rank deficiency, refutation, and pinning; these are not
//! certified first fall degree, F4/F5 solving degree, or regularity.
//! Targets are unlifted field abscissas. Both random controls have zero
//! as a solution. See research/degree_reporting_20260922/README.md.
//!
//! ```sh
//! cargo run --release --example dreg_sweep -- --n-max 7 --d-max 2 --trials 2
//! ```

use crypto_lib::cryptanalysis::koblitz_bench::{
    dreg_summary, format_dreg_table, subspace_ladder, DregSummary,
};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let flag = |name: &str, default: usize| -> usize {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    };
    let d_max = flag("--d-max", 5) as u32;
    let trials = flag("--trials", 8);
    let n_max = flag("--n-max", 21) as u32;
    let max_ell = flag("--max-ell", 12) as u32;
    let seed = flag("--seed", 0x5EED) as u64;
    let with_control = !args.iter().any(|a| a == "--no-control");
    let only_m = args
        .iter()
        .position(|a| a == "--m")
        .and_then(|i| args.get(i + 1))
        .map(|v| {
            v.split(',')
                .filter_map(|x| x.trim().parse::<usize>().ok())
                .collect::<Vec<_>>()
        })
        .unwrap_or_else(|| vec![2, 3]);

    println!("# Bounded Macaulay rank and resolution diagnostics");
    println!();
    println!("d_max = {d_max}, trials = {trials}, seed = {seed:#x}");
    println!();

    let mut rows: Vec<DregSummary> = Vec::new();
    for (n, ell) in subspace_ladder(n_max, max_ell) {
        for &m in &only_m {
            let started = std::time::Instant::now();
            if let Some(r) = dreg_summary(n, 0, m, d_max, trials, seed, with_control) {
                eprintln!(
                    "n={n} ell={ell} m={m}: vars={} refute={:?} ({:.1}s)",
                    r.n_vars,
                    r.refute_mean,
                    started.elapsed().as_secs_f64()
                );
                rows.push(r);
            }
        }
    }

    println!("{}", format_dreg_table(&rows));
    println!(
        "`D_rank_proxy` is rank deficiency, not certified first fall degree.\n\
         `D_refute` is bounded-Macaulay refutation; `delta_means` subtracts\n\
         conditional means over potentially different draws. Both controls\n\
         always have the all-zero root; their degrees measure pinning.\n\
         `unres` combines no refutation/pinning with resource caps; it is\n\
         not a solving-degree lower bound. Targets are not point-lifted.\n\
         These measurements do not certify a complete basis or regularity."
    );
}
