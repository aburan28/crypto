//! Averaged refutation-degree (last-fall / PC-degree) sweep.
//!
//! Where `pc_degree_demo` reports one noisy `D*` per `n`, this demo
//! reports the `D*` *distribution* (min / mean / max + histogram) over a
//! batch of non-decomposable (unsatisfiable) targets per `n` — the stable
//! form that the proposal's prediction #1 actually talks about
//! (`RESEARCH_FFD_PROOF_COMPLEXITY.md` §6): mean `D*` should climb with
//! `n` while the first-fall mean stays roughly flat.
//!
//! ```bash
//! cargo run --example pc_degree_avg_demo --release
//! ```

use crypto_lib::cryptanalysis::pc_degree_avg::{print_pc_avg_sweep, run_pc_avg_sweep};

fn main() {
    // 64 non-decomposable targets per n, degree budget 12.
    let rows = run_pc_avg_sweep(4..=8, 12, 64, 0x_FFD_A6_DE);
    print_pc_avg_sweep(&rows);
}
