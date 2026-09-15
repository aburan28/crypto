//! Timed Semaev 3-decomposition Buchberger + F_p FFD probe (gbrl).
//!
//! ```bash
//! cargo run --release --bin m3_decomp_bench
//! M3_FFD_DMAX=24 cargo run --release --bin m3_decomp_bench
//! ```
//!
//! Public synthetic only. Field is fixed F_1000003 (~20-bit), not a
//! ≥14-bit curve-order ledger claim.

use gbrl::buchberger::BuchbergerState;
use gbrl::ffd::first_fall_degree;
use gbrl::semaev::semaev_decomposition_3_instance;
use gbrl::strategy::{Strategy, SugarStrategy};
use std::time::Instant;

fn main() {
    println!("kind=m3_decomp_bench field=F_1000003 oracle=buchberger_sugar+fp_ffd");
    let d_max = std::env::var("M3_FFD_DMAX")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(24u32);
    let mut ffds: Vec<u32> = Vec::new();
    for fb in [3usize, 4, 5, 6] {
        let sys = semaev_decomposition_3_instance(fb);
        let nvars = sys.first().map(|p| p.nvars).unwrap_or(0);
        let deg = sys.iter().map(|p| p.total_degree()).max().unwrap_or(0);
        let eqs = sys.len();
        let eq_var_ratio = if nvars == 0 {
            0.0
        } else {
            eqs as f64 / nvars as f64
        };

        let t_ffd = Instant::now();
        let (ffd, profiles) = first_fall_degree(&sys, nvars, d_max);
        let ffd_ms = t_ffd.elapsed().as_secs_f64() * 1e3;
        if let Some(d) = ffd {
            ffds.push(d);
        }

        let t0 = Instant::now();
        let mut st = BuchbergerState::new(sys);
        let mut strategy: Box<dyn Strategy> = Box::new(SugarStrategy);
        while !st.is_done() && st.step_count < 200_000 {
            let idx = strategy.select(&st);
            st.step(idx);
        }
        let ms = t0.elapsed().as_secs_f64() * 1e3;
        let last = profiles.last();
        println!(
            "{}",
            serde_json::json!({
                "fb_size": fb,
                "nvars": nvars,
                "eqs": eqs,
                "eq_var_ratio": eq_var_ratio,
                "system_degree": deg,
                "buchberger_steps": st.step_count,
                "basis_len": st.basis.len(),
                "arith_ops": st.arith_ops,
                "done": st.is_done(),
                "wall_ms": ms,
                "ffd_status": if ffd.is_some() { "measured" } else { "not_found_within_dmax" },
                "ffd": ffd,
                "ffd_d_max": d_max,
                "ffd_wall_ms": ffd_ms,
                "macaulay_profiles": profiles.iter().map(|p| serde_json::json!({
                    "degree": p.degree,
                    "rows": p.rows,
                    "cols": p.cols,
                    "rank": p.rank,
                    "syzygies": p.syzygies(),
                })).collect::<Vec<_>>(),
                "macaulay_last": last.map(|p| serde_json::json!({
                    "degree": p.degree,
                    "rows": p.rows,
                    "cols": p.cols,
                    "rank": p.rank,
                })),
                "bits_note": "field is fixed F_1000003 (~20-bit prime); not a 14-bit curve order claim",
            })
        );
    }
    if !ffds.is_empty() {
        let min = *ffds.iter().min().unwrap();
        let max = *ffds.iter().max().unwrap();
        let mean = ffds.iter().sum::<u32>() as f64 / ffds.len() as f64;
        println!(
            "{}",
            serde_json::json!({
                "kind": "ffd_summary",
                "draws": ffds.len(),
                "ffd_min": min,
                "ffd_max": max,
                "ffd_mean": mean,
                "ffds": ffds,
            })
        );
    }
}
