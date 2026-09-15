//! Bounded public-synthetic Macaulay matrix probe for the m=3 ladder.
//!
//! This measures the initial S-polynomial matrix after symbolic preprocessing.
//! It is not a complete F4 Gröbner solve and does not claim an FFD/DoR result.

use gbrl::f4::F4State;
use gbrl::f4::{matrix_reduce, symbolic_preprocess};
use gbrl::semaev::semaev_decomposition_3_instance;
use gbrl::spoly::spoly;
use std::time::Instant;

fn main() {
    println!("kind=m3_f4_matrix_probe field=F_1000003 stage=initial_macaulay_only");
    let max_fb = std::env::var("M3_F4_MAX_FB")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3usize);
    let round_cap = std::env::var("M3_F4_ROUND_CAP")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(100usize);
    let matrix_row_cap = std::env::var("M3_F4_MATRIX_ROW_CAP")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2_000usize);
    let matrix_col_cap = std::env::var("M3_F4_MATRIX_COL_CAP")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2_000usize);
    for fb in [3usize, 4, 5, 6].into_iter().filter(|fb| *fb <= max_fb) {
        let basis = semaev_decomposition_3_instance(fb);
        let mut s_polys = Vec::new();
        for i in 0..basis.len() {
            for j in (i + 1)..basis.len() {
                let (Some(lmi), Some(lmj)) = (basis[i].lm(), basis[j].lm()) else {
                    continue;
                };
                if lmi.gcd_is_one(lmj) {
                    continue;
                }
                s_polys.push(spoly(&basis[i], &basis[j]));
            }
        }
        let t0 = Instant::now();
        let (rows, columns) = symbolic_preprocess(s_polys, &basis);
        let matrix_within_cap = rows.len() <= matrix_row_cap;
        let reduced = if matrix_within_cap {
            matrix_reduce(rows.clone(), &columns)
        } else {
            Vec::new()
        };
        let mut f4 = F4State::new(basis.clone());
        let mut rounds = 0usize;
        let cap = round_cap;
        let mut f4_resource_exhausted = false;
        while !f4.is_done() && rounds < cap {
            let result = f4.step_bounded(0, matrix_row_cap, matrix_col_cap);
            rounds += 1;
            if result.resource_exhausted {
                f4_resource_exhausted = true;
                break;
            }
            if result.matrix_rows == 0 && result.added == 0 {
                break;
            }
        }
        let wall_ms = t0.elapsed().as_secs_f64() * 1e3;
        println!(
            "{}",
            serde_json::json!({
                "fb_size": fb,
                "nvars": basis.first().map(|p| p.nvars).unwrap_or(0),
                "equations": basis.len(),
                "initial_s_polynomials": rows.len(),
                "matrix_rows": rows.len(),
                "matrix_columns": columns.len(),
                "reduced_nonzero_rows": reduced.len(),
                "f4_rounds": rounds,
                "f4_basis_len": f4.basis.len(),
                "f4_pairs_remaining": f4.pairs.len(),
                "f4_complete": f4.is_done(),
                "f4_round_cap": cap,
                "f4_status": if !matrix_within_cap || f4_resource_exhausted { "resource_exhaustion_matrix_cap" } else if f4.is_done() { "completed" } else { "resource_exhaustion_round_cap" },
                "matrix_row_cap": matrix_row_cap,
                "matrix_col_cap": matrix_col_cap,
                "matrix_within_cap": matrix_within_cap,
                "max_input_degree": basis.iter().map(|p| p.total_degree()).max().unwrap_or(0),
                "matrix_wall_ms": wall_ms,
                "ffd_status": "bounded_initial_matrix_only",
                "claim_boundary": "Public synthetic matrix-shape probe only; not a complete F4 solve, not FFD/DoR, not key recovery, and not ledger promotion."
            })
        );
    }
}
