//! Wall-clock benchmark for the Boolean-ring F4 (`cryptanalysis::pq_f4_f2`)
//! and matrix-F5 steps (`cryptanalysis::matrix_f5_f2`)
//! on random quadratic systems with a planted solution, fixed seeds.
//!
//! Prints one JSON line per case: median wall time over `repeats` runs, the
//! build/eliminate split, the counted units, and a fingerprint of the
//! reduced basis so two revisions can be compared for identical output.
//!
//! ```text
//! cargo run --release --example f4_f2_bench -- [repeats] [max_n] [families] [seed_xor_hex]
//! ```
//!
//! `families` is a comma-separated subset of `f4`, `bb` (Buchberger) and
//! `f5` (matrix-F5 steps), all three by default: `-- 5 24 f5` profiles the
//! F5 steps alone.  The F5 lines also split their wall time into the
//! criterion, the row construction, the elimination and the unpacking.
//! The optional hexadecimal seed XOR selects a fresh planted system;
//! omission preserves the original fixed suite exactly.

use crypto_lib::cryptanalysis::matrix_f5_f2::{
    canonical_row_space_fingerprint, matrix_f5_f2_with_form_timed, F5OutputForm,
};
use crypto_lib::cryptanalysis::pq_f4_f2::groebner_basis_f4;
use crypto_lib::cryptanalysis::pq_groebner_f2::{groebner_basis_f2_stats, F2BoolMono, F2BoolPoly};
use serde_json::json;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

fn system(n: usize, m: usize, seed: u64) -> Vec<F2BoolPoly> {
    let mut s = seed | 1;
    let mut rnd = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    let planted = rnd() & ((1u64 << n) - 1);
    let mut monos: Vec<u64> = vec![0];
    for i in 0..n {
        monos.push(1 << i);
        for j in i + 1..n {
            monos.push((1 << i) | (1 << j));
        }
    }
    (0..m)
        .map(|_| {
            let chosen: Vec<u64> = monos.iter().copied().filter(|_| rnd() & 1 == 1).collect();
            // fix the constant so the planted point is a zero
            let val = chosen.iter().filter(|&&t| t & planted == t).count() & 1;
            let mut terms: Vec<F2BoolMono> =
                chosen.into_iter().map(F2BoolMono::from_mask).collect();
            if val == 1 {
                // toggles the constant term: `from_monos` cancels pairs
                terms.push(F2BoolMono::from_mask(0));
            }
            F2BoolPoly::from_monos(terms, n)
        })
        .collect()
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.total_cmp(b));
    v[v.len() / 2]
}

fn main() {
    let repeats: usize = std::env::args()
        .nth(1)
        .map(|s| s.parse().expect("repeats"))
        .unwrap_or(3);
    // optional: skip cases with more than this many variables
    let max_n: usize = std::env::args()
        .nth(2)
        .map(|s| s.parse().expect("max_n"))
        .unwrap_or(usize::MAX);
    let families: Vec<String> = std::env::args()
        .nth(3)
        .unwrap_or_else(|| "f4,bb,f5".into())
        .split(',')
        .map(str::to_owned)
        .collect();
    let seed_xor = std::env::args()
        .nth(4)
        .map(|s| u64::from_str_radix(s.trim_start_matches("0x"), 16).expect("seed_xor_hex"))
        .unwrap_or(0);
    let run = |family: &str| families.iter().any(|f| f == family);
    for &(n, m) in &[
        (12usize, 12usize),
        (14, 14),
        (16, 16),
        (16, 24),
        (18, 18),
        (20, 20),
        (20, 30),
    ] {
        if n > max_n || !run("f4") {
            continue;
        }
        let sys = system(
            n,
            m,
            0x2545_f491_4f6c_dd1d ^ ((n * 64 + m) as u64) ^ seed_xor,
        );
        let mut walls = Vec::new();
        let mut last = None;
        for _ in 0..repeats {
            let (gb, st) = groebner_basis_f4(sys.clone(), n, None);
            walls.push(st.wall_ns as f64 / 1e6);
            last = Some((gb, st));
        }
        let (gb, st) = last.unwrap();
        let mut h = DefaultHasher::new();
        for p in &gb {
            for t in &p.terms {
                t.mask.hash(&mut h);
            }
            u64::MAX.hash(&mut h);
        }
        println!(
            "{}",
            json!({
                "case": format!("n{n}_m{m}"), "wall_ms": median(walls),
                "build_ms": st.build_ns as f64 / 1e6, "eliminate_ms": st.eliminate_ns as f64 / 1e6,
                "other_ms": st.wall_ns.saturating_sub(st.build_ns.saturating_add(st.eliminate_ns)) as f64 / 1e6,
                "word_xors": st.word_xors, "word_xors_performed": st.word_xors_performed, "divisor_tests": st.divisor_tests,
                "steps": st.steps, "basis_len": st.basis_len, "basis_fp": format!("{:016x}", h.finish()),
                "matrix_rows_max": st.matrix_rows_max, "matrix_cols_max": st.matrix_cols_max,
                "peak_matrix_bytes": st.peak_matrix_bytes, "peak_table_bytes": st.peak_table_bytes,
                "pairs_reduced": st.pairs_reduced, "field_pairs_reduced": st.field_pairs_reduced,
                "pairs_chain_skipped": st.pairs_chain_skipped, "pairs_product_skipped": st.pairs_product_skipped,
                "pairs_left": st.pairs_left, "new_elements": st.new_elements,
            })
        );
    }
    // Buchberger (the reference engine) on the smaller systems
    for &(n, m) in &[(8usize, 8usize), (10, 10), (12, 12), (12, 16), (14, 14)] {
        if n > max_n || !run("bb") {
            continue;
        }
        let sys = system(
            n,
            m,
            0x2545_f491_4f6c_dd1d ^ ((n * 64 + m) as u64) ^ seed_xor,
        );
        let mut walls = Vec::new();
        let mut last = None;
        for _ in 0..repeats {
            let (gb, st) = groebner_basis_f2_stats(sys.clone(), n);
            walls.push(st.wall_ns as f64 / 1e6);
            last = Some((gb, st));
        }
        let (gb, st) = last.unwrap();
        let mut h = DefaultHasher::new();
        for p in &gb {
            for t in &p.terms {
                t.mask.hash(&mut h);
            }
            u64::MAX.hash(&mut h);
        }
        println!(
            "{}",
            json!({
                "case": format!("bb_n{n}_m{m}"), "wall_ms": median(walls),
                "mono_ops": st.mono_ops, "reduction_steps": st.reduction_steps,
                "spolys": st.spolys, "pairs_chain_skipped": st.pairs_chain_skipped,
                "field_pairs": st.field_pairs, "basis_len": st.basis_len,
                "basis_fp": format!("{:016x}", h.finish()),
            })
        );
    }
    // matrix-F5 steps (criterion + elimination) on the same systems
    for &(n, m, degree) in &[
        (12usize, 12usize, 4u32),
        (16, 16, 3),
        (16, 16, 4),
        (20, 20, 3),
        (20, 20, 4),
        (24, 24, 3),
        (24, 24, 4),
    ] {
        if n > max_n || !run("f5") {
            continue;
        }
        let sys = system(
            n,
            m,
            0x2545_f491_4f6c_dd1d ^ ((n * 64 + m) as u64) ^ seed_xor,
        );
        let mut walls = Vec::new();
        let mut last = None;
        for _ in 0..repeats {
            let t = std::time::Instant::now();
            let form = match std::env::var("KIC_F5_ECHELON").as_deref() {
                Ok("1") => F5OutputForm::Echelon,
                Ok("2") => F5OutputForm::SelectiveEchelon,
                _ => F5OutputForm::Reduced,
            };
            let r =
                matrix_f5_f2_with_form_timed(&sys, n, degree, form).expect("within size limits");
            walls.push(t.elapsed().as_secs_f64() * 1e3);
            last = Some(r);
        }
        let (rows, rep, phases) = last.unwrap();
        let output_terms: usize = rows.iter().map(|p| p.terms.len()).sum();
        let row_space_fp = canonical_row_space_fingerprint(&rows)
            .expect("returned F5 rows have a valid column set");
        let mut h = DefaultHasher::new();
        for p in &rows {
            for t in &p.terms {
                t.mask.hash(&mut h);
            }
            u64::MAX.hash(&mut h);
        }
        println!(
            "{}",
            json!({
                "case": format!("f5_n{n}_m{m}_d{degree}"), "wall_ms": median(walls),
                "criterion_word_ops": rep.criterion_word_ops, "reduce_word_ops": rep.reduce_word_ops,
                "rows_f4": rep.rows_f4, "rows_built": rep.rows_built, "cols": rep.cols,
                "rows_pruned": rep.rows_pruned, "rank": rep.rank, "rows_fp": format!("{:016x}", h.finish()),
                "output_terms": output_terms,
                "direct_pack_used": phases.direct_pack_used,
                "direct_unpack_used": phases.direct_unpack_used,
                "compact_unpack_used": phases.compact_unpack_used,
                "row_space_fp": format!("{row_space_fp:016x}"),
                "criterion_ms": phases.criterion_ns as f64 / 1e6, "f5_build_ms": phases.build_ns as f64 / 1e6,
                "reduce_ms": phases.reduce_ns as f64 / 1e6, "unpack_ms": phases.unpack_ns as f64 / 1e6,
            })
        );
    }
}
