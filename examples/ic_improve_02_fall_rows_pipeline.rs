//! IC improvement item 02, lead probe (not the item's candidate): does the
//! attribution arm F — the root's degree-3 Macaulay degree-fall rows appended
//! to the generators — also cut the incumbent oracle's cost on the
//! pipeline's own target-fixed systems (`build_decomposition_system`)?
//!
//! Same solver and options as the pipeline (`SolverEngine::default()`,
//! `split_rule_default()`); m = 2 or 3; subspace V = span(1, z, …, z^{ℓ−1}).
//! Gate: identical root sets.
//!
//! ```bash
//! cargo run --release --example ic_improve_02_fall_rows_pipeline -- --n 17 --m 3 --ell 6 --targets 8
//! ```
use crypto_lib::binary_ecc::f2m::F2mElement;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, f4_word_ops_thread, matrix_f4_f2_counted_with,
    solve_boolean_system, split_rule_default, FieldStructure, RowCriterion, SolveOptions,
    SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible;
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn tcpu() -> f64 {
    let mut ts = libc::timespec { tv_sec: 0, tv_nsec: 0 };
    unsafe {
        libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut ts);
    }
    ts.tv_sec as f64 + ts.tv_nsec as f64 * 1e-9
}

fn solve(sys: &[F2BoolPoly], nv: usize) -> (Vec<u64>, f64, u64, usize, bool) {
    let opts = SolveOptions {
        engine: SolverEngine::default(),
        max_solutions: usize::MAX,
        node_budget: usize::MAX,
        split_rule: split_rule_default(),
    };
    let w0 = f4_word_ops_thread();
    let c0 = tcpu();
    let (mut roots, st) = solve_boolean_system(sys, nv, &opts);
    let cpu = tcpu() - c0;
    roots.sort_unstable();
    (roots, cpu, f4_word_ops_thread() - w0, st.splits, st.exhausted)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let get = |k: &str, d: u64| -> u64 {
        a.iter().position(|x| x == k).map(|i| a[i + 1].parse().unwrap()).unwrap_or(d)
    };
    let (n, m, ell, targets) = (get("--n", 13) as u32, get("--m", 3) as usize, get("--ell", 5) as u32, get("--targets", 4));
    let irr = find_irreducible(n).unwrap();
    let st = FieldStructure::new(n, &irr);
    let basis: Vec<F2mElement> = (0..ell).map(|k| F2mElement::from_bit_positions(&[k], n)).collect();
    let b = F2mElement::from_bit_positions(&[0], n);
    for t in 0..targets {
        let mut rng = StdRng::seed_from_u64(20261014 ^ (n as u64 * 1000 + t) ^ 0xfa11);
        let xr_bits: u64 = rng.gen_range(1..(1u64 << n));
        let pos: Vec<u32> = (0..n).filter(|k| (xr_bits >> k) & 1 == 1).collect();
        let xr = F2mElement::from_bit_positions(&pos, n);
        let sys = build_decomposition_system(&basis, &xr, &b, m, &st).expect("system");
        let nv = sys.n_vars;
        let (i_roots, i_cpu, i_ops, i_splits, i_exh) = solve(&sys.equations, nv);
        let c0 = tcpu();
        let (rows, step_ops) = matrix_f4_f2_counted_with(&sys.equations, nv, 3, RowCriterion::None).expect("matrix");
        let mut fsys = sys.equations.clone();
        let extra: Vec<F2BoolPoly> = rows.into_iter().filter(|p| p.terms.iter().all(|t| t.mask.count_ones() < 3)).collect();
        let n_extra = extra.len();
        fsys.extend(extra);
        let f_step = tcpu() - c0;
        let (f_roots, f_solve, f_ops, f_splits, f_exh) = solve(&fsys, nv);
        let (i2_roots, i2_cpu, _, _, _) = solve(&sys.equations, nv);
        let gate = !i_exh && !f_exh && i_roots == f_roots && i2_roots == i_roots;
        println!(
            r#"{{"probe":"fall_rows_pipeline","n":{n},"m":{m},"ell":{ell},"target":{t},"x_r":{xr_bits},"n_vars":{nv},"roots":{},"i_cpu":{i_cpu:.5},"i2_cpu":{i2_cpu:.5},"i_word_ops":{i_ops},"i_splits":{i_splits},"f_cpu":{:.5},"f_step_cpu":{f_step:.5},"f_step_word_ops":{step_ops},"f_word_ops":{f_ops},"f_extra_rows":{n_extra},"f_splits":{f_splits},"gate":{gate}}}"#,
            i_roots.len(),
            f_step + f_solve
        );
    }
}
