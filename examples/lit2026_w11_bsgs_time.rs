//! lit2026 W11 (a): wall time of the single-word BSGS
//! ([`crypto_lib::cryptanalysis::bsgs_fast`]) on one curve, full group,
//! negation map, for the interleaved original / previous / candidate arms.
//!
//! Uses only API present at the program baseline (`4b6e09b78`), so the same
//! source builds there for the `original` arm.  On this branch the step
//! kernel comes from `BSGS_STEP` (`chain` = previous, `paired`, `coz`) and
//! `BSGS_PAIR_K`; at the baseline those variables are ignored.
//!
//! ```sh
//! BSGS_STEP=paired cargo run --release --example lit2026_w11_bsgs_time -- \
//!     <p> <a> <b> <gx> <gy> <n> <targets> <seed> <label>
//! ```
//!
//! Prints one JSON line: build and search seconds, step counts, and a digest
//! of the verified answers (identical across arms by construction).

use std::time::Instant;

use num_bigint::BigUint;

use crypto_lib::cryptanalysis::bsgs_fast::{BsgsFast, BsgsFastPlan, FastCurve, FastPoint};
use crypto_lib::ecc::curve::CurveParams;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    assert!(a.len() >= 9, "usage: p a b gx gy n targets seed label");
    let u = |i: usize| -> u64 { a[i].parse().expect("integer") };
    let params = CurveParams {
        name: "lit2026-w11",
        p: BigUint::from(u(0)),
        a: BigUint::from(u(1)),
        b: BigUint::from(u(2)),
        gx: BigUint::from(u(3)),
        gy: BigUint::from(u(4)),
        n: BigUint::from(u(5)),
        h: 1,
    };
    let targets = u(6) as usize;
    let mut seed = u(7) | 1;
    let label = &a[8];
    let curve = FastCurve::from_params(&params).expect("single-word curve on which G lies");
    let threads = rayon::current_num_threads();
    let plan = BsgsFastPlan::new(curve.n, 0, true, threads, 8);
    let (m, table_bits) = (plan.m, plan.table_bits);

    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    let secrets: Vec<u64> = (0..targets).map(|_| next() % curve.n).collect();
    let qs: Vec<FastPoint> = secrets
        .iter()
        .map(|&x| curve.scalar_mul(curve.g, x))
        .collect();

    let t0 = Instant::now();
    let mut solver = BsgsFast::build(curve, plan);
    let build_s = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let got = solver.solve_many(&qs);
    let search_s = t1.elapsed().as_secs_f64();
    let solved = got
        .iter()
        .zip(&secrets)
        .filter(|(g, &x)| **g == Some(x))
        .count();
    assert_eq!(solved, targets, "every target must verify");
    let digest = secrets.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, &x| {
        (h ^ x).wrapping_mul(0x100_0000_01b3)
    });
    let st = solver.stats();
    let kernel = std::env::var("BSGS_STEP").unwrap_or_else(|_| "unset".into());
    let pair_k = std::env::var("BSGS_PAIR_K").unwrap_or_else(|_| "unset".into());
    println!(
        "{{\"label\":\"{label}\",\"bsgs_step_env\":\"{kernel}\",\"bsgs_pair_k_env\":\"{pair_k}\",\"n\":{},\"bits\":{},\"threads\":{threads},\"m\":{m},\"table_bits\":{table_bits},\"targets\":{targets},\"solved\":{solved},\"answers_digest\":\"{digest:016x}\",\"table_entries\":{},\"baby_steps\":{},\"giant_steps\":{},\"candidates\":{},\"false_candidates\":{},\"build_s\":{build_s:.4},\"search_s\":{search_s:.4},\"total_s\":{:.4}}}",
        curve.n,
        64 - curve.n.leading_zeros(),
        st.table_entries,
        st.baby_steps,
        st.giant_steps,
        st.candidates,
        st.false_candidates,
        build_s + search_s
    );
}
