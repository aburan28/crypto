//! Planted Semaev m=3 decomposition with group-sum witness check + F_p FFD.
//!
//! Public synthetic on y²=x³+x+1 / F_1000003 (|E|≈2^20 ≥ 14 bits).
//! Not a ledger promotion.
//!
//! ```bash
//! cargo run --release --bin m3_witness_ffd
//! ```

use gbrl::buchberger::BuchbergerState;
use gbrl::curve::{add, filter_decomp_3, group_order, point_at_x, Point};
use gbrl::field::{Fp, P};
use gbrl::ffd::first_fall_degree;
use gbrl::monomial::Monomial;
use gbrl::poly::{Poly, Term};
use gbrl::semaev::{decomposition_system_3, Curve};
use gbrl::strategy::{Strategy, SugarStrategy};
use std::time::Instant;

fn factor_base_poly(values: &[Fp]) -> Poly {
    let mut acc = Poly::from_terms(
        vec![Term {
            coef: Fp::one(),
            mono: Monomial::one(1),
        }],
        1,
    );
    for &v in values {
        let factor = Poly::from_terms(
            vec![
                Term {
                    coef: Fp::one(),
                    mono: Monomial { exp: vec![1] },
                },
                Term {
                    coef: -v,
                    mono: Monomial { exp: vec![0] },
                },
            ],
            1,
        );
        acc = acc.mul(&factor);
    }
    acc
}

fn first_n_curve_xs(curve: &Curve, n: usize) -> Vec<Fp> {
    let mut out = Vec::new();
    let mut xi = 0u64;
    while out.len() < n && xi < 50_000 {
        let x = Fp::new(xi as i64);
        if point_at_x(curve, x).is_some() {
            out.push(x);
        }
        xi += 1;
    }
    out
}

fn x_of(p: Point) -> Fp {
    match p {
        Point::Affine(x, _) => x,
        Point::Infinity => panic!("expected affine"),
    }
}

fn main() {
    let curve = Curve {
        a: Fp::new(1),
        b: Fp::new(1),
    };
    let order = group_order(&curve);
    let order_bits = 64 - order.leading_zeros();
    let fb_size: usize = std::env::var("M3_WITNESS_FB")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(4);
    let d_max: u32 = std::env::var("M3_FFD_DMAX")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(16);

    let fb_xs = first_n_curve_xs(&curve, fb_size);
    let fb_pts: Vec<Point> = fb_xs
        .iter()
        .map(|&x| point_at_x(&curve, x).expect("fb x on curve"))
        .collect();
    let fb_poly = factor_base_poly(&fb_xs);

    assert!(fb_pts.len() >= 3);
    let planted = (fb_pts[0], fb_pts[1], fb_pts[2]);
    let r = add(&curve, add(&curve, planted.0, planted.1), planted.2);
    let x_r = x_of(r);
    let planted_xs = [x_of(planted.0), x_of(planted.1), x_of(planted.2)];

    let sys = decomposition_system_3(&curve, x_r, &fb_poly);
    let nvars = 3usize;
    let eqs = sys.len();
    let system_degree = sys.iter().map(|p| p.total_degree()).max().unwrap_or(0);
    let eq_var_ratio = eqs as f64 / nvars as f64;

    let t_ffd = Instant::now();
    let (ffd, profiles) = first_fall_degree(&sys, nvars, d_max);
    let ffd_ms = t_ffd.elapsed().as_secs_f64() * 1e3;

    let t_gb = Instant::now();
    let mut st = BuchbergerState::new(sys.clone());
    let mut strategy: Box<dyn Strategy> = Box::new(SugarStrategy);
    while !st.is_done() && st.step_count < 200_000 {
        let idx = strategy.select(&st);
        st.step(idx);
    }
    let gb_ms = t_gb.elapsed().as_secs_f64() * 1e3;

    let mut algebraic_zeros = 0usize;
    let mut group_witnesses = 0usize;
    let mut recovered_planted = false;
    for &x1 in &fb_xs {
        for &x2 in &fb_xs {
            for &x3 in &fb_xs {
                let pt = [x1, x2, x3];
                if !st.basis.iter().all(|p| p.eval(&pt).is_zero()) {
                    continue;
                }
                algebraic_zeros += 1;
                if filter_decomp_3(&curve, x1, x2, x3, r).is_some() {
                    group_witnesses += 1;
                    let mut xs = [x1, x2, x3];
                    xs.sort_by_key(|v| v.0);
                    let mut planted_sorted = planted_xs;
                    planted_sorted.sort_by_key(|v| v.0);
                    if xs == planted_sorted {
                        recovered_planted = true;
                    }
                }
            }
        }
    }

    let planted_vanishes_input = sys.iter().all(|p| p.eval(&planted_xs).is_zero());
    let planted_group_ok =
        filter_decomp_3(&curve, planted_xs[0], planted_xs[1], planted_xs[2], r).is_some();

    println!(
        "{}",
        serde_json::json!({
            "kind": "m3_witness_ffd",
            "field_p": P,
            "curve": "y^2=x^3+x+1",
            "group_order": order,
            "group_order_bits": order_bits,
            "bits_gate_ge14": order_bits >= 14,
            "fb_size": fb_size,
            "nvars": nvars,
            "eqs": eqs,
            "eq_var_ratio": eq_var_ratio,
            "system_degree": system_degree,
            "buchberger_steps": st.step_count,
            "basis_len": st.basis.len(),
            "arith_ops": st.arith_ops,
            "gb_done": st.is_done(),
            "gb_wall_ms": gb_ms,
            "ffd_status": if ffd.is_some() { "measured" } else { "not_found_within_dmax" },
            "ffd": ffd,
            "ffd_d_max": d_max,
            "ffd_wall_ms": ffd_ms,
            "macaulay_last": profiles.last().map(|p| serde_json::json!({
                "degree": p.degree, "rows": p.rows, "cols": p.cols, "rank": p.rank
            })),
            "planted_vanishes_input": planted_vanishes_input,
            "planted_group_sum_ok": planted_group_ok,
            "algebraic_zeros_over_fb": algebraic_zeros,
            "group_verified_witnesses": group_witnesses,
            "recovered_planted_up_to_perm": recovered_planted,
            "witnesses_sum_in_group": group_witnesses > 0 && planted_group_ok,
            "claim_boundary": "Public synthetic planted m=3 on toy curve |E|~2^20; not key recovery; not ledger promotion.",
        })
    );
}
