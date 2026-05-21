// End-to-end demo: pick a small curve, build a tiny factor base, deliberately
// construct R = P_a + P_b for known (a, b), then build the 2-decomposition
// polynomial system and run Buchberger on it. Filter spurious algebraic roots
// via the actual group law and confirm the planted decomposition is recovered.
//
//     cargo run --release --bin decomp_demo

use gbrl::buchberger::BuchbergerState;
use gbrl::curve::{add, filter_decomp, point_at_x, Point};
use gbrl::field::{Fp, P};
use gbrl::monomial::Monomial;
use gbrl::poly::{Poly, Term};
use gbrl::semaev::{decomposition_system_2, Curve};
use gbrl::strategy::{Strategy, SugarStrategy};

// Univariate factor-base polynomial F(x) = prod (x - v) for v in V.
fn factor_base_poly(values: &[Fp]) -> Poly {
    let mut acc = Poly::from_terms(
        vec![Term { coef: Fp::one(), mono: Monomial::one(1) }],
        1,
    );
    for &v in values {
        let factor = Poly::from_terms(
            vec![
                Term { coef:  Fp::one(), mono: Monomial { exp: vec![1] } },
                Term { coef: -v,         mono: Monomial { exp: vec![0] } },
            ],
            1,
        );
        acc = acc.mul(&factor);
    }
    acc
}

fn main() {
    let curve = Curve { a: Fp::new(1), b: Fp::new(1) };
    println!("curve: y^2 = x^3 + {}*x + {} over F_{}", curve.a, curve.b, P);

    let factor_xs: Vec<u64> = vec![0, 5, 11, 18];
    let fb_points: Vec<(Fp, Point)> = factor_xs.iter()
        .filter_map(|&xi| point_at_x(&curve, Fp::new(xi as i64)).map(|p| (Fp::new(xi as i64), p)))
        .collect();
    println!("\nfactor base ({} pts):", fb_points.len());
    for (x, pt) in &fb_points {
        println!("  x = {}  →  {:?}", x, pt);
    }
    let v_values: Vec<Fp> = fb_points.iter().map(|(x, _)| *x).collect();
    let fb_poly = factor_base_poly(&v_values);

    let a_idx = 0usize;
    let b_idx = 2usize;
    let p_a = fb_points[a_idx].1;
    let p_b = fb_points[b_idx].1;
    let r = add(&curve, p_a, p_b);
    let x_r = match r { Point::Affine(x, _) => x, _ => panic!("R = O") };
    println!("\nchosen: a_idx={} (x={}), b_idx={} (x={})",
             a_idx, fb_points[a_idx].0, b_idx, fb_points[b_idx].0);
    println!("R = P_a + P_b has x_R = {}", x_r);

    let sys = decomposition_system_2(&curve, x_r, &fb_poly);
    println!("\ndecomposition system:");
    println!("  S_3(x_R, x_1, x_2): {} terms, deg = {}", sys[0].nterms(), sys[0].total_degree());
    println!("  F(x_1):            {} terms, deg = {}", sys[1].nterms(), sys[1].total_degree());
    println!("  F(x_2):            {} terms, deg = {}", sys[2].nterms(), sys[2].total_degree());

    let mut state = BuchbergerState::new(sys);
    let mut strategy: Box<dyn Strategy> = Box::new(SugarStrategy);
    while !state.is_done() {
        let idx = strategy.select(&state);
        state.step(idx);
    }
    println!("\nBuchberger: steps={}, final basis={}, arith_ops={}",
             state.step_count, state.basis.len(), state.arith_ops);

    println!("\nbrute-force scan over V × V, classify via group-law filter:");
    let mut real_found = 0;
    for x1 in &v_values {
        for x2 in &v_values {
            let all_zero = state.basis.iter().all(|p| p.eval(&[*x1, *x2]).is_zero());
            if !all_zero { continue; }
            match filter_decomp(&curve, *x1, *x2, r) {
                Some((q1, q2)) => {
                    println!("  ({}, {})  → REAL: {:?} + {:?} = R (or -R)", x1, x2, q1, q2);
                    real_found += 1;
                }
                None => {
                    println!("  ({}, {})  → spurious (no sign choice yields ±R)", x1, x2);
                }
            }
        }
    }
    assert!(real_found >= 1, "should find at least one real decomposition");
    println!("\n{} real decomposition(s) found.", real_found);
}
