//! Seal the inverse-Velu cubic identity with a native replayable certificate.

use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr, IdentityStatement, ProofRef, Verdict,
};

fn v(name: &str) -> Expr {
    Expr::var(name)
}
fn c(value: i64) -> Expr {
    Expr::Const(value)
}
fn add(terms: Vec<Expr>) -> Expr {
    Expr::sum(terms)
}
fn mul(terms: Vec<Expr>) -> Expr {
    Expr::product(terms)
}
fn neg(value: Expr) -> Expr {
    Expr::minus(value)
}

fn cubic() -> IdentityStatement {
    let x3 = || mul(vec![v("x"), v("x"), v("x")]);
    let lhs = add(vec![
        x3(),
        neg(mul(vec![add(vec![v("u"), c(6)]), v("x"), v("x")])),
        mul(vec![add(vec![mul(vec![c(6), v("u")]), c(21)]), v("x")]),
        neg(add(vec![mul(vec![c(9), v("u")]), c(32)])),
    ]);
    let cross = mul(vec![
        add(vec![v("u"), v("x"), c(12)]),
        add(vec![v("x"), c(-3)]),
        add(vec![v("x"), c(-3)]),
    ]);
    let curve_rhs = add(vec![x3(), neg(mul(vec![c(21), v("x")])), c(38)]);
    let rhs = neg(add(vec![cross, neg(mul(vec![c(2), curve_rhs]))]));
    IdentityStatement {
        name: "inverse-Velu cubic for kernel (3,1)".to_owned(),
        source: "research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/inverse_velu/PROTOCOL_V2.md; Lean cubicIdentity".to_owned(),
        variables: vec!["x".to_owned(), "u".to_owned()],
        degree_bound: 3,
        lhs,
        rhs,
    }
}

fn main() {
    let cert = issue(
        cubic(),
        0x2026_1016_1201_0002,
        32,
        "issued by native inverse-Velu example; checked by identity_certificate::check_with_extra on a separate 32-point stream; integer algebra proved by Lean cubicIdentity; finite-field curve bridge separately checked by Rust and remains unformalized",
        Some(ProofRef {
            kind: "lean4".to_owned(),
            path: "research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/inverse_velu/lean/InverseVelu.lean".to_owned(),
            sha256: "cac8b9f311f48353df1476f374dd75e50a2d0cbbbaf69ba7a0c0c7dee8d3be20".to_owned(),
        }),
    )
    .expect("integer cubic identity must issue");
    match check_with_extra(&cert, Some((0xfeed_1201_2026_1016, 32))).expect("valid certificate") {
        Verdict::Accept {
            points: 32,
            extra_points: 32,
            ..
        } => {}
        other => panic!("identity replay failed: {other:?}"),
    }
    println!("{}", serde_json::to_string_pretty(&cert).unwrap());
}
