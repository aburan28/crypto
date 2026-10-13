//! Emit the identity certificate for the weak-twist normal form used in
//! `research/weak_subfield_hyperelliptic_decomposition_20261012/` as JSON on
//! stdout.
//!
//! ```sh
//! cargo run --release --example identity_certificate_weak_twist \
//!   > research/weak_subfield_hyperelliptic_decomposition_20261012/weak-twist-normal-form.certificate.json
//! ```
//!
//! Statement: for `f(x) = x^5 + sum_{k<5} g_k beta^(5-k) x^k`,
//! `f(beta t) = beta^5 (t^5 + sum_{k<5} g_k t^k)`. A curve `y^2 = f(x)` whose
//! coefficients have this shape with `g_k` in a subfield is therefore, on
//! the line `x = beta t`, a scalar multiple of a polynomial over that
//! subfield.

use crypto_lib::cryptanalysis::identity_certificate::{
    check, issue, Expr, IdentityStatement, Verdict,
};

fn pow(base: &str, e: usize) -> Expr {
    Expr::product((0..e).map(|_| Expr::var(base)).collect())
}

fn statement() -> IdentityStatement {
    let g = |k: usize| Expr::var(&format!("g{k}"));
    // lhs: (beta t)^5 + sum_k (g_k beta^(5-k)) (beta t)^k
    let bt = |k: usize| Expr::product((0..k).map(|_| Expr::times(Expr::var("beta"), Expr::var("t"))).collect());
    let mut lhs = vec![bt(5)];
    let mut inner = vec![pow("t", 5)];
    for k in 0..5 {
        lhs.push(Expr::product(vec![g(k), pow("beta", 5 - k), bt(k)]));
        inner.push(Expr::times(g(k), pow("t", k)));
    }
    IdentityStatement {
        name: "weak-twist normal form: f(beta t) = beta^5 g(t)".to_string(),
        source: "research/weak_subfield_hyperelliptic_decomposition_20261012/README.md, Lemma 1; \
                 Bhati-Singh, A new perspective on the decomposition in the Jacobian of small genus \
                 hyperelliptic curve, Eq. (11) and (16)"
            .to_string(),
        variables: ["beta", "t", "g0", "g1", "g2", "g3", "g4"].iter().map(|s| s.to_string()).collect(),
        degree_bound: 10,
        lhs: Expr::sum(lhs),
        rhs: Expr::times(pow("beta", 5), Expr::sum(inner)),
    }
}

fn main() {
    let cert = issue(
        statement(),
        0x2026_1012_0000_2193,
        32,
        "issued by examples/identity_certificate_weak_twist.rs; checked by \
         cryptanalysis::identity_certificate::check replaying from the sealed id and seed",
        None,
    )
    .expect("the normal form holds");
    match check(&cert).expect("well-formed") {
        Verdict::Accept { .. } => {}
        other => panic!("self-check failed: {other:?}"),
    }
    println!("{}", serde_json::to_string_pretty(&cert).expect("serialises"));
}
