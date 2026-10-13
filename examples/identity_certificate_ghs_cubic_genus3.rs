//! Emit the two identity certificates of `research/ghs_cubic_genus3_build_20261012/`:
//! the polynomial (ideal-membership) forms of `DescentCycle.third_cover` and
//! `DescentCycle.cycle_closes`.
//!
//! ```sh
//! cargo run --example identity_certificate_ghs_cubic_genus3 -- research/ghs_cubic_genus3_build_20261012/certificates
//! ```
//!
//! The Lean lemmas derive the conclusion from the hypotheses `y0^2 = u u''`,
//! `y1^2 = u u'`, `u y2 = y0 y1`. A v1 certificate states a pure polynomial
//! identity, so each is stated as the explicit witness
//! `conclusion difference = sum (cofactor) * (hypothesis)`; the identity holds
//! for all `y0, y1, y2, u, u1, u2`, and it is the fact the Lean proof uses.

use crypto_lib::cryptanalysis::identity_certificate::{
    check, issue, Expr, IdentityStatement, ProofRef, Verdict,
};
use crypto_lib::hash::sha256::sha256;
use Expr::{Add, Mul, Neg};

fn v(s: &str) -> Expr {
    Expr::var(s)
}
fn sub(a: Expr, b: Expr) -> Expr {
    Add(vec![a, Neg(Box::new(b))])
}

fn main() {
    let dir = std::env::args().nth(1).expect("output directory");
    let vars: Vec<String> = ["u", "u1", "u2", "y0", "y1", "y2"].iter().map(|s| s.to_string()).collect();
    let (u, u1, u2, y0, y1, y2) = (v("u"), v("u1"), v("u2"), v("y0"), v("y1"), v("y2"));
    let uy2 = Mul(vec![u.clone(), y2.clone()]);
    let y0y1 = Mul(vec![y0.clone(), y1.clone()]);
    let h0 = sub(Mul(vec![y0.clone(), y0.clone()]), Mul(vec![u.clone(), u2.clone()])); // y0^2 - u u''
    let h1 = sub(Mul(vec![y1.clone(), y1.clone()]), Mul(vec![u.clone(), u1.clone()])); // y1^2 - u u'
    let h2 = sub(uy2.clone(), y0y1.clone()); // u y2 - y0 y1

    // third_cover: (u y2)^2 - u^2 u' u'' = (u y2 - y0 y1)(u y2 + y0 y1) + y1^2 (y0^2 - u u'') + u u'' (y1^2 - u u')
    let third = IdentityStatement {
        name: "third conjugate cover: (u y2)^2 - u^2 u1 u2 as a combination of the three cover relations".into(),
        source: "research/ghs_cubic_genus3_build_20261012 (README construction step 1); Lean DescentCycle.third_cover".into(),
        variables: vars.clone(),
        degree_bound: 4,
        lhs: sub(Mul(vec![uy2.clone(), uy2.clone()]), Mul(vec![u.clone(), u.clone(), u1.clone(), u2.clone()])),
        rhs: Add(vec![
            Mul(vec![h2.clone(), Add(vec![uy2.clone(), y0y1.clone()])]),
            Mul(vec![y1.clone(), y1.clone(), h0.clone()]),
            Mul(vec![u.clone(), u2.clone(), h1.clone()]),
        ]),
    };
    // cycle_closes: u y1 y2 - u u' y0 = y1 (u y2 - y0 y1) + y0 (y1^2 - u u')
    let cycle = IdentityStatement {
        name: "descent cycle closes: u y1 y2 - u u1 y0 as a combination of two cover relations".into(),
        source: "research/ghs_cubic_genus3_build_20261012 (README construction step 2); Lean DescentCycle.cycle_closes".into(),
        variables: vars,
        degree_bound: 3,
        lhs: sub(Mul(vec![u.clone(), y1.clone(), y2.clone()]), Mul(vec![u.clone(), u1.clone(), y0.clone()])),
        rhs: Add(vec![Mul(vec![y1.clone(), h2]), Mul(vec![y0.clone(), h1])]),
    };

    let lean_sha = {
        let bytes = std::fs::read("research/ghs_cubic_genus3_build_20261012/lean/DescentCycle.lean")
            .expect("run from the repository root");
        hex::encode(sha256(&bytes))
    };
    let proof = || {
        Some(ProofRef {
            kind: "lean4".into(),
            path: "research/ghs_cubic_genus3_build_20261012/lean/DescentCycle.lean".into(),
            sha256: lean_sha.clone(),
        })
    };
    for (file, stmt, seed) in [("third_cover.json", third, 0x6753_0000_0000_0001u64), ("cycle_closes.json", cycle, 0x6753_0000_0000_0002u64)] {
        let cert = issue(
            stmt,
            seed,
            32,
            "issued by examples/identity_certificate_ghs_cubic_genus3.rs; checked by cryptanalysis::identity_certificate::check replaying from the sealed id and seed. \
             The proof pointer is the Lean file whose lemma derives the same conclusion from the hypotheses; the checker records it and does not run Lean.",
            proof(),
        )
        .expect("identity holds");
        match check(&cert).expect("well-formed") {
            Verdict::Accept { .. } => {}
            other => panic!("self-check failed: {other:?}"),
        }
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(format!("{dir}/{file}"), serde_json::to_string_pretty(&cert).unwrap() + "\n").unwrap();
        println!("{file}: {}", cert.id);
    }
}
