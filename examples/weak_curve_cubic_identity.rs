//! Issue or check the denominator-cleared cubic cocycle identity.
//!
//! For nonzero a, b, c in a commutative field, the identity
//! (b/a)(c/b)(a/c)=1 clears denominators to (bc)a=(ab)c.
//! In the cubic norm-one construction b=σ(a) and c=σ²(a).
//!
//! `cargo run --release --example weak_curve_cubic_identity -- issue` emits
//! the certificate. `... -- check path.json` replays it with fresh points.

use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr, IdentityCertificate, IdentityStatement, Verdict,
};
use std::{env, fs};

fn statement() -> IdentityStatement {
    IdentityStatement {
        name: "Cubic cocycle denominator-cleared product".into(),
        source: "research/weak_curve_database_20261012/lean/WeakCurveDatabase.lean; forward Hilbert-90 identity".into(),
        variables: vec!["a".into(), "b".into(), "c".into()],
        degree_bound: 3,
        lhs: Expr::product(vec![Expr::var("b"), Expr::var("c"), Expr::var("a")]),
        rhs: Expr::product(vec![Expr::var("a"), Expr::var("b"), Expr::var("c")]),
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    match args.as_slice() {
        [_, mode] if mode == "issue" => {
            let cert = issue(
                statement(),
                0x2026_1012_0901,
                32,
                "issued by examples/weak_curve_cubic_identity.rs using the native certificate module; independently replayable from the sealed statement id and seed",
                None,
            )
            .expect("polynomial identity must hold");
            println!("{}", serde_json::to_string_pretty(&cert).unwrap());
        }
        [_, mode, path] if mode == "check" => {
            let cert: IdentityCertificate =
                serde_json::from_slice(&fs::read(path).expect("certificate file")).unwrap();
            assert_eq!(cert.statement, statement(), "wrong identity statement");
            match check_with_extra(&cert, Some((0x71ac_2026_1012, 32))).unwrap() {
                Verdict::Accept { .. } => println!("CHECKED\t{}", cert.id),
                other => panic!("identity certificate rejected: {other:?}"),
            }
        }
        _ => panic!("usage: weak_curve_cubic_identity issue | check CERTIFICATE.json"),
    }
}
