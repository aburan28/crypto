//! Replay the polynomial identities used in the ISO-1 obstruction proof.
//! These records cover polynomial algebra, not finite-field census labels.

use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr, IdentityCertificate, IdentityStatement, Verdict,
};

fn square(value: Expr) -> Expr {
    Expr::times(value.clone(), value)
}

fn statements() -> Vec<IdentityStatement> {
    let x = Expr::var("x");
    let s = Expr::var("s");
    let a = Expr::var("a");
    let b = Expr::var("b");
    let u = Expr::sum(vec![
        square(x.clone()),
        Expr::times(a.clone(), x.clone()),
        b.clone(),
    ]);
    vec![
        IdentityStatement {
            name: "ISO-1 quartic-parameter quotient factorization".into(),
            source: "research/iso1_weak_classes_20261007/REPORT.md; lambda=s^2, s=mu^2".into(),
            variables: vec!["x".into(), "s".into()],
            degree_bound: 4,
            lhs: Expr::sum(vec![
                square(x.clone()),
                Expr::product(vec![
                    Expr::Const(2),
                    Expr::plus(Expr::Const(1), square(s.clone())),
                    x.clone(),
                ]),
                square(Expr::plus(Expr::Const(1), Expr::minus(square(s.clone())))),
            ]),
            rhs: Expr::times(
                Expr::plus(x.clone(), square(Expr::plus(Expr::Const(1), s.clone()))),
                Expr::plus(
                    x.clone(),
                    square(Expr::plus(Expr::Const(1), Expr::minus(s))),
                ),
            ),
        },
        IdentityStatement {
            name: "ISO-1 rational 2-isogeny cleared equation".into(),
            source: "research/iso1_weak_classes_20261007/REPORT.md; U=x^2+a*x+b".into(),
            variables: vec!["x".into(), "a".into(), "b".into()],
            degree_bound: 4,
            lhs: Expr::sum(vec![
                square(u.clone()),
                Expr::product(vec![Expr::Const(-2), a.clone(), x.clone(), u]),
                Expr::times(
                    Expr::plus(square(a), Expr::times(Expr::Const(-4), b.clone())),
                    square(x.clone()),
                ),
            ]),
            rhs: square(Expr::plus(square(x), Expr::minus(b))),
        },
    ]
}

fn replay(certificates: &[IdentityCertificate]) {
    let expected = statements();
    assert_eq!(certificates.len(), expected.len(), "identity record count");
    for (certificate, statement) in certificates.iter().zip(expected) {
        assert_eq!(certificate.statement, statement, "proof statement changed");
        assert!(matches!(
            check_with_extra(certificate, Some((2026100817, 64))).expect("valid record"),
            Verdict::Accept {
                points: 32,
                extra_points: 64,
                ..
            }
        ));
        eprintln!(
            "{}: 32 recorded + 64 additional evaluations accepted",
            certificate.id
        );
    }
}

fn main() {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    match arguments.as_slice() {
        [] => {
            let certificates: Vec<_> = statements().into_iter().map(|statement| {
                issue(statement, 20261008, 32,
                    "Native issuer and replay use the same repository checker; exact integer-polynomial expansion is checked separately by PARI/GP. No formal theorem-prover artifact is asserted.",
                    None).expect("identity holds")
            }).collect();
            replay(&certificates);
            println!(
                "{}",
                serde_json::to_string_pretty(&certificates).expect("serialize records")
            );
        }
        [path] => {
            let bytes = std::fs::read(path).expect("read certificates");
            let certificates: Vec<IdentityCertificate> =
                serde_json::from_slice(&bytes).expect("parse certificates");
            replay(&certificates);
        }
        _ => panic!("usage: iso1_report_check [identity_certificates.json]"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn committed_identities_replay() {
        let manifest = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let root = if manifest
            .join("src/cryptanalysis/identity_certificate.rs")
            .exists()
        {
            manifest
        } else {
            manifest.ancestors().nth(3).expect("repository root")
        };
        let path = root.join("research/iso1_weak_classes_20261007/identity_certificates.json");
        let certificates: Vec<IdentityCertificate> =
            serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
        replay(&certificates);
    }

    #[test]
    fn changed_polynomial_is_refused() {
        for mut statement in statements() {
            statement.rhs = Expr::plus(statement.rhs, Expr::Const(1));
            assert!(issue(statement, 20261008, 32, "mutation control", None).is_err());
        }
    }
}
