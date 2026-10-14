//! Issue or replay the two denominator-cleared polynomial steps behind the
//! cubic torsion-equivariant coordinate. The elliptic-curve group law and
//! exceptional projective points are checked separately in the study.

use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr as E, IdentityCertificate, IdentityStatement, Verdict,
};

fn c(n: i64) -> E {
    E::Const(n)
}

fn sub(a: E, b: E) -> E {
    E::plus(a, E::minus(b))
}

fn pow(a: E, n: usize) -> E {
    E::product((0..n).map(|_| a.clone()).collect())
}

fn statement(kind: &str) -> IdentityStatement {
    let x = E::var("x");
    let y = E::var("y");
    let ym1 = sub(y.clone(), c(1));
    let yp1 = E::plus(y.clone(), c(1));
    let ym2 = sub(y.clone(), c(2));
    let ym3 = sub(y.clone(), c(3));
    match kind {
        // On y²=x³+1, the generic addition formula for T=(0,1) yields
        // y(P+T) = (y-3)/(y+1). Here its denominator-cleared difference is
        // reduced by the curve equation. Exceptions are not in this record.
        "group_law" => {
            let x3 = pow(x.clone(), 3);
            let curve = sub(sub(pow(y.clone(), 2), x3.clone()), c(1));
            let lhs = E::sum(vec![
                E::product(vec![x3.clone(), yp1.clone(), ym2]),
                E::minus(E::product(vec![pow(ym1.clone(), 3), yp1])),
                E::minus(E::times(x3, ym3)),
                E::times(pow(ym1, 2), curve),
            ]);
            IdentityStatement {
                name: "Cubic torsion translation y-coordinate reduction".into(),
                source: "research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/REPORT.md"
                    .into(),
                variables: vec!["x".into(), "y".into()],
                degree_bound: 5,
                lhs,
                rhs: c(0),
            }
        }
        // Conjugate g(y)=(y-3)/(y+1) by chi(y)=(y+1)/(y-1):
        // chi(g(y)) = 1/(1-chi(y)), after clearing denominators.
        "conjugacy" => IdentityStatement {
            name: "Cubic torsion Mobius conjugacy".into(),
            source: "research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/REPORT.md"
                .into(),
            variables: vec!["y".into()],
            degree_bound: 2,
            lhs: E::times(
                E::plus(ym3.clone(), yp1.clone()),
                sub(ym1.clone(), yp1.clone()),
            ),
            rhs: E::times(sub(ym3, yp1), ym1),
        },
        other => panic!("unknown identity kind: {other}"),
    }
}

fn accepted(cert: &IdentityCertificate) -> Verdict {
    let verdict = check_with_extra(cert, Some((0x2026_1014, 64))).expect("well-formed record");
    assert!(matches!(verdict, Verdict::Accept { .. }), "{verdict:?}");
    verdict
}

fn main() {
    let mut args = std::env::args().skip(1);
    match (args.next().as_deref(), args.next()) {
        (Some("issue"), Some(kind)) => {
            let cert = issue(
                statement(&kind),
                0x2026_1013,
                64,
                "issued by examples/isogeny_torsion_identity.rs; the same native checker replays the id-bound points and a separate checker-seed stream",
                None,
            )
            .expect("polynomial identity must hold");
            accepted(&cert);
            println!("{}", serde_json::to_string_pretty(&cert).unwrap());
        }
        (Some("verify"), Some(path)) => {
            let bytes = std::fs::read(&path).expect("certificate path");
            let cert: IdentityCertificate =
                serde_json::from_slice(&bytes).expect("certificate JSON");
            let verdict = accepted(&cert);
            println!("{} {verdict:?}", cert.id);
        }
        _ => panic!("usage: isogeny_torsion_identity issue <group_law|conjugacy> | verify <path>"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_steps_issue_and_replay() {
        for kind in ["group_law", "conjugacy"] {
            let cert = issue(statement(kind), 0x2026_1013, 64, "test", None).unwrap();
            accepted(&cert);
        }
    }

    #[test]
    fn corrupted_group_law_is_rejected() {
        let mut bad = statement("group_law");
        bad.rhs = c(1);
        assert!(issue(bad, 0x2026_1013, 64, "mutation", None).is_err());
    }
}
