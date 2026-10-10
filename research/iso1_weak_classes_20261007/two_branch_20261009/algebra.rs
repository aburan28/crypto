//! Portable identity certificates for the two-branch parametrization.
use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr, IdentityCertificate, IdentityStatement,
};
fn pow(x: Expr, n: usize) -> Expr {
    Expr::product((0..n).map(|_| x.clone()).collect())
}
fn sum(v: Vec<Expr>) -> Expr {
    Expr::sum(v)
}
fn mul(v: Vec<Expr>) -> Expr {
    Expr::product(v)
}
fn statements() -> Vec<IdentityStatement> {
    let q = Expr::var("q");
    let u = Expr::var("u");
    let v = Expr::var("v");
    let b = Expr::var("b");
    let l = Expr::var("l");
    let am = sum(vec![
        mul(vec![b.clone(), sum(vec![u.clone(), Expr::Const(1)])]),
        Expr::minus(mul(vec![b.clone(), sum(vec![u.clone(), Expr::Const(-1)])])),
    ]);
    let ap = sum(vec![
        mul(vec![b.clone(), sum(vec![u.clone(), Expr::Const(1)])]),
        mul(vec![b.clone(), sum(vec![u.clone(), Expr::Const(-1)])]),
    ]);
    let bp = sum(vec![
        Expr::minus(mul(vec![b.clone(), sum(vec![v.clone(), Expr::Const(1)])])),
        mul(vec![b.clone(), sum(vec![v.clone(), Expr::Const(-1)])]),
    ]);
    let bm = sum(vec![
        Expr::minus(mul(vec![b.clone(), sum(vec![v.clone(), Expr::Const(1)])])),
        Expr::minus(mul(vec![b.clone(), sum(vec![v.clone(), Expr::Const(-1)])])),
    ]);
    let mut out = vec![
        IdentityStatement {
            name: "Alternating torus cardinality".into(),
            source: "two_branch_20261009/REPORT.tex; normalization theorem".into(),
            variables: vec!["q".into()],
            degree_bound: 3,
            lhs: mul(vec![
                sum(vec![q.clone(), Expr::Const(1)]),
                sum(vec![
                    pow(q.clone(), 2),
                    Expr::minus(q.clone()),
                    Expr::Const(1),
                ]),
            ]),
            rhs: sum(vec![pow(q, 3), Expr::Const(1)]),
        },
        IdentityStatement {
            name: "Nonsplit-pair cleared cross-ratio".into(),
            source: "two_branch_20261009/REPORT.tex; alpha=b(u+1)/(u-1), alpha_q=-b(v+1)/(v-1)"
                .into(),
            variables: vec!["u".into(), "v".into(), "b".into()],
            degree_bound: 6,
            lhs: mul(vec![u, v, am, bp]),
            rhs: mul(vec![ap, bm]),
        },
    ];
    let alpha = Expr::var("a");
    let conjugate = Expr::var("b");
    let d = Expr::var("d");
    let s = sum(vec![d.clone(), Expr::minus(pow(alpha.clone(), 2))]);
    let c2 = sum(vec![
        Expr::times(Expr::Const(3), pow(alpha.clone(), 2)),
        mul(vec![Expr::Const(-2), alpha.clone(), conjugate.clone()]),
        Expr::minus(d.clone()),
    ]);
    let c3 = mul(vec![
        sum(vec![pow(alpha.clone(), 2), Expr::minus(d.clone())]),
        sum(vec![alpha.clone(), Expr::minus(conjugate.clone())]),
    ]);
    let c1 = sum(vec![
        Expr::times(Expr::Const(3), alpha.clone()),
        Expr::minus(conjugate.clone()),
    ]);
    out.push(IdentityStatement {
        name: "Quadratic-branch rational 2-kernel square".into(),
        source: "two_branch_20261009/REPORT.tex; translated 2-kernel b=(a^2-d)(b^2-d)".into(),
        variables: vec!["a".into(), "b".into(), "d".into()],
        degree_bound: 4,
        lhs: sum(vec![
            Expr::times(Expr::Const(3), pow(s.clone(), 2)),
            mul(vec![Expr::Const(2), c2, s.clone()]),
            mul(vec![c3, c1]),
        ]),
        rhs: mul(vec![s, sum(vec![d, Expr::minus(pow(conjugate, 2))])]),
    });
    let ln = sum(vec![
        pow(l.clone(), 2),
        Expr::minus(l.clone()),
        Expr::Const(1),
    ]);
    out.push(IdentityStatement {
        name: "Legendre reciprocal j numerator and denominator".into(),
        source: "two_branch_20261009/REPORT.tex; j=256(Y-1)^3/(Y-2), Y=l+l^-1".into(),
        variables: vec!["l".into()],
        degree_bound: 10,
        lhs: mul(vec![
            pow(ln.clone(), 3),
            pow(l.clone(), 2),
            sum(vec![
                pow(l.clone(), 2),
                Expr::times(Expr::Const(-2), l.clone()),
                Expr::Const(1),
            ]),
        ]),
        rhs: mul(vec![
            pow(
                sum(vec![
                    pow(l.clone(), 2),
                    Expr::Const(1),
                    Expr::minus(l.clone()),
                ]),
                3,
            ),
            pow(l.clone(), 2),
            pow(sum(vec![l, Expr::Const(-1)]), 2),
        ]),
    });
    let a = Expr::var("a");
    let b = Expr::var("b");
    let beta = Expr::var("beta");
    out.push(IdentityStatement {
        name: "Quadratic-base norm equation".into(),
        source: "two_branch_20261009/REPORT.tex; generalized Hilbert-90 normalization, beta^2=d"
            .into(),
        variables: vec!["a".into(), "b".into(), "beta".into()],
        degree_bound: 4,
        lhs: mul(vec![
            sum(vec![a.clone(), mul(vec![b.clone(), beta.clone()])]),
            sum(vec![
                a.clone(),
                Expr::minus(mul(vec![b.clone(), beta.clone()])),
            ]),
        ]),
        rhs: sum(vec![
            pow(a, 2),
            Expr::minus(mul(vec![pow(b, 2), pow(beta, 2)])),
        ]),
    });
    let degree = Expr::var("degree");
    let n = Expr::var("n");
    out.push(IdentityStatement{name:"Cleared compositum-cover genus formula".into(),source:"two_branch_20261009/REPORT.tex; Riemann-Hurwitz with degree=2^n and n+2 index-two branch points".into(),variables:vec!["degree".into(),"n".into()],degree_bound:2,lhs:sum(vec![mul(vec![degree.clone(),sum(vec![n.clone(),Expr::Const(2)])]),Expr::times(Expr::Const(-4),degree.clone())]),rhs:mul(vec![degree,sum(vec![n,Expr::Const(-2)])])});
    let l = Expr::var("l");
    out.push(IdentityStatement{name:"Cleared norm of complementary cross-ratio".into(),source:"two_branch_20261009/REPORT.tex; exact degree-six quadratic invariant, lambda^Q=lambda^-1".into(),variables:vec!["l".into()],degree_bound:2,lhs:mul(vec![sum(vec![Expr::Const(1),Expr::minus(l.clone())]),sum(vec![l.clone(),Expr::Const(-1)])]),rhs:sum(vec![Expr::times(Expr::Const(2),l.clone()),Expr::minus(pow(l,2)),Expr::Const(-1)])});
    out
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    if a.first().map(String::as_str) == Some("--issue") {
        let c:Vec<_>=statements().into_iter().map(|s|issue(s,2026100910,32,"Cleared polynomial identities; field descent, completeness and orbit length have deductive proofs.",None).unwrap()).collect();
        for v in &c {
            eprintln!("{}", v.id);
        }
        println!("{}", serde_json::to_string_pretty(&c).unwrap());
        return;
    }
    assert_eq!(a.len(), 1);
    let c: Vec<IdentityCertificate> =
        serde_json::from_slice(&std::fs::read(&a[0]).unwrap()).unwrap();
    assert_eq!(c.len(), 7);
    for (c, s) in c.iter().zip(statements()) {
        assert_eq!(c.statement, s);
        assert!(check_with_extra(c, Some((2026100919, 64))).is_ok());
        println!(
            "{} accepted with 32 saved and 64 independent evaluations",
            c.id
        );
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn accept_and_mutation_refusal() {
        for s in statements() {
            let c = issue(s.clone(), 2026100910, 32, "test", None).unwrap();
            assert!(check_with_extra(&c, Some((2026100919, 64))).is_ok());
            let mut wrong = s;
            wrong.rhs = Expr::plus(wrong.rhs, Expr::Const(1));
            assert!(issue(wrong, 2026100910, 32, "mutation", None).is_err());
        }
    }
}
