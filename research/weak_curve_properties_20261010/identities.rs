//! Portable evaluations and exact integer-coefficient checks for cleared identities.
use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr, IdentityCertificate, IdentityStatement,
};
use serde_json::json;
use std::{collections::BTreeMap, fs, path::Path};
fn pow(e: Expr, n: usize) -> Expr {
    Expr::product((0..n).map(|_| e.clone()).collect())
}
fn mul(v: Vec<Expr>) -> Expr {
    Expr::product(v)
}
fn add(v: Vec<Expr>) -> Expr {
    Expr::sum(v)
}
fn statements() -> Vec<IdentityStatement> {
    let r = Expr::var("r");
    let s = Expr::var("s");
    let a = Expr::minus(add(vec![
        pow(r.clone(), 2),
        mul(vec![r.clone(), s.clone()]),
        pow(s.clone(), 2),
    ]));
    let b = mul(vec![r.clone(), s.clone(), add(vec![r.clone(), s.clone()])]);
    let mut out = vec![IdentityStatement {
        name: "Cubic discriminant from three distinct zero-sum roots".into(),
        source: "weak_curve_properties_20261010/REPORT.md; torsion criterion".into(),
        variables: vec!["r".into(), "s".into()],
        degree_bound: 6,
        lhs: mul(vec![
            pow(add(vec![r.clone(), Expr::minus(s.clone())]), 2),
            pow(
                add(vec![Expr::times(Expr::Const(2), r.clone()), s.clone()]),
                2,
            ),
            pow(
                add(vec![r.clone(), Expr::times(Expr::Const(2), s.clone())]),
                2,
            ),
        ]),
        rhs: add(vec![
            Expr::times(Expr::Const(-4), pow(a, 3)),
            Expr::times(Expr::Const(-27), pow(b, 2)),
        ]),
    }];
    let a = Expr::var("a");
    let b = Expr::var("b");
    let x = Expr::var("x");
    let y = Expr::var("y");
    let c = Expr::var("c");
    let numerator = Expr::times(Expr::Const(6912), pow(a.clone(), 3));
    let denominator = add(vec![
        Expr::times(Expr::Const(4), pow(a.clone(), 3)),
        Expr::times(Expr::Const(27), pow(b.clone(), 2)),
    ]);
    out.push(IdentityStatement {
        name: "Cleared normalized j coefficient".into(),
        source: "weak_curve_properties_20261010/REPORT.md; generic twist normalization".into(),
        variables: vec!["a".into(), "b".into()],
        degree_bound: 6,
        lhs: mul(vec![numerator.clone(), Expr::Const(27), pow(b.clone(), 2)]),
        rhs: mul(vec![
            Expr::Const(4),
            pow(a.clone(), 3),
            add(vec![
                Expr::times(Expr::Const(1728), denominator),
                Expr::minus(numerator),
            ]),
        ]),
    });
    out.push(IdentityStatement {
        name: "Cleared quadratic-twist a coefficient".into(),
        source: "weak_curve_properties_20261010/REPORT.md; D=3b/(2a)".into(),
        variables: vec!["a".into(), "b".into()],
        degree_bound: 5,
        lhs: mul(vec![
            pow(Expr::times(Expr::Const(3), b.clone()), 2),
            Expr::Const(4),
            pow(a.clone(), 3),
        ]),
        rhs: mul(vec![
            a.clone(),
            pow(Expr::times(Expr::Const(2), a.clone()), 2),
            Expr::Const(9),
            pow(b.clone(), 2),
        ]),
    });
    out.push(IdentityStatement {
        name: "Cleared quadratic-twist b coefficient".into(),
        source: "weak_curve_properties_20261010/REPORT.md; D=3b/(2a)".into(),
        variables: vec!["a".into(), "b".into()],
        degree_bound: 6,
        lhs: mul(vec![
            pow(Expr::times(Expr::Const(3), b.clone()), 3),
            Expr::Const(8),
            pow(a.clone(), 3),
        ]),
        rhs: mul(vec![
            b.clone(),
            pow(Expr::times(Expr::Const(2), a), 3),
            Expr::Const(27),
            pow(b, 2),
        ]),
    });
    out.push(IdentityStatement {
        name: "Translated Frobenius norm".into(),
        source: "weak_curve_properties_20261010/REPORT.md; independent integral-basis scan".into(),
        variables: vec!["x".into(), "y".into(), "c".into()],
        degree_bound: 2,
        lhs: mul(vec![
            add(vec![x.clone(), Expr::minus(c.clone())]),
            add(vec![y.clone(), Expr::minus(c.clone())]),
        ]),
        rhs: add(vec![
            mul(vec![x.clone(), y.clone()]),
            Expr::minus(mul(vec![c.clone(), add(vec![x, y])])),
            pow(c, 2),
        ]),
    });
    for n in [1usize, 2, 3, 5, 6, 7, 10, 14] {
        let q = Expr::var("q");
        out.push(IdentityStatement {
            name: format!("Degree-{n} unit-group index"),
            source: "weak_curve_properties_20261010/REPORT.md; finite cyclic twist descent".into(),
            variables: vec!["q".into()],
            degree_bound: n as u32,
            lhs: mul(vec![
                add(vec![q.clone(), Expr::Const(-1)]),
                add((0..n).map(|i| pow(q.clone(), i)).collect()),
            ]),
            rhs: add(vec![pow(q, n), Expr::Const(-1)]),
        });
    }
    out
}
type Poly = BTreeMap<Vec<u32>, i128>;
fn expand(e: &Expr, vars: &[String]) -> Poly {
    let constant = |c| {
        let mut p = Poly::new();
        if c != 0 {
            p.insert(vec![0; vars.len()], c);
        }
        p
    };
    match e {
        Expr::Const(c) => constant(*c as i128),
        Expr::Var(v) => {
            let mut exponent = vec![0; vars.len()];
            exponent[vars.iter().position(|x| x == v).unwrap()] = 1;
            [(exponent, 1)].into_iter().collect()
        }
        Expr::Neg(e) => expand(e, vars).into_iter().map(|(m, c)| (m, -c)).collect(),
        Expr::Add(v) => {
            let mut out = Poly::new();
            for e in v {
                for (m, c) in expand(e, vars) {
                    *out.entry(m).or_default() += c;
                }
            }
            out.retain(|_, c| *c != 0);
            out
        }
        Expr::Mul(v) => {
            let mut out = constant(1);
            for e in v {
                let next = expand(e, vars);
                let mut p = Poly::new();
                for (m, c) in &out {
                    for (n, d) in &next {
                        let exponent = m.iter().zip(n).map(|(a, b)| a + b).collect();
                        *p.entry(exponent).or_default() += c * d;
                    }
                }
                p.retain(|_, c| *c != 0);
                out = p;
            }
            out
        }
    }
}
pub fn issue_all(base: &Path) {
    let mut certs = vec![];
    let mut exact = vec![];
    for s in statements() {
        let lhs = expand(&s.lhs, &s.variables);
        let rhs = expand(&s.rhs, &s.variables);
        assert_eq!(lhs, rhs);
        let c=issue(s,202610102601,64,"Native exact integer expansion plus portable evaluations; field and integrality hypotheses have deductive proofs.",None).unwrap();
        exact.push(json!({"id":c.id,"terms":lhs.iter().map(|(e,c)|json!({"exponents":e,"coefficient":c.to_string()})).collect::<Vec<_>>(),"status":"exact_integer_coefficients_equal"}));
        certs.push(c);
    }
    fs::write(
        base.join("identity_certificates.json"),
        serde_json::to_vec_pretty(&certs).unwrap(),
    )
    .unwrap();
    fs::write(
        base.join("identity_exact.json"),
        serde_json::to_vec_pretty(&exact).unwrap(),
    )
    .unwrap();
}
pub fn check_all(base: &Path) {
    let certs: Vec<IdentityCertificate> =
        serde_json::from_slice(&fs::read(base.join("identity_certificates.json")).unwrap())
            .unwrap();
    let exact: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("identity_exact.json")).unwrap()).unwrap();
    assert_eq!(certs.len(), statements().len());
    for ((c, s), record) in certs
        .iter()
        .zip(statements())
        .zip(exact.as_array().unwrap())
    {
        assert_eq!(c.statement, s);
        let p = expand(&s.lhs, &s.variables);
        assert_eq!(p, expand(&s.rhs, &s.variables));
        assert_eq!(
            record["terms"],
            json!(p
                .iter()
                .map(|(e, c)| json!({"exponents":e,"coefficient":c.to_string()}))
                .collect::<Vec<_>>())
        );
        assert_eq!(record["id"], c.id);
        check_with_extra(c, Some((202610102609, 128))).unwrap();
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn identities_and_mutation() {
        for s in statements() {
            assert_eq!(expand(&s.lhs, &s.variables), expand(&s.rhs, &s.variables));
            let mut bad = s.clone();
            bad.rhs = Expr::plus(bad.rhs, Expr::Const(1));
            assert_ne!(
                expand(&bad.lhs, &bad.variables),
                expand(&bad.rhs, &bad.variables)
            );
            assert!(issue(bad, 202610102601, 64, "mutation", None).is_err());
        }
    }
}
