//! Replay the halving and orbit arithmetic accompanying THEOREM.tex.
use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr, IdentityCertificate, IdentityStatement, Verdict,
};

fn square(v: Expr) -> Expr {
    Expr::times(v.clone(), v)
}

fn statements() -> Vec<IdentityStatement> {
    let a = Expr::var("a");
    let b = Expr::var("b");
    let ab = Expr::times(a.clone(), b.clone());
    let p = Expr::var("p");
    let p2 = square(p);
    let p4 = square(p2.clone());
    vec![
        IdentityStatement {
            name:"ISO-1 rational halving point lies on the translated curve".into(),
            source:"research/iso1_weak_classes_20261007/THEOREM.tex; P=(a*b,a*b*(a+b)) after translating a root to zero".into(),
            variables:vec!["a".into(),"b".into()],degree_bound:6,
            lhs:Expr::times(square(ab.clone()),square(Expr::plus(a.clone(),b.clone()))),
            rhs:Expr::product(vec![ab.clone(),Expr::plus(ab.clone(),square(a)),Expr::plus(ab,square(b))]),
        },
        IdentityStatement {
            name:"ISO-1 absolute-orbit counting numerator".into(),
            source:"research/iso1_weak_classes_20261007/THEOREM.tex; twelve times the total orbit count".into(),
            variables:vec!["p".into()],degree_bound:4,
            lhs:Expr::sum(vec![Expr::plus(p4.clone(),Expr::minus(p2.clone())),Expr::times(Expr::Const(4),Expr::plus(p2.clone(),Expr::Const(-1))),Expr::Const(12)]),
            rhs:Expr::sum(vec![p4,Expr::times(Expr::Const(3),p2),Expr::Const(8)]),
        },
    ]
}

fn replay(records: &[IdentityCertificate]) {
    let expected = statements();
    assert_eq!(records.len(), expected.len());
    for (record, statement) in records.iter().zip(expected) {
        assert_eq!(record.statement, statement);
        assert!(matches!(
            check_with_extra(record, Some((2026100917, 64))).unwrap(),
            Verdict::Accept {
                points: 32,
                extra_points: 64,
                ..
            }
        ));
        eprintln!(
            "{}: 32 recorded + 64 additional evaluations accepted",
            record.id
        );
    }
}

fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.is_empty() {
        let records: Vec<_>=statements().into_iter().map(|s|issue(s,20261009,32,"Native polynomial replay; separate exact PARI expansion. The finite-field and order-theoretic arguments are written proofs.",None).unwrap()).collect();
        replay(&records);
        println!("{}", serde_json::to_string_pretty(&records).unwrap());
    } else {
        assert_eq!(
            args.len(),
            1,
            "usage: iso1_theorem_check [theorem_identity_certificates.json]"
        );
        replay(
            &serde_json::from_slice::<Vec<IdentityCertificate>>(&std::fs::read(&args[0]).unwrap())
                .unwrap(),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn retained_theorem_records_replay() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(3)
            .unwrap();
        let path =
            root.join("research/iso1_weak_classes_20261007/theorem_identity_certificates.json");
        replay(
            &serde_json::from_slice::<Vec<IdentityCertificate>>(&std::fs::read(path).unwrap())
                .unwrap(),
        );
    }
    #[test]
    fn changed_theorem_polynomial_is_refused() {
        for mut s in statements() {
            s.rhs = Expr::plus(s.rhs, Expr::Const(1));
            assert!(issue(s, 20261009, 32, "mutation control", None).is_err());
        }
    }
}
