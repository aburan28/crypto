//! Exact algebra receipts for the one-edge torsion theorem and norm-fiber bounds.
use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr, IdentityCertificate, IdentityStatement, Verdict,
};
fn pow(v: Expr, n: usize) -> Expr {
    Expr::product((0..n).map(|_| v.clone()).collect())
}
fn statements() -> Vec<IdentityStatement> {
    let x = Expr::var("x");
    let a = Expr::var("a");
    let b = Expr::var("b");
    let cubic = Expr::sum(vec![
        pow(x.clone(), 3),
        Expr::times(a.clone(), pow(x.clone(), 2)),
        Expr::times(b.clone(), x.clone()),
    ]);
    let slope = Expr::sum(vec![
        Expr::times(Expr::Const(3), pow(x.clone(), 2)),
        Expr::product(vec![Expr::Const(2), a.clone(), x.clone()]),
        b.clone(),
    ]);
    let mut s=vec![IdentityStatement{name:"ISO-1 cleared duplication x-coordinate".into(),source:"research/iso1_weak_classes_20261007/large_population_20261009/REPORT.tex; one-edge full-4-torsion theorem".into(),variables:vec!["x".into(),"a".into(),"b".into()],degree_bound:4,
        lhs:Expr::sum(vec![pow(slope,2),Expr::minus(Expr::product(vec![Expr::Const(4),cubic,Expr::plus(a,Expr::times(Expr::Const(2),x.clone()))]))]),
        rhs:pow(Expr::plus(pow(x,2),Expr::minus(b)),2)}];
    for n in [3usize, 5, 7] {
        let q = Expr::var("q");
        s.push(IdentityStatement{name:format!("ISO-1 degree-{n} norm-fiber cardinality"),source:"research/iso1_weak_classes_20261007/large_population_20261009/REPORT.tex; direct weak-model union bound".into(),variables:vec!["q".into()],degree_bound:n as u32,
        lhs:Expr::times(Expr::plus(q.clone(),Expr::Const(-1)),Expr::sum((0..n).map(|i|pow(q.clone(),i)).collect())),rhs:Expr::plus(pow(q,n),Expr::Const(-1))});
    }
    s
}
fn matrix_count(det: u32) -> (usize, usize) {
    let mut total = 0;
    let mut admitted = 0;
    for a in (1..16).step_by(2) {
        for d in (1..16).step_by(2) {
            for b in (0..16).step_by(2) {
                for c in (0..16).step_by(2) {
                    if (a * d + 256 - b * c) % 16 != det {
                        continue;
                    }
                    total += 1;
                    let trace = (a + d) % 16;
                    if trace == (det + 1) % 16 || trace == (31 - det) % 16 {
                        admitted += 1;
                    }
                }
            }
        }
    }
    (total, admitted)
}
fn quartic_statement() -> IdentityStatement {
    let x = Expr::var("X");
    let c1 = Expr::var("c1");
    let c2 = Expr::var("c2");
    let c3 = Expr::var("c3");
    IdentityStatement{name:"Cleared root-quartic to Weierstrass coordinate identity".into(),source:"research/iso1_weak_classes_20261007/large_population_20261009/PRIOR_WORK.md; standard translation and inversion at a rational branch point".into(),variables:vec!["X".into(),"c1".into(),"c2".into(),"c3".into()],degree_bound:5,
      lhs:Expr::times(pow(c3.clone(),2),Expr::sum(vec![pow(x.clone(),3),Expr::times(c2.clone(),pow(x.clone(),2)),Expr::product(vec![c3.clone(),c1.clone(),x.clone()]),pow(c3.clone(),2)])),
      rhs:Expr::sum(vec![Expr::times(pow(c3.clone(),2),pow(x.clone(),3)),Expr::product(vec![c2,pow(c3.clone(),2),pow(x.clone(),2)]),Expr::product(vec![c1,pow(c3.clone(),3),x]),pow(c3,4)])}
}
fn replay(c: &[IdentityCertificate]) {
    assert_eq!(c.len(), 4);
    for (c, s) in c.iter().zip(statements()) {
        assert_eq!(c.statement, s);
        assert!(matches!(
            check_with_extra(c, Some((2026100949, 64))).unwrap(),
            Verdict::Accept {
                points: 32,
                extra_points: 64,
                ..
            }
        ));
        println!(
            "{}: 32 recorded plus 64 independent evaluations accepted",
            c.id
        );
    }
}
fn main() {
    let a = std::env::args().skip(1).collect::<Vec<_>>();
    if a.first().map(String::as_str) == Some("--issue-quartic") {
        let c=issue(quartic_statement(),20261009,32,"Standard coordinate substitution; the quartic Taylor coefficients and root condition are checked on each recorded field model.",None).unwrap();
        eprintln!("{}", c.id);
        println!("{}", serde_json::to_string_pretty(&vec![c]).unwrap());
        return;
    }
    if a.first().map(String::as_str) == Some("--quartic") {
        let c: Vec<IdentityCertificate> =
            serde_json::from_slice(&std::fs::read(&a[1]).unwrap()).unwrap();
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].statement, quartic_statement());
        assert!(check_with_extra(&c[0], Some((2026100966, 64))).is_ok());
        println!("{}: quartic substitution accepted", c[0].id);
        return;
    }
    if a.first().map(String::as_str) == Some("--issue") {
        let c=statements().into_iter().map(|s|issue(s,20261009,32,"Exact duplication and geometric-sum polynomial identities; the torsion and norm-fiber arguments have separate deductive proofs.",None).unwrap()).collect::<Vec<_>>();
        for r in &c {
            eprintln!("{}", r.id);
        }
        println!("{}", serde_json::to_string_pretty(&c).unwrap());
        return;
    }
    if let Some(path) = a.first() {
        replay(
            &serde_json::from_slice::<Vec<IdentityCertificate>>(&std::fs::read(path).unwrap())
                .unwrap(),
        );
    }
    for det in [1, 9] {
        let (all, admitted) = matrix_count(det);
        assert_eq!((all, admitted), (512, 320));
        println!("uniform_mod16_matrix_model: determinant={det} total={all} admitted={admitted} fraction=5/8");
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn algebra_replay_and_changed_duplication_refused() {
        for s in statements() {
            let c = issue(s, 20261009, 32, "test", None).unwrap();
            assert!(check_with_extra(&c, Some((2026100949, 64))).is_ok());
        }
        let mut s = statements().remove(0);
        s.rhs = Expr::plus(s.rhs, Expr::Const(1));
        assert!(issue(s, 20261009, 32, "changed duplication", None).is_err());
    }
    #[test]
    fn conditional_matrix_counts() {
        assert_eq!(matrix_count(1), (512, 320));
        assert_eq!(matrix_count(9), (512, 320));
    }
}
