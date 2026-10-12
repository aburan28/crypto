//! Prime odd-degree Frobenius/inversion orbit formula and polynomial certificate.
use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr, IdentityCertificate, IdentityStatement, Verdict,
};
fn values(p: u128, n: u32) -> (u128, u128, u128, u128, u128) {
    let big = (p.pow(2 * n) - 1) / (p * p - 1);
    let a = (p.pow(n) - 1) / (p - 1);
    let b = (p.pow(n) + 1) / (p + 1);
    let e = u128::from((p * p) % u128::from(n) == 1);
    let numerator = big + a + b - 3 + 2 * e * u128::from(n - 1).pow(2);
    assert_eq!(numerator % (4 * u128::from(n)), 0);
    (big, a, b, e, numerator / (4 * u128::from(n)))
}
fn orbit(x: u128, p: u128, n: u32, big: u128) -> Vec<u128> {
    let mut items = Vec::with_capacity(4 * n as usize);
    let mut y = x;
    for _ in 0..2 * n {
        items.push(y);
        items.push((big - y) % big);
        y = y * p % big;
    }
    items.sort_unstable();
    items.dedup();
    items
}
fn exhaustive(p: u128, n: u32) -> (usize, usize, usize) {
    let (big, a, b, e, total) = values(p, n);
    assert!(big < 6_000_000);
    let mut visited = vec![false; big as usize];
    let mut counts = [0usize; 3];
    for x in 1..big {
        if visited[x as usize] {
            continue;
        }
        let o = orbit(x, p, n, big);
        for &v in &o {
            visited[v as usize] = true;
        }
        if o.len() == 2 {
            counts[0] += 1;
        } else if o.len() == 2 * n as usize {
            counts[1] += 1;
        } else {
            assert_eq!(o.len(), 4 * n as usize);
            counts[2] += 1;
        }
    }
    assert_eq!(counts[0] as u128, e * u128::from(n - 1) / 2);
    assert_eq!(
        counts[1] as u128,
        (a + b - 2 - e * u128::from(n - 1)) / (2 * u128::from(n))
    );
    assert_eq!(counts[2] as u128, (big + 1 - a - b) / (4 * u128::from(n)));
    assert_eq!(counts.iter().sum::<usize>() as u128, total);
    (counts[0], counts[1], counts[2])
}
fn gcd(mut x: u128, mut y: u128) -> u128 {
    while y != 0 {
        (x, y) = (y, x % y);
    }
    x
}
fn check_stabilizers(p: u128, n: u32) {
    let (big, a, b, e, _) = values(p, n);
    assert_eq!(gcd(a, b), 1);
    assert_eq!(gcd(big, p * p - 1), e * u128::from(n).max(1) + (1 - e));
    for j in 1..2 * n {
        let power = p.pow(j);
        for plus in [false, true] {
            let g = gcd(big, if plus { power + 1 } else { power - 1 });
            let want = if j == n {
                if plus {
                    b
                } else {
                    a
                }
            } else if e == 1
                && (if plus {
                    (power + 1) % u128::from(n)
                } else {
                    (power - 1) % u128::from(n)
                }) == 0
            {
                u128::from(n)
            } else {
                1
            };
            assert_eq!(g, want, "p={p},n={n},j={j},plus={plus}");
        }
    }
    if e == 1 {
        for k in 1..n {
            assert_eq!(
                orbit(u128::from(k) * big / u128::from(n), p, n, big).len(),
                2
            );
        }
    }
}
fn statement() -> IdentityStatement {
    let n = Expr::var("n");
    let e = Expr::var("e");
    let a = Expr::var("A");
    let b = Expr::var("B");
    let big = Expr::var("N");
    let nm1 = Expr::plus(n.clone(), Expr::Const(-1));
    IdentityStatement{name:"ISO-1 prime-degree cleared orbit-sum numerator".into(),source:"research/iso1_weak_classes_20261007/larger_fields_20261009/REPORT.tex; orbit sizes 2,2*n,4*n".into(),variables:vec!["n".into(),"e".into(),"A".into(),"B".into(),"N".into()],degree_bound:3,
      lhs:Expr::sum(vec![Expr::product(vec![Expr::Const(2),n,e.clone(),nm1.clone()]),Expr::times(Expr::Const(2),Expr::sum(vec![a.clone(),b.clone(),Expr::Const(-2),Expr::minus(Expr::times(e.clone(),nm1.clone()))])),Expr::sum(vec![big.clone(),Expr::Const(1),Expr::minus(a.clone()),Expr::minus(b.clone())])]),
      rhs:Expr::sum(vec![big,a,b,Expr::Const(-3),Expr::product(vec![Expr::Const(2),e,nm1.clone(),nm1])])}
}
fn replay(records: &[IdentityCertificate]) {
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].statement, statement());
    assert!(matches!(
        check_with_extra(&records[0], Some((2026100921, 64))).unwrap(),
        Verdict::Accept {
            points: 32,
            extra_points: 64,
            ..
        }
    ));
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.first().is_some_and(|s| s == "--issue") {
        let r=issue(statement(),20261009,32,"Exact cleared orbit-sum algebra; prime-degree stabilizers have a separate deductive proof and cyclic-group audits.",None).unwrap();
        replay(std::slice::from_ref(&r));
        eprintln!("{} accepted", r.id);
        println!("{}", serde_json::to_string_pretty(&vec![r]).unwrap());
        return;
    }
    if let Some(path) = args.first() {
        replay(
            &serde_json::from_slice::<Vec<IdentityCertificate>>(&std::fs::read(path).unwrap())
                .unwrap(),
        );
    }
    for (p, n) in [
        (3, 3),
        (5, 3),
        (7, 3),
        (13, 3),
        (3, 5),
        (5, 5),
        (7, 5),
        (3, 7),
    ] {
        let (big, _, _, _, total) = values(p, n);
        let c = exhaustive(p, n);
        println!("p={p},n={n},parameters={},orbits={total},size2={},size2n={},size4n={},status=EXACT_MATCH",big-1,c.0,c.1,c.2);
    }
    let mut checked = 0;
    for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] {
        for n in [3, 5, 7] {
            check_stabilizers(p, n);
            checked += 1;
        }
    }
    println!("stabilizer_cases={checked},status=EXACT_MATCH");
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn small_group_orbit_populations_agree() {
        for (p, n) in [(3, 3), (5, 3), (3, 5), (3, 7)] {
            exhaustive(p, n);
        }
    }
    #[test]
    fn exceptional_stabilizers_agree() {
        for p in [3, 5, 7, 11, 13, 19, 29] {
            for n in [3, 5, 7] {
                check_stabilizers(p, n);
            }
        }
    }
    #[test]
    fn changed_numerator_is_refused() {
        let mut s = statement();
        s.rhs = Expr::plus(s.rhs, Expr::Const(1));
        assert!(issue(s, 20261009, 32, "mutation", None).is_err());
    }
}
