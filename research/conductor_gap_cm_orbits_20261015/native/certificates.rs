use super::*;
use crypto_lib::cryptanalysis::identity_certificate::{self as ic, Expr as E, IdentityStatement};
fn mul(a: E, b: E) -> E {
    E::times(a, b)
}
fn sq(a: E) -> E {
    mul(a.clone(), a)
}
fn cube(a: E) -> E {
    mul(sq(a.clone()), a)
}
fn add(a: E, b: E) -> E {
    E::plus(a, b)
}
fn c(k: i64, a: E) -> E {
    mul(E::Const(k), a)
}
pub fn run() {
    let root = Path::new(ROOT).join("certificates");
    fs::create_dir_all(&root).unwrap();
    let a = E::var("a");
    let b = E::var("b");
    let ell = E::var("ell");
    let ab = mul(a.clone(), mul(b.clone(), ell.clone()));
    let bell = sq(mul(b.clone(), ell));
    let norm = add(add(c(4, sq(a.clone())), c(4, ab)), c(20, bell.clone()));
    let norm2 = add(sq(add(c(2, a), mul(b.clone(), E::var("ell")))), c(19, bell));
    let u = E::var("u");
    let v = E::var("v");
    let plus = add(u.clone(), v.clone());
    let minus = add(u.clone(), E::minus(v.clone()));
    let cases = vec![
        ("cm-norm-numerator", vec!["a", "b", "ell"], 4, norm, norm2),
        (
            "paired-second-moment",
            vec!["u", "v"],
            2,
            add(sq(plus.clone()), sq(minus.clone())),
            c(2, add(sq(u.clone()), sq(v.clone()))),
        ),
        (
            "paired-third-moment",
            vec!["u", "v"],
            3,
            add(cube(plus), cube(minus)),
            add(c(2, cube(u.clone())), c(6, mul(u, sq(v)))),
        ),
    ];
    let mut receipts = Vec::new();
    for (name, variables, degree, lhs, rhs) in cases {
        let cert = ic::issue(
            IdentityStatement {
                name: name.into(),
                source: format!(
                    "{ROOT}/native/sqrt_velu.cpp; lean/README.md specifies formalization limits"
                ),
                variables: variables.into_iter().map(str::to_string).collect(),
                degree_bound: degree,
                lhs,
                rhs,
            },
            20261015,
            64,
            "native issuer; domain-separated checker uses an additional seed and 64 unseen points",
            None,
        )
        .unwrap();
        let verdict = ic::check_with_extra(&cert, Some((20261016, 64))).unwrap();
        write(
            root.join(format!("{name}.json")),
            &serde_json::to_value(&cert).unwrap(),
        );
        receipts.push(json!({"id":cert.id,"verdict":format!("{:?}",verdict)}));
    }
    write(root.join("independent-checks.json"), &json!(receipts));
}
