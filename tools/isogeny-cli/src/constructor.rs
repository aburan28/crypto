//! Public-curve kernel construction on one selected Frobenius eigenline.
#[path = "../../../research/p192_p224_large_degree_isogenies_20261009/extension.rs"]
mod extension;
use extension::Extension;
use isogeny_algos::{
    api::{self, PrimeCurve, PrimeFieldInt},
    bigint::Big,
    curve::{padd, pmul_big, random_point_f, Curve, Pt},
    field::{Field, Rng},
    fpm::FpM,
    icv1,
    int::Int,
    json::Json,
    poly,
};
fn array(v: &[Int]) -> Json {
    Json::Arr(v.iter().map(Json::str).collect())
}
fn identity(c: &PrimeCurve, n: &Int) -> Json {
    let v = icv1::prime(&c.p, &c.a, &c.b, n).unwrap();
    Json::obj(vec![
        ("slug", Json::str(v.slug)),
        ("icv1", Json::str(v.icv1)),
        ("model_json", Json::str(v.model_json)),
    ])
}
fn stage(name: &str) {
    eprintln!("kernel-first stage={name}");
}
fn product<F: Field>(f: &F, xs: &[F::E]) -> Vec<F::E> {
    if xs.len() <= 16 {
        poly::from_roots(f, xs)
    } else {
        let middle = xs.len() / 2;
        poly::mul_raw(f, &product(f, &xs[..middle]), &product(f, &xs[middle..]))
    }
}
fn construct<F: PrimeFieldInt, const D: usize>(
    base: F,
    c: &PrimeCurve,
    n: &Int,
    ell: u64,
    lambda: u64,
    other: u64,
) -> Vec<Int> {
    stage("extension_begin");
    let f = Extension::<_, D>::search(base, 1);
    stage("extension_verified");
    let p = &c.p;
    let t = &(p + &Int::one()) - n;
    let mut v0 = Int::from(2i64);
    let mut v1 = t.clone();
    for _ in 2..=D {
        let next = &(&t * &v1) - &(p * &v0);
        v0 = v1;
        v1 = next;
    }
    let full = &(&Int::from_big(&f.q()) + &Int::one()) - &v1;
    assert!(!full.is_neg());
    let mut cofactor = full.mag().clone();
    let mut valuation = 0;
    while cofactor.rem_small(ell) == 0 {
        cofactor = cofactor.divrem_small(ell).0;
        valuation += 1;
    }
    assert!(valuation > 0);
    let mut reduction = Big::from_u64(1);
    for _ in 1..valuation {
        reduction = reduction.mul_small(ell);
    }
    let e = Curve::new(f.embed(f.base.elem(&c.a)), f.embed(f.base.elem(&c.b)));
    let mut rng = Rng::new(1);
    let mut chosen = None;
    stage("torsion_begin");
    for _ in 0..32 {
        let p0 = random_point_f(&f, &e, &mut rng);
        let point = pmul_big(&f, &e, &p0, &cofactor);
        let point = pmul_big(&f, &e, &point, &reduction);
        if point == Pt::Inf {
            continue;
        }
        assert_eq!(pmul_big(&f, &e, &point, &Big::from_u64(ell)), Pt::Inf);
        let Pt::Aff(x, y) = point else { unreachable!() };
        let frob = Pt::Aff(f.pow_big(x, c.p.mag()), f.pow_big(y, c.p.mag()));
        let term = pmul_big(&f, &e, &point, &Big::from_u64(other));
        let neg = match term {
            Pt::Inf => Pt::Inf,
            Pt::Aff(x, y) => Pt::Aff(x, f.neg(y)),
        };
        let point = padd(&f, &e, &frob, &neg);
        if point == Pt::Inf {
            continue;
        }
        let Pt::Aff(x, y) = point else { unreachable!() };
        let frob = Pt::Aff(f.pow_big(x, c.p.mag()), f.pow_big(y, c.p.mag()));
        assert_eq!(frob, pmul_big(&f, &e, &point, &Big::from_u64(lambda)));
        assert_eq!(pmul_big(&f, &e, &point, &Big::from_u64(ell)), Pt::Inf);
        chosen = Some(x);
        break;
    }
    let x = chosen.expect("no selected torsion point in the declared point budget");
    stage("torsion_verified");
    let degree = ((ell - 1) / 2) as usize;
    let xs = isogeny_algos::kernel::xonly::multiples_x_fast(&f, &e, x, degree + 1)
        .expect("kernel orbit construction failed");
    assert_eq!(xs[degree], xs[degree - 1]);
    stage("kernel_product_begin");
    let h = product(&f, &xs[..degree]);
    let h: Vec<_> = h
        .into_iter()
        .map(|a| {
            assert!(
                a[1..].iter().all(|v| *v == f.base.zero()),
                "kernel does not descend to the source field"
            );
            f.base.int(a[0])
        })
        .collect();
    stage("kernel_in_base_field");
    h
}
pub fn selection(name: &str, ell: u64) -> Option<(usize, u64, u64)> {
    if !matches!(name, "p192" | "p224")
        || !(3..=1_000_000).contains(&ell)
        || !isogeny_algos::field::is_prime(ell)
    {
        return None;
    }
    let roots = crate::eigenvalues(name, ell)?;
    let mut choices = [
        (
            crate::multiplicative_order(roots[0], ell),
            roots[0],
            roots[1],
        ),
        (
            crate::multiplicative_order(roots[1], ell),
            roots[1],
            roots[0],
        ),
    ];
    choices.sort();
    let (order, lambda, other) = choices[0];
    (order <= 16).then_some((order.max(2) as usize, lambda, other))
}
fn dispatch<F: PrimeFieldInt>(
    base: F,
    c: &PrimeCurve,
    n: &Int,
    ell: u64,
    d: usize,
    lambda: u64,
    other: u64,
) -> Vec<Int> {
    macro_rules! choose { ($($degree:literal),+) => { match d {
        $($degree => construct::<_, $degree>(base, c, n, ell, lambda, other),)+
        _ => panic!("unsupported extension degree"),
    } }; }
    choose!(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
}
pub fn run(name: &str, ell: u64) -> Json {
    run_model(name, ell, &PrimeCurve::preset(name).unwrap())
}
pub fn run_model(name: &str, ell: u64, c: &PrimeCurve) -> Json {
    let n = crate::standard_order(name);
    assert_eq!(c.p, PrimeCurve::preset(name).unwrap().p);
    let (d, lambda, other) = selection(name, ell).expect("unsupported eigenvalue order");
    let h = match name {
        "p192" => dispatch(FpM::<3>::new(c.p.mag()), c, &n, ell, d, lambda, other),
        "p224" => dispatch(FpM::<4>::new(c.p.mag()), c, &n, ell, d, lambda, other),
        _ => unreachable!(),
    };
    stage("map_begin");
    let r = api::isogeny_from_kernel(c, &h, ell, 1).unwrap();
    assert!(r.verified());
    stage("map_checks_passed");
    let cod = PrimeCurve::new(c.p.clone(), r.codomain.0.clone(), r.codomain.1.clone()).unwrap();
    let checks = Json::Arr(
        r.checks
            .iter()
            .map(|v| {
                Json::obj(vec![
                    ("name", Json::str(&v.name)),
                    ("status", Json::str(v.status.as_str())),
                    ("detail", Json::str(&v.detail)),
                ])
            })
            .collect(),
    );
    let iso = Json::obj(vec![
        ("verified", Json::Bool(true)),
        ("kernel", array(&r.kernel)),
        (
            "x_map",
            Json::obj(vec![("num", array(&r.map_num)), ("den", array(&r.map_den))]),
        ),
        (
            "codomain",
            Json::obj(vec![
                ("a", Json::str(&cod.a)),
                ("b", Json::str(&cod.b)),
                ("j", Json::str(&r.j_codomain)),
                ("icv1", identity(&cod, &n)),
            ]),
        ),
        ("checks", checks),
    ]);
    let doc = Json::obj(vec![
        ("status", Json::str("PASS")),
        ("operation", Json::str("kernel-first")),
        ("scope", Json::str("one_frobenius_eigenline")),
        ("degree_coverage", Json::str("PARTIAL")),
        ("expected_total_eigenlines", Json::Num(2)),
        ("constructed_eigenlines", Json::Num(1)),
        (
            "curve",
            Json::obj(vec![
                ("preset", Json::str(name)),
                ("p", Json::str(&c.p)),
                ("a", Json::str(&c.a)),
                ("b", Json::str(&c.b)),
            ]),
        ),
        (
            "result",
            Json::obj(vec![
                ("ell", Json::Num(ell as i64)),
                ("source_icv1", identity(c, &n)),
                ("isogenies", Json::Arr(vec![iso])),
            ]),
        ),
    ]);
    doc
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn selections_keep_old_routes_and_enable_new_degree_nine_routes() {
        assert_eq!(selection("p192", 10453), Some((3, 270, 8274)));
        assert_eq!(selection("p224", 1471), Some((5, 554, 1178)));
        assert_eq!(selection("p192", 2341), Some((9, 850, 1877)));
        assert_eq!(selection("p224", 1117), Some((9, 591, 698)));
        assert_eq!(selection("p192", 1021), None);
        assert_eq!(selection("p224", 1119), None);
        assert_eq!(selection("p224", u64::MAX), None);
    }
    fn field_laws<const D: usize>() {
        let f = Extension::<_, D>::search(isogeny_algos::field::Zp::new(101), 991);
        let mut rng = Rng::new(2721);
        for _ in 0..8 {
            let (a, b, c) = (f.random(&mut rng), f.random(&mut rng), f.random(&mut rng));
            assert_eq!(f.mul(a, f.add(b, c)), f.add(f.mul(a, b), f.mul(a, c)));
            assert_eq!(f.pow_big(a, &f.q()), a);
            if !f.is_zero(a) {
                assert_eq!(f.mul(a, f.inv(a)), f.one());
            }
        }
    }
    #[test]
    fn every_enabled_extension_passes_field_laws() {
        field_laws::<2>();
        field_laws::<3>();
        field_laws::<4>();
        field_laws::<5>();
        field_laws::<6>();
        field_laws::<7>();
        field_laws::<8>();
        field_laws::<9>();
        field_laws::<10>();
        field_laws::<11>();
        field_laws::<12>();
        field_laws::<13>();
        field_laws::<14>();
        field_laws::<15>();
        field_laws::<16>();
    }
}
