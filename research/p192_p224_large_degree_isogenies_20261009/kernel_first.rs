//! Public-curve kernel construction on one selected Frobenius eigenline.
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
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 9);
    assert_eq!(args[0], "kernel-first");
    let get = |key: &str| args.windows(2).find(|a| a[0] == key).unwrap()[1].as_str();
    let name = get("--curve");
    let ell: u64 = get("--ell").parse().unwrap();
    assert_eq!(get("--seed"), "1");
    let c = PrimeCurve::preset(name).unwrap();
    let n = Int::from_big(&Big::from_dec(get("--order")));
    let h = match (name, ell) {
        ("p192", 10453) => construct::<_, 3>(FpM::<3>::new(c.p.mag()), &c, &n, ell, 270, 8274),
        ("p224", 1471) => construct::<_, 5>(FpM::<4>::new(c.p.mag()), &c, &n, ell, 554, 1178),
        _ => panic!("case outside the frozen kernel protocol"),
    };
    stage("map_begin");
    let r = api::isogeny_from_kernel(&c, &h, ell, 1).unwrap();
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
                ("source_icv1", identity(&c, &n)),
                ("isogenies", Json::Arr(vec![iso])),
            ]),
        ),
    ]);
    println!("{}", doc.dump());
}
