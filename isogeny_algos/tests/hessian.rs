use isogeny_algos::curve::{from_j, Curve};
use isogeny_algos::field::{Field, Rng, Zp};
use isogeny_algos::kernel::hessian::*;

/// All points of a twisted Hessian curve over F_p (p small), by brute force.
fn points(f: &Zp, e: &Hessian<u64>) -> Vec<HPt<u64>> {
    let p = f.p;
    let mut v = vec![];
    // Z = 1
    for x in 0..p {
        for y in 0..p {
            let pt = [x, y, 1];
            if e.on_curve(f, &pt) {
                v.push(pt);
            }
        }
    }
    // Z = 0: a X^3 + Y^3 = 0 with X = 1 or (X, Y) = (0, 1)
    for y in 0..p {
        if e.on_curve(f, &[1, y, 0]) {
            v.push([1, y, 0]);
        }
    }
    if e.on_curve(f, &[0, 1, 0]) {
        v.push([0, 1, 0]);
    }
    v
}

fn weierstrass_count(f: &Zp, c: &Curve<u64>) -> u64 {
    let p = f.p;
    1 + (0..p).map(|x| {
        let r = isogeny_algos::curve::rhs(f, c, x);
        if r == 0 { 1 } else if f.legendre(r) == 1 { 2 } else { 0 }
    }).sum::<u64>()
}

#[test]
fn hessian_group_law_isogenies_and_codomain() {
    let p = 1009u64;
    let f = Zp::new(p);
    let mut rng = Rng::new(5);
    let mut tested = std::collections::HashSet::new();
    let mut attempts = 0;
    while tested.iter().map(|t: &(u64, u64)| t.0).collect::<std::collections::HashSet<_>>().len() < 5 && attempts < 1000 {
        attempts += 1;
        let e = Hessian { a: f.random(&mut rng), d: f.random(&mut rng) };
        // non-singular: a (d^3 - 27 a) != 0
        if e.a == 0 || f.sub(f.mul(e.d, f.sq(e.d)), f.mul(27, e.a)) == 0 || e.d == 0 {
            continue;
        }
        let pts = points(&f, &e);
        let n = pts.len() as u64;
        // the j formula: a Weierstrass curve with this j has n or 2p + 2 - n points
        let j = e.j(&f);
        if j != 0 && j != 1728 {
            let w = weierstrass_count(&f, &from_j(&f, j));
            assert!(w == n || w == 2 * p + 2 - n, "j formula");
        }
        // group law: identity, inverse, associativity
        let o = e.identity(&f);
        for _ in 0..20 {
            let (a, b, c) = (pts[rng.below(n) as usize], pts[rng.below(n) as usize], pts[rng.below(n) as usize]);
            assert!(e.eq(&f, &e.add(&f, &a, &o), &a));
            assert!(e.eq(&f, &e.add(&f, &a, &e.neg(&a)), &o));
            let l = e.add(&f, &e.add(&f, &a, &b), &c);
            let r = e.add(&f, &a, &e.add(&f, &b, &c));
            assert!(e.on_curve(&f, &l) && e.eq(&f, &l, &r));
            assert!(e.eq(&f, &e.mul(&f, &a, n), &o));
        }
        for ell in [5u64, 7, 11, 13, 17] {
            if n % ell != 0 || tested.iter().any(|t| t.0 == ell) {
                continue;
            }
            // a point of order ell
            let k = loop {
                let k = e.mul(&f, &pts[rng.below(n) as usize], n / ell);
                if !e.eq(&f, &k, &o) {
                    break k;
                }
            };
            let iso = hessian_isogeny(&f, &e, &k, ell).unwrap();
            let cod = iso.cod;
            let mut checked = 0;
            assert!(e.eq(&f, &iso.eval(&f, &k), &cod.identity(&f)), "kernel -> identity");
            for _ in 0..20 {
                let (a, b) = (pts[rng.below(n) as usize], pts[rng.below(n) as usize]);
                let (ia, ib) = (iso.eval(&f, &a), iso.eval(&f, &b));
                if [ia, ib].iter().any(|v| v.iter().all(|&c| c == 0)) {
                    continue; // exceptional input of the rotated law
                }
                assert!(cod.on_curve(&f, &ia), "image on the codomain");
                let s = iso.eval(&f, &e.add(&f, &a, &b));
                if s.iter().all(|&c| c == 0) {
                    continue;
                }
                assert!(cod.eq(&f, &s, &cod.add(&f, &ia, &ib)), "homomorphism");
                checked += 1;
            }
            assert!(checked >= 15, "only {checked} non-exceptional checks");
            // isogenous curves have the same number of points
            assert_eq!(points(&f, &cod).len() as u64, n);
            tested.insert((ell, n));
        }
    }
    eprintln!("tested (l, #E): {tested:?}");
    assert!(tested.len() >= 5, "{tested:?}");
}
