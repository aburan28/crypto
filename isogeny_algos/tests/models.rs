use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::kernel::models::*;
use isogeny_algos::kernel::montgomery;
use isogeny_algos::path::csidh::Csidh;

/// A point of exact order l (l | n = #E), robust to a non-cyclic l-part.
fn point_of_order(fq: &Zp, ew: &Curve<u64>, n: u64, ell: u64, rng: &mut Rng) -> Pt<u64> {
    let mut lv = 1u64;
    while n.is_multiple_of(lv * ell) {
        lv *= ell;
    }
    loop {
        let r = ew.random_point(fq, rng);
        let mut t = pmul(fq, ew, &r, (n / lv) as u128);
        if t == Pt::Inf {
            continue;
        }
        loop {
            let nt = pmul(fq, ew, &t, ell as u128);
            if nt == Pt::Inf {
                return t;
            }
            t = nt;
        }
    }
}

/// Twisted Edwards on CSIDH toy curves (p + 1 = 4 prod l): codomain vs Montgomery Velu, images on
/// the codomain, homomorphism, kernel to the identity; explicit/one-variable forms for a = 1.
#[test]
fn edwards_isogenies() {
    let cs = Csidh::with_n_primes(8).unwrap();
    let fp = cs.fp;
    let p = cs.p();
    let mut rng = Rng::new(1400);
    for (idx, &ell) in cs.primes.iter().enumerate().take(6) {
        let a = cs.action(0, &[1, -1, 0, 1], &mut rng);
        let ew = montgomery::to_weierstrass(&fp, a);
        let ed = edwards_from_montgomery(&fp, a);
        let third = fp.div(a, 3);
        let to_ed = |q: &Pt<u64>| match *q {
            Pt::Aff(xw, yw) => Some(mont_to_edwards_point(&fp, fp.sub(xw, third), yw)),
            Pt::Inf => None,
        };
        let k = loop {
            let r = random_point_f(&fp, &ew, &mut rng);
            let kk = pmul(&fp, &ew, &r, ((p + 1) / ell) as u128);
            if kk != Pt::Inf {
                break kk;
            }
        };
        let ke = to_ed(&k).unwrap();
        assert!(ed.on_curve(&fp, ke));
        let iso = edwards_isogeny(&fp, &ed, ke, ell);
        // codomain j vs Montgomery Velu
        let Pt::Aff(xk, _) = k else { unreachable!() };
        let xm = fp.sub(xk, third);
        let ms = montgomery::multiples(
            &fp,
            montgomery::a24(&fp, a),
            (xm, 1),
            ((ell - 1) / 2) as usize,
        );
        let a2 = montgomery::velu_codomain(&fp, a, &ms, ell);
        assert_eq!(
            iso.cod.j(&fp),
            montgomery::j_invariant(&fp, a2),
            "l = {ell}"
        );
        // points
        for _ in 0..5 {
            let (q1, q2) = (
                random_point_f(&fp, &ew, &mut rng),
                random_point_f(&fp, &ew, &mut rng),
            );
            let (Some(e1), Some(e2)) = (to_ed(&q1), to_ed(&q2)) else {
                continue;
            };
            let i1 = iso.eval_def(&fp, e1);
            let i2 = iso.eval_def(&fp, e2);
            assert!(iso.cod.on_curve(&fp, i1), "image on codomain");
            assert_eq!(
                iso.eval_def(&fp, ed.add(&fp, e1, e2)),
                iso.cod.add(&fp, i1, i2),
                "homomorphism"
            );
        }
        assert_eq!(iso.eval_def(&fp, ke), (0, 1), "kernel -> identity");
        let _ = idx;
    }
    // a = 1 (Edwards x^2 + y^2 = 1 + d x^2 y^2): ordinary Montgomery curves over a 30-bit field
    // with A + 2 a square, scaled by sqrt(A + 2)
    let fq = Zp::new(next_prime(1u64 << 30));
    let mut done = 0;
    while done < 4 {
        let a = fq.random(&mut rng);
        let ed = edwards_from_montgomery(&fq, a);
        if ed.d == 0 || ed.a == 0 {
            continue;
        }
        let Some(sa) = fq.sqrt(ed.a) else { continue };
        let ew = montgomery::to_weierstrass(&fq, a);
        if !is_smooth(&fq, &ew) {
            continue;
        }
        let n = order(&fq, &ew, &mut rng);
        let Some(&ell) = [3u64, 5, 7, 11, 13].iter().find(|&&l| n.is_multiple_of(l)) else {
            continue;
        };
        let third = fq.div(a, 3);
        let k = point_of_order(&fq, &ew, n, ell, &mut rng);
        let e1 = Edwards {
            a: 1u64,
            d: fq.div(ed.d, ed.a),
        };
        let Pt::Aff(xw, yw) = k else { unreachable!() };
        let (x0, y0) = mont_to_edwards_point(&fq, fq.sub(xw, third), yw);
        let k1 = (fq.mul(x0, sa), y0);
        assert!(e1.on_curve(&fq, k1));
        let iso = edwards_isogeny(&fq, &e1, k1, ell);
        let mut checked = 0;
        for _ in 0..10 {
            let q = ew.random_point(&fq, &mut rng);
            let Pt::Aff(qx, qy) = q else { continue };
            if qy == 0 {
                continue;
            }
            let (x, y) = mont_to_edwards_point(&fq, fq.sub(qx, third), qy);
            let pt = (fq.mul(x, sa), y);
            let d = iso.eval_def(&fq, pt);
            assert!(iso.cod.on_curve(&fq, d));
            assert_eq!(
                iso.eval_explicit(&fq, pt),
                d,
                "Theorem 2 = definition, l = {ell}"
            );
            assert_eq!(iso.eval_x_only(&fq, pt.0), d.0, "one-variable X, l = {ell}");
            checked += 1;
        }
        assert!(checked > 0);
        done += 1;
    }
}

/// Huff curves: codomain vs Weierstrass Velu (j), images on the codomain, homomorphism.
#[test]
fn huff_isogenies() {
    let mut rng = Rng::new(1401);
    let fp = Zp::new(next_prime(1u64 << 30));
    let mut done = 0;
    while done < 4 {
        let h = Huff {
            a: fp.random(&mut rng),
            b: fp.random(&mut rng),
        };
        if h.a == h.b || h.a == 0 || h.b == 0 {
            continue;
        }
        let ew = h.short_weierstrass(&fp);
        if !is_smooth(&fp, &ew) {
            continue;
        }
        let n = order(&fp, &ew, &mut rng);
        let Some(&ell) = [3u64, 5, 7, 11].iter().find(|&&l| n.is_multiple_of(l)) else {
            continue;
        };
        // a point of order l on y^2 = x^3 + (a+b) x^2 + ab x (shift between models: x_s = x + (a+b)/3)
        let shift = fp.div(fp.add(h.a, h.b), 3);
        let k = point_of_order(&fp, &ew, n, ell, &mut rng);
        let Pt::Aff(xs, ys) = k else { unreachable!() };
        let kh = h.from_weierstrass_point(&fp, (fp.sub(xs, shift), ys));
        assert!(h.on_curve(&fp, kh));
        let Some(iso) = huff_isogeny(&fp, &h, kh, ell) else {
            continue;
        };
        let reference = isogeny_algos::kernel::velu::velu_cyclic(&fp, &ew, &k, ell);
        assert_eq!(
            jinv(&fp, &iso.cod.short_weierstrass(&fp)),
            jinv(&fp, &reference.cod),
            "l = {ell}"
        );
        for _ in 0..5 {
            let (q1, q2) = (
                ew.random_point(&fp, &mut rng),
                ew.random_point(&fp, &mut rng),
            );
            let (Pt::Aff(x1, y1), Pt::Aff(x2, y2)) = (q1, q2) else {
                continue;
            };
            let p1 = h.from_weierstrass_point(&fp, (fp.sub(x1, shift), y1));
            let p2 = h.from_weierstrass_point(&fp, (fp.sub(x2, shift), y2));
            let (i1, i2) = (iso.eval(&fp, p1), iso.eval(&fp, p2));
            assert!(iso.cod.on_curve(&fp, i1));
            if let (Some(s), Some(si)) = (h.add(&fp, p1, p2), iso.cod.add(&fp, i1, i2)) {
                assert_eq!(iso.eval(&fp, s), si, "homomorphism");
            }
        }
        done += 1;
    }
}
