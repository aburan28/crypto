mod common;
use common::*;
use eigen_seta::cm::{chain3_vartime, mul_vartime, point_of_order, random_point};
use eigen_seta::isog3::chain3;
use eigen_seta::mont::{ladder, ladder3pt, recover_y, xtpl, MontCurve, XPoint};
use eigen_seta::scalar::Scalar;
use eigen_seta::weier::{WCurve, WIso, WPoint};
use num_bigint::BigUint;
use rand_core::{OsRng, RngCore};

fn rand_scalar(rng: &mut OsRng, bits: usize) -> (Scalar, BigUint) {
    let limbs = bits.div_ceil(64);
    let mut l = vec![0u64; limbs];
    for x in l.iter_mut() {
        *x = rng.next_u64();
    }
    let s = Scalar { l }.mod_2k(bits);
    let b = s.to_biguint();
    (s, b)
}

#[test]
fn montgomery_against_weierstrass() {
    let mut rng = OsRng;
    let a = F::zero(); // y^2 = x^3 + x, (p+1)^2 points
    let mc = MontCurve::new(a);
    let pc = mc.proj();
    let w = WCurve::from_montgomery(a);
    for _ in 0..10 {
        let p = random_point(&w, &mut rng);
        let q = random_point(&w, &mut rng);
        let (k, kb) = rand_scalar(&mut rng, 200);
        let (r0, r1) = ladder(&XPoint::from_affine(p.x), &k, 200, &pc);
        let kp = mul_vartime(&w, &p, &kb);
        assert_eq!(r0.affine(), kp.x, "ladder [k]P");
        let k1p = w.add(&kp, &p);
        assert_eq!(r1.affine(), k1p.x, "ladder [k+1]P");
        let (x, y, z) = recover_y(&a, &p.x, &p.y, &r0, &r1);
        let zi = z.invert();
        assert_eq!(x * zi, kp.x);
        assert_eq!(y * zi, kp.y, "y recovery");
        let pq = w.sub(&p, &q);
        let m3 = ladder3pt(
            &XPoint::from_affine(p.x),
            &XPoint::from_affine(q.x),
            &XPoint::from_affine(pq.x),
            &k,
            200,
            &pc,
        );
        let expect = w.add(&p, &mul_vartime(&w, &q, &kb));
        assert_eq!(m3.affine(), expect.x, "ladder3pt");
        let t = xtpl(&XPoint::from_affine(p.x), &pc);
        assert_eq!(
            t.affine(),
            mul_vartime(&w, &p, &BigUint::from(3u32)).x,
            "xtpl"
        );
    }
    // 3-isogeny chain: x-only vs full-point Velu, compared through an isomorphism
    let n = 12;
    let k = point_of_order(&w, 3, n, &mut rng);
    let r = random_point(&w, &mut rng);
    let mut xs = [XPoint::from_affine(r.x)];
    let cod = chain3(&pc, &XPoint::from_affine(k.x), n, &mut xs);
    let a2 = cod.affine_a();
    let mut ws = [r];
    let wcod = chain3_vartime(&w, &k, n, &mut ws);
    let mcod = WCurve::from_montgomery(a2);
    assert_eq!(mcod.j_invariant(), wcod.j_invariant(), "3-chain codomain j");
    let iso = WIso::find(&wcod, &mcod);
    let im = iso.eval(&ws[0]);
    assert_eq!(im.x, xs[0].affine(), "3-chain image x");
    let _ = WPoint::<F>::infinity();
}
