use isogeny_algos::api::{self, CountMethod, PrimeCurve, Status};
use isogeny_algos::bigint::Big;
use isogeny_algos::curve::Curve;
use isogeny_algos::field::{is_prime, Rng, Zp};
use isogeny_algos::icv1;
use isogeny_algos::int::Int;
use std::process::Command;

fn dec(s: &str) -> Int {
    Int::from_big(&Big::from_dec(s))
}

#[test]
fn p192_order_and_split_degree_maps_are_supported() {
    let c = PrimeCurve::preset("p192").unwrap();
    let n = dec("6277101735386680763835789423176059013767194773182842284081");
    let count = api::check_order(&c, &n, 1);
    assert!(count.certificate.certified && count.twist_certificate.certified);
    let maps = api::isogenies(&c, 13, 1).unwrap();
    assert_eq!(maps.len(), 2);
    assert!(maps.iter().all(|r| r.verified()));
}

#[test]
fn p224_small_low_limb_does_not_reject_modular_construction() {
    use isogeny_algos::{field::Field, find::modpoly::Phi, fpm::FpM};
    let c = PrimeCurve::preset("p224").unwrap();
    let f = FpM::<4>::new(c.p.mag());
    assert_eq!(f.char(), 1);
    assert!(f.characteristic_exceeds(u64::MAX));
    // A second construction route certifies the small modular polynomial exactly.
    assert_eq!(
        Phi::compute_hecke(&f, 3).c,
        Phi::compute_linear_algebra(&f, 3).c
    );
    let n = dec("26959946667150639794667015087019625940457807714424391721682722368061");
    let t = &(&c.p + &Int::one()) - &n;
    let ell = [3u64, 5, 7, 11, 13, 17, 19, 23, 29, 31]
        .into_iter()
        .find(|&l| {
            let tm = t.mod_u64(l);
            let pm = c.p.mod_u64(l);
            (0..l)
                .filter(|&x| (x * x + l - tm * x % l + pm) % l == 0)
                .count()
                == 2
        })
        .expect("a split small prime for the public standard model");
    let maps = api::isogenies(&c, ell, 1).unwrap();
    assert_eq!(maps.len(), 2);
    assert!(maps.iter().all(|m| m.verified()));
    let small = FpM::<1>::from_u64_modulus(101);
    assert!(small.characteristic_exceeds(100));
    assert!(!small.characteristic_exceeds(101));
}

#[test]
fn field_dispatch_handles_every_limb_boundary() {
    use isogeny_algos::{field::Field, fp2::Fp2, fpm::FpM};

    // A quadratic wrapper must retain its multiword base characteristic.
    let mersenne = FpM::<2>::new(&Big::from_dec("170141183460469231731687303715884105727"));
    let extension = Fp2::new(mersenne);
    assert_eq!(extension.char(), u64::MAX);
    assert!(extension.characteristic_exceeds(u64::MAX));

    struct Multiply(Int, Int);
    impl api::FieldTask for Multiply {
        type Out = Int;
        fn run<F: api::PrimeFieldInt>(self, f: &F) -> Int {
            f.int(f.mul(f.elem(&self.0), f.elem(&self.1)))
        }
    }
    // Montgomery arithmetic also works over these odd composite moduli; this checks
    // dispatch and reduction at both ends of every limb width, not primality.
    for bits in [
        63, 64, 65, 128, 129, 192, 193, 256, 257, 320, 321, 384, 385, 448, 449, 512,
    ] {
        let p = &Int::from_big(&Big::from_u64(1).shl(bits)) - &Int::one();
        let a = &p - &Int::from(2i64);
        let b = &p - &Int::from(3i64);
        assert_eq!(
            api::with_field(&p, Multiply(a, b)),
            Int::from(6i64),
            "bits={bits}"
        );
    }
}

/// Standard curves: the orders come out certified and the ICV1 slugs are the ones the crypto
/// repository's curve registry (`docs/curves/registry.json`) lists for P-256 and secp256k1.
#[test]
fn standard_curves_certified_and_named_as_in_the_registry() {
    let p256 = PrimeCurve::preset("p256").unwrap();
    let n256 =
        dec("115792089210356248762697446949407573529996955224135760342422259061068512044369");
    let pc = api::check_order(&p256, &n256, 1);
    assert!(pc.certificate.certified && pc.twist_certificate.certified);
    let id = icv1::prime(&p256.p, &p256.a, &p256.b, &pc.order).unwrap();
    assert_eq!(
        id.slug,
        "icv1-fp256-t89188191154553853111372247798585809583-f188c491"
    );
    // secp256k1 has j = 0: counted by complex multiplication
    let k1 = PrimeCurve::preset("secp256k1").unwrap();
    let pc = api::count_points(&k1, 1).unwrap();
    assert_eq!(pc.method, CountMethod::Cm);
    assert_eq!(
        pc.order,
        dec("115792089237316195423570985008687907852837564279074904382605163141518161494337")
    );
    assert!(pc.certificate.certified);
    let id = icv1::prime(&k1.p, &k1.a, &k1.b, &pc.order).unwrap();
    assert_eq!(
        id.slug,
        "icv1-fp256-t432420386565659656852420866390673177327-76dadd18"
    );
    // a wrong claimed order is rejected (here: off by one)
    let bad = api::check_order(&p256, &(&n256 + &Int::one()), 1);
    assert!(!bad.certificate.certified);
    assert!(bad
        .certificate
        .checks
        .iter()
        .any(|c| c.status == Status::Fail));
}

/// SEA (48-bit p), CM (j = 0 and 1728) and BSGS give the order the crate's exact BSGS gives.
#[test]
fn counts_match_bsgs() {
    let mut rng = Rng::new(7);
    let mut p = (1u64 << 47) + 1001;
    while !is_prime(p) {
        p += 2;
    }
    let fp = Zp::new(p);
    for (a, b) in [(3i64, 1234567i64), (-3, 98765), (0, 5), (7, 0)] {
        let c = PrimeCurve::new(Int::from(p), Int::from(a), Int::from(b)).unwrap();
        let pc = api::count_points(&c, 3).unwrap();
        let e = Curve::new(
            Int::from(a).modulo(&Int::from(p)).to_i128().unwrap() as u64,
            Int::from(b).modulo(&Int::from(p)).to_i128().unwrap() as u64,
        );
        let n = isogeny_algos::curve::order(&fp, &e, &mut rng);
        assert_eq!(pc.order, Int::from(n), "a = {a}, b = {b}");
        assert_eq!(
            &pc.order + &pc.twist_order,
            &Int::from(2 * p) + &Int::from(2i64)
        );
        assert!(pc.certificate.certified || pc.twist_certificate.certified);
        // a transferred certificate rests on the other one's own point orders
        for (c, other) in [
            (&pc.certificate, &pc.twist_certificate),
            (&pc.twist_certificate, &pc.certificate),
        ] {
            if c.via_twist {
                assert!(c.certified && other.certified && !other.via_twist);
            }
        }
    }
}

/// The rational l-isogenies found are exactly the roots of Phi_l(j, Y) in F_p, each verified;
/// Velu from each kernel polynomial gives the same codomain.
#[test]
fn isogenies_match_modular_roots_and_velu() {
    let p = Int::from(2305843009213693951i64); // 2^61 - 1
    let c = PrimeCurve::new(p.clone(), Int::from(3i64), Int::from(1234567i64)).unwrap();
    let mut seen = 0;
    for ell in [3u64, 5, 7, 11, 13] {
        let isos = api::isogenies(&c, ell, 5).unwrap();
        let (_, roots) = api::modular_polynomial_at(&p, ell as usize, &c.j(), 5).unwrap();
        let mut js: Vec<Int> = isos.iter().map(|r| r.j_codomain.clone()).collect();
        js.sort();
        assert_eq!(js, roots, "l = {ell}");
        for r in &isos {
            assert!(r.verified(), "l = {ell}: {:?}", r.checks);
            let again = api::isogeny_from_kernel(&c, &r.kernel, ell, 9).unwrap();
            assert_eq!(again.codomain, r.codomain);
            seen += 1;
        }
    }
    assert!(seen > 0);
    // l = 2 from the rational 2-torsion
    for r in api::isogenies(&c, 2, 5).unwrap() {
        assert!(r.verified());
        assert_eq!(r.kernel.len(), 2);
    }
}

/// The command line: one JSON object, exit 0 when every check passes, 1 on a failed check,
/// 2 on a usage error.
#[test]
fn cli_records_and_exit_codes() {
    let bin = env!("CARGO_BIN_EXE_isogeny-algos");
    let run = |args: &[&str]| {
        let out = Command::new(bin).args(args).output().unwrap();
        (
            out.status.code().unwrap(),
            String::from_utf8(out.stdout).unwrap(),
        )
    };
    let (code, out) = run(&[
        "count", "--p", "10935329", "--a", "5320418", "--b", "8535318",
    ]);
    assert_eq!(code, 0, "{out}");
    assert!(out.starts_with("{\"schema\":\"isogeny-algos/v1\""));
    assert!(
        out.contains("\"icv1\":\"ICV1:fp-10935329:1577:10933753:2786525:unk:unk:r:773361558ca8\"")
    );
    assert!(out.contains("\"status\":\"PASS\""));
    assert_eq!(out.trim().lines().count(), 1);
    let (code, out) = run(&["count", "--curve", "p256", "--order", "1"]);
    assert_eq!(code, 1, "{out}");
    assert!(out.contains("\"status\":\"FAIL\""));
    let (code, _) = run(&["count", "--p", "15"]);
    assert_eq!(code, 2);
    let (code, out) = run(&["modpoly", "--ell", "3", "--integer"]);
    assert_eq!(code, 0);
    assert!(out.contains("\"1855425871872000000000\""));
    let (code, out) = run(&["--version"]);
    assert_eq!(code, 0);
    assert_eq!(
        out.trim(),
        format!(
            "{{\"schema\":\"isogeny-algos/v1\",\"tool\":{{\"name\":\"isogeny-algos\",\"version\":\"{}\"}}}}",
            env!("CARGO_PKG_VERSION")
        )
    );
}
