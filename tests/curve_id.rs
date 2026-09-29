//! ICV1 identities against the reference implementation.
//!
//! Every vector here was produced by `python3 scripts/curve_id.py`; the
//! Rust module must reproduce the reference byte for byte.

use crypto_lib::cryptanalysis::curve_id::{
    binary, binary_degree, koblitz, normalise, prime, resolve, same_curve,
};
use crypto_lib::cryptanalysis::ic_boundary::{
    koblitz_instance, random_binary_instance, roster_prime_instance,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use num_traits::One;

#[test]
fn koblitz_matches_the_reference() {
    let f41 = (BigUint::one() << 41u32) | BigUint::from(0b1001u32);
    let id = koblitz(0, 41, &f41, &BigUint::from(2_199_025_563_772u64)).unwrap();
    assert_eq!(
        id.icv1,
        "ICV1:f2m-41-39c74c32:-2308219:2199025563772:0x1:-7:unk:r:7f48b14a43b8"
    );
    assert_eq!(id.slug, "icv1-f2m41-tm2308219-7f48b14a");
    assert_eq!(
        id.model_json,
        r#"{"a":"0x0","b":"0x1","field":"f2m-41-39c74c32","form":"y^2+xy=x^3+a*x^2+b","modulus":"0x20000000009","v":"1"}"#
    );
    let id = koblitz(1, 11, &BigUint::from(0x805u32), &BigUint::from(1982u32)).unwrap();
    assert_eq!(id.slug, "icv1-f2m11-t67-05f5aa36");
}

#[test]
fn wide_koblitz_matches_the_reference() {
    let f = (BigUint::one() << 131u32) | BigUint::from((1u32 << 13) | 0b111);
    let order: BigUint = "2722258935367507707729280517973639940516".parse().unwrap();
    let id = koblitz(0, 131, &f, &order).unwrap();
    assert_eq!(
        id.icv1,
        "ICV1:f2m-131-5cc3e086:-22283658519494248867:2722258935367507707729280517973639940516:0x1:-7:unk:r:115e0dc5d9eb"
    );
}

#[test]
fn random_binary_matches_the_reference() {
    let id = binary(
        27,
        &BigUint::from(0x800_0027u32),
        &BigUint::one(),
        &BigUint::from(0x84_5462u32),
        &BigUint::from(134_205_186u32),
        None,
    )
    .unwrap();
    assert_eq!(
        id.icv1,
        "ICV1:f2m-27-070a37b6:12543:134205186:0x5fa276d:unk:unk:r:569dca8b7ccb"
    );
    assert_eq!(id.slug, "icv1-f2m27-t12543-569dca8b");
}

#[test]
fn prime_matches_the_reference() {
    let big = |v: u64| BigUint::from(v);
    let id = prime(
        &big(10_935_329),
        &big(5_320_418),
        &big(8_535_318),
        &big(10_933_753),
    )
    .unwrap();
    assert_eq!(
        id.icv1,
        "ICV1:fp-10935329:1577:10933753:2786525:unk:unk:r:773361558ca8"
    );
    assert_eq!(id.slug, "icv1-fp24-t1577-77336155");
    let id = prime(&big(827), &big(1), &big(15), &big(823)).unwrap();
    assert_eq!(id.slug, "icv1-fp10-t5-192cb216");
}

#[test]
fn a_subgroup_order_is_refused() {
    // 549756390943 is the prime subgroup of K_0 over GF(2^41), not its group.
    let f41 = (BigUint::one() << 41u32) | BigUint::from(0b1001u32);
    assert!(koblitz(0, 41, &f41, &BigUint::from(549_756_390_943u64)).is_none());
}

#[test]
fn the_constructors_name_curves_by_the_registry() {
    // Built curves carry the slug the registry holds for them, and every
    // legacy spelling a frozen report used resolves to that slug.
    let kc = KoblitzCurve::new(0, 41).unwrap();
    assert_eq!(kc.label(), "icv1-f2m41-tm2308219-7f48b14a");
    for legacy in ["K_0 / GF(2^41)", "K_0/GF(2^41)", "k0n41", "`K_0/2^41`", "K₀/GF(2^41)"] {
        assert_eq!(resolve(legacy), Some("icv1-f2m41-tm2308219-7f48b14a"), "{legacy}");
    }
    let inst = koblitz_instance(1, 19).unwrap();
    assert!(same_curve(&inst.name, "K_1 / GF(2^19)"));
    assert_eq!(inst.curve_id().slug, inst.name);

    let bench = roster_prime_instance(10).unwrap();
    assert_eq!(bench.name, "icv1-fp10-t5-192cb216");
    assert!(same_curve(&bench.name, "bench-10bit"));
    assert!(!same_curve(&bench.name, "bench-12bit"));

    let binary = random_binary_instance(15, 123_212_651_130, 8).unwrap();
    assert!(same_curve(&binary.name, "random-binary-n15-b524b"));
}

#[test]
fn names_fold_as_the_reference_folds_them() {
    assert_eq!(normalise(" `K₁ / GF(2^{11})` "), "k1/gf(2^11)");
    assert_eq!(normalise("K_1/GF(2^11)"), normalise("K₁/GF(2^11)"));
    assert_eq!(binary_degree("icv1-f2m27-t12543-569dca8b"), Some(27));
    assert_eq!(binary_degree("random-binary-n13-b1503"), Some(13));
    assert_eq!(binary_degree("K_0 / GF(2^41)"), Some(41));
    assert_eq!(resolve("no such curve"), None);
}
