use endomorphism_search::{
    curves::{self, Counts, Curve, Map, Operation},
    orders::{self, Limits},
    scalar,
};
use serde_json::{json, Value};
use std::process::Command;

#[test]
fn exact_enumeration_against_rectangle() {
    for d in [-3, -4, -7, -8, -11, -12, -16, -27, -619] {
        for bound in [1, 2, 10, 25, 155] {
            let actual = orders::enumerate(d, bound, true, Limits::default()).unwrap();
            let mut expected = Vec::new();
            let radius = 2 * orders::isqrt(bound as u64) as i64 + 2;
            for a in -radius..=radius {
                for b in -radius..=radius {
                    let n = orders::norm(d, a, b);
                    if b != 0 && n > 0 && n <= bound as i128 {
                        expected.push((a, b, n as i64));
                    }
                }
            }
            let mut observed: Vec<_> = actual["candidates"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| {
                    (
                        v["a"].as_i64().unwrap(),
                        v["b"].as_i64().unwrap(),
                        v["degree"].as_i64().unwrap(),
                    )
                })
                .collect();
            expected.sort();
            observed.sort();
            assert_eq!(expected, observed, "{d} {bound}");
            assert_eq!(actual["status"], "complete");
        }
    }
}
#[test]
fn order_619_classes_are_not_self_maps() {
    let v = orders::scan(-619, 1000, Limits::default()).unwrap();
    assert_eq!(v["order"]["minimum_non_scalar_degree"], 155);
    assert_eq!(v["order"]["geometric_unit_count"], 2);
    assert_eq!(v["order"]["class_forms"]["class_number"], 5);
    let leading: Vec<_> = v["order"]["class_forms"]["forms"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v["a"].as_i64().unwrap())
        .collect();
    assert_eq!(leading, [1, 5, 5, 7, 7]);
    assert!(!v["norm_search"]["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .any(|v| v["degree"] == 5 || v["degree"] == 7));
    assert_eq!(v["curve_binding"]["status"], "not_bound_to_a_curve");
}
#[test]
fn nonmaximal_units() {
    for (d, count) in [(-3, 6), (-4, 4), (-12, 2), (-16, 2)] {
        assert_eq!(
            orders::scan(d, 10, Limits::default()).unwrap()["order"]["geometric_unit_count"],
            count
        );
    }
}
#[test]
fn conductor_decomposition() {
    for (d, pair) in [
        (-619, (-619, 1)),
        (-112, (-7, 4)),
        (-12, (-3, 2)),
        (-16, (-4, 2)),
        (-27, (-3, 3)),
    ] {
        assert_eq!(orders::fundamental(d).unwrap(), pair);
    }
}
#[test]
fn incomplete_does_not_certify_absence() {
    let v = orders::scan(
        -3,
        10000,
        Limits {
            work: 1,
            ..Limits::default()
        },
    )
    .unwrap();
    assert_eq!(v["norm_search"]["status"], "incomplete");
    assert_eq!(
        v["decision_tree"]["scalar_multiplication"]["status"],
        "search_incomplete"
    );
    assert!(orders::reduced_forms(-619, 1).unwrap()["class_number"].is_null());
    let v = orders::enumerate(
        -3,
        1000,
        true,
        Limits {
            candidates: 1,
            ..Limits::default()
        },
    )
    .unwrap();
    assert_eq!(v["stop_reason"], "candidate_limit");
}
#[test]
fn rejects_invalid_and_overflowing_inputs() {
    for d in [0, 1, -1, -2, -5, i64::MIN] {
        assert!(orders::scan(d, 1000, Limits::default()).is_err());
    }
    assert!(orders::enumerate(-3, 0, true, Limits::default()).is_err());
    assert!(orders::enumerate(-3, i64::MAX, true, Limits::default()).is_err());
    for seconds in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(Limits {
            seconds,
            ..Limits::default()
        }
        .check()
        .is_err());
    }
    for p in [2, 3, 4, 9, 16383] {
        assert!(Curve::new(p, 2, 2).is_err());
    }
    assert!(Curve::new(17, 0, 0).is_err());
}
#[test]
fn integer_square_root_boundaries() {
    for n in 0..10000 {
        let s = orders::isqrt(n);
        assert!(s * s <= n && (s + 1) * (s + 1) > n);
    }
    for n in [u64::MAX, 1u64 << 63, (1u64 << 63) - 1] {
        let s = orders::isqrt(n) as u128;
        assert!(s * s <= n as u128 && (s + 1) * (s + 1) > n as u128);
    }
}
#[test]
fn exhaustive_group_law() {
    let c = Curve::new(17, 2, 2).unwrap();
    let points = c.points();
    assert_eq!(points.len(), 19);
    for p in &points {
        assert_eq!(c.add(*p, c.neg(*p)), None);
        assert_eq!(c.mul(19, *p), None);
        for q in &points {
            assert!(c.contains(c.add(*p, *q)));
            for r in &points {
                assert_eq!(c.add(c.add(*p, *q), *r), c.add(*p, c.add(*q, *r)));
            }
        }
    }
}
fn quotient(c: Curve, kernel: Vec<curves::Point>) -> (Curve, Map) {
    let target = curves::velu_target(c, &kernel).unwrap();
    let m = Map {
        id: "test".into(),
        degree: kernel.len() as i64,
        op: Operation::Velu {
            kernel,
            target,
            scale: 1,
        },
        traces: vec![],
        eigenvalue: None,
        verification: Value::Null,
    };
    (target, m)
}
#[test]
fn exhaustive_velu_homomorphism() {
    let c = Curve::new(29, 4, 22).unwrap();
    let pts = c.points();
    let g = *pts[1..].iter().find(|p| c.mul(2, **p).is_none()).unwrap();
    let (target, m) = quotient(c, vec![None, g]);
    let mut cost = Counts::default();
    assert_eq!(target.points().len(), pts.len());
    assert_eq!(
        pts.iter()
            .filter(|p| m.evaluate(c, **p, &mut cost).is_none())
            .count(),
        2
    );
    for p in &pts {
        assert!(target.contains(m.evaluate(c, *p, &mut cost)));
        for q in &pts {
            assert_eq!(
                m.evaluate(c, c.add(*p, *q), &mut cost),
                target.add(m.evaluate(c, *p, &mut cost), m.evaluate(c, *q, &mut cost))
            );
        }
    }
}
#[test]
fn sage_published_degree_seven_fixture() {
    let c = Curve::new(11, 1, 1).unwrap();
    let kernel = (0..7).map(|i| c.mul(i, Some((6, 5)))).collect();
    let (target, m) = quotient(c, kernel);
    assert_eq!((target.a, target.b), (7, 8));
    assert_eq!(
        m.evaluate(c, Some((4, 5)), &mut Counts::default()),
        Some((10, 0))
    );
}
#[test]
fn rejects_nonclosed_and_duplicate_kernels() {
    let c = Curve::new(17, 2, 2).unwrap();
    let p = c.points()[1];
    assert!(curves::velu_target(c, &[None, p]).is_err());
    assert!(curves::velu_target(c, &[None, None, p]).is_err());
    assert!(curves::velu_target(c, &[None, Some((1, 1))]).is_err());
}
#[test]
fn degree_two_maps_certify_conductor() {
    let v = curves::probe(Curve::new(29, 16, 2).unwrap(), 1000, 32, 30.0, false).unwrap();
    assert_eq!(
        v["ring_binding"]["candidate_discriminants_before_map_search"],
        json!([-7, -28, -112])
    );
    assert_eq!(
        v["ring_binding"]["candidate_discriminants_after_map_search"],
        json!([-7])
    );
    assert_eq!(v["ring_binding"]["full_ring_discriminant"], -7);
    assert!(v["explicit_map_search"]["maps"]
        .as_array()
        .unwrap()
        .iter()
        .any(|m| m["geometric_degree"] == 2));
}
#[test]
fn no_ring_upgrade_without_map_evidence() {
    let v = curves::probe(Curve::new(29, 16, 2).unwrap(), 1, 32, 30.0, false).unwrap();
    assert_eq!(v["ring_binding"]["status"], "unresolved");
    assert!(v["ring_binding"]["full_ring_discriminant"].is_null());
    assert_eq!(v["explicit_map_search"]["maps_constructed_and_verified"], 0);
}
#[test]
fn frobenius_generator_and_evidence_tree() {
    let v = curves::probe(Curve::new(167, 25, 36).unwrap(), 1000, 32, 30.0, true).unwrap();
    assert_eq!(v["trace"], 7);
    assert_eq!(v["subgroup"]["order"], 23);
    assert_eq!(v["ring_binding"]["full_ring_discriminant"], -619);
    let m = &v["explicit_map_search"]["maps"][0];
    assert_eq!(m["geometric_degree"], 155);
    assert_eq!(m["subgroup_eigenvalue"], 21);
    assert_eq!(v["full_field_frobenius_action"]["eigenvalue"], 1);
    assert!(v["speedup"].is_null());
    assert_eq!(
        v["decision_tree"]["scalar_multiplication"]["status"],
        "verified_map_candidates"
    );
    assert!(v["decision_tree"]["scalar_multiplication"]["measured_gain_scope"].is_null());
    assert_eq!(
        v["decision_tree"]["rho"]["extra_geometric_units"],
        "excluded_for_certified_ring"
    );
    assert_eq!(v["decision_tree"]["rho"]["measured_gain"], "unresolved");
    assert_eq!(v["index_calculus"]["status"], "unresolved");
}
#[test]
fn symmetries_count_overlaps_once() {
    assert_eq!(
        curves::symmetries(&[2, 4, 16], 31)["distinct_action_count_including_negation"],
        10
    );
    assert_eq!(
        curves::symmetries(&[0, 1, 30], 31)["distinct_action_count_including_negation"],
        2
    );
}
#[test]
fn split_including_half_integer_ties() {
    for n in [7, 13, 19, 23, 31, 109] {
        for v in 0..n {
            let basis = scalar::lattice(n, v, &mut Counts::default()).unwrap();
            for k in -n..2 * n {
                let (a, b) = scalar::split(k, n, v, basis, &mut Counts::default()).unwrap();
                assert_eq!((a + v * b - k).rem_euclid(n), 0);
            }
        }
    }
}
#[test]
fn signed_joint_arithmetic() {
    let c = Curve::new(17, 2, 2).unwrap();
    let p = c.points()[1];
    let q = c.mul(7, p);
    for a in -10..=10 {
        for b in -10..=10 {
            assert_eq!(
                scalar::joint(c, a, p, b, q, &mut Counts::default()),
                c.add(c.mul(a, p), c.mul(b, q))
            );
        }
    }
}

#[test]
fn malformed_lattices_are_rejected() {
    for (n, basis) in [
        (0, ((0, 0), (0, 0))),
        (7, ((0, 0), (0, 0))),
        (7, ((7, 0), (2, 1))),
        (7, ((i64::MAX, 0), (0, 1))),
    ] {
        assert!(scalar::split(1, n, 3, basis, &mut Counts::default()).is_err());
    }
}
#[test]
fn counts_include_map_and_per_call_tables() {
    let c = Curve::new(17, 2, 2).unwrap();
    let p = c.points()[1];
    let mut measured = Counts::default();
    assert_eq!(c.wnaf(7, p, 4, &mut measured), c.mul(7, p));
    assert!(measured.group_calls >= 4);
    assert!(measured.field_inverse > 0);
    let a = curves::probe(Curve::new(167, 25, 36).unwrap(), 1000, 32, 30.0, true).unwrap();
    assert_eq!(
        a["scalar_operation_diagnostic"]["maps"][0]["exhaustive_scalar_checks"],
        23
    );
    assert!(a["scalar_operation_diagnostic"]["maps"][0]["speedup"].is_null());
}
fn cli(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_endomorphism-search"))
        .args(args)
        .output()
        .unwrap()
}
#[test]
fn cli_conditional_scan_and_partial_sweep() {
    let p = cli(&["scan", "--discriminant=-619", "--degree-bound", "154"]);
    assert!(p.status.success());
    let v: Value = serde_json::from_slice(&p.stdout).unwrap();
    assert_eq!(v["norm_search"]["candidate_count"], 0);
    assert_eq!(v["curve_binding"]["status"], "not_bound_to_a_curve");
    let p = cli(&["sweep", "--max-orders", "1"]);
    assert!(p.status.success());
    let v: Value = serde_json::from_slice(&p.stdout).unwrap();
    assert_eq!(v["orders_scanned"], 1);
    assert_eq!(v["status"], "incomplete_range");
    assert_eq!(v["stop_reason"], "order_limit");
}
#[test]
fn cli_rejections_emit_no_result() {
    for args in [
        vec!["scan", "--discriminant=-5"],
        vec!["scan", "--discriminant=-619", "--seconds", "NaN"],
        vec!["scan", "--discriminant=-619", "--discriminant=-3"],
        vec!["demo", "--count-ops=false"],
        vec!["scan", "--discriminant=-619", "--unknown", "1"],
    ] {
        let p = cli(&args);
        assert!(!p.status.success());
        assert!(p.stdout.is_empty());
    }
}
