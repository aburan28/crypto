use super::*;

fn candidate(discriminant: i64, a: i64, d: i64, b: [i64; 2]) -> Candidate {
    Candidate {
        label: "test".into(),
        discriminant,
        a,
        d,
        b,
        binding: None,
        evidence_refs: vec![],
    }
}
fn input(c: Candidate) -> Vec<u8> {
    json!({"schema_version":"jacobian-candidates/v1","candidates":[c]})
        .to_string()
        .into_bytes()
}
fn h(c: &Candidate, x: Element, y: Element) -> i128 {
    let t = i128::from(c.discriminant).rem_euclid(2);
    let n = (t * t - i128::from(c.discriminant)) / 4;
    let cross = mul(
        mul([x[0] + t * x[1], -x[1]], c.b.map(i128::from), t, n),
        y,
        t,
        n,
    );
    i128::from(c.a) * norm(x, t, n) + i128::from(c.d) * norm(y, t, n) + 2 * cross[0] + t * cross[1]
}

#[test]
fn discriminant_619_certificate() {
    let c = candidate(-619, 12, 13, [0, 1]);
    let r = analyze(&c, 2_000_000);
    assert_eq!(r["minimum"], 12);
    assert_eq!(r["decomposable_minimum_upper_bound"], 7);
    assert_eq!(
        r["reduced_forms"],
        json!([
            [1, 1, 155],
            [5, -1, 31],
            [5, 1, 31],
            [7, -5, 23],
            [7, 5, 23]
        ])
    );
    assert_eq!(r["geometrically_indecomposable"], true);
    let w: [Element; 2] = serde_json::from_value(r["minimum_witness"].clone()).unwrap();
    assert_eq!(h(&c, w[0], w[1]), 12);
    let report = report(&input(c), None, 2_000_000).unwrap();
    assert_eq!(report["records"][0]["status"], "conditional_existence");
    assert!(report["records"][0]["exists_over_base_field"].is_null());
    assert!(report["records"][0]["binding"].is_null());
}

#[test]
fn fundamental_discriminants_and_complete_class_lists() {
    for d in [
        -3, -4, -7, -8, -11, -15, -19, -20, -23, -24, -31, -40, -163, -619,
    ] {
        assert!(fundamental(d));
    }
    for d in [-1, -2, -9, -12, -16, -27, -36, -100, 0, 5] {
        assert!(!fundamental(d));
    }
    assert_eq!(reduced_forms(-23), vec![[1, 1, 6], [2, -1, 3], [2, 1, 3]]);
    assert_eq!(reduced_forms(-20), vec![[1, 0, 5], [2, 2, 3]]);
    assert_eq!(reduced_forms(-4), vec![[1, 0, 1]]);
}

#[test]
fn independent_direct_quadratic_form_crosscheck() {
    for disc in [-3, -4, -7, -8, -11, -19, -20] {
        for a in 1..=4 {
            for d in 1..=4 {
                for u in -2..=2 {
                    for v in -2..=2 {
                        let c = candidate(disc, a, d, [u, v]);
                        let r = analyze(&c, 2_000_000);
                        if r["principal"] != true {
                            continue;
                        }
                        let mut direct = i128::from(a.min(d));
                        for x0 in -4..=4 {
                            for x1 in -4..=4 {
                                for y0 in -4..=4 {
                                    for y1 in -4..=4 {
                                        if [x0, x1, y0, y1] != [0; 4] {
                                            direct = direct.min(h(&c, [x0, x1], [y0, y1]));
                                        }
                                    }
                                }
                            }
                        }
                        assert_eq!(r["minimum"], json!(direct), "{c:?}");
                    }
                }
            }
        }
    }
}

#[test]
fn nonprincipal_classes_and_budget_never_become_negative_existence() {
    let split = analyze(&candidate(-619, 1, 1, [0, 0]), 2_000_000);
    assert_eq!(split["geometrically_indecomposable"], false);
    // A nonprincipal ideal class defeats the naive 'minimum > 1' test.
    let c = candidate(-20, 2, 3, [0, 1]);
    let r = analyze(&c, 2_000_000);
    assert_eq!(r["minimum"], 2);
    assert_eq!(r["status"], "inconclusive");
    assert!(r["geometrically_indecomposable"].is_null());
    let r = report(&input(c), None, 1).unwrap();
    assert_eq!(r["records"][0]["lattice"]["status"], "inconclusive_budget");
    assert!(r["records"][0]["exists_over_base_field"].is_null());
}

#[test]
fn invalid_candidates_and_overflow_inputs_fail_closed() {
    for c in [
        candidate(i64::MIN, 12, 13, [0, 1]),
        candidate(-619, i64::MAX, 13, [0, 1]),
        candidate(-619, 12, 13, [i64::MIN, 1]),
        candidate(-12, 1, 1, [0, 0]),
    ] {
        let r = analyze(&c, 100);
        assert_eq!(r["status"], "unsupported");
        assert!(r["geometrically_indecomposable"].is_null());
    }
    let r = analyze(&candidate(-619, 12, 14, [0, 1]), 2_000_000);
    assert_eq!(r["status"], "not_principal_candidate");
    let mut bad: Value = serde_json::from_slice(&input(candidate(-619, 12, 13, [0, 1]))).unwrap();
    bad["candidates"][0]["endomorphism_verified"] = json!(true);
    assert!(report(bad.to_string().as_bytes(), None, 100).is_err());
}

#[test]
fn exact_catalog_binding_and_sql_escaping() {
    let uid = format!("urn:ec-record:1:sha256:{}", "a".repeat(64));
    let mut c = candidate(-619, 12, 13, [0, 1]);
    c.label = "quoted ' label".into();
    c.binding = Some(Binding {
        slug: "test-model".into(),
        curve_uid: uid.clone(),
    });
    let mut registry = json!({"schema_version":1,"curves":[{"slug":"test-model","model_json":"{}",
        "representations":[{"curve_uid":uid,"ec1":"test-alias","field":{},"curve":{}}]}]});
    let bytes = input(c.clone());
    assert!(report(&bytes, None, 2_000_000).is_err());
    let r = report(&bytes, Some(registry.to_string().as_bytes()), 2_000_000).unwrap();
    assert_eq!(r["records"][0]["binding"]["curve_uid"], uid);
    assert!(sql(&r).contains("quoted '' label"));
    registry["curves"][0]["representations"][0]["curve_uid"] = json!("wrong");
    assert!(report(&bytes, Some(registry.to_string().as_bytes()), 2_000_000).is_err());
}

#[test]
fn deterministic_replay_binds_input_source_and_claim() {
    let bytes = input(candidate(-619, 12, 13, [0, 1]));
    let a = report(&bytes, None, 2_000_000).unwrap();
    assert_eq!(a, report(&bytes, None, 2_000_000).unwrap());
    let b = report(&input(candidate(-619, 12, 14, [0, 1])), None, 2_000_000).unwrap();
    assert_ne!(
        a["records"][0]["certificate_uid"],
        b["records"][0]["certificate_uid"]
    );
    assert_ne!(a["input_blake3"], b["input_blake3"]);
    assert!(!a["checker_source_blake3"].as_str().unwrap().is_empty());
}
