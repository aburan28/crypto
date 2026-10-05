//! Cross-checks between independent routes to the same object.
//!
//! - Over small fields, kernels built by brute force from rational
//!   `ℓ`-torsion points (plain `u64` arithmetic, not [`Field`]) must pass
//!   the verifier, land on a root of `Φ_ℓ(j, Y)`, and be reproduced by the
//!   Elkies construction from that root.
//! - On P-256, the number of roots of `Φ_ℓ(j, Y)` must equal what the
//!   splitting of `ℓ` in `Z[π]` predicts.
//! - The verifier must reject a perturbed kernel and a kernel mixing two
//!   subgroups.
//! - A short P-224 walk must re-verify from its own records.

use num_bigint::BigUint;

use super::curve::{self, Model, OrderAudit};
use super::field::Field;
use super::kernel::{self, EdgeError};
use super::modpoly;
use super::poly;
use super::walk::{self, RouteIndex, StartCurve, Walk, WalkConfig};

type Pt = Option<(u64, u64)>;

fn inv_mod(a: u64, q: u64) -> u64 {
    let mut r = 1u64;
    let (mut b, mut e) = (a % q, q - 2);
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % q;
        }
        b = b * b % q;
        e >>= 1;
    }
    r
}

fn add(q: u64, a: u64, p1: Pt, p2: Pt) -> Pt {
    let (Some((x1, y1)), Some((x2, y2))) = (p1, p2) else {
        return p1.or(p2);
    };
    let lam = if x1 == x2 {
        if (y1 + y2) % q == 0 {
            return None;
        }
        (3 * x1 % q * x1 % q + a) % q * inv_mod(2 * y1 % q, q) % q
    } else {
        (y2 + q - y1) % q * inv_mod((x2 + q - x1) % q, q) % q
    };
    let x3 = (lam * lam % q + 2 * q - x1 - x2) % q;
    let y3 = (lam * ((x1 + q - x3) % q) % q + q - y1) % q;
    Some((x3, y3))
}

fn mul(q: u64, a: u64, k: u64, p: Pt) -> Pt {
    (0..k).fold(None, |acc, _| add(q, a, acc, p))
}

fn points(q: u64, a: u64, b: u64) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    for x in 0..q {
        let r = (x * x % q * x + a * x + b) % q;
        for y in 0..q {
            if y * y % q == r {
                out.push((x, y));
            }
        }
    }
    out
}

/// The kernel polynomial of `<P>` for a point of odd prime order `ℓ`.
fn kernel_of(f: &Field, q: u64, a: u64, ell: u64, p: (u64, u64)) -> poly::Poly {
    let mut h = vec![f.one()];
    for k in 1..=(ell - 1) / 2 {
        let (x, _) = mul(q, a, k, Some(p)).unwrap();
        h = poly::mul(f, &h, &vec![f.neg(&f.from_u64(x)), f.one()]);
    }
    h
}

#[test]
fn rational_torsion_kernels_agree_with_phi_and_elkies() {
    let q = 1009u64;
    let f = Field::new(&BigUint::from(q)).unwrap();
    for ell in [3u64, 5, 7] {
        let phi = modpoly::modular_polynomial(&f, ell).unwrap();
        let (mut seen, mut elkies_checked) = (0, 0);
        'curves: for a in 1..40u64 {
            for b in 1..40u64 {
                let m = Model {
                    a: f.from_u64(a),
                    b: f.from_u64(b),
                };
                let Some(j) = m.j(&f) else { continue };
                if f.is_zero(&j) || j == f.from_u64(1728) {
                    continue;
                }
                let pts = points(q, a, b);
                let Some(&p) = pts.iter().find(|&&p| mul(q, a, ell, Some(p)).is_none()) else {
                    continue;
                };
                let h = kernel_of(&f, q, a, ell, p);
                let cod =
                    kernel::verify_kernel(&f, &m, ell, &h).expect("brute-force kernel verifies");
                let j2 = cod.j(&f).unwrap();
                let fy = phi.at_x(&f, &j);
                assert!(
                    f.is_zero(&poly::eval(&f, &fy, &j2)),
                    "ℓ={ell} a={a} b={b}: j' is a root"
                );
                seen += 1;
                if f.is_zero(&j2) || j2 == f.from_u64(1728) {
                    continue;
                }
                match kernel::elkies_isogeny(&f, &m, &phi, &j2) {
                    Ok(iso) => {
                        assert_eq!(iso.kernel, h, "ℓ={ell} a={a} b={b}: Elkies kernel");
                        assert_eq!(iso.codomain, cod);
                        elkies_checked += 1;
                    }
                    Err(EdgeError::MultipleRoot) | Err(EdgeError::SpecialJ) => {}
                    Err(e) => panic!("ℓ={ell} a={a} b={b}: {e:?}"),
                }
                if elkies_checked >= 12 {
                    break 'curves;
                }
            }
        }
        assert!(
            seen >= 12 && elkies_checked >= 8,
            "ℓ={ell}: {seen} kernels, {elkies_checked} Elkies"
        );
    }
}

#[test]
fn p256_root_counts_match_the_splitting_of_ell() {
    let start = StartCurve::p256();
    let primes = [3u64, 5, 7, 11, 13];
    let class = walk::ClassInfo::compute(&start, &primes);
    let f = Field::new(&start.p).unwrap();
    let m = Model {
        a: f.from_big(&start.a),
        b: f.from_big(&start.b),
    };
    let j = m.j(&f).unwrap();
    for &l in &primes {
        let phi = modpoly::modular_polynomial(&f, l).unwrap();
        let roots = poly::roots(&f, &phi.at_x(&f, &j), 9);
        assert_eq!(Some(roots.len()), class.expected_neighbours(l), "ℓ = {l}");
        for r in roots {
            let iso = kernel::elkies_isogeny(&f, &m, &phi, &r).unwrap();
            let (canon, _) = curve::canonical_model(&f, &iso.codomain).unwrap();
            assert_eq!(
                curve::audit_order(&f, &canon, &start.order, true, 2, 0),
                OrderAudit::ProvedPrime
            );
        }
    }
}

#[test]
fn verifier_rejects_a_perturbed_and_a_mixed_kernel() {
    // A curve over F_211 (211 ≡ 1 mod 5) with all of E[5] rational.
    let q = 211u64;
    let f = Field::new(&BigUint::from(q)).unwrap();
    let mut found = None;
    'search: for a in 1..q {
        for b in 1..q {
            if (4 * a * a % q * a + 27 * b * b).is_multiple_of(q) {
                continue;
            }
            let pts = points(q, a, b);
            if !(pts.len() + 1).is_multiple_of(25) {
                continue;
            }
            let tors: Vec<_> = pts
                .iter()
                .copied()
                .filter(|&p| mul(q, a, 5, Some(p)).is_none())
                .collect();
            if tors.len() == 24 {
                found = Some((a, b, tors));
                break 'search;
            }
        }
    }
    let (a, b, tors) = found.expect("a curve with full rational 5-torsion");
    let m = Model {
        a: f.from_u64(a),
        b: f.from_u64(b),
    };
    let p1 = tors[0];
    let p2 = *tors
        .iter()
        .find(|&&p| (1..5).all(|k| mul(q, a, k, Some(p1)) != Some(p)))
        .unwrap();
    let h1 = kernel_of(&f, q, a, 5, p1);
    kernel::verify_kernel(&f, &m, 5, &h1).unwrap();
    let mut bad = h1.clone();
    bad[0] = f.add(&bad[0], &f.one());
    assert!(kernel::verify_kernel(&f, &m, 5, &bad).is_err());
    let mixed = poly::mul(
        &f,
        &vec![f.neg(&f.from_u64(p1.0)), f.one()],
        &vec![f.neg(&f.from_u64(p2.0)), f.one()],
    );
    assert_eq!(
        kernel::verify_kernel(&f, &m, 5, &mixed),
        Err(EdgeError::NotSubgroup)
    );
}

#[test]
fn a_short_p224_walk_reverifies_from_its_records() {
    let config = WalkConfig {
        primes: vec![3, 5, 7],
        max_curves: 12,
        ..WalkConfig::default()
    };
    let mut w = Walk::new(StartCurve::p224(), config).unwrap();
    w.run();
    assert!(w.failures.is_empty(), "{:?}", w.failures);
    assert!(w.nodes.len() >= 6);
    assert!(w
        .nodes
        .iter()
        .all(|n| n.audit == Some(OrderAudit::ProvedPrime)));
    let ids: Vec<_> = (0..w.nodes.len()).map(|i| w.node_ids(i)).collect();
    assert_eq!(ids[0].ec1.split('h').next(), Some("EC1P224Cp224"));
    let routes = RouteIndex::build(&w, &ids);
    let json = w.routes_json(&ids, &routes).json();
    let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();
    let (nodes, edges) = walk::verify_routes(&parsed, &w.start, 1).unwrap();
    assert_eq!((nodes, edges), (w.nodes.len(), w.edges.len()));
    // A tampered kernel fails the replay.
    let mut tampered = parsed.clone();
    tampered["edges"][0]["kernel_polynomial"][0] = serde_json::Value::String("1".into());
    assert!(walk::verify_routes(&tampered, &w.start, 1).is_err());
    // So do a relabelled direction and a forged route id.
    let mut tampered = parsed.clone();
    let d = tampered["edges"][0]["direction"]
        .as_str()
        .unwrap()
        .to_string();
    tampered["edges"][0]["direction"] =
        serde_json::Value::String(if d == "h" { "d" } else { "h" }.into());
    assert!(walk::verify_routes(&tampered, &w.start, 1).is_err());
    let mut tampered = parsed.clone();
    tampered["routes"][0]["id"] = serde_json::Value::String("IW1E3h1h000000000000".into());
    assert!(walk::verify_routes(&tampered, &w.start, 1).is_err());
    // P-224 sits on the surface of its depth-1 3-volcano: four rational
    // 3-isogenies, one horizontal and three descending to the floor.
    assert_eq!(w.class.depth(3), 1);
    assert!(w.nodes[0].roots.contains(&(3, 4)));
    assert_eq!(
        w.nodes[0].levels.iter().find(|x| x.0 == 3).map(|x| x.1),
        Some(0)
    );
    let from_root: Vec<char> = w
        .edges
        .iter()
        .filter(|e| e.source == 0 && e.ell == 3)
        .map(|e| e.dir)
        .collect();
    assert_eq!(from_root.iter().filter(|&&c| c == 'h').count(), 1);
    assert_eq!(
        from_root.iter().filter(|&&c| c == 'd').count(),
        3,
        "{from_root:?}"
    );
    let yaml = w.curves_yaml(&ids, &routes);
    assert!(yaml.contains("schema_version: 1\nidentity_rule: sha256_sorted_key_compact_utf8_json_of_field_and_curve\ncurves:\n"));
    assert!(yaml.contains(&format!("  {}:\n", ids[1].slug)));
}

#[test]
fn detectors_class_audits_and_queued_shards_agree() {
    use super::queue;
    let config = WalkConfig {
        primes: vec![3, 11],
        max_curves: 12,
        ..WalkConfig::default()
    };
    let mut w = Walk::new(StartCurve::p224(), config).unwrap();
    w.run();
    w.run_class_audits(3);
    let audits = w.class_audits.clone().unwrap().json();
    let audits: serde_json::Value = serde_json::from_str(&audits).unwrap();
    assert_eq!(audits["ecc_safety"]["all_pass"], true);
    assert_eq!(audits["invariance_check"]["curves_sampled"], 3);
    assert_eq!(audits["invariance_check"]["identical_to_root"], true);
    assert_eq!(audits["pkm"]["solinas_weight"], 3);

    let ids: Vec<_> = (0..w.nodes.len()).map(|i| w.node_ids(i)).collect();
    let routes = RouteIndex::build(&w, &ids);
    let text = w.routes_json(&ids, &routes).json();
    // The root's detectors: P-224 is an a = -3 model with a valid generator.
    let one = queue::run_shard(&text, &w.start, 0, 1, 0).unwrap();
    let root: serde_json::Value = serde_json::from_str(one.jsonl.lines().next().unwrap()).unwrap();
    assert_eq!(root["traits"]["a_minus_3_model"]["value"], true);
    assert_eq!(root["traits"]["generator_valid"]["value"], true);
    assert_eq!(root["traits"]["non_singular"]["value"], true);
    assert_eq!(root["traits"]["coefficient_bits"]["value"]["a"], 2);
    // The YAML records the same detectors.
    let yaml = w.curves_yaml(&ids, &routes);
    for name in [
        "generator_valid",
        "non_singular",
        "coefficient_bits",
        "qr_prefix_64",
    ] {
        assert!(yaml.contains(&format!("      {name}:\n")), "{name}");
    }

    // Three shards collect into exactly the single-shard output.
    let base = std::env::temp_dir().join(format!("isogeny-walk-queue-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&base);
    let dirs: Vec<_> = (0..3)
        .map(|i| {
            let d = base.join(format!("shard{i}"));
            let out = queue::run_shard(&text, &w.start, i, 3, 2).unwrap();
            queue::write_shard(&d, &out).unwrap();
            d
        })
        .collect();
    let (merged, summary, class) = queue::collect(&dirs).unwrap();
    assert_eq!(merged, one.jsonl);
    assert!(class.is_some());
    assert!(summary
        .json()
        .contains(&format!("\"curves\": {}", w.nodes.len())));
    // A missing shard, a repeated shard and an edited shard are refused.
    assert!(queue::collect(&dirs[..2]).is_err());
    assert!(queue::collect(&[
        dirs[0].clone(),
        dirs[1].clone(),
        dirs[1].clone(),
        dirs[2].clone()
    ])
    .is_err());
    let f = dirs[2].join(queue::TRAITS_FILE);
    let edited = std::fs::read_to_string(&f)
        .unwrap()
        .replace("true", "false");
    std::fs::write(&f, edited).unwrap();
    assert!(queue::collect(&dirs).is_err());
    let _ = std::fs::remove_dir_all(&base);

    // Plans: one spec per shard, distinct keys, full shas only.
    let pw = queue::PlanWalk {
        curve_args: vec!["--curve".into(), "p224".into()],
        walk_args: vec!["--primes".into(), "3,11".into()],
    };
    assert!(queue::plan(&pw, "main", "cpu", 2, 60, &queue::PlanStore::default()).is_err());
    let specs = queue::plan(
        &pw,
        &"a".repeat(40),
        "cpu",
        2,
        60,
        &queue::PlanStore::default(),
    )
    .unwrap();
    let keys: Vec<String> = specs
        .iter()
        .map(|s| {
            serde_json::from_str::<serde_json::Value>(&s.json()).unwrap()["idempotency_key"]
                .as_str()
                .unwrap()
                .to_string()
        })
        .collect();
    assert_eq!(keys.len(), 2);
    assert_ne!(keys[0], keys[1]);
}
