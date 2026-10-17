//! Independent arithmetic and record check for the final paired TSV.
//! This deliberately does not import the runner's curve or rho implementation.
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;

type Point = Option<(i64, i64)>;

fn modulo(x: i64, p: i64) -> i64 {
    x.rem_euclid(p)
}

fn inverse(a: i64, p: i64) -> i64 {
    let (mut old_r, mut r) = (p, modulo(a, p));
    let (mut old_t, mut t) = (0, 1);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_t, t) = (t, old_t - q * t);
    }
    assert_eq!(old_r, 1);
    modulo(old_t, p)
}

fn add(a: Point, b: Point, p: i64, c: i64) -> Point {
    match (a, b) {
        (None, point) | (point, None) => point,
        (Some((x1, y1)), Some((x2, y2))) => {
            if x1 == x2 && modulo(y1 + y2, p) == 0 {
                return None;
            }
            let slope = if x1 == x2 {
                modulo(3 * x1 * x1 * inverse(2 * y1, p), p)
            } else {
                modulo((y2 - y1) * inverse(x2 - x1, p), p)
            };
            let x3 = modulo(slope * slope - x1 - x2, p);
            let y3 = modulo(slope * (x1 - x3) - y1, p);
            let point = Some((x3, y3));
            assert!(on_curve(point, p, c));
            point
        }
    }
}

fn multiple(mut n: u64, mut point: Point, p: i64, c: i64) -> Point {
    let mut result = None;
    while n > 0 {
        if n & 1 != 0 {
            result = add(result, point, p, c);
        }
        point = add(point, point, p, c);
        n >>= 1;
    }
    result
}

fn on_curve(point: Point, p: i64, c: i64) -> bool {
    match point {
        None => true,
        Some((x, y)) => modulo(y * y - x * x * x - c, p) == 0,
    }
}

fn map(point: Point, p: i64, b: i64) -> Point {
    match point {
        None => None,
        Some((x, y)) => {
            assert_ne!(x, 0);
            let x2 = modulo(x * x, p);
            let x3 = modulo(x2 * x, p);
            Some((
                modulo((x3 + 4 * b) * inverse(x2, p), p),
                modulo(y * (x3 - 8 * b) * inverse(x3, p), p),
            ))
        }
    }
}

fn point_count(p: i64, b: i64) -> u64 {
    let mut count = 1;
    for x in 0..p {
        for y in 0..p {
            if on_curve(Some((x, y)), p, b) {
                count += 1;
            }
        }
    }
    count
}

fn cubic_action(point: Point, p: i64, zeta: i64) -> Point {
    point.map(|(x, y)| (modulo(zeta * x, p), y))
}

fn check_cubic_automorphism(p: i64, b: i64, zeta: i64) -> usize {
    assert_eq!(modulo(zeta * zeta * zeta, p), 1);
    assert_ne!(zeta, 1);
    let mut points = vec![None];
    for x in 0..p {
        for y in 0..p {
            let point = Some((x, y));
            if on_curve(point, p, b) {
                points.push(point);
            }
        }
    }
    for &point in &points {
        let image = cubic_action(point, p, zeta);
        assert!(on_curve(image, p, b));
        assert_eq!(
            cubic_action(cubic_action(cubic_action(point, p, zeta), p, zeta), p, zeta),
            point
        );
        if point.is_some() {
            assert_ne!(image, point);
        }
    }
    for &left in &points {
        for &right in &points {
            assert_eq!(
                cubic_action(add(left, right, p, b), p, zeta),
                add(
                    cubic_action(left, p, zeta),
                    cubic_action(right, p, zeta),
                    p,
                    b
                )
            );
        }
    }
    points.len() * points.len()
}

fn scale(point: Point, p: i64, u: i64) -> Point {
    point.map(|(x, y)| (modulo(u * u * x, p), modulo(u * u * u * y, p)))
}

fn check_scaling_isomorphism(p: i64, b: i64, target_b: i64, u: i64) -> usize {
    let u2 = modulo(u * u, p);
    let u3 = modulo(u2 * u, p);
    assert_eq!(modulo(u3 * u3 * b, p), target_b);
    let mut points = vec![None];
    for x in 0..p {
        for y in 0..p {
            let point = Some((x, y));
            if on_curve(point, p, b) {
                points.push(point);
            }
        }
    }
    let images: BTreeSet<Point> = points.iter().map(|&point| scale(point, p, u)).collect();
    assert_eq!(images.len(), points.len());
    assert_eq!(point_count(p, target_b) as usize, images.len());
    for &point in &images {
        assert!(on_curve(point, p, target_b));
    }
    for &left in &points {
        for &right in &points {
            assert_eq!(
                scale(add(left, right, p, b), p, u),
                add(scale(left, p, u), scale(right, p, u), p, target_b)
            );
        }
    }
    points.len() * points.len()
}

fn isogeny_scaling_scalar(p: i64, b: i64, target_b: i64, r: u64, u: i64, g: Point) -> u64 {
    let phi_g = map(g, p, b);
    let m = (1..r)
        .find(|&candidate| scale(multiple(candidate, g, p, b), p, u) == phi_g)
        .expect("isogeny and scaling maps agree on a scalar multiple of G");
    for x in 0..p {
        for y in 0..p {
            let point = Some((x, y));
            if on_curve(point, p, b) {
                let image = map(point, p, b);
                assert!(on_curve(image, p, target_b));
                assert_eq!(image, scale(multiple(m, point, p, b), p, u));
            }
        }
    }
    m
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = env::args().nth(1).ok_or("verify RAW.tsv")?;
    let raw = fs::read_to_string(path)?;
    let mut lines = raw.lines();
    let header: Vec<&str> = lines.next().ok_or("empty TSV")?.split('\t').collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|&v| v == name)
            .expect("required column")
    };
    let mut answers: BTreeMap<(u64, u64), BTreeMap<String, u64>> = BTreeMap::new();
    let mut records = 0;
    for (n, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), header.len(), "TSV row {} width", n + 2);
        let number = |key: &str| fields[column(key)].parse::<u64>().unwrap();
        let p = number("p");
        let b = number("b");
        let r = number("r");
        let index = number("target_index");
        let (expected_b, expected_r, g) = match p {
            13 => (2, 19, Some((1, 4))),
            199 => (3, 211, Some((1, 2))),
            _ => panic!("unregistered test field"),
        };
        assert_eq!((b, r), (expected_b, expected_r));
        assert!(index < 8);
        let target_b = modulo(-27 * b as i64, p as i64);
        let q = Some((number("qx") as i64, number("qy") as i64));
        let q_image = Some((number("phiqx") as i64, number("phiqy") as i64));
        assert!(on_curve(q, p as i64, b as i64));
        assert!(on_curve(q_image, p as i64, target_b));
        assert_eq!(map(q, p as i64, b as i64), q_image);
        assert_eq!(point_count(p as i64, b as i64), r);
        assert_eq!(point_count(p as i64, target_b), r);
        let arm = fields[column("arm")];
        let d = number("log");
        assert!(d < r);
        assert_eq!(fields[column("status")], "verified");
        let (check_g, check_q, check_b) = if arm == "transport_relations" {
            (map(g, p as i64, b as i64), q_image, target_b)
        } else {
            assert!(arm == "source_relations" || arm == "rho_negation");
            (g, q, b as i64)
        };
        assert_eq!(multiple(d, check_g, p as i64, check_b), check_q);
        let phases = number("setup_gae")
            + number("collection_gae")
            + number("descent_gae")
            + number("verification_gae")
            + number("rho_gae");
        assert_eq!(number("total_gae"), phases);
        if arm == "transport_relations" {
            assert_eq!(number("map_evaluations"), 2);
            assert_eq!(
                (number("map_muls"), number("map_sqrs"), number("map_invs")),
                (10, 4, 2)
            );
        } else {
            assert_eq!(number("map_evaluations"), 0);
        }
        let old = answers
            .entry((p, index))
            .or_default()
            .insert(arm.to_string(), d);
        assert!(old.is_none(), "duplicate arm row");
        records += 1;
    }
    assert_eq!(records, 48);
    assert_eq!(answers.len(), 16);
    let arms: BTreeSet<&str> = ["source_relations", "transport_relations", "rho_negation"]
        .into_iter()
        .collect();
    for ((p, index), values) in &answers {
        assert_eq!(
            values.keys().map(String::as_str).collect::<BTreeSet<_>>(),
            arms
        );
        assert_eq!(
            values["source_relations"], values["transport_relations"],
            "p={p} target={index}"
        );
        assert_eq!(
            values["source_relations"], values["rho_negation"],
            "p={p} target={index}"
        );
    }
    let mut automorphism_pair_checks = 0;
    for (p, b, zeta) in [(13, 2, 3), (13, 11, 3), (199, 3, 92), (199, 118, 92)] {
        automorphism_pair_checks += check_cubic_automorphism(p, b, zeta);
    }
    let mut scaling_isomorphism_pair_checks = 0;
    for (p, b, target_b, u) in [(13, 2, 11, 2), (199, 3, 118, 14)] {
        scaling_isomorphism_pair_checks += check_scaling_isomorphism(p, b, target_b, u);
    }
    let scalar13 = isogeny_scaling_scalar(13, 2, 11, 19, 2, Some((1, 4)));
    let scalar199 = isogeny_scaling_scalar(199, 3, 118, 211, 14, Some((1, 2)));
    println!("independent_arithmetic_check=pass records={records} source_targets={} compared_arms=3 point_counts=4 map_evaluations_checked=48 phase_sums_checked={records} cubic_automorphism_pair_checks={automorphism_pair_checks} scaling_isomorphism_pair_checks={scaling_isomorphism_pair_checks} isogeny_scaling_scalar_13={scalar13} isogeny_scaling_scalar_199={scalar199}", answers.len());
    Ok(())
}
