//! Independent arithmetic, frozen-target, orbit, and restart check for v3.
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

fn key(point: Point, p: i64) -> i64 {
    match point {
        None => 0,
        Some((x, y)) => 2 * (x + 1) + i64::from(y > p - y),
    }
}

fn six_class_checks(p: i64, b: i64, r: u64, zeta: i64, g: Point) -> (u64, usize) {
    let roots: Vec<u64> = (1..r).filter(|&v| (v * v + v + 1) % r == 0).collect();
    assert_eq!(roots.len(), 2);
    let lambda: Vec<u64> = roots
        .into_iter()
        .filter(|&v| multiple(v, g, p, b) == cubic_action(g, p, zeta))
        .collect();
    assert_eq!(lambda.len(), 1);
    let lambda = lambda[0];
    let mut representatives = BTreeSet::new();
    let mut points = vec![None];
    for x in 0..p {
        for y in 0..p {
            let point = Some((x, y));
            if on_curve(point, p, b) {
                points.push(point);
            }
        }
    }
    assert_eq!(points.len(), r as usize);
    for point in points.iter().copied() {
        let mut current = point;
        let mut scalar = 1u64;
        let mut orbit = Vec::new();
        for _ in 0..3 {
            orbit.push((current, scalar));
            orbit.push((current.map(|(x, y)| (x, modulo(-y, p))), (r - scalar) % r));
            current = cubic_action(current, p, zeta);
            scalar = scalar * lambda % r;
        }
        let &(representative, multiplier) = orbit
            .iter()
            .min_by_key(|(member, _)| key(*member, p))
            .unwrap();
        assert_eq!(multiple(multiplier, point, p, b), representative);
        for &(member, mu) in &orbit {
            assert_eq!(multiple(mu, point, p, b), member);
            assert!(key(representative, p) <= key(member, p));
        }
        representatives.insert(representative);
    }
    assert_eq!(representatives.len(), 1 + (r as usize - 1) / 6);
    (lambda, points.len())
}

#[derive(Default)]
struct Attempt {
    seed: u64,
    success: bool,
    gae: u64,
    steps: u64,
    walks: u64,
    dp: u64,
    canonicalisations: u64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let raw = fs::read_to_string(args.next().ok_or("verify RAW.tsv OLD.tsv RUN.log")?)?;
    let old = fs::read_to_string(args.next().ok_or("missing old TSV")?)?;
    let run_log = fs::read_to_string(args.next().ok_or("missing run log")?)?;
    let mut frozen = BTreeMap::new();
    for line in old.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if f[4] == "source_relations" {
            let key = (f[0].parse::<u64>()?, f[3].parse::<u64>()?);
            let q = Some((f[5].parse::<i64>()?, f[6].parse::<i64>()?));
            assert!(frozen.insert(key, (q, f[10].parse::<u64>()?)).is_none());
        }
    }
    assert_eq!(frozen.len(), 16);
    let mut attempts: BTreeMap<(u64, u64), Vec<Attempt>> = BTreeMap::new();
    for line in run_log
        .lines()
        .filter(|line| line.starts_with("six_attempt "))
    {
        let parts: BTreeMap<&str, &str> = line
            .split_whitespace()
            .skip(1)
            .map(|part| part.split_once('=').unwrap())
            .collect();
        let get = |name: &str| parts[name].parse::<u64>().unwrap();
        let p = get("p");
        let index = get("index");
        let j = get("attempt");
        let seed = get("seed");
        let expected_seed =
            (0x5248_4f4e_4547_0000 ^ p ^ (index << 16)) ^ j.wrapping_mul(0x9e37_79b9_7f4a_7c15);
        assert_eq!(seed, expected_seed);
        let group = attempts.entry((p, index)).or_default();
        assert_eq!(j as usize, group.len());
        group.push(Attempt {
            seed,
            success: parts["status"] == "verified",
            gae: get("gae"),
            steps: get("steps"),
            walks: get("walks"),
            dp: get("dp"),
            canonicalisations: get("canonicalisations"),
        });
    }
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
    let mut total_restarts = 0;
    let mut six_class_point_checks = 0;
    for (n, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), header.len(), "TSV row {} width", n + 2);
        let number = |name: &str| fields[column(name)].parse::<u64>().unwrap();
        let p = number("p");
        let b = number("b");
        let r = number("r");
        let index = number("target_index");
        let (expected_b, expected_r, g, u) = match p {
            13 => (2, 19, Some((1, 4)), 2),
            199 => (3, 211, Some((1, 2)), 14),
            _ => panic!("unregistered test field"),
        };
        assert_eq!((b, r), (expected_b, expected_r));
        assert!(index < 8);
        let target_b = modulo(-27 * b as i64, p as i64);
        let q = Some((number("qx") as i64, number("qy") as i64));
        let mapped = Some((number("mapped_qx") as i64, number("mapped_qy") as i64));
        let d = number("log");
        assert_eq!((q, d), frozen[&(p, index)]);
        assert!(d < r);
        assert_eq!(fields[column("status")], "verified");
        assert!(on_curve(q, p as i64, b as i64));
        assert_eq!(multiple(d, g, p as i64, b as i64), q);
        assert_eq!(point_count(p as i64, b as i64), r);
        assert_eq!(point_count(p as i64, target_b), r);
        let arm = fields[column("arm")];
        let image = map(q, p as i64, b as i64);
        let (check_g, check_q, check_b, map_cost) = match arm {
            "source_relations" | "rho_negation" => {
                assert_eq!(mapped, image); // Legacy companion image, not the working target.
                (g, q, b as i64, (0, 0, 0, 0))
            }
            "transport_relations" => {
                assert_eq!(mapped, image);
                (map(g, p as i64, b as i64), image, target_b, (2, 10, 4, 2))
            }
            "isomorphism_relations" => {
                assert_eq!(mapped, scale(q, p as i64, u));
                (scale(g, p as i64, u), mapped, target_b, (2, 6, 0, 0))
            }
            "rho_six_class" => {
                assert_eq!(mapped, q);
                (g, q, b as i64, (0, 0, 0, 0))
            }
            _ => panic!("unknown arm"),
        };
        assert!(on_curve(check_q, p as i64, check_b));
        assert_eq!(multiple(d, check_g, p as i64, check_b), check_q);
        assert_eq!(
            (
                number("map_evaluations"),
                number("map_muls"),
                number("map_sqrs"),
                number("map_invs")
            ),
            map_cost
        );
        let phases = ["setup", "collection", "descent", "verification", "rho"];
        let mut total = 0;
        for phase in phases {
            let gae = number(&format!("{phase}_gae"));
            assert_eq!(
                gae,
                number(&format!("{phase}_adds")) + number(&format!("{phase}_doubles"))
            );
            total += gae;
        }
        assert_eq!(number("total_gae"), total);
        if arm == "rho_six_class" {
            let group = &attempts[&(p, index)];
            assert!((1..=16).contains(&group.len()));
            assert_eq!(number("rho_attempts") as usize, group.len());
            assert!(group.last().unwrap().success);
            assert!(group[..group.len() - 1].iter().all(|v| !v.success));
            assert_eq!(group.iter().map(|v| v.gae).sum::<u64>(), number("rho_gae"));
            assert_eq!(
                group.iter().map(|v| v.steps).sum::<u64>(),
                number("rho_steps")
            );
            assert_eq!(
                group.iter().map(|v| v.walks).sum::<u64>(),
                number("rho_walks")
            );
            assert_eq!(group.iter().map(|v| v.dp).sum::<u64>(), number("rho_dp"));
            assert_eq!(
                group.iter().map(|v| v.canonicalisations).sum::<u64>(),
                number("rho_unpriced_canonicalisations")
            );
            let (lambda, checked) =
                six_class_checks(p as i64, b as i64, r, if p == 13 { 3 } else { 92 }, g);
            six_class_point_checks += checked;
            assert_eq!(number("rho_lambda"), lambda);
            assert_eq!(number("rho_class_setup_field_muls"), 1);
            assert_eq!(
                3 * number("rho_canon_field_muls"),
                2 * number("rho_canon_point_negations")
            );
            assert!(number("rho_canon_field_muls") <= 2 * number("rho_unpriced_canonicalisations"));
            total_restarts += group.len() - 1;
            assert_ne!(group[0].seed, 0);
        } else {
            for name in [
                "rho_lambda",
                "rho_attempts",
                "rho_canon_field_muls",
                "rho_canon_point_negations",
                "rho_class_setup_field_muls",
            ] {
                assert_eq!(number(name), 0);
            }
        }
        assert!(answers
            .entry((p, index))
            .or_default()
            .insert(arm.to_string(), d)
            .is_none());
        records += 1;
    }
    assert_eq!(records, 80);
    assert_eq!(answers.len(), 16);
    assert_eq!(attempts.len(), 16);
    let arms: BTreeSet<&str> = [
        "source_relations",
        "transport_relations",
        "isomorphism_relations",
        "rho_negation",
        "rho_six_class",
    ]
    .into_iter()
    .collect();
    for values in answers.values() {
        assert_eq!(
            values.keys().map(String::as_str).collect::<BTreeSet<_>>(),
            arms
        );
        assert!(values.values().all(|&d| d == values["source_relations"]));
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
    println!("independent_arithmetic_check=pass records={records} source_targets={} compared_arms=5 total_restarts={total_restarts} point_counts=4 phase_sums_checked={records} cubic_automorphism_pair_checks={automorphism_pair_checks} scaling_isomorphism_pair_checks={scaling_isomorphism_pair_checks} six_class_point_checks={six_class_point_checks} isogeny_scaling_scalar_13={scalar13} isogeny_scaling_scalar_199={scalar199}", answers.len());
    Ok(())
}
