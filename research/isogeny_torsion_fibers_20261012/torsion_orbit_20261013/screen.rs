//! Exact rational three-torsion orbit support screen. See PROTOCOL.md.
use std::collections::{BTreeMap, BTreeSet};

type Point = Option<(u64, u64)>;

fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn pow_mod(mut a: u64, mut n: u64, p: u64) -> u64 {
    let mut out = 1;
    while n > 0 {
        if n & 1 == 1 {
            out = mul_mod(out, a, p);
        }
        a = mul_mod(a, a, p);
        n >>= 1;
    }
    out
}

fn sub_mod(a: u64, b: u64, p: u64) -> u64 {
    (a + p - b % p) % p
}

fn add(a: Point, b: Point, p: u64) -> Point {
    let (Some((x1, y1)), Some((x2, y2))) = (a, b) else {
        return a.or(b);
    };
    if x1 == x2 && (y1 + y2) % p == 0 {
        return None;
    }
    let slope = if x1 == x2 {
        mul_mod(mul_mod(3, mul_mod(x1, x1, p), p), pow_mod(2 * y1, p - 2, p), p)
    } else {
        mul_mod(sub_mod(y2, y1, p), pow_mod(sub_mod(x2, x1, p), p - 2, p), p)
    };
    let x3 = sub_mod(sub_mod(mul_mod(slope, slope, p), x1, p), x2, p);
    let y3 = sub_mod(mul_mod(slope, sub_mod(x1, x3, p), p), y1, p);
    Some((x3, y3))
}

fn neg(a: Point, p: u64) -> Point {
    a.map(|(x, y)| (x, sub_mod(0, y, p)))
}

fn scale(mut a: Point, mut n: u64, p: u64) -> Point {
    let mut out = None;
    while n > 0 {
        if n & 1 == 1 {
            out = add(out, a, p);
        }
        a = add(a, a, p);
        n >>= 1;
    }
    out
}

fn on_curve(pt: Point, p: u64, b: u64) -> bool {
    let Some((x, y)) = pt else { return true; };
    mul_mod(y, y, p) == (mul_mod(mul_mod(x, x, p), x, p) + b) % p
}

fn phi(pt: Point, p: u64) -> Point {
    let Some((x, y)) = pt else { return None; };
    if x == 0 {
        return None;
    }
    let x2 = mul_mod(x, x, p);
    let x3 = mul_mod(x2, x, p);
    let ix = pow_mod(x, p - 2, p);
    let ix2 = mul_mod(ix, ix, p);
    let ix3 = mul_mod(ix2, ix, p);
    Some((
        mul_mod((x3 + 4) % p, ix2, p),
        mul_mod(mul_mod(y, sub_mod(x3, 8, p), p), ix3, p),
    ))
}

fn roots(p: u64) -> Vec<Option<u64>> {
    let mut out = vec![None::<u64>; p as usize];
    for y in 0..p {
        let square = mul_mod(y, y, p) as usize;
        out[square] = Some(out[square].map_or(y, |old| old.min(y)));
    }
    out
}

fn base(p: u64, wanted: usize) -> Vec<Point> {
    let roots = roots(p);
    let mut out = Vec::new();
    let mut seen = BTreeSet::new();
    for x in 0..p {
        let rhs = (mul_mod(mul_mod(x, x, p), x, p) + 1) % p;
        if let Some(y) = roots[rhs as usize] {
            let projected = scale(Some((x, y)), 6, p);
            if projected.is_some() && seen.insert(projected) {
                out.push(projected);
                if out.len() == wanted {
                    break;
                }
            }
        }
    }
    assert_eq!(out.len(), wanted);
    out
}

fn group_support(base: &[Point], p: u64) -> BTreeSet<Point> {
    let mut out = BTreeSet::new();
    for i in 0..base.len() {
        for j in i..base.len() {
            out.insert(add(base[i], base[j], p));
        }
    }
    out
}

fn scalar_support(base: &[Point], point_logs: &BTreeMap<Point, u64>, r: u64) -> BTreeSet<u64> {
    let logs: Vec<u64> = base.iter().map(|point| point_logs[point]).collect();
    let mut out = BTreeSet::new();
    for i in 0..logs.len() {
        for j in i..logs.len() {
            out.insert((logs[i] + logs[j]) % r);
        }
    }
    out
}

fn run(p: u64, r: u64) {
    let t = Some((0, 1));
    let t2 = add(t, t, p);
    assert_eq!(t2, Some((0, p - 1)));
    assert_eq!(add(t2, t, p), None);
    let torsion = [None, t, t2];
    let base36 = base(p, 36);
    let base12 = &base36[..12];
    let g = base12[0];
    let mut point_logs = BTreeMap::new();
    let mut current = None;
    for k in 0..r {
        assert!(point_logs.insert(current, k).is_none());
        current = add(current, g, p);
    }
    assert_eq!(current, None);
    assert_eq!(point_logs.len(), r as usize);
    let destination_b = sub_mod(0, 27, p);
    let mut orbit = Vec::new();
    for &point in base12 {
        assert!(on_curve(point, p, 1));
        assert_eq!(scale(point, r, p), None);
        let mut images = BTreeSet::new();
        for &shift in &torsion {
            let lifted = add(point, shift, p);
            assert!(on_curve(lifted, p, 1));
            assert!(images.insert(lifted));
            assert_eq!(phi(lifted, p), phi(point, p));
            assert!(on_curve(phi(lifted, p), p, destination_b));
            assert_eq!(scale(lifted, 6, p), scale(point, 6, p));
            orbit.push((lifted, shift));
        }
        assert_eq!(images.len(), 3);
    }
    for &point in &base36 {
        assert!(on_curve(point, p, 1));
        assert_eq!(scale(point, r, p), None);
    }
    assert_eq!(orbit.len(), 36);
    let support12 = group_support(base12, p);
    let support36 = group_support(&base36, p);
    let mut orbit_support = BTreeSet::new();
    for i in 0..orbit.len() {
        for j in i..orbit.len() {
            let sum = add(orbit[i].0, orbit[j].0, p);
            let correction = neg(add(orbit[i].1, orbit[j].1, p), p);
            orbit_support.insert(add(sum, correction, p));
        }
    }
    assert_eq!(orbit_support, support12);
    for (points, support) in [(base12, &support12), (&base36, &support36)] {
        let scalars = scalar_support(points, &point_logs, r);
        let via_group: BTreeSet<u64> = support.iter().map(|point| point_logs[point]).collect();
        assert_eq!(scalars, via_group);
    }
    let orbit_scalars: BTreeSet<u64> = orbit_support.iter().map(|point| point_logs[point]).collect();
    assert_eq!(orbit_scalars, scalar_support(base12, &point_logs, r));
    for (arm, reps, unknowns, pairs, support) in [
        ("subgroup12", 12, 12, 78, &support12),
        ("torsion_orbit36", 36, 12, 666, &orbit_support),
        ("distinct_subgroup36", 36, 36, 666, &support36),
    ] {
        println!("{p}\t{r}\t{arm}\t{reps}\t{unknowns}\t{pairs}\t{}\t{r}\t{:.9}",
            support.len(), support.len() as f64 / r as f64);
    }
    eprintln!("verified p={p} r={r} orbit_points=36 quotient_classes=12 scalar_index={r} orbit_equals_base=true");
}

fn main() {
    println!("p\tr\tarm\trepresentatives\tunknown_logs\traw_pairs\tdistinct_two_sum_support\tgroup_size\texact_hit_probability");
    for (p, r) in [(617, 103), (2633, 439), (9677, 1613), (38561, 6427)] {
        run(p, r);
    }
}
