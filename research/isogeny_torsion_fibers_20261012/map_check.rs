// Exact public map check from MAP_PROTOCOL.md. No DLP or timing.
use std::collections::BTreeSet;

type Point = Option<(i64, i64)>;
const P: i64 = 13;

fn m(x: i64) -> i64 {
    x.rem_euclid(P)
}

fn pow(mut a: i64, mut n: u64) -> i64 {
    let mut result = 1;
    while n > 0 {
        if n & 1 == 1 {
            result = m(result * a);
        }
        a = m(a * a);
        n >>= 1;
    }
    result
}

fn inv(x: i64) -> i64 {
    assert_ne!(m(x), 0);
    pow(m(x), (P - 2) as u64)
}

fn points(b: i64) -> Vec<Point> {
    let mut result = vec![None];
    for x in 0..P {
        for y in 0..P {
            if m(y * y) == m(x * x * x + b) {
                result.push(Some((x, y)));
            }
        }
    }
    result
}

fn add(a: Point, b: Point, curve_b: i64) -> Point {
    match (a, b) {
        (None, p) | (p, None) => p,
        (Some((x1, y1)), Some((x2, y2))) => {
            assert_eq!(m(y1 * y1), m(x1 * x1 * x1 + curve_b));
            assert_eq!(m(y2 * y2), m(x2 * x2 * x2 + curve_b));
            if x1 == x2 && m(y1 + y2) == 0 {
                return None;
            }
            let slope = if x1 == x2 {
                m(3 * x1 * x1 * inv(2 * y1))
            } else {
                m((y2 - y1) * inv(x2 - x1))
            };
            let x3 = m(slope * slope - x1 - x2);
            let y3 = m(slope * (x1 - x3) - y1);
            Some((x3, y3))
        }
    }
}

fn phi(point: Point, source_b: i64) -> Point {
    match point {
        None => None,
        Some((x, y)) => {
            if x == 0 {
                return None;
            }
            let x2 = m(x * x);
            let x3 = m(x2 * x);
            Some((
                m((x3 + 4 * source_b) * inv(x2)),
                m(y * (x3 - 8 * source_b) * inv(x3)),
            ))
        }
    }
}

fn main() {
    let source_b = 2;
    let target_b = m(-27 * source_b);
    assert_eq!(target_b, 11);
    let source = points(source_b);
    let target = points(target_b);
    assert_eq!(source.len(), 19);
    assert_eq!(target.len(), 19);
    let target_set: BTreeSet<Point> = target.iter().copied().collect();
    let images: BTreeSet<Point> = source
        .iter()
        .map(|&point| {
            let image = phi(point, source_b);
            assert!(target_set.contains(&image));
            image
        })
        .collect();
    assert_eq!(images.len(), source.len());
    let mut pairs = 0u64;
    for &a in &source {
        for &b in &source {
            assert_eq!(
                phi(add(a, b, source_b), source_b),
                add(phi(a, source_b), phi(b, source_b), target_b)
            );
            pairs += 1;
        }
    }
    println!("p=13 source_b=2 target_b=11 source_order={} target_order={} rational_kernel_size=1 image_count={} homomorphism_pairs={} status=pass", source.len(), target.len(), images.len(), pairs);
}
