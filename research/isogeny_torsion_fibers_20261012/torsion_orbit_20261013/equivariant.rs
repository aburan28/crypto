//! Exhaustive rational-point check of EQUIVARIANT_PROTOCOL.md.
use std::collections::BTreeSet;

type Point = Option<(u64, u64)>;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum Projective {
    Finite(u64),
    Infinity,
}

fn mul(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn pow(mut a: u64, mut n: u64, p: u64) -> u64 {
    let mut out = 1;
    while n > 0 {
        if n & 1 == 1 {
            out = mul(out, a, p);
        }
        a = mul(a, a, p);
        n >>= 1;
    }
    out
}

fn sub(a: u64, b: u64, p: u64) -> u64 {
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
        mul(mul(3, mul(x1, x1, p), p), pow(2 * y1, p - 2, p), p)
    } else {
        mul(sub(y2, y1, p), pow(sub(x2, x1, p), p - 2, p), p)
    };
    let x3 = sub(sub(mul(slope, slope, p), x1, p), x2, p);
    let y3 = sub(mul(slope, sub(x1, x3, p), p), y1, p);
    Some((x3, y3))
}

fn on_curve(point: Point, p: u64, b: u64) -> bool {
    let Some((x, y)) = point else { return true; };
    mul(y, y, p) == (mul(mul(x, x, p), x, p) + b) % p
}

fn all_points(p: u64) -> Vec<Point> {
    let mut roots = vec![Vec::<u64>::new(); p as usize];
    for y in 0..p {
        roots[mul(y, y, p) as usize].push(y);
    }
    let mut out = vec![None];
    for x in 0..p {
        let rhs = (mul(mul(x, x, p), x, p) + 1) % p;
        for &y in &roots[rhs as usize] {
            out.push(Some((x, y)));
        }
    }
    out
}

fn chi(point: Point, p: u64) -> Projective {
    let Some((_, y)) = point else { return Projective::Finite(1); };
    if y == 1 {
        Projective::Infinity
    } else {
        Projective::Finite(mul((y + 1) % p, pow(sub(y, 1, p), p - 2, p), p))
    }
}

fn action(z: Projective, p: u64) -> Projective {
    match z {
        Projective::Infinity => Projective::Finite(0),
        Projective::Finite(1) => Projective::Infinity,
        Projective::Finite(v) => Projective::Finite(pow(sub(1, v, p), p - 2, p)),
    }
}

fn phi(point: Point, p: u64) -> Point {
    let Some((x, y)) = point else { return None; };
    if x == 0 {
        return None;
    }
    let x2 = mul(x, x, p);
    let x3 = mul(x2, x, p);
    let ix = pow(x, p - 2, p);
    let ix2 = mul(ix, ix, p);
    let ix3 = mul(ix2, ix, p);
    Some((
        mul((x3 + 4) % p, ix2, p),
        mul(mul(y, sub(x3, 8, p), p), ix3, p),
    ))
}

fn run(p: u64, r: u64) {
    let t = Some((0, 1));
    let t2 = add(t, t, p);
    assert_eq!(t2, Some((0, p - 1)));
    assert_eq!(add(t2, t, p), None);
    let exceptional = [chi(None, p), chi(t, p), chi(t2, p)];
    assert_eq!(exceptional, [Projective::Finite(1), Projective::Infinity, Projective::Finite(0)]);
    assert_eq!(exceptional.iter().copied().collect::<BTreeSet<_>>().len(), 3);
    assert_eq!(phi(None, p), phi(t, p));
    assert_eq!(phi(t, p), phi(t2, p));
    let mut projective_checked = 0;
    for z in (0..p).map(Projective::Finite).chain([Projective::Infinity]) {
        assert_eq!(action(action(action(z, p), p), p), z);
        projective_checked += 1;
    }
    let points = all_points(p);
    assert_eq!(points.len() as u64, p + 1);
    let destination_b = sub(0, 27, p);
    for &point in &points {
        assert!(on_curve(point, p, 1));
        let plus_t = add(point, t, p);
        let plus_t2 = add(plus_t, t, p);
        assert!(on_curve(plus_t, p, 1));
        assert!(on_curve(plus_t2, p, 1));
        assert_eq!(chi(plus_t, p), action(chi(point, p), p));
        assert_eq!(chi(plus_t2, p), action(action(chi(point, p), p), p));
        let image = phi(point, p);
        assert!(on_curve(image, p, destination_b));
        assert_eq!(phi(plus_t, p), image);
        assert_eq!(phi(plus_t2, p), image);
    }
    println!("{p}\t{r}\t{}\t{projective_checked}\t{}\t3\tpass",
        points.len(), 2 * points.len());
    eprintln!("p={p} r={r} chi(O)=1 chi(T)=infinity chi(2T)=0 nonfactorable=true");
}

fn main() {
    println!("p\tr\trational_points\tprojective_action_points\tcoordinate_identities\tdistinct_kernel_chi\tstatus");
    for (p, r) in [(617, 103), (2633, 439), (9677, 1613), (38561, 6427)] {
        run(p, r);
    }
}
