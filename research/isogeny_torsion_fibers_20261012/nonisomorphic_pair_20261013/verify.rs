//! Independent finite-field, map, subgroup, target, and row checker.
//! Deliberately does not import the runner's curve or relation code.
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;

type Point = Option<(i64, i64)>;
const P: i64 = 617;
const R: u64 = 103;
const RELATION_SEED: u64 = 0x5245_4c41_5449_4f4e;
const DESCENT_SEED: u64 = 0x4445_5343_454e_5400;

fn modp(x: i64) -> i64 {
    x.rem_euclid(P)
}
fn inverse(a: i64) -> i64 {
    let (mut u, mut v, mut s, mut t) = (modp(a), P, 1, 0);
    while v != 0 {
        let q = u / v;
        (u, v) = (v, u - q * v);
        (s, t) = (t, s - q * t);
    }
    assert_eq!(u, 1);
    modp(s)
}
fn on(point: Point, b: i64) -> bool {
    point.is_none_or(|(x, y)| modp(y * y - x * x * x - b) == 0)
}
fn add(left: Point, right: Point, b: i64) -> Point {
    match (left, right) {
        (None, p) | (p, None) => p,
        (Some((x, y)), Some((u, v))) => {
            if x == u && modp(y + v) == 0 {
                return None;
            }
            let slope = if x == u {
                modp(3 * x * x * inverse(2 * y))
            } else {
                modp((v - y) * inverse(u - x))
            };
            let nx = modp(slope * slope - x - u);
            let ny = modp(slope * (x - nx) - y);
            let result = Some((nx, ny));
            assert!(on(result, b));
            result
        }
    }
}
fn times(mut n: u64, mut point: Point, b: i64) -> Point {
    let mut result = None;
    while n != 0 {
        if n & 1 != 0 {
            result = add(result, point, b);
        }
        point = add(point, point, b);
        n >>= 1;
    }
    result
}
fn map(point: Point) -> Point {
    match point {
        None | Some((0, _)) => None,
        Some((x, y)) => {
            let x3 = modp(x * x * x);
            Some((
                modp((x3 + 4) * inverse(x * x)),
                modp(y * (x3 - 8) * inverse(x3)),
            ))
        }
    }
}
fn points(b: i64) -> Vec<Point> {
    let mut out = vec![None];
    for x in 0..P {
        for y in 0..P {
            let point = Some((x, y));
            if on(point, b) {
                out.push(point);
            }
        }
    }
    out
}
fn prime(n: u64) -> bool {
    n >= 2 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0)
}
struct SplitMix64(u64);
impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}
fn targets() -> Vec<Point> {
    let mut rng = SplitMix64(0x2026_1013_0001_0000 + P as u64);
    let mut out = Vec::new();
    let mut seen = BTreeSet::new();
    while out.len() < 8 {
        let x = (rng.next() % P as u64) as i64;
        let sign = rng.next() & 1;
        let root = (0..P).find(|&y| on(Some((x, y)), 1));
        if let Some(y) = root {
            let point = Some((x, if sign == 1 && y != 0 { P - y } else { y }));
            let q = times(6, point, 1);
            if q.is_some() && seen.insert(q) {
                out.push(q);
            }
        }
    }
    out
}
fn factor_base(b: i64) -> (Vec<Point>, usize) {
    let mut seen = BTreeSet::new();
    let mut base = Vec::new();
    for x in 0..P {
        if let Some(y) = (0..P).find(|&y| on(Some((x, y)), b)) {
            let projected = times(6, Some((x, y)), b);
            if projected.is_some() && seen.insert(projected) {
                base.push(projected);
            }
            if base.len() == 12 {
                return (base, x as usize + 1);
            }
        }
    }
    panic!("factor base unavailable")
}
fn factor_base_lifts(b: i64) -> usize {
    factor_base(b).1
}
fn modr(value: i64) -> u64 {
    value.rem_euclid(R as i64) as u64
}
fn relation_replay(b: i64, generator: Point, target: Point, index: usize) -> (u64, u64, u64, u64) {
    let (base, _) = factor_base(b);
    let mut pairs = BTreeMap::<Point, Vec<(usize, usize)>>::new();
    for i in 0..base.len() {
        for j in i..base.len() {
            pairs
                .entry(add(base[i], base[j], b))
                .or_default()
                .push((i, j));
        }
    }
    let key = P as u64 ^ ((index as u64) << 16);
    let mut stream = SplitMix64(RELATION_SEED ^ key);
    let mut pivots = vec![None::<(Vec<u64>, u64)>; 12];
    let mut rank = 0;
    let mut mul_subs = 0;
    let mut collection_trials = 0;
    for trial in 1..=10_000 {
        collection_trials = trial;
        let k = 1 + stream.next() % (R - 1);
        let point = times(k, generator, b);
        if let Some(representations) = pairs.get(&point) {
            for &(i, j) in representations {
                let mut row = vec![0u64; 12];
                row[i] += 1;
                row[j] += 1;
                let mut rhs = k;
                for col in 0..12 {
                    let c = row[col];
                    if c == 0 {
                        continue;
                    }
                    if let Some((old, old_rhs)) = &pivots[col] {
                        for j in col..12 {
                            row[j] = modr(row[j] as i64 - (c * old[j]) as i64);
                            mul_subs += 1;
                        }
                        rhs = modr(rhs as i64 - (c * old_rhs) as i64);
                        mul_subs += 1;
                    } else {
                        let inv = (1..R).find(|&v| c * v % R == 1).unwrap();
                        for value in row.iter_mut().skip(col) {
                            *value = *value * inv % R;
                            mul_subs += 1;
                        }
                        rhs = rhs * inv % R;
                        mul_subs += 1;
                        pivots[col] = Some((row, rhs));
                        rank += 1;
                        break;
                    }
                }
                if rank == 12 {
                    break;
                }
            }
        }
        if rank == 12 {
            break;
        }
    }
    assert_eq!(rank, 12);
    let mut logs = [0u64; 12];
    for i in (0..12).rev() {
        let (row, mut rhs) = pivots[i].clone().unwrap();
        for j in i + 1..12 {
            rhs = modr(rhs as i64 - (row[j] * logs[j]) as i64);
            mul_subs += 1;
        }
        logs[i] = rhs;
        assert_eq!(times(rhs, generator, b), base[i]);
    }
    let mut stream = SplitMix64(DESCENT_SEED ^ key);
    for trial in 1..=10_000 {
        let t = stream.next() % R;
        let probe = add(target, times(t, generator, b), b);
        if let Some(representations) = pairs.get(&probe) {
            let (i, j) = representations[0];
            let d = modr(logs[i] as i64 + logs[j] as i64 - t as i64);
            assert_eq!(times(d, generator, b), target);
            return (d, collection_trials, trial, mul_subs);
        }
    }
    panic!("independent descent cap")
}
fn val<'a>(row: &'a [&str], cols: &BTreeMap<&str, usize>, key: &str) -> &'a str {
    row[*cols.get(key).expect("column")]
}
fn number(row: &[&str], cols: &BTreeMap<&str, usize>, key: &str) -> u64 {
    val(row, cols, key).parse().expect("integer")
}

fn main() {
    let path = env::args()
        .nth(1)
        .expect("usage: verify raw-v2.tsv fibers-v2.tsv");
    let fiber_path = env::args().nth(2).expect("fiber output path");
    assert!(prime(P as u64) && prime(R));
    let first = (101..)
        .filter(|&r| prime(r) && prime(6 * r - 1))
        .next()
        .unwrap();
    assert_eq!(first, R);
    let source = points(1);
    let dest = points(590);
    assert_eq!(source.len(), 618);
    assert_eq!(dest.len(), 618);
    assert!(!(1..P).any(|u| modp(u * u * u * u * u * u) == 590));
    let mut fibers = BTreeMap::<Point, usize>::new();
    for &point in &source {
        let image = map(point);
        assert!(on(image, 590));
        *fibers.entry(image).or_default() += 1;
    }
    assert_eq!(fibers.len(), 206);
    assert!(fibers.values().all(|&n| n == 3));
    assert_eq!(
        source
            .iter()
            .copied()
            .filter(|&point| map(point).is_none())
            .collect::<Vec<_>>(),
        vec![None, Some((0, 1)), Some((0, 616))]
    );
    for &left in &source {
        for &right in &source {
            assert_eq!(map(add(left, right, 1)), add(map(left), map(right), 590));
        }
    }
    let generator = source
        .iter()
        .skip(1)
        .copied()
        .map(|p| times(6, p, 1))
        .find(|p| p.is_some())
        .unwrap();
    assert_eq!(generator, Some((66, 264)));
    assert_eq!(times(R, generator, 1), None);
    assert_eq!(times(R, map(generator), 590), None);
    let frozen = targets();
    let kernel = [None, Some((0, 1)), Some((0, 616))];
    let mut fiber_tsv = String::from("target_index\tqx\tqy\timage_x\timage_y\tpreimage0_x\tpreimage0_y\tpreimage1_x\tpreimage1_y\tpreimage2_x\tpreimage2_y\n");
    for (i, &q) in frozen.iter().enumerate() {
        assert_eq!(times(R, q, 1), None);
        let fiber: BTreeSet<Point> = source
            .iter()
            .copied()
            .filter(|&point| map(point) == map(q))
            .collect();
        let expected: BTreeSet<Point> = kernel.iter().copied().map(|k| add(q, k, 1)).collect();
        assert_eq!(fiber.len(), 3);
        assert_eq!(fiber, expected);
        let (qx, qy) = q.unwrap();
        let (ix, iy) = map(q).unwrap();
        fiber_tsv.push_str(&format!("{i}\t{qx}\t{qy}\t{ix}\t{iy}"));
        for point in fiber {
            let (x, y) = point.unwrap();
            fiber_tsv.push_str(&format!("\t{x}\t{y}"));
        }
        fiber_tsv.push('\n');
    }
    let input = fs::read_to_string(path).expect("raw TSV");
    let mut lines = input.lines();
    let header: Vec<_> = lines.next().expect("header").split('\t').collect();
    let cols: BTreeMap<&str, usize> = header
        .iter()
        .enumerate()
        .map(|(i, &name)| (name, i))
        .collect();
    let mut arms = BTreeSet::new();
    let mut values = BTreeMap::<usize, u64>::new();
    let mut count = 0;
    for line in lines {
        let row: Vec<_> = line.split('\t').collect();
        assert_eq!(row.len(), header.len());
        assert_eq!(number(&row, &cols, "p"), P as u64);
        assert_eq!(number(&row, &cols, "b"), 1);
        assert_eq!(number(&row, &cols, "r"), R);
        let i = number(&row, &cols, "target_index") as usize;
        assert!(i < 8);
        let arm = val(&row, &cols, "arm");
        assert!(["source_relations", "transport_relations", "rho_negation"].contains(&arm));
        assert!(arms.insert((i, arm.to_string())));
        let q = Some((
            number(&row, &cols, "qx") as i64,
            number(&row, &cols, "qy") as i64,
        ));
        let mapped = Some((
            number(&row, &cols, "mapped_qx") as i64,
            number(&row, &cols, "mapped_qy") as i64,
        ));
        assert_eq!(q, frozen[i]);
        assert_eq!(mapped, map(q));
        assert_eq!(val(&row, &cols, "status"), "verified");
        let d = number(&row, &cols, "log");
        assert_eq!(times(d, generator, 1), q);
        if arm == "transport_relations" {
            assert_eq!(times(d, map(generator), 590), mapped);
        }
        if let Some(old) = values.insert(i, d) {
            assert_eq!(old, d);
        }
        let phases = [
            "setup_gae",
            "collection_gae",
            "descent_gae",
            "verification_gae",
            "rho_gae",
        ];
        let total: u64 = phases.iter().map(|key| number(&row, &cols, key)).sum();
        assert_eq!(total, number(&row, &cols, "total_gae"));
        for phase in ["setup", "collection", "descent", "verification", "rho"] {
            let gae = format!("{phase}_gae");
            let adds = format!("{phase}_adds");
            let doubles = format!("{phase}_doubles");
            assert_eq!(
                number(&row, &cols, &gae),
                number(&row, &cols, &adds) + number(&row, &cols, &doubles)
            );
        }
        if arm.ends_with("relations") {
            let transported = arm == "transport_relations";
            let (recovered, relation_trials, descent_trials, row_mul_subs) = relation_replay(
                if transported { 590 } else { 1 },
                if transported {
                    map(generator)
                } else {
                    generator
                },
                if transported { mapped } else { q },
                i,
            );
            assert_eq!(recovered, d);
            assert_eq!(relation_trials, number(&row, &cols, "relation_trials"));
            assert_eq!(descent_trials, number(&row, &cols, "descent_trials"));
            assert_eq!(row_mul_subs, number(&row, &cols, "row_mul_subs"));
            assert_eq!(
                number(&row, &cols, "fb_lifts") as usize,
                factor_base_lifts(if transported { 590 } else { 1 })
            );
            assert_eq!(number(&row, &cols, "pair_entries"), 78);
            assert_eq!(number(&row, &cols, "independent_rows"), 12);
            assert_eq!(number(&row, &cols, "row_inversions"), 12);
            assert_eq!(
                number(&row, &cols, "relation_lookups"),
                number(&row, &cols, "relation_trials")
            );
            assert_eq!(
                number(&row, &cols, "descent_lookups"),
                number(&row, &cols, "descent_trials")
            );
            let evals = number(&row, &cols, "map_evaluations");
            assert_eq!(evals, if transported { 2 } else { 0 });
            assert_eq!(
                number(&row, &cols, "map_muls"),
                if transported { 10 } else { 0 }
            );
            assert_eq!(
                number(&row, &cols, "map_sqrs"),
                if transported { 4 } else { 0 }
            );
            assert_eq!(
                number(&row, &cols, "map_invs"),
                if transported { 2 } else { 0 }
            );
        } else {
            assert!(number(&row, &cols, "rho_steps") > 0);
            assert_eq!(number(&row, &cols, "map_evaluations"), 0);
        }
        count += 1;
    }
    assert_eq!(count, 24);
    assert_eq!(arms.len(), 24);
    assert_eq!(values.len(), 8);
    fs::write(fiber_path, fiber_tsv).expect("fiber TSV");
    println!("verified_records={count} replayed_relation_rows=16 source_targets={} checked_target_fibers=8 homomorphism_pairs={} fiber_size=3 image_size=206 source_fb_lifts={} destination_fb_lifts={}",
        values.len(), source.len()*source.len(), factor_base_lifts(1), factor_base_lifts(590));
}
