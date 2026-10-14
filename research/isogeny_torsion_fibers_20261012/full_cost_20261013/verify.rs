//! Independent arithmetic and row checker. Build: rustc --edition=2021 -O verify.rs -o verify.
//! This deliberately does not import the experiment's curve implementation.
use std::collections::{BTreeMap, BTreeSet};
use std::{env, fs};

type Point = Option<(i128, i128)>;

fn residue(n: i128, p: i128) -> i128 {
    n.rem_euclid(p)
}
fn inverse(a: i128, p: i128) -> i128 {
    let (mut old, mut now, mut x, mut y) = (residue(a, p), p, 1, 0);
    while now != 0 {
        let q = old / now;
        (old, now) = (now, old - q * now);
        (x, y) = (y, x - q * y);
    }
    assert_eq!(old, 1);
    residue(x, p)
}
fn on(q: Point, b: i128, p: i128) -> bool {
    q.is_none_or(|(x, y)| residue(y * y - x * x * x - b, p) == 0)
}
fn add(a: Point, b: Point, curve_b: i128, p: i128) -> Point {
    match (a, b) {
        (None, q) | (q, None) => q,
        (Some((x, y)), Some((u, v))) => {
            if x == u && residue(y + v, p) == 0 {
                return None;
            }
            let slope = if x == u {
                residue(3 * x * x * inverse(2 * y, p), p)
            } else {
                residue((v - y) * inverse(u - x, p), p)
            };
            let q = Some((residue(slope * slope - x - u, p), 0));
            let nx = q.unwrap().0;
            let q = Some((nx, residue(slope * (x - nx) - y, p)));
            assert!(on(q, curve_b, p));
            q
        }
    }
}
fn times(mut n: u64, mut q: Point, b: i128, p: i128) -> Point {
    let mut sum = None;
    while n > 0 {
        if n & 1 == 1 {
            sum = add(sum, q, b, p);
        }
        q = add(q, q, b, p);
        n >>= 1;
    }
    sum
}
fn map(q: Point, p: i128) -> Point {
    match q {
        None | Some((0, _)) => None,
        Some((x, y)) => {
            let x3 = residue(x * x * x, p);
            Some((
                residue((x3 + 4) * inverse(x * x, p), p),
                residue(y * (x3 - 8) * inverse(x3, p), p),
            ))
        }
    }
}
fn square_roots(p: i128) -> Vec<Option<i128>> {
    let mut roots = vec![None; p as usize];
    for y in (0..p).rev() {
        roots[residue(y * y, p) as usize] = Some(y);
    }
    roots
}
fn rhs(x: i128, b: i128, p: i128) -> i128 {
    residue(x * x * x + b, p)
}
fn generator(p: i128, roots: &[Option<i128>]) -> Point {
    for x in 0..p {
        if let Some(y) = roots[rhs(x, 1, p) as usize] {
            let projected = times(6, Some((x, y)), 1, p);
            if projected.is_some() {
                return projected;
            }
        }
    }
    panic!("no generator")
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
fn targets(p: i128, roots: &[Option<i128>]) -> Vec<Point> {
    let mut stream = SplitMix64(0x2026_1013_0001_0000 + p as u64);
    let mut out = Vec::new();
    let mut seen = BTreeSet::new();
    while out.len() < 8 {
        let x = (stream.next() % p as u64) as i128;
        let sign = stream.next() & 1;
        if let Some(y0) = roots[rhs(x, 1, p) as usize] {
            let y = if sign == 1 && y0 != 0 { p - y0 } else { y0 };
            let q = times(6, Some((x, y)), 1, p);
            if q.is_some() && seen.insert(q) {
                out.push(q);
            }
        }
    }
    out
}
fn prime(n: u64) -> bool {
    n >= 2 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0)
}
fn selected() -> Vec<(u64, u64)> {
    [(101, 400), (401, 1600), (1601, 6400), (6401, 25600)]
        .into_iter()
        .map(|(lo, hi)| {
            let r = (lo..=hi).find(|&r| prime(r) && prime(6 * r - 1)).unwrap();
            (6 * r - 1, r)
        })
        .collect()
}
fn point(row: &[&str], col: &BTreeMap<&str, usize>, x: &str, y: &str) -> Point {
    Some((row[col[x]].parse().unwrap(), row[col[y]].parse().unwrap()))
}
fn num(row: &[&str], col: &BTreeMap<&str, usize>, name: &str) -> u64 {
    row[col[name]].parse().unwrap()
}

const RELATION_SEED: u64 = 0x5245_4c41_5449_4f4e;
const DESCENT_SEED: u64 = 0x4445_5343_454e_5400;

#[derive(Default)]
struct RelationReplay {
    lifts: u64,
    collection_trials: u64,
    rank: u64,
    row_mul_subs: u64,
    row_inversions: u64,
    descent_trials: u64,
    recovered: Option<u64>,
}

fn modular(n: i128, r: u64) -> u64 {
    n.rem_euclid(r as i128) as u64
}

fn independent_factor_base(b: i128, p: i128, roots: &[Option<i128>]) -> (Vec<Point>, u64) {
    let mut base = Vec::new();
    let mut seen = BTreeSet::new();
    for x in 0..p {
        if let Some(y) = roots[rhs(x, b, p) as usize] {
            let projected = times(6, Some((x, y)), b, p);
            if projected.is_some() && seen.insert(projected) {
                base.push(projected);
                if base.len() == 12 {
                    return (base, (x + 1) as u64);
                }
            }
        }
    }
    panic!("insufficient independent factor base")
}

fn insert_relation(
    pivots: &mut [Option<(Vec<u64>, u64)>],
    mut coeff: Vec<u64>,
    mut rhs: u64,
    r: u64,
    tally: &mut RelationReplay,
) {
    for col in 0..coeff.len() {
        let c = coeff[col];
        if c == 0 {
            continue;
        }
        if let Some((old, old_rhs)) = &pivots[col] {
            for j in col..coeff.len() {
                coeff[j] = modular(coeff[j] as i128 - c as i128 * old[j] as i128, r);
                tally.row_mul_subs += 1;
            }
            rhs = modular(rhs as i128 - c as i128 * *old_rhs as i128, r);
            tally.row_mul_subs += 1;
        } else {
            let scale = inverse(c as i128, r as i128) as u64;
            tally.row_inversions += 1;
            for item in coeff.iter_mut().skip(col) {
                *item = modular(*item as i128 * scale as i128, r);
                tally.row_mul_subs += 1;
            }
            rhs = modular(rhs as i128 * scale as i128, r);
            tally.row_mul_subs += 1;
            pivots[col] = Some((coeff, rhs));
            tally.rank += 1;
            return;
        }
    }
    assert_eq!(rhs, 0, "inconsistent relation");
}

fn independent_relation_replay(
    p: i128,
    b: i128,
    r: u64,
    generator: Point,
    target: Point,
    index: usize,
    round: u64,
    roots: &[Option<i128>],
) -> RelationReplay {
    let (base, lifts) = independent_factor_base(b, p, roots);
    let mut out = RelationReplay {
        lifts,
        ..RelationReplay::default()
    };
    let mut pairs = BTreeMap::<Point, Vec<(usize, usize)>>::new();
    for i in 0..12 {
        for j in i..12 {
            pairs
                .entry(add(base[i], base[j], b, p))
                .or_default()
                .push((i, j));
        }
    }
    let key = p as u64 ^ ((index as u64) << 16) ^ (round << 32);
    let mut stream = SplitMix64(RELATION_SEED ^ key);
    let mut pivots = vec![None::<(Vec<u64>, u64)>; 12];
    for trial in 1..=100_000 {
        out.collection_trials = trial;
        let k = 1 + stream.next() % (r - 1);
        let point = times(k, generator, b, p);
        if let Some(representations) = pairs.get(&point) {
            for &(i, j) in representations {
                let mut coeff = vec![0u64; 12];
                coeff[i] += 1;
                coeff[j] += 1;
                insert_relation(&mut pivots, coeff, k, r, &mut out);
                if out.rank == 12 {
                    break;
                }
            }
        }
        if out.rank == 12 {
            break;
        }
    }
    assert_eq!(
        out.rank, 12,
        "independent replay failed to collect full rank"
    );
    let mut logs = vec![0u64; 12];
    for i in (0..12).rev() {
        let (coeff, rhs) = pivots[i].as_ref().unwrap();
        let mut value = *rhs;
        for (j, &known) in logs.iter().enumerate().skip(i + 1) {
            value = modular(value as i128 - coeff[j] as i128 * known as i128, r);
            out.row_mul_subs += 1;
        }
        logs[i] = value;
    }
    for (&point, &log) in base.iter().zip(&logs) {
        assert_eq!(
            times(log, generator, b, p),
            point,
            "factor-base log mismatch"
        );
    }
    let mut stream = SplitMix64(DESCENT_SEED ^ key);
    for trial in 1..=100_000 {
        out.descent_trials = trial;
        let t = stream.next() % r;
        let probe = add(target, times(t, generator, b, p), b, p);
        if let Some(representations) = pairs.get(&probe) {
            let (i, j) = representations[0];
            let d = (logs[i] + logs[j] + r - t) % r;
            assert_eq!(times(d, generator, b, p), target, "descent mismatch");
            out.recovered = Some(d);
            return out;
        }
    }
    panic!("independent descent cap")
}

fn main() {
    let path = env::args().nth(1).expect("usage: verify raw.tsv");
    let input = fs::read_to_string(path).unwrap();
    let mut lines = input.lines();
    let header = lines.next().unwrap();
    let names: Vec<_> = header.split('\t').collect();
    let col: BTreeMap<&str, usize> = names.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let mut per = BTreeMap::<(u64, u64, u64), Vec<(String, u64, u64, u64)>>::new();
    let mut count = 0usize;
    let mut measured = 0usize;
    for line in lines {
        let row: Vec<_> = line.split('\t').collect();
        assert_eq!(row.len(), names.len());
        let p = num(&row, &col, "p");
        let r = num(&row, &col, "r");
        assert!(selected().contains(&(p, r)));
        let p128 = p as i128;
        let b = (p128 - 27).rem_euclid(p128);
        let roots = square_roots(p128);
        let g = generator(p128, &roots);
        let index = num(&row, &col, "target_index") as usize;
        let round = num(&row, &col, "round");
        assert!(index < 8 && round <= 5);
        let q = point(&row, &col, "qx", "qy");
        let q2 = point(&row, &col, "mapped_qx", "mapped_qy");
        assert_eq!(q, targets(p128, &roots)[index]);
        assert_eq!(q2, map(q, p128));
        assert!(on(q, 1, p128) && on(q2, b, p128));
        assert_eq!(times(r, g, 1, p128), None);
        let g2 = map(g, p128);
        assert_eq!(times(r, g2, b, p128), None);
        assert_eq!(row[col["status"]], "verified");
        let d = num(&row, &col, "log");
        assert!(d < r);
        assert_eq!(times(d, g, 1, p128), q);
        let arm = row[col["arm"]];
        if arm == "source_relations" || arm == "transport_relations" {
            let (curve_b, generator, target) = if arm == "source_relations" {
                (1, g, q)
            } else {
                (b, g2, q2)
            };
            let rel = independent_relation_replay(
                p128, curve_b, r, generator, target, index, round, &roots,
            );
            assert_eq!(rel.recovered, Some(d));
            assert_eq!(rel.lifts, num(&row, &col, "fb_lifts"));
            assert_eq!(num(&row, &col, "pair_entries"), 78);
            assert_eq!(rel.collection_trials, num(&row, &col, "relation_trials"));
            assert_eq!(rel.collection_trials, num(&row, &col, "relation_lookups"));
            assert_eq!(rel.rank, num(&row, &col, "independent_rows"));
            assert_eq!(rel.row_mul_subs, num(&row, &col, "row_mul_subs"));
            assert_eq!(rel.row_inversions, num(&row, &col, "row_inversions"));
            assert_eq!(rel.descent_trials, num(&row, &col, "descent_trials"));
            assert_eq!(rel.descent_trials, num(&row, &col, "descent_lookups"));
        }
        if arm == "transport_relations" {
            assert_eq!(times(d, g2, b, p128), q2);
            assert_eq!(num(&row, &col, "map_evaluations"), 2);
        } else {
            assert_eq!(num(&row, &col, "map_evaluations"), 0);
        }
        let wall = num(&row, &col, "wall_ns");
        assert!(wall > 0);
        if col.contains_key("cpu_ns") {
            let cpu = num(&row, &col, "cpu_ns");
            assert!(cpu > 0 && cpu <= wall);
        }
        let parts = [
            "map_ns",
            "setup_ns",
            "collection_ns",
            "recovery_ns",
            "descent_ns",
            "verification_ns",
            "rho_ns",
            "other_ns",
        ];
        let phase_sum: u64 = parts.iter().map(|&n| num(&row, &col, n)).sum();
        assert_eq!(wall, phase_sum, "phase clocks must partition whole arm");
        assert_eq!(
            num(&row, &col, "total_gae"),
            [
                "setup_gae",
                "collection_gae",
                "descent_gae",
                "verification_gae",
                "rho_gae"
            ]
            .iter()
            .map(|&n| num(&row, &col, n))
            .sum::<u64>()
        );
        assert_eq!(num(&row, &col, "warmup"), u64::from(round == 0));
        if round > 0 {
            measured += 1;
        }
        per.entry((p, index as u64, round)).or_default().push((
            arm.into(),
            d,
            num(&row, &col, "total_gae"),
            num(&row, &col, "rho_steps"),
        ));
        count += 1;
    }
    assert_eq!(count, 768);
    assert_eq!(measured, 640);
    assert_eq!(per.len(), 192);
    for arms in per.values() {
        assert_eq!(arms.len(), 4);
        let names: BTreeSet<_> = arms.iter().map(|a| a.0.as_str()).collect();
        assert_eq!(
            names,
            BTreeSet::from([
                "source_relations",
                "transport_relations",
                "rho_negation",
                "rho_aa"
            ])
        );
        assert!(arms.iter().all(|a| a.1 == arms[0].1));
        let rho = arms.iter().find(|a| a.0 == "rho_negation").unwrap();
        let aa = arms.iter().find(|a| a.0 == "rho_aa").unwrap();
        assert_eq!((rho.1, rho.2, rho.3), (aa.1, aa.2, aa.3));
    }
    println!("independent_check=pass rows={count} measured={measured} paired_workloads={} four_sizes=pass source_and_destination_scalar_checks=pass relation_collection_rank_recovery_descent=pass phase_partitions=pass aa_count_identity=pass",per.len());
}
