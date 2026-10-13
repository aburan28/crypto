//! Non-isomorphic isogeny and rational kernel control: see nonisomorphic_pair_20261013/PROTOCOL.md.
use crypto_lib::cryptanalysis::ic_boundary::{
    rho_walk_with, CountedGroup, GroupOps, NegationClasses, PrimeCurve, PrimePoint, RhoWalk,
};
use std::collections::{BTreeMap, BTreeSet};

const RELATION_SEED: u64 = 0x5245_4c41_5449_4f4e;
const DESCENT_SEED: u64 = 0x4445_5343_454e_5400;
const RHO_SEED: u64 = 0x5248_4f4e_4547_0000;
const TRIAL_CAP: u64 = 10_000;
const P: u64 = 617;
const B: u64 = 1;
const R: u64 = 103;
const COFACTOR: u64 = 6;
const FB_SIZE: usize = 12;

#[derive(Clone, Copy, Default)]
struct MapOps {
    muls: u64,
    sqrs: u64,
    invs: u64,
    evaluations: u64,
}

#[derive(Default)]
struct Solve {
    setup: GroupOps,
    collection: GroupOps,
    descent: GroupOps,
    verification: GroupOps,
    rho: GroupOps,
    lifts: u64,
    pair_entries: u64,
    relation_trials: u64,
    relation_lookups: u64,
    independent_rows: u64,
    row_mul_subs: u64,
    row_inversions: u64,
    descent_trials: u64,
    descent_lookups: u64,
    recovered: Option<u64>,
    status: &'static str,
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

fn mul(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn pow(mut a: u64, mut n: u64, p: u64) -> u64 {
    let mut v = 1;
    while n != 0 {
        if n & 1 != 0 {
            v = mul(v, a, p);
        }
        a = mul(a, a, p);
        n >>= 1;
    }
    v
}

fn sub(a: u64, b: u64, p: u64) -> u64 {
    (a + p - b % p) % p
}

fn phi(point: PrimePoint, p: u64, b: u64, tally: &mut MapOps) -> PrimePoint {
    if point.infinity {
        tally.evaluations += 1;
        return point;
    }
    if point.x == 0 {
        tally.evaluations += 1;
        return PrimePoint::INFINITY;
    }
    let x2 = mul(point.x, point.x, p);
    let x3 = mul(x2, point.x, p);
    let ix = pow(point.x, p - 2, p);
    let ix2 = mul(ix, ix, p);
    let ix3 = mul(ix2, ix, p);
    let nx = mul((x3 + 4 * b) % p, ix2, p);
    let ny = mul(mul(point.y, sub(x3, (8 * b) % p, p), p), ix3, p);
    tally.muls += 5;
    tally.sqrs += 2;
    tally.invs += 1;
    tally.evaluations += 1;
    PrimePoint::affine(nx, ny)
}

fn affine_points(curve: &PrimeCurve) -> Vec<PrimePoint> {
    let mut pts = Vec::new();
    for x in 0..curve.p {
        let rhs = curve.rhs(x);
        if let Some(y) = curve.sqrt(rhs) {
            pts.push(PrimePoint::affine(x, y));
            if y != 0 {
                pts.push(PrimePoint::affine(x, curve.p - y));
            }
        }
    }
    pts.sort_by_key(|pt| (pt.x, pt.y));
    pts
}

fn prime(n: u64) -> bool {
    n >= 2 && (2..).take_while(|d| d * d <= n).all(|d| n % d != 0)
}

fn check_map(source: &PrimeCurve, target: &PrimeCurve) -> (PrimePoint, u64) {
    assert!(prime(P) && prime(R));
    let first = (101..)
        .filter(|&r| prime(r) && prime(6 * r - 1))
        .map(|r| (6 * r - 1, r))
        .next()
        .expect("bounded candidate exists");
    assert_eq!(first, (P, R));
    let src = affine_points(source);
    let dst = affine_points(target);
    assert_eq!(src.len() as u64 + 1, P + 1);
    assert_eq!(dst.len() as u64 + 1, P + 1);
    assert_eq!(source.p, P);
    assert_eq!(source.b, B);
    assert_eq!(target.b, sub(0, 27 * B, P));
    assert!(!(1..P).any(|u| pow(u, 6, P) == target.b));
    let mut images = BTreeMap::<u64, u64>::new();
    *images.entry(target.key(&PrimePoint::INFINITY)).or_default() += 1;
    let mut rational_kernel = vec![PrimePoint::INFINITY];
    for &pt in &src {
        let q = phi(pt, source.p, source.b, &mut MapOps::default());
        assert!(target.is_on_curve(q));
        *images.entry(target.key(&q)).or_default() += 1;
        if q.infinity {
            rational_kernel.push(pt);
        }
    }
    assert_eq!(
        rational_kernel,
        vec![
            PrimePoint::INFINITY,
            PrimePoint::affine(0, 1),
            PrimePoint::affine(0, P - 1)
        ]
    );
    assert_eq!(images.len() as u64, (P + 1) / 3);
    assert!(images.values().all(|&multiplicity| multiplicity == 3));
    let mut all = vec![PrimePoint::INFINITY];
    all.extend(src.iter().copied());
    let mut checks = 0;
    for &a in &all {
        for &b in &all {
            let mut ops = GroupOps::default();
            let sum = source.add(&mut ops, a, b);
            let left = phi(sum, source.p, source.b, &mut MapOps::default());
            let pa = phi(a, source.p, source.b, &mut MapOps::default());
            let pb = phi(b, source.p, source.b, &mut MapOps::default());
            assert_eq!(left, target.add(&mut ops, pa, pb));
            checks += 1;
        }
    }
    let mut ops = GroupOps::default();
    let generator = src
        .iter()
        .find_map(|&point| {
            let projected = source.mul(&mut ops, point, COFACTOR);
            (!projected.infinity).then_some(projected)
        })
        .expect("nonzero subgroup point");
    assert!(source.mul(&mut ops, generator, R).infinity);
    let image = phi(generator, P, B, &mut MapOps::default());
    assert!(!image.infinity);
    assert!(target.mul(&mut ops, image, R).infinity);
    (generator, checks)
}

fn targets(curve: &PrimeCurve) -> Vec<PrimePoint> {
    let mut stream = SplitMix64(0x2026_1013_0001_0000 + curve.p);
    let mut pts = Vec::new();
    let mut seen = BTreeSet::new();
    let mut ops = GroupOps::default();
    while pts.len() < 8 {
        let x = stream.next() % curve.p;
        let sign = stream.next() & 1;
        if let Some(y0) = curve.sqrt(curve.rhs(x)) {
            let y = y0.min(curve.p - y0);
            let pt = PrimePoint::affine(x, if sign == 1 && y != 0 { curve.p - y } else { y });
            let projected = curve.mul(&mut ops, pt, COFACTOR);
            if !projected.infinity && seen.insert(curve.key(&projected)) {
                pts.push(projected);
            }
        }
    }
    pts
}

fn factor_base(curve: &PrimeCurve, size: usize, out: &mut Solve) -> Vec<PrimePoint> {
    let mut fb = Vec::new();
    let mut seen = BTreeSet::new();
    for x in 0..curve.p {
        out.lifts += 1;
        if let Some(y0) = curve.sqrt(curve.rhs(x)) {
            let y = y0.min(curve.p - y0);
            let projected = curve.mul(&mut out.setup, PrimePoint::affine(x, y), COFACTOR);
            if !projected.infinity && seen.insert(curve.key(&projected)) {
                fb.push(projected);
                if fb.len() == size {
                    return fb;
                }
            }
        }
    }
    panic!("insufficient factor base points")
}

type PairTable = BTreeMap<u64, Vec<(usize, usize)>>;
fn pair_table(curve: &PrimeCurve, fb: &[PrimePoint], out: &mut Solve) -> PairTable {
    let mut table = BTreeMap::new();
    for i in 0..fb.len() {
        for j in i..fb.len() {
            let sum = curve.add(&mut out.setup, fb[i], fb[j]);
            table
                .entry(curve.key(&sum))
                .or_insert_with(Vec::new)
                .push((i, j));
            out.pair_entries += 1;
        }
    }
    table
}

struct Row {
    coeff: Vec<u64>,
    rhs: u64,
}

fn insert_row(rows: &mut [Option<Row>], mut row: Row, p: u64, out: &mut Solve) -> bool {
    for (col, existing) in rows.iter_mut().enumerate() {
        let c = row.coeff[col];
        if c == 0 {
            continue;
        }
        if let Some(pivot) = existing {
            for j in col..row.coeff.len() {
                row.coeff[j] = sub(row.coeff[j], mul(c, pivot.coeff[j], p), p);
                out.row_mul_subs += 1;
            }
            row.rhs = sub(row.rhs, mul(c, pivot.rhs, p), p);
            out.row_mul_subs += 1;
        } else {
            let inverse = pow(c, p - 2, p);
            out.row_inversions += 1;
            for j in col..row.coeff.len() {
                row.coeff[j] = mul(row.coeff[j], inverse, p);
                out.row_mul_subs += 1;
            }
            row.rhs = mul(row.rhs, inverse, p);
            out.row_mul_subs += 1;
            *existing = Some(row);
            out.independent_rows += 1;
            return true;
        }
    }
    assert_eq!(row.rhs, 0, "inconsistent exact relation");
    false
}

fn recover_logs(rows: &[Option<Row>], p: u64, out: &mut Solve) -> Vec<u64> {
    let mut logs = vec![0; rows.len()];
    for i in (0..rows.len()).rev() {
        let row = rows[i].as_ref().expect("full rank");
        let mut rhs = row.rhs;
        for (j, value) in logs.iter().enumerate().skip(i + 1) {
            rhs = sub(rhs, mul(row.coeff[j], *value, p), p);
            out.row_mul_subs += 1;
        }
        logs[i] = rhs;
    }
    logs
}

fn relation_solve(
    curve: &PrimeCurve,
    generator: PrimePoint,
    target: PrimePoint,
    r: u64,
    fb_size: usize,
    index: usize,
) -> Solve {
    let mut out = Solve::default();
    let fb = factor_base(curve, fb_size, &mut out);
    let table = pair_table(curve, &fb, &mut out);
    let mut rows: Vec<Option<Row>> = (0..fb_size).map(|_| None).collect();
    let key = curve.p ^ ((index as u64) << 16);
    let mut stream = SplitMix64(RELATION_SEED ^ key);
    for trial in 0..TRIAL_CAP {
        out.relation_trials = trial + 1;
        let k = 1 + stream.next() % (r - 1);
        let point = curve.mul(&mut out.collection, generator, k);
        out.relation_lookups += 1;
        if let Some(pairs) = table.get(&curve.key(&point)) {
            for &(i, j) in pairs {
                let mut coeff = vec![0; fb_size];
                coeff[i] += 1;
                coeff[j] += 1;
                insert_row(&mut rows, Row { coeff, rhs: k }, r, &mut out);
                if out.independent_rows as usize == fb_size {
                    break;
                }
            }
        }
        if out.independent_rows as usize == fb_size {
            break;
        }
    }
    if out.independent_rows as usize != fb_size {
        out.status = "collection_cap";
        return out;
    }
    let logs = recover_logs(&rows, r, &mut out);
    for (&point, &d) in fb.iter().zip(&logs) {
        if curve.mul(&mut out.verification, generator, d) != point {
            out.status = "factor_base_verification_failed";
            return out;
        }
    }
    let mut stream = SplitMix64(DESCENT_SEED ^ key);
    for trial in 0..TRIAL_CAP {
        out.descent_trials = trial + 1;
        let t = stream.next() % r;
        let tg = curve.mul(&mut out.descent, generator, t);
        let probe = curve.add(&mut out.descent, target, tg);
        out.descent_lookups += 1;
        if let Some(pairs) = table.get(&curve.key(&probe)) {
            let (i, j) = pairs[0];
            let d = sub((logs[i] + logs[j]) % r, t, r);
            if curve.mul(&mut out.verification, generator, d) != target {
                out.status = "target_verification_failed";
                return out;
            }
            out.recovered = Some(d);
            out.status = "verified";
            return out;
        }
    }
    out.status = "descent_cap";
    out
}

fn print_solve(
    p: u64,
    b: u64,
    r: u64,
    i: usize,
    arm: &str,
    source_q: PrimePoint,
    transported_q: PrimePoint,
    solved: &Solve,
    maps: MapOps,
    rho: Option<&crypto_lib::cryptanalysis::ic_boundary::RhoResult>,
    rho_class: (u64, u64, u64, u64, u64),
) {
    let mut total = GroupOps::default();
    for ops in [
        solved.setup,
        solved.collection,
        solved.descent,
        solved.verification,
        solved.rho,
    ] {
        total.merge(ops);
    }
    let rc = |key: &str| rho.and_then(|v| v.counters.get(key)).copied().unwrap_or(0);
    let fields = vec![
        p.to_string(),
        b.to_string(),
        r.to_string(),
        i.to_string(),
        arm.to_string(),
        source_q.x.to_string(),
        source_q.y.to_string(),
        transported_q.x.to_string(),
        transported_q.y.to_string(),
        solved.status.to_string(),
        solved
            .recovered
            .map_or_else(|| "none".to_string(), |d| d.to_string()),
        solved.setup.gae().to_string(),
        solved.collection.gae().to_string(),
        solved.descent.gae().to_string(),
        solved.verification.gae().to_string(),
        solved.rho.gae().to_string(),
        total.gae().to_string(),
        total.field_ops.muls.to_string(),
        total.field_ops.sqrs.to_string(),
        total.field_ops.invs.to_string(),
        solved.lifts.to_string(),
        solved.pair_entries.to_string(),
        solved.relation_trials.to_string(),
        solved.relation_lookups.to_string(),
        solved.independent_rows.to_string(),
        solved.row_mul_subs.to_string(),
        solved.row_inversions.to_string(),
        solved.descent_trials.to_string(),
        solved.descent_lookups.to_string(),
        maps.evaluations.to_string(),
        maps.muls.to_string(),
        maps.sqrs.to_string(),
        maps.invs.to_string(),
        solved.setup.adds.to_string(),
        solved.setup.doubles.to_string(),
        solved.collection.adds.to_string(),
        solved.collection.doubles.to_string(),
        solved.descent.adds.to_string(),
        solved.descent.doubles.to_string(),
        solved.verification.adds.to_string(),
        solved.verification.doubles.to_string(),
        solved.rho.adds.to_string(),
        solved.rho.doubles.to_string(),
        rho.map_or(0, |v| v.steps).to_string(),
        rho.map_or(0, |v| v.walks).to_string(),
        rho.map_or(0, |v| v.distinguished_points).to_string(),
        rc("canonicalisations_uncharged").to_string(),
        rho_class.0.to_string(),
        rho_class.1.to_string(),
        rho_class.2.to_string(),
        rho_class.3.to_string(),
        rho_class.4.to_string(),
    ];
    assert_eq!(
        fields.len(),
        HEADER.split('\t').count(),
        "header/row width mismatch"
    );
    println!("{}", fields.join("\t"));
}

const HEADER: &str = "p\tb\tr\ttarget_index\tarm\tqx\tqy\tmapped_qx\tmapped_qy\tstatus\tlog\tsetup_gae\tcollection_gae\tdescent_gae\tverification_gae\trho_gae\ttotal_gae\tfield_muls\tfield_sqrs\tfield_invs\tfb_lifts\tpair_entries\trelation_trials\trelation_lookups\tindependent_rows\trow_mul_subs\trow_inversions\tdescent_trials\tdescent_lookups\tmap_evaluations\tmap_muls\tmap_sqrs\tmap_invs\tsetup_adds\tsetup_doubles\tcollection_adds\tcollection_doubles\tdescent_adds\tdescent_doubles\tverification_adds\tverification_doubles\trho_adds\trho_doubles\trho_steps\trho_walks\trho_dp\trho_unpriced_canonicalisations\trho_canon_field_muls\trho_canon_point_negations\trho_lambda\trho_attempts\trho_class_setup_field_muls";

fn main() {
    println!("{HEADER}");
    let p = P;
    let b = B;
    let r = R;
    let target_b = sub(0, 27 * b, p);
    let source = PrimeCurve { p, a: 0, b };
    let destination = PrimeCurve {
        p,
        a: 0,
        b: target_b,
    };
    let (generator, checks) = check_map(&source, &destination);
    let generator_image = phi(generator, p, b, &mut MapOps::default());
    eprintln!("map_check p={p} b={b} target_b={target_b} order={} subgroup={r} cofactor={COFACTOR} rational_kernel=3 image_size={} sixth_power_isomorphism=false homomorphism_pairs={checks} generator=({}, {}) image_generator=({}, {}) status=pass", p + 1, (p + 1) / 3, generator.x, generator.y, generator_image.x, generator_image.y);
    let mut failures = 0;
    for (i, q) in targets(&source).into_iter().enumerate() {
        let mut map_ops = MapOps::default();
        let g2 = phi(generator, p, b, &mut map_ops);
        let q2 = phi(q, p, b, &mut map_ops);
        assert_eq!(g2, generator_image);
        assert!(destination.is_on_curve(q2) && !q2.infinity);
        eprintln!(
            "target index={i} source=({}, {}) image=({}, {})",
            q.x, q.y, q2.x, q2.y
        );
        let src = relation_solve(&source, generator, q, r, FB_SIZE, i);
        print_solve(
            p,
            b,
            r,
            i,
            "source_relations",
            q,
            q2,
            &src,
            MapOps::default(),
            None,
            (0, 0, 0, 0, 0),
        );
        let mut dst = relation_solve(&destination, g2, q2, r, FB_SIZE, i);
        if let Some(d) = dst.recovered {
            if source.mul(&mut dst.verification, generator, d) != q {
                dst.status = "source_crosscheck_failed";
            }
        }
        print_solve(
            p,
            b,
            r,
            i,
            "transport_relations",
            q,
            q2,
            &dst,
            map_ops,
            None,
            (0, 0, 0, 0, 0),
        );
        let rho = rho_walk_with(
            &source,
            &NegationClasses { r },
            generator,
            q,
            r,
            RHO_SEED ^ p ^ ((i as u64) << 16),
            100_000,
            RhoWalk::negation(),
        );
        let mut verify = GroupOps::default();
        let valid = rho
            .recovered
            .is_some_and(|d| source.mul(&mut verify, generator, d) == q);
        let rho_solve = Solve {
            verification: verify,
            rho: rho.group_ops,
            recovered: rho.recovered,
            status: if valid { "verified" } else { "rho_failed" },
            ..Solve::default()
        };
        print_solve(
            p,
            b,
            r,
            i,
            "rho_negation",
            q,
            q2,
            &rho_solve,
            MapOps::default(),
            Some(&rho),
            (0, 0, 0, 0, 0),
        );
        if src.status != "verified"
            || dst.status != "verified"
            || !valid
            || src.recovered != dst.recovered
            || src.recovered != rho.recovered
        {
            eprintln!(
                "workload_failed index={i} source={} destination={} rho={} values={:?}/{:?}/{:?}",
                src.status, dst.status, valid, src.recovered, dst.recovered, rho.recovered
            );
            failures += 1;
        }
    }
    eprintln!("workloads=8 arms=24 failed_workloads={failures}");
    if failures != 0 {
        std::process::exit(2);
    }
    eprintln!("all_workloads_verified=8 arms_verified=24");
}
