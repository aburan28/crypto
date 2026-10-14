//! Fully priced growing-size control: see full_cost_20261013/PROTOCOL.md.
use crypto_lib::cryptanalysis::ic_boundary::{
    rho_walk_with, CountedGroup, GroupOps, NegationClasses, PrimeCurve, PrimePoint, RhoWalk,
};
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

const RELATION_SEED: u64 = 0x5245_4c41_5449_4f4e;
const DESCENT_SEED: u64 = 0x4445_5343_454e_5400;
const RHO_SEED: u64 = 0x5248_4f4e_4547_0000;
const TRIAL_CAP: u64 = 100_000;
const B: u64 = 1;
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
    setup_ns: u64,
    collection_ns: u64,
    recovery_ns: u64,
    descent_ns: u64,
    verification_ns: u64,
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

fn selected_pairs() -> Vec<(u64, u64)> {
    let bands = [(101, 400), (401, 1600), (1601, 6400), (6401, 25_600)];
    bands
        .into_iter()
        .map(|(lo, hi)| {
            let r = (lo..=hi)
                .find(|&r| prime(r) && prime(6 * r - 1))
                .expect("candidate in frozen band");
            (6 * r - 1, r)
        })
        .collect()
}

fn check_map(source: &PrimeCurve, target: &PrimeCurve, r: u64) -> (PrimePoint, u64) {
    let p = source.p;
    assert!(prime(p) && prime(r) && p == 6 * r - 1);
    let src = affine_points(source);
    let dst = affine_points(target);
    assert_eq!(src.len() as u64 + 1, p + 1);
    assert_eq!(dst.len() as u64 + 1, p + 1);
    assert_eq!(source.b, B);
    assert_eq!(target.b, sub(0, 27 * B, p));
    assert!(!(1..p).any(|u| pow(u, 6, p) == target.b));
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
            PrimePoint::affine(0, p - 1)
        ]
    );
    assert_eq!(images.len() as u64, (p + 1) / 3);
    assert!(images.values().all(|&multiplicity| multiplicity == 3));
    let mut all = vec![PrimePoint::INFINITY];
    all.extend(src.iter().copied());
    let mut checks = 0;
    let mut pair_stream = SplitMix64(0x5041_4952_4348_4543 ^ p);
    let pairs: Vec<(PrimePoint, PrimePoint)> = if p == 617 {
        all.iter()
            .flat_map(|&a| all.iter().map(move |&b| (a, b)))
            .collect()
    } else {
        (0..256)
            .map(|_| {
                let a = all[(pair_stream.next() % all.len() as u64) as usize];
                let b = all[(pair_stream.next() % all.len() as u64) as usize];
                (a, b)
            })
            .collect()
    };
    for (a, b) in pairs {
        let mut ops = GroupOps::default();
        let sum = source.add(&mut ops, a, b);
        let left = phi(sum, source.p, source.b, &mut MapOps::default());
        let pa = phi(a, source.p, source.b, &mut MapOps::default());
        let pb = phi(b, source.p, source.b, &mut MapOps::default());
        assert_eq!(left, target.add(&mut ops, pa, pb));
        checks += 1;
    }
    let mut ops = GroupOps::default();
    let generator = src
        .iter()
        .find_map(|&point| {
            let projected = source.mul(&mut ops, point, COFACTOR);
            (!projected.infinity).then_some(projected)
        })
        .expect("nonzero subgroup point");
    assert!(source.mul(&mut ops, generator, r).infinity);
    let image = phi(generator, p, B, &mut MapOps::default());
    assert!(!image.infinity);
    assert!(target.mul(&mut ops, image, r).infinity);
    for k in 0..r {
        let source_multiple = source.mul(&mut ops, generator, k);
        let destination_multiple = target.mul(&mut ops, image, k);
        assert_eq!(
            phi(source_multiple, p, B, &mut MapOps::default()),
            destination_multiple
        );
    }
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
    round: u64,
) -> Solve {
    let mut out = Solve::default();
    let phase = Instant::now();
    let fb = factor_base(curve, fb_size, &mut out);
    let table = pair_table(curve, &fb, &mut out);
    out.setup_ns = phase.elapsed().as_nanos() as u64;
    let mut rows: Vec<Option<Row>> = (0..fb_size).map(|_| None).collect();
    let key = curve.p ^ ((index as u64) << 16) ^ (round << 32);
    let mut stream = SplitMix64(RELATION_SEED ^ key);
    let collection_clock = Instant::now();
    let mut row_ns = 0u64;
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
                let row_clock = Instant::now();
                insert_row(&mut rows, Row { coeff, rhs: k }, r, &mut out);
                row_ns += row_clock.elapsed().as_nanos() as u64;
                if out.independent_rows as usize == fb_size {
                    break;
                }
            }
        }
        if out.independent_rows as usize == fb_size {
            break;
        }
    }
    out.collection_ns = (collection_clock.elapsed().as_nanos() as u64).saturating_sub(row_ns);
    out.recovery_ns = row_ns;
    if out.independent_rows as usize != fb_size {
        out.status = "collection_cap";
        return out;
    }
    let recovery_clock = Instant::now();
    let logs = recover_logs(&rows, r, &mut out);
    out.recovery_ns += recovery_clock.elapsed().as_nanos() as u64;
    let verification_clock = Instant::now();
    for (&point, &d) in fb.iter().zip(&logs) {
        if curve.mul(&mut out.verification, generator, d) != point {
            out.status = "factor_base_verification_failed";
            return out;
        }
    }
    out.verification_ns = verification_clock.elapsed().as_nanos() as u64;
    let mut stream = SplitMix64(DESCENT_SEED ^ key);
    let descent_clock = Instant::now();
    for trial in 0..TRIAL_CAP {
        out.descent_trials = trial + 1;
        let t = stream.next() % r;
        let tg = curve.mul(&mut out.descent, generator, t);
        let probe = curve.add(&mut out.descent, target, tg);
        out.descent_lookups += 1;
        if let Some(pairs) = table.get(&curve.key(&probe)) {
            let (i, j) = pairs[0];
            let d = sub((logs[i] + logs[j]) % r, t, r);
            out.descent_ns = descent_clock.elapsed().as_nanos() as u64;
            let verify_clock = Instant::now();
            if curve.mul(&mut out.verification, generator, d) != target {
                out.status = "target_verification_failed";
                return out;
            }
            out.verification_ns += verify_clock.elapsed().as_nanos() as u64;
            out.recovered = Some(d);
            out.status = "verified";
            return out;
        }
    }
    out.descent_ns = descent_clock.elapsed().as_nanos() as u64;
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
    round: u64,
    warmup: bool,
    wall_ns: u64,
    map_ns: u64,
    rho_ns: u64,
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
        round.to_string(),
        u8::from(warmup).to_string(),
        wall_ns.to_string(),
        map_ns.to_string(),
        solved.setup_ns.to_string(),
        solved.collection_ns.to_string(),
        solved.recovery_ns.to_string(),
        solved.descent_ns.to_string(),
        solved.verification_ns.to_string(),
        rho_ns.to_string(),
        wall_ns
            .saturating_sub(
                map_ns
                    + solved.setup_ns
                    + solved.collection_ns
                    + solved.recovery_ns
                    + solved.descent_ns
                    + solved.verification_ns
                    + rho_ns,
            )
            .to_string(),
        max_rss_kib().to_string(),
    ];
    assert_eq!(
        fields.len(),
        HEADER.split('\t').count(),
        "header/row width mismatch"
    );
    println!("{}", fields.join("\t"));
}

const HEADER: &str = "p\tb\tr\ttarget_index\tarm\tqx\tqy\tmapped_qx\tmapped_qy\tstatus\tlog\tsetup_gae\tcollection_gae\tdescent_gae\tverification_gae\trho_gae\ttotal_gae\tfield_muls\tfield_sqrs\tfield_invs\tfb_lifts\tpair_entries\trelation_trials\trelation_lookups\tindependent_rows\trow_mul_subs\trow_inversions\tdescent_trials\tdescent_lookups\tmap_evaluations\tmap_muls\tmap_sqrs\tmap_invs\tsetup_adds\tsetup_doubles\tcollection_adds\tcollection_doubles\tdescent_adds\tdescent_doubles\tverification_adds\tverification_doubles\trho_adds\trho_doubles\trho_steps\trho_walks\trho_dp\trho_unpriced_canonicalisations\trho_canon_field_muls\trho_canon_point_negations\trho_lambda\trho_attempts\trho_class_setup_field_muls\tround\twarmup\twall_ns\tmap_ns\tsetup_ns\tcollection_ns\trecovery_ns\tdescent_ns\tverification_ns\trho_ns\tother_ns\tprocess_peak_rss_kib";

fn max_rss_kib() -> u64 {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    // SAFETY: getrusage initializes the entire output structure on success.
    let ok = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
    assert_eq!(ok, 0, "getrusage failed");
    // SAFETY: the successful call above initialized usage.
    let usage = unsafe { usage.assume_init() };
    #[cfg(target_os = "macos")]
    {
        (usage.ru_maxrss as u64) / 1024
    }
    #[cfg(not(target_os = "macos"))]
    {
        usage.ru_maxrss as u64
    }
}

fn main() {
    println!("{HEADER}");
    let pairs = selected_pairs();
    assert_eq!(pairs[0], (617, 103));
    let mut failures = 0;
    for (p, r) in pairs {
        let preflight_clock = Instant::now();
        let source = PrimeCurve { p, a: 0, b: B };
        let destination = PrimeCurve {
            p,
            a: 0,
            b: sub(0, 27 * B, p),
        };
        let (generator, checks) = check_map(&source, &destination, r);
        let generator_image = phi(generator, p, B, &mut MapOps::default());
        let qs = targets(&source);
        let images: Vec<_> = qs
            .iter()
            .map(|&q| phi(q, p, B, &mut MapOps::default()))
            .collect();
        let preflight_ns = preflight_clock.elapsed().as_nanos() as u64;
        eprintln!("pair p={p} r={r} source_b=1 destination_b={} order={} kernel=3 image_size={} nonisomorphic=true homomorphism_pairs={checks} generator=({}, {}) image_generator=({}, {}) preflight_ns={preflight_ns} status=pass", destination.b, p+1, (p+1)/3, generator.x, generator.y, generator_image.x, generator_image.y);
        for (i, (&q, &q2)) in qs.iter().zip(&images).enumerate() {
            assert!(destination.is_on_curve(q2) && !q2.infinity);
            eprintln!(
                "target p={p} index={i} source=({}, {}) image=({}, {})",
                q.x, q.y, q2.x, q2.y
            );
        }
        for round in 0..=5u64 {
            for (i, (&q, &q2)) in qs.iter().zip(&images).enumerate() {
                let mut recovered = [None; 4];
                let mut verified = [false; 4];
                for offset in 0..4 {
                    let arm = (offset + (round as usize + i) % 4) % 4;
                    let clock = Instant::now();
                    match arm {
                        0 | 1 => {
                            let mut map_ops = MapOps::default();
                            let map_clock = Instant::now();
                            let (curve, g, target, name, map_ns) = if arm == 0 {
                                (&source, generator, q, "source_relations", 0)
                            } else {
                                let g = phi(generator, p, B, &mut map_ops);
                                let target = phi(q, p, B, &mut map_ops);
                                assert_eq!(g, generator_image);
                                assert_eq!(target, q2);
                                (
                                    &destination,
                                    g,
                                    target,
                                    "transport_relations",
                                    map_clock.elapsed().as_nanos() as u64,
                                )
                            };
                            let mut solved = relation_solve(curve, g, target, r, FB_SIZE, i, round);
                            if arm == 1 {
                                if let Some(d) = solved.recovered {
                                    let check_clock = Instant::now();
                                    if source.mul(&mut solved.verification, generator, d) != q {
                                        solved.status = "source_crosscheck_failed";
                                    }
                                    solved.verification_ns +=
                                        check_clock.elapsed().as_nanos() as u64;
                                }
                            }
                            let wall_ns = clock.elapsed().as_nanos() as u64;
                            recovered[arm] = solved.recovered;
                            verified[arm] = solved.status == "verified";
                            print_solve(
                                p,
                                B,
                                r,
                                i,
                                name,
                                q,
                                q2,
                                &solved,
                                map_ops,
                                None,
                                (0, 0, 0, 0, 0),
                                round,
                                round == 0,
                                wall_ns,
                                map_ns,
                                0,
                            );
                        }
                        2 | 3 => {
                            let rho_clock = Instant::now();
                            let rho = rho_walk_with(
                                &source,
                                &NegationClasses { r },
                                generator,
                                q,
                                r,
                                RHO_SEED ^ p ^ ((i as u64) << 16) ^ (round << 32),
                                1_000_000,
                                RhoWalk::negation(),
                            );
                            let rho_ns = rho_clock.elapsed().as_nanos() as u64;
                            let verify_clock = Instant::now();
                            let mut verify = GroupOps::default();
                            let valid = rho
                                .recovered
                                .is_some_and(|d| source.mul(&mut verify, generator, d) == q);
                            let verification_ns = verify_clock.elapsed().as_nanos() as u64;
                            let solved = Solve {
                                verification: verify,
                                verification_ns,
                                rho: rho.group_ops,
                                recovered: rho.recovered,
                                status: if valid { "verified" } else { "rho_failed" },
                                ..Solve::default()
                            };
                            let wall_ns = clock.elapsed().as_nanos() as u64;
                            recovered[arm] = solved.recovered;
                            verified[arm] = valid;
                            print_solve(
                                p,
                                B,
                                r,
                                i,
                                if arm == 2 { "rho_negation" } else { "rho_aa" },
                                q,
                                q2,
                                &solved,
                                MapOps::default(),
                                Some(&rho),
                                (0, 0, 0, 0, 0),
                                round,
                                round == 0,
                                wall_ns,
                                0,
                                rho_ns,
                            );
                        }
                        _ => unreachable!(),
                    }
                }
                if !verified.iter().all(|&v| v) || !recovered.iter().all(|&d| d == recovered[0]) {
                    eprintln!("workload_failed p={p} round={round} index={i} verified={verified:?} values={recovered:?}");
                    failures += 1;
                }
            }
        }
        eprintln!("pair_complete p={p} r={r} measured_rows=160 warmup_rows=32 cumulative_failures={failures}");
    }
    eprintln!("sizes=4 workloads=32 measured_rows=640 warmup_rows=128 failed_workloads={failures}");
    if failures != 0 {
        std::process::exit(2);
    }
    eprintln!("all_workloads_verified=32 measured_arms_verified=640");
}
