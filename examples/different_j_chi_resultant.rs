//! A direct chi-coordinate chord resultant on the frozen different-j cell.
//! This correctness gate does not use a pair-sum table or scan base pairs
//! inside the polynomial oracle. It tests every point of the order-229 group.
#![recursion_limit = "256"]

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;
use std::time::Instant;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ext_curve::{
    take_field_counters, ExtCurve, ExtField, ExtPoint, Fp2,
};
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::ic_boundary::rho_reference_negation;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

type El = [u64; 2];
type Point = ExtPoint<El>;
type Curve = ExtCurve<Fp2>;
type Poly = Vec<El>;

fn trim(f: &Fp2, p: &mut Poly) {
    while p.last().is_some_and(|&c| c == f.zero()) {
        p.pop();
    }
}

fn add(f: &Fp2, a: &Poly, b: &Poly) -> Poly {
    let mut out = vec![f.zero(); a.len().max(b.len())];
    for (i, &v) in a.iter().enumerate() {
        out[i] = f.add(out[i], v);
    }
    for (i, &v) in b.iter().enumerate() {
        out[i] = f.add(out[i], v);
    }
    trim(f, &mut out);
    out
}

fn sub(f: &Fp2, a: &Poly, b: &Poly) -> Poly {
    let mut out = vec![f.zero(); a.len().max(b.len())];
    for (i, &v) in a.iter().enumerate() {
        out[i] = f.add(out[i], v);
    }
    for (i, &v) in b.iter().enumerate() {
        out[i] = f.sub(out[i], v);
    }
    trim(f, &mut out);
    out
}

fn scale(f: &Fp2, a: &Poly, k: El) -> Poly {
    let mut out: Poly = a.iter().map(|&v| f.mul(v, k)).collect();
    trim(f, &mut out);
    out
}

fn mul(f: &Fp2, a: &Poly, b: &Poly) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![f.zero(); a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            out[i + j] = f.add(out[i + j], f.mul(x, y));
        }
    }
    trim(f, &mut out);
    out
}

fn div_rem(f: &Fp2, a: &Poly, b: &Poly) -> (Poly, Poly) {
    assert!(!b.is_empty());
    let mut rem = a.clone();
    let mut quo = vec![f.zero(); a.len().saturating_sub(b.len()) + 1];
    let inv_lead = f.inv(*b.last().unwrap()).unwrap();
    while rem.len() >= b.len() {
        let shift = rem.len() - b.len();
        let coeff = f.mul(*rem.last().unwrap(), inv_lead);
        quo[shift] = coeff;
        for (j, &term) in b.iter().enumerate() {
            rem[shift + j] = f.sub(rem[shift + j], f.mul(coeff, term));
        }
        trim(f, &mut rem);
    }
    trim(f, &mut quo);
    (quo, rem)
}

fn gcd(f: &Fp2, mut a: Poly, mut b: Poly) -> Poly {
    while !b.is_empty() {
        let rem = div_rem(f, &a, &b).1;
        a = b;
        b = rem;
    }
    if let Some(&lead) = a.last() {
        let inverse = f.inv(lead).unwrap();
        for coeff in &mut a {
            *coeff = f.mul(*coeff, inverse);
        }
    }
    a
}

fn mul_mod(f: &Fp2, a: &Poly, b: &Poly, modulus: &Poly) -> Poly {
    div_rem(f, &mul(f, a, b), modulus).1
}

fn pow_mod(f: &Fp2, a: &Poly, mut exponent: u64, modulus: &Poly) -> Poly {
    let mut out = vec![f.one()];
    let mut base = div_rem(f, a, modulus).1;
    while exponent > 0 {
        if exponent & 1 != 0 {
            out = mul_mod(f, &out, &base, modulus);
        }
        exponent >>= 1;
        if exponent > 0 {
            base = mul_mod(f, &base, &base, modulus);
        }
    }
    out
}

fn eval(f: &Fp2, a: &Poly, x: El) -> El {
    a.iter()
        .rev()
        .fold(f.zero(), |acc, &c| f.add(f.mul(acc, x), c))
}

fn root_multiplicity(f: &Fp2, a: &Poly, root: El) -> u64 {
    let mut q = a.clone();
    let mut count = 0;
    let linear = vec![f.neg(root), f.one()];
    while q.len() > 1 {
        let (next, rem) = div_rem(f, &q, &linear);
        if !rem.is_empty() {
            break;
        }
        q = next;
        count += 1;
    }
    count
}

#[derive(Default)]
struct Stats {
    resultants: u64,
    gcds: u64,
    split_attempts: u64,
    exhaustive_fallback_evaluations: u64,
    root_candidates: u64,
    chi_lift_scans: u64,
    exact_checks: u64,
    false_candidates: u64,
    vertical_cases: u64,
    tangent_cases: u64,
    horizontal_cases: u64,
    infinity_cases: u64,
    exceptional_chi_cases: u64,
    multiple_resultant_roots: u64,
}

#[derive(Default)]
struct S3Stats {
    quadratics: u64,
    roots: u64,
    lift_scans: u64,
    exact_checks: u64,
    false_candidates: u64,
    degenerate: u64,
    sqrt_failures: u64,
    infinity_cases: u64,
}

fn s3_coefficients(curve: &Curve, x: El, z: El) -> [El; 3] {
    let f=&curve.f;
    let x2=f.sqr(x);
    let z2=f.sqr(z);
    let leading=f.sqr(f.sub(z,x));
    let middle=f.neg(f.scale(f.add(
        f.add(f.mul(x,z2),f.mul(f.add(x2,curve.a),z)),
        f.add(f.mul(x,curve.a),f.scale(curve.b,2))),2));
    let constant=f.sub(
        f.sub(f.mul(x2,z2),
            f.scale(f.mul(f.add(f.mul(x,curve.a),f.scale(curve.b,2)),z),2)),
        f.sub(f.scale(f.mul(curve.b,x),4),f.sqr(curve.a)));
    [constant,middle,leading]
}

fn s3_pairs(curve: &Curve, signed: &[Point], target: Point,
    ops: &mut GroupOps, stats: &mut S3Stats) -> BTreeSet<(usize,usize)> {
    let f=&curve.f;
    let mut out=BTreeSet::new();
    if target.infinity {
        stats.infinity_cases+=1;
        for (i,&p) in signed.iter().enumerate() {
            for (j,&q) in signed.iter().enumerate().skip(i) {
                stats.exact_checks+=1;
                if curve.add(ops,p,q).infinity { out.insert((i,j)); }
            }
        }
        return out;
    }
    for (i,&p) in signed.iter().enumerate() {
        let [c,b,a]=s3_coefficients(curve,p.x,target.x);
        stats.quadratics+=1;
        let mut roots=Vec::new();
        if a==f.zero() {
            stats.degenerate+=1;
            if b!=f.zero() {
                roots.push(f.mul(f.neg(c),f.inv(b).unwrap()));
            } else if c==f.zero() {
                roots.extend(signed.iter().map(|q|q.x));
            }
        } else {
            let discriminant=f.sub(f.sqr(b),f.scale(f.mul(a,c),4));
            if let Some(square_root)=f.sqrt(discriminant) {
                let inv=f.inv(f.scale(a,2)).unwrap();
                roots.push(f.mul(f.add(f.neg(b),square_root),inv));
                roots.push(f.mul(f.sub(f.neg(b),square_root),inv));
            } else {
                stats.sqrt_failures+=1;
            }
        }
        roots.sort_unstable();
        roots.dedup();
        stats.roots+=roots.len() as u64;
        for root in roots {
            for (j,&q) in signed.iter().enumerate() {
                stats.lift_scans+=1;
                if q.x==root {
                    stats.exact_checks+=1;
                    if curve.add(ops,p,q)==target {
                        out.insert((i.min(j),i.max(j)));
                    } else {
                        stats.false_candidates+=1;
                    }
                }
            }
        }
    }
    out
}

fn all_roots(f: &Fp2, polynomial: &Poly, stats: &mut Stats) -> Vec<El> {
    if polynomial.len() <= 1 {
        return Vec::new();
    }
    let x = vec![f.zero(), f.one()];
    let rational = gcd(
        f,
        polynomial.clone(),
        sub(f, &pow_mod(f, &x, f.order(), polynomial), &x),
    );
    stats.gcds += 1;
    fn split(f: &Fp2, a: Poly, stats: &mut Stats, output: &mut Vec<El>) {
        if a.len() <= 1 {
            return;
        }
        if a.len() == 2 {
            output.push(f.neg(f.mul(a[0], f.inv(a[1]).unwrap())));
            return;
        }
        if a.len() == 3 {
            let four_ac = f.scale(f.mul(a[2], a[0]), 4);
            let disc = f.sub(f.sqr(a[1]), four_ac);
            if let Some(root_disc) = f.sqrt(disc) {
                let inv = f.inv(f.scale(a[2], 2)).unwrap();
                output.push(f.mul(f.add(f.neg(a[1]), root_disc), inv));
                output.push(f.mul(f.sub(f.neg(a[1]), root_disc), inv));
                return;
            }
        }
        if a.len() == 4 {
            for shift in 0..128u64 {
                stats.split_attempts += 1;
                let linear = vec![[shift % f.p, 1], f.one()];
                let character = pow_mod(f, &linear, (f.order() - 1) / 2, &a);
                for test in [
                    sub(f, &character, &vec![f.one()]),
                    add(f, &character, &vec![f.one()]),
                    character,
                ] {
                    let divisor = gcd(f, a.clone(), test);
                    stats.gcds += 1;
                    if divisor.len() > 1 && divisor.len() < a.len() {
                        let (quotient, rem) = div_rem(f, &a, &divisor);
                        assert!(rem.is_empty());
                        split(f, divisor, stats, output);
                        split(f, quotient, stats, output);
                        return;
                    }
                }
            }
        }
        // Complete, explicitly charged fallback for a failed deterministic split.
        for packed in 0..f.order() {
            let z = [packed % f.p, packed / f.p];
            stats.exhaustive_fallback_evaluations += 1;
            if eval(f, &a, z) == f.zero() {
                output.push(z);
            }
        }
    }
    let mut roots = Vec::new();
    split(f, rational, stats, &mut roots);
    roots.sort_unstable();
    roots.dedup();
    assert!(roots.iter().all(|&z| eval(f, polynomial, z) == f.zero()));
    roots
}

fn chi(curve: &Curve, point: Point) -> Option<El> {
    if point.infinity {
        return None;
    }
    let f = &curve.f;
    f.inv(f.sub(point.y, f.one()))
        .map(|inv| f.mul(f.add(point.y, f.one()), inv))
}

fn chi_base(curve: &Curve, signed: &[Point]) -> Poly {
    let f = &curve.f;
    let mut base = vec![f.one()];
    for &point in signed {
        if let Some(z) = chi(curve, point) {
            base = mul(f, &base, &vec![f.neg(z), f.one()]);
        }
    }
    base
}

fn chord_resultant(curve: &Curve, point: Point, slope: El) -> Poly {
    let f = &curve.f;
    let d = vec![f.neg(f.one()), f.one()];
    let n = vec![f.one(), f.one()];
    let m3 = f.mul(f.sqr(slope), slope);
    let l = add(f, &scale(f, &d, f.sub(f.mul(slope, point.x), point.y)), &n);
    let d2 = mul(f, &d, &d);
    let d3 = mul(f, &d2, &d);
    let l3 = mul(f, &mul(f, &l, &l), &l);
    let ax = scale(f, &mul(f, &d2, &l), f.mul(curve.a, f.sqr(slope)));
    let b = scale(f, &d3, f.mul(curve.b, m3));
    let y2 = scale(f, &mul(f, &d, &mul(f, &n, &n)), m3);
    sub(f, &add(f, &add(f, &l3, &ax), &b), &y2)
}

fn parse_point(v: &Value) -> Point {
    if v["infinity"] == true {
        return Point::infinity([0, 0]);
    }
    Point::affine(
        [v["x"][0].as_u64().unwrap(), v["x"][1].as_u64().unwrap()],
        [v["y"][0].as_u64().unwrap(), v["y"][1].as_u64().unwrap()],
    )
}

fn point_tuple(point: Point) -> String {
    format!(
        "{},{},{},{},{}",
        u64::from(point.infinity),
        point.x[0],
        point.x[1],
        point.y[0],
        point.y[1]
    )
}

fn velu(src: &Curve, torsion: Point, point: Point, ops: &mut GroupOps) -> Point {
    if point.infinity || point == torsion || point == src.neg(torsion) {
        return src.identity();
    }
    let plus = src.add(ops, point, torsion);
    let minus = src.add(ops, point, src.neg(torsion));
    assert!(!plus.infinity && !minus.infinity);
    let f = &src.f;
    Point::affine(
        f.sub(
            f.add(f.add(point.x, plus.x), minus.x),
            f.scale(torsion.x, 2),
        ),
        f.add(f.add(point.y, plus.y), minus.y),
    )
}

fn freeze_targets(parent_path: &str, public_path: &str, checker_path: &str) {
    assert!(!Path::new(public_path).exists() && !Path::new(checker_path).exists());
    let parent_bytes = fs::read(parent_path).unwrap();
    let mut public: Value = serde_json::from_slice(&parent_bytes).unwrap();
    assert_eq!(public["p"], 211);
    assert_eq!(public["r"], 229);
    assert_eq!(
        public["target_order_domain"],
        "different-j-order-batch-target-v1:13"
    );
    assert_eq!(public["targets"].as_array().unwrap().len(), 64);
    let f = Fp2::new(211).unwrap();
    let source = Curve {
        f,
        a: [190, 0],
        b: [37, 0],
    };
    let destination = Curve {
        f,
        a: [130, 0],
        b: [179, 0],
    };
    let g = parse_point(&public["generator"]);
    let mapped_g = parse_point(&public["mapped_generator"]);
    let t = Point::affine([3, 0], [1, 0]);
    let mut ops = GroupOps::default();
    assert!(source.is_on_curve(g) && destination.is_on_curve(mapped_g));
    assert_eq!(velu(&source, t, g, &mut ops), mapped_g);
    let mut entries = Vec::new();
    for scalar in 1..229u64 {
        let point = source.mul(&mut ops, g, scalar);
        let mapped = destination.mul(&mut ops, mapped_g, scalar);
        assert_eq!(velu(&source, t, point, &mut ops), mapped);
        let tuple = point_tuple(point);
        let digest = sha256_hex(format!("different-j-order-batch-target-v1:13:{tuple}").as_bytes());
        entries.push((digest, tuple, point, mapped, scalar));
    }
    entries.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
    for (existing, generated) in public["targets"]
        .as_array()
        .unwrap()
        .iter()
        .zip(entries.iter())
    {
        assert_eq!(parse_point(&existing["source"]), generated.2);
        assert_eq!(parse_point(&existing["destination"]), generated.3);
    }
    let held = &entries[64..68];
    let targets: Vec<Value> = held
        .iter()
        .map(|entry| json!({"source":entry.2,"destination":entry.3}))
        .collect();
    let scalars: Vec<u64> = held.iter().map(|entry| entry.4).collect();
    let parent_sha = sha256_hex(&parent_bytes);
    public["targets"] = json!(targets);
    public["target_count"] = json!(4);
    public["schema"] = json!("different_j_chi_heldout_public/v1");
    public["parent_public_sha256"] = json!(parent_sha);
    public["selection_positions_zero_indexed"] = json!([64, 65, 66, 67]);
    assert!(public.get("scalars").is_none() && public.get("answer_scalar").is_none());
    fs::write(
        public_path,
        format!("{}\n", serde_json::to_string_pretty(&public).unwrap()),
    )
    .unwrap();
    let public_sha = sha256_hex(&fs::read(public_path).unwrap());
    let checker = json!({"schema":"different_j_chi_heldout_checker/v1",
        "public_sha256":public_sha,"parent_public_sha256":parent_sha,
        "target_count":4,"scalars":scalars});
    fs::write(
        checker_path,
        format!("{}\n", serde_json::to_string_pretty(&checker).unwrap()),
    )
    .unwrap();
    println!(
        "{}",
        json!({"schema":"different_j_chi_heldout_freeze/v1",
        "public_sha256":public_sha,"checker_sha256":sha256_hex(&fs::read(checker_path).unwrap()),
        "parent_public_sha256":parent_sha,"selection_positions_zero_indexed":[64,65,66,67],
        "target_count":4,"source_image_replays":228,"group_ops":ops})
    );
}

fn direct_pairs(
    curve: &Curve,
    signed: &[Point],
    target: Point,
    ops: &mut GroupOps,
) -> BTreeSet<(usize, usize)> {
    let mut out = BTreeSet::new();
    for (i, &a) in signed.iter().enumerate() {
        for (j, &b) in signed.iter().enumerate().skip(i) {
            if curve.add(ops, a, b) == target {
                out.insert((i, j));
            }
        }
    }
    out
}

fn modmul(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % 229) as u64
}

fn insert_relation(pivots: &mut [Option<Vec<u64>>; 9], mut row: Vec<u64>) -> bool {
    for col in 0..9 {
        if row[col] == 0 {
            continue;
        }
        if let Some(pivot) = &pivots[col] {
            let factor = row[col];
            for j in col..10 {
                row[j] = (row[j] + 229 - modmul(factor, pivot[j])) % 229;
            }
        } else {
            let inverse = crypto_lib::cryptanalysis::residual_walk::pow_mod(row[col], 227, 229);
            for item in row.iter_mut().skip(col) {
                *item = modmul(*item, inverse);
            }
            pivots[col] = Some(row);
            return true;
        }
    }
    assert_eq!(row[9], 0);
    false
}

fn recover(pivots: &[Option<Vec<u64>>; 9]) -> Option<[u64; 9]> {
    if pivots.iter().any(Option::is_none) {
        return None;
    }
    let mut output = [0u64; 9];
    for col in (0..9).rev() {
        let row = pivots[col].as_ref().unwrap();
        let mut rhs = row[9];
        for j in col + 1..9 {
            rhs = (rhs + 229 - modmul(row[j], output[j])) % 229;
        }
        output[col] = rhs;
    }
    Some(output)
}

fn thread_cpu_ns() -> u64 {
    let mut ts = libc::timespec { tv_sec: 0, tv_nsec: 0 };
    assert_eq!(unsafe { libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut ts) }, 0);
    ts.tv_sec as u64 * 1_000_000_000 + ts.tv_nsec as u64
}

fn peak_rss_bytes() -> u64 {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::zeroed();
    assert_eq!(unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) }, 0);
    let raw = unsafe { usage.assume_init() }.ru_maxrss as u64;
    if cfg!(target_os = "macos") { raw } else { raw * 1024 }
}

fn run_heldout(path: &str, round: u64, arm: &str, position: usize) {
    assert!(round < 3 && position < 4);
    assert!(["chi-resultant", "direct-pair", "source-s3", "quotient-s3",
        "source-bsgs", "rho-negation", "rho-negation-aa"].contains(&arm));
    let bytes = fs::read(path).unwrap();
    let public: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(public["schema"], "different_j_chi_heldout_public/v1");
    assert_eq!((public["p"].as_u64(), public["r"].as_u64()), (Some(211), Some(229)));
    assert!(public.get("scalars").is_none() && public.get("answer_scalar").is_none());
    assert_eq!(public["selection_positions_zero_indexed"], json!([64, 65, 66, 67]));
    let f = Fp2::new(211).unwrap();
    let source = Curve { f, a: [190, 0], b: [37, 0] };
    let destination = Curve { f, a: [130, 0], b: [179, 0] };
    let g = parse_point(&public["generator"]);
    let mapped_g = parse_point(&public["mapped_generator"]);
    let t = Point::affine([3, 0], [1, 0]);
    let q = parse_point(&public["targets"][position]["source"]);
    let mapped_q = parse_point(&public["targets"][position]["destination"]);
    let wall = Instant::now();
    let cpu = thread_cpu_ns();
    take_field_counters();
    let mut setup_ops = GroupOps::default();
    assert!(source.is_on_curve(g) && source.is_on_curve(q));
    assert!(destination.is_on_curve(mapped_g) && destination.is_on_curve(mapped_q));
    assert_eq!(velu(&source, t, g, &mut setup_ops), mapped_g);
    assert_eq!(velu(&source, t, q, &mut setup_ops), mapped_q);
    assert!(source.mul(&mut setup_ops, g, 229).infinity);
    assert!(source.mul(&mut setup_ops, q, 229).infinity);
    let relation_arm = ["chi-resultant","direct-pair","source-s3","quotient-s3"].contains(&arm);
    let mut signed = Vec::new();
    let mut mapped_signed=Vec::new();
    if relation_arm {
        for witness in public["full_base_witness_source"].as_array().unwrap() {
            let projected = source.mul(&mut setup_ops, parse_point(witness), 195);
            assert!(!projected.infinity && source.mul(&mut setup_ops, projected, 229).infinity);
            signed.push(projected);
            signed.push(source.neg(projected));
        }
        assert_eq!(signed.len(), 16);
        mapped_signed = signed.iter().map(|&p| velu(&source, t, p, &mut setup_ops)).collect();
        let expected_destination = public["full_base_projected"].as_array().unwrap();
        for index in 0..8 {
            assert_eq!(mapped_signed[2 * index], parse_point(&expected_destination[index]));
        }
    }
    let base = if arm == "chi-resultant" { chi_base(&source, &signed) } else { Vec::new() };
    let mut table = BTreeMap::<u64, Vec<(usize, usize)>>::new();
    if arm == "direct-pair" {
        for i in 0..16 {
            for j in i..16 {
                let sum = source.add(&mut setup_ops, signed[i], signed[j]);
                if !sum.infinity {
                    table.entry(source.key(&sum)).or_default().push((i, j));
                }
            }
        }
    }
    let mut baby = BTreeMap::<u64, u64>::new();
    let mut giant_step = source.identity();
    if arm == "source-bsgs" {
        let mut cursor = source.identity();
        for small in 0..16 {
            baby.insert(source.key(&cursor), small);
            cursor = source.add(&mut setup_ops, cursor, g);
        }
        giant_step = source.neg(source.mul(&mut setup_ops, g, 16));
    }
    let setup_cpu_ns = thread_cpu_ns() - cpu;
    let setup_fp = take_field_counters();
    let seed = (0x2026_1015u64 ^ 211)
        .wrapping_add(round.wrapping_mul(0x9e37_79b9))
        .wrapping_add(((position + 64) as u64).wrapping_mul(0x0123_4567));
    let query_cpu_start = thread_cpu_ns();
    let mut ops = GroupOps::default();
    let mut stats = Stats::default();
    let mut s3_stats=S3Stats::default();
    let mut witnesses = Vec::new();
    let mut trials = 0u64;
    let mut relations = 0u64;
    let mut rank = 0u64;
    let mut recovered = None;
    let mut rho_attempts = Vec::new();
    if relation_arm {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut pivots: [Option<Vec<u64>>; 9] = std::array::from_fn(|_| None);
        let quotient=arm=="quotient-s3";
        let query_curve=if quotient {&destination} else {&source};
        let query_g=if quotient {mapped_g} else {g};
        let query_q=if quotient {mapped_q} else {q};
        for _ in 0..1832 {
            trials += 1;
            let a = rng.gen_range(0..229u64);
            let b = rng.gen_range(0..229u64);
            let a_g = query_curve.mul(&mut ops, query_g, a);
            let b_q = query_curve.mul(&mut ops, query_q, b);
            let candidate = query_curve.add(&mut ops, a_g, b_q);
            let projected = query_curve.mul(&mut ops, candidate, 195);
            if projected.infinity { continue; }
            let pairs = if arm == "chi-resultant" {
                polynomial_pairs(&source, &signed, &base, projected, &mut ops, &mut stats)
            } else if arm=="source-s3" {
                s3_pairs(&source,&signed,projected,&mut ops,&mut s3_stats)
            } else if arm=="quotient-s3" {
                s3_pairs(&destination,&mapped_signed,projected,&mut ops,&mut s3_stats)
            } else {
                let mut out = BTreeSet::new();
                if let Some(entries) = table.get(&source.key(&projected)) {
                    for &(i,j) in entries {
                        if source.add(&mut ops, signed[i], signed[j]) == projected {
                            out.insert((i,j));
                        }
                    }
                }
                out
            };
            if let Some(&(i, j)) = pairs.first() {
                let mut row = vec![0u64; 10];
                for index in [i,j] {
                    let column = index / 2;
                    let coefficient = if index % 2 == 0 { 1 } else { 228 };
                    row[column] = (row[column] + coefficient) % 229;
                }
                row[8] = (229 - modmul(195, b)) % 229;
                row[9] = modmul(195, a);
                witnesses.push(json!({"a":a,"b":b,"i":i,"j":j,"projected":projected}));
                relations += 1;
                if insert_relation(&mut pivots, row) { rank += 1; }
                if rank == 9 { break; }
            }
        }
        recovered = recover(&pivots).map(|sol| sol[8]);
    } else if arm == "source-bsgs" {
        let mut cursor = q;
        for giant in 0..15u64 {
            trials += 1;
            if let Some(&small) = baby.get(&source.key(&cursor)) {
                recovered = Some((giant * 16 + small) % 229);
                break;
            }
            cursor = source.add(&mut ops, cursor, giant_step);
        }
    } else {
        for retry in 0..32u64 {
            let try_seed = seed.wrapping_add(retry * 0x0100_0000);
            let result = rho_reference_negation(&source, g, q, 229, try_seed, 4580);
            ops.merge(result.group_ops);
            rho_attempts.push(json!({"retry":retry,"seed":try_seed,"steps":result.steps,
                "walks":result.walks,"group_ops":result.group_ops,"recovered":result.recovered,
                "verified":result.verified}));
            if result.verified { recovered = result.recovered; break; }
        }
    }
    let query_cpu_ns = thread_cpu_ns() - query_cpu_start;
    let query_fp = take_field_counters();
    let mut verify_ops = GroupOps::default();
    let verified = recovered.is_some_and(|scalar| {
        source.mul(&mut verify_ops, g, scalar) == q
            && destination.mul(&mut verify_ops, mapped_g, scalar) == mapped_q
    });
    let verify_fp = take_field_counters();
    let peak_rss = peak_rss_bytes();
    let result = json!({"schema":"different_j_chi_heldout_run/v3",
        "public_sha256":sha256_hex(&bytes),"round":round,"arm":arm,"position":position,
        "seed":seed,"target":q,"mapped_target":mapped_q,
        "setup_group_ops":setup_ops,"query_group_ops":ops,"verify_group_ops":verify_ops,
        "setup_fp_muls":setup_fp.0,"setup_fp_invs":setup_fp.1,
        "query_fp_muls":query_fp.0,"query_fp_invs":query_fp.1,
        "verify_fp_muls":verify_fp.0,"verify_fp_invs":verify_fp.1,
        "setup_cpu_ns":setup_cpu_ns,"query_cpu_ns":query_cpu_ns,
        "total_cpu_ns":thread_cpu_ns()-cpu,"total_wall_ns":wall.elapsed().as_nanos(),
        "peak_rss_bytes":peak_rss,"trials":trials,"relations":relations,
        "rank":rank,"recovered":recovered,"verified":verified,
        "relation_witnesses":witnesses,"rho_attempts":rho_attempts,
        "chi_stats":{"resultants":stats.resultants,"gcds":stats.gcds,
          "split_attempts":stats.split_attempts,
          "exhaustive_fallback_evaluations":stats.exhaustive_fallback_evaluations,
          "root_candidates":stats.root_candidates,"chi_lift_scans":stats.chi_lift_scans,
          "exact_checks":stats.exact_checks,"false_candidates":stats.false_candidates,
          "vertical_cases":stats.vertical_cases,"tangent_cases":stats.tangent_cases,
          "horizontal_cases":stats.horizontal_cases,"infinity_cases":stats.infinity_cases,
          "exceptional_chi_cases":stats.exceptional_chi_cases,
          "multiple_resultant_roots":stats.multiple_resultant_roots},
        "s3_stats":{"quadratics":s3_stats.quadratics,"roots":s3_stats.roots,
          "lift_scans":s3_stats.lift_scans,"exact_checks":s3_stats.exact_checks,
          "false_candidates":s3_stats.false_candidates,
          "degenerate":s3_stats.degenerate,"sqrt_failures":s3_stats.sqrt_failures,
          "infinity_cases":s3_stats.infinity_cases}});
    println!("{result}");
    assert!(verified && peak_rss <= 512 * 1024 * 1024);
}

fn exhaustive_s3(path: &str) {
    let bytes=fs::read(path).unwrap();
    let public:Value=serde_json::from_slice(&bytes).unwrap();
    assert_eq!(public["p"],211);
    assert_eq!(public["r"],229);
    let f=Fp2::new(211).unwrap();
    let source=Curve{f,a:[190,0],b:[37,0]};
    let destination=Curve{f,a:[130,0],b:[179,0]};
    let g=parse_point(&public["generator"]);
    let mg=parse_point(&public["mapped_generator"]);
    let t=Point::affine([3,0],[1,0]);
    let mut ops=GroupOps::default();
    let mut signed=Vec::new();
    for witness in public["full_base_witness_source"].as_array().unwrap() {
        let projected=source.mul(&mut ops,parse_point(witness),195);
        signed.push(projected);
        signed.push(source.neg(projected));
    }
    let mapped:Vec<Point>=signed.iter().map(|&q|velu(&source,t,q,&mut ops)).collect();
    take_field_counters();
    let start=Instant::now();
    let mut source_stats=S3Stats::default();
    let mut quotient_stats=S3Stats::default();
    let mut zero_source=0u64;
    let mut zero_quotient=0u64;
    let mut pairs_source=0u64;
    let mut pairs_quotient=0u64;
    let mut sq=source.identity();
    let mut dq=destination.identity();
    for scalar_label in 0..229 {
        let expected_source=direct_pairs(&source,&signed,sq,&mut ops);
        let actual_source=s3_pairs(&source,&signed,sq,&mut ops,&mut source_stats);
        assert_eq!(actual_source,expected_source);
        let expected_quotient=direct_pairs(&destination,&mapped,dq,&mut ops);
        let actual_quotient=s3_pairs(&destination,&mapped,dq,&mut ops,&mut quotient_stats);
        assert_eq!(actual_quotient,expected_quotient);
        assert_eq!(actual_source,actual_quotient);
        println!("{}",json!({"schema":"different_j_chi_s3_case/v2",
            "scalar_label":scalar_label,"source_target":sq,
            "quotient_target":dq,
            "source_pairs":actual_source.iter().copied().collect::<Vec<_>>(),
            "quotient_pairs":actual_quotient.iter().copied().collect::<Vec<_>>(),
            "status":"pass"}));
        zero_source+=u64::from(actual_source.is_empty());
        zero_quotient+=u64::from(actual_quotient.is_empty());
        pairs_source+=actual_source.len() as u64;
        pairs_quotient+=actual_quotient.len() as u64;
        sq=source.add(&mut ops,sq,g);
        dq=destination.add(&mut ops,dq,mg);
    }
    assert!(sq.infinity && dq.infinity);
    let fp=take_field_counters();
    println!("{}",json!({"schema":"different_j_chi_s3_exhaustive/v2",
        "public_sha256":sha256_hex(&bytes),"source_cases":229,
        "quotient_cases":229,"source_pairs":pairs_source,
        "quotient_pairs":pairs_quotient,"source_zero_cases":zero_source,
        "quotient_zero_cases":zero_quotient,
        "source_s3":{"quadratics":source_stats.quadratics,
          "roots":source_stats.roots,"lift_scans":source_stats.lift_scans,
          "exact_checks":source_stats.exact_checks,
          "false_candidates":source_stats.false_candidates,
          "degenerate":source_stats.degenerate,
          "sqrt_failures":source_stats.sqrt_failures,
          "infinity_cases":source_stats.infinity_cases},
        "quotient_s3":{"quadratics":quotient_stats.quadratics,
          "roots":quotient_stats.roots,"lift_scans":quotient_stats.lift_scans,
          "exact_checks":quotient_stats.exact_checks,
          "false_candidates":quotient_stats.false_candidates,
          "degenerate":quotient_stats.degenerate,
          "sqrt_failures":quotient_stats.sqrt_failures,
          "infinity_cases":quotient_stats.infinity_cases},
        "group_ops":ops,"fp_muls":fp.0,"fp_invs":fp.1,
        "wall_ns":start.elapsed().as_nanos(),"status":"pass"}));
}

fn polynomial_pairs(
    curve: &Curve,
    signed: &[Point],
    base: &Poly,
    target: Point,
    ops: &mut GroupOps,
    stats: &mut Stats,
) -> BTreeSet<(usize, usize)> {
    let f = &curve.f;
    let mut out = BTreeSet::new();
    for (i, &point) in signed.iter().enumerate() {
        if target.infinity {
            stats.infinity_cases += 1;
            for (j, &candidate) in signed.iter().enumerate() {
                if candidate == curve.neg(point) {
                    out.insert((i.min(j), i.max(j)));
                }
            }
            continue;
        }
        if point == target {
            stats.vertical_cases += 1; // The second summand would be infinity.
            continue;
        }
        let minus_target = curve.neg(target);
        let slope = if point == minus_target {
            stats.tangent_cases += 1;
            let numerator = f.add(f.scale(f.sqr(point.x), 3), curve.a);
            let Some(denominator) = f.inv(f.scale(point.y, 2)) else {
                stats.vertical_cases += 1;
                continue;
            };
            f.mul(numerator, denominator)
        } else {
            let Some(denominator) = f.inv(f.sub(minus_target.x, point.x)) else {
                stats.vertical_cases += 1;
                continue;
            };
            f.mul(f.sub(minus_target.y, point.y), denominator)
        };
        if slope == f.zero() {
            stats.horizontal_cases += 1;
            let z = chi(curve, point);
            for (j, &candidate) in signed.iter().enumerate() {
                stats.chi_lift_scans += 1;
                if chi(curve, candidate) == z {
                    stats.exact_checks += 1;
                    if curve.add(ops, point, candidate) == target {
                        out.insert((i.min(j), i.max(j)));
                    } else {
                        stats.false_candidates += 1;
                    }
                }
            }
            continue;
        }
        let resultant = chord_resultant(curve, point, slope);
        stats.resultants += 1;
        let intersection = gcd(f, resultant.clone(), base.clone());
        stats.gcds += 1;
        for root in all_roots(f, &intersection, stats) {
            stats.root_candidates += 1;
            if root_multiplicity(f, &resultant, root) > 1 {
                stats.multiple_resultant_roots += 1;
            }
            for (j, &candidate) in signed.iter().enumerate() {
                stats.chi_lift_scans += 1;
                if chi(curve, candidate) == Some(root) {
                    stats.exact_checks += 1;
                    if curve.add(ops, point, candidate) == target {
                        out.insert((i.min(j), i.max(j)));
                    } else {
                        stats.false_candidates += 1;
                    }
                }
            }
        }
        for (j, &candidate) in signed.iter().enumerate() {
            if chi(curve, candidate).is_none() {
                stats.exceptional_chi_cases += 1;
                stats.exact_checks += 1;
                if curve.add(ops, point, candidate) == target {
                    out.insert((i.min(j), i.max(j)));
                }
            }
        }
    }
    out
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("freeze-targets") {
        assert_eq!(args.len(), 5);
        freeze_targets(&args[2], &args[3], &args[4]);
        return;
    }
    if args.get(1).map(String::as_str) == Some("run-heldout") {
        assert_eq!(args.len(), 6);
        run_heldout(&args[2], args[3].parse().unwrap(), &args[4], args[5].parse().unwrap());
        return;
    }
    if args.get(1).map(String::as_str) == Some("exhaustive-s3") {
        assert_eq!(args.len(),3);
        exhaustive_s3(&args[2]);
        return;
    }
    assert_eq!(args.get(1).map(String::as_str), Some("exhaustive"));
    let path = std::env::args().nth(2).expect("public input path");
    let bytes = fs::read(&path).unwrap();
    let public: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(public["p"], 211);
    assert_eq!(public["r"], 229);
    assert_eq!(public["cofactor"], 195);
    let f = Fp2::new(211).unwrap();
    let curve = Curve {
        f,
        a: [190, 0],
        b: [37, 0],
    };
    let generator = parse_point(&public["generator"]);
    let mut ops = GroupOps::default();
    assert!(curve.is_on_curve(generator));
    assert!(curve.mul(&mut ops, generator, 229).infinity);
    let mut signed = Vec::new();
    for witness in public["full_base_witness_source"].as_array().unwrap() {
        let raw = parse_point(witness);
        let projected = curve.mul(&mut ops, raw, 195);
        assert!(!projected.infinity && curve.mul(&mut ops, projected, 229).infinity);
        signed.push(projected);
        signed.push(curve.neg(projected));
    }
    assert_eq!(signed.len(), 16);
    let base = chi_base(&curve, &signed);
    take_field_counters();
    let start = Instant::now();
    let mut stats = Stats::default();
    let mut target = curve.identity();
    let mut checked = 0;
    let mut pair_count = 0;
    let mut zero_cases = 0;
    for scalar in 0..229 {
        let expected = direct_pairs(&curve, &signed, target, &mut ops);
        let actual = polynomial_pairs(&curve, &signed, &base, target, &mut ops, &mut stats);
        assert_eq!(
            actual, expected,
            "chi resultant mismatch at scalar {scalar}"
        );
        pair_count += actual.len();
        zero_cases += usize::from(actual.is_empty());
        checked += 1;
        println!(
            "{}",
            json!({"schema":"different_j_chi_resultant_case/v1","scalar_label":scalar,
            "target":target,"pairs":actual.iter().copied().collect::<Vec<_>>(),"status":"pass"})
        );
        target = curve.add(&mut ops, target, generator);
    }
    assert!(target.infinity);
    let torsion = Point::affine(f.from_fp(3), f.one());
    for (label, exceptional) in [
        ("torsion_plus", torsion),
        ("torsion_minus", curve.neg(torsion)),
    ] {
        let expected = direct_pairs(&curve, &signed, exceptional, &mut ops);
        let actual = polynomial_pairs(&curve, &signed, &base, exceptional, &mut ops, &mut stats);
        assert_eq!(actual, expected, "chi resultant mismatch at {label}");
        pair_count += actual.len();
        zero_cases += usize::from(actual.is_empty());
        checked += 1;
        println!(
            "{}",
            json!({"schema":"different_j_chi_resultant_case/v1","scalar_label":label,
            "target":exceptional,"pairs":actual.iter().copied().collect::<Vec<_>>(),"status":"pass"})
        );
    }
    let field = take_field_counters();
    println!(
        "{}",
        json!({"schema":"different_j_chi_resultant_summary/v1",
        "public_sha256":sha256_hex(&bytes),"p":211,"r":229,"factor_base_columns":8,
        "signed_points":signed.len(),"chi_polynomial_degree":base.len()-1,
        "cases":checked,"pair_count":pair_count,"zero_cases":zero_cases,
        "resultants":stats.resultants,"gcds":stats.gcds,"split_attempts":stats.split_attempts,
        "exhaustive_fallback_evaluations":stats.exhaustive_fallback_evaluations,
        "root_candidates":stats.root_candidates,"chi_lift_scans":stats.chi_lift_scans,
        "exact_checks":stats.exact_checks,"false_candidates":stats.false_candidates,
        "vertical_cases":stats.vertical_cases,"tangent_cases":stats.tangent_cases,
        "horizontal_cases":stats.horizontal_cases,"infinity_cases":stats.infinity_cases,
        "exceptional_chi_cases":stats.exceptional_chi_cases,
        "multiple_resultant_roots":stats.multiple_resultant_roots,
        "group_ops":ops,"fp_muls":field.0,"fp_invs":field.1,
        "wall_ns":start.elapsed().as_nanos(),"status":"pass"})
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chord_elimination_vanishes_on_group_sum() {
        let f = Fp2::new(211).unwrap();
        let curve = Curve {
            f,
            a: [190, 0],
            b: [37, 0],
        };
        let generator = Point::affine([170, 0], [0, 109]);
        let mut ops = GroupOps::default();
        let p1 = curve.mul(&mut ops, generator, 17);
        let p2 = curve.mul(&mut ops, generator, 73);
        let target = curve.add(&mut ops, p1, p2);
        let minus_target = curve.neg(target);
        let slope = f.mul(
            f.sub(minus_target.y, p1.y),
            f.inv(f.sub(minus_target.x, p1.x)).unwrap(),
        );
        assert_ne!(slope, f.zero());
        let z = chi(&curve, p2).unwrap();
        let h = chord_resultant(&curve, p1, slope);
        assert_eq!(eval(&f, &h, z), f.zero());
    }

    #[test]
    fn rational_root_split_returns_every_distinct_root() {
        let f = Fp2::new(211).unwrap();
        let roots = [[5, 1], [17, 3], [201, 40]];
        let polynomial = roots.iter().fold(vec![f.one()], |acc, &root| {
            mul(&f, &acc, &vec![f.neg(root), f.one()])
        });
        let mut stats = Stats::default();
        let actual = all_roots(&f, &polynomial, &mut stats);
        let mut expected = roots.to_vec();
        expected.sort_unstable();
        assert_eq!(actual, expected);
    }
}
