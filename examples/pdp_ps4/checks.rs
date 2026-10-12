//! Correctness checks of the field, the curve constants (`a = 0`), the
//! translation by `T`, and the exactness theorem of `ps1_algebra/README.md`
//! §3 re-derived and re-checked for `K_0`.

use crate::curve::{add, lift, neg, on_curve, sorted3, thirds, Pt};
use crate::field::{itoh_tsujii_products, Fld, Onb, Small};
use crate::targets::{enumerate_weight, Rng};
use crate::{arg, arg_s, hex};
use crypto_lib::cryptanalysis::koblitz_rotated_chain::{OnbCurve, OnbPoint};
use std::io::Write;

fn open_out(args: &[String]) -> Box<dyn Write> {
    match arg_s(args, "--out") {
        Some(p) => Box::new(
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(p)
                .unwrap(),
        ),
        None => Box::new(std::io::stdout()),
    }
}

const ONB_DEGREES: [usize; 7] = [11, 23, 29, 33, 35, 39, 41];

/// The field against the library's `Type2Onb`.
pub fn cmd_check_field(args: &[String]) {
    let mut out = open_out(args);
    let samples: usize = arg(args, "--samples").unwrap_or(200_000);
    let mut degrees = vec![3usize, 5, 9];
    degrees.extend(ONB_DEGREES);
    for n in degrees {
        let f = Onb::new(n).unwrap();
        let lib = f.lib;
        let mut rng = Rng(0xC0FFEE ^ n as u64);
        let mask = f.mask;
        let (mut mul_pairs, mut mul_bad) = (0u64, 0u64);
        if n <= 11 {
            for a in 0..=mask {
                for b in 0..=mask {
                    mul_pairs += 1;
                    if f.mul_raw(a, b) as u128 != lib.mul(a as u128, b as u128) {
                        mul_bad += 1;
                    }
                }
            }
        } else {
            for _ in 0..samples {
                let (a, b) = (rng.next() & mask, rng.next() & mask);
                mul_pairs += 1;
                if f.mul_raw(a, b) as u128 != lib.mul(a as u128, b as u128) {
                    mul_bad += 1;
                }
            }
        }
        let (mut sqr_n, mut sqr_bad, mut inv_n, mut inv_bad) = (0u64, 0u64, 0u64, 0u64);
        let inv_samples = if n <= 11 { 1usize << n } else { samples / 10 };
        for k in 0..inv_samples {
            let a = if n <= 11 { k as u64 } else { rng.next() & mask };
            sqr_n += 1;
            if f.sqr_raw(a) as u128 != lib.square(a as u128) {
                sqr_bad += 1;
            }
            if a != 0 {
                inv_n += 1;
                let i = f.inv_c(a);
                if Some(i as u128) != lib.inverse(a as u128) || f.mul_raw(i, a) != f.mask {
                    inv_bad += 1;
                }
            }
        }
        // Counted products of one inversion = the Itoh-Tsujii count.
        let before = f.snapshot();
        let _ = f.inv_c(1 + (rng.next() & (mask >> 1)));
        let d = f.snapshot().sub(&before);
        // Lifts against the library's points_with_x.
        let curve = OnbCurve::new(lib);
        let (mut lift_n, mut lift_bad) = (0u64, 0u64);
        let lift_samples = if n <= 11 { 1usize << n } else { 20_000 };
        for k in 0..lift_samples {
            let x = if n <= 11 { k as u64 } else { rng.next() & mask };
            lift_n += 1;
            let pts = curve.points_with_x(x as u128);
            let mine = lift(&f, x);
            let ok = match (mine, pts.first()) {
                (None, None) => true,
                (Some(y), Some(p)) => p.y as u64 == y && on_curve(&f, x, y),
                _ => false,
            };
            if !ok {
                lift_bad += 1;
            }
        }
        writeln!(
            out,
            "{{\"kind\":\"field_check\",\"n\":{n},\"exhaustive\":{},\"mul_pairs\":{mul_pairs},\"mul_mismatch\":{mul_bad},\"sqr_n\":{sqr_n},\"sqr_mismatch\":{sqr_bad},\"inv_n\":{inv_n},\"inv_mismatch\":{inv_bad},\"inv_products_counted\":{},\"inv_products_formula\":{},\"inv_squarings_counted\":{},\"lift_n\":{lift_n},\"lift_mismatch\":{lift_bad}}}",
            n <= 11,
            d.mul,
            itoh_tsujii_products(n),
            d.sqr + d.sqr_k
        )
        .unwrap();
        out.flush().unwrap();
    }
}

fn koblitz_count(n: usize) -> i128 {
    // K_0 over F_2 has 4 points, t_1 = -1; t_n = t_1 t_{n-1} - 2 t_{n-2}, t_0 = 2.
    let (mut t0, mut t1) = (2i128, -1i128);
    for _ in 2..=n {
        let t2 = -t1 - 2 * t0;
        t0 = t1;
        t1 = t2;
    }
    (1i128 << n) + 1 - t1
}

/// The group-law constants for `a = 0` against the library's `OnbCurve`, the
/// point count, `S₃`, and the translation by `T = (0, 1)`.
pub fn cmd_check_curve(args: &[String]) {
    let mut out = open_out(args);
    let samples: usize = arg(args, "--samples").unwrap_or(200_000);
    let mut degrees = vec![3usize, 5, 9];
    degrees.extend(ONB_DEGREES);
    for n in degrees {
        let f = Onb::new(n).unwrap();
        let curve = OnbCurve::new(f.lib);
        let mask = f.mask;
        let mut rng = Rng(0xFACADE ^ n as u64);
        let t_pt = Pt::Aff(0, f.mask);
        // --- the point T ----------------------------------------------------
        let t_on = on_curve(&f, 0, f.mask);
        let t_twice = add(&f, t_pt, t_pt) == Pt::Inf;
        // --- point count (exhaustive up to n = 23) --------------------------
        let (count, expected) = if n <= 23 {
            let mut c = 2u64; // O and T
            for x in 1..=mask {
                if lift(&f, x).is_some() {
                    c += 2;
                }
            }
            (Some(c as i128), Some(koblitz_count(n)))
        } else {
            (None, None)
        };
        // --- translation by T ------------------------------------------------
        let exhaustive_t = n <= 23;
        let t_samples = if exhaustive_t { mask as usize } else { samples };
        let (mut t_n, mut t_bad, mut t_neg_bad) = (0u64, 0u64, 0u64);
        for k in 0..t_samples {
            let x = if exhaustive_t {
                k as u64 + 1
            } else {
                rng.next() & mask
            };
            if x == 0 {
                continue;
            }
            let Some(y) = lift(&f, x) else { continue };
            let xinv = f.inv_c(x);
            t_n += 1;
            match add(&f, Pt::Aff(x, y), t_pt) {
                Pt::Aff(x2, _) if x2 == xinv => {}
                _ => t_bad += 1,
            }
            match add(&f, neg(Pt::Aff(x, y)), t_pt) {
                Pt::Aff(x2, _) if x2 == xinv => {}
                _ => t_neg_bad += 1,
            }
        }
        // --- the group law against the library --------------------------------
        let (mut g_n, mut g_bad, mut s3_bad, mut s3_n, mut s3_nonzero_random) =
            (0u64, 0u64, 0u64, 0u64, 0u64);
        let pairs = if n <= 9 { 0 } else { samples };
        let mut pts: Vec<(u64, u64)> = Vec::new();
        if n <= 9 {
            for x in 0..=mask {
                if let Some(y) = lift(&f, x) {
                    pts.push((x, y));
                    if x != 0 {
                        pts.push((x, y ^ x));
                    }
                }
            }
        }
        let lib_pt = |p: Pt| match p {
            Pt::Inf => OnbPoint::INFINITY,
            Pt::Aff(x, y) => OnbPoint::affine(x as u128, y as u128),
        };
        let cmp = |p: Pt, q: Pt, g_n: &mut u64, g_bad: &mut u64| {
            let mine = add(&f, p, q);
            let theirs = curve.add(lib_pt(p), lib_pt(q));
            *g_n += 1;
            let same = match mine {
                Pt::Inf => theirs.infinity,
                Pt::Aff(x, y) => !theirs.infinity && theirs.x as u64 == x && theirs.y as u64 == y,
            };
            if !same {
                *g_bad += 1;
            }
        };
        if n <= 9 {
            for &(x1, y1) in &pts {
                for &(x2, y2) in &pts {
                    cmp(Pt::Aff(x1, y1), Pt::Aff(x2, y2), &mut g_n, &mut g_bad);
                }
                cmp(Pt::Aff(x1, y1), Pt::Inf, &mut g_n, &mut g_bad);
            }
        }
        for _ in 0..pairs {
            // random points, with repeated and negated operands mixed in
            let rand_pt = |rng: &mut Rng| loop {
                let x = rng.next() & mask;
                if let Some(y) = lift(&f, x) {
                    return if rng.next() & 1 == 1 && x != 0 {
                        (x, y ^ x)
                    } else {
                        (x, y)
                    };
                }
            };
            let p = rand_pt(&mut rng);
            let q = match rng.next() % 4 {
                0 => p,
                1 => (p.0, p.1 ^ p.0),
                _ => rand_pt(&mut rng),
            };
            cmp(Pt::Aff(p.0, p.1), Pt::Aff(q.0, q.1), &mut g_n, &mut g_bad);
            // S3(x1, x2, x(P1 ± P2)) = (x1x2 + x1x3 + x2x3)² + x1x2x3 + 1 must vanish.
            let s3 = |a: u64, b: u64, c: u64| {
                let s = f.mul_raw(a, b) ^ f.mul_raw(a, c) ^ f.mul_raw(b, c);
                f.sqr_raw(s) ^ f.mul_raw(f.mul_raw(a, b), c) ^ f.mask
            };
            for sign in 0..2 {
                let qq = if sign == 0 { q } else { (q.0, q.1 ^ q.0) };
                if let Pt::Aff(x3, _) = add(&f, Pt::Aff(p.0, p.1), Pt::Aff(qq.0, qq.1)) {
                    s3_n += 1;
                    if s3(p.0, q.0, x3) != 0 {
                        s3_bad += 1;
                    }
                    let rnd = rng.next() & mask;
                    if s3(p.0, q.0, rnd) != 0 {
                        s3_nonzero_random += 1;
                    }
                }
            }
        }
        writeln!(
            out,
            "{{\"kind\":\"curve_check\",\"n\":{n},\"T_on_curve\":{t_on},\"T_plus_T_is_O\":{t_twice},\"points\":{},\"points_expected\":{},\"T_translation_exhaustive\":{exhaustive_t},\"T_translation_points\":{t_n},\"T_translation_mismatch\":{t_bad},\"T_translation_neg_mismatch\":{t_neg_bad},\"group_law_pairs\":{g_n},\"group_law_mismatch_vs_lib\":{g_bad},\"S3_sums\":{s3_n},\"S3_nonzero_on_sums\":{s3_bad},\"S3_nonzero_on_random_third\":{s3_nonzero_random}}}",
            count.map_or("null".to_string(), |c| c.to_string()),
            expected.map_or("null".to_string(), |c| c.to_string()),
        )
        .unwrap();
        out.flush().unwrap();
    }
}

// ----------------------------------------------------------------------
// The exactness theorem, K_0.

struct ThmRow {
    n: usize,
    v_name: String,
    v_size: usize,
    targets: usize,
    degenerate_targets: usize,
    triples: u64,
    roots: u64,
    decomps: u64,
    decomps_not_roots: u64,
    extra_f1: u64,
    extra_f2: u64,
    extra_unexplained: u64,
    extras_at_nondegenerate: u64,
    targets_exact: usize,
    targets_with_extras: usize,
}

fn theorem_row(s: &Small, v: &[u64], v_name: &str, xrs: &[u64]) -> ThmRow {
    let f = &s.f;
    let n = f.n;
    let one = s.one();
    let size = s.size;
    let lifts: Vec<Option<u64>> = (0..size as u64).map(|x| lift(s, x)).collect();
    let mut in_v = vec![false; size];
    for &x in v {
        in_v[x as usize] = true;
    }
    // Pair table of the first quadratic: S3(x1, x2, t) = A t² + B t + C.
    let mut pairs: Vec<(u64, u64, u64, u64, u64)> = Vec::new();
    for (i, &x1) in v.iter().enumerate() {
        for &x2 in &v[i..] {
            let a = s.sqr(x1 ^ x2);
            let bb = s.mul(x1, x2);
            let c = s.sqr(bb) ^ one;
            pairs.push((x1, x2, a, bb, c));
        }
    }
    let mut row = ThmRow {
        n,
        v_name: v_name.to_string(),
        v_size: v.len(),
        targets: 0,
        degenerate_targets: 0,
        triples: 0,
        roots: 0,
        decomps: 0,
        decomps_not_roots: 0,
        extra_f1: 0,
        extra_f2: 0,
        extra_unexplained: 0,
        extras_at_nondegenerate: 0,
        targets_exact: 0,
        targets_with_extras: 0,
    };
    for &xr in xrs {
        let Some(yr) = lifts[xr as usize] else {
            continue;
        };
        if xr == 0 {
            continue;
        }
        row.targets += 1;
        let degenerate = in_v[xr as usize] || in_v[s.inv(xr) as usize];
        if degenerate {
            row.degenerate_targets += 1;
        }
        // D: rational decompositions, from the exact group law.
        let mut d: Vec<[u64; 3]> = Vec::new();
        let liftable: Vec<(u64, u64)> = v
            .iter()
            .filter_map(|&x| lifts[x as usize].map(|y| (x, y)))
            .collect();
        for (i, &p1) in liftable.iter().enumerate() {
            for &p2 in &liftable[i..] {
                for t in thirds(s, p1, p2, (xr, yr)).into_iter().flatten() {
                    if in_v[t as usize] {
                        d.push(sorted3(p1.0, p2.0, t));
                    }
                }
            }
        }
        d.sort_unstable();
        d.dedup();
        // A: the roots of S4(x1, x2, x3, x_R) over V^3, by the resultant of
        // the two quadratics S3(x1, x2, t) and S3(x3, x_R, t).
        let mut a_roots: Vec<[u64; 3]> = Vec::new();
        for &x3 in v {
            let a2 = s.sqr(x3 ^ xr);
            let b2 = s.mul(x3, xr);
            let c2 = s.sqr(b2) ^ one;
            for &(x1, x2, a1, b1, c1) in &pairs {
                row.triples += 1;
                let u = s.mul(a1, c2) ^ s.mul(a2, c1);
                let w = s.mul(a1, b2) ^ s.mul(a2, b1);
                let z = s.mul(b1, c2) ^ s.mul(b2, c1);
                if s.sqr(u) ^ s.mul(w, z) == 0 {
                    a_roots.push(sorted3(x1, x2, x3));
                }
            }
        }
        a_roots.sort_unstable();
        a_roots.dedup();
        row.roots += a_roots.len() as u64;
        row.decomps += d.len() as u64;
        row.decomps_not_roots += d
            .iter()
            .filter(|t| a_roots.binary_search(t).is_err())
            .count() as u64;
        let mut extras = 0u64;
        for t in &a_roots {
            if d.binary_search(t).is_ok() {
                continue;
            }
            extras += 1;
            let non: Vec<usize> = (0..3).filter(|&i| lifts[t[i] as usize].is_none()).collect();
            let mut explained = false;
            if non.len() == 2 {
                let (a, b) = (t[non[0]], t[non[1]]);
                let c = t[3 - non[0] - non[1]];
                if a == b && c == xr {
                    row.extra_f1 += 1;
                    explained = true;
                } else if s.mul(a, b) == one && s.mul(c, xr) == one {
                    row.extra_f2 += 1;
                    explained = true;
                }
            }
            if !explained {
                row.extra_unexplained += 1;
            }
            if !degenerate {
                row.extras_at_nondegenerate += 1;
            }
        }
        if extras == 0 {
            row.targets_exact += 1;
        } else {
            row.targets_with_extras += 1;
        }
    }
    row
}

/// `S₄` is symmetric in its four arguments: a spot check of the resultant form.
fn s4_symmetry_check(s: &Small, samples: usize) -> (u64, u64) {
    let one = s.one();
    let res = |x1: u64, x2: u64, x3: u64, x4: u64| {
        let a1 = s.sqr(x1 ^ x2);
        let b1 = s.mul(x1, x2);
        let c1 = s.sqr(b1) ^ one;
        let a2 = s.sqr(x3 ^ x4);
        let b2 = s.mul(x3, x4);
        let c2 = s.sqr(b2) ^ one;
        let u = s.mul(a1, c2) ^ s.mul(a2, c1);
        let w = s.mul(a1, b2) ^ s.mul(a2, b1);
        let z = s.mul(b1, c2) ^ s.mul(b2, c1);
        s.sqr(u) ^ s.mul(w, z)
    };
    let mut rng = Rng(77);
    let (mut n, mut bad) = (0u64, 0u64);
    for _ in 0..samples {
        let v: Vec<u64> = (0..4).map(|_| rng.next() & s.f.mask).collect();
        let base = res(v[0], v[1], v[2], v[3]);
        for perm in [[0, 2, 1, 3], [0, 3, 2, 1], [1, 2, 0, 3], [2, 3, 0, 1]] {
            n += 1;
            if res(v[perm[0]], v[perm[1]], v[perm[2]], v[perm[3]]) != base {
                bad += 1;
            }
        }
    }
    (n, bad)
}

pub fn cmd_check_theorem(args: &[String]) {
    let mut out = open_out(args);
    let extra_random: usize = arg(args, "--random-targets").unwrap_or(0);
    // (n, name, V)
    let mut cases: Vec<(usize, String, Vec<u64>)> = Vec::new();
    for n in [3usize, 5, 9] {
        cases.push((n, "F".to_string(), (0..(1u64 << n)).collect()));
    }
    for w in [2usize, 3] {
        let v: Vec<u64> = enumerate_weight(11, w)
            .into_iter()
            .map(|x| x as u64)
            .collect();
        cases.push((11, format!("wt<={w}"), v));
    }
    for (n, name, v) in cases {
        let s = Small::new(n).unwrap();
        let (sn, sbad) = s4_symmetry_check(&s, 20_000);
        let all: Vec<u64> = (1..(1u64 << n)).collect();
        let xrs: Vec<u64> = if extra_random > 0 && n > 9 {
            // Every degenerate target and a random sample of the others.
            let mut rng = Rng(5150);
            let mut set: Vec<u64> = all
                .iter()
                .copied()
                .filter(|&x| v.contains(&x) || v.contains(&s.inv(x)))
                .collect();
            for _ in 0..extra_random {
                set.push(1 + rng.next() % ((1u64 << n) - 1));
            }
            set.sort_unstable();
            set.dedup();
            set
        } else {
            all
        };
        let t0 = std::time::Instant::now();
        let r = theorem_row(&s, &v, &name, &xrs);
        writeln!(
            out,
            "{{\"kind\":\"theorem_check\",\"n\":{},\"V\":\"{}\",\"V_size\":{},\"targets_checked\":{},\"degenerate_targets\":{},\"triples_tested\":{},\"S4_roots_total\":{},\"rational_decompositions_total\":{},\"decompositions_that_are_not_roots\":{},\"extra_roots_F1\":{},\"extra_roots_F2\":{},\"extra_roots_unexplained\":{},\"extra_roots_at_nondegenerate_targets\":{},\"targets_where_roots_equal_decompositions\":{},\"targets_with_extra_roots\":{},\"s4_symmetry_checks\":{sn},\"s4_symmetry_failures\":{sbad},\"ms_unisolated\":{:.0}}}",
            r.n, r.v_name, r.v_size, r.targets, r.degenerate_targets, r.triples, r.roots,
            r.decomps, r.decomps_not_roots, r.extra_f1, r.extra_f2, r.extra_unexplained,
            r.extras_at_nondegenerate, r.targets_exact, r.targets_with_extras,
            t0.elapsed().as_secs_f64() * 1e3
        )
        .unwrap();
        out.flush().unwrap();
        let _ = hex(0);
    }
}
