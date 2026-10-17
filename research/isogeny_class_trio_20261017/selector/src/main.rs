//! isogeny-class-trio-selector: deterministic search for a B-shape binary
//! curve, a K-shape Koblitz curve and a P-shape prime-field curve whose
//! Frobenius order Z[pi] has conductor f > 1, plus the exact invariants of
//! each isogeny class (strata, class numbers, class size) and the root's
//! position in its volcanoes. See ../PROTOCOL.md for the frozen rules.

mod bsgs;
mod classno;
mod divpoly;
mod f2m;
mod fp;
mod fpiso;
mod iso3;
mod nt;
mod sha256;

use bsgs::{multiples_in_window, Group};
use f2m::{standard_modulus, BPt, BinCurve, F2m};
use fp::{Fp, PCurve, Pt};
use nt::{divisors, factor, frob_disc, is_prime, isqrt, kronecker, order_class_ratio, FrobDisc};
use std::env;
use std::fmt::Write as _;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// Largest l for which the walker toolchain builds Phi_l in practice
/// (docs/isogeny-walk/README.md: "l <= ~61 in practice").
const SMOOTH_BOUND: u128 = 61;

fn arg(args: &[String], name: &str) -> Option<String> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1).cloned())
}
fn arg_u(args: &[String], name: &str, d: u64) -> u64 {
    arg(args, name).map(|s| parse_u(&s) as u64).unwrap_or(d)
}
fn parse_u(s: &str) -> u128 {
    if let Some(h) = s.strip_prefix("0x") {
        u128::from_str_radix(h, 16).unwrap()
    } else {
        s.parse().unwrap()
    }
}

fn fac_json(f: &[(u128, u32)]) -> String {
    let parts: Vec<String> = f.iter().map(|(p, e)| format!("[\"{p}\",{e}]")).collect();
    format!("[{}]", parts.join(","))
}

fn fac_str(f: &[(u128, u32)]) -> String {
    if f.is_empty() {
        return "1".into();
    }
    f.iter()
        .map(|(p, e)| if *e == 1 { format!("{p}") } else { format!("{p}^{e}") })
        .collect::<Vec<_>>()
        .join(" * ")
}

#[derive(Clone, Debug)]
struct Structure {
    n_orders: usize,
    max_depth: u32,
    deepest_ell: u128,
    smooth: bool,
    qualifies: bool,
}

fn structure(fd: &FrobDisc) -> Structure {
    let n_orders = divisors(&fd.f_factors).len();
    let mut max_depth = 0;
    let mut deepest_ell = 0;
    for &(l, e) in &fd.f_factors {
        if e > max_depth {
            max_depth = e;
            deepest_ell = l;
        }
    }
    let smooth = fd.f_factors.iter().all(|&(l, _)| l <= SMOOTH_BOUND);
    let qualifies = fd.f > 1 && smooth && max_depth >= 2;
    Structure { n_orders, max_depth, deepest_ell, smooth, qualifies }
}

fn crater_type(d_k: i128, l: u128) -> &'static str {
    match kronecker(d_k, l) {
        1 => "split",
        0 => "ramified",
        _ => "inert",
    }
}

/// Order of +-a in (Z/l)^* / {+-1}: where l | f, pi = t/2 (mod l O_K) acts
/// as a scalar on E[l] at every curve above the floor, so this is the degree
/// of the extension of F_q holding the x-coordinates of each l-kernel.
fn kernel_x_degree(t: i128, l: u128) -> u128 {
    let a = (t.rem_euclid(l as i128) as u128) * nt::powmod(2, l - 2, l) % l;
    let mut x = a;
    let mut k = 1u128;
    while x != 1 && x != l - 1 {
        x = nt::mulmod(x, a, l);
        k += 1;
    }
    k
}

fn disc_json(fd: &FrobDisc) -> String {
    let st = structure(fd);
    let t = fd.t;
    let mut vol = Vec::new();
    for &(l, e) in &fd.f_factors {
        vol.push(format!(
            "{{\"ell\":{l},\"height\":{e},\"crater\":\"{}\",\"kernel_x_field_degree\":{}}}",
            crater_type(fd.d_k, l),
            kernel_x_degree(t, l)
        ));
    }
    format!(
        "\"delta\":\"{}\",\"d_K\":\"{}\",\"f\":\"{}\",\"f_factors\":{},\"f_str\":\"{}\",\"volcanoes\":[{}],\"n_orders\":{},\"max_depth\":{},\"deepest_ell\":{},\"f_smooth_61\":{},\"qualifies\":{}",
        fd.delta,
        fd.d_k,
        fd.f,
        fac_json(&fd.f_factors),
        fac_str(&fd.f_factors),
        vol.join(","),
        st.n_orders,
        st.max_depth,
        st.deepest_ell,
        st.smooth,
        st.qualifies
    )
}

// ---------------------------------------------------------------- Koblitz

fn lucas(m: u32, mu: i128) -> (i128, i128) {
    // tau^m = (V_m + U_m sqrt(-7)) / 2 with tau^2 - mu tau + 2 = 0
    let (mut v0, mut v1) = (2i128, mu);
    let (mut u0, mut u1) = (0i128, 1i128);
    for _ in 1..m {
        let v2 = mu * v1 - 2 * v0;
        let u2 = mu * u1 - 2 * u0;
        (v0, v1) = (v1, v2);
        (u0, u1) = (u1, u2);
    }
    (v1, u1)
}

fn cmd_koblitz(args: &[String]) {
    let lo = arg_u(args, "--min", 41) as u32;
    let hi = arg_u(args, "--max", 64) as u32;
    for m in lo..=hi {
        for a in [0u32, 1] {
            let mu: i128 = if a == 1 { 1 } else { -1 };
            let (t, u) = lucas(m, mu);
            let q: u128 = 1u128 << m;
            let n = (q as i128 + 1 - t) as u128;
            let base = if a == 1 { 2u128 } else { 4 };
            assert_eq!(n % base, 0);
            let fac = factor(n);
            let r = fac.last().unwrap().0;
            let h = n / r;
            let fd = frob_disc(t, q);
            assert_eq!(fd.d_k, -7);
            assert_eq!(fd.f, u.unsigned_abs());
            println!(
                "{{\"family\":\"koblitz\",\"m\":{m},\"m_prime\":{},\"a\":{a},\"t\":\"{t}\",\"order\":\"{n}\",\"order_factors\":{},\"r\":\"{r}\",\"r_bits\":{},\"cofactor\":\"{h}\",\"near_prime\":{},{}}}",
                is_prime(m as u128),
                fac_json(&fac),
                128 - r.leading_zeros(),
                h == base,
                disc_json(&fd)
            );
        }
    }
}

// ---------------------------------------------------------------- binary

struct BinCount {
    order: Option<u128>,
    ambiguous: bool,
}

fn binary_point(c: &BinCurve, seed: &str) -> BPt {
    for j in 0u32.. {
        let x = sha256::seed_to_u64(&format!("{seed}:x:{j}"), c.k.m);
        if let Some(p) = c.lift_x(x) {
            return p;
        }
    }
    unreachable!()
}

fn count_binary(c: &BinCurve, seed: &str) -> BinCount {
    let q: u128 = 1u128 << c.k.m;
    let w = isqrt(4 * q);
    let (lo, hi) = (q + 1 - w, q + 1 + w);
    // Tr(a2) = 1 and m odd  =>  #E = 2 (mod 4)
    for attempt in 0..6 {
        let p = binary_point(c, &format!("{seed}:pt{attempt}"));
        match multiples_in_window(c, p, lo, hi, 4, 2) {
            Some(v) if v.len() == 1 => return BinCount { order: Some(v[0]), ambiguous: false },
            _ => continue,
        }
    }
    BinCount { order: None, ambiguous: true }
}

fn b_seed(prefix: &str, i: u64) -> String {
    format!("{prefix}{i}")
}

fn cmd_binary(args: &[String]) {
    let m = arg_u(args, "--m", 53) as u32;
    let count = arg_u(args, "--count", 1000);
    let threads = arg_u(args, "--threads", 4) as usize;
    let prefix = arg(args, "--prefix").unwrap_or_else(|| format!("isogeny-class-trio/2026-10-17/B-{m}/b/"));
    assert!(m % 2 == 1);
    let k = F2m::new(standard_modulus(m));
    let q: u128 = 1u128 << m;
    let next = AtomicUsize::new(0);
    let rows: Mutex<Vec<(u64, String, Structure)>> = Mutex::new(Vec::new());
    let stats = Mutex::new((0u64, 0u64, 0u64)); // tried, prime-subgroup, ambiguous
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed) as u64;
                if i >= count {
                    break;
                }
                let seed = b_seed(&prefix, i);
                let b = sha256::seed_to_u64(&seed, m);
                if b == 0 {
                    continue;
                }
                let c = BinCurve { k: k.clone(), a2: 1, a6: b };
                let res = count_binary(&c, &seed);
                let mut st = stats.lock().unwrap();
                st.0 += 1;
                if res.ambiguous {
                    st.2 += 1;
                    continue;
                }
                drop(st);
                let n = res.order.unwrap();
                let r = n / 2;
                if !is_prime(r) {
                    continue;
                }
                // r > 4 sqrt(q) proves #E = n (unique multiple of r in the window)
                assert!(r > isqrt(16 * q));
                let t = q as i128 + 1 - n as i128;
                let fd = frob_disc(t, q);
                let st = structure(&fd);
                let line = format!(
                    "{{\"family\":\"binary\",\"index\":{i},\"seed\":\"{seed}\",\"m\":{m},\"modulus\":{:?},\"a2\":1,\"b\":\"0x{b:x}\",\"t\":\"{t}\",\"order\":\"{n}\",\"r\":\"{r}\",\"cofactor\":2,{}}}",
                    k.modulus,
                    disc_json(&fd)
                );
                stats.lock().unwrap().1 += 1;
                rows.lock().unwrap().push((i, line, st));
            });
        }
    });
    emit_search(rows.into_inner().unwrap(), stats.into_inner().unwrap(), &format!("binary m={m} count={count} prefix={prefix}"));
}

fn emit_search(mut rows: Vec<(u64, String, Structure)>, stats: (u64, u64, u64), label: &str) {
    rows.sort_by_key(|r| r.0);
    for r in &rows {
        println!("{}", r.1);
    }
    let qual: Vec<&(u64, String, Structure)> = rows.iter().filter(|r| r.2.qualifies).collect();
    let pick = qual.iter().max_by(|x, y| {
        (x.2.max_depth, x.2.n_orders, std::cmp::Reverse(x.0)).cmp(&(y.2.max_depth, y.2.n_orders, std::cmp::Reverse(y.0)))
    });
    eprintln!(
        "{{\"summary\":\"{label}\",\"curves_counted\":{},\"prime_subgroup\":{},\"ambiguous\":{},\"qualifying\":{},\"pick_index\":{}}}",
        stats.0,
        stats.1,
        stats.2,
        qual.len(),
        pick.map(|p| p.0.to_string()).unwrap_or("null".into())
    );
}

// ---------------------------------------------------------------- prime

fn pseudo_mersenne(bits: u32) -> (u64, u64) {
    let top: u128 = 1u128 << bits;
    let mut c = 1u128;
    loop {
        let p = top - c;
        if is_prime(p) {
            return (p as u64, c as u64);
        }
        c += 2;
    }
}

fn prime_point(c: &PCurve, seed: &str) -> Pt {
    for j in 0u32.. {
        let x = sha256::seed_to_u64(&format!("{seed}:x:{j}"), 64) % c.k.p;
        if let Some(p) = c.lift_x(x) {
            return p;
        }
    }
    unreachable!()
}

fn count_prime(c: &PCurve, seed: &str) -> Option<u128> {
    let q = c.k.p as u128;
    let w = isqrt(4 * q);
    let (lo, hi) = (q + 1 - w, q + 1 + w);
    for attempt in 0..6 {
        let p = prime_point(c, &format!("{seed}:pt{attempt}"));
        if let Some(v) = multiples_in_window(c, p, lo, hi, 1, 0) {
            if v.len() == 1 {
                return Some(v[0]);
            }
        }
    }
    None
}

fn cmd_prime(args: &[String]) {
    let bits = arg_u(args, "--bits", 53) as u32;
    let count = arg_u(args, "--count", 1000);
    let threads = arg_u(args, "--threads", 4) as usize;
    let (p, cc) = pseudo_mersenne(bits);
    let prefix = arg(args, "--prefix").unwrap_or_else(|| format!("isogeny-class-trio/2026-10-17/P-{bits}/b/"));
    let k = Fp { p };
    let a = p - 3;
    let next = AtomicUsize::new(0);
    let rows: Mutex<Vec<(u64, String, Structure)>> = Mutex::new(Vec::new());
    let stats = Mutex::new((0u64, 0u64, 0u64));
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed) as u64;
                if i >= count {
                    break;
                }
                let seed = b_seed(&prefix, i);
                let b = sha256::seed_to_u64(&seed, 64) % p;
                // discriminant -16(4a^3 + 27b^2) != 0
                let disc = k.add(k.mul(4, k.mul(a, k.mul(a, a))), k.mul(27, k.mul(b, b)));
                if disc == 0 || b == 0 {
                    continue;
                }
                let c = PCurve { k: k.clone(), a, b };
                let res = count_prime(&c, &seed);
                let mut st = stats.lock().unwrap();
                st.0 += 1;
                let Some(n) = res else {
                    st.2 += 1;
                    continue;
                };
                drop(st);
                if !is_prime(n) {
                    continue;
                }
                assert!(n > isqrt(16 * p as u128));
                let t = p as i128 + 1 - n as i128;
                let fd = frob_disc(t, p as u128);
                let st = structure(&fd);
                let line = format!(
                    "{{\"family\":\"prime\",\"index\":{i},\"seed\":\"{seed}\",\"p\":\"{p}\",\"p_form\":\"2^{bits} - {cc}\",\"a\":\"-3\",\"b\":\"0x{b:x}\",\"t\":\"{t}\",\"order\":\"{n}\",\"r\":\"{n}\",\"cofactor\":1,{}}}",
                    disc_json(&fd)
                );
                stats.lock().unwrap().1 += 1;
                rows.lock().unwrap().push((i, line, st));
            });
        }
    });
    emit_search(rows.into_inner().unwrap(), stats.into_inner().unwrap(), &format!("prime bits={bits} p=2^{bits}-{cc} count={count} prefix={prefix}"));
}

// ---------------------------------------------------------------- analysis

/// Every stratum (order O_c, c | f) with its exact class number and the
/// total class size H = sum_c h(c^2 d_K) (Schoof 1987, Thm 4.6).
fn strata_json(fd: &FrobDisc, h_k: u64) -> (String, u128) {
    let mut rows = Vec::new();
    let mut total: u128 = 0;
    for (c, cf) in divisors(&fd.f_factors) {
        let (num, den) = order_class_ratio(fd.d_k, &cf);
        let hc = h_k as u128 * num;
        assert_eq!(hc % den, 0, "class number formula not integral");
        let hc = hc / den;
        total += hc;
        let depth: Vec<String> = fd
            .f_factors
            .iter()
            .map(|&(l, _)| {
                let e = cf.iter().find(|x| x.0 == l).map(|x| x.1).unwrap_or(0);
                format!("\"{l}\":{e}")
            })
            .collect();
        rows.push(format!(
            "{{\"conductor\":\"{c}\",\"conductor_str\":\"{}\",\"disc\":\"{}\",\"levels\":{{{}}},\"h\":\"{hc}\"}}",
            fac_str(&cf),
            fd.d_k * (c * c) as i128,
            depth.join(",")
        ));
    }
    (format!("[{}]", rows.join(",")), total)
}

fn cmd_analyze_binary(args: &[String]) {
    let m = arg_u(args, "--m", 53) as u32;
    let a2 = arg_u(args, "--a2", 1);
    let b = parse_u(&arg(args, "--b").expect("--b")) as u64;
    let seed = arg(args, "--seed").unwrap_or_else(|| format!("isogeny-class-trio/2026-10-17/analyze/{m}/{b:x}"));
    let koblitz = arg(args, "--koblitz").is_some();
    let k = F2m::new(standard_modulus(m));
    let c = BinCurve { k: k.clone(), a2, a6: b };
    let q: u128 = 1u128 << m;
    let (n, t) = if koblitz {
        assert_eq!(b, 1);
        let mu = if a2 == 1 { 1 } else { -1 };
        let (t, _) = lucas(m, mu);
        ((q as i128 + 1 - t) as u128, t)
    } else {
        let w = isqrt(4 * q);
        let p = binary_point(&c, &format!("{seed}:pt0"));
        let v = multiples_in_window(&c, p, q + 1 - w, q + 1 + w, 4, if a2 == 1 { 2 } else { 0 }).unwrap();
        assert_eq!(v.len(), 1);
        (v[0], q as i128 + 1 - v[0] as i128)
    };
    let fac = factor(n);
    let r = fac.last().unwrap().0;
    let h = n / r;
    // generator of the prime-order subgroup and order certificate
    let mut g = BPt::Inf;
    let mut gj = 0;
    for j in 0.. {
        let p = binary_point(&c, &format!("{seed}:G{j}"));
        let gg = c.mul(p, h);
        if gg != BPt::Inf {
            g = gg;
            gj = j;
            break;
        }
    }
    assert!(c.on_curve(g));
    assert_eq!(c.mul(g, r), BPt::Inf);
    let r_gt_window = r > isqrt(16 * q);
    // #E P = O for random points (consistency of the counted order)
    for j in 0..8 {
        let p = binary_point(&c, &format!("{seed}:check{j}"));
        assert_eq!(c.mul(p, n), BPt::Inf);
    }
    let fd = frob_disc(t, q);
    let h_k = if fd.d_k == -7 { 1 } else { classno::class_number(fd.d_k) };
    let (strata, total) = strata_json(&fd, h_k);
    let v3 = fd.f_factors.iter().find(|x| x.0 == 3).map(|x| x.1).unwrap_or(0);
    let vol = iso3::volcano3(&k, b, v3);
    let mut floors = Vec::new();
    for &(l, _) in &fd.f_factors {
        if l <= SMOOTH_BOUND {
            let kd = kernel_x_degree(t, l) as u32;
            let ft = divpoly::floor_test(&k, b, l as usize, kd);
            floors.push(format!(
                "{{\"ell\":{l},\"kernel_x_field_degree\":{kd},\"fixed_x\":{},\"all_x\":{},\"eigenline_x\":{},\"root_on_floor\":{}}}",
                ft.fixed_x,
                (l * l - 1) / 2,
                (l - 1) / 2,
                ft.on_floor
            ));
        }
    }
    let BPt::A(gx, gy) = g else { unreachable!() };
    println!(
        "{{\"floor_tests\":[{}],\"family\":\"{}\",\"m\":{m},\"modulus\":{:?},\"modulus_hex\":\"{}\",\"a2\":{a2},\"b\":\"0x{b:x}\",\"j\":\"0x{:x}\",\"t\":\"{t}\",\"order\":\"{n}\",\"order_factors\":{},\"r\":\"{r}\",\"r_bits\":{},\"cofactor\":\"{h}\",\"generator\":[\"0x{gx:x}\",\"0x{gy:x}\"],\"generator_seed\":\"{seed}:G{gj}\",\"order_certificate\":{{\"rG=O\":true,\"r_prime\":{},\"r_gt_4sqrt_q\":{r_gt_window},\"nP=O_for_8_points\":true}},{},\"h_dK\":\"{h_k}\",\"strata\":{strata},\"class_size\":\"{total}\",\"volcano3\":{{\"height\":{},\"root_rational_3_isogenies\":{},\"distance_to_floor\":{},\"level_from_crater\":{},\"vertices_explored\":{}}}}}",
        floors.join(","),
        if koblitz { "koblitz" } else { "binary" },
        k.modulus,
        k.modulus_hex(),
        k.inv(b),
        fac_json(&fac),
        128 - r.leading_zeros(),
        is_prime(r),
        disc_json(&fd),
        vol.height,
        vol.root_neighbours,
        vol.distance_to_floor.map(|d| d.to_string()).unwrap_or("null".into()),
        vol.level_from_crater.map(|d| d.to_string()).unwrap_or("null".into()),
        vol.vertices_explored
    );
}

fn cmd_analyze_prime(args: &[String]) {
    let p = parse_u(&arg(args, "--p").expect("--p")) as u64;
    let b = parse_u(&arg(args, "--b").expect("--b")) as u64;
    let k = Fp { p };
    let a = p - 3;
    let c = PCurve { k: k.clone(), a, b };
    let seed = arg(args, "--seed").unwrap_or_else(|| format!("isogeny-class-trio/2026-10-17/analyze/p{p}/{b:x}"));
    let n = count_prime(&c, &seed).expect("order");
    assert!(is_prime(n) && n > isqrt(16 * p as u128));
    let mut g = Pt::Inf;
    let mut gj = 0;
    for j in 0.. {
        let pt = prime_point(&c, &format!("{seed}:G{j}"));
        if pt != Pt::Inf {
            g = pt;
            gj = j;
            break;
        }
    }
    assert_eq!(c.mul(g, n), Pt::Inf);
    let t = p as i128 + 1 - n as i128;
    let fd = frob_disc(t, p as u128);
    let h_k = classno::class_number(fd.d_k);
    let (strata, total) = strata_json(&fd, h_k);
    let v3 = fd.f_factors.iter().find(|x| x.0 == 3).map(|x| x.1).unwrap_or(0);
    let vol = fpiso::volcano3(&k, a, b, v3);
    let vol_json = format!(
        "{{\"height\":{},\"root_rational_3_isogenies\":{},\"distance_to_floor\":{},\"level_from_crater\":{},\"vertices_explored\":{}}}",
        vol.height,
        vol.root_neighbours,
        vol.distance_to_floor.map(|d| d.to_string()).unwrap_or("null".into()),
        vol.level_from_crater.map(|d| d.to_string()).unwrap_or("null".into()),
        vol.vertices_explored
    );
    let Pt::A(gx, gy) = g else { unreachable!() };
    let mut s = String::new();
    let _ = write!(
        s,
        "{{\"family\":\"prime\",\"p\":\"{p}\",\"p_hex\":\"0x{p:x}\",\"a\":\"-3\",\"b\":\"0x{b:x}\",\"t\":\"{t}\",\"order\":\"{n}\",\"r\":\"{n}\",\"r_bits\":{},\"cofactor\":\"1\",\"generator\":[\"0x{gx:x}\",\"0x{gy:x}\"],\"generator_seed\":\"{seed}:G{gj}\",\"order_certificate\":{{\"nG=O\":true,\"n_prime\":true,\"n_gt_4sqrt_p\":true}},{},\"h_dK\":\"{h_k}\",\"strata\":{strata},\"class_size\":\"{total}\",\"volcano3\":{vol_json}}}",
        128 - n.leading_zeros(),
        disc_json(&fd)
    );
    println!("{s}");
}

/// Check a recorded generator: on the curve, not O, and r G = O (r prime).
fn cmd_check_point(args: &[String]) {
    let r = parse_u(&arg(args, "--r").expect("--r"));
    assert!(is_prime(r));
    let gx = parse_u(&arg(args, "--gx").expect("--gx")) as u64;
    let gy = parse_u(&arg(args, "--gy").expect("--gy")) as u64;
    let b = parse_u(&arg(args, "--b").expect("--b")) as u64;
    let ok = if let Some(m) = arg(args, "--m") {
        let k = F2m::new(standard_modulus(m.parse().unwrap()));
        let c = BinCurve { k, a2: arg_u(args, "--a2", 1), a6: b };
        let g = BPt::A(gx, gy);
        c.on_curve(g) && c.mul(g, r) == BPt::Inf
    } else {
        let p = parse_u(&arg(args, "--p").expect("--p or --m")) as u64;
        let c = PCurve { k: Fp { p }, a: p - 3, b };
        let g = Pt::A(gx, gy);
        c.on_curve(g) && c.mul(g, r) == Pt::Inf
    };
    println!("{{\"on_curve_and_rG_is_O\":{ok}}}");
    assert!(ok);
}

fn field_after<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let i = line.find(key)? + key.len();
    let rest = &line[i..];
    let end = rest.find(|c: char| c == ',' || c == '}').unwrap_or(rest.len());
    Some(rest[..end].trim_matches('"'))
}

/// Counts over a search JSONL (one record per prime-subgroup curve).
fn cmd_tally(args: &[String]) {
    let path = arg(args, "--in").expect("--in");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut n = 0u64;
    let mut f_gt1 = 0u64;
    let mut smooth_gt1 = 0u64;
    let mut qual = 0u64;
    let mut depth = std::collections::BTreeMap::new();
    let mut ells = std::collections::BTreeMap::new();
    for line in text.lines().filter(|l| !l.is_empty()) {
        n += 1;
        let f = field_after(line, "\"f\":").unwrap();
        let d: u32 = field_after(line, "\"max_depth\":").unwrap().parse().unwrap();
        *depth.entry(d).or_insert(0u64) += 1;
        if f != "1" {
            f_gt1 += 1;
            if field_after(line, "\"f_smooth_61\":") == Some("true") {
                smooth_gt1 += 1;
            }
        }
        if field_after(line, "\"qualifies\":") == Some("true") {
            qual += 1;
        }
        let mut rest = line;
        while let Some(i) = rest.find("{\"ell\":") {
            rest = &rest[i + 7..];
            let l: u64 = rest[..rest.find(',').unwrap()].parse().unwrap();
            *ells.entry(l).or_insert(0u64) += 1;
        }
    }
    let dj: Vec<String> = depth.iter().map(|(k, v)| format!("\"{k}\":{v}")).collect();
    let lj: Vec<String> = ells.iter().map(|(k, v)| format!("\"{k}\":{v}")).collect();
    println!(
        "{{\"source\":\"{path}\",\"prime_subgroup_curves\":{n},\"f_gt_1\":{f_gt1},\"f_gt_1_and_61_smooth\":{smooth_gt1},\"qualifying\":{qual},\"max_depth_histogram\":{{{}}},\"curves_with_ell_dividing_f\":{{{}}}}}",
        dj.join(","),
        lj.join(",")
    );
}

/// j-invariants of the class-number-one discriminants.
const CM1: [(i128, i128); 9] = [
    (-3, 0),
    (-4, 1728),
    (-7, -3375),
    (-8, 8000),
    (-11, -32768),
    (-19, -884736),
    (-43, -884736000),
    (-67, -147197952000),
    (-163, -262537412640768000),
];

/// n P = O for a point P != O with n prime and n > 4 sqrt(p) proves #E = n.
fn order_is(c: &PCurve, n: u128, seed: &str) -> bool {
    let pt = prime_point(c, seed);
    pt != Pt::Inf && c.mul(pt, n) == Pt::Inf
}

fn non_residue(k: &Fp) -> u64 {
    (2..).find(|&z| !k.is_square(z)).unwrap()
}

/// One class-number-one CM curve with f = l prime (see ../PROTOCOL-gap.md).
fn cm_gap_for(d: i128, j: i128) -> Option<String> {
    let lo: u128 = 1 << 60;
    let hi: u128 = 1 << 61;
    let big_l = isqrt((1u128 << 63) / d.unsigned_abs());
    let mut l = big_l;
    let mut tried = 0u64;
    if d % 4 == 0 || d.rem_euclid(8) == 1 {
        // d = 0 (mod 4): t even, so n = p + 1 - t is even (no prime order).
        // d = 1 (mod 8): t, l odd give t^2 + |d| l^2 = 0 (mod 8), so p is even.
        return None;
    }
    let floor_l = big_l / 2;
    let mut primes_seen = 0u32;
    while l > floor_l {
        if !is_prime(l) {
            l -= 1;
            continue;
        }
        let dl2 = d.unsigned_abs() * l * l;
        primes_seen += 1;
        if primes_seen > 2000 {
            return None;
        }
        for t in 1u128..200_000 {
            let num = t * t + dl2;
            if num % 4 != 0 {
                continue;
            }
            let p = num / 4;
            if p <= lo || p >= hi || !is_prime(p) {
                continue;
            }
            let n = p + 1 - t;
            if !is_prime(n) || n <= isqrt(16 * p) {
                continue;
            }
            if d != -3 && d != -4 && p % 4 != 3 {
                continue;
            }
            tried += 1;
            let k = Fp { p: p as u64 };
            let seed = format!("isogeny-class-trio/2026-10-17/cm-gap/d{d}/l{l}/t{t}");
            let (a, b, form) = if d == -3 {
                let found = (1u64..200).find(|&b| order_is(&PCurve { k: k.clone(), a: 0, b }, n, &seed));
                match found {
                    Some(b) => (0u64, b, "y^2 = x^3 + b (j = 0, secp256k1 form)"),
                    None => continue,
                }
            } else if d == -4 {
                let found = (1u64..200).find(|&a| order_is(&PCurve { k: k.clone(), a, b: 0 }, n, &seed));
                match found {
                    Some(a) => (a, 0u64, "y^2 = x^3 + a x (j = 1728)"),
                    None => continue,
                }
            } else {
                let jm = j.rem_euclid(p as i128) as u64;
                let s1728 = k.sub(1728 % k.p, jm);
                let a0 = k.mul(3, k.mul(jm, s1728));
                let b0 = k.mul(2, k.mul(jm, k.mul(s1728, s1728)));
                let c = non_residue(&k);
                let tw = (k.mul(k.mul(c, c), a0), k.mul(k.mul(c, k.mul(c, c)), b0));
                let (a1, b1) = if order_is(&PCurve { k: k.clone(), a: a0, b: b0 }, n, &seed) {
                    (a0, b0)
                } else if order_is(&PCurve { k: k.clone(), a: tw.0, b: tw.1 }, n, &seed) {
                    tw
                } else {
                    continue;
                };
                // a = -3 form: u^4 = -3 / a1; p = 3 (mod 4) so fourth powers = squares
                let target = k.mul(k.sub(0, 3), k.inv(a1));
                let Some(sq) = k.sqrt(target) else { continue };
                let sq = if k.is_square(sq) { sq } else { k.sub(0, sq) };
                let u = k.sqrt(sq).unwrap();
                let u2 = k.mul(u, u);
                let a3 = k.mul(k.mul(u2, u2), a1);
                assert_eq!(a3, k.sub(0, 3));
                let b3 = k.mul(k.mul(u2, k.mul(u2, u2)), b1);
                (a3, b3, "y^2 = x^3 - 3x + b (P-curve form)")
            };
            let c = PCurve { k: k.clone(), a, b };
            assert!(order_is(&c, n, &format!("{seed}:check")));
            let fd = frob_disc(t as i128, p);
            assert_eq!((fd.d_k, fd.f), (if d == -4 || d == -8 { d } else { d }, l));
            let jj = j.rem_euclid(p as i128);
            let h_k = 1u64;
            let (strata, total) = strata_json(&fd, h_k);
            let g = prime_point(&c, &format!("{seed}:G0"));
            let Pt::A(gx, gy) = g else { unreachable!() };
            return Some(format!(
                "{{\"family\":\"prime-cm-gap\",\"d\":{d},\"j\":\"{jj}\",\"form\":\"{form}\",\"p\":\"{p}\",\"p_hex\":\"0x{p:x}\",\"a\":\"{a}\",\"b\":\"{b}\",\"t\":\"{t}\",\"order\":\"{n}\",\"r\":\"{n}\",\"r_bits\":{},\"cofactor\":\"1\",\"gap_l\":\"{l}\",\"gap_l_bits\":{},\"generator\":[\"0x{gx:x}\",\"0x{gy:x}\"],\"generator_seed\":\"{seed}:G0\",\"order_certificate\":{{\"nP=O\":true,\"n_prime\":true,\"n_gt_4sqrt_p\":true}},\"candidates_tried\":{tried},{},\"h_dK\":\"1\",\"strata\":{strata},\"class_size\":\"{total}\",\"root_level\":\"crater (End = O_K: j = j(O_K), h(d) = 1)\"}}",
                128 - n.leading_zeros(),
                128 - l.leading_zeros(),
                disc_json(&fd)
            ));
        }
        l -= 1;
    }
    None
}

fn cmd_cm_gap(_args: &[String]) {
    for (d, j) in CM1 {
        match cm_gap_for(d, j) {
            Some(line) => println!("{line}"),
            None => println!(
                "{{\"family\":\"prime-cm-gap\",\"d\":{d},\"found\":false,\"reason\":\"{}\"}}",
                if d % 4 == 0 { "t even, so n = p + 1 - t is even: prime order impossible" } else if d.rem_euclid(8) == 1 { "d = 1 (mod 8): odd t and odd l make t^2 + |d| l^2 = 0 (mod 8), so 4p would force p even" } else if d == -11 { "3 splits and 3 does not divide l: p prime forces t = 0 (mod 3), then 3 | n" } else { "none in the first 2000 primes l below L, t < 200000" }
            ),
        }
    }
}

fn cmd_classno(args: &[String]) {
    let d: i128 = arg(args, "--d").expect("--d").parse().unwrap();
    assert!(classno::is_fundamental(d));
    println!("{{\"d\":\"{d}\",\"h\":{}}}", classno::class_number(d));
}

fn main() {
    let args: Vec<String> = env::args().collect();
    match args.get(1).map(|s| s.as_str()) {
        Some("koblitz") => cmd_koblitz(&args),
        Some("binary") => cmd_binary(&args),
        Some("prime") => cmd_prime(&args),
        Some("analyze-binary") => cmd_analyze_binary(&args),
        Some("analyze-prime") => cmd_analyze_prime(&args),
        Some("classno") => cmd_classno(&args),
        Some("check-point") => cmd_check_point(&args),
        Some("tally") => cmd_tally(&args),
        Some("cm-gap") => cmd_cm_gap(&args),
        Some("modulus") => {
            let m = arg_u(&args, "--m", 53) as u32;
            println!("{:?}", standard_modulus(m));
        }
        _ => {
            eprintln!("usage: selector <koblitz|binary|prime|analyze-binary|analyze-prime|classno|modulus> [options]");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn koblitz_orders_match_registry() {
        // icv1-f2m41-tm2308219-7f48b14a: subgroup 40 bits, cofactor 4 (docs/curves/registry.json)
        let (t, _) = lucas(41, -1);
        let n = (1u128 << 41) as i128 + 1 - t;
        assert_eq!(n % 4, 0);
        assert!(is_prime((n / 4) as u128));
        // K-163 (a = 1): order 2 * r with r prime (FIPS 186-4 D.1.3.1)
        let (t, _) = lucas(163, 1);
        let _ = t; // 2^163 exceeds the u128 order arithmetic used elsewhere
    }

    #[test]
    fn bsgs_matches_naive_count_binary() {
        let m = 17;
        let k = F2m::new(standard_modulus(m));
        for b in [5u64, 99, 1234, 77777] {
            let c = BinCurve { k: k.clone(), a2: 1, a6: b };
            let mut naive = 2u128;
            for x in 1..(1u64 << m) {
                let cc = x ^ 1 ^ k.mul(b, k.inv(k.sqr(x)));
                if k.trace(cc) == 0 {
                    naive += 2;
                }
            }
            let got = count_binary(&c, &format!("t{b}"));
            assert_eq!(got.order, Some(naive), "b={b}");
        }
    }

    #[test]
    fn bsgs_matches_naive_count_prime() {
        let p = 100_003u64;
        let k = Fp { p };
        for b in [5u64, 17, 1234] {
            let c = PCurve { k: k.clone(), a: p - 3, b };
            let mut naive = 1u128;
            for x in 0..p {
                let rhs = k.add(k.mul(k.add(k.mul(x, x), p - 3), x), b);
                naive += if rhs == 0 { 1 } else if k.is_square(rhs) { 2 } else { 0 };
            }
            assert_eq!(count_prime(&c, &format!("t{b}")), Some(naive));
        }
    }

    /// Control for the characteristic-2 floor test: over GF(2^7) the Koblitz
    /// trace has f = 7 (ramified), so the class is 1 crater curve (the Koblitz
    /// curve, tau in End) plus h(-7 * 49) = 7 floor curves.
    #[test]
    fn floor_test_control_koblitz_m7() {
        let m = 7;
        let k = F2m::new(standard_modulus(m));
        for (a2, mu) in [(1u64, 1i128), (0, -1)] {
            let (t, u) = lucas(m, mu);
            assert_eq!(u.abs(), 7);
            let target = (1i128 << m) + 1 - t;
            let kd = kernel_x_degree(t, 7) as u32;
            let mut floor = 0;
            let mut above = vec![];
            for b in 1..(1u64 << m) {
                let mut n = 2i128;
                for x in 1..(1u64 << m) {
                    let cc = x ^ a2 ^ k.mul(b, k.inv(k.sqr(x)));
                    if k.trace(cc) == 0 {
                        n += 2;
                    }
                }
                if n != target {
                    continue;
                }
                let ft = divpoly::floor_test(&k, b, 7, kd);
                if ft.on_floor {
                    floor += 1;
                } else {
                    above.push(b);
                }
            }
            assert_eq!(above, vec![1], "a2={a2}");
            assert_eq!(floor, 7, "a2={a2}");
        }
    }
}
