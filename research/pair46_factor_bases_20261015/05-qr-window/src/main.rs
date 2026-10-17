// 05-qr-window: twist-aware quadratic-character factor base (positive control).
// Usage: qr_window <out_dir> <T_targets> <n_toy_F> <n_toy_N2> <boots>
// All results are computed here; Python only plots and lays out the PDF.
#![allow(non_snake_case, dead_code)]
use std::fmt::Write as FmtWrite;
use std::io::Write;
use std::time::Instant;

include!("common.rs");

const SEED: u64 = 46005;
const HSEED: u64 = 0x5EED_0005_0A0A_4605;
const ARMS: [&str; 5] = ["A0_plain", "A1_qr_window_chi_f", "A2_spec_chi_x", "A3_random_half", "A4_twist_only"];

fn splitmix(x: u64) -> u64 { let mut z = x.wrapping_add(0x9E3779B97F4A7C15); z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }

/// Membership predicate of each arm (the single instrument): inputs are chi(f(x)), chi(x), x.
#[inline]
fn member(arm: usize, chif: i32, chix: i32, x: u64) -> bool {
    match arm {
        0 => true,
        1 => chif == 1,
        2 => chix != -1,
        3 => splitmix(x ^ HSEED) & 1 == 0,
        4 => chif == -1,
        _ => unreachable!(),
    }
}

#[derive(Clone, Copy, Default)]
struct WinCount { v: [u64; 5], pts: [u64; 5], chif0: u64, len: u64 }

fn count_window(E: &Curve, lo: u64, len: u64) -> WinCount {
    let p = E.p; let mut w = WinCount::default(); w.len = len;
    for x in lo..lo + len {
        let chif = jacobi(E.rhs(x), p); let chix = jacobi(x, p);
        if chif == 0 { w.chif0 += 1; }
        let is_pt = chif != -1; // x-coordinate of an E-point (incl. 2-torsion if chif==0)
        for a in 0..5 { if member(a, chif, chix, x) { w.v[a] += 1; if is_pt { w.pts[a] += 1; } } }
    }
    w
}

fn wilson(k: u64, n: u64) -> (f64, f64) {
    let z = 1.959963984540054; let n_ = n as f64; let ph = k as f64 / n_;
    let d = 1.0 + z * z / n_; let c = (ph + z * z / (2.0 * n_)) / d;
    let h = z * ((ph * (1.0 - ph) / n_ + z * z / (4.0 * n_ * n_)).sqrt()) / d;
    (c - h, c + h)
}
fn pct(v: &mut Vec<f64>, q: f64) -> f64 { v.sort_by(|a, b| a.partial_cmp(b).unwrap()); let i = ((v.len() - 1) as f64 * q).round() as usize; v[i] }

struct Toy { arm: String, bits: u32, p: u64, N: u64, B: u64, v: [u64; 5], pts: [u64; 5], rel: [u64; 5], dec: [u64; 5], e1: Option<bool>, e2: bool, T: u64 }

fn toy_measure(arm_name: &str, E: &Curve, N: u64, T: u64, rng: &mut Rng) -> Toy {
    let p = E.p;
    let B = ((2.0 * N as f64).sqrt()).ceil() as u64;
    let mut memb = vec![vec![false; B as usize]; 5];
    let (mut v, mut pts) = ([0u64; 5], [0u64; 5]);
    let mut fb: Vec<(u64, u64)> = vec![];
    let mut e2 = true;
    for x in 0..B {
        let chif = jacobi(E.rhs(x), p); let chix = jacobi(x, p);
        for a in 0..5 { let m = member(a, chif, chix, x); memb[a][x as usize] = m; if m { v[a] += 1; if chif != -1 { pts[a] += 1; } } }
        if chif != -1 {
            if memb[0][x as usize] != memb[1][x as usize] { e2 = false; }
            fb.push((x, sqrtm(E.rhs(x), p).unwrap()));
        }
    }
    let (mut rel, mut dec) = ([0u64; 5], [0u64; 5]);
    let G = Some(E.rand_point(rng));
    for _ in 0..T {
        let k = 1 + rng.below(N - 1);
        let R = E.mul(k, G);
        let rx = R.unwrap().0;
        let mut c = [0u64; 5];
        for &(x1, y1) in &fb {
            if x1 == rx { continue; }
            let P = Some((x1, y1));
            for Q in [E.add(R, E.neg(P)), E.add(R, P)] {
                if let Some((x2, _)) = Q {
                    if x2 < B && x2 >= x1 {
                        for a in 0..5 { if memb[a][x1 as usize] && memb[a][x2 as usize] { c[a] += 1; } }
                    }
                }
            }
        }
        if c[0] != c[1] { e2 = false; }
        for a in 0..5 { dec[a] += c[a]; if c[a] > 0 { rel[a] += 1; } }
    }
    let e1 = if p < (1u64 << 23) {
        let mut k1 = 0u64; let mut km = 0u64; let mut k0 = 0u64;
        for x in 0..p { match jacobi(E.rhs(x), p) { 1 => k1 += 1, -1 => km += 1, _ => k0 += 1 } }
        let Nt = 2 * p + 2 - N;
        Some(k0 == 0 && k1 == (N - 1) / 2 && km == (Nt - 1) / 2)
    } else { None };
    Toy { arm: arm_name.to_string(), bits: 64 - p.leading_zeros(), p, N, B, v, pts, rel, dec, e1, e2, T }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out = &args[1]; let T: u64 = args[2].parse().unwrap();
    let nF: usize = args[3].parse().unwrap(); let nN2: usize = args[4].parse().unwrap();
    let boots: usize = args[5].parse().unwrap();
    let t0 = Instant::now();
    let mut rng = Rng(SEED);
    let mut log = std::fs::File::create(format!("{}/run.log", out)).unwrap();
    macro_rules! say { ($($a:tt)*) => {{ let s = format!($($a)*); println!("{}", s); writeln!(log, "{}", s).unwrap(); log.flush().unwrap(); }} }
    say!("qr_window seed={} T={} nF={} nN2={} boots={}", SEED, T, nF, nN2, boots);
    let mut js = String::new();
    write!(js, "{{\"item\":\"05-qr-window\",\"seed\":{},\"arms\":{:?},", SEED, ARMS).unwrap();

    // ---------------- pair 46 ----------------
    let p46: u64 = 5947698030651343; let t46: u64 = 34489701; let N46: u64 = 5947697996161643;
    assert!(p46 == 5 * t46 * t46 + 19 * t46 + 19 && N46 == p46 + 1 - t46 && is_prime(p46) && is_prime(N46));
    let curves = [("crater", Curve { p: p46, a: 5757116383569871, b: 3560173328417442 }),
                  ("floor", Curve { p: p46, a: 2087370704375809, b: 1402887937802154 })];
    let B46 = ((2.0 * N46 as f64).sqrt()).ceil() as u64;
    let Nt46 = 2 * p46 + 2 - N46;
    say!("pair46 B=ceil(sqrt(2N))={} twist order={} global exact gain 2p/(N-1)={:.12}", B46, Nt46, 2.0 * p46 as f64 / (N46 - 1) as f64);
    write!(js, "\"pair46\":{{\"p\":{},\"N\":{},\"twist_order\":{},\"B\":{},\"exact_liftable_x\":{},\"exact_twist_x\":{},\"exact_global_gain\":{:.15},\"curves\":[", p46, N46, Nt46, B46, (N46 - 1) / 2, (Nt46 - 1) / 2, 2.0 * p46 as f64 / (N46 - 1) as f64).unwrap();
    let mut wrng = Rng(SEED ^ 0xA5A5);
    let offs: Vec<u64> = (0..32).map(|_| wrng.below(p46 - 65536)).collect();
    for (ci, (name, E)) in curves.iter().enumerate() {
        let ok = order_is(E, N46, &mut rng);
        say!("{}: N*P=O on 3 random points: {}", name, ok);
        assert!(ok);
        let tc = Instant::now();
        let w = count_window(E, 0, B46);
        say!("{} canonical [0,B): chif0={} |V|={:?} pts={:?} ({:.1}s)", name, w.chif0, w.v, w.pts, tc.elapsed().as_secs_f64());
        let m1: Vec<f64> = (0..5).map(|a| w.pts[a] as f64 / w.v[a] as f64).collect();
        let g: Vec<f64> = (0..5).map(|a| m1[a] / m1[0]).collect();
        // canonical-window interval for G(A1) = B/k1: invert Wilson interval of k1/B
        let (lo0, hi0) = wilson(w.pts[0], w.v[0]);
        let gA1 = (1.0 / hi0, 1.0 / lo0);
        let ci2 = { let (l, h) = wilson(w.pts[2], w.v[2]); (l / hi0, h / lo0) };
        let ci3 = { let (l, h) = wilson(w.pts[3], w.v[3]); (l / hi0, h / lo0) };
        say!("{} canonical M1={:?} G={:?} G(A1) CI={:?} G(A2) CI~{:?} G(A3) CI~{:?}", name, m1, g, gA1, ci2, ci3);
        // random windows
        let ws: Vec<WinCount> = offs.iter().map(|&o| count_window(E, o, 65536)).collect();
        let gain = |idx: &[usize], a: usize| -> f64 {
            let (mut va, mut pa, mut v0, mut p0) = (0u64, 0u64, 0u64, 0u64);
            for &i in idx { va += ws[i].v[a]; pa += ws[i].pts[a]; v0 += ws[i].v[0]; p0 += ws[i].pts[0]; }
            if va == 0 { return f64::NAN; } (pa as f64 / va as f64) / (p0 as f64 / v0 as f64)
        };
        let all: Vec<usize> = (0..32).collect();
        let mut brng = Rng(SEED ^ (ci as u64 + 1));
        let mut wg = vec![];
        for a in 0..5 {
            let pt = gain(&all, a);
            let mut bs: Vec<f64> = (0..boots).map(|_| { let idx: Vec<usize> = (0..32).map(|_| brng.below(32) as usize).collect(); gain(&idx, a) }).collect();
            let (l, h) = (pct(&mut bs, 0.025), pct(&mut bs, 0.975));
            wg.push((pt, l, h));
        }
        let chif0w: u64 = ws.iter().map(|w| w.chif0).sum();
        say!("{} 32 random windows (2^16): G[point,lo,hi]={:?} chif0={}", name, wg, chif0w);
        if ci > 0 { js.push(','); }
        write!(js, "{{\"name\":\"{}\",\"a\":{},\"b\":{},\"order_check\":{},\"canonical\":{{\"B\":{},\"chif0\":{},\"V\":{:?},\"pts\":{:?},\"M1\":{:?},\"G\":{:?},\"G_A1_ci95\":[{},{}],\"G_A2_ci95\":[{},{}],\"G_A3_ci95\":[{},{}]}},",
            name, E.a, E.b, ok, B46, w.chif0, w.v, w.pts, m1, g, gA1.0, gA1.1, ci2.0, ci2.1, ci3.0, ci3.1).unwrap();
        write!(js, "\"windows\":{{\"n\":32,\"len\":65536,\"offsets\":{:?},\"chif0\":{},\"V\":{:?},\"pts\":{:?},\"G\":[", offs, chif0w,
            ws.iter().map(|w| w.v).collect::<Vec<_>>(), ws.iter().map(|w| w.pts).collect::<Vec<_>>()).unwrap();
        for (a, (pt, l, h)) in wg.iter().enumerate() { if a > 0 { js.push(','); } if pt.is_nan() { write!(js, "{{\"point\":null}}").unwrap(); } else { write!(js, "{{\"point\":{},\"ci95\":[{},{}]}}", pt, l, h).unwrap(); } }
        js.push_str("]}}");
    }
    js.push_str("]},");
    say!("pair46 done at {:.1}s", t0.elapsed().as_secs_f64());

    // ---------------- toy ----------------
    let csv = std::fs::read_to_string(format!("{}/input_curves_from_04.csv", out)).unwrap();
    let mut sel: Vec<(String, u64, u64, u64, u64)> = vec![];
    let (mut cf, mut cn) = (0, 0);
    for line in csv.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        let arm = f[0];
        let take = (arm == "F" && cf < nF) || (arm == "N2" && cn < nN2);
        if !take { continue; }
        if arm == "F" { cf += 1 } else { cn += 1 }
        sel.push((arm.to_string(), f[2].parse().unwrap(), f[3].parse().unwrap(), f[6].parse().unwrap(), f[7].parse().unwrap()));
    }
    say!("toy curves selected: F={} N2={}", cf, cn);
    let mut toys: Vec<Toy> = vec![];
    for (i, (arm, p, N, a, b)) in sel.iter().enumerate() {
        let E = Curve { p: *p, a: *a, b: *b };
        let mut trng = Rng(SEED ^ (i as u64));
        assert!(order_is(&E, *N, &mut trng), "toy order check failed");
        let t = toy_measure(arm, &E, *N, T, &mut trng);
        if i % 20 == 0 { say!("toy {} {} bits={} B={} dec={:?} rel={:?} e1={:?} e2={} t={:.1}s", i, arm, t.bits, t.B, t.dec, t.rel, t.e1, t.e2, t0.elapsed().as_secs_f64()); }
        toys.push(t);
    }
    let e2_all = toys.iter().all(|t| t.e2);
    let e1_checked = toys.iter().filter(|t| t.e1.is_some()).count();
    let e1_all = toys.iter().all(|t| t.e1 != Some(false));
    say!("toy E2 (A1 == A0 per target, all curves): {}; E1 exact full-field counts: {}/{} checked, all ok: {}", e2_all, e1_checked, toys.len(), e1_all);
    // pooled ratios with curve bootstrap, per toy arm group and overall
    js.push_str("\"toy\":{");
    write!(js, "\"T\":{},\"E2_all\":{},\"E1_checked\":{},\"E1_all_ok\":{},\"groups\":{{", T, e2_all, e1_checked, e1_all).unwrap();
    let groups = [("all", None), ("F", Some("F")), ("N2", Some("N2"))];
    for (gi, (gname, gf)) in groups.iter().enumerate() {
        let idx: Vec<usize> = (0..toys.len()).filter(|&i| gf.map_or(true, |f| toys[i].arm == f)).collect();
        let stat = |ix: &[usize], kind: u8, a: usize| -> f64 {
            let (mut na, mut n0) = (0f64, 0f64);
            for &i in ix { let t = &toys[i]; match kind {
                0 => { na += t.pts[a] as f64 / t.v[a].max(1) as f64; n0 += t.pts[0] as f64 / t.v[0] as f64; }
                1 => { na += t.rel[a] as f64; n0 += t.rel[0] as f64; }
                _ => { na += t.dec[a] as f64; n0 += t.dec[0] as f64; } } }
            na / n0
        };
        let mut brng = Rng(SEED ^ 0xB00 ^ gi as u64);
        if gi > 0 { js.push(','); }
        write!(js, "\"{}\":{{\"n_curves\":{},", gname, idx.len()).unwrap();
        let ttot: u64 = idx.iter().map(|&i| toys[i].T).sum();
        let rel: Vec<u64> = (0..5).map(|a| idx.iter().map(|&i| toys[i].rel[a]).sum()).collect();
        let dec: Vec<u64> = (0..5).map(|a| idx.iter().map(|&i| toys[i].dec[a]).sum()).collect();
        let m2ci: Vec<(f64, f64)> = (0..5).map(|a| wilson(rel[a], ttot)).collect();
        write!(js, "\"targets\":{},\"rel\":{:?},\"dec\":{:?},\"M2\":{:?},\"M2_ci95\":{:?},\"D_per_target\":{:?},", ttot, rel, dec,
            rel.iter().map(|&r| r as f64 / ttot as f64).collect::<Vec<_>>(), m2ci.iter().map(|c| [c.0, c.1]).collect::<Vec<_>>(),
            dec.iter().map(|&d| d as f64 / ttot as f64).collect::<Vec<_>>()).unwrap();
        for (kind, kname) in [(0u8, "G_M1"), (1, "ratio_M2"), (2, "ratio_D")] {
            write!(js, "\"{}\":[", kname).unwrap();
            for a in 0..5 {
                let pt = stat(&idx, kind, a);
                let mut bs: Vec<f64> = (0..boots).map(|_| { let ix: Vec<usize> = (0..idx.len()).map(|_| idx[brng.below(idx.len() as u64) as usize]).collect(); stat(&ix, kind, a) }).collect();
                let (l, h) = (pct(&mut bs, 0.025), pct(&mut bs, 0.975));
                if a > 0 { js.push(','); }
                write!(js, "{{\"point\":{},\"ci95\":[{},{}]}}", pt, l, h).unwrap();
                say!("toy {} {} {}: {:.4} [{:.4},{:.4}]", gname, kname, ARMS[a], pt, l, h);
            }
            js.push_str("],");
        }
        js.pop(); js.push('}');
    }
    js.push_str("},\"curves\":[");
    for (i, t) in toys.iter().enumerate() {
        if i > 0 { js.push(','); }
        write!(js, "{{\"arm\":\"{}\",\"bits\":{},\"p\":{},\"N\":{},\"B\":{},\"V\":{:?},\"pts\":{:?},\"rel\":{:?},\"dec\":{:?},\"E1\":{},\"E2\":{}}}",
            t.arm, t.bits, t.p, t.N, t.B, t.v, t.pts, t.rel, t.dec, t.e1.map_or("null".to_string(), |b| b.to_string()), t.e2).unwrap();
    }
    write!(js, "]}},\"wall_s\":{:.1}}}", t0.elapsed().as_secs_f64()).unwrap();
    std::fs::write(format!("{}/results.json", out), js).unwrap();
    say!("done wall={:.1}s", t0.elapsed().as_secs_f64());
}
