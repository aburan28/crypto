// 13-xy-windows: joint x/y window factor bases vs x-only windows of matched size.
// Usage: xy_windows <out_dir> <T_targets> <n_toy_F> <n_toy_N2> <boots> <toy_budget_s>
// All results are computed here; Python only plots and lays out the PDF.
#![allow(non_snake_case, dead_code)]
use std::collections::HashMap;
use std::fmt::Write as FmtWrite;
use std::io::Write;
use std::time::Instant;

include!("common.rs");

const SEED: u64 = 46013;
const ARMS: [&str; 11] = ["X", "XY_sign", "XY_par", "XY_box2", "XY_box4", "XY_box6", "XY_bits4", "Xbits4", "R_random", "PC_logAP", "X_half"];
const NA: usize = 11;

fn splitmix(x: u64) -> u64 { let mut z = x.wrapping_add(0x9E3779B97F4A7C15); z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }

fn wilson(k: u64, n: u64) -> (f64, f64) {
    let z = 1.959963984540054; let n_ = n as f64; let ph = k as f64 / n_;
    let d = 1.0 + z * z / n_; let c = (ph + z * z / (2.0 * n_)) / d;
    let h = z * ((ph * (1.0 - ph) / n_ + z * z / (4.0 * n_ * n_)).sqrt()) / d;
    (c - h, c + h)
}
fn pct(v: &mut Vec<f64>, q: f64) -> f64 { v.sort_by(|a, b| a.partial_cmp(b).unwrap()); let i = ((v.len() - 1) as f64 * q).round() as usize; v[i] }

/// Scan-arm predicate (the single membership instrument). Inputs: x, one root y, p.
#[inline]
fn scan_member(arm: usize, x: u64, y: u64, p: u64) -> bool {
    let ny = p - y; let cy = y.min(ny);
    match arm {
        0 => true,
        1 => y <= (p - 1) / 2 || ny <= (p - 1) / 2,
        2 => y % 2 == 0 || ny % 2 == 0,
        3 => cy < ((p + 1) / 2) >> 2,
        4 => cy < ((p + 1) / 2) >> 4,
        5 => cy < ((p + 1) / 2) >> 6,
        6 => (cy >> 1) & 15 == 0,
        7 => x % 16 == 0,
        _ => unreachable!(),
    }
}

struct Toy { arm: String, bits: u32, p: u64, N: u64, K: u64, scan: [u64; NA], rel: [u64; NA], dec: [u64; NA], x2_sets: bool, x2_targets: bool, T: u64 }

fn toy_measure(arm_name: &str, E: &Curve, N: u64, T: u64, rng: &mut Rng) -> Toy {
    let p = E.p;
    let K = ((N as f64 / 2.0).sqrt()).ceil() as u64;
    let mut sets: Vec<Vec<(u64, u64)>> = vec![vec![]; NA];
    let mut scan = [0u64; NA];
    // scan arms 0..8
    let mut x = 0u64;
    loop {
        let mut all = true;
        for a in 0..8 { if (sets[a].len() as u64) < K { all = false; } }
        if all || x >= p { break; }
        let r = E.rhs(x);
        if r != 0 { if let Some(y) = sqrtm(r, p) {
            for a in 0..8 { if (sets[a].len() as u64) < K && scan_member(a, x, y, p) { sets[a].push((x, y)); if sets[a].len() as u64 == K { scan[a] = x + 1; } } }
        } }
        x += 1;
    }
    for a in 0..8 { assert!(sets[a].len() as u64 == K, "scan arm {} under-filled", a); }
    // R: random liftable x
    { let mut seen = std::collections::HashSet::new(); let mut tries = 0u64;
      while (sets[8].len() as u64) < K { let x = rng.below(p); tries += 1; let r = E.rhs(x); if r == 0 || seen.contains(&x) { continue; }
        if let Some(y) = sqrtm(r, p) { seen.insert(x); sets[8].push((x, y)); } }
      scan[8] = tries; }
    let G = Some(E.rand_point(rng));
    // PC: jG, j=1..K
    { let mut Pj = G; for _ in 0..K { let (x, y) = Pj.unwrap(); sets[9].push((x, y)); Pj = E.add(Pj, G); } scan[9] = K; }
    // X_half
    let kh = (K + 1) / 2; sets[10] = sets[0][..kh as usize].to_vec(); scan[10] = sets[0][kh as usize - 1].0 + 1;
    let x2_sets = sets[1] == sets[0] && sets[2] == sets[0];
    let mut mask: HashMap<u64, u16> = HashMap::new();
    let mut ys: HashMap<u64, u64> = HashMap::new();
    for a in 0..NA { for &(x, y) in &sets[a] { *mask.entry(x).or_insert(0) |= 1 << a; ys.entry(x).or_insert(y); } }
    let uni: Vec<(u64, u64, u16)> = mask.iter().map(|(&x, &m)| (x, ys[&x], m)).collect();
    let (mut rel, mut dec) = ([0u64; NA], [0u64; NA]);
    let mut x2_targets = true;
    for _ in 0..T {
        let k = 1 + rng.below(N - 1);
        let R = E.mul(k, G);
        let rx = R.unwrap().0;
        let mut c = [0u64; NA];
        for &(x1, y1, m1) in &uni {
            if x1 == rx { continue; }
            let P = Some((x1, y1));
            for Q in [E.add(R, E.neg(P)), E.add(R, P)] {
                if let Some((x2, _)) = Q {
                    if x2 >= x1 { if let Some(&m2) = mask.get(&x2) { let mm = m1 & m2; if mm != 0 { for a in 0..NA { if mm >> a & 1 == 1 { c[a] += 1; } } } } }
                }
            }
        }
        if c[1] != c[0] || c[2] != c[0] { x2_targets = false; }
        for a in 0..NA { dec[a] += c[a]; if c[a] > 0 { rel[a] += 1; } }
    }
    Toy { arm: arm_name.to_string(), bits: 64 - p.leading_zeros(), p, N, K, scan, rel, dec, x2_sets, x2_targets, T }
}

/// Dispersion of box counts in (x, cy) space. cys[x] = canonical |y| or u64::MAX.
fn dispersion(cys: &[u64], p: u64, B1: u64, plant: bool, boots: usize, seed: u64) -> (f64, f64, f64, f64, u64, u64) {
    let H = 16 * p / B1; // lambda ~ 16 (B1 * H / p)
    let half = (p + 1) / 2; let nb = half / H; // full bands only
    let nblk = cys.len() as u64 / B1;
    let mut cnt = vec![0u32; (nblk * nb) as usize];
    for (x, &cy) in cys.iter().enumerate() {
        if cy == u64::MAX { continue; }
        let band = cy / H; if band >= nb { continue; }
        let blk = x as u64 / B1; cnt[(blk * nb + band) as usize] += 1;
    }
    if plant { for i in 0..cnt.len() { if splitmix(i as u64 ^ 0x9A17_0013) % 10 == 0 { cnt[i] += 8; } } }
    let mut s1 = vec![0f64; nblk as usize]; let mut s2 = vec![0f64; nblk as usize];
    for b in 0..nblk as usize { for j in 0..nb as usize { let c = cnt[b * nb as usize + j] as f64; s1[b] += c; s2[b] += c * c; } }
    let idx = |sel: &mut dyn Iterator<Item = usize>| -> f64 {
        let (mut a, mut q, mut n) = (0f64, 0f64, 0f64);
        for b in sel { a += s1[b]; q += s2[b]; n += nb as f64; }
        let m = a / n; let v = (q / n - m * m) * n / (n - 1.0); v / m
    };
    let pt = idx(&mut (0..nblk as usize));
    let mut rng = Rng(seed);
    let mut bs: Vec<f64> = (0..boots).map(|_| { let v: Vec<usize> = (0..nblk).map(|_| rng.below(nblk) as usize).collect(); idx(&mut v.into_iter()) }).collect();
    let mean = s1.iter().sum::<f64>() / (nblk * nb) as f64;
    (pt, pct(&mut bs, 0.025), pct(&mut bs, 0.975), mean, nblk, nb)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out = &args[1]; let T: u64 = args[2].parse().unwrap();
    let nF: usize = args[3].parse().unwrap(); let nN2: usize = args[4].parse().unwrap();
    let boots: usize = args[5].parse().unwrap(); let budget: f64 = args[6].parse().unwrap();
    let t0 = Instant::now();
    let mut rng = Rng(SEED);
    let mut log = std::fs::File::create(format!("{}/run.log", out)).unwrap();
    macro_rules! say { ($($a:tt)*) => {{ let s = format!($($a)*); println!("{}", s); writeln!(log, "{}", s).unwrap(); log.flush().unwrap(); }} }
    say!("xy_windows seed={} T={} nF={} nN2={} boots={} toy_budget_s={}", SEED, T, nF, nN2, boots, budget);
    let mut js = String::new();
    write!(js, "{{\"item\":\"13-xy-windows\",\"seed\":{},\"arms\":{:?},", SEED, ARMS).unwrap();

    // ---------------- pair 46 ----------------
    let p46: u64 = 5947698030651343; let t46: u64 = 34489701; let N46: u64 = 5947697996161643;
    assert!(p46 == 5 * t46 * t46 + 19 * t46 + 19 && N46 == p46 + 1 - t46 && is_prime(p46) && is_prime(N46));
    let mut curves: Vec<(String, Curve, bool)> = vec![
        ("crater".into(), Curve { p: p46, a: 5757116383569871, b: 3560173328417442 }, true),
        ("floor".into(), Curve { p: p46, a: 2087370704375809, b: 1402887937802154 }, true)];
    let mut crng = Rng(SEED ^ 0xC0FFEE);
    for i in 0..4 { loop { let E = Curve { p: p46, a: crng.below(p46), b: crng.below(p46) }; if E.disc_ok() { curves.push((format!("random{}", i), E, false)); break; } } }
    let L: usize = 1 << 23;
    let shapes = [1u64 << 8, 1 << 12, 1 << 16];
    js.push_str(&format!("\"pair46\":{{\"p\":{},\"N\":{},\"L\":{},\"shapes_B1\":{:?},\"curves\":[", p46, N46, L, shapes));
    for (ci, (name, E, known)) in curves.iter().enumerate() {
        let ok = if *known { let o = order_is(E, N46, &mut rng); assert!(o); o } else { true };
        let tc = Instant::now();
        let mut cys = vec![u64::MAX; L];
        let mut nl = 0u64; let mut r0 = 0u64;
        for x in 0..L as u64 { let r = E.rhs(x); if r == 0 { r0 += 1; continue; } if let Some(y) = sqrtm(r, p46) { cys[x as usize] = y.min(p46 - y); nl += 1; } }
        say!("{}: order_check={} liftable={}/{} rhs0={} ({:.1}s)", name, ok, nl, L, r0, tc.elapsed().as_secs_f64());
        // corner densities
        let mut corner = vec![];
        for k in 0..=16u32 { let th = ((p46 + 1) / 2) >> k; let c = cys.iter().filter(|&&c| c != u64::MAX && c < th).count() as u64;
            let (lo, hi) = wilson(c, nl); let e = 1.0 / (1u64 << k) as f64; corner.push((k, c, c as f64 / nl as f64 / e, lo / e, hi / e)); }
        say!("{} corner obs/exp k=0..16: {:?}", name, corner.iter().map(|c| format!("{}:{:.3}", c.0, c.2)).collect::<Vec<_>>());
        let mut disp = vec![];
        for (si, &B1) in shapes.iter().enumerate() {
            let d = dispersion(&cys, p46, B1, false, boots, SEED ^ (ci as u64 * 16 + si as u64));
            say!("{} B1=2^{} I={:.4} [{:.4},{:.4}] mean={:.3} blocks={} bands={}", name, B1.trailing_zeros(), d.0, d.1, d.2, d.3, d.4, d.5);
            let pl = if ci == 0 { let q = dispersion(&cys, p46, B1, true, boots, SEED ^ 0x77 ^ si as u64); say!("{} PLANTED B1=2^{} I={:.4} [{:.4},{:.4}]", name, B1.trailing_zeros(), q.0, q.1, q.2); Some(q) } else { None };
            disp.push((B1, d, pl));
        }
        if ci > 0 { js.push(','); }
        write!(js, "{{\"name\":\"{}\",\"a\":{},\"b\":{},\"pair46_curve\":{},\"order_check\":{},\"liftable\":{},\"rhs0\":{},\"corner\":[", name, E.a, E.b, known, ok, nl, r0).unwrap();
        for (i, c) in corner.iter().enumerate() { if i > 0 { js.push(','); } write!(js, "{{\"k\":{},\"count\":{},\"obs_over_exp\":{},\"ci95\":[{},{}]}}", c.0, c.1, c.2, c.3, c.4).unwrap(); }
        js.push_str("],\"dispersion\":[");
        for (i, (B1, d, pl)) in disp.iter().enumerate() { if i > 0 { js.push(','); }
            write!(js, "{{\"B1\":{},\"H\":{},\"index\":{},\"ci95\":[{},{}],\"mean\":{},\"blocks\":{},\"bands\":{}", B1, 16 * p46 / B1, d.0, d.1, d.2, d.3, d.4, d.5).unwrap();
            if let Some(q) = pl { write!(js, ",\"planted\":{{\"index\":{},\"ci95\":[{},{}],\"mean\":{}}}", q.0, q.1, q.2, q.3).unwrap(); }
            js.push('}'); }
        js.push_str("]}");
    }
    js.push_str("]},");
    let t_p46 = t0.elapsed().as_secs_f64();
    say!("pair46 done at {:.1}s", t_p46);

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
    let tt0 = Instant::now();
    let mut cut_at: Option<usize> = None;
    for (i, (arm, p, N, a, b)) in sel.iter().enumerate() {
        if tt0.elapsed().as_secs_f64() > budget { cut_at = Some(i); say!("BUDGET GUARD: stopping toy loop at curve index {}", i); break; }
        let E = Curve { p: *p, a: *a, b: *b };
        let mut trng = Rng(SEED ^ (i as u64));
        assert!(order_is(&E, *N, &mut trng), "toy order check failed");
        let t = toy_measure(arm, &E, *N, T, &mut trng);
        if i % 20 == 0 { say!("toy {} {} bits={} K={} scan={:?} dec={:?} rel={:?} x2={}/{} t={:.1}s", i, arm, t.bits, t.K, t.scan, t.dec, t.rel, t.x2_sets, t.x2_targets, t0.elapsed().as_secs_f64()); }
        toys.push(t);
    }
    let x2_all = toys.iter().all(|t| t.x2_sets && t.x2_targets);
    say!("toy X2 (XY_sign, XY_par == X: sets and every target): {} over {} curves", x2_all, toys.len());
    write!(js, "\"toy\":{{\"T\":{},\"n_selected\":{},\"n_done\":{},\"cut_at\":{},\"X2_all\":{},\"groups\":{{", T, sel.len(), toys.len(), cut_at.map_or("null".into(), |c| c.to_string()), x2_all).unwrap();
    let groups = [("all", None), ("F", Some("F")), ("N2", Some("N2"))];
    for (gi, (gname, gf)) in groups.iter().enumerate() {
        let idx: Vec<usize> = (0..toys.len()).filter(|&i| gf.map_or(true, |f| toys[i].arm == f)).collect();
        let stat = |ix: &[usize], kind: u8, a: usize| -> f64 {
            let (mut na, mut n0) = (0f64, 0f64);
            for &i in ix { let t = &toys[i]; match kind {
                0 => { na += t.rel[a] as f64; n0 += t.rel[0] as f64; }
                1 => { na += t.dec[a] as f64; n0 += t.dec[0] as f64; }
                _ => { na += t.scan[a] as f64 / t.K as f64; n0 += 1.0; } } }
            na / n0
        };
        let mut brng = Rng(SEED ^ 0xB00 ^ gi as u64);
        if gi > 0 { js.push(','); }
        write!(js, "\"{}\":{{\"n_curves\":{},", gname, idx.len()).unwrap();
        let ttot: u64 = idx.iter().map(|&i| toys[i].T).sum();
        let rel: Vec<u64> = (0..NA).map(|a| idx.iter().map(|&i| toys[i].rel[a]).sum()).collect();
        let dec: Vec<u64> = (0..NA).map(|a| idx.iter().map(|&i| toys[i].dec[a]).sum()).collect();
        let m2ci: Vec<(f64, f64)> = (0..NA).map(|a| wilson(rel[a], ttot)).collect();
        write!(js, "\"targets\":{},\"rel\":{:?},\"dec\":{:?},\"M2\":{:?},\"M2_ci95\":{:?},\"D_per_target\":{:?},", ttot, rel, dec,
            rel.iter().map(|&r| r as f64 / ttot as f64).collect::<Vec<_>>(), m2ci.iter().map(|c| [c.0, c.1]).collect::<Vec<_>>(),
            dec.iter().map(|&d| d as f64 / ttot as f64).collect::<Vec<_>>()).unwrap();
        for (kind, kname) in [(0u8, "ratio_M2"), (1, "ratio_D"), (2, "scan_per_element")] {
            write!(js, "\"{}\":[", kname).unwrap();
            for a in 0..NA {
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
        write!(js, "{{\"arm\":\"{}\",\"bits\":{},\"p\":{},\"N\":{},\"K\":{},\"scan\":{:?},\"rel\":{:?},\"dec\":{:?},\"X2_sets\":{},\"X2_targets\":{}}}",
            t.arm, t.bits, t.p, t.N, t.K, t.scan, t.rel, t.dec, t.x2_sets, t.x2_targets).unwrap();
    }
    write!(js, "]}},\"wall_pair46_s\":{:.1},\"wall_s\":{:.1}}}", t_p46, t0.elapsed().as_secs_f64()).unwrap();
    std::fs::write(format!("{}/results.json", out), js).unwrap();
    say!("done wall={:.1}s", t0.elapsed().as_secs_f64());
}
