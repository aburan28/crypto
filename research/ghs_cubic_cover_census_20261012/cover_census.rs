// Census of GHS-type covers for elliptic curves over K = GF(q^3) (PROTOCOL.md).
// For every F_q-structure tau = A o sigma on P^1_K it computes the (Z/2)^m compositum of the three
// conjugate double covers by branch data: m, r, genus = 2^(m-2)(r-4)+1.
// Dependency-free: rustc -O cover_census.rs ; ./cover_census <q> [sample]
// One JSON object per line on stdout.
use std::collections::HashMap;

const SEED: u64 = 20261012;
struct Rng(u64);
impl Rng {
    fn new(tag: u64) -> Self { Rng(SEED ^ tag.wrapping_mul(0x9E3779B97F4A7C15)) }
    fn next(&mut self) -> u64 { self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15); let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }
    fn below(&mut self, n: u32) -> u32 { (self.next() % n as u64) as u32 }
}

// GF(q^n): 0 is zero, e >= 1 is alpha^(e-1).
struct Fld { q: u32, size: u32, m: u32, zech: Vec<u32>, idx_of: Vec<u32>, el_of: Vec<u32> }
impl Fld {
    fn build(q: u32, n: u32) -> Fld {
        let size = q.pow(n); let m = size - 1; let nn = n as usize;
        let digits = |mut i: u32| -> Vec<u32> { let mut v = vec![0; nn]; for k in 0..nn { v[k] = i % q; i /= q; } v };
        let index = |v: &[u32]| -> u32 { let mut i = 0; for k in (0..nn).rev() { i = i * q + v[k]; } i };
        for cand in 0..size { let md = digits(cand);
            for s in 0..q {
                let mut exp_idx = vec![0u32; m as usize]; let mut cur = vec![0u32; nn]; cur[0] = 1; let mut ok = true;
                for k in 0..m as usize { let ci = index(&cur); if ci == 0 || (k > 0 && ci == 1) { ok = false; break; } exp_idx[k] = ci;
                    let mut nx = vec![0u32; nn]; let top = cur[nn - 1];
                    for i in 0..nn { let sh = if i > 0 { cur[i - 1] } else { 0 }; nx[i] = (sh + (q - md[i]) * top + s * cur[i]) % q; } cur = nx; }
                if !ok || index(&cur) != 1 { continue; }
                let mut el_of = vec![0u32; size as usize]; let mut idx_of = vec![0u32; size as usize];
                for k in 0..m as usize { el_of[exp_idx[k] as usize] = k as u32 + 1; idx_of[k + 1] = exp_idx[k]; }
                let mut zech = vec![0u32; m as usize];
                for k in 0..m as usize { let mut dv = digits(exp_idx[k]); dv[0] = (dv[0] + 1) % q; zech[k] = el_of[index(&dv) as usize]; }
                return Fld { q, size, m, zech, idx_of, el_of };
            } }
        panic!("no field");
    }
    #[inline] fn mul(&self, a: u32, b: u32) -> u32 { if a == 0 || b == 0 { 0 } else { let s = a + b - 2; (if s >= self.m { s - self.m } else { s }) + 1 } }
    #[inline] fn inv(&self, a: u32) -> u32 { (self.m - (a - 1)) % self.m + 1 }
    #[inline] fn neg(&self, a: u32) -> u32 { if a == 0 { 0 } else { (a - 1 + self.m / 2) % self.m + 1 } }
    #[inline] fn add(&self, a: u32, b: u32) -> u32 { if a == 0 { return b; } if b == 0 { return a; }
        let dlt = if b >= a { b - a } else { b + self.m - a }; let z = self.zech[dlt as usize]; if z == 0 { 0 } else { let s = a + z - 2; (if s >= self.m { s - self.m } else { s }) + 1 } }
    #[inline] fn sub(&self, a: u32, b: u32) -> u32 { self.add(a, self.neg(b)) }
    #[inline] fn frob(&self, a: u32) -> u32 { if a == 0 { 0 } else { (((a - 1) as u64 * self.q as u64) % self.m as u64) as u32 + 1 } }
    #[inline] fn is_sq(&self, a: u32) -> bool { a == 0 || (a - 1) % 2 == 0 }
    fn int(&self, c: i64) -> u32 { self.el_of[c.rem_euclid(self.q as i64) as usize] }
    fn pow(&self, a: u32, e: u64) -> u32 { if a == 0 { 0 } else { (((a - 1) as u64 * (e % self.m as u64)) % self.m as u64) as u32 + 1 } }
}

// Small polynomials on the stack, low -> high, n = length (0 for the zero polynomial).
#[derive(Clone, Copy, PartialEq)]
struct P { c: [u32; 14], n: usize }
fn pz() -> P { P { c: [0; 14], n: 0 } }
fn ptrim(mut a: P) -> P { while a.n > 0 && a.c[a.n - 1] == 0 { a.n -= 1; } a }
fn pmul(f: &Fld, a: &P, b: &P) -> P { let mut r = pz(); if a.n == 0 || b.n == 0 { return r; } r.n = a.n + b.n - 1;
    for i in 0..a.n { if a.c[i] == 0 { continue; } for j in 0..b.n { r.c[i + j] = f.add(r.c[i + j], f.mul(a.c[i], b.c[j])); } } r }
fn pdivrem(f: &Fld, a: &P, b: &P) -> (P, P) { let mut r = *a; let mut qv = pz(); if r.n < b.n { return (qv, r); } qv.n = r.n - b.n + 1; let li = f.inv(b.c[b.n - 1]);
    while r.n >= b.n { let k = r.n - b.n; let c = f.mul(r.c[r.n - 1], li); qv.c[k] = c; for j in 0..b.n { r.c[k + j] = f.sub(r.c[k + j], f.mul(c, b.c[j])); } r = ptrim(r); if r.n == 0 { break; } }
    (qv, r) }
fn pmonic(f: &Fld, mut a: P) -> P { if a.n > 0 { let li = f.inv(a.c[a.n - 1]); for i in 0..a.n { a.c[i] = f.mul(a.c[i], li); } } a }
fn pgcd(f: &Fld, a: &P, b: &P) -> P { let (mut x, mut y) = (*a, *b); while y.n > 0 { let r = pdivrem(f, &x, &y).1; x = y; y = r; } pmonic(f, x) }

// A finite subset of P^1: monic squarefree affine polynomial plus a flag for infinity.
#[derive(Clone, Copy, PartialEq)]
struct Set { p: P, inf: bool }
fn size(s: &Set) -> usize { s.p.n.max(1) - 1 + s.inf as usize }
fn sdiff(f: &Fld, a: &Set, b: &Set) -> Set { let g = pgcd(f, &a.p, &b.p); Set { p: pmul(f, &pdivrem(f, &a.p, &g).0, &pdivrem(f, &b.p, &g).0), inf: a.inf ^ b.inf } }
fn union(f: &Fld, a: &Set, b: &Set) -> Set { let g = pgcd(f, &a.p, &b.p); Set { p: pmul(f, &a.p, &pdivrem(f, &b.p, &g).0), inf: a.inf | b.inf } }
// quartic binary form F = sum c[i] X^i Z^(4-i)  ->  its zero set
fn set_of(f: &Fld, c: &[u32; 5]) -> Set { let mut p = pz(); p.n = 5; p.c[..5].copy_from_slice(c); let p = ptrim(p); Set { inf: p.n < 5, p: pmonic(f, p) } }
type Mat = [u32; 4]; // [a, b, c, d] acting on points as x -> (a x + b)/(c x + d)
// tau(F) for tau = A o sigma: conjugate the coefficients, then move the zeros by A:  G(X, Z) = F^sigma(dX - bZ, -cX + aZ)
fn tau(f: &Fld, c: &[u32; 5], a: &Mat) -> [u32; 5] {
    let l1 = [f.neg(a[1]), a[3]]; let l2 = [a[0], f.neg(a[2])];
    let mut p1 = [[0u32; 5]; 5]; let mut p2 = [[0u32; 5]; 5]; p1[0][0] = 1; p2[0][0] = 1;
    for i in 1..5 { for j in 0..i { p1[i][j] = f.add(p1[i][j], f.mul(p1[i - 1][j], l1[0])); p1[i][j + 1] = f.add(p1[i][j + 1], f.mul(p1[i - 1][j], l1[1]));
        p2[i][j] = f.add(p2[i][j], f.mul(p2[i - 1][j], l2[0])); p2[i][j + 1] = f.add(p2[i][j + 1], f.mul(p2[i - 1][j], l2[1])); } }
    let mut g = [0u32; 5];
    for i in 0..5 { let ci = f.frob(c[i]); if ci == 0 { continue; } for u in 0..=i { if p1[i][u] == 0 { continue; } let t = f.mul(ci, p1[i][u]); for v in 0..=(4 - i) { g[u + v] = f.add(g[u + v], f.mul(t, p2[4 - i][v])); } } }
    g
}
fn mmul(f: &Fld, a: &Mat, b: &Mat) -> Mat { [f.add(f.mul(a[0], b[0]), f.mul(a[1], b[2])), f.add(f.mul(a[0], b[1]), f.mul(a[1], b[3])), f.add(f.mul(a[2], b[0]), f.mul(a[3], b[2])), f.add(f.mul(a[2], b[1]), f.mul(a[3], b[3]))] }
fn cocycles(f: &Fld, threads: usize) -> Vec<Mat> { // all A in PGL2(K) with A A^sigma A^sigma^2 scalar, normalised
    let test = |a: &Mat| -> bool { let s1 = [f.frob(a[0]), f.frob(a[1]), f.frob(a[2]), f.frob(a[3])]; let s2 = [f.frob(s1[0]), f.frob(s1[1]), f.frob(s1[2]), f.frob(s1[3])];
        let n = mmul(f, &mmul(f, a, &s1), &s2); n[1] == 0 && n[2] == 0 && n[0] == n[3] && n[0] != 0 };
    let mut out: Vec<Mat> = vec![];
    std::thread::scope(|sc| { let hs: Vec<_> = (0..threads).map(|t| { let test = &test; sc.spawn(move || { let mut v = vec![];
        for b in (t as u32..f.size).step_by(threads) { for c in 0..f.size { for d in 0..f.size { if d == f.mul(b, c) { continue; } let a = [1, b, c, d]; if test(&a) { v.push(a); } } } } v }) }).collect();
        for h in hs { out.extend(h.join().unwrap()); } });
    for c in 1..f.size { for d in 0..f.size { let a = [0, 1, c, d]; if test(&a) { out.push(a); } } }
    out.sort(); out
}

#[derive(Clone, Default)]
struct Res { by_genus: [u64; 18], m1: u64, min_g: u32, types: HashMap<(u32, u32, i32), u64>, mismatch: u64 }
// census of one curve y^2 = x^3 + A x + B over all structures
fn census(f: &Fld, a: u32, b: u32, cos: &[Mat]) -> Res {
    let b0f = [b, a, 0, 1, 0]; let s0 = set_of(f, &b0f); let mut res = Res { min_g: 99, ..Default::default() };
    for mat in cos { let b1f = tau(f, &b0f, mat); let b2f = tau(f, &b1f, mat); let (s1, s2) = (set_of(f, &b1f), set_of(f, &b2f));
        if s1 == s0 { res.m1 += 1; continue; }
        let d01 = sdiff(f, &s0, &s1); let d012 = sdiff(f, &d01, &s2); let m = if size(&d012) == 0 { 2 } else { 3 };
        let u = union(f, &union(f, &s0, &s1), &s2); let r = size(&u) as u32; let g = (if m == 2 { 1 } else { 2 }) * (r - 4) + 1;
        res.by_genus[g as usize] += 1;
        if g <= 5 || g <= res.min_g { // subcover-sum cross-check and hyperelliptic test
            let d02 = sdiff(f, &s0, &s2); let d12 = sdiff(f, &s1, &s2);
            let gq = |s: &Set| (size(s) as i32 - 2) / 2;
            let sum = if m == 2 { gq(&s0) + gq(&s1) + gq(&s2) } else { gq(&s0) + gq(&s1) + gq(&s2) + gq(&d01) + gq(&d02) + gq(&d12) + gq(&d012) };
            if sum != g as i32 { res.mismatch += 1; }
            let hyp = if m == 2 { -1 } else { let lines = [[&s0, &s1, &d01], [&s0, &s2, &d02], [&s1, &s2, &d12], [&s0, &d12, &d012], [&s1, &d02, &d012], [&s2, &d01, &d012], [&d01, &d02, &d12]];
                lines.iter().any(|l| size(&union(f, &union(f, l[0], l[1]), l[2])) == 3) as i32 };
            if g < res.min_g { res.min_g = g; res.types.clear(); }
            if g == res.min_g { *res.types.entry((m, r, hyp)).or_insert(0) += 1; }
        }
    }
    res
}
fn order(f: &Fld, a: u32, b: u32) -> u64 { let mut c = 1u64; for x in 0..f.size { let rhs = f.add(f.add(f.mul(x, f.mul(x, x)), f.mul(a, x)), b); c += if rhs == 0 { 1 } else if f.is_sq(rhs) { 2 } else { 0 }; } c }
fn nonsingular(f: &Fld, a: u32, b: u32) -> bool { f.add(f.mul(f.int(4), f.mul(a, f.mul(a, a))), f.mul(f.int(27), f.mul(b, b))) != 0 }
// K-rational 2-isogeny neighbours: kernel (e, 0) gives A' = -15 e^2 - 4A, B' = -22 e^3 - 8 A e
fn neighbours(f: &Fld, a: u32, b: u32) -> Vec<(u32, u32)> { let mut v = vec![]; for e in 0..f.size { let e2 = f.mul(e, e);
    if f.add(f.add(f.mul(e, e2), f.mul(a, e)), b) == 0 { v.push((f.sub(f.mul(f.int(-15), e2), f.mul(f.int(4), a)), f.sub(f.mul(f.int(-22), f.mul(e, e2)), f.mul(f.int(8), f.mul(a, e))))); } } v }
fn maxprime(mut n: u64) -> u64 { let mut best = 1; let mut p = 2; while p * p <= n { while n % p == 0 { n /= p; best = p; } p += 1; } if n > 1 { n } else { best } }
fn types_json(r: &Res) -> String { let mut v: Vec<_> = r.types.iter().collect(); v.sort();
    format!("[{}]", v.iter().map(|((m, rr, h), c)| format!("{{\"m\":{},\"r\":{},\"hyperelliptic\":{},\"structures\":{}}}", m, rr, match h { 1 => "true", 0 => "false", _ => "null" }, c)).collect::<Vec<_>>().join(",")) }
fn hist_json(h: &HashMap<u32, u64>) -> String { let mut v: Vec<_> = h.iter().collect(); v.sort(); format!("{{{}}}", v.iter().map(|(g, c)| format!("\"{}\":{}", g, c)).collect::<Vec<_>>().join(",")) }
fn par_census(f: &Fld, curves: &[(u32, u32)], cos: &[Mat], threads: usize) -> Vec<Res> {
    let mut out: Vec<Res> = vec![Res::default(); curves.len()];
    std::thread::scope(|sc| { let hs: Vec<_> = (0..threads).map(|t| sc.spawn(move || { let mut v = vec![]; for i in (t..curves.len()).step_by(threads) { v.push((i, census(f, curves[i].0, curves[i].1, cos))); } v })).collect();
        for h in hs { for (i, r) in h.join().unwrap() { out[i] = r; } } });
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect(); let q: u32 = args[1].parse().unwrap(); let sample: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(0);
    let threads: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(8);
    let f = Fld::build(q, 3); let cos = cocycles(&f, threads); let qq = f.size as u64; let expect = (qq * qq * qq - qq) / (q as u64 * q as u64 * q as u64 - q as u64);
    println!("{{\"rec\":\"field\",\"q\":{},\"n\":3,\"structures\":{},\"structures_expected\":{},\"mode\":\"{}\"}}", q, cos.len(), expect, if sample == 0 { "exhaustive" } else { "sample" });
    let emit = |a: u32, b: u32, n: u64, r: &Res, g2: u32, gc: Option<u32>, role: &str| {
        println!("{{\"rec\":\"curve\",\"q\":{},\"role\":\"{}\",\"A_idx\":{},\"B_idx\":{},\"order\":{},\"maxprime\":{},\"min_genus\":{},\"types_at_min\":{},\"structures_by_genus\":{},\"m1_structures\":{},\"mismatch\":{},\"min_genus_2isog\":{},\"min_genus_class\":{}}}",
            q, role, f.idx_of[a as usize], f.idx_of[b as usize], n, maxprime(n), r.min_g, types_json(r),
            format!("{{{}}}", (0..18).filter(|&g| r.by_genus[g] > 0).map(|g| format!("\"{}\":{}", g, r.by_genus[g])).collect::<Vec<_>>().join(",")), r.m1, r.mismatch, g2, gc.map(|x| x.to_string()).unwrap_or("null".into())); };
    if sample == 0 {
        // every K-isomorphism class: orbits of (A, B) -> (u^4 A, u^6 B)
        let sz = f.size as usize; let mut class = vec![u32::MAX; sz * sz]; let mut reps = vec![];
        for a in 0..f.size { for b in 0..f.size { if !nonsingular(&f, a, b) || class[a as usize * sz + b as usize] != u32::MAX { continue; }
            let id = reps.len() as u32; reps.push((a, b)); for u in 1..f.size { let (a2, b2) = (f.mul(a, f.pow(u, 4)), f.mul(b, f.pow(u, 6))); class[a2 as usize * sz + b2 as usize] = id; } } }
        let res = par_census(&f, &reps, &cos, threads); let ords: Vec<u64> = reps.iter().map(|&(a, b)| order(&f, a, b)).collect();
        let mut cmin: HashMap<u64, u32> = HashMap::new(); for i in 0..reps.len() { let e = cmin.entry(ords[i]).or_insert(99); *e = (*e).min(res[i].min_g); }
        let mut h: HashMap<(&str, u64), HashMap<u32, u64>> = HashMap::new();
        for i in 0..reps.len() { let (a, b) = reps[i]; let mut g2 = res[i].min_g; for (a2, b2) in neighbours(&f, a, b) { g2 = g2.min(res[class[a2 as usize * sz + b2 as usize] as usize].min_g); }
            let gc = cmin[&ords[i]]; emit(a, b, ords[i], &res[i], g2, Some(gc), "class");
            for (k, g) in [("self", res[i].min_g), ("2isog", g2), ("class", gc)] { *h.entry((k, ords[i] % 2)).or_default().entry(g).or_insert(0) += 1; } }
        for ((k, par), hh) in { let mut v: Vec<_> = h.into_iter().collect(); v.sort_by_key(|x| x.0); v } { println!("{{\"rec\":\"summary\",\"q\":{},\"scope\":\"{}\",\"order_parity\":\"{}\",\"curves_by_min_genus\":{}}}", q, k, if par == 0 { "even" } else { "odd" }, hist_json(&hh)); }
        println!("{{\"rec\":\"totals\",\"q\":{},\"classes\":{},\"isogeny_classes\":{},\"mismatches\":{}}}", q, reps.len(), cmin.len(), res.iter().map(|r| r.mismatch).sum::<u64>());
    } else {
        let mut r = Rng::new(q as u64); let mut picks = vec![]; while picks.len() < sample { let (a, b) = (r.below(f.size), r.below(f.size)); if nonsingular(&f, a, b) { picks.push((a, b)); } }
        let mut all = picks.clone(); let nb: Vec<Vec<(u32, u32)>> = picks.iter().map(|&(a, b)| neighbours(&f, a, b)).collect(); for v in &nb { all.extend(v.iter().cloned()); }
        let res = par_census(&f, &all, &cos, threads); let mut k = picks.len(); let mut h: HashMap<(&str, u64), HashMap<u32, u64>> = HashMap::new(); let mut mm = 0;
        for i in 0..picks.len() { let (a, b) = picks[i]; let n = order(&f, a, b); let mut g2 = res[i].min_g; for _ in &nb[i] { g2 = g2.min(res[k].min_g); mm += res[k].mismatch; k += 1; } mm += res[i].mismatch;
            emit(a, b, n, &res[i], g2, None, "sample");
            for (key, g) in [("self", res[i].min_g), ("2isog", g2)] { *h.entry((key, n % 2)).or_default().entry(g).or_insert(0) += 1; } }
        for ((k, par), hh) in { let mut v: Vec<_> = h.into_iter().collect(); v.sort_by_key(|x| x.0); v } { println!("{{\"rec\":\"summary\",\"q\":{},\"scope\":\"{}\",\"order_parity\":\"{}\",\"curves_by_min_genus\":{}}}", q, k, if par == 0 { "even" } else { "odd" }, hist_json(&hh)); }
        println!("{{\"rec\":\"totals\",\"q\":{},\"sampled\":{},\"mismatches\":{}}}", q, picks.len(), mm);
    }
}
