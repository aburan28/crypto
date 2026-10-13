// Build the genus-3 plane quartic over F_q attached to E0 : y0^2 = u u'' over K = GF(q^3) (PROTOCOL.md).
// Dependency-free: rustc --edition 2021 -O genus3_build.rs ; ./genus3_build <q> <split cases> <irreducible cases>
// One JSON object per line on stdout.
use std::collections::BTreeMap;

// GF(q^n): 0 is zero, e >= 1 is alpha^(e-1).
#[allow(dead_code)]
struct Fld { q: u32, n: u32, size: u32, m: u32, zech: Vec<u32>, idx_of: Vec<u32>, el_of: Vec<u32> }
impl Fld {
    fn build(q: u32, n: u32) -> Fld {
        let size = q.pow(n); let m = size - 1; let nn = n as usize;
        let digits = |mut i: u32| -> Vec<u32> { let mut v = vec![0; nn]; for k in 0..nn { v[k] = i % q; i /= q; } v };
        let index = |v: &[u32]| -> u32 { let mut i = 0; for k in (0..nn).rev() { i = i * q + v[k]; } i };
        for cand in 0..size { let md = digits(cand);
            for s in 0..q { if n == 1 && (s == 0 || cand != 0) { continue; }
                let mut exp_idx = vec![0u32; m as usize]; let mut cur = vec![0u32; nn]; cur[0] = 1; let mut ok = true;
                for k in 0..m as usize { let ci = index(&cur); if ci == 0 || (k > 0 && ci == 1) { ok = false; break; } exp_idx[k] = ci;
                    let mut nx = vec![0u32; nn];
                    if n == 1 { nx[0] = cur[0] * s % q; } else { let top = cur[nn - 1]; for i in 0..nn { let sh = if i > 0 { cur[i - 1] } else { 0 }; nx[i] = (sh + (q - md[i]) * top + s * cur[i]) % q; } }
                    cur = nx; }
                if !ok || index(&cur) != 1 { continue; }
                let mut el_of = vec![0u32; size as usize]; let mut idx_of = vec![0u32; size as usize];
                for k in 0..m as usize { el_of[exp_idx[k] as usize] = k as u32 + 1; idx_of[k + 1] = exp_idx[k]; }
                let mut zech = vec![0u32; m as usize];
                for k in 0..m as usize { let mut dv = digits(exp_idx[k]); dv[0] = (dv[0] + 1) % q; zech[k] = el_of[index(&dv) as usize]; }
                return Fld { q, n, size, m, zech, idx_of, el_of };
            } }
        panic!("no field");
    }
    #[inline] fn mul(&self, a: u32, b: u32) -> u32 { if a == 0 || b == 0 { 0 } else { (a - 1 + b - 1) % self.m + 1 } }
    #[inline] fn inv(&self, a: u32) -> u32 { assert!(a != 0); (self.m - (a - 1)) % self.m + 1 }
    #[inline] fn neg(&self, a: u32) -> u32 { if a == 0 { 0 } else { (a - 1 + self.m / 2) % self.m + 1 } }
    #[inline] fn add(&self, a: u32, b: u32) -> u32 { if a == 0 { return b; } if b == 0 { return a; }
        let z = self.zech[((b + self.m - a) % self.m) as usize]; if z == 0 { 0 } else { (a - 1 + z - 1) % self.m + 1 } }
    #[inline] fn sub(&self, a: u32, b: u32) -> u32 { self.add(a, self.neg(b)) }
    #[inline] fn frob(&self, a: u32) -> u32 { if a == 0 { 0 } else { (((a - 1) as u64 * self.q as u64) % self.m as u64) as u32 + 1 } }
    #[inline] fn is_sq(&self, a: u32) -> bool { a == 0 || (a - 1) % 2 == 0 }
    #[inline] fn sqrt(&self, a: u32) -> u32 { if a == 0 { 0 } else { (a - 1) / 2 + 1 } }
    #[inline] fn in_fq(&self, a: u32) -> bool { a == 0 || (a - 1) % (self.m / (self.q - 1)) == 0 }
    fn int(&self, c: u32) -> u32 { self.el_of[(c % self.q) as usize] }
}

type Q2 = [u32; 3]; // quadratic c0 + c1 t + c2 t^2
fn ev(f: &Fld, u: &Q2, t: u32) -> u32 { f.add(f.add(u[0], f.mul(u[1], t)), f.mul(u[2], f.mul(t, t))) }
fn conj(f: &Fld, u: &Q2) -> Q2 { [f.frob(u[0]), f.frob(u[1]), f.frob(u[2])] }
fn roots(f: &Fld, u: &Q2) -> Vec<u32> { (0..f.size).filter(|&t| ev(f, u, t) == 0).collect() }
fn resultant_zero(f: &Fld, a: &Q2, b: &Q2) -> bool { // monic quadratics share a root iff resultant = 0
    let d0 = f.sub(a[0], b[0]); let d1 = f.sub(a[1], b[1]);
    f.sub(f.mul(d0, d0), f.mul(d1, f.sub(f.mul(a[1], b[0]), f.mul(b[1], a[0])))) == 0 }
type M3 = [[u32; 3]; 3];
fn det3(f: &Fld, a: &M3) -> u32 { let m = |x, y| f.mul(x, y);
    f.add(f.sub(m(a[0][0], f.sub(m(a[1][1], a[2][2]), m(a[1][2], a[2][1]))), m(a[0][1], f.sub(m(a[1][0], a[2][2]), m(a[1][2], a[2][0])))), m(a[0][2], f.sub(m(a[1][0], a[2][1]), m(a[1][1], a[2][0])))) }
fn inv3(f: &Fld, a: &M3) -> Option<M3> { let d = det3(f, a); if d == 0 { return None; } let di = f.inv(d); let mut r = [[0u32; 3]; 3];
    for i in 0..3 { for j in 0..3 { let (i1, i2, j1, j2) = ((i + 1) % 3, (i + 2) % 3, (j + 1) % 3, (j + 2) % 3);
        r[j][i] = f.mul(di, f.sub(f.mul(a[i1][j1], a[i2][j2]), f.mul(a[i1][j2], a[i2][j1]))); } } Some(r) }

// ternary polynomials: exponent (i, j, k) of X^i Y^j Z^k -> coefficient
type T = BTreeMap<(u8, u8, u8), u32>;
fn tadd(f: &Fld, a: &T, b: &T) -> T { let mut r = a.clone(); for (k, v) in b { let e = r.entry(*k).or_insert(0); *e = f.add(*e, *v); } r.retain(|_, v| *v != 0); r }
fn tmul(f: &Fld, a: &T, b: &T) -> T { let mut r = T::new(); for (k1, v1) in a { for (k2, v2) in b { let e = r.entry((k1.0 + k2.0, k1.1 + k2.1, k1.2 + k2.2)).or_insert(0); *e = f.add(*e, f.mul(*v1, *v2)); } } r.retain(|_, v| *v != 0); r }
fn tscale(f: &Fld, a: &T, c: u32) -> T { a.iter().map(|(k, v)| (*k, f.mul(*v, c))).filter(|x| x.1 != 0).collect() }
fn tlin(c: [u32; 3]) -> T { let mut r = T::new(); for (i, k) in [(1, 0, 0), (0, 1, 0), (0, 0, 1)].iter().enumerate() { if c[i] != 0 { r.insert(*k, c[i]); } } r }
fn teval(f: &Fld, a: &T, p: [u32; 3]) -> u32 { let pw = |x: u32, e: u8| { let mut r = 1; for _ in 0..e { r = f.mul(r, x); } r };
    let mut s = 0; for (k, v) in a { s = f.add(s, f.mul(*v, f.mul(pw(p[0], k.0), f.mul(pw(p[1], k.1), pw(p[2], k.2))))); } s }
fn tsubst(f: &Fld, a: &T, lin: &[T; 3]) -> T { let mut r = T::new(); let one: T = [((0, 0, 0), 1u32)].into_iter().collect();
    for (k, v) in a { let mut term = tscale(f, &one, *v); for (i, e) in [k.0, k.1, k.2].iter().enumerate() { for _ in 0..*e { term = tmul(f, &term, &lin[i]); } } r = tadd(f, &r, &term); } r }
fn count_proj(f: &Fld, a: &T) -> u64 { let mut c = 0u64; // points of the plane curve over the field f
    for y in 0..f.size { for z in 0..f.size { if teval(f, a, [1, y, z]) == 0 { c += 1; } } }
    for z in 0..f.size { if teval(f, a, [0, 1, z]) == 0 { c += 1; } } if teval(f, a, [0, 0, 1]) == 0 { c += 1; } c }
fn normalize(f: &Fld, p: [u32; 3]) -> [u32; 3] { let l = *p.iter().find(|&&x| x != 0).expect("zero point"); let li = f.inv(l); [f.mul(p[0], li), f.mul(p[1], li), f.mul(p[2], li)] }
fn is_prime(n: u64) -> bool { if n < 2 { return false; } let mut p = 2; while p * p <= n { if n % p == 0 { return false; } p += 1; } true }

fn build(f: &Fld, u: &Q2, kind: &str) -> Option<String> {
    let (u1, u2) = (conj(f, u), conj(f, &conj(f, u)));
    // E0 : y^2 = u u'' ; order from the quartic model (two points at infinity)
    let f0 = |t| f.mul(ev(f, u, t), ev(f, &u2, t)); let f1 = |t| f.mul(ev(f, u, t), ev(f, &u1, t));
    let mut ord = 2u64; for t in 0..f.size { let v = f0(t); ord += if v == 0 { 1 } else if f.is_sq(v) { 2 } else { 0 }; }
    if ord % 4 != 0 || !is_prime(ord / 4) { return None; }
    let rinv = inv3(f, &[*u, u1, u2])?; // rows R_0, R_1, R_2 with (1, t, t^2) = R (u, u', u'')
    let sq = |k: (u8, u8, u8)| -> T { [(k, 1u32)].into_iter().collect() };
    let w = [sq((0, 0, 2)), sq((2, 0, 0)), sq((0, 2, 0))]; // (u : u' : u'') = (Z^2 : X^2 : Y^2)
    let row = |r: &[u32; 3]| -> T { let mut s = T::new(); for i in 0..3 { s = tadd(f, &s, &tscale(f, &w[i], r[i])); } s };
    let (r0, r1, r2) = (row(&rinv[0]), row(&rinv[1]), row(&rinv[2]));
    let g = tadd(f, &tmul(f, &r1, &r1), &tscale(f, &tmul(f, &r0, &r2), f.neg(1)));
    // Ciani smoothness on the K model: a X^4 + b Y^4 + c Z^4 + d X^2Y^2 + e Y^2Z^2 + h Z^2X^2
    let co = |k| *g.get(&k).unwrap_or(&0); let (a, b, c, d, e, h) = (co((4, 0, 0)), co((0, 4, 0)), co((0, 0, 4)), co((2, 2, 0)), co((0, 2, 2)), co((2, 0, 2)));
    let four = f.int(4); let disc = |x, y, z| f.sub(f.mul(z, z), f.mul(four, f.mul(x, y)));
    let conic_det = { let t2 = f.inv(f.int(2)); det3(f, &[[a, f.mul(d, t2), f.mul(h, t2)], [f.mul(d, t2), b, f.mul(e, t2)], [f.mul(h, t2), f.mul(e, t2), c]]) };
    let smooth = a != 0 && b != 0 && c != 0 && disc(a, b, d) != 0 && disc(b, c, e) != 0 && disc(a, c, h) != 0 && conic_det != 0;
    // descent: W_j = w_j X + w_j^s Y + w_j^(s^2) Z with w_j = 1, g, g^2 ; substitute (X, Y, Z) = Omega^-1 W
    let om: M3 = { let mut m = [[0u32; 3]; 3]; for j in 0..3 { let wj = if j == 0 { 1 } else { j as u32 + 1 }; m[j] = [wj, f.frob(wj), f.frob(f.frob(wj))]; } m };
    let omi = inv3(f, &om).expect("basis"); let lin = [tlin(omi[0]), tlin(omi[1]), tlin(omi[2])];
    let gw = tsubst(f, &g, &lin); let lead = *gw.values().next().unwrap(); let gwn = tscale(f, &gw, f.inv(lead));
    let descends = gwn.values().all(|&v| f.in_fq(v));
    let ints: Vec<String> = gwn.iter().map(|(k, v)| format!("[{},{},{},{}]", k.0, k.1, k.2, if f.in_fq(*v) { f.idx_of[*v as usize] as i64 } else { -1 })).collect();
    // point counts of the F_q quartic over F_q, F_q^2, F_q^3 (and of the K model over K as a cross-check)
    let mut n = [0u64; 3];
    if descends { for k in 1..=3u32 { let fk = if k == 3 { None } else { Some(Fld::build(f.q, k)) }; let fr = fk.as_ref().unwrap_or(f);
        let gk: T = gwn.iter().map(|(key, v)| (*key, fr.int(f.idx_of[*v as usize]))).collect(); n[k as usize - 1] = count_proj(fr, &gk); } }
    let n3_kmodel = count_proj(f, &g);
    let q = f.q as i64; let p: Vec<i64> = (0..3).map(|k| q.pow(k as u32 + 1) + 1 - n[k] as i64).collect();
    let e1 = p[0]; let e2 = (e1 * p[0] - p[1]) / 2; let e3 = (e2 * p[0] - e1 * p[1] + p[2]) / 3;
    let jac = 1 - e1 + e2 - e3 + q * e2 - q * q * e1 + q * q * q; let te = q * q * q + 1 - ord as i64;
    // transfer divisors
    let xyz = |y: [u32; 3]| [f.mul(y[1], y[2]), f.mul(y[0], y[2]), f.mul(y[0], y[1])];
    let to_w = |p: [u32; 3]| { let mut r = [0u32; 3]; for j in 0..3 { for i in 0..3 { r[j] = f.add(r[j], f.mul(om[j][i], p[i])); } } r };
    let (mut tested, mut on_curve, mut frob_ok, mut on_e0, mut skipped) = (0u64, 0u64, 0u64, 0u64, 0u64);
    for t in 0..f.size { let (v0, v1) = (f0(t), f1(t)); if v0 == 0 || !f.is_sq(v0) { continue; }
        if v1 == 0 || !f.is_sq(v1) || ev(f, u, t) == 0 { skipped += 2; continue; }
        for s0 in 0..2 { let y0 = if s0 == 0 { f.sqrt(v0) } else { f.neg(f.sqrt(v0)) }; // P = (t, y0) on E0
            for s1 in 0..2 { let y1 = if s1 == 0 { f.sqrt(v1) } else { f.neg(f.sqrt(v1)) }; let y2 = f.mul(f.mul(y0, y1), f.inv(ev(f, u, t)));
                let mut pt = (t, [y0, y1, y2]);
                for _ in 0..3 { tested += 1; let wp = to_w(xyz(pt.1));
                    if teval(f, &gwn, wp) == 0 && teval(f, &g, xyz(pt.1)) == 0 { on_curve += 1; }
                    if f.mul(pt.1[0], pt.1[0]) == f0(pt.0) { on_e0 += 1; }
                    let nx = (f.frob(pt.0), [f.frob(pt.1[2]), f.frob(pt.1[0]), f.frob(pt.1[1])]);
                    let wq = [f.frob(wp[0]), f.frob(wp[1]), f.frob(wp[2])];
                    if normalize(f, wq) == normalize(f, to_w(xyz(nx.1))) { frob_ok += 1; }
                    pt = nx; } } } }
    Some(format!("{{\"rec\":\"case\",\"q\":{},\"kind\":\"{}\",\"u_idx\":[{},{}],\"roots_in_K\":{},\"order_E0\":{},\"ell\":{},\"trace_E0\":{},\"descends\":{},\"smooth\":{},\"quartic_Fq\":[{}],\"N\":[{},{},{}],\"N3_K_model\":{},\"N_predicted\":[{},{},{}],\"e\":[{},{},{}],\"jac_order\":{},\"jac_equals_E0\":{},\"divisor_points_tested\":{},\"on_quartic\":{},\"frobenius_matches\":{},\"pushforward_on_E0\":{},\"points_skipped\":{}}}",
        f.q, kind, f.idx_of[u[1] as usize], f.idx_of[u[0] as usize], roots(f, u).len(), ord, ord / 4, te, descends, smooth, ints.join(","), n[0], n[1], n[2], n3_kmodel,
        q + 1, q * q + 1, q * q * q + 1 - 3 * te, e1, e2, e3, jac, jac == ord as i64, tested, on_curve, frob_ok, on_e0, skipped))
}

fn main() { let a: Vec<String> = std::env::args().collect(); let q: u32 = a[1].parse().unwrap(); let want = [a[2].parse::<usize>().unwrap(), a[3].parse::<usize>().unwrap()];
    let f = Fld::build(q, 3); let mut got = [0usize; 2]; let (mut tried, mut dependent) = (0u64, 0u64);
    'o: for bi in 0..f.size { for ci in 0..f.size { let u: Q2 = [f.el_of[ci as usize], f.el_of[bi as usize], 1]; let (u1, u2) = (conj(&f, &u), conj(&f, &conj(&f, &u)));
        let disc = f.sub(f.mul(u[1], u[1]), f.mul(f.int(4), u[0])); if disc == 0 { continue; }
        if resultant_zero(&f, &u, &u1) || resultant_zero(&f, &u, &u2) || resultant_zero(&f, &u1, &u2) { continue; }
        let k = if f.is_sq(disc) { 0 } else { 1 }; if got[k] >= want[k] { continue; }
        tried += 1; if det3(&f, &[u, u1, u2]) == 0 { dependent += 1; continue; }
        if let Some(s) = build(&f, &u, if k == 0 { "split" } else { "irreducible" }) { println!("{}", s); got[k] += 1; }
        if got[0] >= want[0] && got[1] >= want[1] { break 'o; } } }
    println!("{{\"rec\":\"search\",\"q\":{},\"candidates_examined\":{},\"linearly_dependent_skipped\":{},\"split_found\":{},\"irreducible_found\":{}}}", q, tried, dependent, got[0], got[1]);
}
