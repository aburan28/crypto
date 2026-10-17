// 12-model-small-x: small-coordinate factor bases in Montgomery / twisted Edwards / Hessian
// models vs Weierstrass small-x and random controls, m = 2 coverage instrument.
// Single-threaded, std only. See PROTOCOL.md (frozen before data).
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::time::Instant;

// ---------------------------------------------------------------- arithmetic
#[derive(Clone, Copy)]
struct F {
    p: u64,
}
impl F {
    #[inline]
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p { s - self.p } else { s }
    }
    #[inline]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b { a - b } else { a + self.p - b }
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    #[inline]
    fn neg(&self, a: u64) -> u64 {
        if a == 0 { 0 } else { self.p - a }
    }
    fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "inverse of zero");
        let (mut r0, mut r1) = (self.p as i128, a as i128);
        let (mut t0, mut t1) = (0i128, 1i128);
        while r1 != 0 {
            let q = r0 / r1;
            let r2 = r0 - q * r1;
            r0 = r1;
            r1 = r2;
            let t2 = t0 - q * t1;
            t0 = t1;
            t1 = t2;
        }
        let p = self.p as i128;
        (((t0 % p) + p) % p) as u64
    }
    fn pow(&self, mut a: u64, mut e: u128) -> u64 {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
    fn from_i(&self, v: i64) -> u64 {
        let p = self.p as i64;
        (((v % p) + p) % p) as u64
    }
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}
fn seed_of(parts: &[u64]) -> u64 {
    let mut r = Rng(46012);
    let mut s = r.next();
    for &x in parts {
        let mut t = Rng(s ^ x.wrapping_mul(0x2545F4914F6CDD1D));
        s = t.next();
    }
    s
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            return false;
        }
        d += 1;
    }
    true
}

// ---------------------------------------------------------------- curves
const INF: u64 = u64::MAX;
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Pt {
    x: u64,
    y: u64,
}
const O: Pt = Pt { x: INF, y: 0 };

struct W {
    f: F,
    a: u64,
    b: u64,
}
impl W {
    fn rhs(&self, x: u64) -> u64 {
        let f = &self.f;
        f.add(f.mul(f.add(f.mul(x, x), self.a), x), self.b)
    }
    fn on(&self, p: Pt) -> bool {
        p.x == INF || self.f.mul(p.y, p.y) == self.rhs(p.x)
    }
    fn neg(&self, p: Pt) -> Pt {
        if p.x == INF { p } else { Pt { x: p.x, y: self.f.neg(p.y) } }
    }
    #[inline]
    fn add(&self, p: Pt, q: Pt) -> Pt {
        if p.x == INF {
            return q;
        }
        if q.x == INF {
            return p;
        }
        let f = &self.f;
        let l;
        if p.x == q.x {
            if p.y != q.y || p.y == 0 {
                return O;
            }
            let num = f.add(f.mul(3, f.mul(p.x, p.x)), self.a);
            l = f.mul(num, f.inv(f.add(p.y, p.y)));
        } else {
            l = f.mul(f.sub(q.y, p.y), f.inv(f.sub(q.x, p.x)));
        }
        let x3 = f.sub(f.sub(f.mul(l, l), p.x), q.x);
        let y3 = f.sub(f.mul(l, f.sub(p.x, x3)), p.y);
        Pt { x: x3, y: y3 }
    }
    fn mul(&self, p: Pt, k: u128) -> Pt {
        let mut r = O;
        let mut q = p;
        let mut k = k;
        while k > 0 {
            if k & 1 == 1 {
                r = self.add(r, q);
            }
            q = self.add(q, q);
            k >>= 1;
        }
        r
    }
}

// ---------------------------------------------------------------- Part A: pair 46
fn polmulmod(f: &F, a: &[u64; 3], b: &[u64; 3], m: &[u64; 4]) -> [u64; 3] {
    // m monic cubic: x^3 + m2 x^2 + m1 x + m0 ; m[3]=1
    let mut c = [0u64; 5];
    for i in 0..3 {
        for j in 0..3 {
            c[i + j] = f.add(c[i + j], f.mul(a[i], b[j]));
        }
    }
    for d in (3..5).rev() {
        let t = c[d];
        if t != 0 {
            c[d] = 0;
            for k in 0..3 {
                c[d - 3 + k] = f.sub(c[d - 3 + k], f.mul(t, m[k]));
            }
        }
    }
    [c[0], c[1], c[2]]
}
fn poly_trim(v: &mut Vec<u64>) {
    while v.len() > 0 && *v.last().unwrap() == 0 {
        v.pop();
    }
}
fn poly_gcd(f: &F, a: Vec<u64>, b: Vec<u64>) -> Vec<u64> {
    let (mut a, mut b) = (a, b);
    poly_trim(&mut a);
    poly_trim(&mut b);
    while !b.is_empty() {
        // a mod b
        let ib = f.inv(*b.last().unwrap());
        while a.len() >= b.len() && !a.is_empty() {
            let t = f.mul(*a.last().unwrap(), ib);
            let sh = a.len() - b.len();
            for i in 0..b.len() {
                a[sh + i] = f.sub(a[sh + i], f.mul(t, b[i]));
            }
            poly_trim(&mut a);
        }
        std::mem::swap(&mut a, &mut b);
    }
    a
}

fn part_a(out: &mut String) {
    let p: u64 = 5947698030651343;
    let n: u128 = 5947697996161643;
    let t: u128 = 34489701;
    let v: u128 = t + 2;
    let f = F { p };
    let curves = [
        ("crater", 5757116383569871u64, 3560173328417442u64),
        ("floor", 2087370704375809u64, 1402887937802154u64),
    ];
    let tw = 5 * v * v;
    let mut s = String::new();
    write!(
        s,
        "{{\"p\":{},\"N\":{},\"twist_order\":{},\"twist_order_check_p_plus_1_plus_t\":{},\"N_mod2\":{},\"N_mod3\":{},\"N_mod4\":{},\"twist_mod2\":{},\"twist_mod3\":{},\"v_mod3\":{},\"curves\":[",
        p, n, tw, (p as u128 + 1 + t) == tw, n % 2, n % 3, n % 4, tw % 2, tw % 3, v % 3
    )
    .unwrap();
    let mut rng = Rng(seed_of(&[99]));
    for (ci, (name, a, b)) in curves.iter().enumerate() {
        let w = W { f, a: *a, b: *b };
        // x^p mod (x^3 + a x + b)
        let m = [*b, *a, 0, 1];
        let mut r = [1u64, 0, 0];
        let mut base = [0u64, 1, 0];
        let mut e = p;
        while e > 0 {
            if e & 1 == 1 {
                r = polmulmod(&f, &r, &base, &m);
            }
            base = polmulmod(&f, &base, &base, &m);
            e >>= 1;
        }
        let g = vec![r[0], f.sub(r[1], 1), r[2]];
        let gg = poly_gcd(&f, vec![*b, *a, 0, 1], g);
        let gcd_deg = gg.len() as i64 - 1;
        // N * P = O on random points
        let mut np_ok = 0;
        let mut tries = 0;
        while np_ok < 3 && tries < 100 {
            tries += 1;
            let x = rng.below(p);
            let rh = w.rhs(x);
            if f.pow(rh, ((p - 1) / 2) as u128) != 1 {
                continue;
            }
            // p mod 4 ?
            let y = if p % 4 == 3 { f.pow(rh, ((p + 1) / 4) as u128) } else { 0 };
            if y == 0 || f.mul(y, y) != rh {
                continue;
            }
            let pt = Pt { x, y };
            if w.mul(pt, n) == O {
                np_ok += 1;
            } else {
                np_ok = -100;
                break;
            }
        }
        if ci > 0 {
            s.push(',');
        }
        write!(
            s,
            "{{\"name\":\"{}\",\"gcd_deg_xp_minus_x_cubic\":{},\"has_2torsion_over_Fp\":{},\"N_times_P_is_O_on_3_random_points\":{}}}",
            name,
            gcd_deg,
            gcd_deg > 0,
            np_ok == 3
        )
        .unwrap();
    }
    write!(s, "],\"p_mod4\":{},\"p_mod3\":{}}}", p % 4, p % 3).unwrap();
    out.push_str(&s);
}

// ---------------------------------------------------------------- toy curve data
struct Toy {
    w: W,
    pts: Vec<Pt>,          // all affine W points
    n: u64,                // #E incl. identity
    tags: Vec<(String, Vec<u64>)>, // per model coordinate: tag per point (INF = undefined)
    sym_tag: Vec<bool>,
}

fn sqrt_table(p: u64) -> Vec<u32> {
    let mut t = vec![u32::MAX; p as usize];
    for r in 0..=(p / 2) {
        let v = ((r as u128 * r as u128) % p as u128) as usize;
        t[v] = r as u32;
    }
    t
}

fn enum_w(w: &W, sq: &[u32]) -> Vec<Pt> {
    let p = w.f.p;
    let mut v = Vec::with_capacity(p as usize + 2);
    for x in 0..p {
        let r = w.rhs(x);
        let s = sq[r as usize];
        if s == u32::MAX {
            continue;
        }
        let s = s as u64;
        if s == 0 {
            v.push(Pt { x, y: 0 });
        } else {
            v.push(Pt { x, y: s });
            v.push(Pt { x, y: p - s });
        }
    }
    v
}

fn make_montgomery(f: F, sq: &[u32], rng: &mut Rng) -> Result<Toy, String> {
    let p = f.p;
    let (aa, bb) = loop {
        let aa = rng.below(p);
        let bb = rng.below(p);
        if bb != 0 && f.mul(aa, aa) != 4 {
            break (aa, bb);
        }
    };
    let i3 = f.inv(3);
    let ib = f.inv(bb);
    let a = f.mul(f.sub(3, f.mul(aa, aa)), f.mul(i3, f.mul(ib, ib)));
    let b = f.mul(
        f.sub(f.mul(2, f.mul(aa, f.mul(aa, aa))), f.mul(9, aa)),
        f.mul(f.inv(27), f.mul(ib, f.mul(ib, ib))),
    );
    let w = W { f, a, b };
    let pts = enum_w(&w, sq);
    let n = pts.len() as u64 + 1;
    let a3 = f.mul(aa, i3);
    let ae = f.mul(f.add(aa, 2), ib);
    let de = f.mul(f.sub(aa, 2), ib);
    let mut tu = Vec::with_capacity(pts.len());
    let mut tey = Vec::with_capacity(pts.len());
    let mut tex = Vec::with_capacity(pts.len());
    for pt in &pts {
        let u = f.sub(f.mul(bb, pt.x), a3);
        let v = f.mul(bb, pt.y);
        // verify Montgomery equation B v^2 = u^3 + A u^2 + u
        let lhs = f.mul(bb, f.mul(v, v));
        let rhs = f.add(f.mul(f.add(f.mul(u, u), f.mul(aa, u)), u), u);
        if lhs != rhs {
            return Err("montgomery map verification failed".into());
        }
        tu.push(u);
        let up1 = f.add(u, 1);
        if up1 == 0 || v == 0 {
            // Edwards exceptional points (u = -1 or v = 0): not in Edwards affine chart
            let ye = if up1 == 0 { INF } else { f.mul(f.sub(u, 1), f.inv(up1)) };
            tey.push(ye);
            tex.push(INF);
            continue;
        }
        let ye = f.mul(f.sub(u, 1), f.inv(up1));
        let xe = f.mul(u, f.inv(v));
        // verify twisted Edwards a x^2 + y^2 = 1 + d x^2 y^2
        let x2 = f.mul(xe, xe);
        let y2 = f.mul(ye, ye);
        if f.add(f.mul(ae, x2), y2) != f.add(1, f.mul(de, f.mul(x2, y2))) {
            return Err("edwards map verification failed".into());
        }
        tey.push(ye);
        tex.push(xe);
    }
    let tw: Vec<u64> = pts.iter().map(|q| q.x).collect();
    Ok(Toy {
        w,
        pts,
        n,
        tags: vec![
            ("W-x".into(), tw),
            ("M-u".into(), tu),
            ("E-y".into(), tey),
            ("E-x".into(), tex),
        ],
        sym_tag: vec![true, true, true, false],
    })
}

fn solve4(f: &F, mut m: Vec<[u64; 5]>) -> Option<[u64; 4]> {
    for c in 0..4 {
        let piv = (c..4).find(|&r| m[r][c] != 0)?;
        m.swap(c, piv);
        let iv = f.inv(m[c][c]);
        for k in 0..5 {
            m[c][k] = f.mul(m[c][k], iv);
        }
        for r in 0..4 {
            if r != c && m[r][c] != 0 {
                let t = m[r][c];
                for k in 0..5 {
                    m[r][k] = f.sub(m[r][k], f.mul(t, m[c][k]));
                }
            }
        }
    }
    Some([m[0][4], m[1][4], m[2][4], m[3][4]])
}

fn make_hessian(f: F, sq: &[u32], rng: &mut Rng) -> Result<Toy, String> {
    let p = f.p;
    let d = loop {
        let d = rng.below(p);
        if d != 0 && f.mul(d, f.mul(d, d)) != 1 {
            break d;
        }
    };
    let on_h = |x: u64, y: u64| -> bool {
        f.add(f.add(f.mul(x, f.mul(x, x)), f.mul(y, f.mul(y, y))), 1) == f.mul(3, f.mul(d, f.mul(x, y)))
    };
    let msqrt = |v: u64| -> Option<u64> {
        let s = sq[v as usize];
        if s == u32::MAX { None } else { Some(s as u64) }
    };
    let mut hs: HashSet<(u64, u64)> = HashSet::new();
    // x = 0: y^3 = -1
    for y in [f.neg(1)] {
        hs.insert((0, y));
    }
    // y^2 - y + 1 = 0 -> y = (1 +- sqrt(-3))/2
    let i2 = f.inv(2);
    if let Some(s) = msqrt(f.neg(3)) {
        for ss in [s, f.neg(s)] {
            let y = f.mul(f.add(1, ss), i2);
            if on_h(0, y) {
                hs.insert((0, y));
            }
        }
    }
    // lines y = m x - 1 through (0,-1)
    for m in 0..p {
        let a2 = f.add(1, f.mul(m, f.mul(m, m)));
        let a1 = f.neg(f.mul(3, f.add(f.mul(m, m), f.mul(d, m))));
        let a0 = f.mul(3, f.add(m, d));
        let mut xs: Vec<u64> = Vec::new();
        if a2 == 0 {
            if a1 != 0 {
                xs.push(f.mul(f.neg(a0), f.inv(a1)));
            }
        } else {
            let disc = f.sub(f.mul(a1, a1), f.mul(4, f.mul(a2, a0)));
            if let Some(s) = msqrt(disc) {
                let i = f.inv(f.add(a2, a2));
                xs.push(f.mul(f.add(f.neg(a1), s), i));
                xs.push(f.mul(f.sub(f.neg(a1), s), i));
            }
        }
        for x in xs {
            if x == 0 {
                continue;
            }
            let y = f.sub(f.mul(m, x), 1);
            if !on_h(x, y) {
                return Err("hessian enumeration produced off-curve point".into());
            }
            hs.insert((x, y));
        }
    }
    let mut haff: Vec<(u64, u64)> = hs.into_iter().collect();
    haff.sort();
    let extra_inf = if msqrt(f.neg(3)).is_some() { 2 } else { 0 };
    let n_h = haff.len() as u64 + 1 + extra_inf;
    // u0 = 1/(x+y+d), v0 = (y-x) u0
    let mut uv: Vec<(u64, u64)> = Vec::with_capacity(haff.len());
    for &(x, y) in &haff {
        let den = f.add(f.add(x, y), d);
        if den == 0 {
            return Err("affine point on tangent line at O (unexpected)".into());
        }
        let u0 = f.inv(den);
        uv.push((u0, f.mul(f.sub(y, x), u0)));
    }
    // fit v0^2 = c3 u0^3 + c2 u0^2 + c1 u0 + c0 from 4 points with distinct u0
    let mut rows: Vec<[u64; 5]> = Vec::new();
    let mut used: HashSet<u64> = HashSet::new();
    for &(u0, v0) in &uv {
        if used.insert(u0) {
            rows.push([f.mul(u0, f.mul(u0, u0)), f.mul(u0, u0), u0, 1, f.mul(v0, v0)]);
            if rows.len() == 4 {
                break;
            }
        }
    }
    let c = solve4(&f, rows).ok_or("singular fit")?;
    let (c3, c2, c1, c0) = (c[0], c[1], c[2], c[3]);
    if c3 == 0 {
        return Err("c3 = 0".into());
    }
    for &(u0, v0) in &uv {
        let r = f.add(f.add(f.mul(c3, f.mul(u0, f.mul(u0, u0))), f.mul(c2, f.mul(u0, u0))), f.add(f.mul(c1, u0), c0));
        if f.mul(v0, v0) != r {
            return Err("cubic relation verification failed".into());
        }
    }
    let i3 = f.inv(3);
    let c1c3 = f.mul(c1, c3);
    let a = f.sub(c1c3, f.mul(f.mul(c2, c2), i3));
    let b = f.add(
        f.sub(f.mul(f.mul(2, f.mul(c2, f.mul(c2, c2))), f.inv(27)), f.mul(f.mul(c2, c1c3), i3)),
        f.mul(c0, f.mul(c3, c3)),
    );
    let w = W { f, a, b };
    let c2_3 = f.mul(c2, i3);
    let mut pts = Vec::with_capacity(haff.len());
    let mut tx = Vec::with_capacity(haff.len());
    let mut ts = Vec::with_capacity(haff.len());
    let mut keys: HashSet<(u64, u64)> = HashSet::new();
    for (i, &(u0, v0)) in uv.iter().enumerate() {
        let pt = Pt { x: f.add(f.mul(c3, u0), c2_3), y: f.mul(c3, v0) };
        if !w.on(pt) {
            return Err("weierstrass image off curve".into());
        }
        if !keys.insert((pt.x, pt.y)) {
            return Err("map not injective".into());
        }
        pts.push(pt);
        tx.push(haff[i].0);
        ts.push(f.add(haff[i].0, haff[i].1));
    }
    // count W points independently
    let wall = enum_w(&w, sq);
    if wall.len() as u64 + 1 != n_h {
        return Err(format!("point count mismatch W {} vs H {}", wall.len() + 1, n_h));
    }
    // homomorphism spot check: collinear triples on x_H = c sum to O
    let mut byx: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, &(x, _)) in haff.iter().enumerate() {
        if x < 2000 {
            byx.entry(x).or_default().push(i);
        }
    }
    for (_, v) in byx.iter() {
        if v.len() == 3 {
            let s = w.add(w.add(pts[v[0]], pts[v[1]]), pts[v[2]]);
            if s != O {
                return Err("collinear triple does not sum to O".into());
            }
        }
    }
    // negation check: (x,y)->(y,x) maps to -P
    for (i, &(x, y)) in haff.iter().enumerate().take(50) {
        let j = haff.binary_search(&(y, x)).map_err(|_| "neg partner missing")?;
        if w.neg(pts[i]) != pts[j] {
            return Err("negation check failed".into());
        }
    }
    // add the two non-O points at infinity (u0 = 0) to the W point list (no Hessian tag)
    let mut pts = pts;
    let mut tx = tx;
    let mut ts = ts;
    for pt in wall.iter() {
        if !keys.contains(&(pt.x, pt.y)) {
            pts.push(*pt);
            tx.push(INF);
            ts.push(INF);
        }
    }
    if pts.len() as u64 + 1 != n_h {
        return Err("bijection completion failed".into());
    }
    let tw: Vec<u64> = pts.iter().map(|q| q.x).collect();
    Ok(Toy {
        w,
        pts,
        n: n_h,
        tags: vec![("W-x".into(), tw), ("H-x".into(), tx), ("H-s".into(), ts)],
        sym_tag: vec![true, false, true],
    })
}

// ---------------------------------------------------------------- factor bases
fn small_tag_set(tags: &[u64], k: usize) -> Vec<usize> {
    let mut bound = (8 * k) as u64;
    loop {
        let mut v: Vec<(u64, usize)> = tags
            .iter()
            .enumerate()
            .filter(|(_, &t)| t != INF && t < bound)
            .map(|(i, &t)| (t, i))
            .collect();
        if v.len() >= k || bound > tags.len() as u64 * 4 {
            v.sort();
            v.truncate(k);
            return v.into_iter().map(|(_, i)| i).collect();
        }
        bound *= 2;
    }
}

fn key(p: Pt, pp: u64) -> usize {
    if p.x == INF {
        (2 * pp) as usize
    } else {
        (2 * p.x + if p.y > pp / 2 { 1 } else { 0 }) as usize
    }
}

fn rand_sym(t: &Toy, k: usize, rng: &mut Rng, sq: &[u32], exclude: &HashSet<usize>, out: &mut Vec<Pt>) {
    let p = t.w.f.p;
    let mut seen: HashSet<u64> = HashSet::new();
    while out.len() + 1 < k + 1 && out.len() < k {
        let x = rng.below(p);
        if seen.contains(&x) {
            continue;
        }
        let s = sq[t.w.rhs(x) as usize];
        if s == u32::MAX || s == 0 {
            continue;
        }
        let a = Pt { x, y: s as u64 };
        let b = Pt { x, y: p - s as u64 };
        if exclude.contains(&key(a, p)) || exclude.contains(&key(b, p)) {
            continue;
        }
        seen.insert(x);
        out.push(a);
        if out.len() < k {
            out.push(b);
        }
    }
}

fn random_point(t: &Toy, rng: &mut Rng) -> Pt {
    t.pts[rng.below(t.pts.len() as u64) as usize]
}

fn dl_interval(t: &Toy, half: usize, rng: &mut Rng) -> Vec<Pt> {
    let p = t.w.f.p;
    'outer: loop {
        let g = random_point(t, rng);
        let mut v = Vec::with_capacity(2 * half);
        let mut ks: HashSet<usize> = HashSet::new();
        let mut q = O;
        for _ in 0..half {
            q = t.w.add(q, g);
            if q == O || !ks.insert(key(q, p)) || !ks.insert(key(t.w.neg(q), p)) {
                continue 'outer;
            }
            v.push(q);
            v.push(t.w.neg(q));
        }
        return v;
    }
}

struct Cov {
    cover: u64,
    pairs: u64,
    semaev_checked: u64,
    semaev_disagree: u64,
    semaev_degenerate: u64,
}

fn coverage(t: &Toy, fb: &[Pt], semaev: bool, sq: &[u32], rng: &mut Rng) -> Cov {
    let p = t.w.f.p;
    let nbits = (2 * p + 1) as usize;
    let mut bits = vec![0u64; nbits / 64 + 1];
    let mut cover = 0u64;
    let mut pairs = 0u64;
    for i in 0..fb.len() {
        for j in i..fb.len() {
            let s = t.w.add(fb[i], fb[j]);
            let kk = key(s, p);
            let (wd, b) = (kk / 64, kk % 64);
            if bits[wd] >> b & 1 == 0 {
                bits[wd] |= 1 << b;
                cover += 1;
            }
            pairs += 1;
        }
    }
    let mut c = Cov { cover, pairs, semaev_checked: 0, semaev_disagree: 0, semaev_degenerate: 0 };
    if semaev {
        let f = &t.w.f;
        let (a, b) = (t.w.a, t.w.b);
        let xs: HashSet<u64> = fb.iter().map(|q| q.x).collect();
        let xv: Vec<u64> = xs.iter().cloned().collect();
        for _ in 0..200 {
            let r = random_point(t, rng);
            let kk = key(r, p);
            let kn = key(t.w.neg(r), p);
            let enum_dec = (bits[kk / 64] >> (kk % 64) & 1 == 1) || (bits[kn / 64] >> (kn % 64) & 1 == 1);
            let x3 = r.x;
            let mut sem = false;
            for &x1 in &xv {
                let a2 = f.mul(f.sub(x1, x3), f.sub(x1, x3));
                let a1 = f.neg(f.add(
                    f.add(f.mul(2, f.mul(x1, f.mul(x3, x3))), f.mul(2, f.mul(x3, f.add(f.mul(x1, x1), a)))),
                    f.add(f.mul(2, f.mul(a, x1)), f.mul(4, b)),
                ));
                let a0 = f.sub(
                    f.add(f.mul(f.mul(x1, x1), f.mul(x3, x3)), f.mul(a, a)),
                    f.add(f.mul(2, f.mul(x3, f.add(f.mul(a, x1), f.mul(2, b)))), f.mul(4, f.mul(b, x1))),
                );
                if a2 == 0 {
                    if a1 != 0 {
                        let x2 = f.mul(f.neg(a0), f.inv(a1));
                        if xs.contains(&x2) {
                            sem = true;
                            break;
                        }
                    } else if a0 == 0 {
                        c.semaev_degenerate += 1;
                    }
                    continue;
                }
                let disc = f.sub(f.mul(a1, a1), f.mul(4, f.mul(a2, a0)));
                let s = sq[disc as usize];
                if s == u32::MAX {
                    continue;
                }
                let i = f.inv(f.add(a2, a2));
                for ss in [s as u64, f.neg(s as u64)] {
                    let x2 = f.mul(f.add(f.neg(a1), ss), i);
                    if xs.contains(&x2) {
                        sem = true;
                    }
                }
                if sem {
                    break;
                }
            }
            c.semaev_checked += 1;
            if sem != enum_dec {
                c.semaev_disagree += 1;
            }
        }
    }
    c
}

// ---------------------------------------------------------------- stats
fn bootstrap(v: &[f64], seed: u64, qs: &[f64]) -> (f64, Vec<f64>) {
    let n = v.len();
    let mean = v.iter().sum::<f64>() / n as f64;
    let mut rng = Rng(seed);
    let reps = 20000;
    let mut ms = Vec::with_capacity(reps);
    for _ in 0..reps {
        let mut s = 0.0;
        for _ in 0..n {
            s += v[rng.below(n as u64) as usize];
        }
        ms.push(s / n as f64);
    }
    ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let out = qs.iter().map(|&q| ms[((q * reps as f64) as usize).min(reps - 1)]).collect();
    (mean, out)
}

// ---------------------------------------------------------------- main
fn main() {
    let t0 = Instant::now();
    let args: Vec<String> = std::env::args().collect();
    let sizes: Vec<u32> = if args.len() > 1 {
        args[1].split(',').map(|s| s.parse().unwrap()).collect()
    } else {
        vec![16, 20, 22]
    };
    let budget_s: f64 = if args.len() > 2 { args[2].parse().unwrap() } else { 2100.0 };
    let mut json = String::from("{\"protocol\":\"PROTOCOL.md\",\"part_a\":");
    part_a(&mut json);
    eprintln!("part A done");
    json.push_str(",\"records\":[");
    let mut first = true;
    // rel values per (pop, b, arm)
    let mut rels: Vec<(String, u32, String, String, Vec<f64>)> = Vec::new();
    let mut extra: Vec<String> = Vec::new();
    let mut completed: Vec<String> = Vec::new();
    let mut failures: Vec<String> = Vec::new();
    'sizes: for &bits in &sizes {
        let mut p = (1u64 << bits) - 1;
        while !is_prime(p) {
            p -= 1;
        }
        let f = F { p };
        let sq = sqrt_table(p);
        let ncurves = match bits {
            16 => 60,
            20 => 40,
            _ => 20,
        };
        let k = 2 * ((p as f64).sqrt() / 2.0).round() as usize;
        for pop in ["M", "H"] {
            let mut done = 0;
            for ci in 0..ncurves {
                if t0.elapsed().as_secs_f64() > budget_s {
                    failures.push(format!("budget stop at pop {} b {} curve {}", pop, bits, ci));
                    break 'sizes;
                }
                let popid = if pop == "M" { 1 } else { 2 };
                let mut rng = Rng(seed_of(&[popid, bits as u64, ci as u64, 0]));
                let toy = if pop == "M" { make_montgomery(f, &sq, &mut rng) } else { make_hessian(f, &sq, &mut rng) };
                let toy = match toy {
                    Ok(t) => t,
                    Err(e) => {
                        failures.push(format!("pop {} b {} curve {}: {}", pop, bits, ci, e));
                        continue;
                    }
                };
                // arms
                let mut arms: Vec<(String, Vec<Pt>, bool)> = Vec::new();
                for (ti, (name, tags)) in toy.tags.iter().enumerate() {
                    let idx = small_tag_set(tags, k);
                    arms.push((name.clone(), idx.into_iter().map(|i| toy.pts[i]).collect(), toy.sym_tag[ti]));
                }
                let empty = HashSet::new();
                let mut r1 = Rng(seed_of(&[popid, bits as u64, ci as u64, 1]));
                let mut v = Vec::new();
                rand_sym(&toy, k, &mut r1, &sq, &empty, &mut v);
                arms.push(("R-sym".into(), v, true));
                let mut r2 = Rng(seed_of(&[popid, bits as u64, ci as u64, 2]));
                let mut v = Vec::new();
                rand_sym(&toy, k, &mut r2, &sq, &empty, &mut v);
                arms.push(("R-sym2".into(), v, true));
                let mut r3 = Rng(seed_of(&[popid, bits as u64, ci as u64, 3]));
                let mut chosen: HashSet<usize> = HashSet::new();
                let mut v = Vec::new();
                while v.len() < k {
                    let i = r3.below(toy.pts.len() as u64) as usize;
                    if chosen.insert(i) {
                        v.push(toy.pts[i]);
                    }
                }
                arms.push(("R-asym".into(), v, false));
                let mut r4 = Rng(seed_of(&[popid, bits as u64, ci as u64, 4]));
                arms.push(("PC1".into(), dl_interval(&toy, k / 2, &mut r4), true));
                let mut r5 = Rng(seed_of(&[popid, bits as u64, ci as u64, 5]));
                let mut v = dl_interval(&toy, k / 4, &mut r5);
                let ex: HashSet<usize> = v.iter().map(|q| key(*q, p)).collect();
                rand_sym(&toy, k, &mut r5, &sq, &ex, &mut v);
                arms.push(("PC2".into(), v, true));
                let mut covs: HashMap<String, f64> = HashMap::new();
                for (name, fb, sym) in &arms {
                    assert_eq!(fb.len(), k, "arm {} size", name);
                    let dk: HashSet<usize> = fb.iter().map(|q| key(*q, p)).collect();
                    assert_eq!(dk.len(), k, "arm {} distinct", name);
                    let nsym = fb.iter().filter(|q| dk.contains(&key(toy.w.neg(**q), p))).count();
                    let mut r6 = Rng(seed_of(&[popid, bits as u64, ci as u64, 6]));
                    let c = coverage(&toy, fb, *sym, &sq, &mut r6);
                    let cv = c.cover as f64 / toy.n as f64;
                    covs.insert(name.clone(), cv);
                    if !first {
                        json.push(',');
                    }
                    first = false;
                    write!(
                        json,
                        "{{\"pop\":\"{}\",\"bits\":{},\"p\":{},\"curve\":{},\"N\":{},\"a\":{},\"b\":{},\"K\":{},\"arm\":\"{}\",\"neg_closed_frac\":{:.6},\"cover\":{},\"pairs\":{},\"C\":{:.8},\"collision_excess\":{:.8},\"semaev_checked\":{},\"semaev_disagree\":{},\"semaev_degenerate\":{}}}",
                        pop, bits, p, ci, toy.n, toy.w.a, toy.w.b, k, name, nsym as f64 / k as f64, c.cover, c.pairs, cv,
                        (c.pairs - c.cover) as f64 / c.pairs as f64, c.semaev_checked, c.semaev_disagree, c.semaev_degenerate
                    )
                    .unwrap();
                }
                for (name, _, sym) in &arms {
                    if name == "R-sym" || name == "R-asym" {
                        continue;
                    }
                    let ctrl = if *sym { "R-sym" } else { "R-asym" };
                    let rel = covs[name] / covs[ctrl] - 1.0;
                    if let Some(e) = rels.iter_mut().find(|e| e.0 == pop && e.1 == bits && &e.2 == name) {
                        e.4.push(rel);
                    } else {
                        rels.push((pop.to_string(), bits, name.clone(), ctrl.to_string(), vec![rel]));
                    }
                    // secondary: every arm vs R-asym too (asym arms already)
                    if *sym && name != "R-sym2" {
                        let rel2 = covs[name] / covs["R-asym"] - 1.0;
                        let nm = format!("{}|vsR-asym", name);
                        if let Some(e) = rels.iter_mut().find(|e| e.0 == pop && e.1 == bits && e.2 == nm) {
                            e.4.push(rel2);
                        } else {
                            rels.push((pop.to_string(), bits, nm, "R-asym".into(), vec![rel2]));
                        }
                    }
                }
                done += 1;
            }
            completed.push(format!("{{\"pop\":\"{}\",\"bits\":{},\"curves_completed\":{},\"curves_planned\":{},\"K\":{},\"p\":{}}}", pop, bits, done, ncurves, k, p));
            eprintln!("pop {} b {} done ({} curves) t={:.1}s", pop, bits, done, t0.elapsed().as_secs_f64());
        }
    }
    json.push_str("],\"summary\":[");
    let alpha_b = 0.05 / 21.0;
    for (i, (pop, bits, arm, ctrl, v)) in rels.iter().enumerate() {
        let seed = seed_of(&[7, *bits as u64, i as u64]);
        let (m, q) = bootstrap(v, seed, &[0.025, 0.975, alpha_b / 2.0, 1.0 - alpha_b / 2.0]);
        let verdict = if q[2] >= -0.01 && q[3] <= 0.01 {
            "indistinguishable"
        } else if (q[3] < -0.01) || (q[2] > 0.01) {
            "different"
        } else {
            "small_or_inconclusive"
        };
        let excl0 = q[2] > 0.0 || q[3] < 0.0;
        if i > 0 {
            json.push(',');
        }
        write!(
            json,
            "{{\"pop\":\"{}\",\"bits\":{},\"arm\":\"{}\",\"control\":\"{}\",\"n_curves\":{},\"mean_rel\":{:.6e},\"ci95\":[{:.6e},{:.6e}],\"ci_bonf\":[{:.6e},{:.6e}],\"bonf_excludes_0\":{},\"verdict\":\"{}\"}}",
            pop, bits, arm, ctrl, v.len(), m, q[0], q[1], q[2], q[3], excl0, verdict
        )
        .unwrap();
    }
    let _ = extra;
    write!(
        json,
        "],\"completed\":[{}],\"failures\":[{}],\"wall_s\":{:.1}}}",
        completed.join(","),
        failures.iter().map(|s| format!("\"{}\"", s)).collect::<Vec<_>>().join(","),
        t0.elapsed().as_secs_f64()
    )
    .unwrap();
    std::fs::write("results.json", json).unwrap();
    eprintln!("wrote results.json t={:.1}s", t0.elapsed().as_secs_f64());
}
