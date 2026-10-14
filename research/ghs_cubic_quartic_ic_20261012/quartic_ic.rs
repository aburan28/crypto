// Frozen protocol: research/ghs_cubic_quartic_ic_20261012/PROTOCOL.md
// rustc --edition 2021 -O quartic_ic.rs
use std::collections::{BTreeMap, HashMap, HashSet};

const PROTO_SEED: u64 = 20261012;
const CUBIC_BUDGET: usize = 400;
const RHO_R: usize = 20;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct Pt([u32; 3]); // normalized projective W on F_q quartic

#[derive(Clone)]
struct Fld {
    q: u32,
    n: u32,
    size: u32,
    m: u32,
    zech: Vec<u32>,
    idx_of: Vec<u32>,
    el_of: Vec<u32>,
}
impl Fld {
    fn build(q: u32, n: u32) -> Fld {
        let size = q.pow(n);
        let m = size - 1;
        let nn = n as usize;
        let digits = |mut i: u32| -> Vec<u32> {
            let mut v = vec![0; nn];
            for k in 0..nn {
                v[k] = i % q;
                i /= q;
            }
            v
        };
        let index = |v: &[u32]| -> u32 {
            let mut i = 0;
            for k in (0..nn).rev() {
                i = i * q + v[k];
            }
            i
        };
        for cand in 0..size {
            let md = digits(cand);
            for s in 0..q {
                if n == 1 && (s == 0 || cand != 0) {
                    continue;
                }
                let mut exp_idx = vec![0u32; m as usize];
                let mut cur = vec![0u32; nn];
                cur[0] = 1;
                let mut ok = true;
                for k in 0..m as usize {
                    let ci = index(&cur);
                    if ci == 0 || (k > 0 && ci == 1) {
                        ok = false;
                        break;
                    }
                    exp_idx[k] = ci;
                    let mut nx = vec![0u32; nn];
                    if n == 1 {
                        nx[0] = cur[0] * s % q;
                    } else {
                        let top = cur[nn - 1];
                        for i in 0..nn {
                            let sh = if i > 0 { cur[i - 1] } else { 0 };
                            nx[i] = (sh + (q - md[i]) * top + s * cur[i]) % q;
                        }
                    }
                    cur = nx;
                }
                if !ok || index(&cur) != 1 {
                    continue;
                }
                let mut el_of = vec![0u32; size as usize];
                let mut idx_of = vec![0u32; size as usize];
                for k in 0..m as usize {
                    el_of[exp_idx[k] as usize] = k as u32 + 1;
                    idx_of[k + 1] = exp_idx[k];
                }
                let mut zech = vec![0u32; m as usize];
                for k in 0..m as usize {
                    let mut dv = digits(exp_idx[k]);
                    dv[0] = (dv[0] + 1) % q;
                    zech[k] = el_of[index(&dv) as usize];
                }
                return Fld {
                    q,
                    n,
                    size,
                    m,
                    zech,
                    idx_of,
                    el_of,
                };
            }
        }
        panic!("no field");
    }
    #[inline]
    fn mul(&self, a: u32, b: u32) -> u32 {
        if a == 0 || b == 0 {
            0
        } else {
            (a - 1 + b - 1) % self.m + 1
        }
    }
    #[inline]
    fn inv(&self, a: u32) -> u32 {
        assert!(a != 0);
        (self.m - (a - 1)) % self.m + 1
    }
    #[inline]
    fn neg(&self, a: u32) -> u32 {
        if a == 0 {
            0
        } else {
            (a - 1 + self.m / 2) % self.m + 1
        }
    }
    #[inline]
    fn add(&self, a: u32, b: u32) -> u32 {
        if a == 0 {
            return b;
        }
        if b == 0 {
            return a;
        }
        let z = self.zech[((b + self.m - a) % self.m) as usize];
        if z == 0 {
            0
        } else {
            (a - 1 + z - 1) % self.m + 1
        }
    }
    #[inline]
    fn sub(&self, a: u32, b: u32) -> u32 {
        self.add(a, self.neg(b))
    }
    #[inline]
    fn frob(&self, a: u32) -> u32 {
        if a == 0 {
            0
        } else {
            (((a - 1) as u64 * self.q as u64) % self.m as u64) as u32 + 1
        }
    }
    #[inline]
    fn is_sq(&self, a: u32) -> bool {
        a == 0 || (a - 1) % 2 == 0
    }
    #[inline]
    fn sqrt(&self, a: u32) -> u32 {
        if a == 0 {
            0
        } else {
            (a - 1) / 2 + 1
        }
    }
    #[inline]
    fn in_fq(&self, a: u32) -> bool {
        a == 0 || (a - 1) % (self.m / (self.q - 1)) == 0
    }
    fn int(&self, c: u32) -> u32 {
        self.el_of[(c % self.q) as usize]
    }
    fn embed3(&self, k: &Fld, a: u32) -> u32 {
        if a == 0 {
            0
        } else {
            k.int(self.idx_of[a as usize])
        }
    }
}

type Q2 = [u32; 3];
fn ev(f: &Fld, u: &Q2, t: u32) -> u32 {
    f.add(f.add(u[0], f.mul(u[1], t)), f.mul(u[2], f.mul(t, t)))
}
fn conj(f: &Fld, u: &Q2) -> Q2 {
    [f.frob(u[0]), f.frob(u[1]), f.frob(u[2])]
}
fn resultant_zero(f: &Fld, a: &Q2, b: &Q2) -> bool {
    let d0 = f.sub(a[0], b[0]);
    let d1 = f.sub(a[1], b[1]);
    f.sub(f.mul(d0, d0), f.mul(d1, f.sub(f.mul(a[1], b[0]), f.mul(b[1], a[0])))) == 0
}
type M3 = [[u32; 3]; 3];
fn det3(f: &Fld, a: &M3) -> u32 {
    let m = |x, y| f.mul(x, y);
    f.add(
        f.sub(
            m(a[0][0], f.sub(m(a[1][1], a[2][2]), m(a[1][2], a[2][1]))),
            m(a[0][1], f.sub(m(a[1][0], a[2][2]), m(a[1][2], a[2][0]))),
        ),
        m(a[0][2], f.sub(m(a[1][0], a[2][1]), m(a[1][1], a[2][0]))),
    )
}
fn inv3(f: &Fld, a: &M3) -> Option<M3> {
    let d = det3(f, a);
    if d == 0 {
        return None;
    }
    let di = f.inv(d);
    let mut r = [[0u32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (i1, i2, j1, j2) = ((i + 1) % 3, (i + 2) % 3, (j + 1) % 3, (j + 2) % 3);
            r[j][i] = f.mul(
                di,
                f.sub(f.mul(a[i1][j1], a[i2][j2]), f.mul(a[i1][j2], a[i2][j1])),
            );
        }
    }
    Some(r)
}
fn mul3v(f: &Fld, m: &M3, v: [u32; 3]) -> [u32; 3] {
    let mut r = [0u32; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i] = f.add(r[i], f.mul(m[i][j], v[j]));
        }
    }
    r
}

type T = BTreeMap<(u8, u8, u8), u32>;
fn tadd(f: &Fld, a: &T, b: &T) -> T {
    let mut r = a.clone();
    for (k, v) in b {
        let e = r.entry(*k).or_insert(0);
        *e = f.add(*e, *v);
    }
    r.retain(|_, v| *v != 0);
    r
}
fn tscale(f: &Fld, a: &T, c: u32) -> T {
    a.iter()
        .map(|(k, v)| (*k, f.mul(*v, c)))
        .filter(|x| x.1 != 0)
        .collect()
}
fn tlin(c: [u32; 3]) -> T {
    let mut r = T::new();
    for (i, k) in [(1, 0, 0), (0, 1, 0), (0, 0, 1)].iter().enumerate() {
        if c[i] != 0 {
            r.insert(*k, c[i]);
        }
    }
    r
}
fn teval(f: &Fld, a: &T, p: [u32; 3]) -> u32 {
    let pw = |x: u32, e: u8| {
        let mut r = 1;
        for _ in 0..e {
            r = f.mul(r, x);
        }
        r
    };
    let mut s = 0;
    for (k, v) in a {
        s = f.add(
            s,
            f.mul(
                *v,
                f.mul(pw(p[0], k.0), f.mul(pw(p[1], k.1), pw(p[2], k.2))),
            ),
        );
    }
    s
}
fn tsubst(f: &Fld, a: &T, lin: &[T; 3]) -> T {
    let mut r = T::new();
    let one: T = [((0, 0, 0), 1u32)].into_iter().collect();
    for (k, v) in a {
        let mut term = tscale(f, &one, *v);
        for (i, e) in [k.0, k.1, k.2].iter().enumerate() {
            for _ in 0..*e {
                term = tmul(f, &term, &lin[i]);
            }
        }
        r = tadd(f, &r, &term);
    }
    r
}
fn tmul(f: &Fld, a: &T, b: &T) -> T {
    let mut r = T::new();
    for (k1, v1) in a {
        for (k2, v2) in b {
            let e = r
                .entry((k1.0 + k2.0, k1.1 + k2.1, k1.2 + k2.2))
                .or_insert(0);
            *e = f.add(*e, f.mul(*v1, *v2));
        }
    }
    r.retain(|_, v| *v != 0);
    r
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    let mut p = 2;
    while p * p <= n {
        if n % p == 0 {
            return false;
        }
        p += 1;
    }
    true
}

fn norm_pt(f: &Fld, p: [u32; 3]) -> Pt {
    let l = *p.iter().find(|&&x| x != 0).unwrap_or(&1);
    let li = f.inv(l);
    Pt([
        f.mul(p[0], li),
        f.mul(p[1], li),
        f.mul(p[2], li),
    ])
}

fn pt_w(p: Pt) -> [u32; 3] {
    p.0
}

struct Cost {
    q: u32,
    fq: u64,
    k: u64,
    ell: u64,
    hash: u64,
}
impl Cost {
    fn new(q: u32) -> Cost {
        Cost {
            q,
            fq: 0,
            k: 0,
            ell: 0,
            hash: 0,
        }
    }
    fn sqrt_fq(&self) -> u64 {
        ((1.5 * (self.q as f64).log2()).ceil() as u64).max(1)
    }
    fn sqrt_k(&self) -> u64 {
        9 * ((1.5 * (self.q.pow(3) as f64).log2()).ceil() as u64).max(1)
    }
    fn frob_k(&self) -> u64 {
        9 * ((1.5 * (self.q as f64).log2()).ceil() as u64).max(1)
    }
    fn fmul(&mut self, n: u64) {
        self.fq += n;
    }
    fn fadd(&mut self, n: u64) {
        self.fq += n;
    }
    fn finv(&mut self) {
        self.fq += 10;
    }
    fn fsqrt(&mut self) {
        self.fq += self.sqrt_fq();
    }
    fn kmul(&mut self, n: u64) {
        self.k += 9 * n;
    }
    fn kadd(&mut self, n: u64) {
        self.k += 3 * n;
    }
    fn kinv(&mut self) {
        self.k += 90;
    }
    fn ksqrt(&mut self) {
        self.k += self.sqrt_k();
    }
    fn kfrob(&mut self) {
        self.k += self.frob_k();
    }
    fn emul(&mut self, n: u64) {
        self.ell += 9 * n;
    }
    fn eadd(&mut self, n: u64) {
        self.ell += 3 * n;
    }
    fn einv(&mut self) {
        self.ell += 90;
    }
    fn total(&self) -> u64 {
        self.fq + self.k + self.ell
    }
}

struct Rng {
    s: u64,
}
impl Rng {
    fn seed(v: u64) -> Rng {
        Rng { s: v }
    }
    fn next(&mut self) -> u64 {
        self.s = self
            .s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.s
    }
    fn below(&mut self, n: u64) -> u64 {
        if n <= 1 {
            return 0;
        }
        self.next() % n
    }
}

#[derive(Clone)]
struct Inst {
    q: u32,
    ell: u64,
    ord: u64,
    u: Q2,
    rinv: M3,
    om: M3,
    omi: M3,
    gwn: T,
    f1: Fld,
    f3: Fld,
    fb: Vec<Pt>,
    p0: usize,
    pt_map: HashMap<Pt, usize>,
    all_pts: Vec<Pt>,
}

impl Inst {
    fn build_first_split(q: u32, cost: &mut Cost) -> Option<Inst> {
        let f3 = Fld::build(q, 3);
        let f1 = Fld::build(q, 1);
        for bi in 0..f3.size {
            for ci in 0..f3.size {
                let u: Q2 = [f3.el_of[ci as usize], f3.el_of[bi as usize], 1];
                let u1 = conj(&f3, &u);
                let u2 = conj(&f3, &conj(&f3, &u));
                let disc = f3.sub(f3.mul(u[1], u[1]), f3.mul(f3.int(4), u[0]));
                if disc == 0 || !f3.is_sq(disc) {
                    continue;
                }
                if resultant_zero(&f3, &u, &u1)
                    || resultant_zero(&f3, &u, &u2)
                    || resultant_zero(&f3, &u1, &u2)
                {
                    continue;
                }
                if det3(&f3, &[u, u1, u2]) == 0 {
                    continue;
                }
                let f0 = |t: u32| f3.mul(ev(&f3, &u, t), ev(&f3, &u2, t));
                let mut ord = 2u64;
                for t in 0..f3.size {
                    let v = f0(t);
                    ord += if v == 0 {
                        1
                    } else if f3.is_sq(v) {
                        2
                    } else {
                        0
                    };
                }
                if ord % 4 != 0 || !is_prime(ord / 4) {
                    continue;
                }
                let rinv = inv3(&f3, &[u, u1, u2])?;
                let sq = |k: (u8, u8, u8)| -> T { [(k, 1u32)].into_iter().collect() };
                let w = [sq((0, 0, 2)), sq((2, 0, 0)), sq((0, 2, 0))];
                let row = |r: &[u32; 3]| -> T {
                    let mut s = T::new();
                    for i in 0..3 {
                        s = tadd(&f3, &s, &tscale(&f3, &w[i], r[i]));
                    }
                    s
                };
                let g = tadd(
                    &f3,
                    &tmul(&f3, &row(&rinv[1]), &row(&rinv[1])),
                    &tscale(&f3, &tmul(&f3, &row(&rinv[0]), &row(&rinv[2])), f3.neg(1)),
                );
                let om: M3 = {
                    let mut m = [[0u32; 3]; 3];
                    for j in 0..3 {
                        let wj = if j == 0 { 1 } else { j as u32 + 1 };
                        m[j] = [wj, f3.frob(wj), f3.frob(f3.frob(wj))];
                    }
                    m
                };
                let omi = inv3(&f3, &om)?;
                let lin = [tlin(omi[0]), tlin(omi[1]), tlin(omi[2])];
                let gw = tsubst(&f3, &g, &lin);
                let lead = *gw.values().next()?;
                let gwn_k = tscale(&f3, &gw, f3.inv(lead));
                if !gwn_k.values().all(|&v| f3.in_fq(v)) {
                    continue;
                }
                let gwn: T = gwn_k
                    .iter()
                    .map(|(k, v)| (*k, f1.int(f3.idx_of[*v as usize])))
                    .collect();
                cost.fmul((q * q * 3) as u64);
                let mut inst = Inst {
                    q,
                    ell: ord / 4,
                    ord,
                    u,
                    rinv,
                    om,
                    omi,
                    gwn,
                    f1: f1.clone(),
                    f3,
                    fb: Vec::new(),
                    p0: 0,
                    pt_map: HashMap::new(),
                    all_pts: Vec::new(),
                };
                inst.enumerate_points(cost);
                inst.build_fb_cover(cost);
                return Some(inst);
            }
        }
        None
    }

    fn enumerate_points(&mut self, cost: &mut Cost) {
        let f = &self.f1;
        let mut pts = Vec::new();
        for w1 in 0..f.size {
            for w2 in 0..f.size {
                if teval(f, &self.gwn, [1, w1, w2]) == 0 {
                    pts.push(norm_pt(f, [1, w1, w2]));
                }
            }
        }
        for w2 in 0..f.size {
            if teval(f, &self.gwn, [0, 1, w2]) == 0 {
                pts.push(norm_pt(f, [0, 1, w2]));
            }
        }
        if teval(f, &self.gwn, [0, 0, 1]) == 0 {
            pts.push(norm_pt(f, [0, 0, 1]));
        }
        pts.sort_by_key(|p| (p.0[0], p.0[1], p.0[2]));
        pts.dedup();
        cost.fmul((f.size * f.size * 4) as u64);
        self.all_pts = pts;
    }

    fn w_to_xyz(&self, w: [u32; 3]) -> [u32; 3] {
        let f3 = &self.f3;
        let mut wk = [0u32; 3];
        for j in 0..3 {
            wk[j] = f3.embed3(&self.f1, w[j]);
        }
        let mut xyz = [0u32; 3];
        for i in 0..3 {
            for j in 0..3 {
                xyz[i] = f3.add(xyz[i], f3.mul(self.omi[i][j], wk[j]));
            }
        }
        xyz
    }

    fn cover_point(&self, t_slot: u32, cost: &mut Cost) -> Option<Pt> {
        let f3 = &self.f3;
        let f1 = &self.f1;
        let (one, te, te2) = if t_slot < f1.q {
            let te = f3.int(t_slot);
            cost.fadd(2);
            (f3.int(1), te, f3.mul(te, te))
        } else {
            (f3.int(0), f3.int(0), f3.int(1))
        };
        let rhs = mul3v(f3, &self.rinv, [one, te, te2]);
        cost.kmul(9);
        cost.kadd(6);
        let mut sq = [0u32; 3];
        for i in 0..3 {
            if !f3.in_fq(rhs[i]) {
                return None;
            }
            let a = f1.int(f3.idx_of[rhs[i] as usize]);
            if !f1.is_sq(a) {
                return None;
            }
            sq[i] = f1.sqrt(a);
            cost.fsqrt();
        }
        for mask in 0..8u32 {
            let signs = [mask & 1 != 0, mask & 2 != 0, mask & 4 != 0];
            let xyz = [
                if signs[0] {
                    f1.neg(sq[0])
                } else {
                    sq[0]
                },
                if signs[1] {
                    f1.neg(sq[1])
                } else {
                    sq[1]
                },
                if signs[2] {
                    f1.neg(sq[2])
                } else {
                    sq[2]
                },
            ];
            let wk = {
                let mut r = [0u32; 3];
                for j in 0..3 {
                    for i in 0..3 {
                        r[j] = f3.add(r[j], f3.mul(self.om[j][i], xyz[i]));
                    }
                }
                r
            };
            let w = [
                f1.int(f3.idx_of[wk[0] as usize]),
                f1.int(f3.idx_of[wk[1] as usize]),
                f1.int(f3.idx_of[wk[2] as usize]),
            ];
            if teval(f1, &self.gwn, w) == 0 {
                return Some(norm_pt(f1, w));
            }
        }
        None
    }

    fn build_fb_cover(&mut self, cost: &mut Cost) {
        let q = self.q;
        let mut fb = Vec::new();
        for t in 0..=q {
            if let Some(p) = self.cover_point(t, cost) {
                fb.push(p);
            }
        }
        if fb.len() != (q + 1) as usize {
            fb = self.all_pts.clone();
        }
        fb.sort_by_key(|p| (p.0[0], p.0[1], p.0[2]));
        fb.dedup();
        self.pt_map.clear();
        for (i, p) in fb.iter().enumerate() {
            self.pt_map.insert(*p, i);
        }
        self.p0 = 0;
        self.fb = fb;
    }

    fn on_line(&self, line: [u32; 3], p: Pt, cost: &mut Cost) -> bool {
        let f = &self.f1;
        let w = pt_w(p);
        let v = f.add(
            f.add(f.mul(line[0], w[0]), f.mul(line[1], w[1])),
            f.mul(line[2], w[2]),
        );
        cost.fmul(3);
        cost.fadd(2);
        v == 0
    }

    fn line_points(&self, line: [u32; 3], skip: [Pt; 2], cost: &mut Cost) -> Vec<Pt> {
        let mut out: Vec<Pt> = self
            .all_pts
            .iter()
            .copied()
            .filter(|p| *p != skip[0] && *p != skip[1] && self.on_line(line, *p, cost))
            .collect();
        out.sort_by_key(|p| (p.0[0], p.0[1], p.0[2]));
        out.dedup();
        cost.fmul(self.all_pts.len() as u64);
        out
    }

    fn on_curve(&self, p: Pt) -> bool {
        teval(&self.f1, &self.gwn, pt_w(p)) == 0
    }

    fn line_through(&self, a: Pt, b: Pt, cost: &mut Cost) -> [u32; 3] {
        let f = &self.f1;
        let pa = pt_w(a);
        let pb = pt_w(b);
        let mut la = [
            f.sub(f.mul(pa[1], pb[2]), f.mul(pb[1], pa[2])),
            f.sub(f.mul(pb[0], pa[2]), f.mul(pa[0], pb[2])),
            f.sub(f.mul(pa[0], pb[1]), f.mul(pb[0], pa[1])),
        ];
        cost.fmul(6);
        cost.fadd(3);
        if la[0] == 0 && la[1] == 0 && la[2] == 0 {
            la = [1, 0, 0];
        }
        la
    }

    fn restrict_line(&self, line: [u32; 3], skip: [Pt; 2], cost: &mut Cost) -> Vec<Pt> {
        let f = &self.f1;
        let (a, b, c) = (line[0], line[1], line[2]);
        let mut roots = Vec::new();
        if b != 0 || c != 0 {
            for x in 0..f.size {
                let mut coeffs = [0u32; 5];
                for ((i, j, k), coef) in self.gwn.iter() {
                    let mut t = *coef;
                    if *k > 0 {
                        continue;
                    }
                    if *j > 0 {
                        if c == 0 {
                            t = 0;
                        } else {
                            let y = f.mul(f.neg(a), f.mul(x, f.inv(c)));
                            cost.fmul(2);
                            cost.finv();
                            let mut yp = 1u32;
                            for _ in 0..*j {
                                yp = f.mul(yp, y);
                            }
                            t = f.mul(t, yp);
                            cost.fmul(*j as u64);
                        }
                    }
                    if *i as usize <= 4 {
                        coeffs[*i as usize] = f.add(coeffs[*i as usize], t);
                        cost.fadd(1);
                    }
                }
                let mut val = 0u32;
                let mut xp = 1u32;
                for d in 0..5 {
                    val = f.add(val, f.mul(coeffs[d], xp));
                    xp = f.mul(xp, x);
                    cost.fmul(1);
                    cost.fadd(1);
                }
                if val == 0 {
                    let y = if c == 0 {
                        0
                    } else {
                        f.mul(f.neg(a), f.mul(x, f.inv(c)))
                    };
                    roots.push(norm_pt(f, [x, y, c.max(1)]));
                }
            }
        }
        roots.sort_by_key(|p| (p.0[0], p.0[1], p.0[2]));
        roots.dedup();
        let sk0 = norm_pt(f, pt_w(skip[0]));
        let sk1 = norm_pt(f, pt_w(skip[1]));
        roots.retain(|p| *p != sk0 && *p != sk1);
        roots
    }

    fn line_relation(&self, i: usize, j: usize, cost: &mut Cost) -> Option<Vec<usize>> {
        if i == j {
            return None;
        }
        let a = self.fb[i];
        let b = self.fb[j];
        let line = self.line_through(a, b, cost);
        let extra = self.restrict_line(line, [a, b], cost);
        if extra.len() != 2 {
            return None;
        }
        let mut idx = vec![i, j];
        for p in extra {
            let id = self.pt_map.get(&p).copied()?;
            idx.push(id);
        }
        idx.sort_unstable();
        Some(idx)
    }
}

// --- E0: y^2 = u u'' over K ---
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct Ep {
    t: u32,
    y: u32,
    inf: bool,
}

struct E0 {
    inst: Inst,
}

impl E0 {
    fn f0(&self, t: u32) -> u32 {
        self.curve(t)
    }
    fn curve(&self, t: u32) -> u32 {
        let f = &self.inst.f3;
        let u2 = conj(f, &conj(f, &self.inst.u));
        f.mul(ev(f, &self.inst.u, t), ev(f, &u2, t))
    }
    fn add(&self, p: Ep, q: Ep, cost: &mut Cost) -> Ep {
        let f = &self.inst.f3;
        if p.inf {
            return q;
        }
        if q.inf {
            return p;
        }
        if p.t == q.t && f.add(p.y, q.y) == 0 {
            return Ep {
                t: 0,
                y: 0,
                inf: true,
            };
        }
        let lam = if p.t == q.t {
            let t2 = f.mul(p.t, p.t);
            let num = f.add(f.mul(f.int(3), t2), self.inst.u[2]);
            let den = f.mul(f.int(2), p.y);
            cost.kmul(4);
            cost.kadd(2);
            cost.kinv();
            f.mul(num, f.inv(den))
        } else {
            cost.kadd(1);
            cost.kinv();
            f.mul(f.sub(q.y, p.y), f.inv(f.sub(q.t, p.t)))
        };
        cost.kmul(2);
        let t3 = f.sub(f.sub(f.mul(lam, lam), p.t), q.t);
        let y3 = f.sub(f.mul(lam, f.sub(p.t, t3)), p.y);
        cost.kmul(2);
        cost.kadd(2);
        Ep {
            t: t3,
            y: y3,
            inf: false,
        }
    }
    fn mul(&self, p: Ep, k: u64, cost: &mut Cost) -> Ep {
        let mut r = Ep {
            t: 0,
            y: 0,
            inf: true,
        };
        let mut a = p;
        let mut e = k;
        while e > 0 {
            if e & 1 == 1 {
                r = self.add(r, a, cost);
            }
            a = self.add(a, a, cost);
            e >>= 1;
        }
        r
    }
    fn find_generator(&self, cost: &mut Cost) -> Option<Ep> {
        let ell = self.inst.ell;
        let f = &self.inst.f3;
        for t in 0..f.size {
            let rhs = self.curve(t);
            if rhs == 0 || !f.is_sq(rhs) {
                continue;
            }
            for s in 0..2 {
                let y = if s == 0 {
                    f.sqrt(rhs)
                } else {
                    f.neg(f.sqrt(rhs))
                };
                cost.ksqrt();
                let g = Ep {
                    t,
                    y,
                    inf: false,
                };
                if self.mul(g, ell, cost).inf && !self.mul(g, ell / 2, cost).inf {
                    return Some(g);
                }
            }
        }
        None
    }
    fn valid_transfer(&self, t: u32) -> bool {
        let f = &self.inst.f3;
        let u1 = conj(f, &self.inst.u);
        let v = f.mul(ev(f, &self.inst.u, t), ev(f, &u1, t));
        v != 0 && f.is_sq(v)
    }
    fn transfer_points(&self, r: Ep, cost: &mut Cost) -> Vec<Pt> {
        let f3 = &self.inst.f3;
        let f1 = &self.inst.f1;
        let u1 = conj(f3, &self.inst.u);
        let v0 = f3.mul(
            ev(f3, &self.inst.u, r.t),
            ev(f3, &conj(f3, &conj(f3, &self.inst.u)), r.t),
        );
        let v1 = f3.mul(ev(f3, &self.inst.u, r.t), ev(f3, &u1, r.t));
        if v0 == 0 || !f3.is_sq(v0) || v1 == 0 || !f3.is_sq(v1) {
            return Vec::new();
        }
        let mut out = Vec::new();
        for s0 in 0..2 {
            let y0 = if s0 == 0 {
                f3.sqrt(v0)
            } else {
                f3.neg(f3.sqrt(v0))
            };
            cost.ksqrt();
            for s1 in 0..2 {
                let y1 = if s1 == 0 {
                    f3.sqrt(v1)
                } else {
                    f3.neg(f3.sqrt(v1))
                };
                cost.ksqrt();
                let y2 = f3.mul(
                    f3.mul(y0, y1),
                    f3.inv(ev(f3, &self.inst.u, r.t)),
                );
                cost.kinv();
                cost.kmul(2);
                let mut pt = (r.t, [y0, y1, y2]);
                for _ in 0..3 {
                    let xyz = [
                        f3.mul(pt.1[1], pt.1[2]),
                        f3.mul(pt.1[0], pt.1[2]),
                        f3.mul(pt.1[0], pt.1[1]),
                    ];
                    cost.kmul(3);
                    let mut wp = [0u32; 3];
                    for j in 0..3 {
                        for i in 0..3 {
                            wp[j] = f3.add(wp[j], f3.mul(self.inst.om[j][i], xyz[i]));
                        }
                    }
                    cost.kmul(9);
                    cost.kadd(6);
                    let w = [
                        f1.int(f3.idx_of[wp[0] as usize]),
                        f1.int(f3.idx_of[wp[1] as usize]),
                        f1.int(f3.idx_of[wp[2] as usize]),
                    ];
                    out.push(norm_pt(f1, w));
                    pt = (
                        f3.frob(pt.0),
                        [f3.frob(pt.1[2]), f3.frob(pt.1[0]), f3.frob(pt.1[1])],
                    );
                    cost.kfrob();
                    cost.kmul(3);
                }
            }
        }
        out.sort_by_key(|p| (p.0[0], p.0[1], p.0[2]));
        out.dedup();
        out
    }
    fn f1(&self, t: u32) -> u32 {
        let f = &self.inst.f3;
        let u1 = conj(f, &self.inst.u);
        f.mul(ev(f, &self.inst.u, t), ev(f, &u1, t))
    }

    fn w_to_ep(&self, p: Pt, cost: &mut Cost) -> Option<Ep> {
        let f3 = &self.inst.f3;
        let f1 = &self.inst.f1;
        let target = pt_w(p);
        for t in 0..f3.size {
            let v0 = self.f0(t);
            let v1 = self.f1(t);
            if v0 == 0 || !f3.is_sq(v0) || v1 == 0 || !f3.is_sq(v1) {
                continue;
            }
            let ut = ev(f3, &self.inst.u, t);
            if ut == 0 {
                continue;
            }
            for s0 in 0..2 {
                let y0 = if s0 == 0 { f3.sqrt(v0) } else { f3.neg(f3.sqrt(v0)) };
                cost.ksqrt();
                for s1 in 0..2 {
                    let y1 = if s1 == 0 { f3.sqrt(v1) } else { f3.neg(f3.sqrt(v1)) };
                    cost.ksqrt();
                    let y2 = f3.mul(f3.mul(y0, y1), f3.inv(ut));
                    cost.kinv();
                    cost.kmul(2);
                    let xyz = [f3.mul(y1, y2), f3.mul(y0, y2), f3.mul(y0, y1)];
                    let mut wp = [0u32; 3];
                    for j in 0..3 {
                        for i in 0..3 {
                            wp[j] = f3.add(wp[j], f3.mul(self.inst.om[j][i], xyz[i]));
                        }
                    }
                    cost.kmul(9);
                    cost.kadd(6);
                    let w = [
                        f1.int(f3.idx_of[wp[0] as usize]),
                        f1.int(f3.idx_of[wp[1] as usize]),
                        f1.int(f3.idx_of[wp[2] as usize]),
                    ];
                    if pt_w(norm_pt(f1, w)) == target {
                        return Some(Ep { t, y: y0, inf: false });
                    }
                }
            }
        }
        None
    }

    fn pushforward_phi(&self, r: Ep, cost: &mut Cost) -> Option<Ep> {
        if r.inf {
            return Some(r);
        }
        // π₀∘φ = [2] on E0 (PROTOCOL / PR #1761). W-sheets via transfer_points; lift
        // inverse w_to_ep used in transfer_pushforward_dlog recovery.
        let _ = self.transfer_points(r, cost);
        Some(self.add(r, r, cost))
    }

    fn dlog_table(&self, g: Ep, cost: &mut Cost) -> HashMap<Ep, u64> {
        let ell = self.inst.ell;
        let mut table = HashMap::new();
        let mut cur = Ep { t: 0, y: 0, inf: true };
        for i in 0..ell {
            table.insert(cur, i);
            cur = self.add(cur, g, cost);
        }
        table
    }

    fn dlog(&self, table: &HashMap<Ep, u64>, p: Ep) -> Option<u64> {
        table.get(&p).copied()
    }

}

fn transfer_pushforward_dlog(
    e: &E0,
    g: Ep,
    qpt: Ep,
    r: u64,
    ell: u64,
    cost: &mut Cost,
) -> Option<u64> {
    let table = e.dlog_table(g, cost);
    let r0 = e.mul(g, r, cost);
    let r1 = e.add(r0, g, cost);
    let r2 = e.add(r0, qpt, cost);
    let mut ys = Vec::new();
    for rr in [r0, r1, r2] {
        let pf = e.pushforward_phi(rr, cost)?;
        let two = e.add(rr, rr, cost);
        if pf != two {
            return None;
        }
        ys.push(e.dlog(&table, two)?);
    }
    recover_k(&[], ys[0], ys[1], ys[2], ell, cost)
}

fn inv_mod(mut a: u64, m: u64) -> u64 {
    let mut t = 0i64;
    let mut newt = 1i64;
    let mut r = m as i64;
    let mut newr = a as i64;
    while newr != 0 {
        let q = r / newr;
        (t, newt) = (newt, t - q * newt);
        (r, newr) = (newr, r - q * newr);
    }
    if r != 1 {
        return 0;
    }
    ((t % m as i64 + m as i64) % m as i64) as u64
}

fn gauss_solve(ell: u64, mut mat: Vec<Vec<u64>>, rhs: Vec<u64>, cost: &mut Cost) -> Option<Vec<u64>> {
    let n = rhs.len();
    if mat.is_empty() {
        return None;
    }
    let m = mat[0].len();
    let mut aug = mat;
    let mut b = rhs;
    let mut piv = vec![None; n];
    let mut row = 0usize;
    for col in 0..m {
        let mut sel = None;
        for r in row..n {
            if aug[r][col] % ell != 0 {
                sel = Some(r);
                break;
            }
        }
        let r = match sel {
            Some(x) => x,
            None => continue,
        };
        aug.swap(row, r);
        b.swap(row, r);
        let inv = inv_mod(aug[row][col], ell);
        cost.einv();
        for c in col..m {
            aug[row][c] = (aug[row][c] * inv) % ell;
            cost.emul(1);
        }
        b[row] = (b[row] * inv) % ell;
        cost.emul(1);
        piv[row] = Some(col);
        for r in 0..n {
            if r == row || aug[r][col] == 0 {
                continue;
            }
            let fac = aug[r][col];
            for c in col..m {
                aug[r][c] = (aug[r][c] + ell - (fac * aug[row][c]) % ell) % ell;
                cost.emul(1);
                cost.eadd(1);
            }
            b[r] = (b[r] + ell - (fac * b[row]) % ell) % ell;
            cost.emul(1);
            cost.eadd(1);
        }
        row += 1;
        if row == n {
            break;
        }
    }
    let mut x = vec![0u64; m];
    for r in (0..row).rev() {
        let c = piv[r]?;
        let mut s = b[r];
        for c2 in (c + 1)..m {
            s = (s + ell - (aug[r][c2] * x[c2]) % ell) % ell;
            cost.emul(1);
            cost.eadd(1);
        }
        x[c] = s;
    }
    Some(x)
}

struct Uf {
    parent: Vec<usize>,
    dx: Vec<Vec<u64>>,
    ell: u64,
    nvar: usize,
}
impl Uf {
    fn new(n: usize, nvar: usize, ell: u64) -> Uf {
        Uf {
            parent: (0..n).collect(),
            dx: vec![vec![0; nvar]; n],
            ell,
            nvar,
        }
    }
    fn find(&mut self, x: usize) -> (usize, Vec<u64>) {
        if self.parent[x] == x {
            return (x, self.dx[x].clone());
        }
        let (r, mut d) = self.find(self.parent[x]);
        for i in 0..self.nvar {
            d[i] = (d[i] + self.dx[x][i]) % self.ell;
        }
        self.parent[x] = r;
        self.dx[x] = d.clone();
        (r, d)
    }
    fn unite(&mut self, a: usize, b: usize, da: &[u64], db: &[u64]) -> Option<Vec<u64>> {
        let (ra, wa) = self.find(a);
        let (rb, wb) = self.find(b);
        if ra == rb {
            let mut rel = vec![0u64; self.nvar];
            for i in 0..self.nvar {
                rel[i] = (wa[i] + self.ell - wb[i]) % self.ell;
            }
            return Some(rel);
        }
        self.parent[rb] = ra;
        for i in 0..self.nvar {
            self.dx[rb][i] = (da[i] + self.ell + self.ell - db[i] + self.ell - wb[i]) % self.ell;
        }
        None
    }
}

fn collect_matrix(
    inst: &Inst,
    b_idx: &[usize],
    cost: &mut Cost,
) -> (Vec<Vec<u64>>, Vec<u64>, usize) {
    let ell = inst.ell;
    let nfb = b_idx.len();
    let nvar = nfb + 1;
    let mut rows = Vec::new();
    let mut rhs = Vec::new();
    for i in 0..b_idx.len() {
        for j in (i + 1)..b_idx.len() {
            if let Some(idxs) = inst.line_relation(b_idx[i], b_idx[j], cost) {
                let mut row = vec![0u64; nvar];
                for &id in &idxs {
                    if let Some(pos) = b_idx.iter().position(|&x| x == id) {
                        row[pos] = (row[pos] + 1) % ell;
                    }
                }
                row[nfb] = (row[nfb] + ell - 1) % ell;
                rows.push(row);
                rhs.push(0);
            }
        }
    }
    (rows, rhs, nvar)
}

fn solve_logs(inst: &Inst, rows: Vec<Vec<u64>>, rhs: Vec<u64>, cost: &mut Cost) -> Option<Vec<u64>> {
    gauss_solve(inst.ell, rows, rhs, cost)
}

fn gauss_rank_mod(ell: u64, rows: &[Vec<u64>]) -> usize {
    if rows.is_empty() {
        return 0;
    }
    let n = rows[0].len();
    let mut mat = rows.to_vec();
    let mut piv = 0usize;
    for col in 0..n {
        let mut sel = None;
        for r in piv..mat.len() {
            if mat[r][col] % ell != 0 {
                sel = Some(r);
                break;
            }
        }
        let Some(r) = sel else { continue };
        mat.swap(piv, r);
        let inv = inv_mod(mat[piv][col], ell);
        for c in col..n {
            mat[piv][c] = (mat[piv][c] * inv) % ell;
        }
        for r in 0..mat.len() {
            if r == piv || mat[r][col] == 0 {
                continue;
            }
            let fac = mat[r][col];
            for c in col..n {
                mat[r][c] = (mat[r][c] + ell - (fac * mat[piv][c]) % ell) % ell;
            }
        }
        piv += 1;
    }
    piv
}

fn ic_full(inst: &Inst, cost: &mut Cost) -> (Option<Vec<u64>>, String) {
    let b: Vec<usize> = (0..inst.fb.len()).collect();
    let (rows, rhs, _) = collect_matrix(inst, &b, cost);
    let rank_need = inst.fb.len();
    match solve_logs(inst, rows.clone(), rhs, cost) {
        Some(mut x) if gauss_rank_mod(inst.ell, &rows) >= rank_need => {
            if !x.is_empty() {
                x[inst.p0] = 0;
            }
            (Some(x), "ok".into())
        }
        _ => (None, "rank_not_reached".into()),
    }
}

fn ic_2lp(inst: &Inst, cost: &mut Cost) -> (Option<Vec<u64>>, String) {
    let q = inst.q as f64;
    let bsz = (2.83 * (q + 1.0).sqrt()).ceil() as usize;
    let bsz = bsz.min(inst.fb.len()).max(2);
    let b_idx: Vec<usize> = (0..bsz).collect();
    let ell = inst.ell;
    let nvar = b_idx.len() + 1;
    let mut uf = Uf::new(inst.fb.len() + 5000, nvar, ell);
    let mut lp = inst.fb.len();
    let mut rels: Vec<(Vec<u64>, u64)> = Vec::new();
    for i in 0..b_idx.len() {
        for j in (i + 1)..b_idx.len() {
            if let Some(idxs) = inst.line_relation(b_idx[i], b_idx[j], cost) {
                let mut row = vec![0u64; nvar];
                let mut outside = Vec::new();
                for &id in &idxs {
                    if let Some(pos) = b_idx.iter().position(|&x| x == id) {
                        row[pos] = (row[pos] + 1) % ell;
                    } else {
                        outside.push(id);
                    }
                }
                row[nvar - 1] = (row[nvar - 1] + ell - 1) % ell;
                if outside.is_empty() {
                    rels.push((row, 0));
                } else {
                    let lp_id = lp;
                    lp += 1;
                    if let Some(cyc) = uf.unite(b_idx[i], lp_id, &row, &vec![0; nvar]) {
                        rels.push((cyc, 0));
                    }
                    for o in outside {
                        let _ = uf.unite(lp_id, o.min(lp - 1), &vec![0; nvar], &vec![0; nvar]);
                    }
                }
            }
        }
    }
    let mut rows: Vec<Vec<u64>> = rels.iter().map(|(r, _)| r.clone()).collect();
    let rhs: Vec<u64> = rels.iter().map(|(_, h)| *h).collect();
    if rows.len() < nvar {
        return (None, "rank_not_reached".into());
    }
    match solve_logs(inst, rows, rhs, cost) {
        Some(x) => (Some(x), "ok".into()),
        None => (None, "rank_not_reached".into()),
    }
}

fn cubic_monomials() -> Vec<(u8, u8, u8)> {
    let mut v = Vec::new();
    for i in 0..=3 {
        for j in 0..=(3 - i) {
            let k = 3 - i - j;
            v.push((i as u8, j as u8, k as u8));
        }
    }
    v
}

fn cubic_through(
    inst: &Inst,
    pts: &[Pt],
    cost: &mut Cost,
) -> Option<T> {
    let f = &inst.f1;
    let mons = cubic_monomials();
    let mut eqs = Vec::new();
    for p in pts {
        let w = pt_w(*p);
        let mut row = Vec::new();
        for (i, j, k) in &mons {
            let mut t = 1u32;
            for _ in 0..*i {
                t = f.mul(t, w[0]);
            }
            for _ in 0..*j {
                t = f.mul(t, w[1]);
            }
            for _ in 0..*k {
                t = f.mul(t, w[2]);
            }
            row.push(t);
            cost.fmul((*i + *j + *k) as u64);
        }
        eqs.push(row);
    }
    if eqs.len() < mons.len() {
        return None;
    }
    let n = mons.len();
    let mut mat = eqs;
    for col in 0..n {
        let mut piv = None;
        for r in col..mat.len() {
            if mat[r][col] != 0 {
                piv = Some(r);
                break;
            }
        }
        let r = piv?;
        mat.swap(col, r);
        let inv = f.inv(mat[col][col]);
        cost.finv();
        for c in col..n {
            mat[col][c] = f.mul(mat[col][c], inv);
            cost.fmul(1);
        }
        for r in 0..mat.len() {
            if r == col || mat[r][col] == 0 {
                continue;
            }
            let fac = mat[r][col];
            for c in col..n {
                mat[r][c] = f.sub(mat[r][c], f.mul(fac, mat[col][c]));
                cost.fmul(1);
                cost.fadd(1);
            }
        }
    }
    let mut coeff = vec![0u32; n];
    coeff[n - 1] = 1;
    for r in (0..n.saturating_sub(1)).rev() {
        let mut s = 0u32;
        for c in (r + 1)..n {
            s = f.add(s, f.mul(mat[r][c], coeff[c]));
            cost.fmul(1);
            cost.fadd(1);
        }
        coeff[r] = f.sub(0, s);
        cost.fadd(1);
    }
    let mut t = T::new();
    for (i, c) in coeff.iter().enumerate() {
        if *c != 0 {
            t.insert(mons[i], *c);
        }
    }
    Some(t)
}

fn cubic_residual(
    inst: &Inst,
    cubic: &T,
    known: &HashSet<Pt>,
    cost: &mut Cost,
) -> Vec<Pt> {
    let mut out = Vec::new();
    for p in &inst.all_pts {
        if known.contains(p) {
            continue;
        }
        if teval(&inst.f1, cubic, pt_w(*p)) == 0 {
            out.push(*p);
        }
    }
    cost.fmul(inst.all_pts.len() as u64 * 3);
    out
}

fn decompose_target(
    inst: &Inst,
    logs: &[u64],
    dr: &[Pt],
    cost: &mut Cost,
) -> Option<u64> {
    let ell = inst.ell;
    let nfb = inst.fb.len();
    let mut tries = 0usize;
    let fb_len = inst.fb.len();
    for i in 0..fb_len {
        for j in 0..fb_len {
            for k in 0..fb_len {
                if tries >= CUBIC_BUDGET {
                    return None;
                }
                tries += 1;
                let mut pts: Vec<Pt> = dr.to_vec();
                pts.push(inst.fb[i]);
                pts.push(inst.fb[j]);
                pts.push(inst.fb[k]);
                let known: HashSet<Pt> = pts.iter().copied().collect();
                let cubic = cubic_through(inst, &pts, cost)?;
                let res = cubic_residual(inst, &cubic, &known, cost);
                if res.len() != 3 {
                    continue;
                }
                let mut sum = 0u64;
                for p in &res {
                    let id = inst.pt_map.get(p)?;
                    sum = (sum + logs[*id]) % ell;
                    cost.eadd(1);
                }
                let lhs = (3 * logs[nfb] + ell - sum % ell) % ell;
                cost.emul(1);
                cost.eadd(2);
                return Some(lhs);
            }
        }
    }
    None
}

fn recover_k(_logs: &[u64], y0: u64, y1: u64, y2: u64, ell: u64, cost: &mut Cost) -> Option<u64> {
    let den = (y1 + ell - y0) % ell;
    if den == 0 {
        return None;
    }
    cost.einv();
    cost.emul(2);
    cost.eadd(1);
    Some(((y2 + ell - y0) % ell) * inv_mod(den, ell) % ell)
}

fn rho(inst: &Inst, g: Ep, qpt: Ep, k_true: u64, seed: u64, cost: &mut Cost) -> (Option<u64>, u64) {
    let ell = inst.ell;
    let e = E0 { inst: inst.clone() };
    let mut rng = Rng::seed(seed ^ 0x52484F31);
    let mut mults = Vec::new();
    for _ in 0..RHO_R {
        let a = rng.below(ell);
        let b = rng.below(ell - 1) + 1;
        mults.push((a, b, e.add(e.mul(g, a, cost), e.mul(qpt, b, cost), cost)));
    }
    let mut table: HashMap<Ep, (u64, u64)> = HashMap::new();
    let mut steps = 0u64;
    let mut found = None;
    'outer: loop {
        let mut a = rng.below(ell);
        let mut b = rng.below(ell - 1) + 1;
        let mut l = e.add(e.mul(g, a, cost), e.mul(qpt, b, cost), cost);
        loop {
            steps += 1;
            cost.hash += 1;
            if let Some(&(a2, b2)) = table.get(&l) {
                let db = (b + ell - b2) % ell;
                if db != 0 {
                    let da = (a2 + ell - a) % ell;
                    found = Some((da * inv_mod(db, ell) % ell) as u64);
                    break 'outer;
                }
                break;
            }
            table.insert(l, (a, b));
            let h = l.t as u64 ^ ((l.y as u64).wrapping_mul(0x9E37_79B9));
            let j = (h % RHO_R as u64) as usize;
            l = e.add(l, mults[j].2, cost);
            a = (a + mults[j].0) % ell;
            b = (b + mults[j].1) % ell;
            if steps > 128 * (ell as f64).sqrt() as u64 + 5000 {
                break 'outer;
            }
        }
    }
    (found, steps)
}

// --- divisor class arithmetic (residual intersection) ---
#[derive(Clone)]
struct Div {
    terms: Vec<(Pt, i32)>,
}

impl Div {
    fn zero() -> Div {
        Div { terms: vec![] }
    }
    fn point(p: Pt) -> Div {
        Div {
            terms: vec![(p, 1)],
        }
    }
    fn add(&self, o: &Div) -> Div {
        let mut m: HashMap<Pt, i32> = HashMap::new();
        for (p, c) in self.terms.iter().chain(o.terms.iter()) {
            *m.entry(*p).or_insert(0) += c;
        }
        let terms = m.into_iter().filter(|(_, c)| *c != 0).collect();
        Div { terms }
    }
    fn sub(&self, o: &Div) -> Div {
        let mut m: HashMap<Pt, i32> = HashMap::new();
        for (p, c) in &self.terms {
            *m.entry(*p).or_insert(0) += c;
        }
        for (p, c) in &o.terms {
            *m.entry(*p).or_insert(0) -= c;
        }
        let terms = m.into_iter().filter(|(_, c)| *c != 0).collect();
        Div { terms }
    }
    fn deg(&self) -> i32 {
        self.terms.iter().map(|(_, c)| c).sum()
    }
}

impl Inst {
    fn reduce_div(&self, mut d: Div, cost: &mut Cost) -> Div {
        loop {
            let deg = d.deg();
            if deg <= 3 {
                break;
            }
            let pts: Vec<Pt> = d.terms.iter().map(|(p, _)| *p).collect();
            if pts.len() < 2 {
                break;
            }
            let a = pts[0];
            let b = pts[1];
            let line = self.line_through(a, b, cost);
            let extra = self.restrict_line(line, [a, b], cost);
            let mut nd = d.sub(&Div::point(a).add(&Div::point(b)));
            for p in extra {
                nd = nd.add(&Div::point(p));
            }
            d = nd;
        }
        d
    }
    fn canon(&self, d: &Div, cost: &mut Cost) -> Vec<Pt> {
        let r = self.reduce_div(d.clone(), cost);
        let mut v: Vec<Pt> = r.terms.iter().map(|(p, _)| *p).collect();
        v.sort_by_key(|p| (p.0[0], p.0[1], p.0[2]));
        v
    }
    fn equiv(&self, a: &Div, b: &Div, cost: &mut Cost) -> bool {
        self.canon(a, cost) == self.canon(b, cost)
    }
}

fn check_lin_equiv(q: u32) -> (bool, usize) {
    let mut cost = Cost::new(q);
    let inst = match Inst::build_first_split(q, &mut cost) {
        Some(x) => x,
        None => return (false, 0),
    };
    let e = E0 { inst: inst.clone() };
    let g = match e.find_generator(&mut cost) {
        Some(x) => x,
        None => return (false, 0),
    };
    let n_ws = e.transfer_points(g, &mut cost).len();
    let ell = inst.ell;
    let mut rng = Rng::seed(PROTO_SEED ^ q as u64);
    let mut ok = true;
    let two_g = e.add(g, g, &mut cost);
    if e.pushforward_phi(g, &mut cost) != Some(two_g) || two_g.inf {
        ok = false;
    }
    if !e.mul(g, ell, &mut cost).inf {
        ok = false;
    }
    for _ in 0..6 {
        let a = rng.below(ell - 1) + 1;
        let b = rng.below(ell - 1) + 1;
        let p = e.mul(g, a, &mut cost);
        let qpt = e.mul(g, b, &mut cost);
        if !e.valid_transfer(p.t) || !e.valid_transfer(qpt.t) {
            continue;
        }
        let two_p = e.add(p, p, &mut cost);
        let two_q = e.add(qpt, qpt, &mut cost);
        let pf_p = e.pushforward_phi(p, &mut cost);
        let pf_q = e.pushforward_phi(qpt, &mut cost);
        let pq = e.add(p, qpt, &mut cost);
        let pf_pq = e.pushforward_phi(pq, &mut cost);
        if pf_p != Some(two_p) || pf_q != Some(two_q) {
            ok = false;
        }
        if pf_pq != Some(e.add(two_p, two_q, &mut cost)) {
            ok = false;
        }
    }
    (ok, n_ws)
}

fn run_trial(q: u32, trial: u32, variant: &str) {
    let mut cost = Cost::new(q);
    let inst = match Inst::build_first_split(q, &mut cost) {
        Some(x) => x,
        None => {
            println!(
                "{{\"rec\":\"trial\",\"q\":{},\"trial\":{},\"variant\":\"{}\",\"status\":\"no_instance\"}}",
                q, trial, variant
            );
            return;
        }
    };
    let e = E0 { inst: inst.clone() };
    let g = match e.find_generator(&mut cost) {
        Some(x) => x,
        None => {
            println!(
                "{{\"rec\":\"trial\",\"q\":{},\"trial\":{},\"variant\":\"{}\",\"status\":\"no_generator\"}}",
                q, trial, variant
            );
            return;
        }
    };
    let mut rng = Rng::seed(PROTO_SEED ^ ((q as u64) << 32) ^ trial as u64);
    let k = rng.below(inst.ell - 1) + 1;
    let qpt = e.mul(g, k, &mut cost);
    let ell = inst.ell;
    let sqrt_ell = (ell as f64).sqrt();

    let mut run_ic = |v: &str| {
        let mut c = Cost::new(q);
        let inst = Inst::build_first_split(q, &mut c).unwrap();
        let e = E0 { inst: inst.clone() };
        let (logs, status) = if v == "full" {
            ic_full(&inst, &mut c)
        } else {
            ic_2lp(&inst, &mut c)
        };
        let r_probe = rng.below(ell - 1) + 1;
        let mut answer = None;
        let mut ok = false;
        let mut recovery = "";
        let mut ic_answer = false;
        if let Some(ref lv) = logs {
            if lv.len() >= inst.fb.len() + 1 {
                let r0 = e.mul(g, r_probe, &mut c);
                let r1 = e.add(r0, g, &mut c);
                let r2 = e.add(r0, qpt, &mut c);
                let mut ys = Vec::new();
                for rr in [r0, r1, r2] {
                    let dr = e.transfer_points(rr, &mut c);
                    if dr.len() != 6 {
                        continue;
                    }
                    if let Some(y) = decompose_target(&inst, lv, &dr, &mut c) {
                        ys.push(y);
                    }
                }
                if ys.len() == 3 {
                    if let Some(kk) = recover_k(lv, ys[0], ys[1], ys[2], ell, &mut c) {
                        answer = Some(kk);
                        ok = kk == k;
                        ic_answer = true;
                        recovery = "ic";
                    }
                }
            }
        }
        let ic_ops = c.total();
        let mut transfer_ops = 0u64;
        if answer.is_none() {
            let t0 = c.total();
            if let Some(kk) = transfer_pushforward_dlog(&e, g, qpt, r_probe, ell, &mut c) {
                transfer_ops = c.total() - t0;
                answer = Some(kk);
                ok = kk == k;
                recovery = "transfer_dlog";
            }
        }
        let st = if answer.is_some() {
            if !ok {
                "wrong_answer"
            } else if ic_answer {
                "ok"
            } else if v == "full" {
                "ic_fail_transfer_dlog_ok"
            } else {
                "ok_transfer_dlog"
            }
        } else if status == "ok" {
            "no_decomposition"
        } else {
            &status
        };
        let s = c.total() as f64 / sqrt_ell;
        if recovery == "transfer_dlog" && ic_ops > 0 {
            println!(
                "{{\"rec\":\"trial\",\"q\":{},\"trial\":{},\"variant\":\"{}\",\"status\":\"{}\",\"recovery\":\"{}\",\"k_planted\":{},\"k_found\":{},\"correct\":{},\"ops\":{},\"ops_ic\":{},\"ops_transfer\":{},\"S\":{:.6},\"hash_lookups\":{}}}",
                q, trial, v, st, recovery, k,
                answer.map(|x| x.to_string()).unwrap_or_else(|| "null".into()),
                ok, c.total(), ic_ops, transfer_ops, s, c.hash
            );
        } else {
            println!(
                "{{\"rec\":\"trial\",\"q\":{},\"trial\":{},\"variant\":\"{}\",\"status\":\"{}\",\"recovery\":\"{}\",\"k_planted\":{},\"k_found\":{},\"correct\":{},\"ops\":{},\"S\":{:.6},\"hash_lookups\":{}}}",
                q, trial, v, st, recovery, k,
                answer.map(|x| x.to_string()).unwrap_or_else(|| "null".into()),
                ok, c.total(), s, c.hash
            );
        }
        (s, ok)
    };

    if variant == "rho" || variant == "all" {
        let mut c = Cost::new(q);
        let inst = Inst::build_first_split(q, &mut c).unwrap();
        let e = E0 { inst: inst.clone() };
        let g = e.find_generator(&mut c).unwrap();
        let (found, steps) = rho(&inst, g, qpt, k, PROTO_SEED ^ trial as u64, &mut c);
        let ok = found == Some(k);
        let s = c.total() as f64 / sqrt_ell;
        println!(
            "{{\"rec\":\"trial\",\"q\":{},\"trial\":{},\"variant\":\"rho\",\"status\":\"{}\",\"k_planted\":{},\"k_found\":{},\"correct\":{},\"ops\":{},\"S\":{:.6},\"rho_steps\":{}}}",
            q,
            trial,
            if ok { "ok" } else { "fail" },
            k,
            found.map(|x| x.to_string()).unwrap_or_else(|| "null".into()),
            ok,
            c.total(),
            s,
            steps
        );
    }
    if variant == "full" || variant == "all" {
        run_ic("full");
    }
    if variant == "2lp" || variant == "all" {
        run_ic("2lp");
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|a| a == "--check-lin-equiv") {
        let q: u32 = args
            .iter()
            .find_map(|a| a.parse().ok())
            .unwrap_or(7);
        let (ok, n_ws) = check_lin_equiv(q);
        println!(
            "{{\"rec\":\"check-lin-equiv\",\"q\":{},\"check\":\"pushforward_identity\",\"pass\":{},\"n_transfer\":{}}}",
            q, ok, n_ws
        );
        return;
    }
    if args.len() < 4 {
        eprintln!("usage: quartic_ic <q> <trials> <full|2lp|rho|all>");
        eprintln!("       quartic_ic <q> <trials> --check-lin-equiv");
        std::process::exit(1);
    }
    let q: u32 = args[1].parse().expect("q");
    let trials: u32 = args[2].parse().expect("trials");
    let variant = if args[3] == "--check-lin-equiv" {
        let (ok, n_ws) = check_lin_equiv(q);
        println!(
            "{{\"rec\":\"check-lin-equiv\",\"q\":{},\"check\":\"pushforward_identity\",\"pass\":{},\"n_transfer\":{}}}",
            q, ok, n_ws
        );
        return;
    } else {
        args[3].as_str()
    };
    for t in 0..trials {
        run_trial(q, t, variant);
    }
    println!(
        "{{\"rec\":\"summary\",\"q\":{},\"trials\":{},\"variant\":\"{}\"}}",
        q, trials, variant
    );
}
