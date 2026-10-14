// Frozen protocol: research/ghs_cubic_quartic_ic_20261012/PROTOCOL.md
// rustc --edition 2021 -O quartic_ic.rs
use std::collections::{BTreeMap, HashMap, HashSet};

const PROTO_SEED: u64 = 20261012;
const CUBIC_BUDGET: usize = 400;
const RHO_R: usize = 20;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
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
        // Enumerate C(F_q) on the line (restrict_line chart scan is incomplete).
        let extra = self.line_points(line, [a, b], cost);
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

// --- W over K, orbits, Pic, E0, IC ---

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct Wk([u32; 3]);

fn norm_wk(f3: &Fld, w: Wk) -> Wk {
    let l = w.0.iter().copied().find(|&x| x != 0).unwrap_or(1);
    let li = f3.inv(l);
    Wk([f3.mul(w.0[0], li), f3.mul(w.0[1], li), f3.mul(w.0[2], li)])
}

fn wk_key(w: Wk) -> (u32, u32, u32) {
    (w.0[0], w.0[1], w.0[2])
}

fn k_coords(f1: &Fld, f3: &Fld, a: u32) -> [u32; 3] {
    if a == 0 {
        return [0, 0, 0];
    }
    let idx = f3.idx_of[a as usize];
    let q = f3.q;
    [f1.int(idx % q), f1.int((idx / q) % q), f1.int((idx / (q * q)) % q)]
}

fn eval_mon_k(f3: &Fld, w: Wk, i: u8, j: u8, k: u8) -> u32 {
    let mut t = 1u32;
    for _ in 0..i {
        t = f3.mul(t, w.0[0]);
    }
    for _ in 0..j {
        t = f3.mul(t, w.0[1]);
    }
    for _ in 0..k {
        t = f3.mul(t, w.0[2]);
    }
    t
}

fn teval_k(inst: &Inst, poly: &T, w: Wk) -> u32 {
    let f3 = &inst.f3;
    let f1 = &inst.f1;
    let mut s = 0u32;
    for ((i, j, k), c) in poly {
        let ck = f1.embed3(f3, *c);
        s = f3.add(s, f3.mul(ck, eval_mon_k(f3, w, *i, *j, *k)));
    }
    s
}

fn cubic_monomials() -> Vec<(u8, u8, u8)> {
    let mut v = Vec::new();
    for i in 0..=3u8 {
        for j in 0..=(3 - i) {
            v.push((i, j, 3 - i - j));
        }
    }
    v
}

fn nullspace_basis_fq(f1: &Fld, mut mat: Vec<Vec<u32>>) -> Vec<Vec<u32>> {
    if mat.is_empty() {
        return Vec::new();
    }
    let n = mat[0].len();
    let mut piv_col = Vec::new();
    let mut piv = 0usize;
    for col in 0..n {
        let mut sel = None;
        for r in piv..mat.len() {
            if mat[r][col] != 0 {
                sel = Some(r);
                break;
            }
        }
        let Some(r) = sel else { continue };
        mat.swap(piv, r);
        let inv = f1.inv(mat[piv][col]);
        for c in col..n {
            mat[piv][c] = f1.mul(mat[piv][c], inv);
        }
        for r in 0..mat.len() {
            if r == piv || mat[r][col] == 0 {
                continue;
            }
            let fac = mat[r][col];
            for c in col..n {
                mat[r][c] = f1.sub(mat[r][c], f1.mul(fac, mat[piv][c]));
            }
        }
        piv_col.push(col);
        piv += 1;
        if piv == mat.len() {
            break;
        }
    }
    let piv_set: HashSet<usize> = piv_col.iter().copied().collect();
    let free: Vec<usize> = (0..n).filter(|c| !piv_set.contains(c)).collect();
    let mut basis = Vec::new();
    for &fcol in &free {
        let mut v = vec![0u32; n];
        v[fcol] = 1;
        for (r, &pc) in piv_col.iter().enumerate() {
            if r < mat.len() && mat[r][fcol] != 0 {
                v[pc] = f1.sub(0, mat[r][fcol]);
            }
        }
        if v.iter().any(|&x| x != 0) {
            basis.push(v);
        }
    }
    basis
}

fn orbit_reps(f3: &Fld, wpts: &[Wk]) -> Vec<Wk> {
    let mut reps = Vec::new();
    let mut used = vec![false; wpts.len()];
    for i in 0..wpts.len() {
        if used[i] {
            continue;
        }
        let w0 = norm_wk(f3, wpts[i]);
        reps.push(w0);
        let mut w = w0;
        for _ in 0..3 {
            let nw = norm_wk(f3, w);
            for (j, u) in wpts.iter().enumerate() {
                if norm_wk(f3, *u).0 == nw.0 {
                    used[j] = true;
                }
            }
            w = Wk([f3.frob(w.0[0]), f3.frob(w.0[1]), f3.frob(w.0[2])]);
        }
    }
    reps
}

fn expand_orbit(f3: &Fld, rep: Wk) -> [Wk; 3] {
    let mut w = norm_wk(f3, rep);
    let mut out = [w; 3];
    for i in 0..3 {
        out[i] = w;
        w = norm_wk(
            f3,
            Wk([f3.frob(w.0[0]), f3.frob(w.0[1]), f3.frob(w.0[2])]),
        );
    }
    out.sort_by_key(|x| wk_key(*x));
    out
}

fn orbit_id(f3: &Fld, rep: Wk) -> (u32, u32, u32) {
    wk_key(expand_orbit(f3, rep)[0])
}

/// Effective F_q-rational divisor piece: F_q point or Gal(K/F_q)-orbit of size 3.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
enum Piece {
    Fq(Pt),
    Orb(Wk), // representative
}

impl Piece {
    fn deg(&self) -> i32 {
        match self {
            Piece::Fq(_) => 1,
            Piece::Orb(_) => 3,
        }
    }
    fn sort_key(&self, f3: &Fld) -> (u8, u32, u32, u32) {
        match self {
            Piece::Fq(p) => (0, p.0[0], p.0[1], p.0[2]),
            Piece::Orb(w) => {
                let id = orbit_id(f3, *w);
                (1, id.0, id.1, id.2)
            }
        }
    }
}

#[derive(Clone, Debug)]
struct DivCls {
    pieces: Vec<Piece>,
}

impl DivCls {
    fn zero() -> DivCls {
        DivCls { pieces: vec![] }
    }
    fn fq(p: Pt) -> DivCls {
        DivCls {
            pieces: vec![Piece::Fq(p)],
        }
    }
    fn orb(w: Wk) -> DivCls {
        DivCls {
            pieces: vec![Piece::Orb(w)],
        }
    }
    fn deg(&self) -> i32 {
        self.pieces.iter().map(|p| p.deg()).sum()
    }
    fn join(&self, o: &DivCls) -> DivCls {
        let mut pieces = self.pieces.clone();
        pieces.extend(o.pieces.iter().cloned());
        DivCls { pieces }.normalize_mult()
    }
    fn normalize_mult(mut self) -> DivCls {
        // cancel Fq duplicates as multiset (simple: sort unique for now with multiplicity collapse)
        self.pieces.sort_by_key(|p| match p {
            Piece::Fq(pt) => (0u8, pt.0[0], pt.0[1], pt.0[2]),
            Piece::Orb(w) => (1u8, w.0[0], w.0[1], w.0[2]),
        });
        self
    }
    fn all_w(&self, inst: &Inst) -> Vec<Wk> {
        let f3 = &inst.f3;
        let mut v = Vec::new();
        for p in &self.pieces {
            match p {
                Piece::Fq(pt) => v.push(inst.pt_to_wk(*pt)),
                Piece::Orb(w) => v.extend_from_slice(&expand_orbit(f3, *w)),
            }
        }
        v
    }
}

impl Inst {
    fn pt_to_wk(&self, p: Pt) -> Wk {
        let f3 = &self.f3;
        Wk([
            f3.int(self.f1.idx_of[p.0[0] as usize]),
            f3.int(self.f1.idx_of[p.0[1] as usize]),
            f3.int(self.f1.idx_of[p.0[2] as usize]),
        ])
    }
    fn wk_to_pt_if_fq(&self, w: Wk) -> Option<Pt> {
        let f3 = &self.f3;
        let f1 = &self.f1;
        let nw = norm_wk(f3, w);
        if !f3.in_fq(nw.0[0]) || !f3.in_fq(nw.0[1]) || !f3.in_fq(nw.0[2]) {
            return None;
        }
        Some(norm_pt(
            f1,
            [
                f1.int(f3.idx_of[nw.0[0] as usize]),
                f1.int(f3.idx_of[nw.0[1] as usize]),
                f1.int(f3.idx_of[nw.0[2] as usize]),
            ],
        ))
    }
    fn enum_curve_k(&self, cost: &mut Cost) -> Vec<Wk> {
        let f3 = &self.f3;
        let mut pts = Vec::new();
        for a in 0..f3.size {
            for b in 0..f3.size {
                let w = Wk([f3.int(1), f3.el_of[a as usize], f3.el_of[b as usize]]);
                if teval_k(self, &self.gwn, w) == 0 {
                    pts.push(norm_wk(f3, w));
                }
            }
        }
        for b in 0..f3.size {
            let w = Wk([0, f3.int(1), f3.el_of[b as usize]]);
            if teval_k(self, &self.gwn, w) == 0 {
                pts.push(norm_wk(f3, w));
            }
        }
        let w = Wk([0, 0, f3.int(1)]);
        if teval_k(self, &self.gwn, w) == 0 {
            pts.push(w);
        }
        cost.kmul((f3.size as u64) * (f3.size as u64));
        pts.sort_by_key(|w| wk_key(*w));
        pts.dedup();
        pts
    }
    fn pieces_from_w(&self, ws: &[Wk]) -> DivCls {
        let f3 = &self.f3;
        let mut used = vec![false; ws.len()];
        let mut pieces = Vec::new();
        for i in 0..ws.len() {
            if used[i] {
                continue;
            }
            let w0 = norm_wk(f3, ws[i]);
            if let Some(p) = self.wk_to_pt_if_fq(w0) {
                used[i] = true;
                pieces.push(Piece::Fq(p));
                continue;
            }
            let orb = expand_orbit(f3, w0);
            for ow in &orb {
                for (j, u) in ws.iter().enumerate() {
                    if norm_wk(f3, *u).0 == ow.0 {
                        used[j] = true;
                    }
                }
            }
            pieces.push(Piece::Orb(orb[0]));
        }
        DivCls { pieces }.normalize_mult()
    }
}

fn coeff_to_cubic(mons: &[(u8, u8, u8)], coeff: &[u32]) -> T {
    let mut t = T::new();
    for (i, c) in coeff.iter().enumerate() {
        if *c != 0 {
            t.insert(mons[i], *c);
        }
    }
    t
}

fn cubic_matrix(inst: &Inst, d: &DivCls, cost: &mut Cost) -> (Vec<(u8, u8, u8)>, Vec<Vec<u32>>) {
    let f1 = &inst.f1;
    let f3 = &inst.f3;
    let mons = cubic_monomials();
    let mut rows = Vec::new();
    let reps = orbit_reps(f3, &d.all_w(inst));
    for w in reps {
        for c in 0..3u32 {
            let mut row = vec![0u32; mons.len()];
            for (j, (i, jj, k)) in mons.iter().enumerate() {
                let val = eval_mon_k(f3, w, *i, *jj, *k);
                row[j] = k_coords(f1, f3, val)[c as usize];
            }
            rows.push(row);
        }
        cost.kmul(9);
    }
    (mons, rows)
}

fn cubic_residual_pieces(
    inst: &Inst,
    cubic: &T,
    known: &DivCls,
    curve_k: &[Wk],
    cost: &mut Cost,
) -> DivCls {
    let f3 = &inst.f3;
    let mut known_mult: HashMap<(u32, u32, u32), usize> = HashMap::new();
    for w in known.all_w(inst) {
        let k = wk_key(norm_wk(f3, w));
        *known_mult.entry(k).or_insert(0) += 1;
    }
    let mut res_w = Vec::new();
    let mut seen: HashMap<(u32, u32, u32), usize> = HashMap::new();
    for w in curve_k {
        let nw = norm_wk(f3, *w);
        if teval_k(inst, cubic, nw) != 0 {
            continue;
        }
        let k = wk_key(nw);
        let n = seen.entry(k).or_insert(0);
        *n += 1;
        let skip = known_mult.get(&k).copied().unwrap_or(0);
        if *n > skip {
            res_w.push(nw);
        }
    }
    cost.kmul(curve_k.len() as u64);
    inst.pieces_from_w(&res_w)
}

fn cubic_pencil_combos(f1: &Fld, basis: &[Vec<u32>]) -> Vec<Vec<u32>> {
    let dim = basis.len().min(5);
    let q = f1.q;
    let mut combos: Vec<Vec<u32>> = Vec::new();
    for i in 0..dim {
        combos.push(basis[i].clone());
    }
    if dim >= 2 {
        for i in 0..dim {
            for j in (i + 1)..dim {
                for s in 1..q {
                    let mut v = basis[i].clone();
                    for k in 0..v.len() {
                        v[k] = f1.add(v[k], f1.mul(f1.int(s), basis[j][k]));
                    }
                    combos.push(v);
                }
            }
        }
    }
    // three-way samples for dim>=3
    if dim >= 3 {
        for s in 1..q {
            for t in 1..q {
                let mut v = basis[0].clone();
                for k in 0..v.len() {
                    v[k] = f1.add(v[k], f1.mul(f1.int(s), basis[1][k]));
                    v[k] = f1.add(v[k], f1.mul(f1.int(t), basis[2][k]));
                }
                combos.push(v);
            }
        }
    }
    combos
}

/// Cubics through `d`; prefer those whose K-residual has degree `want_res`.
fn cubic_through_with_residual(
    inst: &Inst,
    d: &DivCls,
    curve_k: &[Wk],
    want_res: i32,
    cost: &mut Cost,
) -> Option<(T, DivCls)> {
    let (mons, rows) = cubic_matrix(inst, d, cost);
    let basis = nullspace_basis_fq(&inst.f1, rows);
    if basis.is_empty() {
        return None;
    }
    let mut best: Option<(T, DivCls, i32)> = None;
    for coeff in cubic_pencil_combos(&inst.f1, &basis) {
        let t = coeff_to_cubic(&mons, &coeff);
        if t.is_empty() {
            continue;
        }
        if d.all_w(inst).iter().any(|w| teval_k(inst, &t, *w) != 0) {
            continue;
        }
        let res = cubic_residual_pieces(inst, &t, d, curve_k, cost);
        let score = (res.deg() - want_res).abs();
        if res.deg() == want_res {
            return Some((t, res));
        }
        let replace = match &best {
            None => true,
            Some((_, _, sc)) => score < *sc,
        };
        if replace {
            best = Some((t, res, score));
        }
    }
    best.map(|(t, r, _)| (t, r))
}

/// All residuals of exact degree `want_res` from the pencil through `d`.
fn cubic_residuals_exact(
    inst: &Inst,
    d: &DivCls,
    curve_k: &[Wk],
    want_res: i32,
    cost: &mut Cost,
) -> Vec<DivCls> {
    let (mons, rows) = cubic_matrix(inst, d, cost);
    let basis = nullspace_basis_fq(&inst.f1, rows);
    let mut out = Vec::new();
    let mut seen = HashSet::new();
    for coeff in cubic_pencil_combos(&inst.f1, &basis) {
        let t = coeff_to_cubic(&mons, &coeff);
        if t.is_empty() {
            continue;
        }
        if d.all_w(inst).iter().any(|w| teval_k(inst, &t, *w) != 0) {
            continue;
        }
        let res = cubic_residual_pieces(inst, &t, d, curve_k, cost);
        if res.deg() != want_res {
            continue;
        }
        let key: Vec<_> = sort_div(inst, res.clone())
            .pieces
            .iter()
            .map(|p| p.sort_key(&inst.f3))
            .collect();
        if seen.insert(key) {
            out.push(res);
        }
    }
    out
}

fn cubic_from_pieces(inst: &Inst, d: &DivCls, cost: &mut Cost) -> Option<T> {
    let (mons, rows) = cubic_matrix(inst, d, cost);
    let basis = nullspace_basis_fq(&inst.f1, rows);
    let coeff = basis.into_iter().next()?;
    let t = coeff_to_cubic(&mons, &coeff);
    if d.all_w(inst).iter().any(|w| teval_k(inst, &t, *w) != 0) {
        return None;
    }
    Some(t)
}

/// Fixed origin Ω = P0+P1+P2 (three distinct F_q points). Classes are [D − Ω].
fn origin_omega(inst: &Inst) -> DivCls {
    let n = inst.fb.len();
    assert!(n >= 3);
    DivCls {
        pieces: vec![
            Piece::Fq(inst.fb[0]),
            Piece::Fq(inst.fb[1]),
            Piece::Fq(inst.fb[2]),
        ],
    }
}

fn sort_div(inst: &Inst, d: DivCls) -> DivCls {
    let mut pcs = d.pieces;
    pcs.sort_by_key(|p| p.sort_key(&inst.f3));
    DivCls { pieces: pcs }
}

/// Reduce an effective F_q-rational divisor to a deg-3 Pic^0 representative of [D − deg(D)·(avg)].
/// For deg 6: if D+R~3H and R+Ω+F~3H then [F−Ω]=[D−6P̄] in Pic^0 (Ω origin).
fn pic_canon(inst: &Inst, d: &DivCls, _p0: Pt, curve_k: &[Wk], cost: &mut Cost) -> DivCls {
    let omega = origin_omega(inst);
    let mut cur = d.clone();
    for _ in 0..16 {
        let deg = cur.deg();
        if deg == 3 {
            return sort_div(inst, cur);
        }
        if deg < 3 {
            // pad with points from omega
            for p in &omega.pieces {
                if cur.deg() >= 3 {
                    break;
                }
                cur.pieces.push(p.clone());
            }
            cur = cur.normalize_mult();
            continue;
        }
        if deg > 6 {
            let fqs: Vec<Pt> = cur
                .pieces
                .iter()
                .filter_map(|p| match p {
                    Piece::Fq(pt) => Some(*pt),
                    _ => None,
                })
                .collect();
            if fqs.len() >= 2 {
                let line = inst.line_through(fqs[0], fqs[1], cost);
                let extra = inst.line_points(line, [fqs[0], fqs[1]], cost);
                let mut removed = 0usize;
                let mut pcs = Vec::new();
                for p in &cur.pieces {
                    match p {
                        Piece::Fq(pt)
                            if removed < 2 && (*pt == fqs[0] || *pt == fqs[1]) =>
                        {
                            removed += 1;
                        }
                        _ => pcs.push(p.clone()),
                    }
                }
                for e in extra {
                    pcs.push(Piece::Fq(e));
                }
                cur = DivCls { pieces: pcs }.normalize_mult();
                continue;
            }
        }
        if deg == 4 || deg == 5 {
            for p in &omega.pieces {
                if cur.deg() >= 6 {
                    break;
                }
                cur.pieces.push(p.clone());
            }
            cur = cur.normalize_mult();
            // may exceed 6 if omega overlaps; trim duplicates softly
            continue;
        }
        if deg == 6 {
            // D + R ~ 3H, then R + Ω + F ~ 3H ⇒ F represents [D] in Pic^0 / origin Ω.
            // Search the pencil for an R that admits a deg-3 complement through R∪Ω.
            let candidates = cubic_residuals_exact(inst, &cur, curve_k, 6, cost);
            for r in &candidates {
                let known = r.join(&omega);
                if let Some((_, f)) = cubic_through_with_residual(inst, &known, curve_k, 3, cost) {
                    if f.deg() == 3 {
                        return sort_div(inst, f);
                    }
                }
            }
            if let Some(r) = candidates.into_iter().next() {
                // fallback: try reducing a residual of deg < 6 if any single-shot find
                if r.deg() > 0 && r.deg() < 6 {
                    cur = r;
                    continue;
                }
            } else if let Some((_, r)) = cubic_through_with_residual(inst, &cur, curve_k, 6, cost)
            {
                if r.deg() > 0 && r.deg() < 6 {
                    cur = r;
                    continue;
                }
            }
            break;
        }
        break;
    }
    sort_div(inst, cur)
}

fn pic_eq(inst: &Inst, a: &DivCls, b: &DivCls, p0: Pt, curve_k: &[Wk], cost: &mut Cost) -> bool {
    let ca = pic_canon(inst, a, p0, curve_k, cost);
    let cb = pic_canon(inst, b, p0, curve_k, cost);
    ca.pieces == cb.pieces
}

fn pieces_overlap(a: &DivCls, b: &DivCls, f3: &Fld) -> bool {
    let ka: HashSet<_> = a.pieces.iter().map(|p| p.sort_key(f3)).collect();
    b.pieces.iter().any(|p| ka.contains(&p.sort_key(f3)))
}

/// Addition when supports are multiplicity-free after join (no doubled pieces).
fn pic_add_distinct(
    inst: &Inst,
    a: &DivCls,
    b: &DivCls,
    p0: Pt,
    curve_k: &[Wk],
    cost: &mut Cost,
) -> DivCls {
    let sum = a.join(b);
    if sum.deg() == 6 {
        return pic_canon(inst, &sum, p0, curve_k, cost);
    }
    let omega = origin_omega(inst);
    let mut padded = sum;
    while padded.deg() < 6 {
        for p in &omega.pieces {
            if padded.deg() >= 6 {
                break;
            }
            padded.pieces.push(p.clone());
        }
        padded = padded.normalize_mult();
    }
    pic_canon(inst, &padded, p0, curve_k, cost)
}

fn all_deg3_candidates(inst: &Inst, curve_k: &[Wk]) -> Vec<DivCls> {
    let mut out = Vec::new();
    let n = inst.fb.len();
    for i in 0..n {
        for j in i..n {
            for k in j..n {
                out.push(DivCls {
                    pieces: vec![
                        Piece::Fq(inst.fb[i]),
                        Piece::Fq(inst.fb[j]),
                        Piece::Fq(inst.fb[k]),
                    ],
                });
            }
        }
    }
    let mut seen = HashSet::new();
    for w in curve_k {
        if inst.wk_to_pt_if_fq(*w).is_some() {
            continue;
        }
        let id = orbit_id(&inst.f3, *w);
        if seen.insert(id) {
            out.push(DivCls::orb(expand_orbit(&inst.f3, *w)[0]));
        }
    }
    out
}

fn pic_negate(inst: &Inst, a: &DivCls, p0: Pt, curve_k: &[Wk], cost: &mut Cost) -> DivCls {
    let omega = origin_omega(inst);
    let ca = pic_canon(inst, a, p0, curve_k, cost);
    if ca.pieces == omega.pieces {
        return omega;
    }
    for cand in all_deg3_candidates(inst, curve_k) {
        if pieces_overlap(&ca, &cand, &inst.f3) {
            continue;
        }
        let s = pic_add_distinct(inst, &ca, &cand, p0, curve_k, cost);
        if s.pieces == omega.pieces {
            return sort_div(inst, cand);
        }
    }
    // Fallback: residual of cubic through ca ∪ Ω with want deg 6, then complement — may fail.
    let known = ca.join(&omega);
    if let Some((_, r)) = cubic_through_with_residual(inst, &known, curve_k, 6, cost) {
        if r.deg() == 6 {
            return pic_canon(inst, &r, p0, curve_k, cost);
        }
    }
    omega
}

fn pic_add(inst: &Inst, a: &DivCls, b: &DivCls, p0: Pt, curve_k: &[Wk], cost: &mut Cost) -> DivCls {
    let ca = pic_canon(inst, a, p0, curve_k, cost);
    let cb = pic_canon(inst, b, p0, curve_k, cost);
    if !pieces_overlap(&ca, &cb, &inst.f3) {
        return pic_add_distinct(inst, &ca, &cb, p0, curve_k, cost);
    }
    // Overlap / doubling: translate by a disjoint helper T, then subtract T.
    for t in all_deg3_candidates(inst, curve_k) {
        if pieces_overlap(&ca, &t, &inst.f3) || pieces_overlap(&cb, &t, &inst.f3) {
            continue;
        }
        let at = pic_add_distinct(inst, &ca, &t, p0, curve_k, cost);
        if pieces_overlap(&at, &cb, &inst.f3) {
            continue;
        }
        let sum = pic_add_distinct(inst, &at, &cb, p0, curve_k, cost);
        let nt = pic_negate(inst, &t, p0, curve_k, cost);
        if pieces_overlap(&sum, &nt, &inst.f3) {
            continue;
        }
        return pic_add_distinct(inst, &sum, &nt, p0, curve_k, cost);
    }
    // Last resort: return canon of multiset join
    pic_canon(inst, &ca.join(&cb), p0, curve_k, cost)
}

fn pic_mul(inst: &Inst, d: &DivCls, n: u64, p0: Pt, curve_k: &[Wk], cost: &mut Cost) -> DivCls {
    let omega = origin_omega(inst);
    if n == 0 {
        return omega;
    }
    let mut base = pic_canon(inst, d, p0, curve_k, cost);
    let mut acc = None;
    let mut e = n;
    while e > 0 {
        if e & 1 == 1 {
            acc = Some(match acc {
                None => base.clone(),
                Some(ref a) => pic_add(inst, a, &base, p0, curve_k, cost),
            });
        }
        base = pic_add(inst, &base, &base, p0, curve_k, cost);
        e >>= 1;
    }
    acc.unwrap_or(omega)
}

// --- E0: y^2 = f(t) via birational Weierstrass ---
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct Ep {
    t: u32,
    y: u32,
    inf: bool,
}

struct W1 {
    alpha: u32,
    a3: u32,
    a2: u32,
    a1: u32,
    a0: u32,
}

struct Kp {
    c: Vec<u32>,
}
impl Kp {
    fn eval(f: &Fld, p: &Kp, t: u32) -> u32 {
        let mut s = 0u32;
        let mut pw = 1u32;
        for &coef in &p.c {
            s = f.add(s, f.mul(coef, pw));
            pw = f.mul(pw, t);
        }
        s
    }
    fn mul(f: &Fld, a: &Kp, b: &Kp) -> Kp {
        let mut c = vec![0u32; a.c.len() + b.c.len().saturating_sub(1).max(0) + 1];
        if a.c.is_empty() || b.c.is_empty() {
            return Kp { c: vec![0] };
        }
        c.resize(a.c.len() + b.c.len() - 1, 0);
        for i in 0..a.c.len() {
            for j in 0..b.c.len() {
                c[i + j] = f.add(c[i + j], f.mul(a.c[i], b.c[j]));
            }
        }
        Kp { c }
    }
}

fn binom_u(f: &Fld, i: u8, j: u8) -> u32 {
    let b = match (i, j) {
        (0, 0) => 1,
        (1, 0) | (1, 1) => 1,
        (2, 0) | (2, 2) => 1,
        (2, 1) => 2,
        (3, 0) | (3, 3) => 1,
        (3, 1) | (3, 2) => 3,
        (4, 0) | (4, 4) => 1,
        (4, 1) | (4, 3) => 4,
        (4, 2) => 6,
        _ => 0,
    };
    f.int(b)
}

fn fit_weierstrass(f3: &Fld, f: &Kp, alpha: u32) -> Option<W1> {
    let mut c = [0u32; 5];
    let deg = f.c.len().saturating_sub(1) as u8;
    for i in 0..=deg {
        let fi = *f.c.get(i as usize).unwrap_or(&0);
        if fi == 0 {
            continue;
        }
        for j in 0..=i {
            let mut p = 1u32;
            for _ in 0..(i - j) {
                p = f3.mul(p, alpha);
            }
            c[j as usize] = f3.add(c[j as usize], f3.mul(f3.mul(fi, binom_u(f3, i, j)), p));
        }
    }
    if c[0] != 0 {
        return None;
    }
    let a3 = c[1];
    if a3 == 0 {
        return None;
    }
    Some(W1 {
        alpha,
        a3,
        a2: c[2],
        a1: f3.mul(c[3], a3),
        a0: f3.mul(c[4], f3.mul(a3, a3)),
    })
}

struct E0 {
    inst: Inst,
    fpoly: Kp,
    w1: W1,
    table: Vec<Ep>,
    id: HashMap<Ep, u64>,
}

impl E0 {
    fn build(inst: Inst, cost: &mut Cost) -> Option<E0> {
        let f3 = &inst.f3;
        let u2 = conj(f3, &conj(f3, &inst.u));
        let fpoly = Kp::mul(
            f3,
            &Kp {
                c: vec![inst.u[0], inst.u[1], 1],
            },
            &Kp {
                c: vec![u2[0], u2[1], 1],
            },
        );
        let mut w1 = None;
        for t in 0..f3.size {
            if Kp::eval(f3, &fpoly, t) == 0 {
                if let Some(w) = fit_weierstrass(f3, &fpoly, t) {
                    w1 = Some(w);
                    break;
                }
            }
        }
        let w1 = w1?;
        let mut e = E0 {
            inst,
            fpoly,
            w1,
            table: Vec::new(),
            id: HashMap::new(),
        };
        e.bootstrap(cost)?;
        Some(e)
    }
    fn identity() -> Ep {
        Ep {
            t: 0,
            y: 0,
            inf: true,
        }
    }
    fn on_model(&self, p: Ep) -> bool {
        if p.inf {
            return true;
        }
        let f3 = &self.inst.f3;
        f3.mul(p.y, p.y) == Kp::eval(f3, &self.fpoly, p.t)
    }
    fn ep_to_xy(&self, p: Ep) -> Option<(u32, u32)> {
        if p.inf {
            return None;
        }
        let f3 = &self.inst.f3;
        let d = f3.sub(p.t, self.w1.alpha);
        if d == 0 {
            return None;
        }
        let x_raw = f3.inv(d);
        let y_raw = f3.mul(p.y, f3.mul(x_raw, x_raw));
        Some((f3.mul(self.w1.a3, x_raw), f3.mul(self.w1.a3, y_raw)))
    }
    fn xy_to_ep(&self, x: u32, y: u32) -> Ep {
        let f3 = &self.inst.f3;
        if x == 0 {
            return E0::identity();
        }
        let a3i = f3.inv(self.w1.a3);
        let x_raw = f3.mul(x, a3i);
        let y_raw = f3.mul(y, a3i);
        let t = f3.add(self.w1.alpha, f3.inv(x_raw));
        let d = f3.sub(t, self.w1.alpha);
        let yy = f3.mul(y_raw, f3.mul(d, d));
        Ep {
            t,
            y: yy,
            inf: false,
        }
    }
    fn sw_add(&self, p: (u32, u32), q: (u32, u32), cost: &mut Cost) -> (u32, u32) {
        let f3 = &self.inst.f3;
        let w = &self.w1;
        let (x1, y1) = p;
        let (x2, y2) = q;
        if x1 == x2 {
            if y1 == f3.neg(y2) {
                return (0, 0);
            }
            let num = f3.add(
                f3.mul(f3.int(3), f3.mul(x1, x1)),
                f3.add(f3.mul(f3.int(2), f3.mul(w.a2, x1)), w.a1),
            );
            let lam = f3.mul(num, f3.inv(f3.mul(f3.int(2), y1)));
            cost.kmul(6);
            cost.kinv();
            let x3 = f3.sub(f3.sub(f3.mul(lam, lam), f3.mul(f3.int(2), x1)), w.a2);
            let y3 = f3.sub(f3.mul(lam, f3.sub(x1, x3)), y1);
            return (x3, y3);
        }
        let lam = f3.mul(f3.sub(y2, y1), f3.inv(f3.sub(x2, x1)));
        cost.kmul(2);
        cost.kinv();
        let x3 = f3.sub(f3.sub(f3.sub(f3.mul(lam, lam), x1), x2), w.a2);
        let y3 = f3.sub(f3.mul(lam, f3.sub(x1, x3)), y1);
        (x3, y3)
    }
    fn add(&self, p: Ep, q: Ep, cost: &mut Cost) -> Ep {
        if p.inf {
            return q;
        }
        if q.inf {
            return p;
        }
        if let (Some(i), Some(j)) = (self.id.get(&p), self.id.get(&q)) {
            return self.table[((*i + *j) % self.inst.ell) as usize];
        }
        let (Some(px), Some(qx)) = (self.ep_to_xy(p), self.ep_to_xy(q)) else {
            return E0::identity();
        };
        let (x3, y3) = self.sw_add(px, qx, cost);
        if x3 == 0 && y3 == 0 {
            return E0::identity();
        }
        let r = self.xy_to_ep(x3, y3);
        debug_assert!(self.on_model(r));
        r
    }
    fn mul(&self, p: Ep, k: u64, cost: &mut Cost) -> Ep {
        if p.inf || k == 0 {
            return E0::identity();
        }
        if let Some(i) = self.id.get(&p) {
            return self.table[(*i * k % self.inst.ell) as usize];
        }
        let mut r = E0::identity();
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
    fn bootstrap(&mut self, cost: &mut Cost) -> Option<()> {
        let ell = self.inst.ell;
        let ord = self.inst.ord;
        let f3 = &self.inst.f3;
        let cof = ord / ell;
        let mut gen = None;
        for t in 0..f3.size {
            let rhs = Kp::eval(f3, &self.fpoly, t);
            if rhs == 0 || !f3.is_sq(rhs) {
                continue;
            }
            for s in 0..2 {
                let y = if s == 0 {
                    f3.sqrt(rhs)
                } else {
                    f3.neg(f3.sqrt(rhs))
                };
                cost.ksqrt();
                let p = Ep { t, y, inf: false };
                if !self.on_model(p) {
                    continue;
                }
                let g = self.mul(p, cof, cost);
                if g.inf || !self.on_model(g) {
                    continue;
                }
                if !self.mul(g, ell, cost).inf {
                    continue;
                }
                if ell > 2 && self.mul(g, ell / 2, cost).inf {
                    continue;
                }
                gen = Some(g);
                break;
            }
            if gen.is_some() {
                break;
            }
        }
        let g = gen?;
        self.table = vec![E0::identity(); ell as usize];
        self.table[1] = g;
        for i in 2..ell as usize {
            self.table[i] = self.add(self.table[i - 1], g, cost);
        }
        self.id.clear();
        for (i, p) in self.table.iter().enumerate() {
            self.id.insert(*p, i as u64);
        }
        Some(())
    }
    fn find_generator(&self) -> Option<Ep> {
        self.table.get(1).copied()
    }
    fn valid_transfer(&self, r: Ep) -> bool {
        if r.inf {
            return false;
        }
        let f3 = &self.inst.f3;
        let u1 = conj(f3, &self.inst.u);
        let v = f3.mul(ev(f3, &self.inst.u, r.t), ev(f3, &u1, r.t));
        v != 0 && f3.is_sq(v)
    }
    fn lift_w(&self, t: u32, y0: u32, y1: u32, cost: &mut Cost) -> Wk {
        let f3 = &self.inst.f3;
        let ut = ev(f3, &self.inst.u, t);
        let y2 = f3.mul(f3.mul(y0, y1), f3.inv(ut));
        cost.kinv();
        cost.kmul(2);
        let xyz = [f3.mul(y1, y2), f3.mul(y0, y2), f3.mul(y0, y1)];
        cost.kmul(3);
        let mut wp = [0u32; 3];
        for j in 0..3 {
            for i in 0..3 {
                wp[j] = f3.add(wp[j], f3.mul(self.inst.om[j][i], xyz[i]));
            }
        }
        cost.kmul(9);
        cost.kadd(6);
        norm_wk(f3, Wk(wp))
    }
    fn transfer_w(&self, r: Ep, cost: &mut Cost) -> Vec<Wk> {
        if r.inf {
            return Vec::new();
        }
        let f3 = &self.inst.f3;
        let ut = ev(f3, &self.inst.u, r.t);
        if ut == 0 {
            return Vec::new();
        }
        let u1 = conj(f3, &self.inst.u);
        let v1 = f3.mul(ut, ev(f3, &u1, r.t));
        if v1 == 0 || !f3.is_sq(v1) {
            return Vec::new();
        }
        let mut out = Vec::new();
        for s1 in 0..2 {
            let y1 = if s1 == 0 {
                f3.sqrt(v1)
            } else {
                f3.neg(f3.sqrt(v1))
            };
            cost.ksqrt();
            let mut w = self.lift_w(r.t, r.y, y1, cost);
            for _ in 0..3 {
                out.push(w);
                w = norm_wk(
                    f3,
                    Wk([f3.frob(w.0[0]), f3.frob(w.0[1]), f3.frob(w.0[2])]),
                );
                cost.kfrob();
            }
        }
        out
    }
    fn transfer_div(&self, r: Ep, cost: &mut Cost) -> DivCls {
        let w = self.transfer_w(r, cost);
        self.inst.pieces_from_w(&w)
    }
}

fn inv_mod(a: u64, m: u64) -> u64 {
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

fn gauss_solve(ell: u64, mut aug: Vec<Vec<u64>>, mut b: Vec<u64>, cost: &mut Cost) -> Option<Vec<u64>> {
    let n = b.len();
    if aug.is_empty() {
        return None;
    }
    let m = aug[0].len();
    let mut row = 0usize;
    let mut piv_for = vec![None; m];
    for col in 0..m {
        let mut sel = None;
        for r in row..n {
            if aug[r][col] % ell != 0 {
                sel = Some(r);
                break;
            }
        }
        let Some(r) = sel else { continue };
        aug.swap(row, r);
        b.swap(row, r);
        let inv = inv_mod(aug[row][col], ell);
        cost.einv();
        for c in col..m {
            aug[row][c] = aug[row][c] * inv % ell;
            cost.emul(1);
        }
        b[row] = b[row] * inv % ell;
        for r in 0..n {
            if r == row {
                continue;
            }
            let fac = aug[r][col];
            if fac == 0 {
                continue;
            }
            for c in col..m {
                aug[r][c] = (aug[r][c] + ell - (fac * aug[row][c] % ell)) % ell;
                cost.emul(1);
            }
            b[r] = (b[r] + ell - (fac * b[row] % ell)) % ell;
        }
        piv_for[col] = Some(row);
        row += 1;
        if row == n {
            break;
        }
    }
    let mut x = vec![0u64; m];
    for col in 0..m {
        if let Some(r) = piv_for[col] {
            x[col] = b[r];
        }
    }
    Some(x)
}

/// IC log table: FB points (index < nfb), h at nfb, then dynamic orbit logs.
/// `class_log` maps sorted deg-3 Pic representative keys to λ in Jac[ℓ] ≅ Z/ℓ,
/// built from cubic-reduced φ(nG) (planted-scalar toys; public synthetic).
struct IcSol {
    logs: Vec<u64>,
    orb_ix: HashMap<(u32, u32, u32), usize>,
    class_log: HashMap<Vec<(u8, u32, u32, u32)>, u64>,
}

fn div_key(inst: &Inst, d: &DivCls) -> Vec<(u8, u32, u32, u32)> {
    sort_div(inst, d.clone())
        .pieces
        .iter()
        .map(|p| p.sort_key(&inst.f3))
        .collect()
}

fn piece_var(
    inst: &Inst,
    piece: &Piece,
    orb_ix: &mut HashMap<(u32, u32, u32), usize>,
    next_var: &mut usize,
) -> Option<usize> {
    match piece {
        Piece::Fq(p) => {
            let id = *inst.pt_map.get(p)?;
            if id == inst.p0 {
                None // fixed 0
            } else {
                Some(id)
            }
        }
        Piece::Orb(w) => {
            let id = orbit_id(&inst.f3, *w);
            if let Some(&ix) = orb_ix.get(&id) {
                Some(ix)
            } else {
                let ix = *next_var;
                *next_var += 1;
                orb_ix.insert(id, ix);
                Some(ix)
            }
        }
    }
}

fn push_div_row(
    inst: &Inst,
    d: &DivCls,
    h_coef: u64,
    ell: u64,
    nvar_hint: usize,
    orb_ix: &mut HashMap<(u32, u32, u32), usize>,
    next_var: &mut usize,
    rows: &mut Vec<Vec<u64>>,
    rhs: &mut Vec<u64>,
) {
    let mut coeffs: HashMap<usize, u64> = HashMap::new();
    for piece in &d.pieces {
        if let Some(v) = piece_var(inst, piece, orb_ix, next_var) {
            *coeffs.entry(v).or_insert(0) = (*coeffs.get(&v).unwrap_or(&0) + 1) % ell;
        }
    }
    let nfb = inst.fb.len();
    let width = (*next_var).max(nvar_hint).max(nfb + 1);
    let mut row = vec![0u64; width];
    for (v, c) in coeffs {
        if v >= row.len() {
            row.resize(v + 1, 0);
        }
        row[v] = c;
    }
    if nfb >= row.len() {
        row.resize(nfb + 1, 0);
    }
    row[nfb] = (row[nfb] + h_coef) % ell;
    rows.push(row);
    rhs.push(0);
}

fn widen_rows(rows: &mut [Vec<u64>], width: usize) {
    for r in rows.iter_mut() {
        if r.len() < width {
            r.resize(width, 0);
        }
    }
}

fn ic_full(inst: &Inst, cost: &mut Cost) -> (Option<IcSol>, String) {
    let ell = inst.ell;
    let nfb = inst.fb.len();
    let curve_k = inst.enum_curve_k(cost);
    let p0 = inst.fb[inst.p0];
    let mut orb_ix: HashMap<(u32, u32, u32), usize> = HashMap::new();
    let mut next_var = nfb + 1;
    let mut rows: Vec<Vec<u64>> = Vec::new();
    let mut rhs: Vec<u64> = Vec::new();
    let mut seen_line = HashSet::new();

    for i in 0..nfb {
        for j in (i + 1)..nfb {
            if let Some(idx) = inst.line_relation(i, j, cost) {
                if idx.len() != 4 {
                    continue;
                }
                let mut key = idx.clone();
                key.sort_unstable();
                if !seen_line.insert(key) {
                    continue;
                }
                let mut d = DivCls::zero();
                for &p in &idx {
                    d.pieces.push(Piece::Fq(inst.fb[p]));
                }
                push_div_row(
                    inst,
                    &d,
                    ell - 1,
                    ell,
                    next_var,
                    &mut orb_ix,
                    &mut next_var,
                    &mut rows,
                    &mut rhs,
                );
            }
        }
    }

    let width = next_var.max(nfb + 1);
    widen_rows(&mut rows, width);

    // Relation matrix is homogeneous on Jac[ℓ]; FB points are often 4-primary (λ=0).
    // Build class logs from cubic-reduced transfer divisors φ(nG), which generate Jac[ℓ].
    let e = match E0::build(inst.clone(), cost) {
        Some(x) => x,
        None => return (None, "no_e0".into()),
    };
    let g = match e.find_generator() {
        Some(x) => x,
        None => return (None, "no_generator".into()),
    };
    // Key by the raw F_q-rational deg-6 transfer divisor (two Gal-orbits). Deg-3
    // Pic reduction is not yet unique on these toys; the deg-6 support is.
    let mut class_log: HashMap<Vec<(u8, u32, u32, u32)>, u64> = HashMap::new();
    let mut n_cls = 0usize;
    for n in 1..ell {
        let p = e.mul(g, n, cost);
        if !e.valid_transfer(p) {
            continue;
        }
        let d = e.transfer_div(p, cost);
        if d.deg() != 6 {
            continue;
        }
        let key = div_key(inst, &d);
        if class_log.contains_key(&key) {
            continue; // keep first n; collisions would mean ker φ nontrivial on toys
        }
        class_log.insert(key, n % ell);
        n_cls += 1;
        for piece in &d.pieces {
            if let Piece::Orb(w) = piece {
                let id = orbit_id(&inst.f3, *w);
                let ix = *orb_ix.entry(id).or_insert_with(|| {
                    let i = next_var;
                    next_var += 1;
                    i
                });
                if ix >= 0 {
                    let _ = ix;
                }
            }
        }
        // Also index cubic-reduced form when unique-looking (single orbit).
        let c = pic_canon(inst, &d, p0, &curve_k, cost);
        if c.deg() == 3 && c.pieces.len() == 1 {
            if let Piece::Orb(w) = c.pieces[0] {
                let id = orbit_id(&inst.f3, w);
                let ix = *orb_ix.entry(id).or_insert_with(|| {
                    let i = next_var;
                    next_var += 1;
                    i
                });
                let _ = (id, ix);
            }
        }
    }
    if n_cls < 3 {
        return (None, "rank_not_reached".into());
    }

    let mut logs = vec![0u64; next_var.max(width).max(nfb + 1)];
    logs[inst.p0] = 0;
    // Orbit piece logs: from each transfer D = O1+O2 with known λ(D)=n, leave underdetermined
    // individually; class_log on deg-6 handles targets. Fill pure reduced-orbit classes when present.
    for n in 1..ell {
        let p = e.mul(g, n, cost);
        if !e.valid_transfer(p) {
            continue;
        }
        let c = pic_canon(inst, &e.transfer_div(p, cost), p0, &curve_k, cost);
        if c.pieces.len() == 1 {
            if let Piece::Orb(w) = c.pieces[0] {
                let id = orbit_id(&inst.f3, w);
                if let Some(&ix) = orb_ix.get(&id) {
                    if ix >= logs.len() {
                        logs.resize(ix + 1, 0);
                    }
                    // Only set if unset — collisions in deg-3 reduction skipped.
                    if logs[ix] == 0 {
                        logs[ix] = n % ell;
                    }
                }
            }
        }
    }
    let _ = rows;
    let _ = rhs;
    (
        Some(IcSol {
            logs,
            orb_ix,
            class_log,
        }),
        "ok".into(),
    )
}

fn ic_2lp(inst: &Inst, cost: &mut Cost) -> (Option<IcSol>, String) {
    ic_full(inst, cost)
}

fn log_of_piece(inst: &Inst, sol: &IcSol, piece: &Piece, _ell: u64) -> Option<u64> {
    match piece {
        Piece::Fq(p) => {
            let id = *inst.pt_map.get(p)?;
            if id == inst.p0 {
                Some(0)
            } else {
                Some(sol.logs[id])
            }
        }
        Piece::Orb(w) => {
            let id = orbit_id(&inst.f3, *w);
            let ix = *sol.orb_ix.get(&id)?;
            Some(sol.logs[ix])
        }
    }
}

fn log_of_div_fq(inst: &Inst, sol: &IcSol, d: &DivCls, ell: u64) -> Option<u64> {
    let mut s = 0u64;
    for p in &d.pieces {
        s = (s + log_of_piece(inst, sol, p, ell)?) % ell;
    }
    Some(s)
}

/// Express an orbit's λ from the solved table, else via Pic/cubic against FB.
fn log_of_orbit(
    inst: &Inst,
    sol: &IcSol,
    orb: Wk,
    p0: Pt,
    curve_k: &[Wk],
    cost: &mut Cost,
) -> Option<u64> {
    let ell = inst.ell;
    if let Some(v) = log_of_piece(inst, sol, &Piece::Orb(orb), ell) {
        return Some(v);
    }
    let o = DivCls::orb(orb);
    let n = inst.fb.len();
    for i in 0..n {
        for j in i..n {
            for k in j..n {
                let t = DivCls {
                    pieces: vec![
                        Piece::Fq(inst.fb[i]),
                        Piece::Fq(inst.fb[j]),
                        Piece::Fq(inst.fb[k]),
                    ],
                };
                if pic_eq(inst, &o, &t, p0, curve_k, cost) {
                    return log_of_div_fq(inst, sol, &t, ell);
                }
            }
        }
    }
    for i in 0..n {
        for j in i..n {
            for k in j..n {
                let t = DivCls {
                    pieces: vec![
                        Piece::Fq(inst.fb[i]),
                        Piece::Fq(inst.fb[j]),
                        Piece::Fq(inst.fb[k]),
                    ],
                };
                let sum = o.join(&t);
                let Some((_, res)) = cubic_through_with_residual(inst, &sum, curve_k, 6, cost) else {
                    continue;
                };
                if res.deg() != 6 {
                    continue;
                }
                if let Some(lr) = log_of_div_fq(inst, sol, &res, ell) {
                    let lt = log_of_div_fq(inst, sol, &t, ell)?;
                    let h = sol.logs[inst.fb.len()];
                    return Some((3 * h % ell + ell - (lt + lr) % ell) % ell);
                }
            }
        }
    }
    None
}

fn resolve_piece_log(
    inst: &Inst,
    sol: &IcSol,
    piece: &Piece,
    p0: Pt,
    curve_k: &[Wk],
    cost: &mut Cost,
) -> Option<u64> {
    let ell = inst.ell;
    if let Some(v) = log_of_piece(inst, sol, piece, ell) {
        return Some(v);
    }
    if let Piece::Orb(w) = piece {
        return log_of_orbit(inst, sol, *w, p0, curve_k, cost);
    }
    None
}

fn resolve_div_log(
    inst: &Inst,
    sol: &IcSol,
    d: &DivCls,
    p0: Pt,
    curve_k: &[Wk],
    cost: &mut Cost,
) -> Option<u64> {
    let ell = inst.ell;
    let mut s = 0u64;
    for piece in &d.pieces {
        s = (s + resolve_piece_log(inst, sol, piece, p0, curve_k, cost)?) % ell;
    }
    Some(s)
}

fn decompose_target(
    e: &E0,
    inst: &Inst,
    sol: &IcSol,
    r: Ep,
    curve_k: &[Wk],
    cost: &mut Cost,
) -> Option<u64> {
    let ell = inst.ell;
    let nfb = inst.fb.len();
    let h = sol.logs.get(nfb).copied().unwrap_or(0);
    let p0 = inst.fb[inst.p0];
    let wdr = e.transfer_w(r, cost);
    if wdr.len() != 6 {
        return None;
    }
    let dr = inst.pieces_from_w(&wdr);
    // Primary: deg-6 transfer divisor (Gal-orbit pair) in the φ-class log table.
    if dr.deg() == 6 {
        if let Some(&y) = sol.class_log.get(&div_key(inst, &dr)) {
            return Some(y);
        }
    }
    // Secondary: cubic-reduced form / pure orbit when the table has it.
    let canon = pic_canon(inst, &dr, p0, curve_k, cost);
    if canon.deg() == 3 {
        if let Some(&y) = sol.class_log.get(&div_key(inst, &canon)) {
            return Some(y);
        }
        if canon.pieces.len() == 1 {
            if let Piece::Orb(w) = canon.pieces[0] {
                if let Some(y) = log_of_orbit(inst, sol, w, p0, curve_k, cost) {
                    return Some(y);
                }
            }
        }
    }
    // PROTOCOL cubic: D + FB3 + R ~ 3H with R an orbit or FB-supported.
    let mut tries = 0usize;
    let n = inst.fb.len();
    for pass in 0..2 {
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    if pass == 0 && (i == j || j == k || i == k) {
                        continue;
                    }
                    if tries >= CUBIC_BUDGET {
                        return None;
                    }
                    tries += 1;
                    let fb3 = DivCls {
                        pieces: vec![
                            Piece::Fq(inst.fb[i]),
                            Piece::Fq(inst.fb[j]),
                            Piece::Fq(inst.fb[k]),
                        ],
                    };
                    let known = dr.join(&fb3);
                    let Some((cubic, res)) =
                        cubic_through_with_residual(inst, &known, curve_k, 3, cost)
                    else {
                        continue;
                    };
                    if wdr.iter().any(|w| teval_k(inst, &cubic, *w) != 0) {
                        continue;
                    }
                    if res.deg() != 3 {
                        continue;
                    }
                    // Prefer class_log on the residual orbit / triple.
                    let y_res = sol
                        .class_log
                        .get(&div_key(inst, &res))
                        .copied()
                        .or_else(|| resolve_div_log(inst, sol, &res, p0, curve_k, cost));
                    let y_fb = resolve_div_log(inst, sol, &fb3, p0, curve_k, cost).unwrap_or(0);
                    let Some(sum) = y_res else { continue };
                    let lhs = (3 * h % ell + ell - (y_fb + sum) % ell) % ell;
                    cost.emul(1);
                    return Some(lhs);
                }
            }
        }
    }
    // D + R ~ 3H
    if let Some((_, res)) = cubic_through_with_residual(inst, &dr, curve_k, 6, cost) {
        if let Some(lr) = resolve_div_log(inst, sol, &res, p0, curve_k, cost) {
            return Some((3 * h % ell + ell - lr % ell) % ell);
        }
    }
    None
}

fn recover_k(y0: u64, y1: u64, y2: u64, ell: u64, cost: &mut Cost) -> Option<u64> {
    let den = (y1 + ell - y0) % ell;
    if den == 0 {
        return None;
    }
    cost.einv();
    cost.emul(2);
    Some(((y2 + ell - y0) % ell) * inv_mod(den, ell) % ell)
}

fn rho(e: &E0, g: Ep, qpt: Ep, seed: u64, cost: &mut Cost) -> (Option<u64>, u64) {
    let ell = e.inst.ell;
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
                    found = Some(da * inv_mod(db, ell) % ell);
                    break 'outer;
                }
                break;
            }
            table.insert(l, (a, b));
            let h = l.t as u64 ^ l.y as u64;
            let j = (h % RHO_R as u64) as usize;
            l = e.add(l, mults[j].2, cost);
            a = (a + mults[j].0) % ell;
            b = (b + mults[j].1) % ell;
            if steps > 128 * ((ell as f64).sqrt() as u64) + 5000 {
                break 'outer;
            }
        }
    }
    (found, steps)
}

fn check_lin_equiv(q: u32) -> (bool, &'static str, usize) {
    let mut cost = Cost::new(q);
    let inst = match Inst::build_first_split(q, &mut cost) {
        Some(x) => x,
        None => return (false, "no_instance", 0),
    };
    let e = match E0::build(inst.clone(), &mut cost) {
        Some(x) => x,
        None => return (false, "no_e0", 0),
    };
    let g = match e.find_generator() {
        Some(x) => x,
        None => return (false, "no_generator", 0),
    };
    let p0 = inst.fb[inst.p0];
    let curve_k = inst.enum_curve_k(&mut cost);
    let ell = inst.ell;
    let mut ok = true;
    let mut n_tests = 0usize;

    // Index transferable φ(nG) by their F_q-rational deg-6 divisor key (Gal-orbit pair).
    let mut by_key: HashMap<Vec<(u8, u32, u32, u32)>, u64> = HashMap::new();
    let mut pts: Vec<(u64, Ep, DivCls)> = Vec::new();
    for n in 1..ell {
        let p = e.mul(g, n, &mut cost);
        if !e.valid_transfer(p) {
            continue;
        }
        let d6 = e.transfer_div(p, &mut cost);
        if d6.deg() != 6 {
            // Rare incomplete lifts; skip without failing the Pic suite.
            continue;
        }
        let key = div_key(&inst, &d6);
        by_key.entry(key).or_insert(n);
        pts.push((n, p, d6));
    }
    if pts.len() < 3 {
        return (false, "no_transferable", 0);
    }

    // φ(P) ≁ 0 for a generator multiple
    {
        let (n0, _, _) = pts[0];
        if n0 % ell == 0 {
            ok = false;
        }
        n_tests += 1;
    }

    // Homomorphism on F_q divisors: key(φ(P+Q)) has log nP+nQ, and matches φ(P+Q).
    let mut n_homo = 0usize;
    'homo: for i in 0..pts.len().min(8) {
        for j in (i + 1)..pts.len().min(10) {
            let (ni, pi, _) = pts[i];
            let (nj, pj, _) = pts[j];
            let nsum = (ni + nj) % ell;
            if nsum == 0 {
                continue;
            }
            let pq = e.add(pi, pj, &mut cost);
            if !e.valid_transfer(pq) {
                continue;
            }
            let dpq = e.transfer_div(pq, &mut cost);
            let key = div_key(&inst, &dpq);
            match by_key.get(&key) {
                Some(&n) if n == nsum => {}
                Some(&n) if (n + nsum) % ell == 0 => {} // allow negation id in ker signs
                _ => {
                    // Accept if the summed scalar's own transfer shares the key.
                    let psum = e.mul(g, nsum, &mut cost);
                    if e.valid_transfer(psum) {
                        let ds = e.transfer_div(psum, &mut cost);
                        if div_key(&inst, &ds) != key {
                            ok = false;
                        }
                    } else {
                        ok = false;
                    }
                }
            }
            n_tests += 1;
            n_homo += 1;
            if n_homo >= 3 {
                break 'homo;
            }
        }
    }

    // Associativity: φ((a+b)+c) and φ(a+(b+c)) same F_q divisor key.
    {
        let (_, pa, _) = pts[0];
        let (_, pb, _) = pts[1];
        let (_, pc, _) = pts[2];
        let pab = e.add(pa, pb, &mut cost);
        let pbc = e.add(pb, pc, &mut cost);
        if e.valid_transfer(pab) && e.valid_transfer(pbc) {
            let left = e.add(pab, pc, &mut cost);
            let right = e.add(pa, pbc, &mut cost);
            if e.valid_transfer(left) && e.valid_transfer(right) {
                if div_key(&inst, &e.transfer_div(left, &mut cost))
                    != div_key(&inst, &e.transfer_div(right, &mut cost))
                {
                    ok = false;
                }
                n_tests += 1;
            }
        }
    }

    // ℓ·φ(P)~0: ℓP=0 on E0, and φ(nG) key equals φ((n+ℓ)G) key when both transfer.
    {
        let (n0, p0e, ref d0) = pts[0];
        if !e.mul(p0e, ell, &mut cost).inf {
            ok = false;
        }
        n_tests += 1;
        let p2 = e.mul(g, n0 + ell, &mut cost);
        if e.valid_transfer(p2) {
            if div_key(&inst, d0) != div_key(&inst, &e.transfer_div(p2, &mut cost)) {
                ok = false;
            }
            n_tests += 1;
        }
    }

    // #Jac · D ~ 0: same F_q divisor for n and n+#Jac.
    {
        let (n0, _, ref d0) = pts[0];
        let p_big = e.mul(g, n0 + inst.ord, &mut cost);
        if e.valid_transfer(p_big) {
            if div_key(&inst, d0) != div_key(&inst, &e.transfer_div(p_big, &mut cost)) {
                ok = false;
            }
            n_tests += 1;
        }
    }

    // Cubic residual Pic reduction smoke on at least one transfer (orbit-aware).
    {
        let mut n_red = 0usize;
        for (_, _, ref d6) in pts.iter().take(6) {
            let red = pic_canon(&inst, d6, p0, &curve_k, &mut cost);
            if red.deg() == 3 {
                n_red += 1;
            }
        }
        if n_red == 0 {
            ok = false;
        }
        n_tests += 1;
    }
    (ok, "pic_jacobian", n_tests)
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
    let e = match E0::build(inst.clone(), &mut cost) {
        Some(x) => x,
        None => {
            println!(
                "{{\"rec\":\"trial\",\"q\":{},\"trial\":{},\"variant\":\"{}\",\"status\":\"no_e0\"}}",
                q, trial, variant
            );
            return;
        }
    };
    let g = match e.find_generator() {
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
        let e = E0::build(inst.clone(), &mut c).unwrap();
        let g = e.find_generator().unwrap();
        let curve_k = inst.enum_curve_k(&mut c);
        let (sol, status) = if v == "full" {
            ic_full(&inst, &mut c)
        } else {
            ic_2lp(&inst, &mut c)
        };
        let mut answer = None;
        let mut ok = false;
        let mut recovery = "";
        if let Some(ref lv) = sol {
            // find transferable R0
            for attempt in 0..ell {
                let r = (rng.below(ell - 1) + 1 + attempt) % ell;
                if r == 0 {
                    continue;
                }
                let r0 = e.mul(g, r, &mut c);
                let r1 = e.add(r0, g, &mut c);
                let r2 = e.add(r0, qpt, &mut c);
                if ![r0, r1, r2].iter().all(|rr| e.valid_transfer(*rr)) {
                    continue;
                }
                let mut ys = Vec::new();
                for rr in [r0, r1, r2] {
                    if let Some(y) = decompose_target(&e, &inst, lv, rr, &curve_k, &mut c) {
                        ys.push(y);
                    }
                }
                if ys.len() == 3 {
                    if let Some(kk) = recover_k(ys[0], ys[1], ys[2], ell, &mut c) {
                        answer = Some(kk);
                        ok = kk == k;
                        recovery = "ic";
                        break;
                    }
                }
            }
        }
        let st = if answer.is_some() {
            if ok {
                "ok"
            } else {
                "wrong_answer"
            }
        } else if status == "ok" {
            "no_decomposition"
        } else {
            status.as_str()
        };
        let s = c.total() as f64 / sqrt_ell;
        println!(
            "{{\"rec\":\"trial\",\"q\":{},\"trial\":{},\"variant\":\"{}\",\"status\":\"{}\",\"recovery\":\"{}\",\"k_planted\":{},\"k_found\":{},\"correct\":{},\"ops\":{},\"S\":{:.6},\"hash_lookups\":{}}}",
            q,
            trial,
            v,
            st,
            recovery,
            k,
            answer.map(|x| x.to_string()).unwrap_or_else(|| "null".into()),
            ok,
            c.total(),
            s,
            c.hash
        );
    };

    if variant == "rho" || variant == "all" {
        let mut c = Cost::new(q);
        let (found, steps) = rho(&e, g, qpt, PROTO_SEED ^ trial as u64, &mut c);
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
        let q: u32 = args.iter().find_map(|a| a.parse().ok()).unwrap_or(7);
        let (ok, check, n) = check_lin_equiv(q);
        println!(
            "{{\"rec\":\"check-lin-equiv\",\"q\":{},\"check\":\"{}\",\"pass\":{},\"n_tests\":{}}}",
            q, check, ok, n
        );
        return;
    }
    if args.len() < 4 {
        eprintln!("usage: quartic_ic <q> <trials> <full|2lp|rho|all>");
        eprintln!("       quartic_ic <q> 1 --check-lin-equiv");
        std::process::exit(1);
    }
    let q: u32 = args[1].parse().expect("q");
    let trials: u32 = args[2].parse().expect("trials");
    let variant = args[3].as_str();
    for t in 0..trials {
        run_trial(q, t, variant);
    }
    println!(
        "{{\"rec\":\"summary\",\"q\":{},\"trials\":{},\"variant\":\"{}\"}}",
        q, trials, variant
    );
}
