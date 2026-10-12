//! Benchmark: 4-dimensional GLV-GLS with an isogeny-loop endomorphism (D = -619 by default)
//! against GLS-2 and width-w NAF, on the quadratic GLS twist E'(F_{p^2}), p of 127 bits.
//!
//! All arms run in one binary on the same arithmetic: `Fp2<FpM<2>>` (Karatsuba, 3 base
//! multiplications), Jacobian doubling `dbl-2007-bl` (general a), mixed addition
//! `madd-2007-bl`, odd-multiple tables normalized to affine with one batched inversion per call.
//! The comparison is between scalar-multiplication strategies, not arithmetic styles.
//!
//! Arms (each timed call computes [k]P from (k, P) and returns a Jacobian point):
//! * `wnaf`: width-w NAF of k (about 251 bits);
//! * `gls2`: k = k0 + k1 mu (Babai against a Gauss-reduced basis), psi applied to P's table,
//!   interleaved w-NAF with shared doublings;
//! * `glvgls4`: k = k0 + k1 mu + k2 lambda + k3 lambda mu (Babai against an LLL-reduced basis),
//!   alpha'(P) by the twisted isogeny loop (affine, three inversions), tables for P and
//!   alpha'(P) with one batched inversion, psi applied to both, interleaved w-NAF;
//! * `wnaf (A/A)`: the baseline timed a second time; its ratio to the baseline is the noise floor;
//! * stage diagnostics (not speeds): the 4-dimensional decomposition alone, alpha'(P) alone.
//!
//! Before any timing, on every (k, P) pair: all arms agree (projective equality), psi(P) = [mu]P,
//! alpha'(P) = [lambda]P, every decomposition is congruent to k mod r; on the first 8 pairs the
//! baseline equals textbook double-and-add.  The run exits non-zero if any check fails.
//!
//! ```bash
//! cargo run --release --example gls_loop_bench -- [--d -619] [--pbits 127] [--seed 20261012] \
//!     [--pairs 256] [--rounds 20] [--w 5] [--json out.json]
//! ```
use isogeny_algos::api::PrimeFieldInt;
use isogeny_algos::bigint::Big;
use isogeny_algos::curve::{pmul, pmul_big, random_point_f, Curve, Pt};
use isogeny_algos::field::{Field, Rng};
use isogeny_algos::fp2::Fp2;
use isogeny_algos::fpm::FpM;
use isogeny_algos::gls_loop::{
    build_curve, cheapest_loop, closing_loops, find_instance, gauss2, lll4, Babai, FieldChain,
};
use isogeny_algos::int::Int;
use isogeny_algos::json::Json;
use std::time::Instant;

type B = FpM<2>;
type G = Fp2<B>;
type E = <G as Field>::E;

#[derive(Clone, Copy)]
struct Jac {
    x: E,
    y: E,
    z: E,
}

struct Ar<'a> {
    g: &'a G,
    a: E,
}

impl Ar<'_> {
    fn inf(&self) -> Jac {
        Jac {
            x: self.g.one(),
            y: self.g.one(),
            z: self.g.zero(),
        }
    }
    fn is_inf(&self, p: &Jac) -> bool {
        self.g.is_zero(p.z)
    }
    fn dbl(&self, p: &Jac) -> Jac {
        let g = self.g;
        if self.is_inf(p) {
            return *p;
        }
        let xx = g.sq(p.x);
        let yy = g.sq(p.y);
        let yyyy = g.sq(yy);
        let zz = g.sq(p.z);
        let t = g.sub(g.sub(g.sq(g.add(p.x, yy)), xx), yyyy);
        let s = g.add(t, t);
        let m = g.add(g.add(g.add(xx, xx), xx), g.mul(self.a, g.sq(zz)));
        let x3 = g.sub(g.sq(m), g.add(s, s));
        let y8 = {
            let y2 = g.add(yyyy, yyyy);
            let y4 = g.add(y2, y2);
            g.add(y4, y4)
        };
        let y3 = g.sub(g.mul(m, g.sub(s, x3)), y8);
        let z3 = g.sub(g.sub(g.sq(g.add(p.y, p.z)), yy), zz);
        Jac {
            x: x3,
            y: y3,
            z: z3,
        }
    }
    /// p + (x2, y2) with (x2, y2) affine.
    fn madd(&self, p: &Jac, x2: E, y2: E) -> Jac {
        let g = self.g;
        if self.is_inf(p) {
            return Jac {
                x: x2,
                y: y2,
                z: g.one(),
            };
        }
        let z1z1 = g.sq(p.z);
        let u2 = g.mul(x2, z1z1);
        let s2 = g.mul(y2, g.mul(p.z, z1z1));
        let h = g.sub(u2, p.x);
        let rr = {
            let t = g.sub(s2, p.y);
            g.add(t, t)
        };
        if g.is_zero(h) {
            return if g.is_zero(rr) {
                self.dbl(p)
            } else {
                self.inf()
            };
        }
        let hh = g.sq(h);
        let i = {
            let t = g.add(hh, hh);
            g.add(t, t)
        };
        let j = g.mul(h, i);
        let v = g.mul(p.x, i);
        let x3 = g.sub(g.sub(g.sq(rr), j), g.add(v, v));
        let yj = g.mul(p.y, j);
        let y3 = g.sub(g.mul(rr, g.sub(v, x3)), g.add(yj, yj));
        let z3 = g.sub(g.sub(g.sq(g.add(p.z, h)), z1z1), hh);
        Jac {
            x: x3,
            y: y3,
            z: z3,
        }
    }
    /// Full Jacobian addition (add-2007-bl), used for the tables.
    fn add(&self, p: &Jac, q: &Jac) -> Jac {
        let g = self.g;
        if self.is_inf(p) {
            return *q;
        }
        if self.is_inf(q) {
            return *p;
        }
        let z1z1 = g.sq(p.z);
        let z2z2 = g.sq(q.z);
        let u1 = g.mul(p.x, z2z2);
        let u2 = g.mul(q.x, z1z1);
        let s1 = g.mul(p.y, g.mul(q.z, z2z2));
        let s2 = g.mul(q.y, g.mul(p.z, z1z1));
        let h = g.sub(u2, u1);
        let rr = {
            let t = g.sub(s2, s1);
            g.add(t, t)
        };
        if g.is_zero(h) {
            return if g.is_zero(rr) {
                self.dbl(p)
            } else {
                self.inf()
            };
        }
        let i = g.sq(g.add(h, h));
        let j = g.mul(h, i);
        let v = g.mul(u1, i);
        let x3 = g.sub(g.sub(g.sq(rr), j), g.add(v, v));
        let sj = g.mul(s1, j);
        let y3 = g.sub(g.mul(rr, g.sub(v, x3)), g.add(sj, sj));
        let z3 = g.mul(g.sub(g.sub(g.sq(g.add(p.z, q.z)), z1z1), z2z2), h);
        Jac {
            x: x3,
            y: y3,
            z: z3,
        }
    }
    /// Odd multiples P, 3P, ..., (2^(w-1) - 1) P of several affine base points, all normalized
    /// to affine with one batched inversion.
    fn tables(&self, bases: &[(E, E)], w: u32) -> Vec<Vec<(E, E)>> {
        let g = self.g;
        let jb: Vec<Jac> = bases
            .iter()
            .map(|&(x, y)| Jac { x, y, z: g.one() })
            .collect();
        self.tables_jac(&jb, w)
    }
    /// The same with Jacobian bases (the base itself is normalized in the same batch).
    fn tables_jac(&self, bases: &[Jac], w: u32) -> Vec<Vec<(E, E)>> {
        let g = self.g;
        let n = 1usize << (w - 2);
        let mut jac: Vec<Jac> = vec![];
        for &p in bases {
            let p2 = self.dbl(&p);
            let mut cur = p;
            jac.push(cur);
            for _ in 1..n {
                cur = self.add(&cur, &p2);
                jac.push(cur);
            }
        }
        // batched inversion of all z
        let mut pref = Vec::with_capacity(jac.len());
        let mut acc = g.one();
        for q in &jac {
            pref.push(acc);
            acc = g.mul(acc, q.z);
        }
        let mut inv = g.inv(acc);
        let mut aff = vec![(g.zero(), g.zero()); jac.len()];
        for i in (0..jac.len()).rev() {
            let zi = g.mul(inv, pref[i]);
            inv = g.mul(inv, jac[i].z);
            let zi2 = g.sq(zi);
            aff[i] = (g.mul(jac[i].x, zi2), g.mul(jac[i].y, g.mul(zi2, zi)));
        }
        aff.chunks(n).map(|c| c.to_vec()).collect()
    }
    fn eq(&self, p: &Jac, q: &Jac) -> bool {
        let g = self.g;
        match (self.is_inf(p), self.is_inf(q)) {
            (true, true) => return true,
            (true, false) | (false, true) => return false,
            _ => {}
        }
        let z1z1 = g.sq(p.z);
        let z2z2 = g.sq(q.z);
        g.mul(p.x, z2z2) == g.mul(q.x, z1z1)
            && g.mul(p.y, g.mul(q.z, z2z2)) == g.mul(q.y, g.mul(p.z, z1z1))
    }
    fn lift_aff(&self, p: &Pt<E>) -> Jac {
        match *p {
            Pt::Inf => self.inf(),
            Pt::Aff(x, y) => Jac {
                x,
                y,
                z: self.g.one(),
            },
        }
    }
}

/// The twisted loop in homogeneous form, no inversion: per step, with V = u Z,
/// N = num(X, V), D = den(X, V) (homogenised), N1, D1 the homogenised derivatives:
/// X' = N, Z' = Z D (x' = u num(z)/den(z), z = x/u), Y' = Y (N1 D - N D1), W' = W D^2
/// (y' = y R'(z)); then the closing isomorphism; returned in Jacobian form.
struct ProjChain {
    steps: Vec<[Vec<E>; 4]>,
    w2_inv: E,
    w3_inv: E,
}

fn homog(g: &G, c: &[E], x: E, vpow: &[E]) -> E {
    let d = c.len() - 1;
    let mut acc = c[d];
    for i in (0..d).rev() {
        acc = g.add(g.mul(acc, x), g.mul(c[i], vpow[d - i]));
    }
    acc
}

impl ProjChain {
    fn new(g: &G, fc: &FieldChain<G>) -> ProjChain {
        let deriv = |c: &Vec<E>| -> Vec<E> {
            (1..c.len())
                .map(|i| g.mul(g.from_u64(i as u64), c[i]))
                .collect()
        };
        let steps = fc
            .maps
            .iter()
            .map(|(n, d)| {
                assert_eq!(
                    n.len(),
                    d.len() + 1,
                    "x-map degrees: num = den + 1 expected"
                );
                [n.clone(), d.clone(), deriv(n), deriv(d)]
            })
            .collect();
        ProjChain {
            steps,
            w2_inv: fc.w2_inv,
            w3_inv: fc.w3_inv,
        }
    }
    fn eval(&self, g: &G, u: E, x: E, y: E) -> Jac {
        let (mut xx, mut zz, mut yy, mut ww) = (x, g.one(), y, g.one());
        for [n, d, n1, d1] in &self.steps {
            let v = g.mul(u, zz);
            let mut vpow = Vec::with_capacity(n.len());
            vpow.push(g.one());
            for i in 1..n.len() {
                vpow.push(g.mul(vpow[i - 1], v));
            }
            let nn = homog(g, n, xx, &vpow);
            let dd = homog(g, d, xx, &vpow);
            let nn1 = homog(g, n1, xx, &vpow);
            let dd1 = homog(g, d1, xx, &vpow);
            yy = g.mul(yy, g.sub(g.mul(nn1, dd), g.mul(nn, dd1)));
            ww = g.mul(ww, g.sq(dd));
            xx = nn;
            zz = g.mul(zz, dd);
        }
        xx = g.mul(xx, self.w2_inv);
        yy = g.mul(yy, self.w3_inv);
        // x = X/Z, y = Y/W  ->  Jacobian with Z_J = Z W
        let zw = g.mul(zz, ww);
        let w2 = g.sq(ww);
        Jac {
            x: g.mul(g.mul(xx, zz), w2),
            y: g.mul(g.mul(yy, g.mul(g.sq(zz), zz)), w2),
            z: zw,
        }
    }
}

/// Width-w NAF digits, least significant first.
fn wnaf_i128(mut k: i128, w: u32) -> Vec<i32> {
    let mut out = vec![];
    let m = 1i128 << w;
    while k != 0 {
        if k & 1 == 1 {
            let mut d = k.rem_euclid(m);
            if d >= m / 2 {
                d -= m;
            }
            out.push(d as i32);
            k -= d;
        } else {
            out.push(0);
        }
        k >>= 1;
    }
    out
}

fn wnaf_big(k: &Int, w: u32) -> Vec<i32> {
    let mut k = k.clone();
    let mut out = vec![];
    let m = 1u64 << w;
    let two = Int::from(2i64);
    while !k.is_zero() {
        if !k.is_even() {
            let mut d = k.mod_u64(m) as i64;
            if d >= (m / 2) as i64 {
                d -= m as i64;
            }
            out.push(d as i32);
            k = &k - &Int::from(d);
        } else {
            out.push(0);
        }
        k = k.div_floor(&two);
    }
    out
}

/// Interleaved multi-scalar w-NAF with shared doublings; `tabs[i]` holds the odd multiples of
/// the i-th base, `neg[i]` negates that base.
fn interleaved(ar: &Ar, digits: &[Vec<i32>], tabs: &[&Vec<(E, E)>], neg: &[bool]) -> Jac {
    let g = ar.g;
    let len = digits.iter().map(|d| d.len()).max().unwrap_or(0);
    let mut r = ar.inf();
    for i in (0..len).rev() {
        r = ar.dbl(&r);
        for (j, d) in digits.iter().enumerate() {
            let Some(&dig) = d.get(i) else { continue };
            if dig == 0 {
                continue;
            }
            let (x, y) = tabs[j][(dig.unsigned_abs() as usize) / 2];
            let flip = (dig < 0) != neg[j];
            r = ar.madd(&r, x, if flip { g.neg(y) } else { y });
        }
    }
    r
}

fn pow_big(g: &G, a: E, e: &Big) -> E {
    let mut r = g.one();
    for i in (0..e.bits()).rev() {
        r = g.mul(r, r);
        if e.bit(i) {
            r = g.mul(r, a);
        }
    }
    r
}

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let d: i64 = arg(&args, "--d", -619);
    let pbits: usize = arg(&args, "--pbits", 127);
    let seed: u64 = arg(&args, "--seed", 20_261_012);
    let pairs: usize = arg(&args, "--pairs", 256);
    let rounds: usize = arg(&args, "--rounds", 20);
    let w: u32 = arg(&args, "--w", 5);
    let json_out: Option<String> = args
        .iter()
        .position(|a| a == "--json")
        .and_then(|i| args.get(i + 1))
        .cloned();

    // ---------------- setup (untimed)
    let (lp, _, _, _) = cheapest_loop(d);
    let avoid: Vec<u64> = lp
        .steps
        .iter()
        .map(|&(l, _)| l)
        .filter(|&l| l > 2)
        .collect();
    let inst =
        find_instance(d, pbits, &avoid, 64, seed ^ (1u64 << 32), 20_000_000).expect("instance");
    let cm = build_curve(&inst, seed ^ (1u64 << 32)).expect("curve");
    let chains = closing_loops(
        &cm.curve,
        &cm.crater,
        &lp.degree_sequence(),
        seed ^ (1u64 << 32),
    )
    .expect("walk");
    assert!(!chains.is_empty());
    let p = &inst.p;
    let r = &inst.r;
    let f = B::from_dec(&p.to_string());
    let g = G::new(f.clone());
    let mut rng = Rng::new(seed ^ 0xBE7C);
    let p2m1_half = p.mag().mul(p.mag()).sub_small(1).shr(1);
    let u: E = loop {
        let c = g.random(&mut rng);
        if !g.is_zero(c) && pow_big(&g, c, &p2m1_half) != g.one() {
            break c;
        }
    };
    let u2 = g.mul(u, u);
    let a_tw = g.mul(u2, g.embed(f.elem(&cm.curve.a)));
    let b_tw = g.mul(g.mul(u2, u), g.embed(f.elem(&cm.curve.b)));
    let ct = Curve::new(a_tw, b_tw);
    let ar = Ar { g: &g, a: a_tw };
    let c1 = g.mul(u, g.inv(g.conj(u)));
    let uh_inv = g.inv(pow_big(&g, u, &p.mag().sub_small(1).shr(1)));
    let c2 = g.mul(g.sq(uh_inv), uh_inv);
    let psi_aff = |(x, y): (E, E)| -> (E, E) { (g.mul(c1, g.conj(x)), g.mul(c2, g.conj(y))) };
    let chain: FieldChain<G> = chains[0].to_field(&g, &|v: &Int| g.embed(f.elem(v)));
    let alpha = |(x, y): (E, E)| -> Option<(E, E)> { chain.eval_twisted(&g, u, x, y) };

    // sample points of order r and determine mu, lambda (sign included) on the first one
    let pts: Vec<(E, E)> = (0..pairs)
        .map(|_| loop {
            if let Pt::Aff(x, y) = pmul(&g, &ct, &random_point_f(&g, &ct, &mut rng), inst.h as u128)
            {
                break (x, y);
            }
        })
        .collect();
    let p0 = Pt::Aff(pts[0].0, pts[0].1);
    let mu0 = (&(p - &Int::one()) * &inst.t.inv_mod(r).unwrap()).modulo(r);
    let mu = [mu0.clone(), (-&mu0).modulo(r)]
        .into_iter()
        .find(|m| {
            let (x, y) = psi_aff(pts[0]);
            pmul_big(&g, &ct, &p0, m.mag()) == Pt::Aff(x, y)
        })
        .expect("psi acts as +-mu");
    let s = (&(&(&Int::from(2i64) * &mu0) - &inst.t) * &inst.v.inv_mod(r).unwrap()).modulo(r);
    let inv2 = Int::from(2i64).inv_mod(r).unwrap();
    let mut lam_cands = vec![];
    for ss in [s.clone(), (-&s).modulo(r)] {
        let l = (&(&Int::from(lp.a) + &(&Int::from(lp.b) * &ss)) * &inv2).modulo(r);
        lam_cands.push(l.clone());
        lam_cands.push((-&l).modulo(r));
    }
    let a0 = alpha(pts[0]).expect("alpha(P) finite");
    let lam = lam_cands
        .into_iter()
        .find(|l| pmul_big(&g, &ct, &p0, l.mag()) == Pt::Aff(a0.0, a0.1))
        .expect("alpha' acts as one of +-lambda(+-s)");

    let b2 = Babai::new(gauss2(r, &mu), r.bits());
    let b4 = Babai::new(lll4(r, &mu, &lam), r.bits());
    let lm = (&lam * &mu).modulo(r);
    let rbig = r.mag().clone();
    let scalars: Vec<Int> = (0..pairs)
        .map(|_| {
            let mut v = Int::zero();
            for _ in 0..5 {
                v = &v.shl(64) + &Int::from(rng.next());
            }
            v.modulo(r)
        })
        .collect();

    // ---------------- the arms
    let arm_wnaf = |k: &Int, pt: (E, E)| -> Jac {
        let t = ar.tables(&[pt], w);
        interleaved(&ar, &[wnaf_big(k, w)], &[&t[0]], &[false])
    };
    let arm_gls2 = |k: &Int, pt: (E, E)| -> Jac {
        let c = b2.decompose(k);
        let t = ar.tables(&[pt], w);
        let tpsi: Vec<(E, E)> = t[0].iter().map(|&q| psi_aff(q)).collect();
        let d0 = c[0].to_i128().unwrap();
        let d1 = c[1].to_i128().unwrap();
        interleaved(
            &ar,
            &[wnaf_i128(d0.abs(), w), wnaf_i128(d1.abs(), w)],
            &[&t[0], &tpsi],
            &[d0 < 0, d1 < 0],
        )
    };
    let arm_glvgls4 = |k: &Int, pt: (E, E)| -> Jac {
        let c = b4.decompose(k);
        let ap = alpha(pt).expect("alpha(P) finite");
        let t = ar.tables(&[pt, ap], w);
        let tpsi: Vec<(E, E)> = t[0].iter().map(|&q| psi_aff(q)).collect();
        let tapsi: Vec<(E, E)> = t[1].iter().map(|&q| psi_aff(q)).collect();
        let ds: Vec<i128> = c.iter().map(|x| x.to_i128().unwrap()).collect();
        interleaved(
            &ar,
            &[
                wnaf_i128(ds[0].abs(), w),
                wnaf_i128(ds[1].abs(), w),
                wnaf_i128(ds[2].abs(), w),
                wnaf_i128(ds[3].abs(), w),
            ],
            &[&t[0], &tpsi, &t[1], &tapsi],
            &[ds[0] < 0, ds[1] < 0, ds[2] < 0, ds[3] < 0],
        )
    };

    let pchain = ProjChain::new(&g, &chain);
    let arm_glvgls4p = |k: &Int, pt: (E, E)| -> Jac {
        let c = b4.decompose(k);
        let ap = pchain.eval(&g, u, pt.0, pt.1);
        let t = ar.tables_jac(
            &[
                Jac {
                    x: pt.0,
                    y: pt.1,
                    z: g.one(),
                },
                ap,
            ],
            w,
        );
        let tpsi: Vec<(E, E)> = t[0].iter().map(|&q| psi_aff(q)).collect();
        let tapsi: Vec<(E, E)> = t[1].iter().map(|&q| psi_aff(q)).collect();
        let ds: Vec<i128> = c.iter().map(|x| x.to_i128().unwrap()).collect();
        interleaved(
            &ar,
            &[
                wnaf_i128(ds[0].abs(), w),
                wnaf_i128(ds[1].abs(), w),
                wnaf_i128(ds[2].abs(), w),
                wnaf_i128(ds[3].abs(), w),
            ],
            &[&t[0], &tpsi, &t[1], &tapsi],
            &[ds[0] < 0, ds[1] < 0, ds[2] < 0, ds[3] < 0],
        )
    };

    // ---------------- correctness, on every pair, before any timing
    let mut ok = true;
    let mut max_bits2 = 0usize;
    let mut max_bits4 = 0usize;
    for (i, (k, &pt)) in scalars.iter().zip(&pts).enumerate() {
        let base = arm_wnaf(k, pt);
        ok &= ar.eq(&base, &arm_gls2(k, pt));
        ok &= ar.eq(&base, &arm_glvgls4(k, pt));
        ok &= ar.eq(&base, &arm_glvgls4p(k, pt));
        ok &= ar.eq(
            &pchain.eval(&g, u, pt.0, pt.1),
            &ar.lift_aff(&pmul_big(&g, &ct, &Pt::Aff(pt.0, pt.1), lam.mag())),
        );
        let c2v = b2.decompose(k);
        let c4v = b4.decompose(k);
        ok &= (&(&c2v[0] + &(&c2v[1] * &mu)) - k).modulo(r).is_zero();
        let s4 = &(&(&c4v[0] + &(&c4v[1] * &mu)) + &(&c4v[2] * &lam)) + &(&c4v[3] * &lm);
        ok &= (&s4 - k).modulo(r).is_zero();
        max_bits2 = max_bits2.max(c2v.iter().map(|x| x.bits()).max().unwrap());
        max_bits4 = max_bits4.max(c4v.iter().map(|x| x.bits()).max().unwrap());
        let pp = Pt::Aff(pt.0, pt.1);
        let ps = psi_aff(pt);
        ok &= pmul_big(&g, &ct, &pp, mu.mag()) == Pt::Aff(ps.0, ps.1);
        let ap = alpha(pt).unwrap();
        ok &= pmul_big(&g, &ct, &pp, lam.mag()) == Pt::Aff(ap.0, ap.1);
        if i < 8 {
            ok &= ar.eq(&base, &ar.lift_aff(&pmul_big(&g, &ct, &pp, k.mag())));
            ok &= pmul_big(&g, &ct, &pp, &rbig) == Pt::Inf;
        }
    }
    eprintln!(
        "correctness on {pairs} pairs: {}",
        if ok { "all pass" } else { "FAIL" }
    );
    if !ok {
        std::process::exit(1);
    }

    // ---------------- timing: R rounds, arm order alternating, ns per call
    let names = [
        "wnaf",
        "gls2",
        "glvgls4 (affine loop)",
        "glvgls4 (projective loop)",
        "wnaf (A/A)",
        "decompose4 (stage)",
        "alpha'(P) affine (stage)",
        "alpha'(P) projective (stage)",
    ];
    let first_stage = 5;
    let mut times: Vec<Vec<f64>> = vec![vec![]; names.len()];
    let mut sink = 0u64;
    for round in 0..rounds {
        let order: Vec<usize> = if round % 2 == 0 {
            (0..names.len()).collect()
        } else {
            (0..names.len()).rev().collect()
        };
        for &arm in &order {
            let t0 = Instant::now();
            for (k, &pt) in scalars.iter().zip(&pts) {
                let v = match arm {
                    0 | 4 => arm_wnaf(k, pt).x.0[0],
                    1 => arm_gls2(k, pt).x.0[0],
                    2 => arm_glvgls4(k, pt).x.0[0],
                    3 => arm_glvgls4p(k, pt).x.0[0],
                    5 => b4.decompose(k)[0].bits() as u64,
                    6 => alpha(pt).unwrap().0 .0[0],
                    _ => pchain.eval(&g, u, pt.0, pt.1).x.0[0],
                };
                sink = sink.wrapping_add(v);
            }
            times[arm].push(t0.elapsed().as_nanos() as f64 / pairs as f64);
        }
    }
    let med: Vec<f64> = times.iter().map(|v| median(&mut v.clone())).collect();
    let min: Vec<f64> = times
        .iter()
        .map(|v| v.iter().cloned().fold(f64::INFINITY, f64::min))
        .collect();

    println!("# GLV-GLS with an isogeny-loop endomorphism: timing on E'(F_(p^2))\n");
    println!(
        "D = {d}, loop ({} + {} sqrt({d}))/2 (N = {}, steps {:?}, model cost {} M), p {} bits, r {} bits, h = {}, w = {w}, {pairs} pairs, {rounds} rounds, seed {seed}\n",
        lp.a, lp.b, lp.norm, lp.degree_sequence(), lp.cost, p.bits(), r.bits(), inst.h
    );
    println!("max decomposition coefficient bits over the pairs: GLS-2 {max_bits2}, GLV-GLS-4 {max_bits4}\n");
    println!("| arm | ns per call, median of {rounds} rounds | min of {rounds} rounds | ratio wnaf/arm (median) | ratio wnaf/arm (min) | ratio gls2/arm (median) |");
    println!("|:--|--:|--:|--:|--:|--:|");
    for i in 0..names.len() {
        let stage = i >= first_stage;
        println!(
            "| {} | {:.0} | {:.0} | {} | {} | {} |",
            names[i],
            med[i],
            min[i],
            if stage {
                "—".into()
            } else {
                format!("{:.3}", med[0] / med[i])
            },
            if stage {
                "—".into()
            } else {
                format!("{:.3}", min[0] / min[i])
            },
            if stage {
                "—".into()
            } else {
                format!("{:.3}", med[1] / med[i])
            },
        );
    }
    println!("\nall correctness checks passed: yes (sink {sink:x})");
    if let Some(path) = json_out {
        let doc = Json::obj(vec![
            ("schema", Json::str("gls_loop_bench/v1")),
            ("command", Json::str(args.join(" "))),
            ("d", Json::Num(d)),
            (
                "loop",
                Json::obj(vec![
                    ("A", Json::Num(lp.a)),
                    ("B", Json::Num(lp.b)),
                    ("norm", Json::str(lp.norm)),
                    ("cost_M", Json::Float(lp.cost)),
                ]),
            ),
            ("p", Json::str(p)),
            ("t", Json::str(&inst.t)),
            ("v", Json::str(&inst.v)),
            ("h", Json::Num(inst.h as i64)),
            ("r", Json::str(r)),
            ("mu", Json::str(&mu)),
            ("lambda", Json::str(&lam)),
            ("w", Json::Num(w as i64)),
            ("pairs", Json::Num(pairs as i64)),
            ("rounds", Json::Num(rounds as i64)),
            ("seed", Json::str(seed)),
            ("max_bits_gls2", Json::Num(max_bits2 as i64)),
            ("max_bits_glvgls4", Json::Num(max_bits4 as i64)),
            (
                "arms",
                Json::Arr(
                    (0..names.len())
                        .map(|i| {
                            Json::obj(vec![
                                ("name", Json::str(names[i])),
                                ("median_ns", Json::Float(med[i])),
                                ("min_ns", Json::Float(min[i])),
                                (
                                    "rounds_ns",
                                    Json::Arr(times[i].iter().map(|&x| Json::Float(x)).collect()),
                                ),
                            ])
                        })
                        .collect(),
                ),
            ),
            ("correctness_all_pass", Json::Bool(ok)),
        ]);
        std::fs::write(&path, doc.dump()).expect("write");
        eprintln!("wrote {path}");
    }
}
