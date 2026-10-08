//! Deuring correspondence on E_0: y^2 = x^3 + x (p = 3 mod 4), End(E_0) = O_0 with
//! i -> (x, y) |-> (-x, sqrt(-1) y) and j -> the p-power Frobenius (so k = ij acts as i o pi).
//! Points live over F_{p^4}, where E_0 has full (p^2 - 1)-torsion; the curves and isogenies are
//! defined over F_{p^2} (j-invariants are projected back and checked).
//!
//! Ideal -> curve: an equivalent ideal J = I deltabar/N(I) whose norm divides the odd part T of
//! p^2 - 1 (all its prime powers have rational torsion over F_{p^4}), the kernel
//! E_0[J] = {P in E_0[N(J)] : b(P) = 0 for b in a basis of J} from discrete logarithms in a
//! torsion basis (2-dimensional Pohlig–Hellman on the l^e parts), then a chain of l-isogenies
//! (Vélu). Kernel -> ideal: I_K = O_0 n + O_0 alpha with alpha(K) = 0 solved as two linear
//! congruences mod n.
use super::*;
use crate::curve::{jinv, padd, pmul, Curve, Pt};
use crate::ext::Fp4;
use crate::field::{factor_u64, Field, Rng};
use crate::kernel::velu::velu_cyclic;
use std::collections::HashMap;

type E2 = ((u64, u64), (u64, u64));
type Zp2 = Fp4;

pub struct TorsionBasis {
    pub ell: u64,
    pub e: u32,
    pub p: Pt<E2>,
    pub q: Pt<E2>,
}

pub struct Deuring {
    pub p: u64,
    /// F_{p^4}
    pub f2: Fp4,
    pub alg: Alg,
    pub o0: Lattice,
    pub e0: Curve<E2>,
    /// odd prime powers dividing p^2 - 1, with bases of E_0[l^e] over F_{p^4}
    pub tors: Vec<TorsionBasis>,
    pub t_odd: u64,
}

fn neg_pt(f: &Zp2, p: &Pt<E2>) -> Pt<E2> {
    match *p {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => Pt::Aff(x, f.neg(y)),
    }
}

fn smul(f: &Zp2, e: &Curve<E2>, p: &Pt<E2>, k: i128, n: u64) -> Pt<E2> {
    let k = k.rem_euclid(n as i128) as u128;
    pmul(f, e, p, k)
}

impl Deuring {
    pub fn new(p: u64, rng: &mut Rng) -> Deuring {
        assert!(p % 4 == 3);
        let f2 = Fp4::new(p);
        let alg = Alg::new(&Int::from(p));
        let o0 = alg.o0();
        let e0 = Curve::new(f2.one(), f2.zero());
        let mut tors = vec![];
        let mut t_odd = 1u64;
        let order = (p as u128 + 1) * (p as u128 - 1);
        let mut fac: Vec<(u64, u32)> = factor_u64(p + 1);
        for (l, e) in factor_u64(p - 1) {
            match fac.iter_mut().find(|x| x.0 == l) {
                Some(x) => x.1 += e,
                None => fac.push((l, e)),
            }
        }
        fac.sort();
        for (ell, e) in fac {
            // 2-power torsion is not used (O_0 has denominators 2); primes above 256 would make
            // the l^2-entry discrete-log tables large
            if ell == 2 || ell > 256 {
                continue;
            }
            let le = ell.pow(e);
            t_odd *= le;
            let cof = order / le as u128;
            let rand_tors = |rng: &mut Rng| loop {
                let r = crate::curve::random_point_f(&f2, &e0, rng);
                let t = pmul(&f2, &e0, &r, cof);
                if pmul(&f2, &e0, &t, (le / ell) as u128) != Pt::Inf {
                    return t;
                }
            };
            let bp = rand_tors(rng);
            let bp_small = pmul(&f2, &e0, &bp, (le / ell) as u128);
            let bq = loop {
                let q = rand_tors(rng);
                let qs = pmul(&f2, &e0, &q, (le / ell) as u128);
                // independent mod l: qs not a multiple of bp_small
                let mut m = Pt::Inf;
                let mut dep = false;
                for _ in 0..ell {
                    if m == qs {
                        dep = true;
                        break;
                    }
                    m = padd(&f2, &e0, &m, &bp_small);
                }
                if !dep {
                    break q;
                }
            };
            tors.push(TorsionBasis { ell, e, p: bp, q: bq });
        }
        Deuring { p, f2, alg, o0, e0, tors, t_odd }
    }

    fn i_map(&self, p: &Pt<E2>) -> Pt<E2> {
        match *p {
            Pt::Inf => Pt::Inf,
            Pt::Aff(x, y) => Pt::Aff(self.f2.neg(x), self.f2.mul(self.f2.embed((0, 1)), y)),
        }
    }
    fn frob(&self, p: &Pt<E2>) -> Pt<E2> {
        match *p {
            Pt::Inf => Pt::Inf,
            Pt::Aff(x, y) => Pt::Aff(self.f2.frob(x), self.f2.frob(y)),
        }
    }

    /// alpha(P) for alpha in O_0 and P of odd order dividing n.
    pub fn act(&self, alpha: &Quat, p: &Pt<E2>, n: u64) -> Pt<E2> {
        let f = &self.f2;
        let e = &self.e0;
        let pi = self.frob(p);
        let imgs = [*p, self.i_map(p), pi, self.i_map(&pi)];
        let mut s = Pt::Inf;
        for t in 0..4 {
            let c = alpha.c[t].mod_u64(n) as i128;
            s = padd(f, e, &s, &smul(f, e, &imgs[t], c, n));
        }
        if alpha.d != Int::one() {
            let dinv = alpha.d.inv_mod(&Int::from(n)).expect("denominator invertible mod n");
            s = smul(f, e, &s, dinv.mod_u64(n) as i128, n);
        }
        s
    }

    /// (a, b) mod l^k with R = a P + b Q, for a basis (P, Q) of E[l^k] (Pohlig–Hellman, digit by
    /// digit, l^2 table lookups per digit).
    pub fn dlog2(&self, ell: u64, k: u32, p: &Pt<E2>, q: &Pt<E2>, r: &Pt<E2>) -> Option<(u64, u64)> {
        let f = &self.f2;
        let e = &self.e0;
        let lk = ell.pow(k);
        let top = ell.pow(k - 1);
        let (p1, q1) = (pmul(f, e, p, top as u128), pmul(f, e, q, top as u128));
        let mut table: HashMap<Pt<E2>, (u64, u64)> = HashMap::new();
        for s in 0..ell {
            for t in 0..ell {
                let v = padd(f, e, &pmul(f, e, &p1, s as u128), &pmul(f, e, &q1, t as u128));
                table.insert(v, (s, t));
            }
        }
        let (mut a, mut b) = (0u64, 0u64);
        let mut lpow = 1u64;
        for d in 0..k {
            // T = [l^{k-1-d}] (R - a P - b Q)
            let rem = padd(f, e, r, &neg_pt(f, &padd(f, e, &pmul(f, e, p, a as u128), &pmul(f, e, q, b as u128))));
            let t = pmul(f, e, &rem, ell.pow(k - 1 - d) as u128);
            let (s, u) = *table.get(&t)?;
            a += s * lpow;
            b += u * lpow;
            lpow *= ell;
        }
        Some((a % lk, b % lk))
    }

    /// Generators (one per prime power l^a || N(J)) of the kernel E_0[J], N(J) | t_odd.
    pub fn kernel_of_ideal(&self, j: &Lattice) -> Option<Vec<(u64, u32, Pt<E2>)>> {
        let n = j.norm_int(&self.o0).to_i128()? as u64;
        let mut out = vec![];
        for (ell, a) in factor_u64(n) {
            let tb = self.tors.iter().find(|t| t.ell == ell && t.e >= a)?;
            let la = ell.pow(a);
            let (pa, qa) = (
                pmul(&self.f2, &self.e0, &tb.p, ell.pow(tb.e - a) as u128),
                pmul(&self.f2, &self.e0, &tb.q, ell.pow(tb.e - a) as u128),
            );
            // rows: (A, B) with b(P) = A1 P + A2 Q, b(Q) = B1 P + B2 Q; condition on (x, y):
            // x b(P) + y b(Q) = 0  <=>  x A1 + y B1 = 0 and x A2 + y B2 = 0 (mod l^a)
            let mut rows: Vec<(u64, u64)> = vec![];
            for b in j.basis() {
                let (a1, a2) = self.dlog2(ell, a, &pa, &qa, &self.act(&b, &pa, la))?;
                let (b1, b2) = self.dlog2(ell, a, &pa, &qa, &self.act(&b, &qa, la))?;
                rows.push((a1, b1));
                rows.push((a2, b2));
            }
            // cyclic kernel of order l^a: (1, y) or (x, 1), coordinates lifted digit by digit
            let solve = |swap: bool| -> Option<(u64, u64)> {
                let mut cands = vec![0u64];
                let mut m = 1u64;
                for _ in 0..a {
                    let m2 = m * ell;
                    let mut next = vec![];
                    for &c in &cands {
                        for d in 0..ell {
                            let v = c + d * m;
                            let ok = rows.iter().all(|&(r0, r1)| {
                                let (x, y) = if swap { (v, 1) } else { (1, v) };
                                ((r0 as u128 * x as u128 + r1 as u128 * y as u128) % m2 as u128) == 0
                            });
                            if ok {
                                next.push(v);
                            }
                        }
                    }
                    cands = next;
                    if cands.is_empty() {
                        return None;
                    }
                    m = m2;
                }
                let v = cands[0];
                Some(if swap { (v, 1) } else { (1, v) })
            };
            let (x, y) = solve(false).or_else(|| solve(true))?;
            let k = padd(&self.f2, &self.e0, &pmul(&self.f2, &self.e0, &pa, x as u128), &pmul(&self.f2, &self.e0, &qa, y as u128));
            out.push((ell, a, k));
        }
        Some(out)
    }

    /// The codomain of the isogeny with kernel generated by the given points (orders l^a), by
    /// chains of l-isogenies (Vélu), pushing the remaining generators through each step.
    pub fn isogeny_from_kernels(&self, ks: &[(u64, u32, Pt<E2>)]) -> Curve<E2> {
        let f = &self.f2;
        let mut cur = self.e0;
        let mut pend: Vec<(u64, u32, Pt<E2>)> = ks.to_vec();
        while let Some((ell, a, k)) = pend.first().cloned() {
            if a == 0 {
                pend.remove(0);
                continue;
            }
            let k1 = pmul(f, &cur, &k, ell.pow(a - 1) as u128);
            let iso = velu_cyclic(f, &cur, &k1, ell);
            for item in pend.iter_mut() {
                item.2 = crate::curve::Isogeny::eval(&iso, f, &item.2);
            }
            pend[0].1 -= 1;
            cur = iso.cod;
        }
        cur
    }

    /// An equivalent ideal J = I deltabar / N(I) with N(J) | t_odd (None if not found).
    pub fn smooth_equivalent(&self, i: &Lattice, _rng: &mut Rng) -> Option<(Lattice, Quat)> {
        let alg = &self.alg;
        let (ni, _) = i.norm(&self.o0);
        let d2 = &i.den * &i.den;
        let t = Int::from(self.t_odd);
        let ok_norm = |nv: &Int| -> bool {
            let (q, r) = nv.divrem_trunc(&(&d2 * &ni));
            r.is_zero() && !q.is_zero() && t.modulo(&q).is_zero()
        };
        let rb = i.reduced_basis(alg);
        // every delta with N(J) | t_odd has Nrd(delta) <= t_odd N(I): grow the bound from the
        // minimum by factors of 4 up to that, stopping at the first admissible norm
        let minn = rb.iter().map(|v| alg.qf(v)).min().unwrap();
        let exact = &(&t * &ni) * &d2;
        let mut bound = &minn * &Int::from(4i64);
        let mut cands: Vec<([Int; 4], Int)> = vec![];
        loop {
            let b = if bound < exact { bound.clone() } else { exact.clone() };
            let sv = i.short_vectors(alg, &b, 200_000);
            cands = sv.into_iter().filter(|(_, nv)| ok_norm(nv)).collect();
            if !cands.is_empty() || b == exact {
                break;
            }
            bound = &bound * &Int::from(4i64);
        }
        let (v, _) = cands.into_iter().min_by(|a, b| a.1.cmp(&b.1))?;
        let delta = Quat::new(v, i.den.clone());
        let xi = delta.conj().scale(&Int::one(), &ni);
        Some((i.rmul(alg, &xi), xi))
    }

    /// Deuring: a curve whose endomorphism ring is the right order of I.
    pub fn ideal_to_curve(&self, i: &Lattice, rng: &mut Rng) -> Option<Curve<E2>> {
        let (j, _) = self.smooth_equivalent(i, rng)?;
        let ks = self.kernel_of_ideal(&j)?;
        Some(self.isogeny_from_kernels(&ks))
    }

    /// j-invariant (in F_{p^2}) of the Deuring image of I.
    pub fn ideal_to_j(&self, i: &Lattice, rng: &mut Rng) -> Option<(u64, u64)> {
        let j = jinv(&self.f2, &self.ideal_to_curve(i, rng)?);
        Some(self.f2.project(j).expect("codomain defined over F_{p^2}"))
    }

    /// The left O_0-ideal of a cyclic kernel <K> of order n | t_odd: O_0 n + O_0 alpha with
    /// alpha(K) = 0. Per prime power l^a: the images v_r = b_r(K) of the basis of O_0 in a
    /// torsion basis; two coordinates of alpha are drawn at random, the other two solved from the
    /// 2x2 system (any pair with invertible determinant); coordinates are combined by CRT.
    pub fn ideal_of_kernel(&self, k: &Pt<E2>, n: u64, rng: &mut Rng) -> Option<Lattice> {
        let basis = self.o0.basis();
        let fac = factor_u64(n);
        let mut images = vec![];
        for &(ell, a) in &fac {
            let la = ell.pow(a);
            let tb = self.tors.iter().find(|t| t.ell == ell && t.e >= a)?;
            let pa = pmul(&self.f2, &self.e0, &tb.p, ell.pow(tb.e - a) as u128);
            let qa = pmul(&self.f2, &self.e0, &tb.q, ell.pow(tb.e - a) as u128);
            let kk = pmul(&self.f2, &self.e0, k, (n / la) as u128);
            let mut v = vec![];
            for b in &basis {
                v.push(self.dlog2(ell, a, &pa, &qa, &self.act(b, &kk, la))?);
            }
            images.push((la, v));
        }
        let inv_mod = |x: i128, m: u64| Int::from(x.rem_euclid(m as i128) as i64).inv_mod(&Int::from(m)).map(|v| v.mod_u64(m) as i128);
        for _ in 0..200 {
            let mut coords = [0u128; 4];
            let mut modulus = 1u128;
            let mut ok = true;
            for (la, v) in &images {
                let m = *la as i128;
                // a pair (r, s) with invertible determinant
                let mut pair = None;
                for r in 0..4 {
                    for s2 in r + 1..4 {
                        let det = v[r].0 as i128 * v[s2].1 as i128 - v[s2].0 as i128 * v[r].1 as i128;
                        if let Some(di) = inv_mod(det, *la) {
                            pair = Some((r, s2, di));
                            break;
                        }
                    }
                    if pair.is_some() {
                        break;
                    }
                }
                let Some((r, s2, di)) = pair else {
                    ok = false;
                    break;
                };
                let mut c = [0i128; 4];
                let (mut rhs0, mut rhs1) = (0i128, 0i128);
                for t in 0..4 {
                    if t != r && t != s2 {
                        c[t] = rng.below(*la) as i128;
                        rhs0 -= c[t] * v[t].0 as i128;
                        rhs1 -= c[t] * v[t].1 as i128;
                    }
                }
                let (rhs0, rhs1) = (rhs0.rem_euclid(m), rhs1.rem_euclid(m));
                // [[v_r.0, v_s.0], [v_r.1, v_s.1]] (c_r, c_s) = (rhs0, rhs1)
                let (a11, a12, a21, a22) = (v[r].0 as i128, v[s2].0 as i128, v[r].1 as i128, v[s2].1 as i128);
                c[r] = ((a22 * rhs0 - a12 * rhs1).rem_euclid(m) * di).rem_euclid(m);
                c[s2] = ((a11 * rhs1 - a21 * rhs0).rem_euclid(m) * di).rem_euclid(m);
                // CRT coordinate-wise
                let minv = inv_mod(modulus as i128, *la).unwrap();
                for t in 0..4 {
                    let tt = ((c[t] - coords[t] as i128).rem_euclid(m) * minv).rem_euclid(m) as u128;
                    coords[t] += modulus * tt;
                }
                modulus *= *la as u128;
            }
            if !ok {
                return None;
            }
            let mut alpha = Quat::zero();
            for (t, b) in basis.iter().enumerate() {
                alpha = alpha.add(&b.scale(&Int::from(coords[t] as i128), &Int::one()));
            }
            let id = super::klpt::ideal_from(&self.alg, &self.o0, &Int::from(n), &alpha);
            if id.norm(&self.o0) == (Int::from(n), Int::one()) {
                return Some(id);
            }
        }
        None
    }
}
