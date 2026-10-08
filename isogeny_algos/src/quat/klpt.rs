//! KLPT (Kohel–Lauter–Petit–Tignol 2014) for left ideals of O_0 in B_{p,inf}, p = 3 mod 4:
//! given I, an equivalent left ideal J = I xi of norm l^e (e ~ 7/2 log_l p).
//!  1. L ~ I of prime norm N: N = Nrd(delta)/N(I) prime for a short delta in I; L = I deltabar/N(I).
//!  2. gamma = a + b i + c j + d k in O_0 of norm N l^{e0}: random (c, d), Cornacchia on
//!     a^2 + b^2 = N l^{e0} - p(c^2 + d^2).
//!  3. mu0 = C j + D k with gamma mu0 in L: a linear condition mod N (two functionals cutting
//!     L / N O_0 out of O_0 / N O_0).
//!  4. Strong approximation: mu = lambda mu0 + N mu1 of norm l^{e1}: lambda from a square root mod
//!     N, (c, d) from a linear condition mod N, then a close vector in a 2-dimensional lattice so
//!     that (l^{e1} - p(x^2 + y^2))/N^2 is a prime = 1 mod 4, split by Cornacchia.
//!  J = L betabar / N with beta = gamma mu (norm N l^{e0+e1}).
use super::*;
use crate::field::Rng;

fn i64i(v: i64) -> Int {
    Int::from(v)
}

/// x^2 + y^2 = m for a prime m = 1 mod 4 (or m = 2, or m a perfect square), by Cornacchia.
pub fn cornacchia_sum_two_squares(m: &Int) -> Option<(Int, Int)> {
    if m.is_neg() {
        return None;
    }
    if m.is_square() {
        return Some((m.isqrt(), Int::zero()));
    }
    if *m == i64i(2) {
        return Some((Int::one(), Int::one()));
    }
    if m.mod_u64(4) != 1 || !m.is_probable_prime() {
        return None;
    }
    let r0 = Int::sqrt_mod_prime(&(m - &Int::one()), m)?;
    let (mut a, mut b) = (m.clone(), r0);
    let s = m.isqrt();
    while b > s {
        let r = a.modulo(&b);
        a = b;
        b = r;
    }
    let y2 = m - &(&b * &b);
    if y2.is_square() {
        Some((b, y2.isqrt()))
    } else {
        None
    }
}

/// Coordinates of x in the lattice basis (None if x is not in the lattice).
pub fn coords(lat: &Lattice, x: &Quat) -> Option<[Int; 4]> {
    let mut v: [Int; 4] = x.c.clone().map(|c| &c * &lat.den);
    let mut u: [Int; 4] = std::array::from_fn(|_| Int::zero());
    for col in 0..4 {
        let piv = &lat.m[col][col] * &x.d;
        let (q, r) = v[col].divrem_trunc(&piv);
        if !r.is_zero() {
            return None;
        }
        for t in col..4 {
            v[t] = &v[t] - &(&(&q * &lat.m[col][t]) * &x.d);
        }
        u[col] = q;
    }
    Some(u)
}

/// Nullspace of a 4x4 matrix over F_n (n prime), rows = vectors; returns a basis of
/// {f : M f = 0}.
fn nullspace_mod(mat: &[[Int; 4]], n: &Int) -> Vec<[Int; 4]> {
    let mut a: Vec<[Int; 4]> = mat.iter().map(|r| r.clone().map(|x| x.modulo(n))).collect();
    let mut pivcols = vec![];
    let mut row = 0;
    for col in 0..4 {
        let Some(pr) = (row..a.len()).find(|&r| !a[r][col].is_zero()) else { continue };
        a.swap(row, pr);
        let inv = a[row][col].inv_mod(n).unwrap();
        a[row] = a[row].clone().map(|x| (&x * &inv).modulo(n));
        for r in 0..a.len() {
            if r != row && !a[r][col].is_zero() {
                let f = a[r][col].clone();
                let pr = a[row].clone();
                for t in 0..4 {
                    a[r][t] = (&a[r][t] - &(&f * &pr[t])).modulo(n);
                }
            }
        }
        pivcols.push(col);
        row += 1;
        if row == a.len() {
            break;
        }
    }
    let free: Vec<usize> = (0..4).filter(|c| !pivcols.contains(c)).collect();
    free.iter()
        .map(|&fc| {
            let mut v: [Int; 4] = std::array::from_fn(|_| Int::zero());
            v[fc] = Int::one();
            for (r, &pc) in pivcols.iter().enumerate() {
                v[pc] = (-&a[r][fc]).modulo(n);
            }
            v
        })
        .collect()
}

pub struct Klpt {
    /// output ideal, N(J) = l^e
    pub j: Lattice,
    pub e: u32,
    /// J = I xi
    pub xi: Quat,
    /// norm of the intermediate prime-norm ideal
    pub n_prime: Int,
    pub e0: u32,
    pub e1: u32,
    pub attempts: usize,
}

/// Left ideal O_0 n + O_0 a.
pub fn ideal_from(alg: &Alg, o0: &Lattice, n: &Int, a: &Quat) -> Lattice {
    let mut gens = vec![];
    for b in o0.basis() {
        gens.push(b.scale(n, &Int::one()));
        gens.push(alg.mul(&b, a));
    }
    Lattice::from_gens(&gens)
}

pub fn klpt(alg: &Alg, i: &Lattice, ell: u64, rng: &mut Rng) -> Option<Klpt> {
    let o0 = alg.o0();
    let p = alg.p.clone();
    let l = Int::from(ell);
    let (ni, di) = i.norm(&o0);
    assert!(di == Int::one(), "integral ideal expected");
    // 1. prime-norm equivalent ideal. Candidates delta: short vectors of I (Fincke-Pohst, a
    // capped number), then random small combinations of the whole reduced basis. (Ideals of
    // small norm in O_0 can have a very short rank-2 sublattice of sums of two squares, which
    // never gives N = 3 mod 4 and can hold millions of vectors below the next useful norm.)
    let d2 = &i.den * &i.den;
    let rb = i.reduced_basis(alg);
    let minn = rb.iter().map(|v| alg.qf(v)).min().unwrap();
    let to_cand = |v: [Int; 4], nv: Int| -> Option<([Int; 4], Int)> {
        let (q, r) = nv.divrem_trunc(&(&d2 * &ni));
        if !r.is_zero() {
            return None;
        }
        if q.mod_u64(4) == 3 && q != l && q != p && q.is_probable_prime() {
            Some((v, q))
        } else {
            None
        }
    };
    let mut cands: Vec<([Int; 4], Int)> = i
        .short_vectors(alg, &(&minn * &Int::from(64i64)), 20_000)
        .into_iter()
        .filter_map(|(v, nv)| to_cand(v, nv))
        .collect();
    let mut bnd = 1u64;
    while cands.len() < 8 && bnd <= 64 {
        for _ in 0..500 {
            let c: [Int; 4] = std::array::from_fn(|_| Int::from(rng.below(2 * bnd + 1) as i64 - bnd as i64));
            let v: [Int; 4] = std::array::from_fn(|t| {
                let mut s = Int::zero();
                for r in 0..4 {
                    s = &s + &(&c[r] * &rb[r][t]);
                }
                s
            });
            if v.iter().all(|x| x.is_zero()) {
                continue;
            }
            let nv = alg.qf(&v);
            if let Some(cd) = to_cand(v, nv) {
                if !cands.iter().any(|x| x.1 == cd.1) {
                    cands.push(cd);
                }
            }
        }
        bnd *= 2;
    }
    cands.sort_by(|a, b| a.1.cmp(&b.1));
    let mut attempts = 0usize;
    for (v, n) in cands.into_iter().take(24) {
        attempts += 1;
        let delta = Quat::new(v, i.den.clone());
        let xi1 = delta.conj().scale(&Int::one(), &ni);
        let lid = i.rmul(alg, &xi1);
        debug_assert_eq!(lid.norm(&o0), (n.clone(), Int::one()));
        if let Some(mut res) = klpt_prime(alg, &o0, &lid, &n, ell, rng) {
            res.xi = alg.mul(&xi1, &res.xi);
            res.attempts = attempts;
            return Some(res);
        }
    }
    None
}

/// Steps 2-4 for an ideal L of prime norm n.
fn klpt_prime(alg: &Alg, o0: &Lattice, lid: &Lattice, n: &Int, ell: u64, rng: &mut Rng) -> Option<Klpt> {
    let p = alg.p.clone();
    let l = Int::from(ell);
    let one = Int::one();
    // functionals cutting L/N O_0 out of O_0/N O_0
    let lrows: Vec<[Int; 4]> = lid.basis().iter().map(|b| coords(o0, b).unwrap()).collect();
    let funcs = nullspace_mod(&lrows, n); // f with sum_t row_t f_t = 0 for each row
    if funcs.len() != 2 {
        return None;
    }
    let eval_f = |f: &[Int; 4], x: &Quat| -> Int {
        let u = coords(o0, x).unwrap();
        let mut s = Int::zero();
        for t in 0..4 {
            s = &s + &(&u[t] * &f[t]);
        }
        s.modulo(n)
    };
    // 2. gamma of norm n l^{e0}
    let mut e0 = 0u32;
    while &(n * &l.pow(e0)) < &(&p * &Int::from(1i64 << 12)) {
        e0 += 1;
    }
    for _try in 0..200 {
        let target = n * &l.pow(e0);
        let bound = target.div_floor(&p).isqrt();
        let bnd = bound.to_i128().unwrap_or(i128::MAX).min(1 << 40) as u64;
        let mut gamma = None;
        for _ in 0..2000 {
            let c = Int::from(rng.below(bnd + 1));
            let d = Int::from(rng.below(bnd + 1));
            let m = &target - &(&p * &(&(&c * &c) + &(&d * &d)));
            if m.is_neg() || m.is_zero() {
                continue;
            }
            if let Some((a, b)) = cornacchia_sum_two_squares(&m) {
                gamma = Some(Quat::new([a, b, c, d], one.clone()));
                break;
            }
        }
        let Some(gamma) = gamma else {
            e0 += 1;
            continue;
        };
        // 3. mu0 = C j + D k with gamma mu0 in L
        let gj = alg.mul(&gamma, &alg.j());
        let gk = alg.mul(&gamma, &alg.k());
        let (u1, v1) = (eval_f(&funcs[0], &gj), eval_f(&funcs[0], &gk));
        let (u2, v2) = (eval_f(&funcs[1], &gj), eval_f(&funcs[1], &gk));
        // C u + D v = 0 for both functionals
        let (cc, dd) = if !u1.is_zero() || !v1.is_zero() { (v1.clone(), (-&u1).modulo(n)) } else { (v2.clone(), (-&u2).modulo(n)) };
        if cc.is_zero() && dd.is_zero() {
            continue;
        }
        if !(&(&cc * &u2) + &(&dd * &v2)).modulo(n).is_zero() || !(&(&cc * &u1) + &(&dd * &v1)).modulo(n).is_zero() {
            continue; // gamma in L or degenerate
        }
        // 4. strong approximation
        let nrm0 = (&p * &(&(&cc * &cc) + &(&dd * &dd))).modulo(n);
        if nrm0.is_zero() {
            continue;
        }
        let n4 = n.pow(4);
        let mut e1 = 0u32;
        while l.pow(e1) < &(&p * &n4) * &Int::from(1i64 << 20) {
            e1 += 1;
        }
        for _par in 0..2 {
            let le1 = l.pow(e1);
            let ratio = (&le1 * &nrm0.inv_mod(n).unwrap()).modulo(n);
            let Some(lambda) = Int::sqrt_mod_prime(&ratio, n) else {
                e1 += 1;
                continue;
            };
            if lambda.is_zero() {
                break;
            }
            if let Some(mu) = strong_approx(alg, n, &lambda, &cc, &dd, &le1, rng) {
                let beta = alg.mul(&gamma, &mu);
                debug_assert!(lid.contains(&beta));
                if !lid.contains(&beta) {
                    return None;
                }
                let jid = lid.rmul(alg, &beta.conj().scale(&one, n));
                let mut xi = beta.conj().scale(&one, n);
                let mut jid = jid;
                let mut e = e0 + e1;
                // remove l O_0 factors (non-primitive output)
                let lo = o0.scale(&l, &one);
                while e >= 2 && lo.contains_lattice(&jid) {
                    jid = jid.scale(&one, &l);
                    xi = xi.scale(&one, &l);
                    e -= 2;
                }
                return Some(Klpt { j: jid, e, xi, n_prime: n.clone(), e0, e1, attempts: 0 });
            }
            e1 += 1;
        }
    }
    None
}

/// mu = N a + N b i + x j + y k with x = lambda C + N c, y = lambda D + N d and
/// Nrd(mu) = target (target = l^{e1}).
fn strong_approx(alg: &Alg, n: &Int, lambda: &Int, cc: &Int, dd: &Int, target: &Int, rng: &mut Rng) -> Option<Quat> {
    let p = &alg.p;
    let one = Int::one();
    let two = Int::from(2i64);
    // c C + d D = t mod n, t = ((target - p lambda^2 (C^2 + D^2))/n) (2 p lambda)^{-1}
    let base = &(lambda * lambda) * &(&(cc * cc) + &(dd * dd));
    let num = target - &(p * &base);
    let (q, r) = num.divrem_trunc(n);
    if !r.is_zero() {
        return None;
    }
    let inv = (&(&two * p) * lambda).inv_mod(n)?;
    let t = (&q * &inv).modulo(n);
    // particular (c0, d0) and the lattice {(c, d): cC + dD = 0 mod n}
    let (c0, d0, b1, b2) = if !cc.modulo(n).is_zero() {
        let ci = cc.inv_mod(n).unwrap();
        ((&t * &ci).modulo(n), Int::zero(), [n.clone(), Int::zero()], [(-&(dd * &ci)).modulo(n), one.clone()])
    } else {
        let di = dd.inv_mod(n).unwrap();
        (Int::zero(), (&t * &di).modulo(n), [Int::zero(), n.clone()], [one.clone(), (-&(cc * &di)).modulo(n)])
    };
    // (x, y) = (lambda C + n c0, lambda D + n d0) + n (s b1 + u b2); lattice basis scaled by n
    let x0 = [&(lambda * cc) + &(n * &c0), &(lambda * dd) + &(n * &d0)];
    let mut v1 = [n * &b1[0], n * &b1[1]];
    let mut v2 = [n * &b2[0], n * &b2[1]];
    // Lagrange-Gauss reduction (Euclidean form)
    let dot = |a: &[Int; 2], b: &[Int; 2]| &(&a[0] * &b[0]) + &(&a[1] * &b[1]);
    loop {
        if dot(&v1, &v1) > dot(&v2, &v2) {
            std::mem::swap(&mut v1, &mut v2);
        }
        let mu = dot(&v1, &v2).div_round(&dot(&v1, &v1));
        if mu.is_zero() {
            break;
        }
        v2 = [&v2[0] - &(&mu * &v1[0]), &v2[1] - &(&mu * &v1[1])];
        if dot(&v2, &v2) >= dot(&v1, &v1) {
            break;
        }
    }
    // Babai: coefficients of -x0 in (v1, v2) by Cramer's rule, rounded
    let det = &(&v1[0] * &v2[1]) - &(&v1[1] * &v2[0]);
    let tx = -&x0[0];
    let ty = -&x0[1];
    let s0 = (&(&tx * &v2[1]) - &(&ty * &v2[0])).div_round(&det);
    let u0 = (&(&v1[0] * &ty) - &(&v1[1] * &tx)).div_round(&det);
    let n2 = n * n;
    // search a window of lattice points around the Babai point
    let mut tries = 0;
    for rad in 0i64..40 {
        for ds in -rad..=rad {
            for du in -rad..=rad {
                if ds.abs().max(du.abs()) != rad {
                    continue;
                }
                tries += 1;
                let s = &s0 + &Int::from(ds);
                let u = &u0 + &Int::from(du);
                let x = &(&x0[0] + &(&s * &v1[0])) + &(&u * &v2[0]);
                let y = &(&x0[1] + &(&s * &v1[1])) + &(&u * &v2[1]);
                let rest = target - &(p * &(&(&x * &x) + &(&y * &y)));
                if rest.is_neg() {
                    continue;
                }
                let (m, r) = rest.divrem_trunc(&n2);
                if !r.is_zero() {
                    continue;
                }
                if let Some((a, b)) = cornacchia_sum_two_squares(&m) {
                    let mu = Quat::new([n * &a, n * &b, x, y], one.clone());
                    debug_assert_eq!(alg.nrd(&mu), (target.clone(), one.clone()));
                    return Some(mu);
                }
            }
        }
        if tries > 4000 {
            break;
        }
    }
    let _ = rng;
    None
}
