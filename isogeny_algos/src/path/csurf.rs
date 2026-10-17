//! CSURF, "CSIDH on the surface" (Castryck–Decru 2019), as used in Bernstein–De Feo–Leroux–Smith
//! 2020 ("Faster computation of isogenies of large prime degree", §6), with √élu for the large odd
//! degrees.
//!
//! p = 8 * l_1 * ... * l_n - 1 is prime (so p = 7 mod 8; l_1 = 3 may appear squared, as in
//! CSURF-512). The curves are the F_p-curves on the surface of the 2-volcano,
//! End_{F_p}(E) = Z[(1 + pi)/2], written in Montgomery^- form E_A^-: y^2 = x^3 + A x^2 - x. Every
//! such curve has a unique A (see `from_plus`), its full 2-torsion {(0,0), (alpha,0), (-1/alpha,0)}
//! (alpha^2 + A alpha - 1 = 0) is rational, and the quadratic twist of E_A^- is E_{-A}^-
//! ((x, y) -> (-x, sqrt(-1) y); -1 is a non-square). The ideal class of l acts as in CSIDH, and
//! twisting inverts it: [a] E^t = ([a]^{-1} E)^t.
//!
//! # The ideal above 2 (derived here; cross-checked in the tests)
//! 2 splits in Q(sqrt(-p)) (-p = 1 mod 8); l_2 = (2, (pi - 1)/2). A rational 2-torsion point T
//! lies in ker((pi - 1)/2) iff T is in 2 E(F_p), i.e. for y^2 = (x - e_1)(x - e_2)(x - e_3) and
//! T = (e_1, 0), iff e_1 - e_2 and e_1 - e_3 are both squares. For (0,0) the product of the two
//! differences is -1 (a non-square): the 2-isogeny with kernel (0,0) is the vertical one (to the
//! floor). For T = (alpha, 0) the differences are alpha and alpha + 1/alpha; for (-1/alpha, 0) they
//! are -1/alpha and -(alpha + 1/alpha). So exactly one of the two is halvable over F_p: the
//! l_2-kernel is (alpha, 0) with alpha and alpha + 1/alpha squares (`alpha_of`), and the other is
//! the l_2-bar kernel (halvable on the twist).
//!
//! The step E_A^- -> E_A^- / <(alpha, 0)>: translating x -> x + alpha gives y^2 = x(x^2 + a x + b)
//! with a = 3 alpha + A = 2 alpha + 1/alpha and b = alpha^2 + 1, and the 2-isogeny with kernel (0,0)
//! (Silverman III.4.5) lands on Y^2 = X(X^2 - 2 a X + a^2 - 4 b), where a^2 - 4b = 1/alpha^2 after
//! using alpha^2 + A alpha = 1. Its non-zero 2-torsion abscissae are X_{1,2} = a ± 2 r with
//! r = sqrt(alpha^2 + 1) (a square root of alpha (alpha + 1/alpha), a square on the l_2 side: this is
//! the radical). Moving X_1 to 0 gives X(X^2 + (2 X_1 - X_2) X + C), C = X_1 (X_1 - X_2) = ±4 r X_1;
//! of the two roots exactly one has -C a square (the product of the two C's is -16 r^2 / alpha^2),
//! and with u = sqrt(-C) chosen to be a square, x = X / u is an F_p-isomorphism onto E_{A'}^-,
//!     A' = (2 X_1 - X_2) / u.
//! The dual's kernel is the image of (0,0) and (-1/alpha,0), i.e. X = 0, at x = -X_1/u; so the
//! next l_2-kernel is the other root, alpha' = u / X_1, and a chain of l_2-steps needs two square
//! roots per step (r and u) and no root-finding of alpha after the first. l_2-bar steps are l_2
//! steps on the twist: A -> -step(-A).
//!
//! # Odd l
//! E_A^- is F_p-isomorphic to a Montgomery^+ curve E^+_M: y^2 = x(x^2 + M x + 1) by moving
//! (alpha, 0) to 0 and scaling by s = sqrt(alpha^2 + 1) chosen to be a square: M = (3 alpha + A)/s
//! (`to_plus`). The odd part of the action runs there exactly as CSIDH does (direction by the
//! quadratic character of x^3 + M x^2 + x, kernel from the cofactor ladder, codomain by x-only Vélu
//! or by √élu, `kernel::sqrt_velu_mont`), and the result is brought back to its unique
//! Montgomery^- coefficient (`from_plus`).
use crate::bigint::Big;
use crate::field::{is_prime, Field, Rng, Zp};
use crate::fpm::FpM;
use crate::kernel::montgomery::*;
use crate::kernel::sqrt_velu_mont::SqrtVeluMont;

pub struct Csurf<F: Field> {
    pub fp: F,
    /// odd primes l_i dividing p + 1 (exponent index i + 1; index 0 is l = 2)
    pub primes: Vec<u64>,
    /// p + 1
    pub p1: Big,
    /// (p - 1) / 2
    half: Big,
    /// odd degrees from this one on use √élu instead of x-only Vélu
    pub sqrt_velu_min: u64,
}

/// Default threshold for √élu in the odd steps (the crossover measured for the Montgomery √élu
/// against projective Vélu in this crate is above the CSURF-512 primes; see the tests).
pub const SQRT_VELU_MIN: u64 = u64::MAX;

impl Csurf<Zp> {
    /// p = 8 prod(primes) - 1 (primes odd, with repetition allowed) if it is a prime below 2^62.
    pub fn new(primes: &[u64]) -> Option<Self> {
        let prod: u128 = primes.iter().map(|&l| l as u128).product();
        let p = 8 * prod - 1;
        if p >= (1u128 << 62) || !is_prime(p as u64) {
            return None;
        }
        let mut odd = primes.to_vec();
        odd.dedup();
        Some(Csurf {
            fp: Zp::new(p as u64),
            primes: odd,
            p1: Big::from_u128(p + 1),
            half: Big::from_u128((p - 1) / 2),
            sqrt_velu_min: SQRT_VELU_MIN,
        })
    }

    /// 3, 5, 7, ... (n - 1 consecutive odd primes) and the smallest larger prime making p prime.
    pub fn with_n_primes(n: usize) -> Option<Self> {
        let odd: Vec<u64> = (3u64..).filter(|&q| is_prime(q)).take(n + 40).collect();
        let mut base = odd[..n - 1].to_vec();
        for &last in &odd[n - 1..] {
            base.push(last);
            if let Some(c) = Csurf::new(&base) {
                return Some(c);
            }
            base.pop();
        }
        None
    }

    pub fn p(&self) -> u64 {
        self.fp.p
    }
}

impl<const N: usize> Csurf<FpM<N>> {
    /// p = 8 prod(primes) - 1 over N limbs (primes ascending, repetition allowed); None unless p
    /// passes the Miller–Rabin test.
    pub fn new_big(primes: &[u64]) -> Option<Self> {
        let mut prod = Big::from_u64(8);
        for &l in primes {
            prod = prod.mul_small(l);
        }
        let p = prod.sub_small(1);
        if p.bits() > 64 * N || p.bits() <= 64 * (N - 1) || !crate::bigint::is_probable_prime(&p) {
            return None;
        }
        let mut odd = primes.to_vec();
        odd.dedup();
        Some(Csurf {
            fp: FpM::<N>::new(&p),
            primes: odd,
            half: p.sub_small(1).shr(1),
            p1: prod,
            sqrt_velu_min: SQRT_VELU_MIN,
        })
    }
}

/// The CSURF-512 prime as printed in BDFLS 2020, §6:
/// p = 8 * 9 * (5 * 7 * ... * 337) * 349 * 353 * (367 * ... * 389) - 1 (66 consecutive primes
/// 5..337, then 349, 353, 367, 373, 379, 383, 389). This reading gives a 513-bit prime p = 7 mod 8
/// (checked by Miller–Rabin in the tests). Caveat: the paper's brace says "6 consecutive primes"
/// under (367 ... 389), which holds only five primes; the alternatives with 359 or 397 added
/// are not prime (tests), so this is the only consistent reading found, but it was not checked
/// against Castryck–Decru's own source.
pub fn csurf512_primes() -> Vec<u64> {
    let mut v = vec![3u64, 3];
    v.extend((5..=337u64).filter(|&q| is_prime(q)));
    v.extend([349, 353, 367, 373, 379, 383, 389]);
    v
}

impl Csurf<FpM<9>> {
    pub fn csurf512() -> Self {
        Csurf::<FpM<9>>::new_big(&csurf512_primes()).expect("CSURF-512 prime")
    }
}

/// j-invariant of E_A^-: j = 256 (A^2 + 3)^3 / (A^2 + 4).
pub fn j_minus<F: Field>(f: &F, a: F::E) -> F::E {
    let a2 = f.sq(a);
    let n = f.add(a2, f.from_u64(3));
    f.div(
        f.mul(f.from_u64(256), f.mul(n, f.sq(n))),
        f.add(a2, f.from_u64(4)),
    )
}

/// Short Weierstrass model of E_A^- (x_W = x + A/3): a = -1 - A^2/3, b = 2A^3/27 + A/3.
pub fn minus_to_weierstrass<F: Field>(f: &F, a: F::E) -> crate::curve::Curve<F::E> {
    let three = f.from_u64(3);
    let a2 = f.sq(a);
    let wa = f.sub(f.neg(f.one()), f.div(a2, three));
    let wb = f.add(
        f.div(f.mul(f.from_u64(2), f.mul(a2, a)), f.from_u64(27)),
        f.div(a, three),
    );
    crate::curve::Curve::new(wa, wb)
}

impl<F: Field> Csurf<F> {
    pub fn is_square(&self, v: F::E) -> bool {
        let f = &self.fp;
        f.is_zero(v) || f.pow_big(v, &self.half) == f.one()
    }

    /// The square root of a square v that is itself a square (exists: -1 is a non-square).
    fn sqrt_sq(&self, v: F::E) -> F::E {
        let f = &self.fp;
        let r = f.sqrt(v).expect("csurf: square root of a non-square");
        if self.is_square(r) {
            r
        } else {
            f.neg(r)
        }
    }

    /// The l_2 kernel abscissa of E_A^-: the root alpha of x^2 + A x - 1 with alpha and
    /// alpha + 1/alpha squares. Panics if E_A^- is not on the surface (A^2 + 4 a non-square).
    pub fn alpha_of(&self, a: F::E) -> F::E {
        let f = &self.fp;
        let d = f
            .sqrt(f.add(f.sq(a), f.from_u64(4)))
            .expect("csurf: A^2 + 4 is not a square (curve not on the surface)");
        let half = f.inv(f.from_u64(2));
        for r in [f.mul(f.sub(d, a), half), f.mul(f.sub(f.neg(d), a), half)] {
            if self.is_square(r) && self.is_square(f.add(r, f.inv(r))) {
                return r;
            }
        }
        panic!("csurf: no halvable 2-torsion point (p = 7 mod 8 surface curves always have one)");
    }

    /// One l_2-step from E_A^- with its l_2 kernel abscissa alpha: (A', alpha') where alpha' is the
    /// l_2 kernel abscissa of E_{A'}^- (module documentation).
    pub fn two_step(&self, alpha: F::E) -> (F::E, F::E) {
        let f = &self.fp;
        let ai = f.inv(alpha);
        let aa = f.add(f.add(alpha, alpha), ai); // 2 alpha + 1/alpha
        let r = f
            .sqrt(f.add(f.sq(alpha), f.one()))
            .expect("csurf: alpha^2 + 1 not a square (alpha is not the l_2 kernel)");
        let r2 = f.add(r, r);
        let (mut x1, mut x2) = (f.add(aa, r2), f.sub(aa, r2));
        let mut c = f.mul(x1, f.sub(x1, x2));
        if !self.is_square(f.neg(c)) {
            std::mem::swap(&mut x1, &mut x2);
            c = f.mul(x1, f.sub(x1, x2));
        }
        let u = self.sqrt_sq(f.neg(c));
        let ui = f.inv(u);
        let a_new = f.mul(f.sub(f.add(x1, x1), x2), ui);
        let alpha_new = f.mul(u, f.inv(x1));
        (a_new, alpha_new)
    }

    /// l_2^e for e > 0, l_2-bar^{|e|} for e < 0.
    pub fn two_power(&self, a: F::E, e: i32) -> F::E {
        let f = &self.fp;
        if e == 0 {
            return a;
        }
        // l_2-bar = l_2 on the twist E_{-A}
        let mut cur = if e > 0 { a } else { f.neg(a) };
        let mut alpha = self.alpha_of(cur);
        for _ in 0..e.unsigned_abs() {
            (cur, alpha) = self.two_step(alpha);
        }
        if e > 0 {
            cur
        } else {
            f.neg(cur)
        }
    }

    /// Montgomery^+ coefficient M of a curve F_p-isomorphic to E_A^- (module documentation).
    pub fn to_plus(&self, a: F::E) -> F::E {
        let f = &self.fp;
        let alpha = self.alpha_of(a);
        let s = self.sqrt_sq(f.add(f.sq(alpha), f.one()));
        f.div(f.add(f.mul(f.from_u64(3), alpha), a), s)
    }

    /// The unique A with E_A^- F_p-isomorphic to E^+_M: y^2 = x(x^2 + M x + 1). The 2-torsion
    /// roots beta, 1/beta of x^2 + M x + 1; moving (beta, 0) to 0 gives
    /// x(x^2 + (2 beta - 1/beta) x + beta^2 - 1), and exactly one of beta, 1/beta has 1 - beta^2 a
    /// square (their product is -(1 - beta^2)^2 / beta^2); then x -> x / u with u = sqrt(1 - beta^2)
    /// a square.
    /// Uniqueness: on E_A^- the only 2-torsion point whose translate has -b a square is (0, 0)
    /// (b = -1; for (alpha, 0), (-1/alpha, 0) b is a square), and then u = 1.
    pub fn from_plus(&self, m: F::E) -> F::E {
        let f = &self.fp;
        let d = f
            .sqrt(f.sub(f.sq(m), f.from_u64(4)))
            .expect("csurf: M^2 - 4 not a square (curve not on the surface)");
        let mut beta = f.mul(f.sub(d, m), f.inv(f.from_u64(2)));
        let mut w = f.sub(f.one(), f.sq(beta));
        if !self.is_square(w) {
            beta = f.inv(beta);
            w = f.sub(f.one(), f.sq(beta));
        }
        let u = self.sqrt_sq(w);
        f.div(f.sub(f.add(beta, beta), f.inv(beta)), u)
    }

    fn rhs_plus(&self, m: F::E, x: F::E) -> F::E {
        let f = &self.fp;
        f.mul(f.add(f.mul(f.add(x, m), x), f.one()), x)
    }

    fn cofactor(&self, divide_by: &[u64]) -> Big {
        let mut k = self.p1.clone();
        for &l in divide_by {
            k = k.divrem_small(l).0;
        }
        k
    }

    /// One odd l-step on E^+_M (positive: kernel in E(F_p); negative: in the twist), with the
    /// kernel generator returned for checking: (M', kernel x-coordinate (X:Z)).
    pub fn odd_step_plus(
        &self,
        m: F::E,
        ell: u64,
        positive: bool,
        rng: &mut Rng,
    ) -> (F::E, XZ<F::E>) {
        let f = &self.fp;
        let e24 = a24(f, m);
        let cof = self.cofactor(&[ell]);
        loop {
            let x = f.random(rng);
            let r = self.rhs_plus(m, x);
            if f.is_zero(r) || self.is_square(r) != positive {
                continue;
            }
            let k = ladder_xz(f, e24, (x, f.one()), &cof);
            if f.is_zero(k.1) {
                continue;
            }
            return (self.codomain_plus(m, k, ell), k);
        }
    }

    /// Codomain coefficient of E^+_M / <K> for K of order l: √élu from `sqrt_velu_min` on,
    /// x-only Vélu below.
    pub fn codomain_plus(&self, m: F::E, k: XZ<F::E>, ell: u64) -> F::E {
        let f = &self.fp;
        if ell >= self.sqrt_velu_min {
            if let Some(sv) = SqrtVeluMont::new(f, m, k, ell) {
                return sv.codomain(f);
            }
        }
        let ms = multiples(f, a24(f, m), k, ((ell - 1) / 2) as usize);
        velu_codomain(f, m, &ms, ell)
    }

    /// [l_2^{e_0} prod l_i^{e_i}] E_A^- one isogeny at a time: odd steps on the Montgomery^+
    /// model, then the 2-steps radically. `exps[0]` is the exponent of l_2, `exps[i]` that of
    /// `primes[i - 1]`.
    pub fn action(&self, a: F::E, exps: &[i32], rng: &mut Rng) -> F::E {
        assert_eq!(exps.len(), self.primes.len() + 1);
        let mut a = a;
        if exps[1..].iter().any(|&e| e != 0) {
            let mut m = self.to_plus(a);
            for (i, &e) in exps[1..].iter().enumerate() {
                for _ in 0..e.unsigned_abs() {
                    m = self.odd_step_plus(m, self.primes[i], e > 0, rng).0;
                }
            }
            a = self.from_plus(m);
        }
        self.two_power(a, exps[0])
    }

    /// As [`Csurf::action`], with the odd part batched as in CLMPR (one sampled point serves every
    /// prime of its sign: clear the cofactor once, peel off one prime at a time, push the point
    /// through each isogeny; √élu for the push and the codomain from `sqrt_velu_min` on).
    pub fn action_batched(&self, a: F::E, exps: &[i32], rng: &mut Rng) -> F::E {
        assert_eq!(exps.len(), self.primes.len() + 1);
        let f = &self.fp;
        let mut e = exps[1..].to_vec();
        let mut a = a;
        if e.iter().any(|&x| x != 0) {
            let mut m = self.to_plus(a);
            while e.iter().any(|&x| x != 0) {
                let x = f.random(rng);
                let r = self.rhs_plus(m, x);
                if f.is_zero(r) {
                    continue;
                }
                let positive = self.is_square(r);
                let set: Vec<usize> = (0..e.len())
                    .filter(|&i| e[i] != 0 && (e[i] > 0) == positive)
                    .collect();
                if set.is_empty() {
                    continue;
                }
                let set_primes: Vec<u64> = set.iter().map(|&i| self.primes[i]).collect();
                let mut q = ladder_xz(f, a24(f, m), (x, f.one()), &self.cofactor(&set_primes));
                let mut remaining = set_primes.clone();
                for &i in set.iter().rev() {
                    let ell = self.primes[i];
                    remaining.retain(|&l| l != ell);
                    if f.is_zero(q.1) {
                        break;
                    }
                    let mut k = Big::from_u64(1);
                    for &l in &remaining {
                        k = k.mul_small(l);
                    }
                    let kp = ladder_xz(f, a24(f, m), q, &k);
                    if f.is_zero(kp.1) {
                        continue;
                    }
                    let sv = if ell >= self.sqrt_velu_min {
                        SqrtVeluMont::new(f, m, kp, ell)
                    } else {
                        None
                    };
                    if let Some(sv) = sv {
                        q = sv.eval_xz(f, q);
                        m = sv.codomain(f);
                    } else {
                        let ms = multiples(f, a24(f, m), kp, ((ell - 1) / 2) as usize);
                        let m_new = velu_codomain(f, m, &ms, ell);
                        q = isog_xz(f, &ms, q);
                        m = m_new;
                    }
                    e[i] -= if positive { 1 } else { -1 };
                }
            }
            a = self.from_plus(m);
        }
        self.two_power(a, exps[0])
    }
}
