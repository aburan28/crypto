//! Isogeny-loop endomorphisms on CM curves and their transport to a GLS twist.
//!
//! Koshelev–Sanso (ePrint 2024/1985): on a curve `E/F_p` with `End(E) = O_K` (CM discriminant
//! `D`), a closed walk of horizontal isogenies is an endomorphism `alpha` of degree
//! `N(alpha) = prod l_i`, evaluated step by step.  On the quadratic twist `E'/F_{p^2}`
//! (Galbraith–Lin–Scott) the same loop, twisted, is an `F_{p^2}`-rational endomorphism
//! `alpha'`, alongside the GLS map `psi` (`psi^2 = -1`).  This module builds everything needed to
//! check that on points:
//!
//! * [`cheapest_loop`]: the minimum-cost primitive `alpha = (A + B sqrt D)/2` (exhaustive up to
//!   the norm bound implied by the cost model `7.5 l` per `l`-step);
//! * [`find_instance`]: `p = (t^2 + |D| v^2)/4` prime, `p = 3 mod 4`, with
//!   `#E'(F_{p^2}) = (p-1)^2 + t^2 = h r`, `r` prime;
//! * [`build_curve`]: a curve with `j` a root of `H_D mod p` and trace `t`;
//! * [`closing_loops`]: every non-backtracking horizontal walk with the loop's degree sequence
//!   that returns to the start (they are `alpha` and its conjugate), with the final isomorphism;
//! * [`Chain::eval_twisted`]: the twisted chain on a point of `E'(F_{p^2})`.
use crate::api::{self, IsogenyRecord, PrimeCurve, PrimeFieldInt};
use crate::cm::hilbert_class_poly;
use crate::field::{Field, Rng};
use crate::int::Int;
use crate::poly;
use std::collections::HashSet;

/// Cost per bit of degree in the projective model: min_l 7.5 l / log2 l (l = 3).
pub fn kappa_projective() -> f64 {
    7.5 * 3.0 / 3f64.log2()
}

fn kronecker_odd(d: i64, l: u64) -> i32 {
    let dm = d.rem_euclid(l as i64) as u128;
    if dm == 0 {
        return 0;
    }
    let (mut r, mut b, mut e) = (1u128, dm, (l - 1) / 2);
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % l as u128;
        }
        b = b * b % l as u128;
        e >>= 1;
    }
    if r == 1 {
        1
    } else {
        -1
    }
}

fn kronecker(d: i64, l: u64) -> i32 {
    if l == 2 {
        match d.rem_euclid(8) {
            1 | 7 => 1,
            3 | 5 => -1,
            _ => 0,
        }
    } else {
        kronecker_odd(d, l)
    }
}

fn factor_u64(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = vec![];
    let mut q = 2u64;
    while q * q <= n {
        if n.is_multiple_of(q) {
            let mut e = 0;
            while n.is_multiple_of(q) {
                n /= q;
                e += 1;
            }
            out.push((q, e));
        }
        q += if q == 2 { 1 } else { 2 };
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

/// A loop endomorphism `alpha = (a + b sqrt d)/2`.
#[derive(Clone, Debug)]
pub struct LoopElement {
    pub d: i64,
    pub a: i64,
    pub b: i64,
    pub norm: u64,
    /// (l, number of steps)
    pub steps: Vec<(u64, u32)>,
    /// projective cost, 7.5 l per l-step (1 for a non-trivial automorphism)
    pub cost: f64,
}

impl LoopElement {
    /// The degree of each step, in walking order.
    pub fn degree_sequence(&self) -> Vec<u64> {
        let mut out = vec![];
        for &(l, e) in &self.steps {
            for _ in 0..e {
                out.push(l);
            }
        }
        out
    }
}

/// d_min of the maximal order of discriminant d.
pub fn d_min(d: i64) -> u64 {
    if d.rem_euclid(4) == 0 {
        (-d / 4) as u64
    } else {
        ((1 - d) / 4) as u64
    }
}

/// The cheapest primitive `alpha` with `B >= 1` and `N(alpha) <= nmax`: all ties, and the
/// number of elements enumerated.
pub fn cheapest_loops_upto(d: i64, nmax: u64) -> (Vec<LoopElement>, u64) {
    let ad = (-d) as u64;
    let even = d.rem_euclid(4) == 0;
    let mut best: Vec<LoopElement> = vec![];
    let mut best_cost = f64::INFINITY;
    let mut enumerated = 0u64;
    let mut b: i64 = 1;
    while (ad as u128) * (b as u128) * (b as u128) <= 4 * nmax as u128 {
        let rest = 4 * nmax as u128 - (ad as u128) * (b as u128) * (b as u128);
        let amax = (rest as f64).sqrt() as i64 + 1;
        for a in -amax..=amax {
            if even {
                if a.rem_euclid(2) != 0 {
                    continue;
                }
            } else if (a - b).rem_euclid(2) != 0 {
                continue;
            }
            let n4 = (a as i128) * (a as i128) + (ad as i128) * (b as i128) * (b as i128);
            if n4 % 4 != 0 {
                continue;
            }
            let n = (n4 / 4) as u64;
            if n == 0 || n > nmax {
                continue;
            }
            let (x, y) = if even { (a / 2, b) } else { ((a - b) / 2, b) };
            if gcd_i64(x, y) != 1 {
                continue;
            }
            enumerated += 1;
            let fac = factor_u64(n);
            if fac.iter().any(|&(l, _)| kronecker(d, l) == -1) {
                continue;
            }
            let cost = if n == 1 {
                1.0
            } else {
                fac.iter().map(|&(l, e)| 7.5 * l as f64 * e as f64).sum()
            };
            if cost < best_cost - 1e-9 {
                best_cost = cost;
                best.clear();
            }
            if (cost - best_cost).abs() < 1e-9 {
                best.push(LoopElement {
                    d,
                    a,
                    b,
                    norm: n,
                    steps: fac,
                    cost,
                });
            }
        }
        b += 1;
    }
    (best, enumerated)
}

fn gcd_i64(mut a: i64, mut b: i64) -> i64 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// The cheapest loop, certified: a first pass gives an upper bound C on the cost, and every
/// loop of cost <= C has norm <= 2^(C / kappa), which the second pass enumerates.
pub fn cheapest_loop(d: i64) -> (LoopElement, usize, u64, u64) {
    let dm = d_min(d).max(1);
    let (first, _) = cheapest_loops_upto(d, dm * 64);
    let ub = first.first().map(|l| l.cost).unwrap_or(f64::INFINITY);
    let nmax = if ub.is_finite() {
        (2f64.powf(ub / kappa_projective())).floor() as u64
    } else {
        dm * 4096
    };
    let nmax = nmax.max(dm * 64);
    let (best, enumerated) = cheapest_loops_upto(d, nmax);
    let ties = best.len();
    (
        best.into_iter().next().expect("a loop exists"),
        ties,
        enumerated,
        nmax,
    )
}

/// A GLS instance for discriminant `d`.
#[derive(Clone, Debug)]
pub struct Instance {
    pub d: i64,
    pub p: Int,
    pub t: Int,
    pub v: Int,
    /// #E'(F_{p^2}) = h r
    pub h: u64,
    pub r: Int,
    pub tries: u64,
}

fn random_below_bits(rng: &mut Rng, bits: usize) -> Int {
    let words = bits.div_ceil(64);
    let mut v = Int::zero();
    for _ in 0..words {
        v = &v.shl(64) + &Int::from(rng.next());
    }
    v.shr(words * 64 - bits)
}

/// Search `t, v` with `p = (t^2 + |d| v^2)/4` a prime of exactly `pbits` bits, `p = 3 mod 4`,
/// `gcd(v, l) = 1` for every `l` in `avoid` (so that only the horizontal `l`-isogenies are
/// rational), and `(p-1)^2 + t^2 = h r` with `r` prime and `h <= hmax`.
pub fn find_instance(
    d: i64,
    pbits: usize,
    avoid: &[u64],
    hmax: u64,
    seed: u64,
    max_tries: u64,
) -> Option<Instance> {
    let mut rng = Rng::new(seed);
    let ad = Int::from(-d);
    let four = Int::from(4i64);
    let lad = ((-d) as f64).log2();
    for tries in 1..=max_tries {
        let vbits = ((pbits as f64 + 2.0 - lad) / 2.0).floor() as usize;
        let tbits = (pbits + 2) / 2;
        let v = random_below_bits(&mut rng, vbits.max(1));
        let mut t = random_below_bits(&mut rng, tbits);
        if v.is_zero() || avoid.iter().any(|&l| v.mod_u64(l) == 0) {
            continue;
        }
        if d.rem_euclid(4) == 0 {
            if !t.is_even() {
                t = &t + &Int::one();
            }
        } else if (&t - &v).mod_u64(2) == 1 {
            t = &t + &Int::one();
        }
        let n4 = &(&t * &t) + &(&ad * &(&v * &v));
        if n4.mod_u64(4) != 0 {
            continue;
        }
        let p = n4.div_floor(&four);
        if p.bits() != pbits || p.mod_u64(4) != 3 || !p.is_probable_prime() {
            continue;
        }
        let pm1 = &p - &Int::one();
        let nt = &(&pm1 * &pm1) + &(&t * &t);
        for h in 1..=hmax {
            if nt.mod_u64(h) != 0 {
                continue;
            }
            let r = nt.div_floor(&Int::from(h));
            if r.is_probable_prime() {
                return Some(Instance {
                    d,
                    p,
                    t,
                    v,
                    h,
                    r,
                    tries,
                });
            }
        }
    }
    None
}

/// A curve over F_p with `j` a root of `H_d mod p` and trace exactly `t` (checked on random
/// points), and the set of all roots of `H_d mod p` (the crater, when it is F_p-rational).
pub struct CmCurve {
    pub curve: PrimeCurve,
    pub crater: HashSet<Int>,
    pub hilbert_degree: usize,
    pub hilbert_precision: usize,
    pub hilbert_error_log2: f64,
}

struct RootsTask<'a> {
    coeffs: &'a [Int],
    seed: u64,
}
impl api::FieldTask for RootsTask<'_> {
    type Out = Vec<Int>;
    fn run<F: PrimeFieldInt>(self, f: &F) -> Vec<Int> {
        let pp: poly::Poly<F> = self.coeffs.iter().map(|c| f.elem(c)).collect();
        let mut rng = Rng::new(self.seed);
        poly::roots(f, &pp, &mut rng)
            .into_iter()
            .map(|e| f.int(e))
            .collect()
    }
}

struct OrderCheck<'a> {
    curve: &'a PrimeCurve,
    order: &'a Int,
    seed: u64,
}
impl api::FieldTask for OrderCheck<'_> {
    type Out = bool;
    fn run<F: PrimeFieldInt>(self, f: &F) -> bool {
        let c = crate::curve::Curve::new(f.elem(&self.curve.a), f.elem(&self.curve.b));
        let mut rng = Rng::new(self.seed);
        (0..8).all(|_| {
            let pt = crate::curve::random_point_f(f, &c, &mut rng);
            crate::curve::pmul_big(f, &c, &pt, self.order.mag()) == crate::curve::Pt::Inf
        })
    }
}

/// `[order]P = O` for eight random points.
pub fn order_kills(curve: &PrimeCurve, order: &Int, seed: u64) -> bool {
    api::with_field(&curve.p, OrderCheck { curve, order, seed })
}

pub fn build_curve(inst: &Instance, seed: u64) -> Result<CmCurve, String> {
    let hp = hilbert_class_poly(inst.d)?;
    let roots = api::with_field(
        &inst.p,
        RootsTask {
            coeffs: &hp.coeffs,
            seed,
        },
    );
    let degree = hp.coeffs.len() - 1;
    if roots.len() != degree {
        return Err(format!(
            "H_D mod p has {} roots in F_p, expected {degree}",
            roots.len()
        ));
    }
    let p = &inst.p;
    let j0 = roots
        .iter()
        .find(|j| !j.is_zero() && **j != Int::from(1728i64))
        .ok_or("only j = 0, 1728")?;
    // y^2 = x^3 + 3 j (1728 - j) x + 2 j (1728 - j)^2
    let k = (&Int::from(1728i64) - j0).modulo(p);
    let a = (&(&Int::from(3i64) * j0) * &k).modulo(p);
    let b = (&(&(&Int::from(2i64) * j0) * &k) * &k).modulo(p);
    let mut e = PrimeCurve::new(p.clone(), a, b)?;
    let order = &(p + &Int::one()) - &inst.t;
    if !order_kills(&e, &order, seed) {
        e = e.twist();
        if !order_kills(&e, &order, seed ^ 1) {
            return Err("neither the curve nor its twist has trace t".into());
        }
    }
    Ok(CmCurve {
        curve: e,
        crater: roots.into_iter().collect(),
        hilbert_degree: degree,
        hilbert_precision: hp.precision,
        hilbert_error_log2: hp.worst_error_log2,
    })
}

/// A closed walk: the isogeny records in order and the final isomorphism back to the start,
/// `(x, y) -> (x / w^2, y / w^3)`.
#[derive(Clone, Debug)]
pub struct Chain {
    pub steps: Vec<IsogenyRecord>,
    pub js: Vec<Int>,
    pub w: Int,
}

/// Every non-backtracking walk of horizontal isogenies with the given degree sequence that
/// returns to `start`'s j-invariant, with the isomorphism closing it.
pub fn closing_loops(
    start: &PrimeCurve,
    crater: &HashSet<Int>,
    degrees: &[u64],
    seed: u64,
) -> Result<Vec<Chain>, String> {
    let p = &start.p;
    let j0 = start.j();
    let mut out = vec![];
    // depth-first over (curve, previous j, records)
    let mut stack: Vec<(PrimeCurve, Option<Int>, Vec<IsogenyRecord>, Vec<Int>)> =
        vec![(start.clone(), None, vec![], vec![j0.clone()])];
    while let Some((cur, prev_j, recs, js)) = stack.pop() {
        let depth = recs.len();
        if depth == degrees.len() {
            if cur.j() == j0 {
                // w^2 = (b_m a_0) / (a_m b_0)
                let num = (&cur.b * &start.a).modulo(p);
                let den = (&cur.a * &start.b).modulo(p);
                let w2 = (&num * &den.inv_mod(p).ok_or("a or b zero")?).modulo(p);
                let w = Int::sqrt_mod_prime(&w2, p).ok_or("closing isomorphism not over F_p")?;
                out.push(Chain { steps: recs, js, w });
            }
            continue;
        }
        let l = degrees[depth];
        let same_degree_as_prev = depth > 0 && degrees[depth - 1] == l;
        for rec in api::isogenies(&cur, l, seed)? {
            if !rec.verified() || !crater.contains(&rec.j_codomain) {
                continue;
            }
            if same_degree_as_prev && prev_j.as_ref() == Some(&rec.j_codomain) {
                continue; // the dual of the previous step
            }
            let next = PrimeCurve::new(p.clone(), rec.codomain.0.clone(), rec.codomain.1.clone())?;
            let mut r2 = recs.clone();
            r2.push(rec.clone());
            let mut j2 = js.clone();
            j2.push(rec.j_codomain.clone());
            stack.push((next, Some(cur.j()), r2, j2));
        }
    }
    Ok(out)
}

fn eval_poly<F: Field>(f: &F, c: &[F::E], x: F::E) -> (F::E, F::E) {
    // value and derivative by Horner
    let mut v = f.zero();
    let mut dv = f.zero();
    for &a in c.iter().rev() {
        dv = f.add(f.mul(dv, x), v);
        v = f.add(f.mul(v, x), a);
    }
    (v, dv)
}

/// A chain with its coefficients converted into a field `G` containing F_p.
pub struct FieldChain<G: Field> {
    /// per step: (num, den) coefficients low to high
    pub maps: Vec<(Vec<G::E>, Vec<G::E>)>,
    pub w2_inv: G::E,
    pub w3_inv: G::E,
}

impl Chain {
    /// Convert into `G`, given the embedding of F_p integers into `G`.
    pub fn to_field<G: Field>(&self, g: &G, embed: &dyn Fn(&Int) -> G::E) -> FieldChain<G> {
        let maps = self
            .steps
            .iter()
            .map(|r| {
                (
                    r.map_num.iter().map(embed).collect(),
                    r.map_den.iter().map(embed).collect(),
                )
            })
            .collect();
        let w = embed(&self.w);
        let w2 = g.mul(w, w);
        let w3 = g.mul(w2, w);
        FieldChain {
            maps,
            w2_inv: g.inv(w2),
            w3_inv: g.inv(w3),
        }
    }
}

impl<G: Field> FieldChain<G> {
    /// The chain on the twist by `u`: each step `(x, y) -> (u R(x/u), y R'(x/u))`, then the
    /// closing isomorphism.  With `u = 1` it is the chain itself.  Returns None at infinity.
    pub fn eval_twisted(&self, g: &G, u: G::E, x: G::E, y: G::E) -> Option<(G::E, G::E)> {
        let u_inv = g.inv(u);
        let (mut x, mut y) = (x, y);
        for (num, den) in &self.maps {
            let z = g.mul(x, u_inv);
            let (n, dn) = eval_poly(g, num, z);
            let (dd, ddd) = eval_poly(g, den, z);
            if g.is_zero(dd) {
                return None;
            }
            let di = g.inv(dd);
            let r = g.mul(n, di);
            // R' = (n' d - n d') / d^2
            let rp = g.mul(g.sub(g.mul(dn, dd), g.mul(n, ddd)), g.mul(di, di));
            x = g.mul(u, r);
            y = g.mul(y, rp);
        }
        Some((g.mul(x, self.w2_inv), g.mul(y, self.w3_inv)))
    }
}

// ---------------------------------------------------------------- binary GLS number theory

pub fn squarefree(mut n: i64) -> bool {
    let mut q = 2i64;
    while q * q <= n {
        if n % (q * q) == 0 {
            return false;
        }
        if n % q == 0 {
            n /= q;
        }
        q += 1;
    }
    true
}

pub fn prime_form(d: i64, l: i64) -> Option<crate::path::relation::Form> {
    // (l, b, c) with b^2 = d mod 4l, b = d mod 2
    let d = d as i128;
    let l = l as i128;
    for b in 0..2 * l {
        if (b - d).rem_euclid(2) == 0 && (b * b - d).rem_euclid(4 * l) == 0 {
            return Some(
                crate::path::relation::Form {
                    a: l,
                    b,
                    c: (b * b - d) / (4 * l),
                }
                .reduce(),
            );
        }
    }
    None
}

pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    let mut q = 2u64;
    while q * q <= n {
        if n.is_multiple_of(q) {
            return false;
        }
        q += 1;
    }
    true
}

/// Square roots of d modulo 2^k (k >= 3, d = 1 mod 8): the four lifts.
pub fn sqrt_mod_2k(d: &Int, k: usize) -> Vec<Int> {
    let m = Int::one().shl(k);
    let dm = d.modulo(&m);
    let mut x = Int::one();
    for j in 3..k {
        // x^2 = d mod 2^j; lift to mod 2^(j+1)
        let mj1 = Int::one().shl(j + 1);
        if !(&(&x * &x) - &dm).modulo(&mj1).is_zero() {
            x = &x + &Int::one().shl(j - 1);
        }
    }
    let half = Int::one().shl(k - 1);
    let mut out = vec![
        x.modulo(&m),
        (-&x).modulo(&m),
        (&x + &half).modulo(&m),
        (&(-&x) + &half).modulo(&m),
    ];
    out.retain(|r| (&(r * r) - &dm).modulo(&m).is_zero());
    out.sort();
    out.dedup();
    out
}

/// All (A, B), B > 0, with A^2 + |d| B^2 = 4n, n = l 2^k (l = 1 or an odd prime with (d/l) = 1),
/// found by Cornacchia's reduction from every square root of d mod 4n; each result is verified.
pub fn cornacchia(d: i64, l: u64, k: usize) -> Vec<(Int, Int)> {
    let ad = Int::from(-d);
    let n = &Int::from(l) * &Int::one().shl(k);
    let four_n = &n * &Int::from(4i64);
    let di = Int::from(d);
    let r2 = sqrt_mod_2k(&di, k + 2);
    let rl: Vec<Int> = if l == 1 {
        vec![Int::zero()]
    } else {
        (0..l)
            .filter(|&x| (x as i128 * x as i128 - d as i128).rem_euclid(l as i128) == 0)
            .map(Int::from)
            .collect()
    };
    let m2 = Int::one().shl(k + 2);
    let ml = Int::from(l);
    let bound = (&n * &Int::from(4i64)).isqrt();
    let mut out: Vec<(Int, Int)> = vec![];
    for a2 in &r2 {
        for al in &rl {
            // CRT: x = a2 mod 2^(k+2), x = al mod l
            let x0 = if l == 1 {
                a2.clone()
            } else {
                let t = (&(al - a2) * &m2.inv_mod(&ml).unwrap()).modulo(&ml);
                (a2 + &(&t * &m2)).modulo(&four_n)
            };
            for start in [x0.clone(), (-&x0).modulo(&four_n)] {
                let (mut a, mut b) = (
                    &n * &Int::from(2i64),
                    start.modulo(&(&n * &Int::from(2i64))),
                );
                let mut tried = 0;
                while tried < 4096 {
                    tried += 1;
                    if b <= bound {
                        let rest = &four_n - &(&b * &b);
                        if !rest.is_neg() && rest.modulo(&ad).is_zero() {
                            let c = rest.div_floor(&ad);
                            if c.is_square() && !c.is_zero() {
                                let bb = c.isqrt();
                                if &(&b * &b) + &(&ad * &(&bb * &bb)) == four_n {
                                    out.push((b.clone(), bb.clone()));
                                    out.push((-&b, bb));
                                }
                            }
                        }
                    }
                    if b.is_zero() {
                        break;
                    }
                    let r = a.modulo(&b);
                    a = b;
                    b = r;
                }
            }
        }
    }
    out.sort();
    out.dedup();
    out
}

pub fn v2(x: &Int) -> usize {
    if x.is_zero() {
        return usize::MAX;
    }
    let mut v = 0;
    let mut y = x.abs();
    while y.is_even() {
        y = y.shr(1);
        v += 1;
    }
    v
}

// ---------------------------------------------------------------- decomposition

/// Babai rounding with precomputed constants: c_j = round(k g_j / 2^m), g_j = round(C_j0 2^m / det).
pub struct Babai {
    pub basis: Vec<Vec<Int>>,
    g: Vec<Int>,
    m: usize,
}

fn det3(m: &[[Int; 3]; 3]) -> Int {
    let t = |a: &Int, b: &Int, c: &Int, d: &Int| &(a * b) - &(c * d);
    &(&(&m[0][0] * &t(&m[1][1], &m[2][2], &m[1][2], &m[2][1]))
        - &(&m[0][1] * &t(&m[1][0], &m[2][2], &m[1][2], &m[2][0])))
        + &(&m[0][2] * &t(&m[1][0], &m[2][1], &m[1][1], &m[2][0]))
}

impl Babai {
    pub fn new(basis: Vec<Vec<Int>>, rbits: usize) -> Babai {
        let n = basis.len();
        let (cof, det): (Vec<Int>, Int) = if n == 2 {
            let cof = vec![basis[1][1].clone(), -&basis[0][1]];
            let det = &(&basis[0][0] * &basis[1][1]) - &(&basis[0][1] * &basis[1][0]);
            (cof, det)
        } else {
            let mut cof = vec![];
            for j in 0..4 {
                let rows: Vec<usize> = (0..4).filter(|&i| i != j).collect();
                let mut mm: [[Int; 3]; 3] =
                    std::array::from_fn(|_| std::array::from_fn(|_| Int::zero()));
                for (ri, &i) in rows.iter().enumerate() {
                    for c in 1..4 {
                        mm[ri][c - 1] = basis[i][c].clone();
                    }
                }
                let d3 = det3(&mm);
                cof.push(if j % 2 == 0 { d3 } else { -&d3 });
            }
            let det = (0..4).fold(Int::zero(), |acc, j| &acc + &(&basis[j][0] * &cof[j]));
            (cof, det)
        };
        let m = 2 * rbits + 16;
        let g = cof.iter().map(|c| c.shl(m).div_round(&det)).collect();
        Babai { basis, g, m }
    }
    pub fn decompose(&self, k: &Int) -> Vec<Int> {
        let n = self.basis.len();
        let mut out = vec![Int::zero(); n];
        out[0] = k.clone();
        let half = Int::one().shl(self.m - 1);
        for j in 0..n {
            let t = &(k * &self.g[j]) + &half;
            let c = t.div_floor(&Int::one().shl(self.m));
            if c.is_zero() {
                continue;
            }
            for i in 0..n {
                out[i] = &out[i] - &(&c * &self.basis[j][i]);
            }
        }
        out
    }
}

/// Lagrange–Gauss-reduced basis of {(k0, k1) : k0 + k1 mu = 0 mod r}.
pub fn gauss2(r: &Int, mu: &Int) -> Vec<Vec<Int>> {
    // Lagrange-Gauss reduction of (r, 0), (-mu, 1)
    let mut a = (r.clone(), Int::zero());
    let mut b = ((-mu).modulo(r), Int::one());
    let n2 = |v: &(Int, Int)| &(&v.0 * &v.0) + &(&v.1 * &v.1);
    loop {
        if n2(&a) < n2(&b) {
            std::mem::swap(&mut a, &mut b);
        }
        let dot = &(&a.0 * &b.0) + &(&a.1 * &b.1);
        let q = dot.div_round(&n2(&b));
        if q.is_zero() {
            break;
        }
        a = (&a.0 - &(&q * &b.0), &a.1 - &(&q * &b.1));
    }
    vec![vec![a.0, a.1], vec![b.0, b.1]]
}

/// Exact-LLL-reduced basis of {k : k0 + k1 mu + k2 lam + k3 lam mu = 0 mod r}.
pub fn lll4(r: &Int, mu: &Int, lam: &Int) -> Vec<Vec<Int>> {
    let lm = (lam * mu).modulo(r);
    let z = Int::zero;
    let o = Int::one;
    let rows: Vec<Vec<Int>> = vec![
        vec![r.clone(), z(), z(), z()],
        vec![(-mu).modulo(r), o(), z(), z()],
        vec![(-lam).modulo(r), z(), o(), z()],
        vec![(-&lm).modulo(r), z(), z(), o()],
    ];
    let mut gram: [[Int; 4]; 4] = std::array::from_fn(|_| std::array::from_fn(|_| Int::zero()));
    for i in 0..4 {
        for j in 0..4 {
            gram[i][j] = (0..4).fold(Int::zero(), |acc, c| &acc + &(&rows[i][c] * &rows[j][c]));
        }
    }
    let h = crate::quat::lll_gram(&gram);
    (0..4)
        .map(|i| {
            (0..4)
                .map(|c| (0..4).fold(Int::zero(), |acc, k| &acc + &(&h[i][k] * &rows[k][c])))
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cheapest_loops_reproduce_known_values() {
        // values certified natively in research/isogeny_loop_glv_20261011 and in the paper
        for (d, cost, norm) in [
            (-23i64, 37.5, 6u64),
            (-71, 60.0, 18),
            (-619, 127.5, 175),
            (-3019, 210.0, 1625),
        ] {
            let (l, _, _, _) = cheapest_loop(d);
            assert_eq!(l.cost, cost, "D = {d}");
            assert_eq!(l.norm, norm, "D = {d}");
            assert_eq!((l.a * l.a - d * l.b * l.b) / 4, norm as i64);
        }
    }

    #[test]
    fn loop_closes_and_maps_points_on_a_small_instance() {
        // D = -23 (h = 3), cheapest loop 2*3: a full build at a 40-bit prime
        let d = -23;
        let (l, _, _, _) = cheapest_loop(d);
        let inst = find_instance(d, 40, &[3], 64, 7, 2_000_000).expect("instance");
        let cm = build_curve(&inst, 11).expect("curve");
        assert_eq!(cm.hilbert_degree, 3);
        let chains = closing_loops(&cm.curve, &cm.crater, &l.degree_sequence(), 3).expect("walk");
        assert!(!chains.is_empty(), "no closing loop");
    }
}
