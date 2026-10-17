//! A multiple of the conductor of `End(E)` from a cycle of horizontal
//! `l`-isogenies, `l` a small prime, for curves whose CM field has a small
//! class group.
//!
//! Take `l` split in `K` (`(d_K/l) = 1`) and prime to `v` (`t² − 4q =
//! v²·d_K`).  Then Frobenius has two distinct eigenvalues on `E[l]`, `E`
//! has exactly two rational `l`-isogenies, both horizontal, and so does
//! every curve they reach: the `l`-isogeny graph through `E` is a cycle,
//! of length the order of a prime above `l` in the class group of
//! `End(E)`.  The walk follows it on `j`-invariants, the next one being
//! the root of `Φ_l(j_i, Y)` in `GF(q)` other than `j_{i−1}`, each root
//! required simple (so it is one isogeny, the previous one's dual when it
//! is `j_{i−1}`).
//!
//! When it is back at `j(E)` after `r` steps the composition is an
//! endomorphism `α` of degree `l^r`, defined over `GF(q)` (an isogenous
//! curve with the same `j` is `E`, not its twist, which has trace `−t`),
//! with cyclic kernel since no step undoes the one before.  So `α` is
//! primitive in `O_K` (`l ∤ v`, so `End(E)` and `O_K` agree at `l`) of norm
//! `l^r`: `(α)` is `𝔩^r` or `𝔩̄^r`, and `α = (x + y√d_K)/2` with
//! `x² − d_K y² = 4l^r` and `l ∤ gcd(x, y)`, which Cornacchia's algorithm
//! finds, `|y|` the same for both ideals.  `Z[α]` is the order of conductor
//! `|y|`, so the conductor `f` of `End(E)` divides `y`.  For `d_K = −3, −4`
//! `α` is one of the unit multiples of that generator, so `f` divides one
//! of their `y`.
//!
//! The walk is the proof; nothing is assumed about the class group.  A
//! walk that has not closed within [`WALK_COST`] or [`MAX_STEPS`] says
//! nothing.

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::points::{multiplication_cost, pgcd, sqrt, BinaryW};
use super::torsion::{monic, ppow, psub, trim, Big2, Field, Fp, FpK, Poly, Small2};
use super::Model;
use crate::cryptanalysis::semaev_decomp::Gf2;

/// `Φ_l(X, Y)` for the primes `3 ≤ l ≤ 19`, from PARI/GP's `polmodular`.
const MODULAR_POLYNOMIALS: &str = include_str!("modular_polynomials.txt");

/// Work allowed for one walk, in base-field word products.
pub const WALK_COST: u64 = 4_000_000_000;

/// Steps allowed for one walk.
pub const MAX_STEPS: u32 = 20_000;

/// A closed cycle of `l`-isogenies through `E`.
#[derive(Clone, Debug)]
pub struct Cycle {
    /// The isogeny degree `l`.
    pub degree: u64,
    /// The cycle length `r`.
    pub steps: u32,
    /// The conductor of `End(E)` divides one of these (one unless
    /// `d_K = −3, −4`); empty for a cycle too long to be worth it.
    pub multiples: Vec<BigUint>,
}

/// The coefficients `(i, j, c)` of `c X^i Y^j` in `Φ_l`, both orders.
fn modular_polynomial(l: u64) -> Vec<(usize, usize, BigInt)> {
    let mut out = Vec::new();
    for line in MODULAR_POLYNOMIALS.lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut it = line.split(' ');
        let (Some(ll), Some(i), Some(j), Some(c)) = (it.next(), it.next(), it.next(), it.next())
        else {
            continue;
        };
        if ll.parse::<u64>().ok() != Some(l) {
            continue;
        }
        let (i, j): (usize, usize) = (i.parse().expect("degree"), j.parse().expect("degree"));
        let c: BigInt = c.parse().expect("coefficient");
        if i != j {
            out.push((j, i, c.clone()));
        }
        out.push((i, j, c));
    }
    out
}

/// The primes `l ≤ 19` split in `K`, prime to `v` and to `q`, ascending.
fn degrees(dk: &BigInt, v: &BigUint, q: &BigUint) -> Vec<u64> {
    [3u64, 5, 7, 11, 13, 17, 19].into_iter().filter(|&l| {
        let d = dk.mod_floor(&BigInt::from(l)).to_u64().expect("small");
        !(v % l).is_zero()
            && !(q % l).is_zero()
            && d != 0
            && (1..l).any(|x| x * x % l == d)
    })
    .collect()
}

/// The walks of up to [`DEGREES_TRIED`] usable degrees `l`, at their
/// default budget; closed ones bound the conductor, open ones exclude every
/// order whose class number they outlast.
#[derive(Clone, Debug, Default)]
pub struct Walks {
    pub closed: Vec<Cycle>,
    pub open: Vec<(u64, u32)>,
}

pub fn walks(model: &Model, dk: &BigInt, v: &BigUint) -> Walks {
    let mut out = Walks::default();
    // A cycle's length is the order of a class in a group of order about
    // √|d_K| or more: past MAX_DISC_BITS it cannot close in MAX_STEPS.
    if dk.bits() > MAX_DISC_BITS {
        return out;
    }
    for l in degrees(dk, v, &model.q()) {
        if out.closed.len() + out.open.len() >= DEGREES_TRIED {
            break;
        }
        let Some(steps) = affordable(model, l, 1) else {
            break;
        };
        match walk_degree(model, dk, l, steps) {
            Some(Walked::Closed(c)) => out.closed.push(c),
            Some(Walked::Open { degree, steps }) => out.open.push((degree, steps)),
            None => {}
        }
    }
    out
}

/// The bound the cycles put on `v_ℓ(f)`: each caps it by the largest
/// valuation among its multiples, together by the smallest such cap.
pub fn cap(cycles: &[Cycle], ell: &BigUint) -> u32 {
    let val = |n: &BigUint| {
        let (mut n, mut e) = (n.clone(), 0u32);
        while !n.is_zero() && (&n % ell).is_zero() {
            n /= ell;
            e += 1;
        }
        e
    };
    cycles
        .iter()
        .filter(|cy| !cy.multiples.is_empty())
        .map(|cy| cy.multiples.iter().map(val).max().unwrap_or(0))
        .min()
        .unwrap_or(u32::MAX)
}

/// "a cycle of r l-isogenies", or "cycles of …" for several.
pub fn describe(cycles: &[Cycle]) -> String {
    let walks: Vec<String> = cycles
        .iter()
        .map(|cy| format!("{} {}-isogenies", cy.steps, cy.degree))
        .collect();
    if walks.len() == 1 {
        format!("a cycle of {}", walks[0])
    } else {
        format!("cycles of {}", walks.join(", "))
    }
}

/// The largest norm `l^r`, in bits, whose generator [`walk_degree`] finds
/// by Cornacchia's algorithm.
pub const CORNACCHIA_BITS: f64 = 4096.0;

/// The largest `|d_K|`, in bits, for which [`walks`]
/// walks.
pub const MAX_DISC_BITS: u64 = 48;

/// Degrees `l` that [`walks`] walks.
pub const DEGREES_TRIED: usize = 3;

/// The integer `c` in the field of characteristic `p`.
fn embed<F: Field>(f: &F, c: &BigInt, p: &BigUint) -> F::E {
    let r = c.mod_floor(&BigInt::from(p.clone())).to_biguint().expect("non-negative");
    let base = f.int(1 << 32);
    let mut acc = f.zero();
    for digit in r.to_u32_digits().iter().rev() {
        acc = f.add(&f.mul(&acc, &base), &f.int(u64::from(*digit)));
    }
    acc
}

/// Roots of a monic quadratic `Y² + g₁Y + g₀`, one of them.
enum Roots2<'a, F: Field> {
    Odd { order: &'a BigUint },
    Binary { group: BinaryW<'a, F> },
}

fn quadratic_root<F: Field>(f: &F, how: &Roots2<'_, F>, g1: &F::E, g0: &F::E) -> Option<F::E> {
    match how {
        Roots2::Odd { order } => {
            let disc = f.sub(&f.mul(g1, g1), &f.mul(&f.int(4), g0));
            let s = sqrt(f, &disc, order)?;
            Some(f.mul(&f.sub(&s, g1), &f.inv(&f.int(2))))
        }
        Roots2::Binary { group } => {
            // Y = g₁z, z² + z = g₀/g₁².
            if f.is_zero(g1) {
                return None;
            }
            let c = f.mul(g0, &f.inv(&f.mul(g1, g1)));
            Some(f.mul(g1, &group.solve_artin_schreier(&c)?))
        }
    }
}

/// Walk the cycle of `l`-isogenies from `j0`: `Ok(r)` when it closes after
/// `r` steps, `Err(n)` when it has not after `n = max_steps` clean steps;
/// `None` at a step that is not a clean pair of simple roots.
fn walk<F: Field>(
    f: &F,
    p: &BigUint,
    q: &BigUint,
    j0: &F::E,
    l: u64,
    how: &Roots2<'_, F>,
    max_steps: u32,
) -> Option<Result<u32, u32>> {
    let phi: Vec<(usize, usize, F::E)> = modular_polynomial(l)
        .into_iter()
        .map(|(i, j, c)| (i, j, embed(f, &c, p)))
        .collect();
    let top = l as usize + 1;
    let y = vec![f.zero(), f.one()];
    let eval = |g: &Poly<F::E>, x: &F::E| {
        g.iter()
            .rev()
            .fold(f.zero(), |acc, c| f.add(&f.mul(&acc, x), c))
    };
    let mut prev: Option<F::E> = None;
    let mut cur = j0.clone();
    for step in 1..=max_steps {
        // Φ_l(cur, Y), monic of degree l + 1.
        let mut powers = vec![f.one()];
        for _ in 0..top {
            let last = powers.last().expect("non-empty");
            powers.push(f.mul(last, &cur));
        }
        let mut phi_y = vec![f.zero(); top + 1];
        for (i, j, c) in &phi {
            phi_y[*j] = f.add(&phi_y[*j], &f.mul(c, &powers[*i]));
        }
        let phi_y = trim(f, phi_y);
        if phi_y.len() != top + 1 {
            return None;
        }
        let phi_y = monic(f, &phi_y);
        // Its distinct roots in GF(q): gcd(Y^q − Y, Φ_l(cur, Y)).
        let y_q = ppow(f, &y, q, &phi_y);
        let g = pgcd(f, &phi_y, &psub(f, &y_q, &y));
        if g.len() != 3 {
            return None;
        }
        let g = monic(f, &g);
        // Both roots simple: gcd(g, ∂Φ/∂Y) = 1.
        let derivative = trim(
            f,
            phi_y
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, c)| f.mul(c, &f.int(i as u64)))
                .collect(),
        );
        if pgcd(f, &g, &derivative).len() != 1 {
            return None;
        }
        let next = match &prev {
            Some(back) => {
                if !f.is_zero(&eval(&g, back)) {
                    return None;
                }
                // The roots sum to −g₁.
                f.sub(&f.neg(&g[1]), back)
            }
            None => quadratic_root(f, how, &g[1], &g[0])?,
        };
        if next == *j0 {
            return Some(Ok(step));
        }
        prev = Some(cur);
        cur = next;
    }
    Some(Err(max_steps))
}

/// [`walk`] over `GF(p)`, `p < 2^63`, in plain word arithmetic: the same
/// steps, without the allocation of the generic polynomials.
fn walk_word(p: u64, j0: u64, l: u64, max_steps: u32) -> Option<Result<u32, u32>> {
    let m = |a: u64, b: u64| (u128::from(a) * u128::from(b) % u128::from(p)) as u64;
    let add = |a: u64, b: u64| ((u128::from(a) + u128::from(b)) % u128::from(p)) as u64;
    let sub = |a: u64, b: u64| add(a, p - b % p);
    let pow = |mut b: u64, mut e: u64| {
        let mut acc = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                acc = m(acc, b);
            }
            b = m(b, b);
            e >>= 1;
        }
        acc
    };
    let inv = |a: u64| pow(a, p - 2);
    let strip = |mut v: Vec<u64>| {
        while v.last() == Some(&0) {
            v.pop();
        }
        v
    };
    // a mod b, b non-zero.
    let rem = |a: &[u64], b: &[u64]| {
        let mut r = a.to_vec();
        let lead = inv(*b.last().expect("non-zero"));
        while r.len() >= b.len() {
            let c = m(*r.last().expect("non-empty"), lead);
            let shift = r.len() - b.len();
            for (i, &x) in b.iter().enumerate() {
                r[shift + i] = sub(r[shift + i], m(c, x));
            }
            r = strip(r);
            if r.is_empty() {
                break;
            }
        }
        r
    };
    let gcd = |a: &[u64], b: &[u64]| {
        let (mut a, mut b) = (strip(a.to_vec()), strip(b.to_vec()));
        while !b.is_empty() {
            let r = rem(&a, &b);
            a = b;
            b = r;
        }
        let lead = inv(*a.last().expect("non-zero"));
        a.iter().map(|&c| m(c, lead)).collect::<Vec<u64>>()
    };
    let mulmod = |a: &[u64], b: &[u64], f: &[u64]| {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let mut t = vec![0u64; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                t[i + j] = add(t[i + j], m(x, y));
            }
        }
        rem(&strip(t), f)
    };
    let phi: Vec<(usize, usize, u64)> = modular_polynomial(l)
        .into_iter()
        .map(|(i, j, c)| {
            let r = c.mod_floor(&BigInt::from(p)).to_u64().expect("reduced");
            (i, j, r)
        })
        .collect();
    let top = l as usize + 1;
    let half = (p - 1) / 2;
    let mut prev: Option<u64> = None;
    let mut cur = j0;
    for step in 1..=max_steps {
        let mut powers = vec![1u64];
        for _ in 0..top {
            let last = *powers.last().expect("non-empty");
            powers.push(m(last, cur));
        }
        let mut f = vec![0u64; top + 1];
        for &(i, j, c) in &phi {
            f[j] = add(f[j], m(c, powers[i]));
        }
        let f = strip(f);
        if f.len() != top + 1 {
            return None;
        }
        // Y^p mod f.
        let mut acc = vec![1u64];
        let base = vec![0u64, 1];
        for i in (0..64 - p.leading_zeros()).rev() {
            acc = mulmod(&acc, &acc, &f);
            if (p >> i) & 1 == 1 {
                acc = mulmod(&acc, &base, &f);
            }
        }
        let mut d = acc;
        d.resize(d.len().max(2), 0);
        d[1] = sub(d[1], 1);
        let g = gcd(&f, &strip(d));
        if g.len() != 3 {
            return None;
        }
        let derivative: Vec<u64> = f
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &c)| m(c, i as u64 % p))
            .collect();
        if gcd(&g, &strip(derivative)).len() != 1 {
            return None;
        }
        let next = match prev {
            Some(back) => {
                if add(add(m(back, back), m(g[1], back)), g[0]) != 0 {
                    return None;
                }
                sub(sub(0, g[1]), back)
            }
            None => {
                // (−g₁ + √(g₁² − 4g₀))/2 by Tonelli–Shanks.
                let disc = sub(m(g[1], g[1]), m(4, g[0]));
                if disc != 0 && pow(disc, half) != 1 {
                    return None;
                }
                let root = {
                    let (mut q, mut s) = (p - 1, 0u32);
                    while q % 2 == 0 {
                        q /= 2;
                        s += 1;
                    }
                    let z = (2..p).find(|&z| pow(z, half) == p - 1)?;
                    let (mut mm, mut c, mut t, mut r) =
                        (s, pow(z, q), pow(disc, q), pow(disc, q.div_ceil(2)));
                    while t > 1 {
                        let mut i = 0;
                        let mut tt = t;
                        while tt != 1 {
                            tt = m(tt, tt);
                            i += 1;
                        }
                        let mut b = c;
                        for _ in 0..(mm - i - 1) {
                            b = m(b, b);
                        }
                        mm = i;
                        c = m(b, b);
                        t = m(t, c);
                        r = m(r, b);
                    }
                    if disc == 0 { 0 } else { r }
                };
                m(sub(root, g[1]), inv(2))
            }
        };
        if next == j0 {
            return Some(Ok(step));
        }
        prev = Some(cur);
        cur = next;
    }
    Some(Err(max_steps))
}

/// `α = (x + y√d)/2` as `(x, y)`.
type Quad = (BigInt, BigInt);

fn quad_mul(d: &BigInt, a: &Quad, b: &Quad) -> Quad {
    let x = (&a.0 * &b.0 + d * &a.1 * &b.1) / 2;
    let y = (&a.0 * &b.1 + &b.0 * &a.1) / 2;
    (x, y)
}

/// A primitive solution of `x² − d y² = 4l^r` (`d < 0` a fundamental
/// discriminant, `l` an odd prime split in it), by Cornacchia's algorithm
/// from a square root of `d` mod `4l^r`; checked before it is returned.
fn primitive_solution(d: &BigInt, l: u64, r: u32) -> Option<Quad> {
    let li = BigInt::from(l);
    let m = li.pow(r);
    let s0 = (1..l).find(|&x| BigInt::from(x * x).mod_floor(&li) == d.mod_floor(&li))?;
    // Hensel, doubling the precision: s² ≡ d mod l^k for k up to r.
    let mut s = BigInt::from(s0);
    let mut k = 1u32;
    while k < r {
        k = (2 * k).min(r);
        let modulus = li.pow(k);
        let fix = (&s * &s - d).mod_floor(&modulus);
        let inv = (BigInt::from(2) * &s).extended_gcd(&modulus);
        if !inv.gcd.is_one() {
            return None;
        }
        s = (&s - fix * inv.x).mod_floor(&modulus);
    }
    for root in [s.clone(), (&m - &s).mod_floor(&m)] {
        // root ≡ d mod 2, so root² ≡ d mod 4m.
        let b0 = if (&root - d).is_odd() {
            &root + &m
        } else {
            root
        };
        let limit = (BigInt::from(4) * &m).sqrt();
        let (mut a, mut b) = (BigInt::from(2) * &m, b0);
        while b > limit {
            let r = a.mod_floor(&b);
            a = b;
            b = r;
        }
        let rest = BigInt::from(4) * &m - &b * &b;
        if rest.is_negative() || !rest.mod_floor(&-d).is_zero() {
            continue;
        }
        let c = rest / -d;
        let y = c.sqrt();
        if &y * &y != c || y.is_zero() {
            continue;
        }
        let ok = &b * &b - d * &y * &y == BigInt::from(4) * &m
            && !(b.mod_floor(&li).is_zero() && y.mod_floor(&li).is_zero());
        if ok {
            return Some((b, y));
        }
    }
    None
}

/// What a walk of `l`-isogenies found.
#[derive(Clone, Debug)]
pub enum Walked {
    /// It closed: a cycle.
    Closed(Cycle),
    /// It took `steps` clean steps without coming back.
    Open { degree: u64, steps: u32 },
}

/// The steps a walk of degree `l` may take within [`WALK_COST`] times
/// `scale`, at most [`MAX_STEPS`] times `scale`.
pub fn affordable(model: &Model, l: u64, scale: u64) -> Option<u32> {
    let step_cost = multiplication_cost(model, 1, model.q().bits())?
        .saturating_mul((l + 2) * (l + 2))
        .div_ceil(3);
    Some(
        (WALK_COST.saturating_mul(scale) / step_cost.max(1))
            .min(u64::from(MAX_STEPS) * scale)
            .min(u64::from(u32::MAX)) as u32,
    )
}

/// Walk `l`-isogenies from the curve (`l` from [`degrees`]) for at most
/// `max_steps`.  `None` when the model has no arithmetic here or a step is
/// not a clean pair of simple roots.
pub fn walk_degree(model: &Model, dk: &BigInt, l: u64, max_steps: u32) -> Option<Walked> {
    if dk.sign() != Sign::Minus || max_steps < 3 {
        return None;
    }
    let q = model.q();
    let result = match model {
        Model::Prime { p, a, b } => {
            if p <= &BigUint::from(3u8) {
                return None;
            }
            if let Some(f) = FpK::new(p, &[BigUint::zero()]) {
                // GF(p) in one word.
                let el = |c: &BigUint| f.int(c.to_u64().expect("below p < 2^63"));
                let j = short_j(&f, &el(a), &el(b))?;
                walk_word(f.p, j[0], l, max_steps)?
            } else {
                let f = Fp(p.clone());
                let j = short_j(&f, a, b)?;
                walk(&f, p, &q, &j, l, &Roots2::Odd { order: &q }, max_steps)?
            }
        }
        Model::Extension {
            p, modulus, a, b, ..
        } => {
            if p <= &BigUint::from(3u8) {
                return None;
            }
            let f = FpK::new(p, modulus)?;
            let (a, b) = (f.element(a)?, f.element(b)?);
            let j = short_j(&f, &a, &b)?;
            walk(&f, p, &q, &j, l, &Roots2::Odd { order: &q }, max_steps)?
        }
        Model::Binary {
            n: deg, irr, b, ..
        } => {
            let two = BigUint::from(2u8);
            if *deg <= 63 {
                let f = Small2(Gf2::new(irr));
                let b = f.0.from_element(b);
                let how = Roots2::Binary {
                    group: BinaryW {
                        f: &f,
                        a: f.zero(),
                        b: f.zero(),
                        m: u64::from(*deg),
                    },
                };
                walk(&f, &two, &q, &f.inv(&b), l, &how, max_steps)?
            } else {
                let f = Big2 {
                    n: *deg,
                    irr: irr.clone(),
                };
                let how = Roots2::Binary {
                    group: BinaryW {
                        f: &f,
                        a: f.zero(),
                        b: f.zero(),
                        m: u64::from(*deg),
                    },
                };
                walk(&f, &two, &q, &f.inv(b), l, &how, max_steps)?
            }
        }
        Model::Other { .. } => return None,
    };
    let steps = match result {
        Ok(r) => r,
        Err(n) => return Some(Walked::Open { degree: l, steps: n }),
    };
    if (l as f64).log2() * f64::from(steps) > CORNACCHIA_BITS {
        // The conductor of Z[α] is near l^{r/2}: no use as a bound, and
        // costly to find.  The cycle length still counts.
        return Some(Walked::Closed(Cycle {
            degree: l,
            steps,
            multiples: Vec::new(),
        }));
    }
    let alpha = primitive_solution(dk, l, steps)?;
    // Unit multiples for d_K = −4 (i) and −3 (a sixth root of unity).
    let unit = match dk.to_i64() {
        Some(-4) => Some(((BigInt::zero(), BigInt::one()), 4)),
        Some(-3) => Some(((BigInt::one(), BigInt::one()), 6)),
        _ => None,
    };
    let mut multiples = vec![alpha.1.magnitude().clone()];
    if let Some((u, order)) = unit {
        let mut a = alpha.clone();
        for _ in 1..order {
            a = quad_mul(dk, &a, &u);
            if !a.1.is_zero() {
                multiples.push(a.1.magnitude().clone());
            }
        }
        multiples.sort();
        multiples.dedup();
    }
    Some(Walked::Closed(Cycle {
        degree: l,
        steps,
        multiples,
    }))
}

/// The class number of the fundamental discriminant `d < 0`, by counting
/// reduced forms `(a, b, c)`, `|b| ≤ a ≤ c`, `b ≥ 0` when `|b| = a` or
/// `a = c`; `None` past [`CLASS_NUMBER_BITS`].
pub fn class_number(d: &BigInt) -> Option<u64> {
    if d.sign() != Sign::Minus || d.bits() > CLASS_NUMBER_BITS {
        return None;
    }
    let n = d.magnitude().to_u64()?;
    let mut h = 0u64;
    let mut a = 1u64;
    while 3 * a * a <= n {
        // b ≡ d (mod 2), −a < b ≤ a.
        let start = if n % 2 == 1 { 1 } else { 0 };
        let mut b = start;
        while b <= a {
            let num = b * b + n;
            if num % (4 * a) == 0 {
                let c = num / (4 * a);
                if c >= a {
                    // (a, b, c), and (a, −b, c) unless b = 0, b = a or a = c.
                    h += if b == 0 || b == a || a == c { 1 } else { 2 };
                }
            }
            b += 2;
        }
        a += 1;
    }
    Some(h)
}

/// The largest `|d_K|`, in bits, [`class_number`] counts for.
pub const CLASS_NUMBER_BITS: u64 = 32;

/// `h(O_f) = h_K·f·Π_{p|f}(1 − (d_K/p)/p) / [O_K^* : O_f^*]` for the order
/// of conductor `f = Π p^e`.
pub fn order_class_number(dk: &BigInt, hk: u64, f: &[(BigUint, u32)]) -> BigUint {
    let mut h = BigUint::from(hk);
    let mut nontrivial = false;
    for (p, e) in f {
        if *e == 0 {
            continue;
        }
        nontrivial = true;
        let chi = kronecker(dk, p);
        h *= p.pow(e - 1);
        h = match chi {
            1 => h * (p - 1u8),
            -1 => h * (p + 1u8),
            _ => h * p,
        };
    }
    if nontrivial {
        match dk.to_i64() {
            Some(-4) => h /= 2u8,
            Some(-3) => h /= 3u8,
            _ => {}
        }
    }
    h
}

/// The Kronecker symbol `(d/p)` for a prime `p`.
fn kronecker(d: &BigInt, p: &BigUint) -> i32 {
    if p == &BigUint::from(2u8) {
        return match d.mod_floor(&BigInt::from(8)).to_u8().expect("small") {
            1 | 7 => 1,
            3 | 5 => -1,
            _ => 0,
        };
    }
    let pi = BigInt::from(p.clone());
    let r = d.mod_floor(&pi);
    if r.is_zero() {
        return 0;
    }
    let e = (&pi - 1u8) >> 1u8;
    if r.modpow(&e, &pi).is_one() {
        1
    } else {
        -1
    }
}

/// `j = 1728·4a³/(4a³ + 27b²)` of `y² = x³ + ax + b`; `None` for `j = 0`
/// or `1728`, whose extra automorphisms the walk does not handle.
fn short_j<F: Field>(f: &F, a: &F::E, b: &F::E) -> Option<F::E> {
    if f.is_zero(a) || f.is_zero(b) {
        return None;
    }
    let a3 = f.mul(&f.int(4), &f.mul(a, &f.mul(a, a)));
    let den = f.add(&a3, &f.mul(&f.int(27), &f.mul(b, b)));
    Some(f.mul(&f.mul(&f.int(1728), &a3), &f.inv(&den)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::curve_traits::torsion::{frobenius_is_scalar, scalar_on_power};

    /// The word-arithmetic walk takes the generic walk's steps.
    #[test]
    fn the_word_walk_matches_the_generic_walk() {
        let mut compared = 0;
        for p in [1009u64, 2003, 10007, 65537] {
            let f = Fp(BigUint::from(p));
            let q = BigUint::from(p);
            for a in 1..8u64 {
                for b in 1..8u64 {
                    let Some(j) = short_j(&f, &BigUint::from(a), &BigUint::from(b)) else {
                        continue;
                    };
                    for l in [3u64, 5, 7] {
                        let generic = walk(&f, &q, &q, &j, l, &Roots2::Odd { order: &q }, 200);
                        let word = walk_word(p, j.to_u64().expect("small"), l, 200);
                        assert_eq!(generic, word, "p={p} a={a} b={b} l={l}");
                        compared += usize::from(generic.is_some());
                    }
                }
            }
        }
        assert!(compared > 100, "{compared}");
    }

    #[test]
    fn class_numbers_match_known_values() {
        for (d, h) in [(-3i64, 1u64), (-4, 1), (-7, 1), (-8, 1), (-15, 2), (-20, 2), (-23, 3),
            (-47, 5), (-71, 7), (-84, 4), (-163, 1), (-199, 9), (-260, 8), (-1999, 27)] {
            assert_eq!(class_number(&BigInt::from(d)), Some(h), "h({d})");
        }
        // h(−3·3²) = 1, h(−4·5²) = 2, h(−7·2²) = 1, h(−23·3²) = 3·(3 − 1) = 6.
        let o = |d: i64, p: u32, e: u32, hk: u64| {
            order_class_number(&BigInt::from(d), hk, &[(BigUint::from(p), e)])
        };
        assert_eq!(o(-3, 3, 1, 1), BigUint::from(1u8));
        assert_eq!(o(-4, 5, 1, 1), BigUint::from(2u8));
        assert_eq!(o(-7, 2, 1, 1), BigUint::from(1u8));
        assert_eq!(o(-23, 3, 1, 3), BigUint::from(6u8));
    }

    /// Cornacchia's solutions have the norm and are primitive at `l`.
    #[test]
    fn primitive_solutions_have_the_norm() {
        let mut found = 0;
        for (d, l) in [(-7i64, 11u64), (-23, 3), (-15, 17), (-104, 5), (-4, 5), (-3, 7), (-7, 2 + 1)] {
            let d = BigInt::from(d);
            for r in 1..10 {
                if let Some((x, y)) = primitive_solution(&d, l, r) {
                    assert_eq!(&x * &x - &d * &y * &y, BigInt::from(4) * BigInt::from(l).pow(r));
                    found += 1;
                }
            }
        }
        // 11 = N((3 + √−7)/2) and 4·11 = 3² + 7·... : x² + 7y² = 44 at (±4, ±2)
        // is primitive at 11.
        assert!(primitive_solution(&BigInt::from(-7), 11, 1).is_some());
        assert!(found > 10, "{found}");
    }

    fn pow_mod(b: i128, e: u32, p: i128) -> i128 {
        let mut acc = 1i128;
        for _ in 0..e {
            acc = acc * b.rem_euclid(p) % p;
        }
        acc
    }

    /// `(j(E), j(E/⟨P⟩))` for a rational point `P` of order `l` on
    /// `y² = x³ + ax + b` over GF(p), by Vélu's formulas, if one exists.
    fn velu_pair(p: u64, a: u64, b: u64, l: u64) -> Option<(i128, i128)> {
        let pi = p as i128;
        let md = |x: i128| x.rem_euclid(pi);
        let inv = |x: i128| pow_mod(md(x), (p - 2) as u32, pi);
        let (a, b) = (a as i128, b as i128);
        let add = |s: Option<(i128, i128)>, t: Option<(i128, i128)>| -> Option<(i128, i128)> {
            let (Some((x1, y1)), Some((x2, y2))) = (s, t) else {
                return s.or(t);
            };
            let lam = if x1 == x2 {
                if md(y1 + y2) == 0 {
                    return None;
                }
                md((3 * x1 * x1 + a) * inv(2 * y1))
            } else {
                md((y2 - y1) * inv(x2 - x1))
            };
            let x3 = md(lam * lam - x1 - x2);
            Some((x3, md(lam * (x1 - x3) - y1)))
        };
        let n = super::super::certify::count_prime(p, a as u64, b as u64);
        if n % l != 0 {
            return None;
        }
        let j = |a: i128, b: i128| {
            let a3 = md(4 * a * a * a);
            let den = md(a3 + 27 * b * b);
            (den != 0).then(|| md(1728 * a3 * inv(den)))
        };
        for x in 0..pi {
            let rhs = md(x * x * x + a * x + b);
            let Some(y) = (0..pi).find(|y| md(y * y) == rhs) else {
                continue;
            };
            let mut acc = None;
            for _ in 0..(n / l) {
                acc = add(acc, Some((x, y)));
            }
            let Some(gen) = acc else {
                continue;
            };
            // Vélu over {P, 2P, …, ((l−1)/2)P}.
            let (mut t, mut w) = (0i128, 0i128);
            let mut q = Some(gen);
            for _ in 0..(l - 1) / 2 {
                let (xq, yq) = q?;
                let tq = md(2 * (3 * xq * xq + a));
                t = md(t + tq);
                w = md(w + 4 * yq * yq + xq * tq);
                q = add(q, Some(gen));
            }
            return Some((j(a, b)?, j(md(a - 5 * t), md(b - 7 * w))?));
        }
        None
    }

    /// Every embedded `Φ_l` vanishes on `l`-isogenous `j`-invariants:
    /// codomains by Vélu's formulas over GF(p).
    #[test]
    fn modular_polynomials_vanish_on_isogenies() {
        for l in [3u64, 5, 7, 11, 13, 17, 19] {
            let phi = modular_polynomial(l);
            assert!(phi
                .iter()
                .any(|(i, j, c)| *i == l as usize + 1 && *j == 0 && c.is_one()));
            let mut checked = 0;
            'curves: for p in [1009u64, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061] {
                for a in 1..30u64 {
                    for b in 1..30u64 {
                        if (4 * a * a * a + 27 * b * b) % p == 0 {
                            continue;
                        }
                        let Some((j1, j2)) = velu_pair(p, a, b, l) else {
                            continue;
                        };
                        let pi = p as i128;
                        let mut sum = 0i128;
                        for (i, k, c) in &phi {
                            let c = c.mod_floor(&BigInt::from(p)).to_i128().expect("small");
                            sum = (sum + c * pow_mod(j1, *i as u32, pi) % pi * pow_mod(j2, *k as u32, pi))
                                % pi;
                        }
                        assert_eq!(sum, 0, "Φ_{l} at p={p} a={a} b={b}");
                        checked += 1;
                        if checked >= 8 {
                            break 'curves;
                        }
                    }
                }
            }
            assert!(checked >= 8, "Φ_{l}: only {checked} isogenies");
        }
    }

    /// `(d_K, v)` with `Δ = v²·d_K`, `d_K` fundamental.
    fn fundamental(delta: i64) -> (i64, i64) {
        let (mut dk, mut v) = (delta, 1i64);
        let mut s = 3i64;
        while s * s <= dk.abs() {
            while dk % (s * s) == 0 {
                dk /= s * s;
                v *= s;
            }
            s += 2;
        }
        while dk % 4 == 0 && (dk / 4).rem_euclid(4) <= 1 && dk / 4 != 0 {
            dk /= 4;
            v *= 2;
        }
        (dk, v)
    }

    /// The walk closes on curves over small prime fields, and the bound
    /// holds: the conductor of `End(E)`, from the exact scalar tests on
    /// `E[ℓ^j]` at every `ℓ | v`, divides one of the multiples.
    #[test]
    fn cycles_bound_the_conductor() {
        let mut closed = 0;
        let mut sharp = 0;
        for p in [1009u64, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 2003, 3001, 4001] {
            for a in 1..10u64 {
                'curve: for b in 1..10u64 {
                    if (4 * a * a * a + 27 * b * b) % p == 0 {
                        continue;
                    }
                    let n = super::super::certify::count_prime(p, a, b) as i64;
                    let t = p as i64 + 1 - n;
                    let (dk, v) = fundamental(t * t - 4 * p as i64);
                    let model = Model::Prime {
                        p: BigUint::from(p),
                        a: BigUint::from(a),
                        b: BigUint::from(b),
                    };
                    let tt = BigInt::from(t);
                    // The conductor from the scalar tests.
                    let mut f = 1u64;
                    let mut rest = v as u64;
                    let mut ell = 2u64;
                    while rest > 1 {
                        if rest % ell != 0 {
                            ell += 1;
                            continue;
                        }
                        let mut depth = 0u32;
                        while rest % ell == 0 {
                            rest /= ell;
                            depth += 1;
                        }
                        let mut s = 0u32;
                        while s < depth {
                            let scalar = if s == 0 {
                                frobenius_is_scalar(&model, &tt, ell)
                            } else {
                                scalar_on_power(&model, &tt, ell, s + 1)
                            };
                            match scalar {
                                Some(true) => s += 1,
                                Some(false) => break,
                                None => continue 'curve,
                            }
                        }
                        f *= ell.pow(depth - s);
                        ell += 1;
                    }
                    let cycles = walks(&model, &BigInt::from(dk), &BigUint::from(v as u64)).closed;
                    if cycles.is_empty() {
                        continue;
                    }
                    closed += 1;
                    // The cycle length is the order of a class of the
                    // class group of End(E): it divides h(O_f).
                    let mut ff = Vec::new();
                    let mut rest = f;
                    let mut d = 2u64;
                    while rest > 1 {
                        let mut e = 0;
                        while rest % d == 0 {
                            rest /= d;
                            e += 1;
                        }
                        if e > 0 {
                            ff.push((BigUint::from(d), e));
                        }
                        d += 1;
                    }
                    let hk = class_number(&BigInt::from(dk)).expect("small");
                    let h = order_class_number(&BigInt::from(dk), hk, &ff);
                    for cycle in &cycles {
                        assert!(
                            (&h % cycle.steps).is_zero(),
                            "p={p} a={a} b={b}: cycle of {} {}-isogenies, h(O_{f}) = {h}",
                            cycle.steps,
                            cycle.degree
                        );
                    }
                    for cycle in &cycles {
                    assert!(
                        cycle.multiples.iter().any(|y| (y % f).is_zero()),
                        "p={p} a={a} b={b}: conductor {f} divides none of {:?}",
                        cycle.multiples
                    );
                    }
                    // Sharp when the bounds together pin the conductor.
                    let joint = cycles.iter().fold(BigUint::from(v as u64), |acc, c| {
                        acc.gcd(c.multiples.iter().max().expect("non-empty"))
                    });
                    if joint == BigUint::from(f) {
                        sharp += 1;
                    }
                }
            }
        }
        assert!(closed > 100 && sharp > 50, "closed {closed}, sharp {sharp}");
    }
}
