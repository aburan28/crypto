//! PKM tower systems over any prime-field backend, for `pkm_tower_pilot`
//! at primes its `u64` builder cannot hold (NIST P-192, P-224, P-256), and
//! at small primes as the builder's own cross-check.
//!
//! Everything here is generic over [`PrimeField`]: the factor-base set `V`,
//! the target curve, the targets, `S₃` and its chain, the tower's square
//! rules, the solve ([`f4_tower_generic`]) and the verification. Random
//! draws are made on canonical integers (`BigUint`), never on a backend's
//! representation, so one seed builds the same system under every backend
//! that holds the prime, and the backends can be compared row for row.
//!
//! Two tower kinds:
//! - `kummer`: `y_{j+1} = y_j²`, `V = g·μ_{2^t}`, needs `2^t | p − 1`
//!   (P-224: `v₂(p − 1) = 96`);
//! - `dickson-torus`: `y_{j+1} = y_j² − 2`, `V = {w + w⁻¹ : w ∈ u·μ_{2^t}}`
//!   with `u`, `μ_{2^t}` in the norm-1 torus of `F_{p²}`, needs
//!   `2^t | p + 1` and `p ≡ 3 mod 4` (P-192: `v₂(p + 1) = 64`; P-256: 96).
//!
//! **Verification.** After each solve, every decomposition of the target
//! over `V^m` is enumerated independently of the engine: the `S₃` links are
//! solved as quadratics in their last argument (backend square roots) and
//! membership in `V` is a hash lookup. The engine's verdict must agree (the
//! unit ideal exactly when there is none). Every decomposition found, as a
//! full point (all tower levels and chain unknowns), is then re-checked in
//! `BigUint` arithmetic against every input equation, every tower
//! equation and every element of the returned basis; the planted point,
//! when there is one, likewise.

use std::collections::{BTreeMap, HashSet};
use std::time::{Duration, Instant};

use crypto_lib::cryptanalysis::f4_fp_tower::{self, StepTrace, TowerF4Options};
use crypto_lib::cryptanalysis::f4_tower_generic::{
    f4_tower_generic, staircase, GPoly, GSquareRule, GTowerRing,
};
use crypto_lib::cryptanalysis::fp_backend::PrimeField;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::Rng;
use serde_json::{json, Value};

// ── Canonical random draws ─────────────────────────────────────────

/// A uniform integer in `[0, p)`, drawn 32 bits at a time with rejection.
pub fn rand_below(p: &BigUint, rng: &mut StdRng) -> BigUint {
    let bits = p.bits();
    let words = bits.div_ceil(32) as usize;
    let top = (bits % 32) as u32;
    loop {
        let mut d: Vec<u32> = (0..words).map(|_| rng.gen::<u32>()).collect();
        if top != 0 {
            d[words - 1] &= (1u32 << top) - 1;
        }
        let x = BigUint::new(d);
        if &x < p {
            return x;
        }
    }
}

/// A uniform integer in `[lo, p)`.
fn rand_from(lo: u64, p: &BigUint, rng: &mut StdRng) -> BigUint {
    loop {
        let x = rand_below(p, rng);
        if x >= BigUint::from(lo) {
            return x;
        }
    }
}

fn elem<F: PrimeField>(f: &F, lo: u64, rng: &mut StdRng) -> F::Elem {
    f.from_biguint(&rand_from(lo, f.modulus(), rng))
}

fn s<F: PrimeField>(f: &F, x: F::Elem) -> String {
    f.to_biguint(x).to_string()
}

// ── Polynomials while building a system ────────────────────────────

/// A polynomial in `n` variables: exponent vectors to non-zero coefficients.
#[derive(Clone, Debug)]
pub struct WPol<E> {
    n: usize,
    terms: BTreeMap<Vec<u32>, E>,
}

impl<E: Copy + PartialEq> WPol<E> {
    fn zero(n: usize) -> Self {
        WPol {
            n,
            terms: BTreeMap::new(),
        }
    }
    fn constant<F: PrimeField<Elem = E>>(f: &F, n: usize, c: E) -> Self {
        let mut p = WPol::zero(n);
        if !f.is_zero(c) {
            p.terms.insert(vec![0; n], c);
        }
        p
    }
    fn var<F: PrimeField<Elem = E>>(f: &F, n: usize, i: usize) -> Self {
        let mut e = vec![0; n];
        e[i] = 1;
        let mut p = WPol::zero(n);
        p.terms.insert(e, f.one());
        p
    }
    fn add<F: PrimeField<Elem = E>>(&self, f: &F, o: &Self) -> Self {
        let mut r = self.clone();
        for (e, &c) in &o.terms {
            let v = r.terms.entry(e.clone()).or_insert(f.zero());
            *v = f.add(*v, c);
        }
        r.terms.retain(|_, c| !f.is_zero(*c));
        r
    }
    fn scale<F: PrimeField<Elem = E>>(&self, f: &F, k: E) -> Self {
        let mut r = WPol::zero(self.n);
        for (e, &c) in &self.terms {
            let v = f.mul(c, k);
            if !f.is_zero(v) {
                r.terms.insert(e.clone(), v);
            }
        }
        r
    }
    fn sub<F: PrimeField<Elem = E>>(&self, f: &F, o: &Self) -> Self {
        self.add(f, &o.scale(f, f.neg(f.one())))
    }
    fn mul<F: PrimeField<Elem = E>>(&self, f: &F, o: &Self) -> Self {
        let mut r = WPol::zero(self.n);
        for (e1, &c1) in &self.terms {
            for (e2, &c2) in &o.terms {
                let e: Vec<u32> = e1.iter().zip(e2).map(|(a, b)| a + b).collect();
                let v = r.terms.entry(e).or_insert(f.zero());
                *v = f.add(*v, f.mul(c1, c2));
            }
        }
        r.terms.retain(|_, c| !f.is_zero(*c));
        r
    }
    fn eval<F: PrimeField<Elem = E>>(&self, f: &F, x: &[E]) -> E {
        let mut acc = f.zero();
        for (e, &c) in &self.terms {
            let mut t = c;
            for (i, &k) in e.iter().enumerate() {
                if k > 0 {
                    t = f.mul(t, f.pow_u64(x[i], u64::from(k)));
                }
            }
            acc = f.add(acc, t);
        }
        acc
    }
    pub fn degree(&self) -> u32 {
        self.terms
            .keys()
            .map(|e| e.iter().sum::<u32>())
            .max()
            .unwrap_or(0)
    }
}

/// `p` at `x`, in `BigUint` arithmetic modulo `m`.
fn eval_big(terms: &[(Vec<u32>, BigUint)], x: &[BigUint], m: &BigUint) -> BigUint {
    let mut acc = BigUint::zero();
    for (e, c) in terms {
        let mut t = c.clone();
        for (i, &k) in e.iter().enumerate() {
            if k > 0 {
                t = t * x[i].modpow(&BigUint::from(k), m) % m;
            }
        }
        acc = (acc + t) % m;
    }
    acc
}

/// Semaev's third summation polynomial for `y² = x³ + ax + b`:
/// `(X1 − X2)²X3² − 2((X1 + X2)(X1X2 + a) + 2b)X3 + (X1X2 − a)² − 4b(X1 + X2)`.
fn s3<F: PrimeField>(
    f: &F,
    x1: &WPol<F::Elem>,
    x2: &WPol<F::Elem>,
    x3: &WPol<F::Elem>,
    a: F::Elem,
    b: F::Elem,
) -> WPol<F::Elem> {
    let n = x1.n;
    let d = x1.sub(f, x2);
    let t1 = d.mul(f, &d).mul(f, x3).mul(f, x3);
    let sm = x1.add(f, x2);
    let pr = x1.mul(f, x2);
    let inner = sm
        .mul(f, &pr.add(f, &WPol::constant(f, n, a)))
        .add(f, &WPol::constant(f, n, f.double(b)));
    let t2 = inner.mul(f, x3).scale(f, f.neg(f.from_u64(2)));
    let q = pr.sub(f, &WPol::constant(f, n, a));
    let t3 = q
        .mul(f, &q)
        .sub(f, &sm.scale(f, f.mul(f.from_u64(4), b)));
    t1.add(f, &t2).add(f, &t3)
}

/// The roots `z` of `S₃(x1, x2, z)`, a quadratic in `z`.
fn s3_roots<F: PrimeField>(f: &F, x1: F::Elem, x2: F::Elem, a: F::Elem, b: F::Elem) -> Vec<F::Elem> {
    let d = f.sub(x1, x2);
    let qa = f.sqr(d);
    let sm = f.add(x1, x2);
    let pr = f.mul(x1, x2);
    let qb = f.neg(f.double(f.add(f.mul(sm, f.add(pr, a)), f.double(b))));
    let q = f.sub(pr, a);
    let qc = f.sub(f.sqr(q), f.mul(f.mul(f.from_u64(4), b), sm));
    if f.is_zero(qa) {
        if f.is_zero(qb) {
            assert!(
                !f.is_zero(qc),
                "S3 vanishes identically in its last argument"
            );
            return Vec::new();
        }
        return vec![f.mul(f.neg(qc), f.inv(qb).unwrap())];
    }
    let disc = f.sub(f.sqr(qb), f.mul(f.mul(f.from_u64(4), qa), qc));
    let Some(r) = f.sqrt(disc) else {
        return Vec::new();
    };
    let inv2a = f.inv(f.double(qa)).unwrap();
    let z1 = f.mul(f.sub(r, qb), inv2a);
    let z2 = f.mul(f.sub(f.neg(r), qb), inv2a);
    if z1 == z2 {
        vec![z1]
    } else {
        vec![z1, z2]
    }
}

// ── Towers ─────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WKind {
    Kummer,
    DicksonTorus,
}

impl WKind {
    pub fn name(self) -> &'static str {
        match self {
            WKind::Kummer => "kummer",
            WKind::DicksonTorus => "dickson-torus",
        }
    }
    /// The pilot's `Kind` discriminant, so cell seeds follow its formula.
    pub fn seed_tag(self) -> u64 {
        match self {
            WKind::Kummer => 0,
            WKind::DicksonTorus => 3,
        }
    }
}

pub struct WTower<E> {
    kind: WKind,
    t: usize,
    v: Vec<E>,
    scale: E,
    c: E,
    info: Value,
}

impl<E: Copy + PartialEq + std::hash::Hash + Eq> WTower<E> {
    fn dickson(&self) -> bool {
        self.kind == WKind::DicksonTorus
    }
    fn step<F: PrimeField<Elem = E>>(&self, f: &F, y: E) -> E {
        if self.dickson() {
            f.sub(f.sqr(y), f.from_u64(2))
        } else {
            f.sqr(y)
        }
    }
    fn values<F: PrimeField<Elem = E>>(&self, f: &F, x: E) -> Vec<E> {
        let mut y = vec![f.zero(); self.t];
        y[0] = f.mul(x, f.inv(self.scale).unwrap());
        for j in 0..self.t - 1 {
            y[j + 1] = self.step(f, y[j]);
        }
        y
    }
    fn equations<F: PrimeField<Elem = E>>(&self, f: &F, n: usize, base: usize) -> Vec<WPol<E>> {
        (0..self.t)
            .map(|j| {
                let yj = WPol::var(f, n, base + j);
                let next = if j + 1 < self.t {
                    WPol::var(f, n, base + j + 1)
                } else {
                    WPol::constant(f, n, self.c)
                };
                let mut sq = yj.mul(f, &yj);
                if self.dickson() {
                    sq = sq.sub(f, &WPol::constant(f, n, f.from_u64(2)));
                }
                next.sub(f, &sq)
            })
            .collect()
    }
    fn square_rules<F: PrimeField<Elem = E>>(&self, f: &F, base: usize) -> Vec<GSquareRule<E>> {
        let two = if self.dickson() {
            f.from_u64(2)
        } else {
            f.zero()
        };
        (0..self.t)
            .map(|j| {
                let top = j + 1 == self.t;
                GSquareRule {
                    next: (!top).then_some((base + j + 1) as u8),
                    a: f.zero(),
                    b: if top { f.zero() } else { f.one() },
                    c: f.zero(),
                    d: if top { f.add(self.c, two) } else { two },
                }
            })
            .collect()
    }
    /// `|V| = 2^t`, distinct, and every level equation holds on `V`.
    fn check<F: PrimeField<Elem = E>>(&self, f: &F) {
        let set: HashSet<E> = self.v.iter().copied().collect();
        assert_eq!(set.len(), 1 << self.t, "{:?}: V has repeats", self.kind);
        let eqs = self.equations(f, self.t, 0);
        for &x in &self.v {
            let y = self.values(f, x);
            for (j, e) in eqs.iter().enumerate() {
                assert!(
                    f.is_zero(e.eval(f, &y)),
                    "{:?}: level {j} fails",
                    self.kind
                );
            }
        }
    }
}

fn kummer<F: PrimeField>(f: &F, t: usize, rng: &mut StdRng) -> WTower<F::Elem> {
    let p = f.modulus();
    let pm1 = p - 1u32;
    assert!(
        pm1.trailing_zeros().unwrap_or(0) >= t as u64,
        "2^{t} does not divide p − 1"
    );
    let zeta = loop {
        let h = elem(f, 2, rng);
        if f.legendre(h) == -1 {
            break f.pow_biguint(h, &(&pm1 >> t));
        }
    };
    assert!(f.is_one(f.pow_u64(zeta, 1 << t)));
    assert!(!f.is_one(f.pow_u64(zeta, 1 << (t - 1))), "ζ of exact order 2^t");
    let g = elem(f, 2, rng);
    let mut v = Vec::with_capacity(1 << t);
    let mut w = g;
    for _ in 0..1u64 << t {
        v.push(w);
        w = f.mul(w, zeta);
    }
    WTower {
        kind: WKind::Kummer,
        t,
        v,
        scale: g,
        c: f.one(),
        info: json!({"g": s(f, g), "zeta": s(f, zeta)}),
    }
}

// `F_{p²} = F_p(i)`, `i² = −1` (`p ≡ 3 mod 4`).
type Fp2<E> = (E, E);

fn fp2_mul<F: PrimeField>(f: &F, x: Fp2<F::Elem>, y: Fp2<F::Elem>) -> Fp2<F::Elem> {
    (
        f.sub(f.mul(x.0, y.0), f.mul(x.1, y.1)),
        f.add(f.mul(x.0, y.1), f.mul(x.1, y.0)),
    )
}

fn fp2_pow<F: PrimeField>(f: &F, x: Fp2<F::Elem>, e: &BigUint) -> Fp2<F::Elem> {
    let mut r = (f.one(), f.zero());
    for i in (0..e.bits()).rev() {
        r = fp2_mul(f, r, r);
        if e.bit(i) {
            r = fp2_mul(f, r, x);
        }
    }
    r
}

fn fp2_inv<F: PrimeField>(f: &F, x: Fp2<F::Elem>) -> Fp2<F::Elem> {
    let n = f.add(f.sqr(x.0), f.sqr(x.1));
    let ni = f.inv(n).unwrap();
    (f.mul(x.0, ni), f.mul(f.neg(x.1), ni))
}

fn torus_element<F: PrimeField>(f: &F, rng: &mut StdRng) -> Fp2<F::Elem> {
    let p = f.modulus();
    loop {
        let z = (
            f.from_biguint(&rand_below(p, rng)),
            f.from_biguint(&rand_below(p, rng)),
        );
        if f.is_zero(z.0) && f.is_zero(z.1) {
            continue;
        }
        return fp2_pow(f, z, &(p - 1u32));
    }
}

fn dickson_torus<F: PrimeField>(f: &F, t: usize, rng: &mut StdRng) -> WTower<F::Elem> {
    let p = f.modulus().clone();
    let pp1 = &p + 1u32;
    assert!(
        (&p % 4u32) == BigUint::from(3u32) && pp1.trailing_zeros().unwrap_or(0) >= t as u64,
        "the torus case needs 2^{t} | p + 1 and p ≡ 3 mod 4"
    );
    let one = (f.one(), f.zero());
    let minus_one = (f.neg(f.one()), f.zero());
    let zeta = loop {
        let h = torus_element(f, rng);
        if fp2_pow(f, h, &(&pp1 >> 1)) == minus_one {
            break fp2_pow(f, h, &(&pp1 >> t));
        }
    };
    let two_t = BigUint::one() << t;
    loop {
        let u = torus_element(f, rng);
        if f.is_zero(u.1) || fp2_pow(f, u, &(&two_t << 1)) == one {
            continue;
        }
        let mut v = Vec::with_capacity(1 << t);
        let mut w = u;
        for _ in 0..1u64 << t {
            v.push(f.double(w.0));
            w = fp2_mul(f, w, zeta);
        }
        let set: HashSet<F::Elem> = v.iter().copied().collect();
        if set.len() != 1 << t {
            continue;
        }
        let ut = fp2_pow(f, u, &two_t);
        let c = f.double(ut.0);
        let tower = WTower {
            kind: WKind::DicksonTorus,
            t,
            v,
            scale: f.one(),
            c,
            info: json!({"case": "torus", "u": [s(f, u.0), s(f, u.1)],
                         "zeta": [s(f, zeta.0), s(f, zeta.1)], "c": s(f, c)}),
        };
        check_torus(f, &tower, u, zeta, &p);
        return tower;
    }
}

/// The pilot's `check_torus` over `F` (note §2.1): `ζ` of exact order
/// `2^t` in the torus, `u` in the torus and not in `F_p`, every `V[k]` is
/// `w + 1/w` by the general inverse, `y² − 4` is a non-residue on `V`, no
/// level is ramified, and `c ≠ ±2`.
fn check_torus<F: PrimeField>(
    f: &F,
    tower: &WTower<F::Elem>,
    u: Fp2<F::Elem>,
    zeta: Fp2<F::Elem>,
    p: &BigUint,
) {
    let t = tower.t;
    let one = (f.one(), f.zero());
    let minus_one = (f.neg(f.one()), f.zero());
    assert_eq!(fp2_pow(f, zeta, &(BigUint::one() << t)), one, "ζ^(2^t) ≠ 1");
    assert_eq!(
        fp2_pow(f, zeta, &(BigUint::one() << (t - 1))),
        minus_one,
        "ζ does not have exact order 2^t"
    );
    assert_eq!(fp2_pow(f, u, &(p + 1u32)), one, "u is not in the torus");
    assert!(!f.is_zero(u.1), "u lies in F_p");
    let mut w = u;
    for (k, &x) in tower.v.iter().enumerate() {
        let wi = fp2_inv(f, w);
        let sum = (f.add(w.0, wi.0), f.add(w.1, wi.1));
        assert!(f.is_zero(sum.1), "w + 1/w ∉ F_p at k = {k}");
        assert!(sum.0 == x, "V[{k}] is not w + 1/w");
        let disc = f.sub(f.sqr(x), f.from_u64(4));
        assert_eq!(f.legendre(disc), -1, "not the torus case at V[{k}]");
        let y = tower.values(f, x);
        assert!(y.iter().all(|&yj| !f.is_zero(yj)), "ramified at V[{k}]");
        w = fp2_mul(f, w, zeta);
    }
    let two = f.from_u64(2);
    assert!(tower.c != two && tower.c != f.neg(two), "c = ±2");
}

pub fn build_tower<F: PrimeField>(
    f: &F,
    kind: WKind,
    t: usize,
    rng: &mut StdRng,
) -> WTower<F::Elem> {
    let tower = match kind {
        WKind::Kummer => kummer(f, t, rng),
        WKind::DicksonTorus => dickson_torus(f, t, rng),
    };
    tower.check(f);
    tower
}

// ── Curve and targets ──────────────────────────────────────────────

pub struct WCurve<E> {
    a: E,
    b: E,
}

impl<E: Copy + PartialEq> WCurve<E> {
    fn rhs<F: PrimeField<Elem = E>>(&self, f: &F, x: E) -> E {
        f.add(f.add(f.mul(f.sqr(x), x), f.mul(self.a, x)), self.b)
    }
    /// Affine addition; `None` is the identity.
    fn add<F: PrimeField<Elem = E>>(&self, f: &F, p: (E, E), q: (E, E)) -> Option<(E, E)> {
        let lambda = if p.0 == q.0 {
            if f.is_zero(f.add(p.1, q.1)) {
                return None;
            }
            let num = f.add(f.mul(f.from_u64(3), f.sqr(p.0)), self.a);
            f.mul(num, f.inv(f.double(p.1)).unwrap())
        } else {
            f.mul(f.sub(q.1, p.1), f.inv(f.sub(q.0, p.0)).unwrap())
        };
        let x = f.sub(f.sub(f.sqr(lambda), p.0), q.0);
        let y = f.sub(f.mul(lambda, f.sub(p.0, x)), p.1);
        Some((x, y))
    }
}

pub struct WTarget<E> {
    planted: bool,
    x_r: E,
    solution: Option<Vec<E>>,
}

fn target_curve<F: PrimeField>(
    f: &F,
    tower: &WTower<F::Elem>,
    need: usize,
    rng: &mut StdRng,
) -> (WCurve<F::Elem>, Vec<F::Elem>) {
    loop {
        let a = elem(f, 0, rng);
        let b = elem(f, 1, rng);
        let disc = f.add(
            f.mul(f.from_u64(4), f.mul(f.sqr(a), a)),
            f.mul(f.from_u64(27), f.sqr(b)),
        );
        if f.is_zero(disc) {
            continue;
        }
        let curve = WCurve { a, b };
        let on: Vec<F::Elem> = tower
            .v
            .iter()
            .copied()
            .filter(|&x| f.legendre(curve.rhs(f, x)) == 1)
            .collect();
        if on.len() >= need {
            return (curve, on);
        }
    }
}

fn plant<F: PrimeField>(
    f: &F,
    curve: &WCurve<F::Elem>,
    tower: &WTower<F::Elem>,
    m: usize,
    on: &[F::Elem],
    rng: &mut StdRng,
) -> Option<WTarget<F::Elem>> {
    let mut pick: Vec<F::Elem> = Vec::with_capacity(m);
    while pick.len() < m {
        let x = on[rng.gen_range(0..on.len())];
        if !pick.contains(&x) {
            pick.push(x);
        }
    }
    let pts: Vec<(F::Elem, F::Elem)> = pick
        .iter()
        .map(|&x| {
            let y = f.sqrt(curve.rhs(f, x)).expect("abscissa on the curve");
            (x, if rng.gen::<bool>() { y } else { f.neg(y) })
        })
        .collect();
    let mut acc = pts[0];
    let mut partial = Vec::new();
    for pt in &pts[1..] {
        acc = curve.add(f, acc, *pt)?;
        partial.push(acc.0);
    }
    let t = tower.t;
    let mut sol = vec![f.zero(); m * t + (m - 2)];
    for (i, &x) in pick.iter().enumerate() {
        sol[i * t..(i + 1) * t].copy_from_slice(&tower.values(f, x));
    }
    for k in 0..m - 2 {
        sol[m * t + k] = partial[k];
    }
    Some(WTarget {
        planted: true,
        x_r: acc.0,
        solution: Some(sol),
    })
}

fn random_target<F: PrimeField>(f: &F, curve: &WCurve<F::Elem>, rng: &mut StdRng) -> WTarget<F::Elem> {
    loop {
        let x = elem(f, 0, rng);
        let r = curve.rhs(f, x);
        if f.is_zero(r) || f.legendre(r) != 1 {
            continue;
        }
        return WTarget {
            planted: false,
            x_r: x,
            solution: None,
        };
    }
}

// ── Systems ────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WControl {
    Tower,
    Null,
}

impl WControl {
    pub fn name(self) -> &'static str {
        match self {
            WControl::Tower => "tower",
            WControl::Null => "null",
        }
    }
}

fn summation<F: PrimeField>(
    f: &F,
    tower: &WTower<F::Elem>,
    curve: &WCurve<F::Elem>,
    m: usize,
    n: usize,
    x_r: F::Elem,
) -> Vec<WPol<F::Elem>> {
    let (a, b) = (curve.a, curve.b);
    let x = |i: usize| WPol::var(f, n, i * tower.t).scale(f, tower.scale);
    let xr = WPol::constant(f, n, x_r);
    if m == 2 {
        return vec![s3(f, &x(0), &x(1), &xr, a, b)];
    }
    let u = |k: usize| WPol::var(f, n, m * tower.t + k);
    let mut eqs = vec![s3(f, &x(0), &x(1), &u(0), a, b)];
    for k in 1..m - 2 {
        eqs.push(s3(f, &u(k - 1), &x(k + 1), &u(k), a, b));
    }
    eqs.push(s3(f, &u(m - 3), &x(m - 1), &xr, a, b));
    eqs
}

type System<E> = (Vec<WPol<E>>, usize, Option<Vec<E>>, std::ops::Range<usize>);

fn system<F: PrimeField>(
    f: &F,
    control: WControl,
    tower: &WTower<F::Elem>,
    curve: &WCurve<F::Elem>,
    m: usize,
    target: &WTarget<F::Elem>,
    rng: &mut StdRng,
) -> System<F::Elem> {
    let t = tower.t;
    let n = m * t + (m - 2);
    let mut eqs = summation(f, tower, curve, m, n, target.x_r);
    if control == WControl::Null {
        for eq in eqs.iter_mut() {
            let mut g = WPol::zero(n);
            for e in eq.terms.keys() {
                g.terms.insert(e.clone(), elem(f, 1, rng));
            }
            if let Some(sol) = &target.solution {
                let c = g.eval(f, sol);
                g = g.sub(f, &WPol::constant(f, n, c));
            }
            *eq = g;
        }
    }
    let from = eqs.len();
    for i in 0..m {
        eqs.extend(tower.equations(f, n, i * t));
    }
    let towers = from..eqs.len();
    (eqs, n, target.solution.clone(), towers)
}

// ── Verification ───────────────────────────────────────────────────

/// Every decomposition `x_R = x(P_1 + … + P_m)` over `V^m` as a full point
/// of the tower system (tower levels, then the chain unknowns), enumerated
/// without the engine. Needs the tower control (the null has no curve).
fn decompositions<F: PrimeField>(
    f: &F,
    tower: &WTower<F::Elem>,
    curve: &WCurve<F::Elem>,
    m: usize,
    x_r: F::Elem,
) -> Vec<Vec<F::Elem>> {
    let (a, b) = (curve.a, curve.b);
    let in_v: HashSet<F::Elem> = tower.v.iter().copied().collect();
    // Chains (x_1 … x_k, u_0 … u_{k−2}) with u the last partial sum's x.
    let mut out = Vec::new();
    let point = |xs: &[F::Elem], us: &[F::Elem]| -> Vec<F::Elem> {
        let mut sol = Vec::with_capacity(m * tower.t + us.len());
        for &x in xs {
            sol.extend(tower.values(f, x));
        }
        sol.extend_from_slice(us);
        sol
    };
    if m == 2 {
        for &x1 in &tower.v {
            for z in s3_roots(f, x1, x_r, a, b) {
                if in_v.contains(&z) {
                    out.push(point(&[x1, z], &[]));
                }
            }
        }
        return out;
    }
    let mut partial: Vec<(Vec<F::Elem>, Vec<F::Elem>)> = Vec::new();
    for &x1 in &tower.v {
        for &x2 in &tower.v {
            for u in s3_roots(f, x1, x2, a, b) {
                partial.push((vec![x1, x2], vec![u]));
            }
        }
    }
    for _ in 1..m - 2 {
        let mut next = Vec::new();
        for (xs, us) in &partial {
            let last = *us.last().unwrap();
            for &x in &tower.v {
                for u in s3_roots(f, last, x, a, b) {
                    let mut xs2 = xs.clone();
                    xs2.push(x);
                    let mut us2 = us.clone();
                    us2.push(u);
                    next.push((xs2, us2));
                }
            }
        }
        partial = next;
    }
    for (xs, us) in &partial {
        let last = *us.last().unwrap();
        for z in s3_roots(f, last, x_r, a, b) {
            if in_v.contains(&z) {
                let mut xs2 = xs.clone();
                xs2.push(z);
                out.push(point(&xs2, us));
            }
        }
    }
    out
}

fn big_terms<F: PrimeField>(f: &F, p: &WPol<F::Elem>) -> Vec<(Vec<u32>, BigUint)> {
    p.terms
        .iter()
        .map(|(e, &c)| (e.clone(), f.to_biguint(c)))
        .collect()
}

fn big_basis_terms<F: PrimeField>(
    f: &F,
    g: &GPoly<F::Elem>,
    n_tower: usize,
    n: usize,
) -> Vec<(Vec<u32>, BigUint)> {
    g.monos
        .iter()
        .zip(&g.coefs)
        .map(|(&mo, &c)| {
            let mut e = vec![0u32; n];
            let mask = f4_fp_tower::mask_of(mo);
            for (j, x) in e.iter_mut().enumerate().take(n_tower) {
                *x = ((mask >> j) & 1) as u32;
            }
            let fr = f4_fp_tower::free_of(mo);
            for k in 0..n - n_tower {
                e[n_tower + k] = u32::from(fr[k]);
            }
            (e, f.to_biguint(c))
        })
        .collect()
}

// ── One cell ───────────────────────────────────────────────────────

pub struct WCell {
    pub kind: WKind,
    pub m: usize,
    pub control: WControl,
    pub t: usize,
}

pub struct WArgs {
    pub seed: u64,
    pub planted: usize,
    pub random: usize,
    pub budget: u64,
    pub max_nnz: Option<u64>,
    pub max_dense: Option<u64>,
    pub full_rank_exit: bool,
    pub trace: Option<(fn(usize, &StepTrace), fn(usize, bool, &StepTrace))>,
    pub targets: Option<Vec<usize>>,
    pub verify: bool,
    /// The prime's label (`P-224`, …) for the rows.
    pub label: String,
}

/// Build and solve one cell's systems over `f`; one JSON row per system.
/// Returns the rows and whether any system finished.
pub fn run_cell<F: PrimeField>(f: &F, cell: &WCell, a: &WArgs) -> (Vec<Value>, bool) {
    let cell_seed = a.seed
        ^ (cell.kind.seed_tag() << 40)
        ^ ((cell.m as u64) << 32)
        ^ ((cell.t as u64) << 16)
        ^ 1u64;
    let mut rng = <StdRng as rand::SeedableRng>::seed_from_u64(cell_seed);
    let tower = build_tower(f, cell.kind, cell.t, &mut rng);
    let (curve, on) = target_curve(f, &tower, cell.m.max(2), &mut rng);
    let mut targets = Vec::new();
    while targets.len() < a.planted {
        if let Some(tg) = plant(f, &curve, &tower, cell.m, &on, &mut rng) {
            targets.push(tg);
        }
    }
    for _ in 0..a.random {
        targets.push(random_target(f, &curve, &mut rng));
    }
    let mut rows = Vec::new();
    let mut any_finished = false;
    for (k, target) in targets.iter().enumerate() {
        let (eqs, n, sol, towers) = system(f, cell.control, &tower, &curve, cell.m, target, &mut rng);
        if a.targets.as_ref().is_some_and(|ts| !ts.contains(&k)) {
            continue;
        }
        let d_in: Vec<u32> = eqs.iter().map(WPol::degree).collect();
        let d_max = d_in.iter().copied().max().unwrap_or(1);
        let max_degree = n as u32 + d_max + 6;
        if let Some(s) = &sol {
            for (i, e) in eqs.iter().enumerate() {
                assert!(f.is_zero(e.eval(f, s)), "equation {i} fails at the plant");
            }
        }
        let mut rules = Vec::with_capacity(cell.m * cell.t);
        for i in 0..cell.m {
            rules.extend(tower.square_rules(f, i * cell.t));
        }
        let ring = GTowerRing {
            field: f.clone(),
            rules,
            n_free: cell.m - 2,
        };
        let raw = |e: &WPol<F::Elem>| -> Vec<(Vec<u32>, F::Elem)> {
            e.terms.iter().map(|(k, v)| (k.clone(), *v)).collect()
        };
        for i in towers.clone() {
            assert!(
                ring.from_raw(&raw(&eqs[i])).is_empty(),
                "tower equation {i} is not zero in the ring"
            );
        }
        let input: Vec<GPoly<F::Elem>> = eqs
            .iter()
            .enumerate()
            .filter(|(i, _)| !towers.contains(i))
            .map(|(_, e)| ring.from_raw(&raw(e)))
            .filter(|g| !g.is_empty())
            .collect();
        let ring_degrees: Vec<u32> = input.iter().map(GPoly::degree).collect();
        let mut opts =
            TowerF4Options::new(max_degree).with_budget(Duration::from_secs(a.budget));
        if let Some(c) = a.max_nnz {
            opts = opts.with_max_nnz(c);
        }
        if let Some(c) = a.max_dense {
            opts = opts.with_max_dense(c);
        }
        opts.full_rank_exit = a.full_rank_exit;
        if let Some((step, stop)) = a.trace {
            opts = opts.on_step(step);
            opts.on_stop = Some(stop);
        }
        let started = Instant::now();
        let rep = f4_tower_generic(&input, &ring, &opts);
        let wall = started.elapsed().as_secs_f64();
        let r = &rep.stats;
        let finished = !(r.timed_out || r.oversize);
        any_finished |= finished;
        let planted_ok = sol.as_ref().map(|s| {
            !r.inconsistent && rep.basis.iter().all(|g| f.is_zero(ring.eval(g, s)))
        });
        let lms: Vec<u128> = rep.basis.iter().filter_map(GPoly::lm).collect();
        let pure = |k: usize| {
            lms.iter().any(|&mo| {
                let fr = f4_fp_tower::free_of(mo);
                f4_fp_tower::mask_of(mo) == 0
                    && (0..f4_fp_tower::MAX_FREE).all(|i| (i == k) == (fr[i] > 0))
            })
        };
        let zero_dim = r.inconsistent || (0..ring.n_free).all(pure);

        // Verification, independent of the engine, and BigUint re-checks.
        let verify = (a.verify && finished && cell.control == WControl::Tower).then(|| {
            let vstart = Instant::now();
            let m_big = f.modulus().clone();
            let sols = decompositions(f, &tower, &curve, cell.m, target.x_r);
            let eq_big: Vec<_> = eqs.iter().map(|e| big_terms(f, e)).collect();
            let basis_big: Vec<_> = rep
                .basis
                .iter()
                .map(|g| big_basis_terms(f, g, ring.n_tower(), n))
                .collect();
            let mut checked = 0usize;
            let mut all_ok = true;
            let mut points: Vec<Vec<F::Elem>> = sols.clone();
            if let Some(s) = &sol {
                points.push(s.clone());
            }
            for pt in &points {
                let x: Vec<BigUint> = pt.iter().map(|&v| f.to_biguint(v)).collect();
                let eqs_ok = eq_big.iter().all(|e| eval_big(e, &x, &m_big).is_zero());
                let basis_ok =
                    !r.inconsistent && basis_big.iter().all(|g| eval_big(g, &x, &m_big).is_zero());
                all_ok &= eqs_ok && basis_ok;
                checked += 1;
            }
            let stair = (!r.inconsistent && zero_dim)
                .then(|| {
                    staircase(&lms, ring.n_tower(), ring.n_free, 1 << 20)
                })
                .flatten();
            let distinct: HashSet<Vec<String>> = sols
                .iter()
                .map(|p| p.iter().map(|&v| s(f, v)).collect())
                .collect();
            let verdict_ok = r.inconsistent == distinct.is_empty();
            let count_ok = stair.is_none_or(|c| c >= distinct.len());
            json!({
                "decompositions": distinct.len(),
                "points_checked_biguint": checked,
                "points_ok": all_ok,
                "verdict_agrees": verdict_ok,
                "staircase": stair,
                "staircase_covers": count_ok,
                "ok": all_ok && verdict_ok && count_ok,
                "seconds": vstart.elapsed().as_secs_f64(),
            })
        });
        let row = json!({
            "kind": cell.kind.name(),
            "m": cell.m,
            "control": cell.control.name(),
            "t": cell.t,
            "N": cell.m * cell.t,
            "g": 1,
            "n_vars": n,
            "n_equations": eqs.len(),
            "input_degrees": d_in,
            "target": if target.planted { "planted" } else { "random" },
            "target_index": k,
            "x_r": s(f, target.x_r),
            "p": f.modulus().to_string(),
            "prime": a.label,
            "p_bits": f.modulus().bits(),
            "v2_p_minus_1": (f.modulus() - 1u32).trailing_zeros(),
            "v2_p_plus_1": (f.modulus() + 1u32).trailing_zeros(),
            "curve": {"a": s(f, curve.a), "b": s(f, curve.b)},
            "tower": tower.info,
            "max_degree_bound": max_degree,
            "solving_degree_max": r.solving_degree_max,
            "last_productive_degree": r.last_productive_degree,
            "degree_reached": r.degree_reached,
            "max_cols_to_solution": r.max_cols_to_solution,
            "steps_to_solution": r.steps_to_solution,
            "inconsistent": r.inconsistent,
            "zero_dim": zero_dim,
            "basis_len": rep.basis.len(),
            "steps": r.steps,
            "max_rows": r.max_rows,
            "max_cols": r.max_cols,
            "ms": r.ms,
            "wall_s": wall,
            "timed_out": r.timed_out || r.oversize,
            "oversize": r.oversize,
            "pairs_above_bound": r.pairs_above_bound,
            "planted_ok": planted_ok,
            "engine": "f4_tower_generic",
            "builder": "wide",
            "backend": f.name(),
            "input_ring_degrees": ring_degrees,
            "muladds": r.muladds,
            "max_nnz": r.max_nnz,
            "max_residual_rows": r.max_residual_rows,
            "max_dense_entries": r.max_dense_entries,
            "full_rank_exit": a.full_rank_exit,
            "rows_skipped_full_rank": r.rows_skipped,
            "critical_pairs_reduced": r.critical_pairs_reduced,
            "tower_pairs_reduced": r.tower_pairs_reduced,
            "reducer_rows": r.reducer_rows,
            "stopped_step": r.stopped_step,
            "verify": verify,
        });
        eprintln!(
            "{:13} {} m={} {:5} t={:2} N={:2} {:7} D={:2} cols={:7} steps={:3} muladds={:>15} {:9.1} ms{}{}",
            cell.kind.name(),
            a.label,
            cell.m,
            cell.control.name(),
            cell.t,
            cell.m * cell.t,
            if target.planted { "planted" } else { "random" },
            r.solving_degree_max,
            r.max_cols_to_solution,
            r.steps_to_solution,
            r.muladds,
            r.ms,
            if finished { "" } else { "  STOPPED" },
            match row["verify"]["ok"].as_bool() {
                Some(true) => "  verify ok",
                Some(false) => "  VERIFY-FAIL",
                None => "",
            },
        );
        rows.push(row);
    }
    (rows, any_finished)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_lib::cryptanalysis::fp_backend::{Fp224, Fp32, Fp64, FpWide};
    use rand::SeedableRng;

    fn p224() -> BigUint {
        (BigUint::one() << 224u32) - (BigUint::one() << 96u32) + 1u32
    }
    fn p192() -> BigUint {
        (BigUint::one() << 192u32) - (BigUint::one() << 64u32) - 1u32
    }

    #[test]
    fn towers_check_at_nist_primes() {
        let mut rng = StdRng::seed_from_u64(1);
        for t in 1..=8 {
            let _ = build_tower(&Fp224::new(&p224()), WKind::Kummer, t, &mut rng);
            let _ = build_tower(&FpWide::<3>::new(&p192()), WKind::DicksonTorus, t, &mut rng);
        }
    }

    /// One seed builds the same system under every backend holding the
    /// prime, and the engine gives the same row under each.
    #[test]
    fn rows_agree_across_backends() {
        let p = 3_221_225_473u64;
        let a = WArgs {
            seed: 7,
            planted: 1,
            random: 1,
            budget: 60,
            max_nnz: None,
            max_dense: None,
            full_rank_exit: true,
            trace: None,
            targets: None,
            verify: true,
            label: "p2".into(),
        };
        let cell = WCell {
            kind: WKind::Kummer,
            m: 2,
            control: WControl::Tower,
            t: 3,
        };
        let strip = |rows: Vec<Value>| -> Vec<Value> {
            rows.into_iter()
                .map(|mut r| {
                    let o = r.as_object_mut().unwrap();
                    for k in ["ms", "wall_s", "backend"] {
                        o.remove(k);
                    }
                    o["verify"].as_object_mut().unwrap().remove("seconds");
                    r
                })
                .collect()
        };
        let want = strip(run_cell(&Fp32::new(p), &cell, &a).0);
        assert!(want.iter().all(|r| r["verify"]["ok"] == json!(true)));
        assert_eq!(strip(run_cell(&Fp64::new(p), &cell, &a).0), want);
        assert_eq!(strip(run_cell(&FpWide::<4>::new(&BigUint::from(p)), &cell, &a).0), want);
    }
}
