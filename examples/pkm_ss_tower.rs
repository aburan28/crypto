//! **Supersingular auxiliary curves for the PKM isogeny tower**
//! (`research/pkm_supersingular_aux_20261015/`).
//!
//! Two subcommands:
//!
//! - `structure`: Q1, an exact computation. At small primes `p ≡ 3 mod 4`
//!   with `2^t | p + 1` it lists *every* tower set of each family:
//!   - `torus`: the Dickson torus sets `{w + w^{−1} : w ∈ u·μ_{2^t}}`,
//!     `u` in the norm-1 torus of `F_{p²}`, `u² ∉ μ_{2^t}`, `|V| = 2^t`;
//!   - `ss-j1728`, `ss-j1728-twist`, `ss-j0`, `ss-j0-twist`: every
//!     `x(R' + ⟨Q⟩)` with `⟨Q⟩` a cyclic subgroup of order `2^t` and
//!     `2R' ∉ ⟨Q⟩`, on the supersingular curves of
//!     `prime_tower::supersingular_aux`;
//!   - `ordinary`: the same on a few ordinary curves with a point of order
//!     `2^t` (the control), found by exhaustive point counting;
//!
//!   and decides, for every pair of families, which sets of one are equal
//!   to a set of the other up to a Möbius map `x ↦ (αx + β)/(γx + δ)` over
//!   `F_p` (affine maps included). The test is exact: a set's canonical
//!   form is the lexicographically least image of the set under the `n(n −
//!   1)(n − 2)` Möbius maps sending an ordered triple of its points to
//!   `(0, 1, ∞)`; two sets are `PGL₂(F_p)`-equivalent iff their forms are
//!   equal. A brute-force check over all of `PGL₂(F_p)` confirms the forms
//!   at the smallest primes (`--brute`).
//!
//! - `f4`: Q2, the F4 solving degree of PKM's `m = 2` (`full`) and `m ≥ 3`
//!   (`chain`) systems with the tower from one arm (`torus`, `split`,
//!   `ordinary`, `ss-j1728`, `ss-j1728-twist`, `ss-j0`, `ss-j0-twist`) and
//!   the `tower` or shape-matched `null` control, on the library engine
//!   `f4_fp_tower` (reduced presentation, full-rank exit, degree bound
//!   `n + d + 6`, as round 8). Every finished row is re-checked
//!   exhaustively over `V^m` (`verify` field).
//!
//! The system builder (`Pol`, `s3`, the tower equations and square rules,
//! the summation and the null control, the Dickson towers) is copied from
//! `examples/pkm_tower_pilot.rs`, whose builders are private to that
//! example, with the minimum changes: the isogeny towers come from the
//! library (`prime_tower::{isogeny_tower, supersingular_tower}`), and the
//! seeds of the cells are this example's own. It is a stage diagnostic of
//! the solver axis (`AGENTS.md` §2): it prices no pipeline and reports no
//! `S`.
//!
//! ```bash
//! cargo run --release --example pkm_ss_tower -- structure --p-max 4000 --t 2,3,4 --out q1.jsonl
//! cargo run --release --example pkm_ss_tower -- f4 --p 14680063 --arm ss-j1728 --m 2 \
//!     --t-min 1 --t-max 9 --targets 2 --out r.jsonl
//! ```

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::time::{Duration, Instant};

use crypto_lib::cryptanalysis::f4_fp_tower::{self, SquareRule, TowerF4Options, TowerRing};
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps, PrimeCurve, PrimePoint};
use crypto_lib::cryptanalysis::ic_framework::prime_tower::{
    isogeny_tower, supersingular_aux, supersingular_tower, AuxKind, IsogenyTower,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

// ── Field arithmetic ───────────────────────────────────────────────

fn mulm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn addm(a: u64, b: u64, p: u64) -> u64 {
    let s = a + b;
    if s >= p {
        s - p
    } else {
        s
    }
}

fn subm(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + p - b
    }
}

fn powm(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64 % p;
    a %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulm(r, a, p);
        }
        a = mulm(a, a, p);
        e >>= 1;
    }
    r
}

fn invm(a: u64, p: u64) -> u64 {
    assert!(!a.is_multiple_of(p), "inverse of zero");
    powm(a, p - 2, p)
}

fn is_prime(n: u64) -> bool {
    n >= 2
        && (2..)
            .take_while(|d| d * d <= n)
            .all(|d| !n.is_multiple_of(d))
}

// ── Polynomials (copied from pkm_tower_pilot) ──────────────────────

#[derive(Clone)]
struct Pol {
    n: usize,
    terms: BTreeMap<Vec<u32>, u64>,
}

impl Pol {
    fn zero(n: usize) -> Self {
        Pol {
            n,
            terms: BTreeMap::new(),
        }
    }

    fn constant(n: usize, c: u64, p: u64) -> Self {
        let mut f = Pol::zero(n);
        if !c.is_multiple_of(p) {
            f.terms.insert(vec![0; n], c % p);
        }
        f
    }

    fn var(n: usize, i: usize) -> Self {
        let mut e = vec![0; n];
        e[i] = 1;
        let mut f = Pol::zero(n);
        f.terms.insert(e, 1);
        f
    }

    fn add(&self, o: &Pol, p: u64) -> Pol {
        let mut f = self.clone();
        for (e, c) in &o.terms {
            let v = f.terms.entry(e.clone()).or_insert(0);
            *v = addm(*v, *c, p);
        }
        f.terms.retain(|_, c| *c != 0);
        f
    }

    fn scale(&self, k: u64, p: u64) -> Pol {
        let mut f = Pol::zero(self.n);
        for (e, c) in &self.terms {
            let v = mulm(*c, k, p);
            if v != 0 {
                f.terms.insert(e.clone(), v);
            }
        }
        f
    }

    fn sub(&self, o: &Pol, p: u64) -> Pol {
        self.add(&o.scale(p - 1, p), p)
    }

    fn mul(&self, o: &Pol, p: u64) -> Pol {
        let mut f = Pol::zero(self.n);
        for (e1, c1) in &self.terms {
            for (e2, c2) in &o.terms {
                let e: Vec<u32> = e1.iter().zip(e2).map(|(a, b)| a + b).collect();
                let v = f.terms.entry(e).or_insert(0);
                *v = addm(*v, mulm(*c1, *c2, p), p);
            }
        }
        f.terms.retain(|_, c| *c != 0);
        f
    }

    fn eval(&self, x: &[u64], p: u64) -> u64 {
        let mut acc = 0u64;
        for (e, c) in &self.terms {
            let mut t = *c;
            for (i, k) in e.iter().enumerate() {
                if *k > 0 {
                    t = mulm(t, powm(x[i], *k as u64, p), p);
                }
            }
            acc = addm(acc, t, p);
        }
        acc
    }

    fn degree(&self) -> u32 {
        self.terms
            .keys()
            .map(|e| e.iter().sum::<u32>())
            .max()
            .unwrap_or(0)
    }
}

/// Semaev's `S₃` for `y² = x³ + ax + b`:
/// `(X1 − X2)²X3² − 2((X1 + X2)(X1X2 + a) + 2b)X3 + (X1X2 − a)² − 4b(X1 + X2)`.
fn s3(x1: &Pol, x2: &Pol, x3: &Pol, a: u64, b: u64, p: u64) -> Pol {
    let n = x1.n;
    let d = x1.sub(x2, p);
    let t1 = d.mul(&d, p).mul(x3, p).mul(x3, p);
    let s = x1.add(x2, p);
    let pr = x1.mul(x2, p);
    let inner = s
        .mul(&pr.add(&Pol::constant(n, a, p), p), p)
        .add(&Pol::constant(n, mulm(2, b, p), p), p);
    let t2 = inner.mul(x3, p).scale(p - 2, p);
    let q = pr.sub(&Pol::constant(n, a, p), p);
    let t3 = q.mul(&q, p).sub(&s.scale(mulm(4, b, p), p), p);
    t1.add(&t2, p).add(&t3, p)
}

/// `S₃(x1, x2, X)` as the quadratic `(A, B, C)` in `X`, numerically.
fn s3_quad(x1: u64, x2: u64, a: u64, b: u64, p: u64) -> (u64, u64, u64) {
    let d = subm(x1, x2, p);
    let s = addm(x1, x2, p);
    let pr = mulm(x1, x2, p);
    let qa = mulm(d, d, p);
    let inner = addm(mulm(s, addm(pr, a, p), p), mulm(2, b, p), p);
    let qb = subm(0, mulm(2, inner, p), p);
    let q = subm(pr, a, p);
    let qc = subm(mulm(q, q, p), mulm(mulm(4, b, p), s, p), p);
    (qa, qb, qc)
}

// ── Towers ─────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
enum Step {
    /// `y_{j+1} = y_j² − 2` (Dickson `D_2`).
    Dickson,
    /// `y_{j+1} = y_j + τ/(y_j − ξ)` (Vélu).
    Velu { xi: u64, tau: u64 },
}

struct Tower {
    t: usize,
    v: Vec<u64>,
    steps: Vec<Step>,
    c: u64,
    info: Value,
}

impl Tower {
    fn apply(&self, j: usize, y: u64, p: u64) -> u64 {
        match self.steps[j] {
            Step::Dickson => subm(mulm(y, y, p), 2, p),
            Step::Velu { xi, tau } => addm(y, mulm(tau, invm(subm(y, xi, p), p), p), p),
        }
    }

    fn values(&self, x: u64, p: u64) -> Vec<u64> {
        let mut y = vec![0u64; self.t];
        y[0] = x;
        for j in 0..self.t - 1 {
            y[j + 1] = self.apply(j, y[j], p);
        }
        y
    }

    fn equations(&self, n: usize, base: usize, p: u64) -> Vec<Pol> {
        (0..self.t)
            .map(|j| {
                let yj = Pol::var(n, base + j);
                let next = if j + 1 < self.t {
                    Pol::var(n, base + j + 1)
                } else {
                    Pol::constant(n, self.c, p)
                };
                let sq = yj.mul(&yj, p);
                match self.steps[j] {
                    Step::Dickson => next.sub(&sq.sub(&Pol::constant(n, 2, p), p), p),
                    Step::Velu { xi, tau } => {
                        let lhs = next.mul(&yj.sub(&Pol::constant(n, xi, p), p), p);
                        let rhs = sq
                            .sub(&yj.scale(xi, p), p)
                            .add(&Pol::constant(n, tau, p), p);
                        lhs.sub(&rhs, p)
                    }
                }
            })
            .collect()
    }

    fn square_rules(&self, base: usize, p: u64) -> Vec<SquareRule> {
        (0..self.t)
            .map(|j| {
                let top = j + 1 == self.t;
                let next = (!top).then_some((base + j + 1) as u8);
                match (self.steps[j], top) {
                    (Step::Dickson, false) => SquareRule {
                        next,
                        a: 0,
                        b: 1,
                        c: 0,
                        d: 2 % p,
                    },
                    (Step::Dickson, true) => SquareRule {
                        next,
                        a: 0,
                        b: 0,
                        c: 0,
                        d: addm(self.c, 2, p),
                    },
                    (Step::Velu { xi, tau }, false) => SquareRule {
                        next,
                        a: 1,
                        b: subm(0, xi, p),
                        c: xi,
                        d: subm(0, tau, p),
                    },
                    (Step::Velu { xi, tau }, true) => SquareRule {
                        next,
                        a: 0,
                        b: 0,
                        c: addm(self.c, xi, p),
                        d: subm(0, addm(mulm(self.c, xi, p), tau, p), p),
                    },
                }
            })
            .collect()
    }

    /// `|V| = 2^t`, distinct, and every level equation holds on `V`.
    fn check(&self, p: u64) {
        let mut sorted = self.v.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), 1 << self.t, "V has repeats");
        let eqs = self.equations(self.t, 0, p);
        for &x in &self.v {
            let y = self.values(x, p);
            for (j, e) in eqs.iter().enumerate() {
                assert_eq!(e.eval(&y, p), 0, "level {j} fails at x = {x}");
            }
        }
    }
}

type Fp2 = (u64, u64);

fn fp2_mul(x: Fp2, y: Fp2, p: u64) -> Fp2 {
    (
        subm(mulm(x.0, y.0, p), mulm(x.1, y.1, p), p),
        addm(mulm(x.0, y.1, p), mulm(x.1, y.0, p), p),
    )
}

fn fp2_pow(mut x: Fp2, mut e: u64, p: u64) -> Fp2 {
    let mut r = (1u64, 0u64);
    while e > 0 {
        if e & 1 == 1 {
            r = fp2_mul(r, x, p);
        }
        x = fp2_mul(x, x, p);
        e >>= 1;
    }
    r
}

fn torus_element(p: u64, rng: &mut StdRng) -> Fp2 {
    loop {
        let z = (rng.gen_range(0..p), rng.gen_range(0..p));
        if z == (0, 0) {
            continue;
        }
        return fp2_pow(z, p - 1, p);
    }
}

/// A primitive `2^t`-th root of unity in the torus `T ⊂ F_{p²}^*`.
fn torus_zeta(p: u64, t: usize, rng: &mut StdRng) -> Fp2 {
    loop {
        let h = torus_element(p, rng);
        if fp2_pow(h, p.div_ceil(2), p) == (p - 1, 0) {
            let z = fp2_pow(h, (p + 1) >> t, p);
            assert_eq!(fp2_pow(z, 1u64 << t, p), (1, 0));
            assert_eq!(fp2_pow(z, 1u64 << (t - 1), p), (p - 1, 0));
            return z;
        }
    }
}

/// The Dickson torus tower (pilot's `dickson_torus`, its checks inline).
fn dickson_torus(p: u64, t: usize, rng: &mut StdRng) -> Tower {
    assert!(p % 4 == 3 && (p + 1).is_multiple_of(1u64 << t));
    let zeta = torus_zeta(p, t, rng);
    loop {
        let u = torus_element(p, rng);
        if u.1 == 0 || fp2_pow(u, 1u64 << (t + 1), p) == (1, 0) {
            continue;
        }
        let mut v = Vec::with_capacity(1 << t);
        let mut w = u;
        for _ in 0..1u64 << t {
            v.push(addm(w.0, w.0, p));
            w = fp2_mul(w, zeta, p);
        }
        let mut s = v.clone();
        s.sort_unstable();
        s.dedup();
        if s.len() != 1 << t {
            continue;
        }
        let ut = fp2_pow(u, 1u64 << t, p);
        let c = addm(ut.0, ut.0, p);
        let tower = Tower {
            t,
            v,
            steps: vec![Step::Dickson; t],
            c,
            info: json!({"case": "torus", "u": [u.0, u.1], "zeta": [zeta.0, zeta.1], "c": c}),
        };
        // Torus case: y² − 4 a non-residue on V; unramified levels.
        for &x in &tower.v {
            let disc = subm(mulm(x, x, p), 4, p);
            assert_eq!(powm(disc, (p - 1) / 2, p), p - 1, "not the torus case");
            assert!(tower.values(x, p).iter().all(|&y| y != 0), "ramified");
        }
        assert!(tower.c != 2 && tower.c != p - 2);
        return tower;
    }
}

/// The split Dickson tower (`2^t | p − 1`), the pilot's `dickson`.
fn dickson_split(p: u64, t: usize, rng: &mut StdRng) -> Tower {
    assert!((p - 1).is_multiple_of(1u64 << t));
    let zeta = loop {
        let h = rng.gen_range(2..p);
        if powm(h, (p - 1) / 2, p) == p - 1 {
            break powm(h, (p - 1) >> t, p);
        }
    };
    loop {
        let u = rng.gen_range(2..p);
        if powm(u, 1u64 << (t + 1), p) == 1 {
            continue;
        }
        let v: Vec<u64> = (0..1u64 << t)
            .map(|k| {
                let w = mulm(u, powm(zeta, k, p), p);
                addm(w, invm(w, p), p)
            })
            .collect();
        let mut s = v.clone();
        s.sort_unstable();
        s.dedup();
        if s.len() != 1 << t {
            continue;
        }
        let ut = powm(u, 1u64 << t, p);
        return Tower {
            t,
            v,
            steps: vec![Step::Dickson; t],
            c: addm(ut, invm(ut, p), p),
            info: json!({"case": "split", "u": u, "zeta": zeta}),
        };
    }
}

fn from_isogeny(tw: &IsogenyTower, kind: &str) -> Tower {
    let ordinary = tw.aux_order != tw.p + 1;
    let info = json!({
        "aux": kind, "aux_a": tw.aux.a, "aux_b": tw.aux.b, "aux_order": tw.aux_order,
        "aux_supersingular": !ordinary, "aux_curves_tried": tw.cost.curves_tried,
        "q": [tw.q.x, tw.q.y], "r": [tw.r.x, tw.r.y], "c": tw.c,
        "xi": tw.steps.iter().map(|s| s.xi).collect::<Vec<_>>(),
        "tau": tw.steps.iter().map(|s| s.tau).collect::<Vec<_>>(),
        "v": &tw.v,
    });
    Tower {
        t: tw.t as usize,
        v: tw.v.clone(),
        steps: tw
            .steps
            .iter()
            .map(|s| Step::Velu {
                xi: s.xi,
                tau: s.tau,
            })
            .collect(),
        c: tw.c,
        info,
    }
}

fn build_arm(arm: &str, p: u64, t: usize, seed: u64, rng: &mut StdRng) -> Tower {
    let tower = match arm {
        "torus" => dickson_torus(p, t, rng),
        "split" => dickson_split(p, t, rng),
        "ordinary" => {
            let tw = isogeny_tower(p, t as u32, seed, &mut GroupOps::default())
                .expect("ordinary isogeny tower");
            assert_ne!(tw.aux_order, p + 1, "the searched E' is supersingular");
            from_isogeny(&tw, "ordinary")
        }
        ss => {
            let kind = AuxKind::parse(ss, p).expect("arm");
            let tw = supersingular_tower(p, t as u32, seed, kind, &mut GroupOps::default())
                .expect("supersingular tower");
            assert_eq!(tw.aux_order, p + 1);
            from_isogeny(&tw, kind.name())
        }
    };
    tower.check(p);
    tower
}

// ── Systems ────────────────────────────────────────────────────────

fn nonsingular(a: u64, b: u64, p: u64) -> bool {
    let d = addm(
        mulm(4, mulm(mulm(a, a, p), a, p), p),
        mulm(27, mulm(b, b, p), p),
        p,
    );
    d != 0
}

fn random_point(curve: &PrimeCurve, rng: &mut StdRng) -> PrimePoint {
    loop {
        let x = rng.gen_range(0..curve.p);
        let r = curve.rhs(x);
        if r == 0 {
            continue;
        }
        if let Some(y) = curve.sqrt(r) {
            let y = if rng.gen::<bool>() { y } else { curve.p - y };
            return PrimePoint::affine(x, y);
        }
    }
}

/// A random target curve with at least `need` abscissae of `V` on it
/// (the pilot's `target_curve`).
fn target_curve(p: u64, tower: &Tower, need: usize, rng: &mut StdRng) -> PrimeCurve {
    loop {
        let a = rng.gen_range(0..p);
        let b = rng.gen_range(1..p);
        if !nonsingular(a, b, p) {
            continue;
        }
        let curve = PrimeCurve { p, a, b };
        let on = tower
            .v
            .iter()
            .filter(|&&x| curve.legendre(curve.rhs(x)) == 1)
            .count();
        if on >= need {
            return curve;
        }
    }
}

fn summation(tower: &Tower, curve: &PrimeCurve, m: usize, n: usize, x_r: u64) -> Vec<Pol> {
    let p = curve.p;
    let (a, b) = (curve.a, curve.b);
    let x = |i: usize| Pol::var(n, i * tower.t);
    let xr = Pol::constant(n, x_r, p);
    if m == 2 {
        return vec![s3(&x(0), &x(1), &xr, a, b, p)];
    }
    let u = |k: usize| Pol::var(n, m * tower.t + k);
    let mut eqs = vec![s3(&x(0), &x(1), &u(0), a, b, p)];
    for k in 1..m - 2 {
        eqs.push(s3(&u(k - 1), &x(k + 1), &u(k), a, b, p));
    }
    eqs.push(s3(&u(m - 3), &x(m - 1), &xr, a, b, p));
    eqs
}

/// The summation equations (`null`: each `S₃` link's support with random
/// coefficients, as the pilot's null and `--null-chain`), and the tower
/// ring. Random targets only, so the null has no planted zero.
fn system(
    null: bool,
    tower: &Tower,
    curve: &PrimeCurve,
    m: usize,
    x_r: u64,
    rng: &mut StdRng,
) -> (Vec<Pol>, usize, TowerRing) {
    let p = curve.p;
    let n = m * tower.t + (m - 2);
    let mut eqs = summation(tower, curve, m, n, x_r);
    if null {
        for eq in eqs.iter_mut() {
            let mut f = Pol::zero(n);
            for e in eq.terms.keys() {
                f.terms.insert(e.clone(), rng.gen_range(1..p));
            }
            *eq = f;
        }
    }
    let mut rules = Vec::with_capacity(m * tower.t);
    for i in 0..m {
        rules.extend(tower.square_rules(i * tower.t, p));
    }
    let ring = TowerRing {
        p,
        rules,
        n_free: m - 2,
    };
    // The square rules are the tower equations.
    for i in 0..m {
        for e in tower.equations(n, i * tower.t, p) {
            assert!(
                ring.from_raw(&raw_terms(&e)).is_empty(),
                "a tower equation is not zero in the ring"
            );
        }
    }
    (eqs, n, ring)
}

fn raw_terms(e: &Pol) -> Vec<(Vec<u32>, u64)> {
    e.terms.iter().map(|(k, v)| (k.clone(), *v)).collect()
}

/// Exhaustive recount over `V^m`: does the target decompose (`m = 2`:
/// `S₃(x1, x2, x_R) = 0`; `m = 3`: a common root `u ∈ F̄_p` of
/// `S₃(x1, x2, u)` and `S₃(u, x3, x_R)`, by their resultant, which is also
/// zero when both leading coefficients vanish)? `None` for `m ≥ 4`.
fn decomposes(tower: &Tower, curve: &PrimeCurve, m: usize, x_r: u64) -> Option<bool> {
    let p = curve.p;
    let (a, b) = (curve.a, curve.b);
    let v = &tower.v;
    match m {
        2 => Some(v.iter().any(|&x1| {
            v.iter().any(|&x2| {
                let (qa, qb, qc) = s3_quad(x1, x2, a, b, p);
                addm(
                    addm(mulm(qa, mulm(x_r, x_r, p), p), mulm(qb, x_r, p), p),
                    qc,
                    p,
                ) == 0
            })
        })),
        3 => Some(v.iter().any(|&x1| {
            v.iter().any(|&x2| {
                let (a1, b1, c1) = s3_quad(x1, x2, a, b, p);
                v.iter().any(|&x3| {
                    let (a2, b2, c2) = s3_quad(x3, x_r, a, b, p);
                    let r1 = subm(mulm(a1, c2, p), mulm(a2, c1, p), p);
                    let r2 = subm(mulm(a1, b2, p), mulm(a2, b1, p), p);
                    let r3 = subm(mulm(b1, c2, p), mulm(b2, c1, p), p);
                    subm(mulm(r1, r1, p), mulm(r2, r3, p), p) == 0
                })
            })
        })),
        _ => None,
    }
}

// ── Q1: Möbius classes of tower sets ───────────────────────────────

/// The least image of `set` under the Möbius maps sending an ordered
/// triple of its points to `(0, 1, ∞)`: a complete `PGL₂(F_p)` invariant
/// of the set (all points finite, `|set| ≥ 3`).
fn canonical(set: &[u64], p: u64) -> Vec<u64> {
    let n = set.len();
    let mut best: Option<Vec<u64>> = None;
    let mut img = Vec::with_capacity(n - 3);
    for i in 0..n {
        for j in 0..n {
            if j == i {
                continue;
            }
            for k in 0..n {
                if k == i || k == j {
                    continue;
                }
                let (a, b, c) = (set[i], set[j], set[k]);
                // μ(x) = (x − a)(b − c) / ((x − c)(b − a))
                let bc = subm(b, c, p);
                let ba = subm(b, a, p);
                img.clear();
                for (l, &x) in set.iter().enumerate() {
                    if l == i || l == j || l == k {
                        continue;
                    }
                    let num = mulm(subm(x, a, p), bc, p);
                    let den = mulm(subm(x, c, p), ba, p);
                    img.push(mulm(num, invm(den, p), p));
                }
                img.sort_unstable();
                if best.as_ref().is_none_or(|bst| img < *bst) {
                    best = Some(img.clone());
                }
            }
        }
    }
    best.unwrap_or_default()
}

/// Brute force over `PGL₂(F_p)`: is there `μ` with `μ(A) = B`? (Small `p`
/// only; `∞` is encoded as `p`.)
fn brute_equivalent(a: &[u64], b: &[u64], p: u64) -> bool {
    let target: HashSet<u64> = b.iter().copied().collect();
    let act = |m: (u64, u64, u64, u64), x: u64| -> u64 {
        let (al, be, ga, de) = m;
        let num = addm(mulm(al, x, p), be, p);
        let den = addm(mulm(ga, x, p), de, p);
        if den == 0 {
            p
        } else {
            mulm(num, invm(den, p), p)
        }
    };
    // Normalise: γ = 0, δ = 1 (affine), or γ = 1.
    for ga in 0..2u64 {
        for al in 0..p {
            for be in 0..p {
                let des: Vec<u64> = if ga == 0 { vec![1] } else { (0..p).collect() };
                for de in des {
                    let det = subm(mulm(al, de, p), mulm(be, ga, p), p);
                    if det == 0 {
                        continue;
                    }
                    if a.iter()
                        .all(|&x| target.contains(&act((al, be, ga, de), x)))
                    {
                        return true;
                    }
                }
            }
        }
    }
    false
}

/// All cosets `x(R' + ⟨Q⟩)` (as sorted sets) on `e` with `#e = n`, over
/// every cyclic subgroup of order `2^t` and every `R'` with `2R' ∉ ⟨Q⟩`.
fn isogeny_sets(e: &PrimeCurve, n: u64, t: usize) -> Vec<Vec<u64>> {
    let p = e.p;
    let mut ops = GroupOps::default();
    let mut points = vec![PrimePoint::INFINITY];
    for x in 0..p {
        let r = e.rhs(x);
        if r == 0 {
            points.push(PrimePoint::affine(x, 0));
        } else if let Some(y) = e.sqrt(r) {
            if e.legendre(r) == 1 {
                points.push(PrimePoint::affine(x, y));
                points.push(PrimePoint::affine(x, p - y));
            }
        }
    }
    assert_eq!(points.len() as u64, n, "point count");
    let two_t = 1u64 << t;
    // Cyclic subgroups of order 2^t, by their point sets.
    let mut subgroups: Vec<Vec<PrimePoint>> = Vec::new();
    let mut seen: HashSet<Vec<u64>> = HashSet::new();
    for &q in &points {
        if q.infinity {
            continue;
        }
        let half = e.mul(&mut ops, q, two_t / 2);
        if half.infinity || !e.mul(&mut ops, q, two_t).infinity {
            continue;
        }
        let mut sub = Vec::with_capacity(two_t as usize);
        let mut acc = PrimePoint::INFINITY;
        for _ in 0..two_t {
            sub.push(acc);
            acc = e.add(&mut ops, acc, q);
        }
        let mut key: Vec<u64> = sub.iter().map(|pt| e.key(pt)).collect();
        key.sort_unstable();
        if seen.insert(key) {
            subgroups.push(sub);
        }
    }
    let mut sets: HashSet<Vec<u64>> = HashSet::new();
    for sub in &subgroups {
        let keys: HashSet<u64> = sub.iter().map(|pt| e.key(pt)).collect();
        for &r in &points {
            if r.infinity {
                continue;
            }
            let r2 = e.double(&mut ops, r);
            if keys.contains(&e.key(&r2)) {
                continue;
            }
            let mut v: Vec<u64> = sub
                .iter()
                .map(|&m| e.add(&mut ops, r, m))
                .map(|pt| {
                    assert!(!pt.infinity);
                    pt.x
                })
                .collect();
            v.sort_unstable();
            let len = v.len();
            v.dedup();
            assert_eq!(v.len(), len, "x not injective on the coset");
            sets.insert(v);
        }
    }
    let mut out: Vec<Vec<u64>> = sets.into_iter().collect();
    out.sort();
    out
}

/// All Dickson torus sets `{2·re(w) : w ∈ u·μ_{2^t}}`, `u ∈ T \ F_p`,
/// `u² ∉ μ_{2^t}`, `|V| = 2^t`.
fn torus_sets(p: u64, t: usize) -> Vec<Vec<u64>> {
    let mut rng = StdRng::seed_from_u64(p ^ 0x7031);
    let zeta = torus_zeta(p, t, &mut rng);
    let mut sets: HashSet<Vec<u64>> = HashSet::new();
    let e = PrimeCurve { p, a: 0, b: 1 };
    // T = {a + bi : a² + b² = 1}.
    for re in 0..p {
        let r = subm(1, mulm(re, re, p), p);
        if r == 0 || e.legendre(r) < 0 {
            continue;
        }
        let im0 = e.sqrt(r).expect("square root");
        for im in [im0, p - im0] {
            let u = (re, im);
            if fp2_pow(u, 1u64 << (t + 1), p) == (1, 0) {
                continue;
            }
            let mut v = Vec::with_capacity(1 << t);
            let mut w = u;
            for _ in 0..1u64 << t {
                v.push(addm(w.0, w.0, p));
                w = fp2_mul(w, zeta, p);
            }
            v.sort_unstable();
            let len = v.len();
            v.dedup();
            if v.len() == len {
                sets.insert(v);
            }
        }
    }
    let mut out: Vec<Vec<u64>> = sets.into_iter().collect();
    out.sort();
    out
}

/// `j = 1728·4a³/(4a³ + 27b²)`.
fn j_invariant(a: u64, b: u64, p: u64) -> u64 {
    let a3 = mulm(4, mulm(mulm(a, a, p), a, p), p);
    let den = addm(a3, mulm(27, mulm(b, b, p), p), p);
    mulm(mulm(1728 % p, a3, p), invm(den, p), p)
}

fn run_structure(args: &HashMap<String, String>) {
    let p_min: u64 = args.get("p-min").map_or(100, |v| v.parse().unwrap());
    let p_max: u64 = args.get("p-max").map_or(4000, |v| v.parse().unwrap());
    let ts: Vec<usize> = args
        .get("t")
        .map_or("2,3,4", |v| v.as_str())
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let n_ordinary: usize = args.get("ordinary").map_or(3, |v| v.parse().unwrap());
    let brute_max: u64 = args.get("brute").map_or(0, |v| v.parse().unwrap());
    let mut out = File::create(args.get("out").expect("--out")).expect("out");
    let kinds = [
        AuxKind::SsJ1728,
        AuxKind::SsJ1728Twist,
        AuxKind::SsJ0,
        AuxKind::SsJ0Twist,
    ];
    for p in (p_min..=p_max).filter(|&p| is_prime(p) && p % 4 == 3) {
        let v2 = (p + 1).trailing_zeros() as usize;
        for &t in &ts {
            // t + 1 ≤ v₂(p + 1) so that every family below has its towers
            // (the j = 1728 twist's 2-Sylow is Z/2 × Z/2^{v−1}).
            if t < 2 || t + 1 > v2 {
                continue;
            }
            let started = Instant::now();
            let mut families: Vec<(String, Vec<Vec<u64>>, Value)> = Vec::new();
            families.push(("torus".into(), torus_sets(p, t), json!(null)));
            for kind in kinds {
                if let Ok(e) = supersingular_aux(p, kind) {
                    let sets = isogeny_sets(&e, p + 1, t);
                    families.push((
                        kind.name().into(),
                        sets,
                        json!({"a": e.a, "b": e.b, "j": j_invariant(e.a, e.b, p)}),
                    ));
                }
            }
            // Positive control: y² = x³ + 16x ≅ y² = x³ + x by (x, y) ↦
            // (4x, 8y), so its sets are the ss-j1728 sets scaled by 4 and
            // every one must be found equivalent.
            let iso = PrimeCurve { p, a: 16, b: 0 };
            families.push((
                "ss-j1728-iso16".into(),
                isogeny_sets(&iso, p + 1, t),
                json!({"a": 16, "b": 0, "j": j_invariant(16, 0, p), "positive_control": true}),
            ));
            // Ordinary controls: the first curves y² = x³ + ax + b (a, b
            // from a fixed seed) whose order has a point of order 2^t.
            let mut rng = StdRng::seed_from_u64(p.wrapping_mul(31) ^ t as u64);
            let mut found = 0;
            let mut tried = 0;
            while found < n_ordinary && tried < 20_000 {
                tried += 1;
                let a = rng.gen_range(0..p);
                let b = rng.gen_range(1..p);
                if !nonsingular(a, b, p) {
                    continue;
                }
                let e = PrimeCurve { p, a, b };
                let n = e.point_count();
                if n == p + 1 || !n.is_multiple_of(1 << t) {
                    continue;
                }
                let sets = isogeny_sets(&e, n, t);
                if sets.is_empty() {
                    continue;
                }
                found += 1;
                families.push((
                    format!("ordinary-{found}"),
                    sets,
                    json!({"a": a, "b": b, "order": n, "j": j_invariant(a, b, p)}),
                ));
            }
            let canon: Vec<HashSet<Vec<u64>>> = families
                .iter()
                .map(|(_, sets, _)| sets.iter().map(|s| canonical(s, p)).collect())
                .collect();
            for (fa, (name_a, sets_a, info_a)) in families.iter().enumerate() {
                for (fb, (name_b, sets_b, _)) in families.iter().enumerate() {
                    let hits = sets_a
                        .iter()
                        .filter(|s| canon[fb].contains(&canonical(s, p)))
                        .count();
                    let classes_a = canon[fa].len();
                    let shared = canon[fa].intersection(&canon[fb]).count();
                    let mut row = json!({
                        "p": p, "t": t, "v2_p_plus_1": v2,
                        "family_a": name_a, "family_b": name_b,
                        "sets_a": sets_a.len(), "sets_b": sets_b.len(),
                        "classes_a": classes_a, "classes_b": canon[fb].len(),
                        "classes_shared": shared,
                        "sets_a_equivalent_to_some_b": hits,
                        "fraction": hits as f64 / sets_a.len().max(1) as f64,
                        "curve_a": info_a,
                    });
                    // Brute-force confirmation of the invariant at small p.
                    if p <= brute_max && fa != fb {
                        let mut agree = 0;
                        let mut checked = 0;
                        for s in sets_a.iter().take(6) {
                            let by_canon = canon[fb].contains(&canonical(s, p));
                            let by_brute = sets_b.iter().any(|sb| brute_equivalent(s, sb, p));
                            checked += 1;
                            if by_canon == by_brute {
                                agree += 1;
                            }
                        }
                        row["brute_checked"] = json!(checked);
                        row["brute_agree"] = json!(agree);
                    }
                    writeln!(out, "{row}").unwrap();
                }
            }
            eprintln!(
                "p = {p} t = {t}: {} families, {:.1} s",
                families.len(),
                started.elapsed().as_secs_f64()
            );
        }
    }
}

// ── Q2: F4 ─────────────────────────────────────────────────────────

fn run_f4(args: &HashMap<String, String>) {
    let p: u64 = args.get("p").expect("--p").parse().unwrap();
    let arm = args.get("arm").expect("--arm").clone();
    let null = args.get("control").is_some_and(|c| c == "null");
    let ms: Vec<usize> = args
        .get("m")
        .map_or("2", |v| v.as_str())
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let t_min: usize = args.get("t-min").map_or(1, |v| v.parse().unwrap());
    let t_max: usize = args.get("t-max").map_or(5, |v| v.parse().unwrap());
    let t_max_m3: usize = args.get("t-max-m3").map_or(3, |v| v.parse().unwrap());
    let t_max_m4: usize = args.get("t-max-m4").map_or(2, |v| v.parse().unwrap());
    let targets: usize = args.get("targets").map_or(2, |v| v.parse().unwrap());
    let budget: u64 = args.get("budget").map_or(3600, |v| v.parse().unwrap());
    let max_nnz: u64 = args
        .get("max-nnz")
        .map_or(2_000_000_000, |v| v.parse().unwrap());
    let seed: u64 = args.get("seed").map_or(0x5353_4157, |v| v.parse().unwrap());
    let tower_seed: u64 = args.get("tower-seed").map_or(0, |v| v.parse().unwrap());
    let mut out = File::create(args.get("out").expect("--out")).expect("out");
    let arm_tag = arm
        .bytes()
        .fold(0u64, |h, b| h.wrapping_mul(131) ^ b as u64);
    for &m in &ms {
        let tmax = match m {
            2 => t_max,
            3 => t_max_m3,
            _ => t_max_m4,
        };
        for t in t_min.max(1)..=tmax {
            // The tower seed is shared by the tower and null controls of an
            // arm, so a null runs on the tower its arm measured.
            let cell_seed =
                seed ^ arm_tag.rotate_left(40) ^ ((m as u64) << 32) ^ ((t as u64) << 16);
            let mut rng = StdRng::seed_from_u64(cell_seed);
            let tower = build_arm(&arm, p, t, tower_seed ^ (t as u64), &mut rng);
            let curve = target_curve(p, &tower, m.max(2), &mut rng);
            let xrs: Vec<u64> = (0..targets)
                .map(|_| random_point(&curve, &mut rng).x)
                .collect();
            let mut null_rng = StdRng::seed_from_u64(cell_seed ^ 0x6E75_6C6C);
            for (k, &x_r) in xrs.iter().enumerate() {
                let (eqs, n, ring) = system(null, &tower, &curve, m, x_r, &mut null_rng);
                let d_in: Vec<u32> = eqs.iter().map(Pol::degree).collect();
                let d_max = d_in.iter().copied().max().unwrap_or(1);
                let max_degree = n as u32 + d_max + 6;
                let input: Vec<_> = eqs
                    .iter()
                    .map(|e| ring.from_raw(&raw_terms(e)))
                    .filter(|f| !f.is_empty())
                    .collect();
                let mut opts = TowerF4Options::new(max_degree)
                    .with_budget(Duration::from_secs(budget))
                    .with_max_nnz(max_nnz);
                opts.full_rank_exit = true;
                let started = Instant::now();
                let r = f4_fp_tower::f4_tower(&input, &ring, &opts);
                let wall = started.elapsed().as_secs_f64();
                let finished = !(r.timed_out || r.oversize);
                // Exhaustive recount (tower control only: the null is not a
                // summation system).
                let verify = if !null && finished {
                    decomposes(&tower, &curve, m, x_r).map(|dec| dec == !r.inconsistent)
                } else {
                    None
                };
                let row = json!({
                    "arm": arm, "kind": arm, "control": if null { "null" } else { "tower" },
                    "p": p, "m": m, "t": t, "N": m * t, "n_vars": n,
                    "input_degrees": d_in, "target_index": k, "x_r": x_r,
                    "curve": {"a": curve.a, "b": curve.b},
                    "tower": tower.info,
                    "max_degree_bound": max_degree,
                    "solving_degree_max": r.solving_degree_max,
                    "last_productive_degree": r.last_productive_degree,
                    "degree_reached": r.degree_reached,
                    "max_cols_to_solution": r.max_cols_to_solution,
                    "steps_to_solution": r.steps_to_solution,
                    "inconsistent": r.inconsistent,
                    "basis_len": r.basis.len(),
                    "steps": r.steps, "max_rows": r.max_rows, "max_cols": r.max_cols,
                    "muladds": r.muladds, "max_nnz": r.max_nnz,
                    "timed_out": r.timed_out, "oversize": r.oversize,
                    "ms": r.ms, "wall_s": wall,
                    "engine": "f4_fp_tower", "full_rank_exit": true,
                    "seed": seed, "tower_seed": tower_seed,
                    "verify": verify,
                });
                eprintln!(
                    "{arm:15} {:5} m={m} t={t:2} N={:2} tgt={k} D={} width={} steps={} inc={} {:.1} ms{}{}",
                    if null { "null" } else { "tower" },
                    m * t,
                    r.solving_degree_max,
                    r.max_cols_to_solution,
                    r.steps_to_solution,
                    r.inconsistent,
                    r.ms,
                    if finished { "" } else { "  STOPPED" },
                    if verify == Some(false) { "  VERIFY-MISMATCH" } else { "" }
                );
                writeln!(out, "{row}").unwrap();
                out.flush().unwrap();
                assert_ne!(
                    verify,
                    Some(false),
                    "F4 disagrees with the exhaustive recount"
                );
            }
        }
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let cmd = argv.first().cloned().unwrap_or_default();
    let mut args = HashMap::new();
    let mut it = argv.iter().skip(1);
    while let Some(k) = it.next() {
        let key = k.trim_start_matches("--").to_string();
        let val = it.next().cloned().unwrap_or_default();
        args.insert(key, val);
    }
    match cmd.as_str() {
        "structure" => run_structure(&args),
        "f4" => run_f4(&args),
        _ => {
            eprintln!("usage: pkm_ss_tower structure|f4 [--flags]");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_forms_agree_with_brute_force() {
        let p = 23u64;
        let a = vec![1u64, 2, 5, 7];
        // An explicit Möbius image of a: x ↦ (3x + 1)/(x + 4).
        let img: Vec<u64> = a
            .iter()
            .map(|&x| mulm(addm(mulm(3, x, p), 1, p), invm(addm(x, 4, p), p), p))
            .collect();
        assert_eq!(canonical(&a, p), canonical(&img, p));
        assert!(brute_equivalent(&a, &img, p));
        let other = vec![1u64, 2, 5, 9];
        assert_eq!(
            canonical(&a, p) == canonical(&other, p),
            brute_equivalent(&a, &other, p)
        );
    }

    #[test]
    fn supersingular_and_torus_sets_are_complete() {
        // p = 127: v₂(p + 1) = 7. Every ss-j1728 coset set is a tower set of
        // the library builder's form, and the torus sets have size 2^t.
        let p = 127u64;
        let t = 3;
        let e = supersingular_aux(p, AuxKind::SsJ1728).unwrap();
        let sets = isogeny_sets(&e, p + 1, t);
        // (p + 1)/2^t cosets, each with its negative giving the same set
        // only when 2R' ∈ ⟨Q⟩ (excluded): (p + 1 − 2^{t+1})/2^{t+1}.
        assert_eq!(sets.len() as u64, (p + 1 - (1 << (t + 1))) / (1 << (t + 1)));
        let mut rng = StdRng::seed_from_u64(1);
        let tw = build_arm("ss-j1728", p, t, 0, &mut rng);
        let mut v = tw.v.clone();
        v.sort_unstable();
        assert!(sets.contains(&v));
        for s in torus_sets(p, t) {
            assert_eq!(s.len(), 1 << t);
        }
    }

    #[test]
    fn f4_rows_agree_with_the_exhaustive_recount() {
        let p = 14_680_063u64; // 7·2^21 − 1
        for arm in ["ss-j1728", "ss-j1728-twist", "torus", "ordinary"] {
            let mut rng = StdRng::seed_from_u64(5);
            let tower = build_arm(arm, p, 3, 0, &mut rng);
            let curve = target_curve(p, &tower, 2, &mut rng);
            // A planted target: x(P1 + P2) for two base points decomposes.
            let mut ops = GroupOps::default();
            let on: Vec<u64> = tower
                .v
                .iter()
                .copied()
                .filter(|&x| curve.legendre(curve.rhs(x)) == 1)
                .collect();
            let pt = |x: u64| PrimePoint::affine(x, curve.sqrt(curve.rhs(x)).unwrap());
            let s = curve.add(&mut ops, pt(on[0]), pt(on[1]));
            for (x_r, expect) in [(s.x, Some(true)), (random_point(&curve, &mut rng).x, None)] {
                let (eqs, _, ring) = system(false, &tower, &curve, 2, x_r, &mut rng);
                let input: Vec<_> = eqs.iter().map(|e| ring.from_raw(&raw_terms(e))).collect();
                let opts = TowerF4Options::new(20);
                let r = f4_fp_tower::f4_tower(&input, &ring, &opts);
                let dec = decomposes(&tower, &curve, 2, x_r).unwrap();
                assert_eq!(dec, !r.inconsistent, "{arm}");
                if let Some(e) = expect {
                    assert_eq!(dec, e);
                }
            }
        }
    }
}
