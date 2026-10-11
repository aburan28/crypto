//! **The isogeny walk to a weak curve, priced** (the one unpriced item of
//! `RESEARCH_COVER_DECOMPOSITION_LEDGER.md` §9).
//!
//! The cover route of [`super::jv_cover`] applies to curves over `F_{q³}`,
//! `q = p²`, with a model `y² = h(x)(x − α)(x − σα)`, `h ∈ F_q[x]` linear.
//! [JV12] §4.1 count `Θ(q²)` of them among `Θ(q³)` curves and estimate the
//! walk from a curve of order divisible by `4` to a weak isogenous one at
//! `≈ q` low-degree isogeny steps (cited, conjectural).  This module makes
//! both measurable at toy sizes:
//!
//! * **The weak-class test.**  A curve with full rational 2-torsion,
//!   `y² = (x − e₁)(x − e₂)(x − e₃)` over `F_{q³}`, has a model of the weak
//!   form with `e₁ ↦ ρ ∈ F_q` and `e₂, e₃ ↦ α, σα` exactly when the
//!   `F_{q³}`-affine change `φ(x) = (x − r)/v` with `φ(e₁) = ρ`,
//!   `φ(e₃) = σ(φ(e₂))` exists.  Eliminating `r, v` gives
//!   `σ(a) − ρ = c·(a − ρ)` with `a = φ(e₂)` and `c = (e₃ − e₁)/(e₂ − e₁)`,
//!   whose only solution is the degenerate `a = ρ` unless the `F_q`-linear
//!   map `a ↦ σ(a) − c·a` is singular, i.e. **`N_{F_{q³}/F_q}(c) = 1`**.  So
//!   the class is the curves with full 2-torsion one of whose three
//!   cross-ratios `(e₃ − e₁)/(e₂ − e₁)` has norm one: an isomorphism
//!   invariant, one `F_q`-condition, `≈ 3/q` of the curves with full
//!   2-torsion ([JV12]'s `Θ(q²)` of `Θ(q³)`).
//! * **The walk.**  From a random curve with full 2-torsion (so of order
//!   divisible by `4`), steps by the rational `2`-isogenies (Vélu on a
//!   2-torsion point; the target keeps full 2-torsion when the product of
//!   the other two roots is a square) and the rational `3`-isogenies (Vélu on
//!   a root of the 3-division polynomial), chosen uniformly, until a weak
//!   curve is reached; counted in steps, in distinct `j`-invariants, and in
//!   `F_p` multiplications.

use std::collections::{HashMap, HashSet, VecDeque};
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use super::gaudry_cubic::{mm, sm};
use super::jv_cover::{EllE, Fld, Fq3, PtE6, E2, E6};
use super::residual_walk::pow_mod;

// ── small polynomial toolkit over F_{q³} (coefficients low to high) ──────

type P6 = Vec<E6>;

fn ptrim(p: &mut P6) {
    while p.last() == Some(&E6::ZERO) {
        p.pop();
    }
}

fn pmul(f: &Fq3, a: &[E6], b: &[E6]) -> P6 {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![E6::ZERO; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = f.add(&out[i + j], &f.mul(x, y));
        }
    }
    ptrim(&mut out);
    out
}

fn psub(f: &Fq3, a: &[E6], b: &[E6]) -> P6 {
    let n = a.len().max(b.len());
    let mut out: P6 = (0..n)
        .map(|i| {
            f.sub(
                &a.get(i).copied().unwrap_or(E6::ZERO),
                &b.get(i).copied().unwrap_or(E6::ZERO),
            )
        })
        .collect();
    ptrim(&mut out);
    out
}

/// Remainder of `a` modulo `m` (any leading coefficient).
fn prem(f: &Fq3, a: &[E6], m: &[E6]) -> P6 {
    let mut r = a.to_vec();
    ptrim(&mut r);
    let dm = m.len() - 1;
    let li = f.inv(m.last().unwrap());
    while r.len() > dm {
        let coef = f.mul(r.last().unwrap(), &li);
        let shift = r.len() - m.len();
        for (i, c) in m.iter().enumerate() {
            r[shift + i] = f.sub(&r[shift + i], &f.mul(&coef, c));
        }
        ptrim(&mut r);
    }
    r
}

fn pmonic(f: &Fq3, a: &[E6]) -> P6 {
    match a.last() {
        None => Vec::new(),
        Some(l) => {
            let li = f.inv(l);
            a.iter().map(|c| f.mul(c, &li)).collect()
        }
    }
}

fn pgcd(f: &Fq3, a: &[E6], b: &[E6]) -> P6 {
    let (mut x, mut y) = (a.to_vec(), b.to_vec());
    ptrim(&mut x);
    ptrim(&mut y);
    while !y.is_empty() {
        let r = prem(f, &x, &y);
        x = y;
        y = r;
    }
    pmonic(f, &x)
}

fn ppowmod(f: &Fq3, a: &[E6], mut e: u128, m: &[E6]) -> P6 {
    let mut base = prem(f, a, m);
    let mut acc = vec![E6::ONE];
    while e > 0 {
        if e & 1 == 1 {
            acc = prem(f, &pmul(f, &acc, &base), m);
        }
        e >>= 1;
        if e > 0 {
            base = prem(f, &pmul(f, &base, &base), m);
        }
    }
    acc
}

/// The distinct roots in `F_{q³}` of `a` (Cantor–Zassenhaus).
pub fn roots_q3(f: &Fq3, a: &[E6], rng: &mut StdRng) -> Vec<E6> {
    let a = pmonic(f, a);
    if a.len() <= 1 {
        return Vec::new();
    }
    let q3 = (f.f.p as u128).pow(6);
    let x = vec![E6::ZERO, E6::ONE];
    let xq = ppowmod(f, &x, q3, &a);
    let g = pgcd(f, &a, &psub(f, &xq, &x));
    let mut out = Vec::new();
    fn split(f: &Fq3, g: &[E6], q3: u128, rng: &mut StdRng, out: &mut Vec<E6>) {
        if g.len() <= 1 {
            return;
        }
        if g.len() == 2 {
            out.push(f.neg(&g[0]));
            return;
        }
        loop {
            let r: P6 = vec![f.random(rng), E6::ONE];
            let w = psub(f, &ppowmod(f, &r, (q3 - 1) / 2, g), &[E6::ONE]);
            let h = pgcd(f, g, &w);
            if h.len() > 1 && h.len() < g.len() {
                let other = {
                    // g / h
                    let mut q: P6 = Vec::new();
                    let mut r = g.to_vec();
                    let li = f.inv(h.last().unwrap());
                    let dh = h.len() - 1;
                    q.resize(g.len() - dh, E6::ZERO);
                    while r.len() > dh {
                        let coef = f.mul(r.last().unwrap(), &li);
                        let shift = r.len() - h.len();
                        q[shift] = coef;
                        for (i, c) in h.iter().enumerate() {
                            r[shift + i] = f.sub(&r[shift + i], &f.mul(&coef, c));
                        }
                        ptrim(&mut r);
                    }
                    q
                };
                split(f, &h, q3, rng, out);
                split(f, &other, q3, rng, out);
                return;
            }
        }
    }
    split(f, &g, q3, rng, &mut out);
    out
}

// ── curves with full 2-torsion over F_{q³} ──────────────────────────────

/// `y² = (x − e₀)(x − e₁)(x − e₂)` with distinct roots in `F_{q³}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Curve2 {
    pub e: [E6; 3],
}

/// `N_{F_{q³}/F_q}(z) = z·σ(z)·σ²(z)`.
pub fn norm_q3(f: &Fq3, z: &E6) -> E6 {
    f.mul(z, &f.mul(&f.sigma(z), &f.sigma(&f.sigma(z))))
}

impl Curve2 {
    /// The weak-class test: some ordering `(e₁; e₂, e₃)` has
    /// `N((e₃ − e₁)/(e₂ − e₁)) = 1`.  Returns the index of the root that
    /// plays `ρ`, if any.
    pub fn weak_root(&self, f: &Fq3) -> Option<usize> {
        (0..3).find(|&i| {
            let e1 = self.e[i];
            let e2 = self.e[(i + 1) % 3];
            let e3 = self.e[(i + 2) % 3];
            let c = f.mul(&f.sub(&e3, &e1), &f.inv(&f.sub(&e2, &e1)));
            norm_q3(f, &c) == E6::ONE
        })
    }
    /// The `j`-invariant (the class of the curve over `F_{q³}`).
    pub fn j(&self, f: &Fq3) -> E6 {
        // Legendre form y² = x(x − 1)(x − λ) with λ = (e₂ − e₀)/(e₁ − e₀):
        // j = 256 (λ² − λ + 1)³ / (λ² (λ − 1)²)
        let lam = f.mul(
            &f.sub(&self.e[2], &self.e[0]),
            &f.inv(&f.sub(&self.e[1], &self.e[0])),
        );
        let l2 = f.sq(&lam);
        let num = f.add(&f.sub(&l2, &lam), &E6::ONE);
        let num3 = f.mul(&num, &f.sq(&num));
        let den = f.mul(&l2, &f.sq(&f.sub(&lam, &E6::ONE)));
        let c256 = f.from_fq(E2([256 % f.f.p, 0]));
        f.mul(&f.mul(&c256, &num3), &f.inv(&den))
    }
    /// Short Weierstrass `y² = x³ + a x + b` of the translate with
    /// `e₀ + e₁ + e₂ = 0`, and the translation `t` (`x_short = x − t`).
    pub fn short(&self, f: &Fq3) -> (E6, E6, E6) {
        let three_inv = f.inv(&f.from_fq(E2([3 % f.f.p, 0])));
        let t = f.mul(
            &f.add(&f.add(&self.e[0], &self.e[1]), &self.e[2]),
            &three_inv,
        );
        let r: Vec<E6> = self.e.iter().map(|e| f.sub(e, &t)).collect();
        // (x − r0)(x − r1)(x − r2) = x³ − (Σ r) x² + (Σ r_i r_j) x − r0 r1 r2
        let a = f.add(
            &f.add(&f.mul(&r[0], &r[1]), &f.mul(&r[0], &r[2])),
            &f.mul(&r[1], &r[2]),
        );
        let b = f.neg(&f.mul(&f.mul(&r[0], &r[1]), &r[2]));
        (a, b, t)
    }
    /// The three rational 2-isogenies' targets that keep full 2-torsion:
    /// the kernel `(e_k, 0)`, the other roots `e_i, e_j`, Vélu gives
    /// `y² = x(x² + 2(u + v)x + (u − v)²)` with `u = e_i − e_k`,
    /// `v = e_j − e_k`, whose 2-torsion is rational iff `u·v` is a square.
    pub fn two_isogenies(&self, f: &Fq3, rng: &mut StdRng) -> Vec<Curve2> {
        let mut out = Vec::new();
        for k in 0..3 {
            let u = f.sub(&self.e[(k + 1) % 3], &self.e[k]);
            let v = f.sub(&self.e[(k + 2) % 3], &self.e[k]);
            let Some(s) = f.sqrt(&f.mul(&u, &v), rng) else {
                continue;
            };
            // roots of x² + 2(u+v)x + (u−v)²: −(u+v) ± 2s
            let m = f.neg(&f.add(&u, &v));
            let two_s = f.add(&s, &s);
            let r1 = f.add(&m, &two_s);
            let r2 = f.sub(&m, &two_s);
            if r1 == r2 || r1 == E6::ZERO || r2 == E6::ZERO {
                continue;
            }
            out.push(Curve2 {
                e: [E6::ZERO, r1, r2],
            });
        }
        out
    }
    /// The rational 3-isogenies' targets: for every root `x₀ ∈ F_{q³}` of
    /// the 3-division polynomial of the short model, Vélu's
    /// `a′ = a − 5t`, `b′ = b − 7w` with `t = 2(3x₀² + a)`,
    /// `w = 4y₀² + x₀·t`; the target's 2-torsion is the cubic's roots.
    pub fn three_isogenies(&self, f: &Fq3, rng: &mut StdRng) -> Vec<Curve2> {
        let (a, b, _) = self.short(f);
        let c = |k: u64| f.from_fq(E2([k % f.f.p, 0]));
        // ψ₃ = 3x⁴ + 6a x² + 12b x − a²
        let psi3: P6 = vec![
            f.neg(&f.sq(&a)),
            f.mul(&c(12), &b),
            f.mul(&c(6), &a),
            E6::ZERO,
            c(3),
        ];
        let mut out = Vec::new();
        for x0 in roots_q3(f, &psi3, rng) {
            let y0sq = f.add(&f.add(&f.mul(&x0, &f.sq(&x0)), &f.mul(&a, &x0)), &b);
            let gx = f.add(&f.mul(&c(3), &f.sq(&x0)), &a);
            let t = f.add(&gx, &gx);
            let w = f.add(&f.mul(&c(4), &y0sq), &f.mul(&x0, &t));
            let a2 = f.sub(&a, &f.mul(&c(5), &t));
            let b2 = f.sub(&b, &f.mul(&c(7), &w));
            let cubic: P6 = vec![b2, a2, E6::ZERO, E6::ONE];
            let rs = roots_q3(f, &cubic, rng);
            if rs.len() == 3 {
                out.push(Curve2 {
                    e: [rs[0], rs[1], rs[2]],
                });
            }
        }
        out
    }
}

/// A random curve with full 2-torsion over `F_{q³}`.
pub fn random_curve(f: &Fq3, rng: &mut StdRng) -> Curve2 {
    loop {
        let e = [f.random(rng), f.random(rng), f.random(rng)];
        if e[0] != e[1] && e[1] != e[2] && e[0] != e[2] {
            return Curve2 { e };
        }
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct WalkTrial {
    pub seed: u64,
    pub start_weak: bool,
    pub steps: u64,
    pub distinct: u64,
    pub two_steps: u64,
    pub three_steps: u64,
    pub stuck_restarts: u64,
    pub found: bool,
    /// The walk saw no new `j` for `50·distinct` steps: the component it can
    /// reach by 2- and 3-isogenies is exhausted without a weak curve.
    pub exhausted: bool,
    pub muls: u64,
    pub ms: f64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct WalkReport {
    pub p: u64,
    pub q: u64,
    pub trials: usize,
    pub cap: u64,
    pub use_three: bool,
    pub found: usize,
    pub exhausted: usize,
    pub capped: usize,
    pub mean_steps: f64,
    pub median_steps: f64,
    pub mean_distinct: f64,
    pub mean_muls_per_step: f64,
    pub mean_muls: f64,
    /// `q/3`, the expectation if the weak curves are spread uniformly over the
    /// class and the walk samples it: the registered comparison line.
    pub q_over_3: f64,
    pub cited_q: f64,
    /// Fraction of random full-2-torsion curves that are weak (direct sampling).
    pub weak_fraction: f64,
    pub weak_fraction_times_q: f64,
    pub sampled: u64,
    pub rows: Vec<WalkTrial>,
    pub wall_ms: f64,
}

/// `trials` walks at `p`, each from a fresh random full-2-torsion curve,
/// stepping by uniformly chosen rational 2- (and, with `use_three`, 3-)
/// isogenies until a weak curve or `cap` steps; a curve with no admissible
/// neighbour restarts from a random curve of the same walk (counted).  Also
/// samples `samples` random curves for the weak fraction.
pub fn run_walk(
    p: u64,
    seed: u64,
    trials: usize,
    cap: u64,
    use_three: bool,
    samples: u64,
) -> WalkReport {
    let start = Instant::now();
    let f = Fq3::new(p);
    let q = p * p;
    let mut rows = Vec::new();
    for t in 0..trials as u64 {
        let mut rng = StdRng::seed_from_u64(seed ^ t.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        f.reset_muls();
        let t0 = Instant::now();
        let mut cur = random_curve(&f, &mut rng);
        let mut row = WalkTrial {
            seed: t,
            start_weak: cur.weak_root(&f).is_some(),
            ..Default::default()
        };
        let mut seen: HashSet<E6> = HashSet::new();
        seen.insert(cur.j(&f));
        let mut last_new = 0u64;
        while row.steps < cap {
            if cur.weak_root(&f).is_some() {
                row.found = true;
                break;
            }
            if row.steps > last_new + 50 * seen.len() as u64 + 50 {
                row.exhausted = true;
                break;
            }
            let mut nb: Vec<(Curve2, bool)> = cur
                .two_isogenies(&f, &mut rng)
                .into_iter()
                .map(|c| (c, false))
                .collect();
            if use_three {
                nb.extend(
                    cur.three_isogenies(&f, &mut rng)
                        .into_iter()
                        .map(|c| (c, true)),
                );
            }
            if nb.is_empty() {
                row.stuck_restarts += 1;
                cur = random_curve(&f, &mut rng);
            } else {
                let (next, three) = nb[rng.gen_range(0..nb.len())];
                cur = next;
                if three {
                    row.three_steps += 1;
                } else {
                    row.two_steps += 1;
                }
            }
            row.steps += 1;
            if seen.insert(cur.j(&f)) {
                last_new = row.steps;
            }
        }
        row.distinct = seen.len() as u64;
        row.muls = f.muls();
        row.ms = t0.elapsed().as_secs_f64() * 1e3;
        rows.push(row);
    }
    // direct sampling of the weak fraction
    let mut rng = StdRng::seed_from_u64(seed ^ 0xA11CE);
    let mut weak = 0u64;
    for _ in 0..samples {
        if random_curve(&f, &mut rng).weak_root(&f).is_some() {
            weak += 1;
        }
    }
    let found = rows.iter().filter(|r| r.found).count();
    let exhausted = rows.iter().filter(|r| r.exhausted).count();
    let capped = rows.iter().filter(|r| !r.found && !r.exhausted).count();
    let mut steps: Vec<u64> = rows.iter().filter(|r| r.found).map(|r| r.steps).collect();
    steps.sort_unstable();
    let mean = |v: &[f64]| {
        if v.is_empty() {
            0.0
        } else {
            v.iter().sum::<f64>() / v.len() as f64
        }
    };
    let total_steps: u64 = rows.iter().map(|r| r.steps).sum();
    let total_muls: u64 = rows.iter().map(|r| r.muls).sum();
    WalkReport {
        p,
        q,
        trials,
        cap,
        use_three,
        found,
        exhausted,
        capped,
        mean_steps: mean(&steps.iter().map(|&s| s as f64).collect::<Vec<_>>()),
        median_steps: if steps.is_empty() {
            0.0
        } else {
            steps[steps.len() / 2] as f64
        },
        mean_distinct: mean(&rows.iter().map(|r| r.distinct as f64).collect::<Vec<_>>()),
        mean_muls_per_step: total_muls as f64 / total_steps.max(1) as f64,
        mean_muls: mean(&rows.iter().map(|r| r.muls as f64).collect::<Vec<_>>()),
        q_over_3: q as f64 / 3.0,
        cited_q: q as f64,
        weak_fraction: weak as f64 / samples.max(1) as f64,
        weak_fraction_times_q: weak as f64 / samples.max(1) as f64 * q as f64,
        sampled: samples,
        rows,
        wall_ms: start.elapsed().as_secs_f64() * 1e3,
    }
}

// ── §17: the walk, rebuilt ───────────────────────────────────────────────
//
// `RESEARCH_COVER_DECOMPOSITION_LEDGER.md` §17, registered before this code:
// the curve is carried as its triple of 2-torsion abscissae and never
// normalised; the weak test is three norms and no inversion; a 2-isogeny
// edge costs one character test by the norm down to `F_p` and one square
// root by the odd-degree descent; the 2-isogeny component is enumerated
// breadth first, each curve met once; and a component without a weak curve
// is left by a rational 3-, 5- or 7-isogeny (Vélu's `x`-map on the three
// 2-torsion points gives the image's triple directly).

/// The field with the constants the rebuilt walk keeps: a fixed non-residue
/// of `F_q` raised to the odd part of `q − 1`, for Tonelli–Shanks in `F_q`
/// without a random search per root.
pub struct WalkCtx {
    pub f: Fq3,
    q: u128,
    /// `q − 1 = 2^s · odd`.
    s: u32,
    odd: u128,
    /// `z^{odd}` for a fixed non-residue `z ∈ F_q`.
    z_odd: E2,
}

impl WalkCtx {
    pub fn new(p: u64) -> WalkCtx {
        let f = Fq3::new(p);
        let q = (p as u128) * (p as u128);
        let q1 = q - 1;
        let s = q1.trailing_zeros();
        let odd = q1 >> s;
        // the smallest non-residue of F_q in the enumeration order of F_p × F_p
        let z = (1..)
            .map(|k: u64| E2([k % p, (k / p) % p]))
            .find(|z| *z != E2::ZERO && !f.f.is_square(z))
            .expect("a non-residue of F_q");
        let z_odd = f.f.pow(&z, odd);
        f.reset_muls();
        WalkCtx {
            f,
            q,
            s,
            odd,
            z_odd,
        }
    }
    /// `N_{F_q/F_p}(a) = a₀² − ω a₁²`: three products.
    fn norm_q_p(&self, a: &E2) -> u64 {
        let p = self.f.f.p;
        self.f.f.count_public(3);
        sm(
            mm(a.0[0], a.0[0], p),
            mm(self.f.f.w, mm(a.0[1], a.0[1], p), p),
            p,
        )
    }
    /// Euler's criterion in `F_p`, counted at one product per exponent bit
    /// and a half (square-and-multiply).
    fn is_square_p(&self, x: u64) -> bool {
        let p = self.f.f.p;
        if x == 0 {
            return true;
        }
        let bits = 64 - p.leading_zeros() as u64;
        self.f.f.count_public(bits + bits / 2);
        pow_mod(x, (p - 1) / 2, p) == 1
    }
    /// `w` is a square in `F_{q³}` iff its norm down to `F_p` is a square in
    /// `F_p` (`N(w)^{(p−1)/2} = w^{(q³−1)/2}`): two norms and one Legendre
    /// symbol instead of an exponentiation in `F_{q³}`.
    pub fn is_square_q3(&self, w: &E6) -> bool {
        if *w == E6::ZERO {
            return true;
        }
        let n = norm_q3(&self.f, w);
        debug_assert!(n.in_fq());
        self.is_square_p(self.norm_q_p(&n.0[0]))
    }
    /// Tonelli–Shanks in `F_q` with the cached non-residue; `a` must be a
    /// non-zero square.
    fn sqrt_q(&self, a: &E2) -> E2 {
        let f = &self.f.f;
        let mut m = self.s;
        let mut cc = self.z_odd;
        let mut t = f.pow(a, self.odd);
        let mut r = f.pow(a, self.odd.div_ceil(2));
        while t != E2::ONE {
            let mut i = 0;
            let mut tt = t;
            while tt != E2::ONE {
                tt = f.sq(&tt);
                i += 1;
            }
            let mut b = cc;
            for _ in 0..(m - i - 1) {
                b = f.sq(&b);
            }
            m = i;
            cc = f.sq(&b);
            t = f.mul(&t, &cc);
            r = f.mul(&r, &b);
        }
        debug_assert_eq!(f.sq(&r), *a);
        r
    }
    /// The square root of a square `w ∈ F_{q³}` by the odd-degree descent:
    /// `y = w · σ(w^{(q+1)/2})` has `y² = N_{q³/q}(w) · w`, so
    /// `√w = y / √N(w)`, one exponentiation of `(q + 1)/2` in `F_{q³}` and
    /// one square root in `F_q`.  `None` when `w` is not a square.
    pub fn sqrt_q3(&self, w: &E6) -> Option<E6> {
        if *w == E6::ZERO {
            return Some(E6::ZERO);
        }
        if !self.is_square_q3(w) {
            return None;
        }
        let f = &self.f;
        let n = norm_q3(f, w);
        let rn = self.sqrt_q(&n.0[0]);
        let y = f.mul(&f.sigma(&f.pow(w, self.q.div_ceil(2))), w);
        let r = f.scale(&y, &f.f.inv(&rn));
        debug_assert_eq!(f.sq(&r), *w);
        Some(r)
    }
}

/// A 2-isogeny edge's image, and the root its dual edge needs (the dual's
/// kernel is the image's root `0`, index `0` of the triple, and `u′v′ =
/// (u − v)²`).
#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub curve: Curve2,
    pub dual_root: E6,
    /// The kernel `(e_k, 0)` of the source this edge was taken by (§18's
    /// transport evaluates the isogeny on points by this index).
    pub k: usize,
}

impl Curve2 {
    /// §17's weak test by norms: `N((e₃ − e₁)/(e₂ − e₁)) = 1` is
    /// `N(e₃ − e₁) = N(e₂ − e₁)`, with `N(−1) = −1` kept in the orderings
    /// that reverse a difference.  Agrees with [`Curve2::weak_root`].
    pub fn weak_by_norms(&self, f: &Fq3) -> bool {
        let n01 = norm_q3(f, &f.sub(&self.e[1], &self.e[0]));
        let n02 = norm_q3(f, &f.sub(&self.e[2], &self.e[0]));
        let n12 = norm_q3(f, &f.sub(&self.e[2], &self.e[1]));
        // (e₀; e₁, e₂): N(e₂−e₀) = N(e₁−e₀); (e₁; e₂, e₀): N(e₀−e₁) = N(e₂−e₁);
        // (e₂; e₀, e₁): N(e₁−e₂) = N(e₀−e₂)
        n02 == n01 || n12 == f.neg(&n01) || n12 == n02
    }
    /// The 2-isogeny edges to curves with full 2-torsion, each with the
    /// root its dual needs; `known` is `(kernel index, root)` of the edge
    /// this curve was reached by, whose root is not recomputed.
    pub fn two_edges(&self, ctx: &WalkCtx, known: Option<(usize, E6)>) -> Vec<Edge> {
        let f = &ctx.f;
        let mut out = Vec::new();
        for k in 0..3 {
            let u = f.sub(&self.e[(k + 1) % 3], &self.e[k]);
            let v = f.sub(&self.e[(k + 2) % 3], &self.e[k]);
            let s = match known {
                Some((kk, r)) if kk == k => r,
                _ => match ctx.sqrt_q3(&f.mul(&u, &v)) {
                    Some(s) => s,
                    None => continue,
                },
            };
            let m = f.neg(&f.add(&u, &v));
            let two_s = f.add(&s, &s);
            let r1 = f.add(&m, &two_s);
            let r2 = f.sub(&m, &two_s);
            if r1 == r2 || r1 == E6::ZERO || r2 == E6::ZERO {
                continue;
            }
            out.push(Edge {
                curve: Curve2 {
                    e: [E6::ZERO, r1, r2],
                },
                dual_root: f.sub(&u, &v),
                k,
            });
        }
        out
    }
}

// polynomial helpers over F_{q³} for the division polynomials
fn pscale6(f: &Fq3, a: &[E6], k: &E6) -> P6 {
    let mut out: P6 = a.iter().map(|x| f.mul(x, k)).collect();
    ptrim(&mut out);
    out
}
fn peval6(f: &Fq3, a: &[E6], x: &E6) -> E6 {
    let mut acc = E6::ZERO;
    for c in a.iter().rev() {
        acc = f.add(&f.mul(&acc, x), c);
    }
    acc
}

/// The division polynomials of `y² = x³ + a₂x² + a₄x` that an odd-degree
/// isogeny of degree `3`, `5` or `7` needs: `F = x³ + a₂x² + a₄x`,
/// `ψ₃`, `ψ̃₄ = ψ₄ / ψ₂`, `ψ₅ = 16F²ψ̃₄ − ψ₃³`, `ψ₇ = ψ₅ψ₃³ − 16F²ψ̃₄³`.
struct DivPolys {
    big_f: P6,
    psi3: P6,
    psi4t: P6,
}

impl DivPolys {
    fn new(f: &Fq3, a2: &E6, a4: &E6) -> DivPolys {
        let c = |k: u64| f.from_fq(E2([k % f.f.p, 0]));
        let a4sq = f.sq(a4);
        let big_f: P6 = vec![E6::ZERO, *a4, *a2, E6::ONE];
        // ψ₃ = 3x⁴ + 4a₂x³ + 6a₄x² − a₄²     (b₂ = 4a₂, b₄ = 2a₄, b₆ = 0, b₈ = −a₄²)
        let psi3: P6 = vec![
            f.neg(&a4sq),
            E6::ZERO,
            f.mul(&c(6), a4),
            f.mul(&c(4), a2),
            c(3),
        ];
        // ψ̃₄ = 2x⁶ + b₂x⁵ + 5b₄x⁴ + 10b₆x³ + 10b₈x² + (b₂b₈ − b₄b₆)x + (b₄b₈ − b₆²)
        //     = 2x⁶ + 4a₂x⁵ + 10a₄x⁴ − 10a₄²x² − 4a₂a₄²x − 2a₄³
        let psi4t: P6 = vec![
            f.neg(&f.mul(&c(2), &f.mul(&a4sq, a4))),
            f.neg(&f.mul(&c(4), &f.mul(a2, &a4sq))),
            f.neg(&f.mul(&c(10), &a4sq)),
            E6::ZERO,
            f.mul(&c(10), a4),
            f.mul(&c(4), a2),
            c(2),
        ];
        DivPolys { big_f, psi3, psi4t }
    }
    fn psi(&self, f: &Fq3, ell: u64) -> P6 {
        let c16 = f.from_fq(E2([16 % f.f.p, 0]));
        match ell {
            3 => self.psi3.clone(),
            5 | 7 => {
                let f2 = pmul(f, &self.big_f, &self.big_f);
                let sixteen_f2 = pscale6(f, &f2, &c16);
                let psi3_3 = pmul(f, &self.psi3, &pmul(f, &self.psi3, &self.psi3));
                let psi5 = psub(f, &pmul(f, &sixteen_f2, &self.psi4t), &psi3_3);
                if ell == 5 {
                    psi5
                } else {
                    let psi4t_3 = pmul(f, &self.psi4t, &pmul(f, &self.psi4t, &self.psi4t));
                    psub(f, &pmul(f, &psi5, &psi3_3), &pmul(f, &sixteen_f2, &psi4t_3))
                }
            }
            _ => panic!("ℓ ∈ {{3, 5, 7}}"),
        }
    }
    /// The abscissae `x([k]T)`, `k = 1 … (ℓ−1)/2`, of the kernel generated
    /// by a point `T` with `x(T) = x1`: `x(2T) = x − ψ₃/(4F)`,
    /// `x(3T) = x − 4F·ψ̃₄/ψ₃²`.
    fn kernel_abscissae(&self, f: &Fq3, x1: &E6, ell: u64) -> Vec<E6> {
        let c4 = f.from_fq(E2([4 % f.f.p, 0]));
        let mut out = vec![*x1];
        if ell >= 5 {
            let fx = f.mul(&c4, &peval6(f, &self.big_f, x1));
            let p3 = peval6(f, &self.psi3, x1);
            out.push(f.sub(x1, &f.mul(&p3, &f.inv(&fx))));
            if ell >= 7 {
                let p4 = peval6(f, &self.psi4t, x1);
                out.push(f.sub(x1, &f.mul(&f.mul(&fx, &p4), &f.inv(&f.sq(&p3)))));
            }
        }
        out
    }
}

impl Curve2 {
    /// The images of the rational `ℓ`-isogenies, `ℓ ∈ {3, 5, 7}`, as
    /// triples: for each root of `ψ_ℓ` over `F_{q³}` on the model with
    /// `e₀ ↦ 0`, the kernel abscissae and Vélu's `x`-map on the three
    /// 2-torsion points.  Kernels met twice (through different roots) are
    /// returned once.
    pub fn ell_targets(&self, ctx: &WalkCtx, ell: u64, rng: &mut StdRng) -> Vec<Curve2> {
        self.ell_targets_with_kernels(ctx, ell, rng)
            .into_iter()
            .map(|(c, _)| c)
            .collect()
    }
    /// [`Curve2::ell_targets`] with each target's kernel abscissae on the
    /// model with `e₀ ↦ 0` (§18's transport evaluates the isogeny on points).
    pub fn ell_targets_with_kernels(
        &self,
        ctx: &WalkCtx,
        ell: u64,
        rng: &mut StdRng,
    ) -> Vec<(Curve2, Vec<E6>)> {
        let f = &ctx.f;
        let u = f.sub(&self.e[1], &self.e[0]);
        let v = f.sub(&self.e[2], &self.e[0]);
        let a2 = f.neg(&f.add(&u, &v));
        let a4 = f.mul(&u, &v);
        let dp = DivPolys::new(f, &a2, &a4);
        let psi = dp.psi(f, ell);
        let c = |k: u64| f.from_fq(E2([k % f.f.p, 0]));
        let mut kernels: Vec<Vec<E6>> = Vec::new();
        let mut out = Vec::new();
        for x1 in roots_q3(f, &psi, rng) {
            let mut ker = dp.kernel_abscissae(f, &x1, ell);
            ker.sort();
            if kernels.contains(&ker) {
                continue;
            }
            kernels.push(ker.clone());
            // Vélu: X(P) = x + Σ_k [ t_k/(x − x_k) + u_k/(x − x_k)² ],
            // t_k = 2(3x_k² + 2a₂x_k + a₄), u_k = 4F(x_k)
            let tu: Vec<(E6, E6)> = ker
                .iter()
                .map(|xk| {
                    let gx = f.add(
                        &f.add(&f.mul(&c(3), &f.sq(xk)), &f.mul(&f.mul(&c(2), &a2), xk)),
                        &a4,
                    );
                    (f.add(&gx, &gx), f.mul(&c(4), &peval6(f, &dp.big_f, xk)))
                })
                .collect();
            let image = |x: &E6| -> E6 {
                let mut acc = *x;
                for (xk, (t, uu)) in ker.iter().zip(tu.iter()) {
                    let d = f.inv(&f.sub(x, xk));
                    acc = f.add(&acc, &f.add(&f.mul(t, &d), &f.mul(uu, &f.sq(&d))));
                }
                acc
            };
            let e = [image(&E6::ZERO), image(&u), image(&v)];
            if e[0] != e[1] && e[1] != e[2] && e[0] != e[2] {
                out.push((Curve2 { e }, ker.clone()));
            }
        }
        out
    }
}

// ── §18: evaluating the walk's isogenies on points ───────────────────────
//
// `RESEARCH_COVER_DECOMPOSITION_LEDGER.md` §18.1 steps 1–2: the walk reaches
// a weak curve; to carry a logarithm over, each step's isogeny is evaluated
// on the base point and the target, on the Weierstrass model the walk uses,
// `y² = (x − e₀)(x − e₁)(x − e₂)`.

/// An affine point `(x, y)` of `y² = (x − e₀)(x − e₁)(x − e₂)`.
pub type Pt2 = (E6, E6);

impl Curve2 {
    /// Whether `(x, y)` lies on this curve.
    pub fn on_curve(&self, f: &Fq3, p: &Pt2) -> bool {
        let rhs = f.mul(
            &f.mul(&f.sub(&p.0, &self.e[0]), &f.sub(&p.0, &self.e[1])),
            &f.sub(&p.0, &self.e[2]),
        );
        f.sq(&p.1) == rhs
    }
}

/// The 2-isogeny with kernel `(e_k, 0)` on an affine point.  On the model
/// translated by `e_k`, `y² = x(x² + a x + b)` with `b = u v`, Vélu gives
/// `(x, y) ↦ (y²/x², y(b − x²)/x²)` onto the triple [`Curve2::two_edges`]
/// returns (`Y² = X(X² + 2(u+v)X + (u−v)²)`).  The two-torsion point in the
/// kernel maps to infinity; §18 never transports such a point (`G`, `Q` have
/// odd order `ℓ`).
pub fn map_two(f: &Fq3, from: &Curve2, k: usize, p: &Pt2) -> Pt2 {
    let u = f.sub(&from.e[(k + 1) % 3], &from.e[k]);
    let v = f.sub(&from.e[(k + 2) % 3], &from.e[k]);
    let x = f.sub(&p.0, &from.e[k]);
    let xi2 = f.sq(&f.inv(&x));
    let b = f.mul(&u, &v);
    let xx = f.mul(&f.sq(&p.1), &xi2);
    let yy = f.mul(&f.mul(&p.1, &f.sub(&b, &f.sq(&x))), &xi2);
    (xx, yy)
}

/// Vélu's odd-degree isogeny on an affine point, on the model with `e₀ ↦ 0`
/// (`a₁ = a₃ = 0`), from the kernel abscissae [`Curve2::ell_targets_with_kernels`]
/// records: `X = x + Σ_Q (v_Q/(x − x_Q) + u_Q/(x − x_Q)²)`,
/// `Y = y·(1 − Σ_Q (2 u_Q/(x − x_Q)³ + v_Q/(x − x_Q)²))`, with
/// `u_Q = 4 F(x_Q)`, `v_Q = 2(3 x_Q² + 2 a₂ x_Q + a₄)`.  The codomain triple
/// is the one the walk recorded for this step.
pub fn map_odd(f: &Fq3, from: &Curve2, ker: &[E6], p: &Pt2) -> Pt2 {
    let c = |k: u64| f.from_fq(E2([k % f.f.p, 0]));
    let u = f.sub(&from.e[1], &from.e[0]);
    let v = f.sub(&from.e[2], &from.e[0]);
    let a2 = f.neg(&f.add(&u, &v));
    let a4 = f.mul(&u, &v);
    // F(x) = x³ + a₂ x² + a₄ x on the e₀ ↦ 0 model
    let big_f = |x: &E6| f.mul(&f.add(&f.add(&f.sq(x), &f.mul(&a2, x)), &a4), x);
    let x = f.sub(&p.0, &from.e[0]);
    let mut xx = x;
    let mut ysum = E6::ZERO;
    for xk in ker {
        let gx = f.add(
            &f.add(&f.mul(&c(3), &f.sq(xk)), &f.mul(&f.mul(&c(2), &a2), xk)),
            &a4,
        );
        let vq = f.add(&gx, &gx);
        let uq = f.mul(&c(4), &big_f(xk));
        let d = f.inv(&f.sub(&x, xk));
        let d2 = f.sq(&d);
        xx = f.add(&xx, &f.add(&f.mul(&vq, &d), &f.mul(&uq, &d2)));
        ysum = f.add(
            &ysum,
            &f.add(
                &f.mul(&f.mul(&c(2), &uq), &f.mul(&d2, &d)),
                &f.mul(&vq, &d2),
            ),
        );
    }
    (xx, f.mul(&p.1, &f.sub(&E6::ONE, &ysum)))
}

/// One recorded step of a path: the source curve and the data to replay its
/// isogeny on a point.
#[derive(Clone, Debug)]
pub enum Step {
    /// The 2-isogeny of `from` with kernel `(e_k, 0)`.
    Two { from: Curve2, k: usize },
    /// The odd-degree isogeny of `from` with the given kernel abscissae.
    Odd {
        from: Curve2,
        ell: u64,
        ker: Vec<E6>,
    },
}

/// Apply one recorded step to an affine point.
pub fn map_step(f: &Fq3, step: &Step, p: &Pt2) -> Pt2 {
    match step {
        Step::Two { from, k } => map_two(f, from, *k, p),
        Step::Odd { from, ker, .. } => map_odd(f, from, ker, p),
    }
}

/// §18.1 step 3: move a weak curve onto the cover's form `y² = x(x−α)(x−σα)`
/// and carry two points over.  `w.weak_root` gives the ordering `(e₁; e₂, e₃)`
/// with `N((e₃−e₁)/(e₂−e₁)) = 1`; translating `e₁ ↦ 0` and scaling by `v`
/// with `(e₂−e₁)/v = α`, `(e₃−e₁)/v = σ(α)` needs `α` with `σ(α) = c·α`,
/// `c = (e₃−e₁)/(e₂−e₁)`.  Since `N(c) = 1`, Hilbert 90 gives such an `α`
/// explicitly: `α = θ + c⁻¹σ(θ) + (c·σ(c))⁻¹σ²(θ)` for a `θ` with `α ∉ F_q`.
/// `v` must be a square in `F_{q³}`; if not, `α` is scaled by a non-square of
/// `F_q` (an `F_q`-isomorphism that keeps `σ(α)=cα`), which flips it.  The
/// point map is `(x, y) ↦ ((x−e₁)/v, y/√(v³))`.  `None` if the cover refuses
/// the resulting `α`.
pub fn model_change(ctx: &WalkCtx, w: &Curve2, pts: &[Pt2]) -> Option<(E6, Vec<Pt2>)> {
    let f = &ctx.f;
    let i1 = w.weak_root(f)?;
    let e1 = w.e[i1];
    let e2 = w.e[(i1 + 1) % 3];
    let e3 = w.e[(i1 + 2) % 3];
    let c = f.mul(&f.sub(&e3, &e1), &f.inv(&f.sub(&e2, &e1)));
    // a non-square of F_q (every F_p element is a square in F_q = F_p², so
    // this needs a nonzero imaginary part), still a non-square in F_{q³}
    // because [F_{q³} : F_q] = 3 is odd
    let p = f.f.p;
    let s_nonsq = (1..p * p)
        .map(|k| E2([k % p, k / p]))
        .find(|s| !f.f.is_square(s))
        .map(|s| f.from_fq(s))?;
    let mut rng = StdRng::seed_from_u64(0xA7E_11 ^ (e1.0[0].0[0]).wrapping_mul(0x9E3779B1));
    for _ in 0..64 {
        let theta = f.random(&mut rng);
        let ci = f.inv(&c);
        let csc = f.inv(&f.mul(&c, &f.sigma(&c)));
        let st = f.sigma(&theta);
        let s2t = f.sigma(&st);
        let mut alpha = f.add(&f.add(&theta, &f.mul(&ci, &st)), &f.mul(&csc, &s2t));
        if alpha == E6::ZERO || alpha.in_fq() {
            continue;
        }
        // σ(α) = c·α by construction; v = (e₂−e₁)/α, square-fixed
        let mut v = f.mul(&f.sub(&e2, &e1), &f.inv(&alpha));
        if !ctx.is_square_q3(&v) {
            alpha = f.mul(&alpha, &s_nonsq);
            v = f.mul(&v, &f.inv(&s_nonsq));
        }
        let Some(sqrt_v) = ctx.sqrt_q3(&v) else {
            continue;
        };
        debug_assert_eq!(f.sigma(&alpha), f.mul(&c, &alpha));
        let w_scale = f.mul(&v, &sqrt_v); // √(v³)
        let inv_v = f.inv(&v);
        let inv_w = f.inv(&w_scale);
        if super::jv_cover::Cover::new(f, &alpha).is_none() {
            continue;
        }
        let mapped: Vec<Pt2> = pts
            .iter()
            .map(|(x, y)| (f.mul(&f.sub(x, &e1), &inv_v), f.mul(y, &inv_w)))
            .collect();
        return Some((alpha, mapped));
    }
    None
}

fn isqrt_u128(n: u128) -> u128 {
    if n < 2 {
        return n;
    }
    let mut x = (n as f64).sqrt() as u128;
    while x * x > n {
        x -= 1;
    }
    while (x + 1) * (x + 1) <= n {
        x += 1;
    }
    x
}

/// `#E(F_{q³})` of the curve `y² = (x − e₀)(x − e₁)(x − e₂)` by baby-step
/// giant-step on two random points (the first's multiple in the Hasse
/// interval, confirmed on the second); counted through the field's counter
/// like everything else.  The attacker knows the order of the curve it is
/// given; the walk pays for it once here so that its cost is in the table.
pub fn curve_order(f: &Fq3, c: &Curve2, rng: &mut StdRng) -> u128 {
    let u = f.sub(&c.e[1], &c.e[0]);
    let v = f.sub(&c.e[2], &c.e[0]);
    let ec = EllE::from_a2_a4(f, f.neg(&f.add(&u, &v)), f.mul(&u, &v));
    ell_order(f, &ec, rng)
}

/// `#E(F_{q³})` of `y² = x³ + a₂x² + a₄x` by baby-step giant-step on two
/// random points (the first's order located in the Hasse interval, confirmed
/// on the second).
pub fn ell_order(f: &Fq3, ec: &EllE, rng: &mut StdRng) -> u128 {
    let p = f.f.p as u128;
    let q3 = p.pow(6);
    let two_sqrt = 2 * isqrt_u128(q3) + 2;
    let lo = q3 + 1 - two_sqrt;
    let width = 2 * two_sqrt;
    let steps = isqrt_u128(width) + 1;
    let order_of = |pt: &PtE6| -> Option<u128> {
        let mut table: HashMap<PtE6, u128> = HashMap::new();
        let mut jp = PtE6::INF;
        for j in 0..steps {
            if !jp.inf {
                table.entry(jp).or_insert(j);
            }
            jp = ec.add(&jp, pt);
        }
        let giant = ec.mul_u128(pt, steps);
        let mut t = ec.mul_u128(pt, lo);
        let mut i = 0u128;
        while i * steps <= width + steps {
            if t.inf {
                return Some(lo + i * steps);
            }
            if let Some(&j) = table.get(&ec.neg(&t)) {
                return Some(lo + i * steps + j);
            }
            t = ec.add(&t, &giant);
            i += 1;
        }
        None
    };
    loop {
        let p1 = ec.random_point(rng);
        let Some(m) = order_of(&p1) else { continue };
        let p2 = ec.random_point(rng);
        if ec.mul_u128(&p2, m).inf {
            return m;
        }
    }
}

/// Kronecker symbol `(D / ℓ)` for an odd prime `ℓ`: `0`, `1` or `−1`.
fn kronecker(d: i128, ell: u64) -> i8 {
    let l = ell as i128;
    let r = ((d % l) + l) % l;
    if r == 0 {
        return 0;
    }
    // Euler's criterion in F_ℓ
    let mut acc: i128 = 1;
    let mut base = r;
    let mut e = (l - 1) / 2;
    while e > 0 {
        if e & 1 == 1 {
            acc = acc * base % l;
        }
        base = base * base % l;
        e >>= 1;
    }
    if acc == 1 {
        1
    } else {
        -1
    }
}

/// One walk of §17.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct Walk2Trial {
    pub seed: u64,
    pub start_weak: bool,
    /// `#E(F_{q³})` of the start curve (the whole class's), its trace
    /// `t = q³ + 1 − #E`, and `(D / ℓ)` for `ℓ = 3, 5, 7` with
    /// `D = t² − 4q³`: a degree with `(D/ℓ) = −1` and `ℓ ∤ D` has no rational
    /// `ℓ`-isogeny anywhere in the class and is not tried.
    pub order: u128,
    pub trace: i128,
    pub kronecker: [i8; 3],
    pub degrees_used: Vec<u64>,
    pub found: bool,
    /// Met `cap` distinct curves without a weak one.
    pub capped: bool,
    /// Distinct curves met (the visited set's size).
    pub curves: u64,
    pub components: u64,
    /// Size of the start curve's own 2-isogeny component, and whether it
    /// held a weak curve (P2).
    pub first_component: u64,
    pub first_component_weak: bool,
    /// Jumps taken by degree 3, 5, 7.
    pub jumps: [u64; 3],
    /// Jumps whose image was a curve already met.
    pub wasted_jumps: u64,
    /// No unmet `ℓ`-neighbour was found from the curves tried: the set the
    /// walk can reach is (heuristically, or exactly in closure mode) closed
    /// and holds no weak curve.  The walk stops; it does not restart, since
    /// a restart would change the instance.
    pub exhausted: bool,
    /// Closure mode: the number of 2-isogeny components the exact
    /// `{2} ∪ jumps` closure of the start curve was found to contain.
    pub closure_components: u64,
    pub muls_order: u64,
    pub muls_enumerate: u64,
    pub muls_jumps: u64,
    pub muls: u64,
    pub ms: f64,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct Walk2Report {
    pub p: u64,
    pub q: u64,
    pub trials: usize,
    pub cap_curves: u64,
    pub jump_degrees: Vec<u64>,
    pub closure_mode: bool,
    pub from_weak_class: bool,
    pub weak_class_moves: usize,
    pub found: usize,
    pub capped: usize,
    pub exhausted: usize,
    pub start_weak: usize,
    /// `found / (trials − start_weak)`: the chance that a walk from a
    /// random full-2-torsion curve reaches a weak one within the budget.
    pub success_fraction: f64,
    /// Over the walks that found a weak curve and did not start on one.
    pub median_curves: f64,
    pub mean_curves: f64,
    pub mean_first_component: f64,
    pub frac_first_component_weak: f64,
    pub mean_components: f64,
    pub total_jumps: [u64; 3],
    pub total_wasted_jumps: u64,
    /// `F_p` multiplications per distinct curve met in the enumeration.
    pub c_curve: f64,
    /// `F_p` multiplications per jump attempted (wasted ones included).
    pub c_jump: f64,
    /// `F_p` multiplications of the point count, mean per walk.
    pub c_order: f64,
    /// Mean cost of a walk that found a weak curve (not starting on one),
    /// every multiplication inside; and the mean over every walk, the
    /// capped and exhausted ones included.
    pub mean_muls: f64,
    pub mean_muls_all: f64,
    /// `1.3 · (p³/2) · 331`, rho on the subgroup the route attacks at `p`.
    pub rho_at_p: f64,
    /// `mean_muls / rho_at_p` and `mean_muls_all / rho_at_p`.
    pub walk_over_rho: f64,
    pub walk_over_rho_all: f64,
    pub q_over_3: f64,
    pub weak_fraction: f64,
    pub weak_fraction_times_q: f64,
    pub sampled: u64,
    pub rows: Vec<Walk2Trial>,
    pub wall_ms: f64,
}

/// Enumerate the 2-isogeny component of `start` among full-2-torsion
/// curves (breadth first, keyed by `j`, `seen` shared across the walk),
/// stopping at a weak curve or when `seen` reaches `cap`.  Returns the
/// curves met in this component, the weak curve if any, and whether the
/// cap stopped it.  `start` must already be in `seen`.
fn enumerate_component(
    ctx: &WalkCtx,
    start: Curve2,
    seen: &mut HashSet<E6>,
    cap: u64,
) -> (Vec<Curve2>, Option<Curve2>, bool) {
    let f = &ctx.f;
    let mut comp = vec![start];
    if start.weak_by_norms(f) {
        return (comp, Some(start), false);
    }
    let mut queue: VecDeque<(Curve2, Option<(usize, E6)>)> = VecDeque::new();
    queue.push_back((start, None));
    while let Some((cur, known)) = queue.pop_front() {
        for e in cur.two_edges(ctx, known) {
            if seen.insert(e.curve.j(f)) {
                comp.push(e.curve);
                if e.curve.weak_by_norms(f) {
                    return (comp, Some(e.curve), false);
                }
                if seen.len() as u64 >= cap {
                    return (comp, None, true);
                }
                queue.push_back((e.curve, Some((0, e.dual_root))));
            }
        }
    }
    (comp, None, false)
}

/// A random weak curve moved away from the weak locus: `≥ min_moves` random moves
/// (a 2-edge, or an `ℓ`-jump for `ℓ` in `jumps` when the curve has one)
/// and then further until the curve is not weak.  Same class as the weak
/// curve it started from.
fn weak_class_start(ctx: &WalkCtx, jumps: &[u64], min_moves: usize, rng: &mut StdRng) -> Curve2 {
    let f = &ctx.f;
    loop {
        let mut cur = loop {
            let rho = f.from_fq(f.f.random(rng));
            let alpha = f.random(rng);
            if alpha.in_fq() {
                continue;
            }
            let sa = f.sigma(&alpha);
            if alpha != rho && sa != rho {
                break Curve2 {
                    e: [rho, alpha, sa],
                };
            }
        };
        let mut moves = 0;
        let mut tries = 0;
        while moves < min_moves || cur.weak_by_norms(f) {
            tries += 1;
            if tries > 25 * min_moves {
                break;
            }
            let pick = rng.gen_range(0..4);
            if pick == 0 && !jumps.is_empty() {
                let ell = jumps[rng.gen_range(0..jumps.len())];
                let t = cur.ell_targets(ctx, ell, rng);
                if !t.is_empty() {
                    cur = t[rng.gen_range(0..t.len())];
                    moves += 1;
                }
            } else {
                let e = cur.two_edges(ctx, None);
                if !e.is_empty() {
                    cur = e[rng.gen_range(0..e.len())].curve;
                    moves += 1;
                }
            }
        }
        if !cur.weak_by_norms(f) {
            return cur;
        }
    }
}

/// `trials` walks of §17 at `p`: enumerate the 2-isogeny component, jump by
/// an `ℓ`-isogeny (`jumps`, in order of preference, the degrees inert in
/// the class's order skipped, one source curve per degree and component)
/// when it holds no weak curve, breadth first over components, until one
/// is met, `cap_curves` distinct curves have been, or the reachable set is
/// closed (`exhausted`).
pub fn run_walk2(
    p: u64,
    seed: u64,
    trials: usize,
    cap_curves: u64,
    jumps: &[u64],
    samples: u64,
    closure_mode: bool,
    from_weak_class: bool,
    weak_class_moves: usize,
) -> Walk2Report {
    let start = Instant::now();
    let ctx = WalkCtx::new(p);
    let f = &ctx.f;
    let q = p * p;
    let jidx = |ell: u64| match ell {
        3 => 0,
        5 => 1,
        7 => 2,
        _ => panic!("ℓ ∈ {{3, 5, 7}}"),
    };
    let mut rows = Vec::new();
    for t in 0..trials as u64 {
        let mut rng = StdRng::seed_from_u64(seed ^ t.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        f.reset_muls();
        let t0 = Instant::now();
        let mut row = Walk2Trial {
            seed: t,
            ..Default::default()
        };
        let mut seen: HashSet<E6> = HashSet::new();
        let mut cur = if from_weak_class {
            // a non-weak curve of a class known to hold a weak curve: from a
            // random weak curve, a random path of 2-steps and ℓ-jumps until
            // a non-weak curve at least `weak_class_moves` moves away; its cost is not
            // the walk's and is subtracted below
            let m0 = f.muls();
            let c = weak_class_start(&ctx, jumps, weak_class_moves, &mut rng);
            f.reset_muls();
            let _ = m0;
            c
        } else {
            random_curve(f, &mut rng)
        };
        seen.insert(cur.j(f));
        row.start_weak = cur.weak_by_norms(f);
        // the class's order, trace and discriminant: degrees inert in the
        // order are never tried
        let m_ord = f.muls();
        row.order = curve_order(f, &cur, &mut rng);
        row.muls_order = f.muls() - m_ord;
        let q3 = (p as i128).pow(6);
        row.trace = q3 + 1 - row.order as i128;
        let disc = row.trace * row.trace - 4 * q3;
        row.kronecker = [kronecker(disc, 3), kronecker(disc, 5), kronecker(disc, 7)];
        let degrees: Vec<u64> = jumps
            .iter()
            .copied()
            .filter(|&ell| row.kronecker[jidx(ell)] != -1)
            .collect();
        row.degrees_used = degrees.clone();
        // Breadth first over components, lazily: `pending` holds landing
        // curves of components met but not yet enumerated; `partial` holds
        // components whose (source curve, degree) pairs are not all tried
        // yet.  A component is expanded only until it yields a new
        // neighbour (closure mode: always fully), and a spent `pending`
        // resumes the most recent partial component, so the walk gives up
        // only when the reachable set is closed.
        let mut pending: VecDeque<Curve2> = VecDeque::new();
        let mut partial: Vec<(Vec<Curve2>, Vec<(usize, u64)>, usize)> = Vec::new();
        loop {
            let m0 = f.muls();
            let (comp, weak, capped) = enumerate_component(&ctx, cur, &mut seen, cap_curves);
            row.muls_enumerate += f.muls() - m0;
            row.components += 1;
            if row.components == 1 {
                row.first_component = comp.len() as u64;
                row.first_component_weak = weak.is_some();
            }
            if weak.is_some() {
                row.found = true;
                break;
            }
            if capped {
                row.capped = true;
                break;
            }
            // the (source, degree) pairs of this component, degree-major,
            // sources in a uniform order (all of them in closure mode, up to
            // eight otherwise)
            let mut order: Vec<usize> = (0..comp.len()).collect();
            for i in (1..order.len()).rev() {
                let k = rng.gen_range(0..=i);
                order.swap(i, k);
            }
            if !closure_mode {
                order.truncate(1);
            }
            let plan: Vec<(usize, u64)> = degrees
                .iter()
                .flat_map(|&ell| order.iter().map(move |&src| (src, ell)))
                .collect();
            partial.push((comp, plan, 0));
            // leave for the next component
            let m1 = f.muls();
            let mut landed: Option<Curve2> = None;
            while landed.is_none() {
                if !closure_mode {
                    if let Some(c) = pending.pop_front() {
                        landed = Some(c);
                        break;
                    }
                }
                let Some((pcomp, plan, cursor)) = partial.last_mut() else {
                    break;
                };
                if *cursor >= plan.len() {
                    partial.pop();
                    if closure_mode {
                        landed = pending.pop_front();
                        if landed.is_some() {
                            break;
                        }
                    }
                    continue;
                }
                let (src, ell) = plan[*cursor];
                *cursor += 1;
                let mut targets = pcomp[src].ell_targets(&ctx, ell, &mut rng);
                for i in (1..targets.len()).rev() {
                    let k = rng.gen_range(0..=i);
                    targets.swap(i, k);
                }
                for tgt in targets {
                    let jj = tgt.j(f);
                    if seen.contains(&jj) {
                        row.wasted_jumps += 1;
                        continue;
                    }
                    seen.insert(jj);
                    row.jumps[jidx(ell)] += 1;
                    pending.push_back(tgt);
                }
                if closure_mode && *cursor >= plan.len() {
                    partial.pop();
                    landed = pending.pop_front();
                    if landed.is_some() {
                        break;
                    }
                    if partial.is_empty() {
                        break;
                    }
                }
            }
            row.muls_jumps += f.muls() - m1;
            match landed {
                Some(c) => cur = c,
                None => {
                    row.exhausted = true;
                    row.closure_components = row.components;
                    break;
                }
            }
        }
        row.curves = seen.len() as u64;
        row.muls = f.muls();
        row.ms = t0.elapsed().as_secs_f64() * 1e3;
        rows.push(row);
    }
    // direct sampling of the weak fraction, by the norm test
    let mut rng = StdRng::seed_from_u64(seed ^ 0xA11CE);
    let mut weak = 0u64;
    for _ in 0..samples {
        if random_curve(f, &mut rng).weak_by_norms(f) {
            weak += 1;
        }
    }
    let mean = |v: &[f64]| {
        if v.is_empty() {
            0.0
        } else {
            v.iter().sum::<f64>() / v.len() as f64
        }
    };
    let found = rows.iter().filter(|r| r.found).count();
    let capped = rows.iter().filter(|r| r.capped).count();
    let exhausted = rows.iter().filter(|r| r.exhausted).count();
    let start_weak = rows.iter().filter(|r| r.start_weak).count();
    let proper: Vec<&Walk2Trial> = rows.iter().filter(|r| r.found && !r.start_weak).collect();
    let mut curves: Vec<f64> = proper.iter().map(|r| r.curves as f64).collect();
    curves.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let total_curves: u64 = rows.iter().map(|r| r.curves).sum();
    let total_enum: u64 = rows.iter().map(|r| r.muls_enumerate).sum();
    let total_jump_muls: u64 = rows.iter().map(|r| r.muls_jumps).sum();
    let mut total_jumps = [0u64; 3];
    for r in &rows {
        for i in 0..3 {
            total_jumps[i] += r.jumps[i];
        }
    }
    let total_wasted: u64 = rows.iter().map(|r| r.wasted_jumps).sum();
    let jump_attempts = total_jumps.iter().sum::<u64>() + total_wasted;
    let rho_at_p = 1.3 * ((p as f64).powi(6) / 4.0).sqrt() * 331.0;
    let mean_muls = mean(&proper.iter().map(|r| r.muls as f64).collect::<Vec<_>>());
    let mean_muls_all = mean(&rows.iter().map(|r| r.muls as f64).collect::<Vec<_>>());
    Walk2Report {
        p,
        q,
        trials,
        cap_curves,
        jump_degrees: jumps.to_vec(),
        closure_mode,
        from_weak_class,
        weak_class_moves,
        found,
        capped,
        exhausted,
        start_weak,
        success_fraction: proper.len() as f64 / (trials - start_weak).max(1) as f64,
        median_curves: if curves.is_empty() {
            0.0
        } else {
            curves[curves.len() / 2]
        },
        mean_curves: mean(&curves),
        mean_first_component: mean(
            &rows
                .iter()
                .map(|r| r.first_component as f64)
                .collect::<Vec<_>>(),
        ),
        frac_first_component_weak: rows.iter().filter(|r| r.first_component_weak).count() as f64
            / rows.len().max(1) as f64,
        mean_components: mean(&rows.iter().map(|r| r.components as f64).collect::<Vec<_>>()),
        total_jumps,
        total_wasted_jumps: total_wasted,
        c_curve: total_enum as f64 / total_curves.max(1) as f64,
        c_jump: total_jump_muls as f64 / jump_attempts.max(1) as f64,
        c_order: mean(&rows.iter().map(|r| r.muls_order as f64).collect::<Vec<_>>()),
        mean_muls,
        mean_muls_all,
        rho_at_p,
        walk_over_rho: mean_muls / rho_at_p,
        walk_over_rho_all: mean_muls_all / rho_at_p,
        q_over_3: q as f64 / 3.0,
        weak_fraction: weak as f64 / samples.max(1) as f64,
        weak_fraction_times_q: weak as f64 / samples.max(1) as f64 * q as f64,
        sampled: samples,
        rows,
        wall_ms: start.elapsed().as_secs_f64() * 1e3,
    }
}

/// §18.1 steps 1–2: a breadth-first walk to a weak curve from `start`, with
/// each discovered curve's producing step recorded, so the path from `start`
/// to the weak curve can be replayed on points.  Degrees whose Kronecker
/// symbol `(t²−4q³ / ℓ)` is `−1` are inert in the class and skipped; `trace`
/// is read from the public order, not point-counted.  Stops at a weak curve
/// or when `cap` distinct curves have been met.  Charges nothing itself; the
/// caller brackets it with the field counter.
pub fn walk_to_weak_record(
    ctx: &WalkCtx,
    start: Curve2,
    jumps: &[u64],
    cap: u64,
    trace: i128,
    rng: &mut StdRng,
) -> Option<(Curve2, Vec<Step>)> {
    let f = &ctx.f;
    let p3 = (f.f.p as i128).pow(6);
    let disc = trace * trace - 4 * p3;
    let degrees: Vec<u64> = jumps
        .iter()
        .copied()
        .filter(|&ell| kronecker(disc, ell) != -1)
        .collect();
    let j0 = start.j(f);
    if start.weak_by_norms(f) {
        return Some((start, Vec::new()));
    }
    let mut seen: HashSet<E6> = HashSet::new();
    seen.insert(j0);
    let mut parent: HashMap<E6, (E6, Step)> = HashMap::new();
    let mut curve_of: HashMap<E6, Curve2> = HashMap::new();
    curve_of.insert(j0, start);
    let mut queue: VecDeque<(Curve2, Option<(usize, E6)>)> = VecDeque::new();
    queue.push_back((start, None));
    let reconstruct =
        |weak: Curve2, weak_j: E6, parent: &HashMap<E6, (E6, Step)>| -> (Curve2, Vec<Step>) {
            let mut steps = Vec::new();
            let mut j = weak_j;
            while let Some((pj, step)) = parent.get(&j) {
                steps.push(step.clone());
                j = *pj;
            }
            steps.reverse();
            (weak, steps)
        };
    while let Some((cur, known)) = queue.pop_front() {
        let cj = cur.j(f);
        // 2-isogeny edges
        for e in cur.two_edges(ctx, known) {
            let tj = e.curve.j(f);
            if !seen.insert(tj) {
                continue;
            }
            parent.insert(tj, (cj, Step::Two { from: cur, k: e.k }));
            curve_of.insert(tj, e.curve);
            if e.curve.weak_by_norms(f) {
                return Some(reconstruct(e.curve, tj, &parent));
            }
            if seen.len() as u64 >= cap {
                return None;
            }
            queue.push_back((e.curve, Some((0, e.dual_root))));
        }
        // odd-degree jumps
        for &ell in &degrees {
            for (tgt, ker) in cur.ell_targets_with_kernels(ctx, ell, rng) {
                let tj = tgt.j(f);
                if !seen.insert(tj) {
                    continue;
                }
                parent.insert(
                    tj,
                    (
                        cj,
                        Step::Odd {
                            from: cur,
                            ell,
                            ker,
                        },
                    ),
                );
                curve_of.insert(tj, tgt);
                if tgt.weak_by_norms(f) {
                    return Some(reconstruct(tgt, tj, &parent));
                }
                if seen.len() as u64 >= cap {
                    return None;
                }
                queue.push_back((tgt, None));
            }
        }
    }
    None
}

impl Curve2 {
    /// `[k]·(x, y)` on `y² = (x−e₀)(x−e₁)(x−e₂)`, through the `e₀ ↦ 0`
    /// Weierstrass model and `EllE`.  `None` for the identity.
    pub fn scalar_mul(&self, f: &Fq3, pt: &Pt2, k: u128) -> Option<Pt2> {
        let u = f.sub(&self.e[1], &self.e[0]);
        let v = f.sub(&self.e[2], &self.e[0]);
        let ec = EllE::from_a2_a4(f, f.neg(&f.add(&u, &v)), f.mul(&u, &v));
        let p6 = PtE6 {
            x: f.sub(&pt.0, &self.e[0]),
            y: pt.1,
            inf: false,
        };
        let r = ec.mul_u128(&p6, k);
        if r.inf {
            None
        } else {
            Some((f.add(&r.x, &self.e[0]), r.y))
        }
    }
}

/// §18.1: the whole route from a curve that is not weak.  A weak instance is
/// built (`generate_spec`), moved `moves` random isogeny steps until it is
/// not weak (the challenge `C`, with `G`, `Q` transported onto it — the
/// instance's construction, uncharged), then attacked: the walk back to a
/// weak curve, the points transported along it, the model change, the sieved
/// route, and verification `[d]·G = Q` on `C`.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct EndToEndReport {
    pub p: u64,
    pub seed: u64,
    pub moves: usize,
    pub l: u64,
    pub bits: f64,
    /// The challenge curve was not weak when handed over.
    pub challenge_non_weak: bool,
    /// The walk found a weak curve; its path length in steps.
    pub walk_found: bool,
    pub path_steps: usize,
    pub path_two: usize,
    pub path_odd: usize,
    /// Charged `F_p` multiplications, by phase.
    pub walk_muls: u64,
    pub transport_muls: u64,
    pub model_muls: u64,
    pub route_muls: u64,
    pub total_muls: u64,
    /// The share of the total spent reaching and reshaping the weak curve.
    pub reach_share: f64,
    pub c_add_e: f64,
    /// The route solved and the recovered scalar verified `[d]G = Q` on `C`.
    pub route_solved: bool,
    pub route_correct: bool,
    pub verified_on_challenge: bool,
    /// `S` end to end, and against the pooled rho reference.
    pub s: f64,
    pub rho_s_ref: f64,
    pub s_over_rho: f64,
    /// The same route with no walk, on the weak curve of the same seed
    /// (paired): its `S/rho`, and the end-to-end total over its total.
    pub route_only_s_over_rho: f64,
    pub e2e_over_route_only: f64,
    pub stage: String,
    pub wall_ms: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn run_end_to_end(
    p: u64,
    seed: u64,
    moves: usize,
    jumps: &[u64],
    cap_mult: u64,
    rho_s_ref: f64,
    margin: f64,
    stop: Option<usize>,
) -> EndToEndReport {
    use super::jv_cover::generate_spec;
    use super::jv_sieve::{run_cover_sieve_dlp, run_cover_sieve_dlp_on};
    let start = Instant::now();
    let ctx = WalkCtx::new(p);
    let f = &ctx.f;
    let spec0 = generate_spec(p, seed);
    let l = spec0.l;
    let w0 = Curve2 {
        e: [E6::ZERO, spec0.alpha, f.sigma(&spec0.alpha)],
    };
    let trace = (p as i128).pow(6) + 1 - 4 * l as i128;
    let mut rng = StdRng::seed_from_u64(seed ^ 0xE2E_18);
    let mut rep = EndToEndReport {
        p,
        seed,
        moves,
        l,
        bits: (l as f64).log2(),
        c_add_e: 0.0,
        rho_s_ref,
        stage: String::from("ok"),
        ..Default::default()
    };

    // ── build the challenge C: move away from weak, transporting G, Q
    let mut cur = w0;
    let mut g = (spec0.g.x, spec0.g.y);
    let mut q = (spec0.q.x, spec0.q.y);
    let mut done = false;
    for i in 0..(moves + 200) {
        // candidate moves: 2-edges and ℓ-targets, each with its codomain
        let mut cands: Vec<(Step, Curve2)> = cur
            .two_edges(&ctx, None)
            .into_iter()
            .map(|e| (Step::Two { from: cur, k: e.k }, e.curve))
            .collect();
        for &ell in jumps {
            for (tgt, ker) in cur.ell_targets_with_kernels(&ctx, ell, &mut rng) {
                cands.push((
                    Step::Odd {
                        from: cur,
                        ell,
                        ker,
                    },
                    tgt,
                ));
            }
        }
        if cands.is_empty() {
            break;
        }
        let (step, next) = cands[rng.gen_range(0..cands.len())].clone();
        g = map_step(f, &step, &g);
        q = map_step(f, &step, &q);
        cur = next;
        if i + 1 >= moves && !cur.weak_by_norms(f) {
            done = true;
            break;
        }
    }
    if !done && cur.weak_by_norms(f) {
        // could not leave the weak locus; report the failure honestly
        rep.stage = String::from("challenge_weak");
        rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
        return rep;
    }
    rep.challenge_non_weak = !cur.weak_by_norms(f);
    let challenge = cur;
    let gc = g;
    let qc = q;
    debug_assert!(challenge.on_curve(f, &gc) && challenge.on_curve(f, &qc));

    // ── attack: walk back to a weak curve (charged)
    f.reset_muls();
    let cap = cap_mult * p * p;
    let walk = walk_to_weak_record(&ctx, challenge, jumps, cap, trace, &mut rng);
    rep.walk_muls = f.muls();
    let Some((weak, path)) = walk else {
        rep.stage = String::from("walk_none");
        rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
        return rep;
    };
    rep.walk_found = true;
    rep.path_steps = path.len();
    rep.path_two = path
        .iter()
        .filter(|s| matches!(s, Step::Two { .. }))
        .count();
    rep.path_odd = rep.path_steps - rep.path_two;

    // ── transport G_C, Q_C along the path (charged)
    f.reset_muls();
    let mut gw = gc;
    let mut qw = qc;
    for step in &path {
        gw = map_step(f, step, &gw);
        qw = map_step(f, step, &qw);
    }
    rep.transport_muls = f.muls();
    debug_assert!(weak.on_curve(f, &gw) && weak.on_curve(f, &qw));

    // ── model change onto the cover's form (charged)
    f.reset_muls();
    let Some((alpha, mapped)) = model_change(&ctx, &weak, &[gw, qw]) else {
        rep.model_muls = f.muls();
        rep.stage = String::from("model_change_none");
        rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
        return rep;
    };
    rep.model_muls = f.muls();
    let gm = PtE6 {
        x: mapped[0].0,
        y: mapped[0].1,
        inf: false,
    };
    let qm = PtE6 {
        x: mapped[1].0,
        y: mapped[1].1,
        inf: false,
    };
    let Some(spec) = super::jv_cover::spec_from_curve(p, alpha, l, gm, qm, spec0.d) else {
        rep.stage = String::from("spec_from_curve_none");
        rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
        return rep;
    };

    // ── the sieved route on the transported instance
    let route = run_cover_sieve_dlp_on(&spec, seed, 0, rho_s_ref, stop, margin, None, true);
    rep.route_muls = route.total_muls;
    rep.route_solved = route.solved;
    rep.route_correct = route.correct;
    rep.c_add_e = route.c_add_e;

    // ── verify the recovered scalar on the challenge curve itself
    if route.correct {
        let got = challenge.scalar_mul(f, &gc, spec0.d as u128);
        rep.verified_on_challenge = got == Some(qc);
    }

    // ── cost: every charged phase
    rep.total_muls = rep.walk_muls + rep.transport_muls + rep.model_muls + rep.route_muls;
    rep.reach_share =
        (rep.walk_muls + rep.transport_muls + rep.model_muls) as f64 / rep.total_muls.max(1) as f64;
    let sqrt_l = (l as f64).sqrt();
    rep.s = rep.total_muls as f64 / (rep.c_add_e.max(1.0) * sqrt_l);
    let rho_muls = rho_s_ref * sqrt_l * rep.c_add_e.max(1.0);
    rep.s_over_rho = rep.total_muls as f64 / rho_muls;

    // ── paired walk-free route on the weak curve of the same seed
    let route_only = run_cover_sieve_dlp(p, seed, 0, rho_s_ref, stop, margin, None, true);
    rep.route_only_s_over_rho = route_only.total_muls as f64 / rho_muls;
    rep.e2e_over_route_only = rep.total_muls as f64 / (route_only.total_muls.max(1)) as f64;
    rep.wall_ms = start.elapsed().as_secs_f64() * 1e3;
    rep
}

/// A diagnostic for §17.5 (post hoc, labelled so): are the weak curves
/// spread over the isogeny classes, or concentrated in a few?  `n` random
/// weak curves (the construction `y² = (x − ρ)(x − α)(x − σα)`) and `n`
/// random full-2-torsion curves, each with its trace by [`curve_order`];
/// the traces are the isogeny classes (Tate).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct TraceCensus {
    pub p: u64,
    pub q: u64,
    pub n: u64,
    /// Distinct traces among the weak and among the random curves.
    pub weak_distinct: u64,
    pub random_distinct: u64,
    /// Fraction of the random curves whose trace is one some weak curve
    /// of the census has: an estimate (from below, the census being finite)
    /// of the chance that a random full-2-torsion curve's class holds a
    /// weak curve.
    pub random_in_weak_classes: f64,
    /// The weak traces by frequency: `(t, weak count, random count)`.
    pub top_weak: Vec<(i128, u64, u64)>,
    /// Residues of the weak traces: fraction with `t ≡ r (mod m)` for
    /// `m = 3, 4, 8`, against the random curves' fractions.
    pub weak_mod: Vec<(u64, Vec<f64>)>,
    pub random_mod: Vec<(u64, Vec<f64>)>,
    pub muls: u64,
    pub wall_ms: f64,
}

pub fn trace_census(p: u64, seed: u64, n: u64) -> TraceCensus {
    let start = Instant::now();
    let f = Fq3::new(p);
    let q = p * p;
    let q3 = (p as i128).pow(6);
    let mut rng = StdRng::seed_from_u64(seed ^ 0xCE1505);
    f.reset_muls();
    let mut weak_t: HashMap<i128, u64> = HashMap::new();
    let mut rand_t: HashMap<i128, u64> = HashMap::new();
    for _ in 0..n {
        // a weak curve: ρ ∈ F_q, α ∈ F_{q³} ∖ F_q with ρ, α, σα distinct
        let c = loop {
            let rho = f.from_fq(f.f.random(&mut rng));
            let alpha = f.random(&mut rng);
            if alpha.in_fq() {
                continue;
            }
            let sa = f.sigma(&alpha);
            if alpha != rho && sa != rho {
                break Curve2 {
                    e: [rho, alpha, sa],
                };
            }
        };
        debug_assert!(c.weak_by_norms(&f));
        let t = q3 + 1 - curve_order(&f, &c, &mut rng) as i128;
        *weak_t.entry(t).or_insert(0) += 1;
        let r = random_curve(&f, &mut rng);
        let t = q3 + 1 - curve_order(&f, &r, &mut rng) as i128;
        *rand_t.entry(t).or_insert(0) += 1;
    }
    let in_weak: u64 = rand_t
        .iter()
        .filter(|(t, _)| weak_t.contains_key(t))
        .map(|(_, c)| *c)
        .sum();
    let mut top: Vec<(i128, u64, u64)> = weak_t
        .iter()
        .map(|(t, c)| (*t, *c, *rand_t.get(t).unwrap_or(&0)))
        .collect();
    top.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    top.truncate(12);
    let residues = |m: &HashMap<i128, u64>| -> Vec<(u64, Vec<f64>)> {
        [3u64, 4, 8]
            .iter()
            .map(|&md| {
                let mut v = vec![0u64; md as usize];
                for (t, c) in m {
                    let r = ((t % md as i128) + md as i128) % md as i128;
                    v[r as usize] += c;
                }
                (md, v.iter().map(|&x| x as f64 / n as f64).collect())
            })
            .collect()
    };
    TraceCensus {
        p,
        q,
        n,
        weak_distinct: weak_t.len() as u64,
        random_distinct: rand_t.len() as u64,
        random_in_weak_classes: in_weak as f64 / n as f64,
        top_weak: top,
        weak_mod: residues(&weak_t),
        random_mod: residues(&rand_t),
        muls: f.muls(),
        wall_ms: start.elapsed().as_secs_f64() * 1e3,
    }
}

/// The exact census (post hoc diagnostic of §17.5): **every** weak curve's
/// isogeny class at a small `p`.  Up to the isomorphisms that keep the
/// weak form (`x ↦ x − ρ` with `ρ ∈ F_q`, `x ↦ s²x` with `s ∈ F_q^×`), a
/// weak curve is `y² = x(x − α)(x − σα)` with `α ∈ F_{q³} ∖ F_q` taken up
/// to `F_q^{×2}`: the representatives have the first non-zero of
/// `(α₁, α₂)` in `{1, w}`, `w` the non-residue of `F_q`, `2q² + 2q` of
/// them.  Each gets a probabilistic trace assignment by [`curve_order`];
/// the observed trace set labels the isogeny classes that hold a weak curve.
/// Then `n` random full-2-torsion curves are drawn.  Their fraction in that
/// set estimates, with sampling uncertainty, the chance for a random curve.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ExactCensus {
    pub p: u64,
    pub q: u64,
    pub weak_representatives: u64,
    pub weak_classes: u64,
    /// Among the `n` random curves: distinct traces, and the fraction in a
    /// weak class.
    pub n: u64,
    pub random_distinct: u64,
    pub random_in_weak_classes: f64,
    /// The weak classes by the number of representatives in them, top 12:
    /// `(t, representatives)`.
    pub top_weak: Vec<(i128, u64)>,
    /// §18.2: every weak trace with its representative count, sorted by `t`,
    /// and the traces of the `n` random full-2-torsion curves, kept for
    /// §18.3's characterization.
    #[serde(default)]
    pub weak_traces: Vec<(i128, u64)>,
    #[serde(default)]
    pub random_traces: Vec<i128>,
    /// §18.5: `n_uniform` uniformly random curves `y² = x³ + ax + b` over
    /// `F_{q³}` (all curves, not only full 2-torsion).  A curve with no
    /// rational 2-torsion point has odd order and so no weak curve in its
    /// class; the others are counted by `ell_order`.  `uniform_in_weak_classes`
    /// is the fraction of **all** sampled curves whose trace is a weak class.
    #[serde(default)]
    pub n_uniform: u64,
    #[serde(default)]
    pub uniform_odd_order: u64,
    #[serde(default)]
    pub uniform_four_divides: u64,
    #[serde(default)]
    pub uniform_in_weak_classes: f64,
    #[serde(default)]
    pub uniform_traces: Vec<i128>,
    pub muls: u64,
    pub wall_ms: f64,
}

/// One uniformly random curve `y² = x³ + ax + b` over `F_{q³}` (non-singular):
/// its order through a rational 2-torsion point moved to `0`, or `None` if
/// the cubic has no root in `F_{q³}` (odd order).
fn uniform_curve_order(f: &Fq3, rng: &mut StdRng) -> Option<u128> {
    let c = |k: u64| f.from_fq(E2([k % f.f.p, 0]));
    loop {
        let a = f.random(rng);
        let b = f.random(rng);
        // 4a³ + 27b² ≠ 0
        let disc = f.add(
            &f.mul(&c(4), &f.mul(&a, &f.sq(&a))),
            &f.mul(&c(27), &f.sq(&b)),
        );
        if disc == E6::ZERO {
            continue;
        }
        let cubic: P6 = vec![b, a, E6::ZERO, E6::ONE];
        let roots = roots_q3(f, &cubic, rng);
        let r = roots.first().copied()?;
        // x = X + r:  y² = X³ + 3r X² + (3r² + a) X
        let a2 = f.mul(&c(3), &r);
        let a4 = f.add(&f.mul(&c(3), &f.sq(&r)), &a);
        let ec = EllE::from_a2_a4(f, a2, a4);
        return Some(ell_order(f, &ec, rng));
    }
}

pub fn exact_census(p: u64, seed: u64, n: u64) -> ExactCensus {
    exact_census_full(p, seed, n, 0)
}

pub fn exact_census_full(p: u64, seed: u64, n: u64, n_uniform: u64) -> ExactCensus {
    use rayon::prelude::*;
    let start = Instant::now();
    let f = Fq3::new(p);
    let q = p * p;
    let q3 = (p as i128).pow(6);
    let e2 = |k: u64| E2([k % p, k / p]);
    // f.f.w is a nonsquare in F_p, hence a square in F_{p²}.  The two
    // normal forms require representatives of both F_{p²} square classes.
    let w = (1..q).map(e2).find(|x| !f.f.is_square(x)).unwrap();
    // Enumerate every weak representative in parallel over a0, each task with
    // its own field context (Fq3 holds a non-Sync counter).  `a1 ∈ {1, w}`
    // covers the generic representatives; `a1 = 0` with `a2 ∈ {1, w}` the
    // ones whose first non-zero imaginary part is in a2.
    let both: Vec<(E2, bool)> = [(E2::ONE, false), (w, false)]
        .into_iter()
        .chain([(E2::ONE, true), (w, true)])
        .collect();
    // (a1_or_a2 value, is_a2_branch)
    let partials: Vec<(HashMap<i128, u64>, u64, u64)> = both
        .par_iter()
        .flat_map_iter(|&(val, is_a2)| (0..q).map(move |a0| (val, is_a2, a0)))
        .map(|(val, is_a2, a0)| {
            let f = Fq3::new(p);
            let mut rng = StdRng::seed_from_u64(seed ^ 0xE8AC7 ^ a0.wrapping_mul(0x9E3779B1));
            let mut local: HashMap<i128, u64> = HashMap::new();
            let mut reps = 0u64;
            let a0e = e2(a0);
            if is_a2 {
                // a1 = 0, a2 = val: one curve per a0
                let alpha = E6([a0e, E2::ZERO, val]);
                let c = Curve2 {
                    e: [E6::ZERO, alpha, f.sigma(&alpha)],
                };
                let t = q3 + 1 - curve_order(&f, &c, &mut rng) as i128;
                *local.entry(t).or_insert(0) += 1;
                reps += 1;
            } else {
                for a2 in 0..q {
                    let alpha = E6([a0e, val, e2(a2)]);
                    let c = Curve2 {
                        e: [E6::ZERO, alpha, f.sigma(&alpha)],
                    };
                    let t = q3 + 1 - curve_order(&f, &c, &mut rng) as i128;
                    *local.entry(t).or_insert(0) += 1;
                    reps += 1;
                }
            }
            (local, reps, f.muls())
        })
        .collect();
    let mut weak_t: HashMap<i128, u64> = HashMap::new();
    let mut reps = 0u64;
    let mut census_muls = 0u64;
    for (local, r, m) in partials {
        reps += r;
        census_muls += m;
        for (t, c) in local {
            *weak_t.entry(t).or_insert(0) += c;
        }
    }
    let f = Fq3::new(p);
    let mut rng = StdRng::seed_from_u64(seed ^ 0xE8AC7);
    let mut rand_t: HashMap<i128, u64> = HashMap::new();
    let mut random_traces: Vec<i128> = Vec::with_capacity(n as usize);
    let mut in_weak = 0u64;
    for _ in 0..n {
        let r = random_curve(&f, &mut rng);
        let t = q3 + 1 - curve_order(&f, &r, &mut rng) as i128;
        *rand_t.entry(t).or_insert(0) += 1;
        random_traces.push(t);
        if weak_t.contains_key(&t) {
            in_weak += 1;
        }
    }
    // uniformly random curves over F_{q³}: the reach over all curves
    let mut urng = StdRng::seed_from_u64(seed ^ 0x0A11_C0E5);
    let mut uniform_traces: Vec<i128> = Vec::new();
    let (mut u_odd, mut u_four, mut u_weak) = (0u64, 0u64, 0u64);
    for _ in 0..n_uniform {
        match uniform_curve_order(&f, &mut urng) {
            None => u_odd += 1,
            Some(ord) => {
                if ord % 4 == 0 {
                    u_four += 1;
                }
                let t = q3 + 1 - ord as i128;
                uniform_traces.push(t);
                if weak_t.contains_key(&t) {
                    u_weak += 1;
                }
            }
        }
    }
    let mut top: Vec<(i128, u64)> = weak_t.iter().map(|(t, c)| (*t, *c)).collect();
    top.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut weak_traces = top.clone();
    weak_traces.sort_by_key(|&(t, _)| t);
    top.truncate(12);
    ExactCensus {
        p,
        q,
        weak_representatives: reps,
        weak_classes: weak_t.len() as u64,
        n,
        random_distinct: rand_t.len() as u64,
        random_in_weak_classes: in_weak as f64 / n.max(1) as f64,
        top_weak: top,
        weak_traces,
        random_traces,
        n_uniform,
        uniform_odd_order: u_odd,
        uniform_four_divides: u_four,
        uniform_in_weak_classes: u_weak as f64 / n_uniform.max(1) as f64,
        uniform_traces,
        muls: census_muls + f.muls(),
        wall_ms: start.elapsed().as_secs_f64() * 1e3,
    }
}

/// §18.3 (exploratory, post hoc): what distinguishes the weak classes.
/// Weak membership depends only on the trace `t`, because the Weil
/// restriction's characteristic polynomial over `F_q` is `T⁶ − tT³ + q³`.
/// From frozen `ExactCensus` files: for the full-2-torsion classes the
/// random sample met, the weak fraction by `v₂(D)`, `v₃(D)`, `t mod 3`,
/// `t mod 8` and a class-size proxy (the class's multiplicity among the
/// uniform all-curve sample), with `D = t² − 4q³`; plus twist symmetry
/// (`t` weak ⟹ `−t` weak).  A feature that separates shows fractions of
/// only 0 and 1.  Labelled exploratory; a rule found here is a candidate,
/// to be registered and tested on a size it was not fitted on.
pub fn characterize_weak(reports: &[ExactCensus]) -> String {
    use std::collections::{BTreeMap, HashSet as HSet};
    use std::fmt::Write;
    fn val(mut n: i128, p: i128) -> u32 {
        n = n.abs();
        if n == 0 {
            return 99;
        }
        let mut k = 0;
        while n % p == 0 {
            n /= p;
            k += 1;
        }
        k
    }
    let mut out = String::new();
    for r in reports {
        let q3 = (r.p as i128).pow(6);
        let weak: HSet<i128> = r.weak_traces.iter().map(|&(t, _)| t).collect();
        let classes: HSet<i128> = r.random_traces.iter().copied().collect();
        let mut umult: HashMap<i128, u64> = HashMap::new();
        for &t in &r.uniform_traces {
            *umult.entry(t).or_insert(0) += 1;
        }
        let twist_ok = weak.iter().filter(|t| weak.contains(&-**t)).count();
        let _ = writeln!(
            out,
            "\n### p = {}: {} weak classes (exact), {} full-2-torsion classes sampled, {:.1} % of them weak; weak set closed under t ↦ −t: {}/{}",
            r.p,
            weak.len(),
            classes.len(),
            100.0 * classes.iter().filter(|t| weak.contains(t)).count() as f64
                / classes.len().max(1) as f64,
            twist_ok,
            weak.len()
        );
        let mut feats: Vec<(&str, BTreeMap<String, (u64, u64)>)> = vec![
            ("v2(D)", BTreeMap::new()),
            ("v3(D)", BTreeMap::new()),
            ("t mod 3", BTreeMap::new()),
            ("t mod 8", BTreeMap::new()),
            ("class size proxy (uniform hits)", BTreeMap::new()),
        ];
        for &t in &classes {
            let d = t * t - 4 * q3;
            let w = weak.contains(&t) as u64;
            let um = *umult.get(&t).unwrap_or(&0);
            let vals = [
                format!("{}", val(d, 2)),
                format!("{}", val(d, 3)),
                format!("{}", t.rem_euclid(3)),
                format!("{}", t.rem_euclid(8)),
                format!("{}", um.min(4)),
            ];
            for (i, v) in vals.into_iter().enumerate() {
                let e = feats[i].1.entry(v).or_insert((0, 0));
                e.0 += w;
                e.1 += 1;
            }
        }
        for (name, m) in &feats {
            let cells: Vec<String> = m
                .iter()
                .map(|(v, (w, n))| format!("{v}: {:.2} ({n})", *w as f64 / *n as f64))
                .collect();
            let separates = m.values().all(|(w, n)| *w == 0 || *w == *n);
            let _ = writeln!(
                out,
                "- {name}{}: {}",
                if separates { " — SEPARATES" } else { "" },
                cells.join(" · ")
            );
        }
    }
    out
}

/// The derived tables of ledger §17.5, printed from frozen `Walk2Report`s
/// (so that every number in the note is printed by this code from the
/// experiment file, never computed by hand): per size, the success rate
/// overall and by the number of jump degrees the class admits, the median
/// curves met, the three constants, and the price against rho; then the
/// fits the predictions P1 and P4 ask for.
pub fn summarize_walk2(reports: &[Walk2Report]) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    let _ = writeln!(out, "| p | q | walks (start weak) | success | success by admitted degrees 0 / 1 / 2 / 3 (walks) | exhausted / capped | curves met, median of found (q/3) | first component mean | c_order | c_curve | c_jump | jumps per found walk | walk / rho, found | walk / rho, all |");
    let _ = writeln!(
        out,
        "|---:|--:|:--|--:|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|"
    );
    for r in reports {
        let proper: Vec<&Walk2Trial> = r.rows.iter().filter(|t| !t.start_weak).collect();
        let by_deg: Vec<String> = (0..=3)
            .map(|k| {
                let rows: Vec<&&Walk2Trial> = proper
                    .iter()
                    .filter(|t| t.degrees_used.len() == k)
                    .collect();
                if rows.is_empty() {
                    "–".to_string()
                } else {
                    let f = rows.iter().filter(|t| t.found).count();
                    format!("{:.2} ({})", f as f64 / rows.len() as f64, rows.len())
                }
            })
            .collect();
        let found: Vec<&&Walk2Trial> = proper.iter().filter(|t| t.found).collect();
        let mean_jumps = if found.is_empty() {
            0.0
        } else {
            found
                .iter()
                .map(|t| t.jumps.iter().sum::<u64>() as f64)
                .sum::<f64>()
                / found.len() as f64
        };
        let _ = writeln!(
            out,
            "| {} | {} | {} ({}) | {:.2} | {} | {} / {} | {:.0} ({:.0}) | {:.1} | {:.2e} | {:.0} | {:.2e} | {:.1} | {:.4} | {:.4} |",
            r.p, r.q, r.trials, r.start_weak, r.success_fraction, by_deg.join(" / "), r.exhausted, r.capped,
            r.median_curves, r.q_over_3, r.mean_first_component, r.c_order, r.c_curve, r.c_jump, mean_jumps,
            r.walk_over_rho, r.walk_over_rho_all
        );
    }
    // fits: c_curve against log p; walk/rho (found) against p, sizes ≥ 53
    let fit = |xs: &[f64], ys: &[f64]| -> (f64, f64) {
        let n = xs.len() as f64;
        let mx = xs.iter().sum::<f64>() / n;
        let my = ys.iter().sum::<f64>() / n;
        let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
        let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
        let a = sxy / sxx;
        let b = my - a * mx;
        (a, b)
    };
    let big: Vec<&Walk2Report> = reports
        .iter()
        .filter(|r| r.p >= 53 && r.found > 0)
        .collect();
    if big.len() >= 3 {
        let xs: Vec<f64> = big.iter().map(|r| (r.p as f64).ln()).collect();
        let ys: Vec<f64> = big.iter().map(|r| r.walk_over_rho.ln()).collect();
        let (a, _) = fit(&xs, &ys);
        let _ = writeln!(
            out,
            "\nwalk / rho (found walks) ∝ p^{a:.2} over p ≥ 53 ({} sizes).",
            big.len()
        );
        let ys2: Vec<f64> = big.iter().map(|r| r.c_curve.ln()).collect();
        let (a2, _) = fit(&xs, &ys2);
        let ys3: Vec<f64> = big.iter().map(|r| r.c_jump.ln()).collect();
        let (a3, _) = fit(&xs, &ys3);
        let _ = writeln!(
            out,
            "c_curve ∝ p^{a2:.2}, c_jump ∝ p^{a3:.2} over the same sizes."
        );
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::jv_cover::generate_spec;

    #[test]
    fn fp_nonsquare_is_not_an_fp2_square_class_representative() {
        for p in [7_u64, 11, 13, 17, 37] {
            let f = Fq3::new(p);
            let old = E2([f.f.w, 0]);
            assert!(
                f.f.is_square(&old),
                "p={p}: F_p element must square in F_p2"
            );
            let q = p * p;
            let e2 = |k: u64| E2([k % p, k / p]);
            let w = (1..q).map(e2).find(|x| !f.f.is_square(x)).unwrap();
            assert!(!f.f.is_square(&w), "p={p}: missing second square class");
        }
    }

    #[test]
    fn the_constructed_weak_curves_pass_the_test_and_random_ones_mostly_fail() {
        for p in [7u64, 11, 13] {
            let f = Fq3::new(p);
            let mut rng = StdRng::seed_from_u64(1);
            // the instances of the cover route: y² = x(x − α)(x − σα)
            for seed in 1..4 {
                let spec = generate_spec(p, seed);
                let c = Curve2 {
                    e: [E6::ZERO, spec.alpha, f.sigma(&spec.alpha)],
                };
                assert_eq!(c.weak_root(&f), Some(0), "p = {p}, seed {seed}");
                // and the class is invariant under x ↦ u²x + r
                let (u, r) = (f.random(&mut rng), f.random(&mut rng));
                let u2 = f.sq(&u);
                let moved = Curve2 {
                    e: [
                        f.add(&f.mul(&u2, &c.e[0]), &r),
                        f.add(&f.mul(&u2, &c.e[1]), &r),
                        f.add(&f.mul(&u2, &c.e[2]), &r),
                    ],
                };
                assert!(moved.weak_root(&f).is_some());
                assert_eq!(moved.j(&f), c.j(&f));
            }
            let n = 600;
            let weak = (0..n)
                .filter(|_| random_curve(&f, &mut rng).weak_root(&f).is_some())
                .count();
            let expect = 3.0 * n as f64 / (p * p) as f64;
            assert!(
                (weak as f64) < 3.0 * expect + 6.0 && (weak as f64) > expect / 3.0 - 3.0,
                "p = {p}: {weak} weak of {n}, expected ≈ {expect:.1}"
            );
        }
    }

    #[test]
    fn isogenous_curves_have_the_same_order_class_and_rational_torsion() {
        // the 2- and 3-isogeny targets are curves with full 2-torsion; a
        // 2-isogeny's target composed back reaches the start's j (the dual)
        let p = 11u64;
        let f = Fq3::new(p);
        let mut rng = StdRng::seed_from_u64(2);
        let mut twos = 0;
        let mut threes = 0;
        for _ in 0..20 {
            let c = random_curve(&f, &mut rng);
            for n in c.two_isogenies(&f, &mut rng) {
                twos += 1;
                assert!(n.e[1] != n.e[2] && n.e[1] != E6::ZERO);
                // the dual: one of n's 2-isogenies has c's j
                let back: Vec<E6> = n
                    .two_isogenies(&f, &mut rng)
                    .iter()
                    .map(|m| m.j(&f))
                    .collect();
                assert!(back.contains(&c.j(&f)), "no dual back to the start");
            }
            for n in c.three_isogenies(&f, &mut rng) {
                threes += 1;
                let back: Vec<E6> = n
                    .three_isogenies(&f, &mut rng)
                    .iter()
                    .map(|m| m.j(&f))
                    .collect();
                assert!(back.contains(&c.j(&f)), "no dual 3-isogeny back");
            }
        }
        assert!(twos >= 10, "{twos}");
        assert!(threes >= 3, "{threes}");
    }

    #[test]
    fn a_short_walk_reaches_a_weak_curve_at_p7() {
        let r = run_walk(7, 1, 20, 6000, true, 2000);
        // at p = 7 the reachable components are small: some walks exhaust
        // theirs without a weak curve; the ones that find one do so in about
        // q/3 steps
        assert!(r.found >= 8, "found {} of {}", r.found, r.trials);
        assert!(r.capped == 0, "walks hit the cap: {}", r.capped);
        assert!(r.mean_steps < 5.0 * r.q_over_3, "{}", r.mean_steps);
        assert!(
            r.weak_fraction_times_q > 1.0 && r.weak_fraction_times_q < 6.0,
            "{}",
            r.weak_fraction_times_q
        );
    }

    #[test]
    fn the_norm_test_agrees_with_the_cross_ratio_test() {
        for p in [7u64, 11, 13, 17] {
            let f = Fq3::new(p);
            let mut rng = StdRng::seed_from_u64(3);
            for seed in 1..4 {
                let spec = generate_spec(p, seed);
                let c = Curve2 {
                    e: [E6::ZERO, spec.alpha, f.sigma(&spec.alpha)],
                };
                assert!(c.weak_by_norms(&f));
                // every ordering of the triple
                let perms = [
                    [0, 1, 2],
                    [1, 2, 0],
                    [2, 0, 1],
                    [0, 2, 1],
                    [1, 0, 2],
                    [2, 1, 0],
                ];
                for pm in perms {
                    let d = Curve2 {
                        e: [c.e[pm[0]], c.e[pm[1]], c.e[pm[2]]],
                    };
                    assert!(d.weak_by_norms(&f), "ordering {pm:?}");
                    assert!(d.weak_root(&f).is_some(), "ordering {pm:?}");
                }
            }
            let mut agree = 0;
            let mut weak = 0;
            for _ in 0..400 {
                let c = random_curve(&f, &mut rng);
                let a = c.weak_root(&f).is_some();
                let b = c.weak_by_norms(&f);
                assert_eq!(a, b, "p = {p}: {c:?}");
                agree += 1;
                weak += b as usize;
            }
            assert_eq!(agree, 400);
            assert!(weak > 0 || p > 11, "p = {p}: no weak curve in 400 samples");
        }
    }

    #[test]
    fn the_fast_character_and_root_agree_with_the_slow_ones() {
        for p in [7u64, 11, 13, 17, 23] {
            let ctx = WalkCtx::new(p);
            let f = &ctx.f;
            let mut rng = StdRng::seed_from_u64(4);
            let mut squares = 0;
            for _ in 0..200 {
                let w = f.random(&mut rng);
                assert_eq!(ctx.is_square_q3(&w), f.is_square(&w), "p = {p}: {w:?}");
                if let Some(r) = ctx.sqrt_q3(&w) {
                    squares += 1;
                    assert_eq!(f.sq(&r), w);
                    assert!(f.sqrt(&w, &mut rng).is_some());
                } else {
                    assert!(f.sqrt(&w, &mut rng).is_none());
                }
            }
            assert!(
                squares > 60 && squares < 140,
                "p = {p}: {squares} squares of 200"
            );
            // the fast root is cheaper than Tonelli–Shanks in F_{q³}
            let w = loop {
                let w = f.random(&mut rng);
                if f.is_square(&w) && w != E6::ZERO {
                    break w;
                }
            };
            f.reset_muls();
            ctx.sqrt_q3(&w).unwrap();
            let fast = f.muls();
            f.reset_muls();
            f.sqrt(&w, &mut rng).unwrap();
            let slow = f.muls();
            assert!(fast * 2 < slow, "p = {p}: fast {fast}, slow {slow}");
        }
    }

    #[test]
    fn two_edges_match_the_old_neighbours_and_carry_the_dual_root() {
        let p = 13u64;
        let ctx = WalkCtx::new(p);
        let f = &ctx.f;
        let mut rng = StdRng::seed_from_u64(5);
        let mut edges = 0;
        for _ in 0..30 {
            let c = random_curve(f, &mut rng);
            let old: Vec<E6> = c
                .two_isogenies(f, &mut rng)
                .iter()
                .map(|m| m.j(f))
                .collect();
            let new = c.two_edges(&ctx, None);
            let mut nj: Vec<E6> = new.iter().map(|e| e.curve.j(f)).collect();
            let mut oj = old.clone();
            nj.sort();
            oj.sort();
            assert_eq!(nj, oj, "{c:?}");
            for e in new {
                edges += 1;
                // the carried root is the root of u′v′ at kernel 0 of the image
                let u = f.sub(&e.curve.e[1], &e.curve.e[0]);
                let v = f.sub(&e.curve.e[2], &e.curve.e[0]);
                assert_eq!(f.sq(&e.dual_root), f.mul(&u, &v));
                // and following it returns to the start's class
                let back = e.curve.two_edges(&ctx, Some((0, e.dual_root)));
                assert!(back.iter().any(|b| b.curve.j(f) == c.j(f)), "no dual");
                // with the known root, that edge costs no square root
                f.reset_muls();
                let _ = e.curve.two_edges(&ctx, Some((0, e.dual_root)));
                let with = f.muls();
                f.reset_muls();
                let _ = e.curve.two_edges(&ctx, None);
                let without = f.muls();
                assert!(with < without, "{with} vs {without}");
            }
        }
        assert!(edges >= 20, "{edges}");
    }

    #[test]
    fn ell_targets_are_isogenous_and_three_matches_the_old_code() {
        let p = 11u64;
        let ctx = WalkCtx::new(p);
        let f = &ctx.f;
        let mut rng = StdRng::seed_from_u64(6);
        let mut counts = [0usize; 3];
        for _ in 0..40 {
            let c = random_curve(f, &mut rng);
            // ℓ = 3: the same j-set as the old three_isogenies
            let mut oj: Vec<E6> = c
                .three_isogenies(f, &mut rng)
                .iter()
                .map(|m| m.j(f))
                .collect();
            let mut nj: Vec<E6> = c
                .ell_targets(&ctx, 3, &mut rng)
                .iter()
                .map(|m| m.j(f))
                .collect();
            oj.sort();
            oj.dedup();
            nj.sort();
            nj.dedup();
            assert_eq!(nj, oj, "{c:?}");
            for (i, ell) in [3u64, 5, 7].iter().enumerate() {
                for t in c.ell_targets(&ctx, *ell, &mut rng) {
                    counts[i] += 1;
                    // the dual ℓ-isogeny of the image leads back to the start's class
                    let back: Vec<E6> = t
                        .ell_targets(&ctx, *ell, &mut rng)
                        .iter()
                        .map(|m| m.j(f))
                        .collect();
                    assert!(
                        back.contains(&c.j(f)),
                        "ℓ = {ell}: no dual back, {c:?} -> {t:?}"
                    );
                }
            }
        }
        assert!(
            counts[0] >= 5 && counts[1] >= 3 && counts[2] >= 2,
            "{counts:?}"
        );
    }

    #[test]
    fn the_rebuilt_walk_reaches_weak_curves_at_p7_and_p13() {
        for (p, trials) in [(7u64, 20usize), (13, 12)] {
            let r = run_walk2(p, 1, trials, 3 * p * p, &[3, 5, 7], 1000, false, false, 8);
            // every walk ends found, capped or exhausted; at these sizes many
            // closures hold no weak curve, so exhaustion is common
            assert_eq!(
                r.found + r.capped + r.exhausted + r.start_weak
                    - r.rows.iter().filter(|t| t.start_weak && t.found).count(),
                trials,
                "{r:?}"
            );
            assert!(r.found >= 2, "p = {p}: found {} of {}", r.found, r.trials);
            assert!(
                r.weak_fraction_times_q > 1.0 && r.weak_fraction_times_q < 6.0,
                "{}",
                r.weak_fraction_times_q
            );
            assert!(r.c_curve < 20_000.0, "p = {p}: c_curve {}", r.c_curve);
            for t in &r.rows {
                assert!(
                    t.order % 4 == 0 && t.trace.unsigned_abs() <= 2 * (p as u128).pow(3),
                    "{t:?}"
                );
            }
        }
    }

    #[test]
    fn the_point_count_is_the_order_and_matches_the_cover_instances() {
        for p in [7u64, 11, 13] {
            let f = Fq3::new(p);
            let mut rng = StdRng::seed_from_u64(9);
            for seed in 1..3 {
                let spec = generate_spec(p, seed);
                let c = Curve2 {
                    e: [E6::ZERO, spec.alpha, f.sigma(&spec.alpha)],
                };
                let m = curve_order(&f, &c, &mut rng);
                assert_eq!(m % 4, 0);
                assert_eq!(m, 4 * spec.l as u128, "p = {p} seed {seed}");
                // every 2-isogenous curve has the same order
                let ctx = WalkCtx::new(p);
                for e in c.two_edges(&ctx, None) {
                    assert_eq!(curve_order(&f, &e.curve, &mut rng), m);
                }
                for t in c.ell_targets(&ctx, 3, &mut rng) {
                    assert_eq!(curve_order(&f, &t, &mut rng), m);
                }
            }
        }
    }

    #[test]
    fn point_maps_preserve_the_curve_and_the_logarithm() {
        // On a weak curve, a 2-isogeny and an ℓ-isogeny each carry an
        // order-ℓ point to its image curve and commute with scalar mult.
        for p in [7u64, 11, 13] {
            let ctx = WalkCtx::new(p);
            let f = &ctx.f;
            let mut rng = StdRng::seed_from_u64(7);
            let spec = generate_spec(p, 1);
            let cur = Curve2 {
                e: [E6::ZERO, spec.alpha, f.sigma(&spec.alpha)],
            };
            let g: Pt2 = (spec.g.x, spec.g.y);
            assert!(cur.on_curve(f, &g), "p = {p}");
            // a 2-edge
            let e = cur.two_edges(&ctx, None)[0];
            let gi = map_two(f, &cur, e.k, &g);
            assert!(e.curve.on_curve(f, &gi), "p = {p}: 2-image off curve");
            // [2]·image == image of [2]·g, so the map is a homomorphism
            let two_g = cur.scalar_mul(f, &g, 2).unwrap();
            assert_eq!(
                e.curve.scalar_mul(f, &gi, 2),
                Some(map_two(f, &cur, e.k, &two_g)),
                "p = {p}: 2-isogeny not a homomorphism"
            );
            // an ℓ-isogeny target with its kernel
            if let Some((tgt, ker)) = cur.ell_targets_with_kernels(&ctx, 3, &mut rng).pop() {
                let gi = map_odd(f, &cur, &ker, &g);
                assert!(tgt.on_curve(f, &gi), "p = {p}: 3-image off curve");
                let two_g = cur.scalar_mul(f, &g, 2).unwrap();
                assert_eq!(
                    tgt.scalar_mul(f, &gi, 2),
                    Some(map_odd(f, &cur, &ker, &two_g)),
                    "p = {p}: 3-isogeny not a homomorphism"
                );
            }
        }
    }

    #[test]
    fn end_to_end_from_a_non_weak_curve_recovers_and_verifies() {
        // A whole-method run (§18.1) at a tiny size: the challenge curve is
        // not weak, the walk reaches a weak curve, the route solves, and the
        // recovered scalar verifies on the challenge curve.
        let r = run_end_to_end(101, 1, 8, &[3, 5, 7], 3, 1.3, 1.25, Some(64));
        assert!(r.challenge_non_weak, "challenge was weak: {r:?}");
        assert!(r.walk_found, "walk found no weak curve: {r:?}");
        assert!(r.route_solved && r.route_correct, "route failed: {r:?}");
        assert!(r.verified_on_challenge, "not verified on C: {r:?}");
        assert!(r.total_muls > 0 && r.reach_share >= 0.0 && r.reach_share <= 1.0);
    }
}
