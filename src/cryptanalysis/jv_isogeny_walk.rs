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

use std::collections::HashSet;
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use super::jv_cover::{Fld, Fq3, E2, E6};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::jv_cover::generate_spec;

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
}
