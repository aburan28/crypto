//! # What it costs to decide a target by a Semaev S₃ solve.
//!
//! The yield of a small-x base ([`super::prime_small_x`],
//! [`super::binary_small_x`]) counts how often a random subgroup element `T`
//! is a signed sum of two base points. It says nothing of what one *trial*
//! costs. Deciding a trial without the `B²`-entry pair-sum table means
//! solving Semaev's summation polynomial `S₃(x₁, x₂, X) = 0` with `X = x(T)`
//! for every base abscissa `x₁` and asking whether a root is a base abscissa.
//! The polynomial has degree two in `x₂`, so each `x₁` costs one quadratic.
//!
//! * Prime fields `y² = x³ + a·x + b`:
//!   `S₃ = (x₁ − x₂)²X² − 2((x₁ + x₂)(x₁x₂ + a) + 2b)X + (x₁x₂ − a)² − 4b(x₁ + x₂)`.
//!   In `x₂` with `u = x₁`, `w = X`: `A·x₂² − 2t·x₂ + C`, `A = (w − u)²`,
//!   `t = uw(u + w) + a(u + w) + 2b`, `C = u²w² − 2w(ua + 2b) + a² − 4bu`, and
//!   the roots are `(t ± √(t² − A·C))/A`. A square root per `x₁`.
//! * Binary fields `y² + xy = x³ + a·x² + b`:
//!   `S₃ = x₁²x₂² + x₁²X² + x₂²X² + x₁x₂X + b`, independent of `a`
//!   (the resultant of the group law over `F₂`, computed with a computer
//!   algebra system and checked against the group law in the tests below).
//!   In `x₂`: `A·x₂² + B·x₂ + C`, `A = (u + w)²`, `B = uw`, `C = u²w² + b`;
//!   with `x₂ = (B/A)·t` it is `t² + t = AC/B²`, solvable iff the trace is
//!   zero, and the root is the half-trace. Trace and half-trace are
//!   `F₂`-linear, so each is a precomputed table applied to the bits.
//!
//! The decision agrees with the pair-sum table on every element of small
//! fields (tests), so the two are interchangeable and the S₃ time is the
//! cost of a trial *without* the table. [`S3Cost`] reports it, the time to
//! collect `B + 1` relations at the exact yield, and Pollard rho on the same
//! machine, both in seconds.

use std::collections::HashSet;
use std::time::Instant;

use num_bigint::BigUint;
use rand::{rngs::StdRng, Rng, SeedableRng};

use super::binary_small_x::{self as bin, BinCurve, Elem, Pt as BPt};
use super::prime_small_x::{
    self as pr, check_memory, Curve, Options, PrimeField, Pt as PPt, S3Cost, YieldReport,
};

// ── Binary fields ──────────────────────────────────────────────────

/// The S₃ decision for a binary curve and a fixed base.
pub struct BinaryS3<'a, const L: usize> {
    c: &'a BinCurve<L>,
    xs: Vec<Elem<L>>,
    inv_xs: Vec<Elem<L>>,
    set: HashSet<Elem<L>>,
    trace_mask: Elem<L>,
    half_trace: Vec<Elem<L>>,
}

impl<'a, const L: usize> BinaryS3<'a, L> {
    pub fn new(c: &'a BinCurve<L>, base: &[BPt<L>]) -> Self {
        let f = &c.f;
        let xs: Vec<Elem<L>> = base.iter().map(|p| p.x).collect();
        let inv_xs = xs.iter().map(|x| f.inv(x)).collect();
        let set = xs.iter().copied().collect();
        let m = f.m as usize;
        let mut trace_mask = f.zero();
        let mut half_trace = Vec::with_capacity(m);
        for i in 0..m {
            let mut e = f.zero();
            e[i / 64] |= 1u64 << (i % 64);
            if f.trace(&e) == 1 {
                trace_mask[i / 64] |= 1u64 << (i % 64);
            }
            // HT(e) = Σ_{j ≤ (m−1)/2} e^(4^j), the root of t² + t = e + Tr(e).
            let (mut t, mut acc) = (e, e);
            for _ in 0..(m - 1) / 2 {
                t = f.sqr(&f.sqr(&t));
                acc = f.add(&acc, &t);
            }
            half_trace.push(acc);
        }
        Self {
            c,
            xs,
            inv_xs,
            set,
            trace_mask,
            half_trace,
        }
    }

    fn trace_is_zero(&self, z: &Elem<L>) -> bool {
        let ones: u32 = (0..L)
            .map(|i| (z[i] & self.trace_mask[i]).count_ones())
            .sum();
        ones % 2 == 0
    }

    fn half_trace_of(&self, z: &Elem<L>) -> Elem<L> {
        let mut acc = self.c.f.zero();
        for (i, h) in self.half_trace.iter().enumerate() {
            if (z[i / 64] >> (i % 64)) & 1 == 1 {
                for k in 0..L {
                    acc[k] ^= h[k];
                }
            }
        }
        acc
    }

    /// Whether `X = w ≠ 0` is `x(±a ± b)` for base points `a, b`, by solving
    /// `S₃(u, ·, w)` for every base abscissa `u`.
    pub fn decide(&self, w: &Elem<L>) -> bool {
        let f = &self.c.f;
        let iw = f.inv(w);
        let n = self.xs.len();
        // 1/(u + w) for all u by one inversion (zeros, u = w, are skipped).
        let dens: Vec<Elem<L>> = self.xs.iter().map(|u| f.add(u, w)).collect();
        let mut prefix: Vec<Elem<L>> = Vec::with_capacity(n);
        let mut acc = f.one();
        for d in &dens {
            if !f.is_zero(d) {
                acc = f.mul(&acc, d);
            }
            prefix.push(acc);
        }
        let mut inv_acc = f.inv(&acc);
        let mut invs = vec![f.zero(); n];
        for k in (0..n).rev() {
            if f.is_zero(&dens[k]) {
                continue;
            }
            let before = if k == 0 { f.one() } else { prefix[k - 1] };
            invs[k] = f.mul(&inv_acc, &before);
            inv_acc = f.mul(&inv_acc, &dens[k]);
        }
        for k in 0..n {
            let u = &self.xs[k];
            let uw = f.mul(u, w);
            let c0 = f.add(&f.sqr(&uw), &self.c.b);
            let d = &dens[k];
            if f.is_zero(d) {
                // u = w: A = 0, the equation is linear: x₂ = C/B.
                let v = f.mul(&c0, &f.mul(&self.inv_xs[k], &iw));
                if self.set.contains(&v) {
                    return true;
                }
                continue;
            }
            // z = A·C/B² = (d/(uw))²·C.
            let r = f.mul(&f.mul(d, &self.inv_xs[k]), &iw);
            let z = f.mul(&f.sqr(&r), &c0);
            if !self.trace_is_zero(&z) {
                continue;
            }
            let t = self.half_trace_of(&z);
            // x₂ = (B/A)·t with B/A = uw/(u + w)², and the other root adds B/A.
            let ratio = f.mul(&uw, &f.sqr(&invs[k]));
            let v1 = f.mul(&ratio, &t);
            if self.set.contains(&v1) || self.set.contains(&f.add(&v1, &ratio)) {
                return true;
            }
        }
        false
    }
}

/// Measure the S₃ cost on a binary curve; see [`S3Cost`].
pub fn s3_cost_binary<const L: usize>(
    c: &BinCurve<L>,
    opts: &Options,
    rep: &YieldReport,
    trials: u64,
) -> Result<S3Cost, String> {
    check_memory(opts.size, 8 * L)?;
    let (base, _) = bin::small_x_base(c, opts.size)?;
    let keys = bin::pair_sum_keys(c, &base);
    let ctx = BinaryS3::new(c, &base);
    let mut rng = StdRng::seed_from_u64(opts.seed ^ 0x5333_434f_5354);
    let targets: Vec<BPt<L>> = (0..trials)
        .map(|_| c.random_subgroup_point(&mut rng))
        .collect();

    let t0 = Instant::now();
    let s3_hits: Vec<bool> = targets.iter().map(|p| ctx.decide(&p.x)).collect();
    let s3_secs = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let table_hits: Vec<bool> = targets
        .iter()
        .map(|p| keys.binary_search(&p.x).is_ok())
        .collect();
    let table_secs = t1.elapsed().as_secs_f64();
    let disagreements = s3_hits
        .iter()
        .zip(&table_hits)
        .filter(|(a, b)| a != b)
        .count() as u64;

    // Planted decompositions must be found.
    let (mut planted, mut found) = (0usize, 0usize);
    while planted < 16 {
        let a = base[rng.gen_range(0..base.len())];
        let mut b = base[rng.gen_range(0..base.len())];
        if rng.gen::<bool>() {
            b = c.negate(&b);
        }
        let Some(t) = c.add_points(&a, &b) else {
            continue;
        };
        planted += 1;
        found += usize::from(ctx.decide(&t.x));
    }

    let f = &c.f;
    let mut x = base[0].x;
    let y = base[1].x;
    let t2 = Instant::now();
    for _ in 0..200_000 {
        x = f.mul(&x, &y);
    }
    let mul_ns = t2.elapsed().as_secs_f64() * 1e9 / 200_000.0;
    std::hint::black_box(x);
    let t3 = Instant::now();
    let adds = 400usize;
    for i in 0..adds {
        std::hint::black_box(c.add_points(&base[i % base.len()], &base[(i + 1) % base.len()]));
    }
    let add_us = t3.elapsed().as_secs_f64() * 1e6 / adds as f64;
    Ok(finish(
        opts.size,
        trials,
        s3_secs,
        table_secs,
        mul_ns,
        add_us,
        disagreements,
        found,
        planted,
        rep,
    ))
}

// ── Prime fields ───────────────────────────────────────────────────

/// The S₃ decision for a prime-field curve and a fixed base.
pub struct PrimeS3<'a, F: PrimeField> {
    c: &'a Curve<F>,
    xs: Vec<F::E>,
    set: HashSet<F::K>,
}

impl<'a, F: PrimeField> PrimeS3<'a, F> {
    pub fn new(c: &'a Curve<F>, base: &[PPt<F::E>]) -> Self {
        let f = &c.field;
        Self {
            c,
            xs: base.iter().map(|p| p.x).collect(),
            set: base.iter().map(|p| f.key(p.x)).collect(),
        }
    }

    /// Whether the abscissa `w` is `x(±a ± b)` for base points `a, b`.
    pub fn decide(&self, f: &mut F, w: F::E) -> bool {
        let n = self.xs.len();
        let (a, b) = (self.c.a, self.c.b);
        let two = f.embed_u64(2);
        let four = f.embed_u64(4);
        let two_b = f.mul(two, b);
        let four_b = f.mul(four, b);
        let a2 = f.sqr(a);
        let w2 = f.sqr(w);
        // Per u: A = (w − u)², t, and t² − A·C.
        let mut cap_a = Vec::with_capacity(n);
        let mut ts = Vec::with_capacity(n);
        let mut disc = Vec::with_capacity(n);
        for &u in &self.xs {
            let d = f.sub(w, u);
            let big_a = f.sqr(d);
            let uw = f.mul(u, w);
            let sum = f.add(u, w);
            let t = {
                let p = f.mul(uw, sum);
                let q = f.mul(a, sum);
                let r = f.add(p, q);
                f.add(r, two_b)
            };
            let c0 = {
                let ua = f.mul(u, a);
                let inner = f.add(ua, two_b);
                let mid = f.mul(w, inner);
                let two_mid = f.mul(two, mid);
                let u2w2 = f.sqr(uw);
                let first = f.sub(u2w2, two_mid);
                let first = f.add(first, a2);
                let fbu = f.mul(four_b, u);
                f.sub(first, fbu)
            };
            let ac = f.mul(big_a, c0);
            let t2 = f.sqr(t);
            cap_a.push(big_a);
            ts.push(t);
            disc.push(f.sub(t2, ac));
        }
        let _ = w2;
        // 1/A for all u with A ≠ 0 by one inversion.
        let mut prefix = Vec::with_capacity(n);
        let mut acc = f.one();
        for &big_a in &cap_a {
            if !f.is_zero(big_a) {
                acc = f.mul(acc, big_a);
            }
            prefix.push(acc);
        }
        let mut inv_acc = f.inv(acc).expect("a product of nonzero elements");
        let mut invs = vec![f.zero(); n];
        for k in (0..n).rev() {
            if f.is_zero(cap_a[k]) {
                continue;
            }
            let before = if k == 0 { f.one() } else { prefix[k - 1] };
            invs[k] = f.mul(inv_acc, before);
            inv_acc = f.mul(inv_acc, cap_a[k]);
        }
        for k in 0..n {
            if f.is_zero(cap_a[k]) {
                // u = w: A = 0, the equation is linear, −2t·x₂ + C = 0, and
                // t² − A·C = t²: x₂ = C/(2t) with C = (t² − disc)/A undefined,
                // so recompute C directly.
                let u = self.xs[k];
                let ua = f.mul(u, a);
                let inner = f.add(ua, two_b);
                let mid = f.mul(w, inner);
                let two_mid = f.mul(two, mid);
                let uw = f.mul(u, w);
                let u2w2 = f.sqr(uw);
                let first = f.sub(u2w2, two_mid);
                let first = f.add(first, a2);
                let fbu = f.mul(four_b, u);
                let cc = f.sub(first, fbu);
                let two_t = f.mul(two, ts[k]);
                if let Some(inv) = f.inv(two_t) {
                    let v = f.mul(cc, inv);
                    if self.set.contains(&f.key(v)) {
                        return true;
                    }
                }
                continue;
            }
            let Some(s) = self.c.sqrt(f, disc[k]) else {
                continue;
            };
            let plus = f.add(ts[k], s);
            let minus = f.sub(ts[k], s);
            let v1 = f.mul(plus, invs[k]);
            let v2 = f.mul(minus, invs[k]);
            if self.set.contains(&f.key(v1)) || self.set.contains(&f.key(v2)) {
                return true;
            }
        }
        false
    }
}

/// Measure the S₃ cost on a prime-field curve; see [`S3Cost`].
pub fn s3_cost_prime<F: PrimeField>(
    c: &Curve<F>,
    opts: &Options,
    rep: &YieldReport,
    trials: u64,
) -> Result<S3Cost, String> {
    check_memory(opts.size, std::mem::size_of::<F::K>())?;
    let base = pr::small_x_base(c, opts.size)?;
    let keys = pr::pair_sum_keys(c, &base);
    let ctx = PrimeS3::new(c, &base);
    let mut f = c.field.clone();
    let exp = (f.modulus() - 1u8) >> 1;
    let mut rng = StdRng::seed_from_u64(opts.seed ^ 0x5333_434f_5354);
    let targets: Vec<PPt<F::E>> = (0..trials)
        .map(|_| c.random_subgroup_point(&mut f, &mut rng, &exp))
        .collect();

    let t0 = Instant::now();
    let s3_hits: Vec<bool> = targets.iter().map(|p| ctx.decide(&mut f, p.x)).collect();
    let s3_secs = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let table_hits: Vec<bool> = targets
        .iter()
        .map(|p| keys.binary_search(&f.key(p.x)).is_ok())
        .collect();
    let table_secs = t1.elapsed().as_secs_f64();
    let disagreements = s3_hits
        .iter()
        .zip(&table_hits)
        .filter(|(a, b)| a != b)
        .count() as u64;

    let (mut planted, mut found) = (0usize, 0usize);
    while planted < 16 {
        let a = base[rng.gen_range(0..base.len())];
        let mut b = base[rng.gen_range(0..base.len())];
        if rng.gen::<bool>() {
            b = c.negate(&mut f, &b);
        }
        let Some(t) = c.add_points(&mut f, &a, &b) else {
            continue;
        };
        planted += 1;
        found += usize::from(ctx.decide(&mut f, t.x));
    }

    let mut x = base[0].x;
    let y = base[1].x;
    let t2 = Instant::now();
    for _ in 0..200_000 {
        x = f.mul(x, y);
    }
    let mul_ns = t2.elapsed().as_secs_f64() * 1e9 / 200_000.0;
    std::hint::black_box(&x);
    let t3 = Instant::now();
    let adds = 400usize;
    for i in 0..adds {
        std::hint::black_box(c.add_points(
            &mut f,
            &base[i % base.len()],
            &base[(i + 1) % base.len()],
        ));
    }
    let add_us = t3.elapsed().as_secs_f64() * 1e6 / adds as f64;
    Ok(finish(
        opts.size,
        trials,
        s3_secs,
        table_secs,
        mul_ns,
        add_us,
        disagreements,
        found,
        planted,
        rep,
    ))
}

#[allow(clippy::too_many_arguments)]
fn finish(
    size: usize,
    trials: u64,
    s3_secs: f64,
    table_secs: f64,
    mul_ns: f64,
    add_us: f64,
    disagreements: u64,
    found: usize,
    planted: usize,
    rep: &YieldReport,
) -> S3Cost {
    let n = trials.max(1) as f64;
    let s3_us = s3_secs * 1e6 / n;
    let s3_full = rep.trials_for_full_system_log2 + (s3_us * 1e-6).log2();
    let rho = rep.rho_log2 + (add_us * 1e-6).log2();
    S3Cost {
        base_abscissae: size,
        trials,
        s3_us_per_trial: s3_us,
        table_ns_per_trial: table_secs * 1e9 / n,
        field_mul_ns: mul_ns,
        s3_field_muls_per_trial: s3_us * 1000.0 / mul_ns,
        group_add_us: add_us,
        disagreements,
        planted_found: found,
        planted,
        s3_full_system_log2_s: s3_full,
        rho_log2_s: rho,
        bits_slower_than_rho: s3_full - rho,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::secg_standards::secg;

    /// `S₃` evaluated at three prime-field abscissae.
    fn prime_s3<F: PrimeField>(c: &Curve<F>, f: &mut F, x1: F::E, x2: F::E, x3: F::E) -> F::E {
        // (x₁−x₂)²x₃² − 2((x₁+x₂)(x₁x₂+a)+2b)x₃ + (x₁x₂−a)² − 4b(x₁+x₂)
        let d = f.sub(x1, x2);
        let d2 = f.sqr(d);
        let x3s = f.sqr(x3);
        let t1 = f.mul(d2, x3s);
        let sum = f.add(x1, x2);
        let p = f.mul(x1, x2);
        let inner = f.add(p, c.a);
        let big = f.mul(sum, inner);
        let two = f.embed_u64(2);
        let four = f.embed_u64(4);
        let two_b = f.mul(two, c.b);
        let big = f.add(big, two_b);
        let tb = f.mul(two, big);
        let t2 = f.mul(tb, x3);
        let pa = f.sub(p, c.a);
        let t3 = f.sqr(pa);
        let fb = f.mul(four, c.b);
        let t4 = f.mul(fb, sum);
        let v = f.sub(t1, t2);
        let v = f.add(v, t3);
        f.sub(v, t4)
    }

    #[test]
    fn prime_s3_polynomial_vanishes_on_group_law_sums() {
        let c = pr::toy_curve(16, 5).unwrap();
        let mut f = c.field.clone();
        let base = pr::small_x_base(&c, 30).unwrap();
        for i in 0..base.len() {
            for j in 0..base.len() {
                for neg in [false, true] {
                    let q = if neg {
                        c.negate(&mut f, &base[j])
                    } else {
                        base[j]
                    };
                    let Some(s) = c.add_points(&mut f, &base[i], &q) else {
                        continue;
                    };
                    let v = prime_s3(&c, &mut f, base[i].x, base[j].x, s.x);
                    assert!(f.is_zero(v), "S3 != 0 at ({i},{j},{neg})");
                }
            }
        }
    }

    #[test]
    fn binary_s3_polynomial_vanishes_on_group_law_sums_for_a_zero_and_one() {
        for h in [2u64, 4] {
            // h = 2 uses a = 1 and h = 4 uses a = 0: S₃ must not depend on a.
            let c = bin::toy_curve_with_cofactor::<1>(13, 4, h).unwrap();
            let f = &c.f;
            let (base, _) = bin::small_x_base(&c, 30).unwrap();
            for i in 0..base.len() {
                for j in 0..base.len() {
                    for neg in [false, true] {
                        let q = if neg { c.negate(&base[j]) } else { base[j] };
                        let Some(s) = c.add_points(&base[i], &q) else {
                            continue;
                        };
                        let (x1, x2, x3) = (base[i].x, base[j].x, s.x);
                        let (a, b, cc) = (f.mul(&x1, &x2), f.mul(&x1, &x3), f.mul(&x2, &x3));
                        let sq = f.sqr(&f.add(&f.add(&a, &b), &cc));
                        let v = f.add(&f.add(&sq, &f.mul(&a, &x3)), &c.b);
                        assert!(f.is_zero(&v), "h={h}: S3 != 0 at ({i},{j},{neg})");
                    }
                }
            }
        }
    }

    #[test]
    fn binary_decision_equals_the_pair_table_on_every_field_element() {
        for (h, seed) in [(2u64, 3u64), (4, 5)] {
            let c = bin::toy_curve_with_cofactor::<1>(13, seed, h).unwrap();
            let (base, _) = bin::small_x_base(&c, 40).unwrap();
            let keys = bin::pair_sum_keys(&c, &base);
            let ctx = BinaryS3::new(&c, &base);
            let mut hits = 0;
            for x in 1..(1u64 << 13) {
                let w = c.f.embed_u64(x);
                let table = keys.binary_search(&w).is_ok();
                assert_eq!(ctx.decide(&w), table, "h={h}: x = {x}");
                hits += usize::from(table);
            }
            assert!(hits > 100, "h={h}: only {hits} covered abscissae");
        }
    }

    #[test]
    fn prime_decision_equals_the_pair_table_on_every_field_element() {
        for h in [1u64, 2, 4] {
            let c = if h == 1 {
                pr::toy_curve(13, 3).unwrap()
            } else {
                pr::toy_curve_with_cofactor(13, 3, h).unwrap()
            };
            let mut f = c.field.clone();
            let base = pr::small_x_base(&c, 30).unwrap();
            let keys = pr::pair_sum_keys(&c, &base);
            let ctx = PrimeS3::new(&c, &base);
            let p = f.modulus().to_u64_digits()[0];
            let mut hits = 0;
            for x in 0..p {
                let w = f.embed_u64(x);
                let table = keys.binary_search(&f.key(w)).is_ok();
                assert_eq!(ctx.decide(&mut f, w), table, "h={h}: x = {x}");
                hits += usize::from(table);
            }
            assert!(hits > 100, "h={h}: only {hits} covered abscissae");
        }
    }

    #[test]
    fn trace_and_half_trace_tables_match_the_field() {
        let c = BinCurve::<2>::from_secg(secg("sect113r1").unwrap()).unwrap();
        let (base, _) = bin::small_x_base(&c, 4).unwrap();
        let ctx = BinaryS3::new(&c, &base);
        let f = &c.f;
        let mut rng = StdRng::seed_from_u64(2);
        for _ in 0..50 {
            let z = f.random(&mut rng);
            assert_eq!(ctx.trace_is_zero(&z), f.trace(&z) == 0);
            if f.trace(&z) == 0 {
                let t = ctx.half_trace_of(&z);
                assert_eq!(f.add(&f.sqr(&t), &t), z);
            }
        }
    }

    #[test]
    fn s3_cost_on_real_curves_agrees_with_the_table_and_finds_planted_sums() {
        let opts = Options {
            size: 48,
            trials: 20,
            controls: 4,
            seed: 3,
        };
        let c = BinCurve::<2>::from_secg(secg("sect113r1").unwrap()).unwrap();
        let rep = bin::measure(&c, &opts).unwrap();
        let s = s3_cost_binary(&c, &opts, &rep, 3).unwrap();
        assert_eq!((s.disagreements, s.planted_found, s.planted), (0, 16, 16));
        assert!(s.s3_us_per_trial > 0.0 && s.bits_slower_than_rho > 0.0);
        let c = Curve::<pr::MontField<3>>::from_secg(secg("secp160r1").unwrap()).unwrap();
        let rep = pr::measure(&c, &opts).unwrap();
        let s = s3_cost_prime(&c, &opts, &rep, 3).unwrap();
        assert_eq!((s.disagreements, s.planted_found, s.planted), (0, 16, 16));
    }
}
