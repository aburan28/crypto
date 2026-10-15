//! `prime-kummer:t=T[,g=G]`: Petit–Kosters–Messeng's factor base for prime
//! curves (PKC 2016, §3) in its Kummer form, as an `ic_framework` plug-in.
//!
//! `L(x) = (x/g)^{2^t} − 1` is the composition of `t` squarings, so its
//! zero set `V = g·μ_{2^t}` has exactly `2^t` elements in `F_p` when
//! `2^t | p − 1` (condition 1 of PKM §3.1) and is the zero set of
//! low-degree maps (condition 2).  The base is every point of the curve
//! whose abscissa is in `V`, both signs, one column an abscissa.  `g`
//! defaults to 1 (`V = μ_{2^t}`).
//!
//! Building it costs `2^t` Legendre symbols and square roots; that work is
//! counted in the factor-base phase (`legendre_symbols`, `sqrt_solves`),
//! as `prime-abscissa`'s is.
//!
//! `prime-isogeny:t=T[,seed=S]` is the same paper's isogeny form, which
//! needs no condition on `p` (PKM note §2.1): `V = x(R' + ⟨Q⟩)` for a
//! point `Q` of order `2^t` on an auxiliary curve `E'/F_p` found by a
//! public-seeded search, the zeros of a chain of `t` Vélu 2-isogeny
//! x-maps.  The search (a residue-class baby-step giant-step filter on
//! every curve tried, Mestre's baby-step giant-step order computation on
//! the survivors) and the build of `⟨Q⟩`, `V` and the chain are group
//! operations on `E'`, a curve over the same field, and are charged to the
//! factor-base phase in the harness's unit; the split by stage and
//! `aux_curves_tried` are recorded as counters.

use super::stages::{FactorBaseBuilder, InstanceCtx, Params};
use crate::cryptanalysis::ic_boundary::{
    CountedGroup, FactorBase, GroupOps, PrimeCurve, PrimeInstance, PrimePoint,
};
use crate::cryptanalysis::residual_walk::pow_mod;

fn mul(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

/// A primitive `2^t`-th root of unity in `F_p`, from the first small
/// base that gives one; `None` when `2^t ∤ p − 1`.
pub fn primitive_two_power_root(p: u64, t: u32) -> Option<u64> {
    if t == 0 || t >= 63 || !(p - 1).is_multiple_of(1u64 << t) {
        return None;
    }
    let e = (p - 1) >> t;
    (2..p.min(1 << 20)).find_map(|h| {
        let z = pow_mod(h, e, p);
        (pow_mod(z, 1u64 << (t - 1), p) != 1).then_some(z)
    })
}

/// The Kummer tower base.
pub struct PrimeKummerBase<'i> {
    pub instance: &'i PrimeInstance,
}

impl FactorBaseBuilder<PrimeCurve> for PrimeKummerBase<'_> {
    fn name(&self) -> &str {
        "prime-kummer"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "PKM's Kummer tower: every point with x ∈ g·μ_{{2^t}}, the zeros of (x/g)^(2^t) − 1, t = {}, \
             g = {}, one column an abscissa",
            params.get("t").unwrap_or("?"),
            params.get("g").unwrap_or("1")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("t", "tower height: |V| = 2^t; needs 2^t | p − 1"),
            ("g", "the coset representative, default 1"),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<PrimeCurve>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<PrimePoint>, String> {
        let curve = &self.instance.curve;
        let p = curve.p;
        let t = params.u64("t")? as u32;
        let g = params.u64_or("g", 1)? % p;
        if g == 0 {
            return Err("g must be non-zero".into());
        }
        let zeta = primitive_two_power_root(p, t)
            .ok_or_else(|| format!("2^{t} does not divide p − 1 = {}", p - 1))?;
        let mut xs: Vec<u64> = Vec::with_capacity(1 << t);
        let mut z = g;
        for _ in 0..(1u64 << t) {
            xs.push(z);
            z = mul(z, zeta, p);
        }
        xs.sort_unstable();
        lift_abscissae(self.instance, &xs, |columns| {
            format!("Kummer tower t = {t}, g = {g}, {columns} of 2^{t} abscissae on the curve")
        })
    }
}

/// The base over a sorted abscissa set: every point of the curve whose
/// abscissa is in `xs`, both signs, one column an abscissa.  The
/// Legendre symbols and square roots it takes are counted in the base's
/// phase (`legendre_symbols`, `sqrt_solves`).
fn lift_abscissae(
    instance: &PrimeInstance,
    xs: &[u64],
    describe: impl Fn(usize) -> String,
) -> Result<FactorBase<PrimePoint>, String> {
    let curve = &instance.curve;
    let p = curve.p;
    let (mut points, mut col_of, mut coef_of) = (Vec::new(), Vec::new(), Vec::new());
    let (mut legendre, mut sqrts) = (0u64, 0u64);
    let mut columns = 0usize;
    for &x in xs {
        let rhs = curve.rhs(x);
        legendre += 1;
        if curve.legendre(rhs) < 0 {
            continue;
        }
        sqrts += 1;
        let Some(y) = curve.sqrt(rhs) else { continue };
        let pt = PrimePoint::affine(x, y.min(if y == 0 { 0 } else { p - y }));
        let neg = curve.neg(pt);
        points.push(pt);
        col_of.push(columns);
        coef_of.push(1);
        if neg != pt {
            points.push(neg);
            col_of.push(columns);
            coef_of.push(instance.r - 1);
        }
        columns += 1;
    }
    if columns == 0 {
        return Err("no abscissa of the tower carries a point".into());
    }
    let mut fb = FactorBase::from_column_map(
        describe(columns),
        points,
        col_of,
        coef_of,
        columns,
        |q| curve.key(q),
        |q| curve.key(&curve.neg(*q)),
        |q| q.x,
    )?;
    fb.cost.count("legendre_symbols", legendre);
    fb.cost.count("sqrt_solves", sqrts);
    Ok(fb)
}

// ── `prime-isogeny:t=T[,seed=S]` ────────────────────────────────────

/// SplitMix64: a fixed, version-independent generator, so the
/// auxiliary-curve search is a function of `(p, t, seed)` alone.
struct SplitMix(u64);

impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in `[lo, hi)` by rejection (no modulo bias).
    fn range(&mut self, lo: u64, hi: u64) -> u64 {
        let n = hi - lo;
        let zone = u64::MAX - u64::MAX % n;
        loop {
            let v = self.next();
            if v < zone {
                return lo + v % n;
            }
        }
    }
}

/// The work of the auxiliary-curve search and the tower build, by kind.
/// Every group operation on `E'` is also on the caller's ledger (and so
/// in the factor-base phase's `S`); these are the same counts split by
/// stage, plus the field work outside the group law.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IsogenyTowerCost {
    /// Curves `(a', b')` drawn, singular ones included.
    pub curves_tried: u64,
    pub singular: u64,
    /// Curves the residue-class baby-step giant-step proved to have
    /// `2^t ∤ #E'`.
    pub filtered_out: u64,
    /// Curves that passed the filter and had `#E'` computed (Mestre).
    pub orders_computed: u64,
    /// Passed the filter but `2^t ∤ #E'` (the filter's false positive).
    pub filter_false_positives: u64,
    /// `#E'` not determined within the draw limit (rejected).
    pub orders_undetermined: u64,
    /// `2^t | #E'` but no point of order `2^t` found (2-Sylow not cyclic
    /// enough, or unlucky draws): rejected.
    pub two_sylow_rejected: u64,
    /// Group additions plus doublings of the search (filter and order
    /// computations, every curve tried).
    pub search_gae: u64,
    /// Of which in the order computations of the survivors.
    pub count_gae: u64,
    /// Group additions plus doublings of the build on the accepted curve:
    /// `Q`, `⟨Q⟩`, `R'`, `V`, the kernels of the chain.
    pub build_gae: u64,
    /// Points drawn on auxiliary curves: Legendre symbols and square
    /// roots behind them.
    pub legendre_symbols: u64,
    pub sqrt_solves: u64,
    /// Vélu point evaluations (one inversion each), outside the group law.
    pub velu_evaluations: u64,
    /// Baby-step table entries stored by the searches (memory, not time).
    pub bsgs_table_entries: u64,
}

/// One 2-isogeny step: kernel abscissa `ξ_j`, `τ_j = 3ξ_j² + a'_j`, and
/// the domain's coefficients.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct VeluStep {
    pub xi: u64,
    pub tau: u64,
    pub a: u64,
    pub b: u64,
}

/// The isogeny tower: `V = x(R' + ⟨Q⟩)` on `E' : y² = x³ + a'x + b'`
/// with `Q` of order `2^t`, the chain of Vélu steps, and `c = x(Φ(R'))`.
#[derive(Clone, Debug)]
pub struct IsogenyTower {
    pub p: u64,
    pub t: u32,
    pub aux: PrimeCurve,
    pub aux_order: u64,
    pub q: PrimePoint,
    pub r: PrimePoint,
    /// `x(R' + kQ)` for `k = 0 … 2^t − 1`, in that order.
    pub v: Vec<u64>,
    pub steps: Vec<VeluStep>,
    pub c: u64,
    pub cost: IsogenyTowerCost,
}

impl IsogenyTower {
    /// The composed x-map `L = L_t ∘ … ∘ L_1` at `x`, or `None` at a pole.
    pub fn eval(&self, x: u64) -> Option<u64> {
        let p = self.p;
        let mut y = x;
        for s in &self.steps {
            let d = (y + p - s.xi) % p;
            if d == 0 {
                return None;
            }
            y = (y + mul(s.tau, pow_mod(d, p - 2, p), p)) % p;
        }
        Some(y)
    }
}

/// Vélu's 2-isogeny with kernel `(ξ, 0)` on a point.
fn velu_map(pt: PrimePoint, xi: u64, tau: u64, p: u64) -> PrimePoint {
    if pt.infinity || pt.x == xi {
        return PrimePoint::INFINITY;
    }
    let d = pow_mod((pt.x + p - xi) % p, p - 2, p);
    let x = (pt.x + mul(tau, d, p)) % p;
    let y = mul(pt.y, (1 + p - mul(tau, mul(d, d, p), p)) % p, p);
    PrimePoint::affine(x, y)
}

fn random_point(e: &PrimeCurve, rng: &mut SplitMix, cost: &mut IsogenyTowerCost) -> PrimePoint {
    loop {
        let x = rng.range(0, e.p);
        let rhs = e.rhs(x);
        cost.legendre_symbols += 1;
        if e.legendre(rhs) < 0 {
            continue;
        }
        cost.sqrt_solves += 1;
        if let Some(y) = e.sqrt(rhs) {
            let y = if rng.next() & 1 == 1 && y != 0 {
                e.p - y
            } else {
                y
            };
            return PrimePoint::affine(x, y);
        }
    }
}

/// The multiples `k ≡ 0 (mod m)` in `[lo, hi]` with `[k]P = O`, found by
/// baby-step giant-step on `P_m = [m]P` over the index range: `(count,
/// first, spacing)` with `spacing = ord(P_m)` when `count ≥ 2`.  Every
/// group operation is on `ops`.
fn class_hits(
    e: &PrimeCurve,
    pt: PrimePoint,
    m: u64,
    lo: u64,
    hi: u64,
    ops: &mut GroupOps,
    table_entries: &mut u64,
) -> (u64, Option<u64>, Option<u64>) {
    let i0 = lo.div_ceil(m);
    let i1 = hi / m;
    if i1 < i0 {
        return (0, None, None);
    }
    let len = i1 - i0 + 1;
    let pm = e.mul(ops, pt, m);
    if e.is_identity(&pm) {
        return (len, Some(i0 * m), (len > 1).then_some(1));
    }
    // Baby steps j·P_m, j = 1 … w − 1; a repeat (identity) gives ord(P_m)
    // outright, and then the hits are its multiples.
    let w = (len as f64).sqrt().ceil().max(1.0) as u64;
    let mut table: std::collections::HashMap<u64, u64> = std::collections::HashMap::new();
    let mut acc = pm;
    let mut order = None;
    for j in 1..w {
        if e.is_identity(&acc) {
            order = Some(j);
            break;
        }
        table.entry(e.key(&acc)).or_insert(j);
        acc = e.add(ops, acc, pm);
    }
    *table_entries += table.len() as u64;
    let hits_of = |o: u64| {
        let first = i0.div_ceil(o) * o;
        if first > i1 {
            (0, None, None)
        } else {
            let n = (i1 - first) / o + 1;
            (n, Some(first * m), (n > 1).then_some(o))
        }
    };
    if let Some(o) = order {
        return hits_of(o);
    }
    // `acc` is now w·P_m.  Giant steps S_g = (i0 + g·w)·P_m.
    let giant = acc;
    let mut s = e.mul(ops, pm, i0);
    let mut hits: Vec<u64> = Vec::new();
    let mut base = i0;
    while base <= i1 && hits.len() < 2 {
        // S_g + j·P_m = O for j ∈ [0, w): S_g = O, or −S_g = j·P_m.
        if e.is_identity(&s) {
            hits.push(base);
        } else if let Some(&j) = table.get(&e.key(&e.neg(s))) {
            if base + j <= i1 {
                hits.push(base + j);
            }
        }
        base += w;
        if base <= i1 {
            s = e.add(ops, s, giant);
        }
    }
    match hits.as_slice() {
        [] => (0, None, None),
        [h] => (1, Some(h * m), None),
        [h1, h2, ..] => {
            // Two hits: ord(P_m) divides their difference (below `len`),
            // and is what is left after stripping every prime factor q
            // while [o/q]·P_m = O.
            let mut o = h2 - h1;
            let mut rest = o;
            let mut primes = Vec::new();
            let mut q = 2;
            while q * q <= rest {
                if rest.is_multiple_of(q) {
                    primes.push(q);
                    while rest.is_multiple_of(q) {
                        rest /= q;
                    }
                }
                q += 1;
            }
            if rest > 1 {
                primes.push(rest);
            }
            for q in primes {
                while o.is_multiple_of(q) && e.is_identity(&e.mul(ops, pm, o / q)) {
                    o /= q;
                }
            }
            hits_of(o)
        }
    }
}

/// `#E'(F_p)` by Mestre's baby-step giant-step: each point narrows the
/// class of the order modulo `M` until one multiple of `M` is left in
/// the Hasse interval.  `None` after `draws` points.
fn mestre_order(
    e: &PrimeCurve,
    rng: &mut SplitMix,
    lo: u64,
    hi: u64,
    draws: u32,
    ops: &mut GroupOps,
    cost: &mut IsogenyTowerCost,
) -> Option<u64> {
    let mut m = 1u64;
    for _ in 0..draws {
        let pt = random_point(e, rng, cost);
        let (n, first, spacing) = class_hits(e, pt, m, lo, hi, ops, &mut cost.bsgs_table_entries);
        match (n, first, spacing) {
            (1, Some(k), _) => return Some(k),
            (_, _, Some(o)) if n >= 2 => m *= o,
            _ => return None, // cannot happen: #E' is always a hit
        }
    }
    None
}

/// Search, from the public seed, an auxiliary curve `E'/F_p` with a
/// point `Q` of order `2^t`, and build the isogeny tower on it.  Every
/// group operation on `E'` goes on `ops`; the split by stage and the
/// field work are in the returned cost.
///
/// The search draws `(a', b')`, and on each curve
/// 1. *filter*: one random point `P` and a baby-step giant-step over the
///    multiples of `2^t` in the Hasse interval for `[k]P = O`.  No hit
///    proves `2^t ∤ #E'` (`#E'` would be one), and the curve is dropped;
/// 2. *count*: on a hit, `#E'` by Mestre's baby-step giant-step;
/// 3. `Q = [#E'/2^t]P'` for random `P'` until `[2^{t−1}]Q ≠ O`.
///
/// Then `⟨Q⟩` by `2^t − 1` additions, `R'` with `2R' ∉ ⟨Q⟩`,
/// `V = x(R' + ⟨Q⟩)` by `2^t` additions, and the chain: `K_j =
/// [2^{t−1−j}]Q_j` by doublings, `ξ_j = x(K_j)`, `τ_j = 3ξ_j² + a'_j`,
/// `a'_{j+1} = a'_j − 5τ_j`, `b'_{j+1} = b'_j − 7ξ_jτ_j`, `Q_{j+1}`,
/// `R'_{j+1}` by Vélu, and `c = x(Φ(R'))`.
pub fn isogeny_tower(
    p: u64,
    t: u32,
    seed: u64,
    ops: &mut GroupOps,
) -> Result<IsogenyTower, String> {
    if !(1..=20).contains(&t) {
        return Err(format!("t = {t} is outside 1 … 20"));
    }
    if p < 64 || (1u64 << t) > p / 4 {
        return Err(format!("p = {p} is too small for t = {t}"));
    }
    let two_t = 1u64 << t;
    let w = (4 * p).isqrt();
    let (lo, hi) = (p + 1 - w, p + 1 + w);
    let mut rng = SplitMix(
        p.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ (u64::from(t) << 56) ^ seed.rotate_left(17),
    );
    let mut cost = IsogenyTowerCost::default();
    let gae = |o: &GroupOps| o.adds + o.doubles;
    let mut tried_limit = 64 * two_t + 4096;
    let (aux, n, q) = loop {
        if tried_limit == 0 {
            return Err(format!(
                "no auxiliary curve with a point of order 2^{t} found"
            ));
        }
        tried_limit -= 1;
        cost.curves_tried += 1;
        let a = rng.range(0, p);
        let b = rng.range(1, p);
        let disc = (mul(4, mul(mul(a, a, p), a, p), p) + mul(27, mul(b, b, p), p)) % p;
        if disc == 0 {
            cost.singular += 1;
            continue;
        }
        let e = PrimeCurve { p, a, b };
        let before = gae(ops);
        let pt = random_point(&e, &mut rng, &mut cost);
        let (hits, _, _) = class_hits(&e, pt, two_t, lo, hi, ops, &mut cost.bsgs_table_entries);
        if hits == 0 {
            cost.filtered_out += 1;
            cost.search_gae += gae(ops) - before;
            continue;
        }
        let count_before = gae(ops);
        cost.orders_computed += 1;
        let order = mestre_order(&e, &mut rng, lo, hi, 16, ops, &mut cost);
        cost.count_gae += gae(ops) - count_before;
        cost.search_gae += gae(ops) - before;
        let Some(n) = order else {
            cost.orders_undetermined += 1;
            continue;
        };
        if !n.is_multiple_of(two_t) {
            cost.filter_false_positives += 1;
            continue;
        }
        let build_before = gae(ops);
        let mut q = None;
        for _ in 0..8 {
            let cand = e.mul(ops, random_point(&e, &mut rng, &mut cost), n / two_t);
            let mut half = cand;
            for _ in 0..t - 1 {
                half = e.double(ops, half);
            }
            if !e.is_identity(&half) {
                q = Some(cand);
                break;
            }
        }
        cost.build_gae += gae(ops) - build_before;
        match q {
            Some(q) => break (e, n, q),
            None => cost.two_sylow_rejected += 1,
        }
    };
    let build_before = gae(ops);
    // ⟨Q⟩ by successive additions.
    let mut multiples = Vec::with_capacity(two_t as usize);
    multiples.push(PrimePoint::INFINITY);
    let mut acc = q;
    for _ in 1..two_t {
        multiples.push(acc);
        acc = aux.add(ops, acc, q);
    }
    if !aux.is_identity(&acc) {
        return Err("Q does not have order 2^t".into());
    }
    let keys: std::collections::HashSet<u64> = multiples.iter().map(|m| aux.key(m)).collect();
    let r0 = loop {
        let r = random_point(&aux, &mut rng, &mut cost);
        let r2 = aux.double(ops, r);
        if !keys.contains(&aux.key(&r2)) {
            break r;
        }
    };
    let v: Vec<u64> = multiples
        .iter()
        .map(|&m| aux.add(ops, r0, m))
        .map(|pt| pt.x)
        .collect();
    // The chain.
    let (mut aj, mut bj) = (aux.a, aux.b);
    let (mut qj, mut rj) = (q, r0);
    let mut steps = Vec::with_capacity(t as usize);
    for j in 0..t {
        let cj = PrimeCurve { p, a: aj, b: bj };
        let mut kj = qj;
        for _ in 0..(t - 1 - j) {
            kj = cj.double(ops, kj);
        }
        if kj.infinity || kj.y != 0 {
            return Err(format!("step {j}: the kernel point is not 2-torsion"));
        }
        let xi = kj.x;
        let tau = (mul(3, mul(xi, xi, p), p) + aj) % p;
        steps.push(VeluStep {
            xi,
            tau,
            a: aj,
            b: bj,
        });
        qj = velu_map(qj, xi, tau, p);
        rj = velu_map(rj, xi, tau, p);
        cost.velu_evaluations += 2;
        aj = (aj + p - mul(5, tau, p)) % p;
        bj = (bj + p - mul(7, mul(xi, tau, p), p)) % p;
        let next = PrimeCurve { p, a: aj, b: bj };
        if !next.is_on_curve(qj) || !next.is_on_curve(rj) {
            return Err(format!("step {j}: a Vélu image is off the codomain"));
        }
    }
    cost.build_gae += gae(ops) - build_before;
    let tower = IsogenyTower {
        p,
        t,
        aux,
        aux_order: n,
        q,
        r: r0,
        v,
        steps,
        c: rj.x,
        cost,
    };
    let mut distinct = tower.v.clone();
    distinct.sort_unstable();
    distinct.dedup();
    if distinct.len() != two_t as usize {
        return Err(format!("|V| = {} ≠ 2^{t}", distinct.len()));
    }
    Ok(tower)
}

/// The isogeny tower base: every point of the curve whose abscissa is in
/// `V = x(R' + ⟨Q⟩)`, both signs, one column an abscissa.
pub struct PrimeIsogenyBase<'i> {
    pub instance: &'i PrimeInstance,
}

impl FactorBaseBuilder<PrimeCurve> for PrimeIsogenyBase<'_> {
    fn name(&self) -> &str {
        "prime-isogeny"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "PKM's isogeny tower: every point with x ∈ x(R' + ⟨Q⟩) for Q of order 2^t on a searched \
             auxiliary curve, the zeros of a chain of t Vélu 2-isogeny x-maps, t = {}, seed = {}; \
             the search is charged to the factor-base phase",
            params.get("t").unwrap_or("?"),
            params.get("seed").unwrap_or("0")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("t", "tower height: |V| = 2^t; no condition on p"),
            (
                "seed",
                "public seed of the auxiliary-curve search, default 0",
            ),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<PrimeCurve>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<FactorBase<PrimePoint>, String> {
        let p = self.instance.curve.p;
        let t = params.u64("t")? as u32;
        let seed = params.u64_or("seed", 0)?;
        let mut aux_ops = GroupOps::default();
        let tower = isogeny_tower(p, t, seed, &mut aux_ops)?;
        let mut xs = tower.v.clone();
        xs.sort_unstable();
        let (aa, ab, an, tried) = (
            tower.aux.a,
            tower.aux.b,
            tower.aux_order,
            tower.cost.curves_tried,
        );
        let mut fb = lift_abscissae(self.instance, &xs, |columns| {
            format!(
                "isogeny tower t = {t}, seed = {seed}, E' = (a' {aa}, b' {ab}), #E' = {an}, \
                 {tried} auxiliary curves tried, {columns} of 2^{t} abscissae on the curve"
            )
        })?;
        let c = &tower.cost;
        fb.cost.count("legendre_symbols", c.legendre_symbols);
        fb.cost.count("sqrt_solves", c.sqrt_solves);
        for (name, v) in [
            ("aux_curves_tried", c.curves_tried),
            ("aux_curves_singular", c.singular),
            ("aux_curves_filtered_out", c.filtered_out),
            ("aux_orders_computed", c.orders_computed),
            ("aux_filter_false_positives", c.filter_false_positives),
            ("aux_orders_undetermined", c.orders_undetermined),
            ("aux_two_sylow_rejected", c.two_sylow_rejected),
            ("aux_search_gae", c.search_gae),
            ("aux_count_gae", c.count_gae),
            ("aux_tower_build_gae", c.build_gae),
            ("aux_legendre_symbols", c.legendre_symbols),
            ("aux_sqrt_solves", c.sqrt_solves),
            ("aux_bsgs_table_entries", c.bsgs_table_entries),
            ("velu_evaluations_uncharged", c.velu_evaluations),
        ] {
            fb.cost.count(name, v);
        }
        ops.merge(aux_ops);
        Ok(fb)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::prime_instance_for;

    #[test]
    fn the_tower_is_the_two_power_roots_of_unity() {
        // 2013265921 = 15·2^27 + 1, the PKM rounds' p₁.
        let p = 2_013_265_921u64;
        let z = primitive_two_power_root(p, 10).unwrap();
        assert_eq!(pow_mod(z, 1 << 10, p), 1);
        assert_ne!(pow_mod(z, 1 << 9, p), 1);
        assert!(primitive_two_power_root(p, 28).is_none());
    }

    #[test]
    fn the_base_holds_exactly_the_tower_points() {
        // Any small prime curve with 2^t | p − 1 for some t ≥ 2.
        let mut seed = 1;
        let (inst, t) = loop {
            let inst = prime_instance_for(16, seed);
            let t = (inst.curve.p - 1).trailing_zeros().min(6);
            if t >= 3 {
                break (inst, t);
            }
            seed += 1;
        };
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: inst.generator_point(),
            target: inst.generator_point(),
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        };
        let (_, params) = Params::parse_spec(&format!("prime-kummer:t={t}")).unwrap();
        let fb = PrimeKummerBase { instance: &inst }
            .build(&ctx, &params, &mut GroupOps::default())
            .unwrap();
        let p = inst.curve.p;
        for q in &fb.points {
            assert!(inst.curve.is_on_curve(*q));
            assert_eq!(pow_mod(q.x, 1 << t, p), 1, "abscissa outside μ_2^t");
        }
        assert!(fb.abscissae <= 1 << t);
    }

    #[test]
    fn the_isogeny_tower_has_exactly_two_to_the_t_roots() {
        // A prime with v₂(p − 1) = 1, where no Kummer tower exists, and
        // one with deeper 2-adicity: the isogeny tower needs neither.
        for (p, t) in [(16_411u64, 5u32), (40_961, 6), (65_519, 4)] {
            let mut ops = GroupOps::default();
            let tower = isogeny_tower(p, t, 7, &mut ops).unwrap();
            assert_eq!(tower.v.len(), 1 << t);
            let mut v = tower.v.clone();
            v.sort_unstable();
            v.dedup();
            assert_eq!(v.len(), 1 << t, "|V| = 2^t, distinct");
            for &x in &tower.v {
                assert_eq!(tower.eval(x), Some(tower.c), "x = {x} is not a root");
            }
            // And nothing else is: L(x) = c has exactly the 2^t roots V.
            let roots = (0..p).filter(|&x| tower.eval(x) == Some(tower.c)).count();
            assert_eq!(roots, 1 << t, "p = {p}");
            assert_eq!(tower.steps.len(), t as usize);
            assert_eq!(tower.aux_order, tower.aux.point_count(), "Mestre order");
            assert!(tower.aux_order.is_multiple_of(1 << t));
            // The charge: every group operation of the search and build is
            // on the ledger, and the split adds up to it.
            let c = &tower.cost;
            assert_eq!(ops.adds + ops.doubles, c.search_gae + c.build_gae);
            assert!(c.curves_tried >= 1 && c.search_gae > 0 && c.build_gae > 0);
            assert_eq!(
                c.curves_tried,
                1 + c.singular
                    + c.filtered_out
                    + c.filter_false_positives
                    + c.orders_undetermined
                    + c.two_sylow_rejected
            );
            // Deterministic in (p, t, seed); another seed, another curve.
            let again = isogeny_tower(p, t, 7, &mut GroupOps::default()).unwrap();
            assert_eq!(again.v, tower.v);
            assert_eq!(again.cost, tower.cost);
            let other = isogeny_tower(p, t, 8, &mut GroupOps::default()).unwrap();
            assert_ne!(
                (other.aux.a, other.aux.b, other.r),
                (tower.aux.a, tower.aux.b, tower.r)
            );
        }
    }

    #[test]
    fn mestre_orders_agree_with_point_counting() {
        let p = 10_007u64;
        let mut rng = SplitMix(3);
        let w = (4 * p).isqrt();
        let mut cost = IsogenyTowerCost::default();
        let mut done = 0;
        while done < 40 {
            let e = PrimeCurve {
                p,
                a: rng.range(0, p),
                b: rng.range(1, p),
            };
            let disc =
                (mul(4, mul(mul(e.a, e.a, p), e.a, p), p) + mul(27, mul(e.b, e.b, p), p)) % p;
            if disc == 0 {
                continue;
            }
            let mut ops = GroupOps::default();
            let n = e.point_count();
            if let Some(m) =
                mestre_order(&e, &mut rng, p + 1 - w, p + 1 + w, 16, &mut ops, &mut cost)
            {
                assert_eq!(m, n);
            }
            // The filter never drops a curve whose order is a multiple.
            for t in 1..6u32 {
                let pt = random_point(&e, &mut rng, &mut cost);
                let (hits, _, _) =
                    class_hits(&e, pt, 1 << t, p + 1 - w, p + 1 + w, &mut ops, &mut 0);
                if n.is_multiple_of(1 << t) {
                    assert!(hits >= 1);
                }
            }
            done += 1;
        }
    }

    #[test]
    fn the_isogeny_base_holds_exactly_the_tower_points_that_lift() {
        let inst = prime_instance_for(16, 3);
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: inst.generator_point(),
            target: inst.generator_point(),
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        };
        let t = 5;
        let (_, params) = Params::parse_spec(&format!("prime-isogeny:t={t},seed=2")).unwrap();
        let mut ops = GroupOps::default();
        let fb = PrimeIsogenyBase { instance: &inst }
            .build(&ctx, &params, &mut ops)
            .unwrap();
        let p = inst.curve.p;
        let tower = isogeny_tower(p, t, 2, &mut GroupOps::default()).unwrap();
        let lifting: Vec<u64> = {
            let mut v: Vec<u64> = tower
                .v
                .iter()
                .copied()
                .filter(|&x| inst.curve.legendre(inst.curve.rhs(x)) >= 0)
                .collect();
            v.sort_unstable();
            v
        };
        assert_eq!(
            fb.abscissa_list, lifting,
            "the base's x-set is V ∩ x(E(F_p))"
        );
        for q in &fb.points {
            assert!(inst.curve.is_on_curve(*q));
        }
        assert_eq!(fb.columns, lifting.len());
        // The search is charged on the build's ledger and recorded apart.
        assert_eq!(
            ops.adds + ops.doubles,
            fb.cost.get("aux_search_gae") + fb.cost.get("aux_tower_build_gae")
        );
        assert!(fb.cost.get("aux_curves_tried") >= 1);
        assert!(PrimeIsogenyBase { instance: &inst }
            .build(
                &ctx,
                &Params::parse_spec("prime-isogeny:t=40").unwrap().1,
                &mut ops
            )
            .is_err());
    }
}
