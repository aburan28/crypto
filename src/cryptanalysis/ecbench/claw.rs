//! The **four-point signed-Frobenius pair claw** of
//! `aburan28/cryptanalysis#175` (`experiments/koblitz-pair-claw-20260929`,
//! the `PDP4qtable … LAnone` candidate), ported onto the ecbench ledger.
//!
//! The method, as that experiment runs it on `EC1N83Ckb1h876c2921cb64`:
//!
//! 1. **Known-log base.**  `S` seed scalars `a_i`, seed points `[a_i]G`,
//!    and every signed Frobenius conjugate `σφ^t[a_i]G`, whose logarithm
//!    is `σλ^t a_i` because `φ = [λ]` on the subgroup.  `B = 2nS` points,
//!    `S` folded columns.
//! 2. **Table** (target-independent).  Zero-pair quotient descriptors
//!    `P_i + σφ^t P_j` (`i ≤ j`, first point the unrotated seed), visited
//!    in a fixed permutation without repeats, keyed by the least
//!    normal-basis rotation of the sum's abscissa — the `{±φ^u}` class.
//! 3. **Queries** (online).  Unordered pairs `{F_k, F_l}` of base points,
//!    visited in a fixed permutation without repeats; each costs
//!    `F_k + F_l` and `Q − (F_k + F_l)`, then one canonicalisation and one
//!    probe.  A hit is a natural four-point relation
//!    `Q = F_k + F_l ± φ^u(P_i + σφ^t P_j)`, and every term's logarithm is
//!    known, so the scalar follows with no linear algebra.
//!
//! Because every base logarithm is known, this is a **generic** algorithm:
//! a randomised baby-step giant-step on signed Frobenius classes.  With a
//! table of `M` classes each query lands in it with probability
//! `2nM / r`, so with `M = c·√(r/n)` the expected cost is
//! `M + 2·r/(2nM) = (c + 1/c)·√(r/n)` additions, `S ≈ (c + 1/c)/√n`,
//! minimised at `c = 1` (`2/√n`), against the signed-Frobenius rho
//! constant `√(π/4n)`: about `2.26×` rho at best.  `expected_s` stays
//! `None` because the constant depends on `c`; the protocol states it.
//!
//! Charged: the seeds' scalar multiplications, every table and query
//! addition, and the addition that rebuilds the hit's table sum.
//! Counted, not charged: Frobenius maps, canonicalisations, table inserts
//! and probes (the memory the table holds is `M` entries, as BSGS's is).
//! The planted scalar is never seen; the candidate is returned unverified.

use std::collections::{BTreeMap, HashMap};

use num_traits::ToPrimitive;

use crate::cryptanalysis::ecbench::canonical::derive_u64;
use crate::cryptanalysis::ic_boundary::{BinaryGroup, BinaryInstance, CountedGroup, GroupOps};
use crate::cryptanalysis::ic_measurement::{self as measurement, Phase};
use crate::cryptanalysis::koblitz_fast::{FastPoint, FrobeniusCanon, FrobeniusPowers};

/// The two shape parameters, as multiples of their natural scales.
#[derive(Clone, Copy, Debug)]
pub struct ClawShape {
    /// Seed orbits `S = ⌈base_scale · r^{1/4} n^{−3/4}⌉`.  The PR's n=83
    /// base (48,194 seeds) is `base_scale ≈ 1.06`.
    pub base_scale: f64,
    /// Table classes `M = ⌈table_scale · √(r/n)⌉`.  `1` minimises the
    /// expected additions; the PR's memory-bound `2^31`–`2^33` tables at
    /// n=83 are `table_scale ≈ 1/32`.
    pub table_scale: f64,
    /// Total charged additions allowed before the run is exhausted.
    pub max_adds: u64,
}

impl ClawShape {
    pub fn seeds(&self, r: u64, n: u32) -> u64 {
        let s = self.base_scale * (r as f64).powf(0.25) * (n as f64).powf(-0.75);
        (s.ceil() as u64).max(2)
    }
    pub fn table(&self, r: u64, n: u32) -> u64 {
        let m = self.table_scale * (r as f64 / n as f64).sqrt();
        (m.ceil() as u64).max(1)
    }
}

/// What one claw solve did, phase by phase.
#[derive(Clone, Debug, Default)]
pub struct ClawOutcome {
    pub recovered: Option<u64>,
    pub base: GroupOps,
    pub table: GroupOps,
    pub search: GroupOps,
    pub recover: GroupOps,
    pub counters: BTreeMap<String, u64>,
    pub exhausted: bool,
}

impl ClawOutcome {
    fn count(&mut self, name: &str, by: u64) {
        *self.counters.entry(name.to_string()).or_insert(0) += by;
    }
}

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn addmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + b as u128) % m as u128) as u64
}

fn negmod(a: u64, m: u64) -> u64 {
    if a % m == 0 {
        0
    } else {
        m - a % m
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// A fixed permutation of `[0, d)`: `i ↦ (a·i + b) mod d` with
/// `gcd(a, d) = 1`, so a prefix visits distinct indices — the PR's
/// "unique schedule" without storing it.
struct Schedule {
    d: u64,
    a: u64,
    b: u64,
}

impl Schedule {
    fn new(d: u64, seed: u64, domain: &str) -> Self {
        let mut a = derive_u64(domain, &[seed, 0]) % d.max(1);
        while d > 1 && (a == 0 || gcd(a, d) != 1) {
            a = (a + 1) % d;
        }
        let b = derive_u64(domain, &[seed, 1]) % d.max(1);
        Self { d, a: a.max(1), b }
    }
    fn at(&self, i: u64) -> u64 {
        ((self.a as u128 * i as u128 + self.b as u128) % self.d as u128) as u64
    }
}

/// The class of a point under `{±φ^u}`: the orbit name of its abscissa
/// and the shift `u` that carries the abscissa to the least rotation.
fn class_of(canon: &FrobeniusCanon, p: &FastPoint) -> (u64, u32) {
    let (key, t) = canon.canon_with_shift(p.x);
    // Offset by one so that no finite point shares a key with O.
    (key + 1, t)
}

/// Run the claw on a Koblitz instance.  `None` for any other curve.
pub fn pair_claw(
    inst: &BinaryInstance,
    target: FastPoint,
    seed: u64,
    shape: ClawShape,
) -> Option<ClawOutcome> {
    let kc = inst.koblitz.as_ref()?;
    if kc.k != 1 {
        return None;
    }
    let g = BinaryGroup(&inst.fast);
    let (n, r) = (inst.n, inst.r);
    let lambda = (&kc.lambda % num_bigint::BigUint::from(r)).to_u64()?;
    let canon = FrobeniusCanon::new(&inst.fast.field, n)?;
    let powers = FrobeniusPowers::new(&inst.fast.field, n);
    let mut lambda_pow = Vec::with_capacity(n as usize);
    let mut cur = 1u64;
    for _ in 0..n {
        lambda_pow.push(cur);
        cur = mulmod(cur, lambda, r);
    }
    let frob = |t: u32, p: FastPoint| -> FastPoint {
        if t == 0 || p.infinity {
            p
        } else {
            FastPoint::affine(powers.apply(t, p.x), powers.apply(t, p.y))
        }
    };

    let mut out = ClawOutcome::default();
    let s = shape.seeds(r, n);
    let m = shape.table(r, n);
    out.count("seed_orbits", s);
    out.count("table_target", m);

    // 1. Known-log base: B = 2nS points, index b ↦ (i, t, σ).
    let two_n = 2 * n as u64;
    let b_count = two_n * s;
    out.count("base_points", b_count);
    let mut seeds: Vec<u64> = Vec::with_capacity(s as usize);
    let mut base: Vec<FastPoint> = Vec::with_capacity(b_count as usize);
    for i in 0..s {
        let a = 1 + derive_u64("ecbench.claw.seed", &[seed, i]) % (r - 1);
        seeds.push(a);
        let p = g.mul(&mut out.base, inst.generator, a);
        for t in 0..n {
            let q = frob(t, p);
            base.push(q);
            base.push(g.neg(q));
        }
    }
    out.count("frobenius_maps_uncharged", s * (n as u64 - 1));
    let log_of = |b: u64| -> u64 {
        let i = (b / two_n) as usize;
        let t = ((b % two_n) / 2) as usize;
        let l = mulmod(seeds[i], lambda_pow[t], r);
        if b % 2 == 1 {
            negmod(l, r)
        } else {
            l
        }
    };
    let charged = |o: &ClawOutcome| o.base.gae() + o.table.gae() + o.search.gae();

    // 2. Table: descriptor d ∈ [0, S·B) ↦ (i = d / B, b = d mod B), the sum
    // P_i + F_b with i ≤ seed(b); keyed by class, first descriptor wins.
    let d_table = s * b_count;
    let sched = Schedule::new(d_table, seed, "ecbench.claw.table");
    let mut table: HashMap<u64, u64> = HashMap::with_capacity(m as usize);
    let mut visited = 0u64;
    while (table.len() as u64) < m && visited < d_table {
        let d = sched.at(visited);
        visited += 1;
        let (i, b) = (d / b_count, d % b_count);
        if i > b / two_n {
            out.count("table_descriptors_skipped", 1);
            continue;
        }
        let sum = g.add(&mut out.table, base[(i * two_n) as usize], base[b as usize]);
        if sum.infinity {
            out.count("table_identity_sums", 1);
            continue;
        }
        let (key, _) = class_of(&canon, &sum);
        out.count("canonicalisations_uncharged", 1);
        out.count("inserts_uncharged", 1);
        // The first descriptor of a class is kept, deterministically.
        match table.entry(key) {
            std::collections::hash_map::Entry::Occupied(_) => {
                out.count("table_duplicate_classes", 1);
            }
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert(d);
            }
        }
        if charged(&out) as u64 >= shape.max_adds {
            break;
        }
    }
    out.count("table_classes", table.len() as u64);
    out.count("table_descriptors_visited", visited);

    // 3. Queries: unordered pairs k < l from [0, B)², online.
    let d_query = b_count * b_count;
    let qsched = Schedule::new(d_query, seed, "ecbench.claw.query");
    let mut qi = 0u64;
    // The online window: from the first target-dependent addition to the
    // returned candidate.  The base and the table are reusable set-up.
    measurement::begin_online(Phase::RhoSolve);
    while qi < d_query {
        if charged(&out) as u64 >= shape.max_adds {
            out.exhausted = true;
            out.count("queries_visited", qi);
            measurement::end_online();
            return Some(out);
        }
        let d = qsched.at(qi);
        qi += 1;
        let (k, l) = (d / b_count, d % b_count);
        if k >= l {
            out.count("query_pairs_skipped", 1);
            continue;
        }
        out.count("queries", 1);
        let pair = g.add(&mut out.search, base[k as usize], base[l as usize]);
        let y = g.add(&mut out.search, target, g.neg(pair));
        let pair_log = addmod(log_of(k), log_of(l), r);
        if y.infinity {
            out.count("queries_visited", qi);
            out.count("direct_pair_hits", 1);
            out.recovered = Some(pair_log);
            measurement::end_online();
            return Some(out);
        }
        let (key, ty) = class_of(&canon, &y);
        out.count("canonicalisations_uncharged", 1);
        out.count("lookups_uncharged", 1);
        let Some(&d) = table.get(&key) else { continue };
        // 4. Recover: rebuild the table sum T, align the two conjugates and
        // read the sign off the ordinate:  φ^{tT}T = σ·φ^{ty}Y.
        let (i, b) = (d / b_count, d % b_count);
        let t_sum = g.add(
            &mut out.recover,
            base[(i * two_n) as usize],
            base[b as usize],
        );
        let (_, tt) = class_of(&canon, &t_sum);
        out.count("canonicalisations_uncharged", 1);
        let (ct, cy) = (frob(tt, t_sum), frob(ty, y));
        out.count("frobenius_maps_uncharged", 4);
        let sign_neg = if ct == cy {
            false
        } else if ct == g.neg(cy) {
            true
        } else {
            out.count("key_mismatch", 1);
            continue;
        };
        // Y = σ φ^{tT − tY} T, so log Y = σ λ^{(tT − tY) mod n} log T.
        let u = ((tt + n - ty) % n) as usize;
        let t_log = addmod(log_of(i * two_n), log_of(b), r);
        let mut y_log = mulmod(lambda_pow[u], t_log, r);
        if sign_neg {
            y_log = negmod(y_log, r);
        }
        out.count("queries_visited", qi);
        out.count("relations", 1);
        out.recovered = Some(addmod(pair_log, y_log, r));
        measurement::end_online();
        return Some(out);
    }
    measurement::end_online();
    out.count("queries_visited", qi);
    out.count("query_domain_exhausted", 1);
    out.exhausted = true;
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::koblitz_instance;

    fn shape(base_scale: f64, table_scale: f64) -> ClawShape {
        ClawShape {
            base_scale,
            table_scale,
            max_adds: u64::MAX,
        }
    }

    #[test]
    fn claw_recovers_planted_logs_on_koblitz_curves() {
        for (a, n) in [(0u8, 13u32), (0, 41)] {
            let inst = koblitz_instance(a, n).expect("Koblitz instance");
            let g = BinaryGroup(&inst.fast);
            for (j, k) in [3u64, inst.r - 2, inst.r / 5 + 11].into_iter().enumerate() {
                let mut ops = GroupOps::default();
                let q = g.mul(&mut ops, inst.generator, k);
                for sh in [shape(4.0, 1.0), shape(1.06, 1.0 / 32.0)] {
                    let o = pair_claw(&inst, q, 7 + j as u64, sh).expect("koblitz");
                    if n == 13 && o.exhausted {
                        continue; // a 13-bit group can exhaust a tiny base
                    }
                    assert_eq!(o.recovered, Some(k), "n={n} k={k} {sh:?} {:?}", o.counters);
                }
            }
        }
    }

    #[test]
    fn claw_is_deterministic_in_seed() {
        let inst = koblitz_instance(0, 41).unwrap();
        let g = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        let q = g.mul(&mut ops, inst.generator, 123_456_789);
        let a = pair_claw(&inst, q, 99, shape(4.0, 1.0)).unwrap();
        let b = pair_claw(&inst, q, 99, shape(4.0, 1.0)).unwrap();
        assert_eq!(a.counters, b.counters);
        assert_eq!(a.search.adds, b.search.adds);
    }
}
