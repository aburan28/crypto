//! Interval solvers: `Q = [k]G` with `k ∈ [lo, lo + width)`.
//!
//! * [`bsgs`] — baby-step giant-step with the negation map: `⌈√(width/2)⌉`
//!   baby steps stored as 64-bit keys, giant stride `2m + 1`, every hit
//!   verified by a scalar multiplication so a key collision cannot produce
//!   a wrong answer.  `O(√width)` time and memory.
//! * [`kangaroo`] — van Oorschot–Wiener parallel kangaroos with
//!   distinguished points: tame herd from the interval midpoint, wild herd
//!   from `Q`, powers-of-two jumps whose mean is `≈ N√width / 4` for `N`
//!   kangaroos.  `≈ 2√width` group additions and only the distinguished
//!   points in memory.
//!
//! These run on the real NIST curves: every parameter set, at full field
//! width, solves a planted scalar in a `2^20`-wide interval in well under a
//! second.  That exercises the arithmetic and the bookkeeping; it says
//! nothing about full-width keys, whose cost is the whole-group rho
//! estimate in the audit.

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::group::{mix64, EcdlpGroup};
use super::rho::biguint_to_f64;

/// What an interval solve did.
#[derive(Clone, Debug)]
pub struct IntervalReport {
    pub method: &'static str,
    /// `k` with `Q = [k]G`, verified by a scalar multiplication.
    pub scalar: Option<BigUint>,
    /// Group additions performed (scalar multiplications for setup and
    /// verification counted as their double-and-add length).
    pub group_ops: u64,
    /// Entries held (baby steps, or distinguished points).
    pub table_size: usize,
    pub elapsed_ms: u128,
    /// The textbook expectation for the method at this width.
    pub expected_ops: f64,
    pub threads: usize,
}

fn isqrt_ceil(n: &BigUint) -> BigUint {
    if n.is_zero() {
        return BigUint::zero();
    }
    let s = n.sqrt();
    if &s * &s < *n {
        s + BigUint::one()
    } else {
        s
    }
}

/// Baby-step giant-step on `[lo, lo + width)` with the negation map.
pub fn bsgs<G: EcdlpGroup>(
    g: &G,
    target: &G::Elt,
    lo: &BigUint,
    width: &BigUint,
) -> IntervalReport {
    let t0 = Instant::now();
    let n = g.order();
    let gen = g.generator();
    let width_f = biguint_to_f64(width);
    let mut report = IntervalReport {
        method: "bsgs",
        scalar: None,
        group_ops: 0,
        table_size: 0,
        elapsed_ms: 0,
        expected_ops: 2.0 * (width_f / 2.0).sqrt(),
        threads: 1,
    };
    if width.is_zero() {
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }
    let verify = |k: &BigUint| -> bool { g.mul(&gen, k) == *target };

    // T = Q − [lo]G = [k']G with 0 ≤ k' < width.
    let shift = g.mul(&gen, &(lo % n));
    let t_pt = g.add(target, &g.neg(&shift));
    report.group_ops += 2 * lo.bits();
    if g.is_identity(&t_pt) {
        report.scalar = Some(lo.clone());
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }

    // Baby steps: rep(±[j]G) for j = 0..=m, m = ⌈√(width/2)⌉.
    let m_big = isqrt_ceil(&(width >> 1u32)).max(BigUint::one());
    let m: u64 = m_big
        .iter_u64_digits()
        .next()
        .expect("non-zero")
        .min(u64::MAX / 4);
    let mut table: HashMap<u64, u64> = HashMap::with_capacity(m as usize + 1);
    let mut baby = g.identity();
    for j in 0..=m {
        let (rep, _) = g.canonical_sign(&baby);
        table.entry(g.key(&rep)).or_insert(j);
        baby = g.add(&baby, &gen);
        report.group_ops += 1;
    }
    report.table_size = table.len();

    // Giant steps: T − [i]S, S = [2m+1]G.
    let stride = BigUint::from(2 * m + 1);
    let s_neg = g.neg(&g.mul(&gen, &stride));
    report.group_ops += 2 * stride.bits();
    let giants: u64 = {
        let q = (width + &stride - BigUint::one()) / &stride;
        q.iter_u64_digits().next().unwrap_or(0)
    };
    let mut cur = t_pt;
    for i in 0..=giants {
        let (rep, _) = g.canonical_sign(&cur);
        if let Some(&j) = table.get(&g.key(&rep)) {
            let base = &stride * BigUint::from(i);
            for cand in [&base + BigUint::from(j), (&base + n - BigUint::from(j)) % n] {
                if cand < *width {
                    let k = (lo + &cand) % n;
                    if verify(&k) {
                        report.group_ops += 2 * k.bits();
                        report.scalar = Some(k);
                        report.elapsed_ms = t0.elapsed().as_millis();
                        return report;
                    }
                }
            }
        }
        cur = g.add(&cur, &s_neg);
        report.group_ops += 1;
    }
    report.elapsed_ms = t0.elapsed().as_millis();
    report
}

/// Knobs for [`kangaroo`].
#[derive(Clone, Debug)]
pub struct KangarooOptions {
    pub threads: usize,
    pub seed: u64,
    /// Kangaroos per thread (half tame, half wild); at least 2.
    pub herd: usize,
    /// `None` picks `⌊log₂(√width)/2⌋ − 1` clamped to `[0, 40]`.
    pub dp_bits: Option<u32>,
    /// Abort after this many group additions; `None` means
    /// `64 · (2√width + N · 2^dp_bits)`.
    pub max_ops: Option<u64>,
}

impl Default for KangarooOptions {
    fn default() -> Self {
        Self {
            threads: 1,
            seed: 0x5EED_u64,
            herd: 4,
            dp_bits: None,
            max_ops: None,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Tame,
    Wild,
}

struct Trap {
    kind: Kind,
    /// Tame: `P = [offset]G`.  Wild: `P = Q + [offset]G`.
    offset: BigUint,
}

struct KShared {
    traps: Mutex<HashMap<u64, Trap>>,
    done: AtomicBool,
    ops: AtomicU64,
    dps: AtomicU64,
    result: Mutex<Option<BigUint>>,
}

/// Jump sizes `2^0 … 2^(s−1)` whose mean is at least `target_mean`.
fn jump_sizes(target_mean: f64) -> Vec<u64> {
    let mut s = 1u32;
    while s < 62 {
        let mean = ((1u128 << s) - 1) as f64 / s as f64;
        if mean >= target_mean {
            break;
        }
        s += 1;
    }
    (0..s).map(|j| 1u64 << j).collect()
}

/// Parallel kangaroos on `[lo, lo + width)`.
pub fn kangaroo<G: EcdlpGroup>(
    g: &G,
    target: &G::Elt,
    lo: &BigUint,
    width: &BigUint,
    opts: &KangarooOptions,
) -> IntervalReport {
    let t0 = Instant::now();
    let n = g.order().clone();
    let gen = g.generator();
    let width_f = biguint_to_f64(width);
    let threads = opts.threads.max(1);
    let herd = opts.herd.max(2) & !1usize;
    let total_roos = (threads * herd) as f64;
    let sqrt_w = width_f.sqrt();
    let dp_bits = opts
        .dp_bits
        .unwrap_or_else(|| ((sqrt_w.max(2.0).log2() / 2.0).floor() as i64 - 1).clamp(0, 40) as u32);
    let dp_mask: u64 = (1u64 << dp_bits.min(63)) - 1;
    let budget = opts
        .max_ops
        .unwrap_or_else(|| (64.0 * (2.0 * sqrt_w + total_roos * (1u64 << dp_bits) as f64)) as u64);
    let mut report = IntervalReport {
        method: "kangaroo",
        scalar: None,
        group_ops: 0,
        table_size: 0,
        elapsed_ms: 0,
        expected_ops: 2.0 * sqrt_w + total_roos * (1u64 << dp_bits) as f64,
        threads,
    };
    if width.is_zero() {
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }
    if g.mul(&gen, lo) == *target {
        report.scalar = Some(lo % &n);
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }

    // Mean jump ≈ N√W/4 (Pollard's parallel heuristic).
    let mean = (total_roos * sqrt_w / 4.0).max(1.0);
    let sizes = jump_sizes(mean);
    let jumps: Vec<G::Elt> = sizes
        .iter()
        .map(|d| g.mul(&gen, &BigUint::from(*d)))
        .collect();
    let spacing = BigUint::from(((mean / total_roos).ceil() as u64).max(1));
    let mid = lo + (width >> 1u32);

    let shared = KShared {
        traps: Mutex::new(HashMap::new()),
        done: AtomicBool::new(false),
        ops: AtomicU64::new(0),
        dps: AtomicU64::new(0),
        result: Mutex::new(None),
    };

    std::thread::scope(|scope| {
        for t in 0..threads {
            let shared = &shared;
            let jumps = &jumps;
            let sizes = &sizes;
            let n = &n;
            let gen = &gen;
            let mid = &mid;
            let spacing = &spacing;
            scope.spawn(move || {
                let r = jumps.len() as u64;
                // Kangaroo state: (kind, point, offset).
                let mut roos: Vec<(Kind, G::Elt, BigUint)> = Vec::with_capacity(herd);
                for i in 0..herd {
                    let idx = BigUint::from((t * herd + i) as u64 / 2);
                    let start = &idx * spacing;
                    if i % 2 == 0 {
                        let off = (mid + &start) % n;
                        roos.push((Kind::Tame, g.mul(gen, &off), off));
                    } else {
                        let off = start;
                        roos.push((Kind::Wild, g.add(target, &g.mul(gen, &off)), off));
                    }
                }
                let mut local: u64 = 0;
                let mut restart_jitter: u64 = mix64(opts.seed ^ t as u64);
                while !shared.done.load(Ordering::Relaxed) {
                    for (kind, p, off) in roos.iter_mut() {
                        if g.is_identity(p) {
                            // [off]G = O (tame) or Q = −[off]G (wild).
                            if *kind == Kind::Wild {
                                let k = (n - off.clone() % n) % n;
                                if g.mul(gen, &k) == *target {
                                    let mut res = shared.result.lock().unwrap();
                                    if res.is_none() {
                                        *res = Some(k);
                                    }
                                    shared.done.store(true, Ordering::Relaxed);
                                    return;
                                }
                            }
                            // Hop off the identity.
                            *p = g.add(p, &jumps[0]);
                            *off = (&*off + BigUint::from(sizes[0])) % n;
                            local += 1;
                            continue;
                        }
                        let key = g.key(p);
                        if mix64(key ^ 0x4A4B_0000_0000_0001) & dp_mask == 0 {
                            shared.dps.fetch_add(1, Ordering::Relaxed);
                            let mut traps = shared.traps.lock().unwrap();
                            match traps.get(&key) {
                                Some(tr) if tr.kind != *kind => {
                                    // tame_off ≡ k + wild_off.
                                    let (tame, wild) = if *kind == Kind::Tame {
                                        (off.clone(), tr.offset.clone())
                                    } else {
                                        (tr.offset.clone(), off.clone())
                                    };
                                    let k = (tame + n - wild % n) % n;
                                    if g.mul(gen, &k) == *target {
                                        drop(traps);
                                        let mut res = shared.result.lock().unwrap();
                                        if res.is_none() {
                                            *res = Some(k);
                                        }
                                        shared.done.store(true, Ordering::Relaxed);
                                        return;
                                    }
                                }
                                Some(tr) if tr.offset != *off => {
                                    // Same herd, merged trail: move this one.
                                    drop(traps);
                                    restart_jitter = mix64(restart_jitter);
                                    let jitter = BigUint::from(
                                        restart_jitter
                                            % (sizes.last().copied().unwrap_or(1) * 2 + 1),
                                    );
                                    *off = (&*off + &jitter) % n;
                                    *p = g.add(p, &g.mul(gen, &jitter));
                                    local += 2 * jitter.bits();
                                    continue;
                                }
                                Some(_) => {
                                    drop(traps);
                                }
                                None => {
                                    traps.insert(
                                        key,
                                        Trap {
                                            kind: *kind,
                                            offset: off.clone(),
                                        },
                                    );
                                    drop(traps);
                                }
                            }
                        }
                        let j = (key % r) as usize;
                        *p = g.add(p, &jumps[j]);
                        *off = (&*off + BigUint::from(sizes[j])) % n;
                        local += 1;
                    }
                    if local >= 256 {
                        let total = shared.ops.fetch_add(local, Ordering::Relaxed) + local;
                        local = 0;
                        if total >= budget {
                            shared.done.store(true, Ordering::Relaxed);
                            return;
                        }
                    }
                }
                shared.ops.fetch_add(local, Ordering::Relaxed);
            });
        }
    });

    report.scalar = shared.result.into_inner().unwrap();
    report.group_ops = shared.ops.load(Ordering::Relaxed);
    report.table_size = shared.dps.load(Ordering::Relaxed) as usize;
    report.elapsed_ms = t0.elapsed().as_millis();
    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ecdlp_nist::toy;

    #[test]
    fn bsgs_and_kangaroo_agree_on_toy_curves() {
        let g = toy::prime_mid();
        let lo = BigUint::from(40_000u32);
        let width = BigUint::from(20_000u32);
        for k in [40_000u32, 40_001, 49_999, 59_998, 59_999] {
            let q = g.mul(&g.generator(), &BigUint::from(k));
            let b = bsgs(&g, &q, &lo, &width);
            assert_eq!(b.scalar, Some(BigUint::from(k)), "bsgs k={k}: {b:?}");
            let r = kangaroo(&g, &q, &lo, &width, &KangarooOptions::default());
            assert_eq!(r.scalar, Some(BigUint::from(k)), "kangaroo k={k}: {r:?}");
        }
    }

    #[test]
    fn bsgs_rejects_scalars_outside_the_interval() {
        let g = toy::prime_small();
        let q = g.mul(&g.generator(), &BigUint::from(9000u32));
        let b = bsgs(&g, &q, &BigUint::from(100u32), &BigUint::from(1000u32));
        assert!(b.scalar.is_none());
    }

    #[test]
    fn kangaroo_threads_on_koblitz() {
        let g = toy::koblitz(23, 0).expect("toy");
        let lo = BigUint::from(1u32) << 20u32;
        let width = BigUint::from(1u32) << 18u32;
        let k = &lo + BigUint::from(123_456u32);
        let q = g.mul(&g.generator(), &k);
        let opts = KangarooOptions {
            threads: 3,
            ..KangarooOptions::default()
        };
        let r = kangaroo(&g, &q, &lo, &width, &opts);
        assert_eq!(r.scalar, Some(k));
    }
}
